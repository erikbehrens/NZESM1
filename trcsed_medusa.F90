MODULE trcsed_medusa
   !!======================================================================
   !!                         ***  MODULE trcsed_medusa  ***
   !! TOP :   MEDUSA Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :    -   !  1995-06 (M. Levy)  original code
   !!              -   !  2000-12 (E. Kestenare)  clean up
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90 + simplifications
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!              -   !  2011-04  (A. Yool) updated for ROAM project
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'                                      MEDUSA bio-model
   !!----------------------------------------------------------------------
   !!   trc_sed_medusa        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc
   USE sms_medusa
   !! AXY (10/02/09)
   USE iom
   !! USE trc_nam_dia                   ! JPALM 13-11-2015 -- if iom_use for diag
   !! USE trc_nam_iom_medusa            ! JPALM 13-11-2015 -- if iom_use for diag
   USE fldread                          !  time interpolation
   USE lbclnk
   USE prtctl_trc                       ! Print control for debbuging
   !! JPALM (27-06-2016): add lk_oasis for CO2 and DMS coupling with atm
   USE sbc_oce, ONLY: lk_oasis
   USE oce,     ONLY: Dust_in_cpl
   !! Check Dust dep
# if defined key_debug_medusa   
   !! USE trcrst, ONLY: trc_rst_dia_stat   !! variable stat
# endif


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sed_medusa     ! called in ???
   PUBLIC   trc_sed_medusa_sbc     
   PUBLIC   trc_sed_medusa_dust 

   !! * Module variables
   !! INTEGER ::                   &
     !! ryyss,                     &  !: number of seconds per year
     !! rmtss                         !: number of seconds per month

   !! AXY (10/02/09)
   LOGICAL, PUBLIC  ::   bdustfer  !: boolean for dust input from the atmosphere
   REAL(wp), PUBLIC ::                 &
      sedfeinput = 1.e-9_wp  ,         &
      dustsolub  = 0.014_wp
   REAL(wp), PARAMETER :: Fe_dust_mratio = 0.035   !! Fe:dust mass ratio = 0.035
   INTEGER , PARAMETER ::        nbtimes = 365     !: maximum number of times record in a file
   INTEGER  :: ntimes_dust                         ! number of time steps in a file

   INTEGER ::                          &
      numdust,                         &
      nflx1,  nflx2,                   &
      nflx11, nflx12
   
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_dust      ! structure of input dust

   
   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_sed_medusa( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sed_medusa  ***
      !!
      !! ** Purpose :   compute the now trend due to the vertical sedimentation of
      !!              detritus and add it to the general trend of detritus equations
      !!
      !! ** Method  :   this ROUTINE compute not exactly the advection but the
      !!              transport term, i.e.  dz(wt) and dz(ws)., dz(wtr)
      !!              using an upstream scheme
      !!              the now vertical advection of tracers is given by:
      !!                      dz(trn wn) = 1/bt dk+1( e1t e2t vsed (trn) )
      !!              add this trend now to the general trend of tracer (ta,sa,tra):
      !!                             tra = tra + dz(trn wn)
      !!        
      !!              IF 'key_trc_diabio' is defined, the now vertical advection
      !!              trend of passive tracers is saved for futher diagnostics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !! AXY (10/02/09)
      INTEGER  ::   jnt
      !!
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   ztra
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwork

      !! AXY (10/02/09)
      REAL(wp) ::   rfact2

      CHARACTER (len=25) :: charout
      
      !! JPALM - 26-11-2015 -add iom_use for diagnostic
       REAL(wp), POINTER, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !!
      IF( lk_iomput) THEN  
           IF( med_diag%DSED%dgsave ) THEN
               CALL wrk_alloc( jpi, jpj,      zw2d )
                zw2d(:,:)      = 0.0      !!
           ENDIF
      ENDIF
      
      !! AXY (10/02/09)
      jnt = 1
      rfact2 = 1.0

      ! Number of seconds per year and per month
      !! ryyss = nyear_len(1) * rday
      !! rmtss = ryyss / raamo

      !! AXY (20/11/14): alter this to report on first MEDUSA call
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_sed_medusa: MEDUSA sedimentation'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
	 IF(lwp) WRITE(numout,*) ' kt =',kt
      ENDIF

      ! sedimentation of detrital nitrogen : upstream scheme
      ! ----------------------------------------------------
      !
      zwork(:,:,:) = 0.e0        ! initialisation of sinking variable
      ! for detrital nitrogen sedimentation only - jpdet
      zwork(:,:,1  ) = 0.e0      ! surface value set to zero
      !!    DO ji = 1, jpi
      !!       zirondep(ji,jj,1) = (dustsolub * dust(ji,jj) / (55.85 * rmtss) + 3.e-10 / ryyss) &
      !!       & * rfact2 / fse3t(ji,jj,1)
      !!       zsidep  (ji,jj)   = 8.8 * 0.075 * dust(ji,jj) * rfact2 / &
      !!       & (fse3t(ji,jj,1) * 28.1 * rmtss)
      !!    END DO
      !! END DO

      ! sedimentation of detrital nitrogen : upstream scheme
      ! ----------------------------------------------------
      !
      zwork(:,:,:) = 0.e0        ! initialisation of sinking variable
      ! for detrital nitrogen sedimentation only - jpdet
      zwork(:,:,1  ) = 0.e0      ! surface value set to zero
      zwork(:,:,jpk) = 0.e0      ! bottom value  set to zero
      !
      ! tracer flux at w-point: we use -vsed (downward flux)  with simplification : no e1*e2
      DO jk = 2, jpk
         ! AXY (17/07/14): change "0.d0" to "0."
         ! zwork(:,:,jk) = -vsed * max(trn(:,:,jk-1,jpdet),0.d0) * tmask(:,:,jk-1)
         zwork(:,:,jk) = -vsed * max(trn(:,:,jk-1,jpdet),0.) * tmask(:,:,jk-1)
         !
         ! AXY (16/01/14): stop sinking in upper 10m to reduce model instability 
         !                 in shallower grid cells
         ! if ( jk .lt. 9 ) zwork(:,:,jk) = 0.e0
      END DO
      !
      ! tracer flux divergence at t-point added to the general trend
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               ztra  = - ( zwork(ji,jj,jk) - zwork(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + ztra
               IF( med_diag%DSED%dgsave ) THEN
                   zw2d(ji,jj) = zw2d(ji,jj) + ztra * fse3t(ji,jj,jk) * 86400.
               ENDIF   
                
            END DO
         END DO
      END DO
      !
      IF( med_diag%DSED%dgsave ) THEN
           CALL iom_put( "DSED"  ,  zw2d)
           CALL wrk_dealloc( jpi, jpj,    zw2d  )
      ENDIF
      !!
# if defined key_roam

      ! sedimentation of detrital carbon : upstream scheme
      ! --------------------------------------------------
      !
      zwork(:,:,:) = 0.e0        ! initialisation of sinking variable
      ! for detrital carbon sedimentation only - jpdtc
      zwork(:,:,1  ) = 0.e0      ! surface value set to zero
      zwork(:,:,jpk) = 0.e0      ! bottom value  set to zero
      !
      ! tracer flux at w-point: we use -vsed (downward flux)  with simplification : no e1*e2
      DO jk = 2, jpk
         ! AXY (17/07/14): change "0.d0" to "0."
         ! zwork(:,:,jk) = -vsed * max(trn(:,:,jk-1,jpdtc),0.d0) * tmask(:,:,jk-1)
         zwork(:,:,jk) = -vsed * max(trn(:,:,jk-1,jpdtc),0.) * tmask(:,:,jk-1)
         !
         ! AXY (16/01/14): stop sinking in upper 10m to reduce model instability 
         !                 in shallower grid cells
         ! if ( jk .lt. 9 ) zwork(:,:,jk) = 0.e0
      END DO
      !
      ! tracer flux divergence at t-point added to the general trend
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               ztra  = - ( zwork(ji,jj,jk) - zwork(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               tra(ji,jj,jk,jpdtc) = tra(ji,jj,jk,jpdtc) + ztra
            END DO
         END DO
      END DO
      !

# endif

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sed')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

   END SUBROUTINE trc_sed_medusa

   !! ======================================================================
   !! ======================================================================
   !! ======================================================================

   !! AXY (10/02/09)
   !! JPALM -- 31-03-2016 -- Completely change trc_sed_medusa_sbc.
   !!                     -- We now need to read dust file through a namelist.
   !!                     To be able to use time varying dust depositions from
   !!                     -- copy and adapt the PISCES p4z_sbc_ini subroutine
   !!                     -- Only use the dust related part.      
   SUBROUTINE trc_sed_medusa_sbc(kt)

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sed_medusa_sbc  ***
      !!
      !! ** Purpose :   Read and dust namelist and files.
      !!                The interpolation is done in trc_sed through 
      !!                "CALL fld_read( kt, 1, sf_dust )"
      !!
      !! ** Method  :   Read the sbc namelist, and the adapted dust file, if required
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   -- namelist sbc ref and cfg
      !!                -- external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      INTEGER  :: ji, jj, jk, jm, ifpr
      INTEGER  :: ii0, ii1, ij0, ij1
      INTEGER  :: numdust
      INTEGER  :: ierr 
      INTEGER  :: jfld                ! dummy loop arguments
      INTEGER  :: ios                 ! Local integer output status for namelist read
      INTEGER  :: isrow               ! index for ORCA1 starting row
      REAL(wp) :: ztimes_dust
      REAL(wp), DIMENSION(nbtimes) :: zsteps                 ! times records
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zdust
      !
      CHARACTER(len=100)         ::  cn_dir       ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION(1)  ::   slf_d       ! array of namelist informations on the fields to read
      TYPE(FLD_N)                ::   sn_dust     ! informations about the fields to be read
      !
      NAMELIST/nammedsbc/cn_dir, sn_dust, bdustfer 

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sed_medusa_sbc')
      !
      !                            !* set file information
      REWIND( numnatp_ref )        ! Namelist nammedsbc in reference namelist : MEDUSA external sources of Dust
      READ  ( numnatp_ref, nammedsbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nammedsbc in reference namelist', lwp )

      REWIND( numnatp_cfg )        ! Namelist nammedsbc in configuration namelist : MEDUSA external sources of Dust
      READ  ( numnatp_cfg, nammedsbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nammedsbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nammedsbc )

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist : nammedsbc '
          WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
          WRITE(numout,*) '    dust input from the atmosphere bdustfer     = ', bdustfer
      END IF

      ! dust input from the atmosphere
      ! ------------------------------
      IF( bdustfer ) THEN
         !
         IF(lwp) WRITE(numout,*) '    initialize dust input from atmosphere '
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
         !
         !! already allocated in sms_medusa
         !!ALLOCATE( dust(jpi,jpj) )    ! allocation
         !
         slf_d(1) = sn_dust                          ! store namelist information in an array
         !
         ALLOCATE( sf_dust(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trc_sed_medusa_sbc: unable to allocate sf_dust structure' )
         ALLOCATE( sf_dust(1)%fnow(jpi,jpj,1))
         IF( slf_d(1)%ln_tint )     ALLOCATE( sf_dust(1)%fdta(jpi,jpj,1,2) )
         !
         CALL fld_fill( sf_dust, slf_d, cn_dir, 'trc_sed_medusa_sbc', 'Atmospheric dust deposition', 'nammedsed' )
         !
         CALL fld_read( kt, 1, sf_dust )
         dust(:,:) = sf_dust(1)%fnow(:,:,1)
         !
      ELSEIF (lk_oasis) THEN
         dust = Dust_in_cpl
      ELSE
         dust(:,:) = 0.0
      END IF
      !
      zirondep(:,:) = 0.e0     !! Initialisation of deposition variables
      zirondep(:,:) = dust(:,:) * Fe_dust_mratio / xfe_mass * 1.e6 * 86400.      !! mmol-Fe/m2/d
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sed_medusa_sbc')
      !
   END SUBROUTINE trc_sed_medusa_sbc

   !! ======================================================================
   !! ======================================================================
   !! ======================================================================

   !! AXY & JPALM (28/02/17)

   SUBROUTINE trc_sed_medusa_dust( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sed_medusa_dust  ***
      !!
      !! ** Purpose : compute current dust *before* trc_bio_medusa call
      !!
      !! ** Method  : does what it says on the tin
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      

      !! AXY (20/11/14): alter this to report on first MEDUSA call
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_sed_medusa_dust: MEDUSA dust timestep'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
	 IF(lwp) WRITE(numout,*) ' kt =',kt
      ENDIF

      !! AXY (04/11/13): replace this with a call in trc_ini_medusa
      !! AXY (25/02/10)
      !! call routine for populating CCD array if this is the first time-step
      !! IF( kt == nittrc000 ) CALL medusa_ccd( kt )

      !! AXY (04/11/13): replace this with a call in trc_ini_medusa
      !! AXY (26/01/12)
      !! call routine for populating river arrays if this is the first time-step
      !! IF( kt == nittrc000 ) CALL medusa_river( kt )

      !! AXY (10/02/09)
      !! IF( (jnt == 1) .and. (bdustfer) )  CALL trc_sed_medusa_sbc( kt )

      !! JPALM -- 31-03-2016 -- rewrite trc_sed_medusa_sbc.
      !! IF (kt == nittrc000 ) CALL trc_sed_medusa_sbc 

      !! JPALM -- 20-07-2016 -- adapt dust forcing fields reading and conversion
      !!                     To read dust dep in kg-dust/m2/s instead of g-Fe/m2/month 
      !!                     So all forcings and coupling dust dep are in the same SI units
      !!                     and then convert in mmol-Fe/m2/day

      IF( bdustfer ) THEN
            CALL fld_read( kt, 1, sf_dust )
            dust(:,:) = sf_dust(1)%fnow(:,:,1)
      ELSEIF (lk_oasis) THEN
         dust = Dust_in_cpl
      ELSE
         dust(:,:) = 0.0
      ENDIF
      !!
      zirondep(:,:) = 0.e0     !! Initialisation of deposition variables
      zirondep(:,:) = dust(:,:) * Fe_dust_mratio / xfe_mass * 1.e6 * 86400.  !! mmol-Fe/m2/d
      
      !! JPALM -- 20-07-2016 -- Zirondep and zsidep are not used.
      !!                     So comment out the following lines. but keep them
      !!                     as we may want to used them later on
      !!================================================     
      !!
      !! zirondep(:,:,:) = 0.e0     !! Initialisation of deposition variables
      !! zsidep  (:,:)   = 0.e0
      !!
      !! Iron and Si deposition at the surface
      !! -------------------------------------
      !!
      !! DO jj = 1, jpj
      !!    DO ji = 1, jpi
      !!       zirondep(ji,jj,1) = (dustsolub * dust(ji,jj) / (55.85 * rmtss) + 3.e-10 / ryyss) &
      !!       & * rfact2 / fse3t(ji,jj,1)
      !!       zsidep  (ji,jj)   = 8.8 * 0.075 * dust(ji,jj) * rfact2 / &
      !!       & (fse3t(ji,jj,1) * 28.1 * rmtss)
      !!    END DO
      !! END DO

   END SUBROUTINE trc_sed_medusa_dust

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sed_medusa( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sed_medusa: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sed_medusa
#endif 

   !!======================================================================
END MODULE trcsed_medusa
