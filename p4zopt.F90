MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP - PISCES : Compute the light availability in the water column
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!----------------------------------------------------------------------
#if defined  key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE trc            ! tracer variables
   USE oce_trc        ! tracer-ocean share variables
   USE sms_pisces     ! Source Minus Sink of PISCES
   USE iom            ! I/O manager
   USE fldread         !  time interpolation
   USE prtctl_trc      !  print control for debugging


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_opt        ! called in p4zbio.F90 module
   PUBLIC   p4z_opt_init   ! called in trcsms_pisces.F90 module
   PUBLIC   p4z_opt_alloc

   !! * Shared module variables

   LOGICAL  :: ln_varpar   !: boolean for variable PAR fraction
   REAL(wp) :: parlux      !: Fraction of shortwave as PAR
   REAL(wp) :: xparsw                 !: parlux/3
   REAL(wp) :: xsi0r                 !:  1. /rn_si0

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_par      ! structure of input par
   INTEGER , PARAMETER :: nbtimes = 365  !: maximum number of times record in a file
   INTEGER  :: ntimes_par                ! number of time steps in a file
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: par_varsw    !: PAR fraction of shortwave

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: enano, ediat   !: PAR for phyto, nano and diat 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: etot_ndcy      !: PAR over 24h in case of diurnal cycle
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: emoy           !: averaged PAR in the mixed layer
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ekb, ekg, ekr  !: wavelength (Red-Green-Blue)

   INTEGER  ::   nksrp   ! levels below which the light cannot penetrate ( depth larger than 391 m)

   REAL(wp), DIMENSION(3,61), PUBLIC ::   xkrgb   !: tabulated attenuation coefficients for RGB absorption
   
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_opt( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt  ***
      !!
      !! ** Purpose :   Compute the light availability in the water column
      !!              depending on the depth and the chlorophyll concentration
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   irgb
      REAL(wp) ::   zchl
      REAL(wp) ::   zc0 , zc1 , zc2, zc3, z1_dep
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zdepmoy, zetmp1, zetmp2, zetmp3, zetmp4
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zqsr100, zqsr_corr
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zpar, ze0, ze1, ze2, ze3
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_opt')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      zdepmoy, zetmp1, zetmp2, zetmp3, zetmp4 )
      CALL wrk_alloc( jpi, jpj,      zqsr100, zqsr_corr )
      CALL wrk_alloc( jpi, jpj, jpk, zpar, ze0, ze1, ze2, ze3 )

      IF( knt == 1 .AND. ln_varpar ) CALL p4z_opt_sbc( kt )

      !     Initialisation of variables used to compute PAR
      !     -----------------------------------------------
      ze1(:,:,:) = 0._wp
      ze2(:,:,:) = 0._wp
      ze3(:,:,:) = 0._wp
      !                                        !* attenuation coef. function of Chlorophyll and wavelength (Red-Green-Blue)
      DO jk = 1, jpkm1                         !  --------------------------------------------------------
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zchl = ( trb(ji,jj,jk,jpnch) + trb(ji,jj,jk,jpdch) + rtrn ) * 1.e6
               zchl = MIN(  10. , MAX( 0.05, zchl )  )
               irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )
               !                                                         
               ekb(ji,jj,jk) = xkrgb(1,irgb) * fse3t(ji,jj,jk)
               ekg(ji,jj,jk) = xkrgb(2,irgb) * fse3t(ji,jj,jk)
               ekr(ji,jj,jk) = xkrgb(3,irgb) * fse3t(ji,jj,jk)
            END DO
         END DO
      END DO
      !                                        !* Photosynthetically Available Radiation (PAR)
      !                                        !  --------------------------------------
      IF( l_trcdm2dc ) THEN                     !  diurnal cycle
         !
         zqsr_corr(:,:) = qsr_mean(:,:) / ( 1. - fr_i(:,:) + rtrn )
         !
         CALL p4z_opt_par( kt, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 ) 
         !
         DO jk = 1, nksrp      
            etot_ndcy(:,:,jk) =        ze1(:,:,jk) +        ze2(:,:,jk) +       ze3(:,:,jk)
            enano    (:,:,jk) =  2.1 * ze1(:,:,jk) + 0.42 * ze2(:,:,jk) + 0.4 * ze3(:,:,jk)
            ediat    (:,:,jk) =  1.6 * ze1(:,:,jk) + 0.69 * ze2(:,:,jk) + 0.7 * ze3(:,:,jk)
         END DO
         !
         zqsr_corr(:,:) = qsr(:,:) / ( 1. - fr_i(:,:) + rtrn )
         !
         CALL p4z_opt_par( kt, zqsr_corr, ze1, ze2, ze3 ) 
         !
         DO jk = 1, nksrp      
            etot(:,:,jk) =  ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk)
         END DO
         !
      ELSE
         !
         zqsr_corr(:,:) = qsr(:,:) / ( 1. - fr_i(:,:) + rtrn )
         !
         CALL p4z_opt_par( kt, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 ) 
         !
         DO jk = 1, nksrp      
            etot (:,:,jk) =        ze1(:,:,jk) +        ze2(:,:,jk) +       ze3(:,:,jk)
            enano(:,:,jk) =  2.1 * ze1(:,:,jk) + 0.42 * ze2(:,:,jk) + 0.4 * ze3(:,:,jk)
            ediat(:,:,jk) =  1.6 * ze1(:,:,jk) + 0.69 * ze2(:,:,jk) + 0.7 * ze3(:,:,jk)
         END DO
         etot_ndcy(:,:,:) =  etot(:,:,:) 
      ENDIF


      IF( ln_qsr_bio ) THEN                    !* heat flux accros w-level (used in the dynamics)
         !                                     !  ------------------------
         CALL p4z_opt_par( kt, qsr, ze1, ze2, ze3, pe0=ze0 )
         !
         etot3(:,:,1) =  qsr(:,:) * tmask(:,:,1)
         DO jk = 2, nksrp + 1
            etot3(:,:,jk) =  ( ze0(:,:,jk) + ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk) ) * tmask(:,:,jk)
         END DO
         !                                     !  ------------------------
      ENDIF
      !                                        !* Euphotic depth and level
      neln(:,:) = 1                            !  ------------------------
      heup(:,:) = 300.

      DO jk = 2, nksrp
         DO jj = 1, jpj
           DO ji = 1, jpi
              IF( etot_ndcy(ji,jj,jk) * tmask(ji,jj,jk) >= zqsr100(ji,jj) )  THEN
                 neln(ji,jj) = jk+1                    ! Euphotic level : 1rst T-level strictly below Euphotic layer
                 !                                     ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
                 heup(ji,jj) = fsdepw(ji,jj,jk+1)      ! Euphotic layer depth
              ENDIF
           END DO
        END DO
      END DO
      !
      heup(:,:) = MIN( 300., heup(:,:) )
      !                                        !* mean light over the mixed layer
      zdepmoy(:,:)   = 0.e0                    !  -------------------------------
      zetmp1 (:,:)   = 0.e0
      zetmp2 (:,:)   = 0.e0
      zetmp3 (:,:)   = 0.e0
      zetmp4 (:,:)   = 0.e0

      DO jk = 1, nksrp
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( fsdepw(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  zetmp1 (ji,jj) = zetmp1 (ji,jj) + etot     (ji,jj,jk) * fse3t(ji,jj,jk) ! remineralisation
                  zetmp2 (ji,jj) = zetmp2 (ji,jj) + etot_ndcy(ji,jj,jk) * fse3t(ji,jj,jk) ! production
                  zetmp3 (ji,jj) = zetmp3 (ji,jj) + enano    (ji,jj,jk) * fse3t(ji,jj,jk) ! production
                  zetmp4 (ji,jj) = zetmp4 (ji,jj) + ediat    (ji,jj,jk) * fse3t(ji,jj,jk) ! production
                  zdepmoy(ji,jj) = zdepmoy(ji,jj) + fse3t(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO
      !
      emoy(:,:,:) = etot(:,:,:)       ! remineralisation
      zpar(:,:,:) = etot_ndcy(:,:,:)  ! diagnostic : PAR with no diurnal cycle 
      !
      DO jk = 1, nksrp
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( fsdepw(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                  emoy (ji,jj,jk) = zetmp1(ji,jj) * z1_dep
                  zpar (ji,jj,jk) = zetmp2(ji,jj) * z1_dep
                  enano(ji,jj,jk) = zetmp3(ji,jj) * z1_dep
                  ediat(ji,jj,jk) = zetmp4(ji,jj) * z1_dep
               ENDIF
            END DO
         END DO
      END DO
      !
      IF( lk_iomput ) THEN
        IF( knt == nrdttrc ) THEN
           IF( iom_use( "Heup"  ) ) CALL iom_put( "Heup" , heup(:,:  ) * tmask(:,:,1) )  ! euphotic layer deptht
           IF( iom_use( "PARDM" ) ) CALL iom_put( "PARDM", zpar(:,:,:) * tmask(:,:,:) )  ! Photosynthetically Available Radiation
           IF( iom_use( "PAR"   ) ) CALL iom_put( "PAR"  , emoy(:,:,:) * tmask(:,:,:) )  ! Photosynthetically Available Radiation
        ENDIF
      ELSE
         IF( ln_diatrc ) THEN        ! save output diagnostics
            trc2d(:,:,  jp_pcs0_2d + 10) = heup(:,:  ) * tmask(:,:,1)
            trc3d(:,:,:,jp_pcs0_3d + 3)  = etot(:,:,:) * tmask(:,:,:)
         ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      zdepmoy, zetmp1, zetmp2, zetmp3, zetmp4 )
      CALL wrk_dealloc( jpi, jpj,      zqsr100, zqsr_corr )
      CALL wrk_dealloc( jpi, jpj, jpk, zpar,  ze0, ze1, ze2, ze3 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_opt')
      !
   END SUBROUTINE p4z_opt

   SUBROUTINE p4z_opt_par( kt, pqsr, pe1, pe2, pe3, pe0, pqsr100 ) 
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_opt_par  ***
      !!
      !! ** purpose :   compute PAR of each wavelength (Red-Green-Blue)
      !!                for a given shortwave radiation
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT(in)                                        ::  kt            !   ocean time-step
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in)               ::  pqsr          !   shortwave
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)            ::  pe1 , pe2 , pe3   !  PAR ( R-G-B)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout), OPTIONAL  ::  pe0 
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(out)  , OPTIONAL  ::  pqsr100  
      !! * local variables
      INTEGER    ::   ji, jj, jk     ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj)     ::  zqsr          !   shortwave
      !!----------------------------------------------------------------------

      !  Real shortwave
      IF( ln_varpar ) THEN  ;  zqsr(:,:) = par_varsw(:,:) * pqsr(:,:)
      ELSE                  ;  zqsr(:,:) = xparsw         * pqsr(:,:)
      ENDIF

      !  Light at the euphotic depth 
      IF( PRESENT( pqsr100 ) )  pqsr100(:,:) = 0.01 * 3. * zqsr(:,:)
      !
      IF( PRESENT( pe0 ) ) THEN     !  W-level
         !
         pe0(:,:,1) = pqsr(:,:) - 3. * zqsr(:,:)    !   ( 1 - 3 * alpha ) * q
         pe1(:,:,1) = zqsr(:,:)         
         pe2(:,:,1) = zqsr(:,:)
         pe3(:,:,1) = zqsr(:,:)
         !
         DO jk = 2, nksrp + 1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  pe0(ji,jj,jk) = pe0(ji,jj,jk-1) * EXP( -fse3t(ji,jj,jk-1) * xsi0r )
                  pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -ekb(ji,jj,jk-1 ) )
                  pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -ekg(ji,jj,jk-1 ) )
                  pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -ekr(ji,jj,jk-1 ) )
               END DO
              !
            END DO
            !
         END DO
        !
      ELSE   ! T- level
        !
        pe1(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekb(:,:,1) )
        pe2(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekg(:,:,1) )
        pe3(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekr(:,:,1) )
        !
        DO jk = 2, nksrp      
!CDIR NOVERRCHK
           DO jj = 1, jpj
!CDIR NOVERRCHK
              DO ji = 1, jpi
                 pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -0.5 * ( ekb(ji,jj,jk-1) + ekb(ji,jj,jk) ) )
                 pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -0.5 * ( ekg(ji,jj,jk-1) + ekg(ji,jj,jk) ) )
                 pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -0.5 * ( ekr(ji,jj,jk-1) + ekr(ji,jj,jk) ) )
              END DO
           END DO
        END DO    
        !
      ENDIF
      ! 
   END SUBROUTINE p4z_opt_par


   SUBROUTINE p4z_opt_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_opt_sbc  ***
      !!
      !! ** purpose :   read and interpolate the variable PAR fraction
      !!                of shortwave radiation
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER ,                INTENT(in) ::   kt     ! ocean time step

      !! * local declarations
      INTEGER  :: ji,jj
      REAL(wp) :: zcoef
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_optsbc')
      !
      ! Compute par_varsw at nit000 or only if there is more than 1 time record in par coefficient file
      IF( ln_varpar ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_par > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_par )
            par_varsw(:,:) = ( sf_par(1)%fnow(:,:,1) ) / 3.0
         ENDIF
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_optsbc')
      !
   END SUBROUTINE p4z_opt_sbc

   SUBROUTINE p4z_opt_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!                and of the percentage of PAR in Shortwave
      !!
      !! ** Input   :   external ascii and netcdf files
      !!----------------------------------------------------------------------
      !
      INTEGER :: numpar
      INTEGER :: ierr
      INTEGER :: ios                 ! Local integer output status for namelist read
      REAL(wp), DIMENSION(nbtimes) :: zsteps                 ! times records
      !
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N) ::   sn_par                ! informations about the fields to be read
      !
      NAMELIST/nampisopt/cn_dir, sn_par, ln_varpar, parlux

      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('p4z_opt_init')

      REWIND( numnatp_ref )              ! Namelist nampisopt in reference namelist : Pisces attenuation coef. and PAR
      READ  ( numnatp_ref, nampisopt, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisopt in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisopt in configuration namelist : Pisces attenuation coef. and PAR
      READ  ( numnatp_cfg, nampisopt, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisopt in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisopt )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampisopt '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    PAR as a variable fraction of SW     ln_varpar      = ', ln_varpar
         WRITE(numout,*) '    Default value for the PAR fraction   parlux         = ', parlux
      ENDIF
      !
      xparsw = parlux / 3.0
      xsi0r  = 1.e0 / rn_si0
      !
      ! Variable PAR at the surface of the ocean
      ! ----------------------------------------
      IF( ln_varpar ) THEN
         IF(lwp) WRITE(numout,*) '    initialize variable par fraction '
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !
         ALLOCATE( par_varsw(jpi,jpj) )
         !
         ALLOCATE( sf_par(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_opt_init: unable to allocate sf_par structure' )
         !
         CALL fld_fill( sf_par, (/ sn_par /), cn_dir, 'p4z_opt_init', 'Variable PAR fraction ', 'nampisopt' )
                                   ALLOCATE( sf_par(1)%fnow(jpi,jpj,1)   )
         IF( sn_par%ln_tint )      ALLOCATE( sf_par(1)%fdta(jpi,jpj,1,2) )

         CALL iom_open (  TRIM( sn_par%clname ) , numpar )
         CALL iom_gettime( numpar, zsteps, kntime=ntimes_par)  ! get number of record in file
      ENDIF
      !
      CALL trc_oce_rgb( xkrgb )                  ! tabulated attenuation coefficients
      nksrp = trc_oce_ext_lev( r_si2, 0.33e2 )     ! max level of light extinction (Blue Chl=0.01)
      !
      IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksrp, ' ref depth = ', gdepw_1d(nksrp+1), ' m'
      !
                         ekr      (:,:,:) = 0._wp
                         ekb      (:,:,:) = 0._wp
                         ekg      (:,:,:) = 0._wp
                         etot     (:,:,:) = 0._wp
                         etot_ndcy(:,:,:) = 0._wp
                         enano    (:,:,:) = 0._wp
                         ediat    (:,:,:) = 0._wp
      IF( ln_qsr_bio )   etot3    (:,:,:) = 0._wp
      ! 
      IF( nn_timing == 1 )  CALL timing_stop('p4z_opt_init')
      !
   END SUBROUTINE p4z_opt_init


   INTEGER FUNCTION p4z_opt_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( ekb(jpi,jpj,jpk)      , ekr(jpi,jpj,jpk), ekg(jpi,jpj,jpk),   &
        &       enano(jpi,jpj,jpk)    , ediat(jpi,jpj,jpk), &
        &       etot_ndcy(jpi,jpj,jpk), emoy (jpi,jpj,jpk), STAT=p4z_opt_alloc ) 
         !
      IF( p4z_opt_alloc /= 0 ) CALL ctl_warn('p4z_opt_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_opt_alloc

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE p4z_opt                   ! Empty routine
   END SUBROUTINE p4z_opt
#endif 

   !!======================================================================
END MODULE p4zopt
