MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wscal    !: Calcite and BSi sinking speeds

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes
#if ! defined key_kriest
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes
#endif

   INTEGER  :: ik100

#if  defined key_kriest
   REAL(wp) ::  xkr_sfact    !: Sinking factor
   REAL(wp) ::  xkr_stick    !: Stickiness
   REAL(wp) ::  xkr_nnano    !: Nbr of cell in nano size class
   REAL(wp) ::  xkr_ndiat    !: Nbr of cell in diatoms size class
   REAL(wp) ::  xkr_nmicro   !: Nbr of cell in microzoo size class
   REAL(wp) ::  xkr_nmeso    !: Nbr of cell in mesozoo  size class
   REAL(wp) ::  xkr_naggr    !: Nbr of cell in aggregates  size class

   REAL(wp) ::  xkr_frac 

   REAL(wp), PUBLIC ::  xkr_dnano       !: Size of particles in nano pool
   REAL(wp), PUBLIC ::  xkr_ddiat       !: Size of particles in diatoms pool
   REAL(wp), PUBLIC ::  xkr_dmicro      !: Size of particles in microzoo pool
   REAL(wp), PUBLIC ::  xkr_dmeso       !: Size of particles in mesozoo pool
   REAL(wp), PUBLIC ::  xkr_daggr       !: Size of particles in aggregates pool
   REAL(wp), PUBLIC ::  xkr_wsbio_min   !: min vertical particle speed
   REAL(wp), PUBLIC ::  xkr_wsbio_max   !: max vertical particle speed

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   xnumm   !:  maximum number of particles in aggregates
#endif

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if ! defined key_kriest
   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk, jit
      INTEGER  ::   iiter1, iiter2
      REAL(wp) ::   zagg1, zagg2, zagg3, zagg4
      REAL(wp) ::   zagg , zaggfe, zaggdoc, zaggdoc2, zaggdoc3
      REAL(wp) ::   zfact, zwsmax, zmax, zstep
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zw3d
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., fsdepw(ji,jj,jk+1) - zmax ) / 5000._wp
               wsbio4(ji,jj,jk) = wsbio2 + ( 200.- wsbio2 ) * zfact
            END DO
         END DO
      END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio
      wscal (:,:,:) = wsbio4(:,:,:)
      !
      ! OA This is (I hope) a temporary solution for the problem that may 
      ! OA arise in specific situation where the CFL criterion is broken 
      ! OA for vertical sedimentation of particles. To avoid this, a time
      ! OA splitting algorithm has been coded. A specific maximum
      ! OA iteration number is provided and may be specified in the namelist 
      ! OA This is to avoid very large iteration number when explicit free
      ! OA surface is used (for instance). When niter?max is set to 1, 
      ! OA this computation is skipped. The crude old threshold method is 
      ! OA then applied. This also happens when niter exceeds nitermax.
      IF( MAX( niter1max, niter2max ) == 1 ) THEN
        iiter1 = 1
        iiter2 = 1
      ELSE
        iiter1 = 1
        iiter2 = 1
        DO jk = 1, jpkm1
          DO jj = 1, jpj
             DO ji = 1, jpi
                IF( tmask(ji,jj,jk) == 1) THEN
                   zwsmax =  0.5 * fse3t(ji,jj,jk) / xstep
                   iiter1 =  MAX( iiter1, INT( wsbio3(ji,jj,jk) / zwsmax ) )
                   iiter2 =  MAX( iiter2, INT( wsbio4(ji,jj,jk) / zwsmax ) )
                ENDIF
             END DO
          END DO
        END DO
        IF( lk_mpp ) THEN
           CALL mpp_max( iiter1 )
           CALL mpp_max( iiter2 )
        ENDIF
        iiter1 = MIN( iiter1, niter1max )
        iiter2 = MIN( iiter2, niter2max )
      ENDIF

      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) == 1 ) THEN
                 zwsmax = 0.5 * fse3t(ji,jj,jk) / xstep
                 wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax * FLOAT( iiter1 ) )
                 wsbio4(ji,jj,jk) = MIN( wsbio4(ji,jj,jk), zwsmax * FLOAT( iiter2 ) )
               ENDIF
            END DO
         END DO
      END DO

      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0
      sinkfer2(:,:,:) = 0.e0

      !   Compute the sedimentation term using p4zsink2 for all the sinking particles
      !   -----------------------------------------------------
      DO jit = 1, iiter1
        CALL p4z_sink2( wsbio3, sinking , jppoc, iiter1 )
        CALL p4z_sink2( wsbio3, sinkfer , jpsfe, iiter1 )
      END DO

      DO jit = 1, iiter2
        CALL p4z_sink2( wsbio4, sinking2, jpgoc, iiter2 )
        CALL p4z_sink2( wsbio4, sinkfer2, jpbfe, iiter2 )
        CALL p4z_sink2( wsbio4, sinksil , jpgsi, iiter2 )
        CALL p4z_sink2( wscal , sinkcal , jpcal, iiter2 )
      END DO

      !  Exchange between organic matter compartments due to coagulation/disaggregation
      !  ---------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               zstep = xstep 
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               zfact = zstep * xdiss(ji,jj,jk)
               !  Part I : Coagulation dependent on turbulence
               zagg1 = 25.9  * zfact * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jppoc)
               zagg2 = 4452. * zfact * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpgoc)

               ! Part II : Differential settling

               !  Aggregation of small into large particles
               zagg3 =  47.1 * zstep * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpgoc)
               zagg4 =  3.3  * zstep * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jppoc)

               zagg   = zagg1 + zagg2 + zagg3 + zagg4
               zaggfe = zagg * trb(ji,jj,jk,jpsfe) / ( trb(ji,jj,jk,jppoc) + rtrn )

               ! Aggregation of DOC to POC : 
               ! 1st term is shear aggregation of DOC-DOC
               ! 2nd term is shear aggregation of DOC-POC
               ! 3rd term is differential settling of DOC-POC
               zaggdoc  = ( ( 0.369 * 0.3 * trb(ji,jj,jk,jpdoc) + 102.4 * trb(ji,jj,jk,jppoc) ) * zfact       &
               &            + 2.4 * zstep * trb(ji,jj,jk,jppoc) ) * 0.3 * trb(ji,jj,jk,jpdoc)
               ! transfer of DOC to GOC : 
               ! 1st term is shear aggregation
               ! 2nd term is differential settling 
               zaggdoc2 = ( 3.53E3 * zfact + 0.1 * zstep ) * trb(ji,jj,jk,jpgoc) * 0.3 * trb(ji,jj,jk,jpdoc)
               ! tranfer of DOC to POC due to brownian motion
               zaggdoc3 =  ( 5095. * trb(ji,jj,jk,jppoc) + 114. * 0.3 * trb(ji,jj,jk,jpdoc) ) *zstep * 0.3 * trb(ji,jj,jk,jpdoc)

               !  Update the trends
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zagg + zaggdoc + zaggdoc3
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zagg + zaggdoc2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zaggfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zaggfe
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc - zaggdoc2 - zaggdoc3
               !
            END DO
         END DO
      END DO


     ! Total carbon export per year
     IF( iom_use( "tcexp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
        &   t_oce_co2_exp = glob_sum( ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * e1e2t(:,:) * tmask(:,:,1) )
     !
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          CALL wrk_alloc( jpi, jpj,      zw2d )
          CALL wrk_alloc( jpi, jpj, jpk, zw3d )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPFE100" ) )  THEN
              zw2d(:,:) = ( sinkfer(:,:,ik100) + sinkfer2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of iron at 100m
              CALL iom_put( "EPFE100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSI100" ) )  THEN
              zw2d(:,:) =  sinksil(:,:,ik100) * zfact * tmask(:,:,1) ! Export of bigenic silica at 100m
              CALL iom_put( "EPSI100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPC" ) )  THEN
              zw3d(:,:,:) = ( sinking(:,:,:) + sinking2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPC"  , zw3d )
          ENDIF
          IF( iom_use( "EXPFE" ) )  THEN
              zw3d(:,:,:) = ( sinkfer(:,:,:) + sinkfer2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of iron 
              CALL iom_put( "EXPFE"  , zw3d )
          ENDIF
          IF( iom_use( "EXPCAL" ) )  THEN
              zw3d(:,:,:) = sinkcal(:,:,:) * zfact * tmask(:,:,:) ! Export of calcite 
              CALL iom_put( "EXPCAL"  , zw3d )
          ENDIF
          IF( iom_use( "EXPSI" ) )  THEN
              zw3d(:,:,:) = sinksil(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPSI"  , zw3d )
          ENDIF
          IF( iom_use( "tcexp" ) )  CALL iom_put( "tcexp" , t_oce_co2_exp * zfact )   ! molC/s
          ! 
          CALL wrk_dealloc( jpi, jpj,      zw2d )
          CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
        ENDIF
      ELSE
         IF( ln_diatrc ) THEN
            zfact = 1.e3 * rfact2r
            trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 5) = sinking2(:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 6) = sinkfer (:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 7) = sinkfer2(:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 8) = sinksil (:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 9) = sinkcal (:,:,ik100) * zfact * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink

   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER :: jk

      ik100 = 10        !  last level where depth less than 100 m
      DO jk = jpkm1, 1, -1
         IF( gdept_1d(jk) > 100. )  ik100 = jk - 1
      END DO
      IF (lwp) WRITE(numout,*)
      IF (lwp) WRITE(numout,*) ' Level corresponding to 100m depth ',  ik100 + 1
      IF (lwp) WRITE(numout,*)
      !
      t_oce_co2_exp = 0._wp
      !
   END SUBROUTINE p4z_sink_init

#else
   !!----------------------------------------------------------------------
   !!   'Kriest sinking parameterisation'        key_kriest          ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to
      !!              gravitational sinking - Kriest parameterization
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER  :: ji, jj, jk, jit, niter1, niter2
      REAL(wp) :: zagg1, zagg2, zagg3, zagg4, zagg5, zfract, zaggsi, zaggsh
      REAL(wp) :: zagg , zaggdoc, zaggdoc1, znumdoc
      REAL(wp) :: znum , zeps, zfm, zgm, zsm
      REAL(wp) :: zdiv , zdiv1, zdiv2, zdiv3, zdiv4, zdiv5
      REAL(wp) :: zval1, zval2, zval3, zval4
      REAL(wp) :: zfact
      INTEGER  :: ik1
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: znum3d 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zw3d
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      CALL wrk_alloc( jpi, jpj, jpk, znum3d )
      !
      !     Initialisation of variables used to compute Sinking Speed
      !     ---------------------------------------------------------

      znum3d(:,:,:) = 0.e0
      zval1 = 1. + xkr_zeta
      zval2 = 1. + xkr_zeta + xkr_eta
      zval3 = 1. + xkr_eta

      !     Computation of the vertical sinking speed : Kriest et Evans, 2000
      !     -----------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) /= 0.e0 ) THEN
                  znum = trb(ji,jj,jk,jppoc) / ( trb(ji,jj,jk,jpnum) + rtrn ) / xkr_massp
                  ! -------------- To avoid sinking speed over 50 m/day -------
                  znum  = MIN( xnumm(jk), znum )
                  znum  = MAX( 1.1      , znum )
                  znum3d(ji,jj,jk) = znum
                  !------------------------------------------------------------
                  zeps  = ( zval1 * znum - 1. )/ ( znum - 1. )
                  zfm   = xkr_frac**( 1. - zeps )
                  zgm   = xkr_frac**( zval1 - zeps )
                  zdiv  = MAX( 1.e-4, ABS( zeps - zval2 ) ) * SIGN( 1., ( zeps - zval2 ) )
                  zdiv1 = zeps - zval3
                  wsbio3(ji,jj,jk) = xkr_wsbio_min * ( zeps - zval1 ) / zdiv    &
                     &             - xkr_wsbio_max *   zgm * xkr_eta  / zdiv
                  wsbio4(ji,jj,jk) = xkr_wsbio_min *   ( zeps-1. )    / zdiv1   &
                     &             - xkr_wsbio_max *   zfm * xkr_eta  / zdiv1
                  IF( znum == 1.1)   wsbio3(ji,jj,jk) = wsbio4(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO

      wscal(:,:,:) = MAX( wsbio3(:,:,:), 30._wp )

      !   INITIALIZE TO ZERO ALL THE SINKING ARRAYS
      !   -----------------------------------------

      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0

     !   Compute the sedimentation term using p4zsink2 for all the sinking particles
     !   -----------------------------------------------------

      niter1 = niter1max
      niter2 = niter2max

      DO jit = 1, niter1
        CALL p4z_sink2( wsbio3, sinking , jppoc, niter1 )
        CALL p4z_sink2( wsbio3, sinkfer , jpsfe, niter1 )
        CALL p4z_sink2( wscal , sinksil , jpgsi, niter1 )
        CALL p4z_sink2( wscal , sinkcal , jpcal, niter1 )
      END DO

      DO jit = 1, niter2
        CALL p4z_sink2( wsbio4, sinking2, jpnum, niter2 )
      END DO

     !  Exchange between organic matter compartments due to coagulation/disaggregation
     !  ---------------------------------------------------

      zval1 = 1. + xkr_zeta
      zval2 = 1. + xkr_eta
      zval3 = 3. + xkr_eta
      zval4 = 4. + xkr_eta

      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1,jpi
               IF( tmask(ji,jj,jk) /= 0.e0 ) THEN

                  znum = trb(ji,jj,jk,jppoc)/(trb(ji,jj,jk,jpnum)+rtrn) / xkr_massp
                  !-------------- To avoid sinking speed over 50 m/day -------
                  znum  = min(xnumm(jk),znum)
                  znum  = MAX( 1.1,znum)
                  !------------------------------------------------------------
                  zeps  = ( zval1 * znum - 1.) / ( znum - 1.)
                  zdiv  = MAX( 1.e-4, ABS( zeps - zval3) ) * SIGN( 1., zeps - zval3 )
                  zdiv1 = MAX( 1.e-4, ABS( zeps - 4.   ) ) * SIGN( 1., zeps - 4.    )
                  zdiv2 = zeps - 2.
                  zdiv3 = zeps - 3.
                  zdiv4 = zeps - zval2
                  zdiv5 = 2.* zeps - zval4
                  zfm   = xkr_frac**( 1.- zeps )
                  zsm   = xkr_frac**xkr_eta

                  !    Part I : Coagulation dependant on turbulence
                  !    ----------------------------------------------

                  zagg1 =  0.163 * trb(ji,jj,jk,jpnum)**2               &
                     &            * 2.*( (zfm-1.)*(zfm*xkr_mass_max**3-xkr_mass_min**3)    &
                     &            * (zeps-1)/zdiv1 + 3.*(zfm*xkr_mass_max-xkr_mass_min)    &
                     &            * (zfm*xkr_mass_max**2-xkr_mass_min**2)                  &
                     &            * (zeps-1.)**2/(zdiv2*zdiv3)) 
                  zagg2 =  2*0.163*trb(ji,jj,jk,jpnum)**2*zfm*                       &
                     &                   ((xkr_mass_max**3+3.*(xkr_mass_max**2          &
                     &                    *xkr_mass_min*(zeps-1.)/zdiv2                 &
                     &                    +xkr_mass_max*xkr_mass_min**2*(zeps-1.)/zdiv3)    &
                     &                    +xkr_mass_min**3*(zeps-1)/zdiv1)                  &
                     &                    -zfm*xkr_mass_max**3*(1.+3.*((zeps-1.)/           &
                     &                    (zeps-2.)+(zeps-1.)/zdiv3)+(zeps-1.)/zdiv1))    

                  zagg3 =  0.163*trb(ji,jj,jk,jpnum)**2*zfm**2*8. * xkr_mass_max**3  
                  
                 !    Aggregation of small into large particles
                 !    Part II : Differential settling
                 !    ----------------------------------------------

                  zagg4 =  2.*3.141*0.125*trb(ji,jj,jk,jpnum)**2*                       &
                     &                 xkr_wsbio_min*(zeps-1.)**2                         &
                     &                 *(xkr_mass_min**2*((1.-zsm*zfm)/(zdiv3*zdiv4)      &
                     &                 -(1.-zfm)/(zdiv*(zeps-1.)))-                       &
                     &                 ((zfm*zfm*xkr_mass_max**2*zsm-xkr_mass_min**2)     &
                     &                 *xkr_eta)/(zdiv*zdiv3*zdiv5) )   

                  zagg5 =   2.*3.141*0.125*trb(ji,jj,jk,jpnum)**2                         &
                     &                 *(zeps-1.)*zfm*xkr_wsbio_min                        &
                     &                 *(zsm*(xkr_mass_min**2-zfm*xkr_mass_max**2)         &
                     &                 /zdiv3-(xkr_mass_min**2-zfm*zsm*xkr_mass_max**2)    &
                     &                 /zdiv)  

                  !
                  !     Fractionnation by swimming organisms
                  !     ------------------------------------

                  zfract = 2.*3.141*0.125*trb(ji,jj,jk,jpmes)*12./0.12/0.06**3*trb(ji,jj,jk,jpnum)  &
                    &      * (0.01/xkr_mass_min)**(1.-zeps)*0.1**2  &
                    &      * 10000.*xstep

                  !     Aggregation of DOC to small particles
                  !     --------------------------------------

                  zaggdoc = 0.83 * trb(ji,jj,jk,jpdoc) * xstep * xdiss(ji,jj,jk) * trb(ji,jj,jk,jpdoc)   &
                     &        + 0.005 * 231. * trb(ji,jj,jk,jpdoc) * xstep * trb(ji,jj,jk,jpdoc)
                  zaggdoc1 = 271. * trb(ji,jj,jk,jppoc) * xstep * xdiss(ji,jj,jk) * trb(ji,jj,jk,jpdoc)  &
                     &  + 0.02 * 16706. * trb(ji,jj,jk,jppoc) * xstep * trb(ji,jj,jk,jpdoc)

# if defined key_degrad
                   zagg1   = zagg1   * facvol(ji,jj,jk)                 
                   zagg2   = zagg2   * facvol(ji,jj,jk)                 
                   zagg3   = zagg3   * facvol(ji,jj,jk)                 
                   zagg4   = zagg4   * facvol(ji,jj,jk)                 
                   zagg5   = zagg5   * facvol(ji,jj,jk)                 
                   zaggdoc = zaggdoc * facvol(ji,jj,jk)                 
                   zaggdoc1 = zaggdoc1 * facvol(ji,jj,jk)
# endif
                  zaggsh = ( zagg1 + zagg2 + zagg3 ) * rfact2 * xdiss(ji,jj,jk) / 1000.
                  zaggsi = ( zagg4 + zagg5 ) * xstep / 10.
                  zagg = 0.5 * xkr_stick * ( zaggsh + zaggsi )
                  !
                  znumdoc = trb(ji,jj,jk,jpnum) / ( trb(ji,jj,jk,jppoc) + rtrn )
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zaggdoc + zaggdoc1
                  tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + zfract + zaggdoc / xkr_massp - zagg
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc - zaggdoc1

               ENDIF
            END DO
         END DO
      END DO

     ! Total primary production per year
     t_oce_co2_exp = t_oce_co2_exp + glob_sum( ( sinking(:,:,ik100) * e1e2t(:,:) * tmask(:,:,1) )
     !
     IF( lk_iomput ) THEN
        IF( knt == nrdttrc ) THEN
          CALL wrk_alloc( jpi, jpj,      zw2d )
          CALL wrk_alloc( jpi, jpj, jpk, zw3d )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = sinking(:,:,ik100) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPN100" ) )  THEN
              zw2d(:,:) = sinking2(:,:,ik100) * zfact * tmask(:,:,1) ! Export of number of aggregates ?
              CALL iom_put( "EPN100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSI100" ) )  THEN
              zw2d(:,:) = sinksil(:,:,ik100) * zfact * tmask(:,:,1) ! Export of bigenic silica at 100m
              CALL iom_put( "EPSI100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPC" ) )  THEN
              zw3d(:,:,:) = sinking(:,:,:) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPC"  , zw3d )
          ENDIF
          IF( iom_use( "EXPN" ) )  THEN
              zw3d(:,:,:) = sinking(:,:,:) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPN"  , zw3d )
          ENDIF
          IF( iom_use( "EXPCAL" ) )  THEN
              zw3d(:,:,:) = sinkcal(:,:,:) * zfact * tmask(:,:,:) ! Export of calcite 
              CALL iom_put( "EXPCAL"  , zw3d )
          ENDIF
          IF( iom_use( "EXPSI" ) )  THEN
              zw3d(:,:,:) = sinksil(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPSI"  , zw3d )
          ENDIF
          IF( iom_use( "XNUM" ) )  THEN
              zw3d(:,:,:) =  znum3d(:,:,:) * tmask(:,:,:) !  Number of particles on aggregats
              CALL iom_put( "XNUM"  , zw3d )
          ENDIF
          IF( iom_use( "WSC" ) )  THEN
              zw3d(:,:,:) = wsbio3(:,:,:) * tmask(:,:,:) ! Sinking speed of carbon particles
              CALL iom_put( "WSC"  , zw3d )
          ENDIF
          IF( iom_use( "WSN" ) )  THEN
              zw3d(:,:,:) = wsbio4(:,:,:) * tmask(:,:,:) ! Sinking speed of particles number
              CALL iom_put( "WSN"  , zw3d )
          ENDIF
          !
          CALL wrk_dealloc( jpi, jpj,      zw2d )
          CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
      ELSE
         IF( ln_diatrc ) THEN
            zfact = 1.e3 * rfact2r
            trc2d(:,:  ,jp_pcs0_2d + 4)  = sinking (:,:,ik100)  * zfact * tmask(:,:,1)
            trc2d(:,:  ,jp_pcs0_2d + 5)  = sinking2(:,:,ik100)  * zfact * tmask(:,:,1)
            trc2d(:,:  ,jp_pcs0_2d + 6)  = sinkfer (:,:,ik100)  * zfact * tmask(:,:,1)
            trc2d(:,:  ,jp_pcs0_2d + 7)  = sinksil (:,:,ik100)  * zfact * tmask(:,:,1)
            trc2d(:,:  ,jp_pcs0_2d + 8)  = sinkcal (:,:,ik100)  * zfact * tmask(:,:,1)
            trc3d(:,:,:,jp_pcs0_3d + 11) = sinking (:,:,:)      * zfact * tmask(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 12) = sinking2(:,:,:)      * zfact * tmask(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 13) = sinksil (:,:,:)      * zfact * tmask(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 14) = sinkcal (:,:,:)      * zfact * tmask(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 15) = znum3d  (:,:,:)              * tmask(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 16) = wsbio3  (:,:,:)              * tmask(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 17) = wsbio4  (:,:,:)              * tmask(:,:,:)
         ENDIF
      ENDIF

      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, znum3d )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!
      !! ** Purpose :   Initialization of sinking parameters
      !!                Kriest parameterization only
      !!
      !! ** Method  :   Read the nampiskrs namelist and check the parameters
      !!      called at the first timestep 
      !!
      !! ** input   :   Namelist nampiskrs
      !!----------------------------------------------------------------------
      INTEGER  ::   jk, jn, kiter
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      REAL(wp) ::   znum, zdiv
      REAL(wp) ::   zws, zwr, zwl,wmax, znummax
      REAL(wp) ::   zmin, zmax, zl, zr, xacc
      !
      NAMELIST/nampiskrs/ xkr_sfact, xkr_stick ,  &
         &                xkr_nnano, xkr_ndiat, xkr_nmicro, xkr_nmeso, xkr_naggr
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink_init')
      !

      REWIND( numnatp_ref )              ! Namelist nampiskrs in reference namelist : Pisces sinking Kriest
      READ  ( numnatp_ref, nampiskrs, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiskrs in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampiskrs in configuration namelist : Pisces sinking Kriest
      READ  ( numnatp_cfg, nampiskrs, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiskrs in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiskrs )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampiskrs'
         WRITE(numout,*) '    Sinking factor                           xkr_sfact    = ', xkr_sfact
         WRITE(numout,*) '    Stickiness                               xkr_stick    = ', xkr_stick
         WRITE(numout,*) '    Nbr of cell in nano size class           xkr_nnano    = ', xkr_nnano
         WRITE(numout,*) '    Nbr of cell in diatoms size class        xkr_ndiat    = ', xkr_ndiat
         WRITE(numout,*) '    Nbr of cell in microzoo size class       xkr_nmicro   = ', xkr_nmicro
         WRITE(numout,*) '    Nbr of cell in mesozoo size class        xkr_nmeso    = ', xkr_nmeso
         WRITE(numout,*) '    Nbr of cell in aggregates size class     xkr_naggr    = ', xkr_naggr
      ENDIF


      ! max and min vertical particle speed
      xkr_wsbio_min = xkr_sfact * xkr_mass_min**xkr_eta
      xkr_wsbio_max = xkr_sfact * xkr_mass_max**xkr_eta
      IF (lwp) WRITE(numout,*) ' max and min vertical particle speed ', xkr_wsbio_min, xkr_wsbio_max

      !
      !    effect of the sizes of the different living pools on particle numbers
      !    nano = 2um-20um -> mean size=6.32 um -> ws=2.596 -> xnum=xnnano=2.337
      !    diat and microzoo = 10um-200um -> 44.7 -> 8.732 -> xnum=xndiat=3.718
      !    mesozoo = 200um-2mm -> 632.45 -> 45.14 -> xnum=xnmeso=7.147
      !    aggregates = 200um-10mm -> 1414 -> 74.34 -> xnum=xnaggr=9.877
      !    doc aggregates = 1um
      ! ----------------------------------------------------------

      xkr_dnano = 1. / ( xkr_massp * xkr_nnano )
      xkr_ddiat = 1. / ( xkr_massp * xkr_ndiat )
      xkr_dmicro = 1. / ( xkr_massp * xkr_nmicro )
      xkr_dmeso = 1. / ( xkr_massp * xkr_nmeso )
      xkr_daggr = 1. / ( xkr_massp * xkr_naggr )

      !!---------------------------------------------------------------------
      !!    'key_kriest'                                                  ???
      !!---------------------------------------------------------------------
      !  COMPUTATION OF THE VERTICAL PROFILE OF MAXIMUM SINKING SPEED
      !  Search of the maximum number of particles in aggregates for each k-level.
      !  Bissection Method
      !--------------------------------------------------------------------
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*)'    kriest : Compute maximum number of particles in aggregates'
      ENDIF

      xacc     =  0.001_wp
      kiter    = 50
      zmin     =  1.10_wp
      zmax     = xkr_mass_max / xkr_mass_min
      xkr_frac = zmax

      DO jk = 1,jpk
         zl = zmin
         zr = zmax
         wmax = 0.5 * fse3t(1,1,jk) * rday * float(niter1max) / rfact2
         zdiv = xkr_zeta + xkr_eta - xkr_eta * zl
         znum = zl - 1.
         zwl =  xkr_wsbio_min * xkr_zeta / zdiv &
            & - ( xkr_wsbio_max * xkr_eta * znum * &
            &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
            & - wmax

         zdiv = xkr_zeta + xkr_eta - xkr_eta * zr
         znum = zr - 1.
         zwr =  xkr_wsbio_min * xkr_zeta / zdiv &
            & - ( xkr_wsbio_max * xkr_eta * znum * &
            &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
            & - wmax
iflag:   DO jn = 1, kiter
            IF    ( zwl == 0._wp ) THEN   ;   znummax = zl
            ELSEIF( zwr == 0._wp ) THEN   ;   znummax = zr
            ELSE
               znummax = ( zr + zl ) / 2.
               zdiv = xkr_zeta + xkr_eta - xkr_eta * znummax
               znum = znummax - 1.
               zws =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax
               IF( zws * zwl < 0. ) THEN   ;   zr = znummax
               ELSE                        ;   zl = znummax
               ENDIF
               zdiv = xkr_zeta + xkr_eta - xkr_eta * zl
               znum = zl - 1.
               zwl =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax

               zdiv = xkr_zeta + xkr_eta - xkr_eta * zr
               znum = zr - 1.
               zwr =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax
               !
               IF ( ABS ( zws )  <= xacc ) EXIT iflag
               !
            ENDIF
            !
         END DO iflag

         xnumm(jk) = znummax
         IF (lwp) WRITE(numout,*) '       jk = ', jk, ' wmax = ', wmax,' xnum max = ', xnumm(jk)
         !
      END DO
      !
      ik100 = 10        !  last level where depth less than 100 m
      DO jk = jpkm1, 1, -1
         IF( gdept_1d(jk) > 100. )  iksed = jk - 1
      END DO
      IF (lwp) WRITE(numout,*)
      IF (lwp) WRITE(numout,*) ' Level corresponding to 100m depth ',  ik100 + 1
      IF (lwp) WRITE(numout,*)
      !
      t_oce_co2_exp = 0._wp
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink_init')
      !
  END SUBROUTINE p4z_sink_init

#endif

   SUBROUTINE p4z_sink2( pwsink, psinkflx, jp_tra, kiter )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index      
      INTEGER , INTENT(in   )                         ::   kiter     ! number of iterations for time-splitting 
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztraz, zakz, zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink2')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )

      zstep = rfact2 / FLOAT( kiter ) / 2.

      ztraz(:,:,:) = 0.e0
      zakz (:,:,:) = 0.e0
      ztrb (:,:,:) = trb(:,:,:,jp_tra)

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0
      IF( lk_degrad ) THEN
         zwsink2(:,:,:) = zwsink2(:,:,:) * facvol(:,:,:)
      ENDIF


      ! Vertical advective flux
      DO jn = 1, 2
         !  first guess of the slopes interior values
         DO jk = 2, jpkm1
            ztraz(:,:,jk) = ( trb(:,:,jk-1,jp_tra) - trb(:,:,jk,jp_tra) ) * tmask(:,:,jk)
         END DO
         ztraz(:,:,1  ) = 0.0
         ztraz(:,:,jpk) = 0.0

         ! slopes
         DO jk = 2, jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
                  zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
               END DO
            END DO
         END DO
         
         ! Slopes limitation
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
                     &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
               END DO
            END DO
         END DO
         
         ! vertical advective flux
         DO jk = 1, jpkm1
            DO jj = 1, jpj      
               DO ji = 1, jpi    
                  zigma = zwsink2(ji,jj,jk+1) * zstep / fse3w(ji,jj,jk+1)
                  zew   = zwsink2(ji,jj,jk+1)
                  psinkflx(ji,jj,jk+1) = -zew * ( trb(ji,jj,jk,jp_tra) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
               END DO
            END DO
         END DO
         !
         ! Boundary conditions
         psinkflx(:,:,1  ) = 0.e0
         psinkflx(:,:,jpk) = 0.e0
         
         DO jk=1,jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
                  trb(ji,jj,jk,jp_tra) = trb(ji,jj,jk,jp_tra) + zflx
               END DO
            END DO
         END DO

      ENDDO

      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + 2. * zflx
            END DO
         END DO
      END DO

      trb(:,:,:,jp_tra) = ztrb(:,:,:)
      psinkflx(:,:,:)   = 2. * psinkflx(:,:,:)
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink2')
      !
   END SUBROUTINE p4z_sink2


   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wsbio3 (jpi,jpj,jpk) , wsbio4  (jpi,jpj,jpk) , wscal(jpi,jpj,jpk) ,     &
         &      sinking(jpi,jpj,jpk) , sinking2(jpi,jpj,jpk)                      ,     &                
         &      sinkcal(jpi,jpj,jpk) , sinksil (jpi,jpj,jpk)                      ,     &                
#if defined key_kriest
         &      xnumm(jpk)                                                        ,     &                
#else
         &      sinkfer2(jpi,jpj,jpk)                                             ,     &                
#endif
         &      sinkfer(jpi,jpj,jpk)                                              , STAT=p4z_sink_alloc )                
         !
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_sink_alloc
   
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif 

   !!======================================================================
END MODULE p4zsink
