MODULE trcsbc
   !!==============================================================================
   !!                       ***  MODULE  trcsbc  ***
   !! Ocean passive tracers:  surface boundary condition
   !!======================================================================
   !! History :  8.2  !  1998-10  (G. Madec, G. Roullet, M. Imbard)  Original code
   !!            8.2  !  2001-02  (D. Ludicone)  sea ice and free surface
   !!            8.5  !  2002-06  (G. Madec)  F90: Free form and module
   !!            9.0  !  2004-03  (C. Ethe)  adapted for passive tracers
   !!                 !  2006-08  (C. Deltel) Diagnose ML trends for passive tracers
   !!==============================================================================
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_sbc      : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean  passive tracers variables
   USE prtctl_trc      ! Print control for debbuging
   USE iom
   USE trd_oce
   USE trdtra

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sbc   ! routine called by step.F90

   REAL(wp) ::   r2dt  !  time-step at surface

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sbc ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sbc  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (concentration/dilution effect) and add it to the general 
      !!       trend of tracer equations.
      !!
      !! ** Method :
      !!      * concentration/dilution effect:
      !!            The surface freshwater flux modify the ocean volume
      !!         and thus the concentration of a tracer as :
      !!            tra = tra + emp * trn / e3t   for k=1
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation ) given in kg/m2/s is divided
      !!         by 1035 kg/m3 (density of ocean water) to obtain m/s.
      !!
      !! ** Action  : - Update the 1st level of tra with the trend associated
      !!                with the tracer surface boundary condition 
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jn                                     ! dummy loop indices
      REAL(wp) ::   zse3t, zrtrn, zratio, zfact                    ! temporary scalars
      REAL(wp) ::   zswitch, zftra, zcd, zdtra, ztfx, ztra         ! temporary scalars
      CHARACTER (len=22) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zsfx
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrtrd

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sbc')
      !
      ! Allocate temporary workspace
                      CALL wrk_alloc( jpi, jpj,      zsfx   )
      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrtrd )
      !
      zrtrn = 1.e-15_wp

      SELECT CASE( nn_ice_embd )         ! levitating or embedded sea-ice option
         CASE( 0    )   ;   zswitch = 1  ! (0) standard levitating sea-ice : salt exchange only
         CASE( 1, 2 )   ;   zswitch = 0  ! (1) levitating sea-ice: salt and volume exchange but no pressure effect                                
                                         ! (2) embedded sea-ice : salt and volume fluxes and pressure
      END SELECT

      IF( ln_top_euler) THEN
         r2dt =  rdttrc(1)              ! = rdttrc (use Euler time stepping)
      ELSE
         IF( neuler == 0 .AND. kt == nittrc000 ) THEN     ! at nittrc000
            r2dt = rdttrc(1)           ! = rdttrc (restarting with Euler time stepping)
         ELSEIF( kt <= nittrc000 + nn_dttrc ) THEN          ! at nittrc000 or nittrc000+1
            r2dt = 2. * rdttrc(1)       ! = 2 rdttrc (leapfrog)
         ENDIF
      ENDIF


      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_sbc : Passive tracers surface boundary condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '

         !! JPALM -- 12-01-2016 -- problem after restart, maybe because of this...
         !!                     -- set sbc_trc_b to 0 after restart, first, to check.
         !!------------------------------------------------------------------------------
        ! IF( ln_rsttr .AND. .NOT.ln_top_euler .AND.   &                     ! Restart: read in restart  file
        !    iom_varid( numrtr, 'sbc_'//TRIM(ctrcnm(1))//'_b', ldstop = .FALSE. ) > 0 ) THEN
        !    IF(lwp) WRITE(numout,*) '          nittrc000-nn_dttrc surface tracer content forcing fields red in the restart file'
        !    zfact = 0.5_wp
        !    DO jn = 1, jptra
        !       CALL iom_get( numrtr, jpdom_autoglo, 'sbc_'//TRIM(ctrcnm(jn))//'_b', sbc_trc_b(:,:,jn) )   ! before tracer content sbc
        !    END DO
        ! ELSE                                         ! No restart or restart not found: Euler forward time stepping
           zfact = 1._wp
           sbc_trc_b(:,:,:) = 0._wp
        ! ENDIF
      ELSE                                         ! Swap of forcing fields
         IF( ln_top_euler ) THEN
            zfact = 1._wp
            sbc_trc_b(:,:,:) = 0._wp
         ELSE
            zfact = 0.5_wp
            sbc_trc_b(:,:,:) = sbc_trc(:,:,:)
         ENDIF
         !
      ENDIF

      ! Coupling online : river runoff is added to the horizontal divergence (hdivn) in the subroutine sbc_rnf_div 
      ! one only consider the concentration/dilution effect due to evaporation minus precipitation + freezing/melting of sea-ice
      ! Coupling offline : runoff are in emp which contains E-P-R
      !
      IF( .NOT. lk_offline .AND. lk_vvl ) THEN  ! online coupling with vvl
         zsfx(:,:) = 0._wp
      ELSE                                      ! online coupling free surface or offline with free surface
         zsfx(:,:) = emp(:,:)
      ENDIF

      ! 0. initialization
      DO jn = 1, jptra
         !
         IF( l_trdtrc ) THEN
            ztrtrd(:,:,:) = 0.0
            ztrtrd(:,:,1) = tra(:,:,1,jn)  ! save surface trends
         !                                 ! add the trend to the general tracer trend
         ENDIF

         IF ( nn_ice_tr == -1 ) THEN  ! No tracers in sea ice (null concentration in sea ice)

            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  sbc_trc(ji,jj,jn) = zsfx(ji,jj) * r1_rau0 * trn(ji,jj,1,jn)
               END DO
            END DO

         ELSE

            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zse3t = 1. / fse3t(ji,jj,1)
                  ! tracer flux at the ice/ocean interface (tracer/m2/s)
                  zftra = - trc_i(ji,jj,jn) * fmmflx(ji,jj) ! uptake of tracer in the sea ice
                  zcd   =   trc_o(ji,jj,jn) * fmmflx(ji,jj) ! concentration dilution due to freezing-melting,
                                                               ! only used in the levitating sea ice case
                  ! tracer flux only       : add concentration dilution term in net tracer flux, no F-M in volume flux
                  ! tracer and mass fluxes : no concentration dilution term in net tracer flux, F-M term in volume flux
                  ztfx  = zftra + zswitch * zcd                ! net tracer flux (+C/D if no ice/ocean mass exchange)
   
                  zdtra = r1_rau0 * ( ztfx + zsfx(ji,jj) * trn(ji,jj,1,jn) ) 
                  IF ( zdtra < 0. ) THEN
                     zratio = -zdtra * zse3t * r2dt / ( trn(ji,jj,1,jn) + zrtrn )
                     zdtra = MIN(1.0, zratio) * zdtra ! avoid negative concentrations to arise
                  ENDIF
                  sbc_trc(ji,jj,jn) =  zdtra 
               END DO
            END DO
         ENDIF
         !
         CALL lbc_lnk( sbc_trc(:,:,jn), 'T', 1. )
         !                                       Concentration dilution effect on tracers due to evaporation & precipitation 
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zse3t = zfact / fse3t(ji,jj,1)
               tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ( sbc_trc_b(ji,jj,jn) + sbc_trc(ji,jj,jn) ) * zse3t
            END DO
         END DO
         !
         IF( l_trdtrc ) THEN
            ztrtrd(:,:,1) = tra(:,:,1,jn) - ztrtrd(:,:,1)
            CALL trd_tra( kt, 'TRC', jn, jptra_nsr, ztrtrd )
         END IF
         !                                                       ! ===========
      END DO                                                     ! tracer loop
      !                                                          ! ===========

      !                                           Write in the tracer restar  file
      !                                          *******************************
      IF( lrst_trc .AND. .NOT.ln_top_euler ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : ocean surface tracer content forcing fields written in tracer restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         DO jn = 1, jptra
            CALL iom_rstput( kt, nitrst, numrtw, 'sbc_'//TRIM(ctrcnm(jn))//'_b', sbc_trc(:,:,jn) )
         END DO
      ENDIF
      !
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('sbc ')") ;  CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      ENDIF
                      CALL wrk_dealloc( jpi, jpj,      zsfx   )
      IF( l_trdtrc )  CALL wrk_dealloc( jpi, jpj, jpk, ztrtrd )
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sbc')
      !
   END SUBROUTINE trc_sbc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sbc (kt)              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_sbc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sbc
#endif
   
   !!======================================================================
END MODULE trcsbc
