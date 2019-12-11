MODULE trdken
   !!======================================================================
   !!                       ***  MODULE  trdken  ***
   !! Ocean diagnostics:  compute and output 3D kinetic energy trends
   !!=====================================================================
   !! History :  3.5  !  2012-02  (G. Madec) original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_ken       : compute and output 3D Kinetic energy trends using IOM
   !!   trd_ken_init  : initialisation
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics variables
   USE trd_oce        ! trends: ocean variables
!!gm   USE dynhpg          ! hydrostatic pressure gradient   
   USE zdfbfr         ! bottom friction
   USE ldftra_oce     ! ocean active tracers lateral physics
   USE sbc_oce        ! surface boundary condition: ocean
   USE phycst         ! physical constants
   USE trdvor         ! ocean vorticity trends 
   USE trdglo         ! trends:global domain averaged
   USE trdmxl         ! ocean active mixed layer tracers trends 
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory allocation
   USE ldfslp         ! Isopycnal slopes

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_ken       ! called by trddyn module
   PUBLIC   trd_ken_init  ! called by trdini module

   INTEGER ::   nkstp       ! current time step 

   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   bu, bv   ! volume of u- and v-boxes
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   r1_bt    ! inverse of t-box volume

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
#  include "ldfeiv_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_ken_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION trd_ken_alloc  ***
      !!---------------------------------------------------------------------
      ALLOCATE( bu(jpi,jpj,jpk) , bv(jpi,jpj,jpk) , r1_bt(jpi,jpj,jpk) , STAT= trd_ken_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( trd_ken_alloc )
      IF( trd_ken_alloc /= 0 )   CALL ctl_warn('trd_ken_alloc: failed to allocate arrays')
   END FUNCTION trd_ken_alloc


   SUBROUTINE trd_ken( putrd, pvtrd, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_ken  ***
      !! 
      !! ** Purpose :   output 3D Kinetic Energy trends using IOM
      !!
      !! ** Method  : - apply lbc to the input masked velocity trends 
      !!              - compute the associated KE trend:
      !!          zke = 0.5 * (  mi-1[ un * putrd * bu ] + mj-1[ vn * pvtrd * bv]  ) / bt
      !!      where bu, bv, bt are the volume of u-, v- and t-boxes. 
      !!              - vertical diffusion case (jpdyn_zdf): 
      !!          diagnose separately the KE trend associated with wind stress
      !!              - bottom friction case (jpdyn_bfr):
      !!          explicit case (ln_bfrimp=F): bottom trend put in the 1st level 
      !!                                       of putrd, pvtrd
      !
      !
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   putrd, pvtrd   ! U and V masked trends
      INTEGER                   , INTENT(in   ) ::   ktrd           ! trend index
      INTEGER                   , INTENT(in   ) ::   kt             ! time step
      !
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   ikbu  , ikbv     ! local integers
      INTEGER ::   ikbum1, ikbvm1   !   -       -
      REAL(wp), POINTER, DIMENSION(:,:)   ::   z2dx, z2dy, zke2d   ! 2D workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zke                 ! 3D workspace 
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi, jpj, jpk, zke )
      !
      CALL lbc_lnk( putrd, 'U', -1. )   ;   CALL lbc_lnk( pvtrd, 'V', -1. )      ! lateral boundary conditions
      !
      IF ( lk_vvl .AND. kt /= nkstp ) THEN   ! Variable volume: set box volume at the 1st call of kt time step
         nkstp = kt
         DO jk = 1, jpkm1
            bu   (:,:,jk) =  e1u(:,:) * e2u(:,:) * fse3u_n(:,:,jk)
            bv   (:,:,jk) =  e1v(:,:) * e2v(:,:) * fse3v_n(:,:,jk)
            r1_bt(:,:,jk) = 1._wp / ( e1e2t(:,:) * fse3t_n(:,:,jk) ) * tmask(:,:,jk)
         END DO
      ENDIF
      !
      zke(:,:,jpk) = 0._wp
      zke(1,:, : ) = 0._wp
      zke(:,1, : ) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpj
            DO ji = 2, jpi
               zke(ji,jj,jk) = 0.5_wp * rau0 *( un(ji  ,jj,jk) * putrd(ji  ,jj,jk) * bu(ji  ,jj,jk)  &
                  &                           + un(ji-1,jj,jk) * putrd(ji-1,jj,jk) * bu(ji-1,jj,jk)  &
                  &                           + vn(ji,jj  ,jk) * pvtrd(ji,jj  ,jk) * bv(ji,jj  ,jk)  &
                  &                           + vn(ji,jj-1,jk) * pvtrd(ji,jj-1,jk) * bv(ji,jj-1,jk)  ) * r1_bt(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      SELECT CASE( ktrd )
        CASE( jpdyn_hpg )   ;   CALL iom_put( "ketrd_hpg", zke )    ! hydrostatic pressure gradient
        CASE( jpdyn_spg )   ;   CALL iom_put( "ketrd_spg", zke )    ! surface pressure gradient
        CASE( jpdyn_spgexp );   CALL iom_put( "ketrd_spgexp", zke ) ! surface pressure gradient (explicit)
        CASE( jpdyn_spgflt );   CALL iom_put( "ketrd_spgflt", zke ) ! surface pressure gradient (filter)
        CASE( jpdyn_pvo )   ;   CALL iom_put( "ketrd_pvo", zke )    ! planetary vorticity
        CASE( jpdyn_rvo )   ;   CALL iom_put( "ketrd_rvo", zke )    ! relative  vorticity     (or metric term)
        CASE( jpdyn_keg )   ;   CALL iom_put( "ketrd_keg", zke )    ! Kinetic Energy gradient (or had)
        CASE( jpdyn_zad )   ;   CALL iom_put( "ketrd_zad", zke )    ! vertical   advection
        CASE( jpdyn_ldf )   ;   CALL iom_put( "ketrd_ldf", zke )    ! lateral diffusion
        CASE( jpdyn_zdf )   ;   CALL iom_put( "ketrd_zdf", zke )    ! vertical diffusion 
                                 !                                   ! wind stress trends
                                CALL wrk_alloc( jpi, jpj, z2dx, z2dy, zke2d )
                     z2dx(:,:) = un(:,:,1) * ( utau_b(:,:) + utau(:,:) ) * e1u(:,:) * e2u(:,:) * umask(:,:,1)
                     z2dy(:,:) = vn(:,:,1) * ( vtau_b(:,:) + vtau(:,:) ) * e1v(:,:) * e2v(:,:) * vmask(:,:,1)
                     zke2d(1,:) = 0._wp   ;   zke2d(:,1) = 0._wp
                     DO jj = 2, jpj
                         DO ji = 2, jpi
                           zke2d(ji,jj) = 0.5_wp * (   z2dx(ji,jj) + z2dx(ji-1,jj)   &
                            &                         + z2dy(ji,jj) + z2dy(ji,jj-1)   ) * r1_bt(ji,jj,1)
                         END DO
                     END DO
                                CALL iom_put( "ketrd_tau", zke2d )
                                CALL wrk_dealloc( jpi, jpj     , z2dx, z2dy, zke2d )
        CASE( jpdyn_bfr )   ;   CALL iom_put( "ketrd_bfr", zke )    ! bottom friction (explicit case) 
!!gm TO BE DONE properly
!!gm only valid if ln_bfrimp=F otherwise the bottom stress as to be recomputed at the end of the computation....
!         IF(.NOT. ln_bfrimp) THEN
!            DO jj = 1, jpj    !   
!               DO ji = 1, jpi
!                  ikbu = mbku(ji,jj)         ! deepest ocean u- & v-levels
!                  ikbv = mbkv(ji,jj)   
!                  z2dx(ji,jj) = un(ji,jj,ikbu) * bfrua(ji,jj) * un(ji,jj,ikbu)
!                  z2dy(ji,jj) = vn(ji,jj,ikbu) * bfrva(ji,jj) * vn(ji,jj,ikbv)
!               END DO
!            END DO
!            zke2d(1,:) = 0._wp   ;   zke2d(:,1) = 0._wp
!            DO jj = 2, jpj
!               DO ji = 2, jpi
!                  zke2d(ji,jj) = 0.5_wp * (   z2dx(ji,jj) + z2dx(ji-1,jj)   &
!                     &                      + z2dy(ji,jj) + z2dy(ji,jj-1)   ) * r1_bt(ji,jj,  BEURK!!!
!               END DO
!            END DO
!                              CALL iom_put( "ketrd_bfr", zke2d )    ! bottom friction (explicit case)
!         ENDIF
!!gm end
        CASE( jpdyn_atf )   ;   CALL iom_put( "ketrd_atf", zke )    ! asselin filter trends 
!! a faire !!!!  idee changer dynnxt pour avoir un appel a jpdyn_bfr avant le swap !!!
!! reflechir a une possible sauvegarde du "vrai" un,vn pour le calcul de atf....
!
!         IF( ln_bfrimp ) THEN                                          ! bottom friction (implicit case)
!            DO jj = 1, jpj                                                  ! after velocity known (now filed at this stage)
!               DO ji = 1, jpi
!                  ikbu = mbku(ji,jj)          ! deepest ocean u- & v-levels
!                  ikbv = mbkv(ji,jj)
!                  z2dx(ji,jj) = un(ji,jj,ikbu) * bfrua(ji,jj) * un(ji,jj,ikbu) / fse3u(ji,jj,ikbu)
!                  z2dy(ji,jj) = un(ji,jj,ikbu) * bfrva(ji,jj) * vn(ji,jj,ikbv) / fse3v(ji,jj,ikbv)
!               END DO
!            END DO
!            zke2d(1,:) = 0._wp   ;   zke2d(:,1) = 0._wp
!            DO jj = 2, jpj
!               DO ji = 2, jpi
!                  zke2d(ji,jj) = 0.5_wp * (   z2dx(ji,jj) + z2dx(ji-1,jj)   &
!                     &                      + z2dy(ji,jj) + z2dy(ji,jj-1)   )
!               END DO
!            END DO
!                              CALL iom_put( "ketrd_bfri", zke2d )
!         ENDIF
        CASE( jpdyn_ken )   ;   ! kinetic energy
                    ! called in dynnxt.F90 before asselin time filter
                    ! with putrd=ua and pvtrd=va
                    zke(:,:,:) = 0.5_wp * zke(:,:,:)
                    CALL iom_put( "KE", zke )
                    !
                    CALL ken_p2k( kt , zke )
                      CALL iom_put( "ketrd_convP2K", zke )     ! conversion -rau*g*w
        CASE( jpdyn_eivke )
            ! CMIP6 diagnostic tknebto = tendency of KE from
            ! parameterized mesoscale eddy advection
            ! = vertical_integral( k (N S)^2 ) rho dz
            ! rho = reference density
            ! S = isoneutral slope.
            ! Most terms are on W grid so work on this grid
            CALL wrk_alloc( jpi, jpj, zke2d )
            zke2d(:,:) = 0._wp
            DO jk = 1,jpk
               DO ji = 1,jpi
                  DO jj = 1,jpj
                     zke2d(ji,jj) = zke2d(ji,jj) +  rau0 * fsaeiw(ji, jj, jk)               &
                          &                      * ( wslpi(ji, jj, jk) * wslpi(ji,jj,jk)    &
                          &                      +   wslpj(ji, jj, jk) * wslpj(ji,jj,jk) )  &
                          &                      *   rn2(ji,jj,jk) * fse3w(ji, jj, jk)
                  ENDDO
               ENDDO
            ENDDO
            CALL iom_put("ketrd_eiv", zke2d)
            CALL wrk_dealloc( jpi, jpj, zke2d )
         !
      END SELECT
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zke )
      !
   END SUBROUTINE trd_ken


   SUBROUTINE ken_p2k( kt , pconv )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE ken_p2k  ***
      !!                    
      !! ** Purpose :   compute rate of conversion from potential to kinetic energy
      !!
      !! ** Method  : - compute conv defined as -rau*g*w on T-grid points
      !! 
      !! ** Work only for full steps and partial steps (ln_hpg_zco or ln_hpg_zps)
      !!---------------------------------------------------------------------- 
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(  out) ::   pconv
      !
      INTEGER  ::   ji, jj, jk                       ! dummy loop indices
      INTEGER  ::   iku, ikv                         ! temporary integers
      REAL(wp) ::   zcoef                            ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zconv  ! temporary conv on W-grid
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj,jpk, zconv )
      !
      ! Local constant initialization 
      zcoef = - rau0 * grav * 0.5_wp      
      
      !  Surface value (also valid in partial step case)
      zconv(:,:,1) = zcoef * ( 2._wp * rhd(:,:,1) ) * wn(:,:,1) * fse3w(:,:,1)

      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpk
         zconv(:,:,jk) = zcoef * ( rhd(:,:,jk) + rhd(:,:,jk-1) ) * wn(:,:,jk) * fse3w(:,:,jk)
      END DO

      ! conv value on T-point
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcoef = 0.5_wp / fse3t(ji,jj,jk)
               pconv(ji,jj,jk) = zcoef * ( zconv(ji,jj,jk) + zconv(ji,jj,jk+1) ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zconv )      
      !
   END SUBROUTINE ken_p2k


   SUBROUTINE trd_ken_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_ken_init  ***
      !! 
      !! ** Purpose :   initialisation of 3D Kinetic Energy trend diagnostic
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trd_ken_init : 3D Kinetic Energy trends'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !                           ! allocate box volume arrays
      IF ( trd_ken_alloc() /= 0 )   CALL ctl_stop('trd_ken_alloc: failed to allocate arrays')
      !
!!gm      IF( .NOT. (ln_hpg_zco.OR.ln_hpg_zps) )   &
!!gm         &   CALL ctl_stop('trd_ken_init : only full and partial cells are coded for conversion rate')
      !
      IF ( .NOT.lk_vvl ) THEN     ! constant volume: bu, bv, 1/bt computed one for all
         DO jk = 1, jpkm1
            bu   (:,:,jk) =  e1u(:,:) * e2u(:,:) * fse3u_n(:,:,jk)
            bv   (:,:,jk) =  e1v(:,:) * e2v(:,:) * fse3v_n(:,:,jk)
            r1_bt(:,:,jk) = 1._wp / ( e1e2t(:,:) * fse3t_n(:,:,jk) )
         END DO
      ENDIF
      !
   END SUBROUTINE trd_ken_init

   !!======================================================================
END MODULE trdken
