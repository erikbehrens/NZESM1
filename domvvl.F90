MODULE domvvl
   !!======================================================================
   !!                       ***  MODULE domvvl   ***
   !! Ocean : 
   !!======================================================================
   !! History :  2.0  !  2006-06  (B. Levier, L. Marie)  original code
   !!            3.1  !  2009-02  (G. Madec, M. Leclair, R. Benshila)  pure z* coordinate
   !!            3.3  !  2011-10  (M. Leclair) totally rewrote domvvl:
   !!                                          vvl option includes z_star and z_tilde coordinates
   !!            3.6  !  2014-11  (P. Mathiot) add ice shelf capability
   !!----------------------------------------------------------------------
   !!   'key_vvl'                              variable volume
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   dom_vvl_init     : define initial vertical scale factors, depths and column thickness
   !!   dom_vvl_sf_nxt   : Compute next vertical scale factors
   !!   dom_vvl_sf_swp   : Swap vertical scale factors and update the vertical grid
   !!   dom_vvl_interpol : Interpolate vertical scale factors from one grid point to another
   !!   dom_vvl_rst      : read/write restart file
   !!   dom_vvl_ctl      : Check the vvl options
   !!   dom_vvl_orca_fix : Recompute some area-weighted interpolations of vertical scale factors 
   !!                    : to account for manual changes to e[1,2][u,v] in some Straits 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! ocean surface boundary condition
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O manager library
   USE restart         ! ocean restart
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC  dom_vvl_init       ! called by domain.F90
   PUBLIC  dom_vvl_sf_nxt     ! called by step.F90
   PUBLIC  dom_vvl_sf_swp     ! called by step.F90
   PUBLIC  dom_vvl_interpol   ! called by dynnxt.F90
   PRIVATE dom_vvl_orca_fix   ! called by dom_vvl_interpol

   !!* Namelist nam_vvl
   LOGICAL , PUBLIC                                      :: ln_vvl_zstar = .FALSE.              ! zstar  vertical coordinate
   LOGICAL , PUBLIC                                      :: ln_vvl_ztilde = .FALSE.             ! ztilde vertical coordinate
   LOGICAL , PUBLIC                                      :: ln_vvl_layer = .FALSE.              ! level  vertical coordinate
   LOGICAL , PUBLIC                                      :: ln_vvl_ztilde_as_zstar = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC                                      :: ln_vvl_zstar_at_eqtor = .FALSE.     ! ztilde vertical coordinate
   LOGICAL , PUBLIC                                      :: ln_vvl_kepe = .FALSE.               ! kinetic/potential energy transfer
   !                                                                                            ! conservation: not used yet
   REAL(wp)                                              :: rn_ahe3                   ! thickness diffusion coefficient
   REAL(wp)                                              :: rn_rst_e3t                ! ztilde to zstar restoration timescale [days]
   REAL(wp)                                              :: rn_lf_cutoff              ! cutoff frequency for low-pass filter  [days]
   REAL(wp)                                              :: rn_zdef_max               ! maximum fractional e3t deformation
   LOGICAL , PUBLIC                                      :: ln_vvl_dbg = .FALSE.      ! debug control prints

   !! * Module variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: un_td, vn_td                       ! thickness diffusion transport
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: hdiv_lf                            ! low frequency part of hz divergence
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_b, tilde_e3t_n           ! baroclinic scale factors
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_a, dtilde_e3t_a          ! baroclinic scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_e3t                        ! retoring period for scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_hdv                        ! retoring period for low freq. divergence

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO-Consortium (2010) 
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   INTEGER FUNCTION dom_vvl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION dom_vvl_alloc  ***
      !!----------------------------------------------------------------------
      IF( ln_vvl_zstar ) dom_vvl_alloc = 0
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         ALLOCATE( tilde_e3t_b(jpi,jpj,jpk)  , tilde_e3t_n(jpi,jpj,jpk) , tilde_e3t_a(jpi,jpj,jpk) ,   &
            &      dtilde_e3t_a(jpi,jpj,jpk) , un_td  (jpi,jpj,jpk)     , vn_td  (jpi,jpj,jpk)     ,   &
            &      STAT = dom_vvl_alloc        )
         IF( lk_mpp             )   CALL mpp_sum ( dom_vvl_alloc )
         IF( dom_vvl_alloc /= 0 )   CALL ctl_warn('dom_vvl_alloc: failed to allocate arrays')
         un_td = 0.0_wp
         vn_td = 0.0_wp
      ENDIF
      IF( ln_vvl_ztilde ) THEN
         ALLOCATE( frq_rst_e3t(jpi,jpj) , frq_rst_hdv(jpi,jpj) , hdiv_lf(jpi,jpj,jpk) , STAT= dom_vvl_alloc )
         IF( lk_mpp             )   CALL mpp_sum ( dom_vvl_alloc )
         IF( dom_vvl_alloc /= 0 )   CALL ctl_warn('dom_vvl_alloc: failed to allocate arrays')
      ENDIF

   END FUNCTION dom_vvl_alloc


   SUBROUTINE dom_vvl_init
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_init  ***
      !!                   
      !! ** Purpose :  Initialization of all scale factors, depths
      !!               and water column heights
      !!
      !! ** Method  :  - use restart file and/or initialize
      !!               - interpolate scale factors
      !!
      !! ** Action  : - fse3t_(n/b) and tilde_e3t_(n/b)
      !!              - Regrid: fse3(u/v)_n
      !!                        fse3(u/v)_b       
      !!                        fse3w_n           
      !!                        fse3(u/v)w_b      
      !!                        fse3(u/v)w_n      
      !!                        fsdept_n, fsdepw_n and fsde3w_n
      !!              - h(t/u/v)_0
      !!              - frq_rst_e3t and frq_rst_hdv
      !!
      !! Reference  : Leclair, M., and G. Madec, 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      USE phycst,  ONLY : rpi, rsmall, rad
      !! * Local declarations
      INTEGER ::   ji,jj,jk
      INTEGER ::   ii0, ii1, ij0, ij1
      REAL(wp)::   zcoef
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('dom_vvl_init')

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dom_vvl_init : Variable volume activated'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

      ! choose vertical coordinate (z_star, z_tilde or layer)
      ! ==========================
      CALL dom_vvl_ctl

      ! Allocate module arrays
      ! ======================
      IF( dom_vvl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dom_vvl_init : unable to allocate arrays' )

      ! Read or initialize fse3t_(b/n), tilde_e3t_(b/n) and hdiv_lf (and e3t_a(jpk))
      ! ============================================================================
      CALL dom_vvl_rst( nit000, 'READ' )
      fse3t_a(:,:,jpk) = e3t_0(:,:,jpk)

      ! Reconstruction of all vertical scale factors at now and before time steps
      ! =============================================================================
      ! Horizontal scale factor interpolations
      ! --------------------------------------
      CALL dom_vvl_interpol( fse3t_b(:,:,:), fse3u_b(:,:,:), 'U' )
      CALL dom_vvl_interpol( fse3t_b(:,:,:), fse3v_b(:,:,:), 'V' )
      CALL dom_vvl_interpol( fse3t_n(:,:,:), fse3u_n(:,:,:), 'U' )
      CALL dom_vvl_interpol( fse3t_n(:,:,:), fse3v_n(:,:,:), 'V' )
      CALL dom_vvl_interpol( fse3u_n(:,:,:), fse3f_n(:,:,:), 'F' )
      ! Vertical scale factor interpolations
      ! ------------------------------------
      CALL dom_vvl_interpol( fse3t_n(:,:,:), fse3w_n (:,:,:), 'W'  )
      CALL dom_vvl_interpol( fse3u_n(:,:,:), fse3uw_n(:,:,:), 'UW' )
      CALL dom_vvl_interpol( fse3v_n(:,:,:), fse3vw_n(:,:,:), 'VW' )
      CALL dom_vvl_interpol( fse3t_b(:,:,:), fse3w_b (:,:,:), 'W'  )
      CALL dom_vvl_interpol( fse3u_b(:,:,:), fse3uw_b(:,:,:), 'UW' )
      CALL dom_vvl_interpol( fse3v_b(:,:,:), fse3vw_b(:,:,:), 'VW' )
      ! t- and w- points depth
      ! ----------------------
      ! set the isf depth as it is in the initial step
      fsdept_n(:,:,1) = 0.5_wp * fse3w_n(:,:,1)
      fsdepw_n(:,:,1) = 0.0_wp
      fsde3w_n(:,:,1) = fsdept_n(:,:,1) - sshn(:,:)
      fsdept_b(:,:,1) = 0.5_wp * fse3w_b(:,:,1)
      fsdepw_b(:,:,1) = 0.0_wp

      DO jk = 2, jpk
         DO jj = 1,jpj
            DO ji = 1,jpi
              !    zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))   ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
                                                     ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
                                                     ! 0.5 where jk = mikt  
               zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))
               fsdepw_n(ji,jj,jk) = fsdepw_n(ji,jj,jk-1) + fse3t_n(ji,jj,jk-1)
               fsdept_n(ji,jj,jk) =      zcoef  * ( fsdepw_n(ji,jj,jk  ) + 0.5 * fse3w_n(ji,jj,jk))  &
                   &                + (1-zcoef) * ( fsdept_n(ji,jj,jk-1) +       fse3w_n(ji,jj,jk)) 
               fsde3w_n(ji,jj,jk) = fsdept_n(ji,jj,jk) - sshn(ji,jj)
               fsdepw_b(ji,jj,jk) = fsdepw_b(ji,jj,jk-1) + fse3t_b(ji,jj,jk-1)
               fsdept_b(ji,jj,jk) =      zcoef  * ( fsdepw_b(ji,jj,jk  ) + 0.5 * fse3w_b(ji,jj,jk))  &
                   &                + (1-zcoef) * ( fsdept_b(ji,jj,jk-1) +       fse3w_b(ji,jj,jk)) 
            END DO
         END DO
      END DO

      ! Before depth and Inverse of the local depth of the water column at u- and v- points
      ! -----------------------------------------------------------------------------------
      hu_b(:,:) = 0.
      hv_b(:,:) = 0.
      DO jk = 1, jpkm1
         hu_b(:,:) = hu_b(:,:) + fse3u_b(:,:,jk) * umask(:,:,jk)
         hv_b(:,:) = hv_b(:,:) + fse3v_b(:,:,jk) * vmask(:,:,jk)
      END DO
      hur_b(:,:) = umask_i(:,:) / ( hu_b(:,:) + 1. - umask_i(:,:) )
      hvr_b(:,:) = vmask_i(:,:) / ( hv_b(:,:) + 1. - vmask_i(:,:) )

      ! Restoring frequencies for z_tilde coordinate
      ! ============================================
      IF( ln_vvl_ztilde ) THEN
         ! Values in days provided via the namelist; use rsmall to avoid possible division by zero errors with faulty settings
         frq_rst_e3t(:,:) = 2.0_wp * rpi / ( MAX( rn_rst_e3t  , rsmall ) * 86400.0_wp )
         frq_rst_hdv(:,:) = 2.0_wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.0_wp )
         IF( ln_vvl_ztilde_as_zstar ) THEN
            ! Ignore namelist settings and use these next two to emulate z-star using z-tilde
            frq_rst_e3t(:,:) = 0.0_wp 
            frq_rst_hdv(:,:) = 1.0_wp / rdt
         ENDIF
         IF ( ln_vvl_zstar_at_eqtor ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( ABS(gphit(ji,jj)) >= 6.) THEN
                     ! values outside the equatorial band and transition zone (ztilde)
                     frq_rst_e3t(ji,jj) =  2.0_wp * rpi / ( MAX( rn_rst_e3t  , rsmall ) * 86400.e0_wp )
                     frq_rst_hdv(ji,jj) =  2.0_wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.e0_wp )
                  ELSEIF( ABS(gphit(ji,jj)) <= 2.5) THEN
                     ! values inside the equatorial band (ztilde as zstar)
                     frq_rst_e3t(ji,jj) =  0.0_wp
                     frq_rst_hdv(ji,jj) =  1.0_wp / rdt
                  ELSE
                     ! values in the transition band (linearly vary from ztilde to ztilde as zstar values)
                     frq_rst_e3t(ji,jj) = 0.0_wp + (frq_rst_e3t(ji,jj)-0.0_wp)*0.5_wp   &
                        &            * (  1.0_wp - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
                        &                                          * 180._wp / 3.5_wp ) )
                     frq_rst_hdv(ji,jj) = (1.0_wp / rdt)                                &
                        &            + (  frq_rst_hdv(ji,jj)-(1.e0_wp / rdt) )*0.5_wp   &
                        &            * (  1._wp  - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
                        &                                          * 180._wp / 3.5_wp ) )
                  ENDIF
               END DO
            END DO
            IF( cp_cfg == "orca" .AND. jp_cfg == 3 ) THEN
               ii0 = 103   ;   ii1 = 111        ! Suppress ztilde in the Foxe Basin for ORCA2
               ij0 = 128   ;   ij1 = 135   ;   
               frq_rst_e3t( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) =  0.0_wp
               frq_rst_hdv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) =  1.e0_wp / rdt
            ENDIF
         ENDIF
      ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('dom_vvl_init')

   END SUBROUTINE dom_vvl_init


   SUBROUTINE dom_vvl_sf_nxt( kt, kcall ) 
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_sf_nxt  ***
      !!                   
      !! ** Purpose :  - compute the after scale factors used in tra_zdf, dynnxt,
      !!                 tranxt and dynspg routines
      !!
      !! ** Method  :  - z_star case:  Repartition of ssh INCREMENT proportionnaly to the level thickness.
      !!               - z_tilde_case: after scale factor increment = 
      !!                                    high frequency part of horizontal divergence
      !!                                  + retsoring towards the background grid
      !!                                  + thickness difusion
      !!                               Then repartition of ssh INCREMENT proportionnaly
      !!                               to the "baroclinic" level thickness.
      !!
      !! ** Action  :  - hdiv_lf    : restoring towards full baroclinic divergence in z_tilde case
      !!               - tilde_e3t_a: after increment of vertical scale factor 
      !!                              in z_tilde case
      !!               - fse3(t/u/v)_a
      !!
      !! Reference  : Leclair, M., and Madec, G. 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ze3t
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zht, z_scale, zwu, zwv, zhdiv
      !! * Arguments
      INTEGER, INTENT( in )                  :: kt                    ! time step
      INTEGER, INTENT( in ), OPTIONAL        :: kcall                 ! optional argument indicating call sequence
      !! * Local declarations
      INTEGER                                :: ji, jj, jk            ! dummy loop indices
      INTEGER , DIMENSION(3)                 :: ijk_max, ijk_min      ! temporary integers
      REAL(wp)                               :: z2dt                  ! temporary scalars
      REAL(wp)                               :: z_tmin, z_tmax        ! temporary scalars
      LOGICAL                                :: ll_do_bclinic         ! temporary logical
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('dom_vvl_sf_nxt')
      CALL wrk_alloc( jpi, jpj, zht, z_scale, zwu, zwv, zhdiv )
      CALL wrk_alloc( jpi, jpj, jpk, ze3t                     )

      IF(kt == nit000)   THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_vvl_sf_nxt : compute after scale factors'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      ll_do_bclinic = .TRUE.
      IF( PRESENT(kcall) ) THEN
         IF ( kcall == 2 .AND. ln_vvl_ztilde ) ll_do_bclinic = .FALSE.
      ENDIF

      ! ******************************* !
      ! After acale factors at t-points !
      ! ******************************* !

      !                                                ! --------------------------------------------- !
                                                       ! z_star coordinate and barotropic z-tilde part !
      !                                                ! --------------------------------------------- !

      z_scale(:,:) = ( ssha(:,:) - sshb(:,:) ) * ssmask(:,:) / ( ht_0(:,:) + sshn(:,:) + 1. - ssmask(:,:) )
      DO jk = 1, jpkm1
         ! formally this is the same as fse3t_a = e3t_0*(1+ssha/ht_0)
         fse3t_a(:,:,jk) = fse3t_b(:,:,jk) + fse3t_n(:,:,jk) * z_scale(:,:) * tmask(:,:,jk)
      END DO

      IF( ln_vvl_ztilde .OR. ln_vvl_layer .AND. ll_do_bclinic ) THEN   ! z_tilde or layer coordinate !
         !                                                            ! ------baroclinic part------ !

         ! I - initialization
         ! ==================

         ! 1 - barotropic divergence
         ! -------------------------
         zhdiv(:,:) = 0.
         zht(:,:)   = 0.
         DO jk = 1, jpkm1
            zhdiv(:,:) = zhdiv(:,:) + fse3t_n(:,:,jk) * hdivn(:,:,jk)
            zht  (:,:) = zht  (:,:) + fse3t_n(:,:,jk) * tmask(:,:,jk)
         END DO
         zhdiv(:,:) = zhdiv(:,:) / ( zht(:,:) + 1. - tmask_i(:,:) )

         ! 2 - Low frequency baroclinic horizontal divergence  (z-tilde case only)
         ! --------------------------------------------------
         IF( ln_vvl_ztilde ) THEN
            IF( kt .GT. nit000 ) THEN
               DO jk = 1, jpkm1
                  hdiv_lf(:,:,jk) = hdiv_lf(:,:,jk) - rdt * frq_rst_hdv(:,:)   &
                     &          * ( hdiv_lf(:,:,jk) - fse3t_n(:,:,jk) * ( hdivn(:,:,jk) - zhdiv(:,:) ) )
               END DO
            ENDIF
         END IF

         ! II - after z_tilde increments of vertical scale factors
         ! =======================================================
         tilde_e3t_a(:,:,:) = 0.0_wp  ! tilde_e3t_a used to store tendency terms

         ! 1 - High frequency divergence term
         ! ----------------------------------
         IF( ln_vvl_ztilde ) THEN     ! z_tilde case
            DO jk = 1, jpkm1
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - ( fse3t_n(:,:,jk) * ( hdivn(:,:,jk) - zhdiv(:,:) ) - hdiv_lf(:,:,jk) )
            END DO
         ELSE                         ! layer case
            DO jk = 1, jpkm1
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) -   fse3t_n(:,:,jk) * ( hdivn(:,:,jk) - zhdiv(:,:) ) * tmask(:,:,jk)
            END DO
         END IF

         ! 2 - Restoring term (z-tilde case only)
         ! ------------------
         IF( ln_vvl_ztilde ) THEN
            DO jk = 1, jpk
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - frq_rst_e3t(:,:) * tilde_e3t_b(:,:,jk)
            END DO
         END IF

         ! 3 - Thickness diffusion term
         ! ----------------------------
         zwu(:,:) = 0.0_wp
         zwv(:,:) = 0.0_wp
         ! a - first derivative: diffusive fluxes
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  un_td(ji,jj,jk) = rn_ahe3 * umask(ji,jj,jk) * re2u_e1u(ji,jj) &
                                  & * ( tilde_e3t_b(ji,jj,jk) - tilde_e3t_b(ji+1,jj  ,jk) )
                  vn_td(ji,jj,jk) = rn_ahe3 * vmask(ji,jj,jk) * re1v_e2v(ji,jj) & 
                                  & * ( tilde_e3t_b(ji,jj,jk) - tilde_e3t_b(ji  ,jj+1,jk) )
                  zwu(ji,jj) = zwu(ji,jj) + un_td(ji,jj,jk)
                  zwv(ji,jj) = zwv(ji,jj) + vn_td(ji,jj,jk)
               END DO
            END DO
         END DO
         ! b - correction for last oceanic u-v points
         DO jj = 1, jpj
            DO ji = 1, jpi
               un_td(ji,jj,mbku(ji,jj)) = un_td(ji,jj,mbku(ji,jj)) - zwu(ji,jj)
               vn_td(ji,jj,mbkv(ji,jj)) = vn_td(ji,jj,mbkv(ji,jj)) - zwv(ji,jj)
            END DO
         END DO
         ! c - second derivative: divergence of diffusive fluxes
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  tilde_e3t_a(ji,jj,jk) = tilde_e3t_a(ji,jj,jk) + (   un_td(ji-1,jj  ,jk) - un_td(ji,jj,jk)    &
                     &                                          +     vn_td(ji  ,jj-1,jk) - vn_td(ji,jj,jk)    &
                     &                                            ) * r1_e12t(ji,jj)
               END DO
            END DO
         END DO
         ! d - thickness diffusion transport: boundary conditions
         !     (stored for tracer advction and continuity equation)
         CALL lbc_lnk( un_td , 'U' , -1._wp)
         CALL lbc_lnk( vn_td , 'V' , -1._wp)

         ! 4 - Time stepping of baroclinic scale factors
         ! ---------------------------------------------
         ! Leapfrog time stepping
         ! ~~~~~~~~~~~~~~~~~~~~~~
         IF( neuler == 0 .AND. kt == nit000 ) THEN
            z2dt =  rdt
         ELSE
            z2dt = 2.0_wp * rdt
         ENDIF
         CALL lbc_lnk( tilde_e3t_a(:,:,:), 'T', 1._wp )
         tilde_e3t_a(:,:,:) = tilde_e3t_b(:,:,:) + z2dt * tmask(:,:,:) * tilde_e3t_a(:,:,:)

         ! Maximum deformation control
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ze3t(:,:,jpk) = 0.0_wp
         DO jk = 1, jpkm1
            ze3t(:,:,jk) = tilde_e3t_a(:,:,jk) / e3t_0(:,:,jk) * tmask(:,:,jk) * tmask_i(:,:)
         END DO
         z_tmax = MAXVAL( ze3t(:,:,:) )
         IF( lk_mpp )   CALL mpp_max( z_tmax )                 ! max over the global domain
         z_tmin = MINVAL( ze3t(:,:,:) )
         IF( lk_mpp )   CALL mpp_min( z_tmin )                 ! min over the global domain
         ! - ML - test: for the moment, stop simulation for too large e3_t variations
         IF( ( z_tmax .GT. rn_zdef_max ) .OR. ( z_tmin .LT. - rn_zdef_max ) ) THEN
            IF( lk_mpp ) THEN
               CALL mpp_maxloc( ze3t, tmask, z_tmax, ijk_max(1), ijk_max(2), ijk_max(3) )
               CALL mpp_minloc( ze3t, tmask, z_tmin, ijk_min(1), ijk_min(2), ijk_min(3) )
            ELSE
               ijk_max = MAXLOC( ze3t(:,:,:) )
               ijk_max(1) = ijk_max(1) + nimpp - 1
               ijk_max(2) = ijk_max(2) + njmpp - 1
               ijk_min = MINLOC( ze3t(:,:,:) )
               ijk_min(1) = ijk_min(1) + nimpp - 1
               ijk_min(2) = ijk_min(2) + njmpp - 1
            ENDIF
            IF (lwp) THEN
               WRITE(numout, *) 'MAX( tilde_e3t_a(:,:,:) / e3t_0(:,:,:) ) =', z_tmax
               WRITE(numout, *) 'at i, j, k=', ijk_max
               WRITE(numout, *) 'MIN( tilde_e3t_a(:,:,:) / e3t_0(:,:,:) ) =', z_tmin
               WRITE(numout, *) 'at i, j, k=', ijk_min            
               CALL ctl_warn('MAX( ABS( tilde_e3t_a(:,:,:) ) / e3t_0(:,:,:) ) too high')
            ENDIF
         ENDIF
         ! - ML - end test
         ! - ML - Imposing these limits will cause a baroclinicity error which is corrected for below
         tilde_e3t_a(:,:,:) = MIN( tilde_e3t_a(:,:,:),   rn_zdef_max * e3t_0(:,:,:) )
         tilde_e3t_a(:,:,:) = MAX( tilde_e3t_a(:,:,:), - rn_zdef_max * e3t_0(:,:,:) )

         !
         ! "tilda" change in the after scale factor
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         DO jk = 1, jpkm1
            dtilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - tilde_e3t_b(:,:,jk)
         END DO
         ! III - Barotropic repartition of the sea surface height over the baroclinic profile
         ! ==================================================================================
         ! add ( ssh increment + "baroclinicity error" ) proportionly to e3t(n)
         ! - ML - baroclinicity error should be better treated in the future
         !        i.e. locally and not spread over the water column.
         !        (keep in mind that the idea is to reduce Eulerian velocity as much as possible)
         zht(:,:) = 0.
         DO jk = 1, jpkm1
            zht(:,:)  = zht(:,:) + tilde_e3t_a(:,:,jk) * tmask(:,:,jk)
         END DO
         z_scale(:,:) =  - zht(:,:) / ( ht_0(:,:) + sshn(:,:) + 1. - ssmask(:,:) )
         DO jk = 1, jpkm1
            dtilde_e3t_a(:,:,jk) = dtilde_e3t_a(:,:,jk) + fse3t_n(:,:,jk) * z_scale(:,:) * tmask(:,:,jk)
         END DO

      ENDIF

      IF( ln_vvl_ztilde .OR. ln_vvl_layer )  THEN   ! z_tilde or layer coordinate !
      !                                           ! ---baroclinic part--------- !
         DO jk = 1, jpkm1
            fse3t_a(:,:,jk) = fse3t_a(:,:,jk) + dtilde_e3t_a(:,:,jk) * tmask(:,:,jk)
         END DO
      ENDIF

      IF( ln_vvl_dbg .AND. .NOT. ll_do_bclinic ) THEN   ! - ML - test: control prints for debuging
         !
         IF( lwp ) WRITE(numout, *) 'kt =', kt
         IF ( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
            z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( zht(:,:) ) )
            IF( lk_mpp ) CALL mpp_max( z_tmax )                             ! max over the global domain
            IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(SUM(tilde_e3t_a))) =', z_tmax
         END IF
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + fse3t_n(:,:,jk) * tmask(:,:,jk)
         END DO
         z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( ht_0(:,:) + sshn(:,:) - zht(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ht_0+sshn-SUM(fse3t_n))) =', z_tmax
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + fse3t_a(:,:,jk) * tmask(:,:,jk)
         END DO
         z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( ht_0(:,:) + ssha(:,:) - zht(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ht_0+ssha-SUM(fse3t_a))) =', z_tmax
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + fse3t_b(:,:,jk) * tmask(:,:,jk)
         END DO
         z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( ht_0(:,:) + sshb(:,:) - zht(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ht_0+sshb-SUM(fse3t_b))) =', z_tmax
         !
         z_tmax = MAXVAL( tmask(:,:,1) *  ABS( sshb(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(sshb))) =', z_tmax
         !
         z_tmax = MAXVAL( tmask(:,:,1) *  ABS( sshn(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(sshn))) =', z_tmax
         !
         z_tmax = MAXVAL( tmask(:,:,1) *  ABS( ssha(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ssha))) =', z_tmax
      END IF

      ! *********************************** !
      ! After scale factors at u- v- points !
      ! *********************************** !

      CALL dom_vvl_interpol( fse3t_a(:,:,:), fse3u_a(:,:,:), 'U' )
      CALL dom_vvl_interpol( fse3t_a(:,:,:), fse3v_a(:,:,:), 'V' )

      ! *********************************** !
      ! After depths at u- v points         !
      ! *********************************** !

      hu_a(:,:) = 0._wp                        ! Ocean depth at U-points
      hv_a(:,:) = 0._wp                        ! Ocean depth at V-points
      DO jk = 1, jpkm1
         hu_a(:,:) = hu_a(:,:) + fse3u_a(:,:,jk) * umask(:,:,jk)
         hv_a(:,:) = hv_a(:,:) + fse3v_a(:,:,jk) * vmask(:,:,jk)
      END DO
      !                                        ! Inverse of the local depth
      hur_a(:,:) = 1._wp / ( hu_a(:,:) + 1._wp - umask_i(:,:) ) * umask_i(:,:)
      hvr_a(:,:) = 1._wp / ( hv_a(:,:) + 1._wp - vmask_i(:,:) ) * vmask_i(:,:)

      CALL wrk_dealloc( jpi, jpj, zht, z_scale, zwu, zwv, zhdiv )
      CALL wrk_dealloc( jpi, jpj, jpk, ze3t                     )

      IF( nn_timing == 1 )  CALL timing_stop('dom_vvl_sf_nxt')

   END SUBROUTINE dom_vvl_sf_nxt


   SUBROUTINE dom_vvl_sf_swp( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_sf_swp  ***
      !!                   
      !! ** Purpose :  compute time filter and swap of scale factors 
      !!               compute all depths and related variables for next time step
      !!               write outputs and restart file
      !!
      !! ** Method  :  - swap of e3t with trick for volume/tracer conservation
      !!               - reconstruct scale factor at other grid points (interpolate)
      !!               - recompute depths and water height fields
      !!
      !! ** Action  :  - fse3t_(b/n), tilde_e3t_(b/n) and fse3(u/v)_n ready for next time step
      !!               - Recompute:
      !!                    fse3(u/v)_b       
      !!                    fse3w_n           
      !!                    fse3(u/v)w_b      
      !!                    fse3(u/v)w_n      
      !!                    fsdept_n, fsdepw_n  and fsde3w_n
      !!                    h(u/v) and h(u/v)r
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!              Leclair, M., and G. Madec, 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in )               :: kt       ! time step
      !! * Local declarations
      INTEGER                             :: ji,jj,jk       ! dummy loop indices
      REAL(wp)                            :: zcoef
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('dom_vvl_sf_swp')
      !
      IF( kt == nit000 )   THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_vvl_sf_swp : - time filter and swap of scale factors'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   - interpolate scale factors and compute depths for next time step'
      ENDIF

      !
      ! Time filter and swap of scale factors
      ! =====================================
      ! - ML - fse3(t/u/v)_b are allready computed in dynnxt.
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         IF( neuler == 0 .AND. kt == nit000 ) THEN
            tilde_e3t_b(:,:,:) = tilde_e3t_n(:,:,:)
         ELSE
            tilde_e3t_b(:,:,:) = tilde_e3t_n(:,:,:) & 
            &         + atfp * ( tilde_e3t_b(:,:,:) - 2.0_wp * tilde_e3t_n(:,:,:) + tilde_e3t_a(:,:,:) )
         ENDIF
         tilde_e3t_n(:,:,:) = tilde_e3t_a(:,:,:)
      ENDIF
      fsdept_b(:,:,:) = fsdept_n(:,:,:)
      fsdepw_b(:,:,:) = fsdepw_n(:,:,:)

      fse3t_n(:,:,:) = fse3t_a(:,:,:)
      fse3u_n(:,:,:) = fse3u_a(:,:,:)
      fse3v_n(:,:,:) = fse3v_a(:,:,:)

      ! Compute all missing vertical scale factor and depths
      ! ====================================================
      ! Horizontal scale factor interpolations
      ! --------------------------------------
      ! - ML - fse3u_b and fse3v_b are allready computed in dynnxt
      ! - JC - hu_b, hv_b, hur_b, hvr_b also
      CALL dom_vvl_interpol( fse3u_n(:,:,:), fse3f_n (:,:,:), 'F'  )
      ! Vertical scale factor interpolations
      ! ------------------------------------
      CALL dom_vvl_interpol( fse3t_n(:,:,:), fse3w_n (:,:,:), 'W'  )
      CALL dom_vvl_interpol( fse3u_n(:,:,:), fse3uw_n(:,:,:), 'UW' )
      CALL dom_vvl_interpol( fse3v_n(:,:,:), fse3vw_n(:,:,:), 'VW' )
      CALL dom_vvl_interpol( fse3t_b(:,:,:), fse3w_b (:,:,:), 'W'  )
      CALL dom_vvl_interpol( fse3u_b(:,:,:), fse3uw_b(:,:,:), 'UW' )
      CALL dom_vvl_interpol( fse3v_b(:,:,:), fse3vw_b(:,:,:), 'VW' )
      ! t- and w- points depth
      ! ----------------------
      ! set the isf depth as it is in the initial step
      fsdept_n(:,:,1) = 0.5_wp * fse3w_n(:,:,1)
      fsdepw_n(:,:,1) = 0.0_wp
      fsde3w_n(:,:,1) = fsdept_n(:,:,1) - sshn(:,:)

      DO jk = 2, jpk
         DO jj = 1,jpj
            DO ji = 1,jpi
              !    zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))   ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
                                                                 ! 1 for jk = mikt
               zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))
               fsdepw_n(ji,jj,jk) = fsdepw_n(ji,jj,jk-1) + fse3t_n(ji,jj,jk-1)
               fsdept_n(ji,jj,jk) =      zcoef  * ( fsdepw_n(ji,jj,jk  ) + 0.5 * fse3w_n(ji,jj,jk))  &
                   &                + (1-zcoef) * ( fsdept_n(ji,jj,jk-1) +       fse3w_n(ji,jj,jk)) 
               fsde3w_n(ji,jj,jk) = fsdept_n(ji,jj,jk) - sshn(ji,jj)
            END DO
         END DO
      END DO

      ! Local depth and Inverse of the local depth of the water column at u- and v- points
      ! ----------------------------------------------------------------------------------
      hu (:,:) = hu_a (:,:)
      hv (:,:) = hv_a (:,:)

      ! Inverse of the local depth
      hur(:,:) = hur_a(:,:)
      hvr(:,:) = hvr_a(:,:)

      ! Local depth of the water column at t- points
      ! --------------------------------------------
      ht(:,:) = 0.
      DO jk = 1, jpkm1
         ht(:,:) = ht(:,:) + fse3t_n(:,:,jk) * tmask(:,:,jk)
      END DO

      ! write restart file
      ! ==================
      IF( lrst_oce ) CALL dom_vvl_rst( kt, 'WRITE' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dom_vvl_sf_swp')

   END SUBROUTINE dom_vvl_sf_swp


   SUBROUTINE dom_vvl_interpol( pe3_in, pe3_out, pout )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl__interpol  ***
      !!
      !! ** Purpose :   interpolate scale factors from one grid point to another
      !!
      !! ** Method  :   e3_out = e3_0 + interpolation(e3_in - e3_0)
      !!                - horizontal interpolation: grid cell surface averaging
      !!                - vertical interpolation: simple averaging
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in    ) ::  pe3_in     ! input e3 to be interpolated
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::  pe3_out    ! output interpolated e3
      CHARACTER(LEN=*), INTENT( in )                    ::  pout       ! grid point of out scale factors
      !                                                                !   =  'U', 'V', 'W, 'F', 'UW' or 'VW'
      !! * Local declarations
      INTEGER ::   ji, jj, jk                                          ! dummy loop indices
      LOGICAL ::   l_is_orca                                           ! local logical
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('dom_vvl_interpol')
         !
      l_is_orca = .FALSE.
      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) l_is_orca = .TRUE.      ! ORCA R2 configuration - will need to correct some locations

      SELECT CASE ( pout )
         !               ! ------------------------------------- !
      CASE( 'U' )        ! interpolation from T-point to U-point !
         !               ! ------------------------------------- !
         ! horizontal surface weighted interpolation
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * r1_e12u(ji,jj)                                   &
                     &                       * (   e12t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )     &
                     &                           + e12t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
               END DO
            END DO
         END DO
         !
         IF( l_is_orca ) CALL dom_vvl_orca_fix( pe3_in, pe3_out, pout )
         ! boundary conditions
         CALL lbc_lnk( pe3_out(:,:,:), 'U', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3u_0(:,:,:)
         !               ! ------------------------------------- !
      CASE( 'V' )        ! interpolation from T-point to V-point !
         !               ! ------------------------------------- !
         ! horizontal surface weighted interpolation
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk) * r1_e12v(ji,jj)                                   &
                     &                       * (   e12t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )     &
                     &                           + e12t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
               END DO
            END DO
         END DO
         !
         IF( l_is_orca ) CALL dom_vvl_orca_fix( pe3_in, pe3_out, pout )
         ! boundary conditions
         CALL lbc_lnk( pe3_out(:,:,:), 'V', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3v_0(:,:,:)
         !               ! ------------------------------------- !
      CASE( 'F' )        ! interpolation from U-point to F-point !
         !               ! ------------------------------------- !
         ! horizontal surface weighted interpolation
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk) * r1_e12f(ji,jj)               &
                     &                       * (   e12u(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3u_0(ji,jj  ,jk) )     &
                     &                           + e12u(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3u_0(ji,jj+1,jk) ) )
               END DO
            END DO
         END DO
         !
         IF( l_is_orca ) CALL dom_vvl_orca_fix( pe3_in, pe3_out, pout )
         ! boundary conditions
         CALL lbc_lnk( pe3_out(:,:,:), 'F', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3f_0(:,:,:)
         !               ! ------------------------------------- !
      CASE( 'W' )        ! interpolation from T-point to W-point !
         !               ! ------------------------------------- !
         ! vertical simple interpolation
         pe3_out(:,:,1) = e3w_0(:,:,1) + pe3_in(:,:,1) - e3t_0(:,:,1)
         ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3w_0(:,:,jk) + ( 1.0_wp - 0.5_wp * tmask(:,:,jk) ) * ( pe3_in(:,:,jk-1) - e3t_0(:,:,jk-1) )   &
               &                            +            0.5_wp * tmask(:,:,jk)   * ( pe3_in(:,:,jk  ) - e3t_0(:,:,jk  ) )
         END DO
         !               ! -------------------------------------- !
      CASE( 'UW' )       ! interpolation from U-point to UW-point !
         !               ! -------------------------------------- !
         ! vertical simple interpolation
         pe3_out(:,:,1) = e3uw_0(:,:,1) + pe3_in(:,:,1) - e3u_0(:,:,1)
         ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3uw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * umask(:,:,jk) ) * ( pe3_in(:,:,jk-1) - e3u_0(:,:,jk-1) )   &
               &                             +            0.5_wp * umask(:,:,jk)   * ( pe3_in(:,:,jk  ) - e3u_0(:,:,jk  ) )
         END DO
         !               ! -------------------------------------- !
      CASE( 'VW' )       ! interpolation from V-point to VW-point !
         !               ! -------------------------------------- !
         ! vertical simple interpolation
         pe3_out(:,:,1) = e3vw_0(:,:,1) + pe3_in(:,:,1) - e3v_0(:,:,1)
         ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3vw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * vmask(:,:,jk) ) * ( pe3_in(:,:,jk-1) - e3v_0(:,:,jk-1) )   &
               &                             +            0.5_wp * vmask(:,:,jk)   * ( pe3_in(:,:,jk  ) - e3v_0(:,:,jk  ) )
         END DO
      END SELECT
      !

      IF( nn_timing == 1 )  CALL timing_stop('dom_vvl_interpol')

   END SUBROUTINE dom_vvl_interpol

   SUBROUTINE dom_vvl_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dom_vvl_rst  ***
      !!                     
      !! ** Purpose :   Read or write VVL file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!                if the restart does not contain vertical scale factors,
      !!                they are set to the _0 values
      !!                if the restart does not contain vertical scale factors increments (z_tilde),
      !!                they are set to 0.
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !! * Local declarations
      INTEGER ::   jk
      INTEGER ::   id1, id2, id3, id4, id5     ! local integers
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dom_vvl_rst')
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ===============
         IF( ln_rstart ) THEN                   !* Read the restart file
            CALL rst_read_open                  !  open the restart file if necessary
            CALL iom_get( numror, jpdom_autoglo, 'sshn'   , sshn    )
            !
            id1 = iom_varid( numror, 'fse3t_b', ldstop = .FALSE. )
            id2 = iom_varid( numror, 'fse3t_n', ldstop = .FALSE. )
            id3 = iom_varid( numror, 'tilde_e3t_b', ldstop = .FALSE. )
            id4 = iom_varid( numror, 'tilde_e3t_n', ldstop = .FALSE. )
            id5 = iom_varid( numror, 'hdiv_lf', ldstop = .FALSE. )
            !                             ! --------- !
            !                             ! all cases !
            !                             ! --------- !
            IF( MIN( id1, id2 ) > 0 ) THEN       ! all required arrays exist
               CALL iom_get( numror, jpdom_autoglo, 'fse3t_b', fse3t_b(:,:,:) )
               CALL iom_get( numror, jpdom_autoglo, 'fse3t_n', fse3t_n(:,:,:) )
               ! needed to restart if land processor not computed 
               IF(lwp) write(numout,*) 'dom_vvl_rst : fse3t_b and fse3t_n found in restart files'
               WHERE ( tmask(:,:,:) == 0.0_wp ) 
                  fse3t_n(:,:,:) = e3t_0(:,:,:)
                  fse3t_b(:,:,:) = e3t_0(:,:,:)
               END WHERE
               IF( neuler == 0 ) THEN
                  fse3t_b(:,:,:) = fse3t_n(:,:,:)
               ENDIF
            ELSE IF( id1 > 0 ) THEN
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : fse3t_n not found in restart files'
               IF(lwp) write(numout,*) 'fse3t_n set equal to fse3t_b.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'fse3t_b', fse3t_b(:,:,:) )
               fse3t_n(:,:,:) = fse3t_b(:,:,:)
               neuler = 0
            ELSE IF( id2 > 0 ) THEN
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : fse3t_b not found in restart files'
               IF(lwp) write(numout,*) 'fse3t_b set equal to fse3t_n.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'fse3t_n', fse3t_n(:,:,:) )
               fse3t_b(:,:,:) = fse3t_n(:,:,:)
               neuler = 0
            ELSE
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : fse3t_n not found in restart file'
               IF(lwp) write(numout,*) 'Compute scale factor from sshn'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               DO jk=1,jpk
                  fse3t_n(:,:,jk) =  e3t_0(:,:,jk) * ( ht_0(:,:) + sshn(:,:) ) &
                      &                            / ( ht_0(:,:) + 1._wp - ssmask(:,:) ) * tmask(:,:,jk) &
                      &            + e3t_0(:,:,jk) * (1._wp -tmask(:,:,jk))
               END DO
               fse3t_b(:,:,:) = fse3t_n(:,:,:)
               neuler = 0
            ENDIF
            !                             ! ----------- !
            IF( ln_vvl_zstar ) THEN       ! z_star case !
               !                          ! ----------- !
               IF( MIN( id3, id4 ) > 0 ) THEN
                  CALL ctl_stop( 'dom_vvl_rst: z_star cannot restart from a z_tilde or layer run' )
               ENDIF
               !                          ! ----------------------- !
            ELSE                          ! z_tilde and layer cases !
               !                          ! ----------------------- !
               IF( MIN( id3, id4 ) > 0 ) THEN  ! all required arrays exist
                  CALL iom_get( numror, jpdom_autoglo, 'tilde_e3t_b', tilde_e3t_b(:,:,:) )
                  CALL iom_get( numror, jpdom_autoglo, 'tilde_e3t_n', tilde_e3t_n(:,:,:) )
               ELSE                            ! one at least array is missing
                  tilde_e3t_b(:,:,:) = 0.0_wp
                  tilde_e3t_n(:,:,:) = 0.0_wp
               ENDIF
               !                          ! ------------ !
               IF( ln_vvl_ztilde ) THEN   ! z_tilde case !
                  !                       ! ------------ !
                  IF( id5 > 0 ) THEN  ! required array exists
                     CALL iom_get( numror, jpdom_autoglo, 'hdiv_lf', hdiv_lf(:,:,:) )
                  ELSE                ! array is missing
                     hdiv_lf(:,:,:) = 0.0_wp
                  ENDIF
               ENDIF
            ENDIF
            !
         ELSE                                   !* Initialize at "rest"
            fse3t_b(:,:,:) = e3t_0(:,:,:)
            fse3t_n(:,:,:) = e3t_0(:,:,:)
            sshn(:,:) = 0.0_wp
            IF( ln_vvl_ztilde .OR. ln_vvl_layer) THEN
               tilde_e3t_b(:,:,:) = 0.0_wp
               tilde_e3t_n(:,:,:) = 0.0_wp
               IF( ln_vvl_ztilde ) hdiv_lf(:,:,:) = 0.0_wp
            END IF
         ENDIF

      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! ===================
         IF(lwp) WRITE(numout,*) '---- dom_vvl_rst ----'
         !                                           ! --------- !
         !                                           ! all cases !
         !                                           ! --------- !
         CALL iom_rstput( kt, nitrst, numrow, 'fse3t_b', fse3t_b(:,:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'fse3t_n', fse3t_n(:,:,:) )
         !                                           ! ----------------------- !
         IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN  ! z_tilde and layer cases !
            !                                        ! ----------------------- !
            CALL iom_rstput( kt, nitrst, numrow, 'tilde_e3t_b', tilde_e3t_b(:,:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'tilde_e3t_n', tilde_e3t_n(:,:,:) )
         END IF
         !                                           ! -------------!    
         IF( ln_vvl_ztilde ) THEN                    ! z_tilde case !
            !                                        ! ------------ !
            CALL iom_rstput( kt, nitrst, numrow, 'hdiv_lf', hdiv_lf(:,:,:) )
         ENDIF

      ENDIF
      IF( nn_timing == 1 )  CALL timing_stop('dom_vvl_rst')

   END SUBROUTINE dom_vvl_rst


   SUBROUTINE dom_vvl_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl_ctl  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options
      !!                for vertical coordinate
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio
      INTEGER ::   ios

      NAMELIST/nam_vvl/ ln_vvl_zstar, ln_vvl_ztilde, ln_vvl_layer, ln_vvl_ztilde_as_zstar, &
                      & ln_vvl_zstar_at_eqtor      , rn_ahe3     , rn_rst_e3t            , &
                      & rn_lf_cutoff               , rn_zdef_max , ln_vvl_dbg                ! not yet implemented: ln_vvl_kepe
      !!---------------------------------------------------------------------- 

      REWIND( numnam_ref )              ! Namelist nam_vvl in reference namelist : 
      READ  ( numnam_ref, nam_vvl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_vvl in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_vvl in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, nam_vvl, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_vvl in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_vvl )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_vvl_ctl : choice/control of the variable vertical coordinate'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '           Namelist nam_vvl : chose a vertical coordinate'
         WRITE(numout,*) '              zstar                      ln_vvl_zstar   = ', ln_vvl_zstar
         WRITE(numout,*) '              ztilde                     ln_vvl_ztilde  = ', ln_vvl_ztilde
         WRITE(numout,*) '              layer                      ln_vvl_layer   = ', ln_vvl_layer
         WRITE(numout,*) '              ztilde as zstar   ln_vvl_ztilde_as_zstar  = ', ln_vvl_ztilde_as_zstar
         WRITE(numout,*) '      ztilde near the equator    ln_vvl_zstar_at_eqtor  = ', ln_vvl_zstar_at_eqtor
         ! WRITE(numout,*) '           Namelist nam_vvl : chose kinetic-to-potential energy conservation'
         ! WRITE(numout,*) '                                         ln_vvl_kepe    = ', ln_vvl_kepe
         WRITE(numout,*) '           Namelist nam_vvl : thickness diffusion coefficient'
         WRITE(numout,*) '                                         rn_ahe3        = ', rn_ahe3
         WRITE(numout,*) '           Namelist nam_vvl : maximum e3t deformation fractional change'
         WRITE(numout,*) '                                         rn_zdef_max    = ', rn_zdef_max
         IF( ln_vvl_ztilde_as_zstar ) THEN
            WRITE(numout,*) '           ztilde running in zstar emulation mode; '
            WRITE(numout,*) '           ignoring namelist timescale parameters and using:'
            WRITE(numout,*) '                 hard-wired : z-tilde to zstar restoration timescale (days)'
            WRITE(numout,*) '                                         rn_rst_e3t     =    0.0'
            WRITE(numout,*) '                 hard-wired : z-tilde cutoff frequency of low-pass filter (days)'
            WRITE(numout,*) '                                         rn_lf_cutoff   =    1.0/rdt'
         ELSE
            WRITE(numout,*) '           Namelist nam_vvl : z-tilde to zstar restoration timescale (days)'
            WRITE(numout,*) '                                         rn_rst_e3t     = ', rn_rst_e3t
            WRITE(numout,*) '           Namelist nam_vvl : z-tilde cutoff frequency of low-pass filter (days)'
            WRITE(numout,*) '                                         rn_lf_cutoff   = ', rn_lf_cutoff
         ENDIF
         WRITE(numout,*) '           Namelist nam_vvl : debug prints'
         WRITE(numout,*) '                                         ln_vvl_dbg     = ', ln_vvl_dbg
      ENDIF

      ioptio = 0                      ! Parameter control
      IF( ln_vvl_ztilde_as_zstar ) ln_vvl_ztilde = .true.
      IF( ln_vvl_zstar           )        ioptio = ioptio + 1
      IF( ln_vvl_ztilde          )        ioptio = ioptio + 1
      IF( ln_vvl_layer           )        ioptio = ioptio + 1

      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE vertical coordinate in namelist nam_vvl' )
      IF( .NOT. ln_vvl_zstar .AND. nn_isf .NE. 0) CALL ctl_stop( 'Only vvl_zstar has been tested with ice shelf cavity' )

      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( ln_vvl_zstar           ) WRITE(numout,*) '              zstar vertical coordinate is used'
         IF( ln_vvl_ztilde          ) WRITE(numout,*) '              ztilde vertical coordinate is used'
         IF( ln_vvl_layer           ) WRITE(numout,*) '              layer vertical coordinate is used'
         IF( ln_vvl_ztilde_as_zstar ) WRITE(numout,*) '              to emulate a zstar coordinate'
         ! - ML - Option not developed yet
         ! IF(       ln_vvl_kepe ) WRITE(numout,*) '              kinetic to potential energy transfer : option used'
         ! IF( .NOT. ln_vvl_kepe ) WRITE(numout,*) '              kinetic to potential energy transfer : option not used'
      ENDIF

#if defined key_agrif
      IF (.NOT.Agrif_Root()) CALL ctl_stop( 'AGRIF not implemented with non-linear free surface (key_vvl)' )
#endif

   END SUBROUTINE dom_vvl_ctl

   SUBROUTINE dom_vvl_orca_fix( pe3_in, pe3_out, pout )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dom_vvl_orca_fix  ***
      !!                     
      !! ** Purpose :   Correct surface weighted, horizontally interpolated, 
      !!                scale factors at locations that have been individually
      !!                modified in domhgr. Such modifications break the
      !!                relationship between e12t and e1u*e2u etc.
      !!                Recompute some scale factors ignoring the modified metric.
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in    ) ::  pe3_in     ! input e3 to be interpolated
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::  pe3_out    ! output interpolated e3
      CHARACTER(LEN=*), INTENT( in )                    ::  pout       ! grid point of out scale factors
      !                                                                !   =  'U', 'V', 'W, 'F', 'UW' or 'VW'
      !! * Local declarations
      INTEGER ::   ji, jj, jk                                          ! dummy loop indices
      INTEGER ::   ij0, ij1, ii0, ii1                                  ! dummy loop indices
      INTEGER ::   isrow                                               ! index for ORCA1 starting row
      !! acc
      !! Hmm with the time splitting these "fixes" seem to do more harm than good. Temporarily disabled for
      !! the ORCA2 tests (by changing jp_cfg test from 2 to 3) pending further investigations
      !! 
      !                                                ! =====================
      IF( cp_cfg == "orca" .AND. jp_cfg == 3 ) THEN    ! ORCA R2 configuration
         !                                             ! =====================
      !! acc
         IF( nn_cla == 0 ) THEN
            !
            ii0 = 139   ;   ii1 = 140        ! Gibraltar Strait (e2u was modified)
            ij0 = 102   ;   ij1 = 102
            DO jk = 1, jpkm1
               DO jj = mj0(ij0), mj1(ij1)
                  DO ji = mi0(ii0), mi1(ii1)
                     SELECT CASE ( pout )
                     CASE( 'U' )
                        pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                       &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                       &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                       &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                     CASE( 'F' )
                        pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                       &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                       &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                       &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                     END SELECT
                  END DO
               END DO
            END DO
            !
            ii0 = 160   ;   ii1 = 160        ! Bab el Mandeb (e2u and e1v were modified)
            ij0 =  88   ;   ij1 =  88
            DO jk = 1, jpkm1
               DO jj = mj0(ij0), mj1(ij1)
                  DO ji = mi0(ii0), mi1(ii1)
                     SELECT CASE ( pout )
                     CASE( 'U' )
                        pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                       &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                       &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                       &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                     CASE( 'V' )
                        pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                       &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                       &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                       &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                     CASE( 'F' )
                        pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                       &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                       &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                       &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                     END SELECT
                  END DO
               END DO
            END DO
         ENDIF

         ii0 = 145   ;   ii1 = 146        ! Danish Straits (e2u was modified)
         ij0 = 116   ;   ij1 = 116
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
      ENDIF
      !
         !                                             ! =====================
      IF( cp_cfg == "orca" .AND. jp_cfg == 1 ) THEN    ! ORCA R1 configuration
         !                                             ! =====================
         ! This dirty section will be suppressed by simplification process:
         ! all this will come back in input files
         ! Currently these hard-wired indices relate to configuration with
         ! extend grid (jpjglo=332)
         ! which had a grid-size of 362x292.
         isrow = 332 - jpjglo
         !
         ii0 = 282           ;   ii1 = 283        ! Gibraltar Strait (e2u was modified)
         ij0 = 241 - isrow   ;   ij1 = 241 - isrow
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 = 314           ;   ii1 = 315        ! Bhosporus Strait (e2u was modified)
         ij0 = 248 - isrow   ;   ij1 = 248 - isrow
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &  
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &  
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  44           ;   ii1 =  44        ! Lombok Strait (e1v was modified)
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  48           ;   ii1 =  48        ! Sumba Strait (e1v was modified) [closed from bathy_11 on]
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  53          ;   ii1 =  53        ! Ombai Strait (e1v was modified)
         ij0 = 164 - isrow  ;   ij1 = 165  - isrow  
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  56            ;   ii1 =  56        ! Timor Passage (e1v was modified)
         ij0 = 164 - isrow    ;   ij1 = 165  - isrow  
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  55            ;   ii1 =  55        ! West Halmahera Strait (e1v was modified)
         ij0 = 181 - isrow    ;   ij1 = 182 - isrow  
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  58            ;   ii1 =  58        ! East Halmahera Strait (e1v was modified)
         ij0 = 181 - isrow    ;   ij1 = 182 - isrow  
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
      ENDIF
         !                                             ! =====================
      IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN   ! ORCA R05 configuration
         !                                             ! =====================
         !
         ii0 = 563   ;   ii1 = 564        ! Gibraltar Strait (e2u was modified)
         ij0 = 327   ;   ij1 = 327
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 = 627   ;   ii1 = 628        ! Bosphorus Strait (e2u was modified)
         ij0 = 343   ;   ij1 = 343
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &  
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &  
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  93   ;   ii1 =  94        ! Sumba Strait (e2u was modified)
         ij0 = 232   ;   ij1 = 232
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 = 103   ;   ii1 = 103        ! Ombai Strait (e2u was modified)
         ij0 = 232   ;   ij1 = 232
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  15   ;   ii1 =  15        ! Palk Strait (e2u was modified)
         ij0 = 270   ;   ij1 = 270
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'U' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk)                                        &
                    &                    * (   e1t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) ) &
                    &                    +     e1t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) &
                    &                      ) / e1u(ji,jj)   +   e3u_0(ji,jj,jk)
                  CASE( 'F' )
                     pe3_out(ji,jj,jk) = 0.5_wp * umask(ji,jj,jk) * umask(ji,jj+1,jk)                    &
                    &                    * (   e1u(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3u_0(ji  ,jj,jk) ) &
                    &                    +     e1u(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3u_0(ji+1,jj,jk) ) &
                    &                      ) / e1f(ji,jj)   +   e3f_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 =  87   ;   ii1 =  87        ! Lombok Strait (e1v was modified)
         ij0 = 232   ;   ij1 = 233
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
         !
         ii0 = 662   ;   ii1 = 662        ! Bab el Mandeb (e1v was modified)
         ij0 = 276   ;   ij1 = 276
         DO jk = 1, jpkm1
            DO jj = mj0(ij0), mj1(ij1)
               DO ji = mi0(ii0), mi1(ii1)
                  SELECT CASE ( pout )
                  CASE( 'V' )
                     pe3_out(ji,jj,jk) = 0.5_wp * vmask(ji,jj,jk)                                        &
                    &                    * (   e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) ) &
                    &                    +     e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) &
                    &                      ) / e2v(ji,jj)   +   e3v_0(ji,jj,jk)
                  END SELECT
               END DO
            END DO
         END DO
      ENDIF
   END SUBROUTINE dom_vvl_orca_fix

   !!======================================================================
END MODULE domvvl



