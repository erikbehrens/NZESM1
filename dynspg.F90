MODULE dynspg
   !!======================================================================
   !!                       ***  MODULE  dynspg  ***
   !! Ocean dynamics:  surface pressure gradient control
   !!======================================================================
   !! History :  1.0  ! 2005-12  (C. Talandier, G. Madec, V. Garnier)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg     : update the dynamics trend with the lateral diffusion
   !!   dyn_spg_ctl : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE c1d            ! 1D vertical configuration
   USE phycst         ! physical constants
   USE sbc_oce        ! surface boundary condition: ocean
   USE sbcapr         ! surface boundary condition: atmospheric pressure
   USE dynspg_oce     ! surface pressure gradient variables
   USE dynspg_exp     ! surface pressure gradient     (dyn_spg_exp routine)
   USE dynspg_ts      ! surface pressure gradient     (dyn_spg_ts  routine)
   USE dynspg_flt     ! surface pressure gradient     (dyn_spg_flt routine)
   USE dynadv         ! dynamics: vector invariant versus flux form
   USE dynhpg, ONLY: ln_dynhpg_imp
   USE sbctide
   USE updtide
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE prtctl         ! Print control                     (prt_ctl routine)
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE solver         ! solver initialization
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing


   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_spg        ! routine called by step module
   PUBLIC   dyn_spg_init   ! routine called by opa module

   INTEGER ::   nspg = 0   ! type of surface pressure gradient scheme defined from lk_dynspg_... 

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_spg( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg  ***
      !!
      !! ** Purpose :   achieve the momentum time stepping by computing the
      !!              last trend, the surface pressure gradient including the 
      !!              atmospheric pressure forcing (ln_apr_dyn=T), and performing
      !!              the Leap-Frog integration.
      !!gm              In the current version only the filtered solution provide
      !!gm            the after velocity, in the 2 other (ua,va) are still the trends
      !!
      !! ** Method  :   Three schemes:
      !!              - explicit computation      : the spg is evaluated at now
      !!              - filtered computation      : the Roulet & madec (2000) technique is used
      !!              - split-explicit computation: a time splitting technique is used
      !!
      !!              ln_apr_dyn=T : the atmospheric pressure forcing is applied 
      !!             as the gradient of the inverse barometer ssh:
      !!                apgu = - 1/rau0 di[apr] = 0.5*grav di[ssh_ib+ssh_ibb]
      !!                apgv = - 1/rau0 dj[apr] = 0.5*grav dj[ssh_ib+ssh_ibb]
      !!             Note that as all external forcing a time averaging over a two rdt
      !!             period is used to prevent the divergence of odd and even time step.
      !!
      !! N.B. : When key_esopa is used all the scheme are tested, regardless 
      !!        of the physical meaning of the results. 
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER, INTENT(  out) ::   kindic   ! solver flag
      !
      INTEGER  ::   ji, jj, jk                             ! dummy loop indices
      REAL(wp) ::   z2dt, zg_2, zintp, zgrau0r             ! temporary scalar
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      REAL(wp), POINTER, DIMENSION(:,:)   ::  zpice
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg')
      !

!!gm NOTA BENE : the dynspg_exp and dynspg_ts should be modified so that 
!!gm             they return the after velocity, not the trends (as in trazdf_imp...)
!!gm             In this case, change/simplify dynnxt


      IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF

      IF(      ln_apr_dyn                                                &   ! atmos. pressure
         .OR.  ( .NOT.lk_dynspg_ts .AND. (ln_tide_pot .AND. lk_tide) )   &   ! tide potential (no time slitting)
         .OR.  nn_ice_embd == 2  ) THEN                                      ! embedded sea-ice
         !
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               spgu(ji,jj) = 0._wp
               spgv(ji,jj) = 0._wp
            END DO
         END DO         
         !
         IF( ln_apr_dyn .AND. (.NOT. lk_dynspg_ts) ) THEN                    !==  Atmospheric pressure gradient (added later in time-split case) ==!
            zg_2 = grav * 0.5
            DO jj = 2, jpjm1                          ! gradient of Patm using inverse barometer ssh
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  spgu(ji,jj) = spgu(ji,jj) + zg_2 * (  ssh_ib (ji+1,jj) - ssh_ib (ji,jj)    &
                     &                      + ssh_ibb(ji+1,jj) - ssh_ibb(ji,jj)  ) /e1u(ji,jj)
                  spgv(ji,jj) = spgv(ji,jj) + zg_2 * (  ssh_ib (ji,jj+1) - ssh_ib (ji,jj)    &
                     &                      + ssh_ibb(ji,jj+1) - ssh_ibb(ji,jj)  ) /e2v(ji,jj)
               END DO
            END DO
         ENDIF
         !
         !                                    !==  tide potential forcing term  ==!
         IF( .NOT.lk_dynspg_ts .AND. ( ln_tide_pot .AND. lk_tide )  ) THEN   ! N.B. added directly at sub-time-step in ts-case
            !
            CALL upd_tide( kt )                      ! update tide potential
            !
            DO jj = 2, jpjm1                         ! add tide potential forcing
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  spgu(ji,jj) = spgu(ji,jj) + grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) / e1u(ji,jj)
                  spgv(ji,jj) = spgv(ji,jj) + grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) / e2v(ji,jj)
               END DO 
            END DO
         ENDIF
         !
         IF( nn_ice_embd == 2 ) THEN          !== embedded sea ice: Pressure gradient due to snow-ice mass ==!
            CALL wrk_alloc( jpi, jpj, zpice )
            !                                            
            zintp = REAL( MOD( kt-1, nn_fsbc ) ) / REAL( nn_fsbc )
            zgrau0r     = - grav * r1_rau0
            zpice(:,:) = (  zintp * snwice_mass(:,:) + ( 1.- zintp ) * snwice_mass_b(:,:)  ) * zgrau0r
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  spgu(ji,jj) = spgu(ji,jj) + ( zpice(ji+1,jj) - zpice(ji,jj) ) / e1u(ji,jj)
                  spgv(ji,jj) = spgv(ji,jj) + ( zpice(ji,jj+1) - zpice(ji,jj) ) / e2v(ji,jj)
               END DO
            END DO
            !
            CALL wrk_dealloc( jpi, jpj, zpice )         
         ENDIF
         !
         DO jk = 1, jpkm1                     !== Add all terms to the general trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
                  va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
               END DO
            END DO
         END DO    
         
!!gm add here a call to dyn_trd for ice pressure gradient, the surf pressure trends ????
              
      ENDIF

      SELECT CASE ( nspg )                       ! compute surf. pressure gradient trend and add it to the general trend
      !                                                     
      CASE (  0 )   ;   CALL dyn_spg_exp( kt )              ! explicit
      CASE (  1 )   ;   CALL dyn_spg_ts ( kt )              ! time-splitting
      CASE (  2 )   ;   CALL dyn_spg_flt( kt, kindic )      ! filtered
      !                                                    
      CASE ( -1 )                                ! esopa: test all possibility with control print
                        CALL dyn_spg_exp( kt )
                        CALL prt_ctl( tab3d_1=ua, clinfo1=' spg0 - Ua: ', mask1=umask, &
         &                            tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
                        CALL dyn_spg_ts ( kt )
                        CALL prt_ctl( tab3d_1=ua, clinfo1=' spg1 - Ua: ', mask1=umask, &
         &                           tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
                        CALL dyn_spg_flt( kt, kindic )
                        CALL prt_ctl( tab3d_1=ua, clinfo1=' spg2 - Ua: ', mask1=umask, &
         &                            tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      END SELECT
      !                    
      IF( l_trddyn )   THEN                      ! save the surface pressure gradient trends for further diagnostics
         SELECT CASE ( nspg )
         CASE ( 0, 1 )
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CASE( 2 )
            z2dt = 2. * rdt
            IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt
            ztrdu(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) / z2dt - ztrdu(:,:,:)
            ztrdv(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) / z2dt - ztrdv(:,:,:)
         END SELECT
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_spg, kt )
         !
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg  - Ua: ', mask1=umask, &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg')
      !
   END SUBROUTINE dyn_spg


   SUBROUTINE dyn_spg_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg_init  ***
      !!                
      !! ** Purpose :   Control the consistency between cpp options for 
      !!              surface pressure gradient schemes
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_init')
      !
      IF(lwp) THEN             ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_spg_init : choice of the surface pressure gradient scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '     Explicit free surface                  lk_dynspg_exp = ', lk_dynspg_exp
         WRITE(numout,*) '     Free surface with time splitting       lk_dynspg_ts  = ', lk_dynspg_ts
         WRITE(numout,*) '     Filtered free surface cst volume       lk_dynspg_flt = ', lk_dynspg_flt
      ENDIF

      IF( lk_dynspg_ts ) CALL dyn_spg_ts_init( nit000 )
      ! (do it now, to set nn_baro, used to allocate some arrays later on)
      !                        ! allocate dyn_spg arrays
      IF( lk_dynspg_ts ) THEN
         IF( dynspg_oce_alloc() /= 0 )   CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_oce arrays')
         IF( dyn_spg_ts_alloc() /= 0 )   CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_ts  arrays')
         IF ((neuler/=0).AND.(ln_bt_fw)) CALL ts_rst( nit000, 'READ' )
      ENDIF

      !                        ! Control of surface pressure gradient scheme options
      ioptio = 0
      IF(lk_dynspg_exp)   ioptio = ioptio + 1
      IF(lk_dynspg_ts )   ioptio = ioptio + 1
      IF(lk_dynspg_flt)   ioptio = ioptio + 1
      !
      IF( ( ioptio > 1 .AND. .NOT. lk_esopa ) .OR. ( ioptio == 0 .AND. .NOT. lk_c1d ) )   &
           &   CALL ctl_stop( ' Choose only one surface pressure gradient scheme with a key cpp' )
      IF( ( lk_dynspg_ts .OR. lk_dynspg_exp ) .AND. ln_isfcav )   &
           &   CALL ctl_stop( ' dynspg_ts and dynspg_exp not tested with ice shelf cavity ' )
      !
      IF( lk_esopa     )   nspg = -1
      IF( lk_dynspg_exp)   nspg =  0
      IF( lk_dynspg_ts )   nspg =  1
      IF( lk_dynspg_flt)   nspg =  2
      !
      IF( lk_esopa     )   nspg = -1
      !
      IF(lwp) THEN
         WRITE(numout,*)
         IF( nspg == -1 )   WRITE(numout,*) '     ESOPA test All scheme used'
         IF( nspg ==  0 )   WRITE(numout,*) '     explicit free surface'
         IF( nspg ==  1 )   WRITE(numout,*) '     free surface with time splitting scheme'
         IF( nspg ==  2 )   WRITE(numout,*) '     filtered free surface'
      ENDIF

#if defined key_dynspg_flt || defined key_esopa
      CALL solver_init( nit000 )   ! Elliptic solver initialisation
#endif

      !                        ! Control of timestep choice
      IF( lk_dynspg_ts .OR. lk_dynspg_exp ) THEN
         IF( nn_cla == 1 )   CALL ctl_stop( 'Crossland advection not implemented for this free surface formulation' )
      ENDIF

      !			       ! Control of hydrostatic pressure choice
      IF( lk_dynspg_ts .AND. ln_dynhpg_imp ) THEN
         CALL ctl_stop( 'Semi-implicit hpg not compatible with time splitting' )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_init')
      !
   END SUBROUTINE dyn_spg_init

  !!======================================================================
END MODULE dynspg
