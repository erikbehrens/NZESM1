MODULE sbcice_lim
   !!======================================================================
   !!                       ***  MODULE  sbcice_lim  ***
   !! Surface module :  update the ocean surface boundary condition over ice
   !!       &           covered area using LIM sea-ice model
   !! Sea-Ice model  :  LIM-3 Sea ice model time-stepping
   !!=====================================================================
   !! History :  2.0  ! 2006-12  (M. Vancoppenolle) Original code
   !!            3.0  ! 2008-02  (C. Talandier)  Surface module from icestp.F90
   !!             -   ! 2008-04  (G. Madec)  sltyle and lim_ctl routine
   !!            3.3  ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!            3.4  ! 2011-01  (A Porter)  dynamical allocation
   !!             -   ! 2012-10  (C. Rousset)  add lim_diahsb
   !!            3.6  ! 2014-07  (M. Vancoppenolle, G. Madec, O. Marti) revise coupled interface
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3' :                                  LIM 3.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   sbc_ice_lim  : sea-ice model time-stepping and update ocean sbc over ice-covered area
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ice             ! LIM-3: ice variables
   USE thd_ice         ! LIM-3: thermodynamical variables
   USE dom_ice         ! LIM-3: ice domain

   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice   fields
   USE sbcblk_core     ! Surface boundary condition: CORE bulk
   USE sbcblk_clio     ! Surface boundary condition: CLIO bulk
   USE sbccpl          ! Surface boundary condition: coupled interface
   USE albedo          ! ocean & ice albedo

   USE phycst          ! Define parameters for the routines
   USE eosbn2          ! equation of state
   USE limdyn          ! Ice dynamics
   USE limtrp          ! Ice transport
   USE limhdf          ! Ice horizontal diffusion
   USE limthd          ! Ice thermodynamics
   USE limitd_me       ! Mechanics on ice thickness distribution
   USE limsbc          ! sea surface boundary condition
   USE limdiahsb       ! Ice budget diagnostics
   USE limwri          ! Ice outputs
   USE limrst          ! Ice restarts
   USE limupdate1      ! update of global variables
   USE limupdate2      ! update of global variables
   USE limvar          ! Ice variables switch

   USE limmsh          ! LIM mesh
   USE limistate       ! LIM initial state
   USE limthd_sal      ! LIM ice thermodynamics: salinity

   USE c1d             ! 1D vertical configuration
   USE lbclnk          ! lateral boundary condition - MPP link
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lib_fortran     ! 
   USE limctl

#if defined key_bdy 
   USE bdyice_lim       ! unstructured open boundary data  (bdy_ice_lim routine)
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_ice_lim  ! routine called by sbcmod.F90
   PUBLIC sbc_lim_init ! routine called by sbcmod.F90
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , UCL NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   !!======================================================================

   SUBROUTINE sbc_ice_lim( kt, kblk )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ice_lim  ***
      !!                   
      !! ** Purpose :   update the ocean surface boundary condition via the 
      !!                Louvain la Neuve Sea Ice Model time stepping 
      !!
      !! ** Method  :   ice model time stepping
      !!              - call the ice dynamics routine 
      !!              - call the ice advection/diffusion routine 
      !!              - call the ice thermodynamics routine 
      !!              - call the routine that computes mass and 
      !!                heat fluxes at the ice/ocean interface
      !!              - save the outputs 
      !!              - save the outputs for restart when necessary
      !!
      !! ** Action  : - time evolution of the LIM sea-ice model
      !!              - update all sbc variables below sea-ice:
      !!                utau, vtau, taum, wndm, qns , qsr, emp , sfx 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt      ! ocean time step
      INTEGER, INTENT(in) ::   kblk    ! type of bulk (=3 CLIO, =4 CORE, =5 COUPLED)
      !!
      INTEGER  ::   jl                 ! dummy loop index
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zalb_os, zalb_cs  ! ice albedo under overcast/clear sky
      REAL(wp), POINTER, DIMENSION(:,:  )   ::   zutau_ice, zvtau_ice 
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('sbc_ice_lim')

      !-----------------------!
      ! --- Ice time step --- !
      !-----------------------!
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN

         ! mean surface ocean current at ice velocity point (C-grid dynamics :  U- & V-points as the ocean)
         u_oce(:,:) = ssu_m(:,:) * umask(:,:,1)
         v_oce(:,:) = ssv_m(:,:) * vmask(:,:,1)
         
         ! masked sea surface freezing temperature [Kelvin] (set to rt0 over land)
         CALL eos_fzp( sss_m(:,:) , t_bo(:,:) )
         t_bo(:,:) = ( t_bo(:,:) + rt0 ) * tmask(:,:,1) + rt0 * ( 1._wp - tmask(:,:,1) )
          
         ! Mask sea ice surface temperature (set to rt0 over land)
         DO jl = 1, jpl
            t_su(:,:,jl) = t_su(:,:,jl) * tmask(:,:,1) + rt0 * ( 1._wp - tmask(:,:,1) )
         END DO     
         !
         !------------------------------------------------!                                           
         ! --- Dynamical coupling with the atmosphere --- !                                           
         !------------------------------------------------!
         ! It provides the following fields:
         ! utau_ice, vtau_ice : surface ice stress (U- & V-points)   [N/m2]
         !-----------------------------------------------------------------
         SELECT CASE( kblk )
         CASE( jp_clio    )   ;   CALL blk_ice_clio_tau                         ! CLIO bulk formulation            
         CASE( jp_core    )   ;   CALL blk_ice_core_tau                         ! CORE bulk formulation
         CASE( jp_purecpl )   ;   CALL sbc_cpl_ice_tau( utau_ice , vtau_ice )   ! Coupled   formulation
         END SELECT
         
         IF( ln_mixcpl) THEN   ! Case of a mixed Bulk/Coupled formulation
            CALL wrk_alloc( jpi,jpj    , zutau_ice, zvtau_ice)
            CALL sbc_cpl_ice_tau( zutau_ice , zvtau_ice )
            utau_ice(:,:) = utau_ice(:,:) * xcplmask(:,:,0) + zutau_ice(:,:) * ( 1. - xcplmask(:,:,0) )
            vtau_ice(:,:) = vtau_ice(:,:) * xcplmask(:,:,0) + zvtau_ice(:,:) * ( 1. - xcplmask(:,:,0) )
            CALL wrk_dealloc( jpi,jpj  , zutau_ice, zvtau_ice)
         ENDIF

         !-------------------------------------------------------!
         ! --- ice dynamics and transport (except in 1D case) ---!
         !-------------------------------------------------------!
         numit = numit + nn_fsbc                  ! Ice model time step
         !                                                   
         CALL sbc_lim_bef                         ! Store previous ice values
         CALL sbc_lim_diag0                       ! set diag of mass, heat and salt fluxes to 0
         CALL lim_rst_opn( kt )                   ! Open Ice restart file
         !
         IF( .NOT. lk_c1d ) THEN
            !
            CALL lim_dyn( kt )                    ! Ice dynamics    ( rheology/dynamics )   
            !
            CALL lim_trp( kt )                    ! Ice transport   ( Advection/diffusion )
            !
            IF( nn_monocat /= 2 ) CALL lim_itd_me ! Mechanical redistribution ! (ridging/rafting)
            !
#if defined key_bdy
            CALL bdy_ice_lim( kt )                ! bdy ice thermo 
            IF( ln_icectl )       CALL lim_prt( kt, iiceprt, jiceprt, 1, ' - ice thermo bdy - ' )
#endif
            !
            CALL lim_update1( kt )                ! Corrections
            !
         ENDIF
         
         ! previous lead fraction and ice volume for flux calculations
         CALL sbc_lim_bef                        
         CALL lim_var_glo2eqv                     ! ht_i and ht_s for ice albedo calculation
         CALL lim_var_agg(1)                      ! at_i for coupling (via pfrld) 
         pfrld(:,:)   = 1._wp - at_i(:,:)
         phicif(:,:)  = vt_i(:,:)
         
         !------------------------------------------------------!                                           
         ! --- Thermodynamical coupling with the atmosphere --- !                                           
         !------------------------------------------------------!
         ! It provides the following fields:
         ! qsr_ice , qns_ice  : solar & non solar heat flux over ice   (T-point)         [W/m2]
         ! qla_ice            : latent heat flux over ice              (T-point)         [W/m2]
         ! dqns_ice, dqla_ice : non solar & latent heat sensistivity   (T-point)         [W/m2]
         ! tprecip , sprecip  : total & solid precipitation            (T-point)         [Kg/m2/s]
         ! fr1_i0  , fr2_i0   : 1sr & 2nd fraction of qsr penetration in ice             [%]
         !----------------------------------------------------------------------------------------
         CALL wrk_alloc( jpi,jpj,jpl, zalb_os, zalb_cs )
         CALL albedo_ice( t_su, ht_i, ht_s, zalb_cs, zalb_os ) ! cloud-sky and overcast-sky ice albedos

         SELECT CASE( kblk )
         CASE( jp_clio )                                       ! CLIO bulk formulation
            ! In CLIO the cloud fraction is read in the climatology and the all-sky albedo 
            ! (alb_ice) is computed within the bulk routine
                                 CALL blk_ice_clio_flx( t_su, zalb_cs, zalb_os, alb_ice )
            IF( ln_mixcpl      ) CALL sbc_cpl_ice_flx( p_frld=pfrld, palbi=alb_ice, psst=sst_m, pist=t_su )
            IF( nn_limflx /= 2 ) CALL ice_lim_flx( t_su, alb_ice, qns_ice, qsr_ice, dqns_ice, evap_ice, devap_ice, nn_limflx )
         CASE( jp_core )                                       ! CORE bulk formulation
            ! albedo depends on cloud fraction because of non-linear spectral effects
            alb_ice(:,:,:) = ( 1. - cldf_ice ) * zalb_cs(:,:,:) + cldf_ice * zalb_os(:,:,:)
                                 CALL blk_ice_core_flx( t_su, alb_ice )
            IF( ln_mixcpl      ) CALL sbc_cpl_ice_flx( p_frld=pfrld, palbi=alb_ice, psst=sst_m, pist=t_su )
            IF( nn_limflx /= 2 ) CALL ice_lim_flx( t_su, alb_ice, qns_ice, qsr_ice, dqns_ice, evap_ice, devap_ice, nn_limflx )
         CASE ( jp_purecpl )
            ! albedo depends on cloud fraction because of non-linear spectral effects
            alb_ice(:,:,:) = ( 1. - cldf_ice ) * zalb_cs(:,:,:) + cldf_ice * zalb_os(:,:,:)
                                 CALL sbc_cpl_ice_flx( p_frld=pfrld, palbi=alb_ice, psst=sst_m, pist=t_su )
            IF( nn_limflx == 2 ) CALL ice_lim_flx( t_su, alb_ice, qns_ice, qsr_ice, dqns_ice, evap_ice, devap_ice, nn_limflx )
         END SELECT
         CALL wrk_dealloc( jpi,jpj,jpl, zalb_os, zalb_cs )

         !----------------------------!
         ! --- ice thermodynamics --- !
         !----------------------------!
         CALL lim_thd( kt )                         ! Ice thermodynamics      
         !
         CALL lim_update2( kt )                     ! Corrections
         !
         CALL lim_sbc_flx( kt )                     ! Update surface ocean mass, heat and salt fluxes
         !
         IF(ln_limdiaout) CALL lim_diahsb           ! Diagnostics and outputs 
         !
         CALL lim_wri( 1 )                          ! Ice outputs 
         !
         IF( kt == nit000 .AND. ln_rstart )   &
            &             CALL iom_close( numrir )  ! close input ice restart file
         !
         IF( lrst_ice )   CALL lim_rst_write( kt )  ! Ice restart file 
         !
         IF( ln_icectl )  CALL lim_ctl( kt )        ! alerts in case of model crash
         !
      ENDIF   ! End sea-ice time step only

      !-------------------------!
      ! --- Ocean time step --- !
      !-------------------------!
      ! Update surface ocean stresses (only in ice-dynamic case) otherwise the atm.-ocean stresses are used everywhere
      IF( ln_limdyn )     CALL lim_sbc_tau( kt, ub(:,:,1), vb(:,:,1) )  ! using before instantaneous surf. currents
!!gm   remark, the ocean-ice stress is not saved in ice diag call above .....  find a solution!!!
      !
      IF( nn_timing == 1 ) CALL timing_stop('sbc_ice_lim')
      !
   END SUBROUTINE sbc_ice_lim
   

   SUBROUTINE sbc_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_lim_init  ***
      !!
      !! ** purpose :   Allocate all the dynamic arrays of the LIM-3 modules
      !!----------------------------------------------------------------------
      INTEGER :: ierr
      INTEGER :: ji, jj
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'sbc_ice_lim : update ocean surface boudary condition' 
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   via Louvain la Neuve Ice Model (LIM-3) time stepping'
      !
                                       ! Open the reference and configuration namelist files and namelist output file 
      CALL ctl_opn( numnam_ice_ref, 'namelist_ice_ref',    'OLD',     'FORMATTED', 'SEQUENTIAL', -1, numout, lwp ) 
      CALL ctl_opn( numnam_ice_cfg, 'namelist_ice_cfg',    'OLD',     'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      IF(lwm) CALL ctl_opn( numoni, 'output.namelist.ice', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, 1 )

      CALL ice_run                     ! set some ice run parameters
      !
      !                                ! Allocate the ice arrays
      ierr =        ice_alloc        ()      ! ice variables
      ierr = ierr + dom_ice_alloc    ()      ! domain
      ierr = ierr + sbc_ice_alloc    ()      ! surface forcing
      ierr = ierr + thd_ice_alloc    ()      ! thermodynamics
      ierr = ierr + lim_itd_me_alloc ()      ! ice thickness distribution - mechanics
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'sbc_lim_init : unable to allocate ice arrays')
      !
      !                                ! adequation jpk versus ice/snow layers/categories
      IF( jpl > jpk .OR. (nlay_i+1) > jpk .OR. nlay_s > jpk )   &
         &      CALL ctl_stop( 'STOP',                          &
         &     'sbc_lim_init: the 3rd dimension of workspace arrays is too small.',   &
         &     'use more ocean levels or less ice/snow layers/categories.' )
      !
      CALL lim_itd_init                ! ice thickness distribution initialization
      !
      CALL lim_hdf_init                ! set ice horizontal diffusion computation parameters
      !
      CALL lim_thd_init                ! set ice thermodynics parameters
      !
      CALL lim_thd_sal_init            ! set ice salinity parameters
      !
      CALL lim_msh                     ! ice mesh initialization
      !
      CALL lim_itd_me_init             ! ice thickness distribution initialization for mecanical deformation
      !                                ! Initial sea-ice state
      IF( .NOT. ln_rstart ) THEN              ! start from rest: sea-ice deduced from sst
         numit = 0
         numit = nit000 - 1
         CALL lim_istate
      ELSE                                    ! start from a restart file
         CALL lim_rst_read
         numit = nit000 - 1
      ENDIF
      CALL lim_var_agg(1)
      CALL lim_var_glo2eqv
      !
      CALL lim_sbc_init                 ! ice surface boundary condition   
      !
      fr_i(:,:)     = at_i(:,:)         ! initialisation of sea-ice fraction
      tn_ice(:,:,:) = t_su(:,:,:)       ! initialisation of surface temp for coupled simu
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( gphit(ji,jj) > 0._wp ) THEN  ;  rn_amax_2d(ji,jj) = rn_amax_n  ! NH
            ELSE                             ;  rn_amax_2d(ji,jj) = rn_amax_s  ! SH
            ENDIF
        ENDDO
      ENDDO 
      !
      nstart = numit  + nn_fsbc      
      nitrun = nitend - nit000 + 1 
      nlast  = numit  + nitrun 
      !
      IF( nstock == 0 )   nstock = nlast + 1
      !
   END SUBROUTINE sbc_lim_init


   SUBROUTINE ice_run
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_run ***
      !!                 
      !! ** Purpose :   Definition some run parameter for ice model
      !!
      !! ** Method  :   Read the namicerun namelist and check the parameter 
      !!              values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicerun
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      NAMELIST/namicerun/ jpl, nlay_i, nlay_s, cn_icerst_in, cn_icerst_indir, cn_icerst_out, cn_icerst_outdir,  &
         &                ln_limdyn, rn_amax_n, rn_amax_s, ln_limdiahsb, ln_limdiaout, ln_icectl, iiceprt, jiceprt  
      !!-------------------------------------------------------------------
      !                    
      REWIND( numnam_ice_ref )              ! Namelist namicerun in reference namelist : Parameters for ice
      READ  ( numnam_ice_ref, namicerun, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicerun in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namicerun in configuration namelist : Parameters for ice
      READ  ( numnam_ice_cfg, namicerun, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicerun in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namicerun )
      !
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_run : ice share parameters for dynamics/advection/thermo of sea-ice'
         WRITE(numout,*) ' ~~~~~~'
         WRITE(numout,*) '   number of ice  categories                               = ', jpl
         WRITE(numout,*) '   number of ice  layers                                   = ', nlay_i
         WRITE(numout,*) '   number of snow layers                                   = ', nlay_s
         WRITE(numout,*) '   switch for ice dynamics (1) or not (0)      ln_limdyn   = ', ln_limdyn
         WRITE(numout,*) '   maximum ice concentration for NH                        = ', rn_amax_n 
         WRITE(numout,*) '   maximum ice concentration for SH                        = ', rn_amax_s
         WRITE(numout,*) '   Diagnose heat/salt budget or not          ln_limdiahsb  = ', ln_limdiahsb
         WRITE(numout,*) '   Output   heat/salt budget or not          ln_limdiaout  = ', ln_limdiaout
         WRITE(numout,*) '   control prints in ocean.out for (i,j)=(iiceprt,jiceprt) = ', ln_icectl
         WRITE(numout,*) '   i-index for control prints (ln_icectl=true)             = ', iiceprt
         WRITE(numout,*) '   j-index for control prints (ln_icectl=true)             = ', jiceprt
      ENDIF
      !
      ! sea-ice timestep and inverse
      rdt_ice   = nn_fsbc * rdttra(1)  
      r1_rdtice = 1._wp / rdt_ice 

      ! inverse of nlay_i and nlay_s
      r1_nlay_i = 1._wp / REAL( nlay_i, wp )
      r1_nlay_s = 1._wp / REAL( nlay_s, wp )
      !
#if defined key_bdy
      IF( lwp .AND. ln_limdiahsb )  CALL ctl_warn('online conservation check activated but it does not work with BDY')
#endif
      !
   END SUBROUTINE ice_run


   SUBROUTINE lim_itd_init
      !!------------------------------------------------------------------
      !!                ***  ROUTINE lim_itd_init ***
      !!
      !! ** Purpose :   Initializes the ice thickness distribution
      !! ** Method  :   ...
      !! ** input   :   Namelist namiceitd
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      NAMELIST/namiceitd/ nn_catbnd, rn_himean
      !
      INTEGER  ::   jl                   ! dummy loop index
      REAL(wp) ::   zc1, zc2, zc3, zx1   ! local scalars
      REAL(wp) ::   zhmax, znum, zden, zalpha !
      !!------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )              ! Namelist namiceitd in reference namelist : Parameters for ice
      READ  ( numnam_ice_ref, namiceitd, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceitd in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namiceitd in configuration namelist : Parameters for ice
      READ  ( numnam_ice_cfg, namiceitd, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceitd in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namiceitd )
      !
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_itd : ice cat distribution'
         WRITE(numout,*) ' ~~~~~~'
         WRITE(numout,*) '   shape of ice categories distribution                          nn_catbnd = ', nn_catbnd
         WRITE(numout,*) '   mean ice thickness in the domain (only active if nn_catbnd=2) rn_himean = ', rn_himean
      ENDIF

      !----------------------------------
      !- Thickness categories boundaries 
      !----------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'lim_itd_init : Initialization of ice cat distribution '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

      hi_max(:) = 0._wp

      SELECT CASE ( nn_catbnd  )       
                                   !----------------------
         CASE (1)                  ! tanh function (CICE)
                                   !----------------------
         zc1 =  3._wp / REAL( jpl, wp )
         zc2 = 10._wp * zc1
         zc3 =  3._wp

         DO jl = 1, jpl
            zx1 = REAL( jl-1, wp ) / REAL( jpl, wp )
            hi_max(jl) = hi_max(jl-1) + zc1 + zc2 * (1._wp + TANH( zc3 * (zx1 - 1._wp ) ) )
         END DO

                                   !----------------------
         CASE (2)                  ! h^(-alpha) function
                                   !----------------------
         zalpha = 0.05             ! exponent of the transform function

         zhmax  = 3.*rn_himean

         DO jl = 1, jpl 
            znum = jpl * ( zhmax+1 )**zalpha
            zden = ( jpl - jl ) * ( zhmax+1 )**zalpha + jl
            hi_max(jl) = ( znum / zden )**(1./zalpha) - 1
         END DO

      END SELECT

      DO jl = 1, jpl
         hi_mean(jl) = ( hi_max(jl) + hi_max(jl-1) ) * 0.5_wp
      END DO

      ! Set hi_max(jpl) to a big value to ensure that all ice is thinner than hi_max(jpl)
      hi_max(jpl) = 99._wp

      IF(lwp) WRITE(numout,*) ' Thickness category boundaries '
      IF(lwp) WRITE(numout,*) ' hi_max ', hi_max(0:jpl)
      !
   END SUBROUTINE lim_itd_init

   
   SUBROUTINE ice_lim_flx( ptn_ice, palb_ice, pqns_ice, pqsr_ice, pdqn_ice, pevap_ice, pdevap_ice, k_limflx )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ice_lim_flx  ***
      !!                   
      !! ** Purpose :   update the ice surface boundary condition by averaging and / or
      !!                redistributing fluxes on ice categories                   
      !!
      !! ** Method  :   average then redistribute 
      !!
      !! ** Action  :   
      !!---------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   k_limflx   ! =-1 do nothing; =0 average ; 
                                                                ! =1 average and redistribute ; =2 redistribute
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptn_ice    ! ice surface temperature 
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   palb_ice   ! ice albedo
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pqns_ice   ! non solar flux
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pqsr_ice   ! net solar flux
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdqn_ice   ! non solar flux sensitivity
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pevap_ice  ! sublimation
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdevap_ice ! sublimation sensitivity
      !
      INTEGER  ::   jl      ! dummy loop index
      !
      REAL(wp), POINTER, DIMENSION(:,:) :: zalb_m    ! Mean albedo over all categories
      REAL(wp), POINTER, DIMENSION(:,:) :: ztem_m    ! Mean temperature over all categories
      !
      REAL(wp), POINTER, DIMENSION(:,:) :: z_qsr_m   ! Mean solar heat flux over all categories
      REAL(wp), POINTER, DIMENSION(:,:) :: z_qns_m   ! Mean non solar heat flux over all categories
      REAL(wp), POINTER, DIMENSION(:,:) :: z_evap_m  ! Mean sublimation over all categories
      REAL(wp), POINTER, DIMENSION(:,:) :: z_dqn_m   ! Mean d(qns)/dT over all categories
      REAL(wp), POINTER, DIMENSION(:,:) :: z_devap_m ! Mean d(evap)/dT over all categories
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('ice_lim_flx')
      !
      !
      SELECT CASE( k_limflx )                              !==  averaged on all ice categories  ==!
      CASE( 0 , 1 )
         CALL wrk_alloc( jpi,jpj, z_qsr_m, z_qns_m, z_evap_m, z_dqn_m, z_devap_m)
         !
         z_qns_m  (:,:) = fice_ice_ave ( pqns_ice (:,:,:) )
         z_qsr_m  (:,:) = fice_ice_ave ( pqsr_ice (:,:,:) )
         z_dqn_m  (:,:) = fice_ice_ave ( pdqn_ice (:,:,:) )
         z_evap_m (:,:) = fice_ice_ave ( pevap_ice (:,:,:) )
         z_devap_m(:,:) = fice_ice_ave ( pdevap_ice (:,:,:) )
         DO jl = 1, jpl
            pdqn_ice  (:,:,jl) = z_dqn_m(:,:)
            pdevap_ice(:,:,jl) = z_devap_m(:,:)
         END DO
         !
         DO jl = 1, jpl
            pqns_ice (:,:,jl) = z_qns_m(:,:)
            pqsr_ice (:,:,jl) = z_qsr_m(:,:)
            pevap_ice(:,:,jl) = z_evap_m(:,:)
         END DO
         !
         CALL wrk_dealloc( jpi,jpj, z_qsr_m, z_qns_m, z_evap_m, z_dqn_m, z_devap_m)
      END SELECT

      SELECT CASE( k_limflx )                              !==  redistribution on all ice categories  ==!
      CASE( 1 , 2 )
         CALL wrk_alloc( jpi,jpj, zalb_m, ztem_m )
         !
         zalb_m(:,:) = fice_ice_ave ( palb_ice (:,:,:) ) 
         ztem_m(:,:) = fice_ice_ave ( ptn_ice  (:,:,:) ) 
         DO jl = 1, jpl
            pqns_ice (:,:,jl) = pqns_ice (:,:,jl) + pdqn_ice  (:,:,jl) * ( ptn_ice(:,:,jl) - ztem_m(:,:) )
            pevap_ice(:,:,jl) = pevap_ice(:,:,jl) + pdevap_ice(:,:,jl) * ( ptn_ice(:,:,jl) - ztem_m(:,:) )
            pqsr_ice (:,:,jl) = pqsr_ice (:,:,jl) * ( 1._wp - palb_ice(:,:,jl) ) / ( 1._wp - zalb_m(:,:) ) 
         END DO
         !
         CALL wrk_dealloc( jpi,jpj, zalb_m, ztem_m )
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('ice_lim_flx')
      !
   END SUBROUTINE ice_lim_flx

   SUBROUTINE sbc_lim_bef
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_lim_bef  ***
      !!
      !! ** purpose :  store ice variables at "before" time step 
      !!----------------------------------------------------------------------
      a_i_b  (:,:,:)   = a_i  (:,:,:)     ! ice area
      e_i_b  (:,:,:,:) = e_i  (:,:,:,:)   ! ice thermal energy
      v_i_b  (:,:,:)   = v_i  (:,:,:)     ! ice volume
      v_s_b  (:,:,:)   = v_s  (:,:,:)     ! snow volume 
      e_s_b  (:,:,:,:) = e_s  (:,:,:,:)   ! snow thermal energy
      smv_i_b(:,:,:)   = smv_i(:,:,:)     ! salt content
      oa_i_b (:,:,:)   = oa_i (:,:,:)     ! areal age content
      u_ice_b(:,:)     = u_ice(:,:)
      v_ice_b(:,:)     = v_ice(:,:)
      
   END SUBROUTINE sbc_lim_bef

   SUBROUTINE sbc_lim_diag0
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_lim_diag0  ***
      !!
      !! ** purpose :  set ice-ocean and ice-atm. fluxes to zeros at the beggining
      !!               of the time step
      !!----------------------------------------------------------------------
      sfx    (:,:) = 0._wp   ;
      sfx_bri(:,:) = 0._wp   ; 
      sfx_sni(:,:) = 0._wp   ;   sfx_opw(:,:) = 0._wp
      sfx_bog(:,:) = 0._wp   ;   sfx_dyn(:,:) = 0._wp
      sfx_bom(:,:) = 0._wp   ;   sfx_sum(:,:) = 0._wp
      sfx_res(:,:) = 0._wp   ;   sfx_sub(:,:) = 0._wp
      
      wfx_snw(:,:) = 0._wp   ;   wfx_ice(:,:) = 0._wp
      wfx_sni(:,:) = 0._wp   ;   wfx_opw(:,:) = 0._wp
      wfx_bog(:,:) = 0._wp   ;   wfx_dyn(:,:) = 0._wp
      wfx_bom(:,:) = 0._wp   ;   wfx_sum(:,:) = 0._wp
      wfx_res(:,:) = 0._wp   ;   wfx_sub(:,:) = 0._wp
      wfx_spr(:,:) = 0._wp   ;   
      
      hfx_thd(:,:) = 0._wp   ;   
      hfx_snw(:,:) = 0._wp   ;   hfx_opw(:,:) = 0._wp
      hfx_bog(:,:) = 0._wp   ;   hfx_dyn(:,:) = 0._wp
      hfx_bom(:,:) = 0._wp   ;   hfx_sum(:,:) = 0._wp
      hfx_res(:,:) = 0._wp   ;   hfx_sub(:,:) = 0._wp
      hfx_spr(:,:) = 0._wp   ;   hfx_dif(:,:) = 0._wp 
      hfx_err(:,:) = 0._wp   ;   hfx_err_rem(:,:) = 0._wp
      hfx_err_dif(:,:) = 0._wp
      wfx_err_sub(:,:) = 0._wp
      
      afx_tot(:,:) = 0._wp   ;
      afx_dyn(:,:) = 0._wp   ;   afx_thd(:,:) = 0._wp

      diag_heat(:,:) = 0._wp ;   diag_smvi(:,:) = 0._wp ;
      diag_vice(:,:) = 0._wp ;   diag_vsnw(:,:) = 0._wp ;
      
   END SUBROUTINE sbc_lim_diag0

     
   FUNCTION fice_cell_ave ( ptab )
      !!--------------------------------------------------------------------------
      !! * Compute average over categories, for grid cell (ice covered and free ocean)
      !!--------------------------------------------------------------------------
      REAL (wp), DIMENSION (jpi,jpj) :: fice_cell_ave
      REAL (wp), DIMENSION (jpi,jpj,jpl), INTENT (in) :: ptab
      INTEGER :: jl ! Dummy loop index
      
      fice_cell_ave (:,:) = 0.0_wp
      DO jl = 1, jpl
         fice_cell_ave (:,:) = fice_cell_ave (:,:) + a_i (:,:,jl) * ptab (:,:,jl)
      END DO
      
   END FUNCTION fice_cell_ave
   
   
   FUNCTION fice_ice_ave ( ptab )
      !!--------------------------------------------------------------------------
      !! * Compute average over categories, for ice covered part of grid cell
      !!--------------------------------------------------------------------------
      REAL (kind=wp), DIMENSION (jpi,jpj) :: fice_ice_ave
      REAL (kind=wp), DIMENSION (jpi,jpj,jpl), INTENT(in) :: ptab

      fice_ice_ave (:,:) = 0.0_wp
      WHERE ( at_i (:,:) > 0.0_wp ) fice_ice_ave (:,:) = fice_cell_ave ( ptab (:,:,:)) / at_i (:,:)

   END FUNCTION fice_ice_ave


#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module      NO LIM 3.0 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE sbc_ice_lim ( kt, kblk )     ! Dummy routine
      WRITE(*,*) 'sbc_ice_lim: You should not have seen this print! error?', kt, kblk
   END SUBROUTINE sbc_ice_lim
   SUBROUTINE sbc_lim_init                 ! Dummy routine
   END SUBROUTINE sbc_lim_init
#endif

   !!======================================================================
END MODULE sbcice_lim
