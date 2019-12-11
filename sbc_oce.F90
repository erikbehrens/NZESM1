MODULE sbc_oce
   !!======================================================================
   !!                       ***  MODULE  sbc_oce  ***
   !! Surface module :   variables defined in core memory 
   !!======================================================================
   !! History :  3.0  ! 2006-06  (G. Madec)  Original code
   !!             -   ! 2008-08  (G. Madec)  namsbc moved from sbcmod
   !!            3.3  ! 2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!             -   ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!            3.3  ! 2010-10  (J. Chanut, C. Bricaud)  add the surface pressure forcing
   !!            4.0  ! 2012-05  (C. Rousset) add attenuation coef for use in ice model 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_oce_alloc : allocation of sbc arrays
   !!   sbc_tau2wnd   : wind speed estimated from wind stress
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_oce_alloc   ! routine called in sbcmod.F90
   PUBLIC   sbc_tau2wnd     ! routine called in several sbc modules
   
   !!----------------------------------------------------------------------
   !!           Namelist for the Ocean Surface Boundary Condition
   !!----------------------------------------------------------------------
   !                                   !!* namsbc namelist *
   LOGICAL , PUBLIC ::   ln_ana         !: analytical boundary condition flag
   LOGICAL , PUBLIC ::   ln_flx         !: flux      formulation
   LOGICAL , PUBLIC ::   ln_blk_clio    !: CLIO bulk formulation
   LOGICAL , PUBLIC ::   ln_blk_core    !: CORE bulk formulation
   LOGICAL , PUBLIC ::   ln_blk_mfs     !: MFS  bulk formulation
#if defined key_oasis3
   LOGICAL , PUBLIC ::   lk_oasis = .TRUE.  !: OASIS used
#else
   LOGICAL , PUBLIC ::   lk_oasis = .FALSE. !: OASIS unused
#endif
   LOGICAL , PUBLIC ::   ln_cpl         !: ocean-atmosphere coupled formulation
   LOGICAL , PUBLIC ::   ln_mixcpl      !: ocean-atmosphere forced-coupled mixed formulation
   LOGICAL , PUBLIC ::   ln_dm2dc       !: Daily mean to Diurnal Cycle short wave (qsr)
   LOGICAL , PUBLIC ::   ln_rnf         !: runoffs / runoff mouths
   LOGICAL , PUBLIC ::   ln_ssr         !: Sea Surface restoring on SST and/or SSS      
   LOGICAL , PUBLIC ::   ln_apr_dyn     !: Atmospheric pressure forcing used on dynamics (ocean & ice)
   INTEGER , PUBLIC ::   nn_ice         !: flag for ice in the surface boundary condition (=0/1/2/3)
   INTEGER , PUBLIC ::   nn_isf         !: flag for isf in the surface boundary condition (=0/1/2/3/4) 
   INTEGER , PUBLIC ::   nn_ice_embd    !: flag for levitating/embedding sea-ice in the ocean
   !                                             !: =0 levitating ice (no mass exchange, concentration/dilution effect)
   !                                             !: =1 levitating ice with mass and salt exchange but no presure effect
   !                                             !: =2 embedded sea-ice (full salt and mass exchanges and pressure)
   INTEGER , PUBLIC ::   nn_components  !: flag for sbc module (including sea-ice) coupling mode (see component definition below) 
   INTEGER , PUBLIC ::   nn_limflx      !: LIM3 Multi-category heat flux formulation
   !                                             !: =-1  Use of per-category fluxes
   !                                             !: = 0  Average per-category fluxes
   !                                             !: = 1  Average then redistribute per-category fluxes
   !                                             !: = 2  Redistribute a single flux over categories
   INTEGER , PUBLIC ::   nn_fwb         !: FreshWater Budget: 
   !                                             !:  = 0 unchecked 
   !                                             !:  = 1 global mean of e-p-r set to zero at each nn_fsbc time step
   !                                             !:  = 2 annual global mean of e-p-r set to zero
   LOGICAL , PUBLIC ::   ln_wave        !: true if some coupling with wave model
   LOGICAL , PUBLIC ::   ln_cdgw        !: true if neutral drag coefficient from wave model
   LOGICAL , PUBLIC ::   ln_sdw         !: true if 3d stokes drift from wave model
   !
   LOGICAL , PUBLIC ::   ln_icebergs    !: Icebergs
   !
   INTEGER , PUBLIC ::   nn_lsm         !: Number of iteration if seaoverland is applied
   !!----------------------------------------------------------------------
   !!           switch definition (improve readability)
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC, PARAMETER ::   jp_gyre    = 0        !: GYRE analytical               formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_ana     = 1        !: analytical                    formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_flx     = 2        !: flux                          formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_clio    = 3        !: CLIO bulk                     formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_core    = 4        !: CORE bulk                     formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_purecpl = 5        !: Pure ocean-atmosphere Coupled formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_mfs     = 6        !: MFS  bulk                     formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_none    = 7        !: for OPA when doing coupling via SAS module
   INTEGER , PUBLIC, PARAMETER ::   jp_esopa   = -1       !: esopa test, ALL formulations
   
   !!----------------------------------------------------------------------
   !!           component definition
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC, PARAMETER ::   jp_iam_nemo = 0      !: Initial single executable configuration 
                                                         !  (no internal OASIS coupling)
   INTEGER , PUBLIC, PARAMETER ::   jp_iam_opa  = 1      !: Multi executable configuration - OPA component
                                                         !  (internal OASIS coupling)
   INTEGER , PUBLIC, PARAMETER ::   jp_iam_sas  = 2      !: Multi executable configuration - SAS component
                                                         !  (internal OASIS coupling)
   !!----------------------------------------------------------------------
   !!              Ocean Surface Boundary Condition fields
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC ::  ncpl_qsr_freq            !: qsr coupling frequency per days from atmosphere
   !
   LOGICAL , PUBLIC ::   lhftau = .FALSE.        !: HF tau used in TKE: mean(stress module) - module(mean stress)
   !!                                   !!   now    ! before   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau   , utau_b   !: sea surface i-stress (ocean referential)     [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vtau   , vtau_b   !: sea surface j-stress (ocean referential)     [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   taum              !: module of sea surface stress (at T-point)    [N/m2] 
   !! wndm is used onmpute surface gases exchanges in ice-free ocean or leads
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wndm              !: wind speed module at T-point (=|U10m-Uoce|)  [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr               !: sea heat flux:     solar                     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns    , qns_b    !: sea heat flux: non solar                     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr_tot           !: total     solar heat flux (over sea and ice) [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns_tot           !: total non solar heat flux (over sea and ice) [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp    , emp_b    !: freshwater budget: volume flux               [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx    , sfx_b    !: salt flux                                    [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp_tot           !: total E-P over ocean and ice                 [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fmmflx            !: freshwater budget: freezing/melting          [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rnf    , rnf_b    !: river runoff        [Kg/m2/s]  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fwfisf , fwfisf_b !: ice shelf melting   [Kg/m2/s]  
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  sbc_tsc, sbc_tsc_b  !: sbc content trend                      [K.m/s] jpi,jpj,jpts
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  qsr_hc , qsr_hc_b   !: heat content trend due to qsr flux     [K.m/s] jpi,jpj,jpk
   !!
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tprecip           !: total precipitation                          [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sprecip           !: solid precipitation                          [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fr_i              !: ice fraction = 1 - lead fraction      (between 0 to 1)
#if defined key_cpl_carbon_cycle
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   atm_co2           !: atmospheric pCO2                             [ppm]
#endif
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: xcplmask          !: coupling mask for ln_mixcpl (warning: allocated in sbccpl)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   greenland_icesheet_mass_array, greenland_icesheet_mask
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   antarctica_icesheet_mass_array, antarctica_icesheet_mask

   !!----------------------------------------------------------------------
   !!                     Sea Surface Mean fields
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC                     ::   nn_fsbc   !: frequency of sbc computation (as well as sea-ice model)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssu_m     !: mean (nn_fsbc time-step) surface sea i-current (U-point) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssv_m     !: mean (nn_fsbc time-step) surface sea j-current (V-point) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sst_m     !: mean (nn_fsbc time-step) surface sea temperature     [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sss_m     !: mean (nn_fsbc time-step) surface sea salinity            [psu]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssh_m     !: mean (nn_fsbc time-step) sea surface height                [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e3t_m     !: mean (nn_fsbc time-step) sea surface layer thickness       [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   frq_m     !: mean (nn_fsbc time-step) fraction of solar net radiation absorbed in the 1st T level [-]
   
   !!----------------------------------------------------------------------
   !!  Surface scalars of total ice sheet mass for Greenland and Antarctica, 
   !! passed from atmosphere to be converted to dvol and hence a freshwater 
   !! flux  by using old values. New values are saved in the dump, to become
   !! old values next coupling timestep. Freshwater fluxes split between 
   !! sub iceshelf melting and iceberg calving, scalled to flux per second
   !!----------------------------------------------------------------------
   
   REAL(wp), PUBLIC  :: greenland_icesheet_mass, greenland_icesheet_mass_rate_of_change, greenland_icesheet_timelapsed 
   REAL(wp), PUBLIC  :: antarctica_icesheet_mass, antarctica_icesheet_mass_rate_of_change, antarctica_icesheet_timelapsed

   ! sbccpl namelist parameters associated with icesheet freshwater input code. Included here rather than in sbccpl.F90 to 
   ! avoid circular dependencies.
   INTEGER, PUBLIC     ::   nn_coupled_iceshelf_fluxes     ! =0 : total freshwater input from iceberg calving and ice shelf basal melting 
                                                           ! taken from climatologies used (no action in coupling routines).
                                                           ! =1 :  use rate of change of mass of Greenland and Antarctic icesheets to set the 
                                                           ! combined magnitude of the iceberg calving and iceshelf melting freshwater fluxes.
                                                           ! =2 :  specify constant freshwater inputs in this namelist to set the combined
                                                           ! magnitude of iceberg calving and iceshelf melting freshwater fluxes.
   LOGICAL, PUBLIC     ::   ln_iceshelf_init_atmos         ! If true force ocean to initialise iceshelf masses from atmospheric values rather
                                                           ! than values in ocean restart (applicable if nn_coupled_iceshelf_fluxes=1).
   REAL(wp), PUBLIC    ::   rn_greenland_total_fw_flux    ! Constant total rate of freshwater input (kg/s) for Greenland (if nn_coupled_iceshelf_fluxes=2) 
   REAL(wp), PUBLIC    ::   rn_greenland_calving_fraction  ! Set fraction of total freshwater flux for iceberg calving - remainder goes to iceshelf melting.
   REAL(wp), PUBLIC    ::   rn_antarctica_total_fw_flux   ! Constant total rate of freshwater input (kg/s) for Antarctica (if nn_coupled_iceshelf_fluxes=2) 
   REAL(wp), PUBLIC    ::   rn_antarctica_calving_fraction ! Set fraction of total freshwater flux for iceberg calving - remainder goes to iceshelf melting.
   REAL(wp), PUBLIC    ::   rn_iceshelf_fluxes_tolerance   ! Absolute tolerance for detecting differences in icesheet masses. 

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_oce_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION sbc_oce_alloc  ***
      !!---------------------------------------------------------------------
      INTEGER :: ierr(5)
      !!---------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( utau(jpi,jpj) , utau_b(jpi,jpj) , taum(jpi,jpj) ,     &
         &      vtau(jpi,jpj) , vtau_b(jpi,jpj) , wndm(jpi,jpj) , STAT=ierr(1) ) 
         !
      ALLOCATE( qns_tot(jpi,jpj) , qns  (jpi,jpj) , qns_b(jpi,jpj),        &
         &      qsr_tot(jpi,jpj) , qsr  (jpi,jpj) ,                        &
         &      emp    (jpi,jpj) , emp_b(jpi,jpj) ,                        &
         &      sfx    (jpi,jpj) , sfx_b(jpi,jpj) , emp_tot(jpi,jpj), fmmflx(jpi,jpj), STAT=ierr(2) )
         !
      ALLOCATE( fwfisf  (jpi,jpj), rnf  (jpi,jpj) , sbc_tsc  (jpi,jpj,jpts) , qsr_hc  (jpi,jpj,jpk) ,     &
         &      fwfisf_b(jpi,jpj), rnf_b(jpi,jpj) , sbc_tsc_b(jpi,jpj,jpts) , qsr_hc_b(jpi,jpj,jpk) , STAT=ierr(3) )
         !
      ALLOCATE( tprecip(jpi,jpj) , sprecip(jpi,jpj) , fr_i(jpi,jpj) ,     &
#if defined key_cpl_carbon_cycle
         &      atm_co2(jpi,jpj) ,                                        &
#endif
         &      ssu_m  (jpi,jpj) , sst_m(jpi,jpj) , frq_m(jpi,jpj) ,      &
         &      ssv_m  (jpi,jpj) , sss_m(jpi,jpj) , ssh_m(jpi,jpj) , STAT=ierr(4) )
      ALLOCATE( greenland_icesheet_mass_array(jpi,jpj) , antarctica_icesheet_mass_array(jpi,jpj) )
      ALLOCATE( greenland_icesheet_mask(jpi,jpj) , antarctica_icesheet_mask(jpi,jpj) )
         !
#if defined key_vvl
      ALLOCATE( e3t_m(jpi,jpj) , STAT=ierr(5) )
#endif
         !
      sbc_oce_alloc = MAXVAL( ierr )
      IF( lk_mpp            )   CALL mpp_sum ( sbc_oce_alloc )
      IF( sbc_oce_alloc > 0 )   CALL ctl_warn('sbc_oce_alloc: allocation of arrays failed')
      !
   END FUNCTION sbc_oce_alloc


   SUBROUTINE sbc_tau2wnd
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_tau2wnd  ***
      !!                   
      !! ** Purpose : Estimation of wind speed as a function of wind stress   
      !!
      !! ** Method  : |tau|=rhoa*Cd*|U|^2
      !!---------------------------------------------------------------------
      USE dom_oce         ! ocean space and time domain
      USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, ztau, zcoef ! temporary variables
      INTEGER  ::   ji, jj                ! dummy indices
      !!---------------------------------------------------------------------
      zcoef = 0.5 / ( zrhoa * zcdrag ) 
!CDIR NOVERRCHK
      DO jj = 2, jpjm1
!CDIR NOVERRCHK
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            ztx = utau(ji-1,jj  ) + utau(ji,jj) 
            zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
            ztau = SQRT( ztx * ztx + zty * zty )
            wndm(ji,jj) = SQRT ( ztau * zcoef ) * tmask(ji,jj,1)
         END DO
      END DO
      CALL lbc_lnk( wndm(:,:) , 'T', 1. )
      !
   END SUBROUTINE sbc_tau2wnd

   !!======================================================================
END MODULE sbc_oce
