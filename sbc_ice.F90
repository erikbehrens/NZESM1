MODULE sbc_ice
   !!======================================================================
   !!                 ***  MODULE  sbc_ice  ***
   !! Surface module - LIM-3: parameters & variables defined in memory
   !!======================================================================
   !! History :  3.0  ! 2006-08  (G. Madec)  Surface module
   !!            3.2  ! 2009-06  (S. Masson) merge with ice_oce
   !!            3.3.1! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.4  ! 2011-11  (C. Harris) CICE added as an option
   !!----------------------------------------------------------------------
#if defined key_lim3 || defined key_lim2 || defined key_cice
   !!----------------------------------------------------------------------
   !!   'key_lim2' or 'key_lim3' :             LIM-2 or LIM-3 sea-ice model
   !!----------------------------------------------------------------------
   USE par_oce          ! ocean parameters
   USE sbc_oce          ! surface boundary condition: ocean
# if defined key_lim3
   USE ice              ! LIM-3 parameters
# endif
# if defined key_lim2
   USE par_ice_2        ! LIM-2 parameters
   USE ice_2
# endif
# if defined key_cice
   USE ice_domain_size, only: ncat 
#endif
   USE lib_mpp          ! MPP library
   USE in_out_manager   ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_ice_alloc ! called in iceini(_2).F90

# if defined  key_lim2
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim2    = .TRUE.   !: LIM-2 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim3    = .FALSE.  !: no LIM-3
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .FALSE.  !: no CICE 
#  if defined key_lim2_vp
   CHARACTER(len=1), PUBLIC, PARAMETER ::   cp_ice_msh = 'I'      !: VP : 'I'-grid ice-velocity (B-grid lower left corner)
#  else
   CHARACTER(len=1), PUBLIC, PARAMETER ::   cp_ice_msh = 'C'      !: EVP: 'C'-grid ice-velocity
#  endif
# endif
# if defined  key_lim3
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim2    = .FALSE.  !: no LIM-2
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim3    = .TRUE.   !: LIM-3 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .FALSE.  !: no CICE 
   CHARACTER(len=1), PUBLIC, PARAMETER ::   cp_ice_msh = 'C'      !: 'C'-grid ice-velocity
# endif
# if defined  key_cice
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim2    = .FALSE.  !: no LIM-2
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim3    = .FALSE.  !: no LIM-3
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .TRUE.   !: CICE ice model
   CHARACTER(len=1), PUBLIC            ::   cp_ice_msh = 'F'      !: 'F'-grid ice-velocity
# endif

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qns_ice        !: non solar heat flux over ice                  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice        !: solar heat flux over ice                      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qla_ice        !: latent flux over ice                          [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dqla_ice       !: latent sensibility over ice                 [W/m2/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dqns_ice       !: non solar heat flux over ice (LW+SEN+LA)    [W/m2/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tn_ice         !: ice surface temperature                          [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   alb_ice        !: ice albedo                                       [-]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   utau_ice       !: atmos-ice u-stress. VP: I-pt ; EVP: U,V-pts   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   vtau_ice       !: atmos-ice v-stress. VP: I-pt ; EVP: U,V-pts   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr1_i0         !: Solar surface transmission parameter, thick ice  [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr2_i0         !: Solar surface transmission parameter, thin ice   [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ice        !: sublimation - precip over sea ice          [kg/m2/s]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   topmelt            !: category topmelt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   botmelt            !: category botmelt

#if defined  key_lim3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   evap_ice       !: sublimation                              [kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   devap_ice      !: sublimation sensitivity                [kg/m2/s/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qns_oce        !: non solar heat flux over ocean              [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qsr_oce        !: non solar heat flux over ocean              [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qemp_oce       !: heat flux of precip and evap over ocean     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qemp_ice       !: heat flux of precip and evap over ice       [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qevap_ice      !: heat flux of evap over ice                  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qprec_ice      !: enthalpy of precip over ice                 [J/m3]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_oce        !: evap - precip over ocean                 [kg/m2/s]
#endif
#if defined key_lim3 || defined key_lim2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wndm_ice       !: wind speed module at T-point                 [m/s]
#endif

#if defined key_cice
   !
   ! for consistency with LIM, these are declared with three dimensions
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qlw_ice            !: incoming long-wave
   !
   ! other forcing arrays are two dimensional
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ss_iou             !: x ice-ocean surface stress at NEMO U point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ss_iov             !: y ice-ocean surface stress at NEMO V point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qatm_ice           !: specific humidity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wndi_ice           !: i wind at T point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wndj_ice           !: j wind at T point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   nfrzmlt            !: NEMO frzmlt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr_iu              !: ice fraction at NEMO U point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr_iv              !: ice fraction at NEMO V point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sstfrz             !: sea surface freezing temperature
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tsfc_ice           !: sea-ice surface skin temperature (on categories)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   kn_ice             !: sea-ice surface layer thermal conductivity (on cats)

   ! variables used in the coupled interface
   INTEGER , PUBLIC, PARAMETER ::   jpl = ncat
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   u_ice, v_ice          ! jpi, jpj
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  a_p, ht_p ! Meltpond fraction and depth
   
   !
   
   !
#if defined key_asminc
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ndaice_da          !: NEMO fresh water flux to ocean due to data assim
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   nfresh_da          !: NEMO salt flux to ocean due to data assim
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   nfsalt_da          !: NEMO ice concentration change/second from data assim
#endif
      
#endif
   
#if defined key_lim2 || defined key_cice
   ! already defined in ice.F90 for LIM3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  a_i
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  ht_i, ht_s
#endif

#if defined key_cice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   tatm_ice       !: air temperature [K]
#endif

   REAL(wp), PUBLIC, SAVE ::   cldf_ice = 0.81    !: cloud fraction over sea ice, summer CLIO value   [-]

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_ice_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  FUNCTION sbc_ice_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(5)
      !!----------------------------------------------------------------------
      ierr(:) = 0

#if defined key_lim3 || defined key_lim2
      ALLOCATE( qns_ice (jpi,jpj,jpl) , qsr_ice (jpi,jpj,jpl) ,     &
         &      qla_ice (jpi,jpj,jpl) , dqla_ice(jpi,jpj,jpl) ,     &
         &      dqns_ice(jpi,jpj,jpl) , tn_ice  (jpi,jpj,jpl) , alb_ice (jpi,jpj,jpl) ,   &
         &      utau_ice(jpi,jpj)     , vtau_ice(jpi,jpj)     , wndm_ice(jpi,jpj)     ,   &
         &      fr1_i0  (jpi,jpj)     , fr2_i0  (jpi,jpj)     ,     &
#if defined key_lim2
         &      a_i(jpi,jpj,jpl)      ,                             &
#endif
#if defined key_lim3
         &      evap_ice(jpi,jpj,jpl) , devap_ice(jpi,jpj,jpl) , qprec_ice(jpi,jpj) ,   &
         &      qemp_ice(jpi,jpj)     , qevap_ice(jpi,jpj,jpl) , qemp_oce (jpi,jpj) ,   &
         &      qns_oce (jpi,jpj)     , qsr_oce  (jpi,jpj)     , emp_oce (jpi,jpj)  ,   &
#endif
         &      emp_ice(jpi,jpj)      ,  STAT= ierr(1) )
#endif

#if defined key_cice
      ALLOCATE( qla_ice(jpi,jpj,ncat) , qlw_ice(jpi,jpj,1)    , qsr_ice(jpi,jpj,1)    , &
                wndi_ice(jpi,jpj)     , tatm_ice(jpi,jpj)     , qatm_ice(jpi,jpj)     , &
                wndj_ice(jpi,jpj)     , nfrzmlt(jpi,jpj)      , ss_iou(jpi,jpj)       , &
                ss_iov(jpi,jpj)       , fr_iu(jpi,jpj)        , fr_iv(jpi,jpj)        , &
                a_i(jpi,jpj,ncat)     , topmelt(jpi,jpj,ncat) , botmelt(jpi,jpj,ncat) , &
#if defined key_asminc
                ndaice_da(jpi,jpj)    , nfresh_da(jpi,jpj)    , nfsalt_da(jpi,jpj)    , &
#endif
                sstfrz(jpi,jpj)       , STAT= ierr(1) )
   ! Alex West: Allocating tn_ice with 5 categories.  When NEMO is used with CICE, this variable
   ! represents top layer ice temperature, which is multi-category.
      IF( ln_cpl )   ALLOCATE( u_ice(jpi,jpj)        , fr1_i0(jpi,jpj)       , tn_ice (jpi,jpj,jpl)  , &
         &                     v_ice(jpi,jpj)        , fr2_i0(jpi,jpj)       , alb_ice(jpi,jpj,1)    , &
         &                     emp_ice(jpi,jpj)      , qns_ice(jpi,jpj,1)    , dqns_ice(jpi,jpj,1)   , &
         &                     a_p(jpi,jpj,jpl)      , ht_p(jpi,jpj,jpl)     , tsfc_ice(jpi,jpj,jpl) , &
         &                     kn_ice(jpi,jpj,jpl) ,    STAT=ierr(2) )
      
#endif
         !
#if defined key_cice || defined key_lim2
      IF( ln_cpl )   ALLOCATE( ht_i(jpi,jpj,jpl) , ht_s(jpi,jpj,jpl) , STAT=ierr(5) )
#endif

      sbc_ice_alloc = MAXVAL( ierr )
      IF( lk_mpp            )   CALL mpp_sum ( sbc_ice_alloc )
      IF( sbc_ice_alloc > 0 )   CALL ctl_warn('sbc_ice_alloc: allocation of arrays failed')
   END FUNCTION sbc_ice_alloc

#else
   !!----------------------------------------------------------------------
   !!   Default option                      NO LIM 2.0 or 3.0 or CICE sea-ice model
   !!----------------------------------------------------------------------
   USE in_out_manager   ! I/O manager
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim2    = .FALSE.  !: no LIM-2 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim3    = .FALSE.  !: no LIM-3 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .FALSE.  !: no CICE  ice model
   CHARACTER(len=1), PUBLIC, PARAMETER ::   cp_ice_msh = '-'      !: no grid ice-velocity
   REAL            , PUBLIC, PARAMETER ::   cldf_ice = 0.81       !: cloud fraction over sea ice, summer CLIO value   [-]
   INTEGER         , PUBLIC, PARAMETER ::   jpl = 1 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   u_ice, v_ice,fr1_i0,fr2_i0          ! jpi, jpj
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tn_ice, alb_ice, qns_ice, dqns_ice  ! (jpi,jpj,jpl)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   a_i
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ht_i, ht_s
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   topmelt, botmelt
#endif

   !!======================================================================
END MODULE sbc_ice
