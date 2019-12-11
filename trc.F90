MODULE trc
   !!======================================================================
   !!                      ***  MODULE  trc  ***
   !! Passive tracers   :  module for tracers defined
   !!======================================================================
   !! History :   OPA  !  1996-01  (M. Levy)  Original code
   !!              -   !  2000-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!   NEMO      1.0  !  2004-03  (C. Ethe)  Free form and module
   !!             3.6  !  2016-11  (A. Yool)  Updated diags for CMIP6
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc
   
   IMPLICIT NONE
   PUBLIC

   PUBLIC   trc_alloc   ! called by nemogcm.F90

   !! parameters for the control of passive tracers
   !! ---------------------------------------------   
   INTEGER, PUBLIC                                                 ::   numnat_ref = -1   !: logical unit for the reference passive tracer namelist_top_ref
   INTEGER, PUBLIC                                                 ::   numnat_cfg = -1   !: logical unit for the reference passive tracer namelist_top_cfg
   INTEGER, PUBLIC                                                 ::   numont     = -1   !: logical unit for the reference passive tracer namelist output output.namelist.top
   INTEGER, PUBLIC                                                 ::   numstr     = -1   !: logical unit for tracer statistics
   INTEGER, PUBLIC                                                 ::   numrtr        !: logical unit for trc restart (read )
   INTEGER, PUBLIC                                                 ::   numrtw        !: logical unit for trc restart ( write )

   !! passive tracers fields (before,now,after)
   !! --------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)               ::  trai           !: initial total tracer
# if defined key_medusa && key_roam 
   !! AXY (17/11/2017): elemental cycle initial totals
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)               ::  cycletot       !: initial elemental cycle total
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)               ::  cycletot2      !: initial elemental cycle total excl. halo in mpp_sum
# endif
   REAL(wp), PUBLIC                                                ::  areatot        !: total volume 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  cvol           !: volume correction -degrad option- 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  trn            !: tracer concentration for now time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  tra            !: tracer concentration for next time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)         ::  trb            !: tracer concentration for before time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  sbc_trc_b      !: Before sbc fluxes for tracers
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  sbc_trc        !: Now sbc fluxes for tracers

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  trc_i          !: prescribed tracer concentration in sea ice for SBC
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  )         ::  trc_o          !: prescribed tracer concentration in ocean for SBC
   INTEGER             , PUBLIC                                    ::  nn_ice_tr      !: handling of sea ice tracers

   !! interpolated gradient
   !!--------------------------------------------------  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtru           !: hor. gradient at u-points at bottom ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtrv           !: hor. gradient at v-points at bottom ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtrui          !: hor. gradient at u-points at top    ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)           ::  gtrvi          !: hor. gradient at v-points at top    ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)             ::  qsr_mean        !: daily mean qsr
   
   !! passive tracers  (input and output)
   !! ------------------------------------------  
   LOGICAL             , PUBLIC                                    ::  ln_rsttr       !: boolean term for restart i/o for passive tracers (namelist)
   LOGICAL             , PUBLIC                                    ::  lrst_trc       !: logical to control the trc restart write
   INTEGER             , PUBLIC                                    ::  nn_writetrc    !: time step frequency for concentration outputs (namelist)
   INTEGER             , PUBLIC                                    ::  nutwrs         !: output FILE for passive tracers restart
   INTEGER             , PUBLIC                                    ::  nutrst         !: logical unit for restart FILE for passive tracers
   INTEGER             , PUBLIC                                    ::  nn_rsttr       !: control of the time step ( 0 or 1 ) for pass. tr.
   CHARACTER(len = 80) , PUBLIC                                    ::  cn_trcrst_in   !: suffix of pass. tracer restart name (input)
   CHARACTER(len = 256), PUBLIC                                    ::  cn_trcrst_indir  !: restart input directory
   CHARACTER(len = 80) , PUBLIC                                    ::  cn_trcrst_out  !: suffix of pass. tracer restart name (output)
   CHARACTER(len = 256), PUBLIC                                    ::  cn_trcrst_outdir  !: restart output directory
   REAL(wp)            , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::  rdttrc         !: vertical profile of passive tracer time step
   LOGICAL             , PUBLIC                                    ::  ln_top_euler  !: boolean term for euler integration 
   LOGICAL             , PUBLIC                                    ::  ln_trcdta      !: Read inputs data from files
   LOGICAL             , PUBLIC                                    ::  ln_trcdmp      !: internal damping flag
   LOGICAL             , PUBLIC                                    ::  ln_trcdmp_clo  !: internal damping flag on closed seas
   INTEGER             , PUBLIC                                    ::  nittrc000       !: first time step of passive tracers model
   LOGICAL             , PUBLIC                                    ::  l_trcdm2dc     !: Diurnal cycle for TOP

   !! Information for the ice module for tracers
   !! ------------------------------------------
   TYPE TRC_I_NML                    !--- Ice tracer namelist structure
         REAL(wp)         :: trc_ratio  ! ice-ocean trc ratio
         REAL(wp)         :: trc_prescr ! prescribed ice trc cc
         CHARACTER(len=2) :: ctrc_o     ! choice of ocean trc cc
   END TYPE

   REAL(wp), DIMENSION(jptra), PUBLIC         :: trc_ice_ratio, & ! ice-ocean tracer ratio
                                                 trc_ice_prescr   ! prescribed ice trc cc
   CHARACTER(len=2), DIMENSION(jptra), PUBLIC :: cn_trc_o ! choice of ocean tracer cc

   !! information for outputs
   !! --------------------------------------------------
   TYPE, PUBLIC :: PTRACER                                                            !: Passive tracer type
       CHARACTER(len = 20)  :: clsname  !: short name
       CHARACTER(len = 80)  :: cllname  !: long name
       CHARACTER(len = 20)  :: clunit   !: unit
       LOGICAL              :: llinit   !: read in a file or not
       LOGICAL              :: llsave   !: save the tracer or not
   END TYPE PTRACER
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ctrcnm         !: tracer name 
   CHARACTER(len = 80), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ctrcln         !: trccer field long name
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ctrcun         !: tracer unit
   LOGICAL            , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_wri     !: save the tracer or not

   TYPE, PUBLIC :: DIAG                                                               !: passive trcacer ddditional diagnostic type
      CHARACTER(len = 20)  :: sname    !: short name
      CHARACTER(len = 80)  :: lname    !: long name
      CHARACTER(len = 20)  :: units    !: unit
   END TYPE DIAG

#if defined key_medusa
   TYPE, PUBLIC :: BDIAG
      LOGICAL              :: dgsave
   END TYPE BDIAG
   
   TYPE, PUBLIC :: DIAG_IOM
      TYPE(BDIAG) INVTN, INVTSI, INVTFE, PRN, MPN, PRD, MPD, DSED, OPAL, OPALDISS, GMIPn,            &
                  GMID, MZMI, GMEPN, GMEPD, GMEZMI, GMED, MZME, DEXP, DETN, MDET, AEOLIAN, BENTHIC,  &
                  SCAVENGE, PN_JLIM, PN_NLIM, PN_FELIM, PD_JLIM, PD_NLIM, PD_FELIM, PD_SILIM,        &
                  PDSILIM2, SDT__100, SDT__200, SDT__500, SDT_1000, TOTREG_N, TOTRG_SI, REG__100,    &
                  REG__200, REG__500, REG_1000, FASTN, FASTSI, FASTFE, FASTC, FASTCA, FDT__100,      &
                  FDT__200, FDT__500, FDT_1000, RG__100F, RG__200F, RG__500F, RG_1000F, FDS__100,    &
                  FDS__200, FDS__500, FDS_1000, RGS_100F, RGS_200F, RGS_500F, RGS1000F, REMINN,      &
                  REMINSI, REMINFE, REMINC, REMINCA, SEAFLRN, SEAFLRSI, SEAFLRFE, SEAFLRC, SEAFLRCA, &
                  MED_QSR, MED_XPAR, INTFLX_N, INTFLX_SI, INTFLX_FE, INT_PN, INT_PD, ML_PRN, ML_PRD, &
                  OCAL_CCD, OCAL_LVL, FE_0000, FE_0100, FE_0200, FE_0500, FE_1000, MED_XZE, WIND,    &
                  ATM_PCO2, OCN_PH, OCN_PCO2, OCNH2CO3, OCN_HCO3, OCN_CO3, CO2FLUX, OM_CAL, OM_ARG,  &
                  TCO2, TALK, KW660, ATM_PP0, O2FLUX, O2SAT, CAL_CCD, ARG_CCD, SFR_OCAL, SFR_OARG,   &
                  N_PROD, N_CONS, C_PROD, C_CONS, O2_PROD, O2_CONS, O2_ANOX, RR_0100, RR_0500,       &
                  RR_1000, IBEN_N, IBEN_FE, IBEN_C, IBEN_SI, IBEN_CA, OBEN_N, OBEN_FE, OBEN_C,       &
                  OBEN_SI, OBEN_CA, BEN_N, BEN_FE, BEN_C, BEN_SI, BEN_CA, RUNOFF, RIV_N, RIV_SI,     &
                  RIV_C, RIV_ALK, DETC, SDC__100, SDC__200, SDC__500, SDC_1000, INVTC, INVTALK,      &
                  INVTO2, LYSO_CA, COM_RESP, PN_LLOSS, PD_LLOSS, ZI_LLOSS, ZE_LLOSS, ZI_MES_N,       &
                  ZI_MES_D, ZI_MES_C, ZI_MESDC, ZI_EXCR, ZI_RESP, ZI_GROW, ZE_MES_N, ZE_MES_D,       &
                  ZE_MES_C, ZE_MESDC, ZE_EXCR, ZE_RESP, ZE_GROW, MDETC, GMIDC, GMEDC,                &
                  INT_ZMI, INT_ZME, INT_DET, INT_DTC, DMS_SURF, DMS_ANDR, DMS_SIMO, DMS_ARAN,        &
                  DMS_HALL, DMS_ANDM, ATM_XCO2, OCN_FCO2, ATM_FCO2, OCN_RHOSW, OCN_SCHCO2,           &
                  OCN_KWCO2, OCN_K0, CO2STARAIR, OCN_DPCO2,                                          & ! end of regular 2D
                  TPP3, DETFLUX3, REMIN3N, PH3, OM_CAL3,                                             & ! end of regular 3D
! JPALM (01/09/17): additional UKESM 2D diag
                  CHL_MLD, CHL_CPL,                                                                  &
! AXY (11/11/16): additional CMIP6 2D diagnostics
                  epC100, epCALC100, epN100, epSI100,                                                &
                  FGCO2, INTDISSIC, INTDISSIN, INTDISSISI, INTTALK, O2min, ZO2min,                   &
                  FBDDTALK, FBDDTDIC, FBDDTDIFE, FBDDTDIN, FBDDTDISI,                                & 
! AXY (11/11/16): additional CMIP6 3D diagnostics
                  TPPD3,                                                                             &
                  BDDTALK3, BDDTDIC3, BDDTDIFE3, BDDTDIN3, BDDTDISI3,                                & 
                  FD_NIT3, FD_SIL3, FD_CAR3, FD_CAL3,                                                & 
                  CO33, CO3SATARAG3, CO3SATCALC3, DCALC3,                                            &
                  EXPC3, EXPN3, EXPCALC3, EXPSI3,                                                    &
                  FEDISS3, FESCAV3,                                                                  &
                  MIGRAZP3, MIGRAZD3, MEGRAZP3, MEGRAZD3, MEGRAZZ3,                                  &
                  O2SAT3, PBSI3, PCAL3, REMOC3,                                                      &
                  PNLIMJ3, PNLIMN3, PNLIMFE3, PDLIMJ3, PDLIMN3, PDLIMFE3, PDLIMSI3       
                  !!
                  !! list of all MEDUSA diagnostics that could be called by iom_use
   END TYPE DIAG_IOM  
   !!
   TYPE(DIAG_IOM), PUBLIC :: med_diag  ! define which diagnostics are asked in outputs
# endif                   

   !! information for inputs
   !! --------------------------------------------------
   LOGICAL            , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_ini     !: Initialisation from data input file
   LOGICAL            , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_obc     !: Use open boundary condition data
   LOGICAL            , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_sbc     !: Use surface boundary condition data
   LOGICAL            , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    ::  ln_trc_cbc     !: Use coastal boundary condition data

   !! additional 2D/3D outputs namelist
   !! --------------------------------------------------
   REAL(wp)           , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,  :) ::   trc2d         !: additional 2d outputs array 
   REAL(wp)           , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   trc3d         !: additional 3d outputs array 
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::   ctrc2d        !: 2d field short name
   CHARACTER(len = 80), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::   ctrc2l        !: 2d field long name
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::   ctrc2u        !: 2d field unit
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::   ctrc3d        !: 3d field short name
   CHARACTER(len = 80), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::   ctrc3l        !: 3d field long name
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::   ctrc3u        !: 3d field unit
   LOGICAL            , PUBLIC                                        ::  ln_diatrc      !: boolean term for additional diagnostic
   INTEGER            , PUBLIC                                        ::  nn_writedia    !: frequency of additional outputs

   !! Biological trends
   !! -----------------
   LOGICAL            , PUBLIC                                        ::  ln_diabio      !: boolean term for biological diagnostic
   INTEGER            , PUBLIC                                        ::  nn_writebio    !: frequency of biological outputs
   REAL(wp)           , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  trbio          !: biological trends
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::  ctrbio         !: bio field short name
   CHARACTER(len = 80), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::  ctrbil         !: bio field long name
   CHARACTER(len = 20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::  ctrbiu         !: bio field unit

   !! variables to average over physics over passive tracer sub-steps.
   !! ----------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  un_tm       !: i-horizontal velocity average     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  vn_tm       !: j-horizontal velocity average     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  tsn_tm      !: t/s average     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avt_tm      !: vertical diffusivity coeff. at  w-point   [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  rhop_tm     !: 
# if defined key_zdfddm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avs_tm      !: vertical double diffusivity coeff. at w-point   [m/s]
# endif
#if defined key_ldfslp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  wslpi_tm    !: i-direction slope at u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  wslpj_tm    !: j-direction slope at u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  uslp_tm     !: j-direction slope at u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  vslp_tm     !: j-direction slope at u-, w-points
#endif
#if defined key_trabbl
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  ahu_bbl_tm  !: u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  ahv_bbl_tm  !: j-direction slope at u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  utr_bbl_tm  !: j-direction slope at u-, w-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  vtr_bbl_tm  !: j-direction slope at u-, w-points
#endif
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  sshn_tm     !: average ssh for the now step [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  sshb_hold   !:hold sshb from the beginning of each sub-stepping[m]  

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  rnf_tm     !: river runoff
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  h_rnf_tm   !: depth in metres to the bottom of the relevant grid box
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  hmld_tm    !: mixed layer depth average [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  fr_i_tm    !: average ice fraction     [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  emp_tm     !: freshwater budget: volume flux [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  fmmflx_tm  !: freshwater budget: freezing/melting [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  emp_b_hold !: hold emp from the beginning of each sub-stepping[m]  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  qsr_tm     !: solar radiation average [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  wndm_tm    !: 10m wind average [m]
   !

   ! Temporary physical arrays for sub_stepping
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  tsn_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  un_temp,vn_temp,wn_temp     !: hold current values of avt, un, vn, wn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avt_temp, rhop_temp     !: hold current values of avt, un, vn, wn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  sshn_temp, sshb_temp, ssha_temp, rnf_temp,h_rnf_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  hdivn_temp, rotn_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  hdivb_temp, rotb_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  hmld_temp, qsr_temp, fr_i_temp,wndm_temp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  emp_temp, fmmflx_temp, emp_b_temp
   !
#if defined key_trabbl
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  ahu_bbl_temp, ahv_bbl_temp, utr_bbl_temp, vtr_bbl_temp !: hold current values 
#endif
   !
#if defined key_ldfslp
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  wslpi_temp, wslpj_temp, uslp_temp, vslp_temp    !: hold current values 
#endif
   ! 
# if defined key_zdfddm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  avs_temp      !: salinity vertical diffusivity coeff. at w-point   [m/s]
# endif
   !

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3.1 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE trc_alloc ***
      !!-------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_warn
      !!-------------------------------------------------------------------
      !
      ALLOCATE( trn(jpi,jpj,jpk,jptra), trb(jpi,jpj,jpk,jptra), tra(jpi,jpj,jpk,jptra),       &  
         &      trc_i(jpi,jpj,jptra)  , trc_o(jpi,jpj,jptra)                          ,       &
         &      gtru (jpi,jpj,jptra)  , gtrv (jpi,jpj,jptra)                          ,       &
         &      gtrui(jpi,jpj,jptra)  , gtrvi(jpi,jpj,jptra)                          ,       &
         &      sbc_trc_b(jpi,jpj,jptra), sbc_trc(jpi,jpj,jptra)                      ,       &  
         &      cvol(jpi,jpj,jpk)     , rdttrc(jpk)           , trai(jptra)           ,       &
         &      ctrcnm(jptra)         , ctrcln(jptra)         , ctrcun(jptra)         ,       & 
# if defined key_medusa && defined key_roam
         &      cycletot(6), cycletot2(6)                                             ,       &
# endif
         &      ln_trc_ini(jptra)     , ln_trc_wri(jptra)     , qsr_mean(jpi,jpj)     ,  STAT = trc_alloc  )  

      IF( trc_alloc /= 0 )   CALL ctl_warn('trc_alloc: failed to allocate arrays')

      ! It is known that not intialising SBC_TRC can introduce NaNs
      sbc_trc(:,:,:) = 0.0

      !
   END FUNCTION trc_alloc

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE trc
