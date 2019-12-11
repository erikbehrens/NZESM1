MODULE sms_medusa 
   !!----------------------------------------------------------------------
   !!                     ***  sms_medusa.F90  ***  
   !! TOP :   MEDUSA  Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :    -   !  1999-09 (M. Levy)  original code
   !!              -   !  2000-12 (O. Aumont, E. Kestenare)  add sediment 
   !!             1.0  !  2005-10 (C. Ethe) F90
   !!             1.0  !  2005-03  (A-S Kremeur) add fphylab, fzoolab, fdetlab, fdbod
   !!              -   !  2005-06  (A-S Kremeur) add sedpocb, sedpocn, sedpoca
   !!             2.0  !  2007-04  (C. Deltel, G. Madec)  Free form and modules
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!              -   !  2011-04  (A. Yool) updated for ROAM project
   !!----------------------------------------------------------------------

#if defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'                                         MEDUSA model
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !! NEMO/TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   INTEGER ::   numnatp_ref = -1           !! Logical units for namelist medusa
   INTEGER ::   numnatp_cfg = -1           !! Logical units for namelist medusa
   INTEGER ::   numonp      = -1           !! Logical unit for namelist medusa output

!!----------------------------------------------------------------------
!! Biological parameters
!!----------------------------------------------------------------------
!!
!! Primary production and chl related quantities
   REAL(wp) ::  xxi       !:  conversion factor from gC to mmolN 
   REAL(wp) ::  xaln      !:  Chl-a specific initial slope of P-I curve for non-diatoms 
   REAL(wp) ::  xald      !:  Chl-a specific initial slope of P-I curve for diatoms
   INTEGER  ::  jphy      !:  phytoplankton T-dependent growth switch
   REAL(wp) ::  xvpn      !:  maximum growth rate for non-diatoms
   REAL(wp) ::  xvpd      !:  maximum growth rate for diatoms
   REAL(wp) ::  xthetam   !:  maximum Chl to C ratio for non-diatoms      
   REAL(wp) ::  xthetamd  !:  maximum Chl to C ratio for diatoms    
   REAL(wp) ::  jq10      !:  specific Q10 value (jphy==2)    
!!
!! Diatom silicon parameters
   REAL(wp) ::  xsin0     !:  minimum diatom Si:N ratio
   REAL(wp) ::  xnsi0     !:  minimum diatom N:Si ratio
   REAL(wp) ::  xuif      !:  hypothetical growth ratio at infinite Si:N ratio
!!
!! Nutrient limitation
   INTEGER  ::  jliebig   !:  Liebig nutrient uptake switch
   REAL(wp) ::  xnln      !:  half-sat constant for DIN uptake by non-diatoms 
   REAL(wp) ::  xnld      !:  half-sat constant for DIN uptake by diatoms 
   REAL(wp) ::  xsld      !:  half-sat constant for Si uptake by diatoms 
   REAL(wp) ::  xfln      !:  half-sat constant for Fe uptake by non-diatoms 
   REAL(wp) ::  xfld      !:  half-sat constant for Fe uptake by diatoms  
!!
!! Grazing
   REAL(wp) ::  xgmi      !:  microzoo maximum growth rate 
   REAL(wp) ::  xgme      !:  mesozoo maximum growth rate 
   REAL(wp) ::  xkmi      !:  microzoo grazing half-sat parameter
   REAL(wp) ::  xkme      !:  mesozoo grazing half-sat parameter
   REAL(wp) ::  xphi      !:  micro/mesozoo grazing inefficiency
   REAL(wp) ::  xbetan    !:  micro/mesozoo N assimilation efficiency
   REAL(wp) ::  xbetac    !:  micro/mesozoo C assimilation efficiency
   REAL(wp) ::  xkc       !:  micro/mesozoo net C growth efficiency
   REAL(wp) ::  xpmipn    !:  grazing preference of microzoo for non-diatoms
   REAL(wp) ::  xpmid     !:  grazing preference of microzoo for diatoms
   REAL(wp) ::  xpmepn    !:  grazing preference of mesozoo for non-diatoms 
   REAL(wp) ::  xpmepd    !:  grazing preference of mesozoo for diatoms
   REAL(wp) ::  xpmezmi   !:  grazing preference of mesozoo for microzoo
   REAL(wp) ::  xpmed     !:  grazing preference of mesozoo for detritus
!!
!! Metabolic losses
   REAL(wp) ::  xmetapn   !:  non-diatom metabolic loss rate
   REAL(wp) ::  xmetapd   !:  diatom     metabolic loss rate
   REAL(wp) ::  xmetazmi  !:  microzoo   metabolic loss rate
   REAL(wp) ::  xmetazme  !:  mesozoo    metabolic loss rate
!!
!! Mortality losses
   INTEGER  ::  jmpn      !:  non-diatom mortality functional form
   REAL(wp) ::  xmpn      !:  non-diatom mortality rate
   REAL(wp) ::  xkphn     !:  non-diatom mortality half-sat constant
   INTEGER  ::  jmpd      !:  diatom     mortality functional form
   REAL(wp) ::  xmpd      !:  diatom     mortality rate
   REAL(wp) ::  xkphd     !:  diatom     mortality half-sat constant
   INTEGER  ::  jmzmi     !:  microzoo   mortality functional form
   REAL(wp) ::  xmzmi     !:  microzoo   mortality rate
   REAL(wp) ::  xkzmi     !:  microzoo   mortality half-sat constant
   INTEGER  ::  jmzme     !:  mesozoo    mortality functional form
   REAL(wp) ::  xmzme     !:  mesozoo    mortality rate
   REAL(wp) ::  xkzme     !:  mesozoo    mortality half-sat constant
!!
!! Remineralisation
   INTEGER  ::  jmd       !:  detritus T-dependent remineralisation switch
   INTEGER  ::  jsfd      !:  accelerate seafloor detritus remin. switch
   REAL(wp) ::  xmd       !:  detrital nitrogen remineralisation rate 
   REAL(wp) ::  xmdc      !:  detrital carbon remineralisation rate 
!!
!! Stochiometric ratios
   REAL(wp) ::  xthetapn  !:  non-diatom C:N ratio
   REAL(wp) ::  xthetapd  !:  diatom C:N ratio
   REAL(wp) ::  xthetazmi !:  microzoo C:N ratio
   REAL(wp) ::  xthetazme !:  mesozoo C:N ratio
   REAL(wp) ::  xthetad   !:  detritus C:N ratio
   REAL(wp) ::  xrfn      !:  phytoplankton Fe:N ratio
   REAL(wp) ::  xrsn      !:  diatom Si:N ratio (NOT USED HERE; RETAINED FOR LOBSTER)
!!
!! Iron parameters
   INTEGER  ::  jiron     !:  iron scavenging submodel switch
   REAL(wp) ::  xfe_mass  !:  iron atomic mass
   REAL(wp) ::  xfe_sol   !:  aeolian iron solubility
   REAL(wp) ::  xfe_sed   !:  sediment iron input
   REAL(wp) ::  xLgT      !:  total ligand concentration (umol/m3)
   REAL(wp) ::  xk_FeL    !:  dissociation constant for (Fe + L)
   REAL(wp) ::  xk_sc_Fe  !:  scavenging rate of "free" iron
!!
!! Gravitational sinking      
   REAL(wp) ::  vsed      !:  detritus gravitational sinking rate 
   REAL(wp) ::  xhr       !:  coefficient for Martin et al. (1987) remineralisation
!!
!! Fast-sinking detritus parameters
   INTEGER  ::  jexport   !:  fast detritus remineralisation switch
   INTEGER  ::  jfdfate   !:  fate of fast detritus at seafloor switch
   INTEGER  ::  jrratio   !:  rain ratio switch
   INTEGER  ::  jocalccd  !:  CCD switch
   REAL(wp) ::  xridg_r0  !:  Ridgwell rain ratio coefficient
   REAL(wp) ::  xfdfrac1  !:  fast-sinking fraction of diatom nat. mort. losses
   REAL(wp) ::  xfdfrac2  !:  fast-sinking fraction of mesozooplankton mort. losses
   REAL(wp) ::  xfdfrac3  !:  fast-sinking fraction of diatom silicon grazing losses
   REAL(wp) ::  xcaco3a   !:  polar (high latitude) CaCO3 fraction
   REAL(wp) ::  xcaco3b   !:  equatorial (low latitude) CaCO3 fraction
   REAL(wp) ::  xmassc    !:  organic C mass:mole ratio, C106 H175 O40 N16 P1
   REAL(wp) ::  xmassca   !:  calcium carbonate mass:mole ratio, CaCO3
   REAL(wp) ::  xmasssi   !:  biogenic silicon mass:mole ratio, (H2SiO3)n
   REAL(wp) ::  xprotca   !:  calcium carbonate protection ratio
   REAL(wp) ::  xprotsi   !:  biogenic silicon protection ratio
   REAL(wp) ::  xfastc    !:  organic C remineralisation length scale
   REAL(wp) ::  xfastca   !:  calcium carbonate dissolution length scale
   REAL(wp) ::  xfastsi   !:  biogenic silicon dissolution length scale
!!
!! Benthos parameters
   INTEGER  ::  jorgben   !:  does   organic detritus go to the benthos?
   INTEGER  ::  jinorgben !:  does inorganic detritus go to the benthos?
!!
   REAL(wp) ::  xsedn     !:  organic   nitrogen sediment remineralisation rate 
   REAL(wp) ::  xsedfe    !:  organic   iron     sediment remineralisation rate 
   REAL(wp) ::  xsedsi    !:  inorganic silicon  sediment dissolution      rate 
   REAL(wp) ::  xsedc     !:  organic   carbon   sediment remineralisation rate 
   REAL(wp) ::  xsedca    !:  inorganic carbon   sediment dissolution      rate 
   REAL(wp) ::  xburial   !:  burial rate of seafloor detritus
!!
!! River parameters
   INTEGER  ::  jriver_n  !:  riverine nitrogen?   0 = no, 1 = conc, 2 = flux
   INTEGER  ::  jriver_si !:  riverine silicon?    0 = no, 1 = conc, 2 = flux
   INTEGER  ::  jriver_c  !:  riverine carbon?     0 = no, 1 = conc, 2 = flux
   INTEGER  ::  jriver_alk!:  riverine alkalinity? 0 = no, 1 = conc, 2 = flux
   INTEGER  ::  jriver_dep!:  depth river input added to?  1 = surface, >1 possible
!!
!! Miscellaneous
   REAL(wp) ::  xsdiss    !:  diatom frustule dissolution rate
!!
!! Additional parameters
   INTEGER  ::  jpkb      !:  vertical layer for diagnostic of the vertical flux 
!!
!! UKESM diagnostics
   INTEGER  ::  jdms         !: include DMS diagnostics ? Jpalm (27-08-2014) 
   INTEGER  ::  jdms_input   !: use instant (0) or diel-average (1) inputs (AXY, 08/07/2015)
   INTEGER  ::  jdms_model   !: choice of DMS model passed to atmosphere
!!                              1 = ANDR, 2 = SIMO, 3 = ARAN, 4 = HALL
!! JPALM --19-12-2017 -- UM people need to tune the Anderson DMS scheme    
   REAL(wp) ::  dmsmin       !: Anderson DMS scheme - DMS minimum value
   REAL(wp) ::  dmscut       !: Anderson DMS scheme - DMS cutoff value 
   REAL(wp) ::  dmsslp       !: Anderson DMS scheme - DMS slope 
!! FOR UKESM   
   REAL(wp) ::  scl_chl      !: scaling factor for tuned Chl passed to the UM 
   INTEGER  ::  chl_out      !: select Chl field exported and scaled for UM:
                             !: 1- Surface Chl ; 2- MLD Chl
!!
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   remdmp   !: depth dependent damping coefficient of passive tracers 
!!
!! AXY (27/07/10): add in indices for depth horizons (for sinking flux
!!                 and seafloor iron inputs)
   INTEGER  ::    i0100, i0150, i0200, i0500, i1000, i1100
#if defined key_roam
!!
!! ROAM carbon, alkalinity and oxygen cycle parameters
   REAL(wp) ::  xthetaphy    !:  oxygen evolution/consumption by phytoplankton 
   REAL(wp) ::  xthetazoo    !:  oxygen consumption by zooplankton 
   REAL(wp) ::  xthetanit    !:  oxygen consumption by nitrogen remineralisation
   REAL(wp) ::  xthetarem    !:  oxygen consumption by carbon remineralisation
   REAL(wp) ::  xo2min       !:  oxygen minimum concentration
!! 
!! 3D fields of carbonate system parameters
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: f3_pH       !: 3D pH
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: f3_h2co3    !: 3D carbonic acid
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: f3_hco3     !: 3D bicarbonate
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: f3_co3      !: 3D carbonate
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: f3_omcal    !: 3D omega calcite
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: f3_omarg    !: 3D omega aragonite
!! 
!! 2D fields of calcium carbonate compensation depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: f2_ccd_cal  !: 2D calcite CCD depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: f2_ccd_arg  !: 2D aragonite CCD depth
!!
!! 2D fields of organic and inorganic material sedimented on the seafloor
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_sed_n    !: 2D organic nitrogen   (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_sed_n    !: 2D organic nitrogen   (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_sed_n    !: 2D organic nitrogen   (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_sed_fe   !: 2D organic iron       (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_sed_fe   !: 2D organic iron       (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_sed_fe   !: 2D organic iron       (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_sed_si   !: 2D inorganic silicon  (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_sed_si   !: 2D inorganic silicon  (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_sed_si   !: 2D inorganic silicon  (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_sed_c    !: 2D organic carbon     (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_sed_c    !: 2D organic carbon     (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_sed_c    !: 2D organic carbon     (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_sed_ca   !: 2D inorganic carbon   (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_sed_ca   !: 2D inorganic carbon   (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_sed_ca   !: 2D inorganic carbon   (after)
!!
!! 2D fields of temporally averaged properties for DMS calculations (AXY, 07/07/15)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_dms_chn  !: 2D avg CHN   (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_dms_chn  !: 2D avg CHN   (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_dms_chn  !: 2D avg CHN   (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_dms_chd  !: 2D avg CHD   (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_dms_chd  !: 2D avg CHD   (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_dms_chd  !: 2D avg CHD   (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_dms_mld  !: 2D avg MLD   (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_dms_mld  !: 2D avg MLD   (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_dms_mld  !: 2D avg MLD   (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_dms_qsr  !: 2D avg QSR   (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_dms_qsr  !: 2D avg QSR   (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_dms_qsr  !: 2D avg QSR   (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_dms_din  !: 2D avg DIN   (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_dms_din  !: 2D avg DIN   (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_dms_din  !: 2D avg DIN   (after)
!!
!! 2D fields needing to be knows at first tstp for coupling with atm - UKEMS(Jpalm,14-06-2016)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_co2_flx  !: 2D avg fx co2 (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_co2_flx  !: 2D avg fx co2 (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_co2_flx  !: 2D avg fx co2 (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zb_dms_srf  !: 2D avg sfr dms (before)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_dms_srf  !: 2D avg sfr dms (now)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: za_dms_srf  !: 2D avg srf dms (after)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: zn_chl_srf  !: 2D avg srf chl (now)

#endif

!!----------------------------------------------------------------------
!! CCD parameter
!!----------------------------------------------------------------------
!!
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ocal_ccd  !: CCD depth

!!----------------------------------------------------------------------
!! Dust parameters
!!----------------------------------------------------------------------
!!
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   dust      !: dust parameter 1
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zirondep  !! Fe deposition

!!----------------------------------------------------------------------
!! River parameters
!!----------------------------------------------------------------------
!!
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   riv_n      !: riverine N
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   riv_si     !: riverine Si
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   riv_c      !: riverine C
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   riv_alk    !: riverine alkalinity
   !! AXY (19/07/12): add this to permit river fluxes to be added below top box
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   friver_dep !: where river fluxes added

#if defined key_roam
!!----------------------------------------------------------------------
!! JPALM -- change hist_pco2 init
!!----------------------------------------------------------------------
!!
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)   ::   hist_pco2 !: pCO2
   INTEGER  :: co2_rec
   REAL(wp) :: co2_yinit, co2_yend    !: First and Last year read in the xCO2.atm file
#endif

!!----------------------------------------------------------------------
!! JPALM -- PI CO2 key
!!----------------------------------------------------------------------
!!
#if defined key_axy_pi_co2
   LOGICAL , PUBLIC ::   lk_pi_co2 = .TRUE.  !: PI xCO2 used
#else
   LOGICAL , PUBLIC ::   lk_pi_co2 = .FALSE. !: PI xCO2 unused
#endif

!!----------------------------------------------------------------------
!! Optical parameters
!!----------------------------------------------------------------------
!!
   REAL(wp) ::   xkr0       !: water coefficient absorption in red      (NAMELIST)
   REAL(wp) ::   xkg0       !: water coefficient absorption in green    (NAMELIST)
   REAL(wp) ::   xkrp       !: pigment coefficient absorption in red    (NAMELIST)
   REAL(wp) ::   xkgp       !: pigment coefficient absorption in green  (NAMELIST)
   REAL(wp) ::   xlr        !: exposant for pigment absorption in red   (NAMELIST)
   REAL(wp) ::   xlg        !: exposant for pigment absorption in green (NAMELIST)
   REAL(wp) ::   rpig       !: chla/chla+phea ratio                     (NAMELIST)
                                                        
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   neln    !: number of levels in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   xze     !: euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xpar    !: par (photosynthetic available radiation)

!!----------------------------------------------------------------------
!! Sediment parameters                               
!! 
!! AXY (16/01/12): these parameters were originally part of the pre-
!!                 cursor model on which MEDUSA's code was grounded;
!!                 they do not relate to the sediment/benthos submodel
!!                 added as part of the ROAM project; they have only
!!                 been retained because they are distributed through
!!                 MEDUSA and require a proper clean-up to purge
!!----------------------------------------------------------------------
!!
   REAL(wp) ::   sedlam       !: time coefficient of POC remineralization in sediments
   REAL(wp) ::   sedlostpoc   !: ???
   REAL(wp) ::   areacot      !: ???
                                                        
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: dminl       !: fraction of sinking POC released in sediments
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmin3       !: fraction of sinking POC released at each level
                                                        
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sedpocb     !: mass of POC in sediments
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sedpocn     !: mass of POC in sediments
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sedpoca     !: mass of POC in sediments
                                                        
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fbodf       !: rapid sinking particles
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fbods       !: rapid sinking particles
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fbodn       !: rapid sinking particles
                                                        
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ffln        !: 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: fflf        !: 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ffls        !: 

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: cmask       !: ???

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_medusa_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_medusa_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      INTEGER ::   ierr(8)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
#if defined key_medusa
      !* depth-dependent damping coefficient
      ALLOCATE( remdmp(jpk,jp_medusa),                           STAT=ierr(1) )
# if defined key_roam
      !* 2D and 3D fields of carbonate system parameters
      ALLOCATE( f2_ccd_cal(jpi,jpj)  , f2_ccd_arg(jpi,jpj)  ,       &
         &      f3_pH(jpi,jpj,jpk)   , f3_h2co3(jpi,jpj,jpk),       &
         &      f3_hco3(jpi,jpj,jpk) , f3_co3(jpi,jpj,jpk)  ,       &
         &      f3_omcal(jpi,jpj,jpk), f3_omarg(jpi,jpj,jpk),    STAT=ierr(2) )
      !* 2D fields of organic and inorganic material sedimented on the seafloor
      ALLOCATE( zb_sed_n(jpi,jpj)    , zn_sed_n(jpi,jpj)    ,       &
         &      za_sed_n(jpi,jpj)    ,                              &
         &      zb_sed_fe(jpi,jpj)   , zn_sed_fe(jpi,jpj)   ,       &
         &      za_sed_fe(jpi,jpj)   ,                              &
         &      zb_sed_si(jpi,jpj)   , zn_sed_si(jpi,jpj)   ,       &
         &      za_sed_si(jpi,jpj)   ,                              &
         &      zb_sed_c(jpi,jpj)    , zn_sed_c(jpi,jpj)    ,       &
         &      za_sed_c(jpi,jpj)    ,                              &
         &      zb_sed_ca(jpi,jpj)   , zn_sed_ca(jpi,jpj)   ,       &
         &      za_sed_ca(jpi,jpj)   ,                           STAT=ierr(3) )
      !* 2D fields of temporally averaged properties for DMS calculations (AXY, 07/07/15)
      ALLOCATE( zb_dms_chn(jpi,jpj)  , zn_dms_chn(jpi,jpj)  ,       &
         &      za_dms_chn(jpi,jpj)  ,                              &
         &      zb_dms_chd(jpi,jpj)  , zn_dms_chd(jpi,jpj)  ,       &		  
         &      za_dms_chd(jpi,jpj)  ,                              &
         &      zb_dms_mld(jpi,jpj)  , zn_dms_mld(jpi,jpj)  ,       &		  
         &      za_dms_mld(jpi,jpj)  ,                              &
         &      zb_dms_qsr(jpi,jpj)  , zn_dms_qsr(jpi,jpj)  ,       &		  
         &      za_dms_qsr(jpi,jpj)  ,                              &
         &      zb_dms_din(jpi,jpj)  , zn_dms_din(jpi,jpj)  ,       &		  
         &      za_dms_din(jpi,jpj)  ,                           STAT=ierr(4) )
      !* 2D fields needing to be knows at first tstp for coupling with atm -
      !UKEMSi (Jpalm,14-06-2016) 
      ALLOCATE( zb_co2_flx(jpi,jpj)  , zn_co2_flx(jpi,jpj)  ,       &
         &      za_co2_flx(jpi,jpj)  ,                              &
         &      zb_dms_srf(jpi,jpj)  , zn_dms_srf(jpi,jpj)  ,       &           
         &      za_dms_srf(jpi,jpj)  , zn_chl_srf(jpi,jpj)  ,    STAT=ierr(5) )
# endif
      !* 2D fields of miscellaneous parameters
      ALLOCATE( ocal_ccd(jpi,jpj)    , dust(jpi,jpj)        ,       &
         &      zirondep(jpi,jpj)                           ,       &
         &      riv_n(jpi,jpj)                              ,       &
         &      riv_si(jpi,jpj)      , riv_c(jpi,jpj)       ,       &
         &      riv_alk(jpi,jpj)     , friver_dep(jpk,jpk)  ,    STAT=ierr(6) )
      !* 2D and 3D fields of light parameters
      ALLOCATE( neln(jpi,jpj)        , xze(jpi,jpj)         ,       &
         &      xpar(jpi,jpj,jpk)    ,                           STAT=ierr(7) )
      !* 2D and 3D fields of sediment-associated parameters
      ALLOCATE( dminl(jpi,jpj)       , dmin3(jpi,jpj,jpk)   ,       &
         &      sedpocb(jpi,jpj)     , sedpocn(jpi,jpj)     ,       &
         &      sedpoca(jpi,jpj)     , fbodn(jpi,jpj)       ,       &
         &      fbodf(jpi,jpj)       , fbods(jpi,jpj)       ,       &
         &      ffln(jpi,jpj,jpk)    , fflf(jpi,jpj,jpk)    ,       &
         &      ffls(jpi,jpj,jpk)    , cmask(jpi,jpj)       ,    STAT=ierr(8) ) 
#endif
      !
      sms_medusa_alloc = MAXVAL( ierr )
      !
      IF( sms_medusa_alloc /= 0 )   CALL ctl_warn('sms_medusa_alloc: failed to allocate arrays')
      !
   END FUNCTION sms_medusa_alloc

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     NO MEDUSA model 
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE sms_medusa
