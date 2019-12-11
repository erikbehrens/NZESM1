MODULE trcnam_medusa
   !!======================================================================
   !!                      ***  MODULE trcnam_medusa  ***
   !! TOP :   initialisation of some run parameters for MEDUSA bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!              -   !  2011-04  (A. Yool) updated for ROAM project
   !!              -   !  2013-05  (A. Yool) renamed (from trclsm) for v3.5
   !!              -   !  2015-11  (J. Palmieri) added iom_use for diags
   !!              -   !  2016-11  (A. Yool) updated diags for CMIP6
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'   :                                       MEDUSA model
   !!----------------------------------------------------------------------
   !! trc_nam_medusa      : MEDUSA model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_medusa      ! sms trends
   USE iom             ! I/O manager
   USE sbc_oce, ONLY: lk_oasis
   !!USE trc_nam_dia     ! JPALM 13-11-2015 -- if iom_use for diag

   !! AXY (04/02/14): necessary to find NaNs on HECTOR
!CEB
!$AGRIF_DO_NOT_TREAT
      USE, INTRINSIC :: ieee_arithmetic
!$AGRIF_END_DO_NOT_TREAT
!/CEB

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_medusa       ! called by trcnam.F90 module
   PUBLIC   trc_nam_iom_medusa   ! called by trcnam.F90 module

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
  !CEB!$AGRIF_DO_NOT_TREAT
   SUBROUTINE trc_nam_medusa
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_medusa  ***  
      !!
      !! ** Purpose :   read MEDUSA namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natbio, natopt, and natdbi ("key_trc_diabio")
      !!
      !! ekp: namelist nabio contains ALL parameters of the ecosystem 
      !!      point sourses and sinks PLUS sediment exchange
      !!      dia_bio - used by Lobster to output all point terms 
      !!                (sourses and sinks of bio)
      !!      dia_add - additional diagnostics for biology such as 
      !!                primary production (2d depth integrated field or 3d)
      !!----------------------------------------------------------------------
      !!
      INTEGER            :: ji,jj,jk
      REAL(wp)           :: fthk, fdep, fdep1
      REAL(wp)           :: q1, q2, q3
      !
      NAMELIST/natbio/ xxi,xaln,xald,jphy,xvpn,xvpd,          &
      &    xsin0,xnsi0,xuif,jliebig, jq10,                    &
      &    xthetam,xthetamd,xnln,xnld,xsld,xfln,xfld,         &
      &  xgmi,xgme,xkmi,xkme,xphi,xbetan,xbetac,xkc,          &
      &    xpmipn,xpmid,xpmepn,xpmepd,xpmezmi,xpmed,          &
      &  xmetapn,xmetapd,xmetazmi,xmetazme,                   &
      &  jmpn,xmpn,xkphn,jmpd,xmpd,xkphd,jmzmi,xmzmi,xkzmi,   &
      &    jmzme,xmzme,xkzme,jmd,jsfd,xmd,xmdc,               &
      &  xthetapn,xthetapd,xthetazmi,xthetazme,xthetad,       &
      &    xrfn,xrsn,vsed,xhr,                                &
      &  jiron,xfe_mass,xfe_sol,xfe_sed,xLgT,xk_FeL,xk_sc_Fe, &
      &  jexport,jfdfate,jrratio,jocalccd,xridg_r0,           &
      &    xfdfrac1,xfdfrac2,xfdfrac3,                        &
      &    xcaco3a,xcaco3b,xmassc,xmassca,xmasssi,xprotca,    &
      &    xprotsi,xfastc,xfastca,xfastsi,                    &
      &  jorgben,jinorgben,xsedn,xsedfe,xsedsi,xsedc,xsedca,  &
      &    xburial,                                           &
      &  jriver_n,jriver_si,jriver_c,jriver_alk,jriver_dep,   &
      &  xsdiss,                                              &
      &  sedlam,sedlostpoc,jpkb,jdms,jdms_input,jdms_model,   &
      &  scl_chl, chl_out, dmsmin, dmscut, dmsslp
#if defined key_roam
      NAMELIST/natroam/ xthetaphy,xthetazoo,xthetanit,        &
      &    xthetarem,xo2min 
#endif
      NAMELIST/natopt/xkg0,xkr0,xkgp,xkrp,xlg,xlr,rpig
      INTEGER :: jl, jn
      INTEGER :: ios                 ! Local integer output status for namelist read
      CHARACTER(LEN=32)   ::   clname
      !!
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      clname = 'namelist_medusa'
      IF(lwp) WRITE(numout,*) ' trc_nam_medusa: read MEDUSA namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
# if defined key_debug_medusa
      CALL flush(numout)
# endif


      CALL ctl_opn( numnatp_ref, TRIM( clname )//'_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfg, TRIM( clname )//'_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonp     , 'output.namelist.pis' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

# if defined key_debug_medusa
      CALL flush(numout)
      IF (lwp) write (numout,*) '------------------------------'
      IF (lwp) write (numout,*) 'Jpalm - debug'
      IF (lwp) write (numout,*) 'open namelist_medusa -- OK'
      IF (lwp) write (numout,*) 'Now, read namilists inside :'
      IF (lwp) write (numout,*) ' '
# endif
      !
# if defined key_debug_medusa
      CALL flush(numout)
# endif

      ! 1.4 namelist natbio : biological parameters
      ! -------------------------------------------
      !! Note: the default values below will all be overwritten by the
      !!       input in the namelist natbio.
      
      xxi         = 0.
      xaln        = 0.
      xald        = 0.
      jphy        = 0
      xvpn        = 0.
      xvpd        = 0.
      xthetam     = 0.
      xthetamd    = 0.
!!
      xsin0       = 0.
      xnsi0       = 0.
      xuif        = 0.
!!
      jliebig     = 0
      jq10        = 0.
      xnln        = 0.
      xnld        = 0.
      xsld        = 0.
      xfln        = 0.
      xfld        = 0.
!!
      xgmi        = 0.
      xgme        = 0.
      xkmi        = 0.
      xkme        = 0.
      xphi	  = 0.
      xbetan      = 0.
      xbetac      = 0.
      xkc         = 0.
      xpmipn      = 0.
      xpmid       = 0.
      xpmepn      = 0.
      xpmepd      = 0.
      xpmezmi     = 0.
      xpmed       = 0.
!!
      xmetapn     = 0.
      xmetapd     = 0.
      xmetazmi    = 0.
      xmetazme    = 0.
!!
      jmpn        = 0
      xmpn        = 0.
      xkphn       = 0.
      jmpd        = 0
      xmpd        = 0.
      xkphd       = 0.
      jmzmi       = 0
      xmzmi       = 0.
      xkzmi       = 0.
      jmzme       = 0
      xmzme       = 0.
      xkzme       = 0.
!!
      jmd         = 0
      jsfd        = 0
      xmd         = 0.
      xmdc        = 0.
!!
      xthetapn    = 0.
      xthetapd    = 0.
      xthetazmi   = 0.
      xthetazme   = 0.
      xthetad     = 0.
      xrfn        = 0.
      xrsn        = 0.  !: (NOT USED HERE; RETAINED FOR LOBSTER)
!!
      jiron       = 0
      xfe_mass    = 0.
      xfe_sol     = 0.
      xfe_sed     = 0.
      xLgT        = 0.
      xk_FeL	  = 0.
      xk_sc_Fe    = 0.
!!
      jexport     = 0
      jfdfate     = 0
      jrratio     = 0
      jocalccd    = 0
      xridg_r0    = 0.
      xfdfrac1	  = 0.
      xfdfrac2	  = 0.
      xfdfrac3	  = 0.
      xcaco3a	  = 0.
      xcaco3b	  = 0.
      xmassc	  = 0.
      xmassca	  = 0.
      xmasssi	  = 0.
      xprotca	  = 0.
      xprotsi	  = 0.
      xfastc	  = 0.
      xfastca	  = 0.
      xfastsi	  = 0.
!!
      jorgben     = 0
      jinorgben   = 0
      xsedn       = 0.
      xsedfe      = 0.
      xsedsi      = 0.
      xsedc       = 0.
      xsedca      = 0.
      xburial     = 0.
!!
      jriver_n    = 0
      jriver_si   = 0
      jriver_c    = 0
      jriver_alk  = 0
      jriver_dep  = 1
!!
      xsdiss	  = 0.
!!
      vsed        = 0.
      xhr         = 0.
!!
      sedlam	  = 0.
      sedlostpoc  = 0.
      jpkb	  = 0.
      jdms        = 0
      jdms_input  = 0
      jdms_model  = 0
      scl_chl     = 1.
      chl_out     = 1
      dmsmin      = 2.29 !! Anderson DMS default
      dmscut      = 1.72 !! Anderson DMS default
      dmsslp      = 8.24 !! Anderson DMS default
            
      !REWIND(numnatm)
      !READ(numnatm,natbio)
         ! Namelist natbio
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist natbio in reference namelist : MEDUSA diagnostics
         READ  ( numnatp_ref, natbio, IOSTAT = ios, ERR = 903)
903      IF( ios /= 0 ) CALL ctl_nam ( ios , 'natbio in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist natbio in configuration namelist : MEDUSA diagnostics
         READ  ( numnatp_cfg, natbio, IOSTAT = ios, ERR = 904 )
904      IF( ios /= 0 ) CALL ctl_nam ( ios , 'natbio in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, natbio )

!! Primary production and chl related quantities
!!       xxi         :  conversion factor from gC to mmolN 
!!       xaln        :  Chl-a specific initial slope of P-I curve for non-diatoms 
!!       xald        :  Chl-a specific initial slope of P-I curve for diatoms
!!       jphy        :  phytoplankton T-dependent growth switch
!!       xvpn        :  maximum growth rate for non-diatoms
!!       xvpd        :  maximum growth rate for diatoms
!!       xthetam     :  maximum Chl to C ratio for non-diatoms      
!!       xthetamd    :  maximum Chl to C ratio for diatoms      
!!
!! Diatom silicon parameters
!!       xsin0       :  minimum diatom Si:N ratio
!!       xnsi0       :  minimum diatom N:Si ratio
!!       xuif        :  hypothetical growth ratio at infinite Si:N ratio
!!
!! Nutrient limitation
!!       jliebig     :  Liebig nutrient uptake switch
!!       xnln        :  half-sat constant for DIN uptake by non-diatoms 
!!       xnld        :  half-sat constant for DIN uptake by diatoms 
!!       xsl         :  half-sat constant for Si uptake by diatoms 
!!       xfld        :  half-sat constant for Fe uptake by diatoms  
!!       xfln        :  half-sat constant for Fe uptake by non-datoms 
!!
!! Grazing
!!       xgmi        :  microzoo maximum growth rate 
!!       xgme        :  mesozoo maximum growth rate 
!!       xkmi        :  microzoo grazing half-sat parameter
!!       xkme        :  mesozoo grazing half-sat parameter
!!       xphi        :  micro/mesozoo grazing inefficiency
!!       xbetan      :  micro/mesozoo N assimilation efficiency
!!       xbetac      :  micro/mesozoo C assimilation efficiency
!!       xkc         :  micro/mesozoo net C growth efficiency
!!       xpmipn      :  grazing preference of microzoo for non-diatoms
!!       xpmid       :  grazing preference of microzoo for diatoms
!!       xpmepn      :  grazing preference of mesozoo for non-diatoms 
!!       xpmepd      :  grazing preference of mesozoo for diatoms
!!       xpmezmi     :  grazing preference of mesozoo for microzoo
!!       xpmed       :  grazing preference of mesozoo for detritus
!!
!! Metabolic losses
!!       xmetapn     :  non-diatom metabolic loss rate
!!       xmetapd     :  diatom     metabolic loss rate
!!       xmetazmi    :  microzoo   metabolic loss rate
!!       xmetazme    :  mesozoo    metabolic loss rate
!!
!! Mortality/Remineralisation
!!       jmpn        :  non-diatom mortality functional form
!!       xmpn        :  non-diatom mortality rate
!!       xkphn       :  non-diatom mortality half-sat constant
!!       jmpd        :  diatom     mortality functional form
!!       xmpd        :  diatom mortality rate
!!       xkphd       :  diatom mortality half-sat constant
!!       jmzmi       :  microzoo   mortality functional form
!!       xmzmi       :  microzoo mortality rate
!!       xkzmi       :  microzoo mortality half-sat constant
!!       jmzme       :  mesozoo    mortality functional form
!!       xmzme       :  mesozoo mortality rate
!!       xkzme       :  mesozoo mortality half-sat constant
!!
!! Remineralisation
!!       jmd         :  detritus T-dependent remineralisation switch
!!       jsfd        :  accelerate seafloor detritus remin. switch
!!       xmd         :  detrital nitrogen remineralisation rate 
!!       xmdc        :  detrital carbon remineralisation rate 
!!
!! Stochiometric ratios
!!       xthetapn    :  non-diatom C:N ratio
!!       xthetapd    :  diatom C:N ratio
!!       xthetazmi   :  microzoo C:N ratio
!!       xthetazme   :  mesozoo C:N ratio
!!       xthetad     :  detritus C:N ratio
!!       xrfn        :  phytoplankton Fe:N ratio
!!	 xrsn        :  diatom Si:N ratio (*NOT* used)
!!
!! Iron parameters
!!       jiron       :  iron scavenging submodel switch
!!       xfe_mass    :  iron atomic mass
!!	 xfe_sol     :  aeolian iron solubility
!!	 xfe_sed     :  sediment iron input
!!	 xLgT	     :  total ligand concentration (umol/m3)
!!	 xk_FeL	     :  dissociation constant for (Fe + L)
!!	 xk_sc_Fe    :  scavenging rate of "free" iron
!!	 
!! Fast-sinking detritus parameters
!!       jexport     :  fast detritus remineralisation switch
!!       jfdfate     :  fate of fast detritus at seafloor switch
!!       jrratio     :  rain ratio switch
!!       jocalccd    :  CCD switch
!!       xridg_r0    :  Ridgwell rain ratio coefficient
!!       xfdfrac1    :  fast-sinking fraction of diatom nat. mort. losses
!!       xfdfrac2    :  fast-sinking fraction of meszooplankton mort. losses
!!       xfdfrac3    :  fast-sinking fraction of diatom silicon grazing losses
!!       xcaco3a     :  polar (high latitude) CaCO3 fraction
!!       xcaco3b     :  equatorial (low latitude) CaCO3 fraction
!!       xmassc      :  organic C mass:mole ratio, C106 H175 O40 N16 P1
!!       xmassca     :  calcium carbonate mass:mole ratio, CaCO3
!!       xmasssi     :  biogenic silicon mass:mole ratio, (H2SiO3)n
!!       xprotca     :  calcium carbonate protection ratio
!!       xprotsi     :  biogenic silicon protection ratio
!!       xfastc      :  organic C remineralisation length scale
!!       xfastca     :  calcium carbonate dissolution length scale
!!       xfastsi     :  biogenic silicon dissolution length scale
!!
!! Benthic 
!!       jorgben     :  does   organic detritus go to the benthos?
!!       jinorgben   :  does inorganic detritus go to the benthos?
!!       xsedn       :  organic   nitrogen sediment remineralisation rate 
!!       xsedfe      :  organic   iron     sediment remineralisation rate 
!!       xsedsi      :  inorganic silicon  sediment dissolution      rate 
!!       xsedc       :  organic   carbon   sediment remineralisation rate 
!!       xsedca      :  inorganic carbon   sediment dissolution      rate 
!!       xburial     :  burial rate of seafloor detritus
!!
!! Riverine inputs
!!       jriver_n    :  riverine N          input?
!!       jriver_si   :  riverine Si         input?
!!       jriver_c    :  riverine C          input?
!!       jriver_alk  :  riverine alkalinity input?
!!       jriver_dep  :  depth of riverine   input?
!!
!! Miscellaneous
!!       xsdiss      :  diatom frustule dissolution rate
!!
!! Gravitational sinking      
!!       vsed        :  detritus gravitational sinking rate 
!!       xhr         :  coeff for Martin's remineralisation profile
!!
!! Additional parameters
!!       sedlam      :  time coeff of POC in sediments
!!      sedlostpoc   :  sediment geol loss for POC
!!       jpkb        :  vertical layer for diagnostic of the vertical flux 
!!                      NOTE that in LOBSTER it is a first vertical layers where 
!!                      biology is active  
!!
!! UKESM1 - new diagnostics  !! Jpalm
!!       jdms        :  include dms diagnostics
!!	 jdms_input  :  use instant (0) or diel-avg (1) inputs
!!       jdms_model  :  choice of DMS model passed to atmosphere
!!                      1 = ANDR, 2 = SIMO, 3 = ARAN, 4 = HALL, 5 = ANDM
!!       dmsmin      : DMS minimum value for DMS Anderson (ANDR) sheme ONLY
!!       dmscut      : DMS cutoff value for DMS Anderson (ANDR) sheme ONLY
!!       dmsslp      : DMS slope value for DMS Anderson (ANDR) sheme ONLY
!! UKESM1 - exported Chl to UM
!!       scl_chl     : scaling factor to tune the chl field sent to the UM
!!       chl_out     : select the chl field to send at the UM:
!!                     1- Surf Chl ; 2- MLD Chl 

      IF(lwp) THEN
!!
!! AXY (08/11/13): compilation key notification
         WRITE(numout,*) '=== Compilation keys'
#if defined key_roam
         WRITE(numout,*)     &
         &   ' key_roam                                                               = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_roam                                                               = INACTIVE'
#endif        
#if defined key_axy_carbchem
         WRITE(numout,*)     &
         &   ' key_axy_carbchem                                                       = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_axy_carbchem                                                       = INACTIVE'
#endif        
#if defined key_mocsy
         WRITE(numout,*)     &
         &   ' key_mocsy                                                              = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_mocsy                                                              = INACTIVE'
#endif        
#if defined key_avgqsr_medusa
         WRITE(numout,*)     &
         &   ' key_avgqsr_medusa                                                      = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_avgqsr_medusa                                                      = INACTIVE'
#endif        
#if defined key_bs_axy_zforce
         WRITE(numout,*)     &
         &   ' key_bs_axy_zforce                                                      = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_bs_axy_zforce                                                      = INACTIVE'
#endif        
#if defined key_bs_axy_yrlen
         WRITE(numout,*)     &
         &   ' key_bs_axy_yrlen                                                       = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_bs_axy_yrlen                                                       = INACTIVE'
#endif        
#if defined key_deep_fe_fix
         WRITE(numout,*)     &
         &   ' key_deep_fe_fix                                                        = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_deep_fe_fix                                                        = INACTIVE'
#endif        
#if defined key_axy_nancheck
         WRITE(numout,*)     &
         &   ' key_axy_nancheck                                                       = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_axy_nancheck                                                       = INACTIVE'
#endif        
# if defined key_axy_pi_co2
         WRITE(numout,*)     &
         &   ' key_axy_pi_co2                                                         = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_axy_pi_co2                                                         = INACTIVE'
# endif
# if defined key_debug_medusa
         WRITE(numout,*)     &
         &   ' key_debug_medusa                                                       = ACTIVE'
#else
         WRITE(numout,*)     &
         &   ' key_debug_medusa                                                       = INACTIVE'
# endif
         WRITE(numout,*) ' '

         WRITE(numout,*) 'natbio'
         WRITE(numout,*) ' '
!!
!! Primary production and chl related quantities
         WRITE(numout,*) '=== Primary production'
         WRITE(numout,*)     &
         &   ' conversion factor from gC to mmolN,                        xxi         =', xxi
         WRITE(numout,*)     &
         &   ' Chl-a specific initial slope of P-I curve for non-diatoms, xaln        = ', xaln
         WRITE(numout,*)     &
         &   ' Chl-a specific initial slope of P-I curve for diatoms,     xald        = ', xald
         if (jphy.eq.1) then
            WRITE(numout,*) &
            &   ' phytoplankton growth is *temperature-dependent*            jphy        = ', jphy
         elseif (jphy.eq.2) then
            WRITE(numout,*) &
            &   ' phytoplankton growth is *temperature-dependent(Q10)*       jphy        = ', jphy
         elseif (jphy.eq.0) then
            WRITE(numout,*) &
            &   ' phytoplankton growth is *temperature-independent*          jphy        = ', jphy
         endif
         WRITE(numout,*)     &
         &   ' maximum growth rate for non-diatoms,                       xvpn        = ', xvpn
         WRITE(numout,*)     &
         &   ' maximum growth rate for diatoms,                           xvpn        = ', xvpd
         WRITE(numout,*)     &
         &   ' maximum Chl to C ratio for non-diatoms,                    xthetam     = ', xthetam
         WRITE(numout,*)     &
         &   ' maximum Chl to C ratio for diatoms,                        xthetamd    = ', xthetamd
         WRITE(numout,*)     &
         &   ' specific Q10 value (jphy==2),                                  jq10    = ', jq10
!!
!! Diatom silicon parameters
         WRITE(numout,*) '=== Diatom silicon parameters'
         WRITE(numout,*)     &
         &   ' minimum diatom Si:N ratio,                                 xsin0       = ', xsin0
         WRITE(numout,*)     &
         &   ' minimum diatom N:Si ratio,                                 xnsi0       = ', xnsi0
         WRITE(numout,*)     &
         &   ' hypothetical growth ratio at infinite Si:N ratio,          xuif        = ', xuif
!!
!! Nutrient limitation
         WRITE(numout,*) '=== Nutrient limitation'
         if (jliebig.eq.1) then
            WRITE(numout,*) &
            &   ' nutrient uptake is a Liebig Law (= most limiting) function jliebig     = ', jliebig
         elseif (jliebig.eq.0) then
            WRITE(numout,*) &
            &   ' nutrient uptake is a multiplicative function               jliebig     = ', jliebig
         endif
         WRITE(numout,*)     &
         &   ' half-sat constant for DIN uptake by non-diatoms,           xnln        = ', xnln
         WRITE(numout,*)     &
         &   ' half-sat constant for DIN uptake by diatoms,               xnld        = ', xnld
         WRITE(numout,*)     &
         &   ' half-sat constant for Si uptake by diatoms,                xsld        = ', xsld
         WRITE(numout,*)     &
         &   ' half-sat constant for Fe uptake by non-diatoms,            xfln        = ', xfln
         WRITE(numout,*)     &
         &   ' half-sat constant for Fe uptake by diatoms,                xfld        = ', xfld
!!
!! Grazing
         WRITE(numout,*) '=== Zooplankton grazing'
         WRITE(numout,*)     &
         &   ' microzoo maximum growth rate,                              xgmi        = ', xgmi
         WRITE(numout,*)     &
         &   ' mesozoo maximum growth rate,                               xgme        = ', xgme
         WRITE(numout,*)     &
         &   ' microzoo grazing half-sat parameter,                       xkmi        = ', xkmi
         WRITE(numout,*)     &
         &   ' mesozoo grazing half-sat parameter,                        xkme        = ', xkme
         WRITE(numout,*)     &
         &   ' micro/mesozoo grazing inefficiency,                        xphi        = ', xphi
         WRITE(numout,*)     &
         &   ' micro/mesozoo N assimilation efficiency,                   xbetan      = ', xbetan
         WRITE(numout,*)     &
         &   ' micro/mesozoo C assimilation efficiency,                   xbetac      = ', xbetan
         WRITE(numout,*)     &
         &   ' micro/mesozoo net C growth efficiency,                     xkc         = ', xkc
         WRITE(numout,*)     &
         &   ' grazing preference of microzoo for non-diatoms,            xpmipn      = ', xpmipn
         WRITE(numout,*)     &
         &   ' grazing preference of microzoo for detritus,               xpmid       = ', xpmid
         WRITE(numout,*)     &
         &   ' grazing preference of mesozoo for non-diatoms,             xpmepn      = ', xpmepn
         WRITE(numout,*)     &
         &   ' grazing preference of mesozoo for diatoms,                 xpmepd      = ', xpmepd
         WRITE(numout,*)     &
         &   ' grazing preference of mesozoo for microzoo,                xpmezmi     = ', xpmezmi
         WRITE(numout,*)     &
         &   ' grazing preference of mesozoo for detritus,                xpmed       = ', xpmed
!!
!! Metabolic losses
         WRITE(numout,*) '=== Metabolic losses'
         WRITE(numout,*)     &
         &   ' non-diatom metabolic loss rate,                            xmetapn     = ', xmetapn
         WRITE(numout,*)     &
         &   ' diatom     metabolic loss rate,                            xmetapd     = ', xmetapd
         WRITE(numout,*)     &
         &   ' microzoo   metabolic loss rate,                            xmetazmi    = ', xmetazmi
         WRITE(numout,*)     &
         &   ' mesozoo    metabolic loss rate,                            xmetazme    = ', xmetazme
!!
!! Mortality losses
         WRITE(numout,*) '=== Mortality losses'
         if (jmpn.eq.1) then
            WRITE(numout,*)     &
            &   ' non-diatom mortality functional form,            LINEAR    jmpn        = ', jmpn
         elseif (jmpn.eq.2) then
            WRITE(numout,*)     &
            &   ' non-diatom mortality functional form,         QUADRATIC    jmpn        = ', jmpn
         elseif (jmpn.eq.3) then
            WRITE(numout,*)     &
            &   ' non-diatom mortality functional form,        HYPERBOLIC    jmpn        = ', jmpn
         elseif (jmpn.eq.4) then
            WRITE(numout,*)     &
            &   ' non-diatom mortality functional form,           SIGMOID    jmpn        = ', jmpn
         endif
         WRITE(numout,*)     &
         &   ' non-diatom mortality rate,                                 xmpn        = ', xmpn
         WRITE(numout,*)     &
         &   ' non-diatom mortality half-sat constant                     xkphn       = ', xkphn
         if (jmpd.eq.1) then
            WRITE(numout,*)     &
            &   ' diatom mortality functional form,                LINEAR    jmpd        = ', jmpd
         elseif (jmpd.eq.2) then
            WRITE(numout,*)     &
            &   ' diatom mortality functional form,             QUADRATIC    jmpd        = ', jmpd
         elseif (jmpd.eq.3) then
            WRITE(numout,*)     &
            &   ' diatom mortality functional form,            HYPERBOLIC    jmpd        = ', jmpd
         elseif (jmpd.eq.4) then
            WRITE(numout,*)     &
            &   ' diatom mortality functional form,               SIGMOID    jmpd        = ', jmpd
         endif
         WRITE(numout,*)     &
         &   ' diatom mortality rate,                                     xmpd        = ', xmpd
         WRITE(numout,*)     &
         &   ' diatom mortality half-sat constant                         xkphd       = ', xkphd
         if (jmzmi.eq.1) then
            WRITE(numout,*)     &
            &   ' microzoo mortality functional form,              LINEAR    jmzmi       = ', jmzmi
         elseif (jmzmi.eq.2) then
            WRITE(numout,*)     &
            &   ' microzoo mortality functional form,           QUADRATIC    jmzmi       = ', jmzmi
         elseif (jmzmi.eq.3) then
            WRITE(numout,*)     &
            &   ' microzoo mortality functional form,          HYPERBOLIC    jmzmi       = ', jmzmi
         elseif (jmzmi.eq.4) then
            WRITE(numout,*)     &
            &   ' microzoo mortality functional form,             SIGMOID    jmzmi       = ', jmzmi
         endif
         WRITE(numout,*)     &
         &   ' microzoo mortality rate,                                   xmzmi       = ', xmzmi
         WRITE(numout,*)     &
         &   ' mesozoo mortality half-sat constant,                       xkzmi       = ', xkzmi
         if (jmzme.eq.1) then
            WRITE(numout,*)     &
            &   ' mesozoo mortality functional form,               LINEAR    jmzme       = ', jmzme
         elseif (jmzme.eq.2) then
            WRITE(numout,*)     &
            &   ' mesozoo mortality functional form,            QUADRATIC    jmzme       = ', jmzme
         elseif (jmzme.eq.3) then
            WRITE(numout,*)     &
            &   ' mesozoo mortality functional form,           HYPERBOLIC    jmzme       = ', jmzme
         elseif (jmzme.eq.4) then
            WRITE(numout,*)     &
            &   ' mesozoo mortality functional form,              SIGMOID    jmzme       = ', jmzme
         endif
         WRITE(numout,*)     &
         &   ' mesozoo mortality rate,                                    xmzme       = ', xmzme
         WRITE(numout,*)     &
         &   ' mesozoo mortality half-sat constant,                       xkzme       = ', xkzme
!!
!! Remineralisation
         WRITE(numout,*) '=== Remineralisation'
         if (jmd.eq.1) then
            WRITE(numout,*) &
            &   ' detritus remineralisation is *temperature-dependent*       jmd         = ', jmd
         elseif (jmd.eq.2) then
            WRITE(numout,*) &
            &   ' detritus remineralisation is *temperature-dependent(Q10)*  jmd         = ', jmd
         elseif (jmd.eq.0) then
            WRITE(numout,*) &
            &   ' detritus remineralisation is *temperature-independent*     jmd         = ', jmd
         endif
         if (jsfd.eq.1) then
            WRITE(numout,*) &
            &   ' detritus seafloor remineralisation is *accelerated*        jsfd        = ', jsfd
         else
            WRITE(numout,*) &
            &   ' detritus seafloor remineralisation occurs at same rate     jsfd        = ', jsfd
         endif
         WRITE(numout,*)     &
         &   ' detrital nitrogen remineralisation rate,                   xmd         = ', xmd
         WRITE(numout,*)     &
         &   ' detrital carbon remineralisation rate,                     xmdc        = ', xmdc
!!
!! Stochiometric ratios
         WRITE(numout,*) '=== Stoichiometric ratios'
         WRITE(numout,*)     &
         &   ' non-diatom C:N ratio,                                      xthetapn    = ', xthetapn
         WRITE(numout,*)     &
         &   ' diatom C:N ratio,                                          xthetapd    = ', xthetapd
         WRITE(numout,*)     &
         &   ' microzoo C:N ratio,                                        xthetazmi   = ', xthetazmi
         WRITE(numout,*)     &
         &   ' mesozoo C:N ratio,                                         xthetazme   = ', xthetazme
         WRITE(numout,*)     &
         &   ' detritus C:N ratio,                                        xthetad     = ', xthetad
         WRITE(numout,*)     &
         &   ' phytoplankton Fe:N ratio,                                  xrfn        = ', xrfn
         WRITE(numout,*)     &
         &   ' diatom Si:N ratio,                                         xrsn        = ', xrsn
!!	  
!! Iron parameters
         WRITE(numout,*) '=== Iron parameters'
         if (jiron.eq.1) then
            WRITE(numout,*)     &
            &   ' Dutkiewicz et al. (2005) iron scavenging                   jiron       = ', jiron
         elseif (jiron.eq.2) then
            WRITE(numout,*)     &
            &   ' Moore et al. (2004) iron scavenging                        jiron       = ', jiron
         elseif (jiron.eq.3) then
            WRITE(numout,*)     &
            &   ' Moore et al. (2008) iron scavenging                        jiron       = ', jiron
         elseif (jiron.eq.4) then
            WRITE(numout,*)     &
            &   ' Galbraith et al. (2010) iron scavenging                    jiron       = ', jiron
         else
            WRITE(numout,*)     &
            &   ' There is **no** iron scavenging                            jiron       = ', jiron
         endif
         WRITE(numout,*)     &
         &   ' iron atomic mass,                                          xfe_mass    = ', xfe_mass
         WRITE(numout,*)     &
         &   ' aeolian iron solubility,                                   xfe_sol     = ', xfe_sol
         WRITE(numout,*)     &
         &   ' sediment iron input,                                       xfe_sed     = ', xfe_sed
         WRITE(numout,*)     &
         &   ' total ligand concentration (umol/m3),                      xLgT        = ', xLgT
         WRITE(numout,*)     &
         &   ' dissociation constant for (Fe + L),                        xk_FeL      = ', xk_FeL
         WRITE(numout,*)     &
         &   ' scavenging rate for free iron,                             xk_sc_Fe    = ', xk_sc_Fe
!!
!! Fast-sinking detritus parameters
         WRITE(numout,*) '=== Fast-sinking detritus'
         if (jexport.eq.1) then
            WRITE(numout,*) &
            &   ' fast-detritus remin. uses Dunne et al. (2007; ballast)     jexport     = ', jexport
         elseif (jexport.eq.2) then
            WRITE(numout,*) &
            &   ' fast-detritus remin. uses Martin et al. (1987)             jexport     = ', jexport
         elseif (jexport.eq.2) then
            WRITE(numout,*) &
            &   ' fast-detritus remin. uses Henson et al. (2011)             jexport     = ', jexport
         endif
         if (jfdfate.eq.1) then
            WRITE(numout,*) &
            &   ' fast-detritus reaching seafloor becomes slow-detritus      jfdfate     = ', jfdfate
         elseif (jfdfate.eq.0) then
            WRITE(numout,*) &
            &   ' fast-detritus reaching seafloor instantly remineralised    jfdfate     = ', jfdfate
         endif
#if defined key_roam
         if (jrratio.eq.0) then
            WRITE(numout,*) &
            &   ' Dunne et al. (2005) rain ratio submodel                    jrratio     = ', jrratio
         elseif (jrratio.eq.1) then
            WRITE(numout,*) &
            &   ' Ridgwell et al. (2007) rain ratio submodel (surface omega) jrratio     = ', jrratio
         elseif (jrratio.eq.2) then
            WRITE(numout,*) &
            &   ' Ridgwell et al. (2007) rain ratio submodel (3D omega)      jrratio     = ', jrratio
         endif
#else          
         jrratio = 0
         WRITE(numout,*) &
         &   ' Dunne et al. (2005) rain ratio submodel                    jrratio     = ', jrratio
#endif          
#if defined key_roam
         if (jocalccd.eq.0) then
            WRITE(numout,*) &
            &   ' Default, fixed CCD used                                    jocalccd    = ', jocalccd
         elseif (jocalccd.eq.1) then
            WRITE(numout,*) &
            &   ' Calculated, dynamic CCD used                               jocalccd    = ', jocalccd
         endif
#else          
         jocalccd = 0
         WRITE(numout,*) &
         &   ' Default, fixed CCD used                                    jocalccd    = ', jocalccd
#endif
         WRITE(numout,*)     &
         &   ' Ridgwell rain ratio coefficient,                           xridg_r0    = ', xridg_r0
         WRITE(numout,*)     &
         &   ' fast-sinking fraction of diatom nat. mort. losses,         xfdfrac1    = ', xfdfrac1
         WRITE(numout,*)     &
         &   ' fast-sinking fraction of mesozooplankton mort. losses,     xfdfrac2    = ', xfdfrac2
         WRITE(numout,*)     &
         &   ' fast-sinking fraction of diatom silicon grazing losses,    xfdfrac3    = ', xfdfrac3
         WRITE(numout,*)     &
         &   ' polar (high latitude) CaCO3 fraction,                      xcaco3a     = ', xcaco3a
         WRITE(numout,*)     &
         &   ' equatorial (low latitude) CaCO3 fraction,                  xcaco3b     = ', xcaco3b
         WRITE(numout,*)     &
         &   ' organic C mass:mole ratio, C106 H175 O40 N16 P1,           xmassc      = ', xmassc
         WRITE(numout,*)     &
         &   ' calcium carbonate mass:mole ratio, CaCO3,                  xmassca     = ', xmassca
         WRITE(numout,*)     &
         &   ' biogenic silicon mass:mole ratio, (H2SiO3)n,               xmasssi     = ', xmasssi
         WRITE(numout,*)     &
         &   ' calcium carbonate protection ratio,                        xprotca     = ', xprotca
         WRITE(numout,*)     &
         &   ' biogenic silicon protection ratio,                         xprotsi     = ', xprotsi
         WRITE(numout,*)     &
         &   ' organic C remineralisation length scale,                   xfastc      = ', xfastc
         WRITE(numout,*)     &
         &   ' calcium carbonate dissolution length scale,                xfastca     = ', xfastca
         WRITE(numout,*)     &
         &   ' biogenic silicon dissolution length scale,                 xfastsi     = ', xfastsi
!!
!! Benthos parameters
         WRITE(numout,*) '=== Benthos parameters'
         WRITE(numout,*)     &
         &   ' does   organic detritus go to the benthos?,                jorgben     = ', jorgben
         WRITE(numout,*)     &
         &   ' does inorganic detritus go to the benthos?,                jinorgben   = ', jinorgben
!!
!! Some checks on parameters related to benthos parameters
         if (jorgben.eq.1 .and. jsfd.eq.1) then
            !! slow detritus going to benthos at an accelerated rate
            WRITE(numout,*) '  === WARNING! ==='
            WRITE(numout,*) '  jsfd *and* jorgben are active - please check that you wish this'
            WRITE(numout,*) '  === WARNING! ==='
         endif
         if (jorgben.eq.1 .and. jfdfate.eq.1) then
            !! fast detritus going to benthos but via slow detritus
            WRITE(numout,*) '  === WARNING! ==='
            WRITE(numout,*) '  jfdfate *and* jorgben are active - please check that you wish this'
            WRITE(numout,*) '  === WARNING! ==='
         endif
         if (jorgben.eq.0 .and. jinorgben.eq.1) then
            !! inorganic fast detritus going to benthos but organic fast detritus is not
            WRITE(numout,*) '  === WARNING! ==='
            WRITE(numout,*) '  jinorgben is active but jorgben is not - please check that you wish this'
            WRITE(numout,*) '  === WARNING! ==='
         endif
         WRITE(numout,*)     &
         &   ' organic   nitrogen sediment remineralisation rate,         xsedn       = ', xsedn
         WRITE(numout,*)     &
         &   ' organic   iron     sediment remineralisation rate,         xsedfe      = ', xsedfe
         WRITE(numout,*)     &
         &   ' inorganic silicon  sediment remineralisation rate,         xsedsi      = ', xsedsi
         WRITE(numout,*)     &
         &   ' organic   carbon   sediment remineralisation rate,         xsedc       = ', xsedc
         WRITE(numout,*)     &
         &   ' inorganic carbon   sediment remineralisation rate,         xsedca      = ', xsedca
         WRITE(numout,*)     &
         &   ' burial rate of seafloor detritus,                          xburial     = ', xburial
!!
!! Riverine inputs
         WRITE(numout,*) '=== Riverine inputs'
         if (jriver_n.eq.0) then
            WRITE(numout,*)     &
            &   ' *no* riverine N input,                                     jriver_n    = ', jriver_n
         elseif (jriver_n.eq.1) then
            WRITE(numout,*)     &
            &   ' riverine N concentrations specified,                       jriver_n    = ', jriver_n
         elseif (jriver_n.eq.2) then
            WRITE(numout,*)     &
            &   ' riverine N inputs specified,                               jriver_n    = ', jriver_n
         endif
         if (jriver_si.eq.0) then
            WRITE(numout,*)     &
            &   ' *no* riverine Si input,                                    jriver_si   = ', jriver_si
         elseif (jriver_si.eq.1) then
            WRITE(numout,*)     &
            &   ' riverine Si concentrations specified,                      jriver_si   = ', jriver_si
         elseif (jriver_si.eq.2) then
            WRITE(numout,*)     &
            &   ' riverine Si inputs specified,                              jriver_si   = ', jriver_si
         endif
         if (jriver_c.eq.0) then
            WRITE(numout,*)     &
            &   ' *no* riverine C input,                                     jriver_c    = ', jriver_c
         elseif (jriver_c.eq.1) then
            WRITE(numout,*)     &
            &   ' riverine C concentrations specified,                       jriver_c    = ', jriver_c
         elseif (jriver_c.eq.2) then
            WRITE(numout,*)     &
            &   ' riverine C inputs specified,                               jriver_c    = ', jriver_c
         endif
         if (jriver_alk.eq.0) then
            WRITE(numout,*)     &
            &   ' *no* riverine alkalinity input,                            jriver_alk  = ', jriver_alk
         elseif (jriver_alk.eq.1) then
            WRITE(numout,*)     &
            &   ' riverine alkalinity concentrations specified,              jriver_alk  = ', jriver_alk
         elseif (jriver_alk.eq.2) then
            WRITE(numout,*)     &
            &   ' riverine alkalinity inputs specified,                      jriver_alk  = ', jriver_alk
         endif
         !! AXY (19/07/12): prevent (gross) stupidity on part of user
         if (jriver_dep.lt.1.or.jriver_dep.ge.jpk) then
            jriver_dep = 1
         endif
         WRITE(numout,*)     &
         &   ' riverine input applied to down to depth k = ...            jriver_dep  = ', jriver_dep
!!
!! Miscellaneous
         WRITE(numout,*) '=== Miscellaneous'
         WRITE(numout,*)     &
         &   ' diatom frustule dissolution rate,                          xsdiss      = ', xsdiss
!!
!! Gravitational sinking      
         WRITE(numout,*) '=== Gravitational sinking'
         WRITE(numout,*)     &
         &   ' detritus gravitational sinking rate,                       vsed        = ', vsed
         WRITE(numout,*)     & 
         &   ' coefficient for Martin et al. (1987) remineralisation,     xhr         = ', xhr
!!
!! Non-Medusa parameters
         WRITE(numout,*) '=== Non-Medusa parameters'
         WRITE(numout,*)     & 
         &   ' time coeff of POC in sediments,                            sedlam      = ', sedlam
         WRITE(numout,*)     &
         &   ' Sediment geol loss for POC,                                sedlostpoc  = ', sedlostpoc
         WRITE(numout,*)     &
         &   ' Vert layer for diagnostic of vertical flux,                jpkp        = ', jpkb
!!
!! UKESM1 - new diagnostics  !! Jpalm; AXY (08/07/15)
         WRITE(numout,*) '=== UKESM1-related parameters ==='
         WRITE(numout,*) ' ---- --- ---'

         IF (lk_oasis) THEN
            WRITE(numout,*) '=== UKESM1 --  coupled DMS to the atmosphere'
            WRITE(numout,*)     &
            &   ' include DMS diagnostic?,                                   jdms        = ', jdms
            if (jdms_input .eq. 0) then
               WRITE(numout,*)     &
               &   ' use instant (0) or diel-avg (1) inputs,                    jdms_input  = instantaneous'
            else
               WRITE(numout,*)     &
               &   ' use instant (0) or diel-avg (1) inputs,                    jdms_input  = diel-average'
            endif
   	    if (jdms_model .eq. 1) then
               WRITE(numout,*)     &
               &   ' choice of DMS model passed to atmosphere,                  jdms_model  = Anderson et al. (2001)'
	    elseif (jdms_model .eq. 2) then
               WRITE(numout,*)     &
               &   ' choice of DMS model passed to atmosphere,                  jdms_model  = Simo & Dachs (2002)'
	    elseif (jdms_model .eq. 3) then
               WRITE(numout,*)     &
               &   ' choice of DMS model passed to atmosphere,                  jdms_model  = Aranami & Tsunogai (2004)'
	    elseif (jdms_model .eq. 4) then
               WRITE(numout,*)     &
               &   ' choice of DMS model passed to atmosphere,                  jdms_model  = Halloran et al. (2010)'
	    elseif (jdms_model .eq. 5) then
               WRITE(numout,*)     &
               &   ' choice of DMS model passed to atmosphere,                  jdms_model  = Anderson et al. (2001; modified)'
            endif
            if (jdms_model .eq. 1) then
               WRITE(numout,*)     &
               &   ' Anderson DMS model tuned parameters:                       DMS minimum = ',dmsmin,'. -- Default = 2.29 '
               WRITE(numout,*)     &
               &   ' Anderson DMS model tuned parameters:                       DMS cutoff  = ',dmscut,'. -- Default = 1.72 '
               WRITE(numout,*)     &
               &   ' Anderson DMS model tuned parameters:                       DMS slope   = ',dmsslp,'. -- Default = 8.24 '
            endif

            WRITE(numout,*) '=== UKESM1 --  coupled Chl to the atmosphere'
            WRITE(numout,*)        &
               &   ' Scaling factor to export tuned Chl to the atmosphere       scl_chl  = ', scl_chl
            IF (chl_out .eq. 1) THEN
               WRITE(numout,*)        &
               &   ' Chl field to be scaled and sent to the atmosphere:         chl_out  = Surface Chl field '
            ELSEIF (chl_out .eq. 2) THEN
               WRITE(numout,*)        &
               &   ' Chl field to be scaled and sent to the atmosphere:         chl_out  = MLD Chl field '
            ENDIF
         ENDIF ! IF lk_oasis=true
!!
      ENDIF
!!
!! Key depth positions (with thanks to Andrew Coward for bug-fixing this bit)
      DO jk = 1,jpk
         !! level thickness
         fthk  = e3t_1d(jk)
         !! level depth (top of level)
         fdep  = gdepw_1d(jk)
         !! level depth (bottom of level)
         fdep1 = fdep + fthk
         !!
         if (fdep .lt. 100.0 .AND. fdep1 .gt. 100.0) then        !  100 m
            i0100 = jk
         elseif (fdep .lt. 150.0 .AND. fdep1 .gt. 150.0) then    !  150 m (for BASIN)
            i0150 = jk
         elseif (fdep .lt. 200.0 .AND. fdep1 .gt. 200.0) then    !  200 m
            i0200 = jk
         elseif (fdep .lt. 500.0 .AND. fdep1 .gt. 500.0) then    !  500 m
            i0500 = jk
         elseif (fdep .lt. 1000.0 .AND. fdep1 .gt. 1000.0) then  ! 1000 m
            i1000 = jk
         elseif (fdep1 .lt. 1100.0) then                     ! 1100 m (for Moore et al. sedimentary iron)
            i1100 = jk
         endif
      enddo
      !!
      IF(lwp) THEN
          WRITE(numout,*) '=== Position of key depths'
          WRITE(numout,*)     & 
          &   ' jk position of  100 m horizon                              i0100       = ', i0100
          WRITE(numout,*)     &
          &   ' jk position of  150 m horizon                              i0150       = ', i0150
          WRITE(numout,*)     & 
          &   ' jk position of  200 m horizon                              i0200       = ', i0200
          WRITE(numout,*)     & 
          &   ' jk position of  500 m horizon                              i0500       = ', i0500
          WRITE(numout,*)     & 
          &   ' jk position of 1000 m horizon                              i1000       = ', i1000
          WRITE(numout,*)     & 
          &   ' jk position of 1100 m horizon [*]                          i1100       = ', i1100
          WRITE(numout,*) 'Got here ' , SIZE(friver_dep)
          CALL flush(numout)
      ENDIF

#if defined key_roam

      ! 1.4b namelist natroam : ROAM parameters
      ! ---------------------------------------
      
      xthetaphy = 0.
      xthetazoo = 0.
      xthetanit = 0.
      xthetarem = 0.
      xo2min    = 0.

      !READ(numnatm,natroam)
         ! Namelist natroam
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist natroam in reference namelist : MEDUSA diagnostics
         READ  ( numnatp_ref, natroam, IOSTAT = ios, ERR = 905)
905      IF( ios /= 0 ) CALL ctl_nam ( ios , 'natroam in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist natroam in configuration namelist : MEDUSA diagnostics
         READ  ( numnatp_cfg, natroam, IOSTAT = ios, ERR = 906 )
906      IF( ios /= 0 ) CALL ctl_nam ( ios , 'natroam in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, natroam )

!! ROAM carbon, alkalinity and oxygen cycle parameters
!!       xthetaphy :  oxygen evolution/consumption by phytoplankton
!!       xthetazoo :  oxygen consumption by zooplankton
!!       xthetanit :  oxygen consumption by nitrogen remineralisation
!!       xthetarem :  oxygen consumption by carbon remineralisation
!!       xo2min    :  oxygen minimum concentration

      IF(lwp) THEN
          WRITE(numout,*) 'natroam'
          WRITE(numout,*) ' '
!!
!! ROAM carbon, alkalinity and oxygen cycle parameters
          WRITE(numout,*) '=== ROAM carbon, alkalinity and oxygen cycle parameters'
          WRITE(numout,*)     &
          &   ' oxygen evolution/consumption by phytoplankton              xthetaphy   = ', xthetaphy
          WRITE(numout,*)     &
          &   ' oxygen consumption by zooplankton                          xthetazoo   = ', xthetazoo
          WRITE(numout,*)     &
          &   ' oxygen consumption by nitrogen remineralisation            xthetanit   = ', xthetanit
          WRITE(numout,*)     &
          &   ' oxygen consumption by carbon remineralisation              xthetarem   = ', xthetarem
          WRITE(numout,*)     &
          &   ' oxygen minimum concentration                               xo2min      = ', xo2min
       ENDIF

#endif

      CALL flush(numout)

      ! 1.5 namelist natopt : parameters for optic
      ! ------------------------------------------

      xkg0  = 0.
      xkr0  = 0.
      xkgp  = 0.
      xkrp  = 0.
      xlg   = 0.
      xlr   = 0.
      rpig  = 0.

      !READ(numnatm,natopt)
         ! Namelist natopt
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist natopt in reference namelist : MEDUSA diagnostics
         READ  ( numnatp_ref, natopt, IOSTAT = ios, ERR = 907)
907      IF( ios /= 0 ) CALL ctl_nam ( ios , 'natopt in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist natopt in configuration namelist : MEDUSA diagnostics
         READ  ( numnatp_cfg, natopt, IOSTAT = ios, ERR = 908 )
908      IF( ios /= 0 ) CALL ctl_nam ( ios , 'natopt in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, natopt )

      IF(lwp) THEN
         WRITE(numout,*) 'natopt'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' green   water absorption coeff  xkg0  = ',xkg0
         WRITE(numout,*) ' red water absorption coeff      xkr0  = ',xkr0
         WRITE(numout,*) ' pigment red absorption coeff    xkrp  = ',xkrp
         WRITE(numout,*) ' pigment green absorption coeff  xkgp  = ',xkgp
         WRITE(numout,*) ' green chl exposant              xlg   = ',xlg
         WRITE(numout,*) ' red   chl exposant              xlr   = ',xlr
         WRITE(numout,*) ' chla/chla+phea ratio            rpig  = ',rpig
         WRITE(numout,*) ' '

      ENDIF

      IF(lwp) THEN
         WRITE(numout,*) 'NaN check'
         WRITE(numout,*) ' '
         q1 = -1.
         q2 = 0.
         q3 = log(q1)
         write (numout,*) 'q3 = ', q3
         if ( ieee_is_nan( q3 ) ) then
            write (numout,*) 'NaN detected'
         else
            write (numout,*) 'NaN not detected'
         endif
         WRITE(numout,*) ' '
       ENDIF

   END SUBROUTINE trc_nam_medusa
   
   SUBROUTINE trc_nam_iom_medusa
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam_iom_medusa  ***
      !!
      !! ** Purpose : read all diag requested in iodef file through iom_use
      !!              So it is done only once 
      !!            ** All diagnostic MEDUSA could asked are registered in
      !!            the med_diag type with a boolean value
      !!            So if required, one diagnostic will be true.
      !!
      !!---------------------------------------------------------------------
      !!
      !!
      !!----------------------------------------------------------------------            
      !! Variable conventions
      !!----------------------------------------------------------------------
      !!
      IF (iom_use("INVTN")) THEN 
          med_diag%INVTN%dgsave = .TRUE.
      ELSE 
          med_diag%INVTN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INVTSI")) THEN 
          med_diag%INVTSI%dgsave = .TRUE.
      ELSE 
          med_diag%INVTSI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INVTFE")) THEN 
          med_diag%INVTFE%dgsave = .TRUE.
      ELSE 
          med_diag%INVTFE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PRN")) THEN 
          med_diag%PRN%dgsave = .TRUE.
      ELSE 
          med_diag%PRN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MPN")) THEN 
          med_diag%MPN%dgsave = .TRUE.
      ELSE 
          med_diag%MPN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PRD")) THEN 
          med_diag%PRD%dgsave = .TRUE.
      ELSE 
          med_diag%PRD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MPD")) THEN 
          med_diag%MPD%dgsave = .TRUE.
      ELSE 
          med_diag%MPD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DSED")) THEN 
          med_diag%DSED%dgsave = .TRUE.
      ELSE 
          med_diag%DSED%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OPAL")) THEN 
          med_diag%OPAL%dgsave = .TRUE.
      ELSE 
          med_diag%OPAL%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OPALDISS")) THEN 
          med_diag%OPALDISS%dgsave = .TRUE.
      ELSE 
          med_diag%OPALDISS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMIPn")) THEN 
          med_diag%GMIPn%dgsave = .TRUE.
      ELSE 
          med_diag%GMIPn%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMID")) THEN 
          med_diag%GMID%dgsave = .TRUE.
      ELSE 
          med_diag%GMID%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MZMI")) THEN 
          med_diag%MZMI%dgsave = .TRUE.
      ELSE 
          med_diag%MZMI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMEPN")) THEN 
          med_diag%GMEPN%dgsave = .TRUE.
      ELSE 
          med_diag%GMEPN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMEPD")) THEN 
          med_diag%GMEPD%dgsave = .TRUE.
      ELSE 
          med_diag%GMEPD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMEZMI")) THEN 
          med_diag%GMEZMI%dgsave = .TRUE.
      ELSE 
          med_diag%GMEZMI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMED")) THEN 
          med_diag%GMED%dgsave = .TRUE.
      ELSE 
          med_diag%GMED%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MZME")) THEN 
          med_diag%MZME%dgsave = .TRUE.
      ELSE 
          med_diag%MZME%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DEXP")) THEN 
          med_diag%DEXP%dgsave = .TRUE.
      ELSE 
          med_diag%DEXP%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DETN")) THEN 
          med_diag%DETN%dgsave = .TRUE.
      ELSE 
          med_diag%DETN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MDET")) THEN 
          med_diag%MDET%dgsave = .TRUE.
      ELSE 
          med_diag%MDET%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("AEOLIAN")) THEN 
          med_diag%AEOLIAN%dgsave = .TRUE.
      ELSE 
          med_diag%AEOLIAN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BENTHIC")) THEN 
          med_diag%BENTHIC%dgsave = .TRUE.
      ELSE 
          med_diag%BENTHIC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SCAVENGE")) THEN 
          med_diag%SCAVENGE%dgsave = .TRUE.
      ELSE 
          med_diag%SCAVENGE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PN_JLIM")) THEN 
          med_diag%PN_JLIM%dgsave = .TRUE.
      ELSE 
          med_diag%PN_JLIM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PN_NLIM")) THEN 
          med_diag%PN_NLIM%dgsave = .TRUE.
      ELSE 
          med_diag%PN_NLIM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PN_FELIM")) THEN 
          med_diag%PN_FELIM%dgsave = .TRUE.
      ELSE 
          med_diag%PN_FELIM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PD_JLIM")) THEN 
          med_diag%PD_JLIM%dgsave = .TRUE.
      ELSE 
          med_diag%PD_JLIM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PD_NLIM")) THEN 
          med_diag%PD_NLIM%dgsave = .TRUE.
      ELSE 
          med_diag%PD_NLIM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PD_FELIM")) THEN 
          med_diag%PD_FELIM%dgsave = .TRUE.
      ELSE 
          med_diag%PD_FELIM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PD_SILIM")) THEN 
          med_diag%PD_SILIM%dgsave = .TRUE.
      ELSE 
          med_diag%PD_SILIM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PDSILIM2")) THEN 
          med_diag%PDSILIM2%dgsave = .TRUE.
      ELSE 
          med_diag%PDSILIM2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDT__100")) THEN 
          med_diag%SDT__100%dgsave = .TRUE.
      ELSE 
          med_diag%SDT__100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDT__200")) THEN 
          med_diag%SDT__200%dgsave = .TRUE.
      ELSE 
          med_diag%SDT__200%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDT__500")) THEN 
          med_diag%SDT__500%dgsave = .TRUE.
      ELSE 
          med_diag%SDT__500%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDT_1000")) THEN 
          med_diag%SDT_1000%dgsave = .TRUE.
      ELSE 
          med_diag%SDT_1000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("TOTREG_N")) THEN 
          med_diag%TOTREG_N%dgsave = .TRUE.
      ELSE 
          med_diag%TOTREG_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("TOTRG_SI")) THEN 
          med_diag%TOTRG_SI%dgsave = .TRUE.
      ELSE 
          med_diag%TOTRG_SI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REG__100")) THEN 
          med_diag%REG__100%dgsave = .TRUE.
      ELSE 
          med_diag%REG__100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REG__200")) THEN 
          med_diag%REG__200%dgsave = .TRUE.
      ELSE 
          med_diag%REG__200%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REG__500")) THEN 
          med_diag%REG__500%dgsave = .TRUE.
      ELSE 
          med_diag%REG__500%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REG_1000")) THEN 
          med_diag%REG_1000%dgsave = .TRUE.
      ELSE 
          med_diag%REG_1000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FASTN")) THEN 
          med_diag%FASTN%dgsave = .TRUE.
      ELSE 
          med_diag%FASTN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FASTSI")) THEN 
          med_diag%FASTSI%dgsave = .TRUE.
      ELSE 
          med_diag%FASTSI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FASTFE")) THEN 
          med_diag%FASTFE%dgsave = .TRUE.
      ELSE 
          med_diag%FASTFE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FASTC")) THEN 
          med_diag%FASTC%dgsave = .TRUE.
      ELSE 
          med_diag%FASTC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FASTCA")) THEN 
          med_diag%FASTCA%dgsave = .TRUE.
      ELSE 
          med_diag%FASTCA%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDT__100")) THEN 
          med_diag%FDT__100%dgsave = .TRUE.
      ELSE 
          med_diag%FDT__100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDT__200")) THEN 
          med_diag%FDT__200%dgsave = .TRUE.
      ELSE 
          med_diag%FDT__200%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDT__500")) THEN 
          med_diag%FDT__500%dgsave = .TRUE.
      ELSE 
          med_diag%FDT__500%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDT_1000")) THEN 
          med_diag%FDT_1000%dgsave = .TRUE.
      ELSE 
          med_diag%FDT_1000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RG__100F")) THEN 
          med_diag%RG__100F%dgsave = .TRUE.
      ELSE 
          med_diag%RG__100F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RG__200F")) THEN 
          med_diag%RG__200F%dgsave = .TRUE.
      ELSE 
          med_diag%RG__200F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RG__500F")) THEN 
          med_diag%RG__500F%dgsave = .TRUE.
      ELSE 
          med_diag%RG__500F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RG_1000F")) THEN 
          med_diag%RG_1000F%dgsave = .TRUE.
      ELSE 
          med_diag%RG_1000F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDS__100")) THEN 
          med_diag%FDS__100%dgsave = .TRUE.
      ELSE 
          med_diag%FDS__100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDS__200")) THEN 
          med_diag%FDS__200%dgsave = .TRUE.
      ELSE 
          med_diag%FDS__200%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDS__500")) THEN 
          med_diag%FDS__500%dgsave = .TRUE.
      ELSE 
          med_diag%FDS__500%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FDS_1000")) THEN 
          med_diag%FDS_1000%dgsave = .TRUE.
      ELSE 
          med_diag%FDS_1000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RGS_100F")) THEN 
          med_diag%RGS_100F%dgsave = .TRUE.
      ELSE 
          med_diag%RGS_100F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RGS_200F")) THEN 
          med_diag%RGS_200F%dgsave = .TRUE.
      ELSE 
          med_diag%RGS_200F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RGS_500F")) THEN 
          med_diag%RGS_500F%dgsave = .TRUE.
      ELSE 
          med_diag%RGS_500F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RGS1000F")) THEN 
          med_diag%RGS1000F%dgsave = .TRUE.
      ELSE 
          med_diag%RGS1000F%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REMINN")) THEN 
          med_diag%REMINN%dgsave = .TRUE.
      ELSE 
          med_diag%REMINN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REMINSI")) THEN 
          med_diag%REMINSI%dgsave = .TRUE.
      ELSE 
          med_diag%REMINSI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REMINFE")) THEN 
          med_diag%REMINFE%dgsave = .TRUE.
      ELSE 
          med_diag%REMINFE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REMINC")) THEN 
          med_diag%REMINC%dgsave = .TRUE.
      ELSE 
          med_diag%REMINC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REMINCA")) THEN 
          med_diag%REMINCA%dgsave = .TRUE.
      ELSE 
          med_diag%REMINCA%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SEAFLRN")) THEN 
          med_diag%SEAFLRN%dgsave = .TRUE.
      ELSE 
          med_diag%SEAFLRN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SEAFLRSI")) THEN 
          med_diag%SEAFLRSI%dgsave = .TRUE.
      ELSE 
          med_diag%SEAFLRSI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SEAFLRFE")) THEN 
          med_diag%SEAFLRFE%dgsave = .TRUE.
      ELSE 
          med_diag%SEAFLRFE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SEAFLRC")) THEN 
          med_diag%SEAFLRC%dgsave = .TRUE.
      ELSE 
          med_diag%SEAFLRC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SEAFLRCA")) THEN 
          med_diag%SEAFLRCA%dgsave = .TRUE.
      ELSE 
          med_diag%SEAFLRCA%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MED_QSR")) THEN 
          med_diag%MED_QSR%dgsave = .TRUE.
      ELSE 
          med_diag%MED_QSR%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MED_XPAR")) THEN 
          med_diag%MED_XPAR%dgsave = .TRUE.
      ELSE 
          med_diag%MED_XPAR%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INTFLX_N")) THEN 
          med_diag%INTFLX_N%dgsave = .TRUE.
      ELSE 
          med_diag%INTFLX_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INTFLX_SI")) THEN 
          med_diag%INTFLX_SI%dgsave = .TRUE.
      ELSE 
          med_diag%INTFLX_SI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INTFLX_FE")) THEN 
          med_diag%INTFLX_FE%dgsave = .TRUE.
      ELSE 
          med_diag%INTFLX_FE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INT_PN")) THEN 
          med_diag%INT_PN%dgsave = .TRUE.
      ELSE 
          med_diag%INT_PN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INT_PD")) THEN 
          med_diag%INT_PD%dgsave = .TRUE.
      ELSE 
          med_diag%INT_PD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ML_PRN")) THEN 
          med_diag%ML_PRN%dgsave = .TRUE.
      ELSE 
          med_diag%ML_PRN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ML_PRD")) THEN 
          med_diag%ML_PRD%dgsave = .TRUE.
      ELSE 
          med_diag%ML_PRD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCAL_CCD")) THEN 
          med_diag%OCAL_CCD%dgsave = .TRUE.
      ELSE 
          med_diag%OCAL_CCD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCAL_LVL")) THEN 
          med_diag%OCAL_LVL%dgsave = .TRUE.
      ELSE 
          med_diag%OCAL_LVL%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FE_0000")) THEN 
          med_diag%FE_0000%dgsave = .TRUE.
      ELSE 
          med_diag%FE_0000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FE_0100")) THEN 
          med_diag%FE_0100%dgsave = .TRUE.
      ELSE 
          med_diag%FE_0100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FE_0200")) THEN 
          med_diag%FE_0200%dgsave = .TRUE.
      ELSE 
          med_diag%FE_0200%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FE_0500")) THEN 
          med_diag%FE_0500%dgsave = .TRUE.
      ELSE 
          med_diag%FE_0500%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FE_1000")) THEN 
          med_diag%FE_1000%dgsave = .TRUE.
      ELSE 
          med_diag%FE_1000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MED_XZE")) THEN 
          med_diag%MED_XZE%dgsave = .TRUE.
      ELSE 
          med_diag%MED_XZE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("WIND")) THEN 
          med_diag%WIND%dgsave = .TRUE.
      ELSE 
          med_diag%WIND%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ATM_PCO2")) THEN 
          med_diag%ATM_PCO2%dgsave = .TRUE.
      ELSE 
          med_diag%ATM_PCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_PH")) THEN 
          med_diag%OCN_PH%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_PH%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_PCO2")) THEN 
          med_diag%OCN_PCO2%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_PCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCNH2CO3")) THEN 
          med_diag%OCNH2CO3%dgsave = .TRUE.
      ELSE 
          med_diag%OCNH2CO3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_HCO3")) THEN 
          med_diag%OCN_HCO3%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_HCO3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_CO3")) THEN 
          med_diag%OCN_CO3%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_CO3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("CO2FLUX")) THEN 
          med_diag%CO2FLUX%dgsave = .TRUE.
      ELSE 
          med_diag%CO2FLUX%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OM_CAL")) THEN 
          med_diag%OM_CAL%dgsave = .TRUE.
      ELSE 
          med_diag%OM_CAL%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OM_ARG")) THEN 
          med_diag%OM_ARG%dgsave = .TRUE.
      ELSE 
          med_diag%OM_ARG%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("TCO2")) THEN 
          med_diag%TCO2%dgsave = .TRUE.
      ELSE 
          med_diag%TCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("TALK")) THEN 
          med_diag%TALK%dgsave = .TRUE.
      ELSE 
          med_diag%TALK%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("KW660")) THEN 
          med_diag%KW660%dgsave = .TRUE.
      ELSE 
          med_diag%KW660%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ATM_PP0")) THEN 
          med_diag%ATM_PP0%dgsave = .TRUE.
      ELSE 
          med_diag%ATM_PP0%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("O2FLUX")) THEN 
          med_diag%O2FLUX%dgsave = .TRUE.
      ELSE 
          med_diag%O2FLUX%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("O2SAT")) THEN 
          med_diag%O2SAT%dgsave = .TRUE.
      ELSE 
          med_diag%O2SAT%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("CAL_CCD")) THEN 
          med_diag%CAL_CCD%dgsave = .TRUE.
      ELSE 
          med_diag%CAL_CCD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ARG_CCD")) THEN 
          med_diag%ARG_CCD%dgsave = .TRUE.
      ELSE 
          med_diag%ARG_CCD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SFR_OCAL")) THEN 
          med_diag%SFR_OCAL%dgsave = .TRUE.
      ELSE 
          med_diag%SFR_OCAL%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SFR_OARG")) THEN 
          med_diag%SFR_OARG%dgsave = .TRUE.
      ELSE 
          med_diag%SFR_OARG%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("N_PROD")) THEN 
          med_diag%N_PROD%dgsave = .TRUE.
      ELSE 
          med_diag%N_PROD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("N_CONS")) THEN 
          med_diag%N_CONS%dgsave = .TRUE.
      ELSE 
          med_diag%N_CONS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("C_PROD")) THEN 
          med_diag%C_PROD%dgsave = .TRUE.
      ELSE 
          med_diag%C_PROD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("C_CONS")) THEN 
          med_diag%C_CONS%dgsave = .TRUE.
      ELSE 
          med_diag%C_CONS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("O2_PROD")) THEN 
          med_diag%O2_PROD%dgsave = .TRUE.
      ELSE 
          med_diag%O2_PROD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("O2_CONS")) THEN 
          med_diag%O2_CONS%dgsave = .TRUE.
      ELSE 
          med_diag%O2_CONS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("O2_ANOX")) THEN 
          med_diag%O2_ANOX%dgsave = .TRUE.
      ELSE 
          med_diag%O2_ANOX%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RR_0100")) THEN 
          med_diag%RR_0100%dgsave = .TRUE.
      ELSE 
          med_diag%RR_0100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RR_0500")) THEN 
          med_diag%RR_0500%dgsave = .TRUE.
      ELSE 
          med_diag%RR_0500%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RR_1000")) THEN 
          med_diag%RR_1000%dgsave = .TRUE.
      ELSE 
          med_diag%RR_1000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("IBEN_N")) THEN 
          med_diag%IBEN_N%dgsave = .TRUE.
      ELSE 
          med_diag%IBEN_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("IBEN_FE")) THEN 
          med_diag%IBEN_FE%dgsave = .TRUE.
      ELSE 
          med_diag%IBEN_FE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("IBEN_C")) THEN 
          med_diag%IBEN_C%dgsave = .TRUE.
      ELSE 
          med_diag%IBEN_C%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("IBEN_SI")) THEN 
          med_diag%IBEN_SI%dgsave = .TRUE.
      ELSE 
          med_diag%IBEN_SI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("IBEN_CA")) THEN 
          med_diag%IBEN_CA%dgsave = .TRUE.
      ELSE 
          med_diag%IBEN_CA%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OBEN_N")) THEN 
          med_diag%OBEN_N%dgsave = .TRUE.
      ELSE 
          med_diag%OBEN_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OBEN_FE")) THEN 
          med_diag%OBEN_FE%dgsave = .TRUE.
      ELSE 
          med_diag%OBEN_FE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OBEN_C")) THEN 
          med_diag%OBEN_C%dgsave = .TRUE.
      ELSE 
          med_diag%OBEN_C%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OBEN_SI")) THEN 
          med_diag%OBEN_SI%dgsave = .TRUE.
      ELSE 
          med_diag%OBEN_SI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OBEN_CA")) THEN 
          med_diag%OBEN_CA%dgsave = .TRUE.
      ELSE 
          med_diag%OBEN_CA%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BEN_N")) THEN 
          med_diag%BEN_N%dgsave = .TRUE.
      ELSE 
          med_diag%BEN_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BEN_FE")) THEN 
          med_diag%BEN_FE%dgsave = .TRUE.
      ELSE 
          med_diag%BEN_FE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BEN_C")) THEN 
          med_diag%BEN_C%dgsave = .TRUE.
      ELSE 
          med_diag%BEN_C%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BEN_SI")) THEN 
          med_diag%BEN_SI%dgsave = .TRUE.
      ELSE 
          med_diag%BEN_SI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BEN_CA")) THEN 
          med_diag%BEN_CA%dgsave = .TRUE.
      ELSE 
          med_diag%BEN_CA%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RUNOFF")) THEN 
          med_diag%RUNOFF%dgsave = .TRUE.
      ELSE 
          med_diag%RUNOFF%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RIV_N")) THEN 
          med_diag%RIV_N%dgsave = .TRUE.
      ELSE 
          med_diag%RIV_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RIV_SI")) THEN 
          med_diag%RIV_SI%dgsave = .TRUE.
      ELSE 
          med_diag%RIV_SI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RIV_C")) THEN 
          med_diag%RIV_C%dgsave = .TRUE.
      ELSE 
          med_diag%RIV_C%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("RIV_ALK")) THEN 
          med_diag%RIV_ALK%dgsave = .TRUE.
      ELSE 
          med_diag%RIV_ALK%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DETC")) THEN 
          med_diag%DETC%dgsave = .TRUE.
      ELSE 
          med_diag%DETC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDC__100")) THEN 
          med_diag%SDC__100%dgsave = .TRUE.
      ELSE 
          med_diag%SDC__100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDC__200")) THEN 
          med_diag%SDC__200%dgsave = .TRUE.
      ELSE 
          med_diag%SDC__200%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDC__500")) THEN 
          med_diag%SDC__500%dgsave = .TRUE.
      ELSE 
          med_diag%SDC__500%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("SDC_1000")) THEN 
          med_diag%SDC_1000%dgsave = .TRUE.
      ELSE 
          med_diag%SDC_1000%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INVTC")) THEN 
          med_diag%INVTC%dgsave = .TRUE.
      ELSE 
          med_diag%INVTC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INVTALK")) THEN 
          med_diag%INVTALK%dgsave = .TRUE.
      ELSE 
          med_diag%INVTALK%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INVTO2")) THEN 
          med_diag%INVTO2%dgsave = .TRUE.
      ELSE 
          med_diag%INVTO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("LYSO_CA")) THEN 
          med_diag%LYSO_CA%dgsave = .TRUE.
      ELSE 
          med_diag%LYSO_CA%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("COM_RESP")) THEN 
          med_diag%COM_RESP%dgsave = .TRUE.
      ELSE 
          med_diag%COM_RESP%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PN_LLOSS")) THEN 
          med_diag%PN_LLOSS%dgsave = .TRUE.
      ELSE 
          med_diag%PN_LLOSS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PD_LLOSS")) THEN 
          med_diag%PD_LLOSS%dgsave = .TRUE.
      ELSE 
          med_diag%PD_LLOSS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_LLOSS")) THEN 
          med_diag%ZI_LLOSS%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_LLOSS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_LLOSS")) THEN 
          med_diag%ZE_LLOSS%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_LLOSS%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_MES_N")) THEN 
          med_diag%ZI_MES_N%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_MES_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_MES_D")) THEN 
          med_diag%ZI_MES_D%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_MES_D%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_MES_C")) THEN 
          med_diag%ZI_MES_C%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_MES_C%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_MESDC")) THEN 
          med_diag%ZI_MESDC%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_MESDC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_EXCR")) THEN 
          med_diag%ZI_EXCR%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_EXCR%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_RESP")) THEN 
          med_diag%ZI_RESP%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_RESP%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZI_GROW")) THEN 
          med_diag%ZI_GROW%dgsave = .TRUE.
      ELSE 
          med_diag%ZI_GROW%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_MES_N")) THEN 
          med_diag%ZE_MES_N%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_MES_N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_MES_D")) THEN 
          med_diag%ZE_MES_D%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_MES_D%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_MES_C")) THEN 
          med_diag%ZE_MES_C%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_MES_C%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_MESDC")) THEN 
          med_diag%ZE_MESDC%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_MESDC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_EXCR")) THEN 
          med_diag%ZE_EXCR%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_EXCR%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_RESP")) THEN 
          med_diag%ZE_RESP%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_RESP%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZE_GROW")) THEN 
          med_diag%ZE_GROW%dgsave = .TRUE.
      ELSE 
          med_diag%ZE_GROW%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MDETC")) THEN 
          med_diag%MDETC%dgsave = .TRUE.
      ELSE 
          med_diag%MDETC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMIDC")) THEN 
          med_diag%GMIDC%dgsave = .TRUE.
      ELSE 
          med_diag%GMIDC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("GMEDC")) THEN 
          med_diag%GMEDC%dgsave = .TRUE.
      ELSE 
          med_diag%GMEDC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INT_ZMI")) THEN 
          med_diag%INT_ZMI%dgsave = .TRUE.
      ELSE 
          med_diag%INT_ZMI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INT_ZME")) THEN 
          med_diag%INT_ZME%dgsave = .TRUE.
      ELSE 
          med_diag%INT_ZME%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INT_DET")) THEN 
          med_diag%INT_DET%dgsave = .TRUE.
      ELSE 
          med_diag%INT_DET%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INT_DTC")) THEN 
          med_diag%INT_DTC%dgsave = .TRUE.
      ELSE 
          med_diag%INT_DTC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DMS_SURF")) THEN 
          med_diag%DMS_SURF%dgsave = .TRUE.
      ELSE 
          med_diag%DMS_SURF%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DMS_ANDR")) THEN 
          med_diag%DMS_ANDR%dgsave = .TRUE.
      ELSE 
          med_diag%DMS_ANDR%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DMS_SIMO")) THEN 
          med_diag%DMS_SIMO%dgsave = .TRUE.
      ELSE 
          med_diag%DMS_SIMO%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DMS_ARAN")) THEN 
          med_diag%DMS_ARAN%dgsave = .TRUE.
      ELSE 
          med_diag%DMS_ARAN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DMS_HALL")) THEN 
          med_diag%DMS_HALL%dgsave = .TRUE.
      ELSE 
          med_diag%DMS_HALL%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DMS_ANDM")) THEN 
          med_diag%DMS_ANDM%dgsave = .TRUE.
      ELSE 
          med_diag%DMS_ANDM%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ATM_XCO2")) THEN 
          med_diag%ATM_XCO2%dgsave = .TRUE.
      ELSE 
          med_diag%ATM_XCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_FCO2")) THEN 
          med_diag%OCN_FCO2%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_FCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ATM_FCO2")) THEN 
          med_diag%ATM_FCO2%dgsave = .TRUE.
      ELSE 
          med_diag%ATM_FCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_RHOSW")) THEN 
          med_diag%OCN_RHOSW%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_RHOSW%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_SCHCO2")) THEN 
          med_diag%OCN_SCHCO2%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_SCHCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_KWCO2")) THEN 
          med_diag%OCN_KWCO2%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_KWCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_K0")) THEN 
          med_diag%OCN_K0%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_K0%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("CO2STARAIR")) THEN 
          med_diag%CO2STARAIR%dgsave = .TRUE.
      ELSE 
          med_diag%CO2STARAIR%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OCN_DPCO2")) THEN 
          med_diag%OCN_DPCO2%dgsave = .TRUE.
      ELSE 
          med_diag%OCN_DPCO2%dgsave = .FALSE.
      ENDIF
      !! UKESM additional
      IF  (iom_use("CHL_MLD")) THEN 
          med_diag%CHL_MLD%dgsave = .TRUE.
      ELSE 
          med_diag%CHL_MLD%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("CHL_CPL")) THEN 
          med_diag%CHL_CPL%dgsave = .TRUE.
      ELSE 
          med_diag%CHL_CPL%dgsave = .FALSE.
      ENDIF
      !! 3D
      IF  (iom_use("TPP3")) THEN 
          med_diag%TPP3%dgsave = .TRUE.
      ELSE 
          med_diag%TPP3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DETFLUX3")) THEN 
          med_diag%DETFLUX3%dgsave = .TRUE.
      ELSE 
          med_diag%DETFLUX3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REMIN3N")) THEN 
          med_diag%REMIN3N%dgsave = .TRUE.
      ELSE 
          med_diag%REMIN3N%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PH3")) THEN 
          med_diag%PH3%dgsave = .TRUE.
      ELSE 
          med_diag%PH3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("OM_CAL3")) THEN 
          med_diag%OM_CAL3%dgsave = .TRUE.
      ELSE 
          med_diag%OM_CAL3%dgsave = .FALSE.
      ENDIF
      !!
      !!----------------------------------------------------------------------
      !! AXY (03/11/16): add in additional CMIP6 diagnostics
      !!----------------------------------------------------------------------
      !!
      !! 2D fields
      IF  (iom_use("epC100")) THEN 
          med_diag%epC100%dgsave = .TRUE.
      ELSE 
          med_diag%epC100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("epCALC100")) THEN 
          med_diag%epCALC100%dgsave = .TRUE.
      ELSE 
          med_diag%epCALC100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("epN100")) THEN 
          med_diag%epN100%dgsave = .TRUE.
      ELSE 
          med_diag%epN100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("epSI100")) THEN 
          med_diag%epSI100%dgsave = .TRUE.
      ELSE 
          med_diag%epSI100%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FGCO2")) THEN 
          med_diag%FGCO2%dgsave = .TRUE.
      ELSE 
          med_diag%FGCO2%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INTDISSIC")) THEN 
          med_diag%INTDISSIC%dgsave = .TRUE.
      ELSE 
          med_diag%INTDISSIC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INTDISSIN")) THEN 
          med_diag%INTDISSIN%dgsave = .TRUE.
      ELSE 
          med_diag%INTDISSIN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INTDISSISI")) THEN 
          med_diag%INTDISSISI%dgsave = .TRUE.
      ELSE 
          med_diag%INTDISSISI%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("INTTALK")) THEN 
          med_diag%INTTALK%dgsave = .TRUE.
      ELSE 
          med_diag%INTTALK%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("O2min")) THEN 
          med_diag%O2min%dgsave = .TRUE.
      ELSE 
          med_diag%O2min%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("ZO2min")) THEN 
          med_diag%ZO2min%dgsave = .TRUE.
      ELSE 
          med_diag%ZO2min%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FBDDTALK")) THEN 
          med_diag%FBDDTALK%dgsave = .TRUE.
      ELSE 
          med_diag%FBDDTALK%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FBDDTDIC")) THEN 
          med_diag%FBDDTDIC%dgsave = .TRUE.
      ELSE 
          med_diag%FBDDTDIC%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FBDDTDIFE")) THEN 
          med_diag%FBDDTDIFE%dgsave = .TRUE.
      ELSE 
          med_diag%FBDDTDIFE%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FBDDTDIN")) THEN 
          med_diag%FBDDTDIN%dgsave = .TRUE.
      ELSE 
          med_diag%FBDDTDIN%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FBDDTDISI")) THEN 
          med_diag%FBDDTDISI%dgsave = .TRUE.
      ELSE 
          med_diag%FBDDTDISI%dgsave = .FALSE.
      ENDIF
      !!
      !! 3D
      IF  (iom_use("TPPD3")) THEN 
          med_diag%TPPD3%dgsave = .TRUE.
      ELSE 
          med_diag%TPPD3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BDDTALK3")) THEN 
          med_diag%BDDTALK3%dgsave = .TRUE.
      ELSE 
          med_diag%BDDTALK3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BDDTDIC3")) THEN 
          med_diag%BDDTDIC3%dgsave = .TRUE.
      ELSE 
          med_diag%BDDTDIC3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BDDTDIFE3")) THEN 
          med_diag%BDDTDIFE3%dgsave = .TRUE.
      ELSE 
          med_diag%BDDTDIFE3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BDDTDIN3")) THEN 
          med_diag%BDDTDIN3%dgsave = .TRUE.
      ELSE 
          med_diag%BDDTDIN3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("BDDTDISI3")) THEN 
          med_diag%BDDTDISI3%dgsave = .TRUE.
      ELSE 
          med_diag%BDDTDISI3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FD_NIT3")) THEN 
          med_diag%FD_NIT3%dgsave = .TRUE.
      ELSE 
          med_diag%FD_NIT3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FD_SIL3")) THEN 
          med_diag%FD_SIL3%dgsave = .TRUE.
      ELSE 
          med_diag%FD_SIL3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FD_CAR3")) THEN 
          med_diag%FD_CAR3%dgsave = .TRUE.
      ELSE 
          med_diag%FD_CAR3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FD_CAL3")) THEN 
          med_diag%FD_CAL3%dgsave = .TRUE.
      ELSE 
          med_diag%FD_CAL3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("CO33")) THEN 
          med_diag%CO33%dgsave = .TRUE.
      ELSE 
          med_diag%CO33%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("CO3SATARAG3")) THEN 
          med_diag%CO3SATARAG3%dgsave = .TRUE.
      ELSE 
          med_diag%CO3SATARAG3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("CO3SATCALC3")) THEN 
          med_diag%CO3SATCALC3%dgsave = .TRUE.
      ELSE 
          med_diag%CO3SATCALC3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("DCALC3")) THEN 
          med_diag%DCALC3%dgsave = .TRUE.
      ELSE 
          med_diag%DCALC3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("EXPC3")) THEN 
          med_diag%EXPC3%dgsave = .TRUE.
      ELSE 
          med_diag%EXPC3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("EXPN3")) THEN 
          med_diag%EXPN3%dgsave = .TRUE.
      ELSE 
          med_diag%EXPN3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FEDISS3")) THEN 
          med_diag%FEDISS3%dgsave = .TRUE.
      ELSE 
          med_diag%FEDISS3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("FESCAV3")) THEN 
          med_diag%FESCAV3%dgsave = .TRUE.
      ELSE 
          med_diag%FESCAV3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MIGRAZP3")) THEN 
          med_diag%MIGRAZP3%dgsave = .TRUE.
      ELSE 
          med_diag%MIGRAZP3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MIGRAZD3")) THEN 
          med_diag%MIGRAZD3%dgsave = .TRUE.
      ELSE 
          med_diag%MIGRAZD3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MEGRAZP3")) THEN 
          med_diag%MEGRAZP3%dgsave = .TRUE.
      ELSE 
          med_diag%MEGRAZP3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MEGRAZD3")) THEN 
          med_diag%MEGRAZD3%dgsave = .TRUE.
      ELSE 
          med_diag%MEGRAZD3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("MEGRAZZ3")) THEN 
          med_diag%MEGRAZZ3%dgsave = .TRUE.
      ELSE 
          med_diag%MEGRAZZ3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("O2SAT3")) THEN 
          med_diag%O2SAT3%dgsave = .TRUE.
      ELSE 
          med_diag%O2SAT3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PBSI3")) THEN 
          med_diag%PBSI3%dgsave = .TRUE.
      ELSE 
          med_diag%PBSI3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PCAL3")) THEN 
          med_diag%PCAL3%dgsave = .TRUE.
      ELSE 
          med_diag%PCAL3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("REMOC3")) THEN 
          med_diag%REMOC3%dgsave = .TRUE.
      ELSE 
          med_diag%REMOC3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PNLIMJ3")) THEN 
          med_diag%PNLIMJ3%dgsave = .TRUE.
      ELSE 
          med_diag%PNLIMJ3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PNLIMN3")) THEN 
          med_diag%PNLIMN3%dgsave = .TRUE.
      ELSE 
          med_diag%PNLIMN3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PNLIMFE3")) THEN 
          med_diag%PNLIMFE3%dgsave = .TRUE.
      ELSE 
          med_diag%PNLIMFE3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PDLIMJ3")) THEN 
          med_diag%PDLIMJ3%dgsave = .TRUE.
      ELSE 
          med_diag%PDLIMJ3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PDLIMN3")) THEN 
          med_diag%PDLIMN3%dgsave = .TRUE.
      ELSE 
          med_diag%PDLIMN3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PDLIMFE3")) THEN 
          med_diag%PDLIMFE3%dgsave = .TRUE.
      ELSE 
          med_diag%PDLIMFE3%dgsave = .FALSE.
      ENDIF
      IF  (iom_use("PDLIMSI3")) THEN 
          med_diag%PDLIMSI3%dgsave = .TRUE.
      ELSE 
          med_diag%PDLIMSI3%dgsave = .FALSE.
      ENDIF

   END SUBROUTINE   trc_nam_iom_medusa
   !/CEB !$AGRIF_END_DO_NOT_TREAT   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No MEDUSA
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_medusa                      ! Empty routine
   END  SUBROUTINE  trc_nam_medusa
#endif  

   !!======================================================================
END MODULE trcnam_medusa
