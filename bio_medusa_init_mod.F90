MODULE bio_medusa_init_mod
   !!======================================================================
   !!                         ***  MODULE bio_medusa_init  ***
   !! Initialisation for TRC_BIO_MEDUSA
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!   -   ! 2017-08 (A. Yool)            Add slow-sinking detrius variables
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   bio_medusa_init    ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   !CEB!$AGRIF_DO_NOT_TREAT
   SUBROUTINE bio_medusa_init( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE bio_medusa_init  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - zeroes arrays used in TRC_BIO_MEDUSA
      !!  - sets up the arrays needed for the diagnostics
      !!---------------------------------------------------------------------- 
      USE bio_medusa_mod
      USE par_oce,           ONLY: jpi, jpj, jpk
      USE sms_medusa,        ONLY: jdms
      USE trc,               ONLY: ln_diatrc, med_diag, nittrc000 
      USE in_out_manager,    ONLY: lwp, numout

      USE iom,               ONLY: lk_iomput
      USE trcnam_medusa,     ONLY: trc_nam_iom_medusa

      !! time (integer timestep)
      INTEGER, INTENT( in ) ::    kt

      IF( ln_diatrc ) THEN
         IF (lwp) write (numout,*) 'Diagnostics are now ALL through XIOS (key_xios)'
         IF (lwp) write (numout,*) 'No more key_diatrc anymore.'
      ENDIF

      !!----------------------------------------------------------------------
      !! Zero fast-sinking detritus 2D fields
      !!----------------------------------------------------------------------
      !!
      ffastn(:,:)  = 0.0        !! organic nitrogen
      ffastsi(:,:) = 0.0        !! biogenic silicon
      ffastfe(:,:) = 0.0        !! organic iron
      ffastc(:,:)  = 0.0        !! organic carbon
      ffastca(:,:) = 0.0        !! biogenic calcium carbonate
      !!
      fsedn(:,:)   = 0.0        !! Seafloor flux of N 
      fsedsi(:,:)  = 0.0        !! Seafloor flux of Si
      fsedfe(:,:)  = 0.0        !! Seafloor flux of Fe
      fsedc(:,:)   = 0.0        !! Seafloor flux of C
      fsedca(:,:)  = 0.0        !! Seafloor flux of CaCO3
      !!
      fregenfast(:,:)   = 0.0   !! integrated  N regeneration (fast detritus)
      fregenfastsi(:,:) = 0.0   !! integrated Si regeneration (fast detritus)
# if defined key_roam
      fregenfastc(:,:)  = 0.0   !! integrated  C regeneration (fast detritus)
# endif
      !!
      fccd(:,:)    = 0.0        !! last depth level before CCD

      !!----------------------------------------------------------------------
      !! blank nutrient/flux inventories
      !!----------------------------------------------------------------------
      !!
      fflx_n(:,:)  = 0.0        !! nitrogen flux total
      fflx_si(:,:) = 0.0        !! silicon  flux total
      fflx_fe(:,:) = 0.0        !! iron     flux total
      fifd_n(:,:)  = 0.0        !! nitrogen fast detritus production
      fifd_si(:,:) = 0.0        !! silicon  fast detritus production
      fifd_fe(:,:) = 0.0        !! iron     fast detritus production
      fofd_n(:,:)  = 0.0        !! nitrogen fast detritus remineralisation
      fofd_si(:,:) = 0.0        !! silicon  fast detritus remineralisation
      fofd_fe(:,:) = 0.0        !! iron     fast detritus remineralisation
# if defined key_roam
      fflx_c(:,:)  = 0.0        !! carbon     flux total
      fflx_a(:,:)  = 0.0        !! alkalinity flux total
      fflx_o2(:,:) = 0.0        !! oxygen     flux total
      ftot_c(:,:)  = 0.0        !! carbon     inventory
      ftot_a(:,:)  = 0.0        !! alkalinity inventory
      ftot_o2(:,:) = 0.0        !! oxygen     inventory
      fifd_c(:,:)  = 0.0        !! carbon     fast detritus production
! I don't think fifd_a or fifd_o2 are used - marc 11/5/17
!      fifd_a(:,:)  = 0.0        !! alkalinity fast detritus production
!      fifd_o2(:,:) = 0.0        !! oxygen     fast detritus production
      fofd_c(:,:)  = 0.0        !! carbon     fast detritus remineralisation
! I don't think fofd_a or fofd_o2 are used - marc 11/5/17
!      fofd_a(:,:)  = 0.0        !! alkalinity fast detritus remineralisation
!      fofd_o2(:,:) = 0.0        !! oxygen     fast detritus remineralisation
      !!
      fnit_prod(:,:) = 0.0      !! (organic)   nitrogen production
      fnit_cons(:,:) = 0.0      !! (organic)   nitrogen consumption
      fsil_prod(:,:) = 0.0      !! (inorganic) silicon production
      fsil_cons(:,:) = 0.0      !! (inorganic) silicon consumption
      fcar_prod(:,:) = 0.0      !! (organic)   carbon production
      fcar_cons(:,:) = 0.0      !! (organic)   carbon consumption
      !!
      foxy_prod(:,:) = 0.0      !! oxygen production
      foxy_cons(:,:) = 0.0      !! oxygen consumption
      foxy_anox(:,:) = 0.0      !! unrealised oxygen consumption
      !!
# endif
      ftot_n(:,:)   = 0.0       !! N inventory 
      ftot_si(:,:)  = 0.0       !! Si inventory
      ftot_fe(:,:)  = 0.0       !! Fe inventory
      ftot_pn(:,:)  = 0.0       !! integrated non-diatom phytoplankton
      ftot_pd(:,:)  = 0.0       !! integrated diatom     phytoplankton
      ftot_zmi(:,:) = 0.0       !! integrated microzooplankton
      ftot_zme(:,:) = 0.0       !! integrated mesozooplankton
      ftot_det(:,:) = 0.0       !! integrated slow detritus, nitrogen
      ftot_dtc(:,:) = 0.0       !! integrated slow detritus, carbon
      !!
      fzmi_i(:,:)  = 0.0        !! material grazed by microzooplankton
      fzmi_o(:,:)  = 0.0        !! ... sum of fate of this material
      fzme_i(:,:)  = 0.0        !! material grazed by  mesozooplankton
      fzme_o(:,:)  = 0.0        !! ... sum of fate of this material
      !!
      f_sbenin_n(:,:)  = 0.0    !! slow detritus N  -> benthic pool
      f_sbenin_fe(:,:) = 0.0    !! slow detritus Fe -> benthic pool
      f_sbenin_c(:,:)  = 0.0    !! slow detritus C  -> benthic pool
      f_fbenin_n(:,:)  = 0.0    !! fast detritus N  -> benthic pool
      f_fbenin_fe(:,:) = 0.0    !! fast detritus Fe -> benthic pool
      f_fbenin_si(:,:) = 0.0    !! fast detritus Si -> benthic pool
      f_fbenin_c(:,:)  = 0.0    !! fast detritus C  -> benthic pool
      f_fbenin_ca(:,:) = 0.0    !! fast detritus Ca -> benthic pool
      !!
      f_benout_n(:,:)  = 0.0    !! benthic N  pool  -> dissolved
      f_benout_fe(:,:) = 0.0    !! benthic Fe pool  -> dissolved
      f_benout_si(:,:) = 0.0    !! benthic Si pool  -> dissolved
      f_benout_c(:,:)  = 0.0    !! benthic C  pool  -> dissolved
      f_benout_ca(:,:) = 0.0    !! benthic Ca pool  -> dissolved
      !!
      f_benout_lyso_ca(:,:) = 0.0 !! benthic Ca pool  -> dissolved (when it shouldn't!)
      !!
      f_runoff(:,:)  = 0.0      !! riverine runoff
      f_riv_n(:,:)   = 0.0      !! riverine N   input 
      f_riv_si(:,:)  = 0.0      !! riverine Si  input 
      f_riv_c(:,:)   = 0.0      !! riverine C   input 
      f_riv_alk(:,:) = 0.0      !! riverine alk input 
      !! 
      !! Jpalm -- 06-03-2017 -- Forgotten var to init
      f_omarg(:,:) = 0.0        !!
      f_omcal(:,:) = 0.0 
      xFree(:,:) = 0.0          !! state variables for iron-ligand system
      fcomm_resp(:,:) = 0.0 
      fprn_ml(:,:) = 0.0        !! mixed layer PP diagnostics
      fprd_ml(:,:) = 0.0        !! mixed layer PP diagnostics
      !! AXY (16/08/17)
      fchl_ml(:,:) = 0.0	!! mixed layer chlorophyll diagnostics
      !! 
      fslownflux(:,:) = 0.0 
      fslowcflux(:,:) = 0.0 
      !!
      !! JPALM -- 21-09-2017 -- needed to debug air-sea carb
      f_xco2a(:,:)  = 0.0
      f_pco2w(:,:)  = 0.0
      f_ph(:,:)     = 0.0
      f_kw660(:,:)  = 0.0
      ztmp(:,:)  = 0.0
      zsal(:,:)  = 0.0
      zalk(:,:)  = 0.0
      zdic(:,:)  = 0.0
      zsil(:,:)  = 0.0
# if defined key_mocsy
      ! zpho is only defined if key_mocsy
      ! is active, so we must protect this
      ! initialisation accordingly. 
      zpho(:,:)  = 0.0
# endif
      f_co2flux(:,:)  = 0.0 
      f_pco2atm(:,:)  = 0.0
      f_h2co3(:,:)    = 0.0
      f_hco3(:,:)     = 0.0
      f_co3(:,:)      = 0.0
      f_omarg(:,:)    = 0.0
      f_omcal(:,:)    = 0.0
      !!
      !! AXY (08/08/17): zero slow detritus fluxes
      fslowsink(:,:)  = 0.0
# if defined key_roam
      fslowsinkc(:,:) = 0.0
# endif      
      !!
      !! allocate and initiate 2D diag
      !! -----------------------------
      !! Juju :: add kt condition !!
      IF ( lk_iomput ) THEN 

         !! initialise iom_use test
         IF ( kt == nittrc000 )   CALL trc_nam_iom_medusa 

         ALLOCATE( zw2d(1:jpi, 1:jpj) )
         zw2d(:,:)      = 0.0   !!
         IF ( med_diag%PRN%dgsave ) THEN
            ALLOCATE( fprn2d(1:jpi, 1:jpj) )
            fprn2d(:,:)      = 0.0 !!
         ENDIF
         IF ( med_diag%MPN%dgsave ) THEN
            ALLOCATE( fdpn2d(1:jpi, 1:jpj) )
            fdpn2d(:,:)      = 0.0 !!
         ENDIF
         IF ( med_diag%PRD%dgsave ) THEN
            ALLOCATE( fprd2d(1:jpi, 1:jpj) )
            fprd2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%MPD%dgsave ) THEN
            ALLOCATE( fdpd2d(1:jpi, 1:jpj) )
            fdpd2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OPAL%dgsave ) THEN
            ALLOCATE( fprds2d(1:jpi, 1:jpj) )
            fprds2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OPALDISS%dgsave ) THEN
            ALLOCATE( fsdiss2d(1:jpi, 1:jpj) )
            fsdiss2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMIPn%dgsave ) THEN
            ALLOCATE( fgmipn2d(1:jpi, 1:jpj) )
            fgmipn2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMID%dgsave ) THEN
            ALLOCATE( fgmid2d(1:jpi, 1:jpj) )
            fgmid2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%MZMI%dgsave ) THEN
            ALLOCATE( fdzmi2d(1:jpi, 1:jpj) )
            fdzmi2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMEPN%dgsave ) THEN
            ALLOCATE( fgmepn2d(1:jpi, 1:jpj) )
            fgmepn2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMEPD%dgsave ) THEN
            ALLOCATE( fgmepd2d(1:jpi, 1:jpj) )
            fgmepd2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMEZMI%dgsave ) THEN
            ALLOCATE( fgmezmi2d(1:jpi, 1:jpj) )
            fgmezmi2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMED%dgsave ) THEN
            ALLOCATE( fgmed2d(1:jpi, 1:jpj) )
            fgmed2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%MZME%dgsave ) THEN
            ALLOCATE( fdzme2d(1:jpi, 1:jpj) )
            fdzme2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%DETN%dgsave ) THEN
            ALLOCATE( fslown2d(1:jpi, 1:jpj) )
            fslown2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%MDET%dgsave ) THEN
            ALLOCATE( fdd2d(1:jpi, 1:jpj) )
            fdd2d(:,:)      = 0.0 !!
         ENDIF      
         IF( med_diag%AEOLIAN%dgsave ) THEN
            ALLOCATE( ffetop2d(1:jpi, 1:jpj) )
            ffetop2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%BENTHIC%dgsave ) THEN
            ALLOCATE( ffebot2d(1:jpi, 1:jpj) )
            ffebot2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%SCAVENGE%dgsave ) THEN
            ALLOCATE( ffescav2d(1:jpi, 1:jpj) )
            ffescav2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PN_JLIM%dgsave ) THEN
            ALLOCATE( fjln2d(1:jpi, 1:jpj) )
            fjln2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PN_NLIM%dgsave ) THEN
            ALLOCATE( fnln2d(1:jpi, 1:jpj) )
            fnln2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PN_FELIM%dgsave ) THEN
            ALLOCATE( ffln2d(1:jpi, 1:jpj) )
            ffln2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PD_JLIM%dgsave ) THEN
            ALLOCATE( fjld2d(1:jpi, 1:jpj) )
            fjld2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PD_NLIM%dgsave ) THEN
            ALLOCATE( fnld2d(1:jpi, 1:jpj) )
            fnld2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PD_FELIM%dgsave ) THEN
            ALLOCATE( ffld2d(1:jpi, 1:jpj) )
            ffld2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PD_SILIM%dgsave ) THEN
            ALLOCATE( fsld2d2(1:jpi, 1:jpj) )
            fsld2d2(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PDSILIM2%dgsave ) THEN
            ALLOCATE( fsld2d(1:jpi, 1:jpj) )
            fsld2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip SDT_XXXX diagnostics here
!!
         IF( med_diag%TOTREG_N%dgsave ) THEN
            ALLOCATE( fregen2d(1:jpi, 1:jpj) )
            fregen2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%TOTRG_SI%dgsave ) THEN
            ALLOCATE( fregensi2d(1:jpi, 1:jpj) )
            fregensi2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip REG_XXXX diagnostics here
!!
         IF( med_diag%FASTN%dgsave ) THEN
            ALLOCATE( ftempn2d(1:jpi, 1:jpj) )
            ftempn2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%FASTSI%dgsave ) THEN
            ALLOCATE( ftempsi2d(1:jpi, 1:jpj) )
            ftempsi2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%FASTFE%dgsave ) THEN
            ALLOCATE( ftempfe2d(1:jpi, 1:jpj) )
            ftempfe2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%FASTC%dgsave ) THEN
            ALLOCATE( ftempc2d(1:jpi, 1:jpj) )
            ftempc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%FASTCA%dgsave ) THEN
            ALLOCATE( ftempca2d(1:jpi, 1:jpj) )
            ftempca2d(:,:)      = 0.0 !!
         ENDIF     
!!
!! skip FDT_XXXX, RG_XXXXF, FDS_XXXX, RGS_XXXXF diagnostics here
!!
         IF( med_diag%REMINN%dgsave ) THEN
            ALLOCATE( freminn2d(1:jpi, 1:jpj) )
            freminn2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%REMINSI%dgsave ) THEN
            ALLOCATE( freminsi2d(1:jpi, 1:jpj) )
            freminsi2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%REMINFE%dgsave ) THEN
            ALLOCATE( freminfe2d(1:jpi, 1:jpj) )
            freminfe2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%REMINC%dgsave ) THEN
            ALLOCATE( freminc2d(1:jpi, 1:jpj) )
            freminc2d(:,:)      = 0.0 !! 
         ENDIF
         IF( med_diag%REMINCA%dgsave ) THEN
            ALLOCATE( freminca2d(1:jpi, 1:jpj) )
            freminca2d(:,:)      = 0.0 !!
         ENDIF

# if defined key_roam
!!
!! skip SEAFLRXX, MED_XXXX, INTFLX_XX, INT_XX, ML_XXX, OCAL_XXX, FE_XXXX, MED_XZE, WIND diagnostics here
!!
         IF( med_diag%RR_0100%dgsave ) THEN
            ALLOCATE( ffastca2d(1:jpi, 1:jpj) )
            ffastca2d(:,:)      = 0.0 !!
         ENDIF

         IF( med_diag%ATM_PCO2%dgsave ) THEN
            ALLOCATE( f_pco2a2d(1:jpi, 1:jpj) )
            f_pco2a2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip OCN_PH diagnostic here
!!
         IF( med_diag%OCN_PCO2%dgsave ) THEN
            ALLOCATE( f_pco2w2d(1:jpi, 1:jpj) )
            f_pco2w2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip OCNH2CO3, OCN_HCO3, OCN_CO3 diagnostics here
!!
         IF( med_diag%CO2FLUX%dgsave ) THEN
            ALLOCATE( f_co2flux2d(1:jpi, 1:jpj) )
            f_co2flux2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip OM_XXX diagnostics here
!!
         IF( med_diag%TCO2%dgsave ) THEN
            ALLOCATE( f_TDIC2d(1:jpi, 1:jpj) )
            f_TDIC2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%TALK%dgsave ) THEN
            ALLOCATE( f_TALK2d(1:jpi, 1:jpj) )
            f_TALK2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%KW660%dgsave ) THEN
            ALLOCATE( f_kw6602d(1:jpi, 1:jpj) )
            f_kw6602d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ATM_PP0%dgsave ) THEN
            ALLOCATE( f_pp02d(1:jpi, 1:jpj) )
            f_pp02d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%O2FLUX%dgsave ) THEN
            ALLOCATE( f_o2flux2d(1:jpi, 1:jpj) )
            f_o2flux2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%O2SAT%dgsave ) THEN
            ALLOCATE( f_o2sat2d(1:jpi, 1:jpj) )
            f_o2sat2d(:,:)      = 0.0 !!
         ENDIF 
!!
!! skip XXX_CCD diagnostics here
!! 
         IF( med_diag%SFR_OCAL%dgsave ) THEN
            ALLOCATE( sfr_ocal2d(1:jpi, 1:jpj) )
            sfr_ocal2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%SFR_OARG%dgsave ) THEN
            ALLOCATE( sfr_oarg2d(1:jpi, 1:jpj) )
            sfr_oarg2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip XX_PROD, XX_CONS, O2_ANOX, RR_XXXX diagnostics here
!! 
         IF( med_diag%IBEN_N%dgsave ) THEN
            ALLOCATE( iben_n2d(1:jpi, 1:jpj) )
            iben_n2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%IBEN_FE%dgsave ) THEN
            ALLOCATE( iben_fe2d(1:jpi, 1:jpj) )
            iben_fe2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%IBEN_C%dgsave ) THEN
            ALLOCATE( iben_c2d(1:jpi, 1:jpj) )
            iben_c2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%IBEN_SI%dgsave ) THEN
            ALLOCATE( iben_si2d(1:jpi, 1:jpj) )
            iben_si2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%IBEN_CA%dgsave ) THEN
            ALLOCATE( iben_ca2d(1:jpi, 1:jpj) )
            iben_ca2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OBEN_N%dgsave ) THEN
            ALLOCATE( oben_n2d(1:jpi, 1:jpj) )
            oben_n2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OBEN_FE%dgsave ) THEN
            ALLOCATE( oben_fe2d(1:jpi, 1:jpj) )
            oben_fe2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OBEN_C%dgsave ) THEN
            ALLOCATE( oben_c2d(1:jpi, 1:jpj) )
            oben_c2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OBEN_SI%dgsave ) THEN
            ALLOCATE( oben_si2d(1:jpi, 1:jpj) )
            oben_si2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OBEN_CA%dgsave ) THEN
            ALLOCATE( oben_ca2d(1:jpi, 1:jpj) )
            oben_ca2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip BEN_XX diagnostics here
!!

         IF( med_diag%RIV_N%dgsave ) THEN
            ALLOCATE( rivn2d(1:jpi, 1:jpj) )
            rivn2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%RIV_SI%dgsave ) THEN
            ALLOCATE( rivsi2d(1:jpi, 1:jpj) )
            rivsi2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%RIV_C%dgsave ) THEN
            ALLOCATE( rivc2d(1:jpi, 1:jpj) )
            rivc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%RIV_ALK%dgsave ) THEN
            ALLOCATE( rivalk2d(1:jpi, 1:jpj) )
            rivalk2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%DETC%dgsave ) THEN
            ALLOCATE( fslowc2d(1:jpi, 1:jpj) )
            fslowc2d(:,:)      = 0.0 !!
         ENDIF 
!!
!! skip SDC_XXXX, INVTXXX diagnostics here
!!
         IF( med_diag%LYSO_CA%dgsave ) THEN
            ALLOCATE( lyso_ca2d(1:jpi, 1:jpj) )
            lyso_ca2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip COM_RESP diagnostic here
!!
         IF( med_diag%PN_LLOSS%dgsave ) THEN
            ALLOCATE( fdpn22d(1:jpi, 1:jpj) )
            fdpn22d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%PD_LLOSS%dgsave ) THEN
            ALLOCATE( fdpd22d(1:jpi, 1:jpj) )
            fdpd22d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_LLOSS%dgsave ) THEN
            ALLOCATE( fdzmi22d(1:jpi, 1:jpj) )
            fdzmi22d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZE_LLOSS%dgsave ) THEN
            ALLOCATE( fdzme22d(1:jpi, 1:jpj) )
            fdzme22d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_MES_N%dgsave ) THEN   
            ALLOCATE( zimesn2d(1:jpi, 1:jpj) )
            zimesn2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_MES_D%dgsave ) THEN
            ALLOCATE( zimesd2d(1:jpi, 1:jpj) )
            zimesd2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_MES_C%dgsave ) THEN
            ALLOCATE( zimesc2d(1:jpi, 1:jpj) )
            zimesc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_MESDC%dgsave ) THEN
            ALLOCATE( zimesdc2d(1:jpi, 1:jpj) )
            zimesdc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_EXCR%dgsave ) THEN
            ALLOCATE( ziexcr2d(1:jpi, 1:jpj) )
            ziexcr2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_RESP%dgsave ) THEN
            ALLOCATE( ziresp2d(1:jpi, 1:jpj) )
            ziresp2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZI_GROW%dgsave ) THEN
            ALLOCATE( zigrow2d(1:jpi, 1:jpj) )
            zigrow2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZE_MES_N%dgsave ) THEN
            ALLOCATE( zemesn2d(1:jpi, 1:jpj) )
            zemesn2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZE_MES_D%dgsave ) THEN
            ALLOCATE( zemesd2d(1:jpi, 1:jpj) )
            zemesd2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZE_MES_C%dgsave ) THEN
            ALLOCATE( zemesc2d(1:jpi, 1:jpj) )
            zemesc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZE_MESDC%dgsave ) THEN
            ALLOCATE( zemesdc2d(1:jpi, 1:jpj) )
            zemesdc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZE_EXCR%dgsave ) THEN
            ALLOCATE( zeexcr2d(1:jpi, 1:jpj) )
            zeexcr2d(:,:)      = 0.0 !!
         ENDIF                  
         IF( med_diag%ZE_RESP%dgsave ) THEN
            ALLOCATE( zeresp2d(1:jpi, 1:jpj) )
            zeresp2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ZE_GROW%dgsave ) THEN
            ALLOCATE( zegrow2d(1:jpi, 1:jpj) )
            zegrow2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%MDETC%dgsave ) THEN
            ALLOCATE( mdetc2d(1:jpi, 1:jpj) )
            mdetc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMIDC%dgsave ) THEN
            ALLOCATE( gmidc2d(1:jpi, 1:jpj) )
            gmidc2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%GMEDC%dgsave ) THEN
            ALLOCATE( gmedc2d(1:jpi, 1:jpj) )
            gmedc2d(:,:)      = 0.0 !!
         ENDIF
!!
!! skip INT_XXX diagnostics here
!!
         IF (jdms .eq. 1) THEN
            IF( med_diag%DMS_SURF%dgsave ) THEN
               ALLOCATE( dms_surf2d(1:jpi, 1:jpj) )
               dms_surf2d(:,:)      = 0.0 !!
            ENDIF
            IF( med_diag%DMS_ANDR%dgsave ) THEN
               ALLOCATE( dms_andr2d(1:jpi, 1:jpj) )
               dms_andr2d(:,:)      = 0.0 !!
            ENDIF
            IF( med_diag%DMS_SIMO%dgsave ) THEN
               ALLOCATE( dms_simo2d(1:jpi, 1:jpj) )
               dms_simo2d(:,:)      = 0.0 !!
            ENDIF
            IF( med_diag%DMS_ARAN%dgsave ) THEN
               ALLOCATE( dms_aran2d(1:jpi, 1:jpj) )
               dms_aran2d(:,:)      = 0.0 !!
            ENDIF
            IF( med_diag%DMS_HALL%dgsave ) THEN
               ALLOCATE( dms_hall2d(1:jpi, 1:jpj) )
               dms_hall2d(:,:)      = 0.0 !!
            ENDIF
            IF( med_diag%DMS_ANDM%dgsave ) THEN
               ALLOCATE( dms_andm2d(1:jpi, 1:jpj) )
               dms_andm2d(:,:)      = 0.0 !!
            ENDIF
         ENDIF   
         !!
         !! AXY (24/11/16): extra MOCSY diagnostics, 2D
         IF( med_diag%ATM_XCO2%dgsave ) THEN
            ALLOCATE( f_xco2a_2d(1:jpi, 1:jpj) )
            f_xco2a_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OCN_FCO2%dgsave ) THEN
            ALLOCATE( f_fco2w_2d(1:jpi, 1:jpj) )
            f_fco2w_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%ATM_FCO2%dgsave ) THEN
            ALLOCATE( f_fco2a_2d(1:jpi, 1:jpj) )
            f_fco2a_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OCN_RHOSW%dgsave ) THEN
            ALLOCATE( f_ocnrhosw_2d(1:jpi, 1:jpj) )
            f_ocnrhosw_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OCN_SCHCO2%dgsave ) THEN
            ALLOCATE( f_ocnschco2_2d(1:jpi, 1:jpj) )
            f_ocnschco2_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OCN_KWCO2%dgsave ) THEN
            ALLOCATE( f_ocnkwco2_2d(1:jpi, 1:jpj) )
            f_ocnkwco2_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OCN_K0%dgsave ) THEN
            ALLOCATE( f_ocnk0_2d(1:jpi, 1:jpj) )
            f_ocnk0_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%CO2STARAIR%dgsave ) THEN
            ALLOCATE( f_co2starair_2d(1:jpi, 1:jpj) )
            f_co2starair_2d(:,:)      = 0.0 !!
         ENDIF
         IF( med_diag%OCN_DPCO2%dgsave ) THEN
            ALLOCATE( f_ocndpco2_2d(1:jpi, 1:jpj) )
            f_ocndpco2_2d(:,:)      = 0.0 !!
         ENDIF
# endif  
         IF( med_diag%TPP3%dgsave ) THEN
            ALLOCATE( tpp3d(1:jpi, 1:jpj, 1:jpk) )
            tpp3d(:,:,:)      = 0.0 !! 
         ENDIF
         IF( med_diag%DETFLUX3%dgsave ) THEN
            ALLOCATE( detflux3d(1:jpi, 1:jpj, 1:jpk) )
            detflux3d(:,:,:)      = 0.0 !! 
         ENDIF
         IF( med_diag%REMIN3N%dgsave ) THEN
             ALLOCATE( remin3dn(1:jpi, 1:jpj, 1:jpk) )
             remin3dn(:,:,:)      = 0.0 !! 
         ENDIF
         !! 
         !! AXY (10/11/16): CMIP6 diagnostics, 2D
         !! JPALM -- 17-11-16 -- put fgco2 alloc out of diag request
         !!                   needed for coupling/passed through restart
         !! IF( med_diag%FGCO2%dgsave ) THEN
            ALLOCATE( fgco2(1:jpi, 1:jpj) )
            fgco2(:,:)      = 0.0 !!
         !! ENDIF
         IF( med_diag%INTDISSIC%dgsave ) THEN
            ALLOCATE( intdissic(1:jpi, 1:jpj) )
            intdissic(:,:)  = 0.0 !!
         ENDIF          
         IF( med_diag%INTDISSIN%dgsave ) THEN
            ALLOCATE( intdissin(1:jpi, 1:jpj) )
            intdissin(:,:)  = 0.0 !!
         ENDIF          
         IF( med_diag%INTDISSISI%dgsave ) THEN
            ALLOCATE( intdissisi(1:jpi, 1:jpj) )
            intdissisi(:,:)  = 0.0 !!
         ENDIF          
         IF( med_diag%INTTALK%dgsave ) THEN
            ALLOCATE( inttalk(1:jpi, 1:jpj) )
            inttalk(:,:)  = 0.0 !!
         ENDIF          
         IF( med_diag%O2min%dgsave ) THEN
            ALLOCATE( o2min(1:jpi, 1:jpj) )
            o2min(:,:)  = 1.e3 !! set to high value as we're looking for min(o2)
         ENDIF          
         IF( med_diag%ZO2min%dgsave ) THEN
            ALLOCATE( zo2min(1:jpi, 1:jpj) )
            zo2min(:,:)  = 0.0 !!
         ENDIF          
         IF( med_diag%FBDDTALK%dgsave  ) THEN
            ALLOCATE( fbddtalk(1:jpi, 1:jpj) )
            fbddtalk(:,:)  = 0.0 !! 
         ENDIF
         IF( med_diag%FBDDTDIC%dgsave  ) THEN
            ALLOCATE( fbddtdic(1:jpi, 1:jpj) )
            fbddtdic(:,:)  = 0.0 !! 
         ENDIF
         IF( med_diag%FBDDTDIFE%dgsave ) THEN
            ALLOCATE( fbddtdife(1:jpi, 1:jpj) )
            fbddtdife(:,:) = 0.0 !! 
         ENDIF
         IF( med_diag%FBDDTDIN%dgsave  ) THEN
            ALLOCATE( fbddtdin(1:jpi, 1:jpj) )
            fbddtdin(:,:)  = 0.0 !! 
         ENDIF
         IF( med_diag%FBDDTDISI%dgsave ) THEN
            ALLOCATE( fbddtdisi(1:jpi, 1:jpj) )
            fbddtdisi(:,:) = 0.0 !! 
         ENDIF
         !! 
         !! AXY (10/11/16): CMIP6 diagnostics, 3D
         IF( med_diag%TPPD3%dgsave     ) THEN
            ALLOCATE( tppd3(1:jpi, 1:jpj, 1:jpk) )
            tppd3(:,:,:)     = 0.0 !! 
         ENDIF
         IF( med_diag%BDDTALK3%dgsave  ) THEN
            ALLOCATE( bddtalk3(1:jpi, 1:jpj, 1:jpk) )
            bddtalk3(:,:,:)  = 0.0 !! 
         ENDIF
         IF( med_diag%BDDTDIC3%dgsave  ) THEN
            ALLOCATE( bddtdic3(1:jpi, 1:jpj, 1:jpk) )
            bddtdic3(:,:,:)  = 0.0 !! 
         ENDIF
         IF( med_diag%BDDTDIFE3%dgsave ) THEN
            ALLOCATE( bddtdife3(1:jpi, 1:jpj, 1:jpk) )
            bddtdife3(:,:,:) = 0.0 !! 
         ENDIF
         IF( med_diag%BDDTDIN3%dgsave  ) THEN
            ALLOCATE( bddtdin3(1:jpi, 1:jpj, 1:jpk) )
            bddtdin3(:,:,:)  = 0.0 !! 
         ENDIF
         IF( med_diag%BDDTDISI3%dgsave ) THEN
            ALLOCATE( bddtdisi3(1:jpi, 1:jpj, 1:jpk) )
            bddtdisi3(:,:,:) = 0.0 !! 
         ENDIF
         IF( med_diag%FD_NIT3%dgsave   ) THEN
            ALLOCATE( fd_nit3(1:jpi, 1:jpj, 1:jpk) )
            fd_nit3(:,:,:)   = 0.0 !! 
         ENDIF
         IF( med_diag%FD_SIL3%dgsave   ) THEN
            ALLOCATE( fd_sil3(1:jpi, 1:jpj, 1:jpk) )
            fd_sil3(:,:,:)   = 0.0 !! 
         ENDIF
         IF( med_diag%FD_CAR3%dgsave   ) THEN
            ALLOCATE( fd_car3(1:jpi, 1:jpj, 1:jpk) )
            fd_car3(:,:,:)   = 0.0 !! 
         ENDIF
         IF( med_diag%FD_CAL3%dgsave   ) THEN
            ALLOCATE( fd_cal3(1:jpi, 1:jpj, 1:jpk) )
            fd_cal3(:,:,:)   = 0.0 !! 
         ENDIF
         IF( med_diag%DCALC3%dgsave    ) THEN
            ALLOCATE( dcalc3(1:jpi, 1:jpj, 1:jpk) )
            dcalc3(:,:,: )   = 0.0 !! 
         ENDIF
         IF( med_diag%EXPC3%dgsave     ) THEN
            ALLOCATE( expc3(1:jpi, 1:jpj, 1:jpk) )
            expc3(:,:,: )    = 0.0 !! 
         ENDIF
         IF( med_diag%EXPN3%dgsave     ) THEN
            ALLOCATE( expn3(1:jpi, 1:jpj, 1:jpk) )
            expn3(:,:,: )    = 0.0 !! 
         ENDIF
         IF( med_diag%FEDISS3%dgsave   ) THEN
            ALLOCATE( fediss3(1:jpi, 1:jpj, 1:jpk) )
            fediss3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%FESCAV3%dgsave   ) THEN
            ALLOCATE( fescav3(1:jpi, 1:jpj, 1:jpk) )
            fescav3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%MIGRAZP3%dgsave   ) THEN
            ALLOCATE( migrazp3(1:jpi, 1:jpj, 1:jpk) )
            migrazp3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%MIGRAZD3%dgsave   ) THEN
            ALLOCATE( migrazd3(1:jpi, 1:jpj, 1:jpk) )
            migrazd3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%MEGRAZP3%dgsave   ) THEN
            ALLOCATE( megrazp3(1:jpi, 1:jpj, 1:jpk) )
            megrazp3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%MEGRAZD3%dgsave   ) THEN
            ALLOCATE( megrazd3(1:jpi, 1:jpj, 1:jpk) )
            megrazd3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%MEGRAZZ3%dgsave   ) THEN
            ALLOCATE( megrazz3(1:jpi, 1:jpj, 1:jpk) )
            megrazz3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%O2SAT3%dgsave     ) THEN
            ALLOCATE( o2sat3(1:jpi, 1:jpj, 1:jpk) )
            o2sat3(:,:,: )    = 0.0 !! 
         ENDIF
         IF( med_diag%PBSI3%dgsave      ) THEN
            ALLOCATE( pbsi3(1:jpi, 1:jpj, 1:jpk) )
            pbsi3(:,:,: )     = 0.0 !! 
         ENDIF
         IF( med_diag%PCAL3%dgsave      ) THEN
            ALLOCATE( pcal3(1:jpi, 1:jpj, 1:jpk) )
            pcal3(:,:,: )     = 0.0 !! 
         ENDIF
         IF( med_diag%REMOC3%dgsave     ) THEN
            ALLOCATE( remoc3(1:jpi, 1:jpj, 1:jpk) )
            remoc3(:,:,: )    = 0.0 !! 
         ENDIF
         IF( med_diag%PNLIMJ3%dgsave    ) THEN
            ALLOCATE( pnlimj3(1:jpi, 1:jpj, 1:jpk) )
            pnlimj3(:,:,: )   = 0.0 !! 
         ENDIF
         IF( med_diag%PNLIMN3%dgsave    ) THEN
            ALLOCATE( pnlimn3(1:jpi, 1:jpj, 1:jpk) )
            pnlimn3(:,:,: )   = 0.0 !! 
         ENDIF
         IF( med_diag%PNLIMFE3%dgsave   ) THEN
            ALLOCATE( pnlimfe3(1:jpi, 1:jpj, 1:jpk) )
            pnlimfe3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%PDLIMJ3%dgsave    ) THEN
            ALLOCATE( pdlimj3(1:jpi, 1:jpj, 1:jpk) )
            pdlimj3(:,:,: )   = 0.0 !! 
         ENDIF
         IF( med_diag%PDLIMN3%dgsave    ) THEN
            ALLOCATE( pdlimn3(1:jpi, 1:jpj, 1:jpk) )
            pdlimn3(:,:,: )   = 0.0 !! 
         ENDIF
         IF( med_diag%PDLIMFE3%dgsave   ) THEN
            ALLOCATE( pdlimfe3(1:jpi, 1:jpj, 1:jpk) )
            pdlimfe3(:,:,: )  = 0.0 !! 
         ENDIF
         IF( med_diag%PDLIMSI3%dgsave   ) THEN
            ALLOCATE( pdlimsi3(1:jpi, 1:jpj, 1:jpk) )
            pdlimsi3(:,:,: )  = 0.0 !! 
         ENDIF
      ENDIF
      !! lk_iomput

   END SUBROUTINE bio_medusa_init
   !CEB !$AGRIF_END_DO_NOT_TREAT
#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE bio_medusa_init( )                   ! Empty routine
      IMPLICIT NONE
      WRITE(*,*) 'bio_medusa_init: You should not have seen this print! error?'
   END SUBROUTINE bio_medusa_init
#endif 

   !!======================================================================
END MODULE bio_medusa_init_mod
