MODULE bio_medusa_mod
   !!======================================================================
   !!                      ***  MODULE  bio_medusa_mod  ***
   !! MEDUSA variables   :  module for MEDUSA variables which are shared
   !!                       across subroutines 
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!   -   ! 2017-08 (A. Yool)            Slow detritus, ML-avg chl variables
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'                                           MEDUSA
   !!----------------------------------------------------------------------
   !! Variable conventions
   !!----------------------------------------------------------------------
   !!
   !! names: z*** - state variable 
   !!        f*** - function (or temporary variable used in part of a function)
   !!        b*** - right-hand part (sources and sinks)
   !!        i*** - integer variable (usually used in yes/no flags)
   !!----------------------------------------------------------------------
   USE par_kind,          ONLY: wp
   !CEB
   USE trc 
   IMPLICIT NONE
   PUBLIC

   PUBLIC   bio_medusa_alloc     ! called by trcini.F90

   !! model state variables
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zchn,zchd,zphn,zphd,zpds,zzmi
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zzme,zdet,zdtc,zdin,zsil,zfer
# if defined key_roam
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zdic, zalk, zoxy
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztmp, zsal
# endif
# if defined key_mocsy
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zpho
# endif

   !! integrated source and sink terms
   REAL(wp) ::    b0

   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fthetan,fprn,frn
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fthetad,fprd,frd

   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fjlim_pn,fjlim_pd
   !! AXY (16/07/09): add in Eppley curve functionality
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fun_T,xvpnT,xvpdT

   !! AXY (16/05/11): per Katya's prompting, add in new T-dependence
   !!                 for phytoplankton growth only (i.e. no change
   !!                 for remineralisation)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fun_Q10
   !! AXY (01/03/10): add in mixed layer PP diagnostics
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fprn_ml,fprd_ml
   !! AXY (16/08/17): add in mixed layer chlorophyll diagnostic
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fchl_ml
   !!
   !! nutrient limiting factors
   !! N and Fe (renaming ffln to ffln2 to avoid conflict with
   !! ffln in module sms_medusa - marc 25/4/17)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fnln,ffln2
   !! N, Fe and Si
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fnld,ffld,fsld,fsld2
   !!
   !! silicon cycle
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fsin,fprds,fsdiss

   !! iron cycle; includes parameters for Parekh et al. (2005) iron scheme
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ffetop,ffebot,ffescav
   !! Variable for iron-ligand system
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: xFree

   !! Microzooplankton grazing
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fgmipn,fgmid,fgmidc
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: finmi,ficmi
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fmigrow,fmiexcr,fmiresp
   !!
   !! Mesozooplankton grazing
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fgmepn,fgmepd
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fgmepds,fgmezmi,fgmed,fgmedc
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: finme,ficme
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fmegrow,fmeexcr,fmeresp
   !!
   !! mortality/Remineralisation
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fdpn,fdpd,fdpds,fdzmi,fdzme,fdd
# if defined key_roam
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fddc
# endif
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fdpn2,fdpd2,fdpds2,fdzmi2,fdzme2
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fslown, fslowc
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fslownflux, fslowcflux
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fregen,fregensi
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fregenfast,fregenfastsi
# if defined key_roam
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fregenfastc
# endif
   !!
   !! AXY (08/08/17): sinking of detritus moved here
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::    fslowsink, fslowgain, fslowloss
# if defined key_roam
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::    fslowsinkc, fslowgainc, fslowlossc
# endif
   !!
   !! Particle flux
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fdep1
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftempn,ftempsi,ftempfe
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftempc,ftempca
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: freminn,freminsi,freminfe
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: freminc,freminca
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ffastn,ffastsi,ffastfe
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ffastc,ffastca
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fsedn,fsedsi,fsedfe,fsedc,fsedca
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fccd

   !! AXY (08/07/11): fate of fast detritus reaching the seafloor
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ffast2slown,ffast2slowc

   !! water column nutrient and flux integrals
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftot_n,ftot_si,ftot_fe
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fflx_n,fflx_si,fflx_fe
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fifd_n,fifd_si,fifd_fe
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fofd_n,fofd_si,fofd_fe
# if defined key_roam
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftot_c,ftot_a,ftot_o2
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fflx_c,fflx_a,fflx_o2
! I don't think fifd_a, fifd_o2, fofd_a or fofd_o2 are used - marc 11/5/17
!   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fifd_c,fifd_a,fifd_o2
!   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fofd_c,fofd_a,fofd_o2
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fifd_c
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fofd_c
# endif

   !! Zooplankton grazing integrals
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fzmi_i,fzmi_o,fzme_i,fzme_o

   !! Limitation term temporary variables
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftot_pn,ftot_pd
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftot_zmi,ftot_zme,ftot_det,ftot_dtc

   !! use biological fluxes (1) or not (0)
   INTEGER  ::    ibio_switch
   !!
   !! diagnose fluxes (should only be used in 1D runs)
   INTEGER                               :: idf, idfval

   !! Nitrogen and silicon production and consumption
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fnit_prod,fnit_cons
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fsil_prod,fsil_cons

# if defined key_roam
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_xco2a
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_ph,f_pco2w,f_h2co3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_hco3,f_co3,f_co2flux
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_TDIC,f_TALK,f_dcf,f_henry
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_pp0
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_kw660,f_o2flux,f_o2sat
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_omcal,f_omarg

   !! AXY (23/06/15): additional diagnostics for MOCSY and oxygen
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_pco2atm

   !! Carbon, alkalinity production and consumption
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fcomm_resp
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fcar_prod,fcar_cons

   !! Oxygen production and consumption (and non-consumption)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: foxy_prod,foxy_cons,foxy_anox

   !! Add DMS in MEDUSA for UKESM1 model
   REAL(wp)                              :: dms_surf,dms_andm
   !! AXY (13/03/15): add in other DMS calculations
   REAL(wp)                              :: dms_andr,dms_simo,dms_aran,dms_hall
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: dms_nlim, dms_wtkn
# endif

   !! Benthic fluxes
   INTEGER  ::    ibenthic
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_sbenin_n,f_sbenin_fe,f_sbenin_c
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_fbenin_n,f_fbenin_fe,f_fbenin_si
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_fbenin_c,f_fbenin_ca
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_benout_n,f_benout_fe,f_benout_si 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_benout_c,f_benout_ca

   !! Benthic fluxes of CaCO3 that shouldn't happen because of lysocline
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_benout_lyso_ca

   !! riverine fluxes
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_runoff,f_riv_n,f_riv_si
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_riv_c,f_riv_alk
   !! AXY (19/07/12): variables for local riverine fluxes to handle 
   !! inputs below surface
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_riv_loc_n,f_riv_loc_si
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_riv_loc_c, f_riv_loc_alk

   !! Jpalm -- 11-10-2015 -- adapt diag to iom_use
   !! 2D var for diagnostics.
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fprn2d, fdpn2d, fprd2d, fdpd2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fprds2d, fsdiss2d, fgmipn2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fgmid2d, fdzmi2d, fgmepn2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fgmepd2d, fgmezmi2d, fgmed2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fdzme2d, fslown2d, fdd2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ffetop2d, ffebot2d, ffescav2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fjln2d, fnln2d, ffln2d, fjld2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fnld2d, ffld2d, fsld2d2
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fsld2d, fregen2d, fregensi2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftempn2d, ftempsi2d, ftempfe2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ftempc2d, ftempca2d, freminn2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: freminsi2d, freminfe2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: freminc2d, freminca2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
# if defined key_roam
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ffastca2d, rivn2d, rivsi2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: rivc2d, rivalk2d, fslowc2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fdpn22d, fdpd22d, fdzmi22d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fdzme22d, zimesn2d, zimesd2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zimesc2d, zimesdc2d, ziexcr2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ziresp2d, zigrow2d, zemesn2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zemesd2d, zemesc2d, zemesdc2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zeexcr2d, zeresp2d, zegrow2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: mdetc2d, gmidc2d, gmedc2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_pco2a2d, f_pco2w2d, f_co2flux2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_TDIC2d, f_TALK2d, f_kw6602d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_pp02d, f_o2flux2d, f_o2sat2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: dms_andr2d, dms_simo2d, dms_aran2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: dms_hall2d, dms_andm2d, dms_surf2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: iben_n2d, iben_fe2d, iben_c2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: iben_si2d, iben_ca2d, oben_n2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: oben_fe2d, oben_c2d, oben_si2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: oben_ca2d, sfr_ocal2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: sfr_oarg2d, lyso_ca2d 
   !! AXY (23/11/16): extra MOCSY diagnostics
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_xco2a_2d, f_fco2w_2d, f_fco2a_2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_ocnrhosw_2d, f_ocnschco2_2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_ocnkwco2_2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_ocnk0_2d, f_co2starair_2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_ocndpco2_2d
# endif
   !!
   !! 3D var for diagnostics.
   !! CEB REAL(wp), POINTER, DIMENSION(:,:,:) :: tpp3d, detflux3d, remin3dn
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: tpp3d, detflux3d, remin3dn
   !!/CEB
# if defined key_roam
   !! AXY (04/11/16)
   !! 2D var for new CMIP6 diagnostics (behind a key_roam ifdef 
   !! for simplicity)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fgco2,intdissic,intdissin
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: intdissisi,inttalk,o2min,zo2min
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fbddtalk,fbddtdic,fbddtdife
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fbddtdin,fbddtdisi
   !!
   !! 3D var for new CMIP6 diagnostics
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: tppd3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: bddtalk3,bddtdic3,bddtdife3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: bddtdin3, bddtdisi3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: fd_nit3,fd_sil3,fd_car3,fd_cal3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: co33,co3satarag3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: co3satcalc3,dcalc3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: expc3,expn3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: fediss3,fescav3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: migrazp3,migrazd3,megrazp3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: megrazd3, megrazz3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: o2sat3,pbsi3,pcal3,remoc3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: pnlimj3,pnlimn3,pnlimfe3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: pdlimj3,pdlimn3,pdlimfe3,pdlimsi3
# endif
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3.1 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   INTEGER FUNCTION bio_medusa_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE bio_medusa_alloc ***
      !!-------------------------------------------------------------------
      USE lib_mpp,           ONLY: ctl_warn
      USE par_oce,           ONLY: jpi, jpj
      !!-------------------------------------------------------------------
      !CEB change allocation that it pases through AGRIF
      ALLOCATE(zchn(jpi,jpj),zchd(jpi,jpj),zphn(jpi,jpj),             &
               zphd(jpi,jpj),zpds(jpi,jpj),zzmi(jpi,jpj),             &
               zzme(jpi,jpj),zdet(jpi,jpj),zdtc(jpi,jpj),             &
               zdin(jpi,jpj),zsil(jpi,jpj),zfer(jpi,jpj))             
# if defined key_roam
     ALLOCATE( zdic(jpi,jpj),zalk(jpi,jpj),zoxy(jpi,jpj),             &
               ztmp(jpi,jpj),zsal(jpi,jpj))                           
# endif
# if defined key_mocsy
      ALLOCATE(zpho(jpi,jpj))                                         
# endif
      ALLOCATE(fthetan(jpi,jpj),fprn(jpi,jpj),frn(jpi,jpj),           &
               fthetad(jpi,jpj),fprd(jpi,jpj),frd(jpi,jpj),           &
               fjlim_pn(jpi,jpj),fjlim_pd(jpi,jpj),                   &
               fun_T(jpi,jpj),fun_Q10(jpi,jpj),                       &
               fprn_ml(jpi,jpj),fprd_ml(jpi,jpj),fchl_ml(jpi,jpj))    
      ALLOCATE(fnln(jpi,jpj),ffln2(jpi,jpj),                          &
               fnld(jpi,jpj),ffld(jpi,jpj),fsld(jpi,jpj),             &
               fsld2(jpi,jpj),                                        &
               fsin(jpi,jpj),fprds(jpi,jpj),fsdiss(jpi,jpj),          &
               ffetop(jpi,jpj),ffebot(jpi,jpj),ffescav(jpi,jpj),      &
               xFree(jpi,jpj))                                        
      ALLOCATE(fgmipn(jpi,jpj),fgmid(jpi,jpj),fgmidc(jpi,jpj),        &
               finmi(jpi,jpj),ficmi(jpi,jpj),                         &
               fmigrow(jpi,jpj),fmiexcr(jpi,jpj),fmiresp(jpi,jpj),    &
               fgmepn(jpi,jpj),fgmepd(jpi,jpj),                       &
               fgmepds(jpi,jpj),fgmezmi(jpi,jpj),fgmed(jpi,jpj),      &
               fgmedc(jpi,jpj))                                       
      ALLOCATE(finme(jpi,jpj),ficme(jpi,jpj),                         &
               fmegrow(jpi,jpj),fmeexcr(jpi,jpj),fmeresp(jpi,jpj),    &
               fdpn(jpi,jpj),fdpd(jpi,jpj),fdpds(jpi,jpj),            &
               fdzmi(jpi,jpj),fdzme(jpi,jpj),fdd(jpi,jpj))            
# if defined key_roam
         ALLOCATE(fddc(jpi,jpj)) 
# endif
         ALLOCATE(fdpn2(jpi,jpj),fdpd2(jpi,jpj),fdpds2(jpi,jpj),         &
               fdzmi2(jpi,jpj),fdzme2(jpi,jpj),                       &
               fslown(jpi,jpj),fslowc(jpi,jpj),                       &
               fslownflux(jpi,jpj),fslowcflux(jpi,jpj),               &
               fregen(jpi,jpj),fregensi(jpi,jpj),                     &
               fregenfast(jpi,jpj),fregenfastsi(jpi,jpj))             
# if defined key_roam
       ALLOCATE(fregenfastc(jpi,jpj)) 
# endif
       ALLOCATE(fslowsink(jpi,jpj),fslowgain(jpi,jpj),                 &
               fslowloss(jpi,jpj))                                    
# if defined key_roam
       ALLOCATE(fslowsinkc(jpi,jpj),fslowgainc(jpi,jpj),               &
               fslowlossc(jpi,jpj))                                   
# endif
      ALLOCATE(fdep1(jpi,jpj),                                        &
               ftempn(jpi,jpj),ftempsi(jpi,jpj),ftempfe(jpi,jpj),     &
               ftempc(jpi,jpj),ftempca(jpi,jpj),                      &
               freminn(jpi,jpj),freminsi(jpi,jpj),freminfe(jpi,jpj),  &
               freminc(jpi,jpj),freminca(jpi,jpj),                    &
               ffastn(jpi,jpj),ffastsi(jpi,jpj),ffastfe(jpi,jpj),     &
               ffastc(jpi,jpj),ffastca(jpi,jpj))                      
      ALLOCATE(fsedn(jpi,jpj),fsedsi(jpi,jpj),fsedfe(jpi,jpj),        &
               fsedc(jpi,jpj),fsedca(jpi,jpj),                        &
               fccd(jpi,jpj),                                         &
               ffast2slown(jpi,jpj),ffast2slowc(jpi,jpj),             &
               ftot_n(jpi,jpj),ftot_si(jpi,jpj),ftot_fe(jpi,jpj),     &
               fflx_n(jpi,jpj),fflx_si(jpi,jpj),fflx_fe(jpi,jpj),     &
               fifd_n(jpi,jpj),fifd_si(jpi,jpj),fifd_fe(jpi,jpj),     &
               fofd_n(jpi,jpj),fofd_si(jpi,jpj),fofd_fe(jpi,jpj))     
# if defined key_roam
      ALLOCATE(ftot_c(jpi,jpj),ftot_a(jpi,jpj),ftot_o2(jpi,jpj),      &
               fflx_c(jpi,jpj),fflx_a(jpi,jpj),fflx_o2(jpi,jpj),      &
               fifd_c(jpi,jpj), fofd_c(jpi,jpj))                      
# endif
      ALLOCATE(fzmi_i(jpi,jpj),fzmi_o(jpi,jpj),fzme_i(jpi,jpj),       &
               fzme_o(jpi,jpj),                                       &
               ftot_pn(jpi,jpj),ftot_pd(jpi,jpj),                     &
               ftot_zmi(jpi,jpj),ftot_zme(jpi,jpj),ftot_det(jpi,jpj), &
               ftot_dtc(jpi,jpj),                                     &
               fnit_prod(jpi,jpj),fnit_cons(jpi,jpj),                 &
               fsil_prod(jpi,jpj),fsil_cons(jpi,jpj))                 
# if defined key_roam
      ALLOCATE(f_xco2a(jpi,jpj),                                      &
               f_ph(jpi,jpj),f_pco2w(jpi,jpj),f_h2co3(jpi,jpj),       &
               f_hco3(jpi,jpj),f_co3(jpi,jpj),f_co2flux(jpi,jpj),     &
               f_TDIC(jpi,jpj),f_TALK(jpi,jpj),f_dcf(jpi,jpj),        &
               f_henry(jpi,jpj),                                      &
               f_pp0(jpi,jpj))                                        
      ALLOCATE(f_kw660(jpi,jpj),f_o2flux(jpi,jpj),f_o2sat(jpi,jpj),   &
               f_omcal(jpi,jpj),f_omarg(jpi,jpj),                     &
               f_pco2atm(jpi,jpj),                 &
               fcomm_resp(jpi,jpj),                                   &
               fcar_prod(jpi,jpj),fcar_cons(jpi,jpj),                 &
               foxy_prod(jpi,jpj), foxy_cons(jpi,jpj),                &
               foxy_anox(jpi,jpj),                                    &
               dms_nlim(jpi,jpj),dms_wtkn(jpi,jpj))                   
# endif
      ALLOCATE(f_sbenin_n(jpi,jpj),f_sbenin_fe(jpi,jpj),              &
               f_sbenin_c(jpi,jpj),                                   &
               f_fbenin_n(jpi,jpj),f_fbenin_fe(jpi,jpj),              &
               f_fbenin_si(jpi,jpj),                                  &
               f_fbenin_c(jpi,jpj),f_fbenin_ca(jpi,jpj),              &
               f_benout_n(jpi,jpj),f_benout_fe(jpi,jpj),              &
               f_benout_si(jpi,jpj))                                  
      ALLOCATE(f_benout_c(jpi,jpj),f_benout_ca(jpi,jpj),              &
               f_benout_lyso_ca(jpi,jpj),                             &
               f_runoff(jpi,jpj),f_riv_n(jpi,jpj),f_riv_si(jpi,jpj),  &
               f_riv_c(jpi,jpj),f_riv_alk(jpi,jpj),                   &
               f_riv_loc_n(jpi,jpj),f_riv_loc_si(jpi,jpj),            &
               f_riv_loc_c(jpi,jpj),f_riv_loc_alk(jpi,jpj),           &
               STAT = bio_medusa_alloc)

      !! Check that allocation was successful
      IF ( bio_medusa_alloc /= 0 ) THEN  
         CALL ctl_warn('bio_medusa_alloc: failed to allocate arrays')
      END IF
   
   END FUNCTION bio_medusa_alloc
#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                          No MEDUSA
   !!----------------------------------------------------------------------
#endif 

   !!======================================================================
END MODULE bio_medusa_mod
