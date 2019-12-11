MODULE air_sea_mod
   !!======================================================================
   !!                         ***  MODULE air_sea_mod  ***
   !! Calculate the carbon chemistry for the whole ocean
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!   -   ! 2017-08 (A. Yool)            Add air-sea flux kill switch
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   air_sea        ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE air_sea( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE air_sea  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - calculate air-sea gas exchange
      !!  - river inputs
      !!----------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: f_riv_alk, f_riv_c, f_riv_n,           &
                                   f_riv_si, f_runoff,                    & 
                                   fgco2, zphn, zphd,                     &
# if defined key_roam
                                   dms_andr, dms_andr2d, dms_aran,        &
                                   dms_aran2d, dms_hall, dms_hall2d,      &
                                   dms_simo, dms_simo2d, dms_surf,        &
                                   dms_surf2d, dms_andm, dms_andm2d,      &
                                   dms_nlim, dms_wtkn,                    &
                                   f_co2flux, f_co2flux2d,                &
                                   f_co2starair_2d, f_co3,                &
                                   f_dcf, f_fco2a_2d, f_fco2w_2d,         &
                                   f_h2co3, f_hco3, f_henry,              &
                                   f_kw660, f_kw6602d,                    &
                                   f_o2flux, f_o2flux2d, f_o2sat,         &
                                   f_o2sat2d, f_ocndpco2_2d,              &
                                   f_ocnk0_2d, f_ocnkwco2_2d,             &
                                   f_ocnrhosw_2d, f_ocnschco2_2d,         &
                                   f_omarg, f_omcal,                      &
                                   f_pco2a2d, f_pco2atm, f_pco2w,         &
                                   f_pco2w2d, f_ph, f_pp0, f_pp02d,       &
                                   f_TALK, f_TALK2d, f_TDIC, f_TDIC2d,    &
                                   f_xco2a, f_xco2a_2d,                   &
                                   zalk, zdic, zoxy, zsal, ztmp,          &
# endif
# if defined key_mocsy
                                   zpho,                                  &
# endif
                                   zchd, zchn, zdin, zsil
      USE dom_oce,           ONLY: e3t_0, gphit, tmask, mig, mjg
# if defined key_vvl
      USE dom_oce,           ONLY: e3t_n
# endif
      USE iom,               ONLY: lk_iomput
      USE in_out_manager,    ONLY: lwp, numout
      USE par_kind,          ONLY: wp
      USE par_oce,           ONLY: jpi, jpim1, jpj, jpjm1
      USE sbc_oce,           ONLY: fr_i, qsr, wndm
      USE sms_medusa,        ONLY: jdms, jdms_input, jdms_model,          &
                                   jriver_alk, jriver_c,                  &
                                   jriver_n, jriver_si,                   &
                                   riv_alk, riv_c, riv_n, riv_si,         &
                                   zn_dms_chd, zn_dms_chn, zn_dms_din,    &
                                   zn_dms_mld, zn_dms_qsr,                &
                                   xnln, xnld 
      !CEB USE trc,               ONLY: med_diag
      USE trc
      USE zdfmxl,            ONLY: hmld

# if defined key_roam
      !CEB USE gas_transfer_mod,       ONLY: gas_transfer
#  if defined key_mocsy
      USE mocsy_wrapper,     ONLY: mocsy_interface
#  else
      USE trcco2_medusa,     ONLY: trc_co2_medusa
#  endif
      USE trcdms_medusa,     ONLY: trc_dms_medusa
      USE trcoxy_medusa,     ONLY: trc_oxy_medusa
# endif
      USE lib_mpp,           ONLY: ctl_stop
      USE trcstat,           ONLY: trc_rst_dia_stat 

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! time (integer timestep)
      INTEGER, INTENT( in ) :: kt

      !! Loop variables
      INTEGER :: ji, jj

# if defined key_roam
      !! jpalm 14-07-2016: convert CO2flux diag from mmol/m2/d to kg/m2/s
      REAL, PARAMETER :: weight_CO2_mol = 44.0095  !! g / mol
      REAL, PARAMETER :: secs_in_day    = 86400.0  !! s / d
      REAL, PARAMETER :: CO2flux_conv   = (1.e-6 * weight_CO2_mol) / secs_in_day

      INTEGER :: iters

      !! AXY (23/06/15): additional diagnostics for MOCSY and oxygen
      REAL(wp), DIMENSION(jpi,jpj) :: f_fco2w, f_rhosw
      REAL(wp), DIMENSION(jpi,jpj) :: f_fco2atm
      REAL(wp), DIMENSION(jpi,jpj) :: f_schmidtco2, f_kwco2, f_K0
      REAL(wp), DIMENSION(jpi,jpj) :: f_co2starair, f_dpco2
      !! Output arguments from mocsy_interface, which aren't used
      REAL(wp) :: f_BetaD_dum, f_opres_dum
      REAL(wp) :: f_insitut_dum
      REAL(wp) :: f_kwo2_dum
# endif


# if defined key_roam
      !! init
      f_fco2w(:,:)       = 0.0
      f_fco2atm(:,:)     = 0.0
      f_schmidtco2(:,:)  = 0.0
      f_kwco2(:,:)       = 0.0
      f_co2starair(:,:)  = 0.0
      f_dpco2(:,:)       = 0.0
      f_rhosw(:,:)       = 0.0
      f_K0(:,:)          = 0.0
      !! air pressure (atm); ultimately this will use air 
      !! pressure at the base of the UKESM1 atmosphere 
      !!                                     
      f_pp0(:,:)   = 1.0


      !!-----------------------------------------------------------
      !! Air-sea gas exchange
      !!-----------------------------------------------------------

#   if defined key_debug_medusa
               IF (lwp) write (numout,*)                     & 
               'air-sea: gas_transfer kt = ', kt
               CALL flush(numout)
#   endif
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            IF (tmask(ji,jj,1) == 1) then
               !!
               !! AXY (23/06/15): as part of an effort to update the 
               !!                 carbonate chemistry in MEDUSA, the gas 
               !!                 transfer velocity used in the carbon
               !!                 and oxygen cycles has been harmonised 
               !!                 and is calculated by the same function 
               !!                 here; this harmonisation includes
               !!                 changes to the PML carbonate chemistry 
               !!                 scheme so that it too makes use of the 
               !!                 same gas transfer velocity; the
               !!                 preferred parameterisation of this is 
               !!                 Wanninkhof (2014), option 7
               !!
   !CEB               CALL gas_transfer( wndm(ji,jj), 1, 7,         &  ! inputs
   !CEB                              f_kw660(ji,jj) )              ! outputs
   !CEB see Function below
                   f_kw660(ji,jj) =((0.251*wndm(ji,jj)**2) + (0.0*wndm(ji,jj)))/(100. * 3600.)
            ENDIF
         ENDDO
      ENDDO

#   if defined key_debug_medusa
               IF (lwp) write (numout,*)                     &
               'air-sea: carb-chem kt = ', kt
               CALL flush(numout)
               !! JPALM add carb print:
               call trc_rst_dia_stat(f_xco2a(:,:), 'f_xco2a')
               call trc_rst_dia_stat(wndm(:,:), 'wndm')
               call trc_rst_dia_stat(f_kw660(:,:), 'f_kw660')
               call trc_rst_dia_stat(ztmp(:,:), 'ztmp')
               call trc_rst_dia_stat(zsal(:,:), 'zsal')
               call trc_rst_dia_stat(zalk(:,:), 'zalk')
               call trc_rst_dia_stat(zdic(:,:), 'zdic')
               call trc_rst_dia_stat(zsil(:,:), 'zsil')
               call trc_rst_dia_stat(zpho(:,:), 'zpho')
#   endif
#  if defined key_axy_carbchem
#   if defined key_mocsy
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,1) == 1) then
               !!
               !! Jpalm -- 12-09-2017 -- add extra check after reccurent
               !!          carbonate failure in the coupled run.
               !!          must be associated to air-sea flux or air xCO2...
               !!          Check MOCSY inputs
               IF ( (zsal(ji,jj) > 75.0 ).OR.(zsal(ji,jj) < 0.0 ) .OR.        &
                    (ztmp(ji,jj) > 50.0 ).OR.(ztmp(ji,jj) < -20.0 ) .OR.      &
                    (zalk(ji,jj) > 35.0E2 ).OR.(zalk(ji,jj) <= 0.0 ) .OR.     &
                    (zdic(ji,jj) > 35.0E2 ).OR.(zdic(ji,jj) <= 0.0 ) .OR.     &
                    (f_kw660(ji,jj) > 1.0E-2 ).OR.(f_kw660(ji,jj) < 0.0 ) ) THEN
                  IF(lwp) THEN 
                      WRITE(numout,*) ' surface T = ',ztmp(ji,jj)
                      WRITE(numout,*) ' surface S = ',zsal(ji,jj)
                      WRITE(numout,*) ' surface ALK = ',zalk(ji,jj)
                      WRITE(numout,*) ' surface DIC = ',zdic(ji,jj)
                      WRITE(numout,*) ' KW660 = ',f_kw660(ji,jj)
                      WRITE(numout,*) ' atm xCO2 = ',f_xco2a(ji,jj)   
                      WRITE(numout,*) ' surface pco2w  = ',f_pco2w(ji,jj)
                      WRITE(numout,*) ' surface fco2w  = ',f_fco2w(ji,jj)
                      WRITE(numout,*) ' surface fco2a  = ',f_fco2atm(ji,jj)
                      WRITE(numout,*) ' surface co2flx = ',f_co2flux(ji,jj)
                      WRITE(numout,*) ' surface dpco2  = ',f_dpco2(ji,jj)
                      WRITE(numout,*) ' MOCSY input: ji =', mig(ji),' jj = ', mjg(jj),  &
                                       ' kt = ', kt 
                      WRITE(numout,*) 'MEDUSA - Air-Sea INPUT: unrealistic surface Carb. Chemistry'
                  ENDIF     
                  CALL ctl_stop( 'MEDUSA - Air-Sea INPUT: ',             &
                                 'unrealistic surface Carb. Chemistry -- INPUTS' )
               ENDIF     
               !!
               !! AXY (22/06/15): use Orr & Epitalon (2015) MOCSY-2 carbonate
               !!                 chemistry package; note that depth is set to
               !!                 zero in this call
               CALL mocsy_interface(ztmp(ji,jj),zsal(ji,jj),zalk(ji,jj),     &
                                    zdic(ji,jj),zsil(ji,jj),zpho(ji,jj),     &
                                    f_pp0(ji,jj),0.0,                        &
                                    gphit(ji,jj),f_kw660(ji,jj),             &
                                    f_xco2a(ji,jj),1,f_ph(ji,jj),            &
                                    f_pco2w(ji,jj),f_fco2w(ji,jj),           &
                                    f_h2co3(ji,jj),f_hco3(ji,jj),            &
                                    f_co3(ji,jj),f_omarg(ji,jj),             &
                                    f_omcal(ji,jj),f_BetaD_dum,              &
                                    f_rhosw(ji,jj),f_opres_dum,              &
                                    f_insitut_dum,f_pco2atm(ji,jj),          &
                                    f_fco2atm(ji,jj),f_schmidtco2(ji,jj),    &
                                    f_kwco2(ji,jj),f_K0(ji,jj),              &
                                    f_co2starair(ji,jj),f_co2flux(ji,jj),    &
                                    f_dpco2(ji,jj))
               !! mmol / m3 -> umol / kg
               f_TDIC(ji,jj) = (zdic(ji,jj) / f_rhosw(ji,jj)) * 1000.
               !! meq / m3 ->  ueq / kg
               f_TALK(ji,jj) = (zalk(ji,jj) / f_rhosw(ji,jj)) * 1000.
               f_dcf(ji,jj)  = f_rhosw(ji,jj)
               !! Jpalm -- 12-09-2017 -- add extra check after reccurent
               !!          carbonate failure in the coupled run.
               !!          must be associated to air-sea flux or air xCO2...
               !!          Check MOCSY outputs
               !!===================
               !! Jpalm -- 19-02-2018 -- remove the cap - only check MOCSY inputs
               !!       because of specific area in arabic sea where strangely
               !!       with core 2 forcing, ALK is lower than DIC and result in 
               !!       Enormous dpco2 - even if all carb chem caract are OK.
               !!       and this check stops the model.
               !!       --Input checks are already more than enough to stop the
               !!       model if carb chem goes crazy. 
               !!       we remove the mocsy output checks
               !!===================
               !!IF ( (f_pco2w(ji,jj) > 1.E4 ).OR.(f_pco2w(ji,jj) < 0.0 ) .OR.     &
               !!    (f_fco2w(ji,jj) > 1.E4 ).OR.(f_fco2w(ji,jj) < 0.0 ) .OR.     &   
               !!    (f_fco2atm(ji,jj) > 1.E4 ).OR.(f_fco2atm(ji,jj) < 0.0 ) .OR.     &
               !!    (f_co2flux(ji,jj) > 1.E-1 ).OR.(f_co2flux(ji,jj) < -1.E-1 ) .OR.     &
               !!    (f_dpco2(ji,jj) > 1.E4 ).OR.(f_dpco2(ji,jj) < -1.E4 ) ) THEN
               !!  IF(lwp) THEN 
               !!      WRITE(numout,*) ' surface T = ',ztmp(ji,jj)
               !!      WRITE(numout,*) ' surface S = ',zsal(ji,jj)
               !!      WRITE(numout,*) ' surface ALK = ',zalk(ji,jj)
               !!      WRITE(numout,*) ' surface DIC = ',zdic(ji,jj)
               !!      WRITE(numout,*) ' KW660 = ',f_kw660(ji,jj)
               !!      WRITE(numout,*) ' atm xCO2 = ',f_xco2a(ji,jj)   
               !!      WRITE(numout,*) ' surface pco2w  = ',f_pco2w(ji,jj)
               !!      WRITE(numout,*) ' surface fco2w  = ',f_fco2w(ji,jj)
               !!      WRITE(numout,*) ' surface fco2a  = ',f_fco2atm(ji,jj)
               !!      WRITE(numout,*) ' surface co2flx = ',f_co2flux(ji,jj)
               !!      WRITE(numout,*) ' surface dpco2  = ',f_dpco2(ji,jj)
               !!      WRITE(numout,*) ' MOCSY output: ji =', mig(ji),' jj = ', mjg(jj),  &
               !!                        ' kt = ', kt     
               !!      WRITE(numout,*) 'MEDUSA - Air-Sea OUTPUT: unrealistic surface Carb. Chemistry'
               !!  ENDIF     
               !!  CALL ctl_stop( 'MEDUSA - Air-Sea OUTPUT: ',            &
               !!                 'unrealistic surface Carb. Chemistry -- OUTPUTS' )
               !!ENDIF     
            ENDIF
         ENDDO
      ENDDO

#   if defined key_debug_medusa
               !! JPALM add carb print:
               call trc_rst_dia_stat(f_pco2w(:,:), 'f_pco2w')
               call trc_rst_dia_stat(f_fco2w(:,:), 'f_fco2w')
               call trc_rst_dia_stat(f_fco2atm(:,:), 'f_fco2atm')
               call trc_rst_dia_stat(f_schmidtco2(:,:), 'f_schmidtco2')
               call trc_rst_dia_stat(f_kwco2(:,:), 'f_kwco2')
               call trc_rst_dia_stat(f_co2starair(:,:), 'f_co2starair')
               call trc_rst_dia_stat(f_co2flux(:,:), 'f_co2flux')
               call trc_rst_dia_stat(f_dpco2(:,:), 'f_dpco2')
#   endif
#   else   

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,1) == 1) then     
               iters = 0
               !!
               !! carbon dioxide (CO2); Jerry Blackford code (ostensibly 
               !! OCMIP-2, but not)
               CALL trc_co2_medusa(ztmp(ji,jj),zsal(ji,jj),zdic(ji,jj),      &
                                   zalk(ji,jj),0.0,                          &
                                   f_kw660(ji,jj),f_xco2a(ji,jj),            &
                                   f_ph(ji,jj),                              &
                                   f_pco2w(ji,jj),f_h2co3(ji,jj),            &
                                   f_hco3(ji,jj),f_co3(ji,jj),               &
                                   f_omcal(ji,jj),f_omarg(ji,jj),            &
                                   f_co2flux(ji,jj),f_TDIC(ji,jj),           &
                                   f_TALK(ji,jj),f_dcf(ji,jj),               &
                                   f_henry(ji,jj),iters)
               !!
               !! AXY (09/01/14): removed iteration and NaN checks; these have
               !!                 been moved to trc_co2_medusa together with a
               !!                 fudge that amends erroneous values (this is
               !!                 intended to be a temporary fudge!); the
               !!                 output warnings are retained here so that
               !!                 failure position can be determined
               if (iters .eq. 25) then
                  IF(lwp) WRITE(numout,*) 'air-sea: ITERS WARNING, ',       &
                     iters, ' AT (', ji, ', ', jj, ', 1) AT ', kt
               endif
            ENDIF
         ENDDO
      ENDDO

#   endif
#  else

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,1) == 1) then
               !! AXY (18/04/13): switch off carbonate chemistry 
               !!                 calculations; provide quasi-sensible 
               !!                 alternatives
               f_ph(ji,jj)           = 8.1
               f_pco2w(ji,jj)        = f_xco2a(ji,jj)
               f_h2co3(ji,jj)        = 0.005 * zdic(ji,jj)
               f_hco3(ji,jj)         = 0.865 * zdic(ji,jj)
               f_co3(ji,jj)          = 0.130 * zdic(ji,jj)
               f_omcal(ji,jj) = 4.
               f_omarg(ji,jj) = 2.
               f_co2flux(ji,jj)      = 0.
               f_TDIC(ji,jj)         = zdic(ji,jj)
               f_TALK(ji,jj)         = zalk(ji,jj)
               f_dcf(ji,jj)          = 1.026
               f_henry(ji,jj)        = 1.
               !! AXY (23/06/15): add in some extra MOCSY diagnostics
               f_fco2w(ji,jj)        = f_xco2a(ji,jj)
! This doesn't seem to be used - marc 16/5/17
!               f_BetaD(ji,jj)        = 1.
               f_rhosw(ji,jj)        = 1.026
! This doesn't seem to be used - marc 16/5/17
!               f_opres(ji,jj)        = 0.
!               f_insitut(ji,jj)      = ztmp(ji,jj)
               f_pco2atm(ji,jj)      = f_xco2a(ji,jj)
               f_fco2atm(ji,jj)      = f_xco2a(ji,jj)
               f_schmidtco2(ji,jj)   = 660.
               f_kwco2(ji,jj)        = 0.
               f_K0(ji,jj)           = 0.
               f_co2starair(ji,jj)   = f_xco2a(ji,jj)
               f_dpco2(ji,jj)        = 0.
            ENDIF
         ENDDO
      ENDDO
#  endif

#  if defined key_axy_killco2flux
      !! AXY (18/08/17): single kill switch on air-sea CO2 flux for budget checking
      f_co2flux(:,:) = 0.
#  endif

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,1) == 1) then
               !! mmol/m2/s -> mmol/m3/d; correct for sea-ice; divide 
               !! through by layer thickness
               f_co2flux(ji,jj) = (1. - fr_i(ji,jj)) * f_co2flux(ji,jj) *    &
                                  86400. / fse3t(ji,jj,1)
               !!
               !! oxygen (O2); OCMIP-2 code
               !! AXY (23/06/15): amend input list for oxygen to account 
               !!                 for common gas transfer velocity
               CALL trc_oxy_medusa(ztmp(ji,jj),zsal(ji,jj),f_kw660(ji,jj),   &
                                   f_pp0(ji,jj),zoxy(ji,jj),                 &
                                   f_kwo2_dum,f_o2flux(ji,jj),               &
                                   f_o2sat(ji,jj))
               !!
               !! mmol/m2/s -> mol/m3/d; correct for sea-ice; divide 
               !! through by layer thickness
               f_o2flux(ji,jj)  = (1. - fr_i(ji,jj)) * f_o2flux(ji,jj) *     &
                                  86400. / fse3t(ji,jj,1)
            ENDIF
         ENDDO
      ENDDO

      !! Jpalm (08-2014)
      !! DMS surface concentration calculation
      !! initialy added for UKESM1 model.
      !! using MET-OFFICE subroutine.
      !! DMS module only needs Chl concentration and MLD
      !! to get an aproximate value of DMS concentration.
      !! air-sea fluxes are calculated by atmospheric chemitry model
      !! from atm and oc-surface concentrations.
      !!
      !! AXY (13/03/15): this is amended to calculate all of the DMS
      !!                 estimates examined during UKESM1 (see
      !!                 comments in trcdms_medusa.F90)
      !!
      !! AXY (25/05/17): amended to additionally pass DIN limitation as well as [DIN]; 
      !!                 accounts for differences in nutrient half-saturations; changes 
      !!                 also made in trc_dms_medusa; this permits an additional DMS 
      !!                 calculation while retaining the existing Anderson one 
      !! 
      IF (jdms == 1) THEN
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,1) == 1) then
                  !! calculate weighted half-saturation for DIN uptake
                  dms_wtkn(ji,jj) = ((zphn(ji,jj) * xnln) +            &
                                     (zphd(ji,jj) * xnld)) /           &
                                     (zphn(ji,jj) + zphd(ji,jj))        
                  !!
                  !! feed in correct inputs
                  if (jdms_input .eq. 0) then
                     !! use instantaneous inputs
                     dms_nlim(ji,jj) = zdin(ji,jj) / (zdin(ji,jj) + dms_wtkn(ji,jj))
                     !!
                     CALL trc_dms_medusa(zchn(ji,jj),zchd(ji,jj),             &
                                         hmld(ji,jj),qsr(ji,jj),              &
                                         zdin(ji,jj), dms_nlim(ji,jj),        &
                                         dms_andr,dms_simo,dms_aran,dms_hall, & 
                                         dms_andm)
                  else
                     !! use diel-average inputs
                     dms_nlim(ji,jj) = zn_dms_din(ji,jj) /                    &
                                      (zn_dms_din(ji,jj) + dms_wtkn(ji,jj))
                     !!
                     CALL trc_dms_medusa(zn_dms_chn(ji,jj),zn_dms_chd(ji,jj), &
                                         zn_dms_mld(ji,jj),zn_dms_qsr(ji,jj), &
                                         zn_dms_din(ji,jj),dms_nlim(ji,jj),   &
                                         dms_andr,dms_simo,dms_aran,dms_hall, & 
                                         dms_andm)
                  endif
                  !!
                  !! assign correct output to variable passed to atmosphere
                  if (jdms_model .eq. 1) then
                     dms_surf = dms_andr
                  elseif (jdms_model .eq. 2) then
                     dms_surf = dms_simo
                  elseif (jdms_model .eq. 3) then
                     dms_surf = dms_aran
                  elseif (jdms_model .eq. 4) then
                     dms_surf = dms_hall
                  elseif (jdms_model .eq. 5) then
                     dms_surf = dms_andm
                  endif
                  !!
                  !! 2D diag through iom_use
                  IF( med_diag%DMS_SURF%dgsave ) THEN
                     dms_surf2d(ji,jj) = dms_surf
                  ENDIF
                  IF( med_diag%DMS_ANDR%dgsave ) THEN
                     dms_andr2d(ji,jj) = dms_andr
                  ENDIF
                  IF( med_diag%DMS_SIMO%dgsave ) THEN
                     dms_simo2d(ji,jj) = dms_simo
                  ENDIF
                  IF( med_diag%DMS_ARAN%dgsave ) THEN
                     dms_aran2d(ji,jj) = dms_aran
                  ENDIF
                  IF( med_diag%DMS_HALL%dgsave ) THEN
                     dms_hall2d(ji,jj) = dms_hall
                  ENDIF 
                  IF( med_diag%DMS_ANDM%dgsave ) THEN 
                     dms_andm2d(ji,jj) = dms_andm
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
#   if defined key_debug_medusa
         IF (lwp) write (numout,*)                                &
            'air-sea: finish calculating dms kt = ',kt
            CALL flush(numout)
#   endif 
      ENDIF                  !! End IF (jdms == 1)

      !!
      !! store 2D outputs
      !!
      !! JPALM -- 17-11-16 -- put fgco2 out of diag request
      !!       is needed for coupling; pass through restart
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,1) == 1) then
               !! IF( med_diag%FGCO2%dgsave ) THEN
               !! convert from  mol/m2/day to kg/m2/s
               !! mmol-C/m3/d -> kg-CO2/m2/s
               fgco2(ji,jj) = f_co2flux(ji,jj) * fse3t(ji,jj,1) *             &
                              CO2flux_conv
               !! ENDIF
               IF ( lk_iomput ) THEN
                  IF( med_diag%ATM_PCO2%dgsave ) THEN
                     f_pco2a2d(ji,jj) = f_pco2atm(ji,jj)
                  ENDIF
                  IF( med_diag%OCN_PCO2%dgsave ) THEN
                     f_pco2w2d(ji,jj) = f_pco2w(ji,jj)
                  ENDIF
                  IF( med_diag%CO2FLUX%dgsave ) THEN
                     !! mmol/m3/d -> mmol/m2/d
                     f_co2flux2d(ji,jj) = f_co2flux(ji,jj) * fse3t(ji,jj,1)
                  ENDIF
                  IF( med_diag%TCO2%dgsave ) THEN
                     f_TDIC2d(ji,jj) = f_TDIC(ji,jj)
                  ENDIF
                  IF( med_diag%TALK%dgsave ) THEN
                     f_TALK2d(ji,jj) = f_TALK(ji,jj)
                  ENDIF
                  IF( med_diag%KW660%dgsave ) THEN
                      f_kw6602d(ji,jj) = f_kw660(ji,jj)
                  ENDIF
                  IF( med_diag%ATM_PP0%dgsave ) THEN
                      f_pp02d(ji,jj) = f_pp0(ji,jj)
                  ENDIF
                  IF( med_diag%O2FLUX%dgsave ) THEN
                     f_o2flux2d(ji,jj) = f_o2flux(ji,jj)
                  ENDIF
                  IF( med_diag%O2SAT%dgsave ) THEN
                     f_o2sat2d(ji,jj) = f_o2sat(ji,jj)
                  ENDIF
                  !! AXY (24/11/16): add in extra MOCSY diagnostics
                  IF( med_diag%ATM_XCO2%dgsave ) THEN
                     f_xco2a_2d(ji,jj) = f_xco2a(ji,jj)
                  ENDIF
                  IF( med_diag%OCN_FCO2%dgsave ) THEN
                     f_fco2w_2d(ji,jj) = f_fco2w(ji,jj)
                  ENDIF
                  IF( med_diag%ATM_FCO2%dgsave ) THEN
                     f_fco2a_2d(ji,jj) = f_fco2atm(ji,jj)
                  ENDIF
                  IF( med_diag%OCN_RHOSW%dgsave ) THEN
                     f_ocnrhosw_2d(ji,jj) = f_rhosw(ji,jj)
                  ENDIF
                  IF( med_diag%OCN_SCHCO2%dgsave ) THEN
                     f_ocnschco2_2d(ji,jj) = f_schmidtco2(ji,jj)
                  ENDIF
                  IF( med_diag%OCN_KWCO2%dgsave ) THEN
                     f_ocnkwco2_2d(ji,jj) = f_kwco2(ji,jj)
                  ENDIF
                  IF( med_diag%OCN_K0%dgsave ) THEN
                     f_ocnk0_2d(ji,jj) = f_K0(ji,jj)
                  ENDIF
                  IF( med_diag%CO2STARAIR%dgsave ) THEN
                     f_co2starair_2d(ji,jj) = f_co2starair(ji,jj)
                  ENDIF
                  IF( med_diag%OCN_DPCO2%dgsave ) THEN
                     f_ocndpco2_2d(ji,jj) = f_dpco2(ji,jj)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

# endif

      !!-----------------------------------------------------------------
      !! River inputs
      !!-----------------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,1) == 1) then
               !!
               !! runoff comes in as        kg / m2 / s
               !! used and written out as   m3 / m2 / d (= m / d)
               !! where                     1000 kg / m2 / d = 
               !!                             1 m3 / m2 / d = 1 m / d
               !!
               !! AXY (17/07/14): the compiler doesn't like this line for 
               !!                 some reason; as MEDUSA doesn't even use 
               !!                 runoff for riverine inputs, a temporary 
               !!                 solution is to switch off runoff entirely
               !!                 here; again, this change is one of several 
               !!                 that will need revisiting once MEDUSA has 
               !!                 bedded down in UKESM1; particularly so if 
               !!                 the land scheme provides information
               !!                 concerning nutrient fluxes
               !!
               !! f_runoff(ji,jj) = sf_rnf(1)%fnow(ji,jj,1) / 1000. * 60. *   &
               !!                   60. * 24.
               f_runoff(ji,jj) = 0.0
               !!
               !! nutrients are added via rivers to the model in one of 
               !! two ways:
               !!   1. via river concentration; i.e. the average nutrient 
               !!      concentration of a river water is described by a 
               !!      spatial file, and this is multiplied by runoff to 
               !!      give a nutrient flux
               !!   2. via direct river flux; i.e. the average nutrient 
               !!      flux due to rivers is described by a spatial file, 
               !!      and this is simply applied as a direct nutrient 
               !!      flux (i.e. it does not relate or respond to model 
               !!      runoff) nutrient fields are derived from the 
               !!      GlobalNEWS 2 database; carbon and alkalinity are 
               !!      derived from continent-scale DIC estimates (Huang et 
               !!      al., 2012) and some Arctic river alkalinity 
               !!      estimates (Katya?)
               !! 
               !! as of 19/07/12, riverine nutrients can now be spread 
               !! vertically across several grid cells rather than just 
               !! poured into the surface box; this block of code is still 
               !! executed, however, to set up the total amounts of 
               !! nutrient entering via rivers
               !!
               !! nitrogen
               if (jriver_n .eq. 1) then
                  !! river concentration specified; use runoff to 
                  !! calculate input
                  f_riv_n(ji,jj) = f_runoff(ji,jj) * riv_n(ji,jj)
               elseif (jriver_n .eq. 2) then
                  !! river flux specified; independent of runoff
                  f_riv_n(ji,jj) = riv_n(ji,jj)
               endif
               !!
               !! silicon
               if (jriver_si .eq. 1) then
                  !! river concentration specified; use runoff to 
                  !! calculate input
                  f_riv_si(ji,jj) = f_runoff(ji,jj) * riv_si(ji,jj)
               elseif (jriver_si .eq. 2) then
                  !! river flux specified; independent of runoff
                  f_riv_si(ji,jj) = riv_si(ji,jj)
               endif
               !!
               !! carbon
               if (jriver_c .eq. 1) then
                  !! river concentration specified; use runoff to 
                  !! calculate input
                  f_riv_c(ji,jj) = f_runoff(ji,jj) * riv_c(ji,jj)
               elseif (jriver_c .eq. 2) then
                  !! river flux specified; independent of runoff
                  f_riv_c(ji,jj) = riv_c(ji,jj)
               endif
               !!
               !! alkalinity
               if (jriver_alk .eq. 1) then
                  !! river concentration specified; use runoff to 
                  !! calculate input
                  f_riv_alk(ji,jj) = f_runoff(ji,jj) * riv_alk(ji,jj)
               elseif (jriver_alk .eq. 2) then
                  !! river flux specified; independent of runoff
                  f_riv_alk(ji,jj) = riv_alk(ji,jj)
               endif
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE air_sea
!CONTAINS

!   FUNCTION gas_transfer(wind, N, eqn)
!      implicit none
!      INTEGER :: N, eqn
!      real(wp) :: wind
!      real(wp) :: gas_transfer
!      real(wp) :: a(7)
!      real(wp) :: b(7)
!      real(wp) :: tmp_k
!      data a(1) / 0.166 /  ! Liss & Merlivat (1986)    [approximated]
!      data a(2) / 0.3   /  ! Wanninkhof (1992)         [sans enhancement]
!      data a(3) / 0.23  /  ! Nightingale et al. (2000) [good]
!     data a(4) / 0.23  /  ! Nightingale et al. (2000) [better]
!      data a(5) / 0.222 /  ! Nightingale et al. (2000) [best]
!      data a(6) / 0.337 /  ! OCMIP-2                   [sans variability]
!      data a(7) / 0.251 /  ! Wanninkhof (2014)         [assumes 6h avg winds]

!      data b(1) / 0.133 /
!      data b(2) / 0.0   /
!      data b(3) / 0.0   /
!      data b(4) / 0.1   /
!      data b(5) / 0.333 /
!      data b(6) / 0.0   /
!      data b(7) / 0.0   /
!      tmp_k = (a(eqn) * wind**2) + (b(eqn) * wind)
!      gas_transfer = tmp_k / (100. * 3600.)

!      END FUNCTION gas_transfer

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE air_sea( )                    ! Empty routine
      WRITE(*,*) 'air_sea: You should not have seen this print! error?'
   END SUBROUTINE air_sea
#endif 

   !!======================================================================
END MODULE air_sea_mod
