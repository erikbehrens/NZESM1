MODULE bio_medusa_update_mod
   !!======================================================================
   !!                         ***  MODULE bio_medusa_update_mod  ***
   !! Update tracer fields
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!   -   ! 2017-08 (A. Yool)            Amend slow-detritus bug
   !!   -   ! 2017-08 (A. Yool)            Reformatting for clarity
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   bio_medusa_update        ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE bio_medusa_update( kt, jk )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE bio_medusa_update  ***
      !! This called from TRC_BIO_MEDUSA and updates the tracer fields
      !!---------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: b0, bddtalk3, bddtdic3, bddtdife3,        &
                                   bddtdin3, bddtdisi3,                      &
                                   ibenthic, ibio_switch, idf, idfval,       &
                                   f_benout_c, f_benout_ca, f_benout_n,      &
                                   f_benout_si,                              &
                                   f_co2flux, f_o2flux,                      &
                                   f_riv_loc_alk, f_riv_loc_c,               &
                                   f_riv_loc_n, f_riv_loc_si,                & 
                                   f_riv_alk, f_riv_c, f_riv_n, f_riv_si,    &
                                   fbddtalk, fbddtdic, fbddtdife,            &
                                   fbddtdin, fbddtdisi,                      & 
                                   fdd, fdpd, fdpd2, fdpds, fdpds2,          &
                                   fdpn, fdpn2,                              &
                                   fdzme, fdzme2, fdzmi, fdzmi2,             &
                                   ffast2slowc, ffast2slown,                 &
                                   ffebot, ffetop, ffescav,                  &
                                   fflx_fe, fflx_n, fflx_si,                 &
                                   fgmed, fgmepd, fgmedc, fgmepd, fgmepds,   &
                                   fgmepn, fgmezmi,                          &
                                   fgmid, fgmidc, fgmipn,                    &
                                   ficme, ficmi, finme, finmi,               &
                                   fmeexcr, fmegrow, fmeresp,                &
                                   fmiexcr, fmigrow, fmiresp,                &
                                   fnit_cons, fnit_prod,                     &
                                   fprd, fprds, fprn,                        &
                                   frd,                                      &
                                   freminc, freminca, freminn, freminsi,     &
                                   frn,                                      &
                                   fsil_cons, fsil_prod, fsdiss,             &
                                   ftempca, fthetad, fthetan,                &
                                   fslowsink, fslowgain, fslowloss,          & ! AXY (22/08/17)
                                   f_sbenin_n, f_sbenin_c,                   &
# if defined key_roam
                                   fslowsinkc, fslowgainc, fslowlossc,       & ! AXY (22/08/17)
                                   fcar_cons, fcar_prod, fcomm_resp,         &
                                   fddc, fflx_a, fflx_c, fflx_o2, zoxy,      &
                                   foxy_anox, foxy_cons, foxy_prod,          &
# endif
                                   zpds, zphd, zphn
      USE dom_oce,           ONLY: e3t_0, gphit, mbathy, tmask
# if defined key_vvl
      USE dom_oce,           ONLY: e3t_n
# endif
      USE in_out_manager,    ONLY: lwp, numout
      USE lib_mpp,           ONLY: ctl_stop
      USE par_kind,          ONLY: wp
      USE par_medusa,        ONLY: jp_medusa, jp_msa0, jp_msa1,              &
                                   jpalk, jpchd, jpchn, jpdet, jpdic,        &
                                   jpdin, jpdtc, jpfer, jpoxy, jppds,        &
                                   jpphd, jpphn, jpsil, jpzme, jpzmi,        &
                                   jpalk_lc, jpchd_lc, jpchn_lc, jpdet_lc,   & 
                                   jpdic_lc, jpdin_lc, jpdtc_lc, jpfer_lc,   & 
                                   jpoxy_lc, jppds_lc, jpphd_lc, jpphn_lc,   &
                                   jpsil_lc, jpzme_lc, jpzmi_lc
      USE par_oce,           ONLY: jpi, jpim1, jpj, jpjm1, jpk
      USE par_trc,           ONLY: jptra
      USE sms_medusa,        ONLY: friver_dep,                               &
                                   jinorgben, jorgben,                       &
                                   jriver_alk, jriver_c,                     &
                                   jriver_n, jriver_si,                      &
                                   xbetac, xbetan,                           &
                                   xfdfrac1, xfdfrac2, xfdfrac3,             &
                                   xo2min, xphi, xrfn,                       &
                                   xthetanit, xthetapd, xthetapn,            &
                                   xthetarem, xthetazme, xthetazmi,          &
                                   xxi
     !CEB  USE trc,               ONLY: med_diag, tra
     USE trc
     !/CEB

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! time (integer timestep)
      INTEGER, INTENT( in ) :: kt
      !! Level
      INTEGER, INTENT( in ) :: jk

      !! Loop variables
      INTEGER :: ji, jj, jn

      !! AXY (23/08/13): changed from individual variables for each flux to
      !!                 an array that holds all fluxes
      REAL(wp), DIMENSION(jpi,jpj,jp_medusa) :: btra

      !! nitrogen and silicon production and consumption
      REAL(wp) :: fn_prod, fn_cons, fs_prod, fs_cons

      !! carbon, alkalinity production and consumption
      REAL(wp) :: fc_prod, fc_cons, fa_prod, fa_cons

      !! oxygen production and consumption (and non-consumption)
      REAL(wp), DIMENSION(jpi,jpj) :: fo2_prod, fo2_cons
      REAL(wp), DIMENSION(jpi,jpj) :: fo2_ncons, fo2_ccons

      !! temporary variables
      REAL(wp) :: fq0

      !!==========================================================
      !! LOCAL GRID CELL TRENDS
      !!==========================================================
      !!
      !!----------------------------------------------------------
      !! Determination of trends
      !!----------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then
               !!
               !!----------------------------------------------------------
               !! chlorophyll
               btra(ji,jj,jpchn_lc) = b0 * ( ( (frn(ji,jj) * fprn(ji,jj) *      &
                                             zphn(ji,jj) ) -                 &
                                           fgmipn(ji,jj) - fgmepn(ji,jj) -   &
                                           fdpn(ji,jj) - fdpn2(ji,jj) ) *    &
                                          (fthetan(ji,jj) / xxi) )
               btra(ji,jj,jpchd_lc) = b0 * ( ( (frd(ji,jj) * fprd(ji,jj) *      &
                                             zphd(ji,jj) ) -                 &
                                           fgmepd(ji,jj) - fdpd(ji,jj) -     &
                                           fdpd2(ji,jj) ) *                  &
                                          (fthetad(ji,jj) / xxi) )
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !!----------------------------------------------------------
               !! phytoplankton
               btra(ji,jj,jpphn_lc) = b0 * ( (fprn(ji,jj) * zphn(ji,jj)) -      &
                                          fgmipn(ji,jj) - fgmepn(ji,jj) -    &
                                          fdpn(ji,jj) - fdpn2(ji,jj) )
               btra(ji,jj,jpphd_lc) = b0 * ( (fprd(ji,jj) * zphd(ji,jj)) -      &
                                          fgmepd(ji,jj) - fdpd(ji,jj) -      &
                                          fdpd2(ji,jj) )
               btra(ji,jj,jppds_lc) = b0 * ( (fprds(ji,jj) * zpds(ji,jj)) -     &
                                          fgmepds(ji,jj) - fdpds(ji,jj) -    &
                                          fsdiss(ji,jj) - fdpds2(ji,jj) )
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !!----------------------------------------------------------
               !! zooplankton
               btra(ji,jj,jpzmi_lc) = b0 * (fmigrow(ji,jj) - fgmezmi(ji,jj) -   &
                                         fdzmi(ji,jj) - fdzmi2(ji,jj))
               btra(ji,jj,jpzme_lc) = b0 * (fmegrow(ji,jj) - fdzme(ji,jj) -     &
                                         fdzme2(ji,jj))
            ENDIF
         ENDDO
      ENDDO

      !!----------------------------------------------------------
      !! detritus
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               btra(ji,jj,jpdet_lc) = b0 * (                           &
                   fdpn(ji,jj)                                         & ! mort. losses 
                 + ((1.0 - xfdfrac1) * fdpd(ji,jj))                    & ! mort. losses 
                 + fdzmi(ji,jj)                                        & ! mort. losses
                 + ((1.0 - xfdfrac2) * fdzme(ji,jj))                   & ! mort. losses
                 + ((1.0 - xbetan) * (finmi(ji,jj) + finme(ji,jj)))    & ! assim. inefficiency
                 - fgmid(ji,jj) - fgmed(ji,jj)                         & ! grazing
                 - fdd(ji,jj)                                          & ! remin.
                 + fslowgain(ji,jj) - fslowloss(ji,jj)                 & ! slow-sinking
                 - (f_sbenin_n(ji,jj) / fse3t(ji,jj,jk))               & ! slow-sinking loss to seafloor
                 + ffast2slown(ji,jj) )                                  ! seafloor fast->slow 
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! dissolved inorganic nitrogen nutrient
               !! primary production
               fn_cons = - (fprn(ji,jj) * zphn(ji,jj)) -                     &
                           (fprd(ji,jj) * zphd(ji,jj))
               !! 
               fn_prod =                                                     &
                                         ! messy feeding remin.
                         (xphi * (fgmipn(ji,jj) + fgmid(ji,jj))) +           &
                         (xphi * (fgmepn(ji,jj) + fgmepd(ji,jj) +            &
                                  fgmezmi(ji,jj) + fgmed(ji,jj))) +          &
                                         ! excretion and remin.
                         fmiexcr(ji,jj) + fmeexcr(ji,jj) + fdd(ji,jj) +      &
                         freminn(ji,jj) +                                    &
                                         ! metab. losses
                         fdpn2(ji,jj) + fdpd2(ji,jj) + fdzmi2(ji,jj) +       &
                         fdzme2(ji,jj)
               !! 
               !! riverine flux
               if ( jriver_n .gt. 0 ) then
                  f_riv_loc_n(ji,jj) = f_riv_n(ji,jj) *                      &
                     friver_dep(jk,mbathy(ji,jj)) / fse3t(ji,jj,jk)
                  fn_prod = fn_prod + f_riv_loc_n(ji,jj)
               endif
               !!  
               !! benthic remineralisation
               if (jk.eq.mbathy(ji,jj) .and. jorgben.eq.1 .and.              &
                   ibenthic.eq.1) then
                  fn_prod = fn_prod + (f_benout_n(ji,jj) / fse3t(ji,jj,jk))
               endif
               !!
               btra(ji,jj,jpdin_lc) = b0 * ( fn_prod + fn_cons )
               !!
               !! consumption of dissolved nitrogen
               fnit_cons(ji,jj) = fnit_cons(ji,jj) + ( fse3t(ji,jj,jk) *     &
                                                       fn_cons )
               !! production of dissolved nitrogen
               fnit_prod(ji,jj) = fnit_prod(ji,jj) + ( fse3t(ji,jj,jk) *     &
                                                       fn_prod )
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !!----------------------------------------------------------
               !! dissolved silicic acid nutrient
               !! opal production
               fs_cons = - (fprds(ji,jj) * zpds(ji,jj))
               !!
               fs_prod =                                                     &
                             ! opal dissolution
                         fsdiss(ji,jj) +                                     &
                             ! mort. loss
                         ((1.0 - xfdfrac1) * fdpds(ji,jj)) +                 &
                             ! egestion of grazed Si
                         ((1.0 - xfdfrac3) * fgmepds(ji,jj)) +               &
                             ! fast diss. and metab. losses
                         freminsi(ji,jj) + fdpds2(ji,jj)
               !! 
               !! riverine flux
               if ( jriver_si .gt. 0 ) then
                  f_riv_loc_si(ji,jj) = f_riv_si(ji,jj) *                    &
                                        friver_dep(jk,mbathy(ji,jj)) /       &
                                        fse3t(ji,jj,jk)
                  fs_prod = fs_prod + f_riv_loc_si(ji,jj)
               endif
               !!  
               !! benthic remineralisation
               if (jk.eq.mbathy(ji,jj) .and. jinorgben.eq.1 .and.            &
                   ibenthic.eq.1) then
                  fs_prod = fs_prod + (f_benout_si(ji,jj) / fse3t(ji,jj,jk))
               endif
               !!
               btra(ji,jj,jpsil_lc) = b0 * ( &
                 fs_prod + fs_cons )
               !! consumption of dissolved silicon
               fsil_cons(ji,jj) = fsil_cons(ji,jj) + ( fse3t(ji,jj,jk) *     &
                                                       fs_cons )
               !! production of dissolved silicon
               fsil_prod(ji,jj) = fsil_prod(ji,jj) + ( fse3t(ji,jj,jk) *     &
                                                       fs_prod )
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then !!
               !!----------------------------------------------------------
               !! dissolved "iron" nutrient
               btra(ji,jj,jpfer_lc) = b0 * ( (xrfn * btra(ji,jj,jpdin_lc)) +       &
                                          ffetop(ji,jj) + ffebot(ji,jj) -    &
                                          ffescav(ji,jj) )
            ENDIF
         ENDDO
      ENDDO

# if defined key_roam
      !!----------------------------------------------------------
      !! AXY (26/11/08): implicit detrital carbon change
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then 
               !!
               btra(ji,jj,jpdtc_lc) = b0 * (                           &
                   (xthetapn * fdpn(ji,jj))                            & ! mort. losses 
                 + ((1.0 - xfdfrac1) * (xthetapd * fdpd(ji,jj)))       & ! mort. losses 
                 + (xthetazmi * fdzmi(ji,jj))                          & ! mort. losses 
                 + ((1.0 - xfdfrac2) * (xthetazme * fdzme(ji,jj)))     & ! mort. losses 
                 + ((1.0 - xbetac) * (ficmi(ji,jj) + ficme(ji,jj)))    & ! assim. inefficiency
                 - fgmidc(ji,jj) - fgmedc(ji,jj)                       & ! grazing
                 - fddc(ji,jj)                                         & ! remin.
                 + fslowgainc(ji,jj) - fslowlossc(ji,jj)               & ! slow-sinking
                 - (f_sbenin_c(ji,jj) / fse3t(ji,jj,jk))               & ! slow-sinking loss to seafloor
                 + ffast2slowc(ji,jj) )                                  ! seafloor fast->slow
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !!----------------------------------------------------------
               !! dissolved inorganic carbon
               !! primary production
               fc_cons = - (xthetapn * fprn(ji,jj) * zphn(ji,jj)) -          &
                           (xthetapd * fprd(ji,jj) * zphd(ji,jj))
               !!
               fc_prod =                                                     &
                            ! messy feeding remin
                         (xthetapn * xphi * fgmipn(ji,jj)) +                 &
                         (xphi * fgmidc(ji,jj)) +                            &
                         (xthetapn * xphi * fgmepn(ji,jj)) +                 &
                         (xthetapd * xphi * fgmepd(ji,jj)) +                 &
                         (xthetazmi * xphi * fgmezmi(ji,jj)) +               &
                         (xphi * fgmedc(ji,jj)) +                            &
                            ! resp., remin., losses
                         fmiresp(ji,jj) + fmeresp(ji,jj) + fddc(ji,jj) +     &
                         freminc(ji,jj) + (xthetapn * fdpn2(ji,jj)) +        &
                            ! losses
                         (xthetapd * fdpd2(ji,jj)) +                         &
                         (xthetazmi * fdzmi2(ji,jj)) +                       &
                         (xthetazme * fdzme2(ji,jj))
               !! 
               !! riverine flux
               if ( jriver_c .gt. 0 ) then
                  f_riv_loc_c(ji,jj) = f_riv_c(ji,jj) *                      &
                                       friver_dep(jk,mbathy(ji,jj)) /        &
                                       fse3t(ji,jj,jk)
                  fc_prod = fc_prod + f_riv_loc_c(ji,jj)
               endif
               !!  
               !! benthic remineralisation
               if (jk.eq.mbathy(ji,jj) .and. jorgben.eq.1 .and.              &
                   ibenthic.eq.1) then
                  fc_prod = fc_prod + (f_benout_c(ji,jj) / fse3t(ji,jj,jk))
               endif
               if (jk.eq.mbathy(ji,jj) .and. jinorgben.eq.1 .and.            &
                   ibenthic.eq.1) then
                  fc_prod = fc_prod + (f_benout_ca(ji,jj) / fse3t(ji,jj,jk))
               endif
               !!
               !! community respiration (does not include CaCO3 terms - 
               !! obviously!)
               fcomm_resp(ji,jj) = fcomm_resp(ji,jj) + fc_prod
               !!
               !! CaCO3
               fc_prod = fc_prod - ftempca(ji,jj) + freminca(ji,jj)
               !! 
               !! riverine flux
               if ( jk .eq. 1 .and. jriver_c .gt. 0 ) then
                  fc_prod = fc_prod + f_riv_c(ji,jj)
               endif
               !!
               btra(ji,jj,jpdic_lc) = b0 * ( fc_prod + fc_cons )
               !! consumption of dissolved carbon
               fcar_cons(ji,jj) = fcar_cons(ji,jj) + ( fse3t(ji,jj,jk) *     &
                                                       fc_cons )
               !! production of dissolved carbon
               fcar_prod(ji,jj) = fcar_prod(ji,jj) + ( fse3t(ji,jj,jk) *     &
                                                       fc_prod )
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !!----------------------------------------------------------
               !! alkalinity
               !! CaCO3 dissolution
               fa_prod = 2.0 * freminca(ji,jj)
               !! CaCO3 production
               fa_cons = - 2.0 * ftempca(ji,jj)
               !! 
               !! riverine flux
               if ( jriver_alk .gt. 0 ) then
                  f_riv_loc_alk(ji,jj) = f_riv_alk(ji,jj) *                  &
                                         friver_dep(jk,mbathy(ji,jj)) /      &
                                         fse3t(ji,jj,jk)
                  fa_prod = fa_prod + f_riv_loc_alk(ji,jj)
               endif
               !!  
               !! benthic remineralisation
               if (jk.eq.mbathy(ji,jj) .and. jinorgben.eq.1 .and.            &
                   ibenthic.eq.1) then
                  fa_prod = fa_prod + (2.0 * f_benout_ca(ji,jj) /            &
                                       fse3t(ji,jj,jk))
               endif
               !!
               btra(ji,jj,jpalk_lc) = b0 * ( fa_prod + fa_cons )
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !!----------------------------------------------------------
               !! oxygen (has protection at low O2 concentrations; 
               !! OCMIP-2 style)
               fo2_prod(ji,jj) =                                             &
                                     ! Pn primary production, N
                                 (xthetanit * fprn(ji,jj) * zphn(ji,jj)) +   &
                                     ! Pd primary production, N
                                 (xthetanit * fprd(ji,jj) * zphd(ji,jj)) +   &
                                     ! Pn primary production, C
                                 (xthetarem * xthetapn * fprn(ji,jj) *       &
                                  zphn(ji,jj)) +                             &
                                     ! Pd primary production, C
                                  (xthetarem * xthetapd * fprd(ji,jj) *      &
                                   zphd(ji,jj))
               fo2_ncons(ji,jj) =                                            &
                                     ! Pn messy feeding remin., N 
                                   - (xthetanit * xphi * fgmipn(ji,jj))      &
                                     ! D  messy feeding remin., N
                                   - (xthetanit * xphi * fgmid(ji,jj))       &
                                     ! Pn messy feeding remin., N
                                   - (xthetanit * xphi * fgmepn(ji,jj))      &
                                     ! Pd messy feeding remin., N
                                   - (xthetanit * xphi * fgmepd(ji,jj))      &
                                     ! Zi messy feeding remin., N
                                   - (xthetanit * xphi * fgmezmi(ji,jj))     &
                                     ! D  messy feeding remin., N
                                   - (xthetanit * xphi * fgmed(ji,jj))       &
                                     ! microzoo excretion, N
                                   - (xthetanit * fmiexcr(ji,jj))            &
                                     ! mesozoo  excretion, N
                                   - (xthetanit * fmeexcr(ji,jj))            &
                                     ! slow detritus remin., N
                                   - (xthetanit * fdd(ji,jj))                &
                                     ! fast detritus remin., N
                                   - (xthetanit * freminn(ji,jj))            &
                                     ! Pn  losses, N
                                   - (xthetanit * fdpn2(ji,jj))              &
                                     ! Pd  losses, N
                                   - (xthetanit * fdpd2(ji,jj))              &
                                     ! Zmi losses, N
                                   - (xthetanit * fdzmi2(ji,jj))             &
                                     ! Zme losses, N
                                   - (xthetanit * fdzme2(ji,jj))
               !!  
               !! benthic remineralisation
               if (jk.eq.mbathy(ji,jj) .and. jorgben.eq.1 .and.              &
                   ibenthic.eq.1) then
                  fo2_ncons(ji,jj) = fo2_ncons(ji,jj) -                      &
                                     (xthetanit * f_benout_n(ji,jj) /        &
                                      fse3t(ji,jj,jk))
               endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               fo2_ccons(ji,jj) =                                            &
                                     ! Pn messy feeding remin., C
                                  - (xthetarem * xthetapn * xphi *           &
                                     fgmipn(ji,jj))                          &
                                     ! D  messy feeding remin., C
                                  - (xthetarem * xphi * fgmidc(ji,jj))       &
                                     ! Pn messy feeding remin., C
                                  - (xthetarem * xthetapn * xphi *           &
                                     fgmepn(ji,jj))                          &
                                     ! Pd messy feeding remin., C
                                  - (xthetarem * xthetapd * xphi *           &
                                     fgmepd(ji,jj))                          &
                                     ! Zi messy feeding remin., C
                                  - (xthetarem * xthetazmi * xphi *          &
                                     fgmezmi(ji,jj))                         &
                                     ! D  messy feeding remin., C
                                  - (xthetarem * xphi * fgmedc(ji,jj))       &
                                     ! microzoo respiration, C
                                  - (xthetarem * fmiresp(ji,jj))             &
                                     ! mesozoo  respiration, C
                                  - (xthetarem * fmeresp(ji,jj))             &
                                     ! slow detritus remin., C
                                  - (xthetarem * fddc(ji,jj))                &
                                     ! fast detritus remin., C
                                  - (xthetarem * freminc(ji,jj))             &
                                     ! Pn  losses, C
                                  - (xthetarem * xthetapn * fdpn2(ji,jj))    &
                                     ! Pd  losses, C
                                  - (xthetarem * xthetapd * fdpd2(ji,jj))    &
                                     ! Zmi losses, C
                                  - (xthetarem * xthetazmi * fdzmi2(ji,jj))  &
                                     ! Zme losses, C
                                  - (xthetarem * xthetazme * fdzme2(ji,jj))
               !!  
               !! benthic remineralisation
               if (jk.eq.mbathy(ji,jj) .and. jorgben.eq.1 .and.              &
                   ibenthic.eq.1) then
                  fo2_ccons(ji,jj) = fo2_ccons(ji,jj) - (xthetarem *         &
                                                         f_benout_c(ji,jj) / &
                                                         fse3t(ji,jj,jk))
               endif
               fo2_cons(ji,jj) = fo2_ncons(ji,jj) + fo2_ccons(ji,jj)
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! is this a suboxic zone?
               !! deficient O2; production fluxes only
               if (zoxy(ji,jj).lt.xo2min) then
                  btra(ji,jj,jpoxy_lc) = b0 * fo2_prod(ji,jj)
                  foxy_prod(ji,jj) = foxy_prod(ji,jj) + ( fse3t(ji,jj,jk) *  &
                                                          fo2_prod(ji,jj) )
                  foxy_anox(ji,jj) = foxy_anox(ji,jj) + ( fse3t(ji,jj,jk) *  &
                                                          fo2_cons(ji,jj) )
               else
                  !! sufficient O2; production + consumption fluxes
                  btra(ji,jj,jpoxy_lc) = b0 * ( fo2_prod(ji,jj) +               &
                                             fo2_cons(ji,jj) )
                  foxy_prod(ji,jj) = foxy_prod(ji,jj) +                      &
                                     ( fse3t(ji,jj,jk) * fo2_prod(ji,jj) )
                  foxy_cons(ji,jj) = foxy_cons(ji,jj) +                      &
                                     ( fse3t(ji,jj,jk) * fo2_cons(ji,jj) )
               endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! air-sea fluxes (if this is the surface box)
               if (jk.eq.1) then
                  !!
                  !! CO2 flux
                  btra(ji,jj,jpdic_lc) = btra(ji,jj,jpdic_lc) + (b0 *              &
                                                           f_co2flux(ji,jj))
                  !!
                  !! O2 flux (mol/m3/s -> mmol/m3/d)
                  btra(ji,jj,jpoxy_lc) = btra(ji,jj,jpoxy_lc) + (b0 *              &
                                                           f_o2flux(ji,jj))
               endif
            ENDIF
         ENDDO
      ENDDO
# endif

# if defined key_debug_medusa
! I DON'T THIS IS MUCH USE, NOW IT'S BEEN PULLED OUT OF THE MAIN DO LOOP
! - marc 5/5/17
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !! report state variable fluxes (not including 
               !! fast-sinking detritus)
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpchn_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpchn_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpchd_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpchd_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpphn_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpphn_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpphd_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpphd_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jppds_lc)(',jk,')  = ', &
                                            btra(ji,jj,jppds_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpzmi_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpzmi_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpzme_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpzme_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpdet_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpdet_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpdin_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpdin_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpsil_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpsil_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpfer_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpfer_lc)
#  if defined key_roam
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpdtc_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpdtc_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpdic_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpdic_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpalk_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpalk_lc)
                  IF (lwp) write (numout,*) 'btra(ji,jj,jpoxy_lc)(',jk,')  = ', &
                                            btra(ji,jj,jpoxy_lc)
#  endif
               endif
            ENDIF
         ENDDO
      ENDDO
# endif

      !!----------------------------------------------------------
      !! Integrate calculated fluxes for mass balance
      !!----------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! === nitrogen ===
               fflx_n(ji,jj)  = fflx_n(ji,jj) + fse3t(ji,jj,jk) *            &
                                ( btra(ji,jj,jpphn_lc) + btra(ji,jj,jpphd_lc) +    &
                                  btra(ji,jj,jpzmi_lc) + btra(ji,jj,jpzme_lc) +    &
                                  btra(ji,jj,jpdet_lc) + btra(ji,jj,jpdin_lc) )
               !! === silicon ===
               fflx_si(ji,jj) = fflx_si(ji,jj) + fse3t(ji,jj,jk) *           &
                                ( btra(ji,jj,jppds_lc) + btra(ji,jj,jpsil_lc) )
               !! === iron ===
               fflx_fe(ji,jj) = fflx_fe(ji,jj) + fse3t(ji,jj,jk) *           &
                                ( (xrfn *                                    &
                                   (btra(ji,jj,jpphn_lc) + btra(ji,jj,jpphd_lc) +  &
                                    btra(ji,jj,jpzmi_lc) + btra(ji,jj,jpzme_lc) +  &
                                    btra(ji,jj,jpdet_lc))) + btra(ji,jj,jpfer_lc) )
# if defined key_roam
               !! === carbon ===
               fflx_c(ji,jj)  = fflx_c(ji,jj) + fse3t(ji,jj,jk) *            &
                                ( (xthetapn * btra(ji,jj,jpphn_lc)) +           &
                                  (xthetapd * btra(ji,jj,jpphd_lc)) +           &
                                  (xthetazmi * btra(ji,jj,jpzmi_lc)) +          &
                                  (xthetazme * btra(ji,jj,jpzme_lc)) +          &
                                  btra(ji,jj,jpdtc_lc) + btra(ji,jj,jpdic_lc) )
               !! === alkalinity ===
               fflx_a(ji,jj)  = fflx_a(ji,jj) + fse3t(ji,jj,jk) *            &
                                btra(ji,jj,jpalk_lc)
               !! === oxygen ===
               fflx_o2(ji,jj) = fflx_o2(ji,jj) + fse3t(ji,jj,jk) *           &
                                btra(ji,jj,jpoxy_lc)
# endif
            ENDIF
         ENDDO
      ENDDO

      !!----------------------------------------------------------
      !! Apply calculated tracer fluxes
      !!----------------------------------------------------------
      !!
      !! units: [unit of tracer] per second (fluxes are calculated 
      !! above per day)
      !!
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               ibio_switch = 1
# if defined key_gulf_finland
               !! AXY (17/05/13): fudge in a Gulf of Finland correction; 
               !!                 uses longitude-latitude range to 
               !!                 establish if this is a Gulf of Finland 
               !!                 grid cell; if so, then BGC fluxes are 
               !!                 ignored (though still calculated); for 
               !!                 reference, this is meant to be a 
               !!                 temporary fix to see if all of my 
               !!                 problems can be done away with if I 
               !!                 switch off BGC fluxes in the Gulf of 
               !!                 Finland, which currently appears the 
               !!                 source of trouble
               if ( glamt(ji,jj).gt.24.7 .and. glamt(ji,jj).lt.27.8 .and.    &
                    gphit(ji,jj).gt.59.2 .and. gphit(ji,jj).lt.60.2 ) then
                  ibio_switch = 0
               endif
# endif               
               if (ibio_switch.eq.1) then
                  tra(ji,jj,jk,jpchn) = tra(ji,jj,jk,jpchn) +                &
                                        (btra(ji,jj,jpchn_lc) / 86400.)
                  tra(ji,jj,jk,jpchd) = tra(ji,jj,jk,jpchd) +                &
                                        (btra(ji,jj,jpchd_lc) / 86400.)
                  tra(ji,jj,jk,jpphn) = tra(ji,jj,jk,jpphn) +                &
                                        (btra(ji,jj,jpphn_lc) / 86400.)
                  tra(ji,jj,jk,jpphd) = tra(ji,jj,jk,jpphd) +                &
                                        (btra(ji,jj,jpphd_lc) / 86400.)
                  tra(ji,jj,jk,jppds) = tra(ji,jj,jk,jppds) +                &
                                        (btra(ji,jj,jppds_lc) / 86400.)
                  tra(ji,jj,jk,jpzmi) = tra(ji,jj,jk,jpzmi) +                &
                                        (btra(ji,jj,jpzmi_lc) / 86400.)
                  tra(ji,jj,jk,jpzme) = tra(ji,jj,jk,jpzme) +                &
                                        (btra(ji,jj,jpzme_lc) / 86400.)
                  tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) +                &
                                        (btra(ji,jj,jpdet_lc) / 86400.)
                  tra(ji,jj,jk,jpdin) = tra(ji,jj,jk,jpdin) +                &
                                        (btra(ji,jj,jpdin_lc) / 86400.)
                  tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) +                &
                                        (btra(ji,jj,jpsil_lc) / 86400.)
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) +                &
                                        (btra(ji,jj,jpfer_lc) / 86400.)
# if defined key_roam
                  tra(ji,jj,jk,jpdtc) = tra(ji,jj,jk,jpdtc) +                &
                                        (btra(ji,jj,jpdtc_lc) / 86400.)
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +                &
                                        (btra(ji,jj,jpdic_lc) / 86400.)
                  tra(ji,jj,jk,jpalk) = tra(ji,jj,jk,jpalk) +                &
                                        (btra(ji,jj,jpalk_lc) / 86400.)
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) +                &
                                        (btra(ji,jj,jpoxy_lc) / 86400.)
# endif
               endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then

               !! AXY (18/11/16): CMIP6 diagnostics
               IF( med_diag%FBDDTALK%dgsave )  THEN
                  fbddtalk(ji,jj)  =  fbddtalk(ji,jj)  +                     &
                                      (btra(ji,jj,jpalk_lc) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%FBDDTDIC%dgsave )  THEN
                  fbddtdic(ji,jj)  =  fbddtdic(ji,jj)  +                     &
                                      (btra(ji,jj,jpdic_lc) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%FBDDTDIFE%dgsave ) THEN
                  fbddtdife(ji,jj) =  fbddtdife(ji,jj) +                     &
                                      (btra(ji,jj,jpfer_lc) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%FBDDTDIN%dgsave )  THEN
                  fbddtdin(ji,jj)  =  fbddtdin(ji,jj)  +                     &
                                      (btra(ji,jj,jpdin_lc) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%FBDDTDISI%dgsave ) THEN
                  fbddtdisi(ji,jj) =  fbddtdisi(ji,jj) +                     &
                                      (btra(ji,jj,jpsil_lc) * fse3t(ji,jj,jk))
               ENDIF
	       !!
               IF( med_diag%BDDTALK3%dgsave )  THEN
                  bddtalk3(ji,jj,jk)  =  btra(ji,jj,jpalk_lc)
               ENDIF
               IF( med_diag%BDDTDIC3%dgsave )  THEN
                  bddtdic3(ji,jj,jk)  =  btra(ji,jj,jpdic_lc)
               ENDIF
               IF( med_diag%BDDTDIFE3%dgsave ) THEN
                  bddtdife3(ji,jj,jk) =  btra(ji,jj,jpfer_lc)
               ENDIF
               IF( med_diag%BDDTDIN3%dgsave )  THEN
                  bddtdin3(ji,jj,jk)  =  btra(ji,jj,jpdin_lc)
               ENDIF
               IF( med_diag%BDDTDISI3%dgsave ) THEN
                  bddtdisi3(ji,jj,jk) =  btra(ji,jj,jpsil_lc)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

#   if defined key_debug_medusa
      IF (lwp) write (numout,*) '------'
      IF (lwp) write (numout,*) 'bio_medusa_update: end all calculations'
      IF (lwp) write (numout,*) 'bio_medusa_update: now outputs kt = ', kt
      CALL flush(numout)
#   endif

# if defined key_axy_nancheck
      !!----------------------------------------------------------
      !! Check calculated tracer fluxes
      !!----------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               DO jn = 1,jp_medusa
                  fq0 = btra(ji,jj,jn)
                  !! AXY (30/01/14): "isnan" problem on HECTOR
                  !! if (fq0 /= fq0 ) then
                  if ( ieee_is_nan( fq0 ) ) then
                     !! there's a NaN here
                     if (lwp) write(numout,*) 'NAN detected in btra(ji,jj,',  &
                        ji, ',', jj, ',', jk, ',', jn, ') at time', kt
		     CALL ctl_stop( 'trcbio_medusa, NAN in btra field' )
                  endif
               ENDDO
               DO jn = jp_msa0,jp_msa1
                  fq0 = tra(ji,jj,jk,jn)
                  !! AXY (30/01/14): "isnan" problem on HECTOR
                  !! if (fq0 /= fq0 ) then
                  if ( ieee_is_nan( fq0 ) ) then
                     !! there's a NaN here
                     if (lwp) write(numout,*) 'NAN detected in tra(', ji, &
                        ',', jj, ',', jk, ',', jn, ') at time', kt
   		     CALL ctl_stop( 'trcbio_medusa, NAN in tra field' )
                  endif
               ENDDO
               CALL flush(numout)
            ENDIF
         ENDDO
      ENDDO
# endif


   END SUBROUTINE bio_medusa_update

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE bio_medusa_update( )                    ! Empty routine
      WRITE(*,*) 'bio_medusa_update: You should not have seen this print! error?'
   END SUBROUTINE bio_medusa_update
#endif 

   !!======================================================================
END MODULE bio_medusa_update_mod
