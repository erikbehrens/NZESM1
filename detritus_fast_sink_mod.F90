MODULE detritus_fast_sink_mod
   !!======================================================================
   !!                         ***  MODULE detritus_fast_sink_mod  ***
   !! Calculates fast-sinking detritus
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   detritus_fast_sink        ! Called in detritus.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE detritus_fast_sink( jk, iball )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE detritus_fast_sink  ***
      !! This called from DETRITUS and calculates the fast-sinking detritus 
      !!-------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: b0,                                     &
                                   f_benout_c, f_benout_ca, f_benout_fe,   &
                                   f_benout_lyso_ca, f_benout_n,           &
                                   f_benout_si,                            &
                                   f_fbenin_c, f_fbenin_ca, f_fbenin_fe,   &
                                   f_fbenin_n, f_fbenin_si, f_omcal,       &
                                   fccd, fdep1, fdd,                       &
                                   fdpd, fdpd2, fdpds, fdpds2,             &
                                   fdpn, fdpn2,                            &
                                   fdzme, fdzme2, fdzmi, fdzmi2,           &
                                   ffast2slowc, ffast2slown,               &
                                   ffastc, ffastca, ffastfe, ffastn,       &
                                   ffastsi,                                &
                                   fgmed, fgmepd, fgmepds, fgmepn,         &
                                   fgmezmi,                                &
                                   fgmid, fgmipn,                          &
                                   ficme, ficmi,                           &
                                   fifd_fe, fifd_n, fifd_si,               &
                                   finme, finmi,                           &
                                   fmeexcr, fmiexcr,                       &
                                   fofd_fe, fofd_n, fofd_si,               &
                                   fregen, fregenfast, fregenfastsi,       &
                                   fregensi,                               &
                                   freminc, freminca, freminfe,            &
                                   freminn, freminsi,                      &
                                   fsdiss,                                 &
                                   fsedc, fsedca, fsedn, fsedfe, fsedsi,   &
                                   fslowc, fslowcflux,                     &
                                   fslown, fslownflux,                     &
                                   ftempc, ftempca, ftempfe, ftempn,       &
                                   ftempsi,                                &
# if defined key_roam
                                   fifd_c, fofd_c, fregenfastc,            &
# endif
                                   idf, idfval,                            &
                                   zdet, zdtc
      USE dom_oce,           ONLY: e3t_0, gdepw_0, gphit, mbathy, tmask
# if defined key_vvl
      USE dom_oce,           ONLY: e3t_n, gdepw_n
# endif
      USE in_out_manager,    ONLY: lwp, numout
      USE oce,               ONLY: tsn
      USE par_kind,          ONLY: wp
      USE par_oce,           ONLY: jpi, jpim1, jpj, jpjm1
      USE sms_medusa,        ONLY: f2_ccd_cal, f3_omcal,                   &
                                   jexport, jfdfate, jinorgben, jocalccd,  &
                                   jorgben, jp_tem, jrratio,               &
                                   ocal_ccd, vsed,                         &
                                   xbetac, xbetan,                         &
                                   xcaco3a, xcaco3b,                       &
                                   xfastc, xfastca, xfastsi,               &
                                   xfdfrac1, xfdfrac2, xfdfrac3,           &
                                   xmassc, xmassca, xmasssi,               &
                                   xphi, xprotca, xprotsi,                 &
                                   xrfn, xridg_r0,                         &
                                   xsedc, xsedca, xsedfe,xsedn, xsedsi,    &
                                   xthetapd, xthetapn,                     &
                                   xthetazme, xthetazmi,                   &
                                   zn_sed_c, zn_sed_ca, zn_sed_fe,         &
                                   zn_sed_n, zn_sed_si

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! Level
      INTEGER, INTENT( in ) :: jk
      !! Fast detritus ballast scheme (0 = no; 1 = yes)
      INTEGER, INTENT( in ) :: iball

      !! Loop variables
      INTEGER :: ji, jj

      REAL(wp) :: fb_val, fl_sst
      !! Particle flux
      REAL(wp) :: fcaco3
      REAL(wp) :: fprotf
      REAL(wp), DIMENSION(jpi,jpj) :: fccd_dep
      !! temporary variables
      REAL(wp) :: fq0,fq1,fq2,fq3,fq4,fq5,fq6,fq7,fq8

      !! The two sections below, slow detritus creation and Nutrient 
      !! regeneration are moved from just above the CALL to DETRITUS
      !! in trcbio_medusa.F90.
      !!---------------------------------------------------------
      !! Slow detritus creation
      !!---------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !!
               !! this variable integrates the creation of slow sinking 
               !! detritus to allow the split between fast and slow 
               !! detritus to be diagnosed
               fslown(ji,jj)  = fdpn(ji,jj) + fdzmi(ji,jj) +                 &
                                ((1.0 - xfdfrac1) * fdpd(ji,jj)) +           &
                                ((1.0 - xfdfrac2) * fdzme(ji,jj)) +          &
                                ((1.0 - xbetan) *                            &
                                 (finmi(ji,jj) + finme(ji,jj)))
               !!
               !! this variable records the slow detrital sinking flux at 
               !! this particular depth; it is used in the output of this 
               !! flux at standard depths in the diagnostic outputs; 
               !! needs to be adjusted from per second to per day because 
               !! of parameter vsed
               fslownflux(ji,jj) = zdet(ji,jj) * vsed * 86400.
# if defined key_roam
               !!
               !! and the same for detrital carbon
               fslowc(ji,jj)  = (xthetapn * fdpn(ji,jj)) +                   &
                                (xthetazmi * fdzmi(ji,jj)) +                 &
                                (xthetapd * (1.0 - xfdfrac1) *               &
                                 fdpd(ji,jj)) +                              &
                                (xthetazme * (1.0 - xfdfrac2) *              &
                                 fdzme(ji,jj)) +                             &
                                ((1.0 - xbetac) * (ficmi(ji,jj) +            &
                                                   ficme(ji,jj)))
               !!
               !! this variable records the slow detrital sinking flux 
               !! at this particular depth; it is used in the output of 
               !! this flux at standard depths in the diagnostic 
               !! outputs; needs to be adjusted from per second to per 
               !! day because of parameter vsed
               fslowcflux(ji,jj) = zdtc(ji,jj) * vsed * 86400.
# endif
            ENDIF
         ENDDO
      ENDDO

      !!---------------------------------------------------------
      !! Nutrient regeneration
      !! this variable integrates total nitrogen regeneration down the
      !! watercolumn; its value is stored and output as a 2D diagnostic;
      !! the corresponding dissolution flux of silicon (from sources
      !! other than fast detritus) is also integrated; note that,
      !! confusingly, the linear loss terms from plankton compartments
      !! are labelled as fdX2 when one might have expected fdX or fdX1
      !!---------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !!
               !! nitrogen
               fregen(ji,jj) =                                             &
                                     ! messy feeding
                               (((xphi * (fgmipn(ji,jj) + fgmid(ji,jj))) + &
                                 (xphi *                                   &
                                  (fgmepn(ji,jj) + fgmepd(ji,jj) +         &
                                   fgmezmi(ji,jj) + fgmed(ji,jj))) +       &
                                     ! excretion + D remin.
                                 fmiexcr(ji,jj) + fmeexcr(ji,jj) +         &
                                 fdd(ji,jj) +                              &
                                     ! linear mortality
                                 fdpn2(ji,jj) + fdpd2(ji,jj) +             &
                                 fdzmi2(ji,jj) + fdzme2(ji,jj)) *          &
                                fse3t(ji,jj,jk))
               !!
               !! silicon
               fregensi(ji,jj) =                                           &
                                     ! dissolution + non-lin. mortality
                                 ((fsdiss(ji,jj) +                         &
                                   ((1.0 - xfdfrac1) * fdpds(ji,jj)) +     &
                                     ! egestion by zooplankton
                                   ((1.0 - xfdfrac3) * fgmepds(ji,jj)) +   &
                                     ! linear mortality
                                   fdpds2(ji,jj)) * fse3t(ji,jj,jk))
            ENDIF
         ENDDO
      ENDDO

      !!-------------------------------------------------------------------
      !! Fast-sinking detritus terms
      !! "local" variables declared so that conservation can be checked;
      !! the calculated terms are added to the fast-sinking flux later on
      !! only after the flux entering this level has experienced some
      !! remineralisation
      !! note: these fluxes need to be scaled by the level thickness
      !!-------------------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then

               !! nitrogen:   diatom and mesozooplankton mortality
               ftempn(ji,jj)  = b0 * ((xfdfrac1 * fdpd(ji,jj))  +            &
                                      (xfdfrac2 * fdzme(ji,jj)))
               !!
               !! silicon:    diatom mortality and grazed diatoms
               ftempsi(ji,jj) = b0 * ((xfdfrac1 * fdpds(ji,jj)) +            &
                                      (xfdfrac3 * fgmepds(ji,jj)))
               !!
               !! iron:       diatom and mesozooplankton mortality
               ftempfe(ji,jj) = b0 * (((xfdfrac1 * fdpd(ji,jj)) +            &
                                       (xfdfrac2 * fdzme(ji,jj))) * xrfn)
               !!
               !! carbon:     diatom and mesozooplankton mortality
               ftempc(ji,jj)  = b0 * ((xfdfrac1 * xthetapd * fdpd(ji,jj)) +  & 
                                      (xfdfrac2 * xthetazme * fdzme(ji,jj)))
               !!
            ENDIF
         ENDDO
      ENDDO

# if defined key_roam
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               if (jrratio.eq.0) then
                  !! CaCO3:      latitudinally-based fraction of total 
                  !!             primary production
                  !!               0.10 at equator; 0.02 at pole
                  fcaco3 = xcaco3a + ((xcaco3b - xcaco3a) *                  &
                                      ((90.0 - abs(gphit(ji,jj))) / 90.0))
               elseif (jrratio.eq.1) then
                  !! CaCO3:      Ridgwell et al. (2007) submodel, version 1
                  !!             this uses SURFACE omega calcite to regulate 
                  !!             rain ratio
                  if (f_omcal(ji,jj).ge.1.0) then
                     fq1 = (f_omcal(ji,jj) - 1.0)**0.81
                  else
                     fq1 = 0.
                  endif
                  fcaco3 = xridg_r0 * fq1
               elseif (jrratio.eq.2) then
                  !! CaCO3:      Ridgwell et al. (2007) submodel, version 2
                  !!             this uses FULL 3D omega calcite to regulate
                  !!             rain ratio
                  if (f3_omcal(ji,jj,jk).ge.1.0) then
                     fq1 = (f3_omcal(ji,jj,jk) - 1.0)**0.81
                  else
                     fq1 = 0.
                  endif
                  fcaco3 = xridg_r0 * fq1
               endif
# else
               !! CaCO3:      latitudinally-based fraction of total primary
               !!              production
               !!               0.10 at equator; 0.02 at pole
               fcaco3 = xcaco3a + ((xcaco3b - xcaco3a) *                      &
                                   ((90.0 - abs(gphit(ji,jj))) / 90.0))
# endif
               !! AXY (09/03/09): convert CaCO3 production from function of 
               !! primary production into a function of fast-sinking material; 
               !! technically, this is what Dunne et al. (2007) do anyway; they 
               !! convert total primary production estimated from surface 
               !! chlorophyll to an export flux for which they apply conversion 
               !! factors to estimate the various elemental fractions (Si, Ca)
               ftempca(ji,jj) = ftempc(ji,jj) * fcaco3

# if defined key_debug_medusa
               !! integrate total fast detritus production
               if (idf.eq.1) then
                  fifd_n(ji,jj)  = fifd_n(ji,jj)  + (ftempn(ji,jj)  *         &
                                                     fse3t(ji,jj,jk))
                  fifd_si(ji,jj) = fifd_si(ji,jj) + (ftempsi(ji,jj) *         &
                                                     fse3t(ji,jj,jk))
                  fifd_fe(ji,jj) = fifd_fe(ji,jj) + (ftempfe(ji,jj) *         &
                                                     fse3t(ji,jj,jk))
#  if defined key_roam
                  fifd_c(ji,jj)  = fifd_c(ji,jj)  + (ftempc(ji,jj)  *         &
                                                     fse3t(ji,jj,jk))
#  endif
               endif

               !! report quantities of fast-sinking detritus for each component
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
! These variables are not in this routine - marc 28/4/17
!                  IF (lwp) write (numout,*) 'fdpd(',jk,')    = ', fdpd(ji,jj)
!                  IF (lwp) write (numout,*) 'fdzme(',jk,')   = ', fdzme(ji,jj)
                  IF (lwp) write (numout,*) 'ftempn(',jk,')  = ', ftempn(ji,jj)
                  IF (lwp) write (numout,*) 'ftempsi(',jk,') = ', ftempsi(ji,jj)
                  IF (lwp) write (numout,*) 'ftempfe(',jk,') = ', ftempfe(ji,jj)
                  IF (lwp) write (numout,*) 'ftempc(',jk,')  = ', ftempc(ji,jj)
                  IF (lwp) write (numout,*) 'ftempca(',jk,') = ', ftempca(ji,jj)
                  IF (lwp) write (numout,*) 'flat(',jk,')    = ',             &
                                            abs(gphit(ji,jj))
                  IF (lwp) write (numout,*) 'fcaco3(',jk,')  = ', fcaco3
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      !!----------------------------------------------------------
      !! This version of MEDUSA offers a choice of three methods for
      !! handling the remineralisation of fast detritus.  All three
      !! do so in broadly the same way:
      !!
      !!   1.  Fast detritus is stored as a 2D array  [ ffastX  ]
      !!   2.  Fast detritus is added level-by-level  [ ftempX  ]
      !!   3.  Fast detritus is not remineralised in the top box 
      !!       [ freminX ]
      !!   4.  Remaining fast detritus is remineralised in the 
      !!       bottom  [ fsedX   ] box
      !!
      !! The three remineralisation methods are:
      !!   
      !!   1.  Ballast model (i.e. that published in Yool et al., 
      !!       2011)
      !!  (1b. Ballast-sans-ballast model)
      !!   2.  Martin et al. (1987)
      !!   3.  Henson et al. (2011)
      !! 
      !! The first of these couples C, N and Fe remineralisation to
      !! the remineralisation of particulate Si and CaCO3, but the 
      !! latter two treat remineralisation of C, N, Fe, Si and CaCO3
      !! completely separately.  At present a switch within the code
      !! regulates which submodel is used, but this should be moved
      !! to the namelist file.
      !! 
      !! The ballast-sans-ballast submodel is an original development
      !! feature of MEDUSA in which the ballast submodel's general
      !! framework and parameterisation is used, but in which there
      !! is no protection of organic material afforded by ballasting
      !! minerals.  While similar, it is not the same as the Martin 
      !! et al. (1987) submodel.
      !!
      !! Since the three submodels behave the same in terms of
      !! accumulating sinking material and remineralising it all at
      !! the seafloor, these portions of the code below are common to
      !! all three.
      !!----------------------------------------------------------
      if (jexport.eq.1) then
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  !!=======================================================
                  !! BALLAST SUBMODEL
                  !!=======================================================
                  !! 
                  !!-------------------------------------------------------
                  !! Fast-sinking detritus fluxes, pt. 1: REMINERALISATION
                  !! aside from explicitly modelled, slow-sinking detritus, the
                  !! model includes an implicit representation of detrital
                  !! particles that sink too quickly to be modelled with
                  !! explicit state variables; this sinking flux is instead
                  !! instantaneously remineralised down the water column using
                  !! the version of Armstrong et al. (2002)'s ballast model
                  !! used by Dunne et al. (2007); the version of this model
                  !! here considers silicon and calcium carbonate ballast
                  !! minerals; this section of the code redistributes the fast
                  !! sinking material generated locally down the water column;
                  !! this differs from Dunne et al. (2007) in that fast sinking
                  !! material is distributed at *every* level below that it is
                  !! generated, rather than at every level below some fixed
                  !! depth; this scheme is also different in that sinking 
                  !! material generated in one level is aggregated with that 
                  !! generated by shallower levels; this should make the 
                  !! ballast model more self-consistent (famous last words)
                  !!-------------------------------------------------------
                  !!
                  if (jk.eq.1) then
                     !! this is the SURFACE OCEAN BOX (no remineralisation)
                     !!
                     freminc(ji,jj)  = 0.0
                     freminn(ji,jj)  = 0.0
                     freminfe(ji,jj) = 0.0
                     freminsi(ji,jj) = 0.0
                     freminca(ji,jj) = 0.0
                  elseif (jk.le.mbathy(ji,jj)) then
                     !! this is an OCEAN BOX (remineralise some material)
                     !!
                     !! set up CCD depth to be used depending on user choice
                     if (jocalccd.eq.0) then
                        !! use default CCD field
                        fccd_dep(ji,jj) = ocal_ccd(ji,jj)
                     elseif (jocalccd.eq.1) then
                        !! use calculated CCD field
                        fccd_dep(ji,jj) = f2_ccd_cal(ji,jj)
                     endif
                     !!
                     !! === organic carbon ===
                     !! how much organic C enters this box        (mol)
                     fq0      = ffastc(ji,jj)
                     if (iball.eq.1) then
                        !! how much it weighs
                        fq1      = (fq0 * xmassc)
                        !! how much CaCO3 enters this box
                        fq2      = (ffastca(ji,jj) * xmassca)
                        !! how much  opal enters this box
                        fq3      = (ffastsi(ji,jj) * xmasssi)
                        !! total protected organic C
                        fq4      = (fq2 * xprotca) + (fq3 * xprotsi)
                        !! This next term is calculated for C but used for
                        !! N and Fe as well
                        !! It needs to be protected in case ALL C is protected
                        if (fq4.lt.fq1) then
                           !! protected fraction of total organic C (non-dim)
                           fprotf   = (fq4 / (fq1 + tiny(fq1)))
                        else
                           !! all organic C is protected (non-dim)
                           fprotf   = 1.0
                        endif
                        !! unprotected fraction of total organic C (non-dim)
                        fq5      = (1.0 - fprotf)
                        !! how much organic C is unprotected (mol)
                        fq6      = (fq0 * fq5)
                        !! how much unprotected C leaves this box (mol)
                        fq7      = (fq6 * exp(-(fse3t(ji,jj,jk) / xfastc)))
                        !! how much total C leaves this box (mol)
                        fq8      = (fq7 + (fq0 * fprotf))
                        !! C remineralisation in this box (mol)
                        freminc(ji,jj)  = (fq0 - fq8) / fse3t(ji,jj,jk)
                        ffastc(ji,jj) = fq8
# if defined key_debug_medusa
                        !! report in/out/remin fluxes of carbon for this level
                           if (idf.eq.1.AND.idfval.eq.1) then
                              IF (lwp) write (numout,*)                       &
                                       '------------------------------'
                              IF (lwp) write (numout,*) 'totalC(',jk,')  = ', &
                                       fq1
                              IF (lwp) write (numout,*) 'prtctC(',jk,')  = ', &
                                       fq4
                              IF (lwp) write (numout,*) 'fprotf(',jk,')  = ', &
                                       fprotf
                              IF (lwp) write (numout,*)                       &
                                       '------------------------------'
                              IF (lwp) write (numout,*) 'IN   C(',jk,')  = ', &
                                       fq0
                              IF (lwp) write (numout,*) 'LOST C(',jk,')  = ', &
                                       freminc(ji,jj) * fse3t(ji,jj,jk)
                              IF (lwp) write (numout,*) 'OUT  C(',jk,')  = ', &
                                       fq8
                              IF (lwp) write (numout,*) 'NEW  C(',jk,')  = ', &
                                       ftempc(ji,jj) * fse3t(ji,jj,jk)
                           endif
# endif
                        else
                        !! how much organic C leaves this box (mol)
                        fq1      = fq0 * exp(-(fse3t(ji,jj,jk) / xfastc))
                        !! C remineralisation in this box (mol)
                        freminc(ji,jj)  = (fq0 - fq1) / fse3t(ji,jj,jk)
                        ffastc(ji,jj)  = fq1
                     endif
                     !!
                     !! === organic nitrogen ===
                     !! how much organic N enters this box (mol)
                     fq0      = ffastn(ji,jj)
                     if (iball.eq.1) then
                        !! unprotected fraction of total organic N (non-dim)
                        fq5      = (1.0 - fprotf)
                        !! how much organic N is unprotected (mol)
                        fq6      = (fq0 * fq5)
                        !! how much unprotected N leaves this box (mol)
                        fq7      = (fq6 * exp(-(fse3t(ji,jj,jk) / xfastc)))
                        !! how much total N leaves this box (mol)
                        fq8      = (fq7 + (fq0 * fprotf))
                        !! N remineralisation in this box (mol)
                        freminn(ji,jj)  = (fq0 - fq8) / fse3t(ji,jj,jk)
                        ffastn(ji,jj)  = fq8
# if defined key_debug_medusa
                        !! report in/out/remin fluxes of carbon for this level
                        if (idf.eq.1.AND.idfval.eq.1) then
                           IF (lwp) write (numout,*)                          &
                                    '------------------------------'
                           IF (lwp) write (numout,*) 'totalN(',jk,')  = ', fq1
                           IF (lwp) write (numout,*) 'prtctN(',jk,')  = ', fq4
                           IF (lwp) write (numout,*) 'fprotf(',jk,')  = ',    &
                                    fprotf
                           IF (lwp) write (numout,*)                          &
                                    '------------------------------'
                           if (freminn(ji,jj) < 0.0) then
                              IF (lwp) write (numout,*) '** FREMIN ERROR **'
                           endif
                           IF (lwp) write (numout,*) 'IN   N(',jk,')  = ', fq0
                           IF (lwp) write (numout,*) 'LOST N(',jk,')  = ',    &
                                    freminn(ji,jj) * fse3t(ji,jj,jk)
                           IF (lwp) write (numout,*) 'OUT  N(',jk,')  = ', fq8
                           IF (lwp) write (numout,*) 'NEW  N(',jk,')  = ',    &
                                    ftempn(ji,jj) * fse3t(ji,jj,jk)
                        endif
# endif
                     else
                        !! how much organic N leaves this box (mol)
                        fq1      = fq0 * exp(-(fse3t(ji,jj,jk) / xfastc))
                        !! N remineralisation in this box (mol)
                        freminn(ji,jj)  = (fq0 - fq1) / fse3t(ji,jj,jk)
                        ffastn(ji,jj)  = fq1
                     endif
                     !!
                     !! === organic iron ===
                     !! how much organic Fe enters this box (mol)
                     fq0      = ffastfe(ji,jj)
                     if (iball.eq.1) then
                        !! unprotected fraction of total organic Fe (non-dim)
                        fq5      = (1.0 - fprotf)
                        !! how much organic Fe is unprotected (mol)
                        fq6      = (fq0 * fq5)
                        !! how much unprotected Fe leaves this box (mol)
                        fq7      = (fq6 * exp(-(fse3t(ji,jj,jk) / xfastc)))
                        !! how much total Fe leaves this box (mol)
                        fq8      = (fq7 + (fq0 * fprotf))
                        !! Fe remineralisation in this box (mol)
                        freminfe(ji,jj) = (fq0 - fq8) / fse3t(ji,jj,jk)
                        ffastfe(ji,jj) = fq8
                     else
                        !! how much total Fe leaves this box (mol)
                        fq1      = fq0 * exp(-(fse3t(ji,jj,jk) / xfastc))
                        !! Fe remineralisation in this box (mol)
                        freminfe(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                        ffastfe(ji,jj) = fq1
                     endif
                     !!
                     !! === biogenic silicon ===
                     !! how much  opal centers this box (mol)
                     fq0      = ffastsi(ji,jj)
                     !! how much  opal leaves this box (mol)
                     fq1      = fq0 * exp(-(fse3t(ji,jj,jk) / xfastsi))
                     !! Si remineralisation in this box (mol)
                     freminsi(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                     ffastsi(ji,jj) = fq1
                     !!
                     !! === biogenic calcium carbonate ===
                     !! how much CaCO3 enters this box (mol)
                     fq0      = ffastca(ji,jj)
                     if (fsdepw(ji,jj,jk).le.fccd_dep(ji,jj)) then
                        !! whole grid cell above CCD
                        !! above lysocline, no Ca dissolves (mol)
                        fq1      = fq0
                        !! above lysocline, no Ca dissolves (mol)
                        freminca(ji,jj) = 0.0
                        !! which is the last level above the CCD?    (#)
                        fccd(ji,jj) = real(jk)
                     elseif (fsdepw(ji,jj,jk).ge.fccd_dep(ji,jj)) then
                        !! whole grid cell below CCD
                        !! how much CaCO3 leaves this box (mol)
                        fq1      = fq0 * exp(-(fse3t(ji,jj,jk) / xfastca))
                        !! Ca remineralisation in this box (mol)
                        freminca(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                     else
                        !! partial grid cell below CCD
                        !! amount of grid cell below CCD (m)
                        fq2      = fdep1(ji,jj) - fccd_dep(ji,jj)
                        !! how much CaCO3 leaves this box (mol)
                        fq1      = fq0 * exp(-(fq2 / xfastca))
                        !! Ca remineralisation in this box (mol)
                        freminca(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                     endif
                     ffastca(ji,jj) = fq1 
                  else
                     !! this is BELOW THE LAST OCEAN BOX (do nothing)
                     freminc(ji,jj)  = 0.0
                     freminn(ji,jj)  = 0.0
                     freminfe(ji,jj) = 0.0
                     freminsi(ji,jj) = 0.0
                     freminca(ji,jj) = 0.0              
                  endif
               ENDIF
            ENDDO
         ENDDO
      elseif (jexport.eq.2.or.jexport.eq.3) then
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  if (jexport.eq.2) then
                     !!====================================================
                     !! MARTIN ET AL. (1987) SUBMODEL
                     !!====================================================
                     !! 
                     !!----------------------------------------------------
                     !! This submodel uses the classic Martin et al. (1987) 
                     !! curve to determine the attenuation of fast-sinking 
                     !! detritus down the water column.  All three organic 
                     !! elements, C, N and Fe, are handled identically, and 
                     !! their quantities in sinking particles attenuate 
                     !! according to a power relationship governed by 
                     !! parameter "b".  This is assigned a canonical value 
                     !! of -0.858.  Biogenic opal and calcium carbonate are
                     !! attentuated using the same function as in the
                     !! ballast submodel
                     !!----------------------------------------------------
                     !!
                     fb_val = -0.858
                  elseif (jexport.eq.3) then
                     !!====================================================
                     !! HENSON ET AL. (2011) SUBMODEL
                     !!====================================================
                     !!
                     !!----------------------------------------------------
                     !! This submodel reconfigures the Martin et al. (1987) 
                     !! curve by allowing the "b" value to vary 
                     !! geographically.  Its value is set, following Henson 
                     !! et al. (2011), as a function of local sea surface 
                     !! temperature:
                     !!   b = -1.06 + (0.024 * SST)
                     !! This means that remineralisation length scales are 
                     !! longer in warm, tropical areas and shorter in cold, 
                     !! polar areas.  This does seem back-to-front (i.e. 
                     !! one would expect GREATER remineralisation in warmer 
                     !! waters), but is an outcome of analysis of sediment 
                     !! trap data, and it may reflect details of ecosystem 
                     !! structure that pertain to particle production
                     !! rather than simply Q10.
                     !!----------------------------------------------------
                     !!
                     fl_sst = tsn(ji,jj,1,jp_tem)
                     fb_val = -1.06 + (0.024 * fl_sst)
                  endif
                  !!   
                  if (jk.eq.1) then
                     !! this is the SURFACE OCEAN BOX (no remineralisation)
                     !!
                     freminc(ji,jj)  = 0.0
                     freminn(ji,jj)  = 0.0
                     freminfe(ji,jj) = 0.0
                     freminsi(ji,jj) = 0.0
                     freminca(ji,jj) = 0.0
                  elseif (jk.le.mbathy(ji,jj)) then
                     !! this is an OCEAN BOX (remineralise some material)
                     !!
                     !! === organic carbon ===
                     !! how much organic C enters this box (mol)
                     fq0      = ffastc(ji,jj)
                     !! how much organic C leaves this box (mol)
                     fq1      = fq0 * ((fdep1(ji,jj)/fsdepw(ji,jj,jk))**fb_val)
                     !! C remineralisation in this box (mol)
                     freminc(ji,jj)  = (fq0 - fq1) / fse3t(ji,jj,jk)
                     ffastc(ji,jj)  = fq1
                     !!
                     !! === organic nitrogen ===
                     !! how much organic N enters this box (mol)
                     fq0      = ffastn(ji,jj)
                     !! how much organic N leaves this box (mol)
                     fq1      = fq0 * ((fdep1(ji,jj)/fsdepw(ji,jj,jk))**fb_val)
                     !! N remineralisation in this box (mol)
                     freminn(ji,jj)  = (fq0 - fq1) / fse3t(ji,jj,jk)
                     ffastn(ji,jj)  = fq1
                     !!
                     !! === organic iron ===
                     !! how much organic Fe enters this box (mol)
                     fq0      = ffastfe(ji,jj)
                     !! how much organic Fe leaves this box (mol)
                     fq1      = fq0 * ((fdep1(ji,jj)/fsdepw(ji,jj,jk))**fb_val)
                     !! Fe remineralisation in this box (mol)
                     freminfe(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                     ffastfe(ji,jj) = fq1
                     !!
                     !! === biogenic silicon ===
                     !! how much  opal centers this box (mol)
                     fq0      = ffastsi(ji,jj)
                     !! how much  opal leaves this box (mol)
                     fq1      = fq0 * exp(-(fse3t(ji,jj,jk) / xfastsi))
                     !! Si remineralisation in this box (mol)
                     freminsi(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                     ffastsi(ji,jj) = fq1
                     !!
                     !! === biogenic calcium carbonate ===
                     !! how much CaCO3 enters this box (mol)
                     fq0      = ffastca(ji,jj)
                     if (fsdepw(ji,jj,jk).le.ocal_ccd(ji,jj)) then
                        !! whole grid cell above CCD
                        !! above lysocline, no Ca dissolves (mol)
                        fq1      = fq0
                        !! above lysocline, no Ca dissolves (mol)
                        freminca(ji,jj) = 0.0
                        !! which is the last level above the CCD?    (#)
                        fccd(ji,jj) = real(jk)
                     elseif (fsdepw(ji,jj,jk).ge.ocal_ccd(ji,jj)) then
                        !! whole grid cell below CCD
                        !! how much CaCO3 leaves this box (mol)
                        fq1      = fq0 * exp(-(fse3t(ji,jj,jk) / xfastca))
                        !! Ca remineralisation in this box (mol)
                        freminca(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                     else
                        !! partial grid cell below CCD
                        !! amount of grid cell below CCD (m)
                        fq2      = fdep1(ji,jj) - ocal_ccd(ji,jj)
                        !! how much CaCO3 leaves this box (mol)
                        fq1      = fq0 * exp(-(fq2 / xfastca))
                        !! Ca remineralisation in this box (mol)
                        freminca(ji,jj) = (fq0 - fq1) / fse3t(ji,jj,jk)
                     endif
                     ffastca(ji,jj) = fq1 
                  else
                     !! this is BELOW THE LAST OCEAN BOX (do nothing)
                     freminc(ji,jj)  = 0.0
                     freminn(ji,jj)  = 0.0
                     freminfe(ji,jj) = 0.0
                     freminsi(ji,jj) = 0.0
                     freminca(ji,jj) = 0.0              
                  endif
               ENDIF
            ENDDO
         ENDDO
      endif

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Fast-sinking detritus fluxes, pt. 2: UPDATE FAST FLUXES
               !! here locally calculated additions to the fast-sinking 
               !! flux are added to the total fast-sinking flux; this is 
               !! done here such that material produced in a particular 
               !! layer is only remineralised below this layer
               !!----------------------------------------------------------
               !!
               !! add sinking material generated in this layer to running 
               !! totals
               !!
               !! === organic carbon === 
               !! (diatom and mesozooplankton mortality)
               ffastc(ji,jj)  = ffastc(ji,jj)  + (ftempc(ji,jj)  *           &
                                                  fse3t(ji,jj,jk))
               !!
               !! === organic nitrogen ===
               !! (diatom and mesozooplankton mortality)
               ffastn(ji,jj)  = ffastn(ji,jj)  + (ftempn(ji,jj)  *           &
                                                  fse3t(ji,jj,jk))
               !!
               !! === organic iron ===
               !! (diatom and mesozooplankton mortality)
               ffastfe(ji,jj) = ffastfe(ji,jj) + (ftempfe(ji,jj) *          &
                                                  fse3t(ji,jj,jk))
               !!
               !! === biogenic silicon ===
               !! (diatom mortality and grazed diatoms)
               ffastsi(ji,jj) = ffastsi(ji,jj) + (ftempsi(ji,jj) *          &
                                                  fse3t(ji,jj,jk))
               !!
               !! === biogenic calcium carbonate ===
               !! (latitudinally-based fraction of total primary production)
               ffastca(ji,jj) = ffastca(ji,jj) + (ftempca(ji,jj) *          &
                                                  fse3t(ji,jj,jk))
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Fast-sinking detritus fluxes, pt. 3: SEAFLOOR
               !! remineralise all remaining fast-sinking detritus to dissolved
               !! nutrients; the sedimentation fluxes calculated here allow the
               !! separation of what's remineralised sinking through the final
               !! ocean box from that which is added to the final box by the
               !! remineralisation of material that reaches the seafloor (i.e.
               !! the model assumes that *all* material that hits the seafloor
               !! is remineralised and that none is permanently buried; hey,
               !! this is a giant GCM model that can't be run for long enough
               !! to deal with burial fluxes!)
               !!
               !! in a change to this process, in part so that MEDUSA behaves
               !! a little more like ERSEM et al., fast-sinking detritus (N, Fe
               !! and C) is converted to slow sinking detritus at the seafloor
               !! instead of being remineralised; the rationale is that in
               !! shallower shelf regions (... that are not fully mixed!) this
               !! allows the detrital material to return slowly to dissolved 
               !! nutrient rather than instantaneously as now; the alternative
               !! would be to explicitly handle seafloor organic material - a
               !! headache I don't wish to experience at this point; note that
               !! fast-sinking Si and Ca detritus is just remineralised as 
               !! per usual
               !! 
               !! AXY (13/01/12)
               !! in a further change to this process, again so that MEDUSA is
               !! a little more like ERSEM et al., material that reaches the
               !! seafloor can now be added to sediment pools and stored for
               !! slow release; there are new 2D arrays for organic nitrogen,
               !! iron and carbon and inorganic silicon and carbon that allow
               !! fast and slow detritus that reaches the seafloor to be held
               !! and released back to the water column more slowly; these
               !! arrays are transferred via the tracer restart files between
               !! repeat submissions of the model
               !!----------------------------------------------------------
               !! 
               ffast2slowc(ji,jj)  = 0.0
               ffast2slown(ji,jj)  = 0.0
! I don't think this is used - marc 10/4/17
!               ffast2slowfe(ji,jj) = 0.0
               !!
               if (jk.eq.mbathy(ji,jj)) then
                  !! this is the BOTTOM OCEAN BOX (remineralise everything)
                  !!
                  !! AXY (17/01/12): tweaked to include benthos pools
                  !! 
                  !! === organic carbon ===
                  if (jfdfate.eq.0 .and. jorgben.eq.0) then
                     !! C remineralisation in this box (mol/m3)
                     freminc(ji,jj)  = freminc(ji,jj) + (ffastc(ji,jj) /     &
                                                         fse3t(ji,jj,jk))
                  elseif (jfdfate.eq.1 .and. jorgben.eq.0) then
                     !! fast C -> slow C (mol/m3)
                     ffast2slowc(ji,jj) = ffastc(ji,jj) / fse3t(ji,jj,jk)
                     fslowc(ji,jj)      = fslowc(ji,jj) + ffast2slowc(ji,jj)
                  elseif (jfdfate.eq.0 .and. jorgben.eq.1) then
                     !! fast C -> benthic C (mol/m2)
                     f_fbenin_c(ji,jj)  = ffastc(ji,jj)
                  endif
                  !! record seafloor C (mol/m2)
                  fsedc(ji,jj)   = ffastc(ji,jj)
                  ffastc(ji,jj)  = 0.0
                  !!
                  !! === organic nitrogen ===
                  if (jfdfate.eq.0 .and. jorgben.eq.0) then
                     !! N remineralisation in this box (mol/m3)
                     freminn(ji,jj)  = freminn(ji,jj) + (ffastn(ji,jj) /     &
                                                         fse3t(ji,jj,jk))
                  elseif (jfdfate.eq.1 .and. jorgben.eq.0) then
                     !! fast N -> slow N (mol/m3)
                     ffast2slown(ji,jj) = ffastn(ji,jj) / fse3t(ji,jj,jk)
                     fslown(ji,jj)      = fslown(ji,jj) + ffast2slown(ji,jj)
                  elseif (jfdfate.eq.0 .and. jorgben.eq.1) then
                     !! fast N -> benthic N (mol/m2)
                     f_fbenin_n(ji,jj)  = ffastn(ji,jj)
                  endif
                  !! record seafloor N (mol/m2)
                  fsedn(ji,jj)   = ffastn(ji,jj)
                  ffastn(ji,jj)  = 0.0
                  !!
                  !! === organic iron ===
                  if (jfdfate.eq.0 .and. jorgben.eq.0) then
                     !! Fe remineralisation in this box (mol/m3)
                     freminfe(ji,jj) = freminfe(ji,jj) + (ffastfe(ji,jj) /   &
                                                          fse3t(ji,jj,jk))
! I don't think ffast2slowfe is used - marc 10/4/17
!                  elseif (jfdfate.eq.1 .and. jorgben.eq.0) then
!                     !! fast Fe -> slow Fe (mol/m3)
!                     ffast2slowfe(ji,jj) = ffastn(ji,jj) / fse3t(ji,jj,jk)
                  elseif (jfdfate.eq.0 .and. jorgben.eq.1) then
                     !! fast Fe -> benthic Fe (mol/m2)
                     f_fbenin_fe(ji,jj) = ffastfe(ji,jj)
                  endif
                  !! record seafloor Fe (mol/m2)
                  fsedfe(ji,jj)  = ffastfe(ji,jj)
                  ffastfe(ji,jj) = 0.0
                  !!
                  !! === biogenic silicon ===
                  if (jinorgben.eq.0) then
                     !! Si remineralisation in this box (mol/m3)
                     freminsi(ji,jj) = freminsi(ji,jj) + (ffastsi(ji,jj) /   &
                                                          fse3t(ji,jj,jk))
                  elseif (jinorgben.eq.1) then
                     !! fast Si -> benthic Si
                     f_fbenin_si(ji,jj) = ffastsi(ji,jj)
                  endif
                  !! record seafloor Si (mol/m2)
                  fsedsi(ji,jj)   = ffastsi(ji,jj)
                  ffastsi(ji,jj) = 0.0
                  !!
                  !! === biogenic calcium carbonate ===
                  if (jinorgben.eq.0) then
                     !! Ca remineralisation in this box (mol/m3)
                     freminca(ji,jj) = freminca(ji,jj) + (ffastca(ji,jj) /   &
                                                          fse3t(ji,jj,jk))
                  elseif (jinorgben.eq.1) then
                     !! fast Ca -> benthic Ca (mol/m2)
                     f_fbenin_ca(ji,jj) = ffastca(ji,jj)
                  endif
                  !! record seafloor Ca (mol/m2)
                  fsedca(ji,jj)   = ffastca(ji,jj)
                  ffastca(ji,jj) = 0.0
               endif

# if defined key_debug_medusa
               if (idf.eq.1) then
                  !!-------------------------------------------------------
                  !! Integrate total fast detritus remineralisation
                  !!-------------------------------------------------------
                  !!
                  fofd_n(ji,jj)  = fofd_n(ji,jj)  + (freminn(ji,jj)  *       &
                                                     fse3t(ji,jj,jk))
                  fofd_si(ji,jj) = fofd_si(ji,jj) + (freminsi(ji,jj) *       &
                                                     fse3t(ji,jj,jk))
                  fofd_fe(ji,jj) = fofd_fe(ji,jj) + (freminfe(ji,jj) *       &
                                                     fse3t(ji,jj,jk))
#  if defined key_roam
                  fofd_c(ji,jj)  = fofd_c(ji,jj)  + (freminc(ji,jj)  *       &
                                                     fse3t(ji,jj,jk))
#  endif
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Sort out remineralisation tally of fast-sinking detritus
               !!----------------------------------------------------------
               !!
               !! update fast-sinking regeneration arrays
               fregenfast(ji,jj)   = fregenfast(ji,jj)   +                  &
                                     (freminn(ji,jj)  * fse3t(ji,jj,jk))
               fregenfastsi(ji,jj) = fregenfastsi(ji,jj) +                  &
                                     (freminsi(ji,jj) * fse3t(ji,jj,jk))
# if defined key_roam
               fregenfastc(ji,jj)  = fregenfastc(ji,jj)  +                  &
                                     (freminc(ji,jj)  * fse3t(ji,jj,jk))
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Benthic remineralisation fluxes
               !!----------------------------------------------------------
               !!
               if (jk.eq.mbathy(ji,jj)) then
                  !!
                  !! organic components
                  if (jorgben.eq.1) then
                     f_benout_n(ji,jj)  = xsedn  * zn_sed_n(ji,jj)
                     f_benout_fe(ji,jj) = xsedfe * zn_sed_fe(ji,jj)
                     f_benout_c(ji,jj)  = xsedc  * zn_sed_c(ji,jj)
                  endif
                  !!
                  !! inorganic components
                  if (jinorgben.eq.1) then
                     f_benout_si(ji,jj) = xsedsi * zn_sed_si(ji,jj)
                     f_benout_ca(ji,jj) = xsedca * zn_sed_ca(ji,jj)
                     !!
                     !! account for CaCO3 that dissolves when it shouldn't
                     if ( fsdepw(ji,jj,jk) .le. fccd_dep(ji,jj) ) then
                        f_benout_lyso_ca(ji,jj) = xsedca * zn_sed_ca(ji,jj)
                     endif
                  endif
               endif
               CALL flush(numout)

            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE detritus_fast_sink

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE detritus_fast_sink( )                    ! Empty routine
      WRITE(*,*) 'detritus_fast_sink: You should not have seen this print! error?'
   END SUBROUTINE detritus_fast_sink
#endif 

   !!======================================================================
END MODULE detritus_fast_sink_mod
