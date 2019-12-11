MODULE bio_med_diag_iomput_mod
   !!======================================================================
   !!                         ***  MODULE bio_med_diag_iomput_mod  ***
   !! Calculates diagnostics
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
      
   PUBLIC   bio_med_diag_iomput        ! Called in bio_medusa_diag.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   SUBROUTINE bio_med_diag_iomput( jk )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE bio_med_diag_iomput  ***
      !! Calculates the diagnostics used with iom_put
      !!-------------------------------------------------------------------
      USE bio_medusa_mod
      USE dom_oce,           ONLY: e3t_0, mbathy, tmask
# if defined key_vvl
      USE dom_oce,           ONLY: e3t_n
# endif
      USE in_out_manager,    ONLY: lwp, numout
      USE par_oce,           ONLY: jpim1, jpjm1
      USE phycst,            ONLY: rsmall
      USE sms_medusa,        ONLY: f3_omarg, f3_omcal,                   &
                                   i0100, i0500, i1000,                  &
                                   xbetac, xbetan, xphi,                 &
                                   xthetapd, xthetapn, xthetazmi
      !CEB USE trc,               ONLY: med_diag
      USE trc
      !/CEB
   !!* Substitution
#  include "domzgr_substitute.h90"

      !! level
      INTEGER, INTENT( in ) :: jk

      !! Loop avariables
      INTEGER :: ji, jj, jn

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            IF (tmask(ji,jj,jk) == 1) THEN
               !!-------------------------------------------------------
               !! Add in XML diagnostics stuff
               !!-------------------------------------------------------
               !!
               !! ** 2D diagnostics
#   if defined key_debug_medusa
               IF (lwp) write (numout,*)                                     &
                  'bio_med_diag_iomput: in ij-jj loop jk = ', jk
               CALL flush(numout)
#   endif
               IF ( med_diag%PRN%dgsave ) THEN
                   fprn2d(ji,jj) = fprn2d(ji,jj) +                           &
                                   (fprn(ji,jj)  * zphn(ji,jj) *             &
                                    fse3t(ji,jj,jk)) 
               ENDIF
               IF ( med_diag%MPN%dgsave ) THEN
                   fdpn2d(ji,jj) = fdpn2d(ji,jj) + (fdpn(ji,jj) *            &
                                                    fse3t(ji,jj,jk))
               ENDIF
               IF ( med_diag%PRD%dgsave ) THEN
                   fprd2d(ji,jj) = fprd2d(ji,jj) +                           &
                                   (fprd(ji,jj)  * zphd(ji,jj) *             &
                                    fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%MPD%dgsave ) THEN
                   fdpd2d(ji,jj) = fdpd2d(ji,jj) + (fdpd(ji,jj) *            &
                                                    fse3t(ji,jj,jk)) 
               ENDIF
               !  IF( med_diag%DSED%dgsave ) THEN
               !      CALL iom_put( "DSED"  , ftot_n )
               !  ENDIF
               IF( med_diag%OPAL%dgsave ) THEN
                   fprds2d(ji,jj) = fprds2d(ji,jj) +                         &
                                    (fprds(ji,jj) * zpds(ji,jj) *            &
                                     fse3t(ji,jj,jk)) 
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%OPALDISS%dgsave ) THEN
                   fsdiss2d(ji,jj) = fsdiss2d(ji,jj) + (fsdiss(ji,jj) *      &
                                                        fse3t(ji,jj,jk))  
               ENDIF
               IF( med_diag%GMIPn%dgsave ) THEN
                   fgmipn2d(ji,jj) = fgmipn2d(ji,jj) +                       &
                                     (fgmipn(ji,jj)  * fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%GMID%dgsave ) THEN
                   fgmid2d(ji,jj) = fgmid2d(ji,jj) + (fgmid(ji,jj) *         &
                                                      fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%MZMI%dgsave ) THEN
                   fdzmi2d(ji,jj) = fdzmi2d(ji,jj) + (fdzmi(ji,jj) *         &
                                                      fse3t(ji,jj,jk)) 
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%GMEPN%dgsave ) THEN
                   fgmepn2d(ji,jj) = fgmepn2d(ji,jj) + (fgmepn(ji,jj) *      &
                                                        fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%GMEPD%dgsave ) THEN
                   fgmepd2d(ji,jj) = fgmepd2d(ji,jj) + (fgmepd(ji,jj) *      &
                                                        fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%GMEZMI%dgsave ) THEN
                   fgmezmi2d(ji,jj) = fgmezmi2d(ji,jj) +                     &
                                      (fgmezmi(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%GMED%dgsave ) THEN
                   fgmed2d(ji,jj) = fgmed2d(ji,jj) +                         &
                                       (fgmed(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%MZME%dgsave ) THEN
                  fdzme2d(ji,jj) = fdzme2d(ji,jj) +                          &
                                    (fdzme(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
               !  IF( med_diag%DEXP%dgsave ) THEN
               !      CALL iom_put( "DEXP"  , ftot_n )
               !  ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%DETN%dgsave ) THEN
                  fslown2d(ji,jj) = fslown2d(ji,jj) +                        &
                                    (fslown(ji,jj) * fse3t(ji,jj,jk))  
               ENDIF
               IF( med_diag%MDET%dgsave ) THEN
                  fdd2d(ji,jj) = fdd2d(ji,jj) +                              &
                                 (fdd(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%AEOLIAN%dgsave ) THEN
                  ffetop2d(ji,jj) = ffetop2d(ji,jj) +                        &
                                    (ffetop(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%BENTHIC%dgsave ) THEN
                  ffebot2d(ji,jj) = ffebot2d(ji,jj) +                        &
                                    (ffebot(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%SCAVENGE%dgsave ) THEN
                  ffescav2d(ji,jj) = ffescav2d(ji,jj) +                      &
                                     (ffescav(ji,jj) * fse3t(ji,jj,jk))  
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%PN_JLIM%dgsave ) THEN
                  ! fjln2d(ji,jj) = fjln2d(ji,jj) +                          &
                  !                 (fjln(ji,jj)  * zphn(ji,jj) *            &
                  !                  fse3t(ji,jj,jk)) 
                  fjln2d(ji,jj) = fjln2d(ji,jj) +                            &
                                  (fjlim_pn(ji,jj) * zphn(ji,jj) *           &
                                   fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%PN_NLIM%dgsave ) THEN
                  fnln2d(ji,jj) = fnln2d(ji,jj) +                            &
                                  (fnln(ji,jj) * zphn(ji,jj) *               &
                                   fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%PN_FELIM%dgsave ) THEN
                  ffln2d(ji,jj) = ffln2d(ji,jj) +                            &
                                  (ffln2(ji,jj) * zphn(ji,jj) *              &
                                   fse3t(ji,jj,jk)) 
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%PD_JLIM%dgsave ) THEN
                   ! fjld2d(ji,jj) = fjld2d(ji,jj) +                          &
                   !                 (fjld(ji,jj)  * zphd(ji,jj) *            &
                   !                  fse3t(ji,jj,jk)) 
                   fjld2d(ji,jj) = fjld2d(ji,jj) +                           &
                                   (fjlim_pd(ji,jj) * zphd(ji,jj) *          &
                                    fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%PD_NLIM%dgsave ) THEN
                   fnld2d(ji,jj) = fnld2d(ji,jj) +                           &
                                   (fnld(ji,jj) * zphd(ji,jj) *              &
                                    fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%PD_FELIM%dgsave ) THEN
                   ffld2d(ji,jj) = ffld2d(ji,jj) +                           &
                                   (ffld(ji,jj) * zphd(ji,jj) *              &
                                    fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%PD_SILIM%dgsave ) THEN
                   fsld2d2(ji,jj) = fsld2d2(ji,jj) +                         &
                                    (fsld2(ji,jj) * zphd(ji,jj) *            &
                                     fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%PDSILIM2%dgsave ) THEN
                   fsld2d(ji,jj) = fsld2d(ji,jj) +                           &
                                   (fsld(ji,jj) * zphd(ji,jj) *              &
                                    fse3t(ji,jj,jk))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !! 
               IF( med_diag%TOTREG_N%dgsave ) THEN
                  fregen2d(ji,jj) = fregen2d(ji,jj) + fregen(ji,jj)
               ENDIF
               IF( med_diag%TOTRG_SI%dgsave ) THEN
                  fregensi2d(ji,jj) = fregensi2d(ji,jj) + fregensi(ji,jj)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !! 
               IF( med_diag%FASTN%dgsave ) THEN
                   ftempn2d(ji,jj) = ftempn2d(ji,jj) +                       &
                                     (ftempn(ji,jj)  * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%FASTSI%dgsave ) THEN
                   ftempsi2d(ji,jj) = ftempsi2d(ji,jj) +                     &
                                      (ftempsi(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%FASTFE%dgsave ) THEN
                   ftempfe2d(ji,jj) = ftempfe2d(ji,jj) +                     &
                                      (ftempfe(ji,jj) * fse3t(ji,jj,jk))  
               ENDIF
               IF( med_diag%FASTC%dgsave ) THEN
                   ftempc2d(ji,jj) = ftempc2d(ji,jj) +                       &
                                     (ftempc(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%FASTCA%dgsave ) THEN
                   ftempca2d(ji,jj) = ftempca2d(ji,jj) +                     &
                                      (ftempca(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !! 
               IF( med_diag%REMINN%dgsave ) THEN
                   freminn2d(ji,jj) = freminn2d(ji,jj) +                     &
                                      (freminn(ji,jj)  * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%REMINSI%dgsave ) THEN
                   freminsi2d(ji,jj) = freminsi2d(ji,jj) +                   &
                                       (freminsi(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%REMINFE%dgsave ) THEN
                   freminfe2d(ji,jj) = freminfe2d(ji,jj) +                   &
                                       (freminfe(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%REMINC%dgsave ) THEN
                   freminc2d(ji,jj) = freminc2d(ji,jj) +                     &
                                      (freminc(ji,jj)  * fse3t(ji,jj,jk)) 
               ENDIF
               IF( med_diag%REMINCA%dgsave ) THEN
                   freminca2d(ji,jj) = freminca2d(ji,jj) +                   &
                                       (freminca(ji,jj) * fse3t(ji,jj,jk)) 
               ENDIF
               !!
            ENDIF
         ENDDO
      ENDDO

# if defined key_roam
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !!
               !! AXY (09/11/16): CMIP6 diagnostics
               IF( med_diag%FD_NIT3%dgsave ) THEN
                  fd_nit3(ji,jj,jk) = ffastn(ji,jj)
               ENDIF
               IF( med_diag%FD_SIL3%dgsave ) THEN
                  fd_sil3(ji,jj,jk) = ffastsi(ji,jj)
               ENDIF
               IF( med_diag%FD_CAR3%dgsave ) THEN
                  fd_car3(ji,jj,jk) = ffastc(ji,jj)
               ENDIF
               IF( med_diag%FD_CAL3%dgsave ) THEN
                  fd_cal3(ji,jj,jk) = ffastca(ji,jj)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF (jk.eq.i0100) THEN
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  IF( med_diag%RR_0100%dgsave ) THEN
                     ffastca2d(ji,jj) =                                      &
                                ffastca(ji,jj)/MAX(ffastc(ji,jj), rsmall)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE IF (jk.eq.i0500) THEN
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  IF( med_diag%RR_0500%dgsave ) THEN
                     ffastca2d(ji,jj) =                                      &
                                ffastca(ji,jj)/MAX(ffastc(ji,jj), rsmall)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE IF (jk.eq.i1000) THEN
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  IF( med_diag%RR_1000%dgsave ) THEN
                     ffastca2d(ji,jj) =                                      &
                                ffastca(ji,jj)/MAX(ffastc(ji,jj), rsmall)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (jk.eq.mbathy(ji,jj)) THEN
                  IF (tmask(ji,jj,jk) == 1) THEN
                     IF( med_diag%IBEN_N%dgsave ) THEN
                        iben_n2d(ji,jj) = f_sbenin_n(ji,jj) +                &
                                          f_fbenin_n(ji,jj)
                     ENDIF
                     IF( med_diag%IBEN_FE%dgsave ) THEN
                        iben_fe2d(ji,jj) = f_sbenin_fe(ji,jj) +              &
                                           f_fbenin_fe(ji,jj)
                     ENDIF
                     IF( med_diag%IBEN_C%dgsave ) THEN
                        iben_c2d(ji,jj) = f_sbenin_c(ji,jj) +                &
                                          f_fbenin_c(ji,jj)
                     ENDIF
                     IF( med_diag%IBEN_SI%dgsave ) THEN
                        iben_si2d(ji,jj) = f_fbenin_si(ji,jj)
                     ENDIF
                     IF( med_diag%IBEN_CA%dgsave ) THEN
                        iben_ca2d(ji,jj) = f_fbenin_ca(ji,jj)
                     ENDIF
                     IF( med_diag%OBEN_N%dgsave ) THEN
                        oben_n2d(ji,jj) = f_benout_n(ji,jj)
                     ENDIF
                     IF( med_diag%OBEN_FE%dgsave ) THEN
                        oben_fe2d(ji,jj) = f_benout_fe(ji,jj)
                     ENDIF
                     IF( med_diag%OBEN_C%dgsave ) THEN
                        oben_c2d(ji,jj) = f_benout_c(ji,jj)
                     ENDIF
                     IF( med_diag%OBEN_SI%dgsave ) THEN
                        oben_si2d(ji,jj) = f_benout_si(ji,jj)
                     ENDIF
                     IF( med_diag%OBEN_CA%dgsave ) THEN
                        oben_ca2d(ji,jj) = f_benout_ca(ji,jj)
                     ENDIF
                     IF( med_diag%SFR_OCAL%dgsave ) THEN
                        sfr_ocal2d(ji,jj) = f3_omcal(ji,jj,jk)
                     ENDIF
                     IF( med_diag%SFR_OARG%dgsave ) THEN
                        sfr_oarg2d(ji,jj) =  f3_omarg(ji,jj,jk)
                     ENDIF
                     IF( med_diag%LYSO_CA%dgsave ) THEN
                        lyso_ca2d(ji,jj) = f_benout_lyso_ca(ji,jj)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      !! end bathy-1 diags

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !!
               IF( med_diag%RIV_N%dgsave ) THEN
                  rivn2d(ji,jj) = rivn2d(ji,jj) +                            &
                                  (f_riv_loc_n(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%RIV_SI%dgsave ) THEN
                  rivsi2d(ji,jj) = rivsi2d(ji,jj) +                          &
                                   (f_riv_loc_si(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%RIV_C%dgsave ) THEN
                  rivc2d(ji,jj) = rivc2d(ji,jj) +                            &
                                  (f_riv_loc_c(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%RIV_ALK%dgsave ) THEN
                  rivalk2d(ji,jj) = rivalk2d(ji,jj) +                        &
                                    (f_riv_loc_alk(ji,jj) *                  &
                                     fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%DETC%dgsave ) THEN
                  fslowc2d(ji,jj) = fslowc2d(ji,jj) +                        &
                                    (fslowc(ji,jj)  * fse3t(ji,jj,jk))   
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !!
               IF( med_diag%PN_LLOSS%dgsave ) THEN
                  fdpn22d(ji,jj) = fdpn22d(ji,jj) +                          &
                                   (fdpn2(ji,jj)  * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%PD_LLOSS%dgsave ) THEN
                  fdpd22d(ji,jj) = fdpd22d(ji,jj) +                          &
                                   (fdpd2(ji,jj)  * fse3t(ji,jj,jk))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%ZI_LLOSS%dgsave ) THEN
                  fdzmi22d(ji,jj) = fdzmi22d(ji,jj) +                        &
                                    (fdzmi2(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZE_LLOSS%dgsave ) THEN
                  fdzme22d(ji,jj) = fdzme22d(ji,jj) +                        &
                                    (fdzme2(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%ZI_MES_N%dgsave ) THEN
                  zimesn2d(ji,jj) = zimesn2d(ji,jj) +                        &
                                    (xphi * (fgmipn(ji,jj) +                 &
                                             fgmid(ji,jj)) *                 &
                                     fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZI_MES_D%dgsave ) THEN
                  zimesd2d(ji,jj) = zimesd2d(ji,jj) +                        & 
                                    ((1. - xbetan) * finmi(ji,jj) *          &
                                     fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZI_MES_C%dgsave ) THEN
                  zimesc2d(ji,jj) = zimesc2d(ji,jj) +                        &
                                    (xphi * ((xthetapn * fgmipn(ji,jj)) +    &
                                             fgmidc(ji,jj)) *                &
                                             fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZI_MESDC%dgsave ) THEN
                  zimesdc2d(ji,jj) = zimesdc2d(ji,jj) +                      &
                                     ((1. - xbetac) * ficmi(ji,jj) *         &
                                      fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZI_EXCR%dgsave ) THEN
                  ziexcr2d(ji,jj) = ziexcr2d(ji,jj) +                        &
                                    (fmiexcr(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZI_RESP%dgsave ) THEN
                  ziresp2d(ji,jj) = ziresp2d(ji,jj) +                        &
                                    (fmiresp(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZI_GROW%dgsave ) THEN
                  zigrow2d(ji,jj) = zigrow2d(ji,jj) +                        &
                                    (fmigrow(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               IF( med_diag%ZE_MES_N%dgsave ) THEN
                  zemesn2d(ji,jj) = zemesn2d(ji,jj) +                        &
                                    (xphi *                                  &
                                     (fgmepn(ji,jj) + fgmepd(ji,jj) +        &
                                      fgmezmi(ji,jj) + fgmed(ji,jj)) *       &
                                     fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZE_MES_D%dgsave ) THEN
                  zemesd2d(ji,jj) = zemesd2d(ji,jj) +                        &
                                    ((1. - xbetan) * finme(ji,jj) *          &
                                     fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZE_MES_C%dgsave ) THEN
                  zemesc2d(ji,jj) = zemesc2d(ji,jj) +                        & 
                                    (xphi *                                  &
                                     ((xthetapn * fgmepn(ji,jj)) +           &
                                      (xthetapd * fgmepd(ji,jj)) +           &
                                      (xthetazmi * fgmezmi(ji,jj)) +         &
                                      fgmedc(ji,jj)) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZE_MESDC%dgsave ) THEN
                  zemesdc2d(ji,jj) = zemesdc2d(ji,jj) +                      &
                                     ((1. - xbetac) * ficme(ji,jj) *         &
                                      fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZE_EXCR%dgsave ) THEN
                  zeexcr2d(ji,jj) = zeexcr2d(ji,jj) +                        &
                                    (fmeexcr(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZE_RESP%dgsave ) THEN
                  zeresp2d(ji,jj) = zeresp2d(ji,jj) +                        &
                                    (fmeresp(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%ZE_GROW%dgsave ) THEN
                  zegrow2d(ji,jj) = zegrow2d(ji,jj) +                        &
                                    (fmegrow(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
              IF( med_diag%MDETC%dgsave ) THEN
                  mdetc2d(ji,jj) = mdetc2d(ji,jj) +                          &
                                   (fddc(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%GMIDC%dgsave ) THEN
                  gmidc2d(ji,jj) = gmidc2d(ji,jj) +                          &
                                   (fgmidc(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
               IF( med_diag%GMEDC%dgsave ) THEN
                  gmedc2d(ji,jj) = gmedc2d(ji,jj) +                          &
                                   (fgmedc(ji,jj) * fse3t(ji,jj,jk))
               ENDIF
            ENDIF
         ENDDO
      ENDDO
# endif                   

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !!
               !! ** 3D diagnostics
               IF( med_diag%TPP3%dgsave ) THEN
                  tpp3d(ji,jj,jk) = (fprn(ji,jj) * zphn(ji,jj)) +            &
                                    (fprd(ji,jj) * zphd(ji,jj))
                  !CALL iom_put( "TPP3"  , tpp3d )
               ENDIF
               IF( med_diag%TPPD3%dgsave ) THEN
                  tppd3(ji,jj,jk) = (fprd(ji,jj) * zphd(ji,jj))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
                  
               IF( med_diag%REMIN3N%dgsave ) THEN
                  !! remineralisation
                  remin3dn(ji,jj,jk) = fregen(ji,jj) +                       &
                                       (freminn(ji,jj) * fse3t(ji,jj,jk))
                  !CALL iom_put( "REMIN3N"  , remin3dn )
               ENDIF
               !! IF( med_diag%PH3%dgsave ) THEN
               !!     CALL iom_put( "PH3"  , f3_pH )
               !! ENDIF
               !! IF( med_diag%OM_CAL3%dgsave ) THEN
               !!     CALL iom_put( "OM_CAL3"  , f3_omcal )
               !! ENDIF
	       !! 
	       !! AXY (09/11/16): CMIP6 diagnostics
	       IF ( med_diag%DCALC3%dgsave   ) THEN
                  dcalc3(ji,jj,jk) = freminca(ji,jj)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
	       IF ( med_diag%FEDISS3%dgsave  ) THEN
                  fediss3(ji,jj,jk) = ffetop(ji,jj)
               ENDIF
	       IF ( med_diag%FESCAV3%dgsave  ) THEN
                  fescav3(ji,jj,jk) = ffescav(ji,jj)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
	       IF ( med_diag%MIGRAZP3%dgsave ) THEN
                  migrazp3(ji,jj,jk) = fgmipn(ji,jj) * xthetapn
               ENDIF
	       IF ( med_diag%MIGRAZD3%dgsave ) THEN
                  migrazd3(ji,jj,jk) = fgmidc(ji,jj)
               ENDIF
	       IF ( med_diag%MEGRAZP3%dgsave ) THEN
                  megrazp3(ji,jj,jk) = (fgmepn(ji,jj) * xthetapn) +          &
                                       (fgmepd(ji,jj) * xthetapd)
               ENDIF
	       IF ( med_diag%MEGRAZD3%dgsave ) THEN
                  megrazd3(ji,jj,jk) = fgmedc(ji,jj)
               ENDIF
	       IF ( med_diag%MEGRAZZ3%dgsave ) THEN
                   megrazz3(ji,jj,jk) = (fgmezmi(ji,jj) * xthetazmi)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
	       IF ( med_diag%PBSI3%dgsave    ) THEN
                  pbsi3(ji,jj,jk)    = (fprds(ji,jj) * zpds(ji,jj))
               ENDIF
	       IF ( med_diag%PCAL3%dgsave    ) THEN
                  pcal3(ji,jj,jk)    = ftempca(ji,jj)
               ENDIF
	       IF ( med_diag%REMOC3%dgsave   ) THEN
                  remoc3(ji,jj,jk)   = freminc(ji,jj)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
	       IF ( med_diag%PNLIMJ3%dgsave  ) THEN
                  ! pnlimj3(ji,jj,jk)  = fjln(ji,jj)
                  pnlimj3(ji,jj,jk)  = fjlim_pn(ji,jj)
               ENDIF
	       IF ( med_diag%PNLIMN3%dgsave  ) THEN
                  pnlimn3(ji,jj,jk)  = fnln(ji,jj)
               ENDIF
	       IF ( med_diag%PNLIMFE3%dgsave ) THEN
                  pnlimfe3(ji,jj,jk) = ffln2(ji,jj)
               ENDIF
	       IF ( med_diag%PDLIMJ3%dgsave  ) THEN
                  ! pdlimj3(ji,jj,jk)  = fjld(ji,jj)
                  pdlimj3(ji,jj,jk)  = fjlim_pd(ji,jj)
               ENDIF
	       IF ( med_diag%PDLIMN3%dgsave  ) THEN
                  pdlimn3(ji,jj,jk)  = fnld(ji,jj)
               ENDIF
	       IF ( med_diag%PDLIMFE3%dgsave ) THEN
                  pdlimfe3(ji,jj,jk) = ffld(ji,jj)
               ENDIF
	       IF ( med_diag%PDLIMSI3%dgsave ) THEN
                  pdlimsi3(ji,jj,jk) = fsld2(ji,jj)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE bio_med_diag_iomput
#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE bio_med_diag_iomput( )                    ! Empty routine
      IMPLICIT NONE
      WRITE(*,*) 'bio_med_diag_iomput: You should not have seen this print! error?'
   END SUBROUTINE bio_med_diag_iomput
#endif 

   !!======================================================================
END MODULE bio_med_diag_iomput_mod
