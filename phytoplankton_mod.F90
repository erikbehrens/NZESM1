MODULE phytoplankton_mod
   !!======================================================================
   !!                         ***  MODULE phytoplankton_mod  ***
   !! Calculates the phytoplankton growth
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!   -   ! 2017-08 (A. Yool)            Mean mixed layer chlorophyll
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   phytoplankton        ! Called in plankton.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE phytoplankton( jk )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE phytoplankton  ***
      !! This called from PLANKTON and calculates the phytoplankton
      !! growth.
      !!----------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: fdep1, ffld, ffln2,                   &
                                   fjlim_pd, fjlim_pn,                   &
                                   fnld, fnln,                           &
                                   fprd, fprd_ml, fprds,                 &
                                   fprn, fprn_ml, frd, frn,              &
                                   fsin, fsld, fsld2, fthetad, fthetan,  & 
                                   ftot_det, ftot_dtc, ftot_pd,          &
                                   ftot_pn, ftot_zme, ftot_zmi,          &
                                   fun_Q10, fun_T, idf, idfval,          &
                                   zchd, zchn, zdet, zdin, zdtc,         &
                                   zfer, zpds, zphd, zphn, zsil,         &
                                   zzme, zzmi, fchl_ml
      USE dom_oce,           ONLY: e3t_0, gdepw_0, tmask
#if defined key_vvl
      USE dom_oce,           ONLY: e3t_n, gdepw_n
#endif
      USE in_out_manager,    ONLY: lwp, numout
      USE oce,               ONLY: tsn
      USE par_kind,          ONLY: wp
      USE par_oce,           ONLY: jp_tem, jpi, jpim1, jpj, jpjm1
      USE phycst,            ONLY: rsmall
      USE sms_medusa,        ONLY: jliebig, jphy, jq10,                  &
                                   xald, xaln, xfld, xfln,               &
                                   xnld, xnln, xnsi0, xpar,              &
                                   xsin0, xsld, xthetam, xthetamd, xuif, &
                                   xvpd, xvpn, xxi
      USE zdfmxl,            ONLY: hmld
      USE lbclnk,            ONLY: lbc_lnk

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! Level
      INTEGER, INTENT( in ) :: jk

      INTEGER :: ji, jj

      REAL(wp), DIMENSION(jpi,jpj) :: faln, fchn, fjln
      REAL(wp), DIMENSION(jpi,jpj) :: fald, fchd, fjld
      REAL(wp)                     :: fchn1, fchd1
      !! AXY (03/02/11): add in Liebig terms
      REAL(wp)                     :: fpnlim, fpdlim
      !! AXY (16/07/09): add in Eppley curve functionality
      REAL(wp)                     :: xvpnT,xvpdT
      !! silicon cycle
      REAL(wp)                     :: fnsi

      REAL(wp)                     :: fsin1, fnsi1, fnsi2
      REAL(wp)                     :: fq0

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Chlorophyll calculations
               !!----------------------------------------------------------
               !!
               !! non-diatoms
	       if (zphn(ji,jj).GT.rsmall) then
                  fthetan(ji,jj) = max(tiny(zchn(ji,jj)),                    &
                                       (zchn(ji,jj) * xxi) /                 &
                                       (zphn(ji,jj) + tiny(zphn(ji,jj))))
                  faln(ji,jj)    = xaln * fthetan(ji,jj)
               else
                  fthetan(ji,jj) = 0.
                  faln(ji,jj)    = 0.
               endif
               !!
               !! diatoms
	       if (zphd(ji,jj).GT.rsmall) then
                  fthetad(ji,jj) = max(tiny(zchd(ji,jj)),                   &
                                       (zchd(ji,jj) * xxi) /                &
                                       (zphd(ji,jj) + tiny(zphd(ji,jj))))
                  fald(ji,jj)    = xald * fthetad(ji,jj)
               else
                  fthetad(ji,jj) = 0.
                  fald(ji,jj)    = 0.
               endif

# if defined key_debug_medusa
               !! report biological calculations
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'faln(',jk,') = ', faln(ji,jj)
                  IF (lwp) write (numout,*) 'fald(',jk,') = ', fald(ji,jj)
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Phytoplankton light limitation
               !!----------------------------------------------------------
               !!
               !! It is assumed xpar is the depth-averaged (vertical layer) PAR 
               !! Light limitation (check self-shading) in W/m2
               !!
               !! Note that there is no temperature dependence in phytoplankton
               !! growth rate or any other function. 
               !! In calculation of Chl/Phy ratio tiny(phyto) is introduced to 
               !! avoid NaNs in case of Phy==0.  
               !!
               !! fthetad and fthetan are Chl:C ratio (gChl/gC) in diat and 
               !! non-diat: 
               !! for 1:1 Chl:P ratio (mgChl/mmolN) theta=0.012
               !!
               !! AXY (16/07/09)
               !! temperature for new Eppley style phytoplankton growth
               fun_T(ji,jj)   = 1.066**(1.0 * tsn(ji,jj,jk,jp_tem))
               !! AXY (16/05/11): add in new Q10 (1.5, not 2.0) for
               !!                 phytoplankton growth; remin. unaffected
               fun_Q10(ji,jj) = jq10**((tsn(ji,jj,jk,jp_tem) - 0.0) / 10.0)
               if (jphy.eq.1) then
                  xvpnT = xvpn * fun_T(ji,jj)
                  xvpdT = xvpd * fun_T(ji,jj)
               elseif (jphy.eq.2) then
                  xvpnT = xvpn * fun_Q10(ji,jj)
                  xvpdT = xvpd * fun_Q10(ji,jj)
               else
                  xvpnT = xvpn
                  xvpdT = xvpd
               endif
               !!
               !! non-diatoms
               fchn1 = (xvpnT * xvpnT) +                                     &
                       (faln(ji,jj) * faln(ji,jj) * xpar(ji,jj,jk) *         &
                        xpar(ji,jj,jk))
               if (fchn1.GT.rsmall) then
                  fchn(ji,jj) = xvpnT / (sqrt(fchn1) + tiny(fchn1))
               else
                  fchn(ji,jj) = 0.
               endif
               !! non-diatom J term
               fjln(ji,jj)     = fchn(ji,jj) * faln(ji,jj) * xpar(ji,jj,jk)
               fjlim_pn(ji,jj) = fjln(ji,jj) / xvpnT
               !!
               !! diatoms
               fchd1 = (xvpdT * xvpdT) +                                     &
                       (fald(ji,jj) * fald(ji,jj) * xpar(ji,jj,jk) *         &
                        xpar(ji,jj,jk))
               if (fchd1.GT.rsmall) then
                  fchd(ji,jj) = xvpdT / (sqrt(fchd1) + tiny(fchd1))
               else
                  fchd(ji,jj) = 0.
               endif
               !! diatom J term
               fjld(ji,jj)    = fchd(ji,jj) * fald(ji,jj) * xpar(ji,jj,jk)
               fjlim_pd(ji,jj) = fjld(ji,jj) / xvpdT
      
# if defined key_debug_medusa
               !! report phytoplankton light limitation
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fchn(',jk,') = ', fchn(ji,jj)
                  IF (lwp) write (numout,*) 'fchd(',jk,') = ', fchd(ji,jj)
                  IF (lwp) write (numout,*) 'fjln(',jk,') = ', fjln(ji,jj)
                  IF (lwp) write (numout,*) 'fjld(',jk,') = ', fjld(ji,jj)
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Phytoplankton nutrient limitation
               !!----------------------------------------------------------
               !!
               !! non-diatoms (N, Fe).
               !! non-diatom Qn term
               fnln(ji,jj)  = zdin(ji,jj) / (zdin(ji,jj) + xnln)
               !! non-diatom Qf term
               ffln2(ji,jj) = zfer(ji,jj) / (zfer(ji,jj) + xfln)
               !!
               !! diatoms (N, Si, Fe).
               !! diatom Qn term
               fnld(ji,jj) = zdin(ji,jj) / (zdin(ji,jj) + xnld)
               !! diatom Qs term
               fsld(ji,jj) = zsil(ji,jj) / (zsil(ji,jj) + xsld)
               !! diatom Qf term
               ffld(ji,jj) = zfer(ji,jj) / (zfer(ji,jj) + xfld)

# if defined key_debug_medusa
               !! report phytoplankton nutrient limitation
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fnln(',jk,') = ', fnln(ji,jj)
                  IF (lwp) write (numout,*) 'fnld(',jk,') = ', fnld(ji,jj)
                  IF (lwp) write (numout,*) 'ffln2(',jk,') = ', ffln2(ji,jj)
                  IF (lwp) write (numout,*) 'ffld(',jk,') = ', ffld(ji,jj)
                  IF (lwp) write (numout,*) 'fsld(',jk,') = ', fsld(ji,jj)
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Primary production (non-diatoms)
               !! (note: still needs multiplying by phytoplankton 
               !! concentration)
               !!----------------------------------------------------------
               !!
               if (jliebig .eq. 0) then
                  !! multiplicative nutrient limitation
                  fpnlim = fnln(ji,jj) * ffln2(ji,jj)
               elseif (jliebig .eq. 1) then
                  !! Liebig Law (= most limiting) nutrient limitation
                  fpnlim = min(fnln(ji,jj), ffln2(ji,jj))
               endif
               fprn(ji,jj) = fjln(ji,jj) * fpnlim
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Primary production (diatoms)
               !! (note: still needs multiplying by phytoplankton 
               !! concentration)
               !!
               !! Production here is split between nitrogen production and 
               !! that of silicon; depending upon the "intracellular" ratio 
               !! of Si:N, model diatoms will uptake nitrogen/silicon 
               !! differentially; this borrows from the diatom model of 
               !! Mongin et al. (2006)
               !!----------------------------------------------------------
               !!
               if (jliebig .eq. 0) then
                  !! multiplicative nutrient limitation
                  fpdlim = fnld(ji,jj) * ffld(ji,jj)
               elseif (jliebig .eq. 1) then
                  !! Liebig Law (= most limiting) nutrient limitation
                  fpdlim = min(fnld(ji,jj), ffld(ji,jj))
               endif
               !!
	       if (zphd(ji,jj).GT.rsmall .AND. zpds(ji,jj).GT.rsmall) then
                  !! "intracellular" elemental ratios
                  ! fsin(ji,jj)  = zpds(ji,jj) / (zphd(ji,jj) +              &
                  !                               tiny(zphd(ji,jj)))
                  ! fnsi         = zphd(ji,jj) / (zpds(ji,jj) +              &
                  !                               tiny(zpds(ji,jj)))
                  fsin(ji,jj) = 0.0
                  IF( zphd(ji,jj) .GT. rsmall) fsin(ji,jj)  = zpds(ji,jj) /  &
                                                              zphd(ji,jj)
                  fnsi = 0.0
                  IF( zpds(ji,jj) .GT. rsmall) fnsi  = zphd(ji,jj) /         &
                                                       zpds(ji,jj)
                  !! AXY (23/02/10): these next variables derive from 
                  !! Mongin et al. (2003)
                  fsin1 = 3.0 * xsin0 !! = 0.6
                  fnsi1 = 1.0 / fsin1 !! = 1.667
                  fnsi2 = 1.0 / xsin0 !! = 5.0
                  !!
                  !! conditionalities based on ratios
                  !! nitrogen (and iron and carbon)
                  if (fsin(ji,jj).le.xsin0) then
                     fprd(ji,jj)  = 0.0
                     fsld2(ji,jj) = 0.0
                  elseif (fsin(ji,jj).lt.fsin1) then
                     fprd(ji,jj)  = xuif * ((fsin(ji,jj) - xsin0) /          &
                                            (fsin(ji,jj) +                   &
                                             tiny(fsin(ji,jj)))) *           &
                                    (fjld(ji,jj) * fpdlim)
                     fsld2(ji,jj) = xuif * ((fsin(ji,jj) - xsin0) /          &
                                            (fsin(ji,jj) +                   &
                                             tiny(fsin(ji,jj))))
                  elseif (fsin(ji,jj).ge.fsin1) then
                     fprd(ji,jj)  = (fjld(ji,jj) * fpdlim)
                     fsld2(ji,jj) = 1.0
                  endif
                  !!
                  !! silicon
                  if (fsin(ji,jj).lt.fnsi1) then
                     fprds(ji,jj) = (fjld(ji,jj) * fsld(ji,jj))
                  elseif (fsin(ji,jj).lt.fnsi2) then
                     fprds(ji,jj) = xuif * ((fnsi - xnsi0) /          &
                                            (fnsi + tiny(fnsi))) *           &
                                    (fjld(ji,jj) * fsld(ji,jj))
                  else
                     fprds(ji,jj) = 0.0
                  endif     
               else
                  fsin(ji,jj)  = 0.0
                  fnsi         = 0.0
                  fprd(ji,jj)  = 0.0
                  fsld2(ji,jj) = 0.0
                  fprds(ji,jj) = 0.0
               endif

# if defined key_debug_medusa
               !! report phytoplankton growth (including diatom silicon 
               !! submodel)
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fsin(',jk,')   = ', fsin(ji,jj)
                  IF (lwp) write (numout,*) 'fnsi(',jk,')   = ', fnsi
                  IF (lwp) write (numout,*) 'fsld2(',jk,')  = ', fsld2(ji,jj)
                  IF (lwp) write (numout,*) 'fprn(',jk,')   = ', fprn(ji,jj)
                  IF (lwp) write (numout,*) 'fprd(',jk,')   = ', fprd(ji,jj)
                  IF (lwp) write (numout,*) 'fprds(',jk,')  = ', fprds(ji,jj)
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Mixed layer primary production
               !! this block calculates the amount of primary production 
               !! that occurs within the upper mixed layer; this allows the 
               !! separate diagnosis of "sub-surface" primary production; it 
               !! does assume that short-term variability in mixed layer 
               !! depth doesn't mess with things though
               !!----------------------------------------------------------
               !!
               if (fdep1(ji,jj).le.hmld(ji,jj)) then
                  !! this level is entirely in the mixed layer
                  fq0 = 1.0
               elseif (fsdepw(ji,jj,jk).ge.hmld(ji,jj)) then
                  !! this level is entirely below the mixed layer
                  fq0 = 0.0
               else
                  !! this level straddles the mixed layer
                  fq0 = (hmld(ji,jj) - fsdepw(ji,jj,jk)) / fse3t(ji,jj,jk)
               endif
               !!
               fprn_ml(ji,jj) = fprn_ml(ji,jj) + (fprn(ji,jj) * zphn(ji,jj) * &
                                                  fse3t(ji,jj,jk) * fq0)
               fprd_ml(ji,jj) = fprd_ml(ji,jj) + (fprd(ji,jj) * zphd(ji,jj) * &
                                                  fse3t(ji,jj,jk) * fq0)
	       !! AXY (16/08/17)
	       fchl_ml(ji,jj) = fchl_ml(ji,jj) + ((zchn(ji,jj) + zchd(ji,jj)) * &
	                                          (fse3t(ji,jj,jk) * fq0) / hmld(ji,jj))
            ENDIF
         ENDDO
      ENDDO
      CALL lbc_lnk(fchl_ml(:,:),'T',1. )

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Vertical Integral --
               !!----------------------------------------------------------
               !! vertical integral non-diatom phytoplankton
               ftot_pn(ji,jj)  = ftot_pn(ji,jj)  + (zphn(ji,jj) *             &
                                                    fse3t(ji,jj,jk))
               !! vertical integral diatom phytoplankton
               ftot_pd(ji,jj)  = ftot_pd(ji,jj)  + (zphd(ji,jj) *             &
                                                    fse3t(ji,jj,jk))
               !! vertical integral microzooplankton
               ftot_zmi(ji,jj) = ftot_zmi(ji,jj) + (zzmi(ji,jj) *             &
                                                    fse3t(ji,jj,jk))
               !! vertical integral mesozooplankton
               ftot_zme(ji,jj) = ftot_zme(ji,jj) + (zzme(ji,jj) *             &
                                                    fse3t(ji,jj,jk))
               !! vertical integral slow detritus, nitrogen
               ftot_det(ji,jj) = ftot_det(ji,jj) + (zdet(ji,jj) *             &
                                                    fse3t(ji,jj,jk))
               !! vertical integral slow detritus, carbon
               ftot_dtc(ji,jj) = ftot_dtc(ji,jj) + (zdtc(ji,jj) *             &
                                                    fse3t(ji,jj,jk))
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! More chlorophyll calculations
               !!----------------------------------------------------------
               !!
               !! frn(ji,jj) = (xthetam / fthetan(ji,jj)) *                   &
               !!              (fprn(ji,jj) / (fthetan(ji,jj) * xpar(ji,jj,jk)))
               !! frd(ji,jj) = (xthetam / fthetad(ji,jj)) *                   &
               !!              (fprd(ji,jj) / (fthetad(ji,jj) * xpar(ji,jj,jk)))
               frn(ji,jj) = (xthetam * fchn(ji,jj) * fnln(ji,jj) *            &
                             ffln2(ji,jj)) / (fthetan(ji,jj) +                &
                                             tiny(fthetan(ji,jj)))
               !! AXY (12/05/09): there's potentially a problem here; fsld, 
               !!   silicic acid limitation, is used in the following line 
               !!   to regulate chlorophyll growth in a manner that is 
               !!   inconsistent with its use in the regulation of biomass 
               !!   growth; the Mongin term term used in growth is more
               !!   complex than the simple multiplicative function used
               !!   below
               !! frd(ji,jj) = (xthetam * fchd(ji,jj) * fnld(ji,jj) *        &
               !!               ffld(ji,jj) * fsld(ji,jj)) /                 &
               !!               (fthetad(ji,jj) + tiny(fthetad(ji,jj)))
               !! AXY (12/05/09): this replacement line uses the new 
               !!   variable, fsld2, to regulate chlorophyll growth
               frd(ji,jj) = (xthetamd * fchd(ji,jj) * fnld(ji,jj) *          &
                             ffld(ji,jj) * fsld2(ji,jj)) /                   &
                             (fthetad(ji,jj) + tiny(fthetad(ji,jj)))

# if defined key_debug_medusa
               !! report chlorophyll calculations
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fthetan(',jk,') = ', fthetan(ji,jj)
                  IF (lwp) write (numout,*) 'fthetad(',jk,') = ', fthetad(ji,jj)
                  IF (lwp) write (numout,*) 'frn(',jk,')     = ', frn(ji,jj)
                  IF (lwp) write (numout,*) 'frd(',jk,')     = ', frd(ji,jj)
               endif
# endif

            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE phytoplankton

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE phytoplankton( )                    ! Empty routine
      WRITE(*,*) 'phytoplankton: You should not have seen this print! error?'
   END SUBROUTINE phytoplankton
#endif 

   !!======================================================================
END MODULE phytoplankton_mod
