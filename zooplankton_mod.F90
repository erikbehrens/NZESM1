MODULE zooplankton_mod
   !!======================================================================
   !!                         ***  MODULE zooplankton_mod  ***
   !! Calculates the zooplankton grazing
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
      
   PUBLIC   zooplankton        ! Called in plankton.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE zooplankton( jk )
      !!------------------------------------------------------------------
      !!                     ***  ROUTINE zooplankton  ***
      !! This called from PLANKTON and calculates the zooplankton 
      !! grazing
      !!------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: fgmed, fgmedc, fgmepd, fgmepds,       &
                                   fgmepn, fgmezmi,                      & 
                                   fgmid, fgmidc, fgmipn,                &
                                   ficme, ficmi, finme, finmi,           &
                                   fmeexcr, fmegrow, fmeresp,            &
                                   fmiexcr, fmigrow, fmiresp,            & 
                                   fsin,                                 &
                                   fzme_i, fzme_o, fzmi_i, fzmi_o,       &
                                   idf, idfval,                          &
                                   zdet, zdtc, zphd, zphn, zzme, zzmi
      USE dom_oce,           ONLY: e3t_0, tmask
#if defined key_vvl
      USE dom_oce,           ONLY: e3t_n
#endif
      USE par_kind,          ONLY: wp
      USE in_out_manager,    ONLY: lwp, numout
      USE par_oce,           ONLY: jpim1, jpjm1
      USE phycst,            ONLY: rsmall
      USE sms_medusa,        ONLY: xbetac, xbetan, xgme, xgmi,           &
                                   xkc, xkme, xkmi, xphi,                &
                                   xpmed, xpmepd, xpmepn, xpmezmi,       &
                                   xpmid, xpmipn,                        &
                                   xthetapd, xthetapn,                   &
                                   xthetazme, xthetazmi

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! Level
      INTEGER, INTENT( in ) :: jk

      INTEGER :: ji, jj

      !! Microzooplankton grazing
      REAL(wp) :: fmi1, fmi
      REAL(wp) :: fstarmi, fmith
      !!
      !! Mesozooplankton grazing
      REAL(wp) :: fme1, fme
      REAL(wp) :: fstarme, fmeth

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then

               !!----------------------------------------------------------
               !! Zooplankton Grazing 
               !! this code supplements the base grazing model with one that
               !! considers the C:N ratio of grazed food and balances this
               !! against the requirements of zooplankton growth; this model 
               !! is derived from that of Anderson & Pondaven (2003)
               !!
               !! The current version of the code assumes a fixed C:N ratio 
               !! for detritus (in contrast to Anderson & Pondaven, 2003), 
               !! though the full equations are retained for future extension
               !!----------------------------------------------------------
               !!
               !!----------------------------------------------------------
               !! Microzooplankton first
               !!----------------------------------------------------------
               !!
               fmi1           = (xkmi * xkmi) + (xpmipn * zphn(ji,jj) *      &
                                                 zphn(ji,jj)) +              &
                                (xpmid * zdet(ji,jj) * zdet(ji,jj))
               fmi            = xgmi * zzmi(ji,jj) / fmi1
               !! grazing on non-diatoms
               fgmipn(ji,jj)  = fmi * xpmipn * zphn(ji,jj) * zphn(ji,jj)
               !! grazing on detrital nitrogen
               fgmid(ji,jj)   = fmi * xpmid  * zdet(ji,jj) * zdet(ji,jj)
# if defined key_roam   
               ! acc            
               fgmidc(ji,jj)  = rsmall
               !! grazing on detrital carbon
               IF ( zdet(ji,jj) .GT. rsmall ) fgmidc(ji,jj)  =               &
                  (zdtc(ji,jj) / (zdet(ji,jj) + tiny(zdet(ji,jj)))) *        &
                  fgmid(ji,jj)
# else
               !! AXY (26/11/08): implicit detrital carbon change
               !! grazing on detrital carbon
               fgmidc(ji,jj)  = xthetad * fgmid(ji,jj)
# endif
# if defined key_debug_medusa
               !! report microzooplankton grazing
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fmi1(',jk,')    = ', fmi1
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! which translates to these incoming N and C fluxes
               finmi(ji,jj)   = (1.0 - xphi) * (fgmipn(ji,jj) + fgmid(ji,jj))
               ficmi(ji,jj)   = (1.0 - xphi) * ((xthetapn * fgmipn(ji,jj)) + &
                                                fgmidc(ji,jj))
               !!
               !! the ideal food C:N ratio for microzooplankton
               !! xbetan = 0.77; xthetaz = 5.625; xbetac = 0.64; xkc = 0.80
               fstarmi = (xbetan * xthetazmi) / (xbetac * xkc)
               !!
               !! process these to determine proportioning of grazed N and C
               !! (since there is no explicit consideration of respiration,
               !! only growth and excretion are calculated here)
               fmith = (ficmi(ji,jj) / (finmi(ji,jj) + tiny(finmi(ji,jj))))
               if (fmith.ge.fstarmi) then
                  fmigrow(ji,jj) = xbetan * finmi(ji,jj)
                  fmiexcr(ji,jj) = 0.0
               else
                  fmigrow(ji,jj) = (xbetac * xkc * ficmi(ji,jj)) / xthetazmi
                  fmiexcr(ji,jj) = ficmi(ji,jj) *                            &
                                   ((xbetan / (fmith + tiny(fmith))) -       &
                                    ((xbetac * xkc) / xthetazmi))
               endif
# if defined key_roam
               fmiresp(ji,jj) = (xbetac * ficmi(ji,jj)) -                    &
                                (xthetazmi * fmigrow(ji,jj))
# endif

# if defined key_debug_medusa
               !! report microzooplankton grazing
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fgmipn(',jk,')  = ', fgmipn(ji,jj)
                  IF (lwp) write (numout,*) 'fgmid(',jk,')   = ', fgmid(ji,jj)
                  IF (lwp) write (numout,*) 'fgmidc(',jk,')  = ', fgmidc(ji,jj)
                  IF (lwp) write (numout,*) 'finmi(',jk,')   = ', finmi(ji,jj)
                  IF (lwp) write (numout,*) 'ficmi(',jk,')   = ', ficmi(ji,jj)
                  IF (lwp) write (numout,*) 'fstarmi(',jk,') = ', fstarmi
                  IF (lwp) write (numout,*) 'fmith(',jk,')   = ', fmith
                  IF (lwp) write (numout,*) 'fmigrow(',jk,') = ', fmigrow(ji,jj)
                  IF (lwp) write (numout,*) 'fmiexcr(',jk,') = ', fmiexcr(ji,jj)
#  if defined key_roam
                  IF (lwp) write (numout,*) 'fmiresp(',jk,') = ', fmiresp(ji,jj)
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
               !! Mesozooplankton second
               !!----------------------------------------------------------
               !!
               fme1           = (xkme * xkme) + (xpmepn * zphn(ji,jj) *       &
                                                 zphn(ji,jj)) +               &
                                (xpmepd * zphd(ji,jj) * zphd(ji,jj)) +        & 
                                (xpmezmi * zzmi(ji,jj) * zzmi(ji,jj)) +       &
                                (xpmed * zdet(ji,jj) * zdet(ji,jj))
               fme            = xgme * zzme(ji,jj) / fme1
               !! grazing on non-diatoms
               fgmepn(ji,jj)  = fme * xpmepn  * zphn(ji,jj) * zphn(ji,jj)
               !! grazing on diatoms
               fgmepd(ji,jj)  = fme * xpmepd  * zphd(ji,jj) * zphd(ji,jj)
               !! grazing on diatom silicon
               fgmepds(ji,jj) = fsin(ji,jj) * fgmepd(ji,jj)
               !! grazing on microzooplankton
               fgmezmi(ji,jj) = fme * xpmezmi * zzmi(ji,jj) * zzmi(ji,jj)
               !! grazing on detrital nitrogen
               fgmed(ji,jj)   = fme * xpmed   * zdet(ji,jj) * zdet(ji,jj)
# if defined key_roam
               !! acc
               fgmedc(ji,jj)  = rsmall
               !! grazing on detrital carbon
               IF ( zdet(ji,jj) .GT. rsmall ) fgmedc(ji,jj)  = (zdtc(ji,jj) / &
                  (zdet(ji,jj) + tiny(zdet(ji,jj)))) * fgmed(ji,jj)
# else
               !! AXY (26/11/08): implicit detrital carbon change
               !! grazing on detrital carbon
               fgmedc(ji,jj)  = xthetad * fgmed(ji,jj)
# endif
               !!
               !! which translates to these incoming N and C fluxes
               finme(ji,jj)   = (1.0 - xphi) *                               &
                                (fgmepn(ji,jj) + fgmepd(ji,jj) +             &
                                 fgmezmi(ji,jj) + fgmed(ji,jj))
               ficme(ji,jj)   = (1.0 - xphi) *                               &
                                ((xthetapn * fgmepn(ji,jj)) +                &
                                (xthetapd * fgmepd(ji,jj)) +                 &
                                (xthetazmi * fgmezmi(ji,jj)) + fgmedc(ji,jj))
# if defined key_debug_medusa
               !! report mesozooplankton grazing
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fme1(',jk,')    = ', fme1
                  IF (lwp) write (numout,*) 'fme(',jk,')     = ', fme
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! the ideal food C:N ratio for mesozooplankton
               !! xbetan = 0.77; xthetaz = 5.625; xbetac = 0.64; xkc = 0.80
               fstarme        = (xbetan * xthetazme) / (xbetac * xkc)
               !!
               !! process these to determine proportioning of grazed N and C
               !! (since there is no explicit consideration of respiration,
               !! only growth and excretion are calculated here)
               fmeth   = (ficme(ji,jj) / (finme(ji,jj) + tiny(finme(ji,jj))))
               if (fmeth.ge.fstarme) then
                  fmegrow(ji,jj) = xbetan * finme(ji,jj)
                  fmeexcr(ji,jj) = 0.0
               else
                  fmegrow(ji,jj) = (xbetac * xkc * ficme(ji,jj)) / xthetazme
                  fmeexcr(ji,jj) = ficme(ji,jj) *                            &
                                   ((xbetan / (fmeth + tiny(fmeth))) -       &
                                    ((xbetac * xkc) / xthetazme))
               endif
# if defined key_roam
               fmeresp(ji,jj) = (xbetac * ficme(ji,jj)) - (xthetazme *       &
                                                           fmegrow(ji,jj))
# endif

# if defined key_debug_medusa
               !! report mesozooplankton grazing
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'fgmepn(',jk,')  = ', fgmepn(ji,jj)
                  IF (lwp) write (numout,*) 'fgmepd(',jk,')  = ', fgmepd(ji,jj)
                  IF (lwp) write (numout,*) 'fgmepds(',jk,') = ', fgmepds(ji,jj)
                  IF (lwp) write (numout,*) 'fgmezmi(',jk,') = ', fgmezmi(ji,jj)
                  IF (lwp) write (numout,*) 'fgmed(',jk,')   = ', fgmed(ji,jj)
                  IF (lwp) write (numout,*) 'fgmedc(',jk,')  = ', fgmedc(ji,jj)
                  IF (lwp) write (numout,*) 'finme(',jk,')   = ', finme(ji,jj)
                  IF (lwp) write (numout,*) 'ficme(',jk,')   = ', ficme(ji,jj)
                  IF (lwp) write (numout,*) 'fstarme(',jk,') = ', fstarme
                  IF (lwp) write (numout,*) 'fmeth(',jk,')   = ', fmeth
                  IF (lwp) write (numout,*) 'fmegrow(',jk,') = ', fmegrow(ji,jj)
                  IF (lwp) write (numout,*) 'fmeexcr(',jk,') = ', fmeexcr(ji,jj)
#  if defined key_roam
                  IF (lwp) write (numout,*) 'fmeresp(',jk,') = ', fmeresp(ji,jj)
#  endif
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               fzmi_i(ji,jj)  = fzmi_i(ji,jj)  + fse3t(ji,jj,jk) *          &
                                ( fgmipn(ji,jj) + fgmid(ji,jj) )
               fzmi_o(ji,jj)  = fzmi_o(ji,jj)  + fse3t(ji,jj,jk) *          &
                                ( fmigrow(ji,jj) +                          &
                                  (xphi * (fgmipn(ji,jj) + fgmid(ji,jj))) + &
                                  fmiexcr(ji,jj) + ((1.0 - xbetan) *        &
                                                    finmi(ji,jj)) )
               fzme_i(ji,jj)  = fzme_i(ji,jj)  + fse3t(ji,jj,jk) *          &
                                ( fgmepn(ji,jj) + fgmepd(ji,jj) +           &
                                  fgmezmi(ji,jj) + fgmed(ji,jj) )
               fzme_o(ji,jj)  = fzme_o(ji,jj)  + fse3t(ji,jj,jk) *          &
                                ( fmegrow(ji,jj) +                          &
                                  (xphi * (fgmepn(ji,jj) + fgmepd(ji,jj) +  &
                                  fgmezmi(ji,jj) + fgmed(ji,jj))) +         &
                                  fmeexcr(ji,jj) + ((1.0 - xbetan) *        &
                                                    finme(ji,jj)) )
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE zooplankton

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE zooplankton( )                    ! Empty routine
      WRITE(*,*) 'zooplankton: You should not have seen this print! error?'
   END SUBROUTINE zooplankton
#endif 

   !!======================================================================
END MODULE zooplankton_mod
