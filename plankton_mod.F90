MODULE plankton_mod
   !!======================================================================
   !!                         ***  MODULE plankton_mod  ***
   !! Calculate the carbon chemistry for the whole ocean
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
      
   PUBLIC   plankton        ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE plankton( jk )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE plankton  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - Calculates phytoplankton growth
      !!  - Zooplankton grazing
      !!  - Plankton losses
      !!-------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: fdpd, fdpd2, fdpds, fdpds2,             &
                                   fdpn, fdpn2, fdzme, fdzme2,             &
                                   fdzmi, fdzmi2, fsdiss, fsin,            &
                                   zphd, zphn, zpds, zzme, zzmi
      USE dom_oce,           ONLY: tmask
      USE par_oce,           ONLY: jpim1, jpjm1
      USE phytoplankton_mod, ONLY: phytoplankton
      USE sms_medusa,        ONLY: jmpd, jmpn, jmzme, jmzmi,               &
                                   xkphd, xkphn, xkzme, xkzmi,             &
                                   xmetapd, xmetapn, xmetazme, xmetazmi,   &
                                   xmpd, xmpn, xmzme, xmzmi, xsdiss
      USE zooplankton_mod,   ONLY: zooplankton

      !! Level
      INTEGER, INTENT( in ) :: jk

      INTEGER :: ji, jj

      !!-------------------------------------------------------------------
      !! Calculate phytoplankton growth
      !!-------------------------------------------------------------------
      CALL phytoplankton( jk )

      !!-------------------------------------------------------------------
      !! Calculate zooplankton grazing
      !!-------------------------------------------------------------------
      CALL zooplankton( jk )

      !!-------------------------------------------------------------------
      !! Miscellaneous plankton losses
      !!-------------------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Plankton metabolic losses
               !! Linear loss processes assumed to be metabolic in origin
               !!----------------------------------------------------------
               !!
               fdpn2(ji,jj)  = xmetapn  * zphn(ji,jj)
               fdpd2(ji,jj)  = xmetapd  * zphd(ji,jj)
               fdpds2(ji,jj) = xmetapd  * zpds(ji,jj)
               fdzmi2(ji,jj) = xmetazmi * zzmi(ji,jj)
               fdzme2(ji,jj) = xmetazme * zzme(ji,jj)
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------
               !! Plankton mortality losses
               !! EKP (26/02/09): phytoplankton hyperbolic mortality term 
               !! introduced 
               !! to improve performance in gyres
               !!----------------------------------------------------------
               !!
               !! non-diatom phytoplankton
               !! linear
               if (jmpn.eq.1) fdpn(ji,jj) = xmpn * zphn(ji,jj)
               !! quadratic
               if (jmpn.eq.2) fdpn(ji,jj) = xmpn * zphn(ji,jj) * zphn(ji,jj)
               !! hyperbolic
               if (jmpn.eq.3) fdpn(ji,jj) = xmpn * zphn(ji,jj) *             &
                                            (zphn(ji,jj) /                   &
                                             (xkphn + zphn(ji,jj)))
               !! sigmoid
               if (jmpn.eq.4) fdpn(ji,jj) = xmpn * zphn(ji,jj) *             &
                                            ((zphn(ji,jj) * zphn(ji,jj)) /   &
                                             (xkphn + (zphn(ji,jj) *         &
                                                       zphn(ji,jj))))
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! diatom phytoplankton
               !! linear
               if (jmpd.eq.1) fdpd(ji,jj) = xmpd * zphd(ji,jj)
               !! quadratic
               if (jmpd.eq.2) fdpd(ji,jj) = xmpd * zphd(ji,jj) * zphd(ji,jj)
               !! hyperbolic
               if (jmpd.eq.3) fdpd(ji,jj) = xmpd * zphd(ji,jj) *             &
                                            (zphd(ji,jj) / (xkphd +          &
                                                            zphd(ji,jj)))
               !! sigmoid
               if (jmpd.eq.4) fdpd(ji,jj) = xmpd * zphd(ji,jj) *             &
                                            ((zphd(ji,jj) * zphd(ji,jj)) /   &
                                             (xkphd + (zphd(ji,jj) *         &
                                                       zphd(ji,jj))))
               fdpds(ji,jj) = fdpd(ji,jj) * fsin(ji,jj)
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! microzooplankton
               !! linear
               if (jmzmi.eq.1) fdzmi(ji,jj) = xmzmi * zzmi(ji,jj)
               !! quadratic
               if (jmzmi.eq.2) fdzmi(ji,jj) = xmzmi * zzmi(ji,jj) *          &
                                              zzmi(ji,jj)
               !! hyperbolic
               if (jmzmi.eq.3) fdzmi(ji,jj) = xmzmi * zzmi(ji,jj) *          &
                                              (zzmi(ji,jj) / (xkzmi +        &
                                                              zzmi(ji,jj)))
               !! sigmoid
               if (jmzmi.eq.4) fdzmi(ji,jj) = xmzmi * zzmi(ji,jj) * &
                  ((zzmi(ji,jj) * zzmi(ji,jj)) / (xkzmi + (zzmi(ji,jj) *     &
                                                           zzmi(ji,jj))))
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! mesozooplankton
               !! linear
               if (jmzme.eq.1) fdzme(ji,jj) = xmzme * zzme(ji,jj)
               !! quadratic
               if (jmzme.eq.2) fdzme(ji,jj) = xmzme * zzme(ji,jj) *          &
                                              zzme(ji,jj)
               !! hyperbolic
               if (jmzme.eq.3) fdzme(ji,jj) = xmzme * zzme(ji,jj) *          &
                                              (zzme(ji,jj) / (xkzme +        &
                                                              zzme(ji,jj)))
               !! sigmoid
               if (jmzme.eq.4) fdzme(ji,jj) = xmzme * zzme(ji,jj) *          &
                                              ((zzme(ji,jj) * zzme(ji,jj)) / &
                                               (xkzme + (zzme(ji,jj) *       &
                                                         zzme(ji,jj))))
            ENDIF
         ENDDO
      ENDDO

      !! diatom frustule dissolution. This section is moved from just
      !! below CALL to iron_chem_scav in trcbio_medusa.F90 - marc 9/5/17
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               fsdiss(ji,jj)  = xsdiss * zpds(ji,jj)
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE plankton

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE plankton( )                    ! Empty routine
      WRITE(*,*) 'plankton: You should not have seen this print! error?'
   END SUBROUTINE plankton
#endif 

   !!======================================================================
END MODULE plankton_mod
