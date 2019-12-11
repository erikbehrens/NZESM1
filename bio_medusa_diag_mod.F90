MODULE bio_medusa_diag_mod
   !!======================================================================
   !!                         ***  MODULE bio_medusa_diag_mod  ***
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
      
   PUBLIC   bio_medusa_diag        ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE bio_medusa_diag( jk )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE bio_medusa_diag  ***
      !! This called from TRC_BIO_MEDUSA and calculates diagnostics
      !!-------------------------------------------------------------------
      USE bio_med_diag_iomput_mod,  ONLY: bio_med_diag_iomput
      USE bio_medusa_mod
      USE dom_oce,                  ONLY: e3t_0, gdepw_0, tmask
# if defined key_vvl
      USE dom_oce,                  ONLY: e3t_n, gdepw_n
# endif
      USE in_out_manager,           ONLY: lwp, numout
      USE iom,                      ONLY: lk_iomput
      USE par_oce,                  ONLY: jpim1, jpjm1
      USE sms_medusa,               ONLY: xrfn, xthetapd, xthetapn,      &
                                          xthetazme, xthetazmi
      USE trc,                      ONLY: med_diag 
# if defined key_roam
      USE trcoxy_medusa,            ONLY: oxy_sato
# endif

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! level
      INTEGER, INTENT( in ) :: jk

      !! Loop avariables
      INTEGER :: ji, jj, jn

# if defined key_trc_diabio
      !!==========================================================
      !! LOCAL GRID CELL DIAGNOSTICS
      !!==========================================================
      !!
      !!----------------------------------------------------------
      !! Full diagnostics key_trc_diabio:
      !! LOBSTER and PISCES support full diagnistics option 
      !! key_trc_diabio which gives an option of FULL output of 
      !! biological sourses and sinks. I cannot see any reason 
      !! for doing this. If needed, it can be done as shown
      !! below.
      !!----------------------------------------------------------
      !!
      IF(lwp) WRITE(numout,*) ' MEDUSA does not support key_trc_diabio'
# endif

      !! The section below, down to calculation of zo2min, was moved 
      !! from before the call to AIR_SEA in trcbio_medusa.F90 - marc 9/5/17 
      IF( lk_iomput ) THEN
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  !! sum tracers for inventory checks
                  IF ( med_diag%INVTN%dgsave )   THEN
                     ftot_n(ji,jj)  = ftot_n(ji,jj) +                        &
                        (fse3t(ji,jj,jk) * (zphn(ji,jj) + zphd(ji,jj) +      &
                                            zzmi(ji,jj) + zzme(ji,jj) +      &
                                            zdet(ji,jj) + zdin(ji,jj)))
                  ENDIF
                  IF ( med_diag%INVTSI%dgsave )  THEN
                     ftot_si(ji,jj) = ftot_si(ji,jj) +                       & 
                       (fse3t(ji,jj,jk) * (zpds(ji,jj) + zsil(ji,jj)))
                  ENDIF
                  IF ( med_diag%INVTFE%dgsave )  THEN
                     ftot_fe(ji,jj) = ftot_fe(ji,jj) +                       & 
                        (fse3t(ji,jj,jk) * (xrfn *                           &
                                            (zphn(ji,jj) + zphd(ji,jj) +     &
                                             zzmi(ji,jj) + zzme(ji,jj) +     &
                                             zdet(ji,jj)) +                  &
                                            zfer(ji,jj)))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

# if defined key_roam
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  IF ( med_diag%INVTC%dgsave )  THEN
                     ftot_c(ji,jj)  = ftot_c(ji,jj) +                        & 
                        (fse3t(ji,jj,jk) * ((xthetapn * zphn(ji,jj)) +       &
                                            (xthetapd * zphd(ji,jj)) +       &
                                            (xthetazmi * zzmi(ji,jj)) +      &
                                            (xthetazme * zzme(ji,jj)) +      &
                                            zdtc(ji,jj) + zdic(ji,jj)))
                  ENDIF
                  IF ( med_diag%INVTALK%dgsave ) THEN
                     ftot_a(ji,jj)  = ftot_a(ji,jj) + (fse3t(ji,jj,jk) *     &
                                                       zalk(ji,jj))
                  ENDIF
                  IF ( med_diag%INVTO2%dgsave )  THEN
                     ftot_o2(ji,jj) = ftot_o2(ji,jj) + (fse3t(ji,jj,jk) *    &
                                                        zoxy(ji,jj))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  IF ( med_diag%INVTC%dgsave )  THEN
                     !!
                     !! AXY (10/11/16): CMIP6 diagnostics
                     IF ( med_diag%INTDISSIC%dgsave ) THEN
                        intdissic(ji,jj) = intdissic(ji,jj) +                &
                                           (fse3t(ji,jj,jk) * zdic(ji,jj))
                     ENDIF
                     IF ( med_diag%INTDISSIN%dgsave ) THEN
                        intdissin(ji,jj) = intdissin(ji,jj) +                &
                                           (fse3t(ji,jj,jk) * zdin(ji,jj))
                     ENDIF
                     IF ( med_diag%INTDISSISI%dgsave ) THEN
                        intdissisi(ji,jj) = intdissisi(ji,jj) +              &
                                            (fse3t(ji,jj,jk) * zsil(ji,jj))
                     ENDIF
                     IF ( med_diag%INTTALK%dgsave ) THEN
                        inttalk(ji,jj) = inttalk(ji,jj) +                    &
                                         (fse3t(ji,jj,jk) * zalk(ji,jj))
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  IF ( med_diag%O2min%dgsave ) THEN
                     if ( zoxy(ji,jj) < o2min(ji,jj) ) then
                        o2min(ji,jj)  = zoxy(ji,jj)
                        IF ( med_diag%ZO2min%dgsave ) THEN
                           !! layer midpoint
                           zo2min(ji,jj) = (fsdepw(ji,jj,jk) +               &
                                            fdep1(ji,jj)) / 2.0
                        ENDIF
                     endif
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
# endif
      ENDIF

# if defined key_roam
      !! This section is moved from just below CALL to AIR_SEA in
      !! trcbio_medusa.F90 - marc 9/5/17
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then

               !! AXY (11/11/16): CMIP6 oxygen saturation 3D diagnostic
               IF ( med_diag%O2SAT3%dgsave ) THEN
! Remove f_o2sat3 - marc 9/5/17
!                  call oxy_sato( ztmp(ji,jj), zsal(ji,jj), f_o2sat3 )
!                  o2sat3(ji, jj, jk) = f_o2sat3
                  call oxy_sato( ztmp(ji,jj), zsal(ji,jj),                   &
                                 o2sat3(ji,jj,jk) )
               ENDIF
            ENDIF
         ENDDO
      ENDDO
# endif

      !!---------------------------------------------------------------
      !! Calculates the diagnostics used with iom_put
      !!---------------------------------------------------------------
      CALL bio_med_diag_iomput( jk )

   END SUBROUTINE bio_medusa_diag

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE bio_medusa_diag( )                    ! Empty routine
      WRITE(*,*) 'bio_medusa_diag: You should not have seen this print! error?'
   END SUBROUTINE bio_medusa_diag
#endif 

   !!======================================================================
END MODULE bio_medusa_diag_mod
