MODULE detritus_mod
   !!======================================================================
   !!                         ***  MODULE detritus_mod  ***
   !! Calculates detritus processes and fast-sinking detritus
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!   -   ! 2017-08 (A. Yool)            Revise slow-sinking of detritus
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   detritus        ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE detritus( jk, iball  )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE detritus  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - Calculates detritus processes
      !!  - Fast-sinking detritus 
      !!-------------------------------------------------------------------
      USE bio_medusa_mod,         ONLY: f_sbenin_c, f_sbenin_fe,           &
                                        f_sbenin_n, fdd,                   &
                                        idf, idfval,                       &   
                                        fslowsink,                         &
                                        fslowgain, fslowloss,              &
# if defined key_roam
                                        fslowsinkc,                        &
                                        fslowgainc, fslowlossc,            &
                                        fddc,                              &
# endif
                                        fun_T, fun_Q10, zdet, zdtc
      USE detritus_fast_sink_mod, ONLY: detritus_fast_sink
      USE dom_oce,                ONLY: mbathy, e3t_0, gphit, tmask
# if defined key_vvl
      USE dom_oce,                ONLY: e3t_n
# endif
      USE in_out_manager,         ONLY: lwp, numout
      USE par_oce,                ONLY: jpim1, jpjm1
      USE sms_medusa,             ONLY: jmd, jorgben, jsfd, vsed,          &
                                        xrfn, xmd, xmdc, xthetad

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! Level
      INTEGER, INTENT( in ) :: jk
      !! Fast detritus ballast scheme (0 = no; 1 = yes)
      INTEGER, INTENT( in ) :: iball

      INTEGER :: ji, jj

      !!------------------------------------------------------------------
      !! Detritus remineralisation
      !! Constant or temperature-dependent
      !!------------------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then
               !!
               if (jmd.eq.1) then
                  !! temperature-dependent
                  fdd(ji,jj)  = xmd  * fun_T(ji,jj) * zdet(ji,jj)
# if defined key_roam
                  fddc(ji,jj) = xmdc * fun_T(ji,jj) * zdtc(ji,jj)
# endif
               elseif (jmd.eq.2) then
                  !! AXY (16/05/13): add in Q10-based parameterisation 
                  !! (def in nmlst)
                  !! temperature-dependent
                  fdd(ji,jj)  = xmd  * fun_Q10(ji,jj) * zdet(ji,jj)
#if defined key_roam
                  fddc(ji,jj) = xmdc * fun_Q10(ji,jj) * zdtc(ji,jj)
#endif
               else
                  !! temperature-independent
                  fdd(ji,jj)  = xmd  * zdet(ji,jj)
# if defined key_roam
                  fddc(ji,jj) = xmdc * zdtc(ji,jj)
# endif
               endif
               !!
               !! AXY (22/07/09): accelerate detrital remineralisation 
               !! in the bottom box
               if ((jk.eq.mbathy(ji,jj)) .and. jsfd.eq.1) then
                  fdd(ji,jj)  = 1.0  * zdet(ji,jj)
# if defined key_roam
                  fddc(ji,jj) = 1.0  * zdtc(ji,jj)
# endif
               endif
               
# if defined key_debug_medusa
               !! report plankton mortality and remineralisation
               if (idf.eq.1.AND.idfval.eq.1) then
                  IF (lwp) write (numout,*) '------------------------------'
! I've removed the lines below, because the variables are not in this
! routine. If these debug prints need to stay, they should probably be
! moved - marc 27/4/17
!                  IF (lwp) write (numout,*) 'fdpn2(',jk,') = ', fdpn2(ji,jj)
!                  IF (lwp) write (numout,*) 'fdpd2(',jk,') = ', fdpd2(ji,jj)
!                  IF (lwp) write (numout,*) 'fdpds2(',jk,')= ', fdpds2(ji,jj)
!                  IF (lwp) write (numout,*) 'fdzmi2(',jk,')= ', fdzmi2(ji,jj)
!                  IF (lwp) write (numout,*) 'fdzme2(',jk,')= ', fdzme2(ji,jj)
!                  IF (lwp) write (numout,*) 'fdpn(',jk,')  = ', fdpn(ji,jj)
!                  IF (lwp) write (numout,*) 'fdpd(',jk,')  = ', fdpd(ji,jj)
!                  IF (lwp) write (numout,*) 'fdpds(',jk,') = ', fdpds(ji,jj)
!                  IF (lwp) write (numout,*) 'fdzmi(',jk,') = ', fdzmi(ji,jj)
!                  IF (lwp) write (numout,*) 'fdzme(',jk,') = ', fdzme(ji,jj)
                  IF (lwp) write (numout,*) 'fdd(',jk,')   = ', fdd(ji,jj)
#  if defined key_roam
                  IF (lwp) write (numout,*) 'fddc(',jk,')  = ', fddc(ji,jj)
#  endif
               endif
# endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!----------------------------------------------------------------------
               !! Detritus sinking (AXY, 08/08/18)
	       !! Replaces slow-sinking done in trcsed_medusa.F90
               !!
               !! Uses the fslowsink variable to carry slow-sinking detritus from one
               !! grid level to the next, variable fslowgain to "add" detritus sinking
               !! from above and variable fslowloss to "subtract" detritus sinking out
               !! to below; these variables appear in the differential equations of
               !! detrital nitrogen and carbon below
               !!----------------------------------------------------------------------
               !!
               fslowgain(ji,jj)  = fslowsink(ji,jj) / fse3t(ji,jj,jk)                  ! = mmol N / m3 / d
               if (jk.lt.mbathy(ji,jj)) then
                  fslowloss(ji,jj)  = (zdet(ji,jj) * vsed * 86400.) / fse3t(ji,jj,jk)  ! = mmol N / m3 / d
               else
                  fslowloss(ji,jj)  = 0.                                               ! = mmol N / m3 / d
               endif
               fslowsink(ji,jj) = fslowloss(ji,jj) * fse3t(ji,jj,jk)                   ! = mmol N / m2 / d
               !!
#  if defined key_roam
               fslowgainc(ji,jj) = fslowsinkc(ji,jj) / fse3t(ji,jj,jk)                 ! = mmol C / m3 / d
               if (jk.lt.mbathy(ji,jj)) then
                  fslowlossc(ji,jj) = (zdtc(ji,jj) * vsed * 86400.) / fse3t(ji,jj,jk)  ! = mmol C / m3 / d
               else
                  fslowlossc(ji,jj) = 0.                                               ! = mmol C / m3 / d
               endif
               fslowsinkc(ji,jj) = fslowlossc(ji,jj) * fse3t(ji,jj,jk)                 ! = mmol C / m2 / d
#  endif
            ENDIF
         ENDDO
      ENDDO

      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            if (tmask(ji,jj,jk) == 1) then
               !!---------------------------------------------------------
               !! Detritus addition to benthos
               !! If activated, slow detritus in the bottom box will enter
               !! the benthic pool
               !!---------------------------------------------------------
               !!
               if ((jk.eq.mbathy(ji,jj)) .and. jorgben.eq.1) then
                  !! this is the BOTTOM OCEAN BOX -> into the benthic pool!
                  !!
                  f_sbenin_n(ji,jj)  = (zdet(ji,jj) * vsed * 86400.)
                  f_sbenin_fe(ji,jj) = (zdet(ji,jj) * vsed * 86400. * xrfn)
# if defined key_roam
                  f_sbenin_c(ji,jj)  = (zdtc(ji,jj) * vsed * 86400.)
# else
                  f_sbenin_c(ji,jj)  = (zdet(ji,jj) * vsed * 86400. * xthetad)
# endif
               endif

            ENDIF
         ENDDO
      ENDDO

      !!------------------------------------------------------------------
      !! Fast-sinking detritus
      !!------------------------------------------------------------------
      CALL detritus_fast_sink( jk, iball )

   END SUBROUTINE detritus

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE detritus( )                    ! Empty routine
      WRITE(*,*) 'detritus: You should not have seen this print! error?'
   END SUBROUTINE detritus
#endif 

   !!======================================================================
END MODULE detritus_mod
