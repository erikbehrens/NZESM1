MODULE trcsms_medusa
   !!======================================================================
   !!                         ***  MODULE trcsms_medusa  ***
   !! TOP :   Main module of the MEDUSA tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!              -   !  2017-08  (A. Yool) amend for slow detritus bug
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'                                       bio tracers
   !!----------------------------------------------------------------------
   !! trc_sms_medusa   : MEDUSA_TRC model main routine 
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE trcbio_medusa
   USE trcopt_medusa
   USE trcsed_medusa
   USE trcavg_medusa
   !! for SMS trends
   USE par_medusa,    ONLY: jp_msa0, jp_msa1, jp_medusa
   USE par_oce,       ONLY: jpi, jpj, jpk
   USE trd_oce,       ONLY: jptra_sms, l_trdtrc
   USE trdtrc


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_medusa   ! called by trcsms.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_sms_medusa( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_medusa  ***  
      !!
      !! ** Purpose :   main routine of MEDUSA_TRC model
      !!
      !! ** Method  : - 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt   ! ocean time-step index
      !! Loop variables
      INTEGER :: jn
      !! trend temporary array:
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztrmed


# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) ' MEDUSA inside trc_sms_medusa'
         CALL flush(numout)
# endif

      IF( kt == nittrc000 ) THEN
       IF(lwp) WRITE(numout,*)
       IF(lwp) WRITE(numout,*) ' trc_sms_medusa:  MEDUSA model'
       IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      ENDIF

      !! MEDUSA SMS trends:
      IF( l_trdtrc ) THEN
          CALL wrk_alloc( jpi, jpj, jpk, jp_medusa, ztrmed )
          ztrmed(:,:,:,:)=0.0  
          DO jn = 1, jp_medusa
            ztrmed(:,:,:,jn) = tra(:,:,:,jp_msa0 + jn - 1)
          END DO
      END IF

      CALL trc_avg_medusa( kt ) ! rolling average module
# if defined key_debug_medusa
      IF(lwp) WRITE(numout,*) ' MEDUSA done trc_avg_medusa'
      CALL flush(numout)
# endif
      
      CALL trc_opt_medusa( kt ) ! optical model
# if defined key_debug_medusa
      IF(lwp) WRITE(numout,*) ' MEDUSA done trc_opt_medusa'
      CALL flush(numout)
# endif

      !! AXY & JPALM (28/02/17): call dust before trc_bio_medusa (because of coupling)
      CALL trc_sed_medusa_dust( kt ) ! dust submodel
# if defined key_debug_medusa
      IF(lwp) WRITE(numout,*) ' MEDUSA done trc_sed_medusa_dust'
      CALL flush(numout)
# endif

# if defined key_kill_medusa
      !! MEDUSA skipped
      IF(lwp) WRITE(numout,*) ' MEDUSA killed at kt =', kt
      CALL flush(numout)
# else
      CALL trc_bio_medusa( kt ) ! biological model
#  if defined key_debug_medusa
      IF(lwp) WRITE(numout,*) ' MEDUSA done trc_bio_medusa'
      CALL flush(numout)
#  endif
      
!! AXY (08/08/2017): remove call to buggy subroutine (now handled by detritus.F90)
!!       CALL trc_sed_medusa( kt ) ! sedimentation model
!! #  if defined key_debug_medusa
!!       IF(lwp) WRITE(numout,*) ' MEDUSA done trc_sed_medusa'
!!       CALL flush(numout)
!! #  endif
# endif

      !! MEDUSA SMS trends:
      IF( l_trdtrc ) THEN
          DO jn = 1, jp_medusa
            ztrmed(:,:,:,jn) = tra(:,:,:,jp_msa0 + jn - 1)-ztrmed(:,:,:,jn)
            CALL trd_trc( ztrmed(:,:,:,jn), jn, jptra_sms, kt )   ! save trends
          END DO
          CALL wrk_dealloc( jpi, jpj, jpk, jp_medusa, ztrmed )
      END IF


   END SUBROUTINE trc_sms_medusa
   
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MEDUSA model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_medusa( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_medusa: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_medusa
#endif

   !!======================================================================
END MODULE trcsms_medusa

