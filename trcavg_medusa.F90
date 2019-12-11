MODULE trcavg_medusa
   !!======================================================================
   !!                         ***  MODULE trcavg_medusa  ***
   !! TOP :   MEDUSA
   !!======================================================================
   !! History :    -   !  2015-07  (A. Yool) Original code
   !!----------------------------------------------------------------------
#if defined key_medusa && defined key_roam
   !!----------------------------------------------------------------------
   !!                                        MEDUSA rolling averages
   !!----------------------------------------------------------------------
   !!   trc_avg_medusa        :  
   !!----------------------------------------------------------------------
      USE oce_trc
      USE trc
      USE sms_medusa
      USE lbclnk
      USE prtctl_trc      ! Print control for debugging
      USE in_out_manager  ! I/O manager
      IMPLICIT NONE
      PRIVATE

      PUBLIC   trc_avg_medusa    ! called in trc_sms_medusa

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

!=======================================================================
!
      SUBROUTINE trc_avg_medusa( kt )
!      
!=======================================================================
      !!
      !! Title  : Calculates rolling averages of variables
      !! Author : Andrew Yool
      !! Date   : 23/07/15
      !!
      !! Calculates and updates rolling averages of properties such
      !! as surface irradiance and mixed layer depth that are used
      !! in functions that require average rather than instantaneous
      !! values.
      !!
      !! This functionality was originally added to support the
      !! calculation of surface DMS concentrations - and was done so
      !! within the trcbio_meduse.F90 routine - but was moved to
      !! this separate module so that its calculations could be used
      !! to inform MEDUSA's submarine irradiance field
      !!
!=======================================================================
!
      IMPLICIT NONE
!
      INTEGER, INTENT( in ) ::   kt   ! index of the time stepping
! 
!=======================================================================
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) ' MEDUSA inside trc_avg_medusa'
         CALL flush(numout)
# endif
      !! AXY (24/07/15): alter this to report on first MEDUSA call
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_avg_medusa: MEDUSA rolling average'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
         IF(lwp) WRITE(numout,*) ' kt =',kt
      ENDIF
      !!
      !!----------------------------------------------------------------------
      !! Process average fields
      !! The empirical formulae used for estimating surface DMS concentrations
      !! require temporally averaged input fields; this block calculates these
      !! averages based on diel averages; note that rdt 
      !!----------------------------------------------------------------------
      !!
      zn_dms_chn(:,:) = ( zb_dms_chn(:,:) * ((86400. - rdt) / 86400.) ) &
      &                  + ( trn(:,:,1,jpchn) * (rdt / 86400.) )
      zb_dms_chn(:,:) = zn_dms_chn(:,:)
      zn_dms_chd(:,:) = ( zb_dms_chd(:,:) * ((86400. - rdt) / 86400.) ) &
      &                  + ( trn(:,:,1,jpchd) * (rdt / 86400.) )
      zb_dms_chd(:,:) = zn_dms_chd(:,:)
      zn_dms_mld(:,:) = ( zb_dms_mld(:,:) * ((86400. - rdt) / 86400.) ) &
      &                  + (        hmld(:,:) * (rdt / 86400.) )
      zb_dms_mld(:,:) = zn_dms_mld(:,:)
      zn_dms_qsr(:,:) = ( zb_dms_qsr(:,:) * ((86400. - rdt) / 86400.) ) &
      &                  + (         qsr(:,:) * (rdt / 86400.) )
      zb_dms_qsr(:,:) = zn_dms_qsr(:,:)
      zn_dms_din(:,:) = ( zb_dms_din(:,:) * ((86400. - rdt) / 86400.) ) &
      &                  + ( trn(:,:,1,jpdin) * (rdt / 86400.) )
      zb_dms_din(:,:) = zn_dms_din(:,:)

  END SUBROUTINE trc_avg_medusa

!=======================================================================
!=======================================================================
!=======================================================================

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================

CONTAINS

!=======================================================================
!
   SUBROUTINE trc_avg_medusa( kt )                                      
!      
!
      INTEGER, INTENT( in ) ::   kt
!

      WRITE(*,*) 'trc_avg_medusa: You should not have seen this print! error?'

   END SUBROUTINE trc_avg_medusa
#endif

   !!======================================================================
END MODULE trcavg_medusa
