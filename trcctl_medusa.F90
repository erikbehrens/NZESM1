MODULE trcctl_medusa
   !!======================================================================
   !!                         ***  trcctl_medusa.F90  ***
   !! TOP :                Control of MEDUSA_TRC biogeochemical model
   !!======================================================================
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-12 (C. Ethe) assign a parameter to name individual tracers
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!----------------------------------------------------------------------

#if defined key_medusa

   USE oce_trc
   USE trc
   USE in_out_manager
   USE par_medusa
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_ctl_medusa     ! called by ???


   !!----------------------------------------------------------------------
   !! NEMO/TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_ctl_medusa
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ctl_medusa  ***
      !!
      !! ** Purpose :   control the cpp options, namelist and files 
      !!----------------------------------------------------------------------

!      IF(lwp) WRITE(numout,*)
!      IF(lwp) WRITE(numout,*) 'use MEDUSA biological model  '

! Check number of tracers
! -----------------------
# if defined key_roam
      IF (jp_medusa /= 15) THEN 
          IF (lwp) THEN 
              WRITE (numout,*) ' ===>>>> : W A R N I N G '
              WRITE (numout,*) ' =======   ============= '
              WRITE (numout,*)                             &
              &   ' STOP, change jp_medusa to 15 in '      &
              &  ,' par_medusa.F90 '  
          END IF 
          STOP 'TRC_CTL'
      END IF 
# else
      IF (jp_medusa /= 11) THEN 
          IF (lwp) THEN 
              WRITE (numout,*) ' ===>>>> : W A R N I N G '
              WRITE (numout,*) ' =======   ============= '
              WRITE (numout,*)                             &
              &   ' STOP, change jp_medusa to 11 in '      &
              &  ,' par_medusa.F90 '  
          END IF 
          STOP 'TRC_CTL'
      END IF 
# endif

   END SUBROUTINE trc_ctl_medusa

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                            No MEDUSA
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ctl_medusa               ! Dummy routine
   END SUBROUTINE trc_ctl_medusa
#endif

   !!======================================================================
END MODULE trcctl_medusa
      
