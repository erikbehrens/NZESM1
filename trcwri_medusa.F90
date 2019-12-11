MODULE trcwri_medusa
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    MEDUSA :   Output of MEDUSA tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!             1.1  !  2013-05 (A. Yool)  converted for MEDUSA
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput && defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'                                           MEDUSA model
   !!----------------------------------------------------------------------
   !! trc_wri_medusa   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE sms_medusa  ! MEDUSA variables
   USE iom         ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_medusa 

CONTAINS

   SUBROUTINE trc_wri_medusa
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      DO jn = jp_msa0, jp_msa1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, trn(:,:,:,jn) )
      END DO
      !
   END SUBROUTINE trc_wri_medusa

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_medusa
CONTAINS
   SUBROUTINE trc_wri_medusa                     ! Empty routine  
   END SUBROUTINE trc_wri_medusa
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_medusa
