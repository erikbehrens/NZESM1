MODULE trcwri_cfc
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    cfc :   Output of cfc tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_cfc && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_cfc'                                           cfc model
   !!----------------------------------------------------------------------
   !! trc_wri_cfc   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_cfc 

#  include "top_substitute.h90"
CONTAINS

   SUBROUTINE trc_wri_cfc
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
      DO jn = jp_cfc0, jp_cfc1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, trn(:,:,:,jn) )
      END DO
      !
   END SUBROUTINE trc_wri_cfc

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_cfc
CONTAINS
   SUBROUTINE trc_wri_cfc                     ! Empty routine  
   END SUBROUTINE trc_wri_cfc
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_cfc
