MODULE trcwri_idtra
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    IDEALIZED Tracer :   Output of IDEALIZED Tracer tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!             1.1  !  2013-05 (A. Yool)  converted for MEDUSA
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput && defined key_idtra
   !!----------------------------------------------------------------------
   !!   'key_idtra'                                           IDEALIZED Tracer model
   !!----------------------------------------------------------------------
   !! trc_wri_idtra   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   ! USE oce_trc         ! Ocean variables
   ! USE par_trc         ! TOP parameters
   USE trc             ! passive tracers common variables
   ! USE trcsms_idtra    ! IDEALIZE TRACER sms trends
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_idtra

CONTAINS

   SUBROUTINE trc_wri_idtra
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
      DO jn = jp_idtra0, jp_idtra1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, trn(:,:,:,jn) )
      END DO
      !
   END SUBROUTINE trc_wri_idtra

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_idtra
CONTAINS
   SUBROUTINE trc_wri_idtra                     ! Empty routine  
   END SUBROUTINE trc_wri_idtra
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_idtra

