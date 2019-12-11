MODULE trcnam_cfc
   !!======================================================================
   !!                         ***  MODULE trcnam_cfc  ***
   !! TOP :   initialisation of some run parameters for CFC chemical model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.cfc.h90
   !!----------------------------------------------------------------------
#if defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_nam_cfc      : CFC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_cfc      ! CFC specific variable
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_cfc   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_cfc
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_cfc  ***
      !!                 
      !! ** Purpose :   Definition some run parameter for CFC model
      !!
      !! ** Method  :   Read the namcfc namelist and check the parameter 
      !!       values called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namcfc
      !!----------------------------------------------------------------------
      INTEGER ::  numnatc_ref = -1   ! Logical unit for reference CFC namelist
      INTEGER ::  numnatc_cfg = -1   ! Logical unit for configuration CFC namelist
      INTEGER ::  numonc      = -1   ! Logical unit for output namelist
      INTEGER :: ios                 ! Local integer output status for namelist read
      INTEGER :: jl, jn
      !!
      NAMELIST/namcfcdate/ ndate_beg, nyear_res, simu_type 
      !!----------------------------------------------------------------------
      !                             ! Open namelist files
      CALL ctl_opn( numnatc_ref, 'namelist_cfc_ref'   ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatc_cfg, 'namelist_cfc_cfg'   ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonc, 'output.namelist.cfc', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      REWIND( numnatc_ref )              ! Namelist namcfcdate in reference namelist : CFC parameters
      READ  ( numnatc_ref, namcfcdate, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcfcdate in reference namelist', lwp )

      REWIND( numnatc_cfg )              ! Namelist namcfcdate in configuration namelist : CFC parameters
      READ  ( numnatc_cfg, namcfcdate, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcfcdate in configuration namelist', lwp )
      IF(lwm) WRITE ( numonc, namcfcdate )

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_nam: Read namdates, namelist for CFC chemical model'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '    initial calendar date (aammjj) for CFC  ndate_beg = ', ndate_beg
         WRITE(numout,*) '    restoring time constant (year)          nyear_res = ', nyear_res
         IF (simu_type==1) THEN
            WRITE(numout,*) ' CFC running on SPIN-UP mode             simu_type = ', simu_type
         ELSEIF (simu_type==2) THEN
            WRITE(numout,*) ' CFC running on HINDCAST/PROJECTION mode simu_type = ', simu_type
         ENDIF
      ENDIF
      nyear_beg = ndate_beg / 10000
      IF(lwp) WRITE(numout,*) '    initial year (aa)                       nyear_beg = ', nyear_beg
      !

   IF(lwm) CALL FLUSH ( numonc )     ! flush output namelist CFC

   END SUBROUTINE trc_nam_cfc
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No CFC
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_cfc                      ! Empty routine
   END  SUBROUTINE  trc_nam_cfc
#endif  

   !!======================================================================
END MODULE trcnam_cfc
