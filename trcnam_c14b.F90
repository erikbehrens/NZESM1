MODULE trcnam_c14b
   !!======================================================================
   !!                         ***  MODULE trcnam_c14b  ***
   !! TOP :   initialisation of some run parameters for C14 chemical model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.cfc.h90
   !!----------------------------------------------------------------------
#if defined key_c14b
   !!----------------------------------------------------------------------
   !!   'key_c14b'                                         C14 bomb tracer
   !!----------------------------------------------------------------------
   !! trc_nam_c14b      : C14 model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_c14b     ! C14b specific variable
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_c14b   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_c14b
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_c14b  ***
      !!                 
      !! ** Purpose :   Definition some run parameter for C14 model
      !!
      !! ** Method  :   Read the namc14 namelist and check the parameter 
      !!       values called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namelist_c14b
      !!----------------------------------------------------------------------
      INTEGER ::  numnatb_ref = -1   ! Logical unit for reference c14b namelist
      INTEGER ::  numnatb_cfg = -1   ! Logical unit for configuration c14b namelist
      INTEGER ::  numonb      = -1   ! Logical unit for output namelist
      INTEGER :: ios                 ! Local integer output status for namelist read

      ! definition of additional diagnostic as a structure
      INTEGER :: jl, jn
      !!
      NAMELIST/namc14date/ ndate_beg_b, nyear_res_b
      !!-------------------------------------------------------------------
      !                             ! Open namelist file
      CALL ctl_opn( numnatb_ref, 'namelist_c14b_ref'  ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatb_cfg, 'namelist_c14b_cfg'  ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )   
      IF(lwm) CALL ctl_opn( numonb, 'output.namelist.c14', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )     
      REWIND( numnatb_ref )              ! Namelist namc14date in reference namelist : c14b parameters
      READ  ( numnatb_ref, namc14date, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namc14date in reference namelist', lwp )

      REWIND( numnatb_cfg )              ! Namelist namc14date in configuration namelist : c14b parameters
      READ  ( numnatb_cfg, namc14date, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namc14date in configuration namelist', lwp )
      IF(lwm) WRITE ( numonb, namc14date )
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_nam: Read namdates, namelist for C14 chemical model'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '    initial calendar date (aammjj) for C14  ndate_beg_b = ', ndate_beg_b
         WRITE(numout,*) '    restoring time constant (year)          nyear_res_b = ', nyear_res_b
      ENDIF
      nyear_beg_b = ndate_beg_b / 10000
      IF(lwp) WRITE(numout,*) '    initial year (aa)                  nyear_beg_b = ', nyear_beg_b
      !

   IF(lwm) CALL FLUSH ( numonb )     ! flush output namelist C14b

   END SUBROUTINE trc_nam_c14b
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No 14C
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_c14b                      ! Empty routine
   END  SUBROUTINE  trc_nam_c14b
#endif  

   !!======================================================================
END MODULE trcnam_c14b
