MODULE trcnam_age
   !!======================================================================
   !!                         ***  MODULE trcnam_age  ***
   !! TOP :   initialisation of some run parameters for Age tracer
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) 
   !!----------------------------------------------------------------------
#if defined key_age
   !!----------------------------------------------------------------------
   !!   'key_age'                                               AGE tracers
   !!----------------------------------------------------------------------
   !! trc_nam_age      : AGE  tracer initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trcsms_age      ! AGE specific variable

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_age   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_age
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_age  ***
      !!                 
      !! ** Purpose :   Definition some run parameter for AGE model
      !!
      !! ** input   :   Namelist namage
      !!----------------------------------------------------------------------
      INTEGER ::  numnatg_ref = -1   ! Logical unit for reference AGE namelist
      INTEGER ::  numnatg_cfg = -1   ! Logical unit for configuration AGE namelist
      INTEGER ::  numong      = -1   ! Logical unit for output namelist
      INTEGER :: ios                 ! Local integer output status for namelist read
      INTEGER :: jl, jn
      !!
      NAMELIST/namage/ rn_age_depth, rn_age_kill_rate 
      !!----------------------------------------------------------------------
      !                             ! Open namelist files
      CALL ctl_opn( numnatg_ref, 'namelist_age_ref'   ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatg_cfg, 'namelist_age_cfg'   ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numong, 'output.namelist.age', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      REWIND( numnatg_ref )              ! Namelist namagedate in reference namelist : AGE parameters
      READ  ( numnatg_ref, namage, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namage in reference namelist', lwp )

      REWIND( numnatg_cfg )              ! Namelist namagedate in configuration namelist : AGE parameters
      READ  ( numnatg_cfg, namage, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namage in configuration namelist', lwp )
      IF(lwm) WRITE ( numong, namage )

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_nam_age: Read namage, namelist for Age passive tracer'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '  depth over which age tracer reset to zero                              rn_age_depth      = ', rn_age_depth 
         WRITE(numout,*) '  recip of relax. timescale (s) for age tracer shallower than age_depth  rn_age_kill_rate  = ', rn_age_kill_rate 
      ENDIF

      IF(lwm) CALL FLUSH ( numong )     ! flush output namelist

   END SUBROUTINE trc_nam_age
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No AGE
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_age                      ! Empty routine
   END  SUBROUTINE  trc_nam_age
#endif  

   !!======================================================================
END MODULE trcnam_age
