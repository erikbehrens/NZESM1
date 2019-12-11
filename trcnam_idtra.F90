MODULE trcnam_idtra
   !!======================================================================
   !!                         ***  MODULE trcnam_idtra  ***
   !! TOP :   initialisation of some run parameters for IDEAL-TRACER chemical model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.idtra.h90
   !!----------------------------------------------------------------------
#if defined key_idtra
   !!----------------------------------------------------------------------
   !!   'key_idtra'                                               IDEAL-TRACER tracers
   !!----------------------------------------------------------------------
   !! trc_nam_idtra      : IDEAL-TRACER model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_idtra    ! IDEAL-TRACER specific variable
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_idtra   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_idtra
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_idtra  ***
      !!
      !! ** Purpose :   Definition some run parameter for IDEAL-TRACER model
      !!
      !! ** Method  :   Read the namidtra namelist and check the parameter
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namidtra
      !!----------------------------------------------------------------------
      INTEGER  :: numnatm_ref = -1   ! Logical unit for reference ID-TRA namelist
      INTEGER  :: numnatm_cfg = -1   ! Logical unit for configuration ID-TRA namelist
      INTEGER  :: numonc      = -1   ! Logical unit for output namelist
      INTEGER  :: ios                 ! Local integer output status for namelist read
      REAL(wp) :: tmp_decay          !! Years ; half time decay of our idealize tracer
      REAL(wp) :: TDECyr, TDEC   
      !! ----------------------------------------------------------------
      NAMELIST/namidtra/tmp_decay
      !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !! Jpalm -- 4-11-2014
      !! namelist for idealize tracer
      !! only thing in namelist is the chosen half time decay
      !! no atmospheric conditions, cause we do impose a surface concentration of 1,
      !! and no additionnal diagnostics, 
      !! because the only thing we are interested in is the water mass concentration on this tracer.
      !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_idtra: read IDEAL-TRACER namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      !!
      !! Open the namelist file :
      !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL ctl_opn( numnatm_ref, 'namelist_idtra_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatm_cfg, 'namelist_idtra_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonc, 'output.namelist.idtra', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !! Read the namelists :
      !!~~~~~~~~~~~~~~~~~~~~~~~
      !! First namelist of our idealize tracer :
      !! read the decay 1/2 time of our tracer, to define in the namelist.
      !! tmp_decay = 1y ; 10y ; 100y or 1000y depending of which water mass you want to track
      !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REWIND( numnatm_ref )              ! Namelist namidtra in reference namelist : IDTRA parameters
      READ  ( numnatm_ref, namidtra, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namidtra in reference namelist', lwp )

      REWIND( numnatm_cfg )              ! Namelist namidtra in configuration namelist : IDTRA parameters
      READ  ( numnatm_cfg, namidtra, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namidtra in configuration namelist', lwp )
      IF(lwm) WRITE ( numonc, namidtra )

      IF(lwp) WRITE(numout,*) '   -  half time decay of our idealize tracer : ', tmp_decay

      ! decroissance radioactive du traceur ideal
      ! ---------------------------------------
      ! TDECyr = 12.43/LOG(2.)             !! Tricium as example
       TDECyr = tmp_decay/LOG(2.)          !! Idealise tracer -- with tmp_decay given in the idtracer namelist
       TDEC = TDECyr*365.*24.*60.*60.      !! translate in second
       FDEC = EXP( -rdt/TDEC )


!! #if defined key_trc_diaadd  && ! defined key_iomput
      !!
      !!  -Here you can add tracers names to be read
      !! in a namelist.
      !!  -But this is not necessary with the iomput module
      !! cause names are written in the Iodef file.
      !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!! #endif

   END SUBROUTINE trc_nam_idtra

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No IDEAL-TRACER
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_idtra                      ! Empty routine
   END  SUBROUTINE  trc_nam_idtra
#endif

   !!======================================================================
END MODULE trcnam_idtra






