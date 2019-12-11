MODULE trcice_c14b
   !!======================================================================
   !!                         ***  MODULE trcice_c14b  ***
   !!======================================================================
#if defined key_c14b
   !!----------------------------------------------------------------------
   !!   'key_c14b'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ice_c14b       : MY_TRC model main routine
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ice_ini_c14b       ! called by trcice.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ice_ini_c14b
      !!----------------------------------------------------------------------
      !!                     ***  trc_ice_c14b  ***
      !!
      !!----------------------------------------------------------------------
      !
      !
   END SUBROUTINE trc_ice_ini_c14b


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ice_ini_c14b             ! Empty routine
   END SUBROUTINE trc_ice_ini_c14b
#endif

   !!======================================================================
END MODULE trcice_c14b
