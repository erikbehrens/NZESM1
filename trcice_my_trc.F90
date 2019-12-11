MODULE trcice_my_trc
   !!======================================================================
   !!                         ***  MODULE trcice_my_trc  ***
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ice_my_trc       : MY_TRC model main routine
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ice_ini_my_trc       ! called by trcice.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ice_ini_my_trc
      !!----------------------------------------------------------------------
      !!                     ***  trc_ice_my_trc  ***
      !!
      !!----------------------------------------------------------------------
      !
      !
   END SUBROUTINE trc_ice_ini_my_trc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ice_ini_my_trc             ! Empty routine
   END SUBROUTINE trc_ice_ini_my_trc
#endif

   !!======================================================================
END MODULE trcice_my_trc
