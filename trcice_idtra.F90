MODULE trcice_idtra
   !!======================================================================
   !!                         ***  MODULE trcice_idtra  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_idtra
   !!----------------------------------------------------------------------
   !!   'key_idtra'                                    IDEAL TRACER tracers
   !!----------------------------------------------------------------------
   !! trc_ice_idtra       : MY_TRC model main routine
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ice_ini_idtra       ! called by trcice.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ice_ini_idtra
      !!----------------------------------------------------------------------
      !!                     ***  trc_ice_idtra  ***
      !!
      !! ** Purpose :   main routine of MY_TRC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      !
   END SUBROUTINE trc_ice_ini_idtra


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ice_ini_idtra             ! Empty routine
   END SUBROUTINE trc_ice_ini_idtra
#endif

   !!======================================================================
END MODULE trcice_idtra
