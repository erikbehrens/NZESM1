MODULE obs_const
   !!=====================================================================
   !!                       ***  MODULE obs_const  ***
   !! Observation diagnostics: Constants used by many modules
   !!===================================================================== 
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   !! * Modules used
   USE par_kind, ONLY : & ! Precision variables
      & sp         
   IMPLICIT NONE

   !! * Routine/type accessibility
   PUBLIC

   REAL(kind=sp), PARAMETER :: obfillflt=99999.

END MODULE obs_const

