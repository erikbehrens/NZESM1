MODULE bdy_par
   !!======================================================================
   !!                      ***  MODULE bdy_par   ***
   !! Unstructured Open Boundary Cond. :   define related parameters
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D. Storkey and E. O'Dea) update for Shelf configurations
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
#if defined   key_bdy
   !!----------------------------------------------------------------------
   !!   'key_bdy' :                    Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PUBLIC

# if ! defined key_agrif
   LOGICAL, PUBLIC, PARAMETER ::   lk_bdy  = .TRUE.   !: Unstructured Ocean Boundary Condition flag
# else
   LOGICAL, PUBLIC            ::   lk_bdy  = .TRUE.   !: Unstructured Ocean Boundary Condition flag
# endif
   INTEGER, PUBLIC, PARAMETER ::   jp_bdy  = 10       !: Maximum number of bdy sets
   INTEGER, PUBLIC, PARAMETER ::   jpbgrd  = 3	      !: Number of horizontal grid types used  (T, U, V)

#else
   !!----------------------------------------------------------------------
   !!   Default option :            NO Unstructured open boundary condition
   !!----------------------------------------------------------------------
# if ! defined key_agrif
   LOGICAL, PUBLIC, PARAMETER ::   lk_bdy  = .FALSE.   !: Unstructured Ocean Boundary Condition flag
# else
   LOGICAL, PUBLIC            ::   lk_bdy  = .FALSE.   !: Unstructured Ocean Boundary Condition flag
# endif
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE bdy_par
