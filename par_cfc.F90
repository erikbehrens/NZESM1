MODULE par_cfc
   !!======================================================================
   !!                        ***  par_cfc  ***
   !! TOP :   set the CFC parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!                  !  2017-04  (A. Yool)  add SF6
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   USE par_pisces , ONLY : jp_pisces       !: number of tracers in PISCES
   USE par_pisces , ONLY : jp_pisces_2d    !: number of 2D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_3d    !: number of 3D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_trd   !: number of biological diag in PISCES

   USE par_medusa , ONLY : jp_medusa       !: number of tracers in MEDUSA
   USE par_medusa , ONLY : jp_medusa_2d    !: number of 2D diag in MEDUSA
   USE par_medusa , ONLY : jp_medusa_3d    !: number of 3D diag in MEDUSA
   USE par_medusa , ONLY : jp_medusa_trd   !: number of biological diag in MEDUSA

   USE par_idtra  , ONLY : jp_idtra        !: number of tracers in ideal tracer
   USE par_idtra  , ONLY : jp_idtra_2d     !: number of tracers in ideal tracer
   USE par_idtra  , ONLY : jp_idtra_3d     !: number of tracers in ideal tracer
   USE par_idtra  , ONLY : jp_idtra_trd    !: number of tracers in ideal tracer

   IMPLICIT NONE

   INTEGER, PARAMETER ::   jp_lc      =  jp_pisces     + jp_medusa     + &
                      jp_idtra     !: cumulative number of passive tracers
   INTEGER, PARAMETER ::   jp_lc_2d   =  jp_pisces_2d  + jp_medusa_2d  + &
                      jp_idtra_2d !:
   INTEGER, PARAMETER ::   jp_lc_3d   =  jp_pisces_3d  + jp_medusa_3d  + &
                      jp_idtra_3d !:
   INTEGER, PARAMETER ::   jp_lc_trd  =  jp_pisces_trd + jp_medusa_trd + &
                      jp_idtra_trd !:
   
#if defined key_cfc
   !!---------------------------------------------------------------------
   !!   'key_cfc'   :                                          CFC tracers
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_cfc     = .TRUE.      !: CFC flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc     =  3          !: number of passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc_2d  =  6          !: additional 2d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc_3d  =  0          !: additional 3d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc_trd =  0          !: number of sms trends for CFC
   
   ! assign an index in trc arrays for each CFC prognostic variables
   INTEGER, PUBLIC, PARAMETER ::   jpc11       = jp_lc + 1   !: CFC-11 
   INTEGER, PUBLIC, PARAMETER ::   jpc12       = jp_lc + 2   !: CFC-12 (priority tracer for CMIP6)
   INTEGER, PUBLIC, PARAMETER ::   jpsf6       = jp_lc + 3   !: SF6
#else
   !!---------------------------------------------------------------------
   !!   Default     :                                       No CFC tracers
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_cfc     = .FALSE.     !: CFC flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc     =  0          !: No CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc_2d  =  0          !: No CFC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc_3d  =  0          !: No CFC additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc_trd =  0          !: number of sms trends for CFC
#endif

   ! Starting/ending CFC do-loop indices (N.B. no CFC : jp_cfc0 > jp_cfc1 the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc0     = jp_lc + 1              !: First index of CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc1     = jp_lc + jp_cfc         !: Last  index of CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc0_2d  = jp_lc_2d  + 1          !: First index of CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc1_2d  = jp_lc_2d  + jp_cfc_2d  !: Last  index of CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc0_3d  = jp_lc_3d  + 1          !: First index of CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc1_3d  = jp_lc_3d  + jp_cfc_3d  !: Last  index of CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc0_trd = jp_lc_trd + 1          !: First index of CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cfc1_trd = jp_lc_trd + jp_cfc_trd !: Last  index of CFC tracers

   !!======================================================================
END MODULE par_cfc
