MODULE par_idtra
   !!======================================================================
   !!                        ***  par_idtra  ***
   !! TOP :   set the IDEAL-TRACER parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   USE par_pisces , ONLY : jp_pisces       !: number of tracers in PISCES
   USE par_pisces , ONLY : jp_pisces_2d    !: number of 2D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_3d    !: number of 3D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_trd   !: number of biological diag in PISCES

   USE par_medusa , ONLY : jp_medusa       !: number of tracers in MEDUSA
   USE par_medusa , ONLY : jp_medusa_2d    !: number of 2D diag in MEDUSA
   USE par_medusa , ONLY : jp_medusa_3d    !: number of 3D diag in MEDUSA
   USE par_medusa , ONLY : jp_medusa_trd   !: number of biological diag in MEDUSA

   IMPLICIT NONE

   INTEGER, PARAMETER ::   jp_lp      =  jp_pisces     +  jp_medusa     !: cumulative number of passive tracers
   INTEGER, PARAMETER ::   jp_lp_2d   =  jp_pisces_2d  +  jp_medusa_2d  !:
   INTEGER, PARAMETER ::   jp_lp_3d   =  jp_pisces_3d  +  jp_medusa_3d  !:
   INTEGER, PARAMETER ::   jp_lp_trd  =  jp_pisces_trd +  jp_medusa_trd !:

#if defined key_idtra
   !!---------------------------------------------------------------------
   !!   'key_idtra'   :                                          Ideal tracers
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_idtra     = .TRUE.      !: IDEAL-TRACER flag
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra     =  1          !: number of passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra_2d  =  3          !: additional 2d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra_3d  =  0          !: additional 3d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra_trd =  0          !: number of sms trends for IDEAL-TRACER

   ! assign an index in trc arrays for each IDEAL-TRACER prognostic variables
   INTEGER, PUBLIC, PARAMETER ::   jpidtra       = jp_lp + 1   !: IDEAL-TRACER
#else
   !!---------------------------------------------------------------------
   !!   Default     :                                       No IDEAL-TRACER tracers
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_idtra     = .FALSE.     !: IDEAL-TRACER flag
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra     =  0          !: No IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra_2d  =  0          !: No IDEAL-TRACER additional 2d output arrays
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra_3d  =  0          !: No IDEAL-TRACER additional 3d output arrays
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra_trd =  0          !: number of sms trends for IDEAL-TRACER
#endif

   ! Starting/ending IDEAL-TRACER do-loop indices (N.B. no IDEAL-TRACER : jp_idtra0 > jp_idtra1 the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra0     = jp_lp + 1                !: First index of IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra1     = jp_lp + jp_idtra         !: Last  index of IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra0_2d  = jp_lp_2d  + 1            !: First index of IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra1_2d  = jp_lp_2d  + jp_idtra_2d  !: Last  index of IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra0_3d  = jp_lp_3d  + 1            !: First index of IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra1_3d  = jp_lp_3d  + jp_idtra_3d  !: Last  index of IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra0_trd = jp_lp_trd + 0            !: First index of IDEAL-TRACER tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_idtra1_trd = jp_lp_trd + jp_idtra_trd !: Last  index of IDEAL-TRACER tracers

   !!======================================================================
END MODULE par_idtra



