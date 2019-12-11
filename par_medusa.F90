MODULE par_medusa
   !!======================================================================
   !!                        ***  par_medusa  ***
   !! TOP :   set the MEDUSA parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!              -   !  2011-04  (A. Yool) updated for ROAM project
   !!		   -   !  2013-03  (A. Yool) updated for v3.5 NEMO
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   USE par_pisces , ONLY : jp_pisces       !: number of tracers in PISCES
   USE par_pisces , ONLY : jp_pisces_2d    !: number of 2D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_3d    !: number of 3D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_trd   !: number of biological diag in PISCES

   IMPLICIT NONE

   INTEGER, PARAMETER ::   jp_lm      =  jp_pisces      !: 
   INTEGER, PARAMETER ::   jp_lm_2d   =  jp_pisces_2d   !:
   INTEGER, PARAMETER ::   jp_lm_3d   =  jp_pisces_3d   !:
   INTEGER, PARAMETER ::   jp_lm_trd  =  jp_pisces_trd  !:

#if defined key_medusa
   !!---------------------------------------------------------------------
   !!   'key_medusa'                     user defined tracers (MEDUSA)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_medusa     = .TRUE.   !: PTS flag 
# if defined key_roam
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa     =  15      !: number of PTS tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_2d  =  225     !: additional 2d output arrays (used if ln_diatrc=T)
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_3d  =  5       !: additional 3d output arrays (used if ln_diatrc=T)
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_trd =  0       !: number of sms trends for MEDUSA
# else
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa     =  11      !: number of PTS tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_2d  =  90      !: additional 2d output arrays (used if ln_diatrc=T)
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_3d  =  4       !: additional 3d output arrays (used if ln_diatrc=T)
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_trd =  0       !: number of sms trends for MEDUSA
# endif

   ! assign an index in trc arrays for each PTS prognostic variables
   INTEGER, PUBLIC, PARAMETER ::   jpchn  =  jp_lm + 1      !: non-diatom chlorophyll concentration
   INTEGER, PUBLIC, PARAMETER ::   jpchd  =  jp_lm + 2      !: diatom     chlorophyll concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphn  =  jp_lm + 3      !: non-diatom concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphd  =  jp_lm + 4      !: diatom     concentration
   INTEGER, PUBLIC, PARAMETER ::   jpzmi  =  jp_lm + 5      !: microzooplankton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpzme  =  jp_lm + 6      !: mesozooplankton  concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdin  =  jp_lm + 7      !: dissolved inorganic nitrogen concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil  =  jp_lm + 8      !: silicic acid concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer  =  jp_lm + 9      !: total iron concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdet  =  jp_lm + 10     !: slow-sinking detritus concentration
   INTEGER, PUBLIC, PARAMETER ::   jppds  =  jp_lm + 11     !: diatom silicon concentration
# if defined key_roam
   INTEGER, PUBLIC, PARAMETER ::   jpdtc  =  jp_lm + 12     !: slow-sinking detritus carbon concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdic  =  jp_lm + 13     !: dissolved inorganic carbon concentration
   INTEGER, PUBLIC, PARAMETER ::   jpalk  =  jp_lm + 14     !: alkalinity
   INTEGER, PUBLIC, PARAMETER ::   jpoxy  =  jp_lm + 15     !: dissolved oxygen concentration
# endif

   ! assign an index in trc arrays for each PTS prognostic variables
   INTEGER, PUBLIC, PARAMETER ::   jpchn_lc  =  1      !: non-diatom chlorophyll concentration
   INTEGER, PUBLIC, PARAMETER ::   jpchd_lc  =  2      !: diatom     chlorophyll concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphn_lc  =  3      !: non-diatom concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphd_lc  =  4      !: diatom     concentration
   INTEGER, PUBLIC, PARAMETER ::   jpzmi_lc  =  5      !: microzooplankton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpzme_lc  =  6      !: mesozooplankton  concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdin_lc  =  7      !: dissolved inorganic nitrogen concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil_lc  =  8      !: silicic acid concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer_lc  =  9      !: total iron concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdet_lc  =  10     !: slow-sinking detritus concentration
   INTEGER, PUBLIC, PARAMETER ::   jppds_lc  =  11     !: diatom silicon concentration
# if defined key_roam
   INTEGER, PUBLIC, PARAMETER ::   jpdtc_lc  =  12     !: slow-sinking detritus carbon concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdic_lc  =  13     !: dissolved inorganic carbon concentration
   INTEGER, PUBLIC, PARAMETER ::   jpalk_lc  =  14     !: alkalinity
   INTEGER, PUBLIC, PARAMETER ::   jpoxy_lc  =  15     !: dissolved oxygen concentration
# endif

#else
   !!---------------------------------------------------------------------
   !!   Default                           No user defined tracers (MEDUSA)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_medusa     = .FALSE.  !: MEDUSA flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa     =  0       !: No MEDUSA tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_2d  =  0       !: No MEDUSA additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_3d  =  0       !: No MEDUSA additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_medusa_trd =  0       !: number of sms trends for MEDUSA
#endif

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_msa0     = jp_lm     + 1              !: First index of MEDUSA passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_msa1     = jp_lm     + jp_medusa      !: Last  index of MEDUSA passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_msa0_2d  = jp_lm_2d  + 1              !: First index of MEDUSA passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_msa1_2d  = jp_lm_2d  + jp_medusa_2d   !: Last  index of MEDUSA passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_msa0_3d  = jp_lm_3d  + 1              !: First index of MEDUSA passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_msa1_3d  = jp_lm_3d  + jp_medusa_3d   !: Last  index of MEDUSA passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_msa0_trd = jp_lm_trd + 1              !: First index of MEDUSA passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_msa1_trd = jp_lm_trd + jp_medusa_trd  !: Last  index of MEDUSA passive tracers

   !!======================================================================
END MODULE par_medusa
