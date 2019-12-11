MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the PISCES parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

#if defined key_pisces_reduced
   !!---------------------------------------------------------------------
   !!   'key_pisces_reduced'   :                                LOBSTER bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .TRUE.  !: PISCES flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_p4z        = .FALSE. !: p4z flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     =  6      !: number of passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  =  19     !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  =   3     !: additional 3d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =   17    !: number of sms trends for PISCES

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   INTEGER, PUBLIC, PARAMETER ::   jpdet     =  1        !: detritus                    [mmoleN/m3]
   INTEGER, PUBLIC, PARAMETER ::   jpzoo     =  2        !: zooplancton concentration   [mmoleN/m3]
   INTEGER, PUBLIC, PARAMETER ::   jpphy     =  3        !: phytoplancton concentration [mmoleN/m3]
   INTEGER, PUBLIC, PARAMETER ::   jpno3     =  4        !: nitrate concentration       [mmoleN/m3]
   INTEGER, PUBLIC, PARAMETER ::   jpnh4     =  5        !: ammonium concentration      [mmoleN/m3]
   INTEGER, PUBLIC, PARAMETER ::   jpdom     =  6        !: dissolved organic matter    [mmoleN/m3]

   ! productive layer depth
   INTEGER, PUBLIC, PARAMETER ::   jpkb      = 12        !: first vertical layers where biology is active
   INTEGER, PUBLIC, PARAMETER ::   jpkbm1    = jpkb - 1  !: first vertical layers where biology is active

#elif defined key_pisces  &&  defined key_kriest
   !!---------------------------------------------------------------------
   !!   'key_pisces' & 'key_kriest'                 PISCES bio-model + ???
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .TRUE.  !: PISCES flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_p4z        = .TRUE. !: p4z flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_kriest     = .TRUE.  !: Kriest flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     =  23     !: number of passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  =  13     !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  =  18     !: additional 3d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =   1     !: number of sms trends for PISCES

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   !    WARNING: be carefull about the order when reading the restart
        !   !!gm  this warning should be obsolet with IOM
   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpcal =  4    !: calcite  concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppo4 =  5    !: phosphate concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  6    !: small particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil =  7    !: silicate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  8    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  9    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdoc = 10    !: dissolved organic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpdia = 11    !: Diatoms Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpmes = 12    !: Mesozooplankton Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdsi = 13    !: Diatoms Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer = 14    !: Iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnum = 15    !: Big iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsfe = 16    !: number of particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdfe = 17    !: Diatoms iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgsi = 18    !: (big) Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnfe = 19    !: Nano iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnch = 20    !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdch = 21    !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpno3 = 22    !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnh4 = 23    !: Ammonium Concentration

#elif defined key_pisces
   !!---------------------------------------------------------------------
   !!   'key_pisces'   :                         standard PISCES bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .TRUE.  !: PISCES flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_p4z        = .TRUE.  !: p4z flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_kriest     = .FALSE. !: Kriest flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 24      !: number of PISCES passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 13      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 11      !: additional 3d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =  1      !: number of sms trends for PISCES

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   !    WARNING: be carefull about the order when reading the restart
        !   !!gm  this warning should be obsolet with IOM
   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpcal =  4    !: calcite  concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppo4 =  5    !: phosphate concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  6    !: small particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil =  7    !: silicate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  8    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  9    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdoc = 10    !: dissolved organic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpdia = 11    !: Diatoms Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpmes = 12    !: Mesozooplankton Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdsi = 13    !: Diatoms Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer = 14    !: Iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpbfe = 15    !: Big iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgoc = 16    !: big particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsfe = 17    !: Small iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdfe = 18    !: Diatoms iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgsi = 19    !: (big) Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnfe = 20    !: Nano iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnch = 21    !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdch = 22    !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpno3 = 23    !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnh4 = 24    !: Ammonium Concentration

#else
   !!---------------------------------------------------------------------
   !!   Default                                   No CFC geochemical model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .FALSE.  !: PISCES flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_p4z        = .FALSE.  !: p4z flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     =  0       !: No CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  =  0       !: No CFC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  =  0       !: No CFC additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =  0       !: number of sms trends for PISCES
#endif

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0     = 1                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1     = jp_pisces          !: Last  index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_2d  = 1               !: First index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_2d  = jp_pisces_2d    !: Last  index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_3d  = 1               !: First index of 3D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_3d  = jp_pisces_3d    !: Last  index of 3d diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_trd = 1              !: First index of bio diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_trd = jp_pisces_trd  !: Last  index of bio diag


   !!======================================================================
END MODULE par_pisces
