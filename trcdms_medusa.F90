MODULE trcdms_medusa
   !!======================================================================
   !!                         ***  MODULE trcdms_medusa  ***
   !! TOP :   MEDUSA
   !!======================================================================
   !! History :
   !!  -   !  2014-08  (J. Palmieri - A. Yool)    added for UKESM1 project
   !!  -   !  2017-05  (A. Yool)                  add extra Anderson scheme
   !!----------------------------------------------------------------------
#if defined key_medusa && defined key_roam
   !!----------------------------------------------------------------------
   !!                                        MEDUSA DMS surface concentration
   !!----------------------------------------------------------------------
   !!   trc_dms_medusa        :  
   !!----------------------------------------------------------------------
      USE oce_trc
      USE trc
      USE sms_medusa
      USE lbclnk
      USE prtctl_trc      ! Print control for debugging
      USE in_out_manager  ! I/O manager

      IMPLICIT NONE
      PRIVATE

      PUBLIC   trc_dms_medusa    ! called in trc_bio_medusa

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------


CONTAINS

!=======================================================================
!
   SUBROUTINE trc_dms_medusa( chn, chd, mld, xqsr, xdin, xlim,  &  !! inputs
     &  dms_andr, dms_simo, dms_aran, dms_hall, dms_andm)           !! outputs
!      
!=======================================================================
      !!
      !! Title  : Calculates DMS ocean surface concentration
      !! Author : Julien Palmieri and Andrew Yool
      !! Date   : 08/08/14 
      !!
      !! DMS module is called in trc_bio's huge jk,jj,ji loop
      !! --> DMS concentration is calculated in a specific cell 
      !! (no need of ji,jj,jk)
      !!
      !! AXY (13/03/15): amend to include all four schemes tested
      !!                 during winter/spring 2015; these are:
      !!
      !!                 1. Anderson et al. (2001); this uses fields
      !!                    of surface chl, irradiance and nutrients
      !!                    to empirically estimate DMS via a broken
      !!                    stick approach
      !!
      !!                 2. Simo & Dachs (2002); this uses fields of
      !!                    surface chl and mixed layer depth
      !!
      !!                 3. Aranami & Tsunogai (2004); this is an
      !!                    embellishment of Simo & Dachs
      !!
      !!                 4. Halloran et al. (2010); this is an
      !!                    alternative embellishment of Sim & Dachs
      !!                    and is included because it is formally
      !!                    published (and different from the above)
      !!
      !! AXY (25/05/17): add extra "corrected" Anderson scheme
      !!
      !!                 5. As Anderson et al. (2001) but modified to
      !!                    more accurately reflect nutrient limitation
      !!                    status of phytoplankton community
      !!
      !! AXY (08/07/15): amend to remove Julien's original calculation
      !!                 as this is now superfluous; the four schemes 
      !!                 are calculated and one is chosen to be passed
      !!                 to the atmosphere in trc_bio_medusa
      !!
!=======================================================================

      IMPLICIT NONE
!
      REAL(wp), INTENT( in )    :: chn                  !! non-diatom chlorophyll    (mg/m3)
      REAL(wp), INTENT( in )    :: chd                  !! diatom chlorophyll        (mg/m3)
      REAL(wp), INTENT( in )    :: mld                  !! mix layer depth           (m)
      REAL(wp), INTENT( in )    :: xqsr                 !! surface irradiance        (W/m2)
      REAL(wp), INTENT( in )    :: xdin                 !! surface DIN               (mmol N/m3)
      REAL(wp), INTENT( in )    :: xlim                 !! surface DIN limitation    (mmol N/m3)
      REAL(wp), INTENT( inout ) :: dms_andr             !! DMS surface concentration (nmol/L) 
      REAL(wp), INTENT( inout ) :: dms_simo             !! DMS surface concentration (nmol/L) 
      REAL(wp), INTENT( inout ) :: dms_aran             !! DMS surface concentration (nmol/L) 
      REAL(wp), INTENT( inout ) :: dms_hall             !! DMS surface concentration (nmol/L) 
      REAL(wp), INTENT( inout ) :: dms_andm             !! DMS surface concentration (nmol/L) 
!
      REAL(wp) :: CHL, cmr, sw_dms
      REAL(wp) :: Jterm, Qterm
      !! temporary variables
      REAL(wp) ::    fq1,fq2,fq3
! 
!=======================================================================
!
! AXY (13/03/15): per remarks above, the following calculations estimate
!                 DMS using all of the schemes examined for UKESM1
!
      CHL = 0.0
      CHL = chn+chd                                 !! mg/m3 
      cmr = CHL / mld
!
! AXY (13/03/15): Anderson et al. (2001)
!! JPALM --19-12-2017-- Tunable through the namelist
!!                      within dmsmin - dmscut - dmsslp
        Jterm = xqsr + 1.0e-6
        !! this next line makes a hard-coded assumption about the 
        !! half-saturation constant of MEDUSA (which should be
        !! done properly; perhaps even scaled with the proportion
        !! of diatoms and non-diatoms)
        Qterm = xdin / (xdin + 0.5)
        fq1 = log10(CHL * Jterm * Qterm)
        if (fq1 > dmscut) then
           dms_andr = (dmsslp * (fq1 - dmscut)) + dmsmin
        else
           dms_andr = dmsmin
        endif
!
! AXY (13/03/15): Simo & Dachs (2002)
        fq1 = (-1.0 * log(mld)) + 5.7
        fq2 = (55.8 * cmr) + 0.6
        if (cmr < 0.02) then
           dms_simo = fq1
        else
           dms_simo = fq2
        endif
!           
! AXY (13/03/15): Aranami & Tsunogai (2004)
        fq1 = 60.0 / mld
        fq2 = (55.8 * cmr) + 0.6
        if (cmr < 0.02) then
           dms_aran = fq1
        else
           dms_aran = fq2
        endif
!        
! AXY (13/03/15): Halloran et al. (2010)
        fq1 = (-1.0 * log(mld)) + 5.7
        fq2 = (55.8 * cmr) + 0.6
        fq3 = (90.0 / mld)
        if (cmr < 0.02) then
           dms_hall = fq1
        else
           dms_hall = fq2
        endif
        if (mld > 182.5) then
           dms_hall = fq3
        endif
!
! AXY (25/05/17): modified Anderson et al. (2001)
        Jterm = xqsr + 1.0e-6
        !! this version fixes the hard-coded assumption above
        Qterm = xlim
        fq1 = log10(CHL * Jterm * Qterm)
        if (fq1 > 1.72) then
           dms_andm = (8.24 * (fq1 - 1.72)) + 2.29
        else
           dms_andm = 2.29
        endif

  END SUBROUTINE trc_dms_medusa


!=======================================================================
!=======================================================================
!=======================================================================

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================

CONTAINS

!=======================================================================
!
   SUBROUTINE trc_dms_medusa( kt )                                        
!      
!
      INTEGER, INTENT( in ) ::   kt
!

      WRITE(*,*) 'trc_dms_medusa: You should not have seen this print! error?'

   END SUBROUTINE trc_dms_medusa
#endif

   !!======================================================================
END MODULE trcdms_medusa
