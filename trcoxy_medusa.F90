MODULE trcoxy_medusa
   !!======================================================================
   !!                         ***  MODULE trcoxy_medusa  ***
   !! TOP :   MEDUSA
   !!======================================================================
   !! History :
   !!  -   !  2011-07  (A. Yool)             added for ROAM project
   !!----------------------------------------------------------------------
#if defined key_medusa && defined key_roam
   !!----------------------------------------------------------------------
   !!                                        MEDUSA oxygen cycle
   !!----------------------------------------------------------------------
   !!   trc_oxy_medusa        :  
   !!----------------------------------------------------------------------
      USE oce_trc
      USE trc
      USE sms_medusa
      USE lbclnk
      USE prtctl_trc      ! Print control for debugging

      IMPLICIT NONE
      PRIVATE
      
      PUBLIC   trc_oxy_medusa    ! called in trc_bio_medusa
      PUBLIC   oxy_sato          ! called in trc_bio_medusa

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

! The following is a map of the subroutines contained within this module
! - trc_oxy_medusa
!      - CALLS oxy_schmidt
!      - CALLS oxy_sato

CONTAINS

!=======================================================================
!
   SUBROUTINE trc_oxy_medusa( pt, ps, kw660, pp0, o2,  &  !! inputs
      kwo2, o2flux, o2sat )                               !! outputs
!      
!=======================================================================
      !!
      !! Title  : Calculates O2 change due to air-sea gas exchange
      !! Author : Andrew Yool
      !! Date   : 15/10/04 (revised 08/07/11)
      !!
      !! This subroutine calculates oxygen air-sea gas exchange in the
      !! surface layer of the ocean.  The formulation is taken from one 
      !! written by Ray Najjar for the OCMIP-2 project.  The routine
      !! calls two other subroutines, oxy_schmidt.f (calculates Schmidt 
      !! number of oxygen) and oxy_sato.f (calculates oxygen saturation 
      !! concentration at 1 atm).
      !!
      !! AXY (23/06/15): revised to allow common gas transfer velocity
      !!                 to be used for CO2 and O2; outputs of this
      !!                 routine amended to mmol/m3 from mol/m3
      !! 
      !! Function inputs are (in order) : 
      !!     pt      temperature                     (degrees C)
      !!     ps      salinity                        (PSU)
      !!     kw660   gas transfer velocity           (m/s)
      !!     pp0     surface pressure                (divided by 1 atm)
      !!     o2      surface O2 concentration        (mmol/m3)
      !! (+) kwo2    gas transfer velocity for O2    (m/s)
      !! (*) o2flux  exchange rate of oxygen         (mmol/m2/s)
      !! (+) o2sat   oxygen saturation concentration (mmol/m3)
      !! 
      !! Where (*) is the function output (note its units).  
      !!
!=======================================================================

      implicit none
!
      REAL(wp), INTENT( in )    :: pt
      REAL(wp), INTENT( in )    :: ps
      REAL(wp), INTENT( in )    :: kw660
      REAL(wp), INTENT( in )    :: pp0
      REAL(wp), INTENT( in )    :: o2
      REAL(wp), INTENT( out )   :: kwo2, o2flux, o2sat
!
      REAL(wp) :: o2schmidt, o2sato, mol_o2
!
! Oxygen to mol / m3
!
      mol_o2 = o2 / 1000.
!
! Calculate oxygen Schmidt number
! 
      call oxy_schmidt(pt, o2schmidt)
!
! Calculate the transfer velocity for O2 (m/s)
!
      kwo2 = kw660 * (660 / o2schmidt)**0.5
!
! Calculate the saturation concentration for oxygen (mol/m3)
!
      call oxy_sato(pt, ps, o2sato)
      o2sat = o2sato * pp0
!
! Calculate time rate of change of O2 due to gas exchange (mol/m3/s)
!
      o2flux = kwo2 * (o2sat - mol_o2)
!
! Oxygen flux and saturation to mmol / m3
!
      o2sat  =  o2sat * 1000.
      o2flux = o2flux * 1000.
!
      END SUBROUTINE trc_oxy_medusa

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
   SUBROUTINE oxy_schmidt( pt, &  !! input
      o2_schmidt )                !! output
!      
!=======================================================================
      !!
      !! Title  : Calculates Schmidt number for ocean uptake of O2
      !! Author : Andrew Yool
      !! Date   : 14/10/04 (revised 08/07/11)
      !! 
      !! This subroutine calculates the Schmidt number for O2 using sea 
      !! surface temperature.  The code is based upon that developed as 
      !! part of the OCMIP-2 project (1998-2000).  The coefficients used 
      !! are taken from Keeling et al. (1998, GBC, 12, 141-163).
      !! 
      !! AXY (23/06/2015)
      !! UPDATED: revised formulation from Wanninkhof (2014) for
      !! consistency with MOCSY
      !!
      !! Winninkhof, R. (2014). Relationship between wind speed and gas
      !! exchange over the ocean revisited. LIMNOLOGY AND OCEANOGRAPHY-METHODS
      !! 12, 351-362, doi:10.4319/lom.2014.12.351
      !!
      !! Function inputs are (in order) : 
      !!     t           temperature (degrees C)
      !! (*) o2_schmidt  oxygen Schmidt number
      !! 
      !! Where (*) is the function output.
      !!
!=======================================================================

      implicit none
!
      REAL(wp) :: pt, o2_schmidt
      REAL(wp) :: a0, a1, a2, a3, a4
!
! AXY (23/06/15): OCMIP-2 coefficients
!     data a0 /    1638.0 /
!     data a1 /    -81.83 /
!     data a2 /     1.483 /
!     data a3 / -0.008004 /
!
! AXY (23/06/15): Wanninkhof (2014) coefficients
      data a0 /     1920.4 /
      data a1 /     -135.6 /
      data a2 /     5.2121 /
      data a3 /   -0.10939 /
      data a4 / 0.00093777 /
!
!     o2_schmidt = a0 + pt*(a1 + pt*(a2 + pt*a3))
      o2_schmidt = a0 + pt*(a1 + pt*(a2 + pt*(a3 + pt*a4)))
!
      END SUBROUTINE oxy_schmidt

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
   SUBROUTINE oxy_sato( pt, ps, &  !! inputs 
      o2_sato )                    !! output
!      
!=======================================================================
      !!
      !! Title  : Calculates O2 saturation at 1 atm pressure
      !! Author : Andrew Yool
      !! Date   : 14/10/04 (revised 08/07/11)
      !! 
      !! This subroutine calculates the oxygen saturation concentration
      !! at 1 atmosphere pressure in mol/m3 given ambient temperature
      !! and salinity.  This formulation is (ostensibly) taken from 
      !! Garcia & Gordon (1992, L&O, 1307-1312).  The function works
      !! in the range -1.9 <= T <= 40, 0 <= S <= 42.
      !!
      !! Function inputs are (in order) : 
      !!     pt       temperature (degrees C)
      !!     ps       salinity (PSU)
      !! (*) o2_sato  oxygen saturation (mol/m3)
      !!
      !! Where (*) is the function output (note its units).  
      !!
      !! Check value : T = 10, S = 35, oxy_sato = 0.282015 mol/m3
      !!
!=======================================================================

      implicit none
!
      REAL(wp) :: pt, ps, o2_sato
!
      REAL(wp) :: a0, a1, a2, a3, a4, a5
      REAL(wp) :: b0, b1, b2, b3
      REAL(wp) :: c0
!
      REAL(wp) :: tt, tk, ts, ts2, ts3, ts4, ts5
      REAL(wp) :: ans1, ans2
!
      data a0 /  2.00907    /
      data a1 /  3.22014    /
      data a2 /  4.05010    /
      data a3 /  4.94457    /
      data a4 / -2.56847E-1 /
      data a5 /  3.88767    /
!
      data b0 / -6.24523E-3 /
      data b1 / -7.37614E-3 /
      data b2 / -1.03410E-2 /
      data b3 / -8.17083E-3 /
!
      data c0 / -4.88682E-7 /
!     
      tt   = 298.15 - pt
      tk   = 273.15 + pt
      ts   = log(tt / tk)
      ts2  = ts**2
      ts3  = ts**3
      ts4  = ts**4
      ts5  = ts**5
      ans1 = a0 + a1*ts + a2*ts2 + a3*ts3 + a4*ts4 + a5*ts5  &
           + ps*(b0 + b1*ts + b2*ts2 + b3*ts3)               &
           + c0*(ps*ps)
      ans2 = exp(ans1)
!
!  Convert from ml/l to mol/m3
!
      o2_sato = (ans2 / 22391.6) * 1000.0
!
      END SUBROUTINE oxy_sato

!=======================================================================
!=======================================================================
!=======================================================================

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================

CONTAINS

   SUBROUTINE trc_oxy_medusa( pt, ps, kw660, pp0, o2,  &  !! inputs
      o2flux, o2sat )                                     !! outputs
      USE par_kind

      REAL(wp), INTENT( in )    :: pt
      REAL(wp), INTENT( in )    :: ps
      REAL(wp), INTENT( in )    :: kw660
      REAL(wp), INTENT( in )    :: pp0
      REAL(wp), INTENT( in )    :: o2
      REAL(wp), INTENT( inout ) :: o2flux, o2sat

      WRITE(*,*) 'trc_oxy_medusa: You should not have seen this print! error?', kt

   END SUBROUTINE trc_oxy_medusa
#endif

   !!======================================================================
END MODULE trcoxy_medusa
