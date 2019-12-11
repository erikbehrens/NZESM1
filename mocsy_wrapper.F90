MODULE mocsy_wrapper
   !!======================================================================
   !!                         ***  MODULE mocsy_wrapper  ***
   !! TOP :   MEDUSA
   !!======================================================================
   !! History :
   !!  -   !  2015-06  (A. Yool)             added for UKESM project
   !!  -   !  2017-04  (A. Yool)             alter optK1K2 to 'w14' option
   !!----------------------------------------------------------------------
#if defined key_medusa && defined key_roam
   !!----------------------------------------------------------------------
   !!                                        MEDUSA carbonate chemistry
   !!----------------------------------------------------------------------
   !!   mocsy_wrapper
   !!----------------------------------------------------------------------

   USE mocsy_mainmod
   USE mocsy_gasflux
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC mocsy_interface  ! called in trc_bio_medusa
   PUBLIC mocsy_carbchem   ! called in trc_bio_medusa

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

! The following is a map of the subroutines contained within this module
! - mocsy_interface
!    - CALLS mocsy_carbchem
!
! - mocsy_carbchem
!    - CALLS vars
!      - CALLS p80
!      - CALLS constants
!      - CALLS rho
!      - CALLS varsolver
!        - CALLS [phsolvers routines]
!      - CALLS varsolver (again)
!    - CALLS x2pCO2atm
!    - CALLS p2fCO2
!    - CALLS schmidt_co2
!    - CALLS surface_K0

CONTAINS

!=======================================================================
!
      SUBROUTINE mocsy_interface( temp, sal, alk, dic, sil, phos, &
           Patm, depth, lat, kw660, xco2, N,                      &
           ph, pco2, fco2, co2, hco3, co3, OmegaA,                &
           OmegaC, BetaD, rhoSW, p, tempis,                       &
           pco2atm, fco2atm, schmidtco2, kwco2, K0,               &
           co2starair, co2flux, dpco2 )
!
!=======================================================================
!
! AXY (26/06/15): to preserve both MEDUSA's scalar variables and
!                 MOCSY's vector variables on code traceability
!                 grounds, this additional wrapper does a rather
!                 superfluous conversion between the two; in the
!                 fullness of time, this should be dispensed with
!                 and MOCSY used in its full vector form
!
        USE mocsy_singledouble
        IMPLICIT NONE

!> ======================================================================
!  VARIABLES
!> ======================================================================
!
   INTEGER, INTENT(in) :: N
!
!  MEDUSA-side
!  Input variables
   REAL(kind=wp), INTENT(in)  :: temp, sal, alk, dic, sil, phos
   REAL(kind=wp), INTENT(in)  :: Patm, depth, lat, kw660, xco2
!
!  Output variables
   REAL(kind=wp), INTENT(out) :: ph, pco2, fco2, co2, hco3, co3, OmegaA
   REAL(kind=wp), INTENT(out) :: OmegaC, BetaD, rhoSW, p, tempis
   REAL(kind=wp), INTENT(out) :: pco2atm, fco2atm, schmidtco2, kwco2, K0
   REAL(kind=wp), INTENT(out) :: co2starair, co2flux, dpco2
!
!  MOCSY-side
!  Input variables
   REAL(kind=wp), DIMENSION(N) :: mtemp, msal, malk, mdic, msil, mphos
   REAL(kind=wp), DIMENSION(N) :: mPatm, mdepth, mlat, mkw660, mxco2
!
!  Output variables
   REAL(kind=wp), DIMENSION(N) :: mph, mpco2, mfco2, mco2, mhco3, mco3, mOmegaA
   REAL(kind=wp), DIMENSION(N) :: mOmegaC, mBetaD, mrhoSW, mp, mtempis
   REAL(kind=wp), DIMENSION(N) :: mpco2atm, mfco2atm, mschmidtco2, mkwco2, mK0
   REAL(kind=wp), DIMENSION(N) :: mco2starair, mco2flux, mdpco2
!
!> ----------------------------------------------------------------------
!  Set MOCSY inputs to equal MEDUSA inputs (amend units here)
!> ----------------------------------------------------------------------
!
      mtemp(1)  = temp         ! degrees C
      msal(1)   = sal          ! PSU
      malk(1)   = alk / 1000.  ! meq  / m3 -> eq  / m3
      mdic(1)   = dic / 1000.  ! mmol / m3 -> mol / m3
      msil(1)   = sil / 1000.  ! mmol / m3 -> mol / m3
      mphos(1)  = phos / 1000. ! mmol / m3 -> mol / m3
      mPatm(1)  = Patm         ! atm
      mdepth(1) = depth        ! m
      mlat(1)   = lat          ! degrees N
      mkw660(1) = kw660        ! m / s
      mxco2(1)  = xco2         ! ppm
!
!> ----------------------------------------------------------------------
!  Call MOCSY
!> ----------------------------------------------------------------------
!
      CALL mocsy_carbchem( mtemp, msal, malk, mdic, msil, mphos, & ! inputs
      mPatm, mdepth, mlat, mkw660, mxco2, 1,                     & ! inputs
      mph, mpco2, mfco2, mco2, mhco3, mco3, mOmegaA,             & ! outputs
      mOmegaC, mBetaD, mrhoSW, mp, mtempis,                      & ! outputs
      mpco2atm, mfco2atm, mschmidtco2, mkwco2, mK0,              & ! outputs
      mco2starair, mco2flux, mdpco2 ) ! outputs
!
!> ----------------------------------------------------------------------
!  Set MOCSY outputs to equal MEDUSA outputs (amend units here)
!> ----------------------------------------------------------------------
!
      ph         = mph(1)                 ! standard units
      pco2       = mpco2(1)               ! uatm
      fco2       = mfco2(1)               ! uatm
      co2        = mco2(1) * 1000.        ! mol / m3 -> mmol / m3
      hco3       = mhco3(1) * 1000.       ! mol / m3 -> mmol / m3
      co3        = mco3(1) * 1000.        ! mol / m3 -> mmol / m3
      OmegaA     = mOmegaA(1)             ! dimensionless
      OmegaC     = mOmegaC(1)             ! dimensionless
      BetaD      = mBetaD(1)              ! dimensionless
      rhoSW      = mrhoSW(1)              ! kg / m3
      p          = mp(1)                  ! db
      tempis     = mtempis(1)             ! degrees C
      pco2atm    = mpco2atm(1)            ! uatm
      fco2atm    = mfco2atm(1)            ! uatm
      schmidtco2 = mschmidtco2(1)         ! dimensionless
      kwco2      = mkwco2(1)              ! m / s
      K0         = mK0(1)                 ! (mol/kg) / atm
      co2starair = mco2starair(1) * 1000. ! mol / m3 -> mmol / m3 
      co2flux    = mco2flux(1) * 1000.    ! mol / m2 / s -> mmol / m2 / s 
      dpco2      = mdpco2(1)              ! uatm

  RETURN

  END SUBROUTINE

!-----------------------------------------------------------------------

!=======================================================================
!
      SUBROUTINE mocsy_carbchem( temp, sal, alk, dic, sil, phos,  &
           Patm, depth, lat, kw660, xco2, N,                      &
           ph, pco2, fco2, co2, hco3, co3, OmegaA,                &
           OmegaC, BetaD, rhoSW, p, tempis,                       &
           pco2atm, fco2atm, schmidtco2, kwco2, K0,               &
           co2starair, co2flux, dpco2 )
!
!=======================================================================
! 
! AXY (23/06/15): MOCSY introduced to MEDUSA to update carbonate
!                 chemistry for UKESM1 project
!
! Orr, J. C. and Epitalon, J.-M.: Improved routines to model the 
! ocean carbonate system: mocsy 2.0, Geosci. Model Dev., 8, 485-499, 
! doi:10.5194/gmd-8-485-2015, 2015.
!
! "mocsy is a Fortran 95 package designed to compute all ocean
! carbonate system variables from DIC and total Alk, particularly
! from models.  It updates previous OCMIP code, avoids 3 common 
! model approximations, and offers the best-practice constants as 
! well as more recent options. Its results agree with those from 
! CO2SYS-MATLAB within 0.005%."
!
! Where possible the code remains identical to that published by
! Orr & Epitalon (2015; henceforth OE15); some consolidation of 
! MOCSY files has taken place, and the resulting package is 
! comprised of four modules:
!                 
! 1. mocsy_wrapped.f90
!    - This module contains the interface between MEDUSA and the
!      main MOCSY routines; it draws from (but is very different
!      from) the test_mocsy.f90 routine provided by OE15
!
! 2. mocsy_singledouble.f90
!    - This module contains only precision definitions; it is
!      based on the singledouble.f90 routine provided by OE15
!
! 3. mocsy_phsolvers.f90
!    - This module contains a suite of solvers derived from
!      Munhoven (2013) and modified by OE15; it is based on the
!      phsolver.f90 routine provided by OE15
!
! 4. mocsy_gasflux.f90
!    - This module contains a series of subroutines used to 
!      calculate the air-sea flux of CO2; it consolidates four
!      MOCSY routines with an OCMIP-2 Schmidt number routine
!      updated with the new parameterisations of Wanninkhof (2014)
!
! 5. mocsy_mainmod.f90
!    - This module consolidates all of the remaining MOCSY 
!      functions and subroutines from OE15 into a single file
!      for convenience
!
! NOTE: it is still possible to run PML's carbonate chemistry
! routine in MEDUSA; at present this remains the default, if
! less-preferred, routine (i.e. is used if key_mocsy is not
! present)
!
!=======================================================================
!
! AXY (05/04/17): alter options to include optK1K2 = 'w14'
!
! In conversation with Jim Orr, Waters (2014) is now the 
! preferred option for optK1K2 in the case of global scale 
! simulations since this formulation works over broader ranges
! of temperature (0 < T < 50) and salinity (1 < S < 50).
!
! NOTE: *contrary* to the notice above, MEDUSA has now been
! revised to remove the PML carbonate chemistry routine as part
! of a wider code review. This routine is not anticipated to be 
! used again, and has several out of date parameterisations 
! relative to MOCSY. As a result, key_mocsy is no longer required.
!
        USE mocsy_singledouble
        IMPLICIT NONE

!> ======================================================================
!  VARIABLES
!> ======================================================================
!
!  Input variables
   REAL(kind=wp), INTENT(in),  DIMENSION(N) :: temp, sal, alk, dic, sil, phos
   REAL(kind=wp), INTENT(in),  DIMENSION(N) :: Patm, depth, lat, kw660, xco2
   INTEGER, INTENT(in) :: N
!
!  Output variables
   REAL(kind=wp), INTENT(out), DIMENSION(N) :: ph, pco2, fco2, co2, hco3, co3, OmegaA
   REAL(kind=wp), INTENT(out), DIMENSION(N) :: OmegaC, BetaD, rhoSW, p, tempis
   REAL(kind=wp), INTENT(out), DIMENSION(N) :: pco2atm, fco2atm, schmidtco2, kwco2, K0
   REAL(kind=wp), INTENT(out), DIMENSION(N) :: co2starair, co2flux, dpco2
!
!  Local variables
   REAL(kind=wp), DIMENSION(N) :: depth0, co2star
   INTEGER :: i, totdat
!
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2
!   
!  initialise depth0 to 0
   depth0 = 0.0
   
!> ======================================================================
!  CONFIGURE OPTIONS
!> ======================================================================
!> OPTIONS: see complete documentation in 'vars' subroutine (in mocsy_mainmod.F90)
!> AXY (05/04/17): optK1K2 switched from 'm10' to 'w14'
!> Typical options for MODELS
   optCON  = 'mol/m3' ! input concentrations are in MOL/M3
   optT    = 'Tpot'   ! input temperature, variable 'temp', is POTENTIAL temp [°C]
   optP    = 'm'      ! input variable 'depth' is in METERS
   optB    = 'l10'    ! Lee et al. (2010) formulation for total boron
   optK1K2 = 'w14'    ! Lueker et al. (2000) formulations for K1 & K2 (best practices)
   optKf   = 'dg'     ! Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)
!> optK1K2 = 'l'      ! Lueker et al. (2000) formulations for K1 & K2 (best practices)
!> optKf   = 'dg'     ! Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)

!> ======================================================================
!  CARBONATE CHEMISTRY CALCULATIONS
!> ======================================================================
!> Call mocsy's main subroutine to compute carbonate system variables: 
!> pH, pCO2, fCO2, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R
!> FROM temperature, salinity, total alkalinity, dissolved inorganic 
!> carbon, silica, phosphate, depth (or pressure) (1-D arrays)
   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                      &  ! INPUT
             optCON, optT, optP)
!            optCON, optT, optP, optB, optK1K2, optKf)                                    ! INPUT OPTIONS

!> ======================================================================
!  GAS EXCHANGE CALCULATIONS
!> ======================================================================
!>
!> Only calculate gas exchange fields if depth = 0
   if (depth(1) .eq. 0.0) then
!
!     Compute pCO2atm [uatm] from xCO2 [ppm], atmospheric pressure [atm], & vapor pressure of seawater
!     pCO2atm = (Patm - pH20(i)) * xCO2,   where pH20 is the vapor pressure of seawater [atm]
      CALL x2pCO2atm(xco2, temp, sal, Patm, N, & ! INPUT
      pco2atm)                                   ! OUTPUT

!     Compute fCO2atm [uatm] from pCO2atm [uatm] & fugacity coefficient [unitless]
!     fCO2atm = pCO2atm * fugcoeff,   where fugcoeff= exp(Patm*(B + 2.0*xc2*Del)/(R*tk) )
      CALL p2fCO2(pco2atm, temp, Patm, depth0, N, & ! INPUT
      fco2atm)                                      ! OUTPUT

!     Compute Schmidt number for CO2 from potential temperature
      CALL schmidt_co2(temp, N, & ! INPUT
      schmidtco2)                 ! OUTPUT

!     Compute transfer velocity for CO2 in m/s (see equation [4] in OCMIP2 design document & OCMIP2 Abiotic HOWTO)
      kwco2 = kw660 * (660./schmidtco2)**0.5

!     Surface K0 [(mol/kg) / atm] at T, S of surface water
      CALL surface_K0(temp, sal, N, & ! INPUT
      K0)                             ! OUTPUT

!     "Atmospheric" [CO2*], air-sea CO2 flux, sfc DIC rate of change, & Delta pCO2
!     all "lifted" from the gasx.f90 function of MOCSY
      co2starair = K0 * fco2atm * 1.0e-6_wp * rhoSW ! Equil. [CO2*] for atm CO2 at Patm & sfc-water T,S [mol/m3]
      co2star    = co2                              ! Oceanic [CO2*] in [mol/m3] from vars.f90
      co2flux    = kwco2 * (co2starair - co2star)   ! Air-sea CO2 flux [mol/(m2 * s)]
!     co2ex      = co2flux / dz1                    ! Change in sfc DIC due to gas exchange [mol/[m3 * s)]
      dpco2      = pco2 - pco2atm                   ! Delta pCO2 (oceanic - atmospheric pCO2) [uatm]
!     
   endif

  RETURN

  END SUBROUTINE

!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

#  else
   !!======================================================================
   !!  Dummy module :                         No MOCSY carbonate chemistry
   !!======================================================================

CONTAINS

   SUBROUTINE mocsy_interface( kt )                   ! Empty routine

      INTEGER, INTENT( in ) ::   kt

      WRITE(*,*) 'mocsy_interface: You should not have seen this print! error?', kt

   END SUBROUTINE mocsy_interface
#endif 

   !!======================================================================

END MODULE mocsy_wrapper
