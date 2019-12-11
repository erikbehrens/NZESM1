MODULE trcco2_medusa
   !!======================================================================
   !!                         ***  MODULE trcco2_medusa  ***
   !! TOP :   MEDUSA
   !!======================================================================
   !! History :
   !!  -   !  2010-12  (A. Yool)             added for ROAM project
   !!----------------------------------------------------------------------
#if defined key_medusa && defined key_roam
   !!----------------------------------------------------------------------
   !!                                        MEDUSA carbonate chemistry
   !!----------------------------------------------------------------------
   !!   trc_co2_medusa        :  
   !!----------------------------------------------------------------------
      USE oce_trc
      USE trc
      USE sms_medusa
      USE lbclnk
      USE prtctl_trc      ! Print control for debugging
      USE in_out_manager  ! I/O manager

      IMPLICIT NONE
      PRIVATE
      
      PUBLIC   trc_co2_medusa    ! called in trc_bio_medusa

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

! The following is a map of the subroutines contained within this module
! - trc_co2_medusa
!      - CALLS CO2_dynamics
!           - CALLS CO2DYN
!		 - CALLS POLYCO
!		      - CALLS CO2SET
!		      - CALLS CO2CLC
!           - CALLS CaCO3_Saturation
!      - CALLS Air_sea_exchange

CONTAINS

!=======================================================================
!
      SUBROUTINE trc_co2_medusa( Temp, Sal, DIC, ALK, Depth, xkw, pCO2a, &
      pH, pCO2w, h2co3, hco3, co3, om_cal, om_arg, co2flux, TDIC, TALK,  &
      dcf, henry, iters )
!      
!=======================================================================
!  
! This file contains a set of FORTRAN subroutines that calculate the carbonate system 
! at any given point in marine space time, given values for 
! temperature, salinity, DIC, depth (pressure). 
! This is essentially an implimentation of the Haltafall speciation code
! (Ingri et al 1967, Talanta 14, 1261 - if it ain't broke don't fix it)
! Another routine calulates the air sea exchange of CO2 given wind speed and atmospheric pCO2.
! Code developed by Jerry blackford and others at PML, based on pre-existing code.
! We accept no liability for errors or inaccuracies.
! See Zeebe & Wolf-Gladrow, 2001. CO2 in seawater: equilibrium, kinetics and isotopes. 
! Elsevier Oceanography Series 65, 346. for a reasonable overview. 
! Many other packages exist, replicating the same functionality in different languages.
! See http://cdiac.ornl.gov/oceans/co2rprt.html (CO2sys)
! or  http://neon.otago.ac.nz/research/mfc/people/keith_hunter/software/swco2/
! reference for prior usage of this code: Blackford & Gilbert, 2007. J Mar Sys 64, 229-241.
! 
! Modifications
! 17/02/2010. Added conversion factor from per m3 to per kg (line 108-133)
! 17/02/2010. Update calculation of K1, K2, Kb to make consistant with the OCMIP protocols.
! 29/07/2011. Merged into MEDUSA with a raft of changes to this subroutine; less elsewhere
! 23/06/2015. Modified to take gas transfer velocity as an input (rather than wind speed); 
!             alter CO2 flux to /s rather than /d for consistency with other schemes
!
! Changes for MEDUSA include: 
!     - making the program a module
!     - passing in input variables (obvious given the last point)
!     - making alkalinity a state variable (rather than a function of salinity)
!
      IMPLICIT NONE

      REAL(wp), INTENT( in )    :: Temp       ! degrees C
      REAL(wp), INTENT( in )    :: Sal        ! PSU
      REAL(wp), INTENT( in )    :: DIC        ! mmol / m3
      REAL(wp), INTENT( in )    :: ALK        ! meq  / m3
      REAL(wp), INTENT( in )    :: Depth      ! m
!     REAL(wp), INTENT( in )    :: Wnd        ! m / s
      REAL(wp), INTENT( in )    :: xkw        ! m / s
      REAL(wp), INTENT( in )    :: pCO2a      ! uatm
!----------------------------------------------------------------------
      REAL(wp), INTENT( inout ) :: pH         ! "normal" pH units
      REAL(wp), INTENT( inout ) :: pCO2w      ! uatm
      REAL(wp), INTENT( inout ) :: h2co3      ! mmol / m3
      REAL(wp), INTENT( inout ) :: hco3       ! mmol / m3
      REAL(wp), INTENT( inout ) :: co3        ! mmol / m3
      REAL(wp), INTENT( inout ) :: om_cal     ! normalised
      REAL(wp), INTENT( inout ) :: om_arg     ! normalised
      REAL(wp), INTENT( inout ) :: co2flux    ! mmol / m2 / s
      REAL(wp), INTENT( inout ) :: TDIC       ! umol / kg
      REAL(wp), INTENT( inout ) :: TALK       ! ueq  / kg
      REAL(wp), INTENT( inout ) :: dcf        ! m3 / kg
      REAL(wp), INTENT( inout ) :: henry      ! ?
      INTEGER,  INTENT( inout ) :: iters      ! # iterations to convergence
!----------------------------------------------------------------------
!     REAL(wp) :: cco2,bicarb,carb,henry
      REAL(wp) :: cco2,bicarb,carb
! ======================================================================

! write inputs to screen
!  WRITE(*,*) " "
!  WRITE(*,'(A28)')       "    .........Inputs........."
!  WRITE(*,'(A24,F10.3)') "    Temperature   (C) = ", Temp
!  WRITE(*,'(A24,F10.3)') "    Salinity    (psu) = ", Sal
!  WRITE(*,'(A24,F10.3)') "    Depth         (m) = ", Depth
!  WRITE(*,'(A24,F10.3)') "    DIC     (mmol/m3) = ", DIC
!  WRITE(*,'(A24,F10.3)') "    ALK      (ueq/m3) = ", ALK
!  WRITE(*,'(A24,F10.3)') "    Wind speed  (m/s) = ", Wnd
!  WRITE(*,'(A24,F10.3)') "    pCO2 atmos (uatm) = ", pCO2a

! AXY (07/05/13) ==================================================
! set total number of iterations to zero
      iters = 0
! AXY (07/05/13) ==================================================
  
   call CO2_dynamics( Temp, Sal, Depth, DIC, ALK, pCO2a, &                  ! inputs
      pco2w, ph, h2co3, hco3, co3, henry, om_cal, om_arg, TDIC, TALK, &     ! outputs
      dcf, iters)                                                           ! outputs

   ! AXY (18/08/11): only do air-sea exchange calculation if depth = 0 
   !                 (i.e. surface calculations being performed)
   if (Depth .eq. 0.0) then
      call Air_sea_exchange( Temp, xkw, pCO2w, pCO2a, henry, dcf, &         ! inputs
         co2flux )                                                          ! output
   else
      co2flux = 0.0
   endif      

! AXY (09/01/14) ==================================================
! check for non-convergence and NaNs
   !!
   if (iters .eq. 25) then
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: ITERS WARNING, ', &
      iters, ' AT depth =', Depth
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: ztmp    =', Temp
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: zsal    =', Sal
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: zdic    =', DIC
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: zalk    =', ALK
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: f_kw660 =', xkw
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: f_ph    =', ph
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: f_pco2w =', pCO2w
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: f_h2co3 =', h2co3
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: f_hco3  =', hco3
      IF(lwp) WRITE(numout,*) ' trc_co2_medusa: f_co3   =', co3
   endif
   !!
   !! AXY (29/10/13): NaN checks
   if (co2flux /= co2flux) then
      !! if (lwp) write (numout,'(a,1pe10.2,4i6)') 'CO2FLUX-NAN', &
      !! &        tmask(ji,jj,jk), kt, ji, jj, jk
      if (lwp) write (numout,'(a,a,f10.3,a,f10.3)') 'CO2FLUX-NAN', &
      &        ' TMP', Temp, ' SAL', Sal
      if (lwp) write (numout,'(a,a,f10.3,a,f10.3)') 'CO2FLUX-NAN', &
      &        ' DIC', DIC, ' ALK', ALK
      if (lwp) write (numout,'(a,a,f10.3,a,f10.3)') 'CO2FLUX-NAN', &
      &        ' XKW', xkw, ' PH ', ph
      if (lwp) write (numout,'(a,a,i6)') 'CO2FLUX-NAN', &
      &        ' ITERS', iters
   endif
   !!
   !! AXY (09/01/14): NaN fudges
   if (co2flux /= co2flux) then
      ph      = 8.1
      pCO2w   = pCO2a
      h2co3   = 0.005 * DIC
      hco3    = 0.865 * DIC
      co3     = 0.130 * DIC
      om_cal  = 4.
      om_arg  = 2.
      co2flux = 0.
      TDIC    = DIC
      TALK    = ALK
      dcf     = 1.
      henry   = 1.
   endif
! AXY (09/01/14) ==================================================

! write outputs to screen
!  WRITE(*,*) " "
!  WRITE(*,'(A32,F10.3)') "    ..........Outputs..........."
!  WRITE(*,'(A27,F10.3)') "    pH                 (pH) = ", pH
!  WRITE(*,'(A27,F10.3)') "    DIC           (umol/kg) = ", TDIC
!  WRITE(*,'(A27,F10.3)') "    TALK           (ueq/kg) = ", TALK
!  WRITE(*,'(A27,F10.3)') "    pco2w            (uatm) = ", pco2w
!  WRITE(*,'(A27,F10.3)') "    carbonic acid (mmol/m3) = ", h2co3
!  WRITE(*,'(A27,F10.3)') "    bicarbonate   (mmol/m3) = ", bicarb
!  WRITE(*,'(A27,F10.3)') "    carbonate     (mmol/m3) = ", carb
!  WRITE(*,'(A27,F10.3)') "    Omega calcite       (~) = ", om_cal
!  WRITE(*,'(A27,F10.3)') "    Omega aragonite     (~) = ", om_arg
!  WRITE(*,'(A27,F10.3)') "    air sea flux(mmol/m2/s) = ", flux
!  WRITE(*,*) " "

   END SUBROUTINE trc_co2_medusa
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
      subroutine CO2_dynamics( T, S, Z, DIC, TALK, pco2a, &
      pco2w, ph, h2co3, bicarb, carb, henry, om_cal, om_arg, TCO2, TA, &
      dcf, iters )
!      
!=======================================================================
!  
   IMPLICIT NONE

      REAL(wp), INTENT( in )    :: T          ! temperature (C)
      REAL(wp), INTENT( in )    :: S          ! salinity (psu)
      REAL(wp), INTENT( in )    :: Z          ! depth (metres)
      REAL(wp), INTENT( in )    :: DIC        ! total dissolved inorganic carbon (mmol.m-3)
      REAL(wp), INTENT( in )    :: TALK       ! total alkalinity (meq.m-3)
      REAL(wp), INTENT( in )    :: pco2a      ! atmospheric pCO2
!----------------------------------------------------------------------
      REAL(wp), INTENT( inout ) :: pco2w      ! seawater pCO2
      REAL(wp), INTENT( inout ) :: ph         ! seawater pH
      REAL(wp), INTENT( inout ) :: h2co3      ! seawater carbonic acid (H2CO3)
      REAL(wp), INTENT( inout ) :: bicarb     ! seawater bicarbonate ion (HCO3)
      REAL(wp), INTENT( inout ) :: carb       ! seawater carbonate ion (CO3)
      REAL(wp), INTENT( inout ) :: henry      ! Henry constant
      REAL(wp), INTENT( inout ) :: om_cal     ! Omega calcite
      REAL(wp), INTENT( inout ) :: om_arg     ! Omega aragonite
      REAL(wp), INTENT( inout ) :: TCO2       ! total dissolved inorganic carbon (mol.kg-1)
      REAL(wp), INTENT( inout ) :: TA         ! total alkalinity (eq.kg-1)
      REAL(wp), INTENT( inout ) :: dcf        ! density conversion factor
      INTEGER,  INTENT( inout ) :: iters      ! # iterations to convergence
!----------------------------------------------------------------------
      REAL(wp) :: a, b, c
      REAL(wp) :: ca, bc, cb
      REAL(wp) :: pco2water, fairco2

! ======================================================================

! Adjust to correct units.
! Haltafall uses mol/kg rather than umol/kg necessitating a scaling factor of /1.0D6
! DIC (mmol/m3) needs to be converted to umol/kg via the calculation of water density at prevailing T&S
! sea-water density (Millero & Poisson, Deep-Sea Research, 1981, also known as UNESCO, 1981)
! with T: Temperature in degree Celsius; S: Salinity in practical units, density in kg/m3
! valid for 0<T<40 and 0.5<S<43. Density Conversion factor (dcf) = * density * 1.0D3

          a   =  8.24493d-1 - 4.0899d-3*T +  7.6438d-5*T**2 - 8.2467d-7*T**3 + 5.3875d-9*T**4
          b   = -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T**2 
          c   = 4.8314d-4
          dcf = (999.842594 + 6.793952d-2*T- 9.095290d-3*T**2 + 1.001685d-4*T**3 &
                - 1.120083d-6*T**4 + 6.536332d-9*T**5+a*S+b*S**1.5+c*S**2)/1.0D3

          TA    = TALK / (1.0D6*dcf)
          TCO2  = DIC  / (1.0D6*dcf)

! Call the parent routine for the carbonate system
          CALL CO2DYN ( TCO2, TA, T, S, pco2a, &     ! inputs
          pco2water, pH, HENRY, ca, bc, cb, iters )  ! outputs

! Adjust outputs back to units used in the parent model code (e.g. mmol/m3) if appropriate

          pco2w  = pco2water * (1.0D6) 	! partial pressure of co2 in water 
          TA     = TA * (1.0D6)         ! total alkalinity (umol/kg)
          h2co3  = ca * (1.0D6*dcf)     ! carbonic acid concentration (mmol/m3)
          bicarb = bc * (1.0D6*dcf)     ! bicarbonate ion concentration (mmol/m3)
          carb   = cb * (1.0D6*dcf)     ! carbonate ion concentration (mmol/m3)
          TCO2   = TCO2 * (1.0D6)       ! total C or DIC in units of umol/kg

! Call carbonate saturation state subroutine to calculate calcite and aragonite calcification states

          CALL CaCO3_Saturation ( T, S, Z, cb, &  ! inputs
          om_cal, om_arg )                        ! outputs

   RETURN

   END SUBROUTINE CO2_dynamics
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
      SUBROUTINE Air_sea_exchange( T, xkw, pco2w, pco2a, henry, dcf, &
      flux )
!      
!=======================================================================
!  
!  This routine should be called for the surface box only (like, duh)
!  Uses the Nightingale and Liss parameterisation (GBC 14, 373-388 (2000).
!  SC is the Schmidt number from Wanninkhof (1992, JGR 97,7373-7382), the Nightingale et al (2000)
!  transfer velocity being for Sc=600 and the dependence of tranfer velocity on Sc coming 
!  from Jahne et al (1987, J. Geophys Res 92, 1937-1949).

!  Inputs
!  pCO2w    partial pressure of CO2 in the water (from carbonate system subroutine call)
!  pCO2a    partial pressure of CO2 in the atmosphere (usually external forcing).
!  T        temperature (C)
!  Wnd      wind speed, metres (DELETED)
!  xkw      gas transfer velocity
!  Henry    henry's constant
!  density  the density of water for conversion between mmol/m3 and umol/kg

!  Outputs are
!  flux     flux of CO2 in mmol C /m2/d 
!           +ve is in-gassing (air to sea), -ve is outgassing (sea to air).

   IMPLICIT NONE

      REAL(wp), INTENT( in )    :: T, xkw, pco2w, pco2a, henry, dcf ! INPUT PARAMETERS:
!-----------------------------------------------------------------------
      REAL(wp), INTENT( inout ) :: flux                             ! OUTPUT Variables
!-----------------------------------------------------------------------
      REAL(wp)                  :: sc, fwind                        ! LOCAL VARIABLES:

! calculate the Schmidt number and unit conversions
          sc    = 2073.1-125.62*T+3.6276*T**2.0-0.0432190*T**3.0
!         fwind = (0.222d0 * wnd**2d0 + 0.333d0 * wnd)*(sc/660.d0)**(-0.5)
          fwind = xkw * (sc/660.d0)**(-0.5)
          fwind = fwind*24.d0/100.d0   ! convert to m/day

! flux depends on the difference in partial pressures, wind and henry
! here it is rescaled to mmol/m2/d
          flux = fwind * henry * ( pco2a - pco2w ) * dcf

! AXY (23/06/15): let's get it from /d to /s
          flux = flux / ( 86400. )

  RETURN 

  END SUBROUTINE Air_sea_exchange
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
      SUBROUTINE CO2dyn ( TCO2, TA, T, S, PCO2A, PCO2, PH, HENRY, ca, bc, cb, iters )
!      
!=======================================================================
!
!     This subroutine acts as an interface to the Haltafall iteration, setting options etc.

      REAL(wp), INTENT( in )    :: TCO2, TA, T, S, PCO2A
!-----------------------------------------------------------------------
      REAL(wp), INTENT( inout ) :: PCO2, PH, HENRY, ca, bc, cb
      INTEGER,  INTENT( inout ) :: iters
!-----------------------------------------------------------------------
      REAL(wp) :: PRSS
      INTEGER  :: MCONC, MKVAL, ICONST, ICALC

      PARAMETER ( MCONC = 9,MKVAL = 4 )

      REAL(wp), DIMENSION(MKVAL) :: AKVAL
      REAL(wp), DIMENSION(MCONC) :: CONCS

      ICONST = 6
      PRSS = 1.0D0
      CONCS(1) = TCO2
      CONCS(2) = TA
      ICALC = 1
      
      CALL POLYCO(PRSS,T,S,CONCS,MCONC,AKVAL,MKVAL,ICALC,ICONST,iters)

      if(iters.eq.25) then
         ! AXY (13/05/13): this location has failed to converge; use
         !                 some standard values to plug the hole
         CONCS(3) = PCO2A * 1E-06 ! atmospheric value
         CONCS(4) = 8.1           ! global average pH
         CONCS(5) = TCO2 * 0.005  ! some "standard" values plucked
         CONCS(6) = TCO2 * 0.865  ! from pages 5-6 of Zeebe & Wolf-
         CONCS(7) = TCO2 * 0.130  ! Gladow (2001)
         AKVAL(1) = 0.5E-01       ! trial-and-error value
         IF(lwp) THEN
            WRITE(numout,*) ' CO2DYN : reset failed outputs'
         ENDIF
      endif

      PCO2  = CONCS(3)
      PH    = CONCS(4)
      ca    = CONCS(5)
      bc    = CONCS(6)
      cb    = CONCS(7)
      HENRY = AKVAL(1)
    
      RETURN

      END SUBROUTINE
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
      SUBROUTINE POLYCO(PD,TD,SD,CONCS,NCONC,AKVAL,NKVAL,ICALC,ICONST,iters)
!      
!=======================================================================
!
! MASTER SUBROUTINE FOR CALCULATION OF THE CO2 SYSTEM THERMODYNAMICS
!
! EXPLANATION OF POLYCO PARAMETERS
!       AXY (16/08/11): Yuri change for non-surface carbonate chemistry 
!       P - PRESSURE IN ATMOSPHERES (P<>1 CODED only for ICONST=3 and ICONST=6)
!       T - TEMPERATURE IN DEG.C
!       S - SALINITY IN PPT
!      CONCS(1) - TOTAL C (MOL/KG)
!      CONCS(2) - TOTAL ALKALINITY (MOL/KG)
!      CONCS(3) - PCO2 (ATM)
!      CONCS(4) - PH
!      CONCS(5) - {H2CO3} (MOL/KG)
!      CONCS(6) - {HCO3} (MOL/KG)
!      CONCS(7) - {CO3} (MOL/KG)
!        CONCS(8) - CARBONATE ALKALINITY  ) FOR ICONST = 4,5,6
!        CONCS(9) - BORATE ALKALINITY     )       ONLY
!         NCONC - SIZE OF CONCS ARRAY (7 FOR ICONST=1,2,3; 9 FOR ICONST
!     AKVAL(1) - KP (HENRY'S LAW CONSTANT) (MOL/KG/ATM)
!     AKVAL(2) - K1C (H2CO3 DISSOCIATION) (MOL/KG)
!       AKVAL(3) - K2C (HCO3 DISSOCIATION) (MOL/KG)
!       AKVAL(4) - KB (B(OH)3 DISSOCIATION) (MOL/KG)  FOR ICONST=4,5,6
!        NKVAL - SIZE OF AKVAL ARRAY (3 FOR ICONST=1,2,3; 4 FOR ICONST=
!        ICALC - SELECTION OF THE TWO INPUT PARAMETERS:
!       ICALC = 1  TOTAL C AND ALKALINITY
!       ICALC = 2  TOTAL C AND PCO2
!       ICALC = 3  TOTAL C AND PH
!       ICALC = 4  ALKALINITY AND PCO2
!       ICALC = 5  ALKALINITY AND PH
!       ICALC = 6  PCO2 AND PH
!       ICALC = 7  CALCULATE CONSTANTS AKVAL ONLY
!       ICONST - SELECTION OF PH SCALE AND COMPONENTS:
!       ICONST = 1  NBS PH SCALE
!       ICONST = 2  HANSSON'S SCALE (SWS WITHOUT FLUORIDE)
!       ICONST = 3  SWS PH SCALE
!       ICONST = 4  AS 1 BUT INCLUDING BORATE IN THE CALCULATION
!       ICONST = 5  AS 2 BUT INCLUDING BORATE IN THE CALCULATION
!       ICONST = 6  AS 3 BUT INCLUDING BORATE IN THE CALCULATION

!  NOTE: FOR ICONST=1,2,3 CONCS(2) REPRESENTS CARBONATE ALKALINITY SINC
!        BORATE IS NOT INCLUDED IN THE CALCULATION. FOR ICONST=4,5,6 CO
!        REPRESENTS TOTAL ALKALINITY (CARBONATE + BORATE), THE COMPONEN
!        WHICH ARE GIVEN IN CONCS(8) AND CONCS(9)

      REAL(wp) :: PMIN, PMAX, SMIN, SMAX, TMIN, TMAX, & 
     &     PD, TD, SD, P, T, S, BTOT
      INTEGER  :: MINJC, MAXJC, MINJK, MAXJK, MINCAL, MAXCAL, MINCON,  &
     &     MAXCON, NCONC, NKVAL, ICALC, ICONST, IC, iters
      LOGICAL  :: BORON

      PARAMETER(MINJC=7,MAXJC=9,MINJK=3,MAXJK=4)
      PARAMETER(MINCAL=1,MAXCAL=7,MINCON=1,MAXCON=6)
      !! AXY (11/08/11): TMIN changed to -2 to stop error messages in polar regions
      !! AXY (09/01/14): TMIN changed to -3 to stop error messages in polar regions
      !! AXY (03/03/14): SMAX changed to 42 to stop error messages in salty regions (the Great Answer ...)
      !! AXY (05/03/14): SMAX changed to 45 to stop error messages in salty regions
      PARAMETER(PMIN=0.99999D0,PMAX=1.00001D0,SMIN=0.0D0,  &
     &  SMAX=45.0D0,TMIN=-3.0D0,TMAX=40.0D0)
      REAL(wp), DIMENSION(NKVAL) :: AKVAL
      REAL(wp), DIMENSION(NCONC) :: CONCS

      P = PD
      S = SD
      T = TD
      
      !! AXY (09/08/11): STOP statements commented out and WARNING messages added
      IF(lwp) THEN
         IF(T.LT.TMIN.OR.T.GT.TMAX) WRITE(numout,*) ' trc_co2_medusa: T WARNING, ', T, TMIN, TMAX, S, P
         IF(S.LT.SMIN.OR.S.GT.SMAX) WRITE(numout,*) ' trc_co2_medusa: S WARNING, ', S, SMIN, SMAX, T, P
         IF(P.LT.PMIN.OR.P.GT.PMAX) WRITE(numout,*) ' trc_co2_medusa: P WARNING, ', P, PMIN, PMAX, T, S
      ENDIF
      ! IF(T.LT.TMIN.OR.T.GT.TMAX)WRITE (*,*) P, S, T, TMIN, TMAX
      ! IF(P.LT.PMIN.OR.P.GT.PMAX) STOP'POLYCO - PRESSURE OUT OF RANGE'
      ! IF(S.LT.SMIN.OR.S.GT.SMAX) STOP'POLYCO - SALINITY OUT OF RANGE'
      ! IF(T.LT.TMIN.OR.T.GT.TMAX) STOP'POLYCO - TEMP. OUT OF RANGE'

      !! AXY (17/04/13): iMARNET climate change simulation appears to be compromised
      !!                 by excessive Caspian Sea salinity (> 45 PSU); this may be
      !!                 a runoff glitch that doesn't manifest on local instance of
      !!                 NEMO; fudging this here by replacing excessively high
      !!                 salinity values with maximum salinity value for the purpose
      !!                 of carbonate chemistry calculations; if this works, a more
      !!                 sensible solution might be to kill Caspian, etc.
      IF(lwp) THEN
         IF(S.LT.SMIN) THEN
            WRITE(numout,*) ' trc_co2_medusa: S RESET, ', S, '->', SMIN
            S = SMIN
         ENDIF
         IF(S.GT.SMAX) THEN
            WRITE(numout,*) ' trc_co2_medusa: S RESET, ', S, '->', SMAX
            S = SMAX
         ENDIF
      ENDIF
      
      !! AXY (28/02/13): as above, but for T; this has been done for 1/4-degree
      !!                 simulations to provide a temporary fix to a strange
      !!                 glitch in ocean temperature; basically, T has been found
      !!                 to spike upwards arbitrarily (e.g. 24C -> 53C); this 
      !!                 causes the carbonate chemistry calculations to fail;
      !!                 this fix aims to stop the problem so that we can estimate
      !!                 how frequent a problem this is; it's suggestive of a
      !!                 memory leak
      IF(lwp) THEN
         IF(T.LT.TMIN) THEN
            WRITE(numout,*) ' trc_co2_medusa: T RESET, ', T, '->', TMIN
            T = TMIN
         ENDIF
         IF(T.GT.TMAX) THEN
            WRITE(numout,*) ' trc_co2_medusa: T RESET, ', T, '->', TMAX
            T = TMAX
         ENDIF
      ENDIF
     
      IF(ICALC.LT.MINCAL.OR.ICALC.GT.MAXCAL) STOP 'POLYCO - ICALC OUT OR RANGE'
      IF(ICONST.LT.MINCON.OR.ICONST.GT.MAXCON) STOP 'POLYCO - ICONST OUT OF RANGE'
      BORON=(ICONST.GT.3)
      IF(BORON) THEN
        IC=ICONST-3
        BTOT=0.0004128D0*S/35.0D0
        IF(NCONC.NE.MAXJC) STOP 'POLYCO - WRONG NCONC VALUE'
        IF(NKVAL.NE.MAXJK) STOP 'POLYCO - WRONG NKVAL VALUE'
      ELSE
        IC=ICONST
        IF(NCONC.NE.MINJC) STOP 'POLYCO - WRONG NCONC VALUE'
        IF(NKVAL.NE.MINJK) STOP 'POLYCO - WRONG NKVAL VALUE'
      ENDIF

      CALL CO2SET(P,T,S,AKVAL,NKVAL,IC)

      IF(ICALC.LT.MAXCAL)  &
     & CALL CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT,iters)

      RETURN

      END SUBROUTINE
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
      SUBROUTINE CO2SET(P,T,S,AKVAL,NKVAL,IC)
!      
!=======================================================================
!
! Routine to calculate CO2 system constants under the conditions set by
! P,S,T     (NOTE: PRESSURE <> 1ATM IS NOT YET CODED)

! I. Calculate constants at P=1 and S=0 using

!      ln K0  =  A + B/TK + C ln TK
!                                   (where TK is in Kelvin)

! II. Calculate constants at P=1 and salinity S using

!    ln K  =  ln K0 + (a0 + a1/TK + a2 ln TK) S**1/2
!                 + (b0 + b1TK + b2TK**2) S

! The sources of the coefficients are as follows:

!  IC=                  1                    2                  3
!               (NBS pH scale)        (SWS pH scale       (SWS pH scale
!                                       with no F)

!  KP            WEISS (1974)           WEISS(1974)          WEISS(1974

!  K1C )      MEHRBACH ACC. TO       HANSSON ACC. TO    HANSSON AND MEH
!  K2C )       MILLERO (1979)         MILLERO (1979)      ACC. TO DICKS
!   KB )                                                 AND MILLERO (1
!                                                         (K1C AND K2C
!                                                           HANSSON ACC
!                                                            MILLERO (1
!                                                              (KB ONLY

! ***
!      IMPLICIT real*8 (A-H,O-Z)

!     Modified by jcb 17/02/10 to use OCMIP calculations of K1, K2, Kb. 
!     Differences are subtle rather than significant

      INTEGER             :: MAXK, MAXCON, NKVAL, ICON, IC, IK
! ***
      PARAMETER(MAXK=4,MAXCON=3)
      REAL(wp), DIMENSION(MAXK)        :: A, B, C
      REAL(wp), DIMENSION(MAXK,MAXCON) :: A0, A1, A2
      REAL(wp), DIMENSION(MAXK,MAXCON) :: B0, B1, B2
      REAL(wp), DIMENSION(NKVAL)       :: AKVAL
!     AXY (16/08/11): Yuri change for non-surface carbonate chemistry 
!     REAL(wp)            :: P,T,S,VAL,TK
      REAL(wp)            :: P,T,S,VAL,TK, delta, kappa, Rgas
      REAL(wp)            :: dlogTK, S2, sqrtS, S15, k1, k2, kb
! ***
!     AXY (16/08/11): Yuri change for non-surface carbonate chemistry 
      DATA Rgas/83.131/

      DATA A/-167.8108D0, 290.9097D0, 207.6548D0, 148.0248D0/
      DATA B/9345.17D0, -14554.21D0, -11843.79D0, -8966.9D0/
      DATA C/23.3585D0, -45.0575D0, -33.6485D0, -24.4344D0/
      DATA (A0(1,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (A0(2,ICON),ICON=1,MAXCON) /0.0221D0, 0.5709D0, -45.8076D0/
      DATA (A0(3,ICON),ICON=1,MAXCON) /0.9805D0, 1.4853D0, -39.5492D0/
      DATA (A0(4,ICON),ICON=1,MAXCON) /0.0473D0, 0.5998D0, 0.5998D0/
      DATA (A1(1,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (A1(2,ICON),ICON=1,MAXCON) /34.02D0, -84.25D0, 1935.07D0/
      DATA (A1(3,ICON),ICON=1,MAXCON) /-92.65D0, -192.69D0, 1590.14D0/
      DATA (A1(4,ICON),ICON=1,MAXCON) /49.10D0, -75.25D0, -75.25D0/
      DATA (A2(1,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (A2(2,ICON),ICON=1,MAXCON) /2*0.0D0,6.9513D0/
      DATA (A2(3,ICON),ICON=1,MAXCON) /2*0.0D0,6.1523D0/
      DATA (A2(4,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B0(1,ICON),ICON=1,MAXCON) /3*0.023517D0/
      DATA (B0(2,ICON),ICON=1,MAXCON) /0.0D0,-0.01632D0,-0.01566D0/
      DATA (B0(3,ICON),ICON=1,MAXCON) /-0.03294D0,-0.05058D0,-0.04997D0/
      DATA (B0(4,ICON),ICON=1,MAXCON) /0.0D0, -0.01767D0, -0.01767D0/
      DATA (B1(1,ICON),ICON=1,MAXCON) /3*-2.3656D-4/
      DATA (B1(2,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B1(3,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B1(4,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B2(1,ICON),ICON=1,MAXCON) /3*4.7036D-7/
      DATA (B2(2,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B2(3,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B2(4,ICON),ICON=1,MAXCON) /3*0.0D0/

      TK=T+273.15D0
      DO 100 IK=1,NKVAL
        VAL=A(IK) + B(IK)/TK + C(IK)*LOG(TK)
        VAL=VAL + (A0(IK,IC) + A1(IK,IC)/TK + A2(IK,IC)*LOG(TK))*SQRT(S)
        VAL=VAL + (B0(IK,IC) + B1(IK,IC)*TK + B2(IK,IC)*TK*TK)*S
        AKVAL(IK)=EXP(VAL)
100    CONTINUE

      IF (IC .EQ. 3) THEN
!  Calculation of constants as used in the OCMIP process for ICONST = 3 or 6
!  see http://www.ipsl.jussieu.fr/OCMIP/
!  added jcb 17/02/10

!  Derive simple terms used more than once
	dlogTK = log(TK)
	S2 = S*S
	sqrtS = sqrt(S)
	S15 = S**1.5
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
	k1=10**(-1*(3670.7/TK - 62.008 + 9.7944*dlogTK - &
     &		0.0118 * S + 0.000116*S2))
	k2=10**(-1*(1394.7/TK + 4.777 - &
     &		0.0184*S + 0.000118*S2))
! AXY (16/08/11): Yuri change for non-surface carbonate chemistry 
! Correction for high pressure (from Millero 1995)
! added YA 04/10/2010
        delta=-25.5+0.1271*T
        kappa=(-3.08+0.0877*T)/1000.0
        k1=k1*exp((-delta+0.5*kappa*P)*P/(Rgas*TK))
        delta=-15.82-0.0219*T
        kappa=(1.13-0.1475*T)/1000.0
        k2=k2*exp((-delta+0.5*kappa*P)*P/(Rgas*TK))

! kb = [H][BO2]/[HBO2]
! Millero p.669 (1995) using data from Dickson (1990)
	kb=exp((-8966.90 - 2890.53*sqrtS - 77.942*S + &
     &		1.728*S15 - 0.0996*S2)/TK + &
     &		(148.0248 + 137.1942*sqrtS + 1.62142*S) + &
     &		(-24.4344 - 25.085*sqrtS - 0.2474*S) * &
     &		dlogTK + 0.053105*sqrtS*TK)
! AXY (16/08/11): Yuri change for non-surface carbonate chemistry 
! Correction for high pressure (from Millero, 1995)
! added YA 04/10/2010
        delta=-29.48+0.1622*T-0.002608*T**2.0
        kappa=-2.84/1000.0
        kb=kb*exp((-delta+0.5*kappa*P)*P/(Rgas*TK))
! Correction for high pressure of Kw (from Millero 1995)
! added YA 04/10/2010
        delta=-25.60+0.2324*T-0.0036246*T**2
        kappa=(-5.13+0.0794*T)/1000.0
        AKVAL(1)=AKVAL(1)*exp((-delta+0.5*kappa*P)*P/(Rgas*TK))

! replace haltafall calculations with OCMIP calculations
      AKVAL(2) = k1
      AKVAL(3) = k2
      AKVAL(4) = kb
      END IF ! section implimenting OCMIP coefficients

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
      SUBROUTINE CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT,iters)
!      
!=======================================================================
! 
! ROUTINE TO CARRY OUT CO2 CALCULATIONS WITH 2 FIXED PARAMETERS ACCORDI
! THE EQUATIONS GIVEN BY PARKS(1969) AND SKIRROW (1975)
! WITH ADDITIONS FOR INCLUDING BORON IF BORON=.TRUE.


!      IMPLICIT real*8 (A-H,O-Z)
      INTEGER             :: NCONC, NKVAL, ICALC, II, KARL, LQ, iters
! AXY (07/05/13) ==================================================
! put counter in to check duration in convergence loop
      INTEGER             :: COUNTER,C_CHECK,C_SW,III
! AXY (07/05/13) ==================================================
      REAL(wp)            :: CTOT,ALK,PCO2,PH,H2CO3,HCO3,CO3,ALKC
      REAL(wp)            :: ALKB,AKP,AK1C,AK2C,AKB,BTOT
      REAL(wp)            :: AKR,AHPLUS
      REAL(wp)            :: PROD,tol1,tol2,tol3,tol4,steg,fak
      REAL(wp)            :: STEGBY,Y,X,W,X1,Y1,X2,Y2,FACTOR,TERM,Z
      REAL(wp), DIMENSION(NCONC) :: CONCS
      REAL(wp), DIMENSION(NKVAL) :: AKVAL
      REAL(wp), DIMENSION(9)     :: CONCS2
      REAL(wp), DIMENSION(4)     :: AKVAL2

! AXY (07/05/13) ==================================================
! setup array to store old values of concs
      real(wp),DIMENSION(:) :: old_CONCS(NCONC)
! AXY (07/05/13) ==================================================

      EQUIVALENCE (CTOT  , CONCS2(1)), (ALK   , CONCS2(2)),  &
     &            (PCO2  , CONCS2(3)), (PH    , CONCS2(4)),  &
     &            (H2CO3 , CONCS2(5)), (HCO3  , CONCS2(6)),  &
     &            (CO3   , CONCS2(7)), (ALKC  , CONCS2(8)),  &
     &            (ALKB  , CONCS2(9)),                       &
     &            (AKP   , AKVAL2(1)), (AK1C  , AKVAL2(2)),  &
     &            (AK2C  , AKVAL2(3)), (AKB   , AKVAL2(4))
      LOGICAL             :: BORON,DONE

! AXY (07/05/13) ==================================================
!  DERIVING PH REQUIRES FOLLOWING LOOP TO CONVERGE.
!  THIS SUBROUTINE RELIES ON CONVERGENCE.  IF THE ENVIRONMENTAL
!  CONDITIONS DO NOT ALLOW FOR CONVERGENCE (IN 3D MODEL THIS IS
!  LIKELY TO OCCUR NEAR LOW SALINITY REGIONS) THE MODEL WILL 
!  BE STUCK IN THE LOOP.  TO AVOID THIS A CONVERGENCE CONDITION
!  IS PUT IN PLACE TO SET A FLAGG OF -99 IN THE PH VAR FOR NON CONVEGENCE. 
!  THE MODEL IS THEN ALLOWED TO CONTINUE. 'COUNTER, C_SW,C_CHECK' ARE 
!  THE LOCAL VARS USED.
! C_SW = condition of convergence 0=yes, 1= no
! COUNTER = number of iterations
! C_CHECK = maximum number of iterations

! SET COUNTER AND SWITCH TO ZERO AND OFF
      COUNTER=0
      C_SW=0
! FROM EXPERIENCE IF THE ITERATIONS IN THE FOLLOWING DO LOOP
! EXCEEDS 15 CONVERGENCE WILL NOT OCCUR.  THE OVERHEAD OF 25 ITERATIONS 
! IS OK FOR SMALL DOMAINS WITH 1/10 AND 1/15 DEG RESOLUTION.
! I RECOMMEND A LOWER VALUE OF 15 FOR HIGHER RESOLUTION OR LARGER DOMAINS.
      C_CHECK=25
! AXY (07/05/13) ==================================================

      DO 100 II=1,NCONC
        CONCS2(II)=CONCS(II)

! AXY (07/05/13) ==================================================
! IF CONVERGENCE IS NOT ACHIEVED THE LOCAL ARRAY CONCS MUST BE STORED TO 
! ALLOW THE MODEL TO CONTINUE. THEREFORE ....
! UPDATE OLD_CONCS
        old_CONCS(II)=CONCS(II)
! AXY (07/05/13) ==================================================

100    CONTINUE
      DO 110 II=1,NKVAL
        AKVAL2(II)=AKVAL(II)
110    CONTINUE
      AKR = AK1C/AK2C
      AHPLUS=10.0D0**(-PH)
      PROD=AKR*AKP*PCO2
     
      IF(BORON) THEN

        IF(ICALC.EQ.1.OR.ICALC.EQ.4) THEN
!         *** ALK, BTOT AND CTOT OR PCO2 FIXED ***
!         *** ITERATIVE CALCULATION NECESSARY HERE

!         SET INITIAL GUESSES AND TOLERANCE
          H2CO3=PCO2*AKP
          CO3=ALK/10.0D0
          AHPLUS=1.0D-8
          ALKB=BTOT
          TOL1=ALK/1.0D5
          TOL2=H2CO3/1.0D5
          TOL3=CTOT/1.0D5
          TOL4=BTOT/1.0D5

!         HALTAFALL iteration to determine CO3, ALKB, AHPLUS
          KARL=1
          STEG=2.0D0
          FAK=1.0D0
          STEGBY=0.4D0
10        DONE=.TRUE.

! AXY (07/05/13) ==================================================
! SET COUNTER UPDATE. FLAG 10 IS THE POINT OF RETURN FOR 
! THE CONVERGENCE CONDITION
          COUNTER=COUNTER+1
	  iters  = COUNTER
! CHECK IF CONVERGENCE HAS OCCURED IN THE NUMBER OF 
! ACCEPTABLE ITTERATIONS.  SET C_SW TO LET MODEL KNOW
! WHAT TO DO AT THE END OF THE SUBROUTINE
          if(counter.ge.c_check)then
             IF(lwp) THEN
                WRITE(numout,*) ' CO2CLC : ITERS WARNING, ', iters
             ENDIF
             c_sw=1
             GOTO 123
          endif
! AXY (07/05/13) ==================================================

          IF(ICALC.EQ.4) THEN
!         *** PCO2 IS FIXED ***
            Y=AHPLUS*AHPLUS*CO3/(AK1C*AK2C)
            IF(ABS(Y-H2CO3).GT.TOL2) THEN
              CO3=CO3*H2CO3/Y
              DONE=.FALSE.
            ENDIF
          ELSEIF(ICALC.EQ.1) THEN
!           *** CTOT IS FIXED ***
            Y=CO3*(1.0D0+AHPLUS/AK2C+AHPLUS*AHPLUS/(AK1C*AK2C))
            IF(ABS(Y-CTOT).GT.TOL3) THEN
              CO3=CO3*CTOT/Y
              DONE=.FALSE.
            ENDIF
          ENDIF
          Y=ALKB*(1.0D0+AHPLUS/AKB)
          IF(ABS(Y-BTOT).GT.TOL4) THEN
            ALKB=ALKB*BTOT/Y
            DONE=.FALSE.
          ENDIF

! Alkalinity is equivalent to -(total H+), so the sign of W is opposite
! to that normally used

          Y=CO3*(2.0D0+AHPLUS/AK2C)+ALKB
          IF(ABS(Y-ALK).GT.TOL1) THEN
            DONE=.FALSE.
            X=LOG(AHPLUS)
!           W=SIGN(1.0D0,Y-ALK)
            IF ( Y-ALK .GE. 0.0D0 ) THEN
            W=1.0D0
            ELSE
            W=-1.0D0
            ENDIF
            IF(W.GE.0.0D0) THEN
              X1=X
              Y1=Y
            ELSE
              X2=X
              Y2=Y
            ENDIF
            LQ=KARL
            IF(LQ.EQ.1) THEN
              KARL=2*NINT(W)
            ELSEIF(IABS(LQ).EQ.2.AND.(LQ*W).LT.0.) THEN
              FAK=0.5D0
              KARL=3
            ENDIF
            IF(KARL.EQ.3.AND.STEG.LT.STEGBY) THEN
              W=(X2-X1)/(Y2-Y1)
              X=X1+W*(ALK-Y1)
            ELSE
              STEG=STEG*FAK
              X=X+STEG*W
            ENDIF
            AHPLUS=EXP(X)  
          ENDIF
          IF(.NOT.DONE) GOTO 10
          
          HCO3=CO3*AHPLUS/AK2C
          IF(ICALC.EQ.4) THEN
            CTOT=H2CO3+HCO3+CO3
          ELSEIF(ICALC.EQ.1) THEN
            H2CO3=HCO3*AHPLUS/AK1C
            PCO2=H2CO3/AKP
          ENDIF
          PH=-LOG10(AHPLUS)
          ALKC=ALK-ALKB
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT, PCO2, AND BTOT FIXED ***
          Y=SQRT(PROD*(PROD-4.0D0*AKP*PCO2+4.0D0*CTOT))
          H2CO3=PCO2*AKP
          HCO3=(Y-PROD)/2.0D0
          CO3=CTOT-H2CO3-HCO3
          ALKC=HCO3+2.0D0*CO3
          AHPLUS=AK1C*H2CO3/HCO3
          PH=-LOG10(AHPLUS)
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          ALK=ALKC+ALKB
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT, PH AND BTOT FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
          CO3=FACTOR*AK1C*AK2C
          HCO3=FACTOR*AK1C*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/AKP
          ALKC=HCO3+2.0D0*CO3
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          ALK=ALKC+ALKB
        ELSEIF(ICALC.EQ.5) THEN
!         *** ALK, PH AND BTOT FIXED ***
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          ALKC=ALK-ALKB
          HCO3=ALKC/(1.0D0+2.0D0*AK2C/AHPLUS)
          CO3=HCO3*AK2C/AHPLUS
          H2CO3=HCO3*AHPLUS/AK1C
          PCO2=H2CO3/AKP
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2, PH AND BTOT FIXED ***
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          H2CO3=PCO2*AKP
          HCO3=H2CO3*AK1C/AHPLUS
          CO3=HCO3*AK2C/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          ALKC=HCO3+2.0D0*CO3
          ALK=ALKC+ALKB
        ENDIF
      ELSE
        IF(ICALC.EQ.1) THEN
!         *** CTOT AND ALK FIXED ***
          TERM=4.0D0*ALK+CTOT*AKR-ALK*AKR
          Z=SQRT(TERM*TERM+4.0D0*(AKR-4.0D0)*ALK*ALK)
          CO3=(ALK*AKR-CTOT*AKR-4.0D0*ALK+Z)/(2.0D0*(AKR-4.0D0))
          HCO3=(CTOT*AKR-Z)/(AKR-4.0D0)
          H2CO3=CTOT-ALK+CO3
          PCO2=H2CO3/AKP
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT AND PCO2 FIXED ***
          Y=SQRT(PROD*(PROD-4.0D0*AKP*PCO2+4.0D0*CTOT))
          H2CO3=PCO2*AKP
          HCO3=(Y-PROD)/2.0D0
          CO3=CTOT-H2CO3-HCO3
          ALK=HCO3+2.0D0*CO3
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT AND PH FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
          CO3=FACTOR*AK1C*AK2C
          HCO3=FACTOR*AK1C*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/AKP
          ALK=HCO3+2.0D0*CO3
        ELSEIF(ICALC.EQ.4) THEN
!         *** ALK AND PCO2 FIXED ***
          TERM=SQRT((8.0D0*ALK+PROD)*PROD)
          CO3=ALK/2.0D0+PROD/8.0D0-TERM/8.0D0
          HCO3=-PROD/4.0D0+TERM/4.0D0
          H2CO3=PCO2*AKP
          CTOT=CO3+HCO3+H2CO3
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.5) THEN
!         *** ALK AND PH FIXED ***
          HCO3=ALK/(1.0D0+2.0D0*AK2C/AHPLUS)
          CO3=HCO3*AK2C/AHPLUS
          H2CO3=HCO3*AHPLUS/AK1C
          PCO2=H2CO3/AKP
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2 AND PH FIXED ***
          H2CO3=PCO2*AKP
          HCO3=H2CO3*AK1C/AHPLUS
          CO3=HCO3*AK2C/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          ALK=HCO3+2.0D0*CO3
        ENDIF
      ENDIF

      DO 120 II=1,NCONC
        CONCS(II)=CONCS2(II)
120    CONTINUE

! AXY (07/05/13) ==================================================
! C_SW IS SET AT 0 TO START
! THEN IF NON CONVERGENCE C_SW SET TO 1
123   IF(C_SW.EQ.1)THEN
! IF NON CONVERGENCE, THE MODEL REQUIRES CONCS TO CONTAIN USABLE VALUES.
! BEST OFFER BEING THE OLD CONCS VALUES WHEN CONVERGENCE HAS BEEN 
! ACHIEVED
        DO II=1,NCONC
           CONCS(II)=OLD_CONCS(II)
        END DO

! SPECIFIC CARBONATE VALUES TO PUSH CODE ON THROUGH THE
! NON CONVERGENCE CONDITIONS
! PCO2W = 0 SO CO2 UPTAKE WILL BE ENCOURAGED
!       CONCS(3)=O3C(III)*0.005_fp8/1.e6_fp8
! -99 IS A FLAG TO SHOW WHEN AND WHERE NON CONVERGENCE OCCURED
! CONCS(4)=PH, OUTPUT OF -99 MEANS NON-CONVERGENCE
!       CONCS(4)=-99._fp8
! 
! AXY (10/05/13): remove this -99 flag; this is handled instead up in
!                 trcbio_medusa.F90 where the iters variable is both 
!                 output and may be used to trigger action
!
! RESET SWITCH FOR NEXT CALL TO THIS SUBROUTINE
        C_SW=0
      ENDIF
! AXY (07/05/13) ==================================================


      RETURN

      END SUBROUTINE
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!
      SUBROUTINE CaCO3_Saturation (Tc, S, D, CO3, &
      Om_cal, Om_arg)
!      
!=======================================================================
!
! Routine to calculate the saturation state of calcite and aragonite    
! Inputs:                                                               
!               Tc      Temperature (C)                                 
!               S       Salinity                                        
!               D       Depth (m)                                       
!               CO3     Carbonate ion concentration (mol.kg-1 ie /1D6)
!                                                                       
! Outputs                                                               
!               Om_cal  Calite saturation                               
!               Om_arg  Aragonite saturation                            
!                                                                       
! Intermediates                                                         
!               K_cal   Stoichiometric solubility product for calcite   
!               K_arg   Stoichiometric solubility product for aragonite 
!               Ca      Calcium 2+ concentration (mol.kg-1) 
!               P       Pressure (bars)
!                                                                       
! Source                                                                
!       Zeebe & Wolf-Gladrow 2001 following Mucci (1983)                
!       with pressure corrections from Millero (1995)                   
!       Code tested against reference values given in Z & W-G
!       Built Jerry Blackford, 2008
!=======================================================================

        IMPLICIT None
        REAL(wp), INTENT( in )    :: Tc, S, D, CO3
        REAL(wp), INTENT( inout ) :: Om_cal, Om_arg
        REAL(wp) :: Tk, Kelvin, Ca
        REAL(wp) :: logKspc, Kspc
        REAL(wp) :: logKspa, Kspa
        REAL(wp) :: tmp1, tmp2, tmp3
        REAL(wp) :: dV, dK, P, R
               
! setup
        Kelvin = 273.15
        Tk = Tc + Kelvin 
        Ca = 0.01028    ! Currently oceanic mean value at S=25, needs refining)
        Ca = 0.010279 * (S / 35.0)  ! Ca varies with salinity (cf. Feeley et al., 2004; Yool et al., 2010)
        R = 83.131      !(cm3.bar.mol-1.K-1)
        P = D / 10.0    !pressure in bars

! calculate K for calcite
        tmp1 = -171.9065 - (0.077993*Tk) + (2839.319/Tk) + 71.595*log10(Tk) 
        tmp2 = + (-0.77712 + (0.0028426*Tk) + (178.34/Tk))*SQRT(S) 
        tmp3 = - (0.07711*S) + (0.0041249*(S**1.5))
        logKspc = tmp1 + tmp2 + tmp3
        Kspc = 10.0**logKspc
      
! correction for pressure for calcite
        IF ( D .GT. 0) THEN
          dV = -48.76 + 0.5304*Tc
          dK = -11.76/1.0D3 + (0.3692/1.0D3) * Tc
          tmp1 = -(dV/(R*Tk))*P + (0.5*dK/(R*Tk))*P*P
          Kspc = Kspc*exp(tmp1)
          logKspc = log10(Kspc)
        END IF
        
! calculate K for aragonite
        tmp1 = -171.945 - 0.077993*Tk + 2903.293 / Tk + 71.595* log10(Tk)
        tmp2 = + (-0.068393 + 0.0017276*Tk + 88.135/Tk)*SQRT(S)    
        tmp3 = - 0.10018*S + 0.0059415*S**1.5
        logKspa = tmp1 + tmp2 + tmp3
        Kspa = 10.0**logKspa
     
! correction for pressure for aragonite
        IF ( D .GT. 0) THEN
          dV = -46.00 + 0.5304*Tc
          dK = -11.76/1.0D3 + (0.3692/1.0D3) * Tc
          tmp1 = -(dV/(R*Tk))*P + (0.5*dK/(R*Tk))*P*P
          Kspa = Kspa*exp(tmp1)
          logKspa = log10(Kspa)
        END IF
        
! calculate saturation states
        Om_cal = (CO3 * Ca) / Kspc
        Om_arg = (CO3 * Ca) / Kspa
      
      RETURN

      END SUBROUTINE
!-----------------------------------------------------------------------

!=======================================================================
!=======================================================================
!=======================================================================

#  else
   !!======================================================================
   !!  Dummy module :                         No MEDUSA carbonate chemistry
   !!======================================================================

CONTAINS

   SUBROUTINE trc_co2_medusa( kt )                   ! Empty routine

      INTEGER, INTENT( in ) ::   kt

      WRITE(*,*) 'trc_co2_medusa: You should not have seen this print! error?', kt

   END SUBROUTINE trc_co2_medusa
#endif 

   !!======================================================================
END MODULE  trcco2_medusa
