MODULE mocsy_mainmod

      USE in_out_manager  ! I/O manager

CONTAINS

! ----------------------------------------------------------------------
!  SW_ADTG
! ----------------------------------------------------------------------
!
!> \file sw_adtg.f90
!! \BRIEF 
!> Module with sw_adtg function - compute adiabatic temp. gradient from S,T,P
!>  Function to calculate adiabatic temperature gradient as per UNESCO 1983 routines.
FUNCTION sw_adtg  (s,t,p)

  !     ==================================================================
  !     Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================
  USE mocsy_singledouble
  IMPLICIT NONE
  !> salinity [psu (PSU-78)]
  REAL(kind=wp) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=wp) :: t
  !> pressure [db]
  REAL(kind=wp) :: p

  REAL(kind=wp) :: a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2
  REAL(kind=wp) :: sref

  REAL(kind=wp) :: sw_adtg

  sref = 35.d0
  a0 =  3.5803d-5
  a1 = +8.5258d-6
  a2 = -6.836d-8
  a3 =  6.6228d-10

  b0 = +1.8932d-6
  b1 = -4.2393d-8

  c0 = +1.8741d-8
  c1 = -6.7795d-10
  c2 = +8.733d-12
  c3 = -5.4481d-14

  d0 = -1.1351d-10
  d1 =  2.7759d-12

  e0 = -4.6206d-13
  e1 = +1.8676d-14
  e2 = -2.1687d-16

  sw_adtg =  a0 + (a1 + (a2 + a3*T)*T)*T &
       + (b0 + b1*T)*(S-sref) &
       + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P &
       + (  e0 + (e1 + e2*T)*T )*P*P

END FUNCTION sw_adtg

! ----------------------------------------------------------------------
!  SW_ADTG
! ----------------------------------------------------------------------
!
!> \file sw_ptmp.f90
!! \BRIEF 
!> Module with sw_ptmp function - compute potential T from in-situ T
!> Function to calculate potential temperature [C] from in-situ temperature
FUNCTION sw_ptmp  (s,t,p,pr)

  !     ==================================================================
  !     Calculates potential temperature [C] from in-situ Temperature [C]
  !     From UNESCO 1983 report.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = temperature         [degree C (IPTS-68)]
  !     p  = pressure            [db]
  !     pr = reference pressure  [db]

  USE mocsy_singledouble
  IMPLICIT NONE

! Input arguments
  !> salinity [psu (PSS-78)]
  REAL(kind=wp) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=wp) :: t
  !> pressure [db]
  REAL(kind=wp) :: p
  !> reference pressure  [db]  
  REAL(kind=wp) :: pr

! local arguments
  REAL(kind=wp) :: del_P ,del_th, th, q
  REAL(kind=wp) :: onehalf, two, three
  PARAMETER (onehalf = 0.5d0, two = 2.d0, three = 3.d0 )

! REAL(kind=wp) :: sw_adtg
! EXTERNAL sw_adtg

! Output 
  REAL(kind=wp) :: sw_ptmp

  ! theta1
  del_P  = PR - P
  del_th = del_P*sw_adtg(S,T,P)
  th     = T + onehalf*del_th
  q      = del_th

  ! theta2
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 - 1.d0/SQRT(two))*(del_th - q)
  q      = (two-SQRT(two))*del_th + (-two+three/SQRT(two))*q

  ! theta3
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 + 1.d0/SQRT(two))*(del_th - q)
  q      = (two + SQRT(two))*del_th + (-two-three/SQRT(two))*q

  ! theta4
  del_th = del_P*sw_adtg(S,th,P+del_P)
  sw_ptmp     = th + (del_th - two*q)/(two*three)

  RETURN
END FUNCTION sw_ptmp


! ----------------------------------------------------------------------
!  SW_TEMP
! ----------------------------------------------------------------------
!
!> \file sw_temp.f90
!! \BRIEF 
!> Module with sw_temp function - compute in-situ T from potential T
!> Function to compute in-situ temperature [C] from potential temperature [C]
FUNCTION sw_temp( s, t, p, pr )
  !     =============================================================
  !     SW_TEMP
  !     Computes in-situ temperature [C] from potential temperature [C]
  !     Routine available in seawater.f (used for MIT GCM)
  !     Downloaded seawater.f (on 17 April 2009) from
  !     http://ecco2.jpl.nasa.gov/data1/beaufort/MITgcm/bin/
  !     =============================================================

  !     REFERENCES:
  !     Fofonoff, P. and Millard, R.C. Jr
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !     Eqn.(31) p.39

  !     Bryden, H. 1973.
  !     "New Polynomials for thermal expansion, adiabatic temperature gradient
  !     and potential temperature of sea water."
  !     DEEP-SEA RES., 1973, Vol20,401-408.
  !     =============================================================

  !     Simple modifications: J. C. Orr, 16 April 2009
  !     - combined fortran code from MITgcm site & simplification in
  !       CSIRO code (matlab equivalent) from Phil Morgan

  USE mocsy_singledouble
  IMPLICIT NONE

  !     Input arguments:
  !     -----------------------------------------------
  !     s  = salinity              [psu      (PSS-78) ]
  !     t  = potential temperature [degree C (IPTS-68)]
  !     p  = pressure              [db]
  !     pr = reference pressure    [db]

  !> salinity [psu (PSS-78)]
  REAL(kind=wp) ::   s
  !> potential temperature [degree C (IPTS-68)]
  REAL(kind=wp) ::   t
  !> pressure [db]
  REAL(kind=wp) ::   p
  !> reference pressure [db]
  REAL(kind=wp) ::   pr

  REAL(kind=wp) ::  ds, dt, dp, dpr
  REAL(kind=wp) :: dsw_temp

  REAL(kind=wp) ::   sw_temp
! EXTERNAL sw_ptmp
! REAL(kind=wp) ::   sw_ptmp

  ds = DBLE(s)
  dt = DBLE(t)
  dp = DBLE(p)
  dpr = DBLE(pr)

  !    Simple solution
  !    (see https://svn.mpl.ird.fr/us191/oceano/tags/V0/lib/matlab/seawater/sw_temp.m)
  !    Carry out inverse calculation by swapping P_ref (pr) and Pressure (p)
  !    in routine that is normally used to compute potential temp from temp
  dsw_temp = sw_ptmp(ds, dt, dpr, dp)
  sw_temp = REAL(dsw_temp)

  !    The above simplification works extremely well (compared to Table in 1983 report)
  !    whereas the sw_temp routine from MIT GCM site does not seem to work right

  RETURN
END FUNCTION sw_temp

! ----------------------------------------------------------------------
!  TPOT
! ----------------------------------------------------------------------
!
!> \file tpot.f90
!! \BRIEF 
!>    Module with tpot subroutine - compute potential T from in situ T,S,P
!>    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
!!    This subroutine is needed because sw_ptmp is a function (using scalars not arrays)
SUBROUTINE tpot(salt, tempis, press, pressref, N, tempot)
  !    Purpose:
  !    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_ptmp is a function

  USE mocsy_singledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> salinity [psu]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: salt
  !> in situ temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: tempis
  !> pressure [db]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: press
!f2py optional , depend(salt) :: n=len(salt)
  !> pressure reference level [db]
  REAL(kind=wp), INTENT(in) :: pressref

! OUTPUT variables:
  !> potential temperature [C] for pressref
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: tempot

  REAL(kind=wp) :: dsalt, dtempis, dpress, dpressref
  REAL(kind=wp) :: dtempot

  INTEGER :: i

! REAL(kind=wp) :: sw_ptmp
! EXTERNAL sw_ptmp

  DO i = 1,N
     dsalt     = DBLE(salt(i))
     dtempis   = DBLE(tempis(i))
     dpress    = DBLE(press(i))
     dpressref = DBLE(pressref)

     dtempot   = sw_ptmp(dsalt, dtempis, dpress, dpressref)

     tempot(i) = REAL(dtempot)
  END DO

  RETURN
END SUBROUTINE tpot

! ----------------------------------------------------------------------
!  TIS
! ----------------------------------------------------------------------
!
!> \file tis.f90
!! \BRIEF 
!>    Module with tis subroutine - compute in situ T from S,T,P
!>    Compute in situ temperature from arrays of potential temp, salinity, and pressure.
!!    This subroutine is needed because sw_temp is a function (using scalars not arrays)
SUBROUTINE tis(salt, tempot, press, pressref, N, tempis)
  !    Purpose:
  !    Compute in situ temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_temp is a function

  USE mocsy_singledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> salinity [psu]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: salt
  !> potential temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: tempot
  !> pressure [db]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: press
!f2py optional , depend(salt) :: n=len(salt)
  !> pressure reference level [db]
  REAL(kind=wp), INTENT(in) :: pressref

! OUTPUT variables:
  !> in situ temperature [C] 
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: tempis

! REAL(kind=wp) :: dsalt, dtempis, dpress, dpressref
! REAL(kind=wp) :: dtempot

  INTEGER :: i

! REAL(kind=wp) :: sw_temp
! REAL(kind=wp) :: sw_temp
! EXTERNAL sw_temp

  DO i = 1,N
    !dsalt     = DBLE(salt(i))
    !dtempot   = DBLE(tempot(i))
    !dpress    = DBLE(press(i))
    !dpressref = DBLE(pressref)
    !dtempis   = sw_temp(dsalt, dtempot, dpress, dpressref)
    !tempis(i) = REAL(dtempis)

     tempis   = sw_temp(salt(i), tempot(i), press(i), pressref)
  END DO

  RETURN
END SUBROUTINE tis

! ----------------------------------------------------------------------
!  P80
! ----------------------------------------------------------------------
!
!> \file p80.f90
!! \BRIEF 
!> Module with p80 function - compute pressure from depth
!>     Function to compute pressure from depth using Saunder's (1981) formula with eos80.
FUNCTION p80(dpth,xlat)

  !     Compute Pressure from depth using Saunder's (1981) formula with eos80.

  !     Reference:
  !     Saunders, Peter M. (1981) Practical conversion of pressure
  !     to depth, J. Phys. Ooceanogr., 11, 573-574, (1981)

  !     Coded by:
  !     R. Millard
  !     March 9, 1983
  !     check value: p80=7500.004 dbars at lat=30 deg., depth=7321.45 meters

  !     Modified (slight format changes + added ref. details):
  !     J. Orr, 16 April 2009

  USE mocsy_singledouble
  IMPLICIT NONE

! Input variables:
  !> depth [m]
  REAL(kind=wp), INTENT(in) :: dpth
  !> latitude [degrees]
  REAL(kind=wp), INTENT(in) :: xlat

! Output variable:
  !> pressure [db]
  REAL(kind=wp) :: p80

! Local variables:
  REAL(kind=wp) :: pi
  REAL(kind=wp) :: plat, d, c1

  pi=3.141592654

  plat = ABS(xlat*pi/180.)
  d  = SIN(plat)
  c1 = 5.92e-3+d**2 * 5.25e-3

  p80 = ((1-c1)-SQRT(((1-c1)**2)-(8.84e-6*dpth))) / 4.42e-6

  RETURN
END FUNCTION p80

! ----------------------------------------------------------------------
!  RHO
! ----------------------------------------------------------------------
!
!> \file rho.f90
!! \BRIEF 
!> Module with rho function - computes in situ density from S, T, P
!> Function to compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)
FUNCTION rho(salt, temp, pbar)

  ! Compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)

  USE mocsy_singledouble
  IMPLICIT NONE

  !> salinity [psu]
  REAL(kind=wp) :: salt
  !> in situ temperature (C)
  REAL(kind=wp) :: temp
  !> pressure (bar) [Note units: this is NOT in decibars]
  REAL(kind=wp) :: pbar

  REAL(kind=wp) :: s, t, p
! REAL(kind=wp) :: t68
  REAL(kind=wp) :: X
  REAL(kind=wp) :: rhow, rho0
  REAL(kind=wp) :: a, b, c
  REAL(kind=wp) :: Ksbmw, Ksbm0, Ksbm
  REAL(kind=wp) :: drho

  REAL(kind=wp) :: rho

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = in situ temperature [degree C (IPTS-68)]
  !     p  = pressure            [bar] !!!!  (not in [db]

  s = DBLE(salt)
  t = DBLE(temp)
  p = DBLE(pbar)

! Convert the temperature on today's "ITS 90" scale to older IPTS 68 scale
! (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
! According to Best-Practices guide, line above should be commented & 2 lines below should be uncommented
! Guide's answer of rho (s=35, t=25, p=0) = 1023.343 is for temperature given on ITPS-68 scale
! t68 = (T - 0.0002) / 0.99975
! X = t68
! Finally, don't do the ITS-90 to IPTS-68 conversion (T input var now already on IPTS-68 scale)
  X = T

! Density of pure water
  rhow = 999.842594d0 + 6.793952e-2_wp*X          &
       -9.095290e-3_wp*X*X + 1.001685e-4_wp*X**3  &
       -1.120083e-6_wp*X**4 + 6.536332e-9_wp*X**5

! Density of seawater at 1 atm, P=0
  A = 8.24493e-1_wp - 4.0899e-3_wp*X                         &
       + 7.6438e-5_wp*X*X - 8.2467e-7_wp*X**3 + 5.3875e-9_wp*X**4
  B = -5.72466e-3_wp + 1.0227e-4_wp*X - 1.6546e-6_wp*X*X
  C = 4.8314e-4_wp

  rho0 = rhow + A*S + B*S*SQRT(S) + C*S**2.0d0

! Secant bulk modulus of pure water
! The secant bulk modulus is the average change in pressure
! divided by the total change in volume per unit of initial volume.
  Ksbmw = 19652.21d0 + 148.4206d0*X - 2.327105d0*X*X &
       + 1.360477e-2_wp*X**3 - 5.155288e-5_wp*X**4

! Secant bulk modulus of seawater at 1 atm
  Ksbm0 = Ksbmw + S*( 54.6746d0 - 0.603459d0*X + 1.09987e-2_wp*X**2 &
       - 6.1670e-5_wp*X**3) &
       + S*SQRT(S)*( 7.944e-2_wp + 1.6483e-2_wp*X - 5.3009e-4_wp*X**2)

! Secant bulk modulus of seawater at S,T,P
  Ksbm = Ksbm0 &
       + P*(3.239908d0 + 1.43713e-3_wp*X + 1.16092e-4_wp*X**2 - 5.77905e-7_wp*X**3) &
       + P*S*(2.2838e-3_wp - 1.0981e-5_wp*X - 1.6078e-6_wp*X**2) &
       + P*S*SQRT(S)*1.91075e-4_wp &
       + P*P*(8.50935e-5_wp - 6.12293e-6_wp*X + 5.2787e-8_wp*X**2) &
       + P*P*S*(-9.9348e-7_wp + 2.0816e-8_wp*X + 9.1697e-10_wp*X**2)

! Density of seawater at S,T,P
  drho = rho0/(1.0d0 - P/Ksbm)
  rho = REAL(drho)

  RETURN
END FUNCTION rho

! ----------------------------------------------------------------------
!  RHOINSITU
! ----------------------------------------------------------------------
!
!> \file rhoinsitu.f90
!! \BRIEF 
!> Module with rhoinsitu subroutine - compute in situ density from S, Tis, P
!>     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db).
!!     This subroutine is needed because rho is a function (using scalars not arrays)
SUBROUTINE rhoinsitu(salt, tempis, pdbar, N, rhois)

  !     Purpose:
  !     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db)
  !     Needed because rho is a function

  USE mocsy_singledouble
  IMPLICIT NONE

  INTEGER :: N

! INPUT variables
  ! salt   = salinity [psu]
  ! tempis = in situ temperature [C]
  ! pdbar  = pressure [db]

  !> salinity [psu]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: salt
  !> in situ temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: tempis
  !> pressure [db]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: pdbar
!f2py optional , depend(salt) :: n=len(salt)

! OUTPUT variables:
  ! rhois  = in situ density

  !> in situ density [kg/m3]
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: rhois

! Local variables
  INTEGER :: i

! REAL(kind=wp) ::  rho
! EXTERNAL rho

  DO i = 1,N
     rhois(i) = rho(salt(i), tempis(i), pdbar(i)/10.)
  END DO

  RETURN
END SUBROUTINE rhoinsitu

! ----------------------------------------------------------------------
!  DEPTH2PRESS
! ----------------------------------------------------------------------
!
!> \file depth2press.f90
!! \BRIEF 
!> Module with depth2press subroutine - converts depth to pressure
!! with Saunders (1981) formula
!>     Compute pressure [db] from depth [m] & latitude [degrees north].
!!     This subroutine is needed because p80 is a function (using scalars not arrays)
SUBROUTINE depth2press(depth, lat, pdbar, N)

  !     Purpose:
  !     Compute pressure [db] from depth [m] & latitude [degrees north].
  !     Needed because p80 is a function 

  USE mocsy_singledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> depth [m]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: depth
  !> latitude [degrees]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: lat
!f2py optional , depend(depth) :: n=len(depth)

! OUTPUT variables:
  !> pressure [db]
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: pdbar

  !     Local variables
  INTEGER :: i

! REAL(kind=wp) ::  p80
! EXTERNAL p80

  DO i = 1,N
     pdbar(i) = p80(depth(i), lat(i))
  END DO

  RETURN
END SUBROUTINE depth2press

! ----------------------------------------------------------------------
!  CONSTANTS
! ----------------------------------------------------------------------
!
!> \file constants.f90
!! \BRIEF 
!> Module with contants subroutine - computes carbonate system constants
!! from S,T,P 
!> Compute thermodynamic constants
!! FROM temperature, salinity, and pressure (1D arrays)
SUBROUTINE constants(K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,  &
                     K1p, K2p, K3p, Ksi,                      &
                     St, Ft, Bt,                              &
                     temp, sal, Patm,                         &
                     depth, lat, N,                           &
                     optT, optP, optB, optK1K2, optKf, optGAS)

  !   Purpose:
  !     Compute thermodynamic constants
  !     FROM: temperature, salinity, and pressure (1D arrays)

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
  !     ---------
  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (p = depth)
  !     ---------
  !     optB:
  !     ---------
  !       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
  !       -> 'l10' means use newer   formulation of Lee et al. (2010) for total Boron
  !     ---------
  !     optK1K2:
  !     ---------
  !       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
  !                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
  !       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
  !       -> 'w14' means use Waters (2014) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
  !     -----------
  !     optKf:
  !     ----------
  !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
  !               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
  !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2
  !     ---------
  !       PRESSURE corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  !     St, Ft, Bt

  USE mocsy_singledouble
  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> in <b>situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) [degree C]
  REAL(kind=wp), INTENT(in),    DIMENSION(N) :: temp
  !> depth in <b>meters</b> (when optP='m') or <b>decibars</b> (when optP='db')
  REAL(kind=wp), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=wp), INTENT(in),    DIMENSION(N) :: lat
  !> salinity <b>[psu]</b>
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: sal
!f2py optional , depend(sal) :: n=len(sal)

  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: Patm

  !> for temp input, choose \b 'Tinsitu' for in situ Temp or 
  !! \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
  CHARACTER(7), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(2), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20),
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
!f2py character*3 optional, intent(in) :: optB='l10'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) or \b 'w14' (Waters et al., 2014)
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Ouput variables
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: K2
  !> equilibrium constant for dissociation of boric acid 
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Ksi
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Ft
  !> total boron
  !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: Bt

! Local variables
  REAL(kind=wp) :: ssal
  REAL(kind=wp) :: p
  REAL(kind=wp) :: tempot, tempis68, tempot68
  REAL(kind=wp) :: tempis
  REAL(kind=wp) :: is, invtk, dlogtk, is2, s2, sqrtis
  REAL(kind=wp) :: Ks_0p, Kf_0p
  REAL(kind=wp) :: total2free, free2SWS, total2SWS, SWS2total
  REAL(kind=wp) :: total2free_0p, free2SWS_0p, total2SWS_0p
! REAL(kind=wp) :: free2SWS, free2SWS_0p

  REAL(kind=wp) :: dtempot, dtempot68
  REAL(kind=wp) :: R

  REAL(kind=wp) :: pK1o, ma1, mb1, mc1, pK1
  REAL(kind=wp) :: pK2o, ma2, mb2, mc2, pK2

  REAL(kind=wp), DIMENSION(12) :: a0, a1, a2, b0, b1, b2
  REAL(kind=wp), DIMENSION(12) :: deltav, deltak, lnkpok0
  REAL(kind=wp) :: tmp, nK0we74

  INTEGER :: i, icount, ipc

  REAL(kind=wp) :: t, tk, tk0, prb
  REAL(kind=wp) :: s, sqrts, s15, scl

  REAL(kind=wp) :: Phydro_atm, Patmd, Ptot, Rgas_atm, vbarCO2

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

  ! CONSTANTS
  ! =========
  ! Constants in formulation for Pressure effect on K's (Millero, 95)
  ! with corrected coefficients for Kb, Kw, Ksi, etc.

  ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
  !            9) K1P, 10) K2P, 11) K3P, 12) Ksi

  DATA a0 /-25.5_wp, -15.82_wp, -29.48_wp, -20.02_wp, &
          -18.03_wp,  -9.78_wp, -48.76_wp, -45.96_wp, &
          -14.51_wp, -23.12_wp, -26.57_wp, -29.48_wp/
  DATA a1 /0.1271_wp, -0.0219_wp, 0.1622_wp, 0.1119_wp, &
           0.0466_wp, -0.0090_wp, 0.5304_wp, 0.5304_wp, &
           0.1211_wp, 0.1758_wp, 0.2020_wp, 0.1622_wp/
  DATA a2 /     0.0_wp,       0.0_wp, -2.608e-3_wp, -1.409e-3_wp, &
           0.316e-3_wp, -0.942e-3_wp,  0.0_wp,       0.0_wp, &
          -0.321e-3_wp, -2.647e-3_wp, -3.042e-3_wp, -2.6080e-3_wp/
  DATA b0 /-3.08e-3_wp, 1.13e-3_wp,  -2.84e-3_wp,   -5.13e-3_wp, &
           -4.53e-3_wp, -3.91e-3_wp, -11.76e-3_wp, -11.76e-3_wp, &
           -2.67e-3_wp, -5.15e-3_wp,  -4.08e-3_wp,  -2.84e-3_wp/
  DATA b1 /0.0877e-3_wp, -0.1475e-3_wp, 0.0_wp,       0.0794e-3_wp, &
           0.09e-3_wp,    0.054e-3_wp,  0.3692e-3_wp, 0.3692e-3_wp, &
           0.0427e-3_wp,  0.09e-3_wp,   0.0714e-3_wp, 0.0_wp/
  DATA b2 /12*0.0_wp/

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
    opB = optB
  ELSE
    opB = 'l10'
  ENDIF
  IF (PRESENT(optKf)) THEN
    opKf = optKf
  ELSE
    opKf = 'pf'
  ENDIF
  IF (PRESENT(optK1K2)) THEN
    opK1K2 = optK1K2
  ELSE
    opK1K2 = 'l'
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

  R = 83.14472_wp

  icount = 0
  DO i = 1, N
     icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Model temperature tracer is usually "potential temperature"
!    * Model vertical grid is usually in meters
!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry, if appropriate,
!    convert these 2 model vars (input to this routine)
!     - depth [m] => convert to pressure [db]
!     - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db' ) THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p = depth(i)
     ELSE
        PRINT *,"optP must be 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = temp(i)
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002) / 0.99975
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p, 0.0_wp )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis = 0.99975*tempis68 + 0.0002
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis    = temp(i)
        tempis68  = (temp(i) - 0.0002) / 0.99975
!       dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0_wp)
        dtempot68 = sw_ptmp(sal(i), tempis68, p, 0.0_wp)
        dtempot   = 0.99975*dtempot68 + 0.0002
     ELSE
        PRINT *,"optT must be either 'Tpot' or 'Tinsitu'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF

!    Compute constants:
     IF (temp(i) >= -5. .AND. temp(i) < 1.0e+2) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (      sal(i) < 0.  .OR.  sal(i) > 1e+3) THEN
           PRINT *, 'i, icount, temp, sal =', i, icount, temp(i), sal(i)
        ENDIF
!       Zero out negative salinity (prev case for OCMIP2 model w/ slightly negative S in some coastal cells)
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF

!       Absolute temperature (Kelvin) and related values
        t = DBLE(tempis)
        tk = 273.15d0 + t
        invtk=1.0d0/tk
        dlogtk=LOG(tk)

!       Atmospheric pressure
        Patmd = DBLE(Patm(i))

!       Hydrostatic pressure (prb is in bars)
        prb = DBLE(p) / 10.0d0

!       Salinity and simply related values
        s = DBLE(ssal)
        s2=s*s
        sqrts=SQRT(s)
        s15=s**1.5d0
        scl=s/1.80655d0

!       Ionic strength:
        is = 19.924d0*s/(1000.0d0 - 1.005d0*s)
        is2 = is*is
        sqrtis = SQRT(is)

!       Total concentrations for sulfate, fluoride, and boron

!       Sulfate: Morris & Riley (1966)
        St(i) = 0.14d0 * scl/96.062d0

!       Fluoride:  Riley (1965)
        Ft(i) = 0.000067d0 * scl/18.9984d0

!       Boron:
        IF (trim(opB) == 'l10') THEN
!          New formulation from Lee et al (2010)
           Bt(i) = 0.0002414d0 * scl/10.811d0
        ELSEIF (trim(opB) == 'u74') THEN
!          Classic formulation from Uppström (1974)
           Bt(i) = 0.000232d0  * scl/10.811d0
        ELSE
           PRINT *,"optB must be 'l10' or 'u74'"
           STOP
        ENDIF

!       K0 (K Henry)
!       CO2(g) <-> CO2(aq.)
!       K0  = [CO2]/ fCO2
!       Weiss (1974)   [mol/kg/atm]
        IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
           tk0 = tk                   !in situ temperature (K) for K0 calculation
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           tk0 = dtempot + 273.15d0   !potential temperature (K) for K0 calculation as needed for potential fCO2 & pCO2
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           tk0 = tk                     !in situ temperature (K) for K0 calculation
           Phydro_atm = prb / 1.01325d0 !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
           Ptot = Patmd + Phydro_atm    !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF
        tmp = 9345.17d0/tk0 - 60.2409d0 + 23.3585d0 * LOG(tk0/100.0d0)
        nK0we74 = tmp + s*(0.023517d0 - 0.00023656d0*tk0 + 0.0047036e-4_wp*tk0*tk0)
        K0(i) = EXP(nK0we74)

!       K1 = [H][HCO3]/[H2CO3]
!       K2 = [H][CO3]/[HCO3]
        IF (trim(opK1K2) == 'l') THEN
!         Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total pH scale)
          K1(i) = 10.0d0**(-1.0d0*(3633.86d0*invtk - 61.2172d0 + 9.6777d0*dlogtk  &
                  - 0.011555d0*s + 0.0001152d0*s2))
          K2(i) = 10.0d0**(-1*(471.78d0*invtk + 25.9290d0 - 3.16967d0*dlogtk      &
                  - 0.01781d0*s + 0.0001122d0*s2))
        ELSEIF (trim(opK1K2) == 'm10') THEN
!         Millero (2010, Mar. Fresh Wat. Res.) (seawater pH scale)
          pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
          ma1 = 13.4038d0*sqrts + 0.03206d0*s - (5.242e-5)*s2
          mb1 = -530.659d0*sqrts - 5.8210d0*s
          mc1 = -2.0664d0*sqrts
          pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
          K1(i) = 10.0d0**(-pK1) 

          pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
          ma2 = 21.3728d0*sqrts + 0.1218d0*s - (3.688e-4)*s2
          mb2 = -788.289d0*sqrts - 19.189d0*s
          mc2 = -3.374d0*sqrts
          pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
          K2(i) = 10.0d0**(-pK2)
        ELSEIF (trim(opK1K2) == 'w14') THEN
!         Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014) (seawater scale)
          pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
          ma1 = 13.409160d0*sqrts + 0.031646d0*s - (5.1895e-5)*s2
          mb1 = -531.3642d0*sqrts - 5.713d0*s
          mc1 = -2.0669166d0*sqrts
          pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
          K1(i) = 10.0d0**(-pK1) 

          pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
          ma2 = 21.225890d0*sqrts + 0.12450870d0*s - (3.7243e-4_r8)*s2
          mb2 = -779.3444d0*sqrts - 19.91739d0*s
          mc2 = -3.3534679d0*sqrts
          pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
          K2(i) = 10.0d0**(-pK2)
        ELSE
           PRINT *, "optK1K2 must be either 'l' or 'm10', or 'w14'"
           STOP
        ENDIF

!       Kb = [H][BO2]/[HBO2]
!       (total scale)
!       Millero p.669 (1995) using data from Dickson (1990)
        Kb(i) = EXP((-8966.90d0 - 2890.53d0*sqrts - 77.942d0*s +  &
                1.728d0*s15 - 0.0996d0*s2)*invtk +              &
                (148.0248d0 + 137.1942d0*sqrts + 1.62142d0*s) +   &
                (-24.4344d0 - 25.085d0*sqrts - 0.2474d0*s) *      &
                dlogtk + 0.053105d0*sqrts*tk)

!       K1p = [H][H2PO4]/[H3PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 65
!       Use Millero equation's 115.540 constant instead of 115.525 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K1p(i) = EXP(-4576.752d0*invtk + 115.540d0 - 18.453d0*dlogtk +  &
                 (-106.736d0*invtk + 0.69171d0) * sqrts +             &
                 (-0.65643d0*invtk - 0.01844d0) * s)

!       K2p = [H][HPO4]/[H2PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!       Millero (1995), p.670, eq. 66
!       Use Millero equation's 172.1033 constant instead of 172.0833 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K2p(i) = EXP(-8814.715d0*invtk + 172.1033d0 - 27.927d0*dlogtk +  &
                 (-160.340d0*invtk + 1.3566d0)*sqrts +                 &
                 (0.37335d0*invtk - 0.05778d0)*s)

!       K3p = [H][PO4]/[HPO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 67
!       Use Millero equation's 18.126 constant instead of 18.141 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K3p(i) = EXP(-3070.75d0*invtk - 18.126d0 +            &
                 (17.27039d0*invtk + 2.81197d0) *             &
                 sqrts + (-44.99486d0*invtk - 0.09984d0) * s)

!       Ksi = [H][SiO(OH)3]/[Si(OH)4]
!       (seawater scale)
!       Millero (1995), p.671, eq. 72
!       Use Millero equation's 117.400 constant instead of 117.385 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Ksi(i) = EXP(-8904.2d0*invtk  + 117.400d0 - 19.334d0*dlogtk +  &
                 (-458.79d0*invtk + 3.5913d0) * sqrtis +             &
                 (188.74d0*invtk - 1.5998d0) * is +                  &
                 (-12.1652d0*invtk + 0.07871d0) * is2 +              &
                 LOG(1.0 - 0.001005d0*s))

!       Kw = [H][OH]
!       (seawater scale)
!       Millero (1995) p.670, eq. 63 from composite data
!       Use Millero equation's 148.9802 constant instead of 148.9652 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Kw(i) = EXP(-13847.26d0*invtk + 148.9802d0 - 23.6521d0*dlogtk +  &
               (118.67d0*invtk - 5.977d0 + 1.0495d0 * dlogtk) *          &
               sqrts - 0.01615d0 * s)

!       Ks = [H][SO4]/[HSO4]
!       (free scale)
!       Dickson (1990, J. chem. Thermodynamics 22, 113)
        Ks_0p = EXP(-4276.1d0*invtk + 141.328d0 - 23.093d0*dlogtk          &
                + (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk) * sqrtis  &
                + (35474.d0*invtk - 771.54 + 114.723d0*dlogtk) * is      &
                - 2698.d0*invtk*is**1.5 + 1776.d0*invtk*is2              &
                + LOG(1.0d0 - 0.001005d0*s))

!       Kf = [H][F]/[HF]
!       (total scale)
        IF (trim(opKf) == 'dg') THEN
!          Dickson and Riley (1979) -- change pH scale to total (following Dickson & Goyet, 1994)
           Kf_0p = EXP(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis +  &
                   LOG(1.0d0 - 0.001005d0*s) +                     &
                   LOG(1.0d0 + St(i)/Ks_0p))
        ELSEIF (trim(opKf) == 'pf') THEN
!          Perez and Fraga (1987) - Already on Total scale (no need for last line above)
!          Formulation as given in Dickson et al. (2007)
           Kf_0p = EXP(874.d0*invtk - 9.68d0 + 0.111d0*sqrts)
        ELSE
           PRINT *, "optKf must be either 'dg' or 'pf'"
           STOP
        ENDIF

!       Kspc (calcite) - apparent solubility product of calcite
!       (no scale)
!       Kspc = [Ca2+] [CO32-] when soln is in equilibrium w/ calcite
!       Mucci 1983 mol/kg-soln
        Kspc(i) = 10d0**(-171.9065d0 - 0.077993d0*tk + 2839.319d0/tk    &
                 + 71.595d0*LOG10(tk)                             &
                 + (-0.77712d0 + 0.0028426d0*tk + 178.34d0/tk)*sqrts  &
                 -0.07711d0*s + 0.0041249d0*s15 )


!       Kspa (aragonite) - apparent solubility product of aragonite
!       (no scale)
!       Kspa = [Ca2+] [CO32-] when soln is in equilibrium w/ aragonite
!       Mucci 1983 mol/kg-soln
        Kspa(i) = 10.d0**(-171.945d0 - 0.077993d0*tk + 2903.293d0/tk &
             +71.595d0*LOG10(tk) &
             +(-0.068393d0 + 0.0017276d0*tk + 88.135d0/tk)*sqrts &
             -0.10018d0*s + 0.0059415d0*s15 )

!       Pressure effect on K0 based on Weiss (1974, equation 5)
        Rgas_atm = 82.05736_wp      ! (cm3 * atm) / (mol * K)  CODATA (2006)
        vbarCO2 = 32.3_wp           ! partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
        K0(i) = K0(i) * exp( ((1-Ptot)*vbarCO2)/(Rgas_atm*tk0) )   ! Weiss (1974, equation 5)

!       Pressure effect on all other K's (based on Millero, (1995)
!           index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
!                  K1p(9), K2p(10), K3p(11), Ksi(12)
        DO ipc = 1, 12
           deltav(ipc)  =  a0(ipc) + a1(ipc) *t + a2(ipc) *t*t
           deltak(ipc)   = (b0(ipc)  + b1(ipc) *t + b2(ipc) *t*t)
           lnkpok0(ipc)  = (-(deltav(ipc)) &
                +(0.5d0*deltak(ipc) * prb) &
                )                         * prb/(R*tk)
        END DO

!       Pressure correction on Ks (Free scale)
        Ks(i) = Ks_0p*EXP(lnkpok0(5))
!       Conversion factor total -> free scale
        total2free     = 1.d0/(1.d0 + St(i)/Ks(i))   ! Kfree = Ktotal*total2free
!       Conversion factor total -> free scale at pressure zero
        total2free_0p  = 1.d0/(1.d0 + St(i)/Ks_0p)   ! Kfree = Ktotal*total2free

!       Pressure correction on Kf
!       Kf must be on FREE scale before correction
        Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
        Kf(i) = Kf_0p * EXP(lnkpok0(6)) !Pressure correction (on Free scale)
        Kf(i) = Kf(i)/total2free        !Convert back from Free to Total scale

!       Convert between seawater and total hydrogen (pH) scales
        free2SWS  = 1.d0 + St(i)/Ks(i) + Ft(i)/(Kf(i)*total2free)  ! using Kf on free scale
        total2SWS = total2free * free2SWS                          ! KSWS = Ktotal*total2SWS
        SWS2total = 1.d0 / total2SWS
!       Conversion at pressure zero
        free2SWS_0p  = 1.d0 + St(i)/Ks_0p + Ft(i)/(Kf_0p)  ! using Kf on free scale
        total2SWS_0p = total2free_0p * free2SWS_0p         ! KSWS = Ktotal*total2SWS

!       Convert from Total to Seawater scale before pressure correction
!       Must change to SEAWATER scale: K1, K2, Kb
        IF (trim(optK1K2) == 'l') THEN
          K1(i)  = K1(i)*total2SWS_0p
          K2(i)  = K2(i)*total2SWS_0p
          !This conversion is unnecessary for the K1,K2 from Millero (2010),
          !since we use here the formulation already on the seawater scale
        ENDIF
        Kb(i)  = Kb(i)*total2SWS_0p

!       Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw

!       Other contants (keep on another scale):
!          - K0         (independent of pH scale, already pressure corrected)
!          - Ks         (already on Free scale;   already pressure corrected)
!          - Kf         (already on Total scale;  already pressure corrected)
!          - Kspc, Kspa (independent of pH scale; pressure-corrected below)

!       Perform actual pressure correction (on seawater scale)
        K1(i)   = K1(i)*EXP(lnkpok0(1))
        K2(i)   = K2(i)*EXP(lnkpok0(2))
        Kb(i)   = Kb(i)*EXP(lnkpok0(3))
        Kw(i)   = Kw(i)*EXP(lnkpok0(4))
        Kspc(i) = Kspc(i)*EXP(lnkpok0(7))
        Kspa(i) = Kspa(i)*EXP(lnkpok0(8))
        K1p(i)  = K1p(i)*EXP(lnkpok0(9))
        K2p(i)  = K2p(i)*EXP(lnkpok0(10))
        K3p(i)  = K3p(i)*EXP(lnkpok0(11))
        Ksi(i)  = Ksi(i)*EXP(lnkpok0(12))

!       Convert back to original total scale:
        K1(i)  = K1(i) *SWS2total
        K2(i)  = K2(i) *SWS2total
        K1p(i) = K1p(i)*SWS2total
        K2p(i) = K2p(i)*SWS2total
        K3p(i) = K3p(i)*SWS2total
        Kb(i)  = Kb(i) *SWS2total
        Ksi(i) = Ksi(i)*SWS2total
        Kw(i)  = Kw(i) *SWS2total

     ELSE

        K0(i)   = 1.e20_wp
        K1(i)   = 1.e20_wp
        K2(i)   = 1.e20_wp
        Kb(i)   = 1.e20_wp
        Kw(i)   = 1.e20_wp
        Ks(i)   = 1.e20_wp
        Kf(i)   = 1.e20_wp
        Kspc(i) = 1.e20_wp
        Kspa(i) = 1.e20_wp
        K1p(i)  = 1.e20_wp
        K2p(i)  = 1.e20_wp
        K3p(i)  = 1.e20_wp
        Ksi(i)  = 1.e20_wp
        Bt(i)   = 1.e20_wp
        Ft(i)   = 1.e20_wp
        St(i)   = 1.e20_wp

     ENDIF

  END DO

  RETURN
END SUBROUTINE constants

! ----------------------------------------------------------------------
!  VARSOLVER
! ----------------------------------------------------------------------
!
!> \file varsolver.f90
!! \BRIEF 
!> Module with varsolver subroutine - solve for pH and other carbonate system variables
!>    Solve for pH and other carbonate system variables (with input from vars routine)
SUBROUTINE varsolver(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC,             &
                    temp, salt, ta, tc, pt, sit,                                 &
                    Bt, St, Ft,                                                  &
                    K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,  & 
                    Patm, Phydro_bar, rhodum, optGAS                             )

  !   Purpose: Solve for pH and other carbonate system variables (with input from vars routine)

  !     INPUT variables:
  !     ================
  !     temp    = in situ temperature [degrees C]
  !     ta      = total alkalinity                     in [eq/m^3] or [eq/kg]   based on optCON in calling routine (vars)
  !     tc      = dissolved inorganic carbon           in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     pt      = total dissolved inorganic phosphorus in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     sit     = total dissolved inorganic silicon    in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     Bt      = total dissolved inorganic boron      computed in calling routine (vars)
  !     St      = total dissolved inorganic sulfur     computed in calling routine (vars)
  !     Ft      = total dissolved inorganic fluorine   computed in calling routine (vars)
  !     K's     = K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi 
  !     Patm    = atmospheric pressure [atm]
  !     Phydro_bar = hydrostatic pressure [bar]
  !     rhodum  = density factor as computed in calling routine  (vars)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2 (default optGAS = 'Pinsitu')
  !     ---------
  !       PRESSURE & T corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     ph   = pH on total scale
  !     pco2 = CO2 partial pressure (uatm)
  !     fco2 = CO2 fugacity (uatm)
  !     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] determined by rhodum (depends on optCON in calling routine)
  !     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state

  USE mocsy_singledouble
  USE mocsy_phsolvers

  IMPLICIT NONE

! Input variables
  !> <b>in situ temperature</b> [degrees C]
  REAL(kind=wp), INTENT(in) :: temp
  !> <b>salinity</b> [on the practical salinity scale, dimensionless]
  REAL(kind=wp), INTENT(in) :: salt
  !> total alkalinity in <b>[eq/m^3]</b> OR in <b>[eq/kg]</b>, depending on optCON in calling routine
  REAL(kind=wp), INTENT(in) :: ta
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=wp), INTENT(in) :: tc
  !> phosphate concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=wp), INTENT(in) :: pt
  !> total dissolved inorganic silicon concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=wp), INTENT(in) :: sit
  !> total boron from either Uppstrom (1974) or Lee et al. (2010), depending on optB in calling routine
  REAL(kind=wp), INTENT(in) :: Bt
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=wp), INTENT(in) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=wp), INTENT(in) :: Ft
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=wp), INTENT(in) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=wp), INTENT(in) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=wp), INTENT(in) :: K2
  !> equilibrium constant for dissociation of boric acid 
  REAL(kind=wp), INTENT(in) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=wp), INTENT(in) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=wp), INTENT(in) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! from Dickson and Riley (1979) or Perez and Fraga (1987), depending on optKf
  REAL(kind=wp), INTENT(in) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=wp), INTENT(in) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=wp), INTENT(in) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=wp), INTENT(in) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=wp), INTENT(in) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=wp), INTENT(in) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=wp), INTENT(in) :: Ksi
  !> total atmospheric pressure <b>[atm]</b>
  REAL(kind=wp), INTENT(in) :: Patm
  !> total hydrostatic pressure <b>[bar]</b>
  REAL(kind=wp), INTENT(in) :: Phydro_bar
  !> density factor as computed incalling routine  (vars)
  REAL(kind=wp), INTENT(in) :: rhodum
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(kind=wp), INTENT(out) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=wp), INTENT(out) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=wp), INTENT(out) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=wp), INTENT(out) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=wp), INTENT(out) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=wp), INTENT(out) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=wp), INTENT(out) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=wp), INTENT(out) :: OmegaC

! Local variables
  REAL(kind=wp) :: Phydro_atm, Ptot
  REAL(kind=wp) :: Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=wp) :: tk, tk0
  real(kind=wp) :: temp68, tempot, tempot68
  REAL(kind=wp) :: Hinit, H
  REAL(kind=wp) :: HSO4, HF, HSI, HPO4
  REAL(kind=wp) :: ab, aw, ac
  REAL(kind=wp) :: cu, cb, cc
  REAL(kind=wp) :: Ca
! Array to pass optional arguments
  CHARACTER(7) :: opGAS

  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

! Compute pH from constants and total concentrations
! - use SolveSAPHE v1.0.1 routines from Munhoven (2013, GMD) modified to use mocsy's Ks instead of its own
! 1) Compute best starting point for H+ calculation
  call ahini_for_at(ta, tc, Bt, K1, K2, Kb, Hinit)
! 2) Solve for H+ using above result as the initial H+ value
  H = solve_at_general(ta, tc, Bt,                                         & 
                       pt,     sit,                                        &
                       St, Ft,                                             &
                       K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,     &
                       Hinit)
! 3) Calculate pH from from H+ concentration (mol/kg)
  IF (H > 0.d0) THEN
     pH = -1.*LOG10(H)
  ELSE
     pH = 1.e20_wp
  ENDIF

! Compute carbonate Alk (Ac) by difference: from total Alk and other Alk components
  HSO4 = St/(1.0d0 + Ks/(H/(1.0d0 + St/Ks)))
  HF = 1.0d0/(1.0d0 + Kf/H)
  HSI = 1.0d0/(1.0d0 + H/Ksi)
  HPO4 = (K1p*K2p*(H + 2.*K3p) - H**3) /                &
  (H**3 + K1p*H**2 + K1p*K2p*H + K1p*K2p*K3p)
  ab = Bt/(1.0d0 + H/Kb)
  aw = Kw/H - H/(1.0d0 + St/Ks)
  ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4

! Calculate CO2*, HCO3-, & CO32- (in mol/kg soln) from Ct, Ac, H+, K1, & K2
  cu = (2.0d0 * tc - ac) / (2.0d0 + K1 / H)
  cb = K1 * cu / H
  cc = K2 * cb / H

! When optCON = 'mol/m3' in calling routine (vars), then:
! convert output var concentrations from mol/kg to mol/m^3
! e.g., for case when drho = 1028, multiply by [1.028 kg/L  x  1000 L/m^3])
  co2  = cu * rhodum
  hco3 = cb * rhodum
  co3  = cc * rhodum

! Determine CO2 fugacity [uatm]
! NOTE: equation just below requires CO2* in mol/kg
  fCO2 = cu * 1.e6_wp/K0

! Determine CO2 partial pressure from CO2 fugacity [uatm]
  tk = 273.15d0 + temp
  !Compute EITHER "potential pCO2" OR "in situ pCO2" (T and P used for calculations will differ)
  IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
     tk0 = tk                 !in situ temperature (K) for K0 calculation
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Ppot' .OR. trim(opGAS) == 'ppot') THEN
     !Use potential temperature and atmospheric pressure (water parcel adiabatically brought back to surface)
     !temp68 = (temp - 0.0002d0) / 0.99975d0          !temp = in situ T; temp68 is same converted to ITPS-68 scale
     !tempot68 = sw_ptmp(salt, temp68, Phydro_bar*10d0, 0.0d0) !potential temperature (C)
     !tempot   = 0.99975*tempot68 + 0.0002
     !tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     tempot = sw_ptmp(salt, temp, Phydro_bar*10._wp, 0.0_wp) !potential temperature (C)
     tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
     !Use in situ temperature and total pressure 
     tk0 = tk                             !in situ temperature (K) for fugacity coefficient calculation
     Phydro_atm = Phydro_bar / 1.01325d0  !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
     Ptot = Patm + Phydro_atm            !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
  ELSE
     PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
     STOP
  ENDIF

! Now that we have T and P in the right form, continue with calculation of fugacity coefficient (and pCO2)
  Rgas_atm = 82.05736_wp      ! (cm3 * atm) / (mol * K)  CODATA (2006)
! To compute fugcoeff, we need 3 other terms (B, Del, xc2) in addition to 3 others above (tk, Ptot, Rgas_atm)
  B = -1636.75d0 + 12.0408d0*tk0 - 0.0327957d0*(tk0*tk0) + 0.0000316528d0*(tk0*tk0*tk0)
  Del = 57.7d0 - 0.118d0*tk0
! "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
! x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
! Let's assume that xCO2 = fCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
  xCO2approx = fCO2 * 1.e-6_wp
  IF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
!    xCO2approx = 400.0e-6_wp      !a simple test (gives about same result as seacarb for pCO2insitu)
!    approximate surface xCO2 ~ surface fCO2 (i.e., in situ fCO2 d by exponential pressure correction)
     xCO2approx = xCO2approx * exp( ((1-Ptot)*32.3_wp)/(82.05736_wp*tk0) )   ! of K0 press. correction, see Weiss (1974, equation 5)
  ENDIF
  xc2 = (1.0d0 - xCO2approx)**2 
  fugcoeff = exp( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk0) )
  pCO2 = fCO2 / fugcoeff

! Determine Omega Calcite et Aragonite
! OmegaA = ((0.01028d0*salt/35.0d0)*cc) / Kspa
! OmegaC = ((0.01028d0*salt/35.0d0)*cc) / Kspc
! - see comments from Munhoven on the best value "0.02128" which differs slightly from the best practices guide (0.02127)
  Ca = (0.02128d0/40.078d0) * salt/1.80655d0
  OmegaA = (Ca*cc) / Kspa
  OmegaC = (Ca*cc) / Kspc

  RETURN
END SUBROUTINE varsolver

! ----------------------------------------------------------------------
!  VARS
! ----------------------------------------------------------------------
!
!> \file vars.f90
!! \BRIEF 
!> Module with vars subroutine - compute carbonate system vars from DIC,Alk,T,S,P,nuts
!>    Computes standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!!    as 1D arrays FROM
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations (all 1-D arrays)
SUBROUTINE vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                      &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS                          )

  !   Purpose:
  !     Computes other standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
  !     as 1D arrays
  !     FROM:
  !     temperature, salinity, pressure,
  !     total alkalinity (ALK), dissolved inorganic carbon (DIC),
  !     silica and phosphate concentrations (all 1-D arrays)

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
  !     alk     = total alkalinity in [eq/m^3] with optCON = 'mol/m3'
  !             =               [eq/kg]  with optCON = 'mol/kg'
  !     dic     = dissolved inorganic carbon [mol/m^3] with optCON = 'mol/m3'
  !             =                            [mol/kg]  with optCON = 'mol/kg'
  !     sil     = silica    [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     phos    = phosphate [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     INPUT options:
  !     ==============
  !     -----------
  !     optCON: choose input & output concentration units - mol/kg (data) vs. mol/m^3 (models)
  !     -----------
  !       -> 'mol/kg' for DIC, ALK, sil, & phos given on mokal scale, i.e., in mol/kg  (std DATA units)
  !       -> 'mol/m3' for DIC, ALK, sil, & phos given in mol/m^3 (std MODEL units)
  !     -----------
  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (thus p = depth)
  !     ---------
  !     optB: choose total boron formulation - Uppström (1974) vs. Lee et al. (2010)
  !     ---------
  !       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
  !       -> 'l10' means use newer   formulation of Lee et al. (2010) for total Boron
  !     ---------
  !     optK1K2:
  !     ---------
  !       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
  !                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
  !       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
  !     ----------
  !     optKf:
  !     ----------
  !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
  !               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
  !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2
  !     ---------
  !       PRESSURE corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     ph   = pH on total scale
  !     pco2 = CO2 partial pressure (uatm)
  !     fco2 = CO2 fugacity (uatm)
  !     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state
  !     BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  !     rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  !     p = pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  !     tempis  = in-situ temperature [degrees C]

  USE mocsy_singledouble

  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=wp), INTENT(in),    DIMENSION(N) :: temp
  !> salinity <b>[psu]</b>
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: sal
  !> total alkalinity in <b>[eq/m^3]</b> (when optCON = 'mol/m3') OR in <b>[eq/kg]</b>  (when optCON = 'mol/kg')
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: alk
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: dic
  !> SiO2 concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: sil
  !> phosphate concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: phos
!f2py optional , depend(sal) :: n=len(sal)
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: Patm
  !> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(kind=wp), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=wp), INTENT(in),    DIMENSION(N) :: lat

  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
  CHARACTER(6), INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
  CHARACTER(7), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(2), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20), 
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
!f2py character*3 optional, intent(in) :: optB='l10'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: OmegaC
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: tempis

! Local variables
  REAL(kind=wp) :: ssal, salk, sdic, ssil, sphos
  REAL(kind=wp) :: tempot, tempis68, tempot68
  REAL(kind=wp) :: drho

  REAL(kind=wp) :: K0, K1, K2, Kb, Kw, Ks, Kf, Kspc
  REAL(kind=wp) :: Kspa, K1p, K2p, K3p, Ksi
  REAL(kind=wp) :: St, Ft, Bt

  REAL(kind=wp), DIMENSION(1) :: aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc
  REAL(kind=wp), DIMENSION(1) :: aKspa, aK1p, aK2p, aK3p, aKsi
  REAL(kind=wp), DIMENSION(1) :: aSt, aFt, aBt

  REAL(kind=wp) :: Patmd, Ptot, Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=wp) :: Phydro_atm

  INTEGER :: i, icount

  REAL(kind=wp) :: t, tk, prb
  REAL(kind=wp) :: s
  REAL(kind=wp) :: tc, ta
  REAL(kind=wp) :: sit, pt
  REAL(kind=wp) :: Hinit
  REAL(kind=wp) :: ah1

  REAL(kind=wp) :: HSO4, HF, HSI, HPO4
  REAL(kind=wp) :: ab, aw, ac, ah2, erel

  REAL(kind=wp) :: cu, cb, cc

  REAL(kind=wp), DIMENSION(2) :: dicdel, pco2del
  REAL(kind=wp) :: dx, Rf
  REAL(kind=wp) :: dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC

  INTEGER :: kcomp
  INTEGER :: j, minusplus

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
!   print *,"optB present:"
!   print *,"optB = ", optB 
    opB = optB
  ELSE
!   Default is Lee et al (2010) for total boron
!   print *,"optB NOT present:"
    opB = 'l10'
!   print *,"opB = ", opB 
  ENDIF
  IF (PRESENT(optKf)) THEN
!   print *,"optKf = ", optKf
    opKf = optKf
  ELSE
!   print *,"optKf NOT present:"
!   Default is Perez & Fraga (1987) for Kf
    opKf = 'pf'
!   print *,"opKf = ", opKf
  ENDIF
  IF (PRESENT(optK1K2)) THEN
!   print *,"optK1K2 = ", optK1K2
    opK1K2 = optK1K2
  ELSE
!   print *,"optK1K2 NOT present:"
!   Default is Lueker et al. 2000) for K1 & K2
    opK1K2 = 'l'
!   print *,"opK1K2 = ", opK1K2
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

  icount = 0
  DO i = 1, N
     icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Model temperature tracer is usually "potential temperature"
!    * Model vertical grid is usually in meters
!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry, if appropriate,
!    convert these 2 model vars (input to this routine)
!    - depth [m] => convert to pressure [db]
!    - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     !print *,"optP =", optP, "end"
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p(i) = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db') THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p(i) = depth(i)
     ELSE
        !print *,"optP =", optP, "end"
        PRINT *,"optP must be either 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = temp(i)
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002) / 0.99975
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p(i), 0.0_wp )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis(i) = 0.99975*tempis68 + 0.0002
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis(i) = temp(i)
        tempis68  = (temp(i) - 0.0002) / 0.99975
!       dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0d0)
!       dtempot   = 0.99975*dtempot68 + 0.0002
     ELSE
        PRINT *,"optT must be either 'Tpot' or 'Tinsitu'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF

!    ================================================================
!    Carbonate chemistry computations
!    ================================================================
     IF (dic(i) > 0. .AND. dic(i) < 1.0e+4) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (       sal(i) < 0.   &
             .OR.  alk(i) < 0.   &
             .OR.  dic(i) < 0.   &
             .OR.  sil(i) < 0.   &
             .OR. phos(i) < 0.   &
             .OR.  sal(i) > 1e+3 &
             .OR.  alk(i) > 1e+3 &
             .OR.  dic(i) > 1e+3 &
             .OR.  sil(i) > 1e+3 &
             .OR. phos(i) > 1e+3) THEN
           PRINT *, 'i, icount, tempot, sal,    alk,    dic,    sil,    phos =', &
                     i, icount, tempot, sal(i), alk(i), dic(i), sil(i), phos(i)
        ENDIF
!       Zero out any negative salinity, phosphate, silica, dic, and alk
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF
        IF (phos(i) < 0.0) THEN
           sphos = 0.0
        ELSE
           sphos = phos(i)
        ENDIF
        IF (sil(i) < 0.0) THEN
           ssil = 0.0
        ELSE
           ssil = sil(i)
        ENDIF
        IF (dic(i) < 0.0) THEN
          sdic = 0.0
        ELSE
          sdic = dic(i)
        ENDIF
        IF (alk(i) < 0.0) THEN
          salk = 0.0
        ELSE
          salk = alk(i)
        ENDIF

!       Absolute temperature (Kelvin) & related variables
        t  = DBLE(tempis(i))
        tk = 273.15d0 + t

!       Atmospheric pressure
        Patmd = DBLE(Patm(i))
!       Hydrostatic pressure (prb is in bars)
        prb = DBLE(p(i)) / 10.0d0
        Phydro_atm = prb / 1.01325d0  ! convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
!       Total pressure [atm]
        IF     (trim(opGAS) == 'Pzero'     .OR. trim(opGAS) == 'pzero') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           Ptot = Patmd + Phydro_atm   ! total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF

!       Salinity (equivalent array in double precision)
        s = DBLE(ssal)

!       Get all equilibrium constants and total concentrations of SO4, F, B
        CALL constants(aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc, aKspa,  &
                       aK1p, aK2p, aK3p, aKsi,                           &
                       aSt, aFt, aBt,                                    &
                       temp(i), sal(i), Patm(i),                         &
                       depth(i), lat(i), 1,                              &
                       optT, optP, opB, opK1K2, opKf, opGAS              )

!       Unlike f77, in F90 we can't assign an array (dimen=1) to a scalar in a routine argument
!       Thus, set scalar constants equal to array (dimension=1) values required as arguments
        K0 = aK0(1) ; K1 = aK1(1) ; K2 = aK2(1) ; Kb = aKb(1) ; Kw = aKw(1) 
        Ks = aKs(1) ; Kf = aKs(1) ; Kspc = aKspc(1) ; Kspa = aKspa(1) 
        K1p = aK1p(1) ; K2p = aK2p(1) ; K3p = aK3p(1) ; Ksi = aKsi(1)
        St = aSt(1) ; Ft = aFt(1) ; Bt = aBt(1)

!       Compute in-situ density [kg/m^3]
        rhoSW(i) =  rho(ssal, tempis68, prb)

!       Either convert units of DIC and ALK (MODEL case) or not (DATA case)
        IF     (trim(optCON) == 'mol/kg') THEN
!          No conversion:
!          print *,'DIC and ALK already given in mol/kg (std DATA units)'
           drho = 1.
        ELSEIF (trim(optCON) == 'mol/m3') THEN
!          Do conversion:
!          print *,"DIC and ALK given in mol/m^3 (std MODEL units)"
           drho = DBLE(rhoSW(i))
        ELSE
           PRINT *,"optCON must be either 'mol/kg' or 'mol/m3'"
           STOP
        ENDIF

        tc  = DBLE(sdic)/drho
        ta  = DBLE(salk)/drho
        sit = DBLE(ssil)/drho
        pt  = DBLE(sphos)/drho

!       Solve for pH and all other variables
!       ------------------------------------
        CALL varsolver(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC,     &
                      t, s, ta, tc, pt, sit,                                       &
                      Bt, St, Ft,                                                  &
                      K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,  & 
                      Patmd, prb, drho, opGAS                                     )

!       Convert all output variables from double to single precision
        pH(i)     = REAL(dph)
        co2(i)    = REAL(dco2)
        hco3(i)   = REAL(dhco3)
        co3(i)    = REAL(dco3)
        fCO2(i)   = REAL(dfCO2)
        pCO2(i)   = REAL(dpCO2)
        OmegaA(i) = REAL(dOmegaA)
        OmegaC(i) = REAL(dOmegaC)

!       Compute Revelle factor numerically (derivative using centered-difference scheme)
        DO j=1,2
           minusplus = (-1)**j
           dx = 0.1 * 1e-6         ! Numerical tests found for DIC that optimal dx = 0.1 umol/kg (0.1e-6 mol/kg)
           dicdel(j) = tc + DBLE(minusplus)*dx/2.0d0
            CALL varsolver(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC, &
               t, s, ta, dicdel(j), pt, sit,                                       &
               Bt, St, Ft,                                                         &
               K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,         & 
               Patmd, prb, drho, optGAS                                            )
            pco2del(j) = dpco2
        END DO
       !Classic finite centered difference formula for derivative (2nd order accurate)
        Rf = (pco2del(2) - pco2del(1)) / (dicdel(2) - dicdel(1))       ! dpCO2/dDIC
       !Rf = (pco2del(2) - pco2del(1)) / (dx)                          ! dpCO2/dDIC (same as just above)
        Rf = Rf * tc / dpco2                                           ! R = (dpCO2/dDIC) * (DIC/pCO2)

        BetaD(i) = REAL(Rf)

     ELSE

        ph(i)     = 1.e20_wp
        pco2(i)   = 1.e20_wp
        fco2(i)   = 1.e20_wp
        co2(i)    = 1.e20_wp
        hco3(i)   = 1.e20_wp
        co3(i)    = 1.e20_wp
        OmegaA(i) = 1.e20_wp
        OmegaC(i) = 1.e20_wp
        BetaD(i)  = 1.e20_wp
        rhoSW(i)  = 1.e20_wp
        p(i)      = 1.e20_wp
        tempis(i) = 1.e20_wp

     ENDIF

  END DO

  RETURN
END SUBROUTINE vars

! ----------------------------------------------------------------------
!  P2FCO2
! ----------------------------------------------------------------------
!
!> \file p2fCO2.f90
!! \BRIEF 
!>    Module with p2fCO2 subroutine - compute fCO2 from pCO2, in situ T, atm pressure, hydrostatic pressure
!>    Compute fCO2 from arrays of pCO2, in situ temp, atm pressure, & hydrostatic pressure
SUBROUTINE p2fCO2(pCO2, temp, Patm, p, N, fCO2)
  !    Purpose:
  !    Compute fCO2 from arrays of pCO2, in situ temp, atm pressure, & hydrostatic pressure

  USE mocsy_singledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, INTENT(in) :: N

! INPUT variables
  !> oceanic partial pressure of CO2 [uatm]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: pCO2
  !> in situ temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: temp
  !> atmospheric pressure [atm]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: Patm
  !> hydrostatic pressure [db]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: p

! OUTPUT variables:
  !> fugacity of CO2 [uatm] 
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: fCO2

! LOCAL variables:
  REAL(kind=wp) :: dpCO2, dtemp, tk, dPatm, prb
  REAL(kind=wp) :: Ptot, Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=wp) :: dfCO2

  INTEGER :: i

! REAL(kind=wp) :: sw_ptmp
! EXTERNAL sw_ptmp

  DO i = 1,N
     dpCO2     = DBLE(pCO2(i))
     dtemp     = DBLE(temp(i))
     dPatm     = DBLE(Patm(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     prb = DBLE(p(i)) / 10.0d0         !Pressure effect (prb is in bars)
     Ptot = dPatm + prb/1.01325d0      !Total pressure (atmospheric + hydrostatic) [atm]
     Rgas_atm = 82.05736_wp            !R in (cm3 * atm) / (mol * K)  from CODATA (2006)
!    To compute fugcoeff, we need 3 other terms (B, Del, xc2) as well as 3 others above (tk, Ptot, Rgas_atm)
     B = -1636.75d0 + 12.0408d0*tk - 0.0327957d0*(tk*tk) + 0.0000316528d0*(tk*tk*tk)
     Del = 57.7d0 - 0.118d0*tk
!    "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
!    x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
!    Let's assume that xCO2 = pCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
     xCO2approx = dpCO2 * 1.e-6_wp
     xc2 = (1.0d0 - xCO2approx)**2 
     fugcoeff = EXP( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk) )
     dfCO2 = dpCO2 * fugcoeff
     fCO2(i) = REAL(dfCO2)
  END DO

  RETURN
END SUBROUTINE p2fCO2

! ----------------------------------------------------------------------
!  P2FCO2
! ----------------------------------------------------------------------
!
!> \file f2pCO2.f90
!! \BRIEF 
!>    Module with f2pCO2 subroutine - compute pCO2 from fCO2, in situ T, atm pressure, hydrostatic pressure
!>    Compute pCO2 from arrays of fCO2, in situ temp, atm pressure, & hydrostatic pressure
SUBROUTINE f2pCO2(fCO2, temp, Patm, p, N, pCO2)
  !    Purpose:
  !    Compute pCO2 from arrays of fCO2, in situ temp, atm pressure, & hydrostatic pressure

  USE mocsy_singledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> oceanic fugacity of CO2 [uatm]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: fCO2
  !> in situ temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: temp
  !> atmospheric pressure [atm]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: Patm
  !> hydrostatic pressure [db]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: p

! OUTPUT variables:
  !> oceanic partial pressure of CO2 [uatm] 
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: pCO2

! LOCAL variables:
  REAL(kind=wp) :: dfCO2, dtemp, tk, dPatm, prb
  REAL(kind=wp) :: Ptot, Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=wp) :: dpCO2

  INTEGER :: i

! REAL(kind=wp) :: sw_ptmp
! EXTERNAL sw_ptmp

  DO i = 1,N
     dfCO2     = DBLE(fCO2(i))
     dtemp     = DBLE(temp(i))
     dPatm     = DBLE(Patm(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     prb = DBLE(p(i)) / 10.0d0         !Pressure effect (prb is in bars)
     Ptot = dPatm + prb/1.01325d0      !Total pressure (atmospheric + hydrostatic) [atm]
     Rgas_atm = 82.05736_wp            !R in (cm3 * atm) / (mol * K)  from CODATA (2006)
!    To compute fugcoeff, we need 3 other terms (B, Del, xc2) as well as 3 others above (tk, Ptot, Rgas_atm)
     B = -1636.75d0 + 12.0408d0*tk - 0.0327957d0*(tk*tk) + 0.0000316528d0*(tk*tk*tk)
     Del = 57.7d0 - 0.118d0*tk
!    "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
!    x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
!    Let's assume that xCO2 = fCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
     xCO2approx = dfCO2 * 1.e-6_wp
     xc2 = (1.0d0 - xCO2approx)**2 
     fugcoeff = exp( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk) )
     dpCO2 = dfCO2 / fugcoeff
     pCO2(i) = REAL(dpCO2)
  END DO

  RETURN
END SUBROUTINE f2pCO2

END MODULE mocsy_mainmod
