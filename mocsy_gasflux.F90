MODULE mocsy_gasflux

CONTAINS

! --------------------------------------------------------------------
!  Schmidt CO2 number
! --------------------------------------------------------------------
!
! Title  : Calculates Schmidt number for ocean uptake of CO2
! Author : Andrew Yool
! Date   : 14/10/04
!
! This function calculates the Schmidt number for CO2 using sea surface
! temperature.  The code is based upon that developed as part of the
! OCMIP-2 project (1998-2000).  The coefficients used are taken from
! Wanninkhof (1992, JGR, 97, 7373-7382).
!
! AXY (12/06/2015)
! UPDATED: coefficients used below are now those from Wanninkhof (2014)
! update to original 1992 paper. Full reference is:
!
! Winninkhof, R. (2014). Relationship between wind speed and gas 
! exchange over the ocean revisited. LIMNOLOGY AND OCEANOGRAPHY-METHODS
! 12, 351-362, doi:10.4319/lom.2014.12.351
!
! Check answer for the function at 20 degrees C is 668.
!
! Function inputs are (in order) :
!     t            temperature (degrees C)
! (*) co2_schmidt  carbon dioxide Schmidt number
!
! Where (*) is the function output.
!
      subroutine schmidt_co2(pt, N, co2_schmidt)

      USE mocsy_singledouble

      implicit none
!
      INTEGER, INTENT(in) :: N
      real(kind=wp), INTENT(in),  DIMENSION(N) :: pt
      real(kind=wp), INTENT(out), DIMENSION(N) :: co2_schmidt
!
      real(kind=wp)              :: a0, a1, a2, a3, a4
!
!     data a0 /    2073.1 /
!     data a1 /   -125.62 /
!     data a2 /    3.6276 /
!     data a3 / -0.043219 /
!
      data a0 /    2116.8 /
      data a1 /   -136.25 /
      data a2 /    4.7353 /
      data a3 / -0.092307 /
      data a4 / 0.0007555 /
!
! Wanninkhof (1992)
!     co2_schmidt = a0 + pt*(a1 + pt*(a2 + pt*a3))
!
! Wanninkhof (2014) adds in an extra term
      co2_schmidt = a0 + pt*(a1 + pt*(a2 + pt*(a3 + pt*a4)))
!
      return

      end subroutine schmidt_co2

! --------------------------------------------------------------------
!  Surface K0
! --------------------------------------------------------------------
!
! Title  : Calculates surface K0 from surface T & S
! Author : Andrew Yool
! Date   : 18/06/15
!
! This function is derived from code included in the MOCSY package
! produced by Jim Orr.
!
      subroutine surface_K0(ptmp, saln, N, K0)

      USE mocsy_singledouble

      implicit none
!
      INTEGER, INTENT(in) :: N
      real(kind=wp), INTENT(in),  DIMENSION(N) :: ptmp, saln
      real(kind=wp), INTENT(out), DIMENSION(N) :: K0
!
      real(kind=wp), DIMENSION(N) :: tk, invtk, tmp
      real(kind=wp)               :: a0, a1, a2, a3, a4
!
      tk    = ptmp + 273.15d0
      invtk = 1.0d0 / tk
      tmp = (9345.17d0*invtk) - 60.2409d0 + (23.3585d0 * LOG(tk/100.0d0))
      K0 = EXP( tmp + saln*(0.023517d0 - (0.00023656d0*tk) + (0.0047036e-4_wp*tk*tk)) )
!
      return

      end subroutine surface_K0

! --------------------------------------------------------------------
!  Calculate xCO2
! --------------------------------------------------------------------
!
!>    Compute xCO2 from arrays of pCO2atm, in situ T, S, & atm pressure
SUBROUTINE pCO2atm2xCO2(pCO2atm, temp, salt, Patm, N, xCO2)
  !    Purpose:
  !    Compute xCO2 from arrays of pCO2atm, in situ T, S, & atm pressure

  USE mocsy_singledouble

  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> atmospheric partial pressure of CO2 [uatm] 
  ! AXY (22/06/15): amended this next line to "in" as that's what it should be!
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: pCO2atm
  !> in situ temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: salt
  !> atmospheric pressure [atm]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: Patm
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> mole fraction of CO2 [ppm]
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: xCO2

! LOCAL variables:
  REAL(kind=wp) :: dpCO2atm, dPatm
  REAL(kind=wp), DIMENSION(N) :: pH20
  REAL(kind=wp) :: dxCO2

  INTEGER :: i

  call vapress(temp, salt, N, pH20)

  DO i = 1,N
     dpCO2atm  = DBLE(pCO2atm(i))
     dPatm     = DBLE(Patm(i))
     dxCO2     = dpCO2atm / (dPatm - pH20(i))
     xCO2(i) = REAL(dxCO2)
  END DO

  RETURN
END SUBROUTINE pCO2atm2xCO2

! --------------------------------------------------------------------
!  Calculate pCO2atm
! --------------------------------------------------------------------
!
!>    Compute pCO2atm from arrays of xCO2, in situ T, S, & atm pressure
SUBROUTINE x2pCO2atm(xCO2, temp, salt, Patm, N, pCO2atm)
  !    Purpose:
  !    Compute pCO2atm from arrays of xCO2, in situ T, S, & atm pressure

  USE mocsy_singledouble

  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> mole fraction of CO2 [ppm]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: xCO2
  !> in situ temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: salt
  !> atmospheric pressure [atm]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: Patm
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> oceanic partial pressure of CO2 [uatm] 
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: pCO2atm

! LOCAL variables:
  REAL(kind=wp) :: dxCO2, dPatm
  REAL(kind=wp), DIMENSION(N) :: pH20
  REAL(kind=wp) :: dpCO2atm

  INTEGER :: i

! Compute vapor pressure of seawater [in atm]
  call vapress(temp, salt, N, pH20)

  DO i = 1,N
     dxCO2     = DBLE(xCO2(i))
     dPatm     = DBLE(Patm(i))
     dpCO2atm = (dPatm - pH20(i)) * dxCO2
     pCO2atm(i) = REAL(dpCO2atm)
  END DO

  RETURN
END SUBROUTINE x2pCO2atm

! --------------------------------------------------------------------
!  Calculate seawater vapor pressure
! --------------------------------------------------------------------
!
!>    Compute vapor pressure of seawater (atm) following preocedure from Weiss & Price (1980)
SUBROUTINE vapress(temp, salt, N, vpsw)
  !    Purpose:
  !    Compute vapor pressure of seawater (atm) following preocedure from Weiss & Price (1980)

  USE mocsy_singledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> in situ temperature [C]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=wp), INTENT(in), DIMENSION(N) :: salt
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> vapor pressure of seawater [atm] 
  REAL(kind=wp), INTENT(out), DIMENSION(N) :: vpsw

! LOCAL variables:
  REAL(kind=wp) :: tk, dsalt

  INTEGER :: i

  DO i = 1,N
     dsalt = DBLE(salt(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     vpsw(i) = exp(24.4543d0 - 67.4509d0*(100.0d0/tk) - 4.8489d0*log(tk/100) - 0.000544d0*dsalt)
  END DO

  RETURN
END SUBROUTINE vapress

END MODULE mocsy_gasflux
