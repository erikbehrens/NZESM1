MODULE insitu_tem

  USE dom_oce  ! ocean space and time domain
  USE oce,            ONLY: tsn
  USE par_oce,        ONLY: jpi, jpj, jpk, jpkm1
  USE lbclnk,         ONLY: lbc_lnk
  USE lib_mpp         ! MPP library

  IMPLICIT NONE
  PRIVATE

  PUBLIC theta2t

  !! * Accessibility
  PUBLIC insitu_tem_alloc          ! routines called by step.F90

  REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: insitu_t

   !! * Substitutions
#  include "domzgr_substitute.h90"

CONTAINS

   REAL FUNCTION insitu_tem_alloc()
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ierr
      !!----------------------------------------------------------------------
      !
      ierr = 0
      !
      ALLOCATE( insitu_t(jpi,jpj,jpk), STAT=ierr(1) )
         !
      insitu_tem_alloc = MAXVAL(ierr)
      IF( lk_mpp )   CALL mpp_sum( insitu_tem_alloc )
      !
  END FUNCTION insitu_tem_alloc

!-----------------------------------------------------------------------------------
!
! Calculate the insitu temperature by integrating the adiabatic lapse rate from zero 
! to pressure at depth of tracer level. Based on UM subroutine POTTEM and function ATG
!
! Initial version: D. Acreman June 2006
!
!-----------------------------------------------------------------------------------

  SUBROUTINE theta2t()

    INTEGER, PARAMETER :: num_steps=10                ! number of steps in integration
    INTEGER            :: step                        ! iteration counter
    INTEGER            :: ji, jj, jk                  ! loop indices
    REAL(wp), DIMENSION(jpi,jpj,jpk) :: zP            ! pressure (decibars)
    REAL(wp), DIMENSION(jpi,jpj,jpk) :: zT            ! temperature at pressure p
    REAL(wp), DIMENSION(jpi,jpj,jpk) :: zTB           ! temperature at p-dp
    REAL(wp), DIMENSION(jpi,jpj,jpk) :: zTA           ! temperature at p+dp
    REAL(wp), DIMENSION(jpi,jpj,jpk) :: zDP           ! pressure step
    REAL(wp), DIMENSION(jpi,jpj,jpk) :: zSS           ! salinity(PSU) - 35.0
    REAL(wp), DIMENSION(jpi,jpj,jpk) :: zLAPSE        ! adiabatic lapse rate

!CDIR IEXPAND (ATG)

! Set the pressure interval for the integration. The integration is carried out from 
! zero (pressure at the surface) to the pressure at the depth of the tracer point. The
! pressure in decibars is represented by the depth in metres. Pressures are "Oceanographic"
! pressures equal to absolute pressure minus one atmosphere
     zDP(:,:,:) = 0.0
     DO jk = 1, jpkm1
        DO jj = 1, jpj
           DO ji = 1, jpi
              ! These loops expanded for case where fsdept may be 1D
              zDP(ji,jj,jk) = fsdept(ji,jj,jk) / real(num_steps)
           END DO
        END DO
     END DO

! Salinity at each point
     zSS(:,:,:)      = tsn(:,:,:,jp_sal) - 35.0 

! Set initial values of temperature and pressure. zT is the temperature at pressure zP, 
! zTB is the temperature at pressure zP-zdP and zTA is the temperature at pressure zP+zdP
     zT(:,:,:)  = tsn(:,:,:,jp_tem)
     zP(:,:,:)  = 0.0         ! Pressure at surface
     CALL ATG(zP, zT, zSS, zLAPSE)
     zTB(:,:,:) = zT(:,:,:)  - zLAPSE(:,:,:) * zDP(:,:,:)

     interation: DO step=1, num_steps
        ! Calculate lapse rate (dT/dP) and hence TA
        CALL ATG(zP, zT, zSS, zLAPSE)
        zTA(:,:,:) = zTB(:,:,:) + 2.0 * zLAPSE(:,:,:) * zDP(:,:,:)
        ! Have calculated TB, T and TA for this pressure, now advance solution by dP
        zP(:,:,:)  = zP(:,:,:)  + zDP
        zTB(:,:,:) = zT(:,:,:)
        zT(:,:,:)  = zTA(:,:,:)

     END DO interation

     insitu_t(:,:,:) = zT(:,:,:) * tmask(:,:,:)
     CALL lbc_lnk( insitu_t,  'T', 1.0)

   END SUBROUTINE theta2t

   SUBROUTINE ATG(P,T,DS,zLAPSE)

     REAL, INTENT(IN   ) :: P(jpi,jpj,jpk)        ! PRESSURE (DECIBARS)
     REAL, INTENT(IN   ) :: T(jpi,jpj,jpk)        ! TEMPERATURE (DEG C)
     REAL, INTENT(IN   ) :: DS(jpi,jpj,jpk)       ! SALINITY (PSU) -35.0
     REAL, INTENT(  OUT) :: zLAPSE(jpi,jpj,jpk)   ! LAPSE RATE

       zLAPSE = ((( -2.1687E-16*T+1.8676E-14)*T-4.6206E-13)*P &
           + ((2.7759E-12*T-1.1351E-10)*DS+((-5.4481E-14*T    &
           + 8.733E-12)*T-6.7795E-10)*T+1.8741E-8))*P         &
           + (-4.2393E-8*T+1.8932E-6)*DS                      &
           + ((6.6228E-10*T-6.836E-8)*T+8.5258E-6)*T+3.5803E-5 

   END SUBROUTINE ATG

END MODULE insitu_tem
