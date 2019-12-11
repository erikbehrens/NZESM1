MODULE sbctide
   !!======================================================================
   !!                       ***  MODULE  sbctide  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  9.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE phycst
   USE daymod
   USE dynspg_oce
   USE tideini
   !
   USE iom
   USE in_out_manager  ! I/O units
   USE ioipsl          ! NetCDF IPSL library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PUBLIC

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   pot_astro   !

#if defined key_tide
   !!----------------------------------------------------------------------
   !!   'key_tide' :                                        tidal potential
   !!----------------------------------------------------------------------
   !!   sbc_tide            : 
   !!   tide_init_potential :
   !!----------------------------------------------------------------------

   LOGICAL, PUBLIC, PARAMETER ::   lk_tide  = .TRUE.
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   amp_pot, phi_pot

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.5 , NEMO Consortium (2013)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_tide( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE sbc_tide  ***
      !!----------------------------------------------------------------------      
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step
      INTEGER               ::   jk     ! dummy loop index
      !!----------------------------------------------------------------------

      IF( nsec_day == NINT(0.5_wp * rdttra(1)) ) THEN      ! start a new day
         !
         IF( kt == nit000 ) THEN
            ALLOCATE( amp_pot(jpi,jpj,nb_harmo),                      &
               &      phi_pot(jpi,jpj,nb_harmo), pot_astro(jpi,jpj)   )
         ENDIF
         !
         amp_pot(:,:,:) = 0._wp
         phi_pot(:,:,:) = 0._wp
         pot_astro(:,:) = 0._wp
         !
         CALL tide_harmo( omega_tide, v0tide, utide, ftide, ntide, nb_harmo )
         !
         kt_tide = kt
         !
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_tide : Update of the components and (re)Init. the potential at kt=', kt
            WRITE(numout,*) '~~~~~~~~ '
            DO jk = 1, nb_harmo
               WRITE(numout,*) Wave(ntide(jk))%cname_tide, utide(jk), ftide(jk), v0tide(jk), omega_tide(jk)
            END DO
         ENDIF
         !
         IF( ln_tide_pot )   CALL tide_init_potential
         !
      ENDIF
      !
   END SUBROUTINE sbc_tide


   SUBROUTINE tide_init_potential
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_potential  ***
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcons, ztmp1, ztmp2, zlat, zlon, ztmp, zamp, zcs   ! local scalar
      !!----------------------------------------------------------------------

      DO jk = 1, nb_harmo
         zcons = 0.7_wp * Wave(ntide(jk))%equitide * ftide(jk)
         DO ji = 1, jpi
            DO jj = 1, jpj
               ztmp1 =  amp_pot(ji,jj,jk) * COS( phi_pot(ji,jj,jk) )
               ztmp2 = -amp_pot(ji,jj,jk) * SIN( phi_pot(ji,jj,jk) )
               zlat = gphit(ji,jj)*rad !! latitude en radian
               zlon = glamt(ji,jj)*rad !! longitude en radian
               ztmp = v0tide(jk) + utide(jk) + Wave(ntide(jk))%nutide * zlon
               ! le potentiel est composé des effets des astres:
               IF    ( Wave(ntide(jk))%nutide == 1 )  THEN  ;  zcs = zcons * SIN( 2._wp*zlat )
               ELSEIF( Wave(ntide(jk))%nutide == 2 )  THEN  ;  zcs = zcons * COS( zlat )**2
               ELSE                                         ;  zcs = 0._wp
               ENDIF
               ztmp1 = ztmp1 + zcs * COS( ztmp )
               ztmp2 = ztmp2 - zcs * SIN( ztmp )
               zamp = SQRT( ztmp1*ztmp1 + ztmp2*ztmp2 )
               amp_pot(ji,jj,jk) = zamp
               phi_pot(ji,jj,jk) = ATAN2( -ztmp2 / MAX( 1.e-10_wp , zamp ) ,   &
                  &                        ztmp1 / MAX( 1.e-10_wp,  zamp )   )
            END DO
         END DO
      END DO
      !
   END SUBROUTINE tide_init_potential

#else
  !!----------------------------------------------------------------------
  !!   Default case :   Empty module
  !!----------------------------------------------------------------------
  LOGICAL, PUBLIC, PARAMETER ::   lk_tide = .FALSE.
CONTAINS
  SUBROUTINE sbc_tide( kt )      ! Empty routine
    INTEGER         , INTENT(in) ::   kt         ! ocean time-step
    WRITE(*,*) 'sbc_tide: You should not have seen this print! error?', kt
  END SUBROUTINE sbc_tide
#endif

  !!======================================================================
END MODULE sbctide
