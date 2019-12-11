MODULE trcopt_medusa
   !!======================================================================
   !!                         ***  MODULE trcopt_medusa  ***
   !! TOP :   MEDUSA Compute the light availability in the water column
   !!======================================================================
   !! History :    -   !  1995-05  (M. Levy) Original code
   !!              -   !  1999-09  (J.-M. Andre, M. Levy) 
   !!              -   !  1999-11  (C. Menkes, M.-A. Foujols) itabe initial
   !!              -   !  2000-02  (M.A. Foujols) change x**y par exp(y*log(x))
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'                                      MEDUSA bio-model
   !!----------------------------------------------------------------------
   !!   trc_opt_medusa        :   Compute the light availability in the water column
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc
   USE prtctl_trc      ! Print control for debbuging
   USE sms_medusa

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_opt_medusa   ! called in trcprg.F90

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_opt_medusa( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_opt_medusa  ***
      !!
      !! ** Purpose :   computes the light propagation in the water column
      !!              and the euphotic layer depth
      !!
      !! ** Method  :   local par is computed in w layers using light propagation
      !!              mean par in t layers are computed by integration
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! index of the time stepping
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zpig                                    ! total pigment
      REAL(wp) ::   zkr                                     ! total absorption coefficient in red
      REAL(wp) ::   zkg                                     ! total absorption coefficient in green
      REAL(wp) ::   totchl                                  ! total Chl concentreation
      REAL(wp), DIMENSION(jpi,jpj)     ::   zpar100         ! irradiance at euphotic layer depth
      REAL(wp), DIMENSION(jpi,jpj)     ::   zpar0m          ! irradiance just below the surface
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zparr, zparg    ! red and green compound of par

      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      !! AXY (20/11/14): alter this to report on first MEDUSA call
      !! IF( kt == nit000 ) THEN
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_opt_medusa: MEDUSA optic-model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
	 IF(lwp) WRITE(numout,*) ' kt =',kt
      ENDIF

      ! determination of surface irradiance
      ! -----------------------------------
      ! AXY (23/07/15): the inclusion of empirical DMS calculations requires
      !                 daily averages of a series of properties that are
      !                 used as inputs; these include surface irradiance; 
      !                 here, this is taken advantage of to allow MEDUSA to
      !                 base its submarine light field on daily average
      !                 rather than "instantaneous" irradiance; largely
      !                 because MEDUSA was originally formulated to work
      !                 with diel average irradiance rather than a diel
      !                 cycle; using key_avgqsr_medusa activates this
      !                 functionality, while its absence gives default
      !                 MEDUSA (which is whatever is supplied by NEMO)
# if defined key_avgqsr_medusa
      ! surface irradiance input is rolling average irradiance
      zpar0m (:,:)   = zn_dms_qsr(:,:) * 0.43
# else      
      ! surface irradiance input is   instantaneous irradiance
      zpar0m (:,:)   =        qsr(:,:) * 0.43
# endif
      ! AXY (22/08/14): when zpar0m = 0, zpar100 is also zero and calculating 
      !                 euphotic depth is not possible (cf. the Arctic Octopus); 
      !                 a "solution" to this is to set zpar0m to some minimal
      !                 value such that zpar100 also has a non-zero value and
      !                 euphotic depth can be calculated properly; note that,
      !                 in older, non-diurnal versions of NEMO, this was much
      !                 less of a problem; note also that, if pushed, I will
      !                 claim that my minimal value of zpar0m refers to light
      !                 from stars
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( zpar0m(ji,jj) <= 0.0 ) zpar0m(ji,jj) = 0.001  ! = 1 mW/m2
         ENDDO
      ENDDO
      zpar100(:,:)   = zpar0m(:,:) * 0.01
      xpar   (:,:,1) = zpar0m(:,:)
      zparr  (:,:,1) = 0.5 * zpar0m(:,:)
      zparg  (:,:,1) = 0.5 * zpar0m(:,:)

      ! determination of xpar
      ! ---------------------

      DO jk = 2, jpk                     ! determination of local par in w levels
         DO jj = 1, jpj
            DO ji = 1, jpi
               totchl =trn(ji,jj,jk-1,jpchn)+trn(ji,jj,jk-1,jpchd)
               zpig = MAX( TINY(0.), totchl/rpig) 
               zkr  = xkr0 + xkrp * EXP( xlr * LOG( zpig ) )
               zkg  = xkg0 + xkgp * EXP( xlg * LOG( zpig ) )
               zparr(ji,jj,jk) = zparr(ji,jj,jk-1) * EXP( -zkr * fse3t(ji,jj,jk-1) )
               zparg(ji,jj,jk) = zparg(ji,jj,jk-1) * EXP( -zkg * fse3t(ji,jj,jk-1) )
            END DO
        END DO
      END DO

      DO jk = 1, jpkm1                   ! mean par in t levels
         DO jj = 1, jpj
            DO ji = 1, jpi
               totchl =trn(ji,jj,jk  ,jpchn)+trn(ji,jj,jk  ,jpchd)
               zpig = MAX( TINY(0.), totchl/rpig) 
               zkr  = xkr0 + xkrp * EXP( xlr * LOG( zpig ) )
               zkg  = xkg0 + xkgp * EXP( xlg * LOG( zpig ) )
               zparr(ji,jj,jk)    = zparr(ji,jj,jk) / zkr / fse3t(ji,jj,jk) * ( 1 - EXP( -zkr*fse3t(ji,jj,jk) ) )
               zparg(ji,jj,jk)    = zparg(ji,jj,jk) / zkg / fse3t(ji,jj,jk) * ( 1 - EXP( -zkg*fse3t(ji,jj,jk) ) )
               xpar (ji,jj,jk) = MAX( zparr(ji,jj,jk) + zparg(ji,jj,jk), 1.e-15 )
            END DO
         END DO
      END DO

      ! 3. Determination of euphotic layer depth
      ! ----------------------------------------

      ! Euphotic layer bottom level
      neln(:,:) = 1                                           ! initialisation of EL level
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji = 1, jpi
              IF( xpar(ji,jj,jk) >= zpar100(ji,jj) )   neln(ji,jj) = jk+1 ! 1rst T-level strictly below EL bottom
              !                                                  ! nb. this is to ensure compatibility with
              !                                                  ! nmld_trc definition in trd_mld_trc_zint
           END DO
         END DO
      ENDDO

      ! Euphotic layer depth
      !! Jpalm -- 06-03-2017 -- add init xze, to avoid halo problems within the
      !!                        writing process
      xze(:,:) = 0.0
      DO jj = 1, jpj
         DO ji = 1, jpi
            xze(ji,jj) = fsdepw( ji, jj, neln(ji,jj) )            ! exact EL depth
         END DO
      ENDDO 

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('opt')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF

   END SUBROUTINE trc_opt_medusa

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_opt_medusa( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_opt_medusa: You should not have seen this print! error?', kt
   END SUBROUTINE trc_opt_medusa
#endif 

   !!======================================================================
END MODULE trcopt_medusa
