MODULE zdfbfr
   !!======================================================================
   !!                       ***  MODULE  zdfbfr  ***
   !! Ocean physics: Bottom friction
   !!======================================================================
   !! History :  OPA  ! 1997-06  (G. Madec, A.-M. Treguier)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-09  (A.C.Coward)  Correction to include barotropic contribution
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.4  ! 2011-11  (H. Liu) implementation of semi-implicit bottom friction option
   !!                 ! 2012-06  (H. Liu) implementation of Log Layer bottom friction option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_bfr      : update bottom friction coefficient (non-linear bottom friction only)
   !!   zdf_bfr_init : read in namelist and control the bottom friction parameters.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics variables
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE prtctl          ! Print control
   USE timing          ! Timing
   USE wrk_nemo        ! Memory Allocation
   USE phycst, ONLY: vkarmn

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_bfr         ! called by step.F90
   PUBLIC   zdf_bfr_init    ! called by nemogcm.F90

   !                                 !!* Namelist nambfr: bottom friction namelist *
   INTEGER , PUBLIC ::   nn_bfr       ! = 0/1/2/3 type of bottom friction  (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_bfri1     ! bottom drag coefficient (linear case)  (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_bfri2     ! bottom drag coefficient (non linear case) (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_bfri2_max ! Maximum bottom drag coefficient (non linear case and ln_loglayer=T) (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_bfeb2     ! background bottom turbulent kinetic energy  [m2/s2] (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_bfrien    ! local factor to enhance coefficient bfri (PUBLIC for TAM)
   LOGICAL , PUBLIC ::   ln_bfr2d     ! logical switch for 2D enhancement (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_tfri1     ! top drag coefficient (linear case)  (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_tfri2     ! top drag coefficient (non linear case) (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_tfri2_max ! Maximum top drag coefficient (non linear case and ln_loglayer=T) (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_tfeb2     ! background top turbulent kinetic energy  [m2/s2] (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_tfrien    ! local factor to enhance coefficient tfri (PUBLIC for TAM)
   LOGICAL , PUBLIC ::   ln_tfr2d     ! logical switch for 2D enhancement (PUBLIC for TAM)

   LOGICAL , PUBLIC ::   ln_loglayer  ! switch for log layer bfr coeff. (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_bfrz0     ! bottom roughness for loglayer bfr coeff (PUBLIC for TAM)
   REAL(wp), PUBLIC ::   rn_tfrz0     ! bottom roughness for loglayer bfr coeff (PUBLIC for TAM)
   LOGICAL , PUBLIC ::   ln_bfrimp    ! logical switch for implicit bottom friction
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::  bfrcoef2d, tfrcoef2d   ! 2D bottom/top drag coefficient (PUBLIC for TAM)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_bfr_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_bfr_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( bfrcoef2d(jpi,jpj), tfrcoef2d(jpi,jpj), STAT=zdf_bfr_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( zdf_bfr_alloc )
      IF( zdf_bfr_alloc /= 0 )   CALL ctl_warn('zdf_bfr_alloc: failed to allocate arrays.')
   END FUNCTION zdf_bfr_alloc


   SUBROUTINE zdf_bfr( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_bfr  ***
      !!
      !! ** Purpose :   compute the bottom friction coefficient.
      !!
      !! ** Method  :   Calculate and store part of the momentum trend due
      !!              to bottom friction following the chosen friction type
      !!              (free-slip, linear, or quadratic). The component
      !!              calculated here is multiplied by the bottom velocity in
      !!              dyn_bfr to provide the trend term.
      !!                The coefficients are updated at each time step only
      !!              in the quadratic case.
      !!
      !! ** Action  :   bfrua , bfrva   bottom friction coefficients
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj                       ! dummy loop indices
      INTEGER  ::   ikbt, ikbu, ikbv             ! local integers
      REAL(wp) ::   zvu, zuv, zecu, zecv, ztmp   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::  zbfrt, ztfrt
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_bfr')
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_bfr : Set bottom friction coefficient (non-linear case)'
         WRITE(numout,*) '~~~~~~~~'
      ENDIF
      !
      IF( nn_bfr == 2 ) THEN                 ! quadratic bottom friction only
         !
         CALL wrk_alloc( jpi, jpj, zbfrt, ztfrt )

         IF ( ln_loglayer.AND.lk_vvl ) THEN ! "log layer" bottom friction coefficient

            DO jj = 1, jpj
               DO ji = 1, jpi
                  ikbt = mbkt(ji,jj)
!! JC: possible WAD implementation should modify line below if layers vanish
                  ztmp = tmask(ji,jj,ikbt) * ( vkarmn / LOG( 0.5_wp * fse3t_n(ji,jj,ikbt) / rn_bfrz0 ))**2._wp
                  zbfrt(ji,jj) = MAX(bfrcoef2d(ji,jj), ztmp)
                  zbfrt(ji,jj) = MIN(zbfrt(ji,jj), rn_bfri2_max)
               END DO
            END DO
! (ISF)
            IF ( ln_isfcav ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     ikbt = mikt(ji,jj)
! JC: possible WAD implementation should modify line below if layers vanish
                     ztmp = (1-tmask(ji,jj,1)) * ( vkarmn / LOG( 0.5_wp * fse3t_n(ji,jj,ikbt) / rn_bfrz0 ))**2._wp
                     ztfrt(ji,jj) = MAX(tfrcoef2d(ji,jj), ztmp)
                     ztfrt(ji,jj) = MIN(ztfrt(ji,jj), rn_tfri2_max)
                  END DO
               END DO
            END IF
         !   
         ELSE
            zbfrt(:,:) = bfrcoef2d(:,:)
            ztfrt(:,:) = tfrcoef2d(:,:)
         ENDIF

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points
               ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
               !
               zvu  = 0.25 * (  vn(ji,jj  ,ikbu) + vn(ji+1,jj  ,ikbu)     &
                  &           + vn(ji,jj-1,ikbu) + vn(ji+1,jj-1,ikbu)  )
               zuv  = 0.25 * (  un(ji,jj  ,ikbv) + un(ji-1,jj  ,ikbv)     &
                  &           + un(ji,jj+1,ikbv) + un(ji-1,jj+1,ikbv)  )
               !
               zecu = SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_bfeb2  )
               zecv = SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_bfeb2  )
               !
               bfrua(ji,jj) = - 0.5_wp * ( zbfrt(ji,jj) + zbfrt(ji+1,jj  ) ) * zecu
               bfrva(ji,jj) = - 0.5_wp * ( zbfrt(ji,jj) + zbfrt(ji  ,jj+1) ) * zecv
               !
               ! in case of 2 cell water column, we assume each cell feels the top and bottom friction
               IF ( ln_isfcav ) THEN
                  IF ( miku(ji,jj) + 1 .GE. mbku(ji,jj) ) THEN
                     bfrua(ji,jj) = - 0.5_wp * ( ( zbfrt(ji,jj) + zbfrt(ji+1,jj  ) )   &
                                  &            + ( ztfrt(ji,jj) + ztfrt(ji+1,jj  ) ) ) &
                                  &          * zecu * (1._wp - umask(ji,jj,1))
                  END IF
                  IF ( mikv(ji,jj) + 1 .GE. mbkv(ji,jj) ) THEN
                     bfrva(ji,jj) = - 0.5_wp * ( ( zbfrt(ji,jj) + zbfrt(ji  ,jj+1) )   &
                                  &            + ( ztfrt(ji,jj) + ztfrt(ji  ,jj+1) ) ) &
                                  &          * zecv * (1._wp - vmask(ji,jj,1))
                  END IF
               END IF
            END DO
         END DO
         CALL lbc_lnk( bfrua, 'U', 1. )   ;   CALL lbc_lnk( bfrva, 'V', 1. )      ! Lateral boundary condition

         IF ( ln_isfcav ) THEN
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ! (ISF) ========================================================================
                  ikbu = miku(ji,jj)         ! ocean top level at u- and v-points 
                  ikbv = mikv(ji,jj)         ! (1st wet ocean u- and v-points)
                  !
                  zvu  = 0.25 * (  vn(ji,jj  ,ikbu) + vn(ji+1,jj  ,ikbu)     &
                     &           + vn(ji,jj-1,ikbu) + vn(ji+1,jj-1,ikbu)  )
                  zuv  = 0.25 * (  un(ji,jj  ,ikbv) + un(ji-1,jj  ,ikbv)     &
                     &           + un(ji,jj+1,ikbv) + un(ji-1,jj+1,ikbv)  )
              !
                  zecu = SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_tfeb2 )
                  zecv = SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_tfeb2 )
              !
                  tfrua(ji,jj) = - 0.5_wp * ( ztfrt(ji,jj) + ztfrt(ji+1,jj  ) ) * zecu * (1._wp - umask(ji,jj,1))
                  tfrva(ji,jj) = - 0.5_wp * ( ztfrt(ji,jj) + ztfrt(ji  ,jj+1) ) * zecv * (1._wp - vmask(ji,jj,1))
              ! (ISF) END ====================================================================
              ! in case of 2 cell water column, we assume each cell feels the top and bottom friction
                  IF ( miku(ji,jj) + 1 .GE. mbku(ji,jj) ) THEN
                     tfrua(ji,jj) = - 0.5_wp * ( ( ztfrt(ji,jj) + ztfrt(ji+1,jj  ) )   &
                                  &            + ( zbfrt(ji,jj) + zbfrt(ji+1,jj  ) ) ) &
                                  &          * zecu * (1._wp - umask(ji,jj,1))
                  END IF
                  IF ( mikv(ji,jj) + 1 .GE. mbkv(ji,jj) ) THEN
                     tfrva(ji,jj) = - 0.5_wp * ( ( ztfrt(ji,jj) + ztfrt(ji  ,jj+1) )   &
                                  &            + ( zbfrt(ji,jj) + zbfrt(ji  ,jj+1) ) ) &
                                  &          * zecv * (1._wp - vmask(ji,jj,1))
                  END IF
               END DO
            END DO
            CALL lbc_lnk( tfrua, 'U', 1. )   ;   CALL lbc_lnk( tfrva, 'V', 1. )      ! Lateral boundary condition
         END IF
         !
         !
         IF(ln_ctl)   CALL prt_ctl( tab2d_1=bfrua, clinfo1=' bfr  - u: ', mask1=umask,        &
            &                       tab2d_2=bfrva, clinfo2=       ' v: ', mask2=vmask,ovlap=1 )
         CALL wrk_dealloc( jpi,jpj,zbfrt, ztfrt )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_bfr')
      !
   END SUBROUTINE zdf_bfr


   SUBROUTINE zdf_bfr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_bfr_init  ***
      !!
      !! ** Purpose :   Initialization of the bottom friction
      !!
      !! ** Method  :   Read the nambfr namelist and check their consistency
      !!                called at the first timestep (nit000)
      !!----------------------------------------------------------------------
      USE iom   ! I/O module for ehanced bottom friction file
      !!
      INTEGER   ::   inum         ! logical unit for enhanced bottom friction file
      INTEGER   ::   ji, jj       ! dummy loop indexes
      INTEGER   ::   ikbt, ikbu, ikbv   ! temporary integers
      INTEGER   ::   ictu, ictv         !    -          -
      INTEGER   ::   ios
      REAL(wp)  ::   zminbfr, zmaxbfr   ! temporary scalars
      REAL(wp)  ::   zmintfr, zmaxtfr   ! temporary scalars
      REAL(wp)  ::   ztmp, zfru, zfrv   !    -         -
      !!
      NAMELIST/nambfr/ nn_bfr, rn_bfri1, rn_bfri2, rn_bfri2_max, rn_bfeb2, rn_bfrz0, ln_bfr2d, &
                    &          rn_tfri1, rn_tfri2, rn_tfri2_max, rn_tfeb2, rn_tfrz0, ln_tfr2d, &
                    &  rn_bfrien, rn_tfrien, ln_bfrimp, ln_loglayer
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_bfr_init')
      !
      !                              !* Allocate zdfbfr arrays
      IF( zdf_bfr_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_bfr_init : unable to allocate arrays' )
      !
      !                              !* Parameter control and print
      !
      REWIND( numnam_ref )              ! Namelist nambfr in reference namelist : Bottom momentum boundary condition
      READ  ( numnam_ref, nambfr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambfr in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nambfr in configuration namelist : Bottom momentum boundary condition
      READ  ( numnam_cfg, nambfr, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambfr in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nambfr )
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zdf_bfr_init : momentum bottom friction'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '   Namelist nam_bfr : set bottom friction parameters'
      !
      SELECT CASE (nn_bfr)
      !
      CASE( 0 )
         IF(lwp) WRITE(numout,*) '      free-slip '
         bfrua(:,:) = 0.e0
         bfrva(:,:) = 0.e0
         tfrua(:,:) = 0.e0
         tfrva(:,:) = 0.e0
         !
      CASE( 1 )
         IF(lwp) WRITE(numout,*) '      linear botton friction'
         IF(lwp) WRITE(numout,*) '      bottom friction coef.   rn_bfri1  = ', rn_bfri1
         IF( ln_bfr2d ) THEN
            IF(lwp) WRITE(numout,*) '      read coef. enhancement distribution from file   ln_bfr2d  = ', ln_bfr2d
            IF(lwp) WRITE(numout,*) '      coef rn_bfri2 enhancement factor                rn_bfrien  = ',rn_bfrien
         ENDIF
         IF ( ln_isfcav ) THEN
            IF(lwp) WRITE(numout,*) '      top    friction coef.   rn_bfri1  = ', rn_tfri1
            IF( ln_tfr2d ) THEN
               IF(lwp) WRITE(numout,*) '      read coef. enhancement distribution from file   ln_tfr2d  = ', ln_tfr2d
               IF(lwp) WRITE(numout,*) '      coef rn_tfri2 enhancement factor                rn_tfrien  = ',rn_tfrien
            ENDIF
         END IF
         !
         IF(ln_bfr2d) THEN
            ! bfr_coef is a coefficient in [0,1] giving the mask where to apply the bfr enhancement
            CALL iom_open('bfr_coef.nc',inum)
            CALL iom_get (inum, jpdom_data, 'bfr_coef',bfrcoef2d,1) ! bfrcoef2d is used as tmp array
            CALL iom_close(inum)
            bfrcoef2d(:,:) = rn_bfri1 * ( 1 + rn_bfrien * bfrcoef2d(:,:) )
         ELSE
            bfrcoef2d(:,:) = rn_bfri1  ! initialize bfrcoef2d to the namelist variable
         ENDIF
         !
         bfrua(:,:) = - bfrcoef2d(:,:)
         bfrva(:,:) = - bfrcoef2d(:,:)
         !
         IF ( ln_isfcav ) THEN
            IF(ln_tfr2d) THEN
               ! tfr_coef is a coefficient in [0,1] giving the mask where to apply the bfr enhancement
               CALL iom_open('tfr_coef.nc',inum)
               CALL iom_get (inum, jpdom_data, 'tfr_coef',tfrcoef2d,1) ! tfrcoef2d is used as tmp array
               CALL iom_close(inum)
               tfrcoef2d(:,:) = rn_tfri1 * ( 1 + rn_tfrien * tfrcoef2d(:,:) )
            ELSE
               tfrcoef2d(:,:) = rn_tfri1  ! initialize tfrcoef2d to the namelist variable
            ENDIF
            !
            tfrua(:,:) = - tfrcoef2d(:,:)
            tfrva(:,:) = - tfrcoef2d(:,:)
         END IF
         !
      CASE( 2 )
         IF(lwp) WRITE(numout,*) '      quadratic bottom friction'
         IF(lwp) WRITE(numout,*) '      friction coef.   rn_bfri2  = ', rn_bfri2
         IF(lwp) WRITE(numout,*) '      Max. coef. (log case)   rn_bfri2_max  = ', rn_bfri2_max
         IF(lwp) WRITE(numout,*) '      background tke   rn_bfeb2  = ', rn_bfeb2
         IF(lwp) WRITE(numout,*) '      log formulation   ln_bfr2d = ', ln_loglayer
         IF(lwp) WRITE(numout,*) '      bottom roughness  rn_bfrz0 [m] = ', rn_bfrz0
         IF( rn_bfrz0<=0.e0 ) THEN
            WRITE(ctmp1,*) '      bottom roughness must be strictly positive'
            CALL ctl_stop( ctmp1 )
         ENDIF
         IF( ln_bfr2d ) THEN
            IF(lwp) WRITE(numout,*) '      read coef. enhancement distribution from file   ln_bfr2d  = ', ln_bfr2d
            IF(lwp) WRITE(numout,*) '      coef rn_bfri2 enhancement factor                rn_bfrien  = ',rn_bfrien
         ENDIF
         IF ( ln_isfcav ) THEN
            IF(lwp) WRITE(numout,*) '      quadratic top    friction'
            IF(lwp) WRITE(numout,*) '      friction coef.    rn_tfri2     = ', rn_tfri2
            IF(lwp) WRITE(numout,*) '      Max. coef. (log case)   rn_tfri2_max  = ', rn_tfri2_max
            IF(lwp) WRITE(numout,*) '      background tke    rn_tfeb2     = ', rn_tfeb2
            IF(lwp) WRITE(numout,*) '      log formulation   ln_tfr2d     = ', ln_loglayer
            IF(lwp) WRITE(numout,*) '      top roughness     rn_tfrz0 [m] = ', rn_tfrz0
            IF( rn_tfrz0<=0.e0 ) THEN
               WRITE(ctmp1,*) '      top roughness must be strictly positive'
               CALL ctl_stop( ctmp1 )
            ENDIF
            IF( ln_tfr2d ) THEN
               IF(lwp) WRITE(numout,*) '      read coef. enhancement distribution from file   ln_tfr2d  = ', ln_tfr2d
               IF(lwp) WRITE(numout,*) '      coef rn_tfri2 enhancement factor                rn_tfrien  = ',rn_tfrien
            ENDIF
         END IF
         !
         IF(ln_bfr2d) THEN
            ! bfr_coef is a coefficient in [0,1] giving the mask where to apply the bfr enhancement
            CALL iom_open('bfr_coef.nc',inum)
            CALL iom_get (inum, jpdom_data, 'bfr_coef',bfrcoef2d,1) ! bfrcoef2d is used as tmp array
            CALL iom_close(inum)
            !
            bfrcoef2d(:,:) = rn_bfri2 * ( 1 + rn_bfrien * bfrcoef2d(:,:) )
         ELSE
            bfrcoef2d(:,:) = rn_bfri2  ! initialize bfrcoef2d to the namelist variable
         ENDIF
         
         IF ( ln_isfcav ) THEN
            IF(ln_tfr2d) THEN
               ! tfr_coef is a coefficient in [0,1] giving the mask where to apply the bfr enhancement
               CALL iom_open('tfr_coef.nc',inum)
               CALL iom_get (inum, jpdom_data, 'tfr_coef',tfrcoef2d,1) ! tfrcoef2d is used as tmp array
               CALL iom_close(inum)
               !
               tfrcoef2d(:,:) = rn_tfri2 * ( 1 + rn_tfrien * tfrcoef2d(:,:) )
            ELSE
               tfrcoef2d(:,:) = rn_tfri2  ! initialize tfrcoef2d to the namelist variable
            ENDIF
         END IF
         !
         IF ( ln_loglayer.AND.(.NOT.lk_vvl) ) THEN ! set "log layer" bottom friction once for all
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ikbt = mbkt(ji,jj)
                  ztmp = tmask(ji,jj,ikbt) * ( vkarmn / LOG( 0.5_wp * fse3t_n(ji,jj,ikbt) / rn_bfrz0 ))**2._wp
                  bfrcoef2d(ji,jj) = MAX(bfrcoef2d(ji,jj), ztmp)
                  bfrcoef2d(ji,jj) = MIN(bfrcoef2d(ji,jj), rn_bfri2_max)
               END DO
            END DO
            IF ( ln_isfcav ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     ikbt = mikt(ji,jj)
                     ztmp = tmask(ji,jj,ikbt) * ( vkarmn / LOG( 0.5_wp * fse3t_n(ji,jj,ikbt) / rn_tfrz0 ))**2._wp
                     tfrcoef2d(ji,jj) = MAX(tfrcoef2d(ji,jj), ztmp)
                     tfrcoef2d(ji,jj) = MIN(tfrcoef2d(ji,jj), rn_tfri2_max)
                  END DO
               END DO
            END IF
         ENDIF
         !
      CASE DEFAULT
         IF(lwp) WRITE(ctmp1,*) '         bad flag value for nn_bfr = ', nn_bfr
         CALL ctl_stop( ctmp1 )
         !
      END SELECT
      !
      IF(lwp) WRITE(numout,*) '      implicit bottom friction switch                ln_bfrimp  = ', ln_bfrimp
      !
      !                              ! Make sure ln_zdfexp=.false. when use implicit bfr
      IF( ln_bfrimp .AND. ln_zdfexp ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'Implicit bottom friction can only be used when ln_zdfexp=.false.'
            WRITE(numout,*) '         but you set: ln_bfrimp=.true. and ln_zdfexp=.true.!!!!'
            WRITE(ctmp1,*)  '         set either ln_zdfexp = .false or ln_bfrimp = .false.'
            CALL ctl_stop( ctmp1 )
         END IF
      END IF
      !
      ! Basic stability check on bottom friction coefficient
      !
      ictu = 0               ! counter for stability criterion breaches at U-pts
      ictv = 0               ! counter for stability criterion breaches at V-pts
      zminbfr =  1.e10_wp    ! initialise tracker for minimum of bottom friction coefficient
      zmaxbfr = -1.e10_wp    ! initialise tracker for maximum of bottom friction coefficient
      zmintfr =  1.e10_wp    ! initialise tracker for minimum of bottom friction coefficient
      zmaxtfr = -1.e10_wp    ! initialise tracker for maximum of bottom friction coefficient
      !
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
             ikbu = mbku(ji,jj)       ! deepest ocean level at u- and v-points
             ikbv = mbkv(ji,jj)
             zfru = 0.5 * fse3u(ji,jj,ikbu) / rdt
             zfrv = 0.5 * fse3v(ji,jj,ikbv) / rdt
             IF( ABS( bfrcoef2d(ji,jj) ) > zfru ) THEN
                IF( ln_ctl ) THEN
                   WRITE(numout,*) 'BFR ', narea, nimpp+ji, njmpp+jj, ikbu
                   WRITE(numout,*) 'BFR ', ABS( bfrcoef2d(ji,jj) ), zfru
                ENDIF
                ictu = ictu + 1
             ENDIF
             IF( ABS( bfrcoef2d(ji,jj) ) > zfrv ) THEN
                 IF( ln_ctl ) THEN
                     WRITE(numout,*) 'BFR ', narea, nimpp+ji, njmpp+jj, ikbv
                     WRITE(numout,*) 'BFR ', bfrcoef2d(ji,jj), zfrv
                 ENDIF
                 ictv = ictv + 1
             ENDIF
             zminbfr = MIN(  zminbfr, MIN( zfru, ABS( bfrcoef2d(ji,jj) ) )  )
             zmaxbfr = MAX(  zmaxbfr, MIN( zfrv, ABS( bfrcoef2d(ji,jj) ) )  )
! (ISF)
             IF ( ln_isfcav ) THEN
                ikbu = miku(ji,jj)       ! 1st wet ocean level at u- and v-points
                ikbv = mikv(ji,jj)
                zfru = 0.5 * fse3u(ji,jj,ikbu) / rdt
                zfrv = 0.5 * fse3v(ji,jj,ikbv) / rdt
                IF( ABS( tfrcoef2d(ji,jj) ) > zfru ) THEN
                   IF( ln_ctl ) THEN
                      WRITE(numout,*) 'TFR ', narea, nimpp+ji, njmpp+jj, ikbu
                      WRITE(numout,*) 'TFR ', ABS( tfrcoef2d(ji,jj) ), zfru
                   ENDIF
                   ictu = ictu + 1
                ENDIF
                IF( ABS( tfrcoef2d(ji,jj) ) > zfrv ) THEN
                   IF( ln_ctl ) THEN
                      WRITE(numout,*) 'TFR ', narea, nimpp+ji, njmpp+jj, ikbv
                      WRITE(numout,*) 'TFR ', tfrcoef2d(ji,jj), zfrv
                   ENDIF
                   ictv = ictv + 1
                ENDIF
                zmintfr = MIN(  zmintfr, MIN( zfru, ABS( tfrcoef2d(ji,jj) ) )  )
                zmaxtfr = MAX(  zmaxtfr, MIN( zfrv, ABS( tfrcoef2d(ji,jj) ) )  )
             END IF
! END ISF
         END DO
      END DO
      IF( lk_mpp ) THEN
         CALL mpp_sum( ictu )
         CALL mpp_sum( ictv )
         CALL mpp_min( zminbfr )
         CALL mpp_max( zmaxbfr )
         IF ( ln_isfcav) CALL mpp_min( zmintfr )
         IF ( ln_isfcav) CALL mpp_max( zmaxtfr )
      ENDIF
      IF( .NOT.ln_bfrimp) THEN
      IF( lwp .AND. ictu + ictv > 0 ) THEN
         WRITE(numout,*) ' Bottom/Top friction stability check failed at ', ictu, ' U-points '
         WRITE(numout,*) ' Bottom/Top friction stability check failed at ', ictv, ' V-points '
         WRITE(numout,*) ' Bottom friction coefficient now ranges from: ', zminbfr, ' to ', zmaxbfr
         IF ( ln_isfcav ) WRITE(numout,*) ' Top friction coefficient now ranges from: ', zmintfr, ' to ', zmaxtfr
         WRITE(numout,*) ' Bottom/Top friction coefficient will be reduced where necessary'
      ENDIF
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_bfr_init')
      !
   END SUBROUTINE zdf_bfr_init

   !!======================================================================
END MODULE zdfbfr
