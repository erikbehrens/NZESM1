MODULE dynspg_ts
   !!======================================================================
   !! History :   1.0  ! 2004-12  (L. Bessieres, G. Madec)  Original code
   !!              -   ! 2005-11  (V. Garnier, G. Madec)  optimization
   !!              -   ! 2006-08  (S. Masson)  distributed restart using iom
   !!             2.0  ! 2007-07  (D. Storkey) calls to BDY routines
   !!              -   ! 2008-01  (R. Benshila)  change averaging method
   !!             3.2  ! 2009-07  (R. Benshila, G. Madec) Complete revisit associated to vvl reactivation
   !!             3.3  ! 2010-09  (D. Storkey, E. O'Dea) update for BDY for Shelf configurations
   !!             3.3  ! 2011-03  (R. Benshila, R. Hordoir, P. Oddo) update calculation of ub_b
   !!             3.5  ! 2013-07  (J. Chanut) Switch to Forward-backward time stepping
   !!             3.6  ! 2013-11  (A. Coward) Update for z-tilde compatibility
   !!---------------------------------------------------------------------
#if defined key_dynspg_ts   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dynspg_ts'         split explicit free surface
   !!----------------------------------------------------------------------
   !!   dyn_spg_ts  : compute surface pressure gradient trend using a time-
   !!                 splitting scheme and add to the general trend 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcisf          ! ice shelf variable (fwfisf)
   USE dynspg_oce      ! surface pressure gradient variables
   USE phycst          ! physical constants
   USE dynvor          ! vorticity term
   USE bdy_par         ! for lk_bdy
   USE bdytides        ! open boundary condition data     
   USE bdydyn2d        ! open boundary conditions on barotropic variables
   USE sbctide         ! tides
   USE updtide         ! tide potential
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE in_out_manager  ! I/O manager
   USE iom             ! IOM library
   USE restart         ! only for lrst_oce
   USE zdf_oce         ! Bottom friction coefts
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing    
   USE sbcapr          ! surface boundary condition: atmospheric pressure
   USE dynadv, ONLY: ln_dynadv_vec
#if defined key_agrif
   USE agrif_opa_interp ! agrif
#endif
#if defined key_asminc   
   USE asminc          ! Assimilation increment
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_spg_ts        ! routine called in dynspg.F90 
   PUBLIC dyn_spg_ts_alloc  !    "      "     "    "
   PUBLIC dyn_spg_ts_init   !    "      "     "    "
   PUBLIC ts_rst            !    "      "     "    "

   INTEGER, SAVE :: icycle  ! Number of barotropic sub-steps for each internal step nn_baro <= 2.5 nn_baro
   REAL(wp),SAVE :: rdtbt   ! Barotropic time step

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) :: &
                    wgtbtp1, &              ! Primary weights used for time filtering of barotropic variables
                    wgtbtp2                 ! Secondary weights used for time filtering of barotropic variables

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  zwz          ! ff/h at F points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ftnw, ftne   ! triad of coriolis parameter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ftsw, ftse   ! (only used with een vorticity scheme)

   ! Arrays below are saved to allow testing of the "no time averaging" option
   ! If this option is not retained, these could be replaced by temporary arrays
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sshbb_e, sshb_e, & ! Instantaneous barotropic arrays
                                                   ubb_e, ub_e,     &
                                                   vbb_e, vb_e

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.5 , NEMO Consortium (2013)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_spg_ts_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!----------------------------------------------------------------------
      ierr(:) = 0

      ALLOCATE( sshb_e(jpi,jpj), sshbb_e(jpi,jpj), &
         &      ub_e(jpi,jpj)  , vb_e(jpi,jpj)   , &
         &      ubb_e(jpi,jpj) , vbb_e(jpi,jpj)  , STAT= ierr(1) )

      ALLOCATE( wgtbtp1(3*nn_baro), wgtbtp2(3*nn_baro), zwz(jpi,jpj), STAT= ierr(2) )

      IF( ln_dynvor_een .or. ln_dynvor_een_old ) ALLOCATE( ftnw(jpi,jpj) , ftne(jpi,jpj) , & 
                                                    &      ftsw(jpi,jpj) , ftse(jpi,jpj) , STAT=ierr(3) )

      dyn_spg_ts_alloc = MAXVAL(ierr(:))

      IF( lk_mpp                )   CALL mpp_sum( dyn_spg_ts_alloc )
      IF( dyn_spg_ts_alloc /= 0 )   CALL ctl_warn('dynspg_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION dyn_spg_ts_alloc

   SUBROUTINE dyn_spg_ts( kt )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :   
      !!      -Compute the now trend due to the explicit time stepping
      !!      of the quasi-linear barotropic system. 
      !!
      !! ** Method  :  
      !!      Barotropic variables are advanced from internal time steps
      !!      "n"   to "n+1" if ln_bt_fw=T
      !!      or from 
      !!      "n-1" to "n+1" if ln_bt_fw=F
      !!      thanks to a generalized forward-backward time stepping (see ref. below).
      !!
      !! ** Action :
      !!      -Update the filtered free surface at step "n+1"      : ssha
      !!      -Update filtered barotropic velocities at step "n+1" : ua_b, va_b
      !!      -Compute barotropic advective velocities at step "n" : un_adv, vn_adv
      !!      These are used to advect tracers and are compliant with discrete
      !!      continuity equation taken at the baroclinic time steps. This 
      !!      ensures tracers conservation.
      !!      -Update 3d trend (ua, va) with barotropic component.
      !!
      !! References : Shchepetkin, A.F. and J.C. McWilliams, 2005: 
      !!              The regional oceanic modeling system (ROMS): 
      !!              a split-explicit, free-surface,
      !!              topography-following-coordinate oceanic model. 
      !!              Ocean Modelling, 9, 347-404. 
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in)  ::   kt   ! ocean time-step index
      !
      LOGICAL  ::   ll_fw_start			! if true, forward integration 
      LOGICAL  ::   ll_init			    ! if true, special startup of 2d equations
      INTEGER  ::   ji, jj, jk, jn   		! dummy loop indices
      INTEGER  ::   ikbu, ikbv, noffset    	! local integers
      REAL(wp) ::   zraur, z1_2dt_b, z2dt_bf  	! local scalars
      REAL(wp) ::   zx1, zy1, zx2, zy2         !   -      -
      REAL(wp) ::   z1_12, z1_8, z1_4, z1_2 	  !   -      -
      REAL(wp) ::   zu_spg, zv_spg             !   -      -
      REAL(wp) ::   zhura, zhvra         	     !   -      -
      REAL(wp) ::   za0, za1, za2, za3		       !   -      -
      !
      REAL(wp), POINTER, DIMENSION(:,:) :: zun_e, zvn_e, zsshp2_e
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_trd, zv_trd, zu_frc, zv_frc, zssh_frc
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_sum, zv_sum, zwx, zwy, zhdiv
      REAL(wp), POINTER, DIMENSION(:,:) :: zhup2_e, zhvp2_e, zhust_e, zhvst_e
      REAL(wp), POINTER, DIMENSION(:,:) :: zsshu_a, zsshv_a
      REAL(wp), POINTER, DIMENSION(:,:) :: zhf
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_ts')
      !
      !                                         !* Allocate temporary arrays
      CALL wrk_alloc( jpi, jpj, zsshp2_e, zhdiv )
      CALL wrk_alloc( jpi, jpj, zu_trd, zv_trd, zun_e, zvn_e  )
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zu_sum, zv_sum, zssh_frc, zu_frc, zv_frc)
      CALL wrk_alloc( jpi, jpj, zhup2_e, zhvp2_e, zhust_e, zhvst_e)
      CALL wrk_alloc( jpi, jpj, zsshu_a, zsshv_a                                   )
      CALL wrk_alloc( jpi, jpj, zhf )
      !
      !                                         !* Local constant initialization
      z1_12 = 1._wp / 12._wp 
      z1_8  = 0.125_wp                                   
      z1_4  = 0.25_wp
      z1_2  = 0.5_wp     
      zraur = 1._wp / rau0
      !
      IF( kt == nit000 .AND. neuler == 0 ) THEN  	! reciprocal of baroclinic time step 
        z2dt_bf = rdt
      ELSE
        z2dt_bf = 2.0_wp * rdt
      ENDIF
      z1_2dt_b = 1.0_wp / z2dt_bf 
      !
      ll_init = ln_bt_av                       		! if no time averaging, then no specific restart 
      ll_fw_start = .FALSE.
      !
                                     	                ! time offset in steps for bdy data update
      IF (.NOT.ln_bt_fw) THEN ; noffset=-nn_baro ; ELSE ;  noffset = 0 ; ENDIF
      !
      IF( kt == nit000 ) THEN             	!* initialisation
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
         IF(lwp) WRITE(numout,*)
         !
         IF (neuler==0) ll_init=.TRUE.
         !
         IF (ln_bt_fw.OR.(neuler==0)) THEN
           ll_fw_start=.TRUE.
           noffset = 0
         ELSE
           ll_fw_start=.FALSE.
         ENDIF
         !
         ! Set averaging weights and cycle length:
         CALL ts_wgt(ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2)
         !
         !
      ENDIF
      !
      ! Set arrays to remove/compute coriolis trend.
      ! Do it once at kt=nit000 if volume is fixed, else at each long time step.
      ! Note that these arrays are also used during barotropic loop. These are however frozen
      ! although they should be updated in the variable volume case. Not a big approximation.
      ! To remove this approximation, copy lines below inside barotropic loop
      ! and update depths at T-F points (ht and zhf resp.) at each barotropic time step
      !
      IF ( kt == nit000 .OR. lk_vvl ) THEN
         IF ( ln_dynvor_een_old ) THEN
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zwz(ji,jj) =   ( ht(ji  ,jj+1) + ht(ji+1,jj+1) +                    &
                        &          ht(ji  ,jj  ) + ht(ji+1,jj  )   ) / 4._wp  
                  IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = 1._wp / zwz(ji,jj)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1._wp )
            zwz(:,:) = ff(:,:) * zwz(:,:)

            ftne(1,:) = 0._wp ; ftnw(1,:) = 0._wp ; ftse(1,:) = 0._wp ; ftsw(1,:) = 0._wp
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  ftne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
                  ftnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
                  ftse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
                  ftsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
               END DO
            END DO
         ELSE IF ( ln_dynvor_een ) THEN
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zwz(ji,jj) =   ( ht(ji  ,jj+1) + ht(ji+1,jj+1) +                     &
                        &          ht(ji  ,jj  ) + ht(ji+1,jj  )   )                   &
                        &      / ( MAX( 1.0_wp, tmask(ji  ,jj+1, 1) + tmask(ji+1,jj+1, 1) +    &
                        &                       tmask(ji  ,jj  , 1) + tmask(ji+1,jj  , 1) ) )
                  IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = 1._wp / zwz(ji,jj)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1._wp )
            zwz(:,:) = ff(:,:) * zwz(:,:)

            ftne(1,:) = 0._wp ; ftnw(1,:) = 0._wp ; ftse(1,:) = 0._wp ; ftsw(1,:) = 0._wp
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  ftne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
                  ftnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
                  ftse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
                  ftsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
               END DO
            END DO
         ELSE
            zwz(:,:) = 0._wp
            zhf(:,:) = 0.
            IF ( .not. ln_sco ) THEN
! JC: It not clear yet what should be the depth at f-points over land in z-coordinate case
!     Set it to zero for the time being 
!              IF( rn_hmin < 0._wp ) THEN    ;   jk = - INT( rn_hmin )                                      ! from a nb of level
!              ELSE                          ;   jk = MINLOC( gdepw_0, mask = gdepw_0 > rn_hmin, dim = 1 )  ! from a depth
!              ENDIF
!              zhf(:,:) = gdepw_0(:,:,jk+1)
            ELSE
               zhf(:,:) = hbatf(:,:)
            END IF

            DO jj = 1, jpjm1
               zhf(:,jj) = zhf(:,jj)*(1._wp- umask(:,jj,1) * umask(:,jj+1,1))
            END DO

            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  zhf(:,jj) = zhf(:,jj) + fse3f_n(:,jj,jk) * umask(:,jj,jk) * umask(:,jj+1,jk)
               END DO
            END DO
            CALL lbc_lnk( zhf, 'F', 1._wp )
            ! JC: TBC. hf should be greater than 0 
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( zhf(ji,jj) /= 0._wp )   zwz(ji,jj) = 1._wp / zhf(ji,jj) ! zhf is actually hf here but it saves an array
               END DO
            END DO
            zwz(:,:) = ff(:,:) * zwz(:,:)
         ENDIF
      ENDIF
      !
      ! If forward start at previous time step, and centered integration, 
      ! then update averaging weights:
      IF ((.NOT.ln_bt_fw).AND.((neuler==0).AND.(kt==nit000+1))) THEN
         ll_fw_start=.FALSE.
         CALL ts_wgt(ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2)
      ENDIF
                          
      ! -----------------------------------------------------------------------------
      !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
      ! -----------------------------------------------------------------------------
      !      
      !
      !                                   !* e3*d/dt(Ua) (Vertically integrated)
      !                                   ! --------------------------------------------------
      zu_frc(:,:) = 0._wp
      zv_frc(:,:) = 0._wp
      !
      DO jk = 1, jpkm1
         zu_frc(:,:) = zu_frc(:,:) + fse3u_n(:,:,jk) * ua(:,:,jk) * umask(:,:,jk)
         zv_frc(:,:) = zv_frc(:,:) + fse3v_n(:,:,jk) * va(:,:,jk) * vmask(:,:,jk)         
      END DO
      !
      zu_frc(:,:) = zu_frc(:,:) * hur(:,:)
      zv_frc(:,:) = zv_frc(:,:) * hvr(:,:)
      !
      !
      !                                   !* baroclinic momentum trend (remove the vertical mean trend)
      DO jk = 1, jpkm1                    ! -----------------------------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - zu_frc(ji,jj) * umask(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - zv_frc(ji,jj) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !                                   !* barotropic Coriolis trends (vorticity scheme dependent)
      !                                   ! --------------------------------------------------------
      zwx(:,:) = un_b(:,:) * hu(:,:) * e2u(:,:)        ! now fluxes 
      zwy(:,:) = vn_b(:,:) * hv(:,:) * e1v(:,:)
      !
      IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
               zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
               zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               ! energy conserving formulation for planetary vorticity term
               zu_trd(ji,jj) = z1_4 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               zv_trd(ji,jj) =-z1_4 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
            END DO
         END DO
         !
      ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 =   z1_8 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) &
                 &            + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 = - z1_8 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) &
                 &            + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               zu_trd(ji,jj)  = zy1 * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
               zv_trd(ji,jj)  = zx1 * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
            END DO
         END DO
         !
      ELSEIF ( ln_dynvor_een .or. ln_dynvor_een_old ) THEN  ! enstrophy and energy conserving scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zu_trd(ji,jj) = + z1_12 / e1u(ji,jj) * (  ftne(ji,jj  ) * zwy(ji  ,jj  ) &
                &                                      + ftnw(ji+1,jj) * zwy(ji+1,jj  ) &
                &                                      + ftse(ji,jj  ) * zwy(ji  ,jj-1) &
                &                                      + ftsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zv_trd(ji,jj) = - z1_12 / e2v(ji,jj) * (  ftsw(ji,jj+1) * zwx(ji-1,jj+1) &
                &                                      + ftse(ji,jj+1) * zwx(ji  ,jj+1) &
                &                                      + ftnw(ji,jj  ) * zwx(ji-1,jj  ) &
                &                                      + ftne(ji,jj  ) * zwx(ji  ,jj  ) )
            END DO
         END DO
         !
      ENDIF 
      !
      !                                   !* Right-Hand-Side of the barotropic momentum equation
      !                                   ! ----------------------------------------------------
      IF( lk_vvl ) THEN                         ! Variable volume : remove surface pressure gradient
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zu_trd(ji,jj) = zu_trd(ji,jj) - grav * (  sshn(ji+1,jj  ) - sshn(ji  ,jj  )  ) / e1u(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) - grav * (  sshn(ji  ,jj+1) - sshn(ji  ,jj  )  ) / e2v(ji,jj)
            END DO
         END DO
      ENDIF

      DO jj = 2, jpjm1                          ! Remove coriolis term (and possibly spg) from barotropic trend
         DO ji = fs_2, fs_jpim1
             zu_frc(ji,jj) = zu_frc(ji,jj) - zu_trd(ji,jj) * umask(ji,jj,1)
             zv_frc(ji,jj) = zv_frc(ji,jj) - zv_trd(ji,jj) * vmask(ji,jj,1)
          END DO
      END DO 
      !
      !						! Add bottom stress contribution from baroclinic velocities:      
      IF (ln_bt_fw) THEN
         DO jj = 2, jpjm1                          
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ikbu = mbku(ji,jj)       
               ikbv = mbkv(ji,jj)    
               zwx(ji,jj) = un(ji,jj,ikbu) - un_b(ji,jj) ! NOW bottom baroclinic velocities
               zwy(ji,jj) = vn(ji,jj,ikbv) - vn_b(ji,jj)
            END DO
         END DO
      ELSE
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ikbu = mbku(ji,jj)       
               ikbv = mbkv(ji,jj)    
               zwx(ji,jj) = ub(ji,jj,ikbu) - ub_b(ji,jj) ! BEFORE bottom baroclinic velocities
               zwy(ji,jj) = vb(ji,jj,ikbv) - vb_b(ji,jj)
            END DO
         END DO
      ENDIF
      !
      ! Note that the "unclipped" bottom friction parameter is used even with explicit drag
      zu_frc(:,:) = zu_frc(:,:) + hur(:,:) * bfrua(:,:) * zwx(:,:)
      zv_frc(:,:) = zv_frc(:,:) + hvr(:,:) * bfrva(:,:) * zwy(:,:)
      !       
      IF (ln_bt_fw) THEN                        ! Add wind forcing
         zu_frc(:,:) =  zu_frc(:,:) + zraur * utau(:,:) * hur(:,:)
         zv_frc(:,:) =  zv_frc(:,:) + zraur * vtau(:,:) * hvr(:,:)
      ELSE
         zu_frc(:,:) =  zu_frc(:,:) + zraur * z1_2 * ( utau_b(:,:) + utau(:,:) ) * hur(:,:)
         zv_frc(:,:) =  zv_frc(:,:) + zraur * z1_2 * ( vtau_b(:,:) + vtau(:,:) ) * hvr(:,:)
      ENDIF  
      !
      IF ( ln_apr_dyn ) THEN                    ! Add atm pressure forcing
         IF (ln_bt_fw) THEN
            DO jj = 2, jpjm1              
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_spg =  grav * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj) ) /e1u(ji,jj)
                  zv_spg =  grav * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj) ) /e2v(ji,jj)
                  zu_frc(ji,jj) = zu_frc(ji,jj) + zu_spg
                  zv_frc(ji,jj) = zv_frc(ji,jj) + zv_spg
               END DO
            END DO
         ELSE
            DO jj = 2, jpjm1              
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_spg =  grav * z1_2 * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj)    &
                      &                    + ssh_ibb(ji+1,jj  ) - ssh_ibb(ji,jj)  ) /e1u(ji,jj)
                  zv_spg =  grav * z1_2 * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj)    &
                      &                    + ssh_ibb(ji  ,jj+1) - ssh_ibb(ji,jj)  ) /e2v(ji,jj)
                  zu_frc(ji,jj) = zu_frc(ji,jj) + zu_spg
                  zv_frc(ji,jj) = zv_frc(ji,jj) + zv_spg
               END DO
            END DO
         ENDIF 
      ENDIF
      !                                   !* Right-Hand-Side of the barotropic ssh equation
      !                                   ! -----------------------------------------------
      !                                         ! Surface net water flux and rivers
      IF (ln_bt_fw) THEN
         zssh_frc(:,:) = zraur * ( emp(:,:) - rnf(:,:) + fwfisf(:,:) )
      ELSE
         zssh_frc(:,:) = zraur * z1_2 * (  emp(:,:) + emp_b(:,:) - rnf(:,:) - rnf_b(:,:)   &
                &                        + fwfisf(:,:) + fwfisf_b(:,:)                     )
      ENDIF
#if defined key_asminc
      !                                         ! Include the IAU weighted SSH increment
      IF( lk_asminc .AND. ln_sshinc .AND. ln_asmiau ) THEN
         zssh_frc(:,:) = zssh_frc(:,:) - ssh_iau(:,:)
      ENDIF
#endif
      !                                   !* Fill boundary data arrays for AGRIF
      !                                   ! ------------------------------------
#if defined key_agrif
         IF( .NOT.Agrif_Root() ) CALL agrif_dta_ts( kt )
#endif
      !
      ! -----------------------------------------------------------------------
      !  Phase 2 : Integration of the barotropic equations 
      ! -----------------------------------------------------------------------
      !
      !                                             ! ==================== !
      !                                             !    Initialisations   !
      !                                             ! ==================== !  
      ! Initialize barotropic variables:      
      IF( ll_init )THEN
         sshbb_e(:,:) = 0._wp
         ubb_e  (:,:) = 0._wp
         vbb_e  (:,:) = 0._wp
         sshb_e (:,:) = 0._wp
         ub_e   (:,:) = 0._wp
         vb_e   (:,:) = 0._wp
      ENDIF
      !
      IF (ln_bt_fw) THEN                  ! FORWARD integration: start from NOW fields                    
         sshn_e(:,:) = sshn (:,:)            
         zun_e (:,:) = un_b (:,:)            
         zvn_e (:,:) = vn_b (:,:)
         !
         hu_e  (:,:) = hu   (:,:)       
         hv_e  (:,:) = hv   (:,:) 
         hur_e (:,:) = hur  (:,:)    
         hvr_e (:,:) = hvr  (:,:)
      ELSE                                ! CENTRED integration: start from BEFORE fields
         sshn_e(:,:) = sshb (:,:)
         zun_e (:,:) = ub_b (:,:)         
         zvn_e (:,:) = vb_b (:,:)
         !
         hu_e  (:,:) = hu_b (:,:)       
         hv_e  (:,:) = hv_b (:,:) 
         hur_e (:,:) = hur_b(:,:)    
         hvr_e (:,:) = hvr_b(:,:)
      ENDIF
      !
      !
      !
      ! Initialize sums:
      ua_b  (:,:) = 0._wp       ! After barotropic velocities (or transport if flux form)          
      va_b  (:,:) = 0._wp
      ssha  (:,:) = 0._wp       ! Sum for after averaged sea level
      zu_sum(:,:) = 0._wp       ! Sum for now transport issued from ts loop
      zv_sum(:,:) = 0._wp
      !                                             ! ==================== !
      DO jn = 1, icycle                             !  sub-time-step loop  !
         !                                          ! ==================== !
         !                                                !* Update the forcing (BDY and tides)
         !                                                !  ------------------
         ! Update only tidal forcing at open boundaries
#if defined key_tide
         IF ( lk_bdy .AND. lk_tide ) CALL bdy_dta_tides( kt, kit=jn, time_offset=(noffset+1) )
         IF ( ln_tide_pot .AND. lk_tide ) CALL upd_tide( kt, kit=jn, time_offset=noffset )
#endif
         !
         ! Set extrapolation coefficients for predictor step:
         IF ((jn<3).AND.ll_init) THEN      ! Forward           
           za1 = 1._wp                                          
           za2 = 0._wp                        
           za3 = 0._wp                        
         ELSE                              ! AB3-AM4 Coefficients: bet=0.281105 
           za1 =  1.781105_wp              ! za1 =   3/2 +   bet
           za2 = -1.06221_wp               ! za2 = -(1/2 + 2*bet)
           za3 =  0.281105_wp              ! za3 = bet
         ENDIF

         ! Extrapolate barotropic velocities at step jit+0.5:
         ua_e(:,:) = za1 * zun_e(:,:) + za2 * ub_e(:,:) + za3 * ubb_e(:,:)
         va_e(:,:) = za1 * zvn_e(:,:) + za2 * vb_e(:,:) + za3 * vbb_e(:,:)

         IF( lk_vvl ) THEN                                !* Update ocean depth (variable volume case only)
            !                                             !  ------------------
            ! Extrapolate Sea Level at step jit+0.5:
            zsshp2_e(:,:) = za1 * sshn_e(:,:)  + za2 * sshb_e(:,:) + za3 * sshbb_e(:,:)
            !
            DO jj = 2, jpjm1                                    ! Sea Surface Height at u- & v-points
               DO ji = 2, fs_jpim1   ! Vector opt.
                  zwx(ji,jj) = z1_2 * umask(ji,jj,1)  * r1_e12u(ji,jj)     &
                     &              * ( e12t(ji  ,jj) * zsshp2_e(ji  ,jj)  &
                     &              +   e12t(ji+1,jj) * zsshp2_e(ji+1,jj) )
                  zwy(ji,jj) = z1_2 * vmask(ji,jj,1)  * r1_e12v(ji,jj)     &
                     &              * ( e12t(ji,jj  ) * zsshp2_e(ji,jj  )  &
                     &              +   e12t(ji,jj+1) * zsshp2_e(ji,jj+1) )
               END DO
            END DO
            CALL lbc_lnk_multi( zwx, 'U', 1._wp, zwy, 'V', 1._wp )
            !
            zhup2_e (:,:) = hu_0(:,:) + zwx(:,:)                ! Ocean depth at U- and V-points
            zhvp2_e (:,:) = hv_0(:,:) + zwy(:,:)
         ELSE
            zhup2_e (:,:) = hu(:,:)
            zhvp2_e (:,:) = hv(:,:)
         ENDIF
         !                                                !* after ssh
         !                                                !  -----------
         ! One should enforce volume conservation at open boundaries here
         ! considering fluxes below:
         !
         zwx(:,:) = e2u(:,:) * ua_e(:,:) * zhup2_e(:,:)         ! fluxes at jn+0.5
         zwy(:,:) = e1v(:,:) * va_e(:,:) * zhvp2_e(:,:)
         !
#if defined key_agrif
         ! Set fluxes during predictor step to ensure 
         ! volume conservation
         IF( (.NOT.Agrif_Root()).AND.ln_bt_fw ) THEN
            IF((nbondi == -1).OR.(nbondi == 2)) THEN
               DO jj=1,jpj
                  zwx(2,jj) = ubdy_w(jj) * e2u(2,jj)
               END DO
            ENDIF
            IF((nbondi ==  1).OR.(nbondi == 2)) THEN
               DO jj=1,jpj
                  zwx(nlci-2,jj) = ubdy_e(jj) * e2u(nlci-2,jj)
               END DO
            ENDIF
            IF((nbondj == -1).OR.(nbondj == 2)) THEN
               DO ji=1,jpi
                  zwy(ji,2) = vbdy_s(ji) * e1v(ji,2)
               END DO
            ENDIF
            IF((nbondj ==  1).OR.(nbondj == 2)) THEN
               DO ji=1,jpi
                  zwy(ji,nlcj-2) = vbdy_n(ji) * e1v(ji,nlcj-2)
               END DO
            ENDIF
         ENDIF
#endif
         !
         ! Sum over sub-time-steps to compute advective velocities
         za2 = wgtbtp2(jn)
         zu_sum  (:,:) = zu_sum  (:,:) + za2 * zwx  (:,:) / e2u  (:,:)
         zv_sum  (:,:) = zv_sum  (:,:) + za2 * zwy  (:,:) / e1v  (:,:)
         !
         ! Set next sea level:
         DO jj = 2, jpjm1                                 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zhdiv(ji,jj) = (   zwx(ji,jj) - zwx(ji-1,jj)   &
                  &             + zwy(ji,jj) - zwy(ji,jj-1)   ) * r1_e12t(ji,jj)
            END DO
         END DO
         ssha_e(:,:) = (  sshn_e(:,:) - rdtbt * ( zssh_frc(:,:) + zhdiv(:,:) )  ) * tmask(:,:,1)
         CALL lbc_lnk( ssha_e, 'T',  1._wp )

#if defined key_bdy
         ! Duplicate sea level across open boundaries (this is only cosmetic if lk_vvl=.false.)
         IF (lk_bdy) CALL bdy_ssh( ssha_e )
#endif
#if defined key_agrif
         IF( .NOT.Agrif_Root() ) CALL agrif_ssh_ts( jn )
#endif
         !  
         ! Sea Surface Height at u-,v-points (vvl case only)
         IF ( lk_vvl ) THEN                                
            DO jj = 2, jpjm1
               DO ji = 2, jpim1      ! NO Vector Opt.
                  zsshu_a(ji,jj) = z1_2 * umask(ji,jj,1)  * r1_e12u(ji,jj)  &
                     &              * ( e12t(ji  ,jj  ) * ssha_e(ji  ,jj  ) &
                     &              +   e12t(ji+1,jj  ) * ssha_e(ji+1,jj  ) )
                  zsshv_a(ji,jj) = z1_2 * vmask(ji,jj,1)  * r1_e12v(ji,jj)  &
                     &              * ( e12t(ji  ,jj  ) * ssha_e(ji  ,jj  ) &
                     &              +   e12t(ji  ,jj+1) * ssha_e(ji  ,jj+1) )
               END DO
            END DO
            CALL lbc_lnk_multi( zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp )
         ENDIF   
         !                                 
         ! Half-step back interpolation of SSH for surface pressure computation:
         !----------------------------------------------------------------------
         IF ((jn==1).AND.ll_init) THEN
           za0=1._wp                        ! Forward-backward
           za1=0._wp                           
           za2=0._wp
           za3=0._wp
         ELSEIF ((jn==2).AND.ll_init) THEN  ! AB2-AM3 Coefficients; bet=0 ; gam=-1/6 ; eps=1/12
           za0= 1.0833333333333_wp          ! za0 = 1-gam-eps
           za1=-0.1666666666666_wp          ! za1 = gam
           za2= 0.0833333333333_wp          ! za2 = eps
           za3= 0._wp              
         ELSE                               ! AB3-AM4 Coefficients; bet=0.281105 ; eps=0.013 ; gam=0.0880 
           za0=0.614_wp                     ! za0 = 1/2 +   gam + 2*eps    
           za1=0.285_wp                     ! za1 = 1/2 - 2*gam - 3*eps 
           za2=0.088_wp                     ! za2 = gam
           za3=0.013_wp                     ! za3 = eps
         ENDIF

         zsshp2_e(:,:) = za0 *  ssha_e(:,:) + za1 *  sshn_e (:,:) &
          &            + za2 *  sshb_e(:,:) + za3 *  sshbb_e(:,:)

         !
         ! Compute associated depths at U and V points:
         IF ( lk_vvl.AND.(.NOT.ln_dynadv_vec) ) THEN      
            !                                        
            DO jj = 2, jpjm1                            
               DO ji = 2, jpim1
                  zx1 = z1_2 * umask(ji  ,jj,1) *  r1_e12u(ji  ,jj)    &
                     &      * ( e12t(ji  ,jj  ) * zsshp2_e(ji  ,jj)    &
                     &      +   e12t(ji+1,jj  ) * zsshp2_e(ji+1,jj  ) )
                  zy1 = z1_2 * vmask(ji  ,jj,1) *  r1_e12v(ji  ,jj  )  &
                     &       * ( e12t(ji ,jj  ) * zsshp2_e(ji  ,jj  )  &
                     &       +   e12t(ji ,jj+1) * zsshp2_e(ji  ,jj+1) )
                  zhust_e(ji,jj) = hu_0(ji,jj) + zx1 
                  zhvst_e(ji,jj) = hv_0(ji,jj) + zy1
               END DO
            END DO
         ENDIF
         !
         ! Add Coriolis trend:
         ! zwz array below or triads normally depend on sea level with key_vvl and should be updated
         ! at each time step. We however keep them constant here for optimization.
         ! Recall that zwx and zwy arrays hold fluxes at this stage:
         ! zwx(:,:) = e2u(:,:) * ua_e(:,:) * zhup2_e(:,:)   ! fluxes at jn+0.5
         ! zwy(:,:) = e1v(:,:) * va_e(:,:) * zhvp2_e(:,:)
         !
         IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      !==  energy conserving or mixed scheme  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zy1 = ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
                  zy2 = ( zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
                  zx1 = ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
                  zx2 = ( zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
                  zu_trd(ji,jj) = z1_4 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
                  zv_trd(ji,jj) =-z1_4 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
               END DO
            END DO
            !
         ELSEIF ( ln_dynvor_ens ) THEN                    !==  enstrophy conserving scheme  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zy1 =   z1_8 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) &
                   &             + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
                  zx1 = - z1_8 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) &
                   &             + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
                  zu_trd(ji,jj)  = zy1 * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
                  zv_trd(ji,jj)  = zx1 * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
               END DO
            END DO
            !
         ELSEIF ( ln_dynvor_een .or. ln_dynvor_een_old ) THEN !==  energy and enstrophy conserving scheme  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_trd(ji,jj) = + z1_12 / e1u(ji,jj) * (  ftne(ji,jj  ) * zwy(ji  ,jj  ) &
                     &                                    + ftnw(ji+1,jj) * zwy(ji+1,jj  ) &
                     &                                    + ftse(ji,jj  ) * zwy(ji  ,jj-1) & 
                     &                                    + ftsw(ji+1,jj) * zwy(ji+1,jj-1) )
                  zv_trd(ji,jj) = - z1_12 / e2v(ji,jj) * (  ftsw(ji,jj+1) * zwx(ji-1,jj+1) & 
                     &                                    + ftse(ji,jj+1) * zwx(ji  ,jj+1) &
                     &                                    + ftnw(ji,jj  ) * zwx(ji-1,jj  ) & 
                     &                                    + ftne(ji,jj  ) * zwx(ji  ,jj  ) )
               END DO
            END DO
            ! 
         ENDIF
         !
         ! Add tidal astronomical forcing if defined
         IF ( lk_tide.AND.ln_tide_pot ) THEN
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_spg = grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) / e1u(ji,jj)
                  zv_spg = grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) / e2v(ji,jj)
                  zu_trd(ji,jj) = zu_trd(ji,jj) + zu_spg
                  zv_trd(ji,jj) = zv_trd(ji,jj) + zv_spg
               END DO
            END DO
         ENDIF
         !
         ! Add bottom stresses:
         zu_trd(:,:) = zu_trd(:,:) + bfrua(:,:) * zun_e(:,:) * hur_e(:,:)
         zv_trd(:,:) = zv_trd(:,:) + bfrva(:,:) * zvn_e(:,:) * hvr_e(:,:)
         !
         ! Surface pressure trend:
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! Add surface pressure gradient
               zu_spg = - grav * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) / e1u(ji,jj)
               zv_spg = - grav * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) / e2v(ji,jj)
               zwx(ji,jj) = zu_spg
               zwy(ji,jj) = zv_spg
            END DO
         END DO
         !
         ! Set next velocities:
         IF( ln_dynadv_vec .OR. (.NOT. lk_vvl) ) THEN 	! Vector form
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua_e(ji,jj) = (                                zun_e(ji,jj)   & 
                            &     + rdtbt * (                      zwx(ji,jj)   &
                            &                                 + zu_trd(ji,jj)   &
                            &                                 + zu_frc(ji,jj) ) & 
                            &   ) * umask(ji,jj,1)

                  va_e(ji,jj) = (                                zvn_e(ji,jj)   &
                            &     + rdtbt * (                      zwy(ji,jj)   &
                            &                                 + zv_trd(ji,jj)   &
                            &                                 + zv_frc(ji,jj) ) &
                            &   ) * vmask(ji,jj,1)
               END DO
            END DO

         ELSE					 	! Flux form
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.

                  zhura = umask(ji,jj,1)/(hu_0(ji,jj) + zsshu_a(ji,jj) + 1._wp - umask(ji,jj,1))
                  zhvra = vmask(ji,jj,1)/(hv_0(ji,jj) + zsshv_a(ji,jj) + 1._wp - vmask(ji,jj,1))

                  ua_e(ji,jj) = (                hu_e(ji,jj)  *  zun_e(ji,jj)   & 
                            &     + rdtbt * ( zhust_e(ji,jj)  *    zwx(ji,jj)   & 
                            &               + zhup2_e(ji,jj)  * zu_trd(ji,jj)   &
                            &               +      hu(ji,jj)  * zu_frc(ji,jj) ) &
                            &   ) * zhura

                  va_e(ji,jj) = (                hv_e(ji,jj)  *  zvn_e(ji,jj)   &
                            &     + rdtbt * ( zhvst_e(ji,jj)  *    zwy(ji,jj)   &
                            &               + zhvp2_e(ji,jj)  * zv_trd(ji,jj)   &
                            &               +      hv(ji,jj)  * zv_frc(ji,jj) ) &
                            &   ) * zhvra
               END DO
            END DO
         ENDIF
         !
         IF( lk_vvl ) THEN                             !* Update ocean depth (variable volume case only)
            !                                          !  ----------------------------------------------        
            hu_e (:,:) = hu_0(:,:) + zsshu_a(:,:)
            hv_e (:,:) = hv_0(:,:) + zsshv_a(:,:)
            hur_e(:,:) = umask(:,:,1) / ( hu_e(:,:) + 1._wp - umask(:,:,1) )
            hvr_e(:,:) = vmask(:,:,1) / ( hv_e(:,:) + 1._wp - vmask(:,:,1) )
            !
         ENDIF
         !                                                 !* domain lateral boundary
         !                                                 !  -----------------------
         !
         CALL lbc_lnk_multi( ua_e, 'U', -1._wp, va_e , 'V', -1._wp )

#if defined key_bdy  
                                                           ! open boundaries
         IF( lk_bdy ) CALL bdy_dyn2d( jn, ua_e, va_e, zun_e, zvn_e, hur_e, hvr_e, ssha_e )
#endif
#if defined key_agrif                                                           
         IF( .NOT.Agrif_Root() )  CALL agrif_dyn_ts( jn )  ! Agrif
#endif
         !                                             !* Swap
         !                                             !  ----
         ubb_e  (:,:) = ub_e  (:,:)
         ub_e   (:,:) = zun_e (:,:)
         zun_e  (:,:) = ua_e  (:,:)
         !
         vbb_e  (:,:) = vb_e  (:,:)
         vb_e   (:,:) = zvn_e (:,:)
         zvn_e  (:,:) = va_e  (:,:)
         !
         sshbb_e(:,:) = sshb_e(:,:)
         sshb_e (:,:) = sshn_e(:,:)
         sshn_e (:,:) = ssha_e(:,:)

         !                                             !* Sum over whole bt loop
         !                                             !  ----------------------
         za1 = wgtbtp1(jn)                                    
         IF (( ln_dynadv_vec ).OR. (.NOT. lk_vvl)) THEN    ! Sum velocities
            ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) 
            va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) 
         ELSE                     			   ! Sum transports
            ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) * hu_e (:,:)
            va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) * hv_e (:,:)
         ENDIF
         !                       			   ! Sum sea level
         ssha(:,:) = ssha(:,:) + za1 * ssha_e(:,:)
         !                                                 ! ==================== !
      END DO                                               !        end loop      !
      !                                                    ! ==================== !
      ! -----------------------------------------------------------------------------
      ! Phase 3. update the general trend with the barotropic trend
      ! -----------------------------------------------------------------------------
      !
      ! At this stage ssha holds a time averaged value
      !                                                ! Sea Surface Height at u-,v- and f-points
      IF( lk_vvl ) THEN                                ! (required only in key_vvl case)
         DO jj = 1, jpjm1
            DO ji = 1, jpim1      ! NO Vector Opt.
               zsshu_a(ji,jj) = z1_2 * umask(ji,jj,1)  * r1_e12u(ji,jj) &
                  &              * ( e12t(ji  ,jj) * ssha(ji  ,jj)    &
                  &              +   e12t(ji+1,jj) * ssha(ji+1,jj) )
               zsshv_a(ji,jj) = z1_2 * vmask(ji,jj,1)  * r1_e12v(ji,jj) &
                  &              * ( e12t(ji,jj  ) * ssha(ji,jj  )    &
                  &              +   e12t(ji,jj+1) * ssha(ji,jj+1) )
            END DO
         END DO
         CALL lbc_lnk_multi( zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp ) ! Boundary conditions
      ENDIF
      !
      ! Set advection velocity correction:
      IF (((kt==nit000).AND.(neuler==0)).OR.(.NOT.ln_bt_fw)) THEN     
         un_adv(:,:) = zu_sum(:,:)*hur(:,:)
         vn_adv(:,:) = zv_sum(:,:)*hvr(:,:)
      ELSE
         un_adv(:,:) = z1_2 * ( ub2_b(:,:) + zu_sum(:,:)) * hur(:,:)
         vn_adv(:,:) = z1_2 * ( vb2_b(:,:) + zv_sum(:,:)) * hvr(:,:)
      END IF

      IF (ln_bt_fw) THEN ! Save integrated transport for next computation
         ub2_b(:,:) = zu_sum(:,:)
         vb2_b(:,:) = zv_sum(:,:)
      ENDIF
      !
      ! Update barotropic trend:
      IF (( ln_dynadv_vec ).OR. (.NOT. lk_vvl)) THEN
         DO jk=1,jpkm1
            ua(:,:,jk) = ua(:,:,jk) + ( ua_b(:,:) - ub_b(:,:) ) * z1_2dt_b
            va(:,:,jk) = va(:,:,jk) + ( va_b(:,:) - vb_b(:,:) ) * z1_2dt_b
         END DO
      ELSE
         DO jk=1,jpkm1
            ua(:,:,jk) = ua(:,:,jk) + hur(:,:) * ( ua_b(:,:) - ub_b(:,:) * hu_b(:,:) ) * z1_2dt_b
            va(:,:,jk) = va(:,:,jk) + hvr(:,:) * ( va_b(:,:) - vb_b(:,:) * hv_b(:,:) ) * z1_2dt_b
         END DO
         ! Save barotropic velocities not transport:
         ua_b  (:,:) =  ua_b(:,:) / ( hu_0(:,:) + zsshu_a(:,:) + 1._wp - umask(:,:,1) )
         va_b  (:,:) =  va_b(:,:) / ( hv_0(:,:) + zsshv_a(:,:) + 1._wp - vmask(:,:,1) )
      ENDIF
      !
      DO jk = 1, jpkm1
         ! Correct velocities:
         un(:,:,jk) = ( un(:,:,jk) + un_adv(:,:) - un_b(:,:) )*umask(:,:,jk)
         vn(:,:,jk) = ( vn(:,:,jk) + vn_adv(:,:) - vn_b(:,:) )*vmask(:,:,jk)
         !
      END DO
      !
#if defined key_agrif
      ! Save time integrated fluxes during child grid integration
      ! (used to update coarse grid transports at next time step)
      !
      IF ( (.NOT.Agrif_Root()).AND.(ln_bt_fw) ) THEN
         IF ( Agrif_NbStepint().EQ.0 ) THEN
            ub2_i_b(:,:) = 0.e0
            vb2_i_b(:,:) = 0.e0
         END IF
         !
         za1 = 1._wp / REAL(Agrif_rhot(), wp)
         ub2_i_b(:,:) = ub2_i_b(:,:) + za1 * ub2_b(:,:)
         vb2_i_b(:,:) = vb2_i_b(:,:) + za1 * vb2_b(:,:)
      ENDIF
      !
      !
#endif      
      !
      !                                   !* write time-spliting arrays in the restart
      IF(lrst_oce .AND.ln_bt_fw)   CALL ts_rst( kt, 'WRITE' )
      !
      CALL wrk_dealloc( jpi, jpj, zsshp2_e, zhdiv )
      CALL wrk_dealloc( jpi, jpj, zu_trd, zv_trd, zun_e, zvn_e )
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zu_sum, zv_sum, zssh_frc, zu_frc, zv_frc )
      CALL wrk_dealloc( jpi, jpj, zhup2_e, zhvp2_e, zhust_e, zhvst_e )
      CALL wrk_dealloc( jpi, jpj, zsshu_a, zsshv_a                                   )
      CALL wrk_dealloc( jpi, jpj, zhf )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_ts')
      !
   END SUBROUTINE dyn_spg_ts

   SUBROUTINE ts_wgt( ll_av, ll_fw, jpit, zwgt1, zwgt2)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_wgt  ***
      !!
      !! ** Purpose : Set time-splitting weights for temporal averaging (or not)
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in) ::   ll_av      ! temporal averaging=.true.
      LOGICAL, INTENT(in) ::   ll_fw      ! forward time splitting =.true.
      INTEGER, INTENT(inout) :: jpit      ! cycle length    
      REAL(wp), DIMENSION(3*nn_baro), INTENT(inout) ::   zwgt1, & ! Primary weights
                                                         zwgt2    ! Secondary weights
      
      INTEGER ::  jic, jn, ji                      ! temporary integers
      REAL(wp) :: za1, za2
      !!----------------------------------------------------------------------

      zwgt1(:) = 0._wp
      zwgt2(:) = 0._wp

      ! Set time index when averaged value is requested
      IF (ll_fw) THEN 
         jic = nn_baro
      ELSE
         jic = 2 * nn_baro
      ENDIF

      ! Set primary weights:
      IF (ll_av) THEN
           ! Define simple boxcar window for primary weights 
           ! (width = nn_baro, centered around jic)     
         SELECT CASE ( nn_bt_flt )
              CASE( 0 )  ! No averaging
                 zwgt1(jic) = 1._wp
                 jpit = jic

              CASE( 1 )  ! Boxcar, width = nn_baro
                 DO jn = 1, 3*nn_baro
                    za1 = ABS(float(jn-jic))/float(nn_baro) 
                    IF (za1 < 0.5_wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO

              CASE( 2 )  ! Boxcar, width = 2 * nn_baro
                 DO jn = 1, 3*nn_baro
                    za1 = ABS(float(jn-jic))/float(nn_baro) 
                    IF (za1 < 1._wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO
              CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt' )
         END SELECT

      ELSE ! No time averaging
         zwgt1(jic) = 1._wp
         jpit = jic
      ENDIF
    
      ! Set secondary weights
      DO jn = 1, jpit
        DO ji = jn, jpit
             zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
        END DO
      END DO

      ! Normalize weigths:
      za1 = 1._wp / SUM(zwgt1(1:jpit))
      za2 = 1._wp / SUM(zwgt2(1:jpit))
      DO jn = 1, jpit
        zwgt1(jn) = zwgt1(jn) * za1
        zwgt2(jn) = zwgt2(jn) * za2
      END DO
      !
   END SUBROUTINE ts_wgt

   SUBROUTINE ts_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_rst  ***
      !!
      !! ** Purpose : Read or write time-splitting arrays in restart file
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN
         CALL iom_get( numror, jpdom_autoglo, 'ub2_b'  , ub2_b  (:,:) )   
         CALL iom_get( numror, jpdom_autoglo, 'vb2_b'  , vb2_b  (:,:) ) 
         IF( .NOT.ln_bt_av ) THEN
            CALL iom_get( numror, jpdom_autoglo, 'sshbb_e'  , sshbb_e(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'ubb_e'    ,   ubb_e(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'vbb_e'    ,   vbb_e(:,:) )
            CALL iom_get( numror, jpdom_autoglo, 'sshb_e'   ,  sshb_e(:,:) ) 
            CALL iom_get( numror, jpdom_autoglo, 'ub_e'     ,    ub_e(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'vb_e'     ,    vb_e(:,:) )
         ENDIF
#if defined key_agrif
         ! Read time integrated fluxes
         IF ( .NOT.Agrif_Root() ) THEN
            CALL iom_get( numror, jpdom_autoglo, 'ub2_i_b'  , ub2_i_b(:,:) )   
            CALL iom_get( numror, jpdom_autoglo, 'vb2_i_b'  , vb2_i_b(:,:) )
         ENDIF
#endif
      !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         CALL iom_rstput( kt, nitrst, numrow, 'ub2_b'   , ub2_b  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'vb2_b'   , vb2_b  (:,:) )
         !
         IF (.NOT.ln_bt_av) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'sshbb_e'  , sshbb_e(:,:) ) 
            CALL iom_rstput( kt, nitrst, numrow, 'ubb_e'    ,   ubb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vbb_e'    ,   vbb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'sshb_e'   ,  sshb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'ub_e'     ,    ub_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb_e'     ,    vb_e(:,:) )
         ENDIF
#if defined key_agrif
         ! Save time integrated fluxes
         IF ( .NOT.Agrif_Root() ) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'ub2_i_b'  , ub2_i_b(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb2_i_b'  , vb2_i_b(:,:) )
         ENDIF
#endif
      ENDIF
      !
   END SUBROUTINE ts_rst

   SUBROUTINE dyn_spg_ts_init( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_spg_ts_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      !
      INTEGER  :: ji ,jj
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      REAL(wp) :: zxr2, zyr2, zcmax
      REAL(wp), POINTER, DIMENSION(:,:) :: zcu
      !!
      NAMELIST/namsplit/ ln_bt_fw, ln_bt_av, ln_bt_nn_auto, &
      &                  nn_baro, rn_bt_cmax, nn_bt_flt
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namsplit in reference namelist : time splitting parameters
      READ  ( numnam_ref, namsplit, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsplit in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsplit in configuration namelist : time splitting parameters
      READ  ( numnam_cfg, namsplit, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsplit in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsplit )
      !
      !         ! Max courant number for ext. grav. waves
      !
      CALL wrk_alloc( jpi, jpj, zcu )
      !
      IF (lk_vvl) THEN 
         DO jj = 1, jpj
            DO ji =1, jpi
               zxr2 = 1./(e1t(ji,jj)*e1t(ji,jj))
               zyr2 = 1./(e2t(ji,jj)*e2t(ji,jj))
               zcu(ji,jj) = sqrt(grav*ht_0(ji,jj)*(zxr2 + zyr2) )
            END DO
         END DO
      ELSE
         DO jj = 1, jpj
            DO ji =1, jpi
               zxr2 = 1./(e1t(ji,jj)*e1t(ji,jj))
               zyr2 = 1./(e2t(ji,jj)*e2t(ji,jj))
               zcu(ji,jj) = sqrt(grav*ht(ji,jj)*(zxr2 + zyr2) )
            END DO
         END DO
      ENDIF

      zcmax = MAXVAL(zcu(:,:))
      IF( lk_mpp )   CALL mpp_max( zcmax )

      ! Estimate number of iterations to satisfy a max courant number= rn_bt_cmax
      IF (ln_bt_nn_auto) nn_baro = CEILING( rdt / rn_bt_cmax * zcmax)
      
      rdtbt = rdt / FLOAT(nn_baro)
      zcmax = zcmax * rdtbt
							! Print results
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dyn_spg_ts : split-explicit free surface'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      IF( ln_bt_nn_auto ) THEN
         IF(lwp) WRITE(numout,*) '     ln_ts_nn_auto=.true. Automatically set nn_baro '
         IF(lwp) WRITE(numout,*) '     Max. courant number allowed: ', rn_bt_cmax
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_ts_nn_auto=.false.: Use nn_baro in namelist '
      ENDIF

      IF(ln_bt_av) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_av=.true.  => Time averaging over nn_baro time steps is on '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_av=.false. => No time averaging of barotropic variables '
      ENDIF
      !
      !
      IF(ln_bt_fw) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_fw=.true.  => Forward integration of barotropic variables '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_fw =.false.=> Centred integration of barotropic variables '
      ENDIF
      !
#if defined key_agrif
      ! Restrict the use of Agrif to the forward case only
      IF ((.NOT.ln_bt_fw ).AND.(.NOT.Agrif_Root())) CALL ctl_stop( 'AGRIF not implemented if ln_bt_fw=.FALSE.' )
#endif
      !
      IF(lwp) WRITE(numout,*)    '     Time filter choice, nn_bt_flt: ', nn_bt_flt
      SELECT CASE ( nn_bt_flt )
           CASE( 0 )  ;   IF(lwp) WRITE(numout,*) '           Dirac'
           CASE( 1 )  ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = nn_baro'
           CASE( 2 )  ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = 2*nn_baro' 
           CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt: should 0,1,2' )
      END SELECT
      !
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '     nn_baro = ', nn_baro
      IF(lwp) WRITE(numout,*) '     Barotropic time step [s] is :', rdtbt
      IF(lwp) WRITE(numout,*) '     Maximum Courant number is   :', zcmax
      !
      IF ((.NOT.ln_bt_av).AND.(.NOT.ln_bt_fw)) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: No time averaging => only forward integration is possible' )
      ENDIF
      IF ( zcmax>0.9_wp ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: Maximum Courant number is greater than 0.9: Inc. nn_baro !' )          
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zcu )
      !
   END SUBROUTINE dyn_spg_ts_init

#else
   !!---------------------------------------------------------------------------
   !!   Default case :   Empty module   No split explicit free surface
   !!---------------------------------------------------------------------------
CONTAINS
   INTEGER FUNCTION dyn_spg_ts_alloc()    ! Dummy function
      dyn_spg_ts_alloc = 0
   END FUNCTION dyn_spg_ts_alloc
   SUBROUTINE dyn_spg_ts( kt )            ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'dyn_spg_ts: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_ts
   SUBROUTINE ts_rst( kt, cdrw )          ! Empty routine
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      WRITE(*,*) 'ts_rst    : You should not have seen this print! error?', kt, cdrw
   END SUBROUTINE ts_rst  
   SUBROUTINE dyn_spg_ts_init( kt )       ! Empty routine
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      WRITE(*,*) 'dyn_spg_ts_init   : You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_ts_init
#endif
   
   !!======================================================================
END MODULE dynspg_ts



