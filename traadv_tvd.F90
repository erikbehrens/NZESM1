MODULE traadv_tvd
   !!==============================================================================
   !!                       ***  MODULE  traadv_tvd  ***
   !! Ocean  tracers:  horizontal & vertical advective trend
   !!==============================================================================
   !! History :  OPA  !  1995-12  (L. Mortier)  Original code
   !!                 !  2000-01  (H. Loukos)  adapted to ORCA 
   !!                 !  2000-10  (MA Foujols E.Kestenare)  include file not routine
   !!                 !  2000-12  (E. Kestenare M. Levy)  fix bug in trtrd indexes
   !!                 !  2001-07  (E. Durand G. Madec)  adaptation to ORCA config
   !!            8.5  !  2002-06  (G. Madec)  F90: Free form and module
   !!    NEMO    1.0  !  2004-01  (A. de Miranda, G. Madec, J.M. Molines ): advective bbl
   !!            2.0  !  2008-04  (S. Cravatte) add the i-, j- & k- trends computation
   !!             -   !  2009-11  (V. Garnier) Surface pressure gradient organization
   !!            3.3  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_tvd   : update the tracer trend with the 3D advection trends using a TVD scheme
   !!   nonosc        : compute monotonic tracer fluxes by a non-oscillatory algorithm 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE trc_oce        ! share passive tracers/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! tracers trends
   USE dynspg_oce     ! choice/control of key cpp for surface pressure gradient
   USE diaptr         ! poleward transport diagnostics
   !
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary condition (or mpp link) 
   USE in_out_manager ! I/O manager
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_tvd        ! routine called by traadv.F90
   PUBLIC   tra_adv_tvd_zts    ! routine called by traadv.F90

   LOGICAL ::   l_trd   ! flag to compute trends
   LOGICAL ::   l_trans   ! flag to output vertically integrated transports

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_tvd ( kt, kit000, cdtype, p2dt, pun, pvn, pwn,      &
      &                                       ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_tvd  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!       tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
      !!       corrected flux (monotonic correction)
      !!       note: - this advection scheme needs a leap-frog time scheme
      !!
      !! ** Action : - update (pta) with the now advective tracer trends
      !!             - save the trends 
      !!----------------------------------------------------------------------
      USE oce     , ONLY:   zwx => ua        , zwy => va          ! (ua,va) used as workspace
      !
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt            ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      INTEGER  ::   ik  
      REAL(wp) ::   z2dtt, zbtr, ztra        ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj, zfp_wk   !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj, zfm_wk   !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwi, zwz
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrdx, ztrdy, ztrdz, zptry
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: z2d
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_tvd')
      !
      ALLOCATE(zwi(1:jpi, 1:jpj, 1:jpk))
      ALLOCATE(zwz(1:jpi, 1:jpj, 1:jpk))

      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_tvd : TVD advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         !
      ENDIF
      l_trd = .FALSE.
      l_trans = .FALSE.
      IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) ) l_trd = .TRUE.
      IF( cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") ) ) l_trans = .TRUE.
      !
      IF( l_trd .OR. l_trans )  THEN
         ALLOCATE(ztrdx(1:jpi, 1:jpj, 1:jpk))
         ALLOCATE(ztrdy(1:jpi, 1:jpj, 1:jpk))
         ALLOCATE(ztrdz(1:jpi, 1:jpj, 1:jpk))
         ztrdx(:,:,:) = 0.e0   ;    ztrdy(:,:,:) = 0.e0   ;   ztrdz(:,:,:) = 0.e0
         ALLOCATE(z2d(1:jpi, 1:jpj))
      ENDIF
      !
      IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN  
         ALLOCATE(zptry(1:jpi, 1:jpj, 1:jpk))
         zptry(:,:,:) = 0._wp
      ENDIF
      !
      zwi(:,:,:) = 0.e0 ; 
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         ! 1. Bottom and k=1 value : flux set to zero
         ! ----------------------------------
         zwx(:,:,jpk) = 0.e0    ;    zwz(:,:,jpk) = 0.e0
         zwy(:,:,jpk) = 0.e0    ;    zwi(:,:,jpk) = 0.e0
          
         zwz(:,:,1  ) = 0._wp
         ! 2. upstream advection with initial mass fluxes & intermediate update
         ! --------------------------------------------------------------------
         ! upstream tracer flux in the i and j direction
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! upstream scheme
                  zfp_ui = pun(ji,jj,jk) + ABS( pun(ji,jj,jk) )
                  zfm_ui = pun(ji,jj,jk) - ABS( pun(ji,jj,jk) )
                  zfp_vj = pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) )
                  zfm_vj = pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) )
                  zwx(ji,jj,jk) = 0.5 * ( zfp_ui * ptb(ji,jj,jk,jn) + zfm_ui * ptb(ji+1,jj  ,jk,jn) )
                  zwy(ji,jj,jk) = 0.5 * ( zfp_vj * ptb(ji,jj,jk,jn) + zfm_vj * ptb(ji  ,jj+1,jk,jn) )
               END DO
            END DO
         END DO

         ! upstream tracer flux in the k direction
         ! Interior value
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfp_wk = pwn(ji,jj,jk) + ABS( pwn(ji,jj,jk) )
                  zfm_wk = pwn(ji,jj,jk) - ABS( pwn(ji,jj,jk) )
                  zwz(ji,jj,jk) = 0.5 * ( zfp_wk * ptb(ji,jj,jk,jn) + zfm_wk * ptb(ji,jj,jk-1,jn) ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         ! Surface value
         IF( lk_vvl ) THEN   
            IF ( ln_isfcav ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = 0.e0          ! volume variable
                  END DO
               END DO
            ELSE
               zwz(:,:,1) = 0.e0          ! volume variable
            END IF
         ELSE                
            IF ( ln_isfcav ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptb(ji,jj,mikt(ji,jj),jn)   ! linear free surface 
                  END DO
               END DO   
            ELSE
               zwz(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)   ! linear free surface
            END IF
         ENDIF

         ! total advective trend
         DO jk = 1, jpkm1
            z2dtt = p2dt(jk)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! total intermediate advective trends
                  ztra = - (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                     &      + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                     &      + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) ) / e1e2t(ji,jj)
                  ! update and guess with monotonic sheme
                  pta(ji,jj,jk,jn) =                       pta(ji,jj,jk,jn) +         ztra   / fse3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                  zwi(ji,jj,jk)    = ( fse3t_b(ji,jj,jk) * ptb(ji,jj,jk,jn) + z2dtt * ztra ) / fse3t_a(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !                             ! Lateral boundary conditions on zwi  (unchanged sign)
         CALL lbc_lnk( zwi, 'T', 1. )  

         !                                 ! trend diagnostics (contribution of upstream fluxes)
         IF( l_trd .OR. l_trans )  THEN 
            ! store intermediate advective trends
            ztrdx(:,:,:) = zwx(:,:,:)   ;    ztrdy(:,:,:) = zwy(:,:,:)  ;   ztrdz(:,:,:) = zwz(:,:,:)
         END IF
         !                                 ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( cdtype == 'TRA' .AND. ln_diaptr )    zptry(:,:,:) = zwy(:,:,:) 

         ! 3. antidiffusive flux : high order minus low order
         ! --------------------------------------------------
         ! antidiffusive flux on i and j
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zwx(ji,jj,jk) = 0.5 * pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj,jk,jn) ) - zwx(ji,jj,jk)
                  zwy(ji,jj,jk) = 0.5 * pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj+1,jk,jn) ) - zwy(ji,jj,jk)
               END DO
            END DO
         END DO
      
         ! antidiffusive flux on k
         ! Interior value
         DO jk = 2, jpkm1                    
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zwz(ji,jj,jk) = 0.5 * pwn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj,jk-1,jn) ) - zwz(ji,jj,jk)
               END DO
            END DO
         END DO
         ! surface value
         IF ( ln_isfcav ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zwz(ji,jj,mikt(ji,jj)) = 0.e0
               END DO
            END DO
         ELSE
            zwz(:,:,1) = 0.e0
         END IF
         CALL lbc_lnk( zwx, 'U', -1. )   ;   CALL lbc_lnk( zwy, 'V', -1. )         ! Lateral bondary conditions
         CALL lbc_lnk( zwz, 'W',  1. )

         ! 4. monotonicity algorithm
         ! -------------------------
         CALL nonosc( ptb(:,:,:,jn), zwx, zwy, zwz, zwi, p2dt )


         ! 5. final trend with corrected fluxes
         ! ------------------------------------
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.  
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ! total advective trends
                  ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                     &             + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                     &             + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) )
                  ! add them to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra * tmask(ji,jj,jk)
               END DO
            END DO
         END DO

         !                                 ! trend diagnostics (contribution of upstream fluxes)
         IF( l_trd .OR. l_trans )  THEN 
            ztrdx(:,:,:) = ztrdx(:,:,:) + zwx(:,:,:)  ! <<< Add to previously computed
            ztrdy(:,:,:) = ztrdy(:,:,:) + zwy(:,:,:)  ! <<< Add to previously computed
            ztrdz(:,:,:) = ztrdz(:,:,:) + zwz(:,:,:)  ! <<< Add to previously computed
         ENDIF
         
         IF( l_trd ) THEN 
            CALL trd_tra( kt, cdtype, jn, jptra_xad, ztrdx, pun, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_yad, ztrdy, pvn, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_zad, ztrdz, pwn, ptn(:,:,:,jn) )
         END IF

         IF( l_trans .AND. jn==jp_tem ) THEN
            z2d(:,:) = 0._wp 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     z2d(ji,jj) = z2d(ji,jj) + ztrdx(ji,jj,jk) 
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z2d, 'U', -1. )
            CALL iom_put( "uadv_heattr", rau0_rcp * z2d )       ! heat transport in i-direction
              !
            z2d(:,:) = 0._wp 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     z2d(ji,jj) = z2d(ji,jj) + ztrdy(ji,jj,jk) 
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z2d, 'V', -1. )
            CALL iom_put( "vadv_heattr", rau0_rcp * z2d )       ! heat transport in j-direction
         ENDIF
         ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN  
            zptry(:,:,:) = zptry(:,:,:) + zwy(:,:,:)  ! <<< Add to previously computed
            CALL dia_ptr_ohst_components( jn, 'adv', zptry(:,:,:) )
         ENDIF
         !
      END DO
      !
      DEALLOCATE( zwi )
      DEALLOCATE( zwz )
      IF( l_trd .OR. l_trans )  THEN 
         DEALLOCATE( ztrdx )
         DEALLOCATE( ztrdy )
         DEALLOCATE( ztrdz )
         DEALLOCATE( z2d )
      ENDIF
      IF( cdtype == 'TRA' .AND. ln_diaptr ) DEALLOCATE( zptry )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_tvd')
      !
   END SUBROUTINE tra_adv_tvd

   SUBROUTINE tra_adv_tvd_zts ( kt, kit000, cdtype, p2dt, pun, pvn, pwn,      &
      &                                       ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_tvd_zts  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!       tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   TVD ZTS scheme, i.e. 2nd order centered scheme with
      !!       corrected flux (monotonic correction). This version use sub-
      !!       timestepping for the vertical advection which increases stability
      !!       when vertical metrics are small.
      !!       note: - this advection scheme needs a leap-frog time scheme
      !!
      !! ** Action : - update (pta) with the now advective tracer trends
      !!             - save the trends 
      !!----------------------------------------------------------------------
      USE oce     , ONLY:   zwx => ua        , zwy => va          ! (ua,va) used as workspace
      !
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt            ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      REAL(wp), DIMENSION( jpk )                           ::   zts             ! length of sub-timestep for vertical advection
      REAL(wp), DIMENSION( jpk )                           ::   zr_p2dt         ! reciprocal of tracer timestep
      INTEGER  ::   ji, jj, jk, jl, jn       ! dummy loop indices  
      INTEGER  ::   jnzts = 5       ! number of sub-timesteps for vertical advection
      INTEGER  ::   jtb, jtn, jta   ! sub timestep pointers for leap-frog/euler forward steps
      INTEGER  ::   jtaken          ! toggle for collecting appropriate fluxes from sub timesteps
      REAL(wp) ::   z_rzts          ! Fractional length of Euler forward sub-timestep for vertical advection
      REAL(wp) ::   z2dtt, zbtr, ztra        ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj, zfp_wk   !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj, zfm_wk   !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zwx_sav , zwy_sav
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwi, zwz, zhdiv, zwz_sav, zwzts
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrdx, ztrdy, ztrdz
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zptry
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: ztrs
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_tvd_zts')
      !
      ALLOCATE(zwx_sav(1:jpi, 1:jpj))
      ALLOCATE(zwy_sav(1:jpi, 1:jpj))
      ALLOCATE(zwi(1:jpi, 1:jpj, 1:jpk))
      ALLOCATE(zwz(1:jpi, 1:jpj, 1:jpk))        
      ALLOCATE(zhdiv(1:jpi, 1:jpj, 1:jpk))       
      ALLOCATE(zwz_sav(1:jpi, 1:jpj, 1:jpk))       
      ALLOCATE(zwzts(1:jpi, 1:jpj, 1:jpk)) 
      ALLOCATE(ztrs(1:jpi, 1:jpj, 1:jpk, 1:kjpt+1))
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_tvd_zts : TVD ZTS advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      l_trd = .FALSE.
      IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) ) l_trd = .TRUE.
      !
      IF( l_trd )  THEN
         ALLOCATE(ztrdx(1:jpi, 1:jpj, 1:jpk))       
         ALLOCATE(ztrdy(1:jpi, 1:jpj, 1:jpk))       
         ALLOCATE(ztrdz(1:jpi, 1:jpj, 1:jpk))       
         ztrdx(:,:,:) = 0._wp  ;    ztrdy(:,:,:) = 0._wp  ;   ztrdz(:,:,:) = 0._wp
      ENDIF
      !
      IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN  
         ALLOCATE(zptry(1:jpi, 1:jpj, 1:jpk))       
         zptry(:,:,:) = 0._wp
      ENDIF
      !
      zwi(:,:,:) = 0._wp
      z_rzts = 1._wp / REAL( jnzts, wp )
      zr_p2dt(:) = 1._wp / p2dt(:)
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         ! 1. Bottom value : flux set to zero
         ! ----------------------------------
         zwx(:,:,jpk) = 0._wp   ;    zwz(:,:,jpk) = 0._wp
         zwy(:,:,jpk) = 0._wp   ;    zwi(:,:,jpk) = 0._wp

         ! 2. upstream advection with initial mass fluxes & intermediate update
         ! --------------------------------------------------------------------
         ! upstream tracer flux in the i and j direction
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! upstream scheme
                  zfp_ui = pun(ji,jj,jk) + ABS( pun(ji,jj,jk) )
                  zfm_ui = pun(ji,jj,jk) - ABS( pun(ji,jj,jk) )
                  zfp_vj = pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) )
                  zfm_vj = pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) )
                  zwx(ji,jj,jk) = 0.5_wp * ( zfp_ui * ptb(ji,jj,jk,jn) + zfm_ui * ptb(ji+1,jj  ,jk,jn) )
                  zwy(ji,jj,jk) = 0.5_wp * ( zfp_vj * ptb(ji,jj,jk,jn) + zfm_vj * ptb(ji  ,jj+1,jk,jn) )
               END DO
            END DO
         END DO

         ! upstream tracer flux in the k direction
         ! Interior value
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfp_wk = pwn(ji,jj,jk) + ABS( pwn(ji,jj,jk) )
                  zfm_wk = pwn(ji,jj,jk) - ABS( pwn(ji,jj,jk) )
                  zwz(ji,jj,jk) = 0.5_wp * ( zfp_wk * ptb(ji,jj,jk,jn) + zfm_wk * ptb(ji,jj,jk-1,jn) )
               END DO
            END DO
         END DO
         ! Surface value
         IF( lk_vvl ) THEN
            IF ( ln_isfcav ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = 0.e0          ! volume variable +    isf
                  END DO
               END DO
            ELSE
               zwz(:,:,1) = 0.e0                              ! volume variable + no isf
            END IF
         ELSE
            IF ( ln_isfcav ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptb(ji,jj,mikt(ji,jj),jn)   ! linear free surface +    isf
                  END DO
               END DO
            ELSE
               zwz(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)                                               ! linear free surface + no isf
            END IF
         ENDIF

         ! total advective trend
         DO jk = 1, jpkm1
            z2dtt = p2dt(jk)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! total intermediate advective trends
                  ztra = - (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                     &      + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                     &      + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) ) / e1e2t(ji,jj)
                  ! update and guess with monotonic sheme
                  pta(ji,jj,jk,jn) =                       pta(ji,jj,jk,jn) +         ztra   / fse3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                  zwi(ji,jj,jk)    = ( fse3t_b(ji,jj,jk) * ptb(ji,jj,jk,jn) + z2dtt * ztra ) / fse3t_a(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !                             ! Lateral boundary conditions on zwi  (unchanged sign)
         CALL lbc_lnk( zwi, 'T', 1. )  

         !                                 ! trend diagnostics (contribution of upstream fluxes)
         IF( l_trd )  THEN 
            ! store intermediate advective trends
            ztrdx(:,:,:) = zwx(:,:,:)   ;    ztrdy(:,:,:) = zwy(:,:,:)  ;   ztrdz(:,:,:) = zwz(:,:,:)
         END IF
         !                                 ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( cdtype == 'TRA' .AND. ln_diaptr )  zptry(:,:,:) = zwy(:,:,:)

         ! 3. antidiffusive flux : high order minus low order
         ! --------------------------------------------------
         ! antidiffusive flux on i and j
         !
         DO jk = 1, jpkm1
            !
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zwx_sav(ji,jj) = zwx(ji,jj,jk)
                  zwy_sav(ji,jj) = zwy(ji,jj,jk)

                  zwx(ji,jj,jk) = 0.5_wp * pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj,jk,jn) )
                  zwy(ji,jj,jk) = 0.5_wp * pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj+1,jk,jn) )
               END DO
            END DO

            DO jj = 2, jpjm1         ! partial horizontal divergence
               DO ji = fs_2, fs_jpim1
                  zhdiv(ji,jj,jk) = (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk)   &
                     &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk)  )
               END DO
            END DO

            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zwx(ji,jj,jk) = zwx(ji,jj,jk)  - zwx_sav(ji,jj)
                  zwy(ji,jj,jk) = zwy(ji,jj,jk)  - zwy_sav(ji,jj)
               END DO
            END DO
         END DO
      
         ! antidiffusive flux on k
         zwz(:,:,1) = 0._wp        ! Surface value
         zwz_sav(:,:,:) = zwz(:,:,:)
         !
         ztrs(:,:,:,1) = ptb(:,:,:,jn)
         ztrs(:,:,1,2) = ptb(:,:,1,jn)
         ztrs(:,:,1,3) = ptb(:,:,1,jn)
         zwzts(:,:,:) = 0._wp

         DO jl = 1, jnzts                   ! Start of sub timestepping loop

            IF( jl == 1 ) THEN              ! Euler forward to kick things off
              jtb = 1   ;   jtn = 1   ;   jta = 2
              zts(:) = p2dt(:) * z_rzts
              jtaken = MOD( jnzts + 1 , 2)  ! Toggle to collect every second flux
                                            ! starting at jl =1 if jnzts is odd; 
                                            ! starting at jl =2 otherwise
            ELSEIF( jl == 2 ) THEN          ! First leapfrog step
              jtb = 1   ;   jtn = 2   ;   jta = 3
              zts(:) = 2._wp * p2dt(:) * z_rzts
            ELSE                            ! Shuffle pointers for subsequent leapfrog steps
              jtb = MOD(jtb,3) + 1
              jtn = MOD(jtn,3) + 1
              jta = MOD(jta,3) + 1
            ENDIF
            DO jk = 2, jpkm1          ! Interior value
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zwz(ji,jj,jk) = 0.5_wp * pwn(ji,jj,jk) * ( ztrs(ji,jj,jk,jtn) + ztrs(ji,jj,jk-1,jtn) )
                     IF( jtaken == 0 ) zwzts(ji,jj,jk) = zwzts(ji,jj,jk) + zwz(ji,jj,jk)*zts(jk)           ! Accumulate time-weighted vertcal flux
                  END DO
               END DO
            END DO

            jtaken = MOD( jtaken + 1 , 2 )

            DO jk = 2, jpkm1          ! Interior value
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zbtr = 1._wp / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                     ! total advective trends
                     ztra = - zbtr * (  zhdiv(ji,jj,jk) + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) )
                     ztrs(ji,jj,jk,jta) = ztrs(ji,jj,jk,jtb) + zts(jk) * ztra
                  END DO
               END DO
            END DO

         END DO

         DO jk = 2, jpkm1          ! Anti-diffusive vertical flux using average flux from the sub-timestepping
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zwz(ji,jj,jk) = zwzts(ji,jj,jk) * zr_p2dt(jk) - zwz_sav(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( zwx, 'U', -1. )   ;   CALL lbc_lnk( zwy, 'V', -1. )         ! Lateral bondary conditions
         CALL lbc_lnk( zwz, 'W',  1. )

         ! 4. monotonicity algorithm
         ! -------------------------
         CALL nonosc( ptb(:,:,:,jn), zwx, zwy, zwz, zwi, p2dt )


         ! 5. final trend with corrected fluxes
         ! ------------------------------------
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.  
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ! total advective trends
                  ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                     &             + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                     &             + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) )
                  ! add them to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO

         !                                 ! trend diagnostics (contribution of upstream fluxes)
         IF( l_trd )  THEN 
            ztrdx(:,:,:) = ztrdx(:,:,:) + zwx(:,:,:)  ! <<< Add to previously computed
            ztrdy(:,:,:) = ztrdy(:,:,:) + zwy(:,:,:)  ! <<< Add to previously computed
            ztrdz(:,:,:) = ztrdz(:,:,:) + zwz(:,:,:)  ! <<< Add to previously computed
            
            CALL trd_tra( kt, cdtype, jn, jptra_xad, ztrdx, pun, ptn(:,:,:,jn) )   
            CALL trd_tra( kt, cdtype, jn, jptra_yad, ztrdy, pvn, ptn(:,:,:,jn) )  
            CALL trd_tra( kt, cdtype, jn, jptra_zad, ztrdz, pwn, ptn(:,:,:,jn) ) 
         END IF
         !                                 ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN  
            zptry(:,:,:) = zptry(:,:,:) + zwy(:,:,:) 
            CALL dia_ptr_ohst_components( jn, 'adv', zptry(:,:,:) )
         ENDIF
         !
      END DO
      !
      DEALLOCATE(zwi) 
      DEALLOCATE(zwz) 
      DEALLOCATE(zhdiv) 
      DEALLOCATE(zwz_sav) 
      DEALLOCATE(zwzts)
      DEALLOCATE(ztrs )
      DEALLOCATE(zwx_sav) 
      DEALLOCATE(zwy_sav )

      IF( l_trd )  THEN
          DEALLOCATE(ztrdx) 
          DEALLOCATE(ztrdy) 
          DEALLOCATE(ztrdz)
      END IF

      IF( cdtype == 'TRA' .AND. ln_diaptr ) DEALLOCATE(zptry )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_tvd_zts')
      !
   END SUBROUTINE tra_adv_tvd_zts


   SUBROUTINE nonosc( pbef, paa, pbb, pcc, paft, p2dt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc  ***
      !!     
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpk)         , INTENT(in   ) ::   p2dt            ! vertical profile of tracer time-step
      REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT(in   ) ::   pbef, paft      ! before & after field
      REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT(inout) ::   paa, pbb, pcc   ! monotonic fluxes in the 3 directions
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikm1         ! local integer
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn, z2dtt   ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zup, zdo            !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zbetup, zbetdo, zbup, zbdo
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('nonosc')
      !
      ALLOCATE(zbetup(1:jpi, 1:jpj, 1:jpk))
      ALLOCATE(zbetdo(1:jpi, 1:jpj, 1:jpk))
      ALLOCATE(zbup(1:jpi, 1:jpj, 1:jpk))
      ALLOCATE(zbdo(1:jpi, 1:jpj, 1:jpk))
      !
      zbig  = 1.e+40_wp
      zrtrn = 1.e-15_wp
      zbetup(:,:,:) = 0._wp   ;   zbetdo(:,:,:) = 0._wp

      ! Search local extrema
      ! --------------------
      ! max/min of pbef & paft with large negative/positive value (-/+zbig) inside land
      zbup = MAX( pbef * tmask - zbig * ( 1._wp - tmask ),   &
         &        paft * tmask - zbig * ( 1._wp - tmask )  )
      zbdo = MIN( pbef * tmask + zbig * ( 1._wp - tmask ),   &
         &        paft * tmask + zbig * ( 1._wp - tmask )  )

      DO jk = 1, jpkm1
         ikm1 = MAX(jk-1,1)
         z2dtt = p2dt(jk)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.

               ! search maximum in neighbourhood
               zup = MAX(  zbup(ji  ,jj  ,jk  ),   &
                  &        zbup(ji-1,jj  ,jk  ), zbup(ji+1,jj  ,jk  ),   &
                  &        zbup(ji  ,jj-1,jk  ), zbup(ji  ,jj+1,jk  ),   &
                  &        zbup(ji  ,jj  ,ikm1), zbup(ji  ,jj  ,jk+1)  )

               ! search minimum in neighbourhood
               zdo = MIN(  zbdo(ji  ,jj  ,jk  ),   &
                  &        zbdo(ji-1,jj  ,jk  ), zbdo(ji+1,jj  ,jk  ),   &
                  &        zbdo(ji  ,jj-1,jk  ), zbdo(ji  ,jj+1,jk  ),   &
                  &        zbdo(ji  ,jj  ,ikm1), zbdo(ji  ,jj  ,jk+1)  )

               ! positive part of the flux
               zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   &
                  & + MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )

               ! negative part of the flux
               zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )   &
                  & + MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )

               ! up & down beta terms
               zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) / z2dtt
               zbetup(ji,jj,jk) = ( zup            - paft(ji,jj,jk) ) / ( zpos + zrtrn ) * zbt
               zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zdo            ) / ( zneg + zrtrn ) * zbt
            END DO
         END DO
      END DO
      CALL lbc_lnk( zbetup, 'T', 1. )   ;   CALL lbc_lnk( zbetdo, 'T', 1. )   ! lateral boundary cond. (unchanged sign)

      ! 3. monotonic flux in the i & j direction (paa & pbb)
      ! ----------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zau = MIN( 1._wp, zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
               zbu = MIN( 1._wp, zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
               zcu =       ( 0.5  + SIGN( 0.5 , paa(ji,jj,jk) ) )
               paa(ji,jj,jk) = paa(ji,jj,jk) * ( zcu * zau + ( 1._wp - zcu) * zbu )

               zav = MIN( 1._wp, zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
               zbv = MIN( 1._wp, zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
               zcv =       ( 0.5  + SIGN( 0.5 , pbb(ji,jj,jk) ) )
               pbb(ji,jj,jk) = pbb(ji,jj,jk) * ( zcv * zav + ( 1._wp - zcv) * zbv )

      ! monotonic flux in the k direction, i.e. pcc
      ! -------------------------------------------
               za = MIN( 1., zbetdo(ji,jj,jk+1), zbetup(ji,jj,jk) )
               zb = MIN( 1., zbetup(ji,jj,jk+1), zbetdo(ji,jj,jk) )
               zc =       ( 0.5  + SIGN( 0.5 , pcc(ji,jj,jk+1) ) )
               pcc(ji,jj,jk+1) = pcc(ji,jj,jk+1) * ( zc * za + ( 1._wp - zc) * zb )
            END DO
         END DO
      END DO
      CALL lbc_lnk( paa, 'U', -1. )   ;   CALL lbc_lnk( pbb, 'V', -1. )   ! lateral boundary condition (changed sign)
      !
      DEALLOCATE(zbetup)
      DEALLOCATE(zbetdo) 
      DEALLOCATE(zbup)
      DEALLOCATE(zbdo)
      !
      IF( nn_timing == 1 )  CALL timing_stop('nonosc')
      !
   END SUBROUTINE nonosc

   !!======================================================================
END MODULE traadv_tvd
