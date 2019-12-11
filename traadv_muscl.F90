MODULE traadv_muscl
   !!======================================================================
   !!                       ***  MODULE  traadv_muscl  ***
   !! Ocean  tracers:  horizontal & vertical advective trend
   !!======================================================================
   !! History :       !  2000-06  (A.Estublier)  for passive tracers
   !!                 !  2001-08  (E.Durand, G.Madec)  adapted for T & S
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!            3.4  !  2012-06  (P. Oddo, M. Vichi) include the upstream where needed
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_muscl : update the tracer trend with the horizontal
   !!                   and vertical advection trends using MUSCL scheme
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE trc_oce        ! share passive tracers/Ocean variables
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! tracers trends manager
   USE dynspg_oce     ! choice/control of key cpp for surface pressure gradient
   USE sbcrnf         ! river runoffs
   USE diaptr         ! poleward transport diagnostics
   !
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing
   USE lbclnk         ! ocean lateral boundary condition (or mpp link) 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_muscl   ! routine called by traadv.F90
   
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   upsmsk   !: mixed upstream/centered scheme near some straits
   !                                                           !  and in closed seas (orca 2 and 4 configurations)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xind     !: mixed upstream/centered index
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_muscl( kt, kit000, cdtype, p2dt, pun, pvn, pwn,   &
      &                                                ptb, pta, kjpt, ld_msc_ups )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tra_adv_muscl  ***
      !!
      !! ** Purpose :   Compute the now trend due to total advection of T and 
      !!      S using a MUSCL scheme (Monotone Upstream-centered Scheme for
      !!      Conservation Laws) and add it to the general tracer trend.
      !!
      !! ** Method  : MUSCL scheme plus centered scheme at ocean boundaries
      !!
      !! ** Action  : - update (ta,sa) with the now advective tracer trends
      !!              - save trends 
      !!
      !! References : Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
      !!              IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      LOGICAL                              , INTENT(in   ) ::   ld_msc_ups      ! use upstream scheme within muscl
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt            ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb             ! before tracer field
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn            ! dummy loop indices
      INTEGER  ::   ierr                      ! local integer
      REAL(wp) ::   zu, z0u, zzwx, zw         ! local scalars
      REAL(wp) ::   zv, z0v, zzwy, z0w        !   -      -
      REAL(wp) ::   ztra, zbtr, zdt, zalpha   !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zslpx, zslpy   ! 3D workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zwx  , zwy     ! -      - 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_muscl')
      !
      ALLOCATE( zslpx(1:jpi, 1:jpj, 1:jpk) )
      ALLOCATE( zslpy(1:jpi, 1:jpj, 1:jpk) )
      ALLOCATE( zwx  (1:jpi, 1:jpj, 1:jpk) )
      ALLOCATE( zwy  (1:jpi, 1:jpj, 1:jpk) )
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv : MUSCL advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '        : mixed up-stream           ', ld_msc_ups
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         IF(lwp) WRITE(numout,*)
         !
         !
         IF( ld_msc_ups ) THEN
            IF( .NOT. ALLOCATED( upsmsk ) )  THEN
                ALLOCATE( upsmsk(jpi,jpj), STAT=ierr )
                IF( ierr /= 0 )   CALL ctl_stop('STOP', 'tra_adv_muscl: unable to allocate upsmsk array')
            ENDIF
            upsmsk(:,:) = 0._wp                             ! not upstream by default
         ENDIF

         IF( .NOT. ALLOCATED( xind ) ) THEN
             ALLOCATE( xind(jpi,jpj,jpk), STAT=ierr )
             IF( ierr /= 0 )   CALL ctl_stop('STOP', 'tra_adv_muscl: unable to allocate zind array')
         ENDIF
         !

         !
         ! Upstream / MUSCL scheme indicator
         ! ------------------------------------
!!gm  useless
         xind(:,:,:) = 1._wp                             ! set equal to 1 where up-stream is not needed
!!gm
         !
         IF( ld_msc_ups )  THEN
            DO jk = 1, jpkm1
               xind(:,:,jk) = 1._wp                              &                 ! =>1 where up-stream is not needed
                  &         - MAX ( rnfmsk(:,:) * rnfmsk_z(jk),  &                 ! =>0 near runoff mouths (& closed sea outflows)
                  &                 upsmsk(:,:)                ) * tmask(:,:,jk)   ! =>0 near some straits
            END DO
         ENDIF 
         !
      ENDIF 
      !      
      !                                                     ! ===========
      DO jn = 1, kjpt                                       ! tracer loop
         !                                                  ! ===========
         ! I. Horizontal advective fluxes
         ! ------------------------------
         ! first guess of the slopes
         zwx(:,:,jpk) = 0.e0   ;   zwy(:,:,jpk) = 0.e0        ! bottom values
         ! interior values
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1      
               DO ji = 1, fs_jpim1   ! vector opt.
                  zwx(ji,jj,jk) = umask(ji,jj,jk) * ( ptb(ji+1,jj,jk,jn) - ptb(ji,jj,jk,jn) )
                  zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( ptb(ji,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
           END DO
         END DO
         !
         CALL lbc_lnk( zwx, 'U', -1. )                        ! lateral boundary conditions on zwx, zwy   (changed sign)
         CALL lbc_lnk( zwy, 'V', -1. )
         !                                             !-- Slopes of tracer
         zslpx(:,:,jpk) = 0.e0   ;   zslpy(:,:,jpk) = 0.e0    ! bottom values
         DO jk = 1, jpkm1                                     ! interior values
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji-1,jj  ,jk) )   &
                     &            * ( 0.25 + SIGN( 0.25, zwx(ji,jj,jk) * zwx(ji-1,jj  ,jk) ) )
                  zslpy(ji,jj,jk) =                    ( zwy(ji,jj,jk) + zwy(ji  ,jj-1,jk) )   &
                     &            * ( 0.25 + SIGN( 0.25, zwy(ji,jj,jk) * zwy(ji  ,jj-1,jk) ) )
               END DO
            END DO
         END DO
         !
         DO jk = 1, jpkm1                                     ! Slopes limitation
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  zslpx(ji,jj,jk) = SIGN( 1., zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji  ,jj,jk) ),   &
                     &                                                 2.*ABS( zwx  (ji-1,jj,jk) ),   &
                     &                                                 2.*ABS( zwx  (ji  ,jj,jk) ) )
                  zslpy(ji,jj,jk) = SIGN( 1., zslpy(ji,jj,jk) ) * MIN(    ABS( zslpy(ji,jj  ,jk) ),   &
                     &                                                 2.*ABS( zwy  (ji,jj-1,jk) ),   &
                     &                                                 2.*ABS( zwy  (ji,jj  ,jk) ) )
               END DO
           END DO
         END DO             ! interior values

         !                                             !-- MUSCL horizontal advective fluxes
         DO jk = 1, jpkm1                                     ! interior values
            zdt  = p2dt(jk)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! MUSCL fluxes
                  z0u = SIGN( 0.5, pun(ji,jj,jk) )
                  zalpha = 0.5 - z0u
                  zu  = z0u - 0.5 * pun(ji,jj,jk) * zdt / ( e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk) )
                  zzwx = ptb(ji+1,jj,jk,jn) + xind(ji,jj,jk) * zu * zslpx(ji+1,jj,jk)
                  zzwy = ptb(ji  ,jj,jk,jn) + xind(ji,jj,jk) * zu * zslpx(ji  ,jj,jk)
                  zwx(ji,jj,jk) = pun(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
                  !
                  z0v = SIGN( 0.5, pvn(ji,jj,jk) )
                  zalpha = 0.5 - z0v
                  zv  = z0v - 0.5 * pvn(ji,jj,jk) * zdt / ( e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk) )
                  zzwx = ptb(ji,jj+1,jk,jn) + xind(ji,jj,jk) * zv * zslpy(ji,jj+1,jk)
                  zzwy = ptb(ji,jj  ,jk,jn) + xind(ji,jj,jk) * zv * zslpy(ji,jj  ,jk)
                  zwy(ji,jj,jk) = pvn(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
               END DO
            END DO
         END DO
         !                                                    ! lateral boundary conditions on zwx, zwy   (changed sign)
         CALL lbc_lnk( zwx, 'U', -1. )   ;   CALL lbc_lnk( zwy, 'V', -1. )
         !
         ! Tracer flux divergence at t-point added to the general trend
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ! horizontal advective trends
                  ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                  &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
                  ! add it to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
           END DO
         END DO        
         !                                 ! trend diagnostics (contribution of upstream fluxes)
         IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR.   &
            &( cdtype == 'TRC' .AND. l_trdtrc )      )  THEN
            CALL trd_tra( kt, cdtype, jn, jptra_xad, zwx, pun, ptb(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_yad, zwy, pvn, ptb(:,:,:,jn) )
         END IF
         !                                 ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( cdtype == 'TRA' .AND. ln_diaptr )  CALL dia_ptr_ohst_components( jn, 'adv', zwy(:,:,:)  )

         ! II. Vertical advective fluxes
         ! -----------------------------
         !                                             !-- first guess of the slopes
         zwx (:,:, 1 ) = 0.e0    ;    zwx (:,:,jpk) = 0.e0    ! surface & bottom boundary conditions
         DO jk = 2, jpkm1                                     ! interior values
            zwx(:,:,jk) = tmask(:,:,jk) * ( ptb(:,:,jk-1,jn) - ptb(:,:,jk,jn) )
         END DO

         !                                             !-- Slopes of tracer
         zslpx(:,:,1) = 0.e0                                  ! surface values
         DO jk = 2, jpkm1                                     ! interior value
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
                     &            * ( 0.25 + SIGN( 0.25, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )
               END DO
            END DO
         END DO
         !                                             !-- Slopes limitation
         DO jk = 2, jpkm1                                     ! interior values
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zslpx(ji,jj,jk) = SIGN( 1., zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji,jj,jk  ) ),   &
                     &                                                 2.*ABS( zwx  (ji,jj,jk+1) ),   &
                     &                                                 2.*ABS( zwx  (ji,jj,jk  ) )  )
               END DO
            END DO
         END DO
         !                                             !-- vertical advective flux
         !                                                    ! surface values  (bottom already set to zero)
         IF( lk_vvl ) THEN    ;   zwx(:,:, 1 ) = 0.e0                      !  variable volume
         ELSE                 ;   zwx(:,:, 1 ) = pwn(:,:,1) * ptb(:,:,1,jn)   ! linear free surface
         ENDIF 
         !
         DO jk = 1, jpkm1                                     ! interior values
            zdt  = p2dt(jk)
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3w(ji,jj,jk+1) )
                  z0w = SIGN( 0.5, pwn(ji,jj,jk+1) )
                  zalpha = 0.5 + z0w
                  zw  = z0w - 0.5 * pwn(ji,jj,jk+1) * zdt * zbtr 
                  zzwx = ptb(ji,jj,jk+1,jn) + xind(ji,jj,jk) * zw * zslpx(ji,jj,jk+1)
                  zzwy = ptb(ji,jj,jk  ,jn) + xind(ji,jj,jk) * zw * zslpx(ji,jj,jk  )
                  zwx(ji,jj,jk+1) = pwn(ji,jj,jk+1) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
               END DO 
            END DO
         END DO

         DO jk = 1, jpkm1                    ! Compute & add the vertical advective trend
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1. / ( e1e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ! vertical advective trends 
                  ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
                  ! add it to the general tracer trends
                  pta(ji,jj,jk,jn) =  pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO
         !                                 ! Save the vertical advective trends for diagnostic
         IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR.     &
            &( cdtype == 'TRC' .AND. l_trdtrc )      )   &
            CALL trd_tra( kt, cdtype, jn, jptra_zad, zwx, pwn, ptb(:,:,:,jn) )
         !
      END DO
      !
      DEALLOCATE( zslpx )
      DEALLOCATE( zslpy )
      DEALLOCATE( zwx   )
      DEALLOCATE( zwy   )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_muscl')
      !
   END SUBROUTINE tra_adv_muscl

   !!======================================================================
END MODULE traadv_muscl
