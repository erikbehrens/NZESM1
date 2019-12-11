MODULE dynzdf_imp
   !!======================================================================
   !!                    ***  MODULE  dynzdf_imp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            8.0  !  1997-05  (G. Madec)  vertical component of isopycnal
   !!   NEMO     0.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.4  !  2012-01  (H. Liu) Semi-implicit bottom friction
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf_imp  : update the momentum trend with the vertical diffusion using a implicit time-stepping
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE domvvl          ! variable volume
   USE sbc_oce         ! surface boundary condition: ocean
   USE zdf_oce         ! ocean vertical physics
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE zdfbfr          ! Bottom friction setup
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE dynadv          ! dynamics: vector invariant versus flux form
   USE dynspg_oce, ONLY: lk_dynspg_ts

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf_imp   ! called by step.F90

   REAL(wp) ::  r_vvl     ! variable volume indicator, =1 if lk_vvl=T, =0 otherwise 

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zdf_imp( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_imp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!      and the surface forcing, and add it to the general trend of 
      !!      the momentum equations.
      !!
      !! ** Method  :   The vertical momentum mixing trend is given by :
      !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
      !!      backward time stepping
      !!      Surface boundary conditions: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) arrays with the after vertical diffusive mixing trend.
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in) ::  kt     ! ocean time-step index
      REAL(wp), INTENT(in) ::  p2dt   ! vertical profile of tracer time-step
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikbu, ikbv   ! local integers
      REAL(wp) ::   z1_p2dt, zcoef, zzwi, zzws, zrhs   ! local scalars
      REAL(wp) ::   ze3ua, ze3va
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwi, zwd, zws
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_imp')
      !
      CALL wrk_alloc( jpi,jpj,jpk, zwi, zwd, zws ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         !
         IF( lk_vvl ) THEN   ;    r_vvl = 1._wp       ! Variable volume indicator
         ELSE                ;    r_vvl = 0._wp       
         ENDIF
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      z1_p2dt = 1._wp / p2dt      ! inverse of the timestep

      ! 1. Apply semi-implicit bottom friction
      ! --------------------------------------
      ! Only needed for semi-implicit bottom friction setup. The explicit
      ! bottom friction has been included in "u(v)a" which act as the R.H.S
      ! column vector of the tri-diagonal matrix equation
      !

      IF( ln_bfrimp ) THEN
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = mbku(ji,jj)       ! ocean bottom level at u- and v-points 
               ikbv = mbkv(ji,jj)       ! (deepest ocean u- and v-points)
               avmu(ji,jj,ikbu+1) = -bfrua(ji,jj) * fse3uw(ji,jj,ikbu+1)
               avmv(ji,jj,ikbv+1) = -bfrva(ji,jj) * fse3vw(ji,jj,ikbv+1)
            END DO
         END DO
         IF ( ln_isfcav ) THEN
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ikbu = miku(ji,jj)       ! ocean top level at u- and v-points 
                  ikbv = mikv(ji,jj)       ! (first wet ocean u- and v-points)
                  IF (ikbu .GE. 2) avmu(ji,jj,ikbu) = -tfrua(ji,jj) * fse3uw(ji,jj,ikbu)
                  IF (ikbv .GE. 2) avmv(ji,jj,ikbv) = -tfrva(ji,jj) * fse3vw(ji,jj,ikbv)
               END DO
            END DO
         END IF
      ENDIF

#if defined key_dynspg_ts
      IF( ln_dynadv_vec .OR. .NOT. lk_vvl ) THEN      ! applied on velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = ( ub(:,:,jk) + p2dt * ua(:,:,jk) ) * umask(:,:,jk)
            va(:,:,jk) = ( vb(:,:,jk) + p2dt * va(:,:,jk) ) * vmask(:,:,jk)
         END DO
      ELSE                                            ! applied on thickness weighted velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = (          ub(:,:,jk) * fse3u_b(:,:,jk)      &
               &           + p2dt * ua(:,:,jk) * fse3u_n(:,:,jk)  )   &
               &                               / fse3u_a(:,:,jk) * umask(:,:,jk)
            va(:,:,jk) = (          vb(:,:,jk) * fse3v_b(:,:,jk)      &
               &           + p2dt * va(:,:,jk) * fse3v_n(:,:,jk)  )   &
               &                               / fse3v_a(:,:,jk) * vmask(:,:,jk)
         END DO
      ENDIF

      IF ( ln_bfrimp .AND.lk_dynspg_ts ) THEN
         ! remove barotropic velocities:
         DO jk = 1, jpkm1
            ua(:,:,jk) = (ua(:,:,jk) - ua_b(:,:)) * umask(:,:,jk)
            va(:,:,jk) = (va(:,:,jk) - va_b(:,:)) * vmask(:,:,jk)
         END DO
         ! Add bottom/top stress due to barotropic component only:
         DO jj = 2, jpjm1        
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points 
               ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
               ze3ua =  ( 1._wp - r_vvl ) * fse3u_n(ji,jj,ikbu) + r_vvl   * fse3u_a(ji,jj,ikbu)
               ze3va =  ( 1._wp - r_vvl ) * fse3v_n(ji,jj,ikbv) + r_vvl   * fse3v_a(ji,jj,ikbv)
               ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + p2dt * bfrua(ji,jj) * ua_b(ji,jj) / ze3ua
               va(ji,jj,ikbv) = va(ji,jj,ikbv) + p2dt * bfrva(ji,jj) * va_b(ji,jj) / ze3va
            END DO
         END DO
         IF ( ln_isfcav ) THEN
            DO jj = 2, jpjm1        
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ikbu = miku(ji,jj)         ! top ocean level at u- and v-points 
                  ikbv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
                  ze3ua =  ( 1._wp - r_vvl ) * fse3u_n(ji,jj,ikbu) + r_vvl   * fse3u_a(ji,jj,ikbu)
                  ze3va =  ( 1._wp - r_vvl ) * fse3v_n(ji,jj,ikbv) + r_vvl   * fse3v_a(ji,jj,ikbv)
                  ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + p2dt * tfrua(ji,jj) * ua_b(ji,jj) / ze3ua
                  va(ji,jj,ikbv) = va(ji,jj,ikbv) + p2dt * tfrva(ji,jj) * va_b(ji,jj) / ze3va
               END DO
            END DO
         END IF
      ENDIF
#endif

      ! 2. Vertical diffusion on u
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmu can take
      ! non zero value at the ocean bottom depending on the bottom friction used.
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ze3ua =  ( 1._wp - r_vvl ) * fse3u_n(ji,jj,jk) + r_vvl   * fse3u_a(ji,jj,jk)   ! after scale factor at T-point
               zcoef = - p2dt / ze3ua      
               zzwi          = zcoef * avmu  (ji,jj,jk  ) / fse3uw(ji,jj,jk  )
               zwi(ji,jj,jk) = zzwi  * wumask(ji,jj,jk  )
               zzws          = zcoef * avmu  (ji,jj,jk+1) / fse3uw(ji,jj,jk+1) 
               zws(ji,jj,jk) = zzws  * wumask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion starting from the first level
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and a lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (the after velocity) is in ua
      !-----------------------------------------------------------------------
      !
      !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_dynspg_ts
            ze3ua =  ( 1._wp - r_vvl ) * fse3u_n(ji,jj,1) + r_vvl   * fse3u_a(ji,jj,1) 
            ua(ji,jj,1) = ua(ji,jj,1) + p2dt * 0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
               &                                      / ( ze3ua * rau0 ) * umask(ji,jj,1) 
#else
            ua(ji,jj,1) = ub(ji,jj,1) &
               &                   + p2dt *(ua(ji,jj,1) +  0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
               &                                      / ( fse3u(ji,jj,1) * rau0     ) * umask(ji,jj,1) ) 
#endif
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
#if defined key_dynspg_ts
               zrhs = ua(ji,jj,jk)   ! zrhs=right hand side
#else
               zrhs = ub(ji,jj,jk) + p2dt * ua(ji,jj,jk)
#endif
               ua(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

#if ! defined key_dynspg_ts
      ! Normalization to obtain the general momentum trend ua
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - ub(ji,jj,jk) ) * z1_p2dt
            END DO
         END DO
      END DO
#endif

      ! 3. Vertical diffusion on v
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmv can take
      ! non zero value at the ocean bottom depending on the bottom friction used
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ze3va =  ( 1._wp - r_vvl ) * fse3v_n(ji,jj,jk)  + r_vvl * fse3v_a(ji,jj,jk)   ! after scale factor at T-point
               zcoef = - p2dt / ze3va
               zzwi          = zcoef * avmv (ji,jj,jk  ) / fse3vw(ji,jj,jk  )
               zwi(ji,jj,jk) =  zzwi * wvmask(ji,jj,jk)
               zzws          = zcoef * avmv (ji,jj,jk+1) / fse3vw(ji,jj,jk+1)
               zws(ji,jj,jk) =  zzws * wvmask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (after velocity) is in 2d array va
      !-----------------------------------------------------------------------
      !
      !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
      DO jk = 2, jpkm1        
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_dynspg_ts            
            ze3va =  ( 1._wp - r_vvl ) * fse3v_n(ji,jj,1) + r_vvl   * fse3v_a(ji,jj,1) 
            va(ji,jj,1) = va(ji,jj,1) + p2dt * 0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
               &                                      / ( ze3va * rau0 ) * vmask(ji,jj,1)
#else
            va(ji,jj,1) = vb(ji,jj,1) &
               &                   + p2dt *(va(ji,jj,1) +  0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
               &                                      / ( fse3v(ji,jj,1) * rau0     ) * vmask(ji,jj,1) )
#endif
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_dynspg_ts
               zrhs = va(ji,jj,jk)   ! zrhs=right hand side
#else
               zrhs = vb(ji,jj,jk) + p2dt * va(ji,jj,jk)
#endif
               va(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  third recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
            va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               va(ji,jj,jk) = ( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Normalization to obtain the general momentum trend va
#if ! defined key_dynspg_ts
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               va(ji,jj,jk) = ( va(ji,jj,jk) - vb(ji,jj,jk) ) * z1_p2dt
            END DO
         END DO
      END DO
#endif

      ! J. Chanut: Lines below are useless ?
      !! restore bottom layer avmu(v) 
      IF( ln_bfrimp ) THEN
        DO jj = 2, jpjm1
           DO ji = 2, jpim1
              ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points 
              ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
              avmu(ji,jj,ikbu+1) = 0.e0
              avmv(ji,jj,ikbv+1) = 0.e0
           END DO
        END DO
        IF (ln_isfcav) THEN
           DO jj = 2, jpjm1
              DO ji = 2, jpim1
                 ikbu = miku(ji,jj)         ! ocean top level at u- and v-points 
                 ikbv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
                 IF (ikbu > 1) avmu(ji,jj,ikbu) = 0.e0
                 IF (ikbv > 1) avmv(ji,jj,ikbv) = 0.e0
              END DO
           END DO
        END IF
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zwi, zwd, zws) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf_imp')
      !
   END SUBROUTINE dyn_zdf_imp

   !!==============================================================================
END MODULE dynzdf_imp
