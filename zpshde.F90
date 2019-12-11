MODULE zpshde
   !!======================================================================
   !!                       ***  MODULE zpshde   ***
   !! z-coordinate + partial step : Horizontal Derivative at ocean bottom level
   !!======================================================================
   !! History :  OPA  !  2002-04  (A. Bozec)  Original code
   !!   NEMO     1.0  !  2002-08  (G. Madec E. Durand)  Optimization and Free form
   !!             -   !  2004-03  (C. Ethe)  adapted for passive tracers
   !!            3.3  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!            3.6  !  2014-11  (P. Mathiot) Add zps_hde_isf (needed to open a cavity)
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !!   zps_hde      :  Horizontal DErivative of T, S and rd at the last
   !!                   ocean level (Z-coord. with Partial Steps)
   !!----------------------------------------------------------------------
   USE oce             ! ocean: dynamics and tracers variables
   USE dom_oce         ! domain: ocean variables
   USE phycst          ! physical constants
   USE eosbn2          ! ocean equation of state
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zps_hde     ! routine called by step.F90
   PUBLIC   zps_hde_isf ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zps_hde( kt, kjpt, pta, pgtu, pgtv,   &
      &                          prd, pgru, pgrv    )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE zps_hde  ***
      !!                    
      !! ** Purpose :   Compute the horizontal derivative of T, S and rho
      !!      at u- and v-points with a linear interpolation for z-coordinate
      !!      with partial steps.
      !!
      !! ** Method  :   In z-coord with partial steps, scale factors on last 
      !!      levels are different for each grid point, so that T, S and rd 
      !!      points are not at the same depth as in z-coord. To have horizontal
      !!      gradients again, we interpolate T and S at the good depth : 
      !!      Linear interpolation of T, S   
      !!         Computation of di(tb) and dj(tb) by vertical interpolation:
      !!          di(t) = t~ - t(i,j,k) or t(i+1,j,k) - t~
      !!          dj(t) = t~ - t(i,j,k) or t(i,j+1,k) - t~
      !!         This formulation computes the two cases:
      !!                 CASE 1                   CASE 2  
      !!         k-1  ___ ___________   k-1   ___ ___________
      !!                    Ti  T~                  T~  Ti+1
      !!                  _____                        _____
      !!         k        |   |Ti+1     k           Ti |   |
      !!                  |   |____                ____|   |
      !!              ___ |   |   |           ___  |   |   |
      !!                  
      !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
      !!          t~ = t(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(Ti+1)/e3w(i+1)
      !!        ( t~ = t(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(Tj+1)/e3w(j+1)  )
      !!          or
      !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
      !!          t~ = t(i,j,k) + (e3w(i) - e3w(i+1)) * dk(Ti)/e3w(i )
      !!        ( t~ = t(i,j,k) + (e3w(j) - e3w(j+1)) * dk(Tj)/e3w(j ) )
      !!          Idem for di(s) and dj(s)          
      !!
      !!      For rho, we call eos which will compute rd~(t~,s~) at the right
      !!      depth zh from interpolated T and S for the different formulations
      !!      of the equation of state (eos).
      !!      Gradient formulation for rho :
      !!          di(rho) = rd~ - rd(i,j,k)   or   rd(i+1,j,k) - rd~
      !!
      !! ** Action  : compute for top interfaces
      !!              - pgtu, pgtv: horizontal gradient of tracer at u- & v-points
      !!              - pgru, pgrv: horizontal gradient of rho (if present) at u- & v-points
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   )           ::  kt          ! ocean time-step index
      INTEGER                              , INTENT(in   )           ::  kjpt        ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   )           ::  pta         ! 4D tracers fields
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(  out)           ::  pgtu, pgtv  ! hor. grad. of ptra at u- & v-pts 
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ), OPTIONAL ::  prd         ! 3D density anomaly fields
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pgru, pgrv  ! hor. grad of prd at u- & v-pts (bottom)
      !
      INTEGER  ::   ji, jj, jn      ! Dummy loop indices
      INTEGER  ::   iku, ikv, ikum1, ikvm1   ! partial step level (ocean bottom level) at u- and v-points
      REAL(wp) ::  ze3wu, ze3wv, zmaxu, zmaxv  ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj)      ::  zri, zrj, zhi, zhj   ! NB: 3rd dim=1 to use eos
      REAL(wp), DIMENSION(jpi,jpj,kjpt) ::  zti, ztj             ! 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'zps_hde')
      !
      pgtu(:,:,:)=0.0_wp ; pgtv(:,:,:)=0.0_wp ;
      zti (:,:,:)=0.0_wp ; ztj (:,:,:)=0.0_wp ;
      zhi (:,:  )=0.0_wp ; zhj (:,:  )=0.0_wp ;
      !
      DO jn = 1, kjpt      !==   Interpolation of tracers at the last ocean level   ==!
         !
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)   ;   ikum1 = MAX( iku - 1 , 1 )    ! last and before last ocean level at u- & v-points
               ikv = mbkv(ji,jj)   ;   ikvm1 = MAX( ikv - 1 , 1 )    ! if level first is a p-step, ik.m1=1
               ze3wu = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
               ze3wv = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
               !
               ! i- direction
               IF( ze3wu >= 0._wp ) THEN      ! case 1
                  zmaxu =  ze3wu / fse3w(ji+1,jj,iku)
                  ! interpolated values of tracers
                  zti (ji,jj,jn) = pta(ji+1,jj,iku,jn) + zmaxu * ( pta(ji+1,jj,ikum1,jn) - pta(ji+1,jj,iku,jn) )
                  ! gradient of  tracers
                  pgtu(ji,jj,jn) = umask(ji,jj,1) * ( zti(ji,jj,jn) - pta(ji,jj,iku,jn) )
               ELSE                           ! case 2
                  zmaxu = -ze3wu / fse3w(ji,jj,iku)
                  ! interpolated values of tracers
                  zti (ji,jj,jn) = pta(ji,jj,iku,jn) + zmaxu * ( pta(ji,jj,ikum1,jn) - pta(ji,jj,iku,jn) )
                  ! gradient of tracers
                  pgtu(ji,jj,jn) = umask(ji,jj,1) * ( pta(ji+1,jj,iku,jn) - zti(ji,jj,jn) )
               ENDIF
               !
               ! j- direction
               IF( ze3wv >= 0._wp ) THEN      ! case 1
                  zmaxv =  ze3wv / fse3w(ji,jj+1,ikv)
                  ! interpolated values of tracers
                  ztj (ji,jj,jn) = pta(ji,jj+1,ikv,jn) + zmaxv * ( pta(ji,jj+1,ikvm1,jn) - pta(ji,jj+1,ikv,jn) )
                  ! gradient of tracers
                  pgtv(ji,jj,jn) = vmask(ji,jj,1) * ( ztj(ji,jj,jn) - pta(ji,jj,ikv,jn) )
               ELSE                           ! case 2
                  zmaxv =  -ze3wv / fse3w(ji,jj,ikv)
                  ! interpolated values of tracers
                  ztj (ji,jj,jn) = pta(ji,jj,ikv,jn) + zmaxv * ( pta(ji,jj,ikvm1,jn) - pta(ji,jj,ikv,jn) )
                  ! gradient of tracers
                  pgtv(ji,jj,jn) = vmask(ji,jj,1) * ( pta(ji,jj+1,ikv,jn) - ztj(ji,jj,jn) )
               ENDIF
            END DO
         END DO
         CALL lbc_lnk( pgtu(:,:,jn), 'U', -1. )   ;   CALL lbc_lnk( pgtv(:,:,jn), 'V', -1. )   ! Lateral boundary cond.
         !
      END DO

      ! horizontal derivative of density anomalies (rd)
      IF( PRESENT( prd ) ) THEN         ! depth of the partial step level
         pgru(:,:)=0.0_wp   ; pgrv(:,:)=0.0_wp ; 
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)
               ikv = mbkv(ji,jj)
               ze3wu  = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
               ze3wv  = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
               IF( ze3wu >= 0._wp ) THEN   ;   zhi(ji,jj) = fsdept(ji  ,jj,iku)     ! i-direction: case 1
               ELSE                        ;   zhi(ji,jj) = fsdept(ji+1,jj,iku)     ! -     -      case 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN   ;   zhj(ji,jj) = fsdept(ji,jj  ,ikv)     ! j-direction: case 1
               ELSE                        ;   zhj(ji,jj) = fsdept(ji,jj+1,ikv)     ! -     -      case 2
               ENDIF
            END DO
         END DO

         ! Compute interpolated rd from zti, ztj for the 2 cases at the depth of the partial
         ! step and store it in  zri, zrj for each  case
         CALL eos( zti, zhi, zri )  
         CALL eos( ztj, zhj, zrj )

         ! Gradient of density at the last level 
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)
               ikv = mbkv(ji,jj)
               ze3wu  = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
               ze3wv  = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
               IF( ze3wu >= 0._wp ) THEN   ;   pgru(ji,jj) = umask(ji,jj,1) * ( zri(ji  ,jj    ) - prd(ji,jj,iku) )   ! i: 1
               ELSE                        ;   pgru(ji,jj) = umask(ji,jj,1) * ( prd(ji+1,jj,iku) - zri(ji,jj    ) )   ! i: 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN   ;   pgrv(ji,jj) = vmask(ji,jj,1) * ( zrj(ji,jj      ) - prd(ji,jj,ikv) )   ! j: 1
               ELSE                        ;   pgrv(ji,jj) = vmask(ji,jj,1) * ( prd(ji,jj+1,ikv) - zrj(ji,jj    ) )   ! j: 2
               ENDIF
            END DO
         END DO
         CALL lbc_lnk( pgru , 'U', -1. )   ;   CALL lbc_lnk( pgrv , 'V', -1. )   ! Lateral boundary conditions
         !
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'zps_hde')
      !
   END SUBROUTINE zps_hde
   !
   SUBROUTINE zps_hde_isf( kt, kjpt, pta, pgtu, pgtv,   &
      &                          prd, pgru, pgrv, pmru, pmrv, pgzu, pgzv, pge3ru, pge3rv,  &
      &                   pgtui, pgtvi, pgrui, pgrvi, pmrui, pmrvi, pgzui, pgzvi, pge3rui, pge3rvi )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE zps_hde  ***
      !!                    
      !! ** Purpose :   Compute the horizontal derivative of T, S and rho
      !!      at u- and v-points with a linear interpolation for z-coordinate
      !!      with partial steps.
      !!
      !! ** Method  :   In z-coord with partial steps, scale factors on last 
      !!      levels are different for each grid point, so that T, S and rd 
      !!      points are not at the same depth as in z-coord. To have horizontal
      !!      gradients again, we interpolate T and S at the good depth : 
      !!      Linear interpolation of T, S   
      !!         Computation of di(tb) and dj(tb) by vertical interpolation:
      !!          di(t) = t~ - t(i,j,k) or t(i+1,j,k) - t~
      !!          dj(t) = t~ - t(i,j,k) or t(i,j+1,k) - t~
      !!         This formulation computes the two cases:
      !!                 CASE 1                   CASE 2  
      !!         k-1  ___ ___________   k-1   ___ ___________
      !!                    Ti  T~                  T~  Ti+1
      !!                  _____                        _____
      !!         k        |   |Ti+1     k           Ti |   |
      !!                  |   |____                ____|   |
      !!              ___ |   |   |           ___  |   |   |
      !!                  
      !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
      !!          t~ = t(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(Ti+1)/e3w(i+1)
      !!        ( t~ = t(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(Tj+1)/e3w(j+1)  )
      !!          or
      !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
      !!          t~ = t(i,j,k) + (e3w(i) - e3w(i+1)) * dk(Ti)/e3w(i )
      !!        ( t~ = t(i,j,k) + (e3w(j) - e3w(j+1)) * dk(Tj)/e3w(j ) )
      !!          Idem for di(s) and dj(s)          
      !!
      !!      For rho, we call eos which will compute rd~(t~,s~) at the right
      !!      depth zh from interpolated T and S for the different formulations
      !!      of the equation of state (eos).
      !!      Gradient formulation for rho :
      !!          di(rho) = rd~ - rd(i,j,k)   or   rd(i+1,j,k) - rd~
      !!
      !! ** Action  : compute for top and bottom interfaces
      !!              - pgtu, pgtv, pgtui, pgtvi: horizontal gradient of tracer at u- & v-points
      !!              - pgru, pgrv, pgrui, pgtvi: horizontal gradient of rho (if present) at u- & v-points
      !!              - pmru, pmrv, pmrui, pmrvi: horizontal sum of rho at u- & v- point (used in dynhpg with vvl)
      !!              - pgzu, pgzv, pgzui, pgzvi: horizontal gradient of z at u- and v- point (used in dynhpg with vvl)
      !!              - pge3ru, pge3rv, pge3rui, pge3rvi: horizontal gradient of rho weighted by local e3w at u- & v-points 
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   )           ::  kt          ! ocean time-step index
      INTEGER                              , INTENT(in   )           ::  kjpt        ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   )           ::  pta         ! 4D tracers fields
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(  out)           ::  pgtu, pgtv  ! hor. grad. of ptra at u- & v-pts 
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(  out)           ::  pgtui, pgtvi  ! hor. grad. of stra at u- & v-pts (ISF)
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ), OPTIONAL ::  prd         ! 3D density anomaly fields
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pgru, pgrv      ! hor. grad of prd at u- & v-pts (bottom)
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pmru, pmrv      ! hor. sum  of prd at u- & v-pts (bottom)
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pgzu, pgzv      ! hor. grad of z   at u- & v-pts (bottom)
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pge3ru, pge3rv  ! hor. grad of prd weighted by local e3w at u- & v-pts (bottom)
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pgrui, pgrvi      ! hor. grad of prd at u- & v-pts (top)
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pmrui, pmrvi      ! hor. sum  of prd at u- & v-pts (top)
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pgzui, pgzvi      ! hor. grad of z   at u- & v-pts (top)
      REAL(wp), DIMENSION(jpi,jpj         ), INTENT(  out), OPTIONAL ::  pge3rui, pge3rvi  ! hor. grad of prd weighted by local e3w at u- & v-pts (top)
      !
      INTEGER  ::   ji, jj, jn      ! Dummy loop indices
      INTEGER  ::   iku, ikv, ikum1, ikvm1,ikup1, ikvp1   ! partial step level (ocean bottom level) at u- and v-points
      REAL(wp) ::  ze3wu, ze3wv, zmaxu, zmaxv, zdzwu, zdzwv, zdzwuip1, zdzwvjp1  ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj)      ::  zri, zrj, zhi, zhj   ! NB: 3rd dim=1 to use eos
      REAL(wp), DIMENSION(jpi,jpj,kjpt) ::  zti, ztj             ! 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'zps_hde_isf')
      !
      pgtu(:,:,:)=0.0_wp ; pgtv(:,:,:)=0.0_wp ;
      pgtui(:,:,:)=0.0_wp ; pgtvi(:,:,:)=0.0_wp ;
      zti (:,:,:)=0.0_wp ; ztj (:,:,:)=0.0_wp ;
      zhi (:,:  )=0.0_wp ; zhj (:,:  )=0.0_wp ;
      !
      DO jn = 1, kjpt      !==   Interpolation of tracers at the last ocean level   ==!
         !
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)   ;   ikum1 = MAX( iku - 1 , 1 )    ! last and before last ocean level at u- & v-points
               ikv = mbkv(ji,jj)   ;   ikvm1 = MAX( ikv - 1 , 1 )    ! if level first is a p-step, ik.m1=1
               ! (ISF) case partial step top and bottom in adjacent cell in vertical
               ! cannot used e3w because if 2 cell water column, we have ps at top and bottom
               ! in this case e3w(i,j) - e3w(i,j+1) is not the distance between Tj~ and Tj
               ! the only common depth between cells (i,j) and (i,j+1) is gdepw_0
               ze3wu  = (gdept_0(ji+1,jj,iku) - gdepw_0(ji+1,jj,iku)) - (gdept_0(ji,jj,iku) - gdepw_0(ji,jj,iku))
               ze3wv  = (gdept_0(ji,jj+1,ikv) - gdepw_0(ji,jj+1,ikv)) - (gdept_0(ji,jj,ikv) - gdepw_0(ji,jj,ikv))
               !
               ! i- direction
               IF( ze3wu >= 0._wp ) THEN      ! case 1
                  zmaxu =  ze3wu / fse3w(ji+1,jj,iku)
                  ! interpolated values of tracers
                  zti (ji,jj,jn) = pta(ji+1,jj,iku,jn) + zmaxu * ( pta(ji+1,jj,ikum1,jn) - pta(ji+1,jj,iku,jn) )
                  ! gradient of  tracers
                  pgtu(ji,jj,jn) = umask(ji,jj,iku) * ( zti(ji,jj,jn) - pta(ji,jj,iku,jn) )
               ELSE                           ! case 2
                  zmaxu = -ze3wu / fse3w(ji,jj,iku)
                  ! interpolated values of tracers
                  zti (ji,jj,jn) = pta(ji,jj,iku,jn) + zmaxu * ( pta(ji,jj,ikum1,jn) - pta(ji,jj,iku,jn) )
                  ! gradient of tracers
                  pgtu(ji,jj,jn) = umask(ji,jj,iku) * ( pta(ji+1,jj,iku,jn) - zti(ji,jj,jn) )
               ENDIF
               !
               ! j- direction
               IF( ze3wv >= 0._wp ) THEN      ! case 1
                  zmaxv =  ze3wv / fse3w(ji,jj+1,ikv)
                  ! interpolated values of tracers
                  ztj (ji,jj,jn) = pta(ji,jj+1,ikv,jn) + zmaxv * ( pta(ji,jj+1,ikvm1,jn) - pta(ji,jj+1,ikv,jn) )
                  ! gradient of tracers
                  pgtv(ji,jj,jn) = vmask(ji,jj,ikv) * ( ztj(ji,jj,jn) - pta(ji,jj,ikv,jn) )
               ELSE                           ! case 2
                  zmaxv =  -ze3wv / fse3w(ji,jj,ikv)
                  ! interpolated values of tracers
                  ztj (ji,jj,jn) = pta(ji,jj,ikv,jn) + zmaxv * ( pta(ji,jj,ikvm1,jn) - pta(ji,jj,ikv,jn) )
                  ! gradient of tracers
                  pgtv(ji,jj,jn) = vmask(ji,jj,ikv) * ( pta(ji,jj+1,ikv,jn) - ztj(ji,jj,jn) )
               ENDIF
            END DO
         END DO
         CALL lbc_lnk( pgtu(:,:,jn), 'U', -1. )   ;   CALL lbc_lnk( pgtv(:,:,jn), 'V', -1. )   ! Lateral boundary cond.
         !
      END DO

      ! horizontal derivative of density anomalies (rd)
      IF( PRESENT( prd ) ) THEN         ! depth of the partial step level
         pgru(:,:)=0.0_wp   ; pgrv(:,:)=0.0_wp ; 
         pgzu(:,:)=0.0_wp   ; pgzv(:,:)=0.0_wp ;
         pmru(:,:)=0.0_wp   ; pmru(:,:)=0.0_wp ;
         pge3ru(:,:)=0.0_wp ; pge3rv(:,:)=0.0_wp ;
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)
               ikv = mbkv(ji,jj)
               ze3wu  = (gdept_0(ji+1,jj,iku) - gdepw_0(ji+1,jj,iku)) - (gdept_0(ji,jj,iku) - gdepw_0(ji,jj,iku))
               ze3wv  = (gdept_0(ji,jj+1,ikv) - gdepw_0(ji,jj+1,ikv)) - (gdept_0(ji,jj,ikv) - gdepw_0(ji,jj,ikv))

               IF( ze3wu >= 0._wp ) THEN   ;   zhi(ji,jj) = fsdept(ji+1,jj,iku) - ze3wu     ! i-direction: case 1
               ELSE                        ;   zhi(ji,jj) = fsdept(ji  ,jj,iku) + ze3wu    ! -     -      case 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN   ;   zhj(ji,jj) = fsdept(ji,jj+1,ikv) - ze3wv    ! j-direction: case 1
               ELSE                        ;   zhj(ji,jj) = fsdept(ji,jj  ,ikv) + ze3wv    ! -     -      case 2
               ENDIF
            END DO
         END DO
         
         ! Compute interpolated rd from zti, ztj for the 2 cases at the depth of the partial
         ! step and store it in  zri, zrj for each  case
         CALL eos( zti, zhi, zri )  
         CALL eos( ztj, zhj, zrj )

         ! Gradient of density at the last level 
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj) ; ikum1 = MAX( iku - 1 , 1 )    ! last and before last ocean level at u- & v-points
               ikv = mbkv(ji,jj) ; ikvm1 = MAX( ikv - 1 , 1 )    ! last and before last ocean level at u- & v-points
               ze3wu  = (gdept_0(ji+1,jj,iku) - gdepw_0(ji+1,jj,iku)) - (gdept_0(ji,jj,iku) - gdepw_0(ji,jj,iku))
               ze3wv  = (gdept_0(ji,jj+1,ikv) - gdepw_0(ji,jj+1,ikv)) - (gdept_0(ji,jj,ikv) - gdepw_0(ji,jj,ikv))
               IF( ze3wu >= 0._wp ) THEN 
                  pgzu(ji,jj) = (fsde3w(ji+1,jj,iku) - ze3wu) - fsde3w(ji,jj,iku)
                  pgru(ji,jj) = umask(ji,jj,iku) * ( zri(ji  ,jj) - prd(ji,jj,iku) )   ! i: 1
                  pmru(ji,jj) = umask(ji,jj,iku) * ( zri(ji  ,jj) + prd(ji,jj,iku) )   ! i: 1 
                  pge3ru(ji,jj) = umask(ji,jj,iku)                                                                  &
                                * ( (fse3w(ji+1,jj,iku) - ze3wu )* ( zri(ji  ,jj    ) + prd(ji+1,jj,ikum1) + 2._wp) &
                                   - fse3w(ji  ,jj,iku)          * ( prd(ji  ,jj,iku) + prd(ji  ,jj,ikum1) + 2._wp) )  ! j: 2
               ELSE  
                  pgzu(ji,jj) = fsde3w(ji+1,jj,iku) - (fsde3w(ji,jj,iku) + ze3wu)
                  pgru(ji,jj) = umask(ji,jj,iku) * ( prd(ji+1,jj,iku) - zri(ji,jj) )   ! i: 2
                  pmru(ji,jj) = umask(ji,jj,iku) * ( prd(ji+1,jj,iku) + zri(ji,jj) )   ! i: 2
                  pge3ru(ji,jj) = umask(ji,jj,iku)                                                                  &
                                * (  fse3w(ji+1,jj,iku)          * ( prd(ji+1,jj,iku) + prd(ji+1,jj,ikum1) + 2._wp) &
                                   -(fse3w(ji  ,jj,iku) + ze3wu) * ( zri(ji  ,jj    ) + prd(ji  ,jj,ikum1) + 2._wp) )  ! j: 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN
                  pgzv(ji,jj) = (fsde3w(ji,jj+1,ikv) - ze3wv) - fsde3w(ji,jj,ikv) 
                  pgrv(ji,jj) = vmask(ji,jj,ikv) * ( zrj(ji,jj  ) - prd(ji,jj,ikv) )   ! j: 1
                  pmrv(ji,jj) = vmask(ji,jj,ikv) * ( zrj(ji,jj  ) + prd(ji,jj,ikv) )   ! j: 1
                  pge3rv(ji,jj) = vmask(ji,jj,ikv)                                                                  &
                                * ( (fse3w(ji,jj+1,ikv) - ze3wv )* ( zrj(ji,jj      ) + prd(ji,jj+1,ikvm1) + 2._wp) &
                                   - fse3w(ji,jj  ,ikv)          * ( prd(ji,jj  ,ikv) + prd(ji,jj  ,ikvm1) + 2._wp) )  ! j: 2
               ELSE 
                  pgzv(ji,jj) = fsde3w(ji,jj+1,ikv) - (fsde3w(ji,jj,ikv) + ze3wv)
                  pgrv(ji,jj) = vmask(ji,jj,ikv) * ( prd(ji,jj+1,ikv) - zrj(ji,jj) )   ! j: 2
                  pmrv(ji,jj) = vmask(ji,jj,ikv) * ( prd(ji,jj+1,ikv) + zrj(ji,jj) )   ! j: 2
                  pge3rv(ji,jj) = vmask(ji,jj,ikv)                                                                  &
                                * (  fse3w(ji,jj+1,ikv)          * ( prd(ji,jj+1,ikv) + prd(ji,jj+1,ikvm1) + 2._wp) &
                                   -(fse3w(ji,jj  ,ikv) + ze3wv) * ( zrj(ji,jj      ) + prd(ji,jj  ,ikvm1) + 2._wp) )  ! j: 2
               ENDIF
            END DO
         END DO
         CALL lbc_lnk( pgru   , 'U', -1. )   ;   CALL lbc_lnk( pgrv   , 'V', -1. )   ! Lateral boundary conditions
         CALL lbc_lnk( pmru   , 'U',  1. )   ;   CALL lbc_lnk( pmrv   , 'V',  1. )   ! Lateral boundary conditions
         CALL lbc_lnk( pgzu   , 'U', -1. )   ;   CALL lbc_lnk( pgzv   , 'V', -1. )   ! Lateral boundary conditions
         CALL lbc_lnk( pge3ru , 'U', -1. )   ;   CALL lbc_lnk( pge3rv , 'V', -1. )   ! Lateral boundary conditions
         !
      END IF
         ! (ISH)  compute grui and gruvi
      DO jn = 1, kjpt      !==   Interpolation of tracers at the last ocean level   ==!            !
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = miku(ji,jj)   ;  ikup1 = miku(ji,jj) + 1
               ikv = mikv(ji,jj)   ;  ikvp1 = mikv(ji,jj) + 1
               !
               ! (ISF) case partial step top and bottom in adjacent cell in vertical
               ! cannot used e3w because if 2 cell water column, we have ps at top and bottom
               ! in this case e3w(i,j) - e3w(i,j+1) is not the distance between Tj~ and Tj
               ! the only common depth between cells (i,j) and (i,j+1) is gdepw_0
               ze3wu  = (gdepw_0(ji+1,jj,iku+1) - gdept_0(ji+1,jj,iku)) - (gdepw_0(ji,jj,iku+1) - gdept_0(ji,jj,iku)) 
               ze3wv  = (gdepw_0(ji,jj+1,ikv+1) - gdept_0(ji,jj+1,ikv)) - (gdepw_0(ji,jj,ikv+1) - gdept_0(ji,jj,ikv))
               ! i- direction
               IF( ze3wu >= 0._wp ) THEN      ! case 1
                  zmaxu = ze3wu / fse3w(ji+1,jj,iku+1)
                  ! interpolated values of tracers
                  zti(ji,jj,jn) = pta(ji+1,jj,iku,jn) + zmaxu * ( pta(ji+1,jj,iku+1,jn) - pta(ji+1,jj,iku,jn) )
                  ! gradient of tracers
                  pgtui(ji,jj,jn) = umask(ji,jj,iku) * ( zti(ji,jj,jn) - pta(ji,jj,iku,jn) )
               ELSE                           ! case 2
                  zmaxu = - ze3wu / fse3w(ji,jj,iku+1)
                  ! interpolated values of tracers
                  zti(ji,jj,jn) = pta(ji,jj,iku,jn) + zmaxu * ( pta(ji,jj,iku+1,jn) - pta(ji,jj,iku,jn) )
                  ! gradient of  tracers
                  pgtui(ji,jj,jn) = umask(ji,jj,iku) * ( pta(ji+1,jj,iku,jn) - zti(ji,jj,jn) )
               ENDIF
               !
               ! j- direction
               IF( ze3wv >= 0._wp ) THEN      ! case 1
                  zmaxv =  ze3wv / fse3w(ji,jj+1,ikv+1)
                  ! interpolated values of tracers
                  ztj(ji,jj,jn) = pta(ji,jj+1,ikv,jn) + zmaxv * ( pta(ji,jj+1,ikv+1,jn) - pta(ji,jj+1,ikv,jn) )
                  ! gradient of tracers
                  pgtvi(ji,jj,jn) = vmask(ji,jj,ikv) * ( ztj(ji,jj,jn) - pta(ji,jj,ikv,jn) )
               ELSE                           ! case 2
                  zmaxv =  - ze3wv / fse3w(ji,jj,ikv+1)
                  ! interpolated values of tracers
                  ztj(ji,jj,jn) = pta(ji,jj,ikv,jn) + zmaxv * ( pta(ji,jj,ikv+1,jn) - pta(ji,jj,ikv,jn) )
                  ! gradient of tracers
                  pgtvi(ji,jj,jn) = vmask(ji,jj,ikv) * ( pta(ji,jj+1,ikv,jn) - ztj(ji,jj,jn) )
               ENDIF
            END DO!!
         END DO!!
         CALL lbc_lnk( pgtui(:,:,jn), 'U', -1. )   ;   CALL lbc_lnk( pgtvi(:,:,jn), 'V', -1. )   ! Lateral boundary cond.
         !
      END DO

      ! horizontal derivative of density anomalies (rd)
      IF( PRESENT( prd ) ) THEN         ! depth of the partial step level
         pgrui(:,:)  =0.0_wp ; pgrvi(:,:)  =0.0_wp ;
         pgzui(:,:)  =0.0_wp ; pgzvi(:,:)  =0.0_wp ;
         pmrui(:,:)  =0.0_wp ; pmrui(:,:)  =0.0_wp ;
         pge3rui(:,:)=0.0_wp ; pge3rvi(:,:)=0.0_wp ;

         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = miku(ji,jj)
               ikv = mikv(ji,jj)
               ze3wu  = (gdepw_0(ji+1,jj,iku+1) - gdept_0(ji+1,jj,iku)) - (gdepw_0(ji,jj,iku+1) - gdept_0(ji,jj,iku))
               ze3wv  = (gdepw_0(ji,jj+1,ikv+1) - gdept_0(ji,jj+1,ikv)) - (gdepw_0(ji,jj,ikv+1) - gdept_0(ji,jj,ikv))

               IF( ze3wu >= 0._wp ) THEN   ;   zhi(ji,jj) = fsdept(ji+1,jj,iku) + ze3wu    ! i-direction: case 1
               ELSE                        ;   zhi(ji,jj) = fsdept(ji  ,jj,iku) - ze3wu    ! -     -      case 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN   ;   zhj(ji,jj) = fsdept(ji,jj+1,ikv) + ze3wv    ! j-direction: case 1
               ELSE                        ;   zhj(ji,jj) = fsdept(ji,jj  ,ikv) - ze3wv    ! -     -      case 2
               ENDIF
            END DO
         END DO

         ! Compute interpolated rd from zti, ztj for the 2 cases at the depth of the partial
         ! step and store it in  zri, zrj for each  case
         CALL eos( zti, zhi, zri )  
         CALL eos( ztj, zhj, zrj )

         ! Gradient of density at the last level 
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = miku(ji,jj) ; ikup1 = miku(ji,jj) + 1
               ikv = mikv(ji,jj) ; ikvp1 = mikv(ji,jj) + 1
               ze3wu  = (gdepw_0(ji+1,jj,iku+1) - gdept_0(ji+1,jj,iku)) - (gdepw_0(ji,jj,iku+1) - gdept_0(ji,jj,iku))
               ze3wv  = (gdepw_0(ji,jj+1,ikv+1) - gdept_0(ji,jj+1,ikv)) - (gdepw_0(ji,jj,ikv+1) - gdept_0(ji,jj,ikv))
               IF( ze3wu >= 0._wp ) THEN
                 pgzui  (ji,jj) = (fsde3w(ji+1,jj,iku) + ze3wu) - fsde3w(ji,jj,iku)
                 pgrui  (ji,jj) = umask(ji,jj,iku)   * ( zri(ji,jj) - prd(ji,jj,iku) )          ! i: 1
                 pmrui  (ji,jj) = umask(ji,jj,iku)   * ( zri(ji,jj) + prd(ji,jj,iku) )          ! i: 1 
                 pge3rui(ji,jj) = umask(ji,jj,iku+1)                                                                  &
                                * ( (fse3w(ji+1,jj,iku+1) - ze3wu) * (zri(ji,jj    ) + prd(ji+1,jj,iku+1) + 2._wp)   &
                                   - fse3w(ji  ,jj,iku+1)          * (prd(ji,jj,iku) + prd(ji  ,jj,iku+1) + 2._wp)   ) ! i: 1
               ELSE
                 pgzui  (ji,jj) = fsde3w(ji+1,jj,iku) - (fsde3w(ji,jj,iku) - ze3wu)
                 pgrui  (ji,jj) = umask(ji,jj,iku)   * ( prd(ji+1,jj,iku) - zri(ji,jj) )      ! i: 2
                 pmrui  (ji,jj) = umask(ji,jj,iku)   * ( prd(ji+1,jj,iku) + zri(ji,jj) )      ! i: 2
                 pge3rui(ji,jj) = umask(ji,jj,iku+1)                                                                   &
                                * (  fse3w(ji+1,jj,iku+1)          * (prd(ji+1,jj,iku) + prd(ji+1,jj,iku+1) + 2._wp)  &
                                   -(fse3w(ji  ,jj,iku+1) + ze3wu) * (zri(ji,jj      ) + prd(ji  ,jj,iku+1) + 2._wp)  )     ! i: 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN
                 pgzvi  (ji,jj) = (fsde3w(ji,jj+1,ikv) + ze3wv) - fsde3w(ji,jj,ikv) 
                 pgrvi  (ji,jj) = vmask(ji,jj,ikv)   * ( zrj(ji,jj  ) - prd(ji,jj,ikv) )        ! j: 1
                 pmrvi  (ji,jj) = vmask(ji,jj,ikv)   * ( zrj(ji,jj  ) + prd(ji,jj,ikv) )        ! j: 1
                 pge3rvi(ji,jj) = vmask(ji,jj,ikv+1)                                                                  & 
                                * ( (fse3w(ji,jj+1,ikv+1) - ze3wv) * ( zrj(ji,jj    ) + prd(ji,jj+1,ikv+1) + 2._wp)  &
                                   - fse3w(ji,jj  ,ikv+1)          * ( prd(ji,jj,ikv) + prd(ji,jj  ,ikv+1) + 2._wp)  ) ! j: 1
                                  ! + 2 due to the formulation in density and not in anomalie in hpg sco
               ELSE
                 pgzvi  (ji,jj) = fsde3w(ji,jj+1,ikv) - (fsde3w(ji,jj,ikv) - ze3wv)
                 pgrvi  (ji,jj) = vmask(ji,jj,ikv)   * ( prd(ji,jj+1,ikv) - zrj(ji,jj) )     ! j: 2
                 pmrvi  (ji,jj) = vmask(ji,jj,ikv)   * ( prd(ji,jj+1,ikv) + zrj(ji,jj) )     ! j: 2
                 pge3rvi(ji,jj) = vmask(ji,jj,ikv+1)                                                                   &
                                * (  fse3w(ji,jj+1,ikv+1)          * ( prd(ji,jj+1,ikv) + prd(ji,jj+1,ikv+1) + 2._wp) &
                                   -(fse3w(ji,jj  ,ikv+1) + ze3wv) * ( zrj(ji,jj      ) + prd(ji,jj  ,ikv+1) + 2._wp) )  ! j: 2
               ENDIF
            END DO
         END DO
         CALL lbc_lnk( pgrui   , 'U', -1. )   ;   CALL lbc_lnk( pgrvi   , 'V', -1. )   ! Lateral boundary conditions
         CALL lbc_lnk( pmrui   , 'U',  1. )   ;   CALL lbc_lnk( pmrvi   , 'V',  1. )   ! Lateral boundary conditions
         CALL lbc_lnk( pgzui   , 'U', -1. )   ;   CALL lbc_lnk( pgzvi   , 'V', -1. )   ! Lateral boundary conditions
         CALL lbc_lnk( pge3rui , 'U', -1. )   ;   CALL lbc_lnk( pge3rvi , 'V', -1. )   ! Lateral boundary conditions
         !
      END IF  
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'zps_hde_isf')
      !
   END SUBROUTINE zps_hde_isf
   !!======================================================================
END MODULE zpshde
