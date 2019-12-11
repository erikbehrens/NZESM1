MODULE dynspg_oce
   !!======================================================================
   !!                       ***  MODULE dynspg_oce  ***
   !!       
   !! Ocean dynamics: Define in memory some surface pressure gradient variables
   !!======================================================================
   !! History :  1.0  ! 2005-12  (C. Talandier, G. Madec)  Original code
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PUBLIC           

   PUBLIC   dynspg_oce_alloc   ! called in dynspg.F90
   
   !                                                       !!! Surface pressure gradient logicals
#if   defined key_dynspg_exp  ||  defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .TRUE.   !: Explicit free surface flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .FALSE.  !: Explicit free surface flag
#endif
#if   defined key_dynspg_ts   ||  defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .TRUE.   !: Free surface with time splitting flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .FALSE.  !: Free surface with time splitting flag
#endif
#if   defined key_dynspg_flt  ||  defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .TRUE.   !: Filtered free surface cst volume flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .FALSE.  !: Filtered free surface cst volume flag
#endif

  !                                                                         !!! Time splitting scheme (key_dynspg_ts) 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshn_e, ssha_e   ! sea surface heigth (now, after)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ua_e  , va_e     ! barotropic velocities (after)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hu_e  , hv_e     ! now ocean depth ( = Ho+sshn_e )
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hur_e , hvr_e    ! inverse of hu_e and hv_e
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   un_adv, vn_adv   ! Advection vel. at "now" barocl. step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ub2_b,  vb2_b    ! Half step fluxes (ln_bt_fw=T)
#if defined key_agrif
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ub2_i_b,  vb2_i_b! Half step time integrated fluxes 
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dynspg_oce_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dynspg_oce_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( sshn_e(jpi,jpj) , ua_e(jpi,jpj) , hu_e(jpi,jpj) , hur_e(jpi,jpj) ,      &
         &      ssha_e(jpi,jpj) , va_e(jpi,jpj) , hv_e(jpi,jpj) , hvr_e(jpi,jpj) ,      &
         &      ub2_b(jpi,jpj)  , vb2_b(jpi,jpj)                                 ,      &
#if defined key_agrif
         &      ub2_i_b(jpi,jpj), vb2_i_b(jpi,jpj)                               ,      &
#endif
         &      un_adv(jpi,jpj) , vn_adv(jpi,jpj)                                , STAT = dynspg_oce_alloc )
         !
      IF( lk_mpp                )   CALL mpp_sum ( dynspg_oce_alloc )
      IF( dynspg_oce_alloc /= 0 )   CALL ctl_warn('dynspg_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION dynspg_oce_alloc

   !!======================================================================
END MODULE dynspg_oce
