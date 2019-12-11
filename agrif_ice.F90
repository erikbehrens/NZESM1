MODULE agrif_ice
   !!======================================================================
   !!                       ***  MODULE agrif_ice  ***
   !! AGRIF :   define in memory AGRIF variables for sea-ice
   !!----------------------------------------------------------------------
   !! History :  3.4  ! 2012-08  (R. Benshila)  Original code
   !!----------------------------------------------------------------------
#if defined key_agrif && defined key_lim2
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!----------------------------------------------------------------------
   USE par_oce      ! ocean parameters
   
   IMPLICIT NONE
   PRIVATE 

   PUBLIC agrif_ice_alloc ! routine called by nemo_init in nemogcm.F90

   INTEGER, PUBLIC :: u_ice_id, v_ice_id, adv_ice_id
   REAL(wp), PUBLIC :: lim_nbstep = 0.    ! child time position in sea-ice model
#if defined key_lim2_vp
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)     :: u_ice_nst, v_ice_nst   
#else
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   :: u_ice_oe, u_ice_sn     !: boundaries arrays
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   :: v_ice_oe, v_ice_sn     !:  "          " 
#endif
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) :: adv_ice_oe, adv_ice_sn !:  "          "

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3.4 , NEMO Consortium (2012)
   !! $Id: agrif_ice.F90 3680 2012-11-27 14:42:24Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS 

   INTEGER FUNCTION agrif_ice_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION agrif_ice_alloc  ***
      !!----------------------------------------------------------------------
#if defined key_lim2_vp
      ALLOCATE( u_ice_nst(jpi,jpj), v_ice_nst(jpi,jpj) ,   &
#else
      ALLOCATE( u_ice_oe(4,jpj,2) , v_ice_oe(4,jpj,2) ,    &
         &      u_ice_sn(jpi,4,2) , v_ice_sn(jpi,4,2) ,    &
#endif
         &      adv_ice_oe (4,jpj,7,2) , adv_ice_sn (jpi,4,7,2) ,   &
         &      STAT = agrif_ice_alloc)

#if ! defined key_lim2_vp
      u_ice_oe(:,:,:) =  0.e0
      v_ice_oe(:,:,:) =  0.e0
      u_ice_sn(:,:,:) =  0.e0
      v_ice_sn(:,:,:) =  0.e0
#endif
      adv_ice_oe (:,:,:,:) = 0.e0 
      adv_ice_sn (:,:,:,:) = 0.e0 
      !
   END FUNCTION agrif_ice_alloc

#endif
   !!======================================================================
END MODULE agrif_ice
