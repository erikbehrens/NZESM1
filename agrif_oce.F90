MODULE agrif_oce
   !!======================================================================
   !!                       ***  MODULE agrif_oce  ***
   !! AGRIF :   define in memory AGRIF variables
   !!----------------------------------------------------------------------
   !! History :  2.0  ! 2007-12  (R. Benshila)  Original code
   !!----------------------------------------------------------------------
#if defined key_agrif
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!----------------------------------------------------------------------
   USE par_oce      ! ocean parameters
   USE dom_oce      ! domain parameters

   IMPLICIT NONE
   PRIVATE 

   PUBLIC agrif_oce_alloc ! routine called by nemo_init in nemogcm.F90

   !                                              !!* Namelist namagrif: AGRIF parameters
   LOGICAL , PUBLIC ::   ln_spc_dyn    = .FALSE.   !:
   INTEGER , PUBLIC ::   nn_cln_update = 3         !: update frequency 
   INTEGER , PUBLIC, PARAMETER ::   nn_sponge_len = 2  !: Sponge width (in number of parent grid points)
   REAL(wp), PUBLIC ::   rn_sponge_tra = 2800.     !: sponge coeff. for tracers
   REAL(wp), PUBLIC ::   rn_sponge_dyn = 2800.     !: sponge coeff. for dynamics
   LOGICAL , PUBLIC ::   ln_chk_bathy  = .FALSE.   !: check of parent bathymetry 

   !                                              !!! OLD namelist names
   INTEGER , PUBLIC ::   nbcline = 0               !: update counter
   INTEGER , PUBLIC ::   nbclineupdate             !: update frequency 
   REAL(wp), PUBLIC ::   visc_tra                  !: sponge coeff. for tracers
   REAL(wp), PUBLIC ::   visc_dyn                  !: sponge coeff. for dynamics

   LOGICAL , PUBLIC :: spongedoneT = .FALSE.       !: tracer   sponge layer indicator
   LOGICAL , PUBLIC :: spongedoneU = .FALSE.       !: dynamics sponge layer indicator
   LOGICAL , PUBLIC :: lk_agrif_fstep = .TRUE.     !: if true: first step
   LOGICAL , PUBLIC :: lk_agrif_doupd = .TRUE.     !: if true: send update from current grid
   LOGICAL , PUBLIC :: lk_agrif_debug = .FALSE.    !: if true: print debugging info

   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_tsn
# if defined key_top
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_trn
# endif
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_u
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_v
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,  DIMENSION(:,:) :: fsaht_spu, fsaht_spv !: sponge diffusivities
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,  DIMENSION(:,:) :: fsahm_spt, fsahm_spf !: sponge viscosities

   ! Barotropic arrays used to store open boundary data during
   ! time-splitting loop:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_w, vbdy_w, hbdy_w
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_e, vbdy_e, hbdy_e
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_n, vbdy_n, hbdy_n
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_s, vbdy_s, hbdy_s

   INTEGER :: tsn_id                                                  ! AGRIF profile for tracers interpolation and update
   INTEGER :: un_interp_id, vn_interp_id                              ! AGRIF profiles for interpolations
   INTEGER :: un_update_id, vn_update_id                              ! AGRIF profiles for udpates
   INTEGER :: tsn_sponge_id, un_sponge_id, vn_sponge_id               ! AGRIF profiles for sponge layers
# if defined key_top
   INTEGER :: trn_id, trn_sponge_id
# endif  
   INTEGER :: unb_id, vnb_id, ub2b_interp_id, vb2b_interp_id
   INTEGER :: ub2b_update_id, vb2b_update_id
   INTEGER :: e3t_id, e1u_id, e2v_id, sshn_id
   INTEGER :: scales_t_id
# if defined key_zdftke
   INTEGER :: avt_id, avm_id, en_id
# endif  
   INTEGER :: umsk_id, vmsk_id
   INTEGER :: kindic_agr

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3.1 , NEMO Consortium (2011)
   !! $Id: agrif_oce.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS 

   INTEGER FUNCTION agrif_oce_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION agrif_oce_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ierr
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( fsaht_spu(jpi,jpj), fsaht_spv(jpi,jpj),   &
         &      fsahm_spt(jpi,jpj), fsahm_spf(jpi,jpj),   &
         &      tabspongedone_tsn(jpi,jpj),           &
# if defined key_top         
         &      tabspongedone_trn(jpi,jpj),           &
# endif         
         &      tabspongedone_u  (jpi,jpj),           &
         &      tabspongedone_v  (jpi,jpj), STAT = ierr(1) )

      ALLOCATE( ubdy_w(jpj), vbdy_w(jpj), hbdy_w(jpj),   &
         &      ubdy_e(jpj), vbdy_e(jpj), hbdy_e(jpj),   & 
         &      ubdy_n(jpi), vbdy_n(jpi), hbdy_n(jpi),   & 
         &      ubdy_s(jpi), vbdy_s(jpi), hbdy_s(jpi), STAT = ierr(2) )

      agrif_oce_alloc = MAXVAL(ierr)
      !
   END FUNCTION agrif_oce_alloc

#endif
   !!======================================================================
END MODULE agrif_oce
