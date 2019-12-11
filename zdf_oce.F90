MODULE zdf_oce
   !!======================================================================
   !!              ***  MODULE  zdf_oce  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================
   !! history :  1.0  !  2002-06  (G. Madec) Original code
   !!            3.2  !  2009-07  (G.Madec) addition of avm
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC  zdf_oce_alloc    ! Called in nemogcm.F90

#if defined key_zdfcst   ||   defined key_esopa
   LOGICAL, PARAMETER, PUBLIC ::   lk_zdfcst        = .TRUE.         !: constant vertical mixing flag
#else
   LOGICAL, PARAMETER, PUBLIC ::   lk_zdfcst        = .FALSE.        !: constant vertical mixing flag
#endif

   !                                 !!* namelist namzdf: vertical diffusion *
   REAL(wp), PUBLIC ::   rn_avm0     !: vertical eddy viscosity (m2/s)
   REAL(wp), PUBLIC ::   rn_avt0     !: vertical eddy diffusivity (m2/s)
   INTEGER , PUBLIC ::   nn_avb      !: constant or profile background on avt (=0/1)
   INTEGER , PUBLIC ::   nn_havtb    !: horizontal shape or not for avtb (=0/1)
   LOGICAL , PUBLIC ::   ln_zdfexp   !: explicit vertical diffusion scheme flag
   INTEGER , PUBLIC ::   nn_zdfexp   !: number of sub-time step (explicit time stepping)
   LOGICAL , PUBLIC ::   ln_zdfevd   !: convection: enhanced vertical diffusion flag
   INTEGER , PUBLIC ::   nn_evdm     !: =0/1 flag to apply enhanced avm or not
   REAL(wp), PUBLIC ::   rn_avevd    !: vertical eddy coeff. for enhanced vert. diff. (m2/s)
   LOGICAL , PUBLIC ::   ln_zdfnpc   !: convection: non-penetrative convection flag
   INTEGER , PUBLIC ::   nn_npc      !: non penetrative convective scheme call  frequency
   INTEGER , PUBLIC ::   nn_npcp     !: non penetrative convective scheme print frequency


   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:)     ::   avmb , avtb    !: background profile of avm and avt
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   avtb_2d        !: horizontal shape of background Kz profile
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   bfrua, bfrva   !: Bottom friction coefficients set in zdfbfr
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   tfrua, tfrva   !: top friction coefficients set in zdfbfr
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avmu , avmv    !: vertical viscosity coef at uw- & vw-pts       [m2/s]
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avm  , avt     !: vertical viscosity & diffusivity coef at w-pt [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   avt_k , avm_k  ! not enhanced Kz
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   avmu_k, avmv_k ! not enhanced Kz
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   en              !: now turbulent kinetic energy   [m2/s2]
 
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_oce_alloc()
      !!----------------------------------------------------------------------
      !!            *** FUNCTION zdf_oce_alloc ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE(avmb(jpk) , bfrua(jpi,jpj) ,                         &
         &     avtb(jpk) , bfrva(jpi,jpj) , avtb_2d(jpi,jpj) ,      &
         &     tfrua(jpi, jpj), tfrva(jpi, jpj)              ,      &
         &     avmu(jpi,jpj,jpk), avm(jpi,jpj,jpk)           ,      &
         &     avmv  (jpi,jpj,jpk), avt   (jpi,jpj,jpk)      ,      &
         &     avt_k (jpi,jpj,jpk), avm_k (jpi,jpj,jpk)      ,      & 
         &     avmu_k(jpi,jpj,jpk), avmv_k(jpi,jpj,jpk)      ,      &
         &     en    (jpi,jpj,jpk), STAT = zdf_oce_alloc )
         !
      IF( zdf_oce_alloc /= 0 )   CALL ctl_warn('zdf_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION zdf_oce_alloc

   !!======================================================================
END MODULE zdf_oce
