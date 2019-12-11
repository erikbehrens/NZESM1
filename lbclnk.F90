MODULE lbclnk
   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! Ocean        : lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)     Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)     F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment  
   !!            3.5  ! 2012     (S.Mocavero, I. Epicoco) Add 'lbc_bdy_lnk' 
   !!                            and lbc_obc_lnk' routine to optimize  
   !!                            the BDY/OBC communications
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  add a C1D case  
   !!            3.6  ! 2015-06  (O. Tintó and M. Castrillo)  add lbc_lnk_multi  
   !!----------------------------------------------------------------------
#if defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!   lbc_lnk      : generic interface for mpp_lnk_3d and mpp_lnk_2d routines defined in lib_mpp
   !!   lbc_lnk_e    : generic interface for mpp_lnk_2d_e routine defined in lib_mpp
   !!   lbc_bdy_lnk  : generic interface for mpp_lnk_bdy_2d and mpp_lnk_bdy_3d routines defined in lib_mpp
   !!----------------------------------------------------------------------
   USE lib_mpp          ! distributed memory computing library


   INTERFACE lbc_lnk_multi
      MODULE PROCEDURE mpp_lnk_2d_9, mpp_lnk_2d_multiple
   END INTERFACE

   INTERFACE lbc_lnk
      MODULE PROCEDURE mpp_lnk_3d_gather, mpp_lnk_3d, mpp_lnk_2d
   END INTERFACE

   INTERFACE lbc_bdy_lnk
      MODULE PROCEDURE mpp_lnk_bdy_2d, mpp_lnk_bdy_3d
   END INTERFACE

   INTERFACE lbc_lnk_e
      MODULE PROCEDURE mpp_lnk_2d_e
   END INTERFACE

   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE mpp_lnk_2d_icb
   END INTERFACE

   PUBLIC lbc_lnk       ! ocean lateral boundary conditions
   PUBLIC lbc_lnk_multi ! modified ocean lateral boundary conditions
   PUBLIC lbc_lnk_e
   PUBLIC lbc_bdy_lnk   ! ocean lateral BDY boundary conditions
   PUBLIC lbc_lnk_icb

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

#else
   !!----------------------------------------------------------------------
   !!   Default option                              shared memory computing
   !!----------------------------------------------------------------------
   !!   lbc_lnk      : generic interface for lbc_lnk_3d and lbc_lnk_2d
   !!   lbc_lnk_3d   : set the lateral boundary condition on a 3D variable on ocean mesh
   !!   lbc_lnk_2d   : set the lateral boundary condition on a 2D variable on ocean mesh
   !!   lbc_bdy_lnk  : set the lateral BDY boundary condition
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE in_out_manager  ! I/O manager
   USE lbcnfd          ! north fold

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_lnk
      MODULE PROCEDURE lbc_lnk_3d_gather, lbc_lnk_3d, lbc_lnk_2d
   END INTERFACE

   INTERFACE lbc_lnk_e
      MODULE PROCEDURE lbc_lnk_2d_e
   END INTERFACE

   INTERFACE lbc_lnk_multi
      MODULE PROCEDURE lbc_lnk_2d_9, lbc_lnk_2d_multiple
   END INTERFACE

   INTERFACE lbc_bdy_lnk
      MODULE PROCEDURE lbc_bdy_lnk_2d, lbc_bdy_lnk_3d
   END INTERFACE

   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE lbc_lnk_2d_e
   END INTERFACE
   
   TYPE arrayptr
      REAL , DIMENSION (:,:),  POINTER :: pt2d
   END TYPE arrayptr
   PUBLIC   arrayptr

   PUBLIC   lbc_lnk       ! ocean/ice  lateral boundary conditions
   PUBLIC   lbc_lnk_e 
   PUBLIC   lbc_lnk_multi ! modified ocean lateral boundary conditions
   PUBLIC   lbc_bdy_lnk   ! ocean lateral BDY boundary conditions
   PUBLIC   lbc_lnk_icb
   
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

# if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'                                          1D configuration
   !!----------------------------------------------------------------------

   SUBROUTINE lbc_lnk_3d_gather( pt3d1, cd_type1, pt3d2, cd_type2, psgn )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d_gather  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on two 3D arrays (C1D case)
      !!
      !! ** Method  :   call lbc_lnk_3d on pt3d1 and pt3d2
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   ) ::   cd_type1, cd_type2   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pt3d1   , pt3d2      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   ) ::   psgn                 ! control of the sign 
      !!----------------------------------------------------------------------
      !
      CALL lbc_lnk_3d( pt3d1, cd_type1, psgn)
      CALL lbc_lnk_3d( pt3d2, cd_type2, psgn)
      !
   END SUBROUTINE lbc_lnk_3d_gather


   SUBROUTINE lbc_lnk_3d( pt3d, cd_type, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 3D array (C1D case)
      !!
      !! ** Method  :   1D case, the central water column is set everywhere
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)           ::   pt3d      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   )           ::   psgn      ! control of the sign 
      CHARACTER(len=3)                , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      REAL(wp)                        , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !
      INTEGER  ::   jk     ! dummy loop index
      REAL(wp) ::   ztab   ! local scalar
      !!----------------------------------------------------------------------
      !
      DO jk = 1, jpk
         ztab = pt3d(2,2,jk)
         pt3d(:,:,jk) = ztab
      END DO
      !
   END SUBROUTINE lbc_lnk_3d


   SUBROUTINE lbc_lnk_2d( pt2d, cd_type, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 2D array (non mpp case)
      !!
      !! ** Method  :   1D case, the central water column is set everywhere
      !!----------------------------------------------------------------------
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign 
      CHARACTER(len=3)            , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      REAL(wp)                    , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !
      REAL(wp) ::   ztab   ! local scalar
      !!----------------------------------------------------------------------
      !
      ztab = pt2d(2,2)
      pt2d(:,:) = ztab
      !
   END SUBROUTINE lbc_lnk_2d
   
   SUBROUTINE lbc_lnk_2d_multiple( pt2d_array , type_array , psgn_array , num_fields )
      !!
      INTEGER :: num_fields
      TYPE( arrayptr ), DIMENSION(:) :: pt2d_array
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   type_array   ! define the nature of ptab array grid-points
      !                                                               ! = T , U , V , F , W and I points
      REAL(wp)        , DIMENSION(:), INTENT(in   ) ::   psgn_array   ! =-1 the sign change across the north fold boundary
      !                                                               ! =  1. , the sign is kept
      !
      INTEGER  ::   ii    !!MULTI SEND DUMMY LOOP INDICES
      !
      DO ii = 1, num_fields
        CALL lbc_lnk_2d( pt2d_array(ii)%pt2d, type_array(ii), psgn_array(ii) )
      END DO     
      !
   END SUBROUTINE lbc_lnk_2d_multiple

   SUBROUTINE lbc_lnk_2d_9( pt2dA, cd_typeA, psgnA, pt2dB, cd_typeB, psgnB, pt2dC, cd_typeC, psgnC   &
      &                   , pt2dD, cd_typeD, psgnD, pt2dE, cd_typeE, psgnE, pt2dF, cd_typeF, psgnF   &
      &                   , pt2dG, cd_typeG, psgnG, pt2dH, cd_typeH, psgnH, pt2dI, cd_typeI, psgnI, cd_mpp, pval)
      !!---------------------------------------------------------------------
      ! Second 2D array on which the boundary condition is applied
      REAL(wp), DIMENSION(jpi,jpj), TARGET          , INTENT(inout) ::   pt2dA
      REAL(wp), DIMENSION(jpi,jpj), TARGET, OPTIONAL, INTENT(inout) ::   pt2dB , pt2dC , pt2dD , pt2dE
      REAL(wp), DIMENSION(jpi,jpj), TARGET, OPTIONAL, INTENT(inout) ::   pt2dF , pt2dG , pt2dH , pt2dI
      ! define the nature of ptab array grid-points
      CHARACTER(len=1)                              , INTENT(in   ) ::   cd_typeA
      CHARACTER(len=1)                    , OPTIONAL, INTENT(in   ) ::   cd_typeB , cd_typeC , cd_typeD , cd_typeE
      CHARACTER(len=1)                    , OPTIONAL, INTENT(in   ) ::   cd_typeF , cd_typeG , cd_typeH , cd_typeI
      ! =-1 the sign change across the north fold boundary
      REAL(wp)                                      , INTENT(in   ) ::   psgnA
      REAL(wp)                            , OPTIONAL, INTENT(in   ) ::   psgnB , psgnC , psgnD , psgnE
      REAL(wp)                            , OPTIONAL, INTENT(in   ) ::   psgnF , psgnG , psgnH , psgnI
      CHARACTER(len=3)                    , OPTIONAL, INTENT(in   ) ::   cd_mpp   ! fill the overlap area only
      REAL(wp)                            , OPTIONAL, INTENT(in   ) ::   pval     ! background value (used at closed boundaries)
      !!
      !!---------------------------------------------------------------------

      !!The first array
      CALL lbc_lnk( pt2dA, cd_typeA, psgnA ) 

      !! Look if more arrays to process
      IF(PRESENT (psgnB) )CALL lbc_lnk( pt2dA, cd_typeA, psgnA ) 
      IF(PRESENT (psgnC) )CALL lbc_lnk( pt2dC, cd_typeC, psgnC ) 
      IF(PRESENT (psgnD) )CALL lbc_lnk( pt2dD, cd_typeD, psgnD ) 
      IF(PRESENT (psgnE) )CALL lbc_lnk( pt2dE, cd_typeE, psgnE ) 
      IF(PRESENT (psgnF) )CALL lbc_lnk( pt2dF, cd_typeF, psgnF ) 
      IF(PRESENT (psgnG) )CALL lbc_lnk( pt2dG, cd_typeG, psgnG ) 
      IF(PRESENT (psgnH) )CALL lbc_lnk( pt2dH, cd_typeH, psgnH ) 
      IF(PRESENT (psgnI) )CALL lbc_lnk( pt2dI, cd_typeI, psgnI ) 

   END SUBROUTINE lbc_lnk_2d_9





#else
   !!----------------------------------------------------------------------
   !!   Default option                           3D shared memory computing
   !!----------------------------------------------------------------------

   SUBROUTINE lbc_lnk_3d_gather( pt3d1, cd_type1, pt3d2, cd_type2, psgn )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d_gather  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on two 3D arrays (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   ) ::   cd_type1, cd_type2   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pt3d1   , pt3d2      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   ) ::   psgn                 ! control of the sign 
      !!----------------------------------------------------------------------
      !
      CALL lbc_lnk_3d( pt3d1, cd_type1, psgn)
      CALL lbc_lnk_3d( pt3d2, cd_type2, psgn)
      !
   END SUBROUTINE lbc_lnk_3d_gather


   SUBROUTINE lbc_lnk_3d( pt3d, cd_type, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 3D array (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)           ::   pt3d      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   )           ::   psgn      ! control of the sign 
      CHARACTER(len=3)                , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      REAL(wp)                        , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !!
      REAL(wp) ::   zland
      !!----------------------------------------------------------------------

      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value (zero by default)
      ELSE                         ;   zland = 0._wp
      ENDIF


      IF( PRESENT( cd_mpp ) ) THEN
         ! only fill the overlap area and extra allows 
         ! this is in mpp case. In this module, just do nothing
      ELSE
         !
         !                                     !  East-West boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )
         !
         CASE ( 1 , 4 , 6 )                       !**  cyclic east-west
            pt3d( 1 ,:,:) = pt3d(jpim1,:,:)            ! all points
            pt3d(jpi,:,:) = pt3d(  2  ,:,:)
            !
         CASE DEFAULT                             !**  East closed  --  West closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt3d( 1 ,:,:) = zland
               pt3d(jpi,:,:) = zland
            CASE ( 'F' )                               ! F-point
               pt3d(jpi,:,:) = zland
            END SELECT
            !
         END SELECT
         !
         !                                     ! North-South boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )
         !
         CASE ( 2 )                               !**  South symmetric  --  North closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
               pt3d(:, 1 ,:) = pt3d(:,3,:)
               pt3d(:,jpj,:) = zland
            CASE ( 'V' , 'F' )                         ! V-, F-points
               pt3d(:, 1 ,:) = psgn * pt3d(:,2,:)
               pt3d(:,jpj,:) = zland
            END SELECT
            !
         CASE ( 3 , 4 , 5 , 6 )                   !**  North fold  T or F-point pivot  --  South closed
            SELECT CASE ( cd_type )                    ! South : closed
            CASE ( 'T' , 'U' , 'V' , 'W' , 'I' )             ! all points except F-point
               pt3d(:, 1 ,:) = zland
            END SELECT
            !                                          ! North fold
            CALL lbc_nfd( pt3d(:,:,:), cd_type, psgn )
            !
         CASE DEFAULT                             !**  North closed  --  South closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt3d(:, 1 ,:) = zland
               pt3d(:,jpj,:) = zland
            CASE ( 'F' )                               ! F-point
               pt3d(:,jpj,:) = zland
            END SELECT
            !
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE lbc_lnk_3d

   SUBROUTINE lbc_lnk_2d( pt2d, cd_type, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 2D array (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign 
      CHARACTER(len=3)            , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      REAL(wp)                    , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !!
      REAL(wp) ::   zland
      !!----------------------------------------------------------------------

      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value (zero by default)
      ELSE                         ;   zland = 0._wp
      ENDIF

      IF (PRESENT(cd_mpp)) THEN
         ! only fill the overlap area and extra allows 
         ! this is in mpp case. In this module, just do nothing
      ELSE      
         !
         !                                     ! East-West boundaries
         !                                     ! ====================
         SELECT CASE ( nperio )
         !
         CASE ( 1 , 4 , 6 )                       !** cyclic east-west
            pt2d( 1 ,:) = pt2d(jpim1,:)               ! all points
            pt2d(jpi,:) = pt2d(  2  ,:)
            !
         CASE DEFAULT                             !** East closed  --  West closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )            ! T-, U-, V-, W-points
               pt2d( 1 ,:) = zland
               pt2d(jpi,:) = zland
            CASE ( 'F' )                              ! F-point
               pt2d(jpi,:) = zland
            END SELECT
            !
         END SELECT
         !
         !                                     ! North-South boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )
         !
         CASE ( 2 )                               !**  South symmetric  --  North closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
               pt2d(:, 1 ) = pt2d(:,3)
               pt2d(:,jpj) = zland
            CASE ( 'V' , 'F' )                         ! V-, F-points
               pt2d(:, 1 ) = psgn * pt2d(:,2)
               pt2d(:,jpj) = zland
            END SELECT
            !
         CASE ( 3 , 4 , 5 , 6 )                   !**  North fold  T or F-point pivot  --  South closed
            SELECT CASE ( cd_type )                    ! South : closed
            CASE ( 'T' , 'U' , 'V' , 'W' , 'I' )             ! all points except F-point
               pt2d(:, 1 ) = zland
            END SELECT
            !                                          ! North fold
            CALL lbc_nfd( pt2d(:,:), cd_type, psgn )
            !
         CASE DEFAULT                             !**  North closed  --  South closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt2d(:, 1 ) = zland
               pt2d(:,jpj) = zland
            CASE ( 'F' )                               ! F-point
               pt2d(:,jpj) = zland
            END SELECT
            !
         END SELECT
         !
      ENDIF
      !    
   END SUBROUTINE lbc_lnk_2d
   
   SUBROUTINE lbc_lnk_2d_multiple( pt2d_array , type_array , psgn_array , num_fields )
      !!
      INTEGER :: num_fields
      TYPE( arrayptr ), DIMENSION(:) :: pt2d_array
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   type_array   ! define the nature of ptab array grid-points
      !                                                               ! = T , U , V , F , W and I points
      REAL(wp)        , DIMENSION(:), INTENT(in   ) ::   psgn_array   ! =-1 the sign change across the north fold boundary
      !                                                               ! =  1. , the sign is kept
      !
      INTEGER  ::   ii    !!MULTI SEND DUMMY LOOP INDICES
      !
      DO ii = 1, num_fields
        CALL lbc_lnk_2d( pt2d_array(ii)%pt2d, type_array(ii), psgn_array(ii) )
      END DO     
      !
   END SUBROUTINE lbc_lnk_2d_multiple

   SUBROUTINE lbc_lnk_2d_9( pt2dA, cd_typeA, psgnA, pt2dB, cd_typeB, psgnB, pt2dC, cd_typeC, psgnC   &
      &                   , pt2dD, cd_typeD, psgnD, pt2dE, cd_typeE, psgnE, pt2dF, cd_typeF, psgnF   &
      &                   , pt2dG, cd_typeG, psgnG, pt2dH, cd_typeH, psgnH, pt2dI, cd_typeI, psgnI, cd_mpp, pval)
      !!---------------------------------------------------------------------
      ! Second 2D array on which the boundary condition is applied
      REAL(wp), DIMENSION(jpi,jpj), TARGET          , INTENT(inout) ::   pt2dA
      REAL(wp), DIMENSION(jpi,jpj), TARGET, OPTIONAL, INTENT(inout) ::   pt2dB , pt2dC , pt2dD , pt2dE
      REAL(wp), DIMENSION(jpi,jpj), TARGET, OPTIONAL, INTENT(inout) ::   pt2dF , pt2dG , pt2dH , pt2dI
      ! define the nature of ptab array grid-points
      CHARACTER(len=1)                              , INTENT(in   ) ::   cd_typeA
      CHARACTER(len=1)                    , OPTIONAL, INTENT(in   ) ::   cd_typeB , cd_typeC , cd_typeD , cd_typeE
      CHARACTER(len=1)                    , OPTIONAL, INTENT(in   ) ::   cd_typeF , cd_typeG , cd_typeH , cd_typeI
      ! =-1 the sign change across the north fold boundary
      REAL(wp)                                      , INTENT(in   ) ::   psgnA
      REAL(wp)                            , OPTIONAL, INTENT(in   ) ::   psgnB , psgnC , psgnD , psgnE
      REAL(wp)                            , OPTIONAL, INTENT(in   ) ::   psgnF , psgnG , psgnH , psgnI
      CHARACTER(len=3)                    , OPTIONAL, INTENT(in   ) ::   cd_mpp   ! fill the overlap area only
      REAL(wp)                            , OPTIONAL, INTENT(in   ) ::   pval     ! background value (used at closed boundaries)
      !!
      !!---------------------------------------------------------------------

      !!The first array
      CALL lbc_lnk( pt2dA, cd_typeA, psgnA ) 

      !! Look if more arrays to process
      IF(PRESENT (psgnB) )CALL lbc_lnk( pt2dA, cd_typeA, psgnA ) 
      IF(PRESENT (psgnC) )CALL lbc_lnk( pt2dC, cd_typeC, psgnC ) 
      IF(PRESENT (psgnD) )CALL lbc_lnk( pt2dD, cd_typeD, psgnD ) 
      IF(PRESENT (psgnE) )CALL lbc_lnk( pt2dE, cd_typeE, psgnE ) 
      IF(PRESENT (psgnF) )CALL lbc_lnk( pt2dF, cd_typeF, psgnF ) 
      IF(PRESENT (psgnG) )CALL lbc_lnk( pt2dG, cd_typeG, psgnG ) 
      IF(PRESENT (psgnH) )CALL lbc_lnk( pt2dH, cd_typeH, psgnH ) 
      IF(PRESENT (psgnI) )CALL lbc_lnk( pt2dI, cd_typeI, psgnI ) 

   END SUBROUTINE lbc_lnk_2d_9


#endif


   SUBROUTINE lbc_bdy_lnk_3d( pt3d, cd_type, psgn, ib_bdy )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_bdy_lnk  ***
      !!
      !! ** Purpose :   wrapper rountine to 'lbc_lnk_3d'. This wrapper is used
      !!                to maintain the same interface with regards to the mpp
      !case
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)           ::   pt3d      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   )           ::   psgn      ! control of the sign 
      INTEGER                                                   ::   ib_bdy    ! BDY boundary set
      !!
      CALL lbc_lnk_3d( pt3d, cd_type, psgn)

   END SUBROUTINE lbc_bdy_lnk_3d

   SUBROUTINE lbc_bdy_lnk_2d( pt2d, cd_type, psgn, ib_bdy )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_bdy_lnk  ***
      !!
      !! ** Purpose :   wrapper rountine to 'lbc_lnk_3d'. This wrapper is used
      !!                to maintain the same interface with regards to the mpp
      !case
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(inout)           ::   pt2d      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   )           ::   psgn      ! control of the sign 
      INTEGER                                                   ::   ib_bdy    ! BDY boundary set
      !!
      CALL lbc_lnk_2d( pt2d, cd_type, psgn)

   END SUBROUTINE lbc_bdy_lnk_2d


   SUBROUTINE lbc_lnk_2d_e( pt2d, cd_type, psgn, jpri, jprj )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 2D array (non mpp case)
      !!                special dummy routine to allow for use of halo indexing in mpp case
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign 
      INTEGER                     , INTENT(in   )           ::   jpri      ! size of extra halo (not needed in non-mpp)
      INTEGER                     , INTENT(in   )           ::   jprj      ! size of extra halo (not needed in non-mpp)
      !!----------------------------------------------------------------------

      CALL lbc_lnk_2d( pt2d, cd_type, psgn )
      !    
   END SUBROUTINE lbc_lnk_2d_e

#endif

   !!======================================================================
END MODULE lbclnk

