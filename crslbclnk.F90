MODULE crslbclnk

   !!======================================================================
   !!                       ***  MODULE  crslbclnk  ***
   !!               A temporary solution for lbclnk for coarsened grid.
   !! Ocean        : lateral boundary conditions for grid coarsening
   !!=====================================================================
   !! History :   ! 2012-06  (J. Simeon, G. Madec, C. Ethe, C. Calone)     Original code

   USE dom_oce
   USE crs
   USE lbclnk
   USE par_kind, ONLY: wp
   USE in_out_manager

   
   
   INTERFACE crs_lbc_lnk
      MODULE PROCEDURE crs_lbc_lnk_3d, crs_lbc_lnk_3d_gather, crs_lbc_lnk_2d
   END INTERFACE
   
   PUBLIC crs_lbc_lnk
   
   !! $Id$
CONTAINS

   SUBROUTINE crs_lbc_lnk_3d( pt3d1, cd_type1, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                  ***  SUBROUTINE crs_lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions for coarsened grid
      !!
      !! ** Method  :   Swap domain indices from full to coarse domain
      !!                before arguments are passed directly to lbc_lnk.
      !!                Upon exiting, switch back to full domain indices.
      !!----------------------------------------------------------------------
      !! Arguments
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type1 ! grid type
      REAL(wp)                        , INTENT(in   )           ::   psgn     ! control of the sign

      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout)   ::   pt3d1 ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   ), OPTIONAL ::   pval     ! valeur sur les halo
      CHARACTER(len=3)                , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      
      !! local vairables
      LOGICAL                                                   ::   ll_grid_crs
      REAL(wp)                                                  ::   zval     ! valeur sur les halo

      !!----------------------------------------------------------------------
      
      ll_grid_crs = ( jpi == jpi_crs )
      
      IF( PRESENT(pval) ) THEN  ;  zval = pval
      ELSE                      ;  zval = 0.0
      ENDIF
      
      IF( .NOT. ll_grid_crs )   CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain

      IF( PRESENT( cd_mpp ) ) THEN ; CALL lbc_lnk( pt3d1, cd_type1, psgn, cd_mpp, pval=zval  )
      ELSE                         ; CALL lbc_lnk( pt3d1, cd_type1, psgn, pval=zval  )
      ENDIF

      IF( .NOT. ll_grid_crs )   CALL dom_grid_glo   ! Return to parent grid domain

   END SUBROUTINE crs_lbc_lnk_3d
   
   SUBROUTINE crs_lbc_lnk_3d_gather( pt3d1, cd_type1, pt3d2, cd_type2, psgn )
      !!---------------------------------------------------------------------
      !!                  ***  SUBROUTINE crs_lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions for coarsened grid
      !!
      !! ** Method  :   Swap domain indices from full to coarse domain
      !!                before arguments are passed directly to lbc_lnk.
      !!                Upon exiting, switch back to full domain indices.
      !!----------------------------------------------------------------------
      !! Arguments
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type1,cd_type2 ! grid type
      REAL(wp)                        , INTENT(in   )           ::   psgn     ! control of the sign

      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout)   ::   pt3d1,pt3d2 ! 3D array on which the lbc is applied
      
      !! local vairables
      LOGICAL                                                   ::   ll_grid_crs
      !!----------------------------------------------------------------------
      
      ll_grid_crs = ( jpi == jpi_crs )
      
      IF( .NOT. ll_grid_crs )   CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain

      CALL lbc_lnk( pt3d1, cd_type1, pt3d2, cd_type2, psgn  )

      IF( .NOT. ll_grid_crs )   CALL dom_grid_glo   ! Return to parent grid domain

   END SUBROUTINE crs_lbc_lnk_3d_gather

   
   
   SUBROUTINE crs_lbc_lnk_2d(pt2d, cd_type, psgn, cd_mpp, pval)
      !!---------------------------------------------------------------------
      !!                  ***  SUBROUTINE crs_lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions for coarsened grid
      !!
      !! ** Method  :   Swap domain indices from full to coarse domain
      !!                before arguments are passed directly to lbc_lnk.
      !!                Upon exiting, switch back to full domain indices.
      !!----------------------------------------------------------------------
      !! Arguments
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type  ! grid type
      REAL(wp)                        , INTENT(in   )           ::   psgn     ! control of the sign

      REAL(wp), DIMENSION(jpi_crs,jpj_crs),     INTENT(inout)   ::   pt2d     ! 2D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   ), OPTIONAL ::   pval     ! valeur sur les halo
      CHARACTER(len=3)                , INTENT(in   ), OPTIONAL ::   cd_mpp   ! MPP only (here do nothing)
      !! local variables
      
      LOGICAL                                                   ::   ll_grid_crs
      REAL(wp)                                                  ::   zval     ! valeur sur les halo

      !!----------------------------------------------------------------------
      
      ll_grid_crs = ( jpi == jpi_crs )
      
      IF( PRESENT(pval) ) THEN  ;  zval = pval
      ELSE                      ;  zval = 0.0
      ENDIF
      
      IF( .NOT. ll_grid_crs )   CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain

      IF( PRESENT( cd_mpp ) ) THEN ; CALL lbc_lnk( pt2d, cd_type, psgn, cd_mpp, pval=zval  )
      ELSE                         ; CALL lbc_lnk( pt2d, cd_type, psgn, pval=zval  )
      ENDIF

      IF( .NOT. ll_grid_crs )   CALL dom_grid_glo   ! Return to parent grid domain

   END SUBROUTINE crs_lbc_lnk_2d


END MODULE crslbclnk
