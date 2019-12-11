#define TWO_WAY

MODULE agrif_lim2_update
   !!======================================================================
   !!                       ***  MODULE agrif_lim2_update ***
   !! Nesting module :  update surface ocean boundary condition over ice
   !!                   from a child grif
   !! Sea-Ice model  :  LIM 2.0 Sea ice model time-stepping
   !!======================================================================
   !! History :  2.0   !  04-2008  (F. Dupont)  initial version
   !!            3.4   !  08-2012  (R. Benshila, C. Herbaut) update and EVP
   !!----------------------------------------------------------------------
#if defined key_agrif && defined key_lim2
   !!----------------------------------------------------------------------
   !!   'key_lim2'  :                                 LIM 2.0 sea-ice model
   !!   'key_agrif' :                                 AGRIF library 
   !!----------------------------------------------------------------------
   !!   agrif_update_lim2  : update sea-ice model on boundaries or total
   !!                        sea-ice area for velocities and ice properties
   !!   update_adv_ice     : sea-ice properties
   !!   update_u_ice       : zonal ice velocity
   !!   update_v_ice       : meridional ice velocity
   !!----------------------------------------------------------------------
   USE ice_2
   USE dom_ice_2
   USE sbc_oce
   USE dom_oce
   USE agrif_oce
   USE agrif_ice 

   IMPLICIT NONE
   PRIVATE

   PUBLIC agrif_update_lim2

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.4 , LOCEAN-IPSL (2012)
   !! $Id: agrif_lim2_update.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE agrif_update_lim2 ( kt )
      !!----------------------------------------------------------------------
      !!                     *** ROUTINE agrif_update_lim2 ***
      !! ** Method  :   Call the hydrostaticupdate pressure at the boundary or
      !!                the entire domain 
      !!
      !! ** Action : - Update (u_ice,v_ice) and ice tracers
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt
      !!
      !!----------------------------------------------------------------------
      !
      IF((Agrif_NbStepint() .NE. (Agrif_irhot()-1)).AND.(kt /= 0)) RETURN

      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid = 0.
# if defined TWO_WAY
      IF( MOD(nbcline,nbclineupdate) == 0) THEN
         CALL Agrif_Update_Variable( adv_ice_id , procname = update_adv_ice  )
         CALL Agrif_Update_Variable( u_ice_id   , procname = update_u_ice    )
         CALL Agrif_Update_Variable( v_ice_id   , procname = update_v_ice    )
      ELSE
         CALL Agrif_Update_Variable( adv_ice_id , locupdate=(/0,2/), procname = update_adv_ice  )
         CALL Agrif_Update_Variable( u_ice_id   , locupdate=(/0,1/), procname = update_u_ice    )
         CALL Agrif_Update_Variable( v_ice_id   , locupdate=(/0,1/), procname = update_v_ice    )
      ENDIF
# endif
      !
   END SUBROUTINE agrif_update_lim2


   SUBROUTINE update_adv_ice( tabres, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                        *** ROUTINE update_adv_ice ***
      !! ** Method  : Compute the mass properties on the fine grid and recover
      !!              the properties per mass on the coarse grid
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2,7), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      REAL(wp) :: zrhox, zrhoy
      REAL(wp) :: z1_area
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         zrhox = Agrif_Rhox()
         zrhoy = Agrif_Rhoy()
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj, 1) = frld  (ji,jj  ) * area(ji,jj)
               tabres(ji,jj, 2) = hicif (ji,jj  ) * area(ji,jj)
               tabres(ji,jj, 3) = hsnif (ji,jj  ) * area(ji,jj)
               tabres(ji,jj, 4) = tbif  (ji,jj,1) * area(ji,jj)
               tabres(ji,jj, 5) = tbif  (ji,jj,2) * area(ji,jj)
               tabres(ji,jj, 6) = tbif  (ji,jj,3) * area(ji,jj)
               tabres(ji,jj, 7) = qstoif(ji,jj  ) * area(ji,jj)
            END DO
         END DO
         tabres = zrhox * zrhoy * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               z1_area = 1. / area(ji,jj) * tms(ji,jj)
               frld  (ji,jj)   = tabres(ji,jj, 1) * z1_area
               hicif (ji,jj)   = tabres(ji,jj, 2) * z1_area
               hsnif (ji,jj)   = tabres(ji,jj, 3) * z1_area
               tbif  (ji,jj,1) = tabres(ji,jj, 4) * z1_area
               tbif  (ji,jj,2) = tabres(ji,jj, 5) * z1_area
               tbif  (ji,jj,3) = tabres(ji,jj, 6) * z1_area
               qstoif(ji,jj)   = tabres(ji,jj, 7) * z1_area
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE update_adv_ice


# if defined key_lim2_vp
   SUBROUTINE update_u_ice( tabres, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                        *** ROUTINE update_u_ice ***
      !! ** Method  : Update the fluxes and recover the properties (B-grid)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      REAL(wp) :: zrhoy
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         zrhoy = Agrif_Rhoy()
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               tabres(ji,jj) = e2f(ji-1,jj-1) * u_ice(ji,jj)
            END DO
         END DO
         tabres = zrhoy * tabres
      ELSE
         DO jj= MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               u_ice(ji,jj) = tabres(ji,jj) / (e2f(ji-1,jj-1))
               u_ice(ji,jj) = u_ice(ji,jj) * tmu(ji,jj)
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE update_u_ice


   SUBROUTINE update_v_ice( tabres, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE update_v_ice ***
      !! ** Method  : Update the fluxes and recover the properties (B-grid)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2),  INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      REAL(wp) :: zrhox
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         zrhox = Agrif_Rhox()
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               tabres(ji,jj) = e1f(ji-1,jj-1) * v_ice(ji,jj)
            END DO
         END DO
         tabres = zrhox * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               v_ice(ji,jj) = tabres(ji,jj) / (e1f(ji-1,jj-1))
               v_ice(ji,jj) = v_ice(ji,jj) * tmu(ji,jj)
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE update_v_ice
# else
   SUBROUTINE update_u_ice( tabres, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                        *** ROUTINE update_u_ice ***
      !! ** Method  : Update the fluxes and recover the properties (C-grid)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      REAL(wp) :: zrhoy
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         zrhoy = Agrif_Rhoy()
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               tabres(ji,jj) = e2u(ji,jj) * u_ice(ji,jj)
            END DO
         END DO
         tabres = zrhoy * tabres
      ELSE
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               u_ice(ji,jj) = tabres(ji,jj) / (e2u(ji,jj))
               u_ice(ji,jj) = u_ice(ji,jj) * tmu(ji,jj)
            END DO
         END DO
      ENDIF
      ! 
   END SUBROUTINE update_u_ice


   SUBROUTINE update_v_ice( tabres, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE update_v_ice ***
      !! ** Method  : Update the fluxes and recover the properties (C-grid)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2),  INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      REAL(wp) :: zrhox
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         zrhox = Agrif_Rhox()
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               tabres(ji,jj) = e1v(ji,jj) * v_ice(ji,jj)
            END DO
         END DO
         tabres = zrhox * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               v_ice(ji,jj) = tabres(ji,jj) / (e1v(ji,jj))
               v_ice(ji,jj) = v_ice(ji,jj) * tmv(ji,jj)
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE update_v_ice
# endif

#else
CONTAINS
   SUBROUTINE agrif_lim2_update_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_lim2_update_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_lim2_update : You should not have seen this print! error?'
   END SUBROUTINE agrif_lim2_update_empty
#endif
END MODULE agrif_lim2_update
