#define TWO_WAY        /* TWO WAY NESTING */
#undef DECAL_FEEDBACK /* SEPARATION of INTERFACES*/
 
MODULE agrif_opa_update
#if defined key_agrif  && ! defined key_offline
   USE par_oce
   USE oce
   USE dom_oce
   USE agrif_oce
   USE in_out_manager  ! I/O manager
   USE lib_mpp
   USE wrk_nemo  
   USE dynspg_oce
   USE zdf_oce        ! vertical physics: ocean variables 

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Update_Tra, Agrif_Update_Dyn,Update_Scales
# if defined key_zdftke
   PUBLIC Agrif_Update_Tke
# endif
   !!----------------------------------------------------------------------
   !! NEMO/NST 3.6 , NEMO Consortium (2010)
   !! $Id: agrif_opa_update.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   RECURSIVE SUBROUTINE Agrif_Update_Tra( )
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Update_Tra ***
      !!---------------------------------------------
      ! 
      IF (Agrif_Root()) RETURN
      !
#if defined TWO_WAY  
      IF (lwp.AND.lk_agrif_debug) Write(*,*) 'Update tracers  from grid Number',Agrif_Fixed(), 'nbcline', nbcline

      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid = 0.
      ! 
      IF (MOD(nbcline,nbclineupdate) == 0) THEN
# if ! defined DECAL_FEEDBACK
         CALL Agrif_Update_Variable(tsn_id, procname=updateTS)
# else
         CALL Agrif_Update_Variable(tsn_id, locupdate=(/1,0/),procname=updateTS)
# endif
      ELSE
# if ! defined DECAL_FEEDBACK
         CALL Agrif_Update_Variable(tsn_id,locupdate=(/0,2/), procname=updateTS)
# else
         CALL Agrif_Update_Variable(tsn_id,locupdate=(/1,2/), procname=updateTS)
# endif
      ENDIF
      !
      Agrif_UseSpecialValueInUpdate = .FALSE.
      !
      IF ( lk_agrif_doupd ) THEN ! Initialisation: do recursive update:
         CALL Agrif_ChildGrid_To_ParentGrid()
         CALL Agrif_Update_Tra()
         CALL Agrif_ParentGrid_To_ChildGrid()
      ENDIF
      !
#endif
      !
   END SUBROUTINE Agrif_Update_Tra

   RECURSIVE SUBROUTINE Agrif_Update_Dyn( )
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Update_Dyn ***
      !!---------------------------------------------
      ! 
      IF (Agrif_Root()) RETURN
      !
#if defined TWO_WAY
      IF (lwp.AND.lk_agrif_debug) Write(*,*) 'Update momentum from grid Number',Agrif_Fixed(), 'nbcline', nbcline

      Agrif_UseSpecialValueInUpdate = .FALSE.
      Agrif_SpecialValueFineGrid = 0.
      !     
      IF (mod(nbcline,nbclineupdate) == 0) THEN
# if ! defined DECAL_FEEDBACK
         CALL Agrif_Update_Variable(un_update_id,procname = updateU)
         CALL Agrif_Update_Variable(vn_update_id,procname = updateV)
# else
         CALL Agrif_Update_Variable(un_update_id,locupdate1=(/0,-1/),locupdate2=(/1,-2/),procname = updateU)
         CALL Agrif_Update_Variable(vn_update_id,locupdate1=(/1,-2/),locupdate2=(/0,-1/),procname = updateV)
# endif
      ELSE
# if ! defined DECAL_FEEDBACK
         CALL Agrif_Update_Variable(un_update_id,locupdate=(/0,1/),procname = updateU)
         CALL Agrif_Update_Variable(vn_update_id,locupdate=(/0,1/),procname = updateV)         
# else
         CALL Agrif_Update_Variable(un_update_id,locupdate1=(/0,1/),locupdate2=(/1,1/),procname = updateU)
         CALL Agrif_Update_Variable(vn_update_id,locupdate1=(/1,1/),locupdate2=(/0,1/),procname = updateV)
# endif
      ENDIF

# if ! defined DECAL_FEEDBACK
      CALL Agrif_Update_Variable(e1u_id,procname = updateU2d)
      CALL Agrif_Update_Variable(e2v_id,procname = updateV2d)  
# else
      CALL Agrif_Update_Variable(e1u_id,locupdate1=(/0,-1/),locupdate2=(/1,-2/),procname = updateU2d)
      CALL Agrif_Update_Variable(e2v_id,locupdate1=(/1,-2/),locupdate2=(/0,-1/),procname = updateV2d)  
# endif

# if defined key_dynspg_ts
      IF (ln_bt_fw) THEN
         ! Update time integrated transports
         IF (mod(nbcline,nbclineupdate) == 0) THEN
#  if ! defined DECAL_FEEDBACK
            CALL Agrif_Update_Variable(ub2b_update_id,procname = updateub2b)
            CALL Agrif_Update_Variable(vb2b_update_id,procname = updatevb2b)
#  else
            CALL Agrif_Update_Variable(ub2b_update_id,locupdate1=(/0,-1/),locupdate2=(/1,-2/),procname = updateub2b)
            CALL Agrif_Update_Variable(vb2b_update_id,locupdate1=(/1,-2/),locupdate2=(/0,-1/),procname = updatevb2b)
#  endif
         ELSE
#  if ! defined DECAL_FEEDBACK
            CALL Agrif_Update_Variable(ub2b_update_id,locupdate=(/0,1/),procname = updateub2b)
            CALL Agrif_Update_Variable(vb2b_update_id,locupdate=(/0,1/),procname = updatevb2b)
#  else
            CALL Agrif_Update_Variable(ub2b_update_id,locupdate1=(/0,1/),locupdate2=(/1,1/),procname = updateub2b)
            CALL Agrif_Update_Variable(vb2b_update_id,locupdate1=(/1,1/),locupdate2=(/0,1/),procname = updatevb2b)
#  endif
         ENDIF
      END IF
# endif
      !
      nbcline = nbcline + 1
      !
      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid = 0.
# if ! defined DECAL_FEEDBACK
      CALL Agrif_Update_Variable(sshn_id,procname = updateSSH)
# else
      CALL Agrif_Update_Variable(sshn_id,locupdate=(/1,0/),procname = updateSSH)
# endif
      Agrif_UseSpecialValueInUpdate = .FALSE.
      ! 
#endif
      !
      ! Do recursive update:
      IF ( lk_agrif_doupd ) THEN ! Initialisation: do recursive update:
         CALL Agrif_ChildGrid_To_ParentGrid()
         CALL Agrif_Update_Dyn()
         CALL Agrif_ParentGrid_To_ChildGrid()
      ENDIF
      !
   END SUBROUTINE Agrif_Update_Dyn

# if defined key_zdftke
   SUBROUTINE Agrif_Update_Tke( kt )
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Update_Tke ***
      !!---------------------------------------------
      !!
      INTEGER, INTENT(in) :: kt
      !       
      IF( ( Agrif_NbStepint() .NE. 0 ) .AND. (kt /= 0) ) RETURN
#  if defined TWO_WAY

      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid = 0.

      CALL Agrif_Update_Variable( en_id, locupdate=(/0,0/), procname=updateEN  )
      CALL Agrif_Update_Variable(avt_id, locupdate=(/0,0/), procname=updateAVT )
      CALL Agrif_Update_Variable(avm_id, locupdate=(/0,0/), procname=updateAVM )

      Agrif_UseSpecialValueInUpdate = .FALSE.

#  endif
      
   END SUBROUTINE Agrif_Update_Tke
# endif /* key_zdftke */

   SUBROUTINE updateTS( tabres, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updateT ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji,jj,jk,jn
      !!---------------------------------------------
      !
      IF (before) THEN
         DO jn = n1,n2
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                     tabres(ji,jj,jk,jn) = tsn(ji,jj,jk,jn)
                  END DO
               END DO
            END DO
         END DO
      ELSE
         IF (.NOT.(lk_agrif_fstep.AND.(neuler==0))) THEN
            ! Add asselin part
            DO jn = n1,n2
               DO jk=k1,k2
                  DO jj=j1,j2
                     DO ji=i1,i2
                        IF( tabres(ji,jj,jk,jn) .NE. 0. ) THEN
                           tsb(ji,jj,jk,jn) = tsb(ji,jj,jk,jn) & 
                                 & + atfp * ( tabres(ji,jj,jk,jn) &
                                 &             - tsn(ji,jj,jk,jn) ) * tmask(ji,jj,jk)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         DO jn = n1,n2
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                     IF( tabres(ji,jj,jk,jn) .NE. 0. ) THEN 
                        tsn(ji,jj,jk,jn) = tabres(ji,jj,jk,jn) * tmask(ji,jj,jk)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      ENDIF
      ! 
   END SUBROUTINE updateTS

   SUBROUTINE updateu( tabres, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updateu ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"
      !!
      INTEGER, INTENT(in) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !! 
      INTEGER :: ji, jj, jk
      REAL(wp) :: zrhoy
      !!---------------------------------------------
      ! 
      IF (before) THEN
         zrhoy = Agrif_Rhoy()
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
                  tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3u_n(ji,jj,jk)
               END DO
            END DO
         END DO
         tabres = zrhoy * tabres
      ELSE
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk) = tabres(ji,jj,jk) / e2u(ji,jj) / fse3u_n(ji,jj,jk)
                  !
                  IF (.NOT.(lk_agrif_fstep.AND.(neuler==0))) THEN ! Add asselin part
                     ub(ji,jj,jk) = ub(ji,jj,jk) & 
                           & + atfp * ( tabres(ji,jj,jk) - un(ji,jj,jk) ) * umask(ji,jj,jk)
                  ENDIF
                  !
                  un(ji,jj,jk) = tabres(ji,jj,jk) * umask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF
      ! 
   END SUBROUTINE updateu

   SUBROUTINE updatev( tabres, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updatev ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"
      !!
      INTEGER :: i1,i2,j1,j2,k1,k2
      INTEGER :: ji,jj,jk
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) :: tabres
      LOGICAL :: before
      !!
      REAL(wp) :: zrhox
      !!---------------------------------------------      
      !
      IF (before) THEN
         zrhox = Agrif_Rhox()
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
                  tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3v_n(ji,jj,jk)
               END DO
            END DO
         END DO
         tabres = zrhox * tabres
      ELSE
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk) = tabres(ji,jj,jk) / e1v(ji,jj) / fse3v_n(ji,jj,jk)
                  !
                  IF (.NOT.(lk_agrif_fstep.AND.(neuler==0))) THEN ! Add asselin part
                     vb(ji,jj,jk) = vb(ji,jj,jk) & 
                           & + atfp * ( tabres(ji,jj,jk) - vn(ji,jj,jk) ) * vmask(ji,jj,jk)
                  ENDIF
                  !
                  vn(ji,jj,jk) = tabres(ji,jj,jk) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF
      ! 
   END SUBROUTINE updatev

   SUBROUTINE updateu2d( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updateu2d ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"
      !!
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !! 
      INTEGER :: ji, jj, jk
      REAL(wp) :: zrhoy
      REAL(wp) :: zcorr
      !!---------------------------------------------
      !
      IF (before) THEN
         zrhoy = Agrif_Rhoy()
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = un_b(ji,jj) * hu(ji,jj) * e2u(ji,jj)
            END DO
         END DO
         tabres = zrhoy * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) =  tabres(ji,jj) * hur(ji,jj) / e2u(ji,jj)  
               !    
               ! Update "now" 3d velocities:
               spgu(ji,jj) = 0.e0
               DO jk=1,jpkm1
                  spgu(ji,jj) = spgu(ji,jj) + fse3u_n(ji,jj,jk) * un(ji,jj,jk)
               END DO
               spgu(ji,jj) = spgu(ji,jj) * hur(ji,jj)
               !
               zcorr = tabres(ji,jj) - spgu(ji,jj)
               DO jk=1,jpkm1              
                  un(ji,jj,jk) = un(ji,jj,jk) + zcorr * umask(ji,jj,jk)           
               END DO
               !
               ! Update barotropic velocities:
#if defined key_dynspg_ts
               IF (.NOT.(lk_agrif_fstep.AND.(neuler==0))) THEN ! Add asselin part
                  zcorr = tabres(ji,jj) - un_b(ji,jj)
                  ub_b(ji,jj) = ub_b(ji,jj) + atfp * zcorr * umask(ji,jj,1)
               END IF
#endif               
               un_b(ji,jj) = tabres(ji,jj) * umask(ji,jj,1)
               !       
               ! Correct "before" velocities to hold correct bt component:
               spgu(ji,jj) = 0.e0
               DO jk=1,jpkm1
                  spgu(ji,jj) = spgu(ji,jj) + fse3u_b(ji,jj,jk) * ub(ji,jj,jk)
               END DO
               spgu(ji,jj) = spgu(ji,jj) * hur_b(ji,jj)
               !
               zcorr = ub_b(ji,jj) - spgu(ji,jj)
               DO jk=1,jpkm1              
                  ub(ji,jj,jk) = ub(ji,jj,jk) + zcorr * umask(ji,jj,jk)           
               END DO
               !
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE updateu2d

   SUBROUTINE updatev2d( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updatev2d ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !! 
      INTEGER :: ji, jj, jk
      REAL(wp) :: zrhox
      REAL(wp) :: zcorr
      !!---------------------------------------------
      !
      IF (before) THEN
         zrhox = Agrif_Rhox()
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = vn_b(ji,jj) * hv(ji,jj) * e1v(ji,jj) 
            END DO
         END DO
         tabres = zrhox * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) =  tabres(ji,jj) * hvr(ji,jj) / e1v(ji,jj)  
               !    
               ! Update "now" 3d velocities:
               spgv(ji,jj) = 0.e0
               DO jk=1,jpkm1
                  spgv(ji,jj) = spgv(ji,jj) + fse3v_n(ji,jj,jk) * vn(ji,jj,jk)
               END DO
               spgv(ji,jj) = spgv(ji,jj) * hvr(ji,jj)
               !
               zcorr = tabres(ji,jj) - spgv(ji,jj)
               DO jk=1,jpkm1              
                  vn(ji,jj,jk) = vn(ji,jj,jk) + zcorr * vmask(ji,jj,jk)           
               END DO
               !
               ! Update barotropic velocities:
#if defined key_dynspg_ts
               IF (.NOT.(lk_agrif_fstep.AND.(neuler==0))) THEN ! Add asselin part
                  zcorr = tabres(ji,jj) - vn_b(ji,jj)
                  vb_b(ji,jj) = vb_b(ji,jj) + atfp * zcorr * vmask(ji,jj,1)
               END IF
#endif               
               vn_b(ji,jj) = tabres(ji,jj) * vmask(ji,jj,1)
               !       
               ! Correct "before" velocities to hold correct bt component:
               spgv(ji,jj) = 0.e0
               DO jk=1,jpkm1
                  spgv(ji,jj) = spgv(ji,jj) + fse3v_b(ji,jj,jk) * vb(ji,jj,jk)
               END DO
               spgv(ji,jj) = spgv(ji,jj) * hvr_b(ji,jj)
               !
               zcorr = vb_b(ji,jj) - spgv(ji,jj)
               DO jk=1,jpkm1              
                  vb(ji,jj,jk) = vb(ji,jj,jk) + zcorr * vmask(ji,jj,jk)           
               END DO
               !
            END DO
         END DO
      ENDIF
      ! 
   END SUBROUTINE updatev2d


   SUBROUTINE updateSSH( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updateSSH ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      !!---------------------------------------------
      ! 
      IF (before) THEN
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = sshn(ji,jj)
            END DO
         END DO
      ELSE
#if ! defined key_dynspg_ts
         IF (.NOT.(lk_agrif_fstep.AND.(neuler==0))) THEN
            DO jj=j1,j2
               DO ji=i1,i2
                  sshb(ji,jj) =   sshb(ji,jj) &
                        & + atfp * ( tabres(ji,jj) - sshn(ji,jj) ) * tmask(ji,jj,1)
               END DO
            END DO
         ENDIF
#endif
         DO jj=j1,j2
            DO ji=i1,i2
               sshn(ji,jj) = tabres(ji,jj) * tmask(ji,jj,1)
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE updateSSH

   SUBROUTINE updateub2b( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updateub2b ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      REAL(wp) :: zrhoy
      !!---------------------------------------------
      !
      IF (before) THEN
         zrhoy = Agrif_Rhoy()
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = ub2_i_b(ji,jj) * e2u(ji,jj)
            END DO
         END DO
         tabres = zrhoy * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               ub2_b(ji,jj) = tabres(ji,jj) / e2u(ji,jj)
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE updateub2b

   SUBROUTINE updatevb2b( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updatevb2b ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      REAL(wp) :: zrhox
      !!---------------------------------------------
      !
      IF (before) THEN
         zrhox = Agrif_Rhox()
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = vb2_i_b(ji,jj) * e1v(ji,jj) 
            END DO
         END DO
         tabres = zrhox * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               vb2_b(ji,jj) = tabres(ji,jj) / e1v(ji,jj)
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE updatevb2b


   SUBROUTINE update_scales( tabres, i1, i2, j1, j2, k1, k2, n1,n2, before )
      ! currently not used
      !!---------------------------------------------
      !!           *** ROUTINE updateT ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"

      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: tabres
      LOGICAL, iNTENT(in) :: before

      INTEGER :: ji,jj,jk
      REAL(wp) :: ztemp

      IF (before) THEN
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk,1) = e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,jk)
                  tabres(ji,jj,jk,2) = e1t(ji,jj)*tmask(ji,jj,jk)
                  tabres(ji,jj,jk,3) = e2t(ji,jj)*tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         tabres(:,:,:,1)=tabres(:,:,:,1)*Agrif_Rhox()*Agrif_Rhoy()
         tabres(:,:,:,2)=tabres(:,:,:,2)*Agrif_Rhox()
         tabres(:,:,:,3)=tabres(:,:,:,3)*Agrif_Rhoy()
      ELSE
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  IF( tabres(ji,jj,jk,1) .NE. 0. ) THEN 
                     print *,'VAL = ',ji,jj,jk,tabres(ji,jj,jk,1),e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,jk)
                     print *,'VAL2 = ',ji,jj,jk,tabres(ji,jj,jk,2),e1t(ji,jj)*tmask(ji,jj,jk)
                     print *,'VAL3 = ',ji,jj,jk,tabres(ji,jj,jk,3),e2t(ji,jj)*tmask(ji,jj,jk)
                     ztemp = sqrt(tabres(ji,jj,jk,1)/(tabres(ji,jj,jk,2)*tabres(ji,jj,jk,3)))
                     print *,'CORR = ',ztemp-1.
                     print *,'NEW VALS = ',tabres(ji,jj,jk,2)*ztemp,tabres(ji,jj,jk,3)*ztemp, &
                           tabres(ji,jj,jk,2)*ztemp*tabres(ji,jj,jk,3)*ztemp
                     e1t(ji,jj) = tabres(ji,jj,jk,2)*ztemp
                     e2t(ji,jj) = tabres(ji,jj,jk,3)*ztemp
                  END IF
               END DO
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE update_scales

# if defined key_zdftke
   SUBROUTINE updateEN( ptab, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updateen ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !!---------------------------------------------
      !
      IF (before) THEN
         ptab (i1:i2,j1:j2,k1:k2) = en(i1:i2,j1:j2,k1:k2)
      ELSE
         en(i1:i2,j1:j2,k1:k2) = ptab (i1:i2,j1:j2,k1:k2) 
      ENDIF
      !
   END SUBROUTINE updateEN


   SUBROUTINE updateAVT( ptab, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updateavt ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !!---------------------------------------------
      !
      IF (before) THEN
         ptab (i1:i2,j1:j2,k1:k2) = avt_k(i1:i2,j1:j2,k1:k2)
      ELSE
         avt_k(i1:i2,j1:j2,k1:k2) = ptab (i1:i2,j1:j2,k1:k2) 
      ENDIF
      !
   END SUBROUTINE updateAVT


   SUBROUTINE updateAVM( ptab, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updateavm ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !!---------------------------------------------
      !
      IF (before) THEN
         ptab (i1:i2,j1:j2,k1:k2) = avm_k(i1:i2,j1:j2,k1:k2)
      ELSE
         avm_k(i1:i2,j1:j2,k1:k2) = ptab (i1:i2,j1:j2,k1:k2) 
      ENDIF
      !
   END SUBROUTINE updateAVM

# endif /* key_zdftke */ 

#else
CONTAINS
   SUBROUTINE agrif_opa_update_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_opa_update_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_opa_update : You should not have seen this print! error?'
   END SUBROUTINE agrif_opa_update_empty
#endif
END MODULE agrif_opa_update

