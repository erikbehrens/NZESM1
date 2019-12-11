MODULE agrif_top_interp
#if defined key_agrif && defined key_top
   USE par_oce
   USE oce
   USE dom_oce      
   USE sol_oce
   USE agrif_oce
   USE agrif_top_sponge
   USE par_trc
   USE trc
   USE lib_mpp
   USE wrk_nemo  

   IMPLICIT NONE
   PRIVATE

  !CEB PUBLIC Agrif_trc, interptrn
    PUBLIC Agrif_trc_ceb, interptrn
#  include "domzgr_substitute.h90"  
#  include "vectopt_loop_substitute.h90"
  !!----------------------------------------------------------------------
   !! NEMO/NST 3.6 , NEMO Consortium (2010)
   !! $Id: agrif_top_interp.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   CONTAINS

   !CEB SUBROUTINE Agrif_trc
   SUBROUTINE Agrif_trc_ceb
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_trc  ***
      !!----------------------------------------------------------------------
      !
      IF( Agrif_Root() )   RETURN

      Agrif_SpecialValue    = 0.e0
      Agrif_UseSpecialValue = .TRUE.

      CALL Agrif_Bc_variable( trn_id, procname=interptrn )
      Agrif_UseSpecialValue = .FALSE.
      !
   END SUBROUTINE Agrif_trc_ceb
   !/CEB END SUBROUTINE Agrif_trc
   SUBROUTINE interptrn(ptab,i1,i2,j1,j2,k1,k2,n1,n2,before,nb,ndir)
      !!---------------------------------------------
      !!   *** ROUTINE interptrn ***
      !!---------------------------------------------
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: ptab
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER :: imin, imax, jmin, jmax
      REAL(wp) ::   zrhox , zalpha1, zalpha2, zalpha3
      REAL(wp) ::   zalpha4, zalpha5, zalpha6, zalpha7
      LOGICAL :: western_side, eastern_side,northern_side,southern_side

      IF (before) THEN         
         ptab(i1:i2,j1:j2,k1:k2,n1:n2) = trn(i1:i2,j1:j2,k1:k2,n1:n2)
      ELSE
         !
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
         !
         zrhox = Agrif_Rhox()
         ! 
         zalpha1 = ( zrhox - 1. ) * 0.5
         zalpha2 = 1. - zalpha1
         ! 
         zalpha3 = ( zrhox - 1. ) / ( zrhox + 1. )
         zalpha4 = 1. - zalpha3
         ! 
         zalpha6 = 2. * ( zrhox - 1. ) / ( zrhox + 1. )
         zalpha7 =    - ( zrhox - 1. ) / ( zrhox + 3. )
         zalpha5 = 1. - zalpha6 - zalpha7
         !
         imin = i1
         imax = i2
         jmin = j1
         jmax = j2
         ! 
         ! Remove CORNERS
         IF((nbondj == -1).OR.(nbondj == 2)) jmin = 3
         IF((nbondj == +1).OR.(nbondj == 2)) jmax = nlcj-2
         IF((nbondi == -1).OR.(nbondi == 2)) imin = 3
         IF((nbondi == +1).OR.(nbondi == 2)) imax = nlci-2        
         !
         IF( eastern_side) THEN
            DO jn = 1, jptra
               tra(nlci,j1:j2,k1:k2,jn) = zalpha1 * ptab(nlci,j1:j2,k1:k2,jn) + zalpha2 * ptab(nlci-1,j1:j2,k1:k2,jn)
               DO jk = 1, jpkm1
                  DO jj = jmin,jmax
                     IF( umask(nlci-2,jj,jk) == 0.e0 ) THEN
                        tra(nlci-1,jj,jk,jn) = tra(nlci,jj,jk,jn) * tmask(nlci-1,jj,jk)
                     ELSE
                        tra(nlci-1,jj,jk,jn)=(zalpha4*tra(nlci,jj,jk,jn)+zalpha3*tra(nlci-2,jj,jk,jn))*tmask(nlci-1,jj,jk)
                        IF( un(nlci-2,jj,jk) > 0.e0 ) THEN
                           tra(nlci-1,jj,jk,jn)=( zalpha6*tra(nlci-2,jj,jk,jn)+zalpha5*tra(nlci,jj,jk,jn) & 
                                 + zalpha7*tra(nlci-3,jj,jk,jn) ) * tmask(nlci-1,jj,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            ENDDO
         ENDIF
         ! 
         IF( northern_side ) THEN            
            DO jn = 1, jptra
               tra(i1:i2,nlcj,k1:k2,jn) = zalpha1 * ptab(i1:i2,nlcj,k1:k2,jn) + zalpha2 * ptab(i1:i2,nlcj-1,k1:k2,jn)
               DO jk = 1, jpkm1
                  DO ji = imin,imax
                     IF( vmask(ji,nlcj-2,jk) == 0.e0 ) THEN
                        tra(ji,nlcj-1,jk,jn) = tra(ji,nlcj,jk,jn) * tmask(ji,nlcj-1,jk)
                     ELSE
                        tra(ji,nlcj-1,jk,jn)=(zalpha4*tra(ji,nlcj,jk,jn)+zalpha3*tra(ji,nlcj-2,jk,jn))*tmask(ji,nlcj-1,jk)        
                        IF (vn(ji,nlcj-2,jk) > 0.e0 ) THEN
                           tra(ji,nlcj-1,jk,jn)=( zalpha6*tra(ji,nlcj-2,jk,jn)+zalpha5*tra(ji,nlcj,jk,jn)  &
                                 + zalpha7*tra(ji,nlcj-3,jk,jn) ) * tmask(ji,nlcj-1,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            ENDDO
         ENDIF
         !
         IF( western_side) THEN            
            DO jn = 1, jptra
               tra(1,j1:j2,k1:k2,jn) = zalpha1 * ptab(1,j1:j2,k1:k2,jn) + zalpha2 * ptab(2,j1:j2,k1:k2,jn)
               DO jk = 1, jpkm1
                  DO jj = jmin,jmax
                     IF( umask(2,jj,jk) == 0.e0 ) THEN
                        tra(2,jj,jk,jn) = tra(1,jj,jk,jn) * tmask(2,jj,jk)
                     ELSE
                        tra(2,jj,jk,jn)=(zalpha4*tra(1,jj,jk,jn)+zalpha3*tra(3,jj,jk,jn))*tmask(2,jj,jk)        
                        IF( un(2,jj,jk) < 0.e0 ) THEN
                           tra(2,jj,jk,jn)=(zalpha6*tra(3,jj,jk,jn)+zalpha5*tra(1,jj,jk,jn)+zalpha7*tra(4,jj,jk,jn))*tmask(2,jj,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
         !
         IF( southern_side ) THEN           
            DO jn = 1, jptra
               tra(i1:i2,1,k1:k2,jn) = zalpha1 * ptab(i1:i2,1,k1:k2,jn) + zalpha2 * ptab(i1:i2,2,k1:k2,jn)
               DO jk=1,jpk      
                  DO ji=imin,imax
                     IF( vmask(ji,2,jk) == 0.e0 ) THEN
                        tra(ji,2,jk,jn)=tra(ji,1,jk,jn) * tmask(ji,2,jk)
                     ELSE
                        tra(ji,2,jk,jn)=(zalpha4*tra(ji,1,jk,jn)+zalpha3*tra(ji,3,jk,jn))*tmask(ji,2,jk)
                        IF( vn(ji,2,jk) < 0.e0 ) THEN
                           tra(ji,2,jk,jn)=(zalpha6*tra(ji,3,jk,jn)+zalpha5*tra(ji,1,jk,jn)+zalpha7*tra(ji,4,jk,jn))*tmask(ji,2,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            ENDDO
         ENDIF
         !
         ! Treatment of corners
         ! 
         ! East south
         IF ((eastern_side).AND.((nbondj == -1).OR.(nbondj == 2))) THEN
            tra(nlci-1,2,:,:) = ptab(nlci-1,2,:,:)
         ENDIF
         ! East north
         IF ((eastern_side).AND.((nbondj == 1).OR.(nbondj == 2))) THEN
            tra(nlci-1,nlcj-1,:,:) = ptab(nlci-1,nlcj-1,:,:)
         ENDIF
         ! West south
         IF ((western_side).AND.((nbondj == -1).OR.(nbondj == 2))) THEN
            tra(2,2,:,:) = ptab(2,2,:,:)
         ENDIF
         ! West north
         IF ((western_side).AND.((nbondj == 1).OR.(nbondj == 2))) THEN
            tra(2,nlcj-1,:,:) = ptab(2,nlcj-1,:,:)
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE interptrn

#else
CONTAINS
   SUBROUTINE Agrif_TOP_Interp_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_Top_Interp_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_top_interp : You should not have seen this print! error?'
   END SUBROUTINE Agrif_TOP_Interp_empty
#endif
END MODULE agrif_top_interp
