MODULE agrif_lim2_interp
   !!======================================================================
   !!                       ***  MODULE agrif_lim2_update ***
   !! Nesting module :  update surface ocean boundary condition over ice
   !!                   from a child grif
   !! Sea-Ice model  :  LIM 2.0 Sea ice model time-stepping
   !!======================================================================
   !! History :  2.0   !  04-2008  (F. Dupont)  initial version
   !!            3.4   !  09-2012  (R. Benshila, C. Herbaut) update and EVP
   !!----------------------------------------------------------------------
#if defined key_agrif && defined key_lim2 
   !!----------------------------------------------------------------------
   !!   'key_lim2'  :                                 LIM 2.0 sea-ice model
   !!   'key_agrif' :                                 AGRIF library
   !!----------------------------------------------------------------------
   !!   agrif_interp_lim2   : update sea-ice model on boundaries or total
   !!                         sea-ice area
   !!  agrif_rhg_lim2_load  : interpolcation of ice velocities using Agrif
   !!  agrif_rhg_lim2       : sub-interpolation of ice velocities for both 
   !!                         splitting time and sea-ice time step
   !!  agrif_interp_u_ice   : atomic routine to interpolate u_ice 
   !!  agrif_interp_u_ice   : atomic routine to interpolate v_ice 
   !!  agrif_trp_lim2_load  : interpolcation of ice properties using Agrif
   !!  agrif_trp_lim2       : sub-interpolation of ice properties for  
   !!                         sea-ice time step
   !!  agrif_interp_u_ice   : atomic routine to interpolate ice properties 
   !!----------------------------------------------------------------------
   USE par_oce
   USE dom_oce
   USE sbc_oce
   USE ice_2
   USE dom_ice_2
   USE agrif_ice

   IMPLICIT NONE
   PRIVATE

   PUBLIC agrif_rhg_lim2_load, agrif_rhg_lim2
   PUBLIC agrif_trp_lim2_load, agrif_trp_lim2
   PUBLIC interp_u_ice, interp_v_ice
   PUBLIC interp_adv_ice

   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PRIVATE :: uice_agr, vice_agr
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: tabice_agr 


   !!----------------------------------------------------------------------
   !! NEMO/NST 3.4 , NEMO Consortium (2012)
   !! $Id: agrif_lim2_interp.F90 6794 2016-07-06 13:10:40Z jchanut $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

# if defined key_lim2_vp
   SUBROUTINE agrif_rhg_lim2_load
      !!-----------------------------------------------------------------------
      !!              *** ROUTINE agrif_rhg_lim2_load ***
      !!
      !!  ** Method  : need a special routine for dealing with exchanging data
      !! between the child and parent grid during ice step
      !!
      !!-----------------------------------------------------------------------
      !
      IF (Agrif_Root()) RETURN

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .FALSE.
      u_ice_nst(:,:) = 0.
      v_ice_nst(:,:) = 0.
      CALL Agrif_Bc_variable( u_ice_id ,procname=interp_u_ice, calledweight=1. )
      CALL Agrif_Bc_variable( v_ice_id ,procname=interp_v_ice, calledweight=1. )
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .FALSE.
      !
   END SUBROUTINE agrif_rhg_lim2_load


   SUBROUTINE agrif_rhg_lim2(pu_n,pv_n)
      !!-----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_rhg_lim2 ***
      !!
      !!  ** Method  : we feel the boundaries with values stored above
      !!-----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,0:jpj+1), INTENT(inout) :: pu_n, pv_n
      !!
      REAL(wp) :: zrhox, zrhoy
      INTEGER :: ji,jj
      !!-----------------------------------------------------------------------
      !
      IF (Agrif_Root()) RETURN

      zrhox = Agrif_Rhox()
      zrhoy = Agrif_Rhoy()

      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=2,jpj
            pu_n(3,jj) = u_ice_nst(3,jj)/(zrhoy*e2f(2,jj-1))*tmu(3,jj)
         END DO
         DO jj=2,jpj
            pv_n(3,jj) = v_ice_nst(3,jj)/(zrhox*e1f(2,jj-1))*tmu(3,jj)
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=2,jpj
            pu_n(nlci-1,jj) = u_ice_nst(nlci-1,jj)/(zrhoy*e2f(nlci-2,jj-1))*tmu(nlci-1,jj)
         END DO
         DO jj=2,jpj
            pv_n(nlci-1,jj) = v_ice_nst(nlci-1,jj)/(zrhox*e1f(nlci-2,jj-1))*tmu(nlci-1,jj)
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=2,jpi
            pv_n(ji,3) = v_ice_nst(ji,3)/(zrhox*e1f(ji-1,2))*tmu(ji,3)
         END DO
         DO ji=2,jpi
            pu_n(ji,3) = u_ice_nst(ji,3)/(zrhoy*e2f(ji-1,2))*tmu(ji,3)
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=2,jpi
            pv_n(ji,nlcj-1) = v_ice_nst(ji,nlcj-1)/(zrhox*e1f(ji-1,nlcj-2))*tmu(ji,nlcj-1)
         END DO
         DO ji=2,jpi
            pu_n(ji,nlcj-1) = u_ice_nst(ji,nlcj-1)/(zrhoy*e2f(ji-1,nlcj-2))*tmu(ji,nlcj-1)
         END DO
      ENDIF
      !
   END SUBROUTINE agrif_rhg_lim2

#else
   SUBROUTINE agrif_rhg_lim2_load
      !!-----------------------------------------------------------------------
      !!              *** ROUTINE agrif_rhg_lim2_load ***
      !!
      !!  ** Method  : need a special routine for dealing with exchanging data
      !!  between the child and parent grid during ice step
      !!               we interpolate and store the boundary if needed, ie if
      !!  we are in inside a new parent ice time step
      !!-----------------------------------------------------------------------
      INTEGER :: ji,jj
      REAL(wp) :: zrhox, zrhoy
      !!-----------------------------------------------------------------------
      !
      IF (Agrif_Root()) RETURN

      IF( lim_nbstep == 1. ) THEN
         !
         ! switch old values by hand
         u_ice_oe(:,:,1) =  u_ice_oe(:,:,2)
         v_ice_oe(:,:,1) =  v_ice_oe(:,:,2)
         u_ice_sn(:,:,1) =  u_ice_sn(:,:,2)
         v_ice_sn(:,:,1) =  v_ice_sn(:,:,2)
         ! interpolation of boundaries (called weight prevents AGRIF interpolation)
         Agrif_SpecialValue=-9999.
         Agrif_UseSpecialValue = .TRUE.
         IF( .NOT. ALLOCATED(uice_agr) )THEN
            ALLOCATE(uice_agr(jpi,jpj), vice_agr(jpi,jpj))
         ENDIF
         uice_agr = 0.
         vice_agr = 0.
         CALL Agrif_Bc_variable(u_ice_id,procname=interp_u_ice, calledweight=1.)
         CALL Agrif_Bc_variable(v_ice_id,procname=interp_v_ice, calledweight=1.)
         Agrif_SpecialValue=0.
         Agrif_UseSpecialValue = .FALSE.
         !  
         zrhox = agrif_rhox() ;    zrhoy = agrif_rhoy()      
         uice_agr(:,:) =  uice_agr(:,:)/(zrhoy*e2u(:,:))*umask(:,:,1)
         vice_agr(:,:) =  vice_agr(:,:)/(zrhox*e1v(:,:))*vmask(:,:,1)
         ! fill  boundaries
         DO jj = 1, jpj
            DO ji = 1, 2
               u_ice_oe(ji,  jj,2) = uice_agr(ji       ,jj) 
               u_ice_oe(ji+2,jj,2) = uice_agr(nlci+ji-3,jj)
            END DO
         END DO
         DO jj = 1, jpj
            v_ice_oe(2,jj,2) = vice_agr(2     ,jj) 
            v_ice_oe(4,jj,2) = vice_agr(nlci-1,jj)
         END DO
         DO ji = 1, jpi
            u_ice_sn(ji,2,2) = uice_agr(ji,2     ) 
            u_ice_sn(ji,4,2) = uice_agr(ji,nlcj-1)
         END DO
         DO jj = 1, 2
            DO ji = 1, jpi
               v_ice_sn(ji,jj  ,2) = vice_agr(ji,jj       ) 
               v_ice_sn(ji,jj+2,2) = vice_agr(ji,nlcj+jj-3)
            END DO
         END DO
         !
      ENDIF
      !
   END SUBROUTINE agrif_rhg_lim2_load


   SUBROUTINE agrif_rhg_lim2( kiter, kitermax, cd_type)
      !!-----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_rhg_lim2  ***
      !!
      !!  ** Method  : simple call to atomic routines using stored values to
      !!  fill the boundaries depending of the position of the point and 
      !!  computing factor for time interpolation 
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kiter, kitermax
      CHARACTER(len=1), INTENT( in ) :: cd_type
      !!    
      REAL(wp) :: zalpha, zbeta
      !!-----------------------------------------------------------------------
      !
      IF (Agrif_Root()) RETURN
      zalpha = REAL(lim_nbstep,wp) / (Agrif_Rhot()*Agrif_PArent(nn_fsbc)/REAL(nn_fsbc))
      zbeta  = REAL(kiter,wp) / kitermax
      zbeta = zalpha * zbeta
      SELECT CASE(cd_type)
      CASE('U')
         CALL ParcoursU( zbeta )
      CASE('V')
         CALL ParcoursV( zbeta )
      END SELECT
      !
   END SUBROUTINE agrif_rhg_lim2


   SUBROUTINE ParcoursU( pbeta )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE parcoursU ***
      !!
      !!  ** Method  : time and spatial interpolation for U-point using values
      !!  interpolated from the coarse grid and inside dvalues     
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pbeta
      !!
      INTEGER :: ji, jj
      !!-----------------------------------------------------------------------
      !
      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            DO ji=1,2
               u_ice(ji,jj) = (1-pbeta) * u_ice_oe(ji,jj,1) + pbeta * u_ice_oe(ji,jj,2)
            END DO
         END DO
         DO jj=1,jpj
            u_ice(2,jj) = 0.25*(u_ice(1,jj)+2.*u_ice(2,jj)+u_ice(3,jj))
            u_ice(2,jj) = u_ice(2,jj) * umask(2,jj,1)
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            DO ji=1,2
               u_ice(nlci+ji-3,jj) = (1-pbeta) * u_ice_oe(ji+2,jj,1) + pbeta * u_ice_oe(ji+2,jj,2)
            END DO
         END DO
         DO jj=1,jpj
            u_ice(nlci-2,jj) = 0.25*(u_ice(nlci-3,jj)+2.*u_ice(nlci-2,jj)+u_ice(nlci-1,jj))
            u_ice(nlci-2,jj) = u_ice(nlci-2,jj) * umask(nlci-2,jj,1)
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            u_ice(ji,2) = (1-pbeta) * u_ice_sn(ji,2,1) + pbeta * u_ice_sn(ji,2,2)
            u_ice(ji,2) = u_ice(ji,2)*umask(ji,2,1)
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            u_ice(ji,nlcj-1) = (1-pbeta) * u_ice_sn(ji,4,1) + pbeta * u_ice_sn(ji,4,2)
            u_ice(ji,nlcj-1) = u_ice(ji,nlcj-1)*umask(ji,nlcj-1,1)
         END DO
      ENDIF
      !
   END SUBROUTINE ParcoursU


   SUBROUTINE ParcoursV( pbeta )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE parcoursV ***
      !! 
      !!  ** Method  : time and spatial interpolation for V-point using values
      !!  interpolated from the coarse grid and inside dvalues     
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pbeta
      !!
      INTEGER :: ji, jj
      !!-----------------------------------------------------------------------
      !
      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            v_ice(2,jj) = (1-pbeta) * v_ice_oe(2,jj,1) + pbeta * v_ice_oe(2,jj,2)
            v_ice(2,jj) = v_ice(2,jj) * vmask(2,jj,1)
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            v_ice(nlci-1,jj) = (1-pbeta) * v_ice_oe(4,jj,1) + pbeta * v_ice_oe(4,jj,2)
            v_ice(nlci-1,jj) = v_ice(nlci-1,jj)*vmask(nlci-1,jj,1)
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO jj=1,2
            DO ji=1,jpi
               v_ice(ji,jj) = (1-pbeta) * v_ice_sn(ji,jj,1) + pbeta * v_ice_sn(ji,jj,2)
            END DO
         END DO
         DO ji=1,jpi
            v_ice(ji,2)=0.25*(v_ice(ji,1)+2.*v_ice(ji,2)+v_ice(ji,3))
            v_ice(ji,2)=v_ice(ji,2)*vmask(ji,2,1)
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO jj=1,2
            DO ji=1,jpi
               v_ice(ji,nlcj+jj-3) = (1-pbeta) * v_ice_sn(ji,jj+2,1) + pbeta * v_ice_sn(ji,jj+2,2)
            END DO
         END DO
         DO ji=1,jpi
            v_ice(ji,nlcj-2)=0.25*(v_ice(ji,nlcj-3)+2.*v_ice(ji,nlcj-2)+v_ice(ji,nlcj-1))
            v_ice(ji,nlcj-2) = v_ice(ji,nlcj-2) * vmask(ji,nlcj-2,1)
         END DO
      ENDIF
      !
   END SUBROUTINE ParcoursV
# endif
   SUBROUTINE agrif_trp_lim2_load
      !!-----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_trp_lim2_load ***
      !!
      !!  ** Method  : need a special routine for dealing with exchanging data
      !!  between the child and parent grid during ice step
      !!               we interpolate and store the boundary if needed, ie if
      !!  we are in inside a new parent ice time step
     !!-----------------------------------------------------------------------
      INTEGER :: ji,jj,jn
      !!-----------------------------------------------------------------------
      !
      IF (Agrif_Root()) RETURN
      IF( lim_nbstep == 1. ) THEN
         !
         ! switch old values
         adv_ice_oe(:,:,:,1) =  adv_ice_oe(:,:,:,2)
         adv_ice_sn(:,:,:,1) =  adv_ice_sn(:,:,:,2)
         ! interpolation of boundaries
         IF(.NOT.ALLOCATED(tabice_agr))THEN
            ALLOCATE(tabice_agr(jpi,jpj,7))   
         ENDIF
         tabice_agr(:,:,:) = 0.
         Agrif_SpecialValue=-9999.
         Agrif_UseSpecialValue = .TRUE.
         CALL Agrif_Bc_variable( adv_ice_id ,procname=interp_adv_ice,calledweight=1. )
         Agrif_SpecialValue=0.
         Agrif_UseSpecialValue = .FALSE.
         !  
         ! fill  boundaries
         DO jn =1,7
            DO jj = 1, jpj
               DO ji=1,2
                  adv_ice_oe(ji  ,jj,jn,2) = tabice_agr(ji       ,jj,jn) 
                  adv_ice_oe(ji+2,jj,jn,2) = tabice_agr(nlci-2+ji,jj,jn)
               END DO
            END DO
         END DO

         Do jn =1,7
            Do jj =1,2
               DO ji = 1, jpi
                  adv_ice_sn(ji,jj  ,jn,2) = tabice_agr(ji,jj       ,jn) 
                  adv_ice_sn(ji,jj+2,jn,2) = tabice_agr(ji,nlcj-2+jj,jn)
               END DO
            END DO
         END DO
         !
      ENDIF
      !
   END SUBROUTINE agrif_trp_lim2_load


   SUBROUTINE agrif_trp_lim2
      !!-----------------------------------------------------------------------
      !!                  *** ROUTINE agrif_trp_lim2 ***
      !!
      !!  ** Method  : time coefficient and call to atomic routines
      !!-----------------------------------------------------------------------
      INTEGER :: ji,jj,jn
      REAL(wp) :: zalpha
      !!-----------------------------------------------------------------------      
      !
      IF (Agrif_Root()) RETURN

      zalpha = REAL(lim_nbstep,wp) / (Agrif_Rhot()*Agrif_PArent(nn_fsbc)/REAL(nn_fsbc))
      !
      tabice_agr(:,:,:) = 0.e0
      DO jn =1,7
         DO jj =1,2
            DO ji = 1, jpi
               tabice_agr(ji,jj        ,jn) = (1-zalpha)*adv_ice_sn(ji,jj  ,jn,1) + zalpha*adv_ice_sn(ji,jj  ,jn,2) 
               tabice_agr(ji,nlcj-2+jj ,jn) = (1-zalpha)*adv_ice_sn(ji,jj+2,jn,1) + zalpha*adv_ice_sn(ji,jj+2,jn,2) 
            END DO
         END DO
      END DO

      DO jn =1,7
         DO jj = 1, jpj
            DO ji=1,2
               tabice_agr(ji       ,jj,jn) = (1-zalpha)*adv_ice_oe(ji  ,jj,jn,1) + zalpha*adv_ice_oe(ji  ,jj,jn,2) 
               tabice_agr(nlci-2+ji,jj,jn) = (1-zalpha)*adv_ice_oe(ji+2,jj,jn,1) + zalpha*adv_ice_oe(ji+2,jj,jn,2) 
            END DO
         END DO
      END DO
      !
      CALL parcoursT( tabice_agr(:,:, 1), frld  )
      CALL parcoursT( tabice_agr(:,:, 2), hicif )
      CALL parcoursT( tabice_agr(:,:, 3), hsnif )
      CALL parcoursT( tabice_agr(:,:, 4), tbif(:,:,1) )
      CALL parcoursT( tabice_agr(:,:, 5), tbif(:,:,2) )
      CALL parcoursT( tabice_agr(:,:, 6), tbif(:,:,3) )
      CALL parcoursT( tabice_agr(:,:, 7), qstoif )
      !
   END SUBROUTINE agrif_trp_lim2


   SUBROUTINE parcoursT ( pinterp, pfinal )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE parcoursT ***
      !!
      !!  ** Method  : fill boundaries for T points
      !!-----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pinterp
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pfinal
      !!
      REAL(wp) :: zbound, zvbord
      REAL(wp), DIMENSION(jpi,jpj) ::  zui_u, zvi_v
      INTEGER :: ji, jj
      !!-----------------------------------------------------------------------
      !
      zui_u = 0.e0
      zvi_v = 0.e0
      ! zvbord factor between 1 and 2 to take into account slip or no-slip boundary conditions.
      zbound=0.
      zvbord = 1.0 + ( 1.0 - zbound )
#if defined key_lim2_vp
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
            zui_u(ji,jj) = ( u_ice(ji+1,jj  ) + u_ice(ji+1,jj+1) ) / ( MAX( tmu(ji+1,jj  ) + tmu(ji+1,jj+1), zvbord ) )
            zvi_v(ji,jj) = ( v_ice(ji  ,jj+1) + v_ice(ji+1,jj+1) ) / ( MAX( tmu(ji  ,jj+1) + tmu(ji+1,jj+1), zvbord ) )
         END DO
      END DO
#else
      zui_u(:,:) = u_ice(:,:)
      zvi_v(:,:) = v_ice(:,:)
#endif

      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            !            IF (zui_u(2,jj).EQ.0.) THEN
            !               pfinal (2,jj) = pfinal (1,jj) * tms(2,jj)
            !            ELSE
            pfinal(2,jj) = 0.25* pinterp(1,jj) + 0.5 *  pinterp(2,jj) + 0.25 *pfinal(3,jj)
            !            ENDIF
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            !            IF (zvi_v(ji,2).EQ.0.) THEN
            !               pfinal (ji,2) = pfinal (ji,1) * tms(ji,2)
            !            ELSE
            pfinal(ji,2) = 0.25* pinterp(ji,1) + 0.5 *  pinterp(ji,2) + 0.25 *pfinal(ji,3)
            !            ENDIF 
         END DO
      ENDIF


      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            !            IF (zui_u(nlci-2,jj).EQ.0.) THEN
            !               pfinal(nlci-1,jj) = pfinal (nlci,jj) * tms(nlci-1,jj)
            !            ELSE
            pfinal(nlci-1,jj) = 0.25* pinterp(nlci,jj) + 0.5 *  pinterp(nlci-1,jj) + 0.25 *pfinal(nlci-2,jj)
            !           ENDIF
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            !            IF (zvi_v(ji,nlcj-2).EQ.0.) THEN
            !               pfinal (ji,nlcj-1) =  pfinal(ji,nlcj) * tms(ji,nlcj-1)
            !            ELSE
            pfinal(ji,nlcj-1) = 0.25* pinterp(ji,nlcj) + 0.5 *  pinterp(ji,nlcj-1) + 0.25 *pfinal(ji,nlcj-2)
            !            ENDIF
         END DO
      ENDIF


      pfinal (:,:) = pfinal (:,:) * tms(:,:)
      ! 
   END SUBROUTINE parcoursT


   SUBROUTINE interp_u_ice( tabres, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                     *** ROUTINE interp_u_ice ***
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji,jj
      !!-----------------------------------------------------------------------
      !
#if defined key_lim2_vp
      IF( before ) THEN
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               IF( tmu(ji,jj) == 0. ) THEN
                  tabres(ji,jj) = -9999.
               ELSE
                  tabres(ji,jj) = e2f(ji-1,jj-1) * u_ice(ji,jj)
               ENDIF
            END DO
         END DO
      ELSE
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               uice_agr(ji,jj) = tabres(ji,jj)
            END DO
         END DO
      ENDIF
#else
      IF( before ) THEN
         DO jj= j1, j2
            DO ji= i1, i2
               IF( umask(ji,jj,1) == 0. ) THEN
                  tabres(ji,jj) = -9999.
               ELSE
                  tabres(ji,jj) = e2u(ji,jj) * u_ice(ji,jj)
               ENDIF
            END DO
         END DO
      ELSE
         DO jj= j1, j2
            DO ji= i1, i2
               uice_agr(ji,jj) = tabres(ji,jj)
            END DO
         END DO
      ENDIF
#endif
   END SUBROUTINE interp_u_ice


   SUBROUTINE interp_v_ice( tabres, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE interp_v_ice ***
      !!-----------------------------------------------------------------------      
      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj
      !!-----------------------------------------------------------------------
      !
#if defined key_lim2_vp
      IF( before ) THEN
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               IF( tmu(ji,jj) == 0. ) THEN
                  tabres(ji,jj) = -9999.
               ELSE
                  tabres(ji,jj) = e1f(ji-1,jj-1) * v_ice(ji,jj)
               ENDIF
            END DO
         END DO
      ELSE
         DO jj=MAX(j1,2),j2
            DO ji=MAX(i1,2),i2
               vice_agr(ji,jj) = tabres(ji,jj)
            END DO
         END DO
      ENDIF   
#else
      IF( before ) THEN
	      DO jj= j1 ,j2
	         DO ji = i1, i2
	            IF( vmask(ji,jj,1) == 0. ) THEN
	               tabres(ji,jj) = -9999.
	            ELSE
	               tabres(ji,jj) = e1v(ji,jj) * v_ice(ji,jj)
	            ENDIF
	         END DO
	      END DO
      ELSE
	      DO jj= j1 ,j2
	         DO ji = i1, i2
	            vice_agr(ji,jj) = tabres(ji,jj)
	         END DO
	      END DO
      ENDIF
#endif
   END SUBROUTINE interp_v_ice


   SUBROUTINE interp_adv_ice( tabres, i1, i2, j1, j2, k1, k2, before )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE interp_adv_ice ***                           
      !!
      !! ** Purpose : fill an array with  ice variables
      !!              to be advected
      !!              put -9999 where no ice for correct extrapolation             
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj, jk
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         DO jj=j1,j2
	    DO ji=i1,i2
	       IF( tms(ji,jj) == 0. ) THEN
	          tabres(ji,jj,:) = -9999 
	       ELSE
	          tabres(ji,jj, 1) = frld  (ji,jj)
	          tabres(ji,jj, 2) = hicif (ji,jj)
	          tabres(ji,jj, 3) = hsnif (ji,jj)
	          tabres(ji,jj, 4) = tbif  (ji,jj,1)
	          tabres(ji,jj, 5) = tbif  (ji,jj,2)
	          tabres(ji,jj, 6) = tbif  (ji,jj,3)
	          tabres(ji,jj, 7) = qstoif(ji,jj)
	       ENDIF
	    END DO
         END DO
      ELSE
	 DO jj=j1,j2
	    DO ji=i1,i2
               DO jk=k1, k2
	          tabice_agr(ji,jj,jk) = tabres(ji,jj,jk)
               END DO
	    END DO
	 END DO
      ENDIF
      !
   END SUBROUTINE interp_adv_ice

#else
CONTAINS
   SUBROUTINE agrif_lim2_interp_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_lim2_interp_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_lim2_interp : You should not have seen this print! error?'
   END SUBROUTINE agrif_lim2_interp_empty
#endif
END MODULE agrif_lim2_interp
