MODULE agrif_opa_interp
   !!======================================================================
   !!                   ***  MODULE  agrif_opa_interp  ***
   !! AGRIF: interpolation package
   !!======================================================================
   !! History :  2.0  !  2002-06  (XXX)  Original cade
   !!             -   !  2005-11  (XXX) 
   !!            3.2  !  2009-04  (R. Benshila) 
   !!            3.6  !  2014-09  (R. Benshila) 
   !!----------------------------------------------------------------------
#if defined key_agrif && ! defined key_offline
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!   NOT 'key_offline'                               NO off-line tracers
   !!----------------------------------------------------------------------
   !!   Agrif_tra     :
   !!   Agrif_dyn     : 
   !!   interpu       :
   !!   interpv       :
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE dom_oce      
   USE sol_oce
   USE agrif_oce
   USE phycst
   USE in_out_manager
   USE agrif_opa_sponge
   USE lib_mpp
   USE wrk_nemo
   USE dynspg_oce
   USE zdf_oce
 
   IMPLICIT NONE
   PRIVATE

   INTEGER :: bdy_tinterp = 0

   PUBLIC   Agrif_tra, Agrif_dyn, Agrif_ssh, Agrif_dyn_ts, Agrif_ssh_ts, Agrif_dta_ts
   PUBLIC   interpun, interpvn, interpun2d, interpvn2d 
   PUBLIC   interptsn,  interpsshn
   PUBLIC   interpunb, interpvnb, interpub2b, interpvb2b
   PUBLIC   interpe3t, interpumsk, interpvmsk
# if defined key_zdftke
   PUBLIC   Agrif_tke, interpavm
# endif

#  include "domzgr_substitute.h90"  
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/NST 3.6 , NEMO Consortium (2010)
   !! $Id: agrif_opa_interp.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE Agrif_tra
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_tra  ***
      !!----------------------------------------------------------------------
      !
      IF( Agrif_Root() )   RETURN

      Agrif_SpecialValue    = 0.e0
      Agrif_UseSpecialValue = .TRUE.

      CALL Agrif_Bc_variable( tsn_id, procname=interptsn )
      Agrif_UseSpecialValue = .FALSE.
      !
   END SUBROUTINE Agrif_tra


   SUBROUTINE Agrif_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_DYN  ***
      !!----------------------------------------------------------------------  
      !! 
      INTEGER, INTENT(in) ::   kt
      !!
      INTEGER :: ji,jj,jk, j1,j2, i1,i2
      REAL(wp) :: timeref
      REAL(wp) :: z2dt, znugdt
      REAL(wp) :: zrhox, zrhoy
      REAL(wp), POINTER, DIMENSION(:,:) :: spgv1, spgu1
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN

      CALL wrk_alloc( jpi, jpj, spgv1, spgu1 )

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn

      CALL Agrif_Bc_variable(un_interp_id,procname=interpun)
      CALL Agrif_Bc_variable(vn_interp_id,procname=interpvn)

#if defined key_dynspg_flt
      CALL Agrif_Bc_variable(e1u_id,calledweight=1., procname=interpun2d)
      CALL Agrif_Bc_variable(e2v_id,calledweight=1., procname=interpvn2d)
#endif

      Agrif_UseSpecialValue = .FALSE.

      zrhox = Agrif_Rhox()
      zrhoy = Agrif_Rhoy()

      timeref = 1.
      ! time step: leap-frog
      z2dt = 2. * rdt
      ! time step: Euler if restart from rest
      IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
      ! coefficients
      znugdt =  grav * z2dt    

      ! prevent smoothing in ghost cells
      i1=1
      i2=jpi
      j1=1
      j2=jpj
      IF((nbondj == -1).OR.(nbondj == 2)) j1 = 3
      IF((nbondj == +1).OR.(nbondj == 2)) j2 = nlcj-2
      IF((nbondi == -1).OR.(nbondi == 2)) i1 = 3
      IF((nbondi == +1).OR.(nbondi == 2)) i2 = nlci-2


      IF((nbondi == -1).OR.(nbondi == 2)) THEN
#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO jj=j1,j2
               ua(2,jj,jk) = (ua(2,jj,jk) - z2dt * znugdt * laplacu(2,jj))*umask(2,jj,jk)
            END DO
         END DO

         spgu(2,:)=0.

         DO jk=1,jpkm1
            DO jj=1,jpj
               spgu(2,jj)=spgu(2,jj)+fse3u(2,jj,jk)*ua(2,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            IF (umask(2,jj,1).NE.0.) THEN
               spgu(2,jj)=spgu(2,jj)/hu(2,jj)
            ENDIF
         END DO
#else
         spgu(2,:) = ua_b(2,:)
#endif

         DO jk=1,jpkm1
            DO jj=j1,j2
               ua(2,jj,jk) = 0.25*(ua(1,jj,jk)+2.*ua(2,jj,jk)+ua(3,jj,jk))
               ua(2,jj,jk) = ua(2,jj,jk) * umask(2,jj,jk)
            END DO
         END DO

         spgu1(2,:)=0.

         DO jk=1,jpkm1
            DO jj=1,jpj
               spgu1(2,jj)=spgu1(2,jj)+fse3u(2,jj,jk)*ua(2,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            IF (umask(2,jj,1).NE.0.) THEN
               spgu1(2,jj)=spgu1(2,jj)/hu(2,jj)
            ENDIF
         END DO

         DO jk=1,jpkm1
            DO jj=j1,j2
               ua(2,jj,jk) = (ua(2,jj,jk)+spgu(2,jj)-spgu1(2,jj))*umask(2,jj,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgv1(2,:)=0.
         DO jk=1,jpkm1
            DO jj=1,jpj
               spgv1(2,jj)=spgv1(2,jj)+fse3v_a(2,jj,jk)*va(2,jj,jk)
            END DO
         END DO
         DO jj=1,jpj
            spgv1(2,jj)=spgv1(2,jj)*hvr_a(2,jj)
         END DO
         DO jk=1,jpkm1
            DO jj=1,jpj
               va(2,jj,jk) = (va(2,jj,jk)+va_b(2,jj)-spgv1(2,jj))*vmask(2,jj,jk)
            END DO
         END DO
#endif

      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO jj=j1,j2
               ua(nlci-2,jj,jk) = (ua(nlci-2,jj,jk)- z2dt * znugdt * laplacu(nlci-2,jj))*umask(nlci-2,jj,jk)
            END DO
         END DO
         spgu(nlci-2,:)=0.
         DO jk=1,jpkm1
            DO jj=1,jpj
               spgu(nlci-2,jj)=spgu(nlci-2,jj)+fse3u(nlci-2,jj,jk)*ua(nlci-2,jj,jk)
            ENDDO
         ENDDO
         DO jj=1,jpj
            IF (umask(nlci-2,jj,1).NE.0.) THEN
               spgu(nlci-2,jj)=spgu(nlci-2,jj)/hu(nlci-2,jj)
            ENDIF
         END DO
#else
         spgu(nlci-2,:) = ua_b(nlci-2,:)
#endif
         DO jk=1,jpkm1
            DO jj=j1,j2
               ua(nlci-2,jj,jk) = 0.25*(ua(nlci-3,jj,jk)+2.*ua(nlci-2,jj,jk)+ua(nlci-1,jj,jk))

               ua(nlci-2,jj,jk) = ua(nlci-2,jj,jk) * umask(nlci-2,jj,jk)

            END DO
         END DO
         spgu1(nlci-2,:)=0.
         DO jk=1,jpkm1
            DO jj=1,jpj
               spgu1(nlci-2,jj)=spgu1(nlci-2,jj)+fse3u(nlci-2,jj,jk)*ua(nlci-2,jj,jk)*umask(nlci-2,jj,jk)
            END DO
         END DO
         DO jj=1,jpj
            IF (umask(nlci-2,jj,1).NE.0.) THEN
               spgu1(nlci-2,jj)=spgu1(nlci-2,jj)/hu(nlci-2,jj)
            ENDIF
         END DO
         DO jk=1,jpkm1
            DO jj=j1,j2
               ua(nlci-2,jj,jk) = (ua(nlci-2,jj,jk)+spgu(nlci-2,jj)-spgu1(nlci-2,jj))*umask(nlci-2,jj,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgv1(nlci-1,:)=0._wp
         DO jk=1,jpkm1
            DO jj=1,jpj
               spgv1(nlci-1,jj)=spgv1(nlci-1,jj)+fse3v_a(nlci-1,jj,jk)*va(nlci-1,jj,jk)*vmask(nlci-1,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            spgv1(nlci-1,jj)=spgv1(nlci-1,jj)*hvr_a(nlci-1,jj)
         END DO

         DO jk=1,jpkm1
            DO jj=1,jpj
               va(nlci-1,jj,jk) = (va(nlci-1,jj,jk)+va_b(nlci-1,jj)-spgv1(nlci-1,jj))*vmask(nlci-1,jj,jk)
            END DO
         END DO
#endif

      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN

#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,2,jk) = (va(ji,2,jk) - z2dt * znugdt * laplacv(ji,2))*vmask(ji,2,jk)
            END DO
         END DO

         spgv(:,2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv(ji,2)=spgv(ji,2)+fse3v(ji,2,jk)*va(ji,2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,2,1).NE.0.) THEN
               spgv(ji,2)=spgv(ji,2)/hv(ji,2)
            ENDIF
         END DO
#else
         spgv(:,2)=va_b(:,2)
#endif

         DO jk=1,jpkm1
            DO ji=i1,i2
               va(ji,2,jk)=0.25*(va(ji,1,jk)+2.*va(ji,2,jk)+va(ji,3,jk))
               va(ji,2,jk)=va(ji,2,jk)*vmask(ji,2,jk)
            END DO
         END DO

         spgv1(:,2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv1(ji,2)=spgv1(ji,2)+fse3v(ji,2,jk)*va(ji,2,jk)*vmask(ji,2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,2,1).NE.0.) THEN
               spgv1(ji,2)=spgv1(ji,2)/hv(ji,2)
            ENDIF
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,2,jk) = (va(ji,2,jk)+spgv(ji,2)-spgv1(ji,2))*vmask(ji,2,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgu1(:,2)=0._wp
         DO jk=1,jpkm1
            DO ji=1,jpi
               spgu1(ji,2)=spgu1(ji,2)+fse3u_a(ji,2,jk)*ua(ji,2,jk)*umask(ji,2,jk)
            END DO
         END DO

         DO ji=1,jpi
            spgu1(ji,2)=spgu1(ji,2)*hur_a(ji,2)
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               ua(ji,2,jk) = (ua(ji,2,jk)+ua_b(ji,2)-spgu1(ji,2))*umask(ji,2,jk)
            END DO
         END DO
#endif
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN

#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,nlcj-2,jk) = (va(ji,nlcj-2,jk)-z2dt * znugdt * laplacv(ji,nlcj-2))*vmask(ji,nlcj-2,jk)
            END DO
         END DO


         spgv(:,nlcj-2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv(ji,nlcj-2)=spgv(ji,nlcj-2)+fse3v(ji,nlcj-2,jk)*va(ji,nlcj-2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,nlcj-2,1).NE.0.) THEN
               spgv(ji,nlcj-2)=spgv(ji,nlcj-2)/hv(ji,nlcj-2)
            ENDIF
         END DO

#else
         spgv(:,nlcj-2)=va_b(:,nlcj-2)
#endif

         DO jk=1,jpkm1
            DO ji=i1,i2
               va(ji,nlcj-2,jk)=0.25*(va(ji,nlcj-3,jk)+2.*va(ji,nlcj-2,jk)+va(ji,nlcj-1,jk))
               va(ji,nlcj-2,jk) = va(ji,nlcj-2,jk) * vmask(ji,nlcj-2,jk)
            END DO
         END DO

         spgv1(:,nlcj-2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv1(ji,nlcj-2)=spgv1(ji,nlcj-2)+fse3v(ji,nlcj-2,jk)*va(ji,nlcj-2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,nlcj-2,1).NE.0.) THEN
               spgv1(ji,nlcj-2)=spgv1(ji,nlcj-2)/hv(ji,nlcj-2)
            ENDIF
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,nlcj-2,jk) = (va(ji,nlcj-2,jk)+spgv(ji,nlcj-2)-spgv1(ji,nlcj-2))*vmask(ji,nlcj-2,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgu1(:,nlcj-1)=0._wp
         DO jk=1,jpkm1
            DO ji=1,jpi
               spgu1(ji,nlcj-1)=spgu1(ji,nlcj-1)+fse3u_a(ji,nlcj-1,jk)*ua(ji,nlcj-1,jk)
            END DO
         END DO

         DO ji=1,jpi
            spgu1(ji,nlcj-1)=spgu1(ji,nlcj-1)*hur_a(ji,nlcj-1)
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               ua(ji,nlcj-1,jk) = (ua(ji,nlcj-1,jk)+ua_b(ji,nlcj-1)-spgu1(ji,nlcj-1))*umask(ji,nlcj-1,jk)
            END DO
         END DO
#endif

      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, spgv1, spgu1 )
      !
   END SUBROUTINE Agrif_dyn

   SUBROUTINE Agrif_dyn_ts( jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_dyn_ts  ***
      !!----------------------------------------------------------------------  
      !! 
      INTEGER, INTENT(in) ::   jn
      !!
      INTEGER :: ji, jj
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN

      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            va_e(2,jj) = vbdy_w(jj) * hvr_e(2,jj)
            ! Specified fluxes:
            ua_e(2,jj) = ubdy_w(jj) * hur_e(2,jj)
            ! Characteristics method:
            !alt            ua_e(2,jj) = 0.5_wp * ( ubdy_w(jj) * hur_e(2,jj) + ua_e(3,jj) &
            !alt                       &           - sqrt(grav * hur_e(2,jj)) * (sshn_e(3,jj) - hbdy_w(jj)) )
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            va_e(nlci-1,jj) = vbdy_e(jj) * hvr_e(nlci-1,jj)
            ! Specified fluxes:
            ua_e(nlci-2,jj) = ubdy_e(jj) * hur_e(nlci-2,jj)
            ! Characteristics method:
            !alt            ua_e(nlci-2,jj) = 0.5_wp * ( ubdy_e(jj) * hur_e(nlci-2,jj) + ua_e(nlci-3,jj) &
            !alt                            &           + sqrt(grav * hur_e(nlci-2,jj)) * (sshn_e(nlci-2,jj) - hbdy_e(jj)) )
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ua_e(ji,2) = ubdy_s(ji) * hur_e(ji,2)
            ! Specified fluxes:
            va_e(ji,2) = vbdy_s(ji) * hvr_e(ji,2)
            ! Characteristics method:
            !alt            va_e(ji,2) = 0.5_wp * ( vbdy_s(ji) * hvr_e(ji,2) + va_e(ji,3) &
            !alt                       &           - sqrt(grav * hvr_e(ji,2)) * (sshn_e(ji,3) - hbdy_s(ji)) )
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ua_e(ji,nlcj-1) = ubdy_n(ji) * hur_e(ji,nlcj-1)
            ! Specified fluxes:
            va_e(ji,nlcj-2) = vbdy_n(ji) * hvr_e(ji,nlcj-2)
            ! Characteristics method:
            !alt            va_e(ji,nlcj-2) = 0.5_wp * ( vbdy_n(ji) * hvr_e(ji,nlcj-2)  + va_e(ji,nlcj-3) &
            !alt                            &           + sqrt(grav * hvr_e(ji,nlcj-2)) * (sshn_e(ji,nlcj-2) - hbdy_n(ji)) )
         END DO
      ENDIF
      !
   END SUBROUTINE Agrif_dyn_ts

   SUBROUTINE Agrif_dta_ts( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_dta_ts  ***
      !!----------------------------------------------------------------------  
      !! 
      INTEGER, INTENT(in) ::   kt
      !!
      INTEGER :: ji, jj
      LOGICAL :: ll_int_cons
      REAL(wp) :: zrhot, zt
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN

      ll_int_cons = ln_bt_fw ! Assume conservative temporal integration in
      ! the forward case only

      zrhot = Agrif_rhot()

      ! "Central" time index for interpolation:
      IF (ln_bt_fw) THEN
         zt = REAL(Agrif_NbStepint()+0.5_wp,wp) / zrhot
      ELSE
         zt = REAL(Agrif_NbStepint(),wp) / zrhot
      ENDIF

      ! Linear interpolation of sea level
      Agrif_SpecialValue    = 0.e0
      Agrif_UseSpecialValue = .TRUE.
      CALL Agrif_Bc_variable(sshn_id,calledweight=zt, procname=interpsshn )
      Agrif_UseSpecialValue = .FALSE.

      ! Interpolate barotropic fluxes
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn

      IF (ll_int_cons) THEN ! Conservative interpolation
         ! orders matters here !!!!!!
         CALL Agrif_Bc_variable(ub2b_interp_id,calledweight=1._wp, procname=interpub2b) ! Time integrated
         CALL Agrif_Bc_variable(vb2b_interp_id,calledweight=1._wp, procname=interpvb2b)
         bdy_tinterp = 1
         CALL Agrif_Bc_variable(unb_id ,calledweight=1._wp, procname=interpunb) ! After
         CALL Agrif_Bc_variable(vnb_id ,calledweight=1._wp, procname=interpvnb)
         bdy_tinterp = 2
         CALL Agrif_Bc_variable(unb_id ,calledweight=0._wp, procname=interpunb) ! Before
         CALL Agrif_Bc_variable(vnb_id ,calledweight=0._wp, procname=interpvnb)         
      ELSE ! Linear interpolation
         bdy_tinterp = 0
         ubdy_w(:) = 0.e0 ; vbdy_w(:) = 0.e0 
         ubdy_e(:) = 0.e0 ; vbdy_e(:) = 0.e0 
         ubdy_n(:) = 0.e0 ; vbdy_n(:) = 0.e0 
         ubdy_s(:) = 0.e0 ; vbdy_s(:) = 0.e0 
         CALL Agrif_Bc_variable(unb_id,calledweight=zt, procname=interpunb)
         CALL Agrif_Bc_variable(vnb_id,calledweight=zt, procname=interpvnb)
      ENDIF
      Agrif_UseSpecialValue = .FALSE.
      ! 
   END SUBROUTINE Agrif_dta_ts

   SUBROUTINE Agrif_ssh( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_DYN  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   kt
      !!
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN

      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         ssha(2,:)=ssha(3,:)
         sshn(2,:)=sshn(3,:)
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         ssha(nlci-1,:)=ssha(nlci-2,:)
         sshn(nlci-1,:)=sshn(nlci-2,:)
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         ssha(:,2)=ssha(:,3)
         sshn(:,2)=sshn(:,3)
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         ssha(:,nlcj-1)=ssha(:,nlcj-2)
         sshn(:,nlcj-1)=sshn(:,nlcj-2)
      ENDIF

   END SUBROUTINE Agrif_ssh

   SUBROUTINE Agrif_ssh_ts( jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_ssh_ts  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   jn
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            ssha_e(2,jj) = hbdy_w(jj)
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            ssha_e(nlci-1,jj) = hbdy_e(jj)
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ssha_e(ji,2) = hbdy_s(ji)
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ssha_e(ji,nlcj-1) = hbdy_n(ji)
         END DO
      ENDIF

   END SUBROUTINE Agrif_ssh_ts

# if defined key_zdftke
   SUBROUTINE Agrif_tke
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_tke  ***
      !!----------------------------------------------------------------------  
      REAL(wp) ::   zalpha
      !
      zalpha = REAL( Agrif_NbStepint() + Agrif_IRhot() - 1, wp ) / REAL( Agrif_IRhot(), wp )
      IF( zalpha > 1. )   zalpha = 1.
      
      Agrif_SpecialValue    = 0.e0
      Agrif_UseSpecialValue = .TRUE.
      
      CALL Agrif_Bc_variable(avm_id ,calledweight=zalpha, procname=interpavm)       
              
      Agrif_UseSpecialValue = .FALSE.
      !
   END SUBROUTINE Agrif_tke
# endif

   SUBROUTINE interptsn(ptab,i1,i2,j1,j2,k1,k2,n1,n2,before,nb,ndir)
      !!---------------------------------------------
      !!   *** ROUTINE interptsn ***
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
         ptab(i1:i2,j1:j2,k1:k2,n1:n2) = tsn(i1:i2,j1:j2,k1:k2,n1:n2)
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
            DO jn = 1, jpts
               tsa(nlci,j1:j2,k1:k2,jn) = zalpha1 * ptab(nlci,j1:j2,k1:k2,jn) + zalpha2 * ptab(nlci-1,j1:j2,k1:k2,jn)
               DO jk = 1, jpkm1
                  DO jj = jmin,jmax
                     IF( umask(nlci-2,jj,jk) == 0.e0 ) THEN
                        tsa(nlci-1,jj,jk,jn) = tsa(nlci,jj,jk,jn) * tmask(nlci-1,jj,jk)
                     ELSE
                        tsa(nlci-1,jj,jk,jn)=(zalpha4*tsa(nlci,jj,jk,jn)+zalpha3*tsa(nlci-2,jj,jk,jn))*tmask(nlci-1,jj,jk)
                        IF( un(nlci-2,jj,jk) > 0.e0 ) THEN
                           tsa(nlci-1,jj,jk,jn)=( zalpha6*tsa(nlci-2,jj,jk,jn)+zalpha5*tsa(nlci,jj,jk,jn) & 
                                 + zalpha7*tsa(nlci-3,jj,jk,jn) ) * tmask(nlci-1,jj,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            ENDDO
         ENDIF
         ! 
         IF( northern_side ) THEN            
            DO jn = 1, jpts
               tsa(i1:i2,nlcj,k1:k2,jn) = zalpha1 * ptab(i1:i2,nlcj,k1:k2,jn) + zalpha2 * ptab(i1:i2,nlcj-1,k1:k2,jn)
               DO jk = 1, jpkm1
                  DO ji = imin,imax
                     IF( vmask(ji,nlcj-2,jk) == 0.e0 ) THEN
                        tsa(ji,nlcj-1,jk,jn) = tsa(ji,nlcj,jk,jn) * tmask(ji,nlcj-1,jk)
                     ELSE
                        tsa(ji,nlcj-1,jk,jn)=(zalpha4*tsa(ji,nlcj,jk,jn)+zalpha3*tsa(ji,nlcj-2,jk,jn))*tmask(ji,nlcj-1,jk)        
                        IF (vn(ji,nlcj-2,jk) > 0.e0 ) THEN
                           tsa(ji,nlcj-1,jk,jn)=( zalpha6*tsa(ji,nlcj-2,jk,jn)+zalpha5*tsa(ji,nlcj,jk,jn)  &
                                 + zalpha7*tsa(ji,nlcj-3,jk,jn) ) * tmask(ji,nlcj-1,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            ENDDO
         ENDIF
         !
         IF( western_side) THEN            
            DO jn = 1, jpts
               tsa(1,j1:j2,k1:k2,jn) = zalpha1 * ptab(1,j1:j2,k1:k2,jn) + zalpha2 * ptab(2,j1:j2,k1:k2,jn)
               DO jk = 1, jpkm1
                  DO jj = jmin,jmax
                     IF( umask(2,jj,jk) == 0.e0 ) THEN
                        tsa(2,jj,jk,jn) = tsa(1,jj,jk,jn) * tmask(2,jj,jk)
                     ELSE
                        tsa(2,jj,jk,jn)=(zalpha4*tsa(1,jj,jk,jn)+zalpha3*tsa(3,jj,jk,jn))*tmask(2,jj,jk)        
                        IF( un(2,jj,jk) < 0.e0 ) THEN
                           tsa(2,jj,jk,jn)=(zalpha6*tsa(3,jj,jk,jn)+zalpha5*tsa(1,jj,jk,jn)+zalpha7*tsa(4,jj,jk,jn))*tmask(2,jj,jk)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
         !
         IF( southern_side ) THEN           
            DO jn = 1, jpts
               tsa(i1:i2,1,k1:k2,jn) = zalpha1 * ptab(i1:i2,1,k1:k2,jn) + zalpha2 * ptab(i1:i2,2,k1:k2,jn)
               DO jk=1,jpk      
                  DO ji=imin,imax
                     IF( vmask(ji,2,jk) == 0.e0 ) THEN
                        tsa(ji,2,jk,jn)=tsa(ji,1,jk,jn) * tmask(ji,2,jk)
                     ELSE
                        tsa(ji,2,jk,jn)=(zalpha4*tsa(ji,1,jk,jn)+zalpha3*tsa(ji,3,jk,jn))*tmask(ji,2,jk)
                        IF( vn(ji,2,jk) < 0.e0 ) THEN
                           tsa(ji,2,jk,jn)=(zalpha6*tsa(ji,3,jk,jn)+zalpha5*tsa(ji,1,jk,jn)+zalpha7*tsa(ji,4,jk,jn))*tmask(ji,2,jk)
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
            tsa(nlci-1,2,:,:) = ptab(nlci-1,2,:,:)
         ENDIF
         ! East north
         IF ((eastern_side).AND.((nbondj == 1).OR.(nbondj == 2))) THEN
            tsa(nlci-1,nlcj-1,:,:) = ptab(nlci-1,nlcj-1,:,:)
         ENDIF
         ! West south
         IF ((western_side).AND.((nbondj == -1).OR.(nbondj == 2))) THEN
            tsa(2,2,:,:) = ptab(2,2,:,:)
         ENDIF
         ! West north
         IF ((western_side).AND.((nbondj == 1).OR.(nbondj == 2))) THEN
            tsa(2,nlcj-1,:,:) = ptab(2,nlcj-1,:,:)
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE interptsn

   SUBROUTINE interpsshn(ptab,i1,i2,j1,j2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      LOGICAL :: western_side, eastern_side,northern_side,southern_side
      !!----------------------------------------------------------------------  
      !
      IF( before) THEN
         ptab(i1:i2,j1:j2) = sshn(i1:i2,j1:j2)
      ELSE
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
         IF(western_side)  hbdy_w(j1:j2) = ptab(i1,j1:j2) * tmask(i1,j1:j2,1)
         IF(eastern_side)  hbdy_e(j1:j2) = ptab(i1,j1:j2) * tmask(i1,j1:j2,1)
         IF(southern_side) hbdy_s(i1:i2) = ptab(i1:i2,j1) * tmask(i1:i2,j1,1)
         IF(northern_side) hbdy_n(i1:i2) = ptab(i1:i2,j1) * tmask(i1:i2,j1,1)
      ENDIF
      !
   END SUBROUTINE interpsshn

   SUBROUTINE interpun(ptab,i1,i2,j1,j2,k1,k2, before)
      !!---------------------------------------------
      !!   *** ROUTINE interpun ***
      !!---------------------------------------------    
      !!
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji,jj,jk
      REAL(wp) :: zrhoy	
      !!---------------------------------------------    
      !
      IF (before) THEN 
         DO jk=1,jpk
            DO jj=j1,j2
               DO ji=i1,i2
                  ptab(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
                  ptab(ji,jj,jk) = ptab(ji,jj,jk) * fse3u(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE
         zrhoy = Agrif_Rhoy()
         DO jk=1,jpkm1
            DO jj=j1,j2
               ua(i1:i2,jj,jk) = (ptab(i1:i2,jj,jk)/(zrhoy*e2u(i1:i2,jj)))
               ua(i1:i2,jj,jk) = ua(i1:i2,jj,jk) / fse3u(i1:i2,jj,jk)
            END DO
         END DO
      ENDIF
      ! 
   END SUBROUTINE interpun


   SUBROUTINE interpun2d(ptab,i1,i2,j1,j2,before)
      !!---------------------------------------------
      !!   *** ROUTINE interpun ***
      !!---------------------------------------------    
      !
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !
      INTEGER :: ji,jj
      REAL(wp) :: ztref
      REAL(wp) :: zrhoy	
      !!---------------------------------------------    
      !
      ztref = 1.

      IF (before) THEN 
         DO jj=j1,j2
            DO ji=i1,MIN(i2,nlci-1)
               ptab(ji,jj) = e2u(ji,jj) * ((gcx(ji+1,jj) - gcx(ji,jj))/e1u(ji,jj)) 
            END DO
         END DO
      ELSE
         zrhoy = Agrif_Rhoy()
         DO jj=j1,j2
            laplacu(i1:i2,jj) = ztref * (ptab(i1:i2,jj)/(zrhoy*e2u(i1:i2,jj))) !*umask(i1:i2,jj,1)
         END DO
      ENDIF
      ! 
   END SUBROUTINE interpun2d


   SUBROUTINE interpvn(ptab,i1,i2,j1,j2,k1,k2, before)
      !!---------------------------------------------
      !!   *** ROUTINE interpvn ***
      !!---------------------------------------------    
      !
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !
      INTEGER :: ji,jj,jk
      REAL(wp) :: zrhox	
      !!---------------------------------------------    
      !      
      IF (before) THEN          
         !interpv entre 1 et k2 et interpv2d en jpkp1
         DO jk=k1,jpk
            DO jj=j1,j2
               DO ji=i1,i2
                  ptab(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
                  ptab(ji,jj,jk) = ptab(ji,jj,jk) * fse3v(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE          
         zrhox= Agrif_Rhox()
         DO jk=1,jpkm1
            DO jj=j1,j2
               va(i1:i2,jj,jk) = (ptab(i1:i2,jj,jk)/(zrhox*e1v(i1:i2,jj)))
               va(i1:i2,jj,jk) = va(i1:i2,jj,jk) / fse3v(i1:i2,jj,jk)
            END DO
         END DO
      ENDIF
      !        
   END SUBROUTINE interpvn

   SUBROUTINE interpvn2d(ptab,i1,i2,j1,j2,before)
      !!---------------------------------------------
      !!   *** ROUTINE interpvn ***
      !!---------------------------------------------    
      !
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !
      INTEGER :: ji,jj
      REAL(wp) :: zrhox	
      REAL(wp) :: ztref
      !!---------------------------------------------    
      ! 
      ztref = 1.    
      IF (before) THEN 
         !interpv entre 1 et k2 et interpv2d en jpkp1
         DO jj=j1,MIN(j2,nlcj-1)
            DO ji=i1,i2
               ptab(ji,jj) = e1v(ji,jj) * ((gcx(ji,jj+1) - gcx(ji,jj))/e2v(ji,jj)) * vmask(ji,jj,1)
            END DO
         END DO
      ELSE           
         zrhox = Agrif_Rhox()
         DO ji=i1,i2
            laplacv(ji,j1:j2) = ztref * (ptab(ji,j1:j2)/(zrhox*e1v(ji,j1:j2)))
         END DO
      ENDIF
      !      
   END SUBROUTINE interpvn2d

   SUBROUTINE interpunb(ptab,i1,i2,j1,j2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpunb  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !!
      INTEGER :: ji,jj
      REAL(wp) :: zrhoy, zrhot, zt0, zt1, ztcoeff
      LOGICAL :: western_side, eastern_side,northern_side,southern_side
      !!----------------------------------------------------------------------  
      !
      IF (before) THEN 
         DO jj=j1,j2
            DO ji=i1,i2
               ptab(ji,jj) = un_b(ji,jj) * e2u(ji,jj) * hu(ji,jj) 
            END DO
         END DO
      ELSE
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
         zrhoy = Agrif_Rhoy()
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot      
         ! Polynomial interpolation coefficients:
         IF( bdy_tinterp == 1 ) THEN
            ztcoeff = zrhot * (  zt1**2._wp * (       zt1 - 1._wp)        &
                  &      - zt0**2._wp * (       zt0 - 1._wp)        )
         ELSEIF( bdy_tinterp == 2 ) THEN
            ztcoeff = zrhot * (  zt1        * (       zt1 - 1._wp)**2._wp &
                  &      - zt0        * (       zt0 - 1._wp)**2._wp ) 

         ELSE
            ztcoeff = 1
         ENDIF
         !   
         IF(western_side) THEN
            ubdy_w(j1:j2) = ubdy_w(j1:j2) + ztcoeff * ptab(i1,j1:j2)  
         ENDIF
         IF(eastern_side) THEN
            ubdy_e(j1:j2) = ubdy_e(j1:j2) + ztcoeff * ptab(i1,j1:j2)  
         ENDIF
         IF(southern_side) THEN
            ubdy_s(i1:i2) = ubdy_s(i1:i2) + ztcoeff * ptab(i1:i2,j1) 
         ENDIF
         IF(northern_side) THEN
            ubdy_n(i1:i2) = ubdy_n(i1:i2) + ztcoeff * ptab(i1:i2,j1) 
         ENDIF
         !            
         IF( bdy_tinterp == 0 .OR. bdy_tinterp == 2) THEN
            IF(western_side) THEN
               ubdy_w(j1:j2) = ubdy_w(j1:j2) / (zrhoy*e2u(i1,j1:j2))   &
                     &                                  * umask(i1,j1:j2,1)
            ENDIF
            IF(eastern_side) THEN
               ubdy_e(j1:j2) = ubdy_e(j1:j2) / (zrhoy*e2u(i1,j1:j2))   &
                     &                                  * umask(i1,j1:j2,1)
            ENDIF
            IF(southern_side) THEN
               ubdy_s(i1:i2) = ubdy_s(i1:i2) / (zrhoy*e2u(i1:i2,j1))   &
                     &                                  * umask(i1:i2,j1,1)
            ENDIF
            IF(northern_side) THEN
               ubdy_n(i1:i2) = ubdy_n(i1:i2) / (zrhoy*e2u(i1:i2,j1))   &
                     &                                  * umask(i1:i2,j1,1)
            ENDIF
         ENDIF
      ENDIF
      ! 
   END SUBROUTINE interpunb

   SUBROUTINE interpvnb(ptab,i1,i2,j1,j2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpvnb  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !!
      INTEGER :: ji,jj
      REAL(wp) :: zrhox, zrhot, zt0, zt1, ztcoeff   
      LOGICAL :: western_side, eastern_side,northern_side,southern_side
      !!----------------------------------------------------------------------  
      ! 
      IF (before) THEN 
         DO jj=j1,j2
            DO ji=i1,i2
               ptab(ji,jj) = vn_b(ji,jj) * e1v(ji,jj) * hv(ji,jj) 
            END DO
         END DO
      ELSE
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
         zrhox = Agrif_Rhox()
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot      
         IF( bdy_tinterp == 1 ) THEN
            ztcoeff = zrhot * (  zt1**2._wp * (       zt1 - 1._wp)        &
                  &      - zt0**2._wp * (       zt0 - 1._wp)        )
         ELSEIF( bdy_tinterp == 2 ) THEN
            ztcoeff = zrhot * (  zt1        * (       zt1 - 1._wp)**2._wp &
                  &      - zt0        * (       zt0 - 1._wp)**2._wp ) 

         ELSE
            ztcoeff = 1
         ENDIF
         !
         IF(western_side) THEN
            vbdy_w(j1:j2) = vbdy_w(j1:j2) + ztcoeff * ptab(i1,j1:j2)  
         ENDIF
         IF(eastern_side) THEN
            vbdy_e(j1:j2) = vbdy_e(j1:j2) + ztcoeff * ptab(i1,j1:j2)  
         ENDIF
         IF(southern_side) THEN
            vbdy_s(i1:i2) = vbdy_s(i1:i2) + ztcoeff * ptab(i1:i2,j1)
         ENDIF
         IF(northern_side) THEN
            vbdy_n(i1:i2) = vbdy_n(i1:i2) + ztcoeff * ptab(i1:i2,j1) 
         ENDIF
         !            
         IF( bdy_tinterp == 0 .OR. bdy_tinterp == 2) THEN
            IF(western_side) THEN
               vbdy_w(j1:j2) = vbdy_w(j1:j2) / (zrhox*e1v(i1,j1:j2))   &
                     &                                  * vmask(i1,j1:j2,1)
            ENDIF
            IF(eastern_side) THEN
               vbdy_e(j1:j2) = vbdy_e(j1:j2) / (zrhox*e1v(i1,j1:j2))   &
                     &                                  * vmask(i1,j1:j2,1)
            ENDIF
            IF(southern_side) THEN
               vbdy_s(i1:i2) = vbdy_s(i1:i2) / (zrhox*e1v(i1:i2,j1))   &
                     &                                  * vmask(i1:i2,j1,1)
            ENDIF
            IF(northern_side) THEN
               vbdy_n(i1:i2) = vbdy_n(i1:i2) / (zrhox*e1v(i1:i2,j1))   &
                     &                                  * vmask(i1:i2,j1,1)
            ENDIF
         ENDIF
      ENDIF
      !
   END SUBROUTINE interpvnb

   SUBROUTINE interpub2b(ptab,i1,i2,j1,j2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpub2b  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !!
      INTEGER :: ji,jj
      REAL(wp) :: zrhot, zt0, zt1,zat
      LOGICAL :: western_side, eastern_side,northern_side,southern_side
      !!----------------------------------------------------------------------  
      IF( before ) THEN
         DO jj=j1,j2
            DO ji=i1,i2
               ptab(ji,jj) = ub2_b(ji,jj) * e2u(ji,jj)
            END DO
         END DO
      ELSE
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot
         ! Polynomial interpolation coefficients:
         zat = zrhot * (  zt1**2._wp * (-2._wp*zt1 + 3._wp)        &
               &      - zt0**2._wp * (-2._wp*zt0 + 3._wp)        ) 
         ! 
         IF(western_side ) ubdy_w(j1:j2) = zat * ptab(i1,j1:j2)  
         IF(eastern_side ) ubdy_e(j1:j2) = zat * ptab(i1,j1:j2)  
         IF(southern_side) ubdy_s(i1:i2) = zat * ptab(i1:i2,j1) 
         IF(northern_side) ubdy_n(i1:i2) = zat * ptab(i1:i2,j1) 
      ENDIF
      ! 
   END SUBROUTINE interpub2b

   SUBROUTINE interpvb2b(ptab,i1,i2,j1,j2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpvb2b  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !!
      INTEGER :: ji,jj
      REAL(wp) :: zrhot, zt0, zt1,zat
      LOGICAL :: western_side, eastern_side,northern_side,southern_side
      !!----------------------------------------------------------------------  
      !
      IF( before ) THEN
         DO jj=j1,j2
            DO ji=i1,i2
               ptab(ji,jj) = vb2_b(ji,jj) * e1v(ji,jj)
            END DO
         END DO
      ELSE      
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot
         ! Polynomial interpolation coefficients:
         zat = zrhot * (  zt1**2._wp * (-2._wp*zt1 + 3._wp)        &
               &      - zt0**2._wp * (-2._wp*zt0 + 3._wp)        ) 
         !
         IF(western_side ) vbdy_w(j1:j2) = zat * ptab(i1,j1:j2)  
         IF(eastern_side ) vbdy_e(j1:j2) = zat * ptab(i1,j1:j2)  
         IF(southern_side) vbdy_s(i1:i2) = zat * ptab(i1:i2,j1) 
         IF(northern_side) vbdy_n(i1:i2) = zat * ptab(i1:i2,j1) 
      ENDIF
      !      
   END SUBROUTINE interpvb2b

   SUBROUTINE interpe3t(ptab,i1,i2,j1,j2,k1,k2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpe3t  ***
      !!----------------------------------------------------------------------  
      ! 
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL :: before
      INTEGER, INTENT(in) :: nb , ndir
      !
      INTEGER :: ji, jj, jk
      LOGICAL :: western_side, eastern_side, northern_side, southern_side
      REAL(wp) :: ztmpmsk      
      !!----------------------------------------------------------------------  
      !    
      IF (before) THEN
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  ptab(ji,jj,jk) = tmask(ji,jj,jk) * e3t_0(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  ! Get velocity mask at boundary edge points:
                  IF (western_side)  ztmpmsk = umask(ji    ,jj    ,1)
                  IF (eastern_side)  ztmpmsk = umask(nlci-2,jj    ,1)
                  IF (northern_side) ztmpmsk = vmask(ji    ,nlcj-2,1)
                  IF (southern_side) ztmpmsk = vmask(ji    ,2     ,1)

                  IF (ABS(ptab(ji,jj,jk) - tmask(ji,jj,jk) * e3t_0(ji,jj,jk))*ztmpmsk > 1.D-2) THEN
                     IF (western_side) THEN
                        WRITE(numout,*) 'ERROR bathymetry merge at the western border ji,jj,jk ', ji+nimpp-1,jj+njmpp-1,jk
                     ELSEIF (eastern_side) THEN
                        WRITE(numout,*) 'ERROR bathymetry merge at the eastern border ji,jj,jk ', ji+nimpp-1,jj+njmpp-1,jk
                     ELSEIF (southern_side) THEN
                        WRITE(numout,*) 'ERROR bathymetry merge at the southern border ji,jj,jk', ji+nimpp-1,jj+njmpp-1,jk
                     ELSEIF (northern_side) THEN
                        WRITE(numout,*) 'ERROR bathymetry merge at the northen border ji,jj,jk', ji+nimpp-1,jj+njmpp-1,jk
                     ENDIF
                     WRITE(numout,*) '      ptab(ji,jj,jk), fse3t(ji,jj,jk) ', ptab(ji,jj,jk), e3t_0(ji,jj,jk)
                     kindic_agr = kindic_agr + 1
                  ENDIF
               END DO
            END DO
         END DO

      ENDIF
      ! 
   END SUBROUTINE interpe3t

   SUBROUTINE interpumsk(ptab,i1,i2,j1,j2,k1,k2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpumsk  ***
      !!----------------------------------------------------------------------  
      ! 
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL :: before
      INTEGER, INTENT(in) :: nb , ndir
      !
      INTEGER :: ji, jj, jk
      LOGICAL :: western_side, eastern_side   
      !!----------------------------------------------------------------------  
      !    
      IF (before) THEN
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  ptab(ji,jj,jk) = umask(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE

         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                   ! Velocity mask at boundary edge points:
                  IF (ABS(ptab(ji,jj,jk) - umask(ji,jj,jk)) > 1.D-2) THEN
                     IF (western_side) THEN
                        WRITE(numout,*) 'ERROR with umask at the western border ji,jj,jk ', ji+nimpp-1,jj+njmpp-1,jk
                        WRITE(numout,*) '      masks: parent, child ', ptab(ji,jj,jk), umask(ji,jj,jk)
                        kindic_agr = kindic_agr + 1
                     ELSEIF (eastern_side) THEN
                        WRITE(numout,*) 'ERROR with umask at the eastern border ji,jj,jk ', ji+nimpp-1,jj+njmpp-1,jk
                        WRITE(numout,*) '      masks: parent, child ', ptab(ji,jj,jk), umask(ji,jj,jk)
                        kindic_agr = kindic_agr + 1
                     ENDIF
                  ENDIF
               END DO
            END DO
         END DO

      ENDIF
      ! 
   END SUBROUTINE interpumsk

   SUBROUTINE interpvmsk(ptab,i1,i2,j1,j2,k1,k2,before,nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpvmsk  ***
      !!----------------------------------------------------------------------  
      ! 
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL :: before
      INTEGER, INTENT(in) :: nb , ndir
      !
      INTEGER :: ji, jj, jk
      LOGICAL :: northern_side, southern_side     
      !!----------------------------------------------------------------------  
      !    
      IF (before) THEN
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  ptab(ji,jj,jk) = vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE

         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                   ! Velocity mask at boundary edge points:
                  IF (ABS(ptab(ji,jj,jk) - vmask(ji,jj,jk)) > 1.D-2) THEN
                     IF (southern_side) THEN
                        WRITE(numout,*) 'ERROR with vmask at the southern border ji,jj,jk ', ji+nimpp-1,jj+njmpp-1,jk
                        WRITE(numout,*) '      masks: parent, child ', ptab(ji,jj,jk), vmask(ji,jj,jk)
                        kindic_agr = kindic_agr + 1
                     ELSEIF (northern_side) THEN
                        WRITE(numout,*) 'ERROR with vmask at the northern border ji,jj,jk ', ji+nimpp-1,jj+njmpp-1,jk
                        WRITE(numout,*) '      masks: parent, child ', ptab(ji,jj,jk), vmask(ji,jj,jk)
                        kindic_agr = kindic_agr + 1
                     ENDIF
                  ENDIF
               END DO
            END DO
         END DO

      ENDIF
      ! 
   END SUBROUTINE interpvmsk

# if defined key_zdftke

   SUBROUTINE interpavm(ptab,i1,i2,j1,j2,k1,k2,before)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interavm  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !!----------------------------------------------------------------------  
      !      
      IF( before) THEN
         ptab (i1:i2,j1:j2,k1:k2) = avm_k(i1:i2,j1:j2,k1:k2)
      ELSE
         avm_k(i1:i2,j1:j2,k1:k2) = ptab (i1:i2,j1:j2,k1:k2)
      ENDIF
      !
   END SUBROUTINE interpavm

# endif /* key_zdftke */

#else
   !!----------------------------------------------------------------------
   !!   Empty module                                          no AGRIF zoom
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE Agrif_OPA_Interp_empty
      WRITE(*,*)  'agrif_opa_interp : You should not have seen this print! error?'
   END SUBROUTINE Agrif_OPA_Interp_empty
#endif

   !!======================================================================
END MODULE agrif_opa_interp
