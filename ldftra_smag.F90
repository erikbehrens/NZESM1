MODULE ldftra_smag
   !!======================================================================
   !!                     ***  MODULE  ldftrasmag  ***
   !! Ocean physics:  variable eddy induced velocity coefficients
   !!======================================================================
#if   defined key_traldf_smag   &&   defined key_traldf_c3d
   !!----------------------------------------------------------------------
   !!   'key_traldf_smag'      and           smagorinsky  diffusivity
   !!   'key_traldf_c3d'                    3D tracer lateral  mixing coef.
   !!----------------------------------------------------------------------
   !!   ldf_eiv      : compute the eddy induced velocity coefficients
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcrnf          ! river runoffs
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE phycst          ! physical constants
   USE ldfslp          ! iso-neutral slopes
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom
   USE wrk_nemo
   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC ldf_tra_smag               ! routine called by step.F90
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Id$
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------

CONTAINS

 



  !!----------------------------------------------------------------------
  !!                        ***  ldf_tra_smag.F90  ***
  !!----------------------------------------------------------------------


   SUBROUTINE ldf_tra_smag( kt )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_smag  ***
      !!                   
      !! ** Purpose :   initializations of the horizontal ocean physics
      !!
      !! ** Method  :   3D eddy viscosity coef. 
      !!    M.Griffies, R.Hallberg AMS, 2000
      !!  for laplacian:
      !!   Asmag=(C/pi)^2*dx*dy sqrt(D^2), C=1 for tracers, C=3-4 for viscosity
      !!   D^2= rm_smsh*(du/dx-dv/dy)^2+(dv/dx+du/dy)^2 for Cartesian coordinates
      !!  IF rm_smsh = 0 , only shear is used, recommended for tidal flows
      !!  in general case du/dx ==> e2 d(u/e2)/dx;  du/dy ==> e1 d(u/e1)/dy; 
      !!                  dv/dx ==> e2 d(v/e2)/dx;  dv/dy ==> e1 d(v/e1)/dy
      !!  for bilaplacian: now this option is deleted as unstable or non-conservative
      !!   - delta{Bsmag (delta(T)} = -Bsmag* delta{delta(T)} - delta(Bsmag)*delta( T )
      !!  second term is of arbitrary sign on the edge of fronts and can induce instability 
      !!   Bsmag=Asmag*dx*dy/8
      !!  
      !!       laplacian operator   : ahm1, ahm2 defined at T- and F-points
      !!                              ahm3, ahm4 never used
      !!       bilaplacian operator : ahm1, ahm2 never used
      !!                           :  ahm3, ahm4 defined at U- and V-points
      !!       ??? explanation of the default is missing
      !!  last modified : Maria Luneva, October 2012
      !!----------------------------------------------------------------------
      !!
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl
      REAL ( wp), POINTER , DIMENSION (:,:) :: zux, zvx , zuy , zvy 
      REAL ( wp), POINTER , DIMENSION (:,:) :: zue1, zue2 , zve1 , zve2 
      INTEGER, INTENT( in )                 ::   kt                             ! ocean time-step inedx
      !! * Arguments
      INTEGER                               :: ji,jj,jk

      REAL (wp)                             ::     zdeltau, zdeltav, zhsmag ,zsmsh    ! temporary scalars
      
      CALL wrk_alloc (jpi,jpj,zux, zvx , zuy , zvy )
      CALL wrk_alloc (jpi,jpj,zue1, zue2 , zve1 , zve2 )
      !!----------------------------------------------------------------------
      IF(  kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' ldf_tra_smag : 3D eddy smagorinsky diffusivity '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~   --  '
         IF(lwp) WRITE(numout,*) '               Coefficients are computed'
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)
      ENDIF

      zhsmag = rn_chsmag
      zsmsh  = rn_smsh 
      zux(:,:)=0._wp ; zuy(:,:)=0._wp ; zvx(:,:)=0._wp ; zvy(:,:)=0._wp

      ! -------------------
      ahtt(:,:,:) = rn_aht_0
       IF( ln_traldf_bilap ) THEN
        IF( lwp .AND. kt == nit000) WRITE(numout,* )'ldf_tra_smag :no bilaplacian Smagorinsky diffusivity'
        IF( lwp .AND. kt == nit000) WRITE(numout,* )'ldf_tra_smag :bilaplacian diffusivity set to constant'  
       ENDIF



      ! harmonic operator   (U-, V-, W-points)
      ! ----------------- 

      ahtu(:,:,:) = rn_aht_0                  ! set ahtu , ahtv at u- and v-points,
      ahtv(:,:,:) = rn_aht_0                  ! and ahtw at w-point
      ahtw(:,:,:) = rn_aht_0                  ! (here example: no space variation)
      
      IF( ln_traldf_lap ) THEN

         DO jk=1,jpk

           zue2(:,:)=un(:,:,jk)/e2u(:,:)
           zve1(:,:)=vn(:,:,jk)/e1v(:,:)
           zue1(:,:)=un(:,:,jk)/e1u(:,:)
           zve2(:,:)=vn(:,:,jk)/e2v(:,:)


           DO jj=2,jpj
            DO ji=2,jpi
            zux(ji,jj)=(zue2(ji,jj)-zue2(ji-1,jj))/e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,jk) * zsmsh 
            zvy(ji,jj)=(zve1(ji,jj)-zve1(ji,jj-1))/e2t(ji,jj)*e1t(ji,jj)*tmask(ji,jj,jk) * zsmsh 
            ENDDO
           ENDDO

           DO jj=1,jpjm1
            DO ji=1,jpim1
            zuy(ji,jj)=(zue1(ji,jj+1)-zue1(ji,jj))/e2f(ji,jj)*e1f(ji,jj)*fmask(ji,jj,jk)
            zvx(ji,jj)=(zve2(ji+1,jj)-zve2(ji,jj))/e1f(ji,jj)*e2f(ji,jj)*fmask(ji,jj,jk)
            ENDDO
           ENDDO


          DO jj=2,jpjm1
           DO ji=2,jpim1
           zdeltau=2._wp/( e1u(ji,jj)**(-2)+e2u(ji,jj)**(-2) )
           zdeltav=2._wp/( e1v(ji,jj)**(-2)+e2v(ji,jj)**(-2) )

           ahtu(ji,jj,jk)=MAX( rn_aht_0 , (zhsmag/rpi)**2*zdeltau*                                &
                          SQRT(0.25_wp*( zux(ji,jj)+zux(ji+1,jj)-zvy(ji,jj)-zvy(ji+1,jj) )**2+    &
                               0.25_wp*( zuy(ji,jj)+zuy(ji,jj-1)+zvx(ji,jj)+zvx(ji,jj-1)  )**2) )

           ahtv(ji,jj,jk)=MAX( rn_aht_0 ,  (zhsmag/rpi)**2*zdeltav*                               &
                          SQRT(0.25_wp*( zux(ji,jj)+zux(ji,jj+1)-zvy(ji,jj)-zvy(ji,jj+1) )**2+    &
                               0.25_wp*( zuy(ji,jj)+zuy(ji-1,jj)+zvx(ji-1,jj)+zvx(ji,jj)  )**2) )


         !!! stability criteria: aht<delta**2/(4*dt)   dt=2*rdt , positiveness require aht<delta**2/(8*dt)
             ahtu(ji,jj,jk)=MIN(ahtu(ji,jj,jk),zdeltau/(16*rdt) ,rn_aht_m)
             ahtv(ji,jj,jk)=MIN(ahtv(ji,jj,jk),zdeltav/(16*rdt) ,rn_aht_m)
         ! so...


            ENDDO
           ENDDO
         ENDDO
        ENDIF
            ahtu(:,:,jpk) = ahtu(:,:,jpkm1)
            ahtv(:,:,jpk) = ahtv(:,:,jpkm1)
        CALL lbc_lnk( ahtu, 'U', 1. )   ! Lateral boundary conditions
        CALL lbc_lnk( ahtv, 'V', 1. )

IF(  kt == nit000 ) THEN

      IF(lwp ) THEN                    ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: ahtu at k = 1'
         CALL prihre( ahtu(:,:,1), jpi, jpj, 1, jpi, 1,   &
            &                                1, jpj, 1, 1.e-1, numout )
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: ahtv at k = 1'
         CALL prihre( ahtv(:,:,1), jpi, jpj, 1, jpi, 1,   &
            &                                1, jpj, 1, 1.e-1, numout )
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: ahtw at k = 1'
         CALL prihre( ahtw(:,:,1), jpi, jpj, 1, jpi, 1,   &
            &                                1, jpj, 1, 1.e-1, numout )
      ENDIF
ENDIF

      CALL wrk_dealloc ( jpi,jpj,zux, zvx , zuy , zvy     )
      CALL wrk_dealloc ( jpi,jpj,zue1, zue2 , zve1 , zve2 )


END SUBROUTINE ldf_tra_smag
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ldf_tra_smag( kt )       ! Empty routine
      WRITE(*,*) 'ldf_dyn_smag: You should not have seen this print! error? check keys ldf:c3d+smag', kt
   END SUBROUTINE ldf_tra_smag
#endif

END MODULE ldftra_smag
