MODULE ldfdyn_smag
   !!======================================================================
   !!                     ***  MODULE  ldftrasmag  ***
   !! Ocean physics:  variable eddy induced velocity coefficients
   !!======================================================================
#if   defined key_dynldf_smag   &&   defined key_dynldf_c3d
   !!----------------------------------------------------------------------
   !!   'key_dynldf_smag'      and           smagorinsky  diffusivity
   !!   'key_dynldf_c3d'                    3D tracer lateral  mixing coef.
   !!----------------------------------------------------------------------
   !!   ldf_eiv      : compute the eddy induced velocity coefficients
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcrnf          ! river runoffs
   USE ldfdyn_oce      ! ocean tracer   lateral physics
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
   PUBLIC ldf_dyn_smag               ! routine called by step.F90
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
   !!                        ***  ldfdyn_smag.F90  ***
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Id$
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_dynldf_smag'             3D lateral eddy viscosity coefficients
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_dyn_smag( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_smag  ***
      !!                   
      !! ** Purpose :   initializations of the horizontal ocean physics
      !!
      !! ** Method  :   3D eddy viscosity coef. 
      !!    M.Griffies, R.Hallberg AMS, 2000
      !! for laplacian:
      !!   Asmag=(C/pi)^2*dx*dy sqrt(D^2), C=3-4
      !! for bilaplacian:
      !!   Bsmag=Asmag*dx*dy/8
      !!   D^2=(du/dx-dv/dy)^2+(dv/dx+du/dy)^2 for Cartesian coordinates
      !!  in general case du/dx ==> e2 d(u/e2)/dx;  du/dy ==> e1 d(u/e1)/dy; 
      !!                  dv/dx ==> e2 d(v/e2)/dx;  dv/dy ==> e1 d(v/e1)/dy
      !!
      !!       laplacian operator   : ahm1, ahm2 defined at T- and F-points
      !!                              ahm3, ahm4 never used
      !!       bilaplacian operator : ahm1, ahm2 never used
      !!                           :  ahm3, ahm4 defined at U- and V-points
      !!       explanation of the default is missingi
      !!  last modified : Maria Luneva, September 2011
      !!----------------------------------------------------------------------
      !! * Modules used
      !! ahm0 here is a background viscosity

      !! * Arguments

      !! * local variables

      INTEGER              :: kt                   ! timestep

      INTEGER  ::   ji, jj, jk                     ! dummy loop indices
      REAL (wp):: zdeltat,zdeltaf,zdeltau,zdeltav  ! temporary scalars
      REAL (wp), POINTER, DIMENSION (:,:) ::   zux, zuy , zvx ,zvy, zue1, zue2, zve1, zve2 
      REAL (wp)::  zcmsmag_1, zcmsmag_2 , zcmsh


      !!----------------------------------------------------------------------

      CALL wrk_alloc( jpi,jpj,zux,zuy,zvx,zvy )
      CALL wrk_alloc( jpi,jpj,zue1,zue2,zve1,zve2 )


      IF(  kt == nit000 ) THEN


      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ldf_dyn_smag : 3D lateral eddy viscosity coefficient'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
     
      ENDIF
      
      zcmsmag_1 = rn_cmsmag_1
      zcmsmag_2 = rn_cmsmag_2
      zcmsh     = rn_cmsh

   


      ! Set ahm1 and ahm2  ( T- and F- points) (used for laplacian operators
      ! =================                       whatever its orientation is)
      IF( ln_dynldf_lap ) THEN
         ! define ahm1 and ahm2 at the right grid point position
         ! (USER: modify ahm1 and ahm2 following your desiderata)
         
         DO jk=1,jpk
           zue2(:,:)=un(:,:,jk)/e2u(:,:)
           zve1(:,:)=vn(:,:,jk)/e1v(:,:)
           zue1(:,:)=un(:,:,jk)/e1u(:,:)
           zve2(:,:)=vn(:,:,jk)/e2v(:,:)

 
           DO jj=2,jpj
            DO ji=2,jpi
            zux(ji,jj)=(zue2(ji,jj)-zue2(ji-1,jj))/e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,jk) * zcmsh
            zvy(ji,jj)=(zve1(ji,jj)-zve1(ji,jj-1))/e2t(ji,jj)*e1t(ji,jj)*tmask(ji,jj,jk) * zcmsh
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

            zdeltat=2._wp /(e1t(ji,jj)**(-2)+e2t(ji,jj)**(-2))
            zdeltaf=2._wp /(e1f(ji,jj)**(-2)+e2f(ji,jj)**(-2))   
            ahm1(ji,jj,jk)=(zcmsmag_1/rpi)**2*zdeltat*                                         &
                            sqrt( (zux(ji,jj)-zvy(ji,jj))**2+                                    &
                     0.0625_wp*(zuy(ji,jj)+zuy(ji,jj-1)+zuy(ji-1,jj)+zuy(ji-1,jj-1)+             &
                            zvx(ji,jj)+zvx(ji,jj-1)+zvx(ji-1,jj)+zvx(ji-1,jj-1))**2)

            ahm2(ji,jj,jk)=(zcmsmag_1/rpi)**2*zdeltaf*                                         &
                            sqrt( (zuy(ji,jj)+zvx(ji,jj))**2+                                    &
                     0.0625_wp*(zux(ji,jj)+zux(ji,jj+1)+zux(ji+1,jj)+zux(ji+1,jj+1)-             &
                             zvy(ji,jj)-zvy(ji,jj+1)-zvy(ji+1,jj)-zvy(ji+1,jj+1))**2)

            ahm1(ji,jj,jk)=MAX(ahm1(ji,jj,jk),rn_ahm_0_lap)
            ahm2(ji,jj,jk)=MAX(ahm2(ji,jj,jk),rn_ahm_0_lap)

! stability criteria or upper limit set from namelist
            ahm1(ji,jj,jk)=MIN(ahm1(ji,jj,jk),zdeltat / (16_wp*rdt),rn_ahm_m_lap)
            ahm2(ji,jj,jk)=MIN(ahm2(ji,jj,jk),zdeltaf / (16_wp*rdt),rn_ahm_m_lap)

            ENDDO
           ENDDO

         ENDDO ! jpk
            ahm1(:,:,jpk) = ahm1(:,:,jpkm1)
            ahm2(:,:,jpk) = ahm2(:,:,jpkm1)
            IF(lwp.and.kt==nit000) WRITE(numout,'(36x," ahm ", 7x)')
            DO jk = 1, jpk

               IF(lwp.and.kt==nit000) WRITE(numout,'(30x,E10.2,8x,i3)') ahm1(jpi/2,jpj/2,jk), jk
            END DO
      CALL lbc_lnk( ahm1, 'T', 1. )   ! Lateral boundary conditions on ( ahtt )
      CALL lbc_lnk( ahm2, 'F', 1. )   ! Lateral boundary conditions on ( ahtt )

      ENDIF    ! ln_dynldf_lap
      


      ! ahm3 and ahm4 at U- and V-points (used for bilaplacian operator
      ! ================================  whatever its orientation is)
      ! Here: ahm is proportional to the cube of the maximum of the grid spacing
      !       in the to horizontal direction

      IF( ln_dynldf_bilap ) THEN
         DO jk=1,jpk
           zue2(:,:) = un(:,:,jk)/e2u(:,:)
           zve1(:,:) = vn(:,:,jk)/e1v(:,:)
           zue1(:,:) = un(:,:,jk)/e1u(:,:)
           zve2(:,:) = vn(:,:,jk)/e2v(:,:)


           DO jj=2,jpj
            DO ji=2,jpi
            zux(ji,jj) = (zue2(ji,jj)-zue2(ji-1,jj))/e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,jk)
            zvy(ji,jj) = (zve1(ji,jj)-zve1(ji,jj-1))/e2t(ji,jj)*e1t(ji,jj)*tmask(ji,jj,jk)
            ENDDO
           ENDDO

           DO jj=1,jpjm1
            DO ji=1,jpim1
            zuy(ji,jj) = (zue1(ji,jj+1)-zue1(ji,jj))/e2f(ji,jj)*e1f(ji,jj)*fmask(ji,jj,jk)
            zvx(ji,jj) = (zve2(ji+1,jj)-zve2(ji,jj))/e1f(ji,jj)*e2f(ji,jj)*fmask(ji,jj,jk)
            ENDDO
           ENDDO


          DO jj=2,jpjm1
           DO ji=2,jpim1
            zdeltau = 2._wp/(e1u(ji,jj)**(-2)+e2u(ji,jj)**(-2))
            zdeltav = 2._wp/(e1v(ji,jj)**(-2)+e2v(ji,jj)**(-2))

            ahm3(ji,jj,jk) = -(zcmsmag_2/rpi)**2/8.0_wp*zdeltau**2*                          &

                         sqrt(0.25_wp*(zux(ji,jj)+zux(ji+1,jj)-zvy(ji,jj)-zvy(ji+1,jj))**2+    &
                              0.25_wp*(zuy(ji,jj)+zuy(ji,jj-1)+zvx(ji,jj)+zvx(ji,jj-1))**2)

            ahm4(ji,jj,jk) = -(zcmsmag_2/rpi)**2/8.0_wp*zdeltav**2*                           &

                         sqrt(0.25_wp*(zux(ji,jj)+zux(ji,jj+1)-zvy(ji,jj)-zvy(ji,jj+1))**2+    &
                              0.25_wp*(zuy(ji,jj)+zuy(ji-1,jj)+zvx(ji-1,jj)+zvx(ji,jj))**2)

            ahm3(ji,jj,jk) = MIN (rn_ahm_0_blp , ahm3(ji,jj,jk) )
            ahm4(ji,jj,jk) = MIN (rn_ahm_0_blp , ahm4(ji,jj,jk) ) 

! stability criteria or upper limit set in namelist

            ahm3(ji,jj,jk) = MAX( ahm3(ji,jj,jk),-zdeltau**2/( 128._wp*rdt ),rn_ahm_m_blp )
            ahm4(ji,jj,jk) = MAX( ahm4(ji,jj,jk),-zdeltav**2/( 128._wp*rdt ),rn_ahm_m_blp )
            

            ENDDO
           ENDDO

         ENDDO
            ahm3(:,:,jpk) = ahm3(:,:,jpkm1)
            ahm4(:,:,jpk) = ahm4(:,:,jpkm1)

      DO jk = 1, jpk
      IF(  kt == nit000 ) THEN

               IF(lwp) WRITE(numout,'(30x,E10.2,8x,i3)') ahm3(jpi/2,jpj/2,jk), jk
      ENDIF   
      END DO
      CALL lbc_lnk( ahm3, 'U', 1. )   ! Lateral boundary conditions
      CALL lbc_lnk( ahm4, 'V', 1. )
   ENDIF

      CALL wrk_dealloc( jpi,jpj,zux,zuy,zvx,zvy )
      CALL wrk_dealloc( jpi,jpj,zue1,zue2,zve1,zve2 )
   !    zumax = MAXVAL( ABS( ahm3(:,:,:) ) )                ! slower than the following loop on NEC SX5
      zdeltat = 0._wp
   If(ln_dynldf_lap)THEN
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdeltat = MAX(zdeltat,ABS(ahm1(ji,jj,jk)),ABS(ahm2(ji,jj,jk)) )
          END DO
        END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( zdeltat )                 ! max over the global domain
      !
      IF( MOD( kt, nwrite ) == 1 .AND. lwp )   WRITE(numout,*) ' ==>> time-step= ',kt,'dynlap:  abs(ahm) max: ', zdeltat
    ENDIF
    If(ln_dynldf_bilap)THEN
     zdeltat = 0._wp
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdeltat = MAX(zdeltat,ABS(ahm3(ji,jj,jk)),ABS(ahm3(ji,jj,jk)) )
          END DO
        END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( zdeltat )                 ! max over the global domain
      !
      IF( MOD( kt, nwrite ) == 1 .AND. lwp )   WRITE(numout,*) ' ==>> time-step= ',kt,'dyn_bilap abs(ahm) max: ', zdeltat
      !
   ENDIF
      !

END SUBROUTINE ldf_dyn_smag
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ldf_dyn_smag( kt )       ! Empty routine
      WRITE(*,*) 'ldf_dyn_smag: You should not have seen this print! error? check keys ldf:c3d+smag', kt
   END SUBROUTINE ldf_dyn_smag
#endif

   END MODULE ldfdyn_smag

