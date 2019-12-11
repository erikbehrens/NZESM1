MODULE p2zexp
   !!======================================================================
   !!                         ***  MODULE p2zsed  ***
   !! TOP :   LOBSTER Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :    -   !  1999    (O. Aumont, C. Le Quere)  original code
   !!              -   !  2001-05 (O. Aumont, E. Kestenare) add sediment computations
   !!             1.0  !  2005-06 (A.-S. Kremeur) new temporal integration for sedpoc
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!             3.5  !  2012-03  (C. Ethe)  Merge PISCES-LOBSTER
   !!----------------------------------------------------------------------
#if defined key_pisces_reduced
   !!----------------------------------------------------------------------
   !!   'key_pisces_reduced'                                     LOBSTER bio-model
   !!----------------------------------------------------------------------
   !!   p2z_exp        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc
   USE sms_pisces
   USE p2zsed
   USE lbclnk
   USE prtctl_trc      ! Print control for debbuging
   USE trd_oce
   USE trdtrc
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_exp    
   PUBLIC   p2z_exp_init 
   PUBLIC   p2z_exp_alloc

   !
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   dminl     !: fraction of sinking POC released in sediments
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   dmin3     !: fraction of sinking POC released at each level
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   sedpocb   !: mass of POC in sediments
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   sedpocn   !: mass of POC in sediments
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   cmask     !: Coastal mask area
   REAL(wp)                                ::   areacot   !: surface coastal area

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p2z_exp( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_exp  ***
      !!
      !! ** Purpose :   MODELS EXPORT OF BIOGENIC MATTER (POC ''SOFT
      !!              TISSUE'') AND ITS DISTRIBUTION IN WATER COLUMN
      !!
      !! ** Method  : - IN THE SURFACE LAYER POC IS PRODUCED ACCORDING TO
      !!              NURTRIENTS AVAILABLE AND GROWTH CONDITIONS. NUTRIENT UPTAKE
      !!              KINETICS FOLLOW MICHAELIS-MENTON FORMULATION. 
      !!              THE TOTAL PARTICLE AMOUNT PRODUCED, IS DISTRIBUTED IN THE WATER
      !!              COLUMN BELOW THE SURFACE LAYER.
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER  ::   ji, jj, jk, jl, ikt
      REAL(wp) ::   zgeolpoc, zfact, zwork, ze3t, zsedpocd, zmaskt
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrbio
      REAL(wp), POINTER, DIMENSION(:,:)   ::  zsedpoca
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p2z_exp')
      !
      IF( kt == nittrc000 )   CALL p2z_exp_init

      CALL wrk_alloc( jpi, jpj, zsedpoca ) 
      zsedpoca(:,:) = 0.

      IF( l_trdtrc )  THEN
         CALL wrk_alloc( jpi, jpj, jpk, ztrbio )   ! temporary save of trends
         ztrbio(:,:,:) = tra(:,:,:,jpno3)
      ENDIF

      ! VERTICAL DISTRIBUTION OF NEWLY PRODUCED BIOGENIC
      ! POC IN THE WATER COLUMN
      ! (PARTS OF NEWLY FORMED MATTER REMAINING IN THE DIFFERENT
      ! LAYERS IS DETERMINED BY DMIN3 DEFINED IN sms_p2z.F90
      ! ----------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ze3t = 1. / fse3t(ji,jj,jk)
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + ze3t * dmin3(ji,jj,jk) * xksi(ji,jj)
            END DO
         END DO
      END DO

      ! Find the last level of the water column
      ! Compute fluxes due to sinking particles (slow)
   

      zgeolpoc = 0.e0         !     Initialization
      ! Release of nutrients from the "simple" sediment
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            ikt = mbkt(ji,jj) 
            tra(ji,jj,ikt,jpno3) = tra(ji,jj,ikt,jpno3) + sedlam * sedpocn(ji,jj) / fse3t(ji,jj,ikt) 
            ! Deposition of organic matter in the sediment
            zwork = vsed * trn(ji,jj,ikt,jpdet)
            zsedpoca(ji,jj) = ( zwork + dminl(ji,jj) * xksi(ji,jj)   &
               &           - sedlam * sedpocn(ji,jj) - sedlostpoc * sedpocn(ji,jj) ) * rdt
            zgeolpoc = zgeolpoc + sedlostpoc * sedpocn(ji,jj) * e1e2t(ji,jj)
         END DO
      END DO

      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            tra(ji,jj,1,jpno3) = tra(ji,jj,1,jpno3) + zgeolpoc * cmask(ji,jj) / areacot / fse3t(ji,jj,1)
         END DO
      END DO

      CALL lbc_lnk( sedpocn, 'T', 1. )
 
      ! Oa & Ek: diagnostics depending on jpdia2d !          left as example
      IF( lk_iomput ) THEN  
         CALL iom_put( "SEDPOC" , sedpocn )
      ELSE
         IF( ln_diatrc )           trc2d(:,:,jp_pcs0_2d + 18) = sedpocn(:,:)
      ENDIF

      
      ! Time filter and swap of arrays
      ! ------------------------------
      IF( neuler == 0 .AND. kt == nittrc000 ) THEN        ! Euler time-stepping at first time-step
        !                                             ! (only swap)
        sedpocn(:,:) = zsedpoca(:,:)
        !                                              
      ELSE
        !
        DO jj = 1, jpj
           DO ji = 1, jpi
              zsedpocd = zsedpoca(ji,jj) - 2. * sedpocn(ji,jj) + sedpocb(ji,jj)      ! time laplacian on tracers
              sedpocb(ji,jj) = sedpocn(ji,jj) + atfp * zsedpocd                     ! sedpocb <-- filtered sedpocn
              sedpocn(ji,jj) = zsedpoca(ji,jj)                                       ! sedpocn <-- sedpoca
           END DO
        END DO
        ! 
      ENDIF
      !
      IF( lrst_trc ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'p2z_exp : POC in sediment fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrtw, 'SEDB'//ctrcnm(jpdet), sedpocb(:,:) )
         CALL iom_rstput( kt, nitrst, numrtw, 'SEDN'//ctrcnm(jpdet), sedpocn(:,:) )
      ENDIF
      !
      IF( l_trdtrc ) THEN
         ztrbio(:,:,:) = tra(:,:,:,jpno3) - ztrbio(:,:,:)
         jl = jp_pcs0_trd + 16
         CALL trd_trc( ztrbio, jl, kt )   ! handle the trend
         CALL wrk_dealloc( jpi, jpj, jpk, ztrbio )   ! temporary save of trends
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zsedpoca)   ! temporary save of trends

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('exp')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p2z_exp')
      !
   END SUBROUTINE p2z_exp

   SUBROUTINE p2z_exp_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE p4z_exp_init  ***
      !! ** purpose :   specific initialisation for export
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zmaskt, zfluo, zfluu
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zrro
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdm0
      !!---------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' p2z_exp: LOBSTER export'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '  compute remineralisation-damping arrays for tracers'
      ENDIF
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      zrro )
      CALL wrk_alloc( jpi, jpj, jpk, zdm0 )


      ! Calculate vertical distribution of newly formed biogenic poc
      ! in the water column in the case of max. possible bottom depth
      ! ------------------------------------------------------------
      zdm0 = 0._wp
      zrro = 1._wp
      DO jk = jpkb, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zfluo = ( fsdepw(ji,jj,jk  ) / fsdepw(ji,jj,jpkb) )**xhr
               zfluu = ( fsdepw(ji,jj,jk+1) / fsdepw(ji,jj,jpkb) )**xhr
               IF( zfluo.GT.1. )   zfluo = 1._wp
               zdm0(ji,jj,jk) = zfluo - zfluu
               IF( jk <= jpkb-1 )   zdm0(ji,jj,jk) = 0._wp
               zrro(ji,jj) = zrro(ji,jj) - zdm0(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      zdm0(:,:,jpk) = zrro(:,:)

      ! Calculate vertical distribution of newly formed biogenic poc
      ! in the water column with realistic topography (first "dry" layer
      ! contains total fraction, which has passed to the upper layers)
      ! ----------------------------------------------------------------------
      dminl(:,:)   = 0._wp
      dmin3(:,:,:) = zdm0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) == 0._wp ) THEN
                  dminl(ji,jj) = dminl(ji,jj) + dmin3(ji,jj,jk)
                  dmin3(ji,jj,jk) = 0._wp
               ENDIF
            END DO
         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( tmask(ji,jj,1) == 0 )   dmin3(ji,jj,1) = 0._wp
         END DO
      END DO

      ! Coastal mask 
      cmask(:,:) = 0._wp
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            IF( tmask(ji,jj,1) /= 0. ) THEN
               zmaskt = tmask(ji+1,jj,1) * tmask(ji-1,jj,1) * tmask(ji,jj+1,1) * tmask(ji,jj-1,1) 
               IF( zmaskt == 0. )   cmask(ji,jj) = 1._wp
            END IF
         END DO
      END DO
      CALL lbc_lnk( cmask , 'T', 1. )      ! lateral boundary conditions on cmask   (sign unchanged)
      areacot = glob_sum( e1e2t(:,:) * cmask(:,:) )
      !
      IF( ln_rsttr ) THEN
         CALL iom_get( numrtr, jpdom_autoglo, 'SEDB'//ctrcnm(jpdet), sedpocb(:,:) )
         CALL iom_get( numrtr, jpdom_autoglo, 'SEDN'//ctrcnm(jpdet), sedpocn(:,:) )
      ELSE
         sedpocb(:,:) = 0._wp
         sedpocn(:,:) = 0._wp
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      zrro )
      CALL wrk_dealloc( jpi, jpj, jpk, zdm0 )
      !
   END SUBROUTINE p2z_exp_init

   INTEGER FUNCTION p2z_exp_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_exp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( cmask(jpi,jpj) , dminl(jpi,jpj) , dmin3(jpi,jpj,jpk), &
         &      sedpocb(jpi,jpj) , sedpocn(jpi,jpj),   STAT=p2z_exp_alloc )
      IF( p2z_exp_alloc /= 0 ) CALL ctl_warn('p2z_exp_alloc : failed to allocate arrays.')
      !
   END FUNCTION p2z_exp_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_exp( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_exp: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_exp
#endif 

   !!======================================================================
END MODULE  p2zexp
