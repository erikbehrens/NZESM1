MODULE traqsr
   !!======================================================================
   !!                       ***  MODULE  traqsr  ***
   !! Ocean physics: solar radiation penetration in the top ocean levels
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1996-01  (G. Madec)  s-coordinates
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!             -   !  2005-11  (G. Madec) zco, zps, sco coordinate
   !!            3.2  !  2009-04  (G. Madec & NEMO team) 
   !!            3.4  !  2012-05  (C. Rousset) store attenuation coef for use in ice model 
   !!            3.6  !  2015-12  (O. Aumont, J. Jouanno, C. Ethe) use vertical profile of chlorophyll
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_qsr      : trend due to the solar radiation penetration
   !!   tra_qsr_init : solar radiation penetration initialization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE trc_oce         ! share SMS/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE prtctl          ! Print control
   USE iom             ! I/O manager
   USE fldread         ! read input fields
   USE restart         ! ocean restart
   USE lib_mpp         ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_qsr       ! routine called by step.F90 (ln_traqsr=T)
   PUBLIC   tra_qsr_init  ! routine called by nemogcm.F90

   !                                 !!* Namelist namtra_qsr: penetrative solar radiation
   LOGICAL , PUBLIC ::   ln_traqsr    !: light absorption (qsr) flag
   LOGICAL , PUBLIC ::   ln_qsr_rgb   !: Red-Green-Blue light absorption flag  
   LOGICAL , PUBLIC ::   ln_qsr_2bd   !: 2 band         light absorption flag
   LOGICAL , PUBLIC ::   ln_qsr_bio   !: bio-model      light absorption flag
   LOGICAL , PUBLIC ::   ln_qsr_ice   !: light penetration for ice-model LIM3 (clem)
   INTEGER , PUBLIC ::   nn_chldta    !: use Chlorophyll data (=1) or not (=0)
   REAL(wp), PUBLIC ::   rn_abs       !: fraction absorbed in the very near surface (RGB & 2 bands)
   REAL(wp), PUBLIC ::   rn_si0       !: very near surface depth of extinction      (RGB & 2 bands)
   REAL(wp), PUBLIC ::   rn_si1       !: deepest depth of extinction (water type I)       (2 bands)
 
   ! Module variables
   REAL(wp) ::   xsi0r                           !: inverse of rn_si0
   REAL(wp) ::   xsi1r                           !: inverse of rn_si1
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_chl   ! structure of input Chl (file informations, fields read)
   INTEGER, PUBLIC ::   nksr              ! levels below which the light cannot penetrate ( depth larger than 391 m)
   REAL(wp), DIMENSION(3,61) ::   rkrgb   !: tabulated attenuation coefficients for RGB absorption

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_qsr( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr  ***
      !!
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!      penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 3 wavebands (RGB) and a ratio rn_abs
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration 
      !!         is given by : zta = 1/e3t dk[ I ] / (rau0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         In z-coordinate case, the computation is only done down to the
      !!      level where I(k) < 1.e-15 W/m2. In addition, the coefficients 
      !!      used for the computation are calculated one for once as they
      !!      depends on k only.
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - save the trend in ttrd ('key_trdtra')
      !!
      !! Reference  : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!              Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!              Morel, A. et Berthon, JF, 1989, Limnol Oceanogr 34(8), 1545-1562
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      INTEGER  ::   irgb                 ! local integers
      REAL(wp) ::   zchl, zcoef, zfact   ! local scalars
      REAL(wp) ::   zc0, zc1, zc2, zc3   !    -         -
      REAL(wp) ::   zz0, zz1, z1_e3t     !    -         -
      REAL(wp) ::   zCb, zCmax, zze, zpsi, zpsimax, zdelpsi, zCtot, zCze
      REAL(wp) ::   zlogc, zlogc2, zlogc3 
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zekb, zekg, zekr
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ze0, ze1, ze2, ze3, zea, ztrdt, zchl3d
      !!--------------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_qsr')
      !
      CALL wrk_alloc( jpi, jpj,      zekb, zekg, zekr        ) 
      CALL wrk_alloc( jpi, jpj, jpk, ze0, ze1, ze2, ze3, zea, zchl3d ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_qsr : penetration of the surface solar radiation'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         IF( .NOT.ln_traqsr )   RETURN
      ENDIF

      IF( l_trdtra ) THEN      ! Save ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt ) 
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
      ENDIF

      !                                        Set before qsr tracer content field
      !                                        ***********************************
      IF( kt == nit000 ) THEN                     ! Set the forcing field at nit000 - 1
         !                                        ! -----------------------------------
         qsr_hc(:,:,:) = 0.e0
         !
         IF( ln_rstart .AND.    &                    ! Restart: read in restart file
              & iom_varid( numror, 'qsr_hc_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nit000-1 qsr tracer content forcing field red in the restart file'
            zfact = 0.5e0
            CALL iom_get( numror, jpdom_autoglo, 'qsr_hc_b', qsr_hc_b )   ! before heat content trend due to Qsr flux
         ELSE                                           ! No restart or restart not found: Euler forward time stepping
            zfact = 1.e0
            qsr_hc_b(:,:,:) = 0.e0
         ENDIF
      ELSE                                        ! Swap of forcing field
         !                                        ! ---------------------
         zfact = 0.5e0
         qsr_hc_b(:,:,:) = qsr_hc(:,:,:)
      ENDIF
      !                                        Compute now qsr tracer content field
      !                                        ************************************
      
      !                                           ! ============================================== !
      IF( lk_qsr_bio .AND. ln_qsr_bio ) THEN      !  bio-model fluxes  : all vertical coordinates  !
         !                                        ! ============================================== !
         DO jk = 1, jpkm1
            qsr_hc(:,:,jk) = r1_rau0_rcp * ( etot3(:,:,jk) - etot3(:,:,jk+1) )
         END DO
         !                                        Add to the general trend
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1 
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z1_e3t = zfact / fse3t(ji,jj,jk)
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + ( qsr_hc_b(ji,jj,jk) + qsr_hc(ji,jj,jk) ) * z1_e3t
               END DO
            END DO
         END DO
         CALL iom_put( 'qsr3d', etot3 )   ! Shortwave Radiation 3D distribution
         ! clem: store attenuation coefficient of the first ocean level
         IF ( ln_qsr_ice ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF ( qsr(ji,jj) /= 0._wp ) THEN
                     fraqsr_1lev(ji,jj) = ( qsr_hc(ji,jj,1) / ( r1_rau0_rcp * qsr(ji,jj) ) )
                  ELSE
                     fraqsr_1lev(ji,jj) = 1.
                  ENDIF
               END DO
            END DO
         ENDIF
         !                                        ! ============================================== !
      ELSE                                        !  Ocean alone : 
         !                                        ! ============================================== !
         !
         !                                                ! ------------------------- !
         IF( ln_qsr_rgb) THEN                             !  R-G-B  light penetration !
            !                                             ! ------------------------- !
            ! Set chlorophyl concentration
            IF( nn_chldta == 1 .OR. nn_chldta == 2 .OR. lk_vvl ) THEN    !*  Variable Chlorophyll or ocean volume
               !
               IF( nn_chldta == 1 ) THEN        !*  2D Variable Chlorophyll
                  !
                  CALL fld_read( kt, 1, sf_chl )            ! Read Chl data and provides it at the current time step
                  DO jk = 1, nksr + 1
                     zchl3d(:,:,jk) = sf_chl(1)%fnow(:,:,1) 
                  ENDDO
                  !
               ELSE IF( nn_chldta == 2 ) THEN    !*   -3-D Variable Chlorophyll
                  !
                  CALL fld_read( kt, 1, sf_chl )            ! Read Chl data and provides it at the current time step
!CDIR NOVERRCHK   !
                  DO jj = 1, jpj
!CDIR NOVERRCHK
                     DO ji = 1, jpi
                        zchl    = sf_chl(1)%fnow(ji,jj,1)
                        zCtot   = 40.6  * zchl**0.459
                        zze     = 568.2 * zCtot**(-0.746)
                        IF( zze > 102. ) zze = 200.0 * zCtot**(-0.293)
                        zlogc   = LOG( zchl )
                        zlogc2  = zlogc * zlogc
                        zlogc3  = zlogc * zlogc * zlogc
                        zCb     = 0.768 + 0.087 * zlogc - 0.179 * zlogc2 - 0.025 * zlogc3
                        zCmax   = 0.299 - 0.289 * zlogc + 0.579 * zlogc2
                        zpsimax = 0.6   - 0.640 * zlogc + 0.021 * zlogc2 + 0.115 * zlogc3
                        zdelpsi = 0.710 + 0.159 * zlogc + 0.021 * zlogc2
                        zCze    = 1.12  * (zchl)**0.803 
                        DO jk = 1, nksr + 1
                           zpsi = fsdept(ji,jj,jk) / zze
                           zchl3d(ji,jj,jk) = zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) / zdelpsi )**2 ) )
                        END DO
                        !
                      END DO
                   END DO
                     !
               ELSE                              !* Variable ocean volume but constant chrlorophyll
                  DO jk = 1, nksr + 1
                     zchl3d(:,:,jk) = 0.05 
                  ENDDO
               ENDIF
               !
               zcoef  = ( 1. - rn_abs ) / 3.e0                        !  equi-partition in R-G-B
               ze0(:,:,1) = rn_abs  * qsr(:,:)
               ze1(:,:,1) = zcoef * qsr(:,:)
               ze2(:,:,1) = zcoef * qsr(:,:)
               ze3(:,:,1) = zcoef * qsr(:,:)
               zea(:,:,1) =         qsr(:,:)
               !
               DO jk = 2, nksr+1
                  !
                  DO jj = 1, jpj                                         ! Separation in R-G-B depending of vertical profile of Chl
!CDIR NOVERRCHK
                     DO ji = 1, jpi
                        zchl = MIN( 10. , MAX( 0.03, zchl3d(ji,jj,jk) ) )
                        irgb = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )
                        zekb(ji,jj) = rkrgb(1,irgb)
                        zekg(ji,jj) = rkrgb(2,irgb)
                        zekr(ji,jj) = rkrgb(3,irgb)
                     END DO
                  END DO
!CDIR NOVERRCHK
                  DO jj = 1, jpj
!CDIR NOVERRCHK   
                     DO ji = 1, jpi
                        zc0 = ze0(ji,jj,jk-1) * EXP( - fse3t(ji,jj,jk-1) * xsi0r     )
                        zc1 = ze1(ji,jj,jk-1) * EXP( - fse3t(ji,jj,jk-1) * zekb(ji,jj) )
                        zc2 = ze2(ji,jj,jk-1) * EXP( - fse3t(ji,jj,jk-1) * zekg(ji,jj) )
                        zc3 = ze3(ji,jj,jk-1) * EXP( - fse3t(ji,jj,jk-1) * zekr(ji,jj) )
                        ze0(ji,jj,jk) = zc0
                        ze1(ji,jj,jk) = zc1
                        ze2(ji,jj,jk) = zc2
                        ze3(ji,jj,jk) = zc3
                        zea(ji,jj,jk) = ( zc0 + zc1 + zc2 + zc3 ) * tmask(ji,jj,jk)
                     END DO
                  END DO
               END DO
               !
               DO jk = 1, nksr                                        ! compute and add qsr trend to ta
                  qsr_hc(:,:,jk) = r1_rau0_rcp * ( zea(:,:,jk) - zea(:,:,jk+1) )
               END DO
               zea(:,:,nksr+1:jpk) = 0.e0     ! below 400m set to zero
               CALL iom_put( 'qsr3d', zea )   ! Shortwave Radiation 3D distribution
               !
               IF ( ln_qsr_ice ) THEN    ! store attenuation coefficient of the first ocean level
!CDIR NOVERRCHK
                  DO jj = 1, jpj                                         ! Separation in R-G-B depending of the surface Chl
!CDIR NOVERRCHK
                     DO ji = 1, jpi
                        zchl = MIN( 10. , MAX( 0.03, zchl3d(ji,jj,1) ) )
                        irgb = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )
                        zekb(ji,jj) = rkrgb(1,irgb)
                        zekg(ji,jj) = rkrgb(2,irgb)
                        zekr(ji,jj) = rkrgb(3,irgb)
                     END DO
                  END DO
                  ! 
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        zc0 = rn_abs * EXP( - fse3t(ji,jj,1) * xsi0r     )
                        zc1 = zcoef  * EXP( - fse3t(ji,jj,1) * zekb(ji,jj) )
                        zc2 = zcoef  * EXP( - fse3t(ji,jj,1) * zekg(ji,jj) )
                        zc3 = zcoef  * EXP( - fse3t(ji,jj,1) * zekr(ji,jj) )
                        fraqsr_1lev(ji,jj) = 1.0 - ( zc0 + zc1 + zc2  + zc3  ) * tmask(ji,jj,2) 
                     END DO
                  END DO
                  !
               ENDIF
               !
            ELSE                                                 !*  Constant Chlorophyll
               DO jk = 1, nksr
                  qsr_hc(:,:,jk) =  etot3(:,:,jk) * qsr(:,:)
               END DO
               !  store attenuation coefficient of the first ocean level
               IF( ln_qsr_ice ) THEN
                  fraqsr_1lev(:,:) = etot3(:,:,1) / r1_rau0_rcp
               ENDIF
           ENDIF

         ENDIF
         !                                                ! ------------------------- !
         IF( ln_qsr_2bd ) THEN                            !  2 band light penetration !
            !                                             ! ------------------------- !
            !
            IF( lk_vvl ) THEN                                  !* variable volume
               zz0   =        rn_abs   * r1_rau0_rcp
               zz1   = ( 1. - rn_abs ) * r1_rau0_rcp
               DO jk = 1, nksr                    ! solar heat absorbed at T-point in the top 400m 
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        zc0 = zz0 * EXP( -fsdepw(ji,jj,jk  )*xsi0r ) + zz1 * EXP( -fsdepw(ji,jj,jk  )*xsi1r )
                        zc1 = zz0 * EXP( -fsdepw(ji,jj,jk+1)*xsi0r ) + zz1 * EXP( -fsdepw(ji,jj,jk+1)*xsi1r )
                        qsr_hc(ji,jj,jk) = qsr(ji,jj) * ( zc0*tmask(ji,jj,jk) - zc1*tmask(ji,jj,jk+1) ) 
                     END DO
                  END DO
               END DO
               ! clem: store attenuation coefficient of the first ocean level
               IF ( ln_qsr_ice ) THEN
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        zc0 = zz0 * EXP( -fsdepw(ji,jj,1)*xsi0r ) + zz1 * EXP( -fsdepw(ji,jj,1)*xsi1r )
                        zc1 = zz0 * EXP( -fsdepw(ji,jj,2)*xsi0r ) + zz1 * EXP( -fsdepw(ji,jj,2)*xsi1r )
                        fraqsr_1lev(ji,jj) = ( zc0*tmask(ji,jj,1) - zc1*tmask(ji,jj,2) ) / r1_rau0_rcp
                     END DO
                  END DO
               ENDIF
            ELSE                                               !* constant volume: coef. computed one for all
               DO jk = 1, nksr
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt.
                        ! (ISF) no light penetration below the ice shelves         
                        qsr_hc(ji,jj,jk) =  etot3(ji,jj,jk) * qsr(ji,jj) * tmask(ji,jj,1)
                     END DO
                  END DO
               END DO
               ! clem: store attenuation coefficient of the first ocean level
               IF ( ln_qsr_ice ) THEN
                  fraqsr_1lev(:,:) = etot3(:,:,1) / r1_rau0_rcp
               ENDIF
               !
            ENDIF
            !
         ENDIF
         !
         !                                        Add to the general trend
         DO jk = 1, nksr
            DO jj = 2, jpjm1 
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z1_e3t = zfact / fse3t(ji,jj,jk)
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + ( qsr_hc_b(ji,jj,jk) + qsr_hc(ji,jj,jk) ) * z1_e3t
               END DO
            END DO
         END DO
         !
      ENDIF
      !
      IF( lrst_oce ) THEN   !                  Write in the ocean restart file
         !                                     *******************************
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'qsr tracer content forcing field written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrow, 'qsr_hc_b'   , qsr_hc      )
         CALL iom_rstput( kt, nitrst, numrow, 'fraqsr_1lev', fraqsr_1lev )   ! default definition in sbcssm 
         !
      ENDIF

      IF( l_trdtra ) THEN     ! qsr tracers trends saved for diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_qsr, ztrdt )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt ) 
      ENDIF
      !                       ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' qsr  - Ta: ', mask1=tmask, clinfo3='tra-ta' )
      !
      CALL wrk_dealloc( jpi, jpj,      zekb, zekg, zekr        ) 
      CALL wrk_dealloc( jpi, jpj, jpk, ze0, ze1, ze2, ze3, zea, zchl3d ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_qsr')
      !
   END SUBROUTINE tra_qsr


   SUBROUTINE tra_qsr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_init  ***
      !!
      !! ** Purpose :   Initialization for the penetrative solar radiation
      !!
      !! ** Method  :   The profile of solar radiation within the ocean is set
      !!      from two length scale of penetration (rn_si0,rn_si1) and a ratio
      !!      (rn_abs). These parameters are read in the namtra_qsr namelist. The
      !!      default values correspond to clear water (type I in Jerlov' 
      !!      (1968) classification.
      !!         called by tra_qsr at the first timestep (nit000)
      !!
      !! ** Action  : - initialize rn_si0, rn_si1 and rn_abs
      !!
      !! Reference : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!----------------------------------------------------------------------
      !
      INTEGER  ::   ji, jj, jk                   ! dummy loop indices
      INTEGER  ::   irgb, ierror, ioptio, nqsr   ! local integer
      INTEGER  ::   ios                          ! Local integer output status for namelist read
      REAL(wp) ::   zz0, zc0  , zc1, zcoef       ! local scalars
      REAL(wp) ::   zz1, zc2  , zc3, zchl        !   -      -
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zekb, zekg, zekr
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ze0, ze1, ze2, ze3, zea
      !
      CHARACTER(len=100) ::   cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::   sn_chl   ! informations about the chlorofyl field to be read
      !!
      NAMELIST/namtra_qsr/  sn_chl, cn_dir, ln_traqsr, ln_qsr_rgb, ln_qsr_2bd, ln_qsr_bio, ln_qsr_ice,  &
         &                  nn_chldta, rn_abs, rn_si0, rn_si1
      !!----------------------------------------------------------------------

      !
      IF( nn_timing == 1 )  CALL timing_start('tra_qsr_init')
      !
      CALL wrk_alloc( jpi, jpj,      zekb, zekg, zekr        ) 
      CALL wrk_alloc( jpi, jpj, jpk, ze0, ze1, ze2, ze3, zea ) 
      !

      REWIND( numnam_ref )              ! Namelist namtra_qsr in reference namelist : Ratio and length of penetration
      READ  ( numnam_ref, namtra_qsr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_qsr in reference namelist', lwp )

      REWIND( numnam_cfg )              !  Namelist namtra_qsr in configuration namelist : Ratio and length of penetration
      READ  ( numnam_cfg, namtra_qsr, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_qsr in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_qsr )
      !
      IF(lwp) THEN                ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_qsr_init : penetration of the surface solar radiation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_qsr : set the parameter of penetration'
         WRITE(numout,*) '      Light penetration (T) or not (F)         ln_traqsr  = ', ln_traqsr
         WRITE(numout,*) '      RGB (Red-Green-Blue) light penetration   ln_qsr_rgb = ', ln_qsr_rgb
         WRITE(numout,*) '      2 band               light penetration   ln_qsr_2bd = ', ln_qsr_2bd
         WRITE(numout,*) '      bio-model            light penetration   ln_qsr_bio = ', ln_qsr_bio
         WRITE(numout,*) '      light penetration for ice-model LIM3     ln_qsr_ice = ', ln_qsr_ice
         WRITE(numout,*) '      RGB : Chl data (=1/2) or cst value (=0)  nn_chldta  = ', nn_chldta
         WRITE(numout,*) '      RGB & 2 bands: fraction of light (rn_si1)    rn_abs = ', rn_abs
         WRITE(numout,*) '      RGB & 2 bands: shortess depth of extinction  rn_si0 = ', rn_si0
         WRITE(numout,*) '      2 bands: longest depth of extinction         rn_si1 = ', rn_si1
      ENDIF

      IF( ln_traqsr ) THEN     ! control consistency
         !                      
         IF( .NOT.lk_qsr_bio .AND. ln_qsr_bio )   THEN
            CALL ctl_warn( 'No bio model : force ln_qsr_bio = FALSE ' )
            ln_qsr_bio = .FALSE.
         ENDIF
         !
         ioptio = 0                      ! Parameter control
         IF( ln_qsr_rgb  )   ioptio = ioptio + 1
         IF( ln_qsr_2bd  )   ioptio = ioptio + 1
         IF( ln_qsr_bio  )   ioptio = ioptio + 1
         !
         IF( ioptio /= 1 ) &
            CALL ctl_stop( '          Choose ONE type of light penetration in namelist namtra_qsr',  &
            &              ' 2 bands, 3 RGB bands or bio-model light penetration' )
         !
         IF( ln_qsr_rgb .AND. nn_chldta == 0 )   nqsr =  1 
         IF( ln_qsr_rgb .AND. nn_chldta == 1 )   nqsr =  2
         IF( ln_qsr_rgb .AND. nn_chldta == 2 )   nqsr =  3
         IF( ln_qsr_2bd                      )   nqsr =  4
         IF( ln_qsr_bio                      )   nqsr =  5
         !
         IF(lwp) THEN                   ! Print the choice
            WRITE(numout,*)
            IF( nqsr ==  1 )   WRITE(numout,*) '         R-G-B   light penetration - Constant Chlorophyll'
            IF( nqsr ==  2 )   WRITE(numout,*) '         R-G-B   light penetration - 2D Chl data '
            IF( nqsr ==  3 )   WRITE(numout,*) '         R-G-B   light penetration - 3D Chl data '
            IF( nqsr ==  4 )   WRITE(numout,*) '         2 bands light penetration'
            IF( nqsr ==  5 )   WRITE(numout,*) '         bio-model light penetration'
         ENDIF
         !
      ENDIF
      !                          ! ===================================== !
      IF( ln_traqsr  ) THEN      !  Initialisation of Light Penetration  !  
         !                       ! ===================================== !
         !
         xsi0r = 1.e0 / rn_si0
         xsi1r = 1.e0 / rn_si1
         !                                ! ---------------------------------- !
         IF( ln_qsr_rgb ) THEN            !  Red-Green-Blue light penetration  !
            !                             ! ---------------------------------- !
            !
            CALL trc_oce_rgb( rkrgb )           !* tabulated attenuation coef.
            !
            !                                   !* level of light extinction
            IF(  ln_sco ) THEN   ;   nksr = jpkm1
            ELSE                 ;   nksr = trc_oce_ext_lev( r_si2, 0.33e2 )
            ENDIF

            IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksr, ' ref depth = ', gdepw_1d(nksr+1), ' m'
            !
            IF( nn_chldta == 1  .OR. nn_chldta == 2 ) THEN           !* Chl data : set sf_chl structure
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '        Chlorophyll read in a file'
               ALLOCATE( sf_chl(1), STAT=ierror )
               IF( ierror > 0 ) THEN
                  CALL ctl_stop( 'tra_qsr_init: unable to allocate sf_chl structure' )   ;   RETURN
               ENDIF
               ALLOCATE( sf_chl(1)%fnow(jpi,jpj,1)   )
               IF( sn_chl%ln_tint )ALLOCATE( sf_chl(1)%fdta(jpi,jpj,1,2) )
               !                                        ! fill sf_chl with sn_chl and control print
               CALL fld_fill( sf_chl, (/ sn_chl /), cn_dir, 'tra_qsr_init',   &
                  &                                         'Solar penetration function of read chlorophyll', 'namtra_qsr' )
               !
            ELSE                                !* constant Chl : compute once for all the distribution of light (etot3)
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '        Constant Chlorophyll concentration = 0.05'
               IF( lk_vvl ) THEN                   ! variable volume
                  IF(lwp) WRITE(numout,*) '        key_vvl: light distribution will be computed at each time step'
               ELSE                                ! constant volume: computes one for all
                  IF(lwp) WRITE(numout,*) '        fixed volume: light distribution computed one for all'
                  !
                  zchl = 0.05                                 ! constant chlorophyll
                  irgb = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )
                  zekb(:,:) = rkrgb(1,irgb)                   ! Separation in R-G-B depending of the chlorophyll 
                  zekg(:,:) = rkrgb(2,irgb)
                  zekr(:,:) = rkrgb(3,irgb)
                  !
                  zcoef = ( 1. - rn_abs ) / 3.e0              ! equi-partition in R-G-B
                  ze0(:,:,1) = rn_abs
                  ze1(:,:,1) = zcoef
                  ze2(:,:,1) = zcoef 
                  ze3(:,:,1) = zcoef
                  zea(:,:,1) = tmask(:,:,1)                   ! = ( ze0+ze1+z2+ze3 ) * tmask
               
                  DO jk = 2, nksr+1
!CDIR NOVERRCHK
                     DO jj = 1, jpj
!CDIR NOVERRCHK   
                        DO ji = 1, jpi
                           zc0 = ze0(ji,jj,jk-1) * EXP( - e3t_0(ji,jj,jk-1) * xsi0r     )
                           zc1 = ze1(ji,jj,jk-1) * EXP( - e3t_0(ji,jj,jk-1) * zekb(ji,jj) )
                           zc2 = ze2(ji,jj,jk-1) * EXP( - e3t_0(ji,jj,jk-1) * zekg(ji,jj) )
                           zc3 = ze3(ji,jj,jk-1) * EXP( - e3t_0(ji,jj,jk-1) * zekr(ji,jj) )
                           ze0(ji,jj,jk) = zc0
                           ze1(ji,jj,jk) = zc1
                           ze2(ji,jj,jk) = zc2
                           ze3(ji,jj,jk) = zc3
                           zea(ji,jj,jk) = ( zc0 + zc1 + zc2 + zc3 ) * tmask(ji,jj,jk)
                        END DO
                     END DO
                  END DO 
                  !
                  DO jk = 1, nksr
                     ! (ISF) no light penetration below the ice shelves
                     etot3(:,:,jk) = r1_rau0_rcp * ( zea(:,:,jk) - zea(:,:,jk+1) ) * tmask(:,:,1)
                  END DO
                  etot3(:,:,nksr+1:jpk) = 0.e0                ! below 400m set to zero
               ENDIF
            ENDIF
            !
         ENDIF
            !                             ! ---------------------------------- !
         IF( ln_qsr_2bd ) THEN            !    2 bands    light penetration    !
            !                             ! ---------------------------------- !
            !
            !                                ! level of light extinction
            nksr = trc_oce_ext_lev( rn_si1, 1.e2 )
            IF(lwp) THEN
               WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksr, ' ref depth = ', gdepw_1d(nksr+1), ' m'
            ENDIF
            !
            IF( lk_vvl ) THEN                   ! variable volume
               IF(lwp) WRITE(numout,*) '        key_vvl: light distribution will be computed at each time step'
            ELSE                                ! constant volume: computes one for all
               zz0 =        rn_abs   * r1_rau0_rcp
               zz1 = ( 1. - rn_abs ) * r1_rau0_rcp
               DO jk = 1, nksr                    !*  solar heat absorbed at T-point computed once for all
                  DO jj = 1, jpj                              ! top 400 meters
                     DO ji = 1, jpi
                        zc0 = zz0 * EXP( -fsdepw(ji,jj,jk  )*xsi0r ) + zz1 * EXP( -fsdepw(ji,jj,jk  )*xsi1r )
                        zc1 = zz0 * EXP( -fsdepw(ji,jj,jk+1)*xsi0r ) + zz1 * EXP( -fsdepw(ji,jj,jk+1)*xsi1r )
                        etot3(ji,jj,jk) = (  zc0 * tmask(ji,jj,jk) - zc1 * tmask(ji,jj,jk+1)  ) * tmask(ji,jj,1) 
                     END DO
                  END DO
               END DO
               etot3(:,:,nksr+1:jpk) = 0.e0                   ! below 400m set to zero
               !
            ENDIF
         ENDIF
         !                       ! ===================================== !
      ELSE                       !        No light penetration           !                   
         !                       ! ===================================== !
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'tra_qsr_init : NO solar flux penetration'
            WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      ! initialisation of fraqsr_1lev used in sbcssm
      IF( iom_varid( numror, 'fraqsr_1lev', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( numror, jpdom_autoglo, 'fraqsr_1lev'  , fraqsr_1lev  )
      ELSE
         fraqsr_1lev(:,:) = 1._wp   ! default definition
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      zekb, zekg, zekr        ) 
      CALL wrk_dealloc( jpi, jpj, jpk, ze0, ze1, ze2, ze3, zea ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_qsr_init')
      !
   END SUBROUTINE tra_qsr_init

   !!======================================================================
END MODULE traqsr
