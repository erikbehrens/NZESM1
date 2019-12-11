MODULE sbcisf
   !!======================================================================
   !!                       ***  MODULE  sbcisf  ***
   !! Surface module :  update surface ocean boundary condition under ice
   !!                   shelf
   !!======================================================================
   !! History :  3.2   !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X   !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4   !  2013-03  (P. Mathiot) Merging + parametrization
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_isf        : update sbc under ice shelf
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE eosbn2          ! equation of state
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE lbclnk          !
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing
   USE lib_fortran     ! glob_sum
   USE zdfbfr
   USE fldread         ! read input field at current time step
   USE lib_fortran, ONLY: glob_sum

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_isf, sbc_isf_div, sbc_isf_alloc  ! routine called in sbcmod and divcur

   ! public in order to be able to output then 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   risf_tsc_b, risf_tsc   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qisf              !: net heat flux from ice shelf
   REAL(wp), PUBLIC ::   rn_hisf_tbl                 !: thickness of top boundary layer [m]
   LOGICAL , PUBLIC ::   ln_divisf                   !: flag to correct divergence 
   INTEGER , PUBLIC ::   nn_isfblk                   !: 
   INTEGER , PUBLIC ::   nn_gammablk                 !:
   LOGICAL , PUBLIC ::   ln_conserve                 !:
   REAL(wp), PUBLIC ::   rn_gammat0                  !: temperature exchange coeficient
   REAL(wp), PUBLIC ::   rn_gammas0                  !: salinity    exchange coeficient 
   REAL(wp), PUBLIC ::   rdivisf                     !: flag to test if fwf apply on divergence

   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  rzisf_tbl              !:depth of calving front (shallowest point) nn_isf ==2/3
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  rhisf_tbl, rhisf_tbl_0 !:thickness of tbl
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  r1_hisf_tbl            !:1/thickness of tbl
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  ralpha                 !:proportion of bottom cell influenced by tbl 
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  risfLeff               !:effective length (Leff) BG03 nn_isf==2
   REAL(wp)   , PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  ttbl, stbl, utbl, vtbl !:top boundary layer variable at T point
   INTEGER,    PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)     ::  misfkt, misfkb         !:Level of ice shelf base


   REAL(wp), PUBLIC, SAVE ::   rcpi   = 2000.0_wp     ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   kappa  = 1.54e-6_wp    ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   rhoisf = 920.0_wp      ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   tsurf  = -20.0_wp      ! phycst ?
   REAL(wp), PUBLIC, SAVE ::   lfusisf= 0.334e6_wp    ! phycst ?

!: Variable used in fldread to read the forcing file (nn_isf == 4 .OR. nn_isf == 3)
   CHARACTER(len=100), PUBLIC ::   cn_dirisf  = './'    !: Root directory for location of ssr files
   TYPE(FLD_N)       , PUBLIC ::   sn_qisf, sn_fwfisf     !: information about the runoff file to be read
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_qisf, sf_fwfisf
   TYPE(FLD_N)       , PUBLIC ::   sn_rnfisf              !: information about the runoff file to be read
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_rnfisf           
   TYPE(FLD_N)       , PUBLIC ::   sn_depmax_isf, sn_depmin_isf, sn_Leff_isf     !: information about the runoff file to be read
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.0 , LOCEAN-IPSL (2008)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
 
  SUBROUTINE sbc_isf(kt)
    INTEGER, INTENT(in)          ::   kt         ! ocean time step
    INTEGER                      ::   ji, jj, jk, ijkmin, inum, ierror
    INTEGER                      ::   ikt, ikb   ! top and bottom level of the isf boundary layer
    REAL(wp)                     ::   zgreenland_fwfisf_sum, zantarctica_fwfisf_sum
    REAL(wp)                     ::   rmin
    REAL(wp)                     ::   zhk
    REAL(wp)                     ::   zt_frz, zpress
    CHARACTER(len=256)           ::   cfisf , cvarzisf, cvarhisf   ! name for isf file
    CHARACTER(LEN=256)           :: cnameis                     ! name of iceshelf file
    CHARACTER (LEN=32)           :: cvarLeff                    ! variable name for efficient Length scale
    INTEGER           ::   ios           ! Local integer output status for namelist read

    REAL(wp), DIMENSION(:,:,:), POINTER :: zfwfisf3d, zqhcisf3d, zqlatisf3d
    REAL(wp), DIMENSION(:,:  ), POINTER :: zqhcisf2d
      !
      !!---------------------------------------------------------------------
      NAMELIST/namsbc_isf/ nn_isfblk, rn_hisf_tbl, ln_divisf, ln_conserve, rn_gammat0, rn_gammas0, nn_gammablk, &
                         & sn_fwfisf, sn_qisf, sn_rnfisf, sn_depmax_isf, sn_depmin_isf, sn_Leff_isf
      !
      !
      !                                         ! ====================== !
      IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !
         REWIND( numnam_ref )              ! Namelist namsbc_rnf in reference namelist : Runoffs 
         READ  ( numnam_ref, namsbc_isf, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_isf in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_rnf in configuration namelist : Runoffs
         READ  ( numnam_cfg, namsbc_isf, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_isf in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_isf )


         IF ( lwp ) WRITE(numout,*)
         IF ( lwp ) WRITE(numout,*) 'sbc_isf: heat flux of the ice shelf'
         IF ( lwp ) WRITE(numout,*) '~~~~~~~~~'
         IF ( lwp ) WRITE(numout,*) 'sbcisf :' 
         IF ( lwp ) WRITE(numout,*) '~~~~~~~~'
         IF ( lwp ) WRITE(numout,*) '        nn_isf      = ', nn_isf
         IF ( lwp ) WRITE(numout,*) '        nn_isfblk   = ', nn_isfblk
         IF ( lwp ) WRITE(numout,*) '        rn_hisf_tbl = ', rn_hisf_tbl
         IF ( lwp ) WRITE(numout,*) '        ln_divisf   = ', ln_divisf 
         IF ( lwp ) WRITE(numout,*) '        nn_gammablk = ', nn_gammablk 
         IF ( lwp ) WRITE(numout,*) '        rn_tfri2    = ', rn_tfri2 
         IF (ln_divisf) THEN       ! keep it in the namelist ??? used true anyway as for runoff ? (PM)
            rdivisf = 1._wp
         ELSE
            rdivisf = 0._wp
         END IF
         !
         ! Allocate public variable
         IF ( sbc_isf_alloc()  /= 0 )         CALL ctl_stop( 'STOP', 'sbc_isf : unable to allocate arrays' )
         !
         ! initialisation
         qisf(:,:)        = 0._wp  ; fwfisf(:,:) = 0._wp
         risf_tsc(:,:,:)  = 0._wp
         !
         ! define isf tbl tickness, top and bottom indice
         IF      (nn_isf == 1) THEN
            rhisf_tbl(:,:) = rn_hisf_tbl
            misfkt(:,:)    = mikt(:,:)         ! same indice for bg03 et cav => used in isfdiv
         ELSE IF ((nn_isf == 3) .OR. (nn_isf == 2)) THEN
            ALLOCATE( sf_rnfisf(1), STAT=ierror )
            ALLOCATE( sf_rnfisf(1)%fnow(jpi,jpj,1), sf_rnfisf(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_rnfisf, (/ sn_rnfisf /), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf' )

            !: read effective lenght (BG03)
            IF (nn_isf == 2) THEN
               ! Read Data and save some integral values
               CALL iom_open( sn_Leff_isf%clname, inum )
               cvarLeff  = 'soLeff'               !: variable name for Efficient Length scale
               CALL iom_get( inum, jpdom_data, cvarLeff, risfLeff , 1)
               CALL iom_close(inum)
               !
               risfLeff = risfLeff*1000           !: convertion in m
            END IF

           ! read depth of the top and bottom of the isf top boundary layer (in this case, isf front depth and grounding line depth)
            CALL iom_open( sn_depmax_isf%clname, inum )
            cvarhisf = TRIM(sn_depmax_isf%clvar)
            CALL iom_get( inum, jpdom_data, cvarhisf, rhisf_tbl, 1) !: depth of deepest point of the ice shelf base
            CALL iom_close(inum)
            !
            CALL iom_open( sn_depmin_isf%clname, inum )
            cvarzisf = TRIM(sn_depmin_isf%clvar)
            CALL iom_get( inum, jpdom_data, cvarzisf, rzisf_tbl, 1) !: depth of shallowest point of the ice shelves base
            CALL iom_close(inum)
            !
            rhisf_tbl(:,:) = rhisf_tbl(:,:) - rzisf_tbl(:,:)        !: tickness isf boundary layer

           !! compute first level of the top boundary layer
           DO ji = 1, jpi
              DO jj = 1, jpj
                  jk = 2
                  DO WHILE ( jk .LE. mbkt(ji,jj) .AND. gdepw_0(ji,jj,jk) < rzisf_tbl(ji,jj) ) ;  jk = jk + 1 ;  END DO
                  misfkt(ji,jj) = jk-1
               END DO
            END DO

         ELSE IF ( nn_isf == 4 ) THEN
            ! as in nn_isf == 1
            rhisf_tbl(:,:) = rn_hisf_tbl
            misfkt(:,:)    = mikt(:,:)         ! same indice for bg03 et cav => used in isfdiv
            
            ! load variable used in fldread (use for temporal interpolation of isf fwf forcing)
            ALLOCATE( sf_fwfisf(1), sf_qisf(1), STAT=ierror )
            ALLOCATE( sf_fwfisf(1)%fnow(jpi,jpj,1), sf_fwfisf(1)%fdta(jpi,jpj,1,2) )
            ALLOCATE( sf_qisf(1)%fnow(jpi,jpj,1), sf_qisf(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_fwfisf, (/ sn_fwfisf /), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf' )
            !CALL fld_fill( sf_qisf  , (/ sn_qisf   /), cn_dirisf, 'sbc_isf_init', 'read heat flux isf data'       , 'namsbc_isf' )
         END IF
         
         ! save initial top boundary layer thickness         
         rhisf_tbl_0(:,:) = rhisf_tbl(:,:)

      END IF

      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         fwfisf_b  (:,:  ) = fwfisf  (:,:  )               ! Swap the ocean forcing fields except at nit000
         risf_tsc_b(:,:,:) = risf_tsc(:,:,:)               ! where before fields are set at the end of the routine
         !
      ENDIF

      IF( MOD( kt-1, nn_fsbc) == 0 ) THEN

         ! compute bottom level of isf tbl and thickness of tbl below the ice shelf
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkt(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               rhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), fse3t_n(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt, mbkt(ji,jj)
                  IF ( (SUM(fse3t_n(ji,jj,ikt:jk-1)) .LT. rhisf_tbl(ji,jj)) .AND. (tmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               rhisf_tbl(ji,jj) = MIN(rhisf_tbl(ji,jj), SUM(fse3t_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.
               misfkb(ji,jj) = ikb                                                  ! last wet level of the tbl
               r1_hisf_tbl(ji,jj) = 1._wp / rhisf_tbl(ji,jj)

               zhk           = SUM( fse3t(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)  ! proportion of tbl cover by cell from ikt to ikb - 1
               ralpha(ji,jj) = rhisf_tbl(ji,jj) * (1._wp - zhk ) / fse3t(ji,jj,ikb)  ! proportion of bottom cell influenced by boundary layer
            END DO
         END DO

         ! compute salf and heat flux
         IF (nn_isf == 1) THEN
            ! realistic ice shelf formulation
            ! compute T/S/U/V for the top boundary layer
            CALL sbc_isf_tbl(tsn(:,:,:,jp_tem),ttbl(:,:),'T')
            CALL sbc_isf_tbl(tsn(:,:,:,jp_sal),stbl(:,:),'T')
            CALL sbc_isf_tbl(un(:,:,:),utbl(:,:),'U')
            CALL sbc_isf_tbl(vn(:,:,:),vtbl(:,:),'V')
            ! iom print
            CALL iom_put('ttbl',ttbl(:,:))
            CALL iom_put('stbl',stbl(:,:))
            CALL iom_put('utbl',utbl(:,:))
            CALL iom_put('vtbl',vtbl(:,:))
            ! compute fwf and heat flux
            CALL sbc_isf_cav (kt)

         ELSE IF (nn_isf == 2) THEN
            ! Beckmann and Goosse parametrisation 
            stbl(:,:)   = soce
            CALL sbc_isf_bg03(kt)

         ELSE IF (nn_isf == 3) THEN
            ! specified runoff in depth (Mathiot et al., XXXX in preparation)
            CALL fld_read ( kt, nn_fsbc, sf_rnfisf   )
            fwfisf(:,:) = - sf_rnfisf(1)%fnow(:,:,1)         ! fresh water flux from the isf (fwfisf <0 mean melting) 

            IF( lk_oasis) THEN
            ! nn_coupled_iceshelf_fluxes uninitialised unless lk_oasis=true
            IF( nn_coupled_iceshelf_fluxes .gt. 0 ) THEN

              ! Adjust total iceshelf melt rates so that sum of iceberg calving and iceshelf melting in the northern
              ! and southern hemispheres equals rate of increase of mass of greenland and antarctic ice sheets
              ! to preserve total freshwater conservation in coupled models without an active ice sheet model.
 
              ! All related global sums must be done bit reproducibly
               zgreenland_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * greenland_icesheet_mask(:,:) )

               ! use ABS function because we need to preserve the sign of fwfisf
               WHERE( greenland_icesheet_mask(:,:) == 1.0 )                                                                  &
              &    fwfisf(:,:) = fwfisf(:,:)  * ABS( greenland_icesheet_mass_rate_of_change * (1.0-rn_greenland_calving_fraction) &
              &                           / ( zgreenland_fwfisf_sum + 1.0e-10_wp ) )

               ! check
               IF(lwp) WRITE(numout, *) 'Greenland iceshelf melting climatology (kg/s) : ',zgreenland_fwfisf_sum

               zgreenland_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * greenland_icesheet_mask(:,:) )

               IF(lwp) WRITE(numout, *) 'Greenland iceshelf melting adjusted value (kg/s) : ',zgreenland_fwfisf_sum

               zantarctica_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * antarctica_icesheet_mask(:,:) )

               ! use ABS function because we need to preserve the sign of fwfisf
               WHERE( antarctica_icesheet_mask(:,:) == 1.0 ) &
              &    fwfisf(:,:) = fwfisf(:,:)  * ABS( antarctica_icesheet_mass_rate_of_change * (1.0-rn_antarctica_calving_fraction) &
              &                           / ( zantarctica_fwfisf_sum + 1.0e-10_wp ) )
      
               ! check
               IF(lwp) WRITE(numout, *) 'Antarctica iceshelf melting climatology (kg/s) : ',zantarctica_fwfisf_sum

               zantarctica_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * antarctica_icesheet_mask(:,:) )

               IF(lwp) WRITE(numout, *) 'Antarctica iceshelf melting adjusted value (kg/s) : ',zantarctica_fwfisf_sum

            ENDIF
            ENDIF

            qisf(:,:)   = fwfisf(:,:) * lfusisf              ! heat        flux
            stbl(:,:)   = soce

         ELSE IF (nn_isf == 4) THEN
            ! specified fwf and heat flux forcing beneath the ice shelf
            CALL fld_read ( kt, nn_fsbc, sf_fwfisf   )
            !CALL fld_read ( kt, nn_fsbc, sf_qisf   )
            fwfisf(:,:) = sf_fwfisf(1)%fnow(:,:,1)            ! fwf

            IF( lk_oasis) THEN
            ! nn_coupled_iceshelf_fluxes uninitialised unless lk_oasis=true
            IF( nn_coupled_iceshelf_fluxes .gt. 0 ) THEN

              ! Adjust total iceshelf melt rates so that sum of iceberg calving and iceshelf melting in the northern
              ! and southern hemispheres equals rate of increase of mass of greenland and antarctic ice sheets
              ! to preserve total freshwater conservation in coupled models without an active ice sheet model.

              ! All related global sums must be done bit reproducibly
               zgreenland_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * greenland_icesheet_mask(:,:) )

               ! use ABS function because we need to preserve the sign of fwfisf
               WHERE( greenland_icesheet_mask(:,:) == 1.0 )                                                                  &
              &    fwfisf(:,:) = fwfisf(:,:)  * ABS( greenland_icesheet_mass_rate_of_change * (1.0-rn_greenland_calving_fraction) &
              &                           / ( zgreenland_fwfisf_sum + 1.0e-10_wp ) )

               ! check
               IF(lwp) WRITE(numout, *) 'Greenland iceshelf melting climatology (kg/s) : ',zgreenland_fwfisf_sum

               zgreenland_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * greenland_icesheet_mask(:,:) )

               IF(lwp) WRITE(numout, *) 'Greenland iceshelf melting adjusted value (kg/s) : ',zgreenland_fwfisf_sum

               zantarctica_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * antarctica_icesheet_mask(:,:) )

               ! use ABS function because we need to preserve the sign of fwfisf
               WHERE( antarctica_icesheet_mask(:,:) == 1.0 ) &
              &    fwfisf(:,:) = fwfisf(:,:)  * ABS( antarctica_icesheet_mass_rate_of_change * (1.0-rn_antarctica_calving_fraction) &
              &                           / ( zantarctica_fwfisf_sum + 1.0e-10_wp ) )
      
               ! check
               IF(lwp) WRITE(numout, *) 'Antarctica iceshelf melting climatology (kg/s) : ',zantarctica_fwfisf_sum

               zantarctica_fwfisf_sum = glob_sum( fwfisf(:,:) * e1t(:,:) * e2t(:,:) * antarctica_icesheet_mask(:,:) )

               IF(lwp) WRITE(numout, *) 'Antarctica iceshelf melting adjusted value (kg/s) : ',zantarctica_fwfisf_sum

            ENDIF
            ENDIF

            qisf(:,:)   = fwfisf(:,:) * lfusisf              ! heat        flux
            !qisf(:,:)   = sf_qisf(1)%fnow(:,:,1)              ! heat flux
            stbl(:,:)   = soce

         END IF
         ! compute tsc due to isf
         ! WARNING water add at temp = 0C, correction term is added, maybe better here but need a 3D variable).
!         zpress = grav*rau0*fsdept(ji,jj,jk)*1.e-04
         zt_frz = -1.9 !eos_fzp( tsn(ji,jj,jk,jp_sal), zpress )
         risf_tsc(:,:,jp_tem) = qisf(:,:) * r1_rau0_rcp - rdivisf * fwfisf(:,:) * zt_frz * r1_rau0 !
         
         ! salt effect already take into account in vertical advection
         risf_tsc(:,:,jp_sal) = (1.0_wp-rdivisf) * fwfisf(:,:) * stbl(:,:) * r1_rau0

         ! output
         IF( iom_use('qlatisf' ) )   CALL iom_put('qlatisf', qisf)
         IF( iom_use('fwfisf'  ) )   CALL iom_put('fwfisf' , fwfisf * stbl(:,:) / soce )

         ! if apply only on the trend and not as a volume flux (rdivisf = 0), fwfisf have to be set to 0 now
         fwfisf(:,:) = rdivisf * fwfisf(:,:)         
 
         ! lbclnk
         CALL lbc_lnk(risf_tsc(:,:,jp_tem),'T',1.)
         CALL lbc_lnk(risf_tsc(:,:,jp_sal),'T',1.)
         CALL lbc_lnk(fwfisf(:,:)   ,'T',1.)
         CALL lbc_lnk(qisf(:,:)     ,'T',1.)

!=============================================================================================================================================
         IF ( iom_use('fwfisf3d') .OR. iom_use('qlatisf3d') .OR. iom_use('qhcisf3d') .OR. iom_use('qhcisf')) THEN
            CALL wrk_alloc( jpi,jpj,jpk, zfwfisf3d, zqhcisf3d, zqlatisf3d )
            CALL wrk_alloc( jpi,jpj,     zqhcisf2d                        )

            zfwfisf3d(:,:,:) = 0.0_wp                         ! 3d ice shelf melting (kg/m2/s)
            zqhcisf3d(:,:,:) = 0.0_wp                         ! 3d heat content flux (W/m2)
            zqlatisf3d(:,:,:)= 0.0_wp                         ! 3d ice shelf melting latent heat flux (W/m2)
            zqhcisf2d(:,:)   = fwfisf(:,:) * zt_frz * rcp     ! 2d heat content flux (W/m2)

            DO jj = 1,jpj
               DO ji = 1,jpi
                  ikt = misfkt(ji,jj)
                  ikb = misfkb(ji,jj)
                  DO jk = ikt, ikb - 1
                     zfwfisf3d (ji,jj,jk) = zfwfisf3d (ji,jj,jk) + fwfisf   (ji,jj) * r1_hisf_tbl(ji,jj) * fse3t(ji,jj,jk)
                     zqhcisf3d (ji,jj,jk) = zqhcisf3d (ji,jj,jk) + zqhcisf2d(ji,jj) * r1_hisf_tbl(ji,jj) * fse3t(ji,jj,jk)
                     zqlatisf3d(ji,jj,jk) = zqlatisf3d(ji,jj,jk) + qisf     (ji,jj) * r1_hisf_tbl(ji,jj) * fse3t(ji,jj,jk)
                  END DO
                  zfwfisf3d (ji,jj,jk) = zfwfisf3d (ji,jj,jk) + fwfisf   (ji,jj) * r1_hisf_tbl(ji,jj) * ralpha(ji,jj) * fse3t(ji,jj,jk)
                  zqhcisf3d (ji,jj,jk) = zqhcisf3d (ji,jj,jk) + zqhcisf2d(ji,jj) * r1_hisf_tbl(ji,jj) * ralpha(ji,jj) * fse3t(ji,jj,jk)
                  zqlatisf3d(ji,jj,jk) = zqlatisf3d(ji,jj,jk) + qisf     (ji,jj) * r1_hisf_tbl(ji,jj) * ralpha(ji,jj) * fse3t(ji,jj,jk)
               END DO
            END DO

            CALL iom_put('fwfisf3d' , zfwfisf3d (:,:,:))
            CALL iom_put('qlatisf3d', zqlatisf3d(:,:,:))
            CALL iom_put('qhcisf3d' , zqhcisf3d (:,:,:))
            CALL iom_put('qhcisf'   , zqhcisf2d (:,:  ))

            CALL wrk_dealloc( jpi,jpj,jpk, zfwfisf3d, zqhcisf3d, zqlatisf3d )
            CALL wrk_dealloc( jpi,jpj,     zqhcisf2d                        )
         END IF
!=============================================================================================================================================

         IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
            IF( ln_rstart .AND.    &                     ! Restart: read in restart file
                 & iom_varid( numror, 'fwf_isf_b', ldstop = .FALSE. ) > 0 ) THEN
               IF(lwp) WRITE(numout,*) '          nit000-1 isf tracer content forcing fields read in the restart file'
               CALL iom_get( numror, jpdom_autoglo, 'fwf_isf_b', fwfisf_b(:,:) )   ! before salt content isf_tsc trend
               CALL iom_get( numror, jpdom_autoglo, 'isf_sc_b', risf_tsc_b(:,:,jp_sal) )   ! before salt content isf_tsc trend
               CALL iom_get( numror, jpdom_autoglo, 'isf_hc_b', risf_tsc_b(:,:,jp_tem) )   ! before salt content isf_tsc trend
            ELSE
               fwfisf_b(:,:)    = fwfisf(:,:)
               risf_tsc_b(:,:,:)= risf_tsc(:,:,:)
            END IF
         ENDIF
         ! 
      END IF
  
  END SUBROUTINE sbc_isf

  INTEGER FUNCTION sbc_isf_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION sbc_isf_rnf_alloc  ***
      !!----------------------------------------------------------------------
      sbc_isf_alloc = 0       ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( qisf ) ) THEN
         ALLOCATE(  risf_tsc(jpi,jpj,jpts), risf_tsc_b(jpi,jpj,jpts), qisf(jpi,jpj)   , &
               &    rhisf_tbl(jpi,jpj)    , r1_hisf_tbl(jpi,jpj), rzisf_tbl(jpi,jpj)  , &
               &    ttbl(jpi,jpj)         , stbl(jpi,jpj)       , utbl(jpi,jpj)       , &
               &    vtbl(jpi, jpj)        , risfLeff(jpi,jpj)   , rhisf_tbl_0(jpi,jpj), &
               &    ralpha(jpi,jpj)       , misfkt(jpi,jpj)     , misfkb(jpi,jpj)     , &
               &    STAT= sbc_isf_alloc )
         !
         IF( lk_mpp                  )   CALL mpp_sum ( sbc_isf_alloc )
         IF( sbc_isf_alloc /= 0 )   CALL ctl_warn('sbc_isf_alloc: failed to allocate arrays.')
         !
      ENDIF
  END FUNCTION

  SUBROUTINE sbc_isf_bg03(kt)
   !!==========================================================================
   !!                 *** SUBROUTINE sbcisf_bg03  ***
   !! add net heat and fresh water flux from ice shelf melting
   !! into the adjacent ocean using the parameterisation by
   !! Beckmann and Goosse (2003), "A parameterization of ice shelf-ocean
   !!     interaction for climate models", Ocean Modelling 5(2003) 157-170.
   !!  (hereafter BG)
   !!==========================================================================
   !!----------------------------------------------------------------------
   !!   sbc_isf_bg03      : routine called from sbcmod
   !!----------------------------------------------------------------------
   !!
   !! ** Purpose   :   Add heat and fresh water fluxes due to ice shelf melting
   !! ** Reference :   Beckmann et Goosse, 2003, Ocean Modelling
   !!
   !! History :
   !!      !  06-02  (C. Wang) Original code
   !!----------------------------------------------------------------------

    INTEGER, INTENT ( in ) :: kt

    INTEGER :: ji, jj, jk, jish  !temporary integer
    INTEGER :: ijkmin
    INTEGER :: ii, ij, ik 
    INTEGER :: inum

    REAL(wp) :: zt_sum      ! sum of the temperature between 200m and 600m
    REAL(wp) :: zt_ave      ! averaged temperature between 200m and 600m
    REAL(wp) :: zt_frz      ! freezing point temperature at depth z
    REAL(wp) :: zpress      ! pressure to compute the freezing point in depth
    
    !!----------------------------------------------------------------------
    IF ( nn_timing == 1 ) CALL timing_start('sbc_isf_bg03')
     !

    ! This test is false only in the very first time step of a run (JMM ???- Initialy build to skip 1rst year of run )
    DO ji = 1, jpi
       DO jj = 1, jpj
          ik = misfkt(ji,jj)
          !! Initialize arrays to 0 (each step)
          zt_sum = 0.e0_wp
          IF ( ik .GT. 1 ) THEN
    ! 3. -----------the average temperature between 200m and 600m ---------------------
             DO jk = misfkt(ji,jj),misfkb(ji,jj)
             ! freezing point temperature  at ice shelf base BG eq. 2 (JMM sign pb ??? +7.64e-4 !!!)
             ! after verif with UNESCO, wrong sign in BG eq. 2
             ! Calculate freezing temperature
                zpress = grav*rau0*fsdept(ji,jj,ik)*1.e-04 
                CALL eos_fzp(tsb(ji,jj,ik,jp_sal), zt_frz, zpress) 
                zt_sum = zt_sum + (tsn(ji,jj,ik,jp_tem)-zt_frz) * fse3t(ji,jj,ik) * tmask(ji,jj,ik)  ! sum temp
             ENDDO
             zt_ave = zt_sum/rhisf_tbl(ji,jj) ! calcul mean value
    
    ! 4. ------------Net heat flux and fresh water flux due to the ice shelf
          ! For those corresponding to zonal boundary    
             qisf(ji,jj) = - rau0 * rcp * rn_gammat0 * risfLeff(ji,jj) * e1t(ji,jj) * zt_ave  &
                         & / (e1t(ji,jj) * e2t(ji,jj)) * tmask(ji,jj,ik) 
             
             fwfisf(ji,jj) = qisf(ji,jj) / lfusisf          !fresh water flux kg/(m2s)                  
             fwfisf(ji,jj) = fwfisf(ji,jj) * ( soce / stbl(ji,jj) )
             !add to salinity trend
          ELSE
             qisf(ji,jj) = 0._wp ; fwfisf(ji,jj) = 0._wp
          END IF
       ENDDO
    ENDDO
    !
    IF( nn_timing == 1 )  CALL timing_stop('sbc_isf_bg03')
  END SUBROUTINE sbc_isf_bg03

   SUBROUTINE sbc_isf_cav( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_isf_cav  ***
      !!
      !! ** Purpose :   handle surface boundary condition under ice shelf
      !!
      !! ** Method  : -
      !!
      !! ** Action  :   utau, vtau : remain unchanged
      !!                taum, wndm : remain unchanged
      !!                qns        : update heat flux below ice shelf
      !!                emp, emps  : update freshwater flux below ice shelf
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)          ::   kt         ! ocean time step
      !
      LOGICAL :: ln_isomip = .true.
      REAL(wp), DIMENSION(:,:), POINTER       ::   zfrz,zpress,zti
      REAL(wp), DIMENSION(:,:), POINTER       ::   zgammat2d, zgammas2d 
      !REAL(wp), DIMENSION(:,:), POINTER ::   zqisf, zfwfisf
      REAL(wp) ::   zlamb1, zlamb2, zlamb3
      REAL(wp) ::   zeps1,zeps2,zeps3,zeps4,zeps6,zeps7
      REAL(wp) ::   zaqe,zbqe,zcqe,zaqer,zdis,zsfrz,zcfac
      REAL(wp) ::   zfwflx, zhtflx, zhtflx_b
      REAL(wp) ::   zgammat, zgammas
      REAL(wp) ::   zeps   =  -1.e-20_wp        !==   Local constant initialization   ==!
      INTEGER  ::   ji, jj     ! dummy loop indices
      INTEGER  ::   ii0, ii1, ij0, ij1   ! temporary integers
      INTEGER  ::   ierror     ! return error code
      LOGICAL  ::   lit=.TRUE.
      INTEGER  ::   nit
      !!---------------------------------------------------------------------
      !
      ! coeficient for linearisation of tfreez
      zlamb1=-0.0575
      zlamb2=0.0901
      zlamb3=-7.61e-04
      IF( nn_timing == 1 )  CALL timing_start('sbc_isf_cav')
      !
      CALL wrk_alloc( jpi,jpj, zfrz,zpress,zti, zgammat2d, zgammas2d )

      zcfac=0.0_wp 
      IF (ln_conserve)  zcfac=1.0_wp
      zpress(:,:)=0.0_wp
      zgammat2d(:,:)=0.0_wp
      zgammas2d(:,:)=0.0_wp
      !
      !
!CDIR COLLAPSE
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Crude approximation for pressure (but commonly used)
            ! 1e-04 to convert from Pa to dBar
            zpress(ji,jj)=grav*rau0*fsdepw(ji,jj,mikt(ji,jj))*1.e-04
            !
         END DO
      END DO

! Calculate in-situ temperature (ref to surface)
      zti(:,:)=tinsitu( ttbl, stbl, zpress )
! Calculate freezing temperature
      CALL eos_fzp( sss_m(:,:), zfrz(:,:), zpress )

      
      zhtflx=0._wp ; zfwflx=0._wp
      IF (nn_isfblk == 1) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF (mikt(ji,jj) > 1 ) THEN
                  nit = 1; lit = .TRUE.; zgammat=rn_gammat0; zgammas=rn_gammas0; zhtflx_b=0._wp
                  DO WHILE ( lit )
! compute gamma
                     CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx, ji, jj, lit)
! zhtflx is upward heat flux (out of ocean)
                     zhtflx = zgammat*rcp*rau0*(zti(ji,jj)-zfrz(ji,jj))
! zwflx is upward water flux
                     zfwflx = - zhtflx/lfusisf
! test convergence and compute gammat
                     IF ( (zhtflx - zhtflx_b) .LE. 0.01 ) lit = .FALSE.

                     nit = nit + 1
                     IF (nit .GE. 100) CALL ctl_stop( 'STOP', 'sbc_isf_hol99 : too many iteration ...' )

! save gammat and compute zhtflx_b
                     zgammat2d(ji,jj)=zgammat
                     zhtflx_b = zhtflx
                  END DO

                  qisf(ji,jj) = - zhtflx
! For genuine ISOMIP protocol this should probably be something like
                  fwfisf(ji,jj) = zfwflx  * ( soce / MAX(stbl(ji,jj),zeps))
               ELSE
                  fwfisf(ji,jj) = 0._wp
                  qisf(ji,jj)   = 0._wp
               END IF
            !
            END DO
         END DO

      ELSE IF (nn_isfblk == 2 ) THEN

! More complicated 3 equation thermodynamics as in MITgcm
!CDIR COLLAPSE
         DO jj = 2, jpj
            DO ji = 2, jpi
               IF (mikt(ji,jj) > 1 ) THEN
                  nit=1; lit=.TRUE.; zgammat=rn_gammat0; zgammas=rn_gammas0; zhtflx_b=0._wp; zhtflx=0._wp
                  DO WHILE ( lit )
                     CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx, ji, jj, lit)

                     zeps1=rcp*rau0*zgammat
                     zeps2=lfusisf*rau0*zgammas
                     zeps3=rhoisf*rcpi*kappa/risfdep(ji,jj)
                     zeps4=zlamb2+zlamb3*risfdep(ji,jj)
                     zeps6=zeps4-zti(ji,jj)
                     zeps7=zeps4-tsurf
                     zaqe=zlamb1 * (zeps1 + zeps3)
                     zaqer=0.5/zaqe
                     zbqe=zeps1*zeps6+zeps3*zeps7-zeps2
                     zcqe=zeps2*stbl(ji,jj)
                     zdis=zbqe*zbqe-4.0*zaqe*zcqe               
! Presumably zdis can never be negative because gammas is very small compared to gammat
                     zsfrz=(-zbqe-SQRT(zdis))*zaqer
                     IF (zsfrz .lt. 0.0) zsfrz=(-zbqe+SQRT(zdis))*zaqer
                     zfrz(ji,jj)=zeps4+zlamb1*zsfrz
  
! zfwflx is upward water flux
                     zfwflx= rau0 * zgammas * ( (zsfrz-stbl(ji,jj)) / zsfrz )
! zhtflx is upward heat flux (out of ocean)
! If non conservative we have zcfac=0.0 so zhtflx is as ISOMIP but with different zfrz value
                     zhtflx = ( zgammat*rau0 - zcfac*zfwflx ) * rcp * (zti(ji,jj) - zfrz(ji,jj) ) 
! zwflx is upward water flux
! If non conservative we have zcfac=0.0 so what follows is then zfwflx*sss_m/zsfrz
                     zfwflx = ( zgammas*rau0 - zcfac*zfwflx ) * (zsfrz - stbl(ji,jj)) / stbl(ji,jj)
! test convergence and compute gammat
                     IF (( zhtflx - zhtflx_b) .LE. 0.01 ) lit = .FALSE.

                     nit = nit + 1
                     IF (nit .GE. 51) THEN
                        WRITE(numout,*) "sbcisf : too many iteration ... ", &
                            &  zhtflx, zhtflx_b, zgammat, zgammas, nn_gammablk, ji, jj, mikt(ji,jj), narea
                        CALL ctl_stop( 'STOP', 'sbc_isf_hol99 : too many iteration ...' )
                     END IF
! save gammat and compute zhtflx_b
                     zgammat2d(ji,jj)=zgammat
                     zgammas2d(ji,jj)=zgammas
                     zhtflx_b = zhtflx

                  END DO
! If non conservative we have zcfac=0.0 so zhtflx is as ISOMIP but with different zfrz value
                  qisf(ji,jj) = - zhtflx 
! If non conservative we have zcfac=0.0 so what follows is then zfwflx*sss_m/zsfrz
                  fwfisf(ji,jj) = zfwflx 
               ELSE
                  fwfisf(ji,jj) = 0._wp
                  qisf(ji,jj)   = 0._wp
               ENDIF
               !
            END DO
         END DO
      ENDIF
      ! lbclnk
      CALL lbc_lnk(zgammas2d(:,:),'T',1.)
      CALL lbc_lnk(zgammat2d(:,:),'T',1.)
      ! output
      CALL iom_put('isfgammat', zgammat2d)
      CALL iom_put('isfgammas', zgammas2d)
      !
      CALL wrk_dealloc( jpi,jpj, zfrz,zpress,zti, zgammat2d, zgammas2d )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_isf_cav')

   END SUBROUTINE sbc_isf_cav

   SUBROUTINE sbc_isf_gammats(gt, gs, zqhisf, zqwisf, ji, jj, lit )
      !!----------------------------------------------------------------------
      !! ** Purpose    : compute the coefficient echange for heat flux
      !!
      !! ** Method     : gamma assume constant or depends of u* and stability
      !!
      !! ** References : Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
      !!                Jenkins et al., 2010, JPO, p2298-2312
      !!---------------------------------------------------------------------
      REAL(wp), INTENT(inout) :: gt, gs, zqhisf, zqwisf
      INTEGER , INTENT(in)    :: ji,jj
      LOGICAL , INTENT(inout) :: lit

      INTEGER  :: ikt                 ! loop index
      REAL(wp) :: zut, zvt, zustar           ! U, V at T point and friction velocity
      REAL(wp) :: zdku, zdkv                 ! U, V shear 
      REAL(wp) :: zPr, zSc, zRc              ! Prandtl, Scmidth and Richardson number 
      REAL(wp) :: zmob, zmols                ! Monin Obukov length, coriolis factor at T point
      REAL(wp) :: zbuofdep, zhnu             ! Bouyancy length scale, sublayer tickness
      REAL(wp) :: zhmax                      ! limitation of mol
      REAL(wp) :: zetastar                   ! stability parameter
      REAL(wp) :: zgmolet, zgmoles, zgturb   ! contribution of modelecular sublayer and turbulence 
      REAL(wp) :: zcoef                      ! temporary coef
      REAL(wp) :: zdep
      REAL(wp), PARAMETER :: zxsiN = 0.052   ! dimensionless constant
      REAL(wp), PARAMETER :: epsln = 1.0e-20 ! a small positive number
      REAL(wp), PARAMETER :: znu   = 1.95e-6 ! kinamatic viscosity of sea water (m2.s-1)
      REAL(wp) ::   rcs      = 1.0e-3_wp        ! conversion: mm/s ==> m/s
      REAL(wp), DIMENSION(2) :: zts, zab
      !!---------------------------------------------------------------------
      !
      IF( nn_gammablk == 0 ) THEN
      !! gamma is constant (specified in namelist)
         gt = rn_gammat0
         gs = rn_gammas0
         lit = .FALSE.
      ELSE IF ( nn_gammablk == 1 ) THEN
      !! gamma is assume to be proportional to u* 
      !! WARNING in case of Losh 2008 tbl parametrization, 
      !! you have to used the mean value of u in the boundary layer) 
      !! not yet coded
      !! Jenkins et al., 2010, JPO, p2298-2312
         ikt = mikt(ji,jj)
      !! Compute U and V at T points
   !      zut = 0.5 * ( utbl(ji-1,jj  ) + utbl(ji,jj) )
   !      zvt = 0.5 * ( vtbl(ji  ,jj-1) + vtbl(ji,jj) )
          zut = utbl(ji,jj)
          zvt = vtbl(ji,jj)

      !! compute ustar
         zustar = SQRT( rn_tfri2 * (zut * zut + zvt * zvt) )
      !! Compute mean value over the TBL

      !! Compute gammats
         gt = zustar * rn_gammat0
         gs = zustar * rn_gammas0
         lit = .FALSE.
      ELSE IF ( nn_gammablk == 2 ) THEN
      !! gamma depends of stability of boundary layer
      !! WARNING in case of Losh 2008 tbl parametrization, 
      !! you have to used the mean value of u in the boundary layer) 
      !! not yet coded
      !! Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
      !! as MOL depends of flux and flux depends of MOL, best will be iteration (TO DO)
               ikt = mikt(ji,jj)

      !! Compute U and V at T points
               zut = 0.5 * ( utbl(ji-1,jj  ) + utbl(ji,jj) )
               zvt = 0.5 * ( vtbl(ji  ,jj-1) + vtbl(ji,jj) )

      !! compute ustar
               zustar = SQRT( rn_tfri2 * (zut * zut + zvt * zvt) )
               IF (zustar == 0._wp) THEN           ! only for kt = 1 I think
                 gt = rn_gammat0
                 gs = rn_gammas0
               ELSE
      !! compute Rc number (as done in zdfric.F90)
               zcoef = 0.5 / fse3w(ji,jj,ikt)
               !                                            ! shear of horizontal velocity
               zdku = zcoef * (  un(ji-1,jj  ,ikt  ) + un(ji,jj,ikt  )   &
                  &             -un(ji-1,jj  ,ikt+1) - un(ji,jj,ikt+1)  )
               zdkv = zcoef * (  vn(ji  ,jj-1,ikt  ) + vn(ji,jj,ikt  )   &
                  &             -vn(ji  ,jj-1,ikt+1) - vn(ji,jj,ikt+1)  )
               !                                            ! richardson number (minimum value set to zero)
               zRc = rn2(ji,jj,ikt+1) / ( zdku*zdku + zdkv*zdkv + 1.e-20 )

      !! compute Pr and Sc number (can be improved)
               zPr =   13.8
               zSc = 2432.0

      !! compute gamma mole
               zgmolet = 12.5 * zPr ** (2.0/3.0) - 6.0
               zgmoles = 12.5 * zSc ** (2.0/3.0) -6.0

      !! compute bouyancy 
               zts(jp_tem) = ttbl(ji,jj)
               zts(jp_sal) = stbl(ji,jj)
               zdep        = fsdepw(ji,jj,ikt)
               !
               CALL eos_rab( zts, zdep, zab )
                  !
      !! compute length scale 
               zbuofdep = grav * ( zab(jp_tem) * zqhisf - zab(jp_sal) * zqwisf )  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! compute Monin Obukov Length
               ! Maximum boundary layer depth
               zhmax = fsdept(ji,jj,mbkt(ji,jj)) - fsdepw(ji,jj,mikt(ji,jj)) -0.001
               ! Compute Monin obukhov length scale at the surface and Ekman depth:
               zmob   = zustar ** 3 / (vkarmn * (zbuofdep + epsln))
               zmols  = SIGN(1._wp, zmob) * MIN(ABS(zmob), zhmax) * tmask(ji,jj,ikt)

      !! compute eta* (stability parameter)
               zetastar = 1 / ( SQRT(1 + MAX(zxsiN * zustar / ( ABS(ff(ji,jj)) * zmols * zRc ), 0.0)))

      !! compute the sublayer thickness
               zhnu = 5 * znu / zustar
      !! compute gamma turb
               zgturb = 1/vkarmn * LOG(zustar * zxsiN * zetastar * zetastar / ( ABS(ff(ji,jj)) * zhnu )) &
               &      + 1 / ( 2 * zxsiN * zetastar ) - 1/vkarmn

      !! compute gammats
               gt = zustar / (zgturb + zgmolet)
               gs = zustar / (zgturb + zgmoles)
               END IF
      END IF

   END SUBROUTINE

   SUBROUTINE sbc_isf_tbl( varin, varout, cptin )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE sbc_isf_tbl  ***
      !!
      !! ** Purpose : compute mean T/S/U/V in the boundary layer 
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in) :: varin
      REAL(wp), DIMENSION(:,:)  , INTENT(out):: varout
      
      CHARACTER(len=1), INTENT(in) :: cptin ! point of variable in/out

      REAL(wp) :: ze3, zhk
      REAL(wp), DIMENSION(:,:), POINTER :: zikt

      INTEGER :: ji,jj,jk
      INTEGER :: ikt,ikb
      INTEGER, DIMENSION(:,:), POINTER :: mkt, mkb

      CALL wrk_alloc( jpi,jpj, mkt, mkb  )
      CALL wrk_alloc( jpi,jpj, zikt )

      ! get first and last level of tbl
      mkt(:,:) = misfkt(:,:)
      mkb(:,:) = misfkb(:,:)

      varout(:,:)=0._wp
      DO jj = 2,jpj
         DO ji = 2,jpi
            IF (ssmask(ji,jj) == 1) THEN
               ikt = mkt(ji,jj)
               ikb = mkb(ji,jj)

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = fse3t_n(ji,jj,jk)
                  IF (cptin == 'T' ) varout(ji,jj) = varout(ji,jj) + varin(ji,jj,jk) * r1_hisf_tbl(ji,jj) * ze3
                  IF (cptin == 'U' ) varout(ji,jj) = varout(ji,jj) + 0.5_wp * (varin(ji,jj,jk) + varin(ji-1,jj,jk)) &
                     &                                                       * r1_hisf_tbl(ji,jj) * ze3
                  IF (cptin == 'V' ) varout(ji,jj) = varout(ji,jj) + 0.5_wp * (varin(ji,jj,jk) + varin(ji,jj-1,jk)) &
                     &                                                       * r1_hisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( fse3t_n(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)
               IF (cptin == 'T') &
                   &  varout(ji,jj) = varout(ji,jj) + varin(ji,jj,ikb) * (1._wp - zhk)
               IF (cptin == 'U') &
                   &  varout(ji,jj) = varout(ji,jj) + 0.5_wp * (varin(ji,jj,ikb) + varin(ji-1,jj,ikb)) * (1._wp - zhk)
               IF (cptin == 'V') &
                   &  varout(ji,jj) = varout(ji,jj) + 0.5_wp * (varin(ji,jj,ikb) + varin(ji,jj-1,ikb)) * (1._wp - zhk)
            END IF
         END DO
      END DO

      CALL wrk_dealloc( jpi,jpj, mkt, mkb )      
      CALL wrk_dealloc( jpi,jpj, zikt ) 

      IF (cptin == 'T') CALL lbc_lnk(varout,'T',1.)
      IF (cptin == 'U' .OR. cptin == 'V') CALL lbc_lnk(varout,'T',-1.)

   END SUBROUTINE sbc_isf_tbl
      

   SUBROUTINE sbc_isf_div( phdivn )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE sbc_isf_div  ***
      !!       
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :   
      !!                CAUTION : risf_tsc(:,:,jp_sal) is negative (outflow) increase the 
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phdivn   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikt, ikb 
      INTEGER  ::   nk_isf
      REAL(wp)     ::   zhk, z1_hisf_tbl, zhisf_tbl
      REAL(wp)     ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact   = 0.5_wp
      !
      IF (lk_vvl) THEN     ! need to re compute level distribution of isf fresh water
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkt(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               rhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), fse3t(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               ! test on tmask useless ?????
               DO jk = ikt, mbkt(ji,jj)
                  IF ( (SUM(fse3t(ji,jj,ikt:jk-1)) .LT. rhisf_tbl(ji,jj)) .AND. (tmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               rhisf_tbl(ji,jj) = MIN(rhisf_tbl(ji,jj), SUM(fse3t(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.
               misfkb(ji,jj) = ikb                                                  ! last wet level of the tbl
               r1_hisf_tbl(ji,jj) = 1._wp / rhisf_tbl(ji,jj)

               zhk           = SUM( fse3t(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)  ! proportion of tbl cover by cell from ikt to ikb - 1
               ralpha(ji,jj) = rhisf_tbl(ji,jj) * (1._wp - zhk ) / fse3t(ji,jj,ikb)  ! proportion of bottom cell influenced by boundary layer
            END DO
         END DO
      END IF ! vvl case
      !
      DO jj = 1,jpj
         DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)
               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) + ( fwfisf(ji,jj) + fwfisf_b(ji,jj) ) &
                    &               * r1_hisf_tbl(ji,jj) * r1_rau0 * zfact
               END DO
               ! level partially include in ice shelf boundary layer 
               phdivn(ji,jj,ikb) = phdivn(ji,jj,ikb) + ( fwfisf(ji,jj) &
                  &             + fwfisf_b(ji,jj) ) * r1_hisf_tbl(ji,jj) * r1_rau0 * zfact * ralpha(ji,jj) 
            !==   ice shelf melting mass distributed over several levels   ==!
         END DO
      END DO
      !
   END SUBROUTINE sbc_isf_div
                        
   FUNCTION tinsitu( ptem, psal, ppress ) RESULT( pti )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE eos_init  ***
      !!
      !! ** Purpose :   Compute the in-situ temperature [Celcius]
      !!
      !! ** Method  :   
      !!
      !! Reference  :   Bryden,h.,1973,deep-sea res.,20,401-408
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ptem   ! potential temperature [Celcius]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   psal   ! salinity             [psu]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ppress ! pressure             [dBar]
      REAL(wp), DIMENSION(:,:), POINTER           ::   pti    ! in-situ temperature [Celcius]
!      REAL(wp) :: fsatg
!      REAL(wp) :: pfps, pfpt, pfphp 
      REAL(wp) :: zt, zs, zp, zh, zq, zxk
      INTEGER  :: ji, jj            ! dummy loop indices
      !
      CALL wrk_alloc( jpi,jpj, pti  )
      ! 
      DO jj=1,jpj
         DO ji=1,jpi
            zh = ppress(ji,jj)
! Theta1
            zt = ptem(ji,jj)
            zs = psal(ji,jj)
            zp = 0.0
            zxk= zh * fsatg( zs, zt, zp )
            zt = zt + 0.5 * zxk
            zq = zxk
! Theta2
            zp = zp + 0.5 * zh
            zxk= zh*fsatg( zs, zt, zp )
            zt = zt + 0.29289322 * ( zxk - zq )
            zq = 0.58578644 * zxk + 0.121320344 * zq
! Theta3
            zxk= zh * fsatg( zs, zt, zp )
            zt = zt + 1.707106781 * ( zxk - zq )
            zq = 3.414213562 * zxk - 4.121320344 * zq
! Theta4
            zp = zp + 0.5 * zh
            zxk= zh * fsatg( zs, zt, zp )
            pti(ji,jj) = zt + ( zxk - 2.0 * zq ) / 6.0
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj, pti  )
      !
   END FUNCTION tinsitu
   !
   FUNCTION fsatg( pfps, pfpt, pfphp )
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION fsatg  ***
      !!
      !! ** Purpose    :   Compute the Adiabatic laspse rate [Celcius].[decibar]^-1
      !!
      !! ** Reference  :   Bryden,h.,1973,deep-sea res.,20,401-408
      !! 
      !! ** units      :   pressure        pfphp    decibars
      !!                   temperature     pfpt     deg celsius (ipts-68)
      !!                   salinity        pfps     (ipss-78)
      !!                   adiabatic       fsatg    deg. c/decibar
      !!----------------------------------------------------------------------
      REAL(wp) :: pfps, pfpt, pfphp 
      REAL(wp) :: fsatg
      !
      fsatg = (((-2.1687e-16*pfpt+1.8676e-14)*pfpt-4.6206e-13)*pfphp         &
        &    +((2.7759e-12*pfpt-1.1351e-10)*(pfps-35.)+((-5.4481e-14*pfpt    &
        &    +8.733e-12)*pfpt-6.7795e-10)*pfpt+1.8741e-8))*pfphp             &
        &    +(-4.2393e-8*pfpt+1.8932e-6)*(pfps-35.)                         &
        &    +((6.6228e-10*pfpt-6.836e-8)*pfpt+8.5258e-6)*pfpt+3.5803e-5
      !
    END FUNCTION fsatg
    !!======================================================================
END MODULE sbcisf
