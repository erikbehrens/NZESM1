MODULE diahsb
   !!======================================================================
   !!                       ***  MODULE  diahsb  ***
   !! Ocean diagnostics: Heat, salt and volume budgets
   !!======================================================================
   !! History :  3.3  ! 2010-09  (M. Leclair)  Original code 
   !!                 ! 2012-10  (C. Rousset)  add iom_put
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_hsb       : Diagnose the conservation of ocean heat and salt contents, and volume
   !!   dia_hsb_rst   : Read or write DIA file in restart file
   !!   dia_hsb_init  : Initialization of the conservation diagnostic
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! surface thermohaline fluxes
   USE sbcrnf          ! river runoff
   USE sbcisf          ! ice shelves
   USE domvvl          ! vertical scale factors
   USE traqsr          ! penetrative solar radiation
   USE trabbc          ! bottom boundary condition 
   USE trabbc          ! bottom boundary condition
   USE bdy_par         ! (for lk_bdy)
   USE restart         ! ocean restart
   !
   USE iom             ! I/O manager
   USE in_out_manager  ! I/O manager
   USE lib_fortran     ! glob_sum
   USE lib_mpp         ! distributed memory computing library
   USE timing          ! preformance summary
   USE wrk_nemo        ! work arrays

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hsb        ! routine called by step.F90
   PUBLIC   dia_hsb_init   ! routine called by nemogcm.F90
   PUBLIC   dia_hsb_rst    ! routine called by step.F90

   LOGICAL, PUBLIC ::   ln_diahsb   !: check the heat and salt budgets

   REAL(wp) ::   surf_tot              ! ocean surface
   REAL(wp) ::   frc_t, frc_s, frc_v   ! global forcing trends
   REAL(wp) ::   frc_wn_t, frc_wn_s    ! global forcing trends
   !
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   surf          , ssh_ini          !
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   ssh_hc_loc_ini, ssh_sc_loc_ini   !
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   hc_loc_ini, sc_loc_ini, e3t_ini  !

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_hsb( kt )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!     
      !! ** Purpose: Compute the ocean global heat content, salt content and volume conservation
      !!	
      !! ** Method : - Compute the deviation of heat content, salt content and volume
      !!	            at the current time step from their values at nit000
      !!	            - Compute the contribution of forcing and remove it from these deviations
      !!
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER    ::   ji, jj, jk                  ! dummy loop indice
      REAL(wp)   ::   zdiff_hc    , zdiff_sc      ! heat and salt content variations
      REAL(wp)   ::   zdiff_hc1   , zdiff_sc1     !  -         -     -        - 
      REAL(wp)   ::   zdiff_v1    , zdiff_v2      ! volume variation
      REAL(wp)   ::   zerr_hc1    , zerr_sc1      ! heat and salt content misfit
      REAL(wp)   ::   zvol_tot                    ! volume
      REAL(wp)   ::   z_frc_trd_t , z_frc_trd_s   !    -     -
      REAL(wp)   ::   z_frc_trd_v                 !    -     -
      REAL(wp)   ::   z_wn_trd_t , z_wn_trd_s     !    -     -
      REAL(wp)   ::   z_ssh_hc , z_ssh_sc         !    -     -
      REAL(wp), DIMENSION(:,:), POINTER ::   z2d0, z2d1
      !!---------------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_hsb')      
      CALL wrk_alloc( jpi,jpj,   z2d0, z2d1 )
      !
      tsn(:,:,:,1) = tsn(:,:,:,1) * tmask(:,:,:) ; tsb(:,:,:,1) = tsb(:,:,:,1) * tmask(:,:,:) ;
      tsn(:,:,:,2) = tsn(:,:,:,2) * tmask(:,:,:) ; tsb(:,:,:,2) = tsb(:,:,:,2) * tmask(:,:,:) ;
      ! ------------------------- !
      ! 1 - Trends due to forcing !
      ! ------------------------- !
      z_frc_trd_v = r1_rau0 * glob_sum( - ( emp(:,:) - rnf(:,:) + fwfisf(:,:) ) * surf(:,:) ) ! volume fluxes
      z_frc_trd_t =           glob_sum( sbc_tsc(:,:,jp_tem) * surf(:,:) )                               ! heat fluxes
      z_frc_trd_s =           glob_sum( sbc_tsc(:,:,jp_sal) * surf(:,:) )                               ! salt fluxes
      ! Add runoff    heat & salt input
      IF( ln_rnf    )   z_frc_trd_t = z_frc_trd_t + glob_sum( rnf_tsc(:,:,jp_tem) * surf(:,:) )
      IF( ln_rnf_sal)   z_frc_trd_s = z_frc_trd_s + glob_sum( rnf_tsc(:,:,jp_sal) * surf(:,:) )
      ! Add ice shelf heat & salt input
      IF( nn_isf .GE. 1 )  THEN
          z_frc_trd_t = z_frc_trd_t + glob_sum( risf_tsc(:,:,jp_tem) * surf(:,:) )
          z_frc_trd_s = z_frc_trd_s + glob_sum( risf_tsc(:,:,jp_sal) * surf(:,:) )
      ENDIF

      ! Add penetrative solar radiation
      IF( ln_traqsr )   z_frc_trd_t = z_frc_trd_t + r1_rau0_rcp * glob_sum( qsr     (:,:) * surf(:,:) )
      ! Add geothermal heat flux
      IF( ln_trabbc )   z_frc_trd_t = z_frc_trd_t +               glob_sum( qgh_trd0(:,:) * surf(:,:) )
      !
      IF( .NOT. lk_vvl ) THEN
         IF ( ln_isfcav ) THEN
            DO ji=1,jpi
               DO jj=1,jpj
                  z2d0(ji,jj) = surf(ji,jj) * wn(ji,jj,mikt(ji,jj)) * tsb(ji,jj,mikt(ji,jj),jp_tem)
                  z2d1(ji,jj) = surf(ji,jj) * wn(ji,jj,mikt(ji,jj)) * tsb(ji,jj,mikt(ji,jj),jp_sal)
               ENDDO
            ENDDO
         ELSE
            z2d0(:,:) = surf(:,:) * wn(:,:,1) * tsb(:,:,1,jp_tem)
            z2d1(:,:) = surf(:,:) * wn(:,:,1) * tsb(:,:,1,jp_sal)
         END IF
         z_wn_trd_t = - glob_sum( z2d0 ) 
         z_wn_trd_s = - glob_sum( z2d1 )
      ENDIF

      frc_v = frc_v + z_frc_trd_v * rdt
      frc_t = frc_t + z_frc_trd_t * rdt
      frc_s = frc_s + z_frc_trd_s * rdt
      !                                          ! Advection flux through fixed surface (z=0)
      IF( .NOT. lk_vvl ) THEN
         frc_wn_t = frc_wn_t + z_wn_trd_t * rdt
         frc_wn_s = frc_wn_s + z_wn_trd_s * rdt
      ENDIF

      ! ------------------------ !
      ! 2 -  Content variations !
      ! ------------------------ !
      zdiff_v2 = 0._wp
      zdiff_hc = 0._wp
      zdiff_sc = 0._wp

      ! volume variation (calculated with ssh)
      zdiff_v1 = glob_sum( surf(:,:) * ( sshn(:,:) - ssh_ini(:,:) ) )

      ! heat & salt content variation (associated with ssh)
      IF( .NOT. lk_vvl ) THEN
         IF ( ln_isfcav ) THEN
            DO ji = 1, jpi
               DO jj = 1, jpj
                  z2d0(ji,jj) = surf(ji,jj) * ( tsn(ji,jj,mikt(ji,jj),jp_tem) * sshn(ji,jj) - ssh_hc_loc_ini(ji,jj) ) 
                  z2d1(ji,jj) = surf(ji,jj) * ( tsn(ji,jj,mikt(ji,jj),jp_sal) * sshn(ji,jj) - ssh_sc_loc_ini(ji,jj) ) 
               END DO
            END DO
         ELSE
            z2d0(:,:) = surf(:,:) * ( tsn(:,:,1,jp_tem) * sshn(:,:) - ssh_hc_loc_ini(:,:) ) 
            z2d1(:,:) = surf(:,:) * ( tsn(:,:,1,jp_sal) * sshn(:,:) - ssh_sc_loc_ini(:,:) ) 
         END IF
         z_ssh_hc = glob_sum( z2d0 ) 
         z_ssh_sc = glob_sum( z2d1 ) 
      ENDIF

      DO jk = 1, jpkm1
         ! volume variation (calculated with scale factors)
         zdiff_v2 = zdiff_v2 + glob_sum( surf(:,:) * tmask(:,:,jk) &
            &                           * ( fse3t_n(:,:,jk) - e3t_ini(:,:,jk) ) )
         ! heat content variation
         zdiff_hc = zdiff_hc + glob_sum(  surf(:,:) * tmask(:,:,jk) & 
            &                           * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_tem) - hc_loc_ini(:,:,jk) ) )
         ! salt content variation
         zdiff_sc = zdiff_sc + glob_sum(  surf(:,:) * tmask(:,:,jk)   &
            &                           * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_sal) - sc_loc_ini(:,:,jk) ) )
      ENDDO

      ! Substract forcing from heat content, salt content and volume variations
      zdiff_v1 = zdiff_v1 - frc_v
      IF( lk_vvl )   zdiff_v2 = zdiff_v2 - frc_v
      zdiff_hc = zdiff_hc - frc_t
      zdiff_sc = zdiff_sc - frc_s
      IF( .NOT. lk_vvl ) THEN
         zdiff_hc1 = zdiff_hc + z_ssh_hc 
         zdiff_sc1 = zdiff_sc + z_ssh_sc
         zerr_hc1  = z_ssh_hc - frc_wn_t
         zerr_sc1  = z_ssh_sc - frc_wn_s
      ENDIF

      ! ----------------------- !
      ! 3 - Diagnostics writing !
      ! ----------------------- !
      zvol_tot = 0._wp                    ! total ocean volume (calculated with scale factors)
      DO jk = 1, jpkm1
         zvol_tot  = zvol_tot + glob_sum( surf(:,:) * tmask(:,:,jk) * fse3t_n(:,:,jk) )
      END DO

!!gm to be added ?
!      IF( .NOT. lk_vvl ) THEN            ! fixed volume, add the ssh contribution
!        zvol_tot = zvol_tot + glob_sum( surf(:,:) * sshn(:,:) )
!      ENDIF
!!gm end

      IF( lk_vvl ) THEN
        CALL iom_put( 'bgtemper' , zdiff_hc / zvol_tot )              ! Temperature variation (C) 
        CALL iom_put( 'bgsaline' , zdiff_sc / zvol_tot )              ! Salinity    variation (psu)
        CALL iom_put( 'bgheatco' , zdiff_hc * 1.e-20 * rau0 * rcp )   ! Heat content variation (1.e20 J) 
        CALL iom_put( 'bgsaltco' , zdiff_sc * 1.e-9    )              ! Salt content variation (psu*km3)
        CALL iom_put( 'bgvolssh' , zdiff_v1 * 1.e-9    )              ! volume ssh variation (km3)  
        CALL iom_put( 'bgvole3t' , zdiff_v2 * 1.e-9    )              ! volume e3t variation (km3)  
        CALL iom_put( 'bgfrcvol' , frc_v    * 1.e-9    )              ! vol - surface forcing (km3) 
        CALL iom_put( 'bgfrctem' , frc_t / zvol_tot    )              ! hc  - surface forcing (C) 
        CALL iom_put( 'bgfrcsal' , frc_s / zvol_tot    )              ! sc  - surface forcing (psu) 
      ELSE
        CALL iom_put( 'bgtemper' , zdiff_hc1 / zvol_tot)              ! Heat content variation (C) 
        CALL iom_put( 'bgsaline' , zdiff_sc1 / zvol_tot)              ! Salt content variation (psu)
        CALL iom_put( 'bgheatco' , zdiff_hc1 * 1.e-20 * rau0 * rcp )  ! Heat content variation (1.e20 J) 
        CALL iom_put( 'bgsaltco' , zdiff_sc1 * 1.e-9    )             ! Salt content variation (psu*km3)
        CALL iom_put( 'bgvolssh' , zdiff_v1 * 1.e-9    )              ! volume ssh variation (km3)  
        CALL iom_put( 'bgfrcvol' , frc_v    * 1.e-9    )              ! vol - surface forcing (km3) 
        CALL iom_put( 'bgfrctem' , frc_t / zvol_tot    )              ! hc  - surface forcing (C) 
        CALL iom_put( 'bgfrcsal' , frc_s / zvol_tot    )              ! sc  - surface forcing (psu) 
        CALL iom_put( 'bgmistem' , zerr_hc1 / zvol_tot )              ! hc  - error due to free surface (C)
        CALL iom_put( 'bgmissal' , zerr_sc1 / zvol_tot )              ! sc  - error due to free surface (psu)
      ENDIF
      !
      IF( lrst_oce )   CALL dia_hsb_rst( kt, 'WRITE' )

      CALL wrk_dealloc( jpi,jpj,   z2d0, z2d1 )

      IF( nn_timing == 1 )   CALL timing_stop('dia_hsb')
      !
   END SUBROUTINE dia_hsb


   SUBROUTINE dia_hsb_rst( kt, cdrw )
     !!---------------------------------------------------------------------
     !!                   ***  ROUTINE limdia_rst  ***
     !!                     
     !! ** Purpose :   Read or write DIA file in restart file
     !!
     !! ** Method  :   use of IOM library
     !!----------------------------------------------------------------------
     INTEGER         , INTENT(in) ::   kt     ! ocean time-step
     CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
     !
     INTEGER ::   ji, jj, jk   ! dummy loop indices
     INTEGER ::   id1          ! local integers
     !!----------------------------------------------------------------------
     !
     IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
        IF( ln_rstart ) THEN                   !* Read the restart file
           !id1 = iom_varid( numror, 'frc_vol'  , ldstop = .FALSE. )
           !
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           IF(lwp) WRITE(numout,*) ' dia_hsb_rst at it= ', kt,' date= ', ndastp
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           CALL iom_get( numror, 'frc_v', frc_v )
           CALL iom_get( numror, 'frc_t', frc_t )
           CALL iom_get( numror, 'frc_s', frc_s )
           IF( .NOT. lk_vvl ) THEN
              CALL iom_get( numror, 'frc_wn_t', frc_wn_t )
              CALL iom_get( numror, 'frc_wn_s', frc_wn_s )
           ENDIF
           CALL iom_get( numror, jpdom_autoglo, 'ssh_ini', ssh_ini )
           CALL iom_get( numror, jpdom_autoglo, 'e3t_ini', e3t_ini )
           CALL iom_get( numror, jpdom_autoglo, 'hc_loc_ini', hc_loc_ini )
           CALL iom_get( numror, jpdom_autoglo, 'sc_loc_ini', sc_loc_ini )
           IF( .NOT. lk_vvl ) THEN
              CALL iom_get( numror, jpdom_autoglo, 'ssh_hc_loc_ini', ssh_hc_loc_ini )
              CALL iom_get( numror, jpdom_autoglo, 'ssh_sc_loc_ini', ssh_sc_loc_ini )
           ENDIF
       ELSE
          IF(lwp) WRITE(numout,*) '~~~~~~~'
          IF(lwp) WRITE(numout,*) ' dia_hsb at initial state '
          IF(lwp) WRITE(numout,*) '~~~~~~~'
          ssh_ini(:,:) = sshn(:,:)                                       ! initial ssh
          DO jk = 1, jpk
             e3t_ini   (:,:,jk) = fse3t_n(:,:,jk)                        ! initial vertical scale factors
             hc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_tem) * fse3t_n(:,:,jk)   ! initial heat content
             sc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_sal) * fse3t_n(:,:,jk)   ! initial salt content
          END DO
          frc_v = 0._wp                                           ! volume       trend due to forcing
          frc_t = 0._wp                                           ! heat content   -    -   -    -   
          frc_s = 0._wp                                           ! salt content   -    -   -    -        
          IF( .NOT. lk_vvl ) THEN
             IF ( ln_isfcav ) THEN
                DO ji=1,jpi
                   DO jj=1,jpj
                      ssh_hc_loc_ini(ji,jj) = tsn(ji,jj,mikt(ji,jj),jp_tem) * sshn(ji,jj)   ! initial heat content in ssh
                      ssh_sc_loc_ini(ji,jj) = tsn(ji,jj,mikt(ji,jj),jp_sal) * sshn(ji,jj)   ! initial salt content in ssh
                   ENDDO
                ENDDO
             ELSE
                ssh_hc_loc_ini(:,:) = tsn(:,:,1,jp_tem) * sshn(:,:)   ! initial heat content in ssh
                ssh_sc_loc_ini(:,:) = tsn(:,:,1,jp_sal) * sshn(:,:)   ! initial salt content in ssh
             END IF
             frc_wn_t = 0._wp                                       ! initial heat content misfit due to free surface
             frc_wn_s = 0._wp                                       ! initial salt content misfit due to free surface
          ENDIF
       ENDIF

     ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
        !                                   ! -------------------
        IF(lwp) WRITE(numout,*) '~~~~~~~'
        IF(lwp) WRITE(numout,*) ' dia_hsb_rst at it= ', kt,' date= ', ndastp
        IF(lwp) WRITE(numout,*) '~~~~~~~'

        CALL iom_rstput( kt, nitrst, numrow, 'frc_v'   , frc_v     )
        CALL iom_rstput( kt, nitrst, numrow, 'frc_t'   , frc_t     )
        CALL iom_rstput( kt, nitrst, numrow, 'frc_s'   , frc_s     )
        IF( .NOT. lk_vvl ) THEN
           CALL iom_rstput( kt, nitrst, numrow, 'frc_wn_t', frc_wn_t )
           CALL iom_rstput( kt, nitrst, numrow, 'frc_wn_s', frc_wn_s )
        ENDIF
        CALL iom_rstput( kt, nitrst, numrow, 'ssh_ini', ssh_ini )
        CALL iom_rstput( kt, nitrst, numrow, 'e3t_ini', e3t_ini )
        CALL iom_rstput( kt, nitrst, numrow, 'hc_loc_ini', hc_loc_ini )
        CALL iom_rstput( kt, nitrst, numrow, 'sc_loc_ini', sc_loc_ini )
        IF( .NOT. lk_vvl ) THEN
           CALL iom_rstput( kt, nitrst, numrow, 'ssh_hc_loc_ini', ssh_hc_loc_ini )
           CALL iom_rstput( kt, nitrst, numrow, 'ssh_sc_loc_ini', ssh_sc_loc_ini )
        ENDIF
        !
     ENDIF
     !
   END SUBROUTINE dia_hsb_rst


   SUBROUTINE dia_hsb_init
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!     
      !! ** Purpose: Initialization for the heat salt volume budgets
      !!	
      !! ** Method : Compute initial heat content, salt content and volume
      !!
      !! ** Action : - Compute initial heat content, salt content and volume
      !!             - Initialize forcing trends
      !!             - Compute coefficients for conversion
      !!---------------------------------------------------------------------------
      INTEGER ::   jk       ! dummy loop indice
      INTEGER ::   ierror   ! local integer
      INTEGER ::   ios
      !
      NAMELIST/namhsb/ ln_diahsb
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_hsb_init : check the heat and salt budgets'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF

      REWIND( numnam_ref )              ! Namelist namhsb in reference namelist
      READ  ( numnam_ref, namhsb, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namhsb in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namhsb in configuration namelist
      READ  ( numnam_cfg, namhsb, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namhsb in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namhsb )

      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_hsb_init : check the heat and salt budgets'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namhsb : set hsb parameters'
         WRITE(numout,*) '      Switch for hsb diagnostic (T) or not (F)  ln_diahsb  = ', ln_diahsb
         WRITE(numout,*)
      ENDIF

      IF( .NOT. ln_diahsb )   RETURN
         !      IF( .NOT. lk_mpp_rep ) &
         !        CALL ctl_stop (' Your global mpp_sum if performed in single precision - 64 bits -', &
         !             &         ' whereas the global sum to be precise must be done in double precision ',&
         !             &         ' please add key_mpp_rep')

      ! ------------------- !
      ! 1 - Allocate memory !
      ! ------------------- !
      ALLOCATE( hc_loc_ini(jpi,jpj,jpk), sc_loc_ini(jpi,jpj,jpk), &
         &      e3t_ini(jpi,jpj,jpk), surf(jpi,jpj),  ssh_ini(jpi,jpj), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb: unable to allocate hc_loc_ini' )   ;   RETURN
      ENDIF

      IF(.NOT. lk_vvl ) ALLOCATE( ssh_hc_loc_ini(jpi,jpj), ssh_sc_loc_ini(jpi,jpj),STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb: unable to allocate hc_loc_ini' )   ;   RETURN
      ENDIF

      ! ----------------------------------------------- !
      ! 2 - Time independant variables and file opening !
      ! ----------------------------------------------- !
      IF(lwp) WRITE(numout,*) "dia_hsb: heat salt volume budgets activated"
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      surf(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)      ! masked surface grid cell area
      surf_tot  = glob_sum( surf(:,:) )                                       ! total ocean surface area

      IF( lk_bdy ) CALL ctl_warn( 'dia_hsb does not take open boundary fluxes into account' )         
      !
      ! ---------------------------------- !
      ! 4 - initial conservation variables !
      ! ---------------------------------- !
      CALL dia_hsb_rst( nit000, 'READ' )  !* read or initialize all required files
      !
   END SUBROUTINE dia_hsb_init

   !!======================================================================
END MODULE diahsb
