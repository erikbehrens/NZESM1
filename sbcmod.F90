MODULE sbcmod
   !!======================================================================
   !!                       ***  MODULE  sbcmod  ***
   !! Surface module :  provide to the ocean its surface boundary condition
   !!======================================================================
   !! History :  3.0  ! 2006-07  (G. Madec)  Original code
   !!            3.1  ! 2008-08  (S. Masson, A. Caubel, E. Maisonnave, G. Madec) coupled interface
   !!            3.3  ! 2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.3  ! 2010-10  (S. Masson)  add diurnal cycle
   !!            3.3  ! 2010-09  (D. Storkey) add ice boundary conditions (BDY)
   !!             -   ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!             -   ! 2010-10  (J. Chanut, C. Bricaud, G. Madec)  add the surface pressure forcing
   !!            3.4  ! 2011-11  (C. Harris) CICE added as an option
   !!            3.5  ! 2012-11  (A. Coward, G. Madec) Rethink of heat, mass and salt surface fluxes
   !!            3.6  ! 2014-11  (P. Mathiot, C. Harris) add ice shelves melting                    
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_init       : read namsbc namelist
   !!   sbc            : surface ocean momentum, heat and freshwater boundary conditions
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers
   USE dom_oce          ! ocean space and time domain
   USE phycst           ! physical constants
   USE sbc_oce          ! Surface boundary condition: ocean fields
   USE trc_oce          ! shared ocean-passive tracers variables
   USE sbc_ice          ! Surface boundary condition: ice fields
   USE sbcdcy           ! surface boundary condition: diurnal cycle
   USE sbcssm           ! surface boundary condition: sea-surface mean variables
   USE sbcapr           ! surface boundary condition: atmospheric pressure
   USE sbcana           ! surface boundary condition: analytical formulation
   USE sbcflx           ! surface boundary condition: flux formulation
   USE sbcblk_clio      ! surface boundary condition: bulk formulation : CLIO
   USE sbcblk_core      ! surface boundary condition: bulk formulation : CORE
   USE sbcblk_mfs       ! surface boundary condition: bulk formulation : MFS
   USE sbcice_if        ! surface boundary condition: ice-if sea-ice model
   USE sbcice_lim       ! surface boundary condition: LIM 3.0 sea-ice model
   USE sbcice_lim_2     ! surface boundary condition: LIM 2.0 sea-ice model
   USE sbcice_cice      ! surface boundary condition: CICE    sea-ice model
   USE sbccpl           ! surface boundary condition: coupled florulation
   USE cpl_oasis3       ! OASIS routines for coupling
   USE sbcssr           ! surface boundary condition: sea surface restoring
   USE sbcrnf           ! surface boundary condition: runoffs
   USE sbcisf           ! surface boundary condition: ice shelf
   USE sbcfwb           ! surface boundary condition: freshwater budget
   USE closea           ! closed sea
   USE icbstp           ! Icebergs!

   USE prtctl           ! Print control                    (prt_ctl routine)
   USE iom              ! IOM library
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   USE sbcwave          ! Wave module
   USE bdy_par          ! Require lk_bdy

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc        ! routine called by step.F90
   PUBLIC   sbc_init   ! routine called by opa.F90
   
   INTEGER ::   nsbc   ! type of surface boundary condition (deduced from namsbc informations)
      
   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO-consortium (2011) 
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_init ***
      !!
      !! ** Purpose :   Initialisation of the ocean surface boundary computation
      !!
      !! ** Method  :   Read the namsbc namelist and set derived parameters
      !!                Call init routines for all other SBC modules that have one
      !!
      !! ** Action  : - read namsbc parameters
      !!              - nsbc: type of sbc
      !!----------------------------------------------------------------------
      INTEGER ::   icpt   ! local integer
      !!
      NAMELIST/namsbc/ nn_fsbc   , ln_ana    , ln_flx, ln_blk_clio, ln_blk_core, ln_mixcpl,   &
         &             ln_blk_mfs, ln_apr_dyn, nn_ice, nn_ice_embd, ln_dm2dc   , ln_rnf   ,   &
         &             ln_ssr    , nn_isf    , nn_fwb, ln_cdgw    , ln_wave    , ln_sdw   ,   &
         &             nn_lsm    , nn_limflx , nn_components, ln_cpl
      INTEGER  ::   ios
      INTEGER  ::   ierr, ierr0, ierr1, ierr2, ierr3, jpm
      LOGICAL  ::   ll_purecpl
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_init : surface boundary condition setting'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF

      REWIND( numnam_ref )              ! Namelist namsbc in reference namelist : Surface boundary
      READ  ( numnam_ref, namsbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsbc in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namsbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsbc )

      !                          ! overwrite namelist parameter using CPP key information
      IF( Agrif_Root() ) THEN                ! AGRIF zoom
        IF( lk_lim2 )   nn_ice      = 2
        IF( lk_lim3 )   nn_ice      = 3
        IF( lk_cice )   nn_ice      = 4
      ENDIF
      IF( cp_cfg == 'gyre' ) THEN            ! GYRE configuration
          ln_ana      = .TRUE.   
          nn_ice      =   0
      ENDIF
      !CEB  to avoid coupling of nest, since it crashes, since it needs specified in namecouple
            IF( .NOT. Agrif_Root() ) THEN
            lk_oasis=.FALSE.
            ENDIF

      IF(lwp) THEN               ! Control print
         WRITE(numout,*) '        Namelist namsbc (partly overwritten with CPP key setting)'
         WRITE(numout,*) '           frequency update of sbc (and ice)             nn_fsbc     = ', nn_fsbc
         WRITE(numout,*) '           Type of sbc : '
         WRITE(numout,*) '              analytical formulation                     ln_ana      = ', ln_ana
         WRITE(numout,*) '              flux       formulation                     ln_flx      = ', ln_flx
         WRITE(numout,*) '              CLIO bulk  formulation                     ln_blk_clio = ', ln_blk_clio
         WRITE(numout,*) '              CORE bulk  formulation                     ln_blk_core = ', ln_blk_core
         WRITE(numout,*) '              MFS  bulk  formulation                     ln_blk_mfs  = ', ln_blk_mfs
         WRITE(numout,*) '              ocean-atmosphere coupled formulation       ln_cpl      = ', ln_cpl
         WRITE(numout,*) '              forced-coupled mixed formulation           ln_mixcpl   = ', ln_mixcpl
         WRITE(numout,*) '              OASIS coupling (with atm or sas)           lk_oasis    = ', lk_oasis
         WRITE(numout,*) '              components of your executable              nn_components = ', nn_components
         WRITE(numout,*) '              Multicategory heat flux formulation (LIM3) nn_limflx   = ', nn_limflx
         WRITE(numout,*) '           Misc. options of sbc : '
         WRITE(numout,*) '              Patm gradient added in ocean & ice Eqs.    ln_apr_dyn  = ', ln_apr_dyn
         WRITE(numout,*) '              ice management in the sbc (=0/1/2/3)       nn_ice      = ', nn_ice 
         WRITE(numout,*) '              ice-ocean embedded/levitating (=0/1/2)     nn_ice_embd = ', nn_ice_embd
         WRITE(numout,*) '              daily mean to diurnal cycle qsr            ln_dm2dc    = ', ln_dm2dc 
         WRITE(numout,*) '              runoff / runoff mouths                     ln_rnf      = ', ln_rnf
         WRITE(numout,*) '              iceshelf formulation                       nn_isf      = ', nn_isf
         WRITE(numout,*) '              Sea Surface Restoring on SST and/or SSS    ln_ssr      = ', ln_ssr
         WRITE(numout,*) '              FreshWater Budget control  (=0/1/2)        nn_fwb      = ', nn_fwb
         WRITE(numout,*) '              closed sea (=0/1) (set in namdom)          nn_closea   = ', nn_closea
         WRITE(numout,*) '              n. of iterations if land-sea-mask applied  nn_lsm      = ', nn_lsm
      ENDIF

      ! LIM3 Multi-category heat flux formulation
      SELECT CASE ( nn_limflx)
      CASE ( -1 )
         IF(lwp) WRITE(numout,*) '              Use of per-category fluxes (nn_limflx = -1) '
      CASE ( 0  )
         IF(lwp) WRITE(numout,*) '              Average per-category fluxes (nn_limflx = 0) ' 
      CASE ( 1  )
         IF(lwp) WRITE(numout,*) '              Average then redistribute per-category fluxes (nn_limflx = 1) '
      CASE ( 2  )
         IF(lwp) WRITE(numout,*) '              Redistribute a single flux over categories (nn_limflx = 2) '
      END SELECT
      !
      IF ( nn_components /= jp_iam_nemo .AND. .NOT. lk_oasis )   &
         &      CALL ctl_stop( 'STOP', 'sbc_init : OPA-SAS coupled via OASIS, but key_oasis3 disabled' )
      IF ( nn_components == jp_iam_opa .AND. ln_cpl )   &
         &      CALL ctl_stop( 'STOP', 'sbc_init : OPA-SAS coupled via OASIS, but ln_cpl = T in OPA' )
      IF ( nn_components == jp_iam_opa .AND. ln_mixcpl )   &
         &      CALL ctl_stop( 'STOP', 'sbc_init : OPA-SAS coupled via OASIS, but ln_mixcpl = T in OPA' )
      IF ( ln_cpl .AND. .NOT. lk_oasis )    &
         &      CALL ctl_stop( 'STOP', 'sbc_init : OASIS-coupled atmosphere model, but key_oasis3 disabled' )
      IF( ln_mixcpl .AND. .NOT. lk_oasis )    &
         &      CALL ctl_stop( 'the forced-coupled mixed mode (ln_mixcpl) requires the cpp key key_oasis3' )
      IF( ln_mixcpl .AND. .NOT. ln_cpl )    &
         &      CALL ctl_stop( 'the forced-coupled mixed mode (ln_mixcpl) requires ln_cpl = T' )
      IF( ln_mixcpl .AND. nn_components /= jp_iam_nemo )    &
         &      CALL ctl_stop( 'the forced-coupled mixed mode (ln_mixcpl) is not yet working with sas-opa coupling via oasis' )

      !                              ! allocate sbc arrays
      IF( sbc_oce_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_init : unable to allocate sbc_oce arrays' )

      !                          ! Checks:
      IF( nn_isf .EQ. 0 ) THEN                      ! variable initialisation if no ice shelf 
         IF( sbc_isf_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_init : unable to allocate sbc_isf arrays' )
         fwfisf  (:,:)   = 0.0_wp ; fwfisf_b  (:,:)   = 0.0_wp
         risf_tsc(:,:,:) = 0.0_wp ; risf_tsc_b(:,:,:) = 0.0_wp
         rdivisf       = 0.0_wp
      END IF
      IF( nn_ice == 0 .AND. nn_components /= jp_iam_opa )   fr_i(:,:) = 0.e0 ! no ice in the domain, ice fraction is always zero

      sfx(:,:) = 0.0_wp                            ! the salt flux due to freezing/melting will be computed (i.e. will be non-zero) 
                                                   ! only if sea-ice is present
 
      fmmflx(:,:) = 0.0_wp                        ! freezing-melting array initialisation
      
      taum(:,:) = 0.0_wp                           ! Initialise taum for use in gls in case of reduced restart

      !                                            ! restartability   
      IF( ( nn_ice == 2 .OR. nn_ice ==3 ) .AND. .NOT.( ln_blk_clio .OR. ln_blk_core .OR. ln_cpl ) )   &
         &   CALL ctl_stop( 'LIM sea-ice model requires a bulk formulation or coupled configuration' )
      IF( nn_ice == 4 .AND. .NOT.( ln_blk_core .OR. ln_cpl ) )   &
         &   CALL ctl_stop( 'CICE sea-ice model requires ln_blk_core or ln_cpl' )
      IF( nn_ice == 4 .AND. lk_agrif )   &
         &   WRITE(numout,*) 'CICE + AGRIF CEB' 
      IF( ( nn_ice == 3 .OR. nn_ice == 4 ) .AND. nn_ice_embd == 0 )   &
         &   CALL ctl_stop( 'LIM3 and CICE sea-ice models require nn_ice_embd = 1 or 2' )
      IF( ( nn_ice /= 3 ) .AND. ( nn_limflx >= 0 ) )   &
         &   WRITE(numout,*) 'The nn_limflx>=0 option has no effect if sea ice model is not LIM3'
      IF( ( nn_ice == 3 ) .AND. ( ln_cpl ) .AND. ( ( nn_limflx == -1 ) .OR. ( nn_limflx == 1 ) ) )   &
         &   CALL ctl_stop( 'The chosen nn_limflx for LIM3 in coupled mode must be 0 or 2' )
      IF( ( nn_ice == 3 ) .AND. ( .NOT. ln_cpl ) .AND. ( nn_limflx == 2 ) )   &
         &   CALL ctl_stop( 'The chosen nn_limflx for LIM3 in forced mode cannot be 2' )

      IF( ln_dm2dc )   nday_qsr = -1   ! initialisation flag

      IF( ln_dm2dc .AND. .NOT.( ln_flx .OR. ln_blk_core ) .AND. nn_components /= jp_iam_opa )   &
         &   CALL ctl_stop( 'diurnal cycle into qsr field from daily values requires a flux or core-bulk formulation' )
      
      IF ( ln_wave ) THEN
      !Activated wave module but neither drag nor stokes drift activated
         IF ( .NOT.(ln_cdgw .OR. ln_sdw) )   THEN
            CALL ctl_warn( 'Ask for wave coupling but nor drag coefficient (ln_cdgw=F) neither stokes drift activated (ln_sdw=F)' )
      !drag coefficient read from wave model definable only with mfs bulk formulae and core 
         ELSEIF (ln_cdgw .AND. .NOT.(ln_blk_mfs .OR. ln_blk_core) )       THEN       
             CALL ctl_stop( 'drag coefficient read from wave model definable only with mfs bulk formulae and core')
         ENDIF
      ELSE
      IF ( ln_cdgw .OR. ln_sdw  )                                         & 
         &   CALL ctl_stop('Not Activated Wave Module (ln_wave=F) but     &
         & asked coupling with drag coefficient (ln_cdgw =T) or Stokes drift (ln_sdw=T) ')
      ENDIF 
      !                          ! Choice of the Surface Boudary Condition (set nsbc)
      ll_purecpl = ln_cpl .AND. .NOT. ln_mixcpl
      !
      icpt = 0
      IF( ln_ana          ) THEN   ;   nsbc = jp_ana     ; icpt = icpt + 1   ;   ENDIF       ! analytical           formulation
      IF( ln_flx          ) THEN   ;   nsbc = jp_flx     ; icpt = icpt + 1   ;   ENDIF       ! flux                 formulation
      IF( ln_blk_clio     ) THEN   ;   nsbc = jp_clio    ; icpt = icpt + 1   ;   ENDIF       ! CLIO bulk            formulation
      IF( ln_blk_core     ) THEN   ;   nsbc = jp_core    ; icpt = icpt + 1   ;   ENDIF       ! CORE bulk            formulation
      IF( ln_blk_mfs      ) THEN   ;   nsbc = jp_mfs     ; icpt = icpt + 1   ;   ENDIF       ! MFS  bulk            formulation
      IF( ll_purecpl      ) THEN   ;   nsbc = jp_purecpl ; icpt = icpt + 1   ;   ENDIF       ! Pure Coupled         formulation
      IF( cp_cfg == 'gyre') THEN   ;   nsbc = jp_gyre                        ;   ENDIF       ! GYRE analytical      formulation
      IF( nn_components == jp_iam_opa )   &
         &                  THEN   ;   nsbc = jp_none    ; icpt = icpt + 1   ;   ENDIF       ! opa coupling via SAS module
      IF( lk_esopa        )            nsbc = jp_esopa                                       ! esopa test, ALL formulations
      !
      IF( icpt /= 1 .AND. .NOT.lk_esopa ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '           E R R O R in setting the sbc, one and only one namelist/CPP key option '
         WRITE(numout,*) '                     must be choosen. You choose ', icpt, ' option(s)'
         WRITE(numout,*) '                     We stop'
         nstop = nstop + 1
      ENDIF
      IF(lwp) THEN
         WRITE(numout,*)
         IF( nsbc == jp_esopa   )   WRITE(numout,*) '              ESOPA test All surface boundary conditions'
         IF( nsbc == jp_gyre    )   WRITE(numout,*) '              GYRE analytical formulation'
         IF( nsbc == jp_ana     )   WRITE(numout,*) '              analytical formulation'
         IF( nsbc == jp_flx     )   WRITE(numout,*) '              flux formulation'
         IF( nsbc == jp_clio    )   WRITE(numout,*) '              CLIO bulk formulation'
         IF( nsbc == jp_core    )   WRITE(numout,*) '              CORE bulk formulation'
         IF( nsbc == jp_purecpl )   WRITE(numout,*) '              pure coupled formulation'
         IF( nsbc == jp_mfs     )   WRITE(numout,*) '              MFS Bulk formulation'
         IF( nsbc == jp_none    )   WRITE(numout,*) '              OPA coupled to SAS via oasis'
         IF( ln_mixcpl          )   WRITE(numout,*) '              + forced-coupled mixed formulation'
         IF( nn_components/= jp_iam_nemo )  &
            &                       WRITE(numout,*) '              + OASIS coupled SAS'
      ENDIF
      !
      IF( lk_oasis ) THEN
         IF( sbc_cpl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_cpl_alloc : unable to allocate arrays' )         
         CALL sbc_cpl_init (nn_ice)   ! OASIS initialisation. must be done before: (1) first time step
                                      !                                            (2) the use of nn_fsbc
      ENDIF

!     nn_fsbc initialization if OPA-SAS coupling via OASIS
!     sas model time step has to be declared in OASIS (mandatory) -> nn_fsbc has to be modified accordingly
      IF ( nn_components /= jp_iam_nemo ) THEN

         IF ( nn_components == jp_iam_opa ) nn_fsbc = cpl_freq('O_SFLX') / NINT(rdt)
         IF ( nn_components == jp_iam_sas ) nn_fsbc = cpl_freq('I_SFLX') / NINT(rdt)
         !
         IF(lwp)THEN
            WRITE(numout,*)
            WRITE(numout,*)"   OPA-SAS coupled via OASIS : nn_fsbc re-defined from OASIS namcouple ", nn_fsbc
            WRITE(numout,*)
         ENDIF
      ENDIF

      IF( MOD( nitend - nit000 + 1, nn_fsbc) /= 0 .OR.   &
          MOD( nstock             , nn_fsbc) /= 0 ) THEN 
         WRITE(ctmp1,*) 'experiment length (', nitend - nit000 + 1, ') or nstock (', nstock,   &
            &           ' is NOT a multiple of nn_fsbc (', nn_fsbc, ')'
         CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
      ENDIF
      !
      IF( MOD( rday, REAL(nn_fsbc, wp) * rdt ) /= 0 )   &
         &  CALL ctl_warn( 'nn_fsbc is NOT a multiple of the number of time steps in a day' )
      !
      IF( ln_dm2dc .AND. ( ( NINT(rday) / ( nn_fsbc * NINT(rdt) ) )  < 8 ) )   &
         &   CALL ctl_warn( 'diurnal cycle for qsr: the sampling of the diurnal cycle is too small...' )

                               CALL sbc_ssm_init               ! Sea-surface mean fields initialisation
      !
      IF( ln_ssr           )   CALL sbc_ssr_init               ! Sea-Surface Restoring initialisation
      !
                               CALL sbc_rnf_init               ! Runof initialisation
      !
      IF( nn_ice == 3      )   CALL sbc_lim_init               ! LIM3 initialisation

      IF( nn_ice == 4      )   CALL cice_sbc_init( nsbc )      ! CICE initialisation
      
   END SUBROUTINE sbc_init


   SUBROUTINE sbc( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!                condition (momentum, heat and freshwater fluxes)
      !!
      !! ** Method  :   blah blah  to be written ????????? 
      !!                CAUTION : never mask the surface stress field (tke sbc)
      !!
      !! ** Action  : - set the ocean surface boundary condition at before and now 
      !!                time step, i.e.  
      !!                utau_b, vtau_b, qns_b, qsr_b, emp_n, sfx_b, qrp_b, erp_b
      !!                utau  , vtau  , qns  , qsr  , emp  , sfx  , qrp  , erp
      !!              - updte the ice fraction : fr_i
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc')
      !
      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         utau_b(:,:) = utau(:,:)                         ! Swap the ocean forcing fields
         vtau_b(:,:) = vtau(:,:)                         ! (except at nit000 where before fields
         qns_b (:,:) = qns (:,:)                         !  are set at the end of the routine)
         ! The 3D heat content due to qsr forcing is treated in traqsr
         ! qsr_b (:,:) = qsr (:,:)
         emp_b(:,:) = emp(:,:)
         sfx_b(:,:) = sfx(:,:)
         IF ( ln_rnf ) THEN
            rnf_b    (:,:  ) = rnf    (:,:  )
            rnf_tsc_b(:,:,:) = rnf_tsc(:,:,:)
         ENDIF
      ENDIF
      !                                            ! ---------------------------------------- !
      !                                            !        forcing field computation         !
      !                                            ! ---------------------------------------- !
      !
      IF ( .NOT. lk_bdy ) then
         IF( ln_apr_dyn ) CALL sbc_apr( kt )                ! atmospheric pressure provided at kt+0.5*nn_fsbc
      ENDIF
                                                         ! (caution called before sbc_ssm)
      !
      IF( nn_components /= jp_iam_sas )   CALL sbc_ssm( kt )   ! ocean sea surface variables (sst_m, sss_m, ssu_m, ssv_m)
      !                                                        ! averaged over nf_sbc time-step

      IF (ln_wave) CALL sbc_wave( kt )
                                                   !==  sbc formulation  ==!
                                                            
      SELECT CASE( nsbc )                                ! Compute ocean surface boundary condition
      !                                                  ! (i.e. utau,vtau, qns, qsr, emp, sfx)
      CASE( jp_gyre  )   ;   CALL sbc_gyre    ( kt )                    ! analytical formulation : GYRE configuration
      CASE( jp_ana   )   ;   CALL sbc_ana     ( kt )                    ! analytical formulation : uniform sbc
      CASE( jp_flx   )   ;   CALL sbc_flx     ( kt )                    ! flux formulation
      CASE( jp_clio  )   ;   CALL sbc_blk_clio( kt )                    ! bulk formulation : CLIO for the ocean
      CASE( jp_core  )   
         IF( nn_components == jp_iam_sas ) &
            &                CALL sbc_cpl_rcv ( kt, nn_fsbc, nn_ice )   ! OPA-SAS coupling: SAS receiving fields from OPA 
                             CALL sbc_blk_core( kt )                    ! bulk formulation : CORE for the ocean
                                                                        ! from oce: sea surface variables (sst_m, sss_m,  ssu_m,  ssv_m)
      CASE( jp_purecpl )  ;  CALL sbc_cpl_rcv ( kt, nn_fsbc, nn_ice )   ! pure coupled formulation
                                                                        !
      CASE( jp_mfs   )   ;   CALL sbc_blk_mfs ( kt )                    ! bulk formulation : MFS for the ocean
      CASE( jp_none  ) 
         IF( nn_components == jp_iam_opa ) &
                             CALL sbc_cpl_rcv ( kt, nn_fsbc, nn_ice )   ! OPA-SAS coupling: OPA receiving fields from SAS
      CASE( jp_esopa )                                
                             CALL sbc_ana     ( kt )                    ! ESOPA, test ALL the formulations
                             CALL sbc_gyre    ( kt )                    !
                             CALL sbc_flx     ( kt )                    !
                             CALL sbc_blk_clio( kt )                    !
                             CALL sbc_blk_core( kt )                    !
                             CALL sbc_cpl_rcv ( kt, nn_fsbc, nn_ice )   !
      END SELECT

      IF( ln_mixcpl )        CALL sbc_cpl_rcv ( kt, nn_fsbc, nn_ice )   ! forced-coupled mixed formulation after forcing


      !                                            !==  Misc. Options  ==!
      
      SELECT CASE( nn_ice )                                       ! Update heat and freshwater fluxes over sea-ice areas
      CASE(  1 )   ;         CALL sbc_ice_if   ( kt )                ! Ice-cover climatology ("Ice-if" model)
      CASE(  2 )   ;         CALL sbc_ice_lim_2( kt, nsbc )          ! LIM-2 ice model
      CASE(  3 )   ;         CALL sbc_ice_lim  ( kt, nsbc )          ! LIM-3 ice model
      CASE(  4 )   ;         CALL sbc_ice_cice ( kt, nsbc )          ! CICE ice model
      END SELECT                                              

      IF( ln_icebergs    )   CALL icb_stp( kt )                   ! compute icebergs

      IF( nn_isf   /= 0  )   CALL sbc_isf( kt )                    ! compute iceshelves

      IF( ln_rnf         )   CALL sbc_rnf( kt )                   ! add runoffs to fresh water fluxes
 
      IF( ln_ssr         )   CALL sbc_ssr( kt )                   ! add SST/SSS damping term

      IF( nn_fwb    /= 0 )   CALL sbc_fwb( kt, nn_fwb, nn_fsbc )  ! control the freshwater budget

      IF( nn_closea == 1 )   CALL sbc_clo( kt )                   ! treatment of closed sea in the model domain 
      !                                                           ! (update freshwater fluxes)
!RBbug do not understand why see ticket 667
!clem: it looks like it is necessary for the north fold (in certain circumstances). Don't know why.
      CALL lbc_lnk( emp, 'T', 1. )
      !
      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         IF( ln_rstart .AND.    &                               !* Restart: read in restart file
            & iom_varid( numror, 'utau_b', ldstop = .FALSE. ) > 0 ) THEN 
            IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields red in the restart file'
            CALL iom_get( numror, jpdom_autoglo, 'utau_b', utau_b )   ! before i-stress  (U-point)
            CALL iom_get( numror, jpdom_autoglo, 'vtau_b', vtau_b )   ! before j-stress  (V-point)
            CALL iom_get( numror, jpdom_autoglo, 'qns_b' , qns_b  )   ! before non solar heat flux (T-point)
            ! The 3D heat content due to qsr forcing is treated in traqsr
            ! CALL iom_get( numror, jpdom_autoglo, 'qsr_b' , qsr_b  ) ! before     solar heat flux (T-point)
            CALL iom_get( numror, jpdom_autoglo, 'emp_b', emp_b  )    ! before     freshwater flux (T-point)
            ! To ensure restart capability with 3.3x/3.4 restart files    !! to be removed in v3.6
            IF( iom_varid( numror, 'sfx_b', ldstop = .FALSE. ) > 0 ) THEN
               CALL iom_get( numror, jpdom_autoglo, 'sfx_b', sfx_b )  ! before salt flux (T-point)
            ELSE
               sfx_b (:,:) = sfx(:,:)
            ENDIF
         ELSE                                                   !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields set to nit000'
            utau_b(:,:) = utau(:,:) 
            vtau_b(:,:) = vtau(:,:)
            qns_b (:,:) = qns (:,:)
            emp_b (:,:) = emp(:,:)
            sfx_b (:,:) = sfx(:,:)
         ENDIF
      ENDIF
      !                                                ! ---------------------------------------- !
      IF( lrst_oce ) THEN                              !      Write in the ocean restart file     !
         !                                             ! ---------------------------------------- !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : ocean surface forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrow, 'utau_b' , utau )
         CALL iom_rstput( kt, nitrst, numrow, 'vtau_b' , vtau )
         CALL iom_rstput( kt, nitrst, numrow, 'qns_b'  , qns  )
         ! The 3D heat content due to qsr forcing is treated in traqsr
         ! CALL iom_rstput( kt, nitrst, numrow, 'qsr_b'  , qsr  )
         CALL iom_rstput( kt, nitrst, numrow, 'emp_b'  , emp  )
         CALL iom_rstput( kt, nitrst, numrow, 'sfx_b'  , sfx  )
      ENDIF

      !                                                ! ---------------------------------------- !
      !                                                !        Outputs and control print         !
      !                                                ! ---------------------------------------- !
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN
         CALL iom_put( "empmr"  , emp    - rnf )                ! upward water flux
         CALL iom_put( "empbmr" , emp_b  - rnf )                ! before upward water flux ( needed to recalculate the time evolution of ssh in offline )
         CALL iom_put( "saltflx", sfx  )                        ! downward salt flux  
                                                                ! (includes virtual salt flux beneath ice 
                                                                ! in linear free surface case)
         CALL iom_put( "fmmflx", fmmflx  )                      ! Freezing-melting water flux
         CALL iom_put( "qt"    , qns  + qsr )                   ! total heat flux 
         CALL iom_put( "qns"   , qns        )                   ! solar heat flux
         CALL iom_put( "qsr"   ,       qsr  )                   ! solar heat flux
         IF( nn_ice > 0 .OR. nn_components == jp_iam_opa )   CALL iom_put( "ice_cover", fr_i )   ! ice fraction 
         CALL iom_put( "taum"  , taum       )                   ! wind stress module 
         CALL iom_put( "wspd"  , wndm       )                   ! wind speed  module over free ocean or leads in presence of sea-ice
      ENDIF
      !
      CALL iom_put( "utau", utau )   ! i-wind stress   (stress can be updated at 
      CALL iom_put( "vtau", vtau )   ! j-wind stress    each time step in sea-ice)
      !
      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab2d_1=fr_i              , clinfo1=' fr_i     - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=(emp-rnf + fwfisf), clinfo1=' emp-rnf  - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=(sfx-rnf + fwfisf), clinfo1=' sfx-rnf  - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=qns              , clinfo1=' qns      - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=qsr              , clinfo1=' qsr      - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab3d_1=tmask            , clinfo1=' tmask    - : ', mask1=tmask, ovlap=1, kdim=jpk )
         CALL prt_ctl(tab3d_1=tsn(:,:,:,jp_tem), clinfo1=' sst      - : ', mask1=tmask, ovlap=1, kdim=1   )
         CALL prt_ctl(tab3d_1=tsn(:,:,:,jp_sal), clinfo1=' sss      - : ', mask1=tmask, ovlap=1, kdim=1   )
         CALL prt_ctl(tab2d_1=utau             , clinfo1=' utau     - : ', mask1=umask,                      &
            &         tab2d_2=vtau             , clinfo2=' vtau     - : ', mask2=vmask, ovlap=1 )
      ENDIF

      IF( kt == nitend )   CALL sbc_final         ! Close down surface module if necessary
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc')
      !
   END SUBROUTINE sbc


   SUBROUTINE sbc_final
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_final  ***
      !!
      !! ** Purpose :   Finalize CICE (if used)
      !!---------------------------------------------------------------------
      !
      IF( nn_ice == 4 )   CALL cice_sbc_final
      !
   END SUBROUTINE sbc_final

   !!======================================================================
END MODULE sbcmod
