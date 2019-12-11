MODULE sbcblk_clio
   !!======================================================================
   !!                   ***  MODULE  sbcblk_clio  ***
   !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
   !!=====================================================================
   !! History :  OPA  !  1997-06 (Louvain-La-Neuve)  Original code
   !!                 !  2001-04 (C. Ethe) add flx_blk_declin
   !!   NEMO     2.0  !  2002-08 (C. Ethe, G. Madec) F90: Free form and module
   !!            3.0  !  2008-03 (C. Talandier, G. Madec) surface module + LIM3
   !!            3.2  !  2009-04 (B. Lemaire) Introduce iom_put
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_clio     : CLIO bulk formulation: read and update required input fields
   !!   blk_clio_oce     : ocean CLIO bulk formulea: compute momentum, heat and freswater fluxes for the ocean
   !!   blk_ice_clio     : ice   CLIO bulk formulea: compute momentum, heat and freswater fluxes for the sea-ice
   !!   blk_clio_qsr_oce : shortwave radiation for ocean computed from the cloud cover
   !!   blk_clio_qsr_ice : shortwave radiation for ice   computed from the cloud cover
   !!   flx_blk_declin   : solar declination
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE fldread        ! read input fields
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE wrk_nemo       ! work arrays
   USE timing         ! Timing
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   USE albedo
   USE prtctl          ! Print control
#if defined key_lim3 
   USE ice
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE limthd_dh       ! for CALL lim_thd_snwblow
#elif defined key_lim2
   USE ice_2
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE par_ice_2       ! Surface boundary condition: ice fields
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_blk_clio        ! routine called by sbcmod.F90 
#if defined key_lim2 || defined key_lim3
   PUBLIC blk_ice_clio_tau    ! routine called by sbcice_lim.F90 
   PUBLIC blk_ice_clio_flx    ! routine called by sbcice_lim.F90 
#endif

   INTEGER , PARAMETER ::   jpfld   = 7           ! maximum number of files to read 
   INTEGER , PARAMETER ::   jp_utau = 1           ! index of wind stress (i-component)      (N/m2)    at U-point
   INTEGER , PARAMETER ::   jp_vtau = 2           ! index of wind stress (j-component)      (N/m2)    at V-point
   INTEGER , PARAMETER ::   jp_wndm = 3           ! index of 10m wind module                 (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_humi = 4           ! index of specific humidity               ( % )
   INTEGER , PARAMETER ::   jp_ccov = 5           ! index of cloud cover                     ( % )
   INTEGER , PARAMETER ::   jp_tair = 6           ! index of 10m air temperature             (Kelvin)
   INTEGER , PARAMETER ::   jp_prec = 7           ! index of total precipitation (rain+snow) (Kg/m2/s)

   TYPE(FLD),ALLOCATABLE,DIMENSION(:) :: sf  ! structure of input fields (file informations, fields read)

   INTEGER, PARAMETER  ::   jpintsr = 24          ! number of time step between sunrise and sunset
   !                                              ! uses for heat flux computation
   LOGICAL ::   lbulk_init = .TRUE.               ! flag, bulk initialization done or not)

   REAL(wp) ::   cai = 1.40e-3 ! best estimate of atm drag in order to get correct FS export in ORCA2-LIM
   REAL(wp) ::   cao = 1.00e-3 ! chosen by default  ==> should depends on many things...  !!gmto be updated

   REAL(wp) ::   rdtbs2      !:   
   
   REAL(wp), DIMENSION(19)  ::  budyko            ! BUDYKO's coefficient (cloudiness effect on LW radiation)
   DATA budyko / 1.00, 0.98, 0.95, 0.92, 0.89, 0.86, 0.83, 0.80, 0.78, 0.75,   &
      &          0.72, 0.69, 0.67, 0.64, 0.61, 0.58, 0.56, 0.53, 0.50 /
   REAL(wp), DIMENSION(20)  :: tauco              ! cloud optical depth coefficient
   DATA tauco / 6.6, 6.6, 7.0, 7.2, 7.1, 6.8, 6.5, 6.6, 7.1, 7.6,   &
      &         6.6, 6.1, 5.6, 5.5, 5.8, 5.8, 5.6, 5.6, 5.6, 5.6 /
   !!
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sbudyko      ! cloudiness effect on LW radiation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   stauc        ! cloud optical depth 
   
   REAL(wp) ::   eps20  = 1.e-20   ! constant values
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_blk_clio( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_clio  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!      (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  : (1) READ each fluxes in NetCDF files:
      !!      the i-component of the stress                (N/m2)
      !!      the j-component of the stress                (N/m2)
      !!      the 10m wind speed module                    (m/s)
      !!      the 10m air temperature                      (Kelvin)
      !!      the 10m specific humidity                    (%)
      !!      the cloud cover                              (%)
      !!      the total precipitation (rain+snow)          (Kg/m2/s)
      !!              (2) CALL blk_oce_clio
      !!
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point over free ocean or leads in presence of sea-ice
      !!              - qns         non-solar heat flux including latent heat of solid 
      !!                            precip. melting and emp heat content
      !!              - qsr         solar heat flux
      !!              - emp         upward mass flux (evap. - precip)
      !!              - sfx         salt flux; set to zero at nit000 but possibly non-zero
      !!                            if ice is present (computed in limsbc(_2).F90)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ifpr, jfpr                   ! dummy indices
      INTEGER  ::   ierr0, ierr1, ierr2, ierr3   ! return error code
      INTEGER  ::   ios                          ! Local integer output status for namelist read
      !!
      CHARACTER(len=100) ::  cn_dir                            !   Root directory for location of CLIO files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                 ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_utau, sn_vtau, sn_wndm, sn_tair      ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_humi, sn_ccov, sn_prec               !   "                                 "
      !!
      NAMELIST/namsbc_clio/ cn_dir, sn_utau, sn_vtau, sn_wndm, sn_humi,   &
         &                          sn_ccov, sn_tair, sn_prec
      !!---------------------------------------------------------------------

      !                                         ! ====================== !
      IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !

         REWIND( numnam_ref )              ! Namelist namsbc_clio in reference namelist : CLIO files
         READ  ( numnam_ref, namsbc_clio, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_clio in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_clio in configuration namelist : CLIO files
         READ  ( numnam_cfg, namsbc_clio, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_clio in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_clio )

         ! store namelist information in an array
         slf_i(jp_utau) = sn_utau   ;   slf_i(jp_vtau) = sn_vtau   ;   slf_i(jp_wndm) = sn_wndm
         slf_i(jp_tair) = sn_tair   ;   slf_i(jp_humi) = sn_humi
         slf_i(jp_ccov) = sn_ccov   ;   slf_i(jp_prec) = sn_prec
         
         ! set sf structure
         ALLOCATE( sf(jpfld), STAT=ierr0 )
         IF( ierr0 > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_clio: unable to allocate sf structure' )
         DO ifpr= 1, jpfld
            ALLOCATE( sf(ifpr)%fnow(jpi,jpj,1) , STAT=ierr1)
            IF( slf_i(ifpr)%ln_tint ) ALLOCATE( sf(ifpr)%fdta(jpi,jpj,1,2) , STAT=ierr2 )
         END DO
         IF( ierr1+ierr2 > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_clio: unable to allocate sf array structure' )
         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_clio', 'flux formulation for ocean surface boundary condition', 'namsbc_clio' )
         
         ! allocate sbcblk clio arrays
         ALLOCATE( sbudyko(jpi,jpj) , stauc(jpi,jpj), STAT=ierr3 )
         IF( ierr3 > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_clio: unable to allocate arrays' )
         !
         sfx(:,:) = 0._wp                       ! salt flux; zero unless ice is present (computed in limsbc(_2).F90)
         !
      ENDIF
      !                                         ! ====================== !
      !                                         !    At each time-step   !
      !                                         ! ====================== !
      !
      CALL fld_read( kt, nn_fsbc, sf )                ! input fields provided at the current time-step
      !
      IF( MOD( kt - 1, nn_fsbc ) == 0 )   CALL blk_oce_clio( sf, sst_m )
      !
   END SUBROUTINE sbc_blk_clio


   SUBROUTINE blk_oce_clio( sf, pst )
      !!---------------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce_clio  ***
      !!                 
      !!  ** Purpose :   Compute momentum, heat and freshwater fluxes at ocean surface
      !!               using CLIO bulk formulea
      !!         
      !!  ** Method  :   The flux of heat at the ocean surfaces are derived
      !!       from semi-empirical ( or bulk ) formulae which relate the flux to 
      !!       the properties of the surface and of the lower atmosphere. Here, we
      !!       follow the work of Oberhuber, 1988   
      !!               - momentum flux (stresses) directly read in files at U- and V-points
      !!               - compute ocean/ice albedos (call albedo_oce/albedo_ice)  
      !!               - compute shortwave radiation for ocean (call blk_clio_qsr_oce)
      !!               - compute long-wave radiation for the ocean
      !!               - compute the turbulent heat fluxes over the ocean
      !!               - deduce the evaporation over the ocean
      !!  ** Action  :   Fluxes over the ocean:
      !!               - utau, vtau  i- and j-component of the wind stress
      !!               - taum        wind stress module at T-point
      !!               - wndm        10m wind module at T-point over free ocean or leads in presence of sea-ice
      !!               - qns         non-solar heat flux including latent heat of solid 
      !!                             precip. melting and emp heat content
      !!               - qsr         solar heat flux
      !!               - emp         suface mass flux (evap.-precip.)
      !!  ** Nota    :   sf has to be a dummy argument for AGRIF on NEC
      !!----------------------------------------------------------------------
      TYPE(fld), INTENT(in), DIMENSION(:)       ::   sf    ! input data
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj) ::   pst   ! surface temperature                      [Celcius]
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!
      REAL(wp) ::   zrhova, zcsho, zcleo, zcldeff               ! temporary scalars
      REAL(wp) ::   zqsato, zdteta, zdeltaq, ztvmoy, zobouks    !    -         -
      REAL(wp) ::   zpsims, zpsihs, zpsils, zobouku, zxins, zpsimu   !    -         -
      REAL(wp) ::   zpsihu, zpsilu, zstab,zpsim, zpsih, zpsil   !    -         -
      REAL(wp) ::   zvatmg, zcmn, zchn, zcln, zcmcmn, zdenum    !    -         -
      REAL(wp) ::   zdtetar, ztvmoyr, zlxins, zchcm, zclcm      !    -         -
      REAL(wp) ::   zmt1, zmt2, zmt3, ztatm3, ztamr, ztaevbk    !    -         -
      REAL(wp) ::   zsst, ztatm, zcco1, zpatm, zcmax, zrmax     !    -         -
      REAL(wp) ::   zrhoa, zev, zes, zeso, zqatm, zevsqr        !    -         -
      REAL(wp) ::   ztx2, zty2, zcevap, zcprec                  !    -         -
      REAL(wp), POINTER, DIMENSION(:,:) ::   zqlw        ! long-wave heat flux over ocean
      REAL(wp), POINTER, DIMENSION(:,:) ::   zqla        ! latent heat flux over ocean
      REAL(wp), POINTER, DIMENSION(:,:) ::   zqsb        ! sensible heat flux over ocean
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_oce_clio')
      !
      CALL wrk_alloc( jpi,jpj, zqlw, zqla, zqsb )

      zpatm = 101000._wp      ! atmospheric pressure  (assumed constant here)

      !------------------------------------!
      !   momentum fluxes  (utau, vtau )   !
      !------------------------------------!
!CDIR COLLAPSE
      utau(:,:) = sf(jp_utau)%fnow(:,:,1)
!CDIR COLLAPSE
      vtau(:,:) = sf(jp_vtau)%fnow(:,:,1)

      !------------------------------------!
      !   wind stress module (taum )       !
      !------------------------------------!
!CDIR NOVERRCHK
      DO jj = 2, jpjm1
!CDIR NOVERRCHK
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ztx2 = utau(ji-1,jj  ) + utau(ji,jj)
            zty2 = vtau(ji  ,jj-1) + vtau(ji,jj)
            taum(ji,jj) = 0.5 * SQRT( ztx2 * ztx2 + zty2 * zty2 )
         END DO
      END DO
      utau(:,:) = utau(:,:) * umask(:,:,1)
      vtau(:,:) = vtau(:,:) * vmask(:,:,1)
      taum(:,:) = taum(:,:) * tmask(:,:,1)
      CALL lbc_lnk( taum, 'T', 1. )

      !------------------------------------!
      !   store the wind speed  (wndm )    !
      !------------------------------------!
!CDIR COLLAPSE
      wndm(:,:) = sf(jp_wndm)%fnow(:,:,1)
      wndm(:,:) = wndm(:,:) * tmask(:,:,1)

      !------------------------------------------------!
      !   Shortwave radiation for ocean and snow/ice   !
      !------------------------------------------------!
      
      CALL blk_clio_qsr_oce( qsr )
      qsr(:,:) = qsr(:,:) * tmask(:,:,1) ! no shortwave radiation into the ocean beneath ice shelf
      !------------------------!
      !   Other ocean fluxes   !
      !------------------------!
!CDIR NOVERRCHK
!CDIR COLLAPSE
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            !
            zsst  = pst(ji,jj)              + rt0           ! converte Celcius to Kelvin the SST
            ztatm = sf(jp_tair)%fnow(ji,jj,1)               ! and set minimum value far above 0 K (=rt0 over land)
            zcco1 = 1.0 - sf(jp_ccov)%fnow(ji,jj,1)         ! fraction of clear sky ( 1 - cloud cover)
            zrhoa = zpatm / ( 287.04 * ztatm )              ! air density (equation of state for dry air) 
            ztamr = ztatm - rtt                             ! Saturation water vapour
            zmt1  = SIGN( 17.269,  ztamr )                  !           ||
            zmt2  = SIGN( 21.875,  ztamr )                  !          \  /
            zmt3  = SIGN( 28.200, -ztamr )                  !           \/
            zes   = 611.0 * EXP(  ABS( ztamr ) * MIN ( zmt1, zmt2 ) / ( ztatm - 35.86  + MAX( 0.e0, zmt3 ) )  )
            zev    = sf(jp_humi)%fnow(ji,jj,1) * zes        ! vapour pressure  
            zevsqr = SQRT( zev * 0.01 )                     ! square-root of vapour pressure
            zqatm = 0.622 * zev / ( zpatm - 0.378 * zev )   ! specific humidity 

            !--------------------------------------!
            !  long-wave radiation over the ocean  !  ( Berliand 1952 ; all latitudes )
            !--------------------------------------!
            ztatm3  = ztatm * ztatm * ztatm
            zcldeff = 1.0 - sbudyko(ji,jj) * sf(jp_ccov)%fnow(ji,jj,1) * sf(jp_ccov)%fnow(ji,jj,1)    
            ztaevbk = ztatm * ztatm3 * zcldeff * ( 0.39 - 0.05 * zevsqr ) 
            !
            zqlw(ji,jj) = - emic * stefan * ( ztaevbk + 4. * ztatm3 * ( zsst - ztatm ) ) 

            !--------------------------------------------------
            !  Latent and sensible heat fluxes over the ocean 
            !--------------------------------------------------
            !                                                          ! vapour pressure at saturation of ocean
            zeso =  611.0 * EXP ( 17.2693884 * ( zsst - rtt ) * tmask(ji,jj,1) / ( zsst - 35.86 ) )

            zqsato = ( 0.622 * zeso ) / ( zpatm - 0.378 * zeso )       ! humidity close to the ocean surface (at saturation)

            ! Drag coefficients from Large and Pond (1981,1982)
            !                                                          ! Stability parameters
            zdteta  = zsst - ztatm
            zdeltaq = zqatm - zqsato
            ztvmoy  = ztatm * ( 1. + 2.2e-3 * ztatm * zqatm )
            zdenum  = MAX( sf(jp_wndm)%fnow(ji,jj,1) * sf(jp_wndm)%fnow(ji,jj,1) * ztvmoy, eps20 )
            zdtetar = zdteta / zdenum
            ztvmoyr = ztvmoy * ztvmoy * zdeltaq / zdenum
            !                                                          ! case of stable atmospheric conditions
            zobouks = -70.0 * 10. * ( zdtetar + 3.2e-3 * ztvmoyr )
            zobouks = MAX( 0.e0, zobouks )
            zpsims = -7.0 * zobouks
            zpsihs =  zpsims
            zpsils =  zpsims
            !                                                          ! case of unstable atmospheric conditions
            zobouku = MIN(  0.e0, -100.0 * 10.0 * ( zdtetar + 2.2e-3 * ztvmoyr )  )
            zxins   = ( 1. - 16. * zobouku )**0.25
            zlxins  = LOG( ( 1. + zxins * zxins ) / 2. )
            zpsimu  = 2. * LOG( ( 1 + zxins ) * 0.5 )  + zlxins - 2. * ATAN( zxins ) + rpi * 0.5
            zpsihu  = 2. * zlxins
            zpsilu  = zpsihu
            !                                                          ! intermediate values
            zstab   = MAX( 0.e0, SIGN( 1.e0, zdteta ) )
            zpsim   = zstab * zpsimu + ( 1.0 - zstab ) * zpsims
            zpsih   = zstab * zpsihu + ( 1.0 - zstab ) * zpsihs
            zpsil   = zpsih
            
            zvatmg         = MAX( 0.032 * 1.5e-3 * sf(jp_wndm)%fnow(ji,jj,1) * sf(jp_wndm)%fnow(ji,jj,1) / grav, eps20 )
            zcmn           = vkarmn / LOG ( 10. / zvatmg )
            zchn           = 0.0327 * zcmn
            zcln           = 0.0346 * zcmn
            zcmcmn         = 1. / ( 1. - zcmn * zpsim / vkarmn )
            ! sometimes the ratio zchn * zpsih / ( vkarmn * zcmn ) is too close to 1 and zchcm becomes very very big
            zcmax = 0.1               ! choice for maximum value of the heat transfer coefficient, guided by my intuition
            zrmax = 1 - 3.e-4 / zcmax ! maximum value of the ratio
            zchcm = zcmcmn / ( 1. - MIN ( zchn * zpsih / ( vkarmn * zcmn ) , zrmax ) )
            zclcm          = zchcm
            !                                                          ! transfert coef. (Large and Pond 1981,1982)
            zcsho          = zchn * zchcm                                
            zcleo          = zcln * zclcm 

            zrhova         = zrhoa * sf(jp_wndm)%fnow(ji,jj,1)

            ! sensible heat flux
            zqsb(ji,jj) = zrhova * zcsho * 1004.0  * ( zsst - ztatm )  
         
            ! latent heat flux (bounded by zero)
            zqla(ji,jj) = MAX(  0.e0, zrhova * zcleo * 2.5e+06 * ( zqsato - zqatm )  )
            !               
         END DO
      END DO
      
      ! ----------------------------------------------------------------------------- !
      !     III    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      zcevap = rcp /  cevap    ! convert zqla ==> evap (Kg/m2/s) ==> m/s ==> W/m2
      zcprec = rcp /  rday     ! convert prec ( mm/day ==> m/s)  ==> W/m2

!CDIR COLLAPSE
      emp(:,:) = zqla(:,:) / cevap                                        &   ! freshwater flux
         &     - sf(jp_prec)%fnow(:,:,1) / rday * tmask(:,:,1)
      !
!CDIR COLLAPSE
      qns(:,:) = zqlw(:,:) - zqsb(:,:) - zqla(:,:)                        &   ! Downward Non Solar flux
         &     - zqla(:,:)             * pst(:,:) * zcevap                &   ! remove evap.   heat content at SST in Celcius
         &     + sf(jp_prec)%fnow(:,:,1) * sf(jp_tair)%fnow(:,:,1) * zcprec   ! add    precip. heat content at Tair in Celcius
      qns(:,:) = qns(:,:) * tmask(:,:,1)
#if defined key_lim3
      qns_oce(:,:) = zqlw(:,:) - zqsb(:,:) - zqla(:,:)
      qsr_oce(:,:) = qsr(:,:)
#endif
      ! NB: if sea-ice model, the snow precip are computed and the associated heat is added to qns (see blk_ice_clio)

      IF ( nn_ice == 0 ) THEN
         CALL iom_put( "qlw_oce" ,   zqlw )                 ! output downward longwave  heat over the ocean
         CALL iom_put( "qsb_oce" , - zqsb )                 ! output downward sensible  heat over the ocean
         CALL iom_put( "qla_oce" , - zqla )                 ! output downward latent    heat over the ocean
         CALL iom_put( "qemp_oce",   qns-zqlw+zqsb+zqla )   ! output downward heat content of E-P over the ocean
         CALL iom_put( "qns_oce" ,   qns  )                 ! output downward non solar heat over the ocean
         CALL iom_put( "qsr_oce" ,   qsr  )                 ! output downward solar heat over the ocean
         CALL iom_put( "qt_oce"  ,   qns+qsr )              ! output total downward heat over the ocean
      ENDIF

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=zqsb , clinfo1=' blk_oce_clio: zqsb   : ', tab2d_2=zqlw , clinfo2=' zqlw  : ')
         CALL prt_ctl(tab2d_1=zqla , clinfo1=' blk_oce_clio: zqla   : ', tab2d_2=qsr  , clinfo2=' qsr   : ')
         CALL prt_ctl(tab2d_1=pst  , clinfo1=' blk_oce_clio: pst    : ', tab2d_2=emp  , clinfo2=' emp   : ')
         CALL prt_ctl(tab2d_1=utau , clinfo1=' blk_oce_clio: utau   : ', mask1=umask,   &
            &         tab2d_2=vtau , clinfo2=' vtau : ', mask2=vmask )
      ENDIF

      CALL wrk_dealloc( jpi,jpj, zqlw, zqla, zqsb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('blk_oce_clio')
      !
   END SUBROUTINE blk_oce_clio

# if defined key_lim2 || defined key_lim3
   SUBROUTINE blk_ice_clio_tau
      !!---------------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_clio_tau  ***
      !!                 
      !!  ** Purpose :   Computation momentum flux at the ice-atm interface  
      !!         
      !!  ** Method  :   Read utau from a forcing file. Rearrange if C-grid
      !!
      !!----------------------------------------------------------------------
      REAL(wp) ::   zcoef
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_ice_clio_tau')

      SELECT CASE( cp_ice_msh )

      CASE( 'C' )                          ! C-grid ice dynamics

         zcoef  = cai / cao                         ! Change from air-sea stress to air-ice stress
         utau_ice(:,:) = zcoef * utau(:,:)
         vtau_ice(:,:) = zcoef * vtau(:,:)

      CASE( 'I' )                          ! I-grid ice dynamics:  I-point (i.e. F-point lower-left corner)

         zcoef  = 0.5_wp * cai / cao                ! Change from air-sea stress to air-ice stress
         DO jj = 2, jpj         ! stress from ocean U- and V-points to ice U,V point
            DO ji = 2, jpi   ! I-grid : no vector opt.
               utau_ice(ji,jj) = zcoef * ( utau(ji-1,jj  ) + utau(ji-1,jj-1) )
               vtau_ice(ji,jj) = zcoef * ( vtau(ji  ,jj-1) + vtau(ji-1,jj-1) )
            END DO
         END DO

         CALL lbc_lnk( utau_ice(:,:), 'I', -1. )   ;   CALL lbc_lnk( vtau_ice(:,:), 'I', -1. )   ! I-point

      END SELECT

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=utau_ice , clinfo1=' blk_ice_clio: utau_ice : ', tab2d_2=vtau_ice , clinfo2=' vtau_ice : ')
      ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('blk_ice_clio_tau')

   END SUBROUTINE blk_ice_clio_tau
#endif

# if defined key_lim2 || defined key_lim3
   SUBROUTINE blk_ice_clio_flx(  ptsu , palb_cs, palb_os, palb )
      !!---------------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_clio_flx ***
      !!                 
      !!  ** Purpose :   Computation of the heat fluxes at ocean and snow/ice
      !!       surface the solar heat at ocean and snow/ice surfaces and the 
      !!       sensitivity of total heat fluxes to the SST variations
      !!         
      !!  ** Method  :   The flux of heat at the ice and ocean surfaces are derived
      !!       from semi-empirical ( or bulk ) formulae which relate the flux to 
      !!       the properties of the surface and of the lower atmosphere. Here, we
      !!       follow the work of Oberhuber, 1988   
      !!
      !!  ** Action  :   call albedo_oce/albedo_ice to compute ocean/ice albedo 
      !!               - snow precipitation
      !!               - solar flux at the ocean and ice surfaces
      !!               - the long-wave radiation for the ocean and sea/ice
      !!               - turbulent heat fluxes over water and ice
      !!               - evaporation over water
      !!               - total heat fluxes sensitivity over ice (dQ/dT)
      !!               - latent heat flux sensitivity over ice (dQla/dT)
      !!               - qns  :  modified the non solar heat flux over the ocean
      !!                         to take into account solid precip latent heat flux
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:)   ::   ptsu      ! ice surface temperature                   [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:)   ::   palb_cs  ! ice albedo (clear    sky) (alb_ice_cs)         [-]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:)   ::   palb_os  ! ice albedo (overcast sky) (alb_ice_os)         [-]
      REAL(wp), INTENT(  out), DIMENSION(:,:,:)   ::   palb     ! ice albedo (actual value)                      [-]
      !!
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      !!
      REAL(wp) ::   zmt1, zmt2, zmt3, ztatm3                    ! temporary scalars
      REAL(wp) ::   ztaevbk, zind1, zind2, zind3, ztamr         !    -         -
      REAL(wp) ::   zesi, zqsati, zdesidt                       !    -         -
      REAL(wp) ::   zdqla, zcldeff, zev, zes, zpatm, zrhova     !    -         -
      REAL(wp) ::   zcshi, zclei, zrhovaclei, zrhovacshi        !    -         -
      REAL(wp) ::   ztice3, zticemb, zticemb2, zdqlw, zdqsb     !    -         -
      REAL(wp) ::   z1_lsub                                     !    -         -
      !!
      REAL(wp), DIMENSION(:,:)  , POINTER ::   ztatm   ! Tair in Kelvin
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zqatm   ! specific humidity
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zevsqr  ! vapour pressure square-root
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zrhoa   ! air density
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_qlw, z_qsb
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zevap, zsnw
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_ice_clio_flx')
      !
      CALL wrk_alloc( jpi,jpj, ztatm, zqatm, zevsqr, zrhoa )
      CALL wrk_alloc( jpi,jpj, jpl, z_qlw, z_qsb )

      zpatm = 101000.                        ! atmospheric pressure  (assumed constant  here)
      !--------------------------------------------------------------------------------
      !  Determine cloud optical depths as a function of latitude (Chou et al., 1981).
      !  and the correction factor for taking into account  the effect of clouds 
      !--------------------------------------------------------------------------------

!CDIR NOVERRCHK
!CDIR COLLAPSE
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            ztatm (ji,jj) = sf(jp_tair)%fnow(ji,jj,1)                ! air temperature in Kelvins 
      
            zrhoa(ji,jj) = zpatm / ( 287.04 * ztatm(ji,jj) )         ! air density (equation of state for dry air) 
      
            ztamr = ztatm(ji,jj) - rtt                               ! Saturation water vapour
            zmt1  = SIGN( 17.269,  ztamr )
            zmt2  = SIGN( 21.875,  ztamr )
            zmt3  = SIGN( 28.200, -ztamr )
            zes   = 611.0 * EXP(  ABS( ztamr ) * MIN ( zmt1, zmt2 )   &
               &                / ( ztatm(ji,jj) - 35.86  + MAX( 0.e0, zmt3 ) )  )

            zev = sf(jp_humi)%fnow(ji,jj,1) * zes                    ! vapour pressure  
            zevsqr(ji,jj) = SQRT( zev * 0.01 )                       ! square-root of vapour pressure
            zqatm(ji,jj) = 0.622 * zev / ( zpatm - 0.378 * zev )     ! specific humidity 

            !----------------------------------------------------
            !   Computation of snow precipitation (Ledley, 1985) |
            !----------------------------------------------------
            zmt1  =   253.0 - ztatm(ji,jj)            ;   zind1 = MAX( 0.e0, SIGN( 1.e0, zmt1 ) )
            zmt2  = ( 272.0 - ztatm(ji,jj) ) / 38.0   ;   zind2 = MAX( 0.e0, SIGN( 1.e0, zmt2 ) )
            zmt3  = ( 281.0 - ztatm(ji,jj) ) / 18.0   ;   zind3 = MAX( 0.e0, SIGN( 1.e0, zmt3 ) )
            sprecip(ji,jj) = sf(jp_prec)%fnow(ji,jj,1) / rday   &      ! rday = converte mm/day to kg/m2/s
               &         * (          zind1      &                   ! solid  (snow) precipitation [kg/m2/s]
               &            + ( 1.0 - zind1 ) * (          zind2   * ( 0.5 + zmt2 )   &
               &                                 + ( 1.0 - zind2 ) *  zind3 * zmt3  )   ) 

            !----------------------------------------------------!
            !  fraction of net penetrative shortwave radiation   !
            !----------------------------------------------------!
            ! fraction of qsr_ice which is NOT absorbed in the thin surface layer
            ! and thus which penetrates inside the ice cover ( Maykut and Untersteiner, 1971 ; Elbert anbd Curry, 1993 )
            fr1_i0(ji,jj) = 0.18  * ( 1.e0 - sf(jp_ccov)%fnow(ji,jj,1) ) + 0.35 * sf(jp_ccov)%fnow(ji,jj,1) 
            fr2_i0(ji,jj) = 0.82  * ( 1.e0 - sf(jp_ccov)%fnow(ji,jj,1) ) + 0.65 * sf(jp_ccov)%fnow(ji,jj,1)
         END DO
      END DO
      CALL iom_put( 'snowpre', sprecip )   ! Snow precipitation 
      
      !-----------------------------------------------------------!
      !  snow/ice Shortwave radiation   (abedo already computed)  !
      !-----------------------------------------------------------!
      CALL blk_clio_qsr_ice( palb_cs, palb_os, qsr_ice )
      
      DO jl = 1, jpl
         palb(:,:,jl) = ( palb_cs(:,:,jl) * ( 1.e0 - sf(jp_ccov)%fnow(:,:,1) )   &
            &         +   palb_os(:,:,jl) * sf(jp_ccov)%fnow(:,:,1) )
      END DO

      !                                     ! ========================== !
      DO jl = 1, jpl                       !  Loop over ice categories  !
         !                                  ! ========================== !
!CDIR NOVERRCHK
!CDIR COLLAPSE
         DO jj = 1 , jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               !-------------------------------------------!
               !  long-wave radiation over ice categories  !  ( Berliand 1952 ; all latitudes )
               !-------------------------------------------!
               ztatm3  = ztatm(ji,jj) * ztatm(ji,jj) * ztatm(ji,jj)
               zcldeff = 1.0 - sbudyko(ji,jj) * sf(jp_ccov)%fnow(ji,jj,1) * sf(jp_ccov)%fnow(ji,jj,1)    
               ztaevbk = ztatm3 * ztatm(ji,jj) * zcldeff * ( 0.39 - 0.05 * zevsqr(ji,jj) ) 
               !
               z_qlw(ji,jj,jl) = - emic * stefan * ( ztaevbk + 4. * ztatm3 * ( ptsu(ji,jj,jl) - ztatm(ji,jj) ) ) 

               !----------------------------------------
               !  Turbulent heat fluxes over snow/ice     ( Latent and sensible ) 
               !----------------------------------------        

               ! vapour pressure at saturation of ice (tmask to avoid overflow in the exponential)
               zesi =  611.0 * EXP( 21.8745587 * tmask(ji,jj,1) * ( ptsu(ji,jj,jl) - rtt )/ ( ptsu(ji,jj,jl) - 7.66 ) )
               ! humidity close to the ice surface (at saturation)
               zqsati   = ( 0.622 * zesi ) / ( zpatm - 0.378 * zesi )
               
               !  computation of intermediate values
               zticemb  = ptsu(ji,jj,jl) - 7.66
               zticemb2 = zticemb * zticemb  
               ztice3   = ptsu(ji,jj,jl) * ptsu(ji,jj,jl) * ptsu(ji,jj,jl)
               zdesidt  = zesi * ( 9.5 * LOG( 10.0 ) * ( rtt - 7.66 )  / zticemb2 )
               
               !  Transfer cofficients assumed to be constant (Parkinson 1979 ; Maykut 1982)
               zcshi    = 1.75e-03
               zclei    = zcshi
               
               !  sensible and latent fluxes over ice
               zrhova     = zrhoa(ji,jj) * sf(jp_wndm)%fnow(ji,jj,1)      ! computation of intermediate values
               zrhovaclei = zrhova * zcshi * 2.834e+06
               zrhovacshi = zrhova * zclei * 1004.0
            
               !  sensible heat flux
               z_qsb(ji,jj,jl) = zrhovacshi * ( ptsu(ji,jj,jl) - ztatm(ji,jj) )
            
               !  latent heat flux 
               qla_ice(ji,jj,jl) = MAX(  0.e0, zrhovaclei * ( zqsati - zqatm(ji,jj) )  )
              
               !  sensitivity of non solar fluxes (dQ/dT) (long-wave, sensible and latent fluxes)
               zdqlw = 4.0 * emic * stefan * ztice3
               zdqsb = zrhovacshi
               zdqla = zrhovaclei * ( zdesidt * ( zqsati * zqsati / ( zesi * zesi ) ) * ( zpatm / 0.622 ) )   
               !
               dqla_ice(ji,jj,jl) = zdqla                           ! latent flux sensitivity
               dqns_ice(ji,jj,jl) = -( zdqlw + zdqsb + zdqla )      !  total non solar sensitivity
            END DO
            !
         END DO
         !
      END DO
      !
      ! ----------------------------------------------------------------------------- !
      !    Total FLUXES                                                               !
      ! ----------------------------------------------------------------------------- !
      !
!CDIR COLLAPSE
      qns_ice(:,:,:) = z_qlw (:,:,:) - z_qsb (:,:,:) - qla_ice (:,:,:)      ! Downward Non Solar flux
!CDIR COLLAPSE
      tprecip(:,:)   = sf(jp_prec)%fnow(:,:,1) / rday                     ! total precipitation [kg/m2/s]
      !
      ! ----------------------------------------------------------------------------- !
      !    Correct the OCEAN non solar flux with the existence of solid precipitation !
      ! ---------------=====--------------------------------------------------------- !
!CDIR COLLAPSE
      qns(:,:) = qns(:,:)                                                           &   ! update the non-solar heat flux with:
         &     - sprecip(:,:) * lfus                                                  &   ! remove melting solid precip
         &     + sprecip(:,:) * MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow - rt0 ) * cpic &   ! add solid P at least below melting
         &     - sprecip(:,:) * sf(jp_tair)%fnow(:,:,1)                        * rcp      ! remove solid precip. at Tair

#if defined key_lim3
      ! ----------------------------------------------------------------------------- !
      !    Distribute evapo, precip & associated heat over ice and ocean
      ! ---------------=====--------------------------------------------------------- !
      CALL wrk_alloc( jpi,jpj, zevap, zsnw ) 

      ! --- evaporation --- !
      z1_lsub = 1._wp / Lsub
      evap_ice (:,:,:) = qla_ice (:,:,:) * z1_lsub ! sublimation
      devap_ice(:,:,:) = dqla_ice(:,:,:) * z1_lsub
      zevap    (:,:)   = emp(:,:) + tprecip(:,:)   ! evaporation over ocean

      ! --- evaporation minus precipitation --- !
      zsnw(:,:) = 0._wp
      CALL lim_thd_snwblow( pfrld, zsnw )          ! snow redistribution by wind
      emp_oce(:,:) = pfrld(:,:) * zevap(:,:) - ( tprecip(:,:) - sprecip(:,:) ) - sprecip(:,:) * ( 1._wp - zsnw )
      emp_ice(:,:) = SUM( a_i_b(:,:,:) * evap_ice(:,:,:), dim=3 ) - sprecip(:,:) * zsnw
      emp_tot(:,:) = emp_oce(:,:) + emp_ice(:,:)

      ! --- heat flux associated with emp --- !
      qemp_oce(:,:) = - pfrld(:,:) * zevap(:,:) * sst_m(:,:) * rcp                               & ! evap
         &          + ( tprecip(:,:) - sprecip(:,:) ) * ( sf(jp_tair)%fnow(:,:,1) - rt0 ) * rcp  & ! liquid precip
         &          +   sprecip(:,:) * ( 1._wp - zsnw ) *                                        & ! solid precip
         &              ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow ) - rt0 ) * cpic * tmask(:,:,1) - lfus )
      qemp_ice(:,:) =   sprecip(:,:) * zsnw *                                                    & ! solid precip (only)
         &              ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow ) - rt0 ) * cpic * tmask(:,:,1) - lfus )

      ! --- total solar and non solar fluxes --- !
      qns_tot(:,:) = pfrld(:,:) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 ) + qemp_ice(:,:) + qemp_oce(:,:)
      qsr_tot(:,:) = pfrld(:,:) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      ! --- heat content of precip over ice in J/m3 (to be used in 1D-thermo) --- !
      qprec_ice(:,:) = rhosn * ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow ) - rt0 ) * cpic * tmask(:,:,1) - lfus )

      ! --- heat content of evap over ice in W/m2 (to be used in 1D-thermo) --- !
      DO jl = 1, jpl
         qevap_ice(:,:,jl) = 0._wp ! should be -evap_ice(:,:,jl)*( ( Tice - rt0 ) * cpic * tmask(:,:,1) - lfus )
                                   ! but then qemp_ice should also include sublimation 
      END DO

      CALL wrk_dealloc( jpi,jpj, zevap, zsnw ) 
#endif

!!gm : not necessary as all input data are lbc_lnk...
      CALL lbc_lnk( fr1_i0  (:,:) , 'T', 1. )
      CALL lbc_lnk( fr2_i0  (:,:) , 'T', 1. )
      DO jl = 1, jpl
         CALL lbc_lnk( qns_ice (:,:,jl) , 'T', 1. )
         CALL lbc_lnk( dqns_ice(:,:,jl) , 'T', 1. )
         CALL lbc_lnk( qla_ice (:,:,jl) , 'T', 1. )
         CALL lbc_lnk( dqla_ice(:,:,jl) , 'T', 1. )
      END DO

!!gm : mask is not required on forcing
      DO jl = 1, jpl
         qns_ice (:,:,jl) = qns_ice (:,:,jl) * tmask(:,:,1)
         qla_ice (:,:,jl) = qla_ice (:,:,jl) * tmask(:,:,1)
         dqns_ice(:,:,jl) = dqns_ice(:,:,jl) * tmask(:,:,1)
         dqla_ice(:,:,jl) = dqla_ice(:,:,jl) * tmask(:,:,1)
      END DO

      CALL wrk_dealloc( jpi,jpj, ztatm, zqatm, zevsqr, zrhoa )
      CALL wrk_dealloc( jpi,jpj, jpl  , z_qlw, z_qsb )

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=z_qsb  , clinfo1=' blk_ice_clio: z_qsb  : ', tab3d_2=z_qlw  , clinfo2=' z_qlw  : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=qla_ice  , clinfo1=' blk_ice_clio: z_qla  : ', tab3d_2=qsr_ice  , clinfo2=' qsr_ice  : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=dqns_ice , clinfo1=' blk_ice_clio: dqns_ice : ', tab3d_2=qns_ice  , clinfo2=' qns_ice  : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=dqla_ice , clinfo1=' blk_ice_clio: dqla_ice : ', tab3d_2=ptsu    , clinfo2=' ptsu    : ', kdim=jpl)
         CALL prt_ctl(tab2d_1=tprecip  , clinfo1=' blk_ice_clio: tprecip  : ', tab2d_2=sprecip  , clinfo2=' sprecip  : ')
      ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('blk_ice_clio_flx')
      !
   END SUBROUTINE blk_ice_clio_flx

#endif

   SUBROUTINE blk_clio_qsr_oce( pqsr_oce )
      !!---------------------------------------------------------------------------
      !!                     ***  ROUTINE blk_clio_qsr_oce  ***
      !!                 
      !!  ** Purpose :   Computation of the shortwave radiation at the ocean and the
      !!               snow/ice surfaces. 
      !!         
      !!  ** Method  : - computed qsr from the cloud cover for both ice and ocean 
      !!               - also initialise sbudyko and stauc once for all 
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj)     ::   pqsr_oce    ! shortwave radiation  over the ocean
      !!
      INTEGER, PARAMETER  ::   jp24 = 24   ! sampling of the daylight period (sunrise to sunset) into 24 equal parts
      !!      
      INTEGER  ::   ji, jj, jt    ! dummy loop indices 
      INTEGER  ::   indaet            !  = -1, 0, 1 for odd, normal and leap years resp.
      INTEGER  ::   iday              ! integer part of day
      INTEGER  ::   indxb, indxc      ! index for cloud depth coefficient

      REAL(wp)  ::   zalat , zclat, zcmue, zcmue2    ! local scalars 
      REAL(wp)  ::   zmt1, zmt2, zmt3                ! 
      REAL(wp)  ::   zdecl, zsdecl , zcdecl          ! 
      REAL(wp)  ::   za_oce, ztamr                   !

      REAL(wp) ::   zdl, zlha                        ! local scalars
      REAL(wp) ::   zlmunoon, zcldcor, zdaycor       !   
      REAL(wp) ::   zxday, zdist, zcoef, zcoef1      !
      REAL(wp) ::   zes
      
      REAL(wp), DIMENSION(:,:), POINTER ::   zev          ! vapour pressure
      REAL(wp), DIMENSION(:,:), POINTER ::   zdlha, zlsrise, zlsset     ! 2D workspace
      REAL(wp), DIMENSION(:,:), POINTER ::   zps, zpc   ! sine (cosine) of latitude per sine (cosine) of solar declination 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_clio_qsr_oce')
      !
      CALL wrk_alloc( jpi,jpj, zev, zdlha, zlsrise, zlsset, zps, zpc )

      IF( lbulk_init ) THEN             !   Initilization at first time step only
         rdtbs2 = nn_fsbc * rdt * 0.5
         ! cloud optical depths as a function of latitude (Chou et al., 1981).
         ! and the correction factor for taking into account  the effect of clouds 
         DO jj = 1, jpj
            DO ji = 1 , jpi
               zalat          = ( 90.e0 - ABS( gphit(ji,jj) ) ) /  5.e0
               zclat          = ( 95.e0 -      gphit(ji,jj)   ) / 10.e0
               indxb          = 1 + INT( zalat )
               indxc          = 1 + INT( zclat )
               zdl            = zclat - INT( zclat )
               !  correction factor to account for the effect of clouds
               sbudyko(ji,jj) = budyko(indxb)
               stauc  (ji,jj) = ( 1.e0 - zdl ) * tauco( indxc ) + zdl * tauco( indxc + 1 )
            END DO
         END DO
         lbulk_init = .FALSE.
      ENDIF


      ! Saturated water vapour and vapour pressure
      ! ------------------------------------------
!CDIR NOVERRCHK
!CDIR COLLAPSE
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            ztamr = sf(jp_tair)%fnow(ji,jj,1) - rtt
            zmt1  = SIGN( 17.269,  ztamr )
            zmt2  = SIGN( 21.875,  ztamr )
            zmt3  = SIGN( 28.200, -ztamr )
            zes = 611.0 * EXP(  ABS( ztamr ) * MIN ( zmt1, zmt2 )   &              ! Saturation water vapour
               &                     / ( sf(jp_tair)%fnow(ji,jj,1) - 35.86  + MAX( 0.e0, zmt3 ) )  )
            zev(ji,jj) = sf(jp_humi)%fnow(ji,jj,1) * zes * 1.0e-05                 ! vapour pressure  
         END DO
      END DO

      !-----------------------------------!
      !  Computation of solar irradiance  !
      !-----------------------------------!
!!gm : hard coded  leap year ???
      indaet   = 1                                    ! = -1, 0, 1 for odd, normal and leap years resp.
      zxday = nday_year + rdtbs2 / rday               ! day of the year at which the fluxes are calculated
      iday  = INT( zxday )                            ! (centred at the middle of the ice time step)
      CALL flx_blk_declin( indaet, iday, zdecl )      ! solar declination of the current day
      zsdecl = SIN( zdecl * rad )                     ! its sine
      zcdecl = COS( zdecl * rad )                     ! its cosine


      !  correction factor added for computation of shortwave flux to take into account the variation of
      !  the distance between the sun and the earth during the year (Oberhuber 1988)
      zdist    = zxday * 2. * rpi / REAL(nyear_len(1), wp)
      zdaycor  = 1.0 + 0.0013 * SIN( zdist ) + 0.0342 * COS( zdist )

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            !  product of sine (cosine) of latitude and sine (cosine) of solar declination
            zps(ji,jj) = SIN( gphit(ji,jj) * rad ) * zsdecl
            zpc(ji,jj) = COS( gphit(ji,jj) * rad ) * zcdecl
            !  computation of the both local time of sunrise and sunset
            zlsrise(ji,jj) = ACOS( - SIGN( 1.e0, zps(ji,jj) )    &
               &                   * MIN(  1.e0, SIGN( 1.e0, zps(ji,jj) ) * ( zps(ji,jj) / zpc(ji,jj) )  )   )
            zlsset (ji,jj) = - zlsrise(ji,jj)
            !  dividing the solar day into jp24 segments of length zdlha
            zdlha  (ji,jj) = ( zlsrise(ji,jj) - zlsset(ji,jj) ) / REAL( jp24, wp )
         END DO
      END DO


      !---------------------------------------------!
      !  shortwave radiation absorbed by the ocean  !
      !---------------------------------------------!
      pqsr_oce(:,:)   = 0.e0      ! set ocean qsr to zero      

      ! compute and sum ocean qsr over the daylight (i.e. between sunrise and sunset)
!CDIR NOVERRCHK   
      DO jt = 1, jp24
         zcoef = FLOAT( jt ) - 0.5
!CDIR NOVERRCHK     
!CDIR COLLAPSE
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zlha = COS(  zlsrise(ji,jj) - zcoef * zdlha(ji,jj)  )                  ! local hour angle
               zcmue              = MAX( 0.e0 ,   zps(ji,jj) + zpc(ji,jj) * zlha  )   ! cos of local solar altitude
               zcmue2             = 1368.0 * zcmue * zcmue

               ! ocean albedo depending on the cloud cover (Payne, 1972)
               za_oce     = ( 1.0 - sf(jp_ccov)%fnow(ji,jj,1) ) * 0.05 / ( 1.1 * zcmue**1.4 + 0.15 )   &   ! clear sky
                  &       +         sf(jp_ccov)%fnow(ji,jj,1)   * 0.06                                     ! overcast

                  ! solar heat flux absorbed by the ocean (Zillman, 1972)
               pqsr_oce(ji,jj) = pqsr_oce(ji,jj)                                         &
                  &            + ( 1.0 - za_oce ) * zdlha(ji,jj) * zcmue2                &
                  &            / ( ( zcmue + 2.7 ) * zev(ji,jj) + 1.085 * zcmue +  0.10 )
            END DO
         END DO
      END DO
      ! Taking into account the ellipsity of the earth orbit, the clouds AND masked if sea-ice cover > 0%
      zcoef1 = srgamma * zdaycor / ( 2. * rpi )
!CDIR COLLAPSE
      DO jj = 1, jpj
         DO ji = 1, jpi
            zlmunoon = ASIN( zps(ji,jj) + zpc(ji,jj) ) / rad                         ! local noon solar altitude
            zcldcor  = MIN(  1.e0, ( 1.e0 - 0.62 * sf(jp_ccov)%fnow(ji,jj,1)   &     ! cloud correction (Reed 1977)
               &                          + 0.0019 * zlmunoon )                 )
            pqsr_oce(ji,jj) = zcoef1 * zcldcor * pqsr_oce(ji,jj) * tmask(ji,jj,1)    ! and zcoef1: ellipsity
         END DO
      END DO

      CALL wrk_dealloc( jpi,jpj, zev, zdlha, zlsrise, zlsset, zps, zpc )
      !
      IF( nn_timing == 1 )  CALL timing_stop('blk_clio_qsr_oce')
      !
   END SUBROUTINE blk_clio_qsr_oce


   SUBROUTINE blk_clio_qsr_ice( pa_ice_cs, pa_ice_os, pqsr_ice )
      !!---------------------------------------------------------------------------
      !!                     ***  ROUTINE blk_clio_qsr_ice  ***
      !!                 
      !!  ** Purpose :   Computation of the shortwave radiation at the ocean and the
      !!               snow/ice surfaces. 
      !!         
      !!  ** Method  : - computed qsr from the cloud cover for both ice and ocean 
      !!               - also initialise sbudyko and stauc once for all 
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:) ::   pa_ice_cs   ! albedo of ice under clear sky
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:) ::   pa_ice_os   ! albedo of ice under overcast sky
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pqsr_ice    ! shortwave radiation over the ice/snow
      !!
      INTEGER, PARAMETER  ::   jp24 = 24   ! sampling of the daylight period (sunrise to sunset) into 24 equal parts
      !!
      INTEGER  ::   ji, jj, jl, jt    ! dummy loop indices
      INTEGER  ::   ijpl              ! number of ice categories (3rd dim of pqsr_ice)
      INTEGER  ::   indaet            !  = -1, 0, 1 for odd, normal and leap years resp.
      INTEGER  ::   iday              ! integer part of day
      !!
      REAL(wp) ::   zcmue, zcmue2, ztamr          ! temporary scalars 
      REAL(wp) ::   zmt1, zmt2, zmt3              !    -         -
      REAL(wp) ::   zdecl, zsdecl, zcdecl         !    -         -
      REAL(wp) ::   zlha, zdaycor, zes            !    -         -
      REAL(wp) ::   zxday, zdist, zcoef, zcoef1   !    -         -
      REAL(wp) ::   zqsr_ice_cs, zqsr_ice_os      !    -         -

      REAL(wp), DIMENSION(:,:), POINTER ::   zev                      ! vapour pressure
      REAL(wp), DIMENSION(:,:), POINTER ::   zdlha, zlsrise, zlsset   ! 2D workspace
      REAL(wp), DIMENSION(:,:), POINTER ::   zps, zpc   ! sine (cosine) of latitude per sine (cosine) of solar declination 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_clio_qsr_ice')
      !
      CALL wrk_alloc( jpi,jpj, zev, zdlha, zlsrise, zlsset, zps, zpc )

      ijpl = SIZE(pqsr_ice, 3 )      ! number of ice categories
      
      ! Saturated water vapour and vapour pressure
      ! ------------------------------------------
!CDIR NOVERRCHK
!CDIR COLLAPSE
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi           
            ztamr = sf(jp_tair)%fnow(ji,jj,1) - rtt           
            zmt1  = SIGN( 17.269,  ztamr )
            zmt2  = SIGN( 21.875,  ztamr )
            zmt3  = SIGN( 28.200, -ztamr )
            zes = 611.0 * EXP(  ABS( ztamr ) * MIN ( zmt1, zmt2 )   &              ! Saturation water vapour
               &                     / ( sf(jp_tair)%fnow(ji,jj,1) - 35.86  + MAX( 0.e0, zmt3 ) )  )
            zev(ji,jj) = sf(jp_humi)%fnow(ji,jj,1) * zes * 1.0e-05                 ! vapour pressure  
         END DO
      END DO

      !-----------------------------------!
      !  Computation of solar irradiance  !
      !-----------------------------------!
!!gm : hard coded  leap year ???
      indaet   = 1                                    ! = -1, 0, 1 for odd, normal and leap years resp.
      zxday = nday_year + rdtbs2 / rday               ! day of the year at which the fluxes are calculated
      iday  = INT( zxday )                            ! (centred at the middle of the ice time step)
      CALL flx_blk_declin( indaet, iday, zdecl )      ! solar declination of the current day
      zsdecl = SIN( zdecl * rad )                     ! its sine
      zcdecl = COS( zdecl * rad )                     ! its cosine

      
      !  correction factor added for computation of shortwave flux to take into account the variation of
      !  the distance between the sun and the earth during the year (Oberhuber 1988)
      zdist    = zxday * 2. * rpi / REAL(nyear_len(1), wp)
      zdaycor  = 1.0 + 0.0013 * SIN( zdist ) + 0.0342 * COS( zdist )

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            !  product of sine (cosine) of latitude and sine (cosine) of solar declination
            zps(ji,jj) = SIN( gphit(ji,jj) * rad ) * zsdecl
            zpc(ji,jj) = COS( gphit(ji,jj) * rad ) * zcdecl
            !  computation of the both local time of sunrise and sunset
            zlsrise(ji,jj) = ACOS( - SIGN( 1.e0, zps(ji,jj) )    &
               &                   * MIN(  1.e0, SIGN( 1.e0, zps(ji,jj) ) * ( zps(ji,jj) / zpc(ji,jj) )  )   ) 
            zlsset (ji,jj) = - zlsrise(ji,jj)
            !  dividing the solar day into jp24 segments of length zdlha
            zdlha  (ji,jj) = ( zlsrise(ji,jj) - zlsset(ji,jj) ) / REAL( jp24, wp )
         END DO
      END DO


      !---------------------------------------------!
      !  shortwave radiation absorbed by the ice    !
      !---------------------------------------------!
      ! compute and sum ice qsr over the daylight for each ice categories
      pqsr_ice(:,:,:) = 0.e0
      zcoef1 = zdaycor / ( 2. * rpi )       ! Correction for the ellipsity of the earth orbit
      
      !                    !----------------------------! 
      DO jl = 1, ijpl      !  loop over ice categories  !
         !                 !----------------------------! 
!CDIR NOVERRCHK   
         DO jt = 1, jp24   
            zcoef = FLOAT( jt ) - 0.5
!CDIR NOVERRCHK     
!CDIR COLLAPSE
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  zlha = COS(  zlsrise(ji,jj) - zcoef * zdlha(ji,jj)  )                  ! local hour angle
                  zcmue              = MAX( 0.e0 ,   zps(ji,jj) + zpc(ji,jj) * zlha  )   ! cos of local solar altitude
                  zcmue2             = 1368.0 * zcmue * zcmue
                  
                  !  solar heat flux absorbed by the ice/snow system (Shine and Crane 1984 adapted to high albedo) 
                  zqsr_ice_cs =  ( 1.0 - pa_ice_cs(ji,jj,jl) ) * zdlha(ji,jj) * zcmue2        &   ! clear sky
                     &        / ( ( 1.0 + zcmue ) * zev(ji,jj) + 1.2 * zcmue + 0.0455 )
                  zqsr_ice_os = zdlha(ji,jj) * SQRT( zcmue )                                  &   ! overcast sky
                     &        * ( 53.5 + 1274.5 * zcmue )      * ( 1.0 - 0.996  * pa_ice_os(ji,jj,jl) )    &
                     &        / (  1.0 + 0.139  * stauc(ji,jj) * ( 1.0 - 0.9435 * pa_ice_os(ji,jj,jl) ) )       
             
                  pqsr_ice(ji,jj,jl) = pqsr_ice(ji,jj,jl) + (  ( 1.0 - sf(jp_ccov)%fnow(ji,jj,1) ) * zqsr_ice_cs    &
                     &                                       +         sf(jp_ccov)%fnow(ji,jj,1)   * zqsr_ice_os  )
               END DO
            END DO
         END DO
         !
         ! Correction : Taking into account the ellipsity of the earth orbit
         pqsr_ice(:,:,jl) = pqsr_ice(:,:,jl) * zcoef1 * tmask(:,:,1)
         !
         !                 !--------------------------------! 
      END DO               !  end loop over ice categories  !
      !                    !--------------------------------! 


!!gm  : this should be suppress as input data have been passed through lbc_lnk
      DO jl = 1, ijpl
         CALL lbc_lnk( pqsr_ice(:,:,jl) , 'T', 1. )
      END DO
      !
      CALL wrk_dealloc( jpi,jpj, zev, zdlha, zlsrise, zlsset, zps, zpc )
      !
      IF( nn_timing == 1 )  CALL timing_stop('blk_clio_qsr_ice')
      !
   END SUBROUTINE blk_clio_qsr_ice


   SUBROUTINE flx_blk_declin( ky, kday, pdecl )
      !!---------------------------------------------------------------------------
      !!               ***  ROUTINE flx_blk_declin  ***
      !!          
      !! ** Purpose :   Computation of the solar declination for the day
      !!       
      !! ** Method  :   ???
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in   ) ::   ky      ! = -1, 0, 1 for odd, normal and leap years resp.
      INTEGER , INTENT(in   ) ::   kday    ! day of the year ( kday = 1 on january 1)
      REAL(wp), INTENT(  out) ::   pdecl   ! solar declination
      !!
      REAL(wp) ::   a0  =  0.39507671      ! coefficients for solar declinaison computation
      REAL(wp) ::   a1  = 22.85684301      !     "              ""                 "
      REAL(wp) ::   a2  = -0.38637317      !     "              ""                 "
      REAL(wp) ::   a3  =  0.15096535      !     "              ""                 "
      REAL(wp) ::   a4  = -0.00961411      !     "              ""                 "
      REAL(wp) ::   b1  = -4.29692073      !     "              ""                 "
      REAL(wp) ::   b2  =  0.05702074      !     "              ""                 "
      REAL(wp) ::   b3  = -0.09028607      !     "              ""                 "
      REAL(wp) ::   b4  =  0.00592797
      !!
      REAL(wp) ::   zday   ! corresponding day of type year (cf. ky)
      REAL(wp) ::   zp     ! temporary scalars
      !!---------------------------------------------------------------------
            
      IF    ( ky == 1 )  THEN   ;   zday = REAL( kday, wp ) - 0.5
      ELSEIF( ky == 3 )  THEN   ;   zday = REAL( kday, wp ) - 1.
      ELSE                      ;   zday = REAL( kday, wp )
      ENDIF
      
      zp = rpi * ( 2.0 * zday - 367.0 ) / REAL(nyear_len(1), wp)
      
      pdecl  = a0                                                                      &
         &   + a1 * COS( zp ) + a2 * COS( 2. * zp ) + a3 * COS( 3. * zp ) + a4 * COS( 4. * zp )   &
         &   + b1 * SIN( zp ) + b2 * SIN( 2. * zp ) + b3 * SIN( 3. * zp ) + b4 * SIN( 4. * zp )
      !
   END SUBROUTINE flx_blk_declin

   !!======================================================================
END MODULE sbcblk_clio
