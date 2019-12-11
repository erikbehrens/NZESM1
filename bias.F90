MODULE bias 
   !! Is used by OPA and STEP
   !!======================================================================
   !!                 *** Module bias ***
   !! Code to estimate and apply bias correction.
   !! The bias is in T/S and Pressure. It follows the generalized
   !! bias algorithm presented in Balmaseda et al 2007.
   !!
   !! It can be read from a file offline, estimated online from relaxation
   !! terms or from assimilation increments (this latter estimtd in inner loop)
   !!
   !! The parameters for the time evolution of the bias can be different
   !! from the relaxation terms and for the assim increments. Only the 
   !! parameter for the relaxtion terms are considered here.
   !! 
   !! The offline bias term can contain the seasonal cycle.
   !!
   !! The time evolution of the bias for relaxtion is estimated as followed
   !!   bias_rlx(t+1)=t_rlx_mem*bias_rlx(t)+t_rlx_upd*(t/s)trdmp.
   !!
   !! The total bias in T/S is partion between the correction to T/S only
   !!  (fb_t) and the correction applied to the pressure gradient (fb_p).
   !!  We impose that (fb_p=1.-fb_t). These factors can be different for the
   !!  different terms(fb_t_rxl,fb_t_asm,fb_t_ofl)
   !!
   !!    (t/s)bias =  fb_t_ofl * (t/s)bias_ofl +
   !!                 fb_t_rlx * (t/s)bias_rlx +
   !!                 fb_t_asm * (t/s)bias_asm
   !!    (t/s)bias_p =fb_p_ofl * (t/s)bias_ofl+ 
   !!                 fb_p_rlx * (t/s)bias_rlx_p +
   !!                 fb_p_asm * (t/s)bias_asm_p
   !!    (t/s)bias is applied directely to correct T and S
   !!    (t/s)bias_p is used to compute rhd_pc and gru/v_pc
   !!
   !!  Note: the above is an adhoc /simple way of making the partition
   !!        between bias in pressure and bias in T/S. It would be better
   !!        if the partition is done at the time of estimating the time
   !!        evolution of the bias. That would mean doubling the number of
   !!        3D arrays.
   !!
   !!  New addtion: (not in Balmaseda et al 2007):
   !!  A more physical alternative to the partition of the bias can be
   !!  done in terms of the inertial frequency: when the time scales of 
   !!  adjustment are shorter that >1/f (Equator), the bias correction should
   !!  be in the the pressure term. Otherwise, it can act directly in T/S.
   !!  NOTE that the bias correction in the pressure term here (following
   !!  (Bell et al 2007) is just the "opposite" as the semi-prognostic method 
   !!  in Greatbatch et al 2004. 
   !!  The use of this partition is controlled by ln_inertial=.true.
   !!  
   !!
   !!        2009-03              (M.A. Balmaseda ECMWF)
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   bias_init  : Read in the bias namelist and the bias arrays
   !!   tra_bias   : Apply the bias fields on T/S directly
   !!   dyn_bias   : Compute density correction for dynamic hpg
   !!   bias_opn   : open bias files for restart capabilities
   !!   bias_wrt   : write bias fies "     "          "
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY: &
      & wp
   USE par_oce, ONLY: &
      & jpi, &
      & jpj, &
      & jpk
   USE dom_oce, ONLY: &
      & rdt,          &
      & ln_zps,       &
      & gphit
   USE phycst,  ONLY: &
      & rday,         &
      & rad
   USE oce, ONLY: &
      & tsb, tsn, tsa, &
      & rhop,  &
      & gtsu, gtsv
   USE dynhpg, ONLY:   &
      & ln_dynhpg_imp      
   USE tradmp
   USE dtatsd, ONLY: &
      & ln_tsd_tradmp
   USE in_out_manager, ONLY: &
      & lwp,          &
      & numnam_ref,   &
      & numnam_cfg,   &
      & numond,       &
      & numout,       &
      & lrst_oce,     &
      & nit000
   USE iom
   USE eosbn2
   USE zpshde          ! partial step: hor. derivative (zps_hde routine)
   USE biaspar
   USE fldread         ! read input fields
   USE lbclnk          ! lateral boundary conditions (or mpp link)
   USE asmpar
   USE asminc
   USE lib_mpp, ONLY: &
      & ctl_stop, &
      & ctl_nam

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE   
   PUBLIC bias_init,   &   !: Read in the bias namelist and the bias arrays
      &   tra_bias,    &   !: Estimate/apply bias on traces
      &   dyn_bias,    &   !: " density correction for pressure gradient.
      &   bias_opn,    &
      &   bias_wrt

   !! * Substitutions needed to have the variable fsdept_n
#  include "domzgr_substitute.h90"
   !! * Shared variables
   !! * Private module variables
   REAL(wp), PRIVATE :: &
      & bias_time_unit_asm,   &  !: bias_asm units in s ( per day = 86400 s)   
      & bias_time_unit_rlx,   &  !: bias_rlx units in s ( 1 second) 
      & bias_time_unit_ofl,   &  !: bias_ofl units in s ( 1 second) 
      & t_rlx_mem,            &  !: time param for mem in bias_rlx model
      & t_rlx_upd,            &  !: time param for update in bias_rlx model
                                 !: (pct of incr for computation of bias)
      & t_asm_mem,            &  !: time param for mem in bias_asm model
      & t_asm_upd,            &  !: time param for update in bias_asm model
                                 !: (pct of incr for computation of bias)
      & fb_t_rlx,             &  !: parition of bias in T for rlx bias term
      & fb_t_asm,             &  !: parition of bias in T for asm bias term
      & fb_t_ofl,             &  !: parition of bias in T for ofl bias term
      & fb_p_rlx,             &  !: parition of bias in P for rlx bias term
      & fb_p_asm,             &  !: parition of bias in P for asm bias term
      & fb_p_ofl,             &  !: parition of bias in P for ofl bias term
      & fctamp,               &  !: amplification factor for T if inertial
      & rn_maxlat_bias,       &  !: Max lat for latitudinal ramp
      & rn_minlat_bias,       &  !: Min lat for latitudinal ramp
      & zwgt,                 &  !: weight for IPC
      & ztscale                  !: decay rate for IPC

   LOGICAL,  PRIVATE :: lalloc
   REAL(wp), PRIVATE, DIMENSION(:,:,:), ALLOCATABLE :: &
      & tbias_asm, &       !: Temperature bias field
      & sbias_asm, &       !: Salinity bias field
      & tbias_rlx, &       !: Temperature bias field
      & sbias_rlx, &       !: Salinity bias field
      & tbias_asm_out, &   !: Output temperature bias field
      & sbias_asm_out, &   !: Output salinity bias field
      & tbias_rlx_out, &   !: Output temperature bias field
      & sbias_rlx_out, &   !: Output salinity bias field
      & tbias_p_out,   &   !: Output temperature bias field for P correction
      & sbias_p_out,   &   !: Output salinity bias field for P correction
      & tbias_i_out,   &   !: Output temperature bias field for incremental P correction
      & sbias_i_out,   &   !: Output salinity bias field for incremental P correction
      & tbias_asm_stscale, &   !: Short time scale temperature bias field
      & sbias_asm_stscale, &   !: Short time scale salinity bias field
      & tbias_asm_stscale_out, &   !: Short time scale temperature bias output field
      & sbias_asm_stscale_out   !: Short time scale salinity bias output field

   INTEGER, PRIVATE :: nn_lat_ramp     ! choice of latitude dependent ramp
                                       ! for the pressure correction.
				       ! 1:inertial ramp, 2:linear ramp, else:no ramp
   LOGICAL, PRIVATE :: ln_bsyncro      ! syncronous or assincrous bias correction   
   LOGICAL, PRIVATE :: ln_itdecay      ! evolve bias correction at every time step.  
   LOGICAL, PRIVATE :: ln_incpc        ! incremental pressure correction plus pressure correction
   LOGICAL, PRIVATE :: ln_incpc_only        ! incremental pressure correction only

   REAL(wp), PRIVATE, DIMENSION(:,:), ALLOCATABLE :: fbcoef
   REAL(wp), PRIVATE, DIMENSION(:,:), ALLOCATABLE :: fbcoef_stscale

   INTEGER, PRIVATE ::  &
      & numbias_asm, &    ! netcdf id of bias file from assim
      & numbias_tot, &    ! netcdf id of bias file with total bias
      & nn_bias_itwrt     ! time step for outputting bias pressure corr
         
   CHARACTER(LEN=128), PRIVATE :: &
      & cn_bias_asm,   &  ! name of bias file from assim
      & cn_bias_tot       ! name of bias with total/rlx bias

   ! Structure of input T and S bias offline (file informations, fields read)
   TYPE(FLD), PRIVATE, ALLOCATABLE, DIMENSION(:) ::   sf_tbias_ofl  
   TYPE(FLD), PRIVATE, ALLOCATABLE, DIMENSION(:) ::   sf_sbias_ofl

   TYPE(FLD_N), PRIVATE ::&   ! information about the fields to be read
      &  sn_tbias_ofl, sn_sbias_ofl  
   
CONTAINS

   SUBROUTINE bias_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE bias_init  ***
      !!
      !! ** Purpose : Read in the bias namelist and read in the bias fields.
      !!
      !! ** Method  : Read in the bias namelist and read in the bias fields.
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  08-05  (D. Lea)    Initial version
      !!        !  08-10  (M. Martin) Tidy
      !!        !  09-03  (M. Balmaseda). Generalize to estimate the bias
      !!                                  from relax and offline bias term.
      !!                                  Introduce parameters to control the
      !!                                  model for the bias 
      !!                                  (variables and time evolution)
      !!----------------------------------------------------------------------

      IMPLICIT NONE
      
      !! * Local declarations

      CHARACTER(len=100) ::  cn_dir          ! dir for location ofline bias
      INTEGER            ::  ierror
      INTEGER            ::  ios             ! Local integer output status for namelist read
      REAL(wp)           ::  eft_rlx,  &     ! efolding time (bias memory) in days
         &                   eft_asm,  &     !     "      "
         &                   log2,     &
         &                   lenscl_bias, &  ! lengthscale pressure bias decay between minlat and maxlat.
         &                   minlat_bias, &  ! used in ipc
         &                   maxlat_bias     ! used in ipc
      
      NAMELIST/nambias/ ln_bias, ln_bias_asm, ln_bias_rlx, ln_bias_ofl,   &
         & ln_bias_ts_app, ln_bias_pc_app,                                &
         & fb_t_asm, fb_t_rlx, fb_t_ofl, fb_p_asm, fb_p_rlx, fb_p_ofl,    &
         & eft_rlx, t_rlx_upd, eft_asm, t_asm_upd, nn_lat_ramp,           &
         & bias_time_unit_asm, bias_time_unit_rlx, bias_time_unit_ofl,    &
         & cn_bias_tot, cn_bias_asm, cn_dir, sn_tbias_ofl, sn_sbias_ofl,  &
         & ln_bsyncro, fctamp, rn_maxlat_bias, rn_minlat_bias,            &
         & nn_bias_itwrt, ln_itdecay, ln_incpc,ln_incpc_only, ztscale, zwgt
 

      !-----------------------------------------------------------------------
      ! Read Namelist : bias interface
      !-----------------------------------------------------------------------


      REWIND( numnam_ref )              ! Namelist nambias in reference namelist : Bias pressure correction
      READ  ( numnam_ref, nambias, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambias in reference namelist', lwp )


      ! Set additional default values (note that most values are set in the reference namelist)
      
      IF ( ln_asmiau ) nn_bias_itwrt = nitiaufin
      
      ! ... default values (NB: frequency positive => hours, negative => months)
      !            !   file    ! frequency !  variable  ! time intep !  clim   ! 'yearly' or !
      !            !   name    !  (hours)  !   name     !   (T/F)    !  (T/F)  !  'monthly'  !
      sn_tbias_ofl = FLD_N( 'tbias_ofl'    ,    -1.    ,  'tbias'     ,  .TRUE.   , .FALSE. ,   'yearly', '', '', ''  )
      sn_sbias_ofl = FLD_N( 'sbias_ofl'    ,    -1.    ,  'sbias'     ,  .TRUE.   , .FALSE. ,   'yearly', '', '', ''  )


      REWIND( numnam_cfg )              ! Namelist nambias in configuration namelist : Bias pressure correction
      READ  ( numnam_cfg, nambias, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambias in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nambias )


      IF ( ( .NOT. ln_bias_asm ) .AND. ( .NOT. ln_bias_ofl ) .AND. ( .NOT. ln_bias_rlx ) ) THEN
         ln_bias_ts_app = .FALSE.
         ln_bias_pc_app = .FALSE.
         ln_bias        = .FALSE.
      ENDIF
      
      ! set up decay scales
      log2           = LOG( 2.0_wp )
      t_rlx_mem      = EXP( - log2 * rdt / ( eft_rlx * rday ) )
      t_asm_mem      = EXP( - log2 * bias_time_unit_asm/ ( eft_asm * rday ) )
      
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'bias_init : '
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '  Namelist nambias : '

         WRITE(numout,*) '  Bias switches/options/variables '
         WRITE(numout,*) '     bias main switch               ln_bias        = ',ln_bias
         WRITE(numout,*) '     bias from assim                ln_bias_asm    = ',ln_bias_asm
         WRITE(numout,*) '     bias from relax                ln_bias_rlx    = ',ln_bias_rlx
         WRITE(numout,*) '     bias from offln                ln_bias_ofl    = ',ln_bias_ofl
         WRITE(numout,*) '     bias T and S apply             ln_bias_ts_app = ',ln_bias_ts_app
         WRITE(numout,*) '     bias pressure correctn apply   ln_bias_pc_app = ',ln_bias_pc_app
         WRITE(numout,*) '     bias pressure correctn apply   ln_bias_pc_app = ',ln_bias_pc_app
         WRITE(numout,*) '     lat ramp for bias correction   nn_lat_ramp    = ',nn_lat_ramp
         WRITE(numout,*) '     time step for writing bias fld nn_bias_itwrt  = ',nn_bias_itwrt
	 WRITE(numout,*) '     evolve pcbias at each timestep ln_itdecay     = ',ln_itdecay
	 WRITE(numout,*) '     incremental press. correction  ln_incpc       = ',ln_incpc
	 WRITE(numout,*) '     incremental press. correction only  ln_incpc_only       = ',ln_incpc_only
	 WRITE(numout,*) '     incremental press. correction  zwgt           = ',zwgt
	 WRITE(numout,*) '     incremental press. correction  ztscale        = ',ztscale
         WRITE(numout,*) '  Parameters for parition of bias term '
         WRITE(numout,*) '                                    fb_t_rlx       = ',fb_t_rlx
         WRITE(numout,*) '                                    fb_t_asm       = ',fb_t_asm
         WRITE(numout,*) '                                    fb_t_ofl       = ',fb_t_ofl
         WRITE(numout,*) '                                    fb_p_rlx       = ',fb_p_rlx
         WRITE(numout,*) '                                    fb_p_asm       = ',fb_p_asm
         WRITE(numout,*) '                                    fb_p_ofl       = ',fb_p_ofl
         WRITE(numout,*) '  Parameters for time evolution of bias '
         WRITE(numout,*) '  Rlx   efolding time (mem)         eft_rlx,t_rlx_mem = ', eft_rlx, t_rlx_mem, 1. - log2 * rdt / (eft_rlx * rday)
         WRITE(numout,*) '        uptdate factor              t_rlx_upd         = ',t_rlx_upd
         WRITE(numout,*) '  Asm   efolding time (mem)         eft_asm,t_asm_mem = ', eft_asm, t_asm_mem, 1. - log2 * rdt / (eft_asm * rday)
         WRITE(numout,*) '        uptdate factor              t_asm_upd         = ',t_asm_upd
         WRITE(numout,*) '  Filenames and input structures'
         WRITE(numout,*) '     bias_tot filename              cn_bias_to     = ',cn_bias_tot
         WRITE(numout,*) '     bias_asm filename              cn_bias_asm    = ',cn_bias_asm
         WRITE(numout,*) '     bias_asm time unit (secs)  bias_time_unit_asm = ',bias_time_unit_asm
         WRITE(numout,*) '     structure Tem bias ofl         sn_tbias_ofl   = ',sn_tbias_ofl
         WRITE(numout,*) '     structure Sal bias ofl         sn_sbias_ofl   = ',sn_sbias_ofl
         
         IF ( ( (.NOT. ln_tsd_tradmp) .OR. (.NOT. ln_tradmp) ) .AND. ln_bias_rlx ) &
            &   CALL ctl_stop (' lk_dtatem, lk_dtasal and lk_tradmp need to be true with ln_bias_rlx' )

         IF ( (.NOT. ln_itdecay) .AND. ln_incpc) &   
            &   CALL ctl_stop (' if you set ln_incpc to .true. then you need to set ln_itdecay to .true. as well' )

         IF ( (.NOT. ln_incpc) .AND. ln_incpc_only) &   
            &   CALL ctl_stop (' if you set ln_incpc_only to .true. then you need to set ln_incpc to .true. as well' )
         
         WRITE(numout,*) '     time step is    = ',rdt,'you choose to write pcbias at nn_bias_itwrt  = ',nn_bias_itwrt,'and end of iau is rday/rdt=',rday/rdt 
      ENDIF
      IF( .NOT. ln_bias ) RETURN

      IF( .NOT. lalloc ) THEN

         ALLOCATE( tbias(jpi,jpj,jpk)  , &
            &      sbias(jpi,jpj,jpk)  , &
            &      tbias_p(jpi,jpj,jpk), &
            &      sbias_p(jpi,jpj,jpk), &
            &      tbias_i(jpi,jpj,jpk), &
            &      sbias_i(jpi,jpj,jpk), &
            &      rhd_pc(jpi,jpj,jpk) , &
            &      gru_pc(jpi,jpj)     , &
            &      grv_pc(jpi,jpj)       )

         ALLOCATE( fbcoef(jpi,jpj), fbcoef_stscale(jpi,jpj) )
         
         IF( ln_bias_asm ) ALLOCATE(  tbias_asm(jpi,jpj,jpk),     &
            &                         sbias_asm(jpi,jpj,jpk),     &
                                      tbias_asm_out(jpi,jpj,jpk), &  
                                      sbias_asm_out(jpi,jpj,jpk), &  
                                      tbias_p_out(jpi,jpj,jpk),   & 			      
                                      sbias_p_out(jpi,jpj,jpk)    )

         IF( ln_bias_rlx ) ALLOCATE(  tbias_rlx(jpi,jpj,jpk),     &
            &                         sbias_rlx(jpi,jpj,jpk),     &
                                      tbias_rlx_out(jpi,jpj,jpk), &  
                                      sbias_rlx_out(jpi,jpj,jpk)  )

         IF( ln_incpc )    ALLOCATE(  tbias_asm_stscale(jpi,jpj,jpk), &
            &                         sbias_asm_stscale(jpi,jpj,jpk), &
                                      tbias_asm_stscale_out(jpi,jpj,jpk), &
                                      sbias_asm_stscale_out(jpi,jpj,jpk), &
                                      tbias_i_out(jpi,jpj,jpk),   &  
                                      sbias_i_out(jpi,jpj,jpk)   )

         lalloc = .TRUE.

      ENDIF

      IF( ln_bias_ofl ) THEN      ! set sf_tbias_ofl and sf_sbias_ofl strctrs
         !
         ! tbias
         !
         ALLOCATE( sf_tbias_ofl(1), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'bias_init: unable to allocate sf_tbias_ofl structure' )   ;    RETURN
         ENDIF
         ALLOCATE( sf_tbias_ofl(1)%fnow(jpi,jpj,jpk) )
         ALLOCATE( sf_tbias_ofl(1)%fdta(jpi,jpj,jpk,2) )
         
         ! fill structure with values and control print
         CALL fld_fill( sf_tbias_ofl, (/ sn_tbias_ofl /), cn_dir, 'bias_init', 'Offline T bias term ', 'nam_tbias_ofl' )
         !
         ! salinity bias
         ! 
         ALLOCATE( sf_sbias_ofl(1), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'bias_init: unable to allocate sf_sbias_ofl structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_sbias_ofl(1)%fnow(jpi,jpj,jpk) )
         ALLOCATE( sf_sbias_ofl(1)%fdta(jpi,jpj,jpk,2) )
         
         ! fill structure with values and control print
         CALL fld_fill( sf_sbias_ofl, (/ sn_sbias_ofl /), cn_dir, 'bias_init', 'Offline S bias term ', 'nam_sbias_ofl' )
         
      ENDIF

      ! Read total bias
      IF ( ln_bias ) THEN
         tbias(:,:,:)   = 0.0_wp
         sbias(:,:,:)   = 0.0_wp
         tbias_p(:,:,:) = 0.0_wp
         sbias_p(:,:,:) = 0.0_wp
         tbias_i(:,:,:) = 0.0_wp
         sbias_i(:,:,:) = 0.0_wp
         gru_pc(:,:)    = 0.0_wp
         grv_pc(:,:)    = 0.0_wp
         
         IF ( ln_bias_rlx ) THEN
            tbias_rlx(:,:,:) = 0.0_wp
            sbias_rlx(:,:,:) = 0.0_wp
         ENDIF

         IF ( ln_bias_asm ) THEN   !now rlx and asm bias in same file
            tbias_asm(:,:,:) = 0.0_wp
            sbias_asm(:,:,:) = 0.0_wp
            tbias_asm_out(:,:,:) = 0.0_wp
            sbias_asm_out(:,:,:) = 0.0_wp
         ENDIF

         IF ( ln_incpc ) THEN   !incr pressure correction
            tbias_asm_stscale(:,:,:) = 0.0_wp
            sbias_asm_stscale(:,:,:) = 0.0_wp
            tbias_asm_stscale_out(:,:,:) = 0.0_wp
            sbias_asm_stscale_out(:,:,:) = 0.0_wp
         ENDIF


         numbias_tot    = 0
         ! Get bias from file and prevent fail if the file does not exist
         IF(lwp) WRITE(numout,*) 'Opening ',TRIM( cn_bias_tot ) 
         CALL iom_open( cn_bias_tot, numbias_tot, ldstop=.FALSE. )
         
         IF ( numbias_tot > 0 ) THEN      
            ! Could check validity time of bias fields here...
            ! Get the T and S bias data
            IF(lwp) WRITE(numout,*) 'Reading bias fields from tot...'

            !Search for bias from relaxation term if needed. Use same file
            IF ( ln_bias_rlx ) THEN
               IF(lwp) WRITE(numout,*) 'Reading bias fields for bias rlx from file ',cn_bias_tot
               IF( iom_varid( numbias_tot, 'tbias_rlx' ) > 0 ) THEN
                  ! Get the T and S bias data
                  CALL iom_get( numbias_tot, jpdom_autoglo, 'tbias_rlx', tbias_rlx )
                  CALL iom_get( numbias_tot, jpdom_autoglo, 'sbias_rlx', sbias_rlx )
               ELSE
                  CALL ctl_stop( 'Bias relaxation variables not found in ',cn_bias_tot )
               ENDIF
            ENDIF


            !Search for bias from assim term if needed. Use same file
            IF ( ln_bias_asm .and. .not.ln_bsyncro ) THEN
               IF(lwp) WRITE(numout,*) 'Reading a-syncro bias fields for bias asm from file ',cn_bias_tot
               IF( iom_varid( numbias_tot, 'tbias_asm' ) > 0 ) THEN
                  ! Get the T and S bias data
                  CALL iom_get( numbias_tot, jpdom_autoglo, 'tbias_asm', tbias_asm )
                  CALL iom_get( numbias_tot, jpdom_autoglo, 'sbias_asm', sbias_asm )
               ELSE
                  CALL ctl_stop( 'Bias assim variables not found in ',cn_bias_tot )
               ENDIF
            ENDIF

            IF ( ln_incpc .and. .not.ln_bsyncro ) THEN
               IF(lwp) WRITE(numout,*) 'Reading short time scale bias correction fields for bias asm from file ',cn_bias_tot
               IF( iom_varid( numbias_tot, 'tbias_asm_stscale' ) > 0 ) THEN
                  ! Get the T and S bias data
                  CALL iom_get( numbias_tot, jpdom_autoglo, 'tbias_asm_stscale', tbias_asm_stscale )
                  CALL iom_get( numbias_tot, jpdom_autoglo, 'sbias_asm_stscale', sbias_asm_stscale )
               ELSE
                  CALL ctl_stop( 'Short time scale bias assim variables not found in ',cn_bias_tot )
               ENDIF
            ENDIF  



            ! Close the file
            CALL iom_close(numbias_tot)

         ELSE 
            IF(lwp) WRITE(numout,*) 'No bias file found so T and S bias fields are set to zero'
         ENDIF

      ENDIF

     ! for the time being, the bias_asm is read in the same file as
     ! bias_rlx
     ! Implications: bias_asm is estimated/evolved in time in the second outer 
     !               loop only, when the assimilation increments are ready.
     !               bias_asm is kept cte during the first outer loop.
     !              => Assyncronous bias correction.
     ! Alternative: Syncronous bias correction:
     !              bias_asm estimated/evolved in the first outer loop
     !              with the asm incrments of the previous cycle.
     !              bias_asm kept cte during the second outer loop.
     !              Implication: bias_asm should be estimated really in the
     !              inner loop.
      IF ( ln_bsyncro ) THEN
      ! Read bias from assimilation from a separate file
      IF ( ln_bias_asm ) THEN
         tbias_asm(:,:,:) = 0.0_wp
         sbias_asm(:,:,:) = 0.0_wp
         numbias_asm      = 0
         ! Get bias from file and prevent fail if the file does not exist
         IF(lwp) WRITE(numout,*) 'Opening file for syncro assim bias ',TRIM( cn_bias_asm ) 
         CALL iom_open( cn_bias_asm, numbias_asm, ldstop=.FALSE. )
         
         IF ( numbias_asm > 0 ) THEN      
            ! Could check validity time of bias fields here...
            
            ! Get the T and S bias data
            IF(lwp) WRITE(numout,*) 'Reading syncro bias fields from asm from file ',cn_bias_asm
            CALL iom_get( numbias_asm, jpdom_autoglo, 'tbias_asm', tbias_asm )
            CALL iom_get( numbias_asm, jpdom_autoglo, 'sbias_asm', sbias_asm )
           
!  this is only applicable if tbias_asm were to be calculated in the inner loop
            tbias_asm(:,:,:) = tbias_asm(:,:,:) * rdt / bias_time_unit_asm
            sbias_asm(:,:,:) = sbias_asm(:,:,:) * rdt / bias_time_unit_asm
            
            ! Close the file
            CALL iom_close(numbias_asm)
            
         ELSE 
            IF(lwp) WRITE(numout,*) 'No bias file found from asm so T and S bias fields are set to zero'
         ENDIF
         
      ENDIF

      ENDIF

      !latitudinal dependence of partition coeficients. Adhoc
      IF ( nn_lat_ramp == 1 ) THEN
         ! Use the inertial ramp.
         lenscl_bias = ( rn_maxlat_bias - rn_minlat_bias )*2._wp
         WHERE ( abs( gphit(:,:) ) <= rn_minlat_bias )         
            fbcoef(:,:) = 0._wp          
         ELSEWHERE ( abs( gphit(:,:) ) >= rn_maxlat_bias )         
            fbcoef(:,:) = 1._wp                    
         ELSEWHERE        
            fbcoef(:,:) = 1._wp - exp( -( abs( gphit(:,:) ) - rn_minlat_bias ) &
                           * ( abs( gphit(:,:) ) - rn_minlat_bias ) / lenscl_bias )                         
         ENDWHERE 
      ELSEIF ( nn_lat_ramp == 2 ) THEN   
         ! Use a linear ramp consist with the geostrophic velocity balance ramp in NEMOVAR
     
         WHERE ( abs( gphit(:,:) ) <= rn_minlat_bias )
            fbcoef(:,:) = 0._wp
         ELSEWHERE ( abs( gphit(:,:) ) >= rn_maxlat_bias ) 
            fbcoef(:,:) = 1._wp
         ELSEWHERE
            fbcoef(:,:) = 1._wp - ((rn_maxlat_bias - abs( gphit(:,:)))/(rn_maxlat_bias - rn_minlat_bias))
         ENDWHERE
      ELSE
         fbcoef(:,:) = 0.0_wp
         fctamp      = 0.0_wp
         fbcoef_stscale(:,:) = 0.0_wp
      ENDIF


      IF ( ln_incpc) THEN
         minlat_bias = 3.0_wp
         maxlat_bias = 8.0_wp  
         WHERE ( abs( gphit(:,:) ) <= minlat_bias )
            fbcoef_stscale(:,:)=0._wp
         ELSEWHERE ( abs( gphit(:,:) ) >= maxlat_bias ) 
            fbcoef_stscale(:,:)=1._wp
         ELSEWHERE
            fbcoef_stscale(:,:)=1._wp - ((maxlat_bias - abs( gphit(:,:)))/(maxlat_bias-minlat_bias))
         ENDWHERE   
      ENDIF 


   END SUBROUTINE bias_init

   SUBROUTINE tra_bias ( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tra_bias  ***
      !!
      !! ** Purpose : Update bias field every time step
      !!
      !! ** Method  : add contributions to bias from 3 terms
      !!
      !! ** Action  : Bias from assimilation (read in bias_init) 
      !!              Bias from relaxation term is estimated according to
      !!              the prescribed time evolution of the bias
      !!              Bias from ofl term is read from external file
      !!              The difference contributions are added and the partition
      !!              into direct bias in T/S and pressure perfomed.
      !!
      !! History :  09-03  (M. Balmaseda) 
      !!----------------------------------------------------------------------
      !! called every timestep after dta_sst if ln_bias true.

      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(IN) ::   kt             ! ocean time-step index
      !! * Local variables
      INTEGER             ::   ji,jj,jk, it   ! local loop index
      REAL(wp)            ::   tsclf          ! time saling factor 
      REAL(wp)            ::   fb_t_asm_max, fb_t_rlx_max, fb_t_ofl_max
      REAL(wp)            ::   ztfrac, ztsday
      REAL(wp)            ::   zfrac, zfrac1 ! temporal weights for inst pcbias (names could be changed)
      REAL(wp)            ::   zdecay        ! used in inst pcorr
      REAL(wp), DIMENSION(jpi,jpj) :: zcof1, zcof2

      IF ( .NOT. ln_bias ) RETURN
      fb_t_rlx_max   = MIN(fb_t_rlx*fctamp,1.0_wp)
      fb_t_asm_max   = MIN(fb_t_asm*fctamp,1.0_wp)
      fb_t_ofl_max   = MIN(fb_t_ofl*fctamp,1.0_wp)

      tbias(:,:,:)   = 0.0_wp
      sbias(:,:,:)   = 0.0_wp
      tbias_p(:,:,:) = 0.0_wp
      sbias_p(:,:,:) = 0.0_wp
      tbias_i(:,:,:) = 0.0_wp
      sbias_i(:,:,:) = 0.0_wp

      IF ( ln_bias_asm ) THEN
      
         tsclf = 1
         IF ( .NOT.ln_bsyncro ) tsclf = rdt / bias_time_unit_asm 
         zcof1(:,:) = tsclf * ( ( 1.0_wp - fbcoef(:,:) ) * fb_t_asm + &
            &                              fbcoef(:,:)   * fb_t_asm_max )
         zcof2(:,:) = ( 1.0_wp - fbcoef(:,:) ) * fb_p_asm

         IF ( ln_itdecay ) THEN   !decay the pressure correction at each time step
	 
	    ztsday  = rday / REAL(rdt,wp)

            zdecay = (1-ztscale)**(1/REAL(ztsday,wp)) ! used in ipc
            zfrac1 = zdecay**REAL(kt,wp) ! used in ipc
            IF ( zfrac1 <= 0.0 ) zfrac1 = 0.0_wp

            IF( ln_asmiau .and. ln_trainc ) THEN  ! nudge in increments and decay historical contribution
               
               IF ( kt <= nitiaufin ) THEN  ! During IAU calculate the fraction of increments to apply at each time step

                  ztfrac = REAL(kt,wp) / REAL(nitiaufin,wp)  ! nudging factor during the IAU
	       
                  IF (lwp) THEN
                     WRITE(numout,*) 'tra_bias : bias weights'
                     WRITE(numout,*) '~~~~~~~~~~~~'
                     WRITE(numout,* ) "proportion of  increment applied in pcbias ",ztfrac
                     WRITE(numout,* ) "proportion of  historical pcbias applied ",t_asm_mem**(REAL(kt,wp)/ztsday)
                  ENDIF

                  DO jk = 1, jpkm1  
                     tbias(:,:,jk) = tbias(:,:,jk) +                            &
                     &                ( t_asm_mem**(REAL(kt,wp)/ztsday) * tbias_asm(:,:,jk)  +                    &
                     &                ztfrac * t_asm_upd * t_bkginc(:,:,jk) * tmask(:,:,jk) ) * zcof1(:,:)
                     sbias(:,:,jk) = sbias(:,:,jk) +                            &
                     &               ( t_asm_mem**(REAL(kt,wp)/ztsday) * sbias_asm(:,:,jk)  +                     &
                     &               ztfrac * t_asm_upd * s_bkginc(:,:,jk) * tmask(:,:,jk) ) * zcof1(:,:)

                     tbias_p(:,:,jk) = tbias_p(:,:,jk) +                        &
                     &               ( t_asm_mem**(REAL(kt,wp)/ztsday) * tbias_asm(:,:,jk)  +                     &
                     &               ztfrac * t_asm_upd * t_bkginc(:,:,jk) * tmask(:,:,jk) ) * zcof2(:,:)
                     sbias_p(:,:,jk) = sbias_p(:,:,jk) +                        &  
                     &               ( t_asm_mem**(REAL(kt,wp)/ztsday) * sbias_asm(:,:,jk)  +                     &
                     &               ztfrac * t_asm_upd * s_bkginc(:,:,jk) * tmask(:,:,jk) ) * zcof2(:,:)
                  ENDDO

                  IF (ln_incpc) THEN

                     IF (lwp) THEN
                       WRITE(numout,*) 'tra_bias : bias weights'
                       WRITE(numout,*) '~~~~~~~~~~~~'
                       WRITE(numout,* ) "IPC - proportion of  increment applied in pcbias ",ztfrac
                       WRITE(numout,* ) "IPC - proportion of  historical pcbias applied ",zfrac1
                     ENDIF

                     DO jk = 1, jpkm1

                        tbias_i(:,:,jk) =  ( t_bkginc(:,:,jk) * tmask(:,:,jk)* zwgt * ztfrac * (1.0 - fbcoef_stscale(:,:)) )         &
                        &                + ( tbias_asm_stscale(:,:,jk) * zfrac1 * (1.0 - fbcoef_stscale(:,:)) )
                        sbias_i(:,:,jk) =  ( s_bkginc(:,:,jk) * tmask(:,:,jk)* zwgt * ztfrac * (1.0 - fbcoef_stscale(:,:)) )         &
                        &                + ( sbias_asm_stscale(:,:,jk) * zfrac1 * (1.0 - fbcoef_stscale(:,:)) )

                     ENDDO

                  ENDIF

                  IF ( .not.ln_bsyncro ) THEN 

                     IF ( kt == nn_bias_itwrt ) THEN

                        DO jk = 1, jpkm1
                           tbias_asm_out(:,:,jk) =  t_asm_mem**(REAL(kt,wp)/ztsday) * tbias_asm(:,:,jk)  +       &
                           &                     ztfrac * t_asm_upd * t_bkginc(:,:,jk) * tmask(:,:,jk)
                           sbias_asm_out(:,:,jk) =  t_asm_mem**(REAL(kt,wp)/ztsday) * sbias_asm(:,:,jk)  +       &
                           &                     ztfrac * t_asm_upd * s_bkginc(:,:,jk) * tmask(:,:,jk)
                         END DO

                        IF ( ln_incpc) THEN
                           DO jk = 1, jpkm1
                              tbias_asm_stscale_out(:,:,jk) = ( t_bkginc(:,:,jk) * tmask(:,:,jk) *  zwgt * ztfrac ) + ( tbias_asm_stscale(:,:,jk) * zfrac1 )
                              sbias_asm_stscale_out(:,:,jk) = ( s_bkginc(:,:,jk) * tmask(:,:,jk) *  zwgt * ztfrac ) + ( sbias_asm_stscale(:,:,jk) * zfrac1 )
                           ENDDO
                        ENDIF

                     ENDIF

                  ENDIF

                  ! update the historical component with the increments at the end of the IAU
                  IF ( kt == nitiaufin ) THEN
                     DO jk = 1, jpkm1
                        tbias_asm(:,:,jk) =  t_asm_mem**(REAL(kt,wp)/ztsday) * tbias_asm(:,:,jk)  +       &
                        &                     ztfrac * t_asm_upd * t_bkginc(:,:,jk) * tmask(:,:,jk)
                        sbias_asm(:,:,jk) =  t_asm_mem**(REAL(kt,wp)/ztsday) * sbias_asm(:,:,jk)  +       &
                        &                     ztfrac * t_asm_upd * s_bkginc(:,:,jk) * tmask(:,:,jk)
                     END DO

                     IF (ln_incpc) THEN
                        DO jk = 1, jpkm1
                           tbias_asm_stscale(:,:,jk) = ( t_bkginc(:,:,jk) * tmask(:,:,jk) * zwgt * ztfrac ) + ( tbias_asm_stscale(:,:,jk) * zfrac1 )
                           sbias_asm_stscale(:,:,jk) = ( s_bkginc(:,:,jk) * tmask(:,:,jk) * zwgt * ztfrac ) + ( sbias_asm_stscale(:,:,jk) * zfrac1 )
                        ENDDO
                     ENDIF
     
                  ENDIF
               
               ELSE ! decay pressure correction from combined historical component and increments after IAU

                  ztfrac=t_asm_mem**(REAL(kt-nitiaufin,wp)/REAL(nitiaufin,wp))  ! decay from end of IAU
                  
                  DO jk = 1, jpkm1
                     tbias(:,:,jk) = tbias(:,:,jk) +                            &
                     &                ( ztfrac * tbias_asm(:,:,jk) ) * zcof1(:,:)
                     sbias(:,:,jk) = sbias(:,:,jk) +                            &
                     &               (  ztfrac * sbias_asm(:,:,jk) ) * zcof1(:,:)
                     tbias_p(:,:,jk) = tbias_p(:,:,jk) +                        &
                     &               (  ztfrac * tbias_asm(:,:,jk) ) * zcof2(:,:)
                     sbias_p(:,:,jk) = sbias_p(:,:,jk) +                        &  
                     &               ( ztfrac * sbias_asm(:,:,jk) ) * zcof2(:,:)
                  ENDDO

                 IF (ln_incpc) THEN

                   zfrac  = zdecay**REAL(kt-nitiaufin,wp) 
                   IF ( zfrac <= 0.0 ) zfrac = 0.0_wp
   
                   DO jk = 1, jpkm1
                      tbias_i(:,:,jk) =  tbias_asm_stscale(:,:,jk) * zfrac * (1.0 - fbcoef_stscale(:,:))
                      sbias_i(:,:,jk) =  sbias_asm_stscale(:,:,jk) * zfrac * (1.0 - fbcoef_stscale(:,:))
                   ENDDO

                   IF (lwp) THEN
                      WRITE(numout,*) 'tra_bias : bias weights'
                      WRITE(numout,*) '~~~~~~~~~~~~'
                      WRITE(numout,* ) "IPC - proportion of increments + historical pcbias applied ",zfrac
                   ENDIF

                 ENDIF

                 IF (lwp) THEN
                    WRITE(numout,*) 'tra_bias : bias weights'
                    WRITE(numout,*) '~~~~~~~~~~~~'
                    WRITE(numout,* ) "proportion of increments + historical pcbias applied ",ztfrac
                 ENDIF

                 IF ( ln_trainc .and. .not.ln_bsyncro ) THEN 
                    IF ( kt == nn_bias_itwrt ) THEN
                       DO jk = 1, jpkm1
                          tbias_asm_out(:,:,jk) =  ztfrac * tbias_asm(:,:,jk) 
                          sbias_asm_out(:,:,jk) =  ztfrac * sbias_asm(:,:,jk) 
                       END DO

                       IF ( ln_incpc ) THEN
                          IF (lwp) WRITE(numout,*) 'after end of IAU - IPC - saving s/tbias_asm_stscale'
                             DO jk = 1, jpkm1
                                tbias_asm_stscale_out(:,:,jk) =  tbias_asm_stscale(:,:,jk) * zfrac
                                sbias_asm_stscale_out(:,:,jk) =  sbias_asm_stscale(:,:,jk) * zfrac
                             ENDDO
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF

            ELSE ! no assimilation increments, simply decay pressure correction (e.g for forecasts)
                 ! this occurs at obsoper step as well ( ln_asmiau is false ) 

               DO jk = 1, jpkm1
                  tbias(:,:,jk) = tbias(:,:,jk) +                                                         &
                  &               ( t_asm_mem**(REAL(kt,wp)/ztsday) * tbias_asm(:,:,jk) ) * zcof1(:,:)
                  sbias(:,:,jk) = sbias(:,:,jk) +                                                         &
                  &               ( t_asm_mem**(REAL(kt,wp)/ztsday) * sbias_asm(:,:,jk) ) * zcof1(:,:)
                  tbias_p(:,:,jk) = tbias_p(:,:,jk) +                                                     &
                  &               ( t_asm_mem**(REAL(kt,wp)/ztsday) * tbias_asm(:,:,jk) ) * zcof2(:,:)
                  sbias_p(:,:,jk) = sbias_p(:,:,jk) +                                                     &  
                  &               ( t_asm_mem**(REAL(kt,wp)/ztsday) * sbias_asm(:,:,jk) ) * zcof2(:,:)
               ENDDO

               IF (lwp) THEN
                  WRITE(numout,*) 'tra_bias : bias weights'
                  WRITE(numout,*) '~~~~~~~~~~~~'
                  WRITE(numout,* ) "proportion of historical pcbias applied ",t_asm_mem**(REAL(kt,wp)/ztsday)
               ENDIF

               IF (ln_incpc) THEN
                   IF (lwp) WRITE(numout,*) 'obsoper or forecast mode - IPC - computing tbias_i and sbias_i'
                   DO jk = 1, jpkm1
                      tbias_i(:,:,jk) = tbias_i(:,:,jk) + ( tbias_asm_stscale(:,:,jk) * zfrac1 * (1.0 - fbcoef_stscale(:,:)) )
                      sbias_i(:,:,jk) = sbias_i(:,:,jk) + ( sbias_asm_stscale(:,:,jk) * zfrac1 * (1.0 - fbcoef_stscale(:,:)) )
                   ENDDO
                   tbias_i(:,:,:) =tbias_i(:,:,:)*tmask(:,:,:)
                   sbias_i(:,:,:) =sbias_i(:,:,:)*tmask(:,:,:)
               ENDIF
            ENDIF
 
         ELSE ! maintain a constant pressure correction  

            DO jk = 1, jpkm1
               tbias(:,:,jk) = tbias(:,:,jk) + tbias_asm(:,:,jk) * zcof1(:,:)
               sbias(:,:,jk) = sbias(:,:,jk) + sbias_asm(:,:,jk) * zcof1(:,:)
               tbias_p(:,:,jk) = tbias_p(:,:,jk) + tbias_asm(:,:,jk) * zcof2(:,:)
               sbias_p(:,:,jk) = sbias_p(:,:,jk) + sbias_asm(:,:,jk) * zcof2(:,:)
            END DO     

            IF( ln_asmiau .and. ln_trainc .and. .not.ln_bsyncro ) THEN   
            ! if last outer loop (ln_asmiau=true and ln_trainc=true). t/sbias_asm
            ! is updated, only once (end of run) taking into account units.
               IF ( kt == nn_bias_itwrt ) THEN
                 IF(lwp) WRITE(numout,*)' estimating asm bias at timestep: ',kt
                 DO jk = 1, jpkm1
                   tbias_asm_out(:,:,jk) = t_asm_mem * tbias_asm(:,:,jk)  +             &
                   &                      t_asm_upd * t_bkginc(:,:,jk) * tmask(:,:,jk)
                   sbias_asm_out(:,:,jk) = t_asm_mem * sbias_asm(:,:,jk) +             &
                   &                      t_asm_upd * s_bkginc(:,:,jk) * tmask(:,:,jk)               
                 END DO
               ENDIF
             ENDIF
  
         ENDIF

      ENDIF


#if   defined key_tradmp 
      ! Time evolution of bias from relaxation
      IF ( ln_bias_rlx ) THEN
         tbias_rlx(:,:,:) = t_rlx_mem * tbias_rlx(:,:,:) + &
            &               t_rlx_upd * ttrdmp(:,:,:) * rdt
         sbias_rlx(:,:,:) = t_rlx_mem * sbias_rlx(:,:,:) + &
            &               t_rlx_upd * strdmp(:,:,:) * rdt
         zcof1(:,:) = ( 1.0_wp - fbcoef(:,:) ) * fb_t_rlx +      &
            &                    fbcoef(:,:)   * fb_t_rlx_max
         zcof2(:,:) = ( 1.0_wp - fbcoef(:,:) ) * fb_p_rlx
         DO jk = 1, jpkm1
            tbias(:,:,jk)   = tbias(:,:,jk) + tbias_rlx(:,:,jk) * zcof1(:,:)
            sbias(:,:,jk)   = sbias(:,:,jk) + sbias_rlx(:,:,jk) * zcof1(:,:)
            tbias_p(:,:,jk) = tbias_p(:,:,jk) + tbias_rlx(:,:,jk) * zcof2(:,:)
            sbias_p(:,:,jk) = sbias_p(:,:,jk) + sbias_rlx(:,:,jk) * zcof2(:,:)
         ENDDO
         IF ( kt == nn_bias_itwrt ) THEN
            tbias_rlx_out(:,:,:) = tbias_rlx(:,:,:)
            sbias_rlx_out(:,:,:) = sbias_rlx(:,:,:)
         ENDIF
      ENDIF
#endif
      ! ofline bias
      IF ( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*) ' tra_bias: ln_bias_ofl = ',ln_bias_ofl
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~'
      ENDIF
      IF ( ln_bias_ofl ) THEN
         IF(lwp) WRITE(numout,*) 'reading offline bias'
         CALL fld_read( kt, 1, sf_tbias_ofl ) 
         CALL fld_read( kt, 1, sf_sbias_ofl ) 

         zcof1(:,:) = ( 1.0_wp - fbcoef(:,:) ) * fb_t_ofl +           &
            &                    fbcoef(:,:)   * fb_t_ofl_max
         zcof2(:,:) = ( 1.0_wp - fbcoef(:,:) ) * fb_p_ofl
         DO jk = 1, jpkm1
            tbias(:,:,jk)   = tbias(:,:,jk) + sf_tbias_ofl(1)%fnow(:,:,jk) * zcof1(:,:)
            sbias(:,:,jk)   = sbias(:,:,jk) + sf_sbias_ofl(1)%fnow(:,:,jk) * zcof1(:,:)
            tbias_p(:,:,jk) = tbias_p(:,:,jk) + sf_tbias_ofl(1)%fnow(:,:,jk) * zcof2(:,:)
            sbias_p(:,:,jk) = sbias_p(:,:,jk) + sf_sbias_ofl(1)%fnow(:,:,jk) * zcof2(:,:)
         ENDDO
      ENDIF


      !apply bias on tracers if needed      
      IF ( ln_bias_ts_app ) THEN
         
         ! Add the bias directely to T/s
         tsa(:,:,:,jp_tem) = tsa(:,:,:,jp_tem) + tmask(:,:,:) * tbias(:,:,:) / rdt
         tsa(:,:,:,jp_sal) = tsa(:,:,:,jp_sal) + tmask(:,:,:) * sbias(:,:,:) / rdt

         ! lateral boundary conditions (is this needed?)
         CALL lbc_lnk( tsa(:,:,:,jp_tem), 'T', 1.0_wp )
         CALL lbc_lnk( tsa(:,:,:,jp_sal), 'T', 1.0_wp )   

      ENDIF

      IF ( kt == nn_bias_itwrt ) THEN
            tbias_p_out(:,:,:)   = tbias_p(:,:,:)
            sbias_p_out(:,:,:)   = sbias_p(:,:,:)
            IF (ln_incpc) THEN
               tbias_i_out(:,:,:)   = tbias_i(:,:,:)
               sbias_i_out(:,:,:)   = sbias_i(:,:,:)            
            ENDIF
      ENDIF

   END SUBROUTINE tra_bias

   SUBROUTINE dyn_bias( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_bias  ***
      !!
      !! ** Purpose :   Computes rhd_pc, gru/v_pc bias corrected
      !!                for hydrostatic pressure gradient
      !!                depending on time step (semi-implicit/centered)
      !!                If partial steps computes bottom pressure gradient.
      !!                These correction terms will affect only the dynamical
      !!                component (pressure gradient calculation), but not
      !!                the isopycnal calculation for the lateral mixing.
      !! 
      !! ** Method  :   At this stage of the computation, ta and sa are the 
      !!                after temperature and salinity. If semi-implicit, these
      !!                are used to compute rho and bottom pressure gradient.
      !!                If centered, tb,sb are used instead.
      !!                If bias key is activated, the temperature,salinity are
      !!                bias corrected in the calculation of the density fields
      !!                used in the pressure gradient calculation.
      !!
      !!
      !! ** Action  : - rhd_pc ready. rhop will be overwriten later
      !!              - if ln_zps, bottom density gradients gru/v_pc ready.
      !!----------------------------------------------------------------------
      !!
      !! * Arguments
      INTEGER, INTENT(IN) ::   kt    ! ocean time-step index
      !! * Local variables
      REAL(wp) :: tsw(jpi,jpj,jpk,jpts)
      !!
      !!----------------------------------------------------------------------
      !
      ! gtu,gsu,gtv,gsv rhop will be overwritten later in step.
      ! 
      IF( ln_dynhpg_imp  ) THEN                             ! semi-implicit hpg
         tsw(:,:,:,jp_tem) = tsa(:,:,:,jp_tem) - tbias_p(:,:,:)
         tsw(:,:,:,jp_sal) = tsa(:,:,:,jp_sal) - sbias_p(:,:,:)
         IF ( ln_incpc ) THEN
            tsw(:,:,:,jp_tem) = tsa(:,:,:,jp_tem) - tbias_p(:,:,:) - tbias_i(:,:,:)
            tsw(:,:,:,jp_sal) = tsa(:,:,:,jp_sal) - sbias_p(:,:,:) - sbias_i(:,:,:)
            IF ( ln_incpc_only ) THEN
               IF(lwp) WRITE(numout,*) 'ln_incpc_only =', ln_incpc_only, 'tsw updated with IPC only and ln_dynhpg_imp = ',ln_dynhpg_imp
               tsw(:,:,:,jp_tem) = tsa(:,:,:,jp_tem) - tbias_i(:,:,:)
               tsw(:,:,:,jp_sal) = tsa(:,:,:,jp_sal) - sbias_i(:,:,:)
            ENDIF
         ENDIF
      ELSE
         tsw(:,:,:,jp_tem) = tsb(:,:,:,jp_tem) - tbias_p(:,:,:)
         tsw(:,:,:,jp_sal) = tsb(:,:,:,jp_sal) - sbias_p(:,:,:)
         IF ( ln_incpc ) THEN
            tsw(:,:,:,jp_tem) = tsb(:,:,:,jp_tem) - tbias_p(:,:,:) - tbias_i(:,:,:)
            tsw(:,:,:,jp_sal) = tsb(:,:,:,jp_sal) - sbias_p(:,:,:) - sbias_i(:,:,:)
            IF ( ln_incpc_only ) THEN
               IF(lwp) WRITE(numout,*) 'ln_incpc_only =', ln_incpc_only, 'tsw updated with IPC only and ln_dynhpg_imp = ',ln_dynhpg_imp
               tsw(:,:,:,jp_tem) = tsb(:,:,:,jp_tem) - tbias_i(:,:,:)
               tsw(:,:,:,jp_sal) = tsb(:,:,:,jp_sal) - sbias_i(:,:,:)
            ENDIF
         ENDIF
      ENDIF

      IF(lwp) WRITE(numout,*) 'dyn_bias( kt ) calculating rhd_pc, kt =', kt
      CALL eos( tsw, rhd_pc, rhop, fsdept_n(:,:,:) )
      
      CALL lbc_lnk( rhd_pc, 'T', 1.0_wp )


      ! Partial steps: now horizontal gradient of t,s,rd 
      ! at the bottom ocean level

      IF( ln_zps    )  THEN
         CALL zps_hde( kt, jpts, tsw, gtsu, gtsv,  &
            &          rhd_pc, gru_pc , grv_pc  )
      ENDIF

   END SUBROUTINE dyn_bias
   
   SUBROUTINE bias_opn( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE bias_opn  ***
      !!                     
      !! ** Purpose :  open bias restart file following the same logic as the
      !!               standard restarts.
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(IN) ::   kt     ! ocean time-step
      !! * Local variables
      CHARACTER(LEN=20)   ::   clbkt    ! ocean time-step deine as a character
      CHARACTER(LEN=50)   ::   clbias_tot   ! total bias restart file name
      !!----------------------------------------------------------------------
      !
      IF( lrst_oce .AND. .NOT.lrst_bias ) THEN       ! create bias file
         IF( nitend > 999999999 ) THEN   ;   WRITE(clbkt, *       ) nitend
         ELSE                            ;   WRITE(clbkt, '(i8.8)') nitend
         ENDIF
         clbias_tot = TRIM(cexper)//"_"//TRIM(ADJUSTL(clbkt))//"_"//TRIM(cn_bias_tot)
         IF(lwp) THEN
            WRITE(numout,*)
            SELECT CASE ( jprstlib )
            CASE ( jprstdimg )   ;   WRITE(numout,*) '             open tot bias restart binary file: '//clbias_tot
            CASE DEFAULT         ;   WRITE(numout,*) '             open tot bias restart NetCDF file: '//clbias_tot
            END SELECT
         ENDIF
         CALL iom_open( clbias_tot, numbias_tot , ldwrt = .TRUE., kiolib = jprstlib )
         lrst_bias=.TRUE.

      ENDIF
      !
   END SUBROUTINE bias_opn

   SUBROUTINE bias_wrt( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE bias_wrt  ***
      !!                     
      !! ** Purpose :   Write bias restart fields in the format corresponding to jprstlib
      !!
      !! ** Method  :   Write in numbias_tot when kt == nitend in output
      !!                file, save fields which are necessary for restart.
      !!
      !! Changed the timestep for writing to nitend rather than nitrst as this causes a
      !! problem if the nstock list is used to determine the restart output times and
      !! means that the bias is not output at all. M. Martin. 08/2011.
      !! Need to check with Magdalena that this is ok for ECMWF.
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(IN) ::   kt   ! ocean time-step
      !!----------------------------------------------------------------------
      !                                                                     ! the begining of the run [s]

      IF ( ln_bias_rlx ) THEN
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'tbias_rlx' , tbias_rlx_out )   
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'sbias_rlx' , sbias_rlx_out )   
      ENDIF
      
      IF ( ln_bias_asm ) THEN
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'tbias_asm' , tbias_asm_out )   
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'sbias_asm' , sbias_asm_out )   
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'tbias_p'   , tbias_p_out )
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'sbias_p'   , sbias_p_out )         
      ENDIF

      IF ( ln_incpc ) THEN
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'tbias_asm_stscale' , tbias_asm_stscale_out )   
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'sbias_asm_stscale' , sbias_asm_stscale_out ) 
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'tbias_i'   , tbias_i_out )
         CALL iom_rstput( nn_bias_itwrt, nn_bias_itwrt, numbias_tot, 'sbias_i'   , sbias_i_out )   
      ENDIF
      
      IF( kt == nitend ) THEN
         CALL iom_close( numbias_tot )     ! close the restart file (only at last time step)
         lrst_bias = .FALSE.
      ENDIF
      !
   END SUBROUTINE bias_wrt

END MODULE bias
