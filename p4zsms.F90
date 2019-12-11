MODULE p4zsms
   !!======================================================================
   !!                         ***  MODULE p4zsms  ***
   !! TOP :   PISCES Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4zsms         :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE trcdta
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zbio          !  Biological model
   USE p4zche          !  Chemical model
   USE p4zlys          !  Calcite saturation
   USE p4zflx          !  Gas exchange
   USE p4zsbc          !  External source of nutrients
   USE p4zsed          !  Sedimentation
   USE p4zint          !  time interpolation
   USE p4zrem          !  remineralisation
   USE iom             !  I/O manager
   USE trd_oce         !  Ocean trends variables
   USE trdtrc          !  TOP trends variables
   USE sedmodel        !  Sediment model
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sms_init   ! called in p4zsms.F90
   PUBLIC   p4z_sms        ! called in p4zsms.F90

   REAL(wp) :: alkbudget, no3budget, silbudget, ferbudget, po4budget
   REAL(wp) :: xfact1, xfact2, xfact3
   INTEGER ::  numco2, numnut, numnit  !: logical unit for co2 budget

   !!* Array used to indicate negative tracer values
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     !: ???


   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sms  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of PISCES bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER ::   ji, jj, jk, jnt, jn, jl
      REAL(wp) ::  ztra
#if defined key_kriest
      REAL(wp) ::  zcoef1, zcoef2
#endif
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sms')
      !
      IF( kt == nittrc000 ) THEN
        !
        ALLOCATE( xnegtr(jpi,jpj,jpk) )
        !
        CALL p4z_che                              ! initialize the chemical constants
        !
        IF( .NOT. ln_rsttr ) THEN  ;   CALL p4z_ph_ini   !  set PH at kt=nit000 
        ELSE                       ;   CALL p4z_rst( nittrc000, 'READ' )  !* read or initialize all required fields 
        ENDIF
        !
      ENDIF
      !
      IF( ln_pisdmp .AND. MOD( kt - nn_dttrc, nn_pisdmp ) == 0 )   CALL p4z_dmp( kt )      ! Relaxation of some tracers
      !
      !                                                                    !   set time step size (Euler/Leapfrog)
      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN   ;    rfact = rdttrc(1)     !  at nittrc000
      ELSEIF( kt <= nittrc000 + nn_dttrc )                          THEN   ;    rfact = 2. * rdttrc(1)   ! at nittrc000 or nittrc000+nn_dttrc (Leapfrog)
      ENDIF
      !
      IF( ( ln_top_euler .AND. kt == nittrc000 )  .OR. ( .NOT.ln_top_euler .AND. kt <= nittrc000 + nn_dttrc ) ) THEN
         rfactr  = 1. / rfact
         rfact2  = rfact / FLOAT( nrdttrc )
         rfact2r = 1. / rfact2
         xstep = rfact2 / rday         ! Time step duration for biology
         IF(lwp) WRITE(numout,*) 
         IF(lwp) WRITE(numout,*) '    Passive Tracer  time step    rfact  = ', rfact, ' rdt = ', rdttra(1)
         IF(lwp) write(numout,*) '    PISCES  Biology time step    rfact2 = ', rfact2
         IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN
         DO jn = jp_pcs0, jp_pcs1              !   SMS on tracer without Asselin time-filter
            trb(:,:,:,jn) = trn(:,:,:,jn)
         END DO
      ENDIF
      !
      IF( ndayflxtr /= nday_year ) THEN      ! New days
         !
         ndayflxtr = nday_year

         IF(lwp) write(numout,*)
         IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
         IF(lwp) write(numout,*) '~~~~~~'

         CALL p4z_che              ! computation of chemical constants
         CALL p4z_int( kt )        ! computation of various rates for biogeochemistry
         !
      ENDIF

      IF( ll_sbc ) CALL p4z_sbc( kt )   ! external sources of nutrients 

      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL p4z_bio( kt, jnt )   ! Biology
         CALL p4z_lys( kt, jnt )   ! Compute CaCO3 saturation
         CALL p4z_sed( kt, jnt )   ! Surface and Bottom boundary conditions
         CALL p4z_flx( kt, jnt )   ! Compute surface fluxes
         !
         xnegtr(:,:,:) = 1.e0
         DO jn = jp_pcs0, jp_pcs1
            DO jk = 1, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF( ( trb(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN
                        ztra             = ABS( trb(ji,jj,jk,jn) ) / ( ABS( tra(ji,jj,jk,jn) ) + rtrn )
                        xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
                     ENDIF
                 END DO
               END DO
            END DO
         END DO
         !                                ! where at least 1 tracer concentration becomes negative
         !                                ! 
         DO jn = jp_pcs0, jp_pcs1
           trb(:,:,:,jn) = trb(:,:,:,jn) + xnegtr(:,:,:) * tra(:,:,:,jn)
         END DO
        !
         DO jn = jp_pcs0, jp_pcs1
            tra(:,:,:,jn) = 0._wp
         END DO
         !
         IF( ln_top_euler ) THEN
            DO jn = jp_pcs0, jp_pcs1
               trn(:,:,:,jn) = trb(:,:,:,jn)
            END DO
         ENDIF
      END DO

#if defined key_kriest
      ! 
      zcoef1 = 1.e0 / xkr_massp 
      zcoef2 = 1.e0 / xkr_massp / 1.1
      DO jk = 1,jpkm1
         trb(:,:,jk,jpnum) = MAX(  trb(:,:,jk,jpnum), trb(:,:,jk,jppoc) * zcoef1 / xnumm(jk)  )
         trb(:,:,jk,jpnum) = MIN(  trb(:,:,jk,jpnum), trb(:,:,jk,jppoc) * zcoef2              )
      END DO
      !
#endif
      !
      !
      IF( l_trdtrc ) THEN
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn), jn, jptra_sms, kt )   ! save trends
         END DO
      END IF
      !
      IF( lk_sed ) THEN 
         !
         CALL sed_model( kt )     !  Main program of Sediment model
         !
         DO jn = jp_pcs0, jp_pcs1
           CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
         END DO
         !
      ENDIF
      !
      IF( lrst_trc )  CALL p4z_rst( kt, 'WRITE' )  !* Write PISCES informations in restart file 
      !

      IF( lk_iomput .OR. ln_check_mass )  CALL p4z_chk_mass( kt ) ! Mass conservation checking

      IF ( lwm .AND. kt == nittrc000 ) CALL FLUSH    ( numonp )     ! flush output namelist PISCES
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sms')
      !
      !
   END SUBROUTINE p4z_sms

   SUBROUTINE p4z_sms_init
      !!----------------------------------------------------------------------
      !!                     ***  p4z_sms_init  ***  
      !!
      !! ** Purpose :   read PISCES namelist
      !!
      !! ** input   :   file 'namelist.trc.s' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      NAMELIST/nampisbio/ nrdttrc, wsbio, xkmort, ferat3, wsbio2, niter1max, niter2max
#if defined key_kriest
      NAMELIST/nampiskrp/ xkr_eta, xkr_zeta, xkr_ncontent, xkr_mass_min, xkr_mass_max
#endif
      NAMELIST/nampisdmp/ ln_pisdmp, nn_pisdmp
      NAMELIST/nampismass/ ln_check_mass
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampisbio in reference namelist : Pisces variables
      READ  ( numnatp_ref, nampisbio, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisbio in configuration namelist : Pisces variables
      READ  ( numnatp_cfg, nampisbio, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisbio )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' Namelist : nampisbio'
         WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
         WRITE(numout,*) '    POC sinking speed                         wsbio     =', wsbio
         WRITE(numout,*) '    half saturation constant for mortality    xkmort    =', xkmort
         WRITE(numout,*) '    Fe/C in zooplankton                       ferat3    =', ferat3
         WRITE(numout,*) '    Big particles sinking speed               wsbio2    =', wsbio2
         WRITE(numout,*) '    Maximum number of iterations for POC      niter1max =', niter1max
         WRITE(numout,*) '    Maximum number of iterations for GOC      niter2max =', niter2max
      ENDIF

#if defined key_kriest

      !                               ! nampiskrp : kriest parameters
      !                               ! -----------------------------
      REWIND( numnatp_ref )              ! Namelist nampiskrp in reference namelist : Pisces Kriest
      READ  ( numnatp_ref, nampiskrp, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiskrp in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampiskrp in configuration namelist : Pisces Kriest
      READ  ( numnatp_cfg, nampiskrp, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiskrp in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiskrp )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampiskrp'
         WRITE(numout,*) '    Sinking  exponent                        xkr_eta      = ', xkr_eta
         WRITE(numout,*) '    N content exponent                       xkr_zeta     = ', xkr_zeta
         WRITE(numout,*) '    N content factor                         xkr_ncontent = ', xkr_ncontent
         WRITE(numout,*) '    Minimum mass for Aggregates              xkr_mass_min = ', xkr_mass_min
         WRITE(numout,*) '    Maximum mass for Aggregates              xkr_mass_max = ', xkr_mass_max
         WRITE(numout,*)
     ENDIF


     ! Computation of some variables
     xkr_massp = xkr_ncontent * 7.625 * xkr_mass_min**xkr_zeta

#endif

      REWIND( numnatp_ref )              ! Namelist nampisdmp in reference namelist : Pisces damping
      READ  ( numnatp_ref, nampisdmp, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisdmp in configuration namelist : Pisces damping
      READ  ( numnatp_cfg, nampisdmp, IOSTAT = ios, ERR = 906 )
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisdmp )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampisdmp'
         WRITE(numout,*) '    Relaxation of tracer to glodap mean value             ln_pisdmp      =', ln_pisdmp
         WRITE(numout,*) '    Frequency of Relaxation                               nn_pisdmp      =', nn_pisdmp
         WRITE(numout,*) ' '
      ENDIF

      REWIND( numnatp_ref )              ! Namelist nampismass in reference namelist : Pisces mass conservation check
      READ  ( numnatp_ref, nampismass, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismass in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismass in configuration namelist : Pisces mass conservation check 
      READ  ( numnatp_cfg, nampismass, IOSTAT = ios, ERR = 908 )
908   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismass in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismass )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameter for mass conservation checking'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Flag to check mass conservation of NO3/Si/TALK ln_check_mass = ', ln_check_mass
      ENDIF

   END SUBROUTINE p4z_sms_init

   SUBROUTINE p4z_ph_ini
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_ini_ph  ***
      !!
      !!  ** Purpose : Initialization of chemical variables of the carbon cycle
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!---------------------------------------------------------------------

      ! Set PH from  total alkalinity, borat (???), akb3 (???) and ak23 (???)
      ! --------------------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztmas   = tmask(ji,jj,jk)
               ztmas1  = 1. - tmask(ji,jj,jk)
               zcaralk = trb(ji,jj,jk,jptal) - borat(ji,jj,jk) / (  1. + 1.E-8 / ( rtrn + akb3(ji,jj,jk) )  )
               zco3    = ( zcaralk - trb(ji,jj,jk,jpdic) ) * ztmas + 0.5e-3 * ztmas1
               zbicarb = ( 2. * trb(ji,jj,jk,jpdic) - zcaralk )
               hi(ji,jj,jk) = ( ak23(ji,jj,jk) * zbicarb / zco3 ) * ztmas + 1.e-9 * ztmas1
            END DO
         END DO
     END DO
     !
   END SUBROUTINE p4z_ph_ini

   SUBROUTINE p4z_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_rst  ***
      !!
      !!  ** Purpose : Read or write variables in restart file:
      !!
      !!  WRITE(READ) mode:
      !!       kt        : number of time step since the begining of the experiment at the
      !!                   end of the current(previous) run
      !!---------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!---------------------------------------------------------------------

      IF( TRIM(cdrw) == 'READ' ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' p4z_rst : Read specific variables from pisces model '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
         ! 
         IF( iom_varid( numrtr, 'PH', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_autoglo, 'PH' , hi(:,:,:)  )
         ELSE
!            hi(:,:,:) = 1.e-9 
            CALL p4z_ph_ini
         ENDIF
         CALL iom_get( numrtr, jpdom_autoglo, 'Silicalim', xksi(:,:) )
         IF( iom_varid( numrtr, 'Silicamax', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_autoglo, 'Silicamax' , xksimax(:,:)  )
         ELSE
            xksimax(:,:) = xksi(:,:)
         ENDIF
         !
         IF( iom_varid( numrtr, 'tcflxcum', ldstop = .FALSE. ) > 0 ) THEN  ! cumulative total flux of carbon
            CALL iom_get( numrtr, 'tcflxcum' , t_oce_co2_flx_cum  )
         ELSE
            t_oce_co2_flx_cum = 0._wp
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         IF( kt == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'p4z_rst : write pisces restart file  kt =', kt
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         CALL iom_rstput( kt, nitrst, numrtw, 'PH', hi(:,:,:) )
         CALL iom_rstput( kt, nitrst, numrtw, 'Silicalim', xksi(:,:) )
         CALL iom_rstput( kt, nitrst, numrtw, 'Silicamax', xksimax(:,:) )
         CALL iom_rstput( kt, nitrst, numrtw, 'tcflxcum', t_oce_co2_flx_cum )
      ENDIF
      !
   END SUBROUTINE p4z_rst

   SUBROUTINE p4z_dmp( kt )
      !!----------------------------------------------------------------------
      !!                    ***  p4z_dmp  ***
      !!
      !! ** purpose  : Relaxation of some tracers
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in )  ::     kt ! time step
      !
      REAL(wp) ::  alkmean = 2426.     ! mean value of alkalinity ( Glodap ; for Goyet 2391. )
      REAL(wp) ::  po4mean = 2.165     ! mean value of phosphates
      REAL(wp) ::  no3mean = 30.90     ! mean value of nitrate
      REAL(wp) ::  silmean = 91.51     ! mean value of silicate
      !
      REAL(wp) :: zarea, zalksumn, zpo4sumn, zno3sumn, zsilsumn
      REAL(wp) :: zalksumb, zpo4sumb, zno3sumb, zsilsumb
      !!---------------------------------------------------------------------


      IF(lwp)  WRITE(numout,*)
      IF(lwp)  WRITE(numout,*) ' p4z_dmp : Restoring of nutrients at time-step kt = ', kt
      IF(lwp)  WRITE(numout,*)

      IF( cp_cfg == "orca" .AND. .NOT. lk_c1d ) THEN      ! ORCA configuration (not 1D) !
         !                                                    ! --------------------------- !
         ! set total alkalinity, phosphate, nitrate & silicate
         zarea          = 1._wp / glob_sum( cvol(:,:,:) ) * 1e6              

         zalksumn = glob_sum( trn(:,:,:,jptal) * cvol(:,:,:)  ) * zarea
         zpo4sumn = glob_sum( trn(:,:,:,jppo4) * cvol(:,:,:)  ) * zarea * po4r
         zno3sumn = glob_sum( trn(:,:,:,jpno3) * cvol(:,:,:)  ) * zarea * rno3
         zsilsumn = glob_sum( trn(:,:,:,jpsil) * cvol(:,:,:)  ) * zarea
 
         IF(lwp) WRITE(numout,*) '       TALKN mean : ', zalksumn
         trn(:,:,:,jptal) = trn(:,:,:,jptal) * alkmean / zalksumn

         IF(lwp) WRITE(numout,*) '       PO4N  mean : ', zpo4sumn
         trn(:,:,:,jppo4) = trn(:,:,:,jppo4) * po4mean / zpo4sumn

         IF(lwp) WRITE(numout,*) '       NO3N  mean : ', zno3sumn
         trn(:,:,:,jpno3) = trn(:,:,:,jpno3) * no3mean / zno3sumn

         IF(lwp) WRITE(numout,*) '       SiO3N mean : ', zsilsumn
         trn(:,:,:,jpsil) = MIN( 400.e-6,trn(:,:,:,jpsil) * silmean / zsilsumn )
         !
         !
         IF( .NOT. ln_top_euler ) THEN
            zalksumb = glob_sum( trb(:,:,:,jptal) * cvol(:,:,:)  ) * zarea
            zpo4sumb = glob_sum( trb(:,:,:,jppo4) * cvol(:,:,:)  ) * zarea * po4r
            zno3sumb = glob_sum( trb(:,:,:,jpno3) * cvol(:,:,:)  ) * zarea * rno3
            zsilsumb = glob_sum( trb(:,:,:,jpsil) * cvol(:,:,:)  ) * zarea
 
            IF(lwp) WRITE(numout,*) ' '
            IF(lwp) WRITE(numout,*) '       TALKB mean : ', zalksumb
            trb(:,:,:,jptal) = trb(:,:,:,jptal) * alkmean / zalksumb

            IF(lwp) WRITE(numout,*) '       PO4B  mean : ', zpo4sumb
            trb(:,:,:,jppo4) = trb(:,:,:,jppo4) * po4mean / zpo4sumb

            IF(lwp) WRITE(numout,*) '       NO3B  mean : ', zno3sumb
            trb(:,:,:,jpno3) = trb(:,:,:,jpno3) * no3mean / zno3sumb

            IF(lwp) WRITE(numout,*) '       SiO3B mean : ', zsilsumb
            trb(:,:,:,jpsil) = MIN( 400.e-6,trb(:,:,:,jpsil) * silmean / zsilsumb )
        ENDIF
        !
      ENDIF
        !
   END SUBROUTINE p4z_dmp


   SUBROUTINE p4z_chk_mass( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_chk_mass  ***
      !!
      !! ** Purpose :  Mass conservation check 
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      REAL(wp)             ::  zrdenittot, zsdenittot, znitrpottot
      CHARACTER(LEN=100)   ::   cltxt
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zvol
      INTEGER :: jk
      !!----------------------------------------------------------------------

      !
      !!---------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN 
         IF( ln_check_mass .AND. lwp) THEN      !   Open budget file of NO3, ALK, Si, Fer
            CALL ctl_opn( numco2, 'carbon.budget'  , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numnut, 'nutrient.budget', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numnit, 'nitrogen.budget', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            xfact1 = rfact2r * 12. / 1.e15 * ryyss    ! conversion molC/kt --> PgC/yr
            xfact2 = 1.e+3 * rno3 * 14. / 1.e12 * ryyss   ! conversion molC/l/s ----> TgN/m3/yr
            xfact3 = 1.e+3 * rfact2r * rno3   ! conversion molC/l/kt ----> molN/m3/s
            cltxt='time-step   Alkalinity        Nitrate        Phosphorus         Silicate           Iron'
            IF( lwp ) WRITE(numnut,*)  TRIM(cltxt)
            IF( lwp ) WRITE(numnut,*) 
         ENDIF
      ENDIF

      !
      IF( iom_use( "pno3tot" ) .OR. ( ln_check_mass .AND. kt == nitend )  ) THEN
         !   Compute the budget of NO3, ALK, Si, Fer
         no3budget = glob_sum( (   trn(:,:,:,jpno3) + trn(:,:,:,jpnh4)  &
            &                    + trn(:,:,:,jpphy) + trn(:,:,:,jpdia)  &
            &                    + trn(:,:,:,jpzoo) + trn(:,:,:,jpmes)  &
            &                    + trn(:,:,:,jppoc)                     &
#if ! defined key_kriest
            &                    + trn(:,:,:,jpgoc)                     &
#endif
            &                    + trn(:,:,:,jpdoc)                     ) * cvol(:,:,:)  )
         !
         no3budget = no3budget / areatot
         CALL iom_put( "pno3tot", no3budget )
      ENDIF
      !
      IF( iom_use( "ppo4tot" ) .OR. ( ln_check_mass .AND. kt == nitend )  ) THEN
         po4budget = glob_sum( (   trn(:,:,:,jppo4)                     &
            &                    + trn(:,:,:,jpphy) + trn(:,:,:,jpdia)  &
            &                    + trn(:,:,:,jpzoo) + trn(:,:,:,jpmes)  &
            &                    + trn(:,:,:,jppoc)                     &
#if ! defined key_kriest
            &                    + trn(:,:,:,jpgoc)                     &
#endif
            &                    + trn(:,:,:,jpdoc)                     ) * cvol(:,:,:)  )
         po4budget = po4budget / areatot
         CALL iom_put( "ppo4tot", po4budget )
      ENDIF
      !
      IF( iom_use( "psiltot" ) .OR. ( ln_check_mass .AND. kt == nitend )  ) THEN
         silbudget = glob_sum( (   trn(:,:,:,jpsil) + trn(:,:,:,jpgsi)  &
            &                    + trn(:,:,:,jpdsi)                     ) * cvol(:,:,:)  )
         !
         silbudget = silbudget / areatot
         CALL iom_put( "psiltot", silbudget )
      ENDIF
      !
      IF( iom_use( "palktot" ) .OR. ( ln_check_mass .AND. kt == nitend )  ) THEN
         alkbudget = glob_sum( (   trn(:,:,:,jpno3) * rno3              &
            &                    + trn(:,:,:,jptal)                     &
            &                    + trn(:,:,:,jpcal) * 2.                ) * cvol(:,:,:)  )
         !
         alkbudget = alkbudget / areatot
         CALL iom_put( "palktot", alkbudget )
      ENDIF
      !
      IF( iom_use( "pfertot" ) .OR. ( ln_check_mass .AND. kt == nitend )  ) THEN
         ferbudget = glob_sum( (   trn(:,:,:,jpfer) + trn(:,:,:,jpnfe)  &
            &                    + trn(:,:,:,jpdfe)                     &
#if ! defined key_kriest
            &                    + trn(:,:,:,jpbfe)                     &
#endif
            &                    + trn(:,:,:,jpsfe)                     &
            &                    + trn(:,:,:,jpzoo) * ferat3            &
            &                    + trn(:,:,:,jpmes) * ferat3            ) * cvol(:,:,:)  )
         !
         ferbudget = ferbudget / areatot
         CALL iom_put( "pfertot", ferbudget )
      ENDIF
      !

      ! Global budget of N SMS : denitrification in the water column and in the sediment
      !                          nitrogen fixation by the diazotrophs
      ! --------------------------------------------------------------------------------
      IF( iom_use( "tnfix" ) .OR.  ( ln_check_mass .AND. kt == nitend )  ) THEN
         znitrpottot  = glob_sum ( nitrpot(:,:,:) * nitrfix * cvol(:,:,:) )
         CALL iom_put( "tnfix"  , znitrpottot * 1.e+3 * rno3 )  ! Global  nitrogen fixation molC/l  to molN/m3 
      ENDIF
      !
      IF( iom_use( "tdenit" ) .OR.  ( ln_check_mass .AND. kt == nitend )  ) THEN
         zrdenittot   = glob_sum ( denitr(:,:,:) * rdenit * xnegtr(:,:,:) * cvol(:,:,:) )
         CALL iom_put( "tdenit"  , zrdenittot * 1.e+3 * rno3 )  ! Total denitrification molC/l to molN/m3 
      ENDIF
      !
      IF( iom_use( "Sdenit" ) .OR.  ( ln_check_mass .AND. kt == nitend )  ) THEN
         zsdenittot   = glob_sum ( sdenit(:,:) * e1e2t(:,:) )
         CALL iom_put( "Sdenit", sdenit(:,:) * xfact3 * tmask(:,:,1) )  ! Nitrate reduction in the sediments
      ENDIF

      IF( ln_check_mass .AND. kt == nitend ) THEN   ! Compute the budget of NO3, ALK, Si, Fer
         t_atm_co2_flx  = t_atm_co2_flx / glob_sum( e1e2t(:,:) )
         t_oce_co2_flx  = t_oce_co2_flx         * xfact1 * (-1 )
         tpp            = tpp           * 1000. * xfact1
         t_oce_co2_exp  = t_oce_co2_exp * 1000. * xfact1
         IF( lwp ) WRITE(numco2,9000) ndastp, t_atm_co2_flx, t_oce_co2_flx, tpp, t_oce_co2_exp
         IF( lwp ) WRITE(numnut,9100) ndastp, alkbudget        * 1.e+06, &
             &                                no3budget * rno3 * 1.e+06, &
             &                                po4budget * po4r * 1.e+06, &
             &                                silbudget        * 1.e+06, &
             &                                ferbudget        * 1.e+09
         !
         IF( lwp ) WRITE(numnit,9200) ndastp, znitrpottot * xfact2  , &
         &                             zrdenittot  * xfact2  , &
         &                             zsdenittot  * xfact2

      ENDIF
      !
 9000  FORMAT(i8,f10.5,e18.10,f10.5,f10.5)
 9100  FORMAT(i8,5e18.10)
 9200  FORMAT(i8,3f10.5)

       !
   END SUBROUTINE p4z_chk_mass

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_sms: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_sms
#endif 

   !!======================================================================
END MODULE p4zsms 
