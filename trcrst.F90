MODULE trcrst
   !!======================================================================
   !!                         ***  MODULE trcrst  ***
   !! TOP :   Manage the passive tracer restart
   !!======================================================================
   !! History :    -   !  1991-03  ()  original code
   !!             1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!              -   !  2005-10 (C. Ethe) print control
   !!             2.0  !  2005-10 (C. Ethe, G. Madec) revised architecture
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   trc_rst :   Restart for passive tracer
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_rst_opn    : open  restart file
   !!   trc_rst_read   : read  restart file
   !!   trc_rst_wri    : write restart file
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trc
   USE trcnam_trp
   USE iom
   USE ioipsl, ONLY : ju2ymds    ! for calendar
   USE daymod
   !! AXY (05/11/13): need these for MEDUSA to input/output benthic reservoirs
   USE par_medusa
   USE sms_medusa
   USE trcsms_medusa
   !!
#if defined key_idtra
   USE trcsms_idtra
#endif
   !!
#if defined key_cfc
   USE trcsms_cfc
#endif
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE sbc_oce, ONLY: lk_oasis 
   USE oce,     ONLY: CO2Flux_out_cpl, DMS_out_cpl, chloro_out_cpl  !! Coupling variable
   USE trcstat

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_rst_opn       ! called by ???
   PUBLIC   trc_rst_read      ! called by ???
   PUBLIC   trc_rst_wri       ! called by ???
   PUBLIC   trc_rst_cal
   PUBLIC   trc_rst_stat
#if defined key_medusa && defined key_roam
   PUBLIC   trc_rst_conserve
#endif

   !! * Substitutions
#  include "top_substitute.h90"
   
CONTAINS
   
   SUBROUTINE trc_rst_opn( kt )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_opn  ***
      !!
      !! ** purpose  :   output of sea-trc variable in a netcdf file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! number of iteration
      INTEGER             ::   iyear, imonth, iday
      REAL (wp)           ::   zsec
      REAL (wp)           ::   zfjulday
      !
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! trc output restart file name
      CHARACTER(LEN=256)  ::   clpath   ! full path to ocean output restart file
      !!----------------------------------------------------------------------
      !
      IF( lk_offline ) THEN
         IF( kt == nittrc000 ) THEN
            lrst_trc = .FALSE.
            IF( ln_rst_list ) THEN
               nrst_lst = 1
               nitrst = nstocklist( nrst_lst )
            ELSE
               nitrst = nitend
            ENDIF
         ENDIF

         IF( .NOT. ln_rst_list .AND. MOD( kt - 1, nstock ) == 0 ) THEN
            ! we use kt - 1 and not kt - nittrc000 to keep the same periodicity from the beginning of the experiment
            nitrst = kt + nstock - 1                  ! define the next value of nitrst for restart writing
            IF( nitrst > nitend )   nitrst = nitend   ! make sure we write a restart at the end of the run
         ENDIF
      ELSE
         IF( kt == nittrc000 ) lrst_trc = .FALSE.
      ENDIF

      ! to get better performances with NetCDF format:
      ! we open and define the tracer restart file one tracer time step before writing the data (-> at nitrst - 2*nn_dttrc + 1)
      ! except if we write tracer restart files every tracer time step or if a tracer restart file was writen at nitend - 2*nn_dttrc + 1
      IF( kt == nitrst - 2*nn_dttrc .OR. nstock == nn_dttrc .OR. ( kt == nitend - nn_dttrc .AND. .NOT. lrst_trc ) ) THEN
         IF ( ln_rstdate ) THEN
            !! JPALM -- 22-12-2015 -- modif to get the good date on restart trc file name
            !!                     -- the condition to open the rst file is not the same than for the dynamic rst.
            !!                     -- here it - for an obscure reason - is open 2 time-step before the restart writing process
            !!                     instead of 1.
            !!                     -- i am not sure if someone forgot +1 in the if loop condition as
            !!                     it is writen in all comments nitrst -2*nn_dttrc + 1 and the condition is 
            !!                     nitrst - 2*nn_dttrc
            !!                     -- nevertheless we didn't wanted to broke something already working 
            !!                     and just adapted the part we added.
            !!                     -- So instead of calling ju2ymds( fjulday + (rdttra(1)) 
            !!                     we call ju2ymds( fjulday + (2*rdttra(1)) 
            !!--------------------------------------------------------------------      
            zfjulday = fjulday + (2*rdttra(1)) / rday
            IF( ABS(zfjulday - REAL(NINT(zfjulday),wp)) < 0.1 / rday )   zfjulday = REAL(NINT(zfjulday),wp)   ! avoid truncation error
            CALL ju2ymds( zfjulday + (2*rdttra(1)) / rday, iyear, imonth, iday, zsec )
            WRITE(clkt, '(i4.4,2i2.2)') iyear, imonth, iday
         ELSE
            ! beware of the format used to write kt (default is i8.8, that should be large enough)
            IF( nitrst > 1.0e9 ) THEN   ;   WRITE(clkt,*       ) nitrst
            ELSE                        ;   WRITE(clkt,'(i8.8)') nitrst
            ENDIF
         ENDIF
         ! create the file
         IF(lwp) WRITE(numout,*)
         clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_trcrst_out)
         clpath = TRIM(cn_trcrst_outdir)
         IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
         IF(lwp) WRITE(numout,*) &
             '             open trc restart.output NetCDF file: ',TRIM(clpath)//clname
         CALL iom_open( TRIM(clpath)//TRIM(clname), numrtw, ldwrt = .TRUE., kiolib = jprstlib )
         lrst_trc = .TRUE.
      ENDIF
      !
   END SUBROUTINE trc_rst_opn

   SUBROUTINE trc_rst_read
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_opn  ***
      !!
      !! ** purpose  :   read passive tracer fields in restart files
      !!----------------------------------------------------------------------
      INTEGER  ::  jn, jl     
      !! AXY (05/11/13): temporary variables
      REAL(wp) ::    fq0,fq1,fq2

      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_rst_read : read data in the TOP restart file'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

      ! READ prognostic variables and computes diagnostic variable
      DO jn = 1, jptra
         CALL iom_get( numrtr, jpdom_autoglo, 'TRN'//ctrcnm(jn), trn(:,:,:,jn) )
         trn(:,:,:,jn) = trn(:,:,:,jn) * tmask(:,:,:)
      END DO

      DO jn = 1, jptra
         CALL iom_get( numrtr, jpdom_autoglo, 'TRB'//ctrcnm(jn), trb(:,:,:,jn) )
         trb(:,:,:,jn) = trb(:,:,:,jn) * tmask(:,:,:)
      END DO
      !
      !! AXY (09/06/14): the ARCHER version of MEDUSA-2 does not include an equivalent
      !!                 call to MEDUSA-2 at this point; this suggests that the FCM
      !!                 version of NEMO date significantly earlier than the current
      !!                 version

#if defined key_medusa
      !! AXY (13/01/12): check if the restart contains sediment fields;
      !!                 this is only relevant for simulations that include
      !!                 biogeochemistry and are restarted from earlier runs
      !!                 in which there was no sediment component
      !!
      IF( iom_varid( numrtr, 'B_SED_N', ldstop = .FALSE. ) > 0 ) THEN
         !! YES; in which case read them
         !!
         IF(lwp) WRITE(numout,*) ' MEDUSA sediment fields present - reading in ...'
         CALL iom_get( numrtr, jpdom_autoglo, 'B_SED_N',  zb_sed_n(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_SED_N',  zn_sed_n(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_SED_FE', zb_sed_fe(:,:) )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_SED_FE', zn_sed_fe(:,:) )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_SED_SI', zb_sed_si(:,:) )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_SED_SI', zn_sed_si(:,:) )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_SED_C',  zb_sed_c(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_SED_C',  zn_sed_c(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_SED_CA', zb_sed_ca(:,:) )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_SED_CA', zn_sed_ca(:,:) )
      ELSE
         !! NO; in which case set them to zero
         !!
         IF(lwp) WRITE(numout,*) ' MEDUSA sediment fields absent - setting to zero ...'
         zb_sed_n(:,:)  = 0.0   !! organic N
         zn_sed_n(:,:)  = 0.0
         zb_sed_fe(:,:) = 0.0   !! organic Fe
         zn_sed_fe(:,:) = 0.0
         zb_sed_si(:,:) = 0.0   !! inorganic Si
         zn_sed_si(:,:) = 0.0
         zb_sed_c(:,:)  = 0.0   !! organic C
         zn_sed_c(:,:)  = 0.0
         zb_sed_ca(:,:) = 0.0   !! inorganic C
         zn_sed_ca(:,:) = 0.0
      ENDIF
      !!
      !! calculate stats on these fields
      IF(lwp) WRITE(numout,*) ' MEDUSA sediment field stats (min, max, sum) ...'
      call trc_rst_dia_stat(zn_sed_n(:,:), 'Sediment  N')
      call trc_rst_dia_stat(zn_sed_fe(:,:), 'Sediment Fe')
      call trc_rst_dia_stat(zn_sed_si(:,:), 'Sediment Si')
      call trc_rst_dia_stat(zn_sed_c(:,:), 'Sediment C')
      call trc_rst_dia_stat(zn_sed_ca(:,:), 'Sediment Ca')
      !!
      !! AXY (07/07/15): read in temporally averaged fields for DMS
      !!                 calculations
      !!
      IF( iom_varid( numrtr, 'B_DMS_CHN', ldstop = .FALSE. ) > 0 ) THEN
         !! YES; in which case read them
         !!
         IF(lwp) WRITE(numout,*) ' MEDUSA averaged properties for DMS present - reading in ...'
         CALL iom_get( numrtr, jpdom_autoglo, 'B_DMS_CHN',  zb_dms_chn(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_DMS_CHN',  zn_dms_chn(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_DMS_CHD',  zb_dms_chd(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_DMS_CHD',  zn_dms_chd(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_DMS_MLD',  zb_dms_mld(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_DMS_MLD',  zn_dms_mld(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_DMS_QSR',  zb_dms_qsr(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_DMS_QSR',  zn_dms_qsr(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'B_DMS_DIN',  zb_dms_din(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_DMS_DIN',  zn_dms_din(:,:)  )
      ELSE
         !! NO; in which case set them to zero
         !!
         IF(lwp) WRITE(numout,*) ' MEDUSA averaged properties for DMS absent - setting to zero ...'
         zb_dms_chn(:,:)  = 0.0   !! CHN
         zn_dms_chn(:,:)  = 0.0
         zb_dms_chd(:,:)  = 0.0   !! CHD
         zn_dms_chd(:,:)  = 0.0
         zb_dms_mld(:,:)  = 0.0   !! MLD
         zn_dms_mld(:,:)  = 0.0
         zb_dms_qsr(:,:)  = 0.0   !! QSR
         zn_dms_qsr(:,:)  = 0.0
         zb_dms_din(:,:)  = 0.0   !! DIN
         zn_dms_din(:,:)  = 0.0
      ENDIF
      !!  
      !! JPALM 14-06-2016 -- add CO2 flux and DMS surf through the restart
      !!                  -- needed for the coupling with atm
      IF( iom_varid( numrtr, 'N_DMS_srf', ldstop = .FALSE. ) > 0 ) THEN
         IF(lwp) WRITE(numout,*) 'DMS surf concentration - reading in ...'
         CALL iom_get( numrtr, jpdom_autoglo, 'B_DMS_srf',  zb_dms_srf(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_DMS_srf',  zn_dms_srf(:,:)  )
      ELSE
         IF(lwp) WRITE(numout,*) 'DMS surf concentration - setting to zero ...'
         zb_dms_srf(:,:)  = 0.0   !! DMS
         zn_dms_srf(:,:)  = 0.0
      ENDIF
      IF (lk_oasis) THEN
         DMS_out_cpl(:,:) = zn_dms_srf(:,:)        !! Coupling variable
      END IF
      !!
      IF( iom_varid( numrtr, 'B_CO2_flx', ldstop = .FALSE. ) > 0 ) THEN
         IF(lwp) WRITE(numout,*) 'CO2 air-sea flux - reading in ...'
         CALL iom_get( numrtr, jpdom_autoglo, 'B_CO2_flx',  zb_co2_flx(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'N_CO2_flx',  zn_co2_flx(:,:)  )
      ELSE
         IF(lwp) WRITE(numout,*) 'CO2 air-sea flux - setting to zero ...'
         zb_co2_flx(:,:)  = 0.0   !! CO2 flx
         zn_co2_flx(:,:)  = 0.0
      ENDIF
      IF (lk_oasis) THEN
         CO2Flux_out_cpl(:,:) =  zn_co2_flx(:,:)   !! Coupling variable
      END IF
      !!
      !! JPALM 02-06-2017 -- in complement to DMS surf 
      !!                  -- the atm model needs surf Chl 
      !!                     as proxy of org matter from the ocean
      !!                  -- needed for the coupling with atm
      !!       07-12-2017 -- To make things cleaner, we want to store an  
      !!                     unscaled Chl field in the restart and only 
      !!                     scale it when reading it in.

      IF( iom_varid( numrtr, 'N_CHL_srf', ldstop = .FALSE. ) > 0 ) THEN
         IF(lwp) WRITE(numout,*) 'Chl cpl concentration - reading in ... - scale by ', scl_chl
         CALL iom_get( numrtr, jpdom_autoglo, 'N_CHL_srf',  zn_chl_srf(:,:)  )
      ELSE
         IF(lwp) WRITE(numout,*) 'set Chl coupled concentration - scaled by ', scl_chl
         zn_chl_srf(:,:)  = MAX( 0.0, (trn(:,:,1,jpchn) + trn(:,:,1,jpchd)) * 1.E-6 )
      ENDIF
      IF (lk_oasis) THEN
         chloro_out_cpl(:,:) = zn_chl_srf(:,:) * scl_chl        !! Coupling variable
      END IF
      !!
      !! calculate stats on these fields
      IF(lwp) WRITE(numout,*) ' MEDUSA averaged properties for DMS stats (min, max, sum) ...'
      call trc_rst_dia_stat(zn_dms_chn(:,:), 'DMS, CHN')
      call trc_rst_dia_stat(zn_dms_chd(:,:), 'DMS, CHD')
      call trc_rst_dia_stat(zn_dms_mld(:,:), 'DMS, MLD')
      call trc_rst_dia_stat(zn_dms_qsr(:,:), 'DMS, QSR')
      call trc_rst_dia_stat(zn_dms_din(:,:), 'DMS, DIN')
      call trc_rst_dia_stat(zn_dms_srf(:,:), 'DMS surf')
      call trc_rst_dia_stat(zn_co2_flx(:,:), 'CO2 flux')
      IF (lk_oasis) THEN
         call trc_rst_dia_stat(chloro_out_cpl(:,:), 'CHL  cpl')
      END IF
      !!  
      !! JPALM 14-06-2016 -- add Carbonate chenistry variables through the restart
      !!                  -- needed for monthly call of carb-chem routine and better reproducibility
# if defined key_roam
      IF( iom_varid( numrtr, 'pH_3D', ldstop = .FALSE. ) > 0 ) THEN
         IF(lwp) WRITE(numout,*) 'Carbonate chem variable - reading in ...'
         CALL iom_get( numrtr, jpdom_autoglo, 'pH_3D'   ,  f3_pH(:,:,:)     )
         CALL iom_get( numrtr, jpdom_autoglo, 'h2CO3_3D',  f3_h2co3(:,:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'hCO3_3D' ,  f3_hco3(:,:,:)   )
         CALL iom_get( numrtr, jpdom_autoglo, 'CO3_3D'  ,  f3_co3(:,:,:)    )
         CALL iom_get( numrtr, jpdom_autoglo, 'omcal_3D',  f3_omcal(:,:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'omarg_3D',  f3_omarg(:,:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'CCD_CAL' ,  f2_ccd_cal(:,:)  )
         CALL iom_get( numrtr, jpdom_autoglo, 'CCD_ARG' ,  f2_ccd_arg(:,:)  )
         !!
         IF(lwp) WRITE(numout,*) ' MEDUSA averaged Carb-chem stats (min, max, sum) ...'
      call trc_rst_dia_stat( f3_pH(:,:,1)   ,'pH 3D surf')
      call trc_rst_dia_stat( f3_h2co3(:,:,1),'h2CO3 3D surf')
      call trc_rst_dia_stat( f3_hco3(:,:,1) ,'hCO3 3D surf' )
      call trc_rst_dia_stat( f3_co3(:,:,1)  ,'CO3 3D surf' )
      call trc_rst_dia_stat( f3_omcal(:,:,1),'omcal 3D surf')
      call trc_rst_dia_stat( f3_omarg(:,:,1),'omarg 3D surf')
      call trc_rst_dia_stat( f2_ccd_cal(:,:),'CCD_CAL')
      call trc_rst_dia_stat( f2_ccd_arg(:,:),'CCD_ARG')

      ELSE
         IF(lwp) WRITE(numout,*) 'WARNING : No Carbonate-chem variable in the restart.... '
         IF(lwp) WRITE(numout,*) 'Is not a problem if start a month, but may be very problematic if not '
         IF(lwp) WRITE(numout,*) 'Check if   mod(kt*rdt,2592000) == rdt' 
         IF(lwp) WRITE(numout,*) 'Or don t start from uncomplete restart...' 
      ENDIF
# endif


#endif
      !
#if defined key_idtra
      !! JPALM -- 05-01-2016 -- mv idtra and CFC restart reading and 
      !!                        writting here undre their key.
      !!                        problems in CFC restart, maybe because of this...
      !!                        and pb in idtra diag or diad-restart writing.
      !!----------------------------------------------------------------------
      IF( iom_varid( numrtr, 'qint_IDTRA', ldstop = .FALSE. ) > 0 ) THEN
         !! YES; in which case read them
         !!
         IF(lwp) WRITE(numout,*) ' IDTRA averaged properties present - reading in ...'
         CALL iom_get( numrtr, jpdom_autoglo, 'qint_IDTRA',  qint_idtra(:,:,1)  )
      ELSE
         !! NO; in which case set them to zero
         !!
         IF(lwp) WRITE(numout,*) ' IDTRA averaged properties absent - setting to zero ...'
         qint_idtra(:,:,1)  = 0.0   !! CHN
      ENDIF
      !!
      !! calculate stats on these fields
      IF(lwp) WRITE(numout,*) ' IDTRA averaged properties stats (min, max, sum) ...'
      call trc_rst_dia_stat(qint_idtra(:,:,1), 'qint_IDTRA')
#endif
      !
#if defined key_cfc
      DO jl = 1, jp_cfc
         jn = jp_cfc0 + jl - 1
         IF( iom_varid( numrtr, 'qint_'//ctrcnm(jn), ldstop = .FALSE. ) > 0 ) THEN
            !! YES; in which case read them
            !!
            IF(lwp) WRITE(numout,*) ' CFC averaged properties present - reading in ...'
            CALL iom_get( numrtr, jpdom_autoglo, 'qint_'//ctrcnm(jn), qint_cfc(:,:,jl) )
         ELSE
            !! NO; in which case set them to zero
            !!
            IF(lwp) WRITE(numout,*) ' CFC averaged properties absent - setting to zero ...'
            qint_cfc(:,:,jl)  = 0.0   !! CHN
         ENDIF
         !!
         !! calculate stats on these fields
         IF(lwp) WRITE(numout,*) ' CFC averaged properties stats (min, max, sum) ...'
         call trc_rst_dia_stat(qint_cfc(:,:,jl), 'qint_'//ctrcnm(jn))
      END DO
#endif
      !
   END SUBROUTINE trc_rst_read

   SUBROUTINE trc_rst_wri( kt )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_wri  ***
      !!
      !! ** purpose  :   write passive tracer fields in restart files
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index
      !!
      INTEGER  :: jn, jl
      REAL(wp) :: zarak0
      !! AXY (05/11/13): temporary variables
      REAL(wp) ::    fq0,fq1,fq2
      !!----------------------------------------------------------------------
      !
      CALL iom_rstput( kt, nitrst, numrtw, 'rdttrc1', rdttrc(1) )   ! surface passive tracer time step
      ! prognostic variables 
      ! -------------------- 
      DO jn = 1, jptra
         CALL iom_rstput( kt, nitrst, numrtw, 'TRN'//ctrcnm(jn), trn(:,:,:,jn) )
      END DO

      DO jn = 1, jptra
         CALL iom_rstput( kt, nitrst, numrtw, 'TRB'//ctrcnm(jn), trb(:,:,:,jn) )
      END DO

      !! AXY (09/06/14): the ARCHER version of MEDUSA-2 does not include an equivalent
      !!                 call to MEDUSA-2 at this point; this suggests that the FCM
      !!                 version of NEMO date significantly earlier than the current
      !!                 version

#if defined key_medusa
      !! AXY (13/01/12): write out "before" and "now" state of seafloor
      !!                 sediment pools into restart; this happens
      !!                 whether or not the pools are to be used by
      !!                 MEDUSA (which is controlled by a switch in the
      !!                 namelist_top file)
      !!
      IF(lwp) WRITE(numout,*) ' MEDUSA sediment fields - writing out ...'
      CALL iom_rstput( kt, nitrst, numrtw, 'B_SED_N',  zb_sed_n(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_SED_N',  zn_sed_n(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_SED_FE', zb_sed_fe(:,:) )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_SED_FE', zn_sed_fe(:,:) )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_SED_SI', zb_sed_si(:,:) )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_SED_SI', zn_sed_si(:,:) )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_SED_C',  zb_sed_c(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_SED_C',  zn_sed_c(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_SED_CA', zb_sed_ca(:,:) )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_SED_CA', zn_sed_ca(:,:) )
      !!
      !! calculate stats on these fields
      IF(lwp) WRITE(numout,*) ' MEDUSA sediment field stats (min, max, sum) ...'
      call trc_rst_dia_stat(zn_sed_n(:,:), 'Sediment  N')
      call trc_rst_dia_stat(zn_sed_fe(:,:), 'Sediment Fe')
      call trc_rst_dia_stat(zn_sed_si(:,:), 'Sediment Si')
      call trc_rst_dia_stat(zn_sed_c(:,:), 'Sediment C')
      call trc_rst_dia_stat(zn_sed_ca(:,:), 'Sediment Ca')
      !!
      !! AXY (07/07/15): write out temporally averaged fields for DMS
      !!                 calculations
      !!
      IF(lwp) WRITE(numout,*) ' MEDUSA averaged properties for DMS - writing out ...'
      CALL iom_rstput( kt, nitrst, numrtw, 'B_DMS_CHN',  zb_dms_chn(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_DMS_CHN',  zn_dms_chn(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_DMS_CHD',  zb_dms_chd(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_DMS_CHD',  zn_dms_chd(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_DMS_MLD',  zb_dms_mld(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_DMS_MLD',  zn_dms_mld(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_DMS_QSR',  zb_dms_qsr(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_DMS_QSR',  zn_dms_qsr(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_DMS_DIN',  zb_dms_din(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_DMS_DIN',  zn_dms_din(:,:)  )
         !! JPALM 14-06-2016 -- add CO2 flux and DMS surf through the restart
         !!                  -- needed for the coupling with atm
      CALL iom_rstput( kt, nitrst, numrtw, 'B_DMS_srf',  zb_dms_srf(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_DMS_srf',  zn_dms_srf(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'B_CO2_flx',  zb_co2_flx(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'N_CO2_flx',  zn_co2_flx(:,:)  )
      !! JPALM 07-12-2017 -- To make things cleaner, we want to store an  
      !!                     unscaled Chl field in the restart and only 
      !!                     scale it when reading it in.
      CALL iom_rstput( kt, nitrst, numrtw, 'N_CHL_srf',  zn_chl_srf(:,:)  )
      !!
      !! calculate stats on these fields
      IF(lwp) WRITE(numout,*) ' MEDUSA averaged properties for DMS stats (min, max, sum) ...'
      call trc_rst_dia_stat(zn_dms_chn(:,:), 'DMS, CHN')
      call trc_rst_dia_stat(zn_dms_chd(:,:), 'DMS, CHD')
      call trc_rst_dia_stat(zn_dms_mld(:,:), 'DMS, MLD')
      call trc_rst_dia_stat(zn_dms_qsr(:,:), 'DMS, QSR')
      call trc_rst_dia_stat(zn_dms_din(:,:), 'DMS, DIN')
      call trc_rst_dia_stat(zn_dms_srf(:,:), 'DMS surf')
      call trc_rst_dia_stat(zn_co2_flx(:,:), 'CO2 flux')
      call trc_rst_dia_stat(zn_chl_srf(:,:), 'unscaled CHL cpl')
      !!
      IF(lwp) WRITE(numout,*) ' MEDUSA averaged prop. for dust and iron dep.'
      call trc_rst_dia_stat(dust(:,:), 'Dust dep')
      call trc_rst_dia_stat(zirondep(:,:), 'Iron dep')
      !! 
      !!  
      !! JPALM 14-06-2016 -- add Carbonate chenistry variables through the restart
      !!                  -- needed for monthly call of carb-chem routine and better reproducibility
# if defined key_roam
      IF(lwp) WRITE(numout,*) 'Carbonate chem variable - writing out ...'
      CALL iom_rstput( kt, nitrst, numrtw, 'pH_3D'   ,  f3_pH(:,:,:)     )
      CALL iom_rstput( kt, nitrst, numrtw, 'h2CO3_3D',  f3_h2co3(:,:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'hCO3_3D' ,  f3_hco3(:,:,:)   )
      CALL iom_rstput( kt, nitrst, numrtw, 'CO3_3D'  ,  f3_co3(:,:,:)    )
      CALL iom_rstput( kt, nitrst, numrtw, 'omcal_3D',  f3_omcal(:,:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'omarg_3D',  f3_omarg(:,:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'CCD_CAL' ,  f2_ccd_cal(:,:)  )
      CALL iom_rstput( kt, nitrst, numrtw, 'CCD_ARG' ,  f2_ccd_arg(:,:)  )
      !!
      IF(lwp) WRITE(numout,*) ' MEDUSA averaged Carb-chem stats (min, max, sum) ...'
      call trc_rst_dia_stat( f3_pH(:,:,1)   ,'pH 3D surf')
      call trc_rst_dia_stat( f3_h2co3(:,:,1),'h2CO3 3D surf')
      call trc_rst_dia_stat( f3_hco3(:,:,1) ,'hCO3 3D surf' )
      call trc_rst_dia_stat( f3_co3(:,:,1)  ,'CO3 3D surf' )
      call trc_rst_dia_stat( f3_omcal(:,:,1),'omcal 3D surf')
      call trc_rst_dia_stat( f3_omarg(:,:,1),'omarg 3D surf')
      call trc_rst_dia_stat( f2_ccd_cal(:,:),'CCD_CAL')
      call trc_rst_dia_stat( f2_ccd_arg(:,:),'CCD_ARG')
      !!
# endif
!!
#endif
      !
#if defined key_idtra
      !! JPALM -- 05-01-2016 -- mv idtra and CFC restart reading and 
      !!                        writting here undre their key.
      !!                        problems in CFC restart, maybe because of this...
      !!                        and pb in idtra diag or diad-restart writing.
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*) ' IDTRA averaged properties - writing out ...'
      CALL iom_rstput( kt, nitrst, numrtw, 'qint_IDTRA',  qint_idtra(:,:,1) )
      !!
      !! calculate stats on these fields
      IF(lwp) WRITE(numout,*) ' IDTRA averaged properties stats (min, max, sum) ...'
      call trc_rst_dia_stat(qint_idtra(:,:,1), 'qint_IDTRA')
#endif
      !
#if defined key_cfc
      DO jl = 1, jp_cfc
         jn = jp_cfc0 + jl - 1
         IF(lwp) WRITE(numout,*) ' CFC averaged properties - writing out ...'
         CALL iom_rstput( kt, nitrst, numrtw, 'qint_'//ctrcnm(jn), qint_cfc(:,:,jl) )
         !!
         !! calculate stats on these fields
         IF(lwp) WRITE(numout,*) ' CFC averaged properties stats (min, max, sum) ...'
         call trc_rst_dia_stat(qint_cfc(:,:,jl), 'qint_'//ctrcnm(jn))
      END DO
#endif
      !

      IF( kt == nitrst ) THEN
          CALL trc_rst_stat            ! statistics
#if defined key_medusa && defined key_roam
          CALL trc_rst_conserve        ! conservation check
#endif
          CALL iom_close( numrtw )     ! close the restart file (only at last time step)
#if ! defined key_trdmxl_trc
          lrst_trc = .FALSE.
#endif
          IF( lk_offline .AND. ln_rst_list ) THEN
             nrst_lst = nrst_lst + 1
             nitrst = nstocklist( nrst_lst )
          ENDIF
      ENDIF
      !
   END SUBROUTINE trc_rst_wri 


   SUBROUTINE trc_rst_cal( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE trc_rst_cal  ***
      !!
      !!  ** Purpose : Read or write calendar in restart file:
      !!
      !!  WRITE(READ) mode:
      !!       kt        : number of time step since the begining of the experiment at the
      !!                   end of the current(previous) run
      !!       adatrj(0) : number of elapsed days since the begining of the experiment at the
      !!                   end of the current(previous) run (REAL -> keep fractions of day)
      !!       ndastp    : date at the end of the current(previous) run (coded as yyyymmdd integer)
      !!
      !!   According to namelist parameter nrstdt,
      !!       nn_rsttr = 0  no control on the date (nittrc000 is  arbitrary).
      !!       nn_rsttr = 1  we verify that nittrc000 is equal to the last
      !!                   time step of previous run + 1.
      !!       In both those options, the  exact duration of the experiment
      !!       since the beginning (cumulated duration of all previous restart runs)
      !!       is not stored in the restart and is assumed to be (nittrc000-1)*rdt.
      !!       This is valid is the time step has remained constant.
      !!
      !!       nn_rsttr = 2  the duration of the experiment in days (adatrj)
      !!                    has been stored in the restart file.
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      INTEGER  ::  jlibalt = jprstlib
      LOGICAL  ::  llok
      REAL(wp) ::  zkt, zrdttrc1
      REAL(wp) ::  zndastp

      ! Time domain : restart
      ! ---------------------

      IF( TRIM(cdrw) == 'READ' ) THEN

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_rst_cal : read the TOP restart file for calendar'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

         IF ( jprstlib == jprstdimg ) THEN
           ! eventually read netcdf file (monobloc)  for restarting on different number of processors
           ! if {cn_trcrst_in}.nc exists, then set jlibalt to jpnf90 
           INQUIRE( FILE = TRIM(cn_trcrst_indir)//'/'//TRIM(cn_trcrst_in)//'.nc', EXIST = llok )
           IF ( llok ) THEN ; jlibalt = jpnf90  ; ELSE ; jlibalt = jprstlib ; ENDIF
         ENDIF

         IF( ln_rsttr ) THEN
            CALL iom_open( TRIM(cn_trcrst_indir)//'/'//cn_trcrst_in, numrtr, kiolib = jlibalt )
            CALL iom_get ( numrtr, 'kt', zkt )   ! last time-step of previous run

            IF(lwp) THEN
               WRITE(numout,*) ' *** Info read in restart : '
               WRITE(numout,*) '   previous time-step                               : ', NINT( zkt )
               WRITE(numout,*) ' *** restart option'
               SELECT CASE ( nn_rsttr )
               CASE ( 0 )   ;   WRITE(numout,*) ' nn_rsttr = 0 : no control of nittrc000'
               CASE ( 1 )   ;   WRITE(numout,*) ' nn_rsttr = 1 : no control the date at nittrc000 (use ndate0 read in the namelist)'
               CASE ( 2 )   ;   WRITE(numout,*) ' nn_rsttr = 2 : calendar parameters read in restart'
               END SELECT
               WRITE(numout,*)
            ENDIF
            ! Control of date 
            IF( nittrc000  - NINT( zkt ) /= nn_dttrc .AND.  nn_rsttr /= 0 )                                  &
               &   CALL ctl_stop( ' ===>>>> : problem with nittrc000 for the restart',                 &
               &                  ' verify the restart file or rerun with nn_rsttr = 0 (namelist)' )
         ENDIF
         !
         IF( lk_offline ) THEN    
            !                                          ! set the date in offline mode
            IF( ln_rsttr .AND. nn_rsttr == 2 ) THEN
               CALL iom_get( numrtr, 'ndastp', zndastp ) 
               ndastp = NINT( zndastp )
               CALL iom_get( numrtr, 'adatrj', adatrj  )
             ELSE
               ndastp = ndate0 - 1     ! ndate0 read in the namelist in dom_nam
               adatrj = ( REAL( nittrc000-1, wp ) * rdttra(1) ) / rday
               ! note this is wrong if time step has changed during run
            ENDIF
            !
            IF(lwp) THEN
              WRITE(numout,*) ' *** Info used values : '
              WRITE(numout,*) '   date ndastp                                      : ', ndastp
              WRITE(numout,*) '   number of elapsed days since the begining of run : ', adatrj
              WRITE(numout,*)
            ENDIF
            !
            IF( ln_rsttr )  THEN   ;    neuler = 1
            ELSE                   ;    neuler = 0
            ENDIF
            !
            CALL day_init          ! compute calendar
            !
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         !
         IF(  kt == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'trc_wri : write the TOP restart file (NetCDF) at it= ', kt, ' date= ', ndastp
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         CALL iom_rstput( kt, nitrst, numrtw, 'kt'     , REAL( kt    , wp) )   ! time-step
         CALL iom_rstput( kt, nitrst, numrtw, 'ndastp' , REAL( ndastp, wp) )   ! date
         CALL iom_rstput( kt, nitrst, numrtw, 'adatrj' , adatrj            )   ! number of elapsed days since
         !                                                                     ! the begining of the run [s]
      ENDIF

   END SUBROUTINE trc_rst_cal


   SUBROUTINE trc_rst_stat
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_stat  ***
      !!
      !! ** purpose  :   Compute tracers statistics
      !!----------------------------------------------------------------------
      INTEGER  :: jk, jn
      REAL(wp) :: ztraf, zmin, zmax, zmean, zdrift
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zvol
      !!----------------------------------------------------------------------

      IF( lwp ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) '           ----TRACER STAT----             '
         WRITE(numout,*) 
      ENDIF
      !
      DO jk = 1, jpk
         zvol(:,:,jk) = e1e2t(:,:) * fse3t_a(:,:,jk) * tmask(:,:,jk)
      END DO
      !
      DO jn = 1, jptra
         ztraf = glob_sum( trn(:,:,:,jn) * zvol(:,:,:) )
         zmin  = MINVAL( trn(:,:,:,jn), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
         zmax  = MAXVAL( trn(:,:,:,jn), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
         IF( lk_mpp ) THEN
            CALL mpp_min( zmin )      ! min over the global domain
            CALL mpp_max( zmax )      ! max over the global domain
         END IF
         zmean  = ztraf / areatot
         zdrift = ( ( ztraf - trai(jn) ) / ( trai(jn) + 1.e-12 )  ) * 100._wp
         IF(lwp) WRITE(numout,9000) jn, TRIM( ctrcnm(jn) ), zmean, zmin, zmax, zdrift
      END DO
      IF(lwp) WRITE(numout,*) 
9000  FORMAT(' tracer nb :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, &
      &      '    max :',e18.10,'    drift :',e18.10, ' %')
      !
   END SUBROUTINE trc_rst_stat


# if defined key_medusa && defined key_roam
   SUBROUTINE trc_rst_conserve
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_conserve  ***
      !!
      !! ** purpose  :   Compute tracers conservation statistics
      !!
      !! AXY (17/11/2017)
      !! This routine calculates the "now" inventories of the elemental 
      !! cycles of MEDUSA and compares them to those calculate when the
      !! model was initialised / restarted; the cycles calculated are:
      !!    nitrogen, silicon, iron, carbon, alkalinity and oxygen
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj, jk, jn
      REAL(wp) :: zsum3d, zsum2d, zinvt, zdelta, zratio
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: z3d, zvol
      REAL(wp), DIMENSION(jpi,jpj)     :: z2d, zarea
      REAL(wp), DIMENSION(6)           :: loc_cycletot3, loc_cycletot2
      !!----------------------------------------------------------------------
      !
      IF( lwp ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) '           ----TRACER CONSERVATION----             '
         WRITE(numout,*) 
      ENDIF
      !
      ! ocean volume
      DO jk = 1, jpk
         zvol(:,:,jk) = e1e2t(:,:) * fse3t_a(:,:,jk) * tmask(:,:,jk)
      END DO
      !
      ! ocean area (for sediments)
      zarea(:,:)      = e1e2t(:,:) * tmask(:,:,1)
      !
      !----------------------------------------------------------------------
      ! nitrogen
      z3d(:,:,:) = trn(:,:,:,jpphn) + trn(:,:,:,jpphd) + trn(:,:,:,jpzmi) + &
                   trn(:,:,:,jpzme) + trn(:,:,:,jpdet) + trn(:,:,:,jpdin)
      z2d(:,:)   = zn_sed_n(:,:)
      zsum3d     = glob_sum( z3d(:,:,:) * zvol(:,:,:) )
      zsum2d     = glob_sum( z2d(:,:) * zarea(:,:) )
      ! total tracer, and delta
      zinvt      = zsum3d + zsum2d
      zdelta     = zinvt - cycletot(1)
      zratio     = 1.0e2 * zdelta / cycletot(1)
      !
      IF ( lwp ) WRITE(numout,9010) 'nitrogen', zsum3d, zsum2d, zinvt,   &
         cycletot(1), zdelta, zratio
      IF ( lwp ) WRITE(numout,*) 
      !
      !----------------------------------------------------------------------
      ! silicon
      z3d(:,:,:) = trn(:,:,:,jppds) + trn(:,:,:,jpsil)
      z2d(:,:)   = zn_sed_si(:,:)
      zsum3d     = glob_sum( z3d(:,:,:) * zvol(:,:,:) )
      zsum2d     = glob_sum( z2d(:,:) * zarea(:,:) )
      ! total tracer, and delta
      zinvt      = zsum3d + zsum2d
      zdelta     = zinvt - cycletot(2)
      zratio     = 1.0e2 * zdelta / cycletot(2)
      !
      IF ( lwp ) WRITE(numout,9010) 'silicon', zsum3d, zsum2d, zinvt,    &
         cycletot(2), zdelta, zratio
      IF ( lwp ) WRITE(numout,*) 
      !
      !----------------------------------------------------------------------
      ! iron
      z3d(:,:,:) = ((trn(:,:,:,jpphn) + trn(:,:,:,jpphd) + trn(:,:,:,jpzmi) + &
            trn(:,:,:,jpzme) + trn(:,:,:,jpdet)) * xrfn) + trn(:,:,:,jpfer)
      z2d(:,:)   = zn_sed_fe(:,:)
      zsum3d     = glob_sum( z3d(:,:,:) * zvol(:,:,:) )
      zsum2d     = glob_sum( z2d(:,:) * zarea(:,:) )
      ! total tracer, and delta
      zinvt      = zsum3d + zsum2d
      zdelta     = zinvt - cycletot(3)
      zratio     = 1.0e2 * zdelta / cycletot(3)
      !
      IF ( lwp ) WRITE(numout,9010) 'iron', zsum3d, zsum2d, zinvt,       &
         cycletot(3), zdelta, zratio
      IF ( lwp ) WRITE(numout,*) 
      !
      !----------------------------------------------------------------------
      ! carbon
      z3d(:,:,:) = (trn(:,:,:,jpphn) * xthetapn)  + (trn(:,:,:,jpphd) * xthetapd)  + &
                   (trn(:,:,:,jpzmi) * xthetazmi) + (trn(:,:,:,jpzme) * xthetazme) + &
                   trn(:,:,:,jpdtc) + trn(:,:,:,jpdic)
      z2d(:,:)   = zn_sed_c(:,:) + zn_sed_ca(:,:)
      zsum3d     = glob_sum( z3d(:,:,:) * zvol(:,:,:) )
      zsum2d     = glob_sum( z2d(:,:) * zarea(:,:) )
      ! total tracer, and delta
      zinvt      = zsum3d + zsum2d
      zdelta     = zinvt - cycletot(4)
      zratio     = 1.0e2 * zdelta / cycletot(4)
      !
      IF ( lwp ) WRITE(numout,9010) 'carbon', zsum3d, zsum2d, zinvt,     &
         cycletot(4), zdelta, zratio
      IF ( lwp ) WRITE(numout,*) 
      !
      !----------------------------------------------------------------------
      ! alkalinity
      z3d(:,:,:) = trn(:,:,:,jpalk)
      z2d(:,:)   = zn_sed_ca(:,:) * 2.0
      zsum3d     = glob_sum( z3d(:,:,:) * zvol(:,:,:) )
      zsum2d     = glob_sum( z2d(:,:) * zarea(:,:) )
      ! total tracer, and delta
      zinvt      = zsum3d + zsum2d
      zdelta     = zinvt - cycletot(5)
      zratio     = 1.0e2 * zdelta / cycletot(5)
      !
      IF ( lwp ) WRITE(numout,9010) 'alkalinity', zsum3d, zsum2d, zinvt, &
         cycletot(5), zdelta, zratio
      IF ( lwp ) WRITE(numout,*) 
      !
      !----------------------------------------------------------------------
      ! oxygen
      z3d(:,:,:) = trn(:,:,:,jpoxy)
      z2d(:,:)   = 0.0
      zsum3d     = glob_sum( z3d(:,:,:) * zvol(:,:,:) )
      zsum2d     = glob_sum( z2d(:,:) * zarea(:,:) )
      ! total tracer, and delta
      zinvt      = zsum3d + zsum2d
      zdelta     = zinvt - cycletot(6)
      zratio     = 1.0e2 * zdelta / cycletot(6)
      !
      IF ( lwp ) WRITE(numout,9010) 'oxygen', zsum3d, zsum2d, zinvt,     &
         cycletot(6), zdelta, zratio
      !
      !----------------------------------------------------------------------
      ! Check 
      zsum3d        = glob_sum( zvol(:,:,:) )
      zsum2d        = glob_sum( zarea(:,:) )
      IF ( lwp ) THEN 
         WRITE(numout,*)
         WRITE(numout,*) ' check : cvol    : ', zsum3d
         WRITE(numout,*) ' check : carea   : ', zsum2d
         WRITE(numout,*)
      ENDIF
      !
9010  FORMAT(' element:',a10,                     &
             ' 3d sum:',e18.10,' 2d sum:',e18.10, &
             ' total:',e18.10,' initial:',e18.10, &
             ' delta:',e18.10,' %:',e18.10)
      !
   END SUBROUTINE trc_rst_conserve 
# endif


#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rst_read                      ! Empty routines
   END SUBROUTINE trc_rst_read
   SUBROUTINE trc_rst_wri( kt )
      INTEGER, INTENT ( in ) :: kt
      WRITE(*,*) 'trc_rst_wri: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rst_wri   
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcrst
