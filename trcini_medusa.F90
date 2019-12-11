MODULE trcini_medusa
   !!======================================================================
   !!                         ***  MODULE trcini_medusa  ***
   !! TOP :   initialisation of the MEDUSA tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!              -   !  2008-08  (K. Popova) adaptation for MEDUSA
   !!              -   !  2008-11  (A. Yool) continuing adaptation for MEDUSA
   !!              -   !  2010-03  (A. Yool) updated for branch inclusion
   !!              -   !  2011-04  (A. Yool) updated for ROAM project
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!   'key_medusa'                                         MEDUSA tracers
   !!----------------------------------------------------------------------
   !! trc_ini_medusa   : MEDUSA model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE in_out_manager
   !! AXY (04/11/13): add this in for initialisation stuff
   USE iom
   USE par_medusa
   !! AXY (13/01/12): add this in for sediment variables
   USE sms_medusa
   !! AXY (04/11/13): add this in for initialisation stuff
   USE trcsed_medusa
   USE sbc_oce, ONLY: lk_oasis
   USE oce,     ONLY: CO2Flux_out_cpl, DMS_out_cpl, chloro_out_cpl  !! Coupling variable


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_medusa   ! called by trcini.F90 module

   !! AXY (25/02/10)
   LOGICAL, PUBLIC ::                  &
      bocalccd = .TRUE.
   !! JPALM (14/09/15)
   LOGICAL, PUBLIC ::                  &
      ln_ccd = .TRUE.

   INTEGER ::                          &
      numccd

   !! AXY (25/02/10)
   INTEGER ::                          &
      numriv

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_ini_medusa
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_medusa  ***  
      !!
      !! ** Purpose :   initialization for MEDUSA model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------

      !! vertical array index
      INTEGER  ::    jk, ierr
      !! AXY (19/07/12): added jk2 to set up friver_dep array
      INTEGER            :: jk2
      !! AXY (19/07/12): added tfthk to set up friver_dep array
      REAL(wp)           :: fthk, tfthk
      !! AXY (04/11/13): add in temporary variables for checks
      REAL(wp)           :: fq0, fq1, fq2

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: initialisation of MEDUSA model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
# if defined key_debug_medusa
            CALL flush(numout)
# endif

                                           ! Allocate MEDUSA arrays
      ierr =         sms_medusa_alloc()
# if defined key_debug_medusa
            IF (lwp) write (numout,*) '------------------------------'
            IF (lwp) write (numout,*) 'Jpalm - debug'
            IF (lwp) write (numout,*) 'in trc_ini_medusa, just after array allocate'
            IF (lwp) write (numout,*) ' '
            CALL flush(numout)
# endif

!!
!! AXY (19/07/12): setup array to control distribution of river nutrients
      friver_dep(:,:) = 0.
      DO jk = 1,jpk
         tfthk = 0.
         DO jk2 = 1,jriver_dep
            fthk  = e3t_1d(jk2)
            if (jk2 .le. jk) then
               tfthk = tfthk + fthk
               friver_dep(jk2,jk) = fthk
            endif
         ENDDO
         DO jk2 = 1,jriver_dep
            friver_dep(jk2,jk) = friver_dep(jk2,jk) / tfthk
         ENDDO
      ENDDO
!!
!! Have a look at the result of this for a single depth (jriver_dep + 1)
      IF(lwp) THEN
          WRITE(numout,*) '=== River nutrient fraction by depth (for a water column of jpk depth)'
          DO jk = 1,jpk
             WRITE(numout,*)     &
             &   ' cell = ', jk, ', friver_dep value = ', friver_dep(jk,jpk)
          ENDDO
          IF(lwp) CALL flush(numout)
       ENDIF

#if defined key_roam
!! ROAM 3D and 2D carbonate system fields (calculated on first time
!! step, then monthly)
      f3_pH(:,:,:)    = 0.
      f3_h2co3(:,:,:) = 0.
      f3_hco3(:,:,:)  = 0.
      f3_co3(:,:,:)   = 0.
      f3_omcal(:,:,:) = 0.
      f3_omarg(:,:,:) = 0.
!!
      f2_ccd_cal(:,:) = 0.
      f2_ccd_arg(:,:) = 0.
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: carbonate fields initialised to zero'
#endif
      IF(lwp) CALL flush(numout)

      !!----------------------------------------------------------------------
      !! State variable initial conditions (all mmol / m3)
      !!----------------------------------------------------------------------
      !!     
      !! biological and detrital components are initialised to nominal
      !! values above 100 m depth and zero below; the latter condition
      !! is applied since non-linear loss processes allow significant
      !! concentrations of these components to persist at depth
      !!
      trn(:,:,:,jpchn) = 0.
      trn(:,:,:,jpchd) = 0.
      trn(:,:,:,jpphn) = 0.
      trn(:,:,:,jpphd) = 0.
      trn(:,:,:,jppds) = 0.
      trn(:,:,:,jpzmi) = 0.
      trn(:,:,:,jpzme) = 0.
      trn(:,:,:,jpdet) = 0.
      !!
      DO jk = 1,13
         !! non-diatom chlorophyll         (nominal)
         trn(:,:,jk,jpchn) = 0.01
         !!
         !! diatom chlorophyll             (nominal)
         trn(:,:,jk,jpchd) = 0.01
         !!
         !! non-diatom                     (nominal)
         trn(:,:,jk,jpphn) = 0.01
         !!
         !! diatom                         (nominal)
         trn(:,:,jk,jpphd) = 0.01
         !!
         !! diatom silicon                 (nominal)
         trn(:,:,jk,jppds) = 0.01
         !!
         !! microzooplankton               (nominal)
         trn(:,:,jk,jpzmi) = 0.01
         !!
         !! mesozooplankton                (nominal)
         trn(:,:,jk,jpzme) = 0.01
         !!
         !! detrital nitrogen              (nominal)
         trn(:,:,jk,jpdet) = 0.01
      ENDDO
      !!
      !! dissolved inorganic nitrogen     (nominal average value; typically initialised from climatology)
      trn(:,:,:,jpdin) = 30.
      !!
      !! dissolved silicic acid           (nominal average value; typically initialised from climatology)
      trn(:,:,:,jpsil) = 90.
      !!
      !! dissolved "total" iron           (nominal; typically initialised from model-derived climatology)
      trn(:,:,:,jpfer) = 1.0e-4           !! = 0.1 umol Fe / m3
      !!
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: MEDUSA-1 fields initialised to defaults'
# if defined key_roam
      !!
      !! detrital carbon                  (nominal)
      trn(:,:,:,jpdtc) = 0.
      DO jk = 1,13
         trn(:,:,jk,jpdtc) = 0.06625
      ENDDO
      !!
      !! dissolved inorganic carbon (DIC) (nominal average value; typically initialised from climatology)
      trn(:,:,:,jpdic) = 2330.
      !!
      !! total alkalinity                 (nominal average value; typically initialised from climatology)
      trn(:,:,:,jpalk) = 2450.
      !!
      !! dissolved oxygen                 (nominal average value; typically initialised from climatology)
      trn(:,:,:,jpoxy) = 175.
      !!
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: MEDUSA-2 fields initialised to defaults'
# endif
      IF(lwp) CALL flush(numout)

      !!----------------------------------------------------------------------
      !! Sediment pools initial conditions (all mmol / m2)
      !!----------------------------------------------------------------------
      !!     
      !! these pools store biogenic material that has sunk to the seabed,
      !! and act as a temporary reservoir
      zb_sed_n(:,:)  = 0.0  !! organic N
      zn_sed_n(:,:)  = 0.0
      za_sed_n(:,:)  = 0.0
      zb_sed_fe(:,:) = 0.0  !! organic Fe
      zn_sed_fe(:,:) = 0.0
      za_sed_fe(:,:) = 0.0
      zb_sed_si(:,:) = 0.0  !! inorganic Si
      zn_sed_si(:,:) = 0.0
      za_sed_si(:,:) = 0.0
      zb_sed_c(:,:)  = 0.0  !! organic C
      zn_sed_c(:,:)  = 0.0
      za_sed_c(:,:)  = 0.0
      zb_sed_ca(:,:) = 0.0  !! inorganic C
      zn_sed_ca(:,:) = 0.0
      za_sed_ca(:,:) = 0.0
      !!
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: benthic fields initialised to zero'
      IF(lwp) CALL flush(numout)
     
      !!----------------------------------------------------------------------
      !! Averaged properties for DMS calculations (various units)
      !!----------------------------------------------------------------------
      !!     
      !! these store temporally averaged properties for DMS calculations (AXY, 07/07/15)
      zb_dms_chn(:,:)  = 0.0  !! CHN
      zn_dms_chn(:,:)  = 0.0
      za_dms_chn(:,:)  = 0.0
      zb_dms_chd(:,:)  = 0.0  !! CHD
      zn_dms_chd(:,:)  = 0.0
      za_dms_chd(:,:)  = 0.0
      zb_dms_mld(:,:)  = 0.0  !! MLD
      zn_dms_mld(:,:)  = 0.0
      za_dms_mld(:,:)  = 0.0
      zb_dms_qsr(:,:)  = 0.0  !! QSR
      zn_dms_qsr(:,:)  = 0.0
      za_dms_qsr(:,:)  = 0.0
      zb_dms_din(:,:)  = 0.0  !! DIN
      zn_dms_din(:,:)  = 0.0
      za_dms_din(:,:)  = 0.0
      !!
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: average fields for DMS initialised to zero'
      IF(lwp) CALL flush(numout)
      !!
      !!---------------------------------------------------------------------
      !!JPALM (14-06-2016): init dms and co2 flux for coupling with atm (UKESM)
      !!---------------------------------------------------------------------
      !!
      zb_co2_flx(:,:)  = 0.0  !! CO2 flx
      zn_co2_flx(:,:)  = 0.0
      za_co2_flx(:,:)  = 0.0
      zb_dms_srf(:,:)  = 0.0  !! DMS srf
      zn_dms_srf(:,:)  = 0.0
      za_dms_srf(:,:)  = 0.0
      zn_chl_srf(:,:)  = 2.0E-8 !! Chl cpl - set first as surf
      !!
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: DMS and CO2 flux (UKESM) initialised to zero'
      IF(lwp) CALL flush(numout)
      IF (lk_oasis) THEN
         CO2Flux_out_cpl(:,:) =  zn_co2_flx(:,:)   !! Coupling variable
         DMS_out_cpl(:,:)     =  zn_dms_srf(:,:)   !! Coupling variable
         chloro_out_cpl(:,:)  =  zn_chl_srf(:,:) * scl_chl   !! Coupling variable
      END IF
      !!
      !!----------------------------------------------------------------------
      !! AXY (04/11/13): initialise fields previously done by trc_sed_medusa
      !!----------------------------------------------------------------------
      !!     
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: initialising dust deposition fields'
      CALL trc_sed_medusa_sbc( nit000 )
      !!
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: initialising ocean CCD array'
      CALL trc_ini_medusa_ccd( nit000 )
      fq0 = MINVAL(ocal_ccd(:,:))
      fq1 = MAXVAL(ocal_ccd(:,:))
      if (lwp) write (numout,'(a,f10.3,a,f10.3)') & 
         & 'CCD: min ', fq0, ' max ', fq1
      !!
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: initialising riverine nutrient arrays'
      riv_n(:,:)   = 0.0
      riv_si(:,:)  = 0.0
      riv_c(:,:)   = 0.0
      riv_alk(:,:) = 0.0 
      !!
      CALL trc_ini_medusa_river( nit000 )
      fq0 = MINVAL(riv_n(:,:))
      fq1 = MAXVAL(riv_n(:,:))
      if (lwp) write (numout,'(a,f10.3,a,f10.3)') & 
         & 'RIV_N:   min ', fq0, ' max ', fq1
      fq0 = MINVAL(riv_si(:,:))
      fq1 = MAXVAL(riv_si(:,:))
      if (lwp) write (numout,'(a,f10.3,a,f10.3)') & 
         & 'RIV_SI:  min ', fq0, ' max ', fq1
      fq0 = MINVAL(riv_c(:,:))
      fq1 = MAXVAL(riv_c(:,:))
      if (lwp) write (numout,'(a,f10.3,a,f10.3)') & 
         & 'RIV_C:   min ', fq0, ' max ', fq1
      fq0 = MINVAL(riv_alk(:,:))
      fq1 = MAXVAL(riv_alk(:,:))
      if (lwp) write (numout,'(a,f10.3,a,f10.3)') & 
         & 'RIV_ALK: min ', fq0, ' max ', fq1
      IF(lwp) CALL flush(numout)

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa: MEDUSA initialised'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      IF(lwp) CALL flush(numout)
      !!
      !!----------------------------------------------------------------------
      !! JPALM (23-01-2017): new way to initialize CO2-atm for cmip6 
      !!                     initially done in trcsms_medusa
      !!----------------------------------------------------------------------
      !! 
      IF( ( .NOT.lk_oasis ) .AND. ( .NOT.lk_pi_co2 ) ) THEN
         IF(lwp) WRITE(numout,*) ' trc_ini_medusa: initialisating atm CO2 record'
         CALL trc_ini_medusa_co2atm
      ENDIF

   END SUBROUTINE trc_ini_medusa

   !! ======================================================================
   !! ======================================================================
   !! ======================================================================

   !! AXY (25/02/10)
   SUBROUTINE trc_ini_medusa_ccd(kt)

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_ini_medusa_ccd  ***
      !!
      !! ** Purpose :   Read CCD field
      !!
      !! ** Method  :   Read the file
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !!---------------------------------------------------------------------

      !! Open the file
      !! -------------
      !!
      !!!! JPALM -- 14-09-2015 -- 
      !!!!       -- to test on ORCA2 with Christian, no file available, so initiate to 0 
      IF (ln_ccd) THEN
         IF(lwp) WRITE(numout,*) ' '
         IF(lwp) WRITE(numout,*) ' **** Routine trc_ini_medusa_ccd'
         CALL iom_open ( 'ccd_ocal_nemo.nc', numccd )
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_ccd: ccd_ocal_nemo.nc opened'

      !! Read the data
      !! -------------
      !!
         CALL iom_get ( numccd, jpdom_data, 'OCAL_CCD', ocal_ccd )
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_ccd: data read'

      !! Close the file
      !! --------------
      !!
         CALL iom_close ( numccd )
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_ccd: ccd_ocal_nemo.nc closed'
         IF(lwp) CALL flush(numout)
      ELSE
         IF(lwp) WRITE(numout,*) ' '
         IF(lwp) WRITE(numout,*) ' **** Routine trc_ini_medusa_ccd'
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_ccd: do not read ccd_ocal_nemo.nc'
         IF(lwp) WRITE(numout,*) ' **** ln_ccd = FALSE and ocal_ccd = 0.0 ---'
         ocal_ccd(:,:) = 0.0  
      ENDIF
 
   END SUBROUTINE trc_ini_medusa_ccd

   !! ======================================================================
   !! ======================================================================
   !! ======================================================================

   !! AXY (26/01/12)
   SUBROUTINE trc_ini_medusa_river(kt)

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_ini_medusa_river  ***
      !!
      !! ** Purpose :   Read riverine nutrient fields
      !!
      !! ** Method  :   Read the file
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !!---------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' **** Routine trc_ini_medusa_river'
         WRITE(numout,*) ' '
      ENDIF

      !! Open and read the files
      !! -----------------------
      !!
      if (jriver_n.gt.0) then
         if (jriver_n.eq.1) CALL iom_open ( 'river_N_conc_orca100.nc', numriv )
         if (jriver_n.eq.2) CALL iom_open ( 'river_N_flux_orca100.nc', numriv )
         CALL iom_get  ( numriv, jpdom_data, 'RIV_N', riv_n )
         IF(lwp) THEN
            if (jriver_n.eq.1) WRITE(numout,*) ' **** trc_ini_medusa_river: N CONC data read'
            if (jriver_n.eq.2) WRITE(numout,*) ' **** trc_ini_medusa_river: N FLUX data read'
         ENDIF
         CALL iom_close ( numriv )
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_river: river N file closed'
      else
         IF(lwp) THEN
            WRITE(numout,*) ' **** trc_ini_medusa_river: N data NOT read'
         ENDIF
      endif
      !!
      if (jriver_si.gt.0) then
         if (jriver_si.eq.1) CALL iom_open ( 'river_Si_conc_orca100.nc', numriv )
         if (jriver_si.eq.2) CALL iom_open ( 'river_Si_flux_orca100.nc', numriv )
         CALL iom_get  ( numriv, jpdom_data, 'RIV_SI', riv_si )
         IF(lwp) THEN
            if (jriver_si.eq.1) WRITE(numout,*) ' **** trc_ini_medusa_river: Si CONC data read'
            if (jriver_si.eq.2) WRITE(numout,*) ' **** trc_ini_medusa_river: Si FLUX data read'
         ENDIF
         CALL iom_close ( numriv )
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_river: river Si file closed'
      else
         IF(lwp) THEN
            WRITE(numout,*) ' **** trc_ini_medusa_river: Si data NOT read'
         ENDIF
      endif
      !!
      if (jriver_c.gt.0) then
         if (jriver_c.eq.1) CALL iom_open ( 'river_C_conc_orca100.nc', numriv )
         if (jriver_c.eq.2) CALL iom_open ( 'river_C_flux_orca100.nc', numriv )
         CALL iom_get  ( numriv, jpdom_data, 'RIV_C', riv_c )
         IF(lwp) THEN
            if (jriver_c.eq.1) WRITE(numout,*) ' **** trc_ini_medusa_river: C CONC data read'
            if (jriver_c.eq.2) WRITE(numout,*) ' **** trc_ini_medusa_river: C FLUX data read'
         ENDIF
         CALL iom_close ( numriv )
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_river: river C file closed'
      else
         IF(lwp) THEN
            WRITE(numout,*) ' **** trc_ini_medusa_river: C data NOT read'
         ENDIF
      endif
      !!
      if (jriver_alk.gt.0) then
         if (jriver_alk.eq.1) CALL iom_open ( 'river_alk_conc_orca100.nc', numriv )
         if (jriver_alk.eq.2) CALL iom_open ( 'river_alk_flux_orca100.nc', numriv )
         CALL iom_get  ( numriv, jpdom_data, 'RIV_ALK', riv_alk )
         IF(lwp) THEN
            if (jriver_alk.eq.1) WRITE(numout,*) ' **** trc_ini_medusa_river: alkalinity CONC data read'
            if (jriver_alk.eq.2) WRITE(numout,*) ' **** trc_ini_medusa_river: alkalinity FLUX data read'
         ENDIF
         CALL iom_close ( numriv )
         IF(lwp) WRITE(numout,*) ' **** trc_ini_medusa_river: river alkalinity file closed'
      else
         IF(lwp) THEN
            WRITE(numout,*) ' **** trc_ini_medusa_river: alkalinity data NOT read'
         ENDIF
      endif
      IF(lwp) CALL flush(numout)

   END SUBROUTINE trc_ini_medusa_river
   
   SUBROUTINE trc_ini_medusa_co2atm
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_medusa_co2atm  ***  
      !!
      !! ** Purpose :   initialization atmospheric co2 record
      !!
      !! ** Method  : - Read the xco2 file
      !!----------------------------------------------------------------------
      INTEGER                       ::  jn, jm, io, ierr, inum, iostatus
      INTEGER, PARAMETER            ::  iskip = 4   ! number of 1st descriptor lines
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)   ::   zyy !: xCO2 record years
      CHARACTER (len=10), PARAMETER ::  clname = 'xco2.atm'  !! atm CO2 record file
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_medusa_co2atm: initialisation of atm CO2 historical record'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~'


      IF(lwp) WRITE(numout,*) 'read of formatted file xco2.atm'

      CALL ctl_opn( inum, clname, 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !!!
      ! -Compute the number of year in the file
      ! -File starts in co2_yinit, jn represents the record number in the file.
      ! -Remove the file head (iskip lines) to jn
      ! -The year is jn + yinit - 1 
      !! Determine the number of lines in xCO2 input file
      iostatus = 0
      jn = 1
      DO WHILE ( iostatus == 0 )
        READ(inum,'(1x)', IOSTAT=iostatus, END=100)
        jn = jn + 1
      ENDDO
      IF( iostatus .NE. 0 ) THEN
        !! Error while reading xCO2 input file 
        CALL ctl_stop('trc_ini_medusa_co2atm: &
                      & Error on the 1st reading of xco2.atm')
        RETURN
      ENDIF
 100  co2_rec = jn - 1 - iskip
      IF ( lwp) WRITE(numout,*) '    ', co2_rec ,' years read in the file'
      !                                ! Allocate CO2 hist arrays
      ierr = 0 
      ALLOCATE( hist_pco2(co2_rec),zyy(co2_rec), STAT=ierr )
      IF( ierr > 0 ) THEN
         CALL ctl_stop( 'trc_ini_medusa_co2atm: unable to allocate  array' )  
         RETURN
      ENDIF

      REWIND(inum)

      DO jm = 1, iskip        ! Skip over 1st six descriptor lines
         READ(inum,'(1x)')
      END DO
      ! file starts in 1931 do jn represent the year in the century.jhh
      ! Read file till the end
      ! allocate start and end year of the file
      DO jn = 1, co2_rec
        READ(inum,'(F6.1,F12.7)', IOSTAT=io) zyy(jn), hist_pco2(jn)
        IF( io .NE. 0 ) THEN
          !! Error while reading xCO2 input file 
          CALL ctl_stop('trc_ini_medusa_co2atm: &
                        & Error on the 2nd reading of xco2.atm')
          RETURN
        ENDIF

        IF(jn==1) co2_yinit = zyy(jn)
      END DO
      co2_yend = co2_yinit + real(co2_rec) - 1.

      IF(lwp) THEN        ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'CO2 hist start year: ', co2_yinit
         WRITE(numout,*) 'CO2 hist end   year: ', co2_yend
         WRITE(numout,*) ' Year   xCO2 atm '
         DO jn = 1, co2_rec
            WRITE(numout, '(F6.1,F12.7)') zyy(jn), hist_pco2(jn)
         END DO
      ENDIF

   END SUBROUTINE trc_ini_medusa_co2atm


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MEDUSA model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_medusa             ! Empty routine
   END SUBROUTINE trc_ini_medusa
#endif

   !!======================================================================
END MODULE trcini_medusa
