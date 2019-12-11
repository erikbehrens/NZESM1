MODULE trcbio_medusa
   !!======================================================================
   !!                         ***  MODULE trcbio  ***
   !! TOP :   MEDUSA
   !!======================================================================
   !! History :
   !!  -   !  1999-07  (M. Levy)              original code
   !!  -   !  2000-12  (E. Kestenare)         assign parameters to name 
   !!                                         individual tracers
   !!  -   !  2001-03  (M. Levy)              LNO3 + dia2d 
   !! 2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!  -   !  2008-08  (K. Popova)            adaptation for MEDUSA
   !!  -   !  2008-11  (A. Yool)              continuing adaptation for MEDUSA
   !!  -   !  2010-03  (A. Yool)              updated for branch inclusion
   !!  -   !  2011-08  (A. Yool)              updated for ROAM (see below)
   !!  -   !  2013-03  (A. Yool)              updated for iMARNET
   !!  -   !  2013-05  (A. Yool)              updated for v3.5
   !!  -   !  2014-08  (A. Yool, J. Palm)     Add DMS module for UKESM1 model
   !!  -   !  2015-06  (A. Yool)              Update to include MOCSY
   !!  -   !  2015-07  (A. Yool)              Update for rolling averages
   !!  -   !  2015-10  (J. Palm)              Update for diag outputs through 
   !!                                         iom_use  
   !!  -   !  2016-11  (A. Yool)              Updated diags for CMIP6
   !!  -   !  2017-05  (A. Yool)              Added extra DMS calculation
   !!  -   !  2017-11  (J. Palm, A. Yool)     Diagnose tracer excursions
   !!----------------------------------------------------------------------
   !!
#if defined key_roam
   !!----------------------------------------------------------------------
   !! Updates for the ROAM project include:
   !!   - addition of DIC, alkalinity, detrital carbon and oxygen tracers
   !!   - addition of air-sea fluxes of CO2 and oxygen
   !!   - periodic (monthly) calculation of full 3D carbonate chemistry 
   !!   - detrital C:N ratio now free to evolve dynamically
   !!   - benthic storage pools
   !! 
   !! Opportunity also taken to add functionality:
   !!   - switch for Liebig Law (= most-limiting) nutrient uptake
   !!   - switch for accelerated seafloor detritus remineralisation
   !!   - switch for fast -> slow detritus transfer at seafloor
   !!   - switch for ballast vs. Martin vs. Henson fast detritus remin.
   !!   - per GMD referee remarks, xfdfrac3 introduced for grazed PDS
   !!----------------------------------------------------------------------
#endif
   !!
#if defined key_mocsy
   !!----------------------------------------------------------------------
   !! Updates with the addition of MOCSY include:
   !!   - option to use PML or MOCSY carbonate chemistry (the latter is 
   !!     preferred)
   !!   - central calculation of gas transfer velocity, f_kw660; previously
   !!     this was done separately for CO2 and O2 with predictable results
   !!   - distribution of f_kw660 to both PML and MOCSY CO2 air-sea flux 
   !!     calculations and to those for O2 air-sea flux
   !!   - extra diagnostics included for MOCSY
   !!----------------------------------------------------------------------
#endif
   !!
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                        MEDUSA bio-model
   !!----------------------------------------------------------------------
   !!   trc_bio_medusa        :  
   !!----------------------------------------------------------------------
      !! AXY (30/01/14): necessary to find NaNs on HECTOR
!CEB
!$AGRIF_DO_NOT_TREAT
      USE, INTRINSIC :: ieee_arithmetic 
!$AGRIF_END_DO_NOT_TREAT
!/CEB
      USE bio_medusa_mod,             ONLY: b0, fdep1,                      & 
                                            ibenthic, idf, idfval,          &
# if defined key_roam
                                            f_xco2a,                        &
                                            zalk, zdic, zoxy, zsal, ztmp,   &
# endif
# if defined key_mocsy
                                            zpho,                           &
# endif
                                            zchd, zchn, zdet, zdin, zdtc,   &
                                            zfer, zpds, zphd, zphn, zsil,   &
                                            zzme, zzmi
      USE dom_oce,                    ONLY: e3t_0, gdept_0,  gdepw_0,       &
                                            nday_year, nsec_day,            & 
                                            nyear, nyear_len, ndastp,       &
                                            nsec_month,                     &
                                            rdt, tmask, mig, mjg, nimpp,    &
                                            njmpp 
#if defined key_vvl
      USE dom_oce,                    ONLY: e3t_n, gdept_n,  gdepw_n
#endif

      USE in_out_manager,             ONLY: lwp, numout, nn_date0
      USE iom,                        ONLY: lk_iomput
      USE lbclnk,                     ONLY: lbc_lnk
      USE lib_mpp,                    ONLY: mpp_max, mpp_maxloc,            & 
                                            mpp_min, mpp_minloc,            &
                                            ctl_stop, ctl_warn, lk_mpp  
      USE oce,                        ONLY: tsb, tsn, PCO2a_in_cpl
      USE par_kind,                   ONLY: wp
      USE par_medusa,                 ONLY: jpalk, jpchd, jpchn, jpdet,     &
                                            jpdic, jpdin, jpdtc, jpfer,     &
                                            jpoxy, jppds, jpphd, jpphn,     &
                                            jpsil, jpzme, jpzmi
      USE par_oce,                    ONLY: jp_sal, jp_tem, jpi, jpim1,     &
                                            jpj, jpjm1, jpk
      !! JPALM (27-06-2016): add lk_oasis for CO2 and DMS coupling with atm
      USE sbc_oce,                    ONLY: lk_oasis
      USE sms_medusa,                 ONLY: hist_pco2, co2_yinit, co2_yend, &
                                            lk_pi_co2
      USE trc,                        ONLY: ln_rsttr, nittrc000, trn
      USE bio_medusa_init_mod,        ONLY: bio_medusa_init
      USE carb_chem_mod,              ONLY: carb_chem
      USE air_sea_mod,                ONLY: air_sea
      USE plankton_mod,               ONLY: plankton
      USE iron_chem_scav_mod,         ONLY: iron_chem_scav
      USE detritus_mod,               ONLY: detritus
      USE bio_medusa_update_mod,      ONLY: bio_medusa_update
      USE bio_medusa_diag_mod,        ONLY: bio_medusa_diag
      USE bio_medusa_diag_slice_mod,  ONLY: bio_medusa_diag_slice
      USE bio_medusa_fin_mod,         ONLY: bio_medusa_fin
      USE trcstat,                    ONLY: trc_rst_dia_stat

      IMPLICIT NONE
      PRIVATE
      
      PUBLIC   trc_bio_medusa    ! called in trcsms_medusa.F90
      PUBLIC   trc_bio_exceptionnal_fix ! here 

   !!* Substitution
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_bio_medusa( kt )
      !!------------------------------------------------------------------
      !!                     ***  ROUTINE trc_bio  ***
      !!
      !! ** Purpose : compute the now trend due to biogeochemical processes
      !!              and add it to the general trend of passive tracers 
      !!              equations
      !!
      !! ** Method  : each now biological flux is calculated in function of
      !!              now concentrations of tracers.
      !!              depending on the tracer, these fluxes are sources or 
      !!              sinks.
      !!              The total of the sources and sinks for each tracer
      !!              is added to the general trend.
      !!        
      !!                      tra = tra + zf...tra - zftra...
      !!                                     |         |
      !!                                     |         |
      !!                                  source      sink
      !!        
      !!              IF 'key_trc_diabio' defined , the biogeochemical trends
      !!              for passive tracers are saved for futher diagnostics.
      !!------------------------------------------------------------------
      !!
      !!
      !!------------------------------------------------------------------
      !! Variable conventions
      !!------------------------------------------------------------------
      !!
      !! names: z*** - state variable 
      !!        f*** - function (or temporary variable used in part of 
      !!               a function)
      !!        x*** - parameter
      !!        b*** - right-hand part (sources and sinks)
      !!        i*** - integer variable (usually used in yes/no flags)
      !!
      !! time (integer timestep)
      INTEGER, INTENT( in ) ::    kt
      !!
      !! spatial array indices
      INTEGER  ::    ji,jj,jk,jn
      !!
      INTEGER  ::    iball
# if defined key_roam
      !!
      INTEGER  ::    iyr1, iyr2
      !!
# endif
      !! 
      !! temporary variables
      REAL(wp) ::    fq3,fq4
      REAL(wp) ::    this_y, this_d, this_s, fyear
      !!
      !! T and S check temporary variable
      REAL(wp) ::    sumtsn, tsnavg
      INTEGER  ::    summask
      CHARACTER(40) :: charout, charout2, charout3, charout4, charout5
      !!
      !!------------------------------------------------------------------

# if defined key_debug_medusa
      IF (lwp) write (numout,*) 'trc_bio_medusa: variables defined'
      CALL flush(numout)
# endif 

      !! AXY (20/11/14): alter this to report on first MEDUSA call
      !! IF( kt == nit000 ) THEN
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_bio: MEDUSA bio-model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
	 IF(lwp) WRITE(numout,*) ' kt =',kt
      ENDIF

      !! AXY (13/01/12): is benthic model properly interactive? 0 = no, 1 = yes
      ibenthic = 1

      !!------------------------------------------------------------------
      !! b0 is present for debugging purposes; using b0 = 0 sets the tendency
      !! terms of all biological equations to 0.
      !!------------------------------------------------------------------
      !!
      !! AXY (03/09/14): probably not the smartest move ever, but it'll fit
      !!                 the bill for now; another item on the things-to-sort-
      !!		 out-in-the-future list ...
# if defined key_kill_medusa
      b0 = 0.
# else
      b0 = 1.
# endif
      !!------------------------------------------------------------------
      !! fast detritus ballast scheme (0 = no; 1 = yes)
      !! alternative to ballast scheme is same scheme but with no ballast
      !! protection (not dissimilar to Martin et al., 1987)
      !!------------------------------------------------------------------
      !!
      iball = 1

      !!------------------------------------------------------------------
      !! full flux diagnostics (0 = no; 1 = yes); appear in ocean.output
      !! these should *only* be used in 1D since they give comprehensive
      !! output for ecological functions in the model; primarily used in
      !! debugging
      !!------------------------------------------------------------------
      !!
      idf    = 0
      !!
      !! timer mechanism
      if (kt/120*120.eq.kt) then
         idfval = 1
      else
         idfval = 0
      endif

      !!------------------------------------------------------------------
      !! Initialise arrays to zero and set up arrays for diagnostics
      !!------------------------------------------------------------------
      CALL bio_medusa_init( kt )

# if defined key_axy_nancheck
       DO jn = jp_msa0,jp_msa1
          !! fq0 = MINVAL(trn(:,:,:,jn))
          !! fq1 = MAXVAL(trn(:,:,:,jn))
          fq2 = SUM(trn(:,:,:,jn))
          !! if (lwp) write (numout,'(a,2i6,3(1x,1pe15.5))') 'NAN-CHECK',     &
          !!                kt, jn, fq0, fq1, fq2
          !! AXY (30/01/14): much to our surprise, the next line doesn't 
          !!                 work on HECTOR and has been replaced here with 
          !!                 a specialist routine
          !! if (fq2 /= fq2 ) then
          if ( ieee_is_nan( fq2 ) ) then
             !! there's a NaN here
             if (lwp) write(numout,*) 'NAN detected in field', jn,           &
                                      'at time', kt, 'at position:'
             DO jk = 1,jpk
                DO jj = 1,jpj
                   DO ji = 1,jpi
                      !! AXY (30/01/14): "isnan" problem on HECTOR
                      !! if (trn(ji,jj,jk,jn) /= trn(ji,jj,jk,jn)) then
                      if ( ieee_is_nan( trn(ji,jj,jk,jn) ) ) then
                         if (lwp) write (numout,'(a,1pe12.2,4i6)')           &
                            'NAN-CHECK', tmask(ji,jj,jk), ji, jj, jk, jn
                      endif
                   enddo
                enddo
             enddo
             CALL ctl_stop( 'trcbio_medusa, NAN in incoming tracer field' )
          endif
       ENDDO
       CALL flush(numout)
# endif

# if defined key_debug_medusa
      IF (lwp) write (numout,*)                                              &
                     'trc_bio_medusa: variables initialised and checked'
      CALL flush(numout)
# endif 

# if defined key_roam
      !!------------------------------------------------------------------
      !! calculate atmospheric pCO2
      !!------------------------------------------------------------------
      !!
      IF (lk_oasis) THEN 
         !! xCO2 from coupled
         IF ( ( kt == nittrc000 ) .AND. lwp )     &
             WRITE(numout,*) '** MEDUSA Atm xCO2 given by the UM **'
         f_xco2a(:,:) = PCO2a_in_cpl(:,:)
         !! Check the xCO2 from the UM is OK
         !! piece of code moved from air-sea.F90
         !!---
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               !! OPEN wet point IF..THEN loop
               IF (tmask(ji,jj,1) == 1) then
                  !!! Jpalm test on atm xCO2
                  IF ( (f_xco2a(ji,jj) .GT. 10000.0 ).OR.     &
                       (f_xco2a(ji,jj) .LE. 0.0 ) ) THEN
                     IF(lwp) THEN
                        WRITE(numout,*) ' atm xCO2 = ',f_xco2a(ji,jj),       &
                                        ' -- ji =', mig(ji),' jj = ', mjg(jj)
                     ENDIF
                     CALL ctl_stop( 'MEDUSA - trc_bio :',     &
                                    'unrealistic coupled atm xCO2 ' )
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         !!--- 
      ELSEIF (lk_pi_co2) THEN
         !! OCMIP pre-industrial xCO2
         IF ( ( kt == nittrc000 ) .AND. lwp )     &
             WRITE(numout,*) '** MEDUSA Atm xCO2 fixed to pre-industrial value **'
         !! f_xco2a(:,:) = 284.725       !! CMIP5 pre-industrial pCO2
         f_xco2a(:,:) = 284.317          !! CMIP6 pre-industrial pCO2 
      ELSE      
         !! xCO2 from file
         !! AXY - JPALM new interpolation scheme usinf nyear_len
         this_y = real(nyear)
         this_d = real(nday_year)
         this_s = real(nsec_day)
         !!
         fyear = this_y + ((this_d - 1) + (this_s / (60. * 60. * 24.))) / real(nyear_len(1))
         !!
         IF ( ( kt == nittrc000 ) .AND. lwp ) THEN 
            WRITE(numout,*) '** MEDUSA Atm xCO2 from file **'
            WRITE(numout,*) ' MEDUSA year      =', this_y
            WRITE(numout,*) ' Year length      =', real(nyear_len(1))
            WRITE(numout,*) ' MEDUSA nday_year =', this_d
            WRITE(numout,*) ' MEDUSA nsec_day  =', this_s
         ENDIF
         !! 
         !! different case test 
         IF (fyear .LE. co2_yinit) THEN
            !! before first record -- pre-industrial value
            f_xco2a(:,:) = hist_pco2(1)  
         ELSEIF (fyear .GE. co2_yend) THEN
            !! after last record - continue to use the last value
            f_xco2a(:,:) = hist_pco2(int(co2_yend - co2_yinit + 1.) )
         ELSE
            !! just right
            iyr1 = int(fyear - co2_yinit) + 1
            iyr2 = iyr1 + 1 
            fq3 = fyear - real(iyr1) - co2_yinit + 1.
            fq4 = ((1 - fq3) * hist_pco2(iyr1)) + (fq3 * hist_pco2(iyr2))
            f_xco2a(:,:) = fq4
            !!
            IF ( ( kt == nittrc000 ) .AND. lwp ) THEN 
               WRITE(numout,*) ' MEDUSA year1     =', iyr1
               WRITE(numout,*) ' MEDUSA year2     =', iyr2
               WRITE(numout,*) ' xCO2 year1       =', hist_pco2(iyr1)
               WRITE(numout,*) ' xCO2 year2       =', hist_pco2(iyr2)
               WRITE(numout,*) ' Year2 weight     =', fq3
            ENDIF
         ENDIF
      ENDIF
 
      !! Writing xCO2 in output on start and on the 1st tsp of each  month
      IF ( (kt == nittrc000 .AND. .NOT.ln_rsttr) .OR.                        &
           ( nsec_month .LE. INT(rdt) ) )  THEN
           IF ( lwp )  WRITE(numout,*) ' *** Atm xCO2 *** -- kt:', kt,       &
                                       ';  current date:', ndastp
          call trc_rst_dia_stat(f_xco2a(:,:), 'atm xCO2')
      ENDIF
# endif

# if defined key_debug_medusa
      IF (lwp) write (numout,*) 'trc_bio_medusa: ready for carbonate chemistry'
      IF (lwp) write (numout,*) 'trc_bio_medusa: kt = ', kt
      IF (lwp) write (numout,*) 'trc_bio_medusa: nittrc000 = ', nittrc000
      CALL flush(numout)
# endif 

# if defined key_roam
      !! AXY (20/11/14): alter to call on first MEDUSA timestep and then every
      !!                 month (this is hardwired as 960 timesteps but should
      !!                 be calculated and done properly
      !! IF( kt == nit000 .or. mod(kt,1920) == 0 ) THEN
      !! IF( kt == nittrc000 .or. mod(kt,960) == 0 ) THEN 
      !!=============================
      !! Jpalm -- 07-10-2016 -- need to change carb-chem frequency call :
      !!          we don't want to call on the first time-step of all run 
      !!          submission, but only on the very first time-step, and 
      !!          then every month. So we call on nittrc000 if not 
      !!          restarted run, else if one month after last call.
      !!          Assume one month is 30d --> 3600*24*30 : 2592000s
      !!          try to call carb-chem at 1st month's tm-stp : 
      !!          x * 30d + 1*rdt(i.e: mod = rdt)   
      !!          ++ need to pass carb-chem output var through restarts
      !!If ( (kt == nitt8rc000 .AND. .NOT.ln_rsttr) .OR.          &
      !!     ( (mod(kt*rdt,2592000.)) == rdt) THEN
      !!=============================
      !! (Jpalm -- updated for restartability issues)
      !! We want this to be start of month or if starting afresh from  
      !! climatology - marc 20/6/17 
      !!If ( (kt == nittrc000 .AND. .NOT.ln_rsttr) .OR.                         & 
      !!     ((86400*mod(nn_date0,100) + mod(kt*rdt,2592000.)) == rdt) ) THEN 
      !!=============================
      !! Jpalm -- 15-02-2018 -- need to change 3D carb-chem call freq again.
      !! previous call did not work, probably the (86400*mod(nn_date0,100) part
      !! should not be in...
      !! now use the NEMO calendar tool : nsec_month to be sure to call 
      !! at the beginning of a new month .
      IF ( (kt == nittrc000 .AND. .NOT.ln_rsttr) .OR.                        &
           ( nsec_month .LE. INT(rdt) ) )  THEN
           IF ( lwp )  WRITE(numout,*)                                       &         
                              ' *** 3D carb chem call *** -- kt:', kt,       &
                              ';  current date:', ndastp 
         !!---------------------------------------------------------------
         !! Calculate the carbonate chemistry for the whole ocean on the first
         !! simulation timestep and every month subsequently; the resulting 3D
         !! field of omega calcite is used to determine the depth of the CCD
         !!---------------------------------------------------------------
         CALL carb_chem( kt )

      ENDIF
# endif

# if defined key_debug_medusa
      IF (lwp) write (numout,*) 'trc_bio_medusa: ready for full domain calculations'
      CALL flush(numout)
# endif 

      !!------------------------------------------------------------------
      !! MEDUSA has unified equation through the water column
      !! (Diff. from LOBSTER which has two sets: bio- and non-bio layers) 
      !! Statement below in LOBSTER is different: DO jk = 1, jpkbm1          
      !!------------------------------------------------------------------
      !!
      !! NOTE: the ordering of the loops below differs from that of some other
      !! models; looping over the vertical dimension is the outermost loop and
      !! this complicates some calculations (e.g. storage of vertical fluxes
      !! that can otherwise be done via a singular variable require 2D fields
      !! here); however, these issues are relatively easily resolved, but the
      !! loops CANNOT be reordered without potentially causing code efficiency
      !! problems (e.g. array indexing means that reordering the loops would
      !! require skipping between widely-spaced memory location; potentially
      !! outside those immediately cached)
      !!
      !! OPEN vertical loop
      !CEB fix negative values in DIC and ALK in the nested region
     !  WHERE(trn(:,:,:,jpdic) .lt. 500 .AND. tmask(:,:,:) ==1 )
     !  trn(:,:,:,jpdic)=2050
     !  tra(:,:,:,jpdic)=2050
     !  END WHERE
     !  WHERE(trn(:,:,:,jpalk) .lt. 200 .AND. tmask(:,:,:) ==1 )
     !  trn(:,:,:,jpalk)=2350
    !   tra(:,:,:,jpalk)=2350
    !   END WHERE
    !   WHERE(trn(:,:,:,jpdic) .gt. 3500 .AND. tmask(:,:,:) ==1 )
    !   trn(:,:,:,jpdic)=2050
    !   tra(:,:,:,jpdic)=2050
    !   END WHERE
    !   WHERE(trn(:,:,:,jpalk) .gt. 3200 .AND. tmask(:,:,:) ==1 )
    !   trn(:,:,:,jpalk)=2350
    !   tra(:,:,:,jpalk)=2350
    !   END WHERE
     ! DO jk = 1,jpk
     !    !! OPEN horizontal loops
      !   DO jj = 2,jpjm1
       !     DO ji = 2,jpim1
        !      !! OPEN wet point IF..THEN loop
         !     IF (tmask(ji,jj,jk) == 1) then
          !         IF (trn(ji,jj,jk,jpdic) .lt. 500 .OR. trn(ji,jj,jk,jpdic) .gt. 3500 ) trn(ji,jj,jk,jpdic)=2030
           !        IF (trn(ji,jj,jk,jpalk) .lt. 500 .OR. trn(ji,jj,jk,jpalk) .gt. 3500 ) trn(ji,jj,jk,jpalk)=2450
            !   ENDIF
     !     ENDDO
     !   ENDDO
     ! ENDDO

      !/CEB

      DO jk = 1,jpk
         !! OPEN horizontal loops
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               !! OPEN wet point IF..THEN loop
               if (tmask(ji,jj,jk) == 1) then               
                  !!======================================================
                  !! SETUP LOCAL GRID CELL
                  !!======================================================
                  !!
                  !!------------------------------------------------------
                  !! Some notes on grid vertical structure
                  !! - fsdepw(ji,jj,jk) is the depth of the upper surface of 
                  !!   level jk
                  !! - fsde3w(ji,jj,jk) is *approximately* the midpoint of 
                  !!   level jk
                  !! - fse3t(ji,jj,jk)  is the thickness of level jk
                  !!------------------------------------------------------
                  !!
                  !! AXY (01/03/10): set up level depth (bottom of level)
                  fdep1(ji,jj) = fsdepw(ji,jj,jk) + fse3t(ji,jj,jk)
                  !!
                  !! set up model tracers
                  !! negative values of state variables are not allowed to
                  !! contribute to the calculated fluxes
                  !! non-diatom chlorophyll
                  zchn(ji,jj) = max(0.,trn(ji,jj,jk,jpchn))
                  !! diatom chlorophyll
                  zchd(ji,jj) = max(0.,trn(ji,jj,jk,jpchd))
                  !! non-diatoms
                  zphn(ji,jj) = max(0.,trn(ji,jj,jk,jpphn))
                  !! diatoms
                  zphd(ji,jj) = max(0.,trn(ji,jj,jk,jpphd))
                  !! diatom silicon
                  zpds(ji,jj) = max(0.,trn(ji,jj,jk,jppds))
                  !! AXY (28/01/10): probably need to take account of 
                  !! chl/biomass connection
                  if (zchn(ji,jj).eq.0.) zphn(ji,jj) = 0.
                  if (zchd(ji,jj).eq.0.) zphd(ji,jj) = 0.
                  if (zphn(ji,jj).eq.0.) zchn(ji,jj) = 0.
                  if (zphd(ji,jj).eq.0.) zchd(ji,jj) = 0.
	          !! AXY (23/01/14): duh - why did I forget diatom silicon?
	          if (zpds(ji,jj).eq.0.) zphd(ji,jj) = 0.
	          if (zphd(ji,jj).eq.0.) zpds(ji,jj) = 0.
               ENDIF
            ENDDO
         ENDDO

         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  !! microzooplankton
                  zzmi(ji,jj) = max(0.,trn(ji,jj,jk,jpzmi))
                  !! mesozooplankton
                  zzme(ji,jj) = max(0.,trn(ji,jj,jk,jpzme))
                  !! detrital nitrogen
                  zdet(ji,jj) = max(0.,trn(ji,jj,jk,jpdet))
                  !! dissolved inorganic nitrogen
                  zdin(ji,jj) = max(0.,trn(ji,jj,jk,jpdin))
                  !! dissolved silicic acid
                  zsil(ji,jj) = max(0.,trn(ji,jj,jk,jpsil))
                  !! dissolved "iron"
                  zfer(ji,jj) = max(0.,trn(ji,jj,jk,jpfer))
               ENDIF
            ENDDO
         ENDDO

# if defined key_roam
         !! extra MEDUSA-2 tracers
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  !! detrital carbon
                  zdtc(ji,jj) = max(0.,trn(ji,jj,jk,jpdtc))
                  !! dissolved inorganic carbon
                  zdic(ji,jj) = trn(ji,jj,jk,jpdic)
                  !! alkalinity
                  zalk(ji,jj) = trn(ji,jj,jk,jpalk)
                  !! oxygen
                  zoxy(ji,jj) = max(0.,trn(ji,jj,jk,jpoxy))
#  if defined key_axy_carbchem && defined key_mocsy
                  !! phosphate via DIN and Redfield
                  zpho(ji,jj) = max(0.,trn(ji,jj,jk,jpdin)) / 16.0
#  endif
                  !!
                  !! also need physical parameters for gas exchange 
                  !! calculations
                  ztmp(ji,jj) = tsn(ji,jj,jk,jp_tem)
                  zsal(ji,jj) = tsn(ji,jj,jk,jp_sal)
               ENDIF
            ENDDO
         ENDDO
# else
         !! diagnostic MEDUSA-1 detrital carbon tracer
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  !! implicit detrital carbon
                  zdtc(ji,jj) = zdet(ji,jj) * xthetad
               ENDIF
            ENDDO
         ENDDO
# endif

# if defined key_roam
         !! ---------------------------------------------
         !! JPALM -- 14-12-2017 -- Here, before any exeptionnal crazy value is
         !!              removed, we want to tell to the Master Processor that this 
         !!              Exceptionnal value did exist
         !!

         Call trc_bio_check(kt, jk)
         !!================================================================
	 !! AXY (03/11/17): check input fields
         !! tracer values that exceed thresholds can cause carbonate system
         !! failures when passed to MOCSY; temporary temperature excursions
         !! in recent UKESM0.8 runs trigger such failures but are too short
         !! to have physical consequences; this section checks for such
         !! values and amends them using neighbouring values
         !! 
         !! the check on temperature here is also carried out at the end of
         !! each model time step and anomalies are reported in the master
         !! ocean.output file; the error reporting below is strictly local
         !! to the relevant ocean.output_XXXX file so will not be visible
         !! unless all processors are reporting output
         !!================================================================
         !!
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  !! the thresholds for the four tracers are ...
                  IF ( (ztmp(ji,jj) .LT. -3.0) .OR. (ztmp(ji,jj) .GT. 40.0  ) .OR.   &
                       (zsal(ji,jj) .LE.  0.0) .OR. (zsal(ji,jj) .GT. 50.0  ) .OR.   &
                       (zdic(ji,jj) .LE.  0.0) .OR. (zdic(ji,jj) .GT. 4.0E3 ) .OR.   &
                       (zalk(ji,jj) .LE.  0.0) .OR. (zalk(ji,jj) .GT. 4.0E3 ) ) THEN
                     !!
                     !! all tracer values are reported in the event of any excursion
                     IF (lwp) THEN
                       !CEB  WRITE(charout,*)  ' Tmp = ', ztmp(ji,jj)
                       !CEB  WRITE(charout2,*) ' Sal = ', zsal(ji,jj)
                       !CEB  WRITE(charout3,*) ' DIC = ', zdic(ji,jj)
                       !CEB  WRITE(charout4,*) ' Alk = ', zalk(ji,jj)
                       !CEB  WRITE(charout5,*) mig(ji), mjg(jj), jk, kt 
                         CALL ctl_warn( 'trc_bio_medusa: carbonate chemistry WARNING:',  &
                            TRIM(charout),TRIM(charout2),TRIM(charout3),TRIM(charout4),  & 
                            'at i, j, k, kt:', TRIM(charout5) )
                     ENDIF
                     !!
                     !! Detect, report and correct tracer excursions
                     IF ( (ztmp(ji,jj) .LT. -3.0) .OR. (ztmp(ji,jj) .GT.  40.0) )             &
                        CALL trc_bio_exceptionnal_fix(                                         &
                        tsn(ji-1:ji+1,jj-1:jj+1,jk,jp_tem), tmask(ji-1:ji+1,jj-1:jj+1,jk),    &
                        'Tmp', -3.0, 40.0, ztmp(ji,jj) )
                     !!
                     IF ( (zsal(ji,jj) .LE. 0.0) .OR. (zsal(ji,jj) .GT.  50.0)  )             &
                        CALL trc_bio_exceptionnal_fix(                                         &
                        tsn(ji-1:ji+1,jj-1:jj+1,jk,jp_sal), tmask(ji-1:ji+1,jj-1:jj+1,jk),    &
                        'Sal', 1.0, 50.0, zsal(ji,jj) )
                     !!
                     IF ( (zdic(ji,jj) .LE. 0.0) .OR. (zdic(ji,jj) .GT. 4.0E3)  )             &
                        CALL trc_bio_exceptionnal_fix(                                         &
                        trn(ji-1:ji+1,jj-1:jj+1,jk,jpdic), tmask(ji-1:ji+1,jj-1:jj+1,jk),     &
                        'DIC', 100.0, 4.0E3, zdic(ji,jj) )
                     !!
                     IF ( (zalk(ji,jj) .LE. 0.0) .OR. (zalk(ji,jj) .GT. 4.0E3)  )             &
                        CALL trc_bio_exceptionnal_fix(                                         &
                        trn(ji-1:ji+1,jj-1:jj+1,jk,jpalk), tmask(ji-1:ji+1,jj-1:jj+1,jk),     &
                        'Alk', 100.0, 4.0E3, zalk(ji,jj) )
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
# endif

# if defined key_debug_medusa
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,jk) == 1) then
                  if (idf.eq.1) then
                     !! AXY (15/01/10)
                     if (trn(ji,jj,jk,jpdin).lt.0.) then
                        IF (lwp) write (numout,*)                            &
                           '------------------------------'
                        IF (lwp) write (numout,*) 'NEGATIVE DIN ERROR =',    &
                           trn(ji,jj,jk,jpdin)
                        IF (lwp) write (numout,*) 'NEGATIVE DIN ERROR @',    &
                           ji, jj, jk, kt
                     endif
                     if (trn(ji,jj,jk,jpsil).lt.0.) then
                        IF (lwp) write (numout,*)                            &
                           '------------------------------'
                        IF (lwp) write (numout,*) 'NEGATIVE SIL ERROR =',    &
                           trn(ji,jj,jk,jpsil)
                        IF (lwp) write (numout,*) 'NEGATIVE SIL ERROR @',    &
                           ji, jj, jk, kt
                     endif
#  if defined key_roam
                     if (trn(ji,jj,jk,jpdic).lt.0.) then
                        IF (lwp) write (numout,*)                            &
                           '------------------------------'
                        IF (lwp) write (numout,*) 'NEGATIVE DIC ERROR =',    &
                           trn(ji,jj,jk,jpdic)
                        IF (lwp) write (numout,*) 'NEGATIVE DIC ERROR @',    &
                           ji, jj, jk, kt
                     endif
                     if (trn(ji,jj,jk,jpalk).lt.0.) then
                        IF (lwp) write (numout,*)                            &
                           '------------------------------'
                        IF (lwp) write (numout,*) 'NEGATIVE ALK ERROR =',    &
                           trn(ji,jj,jk,jpalk)
                        IF (lwp) write (numout,*) 'NEGATIVE ALK ERROR @',    &
                           ji, jj, jk, kt
                     endif
                     if (trn(ji,jj,jk,jpoxy).lt.0.) then
                        IF (lwp) write (numout,*)                            &
                           '------------------------------'
                        IF (lwp) write (numout,*) 'NEGATIVE OXY ERROR =',    &
                           trn(ji,jj,jk,jpoxy)
                        IF (lwp) write (numout,*) 'NEGATIVE OXY ERROR @',    &
                           ji, jj, jk, kt
                     endif
#  endif
                  endif
               ENDIF
            ENDDO
         ENDDO
# endif
# if defined key_debug_medusa
! I'M NOT SURE THIS IS USEFUL NOW THAT I'VE SPLIT THE DO LOOP - marc 8/5/17
!         if (idf.eq.1.AND.idfval.eq.1) then
!            DO jj = 2,jpjm1
!               DO ji = 2,jpim1
!                  if (tmask(ji,jj,jk) == 1) then
!                     !! report state variable values
!                     IF (lwp) write (numout,*)                               &
!                        '------------------------------'
!                     IF (lwp) write (numout,*) 'fthk(',jk,') = ',            &
!                        fse3t(ji,jj,jk)
!                     IF (lwp) write (numout,*) 'zphn(',jk,') = ', zphn(ji,jj)
!                     IF (lwp) write (numout,*) 'zphd(',jk,') = ', zphd(ji,jj)
!                     IF (lwp) write (numout,*) 'zpds(',jk,') = ', zpds(ji,jj)
!                     IF (lwp) write (numout,*) 'zzmi(',jk,') = ', zzmi(ji,jj)
!                     IF (lwp) write (numout,*) 'zzme(',jk,') = ', zzme(ji,jj)
!                     IF (lwp) write (numout,*) 'zdet(',jk,') = ', zdet(ji,jj)
!                     IF (lwp) write (numout,*) 'zdin(',jk,') = ', zdin(ji,jj)
!                     IF (lwp) write (numout,*) 'zsil(',jk,') = ', zsil(ji,jj)
!                     IF (lwp) write (numout,*) 'zfer(',jk,') = ', zfer(ji,jj)
#  if defined key_roam
!                     IF (lwp) write (numout,*) 'zdtc(',jk,') = ', zdtc(ji,jj)
!                     IF (lwp) write (numout,*) 'zdic(',jk,') = ', zdic(ji,jj)
!                     IF (lwp) write (numout,*) 'zalk(',jk,') = ', zalk(ji,jj)
!                     IF (lwp) write (numout,*) 'zoxy(',jk,') = ', zoxy(ji,jj)
#  endif
!                  ENDIF
!               ENDDO
!            ENDDO
!         endif
# endif

# if defined key_debug_medusa
! I'M NOT SURE THIS IS USEFUL NOW THAT I'VE SPLIT THE DO LOOP - marc 8/5/17
!         if (idf.eq.1.AND.idfval.eq.1.AND.jk.eq.1) then
!            DO jj = 2,jpjm1
!               DO ji = 2,jpim1
!                  if (tmask(ji,jj,jk) == 1) then
!                     IF (lwp) write (numout,*)                               &
!                       '------------------------------'
!                     IF (lwp) write (numout,*) 'dust      = ', dust(ji,jj)
!                  ENDIF
!               ENDDO
!            ENDDO
!         endif
# endif

         !!---------------------------------------------------------------
         !! Calculate air-sea gas exchange and river inputs
         !!---------------------------------------------------------------
         IF ( jk == 1 ) THEN
            CALL air_sea( kt )
         ENDIF

         !!---------------------------------------------------------------
         !! Phytoplankton growth, zooplankton grazing and miscellaneous
         !! plankton losses. 
         !!---------------------------------------------------------------
         CALL plankton( jk )

         !!---------------------------------------------------------------
         !! Iron chemistry and scavenging
         !!---------------------------------------------------------------
         CALL iron_chem_scav( jk )

         !!---------------------------------------------------------------
         !! Detritus processes
         !!---------------------------------------------------------------
         CALL detritus( jk, iball )

         !!---------------------------------------------------------------
         !! Updating tracers
         !!---------------------------------------------------------------
         CALL bio_medusa_update( kt, jk )

         !!---------------------------------------------------------------
         !! Diagnostics
         !!---------------------------------------------------------------
         CALL bio_medusa_diag( jk )

         !!-------------------------------------------------------
         !! 2d specific k level diags
         !!-------------------------------------------------------
         IF( lk_iomput ) THEN
            CALL bio_medusa_diag_slice( jk )
         ENDIF

      !! CLOSE vertical loop
      ENDDO

      !!------------------------------------------------------------------
      !! Final calculations for diagnostics
      !!------------------------------------------------------------------
      CALL bio_medusa_fin( kt )

# if defined key_debug_medusa
       IF(lwp) WRITE(numout,*) ' MEDUSA exiting trc_bio_medusa at kt =', kt
       CALL flush(numout)
# endif

   END SUBROUTINE trc_bio_medusa



   SUBROUTINE trc_bio_exceptionnal_fix(tiny_var, tiny_mask, var_nm, mini, maxi, varout)
      !! JPALM (27/10/17): This function is called only when abnormal values that 
      !! could break the model's carbonate system are fed to MEDUSA
      !! 
      !! The function calculates an average tracer value based on the 3x3 cell
      !! neighbourhood around the abnormal cell, and reports this back
      !! 
      !! Tracer array values are not modified, but MEDUSA uses "corrected" values
      !! in its calculations
      !!
      !! temporary variables
      REAL(wp), INTENT( in ), DIMENSION(3,3) :: tiny_var, tiny_mask
      CHARACTER(3), INTENT( in )            :: var_nm
      REAL(wp), INTENT( in )                 :: mini, maxi
      REAL(wp), INTENT( out )                :: varout
      REAL(wp)      :: sumtsn, tsnavg
      INTEGER       :: summask
      CHARACTER(25) :: charout1, charout2
      CHARACTER(60) :: charout3, charout4
      INTEGER       :: ii,ij
    
      !! point to the center of the 3*3 zoom-grid, to check around
      ii = 2
      ij = 2
      !! Print surounding values to check if isolated Crazy value or 
      !! General error
      IF(lwp) THEN 
          WRITE(numout,*)                                 &
            '----------------------------------------------------------------------'
          WRITE(numout,*)                                 &
            'trc_bio_medusa: 3x3 neighbourhood surrounding abnormal ', TRIM(var_nm)
          WRITE(numout,9100)                              &
            3, tiny_var(ii-1,ij+1), tiny_var(ii  ,ij+1), tiny_var(ii+1,ij+1)
          WRITE(numout,9100)                              &
            2, tiny_var(ii-1,ij  ), tiny_var(ii  ,ij  ), tiny_var(ii+1,ij  )
          WRITE(numout,9100)                              &
            1, tiny_var(ii-1,ij-1), tiny_var(ii  ,ij-1), tiny_var(ii+1,ij-1)
          WRITE(numout,*)                                 &
            'trc_bio_medusa: 3x3 land-sea neighbourhood, tmask'
          WRITE(numout,9100)                              &
            3, tiny_mask(ii-1,ij+1), tiny_mask(ii  ,ij+1), tiny_mask(ii+1,ij+1)
          WRITE(numout,9100)                              &
            2, tiny_mask(ii-1,ij  ), tiny_mask(ii  ,ij  ), tiny_mask(ii+1,ij  )
          WRITE(numout,9100)                              &
            1, tiny_mask(ii-1,ij-1), tiny_mask(ii  ,ij-1), tiny_mask(ii+1,ij-1)
      ENDIF
      !! Correct out of range values
      sumtsn = ( tiny_mask(ii-1,ij+1) * tiny_var(ii-1,ij+1) ) +  &
               ( tiny_mask(ii  ,ij+1) * tiny_var(ii  ,ij+1) ) +  &
               ( tiny_mask(ii+1,ij+1) * tiny_var(ii+1,ij+1) ) +  &
               ( tiny_mask(ii-1,ij  ) * tiny_var(ii-1,ij  ) ) +  &
               ( tiny_mask(ii+1,ij  ) * tiny_var(ii+1,ij  ) ) +  &
               ( tiny_mask(ii-1,ij-1) * tiny_var(ii-1,ij-1) ) +  &
               ( tiny_mask(ii  ,ij-1) * tiny_var(ii  ,ij-1) ) +  &
               ( tiny_mask(ii+1,ij-1) * tiny_var(ii+1,ij-1) )
      !!
      summask = tiny_mask(ii-1,ij+1) + tiny_mask(ii  ,ij+1)   +  &
                tiny_mask(ii+1,ij+1) + tiny_mask(ii-1,ij  )   +  &
                tiny_mask(ii+1,ij  ) + tiny_mask(ii-1,ij-1)   +  &
                tiny_mask(ii  ,ij-1) + tiny_mask(ii+1,ij-1)
      !!
      IF ( summask .GT. 0 ) THEN
         tsnavg = ( sumtsn / summask )
         varout = MAX( MIN( maxi, tsnavg), mini )
      ELSE   
         IF (ztmp(ii,ij) .LT. mini )  varout = mini
         IF (ztmp(ii,ij) .GT. maxi )  varout = maxi
      ENDIF
      !!
      IF (lwp) THEN 
          WRITE(charout1,9200) tiny_var(ii,ij)
          WRITE(charout2,9200) varout
          WRITE(charout3,*) ' ', charout1, ' -> ', charout2
          WRITE(charout4,*) ' Tracer: ', trim(var_nm)
      !!
          WRITE(numout,*) 'trc_bio_medusa: ** EXCEPTIONAL VALUE SWITCHING **'
          WRITE(numout,*) charout4 
          WRITE(numout,*) charout3
          WRITE(numout,*) '----------------------------------------------------------------------'
          WRITE(numout,*) ' '
      ENDIF

9100  FORMAT('Row:', i1, '  ', e16.6, ' ', e16.6, ' ', e16.6)
9200  FORMAT(e16.6)

   END SUBROUTINE trc_bio_exceptionnal_fix 

   SUBROUTINE trc_bio_check(kt, jk)
      !!-----------------------------------
      !! JPALM -- 14-12-2017 -- Still dealing with this micro-boil/carb failure
      !!                     problem. The model is now able to correct a local
      !!                     crazy value. but does it silently.
      !!                     We need to spread the word to the master processor. we
      !!                     don't want the model  to correct values without telling us
      !!                     This module will tell at least when crazy DIC or
      !!                     ALK appears.
      !!-------------------------------------
      REAL(wp)              :: zmax, zmin    ! temporary scalars
      INTEGER               :: ji,jj         ! dummy loop
      INTEGER               :: ii,ij         ! temporary scalars 
      INTEGER, DIMENSION(2) :: ilocs         ! 
      INTEGER, INTENT( in ) :: kt, jk
      !!
      !!==========================
      !! DIC Check
      zmax =  -5.0  ! arbitrary  low maximum value
      zmin =  4.0E4  ! arbitrary high minimum value
      DO jj = 2, jpjm1
         DO ji = 2,jpim1
            IF( tmask(ji,jj,1) == 1) THEN
               zmax = MAX(zmax,zdic(ji,jj))     ! find local maximum
               zmin = MIN(zmin,zdic(ji,jj))     ! find local minimum
            ENDIF
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( zmax )       ! max over the global domain
      IF( lk_mpp )   CALL mpp_min( zmin )       ! min over the global domain
      !
      IF( zmax .GT. 4.0E3) THEN  ! we've got a problem
         IF (lk_mpp) THEN
            CALL mpp_maxloc ( zdic(:,:),tmask(:,:,1), zmax, ii,ij )
         ELSE
            ilocs = MAXLOC( zdic(:,:), mask = tmask(:,:,1) == 1. )
            ii = ilocs(1) + nimpp - 1
            ij = ilocs(2) + njmpp - 1
         ENDIF
         !
         IF(lwp) THEN
            WRITE(numout,*) 'trc_bio:tracer anomaly: *****    WARNING     *****'
            WRITE(numout,*) 'trc_bio:tracer anomaly: DIC concentration > 4000 '
            WRITE(numout,9600) kt, zmax, ii, ij, jk
            WRITE(numout,*) 'trc_bio:tracer anomaly: ***** END OF WARNING *****'
         ENDIF
      ENDIF
      !
      IF( zmin .LE. 0.0) THEN  ! we've got a problem
         IF (lk_mpp) THEN
            CALL mpp_minloc ( zdic(:,:),tmask(:,:,1), zmin, ii,ij )
         ELSE
            ilocs = MINLOC( zdic(:,:), mask = tmask(:,:,1) == 1. )
            ii = ilocs(1) + nimpp - 1
            ij = ilocs(2) + njmpp - 1
         ENDIF
         !
         IF(lwp) THEN
            WRITE(numout,*) 'trc_bio:tracer anomaly: *****    WARNING     *****'
            WRITE(numout,*) 'trc_bio:tracer anomaly: DIC concentration <= 0 '
            WRITE(numout,9700) kt, zmin, ii, ij, jk
            WRITE(numout,*) 'trc_bio:tracer anomaly: ***** END OF WARNING *****'
         ENDIF
      ENDIF
      !!
      !!==========================
      !! ALKALINITY Check
      zmax =  -5.0  ! arbitrary  low maximum value
      zmin =  4.0E4  ! arbitrary high minimum value
      DO jj = 2, jpjm1
         DO ji = 2,jpim1
            IF( tmask(ji,jj,1) == 1) THEN
               zmax = MAX(zmax,zalk(ji,jj))     ! find local maximum
               zmin = MIN(zmin,zalk(ji,jj))     ! find local minimum
            ENDIF
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( zmax )       ! max over the global domain
      IF( lk_mpp )   CALL mpp_min( zmin )       ! min over the global domain
      !
      IF( zmax .GT. 4.0E3) THEN  ! we've got a problem
         IF (lk_mpp) THEN
            CALL mpp_maxloc ( zalk(:,:),tmask(:,:,1), zmax, ii,ij )
         ELSE
            ilocs = MAXLOC( zalk(:,:), mask = tmask(:,:,1) == 1. )
            ii = ilocs(1) + nimpp - 1
            ij = ilocs(2) + njmpp - 1
         ENDIF
         !
         IF(lwp) THEN
            WRITE(numout,*) 'trc_bio:tracer anomaly: *****     WARNING    *****'
            WRITE(numout,*) 'trc_bio:tracer anomaly: ALK concentration > 4000 '
            WRITE(numout,9800) kt, zmax, ii, ij, jk
            WRITE(numout,*) 'trc_bio:tracer anomaly: ***** END OF WARNING *****'
         ENDIF
      ENDIF
      !
      IF( zmin .LE. 0.0) THEN  ! we've got a problem
         IF (lk_mpp) THEN
            CALL mpp_minloc ( zalk(:,:),tmask(:,:,1), zmin, ii,ij )
         ELSE
            ilocs = MINLOC( zalk(:,:), mask = tmask(:,:,1) == 1. )
            ii = ilocs(1) + nimpp - 1
            ij = ilocs(2) + njmpp - 1
         ENDIF
         !
         IF(lwp) THEN
            WRITE(numout,*) 'trc_bio:tracer anomaly: *****    WARNING     *****'
            WRITE(numout,*) 'trc_bio:tracer anomaly:  ALK concentration <= 0 '
            WRITE(numout,9900) kt, zmin, ii, ij, jk
            WRITE(numout,*) 'trc_bio:tracer anomaly: ***** END OF WARNING *****'
         ENDIF
      ENDIF


9600  FORMAT ('trc_bio:tracer anomaly: kt=',i6,' max DIC: ',f16.10,', i j k: ',3i5)
9700  FORMAT ('trc_bio:tracer anomaly: kt=',i6,' min DIC: ',f16.10,', i j k: ',3i5)
9800  FORMAT ('trc_bio:tracer anomaly: kt=',i6,' max ALK: ',f16.10,', i j k: ',3i5)
9900  FORMAT ('trc_bio:tracer anomaly: kt=',i6,' min ALK: ',f16.10,', i j k: ',3i5)

   END SUBROUTINE trc_bio_check


#else
   !!=====================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!=====================================================================
CONTAINS
   SUBROUTINE trc_bio_medusa( kt )                   ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_bio_medusa: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bio_medusa
#endif 

   !!=====================================================================
END MODULE trcbio_medusa
