MODULE bio_medusa_fin_mod
   !!======================================================================
   !!                         ***  MODULE bio_medusa_fin_mod  ***
   !! Finalisation for TRC_BIO_MEDUSA
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!   -   ! 2017-08 (A. Yool)            Amend bethic reservoir updating
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   bio_medusa_fin     ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE bio_medusa_fin( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE bio_medusa_fin  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - ...
      !!----------------------------------------------------------------------
      USE bio_medusa_mod
      USE dom_oce,           ONLY: atfp, atfp1, neuler, rdt, tmask
      USE in_out_manager,    ONLY: lwp, numout
      USE iom,               ONLY: iom_put
      USE lbclnk,            ONLY: lbc_lnk
      USE oce,               ONLY: chloro_out_cpl 
      USE par_medusa,        ONLY: jp_medusa_2d, jp_medusa_3d,          &
                                   jp_medusa_trd, jpchd, jpchn
      USE par_oce,           ONLY: jpi, jpim1, jpj, jpjm1, jpk
      USE phycst,            ONLY: rsmall
      USE sbc_oce,           ONLY: lk_oasis
      USE sms_medusa,        ONLY: jinorgben, jorgben,                  &
                                   f3_co3, f3_h2co3, f3_hco3,           &
                                   f3_omarg, f3_omcal, f3_pH,           &
                                   za_sed_c, za_sed_ca, za_sed_fe,      &
                                   za_sed_n, za_sed_si,                 &
                                   zb_sed_c, zb_sed_ca, zb_sed_fe,      &
                                   zb_sed_n, zb_sed_si,                 &
                                   zn_sed_c, zn_sed_ca, zn_sed_fe,      &
                                   zn_sed_n, zn_sed_si, zn_chl_srf,     &
                                   scl_chl, chl_out
      !CEB USE trc,               ONLY: med_diag, nittrc000, trn 
      USE trc
      !/CEB
      USE trcnam_trp,        ONLY: ln_trcadv_cen2, ln_trcadv_tvd
 
      !! time (integer timestep)
      INTEGER, INTENT( in ) ::    kt

      INTEGER :: ji, jj
      INTEGER :: jn

      REAL(wp) :: fq0,fq1,fq2,fq3
      !CEB!$AGRIF_DO_NOT_TREAT
# if defined key_roam                     
      !!----------------------------------------------------------------------
      !! AXY (09/08/17): fix benthic submodel
      !!----------------------------------------------------------------------
      !! Process benthic in/out fluxes
      !! These can be handled outside of the 3D calculations since the
      !! benthic pools (and fluxes) are 2D in nature; this code was
      !! developed with help from George Nurser (NOC); it cannot be run
      !! in a configuration with variable time-stepping with depth
      !!----------------------------------------------------------------------
      !!
      !! time-step calculation
      if (jorgben.eq.1) then
         za_sed_n(:,:)  = zb_sed_n(:,:)  + ((2. * (rdt / 86400.)) * &
                          ( f_sbenin_n(:,:)  + f_fbenin_n(:,:)  - f_benout_n(:,:)  ))
         za_sed_fe(:,:) = zb_sed_fe(:,:) + ((2. * (rdt / 86400.)) * &
                          ( f_sbenin_fe(:,:) + f_fbenin_fe(:,:) - f_benout_fe(:,:) ))
         za_sed_c(:,:)  = zb_sed_c(:,:)  + ((2. * (rdt / 86400.)) * &
                          ( f_sbenin_c(:,:)  + f_fbenin_c(:,:)  - f_benout_c(:,:)  ))
      endif
      if (jinorgben.eq.1) then
         za_sed_si(:,:) = zb_sed_si(:,:) + ((2. * (rdt / 86400.)) * &
                          ( f_fbenin_si(:,:) - f_benout_si(:,:) ))
         za_sed_ca(:,:) = zb_sed_ca(:,:) + ((2. * (rdt / 86400.)) * &
                          ( f_fbenin_ca(:,:) - f_benout_ca(:,:) ))
      endif
      !!
      !! time-level calculation
      if (jorgben.eq.1) then
         zb_sed_n(:,:)  = zn_sed_n(:,:)  + (atfp * &
                          ( za_sed_n(:,:)  - (2. * zn_sed_n(:,:))  + zb_sed_n(:,:)  ))
         zn_sed_n(:,:)  = za_sed_n(:,:)
         zb_sed_fe(:,:) = zn_sed_fe(:,:) + (atfp * &
                          ( za_sed_fe(:,:) - (2. * zn_sed_fe(:,:)) + zb_sed_fe(:,:) ))
         zn_sed_fe(:,:) = za_sed_fe(:,:)
         zb_sed_c(:,:)  = zn_sed_c(:,:)  + (atfp * &
                          ( za_sed_c(:,:)  - (2. * zn_sed_c(:,:))  + zb_sed_c(:,:)  ))
         zn_sed_c(:,:)  = za_sed_c(:,:)
      endif
      if (jinorgben.eq.1) then
         zb_sed_si(:,:) = zn_sed_si(:,:) + (atfp * &
                          ( za_sed_si(:,:) - (2. * zn_sed_si(:,:)) + zb_sed_si(:,:) ))
         zn_sed_si(:,:) = za_sed_si(:,:)
         zb_sed_ca(:,:) = zn_sed_ca(:,:) + (atfp * &
                          ( za_sed_ca(:,:) - (2. * zn_sed_ca(:,:)) + zb_sed_ca(:,:) ))
         zn_sed_ca(:,:) = za_sed_ca(:,:)
      endif
# endif      

#  if defined key_debug_medusa
         !! AXY (12/07/17)
         !!-----------------------------------------------------------------
         !! Check conservation of MEDUSA's sinks-minus-sources using fflx_X 
         !! diagnostics (i.e. biogeochemical processes only)
         !!   - fflx_X diagnostics *should* include all transfers between
         !!     modelled components
         !!   - they should also include gains / losses due to air-sea 
         !!     fluxes of C and O2, aeolian and seafloor inputs of Fe, and 
         !!     inputs from seafloor "benthic buckets" (N, Si, Fe, C and 
         !!     alkalinity)
         !!   - however, they do not include the transfer of material to 
         !!     "benthic buckets" by sedimenting slow- and fast-sinking
         !!     detritus since these are separate 2D reservoirs
         !!   - consequently, for a given water column, the integrated
         !!     fluxes should sum to the "loss" of material to the "benthic
         !!     buckets"
         !!   - if they do not, this suggests that MEDUSA contains errors
         !!     in its accounting (e.g. processes omitted from calculated
         !!     fluxes)
         !!   - here, the local integrated fluxes and benthic inputs (plus
         !!     air-sea fluxes in the case of C) are reported together with
         !!     the resulting error
         !!   - only N, Si, C and alkalinity inventories considered; Fe and
         !!     O2 overlooked because of wholesale loss (and addition, in 
         !!     the case of O2) of these tracers within the water column
         !!-----------------------------------------------------------------
         !!
         !! nitrogen
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,1) == 1) then
                  fq0 = fflx_n(ji,jj)
                  fq1 = f_sbenin_n(ji,jj) + f_fbenin_n(ji,jj)
                  fq2 = fq0 + fq1
                  fq3 = f_benout_n(ji,jj)
                  if (lwp) write (numout,'a,2i3,a,4f15,5)')                   &
                     'AXY N   cons: (i,j)=',ji,jj,', (flx,ben,err,out)=',      &
                     fq0,fq1,fq2,fq3
               ENDIF
            ENDDO
         ENDDO   
         !! silicon
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,1) == 1) then
                  fq0 = fflx_si(ji,jj)
                  fq1 = f_fbenin_si(ji,jj)
                  fq2 = fq0 + fq1
                  fq3 = f_benout_si(ji,jj)
                  if (lwp) write (numout,'a,2i3,a,4f15,5)')                   &
                     'AXY Si  cons: (i,j)=',ji,jj,', (flx,ben,err,out)=',     &
                     fq0,fq1,fq2,fq3
               ENDIF
            ENDDO
         ENDDO   
         !! carbon
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               if (tmask(ji,jj,1) == 1) then
                  fq0 = fflx_c(ji,jj)
                  fq1 = f_sbenin_c(ji,jj) + f_fbenin_c(ji,jj) + f_fbenin_ca(ji,jj)
                  fq2 = f_co2flux(ji,jj) * fse3t(ji,jj,1)
                  fq3 = fq0 + fq1
                  fq4 = f_benout_c(ji,jj) + f_benout_ca(ji,jj)
                  if (lwp) write (numout,'a,2i3,a,5f15,5)')                   &
                     'AXY C   cons: (i,j)=',ji,jj,', (flx,ben,asf,err,out)=', &
                     fq0,fq1,fq2,fq3,fq4
                ENDIF
             ENDDO
          ENDDO   
          !! alkalinity
          DO jj = 2,jpjm1
             DO ji = 2,jpim1
                if (tmask(ji,jj,1) == 1) then
                   fq0 = fflx_a(ji,jj)
                   fq1 = 2.0 * f_fbenin_ca(ji,jj)
                   fq2 = fq0 + fq1
                   fq3 = 2.0 * f_benout_ca(ji,jj)
                   if (lwp) write (numout,'a,2i3,a,4f15,5)')                   &
                      'AXY alk cons: (i,j)=',ji,jj,', (flx,ben,err,out)=',     &
                      fq0,fq1,fq2,fq3
               ENDIF
            ENDDO
         ENDDO   
#  endif

         !!!---------------------------------------------------------------
         !! Add very last diag calculations 
         !!!---------------------------------------------------------------
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               !!         
               IF( med_diag%PN_JLIM%dgsave ) THEN
                  fjln2d(ji,jj) = fjln2d(ji,jj)   / MAX(ftot_pn(ji,jj), rsmall)
               ENDIF
               IF( med_diag%PN_NLIM%dgsave ) THEN
                  fnln2d(ji,jj) = fnln2d(ji,jj)   / MAX(ftot_pn(ji,jj), rsmall)
               ENDIF
               IF( med_diag%PN_FELIM%dgsave ) THEN
                  ffln2d(ji,jj) = ffln2d(ji,jj)   / MAX(ftot_pn(ji,jj), rsmall)
               ENDIF
               IF( med_diag%PD_JLIM%dgsave ) THEN
                  fjld2d(ji,jj) = fjld2d(ji,jj)   / MAX(ftot_pd(ji,jj), rsmall)
               ENDIF
               IF( med_diag%PD_NLIM%dgsave ) THEN
                  fnld2d(ji,jj) = fnld2d(ji,jj)   / MAX(ftot_pd(ji,jj), rsmall)
               ENDIF
               IF( med_diag%PD_FELIM%dgsave ) THEN
                  ffld2d(ji,jj) = ffld2d(ji,jj)   / MAX(ftot_pd(ji,jj), rsmall)
               ENDIF
               IF( med_diag%PD_SILIM%dgsave ) THEN
                  fsld2d2(ji,jj) = fsld2d2(ji,jj) / MAX(ftot_pd(ji,jj), rsmall)
               ENDIF
               IF( med_diag%PDSILIM2%dgsave ) THEN
                  fsld2d(ji,jj) = fsld2d(ji,jj)   / MAX(ftot_pd(ji,jj), rsmall)
               ENDIF
            ENDDO
         ENDDO

         !!!---------------------------------------------------------------
         !! Calculates Chl diag for UM coupling 
         !!!---------------------------------------------------------------
         !! JPALM -- 02-06-2017 --
         !! add Chl surf coupling
         !! no need to output, just pass to cpl var
         IF (lk_oasis) THEN
            IF (chl_out.eq.1) THEN
               !! export and scale surface chl
               zn_chl_srf(:,:) = MAX( 0.0, (trn(:,:,1,jpchd) + trn(:,:,1,jpchn)) * 1.0E-6 )
                                 !! surf Chl in Kg-chl/m3 as needed for cpl
            ELSEIF (chl_out.eq.2) THEN
               !! export and scale mld chl
               zn_chl_srf(:,:) = MAX( 0.0, fchl_ml(:,:) * 1.0E-6 )
                                 !! mld Chl in Kg-chl/m3 as needed for cpl
            ENDIF
            chloro_out_cpl(:,:) = zn_chl_srf(:,:) * scl_chl        !! Coupling Chl
         END IF

         !!----------------------------------------------------------------
         !! Add in XML diagnostics stuff
         !!----------------------------------------------------------------
         !!
         !! ** 2D diagnostics
#   if defined key_debug_medusa
         IF (lwp) write (numout,*) 'bio_medusa_fin: export all diag kt = ', kt
         CALL flush(numout)
#   endif
         IF ( med_diag%INVTN%dgsave ) THEN
            CALL iom_put( "INVTN"  , ftot_n )
         ENDIF
         IF ( med_diag%INVTSI%dgsave ) THEN
            CALL iom_put( "INVTSI"  , ftot_si )
         ENDIF
         IF ( med_diag%INVTFE%dgsave ) THEN
            CALL iom_put( "INVTFE"  , ftot_fe )
         ENDIF                           
         IF ( med_diag%ML_PRN%dgsave ) THEN
            CALL iom_put( "ML_PRN"  , fprn_ml )
         ENDIF
         IF ( med_diag%ML_PRD%dgsave ) THEN
            CALL iom_put( "ML_PRD"  , fprd_ml )
         ENDIF
         IF ( med_diag%OCAL_LVL%dgsave ) THEN
            CALL iom_put( "OCAL_LVL"  , fccd )
         ENDIF
         IF ( med_diag%CHL_MLD%dgsave ) THEN
            CALL iom_put( "CHL_MLD"  , fchl_ml )
         ENDIF
         IF (lk_oasis) THEN
            IF ( med_diag%CHL_CPL%dgsave ) THEN
               CALL iom_put( "CHL_CPL"  , chloro_out_cpl )
            ENDIF
         ENDIF
         IF ( med_diag%PN_JLIM%dgsave ) THEN
            CALL iom_put( "PN_JLIM"  , fjln2d )
            DEALLOCATE( fjln2d )
         ENDIF
         IF ( med_diag%PN_NLIM%dgsave ) THEN
            CALL iom_put( "PN_NLIM"  , fnln2d )
            DEALLOCATE( fnln2d )
         ENDIF
         IF ( med_diag%PN_FELIM%dgsave ) THEN
            CALL iom_put( "PN_FELIM"  , ffln2d )
            DEALLOCATE( ffln2d )
         ENDIF
         IF ( med_diag%PD_JLIM%dgsave ) THEN
            CALL iom_put( "PD_JLIM"  , fjld2d )
            DEALLOCATE( fjld2d )
         ENDIF
         IF ( med_diag%PD_NLIM%dgsave ) THEN
            CALL iom_put( "PD_NLIM"  , fnld2d )
            DEALLOCATE( fnld2d )
         ENDIF
         IF ( med_diag%PD_FELIM%dgsave ) THEN
            CALL iom_put( "PD_FELIM"  , ffld2d )
            DEALLOCATE( ffld2d )
         ENDIF
         IF ( med_diag%PD_SILIM%dgsave ) THEN
            CALL iom_put( "PD_SILIM"  , fsld2d2 )
            DEALLOCATE( fsld2d2 )
         ENDIF
         IF ( med_diag%PDSILIM2%dgsave ) THEN
            CALL iom_put( "PDSILIM2"  , fsld2d )
            DEALLOCATE( fsld2d )
         ENDIF
         IF ( med_diag%INTFLX_N%dgsave ) THEN
            CALL iom_put( "INTFLX_N"  , fflx_n )
         ENDIF
         IF ( med_diag%INTFLX_SI%dgsave ) THEN
            CALL iom_put( "INTFLX_SI"  , fflx_si )
         ENDIF
         IF ( med_diag%INTFLX_FE%dgsave ) THEN
            CALL iom_put( "INTFLX_FE"  , fflx_fe )
         ENDIF        
         IF ( med_diag%INT_PN%dgsave ) THEN
            CALL iom_put( "INT_PN"  , ftot_pn )
         ENDIF
         IF ( med_diag%INT_PD%dgsave ) THEN
            CALL iom_put( "INT_PD"  , ftot_pd )
         ENDIF         
         IF ( med_diag%INT_ZMI%dgsave ) THEN
            CALL iom_put( "INT_ZMI"  , ftot_zmi )
         ENDIF
         IF ( med_diag%INT_ZME%dgsave ) THEN
            CALL iom_put( "INT_ZME"  , ftot_zme )
         ENDIF
         IF ( med_diag%INT_DET%dgsave ) THEN
            CALL iom_put( "INT_DET"  , ftot_det )
         ENDIF
         IF ( med_diag%INT_DTC%dgsave ) THEN
            CALL iom_put( "INT_DTC"  , ftot_dtc )
         ENDIF
         IF ( med_diag%BEN_N%dgsave ) THEN
            CALL iom_put( "BEN_N"  , za_sed_n )
         ENDIF
         IF ( med_diag%BEN_FE%dgsave ) THEN
            CALL iom_put( "BEN_FE"  , za_sed_fe )
         ENDIF
         IF ( med_diag%BEN_C%dgsave ) THEN
            CALL iom_put( "BEN_C"  , za_sed_c )
         ENDIF
         IF ( med_diag%BEN_SI%dgsave ) THEN
            CALL iom_put( "BEN_SI"  , za_sed_si )
         ENDIF
         IF ( med_diag%BEN_CA%dgsave ) THEN
            CALL iom_put( "BEN_CA"  , za_sed_ca )
         ENDIF
         IF ( med_diag%RUNOFF%dgsave ) THEN
            CALL iom_put( "RUNOFF"  , f_runoff )
         ENDIF 
# if defined key_roam        
         IF ( med_diag%N_PROD%dgsave ) THEN
            CALL iom_put( "N_PROD"  , fnit_prod )
         ENDIF
         IF ( med_diag%N_CONS%dgsave ) THEN
            CALL iom_put( "N_CONS"  , fnit_cons )
         ENDIF
         IF ( med_diag%C_PROD%dgsave ) THEN
            CALL iom_put( "C_PROD"  , fcar_prod )
         ENDIF
         IF ( med_diag%C_CONS%dgsave ) THEN
            CALL iom_put( "C_CONS"  , fcar_cons )
         ENDIF
         IF ( med_diag%O2_PROD%dgsave ) THEN
            CALL iom_put( "O2_PROD"  , foxy_prod )
         ENDIF
         IF ( med_diag%O2_CONS%dgsave ) THEN
            CALL iom_put( "O2_CONS"  , foxy_cons )
         ENDIF
         IF ( med_diag%O2_ANOX%dgsave ) THEN
            CALL iom_put( "O2_ANOX"  , foxy_anox )
         ENDIF
         IF ( med_diag%INVTC%dgsave ) THEN
            CALL iom_put( "INVTC"  , ftot_c )
         ENDIF
         IF ( med_diag%INVTALK%dgsave ) THEN
            CALL iom_put( "INVTALK"  , ftot_a )
         ENDIF
         IF ( med_diag%INVTO2%dgsave ) THEN
            CALL iom_put( "INVTO2"  , ftot_o2 )
         ENDIF
         IF ( med_diag%COM_RESP%dgsave ) THEN
            CALL iom_put( "COM_RESP"  , fcomm_resp )
         ENDIF         
# endif      
         !!
         !! diagnostic filled in the i-j-k main loop
         !!--------------------------------------------
         IF ( med_diag%PRN%dgsave ) THEN
            CALL iom_put( "PRN"  , fprn2d )
            DEALLOCATE( fprn2d )
         ENDIF
         IF ( med_diag%MPN%dgsave ) THEN
            CALL iom_put( "MPN"  ,fdpn2d )
            DEALLOCATE( fdpn2d )
         ENDIF
         IF ( med_diag%PRD%dgsave ) THEN
            CALL iom_put( "PRD"  ,fprd2d )
            DEALLOCATE( fprd2d )
         ENDIF
         IF( med_diag%MPD%dgsave ) THEN
            CALL iom_put( "MPD"  , fdpd2d )
            DEALLOCATE( fdpd2d )
         ENDIF
         !  IF( med_diag%DSED%dgsave ) THEN
         !      CALL iom_put( "DSED"  , ftot_n )
         !  ENDIF
         IF( med_diag%OPAL%dgsave ) THEN
            CALL iom_put( "OPAL"  , fprds2d )
            DEALLOCATE( fprds2d )
         ENDIF
         IF( med_diag%OPALDISS%dgsave ) THEN
            CALL iom_put( "OPALDISS"  , fsdiss2d )
            DEALLOCATE( fsdiss2d )
         ENDIF
         IF( med_diag%GMIPn%dgsave ) THEN
            CALL iom_put( "GMIPn"  , fgmipn2d )
            DEALLOCATE( fgmipn2d )
         ENDIF
         IF( med_diag%GMID%dgsave ) THEN
            CALL iom_put( "GMID"  , fgmid2d )
            DEALLOCATE( fgmid2d )
         ENDIF
         IF( med_diag%MZMI%dgsave ) THEN
            CALL iom_put( "MZMI"  , fdzmi2d )
            DEALLOCATE( fdzmi2d )
         ENDIF
         IF( med_diag%GMEPN%dgsave ) THEN
            CALL iom_put( "GMEPN"  , fgmepn2d )
            DEALLOCATE( fgmepn2d )
         ENDIF
         IF( med_diag%GMEPD%dgsave ) THEN
            CALL iom_put( "GMEPD"  , fgmepd2d )
            DEALLOCATE( fgmepd2d )
         ENDIF
         IF( med_diag%GMEZMI%dgsave ) THEN
            CALL iom_put( "GMEZMI"  , fgmezmi2d )
            DEALLOCATE( fgmezmi2d )
         ENDIF
         IF( med_diag%GMED%dgsave ) THEN
            CALL iom_put( "GMED"  , fgmed2d )
            DEALLOCATE( fgmed2d )
         ENDIF
         IF( med_diag%MZME%dgsave ) THEN
            CALL iom_put( "MZME"  , fdzme2d )
            DEALLOCATE( fdzme2d )
         ENDIF
         !  IF( med_diag%DEXP%dgsave ) THEN
         !      CALL iom_put( "DEXP"  , ftot_n )
         !  ENDIF
         IF( med_diag%DETN%dgsave ) THEN
            CALL iom_put( "DETN"  , fslown2d )
            DEALLOCATE( fslown2d )
         ENDIF
         IF( med_diag%MDET%dgsave ) THEN
            CALL iom_put( "MDET"  , fdd2d )
            DEALLOCATE( fdd2d )
         ENDIF
         IF( med_diag%AEOLIAN%dgsave ) THEN
            CALL iom_put( "AEOLIAN"  , ffetop2d )
            DEALLOCATE( ffetop2d )
         ENDIF
         IF( med_diag%BENTHIC%dgsave ) THEN
            CALL iom_put( "BENTHIC"  , ffebot2d )
            DEALLOCATE( ffebot2d )
         ENDIF
         IF( med_diag%SCAVENGE%dgsave ) THEN
            CALL iom_put( "SCAVENGE"  , ffescav2d )
            DEALLOCATE( ffescav2d )
         ENDIF
         !! 
         IF( med_diag%TOTREG_N%dgsave ) THEN
            CALL iom_put( "TOTREG_N"  , fregen2d )
            DEALLOCATE( fregen2d )
         ENDIF
         IF( med_diag%TOTRG_SI%dgsave ) THEN
            CALL iom_put( "TOTRG_SI"  , fregensi2d )
            DEALLOCATE( fregensi2d )
         ENDIF
         !! 
         IF( med_diag%FASTN%dgsave ) THEN
            CALL iom_put( "FASTN"  , ftempn2d )
            DEALLOCATE( ftempn2d )
         ENDIF
         IF( med_diag%FASTSI%dgsave ) THEN
            CALL iom_put( "FASTSI"  , ftempsi2d )
            DEALLOCATE( ftempsi2d )
         ENDIF
         IF( med_diag%FASTFE%dgsave ) THEN
            CALL iom_put( "FASTFE"  , ftempfe2d )
            DEALLOCATE( ftempfe2d )
         ENDIF
         IF( med_diag%FASTC%dgsave ) THEN
            CALL iom_put( "FASTC"  , ftempc2d )
            DEALLOCATE( ftempc2d )
         ENDIF
         IF( med_diag%FASTCA%dgsave ) THEN
            CALL iom_put( "FASTCA"  , ftempca2d )
            DEALLOCATE( ftempca2d )
         ENDIF
         !! 
         IF( med_diag%REMINN%dgsave ) THEN
            CALL iom_put( "REMINN"  , freminn2d )
            DEALLOCATE( freminn2d )
         ENDIF
         IF( med_diag%REMINSI%dgsave ) THEN
            CALL iom_put( "REMINSI"  , freminsi2d )
            DEALLOCATE( freminsi2d )
         ENDIF
         IF( med_diag%REMINFE%dgsave ) THEN
            CALL iom_put( "REMINFE"  , freminfe2d )
            DEALLOCATE( freminfe2d )
         ENDIF
         IF( med_diag%REMINC%dgsave ) THEN
            CALL iom_put( "REMINC"  , freminc2d )
            DEALLOCATE( freminc2d )
         ENDIF
         IF( med_diag%REMINCA%dgsave ) THEN
            CALL iom_put( "REMINCA"  , freminca2d )
            DEALLOCATE( freminca2d )
         ENDIF
         IF( med_diag%SEAFLRN%dgsave ) THEN
            CALL iom_put( "SEAFLRN"  , fsedn )
         ENDIF
         IF( med_diag%SEAFLRSI%dgsave ) THEN
            CALL iom_put( "SEAFLRSI"  , fsedsi )
         ENDIF
         IF( med_diag%SEAFLRFE%dgsave ) THEN
            CALL iom_put( "SEAFLRFE"  , fsedfe )
         ENDIF
         IF( med_diag%SEAFLRC%dgsave ) THEN
            CALL iom_put( "SEAFLRC"  , fsedc )
         ENDIF
         IF( med_diag%SEAFLRCA%dgsave ) THEN
            CALL iom_put( "SEAFLRCA"  , fsedca )
         ENDIF
         !!
# if defined key_roam            
         !!
         IF( med_diag%RIV_N%dgsave ) THEN
            CALL iom_put( "RIV_N"  , rivn2d )
            DEALLOCATE( rivn2d )
         ENDIF
         IF( med_diag%RIV_SI%dgsave ) THEN
            CALL iom_put( "RIV_SI"  , rivsi2d )
            DEALLOCATE( rivsi2d )
         ENDIF
         IF( med_diag%RIV_C%dgsave ) THEN
            CALL iom_put( "RIV_C"  , rivc2d )
            DEALLOCATE( rivc2d )
         ENDIF
         IF( med_diag%RIV_ALK%dgsave ) THEN
            CALL iom_put( "RIV_ALK"  , rivalk2d )
            DEALLOCATE( rivalk2d )
         ENDIF
         IF( med_diag%DETC%dgsave ) THEN
            CALL iom_put( "DETC"  , fslowc2d )
            DEALLOCATE( fslowc2d )
         ENDIF
         !!
         IF( med_diag%PN_LLOSS%dgsave ) THEN
            CALL iom_put( "PN_LLOSS"  , fdpn22d )
            DEALLOCATE( fdpn22d )
         ENDIF
         IF( med_diag%PD_LLOSS%dgsave ) THEN
            CALL iom_put( "PD_LLOSS"  , fdpd22d )
            DEALLOCATE( fdpd22d )
         ENDIF
         IF( med_diag%ZI_LLOSS%dgsave ) THEN
            CALL iom_put( "ZI_LLOSS"  , fdzmi22d )
             DEALLOCATE( fdzmi22d )
          ENDIF
          IF( med_diag%ZE_LLOSS%dgsave ) THEN
             CALL iom_put( "ZE_LLOSS"  , fdzme22d )
             DEALLOCATE( fdzme22d )
          ENDIF
          IF( med_diag%ZI_MES_N%dgsave ) THEN
             CALL iom_put( "ZI_MES_N"  , zimesn2d )
             DEALLOCATE( zimesn2d )
          ENDIF
          IF( med_diag%ZI_MES_D%dgsave ) THEN
             CALL iom_put( "ZI_MES_D"  , zimesd2d )
             DEALLOCATE( zimesd2d )
          ENDIF
          IF( med_diag%ZI_MES_C%dgsave ) THEN
             CALL iom_put( "ZI_MES_C"  , zimesc2d )
             DEALLOCATE( zimesc2d )
          ENDIF
          IF( med_diag%ZI_MESDC%dgsave ) THEN
             CALL iom_put( "ZI_MESDC"  ,zimesdc2d )
             DEALLOCATE( zimesdc2d )
          ENDIF
          IF( med_diag%ZI_EXCR%dgsave ) THEN
             CALL iom_put( "ZI_EXCR"  , ziexcr2d )
             DEALLOCATE( ziexcr2d )
          ENDIF
          IF( med_diag%ZI_RESP%dgsave ) THEN
             CALL iom_put( "ZI_RESP"  , ziresp2d )
             DEALLOCATE( ziresp2d )
          ENDIF
          IF( med_diag%ZI_GROW%dgsave ) THEN
             CALL iom_put( "ZI_GROW"  , zigrow2d )
             DEALLOCATE( zigrow2d )
          ENDIF
          IF( med_diag%ZE_MES_N%dgsave ) THEN
             CALL iom_put( "ZE_MES_N"  , zemesn2d )
             DEALLOCATE( zemesn2d )
          ENDIF
          IF( med_diag%ZE_MES_D%dgsave ) THEN
             CALL iom_put( "ZE_MES_D"  , zemesd2d )
             DEALLOCATE( zemesd2d )
          ENDIF
          IF( med_diag%ZE_MES_C%dgsave ) THEN
             CALL iom_put( "ZE_MES_C"  , zemesc2d )
             DEALLOCATE( zemesc2d )
          ENDIF
          IF( med_diag%ZE_MESDC%dgsave ) THEN
             CALL iom_put( "ZE_MESDC"  , zemesdc2d )
             DEALLOCATE( zemesdc2d )
          ENDIF
          IF( med_diag%ZE_EXCR%dgsave ) THEN
             CALL iom_put( "ZE_EXCR"  , zeexcr2d )
             DEALLOCATE( zeexcr2d )
          ENDIF
          IF( med_diag%ZE_RESP%dgsave ) THEN
             CALL iom_put( "ZE_RESP"  , zeresp2d )
             DEALLOCATE( zeresp2d )
          ENDIF
          IF( med_diag%ZE_GROW%dgsave ) THEN
             CALL iom_put( "ZE_GROW"  , zegrow2d )
             DEALLOCATE( zegrow2d )
          ENDIF
          IF( med_diag%MDETC%dgsave ) THEN
             CALL iom_put( "MDETC"  , mdetc2d )
             DEALLOCATE( mdetc2d )
          ENDIF
          IF( med_diag%GMIDC%dgsave ) THEN
             CALL iom_put( "GMIDC"  , gmidc2d )
             DEALLOCATE( gmidc2d )
          ENDIF
          IF( med_diag%GMEDC%dgsave ) THEN
             CALL iom_put( "GMEDC"  , gmedc2d )
             DEALLOCATE( gmedc2d )
          ENDIF
          IF( med_diag%IBEN_N%dgsave ) THEN
             CALL iom_put( "IBEN_N"  , iben_n2d )
             DEALLOCATE( iben_n2d )
          ENDIF
          IF( med_diag%IBEN_FE%dgsave ) THEN
             CALL iom_put( "IBEN_FE"  , iben_fe2d )
             DEALLOCATE( iben_fe2d )
          ENDIF
          IF( med_diag%IBEN_C%dgsave ) THEN
             CALL iom_put( "IBEN_C"  , iben_c2d )
             DEALLOCATE( iben_c2d )
          ENDIF
          IF( med_diag%IBEN_SI%dgsave ) THEN
             CALL iom_put( "IBEN_SI"  , iben_si2d )
             DEALLOCATE( iben_si2d )
          ENDIF
          IF( med_diag%IBEN_CA%dgsave ) THEN
             CALL iom_put( "IBEN_CA"  , iben_ca2d )
             DEALLOCATE( iben_ca2d )
          ENDIF
          IF( med_diag%OBEN_N%dgsave ) THEN
             CALL iom_put( "OBEN_N"  , oben_n2d )
             DEALLOCATE( oben_n2d )
          ENDIF
          IF( med_diag%OBEN_FE%dgsave ) THEN
             CALL iom_put( "OBEN_FE"  , oben_fe2d )
             DEALLOCATE( oben_fe2d )
          ENDIF
          IF( med_diag%OBEN_C%dgsave ) THEN
             CALL iom_put( "OBEN_C"  , oben_c2d )
             DEALLOCATE( oben_c2d )
          ENDIF
          IF( med_diag%OBEN_SI%dgsave ) THEN
             CALL iom_put( "OBEN_SI"  , oben_si2d )
             DEALLOCATE( oben_si2d )
          ENDIF
          IF( med_diag%OBEN_CA%dgsave ) THEN
             CALL iom_put( "OBEN_CA"  , oben_ca2d )
             DEALLOCATE( oben_ca2d )
          ENDIF
          IF( med_diag%SFR_OCAL%dgsave ) THEN
             CALL iom_put( "SFR_OCAL"  , sfr_ocal2d )
             DEALLOCATE( sfr_ocal2d )
          ENDIF
          IF( med_diag%SFR_OARG%dgsave ) THEN
             CALL iom_put( "SFR_OARG"  , sfr_oarg2d )
             DEALLOCATE( sfr_oarg2d )
          ENDIF
          IF( med_diag%LYSO_CA%dgsave ) THEN
             CALL iom_put( "LYSO_CA"  , lyso_ca2d )
             DEALLOCATE( lyso_ca2d )
          ENDIF
# endif                   
          !!
          !! ** 3D diagnostics
          IF( med_diag%TPP3%dgsave ) THEN
             CALL iom_put( "TPP3"  , tpp3d )
             DEALLOCATE( tpp3d )
          ENDIF
          IF( med_diag%DETFLUX3%dgsave ) THEN
             CALL iom_put( "DETFLUX3"  , detflux3d )
             DEALLOCATE( detflux3d )
          ENDIF
          IF( med_diag%REMIN3N%dgsave ) THEN
             CALL iom_put( "REMIN3N"  , remin3dn )
             DEALLOCATE( remin3dn )
          ENDIF
# if defined key_roam          
          IF( med_diag%PH3%dgsave ) THEN
             CALL iom_put( "PH3"  , f3_pH )
          ENDIF
          IF( med_diag%OM_CAL3%dgsave ) THEN
             CALL iom_put( "OM_CAL3"  , f3_omcal )
          ENDIF
          !!
          !! AXY (09/11/16): 2D CMIP6 diagnostics
          IF( med_diag%INTDISSIC%dgsave ) THEN
             CALL iom_put( "INTDISSIC"  , intdissic )
             DEALLOCATE( intdissic )
          ENDIF          
          IF( med_diag%INTDISSIN%dgsave ) THEN
             CALL iom_put( "INTDISSIN"  , intdissin )
             DEALLOCATE( intdissin )
          ENDIF          
          IF( med_diag%INTDISSISI%dgsave ) THEN
             CALL iom_put( "INTDISSISI"  , intdissisi )
             DEALLOCATE( intdissisi )
          ENDIF          
          IF( med_diag%INTTALK%dgsave ) THEN
             CALL iom_put( "INTTALK"  , inttalk )
             DEALLOCATE( inttalk )
          ENDIF          
          IF( med_diag%O2min%dgsave ) THEN
             CALL iom_put( "O2min"  , o2min )
             DEALLOCATE( o2min )
          ENDIF          
          IF( med_diag%ZO2min%dgsave ) THEN
             CALL iom_put( "ZO2min"  , zo2min )
             DEALLOCATE( zo2min )
          ENDIF          
          IF( med_diag%FBDDTALK%dgsave ) THEN
             CALL iom_put( "FBDDTALK"  , fbddtalk )
             DEALLOCATE( fbddtalk )
          ENDIF          
          IF( med_diag%FBDDTDIC%dgsave ) THEN
             CALL iom_put( "FBDDTDIC"  , fbddtdic )
             DEALLOCATE( fbddtdic )
          ENDIF          
          IF( med_diag%FBDDTDIFE%dgsave ) THEN
             CALL iom_put( "FBDDTDIFE" , fbddtdife )
             DEALLOCATE( fbddtdife )
          ENDIF          
          IF( med_diag%FBDDTDIN%dgsave ) THEN
             CALL iom_put( "FBDDTDIN"  , fbddtdin )
             DEALLOCATE( fbddtdin )
          ENDIF          
          IF( med_diag%FBDDTDISI%dgsave ) THEN
             CALL iom_put( "FBDDTDISI" , fbddtdisi )
             DEALLOCATE( fbddtdisi )
          ENDIF    
          !!
          !! AXY (09/11/16): 3D CMIP6 diagnostics
          IF( med_diag%TPPD3%dgsave ) THEN
             CALL iom_put( "TPPD3"     , tppd3 )
             DEALLOCATE( tppd3 )
          ENDIF          
          IF( med_diag%BDDTALK3%dgsave ) THEN
             CALL iom_put( "BDDTALK3"  , bddtalk3 )
             DEALLOCATE( bddtalk3 )
          ENDIF          
          IF( med_diag%BDDTDIC3%dgsave ) THEN
             CALL iom_put( "BDDTDIC3"  , bddtdic3 )
             DEALLOCATE( bddtdic3 )
          ENDIF          
          IF( med_diag%BDDTDIFE3%dgsave ) THEN
             CALL iom_put( "BDDTDIFE3" , bddtdife3 )
             DEALLOCATE( bddtdife3 )
          ENDIF          
          IF( med_diag%BDDTDIN3%dgsave ) THEN
             CALL iom_put( "BDDTDIN3"  , bddtdin3 )
             DEALLOCATE( bddtdin3 )
          ENDIF          
          IF( med_diag%BDDTDISI3%dgsave ) THEN
             CALL iom_put( "BDDTDISI3" , bddtdisi3 )
             DEALLOCATE( bddtdisi3 )
          ENDIF    
          IF( med_diag%FD_NIT3%dgsave ) THEN
             CALL iom_put( "FD_NIT3"  , fd_nit3 )
             DEALLOCATE( fd_nit3 )
          ENDIF
          IF( med_diag%FD_SIL3%dgsave ) THEN
             CALL iom_put( "FD_SIL3"  , fd_sil3 )
             DEALLOCATE( fd_sil3 )
          ENDIF
          IF( med_diag%FD_CAL3%dgsave ) THEN
             CALL iom_put( "FD_CAL3"  , fd_cal3 )
             DEALLOCATE( fd_cal3 )
          ENDIF
          IF( med_diag%FD_CAR3%dgsave ) THEN
             CALL iom_put( "FD_CAR3"  , fd_car3 )
             DEALLOCATE( fd_car3 )
          ENDIF
          IF( med_diag%CO33%dgsave ) THEN
             CALL iom_put( "CO33"  , f3_co3 )
          ENDIF          	        
          IF( med_diag%CO3SATARAG3%dgsave ) THEN
             CALL iom_put( "CO3SATARAG3"  , f3_omarg )
          ENDIF          	        
          IF( med_diag%CO3SATCALC3%dgsave ) THEN
             CALL iom_put( "CO3SATCALC3"  , f3_omcal )
          ENDIF          	        
          IF( med_diag%EXPC3%dgsave ) THEN
             CALL iom_put( "EXPC3"  , expc3 )
             DEALLOCATE( expc3 )
          ENDIF          	        
          IF( med_diag%EXPN3%dgsave ) THEN
             CALL iom_put( "EXPN3"  , expn3 )
             DEALLOCATE( expn3 )
          ENDIF          	        
          IF( med_diag%DCALC3%dgsave ) THEN
             CALL iom_put( "DCALC3"  , dcalc3 )
             DEALLOCATE( dcalc3 )
          ENDIF          	        
          IF( med_diag%FEDISS3%dgsave ) THEN
             CALL iom_put( "FEDISS3"  , fediss3 )
             DEALLOCATE( fediss3 )
          ENDIF          	        
          IF( med_diag%FESCAV3%dgsave ) THEN
             CALL iom_put( "FESCAV3"  , fescav3 )
             DEALLOCATE( fescav3 )
          ENDIF          	        
          IF( med_diag%MIGRAZP3%dgsave ) THEN
             CALL iom_put( "MIGRAZP3"  , migrazp3 )
             DEALLOCATE( migrazp3 )
          ENDIF          	        
          IF( med_diag%MIGRAZD3%dgsave ) THEN
             CALL iom_put( "MIGRAZD3"  , migrazd3 )
             DEALLOCATE( migrazd3 )
          ENDIF          	        
          IF( med_diag%MEGRAZP3%dgsave ) THEN
             CALL iom_put( "MEGRAZP3"  , megrazp3 )
             DEALLOCATE( megrazp3 )
          ENDIF          	        
          IF( med_diag%MEGRAZD3%dgsave ) THEN
             CALL iom_put( "MEGRAZD3"  , megrazd3 )
             DEALLOCATE( megrazd3 )
          ENDIF          	        
          IF( med_diag%MEGRAZZ3%dgsave ) THEN
             CALL iom_put( "MEGRAZZ3"  , megrazz3 )
             DEALLOCATE( megrazz3 )
          ENDIF          	        
          IF( med_diag%O2SAT3%dgsave ) THEN
             CALL iom_put( "O2SAT3"  , o2sat3 )
             DEALLOCATE( o2sat3 )
          ENDIF          	        
          IF( med_diag%PBSI3%dgsave ) THEN
             CALL iom_put( "PBSI3"  , pbsi3 )
             DEALLOCATE( pbsi3 )
          ENDIF          	        
          IF( med_diag%PCAL3%dgsave ) THEN
             CALL iom_put( "PCAL3"  , pcal3 )
             DEALLOCATE( pcal3 )
          ENDIF          	        
          IF( med_diag%REMOC3%dgsave ) THEN
             CALL iom_put( "REMOC3"  , remoc3 )
             DEALLOCATE( remoc3 )
          ENDIF          	        
          IF( med_diag%PNLIMJ3%dgsave ) THEN
             CALL iom_put( "PNLIMJ3" , pnlimj3 )
             DEALLOCATE( pnlimj3 )
          ENDIF          	        
          IF( med_diag%PNLIMN3%dgsave ) THEN
             CALL iom_put( "PNLIMN3" , pnlimn3 )
             DEALLOCATE( pnlimn3 )
          ENDIF          	        
          IF( med_diag%PNLIMFE3%dgsave ) THEN
             CALL iom_put( "PNLIMFE3" , pnlimfe3 )
             DEALLOCATE( pnlimfe3 )
          ENDIF          	        
          IF( med_diag%PDLIMJ3%dgsave ) THEN
             CALL iom_put( "PDLIMJ3" , pdlimj3 )
             DEALLOCATE( pdlimj3 )
          ENDIF          	        
          IF( med_diag%PDLIMN3%dgsave ) THEN
             CALL iom_put( "PDLIMN3" , pdlimn3 )
             DEALLOCATE( pdlimn3 )
          ENDIF          	        
          IF( med_diag%PDLIMFE3%dgsave ) THEN
             CALL iom_put( "PDLIMFE3" , pdlimfe3 )
             DEALLOCATE( pdlimfe3 )
          ENDIF          	        
          IF( med_diag%PDLIMSI3%dgsave ) THEN
             CALL iom_put( "PDLIMSI3" , pdlimsi3 )
             DEALLOCATE( pdlimsi3 )
          ENDIF          	        
          
# endif         
    !/CEB!$AGRIF_END_DO_NOT_TREAT
          DEALLOCATE( zw2d )

   END SUBROUTINE bio_medusa_fin

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE bio_medusa_fin( )                    ! Empty routine
      WRITE(*,*) 'bio_medusa_fin: You should not have seen this print! error?'
   END SUBROUTINE bio_medusa_fin
#endif 

   !!======================================================================
END MODULE bio_medusa_fin_mod
