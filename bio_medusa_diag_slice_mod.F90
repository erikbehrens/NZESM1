MODULE bio_medusa_diag_slice_mod
   !!======================================================================
   !!                         ***  MODULE bio_medusa_diag_slice_mod  ***
   !! Diagnostic calculations at different levels
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   bio_medusa_diag_slice     ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE bio_medusa_diag_slice( jk )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE bio_medusa_diag_slice  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - ...
      !!----------------------------------------------------------------------
      USE bio_medusa_mod
      USE par_medusa,        ONLY: jpchd, jpchn      
      USE dom_oce,           ONLY: tmask
      USE in_out_manager,    ONLY: lwp, numout
      USE iom,               ONLY: iom_put
      USE lbclnk,            ONLY: lbc_lnk
      USE oce,               ONLY: CO2Flux_out_cpl, DMS_out_cpl
      USE par_oce,           ONLY: jpi, jpj
      USE sbc_oce,           ONLY: lk_oasis, qsr, wndm
      USE sms_medusa,        ONLY: i0100, i0150, i0200, i0500, i1000,      &
                                   f2_ccd_arg, f2_ccd_cal,                 &
                                   f3_co3, f3_h2co3, f3_hco3, f3_pH,       &
                                   jdms, ocal_ccd, xpar, xze,              &
                                   zb_co2_flx, zb_dms_srf,                 &
                                   zn_co2_flx, zn_dms_srf
      !CEB USE trc,               ONLY: med_diag
      USE trc
      !/CEB
      !! The vertical level
      INTEGER, INTENT( in ) ::    jk
      !!----------------------------------------------------------------------

      !!-----------------------------------------
      !!
      !! 2d specific k level diags
      !!
      !!-----------------------------------------
#   if defined key_debug_medusa
         IF (lwp) write (numout,*) 'bio_medusa_diag_slice: start jk = ', jk
         CALL flush(numout)
#   endif
      !!
      IF (jk.eq.1) THEN
         IF( med_diag%MED_QSR%dgsave ) THEN
            CALL iom_put( "MED_QSR"  , qsr ) !
         ENDIF
         IF( med_diag%MED_XPAR%dgsave ) THEN
            CALL iom_put( "MED_XPAR"  , xpar(:,:,jk) ) !
         ENDIF       
         IF( med_diag%OCAL_CCD%dgsave ) THEN
            CALL iom_put( "OCAL_CCD"  , ocal_ccd ) !
         ENDIF
         IF( med_diag%FE_0000%dgsave ) THEN
            CALL iom_put( "FE_0000"  , xFree ) !
         ENDIF                     
         IF( med_diag%MED_XZE%dgsave ) THEN
            CALL iom_put( "MED_XZE"  , xze ) !
         ENDIF 
# if defined key_roam                     
         IF( med_diag%WIND%dgsave ) THEN
            CALL iom_put( "WIND"  , wndm )
         ENDIF
         IF( med_diag%ATM_PCO2%dgsave ) THEN
            CALL iom_put( "ATM_PCO2"  , f_pco2a2d )
            DEALLOCATE( f_pco2a2d )
         ENDIF
         IF( med_diag%OCN_PH%dgsave ) THEN
            zw2d(:,:) = f3_pH(:,:,jk)
            CALL iom_put( "OCN_PH"  , zw2d )
         ENDIF
         IF( med_diag%OCN_PCO2%dgsave ) THEN
            CALL iom_put( "OCN_PCO2"  , f_pco2w2d )
            DEALLOCATE( f_pco2w2d )
         ENDIF
         IF( med_diag%OCNH2CO3%dgsave ) THEN
            zw2d(:,:) = f3_h2co3(:,:,jk)
            CALL iom_put( "OCNH2CO3"  , zw2d )
         ENDIF
         IF( med_diag%OCN_HCO3%dgsave ) THEN
            zw2d(:,:) = f3_hco3(:,:,jk)
            CALL iom_put( "OCN_HCO3"  , zw2d )
         ENDIF
         IF( med_diag%OCN_CO3%dgsave ) THEN
            zw2d(:,:) = f3_co3(:,:,jk)
            CALL iom_put( "OCN_CO3"  , zw2d )
         ENDIF
         IF( med_diag%CO2FLUX%dgsave ) THEN
            CALL iom_put( "CO2FLUX"  , f_co2flux2d )
            DEALLOCATE( f_co2flux2d )
         ENDIF
         !! 
         !! AXY (10/11/16): repeat CO2 flux diagnostic in UKMO/CMIP6 units;
         !!                 this both outputs the CO2 flux in specified units
         !!                 and sends the resulting field to the coupler
         !! JPALM (17/11/16): put CO2 flux (fgco2) alloc/unalloc/pass to zn 
         !!                   out of diag list request 
         CALL lbc_lnk( fgco2(:,:),'T',1. )
         IF( med_diag%FGCO2%dgsave ) THEN
            CALL iom_put( "FGCO2"  , fgco2 )
         ENDIF
         !! JPALM (17/11/16): should mv this fgco2 part 
         !!                   out of lk_iomput loop
         zb_co2_flx = zn_co2_flx
         zn_co2_flx = fgco2
         IF (lk_oasis) THEN
            CO2Flux_out_cpl = zn_co2_flx
         ENDIF
         DEALLOCATE( fgco2 )
         !! ---
         IF( med_diag%OM_CAL%dgsave ) THEN
            CALL iom_put( "OM_CAL"  , f_omcal )
         ENDIF
         IF( med_diag%OM_ARG%dgsave ) THEN
            CALL iom_put( "OM_ARG"  , f_omarg )
         ENDIF
         IF( med_diag%TCO2%dgsave ) THEN
            CALL iom_put( "TCO2"  , f_TDIC2d )
            DEALLOCATE( f_TDIC2d )
         ENDIF
         IF( med_diag%TALK%dgsave ) THEN
            CALL iom_put( "TALK"  , f_TALK2d )
            DEALLOCATE( f_TALK2d )
         ENDIF
         IF( med_diag%KW660%dgsave ) THEN
            CALL iom_put( "KW660"  , f_kw6602d )
            DEALLOCATE( f_kw6602d )
         ENDIF
         IF( med_diag%ATM_PP0%dgsave ) THEN
            CALL iom_put( "ATM_PP0"  , f_pp02d )
            DEALLOCATE( f_pp02d )
         ENDIF
         IF( med_diag%O2FLUX%dgsave ) THEN
            CALL iom_put( "O2FLUX"  , f_o2flux2d )
            DEALLOCATE( f_o2flux2d )
         ENDIF
         IF( med_diag%O2SAT%dgsave ) THEN
            CALL iom_put( "O2SAT"  , f_o2sat2d )
            DEALLOCATE( f_o2sat2d )
         ENDIF
         IF( med_diag%CAL_CCD%dgsave ) THEN
            CALL iom_put( "CAL_CCD"  , f2_ccd_cal )
         ENDIF
         IF( med_diag%ARG_CCD%dgsave ) THEN
            CALL iom_put( "ARG_CCD"  , f2_ccd_arg )
         ENDIF
         IF (jdms .eq. 1) THEN
            IF( med_diag%DMS_SURF%dgsave ) THEN
               CALL lbc_lnk(dms_surf2d(:,:),'T',1. )
               CALL iom_put( "DMS_SURF"  , dms_surf2d )
               zb_dms_srf = zn_dms_srf
               zn_dms_srf = dms_surf2d
               IF (lk_oasis) THEN
                  DMS_out_cpl = zn_dms_srf
               ENDIF
               DEALLOCATE( dms_surf2d ) 
            ENDIF
            IF( med_diag%DMS_ANDR%dgsave ) THEN
               CALL iom_put( "DMS_ANDR"  , dms_andr2d )
               DEALLOCATE( dms_andr2d )
            ENDIF
            IF( med_diag%DMS_SIMO%dgsave ) THEN
               CALL iom_put( "DMS_SIMO"  , dms_simo2d )
               DEALLOCATE( dms_simo2d )
            ENDIF
            IF( med_diag%DMS_ARAN%dgsave ) THEN
               CALL iom_put( "DMS_ARAN"  , dms_aran2d )
               DEALLOCATE( dms_aran2d )
            ENDIF
            IF( med_diag%DMS_HALL%dgsave ) THEN
               CALL iom_put( "DMS_HALL"  , dms_hall2d )
               DEALLOCATE( dms_hall2d )
            ENDIF
            IF( med_diag%DMS_ANDM%dgsave ) THEN
               CALL iom_put( "DMS_ANDM"  , dms_andm2d )
               DEALLOCATE( dms_andm2d )
            ENDIF
         ENDIF
         !! AXY (24/11/16): extra MOCSY diagnostics
         IF( med_diag%ATM_XCO2%dgsave ) THEN
            CALL iom_put( "ATM_XCO2"  ,   f_xco2a_2d      )
            DEALLOCATE( f_xco2a_2d )
         ENDIF
         IF( med_diag%OCN_FCO2%dgsave ) THEN
            CALL iom_put( "OCN_FCO2"  ,   f_fco2w_2d      )
            DEALLOCATE( f_fco2w_2d )
         ENDIF
         IF( med_diag%ATM_FCO2%dgsave ) THEN
            CALL iom_put( "ATM_FCO2"  ,   f_fco2a_2d      )
            DEALLOCATE( f_fco2a_2d )
         ENDIF
         IF( med_diag%OCN_RHOSW%dgsave ) THEN
            CALL iom_put( "OCN_RHOSW"  ,  f_ocnrhosw_2d   )
            DEALLOCATE( f_ocnrhosw_2d )
         ENDIF
         IF( med_diag%OCN_SCHCO2%dgsave ) THEN
            CALL iom_put( "OCN_SCHCO2"  , f_ocnschco2_2d  )
            DEALLOCATE( f_ocnschco2_2d )
         ENDIF
         IF( med_diag%OCN_KWCO2%dgsave ) THEN
            CALL iom_put( "OCN_KWCO2"  ,  f_ocnkwco2_2d   )
            DEALLOCATE( f_ocnkwco2_2d )
         ENDIF
         IF( med_diag%OCN_K0%dgsave ) THEN
            CALL iom_put( "OCN_K0"  ,     f_ocnk0_2d      )
            DEALLOCATE( f_ocnk0_2d )
         ENDIF
         IF( med_diag%CO2STARAIR%dgsave ) THEN
            CALL iom_put( "CO2STARAIR"  , f_co2starair_2d )
            DEALLOCATE( f_co2starair_2d )
         ENDIF
         IF( med_diag%OCN_DPCO2%dgsave ) THEN
            CALL iom_put( "OCN_DPCO2"  ,  f_ocndpco2_2d   )
            DEALLOCATE( f_ocndpco2_2d )
         ENDIF
# endif                     
      ELSE IF (jk.eq.i0100) THEN 
         IF( med_diag%SDT__100%dgsave ) THEN
            zw2d(:,:) = fslownflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDT__100"  , zw2d )
         ENDIF
         IF( med_diag%REG__100%dgsave ) THEN
            CALL iom_put( "REG__100"  , fregen2d )
         ENDIF
         IF( med_diag%FDT__100%dgsave ) THEN
            CALL iom_put( "FDT__100"  , ffastn )
         ENDIF           
         IF( med_diag%RG__100F%dgsave ) THEN
            CALL iom_put( "RG__100F"  , fregenfast )
         ENDIF
         IF( med_diag%FDS__100%dgsave ) THEN
            CALL iom_put( "FDS__100"  , ffastsi )
         ENDIF         
         IF( med_diag%RGS_100F%dgsave ) THEN
            CALL iom_put( "RGS_100F"  , fregenfastsi )
         ENDIF
         IF( med_diag%FE_0100%dgsave ) THEN
            CALL iom_put( "FE_0100"  , xFree )
         ENDIF
# if defined key_roam                     
         IF( med_diag%RR_0100%dgsave ) THEN
            CALL iom_put( "RR_0100"  , ffastca2d )
         ENDIF                     
         IF( med_diag%SDC__100%dgsave ) THEN
            zw2d(:,:) = fslowcflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDC__100"  , zw2d )
         ENDIF                  
         IF( med_diag%epC100%dgsave    ) THEN
            zw2d(:,:) = (fslowcflux + ffastc) * tmask(:,:,jk)
            CALL iom_put( "epC100"    , zw2d )
         ENDIF		     
         IF( med_diag%epCALC100%dgsave ) THEN
            CALL iom_put( "epCALC100" , ffastca )
         ENDIF		     
         IF( med_diag%epN100%dgsave    ) THEN
            zw2d(:,:) = (fslownflux + ffastn) * tmask(:,:,jk)
            CALL iom_put( "epN100"    , zw2d )
         ENDIF		     
         IF( med_diag%epSI100%dgsave   ) THEN
            CALL iom_put( "epSI100"   , ffastsi )
         ENDIF		     
# endif                     
      ELSE IF (jk.eq.i0200) THEN
         IF( med_diag%SDT__200%dgsave ) THEN
            zw2d(:,:) = fslownflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDT__200"  , zw2d )
         ENDIF
         IF( med_diag%REG__200%dgsave ) THEN
            CALL iom_put( "REG__200"  , fregen2d )
         ENDIF
         IF( med_diag%FDT__200%dgsave ) THEN
            CALL iom_put( "FDT__200"  , ffastn )
         ENDIF
         IF( med_diag%RG__200F%dgsave ) THEN
            CALL iom_put( "RG__200F"  , fregenfast )
         ENDIF
         IF( med_diag%FDS__200%dgsave ) THEN
            CALL iom_put( "FDS__200"  , ffastsi )
         ENDIF
         IF( med_diag%RGS_200F%dgsave ) THEN
            CALL iom_put( "RGS_200F"  , fregenfastsi )
         ENDIF
         IF( med_diag%FE_0200%dgsave ) THEN
            CALL iom_put( "FE_0200"   , xFree )
         ENDIF
# if defined key_roam                     
         IF( med_diag%SDC__200%dgsave ) THEN
            zw2d(:,:) = fslowcflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDC__200"  , zw2d )
         ENDIF
# endif                     
      ELSE IF (jk.eq.i0500) THEN
         IF( med_diag%SDT__500%dgsave ) THEN
            zw2d(:,:) = fslownflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDT__500"  , zw2d )
         ENDIF
         IF( med_diag%REG__500%dgsave ) THEN
            CALL iom_put( "REG__500"  , fregen2d )
         ENDIF      
         IF( med_diag%FDT__500%dgsave ) THEN
            CALL iom_put( "FDT__500"  , ffastn )
         ENDIF
         IF( med_diag%RG__500F%dgsave ) THEN
            CALL iom_put( "RG__500F"  , fregenfast )
         ENDIF
         IF( med_diag%FDS__500%dgsave ) THEN
            CALL iom_put( "FDS__500"  , ffastsi )
         ENDIF
         IF( med_diag%RGS_500F%dgsave ) THEN
            CALL iom_put( "RGS_500F"  , fregenfastsi )
         ENDIF
         IF( med_diag%FE_0500%dgsave ) THEN
            CALL iom_put( "FE_0500"  , xFree )
         ENDIF
# if defined key_roam                     
         IF( med_diag%RR_0500%dgsave ) THEN
            CALL iom_put( "RR_0500"  , ffastca2d )
         ENDIF
         IF( med_diag%SDC__500%dgsave ) THEN
            zw2d(:,:) = fslowcflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDC__500"  , zw2d )
         ENDIF  
# endif                      
      ELSE IF (jk.eq.i1000) THEN
         IF( med_diag%SDT_1000%dgsave ) THEN
            zw2d(:,:) = fslownflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDT_1000"  , zw2d )
         ENDIF
         IF( med_diag%REG_1000%dgsave ) THEN
            CALL iom_put( "REG_1000"  , fregen2d )
         ENDIF  
         IF( med_diag%FDT_1000%dgsave ) THEN
            CALL iom_put( "FDT_1000"  , ffastn )
         ENDIF
         IF( med_diag%RG_1000F%dgsave ) THEN
            CALL iom_put( "RG_1000F"  , fregenfast )
         ENDIF
         IF( med_diag%FDS_1000%dgsave ) THEN
            CALL iom_put( "FDS_1000"  , ffastsi )
         ENDIF
         IF( med_diag%RGS1000F%dgsave ) THEN
            CALL iom_put( "RGS1000F"  , fregenfastsi )
         ENDIF
         IF( med_diag%FE_1000%dgsave ) THEN
            CALL iom_put( "FE_1000"  , xFree )
         ENDIF
# if defined key_roam                     
         IF( med_diag%RR_1000%dgsave ) THEN
            CALL iom_put( "RR_1000"  , ffastca2d )
            DEALLOCATE( ffastca2d )
         ENDIF
         IF( med_diag%SDC_1000%dgsave ) THEN
            zw2d(:,:) = fslowcflux(:,:) * tmask(:,:,jk)
            CALL iom_put( "SDC_1000"  , zw2d )
         ENDIF 
# endif                      
      ENDIF
      !! to do on every k loop :
      IF( med_diag%DETFLUX3%dgsave ) THEN
         !! detrital flux
         detflux3d(:,:,jk) = (fslownflux(:,:) + ffastn(:,:)) * tmask(:,:,jk)
         !CALL iom_put( "DETFLUX3"  , ftot_n )
      ENDIF
# if defined key_roam                     
      IF( med_diag%EXPC3%dgsave ) THEN
         expc3(:,:,jk) = (fslowcflux(:,:) + ffastc(:,:)) * tmask(:,:,jk)
      ENDIF		     
      IF( med_diag%EXPN3%dgsave ) THEN
         expn3(:,:,jk) = (fslownflux(:,:) + ffastn(:,:)) * tmask(:,:,jk)
      ENDIF		     
# endif		     
#   if defined key_debug_medusa
         IF (lwp) write (numout,*) 'bio_medusa_diag_slice: end jk = ', jk
         CALL flush(numout)
#   endif

   END SUBROUTINE bio_medusa_diag_slice

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE bio_medusa_diag_slice( )                  ! Empty routine
      WRITE(*,*) 'bio_medusa_diag_slice: You should not have seen this print! error?'
   END SUBROUTINE bio_medusa_diag_slice
#endif 

   !!======================================================================
END MODULE bio_medusa_diag_slice_mod
