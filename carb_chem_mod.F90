MODULE carb_chem_mod
   !!======================================================================
   !!                         ***  MODULE carb_chem_mod  ***
   !! Calculate the carbon chemistry for the whole ocean
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!----------------------------------------------------------------------
#if defined key_roam
   !!----------------------------------------------------------------------
   !!                                                   key_roam?
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   carb_chem        ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE carb_chem( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE carb_chem  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - ...
      !!----------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: f_co2flux, f_co3, f_dcf,               &
                                   f_h2co3, f_hco3, f_henry,              &
                                   f_kw660, f_omarg, f_omcal,             &
                                   f_pco2atm, f_pco2w, f_ph, f_pp0,       &
                                   f_TALK, f_TDIC, f_xco2a,               &
# if defined key_mocsy
                                   zpho,                                  &
# endif
                                   zalk, zdic, zsal, zsil, ztmp 
      USE dom_oce,           ONLY: gdept_0, gdepw_0, gphit, mbathy, tmask
# if defined key_vvl
      USE dom_oce,           ONLY: gdept_n, gdepw_n
# endif
      USE in_out_manager,    ONLY: lwp, numout
      USE oce,               ONLY: PCO2a_in_cpl, tsb, tsn
      USE par_kind,          ONLY: wp
      USE par_medusa,        ONLY: jpalk, jpdic, jpdin, jpsil 
      USE par_oce,           ONLY: jp_sal, jp_tem, jpi, jpim1, jpj,       &
                                   jpjm1, jpk
      USE sbc_oce,           ONLY: lk_oasis
      USE sms_medusa,        ONLY: f2_ccd_arg, f2_ccd_cal, f3_co3,        &
                                   f3_h2co3, f3_hco3, f3_omarg,           &
                                   f3_omcal, f3_pH
      USE trc,               ONLY: trn

# if defined key_mocsy
      USE mocsy_wrapper,     ONLY: mocsy_interface
# else
      USE trcco2_medusa,     ONLY: trc_co2_medusa
# endif

   !!* Substitution
#  include "domzgr_substitute.h90"
 
      !! time (integer timestep)
      INTEGER, INTENT( in ) ::    kt

      !! Flags to help with calculating the position of the CCD
      INTEGER, DIMENSION(jpi,jpj) ::     i2_omcal,i2_omarg

      !! AXY (23/06/15): additional diagnostics for MOCSY and oxygen
      REAL(wp) :: f_rhosw
      !! Output arguments from mocsy_interface, which aren't used
      REAL(wp) :: f_fco2w_dum, f_BetaD_dum, f_opres_dum
      REAL(wp) :: f_insitut_dum, f_fco2atm_dum
      REAL(wp) :: f_schmidtco2_dum, f_kwco2_dum, f_K0_dum
      REAL(wp) :: f_co2starair_dum, f_dpco2_dum
      !! temporary variables
      REAL(wp) ::    fq0,fq1,fq2,fq3,fq4

      INTEGER :: iters
      !! Loop variables
      INTEGER :: ji, jj, jk

      !!----------------------------------------------------------------------
      !! Calculate the carbonate chemistry for the whole ocean on the first
      !! simulation timestep and every month subsequently; the resulting 3D
      !! field of omega calcite is used to determine the depth of the CCD
      !!----------------------------------------------------------------------
      !!
!      IF(lwp) WRITE(numout,*)                                               &
!              ' MEDUSA calculating all carbonate chemistry at kt =', kt
!         CALL flush(numout)
      !! blank flags
      i2_omcal(:,:) = 0
      i2_omarg(:,:) = 0

      !! Loop over levels
      DO jk = 1,jpk

         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               !! OPEN wet point IF..THEN loop
               IF (tmask(ji,jj,jk).eq.1) THEN
                  !! Do carbonate chemistry
                  !!
                  !! Set up required state variables.
                  !! dissolved inorganic carbon
                  zdic(ji,jj) = max(0.,trn(ji,jj,jk,jpdic))
                  !! alkalinity
                  zalk(ji,jj) = max(0.,trn(ji,jj,jk,jpalk))
                  !! temperature
                  ztmp(ji,jj) = tsn(ji,jj,jk,jp_tem)
                  !! salinity
                  zsal(ji,jj) = tsn(ji,jj,jk,jp_sal)
#  if defined key_mocsy
                  !! silicic acid
                  zsil(ji,jj) = max(0.,trn(ji,jj,jk,jpsil))
                  !! phosphate via DIN and Redfield
                  zpho(ji,jj) = max(0.,trn(ji,jj,jk,jpdin)) / 16.0 
#  endif
		  !!
		  !! AXY (28/02/14): check input fields
		  if (ztmp(ji,jj) .lt. -3.0 .or. ztmp(ji,jj) .gt. 40.0 ) then
                     IF (lwp) WRITE(numout,*) ' carb_chem: T WARNING 3D, ',  &
                        tsb(ji,jj,jk,jp_tem), tsn(ji,jj,jk,jp_tem),          &
                        ' at (', ji, ',', jj, ',', jk, ') at time', kt
	             IF(lwp) WRITE(numout,*)                                 &
                        ' carb_chem: T SWITCHING 3D, ',                      &
			tsn(ji,jj,jk,jp_tem), ' -> ', tsb(ji,jj,jk,jp_tem)
                     ztmp(ji,jj) = tsb(ji,jj,jk,jp_tem)     !! temperature
                  endif
		  if (zsal(ji,jj) .lt. 0.0 .or. zsal(ji,jj) .gt. 45.0 ) then
                     IF (lwp) WRITE(numout,*) ' carb_chem: S WARNING 3D, ',  &
                        tsb(ji,jj,jk,jp_sal), tsn(ji,jj,jk,jp_sal),          &
                        ' at (', ji, ',', jj, ',', jk, ') at time', kt
                  endif
               ENDIF
            ENDDO
         ENDDO

         !! Blank input variables not used at this stage (they 
         !! relate to air-sea flux)
         f_kw660(:,:) = 1.0
         f_pp0(:,:)   = 1.0

         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk).eq.1) THEN
                  !! calculate carbonate chemistry at grid cell midpoint
#  if defined key_mocsy
                  !! AXY (22/06/15): use Orr & Epitalon (2015) MOCSY-2 carbonate
                  !!                 chemistry package
                  CALL mocsy_interface(ztmp(ji,jj),zsal(ji,jj),zalk(ji,jj), &
                                       zdic(ji,jj),zsil(ji,jj),zpho(ji,jj), &
                                       f_pp0(ji,jj),fsdept(ji,jj,jk),       &
                                       gphit(ji,jj),f_kw660(ji,jj),         &
                                       f_xco2a(ji,jj),1,f_ph(ji,jj),        &
                                       f_pco2w(ji,jj),f_fco2w_dum,          &
                                       f_h2co3(ji,jj),f_hco3(ji,jj),        &
                                       f_co3(ji,jj),f_omarg(ji,jj),         &
                                       f_omcal(ji,jj),f_BetaD_dum,          &
                                       f_rhosw,f_opres_dum,                 &
                                       f_insitut_dum,f_pco2atm(ji,jj),      &
                                       f_fco2atm_dum,f_schmidtco2_dum,      &
                                       f_kwco2_dum,f_K0_dum,                &
                                       f_co2starair_dum,f_co2flux(ji,jj),   & 
                                       f_dpco2_dum)
                  !!
                  !! mmol / m3 -> umol / kg
                  f_TDIC(ji,jj) = (zdic(ji,jj) / f_rhosw) * 1000.
                  !! meq / m3 -> ueq / kg
                  f_TALK(ji,jj) = (zalk(ji,jj) / f_rhosw) * 1000.
                  f_dcf(ji,jj)  = f_rhosw
#  else
                  !! AXY (22/06/15): use old PML carbonate chemistry 
                  !! package (the MEDUSA-2 default)
                  CALL trc_co2_medusa(ztmp(ji,jj),zsal(ji,jj),zdic(ji,jj),  &
                                      zalk(ji,jj),fsdept(ji,jj,jk),         &
                                      f_kw660(ji,jj),f_xco2a(ji,jj),        &
                                      f_ph(ji,jj),                          &
                                      f_pco2w(ji,jj),f_h2co3(ji,jj),        &
                                      f_hco3(ji,jj),f_co3(ji,jj),           &
                                      f_omcal(ji,jj),f_omarg(ji,jj),        &
                                      f_co2flux(ji,jj),f_TDIC(ji,jj),       &
                                      f_TALK(ji,jj),f_dcf(ji,jj),           &
                                      f_henry(ji,jj),iters)
                  !! 
                  !! AXY (28/02/14): check output fields
                  IF (iters .eq. 25) THEN
                     IF(lwp) WRITE(numout,*)                                &
                        ' carb_chem: 3D ITERS WARNING, ', iters, ' AT (',   &
                        ji, ', ', jj, ', ', jk, ') AT ', kt
                  ENDIF
#  endif
                  !!
                  !! store 3D outputs
                  f3_pH(ji,jj,jk)    = f_ph(ji,jj)
                  f3_h2co3(ji,jj,jk) = f_h2co3(ji,jj)
                  f3_hco3(ji,jj,jk)  = f_hco3(ji,jj)
                  f3_co3(ji,jj,jk)   = f_co3(ji,jj)
                  f3_omcal(ji,jj,jk) = f_omcal(ji,jj)
                  f3_omarg(ji,jj,jk) = f_omarg(ji,jj)
                  !!
                  !! CCD calculation: calcite
                  if (i2_omcal(ji,jj) == 0 .and. f_omcal(ji,jj) < 1.0) then
                     if (jk .eq. 1) then
                        f2_ccd_cal(ji,jj) = fsdept(ji,jj,jk)
                     else
                        fq0 = f3_omcal(ji,jj,jk-1) - f_omcal(ji,jj)
                        fq1 = f3_omcal(ji,jj,jk-1) - 1.0
                        fq2 = fq1 / (fq0 + tiny(fq0))
                        fq3 = fsdept(ji,jj,jk) - fsdept(ji,jj,jk-1)
                        fq4 = fq2 * fq3
                        f2_ccd_cal(ji,jj) = fsdept(ji,jj,jk-1) + fq4
                     endif
                     i2_omcal(ji,jj)   = 1
                  endif
                  if ( i2_omcal(ji,jj) == 0 .and. jk == mbathy(ji,jj) ) then
                     !! reached seafloor and still no dissolution; set 
                     !! to seafloor (W-point)
                     f2_ccd_cal(ji,jj) = fsdepw(ji,jj,jk+1)
                     i2_omcal(ji,jj)   = 1
                  endif
                  !!
                  !! CCD calculation: aragonite
                  if (i2_omarg(ji,jj) == 0 .and. f_omarg(ji,jj) < 1.0) then
                     if (jk .eq. 1) then
                        f2_ccd_arg(ji,jj) = fsdept(ji,jj,jk)
                     else
                        fq0 = f3_omarg(ji,jj,jk-1) - f_omarg(ji,jj)
                        fq1 = f3_omarg(ji,jj,jk-1) - 1.0
                        fq2 = fq1 / (fq0 + tiny(fq0))
                        fq3 = fsdept(ji,jj,jk) - fsdept(ji,jj,jk-1)
                        fq4 = fq2 * fq3
                        f2_ccd_arg(ji,jj) = fsdept(ji,jj,jk-1) + fq4
                     endif
                     i2_omarg(ji,jj)   = 1
                  endif
                  if ( i2_omarg(ji,jj) == 0 .and. jk == mbathy(ji,jj) ) then
                     !! reached seafloor and still no dissolution; set 
                     !! to seafloor (W-point)
                     f2_ccd_arg(ji,jj) = fsdepw(ji,jj,jk+1)
                     i2_omarg(ji,jj)   = 1
                  endif
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE carb_chem

#else
   !!======================================================================
   !!  Dummy module :                                   None key_roam?
   !!======================================================================
CONTAINS
   SUBROUTINE carb_chem( )                    ! Empty routine
      WRITE(*,*) 'carb_chem: You should not have seen this print! error?'
   END SUBROUTINE carb_chem
#endif 

   !!======================================================================
END MODULE carb_chem_mod
