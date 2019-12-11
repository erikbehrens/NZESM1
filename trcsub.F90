MODULE trcsub
   !!======================================================================
   !!                       ***  MODULE trcsubstp  ***
   !!TOP :   Averages physics variables for TOP substepping. 
   !!======================================================================
   !! History :  1.0  !  2011-10  (K. Edwards)  Original
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   trc_sub    : passive tracer system sub-stepping 
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE trc
   USE prtctl_trc       ! Print control for debbuging
   USE iom
   USE in_out_manager
   USE lbclnk
#if defined key_zdftke
   USE zdftke          ! twice TKE (en)
#endif
#if defined key_zdfgls
   USE zdfgls, ONLY: en
#endif
   USE trabbl
   USE zdf_oce
   USE domvvl
   USE divcur          ! hor. divergence and curl      (div & cur routines)
   USE sbcrnf, ONLY: h_rnf, nk_rnf   ! River runoff 
   USE bdy_oce
#if defined key_agrif
   USE agrif_opa_update
   USE agrif_opa_interp
#endif

   IMPLICIT NONE

   PUBLIC   trc_sub_stp      ! called by trc_stp
   PUBLIC   trc_sub_ini      ! called by trc_ini to initialize substepping arrays.
   PUBLIC   trc_sub_reset    ! called by trc_stp to reset physics variables
   PUBLIC   trc_sub_ssh      ! called by trc_stp to reset physics variables

   !!* Module variables
   REAL(wp)  :: r1_ndttrc     !    1 /  nn_dttrc 
   REAL(wp)  :: r1_ndttrcp1   !    1 / (nn_dttrc+1) 

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sub_stp( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Average variables needed for sub-stepping passive tracers
      !! 
      !! ** Method  : Called every timestep to increment _tm (time mean) variables
      !!              on TOP steps, calculate averages.
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt        ! ocean time-step index
      INTEGER               ::  ji,jj,jk  ! dummy loop indices
      REAL(wp)              ::  z1_ne3t, z1_ne3u, z1_ne3v, z1_ne3w
      !!-------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sub_stp')
      !
      IF( kt == nit000 ) THEN
           IF(lwp) WRITE(numout,*)
           IF(lwp) WRITE(numout,*) 'trc_sub_stp : substepping of the passive tracers'
           IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
           !
           sshb_hold  (:,:) = sshn  (:,:)
           emp_b_hold (:,:) = emp_b (:,:)
           !
           r1_ndttrc        = 1._wp / REAL( nn_dttrc    , wp ) 
           r1_ndttrcp1      = 1._wp / REAL( nn_dttrc + 1, wp )
           !
      ENDIF  

       IF( MOD( kt , nn_dttrc ) /= 0 ) THEN
          !
          un_tm   (:,:,:)        = un_tm   (:,:,:)        + un   (:,:,:)        * fse3u(:,:,:) 
          vn_tm   (:,:,:)        = vn_tm   (:,:,:)        + vn   (:,:,:)        * fse3v(:,:,:) 
          tsn_tm  (:,:,:,jp_tem) = tsn_tm  (:,:,:,jp_tem) + tsn  (:,:,:,jp_tem) * fse3t(:,:,:)  
          tsn_tm  (:,:,:,jp_sal) = tsn_tm  (:,:,:,jp_sal) + tsn  (:,:,:,jp_sal) * fse3t(:,:,:)  
          rhop_tm (:,:,:)        = rhop_tm (:,:,:)        + rhop (:,:,:)        * fse3t(:,:,:)  
          avt_tm  (:,:,:)        = avt_tm  (:,:,:)        + avt  (:,:,:)        * fse3w(:,:,:)  
# if defined key_zdfddm
          avs_tm  (:,:,:)        = avs_tm  (:,:,:)        + avs  (:,:,:)        * fse3w(:,:,:)  
# endif
#if defined key_ldfslp
          wslpi_tm(:,:,:)        = wslpi_tm(:,:,:)        + wslpi(:,:,:)
          wslpj_tm(:,:,:)        = wslpj_tm(:,:,:)        + wslpj(:,:,:)
          uslp_tm (:,:,:)        = uslp_tm (:,:,:)        + uslp (:,:,:)
          vslp_tm (:,:,:)        = vslp_tm (:,:,:)        + vslp (:,:,:)
#endif
# if defined key_trabbl
          IF( nn_bbl_ldf == 1 ) THEN
             ahu_bbl_tm(:,:)     = ahu_bbl_tm(:,:)        + ahu_bbl(:,:) 
             ahv_bbl_tm(:,:)     = ahv_bbl_tm(:,:)        + ahv_bbl(:,:) 
          ENDIF
          IF( nn_bbl_adv == 1 ) THEN
             utr_bbl_tm(:,:)     = utr_bbl_tm(:,:)        + utr_bbl(:,:) 
             vtr_bbl_tm(:,:)     = vtr_bbl_tm(:,:)        + vtr_bbl(:,:) 
          ENDIF
# endif
          !
          sshn_tm  (:,:)         = sshn_tm  (:,:)         + sshn  (:,:) 
          rnf_tm   (:,:)         = rnf_tm   (:,:)         + rnf   (:,:) 
          h_rnf_tm (:,:)         = h_rnf_tm (:,:)         + h_rnf (:,:) 
          hmld_tm  (:,:)         = hmld_tm  (:,:)         + hmld  (:,:)
          fr_i_tm  (:,:)         = fr_i_tm  (:,:)         + fr_i  (:,:)
          emp_tm   (:,:)         = emp_tm   (:,:)         + emp   (:,:) 
          fmmflx_tm(:,:)         = fmmflx_tm(:,:)         + fmmflx(:,:)
          qsr_tm   (:,:)         = qsr_tm   (:,:)         + qsr   (:,:)
          wndm_tm  (:,:)         = wndm_tm  (:,:)         + wndm  (:,:)

      ELSE                           !  It is time to substep 
         !   1. set temporary arrays to hold physics variables
         un_temp    (:,:,:)      = un    (:,:,:)
         vn_temp    (:,:,:)      = vn    (:,:,:)
         wn_temp    (:,:,:)      = wn    (:,:,:)
         tsn_temp   (:,:,:,:)    = tsn   (:,:,:,:)
         rhop_temp  (:,:,:)      = rhop  (:,:,:)    
         avt_temp   (:,:,:)      = avt   (:,:,:)
# if defined key_zdfddm
         avs_temp   (:,:,:)      = avs   (:,:,:)
# endif
#if defined key_ldfslp
         wslpi_temp (:,:,:)      = wslpi (:,:,:)
         wslpj_temp (:,:,:)      = wslpj (:,:,:)
         uslp_temp  (:,:,:)      = uslp  (:,:,:)
         vslp_temp  (:,:,:)      = vslp  (:,:,:)
#endif
# if defined key_trabbl
          IF( nn_bbl_ldf == 1 ) THEN
             ahu_bbl_temp(:,:)   = ahu_bbl(:,:)  
             ahv_bbl_temp(:,:)   = ahv_bbl(:,:) 
          ENDIF
          IF( nn_bbl_adv == 1 ) THEN
             utr_bbl_temp(:,:)   = utr_bbl(:,:) 
             vtr_bbl_temp(:,:)   = vtr_bbl(:,:) 
          ENDIF
# endif
         sshn_temp  (:,:)        = sshn  (:,:)
         sshb_temp  (:,:)        = sshb  (:,:)
         ssha_temp  (:,:)        = ssha  (:,:)
         rnf_temp   (:,:)        = rnf   (:,:)
         h_rnf_temp (:,:)        = h_rnf (:,:)
         hmld_temp  (:,:)        = hmld  (:,:)
         fr_i_temp  (:,:)        = fr_i  (:,:)
         emp_temp   (:,:)        = emp   (:,:)
         emp_b_temp (:,:)        = emp_b (:,:)
         fmmflx_temp(:,:)        = fmmflx(:,:)
         qsr_temp   (:,:)        = qsr   (:,:)
         wndm_temp  (:,:)        = wndm  (:,:)
         !                                    !  Variables reset in trc_sub_ssh
         rotn_temp  (:,:,:)      = rotn  (:,:,:)
         hdivn_temp (:,:,:)      = hdivn (:,:,:)
         rotb_temp  (:,:,:)      = rotb  (:,:,:)
         hdivb_temp (:,:,:)      = hdivb (:,:,:)
         !
         ! 2. Create averages and reassign variables
         un_tm    (:,:,:)        = un_tm   (:,:,:)        + un   (:,:,:)        * fse3u(:,:,:) 
         vn_tm    (:,:,:)        = vn_tm   (:,:,:)        + vn   (:,:,:)        * fse3v(:,:,:) 
         tsn_tm   (:,:,:,jp_tem) = tsn_tm  (:,:,:,jp_tem) + tsn  (:,:,:,jp_tem) * fse3t(:,:,:)  
         tsn_tm   (:,:,:,jp_sal) = tsn_tm  (:,:,:,jp_sal) + tsn  (:,:,:,jp_sal) * fse3t(:,:,:)  
         rhop_tm (:,:,:)         = rhop_tm (:,:,:)        + rhop (:,:,:)        * fse3t(:,:,:)  
         avt_tm   (:,:,:)        = avt_tm  (:,:,:)        + avt  (:,:,:)        * fse3w(:,:,:)  
# if defined key_zdfddm
         avs_tm   (:,:,:)        = avs_tm  (:,:,:)        + avs  (:,:,:)        * fse3w(:,:,:)  
# endif
#if defined key_ldfslp
         wslpi_tm (:,:,:)        = wslpi_tm(:,:,:)        + wslpi(:,:,:) 
         wslpj_tm (:,:,:)        = wslpj_tm(:,:,:)        + wslpj(:,:,:) 
         uslp_tm  (:,:,:)        = uslp_tm (:,:,:)        + uslp (:,:,:) 
         vslp_tm  (:,:,:)        = vslp_tm (:,:,:)        + vslp (:,:,:)
#endif
# if defined key_trabbl
          IF( nn_bbl_ldf == 1 ) THEN
             ahu_bbl_tm(:,:)     = ahu_bbl_tm(:,:)        + ahu_bbl(:,:) 
             ahv_bbl_tm(:,:)     = ahv_bbl_tm(:,:)        + ahv_bbl(:,:) 
          ENDIF
          IF( nn_bbl_adv == 1 ) THEN
             utr_bbl_tm(:,:)     = utr_bbl_tm(:,:)        + utr_bbl(:,:) 
             vtr_bbl_tm(:,:)     = vtr_bbl_tm(:,:)        + vtr_bbl(:,:) 
          ENDIF
# endif
         sshn_tm  (:,:)          = sshn_tm    (:,:)       + sshn  (:,:) 
         rnf_tm   (:,:)          = rnf_tm     (:,:)       + rnf   (:,:) 
         h_rnf_tm (:,:)          = h_rnf_tm   (:,:)       + h_rnf (:,:) 
         hmld_tm  (:,:)          = hmld_tm    (:,:)       + hmld  (:,:)
         fr_i_tm  (:,:)          = fr_i_tm    (:,:)       + fr_i  (:,:)
         emp_tm   (:,:)          = emp_tm     (:,:)       + emp   (:,:) 
         fmmflx_tm(:,:)          = fmmflx_tm  (:,:)       + fmmflx(:,:)
         qsr_tm   (:,:)          = qsr_tm     (:,:)       + qsr   (:,:)
         wndm_tm  (:,:)          = wndm_tm    (:,:)       + wndm  (:,:)
         !
         sshn     (:,:)          = sshn_tm    (:,:) * r1_ndttrcp1 
         sshb     (:,:)          = sshb_hold  (:,:)
         rnf      (:,:)          = rnf_tm     (:,:) * r1_ndttrcp1 
         h_rnf    (:,:)          = h_rnf_tm   (:,:) * r1_ndttrcp1 
         hmld     (:,:)          = hmld_tm    (:,:) * r1_ndttrcp1 
         !  variables that are initialized after averages
         emp_b    (:,:) = emp_b_hold (:,:)
         IF( kt == nittrc000 ) THEN
            wndm  (:,:)          = wndm_tm    (:,:) * r1_ndttrc 
            qsr   (:,:)          = qsr_tm     (:,:) * r1_ndttrc 
            emp   (:,:)          = emp_tm     (:,:) * r1_ndttrc 
            fmmflx(:,:)          = fmmflx_tm  (:,:) * r1_ndttrc 
            fr_i  (:,:)          = fr_i_tm    (:,:) * r1_ndttrc
# if defined key_trabbl
            IF( nn_bbl_ldf == 1 ) THEN
               ahu_bbl(:,:)      = ahu_bbl_tm (:,:) * r1_ndttrc  
               ahv_bbl(:,:)      = ahv_bbl_tm (:,:) * r1_ndttrc 
            ENDIF
            IF( nn_bbl_adv == 1 ) THEN
               utr_bbl(:,:)      = utr_bbl_tm (:,:) * r1_ndttrc  
               vtr_bbl(:,:)      = vtr_bbl_tm (:,:) * r1_ndttrc 
            ENDIF
# endif
         ELSE
            wndm  (:,:)          = wndm_tm    (:,:) * r1_ndttrcp1 
            qsr   (:,:)          = qsr_tm     (:,:) * r1_ndttrcp1 
            emp   (:,:)          = emp_tm     (:,:) * r1_ndttrcp1 
            fmmflx(:,:)          = fmmflx_tm  (:,:) * r1_ndttrcp1 
            fr_i  (:,:)          = fr_i_tm    (:,:) * r1_ndttrcp1 
# if defined key_trabbl
            IF( nn_bbl_ldf == 1 ) THEN
               ahu_bbl(:,:)      = ahu_bbl_tm (:,:) * r1_ndttrcp1  
               ahv_bbl(:,:)      = ahv_bbl_tm (:,:) * r1_ndttrcp1 
            ENDIF
            IF( nn_bbl_adv == 1 ) THEN
               utr_bbl(:,:)      = utr_bbl_tm (:,:) * r1_ndttrcp1  
               vtr_bbl(:,:)      = vtr_bbl_tm (:,:) * r1_ndttrcp1 
            ENDIF
# endif
         ENDIF
         !
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z1_ne3t = r1_ndttrcp1  / fse3t(ji,jj,jk)
                  z1_ne3u = r1_ndttrcp1  / fse3u(ji,jj,jk)
                  z1_ne3v = r1_ndttrcp1  / fse3v(ji,jj,jk)
                  z1_ne3w = r1_ndttrcp1  / fse3w(ji,jj,jk)
                  !
                  un   (ji,jj,jk)        = un_tm   (ji,jj,jk)        * z1_ne3u
                  vn   (ji,jj,jk)        = vn_tm   (ji,jj,jk)        * z1_ne3v
                  tsn  (ji,jj,jk,jp_tem) = tsn_tm  (ji,jj,jk,jp_tem) * z1_ne3t
                  tsn  (ji,jj,jk,jp_sal) = tsn_tm  (ji,jj,jk,jp_sal) * z1_ne3t
                  rhop (ji,jj,jk)        = rhop_tm (ji,jj,jk)        * z1_ne3t
                  avt  (ji,jj,jk)        = avt_tm  (ji,jj,jk)        * z1_ne3w
# if defined key_zdfddm
                  avs  (ji,jj,jk)        = avs_tm  (ji,jj,jk)        * z1_ne3w
# endif
#if defined key_ldfslp
                  wslpi(ji,jj,jk)        = wslpi_tm(ji,jj,jk) 
                  wslpj(ji,jj,jk)        = wslpj_tm(ji,jj,jk)
                  uslp (ji,jj,jk)        = uslp_tm (ji,jj,jk)
                  vslp (ji,jj,jk)        = vslp_tm (ji,jj,jk)
#endif
               ENDDO
            ENDDO
         ENDDO
         !
         CALL trc_sub_ssh( kt )         ! after ssh & vertical velocity
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sub_stp')
      !
   END SUBROUTINE trc_sub_stp

   SUBROUTINE trc_sub_ini
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sub_ini  ***
      !!                      
      !! ** Purpose : Initialize variables needed for sub-stepping passive tracers
      !! 
      !! ** Method  : 
      !!              Compute the averages for sub-stepping
      !!-------------------------------------------------------------------
      INTEGER ::   ierr
      !!-------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sub_ini')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_sub_ini : initial set up of the passive tracers substepping'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      ierr =  trc_sub_alloc    ()
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'top_sub_alloc : unable to allocate standard ocean arrays' )

      un_tm   (:,:,:)        = un   (:,:,:)        * fse3u(:,:,:) 
      vn_tm   (:,:,:)        = vn   (:,:,:)        * fse3v(:,:,:) 
      tsn_tm  (:,:,:,jp_tem) = tsn  (:,:,:,jp_tem) * fse3t(:,:,:)  
      tsn_tm  (:,:,:,jp_sal) = tsn  (:,:,:,jp_sal) * fse3t(:,:,:)  
      rhop_tm (:,:,:)        = rhop (:,:,:)        * fse3t(:,:,:)  
      avt_tm  (:,:,:)        = avt  (:,:,:)        * fse3w(:,:,:)  
# if defined key_zdfddm
      avs_tm  (:,:,:)        = avs  (:,:,:)        * fse3w(:,:,:)  
# endif
#if defined key_ldfslp
      wslpi_tm(:,:,:)        = wslpi(:,:,:)
      wslpj_tm(:,:,:)        = wslpj(:,:,:)
      uslp_tm (:,:,:)        = uslp (:,:,:)
      vslp_tm (:,:,:)        = vslp (:,:,:)
#endif
      sshn_tm  (:,:) = sshn  (:,:) 
      rnf_tm   (:,:) = rnf   (:,:) 
      h_rnf_tm (:,:) = h_rnf (:,:) 
      hmld_tm  (:,:) = hmld  (:,:)

      ! Physics variables that are set after initialization:
      fr_i_tm(:,:) = 0._wp
      emp_tm (:,:) = 0._wp
      fmmflx_tm(:,:)  = 0._wp
      qsr_tm (:,:) = 0._wp
      wndm_tm(:,:) = 0._wp
# if defined key_trabbl
      IF( nn_bbl_ldf == 1 ) THEN
         ahu_bbl_tm(:,:) = 0._wp
         ahv_bbl_tm(:,:) = 0._wp
      ENDIF
      IF( nn_bbl_adv == 1 ) THEN
         utr_bbl_tm(:,:) = 0._wp
         vtr_bbl_tm(:,:) = 0._wp
      ENDIF
# endif
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sub_ini')
      !
   END SUBROUTINE trc_sub_ini

   SUBROUTINE trc_sub_reset( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sub_reset  ***
      !!                      
      !! ** Purpose : Reset physics variables averaged for substepping
      !! 
      !! ** Method  : 
      !!              Compute the averages for sub-stepping
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      INTEGER :: jk                 ! dummy loop indices
      !!-------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sub_reset')
      !
      !   restore physics variables
      un    (:,:,:)   =  un_temp    (:,:,:)
      vn    (:,:,:)   =  vn_temp    (:,:,:)
      wn    (:,:,:)   =  wn_temp    (:,:,:)
      tsn   (:,:,:,:) =  tsn_temp   (:,:,:,:)
      rhop  (:,:,:)   =  rhop_temp  (:,:,:)
      avt   (:,:,:)   =  avt_temp   (:,:,:)
# if defined key_zdfddm
      avs   (:,:,:)   =  avs_temp   (:,:,:)
# endif
#if defined key_ldfslp
      wslpi (:,:,:)   =  wslpi_temp (:,:,:)
      wslpj (:,:,:)   =  wslpj_temp (:,:,:)
      uslp  (:,:,:)   =  uslp_temp  (:,:,:)
      vslp  (:,:,:)   =  vslp_temp  (:,:,:)
#endif
      sshn  (:,:)     =  sshn_temp  (:,:)
      sshb  (:,:)     =  sshb_temp  (:,:)
      ssha  (:,:)     =  ssha_temp  (:,:)
      rnf   (:,:)     =  rnf_temp   (:,:)
      h_rnf (:,:)     =  h_rnf_temp (:,:)
      !
      hmld  (:,:)     =  hmld_temp  (:,:)
      fr_i  (:,:)     =  fr_i_temp  (:,:)
      emp   (:,:)     =  emp_temp   (:,:)
      fmmflx(:,:)     =  fmmflx_temp(:,:)
      emp_b (:,:)     =  emp_b_temp (:,:)
      qsr   (:,:)     =  qsr_temp   (:,:)
      wndm  (:,:)     =  wndm_temp  (:,:)
# if defined key_trabbl
      IF( nn_bbl_ldf == 1 ) THEN
         ahu_bbl(:,:) = ahu_bbl_temp(:,:) 
         ahv_bbl(:,:) = ahv_bbl_temp(:,:) 
      ENDIF
      IF( nn_bbl_adv == 1 ) THEN
         utr_bbl(:,:) = utr_bbl_temp(:,:) 
         vtr_bbl(:,:) = vtr_bbl_temp(:,:) 
      ENDIF
# endif
      !
      hdivn (:,:,:)   =  hdivn_temp (:,:,:)
      rotn  (:,:,:)   =  rotn_temp  (:,:,:)
      hdivb (:,:,:)   =  hdivb_temp (:,:,:)
      rotb  (:,:,:)   =  rotb_temp  (:,:,:)
      !                                      

      ! Start new averages
         un_tm   (:,:,:)        = un   (:,:,:)        * fse3u(:,:,:) 
         vn_tm   (:,:,:)        = vn   (:,:,:)        * fse3v(:,:,:) 
         tsn_tm  (:,:,:,jp_tem) = tsn  (:,:,:,jp_tem) * fse3t(:,:,:)  
         tsn_tm  (:,:,:,jp_sal) = tsn  (:,:,:,jp_sal) * fse3t(:,:,:)  
         rhop_tm (:,:,:)        = rhop (:,:,:)        * fse3t(:,:,:)  
         avt_tm  (:,:,:)        = avt  (:,:,:)        * fse3w(:,:,:)  
# if defined key_zdfddm
         avs_tm  (:,:,:)        = avs  (:,:,:)        * fse3w(:,:,:)  
# endif
#if defined key_ldfslp
         wslpi_tm(:,:,:)        = wslpi(:,:,:) 
         wslpj_tm(:,:,:)        = wslpj(:,:,:)
         uslp_tm (:,:,:)        = uslp (:,:,:)
         vslp_tm (:,:,:)        = vslp (:,:,:)
#endif
      !
      sshb_hold  (:,:) = sshn  (:,:)
      emp_b_hold (:,:) = emp   (:,:)
      sshn_tm    (:,:) = sshn  (:,:) 
      rnf_tm     (:,:) = rnf   (:,:) 
      h_rnf_tm   (:,:) = h_rnf (:,:) 
      hmld_tm    (:,:) = hmld  (:,:)
      fr_i_tm    (:,:) = fr_i  (:,:)
      emp_tm     (:,:) = emp   (:,:)
      fmmflx_tm  (:,:) = fmmflx(:,:)
      qsr_tm     (:,:) = qsr   (:,:)
      wndm_tm    (:,:) = wndm  (:,:)
# if defined key_trabbl
      IF( nn_bbl_ldf == 1 ) THEN
         ahu_bbl_tm(:,:) = ahu_bbl(:,:) 
         ahv_bbl_tm(:,:) = ahv_bbl(:,:) 
      ENDIF
      IF( nn_bbl_adv == 1 ) THEN
         utr_bbl_tm(:,:) = utr_bbl(:,:) 
         vtr_bbl_tm(:,:) = vtr_bbl(:,:) 
      ENDIF
# endif
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sub_reset')
      !
   END SUBROUTINE trc_sub_reset


   SUBROUTINE trc_sub_ssh( kt ) 
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE trc_sub_ssh  ***
      !!                   
      !! ** Purpose :   compute the after ssh (ssha), the now vertical velocity
      !!              and update the now vertical coordinate (lk_vvl=T).
      !!
      !! ** Method  : - Using the incompressibility hypothesis, the vertical 
      !!      velocity is computed by integrating the horizontal divergence  
      !!      from the bottom to the surface minus the scale factor evolution.
      !!        The boundary conditions are w=0 at the bottom (no flux) and.
      !!
      !! ** action  :   ssha    : after sea surface height
      !!                wn      : now vertical velocity
      !!                sshu_a, sshv_a, sshf_a  : after sea surface height (lk_vvl=T)
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! time step
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcoefu, zcoefv, zcoeff, z2dt, z1_2dt, z1_rau0   ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) :: zhdiv
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sub_ssh')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, zhdiv )

      IF( kt == nittrc000 ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_sub_ssh : after sea surface height and now vertical velocity '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         !
         wn(:,:,jpk) = 0._wp                  ! bottom boundary condition: w=0 (set once for all)
         !
      ENDIF
      !
      CALL div_cur( kt )                              ! Horizontal divergence & Relative vorticity
      !
      z2dt = 2._wp * rdt                              ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == nittrc000 )   z2dt = rdt

      !                                           !------------------------------!
      !                                           !   After Sea Surface Height   !
      !                                           !------------------------------!
      zhdiv(:,:) = 0._wp
      DO jk = 1, jpkm1                                 ! Horizontal divergence of barotropic transports
        zhdiv(:,:) = zhdiv(:,:) + fse3t(:,:,jk) * hdivn(:,:,jk)
      END DO
      !                                                ! Sea surface elevation time stepping
      ! In forward Euler time stepping case, the same formulation as in the leap-frog case can be used
      ! because emp_b field is initialized with the vlaues of emp field. Hence, 0.5 * ( emp + emp_b ) = emp
      z1_rau0 = 0.5 / rau0
      ssha(:,:) = (  sshb(:,:) - z2dt * ( z1_rau0 * ( emp_b(:,:) + emp(:,:) ) + zhdiv(:,:) )  ) * tmask(:,:,1)
#if ! defined key_dynspg_ts
      ! These lines are not necessary with time splitting since
      ! boundary condition on sea level is set during ts loop
#if defined key_agrif
      CALL agrif_ssh( kt )
#endif
#if defined key_bdy
      ssha(:,:) = ssha(:,:) * bdytmask(:,:)
      CALL lbc_lnk( ssha, 'T', 1. ) 
#endif
#endif


      !                                           !------------------------------!
      !                                           !     Now Vertical Velocity    !
      !                                           !------------------------------!
      z1_2dt = 1.e0 / z2dt
      DO jk = jpkm1, 1, -1                             ! integrate from the bottom the hor. divergence
         ! - ML - need 3 lines here because replacement of fse3t by its expression yields too long lines otherwise
         wn(:,:,jk) = wn(:,:,jk+1) -   fse3t_n(:,:,jk) * hdivn(:,:,jk)        &
            &                      - ( fse3t_a(:,:,jk) - fse3t_b(:,:,jk) )    &
            &                         * tmask(:,:,jk) * z1_2dt
#if defined key_bdy
         wn(:,:,jk) = wn(:,:,jk) * bdytmask(:,:)
#endif
      END DO

      !
      CALL wrk_dealloc( jpi, jpj, zhdiv )
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sub_ssh')
      !
   END SUBROUTINE trc_sub_ssh

   INTEGER FUNCTION trc_sub_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE trc_sub_alloc ***
      !!-------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_warn
      INTEGER ::  ierr
      !!-------------------------------------------------------------------
      !
      ALLOCATE( un_temp(jpi,jpj,jpk)        ,  vn_temp(jpi,jpj,jpk)  ,   &
         &      wn_temp(jpi,jpj,jpk)        ,  avt_temp(jpi,jpj,jpk) ,   &
         &      rhop_temp(jpi,jpj,jpk)      ,  rhop_tm(jpi,jpj,jpk) ,   &
         &      sshn_temp(jpi,jpj)          ,  sshb_temp(jpi,jpj) ,      &
         &      ssha_temp(jpi,jpj)          ,                           &
#if defined key_ldfslp
         &      wslpi_temp(jpi,jpj,jpk)     ,  wslpj_temp(jpi,jpj,jpk),  &
         &      uslp_temp(jpi,jpj,jpk)      ,  vslp_temp(jpi,jpj,jpk),   &
#endif
#if defined key_trabbl
         &      ahu_bbl_temp(jpi,jpj)       ,  ahv_bbl_temp(jpi,jpj),    &
         &      utr_bbl_temp(jpi,jpj)       ,  vtr_bbl_temp(jpi,jpj),    &
#endif
         &      rnf_temp(jpi,jpj)           ,  h_rnf_temp(jpi,jpj) ,     &
         &      tsn_temp(jpi,jpj,jpk,2)     ,  emp_b_temp(jpi,jpj),      &
         &      emp_temp(jpi,jpj)           ,  fmmflx_temp(jpi,jpj),     &
         &      hmld_temp(jpi,jpj)          ,  qsr_temp(jpi,jpj) ,       &
         &      fr_i_temp(jpi,jpj)          ,  fr_i_tm(jpi,jpj) ,        &
         &      wndm_temp(jpi,jpj)          ,  wndm_tm(jpi,jpj) ,        &
# if defined key_zdfddm
         &      avs_tm(jpi,jpj,jpk)         ,  avs_temp(jpi,jpj,jpk) ,   &
# endif
         &      hdivn_temp(jpi,jpj,jpk)     ,  hdivb_temp(jpi,jpj,jpk),  &
         &      rotn_temp(jpi,jpj,jpk)      ,  rotb_temp(jpi,jpj,jpk),   &
         &      un_tm(jpi,jpj,jpk)          ,  vn_tm(jpi,jpj,jpk)  ,     &
         &      avt_tm(jpi,jpj,jpk)                                ,     &
         &      sshn_tm(jpi,jpj)            ,  sshb_hold(jpi,jpj) ,      &
         &      tsn_tm(jpi,jpj,jpk,2)       ,                            &
         &      emp_tm(jpi,jpj)             ,  fmmflx_tm(jpi,jpj)  ,     &
         &      emp_b_hold(jpi,jpj)         ,                            &
         &      hmld_tm(jpi,jpj)            ,  qsr_tm(jpi,jpj) ,         &
#if defined key_ldfslp
         &      wslpi_tm(jpi,jpj,jpk)       ,  wslpj_tm(jpi,jpj,jpk),    &
         &      uslp_tm(jpi,jpj,jpk)        ,  vslp_tm(jpi,jpj,jpk),     &
#endif
#if defined key_trabbl
         &      ahu_bbl_tm(jpi,jpj)         ,  ahv_bbl_tm(jpi,jpj),      &
         &      utr_bbl_tm(jpi,jpj)         ,  vtr_bbl_tm(jpi,jpj),      &
#endif
         &      rnf_tm(jpi,jpj)             ,  h_rnf_tm(jpi,jpj) ,       &
         &                                    STAT=trc_sub_alloc )  
      IF( trc_sub_alloc /= 0 )   CALL ctl_warn('trc_sub_alloc: failed to allocate arrays')

      !
   END FUNCTION trc_sub_alloc

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sub_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_sub_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sub_stp
   SUBROUTINE trc_sub_ini        ! Empty routine
      WRITE(*,*) 'trc_sub_ini: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sub_ini

#endif

   !!======================================================================
END MODULE trcsub
