#if defined key_agrif
!!----------------------------------------------------------------------
!! NEMO/NST 3.4 , NEMO Consortium (2012)
!! $Id: agrif_user.F90 6204 2016-01-04 13:47:06Z cetlod $
!! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
SUBROUTINE agrif_user
END SUBROUTINE agrif_user

SUBROUTINE agrif_before_regridding
END SUBROUTINE agrif_before_regridding

SUBROUTINE Agrif_InitWorkspace
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_InitWorkspace ***
   !!----------------------------------------------------------------------
   USE par_oce
   USE dom_oce
   USE nemogcm
   !
   IMPLICIT NONE
   !!----------------------------------------------------------------------
   !
   IF( .NOT. Agrif_Root() ) THEN
      jpni = Agrif_Parent(jpni)
      jpnj = Agrif_Parent(jpnj)
      jpnij = Agrif_Parent(jpnij)
      jpiglo  = nbcellsx + 2 + 2*nbghostcells
      jpjglo  = nbcellsy + 2 + 2*nbghostcells
      jpi     = ( jpiglo-2*jpreci + (jpni-1+0) ) / jpni + 2*jpreci
      jpj     = ( jpjglo-2*jprecj + (jpnj-1+0) ) / jpnj + 2*jprecj
! JC: change to allow for different vertical levels
!     jpk is already set
!     keep it jpk possibly different from jpkdta which 
!     hold parent grid vertical levels number (set earlier)
!      jpk     = jpkdta 
      jpim1   = jpi-1 
      jpjm1   = jpj-1 
      jpkm1   = jpk-1                                         
      jpij    = jpi*jpj 
      jpidta  = jpiglo
      jpjdta  = jpjglo
      jpizoom = 1
      jpjzoom = 1
      nperio  = 0
      jperio  = 0
   ENDIF
   !
END SUBROUTINE Agrif_InitWorkspace


SUBROUTINE Agrif_InitValues
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_InitValues ***
   !!
   !! ** Purpose :: Declaration of variables to be interpolated
   !!----------------------------------------------------------------------
   USE Agrif_Util
   USE oce 
   USE dom_oce
   USE nemogcm
   USE tradmp
   USE bdy_par

   IMPLICIT NONE
   !!----------------------------------------------------------------------
   ! 0. Initializations
   !-------------------
   IF( cp_cfg == 'orca' ) THEN
      IF ( jp_cfg == 2 .OR. jp_cfg == 025 .OR. jp_cfg == 05 &
            &                      .OR. jp_cfg == 4 ) THEN
         jp_cfg = -1    ! set special value for jp_cfg on fine grids
         cp_cfg = "default"
      ENDIF
   ENDIF
   ! Specific fine grid Initializations
   ! no tracer damping on fine grids
   ln_tradmp = .FALSE.
   ! no open boundary on fine grids
   lk_bdy = .FALSE.


   CALL nemo_init  ! Initializations of each fine grid

   CALL agrif_nemo_init
    IF(lwp) WRITE(numout,*)'Agrif_after_init'
   CALL Agrif_InitValues_cont_dom
   IF(lwp) WRITE(numout,*)'Agrif_after_initvalues'
# if ! defined key_offline
   CALL Agrif_InitValues_cont
# endif      
 IF(lwp) WRITE(numout,*)'Agrif before top'
# if defined key_top
!CEB   CALL Agrif_InitValues_cont_top
# endif     
 IF(lwp) WRITE(numout,*)'Agrif end'
END SUBROUTINE Agrif_initvalues

SUBROUTINE Agrif_InitValues_cont_dom
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_InitValues_cont ***
   !!
   !! ** Purpose ::   Declaration of variables to be interpolated
   !!----------------------------------------------------------------------
   USE Agrif_Util
   USE oce 
   USE dom_oce
   USE nemogcm
   USE sol_oce
   USE in_out_manager
   USE agrif_opa_update
   USE agrif_opa_interp
   USE agrif_opa_sponge
   !
   IMPLICIT NONE
   !
   !!----------------------------------------------------------------------

   ! Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_var_dom
   !
END SUBROUTINE Agrif_InitValues_cont_dom


SUBROUTINE agrif_declare_var_dom
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE agrif_declare_var ***
   !!
   !! ** Purpose :: Declaration of variables to be interpolated
   !!----------------------------------------------------------------------
   USE agrif_util
   USE par_oce       
   USE oce
   IMPLICIT NONE
   !!----------------------------------------------------------------------

   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_variable((/1,2/),(/2,3/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),e1u_id)
   CALL agrif_declare_variable((/2,1/),(/3,2/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),e2v_id)

   ! 2. Type of interpolation
   !-------------------------
   CALL Agrif_Set_bcinterp(e1u_id,interp1=Agrif_linear,interp2=AGRIF_ppm)
   CALL Agrif_Set_bcinterp(e2v_id,interp1=AGRIF_ppm,interp2=Agrif_linear)

   ! 3. Location of interpolation
   !-----------------------------
   CALL Agrif_Set_bc(e1u_id,(/0,0/))
   CALL Agrif_Set_bc(e2v_id,(/0,0/))

   ! 5. Update type
   !--------------- 
   CALL Agrif_Set_Updatetype(e1u_id,update1 = Agrif_Update_Copy, update2=Agrif_Update_Average)
   CALL Agrif_Set_Updatetype(e2v_id,update1 = Agrif_Update_Average, update2=Agrif_Update_Copy)

! High order updates
!   CALL Agrif_Set_Updatetype(e1u_id,update1 = Agrif_Update_Average,            update2=Agrif_Update_Full_Weighting)
!   CALL Agrif_Set_Updatetype(e2v_id,update1 = Agrif_Update_Full_Weighting,     update2=Agrif_Update_Average)
    !
END SUBROUTINE agrif_declare_var_dom


# if ! defined key_offline

SUBROUTINE Agrif_InitValues_cont
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_InitValues_cont ***
   !!
   !! ** Purpose ::   Declaration of variables to be interpolated
   !!----------------------------------------------------------------------
   USE Agrif_Util
   USE oce 
   USE dom_oce
   USE nemogcm
   USE sol_oce
   USE lib_mpp
   USE in_out_manager
   USE agrif_opa_update
   USE agrif_opa_interp
   USE agrif_opa_sponge
   !
   IMPLICIT NONE
   !
   LOGICAL :: check_namelist
   CHARACTER(len=15) :: cl_check1, cl_check2, cl_check3
   !!----------------------------------------------------------------------

   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_var
   IF(lwp) WRITE(numout,*)'Agrif after declare'
   ! 2. First interpolations of potentially non zero fields
   !-------------------------------------------------------
   Agrif_SpecialValue=0.
   Agrif_UseSpecialValue = .TRUE.
   CALL Agrif_Bc_variable(tsn_id,calledweight=1.,procname=interptsn)
   CALL Agrif_Sponge
   tabspongedone_tsn = .FALSE.
   CALL Agrif_Bc_variable(tsn_sponge_id,calledweight=1.,procname=interptsn_sponge)
   ! reset tsa to zero
   tsa(:,:,:,:) = 0.

   Agrif_UseSpecialValue = ln_spc_dyn
   CALL Agrif_Bc_variable(un_interp_id,calledweight=1.,procname=interpun)
   CALL Agrif_Bc_variable(vn_interp_id,calledweight=1.,procname=interpvn)
   tabspongedone_u = .FALSE.
   tabspongedone_v = .FALSE.
   CALL Agrif_Bc_variable(un_sponge_id,calledweight=1.,procname=interpun_sponge)
   tabspongedone_u = .FALSE.
   tabspongedone_v = .FALSE.
   CALL Agrif_Bc_variable(vn_sponge_id,calledweight=1.,procname=interpvn_sponge)

#if defined key_dynspg_ts
   Agrif_UseSpecialValue = .TRUE.
   CALL Agrif_Bc_variable(sshn_id,calledweight=1., procname=interpsshn )

   Agrif_UseSpecialValue = ln_spc_dyn
   CALL Agrif_Bc_variable(unb_id,calledweight=1.,procname=interpunb)
   CALL Agrif_Bc_variable(vnb_id,calledweight=1.,procname=interpvnb)
   CALL Agrif_Bc_variable(ub2b_interp_id,calledweight=1.,procname=interpub2b)
   CALL Agrif_Bc_variable(vb2b_interp_id,calledweight=1.,procname=interpvb2b)
   ubdy_w(:) = 0.e0 ; vbdy_w(:) = 0.e0 ; hbdy_w(:) =0.e0
   ubdy_e(:) = 0.e0 ; vbdy_e(:) = 0.e0 ; hbdy_e(:) =0.e0 
   ubdy_n(:) = 0.e0 ; vbdy_n(:) = 0.e0 ; hbdy_n(:) =0.e0 
   ubdy_s(:) = 0.e0 ; vbdy_s(:) = 0.e0 ; hbdy_s(:) =0.e0
#endif

   Agrif_UseSpecialValue = .FALSE. 
   ! reset velocities to zero
   ua(:,:,:) = 0.
   va(:,:,:) = 0.
   IF(lwp) WRITE(numout,*)'Agrif after 2'

   ! 3. Some controls
   !-----------------
   check_namelist = .TRUE.

   IF( check_namelist ) THEN 

      ! Check time steps           
      IF( NINT(Agrif_Rhot()) * NINT(rdt) .NE. Agrif_Parent(rdt) ) THEN
         WRITE(cl_check1,*)  NINT(Agrif_Parent(rdt))
         WRITE(cl_check2,*)  NINT(rdt)
         WRITE(cl_check3,*)  NINT(Agrif_Parent(rdt)/Agrif_Rhot())
         CALL ctl_warn( 'incompatible time step between grids',   &
               &               'parent grid value : '//cl_check1    ,   & 
               &               'child  grid value : '//cl_check2    ,   & 
               &               'value on child grid will be changed to : '//cl_check3 )
         rdt=Agrif_Parent(rdt)/Agrif_Rhot()
      ENDIF

      ! Check run length
      IF( Agrif_IRhot() * (Agrif_Parent(nitend)- &
            Agrif_Parent(nit000)+1) .NE. (nitend-nit000+1) ) THEN
         WRITE(cl_check1,*)  (Agrif_Parent(nit000)-1)*Agrif_IRhot() + 1
         WRITE(cl_check2,*)   Agrif_Parent(nitend)   *Agrif_IRhot()
         CALL ctl_warn( 'incompatible run length between grids'               ,   &
               &              ' nit000 on fine grid will be change to : '//cl_check1,   &
               &              ' nitend on fine grid will be change to : '//cl_check2    )
         nit000 = (Agrif_Parent(nit000)-1)*Agrif_IRhot() + 1
         nitend =  Agrif_Parent(nitend)   *Agrif_IRhot()
      ENDIF

      ! Check coordinates
      IF( ln_zps ) THEN
         ! check parameters for partial steps 
         IF( Agrif_Parent(e3zps_min) .NE. e3zps_min ) THEN
            WRITE(*,*) 'incompatible e3zps_min between grids'
            WRITE(*,*) 'parent grid :',Agrif_Parent(e3zps_min)
            WRITE(*,*) 'child grid  :',e3zps_min
            WRITE(*,*) 'those values should be identical'
            STOP
         ENDIF
         IF( Agrif_Parent(e3zps_rat) /= e3zps_rat ) THEN
            WRITE(*,*) 'incompatible e3zps_rat between grids'
            WRITE(*,*) 'parent grid :',Agrif_Parent(e3zps_rat)
            WRITE(*,*) 'child grid  :',e3zps_rat
            WRITE(*,*) 'those values should be identical'                  
            STOP
         ENDIF
      ENDIF
      ! check if masks and bathymetries match
      IF(ln_chk_bathy) THEN
         !
         IF(lwp) WRITE(numout,*) 'AGRIF: Check Bathymetry and masks near bdys. Level: ', Agrif_Level()
         !
         kindic_agr = 0
         ! check if umask agree with parent along western and eastern boundaries:
         CALL Agrif_Bc_variable(umsk_id,calledweight=1.,procname=interpumsk)
         ! check if vmask agree with parent along northern and southern boundaries:
         CALL Agrif_Bc_variable(vmsk_id,calledweight=1.,procname=interpvmsk)
	 ! check if tmask and vertical scale factors agree with parent over first two coarse grid points:
         CALL Agrif_Bc_variable(e3t_id,calledweight=1.,procname=interpe3t)
         !
         IF (lk_mpp) CALL mpp_sum( kindic_agr )
         IF( kindic_agr /= 0 ) THEN             		
            CALL ctl_stop('Child Bathymetry is not correct near boundaries.')
         ELSE
            IF(lwp) WRITE(numout,*) 'Child Bathymetry is ok near boundaries.'
         END IF
      ENDIF
      !
   ENDIF
   IF(lwp) WRITE(numout,*)'Agrif after checks'

   ! 
   ! Do update at initialisation because not done before writing restarts
   ! This would indeed change boundary conditions values at initial time
   ! hence produce restartability issues.
   ! Note that update below is recursive (with lk_agrif_doupd=T):
   ! 
! JC: I am not sure if Agrif_MaxLevel() is the "relative"
!     or the absolute maximum nesting level...TBC                        
   IF ( Agrif_Level().EQ.Agrif_MaxLevel() ) THEN 
      ! NB: Do tracers first, dynamics after because nbcline incremented in dynamics
      CALL Agrif_Update_tra()
      CALL Agrif_Update_dyn()
   ENDIF
   !
# if defined key_zdftke
!   CALL Agrif_Update_tke(0)
# endif
   !
   Agrif_UseSpecialValueInUpdate = .FALSE.
   nbcline = 0
   lk_agrif_doupd = .FALSE.
   !
END SUBROUTINE Agrif_InitValues_cont


SUBROUTINE agrif_declare_var
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE agrif_declarE_var ***
   !!
   !! ** Purpose :: Declaration of variables to be interpolated
   !!----------------------------------------------------------------------
   USE agrif_util
   USE par_oce       !   ONLY : jpts
   USE oce
   USE agrif_oce
   IMPLICIT NONE
   !!----------------------------------------------------------------------

   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_variable((/2,2,0,0/),(/3,3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/nlci,nlcj,jpk,jpts/),tsn_id)
   CALL agrif_declare_variable((/2,2,0,0/),(/3,3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/nlci,nlcj,jpk,jpts/),tsn_sponge_id)

   CALL agrif_declare_variable((/1,2,0/),(/2,3,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),un_interp_id)
   CALL agrif_declare_variable((/2,1,0/),(/3,2,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),vn_interp_id)
   CALL agrif_declare_variable((/1,2,0/),(/2,3,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),un_update_id)
   CALL agrif_declare_variable((/2,1,0/),(/3,2,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),vn_update_id)
   CALL agrif_declare_variable((/1,2,0/),(/2,3,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),un_sponge_id)
   CALL agrif_declare_variable((/2,1,0/),(/3,2,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),vn_sponge_id)

   CALL agrif_declare_variable((/2,2,0/),(/3,3,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),e3t_id)
   CALL agrif_declare_variable((/1,2,0/),(/2,3,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),umsk_id)
   CALL agrif_declare_variable((/2,1,0/),(/3,2,0/),(/'x','y','N'/),(/1,1,1/),(/nlci,nlcj,jpk/),vmsk_id)

   CALL agrif_declare_variable((/2,2,0,0/),(/3,3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/nlci,nlcj,jpk,3/),scales_t_id)

   CALL agrif_declare_variable((/1,2/),(/2,3/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),unb_id)
   CALL agrif_declare_variable((/2,1/),(/3,2/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),vnb_id)
   CALL agrif_declare_variable((/1,2/),(/2,3/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),ub2b_interp_id)
   CALL agrif_declare_variable((/2,1/),(/3,2/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),vb2b_interp_id)
   CALL agrif_declare_variable((/1,2/),(/2,3/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),ub2b_update_id)
   CALL agrif_declare_variable((/2,1/),(/3,2/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),vb2b_update_id)

   CALL agrif_declare_variable((/2,2/),(/3,3/),(/'x','y'/),(/1,1/),(/nlci,nlcj/),sshn_id)

# if defined key_zdftke
   CALL agrif_declare_variable((/2,2,0/),(/3,3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj,jpk/), en_id)
   CALL agrif_declare_variable((/2,2,0/),(/3,3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj,jpk/),avt_id)
   CALL agrif_declare_variable((/2,2,0/),(/3,3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj,jpk/),avm_id)
# endif

   ! 2. Type of interpolation
   !-------------------------
   CALL Agrif_Set_bcinterp(tsn_id,interp=AGRIF_linear)

   CALL Agrif_Set_bcinterp(un_interp_id,interp1=Agrif_linear,interp2=AGRIF_ppm)
   CALL Agrif_Set_bcinterp(vn_interp_id,interp1=AGRIF_ppm,interp2=Agrif_linear)

   CALL Agrif_Set_bcinterp(tsn_sponge_id,interp=AGRIF_linear)

   CALL Agrif_Set_bcinterp(sshn_id,interp=AGRIF_linear)
   CALL Agrif_Set_bcinterp(unb_id,interp1=Agrif_linear,interp2=AGRIF_ppm)
   CALL Agrif_Set_bcinterp(vnb_id,interp1=AGRIF_ppm,interp2=Agrif_linear)
   CALL Agrif_Set_bcinterp(ub2b_interp_id,interp1=Agrif_linear,interp2=AGRIF_ppm)
   CALL Agrif_Set_bcinterp(vb2b_interp_id,interp1=AGRIF_ppm,interp2=Agrif_linear)


   CALL Agrif_Set_bcinterp(un_sponge_id,interp1=Agrif_linear,interp2=AGRIF_ppm)
   CALL Agrif_Set_bcinterp(vn_sponge_id,interp1=AGRIF_ppm,interp2=Agrif_linear)

   CALL Agrif_Set_bcinterp(e3t_id,interp=AGRIF_constant)
   CALL Agrif_Set_bcinterp(umsk_id,interp=AGRIF_constant)
   CALL Agrif_Set_bcinterp(vmsk_id,interp=AGRIF_constant)

# if defined key_zdftke
   CALL Agrif_Set_bcinterp(avm_id ,interp=AGRIF_linear)
# endif


   ! 3. Location of interpolation
   !-----------------------------
   CALL Agrif_Set_bc(tsn_id,(/0,1/))
   CALL Agrif_Set_bc(un_interp_id,(/0,1/))
   CALL Agrif_Set_bc(vn_interp_id,(/0,1/))

!   CALL Agrif_Set_bc(tsn_sponge_id,(/-3*Agrif_irhox(),0/))
!   CALL Agrif_Set_bc(un_sponge_id,(/-2*Agrif_irhox()-1,0/))
!   CALL Agrif_Set_bc(vn_sponge_id,(/-2*Agrif_irhox()-1,0/))
   CALL Agrif_Set_bc(tsn_sponge_id,(/-nn_sponge_len*Agrif_irhox()-1,0/))
   CALL Agrif_Set_bc(un_sponge_id ,(/-nn_sponge_len*Agrif_irhox()-1,0/))
   CALL Agrif_Set_bc(vn_sponge_id ,(/-nn_sponge_len*Agrif_irhox()-1,0/))

   CALL Agrif_Set_bc(sshn_id,(/0,0/))
   CALL Agrif_Set_bc(unb_id ,(/0,0/))
   CALL Agrif_Set_bc(vnb_id ,(/0,0/))
   CALL Agrif_Set_bc(ub2b_interp_id,(/0,0/))
   CALL Agrif_Set_bc(vb2b_interp_id,(/0,0/))

   CALL Agrif_Set_bc(e3t_id,(/-2*Agrif_irhox()-1,0/))   ! if west and rhox=3: column 2 to 9
   CALL Agrif_Set_bc(umsk_id,(/0,0/))
   CALL Agrif_Set_bc(vmsk_id,(/0,0/))

# if defined key_zdftke
   CALL Agrif_Set_bc(avm_id ,(/0,1/))
# endif

   ! 5. Update type
   !--------------- 
   CALL Agrif_Set_Updatetype(tsn_id, update = AGRIF_Update_Average)

   CALL Agrif_Set_Updatetype(scales_t_id, update = AGRIF_Update_Average)

   CALL Agrif_Set_Updatetype(un_update_id,update1 = Agrif_Update_Copy, update2 = Agrif_Update_Average)
   CALL Agrif_Set_Updatetype(vn_update_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy)

   CALL Agrif_Set_Updatetype(sshn_id, update = AGRIF_Update_Average)

   CALL Agrif_Set_Updatetype(ub2b_update_id,update1 = Agrif_Update_Copy, update2 = Agrif_Update_Average)
   CALL Agrif_Set_Updatetype(vb2b_update_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy)

# if defined key_zdftke
   CALL Agrif_Set_Updatetype( en_id, update = AGRIF_Update_Average)
   CALL Agrif_Set_Updatetype(avt_id, update = AGRIF_Update_Average)
   CALL Agrif_Set_Updatetype(avm_id, update = AGRIF_Update_Average)
# endif

! High order updates
!   CALL Agrif_Set_Updatetype(tsn_id, update = Agrif_Update_Full_Weighting)
!   CALL Agrif_Set_Updatetype(un_update_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Full_Weighting)
!   CALL Agrif_Set_Updatetype(vn_update_id,update1 = Agrif_Update_Full_Weighting, update2 = Agrif_Update_Average)
!
!   CALL Agrif_Set_Updatetype(ub2b_update_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Full_Weighting)
!   CALL Agrif_Set_Updatetype(vb2b_update_id,update1 = Agrif_Update_Full_Weighting, update2 = Agrif_Update_Average)
!   CALL Agrif_Set_Updatetype(sshn_id, update = Agrif_Update_Full_Weighting)
 
   !
END SUBROUTINE agrif_declare_var
# endif

#  if defined key_lim2
SUBROUTINE Agrif_InitValues_cont_lim2
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_InitValues_cont_lim2 ***
   !!
   !! ** Purpose :: Initialisation of variables to be interpolated for LIM2
   !!----------------------------------------------------------------------
   USE Agrif_Util
   USE ice_2
   USE agrif_ice
   USE in_out_manager
   USE agrif_lim2_update
   USE agrif_lim2_interp
   USE lib_mpp
   !
   IMPLICIT NONE
   !
   !!----------------------------------------------------------------------

   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_var_lim2

   ! 2. First interpolations of potentially non zero fields
   !-------------------------------------------------------
   Agrif_SpecialValue=-9999.
   Agrif_UseSpecialValue = .TRUE.
   !     Call Agrif_Bc_variable(zadv ,adv_ice_id ,calledweight=1.,procname=interp_adv_ice )
   !     Call Agrif_Bc_variable(zvel ,u_ice_id   ,calledweight=1.,procname=interp_u_ice   )
   !     Call Agrif_Bc_variable(zvel ,v_ice_id   ,calledweight=1.,procname=interp_v_ice   )
   Agrif_SpecialValue=0.
   Agrif_UseSpecialValue = .FALSE.

   ! 3. Some controls
   !-----------------

#   if ! defined key_lim2_vp
   lim_nbstep = 1.
   CALL agrif_rhg_lim2_load
   CALL agrif_trp_lim2_load
   lim_nbstep = 0.
#   endif
   !RB mandatory but why ???
   !      IF( nbclineupdate /= nn_fsbc .AND. nn_ice == 2 )THEN
   !         CALL ctl_warn ('With ice model on child grid, nbclineupdate is set to nn_fsbc')
   !         nbclineupdate = nn_fsbc
   !       ENDIF
   CALL Agrif_Update_lim2(0)
   !
END SUBROUTINE Agrif_InitValues_cont_lim2

SUBROUTINE agrif_declare_var_lim2
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE agrif_declare_var_lim2 ***
   !!
   !! ** Purpose :: Declaration of variables to be interpolated for LIM2
   !!----------------------------------------------------------------------
   USE agrif_util
   USE ice_2

   IMPLICIT NONE
   !!----------------------------------------------------------------------

   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_variable((/2,2,0/),(/3,3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj, 7/),adv_ice_id )
#   if defined key_lim2_vp
   CALL agrif_declare_variable((/1,1/),(/3,3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),u_ice_id)
   CALL agrif_declare_variable((/1,1/),(/3,3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),v_ice_id)
#   else
   CALL agrif_declare_variable((/1,2/),(/2,3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),u_ice_id)
   CALL agrif_declare_variable((/2,1/),(/3,2/),(/'x','y'/),(/1,1/),(/jpi,jpj/),v_ice_id)
#   endif

   ! 2. Type of interpolation
   !-------------------------
   CALL Agrif_Set_bcinterp(adv_ice_id ,interp=AGRIF_linear)
   CALL Agrif_Set_bcinterp(u_ice_id,interp1=Agrif_linear,interp2=AGRIF_ppm)
   CALL Agrif_Set_bcinterp(v_ice_id,interp1=AGRIF_ppm,interp2=Agrif_linear)

   ! 3. Location of interpolation
   !-----------------------------
   CALL Agrif_Set_bc(adv_ice_id ,(/0,1/))
   CALL Agrif_Set_bc(u_ice_id,(/0,1/))
   CALL Agrif_Set_bc(v_ice_id,(/0,1/))

   ! 5. Update type
   !---------------
   CALL Agrif_Set_Updatetype(adv_ice_id , update = AGRIF_Update_Average)
   CALL Agrif_Set_Updatetype(u_ice_id,update1 = Agrif_Update_Copy, update2 = Agrif_Update_Average)
   CALL Agrif_Set_Updatetype(v_ice_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy)
   ! 
END SUBROUTINE agrif_declare_var_lim2
#  endif


# if defined key_top
SUBROUTINE Agrif_InitValues_cont_top
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_InitValues_cont_top ***
   !!
   !! ** Purpose :: Declaration of variables to be interpolated
   !!----------------------------------------------------------------------
   USE Agrif_Util
   USE oce 
   USE dom_oce
   USE nemogcm
   USE par_trc
   USE lib_mpp
   USE trc
   USE in_out_manager
   USE agrif_opa_sponge
   USE agrif_top_update
   USE agrif_top_interp
   USE agrif_top_sponge
   !
   IMPLICIT NONE
   !
   CHARACTER(len=10) :: cl_check1, cl_check2, cl_check3
   LOGICAL :: check_namelist
   !!----------------------------------------------------------------------


   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_var_top

   ! 2. First interpolations of potentially non zero fields
   !-------------------------------------------------------
   Agrif_SpecialValue=0.
   Agrif_UseSpecialValue = .TRUE.
   CALL Agrif_Bc_variable(trn_id,calledweight=1.,procname=interptrn)
   Agrif_UseSpecialValue = .FALSE.
   CALL Agrif_Sponge
   tabspongedone_trn = .FALSE.
   CALL Agrif_Bc_variable(trn_sponge_id,calledweight=1.,procname=interptrn_sponge)
   ! reset tsa to zero
   tra(:,:,:,:) = 0.


   ! 3. Some controls
   !-----------------
   check_namelist = .TRUE.

   IF( check_namelist ) THEN
# if defined key_offline
      ! Check time steps
      IF( NINT(Agrif_Rhot()) * NINT(rdt) .NE. Agrif_Parent(rdt) ) THEN
         WRITE(cl_check1,*)  Agrif_Parent(rdt)
         WRITE(cl_check2,*)  rdt
         WRITE(cl_check3,*)  rdt*Agrif_Rhot()
         CALL ctl_warn( 'incompatible time step between grids',   &
               &               'parent grid value : '//cl_check1    ,   & 
               &               'child  grid value : '//cl_check2    ,   & 
               &               'value on child grid will be changed to  &
               &               :'//cl_check3  )
         rdt=rdt*Agrif_Rhot()
      ENDIF

      ! Check run length
      IF( Agrif_IRhot() * (Agrif_Parent(nitend)- &
            Agrif_Parent(nit000)+1) .NE. (nitend-nit000+1) ) THEN
         WRITE(cl_check1,*)  (Agrif_Parent(nit000)-1)*Agrif_IRhot() + 1
         WRITE(cl_check2,*)   Agrif_Parent(nitend)   *Agrif_IRhot()
         CALL ctl_warn( 'incompatible run length between grids'               ,   &
               &              ' nit000 on fine grid will be change to : '//cl_check1,   &
               &              ' nitend on fine grid will be change to : '//cl_check2    )
         nit000 = (Agrif_Parent(nit000)-1)*Agrif_IRhot() + 1
         nitend =  Agrif_Parent(nitend)   *Agrif_IRhot()
      ENDIF

      ! Check coordinates
      IF( ln_zps ) THEN
         ! check parameters for partial steps 
         IF( Agrif_Parent(e3zps_min) .NE. e3zps_min ) THEN
            WRITE(*,*) 'incompatible e3zps_min between grids'
            WRITE(*,*) 'parent grid :',Agrif_Parent(e3zps_min)
            WRITE(*,*) 'child grid  :',e3zps_min
            WRITE(*,*) 'those values should be identical'
            STOP
         ENDIF
         IF( Agrif_Parent(e3zps_rat) .NE. e3zps_rat ) THEN
            WRITE(*,*) 'incompatible e3zps_rat between grids'
            WRITE(*,*) 'parent grid :',Agrif_Parent(e3zps_rat)
            WRITE(*,*) 'child grid  :',e3zps_rat
            WRITE(*,*) 'those values should be identical'                  
            STOP
         ENDIF
      ENDIF
#  endif         
      ! Check passive tracer cell
      IF( nn_dttrc .NE. 1 ) THEN
         WRITE(*,*) 'nn_dttrc should be equal to 1'
      ENDIF
   ENDIF

   CALL Agrif_Update_trc(0)
   !
   Agrif_UseSpecialValueInUpdate = .FALSE.
   nbcline_trc = 0
   !
END SUBROUTINE Agrif_InitValues_cont_top


SUBROUTINE agrif_declare_var_top
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE agrif_declare_var_top ***
   !!
   !! ** Purpose :: Declaration of TOP variables to be interpolated
   !!----------------------------------------------------------------------
   USE agrif_util
   USE agrif_oce
   USE dom_oce
   USE trc

   IMPLICIT NONE

   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
   CALL agrif_declare_variable((/2,2,0,0/),(/3,3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/nlci,nlcj,jpk,jptra/),trn_id)
   CALL agrif_declare_variable((/2,2,0,0/),(/3,3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/nlci,nlcj,jpk,jptra/),trn_sponge_id)

   ! 2. Type of interpolation
   !-------------------------
   CALL Agrif_Set_bcinterp(trn_id,interp=AGRIF_linear)
   CALL Agrif_Set_bcinterp(trn_sponge_id,interp=AGRIF_linear)

   ! 3. Location of interpolation
   !-----------------------------
   CALL Agrif_Set_bc(trn_id,(/0,1/))
!   CALL Agrif_Set_bc(trn_sponge_id,(/-3*Agrif_irhox(),0/))
   CALL Agrif_Set_bc(trn_sponge_id,(/-nn_sponge_len*Agrif_irhox()-1,0/))

   ! 5. Update type
   !--------------- 
   CALL Agrif_Set_Updatetype(trn_id, update = AGRIF_Update_Average)

!   Higher order update
!   CALL Agrif_Set_Updatetype(tsn_id, update = Agrif_Update_Full_Weighting)

   !
END SUBROUTINE agrif_declare_var_top
# endif

SUBROUTINE Agrif_detect( kg, ksizex )
   !!----------------------------------------------------------------------
   !!   *** ROUTINE Agrif_detect ***
   !!----------------------------------------------------------------------
   !
   INTEGER, DIMENSION(2) :: ksizex
   INTEGER, DIMENSION(ksizex(1),ksizex(2)) :: kg 
   !!----------------------------------------------------------------------
   !
   RETURN
   !
END SUBROUTINE Agrif_detect


SUBROUTINE agrif_nemo_init
   !!----------------------------------------------------------------------
   !!                     *** ROUTINE agrif_init ***
   !!----------------------------------------------------------------------
   USE agrif_oce 
   USE agrif_ice
   USE in_out_manager
   USE lib_mpp
   IMPLICIT NONE
   !
   INTEGER  ::   ios                 ! Local integer output status for namelist read
   INTEGER  ::   iminspon
   NAMELIST/namagrif/ nn_cln_update, rn_sponge_tra, rn_sponge_dyn, ln_spc_dyn, ln_chk_bathy
   !!--------------------------------------------------------------------------------------
   !
   REWIND( numnam_ref )              ! Namelist namagrif in reference namelist : AGRIF zoom
   READ  ( numnam_ref, namagrif, IOSTAT = ios, ERR = 901)
901 IF( ios /= 0 ) CALL ctl_nam ( ios , 'namagrif in reference namelist', lwp )

   REWIND( numnam_cfg )              ! Namelist namagrif in configuration namelist : AGRIF zoom
   READ  ( numnam_cfg, namagrif, IOSTAT = ios, ERR = 902 )
902 IF( ios /= 0 ) CALL ctl_nam ( ios , 'namagrif in configuration namelist', lwp )
   IF(lwm) WRITE ( numond, namagrif )
   !
   IF(lwp) THEN                    ! control print
      WRITE(numout,*)
      WRITE(numout,*) 'agrif_nemo_init : AGRIF parameters'
      WRITE(numout,*) '~~~~~~~~~~~~~~~'
      WRITE(numout,*) '   Namelist namagrif : set AGRIF parameters'
      WRITE(numout,*) '      baroclinic update frequency       nn_cln_update = ', nn_cln_update
      WRITE(numout,*) '      sponge coefficient for tracers    rn_sponge_tra = ', rn_sponge_tra, ' s'
      WRITE(numout,*) '      sponge coefficient for dynamics   rn_sponge_tra = ', rn_sponge_dyn, ' s'
      WRITE(numout,*) '      use special values for dynamics   ln_spc_dyn    = ', ln_spc_dyn
      WRITE(numout,*) '      check bathymetry                  ln_chk_bathy  = ', ln_chk_bathy
      WRITE(numout,*) 
   ENDIF
   !
   ! convert DOCTOR namelist name into OLD names
   nbclineupdate = nn_cln_update
   visc_tra      = rn_sponge_tra
   visc_dyn      = rn_sponge_dyn
   !
   ! Check sponge length:
   iminspon = MIN(FLOOR(REAL(jpiglo-4)/REAL(2*Agrif_irhox())), FLOOR(REAL(jpjglo-4)/REAL(2*Agrif_irhox())) )
   IF (lk_mpp) iminspon = MIN(iminspon,FLOOR(REAL(jpi-2)/REAL(Agrif_irhox())), FLOOR(REAL(jpj-2)/REAL(Agrif_irhox())) )
   IF (nn_sponge_len > iminspon)  CALL ctl_stop('agrif sponge length is too large')
   !
   IF( agrif_oce_alloc()  > 0 )   CALL ctl_warn('agrif agrif_oce_alloc: allocation of arrays failed')
# if defined key_lim2
   IF( agrif_ice_alloc()  > 0 )   CALL ctl_stop('agrif agrif_ice_alloc: allocation of arrays failed')
# endif
   !
END SUBROUTINE agrif_nemo_init

# if defined key_mpp_mpi

SUBROUTINE Agrif_InvLoc( indloc, nprocloc, i, indglob )
   !!----------------------------------------------------------------------
   !!                     *** ROUTINE Agrif_detect ***
   !!----------------------------------------------------------------------
   USE dom_oce
   IMPLICIT NONE
   !
   INTEGER :: indglob, indloc, nprocloc, i
   !!----------------------------------------------------------------------
   !
   SELECT CASE( i )
   CASE(1)   ;   indglob = indloc + nimppt(nprocloc+1) - 1
   CASE(2)   ;   indglob = indloc + njmppt(nprocloc+1) - 1
   CASE DEFAULT
      indglob = indloc
   END SELECT
   !
END SUBROUTINE Agrif_InvLoc

SUBROUTINE Agrif_get_proc_info( imin, imax, jmin, jmax )
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_get_proc_info ***
   !!----------------------------------------------------------------------
   USE par_oce
   IMPLICIT NONE
   !
   INTEGER, INTENT(out) :: imin, imax
   INTEGER, INTENT(out) :: jmin, jmax
   !!----------------------------------------------------------------------
   !
   imin = nimppt(Agrif_Procrank+1)  ! ?????
   jmin = njmppt(Agrif_Procrank+1)  ! ?????
   imax = imin + jpi - 1
   jmax = jmin + jpj - 1
   ! 
END SUBROUTINE Agrif_get_proc_info

SUBROUTINE Agrif_estimate_parallel_cost(imin, imax,jmin, jmax, nbprocs, grid_cost)
   !!----------------------------------------------------------------------
   !!                 *** ROUTINE Agrif_estimate_parallel_cost ***
   !!----------------------------------------------------------------------
   USE par_oce
   IMPLICIT NONE
   !
   INTEGER,  INTENT(in)  :: imin, imax
   INTEGER,  INTENT(in)  :: jmin, jmax
   INTEGER,  INTENT(in)  :: nbprocs
   REAL(wp), INTENT(out) :: grid_cost
   !!----------------------------------------------------------------------
   !
   grid_cost = REAL(imax-imin+1,wp)*REAL(jmax-jmin+1,wp) / REAL(nbprocs,wp)
   !
END SUBROUTINE Agrif_estimate_parallel_cost

# endif

#else
SUBROUTINE Subcalledbyagrif
   !!----------------------------------------------------------------------
   !!                   *** ROUTINE Subcalledbyagrif ***
   !!----------------------------------------------------------------------
   WRITE(*,*) 'Impossible to be here'
END SUBROUTINE Subcalledbyagrif
#endif
