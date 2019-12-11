MODULE oce_trc
   !!======================================================================
   !!                      ***  MODULE  oce_trc  ***
   !! TOP :   variables shared between ocean and passive tracers
   !!======================================================================
   !! History :   1.0  !  2004-03  (C. Ethe)  original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  rewritting
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------

   !* Domain size *
   USE par_oce , ONLY :   jpi      =>   jpi        !: first  dimension of grid --> i 
   USE par_oce , ONLY :   jpj      =>   jpj        !: second dimension of grid --> j  
   USE par_oce , ONLY :   jpk      =>   jpk        !: number of levels  
   USE par_oce , ONLY :   jpim1    =>   jpim1      !: jpi - 1
   USE par_oce , ONLY :   jpjm1    =>   jpjm1      !: jpj - 1 
   USE par_oce , ONLY :   jpkm1    =>   jpkm1      !: jpk - 1  
   USE par_oce , ONLY :   jpij     =>   jpij       !: jpi x jpj
   USE par_oce , ONLY :   lk_esopa =>   lk_esopa   !: flag to activate the all option
   USE par_oce , ONLY :   jp_tem   =>   jp_tem     !: indice for temperature
   USE par_oce , ONLY :   jp_sal   =>   jp_sal     !: indice for salinity

   !* IO manager *
   USE in_out_manager    
 
   !* Memory Allocation *
   USE wrk_nemo      
 
   !* Timing *
   USE timing    
 
   !* MPP library                         
   USE lib_mpp 

   !* Fortran utilities                         
   USE lib_fortran

   !* Lateral boundary conditions                         
   USE lbclnk

   !* physical constants *
   USE phycst            

   !* 1D configuration
   USE c1d                                         

   !* model domain *
   USE dom_oce 

   USE domvvl, ONLY : un_td, vn_td          !: thickness diffusion transport
   USE domvvl, ONLY : ln_vvl_ztilde         !: ztilde vertical coordinate
   USE domvvl, ONLY : ln_vvl_layer          !: level  vertical coordinate

   !* ocean fields: here now and after fields *
   USE oce , ONLY :   ua      =>    ua      !: i-horizontal velocity (m s-1) 
   USE oce , ONLY :   va      =>    va      !: j-horizontal velocity (m s-1)
   USE oce , ONLY :   un      =>    un      !: i-horizontal velocity (m s-1) 
   USE oce , ONLY :   vn      =>    vn      !: j-horizontal velocity (m s-1)
   USE oce , ONLY :   wn      =>    wn      !: vertical velocity (m s-1)  
   USE oce , ONLY :   tsn     =>    tsn     !: 4D array contaning ( tn, sn )
   USE oce , ONLY :   tsb     =>    tsb     !: 4D array contaning ( tb, sb )
   USE oce , ONLY :   tsa     =>    tsa     !: 4D array contaning ( ta, sa )
   USE oce , ONLY :   rhop    =>    rhop    !: potential volumic mass (kg m-3) 
   USE oce , ONLY :   rhd     =>    rhd     !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)
#if defined key_offline
   USE oce , ONLY :   rab_n   =>    rab_n   !: local thermal/haline expension ratio at T-points
#endif
   USE oce , ONLY :   hdivn   =>    hdivn   !: horizontal divergence (1/s)
   USE oce , ONLY :   rotn    =>    rotn    !: relative vorticity    [s-1]
   USE oce , ONLY :   hdivb   =>    hdivb   !: horizontal divergence (1/s)
   USE oce , ONLY :   rotb    =>    rotb    !: relative vorticity    [s-1]
   USE oce , ONLY :   sshn    =>    sshn    !: sea surface height at t-point [m]   
   USE oce , ONLY :   sshb    =>    sshb    !: sea surface height at t-point [m]   
   USE oce , ONLY :   ssha    =>    ssha    !: sea surface height at t-point [m]   
   USE oce , ONLY :   l_traldf_rot => l_traldf_rot  !: rotated laplacian operator for lateral diffusion

   !* surface fluxes *
   USE sbc_oce , ONLY :   utau       =>    utau       !: i-surface stress component
   USE sbc_oce , ONLY :   vtau       =>    vtau       !: j-surface stress component
   USE sbc_oce , ONLY :   wndm       =>    wndm       !: 10m wind speed 
   USE sbc_oce , ONLY :   qsr        =>    qsr        !: penetrative solar radiation (w m-2)
   USE sbc_oce , ONLY :   emp        =>    emp        !: freshwater budget: volume flux               [Kg/m2/s]
   USE sbc_oce , ONLY :   emp_b      =>    emp_b      !: freshwater budget: volume flux               [Kg/m2/s]
   USE sbc_oce , ONLY :   fmmflx     =>    fmmflx     !: freshwater budget: volume flux               [Kg/m2/s]
   USE sbc_oce , ONLY :   rnf        =>    rnf        !: river runoff   [Kg/m2/s]
   USE sbc_oce , ONLY :   ln_dm2dc   =>    ln_dm2dc   !: Diurnal Cycle 
   USE sbc_oce , ONLY :   ncpl_qsr_freq   =>   ncpl_qsr_freq   !: qsr coupling frequency per days from atmospher
   USE sbc_oce , ONLY :   ln_rnf     =>    ln_rnf     !: runoffs / runoff mouths
   USE sbc_oce , ONLY :   fr_i       =>    fr_i       !: ice fraction (between 0 to 1)
   USE sbc_oce , ONLY :   nn_ice_embd => nn_ice_embd  !: flag for  levitating/embedding sea-ice in the ocean
   USE traqsr  , ONLY :   rn_abs     =>    rn_abs     !: fraction absorbed in the very near surface
   USE traqsr  , ONLY :   rn_si0     =>    rn_si0     !: very near surface depth of extinction
   USE traqsr  , ONLY :   ln_qsr_bio =>    ln_qsr_bio !: flag to use or not the biological fluxes for light
   USE sbcrnf  , ONLY :   rnfmsk     =>    rnfmsk     !: mixed adv scheme in runoffs vicinity (hori.) 
   USE sbcrnf  , ONLY :   rnfmsk_z   =>    rnfmsk_z   !: mixed adv scheme in runoffs vicinity (vert.)
   USE sbcrnf  , ONLY :   h_rnf      =>    h_rnf      !: river runoff   [Kg/m2/s]
   USE sbcrnf  , ONLY :   nk_rnf     =>    nk_rnf     !: depth of runoff in model level

   USE trc_oce

   !* lateral diffusivity (tracers) *
   USE ldftra_oce , ONLY :  rldf     =>   rldf        !: multiplicative coef. for lateral diffusivity
   USE ldftra_oce , ONLY :  rn_aht_0 =>   rn_aht_0    !: horizontal eddy diffusivity for tracers (m2/s)
   USE ldftra_oce , ONLY :  aht0     =>   aht0        !: horizontal eddy diffusivity for tracers (m2/s)
   USE ldftra_oce , ONLY :  ahtb0    =>   ahtb0       !: background eddy diffusivity for isopycnal diff. (m2/s)
   USE ldftra_oce , ONLY :  ahtu     =>   ahtu        !: lateral diffusivity coef. at u-points 
   USE ldftra_oce , ONLY :  ahtv     =>   ahtv        !: lateral diffusivity coef. at v-points 
   USE ldftra_oce , ONLY :  ahtw     =>   ahtw        !: lateral diffusivity coef. at w-points 
   USE ldftra_oce , ONLY :  ahtt     =>   ahtt        !: lateral diffusivity coef. at t-points
   USE ldftra_oce , ONLY :  aeiv0    =>   aeiv0       !: eddy induced velocity coefficient (m2/s) 
   USE ldftra_oce , ONLY :  aeiu     =>   aeiu        !: eddy induced velocity coef. at u-points (m2/s)   
   USE ldftra_oce , ONLY :  aeiv     =>   aeiv        !: eddy induced velocity coef. at v-points (m2/s) 
   USE ldftra_oce , ONLY :  aeiw     =>   aeiw        !: eddy induced velocity coef. at w-points (m2/s) 
   USE ldftra_oce , ONLY :  lk_traldf_eiv  =>  lk_traldf_eiv     !: eddy induced velocity flag
   USE ldftra_oce , ONLY :  r_fact_lap     =>  r_fact_lap        !: enhanced zonal diffusivity coefficient

   !* vertical diffusion *
   USE zdf_oce , ONLY :   avt        =>   avt         !: vert. diffusivity coef. at w-point for temp  
# if defined key_zdfddm
   USE zdfddm  , ONLY :   avs        =>   avs         !: salinity vertical diffusivity coeff. at w-point
# endif

   !* mixing & mixed layer depth *
   USE zdfmxl , ONLY :   nmln        =>   nmln        !: number of level in the mixed layer
   USE zdfmxl , ONLY :   hmld        =>   hmld        !: mixing layer depth (turbocline)
   USE zdfmxl , ONLY :   hmlp        =>   hmlp        !: mixed layer depth  (rho=rho0+zdcrit) (m)
   USE zdfmxl , ONLY :   hmlpt       =>   hmlpt       !: mixed layer depth at t-points (m)

   !* direction of lateral diffusion *
   USE ldfslp , ONLY :   lk_ldfslp  =>  lk_ldfslp     !: slopes flag
# if   defined key_ldfslp
   USE ldfslp , ONLY :   uslp       =>   uslp         !: i-direction slope at u-, w-points
   USE ldfslp , ONLY :   vslp       =>   vslp         !: j-direction slope at v-, w-points
   USE ldfslp , ONLY :   wslpi      =>   wslpi        !: i-direction slope at u-, w-points
   USE ldfslp , ONLY :   wslpj      =>   wslpj        !: j-direction slope at v-, w-points
# endif

   USE diaar5 , ONLY :   lk_diaar5  =>   lk_diaar5
#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE oce_trc
