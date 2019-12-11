MODULE trd_oce
   !!======================================================================
   !!                   ***  MODULE trd_oce  ***
   !! Ocean trends :   set tracer and momentum trend variables
   !!======================================================================
   !! History :  1.0  !  2004-08  (C. Talandier) Original code
   !!----------------------------------------------------------------------
   USE par_oce                 ! ocean parameters
   USE trdmxl_oce              ! ocean active mixed layer tracers trends variables
   USE trdvor_oce              ! ocean vorticity trends variables

   IMPLICIT NONE
   PUBLIC

   !                                                   !!* Namelist namtrd:  diagnostics on dynamics/tracer trends *
   LOGICAL , PUBLIC  ::   ln_dyn_trd   = .FALSE.        !: (T) 3D momentum             trends or (F) not
   LOGICAL , PUBLIC  ::   ln_tra_trd   = .FALSE.        !: (T) 3D tracer               trends or (F) not
   LOGICAL , PUBLIC  ::   ln_KE_trd    = .FALSE.        !: (T) 3D Kinetic   Energy     trends or (F) not
   LOGICAL , PUBLIC  ::   ln_PE_trd    = .FALSE.        !: (T) 3D Potential Energy     trends or (F) not
   LOGICAL , PUBLIC  ::   ln_vor_trd   = .FALSE.        !: (T) 3D barotropic vorticity trends or (F) not
   LOGICAL , PUBLIC  ::   ln_glo_trd   = .FALSE.        !: (T) global domain averaged diag for T, T^2, KE, and PE
   LOGICAL , PUBLIC  ::   ln_dyn_mxl   = .FALSE.        !: (T) 2D tracer   trends averaged over the mixed layer 
   LOGICAL , PUBLIC  ::   ln_tra_mxl   = .FALSE.        !: (T) 2D momentum trends averaged over the mixed layer 
   INTEGER , PUBLIC  ::   nn_trd       = 10             !: time step frequency for ln_glo_trd=T only

   LOGICAL , PUBLIC ::   l_trdtra        !: tracers  trend flag (set from namelist in trdini)
   LOGICAL , PUBLIC ::   l_trddyn        !: momentum trend flag (set from namelist in trdini)
   
# if ( defined key_trdtrc && defined key_iomput )  ||  defined key_trdmxl_trc
   LOGICAL , PUBLIC ::   l_trdtrc = .TRUE.        !: tracers  trend flag
# else
   LOGICAL , PUBLIC ::   l_trdtrc = .FALSE.       !: tracers  trend flag
# endif
   !                                                  !!!* Active tracers trends indexes
   INTEGER, PUBLIC, PARAMETER ::   jptot_tra  = 20     !: Total trend nb: change it when adding/removing one indice below
   !                               ===============     !  
   INTEGER, PUBLIC, PARAMETER ::   jptra_xad  =  1     !: x- horizontal advection
   INTEGER, PUBLIC, PARAMETER ::   jptra_yad  =  2     !: y- horizontal advection
   INTEGER, PUBLIC, PARAMETER ::   jptra_zad  =  3     !: z- vertical   advection
   INTEGER, PUBLIC, PARAMETER ::   jptra_sad  =  4     !: z- vertical   advection
   INTEGER, PUBLIC, PARAMETER ::   jptra_totad  =  5   !: total         advection
   INTEGER, PUBLIC, PARAMETER ::   jptra_ldf  =  6     !: lateral       diffusion
   INTEGER, PUBLIC, PARAMETER ::   jptra_iso_x  =  7   !: x-component of isopycnal diffusion
   INTEGER, PUBLIC, PARAMETER ::   jptra_iso_y  =  8   !: y-component of isopycnal diffusion
   INTEGER, PUBLIC, PARAMETER ::   jptra_iso_z1 =  9   !: z-component of isopycnal diffusion
   INTEGER, PUBLIC, PARAMETER ::   jptra_zdf  = 10     !: vertical      diffusion
   INTEGER, PUBLIC, PARAMETER ::   jptra_zdfp = 11     !: "PURE" vert.  diffusion (ln_traldf_iso=T)
   INTEGER, PUBLIC, PARAMETER ::   jptra_evd  = 12     !: EVD term (convection)
   INTEGER, PUBLIC, PARAMETER ::   jptra_bbc  = 13     !: Bottom Boundary Condition (geoth. heating) 
   INTEGER, PUBLIC, PARAMETER ::   jptra_bbl  = 14     !: Bottom Boundary Layer (diffusive and/or advective)
   INTEGER, PUBLIC, PARAMETER ::   jptra_npc  = 15     !: non-penetrative convection treatment
   INTEGER, PUBLIC, PARAMETER ::   jptra_dmp  = 16     !: internal restoring (damping)
   INTEGER, PUBLIC, PARAMETER ::   jptra_qsr  = 17     !: penetrative solar radiation
   INTEGER, PUBLIC, PARAMETER ::   jptra_nsr  = 18     !: non solar radiation / C/D on salinity  (+runoff if ln_rnf=T)
   INTEGER, PUBLIC, PARAMETER ::   jptra_atf  = 19     !: Asselin time filter
   INTEGER, PUBLIC, PARAMETER ::   jptra_tot  = 20     !: Model total trend
   !
   !                                                  !!!* Passive tracers trends indices (use if "key_top" defined)
   INTEGER, PUBLIC, PARAMETER ::   jptra_sms  = 19     !: sources m. sinks
   INTEGER, PUBLIC, PARAMETER ::   jptra_radn = 20     !: corr. trn<0 in trcrad
   INTEGER, PUBLIC, PARAMETER ::   jptra_radb = 21     !: corr. trb<0 in trcrad (like atf)
   !
   !                                                  !!!* Momentum trends indices
   INTEGER, PUBLIC, PARAMETER ::   jptot_dyn  = 16     !: Total trend nb: change it when adding/removing one indice below
   !                               ===============     !  
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_hpg  =  1     !: hydrostatic pressure gradient 
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_spg  =  2     !: surface     pressure gradient
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_keg  =  3     !: kinetic energy gradient  or horizontal advection
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_rvo  =  4     !: relative  vorticity      or metric term
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_pvo  =  5     !: planetary vorticity
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_zad  =  6     !: vertical advection
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_ldf  =  7     !: horizontal diffusion   
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_zdf  =  8     !: vertical   diffusion
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_bfr  =  9     !: bottom  stress 
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_atf  = 10     !: Asselin time filter
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_tau  = 11     !: surface stress
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_bfri = 12     !: implicit bottom friction (ln_bfrimp=.TRUE.)
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_ken  = 13     !: use for calculation of KE
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_spgflt  = 14  !: filter contribution to surface pressure gradient (spg_flt)
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_spgexp  = 15  !: explicit contribution to surface pressure gradient (spg_flt)
   INTEGER, PUBLIC, PARAMETER ::   jpdyn_eivke   = 16  !: K.E trend from Gent McWilliams scheme
   !
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trd_oce
