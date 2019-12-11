MODULE biaspar
   !! Variables relevant to bias module
   !!======================================================================
   !!                 *** Module biaspar ***
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY: &
      & wp
   USE par_oce, ONLY: &
      & jpi, &
      & jpj, &
      & jpk

   IMPLICIT NONE
   PUBLIC

   !! * Shared module variables
   LOGICAL, PUBLIC :: ln_bias        = .FALSE. !: estimate (apply) bias arrays 
   LOGICAL, PUBLIC :: ln_bias_asm    = .FALSE. !: estimate bias from assim incr
   LOGICAL, PUBLIC :: ln_bias_rlx    = .FALSE. !: estimate bias from relaxation 
   LOGICAL, PUBLIC :: ln_bias_ofl    = .FALSE. !: bias estimated offline
   LOGICAL, PUBLIC :: ln_bias_ts_app = .FALSE. !: estimate (apply) bias arrays 
   LOGICAL, PUBLIC :: ln_bias_pc_app = .FALSE. !: estimate bias from assim incr
   LOGICAL, PUBLIC :: lrst_bias      = .FALSE. !: estimate bias from assim incr

   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: &
      & tbias, &       !: Temperature bias field for T correction
      & tbias_p, &     !:  "           "    "    "   P correction
      & tbias_i, &     !:  "           "    "    "   incremental P correction
      & sbias, &       !: Salinity bias field    for S correction
      & sbias_p, &     !:  "           "    "    "   P correction
      & sbias_i, &     !:  "           "    "    "   incremental P correction
      & rhd_pc         !: Press corrtd density from online to use in dyn_hpg


   REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE :: &
      & gru_pc, &      !: Press corrtd bottom pressure gradient (x-dir)
      & grv_pc         !: Press corrtd bottom pressure gradient (y-dir)



END MODULE biaspar
