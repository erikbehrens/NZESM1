MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!----------------------------------------------------------------------
#if defined key_pisces || defined key_pisces_reduced 
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                         PISCES model
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   INTEGER ::   numnatp_ref = -1           !! Logical units for namelist pisces
   INTEGER ::   numnatp_cfg = -1           !! Logical units for namelist pisces
   INTEGER ::   numonp      = -1           !! Logical unit for namelist pisces output

   !!*  Biological fluxes for light : variables shared by pisces & lobster
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  neln  !: number of T-levels + 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup  !: euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot  !: par (photosynthetic available radiation)
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  xksi  !:  LOBSTER : zooplakton closure
   !                                                       !:  PISCES  : silicon dependant half saturation

#if defined key_pisces 
   !!*  Time variables
   INTEGER  ::   nrdttrc           !: ???
   INTEGER  ::   ndayflxtr         !: ???
   REAL(wp) ::   rfact , rfactr    !: ???
   REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstep             !: Time step duration for biology
   REAL(wp) ::   ryyss             !: number of seconds per year 
   REAL(wp) ::   r1_ryyss          !: inverse number of seconds per year 


   !!*  Biological parameters 
   INTEGER  ::   niter1max, niter2max   !: Maximum number of iterations for sinking
   REAL(wp) ::   rno3              !: ???
   REAL(wp) ::   o2ut              !: ???
   REAL(wp) ::   po4r              !: ???
   REAL(wp) ::   rdenit            !: ???
   REAL(wp) ::   rdenita           !: ???
   REAL(wp) ::   o2nit             !: ???
   REAL(wp) ::   wsbio, wsbio2     !: ???
   REAL(wp) ::   xkmort            !: ???
   REAL(wp) ::   ferat3            !: ???

   !!*  diagnostic parameters 
   REAL(wp) ::  tpp                !: total primary production
   REAL(wp) ::  t_oce_co2_exp      !: total carbon export
   REAL(wp) ::  t_oce_co2_flx      !: Total ocean carbon flux
   REAL(wp) ::  t_oce_co2_flx_cum  !: Cumulative Total ocean carbon flux
   REAL(wp) ::  t_atm_co2_flx      !: global mean of atmospheric pco2

   !!* restoring
   LOGICAL  ::  ln_pisdmp          !: restoring or not of nutrients to a mean value
   INTEGER  ::  nn_pisdmp          !: frequency of relaxation or not of nutrients to a mean value

   !!* Mass conservation
   LOGICAL  ::  ln_check_mass      !: Flag to check mass conservation

   !!*  Biological fluxes for primary production
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xksimax    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanono3   !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatno3   !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanonh4   !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatnh4   !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanopo4   !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatpo4   !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphy    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdia    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdfe    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnfe    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnfe    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdfe    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimsi     !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   biron      !: bioavailable fraction of iron


   !!*  SMS for the organic matter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfracal    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xlimbac    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xlimbacl   !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xdiss      !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal    !: Calcite production

   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akb3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akw3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   borat      !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aphscale   !: 


   !!* Temperature dependancy of SMS terms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp. dependancy of various biological rates
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc2   !: Temp. dependancy of mesozooplankton rates

#if defined key_kriest
   !!*  Kriest parameter for aggregation
   REAL(wp) ::   xkr_eta                            !: Sinking  exponent 
   REAL(wp) ::   xkr_zeta                           !:  N content exponent 
   REAL(wp) ::   xkr_ncontent                       !:  N content factor   
   REAL(wp) ::   xkr_massp                          !: 
   REAL(wp) ::   xkr_mass_min, xkr_mass_max         !:  Minimum, Maximum mass for Aggregates 
#endif

#endif
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_pisces_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      INTEGER ::   ierr(5)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light : shared variables for pisces & lobster
      ALLOCATE( etot(jpi,jpj,jpk), neln(jpi,jpj), heup(jpi,jpj), xksi(jpi,jpj), STAT=ierr(1) )
      !
#if defined key_pisces
      !*  Biological fluxes for primary production
      ALLOCATE( xksimax(jpi,jpj)     , biron   (jpi,jpj,jpk),       &
         &      xnanono3(jpi,jpj,jpk), xdiatno3(jpi,jpj,jpk),       &
         &      xnanonh4(jpi,jpj,jpk), xdiatnh4(jpi,jpj,jpk),       &
         &      xnanopo4(jpi,jpj,jpk), xdiatpo4(jpi,jpj,jpk),       &
         &      xlimphy (jpi,jpj,jpk), xlimdia (jpi,jpj,jpk),       &
         &      xlimnfe (jpi,jpj,jpk), xlimdfe (jpi,jpj,jpk),       &
         &      xlimsi  (jpi,jpj,jpk), concdfe (jpi,jpj,jpk),       &
         &      concnfe (jpi,jpj,jpk),                           STAT=ierr(2) ) 
         !
      !*  SMS for the organic matter
      ALLOCATE( xfracal (jpi,jpj,jpk), nitrfac(jpi,jpj,jpk),       &
         &      xlimbac (jpi,jpj,jpk), xdiss  (jpi,jpj,jpk),       & 
         &      xlimbacl(jpi,jpj,jpk), prodcal(jpi,jpj,jpk),     STAT=ierr(3) )

      !* Variable for chemistry of the CO2 cycle
      ALLOCATE( akb3(jpi,jpj,jpk)    , ak13  (jpi,jpj,jpk) ,       &
         &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,       &
         &      akw3(jpi,jpj,jpk)    , borat (jpi,jpj,jpk) ,       &
         &      hi  (jpi,jpj,jpk)    , excess(jpi,jpj,jpk) ,       &
         &      aphscale(jpi,jpj,jpk),                           STAT=ierr(4) )
         !
      !* Temperature dependancy of SMS terms
      ALLOCATE( tgfunc(jpi,jpj,jpk)  , tgfunc2(jpi,jpj,jpk) ,    STAT=ierr(5) )
         !
#endif
      !
      sms_pisces_alloc = MAXVAL( ierr )
      !
      IF( sms_pisces_alloc /= 0 )   CALL ctl_warn('sms_pisces_alloc: failed to allocate arrays') 
      !
   END FUNCTION sms_pisces_alloc

#else
   !!----------------------------------------------------------------------   
   !!  Empty module :                                     NO PISCES model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================   
END MODULE sms_pisces    
