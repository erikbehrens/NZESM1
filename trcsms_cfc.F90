MODULE trcsms_cfc
   !!======================================================================
   !!                      ***  MODULE trcsms_cfc  ***
   !! TOP : CFC main model
   !!======================================================================
   !! History :  OPA  !  1999-10  (JC. Dutay)  original code
   !!  NEMO      1.0  !  2004-03  (C. Ethe) free form + modularity
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  reorganisation
   !!                 !  2016-06  (J. Palmieri)  update for UKESM1
   !!                 !  2017-04  (A. Yool)  update to add SF6, fix coefficients
   !!----------------------------------------------------------------------
#if defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !!   trc_sms_cfc  :  compute and add CFC suface forcing to CFC trends
   !!   cfc_init     :  sets constants for CFC surface forcing computation
   !!----------------------------------------------------------------------
   USE dom_oce       ! ocean space and time domain
   USE oce_trc       ! Ocean variables
   USE par_trc       ! TOP parameters
   USE trc           ! TOP variables
   USE trd_oce
   USE trdtrc
   USE iom           ! I/O library
  !CEB
  !!     USE wrk_nemo        ! work arrays
     !/CEB
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_cfc         ! called in ???    
   PUBLIC   trc_sms_cfc_alloc   ! called in trcini_cfc.F90

   INTEGER , PUBLIC, PARAMETER ::   jphem  =   2   ! parameter for the 2 hemispheres
   INTEGER , PUBLIC            ::   jpyear         ! Number of years read in CFC1112 file
   INTEGER , PUBLIC            ::   ndate_beg      ! initial calendar date (aammjj) for CFC
   INTEGER , PUBLIC            ::   simu_type      ! Kind of simulation: 1- Spin-up 
                                                   !                     2- Hindcast/projection
   INTEGER , PUBLIC            ::   nyear_res      ! restoring time constant (year)
   INTEGER , PUBLIC            ::   nyear_beg      ! initial year (aa) 
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   p_cfc    ! partial hemispheric pressure for CFC
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   xphem    ! spatial interpolation factor for patm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qtr_cfc  ! flux at surface
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qint_cfc ! cumulative flux 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   patm     ! atmospheric function

   REAL(wp), DIMENSION(4,3) ::   soa   ! coefficient for solubility of CFC [mol/l/atm]
   REAL(wp), DIMENSION(3,3) ::   sob   !    "               "
   REAL(wp), DIMENSION(5,3) ::   sca   ! coefficients for schmidt number in degre Celcius
      
   !                          ! coefficients for conversion
   REAL(wp) ::   xconv1 = 1.0          ! conversion from to 
   REAL(wp) ::   xconv2 = 0.01/3600.   ! conversion from cm/h to m/s: 
   REAL(wp) ::   xconv3 = 1.0e+3       ! conversion from mol/l/atm to mol/m3/atm
   REAL(wp) ::   xconv4 = 1.0e-12      ! conversion from mol/m3/atm to mol/m3/pptv 

   !! trend temporary array:
   !CEB   REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrcfc
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrcfc
   !/CEB
   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_cfc( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_cfc  ***
      !!
      !! ** Purpose :   Compute the surface boundary contition on CFC 11 
      !!             passive tracer associated with air-mer fluxes and add it 
      !!             to the general trend of tracers equations.
      !!
      !! ** Method  : - get the atmospheric partial pressure - given in pico -
      !!              - computation of solubility ( in 1.e-12 mol/l then in 1.e-9 mol/m3)
      !!              - computation of transfert speed ( given in cm/hour ----> cm/s )
      !!              - the input function is given by : 
      !!                speed * ( concentration at equilibrium - concentration at surface )
      !!              - the input function is in pico-mol/m3/s and the
      !!                CFC concentration in pico-mol/m3
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jn, jl, jm, js
      INTEGER  ::   iyear_beg, iyear_end, iyear_tmp
      INTEGER  ::   im1, im2, ierr
      REAL(wp) ::   ztap, zdtap        
      REAL(wp) ::   zt1, zt2, zt3, zt4, zv2
      REAL(wp) ::   zsol      ! solubility
      REAL(wp) ::   zsch      ! schmidt number 
      REAL(wp) ::   zpp_cfc   ! atmospheric partial pressure of CFC
      REAL(wp) ::   zca_cfc   ! concentration at equilibrium
      REAL(wp) ::   zak_cfc   ! transfert coefficients
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)  ::   zpatm     ! atmospheric function
      !!----------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_cfc')
      !
      ALLOCATE( zpatm(jphem,jp_cfc), STAT=ierr )
      IF( ierr > 0 ) THEN
         CALL ctl_stop( 'trc_sms_cfc: unable to allocate zpatm array' )   ;   RETURN
      ENDIF

      IF( kt == nittrc000 )   CALL cfc_init

      ! Temporal interpolation
      ! ----------------------
      !! JPALM -- 15-06-2016 -- define 2 kinds of CFC run:
      !!                     1- the SPIN-UP and 2- Hindcast/Projections
      !!                     -- main difference is the way to define the year of
      !!                     simulation, that determine the atm pCFC.
      !!                     1-- Spin-up: our atm forcing is of 30y we cycle on.
      !!                     So we do 90y CFC cycles to be in good
      !!                     correspondence with the atmosphere
      !!                     2-- Hindcast/proj, instead of nyear-1900 we keep
      !!                     the 2 last digit, and enable 3 cycle from 1800 to 2100.  
      !!----------------------------------------------------------------------
      IF (simu_type==1) THEN
         !! 1 -- SPIN-UP
         iyear_tmp = nyear - nyear_res  !! JPALM -- in our spin-up, nyear_res is 1000
         iyear_beg = MOD( iyear_tmp , 90 )
         !! JPALM -- the pCFC file only got 78 years.
         !!       So if iyear_beg > 78 then we set pCFC to 0
         !!             iyear_beg = 0 as well -- must try to avoid obvious problems
         !!             as Pcfc is set to 0.00 up to year 32, let set iyear_beg to year 10
         !!          else, must add 30 to iyear_beg to match with P_cfc indices
         !!---------------------------------------
         IF ((iyear_beg > 77) .OR. (iyear_beg==0)) THEN
            iyear_beg = 10
         ELSE 
            iyear_beg = iyear_beg + 30
         ENDIF
      ELSEIF (simu_type==2) THEN
         !! 2 -- Hindcast/proj
         iyear_beg = MOD(nyear, 100)
         IF (iyear_beg < 20)  iyear_beg = iyear_beg + 100
         !! JPALM -- Same than previously, if iyear_beg is out of P_cfc range,
         !!       we want to set p_CFC to 0.00 --> set iyear_beg = 10
         IF ((iyear_beg < 30) .OR. (iyear_beg > 115)) iyear_beg = 10             
      ENDIF
      !!
      IF ( nmonth <= 6 ) THEN
         iyear_beg = iyear_beg - 1
         im1       =  6 - nmonth + 1
         im2       =  6 + nmonth - 1
      ELSE
         im1       = 12 - nmonth + 7
         im2       =      nmonth - 7
      ENDIF
      iyear_end = iyear_beg + 1

      !                                                  !------------!
      DO jl = 1, jp_cfc                                  !  CFC loop  !
         !                                               !------------!
         jn = jp_cfc0 + jl - 1
         ! time interpolation at time kt
         DO jm = 1, jphem
            zpatm(jm,jl) = (  p_cfc(iyear_beg, jm, jl) * FLOAT (im1)  &
               &           +  p_cfc(iyear_end, jm, jl) * FLOAT (im2) ) / 12.
         END DO
         
         !                                                         !------------!
         DO jj = 1, jpj                                            !  i-j loop  !
            DO ji = 1, jpi                                         !------------!
 
               ! space interpolation
               zpp_cfc  =       xphem(ji,jj)   * zpatm(1,jl)   &
                  &     + ( 1.- xphem(ji,jj) ) * zpatm(2,jl)

               ! Computation of concentration at equilibrium : in picomol/l
               ! coefficient for solubility for CFC-11/12 in  mol/l/atm
               IF( tmask(ji,jj,1) .GE. 0.5 ) THEN
                  ztap  = ( tsn(ji,jj,1,jp_tem) + 273.16 ) * 0.01
                  zdtap = sob(1,jl) + ztap * ( sob(2,jl) + ztap * sob(3,jl) ) 
                  zsol  =  EXP( soa(1,jl) + soa(2,jl) / ztap + soa(3,jl) * LOG( ztap )   &
                     &                    + soa(4,jl) * ztap * ztap + tsn(ji,jj,1,jp_sal) * zdtap ) 
               ELSE
                  zsol  = 0.e0
               ENDIF
               ! conversion from mol/l/atm to mol/m3/atm and from mol/m3/atm to mol/m3/pptv    
               zsol = xconv4 * xconv3 * zsol * tmask(ji,jj,1)  
               ! concentration at equilibrium
               zca_cfc = xconv1 * zpp_cfc * zsol * tmask(ji,jj,1)             
  
               ! Computation of speed transfert
               !    Schmidt number
               zt1  = tsn(ji,jj,1,jp_tem)
               zt2  = zt1 * zt1 
               zt3  = zt1 * zt2
               zt4  = zt1 * zt3
               zsch = sca(1,jl) + sca(2,jl) * zt1 + sca(3,jl) * zt2 + sca(4,jl) * zt3 + sca(5,jl) * zt4

               !    speed transfert : formulae of wanninkhof 1992
               zv2     = wndm(ji,jj) * wndm(ji,jj)
               zsch    = zsch / 660.
               ! AXY (25/04/17): OMIP protocol specifies lower Wanninkhof (2014) value
               ! zak_cfc = ( 0.39 * xconv2 * zv2 / SQRT(zsch) ) * tmask(ji,jj,1)
               zak_cfc = ( 0.251 * xconv2 * zv2 / SQRT(zsch) ) * tmask(ji,jj,1)

               ! Input function  : speed *( conc. at equil - concen at surface )
               ! trn in pico-mol/l idem qtr; ak in en m/a
               qtr_cfc(ji,jj,jl) = -zak_cfc * ( trb(ji,jj,1,jn) - zca_cfc )   &
#if defined key_degrad
                  &                         * facvol(ji,jj,1)                           &
#endif
                  &                         * tmask(ji,jj,1) * ( 1. - fr_i(ji,jj) )
               ! Add the surface flux to the trend
               tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + qtr_cfc(ji,jj,jl) / fse3t(ji,jj,1) 

               ! cumulation of surface flux at each time step
               qint_cfc(ji,jj,jl) = qint_cfc(ji,jj,jl) + qtr_cfc(ji,jj,jl) * rdt
               !                                               !----------------!
            END DO                                             !  end i-j loop  !
         END DO                                                !----------------!
         !                                                  !----------------!
      END DO                                                !  end CFC loop  !
         !
      IF( kt == nittrc000 ) THEN
         DO jl = 1, jp_cfc   
             WRITE(NUMOUT,*) ' '
             WRITE(NUMOUT,*) 'CFC interpolation verification '  !! Jpalm  
             WRITE(NUMOUT,*) '################################## '
             WRITE(NUMOUT,*) ' '
               if (jl.EQ.1) then
                   WRITE(NUMOUT,*) 'Traceur = CFC11: '
               elseif (jl.EQ.2) then
                   WRITE(NUMOUT,*) 'Traceur = CFC12: '
               elseif (jl.EQ.3) then
                   WRITE(NUMOUT,*) 'Traceur = SF6: '
               endif
             WRITE(NUMOUT,*) 'nyear    = ', nyear
             WRITE(NUMOUT,*) 'nmonth   = ', nmonth
             WRITE(NUMOUT,*) 'iyear_beg= ', iyear_beg
             WRITE(NUMOUT,*) 'iyear_end= ', iyear_end
             WRITE(NUMOUT,*) 'p_cfc(iyear_beg)= ',p_cfc(iyear_beg, 1, jl)
             WRITE(NUMOUT,*) 'p_cfc(iyear_end)= ',p_cfc(iyear_end, 1, jl)
             WRITE(NUMOUT,*) 'Im1= ',im1
             WRITE(NUMOUT,*) 'Im2= ',im2
             WRITE(NUMOUT,*) 'zpp_cfc = ',zpp_cfc
             WRITE(NUMOUT,*) ' '
         END DO  
# if defined key_debug_medusa
         CALL flush(numout)
# endif
      ENDIF
        !
      !IF( lrst_trc ) THEN
      !   IF(lwp) WRITE(numout,*)
      !   IF(lwp) WRITE(numout,*) 'trc_sms_cfc : cumulated input function fields written in ocean restart file ',   &
      !      &                    'at it= ', kt,' date= ', ndastp
      !   IF(lwp) WRITE(numout,*) '~~~~'
      !   DO jn = jp_cfc0, jp_cfc1
      !      CALL iom_rstput( kt, nitrst, numrtw, 'qint_'//ctrcnm(jn), qint_cfc(:,:,jn) )
      !   END DO
      !ENDIF                                            
      !
      IF  (iom_use("qtrCFC11"))  CALL iom_put( "qtrCFC11"  , qtr_cfc (:,:,1) )
      IF  (iom_use("qintCFC11")) CALL iom_put( "qintCFC11" , qint_cfc(:,:,1) )
      IF  (iom_use("qtrCFC12"))  CALL iom_put( "qtrCFC12"  , qtr_cfc (:,:,2) )
      IF  (iom_use("qintCFC12")) CALL iom_put( "qintCFC12" , qint_cfc(:,:,2) )
      IF  (iom_use("qtrSF6"))    CALL iom_put( "qtrSF6"    , qtr_cfc (:,:,3) )
      IF  (iom_use("qintSF6"))   CALL iom_put( "qintSF6"   , qint_cfc(:,:,3) )
      !
      IF( l_trdtrc ) THEN
          !CEB 
          !CALL wrk_alloc( jpi, jpj, jpk, ztrcfc )
          ALLOCATE(ztrcfc(jpi,jpj,jpk))
          !/CEB 
          DO jn = jp_cfc0, jp_cfc1
             ztrcfc(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrcfc, jn, jptra_sms, kt )   ! save trends
          END DO
          !CEB
          !CALL wrk_dealloc( jpi, jpj, jpk, ztrcfc )
          DEALLOCATE(ztrcfc)
          !/CEB
      END IF
      !
# if defined key_debug_medusa
      IF(lwp) WRITE(numout,*) '   CFC - Check: nn_timing = ', nn_timing
      CALL flush(numout)
# endif
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_cfc')
      !
   END SUBROUTINE trc_sms_cfc


   SUBROUTINE cfc_init
      !!---------------------------------------------------------------------
      !!                     ***  cfc_init  ***  
      !!
      !! ** Purpose : sets constants for CFC model
      !!---------------------------------------------------------------------
      INTEGER :: jl, jn, iyear_beg, iyear_tmp

      ! coefficient for CFC11 
      !----------------------

      ! Solubility
      soa(1,1) = -229.9261 
      soa(2,1) =  319.6552
      soa(3,1) =  119.4471
      soa(4,1) =   -1.39165

      sob(1,1) = -0.142382
      sob(2,1) =  0.091459
      sob(3,1) = -0.0157274

      ! Schmidt number          AXY (25/04/17)
      sca(1,1) = 3579.2       ! = 3501.8
      sca(2,1) = -222.63      ! = -210.31
      sca(3,1) =    7.5749    ! =    6.1851
      sca(4,1) =   -0.14595   ! =   -0.07513
      sca(5,1) =    0.0011874 ! = absent

      ! coefficient for CFC12 
      !----------------------

      ! Solubility
      soa(1,2) = -218.0971
      soa(2,2) =  298.9702
      soa(3,2) =  113.8049
      soa(4,2) =   -1.39165

      sob(1,2) = -0.143566
      sob(2,2) =  0.091015
      sob(3,2) = -0.0153924

      ! schmidt number         AXY (25/04/17)
      sca(1,2) = 3828.1      ! = 3845.4 
      sca(2,2) = -249.86     ! = -228.95
      sca(3,2) =    8.7603   ! =    6.1908 
      sca(4,2) =   -0.1716   ! =   -0.067430
      sca(5,2) =    0.001408 ! = absent

      ! coefficients for SF6   AXY (25/04/17)
      !---------------------
      
      ! Solubility
      soa(1,3) =  -80.0343
      soa(2,3) =  117.232
      soa(3,3) =   29.5817
      soa(4,3) =    0.0

      sob(1,3) =  0.0335183
      sob(2,3) = -0.0373942
      sob(3,3) =  0.00774862

      ! Schmidt number
      sca(1,3) = 3177.5
      sca(2,3) = -200.57
      sca(3,3) =    6.8865
      sca(4,3) =   -0.13335
      sca(5,3) =    0.0010877

      !!---------------------------------------------
      !! JPALM -- re-initialize CFC fields and diags if restart a CFC cycle,
      !!       Or if out of P_cfc range
      IF (simu_type==1) THEN
         iyear_tmp = nyear - nyear_res  !! JPALM -- in our spin-up, nyear_res is 1000
         iyear_beg = MOD( iyear_tmp , 90 )
         !!---------------------------------------
         IF ((iyear_beg > 77) .OR. (iyear_beg==0)) THEN
            qtr_cfc(:,:,:) = 0._wp
            IF(lwp) THEN
               WRITE(numout,*) 
               WRITE(numout,*) 'restart a CFC cycle or out of P_cfc year bounds zero --'
               WRITE(numout,*) '                          --    set qtr_CFC = 0.00   --'
               WRITE(numout,*) '                          --   set qint_CFC = 0.00   --'
               WRITE(numout,*) '                          --   set trn(CFC) = 0.00   --'
            ENDIF
            qtr_cfc(:,:,:) = 0._wp
            qint_cfc(:,:,:) = 0._wp
            trn(:,:,:,jp_cfc0:jp_cfc1) = 0._wp
            trb(:,:,:,jp_cfc0:jp_cfc1) = 0._wp
         ENDIF
      !!
      !! 2 -- Hindcast/proj
      ELSEIF (simu_type==2) THEN
         iyear_beg = MOD(nyear, 100)
         IF (iyear_beg < 20)  iyear_beg = iyear_beg + 100
         IF ((iyear_beg < 30) .OR. (iyear_beg > 115)) THEN
            qtr_cfc(:,:,:) = 0._wp
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'restart a CFC cycle or out of P_cfc year bounds zero --'
               WRITE(numout,*) '                          --    set qtr_CFC = 0.00   --'
               WRITE(numout,*) '                          --   set qint_CFC = 0.00   --'
               WRITE(numout,*) '                          --   set trn(CFC) = 0.00   --'
            ENDIF
            qtr_cfc(:,:,:) = 0._wp
            qint_cfc(:,:,:) = 0._wp
            trn(:,:,:,jp_cfc0:jp_cfc1) = 0._wp
            trb(:,:,:,jp_cfc0:jp_cfc1) = 0._wp
         ENDIF
      ENDIF

      IF(lwp) WRITE(numout,*)
      !
   END SUBROUTINE cfc_init


   INTEGER FUNCTION trc_sms_cfc_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_cfc_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( xphem   (jpi,jpj)        ,     &
         &      qtr_cfc (jpi,jpj,jp_cfc) ,     &
         &      qint_cfc(jpi,jpj,jp_cfc) , STAT=trc_sms_cfc_alloc )
         !
      IF( trc_sms_cfc_alloc /= 0 ) CALL ctl_warn('trc_sms_cfc_alloc : failed to allocate arrays.')
      !
   END FUNCTION trc_sms_cfc_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                         No CFC tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_cfc( kt )       ! Empty routine
      WRITE(*,*) 'trc_sms_cfc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_cfc
#endif

   !!======================================================================
END MODULE trcsms_cfc
