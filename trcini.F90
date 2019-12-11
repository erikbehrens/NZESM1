MODULE trcini
   !!======================================================================
   !!                         ***  MODULE trcini  ***
   !! TOP :   Manage the passive tracer initialization
   !!======================================================================
   !! History :   -   ! 1991-03 (O. Marti)  original code
   !!            1.0  ! 2005-03 (O. Aumont, A. El Moussaoui) F90
   !!            2.0  ! 2005-10 (C. Ethe, G. Madec) revised architecture
   !!            4.0  ! 2011-01 (A. R. Porter, STFC Daresbury) dynamical allocation
   !!             -   ! 2014-06 (A. Yool, J. Palmieri) adding MEDUSA-2
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_init  :   Initialization for passive tracer
   !!   top_alloc :   allocate the TOP arrays
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables
   USE trcrst          ! passive tracers restart
   USE trcnam          ! Namelist read
   USE trcini_cfc      ! CFC      initialisation
   USE trcini_pisces   ! PISCES   initialisation
   USE trcini_c14b     ! C14 bomb initialisation
   USE trcini_age      ! AGE      initialisation
   USE trcini_my_trc   ! MY_TRC   initialisation
   USE trcini_idtra    ! idealize tracer initialisation
   USE trcini_medusa   ! MEDUSA   initialisation
   USE par_medusa      ! MEDUSA   parameters (needed for elemental cycles)
   USE trcdta          ! initialisation from files
   USE daymod          ! calendar manager
   USE prtctl_trc      ! Print control passive tracers (prt_ctl_trc_init routine)
   USE trcsub          ! variables to substep passive tracers
   USE lib_mpp         ! distribued memory computing library
   USE sbc_oce
   USE trcice          ! tracers in sea ice
   USE sms_medusa      ! MEDUSA   initialisation
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_init   ! called by opa

    !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE trc_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_init  ***
      !!
      !! ** Purpose :   Initialization of the passive tracer fields 
      !!
      !! ** Method  : - read namelist
      !!              - control the consistancy 
      !!              - compute specific initialisations
      !!              - set initial tracer fields (either read restart 
      !!                or read data or analytical formulation
      !!---------------------------------------------------------------------
      INTEGER ::   ji, jj, jk, jn, jl    ! dummy loop indices
# if defined key_medusa && defined key_roam
      !! AXY (23/11/2017)
      REAL(wp)                         :: zsum3d, zsum2d
      REAL(wp)                         :: zq1, zq2, loc_vol, loc_area
      REAL(wp), DIMENSION(6)           :: loc_cycletot3, loc_cycletot2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztot3d
      REAL(wp), DIMENSION(jpi,jpj)     :: ztot2d, carea
# endif
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_init')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_init : initial set up of the passive tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      CALL top_alloc()              ! allocate TOP arrays

      l_trcdm2dc = ln_dm2dc .OR. ( ln_cpl .AND. ncpl_qsr_freq /= 1 )
      l_trcdm2dc = l_trcdm2dc  .AND. .NOT. lk_offline
      IF( l_trcdm2dc .AND. lwp ) &
         &   CALL ctl_warn(' Coupling with passive tracers and used of diurnal cycle. &
         & Computation of a daily mean shortwave for some biogeochemical models) ')

      IF( nn_cla == 1 )   &
         &  CALL ctl_stop( ' Cross Land Advection not yet implemented with passive tracer ; nn_cla must be 0' )

      CALL trc_nam      ! read passive tracers namelists
      !
      IF(lwp) WRITE(numout,*)
      !
      IF( ln_rsttr .AND. .NOT. lk_offline ) CALL trc_rst_cal( nit000, 'READ' )   ! calendar
      !
      IF(lwp) WRITE(numout,*)
                                                              ! masked grid volume
      !                                                              ! masked grid volume
      DO jk = 1, jpk
         cvol(:,:,jk) = e1e2t(:,:) * fse3t(:,:,jk) * tmask(:,:,jk)
      END DO
      IF( lk_degrad ) cvol(:,:,:) = cvol(:,:,:) * facvol(:,:,:)      ! degrad option: reduction by facvol
      !                                                              ! total volume of the ocean 
      areatot = glob_sum( cvol(:,:,:) )
      carea(:,:) = e1e2t(:,:) * tmask(:,:,1) 

      IF( lk_pisces  )       CALL trc_ini_pisces       ! PISCES  bio-model
      IF( lk_cfc     )       CALL trc_ini_cfc          ! CFC     tracers
      IF( lk_c14b    )       CALL trc_ini_c14b         ! C14 bomb  tracer
      IF( lk_age     )       CALL trc_ini_age          ! AGE       tracer
      IF( lk_my_trc  )       CALL trc_ini_my_trc       ! MY_TRC  tracers
      IF( lk_idtra   )       CALL trc_ini_idtra        ! Idealize tracers
      IF( lk_medusa  )       CALL trc_ini_medusa       ! MEDUSA  tracers

      CALL trc_ice_ini                                 ! Tracers in sea ice

      IF( ln_ctl ) THEN
         !
         IF (narea == 1) THEN  
            ! The tracer.stat file only contains global tracer sum values, if 
            ! it contains anything at all. Hence it only needs to be opened 
            ! and written to on the master PE, not on all PEs.  
            CALL ctl_opn( numstr, 'tracer.stat', 'REPLACE','FORMATTED',  & 
                          'SEQUENTIAL', -1, numout, lwp , narea ) 
         ENDIF  
         !
      ENDIF
      IF(lwp) WRITE(numout,*) 'before trc_dta_init'
      IF( ln_trcdta ) THEN
         CALL trc_dta_init(jptra)
      ENDIF
      IF(lwp) WRITE(numout,*) 'after trc_dta_init'
      IF( ln_rsttr ) THEN
        !
        CALL trc_rst_read              ! restart from a file
        !
      ELSE
        !
        IF(lwp) WRITE(numout,*) 'after restart'
        IF( ln_trcdta .AND. nb_trcdta > 0 ) THEN  ! Initialisation of tracer from a file that may also be used for damping
            !
            DO jn = 1, jptra
               IF( ln_trc_ini(jn) ) THEN      ! update passive tracers arrays with input data read from file
                  jl = n_trc_index(jn) 
                  CALL trc_dta( nit000, sf_trcdta(jl), rf_trfac(jl) )   ! read tracer data at nit000
                  trn(:,:,:,jn) = sf_trcdta(jl)%fnow(:,:,:) 
                  IF( .NOT.ln_trcdmp .AND. .NOT.ln_trcdmp_clo ) THEN      !== deallocate data structure   ==!
                     !                                                    (data used only for initialisation)
                     IF(lwp) WRITE(numout,*) 'trc_dta: deallocate data arrays as they are only used to initialize the run'
                                                  DEALLOCATE( sf_trcdta(jl)%fnow )     !  arrays in the structure
                     IF( sf_trcdta(jl)%ln_tint )  DEALLOCATE( sf_trcdta(jl)%fdta )
                     !
                  ENDIF
               ENDIF
            ENDDO
            !
        ENDIF
        !
        trb(:,:,:,:) = trn(:,:,:,:)
        ! 
        IF(lwp) WRITE(numout,*) 'after trc=trn'
      ENDIF
 
      tra(:,:,:,:) = 0._wp
      !
      IF( nn_dttrc /= 1 )        CALL trc_sub_ini      ! Initialize variables for substepping passive tracers
      !
      IF(lwp) WRITE(numout,*) 'after trc_sub_ini'
      trai(:) = 0._wp                                                   ! initial content of all tracers
      !CEB
      IF( Agrif_Root() ) THEN
      DO jn = 1, jptra
         trai(jn) = trai(jn) + glob_sum( trn(:,:,:,jn) * cvol(:,:,:)   )
         IF(lwp) WRITE(numout,*) 'global numbers', jn
      END DO
      ENDIF

      IF(lwp) THEN               ! control print
         WRITE(numout,*)
         WRITE(numout,*)
         WRITE(numout,*) '          *** Total number of passive tracer jptra = ', jptra
         WRITE(numout,*) '          *** Total volume of ocean                = ', areatot
         WRITE(numout,*) '          *** Total inital content of all tracers '
         WRITE(numout,*)
# if defined key_debug_medusa
         CALL flush(numout)
# endif
         !
# if defined key_debug_medusa
         WRITE(numout,*) ' litle check :  ', ctrcnm(1)
         CALL flush(numout)
# endif
         DO jn = 1, jptra
            WRITE(numout,9000) jn, TRIM( ctrcnm(jn) ), trai(jn)
         ENDDO
         WRITE(numout,*)
      ENDIF
      IF(lwp) WRITE(numout,*)
      IF(ln_ctl) THEN            ! print mean trends (used for debugging)
         CALL prt_ctl_trc_init
         WRITE(charout, FMT="('ini ')")
         CALL prt_ctl_trc_info( charout )
         CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF

      !CEB
      IF( Agrif_Root() ) THEN
# if defined key_medusa && defined key_roam
      ! AXY (17/11/2017): calculate initial totals of elemental cycles
      !
      ! This is done in a very hard-wired way here; in future, this could be
      ! replaced with loops and using a 2D array; one dimension would cover
      ! the tracers, the other would be for the elements; each tracer would
      ! have a factor for each element to say how much of that element was
      ! in that tracer; for example, PHN would be 1.0 for N, xrfn for Fe and
      ! xthetapn for C, with the other elements 0.0; the array entry for PHN
      ! would then be (1. 0. xrfn xthetapn 0. 0.) for (N, Si, Fe, C, A, O2);
      ! saving this for the next iteration
      !
      cycletot(:) = 0._wp
      ! report elemental totals at initialisation as we go along
      IF ( lwp ) WRITE(numout,*)
      IF ( lwp ) WRITE(numout,*)    ' Elemental cycle totals: '
      ! nitrogen
      ztot3d(:,:,:) = trn(:,:,:,jpphn) + trn(:,:,:,jpphd) + trn(:,:,:,jpzmi) + &
                      trn(:,:,:,jpzme) + trn(:,:,:,jpdet) + trn(:,:,:,jpdin)
      ztot2d(:,:)   = zn_sed_n(:,:)
      zsum3d        = glob_sum( ztot3d(:,:,:) * cvol(:,:,:) )
      zsum2d        = glob_sum( ztot2d(:,:) * carea(:,:) )
      cycletot(1)   = zsum3d + zsum2d
      IF ( lwp ) WRITE(numout,9010) 'nitrogen', zsum3d, zsum2d, cycletot(1)
      ! silicon
      ztot3d(:,:,:) = trn(:,:,:,jppds) + trn(:,:,:,jpsil)
      ztot2d(:,:)   = zn_sed_si(:,:)
      zsum3d        = glob_sum( ztot3d(:,:,:) * cvol(:,:,:) )
      zsum2d        = glob_sum( ztot2d(:,:) * carea(:,:) )
      cycletot(2)   = zsum3d + zsum2d
      IF ( lwp ) WRITE(numout,9010) 'silicon', zsum3d, zsum2d, cycletot(2)
      ! iron
      ztot3d(:,:,:) = ((trn(:,:,:,jpphn) + trn(:,:,:,jpphd) + trn(:,:,:,jpzmi) + &
                      trn(:,:,:,jpzme) + trn(:,:,:,jpdet)) * xrfn) + trn(:,:,:,jpfer)
      ztot2d(:,:)   = zn_sed_fe(:,:)
      zsum3d        = glob_sum( ztot3d(:,:,:) * cvol(:,:,:) )
      zsum2d        = glob_sum( ztot2d(:,:) * carea(:,:) )
      cycletot(3)   = zsum3d + zsum2d
      IF ( lwp ) WRITE(numout,9010) 'iron', zsum3d, zsum2d, cycletot(3)
      ! carbon (uses fixed C:N ratios on plankton tracers)
      ztot3d(:,:,:) = (trn(:,:,:,jpphn) * xthetapn)  + (trn(:,:,:,jpphd) * xthetapd)  +  &
                      (trn(:,:,:,jpzmi) * xthetazmi) + (trn(:,:,:,jpzme) * xthetazme) +  &
                      trn(:,:,:,jpdtc) + trn(:,:,:,jpdic)
      ztot2d(:,:)   = zn_sed_c(:,:) + zn_sed_ca(:,:)
      zsum3d        = glob_sum( ztot3d(:,:,:) * cvol(:,:,:) )
      zsum2d        = glob_sum( ztot2d(:,:) * carea(:,:) )
      cycletot(4)   = zsum3d + zsum2d
      IF ( lwp ) WRITE(numout,9010) 'carbon', zsum3d, zsum2d, cycletot(4)
      ! alkalinity (note benthic correction)
      ztot3d(:,:,:) = trn(:,:,:,jpalk)
      ztot2d(:,:)   = zn_sed_ca(:,:) * 2._wp
      zsum3d        = glob_sum( ztot3d(:,:,:) * cvol(:,:,:) )
      zsum2d        = glob_sum( ztot2d(:,:) * carea(:,:) )
      cycletot(5)   = zsum3d + zsum2d
      IF ( lwp ) WRITE(numout,9010) 'alkalinity', zsum3d, zsum2d, cycletot(5)
      ! oxygen (note no benthic)
      ztot3d(:,:,:) = trn(:,:,:,jpoxy)
      ztot2d(:,:)   = 0._wp
      zsum3d        = glob_sum( ztot3d(:,:,:) * cvol(:,:,:) )
      zsum2d        = glob_sum( ztot2d(:,:) * carea(:,:) )
      cycletot(6)   = zsum3d + zsum2d
      IF ( lwp ) WRITE(numout,9010) 'oxygen', zsum3d, zsum2d, cycletot(6)
      ! Check
      zsum3d        = glob_sum( cvol(:,:,:) )
      zsum2d        = glob_sum( carea(:,:) )      
      IF ( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' check : cvol    : ', zsum3d
         WRITE(numout,*) ' check : carea   : ', zsum2d
         WRITE(numout,*)
      ENDIF
      !
      ENDIF
      !/CEB
# endif

      IF(lwp) THEN 
          WRITE(numout,*)
          WRITE(numout,*) 'trc_init : passive tracer set up completed'
          WRITE(numout,*) '~~~~~~~'
      ENDIF 
# if defined key_debug_medusa
         CALL trc_rst_stat
         CALL flush(numout)
# endif

9000  FORMAT(' tracer nb : ',i2,'      name :',a10,'      initial content :',e18.10)
9010  FORMAT(' element:',a10,                     &
             ' 3d sum:',e18.10,' 2d sum:',e18.10, &
             ' total:',e18.10)
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_init')
      !
   END SUBROUTINE trc_init


   SUBROUTINE top_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE top_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!----------------------------------------------------------------------
      USE trcadv        , ONLY:   trc_adv_alloc          ! TOP-related alloc routines...
      USE trc           , ONLY:   trc_alloc
      USE trcnxt        , ONLY:   trc_nxt_alloc
      USE trczdf        , ONLY:   trc_zdf_alloc
      USE trdtrc_oce    , ONLY:   trd_trc_oce_alloc
#if defined key_trdmxl_trc 
      USE trdmxl_trc    , ONLY:   trd_mxl_trc_alloc
#endif
# if defined key_medusa
      USE bio_medusa_mod, ONLY:   bio_medusa_alloc
# endif

      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        trc_adv_alloc()          ! Start of TOP-related alloc routines...
      ierr = ierr + trc_alloc    ()
      ierr = ierr + trc_nxt_alloc()
      ierr = ierr + trc_zdf_alloc()
      ierr = ierr + trd_trc_oce_alloc()
#if defined key_trdmxl_trc 
      ierr = ierr + trd_mxl_trc_alloc()
#endif
#if defined key_medusa
      ierr = ierr + bio_medusa_alloc()
#endif
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'top_alloc : unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE top_alloc

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_init                      ! Dummy routine   
   END SUBROUTINE trc_init
#endif

   !!======================================================================
END MODULE trcini
