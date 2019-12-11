MODULE tideini
   !!======================================================================
   !!                       ***  MODULE  tideini  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  1.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE phycst
   USE daymod
   USE dynspg_oce
   USE tide_mod
   !
   USE iom
   USE in_out_manager  ! I/O units
   USE ioipsl          ! NetCDF IPSL library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PUBLIC

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   omega_tide   !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   v0tide       !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   utide        !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ftide        !:

   LOGICAL , PUBLIC ::   ln_tide_pot     !:
   LOGICAL , PUBLIC ::   ln_tide_ramp    !:
   INTEGER , PUBLIC ::   nb_harmo                 !:
   INTEGER , PUBLIC ::   kt_tide                  !:
   REAL(wp), PUBLIC ::   rdttideramp              !:
   
   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ntide   !:

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.5 , NEMO Consortium (2013)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
  SUBROUTINE tide_init ( kt )
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE tide_init  ***
    !!----------------------------------------------------------------------      
    !! * Local declarations
    INTEGER  :: ji, jk
    INTEGER, INTENT( in ) ::   kt     ! ocean time-step
    CHARACTER(LEN=4), DIMENSION(jpmax_harmo) :: clname
    INTEGER  ::   ios                 ! Local integer output status for namelist read
    !
    NAMELIST/nam_tide/ln_tide_pot, ln_tide_ramp, rdttideramp, clname
    !!----------------------------------------------------------------------

    IF ( kt == nit000 ) THEN
       !
       IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'tide_init : Initialization of the tidal components'
          WRITE(numout,*) '~~~~~~~~~ '
       ENDIF
       !
       CALL tide_init_Wave
       !
       ! Read Namelist nam_tide
       REWIND( numnam_ref )              ! Namelist nam_tide in reference namelist : Tides
       READ  ( numnam_ref, nam_tide, IOSTAT = ios, ERR = 901)
901    IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_tide in reference namelist', lwp )

       REWIND( numnam_cfg )              ! Namelist nam_tide in configuration namelist : Tides
       READ  ( numnam_cfg, nam_tide, IOSTAT = ios, ERR = 902 )
902    IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_tide in configuration namelist', lwp )
       IF(lwm) WRITE ( numond, nam_tide )
       !
       nb_harmo=0
       DO jk = 1, jpmax_harmo
          DO ji = 1,jpmax_harmo
             IF( TRIM(clname(jk)) == Wave(ji)%cname_tide )   nb_harmo = nb_harmo + 1
          END DO
       END DO
       !       
       ! Ensure that tidal components have been set in namelist_cfg
       IF( nb_harmo .EQ. 0 ) CALL ctl_stop( 'tide_init : No tidal components set in nam_tide' )
       !
       IF(lwp) THEN
          WRITE(numout,*) '   Namelist nam_tide'
          WRITE(numout,*) '      Apply astronomical potential : ln_tide_pot  =', ln_tide_pot
          WRITE(numout,*) '                                     nb_harmo     = ', nb_harmo
          WRITE(numout,*) '                                     ln_tide_ramp = ', ln_tide_ramp 
          WRITE(numout,*) '                                     rdttideramp  = ', rdttideramp
       ENDIF
       IF( ln_tide_ramp.AND.((nitend-nit000+1)*rdt/rday < rdttideramp) )   &
          &   CALL ctl_stop('rdttideramp must be lower than run duration')
       IF( ln_tide_ramp.AND.(rdttideramp<0.) ) &
          &   CALL ctl_stop('rdttideramp must be positive')
       !
       IF( .NOT. lk_dynspg_ts )   CALL ctl_warn( 'sbc_tide : use of time splitting is recommended' )
       !
       ALLOCATE( ntide(nb_harmo) )
       DO jk = 1, nb_harmo
          DO ji = 1, jpmax_harmo
             IF( TRIM(clname(jk)) .eq. Wave(ji)%cname_tide ) THEN
                ntide(jk) = ji
                EXIT
             END IF
          END DO
       END DO
       !
       ALLOCATE( omega_tide(nb_harmo), v0tide    (nb_harmo),   &
          &      utide     (nb_harmo), ftide     (nb_harmo)  )
       kt_tide = kt
       !
      ENDIF
      !
   END SUBROUTINE tide_init
     
   !!======================================================================
END MODULE tideini
