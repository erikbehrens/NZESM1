MODULE trcice
   !!======================================================================
   !!                         ***  MODULE trcice  ***
   !! TOP :   Manage the communication between TOP and sea ice
   !!======================================================================
   !! History :  3.5  ! 2013    (M. Vancoppenolle, O. Aumont, G. Madec), original code
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_ice   :  Call the appropriate sea ice tracer subroutine
   !!----------------------------------------------------------------------

   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables
   USE trcice_cfc      ! CFC      initialisation
   USE trcice_pisces   ! PISCES   initialisation
   USE trcice_c14b     ! C14 bomb initialisation
   USE trcice_my_trc   ! MY_TRC   initialisation
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_ice_ini ! called by trc_nam

CONTAINS
   
   SUBROUTINE trc_ice_ini
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ice_ini ***
      !!
      !! ** Purpose :   Initialization of the ice module for tracers
      !!
      !! ** Method  : - 
      !!            
      !!---------------------------------------------------------------------
      ! --- Variable declarations --- !

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_ice_ini : Initialize sea ice tracer boundary condition'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF

      IF( nn_timing == 1 )  CALL timing_start('trc_ice_ini')

      !
      trc_i(:,:,:) = 0.0d0 ! by default
      trc_o(:,:,:) = 0.0d0 ! by default

      IF ( nn_ice_tr == 1 ) THEN
         IF( lk_pisces  )    CALL trc_ice_ini_pisces       ! PISCES  bio-model
         IF( lk_cfc     )    CALL trc_ice_ini_cfc          ! CFC     tracers
         IF( lk_c14b    )    CALL trc_ice_ini_c14b         ! C14 bomb  tracer
         IF( lk_my_trc  )    CALL trc_ice_ini_my_trc       ! MY_TRC  tracers
      ENDIF
      !CEB
      WRITE(numout,*) 'after trc ice init' 
      IF( nn_timing == 1 )   CALL timing_stop('trc_ice_ini')
      !
   END SUBROUTINE trc_ice_ini

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ice_ini                   ! Dummy routine   
   END SUBROUTINE trc_ice_ini
#endif

   !!======================================================================
END MODULE trcice
