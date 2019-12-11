MODULE trcini_idtra
   !!======================================================================
   !!                         ***  MODULE trcini_idtra  ***
   !! TOP :   initialisation of the IDEAL-TRACER tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.idtra.h90
   !!----------------------------------------------------------------------
#if defined key_idtra
   !!----------------------------------------------------------------------
   !!   'key_idtra'                                               IDEAL-TRACER tracers
   !!----------------------------------------------------------------------
   !! trc_ini_idtra      : IDEAL-TRACER model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_idtra    ! IDEAL-TRACER sms trends
   ! USE par_idtra       ! IDEAL-TRACER parameters
   ! USE in_out_manager  ! I/O manager
   ! USE lib_mpp
   ! USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_idtra   ! called by trcini.F90 module

   INTEGER  ::   inum                   ! unit number

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_idtra
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_idtra  ***
      !!
      !! ** Purpose :   initialization for idtra model
      !!
      !! ** Method  : - Read the namidtra namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER  ::    jn, jl 
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_idtra: initialisation of Ideal Tracers model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~'

      IF( trc_sms_idtra_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_idtra:unable to allocate CFC arrays' )


      ! Initialization of trn in case of  no restart
      !----------------------------------------------
      qtr_idtra(:,:,:) = 0._wp
      inv_idtra(:,:,:) = 0._wp
      IF( .NOT. ln_rsttr ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'Initialization of id-tracers ; No restart : '
            WRITE(numout,*) '                             ; Init field equal 1 at surface - zero elsewhere'
            WRITE(numout,*) '                             ; qint idtra equal 0 '
         ENDIF
         qint_idtra(:,:,:) = 0._wp
         DO jn = jp_idtra0, jp_idtra1
             trn(:,:,:,jn) = 0.e0
             trn(:,:,1,jn) = 1.0
           IF(lwp) WRITE(numout,*) 'Idealise Tracer initialisation -- jn = ',jn
         END DO
      ENDIF


      !   Ideal traceur do not need any atmospheric concentration.
      ! We consider that sucface concentration is equal to 1,
      ! that it is advectied within the water circulation,
      ! and that it is regularly degraded as if it was a radiactive tracer (tricium for example)
      ! But we can play with tha caracteristic time of
      !--------------------------------------------------------------------



      IF(lwp) WRITE(numout,*) 'Initialization of IDEAL-TRACER tracers done'
      IF(lwp) WRITE(numout,*) ' '

   END SUBROUTINE trc_ini_idtra

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                         No IDEAL-TRACER tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_idtra             ! Empty routine


   END SUBROUTINE trc_ini_idtra
#endif

   !!======================================================================
END MODULE trcini_idtra
                                                



