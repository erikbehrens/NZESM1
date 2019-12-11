MODULE trcsms_age
   !!======================================================================
   !!                         ***  MODULE trcsms_age  ***
   !! TOP :   Main module of the AGE tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_age
   !!----------------------------------------------------------------------
   !!   'key_age'                                               AGE tracer
   !!----------------------------------------------------------------------
   !! trc_sms_age       : AGE model main routine
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_age       ! called by trcsms.F90 module

   INTEGER , PUBLIC :: nl_age             ! T level surrounding age_depth
   INTEGER , PUBLIC :: nla_age            ! T level wholly above age_depth
   INTEGER , PUBLIC :: nlb_age            ! T level wholly below age_depth

   REAL(wp), PUBLIC :: rn_age_depth       ! = 10       depth over which age tracer reset to zero
   REAL(wp), PUBLIC :: rn_age_kill_rate   ! = -1./7200  recip of relaxation timescale (s) for  age tracer shallower than age_depth
   
   REAL(wp), PUBLIC :: rryear          !: recip number of seconds in one year
   REAL(wp), PUBLIC :: frac_kill_age   !: fraction of level nl_age above age_depth where it is relaxed towards zero
   REAL(wp), PUBLIC :: frac_add_age    !: fraction of level nl_age below age_depth where it is incremented


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_age( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_age  ***
      !!
      !! ** Purpose :   main routine of AGE model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::   jn, jk   ! dummy loop index
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrage
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_age')
      !
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_sms_age:  AGE model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      ENDIF

      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrage )

      DO jk = 1, nla_age
         tra(:,:,jk,jpage1) = rn_age_kill_rate * trb(:,:,jk,jpage1)
      ENDDO
      !
      tra(:,:,nl_age,jpage1) = frac_kill_age * rn_age_kill_rate * trb(:,:,nl_age,jpage1)  &
          &                  + frac_add_age  * rryear * tmask(:,:,nl_age)
      !
      DO jk = nlb_age, jpk
         tra(:,:,jk,jpage1) = tmask(:,:,jk) * rryear
      ENDDO
      !
      IF( l_trdtrc ) THEN      ! Save the trends in the mixed layer
          DO jn = jp_age0, jp_age1
            ztrage(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrage, jn, jptra_sms, kt )   ! save trends
          END DO
          CALL wrk_dealloc( jpi, jpj, jpk, ztrage )
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_age')
      !
   END SUBROUTINE trc_sms_age

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No AGE model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_age( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_age: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_age
#endif

   !!======================================================================
END MODULE trcsms_age
