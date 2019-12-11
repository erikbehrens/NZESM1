MODULE trdtrc
   !!======================================================================
   !!                       ***  MODULE  trdtrc  ***
   !! Ocean diagnostics:  mixed layer passive tracer trends 
   !!======================================================================
   !! History :  3.0  !  2010-07  (C. Ethe)  Original code (from trdtrc.F90)
   !!----------------------------------------------------------------------
#if   defined key_top && ( defined key_trdmxl_trc   ||   defined key_trdtrc )
   !!----------------------------------------------------------------------
   !!   'key_trdmxl_trc'                  mixed layer trend diagnostics
   !!   'key_trdtrc'                      3D trend diagnostics
   !!----------------------------------------------------------------------
   !!   trdtrc      : passive tracer trends 
   !!----------------------------------------------------------------------
   USE trc               ! tracer definitions (trn, trb, tra, etc.)
   USE trcnam_trp
   USE trd_oce
   USE trdtrc_oce       ! definition of main arrays used for trends computations
   USE trdmxl_trc        ! Mixed layer trends diag.
   USE iom               ! I/O library
# if defined key_debug_medusa
   USE trcstat,          ONLY: trc_rst_dia_stat     
# endif

   IMPLICIT NONE
   PRIVATE

   INTERFACE trd_trc
      MODULE PROCEDURE trd_trc_trp, trd_trc_bio
   END INTERFACE

   PUBLIC trd_trc

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trd_trc_trp( ptrtrd, kjn, ktrd, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_trc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      INTEGER, INTENT( in )  ::   ktrd                                ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER (len=20) :: cltra
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
!         IF(lwp)WRITE(numout,*)
!         IF(lwp)WRITE(numout,*) 'trd_trc:'
!         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Mixed layer trends for passive tracers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if defined key_trdmxl_trc  
      IF( lk_trdmxl_trc .AND. ln_trdtrc( kjn ) ) THEN
         !
         SELECT CASE ( ktrd )
         CASE ( jptra_xad     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_xad, '3D', kjn )
         CASE ( jptra_yad     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_yad, '3D', kjn )
         CASE ( jptra_zad     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_zad, '3D', kjn )
         CASE ( jptra_ldf     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_ldf, '3D', kjn )
         CASE ( jptra_bbl     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_bbl, '3D', kjn )
         CASE ( jptra_zdf     )
            IF( ln_trcldf_iso ) THEN
               CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_ldf, '3D', kjn )
            ELSE
               CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_zdf, '3D', kjn )
            ENDIF
         CASE ( jptra_dmp     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_dmp , '3D', kjn )
         CASE ( jptra_nsr     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_sbc , '2D', kjn )
         CASE ( jptra_sms     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_sms , '3D', kjn )
         CASE ( jptra_radb    )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_radb, '3D', kjn )
         CASE ( jptra_radn    )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_radn, '3D', kjn )
         CASE ( jptra_atf     )   ;   CALL trd_mxl_trc_zint( ptrtrd, jpmxl_trc_atf , '3D', kjn )
         END SELECT
         !
      END IF
#endif

      IF( lk_trdtrc .AND. ln_trdtrc( kjn ) ) THEN
      !! JPALM -- 17-08-2017 -- modif following trd_tra_iom as suggested by Georges
      !!                     -- add jptra_tot; jptra_totad; jptra_zdfp
      !!                     -- shange to output trends every 2 time-step, except tot.
      !!                     -- move cltra and iomput inside the select case
      !!                     So if an non-wanted case arrives here it will not go
      !!                     through cltra (without value) and break iomput.
      !!                     -- Add iom_use in prevision of not using All trends
      !!                     for All passive tracers (will create a HUGE 3D file otherwise --
      !!                     might be interested in very few of them : SMS and TOT probably)
         !
         SELECT CASE( ktrd )
         !! tot - output every time-step:
         CASE( jptra_tot  )       ;    WRITE (cltra,'("TOT_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         END SELECT
         !
       IF( MOD( kt, 2 ) == 0 ) THEN
         SELECT CASE( ktrd )
         CASE( jptra_xad  )       ;    WRITE (cltra,'("XAD_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_yad  )       ;    WRITE (cltra,'("YAD_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_zad  )       ;    WRITE (cltra,'("ZAD_",4a)')      !! care vvl case
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_totad  )     ;    WRITE (cltra,'("TAD_",4a)')      !! total adv
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_ldf  )       ;    WRITE (cltra,'("LDF_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_bbl  )       ;    WRITE (cltra,'("BBL_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_nsr  )       ;    WRITE (cltra,'("FOR_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_zdf  )       ;    WRITE (cltra,'("ZDF_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_zdfp )       ;    WRITE (cltra,'("ZDP_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_dmp  )       ;    WRITE (cltra,'("DMP_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_sms  )       ;    WRITE (cltra,'("SMS_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_radb )       ;    WRITE (cltra,'("RDB_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         CASE( jptra_radn )       ;    WRITE (cltra,'("RDN_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         END SELECT
       ELSE IF( MOD( kt, 2 ) == 1 ) THEN
         SELECT CASE( ktrd )
         CASE( jptra_atf  )       ;    WRITE (cltra,'("ATF_",4a)')
                           cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
                           CALL trd_trc_iomput( cltra, ptrtrd, kjn, kt )
         END SELECT
       END IF
         !
      END IF

   END SUBROUTINE trd_trc_trp

   SUBROUTINE trd_trc_bio( ptrbio, ktrd, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_bio  ***
      !!----------------------------------------------------------------------

      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   ktrd                                ! bio trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout )  ::   ptrbio  ! Bio trend
      !!----------------------------------------------------------------------

#if defined key_trdmxl_trc  
      CALL trd_mxl_bio_zint( ptrbio, ktrd ) ! Verticaly integrated biological trends
#endif

   END SUBROUTINE trd_trc_bio

   SUBROUTINE trd_trc_iomput( cltra, ptrtrd, kjn, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_trc_iomput  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! timestep
      INTEGER, INTENT( in )  ::   kjn                                 ! biotrend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout )  ::   ptrtrd  ! var trend
      CHARACTER (len=*),INTENT( in ) :: cltra                         ! trend name
      !!----------------------------------------------------------------------


      IF  (iom_use(cltra)) THEN
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) ' TREND stats (min, max,sum) kt = ',kt ,' jn = ',kjn
         CALL trc_rst_dia_stat( ptrtrd(:,:,1), cltra)
# endif
         CALL iom_put( cltra,  ptrtrd(:,:,:) )
# if defined key_debug_medusa
      ELSE
         IF(lwp) WRITE(numout,*) &
                      ' TREND -- No output asked for ',cltra,' kt = ',kt,' jn = ',kjn
         CALL trc_rst_dia_stat( ptrtrd(:,:,1), cltra)
# endif
      ENDIF

   END SUBROUTINE trd_trc_iomput


#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------

   INTERFACE trd_trc
      MODULE PROCEDURE trd_trc_trp, trd_trc_bio
   END INTERFACE

CONTAINS

   SUBROUTINE trd_trc_trp( ptrtrd, kjn, ktrd, kt )
      INTEGER               , INTENT( in )     ::   kt      ! time step
      INTEGER               , INTENT( in )     ::   kjn     ! tracer index
      INTEGER               , INTENT( in )     ::   ktrd    ! tracer trend index
      REAL, DIMENSION(:,:,:), INTENT( inout )  ::   ptrtrd  ! Temperature or U trend
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_trc_trp

   SUBROUTINE trd_trc_bio( ptrbio, ktrd, kt )
      INTEGER               , INTENT( in )     ::   kt      ! time step
      INTEGER               , INTENT( in )     ::   ktrd    ! tracer trend index
      REAL, DIMENSION(:,:,:), INTENT( inout )  ::   ptrbio  ! Temperature or U trend
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrbio(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_trc_bio

#endif
   !!======================================================================
END MODULE trdtrc
