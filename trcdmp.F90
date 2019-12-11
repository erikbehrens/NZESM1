MODULE trcdmp
   !!======================================================================
   !!                       ***  MODULE  trcdmp  ***
   !! Ocean physics: internal restoring trend on passive tracers
   !!======================================================================
   !! History :  OPA  !  1991-03  (O. Marti, G. Madec)  Original code
   !!                 !  1996-01  (G. Madec) statement function for e3
   !!                 !  1997-05  (H. Loukos)  adapted for passive tracers
   !!    NEMO    9.0  !  2004-03  (C. Ethe)    free form + modules
   !!            3.2  !  2007-02  (C. Deltel)  Diagnose ML trends for passive tracers
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) merge TRA-TRC 
   !!----------------------------------------------------------------------
#if  defined key_top 
   !!----------------------------------------------------------------------
   !!   trc_dmp      : update the tracer trend with the internal damping
   !!   trc_dmp_init : initialization, namlist read, parameters control
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE trcdta
   USE tradmp
   USE prtctl_trc      ! Print control for debbuging
   USE trdtra
   USE trd_oce
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_dmp            ! routine called by step.F90
   PUBLIC trc_dmp_clo        ! routine called by step.F90
   PUBLIC trc_dmp_alloc      ! routine called by nemogcm.F90

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   restotr   ! restoring coeff. on tracers (s-1)

   INTEGER, PARAMETER           ::   npncts   = 8        ! number of closed sea
   INTEGER, DIMENSION(npncts)   ::   nctsi1, nctsj1      ! south-west closed sea limits (i,j)
   INTEGER, DIMENSION(npncts)   ::   nctsi2, nctsj2      ! north-east closed sea limits (i,j)

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_dmp_alloc()
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( restotr(jpi,jpj,jpk) , STAT=trc_dmp_alloc )
      !
      IF( trc_dmp_alloc /= 0 )   CALL ctl_warn('trc_dmp_alloc: failed to allocate array')
      !
   END FUNCTION trc_dmp_alloc


   SUBROUTINE trc_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dmp  ***
      !!                  
      !! ** Purpose :   Compute the passive tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards trdta computed 
      !!      and add to the general tracer trends:
      !!                     trn = tra + restotr * (trdta - trb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmptr=0) or in area of weak vertical mixing (nlmdmptr=1) or
      !!      below the well mixed layer (nlmdmptr=2)
      !!
      !! ** Action  : - update the tracer trends tra with the newtonian 
      !!                damping trends.
      !!              - save the trends ('key_trdmxl_trc')
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jn, jl       ! dummy loop indices
      REAL(wp) ::   ztra                 ! temporary scalars
      CHARACTER (len=22) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrtrd
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrcdta   ! 3D  workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_dmp')
      !
      ! 0. Initialization (first time-step only)
      !    --------------
      IF( kt == nittrc000 ) CALL trc_dmp_init

      IF( l_trdtrc )   CALL wrk_alloc( jpi, jpj, jpk, ztrtrd )   ! temporary save of trends
      !
      IF( nb_trcdta > 0 ) THEN  ! Initialisation of tracer from a file that may also be used for damping
         !
         CALL wrk_alloc( jpi, jpj, jpk, ztrcdta )    ! Memory allocation
         !                                                          ! ===========
         DO jn = 1, jptra                                           ! tracer loop
            !                                                       ! ===========
            IF( l_trdtrc ) ztrtrd(:,:,:) = tra(:,:,:,jn)    ! save trends 
            !
            IF( ln_trc_ini(jn) ) THEN      ! update passive tracers arrays with input data read from file
               
               jl = n_trc_index(jn) 
               CALL trc_dta( kt, sf_trcdta(jl), rf_trfac(jl), ztrcdta )   ! read tracer data at nit000

               SELECT CASE ( nn_zdmp_tr )
               !
               CASE( 0 )                !==  newtonian damping throughout the water column  ==!
                  DO jk = 1, jpkm1
                     DO jj = 2, jpjm1
                        DO ji = fs_2, fs_jpim1   ! vector opt.
                           ztra = restotr(ji,jj,jk) * ( ztrcdta(ji,jj,jk) - trb(ji,jj,jk,jn) )
                           tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
                        END DO
                     END DO
                  END DO
               !
               CASE ( 1 )                !==  no damping in the turbocline (avt > 5 cm2/s)  ==!
                  DO jk = 1, jpkm1
                     DO jj = 2, jpjm1
                        DO ji = fs_2, fs_jpim1   ! vector opt.
                           IF( avt(ji,jj,jk) <= 5.e-4_wp )  THEN 
                              ztra = restotr(ji,jj,jk) * ( ztrcdta(ji,jj,jk) - trb(ji,jj,jk,jn) )
                              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
                           ENDIF
                        END DO
                     END DO
                  END DO
               !
               CASE ( 2 )               !==  no damping in the mixed layer   ==! 
                  DO jk = 1, jpkm1
                     DO jj = 2, jpjm1
                        DO ji = fs_2, fs_jpim1   ! vector opt.
                           IF( fsdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                              ztra = restotr(ji,jj,jk) * ( ztrcdta(ji,jj,jk) - trb(ji,jj,jk,jn) )
                              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
                           END IF
                        END DO
                     END DO
                  END DO
               !  
               END SELECT
               ! 
            ENDIF
            !
            IF( l_trdtrc ) THEN
               ztrtrd(:,:,:) = tra(:,:,:,jn) -  ztrtrd(:,:,:)
               CALL trd_tra( kt, 'TRC', jn, jptra_dmp, ztrtrd )
            END IF
            !                                                       ! ===========
         END DO                                                     ! tracer loop
         !                                                          ! ===========
         CALL wrk_dealloc( jpi, jpj, jpk, ztrcdta )
      ENDIF
      !
      IF( l_trdtrc )  CALL wrk_dealloc( jpi, jpj, jpk, ztrtrd )
      !                                          ! print mean trends (used for debugging)
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('dmp ')") ;  CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_dmp')
      !
   END SUBROUTINE trc_dmp

   SUBROUTINE trc_dmp_clo( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc_dmp_clo  ***
      !!
      !! ** Purpose :   Closed sea domain initialization
      !!
      !! ** Method  :   if a closed sea is located only in a model grid point
      !!                we restore to initial data
      !!
      !! ** Action  :   nctsi1(), nctsj1() : south-west closed sea limits (i,j)
      !!                nctsi2(), nctsj2() : north-east Closed sea limits (i,j)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER :: ji , jj, jk, jn, jl, jc                     ! dummy loop indicesa
      INTEGER :: isrow                                      ! local index
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrcdta       ! 3D  workspace

      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         ! initial values
         nctsi1(:) = 1  ;  nctsi2(:) = 1
         nctsj1(:) = 1  ;  nctsj2(:) = 1

         ! set the closed seas (in data domain indices)
         ! -------------------

         IF( cp_cfg == "orca" ) THEN
            !
            SELECT CASE ( jp_cfg )
            !                                           ! =======================
            CASE ( 1 )                                  ! eORCA_R1 configuration
            !                                           ! =======================
            isrow = 332 - jpjglo
            !
                                                        ! Caspian Sea
            nctsi1(1)   = 333  ; nctsj1(1)   = 243 - isrow
            nctsi2(1)   = 342  ; nctsj2(1)   = 274 - isrow
            !                                           ! Lake Superior
            nctsi1(2)   = 198  ; nctsj1(2)   = 258 - isrow
            nctsi2(2)   = 204  ; nctsj2(2)   = 262 - isrow
            !                                           ! Lake Michigan
            nctsi1(3)   = 201  ; nctsj1(3)   = 250 - isrow
            nctsi2(3)   = 203  ; nctsj2(3)   = 256 - isrow
            !                                           ! Lake Huron
            nctsi1(4)   = 204  ; nctsj1(4)   = 252 - isrow
            nctsi2(4)   = 209  ; nctsj2(4)   = 256 - isrow
            !                                           ! Lake Erie
            nctsi1(5)   = 206  ; nctsj1(5)   = 249 - isrow
            nctsi2(5)   = 209  ; nctsj2(5)   = 251 - isrow
            !                                           ! Lake Ontario
            nctsi1(6)   = 210  ; nctsj1(6)   = 252 - isrow
            nctsi2(6)   = 212  ; nctsj2(6)   = 252 - isrow
            !                                           ! Victoria Lake
            nctsi1(7)   = 321  ; nctsj1(7)   = 180 - isrow
            nctsi2(7)   = 322  ; nctsj2(7)   = 189 - isrow
            !                                           ! Baltic Sea
            nctsi1(8)   = 297  ; nctsj1(8)   = 270 - isrow
            nctsi2(8)   = 308  ; nctsj2(8)   = 293 - isrow
            !                                        
            !                                           ! =======================
            CASE ( 2 )                                  !  ORCA_R2 configuration
               !                                        ! =======================
               !                                            ! Caspian Sea
               nctsi1(1)   =  11  ;  nctsj1(1)   = 103
               nctsi2(1)   =  17  ;  nctsj2(1)   = 112
               !                                            ! Great North American Lakes
               nctsi1(2)   =  97  ;  nctsj1(2)   = 107
               nctsi2(2)   = 103  ;  nctsj2(2)   = 111
               !                                            ! Black Sea 1 : west part of the Black Sea
               nctsi1(3)   = 174  ;  nctsj1(3)   = 107
               nctsi2(3)   = 181  ;  nctsj2(3)   = 112
              !                                            ! Black Sea 2 : est part of the Black Sea
               nctsi1(4)   =   2  ;  nctsj1(4)   = 107
               nctsi2(4)   =   6  ;  nctsj2(4)   = 112
               !                                            ! Baltic Sea
               nctsi1(5)   =  145 ;  nctsj1(5)   = 116
               nctsi2(5)   =  150 ;  nctsj2(5)   = 126
               !                                        ! =======================
            CASE ( 4 )                                  !  ORCA_R4 configuration
               !                                        ! =======================
               !                                            ! Caspian Sea
               nctsi1(1)   =  4  ;  nctsj1(1)   = 53
               nctsi2(1)   =  4  ;  nctsj2(1)   = 56
               !                                            ! Great North American Lakes
               nctsi1(2)   = 49  ;  nctsj1(2)   = 55
               nctsi2(2)   = 51  ;  nctsj2(2)   = 56
               !                                            ! Black Sea
               nctsi1(3)   = 88  ;  nctsj1(3)   = 55
               nctsi2(3)   = 91  ;  nctsj2(3)   = 56
               !                                            ! Baltic Sea
               nctsi1(4)   = 75  ;  nctsj1(4)   = 59
               nctsi2(4)   = 76  ;  nctsj2(4)   = 61
               !                                        ! =======================
            CASE ( 025 )                                ! ORCA_R025 configuration
               !                                        ! =======================
                                                     ! Caspian + Aral sea
               nctsi1(1)   = 1330 ; nctsj1(1)   = 645
               nctsi2(1)   = 1400 ; nctsj2(1)   = 795
               !                                        ! Azov Sea
               nctsi1(2)   = 1284 ; nctsj1(2)   = 722
               nctsi2(2)   = 1304 ; nctsj2(2)   = 747
               !
            END SELECT
            !
         ENDIF
         !

         ! convert the position in local domain indices
         ! --------------------------------------------
         DO jc = 1, npncts
            nctsi1(jc)   = mi0( nctsi1(jc) )
            nctsj1(jc)   = mj0( nctsj1(jc) )

            nctsi2(jc)   = mi1( nctsi2(jc) )
            nctsj2(jc)   = mj1( nctsj2(jc) )
         END DO
         !
      ENDIF

      ! Restore close seas values to initial data
      IF( ln_trcdta .AND. nb_trcdta > 0 )  THEN   ! Initialisation of tracer from a file that may also be used for damping
         !
         IF(lwp)  WRITE(numout,*)
         IF(lwp)  WRITE(numout,*) ' trc_dmp_clo : Restoring of nutrients on close seas at time-step kt = ', kt
         IF(lwp)  WRITE(numout,*)
         !
         CALL wrk_alloc( jpi, jpj, jpk, ztrcdta )   ! Memory allocation
         !
         DO jn = 1, jptra
            IF( ln_trc_ini(jn) ) THEN      ! update passive tracers arrays with input data read from file
                jl = n_trc_index(jn)
                CALL trc_dta( kt, sf_trcdta(jl), rf_trfac(jl), ztrcdta )   ! read tracer data at nit000
                DO jc = 1, npncts
                   DO jk = 1, jpkm1
                      DO jj = nctsj1(jc), nctsj2(jc)
                         DO ji = nctsi1(jc), nctsi2(jc)
                            trn(ji,jj,jk,jn) = ztrcdta(ji,jj,jk)
                            trb(ji,jj,jk,jn) = trn(ji,jj,jk,jn)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          CALL wrk_dealloc( jpi, jpj, jpk, ztrcdta )
      ENDIF
      !
   END SUBROUTINE trc_dmp_clo


   SUBROUTINE trc_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the nammbf namelist and check the parameters
      !!              called by trc_dmp at the first timestep (nittrc000)
      !!----------------------------------------------------------------------
      !
      INTEGER :: imask  !local file handle

      IF( nn_timing == 1 )  CALL timing_start('trc_dmp_init')
      !
      !Allocate arrays
      IF( trc_dmp_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_dmp_init: unable to allocate arrays' )

      IF( lzoom )   nn_zdmp_tr = 0           ! restoring to climatology at closed north or south boundaries
      SELECT CASE ( nn_zdmp_tr )
      CASE ( 0 )   ;   IF(lwp) WRITE(numout,*) '   tracer damping throughout the water column'
      CASE ( 1 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the turbocline (avt > 5 cm2/s)'
      CASE ( 2 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixed layer'
      CASE DEFAULT
         WRITE(ctmp1,*) 'bad flag value for nn_zdmp_tr = ', nn_zdmp_tr
         CALL ctl_stop(ctmp1)
      END SELECT

      IF( .NOT. ln_tradmp )   &
         &   CALL ctl_stop( 'passive trace damping need key_tradmp to compute damping coef.' )
      !
      !                          ! Read damping coefficients from file
      !Read in mask from file
      CALL iom_open ( cn_resto_tr, imask)
      CALL iom_get  ( imask, jpdom_autoglo, 'resto', restotr)
      CALL iom_close( imask )
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_dmp_init')
      !
   END SUBROUTINE trc_dmp_init

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_dmp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dmp
#endif


   !!======================================================================
END MODULE trcdmp
