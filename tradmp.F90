MODULE tradmp
   !!======================================================================
   !!                       ***  MODULE  tradmp  ***
   !! Ocean physics: internal restoring trend on active tracers (T and S)
   !!======================================================================
   !! History :  OPA  ! 1991-03  (O. Marti, G. Madec)  Original code
   !!                 ! 1992-06  (M. Imbard)  doctor norme
   !!                 ! 1996-01  (G. Madec)  statement function for e3
   !!                 ! 1997-05  (G. Madec)  macro-tasked on jk-slab
   !!                 ! 1998-07  (M. Imbard, G. Madec) ORCA version
   !!            7.0  ! 2001-02  (M. Imbard)  cofdis, Original code
   !!            8.1  ! 2001-02  (G. Madec, E. Durand)  cleaning
   !!  NEMO      1.0  ! 2002-08  (G. Madec, E. Durand)  free form + modules
   !!            3.2  ! 2009-08  (G. Madec, C. Talandier)  DOCTOR norm for namelist parameter
   !!            3.3  ! 2010-06  (C. Ethe, G. Madec) merge TRA-TRC 
   !!            3.4  ! 2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal + suppression of CPP keys
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_dmp_alloc : allocate tradmp arrays
   !!   tra_dmp       : update the tracer trend with the internal damping
   !!   tra_dmp_init  : initialization, namlist read, parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean: variables
   USE dom_oce        ! ocean: domain variables
   USE c1d            ! 1D vertical configuration
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE zdf_oce        ! ocean: vertical physics
   USE phycst         ! physical constants
   USE dtatsd         ! data: temperature & salinity
   USE zdfmxl         ! vertical physics: mixed layer depth
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE wrk_nemo       ! Memory allocation
   USE timing         ! Timing
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_dmp      ! routine called by step.F90
   PUBLIC   tra_dmp_init ! routine called by opa.F90

   !                               !!* Namelist namtra_dmp : T & S newtonian damping *
   ! nn_zdmp and cn_resto are public as they are used by C1D/dyndmp.F90
   LOGICAL , PUBLIC ::   ln_tradmp   !: internal damping flag
   INTEGER , PUBLIC ::   nn_zdmp     ! = 0/1/2 flag for damping in the mixed layer
   CHARACTER(LEN=200) , PUBLIC :: cn_resto      ! name of netcdf file containing restoration coefficient field
   !


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   strdmp   !: damping salinity trend (psu/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ttrdmp   !: damping temperature trend (Celcius/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   resto    !: restoring coeff. on T and S (s-1)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_dmp_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( strdmp(jpi,jpj,jpk) , ttrdmp(jpi,jpj,jpk), resto(jpi,jpj,jpk), STAT= tra_dmp_alloc )
      !
      IF( lk_mpp            )   CALL mpp_sum ( tra_dmp_alloc )
      IF( tra_dmp_alloc > 0 )   CALL ctl_warn('tra_dmp_alloc: allocation of arrays failed')
      !
   END FUNCTION tra_dmp_alloc


   SUBROUTINE tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp  ***
      !!                  
      !! ** Purpose :   Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - (ta,sa)   tracer trends updated with the damping trend
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zta, zsa             ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  zts_dta 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_dmp')
      !
      CALL wrk_alloc( jpi, jpj, jpk, jpts,  zts_dta )
      !
      !                           !==   input T-S data at kt   ==!
      CALL dta_tsd( kt, zts_dta )            ! read and interpolates T-S data at kt
      !
      SELECT CASE ( nn_zdmp )     !==    type of damping   ==!
      !
      CASE( 0 )                   !==  newtonian damping throughout the water column  ==!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zta = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                  zsa = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + zta
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + zsa
                  strdmp(ji,jj,jk) = zsa           ! save the trend (used in asmtrj)
                  ttrdmp(ji,jj,jk) = zta      
               END DO
            END DO
         END DO
         !
      CASE ( 1 )                  !==  no damping in the turbocline (avt > 5 cm2/s)  ==!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( avt(ji,jj,jk) <= 5.e-4_wp ) THEN
                     zta = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                     zsa = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  ELSE
                     zta = 0._wp
                     zsa = 0._wp  
                  ENDIF
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + zta
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + zsa
                  strdmp(ji,jj,jk) = zsa           ! save the salinity trend (used in asmtrj)
                  ttrdmp(ji,jj,jk) = zta
               END DO
            END DO
         END DO
         !
      CASE ( 2 )                  !==  no damping in the mixed layer   ==!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fsdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     zta = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                     zsa = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  ELSE
                     zta = 0._wp
                     zsa = 0._wp  
                  ENDIF
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + zta
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + zsa
                  strdmp(ji,jj,jk) = zsa           ! save the salinity trend (used in asmtrj)
                  ttrdmp(ji,jj,jk) = zta
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      IF( l_trdtra )   THEN       ! trend diagnostic
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_dmp, ttrdmp )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_dmp, strdmp )
      ENDIF
      !                           ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' dmp  - Ta: ', mask1=tmask,   &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      CALL wrk_dealloc( jpi, jpj, jpk, jpts,  zts_dta )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_dmp')
      !
   END SUBROUTINE tra_dmp


   SUBROUTINE tra_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the namtra_dmp namelist and check the parameters
      !!----------------------------------------------------------------------
      NAMELIST/namtra_dmp/ ln_tradmp, nn_zdmp, cn_resto
      INTEGER ::  ios         ! Local integer for output status of namelist read
      INTEGER :: imask        ! File handle 
      !!
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )   ! Namelist namtra_dmp in reference namelist : T & S relaxation
      READ  ( numnam_ref, namtra_dmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_dmp in reference namelist', lwp )
      !
      REWIND( numnam_cfg )   ! Namelist namtra_dmp in configuration namelist : T & S relaxation
      READ  ( numnam_cfg, namtra_dmp, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_dmp in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_dmp )

      IF(lwp) THEN                 !Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_dmp_init : T and S newtonian relaxation'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_dmp : set relaxation parameters'
         WRITE(numout,*) '      Apply relaxation   or not       ln_tradmp = ', ln_tradmp
         WRITE(numout,*) '      mixed layer damping option      nn_zdmp   = ', nn_zdmp
         WRITE(numout,*) '      Damping file name               cn_resto  = ', cn_resto
         WRITE(numout,*)
      ENDIF

      IF( ln_tradmp) THEN
         !
         !Allocate arrays
         IF( tra_dmp_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'tra_dmp_init: unable to allocate arrays' )

         !Check values of nn_zdmp
         SELECT CASE (nn_zdmp)
         CASE ( 0 )  ; IF(lwp) WRITE(numout,*) '   tracer damping as specified by mask'
         CASE ( 1 )  ; IF(lwp) WRITE(numout,*) '   no tracer damping in the turbocline'
         CASE ( 2 )  ; IF(lwp) WRITE(numout,*) '   no tracer damping in the mixed layer'
         END SELECT

         !TG: Initialisation of dtatsd - Would it be better to have dmpdta routine
         !so can damp to something other than intitial conditions files?
         IF( .NOT.ln_tsd_tradmp ) THEN
            CALL ctl_warn( 'tra_dmp_init: read T-S data not initialized, we force ln_tsd_tradmp=T' )
            CALL dta_tsd_init( ld_tradmp=ln_tradmp )        ! forces the initialisation of T-S data
         ENDIF

         !initialise arrays - Are these actually used anywhere else?
         strdmp(:,:,:) = 0._wp
         ttrdmp(:,:,:) = 0._wp

         !Read in mask from file
         CALL iom_open ( cn_resto, imask)
         CALL iom_get  ( imask, jpdom_autoglo, 'resto', resto)
         CALL iom_close( imask )
       ENDIF

   END SUBROUTINE tra_dmp_init

END MODULE tradmp
