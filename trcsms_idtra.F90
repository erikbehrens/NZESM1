MODULE trcsms_idtra
   !!======================================================================
   !!                      ***  MODULE trcsms_idtra  ***
   !! TOP : TRI main model
   !!======================================================================
   !! History :    -   !  1999-10  (JC. Dutay)  original code
   !!             1.0  !  2004-03 (C. Ethe) free form + modularity
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  reorganisation
   !!----------------------------------------------------------------------
#if defined key_idtra
   !!----------------------------------------------------------------------
   !!   'key_idtra'                                               TRI tracers
   !!----------------------------------------------------------------------
   !!   trc_sms_idtra     :  compute and add TRI suface forcing to TRI trends
   !!   trc_idtra_cst :  sets constants for TRI surface forcing computation
   !!----------------------------------------------------------------------
   USE oce_trc      ! Ocean variables
   USE par_trc      ! TOP parameters
   USE trc          ! TOP variables
   USE trd_oce
   USE trdtrc
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_idtra        ! called in ???
   PUBLIC   trc_sms_idtra_alloc  ! called in ???
   !
   INTEGER , PUBLIC    ::   nyear_res      ! restoring time constant (year)
   INTEGER , PUBLIC    ::   numnatm
   REAL(wp), PUBLIC    ::   FDEC
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qtr_idtra  ! flux at surface
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qint_idtra ! cumulative flux 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   inv_idtra  ! vertic. inventory

   !                          ! coefficients for conversion
   REAL(wp) ::  WTEMP


   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_sms_idtra( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_idtra  ***
      !!
      !! ** Purpose :   Compute the surface boundary contition on TRI 11
      !!             passive tracer associated with air-mer fluxes and add it
      !!             to the general trend of tracers equations.
      !!
      !! ** Method  : - get the atmospheric partial pressure - given in pico -
      !!              - computation of solubility ( in 1.e-12 mol/l then in 1.e-9 mol/m3)
      !!              - computation of transfert speed ( given in cm/hour ----> cm/s )
      !!              - the input function is given by :
      !!                speed * ( concentration at equilibrium - concentration at surface )
      !!              - the input function is in pico-mol/m3/s and the
      !!                TRI concentration in pico-mol/m3
      !!               
      !! *** For Idealized Tracers              
      !!              - no need for any temporal references, 
      !!              nor any atmospheric concentration, nor air -sea fluxes
      !!              - Here we fixe surface concentration to 1.0 Tracer-Unit/m3
      !!              - Then we add a decay (radioactive-like) to this tracer concentration 
      !!              - the Half life deccay is chosen by the user, depending of the experiment. 
      !!              
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt    ! ocean time-step index
      !!
      INTEGER                ::   ji, jj, jn, jl, jk
      REAL(wp)               ::   rlx                 !! relaxation time (1 day)
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_idtra')
      !
      rlx = 10./(60. * 60. * 24.)                              !! relaxation time (1/10 day)
      IF (kt == nittrc000) THEN
         IF(lwp) WRITE(numout,*) '   trcsms_idtra :'
         IF(lwp) WRITE(numout,*) '   ~~~~~~~~~~~~~~~~~'
         IF(lwp) WRITE(numout,*) '   - idtra decay factor : ', FDEC
         IF(lwp) WRITE(numout,*) '   - relaxation time    : ', rlx
# if defined key_debug_medusa
         CALL flush(numout)
# endif
      !   CALL idtra_init
      ENDIF

         !
      inv_idtra(:,:,:) = 0.0                                   !! init the inventory
      qtr_idtra(:,:,:) = 0.0                                   !! init the air-sea flux
      DO jl = 1, jp_idtra
         jn = jp_idtra0 + jl - 1

      !!   DO jj = 1, jpj
      !!      DO ji = 1, jpi
           DO jj = 2,jpjm1
              DO ji = 2,jpim1

         !! First, a crude version. will be much inproved later.
             qtr_idtra(ji,jj,jl)  = rlx * (1. - trb(ji,jj,1,jn)) * tmask(ji,jj,1) *   & 
                                  fse3t(ji,jj,1)                  !! Air-sea Flux

           !! DEBUG-TEST : Set flux equal to 0, see if it induces the pb we see in the MED  
           !!  qtr_idtra(ji,jj,jl)  = 0.0
           ENDDO
         ENDDO
         tra(:,:,1,jn)      = tra(:,:,1,jn) + ( qtr_idtra(:,:,jl) *  &
                            tmask(:,:,1) / fse3t(:,:,1) )
         qint_idtra(:,:,jl) = qint_idtra(:,:,jl) +                   &           
                              qtr_idtra(:,:,jl) * rdt              !! Cumulative Air-sea Flux


         DO jk =1,jpk
            inv_idtra(:,:,jl) = inv_idtra(:,:,jl) +                  &
                     (trn(:,:,jk,jn) * fse3t(:,:,jk) * tmask(:,:,jk))  !! vertical inventory
         ENDDO
!
!DECAY of OUR IDEALIZED TRACER
! ---------------------------------------

         DO  jk =1,jpk
      !!      DO jj=1,jpj
      !!        DO  ji =1,jpi
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
            
                 !! IF (trn(ji,jj,jk,jn) > 0.0) THEN
                    WTEMP = trn(ji,jj,jk,jn) * (1. - FDEC )
                    tra(ji,jj,jk,jn) = (tra(ji,jj,jk,jn) - WTEMP/rdt ) * &
                                     tmask(ji,jj,jk)
                 !! ENDIF 
              ENDDO 
            ENDDO
         ENDDO

      ENDDO
    !! jn loop
!
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) '   IDTRA - calculation part - DONE trc_sms_idtra -- '
      CALL flush(numout)
# endif
        !
        !! restart and diagnostics management -- 
      !IF( lrst_trc ) THEN
      !   IF(lwp) WRITE(numout,*)
      !   IF(lwp) WRITE(numout,*) 'trc_sms_idtra : cumulated input function fields written in ocean restart file ',   &
      !      &                    'at it= ', kt,' date= ', ndastp
      !   IF(lwp) WRITE(numout,*) '~~~~'
      !   !!DO jn = jp_idtra0, jp_idtra1
      !      CALL iom_rstput( kt, nitrst, numrtw, 'qint_IDTRA', qint_idtra(:,:,1) )
      !   !!END DO
 ! if defined key_debug_medusa
      !   IF(lwp) WRITE(numout,*) '   IDTRA - writing diag-restart - DONE trc_sms_idtra -- '
      !   CALL flush(numout)
 ! endif
      !ENDIF
      !
         CALL iom_put( "qtrIDTRA"  , qtr_idtra (:,:,1) )
         CALL iom_put( "qintIDTRA" , qint_idtra(:,:,1) )
         CALL iom_put( "invIDTRA" , inv_idtra(:,:,1) )
      !
# if defined key_debug_medusa
      IF(lwp) WRITE(numout,*) '   IDTRA - writing diag - DONE trc_sms_idtra -- '
      CALL flush(numout)
# endif
      !
      IF( l_trdtrc ) THEN
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) '   IDTRA - writing trends - trc_sms_idtra -- '
         CALL flush(numout)
# endif
          DO jn = jp_idtra0, jp_idtra1
            CALL trd_trc( tra(:,:,:,jn), jn, jptra_sms, kt )   ! save trends
          END DO
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) '   IDTRA - writing trends - DONE trc_sms_idtra -- '
         CALL flush(numout)
# endif
      END IF
      !
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) '   IDTRA - Check: nn_timing = ', nn_timing 
         CALL flush(numout)
# endif
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_idtra')
      !
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) '   IDTRA DONE trc_sms_idtra -- '
      CALL flush(numout)
# endif
      !
   END SUBROUTINE trc_sms_idtra

   SUBROUTINE idtra_init
      !!---------------------------------------------------------------------
      !!                     ***  idtra_init  ***
      !!
      !! ** Purpose : read restart values for IDTRA model
      !!---------------------------------------------------------------------
      INTEGER :: jn

      IF( ln_rsttr ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' Read specific variables from Ideal Tracers model '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
         !
         DO jn = jp_idtra0, jp_idtra1
            CALL iom_get( numrtr, jpdom_autoglo, 'qint_IDTRA', qint_idtra(:,:,jn) )
         END DO
      ENDIF
      IF(lwp) WRITE(numout,*) 'idtra restart variables read -- OK'
      !
   END SUBROUTINE idtra_init

   INTEGER FUNCTION trc_sms_idtra_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_idtra_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( qtr_idtra (jpi,jpj,jp_idtra) ,     &
         &      inv_idtra(jpi,jpj,jp_idtra)  ,     &
         &      qint_idtra(jpi,jpj,jp_idtra) , STAT=trc_sms_idtra_alloc )
         !
      IF( trc_sms_idtra_alloc /= 0 ) CALL ctl_warn('trc_sms_idtra_alloc : failed to allocate arrays.')
      !
   END FUNCTION trc_sms_idtra_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                         No TRI tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_idtra( kt )       ! Empty routine
      WRITE(*,*) 'trc_sms_idtra: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_idtra
#endif

   !!======================================================================
END MODULE trcsms_idtra




