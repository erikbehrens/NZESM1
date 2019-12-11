MODULE trctrp
   !!======================================================================
   !!                       ***  MODULE trctrp  ***
   !! Ocean Physics    : manage the passive tracer transport
   !!======================================================================
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             3.3  !  2010-07 (C. Ethe) Merge TRA-TRC
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_trp        : passive tracer transport
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean passive tracers variables 
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE trabbl          ! bottom boundary layer               (trc_bbl routine)
   USE trcbbl          ! bottom boundary layer               (trc_bbl routine)
   USE zdfkpp          ! KPP non-local tracer fluxes         (trc_kpp routine)
   USE trcdmp          ! internal damping                    (trc_dmp routine)
   USE trcldf          ! lateral mixing                      (trc_ldf routine)
   USE trcadv          ! advection                           (trc_adv routine)
   USE trczdf          ! vertical diffusion                  (trc_zdf routine)
   USE trcnxt          ! time-stepping                       (trc_nxt routine)
   USE trcrad          ! positivity                          (trc_rad routine)
   USE trcsbc          ! surface boundary condition          (trc_sbc routine)
   USE zpshde          ! partial step: hor. derivative       (zps_hde routine)
# if defined key_debug_medusa
   USE trcstat
# endif


#if defined key_agrif
   USE agrif_top_sponge ! tracers sponges
   USE agrif_top_update ! tracers updates
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_trp    ! called by trc_stp

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_trp( kstp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_trp  ***
      !!                      
      !! ** Purpose :   Management of passive tracers transport
      !! 
      !! ** Method  : - Compute the passive tracers trends 
      !!              - Update the passive tracers
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kstp  ! ocean time-step index
      !! ---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_trp')
      !
      IF( .NOT. lk_c1d ) THEN
         !
                                CALL trc_sbc( kstp )            ! surface boundary condition
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) ' MEDUSA trc_trp after trc_sbc at kt =', kstp
         CALL trc_rst_tra_stat
         CALL flush(numout)
# endif
         IF( lk_trabbl )        CALL trc_bbl( kstp )            ! advective (and/or diffusive) bottom boundary layer scheme
         IF( ln_trcdmp )        CALL trc_dmp( kstp )            ! internal damping trends
                                CALL trc_adv( kstp )            ! horizontal & vertical advection 
         IF( ln_zps ) THEN
           IF( ln_isfcav ) THEN ; CALL zps_hde_isf( kstp, jptra, trb, pgtu=gtru, pgtv=gtrv, pgtui=gtrui, pgtvi=gtrvi )  ! both top & bottom
           ELSE                 ; CALL zps_hde    ( kstp, jptra, trb, gtru, gtrv )                                      !  only bottom
           ENDIF
         ENDIF
                                CALL trc_ldf( kstp )            ! lateral mixing
         IF( .NOT. lk_offline .AND. lk_zdfkpp )    &
            &                   CALL trc_kpp( kstp )            ! KPP non-local tracer fluxes
#if defined key_agrif
         IF(.NOT. Agrif_Root()) CALL Agrif_Sponge_trc           ! tracers sponge
#endif
                                CALL trc_zdf( kstp )            ! vertical mixing and after tracer fields
                                CALL trc_nxt( kstp )            ! tracer fields at next time step     
# if defined key_debug_medusa
         IF(lwp) WRITE(numout,*) ' MEDUSA trc_trp after trc_nxt at kt =', kstp
         CALL trc_rst_tra_stat
         CALL flush(numout)
# endif
         IF( ln_trcrad )        CALL trc_rad( kstp )            ! Correct artificial negative concentrations
         IF( ln_trcdmp_clo )    CALL trc_dmp_clo( kstp )        ! internal damping trends on closed seas only

#if defined key_agrif
      IF( .NOT. Agrif_Root())   CALL Agrif_Update_Trc( kstp )   ! Update tracer at AGRIF zoom boundaries : children only
#endif

      ELSE                                               ! 1D vertical configuration
                                CALL trc_sbc( kstp )            ! surface boundary condition
         IF( .NOT. lk_offline .AND. lk_zdfkpp )    &
            &                   CALL trc_kpp( kstp )            ! KPP non-local tracer fluxes
                                CALL trc_zdf( kstp )            ! vertical mixing and after tracer fields
                                CALL trc_nxt( kstp )            ! tracer fields at next time step     
          IF( ln_trcrad )       CALL trc_rad( kstp )            ! Correct artificial negative concentrations
         !
      END IF
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_trp')
      !
9400  FORMAT(a25,i4,D23.16)
      !
   END SUBROUTINE trc_trp

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        No TOP models
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_trp( kstp )              ! Empty routine
      INTEGER, INTENT(in) ::   kstp
      WRITE(*,*) 'trc_trp: You should not have seen this print! error?', kstp
   END SUBROUTINE trc_trp
#endif
   
   !!======================================================================
END MODULE trctrp
