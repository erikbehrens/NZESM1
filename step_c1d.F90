MODULE step_c1d
   !!======================================================================
   !!                       ***  MODULE step_c1d  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping - c1d case
   !!======================================================================
   !! History :   2.0  !  2004-04  (C. Ethe)  adapted from step.F90 for C1D
   !!             3.0  !  2008-04  (G. Madec)  redo the adaptation to include SBC
   !!----------------------------------------------------------------------
#if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'                                       1D Configuration
   !!----------------------------------------------------------------------  
   !!   stp_c1d        : NEMO system time-stepping in c1d case
   !!----------------------------------------------------------------------
   USE step_oce        ! time stepping definition modules 
#if defined key_top
   USE trcstp          ! passive tracer time-stepping      (trc_stp routine)
#endif
   USE dyncor_c1d      ! Coriolis term (c1d case)         (dyn_cor_1d     )
   USE dynnxt_c1d      ! time-stepping                    (dyn_nxt routine)
   USE dyndmp          ! U & V momentum damping           (dyn_dmp routine)
   USE restart         ! restart 

   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_c1d      ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/C1D 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_c1d( kstp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp_c1d  ***
      !!                      
      !! ** Purpose :  - Time stepping of SBC including LIM (dynamic and thermodynamic eqs.)
      !!               - Time stepping of OPA (momentum and active tracer eqs.)
      !!               - Time stepping of TOP (passive tracer eqs.)
      !! 
      !! ** Method  : -1- Update forcings and data  
      !!              -2- Update vertical ocean physics 
      !!              -3- Compute the t and s trends 
      !!              -4- Update t and s 
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
      INTEGER ::   jk       ! dummy loop indice
      INTEGER ::   indic    ! error indicator if < 0
      !! ---------------------------------------------------------------------

                             indic = 0                ! reset to no error condition
      IF( kstp == nit000 )   CALL iom_init( "nemo")   ! iom_put initialization (must be done after nemo_init for AGRIF+XIOS+OASIS)
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1, "nemo" )   ! say to iom that we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update data, open boundaries, surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL sbc    ( kstp )         ! Sea Boundary Condition (including sea-ice)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update                (ua, va, ta, sa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL eos_rab( tsb, rab_b )   ! before local thermal/haline expension ratio at T-points
                         CALL eos_rab( tsn, rab_n )   ! now    local thermal/haline expension ratio at T-points
                         CALL bn2( tsb, rab_b, rn2b ) ! before Brunt-Vaisala frequency
                         CALL bn2( tsn, rab_n, rn2  ) ! now    Brunt-Vaisala frequency
      !  VERTICAL PHYSICS   
                         CALL zdf_bfr( kstp )         ! bottom friction
      !                                               ! Vertical eddy viscosity and diffusivity coefficients
      IF( lk_zdfric  )   CALL zdf_ric( kstp )            ! Richardson number dependent Kz
      IF( lk_zdftke  )   CALL zdf_tke( kstp )            ! TKE closure scheme for Kz
      IF( lk_zdfgls  )   CALL zdf_gls( kstp )            ! GLS closure scheme for Kz
      IF( lk_zdfkpp  )   CALL zdf_kpp( kstp )            ! KPP closure scheme for Kz
      IF( lk_zdfcst  )   THEN                            ! Constant Kz (reset avt, avm[uv] to the background value)
         avt (:,:,:) = rn_avt0 * tmask(:,:,:)
         avmu(:,:,:) = rn_avm0 * umask(:,:,:)
         avmv(:,:,:) = rn_avm0 * vmask(:,:,:)
      ENDIF

      IF( ln_rnf_mouth ) THEN                         ! increase diffusivity at rivers mouths
         DO jk = 2, nkrnf   ;   avt(:,:,jk) = avt(:,:,jk) + 2.e0 * rn_avt_rnf * rnfmsk(:,:)   ;   END DO
      ENDIF
      IF( ln_zdfevd  )   CALL zdf_evd( kstp )         ! enhanced vertical eddy diffusivity

      IF( lk_zdftmx  )   CALL zdf_tmx( kstp )         ! tidal vertical mixing

      IF( lk_zdfddm .AND. .NOT. lk_zdfkpp )   &
         &               CALL zdf_ddm( kstp )         ! double diffusive mixing
         
                         CALL zdf_mxl( kstp )         ! mixed layer depth

                                                      ! write tke information in the restart file
      IF( lrst_oce .AND. lk_zdftke )   CALL tke_rst( kstp, 'WRITE' )
                                                      ! write gls information in the restart file
      IF( lrst_oce .AND. lk_zdfgls )   CALL gls_rst( kstp, 'WRITE' )

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs             (ua, va, ta, sa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL dia_wri( kstp )       ! ocean model: outputs
      IF( lk_diahth  )   CALL dia_hth( kstp )       ! Thermocline depth (20°C)


#if defined key_top
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        CALL trc_stp( kstp )       ! time-stepping
#endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Active tracers                              (ua, va used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        tsa(:,:,:,:) = 0._wp       ! set tracer trends to zero

                        CALL tra_sbc( kstp )       ! surface boundary condition
      IF( ln_traqsr )   CALL tra_qsr( kstp )       ! penetrative solar radiation qsr
      IF( ln_tradmp )   CALL tra_dmp( kstp )       ! internal damping trends- tracers
      IF( lk_zdfkpp )   CALL tra_kpp( kstp )       ! KPP non-local tracer fluxes
                        CALL tra_zdf( kstp )       ! vertical mixing
                        CALL eos( tsn, rhd, rhop, gdept_0(:,:,:) )   ! now potential density for zdfmxl
      IF( ln_zdfnpc )   CALL tra_npc( kstp )       ! applied non penetrative convective adjustment on (t,s)
                        CALL tra_nxt( kstp )       ! tracer fields at next time step



      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Dynamics                                    (ta, sa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        ua(:,:,:) = 0._wp          ! set dynamics trends to zero
                        va(:,:,:) = 0._wp

      IF( ln_dyndmp )   CALL dyn_dmp    ( kstp )   ! internal damping trends- momentum
                        CALL dyn_cor_c1d( kstp )   ! vorticity term including Coriolis
                        CALL dyn_zdf    ( kstp )   ! vertical diffusion
                        CALL dyn_nxt_c1d( kstp )   ! lateral velocity at next time step

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control and restarts
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             CALL stp_ctl( kstp, indic )
      IF( kstp == nit000 )   CALL iom_close( numror )      ! close input  ocean restart file
      IF( lrst_oce       )   CALL rst_write( kstp )        ! write output ocean restart file
      !
#if defined key_iomput
      IF( kstp == nitend .OR. indic < 0 )   CALL xios_context_finalize()   ! needed for XIOS
      !
#endif
   END SUBROUTINE stp_c1d

#else
   !!----------------------------------------------------------------------
   !!   Default key                                            NO 1D Config
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE stp_c1d ( kt )      ! dummy routine
      WRITE(*,*) 'stp_c1d: You should not have seen this print! error?', kt
   END SUBROUTINE stp_c1d
#endif

   !!======================================================================
END MODULE step_c1d
