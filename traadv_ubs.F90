MODULE traadv_ubs
   !!==============================================================================
   !!                       ***  MODULE  traadv_ubs  ***
   !! Ocean active tracers:  horizontal & vertical advective trend
   !!==============================================================================
   !! History :  1.0  !  2006-08  (L. Debreu, R. Benshila)  Original code
   !!            3.3  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_ubs : update the tracer trend with the horizontal
   !!                 advection trends using a third order biaised scheme  
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE trc_oce        ! share passive tracers/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE dynspg_oce     ! choice/control of key cpp for surface pressure gradient
   USE diaptr         ! poleward transport diagnostics
   !
   USE lib_mpp        ! I/O library
   USE lbclnk         ! ocean lateral boundary condition (or mpp link)
   USE in_out_manager ! I/O manager
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_ubs   ! routine called by traadv module

   LOGICAL :: l_trd  ! flag to compute trends or not

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_ubs ( kt, kit000, cdtype, p2dt, pun, pvn, pwn,      &
      &                                       ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_ubs  ***
      !!                 
      !! ** Purpose :   Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method  :   The upstream biased scheme (UBS) is based on a 3rd order
      !!      upstream-biased parabolic interpolation (Shchepetkin and McWilliams 2005)
      !!      It is only used in the horizontal direction.
      !!      For example the i-component of the advective fluxes are given by :
      !!                !  e2u e3u un ( mi(Tn) - zltu(i  ) )   if un(i) >= 0
      !!          ztu = !  or 
      !!                !  e2u e3u un ( mi(Tn) - zltu(i+1) )   if un(i) < 0
      !!      where zltu is the second derivative of the before temperature field:
      !!          zltu = 1/e3t di[ e2u e3u / e1u di[Tb] ]
      !!      This results in a dissipatively dominant (i.e. hyper-diffusive) 
      !!      truncation error. The overall performance of the advection scheme 
      !!      is similar to that reported in (Farrow and Stevens, 1995). 
      !!      For stability reasons, the first term of the fluxes which corresponds
      !!      to a second order centered scheme is evaluated using the now velocity 
      !!      (centered in time) while the second term which is the diffusive part 
      !!      of the scheme, is evaluated using the before velocity (forward in time). 
      !!      Note that UBS is not positive. Do not use it on passive tracers.
      !!                On the vertical, the advection is evaluated using a TVD scheme,
      !!      as the UBS have been found to be too diffusive.
      !!
      !! ** Action : - update (pta) with the now advective tracer trends
      !!
      !! Reference : Shchepetkin, A. F., J. C. McWilliams, 2005, Ocean Modelling, 9, 347-404. 
      !!             Farrow, D.E., Stevens, D.P., 1995, J. Phys. Ocean. 25, 1731Ğ1741. 
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt            ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean transport components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   ztra, zbtr, zcoef, z2dtt                       ! local scalars
      REAL(wp) ::   zfp_ui, zfm_ui, zcenut, ztak, zfp_wk, zfm_wk   !   -      -
      REAL(wp) ::   zfp_vj, zfm_vj, zcenvt, zeeu, zeev, z_hdivn    !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztu, ztv, zltu, zltv, zti, ztw
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_ubs')
      !
      CALL wrk_alloc( jpi, jpj, jpk, ztu, ztv, zltu, zltv, zti, ztw )
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_ubs :  horizontal UBS advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      l_trd = .FALSE.
      IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) ) l_trd = .TRUE.
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         ! 1. Bottom value : flux set to zero
         ! ----------------------------------
         zltu(:,:,jpk) = 0.e0       ;      zltv(:,:,jpk) = 0.e0
         !                                              
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                   
            !  Laplacian
            DO jj = 1, jpjm1            ! First derivative (gradient)
               DO ji = 1, fs_jpim1   ! vector opt.
                  zeeu = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj) * umask(ji,jj,jk)
                  zeev = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj) * vmask(ji,jj,jk)
                  ztu(ji,jj,jk) = zeeu * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zeev * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
            END DO
            DO jj = 2, jpjm1            ! Second derivative (divergence)
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zcoef = 1. / ( 6. * fse3t(ji,jj,jk) )
                  zltu(ji,jj,jk) = (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)  ) * zcoef
                  zltv(ji,jj,jk) = (  ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  ) * zcoef
               END DO
            END DO
            !                                    
         END DO                                           ! End of slab         
         CALL lbc_lnk( zltu, 'T', 1. )   ;    CALL lbc_lnk( zltv, 'T', 1. )   ! Lateral boundary cond. (unchanged sgn)

         !    
         !  Horizontal advective fluxes               
         DO jk = 1, jpkm1                                 ! Horizontal slab
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! upstream transport (x2)
                  zfp_ui = pun(ji,jj,jk) + ABS( pun(ji,jj,jk) )
                  zfm_ui = pun(ji,jj,jk) - ABS( pun(ji,jj,jk) )
                  zfp_vj = pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) )
                  zfm_vj = pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) )
                  ! 2nd order centered advective fluxes (x2)
                  zcenut = pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj  ,jk,jn) )
                  zcenvt = pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji  ,jj+1,jk,jn) )
                  ! UBS advective fluxes
                  ztu(ji,jj,jk) = 0.5 * ( zcenut - zfp_ui * zltu(ji,jj,jk) - zfm_ui * zltu(ji+1,jj,jk) )
                  ztv(ji,jj,jk) = 0.5 * ( zcenvt - zfp_vj * zltv(ji,jj,jk) - zfm_vj * zltv(ji,jj+1,jk) )
               END DO
            END DO
         END DO                                           ! End of slab         

         zltu(:,:,:) = pta(:,:,:,jn)      ! store pta trends

         DO jk = 1, jpkm1                 ! Horizontal advective trends
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn)                        &
                     &             - (  ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk)    &
                     &                + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk)  ) / ( e1e2t(ji,jj) * fse3t(ji,jj,jk) )
               END DO
            END DO
            !                                             
         END DO                                           !   End of slab

         ! Horizontal trend used in tra_adv_ztvd subroutine
         zltu(:,:,:) = pta(:,:,:,jn) - zltu(:,:,:)

         !                
         IF( l_trd ) THEN                  ! trend diagnostics
             CALL trd_tra( kt, cdtype, jn, jptra_xad, ztu, pun, ptn(:,:,:,jn) )
             CALL trd_tra( kt, cdtype, jn, jptra_yad, ztv, pvn, ptn(:,:,:,jn) )
         END IF
         !                                 ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( cdtype == 'TRA' .AND. ln_diaptr ) CALL dia_ptr_ohst_components( jn, 'adv', ztv(:,:,:) )
         
         ! TVD scheme for the vertical direction  
         ! ----------------------
         IF( l_trd )   zltv(:,:,:) = pta(:,:,:,jn)          ! store pta if trend diag.

         !  Bottom value : flux set to zero
         ztw(:,:,jpk) = 0.e0   ;   zti(:,:,jpk) = 0.e0

         ! Surface value
         IF( lk_vvl ) THEN   ;   ztw(:,:,1) = 0.e0                      ! variable volume : flux set to zero
         ELSE                ;   ztw(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)   ! free constant surface 
         ENDIF
         !  upstream advection with initial mass fluxes & intermediate update
         ! -------------------------------------------------------------------
         ! Interior value
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                   zfp_wk = pwn(ji,jj,jk) + ABS( pwn(ji,jj,jk) )
                   zfm_wk = pwn(ji,jj,jk) - ABS( pwn(ji,jj,jk) )
                   ztw(ji,jj,jk) = 0.5 * (  zfp_wk * ptb(ji,jj,jk,jn) + zfm_wk * ptb(ji,jj,jk-1,jn)  )
               END DO
            END DO
         END DO 
         ! update and guess with monotonic sheme
         DO jk = 1, jpkm1
            z2dtt = p2dt(jk)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ztak = - ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) ) * zbtr
                  pta(ji,jj,jk,jn) =   pta(ji,jj,jk,jn) +  ztak 
                  zti(ji,jj,jk)    = ( ptb(ji,jj,jk,jn) + z2dtt * ( ztak + zltu(ji,jj,jk) ) ) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         CALL lbc_lnk( zti, 'T', 1. )      ! Lateral boundary conditions on zti, zsi   (unchanged sign)

         !  antidiffusive flux : high order minus low order
         ztw(:,:,1) = 0.e0       ! Surface value
         DO jk = 2, jpkm1        ! Interior value
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztw(ji,jj,jk) = 0.5 * pwn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj,jk-1,jn) ) - ztw(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         CALL nonosc_z( ptb(:,:,:,jn), ztw, zti, p2dt )      !  monotonicity algorithm

         !  final trend with corrected fluxes
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1 
               DO ji = fs_2, fs_jpim1   ! vector opt.   
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ! k- vertical advective trends  
                  ztra = - zbtr * ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) )
                  ! added to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO

         !  Save the final vertical advective trends
         IF( l_trd )  THEN                        ! vertical advective trend diagnostics
            DO jk = 1, jpkm1                       ! (compute -w.dk[ptn]= -dk[w.ptn] + ptn.dk[w])
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     zbtr = 1.e0 / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                     z_hdivn = (  pwn(ji,jj,jk) - pwn(ji,jj,jk+1)  ) * zbtr
                     zltv(ji,jj,jk) = pta(ji,jj,jk,jn) - zltv(ji,jj,jk) + ptn(ji,jj,jk,jn) * z_hdivn
                  END DO
               END DO
            END DO
            CALL trd_tra( kt, cdtype, jn, jptra_zad, zltv )
         ENDIF
         !
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztu, ztv, zltu, zltv, zti, ztw )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_ubs')
      !
   END SUBROUTINE tra_adv_ubs


   SUBROUTINE nonosc_z( pbef, pcc, paft, p2dt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc_z  ***
      !!     
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ), DIMENSION(jpk)          ::   p2dt   ! vertical profile of tracer time-step
      REAL(wp),                DIMENSION (jpi,jpj,jpk) ::   pbef   ! before field
      REAL(wp), INTENT(inout), DIMENSION (jpi,jpj,jpk) ::   paft   ! after field
      REAL(wp), INTENT(inout), DIMENSION (jpi,jpj,jpk) ::   pcc    ! monotonic flux in the k direction
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikm1         ! local integer
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn, z2dtt   ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zbetup, zbetdo
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('nonosc_z')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zbetup, zbetdo )
      !
      zbig  = 1.e+40_wp
      zrtrn = 1.e-15_wp
      zbetup(:,:,:) = 0._wp   ;   zbetdo(:,:,:) = 0._wp

      ! Search local extrema
      ! --------------------
      ! large negative value (-zbig) inside land
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
      ! search maximum in neighbourhood
      DO jk = 1, jpkm1
         ikm1 = MAX(jk-1,1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbetup(ji,jj,jk) = MAX(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
                  &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
                  &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
            END DO
         END DO
      END DO
      ! large positive value (+zbig) inside land
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
      ! search minimum in neighbourhood
      DO jk = 1, jpkm1
         ikm1 = MAX(jk-1,1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbetdo(ji,jj,jk) = MIN(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
                  &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
                  &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
            END DO
         END DO
      END DO

      ! restore masked values to zero
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:)
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:)


      ! 2. Positive and negative part of fluxes and beta terms
      ! ------------------------------------------------------

      DO jk = 1, jpkm1
         z2dtt = p2dt(jk)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! positive & negative part of the flux
               zpos = MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )
               zneg = MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )
               ! up & down beta terms
               zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) / z2dtt
               zbetup(ji,jj,jk) = ( zbetup(ji,jj,jk) - paft(ji,jj,jk) ) / (zpos+zrtrn) * zbt
               zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zbetdo(ji,jj,jk) ) / (zneg+zrtrn) * zbt
            END DO
         END DO
      END DO
      ! monotonic flux in the k direction, i.e. pcc
      ! -------------------------------------------
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj,jk-1) )
               zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj,jk-1) )
               zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pcc(ji,jj,jk) ) )
               pcc(ji,jj,jk) = pcc(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zbetup, zbetdo )
      !
      IF( nn_timing == 1 )  CALL timing_stop('nonosc_z')
      !
   END SUBROUTINE nonosc_z

   !!======================================================================
END MODULE traadv_ubs
