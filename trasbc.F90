MODULE trasbc
   !!==============================================================================
   !!                       ***  MODULE  trasbc  ***
   !! Ocean active tracers:  surface boundary condition
   !!==============================================================================
   !! History :  OPA  !  1998-10  (G. Madec, G. Roullet, M. Imbard)  Original code
   !!            8.2  !  2001-02  (D. Ludicone)  sea ice and free surface
   !!  NEMO      1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!             -   !  2010-09  (C. Ethe, G. Madec) Merge TRA-TRC
   !!            3.6  !  2014-11  (P. Mathiot) isf melting forcing 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_sbc      : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE sbc_oce         ! surface boundary condition: ocean
   USE dom_oce         ! ocean space domain variables
   USE phycst          ! physical constant
   USE sbcmod          ! ln_rnf  
   USE sbcrnf          ! River runoff  
   USE sbcisf          ! Ice shelf   
   USE traqsr          ! solar radiation penetration
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   !
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE iom
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE eosbn2
#if defined key_asminc   
   USE asminc          ! Assimilation increment
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_sbc    ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_sbc ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sbc  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (flux through the interface, concentration/dilution effect)
      !!      and add it to the general trend of tracer equations.
      !!
      !! ** Method :
      !!      Following Roullet and Madec (2000), the air-sea flux can be divided 
      !!      into three effects: (1) Fext, external forcing; 
      !!      (2) Fwi, concentration/dilution effect due to water exchanged 
      !!         at the surface by evaporation, precipitations and runoff (E-P-R); 
      !!      (3) Fwe, tracer carried with the water that is exchanged. 
      !!            - salinity    : salt flux only due to freezing/melting
      !!            sa = sa +  sfx / rau0 / e3t  for k=1
      !!
      !!      Fext, flux through the air-sea interface for temperature and salt: 
      !!            - temperature : heat flux q (w/m2). If penetrative solar
      !!         radiation q is only the non solar part of the heat flux, the
      !!         solar part is added in traqsr.F routine.
      !!            ta = ta + q /(rau0 rcp e3t)  for k=1
      !!            - salinity    : no salt flux
      !!
      !!      The formulation for Fwb and Fwi vary according to the free 
      !!      surface formulation (linear or variable volume). 
      !!      * Linear free surface
      !!            The surface freshwater flux modifies the ocean volume
      !!         and thus the concentration of a tracer and the temperature.
      !!         First order of the effect of surface freshwater exchange 
      !!         for salinity, it can be neglected on temperature (especially
      !!         as the temperature of precipitations and runoffs is usually
      !!         unknown).
      !!            - temperature : we assume that the temperature of both
      !!         precipitations and runoffs is equal to the SST, thus there
      !!         is no additional flux since in this case, the concentration
      !!         dilution effect is balanced by the net heat flux associated
      !!         to the freshwater exchange (Fwe+Fwi=0):
      !!            (Tp P - Te E) + SST (P-E) = 0 when Tp=Te=SST
      !!            - salinity    : evaporation, precipitation and runoff
      !!         water has a zero salinity  but there is a salt flux due to 
      !!         freezing/melting, thus:
      !!            sa = sa + emp * sn / rau0 / e3t   for k=1
      !!                    + sfx    / rau0 / e3t
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation minus runoff) given in kg/m2/s is divided
      !!         by rau0 (density of sea water) to obtain m/s.    
      !!         Note: even though Fwe does not appear explicitly for 
      !!         temperature in this routine, the heat carried by the water
      !!         exchanged through the surface is part of the total heat flux
      !!         forcing and must be taken into account in the global heat
      !!         balance).
      !!      * nonlinear free surface (variable volume, lk_vvl)
      !!         contrary to the linear free surface case, Fwi is properly 
      !!         taken into account by using the true layer thicknesses to       
      !!         calculate tracer content and advection. There is no need to 
      !!         deal with it in this routine.
      !!           - temperature: Fwe=SST (P-E+R) is added to Fext.
      !!           - salinity:  Fwe = 0, there is no surface flux of salt.
      !!
      !! ** Action  : - Update the 1st level of (ta,sa) with the trend associated
      !!                with the tracer surface boundary condition 
      !!              - send trends to trdtra module (l_trdtra=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices  
      INTEGER  ::   ikt, ikb 
      INTEGER  ::   nk_isf
      REAL(wp) ::   zfact, z1_e3t, zdep
      REAL(wp) ::   zalpha, zhk
      REAL(wp) ::  zt_frz, zpress
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_sbc')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_sbc : TRAcer Surface Boundary Condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      IF( l_trdtra ) THEN                    !* Save ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds ) 
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF

!!gm      IF( .NOT.ln_traqsr )   qsr(:,:) = 0.e0   ! no solar radiation penetration
      IF( .NOT.ln_traqsr ) THEN     ! no solar radiation penetration
         qns(:,:) = qns(:,:) + qsr(:,:)      ! total heat flux in qns
         qsr(:,:) = 0.e0                     ! qsr set to zero
      ENDIF

      !----------------------------------------
      !        EMP, SFX and QNS effects
      !----------------------------------------
      !                                          Set before sbc tracer content fields
      !                                          ************************************
      IF( kt == nit000 ) THEN                      ! Set the forcing field at nit000 - 1
         !                                         ! -----------------------------------
         IF( ln_rstart .AND.    &                     ! Restart: read in restart file
              & iom_varid( numror, 'sbc_hc_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nit000-1 surface tracer content forcing fields red in the restart file'
            zfact = 0.5_wp
            CALL iom_get( numror, jpdom_autoglo, 'sbc_hc_b', sbc_tsc_b(:,:,jp_tem) )   ! before heat content sbc trend
            CALL iom_get( numror, jpdom_autoglo, 'sbc_sc_b', sbc_tsc_b(:,:,jp_sal) )   ! before salt content sbc trend
         ELSE                                         ! No restart or restart not found: Euler forward time stepping
            zfact = 1._wp
            sbc_tsc(:,:,:) = 0._wp
            sbc_tsc_b(:,:,:) = 0._wp
         ENDIF
      ELSE                                         ! Swap of forcing fields
         !                                         ! ----------------------
         zfact = 0.5_wp
         sbc_tsc_b(:,:,:) = sbc_tsc(:,:,:)
      ENDIF
      !                                          Compute now sbc tracer content fields
      !                                          *************************************

                                                   ! Concentration dilution effect on (t,s) due to  
                                                   ! evaporation, precipitation and qns, but not river runoff 
                                               
      IF( lk_vvl ) THEN                            ! Variable Volume case  ==>> heat content of mass flux is in qns
         DO jj = 1, jpj
            DO ji = 1, jpi 
               sbc_tsc(ji,jj,jp_tem) = r1_rau0_rcp * qns(ji,jj)                              ! non solar heat flux
               sbc_tsc(ji,jj,jp_sal) = r1_rau0     * sfx(ji,jj)                              ! salt flux due to freezing/melting
            END DO
         END DO
      ELSE                                         ! Constant Volume case ==>> Concentration dilution effect
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! temperature : heat flux
               sbc_tsc(ji,jj,jp_tem) = r1_rau0_rcp * qns(ji,jj)                          &   ! non solar heat flux
                  &                  + r1_rau0     * emp(ji,jj)  * tsn(ji,jj,1,jp_tem)       ! concent./dilut. effect
               ! salinity    : salt flux + concent./dilut. effect (both in sfx)
               sbc_tsc(ji,jj,jp_sal) = r1_rau0  * (  sfx(ji,jj)                          &   ! salt flux (freezing/melting)
                  &                                + emp(ji,jj) * tsn(ji,jj,1,jp_sal) )      ! concent./dilut. effect
            END DO
         END DO
         IF( iom_use('emp_x_sst') )   CALL iom_put( "emp_x_sst", emp (:,:) * tsn(:,:,1,jp_tem) )   ! c/d term on sst
         IF( iom_use('emp_x_sss') )   CALL iom_put( "emp_x_sss", emp (:,:) * tsn(:,:,1,jp_sal) )   ! c/d term on sss
      ENDIF
      ! Concentration dilution effect on (t,s) due to evapouration, precipitation and qns, but not river runoff  
      DO jn = 1, jpts
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z1_e3t = zfact / fse3t(ji,jj,1)
               tsa(ji,jj,1,jn) = tsa(ji,jj,1,jn) + ( sbc_tsc_b(ji,jj,jn) + sbc_tsc(ji,jj,jn) ) * z1_e3t
            END DO
         END DO
      END DO
      !                                          Write in the ocean restart file
      !                                          *******************************
      IF( lrst_oce ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : ocean surface tracer content forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrow, 'sbc_hc_b', sbc_tsc(:,:,jp_tem) )
         CALL iom_rstput( kt, nitrst, numrow, 'sbc_sc_b', sbc_tsc(:,:,jp_sal) )
      ENDIF
      !
      !
      !----------------------------------------
      !       Ice Shelf effects (ISF)
      !     tbl treated as in Losh (2008) JGR
      !----------------------------------------
      !
      IF( nn_isf > 0 ) THEN
         zfact = 0.5e0
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1
         
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)
   
               ! level fully include in the ice shelf boundary layer
               ! if isfdiv, we have to remove heat flux due to inflow at 0oC (as in rnf when you add rnf at sst)
               ! sign - because fwf sign of evapo (rnf sign of precip)
               DO jk = ikt, ikb - 1
               ! compute tfreez for the temperature correction (we add water at freezing temperature)
               ! compute trend
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)                                          &
                     &           + zfact * (risf_tsc_b(ji,jj,jp_tem) + risf_tsc(ji,jj,jp_tem)) * r1_hisf_tbl(ji,jj)
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)                                          &
                     &           + zfact * (risf_tsc_b(ji,jj,jp_sal) + risf_tsc(ji,jj,jp_sal)) * r1_hisf_tbl(ji,jj)
               END DO
   
               ! level partially include in ice shelf boundary layer 
               ! compute tfreez for the temperature correction (we add water at freezing temperature)
               ! compute trend
               tsa(ji,jj,ikb,jp_tem) = tsa(ji,jj,ikb,jp_tem)                                           &
                  &              + zfact * (risf_tsc_b(ji,jj,jp_tem) + risf_tsc(ji,jj,jp_tem)) * r1_hisf_tbl(ji,jj) * ralpha(ji,jj)
               tsa(ji,jj,ikb,jp_sal) = tsa(ji,jj,ikb,jp_sal)                                           &
                  &              + zfact * (risf_tsc_b(ji,jj,jp_sal) + risf_tsc(ji,jj,jp_sal)) * r1_hisf_tbl(ji,jj) * ralpha(ji,jj) 
            END DO
         END DO
         IF( lrst_oce ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc : isf surface tracer content forcing fields written in ocean restart file ',   &
               &                    'at it= ', kt,' date= ', ndastp
            IF(lwp) WRITE(numout,*) '~~~~'
            CALL iom_rstput( kt, nitrst, numrow, 'fwf_isf_b', fwfisf(:,:)          )
            CALL iom_rstput( kt, nitrst, numrow, 'isf_hc_b' , risf_tsc(:,:,jp_tem) )
            CALL iom_rstput( kt, nitrst, numrow, 'isf_sc_b' , risf_tsc(:,:,jp_sal) )
         ENDIF
      END IF
      !
      !----------------------------------------
      !        River Runoff effects
      !----------------------------------------
      !
      IF( ln_rnf ) THEN         ! input of heat and salt due to river runoff 
         zfact = 0.5_wp
         DO jj = 2, jpj 
            DO ji = fs_2, fs_jpim1
               IF( rnf(ji,jj) /= 0._wp ) THEN
                  zdep = zfact / h_rnf(ji,jj)
                  DO jk = 1, nk_rnf(ji,jj)
                                        tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)   &
                                          &               +  ( rnf_tsc_b(ji,jj,jp_tem) + rnf_tsc(ji,jj,jp_tem) ) * zdep
                     IF( ln_rnf_sal )   tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)   &
                                          &               +  ( rnf_tsc_b(ji,jj,jp_sal) + rnf_tsc(ji,jj,jp_sal) ) * zdep 
                  END DO
               ENDIF
            END DO  
         END DO  
      ENDIF

      IF( iom_use('rnf_x_sst') )   CALL iom_put( "rnf_x_sst", rnf*tsn(:,:,1,jp_tem) )   ! runoff term on sst
      IF( iom_use('rnf_x_sss') )   CALL iom_put( "rnf_x_sss", rnf*tsn(:,:,1,jp_sal) )   ! runoff term on sss

#if defined key_asminc
! WARNING: THIS MAY WELL NOT BE REQUIRED - WE DON'T WANT TO CHANGE T&S BUT THIS MAY COMPENSATE ANOTHER TERM...
! Rate of change in e3t for each level is ssh_iau*e3t_0/ht_0
! Contribution to tsa should be rate of change in level / per m of ocean? (hence the division by fse3t_n)
      IF( ln_sshinc ) THEN         ! input of heat and salt due to assimilation
         DO jj = 2, jpj 
            DO ji = fs_2, fs_jpim1
               zdep = ssh_iau(ji,jj) / ( ht_0(ji,jj) + 1.0 - ssmask(ji, jj) )
               DO jk = 1, jpkm1
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)   &
                                        &            + tsn(ji,jj,jk,jp_tem) * zdep * ( e3t_0(ji,jj,jk) / fse3t_n(ji,jj,jk) )
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)   &
                                        &            + tsn(ji,jj,jk,jp_sal) * zdep * ( e3t_0(ji,jj,jk) / fse3t_n(ji,jj,jk) )
               END DO
            END DO  
         END DO  
      ENDIF
#endif
 
      IF( l_trdtra )   THEN                      ! send trends for further diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal) - ztrds(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_nsr, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_nsr, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds ) 
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' sbc  - Ta: ', mask1=tmask,   &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_sbc')
      !
   END SUBROUTINE tra_sbc

   !!======================================================================
END MODULE trasbc
