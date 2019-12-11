MODULE p4zfechem
   !!======================================================================
   !!                         ***  MODULE p4zfechem  ***
   !! TOP :   PISCES Compute iron chemistry and scavenging
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, A. Tagliabue, C. Ethe) Original code
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_fechem       :  Compute remineralization/scavenging of iron
   !!   p4z_fechem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_fechem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
   USE p4zche          !  chemical model
   USE p4zsbc          !  Boundary conditions from sediments
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_fechem      ! called in p4zbio.F90
   PUBLIC   p4z_fechem_init ! called in trcsms_pisces.F90

   !! * Shared module variables
   LOGICAL          ::  ln_fechem    !: boolean for complex iron chemistry following Tagliabue and voelker
   LOGICAL          ::  ln_ligvar    !: boolean for variable ligand concentration following Tagliabue and voelker
   REAL(wp), PUBLIC ::  xlam1        !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::  xlamdust     !: scavenging rate of Iron by dust 
   REAL(wp), PUBLIC ::  ligand       !: ligand concentration in the ocean 

   REAL(wp) :: kl1, kl2, kb1, kb2, ks, kpr, spd, con, kth

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_fechem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_fechem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of iron
      !!
      !! ** Method  :   2 different chemistry models are available for iron
      !!                (1) The simple chemistry model of Aumont and Bopp (2006)
      !!                    based on one ligand and one inorganic form
      !!                (2) The complex chemistry model of Tagliabue and 
      !!                    Voelker (2009) based on 2 ligands, 2 inorganic forms
      !!                    and one particulate form (ln_fechem)
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk, jic
      REAL(wp) ::   zdep, zlam1a, zlam1b, zlamfac
      REAL(wp) ::   zkeq, zfeequi, zfesatur, zfecoll
      REAL(wp) ::   zdenom1, zscave, zaggdfea, zaggdfeb, zcoag
      REAL(wp) ::   ztrc, zdust
#if ! defined key_kriest
      REAL(wp) ::   zdenom, zdenom2
#endif
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zTL1, zFe3, ztotlig
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zFeL1, zFeL2, zTL2, zFe2, zFeP
      REAL(wp) :: zkox, zkph1, zkph2, zph, zionic, ztligand
      REAL(wp) :: za, zb, zc, zkappa1, zkappa2, za0, za1, za2
      REAL(wp) :: zxs, zfunc, zp, zq, zd, zr, zphi, zfff, zp3, zq2
      REAL(wp) :: ztfe, zoxy
      REAL(wp) :: zstep
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_fechem')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, zFe3, zFeL1, zTL1, ztotlig )
      zFe3 (:,:,:) = 0.
      zFeL1(:,:,:) = 0.
      zTL1 (:,:,:) = 0.
      IF( ln_fechem ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, zFe2, zFeL2, zTL2, zFeP )
         zFe2 (:,:,:) = 0.
         zFeL2(:,:,:) = 0.
         zTL2 (:,:,:) = 0.
         zFeP (:,:,:) = 0.
      ENDIF

      ! Total ligand concentration : Ligands can be chosen to be constant or variable
      ! Parameterization from Tagliabue and Voelker (2011)
      ! -------------------------------------------------
      IF( ln_ligvar ) THEN
         ztotlig(:,:,:) =  0.09 * trb(:,:,:,jpdoc) * 1E6 + ligand * 1E9
         ztotlig(:,:,:) =  MIN( ztotlig(:,:,:), 10. )
      ELSE
         ztotlig(:,:,:) = ligand * 1E9
      ENDIF

      IF( ln_fechem ) THEN
         ! ------------------------------------------------------------
         ! NEW FE CHEMISTRY ROUTINE from Tagliabue and Volker (2009)
         ! This model is based on two ligands, Fe2+, Fe3+ and Fep
         ! Chemistry is supposed to be fast enough to be at equilibrium
         ! ------------------------------------------------------------
!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  ! Calculate ligand concentrations : assume 2/3rd of excess goes to
                  ! strong ligands (L1) and 1/3rd to weak ligands (L2)
                  ztligand       = ztotlig(ji,jj,jk) - ligand * 1E9
                  zTL1(ji,jj,jk) =                0.000001 + 0.67 * ztligand
                  zTL2(ji,jj,jk) = ligand * 1E9 - 0.000001 + 0.33 * ztligand
                  ! ionic strength from Millero et al. 1987
                  zionic = 19.9201 * tsn(ji,jj,jk,jp_sal) / ( 1000. - 1.00488 * tsn(ji,jj,jk,jp_sal) + rtrn )
                  zph    = -LOG10( MAX( hi(ji,jj,jk), rtrn) )
                  zoxy   = trb(ji,jj,jk,jpoxy) * ( rhop(ji,jj,jk) / 1.e3 )
                  ! Fe2+ oxydation rate from Santana-Casiano et al. (2005)
                  zkox   = 35.407 - 6.7109 * zph + 0.5342 * zph * zph - 5362.6 / ( tsn(ji,jj,jk,jp_tem) + 273.15 )  &
                    &    - 0.04406 * SQRT( tsn(ji,jj,jk,jp_sal) ) - 0.002847 * tsn(ji,jj,jk,jp_sal)
                  zkox   = ( 10.** zkox ) * spd
                  zkox   = zkox * MAX( 1.e-6, zoxy) / ( chemo2(ji,jj,jk) + rtrn )
                  ! PHOTOREDUCTION of complexed iron : Tagliabue and Arrigo (2006)
                  zkph2 = MAX( 0., 15. * etot(ji,jj,jk) / ( etot(ji,jj,jk) + 2. ) )
                  zkph1 = zkph2 / 5.
                  ! pass the dfe concentration from PISCES
                  ztfe = trb(ji,jj,jk,jpfer) * 1e9
                  ! ----------------------------------------------------------
                  ! ANALYTICAL SOLUTION OF ROOTS OF THE FE3+ EQUATION
                  ! As shown in Tagliabue and Voelker (2009), Fe3+ is the root of a 3rd order polynom. 
                  ! ----------------------------------------------------------
                  ! calculate some parameters
                  za = 1 + ks / kpr
                  zb = 1 + ( zkph1 + kth ) / ( zkox + rtrn )
                  zc = 1 + zkph2 / ( zkox + rtrn )
                  zkappa1 = ( kb1 + zkph1 + kth ) / kl1
                  zkappa2 = ( kb2 + zkph2 ) / kl2
                  za2 = zTL1(ji,jj,jk) * zb / za + zTL2(ji,jj,jk) * zc / za + zkappa1 + zkappa2 - ztfe / za
                  za1 = zkappa2 * zTL1(ji,jj,jk) * zb / za + zkappa1 * zTL2(ji,jj,jk) * zc / za &
                      & + zkappa1 * zkappa2 - ( zkappa1 + zkappa2 ) * ztfe / za
                  za0 = -zkappa1 * zkappa2 * ztfe / za
                  zp  = za1 - za2 * za2 / 3.
                  zq  = za2 * za2 * za2 * 2. / 27. - za2 * za1 / 3. + za0
                  zp3 = zp / 3.
                  zq2 = zq / 2.
                  zd  = zp3 * zp3 * zp3 + zq2 * zq2
                  zr  = zq / ABS( zq ) * SQRT( ABS( zp ) / 3. )
                  ! compute the roots
                  IF( zp > 0.) THEN
                     ! zphi = ASINH( zq / ( 2. * zr * zr * zr ) )
                     zphi =  zq / ( 2. * zr * zr * zr ) 
                     zphi = LOG( zphi + SQRT( zphi * zphi + 1 ) )  ! asinh(x) = log(x + sqrt(x^2+1))
                     zxs  = -2. * zr * SINH( zphi / 3. ) - za1 / 3.
                  ELSE
                     IF( zd > 0. ) THEN
                        zfff = MAX( 1., zq / ( 2. * zr * zr * zr ) )
                        ! zphi = ACOSH( zfff )
                        zphi = LOG( zfff + SQRT( zfff * zfff - 1 ) )  ! acosh(x) = log(x + sqrt(x^2-1))
                        zxs = -2. * zr * COSH( zphi / 3. ) - za1 / 3.
                     ELSE
                        zfff = MIN( 1., zq / ( 2. * zr * zr * zr ) )
                        zphi = ACOS( zfff )
                        DO jic = 1, 3
                           zfunc = -2 * zr * COS( zphi / 3. + 2. * FLOAT( jic - 1 ) * rpi / 3. ) - za2 / 3.
                           IF( zfunc > 0. .AND. zfunc <= ztfe)  zxs = zfunc
                        END DO
                     ENDIF
                  ENDIF
                  ! solve for the other Fe species
                  zFe3(ji,jj,jk) = MAX( 0., zxs ) 
                  zFep(ji,jj,jk) = MAX( 0., ( ks * zFe3(ji,jj,jk) / kpr ) )
                  zkappa2 = ( kb2 + zkph2 ) / kl2
                  zFeL2(ji,jj,jk) = MAX( 0., ( zFe3(ji,jj,jk) * zTL2(ji,jj,jk) ) / ( zkappa2 + zFe3(ji,jj,jk) ) )
                  zFeL1(ji,jj,jk) = MAX( 0., ( ztfe / zb - za / zb * zFe3(ji,jj,jk) - zc / zb * zFeL2(ji,jj,jk) ) )
                  zFe2 (ji,jj,jk) = MAX( 0., ( ( zkph1 * zFeL1(ji,jj,jk) + zkph2 * zFeL2(ji,jj,jk) ) / zkox ) )
               END DO
            END DO
         END DO
      ELSE
         ! ------------------------------------------------------------
         ! OLD FE CHEMISTRY ROUTINE from Aumont and Bopp (2006)
         ! This model is based on one ligand and Fe' 
         ! Chemistry is supposed to be fast enough to be at equilibrium
         ! ------------------------------------------------------------
!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  zTL1(ji,jj,jk) = ztotlig(ji,jj,jk)
                  zkeq           = fekeq(ji,jj,jk)
                  zfesatur       = zTL1(ji,jj,jk) * 1E-9
                  ztfe           = trb(ji,jj,jk,jpfer) 
                  ! Fe' is the root of a 2nd order polynom
                  zFe3 (ji,jj,jk) = ( -( 1. + zfesatur * zkeq - zkeq * ztfe )               &
                     &             + SQRT( ( 1. + zfesatur * zkeq - zkeq * ztfe )**2       &
                     &               + 4. * ztfe * zkeq) ) / ( 2. * zkeq )
                  zFe3 (ji,jj,jk) = zFe3(ji,jj,jk) * 1E9
                  zFeL1(ji,jj,jk) = MAX( 0., trb(ji,jj,jk,jpfer) * 1E9 - zFe3(ji,jj,jk) )
              END DO
            END DO
         END DO
         !
      ENDIF

      zdust = 0.         ! if no dust available
!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zstep = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               ! Scavenging rate of iron. This scavenging rate depends on the load of particles of sea water. 
               ! This parameterization assumes a simple second order kinetics (k[Particles][Fe]).
               ! Scavenging onto dust is also included as evidenced from the DUNE experiments.
               ! --------------------------------------------------------------------------------------
               IF( ln_fechem ) THEN
                  zfeequi = ( zFe3(ji,jj,jk) + zFe2(ji,jj,jk) + zFeP(ji,jj,jk) ) * 1E-9
                  zfecoll = ( 0.3 * zFeL1(ji,jj,jk) + 0.5 * zFeL2(ji,jj,jk) ) * 1E-9
               ELSE
                  zfeequi = zFe3(ji,jj,jk) * 1E-9 
                  zfecoll = 0.5 * zFeL1(ji,jj,jk) * 1E-9
               ENDIF
#if defined key_kriest
               ztrc   = ( trb(ji,jj,jk,jppoc) + trb(ji,jj,jk,jpcal) + trb(ji,jj,jk,jpgsi) ) * 1.e6 
#else
               ztrc   = ( trb(ji,jj,jk,jppoc) + trb(ji,jj,jk,jpgoc) + trb(ji,jj,jk,jpcal) + trb(ji,jj,jk,jpgsi) ) * 1.e6 
#endif
               IF( ln_dust )  zdust  = dust(ji,jj) / ( wdust / rday ) * tmask(ji,jj,jk) ! dust in kg/m2/s
               zlam1b = 3.e-5 + xlamdust * zdust + xlam1 * ztrc
               zscave = zfeequi * zlam1b * zstep

               ! Compute the different ratios for scavenging of iron
               ! to later allocate scavenged iron to the different organic pools
               ! ---------------------------------------------------------
               zdenom1 = xlam1 * trb(ji,jj,jk,jppoc) / zlam1b
#if ! defined key_kriest
               zdenom2 = xlam1 * trb(ji,jj,jk,jpgoc) / zlam1b
#endif

               !  Increased scavenging for very high iron concentrations found near the coasts 
               !  due to increased lithogenic particles and let say it is unknown processes (precipitation, ...)
               !  -----------------------------------------------------------
               zlamfac = MAX( 0.e0, ( gphit(ji,jj) + 55.) / 30. )
               zlamfac = MIN( 1.  , zlamfac )
               zdep    = MIN( 1., 1000. / fsdept(ji,jj,jk) )
               zlam1b  = xlam1 * MAX( 0.e0, ( trb(ji,jj,jk,jpfer) * 1.e9 - ztotlig(ji,jj,jk) ) )
               zcoag   = zfeequi * zlam1b * zstep + 1E-4 * ( 1. - zlamfac ) * zdep * zstep * trb(ji,jj,jk,jpfer)

               !  Compute the coagulation of colloidal iron. This parameterization 
               !  could be thought as an equivalent of colloidal pumping.
               !  It requires certainly some more work as it is very poorly constrained.
               !  ----------------------------------------------------------------
               zlam1a  = ( 0.369  * 0.3 * trb(ji,jj,jk,jpdoc) + 102.4  * trb(ji,jj,jk,jppoc) ) * xdiss(ji,jj,jk)    &
                   &   + ( 114.   * 0.3 * trb(ji,jj,jk,jpdoc) + 5.09E3 * trb(ji,jj,jk,jppoc) )
               zaggdfea = zlam1a * zstep * zfecoll
#if defined key_kriest
               zaggdfeb = 0.
               !
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfea - zaggdfeb - zcoag
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1 + zaggdfea + zaggdfeb
#else
               zlam1b = 3.53E3 *   trb(ji,jj,jk,jpgoc) * xdiss(ji,jj,jk)
               zaggdfeb = zlam1b * zstep * zfecoll
               !
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfea - zaggdfeb - zcoag
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1 + zaggdfea
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zscave * zdenom2 + zaggdfeb
#endif
            END DO
         END DO
      END DO
      !
      !  Define the bioavailable fraction of iron
      !  ----------------------------------------
      IF( ln_fechem ) THEN
          biron(:,:,:) = MAX( 0., trb(:,:,:,jpfer) - zFeP(:,:,:) * 1E-9 )
      ELSE
          biron(:,:,:) = trb(:,:,:,jpfer) 
      ENDIF

      !  Output of some diagnostics variables
      !     ---------------------------------
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         IF( iom_use("Fe3")    )  CALL iom_put("Fe3"    , zFe3   (:,:,:)       * tmask(:,:,:) )   ! Fe3+
         IF( iom_use("FeL1")   )  CALL iom_put("FeL1"   , zFeL1  (:,:,:)       * tmask(:,:,:) )   ! FeL1
         IF( iom_use("TL1")    )  CALL iom_put("TL1"    , zTL1   (:,:,:)       * tmask(:,:,:) )   ! TL1
         IF( iom_use("Totlig") )  CALL iom_put("Totlig" , ztotlig(:,:,:)       * tmask(:,:,:) )   ! TL
         IF( iom_use("Biron")  )  CALL iom_put("Biron"  , biron  (:,:,:) * 1e9 * tmask(:,:,:) )   ! biron
         IF( ln_fechem ) THEN
            IF( iom_use("Fe2")  ) CALL iom_put("Fe2"    , zFe2   (:,:,:)       * tmask(:,:,:) )   ! Fe2+
            IF( iom_use("FeL2") ) CALL iom_put("FeL2"   , zFeL2  (:,:,:)       * tmask(:,:,:) )   ! FeL2
            IF( iom_use("FeP")  ) CALL iom_put("FeP"    , zFeP   (:,:,:)       * tmask(:,:,:) )   ! FeP
            IF( iom_use("TL2")  ) CALL iom_put("TL2"    , zTL2   (:,:,:)       * tmask(:,:,:) )   ! TL2
         ENDIF
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('fechem')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
                       CALL wrk_dealloc( jpi, jpj, jpk, zFe3, zFeL1, zTL1, ztotlig )
      IF( ln_fechem )  CALL wrk_dealloc( jpi, jpj, jpk, zFe2, zFeL2, zTL2, zFeP )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_fechem')
      !
   END SUBROUTINE p4z_fechem


   SUBROUTINE p4z_fechem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_fechem_init  ***
      !!
      !! ** Purpose :   Initialization of iron chemistry parameters
      !!
      !! ** Method  :   Read the nampisfer namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisfer
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisfer/ ln_fechem, ln_ligvar, xlam1, xlamdust, ligand 
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampisfer in reference namelist : Pisces iron chemistry
      READ  ( numnatp_ref, nampisfer, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisfer in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisfer in configuration namelist : Pisces iron chemistry
      READ  ( numnatp_cfg, nampisfer, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisfer in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisfer )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for Iron chemistry, nampisfer'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    enable complex iron chemistry scheme      ln_fechem =', ln_fechem
         WRITE(numout,*) '    variable concentration of ligand          ln_ligvar =', ln_ligvar
         WRITE(numout,*) '    scavenging rate of Iron                   xlam1     =', xlam1
         WRITE(numout,*) '    scavenging rate of Iron by dust           xlamdust  =', xlamdust
         WRITE(numout,*) '    ligand concentration in the ocean         ligand    =', ligand
      ENDIF
      !
      IF( ln_fechem ) THEN
         ! initialization of some constants used by the complexe chemistry scheme
         ! ----------------------------------------------------------------------
         spd = 3600. * 24.
         con = 1.E9
         ! LIGAND KINETICS (values from Witter et al. 2000)
         ! Weak (L2) ligands
         ! Phaeophytin
         kl2 = 12.2E5  * spd / con
         kb2 = 12.3E-6 * spd
         ! Strong (L1) ligands
         ! Saccharides
         ! kl1 = 12.2E5  * spd / con
         ! kb1 = 12.3E-6 * spd
         ! DFOB-
         kl1 = 19.6e5  * spd / con
         kb1 = 1.5e-6  * spd
         ! pcp and remin of Fe3p
         ks  = 0.075
         kpr = 0.05
         ! thermal reduction of Fe3
         kth = 0.0048 * 24.
         !
      ENDIF
      !
   END SUBROUTINE p4z_fechem_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_fechem                    ! Empty routine
   END SUBROUTINE p4z_fechem
#endif 

   !!======================================================================
END MODULE p4zfechem
