MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_meso       :   Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init  :   Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zint          !  interpolation and computation of various fields
   USE p4zprod         !  production
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xprefc       !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xprefp       !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xprefz       !: mesozoo preference for diatoms
   REAL(wp), PUBLIC ::  xprefpoc     !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: non assimilated fraction of P by mesozoo 
   REAL(wp), PUBLIC ::  unass2       !: Efficicency of mesozoo growth 
   REAL(wp), PUBLIC ::  sigma2       !: Fraction of mesozoo excretion as DOM 
   REAL(wp), PUBLIC ::  epsher2      !: half sturation constant for grazing 2
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam
      REAL(wp) :: zgraze2 , zdenom, zdenom2
      REAL(wp) :: zfact   , zstep, zfood, zfoodlim, zproport
      REAL(wp) :: zmortzgoc, zfrac, zfracfe, zratio, zratio2
      REAL(wp) :: zepshert, zepsherv, zgrarsig, zgraztot, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem2, zgrafer2, zgrapoc2, zprcaca, zmortz2, zgrasrat, zgrasratn
#if defined key_kriest
      REAL znumpoc
#endif
      REAL(wp) :: zrespz2, ztortz2, zgrazd, zgrazz, zgrazpof
      REAL(wp) :: zgrazn, zgrazpoc, zgraznf, zgrazf
      REAL(wp) :: zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zgrazing, zw3d

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_meso')
      !
      IF( lk_iomput ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, zgrazing )
         zgrazing(:,:,:) = 0._wp
      ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompam   = MAX( ( trb(ji,jj,jk,jpmes) - 1.e-9 ), 0.e0 )
# if defined key_degrad
               zstep     = xstep * facvol(ji,jj,jk)
# else
               zstep     = xstep
# endif
               zfact     = zstep * tgfunc2(ji,jj,jk) * zcompam

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz2   = resrat2 * zfact * trb(ji,jj,jk,jpmes) / ( xkmort + trb(ji,jj,jk,jpmes) )  &
                  &      + resrat2 * zfact * 3. * nitrfac(ji,jj,jk)

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation
               !  ---------------------------------------------------------------
               ztortz2   = mzrat2 * 1.e6 * zfact * trb(ji,jj,jk,jpmes)
               !
               zcompadi  = MAX( ( trb(ji,jj,jk,jpdia) - xthresh2dia ), 0.e0 )
               zcompaz   = MAX( ( trb(ji,jj,jk,jpzoo) - xthresh2zoo ), 0.e0 )
               ! Size effect of nanophytoplankton on grazing : the smaller it is, the less prone
               ! it is to predation by mesozooplankton
               ! -------------------------------------------------------------------------------
               zcompaph  = MAX( ( trb(ji,jj,jk,jpphy) - xthresh2phy ), 0.e0 ) &
                  &      * MIN(1., MAX( 0., ( quotan(ji,jj,jk) - 0.2) / 0.3 ) )
               zcompapoc = MAX( ( trb(ji,jj,jk,jppoc) - xthresh2poc ), 0.e0 )

               zfood     = xprefc * zcompadi + xprefz * zcompaz + xprefp * zcompaph + xprefpoc * zcompapoc 
               zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
               zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )
               zgraze2   = grazrat2 * zstep * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpmes) 

               zgrazd    = zgraze2  * xprefc   * zcompadi  * zdenom2 
               zgrazz    = zgraze2  * xprefz   * zcompaz   * zdenom2 
               zgrazn    = zgraze2  * xprefp   * zcompaph  * zdenom2 
               zgrazpoc  = zgraze2  * xprefpoc * zcompapoc * zdenom2 

               zgraznf   = zgrazn   * trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) + rtrn)
               zgrazf    = zgrazd   * trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) + rtrn)
               zgrazpof  = zgrazpoc * trb(ji,jj,jk,jpsfe) / ( trb(ji,jj,jk,jppoc) + rtrn)

               !  Mesozooplankton flux feeding on GOC
               !  ----------------------------------
               !  ----------------------------------
# if ! defined key_kriest
               zgrazffeg = grazflux  * zstep * wsbio4(ji,jj,jk)      &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpgoc) * trb(ji,jj,jk,jpmes)
               zgrazfffg = zgrazffeg * trb(ji,jj,jk,jpbfe) / (trb(ji,jj,jk,jpgoc) + rtrn)
# endif
               zgrazffep = grazflux  * zstep *  wsbio3(ji,jj,jk)     &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpmes)
               zgrazfffp = zgrazffep * trb(ji,jj,jk,jpsfe) / (trb(ji,jj,jk,jppoc) + rtrn)
              !
# if ! defined key_kriest
              zgraztot  = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg
              ! Compute the proportion of filter feeders
              zproport  = (zgrazffep + zgrazffeg)/(rtrn + zgraztot)
              ! Compute fractionation of aggregates. It is assumed that 
              ! diatoms based aggregates are more prone to fractionation
              ! since they are more porous (marine snow instead of fecal pellets)
              zratio    = trb(ji,jj,jk,jpgsi) / ( trb(ji,jj,jk,jpgoc) + rtrn )
              zratio2   = zratio * zratio
              zfrac     = zproport * grazflux  * zstep * wsbio4(ji,jj,jk)      &
               &          * trb(ji,jj,jk,jpgoc) * trb(ji,jj,jk,jpmes)          &
               &          * ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 ) )
              zfracfe   = zfrac * trb(ji,jj,jk,jpbfe) / (trb(ji,jj,jk,jpgoc) + rtrn)

              zgrazffep = zproport * zgrazffep
              zgrazffeg = zproport * zgrazffeg
              zgrazfffp = zproport * zgrazfffp
              zgrazfffg = zproport * zgrazfffg
              zgraztot  = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg
              zgraztotn = zgrazd * quotad(ji,jj,jk) + zgrazz + zgrazn * quotan(ji,jj,jk)   &
              &   + zgrazpoc + zgrazffep + zgrazffeg
              zgraztotf = zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfffp + zgrazfffg
# else
              zgraztot  = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep
              ! Compute the proportion of filter feeders
              zproport  = zgrazffep / ( zgraztot + rtrn )
              zgrazffep = zproport * zgrazffep
              zgrazfffp = zproport * zgrazfffp
              zgraztot  = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep
              zgraztotn = zgrazd * quotad(ji,jj,jk) + zgrazz + zgrazn * quotan(ji,jj,jk) + zgrazpoc + zgrazffep
              zgraztotf = zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfffp
# endif

              ! Total grazing ( grazing by microzoo is already computed in p4zmicro )
              IF( lk_iomput )  zgrazing(ji,jj,jk) = zgraztot

              !    Mesozooplankton efficiency
              !    --------------------------
               zgrasrat  =  ( zgraztotf +rtrn )/ ( zgraztot + rtrn )
               zgrasratn =  ( zgraztotn +rtrn )/ ( zgraztot + rtrn )
               zepshert  = MIN( 1., zgrasratn, zgrasrat / ferat3)
               zepsherv  = zepshert * MIN( epsher2, (1. - unass2) * zgrasrat / ferat3, (1. - unass2) * zgrasratn )
               zgrarem2  = zgraztot * ( 1. - zepsherv - unass2 ) &
                &       + ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 ) * ztortz2
               zgrafer2  = zgraztot * MAX( 0. , ( 1. - unass2 ) * zgrasrat - ferat3 * zepsherv )    &
                &       + ferat3 * ( ( 1. - epsher2 - unass2 ) /( 1. - epsher2 ) * ztortz2 )
               zgrapoc2  = zgraztot * unass2

               !   Update the arrays TRA which contain the biological sources and sinks
               zgrarsig  = zgrarem2 * sigma2
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem2 - zgrarsig
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig              

               zmortz2 = ztortz2 + zrespz2
               zmortzgoc = unass2 / ( 1. - epsher2 ) * ztortz2 + zrespz2
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zmortz2 + zepsherv * zgraztot 
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazd
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazn
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazn * trb(ji,jj,jk,jpnch) / ( trb(ji,jj,jk,jpphy) + rtrn )
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazd * trb(ji,jj,jk,jpdch) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazd * trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazd * trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazf

               ! calcite production
               zprcaca = xfracal(ji,jj,jk) * zgrazn
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part2 * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
#if defined key_kriest
              znumpoc = trb(ji,jj,jk,jpnum) / ( trb(ji,jj,jk,jppoc) + rtrn )
              tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortzgoc - zgrazpoc - zgrazffep + zgrapoc2
              tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) - zgrazpoc * znumpoc + zgrapoc2 * xkr_dmeso      &
                 &   + zmortzgoc * xkr_dmeso - zgrazffep * znumpoc * wsbio4(ji,jj,jk) / ( wsbio3(ji,jj,jk) + rtrn )
              tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortzgoc - zgrazfffp - zgrazpof    &
                 &                 + zgraztotf * unass2
#else
              tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zgrazpoc - zgrazffep + zfrac
              tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zmortzgoc - zgrazffeg + zgrapoc2 - zfrac
              tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zgrazpof - zgrazfffp + zfracfe
              tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ferat3 * zmortzgoc - zgrazfffg     &
                 &                + zgraztotf * unass2 - zfracfe
#endif
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, zw3d )
         IF( iom_use( "GRAZ2" ) ) THEN
            zw3d(:,:,:) = zgrazing(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !   Total grazing of phyto by zooplankton
            CALL iom_put( "GRAZ2", zw3d )
         ENDIF
         IF( iom_use( "PCAL" ) ) THEN
            zw3d(:,:,:) = prodcal(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)   !  Calcite production
            CALL iom_put( "PCAL", zw3d )  
         ENDIF
         CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( lk_iomput )  CALL wrk_dealloc( jpi, jpj, jpk, zgrazing )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso

   SUBROUTINE p4z_meso_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismes/ part2, grazrat2, resrat2, mzrat2, xprefc, xprefp, xprefz,   &
         &                xprefpoc, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2, xkgraz2, epsher2, sigma2, unass2, grazflux
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampismes in reference namelist : Pisces mesozooplankton
      READ  ( numnatp_ref, nampismes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismes in configuration namelist : Pisces mesozooplankton
      READ  ( numnatp_cfg, nampismes, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismes )


      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, nampismes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '    mesozoo preference for phyto                   xprefc       =', xprefc
         WRITE(numout,*) '    mesozoo preference for POC                     xprefp       =', xprefp
         WRITE(numout,*) '    mesozoo preference for zoo                     xprefz       =', xprefz
         WRITE(numout,*) '    mesozoo preference for poc                     xprefpoc     =', xprefpoc
         WRITE(numout,*) '    microzoo feeding threshold  for mesozoo        xthresh2zoo  =', xthresh2zoo
         WRITE(numout,*) '    diatoms feeding threshold  for mesozoo         xthresh2dia  =', xthresh2dia
         WRITE(numout,*) '    nanophyto feeding threshold for mesozoo        xthresh2phy  =', xthresh2phy
         WRITE(numout,*) '    poc feeding threshold for mesozoo              xthresh2poc  =', xthresh2poc
         WRITE(numout,*) '    feeding threshold for mesozooplankton          xthresh2     =', xthresh2
         WRITE(numout,*) '    exsudation rate of mesozooplankton             resrat2      =', resrat2
         WRITE(numout,*) '    mesozooplankton mortality rate                 mzrat2       =', mzrat2
         WRITE(numout,*) '    maximal mesozoo grazing rate                   grazrat2     =', grazrat2
         WRITE(numout,*) '    mesozoo flux feeding rate                      grazflux     =', grazflux
         WRITE(numout,*) '    non assimilated fraction of P by mesozoo       unass2       =', unass2
         WRITE(numout,*) '    Efficicency of Mesozoo growth                  epsher2      =', epsher2
         WRITE(numout,*) '    Fraction of mesozoo excretion as DOM           sigma2       =', sigma2
         WRITE(numout,*) '    half sturation constant for grazing 2          xkgraz2      =', xkgraz2
      ENDIF


   END SUBROUTINE p4z_meso_init


#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_meso                    ! Empty routine
   END SUBROUTINE p4z_meso
#endif 

   !!======================================================================
END MODULE p4zmeso
