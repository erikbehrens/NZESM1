MODULE domc1d
   !!======================================================================
   !!                     ***  MODULE  domc1d  ***
   !! Ocean Domain : 1D column position from lat/lon namelist specification
   !!======================================================================
   !! History :  3.5  !  2013-04  (D. Calvert)  Original code
   !!----------------------------------------------------------------------
#if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'   :                                      1D Configuration
   !!----------------------------------------------------------------------
   !!   dom_c1d     : Determine jpizoom/jpjzoom from a given lat/lon
   !!----------------------------------------------------------------------
   USE phycst                        ! Physical constants (and par_oce)
   USE iom                           ! I/O library (iom_get)
   USE in_out_manager                ! I/O manager (ctmp1)
   USE dom_oce , ONLY : nimpp, njmpp ! Shared/distributed memory setting (mpp_init routine)
   USE wrk_nemo                      ! Memory allocation
   USE timing                        ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_c1d                  ! Routine called in domcfg.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE dom_c1d( plat, plon )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dom_c1d  ***
      !! 
      !! ** Purpose : Recalculate jpizoom/jpjzoom indices from lat/lon point
      !!
      !! ** Method  : Calculate global gphit and glamt as for dom_hgr.
      !!              After, find closest grid point to lat/lon point as for 
      !!              dom_ngb on T grid. From this infer jpizoom and jpjzoom.
      !!
      !! ** Action  : Recalculate jpizoom, jpjzoom (indices of C1D zoom)
      !!----------------------------------------------------------------------
      NAMELIST/namdom/ nn_bathy, rn_bathy , rn_e3zps_min, rn_e3zps_rat, nn_msh, rn_hmin,   &
         &             nn_acc   , rn_atfp     , rn_rdt      , rn_rdtmin ,                  &
         &             rn_rdtmax, rn_rdth     , nn_closea , ln_crs,    &
         &             jphgr_msh, &
         &             ppglam0, ppgphi0, ppe1_deg, ppe2_deg, ppe1_m, ppe2_m, &
         &             ppsur, ppa0, ppa1, ppkth, ppacr, ppdzmin, pphmax, ldbletanh, &
         &             ppa2, ppkth2, ppacr2

      INTEGER  ::  ji, jj                          ! Dummy loop indices
      INTEGER  ::  inum                            ! Coordinate file handle (case 0)
      INTEGER  ::  ijeq                            ! Index of equator T point (case 4)
      INTEGER  ::  ios                             ! Local integer output status for namelist read

      INTEGER , DIMENSION(2) ::   iloc             ! Minloc returned indices

      REAL(wp), INTENT(in) ::  plat                ! Column latitude
      REAL(wp), INTENT(in) ::  plon                ! Column longitude

      REAL(wp) ::  zlon                            ! Wraparound longitude
      REAL(wp) ::  zti, ztj, zarg                  ! Local scalars
      REAL(wp) ::  glam0, gphi0                    ! Variables corresponding to parameters ppglam0 ppgphi0 set in par_oce
      REAL(wp) ::  zlam1, zcos_alpha, ze1, ze1deg  ! Case 5 local scalars
      REAL(wp) ::  zphi1, zsin_alpha, zim05, zjm05 !          "

      REAL(wp) , POINTER, DIMENSION(:,:) ::  gphidta, glamdta, zdist ! Global lat/lon
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )   CALL timing_start('dom_c1d')

      REWIND( numnam_ref )              ! Namelist namdom in reference namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_ref, namdom, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdom in reference namelist', lwp )
  
      !
      REWIND( numnam_cfg )              ! Namelist namdom in configuration namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_cfg, namdom, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdom in configuration namelist', lwp )

      CALL wrk_alloc( jpidta, jpjdta, gphidta, glamdta, zdist )


      ! ============================= !
      !  Code from dom_hgr:           !
      !  Calculate global horizontal  !
      !  mesh, only glamt and gphit   !
      ! ============================= !

      SELECT CASE( jphgr_msh )   ! type of horizontal mesh

      CASE ( 0 )                 !  curvilinear coordinate on the sphere read in coordinate.nc file

         CALL iom_open( 'coordinates', inum )
         CALL iom_get( inum, jpdom_unknown, 'glamt', glamdta ) ! mig, mjg undefined at this point
         CALL iom_get( inum, jpdom_unknown, 'gphit', gphidta ) ! so use jpdom_unknown not jpdom_data
         CALL iom_close ( inum )

      CASE ( 1 )                 ! geographical mesh on the sphere with regular grid-spacing

         DO jj = 1, jpjdta
            DO ji = 1, jpidta
               zti = FLOAT( ji - 1 + nimpp - 1 )
               ztj = FLOAT( jj - 1 + njmpp - 1 )

               glamdta(ji,jj) = ppglam0 + ppe1_deg * zti
               gphidta(ji,jj) = ppgphi0 + ppe2_deg * ztj
            END DO
         END DO

      CASE ( 2:3 )               ! f- or beta-plane with regular grid-spacing
         
         glam0 = 0.e0
         gphi0 = - ppe2_m * 1.e-3

         DO jj = 1, jpjdta
            DO ji = 1, jpidta
               glamdta(ji,jj) = glam0 + ppe1_m * 1.e-3 * FLOAT( ji - 1 + nimpp - 1 )
               gphidta(ji,jj) = gphi0 + ppe2_m * 1.e-3 * FLOAT( jj - 1 + njmpp - 1 )
            END DO
         END DO

      CASE ( 4 )                 ! geographical mesh on the sphere, isotropic MERCATOR type

         IF( ppgphi0 == -90 )   CALL ctl_stop( ' Mercator grid cannot start at south pole !!!! ' )

         zarg = rpi / 4. - rpi / 180. * ppgphi0 / 2.
         ijeq = ABS( 180. / rpi * LOG( COS( zarg ) / SIN( zarg ) ) / ppe1_deg )
         IF( ppgphi0 > 0 )   ijeq = -ijeq

         DO jj = 1, jpjdta
            DO ji = 1, jpidta
               zti = FLOAT( ji - 1    + nimpp - 1 )
               ztj = FLOAT( jj - ijeq + njmpp - 1 )

               glamdta(ji,jj) = ppglam0 + ppe1_deg * zti
               gphidta(ji,jj) = 1. / rad * ASIN ( TANH( ppe1_deg * rad * ztj ) )
            END DO
         END DO

      CASE ( 5 )                 ! beta-plane with regular grid-spacing and rotated domain (GYRE configuration)
   
         zlam1 = -85
         zphi1 = 29
         ze1 = 106000. / FLOAT(jp_cfg)
 
         zsin_alpha = - SQRT( 2. ) / 2.
         zcos_alpha =   SQRT( 2. ) / 2.
         ze1deg = ze1 / (ra * rad)

         glam0 = zlam1 + zcos_alpha * ze1deg * FLOAT( jpjdta-2 ) ! Force global
         gphi0 = zphi1 + zsin_alpha * ze1deg * FLOAT( jpjdta-2 )

         DO jj = 1, jpjdta
            DO ji = 1, jpidta
               zim05 = FLOAT( ji + nimpp - 1 ) - 1.5
               zjm05 = FLOAT( jj + njmpp - 1 ) - 1.5

               glamdta(ji,jj) = glam0 + zim05 * ze1deg * zcos_alpha + zjm05 * ze1deg * zsin_alpha
               gphidta(ji,jj) = gphi0 - zim05 * ze1deg * zsin_alpha + zjm05 * ze1deg * zcos_alpha
            END DO
         END DO

      CASE DEFAULT

         WRITE(ctmp1,*) '          bad flag value for jphgr_msh = ', jphgr_msh
         CALL ctl_stop( ctmp1 )

      END SELECT


      ! ============================== !
      !  Code from dom_ngb:            !
      !  Calculate the nearest grid    !
      !  point to the given lat/lon &  !
      !  update jpizoom and jpjzoom    !
      ! ============================== !

      zlon         = MOD( plon         + 720., 360. )                                      ! plon    between    0 and 360
      glamdta(:,:) = MOD( glamdta(:,:) + 720., 360. )                                      ! glamdta between    0 and 360
      IF( zlon > 270. )   zlon = zlon - 360.                                               ! zlon    between  -90 and 270
      IF( zlon <  90. )   WHERE( glamdta(:,:) > 180. ) glamdta(:,:) = glamdta(:,:) - 360.  ! glamdta between -180 and 180

      glamdta(:,:) = glamdta(:,:) - zlon
      gphidta(:,:) = gphidta(:,:) - plat
      zdist(:,:)   = glamdta(:,:) * glamdta(:,:) + gphidta(:,:) * gphidta(:,:)
      
      iloc(:) = MINLOC( zdist(:,:) ) ! No mask; zoom indices freely defined
      jpizoom = iloc(1) + nimpp - 2  ! Minloc index - 1; want the bottom-left
      jpjzoom = iloc(2) + njmpp - 2  ! corner index of the zoom domain.

      CALL wrk_dealloc( jpidta, jpjdta, gphidta, glamdta, zdist )

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_c1d : compute jpizoom & jpjzoom from global mesh and given coordinates'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '      column i zoom index             jpizoom = ', jpizoom
         WRITE(numout,*) '      column j zoom index             jpjzoom = ', jpjzoom
         WRITE(numout,*)
      ENDIF

      IF( nn_timing == 1 )   CALL timing_stop('dom_c1d')

   END SUBROUTINE dom_c1d

#else
   !!----------------------------------------------------------------------
   !!   Default option                                  NO 1D Configuration
   !!----------------------------------------------------------------------
CONTAINS  
   SUBROUTINE dom_c1d( plat, plon )     ! Empty routine
      REAL :: plat, plon
      WRITE(*,*) 'dom_c1d: You should not have seen this print! error?',plat,plon
   END SUBROUTINE dom_c1d
#endif

   !!======================================================================
END MODULE domc1d
