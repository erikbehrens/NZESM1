MODULE sbcwave
   !!======================================================================
   !!                       ***  MODULE  sbcwave  ***
   !! Wave module 
   !!======================================================================
   !! History :  3.3.1  !   2011-09  (Adani M)  Original code: Drag Coefficient 
   !!         :  3.4    !   2012-10  (Adani M)                 Stokes Drift 
   !!----------------------------------------------------------------------
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE fldread	       ! read input fields
   USE oce
   USE sbc_oce	       ! Surface boundary condition: ocean fields
   USE domvvl

   
   !!----------------------------------------------------------------------
   !!   sbc_wave       : read drag coefficient from wave model in netcdf files 
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_wave    ! routine called in sbc_blk_core or sbc_blk_mfs
   
   INTEGER , PARAMETER ::   jpfld  = 3           ! maximum number of files to read for srokes drift
   INTEGER , PARAMETER ::   jp_usd = 1           ! index of stokes drift  (i-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_vsd = 2           ! index of stokes drift  (j-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_wn  = 3           ! index of wave number                 (1/m)    at T-point
   TYPE(FLD), ALLOCATABLE, DIMENSION(:)  :: sf_cd	  ! structure of input fields (file informations, fields read) Drag Coefficient
   TYPE(FLD), ALLOCATABLE, DIMENSION(:)  :: sf_sd	  ! structure of input fields (file informations, fields read) Stokes Drift
   REAL(wp),PUBLIC,ALLOCATABLE,DIMENSION (:,:)       :: cdn_wave 
   REAL(wp),ALLOCATABLE,DIMENSION (:,:)              :: usd2d,vsd2d,uwavenum,vwavenum 
   REAL(wp),PUBLIC,ALLOCATABLE,DIMENSION (:,:,:)     :: usd3d,vsd3d,wsd3d 

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011) 
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_wave( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_apr  ***
      !!
      !! ** Purpose :   read drag coefficient from wave model  in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - Read Cd_n10 fields in netcdf files 
      !!              - Read stokes drift 2d in netcdf files 
      !!              - Read wave number      in netcdf files 
      !!              - Compute 3d stokes drift using monochromatic
      !! ** action  :   
      !!               
      !!---------------------------------------------------------------------
      USE oce,  ONLY : un,vn,hdivn,rotn
      USE divcur
      USE wrk_nemo
#if defined key_bdy
      USE bdy_oce, ONLY : bdytmask
#endif
      INTEGER, INTENT( in  ) ::  kt       ! ocean time step
      INTEGER                ::  ierror   ! return error code
      INTEGER                ::  ifpr, jj,ji,jk 
      INTEGER                ::   ios     ! Local integer output status for namelist read
      REAL(wp),DIMENSION(:,:,:),POINTER             ::  udummy,vdummy,hdivdummy,rotdummy
      REAL                                          ::  z2dt,z1_2dt
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i     ! array of namelist informations on the fields to read
      CHARACTER(len=100)     ::  cn_dir                          ! Root directory for location of drag coefficient files
      TYPE(FLD_N)            ::  sn_cdg, sn_usd, sn_vsd, sn_wn   ! informations about the fields to be read
      !!---------------------------------------------------------------------
      NAMELIST/namsbc_wave/  sn_cdg, cn_dir, sn_usd, sn_vsd, sn_wn
      !!---------------------------------------------------------------------

      !!----------------------------------------------------------------------
      !
      !
      !                                         ! -------------------- !
      IF( kt == nit000 ) THEN                   ! First call kt=nit000 !
         !                                      ! -------------------- !
         REWIND( numnam_ref )              ! Namelist namsbc_wave in reference namelist : File for drag coeff. from wave model
         READ  ( numnam_ref, namsbc_wave, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_wave in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_wave in configuration namelist : File for drag coeff. from wave model
         READ  ( numnam_cfg, namsbc_wave, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_wave in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_wave )
         !

         IF ( ln_cdgw ) THEN
            ALLOCATE( sf_cd(1), STAT=ierror )           !* allocate and fill sf_wave with sn_cdg
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave: unable to allocate sf_wave structure' )
            !
                                   ALLOCATE( sf_cd(1)%fnow(jpi,jpj,1)   )
            IF( sn_cdg%ln_tint )   ALLOCATE( sf_cd(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_cd, (/ sn_cdg /), cn_dir, 'sbc_wave', 'Wave module ', 'namsbc_wave' )
            ALLOCATE( cdn_wave(jpi,jpj) )
            cdn_wave(:,:) = 0.0
        ENDIF
         IF ( ln_sdw ) THEN
            slf_i(jp_usd) = sn_usd ; slf_i(jp_vsd) = sn_vsd; slf_i(jp_wn) = sn_wn
            ALLOCATE( sf_sd(3), STAT=ierror )           !* allocate and fill sf_wave with sn_cdg
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave: unable to allocate sf_wave structure' )
            !
            DO ifpr= 1, jpfld
               ALLOCATE( sf_sd(ifpr)%fnow(jpi,jpj,1) )
               IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf_sd(ifpr)%fdta(jpi,jpj,1,2) )
            END DO
            CALL fld_fill( sf_sd, slf_i, cn_dir, 'sbc_wave', 'Wave module ', 'namsbc_wave' )
            ALLOCATE( usd2d(jpi,jpj),vsd2d(jpi,jpj),uwavenum(jpi,jpj),vwavenum(jpi,jpj) )
            ALLOCATE( usd3d(jpi,jpj,jpk),vsd3d(jpi,jpj,jpk),wsd3d(jpi,jpj,jpk) )
            usd2d(:,:) = 0.0 ;  vsd2d(:,:) = 0.0 ; uwavenum(:,:) = 0.0 ; vwavenum(:,:) = 0.0
            usd3d(:,:,:) = 0.0 ;vsd3d(:,:,:) = 0.0 ; wsd3d(:,:,:) = 0.0
         ENDIF
      ENDIF
         !
         !
      IF ( ln_cdgw ) THEN
         CALL fld_read( kt, nn_fsbc, sf_cd )      !* read drag coefficient from external forcing
         cdn_wave(:,:) = sf_cd(1)%fnow(:,:,1)
      ENDIF
      IF ( ln_sdw )  THEN
          CALL fld_read( kt, nn_fsbc, sf_sd )      !* read drag coefficient from external forcing

         ! Interpolate wavenumber, stokes drift into the grid_V and grid_V
         !-------------------------------------------------

         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               uwavenum(ji,jj)=0.5 * ( 2. - umask(ji,jj,1) ) * ( sf_sd(3)%fnow(ji,jj,1) * tmask(ji,jj,1) &
               &                                + sf_sd(3)%fnow(ji+1,jj,1) * tmask(ji+1,jj,1) )

               vwavenum(ji,jj)=0.5 * ( 2. - vmask(ji,jj,1) ) * ( sf_sd(3)%fnow(ji,jj,1) * tmask(ji,jj,1) &
               &                                + sf_sd(3)%fnow(ji,jj+1,1) * tmask(ji,jj+1,1) )

               usd2d(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( sf_sd(1)%fnow(ji,jj,1) * tmask(ji,jj,1) &
               &                                + sf_sd(1)%fnow(ji+1,jj,1) * tmask(ji+1,jj,1) )

               vsd2d(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( sf_sd(2)%fnow(ji,jj,1) * tmask(ji,jj,1) &
               &                                + sf_sd(2)%fnow(ji,jj+1,1) * tmask(ji,jj+1,1) )
            END DO
         END DO

          !Computation of the 3d Stokes Drift
          DO jk = 1, jpk
             DO jj = 1, jpj-1
                DO ji = 1, jpi-1
                   usd3d(ji,jj,jk) = usd2d(ji,jj)*exp(2.0*uwavenum(ji,jj)*(-MIN( gdept_0(ji,jj,jk) , gdept_0(ji+1,jj  ,jk))))
                   vsd3d(ji,jj,jk) = vsd2d(ji,jj)*exp(2.0*vwavenum(ji,jj)*(-MIN( gdept_0(ji,jj,jk) , gdept_0(ji  ,jj+1,jk))))
                END DO
             END DO
             usd3d(jpi,:,jk) = usd2d(jpi,:)*exp( 2.0*uwavenum(jpi,:)*(-gdept_0(jpi,:,jk)) )
             vsd3d(:,jpj,jk) = vsd2d(:,jpj)*exp( 2.0*vwavenum(:,jpj)*(-gdept_0(:,jpj,jk)) )
          END DO

          CALL wrk_alloc( jpi,jpj,jpk,udummy,vdummy,hdivdummy,rotdummy)
          
          udummy(:,:,:)=un(:,:,:)
          vdummy(:,:,:)=vn(:,:,:)
          hdivdummy(:,:,:)=hdivn(:,:,:)
          rotdummy(:,:,:)=rotn(:,:,:)
          un(:,:,:)=usd3d(:,:,:)
          vn(:,:,:)=vsd3d(:,:,:)
          CALL div_cur(kt)
      !                                           !------------------------------!
      !                                           !     Now Vertical Velocity    !
      !                                           !------------------------------!
          z2dt = 2._wp * rdt                              ! set time step size (Euler/Leapfrog)

          z1_2dt = 1.e0 / z2dt
          DO jk = jpkm1, 1, -1                             ! integrate from the bottom the hor. divergence
             ! - ML - need 3 lines here because replacement of fse3t by its expression yields too long lines otherwise
             wsd3d(:,:,jk) = wsd3d(:,:,jk+1) -   fse3t_n(:,:,jk) * hdivn(:,:,jk)        &
                &                      - ( fse3t_a(:,:,jk) - fse3t_b(:,:,jk) )    &
                &                         * tmask(:,:,jk) * z1_2dt
#if defined key_bdy
             wsd3d(:,:,jk) = wsd3d(:,:,jk) * bdytmask(:,:)
#endif
          END DO
          hdivn(:,:,:)=hdivdummy(:,:,:)
          rotn(:,:,:)=rotdummy(:,:,:)
          vn(:,:,:)=vdummy(:,:,:)
          un(:,:,:)=udummy(:,:,:)
          CALL wrk_dealloc( jpi,jpj,jpk,udummy,vdummy,hdivdummy,rotdummy)
      ENDIF
   END SUBROUTINE sbc_wave
      
   !!======================================================================
END MODULE sbcwave
