MODULE crsdomwri
   !!======================================================================
   !! Coarse Ocean initialization : write the coarse ocean domain mesh and mask files
   !!======================================================================
   !! History :  OPA  ! 1997-02  (G. Madec)  Original code
   !!            8.1  ! 1999-11  (M. Imbard)  NetCDF FORMAT with IOIPSL
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90 and several file
   !!            3.0  ! 2008-01  (S. Masson) add dom_uniq_crs
   !!            4.0  ! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!                 ! 2012-06  (J. Simeon, C. Calone, C Ethe )  Reduced and modified for coarse grid
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   crs_dom_wri    : create and write mesh and mask file(s)
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE par_kind, ONLY: wp
   USE lib_mpp         ! MPP library
   USE iom_def
   USE iom
   USE crs         ! coarse grid domain
   USE crsdom         ! coarse grid domain
   USE crslbclnk       ! crs mediator to lbclnk
   USE wrk_nemo        ! Working array



   IMPLICIT NONE
   PRIVATE

   PUBLIC crs_dom_wri        ! routine called by crsini.F90

   !! $Id$
CONTAINS

   SUBROUTINE crs_dom_wri
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE crs_dom_wri  ***
      !!
      !! ** Purpose :   Create the NetCDF file(s) which contain(s) all the
      !!      ocean domain informations (mesh and mask arrays). This (these)
      !!      file(s) is (are) used for visualisation (SAXO software) and
      !!      diagnostic computation.
      !!
      !! ** Method  :   Write in a file all the arrays generated in routines
      !!      crsini for meshes and mask. In three separate files: 
      !!      domain size, horizontal grid-point position,
      !!      masks, depth and vertical scale factors
      !!      
      !! ** Output files :   mesh_hgr_crs.nc, mesh_zgr_crs.nc, mesh_mask.nc
      !!----------------------------------------------------------------------
      !!
      INTEGER           ::   inum0    ! temprary units for 'mesh_mask.nc' file
      INTEGER           ::   inum1    ! temprary units for 'mesh.nc'      file
      INTEGER           ::   inum2    ! temprary units for 'mask.nc'      file
      INTEGER           ::   inum3    ! temprary units for 'mesh_hgr.nc'  file
      INTEGER           ::   inum4    ! temprary units for 'mesh_zgr.nc'  file
      INTEGER           ::   iif, iil, ijf, ijl
      CHARACTER(len=21) ::   clnam0   ! filename (mesh and mask informations)
      CHARACTER(len=21) ::   clnam1   ! filename (mesh informations)
      CHARACTER(len=21) ::   clnam2   ! filename (mask informations)
      CHARACTER(len=21) ::   clnam3   ! filename (horizontal mesh informations)
      CHARACTER(len=21) ::   clnam4   ! filename (vertical   mesh informations)
      INTEGER           ::   ji, jj, jk   ! dummy loop indices
      !                                   !  workspaces
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zprt, zprw 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdepu, zdepv
      REAL(wp), POINTER, DIMENSION(:,:  ) :: ze3tp, ze3wp
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('crs_dom_wri')
      !
      CALL wrk_alloc( jpi_crs, jpj_crs,      zprt , zprw  )
      CALL wrk_alloc( jpi_crs, jpj_crs,      ze3tp, ze3wp )
      CALL wrk_alloc( jpi_crs, jpj_crs, jpk, zdepu, zdepv )

      ze3tp(:,:) = 0.0
      ze3wp(:,:) = 0.0

      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'crs_dom_wri : create NetCDF mesh and mask information file(s)'
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      
      clnam0 = 'mesh_mask_crs'  ! filename (mesh and mask informations)
      clnam1 = 'mesh_crs'       ! filename (mesh informations)
      clnam2 = 'mask_crs'       ! filename (mask informations)
      clnam3 = 'mesh_hgr_crs'   ! filename (horizontal mesh informations)
      clnam4 = 'mesh_zgr_crs'   ! filename (vertical   mesh informations)
      

      SELECT CASE ( MOD(nn_msh_crs, 3) )
         !                                  ! ============================
      CASE ( 1 )                            !  create 'mesh_mask.nc' file
         !                                  ! ============================
         CALL iom_open( TRIM(clnam0), inum0, ldwrt = .TRUE., kiolib = jprstlib )
         inum2 = inum0                                            ! put all the informations
         inum3 = inum0                                            ! in unit inum0
         inum4 = inum0
         
         !                                  ! ============================
      CASE ( 2 )                            !  create 'mesh.nc' and 
         !                                  !         'mask.nc' files
         !                                  ! ============================
         CALL iom_open( TRIM(clnam1), inum1, ldwrt = .TRUE., kiolib = jprstlib )
         CALL iom_open( TRIM(clnam2), inum2, ldwrt = .TRUE., kiolib = jprstlib )
         inum3 = inum1                                            ! put mesh informations 
         inum4 = inum1                                            ! in unit inum1 
         !                                  ! ============================
      CASE ( 0 )                            !  create 'mesh_hgr.nc'
         !                                  !         'mesh_zgr.nc' and
         !                                  !         'mask.nc'     files
         !                                  ! ============================
         CALL iom_open( TRIM(clnam2), inum2, ldwrt = .TRUE., kiolib = jprstlib )
         CALL iom_open( TRIM(clnam3), inum3, ldwrt = .TRUE., kiolib = jprstlib )
         CALL iom_open( TRIM(clnam4), inum4, ldwrt = .TRUE., kiolib = jprstlib )
         !
      END SELECT
 
      !========================================================
      !                                                         ! masks (inum2) 
      CALL iom_rstput( 0, 0, inum2, 'tmask', tmask_crs, ktype = jp_i1 )     !    ! land-sea mask
      CALL iom_rstput( 0, 0, inum2, 'umask', umask_crs, ktype = jp_i1 )
      CALL iom_rstput( 0, 0, inum2, 'vmask', vmask_crs, ktype = jp_i1 )
      CALL iom_rstput( 0, 0, inum2, 'fmask', fmask_crs, ktype = jp_i1 )
      
      
      tmask_i_crs(:,:) = tmask_crs(:,:,1)
      iif = jpreci
      iil = nlci_crs - jpreci + 1
      ijf = jpreci
      ijl = nlcj_crs - jprecj + 1
     
      tmask_i_crs( 1:iif ,    :  ) = 0._wp
      tmask_i_crs(iil:jpi_crs,    :  ) = 0._wp
      tmask_i_crs(   :   , 1:ijf ) = 0._wp
      tmask_i_crs(   :   ,ijl:jpj_crs) = 0._wp
      
      
      tpol_crs(1:jpiglo_crs,:) = 1._wp
      fpol_crs(1:jpiglo_crs,:) = 1._wp
      IF( jperio == 3 .OR. jperio == 4 ) THEN
         tpol_crs(jpiglo_crs/2+1:jpiglo_crs,:) = 0._wp
         fpol_crs(       1      :jpiglo_crs,:) = 0._wp
         IF( mjg_crs(nlej_crs) == jpiglo_crs ) THEN
            DO ji = iif+1, iil-1
               tmask_i_crs(ji,nlej_crs-1) = tmask_i_crs(ji,nlej_crs-1) &
               & * tpol_crs(mig_crs(ji),1)
            ENDDO
         ENDIF
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN
         tpol_crs(      1       :jpiglo_crs,:)=0._wp
         fpol_crs(jpiglo_crs/2+1:jpiglo_crs,:)=0._wp
      ENDIF
      
      CALL iom_rstput( 0, 0, inum2, 'tmaskutil', tmask_i_crs, ktype = jp_i1 )
                                   !    ! unique point mask
      CALL dom_uniq_crs( zprw, 'U' )
      zprt = umask_crs(:,:,1) * zprw
      CALL iom_rstput( 0, 0, inum2, 'umaskutil', zprt, ktype = jp_i1 )  
      CALL dom_uniq_crs( zprw, 'V' )
      zprt = vmask_crs(:,:,1) * zprw
      CALL iom_rstput( 0, 0, inum2, 'vmaskutil', zprt, ktype = jp_i1 )  
      CALL dom_uniq_crs( zprw, 'F' )
      zprt = fmask_crs(:,:,1) * zprw
      CALL iom_rstput( 0, 0, inum2, 'fmaskutil', zprt, ktype = jp_i1 )  
      !========================================================
      !                                                         ! horizontal mesh (inum3)
      CALL iom_rstput( 0, 0, inum3, 'glamt', glamt_crs, ktype = jp_r4 )     !    ! latitude
      CALL iom_rstput( 0, 0, inum3, 'glamu', glamu_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum3, 'glamv', glamv_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum3, 'glamf', glamf_crs, ktype = jp_r4 )
      
      CALL iom_rstput( 0, 0, inum3, 'gphit', gphit_crs, ktype = jp_r4 )     !    ! longitude
      CALL iom_rstput( 0, 0, inum3, 'gphiu', gphiu_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum3, 'gphiv', gphiv_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum3, 'gphif', gphif_crs, ktype = jp_r4 )
      
      CALL iom_rstput( 0, 0, inum3, 'e1t', e1t_crs, ktype = jp_r8 )         !    ! e1 scale factors
      CALL iom_rstput( 0, 0, inum3, 'e1u', e1u_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum3, 'e1v', e1v_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum3, 'e1f', e1f_crs, ktype = jp_r8 )
      
      CALL iom_rstput( 0, 0, inum3, 'e2t', e2t_crs, ktype = jp_r8 )         !    ! e2 scale factors
      CALL iom_rstput( 0, 0, inum3, 'e2u', e2u_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum3, 'e2v', e2v_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum3, 'e2f', e2f_crs, ktype = jp_r8 )
      
      CALL iom_rstput( 0, 0, inum3, 'ff', ff_crs, ktype = jp_r8 )           !    ! coriolis factor

      !========================================================
      !                                                         ! vertical mesh (inum4) 
!     ! note that mbkt is set to 1 over land ==> use surface tmask_crs
      zprt(:,:) = tmask_crs(:,:,1) * REAL( mbkt_crs(:,:) , wp )
      CALL iom_rstput( 0, 0, inum4, 'mbathy', zprt, ktype = jp_i2 )     !    ! nb of ocean T-points

      IF( ln_zps ) THEN                       ! z-coordinate - partial steps

            
         IF ( nn_msh_crs <= 6 ) THEN
            CALL iom_rstput( 0, 0, inum4, 'e3t', e3t_crs )      
            CALL iom_rstput( 0, 0, inum4, 'e3w', e3w_crs )      
            CALL iom_rstput( 0, 0, inum4, 'e3u', e3u_crs )      
            CALL iom_rstput( 0, 0, inum4, 'e3v', e3v_crs )      
         ELSE
            DO jj = 1,jpj_crs   
               DO ji = 1,jpi_crs
                  ze3tp(ji,jj) = e3t_crs(ji,jj,mbkt_crs(ji,jj)) * tmask_crs(ji,jj,1)
                  ze3wp(ji,jj) = e3w_crs(ji,jj,mbkt_crs(ji,jj)) * tmask_crs(ji,jj,1)
               END DO
            END DO

            CALL crs_lbc_lnk( ze3tp,'T', 1.0 )
            CALL crs_lbc_lnk( ze3wp,'W', 1.0 )
  
            CALL iom_rstput( 0, 0, inum4, 'e3t_ps', ze3tp )      
            CALL iom_rstput( 0, 0, inum4, 'e3w_ps', ze3wp )
         ENDIF

         IF ( nn_msh_crs <= 3 ) THEN
            CALL iom_rstput( 0, 0, inum4, 'gdept', gdept_crs, ktype = jp_r4 ) 
            DO jk = 1,jpk   
               DO jj = 1, jpj_crsm1   
                  DO ji = 1, jpi_crsm1  ! jes what to do for fs_jpim1??vector opt.
                     zdepu(ji,jj,jk) = MIN( gdept_crs(ji,jj,jk) , gdept_crs(ji+1,jj  ,jk) ) * umask_crs(ji,jj,jk)
                     zdepv(ji,jj,jk) = MIN( gdept_crs(ji,jj,jk) , gdept_crs(ji  ,jj+1,jk) ) * vmask_crs(ji,jj,jk)
                  END DO   
               END DO   
            END DO

            CALL crs_lbc_lnk( zdepu,'U', 1. )   ;   CALL crs_lbc_lnk( zdepv,'V', 1. ) 
            CALL iom_rstput( 0, 0, inum4, 'gdepu', zdepu, ktype = jp_r4 )
            CALL iom_rstput( 0, 0, inum4, 'gdepv', zdepv, ktype = jp_r4 )
            CALL iom_rstput( 0, 0, inum4, 'gdepw', gdepw_crs, ktype = jp_r4 )
         ELSE
            DO jj = 1,jpj_crs   
               DO ji = 1,jpi_crs
                  zprt(ji,jj) = gdept_0(ji,jj,mbkt(ji,jj)  ) * tmask(ji,jj,1)
                  zprw(ji,jj) = gdepw_0(ji,jj,mbkt(ji,jj)+1) * tmask(ji,jj,1)
               END DO
            END DO
            CALL iom_rstput( 0, 0, inum4, 'hdept', zprt, ktype = jp_r4 )     
            CALL iom_rstput( 0, 0, inum4, 'hdepw', zprw, ktype = jp_r4 ) 
         ENDIF

         CALL iom_rstput( 0, 0, inum4, 'gdept_1d', gdept_1d )     !    ! reference z-coord.
         CALL iom_rstput( 0, 0, inum4, 'gdepw_1d', gdepw_1d )
         CALL iom_rstput( 0, 0, inum4, 'e3t_1d'  , e3t_1d   )
         CALL iom_rstput( 0, 0, inum4, 'e3w_1d'  , e3w_1d   )

         CALL iom_rstput(  0, 0, inum4, 'ocean_volume_t', ocean_volume_crs_t ) 
         CALL iom_rstput(  0, 0, inum4, 'facvol_t' , facvol_t  ) 
         CALL iom_rstput(  0, 0, inum4, 'facvol_w' , facvol_w  ) 
         CALL iom_rstput(  0, 0, inum4, 'facsurfu' , facsurfu  ) 
         CALL iom_rstput(  0, 0, inum4, 'facsurfv' , facsurfv  ) 
         CALL iom_rstput(  0, 0, inum4, 'e1e2w_msk', e1e2w_msk ) 
         CALL iom_rstput(  0, 0, inum4, 'e2e3u_msk', e2e3u_msk ) 
         CALL iom_rstput(  0, 0, inum4, 'e1e3v_msk', e1e3v_msk )
         CALL iom_rstput(  0, 0, inum4, 'e1e2w'    , e1e2w_crs ) 
         CALL iom_rstput(  0, 0, inum4, 'e2e3u'    , e2e3u_crs ) 
         CALL iom_rstput(  0, 0, inum4, 'e1e3v'    , e1e3v_crs )
         CALL iom_rstput(  0, 0, inum4, 'bt'       , bt_crs    )
         CALL iom_rstput(  0, 0, inum4, 'r1_bt'    , r1_bt_crs )

         CALL iom_rstput(  0, 0, inum4, 'crs_surfu_wgt', crs_surfu_wgt ) 
         CALL iom_rstput(  0, 0, inum4, 'crs_surfv_wgt', crs_surfv_wgt ) 
         CALL iom_rstput(  0, 0, inum4, 'crs_volt_wgt' , crs_volt_wgt  ) 

      ENDIF
      
     IF( ln_zco ) THEN
         !                                                      ! z-coordinate - full steps
        CALL iom_rstput( 0, 0, inum4, 'gdept_1d', gdept_1d )     !    ! depth
        CALL iom_rstput( 0, 0, inum4, 'gdepw_1d', gdepw_1d )
        CALL iom_rstput( 0, 0, inum4, 'e3t_1d'  , e3t_1d   )     !    ! scale factors
        CALL iom_rstput( 0, 0, inum4, 'e3w_1d'  , e3w_1d   )
     ENDIF
      !                                     ! ============================
      !                                     !        close the files 
      !                                     ! ============================
      SELECT CASE ( MOD(nn_msh_crs, 3) )
      CASE ( 1 )                
         CALL iom_close( inum0 )
      CASE ( 2 )
         CALL iom_close( inum1 )
         CALL iom_close( inum2 )
      CASE ( 0 )
         CALL iom_close( inum2 )
         CALL iom_close( inum3 )
         CALL iom_close( inum4 )
      END SELECT
      !
      CALL wrk_dealloc( jpi_crs, jpj_crs,      zprt , zprw  )
      CALL wrk_dealloc( jpi_crs, jpj_crs,      ze3tp, ze3wp )
      CALL wrk_dealloc( jpi_crs, jpj_crs, jpk, zdepu, zdepv )
      !
      IF( nn_timing == 1 )  CALL timing_stop('crs_dom_wri')
      !
      
   END SUBROUTINE crs_dom_wri


   SUBROUTINE dom_uniq_crs( puniq, cdgrd )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE crs_dom_uniq_crs  ***
      !!                   
      !! ** Purpose :   identify unique point of a grid (TUVF)
      !!
      !! ** Method  :   1) apply crs_lbc_lnk on an array with different values for each element
      !!                2) check which elements have been changed
      !!----------------------------------------------------------------------
      !
      CHARACTER(len=1)        , INTENT(in   ) ::   cdgrd   ! 
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   puniq   ! 
      !
      REAL(wp) ::  zshift   ! shift value link to the process number
      INTEGER  ::  ji       ! dummy loop indices
      LOGICAL, DIMENSION(SIZE(puniq,1),SIZE(puniq,2),1) ::  lldbl  ! store whether each point is unique or not
      REAL(wp), POINTER, DIMENSION(:,:) :: ztstref
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('crs_dom_uniq_crs')
      !
      CALL wrk_alloc( jpi_crs, jpj_crs, ztstref )
      !
      ! build an array with different values for each element 
      ! in mpp: make sure that these values are different even between process
      ! -> apply a shift value according to the process number
      zshift = jpi_crs * jpj_crs * ( narea - 1 )
      ztstref(:,:) = RESHAPE( (/ (zshift + REAL(ji,wp), ji = 1, jpi_crs*jpj_crs) /), (/ jpi_crs, jpj_crs /) )
      !
      puniq(:,:) = ztstref(:,:)                   ! default definition
      CALL crs_lbc_lnk( puniq,cdgrd, 1. )            ! apply boundary conditions
      lldbl(:,:,1) = puniq(:,:) == ztstref(:,:)   ! check which values have been changed 
      !
      puniq(:,:) = 1.                             ! default definition
      ! fill only the inner part of the cpu with llbl converted into real 
      puniq(nldi_crs:nlei_crs,nldj_crs:nlej_crs) = REAL( COUNT( lldbl(nldi_crs:nlei_crs,nldj_crs:nlej_crs,:), dim = 3 ) , wp )
      !
      CALL wrk_dealloc( jpi_crs, jpj_crs, ztstref )
      !
      IF( nn_timing == 1 )  CALL timing_stop('crs_dom_uniq_crs')
      !
      
   END SUBROUTINE dom_uniq_crs

   !!======================================================================

END MODULE crsdomwri


