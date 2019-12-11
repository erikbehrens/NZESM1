MODULE bdytra
   !!======================================================================
   !!                       ***  MODULE  bdytra  ***
   !! Ocean tracers:   Apply boundary conditions for tracers
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!----------------------------------------------------------------------
#if defined key_bdy
   !!----------------------------------------------------------------------
   !!   'key_bdy'                     Unstructured Open Boundary Conditions
   !!----------------------------------------------------------------------
   !!   bdy_tra            : Apply open boundary conditions to T and S
   !!   bdy_tra_frs        : Apply Flow Relaxation Scheme
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE bdy_oce         ! ocean open boundary conditions
   USE bdylib          ! for orlanski library routines
   USE bdydta, ONLY:   bf
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC bdy_tra      ! routine called in tranxt.F90 
   PUBLIC bdy_tra_dmp  ! routine called in step.F90 

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_tra( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_tra  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for temperature and salinity
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: kt     ! Main time step counter
      !!
      INTEGER               :: ib_bdy ! Loop index

      DO ib_bdy=1, nb_bdy

         SELECT CASE( cn_tra(ib_bdy) )
         CASE('none')
            CYCLE
         CASE('frs')
            CALL bdy_tra_frs( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt )
         CASE('specified')
            CALL bdy_tra_spe( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt )
         CASE('neumann')
            CALL bdy_tra_nmn( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt )
         CASE('orlanski')
            CALL bdy_tra_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ll_npo=.false. )
         CASE('orlanski_npo')
            CALL bdy_tra_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ll_npo=.true. )
         CASE('runoff')
            CALL bdy_tra_rnf( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt )
         CASE DEFAULT
            CALL ctl_stop( 'bdy_tra : unrecognised option for open boundaries for T and S' )
         END SELECT
         ! Boundary points should be updated
         CALL lbc_bdy_lnk( tsa(:,:,:,jp_tem), 'T', 1., ib_bdy )
         CALL lbc_bdy_lnk( tsa(:,:,:,jp_sal), 'T', 1., ib_bdy )
      ENDDO
      !

   END SUBROUTINE bdy_tra

   SUBROUTINE bdy_tra_frs( idx, dta, kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_frs  ***
      !!                    
      !! ** Purpose : Apply the Flow Relaxation Scheme for tracers at open boundaries.
      !! 
      !! Reference : Engedahl H., 1995, Tellus, 365-382.
      !!----------------------------------------------------------------------
      INTEGER,         INTENT(in) ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      !! 
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_tra_frs')
      !
      igrd = 1                       ! Everything is at T-points here
      DO ib = 1, idx%nblen(igrd)
         DO ik = 1, jpkm1
            ii = idx%nbi(ib,igrd)
            ij = idx%nbj(ib,igrd)
            zwgt = idx%nbw(ib,igrd)
            tsa(ii,ij,ik,jp_tem) = ( tsa(ii,ij,ik,jp_tem) + zwgt * ( dta%tem(ib,ik) - tsa(ii,ij,ik,jp_tem) ) ) * tmask(ii,ij,ik)         
            tsa(ii,ij,ik,jp_sal) = ( tsa(ii,ij,ik,jp_sal) + zwgt * ( dta%sal(ib,ik) - tsa(ii,ij,ik,jp_sal) ) ) * tmask(ii,ij,ik)
         END DO
      END DO 
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_tra_frs')
      !
   END SUBROUTINE bdy_tra_frs
  
   SUBROUTINE bdy_tra_spe( idx, dta, kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_frs  ***
      !!                    
      !! ** Purpose : Apply a specified value for tracers at open boundaries.
      !! 
      !!----------------------------------------------------------------------
      INTEGER,         INTENT(in) ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      !! 
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_tra_spe')
      !
      igrd = 1                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            tsa(ii,ij,ik,jp_tem) = dta%tem(ib,ik) * tmask(ii,ij,ik)
            tsa(ii,ij,ik,jp_sal) = dta%sal(ib,ik) * tmask(ii,ij,ik)
         END DO
      END DO
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_tra_spe')
      !
   END SUBROUTINE bdy_tra_spe

   SUBROUTINE bdy_tra_nmn( idx, dta, kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_nmn  ***
      !!                    
      !! ** Purpose : Duplicate the value for tracers at open boundaries.
      !! 
      !!----------------------------------------------------------------------
      INTEGER,         INTENT(in) ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      !! 
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij,zcoef, zcoef1,zcoef2, ip, jp   ! 2D addresses
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_tra_nmn')
      !
      igrd = 1                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            ! search the sense of the gradient
            zcoef1 = bdytmask(ii-1,ij  ) +  bdytmask(ii+1,ij  )
            zcoef2 = bdytmask(ii  ,ij-1) +  bdytmask(ii  ,ij+1)
            IF ( zcoef1+zcoef2 == 0) THEN
               ! corner
               zcoef = tmask(ii-1,ij,ik) + tmask(ii+1,ij,ik) +  tmask(ii,ij-1,ik) +  tmask(ii,ij+1,ik)
               tsa(ii,ij,ik,jp_tem) = tsa(ii-1,ij  ,ik,jp_tem) * tmask(ii-1,ij  ,ik) + &
                 &                    tsa(ii+1,ij  ,ik,jp_tem) * tmask(ii+1,ij  ,ik) + &
                 &                    tsa(ii  ,ij-1,ik,jp_tem) * tmask(ii  ,ij-1,ik) + &
                 &                    tsa(ii  ,ij+1,ik,jp_tem) * tmask(ii  ,ij+1,ik)
               tsa(ii,ij,ik,jp_tem) = ( tsa(ii,ij,ik,jp_tem) / MAX( 1, zcoef) ) * tmask(ii,ij,ik)
               tsa(ii,ij,ik,jp_sal) = tsa(ii-1,ij  ,ik,jp_sal) * tmask(ii-1,ij  ,ik) + &
                 &                    tsa(ii+1,ij  ,ik,jp_sal) * tmask(ii+1,ij  ,ik) + &
                 &                    tsa(ii  ,ij-1,ik,jp_sal) * tmask(ii  ,ij-1,ik) + &
                 &                    tsa(ii  ,ij+1,ik,jp_sal) * tmask(ii  ,ij+1,ik)
               tsa(ii,ij,ik,jp_sal) = ( tsa(ii,ij,ik,jp_sal) / MAX( 1, zcoef) ) * tmask(ii,ij,ik)
            ELSE
               ip = bdytmask(ii+1,ij  ) - bdytmask(ii-1,ij  )
               jp = bdytmask(ii  ,ij+1) - bdytmask(ii  ,ij-1)
               tsa(ii,ij,ik,jp_tem) = tsa(ii+ip,ij+jp,ik,jp_tem) * tmask(ii+ip,ij+jp,ik)
               tsa(ii,ij,ik,jp_sal) = tsa(ii+ip,ij+jp,ik,jp_sal) * tmask(ii+ip,ij+jp,ik)
            ENDIF
         END DO
      END DO
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_tra_nmn')
      !
   END SUBROUTINE bdy_tra_nmn
 

   SUBROUTINE bdy_tra_orlanski( idx, dta, ll_npo )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_orlanski  ***
      !!             
      !!              - Apply Orlanski radiation to temperature and salinity. 
      !!              - Wrapper routine for bdy_orlanski_3d
      !! 
      !!
      !! References:  Marchesiello, McWilliams and Shchepetkin, Ocean Modelling vol. 3 (2001)    
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),               INTENT(in) ::   dta  ! OBC external data
      LOGICAL,                      INTENT(in) ::   ll_npo  ! switch for NPO version

      INTEGER  ::   igrd                                    ! grid index
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 ) CALL timing_start('bdy_tra_orlanski')
      !
      igrd = 1      ! Orlanski bc on temperature; 
      !            
      CALL bdy_orlanski_3d( idx, igrd, tsb(:,:,:,jp_tem), tsa(:,:,:,jp_tem), dta%tem, ll_npo )

      igrd = 1      ! Orlanski bc on salinity;
      !  
      CALL bdy_orlanski_3d( idx, igrd, tsb(:,:,:,jp_sal), tsa(:,:,:,jp_sal), dta%sal, ll_npo )
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_tra_orlanski')
      !

   END SUBROUTINE bdy_tra_orlanski


   SUBROUTINE bdy_tra_rnf( idx, dta, kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_rnf  ***
      !!                    
      !! ** Purpose : Apply the runoff values for tracers at open boundaries:
      !!                  - specified to 0.1 PSU for the salinity
      !!                  - duplicate the value for the temperature
      !! 
      !!----------------------------------------------------------------------
      INTEGER,         INTENT(in) ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      !! 
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij, ip, jp ! 2D addresses
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_tra_rnf')
      !
      igrd = 1                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            ip = bdytmask(ii+1,ij  ) - bdytmask(ii-1,ij  )
            jp = bdytmask(ii  ,ij+1) - bdytmask(ii  ,ij-1)
            tsa(ii,ij,ik,jp_tem) = tsa(ii+ip,ij+jp,ik,jp_tem) * tmask(ii,ij,ik)
            tsa(ii,ij,ik,jp_sal) =                        0.1 * tmask(ii,ij,ik)
         END DO
      END DO
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_tra_rnf')
      !
   END SUBROUTINE bdy_tra_rnf

   SUBROUTINE bdy_tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_dmp  ***
      !!                    
      !! ** Purpose : Apply damping for tracers at open boundaries.
      !! 
      !!----------------------------------------------------------------------
      INTEGER,         INTENT(in) ::   kt
      !! 
      REAL(wp) ::   zwgt           ! boundary weight
      REAL(wp) ::   zta, zsa, ztime
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      INTEGER  ::   ib_bdy         ! Loop index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_tra_dmp')
      !
      DO ib_bdy=1, nb_bdy
         IF ( ln_tra_dmp(ib_bdy) ) THEN
            igrd = 1                       ! Everything is at T-points here
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(ib,igrd)
               DO ik = 1, jpkm1
                  zta = zwgt * ( dta_bdy(ib_bdy)%tem(ib,ik) - tsb(ii,ij,ik,jp_tem) ) * tmask(ii,ij,ik)
                  zsa = zwgt * ( dta_bdy(ib_bdy)%sal(ib,ik) - tsb(ii,ij,ik,jp_sal) ) * tmask(ii,ij,ik)
                  tsa(ii,ij,ik,jp_tem) = tsa(ii,ij,ik,jp_tem) + zta
                  tsa(ii,ij,ik,jp_sal) = tsa(ii,ij,ik,jp_sal) + zsa
               END DO
            END DO
         ENDIF
      ENDDO
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_tra_dmp')
      !
   END SUBROUTINE bdy_tra_dmp
 
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_tra(kt)      ! Empty routine
      WRITE(*,*) 'bdy_tra: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_tra

   SUBROUTINE bdy_tra_dmp(kt)      ! Empty routine
      WRITE(*,*) 'bdy_tra_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_tra_dmp

#endif

   !!======================================================================
END MODULE bdytra
