MODULE albedo
   !!======================================================================
   !!                       ***  MODULE  albedo  ***
   !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
   !!=====================================================================
   !! History :  8.0  ! 2001-04  (LIM 1.0)
   !!   NEMO     1.0  ! 2003-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
   !!             -   ! 2004-11  (C. Talandier)  add albedo_init
   !!             -   ! 2001-06  (M. Vancoppenolle) LIM 3.0
   !!             -   ! 2006-08  (G. Madec)  cleaning for surface module
   !!            3.6  ! 2016-01  (C. Rousset) new parameterization for sea ice albedo
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   albedo_ice    : albedo for   ice (clear and overcast skies)
   !!   albedo_oce    : albedo for ocean (clear and overcast skies)
   !!   albedo_init   : initialisation of albedo computation
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   albedo_ice   ! routine called sbcice_lim.F90
   PUBLIC   albedo_oce   ! routine called by ???

   INTEGER  ::   albd_init = 0      !: control flag for initialization
  
   REAL(wp) ::   rmue     = 0.40    !  cosine of local solar altitude
   REAL(wp) ::   ralb_oce = 0.066   ! ocean or lead albedo (Pegau and Paulson, Ann. Glac. 2001)
   REAL(wp) ::   c1       = 0.05    ! snow thickness (only for nn_ice_alb=0)
   REAL(wp) ::   c2       = 0.10    !  "        "
   REAL(wp) ::   rcloud   = 0.06    ! cloud effect on albedo (only-for nn_ice_alb=0)
 
   !                             !!* namelist namsbc_alb
   INTEGER  ::   nn_ice_alb
   REAL(wp) ::   rn_albice

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE albedo_ice( pt_ice, ph_ice, ph_snw, pa_ice_cs, pa_ice_os )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE albedo_ice  ***
      !!          
      !! ** Purpose :   Computation of the albedo of the snow/ice system 
      !!       
      !! ** Method  :   Two schemes are available (from namelist parameter nn_ice_alb)
      !!                  0: the scheme is that of Shine & Henderson-Sellers (JGR 1985) for clear-skies
      !!                  1: the scheme is "home made" (for cloudy skies) and based on Brandt et al. (J. Climate 2005)
      !!                                                                           and Grenfell & Perovich (JGR 2004)
      !!                Description of scheme 1:
      !!                  1) Albedo dependency on ice thickness follows the findings from Brandt et al (2005)
      !!                     which are an update of Allison et al. (JGR 1993) ; Brandt et al. 1999
      !!                     0-5cm  : linear function of ice thickness
      !!                     5-150cm: log    function of ice thickness
      !!                     > 150cm: constant
      !!                  2) Albedo dependency on snow thickness follows the findings from Grenfell & Perovich (2004)
      !!                     i.e. it increases as -EXP(-snw_thick/0.02) during freezing and -EXP(-snw_thick/0.03) during melting
      !!                  3) Albedo dependency on clouds is speculated from measurements of Grenfell and Perovich (2004)
      !!                     i.e. cloudy-clear albedo depend on cloudy albedo following a 2d order polynomial law
      !!                  4) The needed 4 parameters are: dry and melting snow, freezing ice and bare puddled ice
      !!
      !! ** Note    :   The parameterization from Shine & Henderson-Sellers presents several misconstructions:
      !!                  1) ice albedo when ice thick. tends to 0 is different than ocean albedo
      !!                  2) for small ice thick. covered with some snow (<3cm?), albedo is larger 
      !!                     under melting conditions than under freezing conditions
      !!                  3) the evolution of ice albedo as a function of ice thickness shows  
      !!                     3 sharp inflexion points (at 5cm, 100cm and 150cm) that look highly unrealistic
      !!
      !! References :   Shine & Henderson-Sellers 1985, JGR, 90(D1), 2243-2250.
      !!                Brandt et al. 2005, J. Climate, vol 18
      !!                Grenfell & Perovich 2004, JGR, vol 109 
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:) ::   pt_ice      !  ice surface temperature (Kelvin)
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:) ::   ph_ice      !  sea-ice thickness
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:) ::   ph_snw      !  snow thickness
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pa_ice_cs   !  albedo of ice under clear    sky
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pa_ice_os   !  albedo of ice under overcast sky
      !!
      INTEGER  ::   ji, jj, jl         ! dummy loop indices
      INTEGER  ::   ijpl               ! number of ice categories (3rd dim of ice input arrays)
      REAL(wp)            ::   ralb_im, ralb_sf, ralb_sm, ralb_if
      REAL(wp)            ::   zswitch, z1_c1, z1_c2
      REAL(wp)                            ::   zalb_sm, zalb_sf, zalb_st ! albedo of snow melting, freezing, total
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zalb, zalb_it             ! intermediate variable & albedo of ice (snow free)
      !!---------------------------------------------------------------------

      ijpl = SIZE( pt_ice, 3 )                     ! number of ice categories
      
      CALL wrk_alloc( jpi,jpj,ijpl, zalb, zalb_it )

      IF( albd_init == 0 )   CALL albedo_init      ! initialization 

      
      SELECT CASE ( nn_ice_alb )

      !------------------------------------------
      !  Shine and Henderson-Sellers (1985)
      !------------------------------------------
      CASE( 0 )
       
         ralb_sf = 0.80       ! dry snow
         ralb_sm = 0.65       ! melting snow
         ralb_if = 0.72       ! bare frozen ice
         ralb_im = rn_albice  ! bare puddled ice 
         
         !  Computation of ice albedo (free of snow)
         WHERE     ( ph_snw == 0._wp .AND. pt_ice >= rt0_ice )   ;   zalb(:,:,:) = ralb_im
         ELSE WHERE                                              ;   zalb(:,:,:) = ralb_if
         END  WHERE
      
         WHERE     ( 1.5  < ph_ice                     )  ;  zalb_it = zalb
         ELSE WHERE( 1.0  < ph_ice .AND. ph_ice <= 1.5 )  ;  zalb_it = 0.472  + 2.0 * ( zalb - 0.472 ) * ( ph_ice - 1.0 )
         ELSE WHERE( 0.05 < ph_ice .AND. ph_ice <= 1.0 )  ;  zalb_it = 0.2467 + 0.7049 * ph_ice              &
            &                                                                 - 0.8608 * ph_ice * ph_ice     &
            &                                                                 + 0.3812 * ph_ice * ph_ice * ph_ice
         ELSE WHERE                                       ;  zalb_it = 0.1    + 3.6    * ph_ice
         END WHERE
     
         DO jl = 1, ijpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! freezing snow
                  ! no effect of underlying ice layer IF snow thickness > c1. Albedo does not depend on snow thick if > c2
                  !                                        !  freezing snow        
                  zswitch   = 1._wp - MAX( 0._wp , SIGN( 1._wp , - ( ph_snw(ji,jj,jl) - c1 ) ) )
                  zalb_sf   = ( 1._wp - zswitch ) * (  zalb_it(ji,jj,jl)  &
                     &                           + ph_snw(ji,jj,jl) * ( ralb_sf - zalb_it(ji,jj,jl) ) / c1  )   &
                     &        +         zswitch   * ralb_sf  

                  ! melting snow
                  ! no effect of underlying ice layer. Albedo does not depend on snow thick IF > c2
                  zswitch   = MAX( 0._wp , SIGN( 1._wp , ph_snw(ji,jj,jl) - c2 ) )
                  zalb_sm = ( 1._wp - zswitch ) * ( ralb_im + ph_snw(ji,jj,jl) * ( ralb_sm - ralb_im ) / c2 )   &
                      &     +         zswitch   *   ralb_sm 
                  !
                  ! snow albedo
                  zswitch  =  MAX( 0._wp , SIGN( 1._wp , pt_ice(ji,jj,jl) - rt0_snow ) )   
                  zalb_st  =  zswitch * zalb_sm + ( 1._wp - zswitch ) * zalb_sf
               
                  ! Ice/snow albedo
                  zswitch   = 1._wp - MAX( 0._wp , SIGN( 1._wp , - ph_snw(ji,jj,jl) ) )
                  pa_ice_cs(ji,jj,jl) =  zswitch * zalb_st + ( 1._wp - zswitch ) * zalb_it(ji,jj,jl)
                  !
               END DO
            END DO
         END DO

         pa_ice_os(:,:,:) = pa_ice_cs(:,:,:) + rcloud       ! Oberhuber correction for overcast sky

      !------------------------------------------
      !  New parameterization (2016)
      !------------------------------------------
      CASE( 1 ) 

         ralb_im = rn_albice  ! bare puddled ice
! compilation of values from literature
         ralb_sf = 0.85      ! dry snow
         ralb_sm = 0.75      ! melting snow
         ralb_if = 0.60      ! bare frozen ice
! Perovich et al 2002 (Sheba) => the only dataset for which all types of ice/snow were retrieved
!         ralb_sf = 0.85       ! dry snow
!         ralb_sm = 0.72       ! melting snow
!         ralb_if = 0.65       ! bare frozen ice
! Brandt et al 2005 (East Antarctica)
!         ralb_sf = 0.87      ! dry snow
!         ralb_sm = 0.82      ! melting snow
!         ralb_if = 0.54      ! bare frozen ice
! 
         !  Computation of ice albedo (free of snow)
         z1_c1 = 1. / ( LOG(1.5) - LOG(0.05) ) 
         z1_c2 = 1. / 0.05
         WHERE     ( ph_snw == 0._wp .AND. pt_ice >= rt0_ice )   ;   zalb = ralb_im
         ELSE WHERE                                              ;   zalb = ralb_if
         END  WHERE
         
         WHERE     ( 1.5  < ph_ice                     )  ;  zalb_it = zalb
         ELSE WHERE( 0.05 < ph_ice .AND. ph_ice <= 1.5 )  ;  zalb_it = zalb     + ( 0.18 - zalb     ) * z1_c1 *  &
            &                                                                     ( LOG(1.5) - LOG(ph_ice) )
         ELSE WHERE                                       ;  zalb_it = ralb_oce + ( 0.18 - ralb_oce ) * z1_c2 * ph_ice
         END WHERE

         z1_c1 = 1. / 0.02
         z1_c2 = 1. / 0.03
         !  Computation of the snow/ice albedo
         DO jl = 1, ijpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zalb_sf = ralb_sf - ( ralb_sf - zalb_it(ji,jj,jl)) * EXP( - ph_snw(ji,jj,jl) * z1_c1 );
                  zalb_sm = ralb_sm - ( ralb_sm - zalb_it(ji,jj,jl)) * EXP( - ph_snw(ji,jj,jl) * z1_c2 );

                   ! snow albedo
                  zswitch = MAX( 0._wp , SIGN( 1._wp , pt_ice(ji,jj,jl) - rt0_snow ) )   
                  zalb_st = zswitch * zalb_sm + ( 1._wp - zswitch ) * zalb_sf

                  ! Ice/snow albedo   
                  zswitch             = MAX( 0._wp , SIGN( 1._wp , - ph_snw(ji,jj,jl) ) )
                  pa_ice_os(ji,jj,jl) = ( 1._wp - zswitch ) * zalb_st + zswitch *  zalb_it(ji,jj,jl)

              END DO
            END DO
         END DO
         ! Effect of the clouds (2d order polynomial)
         pa_ice_cs = pa_ice_os - ( - 0.1010 * pa_ice_os * pa_ice_os + 0.1933 * pa_ice_os - 0.0148 ); 

      END SELECT
      
      CALL wrk_dealloc( jpi,jpj,ijpl, zalb, zalb_it )
      !
   END SUBROUTINE albedo_ice


   SUBROUTINE albedo_oce( pa_oce_os , pa_oce_cs )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE albedo_oce  ***
      !! 
      !! ** Purpose :   Computation of the albedo of the ocean
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pa_oce_os   !  albedo of ocean under overcast sky
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pa_oce_cs   !  albedo of ocean under clear sky
      !!
      REAL(wp) :: zcoef 
      !!----------------------------------------------------------------------
      !
      zcoef = 0.05 / ( 1.1 * rmue**1.4 + 0.15 )   ! Parameterization of Briegled and Ramanathan, 1982
      pa_oce_cs(:,:) = zcoef 
      pa_oce_os(:,:) = 0.06                       ! Parameterization of Kondratyev, 1969 and Payne, 1972
      !
   END SUBROUTINE albedo_oce


   SUBROUTINE albedo_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE albedo_init  ***
      !!
      !! ** Purpose :   initializations for the albedo parameters
      !!
      !! ** Method  :   Read the namelist namsbc_alb
      !!----------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      NAMELIST/namsbc_alb/ nn_ice_alb, rn_albice 
      !!----------------------------------------------------------------------
      !
      albd_init = 1                     ! indicate that the initialization has been done
      !
      REWIND( numnam_ref )              ! Namelist namsbc_alb in reference namelist : Albedo parameters
      READ  ( numnam_ref, namsbc_alb, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_alb in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsbc_alb in configuration namelist : Albedo parameters
      READ  ( numnam_cfg, namsbc_alb, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_alb in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsbc_alb )
      !
      IF(lwp) THEN                      ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'albedo : set albedo parameters'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namsbc_alb : albedo '
         WRITE(numout,*) '      choose the albedo parameterization                  nn_ice_alb = ', nn_ice_alb
         WRITE(numout,*) '      albedo of bare puddled ice                          rn_albice  = ', rn_albice
      ENDIF
      !
   END SUBROUTINE albedo_init

   !!======================================================================
END MODULE albedo
