MODULE p4zche
   !!======================================================================
   !!                         ***  MODULE p4zche  ***
   !! TOP :   PISCES Sea water chemistry computed following OCMIP protocol
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_che      :  Sea water chemistry computed following OCMIP protocol
   !!----------------------------------------------------------------------
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables
   USE sms_pisces    !  PISCES Source Minus Sink variables
   USE lib_mpp       !  MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_che         !
   PUBLIC   p4z_che_alloc   !

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sio3eq   ! chemistry of Si
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   fekeq    ! chemistry of Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   chemc    ! Solubilities of O2 and CO2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   chemo2   ! Solubilities of O2 and CO2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tempis   ! In situ temperature

   REAL(wp), PUBLIC ::   atcox  = 0.20946         ! units atm

   REAL(wp) ::   salchl = 1. / 1.80655    ! conversion factor for salinity --> chlorinity (Wooster et al. 1969)
   REAL(wp) ::   o2atm  = 1. / ( 1000. * 0.20946 )  

   REAL(wp) ::   rgas   = 83.14472       ! universal gas constants
   REAL(wp) ::   oxyco  = 1. / 22.4144   ! converts from liters of an ideal gas to moles

   REAL(wp) ::   bor1   = 0.00023        ! borat constants
   REAL(wp) ::   bor2   = 1. / 10.82

   REAL(wp) ::   st1    =      0.14     ! constants for calculate concentrations for sulfate
   REAL(wp) ::   st2    =  1./96.062    !  (Morris & Riley 1966)

   REAL(wp) ::   ft1    =    0.000067   ! constants for calculate concentrations for fluorides
   REAL(wp) ::   ft2    = 1./18.9984    ! (Dickson & Riley 1979 )

   !                                    ! volumetric solubility constants for o2 in ml/L  
   REAL(wp) ::   ox0    =  2.00856      ! from Table 1 for Eq 8 of Garcia and Gordon, 1992.
   REAL(wp) ::   ox1    =  3.22400      ! corrects for moisture and fugacity, but not total atmospheric pressure
   REAL(wp) ::   ox2    =  3.99063      !      Original PISCES code noted this was a solubility, but 
   REAL(wp) ::   ox3    =  4.80299      ! was in fact a bunsen coefficient with units L-O2/(Lsw atm-O2)
   REAL(wp) ::   ox4    =  9.78188e-1   ! Hence, need to divide EXP( zoxy ) by 1000, ml-O2 => L-O2
   REAL(wp) ::   ox5    =  1.71069      ! and atcox = 0.20946 to add the 1/atm dimension.
   REAL(wp) ::   ox6    = -6.24097e-3   
   REAL(wp) ::   ox7    = -6.93498e-3 
   REAL(wp) ::   ox8    = -6.90358e-3
   REAL(wp) ::   ox9    = -4.29155e-3 
   REAL(wp) ::   ox10   = -3.11680e-7 

   !                                    ! coeff. for seawater pressure correction : millero 95
   !                                    ! AGRIF doesn't like the DATA instruction
   REAL(wp) :: devk11  = -25.5
   REAL(wp) :: devk12  = -15.82
   REAL(wp) :: devk13  = -29.48
   REAL(wp) :: devk14  = -25.60
   REAL(wp) :: devk15  = -48.76
   !
   REAL(wp) :: devk21  = 0.1271
   REAL(wp) :: devk22  = -0.0219
   REAL(wp) :: devk23  = 0.1622
   REAL(wp) :: devk24  = 0.2324
   REAL(wp) :: devk25  = 0.5304
   !
   REAL(wp) :: devk31  = 0.
   REAL(wp) :: devk32  = 0.
   REAL(wp) :: devk33  = 2.608E-3
   REAL(wp) :: devk34  = -3.6246E-3
   REAL(wp) :: devk35  = 0.
   !
   REAL(wp) :: devk41  = -3.08E-3
   REAL(wp) :: devk42  = 1.13E-3
   REAL(wp) :: devk43  = -2.84E-3
   REAL(wp) :: devk44  = -5.13E-3
   REAL(wp) :: devk45  = -11.76E-3
   !
   REAL(wp) :: devk51  = 0.0877E-3
   REAL(wp) :: devk52  = -0.1475E-3     
   REAL(wp) :: devk53  = 0.
   REAL(wp) :: devk54  = 0.0794E-3      
   REAL(wp) :: devk55  = 0.3692E-3      

   !!* Substitution
#include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_che
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_che  ***
      !!
      !! ** Purpose :   Sea water chemistry computed following OCMIP protocol
      !!
      !! ** Method  : - ...
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   ztkel, zt   , zt2   , zsal  , zsal2 , zbuf1 , zbuf2
      REAL(wp) ::   ztgg , ztgg2, ztgg3 , ztgg4 , ztgg5
      REAL(wp) ::   zpres, ztc  , zcl   , zcpexp, zoxy  , zcpexp2
      REAL(wp) ::   zsqrt, ztr  , zlogt , zcek1, zc1, zplat
      REAL(wp) ::   zis  , zis2 , zsal15, zisqrt, za1  , za2
      REAL(wp) ::   zckb , zck1 , zck2  , zckw  , zak1 , zak2  , zakb , zaksp0, zakw
      REAL(wp) ::   zst  , zft  , zcks  , zckf  , zaksp1
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_che')
      !
      ! Computations of chemical constants require in situ temperature
      ! Here a quite simple formulation is used to convert 
      ! potential temperature to in situ temperature. The errors is less than 
      ! 0.04°C relative to an exact computation
      ! ---------------------------------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zpres = fsdept(ji,jj,jk) / 1000.
               za1 = 0.04 * ( 1.0 + 0.185 * tsn(ji,jj,jk,jp_tem) + 0.035 * (tsn(ji,jj,jk,jp_sal) - 35.0) )
               za2 = 0.0075 * ( 1.0 - tsn(ji,jj,jk,jp_tem) / 30.0 )
               tempis(ji,jj,jk) = tsn(ji,jj,jk,jp_tem) - za1 * zpres + za2 * zpres**2
            END DO
         END DO
      END DO
      !
      ! CHEMICAL CONSTANTS - SURFACE LAYER
      ! ----------------------------------
!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            !                             ! SET ABSOLUTE TEMPERATURE
            ztkel = tempis(ji,jj,1) + 273.15
            zt    = ztkel * 0.01
            zt2   = zt * zt
            zsal  = tsn(ji,jj,1,jp_sal) + ( 1.- tmask(ji,jj,1) ) * 35.
            zsal2 = zsal * zsal
            zlogt = LOG( zt )
            !                             ! LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1980)
            !                             !     AND FOR THE ATMOSPHERE FOR NON IDEAL GAS
            zcek1 = 9345.17/ztkel - 60.2409 + 23.3585 * LOG(zt) + zsal*(0.023517 - 0.00023656*ztkel    &
            &       + 0.0047036e-4*ztkel**2)
            !                             ! SET SOLUBILITIES OF O2 AND CO2 
            chemc(ji,jj,1) = EXP( zcek1 ) * 1.e-6 * rhop(ji,jj,1) / 1000. ! mol/(kg uatm)
            chemc(ji,jj,2) = -1636.75 + 12.0408*ztkel - 0.0327957*ztkel**2 + 0.0000316528*ztkel**3
            chemc(ji,jj,3) = 57.7 - 0.118*ztkel
            !
         END DO
      END DO

      ! OXYGEN SOLUBILITY - DEEP OCEAN
      ! -------------------------------
!CDIR NOVERRCHK
      DO jk = 1, jpk
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
              ztkel = tempis(ji,jj,jk) + 273.15
              zsal  = tsn(ji,jj,jk,jp_sal) + ( 1.- tmask(ji,jj,jk) ) * 35.
              zsal2 = zsal * zsal
              ztgg  = LOG( ( 298.15 - tempis(ji,jj,jk) ) / ztkel )  ! Set the GORDON & GARCIA scaled temperature
              ztgg2 = ztgg  * ztgg
              ztgg3 = ztgg2 * ztgg
              ztgg4 = ztgg3 * ztgg
              ztgg5 = ztgg4 * ztgg
              zoxy  = ox0 + ox1 * ztgg + ox2 * ztgg2 + ox3 * ztgg3 + ox4 * ztgg4 + ox5 * ztgg5   &
                     + zsal * ( ox6 + ox7 * ztgg + ox8 * ztgg2 + ox9 * ztgg3 ) +  ox10 * zsal2
              chemo2(ji,jj,jk) = ( EXP( zoxy ) * o2atm ) * oxyco * atcox     ! mol/(L atm)
            END DO
          END DO
        END DO



      ! CHEMICAL CONSTANTS - DEEP OCEAN
      ! -------------------------------
!CDIR NOVERRCHK
      DO jk = 1, jpk
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi

               ! SET PRESSION ACCORDING TO SAUNDER (1980)
               zplat   = SIN ( ABS(gphit(ji,jj)*3.141592654/180.) )
               zc1 = 5.92E-3 + zplat**2 * 5.25E-3
               zpres = ((1-zc1)-SQRT(((1-zc1)**2)-(8.84E-6*fsdept(ji,jj,jk)))) / 4.42E-6
               zpres = zpres / 10.0

               ! SET ABSOLUTE TEMPERATURE
               ztkel   = tempis(ji,jj,jk) + 273.15
               zsal    = tsn(ji,jj,jk,jp_sal) + ( 1.-tmask(ji,jj,jk) ) * 35.
               zsqrt  = SQRT( zsal )
               zsal15  = zsqrt * zsal
               zlogt  = LOG( ztkel )
               ztr    = 1. / ztkel
               zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
               zis2   = zis * zis
               zisqrt = SQRT( zis )
               ztc     = tempis(ji,jj,jk) + ( 1.- tmask(ji,jj,jk) ) * 20.

               ! CHLORINITY (WOOSTER ET AL., 1969)
               zcl     = zsal * salchl

               ! TOTAL SULFATE CONCENTR. [MOLES/kg soln]
               zst     = st1 * zcl * st2

               ! TOTAL FLUORIDE CONCENTR. [MOLES/kg soln]
               zft     = ft1 * zcl * ft2

               ! DISSOCIATION CONSTANT FOR SULFATES on free H scale (Dickson 1990)
               zcks    = EXP(-4276.1 * ztr + 141.328 - 23.093 * zlogt         &
               &         + (-13856. * ztr + 324.57 - 47.986 * zlogt) * zisqrt &
               &         + (35474. * ztr - 771.54 + 114.723 * zlogt) * zis    &
               &         - 2698. * ztr * zis**1.5 + 1776.* ztr * zis2         &
               &         + LOG(1.0 - 0.001005 * zsal))
               !
               aphscale(ji,jj,jk) = ( 1. + zst / zcks )

               ! DISSOCIATION CONSTANT FOR FLUORIDES on free H scale (Dickson and Riley 79)
               zckf    = EXP( 1590.2*ztr - 12.641 + 1.525*zisqrt   &
               &         + LOG(1.0d0 - 0.001005d0*zsal)            &
               &         + LOG(1.0d0 + zst/zcks))

               ! DISSOCIATION CONSTANT FOR CARBONATE AND BORATE
               zckb=  (-8966.90 - 2890.53*zsqrt - 77.942*zsal        &
               &      + 1.728*zsal15 - 0.0996*zsal*zsal)*ztr         &
               &      + (148.0248 + 137.1942*zsqrt + 1.62142*zsal)   &
               &      + (-24.4344 - 25.085*zsqrt - 0.2474*zsal)      & 
               &      * zlogt + 0.053105*zsqrt*ztkel


               ! DISSOCIATION COEFFICIENT FOR CARBONATE ACCORDING TO 
               ! MEHRBACH (1973) REFIT BY MILLERO (1995), seawater scale
               zck1    = -1.0*(3633.86*ztr - 61.2172 + 9.6777*zlogt  &
                  - 0.011555*zsal + 0.0001152*zsal*zsal)
               zck2    = -1.0*(471.78*ztr + 25.9290 - 3.16967*zlogt      &
                  - 0.01781*zsal + 0.0001122*zsal*zsal)

               ! PKW (H2O) (DICKSON AND RILEY, 1979)
               zckw = -13847.26*ztr + 148.9652 - 23.6521 * zlogt    & 
               &     + (118.67*ztr - 5.977 + 1.0495 * zlogt)        &
               &     * zsqrt - 0.01615 * zsal

               ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
               !       (S=27-43, T=2-25 DEG C) at pres =0 (atmos. pressure) (MUCCI 1983)
               zaksp0  = -171.9065 -0.077993*ztkel + 2839.319*ztr + 71.595*LOG10( ztkel )   &
                  &      + (-0.77712 + 0.00284263*ztkel + 178.34*ztr) * zsqrt  &
                  &      - 0.07711*zsal + 0.0041249*zsal15

               ! K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
               zak1    = 10**(zck1)
               zak2    = 10**(zck2)
               zakb    = EXP( zckb  )
               zakw    = EXP( zckw )
               zaksp1  = 10**(zaksp0)

               ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
               !        (REFERENCE TO CULBERSON & PYTKOQICZ (1968) AS MADE
               !        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
               !        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
               !        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
               !        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
               !        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
               !        & GIESKES (1970), P. 1285-1286 (THE SMALL
               !        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
               !        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON P. 1285))
               zcpexp  = zpres /(rgas*ztkel)
               zcpexp2 = zpres * zpres/(rgas*ztkel)

               ! KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
               !        CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
               !        (CF. BROECKER ET AL., 1982)

               zbuf1  = -     ( devk11 + devk21 * ztc + devk31 * ztc * ztc )
               zbuf2  = 0.5 * ( devk41 + devk51 * ztc )
               ak13(ji,jj,jk) = zak1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk12 + devk22 * ztc + devk32 * ztc * ztc )
               zbuf2  = 0.5 * ( devk42 + devk52 * ztc )
               ak23(ji,jj,jk) = zak2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk13 + devk23 * ztc + devk33 * ztc * ztc )
               zbuf2  = 0.5 * ( devk43 + devk53 * ztc )
               akb3(ji,jj,jk) = zakb * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk14 + devk24 * ztc + devk34 * ztc * ztc )
               zbuf2  = 0.5 * ( devk44 + devk54 * ztc )
               akw3(ji,jj,jk) = zakw * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )


               ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE 
               !        AS FUNCTION OF PRESSURE FOLLOWING MILLERO
               !        (P. 1285) AND BERNER (1976)
               zbuf1  =     - ( devk15 + devk25 * ztc + devk35 * ztc * ztc )
               zbuf2  = 0.5 * ( devk45 + devk55 * ztc )
               aksp(ji,jj,jk) = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               ! TOTAL BORATE CONCENTR. [MOLES/L]
               borat(ji,jj,jk) = bor1 * zcl * bor2

               ! Iron and SIO3 saturation concentration from ...
               sio3eq(ji,jj,jk) = EXP(  LOG( 10.) * ( 6.44 - 968. / ztkel )  ) * 1.e-6
               fekeq (ji,jj,jk) = 10**( 17.27 - 1565.7 / ( 273.15 + ztc ) )

            END DO
         END DO
      END DO
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_che')
      !
   END SUBROUTINE p4z_che


   INTEGER FUNCTION p4z_che_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_che_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( sio3eq(jpi,jpj,jpk), fekeq(jpi,jpj,jpk), chemc(jpi,jpj,3), chemo2(jpi,jpj,jpk),   &
      &         tempis(jpi,jpj,jpk), STAT=p4z_che_alloc )
      !
      IF( p4z_che_alloc /= 0 )   CALL ctl_warn('p4z_che_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_che_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_che( kt )                   ! Empty routine
      INTEGER, INTENT(in) ::   kt
      WRITE(*,*) 'p4z_che: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_che
#endif 

   !!======================================================================
END MODULE p4zche
