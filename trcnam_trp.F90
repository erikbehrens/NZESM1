MODULE trcnam_trp
   !!======================================================================
   !!                       ***  MODULE  trcnam_trp  ***
   !! TOP :   namelist read options for transport
   !!======================================================================
   !! History :   1.0  !  2004-03  (C. Ethe)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_nam_trp  : read the passive tracer namelist for transport
   !!----------------------------------------------------------------------
   USE trc                 ! passive tracers variables
   USE in_out_manager      ! ocean dynamics and active tracers variables
   USE lib_mpp           ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_trp    ! routine called by step module
 
   !                                        !!: ** Advection (namtrc_adv) **
   LOGICAL , PUBLIC ::   ln_trcadv_cen2      ! 2nd order centered scheme flag
   LOGICAL , PUBLIC ::   ln_trcadv_tvd       ! TVD scheme flag
   LOGICAL , PUBLIC ::   ln_trcadv_muscl     ! MUSCL scheme flag
   LOGICAL , PUBLIC ::   ln_trcadv_muscl2    ! MUSCL2 scheme flag
   LOGICAL , PUBLIC ::   ln_trcadv_ubs       ! UBS scheme flag
   LOGICAL , PUBLIC ::   ln_trcadv_qck       ! QUICKEST scheme flag
   LOGICAL , PUBLIC ::   ln_trcadv_msc_ups   ! use upstream scheme within muscl


   !                                        !!: ** lateral mixing namelist (nam_trcldf) **
   LOGICAL , PUBLIC ::   ln_trcldf_lap       !: laplacian operator
   LOGICAL , PUBLIC ::   ln_trcldf_bilap     !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_trcldf_level     !: iso-level direction
   LOGICAL , PUBLIC ::   ln_trcldf_hor       !: horizontal (geopotential) direction
   LOGICAL , PUBLIC ::   ln_trcldf_iso       !: iso-neutral direction
   REAL(wp), PUBLIC ::   rn_ahtrc_0          !: diffusivity coefficient for passive tracer (m2/s)
   REAL(wp), PUBLIC ::   rn_ahtrb_0          !: background diffusivity coefficient for passive tracer (m2/s)
   REAL(wp), PUBLIC ::   rn_fact_lap         !: Enhanced zonal diffusivity coefficent in the equatorial domain

   !                                        !!: ** Treatment of Negative concentrations ( nam_trcrad )
   LOGICAL , PUBLIC ::   ln_trcrad           !: flag to artificially correct negative concentrations

   !                                        !!: ** Vertical diffusion (nam_trczdf) **
   LOGICAL , PUBLIC ::   ln_trczdf_exp       !: explicit vertical diffusion scheme flag
   INTEGER , PUBLIC ::   nn_trczdf_exp       !: number of sub-time step (explicit time stepping)

   !                                                 !!: ** newtonian damping namelist (nam_trcdmp) **
   !                          !!* Namelist namtrc_dmp : passive tracer newtonian damping *
   INTEGER , PUBLIC ::   nn_zdmp_tr    ! = 0/1/2 flag for damping in the mixed layer
   CHARACTER(LEN=200) , PUBLIC :: cn_resto_tr    !File containing restoration coefficient

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_trp
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_trp  ***
      !!                
      !! ** Purpose :   Read Namelist for tracer transport option
      !!----------------------------------------------------------------------
      INTEGER ::  ios                 ! Local integer output status for namelist read
      NAMELIST/namtrc_adv/ ln_trcadv_cen2 , ln_trcadv_tvd   ,    &
         &                 ln_trcadv_muscl, ln_trcadv_muscl2,    &
         &                 ln_trcadv_ubs  , ln_trcadv_qck, ln_trcadv_msc_ups

      NAMELIST/namtrc_ldf/ ln_trcldf_lap  ,     &
         &                 ln_trcldf_bilap, ln_trcldf_level,     &
         &                 ln_trcldf_hor  , ln_trcldf_iso  , rn_ahtrc_0, rn_ahtrb_0,   &
         &                 rn_fact_lap

      NAMELIST/namtrc_zdf/ ln_trczdf_exp  , nn_trczdf_exp
      NAMELIST/namtrc_rad/ ln_trcrad
      NAMELIST/namtrc_dmp/ nn_zdmp_tr , cn_resto_tr
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_trp: read namelist for tracer transport'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~'

      REWIND( numnat_ref )              ! Namelist namtrc_adv in reference namelist : Tracer advection scheme
      READ  ( numnat_ref, namtrc_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_adv in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_adv in configuration namelist : Tracer advection scheme
      READ  ( numnat_cfg, namtrc_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_adv in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_adv )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'trc_adv_ctl : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtrc_adv : chose a advection scheme for tracers'
         WRITE(numout,*) '      2nd order advection scheme     ln_trcadv_cen2   = ', ln_trcadv_cen2
         WRITE(numout,*) '      TVD advection scheme           ln_trcadv_tvd    = ', ln_trcadv_tvd
         WRITE(numout,*) '      MUSCL  advection scheme        ln_trcadv_muscl  = ', ln_trcadv_muscl
         WRITE(numout,*) '      MUSCL2 advection scheme        ln_trcadv_muscl2 = ', ln_trcadv_muscl2
         WRITE(numout,*) '      UBS    advection scheme        ln_trcadv_ubs    = ', ln_trcadv_ubs
         WRITE(numout,*) '      QUICKEST advection scheme      ln_trcadv_qck    = ', ln_trcadv_qck
      ENDIF
      !
      REWIND( numnat_ref )              ! Namelist namtrc_ldf in reference namelist : Tracer lateral diffusive operator
      READ  ( numnat_ref, namtrc_ldf, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_ldf in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_ldf in configuration namelist : Tracer lateral diffusive operator
      READ  ( numnat_cfg, namtrc_ldf, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_ldf in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_ldf )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'trc:ldf_ctl : lateral tracer diffusive operator'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtrc_ldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '      laplacian operator                                 ln_trcldf_lap   = ', ln_trcldf_lap
         WRITE(numout,*) '      bilaplacian operator                               ln_trcldf_bilap = ', ln_trcldf_bilap
         WRITE(numout,*) '      iso-level                                          ln_trcldf_level = ', ln_trcldf_level
         WRITE(numout,*) '      horizontal (geopotential)                          ln_trcldf_hor   = ', ln_trcldf_hor
         WRITE(numout,*) '      iso-neutral                                        ln_trcldf_iso   = ', ln_trcldf_iso
         WRITE(numout,*) '      diffusivity coefficient                                 rn_ahtrc_0 = ', rn_ahtrc_0
         WRITE(numout,*) '      background hor. diffusivity                             rn_ahtrb_0 = ', rn_ahtrb_0
         WRITE(numout,*) '      enhanced zonal diffusivity                             rn_fact_lap = ', rn_fact_lap
      ENDIF

      !                                ! Vertical mixing
      REWIND( numnat_ref )              ! Namelist namtrc_zdf in reference namelist : Tracer vertical mixing
      READ  ( numnat_ref, namtrc_zdf, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_zdf in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_zdf in configuration namelist : Tracer vertical mixing
      READ  ( numnat_cfg, namtrc_zdf, IOSTAT = ios, ERR = 906 )
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_zdf in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_zdf )

      IF(lwp) THEN                     !   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namtrc_zdf : set vertical diffusion parameters'
         WRITE(numout,*) '      time splitting / backward scheme ln_trczdf_exp = ', ln_trczdf_exp
         WRITE(numout,*) '      number of time step              nn_trczdf_exp = ', nn_trczdf_exp
      ENDIF

      !
      REWIND( numnat_ref )              ! Namelist namtrc_rad in reference namelist : Tracer negative concentrations
      READ  ( numnat_ref, namtrc_rad, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_rad in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_rad in configuration namelist : Tracer negative concentrations
      READ  ( numnat_cfg, namtrc_rad, IOSTAT = ios, ERR = 908 )
908   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_rad in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_rad )

      IF(lwp) THEN                     !   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namtrc_rad : treatment of negative concentrations'
         WRITE(numout,*) '      correct artificially negative concen. or not ln_trcrad = ', ln_trcrad
      ENDIF


      REWIND( numnat_ref )              ! Namelist namtrc_dmp in reference namelist : Passive tracers newtonian damping
      READ  ( numnat_ref, namtrc_dmp, IOSTAT = ios, ERR = 909)
909   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_dmp in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_dmp in configuration namelist : Passive tracers newtonian damping
      READ  ( numnat_cfg, namtrc_dmp, IOSTAT = ios, ERR = 910)
910   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_dmp in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_dmp )

      IF(lwp) THEN                       ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'trc_dmp : Passive tracers newtonian damping'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namtrc_dmp : set damping parameter'
         WRITE(numout,*) '      mixed layer damping option     nn_zdmp_tr = ', nn_zdmp_tr, '(zoom: forced to 0)'
         WRITE(numout,*) '      Restoration coeff file    cn_resto_tr = ', cn_resto_tr
      ENDIF
      !
   END SUBROUTINE trc_nam_trp
   
#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                         No TOP model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_trp              ! Empty routine
   END SUBROUTINE trc_nam_trp
#endif

  !!======================================================================
END MODULE trcnam_trp
