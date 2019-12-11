MODULE trcini_age
   !!======================================================================
   !!                         ***  MODULE trcini_age  ***
   !! TOP :   initialisation of the AGE tracer
   !!======================================================================
   !! History :   2.0  !  2007-12  (G. Nurser, G. Madec, C. Ethe ) Original code
   !!----------------------------------------------------------------------
#if defined key_age
   !!----------------------------------------------------------------------
   !!   'key_age'                                               AGE tracer
   !!----------------------------------------------------------------------
   !! trc_ini_age   : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trc
   USE trcsms_age

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_age   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_age
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_age  ***  
      !!
      !! ** Purpose :   initialization for AGE model
      !!
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_age: passive tracer age'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      rryear  = 1._wp / ( nyear_len(1) * rday )    ! recip number of seconds in one year

      !! BUG in s-coordinate this does not work!
      nlb_age = MINLOC( gdepw_1d, mask = gdepw_1d > rn_age_depth, dim = 1 ) ! shallowest W level Below age_depth
                                                                            !  = shallowest T level wholly below age_depth

      nl_age  = nlb_age - 1                                                 ! deepest    W level Above age_depth
                                                                            !  = T level surrounding age_depth

      nla_age = nl_age - 1                                                   ! deepest    T level wholly above age_depth

      frac_kill_age = ( rn_age_depth - gdepw_1d(nl_age) ) / e3t_1d(nl_age)      ! fraction of level nl_age above age_depth
      frac_add_age  = 1._wp -  frac_kill_age                                    ! fraction of level nl_age below age_depth

      
      IF( .NOT. ln_rsttr ) trn(:,:,:,jp_age0:jp_age1) = 0.
      !
   END SUBROUTINE trc_ini_age

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No AGE model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_age             ! Empty routine
   END SUBROUTINE trc_ini_age
#endif

   !!======================================================================
END MODULE trcini_age
