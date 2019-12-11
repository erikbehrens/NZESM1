MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    PISCES :   Output of PISCES tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput && ( defined key_pisces || defined key_pisces_reduced )
   !!----------------------------------------------------------------------
   !!   'key_pisces or key_pisces_reduced'                     PISCES model
   !!----------------------------------------------------------------------
   !! trc_wri_pisces   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE sms_pisces  ! PISCES variables
   USE iom         ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_pisces 

#  include "top_substitute.h90"
CONTAINS

   SUBROUTINE trc_wri_pisces
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)           :: cltra
      REAL(wp)                     :: zfact
      INTEGER                      :: ji, jj, jk, jn
      REAL(wp), DIMENSION(jpi,jpj) :: zdic, zo2min, zdepo2min
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
#if defined key_pisces_reduced
      DO jn = jp_pcs0, jp_pcs1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, trn(:,:,:,jn) )
      END DO
#else
      DO jn = jp_pcs0, jp_pcs1
         zfact = 1.0e+6 
         IF( jn == jpno3 .OR. jn == jpnh4 ) zfact = rno3 * 1.0e+6 
         IF( jn == jppo4  )                 zfact = po4r * 1.0e+6
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         IF( iom_use( cltra ) )  CALL iom_put( cltra, trn(:,:,:,jn) * zfact )
      END DO

      IF( iom_use( "INTDIC" ) ) THEN                     !   DIC content in kg/m2
         zdic(:,:) = 0.
         DO jk = 1, jpkm1
            zdic(:,:) = zdic(:,:) + trn(:,:,jk,jpdic) * fse3t(:,:,jk) * tmask(:,:,jk) * 12.
         ENDDO
         CALL iom_put( 'INTDIC', zdic )     
      ENDIF
      !
      IF( iom_use( "O2MIN" ) .OR. iom_use ( "ZO2MIN" ) ) THEN  ! Oxygen minimum concentration and depth 
         zo2min   (:,:) = trn(:,:,1,jpoxy) * tmask(:,:,1)
         zdepo2min(:,:) = fsdepw(:,:,1)    * tmask(:,:,1)
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( tmask(ji,jj,jk) == 1 ) then
                     IF( trn(ji,jj,jk,jpoxy) < zo2min(ji,jj) ) then
                        zo2min   (ji,jj) = trn(ji,jj,jk,jpoxy)
                        zdepo2min(ji,jj) = fsdepw(ji,jj,jk)
                     ENDIF
                  ENDIF
               END DO
            END DO
         END DO
         !
         CALL iom_put('O2MIN' , zo2min     )                              ! oxygen minimum concentration
         CALL iom_put('ZO2MIN', zdepo2min  )                              ! depth of oxygen minimum concentration
          !
      ENDIF
#endif
      !
   END SUBROUTINE trc_wri_pisces

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_pisces
CONTAINS
   SUBROUTINE trc_wri_pisces                     ! Empty routine  
   END SUBROUTINE trc_wri_pisces
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_pisces
