#define SPONGE && define SPONGE_TOP

MODULE agrif_opa_sponge
#if defined key_agrif  && ! defined key_offline
   USE par_oce
   USE oce
   USE dom_oce
   USE in_out_manager
   USE agrif_oce
   USE wrk_nemo  
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Sponge, Agrif_Sponge_Tra, Agrif_Sponge_Dyn
   PUBLIC interptsn_sponge, interpun_sponge, interpvn_sponge

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3 , NEMO Consortium (2010)
   !! $Id: agrif_opa_sponge.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE Agrif_Sponge_Tra
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Sponge_Tra ***
      !!---------------------------------------------
      !!
      REAL(wp) :: timecoeff

#if defined SPONGE
      timecoeff = REAL(Agrif_NbStepint(),wp)/Agrif_rhot()

      CALL Agrif_Sponge
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      tabspongedone_tsn = .FALSE.

      CALL Agrif_Bc_Variable(tsn_sponge_id,calledweight=timecoeff,procname=interptsn_sponge)

      Agrif_UseSpecialValue = .FALSE.
#endif

   END SUBROUTINE Agrif_Sponge_Tra

   SUBROUTINE Agrif_Sponge_dyn
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Sponge_dyn ***
      !!---------------------------------------------
      !!
      REAL(wp) :: timecoeff

#if defined SPONGE
      timecoeff = REAL(Agrif_NbStepint(),wp)/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn

      tabspongedone_u = .FALSE.
      tabspongedone_v = .FALSE.         
      CALL Agrif_Bc_Variable(un_sponge_id,calledweight=timecoeff,procname=interpun_sponge)

      tabspongedone_u = .FALSE.
      tabspongedone_v = .FALSE.
      CALL Agrif_Bc_Variable(vn_sponge_id,calledweight=timecoeff,procname=interpvn_sponge)

      Agrif_UseSpecialValue = .FALSE.
#endif

   END SUBROUTINE Agrif_Sponge_dyn

   SUBROUTINE Agrif_Sponge
      !!---------------------------------------------
      !!   *** ROUTINE  Agrif_Sponge ***
      !!---------------------------------------------
      INTEGER  :: ji,jj,jk
      INTEGER  :: ispongearea, ilci, ilcj
      LOGICAL  :: ll_spdone
      REAL(wp) :: z1spongearea, zramp
      REAL(wp), POINTER, DIMENSION(:,:) :: ztabramp

#if defined SPONGE || defined SPONGE_TOP
      ll_spdone=.TRUE.
      IF (( .NOT. spongedoneT ).OR.( .NOT. spongedoneU )) THEN
         ! Define ramp from boundaries towards domain interior
         ! at T-points
         ! Store it in ztabramp
         ll_spdone=.FALSE.

         CALL wrk_alloc( jpi, jpj, ztabramp )

         ispongearea  = 2 + nn_sponge_len * Agrif_irhox()
         ilci = nlci - ispongearea
         ilcj = nlcj - ispongearea 
         z1spongearea = 1._wp / REAL( ispongearea - 2 )

         ztabramp(:,:) = 0._wp

         IF( (nbondi == -1) .OR. (nbondi == 2) ) THEN
            DO jj = 1, jpj
               IF ( umask(2,jj,1) == 1._wp ) THEN
                 DO ji = 2, ispongearea                  
                    ztabramp(ji,jj) = ( ispongearea-ji ) * z1spongearea
                 END DO
               ENDIF
            ENDDO
         ENDIF

         IF( (nbondi == 1) .OR. (nbondi == 2) ) THEN
            DO jj = 1, jpj
               IF ( umask(nlci-2,jj,1) == 1._wp ) THEN
                  DO ji = ilci+1,nlci-1
                     zramp = (ji - (ilci+1) ) * z1spongearea
                     ztabramp(ji,jj) = MAX( ztabramp(ji,jj), zramp )
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

         IF( (nbondj == -1) .OR. (nbondj == 2) ) THEN
            DO ji = 1, jpi
               IF ( vmask(ji,2,1) == 1._wp ) THEN
                  DO jj = 2, ispongearea
                     zramp = ( ispongearea-jj ) * z1spongearea
                     ztabramp(ji,jj) = MAX( ztabramp(ji,jj), zramp )
                  END DO
               ENDIF
            ENDDO
         ENDIF

         IF( (nbondj == 1) .OR. (nbondj == 2) ) THEN
            DO ji = 1, jpi
               IF ( vmask(ji,nlcj-2,1) == 1._wp ) THEN
                  DO jj = ilcj+1,nlcj-1
                     zramp = (jj - (ilcj+1) ) * z1spongearea
                     ztabramp(ji,jj) = MAX( ztabramp(ji,jj), zramp )
                  END DO
               ENDIF
            ENDDO
         ENDIF

      ENDIF

      ! Tracers
      IF( .NOT. spongedoneT ) THEN
         fsaht_spu(:,:) = 0._wp
         fsaht_spv(:,:) = 0._wp
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               fsaht_spu(ji,jj) = 0.5_wp * visc_tra * (ztabramp(ji,jj) + ztabramp(ji+1,jj  ))
               fsaht_spv(ji,jj) = 0.5_wp * visc_tra * (ztabramp(ji,jj) + ztabramp(ji  ,jj+1))
            END DO
         END DO

         CALL lbc_lnk( fsaht_spu, 'U', 1. )   ! Lateral boundary conditions
         CALL lbc_lnk( fsaht_spv, 'V', 1. )
         spongedoneT = .TRUE.
      ENDIF

      ! Dynamics
      IF( .NOT. spongedoneU ) THEN
         fsahm_spt(:,:) = 0._wp
         fsahm_spf(:,:) = 0._wp
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               fsahm_spt(ji,jj) = visc_dyn * ztabramp(ji,jj)
               fsahm_spf(ji,jj) = 0.25_wp * visc_dyn * ( ztabramp(ji,jj) + ztabramp(ji  ,jj+1) &
                                                     &  +ztabramp(ji,jj) + ztabramp(ji+1,jj  ) )
            END DO
         END DO

         CALL lbc_lnk( fsahm_spt, 'T', 1. )   ! Lateral boundary conditions
         CALL lbc_lnk( fsahm_spf, 'F', 1. )
         spongedoneU = .TRUE.
      ENDIF
      !
      IF (.NOT.ll_spdone) CALL wrk_dealloc( jpi, jpj, ztabramp )
      !
#endif

   END SUBROUTINE Agrif_Sponge

   SUBROUTINE interptsn_sponge(tabres,i1,i2,j1,j2,k1,k2,n1,n2,before)
      !!---------------------------------------------
      !!   *** ROUTINE interptsn_sponge ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before


      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   iku, ikv
      REAL(wp) :: ztsa, zabe1, zabe2, zbtr
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) :: ztu, ztv
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2) ::tsbdiff
      !
      IF (before) THEN
         tabres(i1:i2,j1:j2,k1:k2,n1:n2) = tsn(i1:i2,j1:j2,k1:k2,n1:n2)
      ELSE   
   
         tsbdiff(:,:,:,:) = tsb(i1:i2,j1:j2,:,:) - tabres(:,:,:,:)    
         DO jn = 1, jpts            
            DO jk = 1, jpkm1
               DO jj = j1,j2-1
                  DO ji = i1,i2-1
                     zabe1 = fsaht_spu(ji,jj) * umask(ji,jj,jk) * re2u_e1u(ji,jj) * fse3u_n(ji,jj,jk)
                     zabe2 = fsaht_spv(ji,jj) * vmask(ji,jj,jk) * re1v_e2v(ji,jj) * fse3v_n(ji,jj,jk)
                     ztu(ji,jj,jk) = zabe1 * ( tsbdiff(ji+1,jj  ,jk,jn) - tsbdiff(ji,jj,jk,jn) ) 
                     ztv(ji,jj,jk) = zabe2 * ( tsbdiff(ji  ,jj+1,jk,jn) - tsbdiff(ji,jj,jk,jn) )
                  ENDDO
               ENDDO

               IF( ln_zps ) THEN      ! set gradient at partial step level
                  DO jj = j1,j2-1
                     DO ji = i1,i2-1
                        ! last level
                        iku = mbku(ji,jj)
                        ikv = mbkv(ji,jj)
                        IF( iku == jk ) THEN
                           ztu(ji,jj,jk) = 0._wp
                        ENDIF
                        IF( ikv == jk ) THEN
                           ztv(ji,jj,jk) = 0._wp
                        ENDIF
                     END DO
                  END DO
               ENDIF
            ENDDO

            DO jk = 1, jpkm1
               DO jj = j1+1,j2-1
                  DO ji = i1+1,i2-1

                     IF (.NOT. tabspongedone_tsn(ji,jj)) THEN 
                        zbtr = r1_e12t(ji,jj) / fse3t_n(ji,jj,jk)
                        ! horizontal diffusive trends
                        ztsa = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk)  )
                        ! add it to the general tracer trends
                        tsa(ji,jj,jk,jn) = tsa(ji,jj,jk,jn) + ztsa
                     ENDIF

                  ENDDO
               ENDDO

            ENDDO
         ENDDO

         tabspongedone_tsn(i1+1:i2-1,j1+1:j2-1) = .TRUE.

      ENDIF

   END SUBROUTINE interptsn_sponge

   SUBROUTINE interpun_sponge(tabres,i1,i2,j1,j2,k1,k2, before)
      !!---------------------------------------------
      !!   *** ROUTINE interpun_sponge ***
      !!---------------------------------------------    
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before

      INTEGER :: ji,jj,jk

      ! sponge parameters 
      REAL(wp) :: ze2u, ze1v, zua, zva, zbtr
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) :: ubdiff
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) :: rotdiff, hdivdiff
      INTEGER :: jmax
      !


      IF (before) THEN
         tabres = un(i1:i2,j1:j2,:)
      ELSE

         ubdiff(i1:i2,j1:j2,:) = (ub(i1:i2,j1:j2,:) - tabres(:,:,:))*umask(i1:i2,j1:j2,:)

         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============

            !                                             ! --------
            ! Horizontal divergence                       !   div
            !                                             ! --------
            DO jj = j1,j2
               DO ji = i1+1,i2   ! vector opt.
                  zbtr = r1_e12t(ji,jj) / fse3t_n(ji,jj,jk) * fsahm_spt(ji,jj)
                  hdivdiff(ji,jj,jk) = (  e2u(ji  ,jj)*fse3u_n(ji  ,jj,jk) * ubdiff(ji  ,jj,jk) &
                                     &   -e2u(ji-1,jj)*fse3u_n(ji-1,jj,jk) * ubdiff(ji-1,jj,jk) ) * zbtr
               END DO
            END DO

            DO jj = j1,j2-1
               DO ji = i1,i2   ! vector opt.
                  zbtr = r1_e12f(ji,jj) * fse3f_n(ji,jj,jk) * fsahm_spf(ji,jj)
                  rotdiff(ji,jj,jk) = (-e1u(ji,jj+1) * ubdiff(ji,jj+1,jk) &
                                       +e1u(ji,jj  ) * ubdiff(ji,jj  ,jk) & 
                                    & ) * fmask(ji,jj,jk) * zbtr 
               END DO
            END DO
         ENDDO

         !



         DO jj = j1+1, j2-1
            DO ji = i1+1, i2-1   ! vector opt.

               IF (.NOT. tabspongedone_u(ji,jj)) THEN
                  DO jk = 1, jpkm1                                 ! Horizontal slab
                     ze2u = rotdiff (ji,jj,jk)
                     ze1v = hdivdiff(ji,jj,jk)
                     ! horizontal diffusive trends
                     zua = - ( ze2u - rotdiff (ji,jj-1,jk)) / ( e2u(ji,jj) * fse3u_n(ji,jj,jk) )   &
                           + ( hdivdiff(ji+1,jj,jk) - ze1v  ) / e1u(ji,jj)

                     ! add it to the general momentum trends
                     ua(ji,jj,jk) = ua(ji,jj,jk) + zua

                  END DO
               ENDIF

            END DO
         END DO

         tabspongedone_u(i1+1:i2-1,j1+1:j2-1) = .TRUE.

         jmax = j2-1
         IF ((nbondj == 1).OR.(nbondj == 2)) jmax = MIN(jmax,nlcj-3)

         DO jj = j1+1, jmax
            DO ji = i1+1, i2   ! vector opt.

               IF (.NOT. tabspongedone_v(ji,jj)) THEN
                  DO jk = 1, jpkm1                                 ! Horizontal slab
                     ze2u = rotdiff (ji,jj,jk)
                     ze1v = hdivdiff(ji,jj,jk)

                     ! horizontal diffusive trends
                     zva = + ( ze2u - rotdiff (ji-1,jj,jk)) / ( e1v(ji,jj) * fse3v_n(ji,jj,jk) )   &
                           + ( hdivdiff(ji,jj+1,jk) - ze1v  ) / e2v(ji,jj)

                     ! add it to the general momentum trends
                     va(ji,jj,jk) = va(ji,jj,jk) + zva
                  END DO
               ENDIF

            END DO
         END DO


         tabspongedone_v(i1+1:i2,j1+1:jmax) = .TRUE.

      ENDIF


   END SUBROUTINE interpun_sponge


   SUBROUTINE interpvn_sponge(tabres,i1,i2,j1,j2,k1,k2, before,nb,ndir)
      !!---------------------------------------------
      !!   *** ROUTINE interpvn_sponge ***
      !!--------------------------------------------- 
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir

      INTEGER :: ji,jj,jk

      REAL(wp) :: ze2u, ze1v, zua, zva, zbtr

      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) :: vbdiff
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) :: rotdiff, hdivdiff
      INTEGER :: imax
      !

      IF (before) THEN 
         tabres = vn(i1:i2,j1:j2,:)
      ELSE

         vbdiff(i1:i2,j1:j2,:) = (vb(i1:i2,j1:j2,:) - tabres(:,:,:))*vmask(i1:i2,j1:j2,:)

         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============

            !                                             ! --------
            ! Horizontal divergence                       !   div
            !                                             ! --------
            DO jj = j1+1,j2
               DO ji = i1,i2   ! vector opt.
                  zbtr = r1_e12t(ji,jj) / fse3t_n(ji,jj,jk) * fsahm_spt(ji,jj)
                  hdivdiff(ji,jj,jk) = ( e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * vbdiff(ji,jj  ,jk)  &
                                     &  -e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * vbdiff(ji,jj-1,jk)  ) * zbtr
               END DO
            END DO
            DO jj = j1,j2
               DO ji = i1,i2-1   ! vector opt.
                  zbtr = r1_e12f(ji,jj) * fse3f_n(ji,jj,jk) * fsahm_spf(ji,jj)
                  rotdiff(ji,jj,jk) = ( e2v(ji+1,jj) * vbdiff(ji+1,jj,jk) & 
                                    &  -e2v(ji  ,jj) * vbdiff(ji  ,jj,jk) &
                                    & ) * fmask(ji,jj,jk) * zbtr
               END DO
            END DO
         ENDDO

         !                                                ! ===============
         !                                                

         imax = i2-1
         IF ((nbondi == 1).OR.(nbondi == 2)) imax = MIN(imax,nlci-3)

         DO jj = j1+1, j2
            DO ji = i1+1, imax   ! vector opt.
               IF (.NOT. tabspongedone_u(ji,jj)) THEN
                  DO jk = 1, jpkm1                                 ! Horizontal slab
                     ze2u = rotdiff (ji,jj,jk)
                     ze1v = hdivdiff(ji,jj,jk)
                     ! horizontal diffusive trends
                     zua = - ( ze2u - rotdiff (ji,jj-1,jk)) / ( e2u(ji,jj) * fse3u_n(ji,jj,jk) ) + ( hdivdiff(ji+1,jj,jk) - ze1v) &
                           / e1u(ji,jj)


                     ! add it to the general momentum trends
                     ua(ji,jj,jk) = ua(ji,jj,jk) + zua
                  END DO

               ENDIF
            END DO
         END DO

         tabspongedone_u(i1+1:imax,j1+1:j2) = .TRUE.

         DO jj = j1+1, j2-1
            DO ji = i1+1, i2-1   ! vector opt.
               IF (.NOT. tabspongedone_v(ji,jj)) THEN
                  DO jk = 1, jpkm1                                 ! Horizontal slab
                     ze2u = rotdiff (ji,jj,jk)
                     ze1v = hdivdiff(ji,jj,jk)
                     ! horizontal diffusive trends

                     zva = + ( ze2u - rotdiff (ji-1,jj,jk)) / ( e1v(ji,jj) * fse3v_n(ji,jj,jk) ) + ( hdivdiff(ji,jj+1,jk) - ze1v) &
                           / e2v(ji,jj)

                     ! add it to the general momentum trends
                     va(ji,jj,jk) = va(ji,jj,jk) + zva
                  END DO
               ENDIF
            END DO
         END DO
         tabspongedone_v(i1+1:i2-1,j1+1:j2-1) = .TRUE.
      ENDIF

   END SUBROUTINE interpvn_sponge

#else
CONTAINS

   SUBROUTINE agrif_opa_sponge_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_OPA_sponge_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_opa_sponge : You should not have seen this print! error?'
   END SUBROUTINE agrif_opa_sponge_empty
#endif

END MODULE agrif_opa_sponge
