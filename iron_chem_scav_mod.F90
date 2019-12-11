MODULE iron_chem_scav_mod
   !!======================================================================
   !!                         ***  MODULE iron_chem_scav_mod  ***
   !! Calculate the iron chemistry and scavenging.
   !!======================================================================
   !! History :
   !!   -   ! 2017-04 (M. Stringer)        Code taken from trcbio_medusa.F90
   !!----------------------------------------------------------------------
#if defined key_medusa
   !!----------------------------------------------------------------------
   !!                                                   MEDUSA bio-model
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE
      
   PUBLIC   iron_chem_scav        ! Called in trcbio_medusa.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE iron_chem_scav( jk )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE iron_chem_scav  ***
      !! This called from TRC_BIO_MEDUSA and 
      !!  - 
      !!-------------------------------------------------------------------
      USE bio_medusa_mod,    ONLY: ffastc, ffastca, ffastsi,              &
                                   ffetop, ffebot, ffescav, xfree,        & 
                                   zdet, zfer, zphd, zphn, zzme, zzmi,    &
                                   idf, idfval                          
      USE dom_oce,           ONLY: e3t_0, gdepw_0, mbathy, tmask
#if defined key_vvl
      USE dom_oce,           ONLY: e3t_n, gdepw_n
#endif
      USE par_kind,          ONLY: wp
      USE in_out_manager,    ONLY: lwp, numout
      USE par_oce,           ONLY: jpi, jpim1, jpj, jpjm1
      USE sms_medusa,        ONLY: i0500, jiron, xfe_sed, xfe_sol,        &
                                   xfe_mass,                              &
                                   xk_FeL, xk_sc_Fe, xLgT,                &
                                   xmassc, xmassca, xmasssi,              &
                                   xthetad, xthetapd, xthetapn,           &
                                   xthetazme, xthetazmi,                  &
                                   zirondep

   !!* Substitution
#  include "domzgr_substitute.h90"

      !! Level
      INTEGER, INTENT( in ) :: jk

      !! iron cycle; includes parameters for Parekh et al. (2005) iron scheme
      !! state variables for iron-ligand system
      REAL(wp), DIMENSION(jpi,jpj) :: xFeT, xFeF, xFeL
      REAL(wp) :: xLgF
      !! iron-ligand parameters
      REAL(wp) :: xb_coef_tmp, xb2M4ac
      !! max Fe' parameters
      REAL(wp) :: xmaxFeF,fdeltaFe
      !!
      !! local parameters for Moore et al. (2004) alternative scavenging 
      !! scheme
      REAL(wp) :: fbase_scav,fscal_sink,fscal_part,fscal_scav
      !!
      !! local parameters for Moore et al. (2008) alternative scavenging 
      !! scheme
      REAL(wp) :: fscal_csink,fscal_sisink,fscal_casink
      !!
      !! local parameters for Galbraith et al. (2010) alternative 
      !! scavenging scheme.
      !! organic portion of scavenging
      REAL(wp) :: xCscav1, xCscav2, xk_org, xORGscav
      !! inorganic portion of scavenging
      REAL(wp) :: xk_inorg, xINORGscav

      INTEGER :: ji, jj

      !!------------------------------------------------------------------
      !! Iron chemistry and fractionation
      !! following the Parekh et al. (2004) scheme adopted by the Met.
      !! Office, Medusa models total iron but considers "free" and
      !! ligand-bound forms for the purposes of scavenging (only "free"
      !! iron can be scavenged
      !!------------------------------------------------------------------
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            !! OPEN wet point IF..THEN loop
            if (tmask(ji,jj,jk) == 1) then
               !!
               !! total iron concentration (mmol Fe / m3 -> umol Fe / m3)
               xFeT(ji,jj) = zfer(ji,jj) * 1.e3
               !!
               !! calculate fractionation (based on Diat-HadOCC; in turn 
               !! based on Parekh et al., 2004)
               xb_coef_tmp = xk_FeL * (xLgT - xFeT(ji,jj)) - 1.0
               xb2M4ac     = max(((xb_coef_tmp * xb_coef_tmp) +              &
                                  (4.0 * xk_FeL * xLgT)), 0.0)
               !!
               !! "free" ligand concentration
               xLgF        = 0.5 * (xb_coef_tmp + (xb2M4ac**0.5)) / xk_FeL
               !!
               !! ligand-bound iron concentration
               xFeL(ji,jj) = xLgT - xLgF
               !!
               !! "free" iron concentration (and convert to mmol Fe / m3)
               xFeF(ji,jj) = (xFeT(ji,jj) - xFeL(ji,jj)) * 1.e-3
               xFree(ji,jj)= xFeF(ji,jj) / (zfer(ji,jj) + tiny(zfer(ji,jj)))
            ENDIF
         ENDDO
      ENDDO


      !!
      !! scavenging of iron (multiple schemes); I'm only really 
      !! happy with the first one at the moment - the others 
      !! involve assumptions (sometimes guessed at by me) that 
      !! are potentially questionable
      !!
      if (jiron.eq.1) then
         !!------------------------------------------------------
         !! Scheme 1: Dutkiewicz et al. (2005)
         !! This scheme includes a single scavenging term based 
         !! solely on a fixed rate and the availablility of 
         !! "free" iron
         !!------------------------------------------------------
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  !! = mmol/m3/d
                  ffescav(ji,jj)     = xk_sc_Fe * xFeF(ji,jj)
                  !!
                  !!------------------------------------------------------
                  !!
                  !! Mick's code contains a further (optional) implicit 
                  !! "scavenging" of iron that sets an upper bound on 
                  !! "free" iron concentration, and essentially caps the 
                  !! concentration of total iron as xFeL + "free" iron; 
                  !! since the former is constrained by a fixed total 
                  !! ligand concentration (= 1.0 umol/m3), and the latter 
                  !! isn't allowed above this upper bound, total iron is 
                  !! constrained to a maximum of ...
                  !!
                  !!    xFeL(ji,jj) + min(xFeF(ji,jj), 0.3 umol/m3) = 1.0 + 0.3 
                  !!                                  = 1.3 umol / m3
                  !! 
                  !! In Mick's code, the actual value of total iron is 
                  !! reset to this sum (i.e. TFe = FeL + Fe'; but 
                  !! Fe' <= 0.3 umol/m3); this isn't our favoured approach 
                  !! to tracer updating here (not least because of the 
                  !! leapfrog), so here the amount scavenged is augmented 
                  !! by an additional amount that serves to drag total 
                  !! iron back towards that expected from this limitation 
                  !! on iron concentration ...
                  !!
                  !! = umol/m3
                  xmaxFeF     = min((xFeF(ji,jj) * 1.e3), 0.3)
                  !!
                  !! Here, the difference between current total Fe and 
                  !! (FeL + Fe') is calculated and added to the scavenging 
                  !! flux already calculated above ...
                  !!
                  !! = mmol/m3
                  fdeltaFe    = (xFeT(ji,jj) - (xFeL(ji,jj) + xmaxFeF)) * 1.e-3
                  !!
                  !! This assumes that the "excess" iron is dissipated 
                  !! with a time-scale of 1 day; seems reasonable to me 
                  !! ... (famous last words)
                  !!
                  !! = mmol/m3/d
                  ffescav(ji,jj)     = ffescav(ji,jj) + fdeltaFe
                  !!
# if defined key_deep_fe_fix
                  !! AXY (17/01/13)
                  !! stop scavenging for iron concentrations below 
                  !! 0.5 umol / m3 at depths greater than 1000 m; this 
                  !! aims to end MEDUSA's continual loss of iron at depth 
                  !! without impacting things at the surface too much; the 
                  !! justification for this is that it appears to be what 
                  !! Mick Follows et al. do in their work (as evidenced by 
                  !! the iron initial condition they supplied me with); to 
                  !! be honest, it looks like Follow et al. do this at
                  !! shallower depths than 1000 m, but I'll stick with this
                  !! for now; I suspect that this seemingly arbitrary
                  !! approach effectively "parameterises" the 
                  !! particle-based scavenging rates that other models use 
                  !! (i.e. at depth there are no sinking particles, so 
                  !! scavenging stops); it might be fun justifying this in 
                  !! a paper though!
                  !!
                  if ((fsdepw(ji,jj,jk).gt.1000.) .and.                       &
                       (xFeT(ji,jj).lt.0.5)) then
                     ffescav(ji,jj) = 0.
                  endif
# endif
               ENDIF
            ENDDO
         ENDDO
      elseif (jiron.eq.2) then
         !!------------------------------------------------------
         !! Scheme 2: Moore et al. (2004)
         !! This scheme includes a single scavenging term that 
         !! accounts for both suspended and sinking particles in 
         !! the water column; this term scavenges total iron rather 
         !! than "free" iron
         !!------------------------------------------------------
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  !!
                  !! total iron concentration (mmol Fe / m3 -> umol Fe / m3)
                  xFeT(ji,jj) = zfer(ji,jj) * 1.e3
                  !!
                  !! this has a base scavenging rate (12% / y) which is 
                  !! modified by local particle concentration and sinking 
                  !! flux (and dust - but I'm ignoring that here for now) 
                  !! and which is accelerated when Fe concentration gets
                  !! 0.6 nM (= 0.6 umol/m3 = 0.0006 mmol/m3), and decreased 
                  !! as concentrations below 0.4 nM (= 0.4 umol/m3 = 
                  !! 0.0004 mmol/m3)
                  !!
                  !! base scavenging rate (0.12 / y)
                  fbase_scav = 0.12 / 365.25
                  !!
                  !! calculate sinking particle part of scaling factor
                  !! this takes local fast sinking carbon (mmol C / m2 / d)
                  !! and gets it into nmol C / cm3 / s ("rdt" below is the 
                  !! number of seconds in a model timestep)
                  !!
                  !! fscal_sink = ffastc(ji,jj) * 1.e2 / (86400.)
                  !!
                  !! ... actually, re-reading Moore et al.'s equations, it 
                  !! looks like he uses his sinking flux directly, without 
                  !! scaling it by time-step or anything, so I'll copy this 
                  !! here ...
                  !!
                  fscal_sink = ffastc(ji,jj) * 1.e2
                  !!
                  !! calculate particle part of scaling factor
                  !! this totals up the carbon in suspended particles 
                  !! (Pn, Pd, Zmi, Zme, D),
                  !! which comes out in mmol C / m3 (= nmol C / cm3), and 
                  !! then multiplies it by a magic factor, 0.002, to get it 
                  !! into nmol C / cm2 / s
                  !!
                  fscal_part = ( (xthetapn * zphn(ji,jj)) +                  &
                                 (xthetapd * zphd(ji,jj)) +                  &
                                 (xthetazmi * zzmi(ji,jj)) +                 &
                                 (xthetazme * zzme(ji,jj)) +                 &
                                 (xthetad * zdet(ji,jj)) ) * 0.002
                  !!
                  !! calculate scaling factor for base scavenging rate
                  !! this uses the (now correctly scaled) sinking flux and 
                  !! standing
                  !! particle concentration, divides through by some sort 
                  !! of reference value (= 0.0066 nmol C / cm2 / s) and 
                  !! then uses this, or not if its too high, to rescale the 
                  !! base scavenging rate
                  !!
                  fscal_scav = fbase_scav *                                  &
                               min(((fscal_sink + fscal_part) / 0.0066), 4.0)
                  !!
                  !! the resulting scavenging rate is then scaled further 
                  !! according to the local iron concentration (i.e. 
                  !! diminished in low iron regions; enhanced in high iron 
                  !! regions; less alone in intermediate iron regions)
                  !!
                  if (xFeT(ji,jj).lt.0.4) then
                     !!
                     !! low iron region
                     !!
                     fscal_scav = fscal_scav * (xFeT(ji,jj) / 0.4)
                     !!
                  elseif (xFeT(ji,jj).gt.0.6) then
                     !!
                     !! high iron region
                     !!
                     fscal_scav = fscal_scav + ((xFeT(ji,jj) / 0.6) *        &
                                                (6.0 / 1.4))
                     !!
                  else
                     !!
                     !! intermediate iron region: do nothing
                     !!
                  endif
                  !! 
                  !! apply the calculated scavenging rate ...
                  !!
                  ffescav(ji,jj) = fscal_scav * zfer(ji,jj)
                  !!
               ENDIF
            ENDDO
         ENDDO
      elseif (jiron.eq.3) then
         !!------------------------------------------------------
         !! Scheme 3: Moore et al. (2008)
         !! This scheme includes a single scavenging term that 
         !! accounts for sinking particles in the water column, 
         !! and includes organic C, biogenic opal, calcium 
         !! carbonate and dust in this (though the latter is 
         !! ignored here until I work out what units the incoming
         !! "dust" flux is in); this term scavenges total iron 
         !! rather than "free" iron
         !!------------------------------------------------------
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  !!
                  !! total iron concentration (mmol Fe / m3 -> umol Fe / m3)
                  xFeT(ji,jj) = zfer(ji,jj) * 1.e3
                  !!
                  !! this has a base scavenging rate which is modified by 
                  !! local particle sinking flux (including dust - but I'm 
                  !! ignoring that here for now) and which is accelerated 
                  !! when Fe concentration is > 0.6 nM (= 0.6 umol/m3 = 
                  !! 0.0006 mmol/m3), and decreased as concentrations < 
                  !! 0.5 nM (= 0.5 umol/m3 = 0.0005 mmol/m3)
                  !!
                  !! base scavenging rate (Fe_b in paper; units may be 
                  !! wrong there)
                  fbase_scav = 0.00384 ! (ng)^-1 cm
                  !!
                  !! calculate sinking particle part of scaling factor; 
                  !! this converts mmol / m2 / d fluxes of organic carbon, 
                  !! silicon and calcium carbonate into ng / cm2 / s 
                  !! fluxes; it is assumed here that the mass conversions 
                  !! simply consider the mass of the main element
                  !! (C, Si and Ca) and *not* the mass of the molecules 
                  !! that they are part of; Moore et al. (2008) is unclear 
                  !! on the conversion that should be used
                  !!
                  !! milli -> nano; mol -> gram; /m2 -> /cm2; /d -> /s
                  !! ng C  / cm2 / s
                  fscal_csink  = (ffastc(ji,jj)  * 1.e6 * xmassc  *          &
                                  1.e-4 / 86400.)
                  !! ng Si / cm2 / s
                  fscal_sisink = (ffastsi(ji,jj) * 1.e6 * xmasssi *          &
                                  1.e-4 / 86400.)
                  !! ng Ca / cm2 / s
                  fscal_casink = (ffastca(ji,jj) * 1.e6 * xmassca *          &
                                  1.e-4 / 86400.)
                  !! 
                  !! sum up these sinking fluxes and convert to ng / cm 
                  !! by dividing through by a sinking rate of 
                  !! 100 m / d = 1.157 cm / s
                  !! ng / cm
                  fscal_sink   = ((fscal_csink * 6.) + fscal_sisink +        &
                                  fscal_casink) / (100. * 1.e3 / 86400)
                  !!
                  !! now calculate the scavenging rate based upon the base 
                  !! rate and this particle flux scaling; according to the 
                  !! published units, the result actually has *no* units, 
                  !! but as it must be expressed per unit time for it to 
                  !! make any sense, I'm assuming a missing "per second"
                  !! / s
                  fscal_scav = fbase_scav * fscal_sink
                  !!
                  !! the resulting scavenging rate is then scaled further 
                  !! according to the local iron concentration (i.e. 
                  !! diminished in low iron regions; enhanced in high iron 
                  !! regions; less alone in intermediate iron regions)
                  !!
                  if (xFeT(ji,jj).lt.0.5) then
                     !!
                     !! low iron region (0.5 instead of the 0.4 in Moore 
                     !! et al., 2004)
                     !!
                     fscal_scav = fscal_scav * (xFeT(ji,jj) / 0.5)
                     !!
                  elseif (xFeT(ji,jj).gt.0.6) then
                     !!
                     !! high iron region (functional form different in 
                     !! Moore et al., 2004)
                     !!
                     fscal_scav = fscal_scav + ((xFeT(ji,jj) - 0.6) * 0.00904)
                     !!
                  else
                     !!
                     !! intermediate iron region: do nothing
                     !!
                  endif
                  !! 
                  !! apply the calculated scavenging rate ...
                  !!
                  ffescav(ji,jj) = fscal_scav * zfer(ji,jj)
               ENDIF
            ENDDO
         ENDDO
      elseif (jiron.eq.4) then
         !!------------------------------------------------------
         !! Scheme 4: Galbraith et al. (2010)
         !! This scheme includes two scavenging terms, one for 
         !! organic, particle-based scavenging, and another for 
         !! inorganic scavenging; both terms scavenge "free" iron 
         !! only
         !!------------------------------------------------------
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  !!
                  !! Galbraith et al. (2010) present a more straightforward 
                  !! outline of the scheme in Parekh et al. (2005) ...
                  !! 
                  !! sinking particulate carbon available for scavenging
                  !! this assumes a sinking rate of 100 m / d (Moore & 
                  !! Braucher, 2008),
                  xCscav1    = (ffastc(ji,jj) * xmassc) / 100. ! = mg C / m3
                  !! 
                  !! scale by Honeyman et al. (1981) exponent coefficient
                  !! multiply by 1.e-3 to express C flux in g C rather than 
                  !! mg C
                  xCscav2    = (xCscav1 * 1.e-3)**0.58
                  !!
                  !! multiply by Galbraith et al. (2010) scavenging rate
                  xk_org     = 0.5 ! ((g C m/3)^-1) / d
                  xORGscav   = xk_org * xCscav2 * xFeF(ji,jj)
                  !!
                  !! Galbraith et al. (2010) also include an inorganic bit ...
                  !! 
                  !! this occurs at a fixed rate, again based on the 
                  !! availability of "free" iron
                  !!
                  !! k_inorg = 1000 d**-1 nmol Fe**-0.5 kg**-0.5
                  !!
                  !! to implement this here, scale xFeF by 1026 to put in 
                  !! units of umol/m3 which approximately equal nmol/kg
                  !!
                  xk_inorg   = 1000. ! ((nmol Fe / kg)^1.5)
                  xINORGscav = (xk_inorg * (xFeF(ji,jj) * 1026.)**1.5) * 1.e-3
                  !!
                  !! sum these two terms together
                  ffescav(ji,jj) = xORGscav + xINORGscav
               ENDIF
            ENDDO
         ENDDO
      else
         !!------------------------------------------------------
         !! No Scheme: you coward!
         !! This scheme puts its head in the sand and eskews any 
         !! decision about how iron is removed from the ocean; 
         !! prepare to get deluged in iron you fool!
         !!------------------------------------------------------
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  ffescav(ji,jj) = 0.
               ENDIF
            ENDDO
         ENDDO
      endif

      !!---------------------------------------------------------
      !! Other iron cycle processes
      !!---------------------------------------------------------
      !!
      !! aeolian iron deposition
      !! zirondep      is in mmol-Fe / m2 / day
      !! ffetop(ji,jj) is in mmol-dissolved-Fe / m3 / day
      if (jk == 1) then
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  ffetop(ji,jj)  = zirondep(ji,jj) * xfe_sol / fse3t(ji,jj,jk) 
               ENDIF
            ENDDO
         ENDDO
      else
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  ffetop(ji,jj)  = 0.0
                ENDIF
            ENDDO
         ENDDO
      endif
      !!
      !! seafloor iron addition
      DO jj = 2,jpjm1
         DO ji = 2,jpim1
            IF (tmask(ji,jj,jk) == 1) THEN
               !! AXY (10/07/12): amended to only apply sedimentary flux up 
               !! to ~500 m down
               !! if (jk.eq.(mbathy(ji,jj)-1).AND.jk.lt.i1100) then
               if ((jk.eq.mbathy(ji,jj)).AND.jk.le.i0500) then
                  !! Moore et al. (2004) cite a coastal California value of 
                  !! 5 umol/m2/d, but adopt a global value of 2 umol/m2/d 
                  !! for all areas < 1100 m; here we use this latter value
                  !! but apply it everywhere
                  !! AXY (21/07/09): actually, let's just apply it below 
                  !! 1100 m (levels 1-37)
                  ffebot(ji,jj)  = (xfe_sed / fse3t(ji,jj,jk))
               else
                  ffebot(ji,jj)  = 0.0
               endif
            ENDIF
         ENDDO
      ENDDO

      !! AXY (16/12/09): remove iron addition/removal processes
      !! For the purposes of the quarter degree run, the iron 
      !! cycle is being pegged to the initial condition supplied 
      !! by Mick Follows via restoration with a 30 day period;
      !! iron addition at the seafloor is still permitted with 
      !! the idea that this extra iron will be removed by the 
      !! restoration away from the source
      !! ffescav(ji,jj) = 0.0
      !! ffetop(ji,jj)  = 0.0
      !! ffebot(ji,jj)  = 0.0

# if defined key_debug_medusa
      !! report miscellaneous calculations
      !! report miscellaneous calculations
      if (idf.eq.1.AND.idfval.eq.1) then
         DO jj = 2,jpjm1
            DO ji = 2,jpim1
               IF (tmask(ji,jj,jk) == 1) THEN
                  IF (lwp) write (numout,*) '------------------------------'
                  IF (lwp) write (numout,*) 'xfe_sol  = ', xfe_sol
                  IF (lwp) write (numout,*) 'xfe_mass = ', xfe_mass
                  IF (lwp) write (numout,*) 'ffetop(',jk,')  = ', ffetop(ji,jj)
                  IF (lwp) write (numout,*) 'ffebot(',jk,')  = ', ffebot(ji,jj)
                  IF (lwp) write (numout,*) 'xFree(',jk,')   = ', xFree(ji,jj)
                  IF (lwp) write (numout,*) 'ffescav(',jk,') = ', ffescav(ji,jj)
               ENDIF
            ENDDO
         ENDDO
      endif
# endif

   END SUBROUTINE iron_chem_scav

#else
   !!======================================================================
   !!  Dummy module :                                   No MEDUSA bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE iron_chem_scav( )                    ! Empty routine
      WRITE(*,*) 'iron_chem_scav: You should not have seen this print! error?'
   END SUBROUTINE iron_chem_scav
#endif 

   !!======================================================================
END MODULE iron_chem_scav_mod
