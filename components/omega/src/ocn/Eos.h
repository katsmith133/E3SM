#ifndef OMEGA_EOS_H
#define OMEGA_EOS_H
//===-- ocn/Eos.h - Equation of State --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating specific volume
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate the specific volume based on the choice of EOS
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

enum class EosType {
   LinearEos,  /// Linear equation of state
   Teos10Eos,  /// Roquet et al. 2015 75 term expansion
   ConstantEos /// Constant specific volume equation of state
};

/// TEOS10 75-term Polynomial Equation of State
class Teos10Eos {
 public:
   /// constructor declaration
   Teos10Eos(const VertCoord *VCoord);

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   specific volume according to the Roquet et al. 2015 75 term expansion.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity,
                                   const Array2DReal &Pressure,
                                   I4 KDisp) const {

      Real SpecVolPCoeffs[6 * VecLength];
      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         /// Calculate the local specific volume polynomial pressure
         /// coefficients with cell center values
         calcPCoeffs(SpecVolPCoeffs, KVec, ConservTemp(ICell, K),
                     AbsSalinity(ICell, K));

         /// Calculate the specific volume at the given pressure
         /// If KDisp is 0, we use the current pressure, otherwise
         /// we use the displaced pressure (K + KDisp)
         /// Note: KDisp is only used for TEOS-10, for Linear EOS it
         /// is always 0.
         if (KDisp == 0) {
            // No displacement
            SpecVol(ICell, K) =
                calcRefProfile(Pressure(ICell, K) * Pa2Db) +
                calcDelta(SpecVolPCoeffs, KVec, Pressure(ICell, K) * Pa2Db);
         } else {
            // Displacement, use the displaced pressure
            I4 KTmp = Kokkos::min(K + KDisp, MaxLayerCell(ICell));
            KTmp    = Kokkos::max(MinLayerCell(ICell), KTmp);
            SpecVol(ICell, K) =
                calcRefProfile(Pressure(ICell, KTmp) * Pa2Db) +
                calcDelta(SpecVolPCoeffs, KVec, Pressure(ICell, KTmp) * Pa2Db);
         }
      }
   }

   /// TEOS-10 helpers
   /// Calculate pressure polynomial coefficients for TEOS-10
   KOKKOS_FUNCTION void calcPCoeffs(Real (&SpecVolPCoeffs)[6 * VecLength],
                                    const I4 KVec, const Real Ct,
                                    const Real Sa) const {
      constexpr Real SaNorm = 40.0 * 35.16504 / 35.0;
      constexpr Real CtNorm = 40.0;
      constexpr Real DeltaS = 24.0;
      Real Ss               = Kokkos::sqrt((Sa + DeltaS) / SaNorm);
      Real Tt               = Ct / CtNorm;

      /// Coefficients for the polynomial expansion
      constexpr Real V000 = 1.0769995862e-03;
      constexpr Real V100 = -3.1038981976e-04;
      constexpr Real V200 = 6.6928067038e-04;
      constexpr Real V300 = -8.5047933937e-04;
      constexpr Real V400 = 5.8086069943e-04;
      constexpr Real V500 = -2.1092370507e-04;
      constexpr Real V600 = 3.1932457305e-05;
      constexpr Real V010 = -1.5649734675e-05;
      constexpr Real V110 = 3.5009599764e-05;
      constexpr Real V210 = -4.3592678561e-05;
      constexpr Real V310 = 3.4532461828e-05;
      constexpr Real V410 = -1.1959409788e-05;
      constexpr Real V510 = 1.3864594581e-06;
      constexpr Real V020 = 2.7762106484e-05;
      constexpr Real V120 = -3.7435842344e-05;
      constexpr Real V220 = 3.5907822760e-05;
      constexpr Real V320 = -1.8698584187e-05;
      constexpr Real V420 = 3.8595339244e-06;
      constexpr Real V030 = -1.6521159259e-05;
      constexpr Real V130 = 2.4141479483e-05;
      constexpr Real V230 = -1.4353633048e-05;
      constexpr Real V330 = 2.2863324556e-06;
      constexpr Real V040 = 6.9111322702e-06;
      constexpr Real V140 = -8.7595873154e-06;
      constexpr Real V240 = 4.3703680598e-06;
      constexpr Real V050 = -8.0539615540e-07;
      constexpr Real V150 = -3.3052758900e-07;
      constexpr Real V060 = 2.0543094268e-07;
      constexpr Real V001 = -1.6784136540e-05;
      constexpr Real V101 = 2.4262468747e-05;
      constexpr Real V201 = -3.4792460974e-05;
      constexpr Real V301 = 3.7470777305e-05;
      constexpr Real V401 = -1.7322218612e-05;
      constexpr Real V501 = 3.0927427253e-06;
      constexpr Real V011 = 1.8505765429e-05;
      constexpr Real V111 = -9.5677088156e-06;
      constexpr Real V211 = 1.1100834765e-05;
      constexpr Real V311 = -9.8447117844e-06;
      constexpr Real V411 = 2.5909225260e-06;
      constexpr Real V021 = -1.1716606853e-05;
      constexpr Real V121 = -2.3678308361e-07;
      constexpr Real V221 = 2.9283346295e-06;
      constexpr Real V321 = -4.8826139200e-07;
      constexpr Real V031 = 7.9279656173e-06;
      constexpr Real V131 = -3.4558773655e-06;
      constexpr Real V231 = 3.1655306078e-07;
      constexpr Real V041 = -3.4102187482e-06;
      constexpr Real V141 = 1.2956717783e-06;
      constexpr Real V051 = 5.0736766814e-07;
      constexpr Real V002 = 3.0623833435e-06;
      constexpr Real V102 = -5.8484432984e-07;
      constexpr Real V202 = -4.8122251597e-06;
      constexpr Real V302 = 4.9263106998e-06;
      constexpr Real V402 = -1.7811974727e-06;
      constexpr Real V012 = -1.1736386731e-06;
      constexpr Real V112 = -5.5699154557e-06;
      constexpr Real V212 = 5.4620748834e-06;
      constexpr Real V312 = -1.3544185627e-06;
      constexpr Real V022 = 2.1305028740e-06;
      constexpr Real V122 = 3.9137387080e-07;
      constexpr Real V222 = -6.5731104067e-07;
      constexpr Real V032 = -4.6132540037e-07;
      constexpr Real V132 = 7.7618888092e-09;
      constexpr Real V042 = -6.3352916514e-08;
      constexpr Real V003 = -3.8088938393e-07;
      constexpr Real V103 = 3.6310188515e-07;
      constexpr Real V203 = 1.6746303780e-08;
      constexpr Real V013 = -3.6527006553e-07;
      constexpr Real V113 = -2.7295696237e-07;
      constexpr Real V023 = 2.8695905159e-07;
      constexpr Real V004 = 8.8302421514e-08;
      constexpr Real V104 = -1.1147125423e-07;
      constexpr Real V014 = 3.1454099902e-07;
      constexpr Real V005 = 4.2369007180e-09;

      SpecVolPCoeffs[5 + 6 * KVec] = V005;
      SpecVolPCoeffs[4 + 6 * KVec] = V014 * Tt + V104 * Ss + V004;
      SpecVolPCoeffs[3 + 6 * KVec] =
          (V023 * Tt + V113 * Ss + V013) * Tt + (V203 * Ss + V103) * Ss + V003;
      SpecVolPCoeffs[2 + 6 * KVec] =
          (((V042 * Tt + V132 * Ss + V032) * Tt + (V222 * Ss + V122) * Ss +
            V022) *
               Tt +
           ((V312 * Ss + V212) * Ss + V112) * Ss + V012) *
              Tt +
          (((V402 * Ss + V302) * Ss + V202) * Ss + V102) * Ss + V002;
      SpecVolPCoeffs[1 + 6 * KVec] =
          ((((V051 * Tt + V141 * Ss + V041) * Tt + (V231 * Ss + V131) * Ss +
             V031) *
                Tt +
            ((V321 * Ss + V221) * Ss + V121) * Ss + V021) *
               Tt +
           (((V411 * Ss + V311) * Ss + V211) * Ss + V111) * Ss + V011) *
              Tt +
          ((((V501 * Ss + V401) * Ss + V301) * Ss + V201) * Ss + V101) * Ss +
          V001;
      SpecVolPCoeffs[0 + 6 * KVec] =
          (((((V060 * Tt + V150 * Ss + V050) * Tt + (V240 * Ss + V140) * Ss +
              V040) *
                 Tt +
             ((V330 * Ss + V230) * Ss + V130) * Ss + V030) *
                Tt +
            (((V420 * Ss + V320) * Ss + V220) * Ss + V120) * Ss + V020) *
               Tt +
           ((((V510 * Ss + V410) * Ss + V310) * Ss + V210) * Ss + V110) * Ss +
           V010) *
              Tt +
          (((((V600 * Ss + V500) * Ss + V400) * Ss + V300) * Ss + V200) * Ss +
           V100) *
              Ss +
          V000;
   }

   /// Evaluate pressure polynomial delta for TEOS-10
   KOKKOS_FUNCTION Real calcDelta(const Real (&SpecVolPCoeffs)[6 * VecLength],
                                  const I4 KVec, const Real P) const {

      constexpr Real PNorm = 1e-4;
      Real Pp              = P * PNorm;

      Real Delta = ((((SpecVolPCoeffs[5 + 6 * KVec] * Pp +
                       SpecVolPCoeffs[4 + 6 * KVec]) *
                          Pp +
                      SpecVolPCoeffs[3 + 6 * KVec]) *
                         Pp +
                     SpecVolPCoeffs[2 + 6 * KVec]) *
                        Pp +
                    SpecVolPCoeffs[1 + 6 * KVec]) *
                       Pp +
                   SpecVolPCoeffs[0 + 6 * KVec];

      return Delta;
   }

   /// Calculate reference profile for TEOS-10
   KOKKOS_FUNCTION Real calcRefProfile(Real P) const {
      constexpr Real PNorm = 1e-4;
      constexpr Real V00   = -4.4015007269e-05;
      constexpr Real V01   = 6.9232335784e-06;
      constexpr Real V02   = -7.5004675975e-07;
      constexpr Real V03   = 1.7009109288e-08;
      constexpr Real V04   = -1.6884162004e-08;
      constexpr Real V05   = 1.9613503930e-09;
      Real Pp              = P * PNorm;

      Real V0 =
          (((((V05 * Pp + V04) * Pp + V03) * Pp + V02) * Pp + V01) * Pp + V00) *
          Pp;
      return V0;
   }

   /// Calculate 2nd derivative of Gibbs wrt pot temp at ref P for TEOS-10
   /// GSW Toolbox function gsw_gibbs_pt0_pt0
   KOKKOS_FUNCTION Real calcGibbsDerivPt0Pt0(Real Sa, Real P) const {
      Real x2 = Sfac * Sa;
      Real x  = Kokkos::sqrt(x2);
      Real y  = P * 0.025;

      Real g03 =
          -24715.571866078 +
          y * (4420.4472249096725 +
               y * (-1778.231237203896 +
                    y * (1160.5182516851419 +
                         y * (-569.531539542516 + y * 128.13429152494615))));

      Real g08 =
          x2 *
          (1760.062705994408 +
           x * (-86.1329351956084 +
                x * (-137.1145018408982 +
                     y * (296.20061691375236 +
                          y * (-205.67709290374563 + 49.9394019139016 * y))) +
                y * (-60.136422517125 + y * 10.50720794170734)) +
           y * (-1351.605895580406 +
                y * (1097.1125373015109 +
                     y * (-433.20648175062206 + 63.905091254154904 * y))));

      return ((g03 + g08) * 0.000625);
   }

   /// Calculate Pot Temmperature from Conservative Temp
   /// GSW Toolbox function gsw_pt_from_ct
   KOKKOS_FUNCTION Real calcPtFromCt(Real Sa, Real Ct) const {
      constexpr Real a0 = -1.446013646344788e-2;
      constexpr Real a1 = -3.305308995852924e-3;
      constexpr Real a2 = 1.062415929128982e-4;
      constexpr Real a3 = 9.477566673794488e-1;
      constexpr Real a4 = 2.166591947736613e-3;
      constexpr Real a5 = 3.828842955039902e-3;
      constexpr Real b0 = 1.000000000000000e0;
      constexpr Real b1 = 6.506097115635800e-4;
      constexpr Real b2 = 3.830289486850898e-3;
      constexpr Real b3 = 1.247811760368034e-6;
      Real a5ct, b3ct, ct_factor, pt_num, pt_recden, ct_diff;
      Real pt, pt_old, ptm, dpt_dct, s1;

      s1 = Sa / Psu2Gpkg;

      a5ct = a5 * Ct;
      b3ct = b3 * Ct;

      ct_factor = (a3 + a4 * s1 + a5ct);
      pt_num    = a0 + s1 * (a1 + a2 * s1) + Ct * ct_factor;
      pt_recden = 1.0 / (b0 + b1 * s1 + Ct * (b2 + b3ct));
      pt        = pt_num * pt_recden;

      dpt_dct = pt_recden * (ct_factor + a5ct - (b2 + b3ct + b3ct) * pt);

      ct_diff = calcCtFromPt(Sa, pt) - Ct;
      pt_old  = pt;
      pt      = pt_old - ct_diff * dpt_dct;
      ptm     = 0.5 * (pt + pt_old);

      dpt_dct = -Cp0Sw / ((ptm + TkFrz) * calcGibbsDerivPt0Pt0(Sa, ptm));

      pt      = pt_old - ct_diff * dpt_dct;
      ct_diff = calcCtFromPt(Sa, pt) - Ct;
      pt_old  = pt;
      return (pt_old - ct_diff * dpt_dct);
   }

   /// Calculate Conservative Temperature from Potential Temp
   /// GSW Toolbox function gsw_ct_from_pt
   KOKKOS_FUNCTION Real calcCtFromPt(Real Sa, Real Pt) const {
      Real x2, x, y, pot_enthalpy;

      x2 = Sfac * Sa;
      x  = Kokkos::sqrt(x2);
      y  = Pt * 0.025e0; /*! normalize for F03 and F08 */
      pot_enthalpy =
          61.01362420681071e0 +
          y * (168776.46138048015e0 +
               y * (-2735.2785605119625e0 +
                    y * (2574.2164453821433e0 +
                         y * (-1536.6644434977543e0 +
                              y * (545.7340497931629e0 +
                                   (-50.91091728474331e0 -
                                    18.30489878927802e0 * y) *
                                       y))))) +
          x2 *
              (268.5520265845071e0 +
               y * (-12019.028203559312e0 +
                    y * (3734.858026725145e0 +
                         y * (-2046.7671145057618e0 +
                              y * (465.28655623826234e0 +
                                   (-0.6370820302376359e0 -
                                    10.650848542359153e0 * y) *
                                       y)))) +
               x * (937.2099110620707e0 +
                    y * (588.1802812170108e0 + y * (248.39476522971285e0 +
                                                    (-3.871557904936333e0 -
                                                     2.6268019854268356e0 * y) *
                                                        y)) +
                    x * (-1687.914374187449e0 +
                         x * (246.9598888781377e0 +
                              x * (123.59576582457964e0 -
                                   48.5891069025409e0 * x)) +
                         y * (936.3206544460336e0 +
                              y * (-942.7827304544439e0 +
                                   y * (369.4389437509002e0 +
                                        (-33.83664947895248e0 -
                                         9.987880382780322e0 * y) *
                                            y))))));

      return (pot_enthalpy / Cp0Sw);
   }

   /// Calculates freezing Conservative Temperature using TEOS-10 polynomial
   /// (polynomial error in [-5e-4, 6e-4] K, from GSW package)
   /// GSW Toolbox function gsw_ct_freezing_poly
   KOKKOS_FUNCTION Real calcCtFreezing(const Real Sa, const Real P,
                                       const Real SaturationFract) const {
      constexpr Real C0  = 0.017947064327968736;
      constexpr Real C1  = -6.076099099929818;
      constexpr Real C2  = 4.883198653547851;
      constexpr Real C3  = -11.88081601230542;
      constexpr Real C4  = 13.34658511480257;
      constexpr Real C5  = -8.722761043208607;
      constexpr Real C6  = 2.082038908808201;
      constexpr Real C7  = -7.389420998107497;
      constexpr Real C8  = -2.110913185058476;
      constexpr Real C9  = 0.2295491578006229;
      constexpr Real C10 = -0.9891538123307282;
      constexpr Real C11 = -0.08987150128406496;
      constexpr Real C12 = 0.3831132432071728;
      constexpr Real C13 = 1.054318231187074;
      constexpr Real C14 = 1.065556599652796;
      constexpr Real C15 = -0.7997496801694032;
      constexpr Real C16 = 0.3850133554097069;
      constexpr Real C17 = -2.078616693017569;
      constexpr Real C18 = 0.8756340772729538;
      constexpr Real C19 = -2.079022768390933;
      constexpr Real C20 = 1.596435439942262;
      constexpr Real C21 = 0.1338002171109174;
      constexpr Real C22 = 1.242891021876471;

      // Note: a = 0.502500117621 / SS0
      constexpr Real A = 0.014289763856964;
      constexpr Real B = 0.057000649899720;

      Real Sar = Sa * 1.0e-2;
      Real X   = Kokkos::sqrt(Sar);
      Real Pr  = P * 1.0e-4;

      Real CtFreez =
          C0 + Sar * (C1 + X * (C2 + X * (C3 + X * (C4 + X * (C5 + C6 * X))))) +
          Pr * (C7 + Pr * (C8 + C9 * Pr)) +
          Sar * Pr *
              (C10 + Pr * (C12 + Pr * (C15 + C21 * Sar)) +
               Sar * (C13 + C17 * Pr + C19 * Sar) +
               X * (C11 + Pr * (C14 + C18 * Pr) +
                    Sar * (C16 + C20 * Pr + C22 * Sar)));

      /* Adjust for the effects of dissolved air */
      CtFreez = CtFreez - SaturationFract * (1e-3) * (2.4 - A * Sa) *
                              (1.0 + B * (1.0 - Sa / SS0));

      return CtFreez;
   }

   /// Calculates the in-situ temperature at which seawater freezes
   /// GSW Toolbox function gsw_t_freezing_poly
   KOKKOS_FUNCTION Real calcISTempFreezing(const Real Salinity,
                                           const Real Pressure,
                                           const Real SaturationFract) const {

      constexpr Real T0  = 0.002519;
      constexpr Real T1  = -5.946302841607319;
      constexpr Real T2  = 4.136051661346983;
      constexpr Real T3  = -1.115150523403847e1;
      constexpr Real T4  = 1.476878746184548e1;
      constexpr Real T5  = -1.088873263630961e1;
      constexpr Real T6  = 2.961018839640730;
      constexpr Real T7  = -7.433320943962606;
      constexpr Real T8  = -1.561578562479883;
      constexpr Real T9  = 4.073774363480365e-2;
      constexpr Real T10 = 1.158414435887717e-2;
      constexpr Real T11 = -4.122639292422863e-1;
      constexpr Real T12 = -1.123186915628260e-1;
      constexpr Real T13 = 5.715012685553502e-1;
      constexpr Real T14 = 2.021682115652684e-1;
      constexpr Real T15 = 4.140574258089767e-2;
      constexpr Real T16 = -6.034228641903586e-1;
      constexpr Real T17 = -1.205825928146808e-2;
      constexpr Real T18 = -2.812172968619369e-1;
      constexpr Real T19 = 1.877244474023750e-2;
      constexpr Real T20 = -1.204395563789007e-1;
      constexpr Real T21 = 2.349147739749606e-1;
      constexpr Real T22 = 2.748444541144219e-3;

      const Real SaR = Salinity * 1.0e-2;
      const Real X   = Kokkos::sqrt(SaR);
      const Real Pr  = Pressure * 1.0e-4;

      Real TFreez =
          T0 + SaR * (T1 + X * (T2 + X * (T3 + X * (T4 + X * (T5 + T6 * X))))) +
          Pr * (T7 + Pr * (T8 + T9 * Pr)) +
          SaR * Pr * (T10 + Pr * (T12 + Pr * (T15 + T21 * SaR))) +
          SaR * (T13 + T17 * Pr + T19 * SaR) +
          X * (T11 + Pr * (T14 + T18 * Pr) +
               SaR * (T16 + T20 * Pr + T22 * SaR));

      TFreez = TFreez - SaturationFract * 1e-3_Real *
                            (2.4_Real - Salinity / (2.0_Real * SS0));

      return TFreez;
   }

   /// Calculates potential enthalpy of ice from potential temperature of ice
   /// using TEOS-10 polynomial The error of this fit ranges between -6e-3 and
   /// 6e-3 J/kg over the potential temperature range of -100 to 2 degC GSW
   /// Toolbox function gsw_pot_enthalpy_from_pt_ice_poly
   KOKKOS_FUNCTION Real
   calcPotEnthalpyIceFromPotTempIce(const Real PotTempIce) const {
      constexpr Real P0 = -3.333601570157700e5;
      constexpr Real P1 = 2.096693916810367e3;
      constexpr Real P2 = 3.687110754043292;
      constexpr Real P3 = 4.559401565980682e-4;
      constexpr Real P4 = -2.516011957758120e-6;
      constexpr Real P5 = -1.040364574632784e-8;
      constexpr Real P6 = -1.701786588412454e-10;
      constexpr Real P7 = -7.667191301635057e-13;

      // Initial estimate of the potential enthalpy
      Real PotEnthalpyIce =
          P0 +
          PotTempIce *
              (P1 +
               PotTempIce *
                   (P2 +
                    PotTempIce *
                        (P3 +
                         PotTempIce *
                             (P4 + PotTempIce *
                                       (P5 + PotTempIce *
                                                 (P6 + P7 * PotTempIce))))));

      Real DPotTempDPotEnth =
          calcPotTempIceFromPotEnthalpyIceDeriv(PotEnthalpyIce);

      Real PotEnthalpyIceOld, PotEnthalpyIceMid, F;

      for (int I = 0; I < 5; ++I) {
         PotEnthalpyIceOld = PotEnthalpyIce;
         F = calcPotTempIceFromPotEnthalpyIce(PotEnthalpyIceOld) - PotTempIce;
         PotEnthalpyIce    = PotEnthalpyIceOld - F / DPotTempDPotEnth;
         PotEnthalpyIceMid = 0.5_Real * (PotEnthalpyIce + PotEnthalpyIceOld);
         DPotTempDPotEnth =
             calcPotTempIceFromPotEnthalpyIceDeriv(PotEnthalpyIceMid);
         PotEnthalpyIce = PotEnthalpyIceOld - F / DPotTempDPotEnth;
      }

      return PotEnthalpyIce;
   }

   /// Calculate potential temperature of ice from the potential enthalpy of ice
   /// The error of this fit ranges between -5e-5 and 2e-4 degC over the
   /// potential temperature range of -100 to 2 degC. GSW Toolbox function
   /// gsw_pt_from_pot_enthalpy_ice_poly
   KOKKOS_FUNCTION Real
   calcPotTempIceFromPotEnthalpyIce(const Real PotEnthalpyIce) const {
      constexpr Real Q0 = 2.533588268773218e2;
      constexpr Real Q1 = 2.594351081876611e-3;
      constexpr Real Q2 = 1.765077810213815e-8;
      constexpr Real Q3 = 7.768070564290540e-14;
      constexpr Real Q4 = 2.034842254277530e-19;
      constexpr Real Q5 = 3.220014531712841e-25;
      constexpr Real Q6 = 2.845172809636068e-31;
      constexpr Real Q7 = 1.094005878892950e-37;

      Real IcePotTemp =
          Q0 +
          PotEnthalpyIce *
              (Q1 +
               PotEnthalpyIce *
                   (Q2 +
                    PotEnthalpyIce *
                        (Q3 +
                         PotEnthalpyIce *
                             (Q4 +
                              PotEnthalpyIce *
                                  (Q5 + PotEnthalpyIce *
                                            (Q6 + Q7 * PotEnthalpyIce))))));

      return IcePotTemp;
   }

   /// Calculate the derivative of potential temperature of ice with respect to
   /// potential enthalpy of ice GSW Toolbox function
   /// gsw_pt_from_pot_enthalpy_ice_poly_derivatives
   KOKKOS_FUNCTION Real
   calcPotTempIceFromPotEnthalpyIceDeriv(const Real PotEnthalpyIce) const {
      constexpr Real Q1 = 2.594351081876611e-3;
      constexpr Real P2 = 3.530155620427630e-8;
      constexpr Real P3 = 2.330421169287162e-13;
      constexpr Real P4 = 8.139369017110120e-19;
      constexpr Real P5 = 1.610007265856420e-24;
      constexpr Real P6 = 1.707103685781641e-30;
      constexpr Real P7 = 7.658041152250651e-37;

      Real DPotTempDPotEnth =
          Q1 +
          PotEnthalpyIce *
              (P2 +
               PotEnthalpyIce *
                   (P3 + PotEnthalpyIce *
                             (P4 + PotEnthalpyIce *
                                       (P5 + PotEnthalpyIce *
                                                 (P6 + P7 * PotEnthalpyIce)))));

      return DPotTempDPotEnth;
   }

   /// Calculate the mass fraction of ice (mass of ice divided by mass of ice
   /// plus seawater), which results from given values of the bulk absolute
   /// salinity, bulk potential enthalpy, and pressure. The final equalibrium
   /// values of absolute salinity and conservative temperature of the
   /// interstitial seawater phase are also returned. This code assumes there is
   /// no dissolved air in the seawater. When the mass fraction is positive, the
   /// seawater-ice mixture is at thermodynamic equilibrium. It also returns 0
   /// for the mass fraction if the input bulk enthalpy is sufficiently large,
   /// meaning there is no ice present in the final state. GSW Toolbox function
   /// gsw_frazil_properties_potential_poly
   KOKKOS_FUNCTION void calcFrazilProperties(const Real BulkAbsSalinity,
                                             const Real BulkPotEnthalpy,
                                             const Real Pressure,
                                             Real &InterstitialAbsSalinity,
                                             Real &InterstitialConservTemp,
                                             Real &IceMassFraction) const {

      constexpr Real F01 = -9.041191886754806e-1;
      constexpr Real F02 = 4.169608567309818e-2;
      constexpr Real F03 = -9.325971761333677e-3;
      constexpr Real F04 = 4.699055851002199e-2;
      constexpr Real F05 = -3.086923404061666e-2;
      constexpr Real F06 = 1.057761186019000e-2;
      constexpr Real F07 = -7.349302346007727e-2;
      constexpr Real F08 = 1.444842576424337e-1;
      constexpr Real F09 = -1.408425967872030e-1;
      constexpr Real F10 = 1.070981398567760e-1;
      constexpr Real F11 = -1.768451760854797e-2;
      constexpr Real F12 = -4.013688314067293e-1;
      constexpr Real F13 = 7.209753205388577e-1;
      constexpr Real F14 = -1.807444462285120e-1;
      constexpr Real F15 = 1.362305015808993e-1;
      constexpr Real F16 = -9.500974920072897e-1;
      constexpr Real F17 = 1.192134856624248;
      constexpr Real F18 = -9.191161283559850e-2;
      constexpr Real F19 = -1.008594411490973;
      constexpr Real F20 = 8.020279271484482e-1;
      constexpr Real F21 = -3.930534388853466e-1;
      constexpr Real F22 = -2.026853316399942e-2;
      constexpr Real F23 = -2.722731069001690e-2;
      constexpr Real F24 = 5.032098120548072e-2;
      constexpr Real F25 = -2.354888890484222e-2;
      constexpr Real F26 = -2.454090179215001e-2;
      constexpr Real F27 = 4.125987229048937e-2;
      constexpr Real F28 = -3.533404753585094e-2;
      constexpr Real F29 = 3.766063025852511e-2;
      constexpr Real F30 = -3.358409746243470e-2;
      constexpr Real F31 = -2.242158862056258e-2;
      constexpr Real F32 = 2.102254738058931e-2;
      constexpr Real F33 = -3.048635435546108e-2;
      constexpr Real F34 = -1.996293091714222e-2;
      constexpr Real F35 = 2.577703068234217e-2;
      constexpr Real F36 = -1.292053030649309e-2;

      constexpr Real G01 = 3.332286683867741e5;
      constexpr Real G02 = 1.416532517833479e4;
      constexpr Real G03 = -1.021129089258645e4;
      constexpr Real G04 = 2.356370992641009e4;
      constexpr Real G05 = -8.483432350173174e3;
      constexpr Real G06 = 2.279927781684362e4;
      constexpr Real G07 = 1.506238790315354e4;
      constexpr Real G08 = 4.194030718568807e3;
      constexpr Real G09 = -3.146939594885272e5;
      constexpr Real G10 = -7.549939721380912e4;
      constexpr Real G11 = 2.790535212869292e6;
      constexpr Real G12 = 1.078851928118102e5;
      constexpr Real G13 = -1.062493860205067e7;
      constexpr Real G14 = 2.082909703458225e7;
      constexpr Real G15 = -2.046810820868635e7;
      constexpr Real G16 = 8.039606992745191e6;
      constexpr Real G17 = -2.023984705844567e4;
      constexpr Real G18 = 2.871769638352535e4;
      constexpr Real G19 = -1.444841553038544e4;
      constexpr Real G20 = 2.261532522236573e4;
      constexpr Real G21 = -2.090579366221046e4;
      constexpr Real G22 = -1.128417003723530e4;
      constexpr Real G23 = 3.222965226084112e3;
      constexpr Real G24 = -1.226388046175992e4;
      constexpr Real G25 = 1.506847628109789e4;
      constexpr Real G26 = -4.584670946447444e4;
      constexpr Real G27 = 1.596119496322347e4;
      constexpr Real G28 = -6.338852410446789e4;
      constexpr Real G29 = 8.951570926106525e4;

      Real SaturationFract = 0.0_Real;
      Real DCtDSa, CtFreezing, DFuncDPotEnthalpy, DFuncDPotEnthalpyMeanPoly,
          DPotEnthDSa;
      Real Func, Func0, PotEnthalpy, InterstitialAbsSalinityTmp, IceMassFracOld,
          IceMassFracTmp, X, Xa, Y, Z;
      I4 Iterations, MaxIterations;

      //  Finding Func0.  This is the value of the function, Func, that would
      // result in the output IceMassFraction being exactly zero.
      Func0 =
          BulkPotEnthalpy -
          Cp0Sw * calcCtFreezing(BulkAbsSalinity, Pressure, SaturationFract);

      //  Setting the three outputs for data points that have Func0
      //  non-negative. When Func0 is zero or positive then the final answer
      //  will contain no frazil ice.
      if (Func0 >= 0.0_Real) {
         InterstitialAbsSalinity = BulkAbsSalinity;
         InterstitialConservTemp = BulkPotEnthalpy / Cp0Sw;
         IceMassFraction         = 0.0_Real;
         return;
      }

      // Begin finding the solution for data points that have Func0 < 0, so that
      // the output will have a positive ice mass fraction.
      // Evaluate a polynomial for IceMassFracTmp in terms of BulkAbsSalinity,
      // Func0, and Pressure.
      X = BulkAbsSalinity * 1.0e-2;
      Y = Func0 / 3.0e5;
      Z = Pressure * 1.0e-4;

      IceMassFracTmp =
          Y * (F01 + X * (F02 + X * (F03 + X * (F04 + X * (F05 + F06 * X)))) +
               Y * (F07 + X * (F08 + X * (F09 + X * (F10 + F11 * X))) +
                    Y * (F12 + X * (F13 + X * (F14 + F15 * X)) +
                         Y * (F16 + X * (F17 + F18 * X) +
                              Y * (F19 + F20 * X + F21 * Y)))) +
               Z * (F22 + X * (F23 + X * (F24 + F25 * X)) +
                    Y * (X * (F26 + F27 * X) + Y * (F28 + F29 * X + F30 * Y)) +
                    Z * (F31 + X * (F32 + F33 * X) +
                         Y * (F34 + F35 * X + F36 * Y))));

      // The ice mass fraction out of this code is restricted to be less than
      // 0.9.
      IceMassFracTmp = Kokkos::min(IceMassFracTmp, 0.9_Real);

      // The initial guess at the absolute salinity of the interstitial seawater
      InterstitialAbsSalinityTmp =
          BulkAbsSalinity / (1.0_Real - IceMassFracTmp);

      // Doing a Newton step with a separate polynomial estimate of the mean
      // derivative DFuncDPotEnthalpyMeanPoly.
      CtFreezing =
          calcCtFreezing(InterstitialAbsSalinityTmp, Pressure, SaturationFract);
      PotEnthalpy =
          calcPotEnthalpyIceFreezing(InterstitialAbsSalinityTmp, Pressure);
      Func = BulkPotEnthalpy -
             (1.0_Real - IceMassFracTmp) * Cp0Sw * CtFreezing -
             IceMassFracTmp * PotEnthalpy;

      Xa = InterstitialAbsSalinityTmp * 1.0e-2;

      DFuncDPotEnthalpyMeanPoly =
          G01 + Xa * (G02 + Xa * (G03 + Xa * (G04 + G05 * Xa))) +
          IceMassFracTmp *
              (Xa * (G06 + Xa * (G07 + G08 * Xa)) +
               IceMassFracTmp *
                   (Xa * (G09 + G10 * Xa) +
                    IceMassFracTmp * Xa *
                        (G11 + G12 * Xa +
                         IceMassFracTmp *
                             (G13 +
                              IceMassFracTmp *
                                  (G14 + IceMassFracTmp *
                                             (G15 + G16 * IceMassFracTmp)))))) +
          Z * (G17 + Xa * (G18 + G19 * Xa) +
               IceMassFracTmp *
                   (G20 + IceMassFracTmp * (G21 + G22 * IceMassFracTmp) +
                    Xa * (G23 + G24 * Xa * IceMassFracTmp)) +
               Z * (G25 + Xa * (G26 + G27 * Xa) +
                    IceMassFracTmp * (G28 + G29 * IceMassFracTmp)));

      IceMassFracOld = IceMassFracTmp;
      IceMassFracTmp = IceMassFracOld - Func / DFuncDPotEnthalpyMeanPoly;
      InterstitialAbsSalinityTmp =
          BulkAbsSalinity / (1.0_Real - IceMassFracTmp);

      // Calculating the estimate of the derivative of Func to be fed into
      // Newton's Method.
      CtFreezing =
          calcCtFreezing(InterstitialAbsSalinityTmp, Pressure, SaturationFract);
      PotEnthalpy =
          calcPotEnthalpyIceFreezing(InterstitialAbsSalinityTmp, Pressure);

      DCtDSa = calcConsTempFreezingFirstDerivPoly(InterstitialAbsSalinityTmp,
                                                  Pressure, SaturationFract);
      DPotEnthDSa = calcPotEnthalpyIceFirstDerivPoly(InterstitialAbsSalinityTmp,
                                                     Pressure);

      DFuncDPotEnthalpy =
          Cp0Sw * CtFreezing - PotEnthalpy -
          InterstitialAbsSalinityTmp *
              (Cp0Sw * DCtDSa +
               IceMassFracTmp * DPotEnthDSa / (1.0_Real - IceMassFracTmp));

      if (IceMassFracTmp >= 0.0_Real && IceMassFracTmp <= 0.2_Real &&
          InterstitialAbsSalinityTmp > 15.0_Real &&
          InterstitialAbsSalinityTmp < 60.0_Real && Pressure <= 3000.0_Real) {
         MaxIterations = 1;
      } else if (IceMassFracTmp >= 0.0_Real && IceMassFracTmp <= 0.85_Real &&
                 InterstitialAbsSalinityTmp > 0.0_Real &&
                 InterstitialAbsSalinityTmp < 120.0_Real &&
                 Pressure <= 3500.0_Real) {
         MaxIterations = 2;
      } else {
         MaxIterations = 3;
      }

      for (Iterations = 0; Iterations < MaxIterations; ++Iterations) {
         if (Iterations > 1) {
            // On the first iteration CtFreezing and PotEnthalpy are both known
            CtFreezing  = calcCtFreezing(InterstitialAbsSalinityTmp, Pressure,
                                         SaturationFract);
            PotEnthalpy = calcPotEnthalpyIceFreezing(InterstitialAbsSalinityTmp,
                                                     Pressure);
         }

         // This is the function, Fun, whose zero we seek
         Func = BulkPotEnthalpy -
                (1.0_Real - IceMassFracTmp) * Cp0Sw * CtFreezing -
                IceMassFracTmp * PotEnthalpy;
         IceMassFracOld = IceMassFracTmp;
         IceMassFracTmp = IceMassFracOld - Func / DFuncDPotEnthalpy;

         // The ice mass fraction out of this code is restricted to be less than
         // 0.9.
         IceMassFracTmp = Kokkos::min(IceMassFracTmp, 0.9_Real);

         InterstitialAbsSalinityTmp =
             BulkAbsSalinity / (1.0_Real - IceMassFracTmp);
      }

      if (IceMassFracTmp < 0.0_Real) {
         InterstitialAbsSalinity = BulkAbsSalinity;
         InterstitialConservTemp = BulkPotEnthalpy / Cp0Sw;
         IceMassFraction         = 0.0_Real;
      } else {
         InterstitialAbsSalinity = InterstitialAbsSalinityTmp;
         InterstitialConservTemp = calcCtFreezing(InterstitialAbsSalinityTmp,
                                                  Pressure, SaturationFract);
         IceMassFraction         = IceMassFracTmp;
      }
   }

   /// Calculates the first derivatives of the conservative temperature at which
   /// seawater freezes, with respect to absolute salinity GSW Toolbox function
   /// gsw_ct_freezing_first_derivatives_poly
   KOKKOS_FUNCTION Real
   calcConsTempFreezingFirstDerivPoly(const Real Salinity, const Real Pressure,
                                      const Real SaturationFract) const {
      constexpr Real A   = 0.014289763856964;
      constexpr Real B   = 0.057000649899720;
      constexpr Real C0  = 0.017947064327968736;
      constexpr Real C1  = -6.076099099929818;
      constexpr Real C2  = 4.883198653547851;
      constexpr Real C3  = -11.88081601230542;
      constexpr Real C4  = 13.34658511480257;
      constexpr Real C5  = -8.722761043208607;
      constexpr Real C6  = 2.082038908808201;
      constexpr Real C10 = -0.9891538123307282;
      constexpr Real C11 = -0.08987150128406496;
      constexpr Real C12 = 0.3831132432071728;
      constexpr Real C13 = 1.054318231187074;
      constexpr Real C14 = 1.065556599652796;
      constexpr Real C15 = -0.7997496801694032;
      constexpr Real C16 = 0.3850133554097069;
      constexpr Real C17 = -2.078616693017569;
      constexpr Real C18 = 0.8756340772729538;
      constexpr Real C19 = -2.079022768390933;
      constexpr Real C20 = 1.596435439942262;
      constexpr Real C21 = 0.1338002171109174;
      constexpr Real C22 = 1.242891021876471;

      Real D = -A - A * B - 2.4_Real * B / SS0;
      Real E = 2.0_Real * A * B / SS0;

      const Real SaR = Salinity * 1.0e-2;
      const Real X   = Kokkos::sqrt(SaR);
      const Real Pr  = Pressure * 1.0e-4;

      return (C1 +
              X * (1.5_Real * C2 +
                   X * (2.0_Real * C3 +
                        X * (2.5_Real * C4 +
                             X * (3.0_Real * C5 + 3.5_Real * C6 * X)))) +
              Pr *
                  (C10 +
                   X * (1.5_Real * C11 +
                        X * (2.0_Real * C13 +
                             X * (2.5_Real * C16 +
                                  X * (3.0_Real * C19 + 3.5_Real * C22 * X)))) +
                   Pr * (C12 +
                         X * (1.5_Real * C14 +
                              X * (2.0_Real * C17 + 2.5_Real * C20 * X)) +
                         Pr * (C15 +
                               X * (1.5_Real * C18 + 2.0_Real * C21 * X))))) *
                 1e-2_Real -
             SaturationFract * 1e-3_Real * (D - Salinity * E);
   }

   /// Calculate the first derivatives of the potential enthalpy of ice at which
   /// ice melts into seawater with respect to absolute salinity GSW Toolbox
   /// function gsw_pot_enthalpy_ice_freezing_first_derivatives_poly
   KOKKOS_FUNCTION Real calcPotEnthalpyIceFirstDerivPoly(
       const Real Salinity, const Real Pressure) const {
      constexpr Real D1  = -1.249490228128056e4;
      constexpr Real D2  = 1.336783910789822e4;
      constexpr Real D3  = -4.811989517774642e4;
      constexpr Real D4  = 8.044864276240987e4;
      constexpr Real D5  = -7.124452125071862e4;
      constexpr Real D6  = 2.280706828014839e4;
      constexpr Real D7  = 0.315423710959628e3;
      constexpr Real D8  = -3.592775732074710e2;
      constexpr Real D9  = 1.644828513129230e3;
      constexpr Real D10 = -4.809640968940840e3;
      constexpr Real D11 = 2.901071777977272e3;
      constexpr Real D12 = -9.218459682855746e2;
      constexpr Real D13 = 0.379377450285737e3;
      constexpr Real D14 = -2.672164989849465e3;
      constexpr Real D15 = 5.044317489422632e3;
      constexpr Real D16 = -2.631711865886377e3;
      constexpr Real D17 = -0.160245473297112e3;
      constexpr Real D18 = 4.029061696035465e2;
      constexpr Real D19 = -3.682950019675760e2;

      constexpr Real F1  = -2.034535061416256e4;
      constexpr Real F2  = 0.315423710959628e3;
      constexpr Real F3  = -0.239518382138314e3;
      constexpr Real F4  = 0.822414256564615e3;
      constexpr Real F5  = -1.923856387576336e3;
      constexpr Real F6  = 0.967023925992424e3;
      constexpr Real F7  = -0.263384562367307e3;
      constexpr Real F8  = -5.051613740291480e3;
      constexpr Real F9  = 7.587549005714740e2;
      constexpr Real F10 = -3.562886653132620e3;
      constexpr Real F11 = 5.044317489422632e3;
      constexpr Real F12 = -2.105369492709102e3;
      constexpr Real F13 = 6.387082316647800e2;
      constexpr Real F14 = -4.807364198913360e2;
      constexpr Real F15 = 8.058123392070929e2;
      constexpr Real F16 = -5.524425029513641e2;

      const Real SaR = Salinity * 1.0e-2;
      const Real X   = Kokkos::sqrt(SaR);
      const Real Pr  = Pressure * 1.0e-4;

      return (D1 + X * (D2 + X * (D3 + X * (D4 + X * (D5 + D6 * X)))) +
              Pr * (D7 + X * (D8 + X * (D9 + X * (D10 + X * (D11 + D12 * X)))) +
                    Pr * (D13 + X * (D14 + X * (D15 + D16 * X)) +
                          Pr * (D17 + X * (D18 + D19 * X))))) *
             1e-2;
   }

   /// Calculates the potential enthalpy of ice at which seawater freezes.
   /// The error of this fit ranges between -2.5 and 1 J/kg with an rms of
   /// 1.07, between SA of 0 and 120 g/kg and Pressure between 0 and 10,000 dbar
   /// The error of the fit is between -0.7 and 0.7 with and rums of 0.3
   /// between SA of 0 and 120 g/kg and Pressure between 0 and 5,000 dbar.
   /// GSW Toolbox function
   KOKKOS_FUNCTION Real calcPotEnthalpyIceFreezing(const Real Salinity,
                                                   const Real Pressure) const {
      constexpr Real C0  = -3.333548730778702e5;
      constexpr Real C1  = -1.249490228128056e4;
      constexpr Real C2  = 0.891189273859881e4;
      constexpr Real C3  = -2.405994758887321e4;
      constexpr Real C4  = 3.217945710496395e4;
      constexpr Real C5  = -2.374817375023954e4;
      constexpr Real C6  = 0.651630522289954e4;
      constexpr Real C7  = -2.034535061416256e4;
      constexpr Real C8  = -0.252580687014574e4;
      constexpr Real C9  = 0.021290274388826e4;
      constexpr Real C10 = 0.315423710959628e3;
      constexpr Real C11 = -0.239518382138314e3;
      constexpr Real C12 = 0.379377450285737e3;
      constexpr Real C13 = 0.822414256564615e3;
      constexpr Real C14 = -1.781443326566310e3;
      constexpr Real C15 = -0.160245473297112e3;
      constexpr Real C16 = -1.923856387576336e3;
      constexpr Real C17 = 2.522158744711316e3;
      constexpr Real C18 = 0.268604113069031e3;
      constexpr Real C19 = 0.967023925992424e3;
      constexpr Real C20 = -1.052684746354551e3;
      constexpr Real C21 = -0.184147500983788e3;
      constexpr Real C22 = -0.263384562367307e3;

      const Real SaR = Salinity * 1.0e-2;
      const Real X   = Kokkos::sqrt(SaR);
      const Real Pr  = Pressure * 1.0e-4;

      return C0 +
             SaR * (C1 + X * (C2 + X * (C3 + X * (C4 + X * (C5 + C6 * X))))) +
             Pr * (C7 + Pr * (C8 + C9 * Pr)) +
             SaR * Pr *
                 (C10 + Pr * (C12 + Pr * (C15 + C21 * SaR)) +
                  SaR * (C13 + C17 * Pr + C19 * SaR) +
                  X * (C11 + Pr * (C14 + C18 * Pr) +
                       SaR * (C16 + C20 * Pr + C22 * SaR)));
   }

   /// Calculates the final absolute salinity, final conservative temperature,
   /// and final ice mass fraction that results when a given mass fraction of
   /// ice melts and is mized into seawater whose properties are (SA, CT, P).
   /// When the mass fraction is calculated as being a positive value, the
   /// seawater-ice mixture is at thermodynamic equilibrium. It also returns 0
   /// for the mass fraction if the input bulk enthalpy is sufficiently large,
   /// meaning there is no ice present in the final state. GSW Toolbox function
   /// gsw_melting_ice_into_seawater
   KOKKOS_FUNCTION void calcMeltingIceIntoSeawater(
       const Real AbsSalinity, const Real ConservTemp, const Real Pressure,
       const Real IceMassFractionInput, const Real IceTemperature,
       Real &AbsSalinityFinal, Real &ConservTempFinal,
       Real &IceMassFractionFinal) const {

      Real SaturationFract = 0.0_Real;
      // Note: GSW Toolbox uses gsw_ct_freezing_exact() here, which uses a
      // Newton-Raphson iteration to find the freezing temperature. This is not
      // necessary for our purposes, so we use the polynomial fit instead,
      // gsw_ct_freezing_poly() [our calcCtFreezing function above]. The error
      // of this fit ranges between -5e-4 K and 6e-4 K compared to the exact
      // freezing temperature.
      // Real CtFreezingWater = calcCtFreezing(AbsSalinity, Pressure,
      // SaturationFract); if (ConservTemp < CtFreezingWater) {
      // The sea water Conservative temperature is below the freezing point
      // Error handling: sea water Conservative temperature is below the
      // freezing point
      //   Error;
      //}

      // Note: GSW Toolbox uses gsw_t_freezing_exact() here, which uses a
      // Newton-Raphson iteration to find the freezing temperature. This is not
      // necessary for our purposes, so we use the polynomial fit instead,
      // gsw_t_freezing_poly() [our calcISTempFreezing function above].
      // Real TFreezingIce = calcISTempFreezing(0.0_Real, Pressure,
      // SaturationFract) - 1e-6_Real; if (IceTemperature > TFreezingIce) {
      // Error handling: ice temperature is above the freezing point
      //   Error;
      //}

      Real BulkAbsSalinity = (1.0_Real - IceMassFractionInput) * AbsSalinity;

      // Note: GSW Toolbox uses gsw_enthalpy_ct_exact() here, which uses a
      // Newton-Raphson iteration to find the enthalpy. This is not necessary
      // for our purposes, so we use the polynomial fit instead, gsw_enthalpy()
      // [our calcSpecEnthalpySW function below].
      Real BulkPotEnthalpy =
          (1.0_Real - IceMassFractionInput) *
              calcSpecEnthalpySW(AbsSalinity, ConservTemp, Pressure) +
          IceMassFractionInput * calcSpecEnthalpyIce(IceTemperature, Pressure);

      calcFrazilProperties(BulkAbsSalinity, BulkPotEnthalpy, Pressure,
                           AbsSalinityFinal, ConservTempFinal,
                           IceMassFractionFinal);
   }

   /// Calculate the specific enthalpy of ice
   /// GSW Toolbox function gsw_enthalpy_ice
   KOKKOS_FUNCTION Real calcSpecEnthalpyIce(const Real ISTemperature,
                                            const Real Pressure) const {
      constexpr Real G00 = -6.32020233335886e5;
      constexpr Real G01 = 6.55022213658955e-1;
      constexpr Real G02 = -1.89369929326131e-8;
      constexpr Real G03 = 3.3974612327105304e-15;
      constexpr Real G04 = -5.564648690589909e-22;
      const Kokkos::complex<Real> T1(3.68017112855051e-2_Real,
                                     5.10878114959572e-2_Real);
      const Kokkos::complex<Real> T2(3.37315741065416e-1_Real,
                                     3.35449415919309e-1_Real);
      const Kokkos::complex<Real> R1(4.47050716285388e1_Real,
                                     6.56876847463481e1_Real);
      const Kokkos::complex<Real> R20(-7.25974574329220e1_Real,
                                      -7.81008427112870e1_Real);
      const Kokkos::complex<Real> R21(-5.57107698030123e-5_Real,
                                      4.64578634580806e-5_Real);
      const Kokkos::complex<Real> R22(2.34801409215913e-11_Real,
                                      -2.85651142904972e-11_Real);

      Real Tau = (ISTemperature + TkFrz) / TkTrip;
      Real DZI = Db2Pa * Pressure / 611.657; // PTrip;
      Real G0  = G00 + DZI * (G01 + DZI * (G02 + DZI * (G03 + G04 * DZI)));
      Kokkos::complex<Real> R2      = R20 + DZI * (R21 + R22 * DZI);
      Kokkos::complex<Real> SqTauT1 = Kokkos::pow(Tau / T1, 2.0_Real);
      Kokkos::complex<Real> SqTauT2 = Kokkos::pow(Tau / T2, 2.0_Real);
      Kokkos::complex<Real> G =
          R1 * T1 * (Kokkos::log(1.0_Real - SqTauT1) + SqTauT2) +
          R2 * T2 * (Kokkos::log(1.0_Real - SqTauT2) + SqTauT1);

      return G0 + TkTrip * G.real(); // Kokkos::abs(G);
   }

   /// Calculate the specific enthalpy of seawater from absolute salinity and
   /// conservative temperature and pressure GSW Toolbox function  gsw_enthalpy
   /// and gsw_dynamic_enthalpy
   KOKKOS_FUNCTION Real calcSpecEnthalpySW(const Real AbsSalinity,
                                           const Real ConservTemp,
                                           const Real Pressure) const {
      constexpr Real H001 = 1.07699958620e-3;
      constexpr Real H002 = -3.03995719050e-5;
      constexpr Real H003 = 3.32853897400e-6;
      constexpr Real H004 = -2.82734035930e-7;
      constexpr Real H005 = 2.10623061600e-8;
      constexpr Real H006 = -2.10787688100e-9;
      constexpr Real H007 = 2.80192913290e-10;
      constexpr Real H011 = -1.56497346750e-5;
      constexpr Real H012 = 9.25288271450e-6;
      constexpr Real H013 = -3.91212891030e-7;
      constexpr Real H014 = -9.13175163830e-8;
      constexpr Real H015 = 6.29081998040e-8;
      constexpr Real H021 = 2.77621064840e-5;
      constexpr Real H022 = -5.85830342650e-6;
      constexpr Real H023 = 7.10167624670e-7;
      constexpr Real H024 = 7.17397628980e-8;
      constexpr Real H031 = -1.65211592590e-5;
      constexpr Real H032 = 3.96398280870e-6;
      constexpr Real H033 = -1.53775133460e-7;
      constexpr Real H042 = -1.70510937410e-6;
      constexpr Real H043 = -2.11176388380e-8;
      constexpr Real H041 = 6.91113227020e-6;
      constexpr Real H051 = -8.05396155400e-7;
      constexpr Real H052 = 2.53683834070e-7;
      constexpr Real H061 = 2.05430942680e-7;
      constexpr Real H101 = -3.10389819760e-4;
      constexpr Real H102 = 1.21312343735e-5;
      constexpr Real H103 = -1.94948109950e-7;
      constexpr Real H104 = 9.07754712880e-8;
      constexpr Real H105 = -2.22942508460e-8;
      constexpr Real H111 = 3.50095997640e-5;
      constexpr Real H112 = -4.78385440780e-6;
      constexpr Real H113 = -1.85663848520e-6;
      constexpr Real H114 = -6.82392405930e-8;
      constexpr Real H121 = -3.74358423440e-5;
      constexpr Real H122 = -1.18391541805e-7;
      constexpr Real H123 = 1.30457956930e-7;
      constexpr Real H131 = 2.41414794830e-5;
      constexpr Real H132 = -1.72793868275e-6;
      constexpr Real H133 = 2.58729626970e-9;
      constexpr Real H141 = -8.75958731540e-6;
      constexpr Real H142 = 6.47835889150e-7;
      constexpr Real H151 = -3.30527589000e-7;
      constexpr Real H201 = 6.69280670380e-4;
      constexpr Real H202 = -1.73962304870e-5;
      constexpr Real H203 = -1.60407505320e-6;
      constexpr Real H204 = 4.18657594500e-9;
      constexpr Real H211 = -4.35926785610e-5;
      constexpr Real H212 = 5.55041738250e-6;
      constexpr Real H213 = 1.82069162780e-6;
      constexpr Real H221 = 3.59078227600e-5;
      constexpr Real H222 = 1.46416731475e-6;
      constexpr Real H223 = -2.19103680220e-7;
      constexpr Real H231 = -1.43536330480e-5;
      constexpr Real H232 = 1.58276530390e-7;
      constexpr Real H241 = 4.37036805980e-6;
      constexpr Real H301 = -8.50479339370e-4;
      constexpr Real H302 = 1.87353886525e-5;
      constexpr Real H303 = 1.64210356660e-6;
      constexpr Real H311 = 3.45324618280e-5;
      constexpr Real H312 = -4.92235589220e-6;
      constexpr Real H313 = -4.51472854230e-7;
      constexpr Real H321 = -1.86985841870e-5;
      constexpr Real H322 = -2.44130696000e-7;
      constexpr Real H331 = 2.28633245560e-6;
      constexpr Real H401 = 5.80860699430e-4;
      constexpr Real H402 = -8.66110930600e-6;
      constexpr Real H403 = -5.93732490900e-7;
      constexpr Real H411 = -1.19594097880e-5;
      constexpr Real H421 = 3.85953392440e-6;
      constexpr Real H412 = 1.29546126300e-6;
      constexpr Real H501 = -2.10923705070e-4;
      constexpr Real H502 = 1.54637136265e-6;
      constexpr Real H511 = 1.38645945810e-6;
      constexpr Real H601 = 3.19324573050e-5;

      Real SFac   = 0.0248826675584615_Real;
      Real Offset = 5.971840214030754e-1_Real;
      Real Xs     = Kokkos::sqrt(AbsSalinity * SFac + Offset);
      Real Ys     = ConservTemp * 0.025_Real;
      Real Z      = Pressure * 1.0e-4_Real;

      return (Z *
              (H001 +
               Xs * (H101 +
                     Xs * (H201 +
                           Xs * (H301 +
                                 Xs * (H401 + Xs * (H501 + H601 * Xs))))) +
               Ys *
                   (H011 +
                    Xs * (H111 +
                          Xs * (H211 + Xs * (H311 + Xs * (H411 + H511 * Xs)))) +
                    Ys * (H021 +
                          Xs * (H121 + Xs * (H221 + Xs * (H321 + H421 * Xs))) +
                          Ys * (H031 + Xs * (H131 + Xs * (H231 + H331 * Xs)) +
                                Ys * (H041 + Xs * (H141 + H241 * Xs) +
                                      Ys * (H051 + H151 * Xs + H061 * Ys))))) +
               Z * (H002 +
                    Xs * (H102 +
                          Xs * (H202 + Xs * (H302 + Xs * (H402 + H502 * Xs)))) +
                    Ys * (H012 +
                          Xs * (H112 + Xs * (H212 + Xs * (H312 + H412 * Xs))) +
                          Ys * (H022 + Xs * (H122 + Xs * (H222 + H322 * Xs)) +
                                Ys * (H032 + Xs * (H132 + H232 * Xs) +
                                      Ys * (H042 + H142 * Xs + H052 * Ys)))) +
                    Z * (H003 +
                         Xs * (H103 + Xs * (H203 + Xs * (H303 + H403 * Xs))) +
                         Ys * (H013 + Xs * (H113 + Xs * (H213 + H313 * Xs)) +
                               Ys * (H023 + Xs * (H123 + H223 * Xs) +
                                     Ys * (H033 + H133 * Xs + H043 * Ys))) +
                         Z * (H004 + Xs * (H104 + H204 * Xs) +
                              Ys * (H014 + H114 * Xs + H024 * Ys) +
                              Z * (H005 + H105 * Xs + H015 * Ys +
                                   Z * (H006 + H007 * Z))))))) *
                 Db2Pa * 1.0e-4_Real +
             Cp0Sw * ConservTemp;
   }

   /// Calculate in-situ temperature from the potential temperature of ice with
   /// reference pressure of 0 dbar and the in-situ pressure GSW Toolbox
   /// function gsw_t_from_pt0_ice
   KOKKOS_FUNCTION Real
   calcISTempFromPotTempIce(const Real ISTemperature) const {
      constexpr Real P1 = -2.259745637898635e-4;
      constexpr Real P2 = 1.486236778150360e-9;
      constexpr Real P3 = 6.257869607978536e-12;
      constexpr Real P4 = -5.253795281359302e-7;
      constexpr Real P5 = 6.752596995671330e-9;
      constexpr Real P6 = 2.082992190070936e-11;

      constexpr Real Q1 = -5.849191185294459e-15;
      constexpr Real Q2 = 9.330347971181604e-11;
      constexpr Real Q3 = 3.415888886921213e-13;
      constexpr Real Q4 = 1.064901553161811e-12;
      constexpr Real Q5 = -1.454060359158787e-10;
      constexpr Real Q6 = -5.323461372791532e-13;

      ///// need to add calculation
      Real Pressure = 0.0_Real;
      Real Dp       = Pressure - PRef;
      Real PotTempIce =
          ISTemperature +
          Dp * (P1 + (Pressure + PRef) * (P2 + P3 * ISTemperature) +
                ISTemperature * (P5 + P6 * ISTemperature));

      PotTempIce = Kokkos::min(PotTempIce, -T0);
      PotTempIce = Kokkos::max(PotTempIce, -TkFrz + 0.15_Real);

      DEntropyDt = calcGibbsIce(2, 9, PotTempIce, PRef);

      EntropyTrue = calcGibbsIcePartT(ISTemperature, Pressure);

      ///.....

      return;
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DReal LatCell;
   Array1DReal LonCell;
};

/// Linear Equation of State
class LinearEos {
 public:
   /// Coefficients for LinearEos (overwritten by config file if set there)
   Real DRhodT  = -0.2;  ///< Thermal expansion coefficient (kg m^-3 degC^-1)
   Real DRhodS  = 0.8;   ///< Haline contraction coefficient (kg m^-3)
   Real RhoT0S0 = RhoFw; ///< Reference density (kg m^-3) at (T,S)=(0,0)

   /// constructor declaration
   LinearEos(const VertCoord *VCoord);

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   linear specific volume.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         SpecVol(ICell, K) =
             1.0_Real / (RhoT0S0 + (DRhodT * ConservTemp(ICell, K) +
                                    DRhodS * AbsSalinity(ICell, K)));
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Constant Equation of State
class ConstantEos {
 public:
   /// constructor declaration
   ConstantEos(const VertCoord *VCoord);

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and returns a constant specific volume
   //   value for all active layers.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));
      (void)ConservTemp;
      (void)AbsSalinity;

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K        = KStart + KVec;
         SpecVol(ICell, K) = 1.0_Real / RhoSw;
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Functor for calculating the squared Brunt-Vaisala frequency using TEOS-10
class Teos10BruntVaisalaFreqSq {
 public:
   /// Constructor for BruntVaisalaFreqSq
   Teos10BruntVaisalaFreqSq(const VertCoord *VCoord);

   //   The functor takes the full arrays of squared Brunt-Vaisala frequency
   //   (inout) the index ICell, and the ocean tracers (conservative)
   //   temperature, (absolute) salinity, pressure, specific volume as inputs,
   //   and outputs the squared Brunt-Vaisala frequency.
   KOKKOS_FUNCTION void operator()(Array2DReal BruntVaisalaFreqSq, I4 ICell,
                                   I4 KChunk, const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity,
                                   const Array2DReal &Pressure,
                                   const Array2DReal &SpecVol) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell) + 1);
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         // Calculate squared Brunt-Vaisala frequency
         Real CtInt =
             0.5_Real * (ConservTemp(ICell, K) + ConservTemp(ICell, K - 1));
         Real SaInt =
             0.5_Real * (AbsSalinity(ICell, K) + AbsSalinity(ICell, K - 1));
         Real PInt  = 0.5_Real * (Pressure(ICell, K) + Pressure(ICell, K - 1));
         Real SpInt = 0.5_Real * (SpecVol(ICell, K) + SpecVol(ICell, K - 1));
         Real AlphaInt = calcAlpha(SaInt, CtInt, PInt * Pa2Db, SpInt);
         Real BetaInt  = calcBeta(SaInt, CtInt, PInt * Pa2Db, SpInt);
         Real DSa      = AbsSalinity(ICell, K) - AbsSalinity(ICell, K - 1);
         Real DCt      = ConservTemp(ICell, K) - ConservTemp(ICell, K - 1);
         Real DP       = Pressure(ICell, K) - Pressure(ICell, K - 1);

         BruntVaisalaFreqSq(ICell, K) = Gravity * Gravity *
                                        (BetaInt * DSa - AlphaInt * DCt) /
                                        (SpInt * DP);
      }
   }

   /// Calculate alpha values for the squared Brunt-Vaisala frequency
   KOKKOS_FUNCTION Real calcAlpha(Real Sa, Real Ct, Real P, Real Sp) const {
      constexpr Real Factor = 0.0248826675584615;
      constexpr Real Offset = 5.971840214030754e-1;
      constexpr Real PNorm  = 1.0e-4;
      Real Ss               = Kokkos::sqrt(Factor * Sa + Offset);
      Real Tt               = 0.025_Real * Ct;
      Real Pp               = P * PNorm;

      constexpr Real A000 = -1.56497346750e-5;
      constexpr Real A001 = 1.85057654290e-5;
      constexpr Real A002 = -1.17363867310e-6;
      constexpr Real A003 = -3.65270065530e-7;
      constexpr Real A004 = 3.14540999020e-7;
      constexpr Real A010 = 5.55242129680e-5;
      constexpr Real A011 = -2.34332137060e-5;
      constexpr Real A012 = 4.26100574800e-6;
      constexpr Real A013 = 5.73918103180e-7;
      constexpr Real A020 = -4.95634777770e-5;
      constexpr Real A021 = 2.37838968519e-5;
      constexpr Real A022 = -1.38397620111e-6;
      constexpr Real A030 = 2.76445290808e-5;
      constexpr Real A031 = -1.36408749928e-5;
      constexpr Real A032 = -2.53411666056e-7;
      constexpr Real A040 = -4.02698077700e-6;
      constexpr Real A041 = 2.53683834070e-6;
      constexpr Real A050 = 1.23258565608e-6;
      constexpr Real A100 = 3.50095997640e-5;
      constexpr Real A101 = -9.56770881560e-6;
      constexpr Real A102 = -5.56991545570e-6;
      constexpr Real A103 = -2.72956962370e-7;
      constexpr Real A110 = -7.48716846880e-5;
      constexpr Real A111 = -4.73566167220e-7;
      constexpr Real A112 = 7.82747741600e-7;
      constexpr Real A120 = 7.24244384490e-5;
      constexpr Real A121 = -1.03676320965e-5;
      constexpr Real A122 = 2.32856664276e-8;
      constexpr Real A130 = -3.50383492616e-5;
      constexpr Real A131 = 5.18268711320e-6;
      constexpr Real A140 = -1.65263794500e-6;
      constexpr Real A200 = -4.35926785610e-5;
      constexpr Real A201 = 1.11008347650e-5;
      constexpr Real A202 = 5.46207488340e-6;
      constexpr Real A210 = 7.18156455200e-5;
      constexpr Real A211 = 5.85666925900e-6;
      constexpr Real A212 = -1.31462208134e-6;
      constexpr Real A220 = -4.30608991440e-5;
      constexpr Real A221 = 9.49659182340e-7;
      constexpr Real A230 = 1.74814722392e-5;
      constexpr Real A300 = 3.45324618280e-5;
      constexpr Real A301 = -9.84471178440e-6;
      constexpr Real A302 = -1.35441856270e-6;
      constexpr Real A310 = -3.73971683740e-5;
      constexpr Real A311 = -9.76522784000e-7;
      constexpr Real A320 = 6.85899736680e-6;
      constexpr Real A400 = -1.19594097880e-5;
      constexpr Real A401 = 2.59092252600e-6;
      constexpr Real A410 = 7.71906784880e-6;
      constexpr Real A500 = 1.38645945810e-6;

      Real Rval =
          A000 +
          Ss * (A100 + Ss * (A200 + Ss * (A300 + Ss * (A400 + A500 * Ss)))) +
          Tt * (A010 + Ss * (A110 + Ss * (A210 + Ss * (A310 + A410 * Ss))) +
                Tt * (A020 + Ss * (A120 + Ss * (A220 + A320 * Ss)) +
                      Tt * (A030 + Ss * (A130 + A230 * Ss) +
                            Tt * (A040 + A140 * Ss + A050 * Tt)))) +
          Pp * (A001 + Ss * (A101 + Ss * (A201 + Ss * (A301 + A401 * Ss))) +
                Tt * (A011 + Ss * (A111 + Ss * (A211 + A311 * Ss)) +
                      Tt * (A021 + Ss * (A121 + A221 * Ss) +
                            Tt * (A031 + A131 * Ss + A041 * Tt))) +
                Pp * (A002 + Ss * (A102 + Ss * (A202 + A302 * Ss)) +
                      Tt * (A012 + Ss * (A112 + A212 * Ss) +
                            Tt * (A022 + A122 * Ss + A032 * Tt)) +
                      Pp * (A003 + A103 * Ss + A013 * Tt + A004 * Pp)));

      return 0.025_Real * Rval / Sp;
   }

   /// Calculate beta values for the squared Brunt-Vaisala frequency
   KOKKOS_FUNCTION Real calcBeta(Real Sa, Real Ct, Real P, Real Sp) const {
      constexpr Real Factor = 0.0248826675584615;
      constexpr Real Offset = 5.971840214030754e-1;
      constexpr Real PNorm  = 1.0e-4;
      Real Ss               = Kokkos::sqrt(Factor * Sa + Offset);
      Real Tt               = 0.025_Real * Ct;
      Real Pp               = P * PNorm;

      constexpr Real B000 = -3.10389819760e-4;
      constexpr Real B003 = 3.63101885150e-7;
      constexpr Real B004 = -1.11471254230e-7;
      constexpr Real B010 = 3.50095997640e-5;
      constexpr Real B013 = -2.72956962370e-7;
      constexpr Real B020 = -3.74358423440e-5;
      constexpr Real B030 = 2.41414794830e-5;
      constexpr Real B040 = -8.75958731540e-6;
      constexpr Real B050 = -3.30527589000e-7;
      constexpr Real B100 = 1.33856134076e-3;
      constexpr Real B103 = 3.34926075600e-8;
      constexpr Real B110 = -8.71853571220e-5;
      constexpr Real B120 = 7.18156455200e-5;
      constexpr Real B130 = -2.87072660960e-5;
      constexpr Real B140 = 8.74073611960e-6;
      constexpr Real B200 = -2.55143801811e-3;
      constexpr Real B210 = 1.03597385484e-4;
      constexpr Real B220 = -5.60957525610e-5;
      constexpr Real B230 = 6.85899736680e-6;
      constexpr Real B300 = 2.32344279772e-3;
      constexpr Real B310 = -4.78376391520e-5;
      constexpr Real B320 = 1.54381356976e-5;
      constexpr Real B400 = -1.05461852535e-3;
      constexpr Real B410 = 6.93229729050e-6;
      constexpr Real B500 = 1.91594743830e-4;
      constexpr Real B001 = 2.42624687470e-5;
      constexpr Real B011 = -9.56770881560e-6;
      constexpr Real B021 = -2.36783083610e-7;
      constexpr Real B031 = -3.45587736550e-6;
      constexpr Real B041 = 1.29567177830e-6;
      constexpr Real B101 = -6.95849219480e-5;
      constexpr Real B111 = 2.22016695300e-5;
      constexpr Real B121 = 5.85666925900e-6;
      constexpr Real B131 = 6.33106121560e-7;
      constexpr Real B201 = 1.12412331915e-4;
      constexpr Real B211 = -2.95341353532e-5;
      constexpr Real B221 = -1.46478417600e-6;
      constexpr Real B301 = -6.92888744480e-5;
      constexpr Real B311 = 1.03636901040e-5;
      constexpr Real B401 = 1.54637136265e-5;
      constexpr Real B002 = -5.84844329840e-7;
      constexpr Real B012 = -5.56991545570e-6;
      constexpr Real B022 = 3.91373870800e-7;
      constexpr Real B032 = 7.76188880920e-9;
      constexpr Real B102 = -9.62445031940e-6;
      constexpr Real B112 = 1.09241497668e-5;
      constexpr Real B122 = -1.31462208134e-6;
      constexpr Real B202 = 1.47789320994e-5;
      constexpr Real B212 = -4.06325568810e-6;
      constexpr Real B302 = -7.12478989080e-6;

      Real Rval =
          B000 +
          Ss * (B100 + Ss * (B200 + Ss * (B300 + Ss * (B400 + B500 * Ss)))) +
          Tt * (B010 + Ss * (B110 + Ss * (B210 + Ss * (B310 + B410 * Ss))) +
                Tt * (B020 + Ss * (B120 + Ss * (B220 + B320 * Ss)) +
                      Tt * (B030 + Ss * (B130 + B230 * Ss) +
                            Tt * (B040 + B140 * Ss + B050 * Tt)))) +
          Pp * (B001 + Ss * (B101 + Ss * (B201 + Ss * (B301 + B401 * Ss))) +
                Tt * (B011 + Ss * (B111 + Ss * (B211 + B311 * Ss)) +
                      Tt * (B021 + Ss * (B121 + B221 * Ss) +
                            Tt * (B031 + B131 * Ss + B041 * Tt))) +
                Pp * (B002 + Ss * (B102 + Ss * (B202 + B302 * Ss)) +
                      Tt * (B012 + Ss * (B112 + B212 * Ss) +
                            Tt * (B022 + B122 * Ss + B032 * Tt)) +
                      Pp * (B003 + B103 * Ss + B013 * Tt + B004 * Pp)));

      return -0.5_Real * Rval * Factor / (Sp * Ss);
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Linear squared Brunt-Vaisala frequency calculator
class LinearBruntVaisalaFreqSq {
 public:
   /// constructor declaration
   LinearBruntVaisalaFreqSq(const VertCoord *VCoord);

   //   The functor takes the full arrays of squared Brunt-Vaisala frequency
   //   (inout), the index ICell, and the specific volume and pseudo-thickness
   //   as inputs, and outputs the squared Brunt-Vaisala frequency.
   KOKKOS_FUNCTION void operator()(Array2DReal BruntVaisalaFreqSq, I4 ICell,
                                   I4 KChunk,
                                   const Array2DReal &SpecVol) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell) + 1);
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         /// Calculate squared Brunt-Vaisala frequency at mid-point between
         /// K-1 and K Do not need to use displaced specific volume here
         /// since only the linear EOS is used with this BVF formulation.
         BruntVaisalaFreqSq(ICell, K) =
             -(Gravity / RhoSw) *
             ((1.0_Real / SpecVol(ICell, K - 1)) -
              (1.0_Real / SpecVol(ICell, K))) /
             (GeomZMid(ICell, K - 1) - GeomZMid(ICell, K));
      }
   }

 private:
   Array2DReal GeomZMid;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Class for Equation of State (EOS) calculations
class Eos {
 public:
   /// Get instance of Eos
   static Eos *getInstance();

   /// Destroy instance (frees Kokkos views)
   static void destroyInstance();

   EosType EosChoice;              ///< Current EOS type in use
   Array2DReal SpecVol;            ///< Specific volume field at level centers
   Array2DReal SpecVolDisplaced;   ///< Displaced specific volume field
   Array2DReal BruntVaisalaFreqSq; ///< Squared Brunt-Vaisala frequency field

   std::string SpecVolFldName; ///< Field name for specific volume
   std::string
       SpecVolDisplacedFldName; ///< Field name for displaced specific volume
   std::string BruntVaisalaFreqSqFldName; ///< Field name for squared
                                          ///< Brunt-Vaisala frequency
   std::string EosGroupName;              ///< EOS group name (for config)
   std::string Name;                      ///< Name of this EOS instance

   /// Compute specific volume for all cells/layers
   void computeSpecVol(const Array2DReal &ConservTemp,
                       const Array2DReal &AbsSalinity,
                       const Array2DReal &Pressure);

   /// Compute displaced specific volume (for vertical displacement)
   void computeSpecVolDisp(const Array2DReal &ConservTemp,
                           const Array2DReal &AbsSalinity,
                           const Array2DReal &Pressure, I4 KDisp);

   /// Compute squared Brunt-Vaisala frequency for all cells/layers
   void computeBruntVaisalaFreqSq(const Array2DReal &ConservTemp,
                                  const Array2DReal &AbsSalinity,
                                  const Array2DReal &Pressure,
                                  const Array2DReal &SpecVol);

   /// Convert Conservative Temperature to potential temperature
   /// GSW Toolbox function gsw_
   Real calcPtFromCt(const Real &Sa, const Real &Ct) const;

   /// Convert potential temperature to Conservative Temperature
   /// GSW Toolbox function gsw_
   Real calcCtFromPt(const Real &Sa, const Real &Pt) const;

   /// Calculate potential enthalpy of ice from potential temperature of ice
   /// GSW Toolbox function gsw_pot_enthalpy_from_pt_ice_poly
   Real calcPotEnthalpyIceFromPotTempIce(const Real &PotTempIce) const;

   /// Calculate frazil properties
   /// GSW Toolbox function gsw_frazil_properties_potential_poly
   void calcFrazilProperties(const Real &BulkAbsSalinity,
                             const Real &BulkPotEnthalpy, const Real &Pressure,
                             Real &InterstitialAbsSalinity,
                             Real &InterstitialConservTemp,
                             Real &IceMassFraction) const;

   /// Calculate melting ice properties
   /// GSW Toolbox function gsw_melting_ice_into_seawater
   void calcMeltingIceIntoSeawater(
       const Real &AbsSalinity, const Real &ConservTemp, const Real &Pressure,
       const Real &IceMassFractionInput, const Real &IceTemperature,
       Real &AbsSalinityFinal, Real &ConservTempFinal,
       Real &IceMassFractionFinal) const;

   /// Calculates the potential temperature of ice from potential enthalpy of
   /// ice GSW Toolbox function gsw_pt_from_pot_enthalpy_ice_poly
   Real calcPotTempIceFromPotEnthalpyIce(const Real &PotEnthalpyIce) const;

   /// Calculates in-situ temperature from the potential tmeperature of ice with
   /// reference pressure of 0 dbar and the in-situ pressure
   /// GSW Toolbox function gsw_t_from_pt0_ice
   Real calcISTempFromPotTempIce(const Real &PotTempIce,
                                 const Real &Pressure) const;

   /// Initialize EOS from config and mesh
   static void init();

 private:
   /// Private constructor
   Eos(const std::string &Name, const HorzMesh *Mesh, const VertCoord *VCoord);

   /// Private destructor
   ~Eos();

   /// Instance pointer
   static Eos *Instance;

   /// Delete copy and move constructors and assignment operators
   Eos(const Eos &)            = delete;
   Eos &operator=(const Eos &) = delete;
   Eos(Eos &&)                 = delete;
   Eos &operator=(Eos &&)      = delete;

   const HorzMesh *Mesh;    ///< Horizontal mesh
   const VertCoord *VCoord; ///< Vertical coordinate

   Teos10Eos ComputeSpecVolTeos10;     ///< TEOS-10 specific volume calculator
   LinearEos ComputeSpecVolLinear;     ///< Linear specific volume calculator
   ConstantEos ComputeSpecVolConstant; ///< Constant specific volume calculator
   Teos10BruntVaisalaFreqSq
       ComputeBruntVaisalaFreqSqTeos10; ///< TEOS-10 squared Brunt-Vaisala
                                        ///< calculator
   LinearBruntVaisalaFreqSq
       ComputeBruntVaisalaFreqSqLinear; ///< Linear squared Brunt-Vaisala
                                        ///< calculator

   // Define fields and metadata
   void defineFields();

}; // End class Eos

} // namespace OMEGA
#endif
