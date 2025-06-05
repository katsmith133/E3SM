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
#include "EosConstants.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include <string>

namespace OMEGA {

enum class EosType {
   Linear,       /// Linear equation of state
   Teos10Poly75t /// Roquet et al. 2015 75 term expansion
};

class Eos {
 public:
   EosType EosChoice;
   Array2DReal SpecVol;
   Array2DReal SpecVolDisplaced;
   Array2DReal SpecVolPCoeffs; // Only used for TEOS-10
   mutable bool PCoeffsInit = false; // Only used for TEOS-10

   // Linear EOS parameters
   Real DRhodT  = {-0.2};   // alpha in kg.m-3 degC-1
   Real DRhodS  = {0.8};    // beta in kg m-3
   Real RhoT0S0 = {1000.0}; // density at (T,S)=(0,0) in kg.m-3

   std::string SpecVolFldName;
   std::string SpecVolDisplacedFldName;
   std::string EosGroupName;
   std::string Name;

   void computeSpecVol(const Array2DReal &SpecVol,
                       const Array2DReal &ConservTemp,
                       const Array2DReal &AbsSalinity,
                       const Array2DReal &Pressure) const;

   // Static management functions
   static I4 init();
   static Eos *create(const std::string &Name, const HorzMesh *Mesh,
                      int NVertLevels) {
      // Check to see if eos of the same name already exist and
      // if so, exit with an error
      if (AllEos.find(Name) != AllEos.end()) {
         LOG_ERROR("Attempted to create Eos with name {} but Eos of "
                   "that name already exists",
                   Name);
         return nullptr;
      }

      // Create new eos on the heap and put it in a map of
      // unique_ptrs, which will manage its lifetime
      auto *NewEos = new Eos(Name, Mesh, NVertLevels);
      AllEos.emplace(Name, NewEos);

      return NewEos; // get(Name);
   }
   ~Eos();
   // Deallocates arrays
   static void clear();
   // Remove Eos object by name
   static void erase(const std::string &Name); ///< [in]
   // Get default eos object
   static Eos *getDefault();
   // Get eos object by name
   static Eos *get(const std::string &Name); ///< [in]

 private:
   I4 NCellsAll;
   I4 NChunks;
   const int NVertLevels;

   // Linear EOS calculation
   KOKKOS_FUNCTION void computeSpecVolLinear(const Array2DReal &SpecVol,
                                             I4 ICell, I4 KChunk,
                                             const Array2DReal &ConservTemp,
                                             const Array2DReal &AbsSalinity) const {
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         SpecVol(ICell, K) =
             1.0 / (RhoT0S0 + (DRhodT * ConservTemp(ICell, K) +
                               DRhodS * AbsSalinity(ICell, K)));
      }
   }

   // TEOS-10 calculation
   KOKKOS_FUNCTION void computeSpecVolTeos10(const Array2DReal &SpecVol,
                                             I4 ICell, I4 KChunk,
                                             const Array2DReal &ConservTemp,
                                             const Array2DReal &AbsSalinity,
                                             const Array2DReal &Pressure) const {
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         calcPCoeffs(SpecVolPCoeffs, KVec, ConservTemp(ICell, K), AbsSalinity(ICell, K));
         PCoeffsInit = true;
         SpecVol(ICell, K) =
             calcRefProfile(Pressure(ICell, K)) +
             calcDelta(SpecVolPCoeffs, KVec, Pressure(ICell, K));
      }
   }

   // Calculate displaced specific volume
   KOKKOS_FUNCTION void calcDisplacedSpecVol(const Array2DReal &SpecVolDisplaced, 
                                             const I4 ICell, const I4 KChunk, 
                                             const Array2DReal &Pressure) const {
      // Checking to make sure calcPCoeffs was called already
      if (!PCoeffsInit) {
         LOG_ERROR("Eos: calcPCoeffs must be called before calcDisplacedSpecVol");
         return;
      }

      OMEGA_SCOPE(LocSpecVolPCoeffs, SpecVolPCoeffs);
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K    = KStart + KVec;
         const I4 KTmp = Kokkos::min(K + 1, NVertLevels);
         SpecVolDisplaced(ICell, K) =
             calcRefProfile(Pressure(ICell, KTmp)) +
             calcDelta(LocSpecVolPCoeffs, KVec, Pressure(ICell, KTmp));
      }
   }

   // TEOS-10 helpers
   // Calculate pressure coefficients
   KOKKOS_FUNCTION void calcPCoeffs(const Array2DReal &SpecVolPCoeffs, 
                                    const I4 K, const Real Ct, 
                                    const Real Sa) const {
      const Real SAu    = 40.0 * 35.16504 / 35.0;
      const Real CTu    = 40.0;
      const Real DeltaS = 24.0;
      GSW_SPECVOL_COEFFICIENTS;
      const Real Ss     = Kokkos::sqrt((Sa + DeltaS) / SAu);
      Real Tt           = Ct / CTu;
      SpecVolPCoeffs(5, K) = V005;
      SpecVolPCoeffs(4, K) = V014 * Tt + V104 * Ss + V004;
      SpecVolPCoeffs(3, K) = 
            (V023 * Tt + V113 * Ss + V013) * Tt + 
            (V203 * Ss + V103) * Ss + V003;
      SpecVolPCoeffs(2, K) = 
            (((V042 * Tt + V132 * Ss + V032) * Tt + 
            (V222 * Ss + V122) * Ss + V022) * Tt +
            ((V312 * Ss + V212) * Ss + V112) * Ss + V012) * Tt +
            (((V402 * Ss + V302) * Ss + V202) * Ss + V102) * Ss + V002;
      SpecVolPCoeffs(1, K) = 
            ((((V051 * Tt + V141 * Ss + V041) * Tt + 
            (V231 * Ss + V131) * Ss + V031) * Tt +
            ((V321 * Ss + V221) * Ss + V121) * Ss + V021) * Tt +
            (((V411 * Ss + V311) * Ss + V211) * Ss + V111) * Ss + V011) * Tt +
            ((((V501 * Ss + V401) * Ss + V301) * Ss + V201) * Ss + V101) * Ss + V001;
      SpecVolPCoeffs(0, K) =
            (((((V060 * Tt + V150 * Ss + V050) * Tt + 
            (V240 * Ss + V140) * Ss + V040) * Tt +
            ((V330 * Ss + V230) * Ss + V130) * Ss + V030) * Tt +
            (((V420 * Ss + V320) * Ss + V220) * Ss + V120) * Ss + V020) * Tt +
            ((((V510 * Ss + V410) * Ss + V310) * Ss + V210) * Ss + V110) * Ss + V010) * Tt +
            (((((V600 * Ss + V500) * Ss + V400) * Ss + V300) * Ss + V200) * Ss + V100) * Ss + V000;

      // could insert a check here (abs(value)> 0 value or <e+33)
   }

   // Calculate delta
   KOKKOS_FUNCTION Real calcDelta(const Array2DReal &SpecVolPCoeffs, 
                                  const I4 K, const Real P) const {
      // Check if the coefficients have been initialized
      if (!PCoeffsInit) {
         LOG_ERROR("Eos: calcPCoeffs must be called before calcDelta");
         return;
      }
      const Real Pu = 1e4;
      Real Pp       = P / Pu;

      Real Delta = ((((SpecVolPCoeffs(5, K) * Pp + SpecVolPCoeffs(4, K)) * Pp +
                      SpecVolPCoeffs(3, K)) * Pp + SpecVolPCoeffs(2, K)) * Pp +
                      SpecVolPCoeffs(1, K)) * Pp + SpecVolPCoeffs(0, K);
      return Delta;
   }

   // Calculate reference profile
   KOKKOS_FUNCTION Real calcRefProfile(Real P) const {
      const Real Pu  = 1e4;
      const Real V00 = -4.4015007269e-05;
      const Real V01 = 6.9232335784e-06;
      const Real V02 = -7.5004675975e-07;
      const Real V03 = 1.7009109288e-08;
      const Real V04 = -1.6884162004e-08;
      const Real V05 = 1.9613503930e-09;
      Real Pp        = P / Pu;

      Real V0 =
          (((((V05 * Pp + V04) * Pp + V03) * Pp + V02) * Pp + V01) * Pp + V00) * Pp;
      return V0;
   }

   // Constructor declaration
   Eos(const std::string &Name, ///< [in] Name for eos object
       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
       int NVertLevels          ///< [in] Number of vertical levels
   );
   // Pointer to default eos
   static Eos *DefaultEos;
   // Map with all eos objects
   static std::map<std::string, std::unique_ptr<Eos>> AllEos;
   void defineFields();
}; // End class Eos

} // namespace OMEGA
#endif