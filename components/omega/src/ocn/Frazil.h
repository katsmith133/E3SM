#ifndef OMEGA_FRAZIL_H
#define OMEGA_FRAZIL_H
//===-- ocn/Frazil.h - Frazil Ice Formation -------------------*- C++ -*-===//
//
// The Frazil class manages frazil tendencies and accumulators.
// This initial implementation only has a basic configuration.
// but carries scaffolding for other implementations.
//
//===----------------------------------------------------------------------===//

#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

enum class FrazilType {
   BasicFrazil ///< MPAS-O style basic frazil option
   // Additional options can be added here in the future.
};

class BasicFrazilFormation {
 public:
   BasicFrazilFormation();

   Real massLimit         = 0.1_Real;  // to do:  remove default
   Real FrazilIceSalinity = IceRefSal; // Global constant
   Real LatFrazil         = LatIce;    // Global constant

   Real FrazilPorosity =
       1.0_Real; // Internal for now; can move to config later.

   KOKKOS_FUNCTION void operator()(const Real SA, const Real CT, const Real PDb,
                                   const Real H, Real &SumIceThickness,
                                   Real &SumSalt, Real &SumEnergy, Real &HTend,
                                   Real &TTend, Real &STend,
                                   const Real Tfrz) const {

      const Real potential      = H * Cp0Sw * RhoSw * (CT - Tfrz);
      const Real freezingEnergy = Kokkos::max(0.0_Real, -potential);

      HTend = 0.0_Real;
      TTend = 0.0_Real;
      STend = 0.0_Real;

      Real newFrzThickness =
          freezingEnergy /
          (LatFrazil * RhoSw); // frazil (ice) mass in pseudo-thickness terms
      newFrzThickness = Kokkos::min(newFrzThickness, H * massLimit);
      const Real newFrzEnergy =
          -newFrzThickness * LatFrazil; // (<0; enthalpy of ice)

      // MANUAL TOGGLE: uncomment line below to use porosity
      // const Real FrazilIceSalinity = FrazilPorosity * SA;

      const Real frazilSalinity = Kokkos::min(FrazilIceSalinity, SA);
      const Real newSaltContent =
          newFrzThickness * frazilSalinity; // in m.(g/kg)

      HTend = -newFrzThickness - newSaltContent;
      // TTend below should include the enthalpy associated with the mass flux
      // this is intentionally not added here to match the mpas-o implementation
      TTend =
          -(newFrzEnergy) / Cp0Sw; // (E< 0 so TTend>0) // scaled to h.CT tend
      STend = -newSaltContent;

      SumIceThickness += newFrzThickness;
      // non-conservation between Sum and HTend by construction in the original
      SumSalt += newSaltContent;
      SumEnergy += newFrzEnergy;
   }
};

// TO-DO: is there a new for salt reconciliation at the surface?
// frazilSalinityTendency(minLevelCell(iCell),iCell) =
// frazilSalinityTendency(minLevelCell(iCell),iCell) + &
//     max(0.0_RKIND,(sumNewThicknessWeightedSaltContent -
//     newThicknessWeightedSaltContent) ) / dt
// accumulatedFrazilIceSalinityNew(iCell) =
// accumulatedFrazilIceSalinityOld(iCell) + newThicknessWeightedSaltContent

class BasicFrazilMelt {
 public:
   BasicFrazilMelt();

   Real massLimit         = 0.1_Real;  // to do:  remove default
   Real FrazilIceSalinity = IceRefSal; // Global constant
   Real LatFrazil         = LatIce;    // Global constant

   KOKKOS_FUNCTION void operator()(const Real SA, const Real CT, const Real PDb,
                                   const Real H, Real &SumIceThickness,
                                   Real &SumSalt, Real &SumEnergy, Real &HTend,
                                   Real &TTend, Real &STend,
                                   const Real Tfrz) const {
      constexpr Real Eps = 1.0e-12_Real;

      if (SumIceThickness <= Eps) { // potential leak if we dont redistribute
         SumIceThickness = 0.0_Real;
         SumSalt         = 0.0_Real;
         HTend           = 0.0_Real;
         TTend           = 0.0_Real;
         STend           = 0.0_Real;
         return;
      }

      const Real potential       = H * Cp0Sw * RhoSw * (CT - Tfrz);
      const Real availableEnergy = Kokkos::max(0.0_Real, potential);

      HTend = 0.0_Real;
      TTend = 0.0_Real;
      STend = 0.0_Real;

      Real meltThickness =
          availableEnergy /
          (LatFrazil * RhoSw); // mass in pseudo-thickness units
      meltThickness = Kokkos::min(meltThickness, SumIceThickness);
      meltThickness = Kokkos::min(meltThickness,
                                  H * massLimit); // also 0.1h lim on added mass
      const Real meltAverageSalinity = SumSalt / SumIceThickness;
      const Real meltingEnergy =
          meltThickness * LatFrazil; // (energy needed to melt)

      HTend = +meltThickness * (1 + meltAverageSalinity); // (>0 so HTend>0)
      TTend =
          -meltingEnergy / Cp0Sw + // (TTend < 0 when melting for phase change)
          meltThickness *
              Tfrz; // (added the enthalpy of the melt water at freezing temp)
      STend = +meltAverageSalinity * meltThickness; // (STend > 0 when melting)

      SumIceThickness -= meltThickness;
      SumSalt -= meltAverageSalinity * meltThickness;
      SumEnergy -=
          meltingEnergy; // necessarily non-conservative by construction
   }
};

class Frazil {
 public:
   static void init();
   /// Creates a new frazil object and stores it in the AllFrazil map.
   static Frazil *create(const std::string &Name);

   /// Retrieve frazil object by name.
   static Frazil *get(const std::string &Name);

   /// Retrieve default frazil object.
   static Frazil *getDefault();

   /// Destructor
   ~Frazil();

   /// Deallocates arrays
   static void clear();

   /// Remove frazil object by name.
   static void erase(std::string InName); ///< [in] name to remove

   Array2DReal FrazilTTend;
   Array2DReal FrazilSTend;
   Array2DReal FrazilHTend;
   Array1DReal AccMIce;
   Array1DReal AccEIce;
   Array1DReal AccMLiq;
   Array1DReal AccELiq;
   Array1DReal AccMSalt;

   void computeFrazil(const Array2DReal &CT, const Array2DReal &SA,
                      const Array2DReal &P, const Array2DReal &H);
   void computeFrazilBasicImpl(const Array2DReal &CT, const Array2DReal &SA,
                               const Array2DReal &P, const Array2DReal &LayerH);

   bool conservationCheck = false;
   Real depthLimit        = -1.0_Real;

 private:
   static Frazil *DefaultFrazil;
   static std::map<std::string, std::unique_ptr<Frazil>> AllFrazil;

   Frazil(const HorzMesh *Mesh, const VertCoord *VCoord);

   // Forbid copy and move construction/assignment.
   Frazil(const Frazil &)            = delete;
   Frazil &operator=(const Frazil &) = delete;
   Frazil(Frazil &&)                 = delete;
   Frazil &operator=(Frazil &&)      = delete;

   FrazilType frazilChoice;
   BasicFrazilFormation computeBasicFrazilFormation;
   BasicFrazilMelt computeBasicFrazilMelt;

   I4 NCellsAll;
   I4 NChunks;

   const HorzMesh *MeshPtr;
   const VertCoord *VCoordPtr;

   void checkColumnConservation() const;
};

} // namespace OMEGA

#endif
