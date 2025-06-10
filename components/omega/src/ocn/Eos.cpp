//===-- ocn/Eos.cpp - Equation of State ------------------*- C++ -*-===//
//
// The Eos class is responsible for managing the equation of state. It 
// has a linear EOS and TEOS-10 EOS option, which is determined at 
// initialization. It contains arrays that store the specific volume and 
// displaced specific volume data.
//
//===----------------------------------------------------------------------===//

#include "Eos.h"
#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {

/// Constructor for Eos
Eos::Eos(const std::string &Name,  ///< [in] Name for eos object
         const HorzMesh *Mesh,     ///< [in] Horizontal mesh
         int NVertLevels           ///< [in] Number of vertical levels
         )
    : 
      SpecVol("SpecVol", Mesh->NCellsAll, NVertLevels),
      SpecVolDisplaced("SpecVolDisplaced", Mesh->NCellsAll, NVertLevels),
      SpecVolPCoeffs("SpecVolPCoeffs", 6, VecLength),
      Name(Name),
      NCellsAll(Mesh->NCellsAll),
      NChunks(NVertLevels / VecLength),
      NVertLevels(NVertLevels)
{
    defineFields();
}

/// Destructor for Eos 
Eos::~Eos() {}

/// Instance management
Eos* Eos::instance_ = nullptr;

/// Get instance of Eos
Eos* Eos::getInstance(const std::string &Name, const HorzMesh *Mesh, int NVertLevels) {
   /// Create instance if it doesn't exist
   if (!instance_) {
       instance_ = new Eos(Name, Mesh, NVertLevels);
   }
   return instance_;
}

/// Destroy instance of Eos
void Eos::destroyInstance() {
   delete instance_;
   instance_ = nullptr;
}

/// Initializes the Eos (Equation of State) class and its options.
/// it ASSUMES that HorzMesh was initialized and initializes the Eos class by
/// using the default mesh, reading the config file, and setting parameters
/// for either a Linear or TEOS-10 equation.
/// Returns 0 on success, or an error code if any required option is missing.
int Eos::init() {

   int Err               = 0;
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   I4 NVertLevels        = DefHorzMesh->NVertLevels;

   /// Create default eos
   Eos* eos = Eos::getInstance("Default", DefHorzMesh, NVertLevels);

   /// Get EosConfig group from Omega config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config EosConfig("Eos");
   Err = OmegaConfig->get(EosConfig);
   if (Err != 0) {
      LOG_CRITICAL("Eos::init: Eos group not found in Config");
      return Err;
   }

   /// Get EosType from EosConfig
   /// and set the EosChoice accordingly
   std::string EosTypeStr;
   Err = EosConfig.get("EosType", EosTypeStr);
   if (Err != 0) {
      LOG_CRITICAL("Eos::init: EosType subgroup not found in EosConfig");
      return Err;
   }

   /// Set EosChoice based on EosTypeStr and get parameters
   if (EosTypeStr == "Linear" or EosTypeStr == "linear") {
      Config EosLinConfig("Linear");
      Err = EosConfig.get(EosLinConfig);
      if (Err != 0) {
         LOG_CRITICAL("Eos::init: Linear subgroup not found in EosConfig");
         return Err;
      }
      eos->EosChoice = EosType::Linear;
      Err                   = EosLinConfig.get("DRhoDT", eos->DRhodT);
      if (Err != 0) {
         LOG_CRITICAL("Eos::init: Parameter Linear:DRhodT not found in EosLinConfig");
         return Err;
      }
      Err = EosLinConfig.get("DRhoDS", eos->DRhodS);
      if (Err != 0) {
         LOG_CRITICAL("Eos::init: Parameter Linear:DRhodS not found in EosLinConfig");
         return Err;
      }
      Err = EosLinConfig.get("RhoT0S0", eos->RhoT0S0);
      if (Err != 0) {
         LOG_CRITICAL("Eos::init: Parameter Linear:RhoT0S0 not found in EosLinConfig");
      }
   } else if ((EosTypeStr == "teos10") or (EosTypeStr == "teos-10") or
              (EosTypeStr == "TEOS-10")) {
      eos->EosChoice = EosType::Teos10Poly75t;
   } else {
      LOG_CRITICAL("Eos::init: Unknown EosType requested");
      Err = -1;
      return Err;
   }

   return Err;
} // end init

/// Compute specific volume for all cells/levels (no displacement)
void Eos::computeSpecVol(Array2DReal SpecVol,
                         const Array2DReal &ConservTemp,
                         const Array2DReal &AbsSalinity,
                         const Array2DReal &Pressure) {
   OMEGA_SCOPE(LocSpecVol, SpecVol); /// Create a local view for computation
   deepCopy(LocSpecVol, 0.0);
   I4 KDisp = 0; /// No displacement in this case
   /// Dispatch to the correct EOS calculation
   if (EosChoice == EosType::Linear) {
      parallelFor(
          "eos-linear", {NCellsAll, NChunks},
          KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
             computeSpecVolLinear(LocSpecVol, ICell, KChunk,
                                     ConservTemp, AbsSalinity);
          });
   } else if (EosChoice == EosType::Teos10Poly75t) {
      parallelFor(
          "eos-teos10", {NCellsAll, NChunks},
          KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
             computeSpecVolTeos10(LocSpecVol, ICell, KChunk,
                                            ConservTemp, AbsSalinity, 
                                            Pressure, KDisp);
          });
   }
   deepCopy(SpecVol, LocSpecVol); /// Copy result back to output view
}

/// Compute displaced specific volume (for vertical displacement)
void Eos::computeSpecVolDisp(Array2DReal SpecVolDisplaced,
                         const Array2DReal &ConservTemp,
                         const Array2DReal &AbsSalinity,
                         const Array2DReal &Pressure,
                         I4 KDisp) {
   OMEGA_SCOPE(LocSpecVolDisplaced, SpecVolDisplaced); /// Local view for computation
   deepCopy(LocSpecVolDisplaced, 0.0);
   if (EosChoice == EosType::Linear) {
      LOG_INFO(
          "Eos::computeSpecVolDisp called with Linear EOS. "
          "SpecVol is independent of pressure/depth, so the "
          "displaced value will be the same as SpecVol.");
      parallelFor(
          "eos-linear", {NCellsAll, NChunks},
          KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
             computeSpecVolLinear(LocSpecVolDisplaced, ICell, KChunk,
                                     ConservTemp, AbsSalinity);
          });
   } else if (EosChoice == EosType::Teos10Poly75t) {
      parallelFor(
          "eos-teos10", {NCellsAll, NChunks},
          KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
             computeSpecVolTeos10(LocSpecVolDisplaced, ICell, KChunk,
                                            ConservTemp, AbsSalinity, 
                                            Pressure, KDisp);
          });
   }
   deepCopy(SpecVolDisplaced, LocSpecVolDisplaced); /// Copy result back to output view
}

/// Define IO fields and metadata for output
void Eos::defineFields() {

   int Err = 0;

   /// Set field names (append Name if not default)
   SpecVolFldName          = "SpecVol";
   SpecVolDisplacedFldName = "SpecVolDisplaced";
   if (Name != "Default") {
      SpecVolFldName.append(Name);
      SpecVolDisplacedFldName.append(Name);
   }

   /// Create fields for state variables
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";

   /// Create and register the specific volume field
   auto SpecVolField =
       Field::create(SpecVolFldName,                   // Field name
                     "Layer-averaged Specific Volume",  // Long Name
                     "m3 kg-1",                        // Units
                     "sea_water_specific_volume",       // CF-ish Name
                     0.0,                              // Min valid value
                     9.99E+30,                         // Max valid value
                     -9.99E+30,                        // Scalar used for undefined entries
                     NDims,                            // Number of dimensions
                     DimNames                          // Dimension names
       );
   /// Create and register the displaced specific volume field
   auto SpecVolDisplacedField =
       Field::create(SpecVolDisplacedFldName,                  // Field name
                     "Specific Volume displaced adiabatically "
                     "to specified layer",                      // long Name
                     "m3 kg-1",                                // Units
                     "sea_water_specific_volume_displaced",     // CF-ish Name
                     0.0,                                      // Min valid value
                     9.99E+30,                                 // Max valid value
                     -9.99E+30,                                // Scalar used for undefined entried
                     NDims,                                    // Number of dimensions
                     DimNames                                  // Dimension names
       );

   // Create a field group for the eos-specific state fields
   EosGroupName = "Eos";
   if (Name != "Default") {
      EosGroupName.append(Name);
   }
   auto EosGroup = FieldGroup::create(EosGroupName);

   // Add fields to the EOS group
   Err = EosGroup->addField(SpecVolDisplacedFldName);
   if (Err != 0)
      LOG_ERROR("Eos::defineFields: Error adding {} to field group {}", 
                SpecVolDisplacedFldName, EosGroupName);
   Err = EosGroup->addField(SpecVolFldName);
   if (Err != 0)
      LOG_ERROR("Eos::defineFields: Error adding {} to field group {}", 
                SpecVolFldName, EosGroupName);

   // Attach Kokkos views to the fields
   Err = SpecVolDisplacedField->attachData<Array2DReal>(SpecVolDisplaced);
   if (Err != 0)
      LOG_ERROR("Eos::defineFields: Error attaching data array to field {}",
                SpecVolDisplacedFldName);
   Err = SpecVolField->attachData<Array2DReal>(SpecVol);
   if (Err != 0)
      LOG_ERROR("Eos::defineFields: Error attaching data array to field {}", SpecVolFldName);

} // end defineIOFields

} // namespace OMEGA