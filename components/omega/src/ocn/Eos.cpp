#include "Eos.h"
#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {

Eos *Eos::DefaultEos = nullptr;
std::map<std::string, std::unique_ptr<Eos>> Eos::AllEos;

Teos10Poly75t::Teos10Poly75t(int NVertLevels) : NVertLevels(NVertLevels) {
   SpecVolPCoeffs = Array2DReal("SpecVolPCoeffs", 6, VecLength);
}

LinearEOS::LinearEOS() {}

Eos::Eos(const std::string &Name_, ///< [in] Name for eos object
         const HorzMesh *Mesh,     ///< [in] Horizontal mesh
         int NVertLevels           ///< [in] Number of vertical levels
         )
    : computeSpecVolTeos10Poly75t(NVertLevels) {
   SpecVol = Array2DReal("SpecVol", Mesh->NCellsAll, NVertLevels);
   SpecVolDisplaced =
       Array2DReal("SpecVolDisplaced", Mesh->NCellsAll, NVertLevels);
   // Array dimension lengths
   NCellsAll = Mesh->NCellsAll;
   NChunks   = NVertLevels / VecLength;
   Name      = Name_;

   // Register fields and metadata for IO
   defineFields();

} // end constructor

// Initializes the Eos (Equation of State) class and its options.
// it ASSUMES that HorzMesh was initialized and initializes the Eos class by
// using the default mesh, reading the config file, and setting parameters
// for either a Linear or TEOS-10 equation.
// Returns 0 on success, or an error code if any required option is missing.
int Eos::init() {

   int Err               = 0;
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   I4 NVertLevels        = DefHorzMesh->NVertLevels;

   // Create default eos
   Eos::DefaultEos = create("Default", DefHorzMesh, NVertLevels);

   // Get EosConfig group
   Config *OmegaConfig = Config::getOmegaConfig();
   Config EosConfig("Eos");
   Err = OmegaConfig->get(EosConfig);
   if (Err != 0) {
      LOG_CRITICAL("Eos: Eos group not found in Config");
      return Err;
   }
   std::string EosTypeStr;
   Err = EosConfig.get("EosType", EosTypeStr);
   if (Err != 0) {
      LOG_CRITICAL("Eos: EosType subgroup not found in EosConfig");
      return Err;
   }

   if (EosTypeStr == "Linear" or EosTypeStr == "linear") {
      Config EosLinConfig("Linear");
      Err = EosConfig.get(EosLinConfig);
      if (Err != 0) {
         LOG_CRITICAL("Eos: Linear subgroup not found in EosConfig");
         return Err;
      }
      DefaultEos->EosChoice = EosType::Linear;
      Err                   = EosLinConfig.get("DRhoDT",
                                            DefaultEos->computeSpecVolLinear.DRhodT);
      if (Err != 0) {
         LOG_CRITICAL("Eos: Parameter Linear:DRhodT not found in EosLinConfig");
         return Err;
      }
      Err = EosLinConfig.get("DRhoDS",
                          DefaultEos->computeSpecVolLinear.DRhodS);
      if (Err != 0) {
         LOG_CRITICAL("Eos: Parameter Linear:DRhodS not found in EosLinConfig");
         return Err;
      }
      Err = EosLinConfig.get("RhoT0S0",
                          DefaultEos->computeSpecVolLinear.RhoT0S0);
      if (Err != 0) {
         LOG_CRITICAL("Eos: Parameter Linear:RhoT0S0 not found in EosLinConfig");
      }
   } else if ((EosTypeStr == "teos10") or (EosTypeStr == "teos-10") or
              (EosTypeStr == "TEOS-10")) {
      DefaultEos->EosChoice = EosType::Teos10Poly75t;
   } else {
      LOG_CRITICAL("Eos: Unknown EosType requested");
      Err = -1;
      return Err;
   }

   return Err;
} // end init

void Eos::computeSpecVol(const Array2DReal &SpecVol,
                         const Array2DReal &ConservTemp,
                         const Array2DReal &AbsSalinity,
                         const Array2DReal &Pressure) const {
   OMEGA_SCOPE(LocSpecVol, SpecVol);
   OMEGA_SCOPE(LocComputeSpecVolLinear, computeSpecVolLinear);
   OMEGA_SCOPE(LocComputeSpecVolTeos10Poly75t, computeSpecVolTeos10Poly75t);
   deepCopy(LocSpecVol, 0.0);
   if (EosChoice == EosType::Linear) {
      parallelFor(
          "eos-linear", {NCellsAll, NChunks},
          KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
             LocComputeSpecVolLinear(LocSpecVol, ICell, KChunk,
                                     ConservTemp, AbsSalinity);
          });
   } else if (EosChoice == EosType::Teos10Poly75t) {
      parallelFor(
          "eos-teos10", {NCellsAll, NChunks},
          KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
             LocComputeSpecVolTeos10Poly75t(LocSpecVol, ICell, KChunk,
                                            ConservTemp,
                                            AbsSalinity, Pressure);
          });
   }
}

//------------------------------------------------------------------------------
// Define IO fields and metadata
void Eos::defineFields() {

   int Err = 0;

   SpecVolFldName          = "SpecVol";
   SpecVolDisplacedFldName = "SpecVolDisplaced";
   if (Name != "Default") {
      SpecVolFldName.append(Name);
      SpecVolDisplacedFldName.append(Name);
   }

   // Create fields for state variables
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";

   auto SpecVolField =
       Field::create(SpecVolFldName,                   // Field name
                     "Layer-averaged Specific Volume", /// long Name
                     "m3 kg-1",                        // units
                     "sea_water_specific_volume",      // CF-ish Name
                     0.0,                              // min valid value
                     9.99E+30,                         // max valid value
                     -9.99E+30, // scalar used for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );
   auto SpecVolDisplacedField =
       Field::create(SpecVolDisplacedFldName, // Field name
                     "Specific Volume displaced adiabatically "
                     "to 1 layer below",                    /// long Name
                     "m3 kg-1",                             // units
                     "sea_water_specific_volume_displaced", // CF-ish Name
                     0.0,                                   // min valid value
                     9.99E+30,                              // max valid value
                     -9.99E+30, // scalar used for undefined entried
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   // Create a field group for the eos-specific state fields
   EosGroupName = "Eos";
   if (Name != "Default") {
      EosGroupName.append(Name);
   }
   auto EosGroup = FieldGroup::create(EosGroupName);

   // Add restart group if needed
   if (!FieldGroup::exists("Restart"))
      auto RestartGroup = FieldGroup::create("Restart");

   Err = EosGroup->addField(SpecVolDisplacedFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", SpecVolDisplacedFldName,
                EosGroupName);
   Err = EosGroup->addField(SpecVolFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", SpecVolFldName,
                EosGroupName);

   Err = FieldGroup::addFieldToGroup(SpecVolDisplacedFldName, "Restart");
   if (Err != 0)
      LOG_ERROR("Error adding {} to Restart field group",
                SpecVolDisplacedFldName);
   Err = FieldGroup::addFieldToGroup(SpecVolFldName, "Restart");
   if (Err != 0)
      LOG_ERROR("Error adding {} to Restart field group", SpecVolFldName);

   // Associate Field with data
   Err = SpecVolDisplacedField->attachData<Array2DReal>(SpecVolDisplaced);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}",
                SpecVolDisplacedFldName);
   Err = SpecVolField->attachData<Array2DReal>(SpecVol);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}", SpecVolFldName);

} // end defineIOFields

//------------------------------------------------------------------------------
// Destroys the eos class
Eos::~Eos() {

   // Kokkos arrays removed when no longer in scope
   int Err;
   Err = FieldGroup::destroy(EosGroupName);
   if (Err != 0)
      LOG_ERROR("Error removing FieldGroup {}", EosGroupName);
   Err = Field::destroy(SpecVolFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", SpecVolFldName);
   Err = Field::destroy(SpecVolDisplacedFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", SpecVolDisplacedFldName);

} // end destructor

//------------------------------------------------------------------------------
// Removes all eos instances before exit
void Eos::clear() { AllEos.clear(); } // end clear

//------------------------------------------------------------------------------
// Removes eos from list by name
void Eos::erase(const std::string &Name) { AllEos.erase(Name); } // end erase

//------------------------------------------------------------------------------
// Get default eos
Eos *Eos::getDefault() { return Eos::DefaultEos; } // end get default

//------------------------------------------------------------------------------
// Get eos by name
Eos *Eos::get(const std::string &Name ///< [in] Name of eos
) {

   auto it = AllEos.find(Name);

   if (it != AllEos.end()) {
      return it->second.get();
   } else {
      LOG_ERROR("Eos::get: Attempt to retrieve non-existent eos:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get eos

} // namespace OMEGA
