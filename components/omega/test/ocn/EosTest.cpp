//===-- Test driver for OMEGA GSW-C library -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA GSW-C external library
///
/// This driver tests that the GSW-C library can be called
/// and returns expected value (as published in Roquet et al 2015)
//
//===-----------------------------------------------------------------------===/

#include "Eos.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "EosConstants.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "mpi.h"

// added for debug
#include "AuxiliaryState.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"

#include <gswteos-10.h>

using namespace OMEGA;

constexpr Geometry Geom   = Geometry::Spherical;
constexpr int NVertLevels = 60;
// Published values (TEOS-10) to test against
const Real TeosExpValueDelta = 0.0009776149797; 
const Real TeosExpValueVol   = 0.0009732819628; 
double Sa                    = 30.0; // Absolute Salinity in g/kg
double Ct                    = 10.0; // Conservative Temperature in degC
double P                     = 1000.0; // Pressure in dbar
const Real RTol              = 1e-10; // Relative tolerance for isApprox checks
// Expected value for Linear eos (and default parameters)
const Real LinearExpValue = 0.0009784735812133072;

//------------------------------------------------------------------------------
// The initialization routine for Eos testing. It calls various
// init routines, including the creation of the default decomposition.
I4 initEosTest(const std::string &mesh) {

   I4 Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_ERROR("Eos: Error reading config file");
      return Err;
   }

   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(mesh);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing default decomposition");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing default mesh");
   }

   const auto &Mesh = HorzMesh::getDefault();
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);
   return Err;
}

int testEosMapping() {
   int Err = 0;
   // test initialization
   int EosErr = Eos::init();
   if (EosErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing default Eos");
   }

   // test retrieval of default
   Eos *DefEos = Eos::getDefault();

   if (DefEos) {
      LOG_INFO("EosTest: Default Eos retrieval PASS");
   } else {
      Err++;
      LOG_INFO("EosTest: Default Eos retrieval FAIL");
      return -1;
   }

   const auto *Mesh = HorzMesh::getDefault();
   // test creation of another Eos
   Eos::create("TestEos", Mesh, NVertLevels);

   if (Eos::get("TestEos")) {
      LOG_INFO("EosTest: Non-default Eos retrieval PASS");
   } else {
      Err++;
      LOG_INFO("EosTest: Non-default Eos retrieval FAIL");
   }

   // test erase
   Eos::erase("TestEos");
   if (Eos::get("TestEos")) {
      Err++;
      LOG_INFO("EosTest: Non-default Eos erase FAIL");
   } else {
      LOG_INFO("EosTest: Non-default Eos erase PASS");
   }

   return Err;
}

int testEosLinear() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   // create Eos to test
   Eos::create("LinearEos", Mesh, NVertLevels);
   Eos *TestEos       = Eos::get("LinearEos");
   TestEos->EosChoice = EosType::Linear;

   // create ocean state array
   Array2DReal SArray = Array2DReal("SArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal TArray = Array2DReal("TArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal PArray = Array2DReal("PArray", Mesh->NCellsAll, NVertLevels);
   // Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   // Key calculation
   TestEos->computeSpecVol(TestEos->SpecVol, TArray, SArray, PArray);

   // check on all array values
   int numMismatches = 0;
   Array2DReal SpecVol = TestEos->SpecVol;
   parallelReduce("CheckSpecVolMatrix-linear", {Mesh->NCellsAll, NVertLevels},
                  KOKKOS_LAMBDA(int i, int j, int &localCount) {
                     if (!isApprox(SpecVol(i, j), LinearExpValue, RTol)) {
                        localCount++;
                     }
                  },
                  numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: SpecVol Linear isApprox FAIL, "
                "expected {}, got {} mismatches",
                LinearExpValue, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalc Linear: PASS");
   }

   Eos::erase("LinearEos");
   return Err;
}

int testEosTeos10() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   // create Eos to test
   Eos::create("TeosEos", Mesh, NVertLevels);
   Eos *TestEos       = Eos::get("TeosEos");
   TestEos->EosChoice = EosType::Teos10Poly75t;

   // create ocean state array
   Array2DReal SArray = Array2DReal("SArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal TArray = Array2DReal("TArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal PArray = Array2DReal("PArray", Mesh->NCellsAll, NVertLevels);
   // Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   // Key calculation
   TestEos->computeSpecVol(TestEos->SpecVol, TArray, SArray, PArray);

   // check on all array values
   int numMismatches = 0;
   Array2DReal SpecVol = TestEos->SpecVol;
   parallelReduce("CheckSpecVolMatrix-Teos", {Mesh->NCellsAll, NVertLevels},
                  KOKKOS_LAMBDA(int i, int j, int &localCount) {
                     if (!isApprox(SpecVol(i, j), TeosExpValueVol, RTol)) {
                        localCount++;
                     }
                  },
                  numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: TEOS SpecVol isApprox FAIL, "
                "expected {}, got {} mismatches",
                TeosExpValueVol, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalc TEOS-10: PASS");
   }

   Eos::erase("TeosEos");
   return Err;
}

void finalizeEosTest() {
   Eos::clear();
   HorzMesh::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

int checkValueGswcSpecVol() {
   int Err         = 0;
   const Real RTol = 1e-10;

   double SpecVol = gsw_specvol(Sa, Ct, P);
   bool Check = isApprox(SpecVol, TeosExpValueVol, RTol);
   if (!Check) {
      Err++;
      LOG_ERROR(
          "checkValueGswcSpecVol: SpecVol isApprox FAIL, expected {}, got {}",
          TeosExpValueVol, SpecVol);
   }
   if (Err == 0) {
      LOG_INFO("checkValueGswcSpecVol: PASS");
   }
   return Err;
}

// intermediate test accessing the coefficients
int fetchCoeff() {
   int Err       = 0;
   double ExpVal = 0.0010769995862;
   GSW_SPECVOL_COEFFICIENTS;

   if (!isApprox(V000, ExpVal, RTol)) {
      Err++;
      LOG_ERROR("EosTest: Coeff V000 isApprox FAIL, expected {}, got {}",
                ExpVal, V000);
   }
   if (Err == 0) {
      LOG_INFO("GswcTeosTest: check PASS");
   }
   return Err;
}

// the main test (all in one to have the same log)
// Single value tests:
// --> one test calls the external GSW-C library
// and compares the specific volume to the published value
// --> next test call the external GSW-C library
// and compares the V000 coefficient to the expected value
// Full array tests:
// --> one tests the initialization/retrieval of Eos
// --> next checks the value on a Eos with linear option
// --> next checks the value on a Eos with TEOS-10 option
int eosTest(const std::string &MeshFile = "OmegaMesh.nc") {
   int Err = initEosTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("EosTest: Error initializing");
   }
   const auto &Mesh = HorzMesh::getDefault();
   
   LOG_INFO("Single value checks:");
   Err += checkValueGswcSpecVol();
   Err += fetchCoeff();

   LOG_INFO("Full array checks:");
   Err += testEosMapping();
   Err += testEosLinear();
   Err += testEosTeos10();

   if (Err == 0) {
      LOG_INFO("EosTest: Successful completion");
   }
   finalizeEosTest();

   return Err;
}

//------------------------------------------------------------------------------
// The test driver for Eos testing
int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   { RetVal += eosTest(); }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
