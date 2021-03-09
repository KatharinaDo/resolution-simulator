// Simon Spannagel (DESY) January 2016

#include "TCanvas.h"
#include "TProfile.h"
#include "TString.h"
#include "TFile.h"

#include "assembly.h"
#include "propagate.h"
#include "materials.h"
#include "constants.h"
#include "log.h"

using namespace std;
using namespace gblsim;
using namespace unilog;

int main(int argc, char* argv[]) {

  /*
   * Telescope resolution simulation for the DATURA telescope at the DESY TB21 beam line
   * Six MIMOSA26 planes with 150mm spacing, intrinsic sensor resolution 3.4um
   * DUT with variable thickness (scan)
   */

  
  Log::ReportingLevel() = Log::FromString("INFO");

  for (int i = 1; i < argc; i++) {
    // Setting verbosity:
    if (std::string(argv[i]) == "-v") { 
      Log::ReportingLevel() = Log::FromString(std::string(argv[++i]));
      continue;
    } 
  }

  TFile * out = TFile::Open("datura-resolution.root","RECREATE");
  gDirectory->pwd();

  TCanvas *c1 = new TCanvas("c1","resolution",700,700);
  TProfile *resolution = new TProfile("resolution"," ",100,0,0.02);

  //----------------------------------------------------------------------------
  // Preparation of the telescope and beam properties:

  // MIMOSA26 telescope planes consist of 50um silicon plus 2x25um Kapton foil only:
  double MIM26 = 55e-3 / X0_Si + 50e-3 / X0_Kapton;
  // The intrinsic resolution has been measured to be around 3.25um:
  double RES = 3.31e-3; // according to paper : 3.24 +- 0.09
  // DUT
  double CLICTDBudget = 0.025;
  // TPX3
  double TPX3Budget = 0.01068;
  double TPX3Res = 12.8e-3;

  // M26  M26  M26      DUT      M26  M26  M26
  //  |    |    |        |        |    |    |
  //  |    |    |        |        |    |    |
  //  |<-->|    |<------>|<------>|    |    |
  //   DIST      DUT_DIST DUT_DIST
  
  // Distance between telescope planes in mm:
  //double DIST = 20;
  // Distance of telescope arms and DUT assembly:
  //double DUT_DIST = 20;

//   double POSITIONS [7] = {0.0, 277.0, 305.0, 351.0, 377.0, 624.0, 669.0}; // A1/-3V/-3V -> no cut-out -> material: 0.0217
//   double DUT_POS = 332.0;
   //double POSITIONS [7] = {0.0, 277.0, 305.0, 364.0, 391.0, 664.0, 691.0}; // A4/-6V/-6V -> with cut-out -> material: 0.0032 
   //double DUT_POS = 344.0;
  // double POSITIONS [7] = {0.0, 277.0, 305.0, 476.0, 504.0, 752.0, 797.0}; // A4/-6V/-6V (rotations sp) -> with cut-out -> material: 0.00032 
  // double DUT_POS = 424.0;
     double POSITIONS [7] = {0.0, 277.0, 305.0, 482.0, 510.0, 758.0, 803.0}; // A4/-6V/-6V (rotations sp) - 2 -> with cut-out -> material: 0.00032
     double DUT_POS = 424.5;
   //double POSITIONS [7] = {0.0, 153.0, 305.0, 345.0, 455.0, 565.0, 621.0};
  // double DUT_POS = 331.0;


  // Beam energy 5 GeV electrons/positrons at DESY:
  double BEAM = 5.4;


  //----------------------------------------------------------------------------
  // Build the trajectory through the telescope device:

  // Build a vector of all telescope planes:
  std::vector<plane> datura;
  double position = 0;
  LOG(logRESULT) << "MIMOSA budget: " << MIM26;

  // Upstream telescope arm:
  for(int i = 0; i < 3; i++) {
    LOG(logRESULT) << "Position of plane: " << i << ": " << POSITIONS[i];
    position = POSITIONS[i];
    datura.push_back(plane(position,MIM26,true,RES));
  }

  // Downstream telescope arm:
  for(int i = 3; i < 6; i++) {
    LOG(logRESULT) << "Position of plane: " << i << ": " << POSITIONS[i];
    position = POSITIONS[i];
    datura.push_back(plane(position,MIM26,true,RES));
  }

  // set the timepix
  LOG(logRESULT) << "Position of plane: " << 7 << ": " << POSITIONS[6];
  datura.push_back(plane(POSITIONS[6],TPX3Budget,true,TPX3Res));

  //for(double dut_x0 = 0.00032; dut_x0 < 0.03; dut_x0 += 0.0005) {
    for(double dut_x0 = 0.010; dut_x0 < 0.0242; dut_x0 += 0.0004) {

    // Prepare the DUT (no measurement, just scatterer
    plane dut(DUT_POS, dut_x0, false);

    // Duplicate the planes vector and add the current DUT:
    std::vector<plane> planes = datura;
    planes.push_back(dut);

    // Build the telescope:
    telescope mytel(planes, BEAM);

    // Get the resolution at plane-vector position (x):
    LOG(logRESULT) << "Track resolution at DUT with " << dut_x0 << "% X0: " << mytel.getResolution(3);
    resolution->Fill(dut_x0,mytel.getResolution(3),1);
  }

    c1->cd();
  resolution->SetTitle("DATURA Track Resolution at DUT;DUT material budget x/X_{0};resolution at DUT #left[#mum#right]");
  resolution->GetYaxis()->SetRangeUser(1.5,4);
  resolution->SetMarkerStyle(0);
  resolution->SetLineColor(kRed+1);
  resolution->SetLineWidth(2);
  resolution->SetMarkerColor(kRed+1);
  
  resolution->Draw();
  c1->Write();

  // Write result to file
  out->Write();
  return 0;
}
