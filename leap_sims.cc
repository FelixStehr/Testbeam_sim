//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file polarisation/Pol01/Pol01.cc
/// \brief Main program of the polarisation/Pol01 example
//
// $Id: Pol01.cc 86418 2014-11-11 10:39:38Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"


#include "G4VisExecutive.hh"



#include "G4UIExecutive.hh"

// include optical physics
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " leap_sims [-m macro ] [-f outFileName] [-t outType] [-v version] [-b beamline]" << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

G4String session;
G4String macro;
G4String outFile;
G4String outType;
G4String version;
G4String beamline;

for ( G4int i=1; i<argc; i=i+2 ) {
  if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
  else if ( G4String(argv[i]) == "-f" ) outFile = argv[i+1];
  else if ( G4String(argv[i]) == "-t" ) outType = argv[i+1];
  else if ( G4String(argv[i]) == "-v" ) version = argv[i+1];
  else if ( G4String(argv[i]) == "-b" ) beamline = argv[i+1];
  else{
    PrintUsage();
    return 1;
  }
}
  if(outFile.empty()){outFile = "result.root";}
  if (outType.empty()){outType = "bunch";}
  if (version.empty()){version = "Cal";}
  if (version!="Cal"){version = "Cal";} // THE CODE IS  AT THE MONENT ONLY WITH CAL FUNCTIONAL
  if (beamline.empty()){beamline = "off";}

  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  long seeds[2];
  time_t systime = time(NULL);
  seeds[0] = (long) systime;
  seeds[1] = (long) (systime*G4UniformRand());
  G4Random::setTheSeeds(seeds);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  DetectorConstruction* det;
  PrimaryGeneratorAction* prim;
  runManager->SetUserInitialization(det = new DetectorConstruction(version,beamline));

  // PhysicsList FTFP_BERT with optical PhysList (like in all optical physics examples for Geant4)
  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  // physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  // G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  // physicsList->RegisterPhysics(opticalPhysics);
  // runManager-> SetUserInitialization(physicsList);

  // Physicslist user defined
  runManager->SetUserInitialization(new PhysicsList(version));
  runManager->SetUserAction(prim = new PrimaryGeneratorAction());


  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();


  // set user action classes
  RunAction* run;
  EventAction* event;
  runManager->SetUserAction(run = new RunAction(det,prim,outFile,outType,version));
  runManager->SetUserAction(event = new EventAction(run,outType,version));
  runManager->SetUserAction(new SteppingAction(det,event,run,outType,version));

  // get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    //
    if ( macro.size() ) {
      // batch mode
      G4String command = "/control/execute ";
      UImanager->ApplyCommand(command+macro);
    }
    else  {
      // interactive mode : define UI session
      UImanager->ApplyCommand("/control/execute init_vis.mac");
      if (ui->IsGUI()) {
        UImanager->ApplyCommand("/control/execute gui.mac");
      }
      ui->SessionStart();
      delete ui;
    }

  // job termination

  delete visManager;

  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
