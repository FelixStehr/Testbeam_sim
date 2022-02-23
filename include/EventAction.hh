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
/// \file polarisation/Pol01/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 98772 2016-08-09 14:25:31Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "RunAction.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction*, G4String outType, G4String version);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void AddVals(G4double Eval, G4double Npart);
  void AddGammaVals(G4double EGammapol, G4double Ngamma);
  void AddEnergyCalo(G4double ECalo);
  void AddPhotonEnergy(G4double EPhoton);
  void AddGammaEnergy(G4double EGamma);

  // here again just the special case for the two cristals
  void AddEcalo0(G4double Ec0);
  void AddEelectron0(G4double Ee0);
  void AddEgamma0(G4double Eg0);
  void AddErest0(G4double Er0);

  void AddEcalo1(G4double Ec1);
  void AddEelectron1(G4double Ee1);
  void AddEgamma1(G4double Eg1);
  void AddErest1(G4double Er1);

  // void AddNelectron0(G4double Ne0);
  // void AddNelectron1(G4double Ne1);
  // void AddNgamma0(G4double Ng0);
  // void AddNgamma1(G4double Ng1);
  // void AddNrest0(G4double Nr0);
  // void AddNrest1(G4double Nr1);



private:
  RunAction* fRunAction;

  G4double fEnergySum;
  G4double fNP;
  G4double fGammaEnergySum;
  G4double fNGamma;
  G4double fEnergyCalo;
  G4double fPhotonEnergySum;
  G4double fGammaEnergyIn;
  G4String outputType;
  G4String versionType;

  //here 2 Chrystals CP 0 and 1
  G4double Ecalo0;
  G4double Eelectron0;
  G4double Egamma0;
  G4double Erest0;

  G4double Ecalo1;
  G4double Eelectron1;
  G4double Egamma1;
  G4double Erest1;


  G4double Nelectron0;
  G4double Ngamma0;
  G4double Nrest0;

  G4double Nelectron1;
  G4double Ngamma1;
  G4double Nrest1;





};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// inline functions
inline void EventAction::AddVals(G4double Eval, G4double Npart) {
  fEnergySum += Eval;
  fNP += Npart;
}
inline void EventAction::AddGammaVals(G4double EGammapol, G4double Ngamma) {
  fGammaEnergySum += EGammapol;
  fNGamma += Ngamma;
}
inline void EventAction::AddEnergyCalo(G4double ECalo) { //deposited energy in the crstals
  fEnergyCalo += ECalo;
}

inline void EventAction::AddPhotonEnergy(G4double EPhoton) { //Photon Energy which comes out of the crystal
  fPhotonEnergySum += EPhoton;
}

inline void EventAction::AddGammaEnergy(G4double EGamma) { //Gamma Energy on the entrace of the crystals
  fGammaEnergyIn+= EGamma;
}


// Here special case of two crystals

// Crystal 0
inline void EventAction::AddEcalo0(G4double Ec0) { //Gamma Energy on the entrace of the crystals
  Ecalo0+= Ec0;
}

inline void EventAction::AddEelectron0(G4double Ee0) { //Gamma Energy on the entrace of the crystals
  Eelectron0+= Ee0;
  Nelectron0+= 1.;
}

inline void EventAction::AddEgamma0(G4double Eg0) { //Gamma Energy on the entrace of the crystals
  Egamma0+= Eg0;
  Ngamma0+= 1.;
}

inline void EventAction::AddErest0(G4double Er0) { //Gamma Energy on the entrace of the crystals
  Erest0+= Er0;
  Nrest0+= 1.;
}


// Crystal 1
inline void EventAction::AddEcalo1(G4double Ec1) { //Gamma Energy on the entrace of the crystals
  Ecalo1+= Ec1;
}

inline void EventAction::AddEelectron1(G4double Ee1) { //Gamma Energy on the entrace of the crystals
  Eelectron1+= Ee1;
  Nelectron1+= 1.;
}

inline void EventAction::AddEgamma1(G4double Eg1) { //Gamma Energy on the entrace of the crystals
  Egamma1+= Eg1;
  Ngamma1+=1.;
}

inline void EventAction::AddErest1(G4double Er1) { //Gamma Energy on the entrace of the crystals
  Erest1+= Er1;
  Nrest1+= 1.;
}


#endif
