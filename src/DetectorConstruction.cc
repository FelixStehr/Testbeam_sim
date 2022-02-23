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
/// \file polarisation/Pol01/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 98772 2016-08-09 14:25:31Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"

#include "G4PolarizationManager.hh"
// #include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4RotationMatrix.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4String version, G4String beamline)
: G4VUserDetectorConstruction(),
  PhysicalWorld(0), PhysicalCore(0), fConvMaterial(0), fWorldMaterial(0), fCaloMaterial(0)
{
  versionType=version;
  beamlineStatus=beamline;
  allMaterials = new Materials();
  allMaterials->DefineMaterials();

  fSizeXY = 50*mm;
  fCoreThick = 75*mm;
  fConvThick = 1.75*mm;
  fWorldSize = 9*m;
  CrystalNumber= "two";
  SFStatus ="false";
  LanexStatus ="false";
  StrawStatus ="false";
  dCalo = 10*cm;
  RCollimator = 2.5*mm;
  CaloXpos = 0.*mm;
  // G4cout <<"TEST TEST TEST"<< G4endl;


  SetConvMaterial("G4_W");
  // if(beamlineStatus=="on"){
    // SetWorldMaterial("Air");}
  // else if(beamlineStatus=="off"){
    // SetWorldMaterial("Galactic");}
  // else {
  //   SetWorldMaterial("Galactic");
  //   G4cout << "NO VALID IMPUT FOR BEAMLINE: beamlineStatus set to off!->all simulations in vacuum"<< G4endl;
  //    }
  SetCaloMaterial("TF101");
  SetWorldMaterial("Air_NoRI");

  fMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4PolarizationManager::GetInstance()->Clean();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Geometry parameters
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //
  // Polarimeter
  //
  G4double  maggap1 = 48.5*mm;
  G4double  maggap2 = 12.5*mm;
  G4double  absrad  = fSizeXY/2.;
  G4double shieldrad=75.0*mm;
  G4double vacthick = 1.0*mm;//4.*mm;
  G4double corethick = fCoreThick;
  G4double convthick = fConvThick;
  G4double  coilthick = corethick + 25.0*mm;
  G4double  shieldthick = corethick - 25.0*mm;
  G4double  conedist = corethick/2. + maggap2;
  G4double  magthick = 2.*(maggap1+convthick+maggap2)+corethick;

  //
  // Calorimeter
  //

  G4double detthick = 45.*cm; // crystal size
  G4double detx = 3.8 *cm;
  G4double dety = 3.8 *cm;


  G4double AirgabCrystalPMTthick =1.*mm;//4.*mm;
  G4double AirgabCrystalPMTxy= 27.*mm; // entrance size of the PMT mount
  G4double PMTGlassthick=0.5 *mm;//4.*mm;
  G4double PMTGlassxy   = 25. *mm;
  G4double DetVacstepthick= 1. *mm; //4.*mm;// here VacStep as detector volume later kathod material
  G4double DetVacstepxy= PMTGlassxy; //18
  G4double PMTmountlength = AirgabCrystalPMTthick+PMTGlassthick+DetVacstepthick;
  G4double PMTmountxy = detx;

  G4double alairgapthick = 0.001 *mm; //4.*mm;   // thickness of the air gap between the aluwrapping and the crystal
  G4double alairgapx = detx + 2*alairgapthick;
  G4double alairgapy = dety + 2*alairgapthick;
  G4double alairgaplength = detthick + alairgapthick + PMTmountlength;

  G4double aluwrapthick =0.01  *mm;//4.*mm;  // wikipedia: alu foil thickness between 0.004 and 0.02 mm
  G4double aluwrapx = alairgapx + 2*aluwrapthick;
  G4double aluwrapy = alairgapy + 2*aluwrapthick;
  G4double aluwraplength = alairgaplength + aluwrapthick;

  G4double venylfoilthick = 0.1 *mm;//4.*mm; // just a guess becouse we have more than 1 layer
  G4double venylfoilxy = aluwrapx + 2* venylfoilthick;
  G4double venylfoillength = aluwraplength + venylfoilthick;

  //defining the size of the Calorimeterzell and the virtual calorimeter (mother volume of the calorimetercells)
  //G4int NbofCalor = 9; //here later free paramter to select numer of crystals
  //G4double calorcellxy = aluwrapx;

  G4double calorcellxy;
  if (venylfoilxy<= 4*cm){
    calorcellxy= 4*cm;}// here I took 4cm becouse in the experiment wie assumed that the crystal centres are 4cm apart
  else{ calorcellxy=venylfoilxy;}
  G4double calorcelllength = venylfoillength + vacthick;
  G4double virtcalorxy;

  if (CrystalNumber == "one"){
    virtcalorxy = calorcellxy;
  }
  else if (CrystalNumber == "nine" || CrystalNumber == "four" ){
    virtcalorxy = 3*calorcellxy;
  }
  else if (CrystalNumber == "two"){
    virtcalorxy = 2*calorcellxy;
  }
  else {
    virtcalorxy = calorcellxy;
    SetCrystalnumber("one");
    G4cout << "NO VALID IMPUT FOR CRYSTALNUMBER: Crystalnumber set to one!"<< G4endl;
  }


  G4double virtcalorlength = calorcelllength;
  // G4double spacePolCal = dCalo;//50. *mm;
  G4double caloZposition;

  if(beamlineStatus=="on"){
    caloZposition = virtcalorlength/2+dCalo;}
  else {caloZposition = (magthick+virtcalorlength)/2+dCalo;}

  //
  //Scintillator fingers dimensions
  //

  G4double SFlength = 60*mm; // Scintillator dimensions just roughly measured
  G4double SFx = 4* mm;
  G4double SFy = 13* mm;

  G4double SFairgapthick = 0.001 *mm;
  G4double SFairgapx = SFx + 2*SFairgapthick;
  G4double SFairgapy = SFy + 2*SFairgapthick;
  G4double SFairgaplength = SFlength + SFairgapthick;

  G4double SFAluthick = 0.01 *mm;
  G4double SFAlux =  SFairgapx + 2*SFAluthick   ;
  G4double SFAluy = SFairgapy + 2*SFAluthick    ;
  G4double SFAlulength = SFairgaplength + SFAluthick;

  G4double SFVacstepthick = 1 *mm;
  G4double SFVacstepx =  SFAlux;
  G4double SFVacstepy =  SFAluy;

  G4double SFVinylthick = 0.1 *mm;
  G4double SFVinylx= SFAlux + 2*SFVinylthick  ;
  G4double SFVinyly= SFAluy + 2*SFVinylthick  ;
  G4double SFVinyllength = SFAlulength + SFVinylthick + SFVacstepthick;

  G4double SFpositonZ = 420*mm - SFVinylx/2;
  G4double Intersectiondis = 10*mm;     // Position were the electron beam hits the SF from its edge
  // G4double SFpositonx = -(SFVinyllength/2.-Intersectiondis)/sqrt(2);
  // G4double SFpositony = SFpositonx;
  G4double SFpositonx = 400*mm;
  G4double SFpositony = 0*mm;

  G4RotationMatrix* SFRotation = new G4RotationMatrix();
  SFRotation->rotateY(90.*deg);
  SFRotation->rotateX(-45.*deg);
  // SFRotation->rotateZ(45.*deg);

  //
  // Lanex screen (from Jon)
  //
  G4double ToyX = 400*mm;

  G4double scz1 = 0.25 * mm;
  G4double scz2 = 0.14 * mm;
  G4double scz3 = 0.006 * mm;
  G4double scx = 10 * cm;
  G4double scy = 10 * cm;
  G4double ZposLanex = dCalo - 5*mm - (scz1+scz2+scz3)/2;
  //
  // Cherencov Straw from Luxe
  //
  G4double ZposStraw = 28* cm;
  G4double LStraw = 200* mm;
  G4double dMylar = 36E-6*m;
  G4double dCuAu = 50E-9*m;
  G4double rStraw = 12/2 *mm;
  G4double rMylar = 10/2 *mm;
  G4double rCuAu = rMylar - dMylar;
  G4double rOil = rCuAu - dCuAu;

  G4RotationMatrix* StrawRotation = new G4RotationMatrix();
  StrawRotation->rotateX(90.*deg);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Get materials
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // G4Material* absMat    = allMaterials->GetMat("G4_W");
  G4Material* absMat    = allMaterials->GetMat("Galactic");
  G4Material* magMat    = allMaterials->GetMat("G4_Fe");
  G4Material* coilMat   = allMaterials->GetMat("G4_Cu");
  G4Material* shieldMat = allMaterials->GetMat("G4_Pb");
  G4Material* Air       = allMaterials->GetMat("Air");
  G4Material* Air_NoRI = allMaterials->GetMat("Air_NoRI");
  G4Material* Al        = allMaterials->GetMat("Aluminium");
  G4Material* Vacuum    = allMaterials->GetMat("Galactic");
  G4Material* Concrete  = allMaterials->GetMat("G4_CONCRETE");
  G4Material* StainlessSteel  = allMaterials->GetMat("G4_STAINLESS-STEEL");
  G4Material* Copper = allMaterials->GetMat("G4_Cu");
  G4Material* Venyl = allMaterials->GetMat("G4_POLYVINYL_CHLORIDE");
  G4Material* PMTGlass = allMaterials->GetMat("Glass");
  //change the SFmat
//  G4Material* SFMat = allMaterials->GetMat("G4_PLASTIC_SC_VINYLTOLUENE"); // not sure if it is the right material
  G4Material* SFMat = allMaterials->GetMat("G4_PARAFFIN");
  G4Material* lanex = allMaterials->GetMat("G4_GADOLINIUM_OXYSULFIDE");
  G4Material* Pstyrene = allMaterials->GetMat("Polystyrene");
  G4Material* CuAu = allMaterials->GetMat("GoldAlloy");
  G4Material* Mylar = allMaterials->GetMat("G4_MYLAR");
  G4Material* StrawMat = allMaterials->GetMat("PLA");
  // G4Material* Parafin = allMaterials->GetMat("G4_PARAFFIN");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //VisAttributes
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //Polarimeter
  //

    G4VisAttributes * MagnetVis= new G4VisAttributes( G4Colour(255/255. ,102/255. ,102/255. ));
    MagnetVis->SetVisibility(true);
    MagnetVis->SetLineWidth(1);

    G4VisAttributes * CopperCoilVis= new G4VisAttributes( G4Colour(255/255. ,0/255. ,255/255. ));
    CopperCoilVis->SetVisibility(true);
    CopperCoilVis->SetLineWidth(1);

    G4VisAttributes * LeadTubeVis= new G4VisAttributes( G4Colour(0/255. ,102/255. ,204/255. ));
    LeadTubeVis->SetVisibility(true);
    LeadTubeVis->SetLineWidth(1);

    G4VisAttributes * ConversionTargetVis= new G4VisAttributes( G4Colour(105/255. ,105/255. ,105/255. ));
    ConversionTargetVis->SetVisibility(true);
    ConversionTargetVis->SetLineWidth(2);
    ConversionTargetVis->SetForceSolid(true);

    G4VisAttributes * IronCoreVis= new G4VisAttributes( G4Colour(51/255. ,51/255. ,255/255. ));
    IronCoreVis->SetVisibility(true);
    IronCoreVis->SetLineWidth(2);
    IronCoreVis->SetForceSolid(true);

    G4VisAttributes * VacStepVis= new G4VisAttributes( G4Colour(255/255. ,165/255. ,0/255. ));
    VacStepVis->SetVisibility(true);
    VacStepVis->SetLineWidth(1);
    VacStepVis->SetForceSolid(true);


  //Calorimeter
  //
    G4VisAttributes * AirVis= new G4VisAttributes( G4Colour(0/255. ,255/255. ,255/255. ));
    AirVis->SetVisibility(true);
    AirVis->SetLineWidth(2);
    // AirVis->SetForceWireframe( true );
    AirVis->SetForceSolid(true);


    G4VisAttributes * AluVis= new G4VisAttributes( G4Colour(87/255. ,87/255. ,87/255. ));
    AluVis->SetVisibility(true);
    AluVis->SetLineWidth(2);
    AluVis->SetForceSolid(true);

    G4VisAttributes * CrystalVis= new G4VisAttributes( G4Colour(224/255. ,255/255. ,255/255. ));
    // G4VisAttributes * CrystalVis= new G4VisAttributes( G4Colour(255/255. ,193/255. ,193/255. ));
    CrystalVis->SetVisibility(true);
    CrystalVis->SetLineWidth(2);
    CrystalVis->SetForceSolid(true);

    G4VisAttributes * VacStep3Vis= new G4VisAttributes( G4Colour(255/255. ,165/255. ,0/255. ));
    VacStep3Vis->SetVisibility(true);
    VacStep3Vis->SetLineWidth(1);
    VacStep3Vis->SetForceSolid(true);

    G4VisAttributes * PMTmountVis= new G4VisAttributes( G4Colour(255/255. ,0/255. ,0/255., 0.95));
    PMTmountVis->SetVisibility(false);
    PMTmountVis->SetLineWidth(1);
    //PMTmountVis->SetForceSolid(true);


    G4VisAttributes * VenyfoilVis= new G4VisAttributes( G4Colour(255/255. ,0/255. ,255/255. ));
    VenyfoilVis->SetVisibility(true);
    VenyfoilVis->SetLineWidth(1);
    VenyfoilVis->SetForceSolid(true);

    G4VisAttributes * GlassVis= new G4VisAttributes( G4Colour(0/255. ,0/255. ,205/255. ));
    GlassVis->SetVisibility(true);
    GlassVis->SetLineWidth(1);
    GlassVis->SetForceSolid(true);




  //
  //Beam Line
  //
    G4VisAttributes * ConcreteVis= new G4VisAttributes( G4Colour(119/255. ,136/255. ,153/255., 0.95 ));
    ConcreteVis->SetVisibility(true);
    ConcreteVis->SetLineWidth(2);
    ConcreteVis->SetForceWireframe( true );
    //ConcreteVis->SetForceSolid(true);

    G4VisAttributes * ConcreteHoleVis= new G4VisAttributes( G4Colour(119/255. ,136/255. ,153/255. ));
    ConcreteHoleVis->SetVisibility(true);
    ConcreteHoleVis->SetLineWidth(2);
    //ConcreteHoleVis->SetForceSolid(false);
    ConcreteHoleVis->SetForceWireframe( true );

    G4VisAttributes * BeamTubeVis= new G4VisAttributes( G4Colour(119/255. ,136/255. ,153/255., 0.95 ));
    BeamTubeVis->SetVisibility(true);
    BeamTubeVis->SetLineWidth(1);
    BeamTubeVis->SetForceSolid(true);

    G4VisAttributes * AluWindowVis= new G4VisAttributes( G4Colour(255/255. ,165/255. ,0/255.));
    AluWindowVis->SetVisibility(true);
    AluWindowVis->SetLineWidth(2);
    AluWindowVis->SetForceSolid(true);

    G4VisAttributes * CollimatorVis= new G4VisAttributes( G4Colour(205/255. ,55/255. ,0/255.));
    CollimatorVis->SetVisibility(true);
    CollimatorVis->SetLineWidth(1);
    CollimatorVis->SetForceSolid(true);

    G4VisAttributes * BeamLineVacVis= new G4VisAttributes( G4Colour(142/255. ,229/255. ,238/255.));
    CollimatorVis->SetVisibility(true);
    CollimatorVis->SetLineWidth(1);
    CollimatorVis->SetForceSolid(true);

    G4VisAttributes * BleiVis= new G4VisAttributes( G4Colour(238/255. ,238/255. ,0/255., 0.4));
    BleiVis->SetVisibility(true);
    BleiVis->SetLineWidth(1);
  //  BleiVis->SetForceSolid(true);


  //
  // Luxe cherencov straw and Lanex scintillator screen
  //

    G4VisAttributes * LanexVis= new G4VisAttributes( G4Colour(0/255. ,255/255. ,0/255.));
    LanexVis->SetVisibility(true);
    LanexVis->SetLineWidth(1);

    G4VisAttributes * PolyVis1= new G4VisAttributes( G4Colour(173/255. ,216/255. ,230/255.));
    PolyVis1->SetVisibility(true);
    PolyVis1->SetLineWidth(1);

    G4VisAttributes * PolyVis2= new G4VisAttributes( G4Colour(135/255. ,206/255. ,235/255.));
    PolyVis2->SetVisibility(true);
    PolyVis2->SetLineWidth(1);

    G4VisAttributes * StrawVis= new G4VisAttributes( G4Colour(0/255. ,0/255. ,0/255.));
    StrawVis->SetVisibility(true);
    StrawVis->SetLineWidth(1);



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // World
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  G4Box*
  SolidWorld = new G4Box("World",                            //name
                   fWorldSize/2,fWorldSize/2,fWorldSize/2); //dimensions

  G4LogicalVolume*
  LogicalWorld = new G4LogicalVolume(SolidWorld,                   //shape
                               fWorldMaterial,           //material
                              "World");                  //name

  PhysicalWorld = new G4PVPlacement(0,                          //no rotation
                             G4ThreeVector(),            //at (0,0,0)
                             LogicalWorld,                     //logical volume
                             "World",                    //name
                             0,                          //mother volume
                             false,                      //no boolean operation
                             0);                         //copy number


  LogicalWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Polarimeter geometry
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (versionType == "Pol" || versionType == "PolCal"){
  // Magnet
  //
  G4double DzArrayMagnet   [] = {-magthick/2.  , -conedist , -corethick/2. , corethick/2.,  conedist, magthick/2.    };
  G4double RminArrayMagnet [] = {36.84308*mm,  absrad,  absrad , absrad,  absrad,  36.84308*mm};
  G4double RmaxArrayMagnet [] = {196.0*mm   ,  196.0*mm, 196.0*mm ,196.0*mm, 196.0*mm, 196.0*mm    };

  G4Polycone *solidMagnet = new G4Polycone("solidMagnet", 	 //its name
            0.0*deg, 		 //its start angle
            360.0*deg,		 //its opening angle
            6, 		         //its nZ
            DzArrayMagnet, 	 //z value
            RminArrayMagnet, 	 //rmin
            RmaxArrayMagnet ); 	 //rmax

  G4LogicalVolume * LogicalMagnet = new G4LogicalVolume(solidMagnet, //its solid
                  magMat, 	 //its material
                   "Magnet" ,		 //its name
                   0,0,0);

  G4VPhysicalVolume *  PhysicalMagnet= new G4PVPlacement(0,	//rotation
  								 G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm),// translation position
  								 LogicalMagnet,      //its logical volume
  						       "PhysicalMagnet",   //its name  (2nd constructor)
  						       LogicalWorld,         //its mother volume
  						       false,              //no boolean operation
  						       0);                 //copy number

  LogicalMagnet->SetVisAttributes(MagnetVis);

  //Copper Coils
  //
  G4Tubs *solidCuTube= new G4Tubs("solidCuTube", //name
                                  shieldrad, // inner radius
                                  170.0*mm,  // outer radius
                                  coilthick/2., // half length in z
                                  0.0*deg,  // starting angle
                                  360.0*deg ); // total angle

  G4LogicalVolume * LogicalCuTube = new G4LogicalVolume(solidCuTube, //its solid
             coilMat,              //its material
             "CuTube" ,		         //its name
             0,0,0);

  G4VPhysicalVolume *PhysicalCuTube= new G4PVPlacement(0,	//rotation
         G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm), // translation position
           LogicalCuTube,      // its logical volume
           "PhysicalCuTube",   // its name  (2nd constructor)
           LogicalMagnet,     // its mother volume
           false,              // no boolean operation
           0);

  LogicalCuTube->SetVisAttributes(CopperCoilVis);

  // Lead Tube
  //
  G4Tubs *solidPbtube= new G4Tubs("solidPbtube", // name
                                  absrad, // inner radius
                                  shieldrad, // outer radius
                                  shieldthick/2., // half length in z
                                  0.0*deg, // start angle
                                  360.0*deg ); // total angle

  G4LogicalVolume * LogicalPbtube = new G4LogicalVolume(solidPbtube, 	 //its solid
  						shieldMat, 		 //its material
  						"Pbtube" ,		 //its name
  						0,0,0);

  G4VPhysicalVolume *PhysicalPbTube= new G4PVPlacement(0,	//rotation
  				G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm),
  				LogicalPbtube,      //its logical volume
  			    "PhysicalPbTube",   //its name  (2nd constructor)
  			    LogicalMagnet,     //its mother volume
  			    false,              //no boolean operation
  			    0);                 //copy number

  LogicalPbtube->SetVisAttributes(LeadTubeVis);

  // Reconversion Target
  //
  auto solidReconversion = new G4Tubs("solidReconversion", // name
                                      0.0*mm, // inner radius
                                      absrad, // outer radius
                                      convthick/2, // half length in z
                                      0.0*deg, // starting angle
                                      360.0*deg ); // total angle

  auto LogicalReconversion = new G4LogicalVolume(solidReconversion, 	 //its solid
              absMat,          //its material
              "solidReconversion" ,	 //its name
              0,0,0);

  auto PhysicalReconversion= new G4PVPlacement(0,	//rotation
              G4ThreeVector(0.0*mm, 0.0*mm, -maggap2-convthick/2-corethick/2),
            LogicalReconversion,         //its logical volume
            "PhysicalReconversion",   //its name  (2nd constructor)
            LogicalWorld,              //its mother volume
            false,                 //no boolean operation
            0);                       //copy number

  LogicalReconversion->SetVisAttributes(ConversionTargetVis);

  // Iron Core
  //
  auto solidCore = new G4Tubs ("Container",                           //its name
                   0.0*mm, absrad*mm, corethick/2, 0.0*deg, 360.0*deg );//its dimensions

  G4LogicalVolume*
  LogicalCore = new G4LogicalVolume(solidCore,                        //its shape
                             magMat,             //its material
                             "IronCore");                   //its name

  PhysicalCore = new G4PVPlacement(0,                             //no rotation
                           G4ThreeVector(),               //at (0,0,0)
                           LogicalCore,                          //its logical volume
                           "IronCorePV",    //its name
                           LogicalWorld,                        //its mother  volume
                           false,                         //no boolean operation
                           0);                            //copy number

  // register logical Volume in PolarizationManager with polarization
  G4PolarizationManager * polMgr = G4PolarizationManager::GetInstance();
  polMgr->SetVolumePolarization(LogicalCore,G4ThreeVector(0.,0.,1.));

  LogicalCore->SetVisAttributes(IronCoreVis);

  //
  //vacuum step 1
  //
  auto VacStepS1 = new G4Tubs("VacStep1",  //Name
                              0.,         // inner radius
                              absrad,     // outer radius
                              vacthick/2., // half length in z
                              0.0*deg,    // starting phi angle
                              360.0*deg); // angle of the segment


  auto  VacStepLV1 = new G4LogicalVolume(VacStepS1,    //its solid
                                         fWorldMaterial,    //its material
                                         "VacStep1");  //its name

  fVacStepPV1 = new G4PVPlacement(0,                 //no rotation
                       G4ThreeVector(0.,0., - corethick/2 -maggap2 +vacthick/2 +1.0*mm),    //its position
                               VacStepLV1,            //its logical volume
                               "VacStep1",                 //its name
                               LogicalWorld,               //its mother
                               false,                     //no boolean operat
                               0);                        //copy number

  VacStepLV1->SetVisAttributes(VacStepVis);

  // vacuum step 2
  //
  auto VacStepS2 = new G4Tubs("VacStep2",  //Name
                               0.,         // inner radius
                               absrad,     // outer radius
                               vacthick/2., // half length in z
                               0.0*deg,    // starting phi angle
                               360.0*deg); // angle of the segment

  auto VacStepLV2 = new G4LogicalVolume(VacStepS2,    //its solid
                                        fWorldMaterial,    //its material
                                        "VacStep1");       //its name

  fVacStepPV2 = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,corethick/2 + vacthick/2 + 10.0*mm),    //its position
                                VacStepLV2,            //its logical volume
                                "VacStep2",                 //its name
                                LogicalWorld,               //its mother
                                false,                     //no boolean operat
                                0);                       //copy number

  VacStepLV2->SetVisAttributes(VacStepVis);

  } //end if-statement polarimeter


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calorimeter geometrys
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (versionType == "Cal" || versionType == "PolCal"){

  // G4RotationMatrix* myRotation = new G4RotationMatrix();
  // myRotation->rotateX(90.*deg);
  // myRotation->rotateY(0.*deg);
  // myRotation->rotateZ(0.*deg);

  // Virtuel calorimeter (mother volume for the hole calorimeter/detector)
  auto fVirtCaloS= new G4Box("virtualCalorimeterS",  //Name
                                virtcalorxy/2.,   // x size
                                virtcalorxy/2.,     // y size
                                virtcalorlength/2.); // z size


  auto fVirtCaloLV = new G4LogicalVolume(fVirtCaloS,    //its solid
                                         fWorldMaterial,    //its material
                                         "virtualCalorimeterLV");       //its name

  fVirtCaloPV = new G4PVPlacement(0,                   //no rotation: 0 or rotation 90 deg around x axis : myRotaion
                         G4ThreeVector(calorcellxy/2+CaloXpos,0.,caloZposition),    //if CaloXpos =0 the crystal with the square PMT is in the centre of the beam
                                 fVirtCaloLV,            //its logical volume
                                 "virtualCalorimeterPV",                 //its name
                                LogicalWorld,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

  fVirtCaloLV->SetVisAttributes(G4VisAttributes::GetInvisible());



  // Calorimeter cells (placing (for now) 9 calorimeter cells in the virtual calorimeter)
  //
  auto fCaloCellS= new G4Box("physicalcalorimeterS",  //Name
                               calorcellxy/2.,   // x size
                               calorcellxy/2.,     // y size
                               calorcelllength/2.); // z size

  auto fCaloCellLV = new G4LogicalVolume(fCaloCellS,    //its solid
                                        fWorldMaterial,    //its material
                                        "physicalcalorimeterLV");       //its name

  if(CrystalNumber=="nine"){
  //the array for the placement of the 9 calorimetercells in the virtual calorimeter
  G4double CalorRX[9]={0,0,calorcellxy,calorcellxy,calorcellxy,0,-calorcellxy,-calorcellxy,-calorcellxy};
  G4double CalorRY[9]={0,calorcellxy,calorcellxy,0,-calorcellxy,-calorcellxy,-calorcellxy,0,calorcellxy};

  for (G4int i=0;i<=8;i++){
  fCaloCellPV = new G4PVPlacement(0,		       //no rotation
               G4ThreeVector(CalorRX[i],CalorRY[i],0),  //its position
               fCaloCellLV,            //its logical volume
              "physicalcalorimeterPV",    //its name
               fVirtCaloLV,               //its mother
               false,                     //no boolean operat
               i);                        //copy number       //copy number
  }
  }

  else if (CrystalNumber == "one"){
  fCaloCellPV = new G4PVPlacement(0,		       //no rotation
               G4ThreeVector(0,0,0),  //its position
               fCaloCellLV,            //its logical volume
              "physicalcalorimeterPV",    //its name
               fVirtCaloLV,               //its mother
               false,                     //no boolean operat
               0);                        //copy number       //copy number
  }

  else if(CrystalNumber=="four"){
  //the array for the placement of the 4 calorimetercells in the virtual calorimeter
  G4double CalorRX4[4]={0,calorcellxy,-calorcellxy,0};
  G4double CalorRY4[4]={0,0,0,calorcellxy};
  for (G4int i=0;i<=3;i++){
  fCaloCellPV = new G4PVPlacement(0,		       //no rotation
               G4ThreeVector(CalorRX4[i],CalorRY4[i],0),  //its position
               fCaloCellLV,            //its logical volume
              "physicalcalorimeterPV",    //its name
               fVirtCaloLV,               //its mother
               false,                     //no boolean operat
               i);                        //copy number       //copy number
  }
  }

  else if(CrystalNumber=="two"){
  //the array for the placement of the 2 calorimetercells in the virtual calorimeter
  G4double CalorRX4[2]={-calorcellxy/2,calorcellxy/2};
  G4double CalorRY4[2]={0,0};
  for (G4int i=0;i<=1;i++){
  fCaloCellPV = new G4PVPlacement(0,		       //no rotation
               G4ThreeVector(CalorRX4[i],CalorRY4[i],0),  //its position
               fCaloCellLV,            //its logical volume
              "physicalcalorimeterPV",    //its name
               fVirtCaloLV,               //its mother
               false,                     //no boolean operat
               i);                        //copy number       //copy number
  }
  }


  fCaloCellLV->SetVisAttributes(G4VisAttributes::GetInvisible());





  auto fVenylS= new G4Box("VenylfoilS",  //Name
                                venylfoilxy/2.,   // x size
                                venylfoilxy/2.,     // y size
                                venylfoillength/2.); // z size


  auto fVenylLV = new G4LogicalVolume(fVenylS,    //its solid
                                         Venyl,    //its material
                                         "VenylfoilLV");       //its name

  fVenylPV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,vacthick/2),    //its position // old 0.,0.,-vacthick/2
                                 fVenylLV,            //its logical volume
                                 "VenylfoilPV",                 //its name
                                 fCaloCellLV,               //its mother
                                 false,                     //no boolean operat
                                     0);                        //copy number

  fVenylLV->SetVisAttributes(VenyfoilVis);

  // Alu-wrapping
  //
  auto fAluWrapS= new G4Box("AluWrappingS",  //Name
                                aluwrapx/2.,   // x size
                                aluwrapy/2.,     // y size
                                aluwraplength/2.); // z size


  auto fAluwrapLV = new G4LogicalVolume(fAluWrapS,    //its solid
                                         Al,    //its material
                                         "AluWrappingLV");       //its name

  fAluwrapPV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,venylfoilthick/2),    // its positon
                                 fAluwrapLV,            //its logical volume
                                 "AluWrappingPV",                 //its name
                                 fVenylLV,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

  fAluwrapLV->SetVisAttributes(AluVis);

  //
  // AirGap
  //
  auto fAlAirGapS= new G4Box("AlAirGapS",  //Name
                               alairgapx/2.,   // x size
                               alairgapy/2.,     // y size
                               alairgaplength/2.); // z size


  auto fAlAirGapLV = new G4LogicalVolume(fAlAirGapS,    //its solid
                                        Air,    //its material
                                        "AlAirGapLV");       //its name

  fAlAirGapPV = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,aluwrapthick/2),    //its position
                                fAlAirGapLV,            //its logical volume
                                "AlAirGapPV",                 //its name
                                fAluwrapLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  // fAlAirGapLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  fAlAirGapLV->SetVisAttributes(AirVis);

  //
  // Crystal (aka Detector)
  //
  auto fDetectorS= new G4Box("DetectorS",  //Name
                               detx/2.,   // x size
                               dety/2.,     // y size
                               detthick/2.); // z size


  fDetectorLV = new G4LogicalVolume(fDetectorS,    //its solid
                                        fCaloMaterial,    //its material
                                        "DetectorLV");       //its name

  fDetectorPV = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,(alairgapthick-PMTmountlength)/2),    //its position
                                fDetectorLV,            //its logical volume
                                "DetectorPV",                 //its name
                                fAlAirGapLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  fDetectorLV->SetVisAttributes(CrystalVis);

  //
  // PMTmount
  //
  auto fPMTmountS= new G4Box("PMTmountS",  //Name
                               alairgapx/2.,  //PMTmountxy/2.,   // x size
                               alairgapx/2.,  //PMTmountxy/2.,     // y size
                               PMTmountlength/2.); // z size


  auto fPMTmountLV = new G4LogicalVolume(fPMTmountS,    //its solid
                                        Venyl,    //its material
                                        "PMTmountLV");       //its name

  auto fPMTmountPV = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,(alairgapthick+detthick)/2),    //its position
                                fPMTmountLV,            //its logical volume
                                "PMTmountPV",                 //its name
                                fAlAirGapLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  fPMTmountLV->SetVisAttributes(PMTmountVis);

  //
  // Space of Air between the Crystal and the PMT
  //

  auto fAirgabCrystalPMTS= new G4Box("AirgabCrystalPMTS",  //Name
                               detx/2,//AirgabCrystalPMTxy/2.,   // x size
                               detx/2,//AirgabCrystalPMTxy/2.,     // y size
                               AirgabCrystalPMTthick/2.); // z size


  auto fAirgabCrystalPMTLV = new G4LogicalVolume(fAirgabCrystalPMTS,    //its solid
                                        Air,    //its material
                                        "AirgabCrystalPMTLV");       //its name

  auto fAirgabCrystalPMTPV = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,-(PMTGlassthick+DetVacstepthick)/2),    //its position
                                fAirgabCrystalPMTLV,            //its logical volume
                                "AirgabCrystalPMTPV",                 //its name
                                fPMTmountLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number
  fAirgabCrystalPMTLV->SetVisAttributes(AirVis);

  //
  // PMTGlass Window
  //

  auto fPMTGlassS= new G4Box("PMTGlassS",  //Name
                                detx/2,//PMTGlassxy/2.,   // x size
                                detx/2,//PMTGlassxy/2.,     // y size
                               PMTGlassthick/2.); // z size


  auto fPMTGlassLV = new G4LogicalVolume(fPMTGlassS,    //its solid
                                        PMTGlass,    //its material
                                        "PMTGlassLV");       //its name

  auto fPMTGlassPV = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,(AirgabCrystalPMTthick-DetVacstepthick)/2),    //its position
                                fPMTGlassLV,            //its logical volume
                                "PMTGlassPV",                 //its name
                                fPMTmountLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  fPMTGlassLV->SetVisAttributes(GlassVis);

  //
  //vacuum step 3 this is the detector behind the PMTGLass // use 4 anodes
  auto fVacStepS3 = new G4Box("VacStep3S",  //Name
                               detx/2,//DetVacstepxy/2.,
                               detx/2,// DetVacstepxy/2,
                               DetVacstepthick/2.);

  auto fVacStepLV3 = new G4LogicalVolume(fVacStepS3,    //its solid
                                        Vacuum,    //its material
                                        "VacStep3LV");       //its name


  fVacStepPV3 = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0,0,(AirgabCrystalPMTthick+PMTGlassthick)/2),    //its position //old 0.,0.,aluwraplength/2
                                fVacStepLV3,            //its logical volume
                                "VacStep3PV",                 //its name
                                fPMTmountLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  fVacStepLV3->SetVisAttributes(VacStep3Vis);

  //
  // vacuum step 4
  // all in front of the crystall

  auto fVacStepS4 = new G4Box("VacStep4S",  //Name
                               calorcellxy/2.,
                               calorcellxy/2,
                               vacthick/2.);

  auto fVacStepLV4 = new G4LogicalVolume(fVacStepS4,    //its solid
                                        fWorldMaterial,    //its material
                                        "VacStep4LV");       //its name

  fVacStepPV4 = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,-venylfoillength/2),    //its position //old 0.,0.,aluwraplength/2
                                fVacStepLV4,            //its logical volume
                                "VacStep4PV",                 //its name
                                fCaloCellLV,               //its mother //old fCaloCellLV
                                false,                     //no boolean operat
                                0);                        //copy number

   fVacStepLV4->SetVisAttributes(VacStep3Vis);
  // fVacStepLV4->SetVisAttributes(G4VisAttributes::GetInvisible());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Scintillator Fingers  (SF)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(SFStatus=="true"){

  //Vinylfoil for the Scintillator fingers
  //
  auto SFVinylS= new G4Box("SFVinylfoilS",  //Name
                                SFVinylx/2.,   // x size
                                SFVinyly/2.,     // y size
                                SFVinyllength/2.); // z size


  auto SFVinylLV = new G4LogicalVolume(SFVinylS,    //its solid
                                         Venyl,    //its material
                                         "SFVinylfoilLV");       //its name

  auto SFVinylPV = new G4PVPlacement(SFRotation,                   //no rotation
                         G4ThreeVector(SFpositonx,SFpositony,SFpositonZ),    //its position // old 0.,0.,-vacthick/2
                                 SFVinylLV,            //its logical volume
                                 "SFVinylfoilPV",                 //its name
                                LogicalWorld,               //its mother
                                 false,                     //no boolean operat
                                     0);                        //copy number

  SFVinylLV->SetVisAttributes(VenyfoilVis);

  // Alu-wrapping for SF
  //
  auto SFAluWrapS= new G4Box("SFAluWrappingS",  //Name
                                SFAlux/2.,   // x size
                                SFAluy/2.,     // y size
                                SFAlulength/2.); // z size


  auto SFAluwrapLV = new G4LogicalVolume(SFAluWrapS,    //its solid
                                         Al,    //its material
                                         "SFAluWrappingLV");       //its name

  auto SFAluwrapPV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,-(SFVacstepthick-SFVinylthick)/2),    // its positon
                                 SFAluwrapLV,            //its logical volume
                                 "SFAluWrappingPV",                 //its name
                                 SFVinylLV,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

  SFAluwrapLV->SetVisAttributes(AluVis);

  //
  // AirGap for SF Material and Aluwrapping
  //
  auto SFAirGapS= new G4Box("SFAirGapS",  //Name
                               SFairgapx/2.,   // x size
                               SFairgapy/2.,     // y size
                               SFairgaplength/2.); // z size


  auto SFAirGapLV = new G4LogicalVolume(SFAirGapS,    //its solid
                                        Air_NoRI,    //its material
                                        "SFAirGapLV");       //its name

  auto SFAirGapPV = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,SFAluthick/2),    //its position
                                SFAirGapLV,            //its logical volume
                                "SFAirGapPV",                 //its name
                                SFAluwrapLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  //SFAirGapLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  SFAirGapLV->SetVisAttributes(AirVis);

  //
  // Scintillator
  //
  auto SFS= new G4Box("ScintillarofingerS",  //Name
                               SFx/2.,   // x size
                               SFy/2.,     // y size
                               SFlength/2.); // z size


  auto SFLV = new G4LogicalVolume(SFS,    //its solid
                                        SFMat,    //its material
                                        "ScintillarofingerLV");       //its name

  auto SFPV = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0.,SFairgapthick/2),    //its position
                                SFLV,            //its logical volume
                                "ScintillarofingerPV",                 //its name
                                SFAirGapLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  SFLV->SetVisAttributes(CrystalVis);


  //vacuum step 5 this is the detector for the Scintillaror fingers
  auto fVacStepS5 = new G4Box("VacStep5Solit",  //Name
                               SFVacstepx/2,//DetVacstepxy/2.,
                               SFVacstepy/2,// DetVacstepxy/2,
                               SFVacstepthick/2.);

  auto fVacStepLV5 = new G4LogicalVolume(fVacStepS5,    //its solid
                                        Vacuum,    //its material
                                        "VacStep5LV");       //its name


  auto fVacStepPV5 = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0,0,(SFAlulength+SFVinylthick)/2),    //its position //old 0.,0.,aluwraplength/2
                                fVacStepLV5,            //its logical volume
                                "VacStep5PV",                 //its name
                                SFVinylLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  fVacStepLV5->SetVisAttributes(VacStep3Vis);
  G4cout<<"Scintillator-Finger is constructed"<<G4endl;
 }
 else {  G4cout<<"There is no Scintillator-Finger in the beam"<<G4endl;}

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Lanex scintillator screen
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 if(LanexStatus == "true"){

  auto scintArmSolid = new G4Box("scintArmSolid", scx/2., scy/2., (scz1+scz2+scz3)/2.);
  auto scintArmLogical = new G4LogicalVolume(scintArmSolid,fWorldMaterial ,"scintArmLogical");

  auto scintBaseSolid = new G4Box("scintBaseBox", scx/2., scy/2., scz1/2.);
  auto scintBaseLogical = new G4LogicalVolume(scintBaseSolid, Pstyrene,"scintBaseLogical");

  auto scintPhosphorSolid = new G4Box("scintPhosphorBox", scx/2., scy/2., scz2/2.);
  auto scintPhosphorLogical = new G4LogicalVolume(scintPhosphorSolid, lanex,"scintPhosphorLogical");

  auto scintFinishSolid = new G4Box("scintFinishBox", scx/2., scy/2., scz3/2.);
  auto scintFinishLogical = new G4LogicalVolume(scintFinishSolid, Pstyrene,"scintFinishLogical");


  new G4PVPlacement(0,G4ThreeVector(0.,0.,-(scz2+scz3)/2.),
                    scintBaseLogical,"scintBasePhysical",scintArmLogical,false,0);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,(scz1 - scz3)/2.),
                    scintPhosphorLogical,"scintPhosphorPhysical",scintArmLogical,false,0);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,(scz1+scz2)/2.),
                    scintFinishLogical,"scintFinishPhysical",scintArmLogical,false,0);

  new G4PVPlacement(0,G4ThreeVector(ToyX,0,ZposLanex),
                    scintArmLogical,"scintArmPhysical",LogicalWorld,false,0);

  scintPhosphorLogical->SetVisAttributes(LanexVis);
  scintFinishLogical->SetVisAttributes(PolyVis1);
  scintBaseLogical->SetVisAttributes(PolyVis2);
  scintArmLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4cout<<"LANEX screen is constructed"<<G4endl;
 }
 else {  G4cout<<"There is no Lanex screen in the beam"<<G4endl;}
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Cherencov Straw for LUXE
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 if(StrawStatus == "true"){

  auto StrawS= new G4Tubs("Straw", 0.,rStraw,LStraw/2., 0.0*deg, 360.0*deg );
  auto StrawLV = new G4LogicalVolume(StrawS,StrawMat,"StawLogical");

  auto MylarS= new G4Tubs("Mylar", 0.,rMylar,LStraw/2., 0.0*deg, 360.0*deg );
  auto MylarLV = new G4LogicalVolume(MylarS,Mylar,"MylarLogical");

  auto CuAuS= new G4Tubs("CuAu", 0.,rCuAu,LStraw/2., 0.0*deg, 360.0*deg );
  auto CuAuLV = new G4LogicalVolume(CuAuS,CuAu,"CuAuLogical");

  auto OilS= new G4Tubs("Oil", 0.,rOil,LStraw/2., 0.0*deg, 360.0*deg );
  auto OilLV = new G4LogicalVolume(OilS,Air,"OilLogical");


  new G4PVPlacement(StrawRotation, G4ThreeVector(ToyX,0,ZposStraw),
                    StrawLV,"StrawPhysical", LogicalWorld, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0,0,0),
                    MylarLV,"MylarPhysical", StrawLV, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0,0,0),
                    CuAuLV,"CuAuPhysical", MylarLV, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0,0,0),
                    OilLV,"OilPhysical", CuAuLV, false, 0);

  StrawLV->SetVisAttributes(StrawVis);
  G4cout<<"cherenkov straw is constructed"<<G4endl;
 }
 else {  G4cout<<"There is no cherenkov straw in the beam"<<G4endl;}

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Beamline geometry
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 if(beamlineStatus=="on"){

   G4double concreteBlocklength = 80. *cm;  //length of the conrete block
   G4double concreteBlockhigth= 120.*cm;   //hight of the conrete block
   G4double concreteBlockwidth= 50.*cm;    //width of the conrete block

   G4double  concreteHolehigth =  20.*cm;     //hight of the hole in the conrete concrete block
   G4double  concreteHolewidth =  9.5*cm;     //width of the hole in the conrete concrete block

   G4double  collimatorwidth   =  10.*cm;
   G4double  collimatorlength  =  20.*cm;     //lenght of the colimator copper block
   G4double  collimatorradius  =  RCollimator;      //radius of the colimator hole


   G4double distColltoCBlock= 2.*cm;         //distance from the Collimator to concrete block
   G4double distoCBlock= 18.*cm;            //distance from the aluwindow of the beam line ot to the conrete block

   G4double Aluwindowthick= 0.2*mm;         // thickness of the aluwindow
   G4double beamlineVaclength=3 *m;      // lenght of vaccum in the beam line
   G4double beamlinelength= Aluwindowthick+ beamlineVaclength;// lenght of the "toy" beamline
   G4double rInner= 1.85 *cm;                 // inner radius of the beam line
   G4double rOuter= 1.925 *cm;               //outer radius of the beam line




   //Concrete Block
   auto fconcreteBlockS= new G4Box("concreteBlock",         // Name
                                concreteBlockwidth/2.,     // x size
                                concreteBlockhigth/2.,     // y size
                                concreteBlocklength/2.);    // z size


   auto fconcreteBlockLV = new G4LogicalVolume(fconcreteBlockS,  //its solid
                                         Concrete,    //its material
                                         "concreteBlock");  //its name

   auto fconcreteBlockPV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,-415,-(collimatorlength+distColltoCBlock+concreteBlocklength/2.)),    //its position
                                 fconcreteBlockLV,            //its logical volume
                                 "concreteBlock",                 //its name
                                LogicalWorld,               //its mother
                                 false,                     //no boolean operat
                                 0);                 //copy number

   fconcreteBlockLV->SetVisAttributes(ConcreteVis);



   // hole in conrete block
   auto fconcreteHoleS= new G4Box("concreteHole",         // Name
                                concreteHolewidth/2.,     // x size
                                concreteHolehigth/2.,     // y size
                                concreteBlocklength/2.);    // z size


   auto fconcreteHoleLV = new G4LogicalVolume(fconcreteHoleS,  //its solid
                                         fWorldMaterial,    //its material
                                         "concreteHole");  //its name

   auto fconcreteHolePV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,415.,0.),    //its position
                                 fconcreteHoleLV,            //its logical volume
                                 "concreteHole",                 //its name
                                 fconcreteBlockLV,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

  fconcreteHoleLV->SetVisAttributes(ConcreteHoleVis);


   if (RCollimator>0){
   // copper collimator
   auto fcollimatorS= new G4Box("collimator",         // Name
                                collimatorwidth/2.,     // x size
                                collimatorwidth/2.,     // y size
                                collimatorlength/2.);    // z size


   auto fcollimatorLV = new G4LogicalVolume(fcollimatorS,  //its solid
                                          Copper,    //its material
                                         "collimator");  //its name

   auto fcollimatorPV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,-collimatorlength/2),    //its position
                                 fcollimatorLV,            //its logical volume
                                 "collimator",                 //its name
                                 LogicalWorld,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

   fcollimatorLV->SetVisAttributes(CollimatorVis);


   // copper collimator hole
   auto fcollimatorHoleS= new G4Tubs("collimatorhole", //name
                                  0., // inner radius
                                  collimatorradius,  // outer radius
                                  collimatorlength/2., // half length in z
                                  0.0*deg,  // starting angle
                                  360.0*deg ); // tota angle

   auto fcollimatorHoleLV = new G4LogicalVolume(fcollimatorHoleS,  //its solid
                                          fWorldMaterial,    //its material
                                         "collimatorhole");  //its name

   auto fcollimatorHolePV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,0.),    //its position
                                 fcollimatorHoleLV,            //its logical volume
                                 "collimatorhole",                 //its name
                                 fcollimatorLV,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

   }


   // Lead shelding 2 Blocks next to the Collimator
   auto fBlei1S= new G4Box("Blei1",         // Name
                              200/2.,     // x size
                              250/2.,     // y size
                              100/2.);    // z size


   auto fBlei1LV = new G4LogicalVolume(fBlei1S,  //its solid
                                         shieldMat,    //its material
                                         "Blei1");  //its name

   auto fBlei1PV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(-160,40,-165),    //its position
                                 fBlei1LV,            //its logical volume
                                 "Blei1",                 //its name
                                LogicalWorld,               //its mother
                                 false,                     //no boolean operat
                                 0);                 //copy number

   auto fBlei2S= new G4Box("Blei2",         // Name
                              200/2.,     // x size
                              300/2.,     // y size
                              100/2.);    // z size


   auto fBlei2LV = new G4LogicalVolume(fBlei2S,  //its solid
                                         shieldMat,    //its material
                                         "Blei2");  //its name

   auto fBlei2PV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(160,65,-165),    //its position
                                 fBlei2LV,            //its logical volume
                                 "Blei2",                 //its name
                                LogicalWorld,               //its mother
                                 false,                     //no boolean operat
                                 0);                 //copy number

   fBlei2LV->SetVisAttributes(BleiVis);
   fBlei1LV->SetVisAttributes(BleiVis);

   //
   // beam line tube
   auto fbeamlineTubeS= new G4Tubs("beamline", //name
                                  0., // inner radius
                                  rOuter,  // outer radius
                                  beamlinelength/2., // half length in z
                                  0.0*deg,  // starting angle
                                  360.0*deg ); // tota angle

   auto fbeamlineTubeLV = new G4LogicalVolume(fbeamlineTubeS,  //its solid
                                          StainlessSteel,    //its material
                                         "beamline");  //its name

   auto fbeamlineTubePV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,-(collimatorlength+distColltoCBlock+concreteBlocklength+distoCBlock+beamlinelength/2.)),    //its position
                                 fbeamlineTubeLV,            //its logical volume
                                 "beamline",                 //its name
                                 LogicalWorld,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

   fbeamlineTubeLV->SetVisAttributes(BeamTubeVis);

   // beam line Vaccuum
   auto fbeamlineVacS= new G4Tubs("beamlineVac", //name
                                  0., // inner radius
                                  rInner,  // outer radius
                                  beamlineVaclength/2., // half length in z
                                  0.0*deg,  // starting angle
                                  360.0*deg ); // tota angle

   auto fbeamlineVacLV = new G4LogicalVolume(fbeamlineVacS,  //its solid
                                          Vacuum,    //its material
                                         "beamlineVac");  //its name

   auto fbeamlineVacPV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,-Aluwindowthick/2.),    //its position
                                 fbeamlineVacLV,            //its logical volume
                                 "beamlineVac",                 //its name
                                 fbeamlineTubeLV,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

   fbeamlineVacLV->SetVisAttributes(BeamLineVacVis);


   // Aluminum Window
   auto fAluWindowS= new G4Tubs("AluWindow", //name
                                  0., // inner radius
                                  rOuter,  // outer radius
                                  Aluwindowthick/2., // half length in z
                                  0.0*deg,  // starting angle
                                  360.0*deg ); // tota angle

   auto fAluWindowLV = new G4LogicalVolume(fAluWindowS,  //its solid
                                          Al,    //its material
                                         "AluWindow");  //its name

   auto fAluWindowPV = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(0.,0.,beamlineVaclength/2.),    //its position
                                 fAluWindowLV,            //its logical volume
                                 "AluWindow",                 //its name
                                 fbeamlineTubeLV,               //its mother
                                 false,                     //no boolean operat
                                 0);                        //copy number

   fAluWindowLV->SetVisAttributes(AluWindowVis);
 }



 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Optical boundary surfaces
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // Al <-> Air optical surface copied from QuaSim
 //
 G4MaterialPropertiesTable *AlWrapProperty = new G4MaterialPropertiesTable();
 const G4int n_AlAir=8;
 // G4double PE_AlWrap[n_AlAir] = {1.38*eV, 1.90*eV, 2.38*eV, 2.48*eV, 5.17*eV, 5.64*eV, 6.20*eV, 7.77*eV };
 // G4double RE_AlWrap[n_AlAir] = {0.82,    0.83,    0.84,    0.85,    0.86,    0.84,    0.81,    0.74 };
 G4double PE_AlWrap[n_AlAir] = {1.38*eV, 1.90*eV, 2.38*eV, 2.48*eV, 5.17*eV, 5.64*eV, 6.20*eV, 6.7*eV };
 G4double RE_AlWrap[n_AlAir] = {0.82,    0.83,    0.84,    0.85,    0.86,    0.84,    0.81,    0.79 };
 AlWrapProperty -> AddProperty("REFLECTIVITY", PE_AlWrap, RE_AlWrap, n_AlAir);

 // optical and logical surface
 G4OpticalSurface* OpAlWrapSurface = new G4OpticalSurface("AirAluSurface");
 OpAlWrapSurface -> SetType(dielectric_metal);
 OpAlWrapSurface -> SetFinish(ground);
 G4double AlPolish = 0; // 1 = smooth, 0 = maximum roughness
 OpAlWrapSurface -> SetPolish(AlPolish);
 OpAlWrapSurface -> SetModel(glisur);
 OpAlWrapSurface -> SetMaterialPropertiesTable(AlWrapProperty);
 // need to attach them to all the physical volumes for the alu
 G4LogicalBorderSurface* AlWrapSurface =
          new G4LogicalBorderSurface("AirAluSurface",fAluwrapPV , fAlAirGapPV, OpAlWrapSurface);
 G4LogicalBorderSurface* AlWrapSurface2 =
          new G4LogicalBorderSurface("AirAluSurface", fAlAirGapPV, fAluwrapPV, OpAlWrapSurface);


  //

  // here we take the defould model dielectric_dielectric with a perfect smooth surface so we do not have to define a border surface!
  /// - Quartz <-> Air: Unified model (taken from QuaSi) /have to check if TF1 <-> Air is different
  // properties table
  // G4MaterialPropertiesTable* AirQSurfProp = new G4MaterialPropertiesTable();
  // const G4int n_AirQ=2;
  // G4double OpAirSpecularlobe[n_AirQ] = {1.0, 1.0};
  // G4double OpAirSpecularspike[n_AirQ] = {0.0, 0.0};
  // G4double OpAirBackscatter[n_AirQ] = {0.0, 0.0};
  //
  // G4double PE_OpAir[n_AirQ] = {1.38*eV, 6.70*eV};
  // G4double RI_OpAir[n_AirQ] = {1.00029, 1.00029};
  // AirQSurfProp -> AddProperty("RINDEX", PE_OpAir, RI_OpAir, n_AirQ);
  // AirQSurfProp -> AddProperty("SPECULARLOBECONSTANT", PE_OpAir, OpAirSpecularlobe, n_AirQ);
  // AirQSurfProp -> AddProperty("SPECULARSPIKECONSTANT", PE_OpAir, OpAirSpecularspike, n_AirQ);
  // AirQSurfProp -> AddProperty("BACKSCATTERCONSTANT", PE_OpAir, OpAirBackscatter, n_AirQ);
  //
  // // optical and logical surface
  // G4OpticalSurface* OpAirSurface = new G4OpticalSurface("QAirSurface");
  // OpAirSurface -> SetType(dielectric_dielectric);
  // OpAirSurface -> SetModel(unified);
  //
  // // from Janacek, Morses, Simulating Scintillator Light ..., IEEE 2010:
  // // polished 1.3 degree, etched 3.8 degree, ground 12 degree
  // //          0.023              0.066              0.21
  // // degree   1      2      3      4      5      6      7      8      9      10     11     12
  // // radians  0.017  0.035  0.052  0.070  0.087  0.105  0.122  0.140  0.157  0.175  0.192  0.209
  // // degree   13     14     15
  // // radians  0.227  0.244  0.262
  //
  // OpAirSurface -> SetFinish(ground);
  // //  G4double sigma_alpha = 0.023;
  // //  G4double sigma_alpha = 0.066;
  // G4double sigma_alpha = 0.209;
  // OpAirSurface -> SetSigmaAlpha(sigma_alpha);
  // OpAirSurface -> SetMaterialPropertiesTable(AirQSurfProp);
  // G4LogicalBorderSurface* AirSurface =
  // 	new G4LogicalBorderSurface("QAirSurface", fDetectorPV, fAlAirGapPV, OpAirSurface);
  // G4LogicalBorderSurface* AirSurface2 =
  // 	new G4LogicalBorderSurface("QAirSurface", fAlAirGapPV, fDetectorPV, OpAirSurface);

  } //end if-statement for the calorimeter geometry

  PrintParameters();

  //always return the root volume
  //
  return PhysicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{ if(versionType=="Pol" || versionType=="PolCal"){
  G4cout << "\n The ConverterTarget is made of " << fConvMaterial->GetName()
          << " , " << G4BestUnit(fSizeXY,"Length")<<  "in diameter and "
         <<  G4BestUnit(fConvThick,"Length") << " thick"
           << G4endl;
 G4cout << "\n The IronCore is"
        <<  G4BestUnit(fCoreThick,"Length") << " thick"
          << G4endl;
  }

  if(versionType=="Cal" || versionType=="PolCal"){
    G4cout << "\n The Calorimeter is made of " << fCaloMaterial->GetName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetConvMaterial(G4String materialChoice)
{
  // search the material by its name
  // G4Material* mat =
  //   G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  G4Material* mat = allMaterials->GetMat(materialChoice);
  if (mat != fConvMaterial) {
    if(mat) {
      fConvMaterial = mat;
      UpdateGeometry();
    } else {
      G4cout << "### Warning!  Converter Target material: <"
           << materialChoice << "> not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  // G4Material* mat =
  //   G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  G4Material* mat = allMaterials->GetMat(materialChoice);
  if (mat != fWorldMaterial) {
    if(mat) {
      fWorldMaterial = mat;
      UpdateGeometry();
    } else {
      G4cout << "### Warning! World material: <"
           << materialChoice << "> not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCaloMaterial(G4String materialChoice)
{
  // search the material by its name
  // G4Material* mat =
  //   G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  G4Material* mat = allMaterials->GetMat(materialChoice);
  if (mat != fCaloMaterial) {
    if(mat) {
      fCaloMaterial = mat;
      UpdateGeometry();
    } else {
      G4cout << "### Warning! Calorimeter material: <"
           << materialChoice << "> not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeXY(G4double value)
{
  fSizeXY = value;
  if (fWorldSize<fSizeXY) fWorldSize = 10*fSizeXY;
  UpdateGeometry();
}

void DetectorConstruction::SetCoreThick(G4double value)
{
  fCoreThick = value;
  if (fWorldSize<fCoreThick) fWorldSize = 10*fCoreThick;
  UpdateGeometry();
}

void DetectorConstruction::SetConvThick(G4double value)
{
  fConvThick = value;
  UpdateGeometry();
}

void DetectorConstruction::SetCrystalnumber(G4String value)
{
  CrystalNumber=value;
  UpdateGeometry();
}

void DetectorConstruction::SetCollimatorRadius(G4double value)
{
  RCollimator = value;
  UpdateGeometry();
}

void DetectorConstruction::SetCaloDistance(G4double value)
{
  dCalo = value;
  UpdateGeometry();
}

void DetectorConstruction::SetSFStatus(G4String value)
{
  SFStatus = value;
  UpdateGeometry();
}

void DetectorConstruction::SetLanexStatus(G4String value)
{
  LanexStatus = value;
  UpdateGeometry();
}

void DetectorConstruction::SetStrawStatus(G4String value)
{
  StrawStatus = value;
  UpdateGeometry();
}

void DetectorConstruction::SetCaloXposition(G4double value)
{
  CaloXpos = value;
  UpdateGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  if (PhysicalWorld){
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
