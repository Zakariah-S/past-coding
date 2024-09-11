#include "construction.hh"
#include "detector.hh"

MyDetectorConstruction::MyDetectorConstruction() {}

MyDetectorConstruction::~MyDetectorConstruction() {}

G4VPhysicalVolume *MyDetectorConstruction::Construct() {

    //Get particle data
    G4NistManager *nist = G4NistManager::Instance();

    //Define SiO2
    /*G4Material *SiO2 = new G4Material("SiO2", 2.201*g/cm3, 2);
    SiO2->AddElement(nist->FindOrBuildElement("Si"), 1);
    SiO2->AddElement(nist->FindOrBuildElement("O"), 2);

    //Define water
    G4Material *H2O = new G4Material("H2O", 1*g/cm3, 2);
    H2O->AddElement(nist->FindOrBuildElement("O"), 1);
    H2O->AddElement(nist->FindOrBuildElement("H"), 2);

    //Get carbon as a usable element
    G4Element *C = nist->FindOrBuildElement("C");

    //Define aerogel (mix of SiO2, water, carbon)
    G4Material *Aerogel = new G4Material("Aerogel", 0.200*g/cm3, 3);
    Aerogel->AddMaterial(SiO2, 62.5*perCent);
    Aerogel->AddMaterial(H2O, 37.4*perCent);
    Aerogel->AddElement(C, 0.1*perCent);

    //Define refractive indices at particular photon energies
    G4double energy[2] = {1.239841939*eV/0.9, 1.239841939*eV/0.2};
    G4double rindexAerogel[2] = {1.1, 1.1};
    G4double rindexWorld[2] = {1., 1.};

    //Implement index n of aerogel
    G4MaterialPropertiesTable *mptAerogel = new G4MaterialPropertiesTable();
    mptAerogel->AddProperty("RINDEX", energy, rindexAerogel, 2);
    Aerogel->SetMaterialPropertiesTable(mptAerogel);*/

    //Get air as a usable material
    G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");

    //Give air appropriate refractive index
    /*G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
    mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
    worldMat->SetMaterialPropertiesTable(mptWorld);*/

    //Define world volume, fill it with air
    G4Box *solidWorld = new G4Box("solidWorld", 0.5*m, 0.5*m, 0.5*m);

    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");

    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

    //Define aerogel box
    /*G4Box *solidRadiator = new G4Box("solidRadiator", 0.4*m, 0.4*m, 0.01*m);

    G4LogicalVolume *logicRadiator = new G4LogicalVolume(solidRadiator, Aerogel, "logicRadiator");

    G4VPhysicalVolume *physRadiator = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.25*m), logicRadiator, "physRadiator", logicWorld, false, 0, true);*/

    /*G4Box *solidDetector = new G4Box("solidDetector", 0.005*m, 0.005*m, 0.01*m);

    logicDetector = new G4LogicalVolume(solidDetector, worldMat, "logicDetector");

    for(G4int i = 0; i < 100; i++)
    {
        for(G4int j = 0; j < 100; j++)
        {
            G4VPhysicalVolume *physDetector = new G4PVPlacement(0, G4ThreeVector(-0.5*m+(i+0.5)*m/100, -0.5*m+(j+0.5)*m/100, 0.49*m), logicDetector, "physDetector", logicWorld, false, j+i*100, true);
        }
    }*/

    G4double aluminiumThickness = 3.0*mm;
	G4double aluminiumBoxX = 2.5*cm;
	G4double aluminiumBoxY = 2.5*cm;
	G4double aluminiumBoxZ = 5.0*cm;

	G4double aluminiumBoxPosX = 0.*cm;
	G4double aluminiumBoxPosY = 0.*cm;
	G4double aluminiumBoxPosZ = 0.*cm;

	G4Box *solidOuterAluminiumBox = new G4Box("solidOuterAluminiumBox", aluminiumBoxX/2., aluminiumBoxY/2., aluminiumBoxZ/2.);
	G4Box *solidInnerAluminiumBox = new G4Box("solidInnerAluminiumBox", (aluminiumBoxX - 2*aluminiumThickness)/2, (aluminiumBoxY - 2*aluminiumThickness)/2, (aluminiumBoxZ - 2*aluminiumThickness)/2);
	G4SubtractionSolid *solidAluminiumBox = new G4SubtractionSolid("solidAluminiumBox", solidOuterAluminiumBox, solidInnerAluminiumBox, 0, G4ThreeVector(0.*mm, 0.*mm, 0.*mm));

	G4LogicalVolume *logicAluminiumBox = new G4LogicalVolume(solidAluminiumBox, nist->FindOrBuildMaterial("G4_Cu"), "logicAluminiumBox");

	G4VPhysicalVolume *physAluminiumBox = new G4PVPlacement(
			0,									// no rotation
			G4ThreeVector(aluminiumBoxPosX, aluminiumBoxPosY, aluminiumBoxPosZ),	// Placed at the center of the world volume
			logicAluminiumBox, 							// it's locical volume
			"AluminiumBox", 							// It's name
			logicWorld,								// logic volume of the mother
			false, 									// no boolean operator
			0, 									// copy number
			false); 								// checking for overlaps

    return physWorld;
}

//Set up the small bois as detectors (see detector.hh/cc for further code)
/*void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");
    logicDetector->SetSensitiveDetector(sensDet);
}*/
