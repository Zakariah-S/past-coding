#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
    //Set up particle gun to shoot 1 at a time
    /*fParticleGun = new G4ParticleGun(1);*/
    fParticleGun = new G4GeneralParticleSource();

    /*Set up particle gun to shoot protons with given energy, 
    momentum, momentum direction, initial position*/
    /*G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();    
    G4ParticleDefinition *particle = particleTable->FindParticle("geantino");*/

    /*ParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));    
    fParticleGun->SetParticleMomentum(0.*GeV);
    fParticleGun->SetParticleDefinition(particle);*/
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
    /*G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();

    if (particle == G4Geantino::Geantino())
    {
        G4int Z = 90;
        G4int A = 232;
        G4double charge = -1.*eplus;
        G4double energy = 0.*keV;

        G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(charge);
        fParticleGun->SetParticleMomentum(0.*GeV);

    }*/

    fParticleGun->GeneratePrimaryVertex(anEvent);
}