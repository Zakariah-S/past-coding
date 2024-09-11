#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
//#include "G4SystemofUnits.hh"
//#include "G4ParticleTable.hh"
//#include "G4Geantino.hh"
//#include "G4IonTable.hh"
#include "G4GeneralParticleSource.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
    MyPrimaryGenerator();
    ~MyPrimaryGenerator();

    virtual void GeneratePrimaries(G4Event*);

private:
    G4GeneralParticleSource *fParticleGun = nullptr;
};

#endif