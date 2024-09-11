#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name) {}

MySensitiveDetector::~MySensitiveDetector() {}

//Do stuff when stuff hits the detectors
G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    //Get track info
    G4Track *track = aStep->GetTrack();

    //kill track (stop photon from travelling to the next detector over)
    track->SetTrackStatus(fStopAndKill);

    //Get info about start and end of track
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    //Get and print position of photon when it enters detector
    G4ThreeVector posPhoton = preStepPoint->GetPosition();
    //G4cout << "Photon position" << posPhoton << G4endl;

    //Get detector object and print its copy number (defined in detector.cc)
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int copyNo = touchable->GetCopyNumber();
    //G4 cout << "Copy Number: " << copyNo << G4endl;

    //Get and print position of the detector-copy that was struck
    G4VPhysicalVolume *physvol = touchable->GetVolume();
    G4ThreeVector posDetector = physvol->GetTranslation();

    G4cout << "Detector position: " << posDetector << G4endl;
}