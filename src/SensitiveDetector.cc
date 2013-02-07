#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <map>
using std::map;

#include "SensitiveDetector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDetector::SensitiveDetector(const G4String& name)
    : G4VSensitiveDetector(name),
      m_hitsCollection(NULL)
{
    G4String cN = name + "_HC";
    collectionName.insert(cN);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDetector::~SensitiveDetector()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDetector::Initialize(G4HCofThisEvent *hcte)
{
    // Create hits collection
    m_hitsCollection = new SDHitsCollection(SensitiveDetectorName, collectionName[0]);

    // Add this collection in hcte
    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hcte->AddHitsCollection(hcID, m_hitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *history)
{
    SDHit *hit = new SDHit();

    G4Track *track = step->GetTrack();
    // energy deposit
    G4double energyDeposited = step->GetTotalEnergyDeposit();

    hit->SetTrackId(track->GetTrackID());
    if(!m_particleTypes.count(track->GetTrackID()))
        m_particleTypes[track->GetTrackID()] = track->GetDefinition()->GetParticleName();
    hit->SetParentId(track->GetParentID());
    hit->SetParticleType(track->GetDefinition()->GetParticleName());
    if(track->GetParentID())
        hit->SetParentType(m_particleTypes[track->GetParentID()]);
    else
        hit->SetParentType(G4String("none"));
    if(track->GetCreatorProcess())
        hit->SetCreatorProcess(track->GetCreatorProcess()->GetProcessName());
    else
        hit->SetCreatorProcess(G4String("Null"));
    hit->SetDepositingProcess(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
    hit->SetPosition(step->GetPostStepPoint()->GetPosition());
    hit->SetEnergyDeposited(energyDeposited);
    hit->SetKineticEnergy(track->GetKineticEnergy());
    hit->SetTime(track->GetGlobalTime());

    m_hitsCollection->insert(hit);

/*
    G4TouchableHistory* touchable
        = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    // Get calorimeter cell id
    G4int layerNumber = touchable->GetReplicaNumber(1);

    // Get hit accounting data for this cell
    SDHit* hit = (*m_hitsCollection)[layerNumber];
    if (!hit) {
        G4cerr << "Cannot access hit " << layerNumber << G4endl;
        exit(1);
    }
*/
/*
    if(track->GetDefinition()->GetParticleName() == "e-") {
        (*m_hitsCollection)[0]->Add(eDep, stepLength);
    } else if(track->GetDefinition()->GetParticleName() == "e+") {
        (*m_hitsCollection)[1]->Add(eDep, stepLength);
    } else if(track->GetDefinition()->GetParticleName() == "gamma") {
        (*m_hitsCollection)[2]->Add(eDep, stepLength);
    } else if(track->GetDefinition()->GetParticleName() == "opticalphoton") {
        (*m_hitsCollection)[3]->Add(eDep, stepLength);
    } else {
        return false;
    }
*/
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDetector::EndOfEvent(G4HCofThisEvent *hcte)
{
    if(verboseLevel>0) { // `/hits/verbose x' to change this verbose level
        G4int nbHits = m_hitsCollection->entries();
        G4cout << "\n-------->Hits Collection: in this event there are " << nbHits
               << " hits in the SD: " << G4endl;
        for (G4int i=0; i<nbHits; i++) (*m_hitsCollection)[i]->Print();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
