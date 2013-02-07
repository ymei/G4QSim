#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"

#include "TrackingAction.hh"
#include "TrackingMessenger.hh"
#include "RunAction.hh"
#include "EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(RunAction* RA, EventAction* EA)
    :m_runAction(RA),m_eventAction(EA),m_fullChainQ(true)
{
    m_trackingMessenger = new TrackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
    delete m_trackingMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
#if 0
    G4ParticleDefinition* particle = track->GetDefinition();
    G4String name   = particle->GetParticleName();
    m_charge = particle->GetPDGCharge();
    m_mass   = particle->GetPDGMass();  
    
    G4double Ekin = track->GetKineticEnergy();
    G4int ID      = track->GetTrackID();
  
    G4bool condition = false;  

    //count particles
    //
    m_runAction->ParticleCount(name, Ekin);
  
    /*
    if(particle == G4Electron::Electron())
    G4cout << "ElectronEk: " << Ekin << " TrackID: " << ID << G4endl;
    */
  
    //fullChain: stop ion and print decay chain
    //
    if (m_charge > 2.) {
        G4Track* tr = (G4Track*) track;
        if (m_fullChainQ) tr->SetTrackStatus(fStopButAlive);
        if (ID == 1) m_eventAction->AddDecayChain(name);
        else         m_eventAction->AddDecayChain(" ---> " + name); 
    }
  
    //example of saving random number seed of this event, under condition
    //
    // if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
#if 0
    //keep only ions
    //
    if (m_charge < 3. ) return;

    //get time
    //   
    G4double time = track->GetGlobalTime();

    //energy and momentum balance (from secondaries)
    //
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();

    size_t nbtrk = (*secondaries).size();
    if (nbtrk) {
        //there are secondaries --> it is a decay
        //
        //force 'single' decay
        G4int ID = track->GetTrackID();
        if ((!m_fullChainQ) && (ID > 1)) G4RunManager::GetRunManager()->AbortEvent();
        //
        //balance
        G4double Ebalance = - track->GetTotalEnergy();
        G4ThreeVector Pbalance = - track->GetMomentum();
        for (size_t itr=0; itr<nbtrk; itr++) {
            G4Track* trk = (*secondaries)[itr];
            Ebalance += trk->GetTotalEnergy();
            //exclude gamma desexcitation from momentum balance
            if (trk->GetDefinition() != G4Gamma::Gamma())     
                Pbalance += trk->GetMomentum();             
        }
        G4double Pbal = Pbalance.mag();  
        m_runAction->Balance(Ebalance,Pbal);  
    }

    //no secondaries --> end of chain    
    //  
    if (!nbtrk) {
        m_runAction->EventTiming(time); //time of life
    }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
