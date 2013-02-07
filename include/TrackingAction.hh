#ifndef __TRACKINGACTION_H__
#define __TRACKINGACTION_H__

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class RunAction;
class EventAction;
class TrackingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction {

public:  
    TrackingAction(RunAction*, EventAction*);
    ~TrackingAction();
   
    void  PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*);
    
    void SetFullChain(G4bool flag) { m_fullChainQ = flag;};
    
private:
    RunAction*          m_runAction;
    EventAction*        m_eventAction;
    TrackingMessenger*  m_trackingMessenger;
    
    G4double m_charge, m_mass;        
    G4bool m_fullChainQ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __TRACKINGACTION_H__ */
