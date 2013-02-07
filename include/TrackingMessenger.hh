#ifndef __TRACKINGMESSENGER_H__
#define __TRACKINGMESSENGER_H__

#include "globals.hh"
#include "G4UImessenger.hh"

class TrackingAction;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingMessenger: public G4UImessenger
{
  public:
    TrackingMessenger(TrackingAction*);
   ~TrackingMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    TrackingAction* m_trackingAction;
    G4UIcmdWithABool* m_trackingCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __TRACKINGMESSENGER_H__ */
