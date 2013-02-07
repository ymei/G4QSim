#ifndef __EVENTMESSENGER_H__
#define __EVENTMESSENGER_H__

#include "G4UImessenger.hh"
#include "globals.hh"

class EventAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventMessenger: public G4UImessenger
{
public:
    EventMessenger(EventAction*);
    ~EventMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    EventAction* m_eventAction;
    
    G4UIdirectory*        m_dir;       
    G4UIdirectory*        m_eventDir;   
    G4UIcmdWithAnInteger* m_printCmd;
    G4UIcmdWithABool* m_storeTrajectoryCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __EVENTMESSENGER_H__ */
