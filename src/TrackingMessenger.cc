#include "TrackingMessenger.hh"
#include "TrackingAction.hh"

#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingMessenger::TrackingMessenger(TrackingAction* trackA)
    :m_trackingAction(trackA)
{
    m_trackingCmd = new G4UIcmdWithABool("/G4QSim/fullChain",this);
    m_trackingCmd->SetGuidance("allow full decay chain");
    m_trackingCmd->SetParameterName("fullChain",true);
    m_trackingCmd->SetDefaultValue(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingMessenger::~TrackingMessenger()
{
    delete m_trackingCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
    if (command == m_trackingCmd) {
        m_trackingAction->SetFullChain(m_trackingCmd->GetNewBoolValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
