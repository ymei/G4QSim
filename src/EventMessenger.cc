#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

#include "EventMessenger.hh"
#include "EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventMessenger::EventMessenger(EventAction* evAct)
    :m_eventAction(evAct)
{   
    m_dir = new G4UIdirectory("/G4QSim/");
    m_dir->SetGuidance("G4QSim top level");

    m_eventDir = new G4UIdirectory("/G4QSim/event/");
    m_eventDir ->SetGuidance("G4QSim event control");
      
    m_printCmd = new G4UIcmdWithAnInteger("/G4QSim/event/printModulo",this);
    m_printCmd->SetGuidance("Print events modulo n");
    m_printCmd->SetParameterName("eventNb",false);
    m_printCmd->SetRange("eventNb>0");
    m_printCmd->AvailableForStates(G4State_Idle);

    m_storeTrajectoryCmd = new G4UIcmdWithABool("/G4QSim/event/storeTrajectory",this);
    m_storeTrajectoryCmd->SetGuidance("Store particle trajectories into the data file");
    m_storeTrajectoryCmd->SetParameterName("storeTrajectory",true);
    m_storeTrajectoryCmd->SetDefaultValue(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventMessenger::~EventMessenger()
{
    delete m_storeTrajectoryCmd;
    delete m_printCmd;
    delete m_eventDir;
    delete m_dir;     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
    if (command == m_printCmd) {
        m_eventAction->SetPrintModulo(m_printCmd->GetNewIntValue(newValue));
    }
    if (command == m_storeTrajectoryCmd) {
        m_eventAction->SetStoreTrajectory(m_storeTrajectoryCmd->GetNewBoolValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
