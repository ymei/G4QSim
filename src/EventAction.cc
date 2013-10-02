#include "EventAction.hh"
#include "EventMessenger.hh"

#include "G4Event.hh"
#include <iomanip>

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(AnalysisManager *analysisManager)
    :m_printModulo(1000),m_storeTrajectory(false),m_eventMessenger(NULL)
{
    m_eventMessenger = new EventMessenger(this);
    m_analysisManager = analysisManager;
    m_analysisManager->SetSaveInteval(m_printModulo*10);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
    delete m_eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
    m_decayChain = " ";

    if(event->GetEventID() % m_printModulo == 0) {
        G4cout << G4endl;
        G4cout << "------ Begin event " << event->GetEventID()
               << "------" << G4endl;
    }
        
    if(m_analysisManager)
        m_analysisManager->BeginOfEvent(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    if(m_analysisManager) {
        if(m_storeTrajectory)
            m_analysisManager->EndOfEventStoreTrajectory(event);
        else
            m_analysisManager->EndOfEvent(event);
    }

    //printing survey
    //
    G4int eventNb = event->GetEventID(); 
    if (eventNb % m_printModulo == 0) 
        G4cout << "------ End of event " << std::setw(6) << eventNb 
               << " :" + m_decayChain << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
