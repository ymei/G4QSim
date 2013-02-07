#ifndef __EVENTACTION_H__
#define __EVENTACTION_H__

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "AnalysisManager.hh"

class EventMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(AnalysisManager *analysisManager);
   ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    void SetPrintModulo(G4int val)   {m_printModulo = val;}
    void SetStoreTrajectory(G4bool val) {m_storeTrajectory = val;}
    void AddDecayChain(G4String val) {m_decayChain += val;}
               
  private:
    G4int            m_printModulo;
    G4bool           m_storeTrajectory;
    G4String         m_decayChain;                   
    EventMessenger*  m_eventMessenger;
    AnalysisManager* m_analysisManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __EVENTACTION_H__ */
