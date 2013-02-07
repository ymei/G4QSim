#ifndef __RUNACTION_H__
#define __RUNACTION_H__

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <map>

class G4Run;
class PrimaryGeneratorAction;
class AnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(PrimaryGeneratorAction*, AnalysisManager*);
   ~RunAction();
   
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
    
    void ParticleCount(G4String, G4double);
    void Balance(G4double,G4double);
    void EventTiming(G4double);
    
  private:
    PrimaryGeneratorAction* m_primaryGeneratorAction;
    AnalysisManager *m_analysisManager;
    
    std::map<G4String,G4int> m_particleCount;
    std::map<G4String,G4double> m_eMean;
    std::map<G4String,G4double> m_eMin;
    std::map<G4String,G4double> m_eMax;
    G4int    m_decayCount;
    G4double m_eBalance[3];
    G4double m_pBalance[3];
    G4double m_eventTime[3];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __RUNACTION_H__ */
