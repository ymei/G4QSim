#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(PrimaryGeneratorAction* primaryGeneratorAction, AnalysisManager* analysisManager)
    :m_primaryGeneratorAction(primaryGeneratorAction)
{
    m_analysisManager = analysisManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{ 
    if(m_analysisManager)
        m_analysisManager->BeginOfRun(run);

    //initialize arrays
    //
    m_decayCount = 0;
    for (G4int i=0; i<3; i++) m_eBalance[i] = m_pBalance[i] = m_eventTime[i] = 0. ;
       
    //inform the runManager to save random number seed
    //
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleCount(G4String name, G4double Ekin)
{
    m_particleCount[name]++;
    m_eMean[name] += Ekin;
    if (m_particleCount[name] == 1) m_eMin[name] = m_eMax[name] = Ekin;
    if (Ekin < m_eMin[name]) m_eMin[name] = Ekin;
    if (Ekin > m_eMax[name]) m_eMax[name] = Ekin;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Balance(G4double Ebal, G4double Pbal)
{
    m_decayCount++;
    m_eBalance[0] += Ebal;
    if (m_decayCount == 1) m_eBalance[1] = m_eBalance[2] = Ebal;
    if (Ebal < m_eBalance[1]) m_eBalance[1] = Ebal;
    if (Ebal > m_eBalance[2]) m_eBalance[2] = Ebal;
  
    m_pBalance[0] += Pbal;
    if (m_decayCount == 1) m_pBalance[1] = m_pBalance[2] = Pbal;  
    if (Pbal < m_pBalance[1]) m_pBalance[1] = Pbal;
    if (Pbal > m_pBalance[2]) m_pBalance[2] = Pbal;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EventTiming(G4double time)
{  
    m_eventTime[0] += time;
    if (m_decayCount == 1) m_eventTime[1] = m_eventTime[2] = time;  
    if (time < m_eventTime[1]) m_eventTime[1] = time;
    if (time > m_eventTime[2]) m_eventTime[2] = time;    
}    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
    if(m_analysisManager)
        m_analysisManager->EndOfRun(run);

    G4int nbEvents = run->GetNumberOfEvent();
    if (nbEvents == 0) return;
#if 0 
    G4ParticleDefinition* particle = m_primaryGeneratorAction->GetParticleGun()
        ->GetParticleDefinition();
    G4String partName = particle->GetParticleName();
    G4double eprimary = m_primaryGeneratorAction->GetParticleGun()->GetParticleEnergy();
 
    G4cout << "\n ======================== run summary ======================";  
    G4cout << "\n The run was " << nbEvents << " " << partName << " of "
           << G4BestUnit(eprimary,"Energy");
    G4cout << "\n ===========================================================\n";
    G4cout << G4endl;

    G4int prec = 4, wid = prec + 2;
    G4int dfprec = G4cout.precision(prec);
      
    //particle count
    //
    G4cout << " Nb of generated particles: \n" << G4endl;
     
    std::map<G4String,G4int>::iterator it;               
    for (it = m_particleCount.begin(); it != m_particleCount.end(); it++) {
        G4String name = it->first;
        G4int count   = it->second;
        G4double eMean = m_eMean[name]/count;
        G4double eMin = m_eMin[name], eMax = m_eMax[name];    
         
        G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
               << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
               << "\t( "  << G4BestUnit(eMin, "Energy")
               << " --> " << G4BestUnit(eMax, "Energy") 
               << ")" << G4endl;       
    }
 
    //energy momentum balance
    //
    G4double Ebmean = m_eBalance[0]/m_decayCount;
    G4double Pbmean = m_pBalance[0]/m_decayCount;
  
    G4cout << "\n Energy and momentum balance : final state - initial state"
           << "\n (excluding gamma desexcitation from momentum balance) : \n"  
           << G4endl;
         
    G4cout 
        << "  Energy:   mean = " << std::setw(wid) << G4BestUnit(Ebmean, "Energy")
        << "\t( "  << G4BestUnit(m_eBalance[1], "Energy")
        << " --> " << G4BestUnit(m_eBalance[2], "Energy")
        << ")" << G4endl;
       
    G4cout 
        << "  Momentum: mean = " << std::setw(wid) << G4BestUnit(Pbmean, "Energy")
        << "\t( "  << G4BestUnit(m_pBalance[1], "Energy")
        << " --> " << G4BestUnit(m_pBalance[2], "Energy")
        << ")" << G4endl;
        
    //time of life
    //
    G4double Tmean = m_eventTime[0]/nbEvents;
    G4double halfLife = Tmean*std::log(2.);
   
    G4cout << "\n Time of life : mean = "
           << std::setw(wid) << G4BestUnit(Tmean, "Time")
           << "  half-life = "
           << std::setw(wid) << G4BestUnit(halfLife, "Time")
           << "   ( "  << G4BestUnit(m_eventTime[1], "Time")
           << " --> "  << G4BestUnit(m_eventTime[2], "Time")
           << ")" << G4endl;
        
    //activity
    //
    G4double molMass = particle->GetAtomicMass()*g/mole;
    G4double nAtoms = Avogadro/molMass;
    G4double ActivPerAtom = 1./Tmean;
    G4double ActivPerMass = ActivPerAtom*nAtoms;
   
    G4cout << "\n Activity = "
           << std::setw(wid) << ActivPerMass*g/becquerel
           << " Bq/g   ("    << ActivPerMass*g/curie
           << " Ci/g) \n" 
           << G4endl;

    // remove all contents in particleCount
    // 
    m_particleCount.clear(); 
    m_eMean.clear();  m_eMin.clear(); m_eMax.clear();

    // restore default precision
    // 
    G4cout.precision(dfprec);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
