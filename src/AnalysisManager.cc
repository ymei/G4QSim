#include <G4SDManager.hh>
#include <G4Run.hh>
#include <G4Event.hh>
#include <G4HCofThisEvent.hh>
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4SystemOfUnits.hh"

#include <numeric>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>

#include "DetectorConstruction.hh"
#include "SDHit.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventData.hh"
#include "AnalysisManager.hh"


AnalysisManager::AnalysisManager(PrimaryGeneratorAction *primaryGeneratorAction)
{
    m_hitsCollectionID = -1;
    m_dataFilename = "events.root";

    m_primaryGeneratorAction = primaryGeneratorAction;
    m_eventData = new EventData();

    gROOT->ProcessLine("#include <vector>");
}

AnalysisManager::~AnalysisManager()
{
}

void
AnalysisManager::BeginOfRun(const G4Run *run)
{
    m_tfile = new TFile(m_dataFilename.c_str(), "RECREATE", "File containing event data from G4QSim");
    m_ttree = new TTree("t1", "Tree containing event data from G4QSim");

    m_ttree->Branch("eventid", &m_eventData->m_eventId, "eventid/I");
    m_ttree->Branch("nbhits", &m_eventData->m_nbHits, "nbhits/I");
    m_ttree->Branch("nbsteps", &m_eventData->m_nbSteps, "nbsteps/I");

    m_ttree->Branch("etot", &m_eventData->m_totalEnergyDeposited, "etot/D");
    m_ttree->Branch("trackid", "vector<int>", &m_eventData->m_trackId);
    m_ttree->Branch("parentid", "vector<int>", &m_eventData->m_parentId);
    m_ttree->Branch("type", "vector<string>", &m_eventData->m_particleType);
    m_ttree->Branch("parenttype", "vector<string>", &m_eventData->m_parentType);
    m_ttree->Branch("creatproc", "vector<string>", &m_eventData->m_creatorProcess);
    m_ttree->Branch("edproc", "vector<string>", &m_eventData->m_depositingProcess);
    m_ttree->Branch("xp", "vector<double>", &m_eventData->m_xp);
    m_ttree->Branch("yp", "vector<double>", &m_eventData->m_yp);
    m_ttree->Branch("zp", "vector<double>", &m_eventData->m_zp);
    m_ttree->Branch("ed", "vector<double>", &m_eventData->m_energyDeposited);
    m_ttree->Branch("ek", "vector<double>", &m_eventData->m_kineticEnergy);
    m_ttree->Branch("time", "vector<double>", &m_eventData->m_time);

    m_ttree->Branch("type_pri", "vector<string>", &m_eventData->m_primaryParticleType);
    m_ttree->Branch("etot_pri", &m_eventData->m_primaryParticleEnergy, "etot_pri/D");

    m_ttree->Branch("xm_pri", &m_eventData->m_primaryXm, "xm_pri/D");
    m_ttree->Branch("ym_pri", &m_eventData->m_primaryYm, "ym_pri/D");
    m_ttree->Branch("zm_pri", &m_eventData->m_primaryZm, "zm_pri/D");

    m_ttree->Branch("xp_pri", &m_eventData->m_primaryX, "xp_pri/D");
    m_ttree->Branch("yp_pri", &m_eventData->m_primaryY, "yp_pri/D");
    m_ttree->Branch("zp_pri", &m_eventData->m_primaryZ, "zp_pri/D");

    m_ttree->SetMaxTreeSize(1e12);
    // m_ttree->SetAutoSave(50000000);
    m_ttree->AutoSave();

    m_nbEventsToSimulateParameter = new TParameter<int>("nbevents", m_nbEventsToSimulate);
    m_nbEventsToSimulateParameter->Write();
}

void
AnalysisManager::EndOfRun(const G4Run *run)
{
    // in case of automatic file change
    m_tfile->CurrentFile()->Write();
    m_tfile->CurrentFile()->Close();
}

void
AnalysisManager::BeginOfEvent(const G4Event *event)
{
    if(m_hitsCollectionID == -1)
    {
        G4SDManager *sdManager = G4SDManager::GetSDMpointer();
        // ``Tracker_HC'' as defined in `DetectorConstruction.cc'
        m_hitsCollectionID = sdManager->GetCollectionID("Tracker_HC");
    } 
}

void
AnalysisManager::EndOfEvent(const G4Event *event)
{

    G4HCofThisEvent* hcte = event->GetHCofThisEvent();
    SDHitsCollection* hitsCollection = NULL;

    G4int nbHits = 0;

    if(hcte) {
        if(m_hitsCollectionID != -1) {
            hitsCollection = (SDHitsCollection *)(hcte->GetHC(m_hitsCollectionID));
            nbHits = (hitsCollection)?(hitsCollection->entries()):(0);
        }
    }

    m_eventData->m_primaryParticleType->push_back(m_primaryGeneratorAction->GetParticleTypeOfPrimary());
    m_eventData->m_primaryParticleEnergy = m_primaryGeneratorAction->GetEnergyOfPrimary()/keV;
    m_eventData->m_primaryXm = m_primaryGeneratorAction->GetMomentumDirectionOfPrimary().x();
    m_eventData->m_primaryYm = m_primaryGeneratorAction->GetMomentumDirectionOfPrimary().y();
    m_eventData->m_primaryZm = m_primaryGeneratorAction->GetMomentumDirectionOfPrimary().z();
    m_eventData->m_primaryX = m_primaryGeneratorAction->GetPositionOfPrimary().x()/mm;
    m_eventData->m_primaryY = m_primaryGeneratorAction->GetPositionOfPrimary().y()/mm;
    m_eventData->m_primaryZ = m_primaryGeneratorAction->GetPositionOfPrimary().z()/mm;

    if(nbHits) {
        m_eventData->m_nbHits = nbHits;
        m_eventData->m_eventId = event->GetEventID();

        // G4cout << "Event " << event->GetEventID() << ": There are " << nbHits << " hits" << G4endl;

        G4int nbSteps = 0;
        G4float totalEnergyDeposited = 0.;

        // SD hits
        for(G4int i=0; i<nbHits; i++) {
            SDHit *hit = (*hitsCollection)[i];

            if(true /* hit->GetParticleType() != "opticalphoton" */) {

                m_eventData->m_trackId->push_back(hit->GetTrackId());
                m_eventData->m_parentId->push_back(hit->GetParentId());

                m_eventData->m_particleType->push_back(hit->GetParticleType());
                m_eventData->m_parentType->push_back(hit->GetParentType());
                m_eventData->m_creatorProcess->push_back(hit->GetCreatorProcess());
                m_eventData->m_depositingProcess->push_back(hit->GetDepositingProcess());

                m_eventData->m_xp->push_back(hit->GetPosition().x()/mm);
                m_eventData->m_yp->push_back(hit->GetPosition().y()/mm);
                m_eventData->m_zp->push_back(hit->GetPosition().z()/mm);

                totalEnergyDeposited += hit->GetEnergyDeposited()/keV;
                m_eventData->m_energyDeposited->push_back(hit->GetEnergyDeposited()/keV);
                m_eventData->m_kineticEnergy->push_back(hit->GetKineticEnergy()/keV);
                m_eventData->m_time->push_back(hit->GetTime()/second);

                nbSteps++;
            }
        }
        m_eventData->m_nbSteps = nbSteps;
        m_eventData->m_totalEnergyDeposited = totalEnergyDeposited;

    }
    m_ttree->Fill();
    if(event->GetEventID() % m_saveInterval == 0)
        m_ttree->AutoSave("SaveSelf");

    m_eventData->Clear();
}

void
AnalysisManager::EndOfEventStoreTrajectory(const G4Event *event)
{
    m_eventData->m_primaryParticleType->push_back(m_primaryGeneratorAction->GetParticleTypeOfPrimary());
    m_eventData->m_primaryParticleEnergy = m_primaryGeneratorAction->GetEnergyOfPrimary()/keV;
    m_eventData->m_primaryXm = m_primaryGeneratorAction->GetMomentumDirectionOfPrimary().x();
    m_eventData->m_primaryYm = m_primaryGeneratorAction->GetMomentumDirectionOfPrimary().y();
    m_eventData->m_primaryZm = m_primaryGeneratorAction->GetMomentumDirectionOfPrimary().z();
    m_eventData->m_primaryX = m_primaryGeneratorAction->GetPositionOfPrimary().x()/mm;
    m_eventData->m_primaryY = m_primaryGeneratorAction->GetPositionOfPrimary().y()/mm;
    m_eventData->m_primaryZ = m_primaryGeneratorAction->GetPositionOfPrimary().z()/mm;

    m_eventData->m_eventId = event->GetEventID();

    G4Trajectory *trajectory;
    G4String particleName;
    G4int nbSteps, nbPointEntries;
    G4ThreeVector pointn, pointn_1, pdiff, initMomentum;
    G4double totalEnergyDeposited;

    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    G4int nbTrajectories = 0;
    if (trajectoryContainer) nbTrajectories = trajectoryContainer->entries();
    m_eventData->m_nbHits = nbTrajectories;

    nbSteps = 0;
    totalEnergyDeposited = 0.0;
    if (nbTrajectories > 0) {
        for(int j=0; j<nbTrajectories; j++) {
            trajectory = (G4Trajectory*)(*trajectoryContainer)[j];
            particleName = trajectory->GetParticleName();
            nbPointEntries = trajectory->GetPointEntries();
            initMomentum = trajectory->GetInitialMomentum();

            m_eventData->m_trackId->push_back(trajectory->GetTrackID());
            m_eventData->m_parentId->push_back(trajectory->GetParentID());

            m_eventData->m_particleType->push_back(particleName);
            m_eventData->m_parentType->push_back("");

            m_eventData->m_creatorProcess->push_back("");
            m_eventData->m_depositingProcess->push_back("");

            // `ed' becomes the norm of the initial momentum
            m_eventData->m_energyDeposited->push_back(trajectory->GetInitialMomentum().mag()/keV);
            // `ek' becomes initial kinetic energy
            m_eventData->m_kineticEnergy->push_back(trajectory->GetInitialKineticEnergy()/keV);
            totalEnergyDeposited += trajectory->GetInitialKineticEnergy()/keV;
            // `time' is used to record the index of new trajectory starting point in m_xp,yp,zp values
            m_eventData->m_time->push_back((G4double)nbSteps);

            for(int i=0; i<nbPointEntries; i++) {
                pointn = trajectory->GetPoint(i)->GetPosition();

                m_eventData->m_xp->push_back(pointn.x()/mm);
                m_eventData->m_yp->push_back(pointn.y()/mm);
                m_eventData->m_zp->push_back(pointn.z()/mm);

                nbSteps++;
            }
        }
    }
    m_eventData->m_nbSteps = nbSteps;
    m_eventData->m_totalEnergyDeposited = totalEnergyDeposited;

    m_ttree->Fill();
    if(event->GetEventID() % m_saveInterval == 0)
        m_ttree->AutoSave("SaveSelf");
    m_eventData->Clear();
}

void
AnalysisManager::Step(const G4Step *step)
{

}
/*
G4bool
AnalysisManager::FilterEvent(EventData *pEventData)
{
    G4double dEnergyDepositedSensitiveRegion = 0.;

    vector<float> *pX = pEventData->m_pX;
    vector<float> *pY = pEventData->m_pY;
    vector<float> *pZ = pEventData->m_pZ;
    vector<float> *pEnergyDeposited = pEventData->m_pEnergyDeposited;

    const G4double dDriftLength = DetectorConstruction::GetGeometryParameter("DriftLength");
    const G4double dRadius = DetectorConstruction::GetGeometryParameter("TeflonCylinderInnerRadius");

    for(G4int i=0; i<pEnergyDeposited->size(); i++)
    {
        if((*pZ)[i] < 0. && (*pZ)[i] > -dDriftLength && std::sqrt((*pX)[i]*(*pX)[i] + (*pY)[i]*(*pY)[i]) < dRadius)
            dEnergyDepositedSensitiveRegion += (*pEnergyDeposited)[i];
    }

//    if(dEnergyDepositedSensitiveRegion > 0. && dEnergyDepositedSensitiveRegion < 100.)
    if(dEnergyDepositedSensitiveRegion > 0.)
        return false;
    else
        return true;
}
*/
