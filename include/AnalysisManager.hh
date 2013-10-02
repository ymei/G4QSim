#ifndef __ANALYSISMANAGER_H__
#define __ANALYSISMANAGER_H__

#include <globals.hh>

#include <TParameter.h>

class G4Run;
class G4Event;
class G4Step;

class TFile;
class TTree;

class PrimaryGeneratorAction;
class EventData;

class AnalysisManager
{
public:
    AnalysisManager(PrimaryGeneratorAction *primaryGeneratorAction);
    virtual ~AnalysisManager();

public:
    virtual void BeginOfRun(const G4Run *run);
    virtual void EndOfRun(const G4Run *run);
    virtual void BeginOfEvent(const G4Event *event);
    virtual void EndOfEvent(const G4Event *event);
    virtual void EndOfEventStoreTrajectory(const G4Event *event);
    virtual void Step(const G4Step *step);

    void SetDataFilename(const G4String &filename) { m_dataFilename = filename; }
    void SetNbEventsToSimulate(G4int nbEventsToSimulate)
        { m_nbEventsToSimulate = nbEventsToSimulate; }
    void SetSaveInteval(G4int val=100000) { m_saveInterval = val; }

private:
//    G4bool FilterEvent(EventData *pEventData);

private:
    G4int m_hitsCollectionID;

    G4String m_dataFilename;
    G4int m_nbEventsToSimulate;

    TFile *m_tfile;
    TTree *m_ttree;
    TParameter<int> *m_nbEventsToSimulateParameter;
    G4int m_saveInterval;

    PrimaryGeneratorAction *m_primaryGeneratorAction;
    EventData *m_eventData;
};

#endif // __ANALYSISMANAGER_H__
