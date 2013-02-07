#ifndef __STACKINGACTION_H__
#define __STACKINGACTION_H__

#include <globals.hh>
#include <G4UserStackingAction.hh>

class AnalysisManager;

class StackingAction: public G4UserStackingAction
{
public:
    StackingAction(AnalysisManager *analysisManager=NULL);
    ~StackingAction();

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();

private:
    AnalysisManager *m_analysisManager;
};

#endif /* __STACKINGACTION_H__ */
