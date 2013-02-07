#include <G4ios.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4Track.hh>
#include <G4Event.hh>
#include <G4VProcess.hh>
#include <G4StackManager.hh>

#include "AnalysisManager.hh"

#include "StackingAction.hh"

StackingAction::StackingAction(AnalysisManager *analysisManager)
{
    m_analysisManager = analysisManager;
}

StackingAction::~StackingAction()
{
}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track *aTrack)
{
    G4ClassificationOfNewTrack trackClassification = fUrgent;

    if(aTrack->GetDefinition()->GetParticleType() == "nucleus" &&
       !aTrack->GetDefinition()->GetPDGStable()) {
        if(aTrack->GetParentID() > 0 && aTrack->GetCreatorProcess()->GetProcessName()
           == "RadioactiveDecay")
            trackClassification = fPostpone;
    }

    return trackClassification;
}

void
StackingAction::NewStage()
{
}

void
StackingAction::PrepareNewEvent()
{ 
}
