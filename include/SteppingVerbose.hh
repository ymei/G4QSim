#ifndef __STEPPINGVERBOSE_H__
#define __STEPPINGVERBOSE_H__

#include "G4SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingVerbose : public G4SteppingVerbose {

public:   

    SteppingVerbose();
    ~SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __STEPPINGVERBOSE_H__ */
