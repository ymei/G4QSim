#ifndef __SENSITIVEDETECTOR_H__
#define __SENSITIVEDETECTOR_H__

#include "G4VSensitiveDetector.hh"
#include "SDHit.hh"

class G4Step;

class SensitiveDetector : public G4VSensitiveDetector
{
public:
    SensitiveDetector(const G4String&);
    ~SensitiveDetector();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

private:

    SDHitsCollection* m_hitsCollection;
    std::map<G4int,G4String> m_particleTypes;
};

#endif /* __SENSITIVEDETECTOR_H__ */
