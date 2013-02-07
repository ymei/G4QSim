#ifndef __HIT_H__
#define __HIT_H__

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
class SDHit : public G4VHit
{
    
public:

    SDHit();
    SDHit(const SDHit&);
    
    virtual ~SDHit();

    // operators
    const SDHit& operator=(const SDHit&);
    G4int operator==(const SDHit&) const;
    
    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {;}
    virtual void Print();

    // methods to handle data
    void SetTrackId(G4int trackId) { m_trackId = trackId; }
    void SetParentId(G4int parentId) { m_parentId = parentId; }
    void SetParticleType(const G4String &particleType) { m_particleType = new G4String(particleType); }
    void SetParentType(const G4String &parentType) { m_parentType = new G4String(parentType); }
    void SetCreatorProcess(const G4String &creatorProcess) { m_creatorProcess = new G4String(creatorProcess); }
    void SetDepositingProcess(const G4String &depositingProcess) { m_depositingProcess = new G4String(depositingProcess); }
    void SetPosition(G4ThreeVector position) { m_position = position; }
    void SetEnergyDeposited(G4double energyDeposited) { m_energyDeposited = energyDeposited; }
    void SetKineticEnergy(G4double kineticEnergy) { m_kineticEnergy = kineticEnergy; }
    void SetTime(G4double time) { m_time = time; }

    // get methods
    G4int GetTrackId() const { return m_trackId; }
    G4int GetParentId() const { return m_parentId; }
    G4String &GetParticleType() const { return *m_particleType; }
    G4String &GetParentType() const { return *m_parentType; }
    G4String &GetCreatorProcess() const { return *m_creatorProcess; }
    G4String &GetDepositingProcess() const { return *m_depositingProcess; }
    G4ThreeVector GetPosition() const { return m_position; }
    G4double GetEnergyDeposited() const { return m_energyDeposited; }
    G4double GetKineticEnergy() const { return m_kineticEnergy; }
    G4double GetTime() const { return m_time; }

private:
    G4int m_trackId;
    G4int m_parentId;
    G4String *m_particleType;
    G4String *m_parentType;
    G4String *m_creatorProcess;
    G4String *m_depositingProcess;
    G4ThreeVector m_position;
    G4double m_energyDeposited;
    G4double m_kineticEnergy;
    G4double m_time;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<SDHit> SDHitsCollection;
extern G4Allocator<SDHit> SDHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* SDHit::operator new(size_t)
{
    void *hit;
    hit = (void *)SDHitAllocator.MallocSingle();
    return hit;
}

inline void SDHit::operator delete(void *hit)
{
    SDHitAllocator.FreeSingle((SDHit*)hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __HIT_H__ */
