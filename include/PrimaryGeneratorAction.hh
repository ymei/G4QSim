#ifndef __PRIMARYGENERATORACTION_H__
#define __PRIMARYGENERATORACTION_H__

#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();    
    ~PrimaryGeneratorAction();

public:
    void GeneratePrimaries(G4Event*);
    G4GeneralParticleSource* GetParticleGun() { return m_particleGun;}

    const G4String &GetParticleTypeOfPrimary() { return m_particleTypeOfPrimary; }
    G4ThreeVector GetPositionOfPrimary() const { return m_positionOfPrimary; }
    G4ThreeVector GetMomentumDirectionOfPrimary() const { return m_momentumDirectionOfPrimary; }
    G4double GetEnergyOfPrimary() const { return m_energyOfPrimary; }

private:
    G4GeneralParticleSource *m_particleGun;

    G4String m_particleTypeOfPrimary;
    G4ThreeVector m_positionOfPrimary;
    G4ThreeVector m_momentumDirectionOfPrimary;
    G4double m_energyOfPrimary;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __PRIMARYGENERATORACTION_H__ */
