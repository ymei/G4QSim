#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "ParticleSource.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4int particleSourceType)
{
    m_particleSourceType = particleSourceType;
    m_particleTypeOfPrimary = "";
    m_positionOfPrimary = G4ThreeVector(0., 0., 0.);
    m_momentumDirectionOfPrimary = G4ThreeVector(0., 0., 0.);
    m_energyOfPrimary = 0.;

    // `gps' that comes with Geant4
    m_gpsGun = new G4GeneralParticleSource();
    // custom defined
    m_particleGun = new ParticleSource();

    // just an example, usually set in macro
    m_gpsGun->SetParticlePosition(m_positionOfPrimary);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete m_gpsGun;
    delete m_particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (m_gpsGun->GetParticleDefinition() == G4Geantino::Geantino()) {
        G4int Z = 10, A = 24;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;

        G4ParticleDefinition* ion
            = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
        m_gpsGun->SetParticleDefinition(ion);
        m_gpsGun->SetParticleCharge(ionCharge);
    }

    //create vertex
    if(m_particleSourceType) {
        m_particleGun->GeneratePrimaryVertex(anEvent);
        m_particleTypeOfPrimary      = m_particleGun->GetParticleDefinition()->GetParticleName();
        m_positionOfPrimary          = m_particleGun->GetParticlePosition();
        m_momentumDirectionOfPrimary = m_particleGun->GetParticleMomentumDirection();
        m_energyOfPrimary            = m_particleGun->GetParticleEnergy();
    } else {
        m_gpsGun->GeneratePrimaryVertex(anEvent);
        m_particleTypeOfPrimary      = m_gpsGun->GetParticleDefinition()->GetParticleName();
        m_positionOfPrimary          = m_gpsGun->GetParticlePosition();
        m_momentumDirectionOfPrimary = m_gpsGun->GetParticleMomentumDirection();
        m_energyOfPrimary            = m_gpsGun->GetParticleEnergy();
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
