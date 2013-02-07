#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    m_particleTypeOfPrimary = "";
    m_positionOfPrimary = G4ThreeVector(0., 0., 0.);
    m_momentumDirectionOfPrimary = G4ThreeVector(0., 0., 0.);
    m_energyOfPrimary = 0.;

    // G4int nParticle = 1;
    // m_particleGun  = new G4ParticleGun(nParticle);
    m_particleGun  = new G4GeneralParticleSource();

    //m_particleGun->SetParticleEnergy(3.0*eV);
    m_particleGun->SetParticlePosition(m_positionOfPrimary);
    //m_particleGun->SetParticleMomentumDirection(G4ThreeVector(1.0, 0.0, 0.0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete m_particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (m_particleGun->GetParticleDefinition() == G4Geantino::Geantino()) {
        G4int Z = 10, A = 24;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;
    
        G4ParticleDefinition* ion
            = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
        m_particleGun->SetParticleDefinition(ion);
        m_particleGun->SetParticleCharge(ionCharge);
    }  
  
    //create vertex
    //   
    m_particleGun->GeneratePrimaryVertex(anEvent);

    m_particleTypeOfPrimary = m_particleGun->GetParticleDefinition()->GetParticleName();
    m_positionOfPrimary = m_particleGun->GetParticlePosition();
    m_momentumDirectionOfPrimary = m_particleGun->GetParticleMomentumDirection();
    m_energyOfPrimary = m_particleGun->GetParticleEnergy();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
