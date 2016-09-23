#ifndef __PARTICLESOURCE_H__
#define __PARTICLESOURCE_H__

#include <G4VPrimaryGenerator.hh>
#include <G4Navigator.hh>
#include <G4ParticleMomentum.hh>
#include <G4ParticleDefinition.hh>
#include <G4Track.hh>
#include <Randomize.hh>

#include "ParticleSourceMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ParticleSource: public G4VPrimaryGenerator
{
public:
    ParticleSource();
    ~ParticleSource();

public:
    void GeneratePrimaryVertex(G4Event *event);

    void SetGunType(G4String gunType) { m_gunType = gunType; }
    void SetQValue(G4double QValue) { m_QValue = QValue; }
    int LoadEHist(G4String EHistFName);
    int LoadDBDEvents(G4String dBDEventsFName);
    void SetFIonPDir(G4ThreeVector fIonPDir) { m_fIonPDir = fIonPDir; }
    void SetFIonPos(G4ThreeVector fIonPos) { m_fIonPos = fIonPos; }
    void SetFIonEk(G4double fIonEk) { m_fIonEk = fIonEk; }
    void SetFIon(G4ParticleDefinition *ionDef, G4double ionCharge)
        { m_fIonDef = ionDef; m_fIonCharge = ionCharge; }
    void SetDIon(G4ParticleDefinition *ionDef, G4double ionCharge)
        { m_dIonDef = ionDef; m_dIonCharge = ionCharge; }
    void SetBeNu_a(G4double beNu_a) { m_beNu_a = beNu_a;}

    G4ParticleDefinition *GetParticleDefinition(){return m_fIonDef;}
    G4ThreeVector GetParticlePosition() {return m_fIonPos;}
    G4ThreeVector GetParticleMomentumDirection() {return m_fIonPDir;}
    G4double GetParticleEnergy() {return m_fIonEk;}

    size_t m_dBDi; /// index of the list of DBD events

private:
    ParticleSourceMessenger *m_sourceMessenger;

    G4String m_gunType;
    G4double m_QValue;
    G4double m_beNu_a;
    G4ThreeVector m_fIonPos;
    G4double m_fIonEk;
    G4ThreeVector m_fIonPDir; /// father ion momentum direction
    G4ParticleDefinition *m_fIonDef, *m_dIonDef;
    G4double m_fIonCharge, m_dIonCharge;
    size_t m_EHistN;
    G4double *m_EHistE, *m_EHistP;
    G4RandGeneral *m_EGen;
    size_t m_dBDN;
    G4ThreeVector *m_dBDb[2]; /// momenta of the two betas from DBD
    G4double m_nuE;
    G4ThreeVector m_nuP;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __PARTICLESOURCE_H__ */
