#ifndef __PARTICLESOURCEMESSENGER_H__
#define __PARTICLESOURCEMESSENGER_H__

#include <G4UImessenger.hh>
#include <globals.hh>

class ParticleSource;

class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ParticleSourceMessenger: public G4UImessenger
{
public:
    ParticleSourceMessenger(ParticleSource *particleSource);
    ~ParticleSourceMessenger();
    /// Called by UI to set new values.
    void SetNewValue(G4UIcommand *command, G4String newValues);
 
private:
    G4ParticleDefinition *ParseIonValues(G4String newValues);

    ParticleSource *m_particleSource; ///< ParticleSource passed into this class.
    G4ParticleTable *m_particleTable;
    G4ParticleDefinition *m_ionDef;
    G4double m_ionCharge;

    // need to be deleted in destructor
    G4UIdirectory *m_uiDirectory;
    // commands
    G4UIcmdWithoutParameter *m_listCmd;
    G4UIcmdWithAString *m_typeCmd;
    G4UIcmdWithADoubleAndUnit *m_QValueCmd;
    G4UIcmdWithAString *m_EHistCmd;
    G4UIcmdWithAString *m_dBDEventsCmd;
    G4UIcmdWith3Vector *m_fIonPDirCmd;
    G4UIcmdWith3VectorAndUnit *m_fIonPosCmd;
    G4UIcmdWithADoubleAndUnit *m_fIonEkCmd;
    G4UIcommand *m_fIonCmd, *m_dIonCmd;
    G4UIcmdWithADouble *m_beNu_aCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __PARTICLESOURCEMESSENGER_H__ */
