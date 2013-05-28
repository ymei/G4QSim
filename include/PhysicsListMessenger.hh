#ifndef __PHYSICSLISTMESSENGER_H__
#define __PHYSICSLISTMESSENGER_H__

#include "G4UImessenger.hh"
#include "globals.hh"

class PhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsListMessenger: public G4UImessenger
{
public:
  
    PhysicsListMessenger(PhysicsList* );
    ~PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:
  
    PhysicsList*               m_physicsList;
    
    G4UIdirectory*             m_physDir;
    G4UIcmdWithADoubleAndUnit* m_gammaCutCmd;
    G4UIcmdWithADoubleAndUnit* m_electronCutCmd;
    G4UIcmdWithADoubleAndUnit* m_positronCutCmd;
    G4UIcmdWithADoubleAndUnit* m_allCutCmd;
    G4UIcmdWithAString*        m_addListCmd;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __PHYSICSLISTMESSENGER_H__ */
