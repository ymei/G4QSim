#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* physicsList)
    :m_physicsList(physicsList)
{ 
    m_physDir = new G4UIdirectory("/G4QSim/phys/");
    m_physDir->SetGuidance("physics list commands");

    m_gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/G4QSim/phys/setGCut",this);  
    m_gammaCutCmd->SetGuidance("Set gamma cut.");
    m_gammaCutCmd->SetParameterName("Gcut",false);
    m_gammaCutCmd->SetUnitCategory("Length");
    m_gammaCutCmd->SetRange("Gcut>0.0");
    m_gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    m_electronCutCmd = new G4UIcmdWithADoubleAndUnit("/G4QSim/phys/setECut",this);  
    m_electronCutCmd->SetGuidance("Set electron cut.");
    m_electronCutCmd->SetParameterName("Ecut",false);
    m_electronCutCmd->SetUnitCategory("Length");
    m_electronCutCmd->SetRange("Ecut>0.0");
    m_electronCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
    m_positronCutCmd = new G4UIcmdWithADoubleAndUnit("/G4QSim/phys/setPCut",this);  
    m_positronCutCmd->SetGuidance("Set positron cut.");
    m_positronCutCmd->SetParameterName("Pcut",false);
    m_positronCutCmd->SetUnitCategory("Length");
    m_positronCutCmd->SetRange("Pcut>0.0");
    m_positronCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

    m_allCutCmd = new G4UIcmdWithADoubleAndUnit("/G4QSim/phys/setCuts",this);  
    m_allCutCmd->SetGuidance("Set cut for all.");
    m_allCutCmd->SetParameterName("cut",false);
    m_allCutCmd->SetUnitCategory("Length");
    m_allCutCmd->SetRange("cut>0.0");
    m_allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    m_addListCmd = new G4UIcmdWithAString("/G4QSim/phys/addPhysics",this);  
    m_addListCmd->SetGuidance("Add modular physics list.");
    m_addListCmd->SetParameterName("PList",false);
    // addPhysics must be used before physics list initialization
    m_addListCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
    delete m_gammaCutCmd;
    delete m_electronCutCmd;
    delete m_positronCutCmd;
    delete m_allCutCmd;
    delete m_addListCmd;
    delete m_physDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                       G4String newValue)
{       
    if( command == m_gammaCutCmd )
    { m_physicsList->SetCutForGamma(m_gammaCutCmd->GetNewDoubleValue(newValue));}
     
    if( command == m_electronCutCmd )
    { m_physicsList->SetCutForElectron(m_electronCutCmd->GetNewDoubleValue(newValue));}
     
    if( command == m_positronCutCmd )
    { m_physicsList->SetCutForPositron(m_positronCutCmd->GetNewDoubleValue(newValue));}

    if( command == m_allCutCmd )
    {
        G4double cut = m_allCutCmd->GetNewDoubleValue(newValue);
        m_physicsList->SetCutForGamma(cut);
        m_physicsList->SetCutForElectron(cut);
        m_physicsList->SetCutForPositron(cut);
    }
    
    if( command == m_addListCmd )
    { m_physicsList->AddPhysicsList(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
