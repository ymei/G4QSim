#ifndef __PHYSICSLIST_H__
#define __PHYSICSLIST_H__

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VUserPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList();

public:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess(); 

    void SetCuts();

    //these methods Construct particles
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructHadrons();

    //these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void ConstructOp();
    
    //for the Messenger 
    void SetVerbose(G4int);
    void SetNbOfPhotonsCerenkov(G4int);

private:
    G4Cerenkov*          m_theCerenkovProcess;
    G4Scintillation*     m_theScintillationProcess;
    G4OpAbsorption*      m_theAbsorptionProcess;
    G4OpRayleigh*        m_theRayleighScatteringProcess;
    G4OpMieHG*           m_theMieHGScatteringProcess;
    G4OpBoundaryProcess* m_theBoundaryProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* __PHYSICSLIST_H__ */
