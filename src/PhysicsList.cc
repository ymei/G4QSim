#include "globals.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4ProcessManager.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLivermorePhysics.hh"

#include "G4VModularPhysicsList.hh"

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
    : G4VUserPhysicsList()
{
    //add new units for radioActive decays
    //
    const G4double minute = 60*second;
    const G4double hour   = 60*minute;
    const G4double day    = 24*hour;
    const G4double year   = 365*day;
    new G4UnitDefinition("minute", "min", "Time", minute);
    new G4UnitDefinition("hour",   "h",   "Time", hour);
    new G4UnitDefinition("day",    "d",   "Time", day);
    new G4UnitDefinition("year",   "y",   "Time", year);

    m_theCerenkovProcess           = NULL;
    m_theScintillationProcess      = NULL;
    m_theAbsorptionProcess         = NULL;
    m_theRayleighScatteringProcess = NULL;
    m_theMieHGScatteringProcess    = NULL;
    m_theBoundaryProcess           = NULL;
    m_theParticleIterator          = theParticleTable->GetIterator();

    m_currentDefaultCut = 1.0*mm;
    m_cutForGamma = m_currentDefaultCut;
    m_cutForElectron = m_currentDefaultCut;
    m_cutForPositron = m_currentDefaultCut;

    m_physicsListMessenger = new PhysicsListMessenger(this);
    m_emName = G4String("EmPenelope");
    m_emPhysicsList = new G4EmPenelopePhysics();

    SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
    // In this method, static member functions should be called
    // for all particles which you want to use.
    // This ensures that objects of these particle types will be
    // created in the program.

    ConstructBosons();
    ConstructLeptons();
    ConstructHadrons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructBosons()
{
    // pseudo-particles
    G4Geantino::GeantinoDefinition();
    G4ChargedGeantino::ChargedGeantinoDefinition();

    // gamma
    G4Gamma::GammaDefinition();

    // optical photon
    G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructLeptons()
{
    // leptons
    //  e+/-
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    // mu+/-
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
    // nu_e
    G4NeutrinoE::NeutrinoEDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    // nu_mu
    G4NeutrinoMu::NeutrinoMuDefinition();
    G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructHadrons()
{
    //  mesons
    G4MesonConstructor mConstructor;
    mConstructor.ConstructParticle();

    //  baryons
    G4BaryonConstructor bConstructor;
    bConstructor.ConstructParticle();

    //  ions
    G4IonConstructor iConstructor;
    iConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
    AddTransportation();

    ConstructGeneral();
    // ConstructEM();
    m_emPhysicsList->ConstructProcess();
    ConstructOp();
    // Optical Physics
    /*
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    RegisterPhysics(opticalPhysics);
    opticalPhysics->SetMaxNumPhotonsPerStep(1000);
    opticalPhysics->SetMaxBetaChangePerStep(10.0);
    opticalPhysics->SetTrackSecondariesFirst(true);
    */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String &name)
{
    if (verboseLevel>0) {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }

    if (name == m_emName) {
        return;
    } else if (name == "EmLivermore") {
        m_emName = name;
        delete m_emPhysicsList;
        m_emPhysicsList = new G4EmLivermorePhysics();
    } else {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
               << " is not defined" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructGeneral()
{
    // Add Decay Process
    G4Decay* m_theDecayProcess = new G4Decay();

    m_theParticleIterator->reset();
    while( (*m_theParticleIterator)() ){
        G4ParticleDefinition* particle = m_theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if (m_theDecayProcess->IsApplicable(*particle)) {
            pmanager ->AddProcess(m_theDecayProcess);
            // set ordering for PostStepDoIt and AtRestDoIt
            pmanager ->SetProcessOrdering(m_theDecayProcess, idxPostStep);
            pmanager ->SetProcessOrdering(m_theDecayProcess, idxAtRest);
        }
    }

    // Add radioactive decay
    G4RadioactiveDecay* m_theRadioactiveDecay = new G4RadioactiveDecay();
    m_theRadioactiveDecay->SetHLThreshold(-1.*s);
    m_theRadioactiveDecay->SetICM(true);
    m_theRadioactiveDecay->SetARM(false);

    const G4IonTable *m_theIonTable =
        G4ParticleTable::GetParticleTable()->GetIonTable();
    for (G4int i=0; i<m_theIonTable->Entries(); i++)
    {
        G4String particleName = m_theIonTable->GetParticle(i)->GetParticleName();
        G4String particleType = m_theIonTable->GetParticle(i)->GetParticleType();

        G4cout << __FUNCTION__ << "i=" << i << " " << particleName << G4endl;

        if (particleName == "GenericIon")
        {
            G4ProcessManager* pmanager =
                m_theIonTable->GetParticle(i)->GetProcessManager();
            pmanager->SetVerboseLevel(0);
            pmanager->AddProcess(m_theRadioactiveDecay);
            pmanager->SetProcessOrdering(m_theRadioactiveDecay, idxPostStep);
            pmanager->SetProcessOrdering(m_theRadioactiveDecay, idxAtRest);
        }
    }

//    G4ProcessManager* pmanager = G4GenericIon::GenericIon()->GetProcessManager();
//    pmanager->AddProcess(m_theRadioactiveDecay, 0, -1, 1);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RayleighScattering.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4PenelopeAnnihilationModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PenelopeComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeRayleighModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructEM()
{
    m_theParticleIterator->reset();
    while( (*m_theParticleIterator)() ){
        G4ParticleDefinition* particle = m_theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (particleName == "gamma") {
            // gamma
            // Construct processes for gamma
            /*
            pmanager->AddDiscreteProcess(new G4GammaConversion());
            pmanager->AddDiscreteProcess(new G4ComptonScattering());
            pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
            */
            // using Penelope model
            G4RayleighScattering *m_theRayleighScattering = new G4RayleighScattering();
            m_theRayleighScattering->SetEmModel(new G4PenelopeRayleighModel());
            pmanager->AddDiscreteProcess(m_theRayleighScattering);

            G4ComptonScattering *m_theComptonScattering = new G4ComptonScattering();
            m_theComptonScattering->SetEmModel(new G4PenelopeComptonModel());
            pmanager->AddDiscreteProcess(m_theComptonScattering);

            G4PhotoElectricEffect *m_thePhotoElectricEffect = new G4PhotoElectricEffect();
            m_thePhotoElectricEffect->SetEmModel(new G4PenelopePhotoElectricModel());
            pmanager->AddDiscreteProcess(m_thePhotoElectricEffect);

            G4GammaConversion *m_theGammaConversion = new G4GammaConversion();
            m_theGammaConversion->SetEmModel(new G4PenelopeGammaConversionModel());
            pmanager->AddDiscreteProcess(m_theGammaConversion);

        } else if (particleName == "e-") {
            //electron
            // Construct processes for electron
            pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
            /*
            pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);
            */
            G4eIonisation *m_theIonisation = new G4eIonisation();
            m_theIonisation->SetEmModel(new G4PenelopeIonisationModel());
            pmanager->AddProcess(m_theIonisation, -1, 2, 2);

            G4eBremsstrahlung *m_theBremsstrahlung = new G4eBremsstrahlung();
            m_theBremsstrahlung->SetEmModel(new G4PenelopeBremsstrahlungModel());
            pmanager->AddProcess(m_theBremsstrahlung, -1, -3, 3);

        } else if (particleName == "e+") {
            //positron
            // Construct processes for positron
            pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
            /*
            pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3, 3);
            pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1, 4);
            */
            G4eIonisation *m_theIonisation = new G4eIonisation();
            m_theIonisation->SetEmModel(new G4PenelopeIonisationModel());
            pmanager->AddProcess(m_theIonisation, -1, 2, 2);

            G4eBremsstrahlung *m_theBremsstrahlung = new G4eBremsstrahlung();
            m_theBremsstrahlung->SetEmModel(new G4PenelopeBremsstrahlungModel());
            pmanager->AddProcess(m_theBremsstrahlung, -1, -3, 3);

            G4eplusAnnihilation *m_theplusAnnihilation = new G4eplusAnnihilation();
            m_theplusAnnihilation->SetEmModel(new G4PenelopeAnnihilationModel());
            pmanager->AddDiscreteProcess(m_theplusAnnihilation);

        } else if( particleName == "mu+" ||
                   particleName == "mu-"    ) {
            //muon
            // Construct processes for muon
            pmanager->AddProcess(new G4MuMultipleScattering(),-1, 1, 1);
            pmanager->AddProcess(new G4MuIonisation(),      -1, 2, 2);
            pmanager->AddProcess(new G4MuBremsstrahlung(),  -1, 3, 3);
            pmanager->AddProcess(new G4MuPairProduction(),  -1, 4, 4);

        } else {
            if ((particle->GetPDGCharge() != 0.0) &&
                (particle->GetParticleName() != "chargedgeantino")) {
                // all others charged particles except geantino
                pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
                pmanager->AddProcess(new G4hIonisation(),       -1,2,2);
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructOp()
{
    m_theCerenkovProcess           = new G4Cerenkov("Cerenkov");
    m_theScintillationProcess      = new G4Scintillation("Scintillation");
    m_theAbsorptionProcess         = new G4OpAbsorption();
    m_theRayleighScatteringProcess = new G4OpRayleigh();
    m_theMieHGScatteringProcess    = new G4OpMieHG();
    m_theBoundaryProcess           = new G4OpBoundaryProcess();

//  m_theCerenkovProcess->DumpPhysicsTable();
//  m_theScintillationProcess->DumpPhysicsTable();
//  m_theRayleighScatteringProcess->DumpPhysicsTable();

    SetVerbose(0);

    m_theCerenkovProcess->SetMaxNumPhotonsPerStep(1000);
    m_theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    m_theCerenkovProcess->SetTrackSecondariesFirst(true);

    m_theScintillationProcess->SetScintillationYieldFactor(1.);
    m_theScintillationProcess->SetTrackSecondariesFirst(true);

    // Use Birks Correction in the Scintillation process

    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    m_theScintillationProcess->AddSaturation(emSaturation);

    // obsolete
    // G4OpticalSurfaceModel themodel = unified;
    // theBoundaryProcess->SetModel(themodel);

    m_theParticleIterator->reset();
    while( (*m_theParticleIterator)() ){
        G4ParticleDefinition* particle = m_theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (m_theCerenkovProcess->IsApplicable(*particle)) {
            //pmanager->AddProcess(m_theCerenkovProcess);
            //pmanager->SetProcessOrdering(m_theCerenkovProcess,idxPostStep);
        }
        /* Disable scintillation
        if (m_theScintillationProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(m_theScintillationProcess);
            pmanager->SetProcessOrderingToLast(m_theScintillationProcess, idxAtRest);
            pmanager->SetProcessOrderingToLast(m_theScintillationProcess, idxPostStep);
        }
        */
        if (particleName == "opticalphoton") {
            G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
            pmanager->AddDiscreteProcess(m_theAbsorptionProcess);
            // pmanager->AddDiscreteProcess(m_theRayleighScatteringProcess);
            // pmanager->AddDiscreteProcess(m_theMieHGScatteringProcess);
            pmanager->AddDiscreteProcess(m_theBoundaryProcess);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVerbose(G4int verbose)
{
    m_theCerenkovProcess->SetVerboseLevel(verbose);
    m_theScintillationProcess->SetVerboseLevel(verbose);
    m_theAbsorptionProcess->SetVerboseLevel(verbose);
    m_theRayleighScatteringProcess->SetVerboseLevel(verbose);
    m_theMieHGScatteringProcess->SetVerboseLevel(verbose);
    m_theBoundaryProcess->SetVerboseLevel(verbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
    m_theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumber);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
    // fixe lower limit for cut
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);

    // "G4VUserPhysicsList::SetCutsWithDefault" method sets
    // the default cut value for all particle types
    // SetCutsWithDefault();

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    SetCutValue(m_cutForGamma, "gamma");
    SetCutValue(m_cutForElectron, "e-");
    SetCutValue(m_cutForPositron, "e+");

    DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
    m_cutForGamma = cut;
    SetParticleCuts(m_cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
    m_cutForElectron = cut;
    SetParticleCuts(m_cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
    m_cutForPositron = cut;
    SetParticleCuts(m_cutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
