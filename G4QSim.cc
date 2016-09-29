/** \file
 * G4QSim main entry.
 *
 * See print_usage() for available arguments.
 */
#include <string>
#include <sstream>
#include <unistd.h>

#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "SensitiveDetector.hh"
#include "StackingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4GDMLParser.hh"

/** print arguments and exit.
 */
void print_usage(void)
{
    G4cerr << "Usage:" << G4endl;
    G4cerr << "       -g [geometry description file name, .gdml]"       << G4endl;
    G4cerr << "       -i [tcsh | qt], interactive with UI type"         << G4endl;
    G4cerr << "       -m [macro (input) file name]"                     << G4endl;
    G4cerr << "       -n [#] of events to simulate"                     << G4endl;
    G4cerr << "       -o [output (data) file name]"                     << G4endl;
    G4cerr << "       -s [random number seed]"                          << G4endl;
    G4cerr << "       -u [particleSource, 0:gps, 1:(double)BetaDecay]"  << G4endl;
    G4cerr << "       -v [vrml | opengl], visual engine selection"      << G4endl;
}

/** main entry.
 */
int
main(int argc, char **argv)
{
    // switches
    int optC = 0;

    std::stringstream optStrStream;
    G4String geometryFilename, macroFilename, dataFilename, uiTypeName, cmdString;
    
    G4bool interactiveQ = false;
    G4bool visualizeQ = false;
    G4bool vrmlVisualizeQ = false;
    G4bool openGLVisualizeQ = false;
    G4int nbEventsToSimulate = 0;
    G4int particleSourceType = 0;

    G4GDMLParser gdmlParser;
    G4long randSeed = time(NULL);

    // parse switches
    while((optC = getopt(argc, argv, "g:i:m:n:o:s:u:v:")) != -1)
    {
        switch(optC)
        {
        case 'g':
            geometryFilename = optarg;
            break;
            
        case 'i':
            interactiveQ = true;
            uiTypeName = optarg;
            break;

        case 'm':
            macroFilename = optarg;
            break;

        case 'n':
            optStrStream.str(optarg);
            optStrStream.clear(); // clear `error control' state, not the content
            optStrStream >> nbEventsToSimulate;
            break;

        case 'o':
            dataFilename = optarg;
            break;

        case 's':
            optStrStream.str(optarg);
            optStrStream.clear();
            optStrStream >> randSeed;
            break;

        case 'u':
            optStrStream.str(optarg);
            optStrStream.clear();
            optStrStream >> particleSourceType;

        case 'v':
            visualizeQ = true;
            optStrStream.str(optarg);
            if(optStrStream.str() == "vrml")
                vrmlVisualizeQ = true;
            else if(optStrStream.str() == "opengl")
                openGLVisualizeQ = true;
            break;

        default:
            print_usage();
            exit(EXIT_FAILURE);
            break;
        }
    }

    // read and setup the geometry
    if(geometryFilename.empty()) {
        G4cerr << "Geometry File not found!" << G4endl;
        print_usage();
        exit(EXIT_FAILURE);
    } else {
        gdmlParser.Read(geometryFilename);
    }

    // Random Number Generator Initialization
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(randSeed);
    CLHEP::HepRandom::showEngineStatus();

    // my Verbose output class
    G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    
    // create the run manager
    G4RunManager *runManager = new G4RunManager;

    // set user-defined initialization classes
    runManager->SetUserInitialization(new DetectorConstruction(&gdmlParser));
    runManager->SetUserInitialization(new PhysicsList);
    /*
    G4VModularPhysicsList* physicsList = new LHEP_EMV;
    runManager->SetUserInitialization(physicsList);
    */

    // create the primary generator action
    PrimaryGeneratorAction *primaryGeneratorAction = new PrimaryGeneratorAction(particleSourceType);
    // create an analysis manager object
    AnalysisManager *analysisManager = new AnalysisManager(primaryGeneratorAction);
    // set user-defined action classes
    runManager->SetUserAction(primaryGeneratorAction);
    runManager->SetUserAction(new StackingAction(analysisManager));
    RunAction* runAction = new RunAction(primaryGeneratorAction, analysisManager);
    runManager->SetUserAction(runAction);
    EventAction *eventAction = new EventAction(analysisManager);
    runManager->SetUserAction(eventAction);
    runManager->SetUserAction(new TrackingAction(runAction, eventAction));

    runManager->Initialize();
    G4UImanager* uiManager = G4UImanager::GetUIpointer();

    if(!dataFilename.empty()) {
        analysisManager->SetDataFilename(dataFilename);
    }
    if(nbEventsToSimulate) {
        analysisManager->SetNbEventsToSimulate(nbEventsToSimulate);
    }
 
#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif

    G4UIExecutive* uiExec = NULL;
    if(interactiveQ) {
#ifdef G4UI_USE
        if(!uiTypeName.empty())
            uiExec = new G4UIExecutive(argc, argv, uiTypeName);
        else // default UI type is tcsh
            uiExec = new G4UIExecutive(argc, argv, "tcsh");
#endif
    }
    if(visualizeQ) {
#ifdef G4VIS_USE
        uiManager->ApplyCommand("/control/verbose 2");
        uiManager->ApplyCommand("/run/verbose 2");
        uiManager->ApplyCommand("/run/initialize");

        if(vrmlVisualizeQ)
            uiManager->ApplyCommand("/vis/open VRML2FILE");
        if(openGLVisualizeQ)
            uiManager->ApplyCommand("/vis/open OGL");

        uiManager->ApplyCommand("/vis/viewer/set/autoRefresh false");
        uiManager->ApplyCommand("/vis/verbose errors");

        uiManager->ApplyCommand("/vis/drawVolume");
        uiManager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 30 -30 deg");
        uiManager->ApplyCommand("/vis/viewer/zoom 1.0");
        uiManager->ApplyCommand("/vis/viewer/set/style wireframe");

        // uiManager->ApplyCommand("/vis/scene/add/axes -1 -1 1 0.3 m");
        uiManager->ApplyCommand("/vis/scene/add/trajectories smooth");
        uiManager->ApplyCommand("/vis/scene/add/hits");
        uiManager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
        
        uiManager->ApplyCommand("/vis/viewer/set/autoRefresh true");
        uiManager->ApplyCommand("/vis/verbose warnings");
        
        uiManager->ApplyCommand("/tracking/verbose 2");
        uiManager->ApplyCommand("/tracking/storeTrajectory 1");

        uiManager->ApplyCommand("/vis/viewer/flush");
#endif
    }

    if(!macroFilename.empty()) {
        cmdString = "/control/execute " + macroFilename;
        uiManager->ApplyCommand(cmdString);
    }
    if(nbEventsToSimulate) {
        optStrStream.str("");
        optStrStream.clear();
        optStrStream << "/run/beamOn " << nbEventsToSimulate;
        uiManager->ApplyCommand(optStrStream.str());
    }

    if(interactiveQ) {
#ifdef G4UI_USE
        uiExec->SessionStart();
        delete uiExec;
#endif
    }

    delete analysisManager;
#ifdef G4VIS_USE
    if(visualizeQ)
        delete visManager;
#endif
    delete runManager;

    return EXIT_SUCCESS;
}
