#include "G4SDManager.hh"

#include <map>
using std::map;

#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"

DetectorConstruction::DetectorConstruction(G4GDMLParser *gdmlParser)
{   
    m_world = NULL;
    m_gdmlParser = gdmlParser;
    if(m_gdmlParser)
        m_world = m_gdmlParser->GetWorldVolume();
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume * DetectorConstruction::Construct()
{
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4String trackerChamberSDname = "Tracker";
    SensitiveDetector* aTrackerSD = new SensitiveDetector(trackerChamberSDname);
    sdManager->AddNewDetector(aTrackerSD);

    ///////////////////////////////////////////////////////////////////////
    //
    // Retrieve Auxiliary Information for sensitive detector
    //
    const G4GDMLAuxMapType* auxmap = m_gdmlParser->GetAuxMap();
    
    G4cout << "Found " << auxmap->size()
           << " volume(s) with auxiliary information."
           << G4endl;
    
    // now we are looking for sensitive detectors
    // and setting them for the volumes

    for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
        iter!=auxmap->end(); iter++) {
        
        G4cout << "Volume ``" << ((*iter).first)->GetName()
               << "'' has the following list of auxiliary information: "
               << G4endl;
        
        for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
             vit!=(*iter).second.end();vit++) {

            if ((*vit).type=="SensDet") {
                G4cout << "Attaching sensitive detector ``" << (*vit).value
                       << "'' to volume ``" << ((*iter).first)->GetName()
                       << "''" <<  G4endl;

                G4VSensitiveDetector* mydet = sdManager->FindSensitiveDetector((*vit).value);
                if(mydet) {
                    G4LogicalVolume* myvol = (*iter).first;
                    myvol->SetSensitiveDetector(mydet);
                } else {
                    G4cout << (*vit).value << " detector not found" << G4endl;
                }
            }
        }
    }

    //
    // End of Auxiliary Information block
    //
    ////////////////////////////////////////////////////////////////////////

    return m_world;
}
