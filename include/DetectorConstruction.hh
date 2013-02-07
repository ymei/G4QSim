#ifndef _DETECTORCONSTRUCTION_H_
#define _DETECTORCONSTRUCTION_H_

#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

    DetectorConstruction(G4GDMLParser *gdmlParser = NULL);
    ~DetectorConstruction();

    G4VPhysicalVolume *Construct();

private:

    G4GDMLParser *m_gdmlParser;
    G4VPhysicalVolume *m_world;
};

#endif /* _DETECTORCONSTRUCTION_H_ */
