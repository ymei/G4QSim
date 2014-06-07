#include "SDHit.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4Allocator<SDHit> SDHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDHit::SDHit()
    : G4VHit(),
      m_trackId(-1), m_parentId(-1), m_energyDeposited(0.), m_kineticEnergy(0.), m_time(0.)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDHit::~SDHit() 
{
    if(m_particleType) delete m_particleType;
    if(m_parentType) delete m_parentType;
    if(m_creatorProcess) delete m_creatorProcess;
    if(m_depositingProcess) delete m_depositingProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDHit::SDHit(const SDHit& right)
    : G4VHit()
{
    m_trackId           = right.m_trackId;
    m_parentId          = right.m_parentId;
    m_particleType      = right.m_particleType;
    m_parentType        = right.m_parentType;
    m_creatorProcess    = right.m_creatorProcess;
    m_depositingProcess = right.m_depositingProcess;
    m_position          = right.m_position;
    m_energyDeposited   = right.m_energyDeposited;
    m_kineticEnergy     = right.m_kineticEnergy;
    m_time              = right.m_time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const SDHit& SDHit::operator=(const SDHit& right)
{
    m_trackId           = right.m_trackId;
    m_parentId          = right.m_parentId;
    m_particleType      = right.m_particleType;
    m_parentType        = right.m_parentType;
    m_creatorProcess    = right.m_creatorProcess;
    m_depositingProcess = right.m_depositingProcess;
    m_position          = right.m_position;
    m_energyDeposited   = right.m_energyDeposited;
    m_kineticEnergy     = right.m_kineticEnergy;
    m_time              = right.m_time;

    return *this;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SDHit::operator==(const SDHit& right) const
{
    return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDHit::Print()
{
    G4cout << "---- SD Hit ----"
           << " trackId: " << m_trackId
           << " parentId: " << m_parentId
           << " particle: " << *m_particleType
           << " parentType: " << *m_parentType << G4endl
           << " CreatorProcess: " << *m_creatorProcess
           << " DepositingProcess: " << *m_depositingProcess << G4endl
           << " Position: " << m_position.x()/mm
           << " " << m_position.y()/mm
           << " " << m_position.z()/mm
           << " mm" << G4endl
           << " EnergyDeposited: " << m_energyDeposited/keV << " keV"
           << " KineticEnergyLeft: " << m_kineticEnergy/keV << " keV"
           << " Time: " << m_time/s << " s" << G4endl;
/*    

        << "EDeposit: "
        << std::setw(7) << G4BestUnit(m_eDeposit,"Energy")
        << ", Track Length: "
        << std::setw(7) << G4BestUnit(m_trackLength,"Length")
        << G4endl;
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
