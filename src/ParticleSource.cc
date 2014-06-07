#include <G4PrimaryParticle.hh>
#include <G4Event.hh>
#include <G4TransportationManager.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4IonTable.hh>
#include <G4Ions.hh>
#include <G4TrackingManager.hh>
#include <G4Track.hh>
#include "G4SystemOfUnits.hh"
#include <Randomize.hh>

#include <G4Geantino.hh>
#include <G4Electron.hh>
#include <G4Positron.hh>

#include "ParticleSource.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParticleSource::ParticleSource()
    : m_EHistN(0), m_EHistE(0), m_EHistP(0), m_EGen(0)
{
    m_sourceMessenger = new ParticleSourceMessenger(this);

    m_gunType = "";
    m_beNu_a = -1.0/3.0;
    m_QValue = 0.0;
    m_fIonPos = G4ThreeVector(0.0,0.0,0.0);
    m_fIonEk  = 0.0;
    m_fIonPDir = G4ThreeVector(1.0,0.0,0.0);
    m_fIonDef = 0; m_dIonDef = 0;
    m_fIonCharge = 0.0; m_dIonCharge = 0.0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParticleSource::~ParticleSource()
{
    delete m_sourceMessenger;
    if(m_EHistE) delete[] m_EHistE;
    if(m_EHistP) delete[] m_EHistP;
    if(m_EGen) delete m_EGen;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int ParticleSource::LoadEHist(G4String EHistFName)
{
    FILE *fp;
    char buf[LINE_MAX];
    size_t i;
    G4double he, hp;

    if(m_EHistE && m_EHistP) { // file has already been loaded
        G4cout << EHistFName << " has been loaded already." << G4endl;
        return 2;
    }
    
    if((fp=fopen(EHistFName.data(), "r")) == NULL) {
        perror(EHistFName.data());
        return 0;
    }
    
    m_EHistN = 0;
    while(fgets(buf, sizeof(buf), fp) && (!feof(fp))) {
        if(buf[0]=='#') continue;
        m_EHistN++;
    }
    rewind(fp);
    m_EHistE = new G4double[m_EHistN];
    m_EHistP = new G4double[m_EHistN];
    i = 0;
    while(fgets(buf, sizeof(buf), fp) && (!feof(fp))) {
        if(buf[0]=='#') continue;
        sscanf(buf, "%lf%lf", &he, &hp);
        m_EHistE[i] = he * MeV;
        m_EHistP[i] = hp;
        i++;
    }
    fclose(fp);
    m_EGen = new G4RandGeneral(m_EHistP, m_EHistN);
    return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// generate random variate with linear pdf y = a * x + 1/2 with x in [-1,1]
static inline G4double inv_ld_cdf(G4double u, G4double a)
{
    G4double p0;
    p0 = 1.0 - a - 2.0 * u;
    return (-2.0 * p0) / (sqrt(1.0 - 4.0 * a * p0) + 1.0);
}

static inline G4double rand_ld(G4double a)
{
    return inv_ld_cdf(G4UniformRand(), a);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// solve for p_nu knowing p_e and Q, assuming m_nu = 0
static G4double solve_p_nu(G4double pe, G4double Q, G4double me, G4double cosTh, G4double md)
{
    G4double tmp, sol;
    
    tmp = 2.0*md + 2.0*cosTh*pe;
    sol = 1.0/2.0*(-2.0*md-2.0*cosTh*pe
                   + std::sqrt(tmp*tmp-4.0*(-2.0*md*me + pe*pe
                                            + 2.0*md*std::sqrt(me*me + pe*pe) - 2.0*md*Q)));
    return sol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticleSource::GeneratePrimaryVertex(G4Event *event)
{
    G4double primaryParticleTime = 0.0;
    G4double tmp, x1, x2;
    G4double energy, mass, charge, boostBeta;
    G4double betaEk, betaP, betaMass, cosTh;
    G4double nuP, dIonMass;
    G4ThreeVector eP3v, eP3vNorm, nuP3v, dIP3v, rotAxisV;
    G4LorentzVector lv;
    G4PrimaryParticle *particle;

    if(m_EGen == 0) {
        G4cerr << "Load EHist first!\n" << G4endl;
        return;
    }

    mass = m_fIonDef->GetPDGMass();
    tmp = mass/(m_fIonEk + mass);
    boostBeta = std::sqrt(1.0 - tmp * tmp); // v/c value for boost

    // create a new vertex
    G4PrimaryVertex *vertex = new G4PrimaryVertex(m_fIonPos, primaryParticleTime);

    // generate beta
    charge = G4Electron::ElectronDefinition()->GetPDGCharge();
    mass = betaMass = G4Electron::ElectronDefinition()->GetPDGMass();
    betaEk = m_EGen->shoot() * m_QValue;
    energy = betaEk + mass;
    betaP  = std::sqrt(energy * energy - mass * mass);
    // uniform sphere point picking
    do {
        x1 = G4UniformRand()*2.0 - 1.0;
        x2 = G4UniformRand()*2.0 - 1.0;
        tmp = x1*x1 + x2*x2;
    } while(tmp>=1.0);
    eP3v.setX(2.0*x1*std::sqrt(1.0 - tmp));
    eP3v.setY(2.0*x2*std::sqrt(1.0 - tmp));
    eP3v.setZ(1.0 - 2.0 * tmp);

/*  // does the same thing but less efficient
    x1 = G4UniformRand()*2.0 - 1.0;
    x2 = std::sqrt(1.0 - x1*x1);
    tmp = 2.0*M_PI*G4UniformRand();
    eP3v.setX(x2 * cos(tmp));
    eP3v.setY(x2 * sin(tmp));
    eP3v.setZ(x1);
*/

    eP3vNorm = eP3v/eP3v.mag();
    eP3v *= betaP/eP3v.mag();
    lv.setVect(eP3v);
    lv.setE(energy);
    lv.boost(m_fIonPDir, boostBeta);

    particle = new G4PrimaryParticle(G4Electron::ElectronDefinition(),
                                     lv.px(), lv.py(), lv.pz());
    particle->SetMass(mass);
    particle->SetCharge(charge);
    // particle->SetPolarization(,,);
    vertex->SetPrimary(particle);

    // generate recoil ion
    charge = m_dIonCharge;
    mass = dIonMass = m_dIonDef->GetPDGMass();
    cosTh = rand_ld(m_beNu_a);
    nuP = solve_p_nu(betaP, m_QValue, betaMass, cosTh, dIonMass);
    rotAxisV = eP3vNorm.cross(G4ThreeVector(0,0,1));
    // if eP3vNorm happens to be z axis:
    if(rotAxisV.mag2() == 0.0) rotAxisV.set(1.0, 0.0, 0.0);
    nuP3v = eP3vNorm;
    nuP3v.rotate(rotAxisV, std::acos(cosTh));
    nuP3v.rotate(eP3vNorm, 2.0*M_PI*G4UniformRand());
    nuP3v *= nuP/nuP3v.mag();
    m_nuP = nuP3v;
    m_nuE = nuP3v.mag();

    lv.setVect(m_nuP);
    lv.setE(m_nuE);
    lv.boost(m_fIonPDir, boostBeta);
    particle = new G4PrimaryParticle(G4Geantino::GeantinoDefinition(),
                                     lv.px(), lv.py(), lv.pz());
    vertex->SetPrimary(particle);

    dIP3v = - (eP3v + nuP3v);
    energy = std::sqrt(dIP3v.mag2() + dIonMass * dIonMass);
    lv.setVect(dIP3v);
    lv.setE(energy);
    lv.boost(m_fIonPDir, boostBeta);
    particle = new G4PrimaryParticle(m_dIonDef, lv.px(), lv.py(), lv.pz());

    particle->SetMass(mass);
    particle->SetCharge(charge);
    vertex->SetPrimary(particle);
/*
    G4cout << "fIonMass: " << m_fIonDef->GetPDGMass() << "boostBeta: " << boostBeta
           << " dIonMass: " << dIonMass << " charge: " << charge
           << " dIP3v.mag(): " << dIP3v.mag() << " energy: " << energy << " lv mag: "
           << std::sqrt(lv.px()*lv.px() + lv.py()*lv.py() + lv.pz()*lv.pz())
           << G4endl;
*/
    event->AddPrimaryVertex(vertex);
}
