#include "EventData.hh"

EventData::EventData()
{
    m_eventId = 0;
    m_nbHits = 0;
    m_nbSteps = 0;
    
    m_totalEnergyDeposited = 0.0;
    m_trackId = new vector<int>;
    m_parentId = new vector<int>;
    m_particleType = new vector<string>;
    m_parentType = new vector<string>;
    m_creatorProcess = new vector<string>;
    m_depositingProcess = new vector<string>;
    m_xp = new vector<double>;
    m_yp = new vector<double>;
    m_zp = new vector<double>;
    m_energyDeposited = new vector<double>;
    m_kineticEnergy = new vector<double>;
    m_time = new vector<double>;

    m_primaryParticleType = new vector<string>;
    m_primaryParticleEnergy = 0.0;
    m_primaryXm = 0.0;
    m_primaryYm = 0.0;
    m_primaryZm = 0.0;
    m_primaryX = 0.0;
    m_primaryY = 0.0;
    m_primaryZ = 0.0;
}

EventData::~EventData()
{
    delete m_trackId;
    delete m_parentId;
    delete m_particleType;
    delete m_parentType;
    delete m_creatorProcess;
    delete m_depositingProcess;
    delete m_xp;
    delete m_yp;
    delete m_zp;
    delete m_energyDeposited;
    delete m_kineticEnergy;
    delete m_time;
    delete m_primaryParticleType;
}

void
EventData::Clear()
{
    m_eventId = 0;
    m_nbHits = 0;
    m_nbSteps = 0;

    m_totalEnergyDeposited = 0.0;
    m_trackId->clear();
    m_parentId->clear();
    m_particleType->clear();
    m_parentType->clear();
    m_creatorProcess->clear();
    m_depositingProcess->clear();
    m_xp->clear();
    m_yp->clear();
    m_zp->clear();
    m_energyDeposited->clear();
    m_kineticEnergy->clear();
    m_time->clear();

    m_primaryParticleType->clear();
    m_primaryParticleEnergy = 0.0;
    m_primaryXm = 0.0;
    m_primaryYm = 0.0;
    m_primaryZm = 0.0;
    m_primaryX = 0.0;
    m_primaryY = 0.0;
    m_primaryZ = 0.0;
}
