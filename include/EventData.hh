#ifndef __EVENTDATA_H__
#define __EVENTDATA_H__

#include <string>
#include <vector>

using std::string;
using std::vector;

class EventData
{
public:

    EventData();
    ~EventData();
    
public:

    void Clear();
    
public:

    int m_eventId;
    int m_nbHits;
    int m_nbSteps;

    double m_totalEnergyDeposited;
    vector<int> *m_trackId;
    vector<int> *m_parentId;
    vector<string> *m_particleType;
    vector<string> *m_parentType;
    vector<string> *m_creatorProcess;
    vector<string> *m_depositingProcess;
    vector<double> *m_xp;               // position x,y,z
    vector<double> *m_yp;
    vector<double> *m_zp;
    vector<double> *m_energyDeposited;
    vector<double> *m_kineticEnergy;
    vector<double> *m_time;

    vector<string> *m_primaryParticleType;
    double m_primaryParticleEnergy;
    double m_primaryXm;
    double m_primaryYm;
    double m_primaryZm;
    double m_primaryX;
    double m_primaryY;
    double m_primaryZ;
};


#endif /* __EVENTDATA_H__ */
