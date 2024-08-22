#include "userinfo.h"
#include <vector>
#include <map>
#include <string> 
// For Interface to Set UserInfo Class

using namespace Belle;

namespace Belle{
void setUserInfo(Particle &p, map <string, int> ch)
{ 
  if(! &p.userInfo() )// return;
    p.userInfo((UserInfo(ch))); 
}

void setUserInfo(vector<Particle> &p, map <string, int> ch)
{ for(int i=0;i<p.size();++i) setUserInfo(p[i],ch); }

// UserInfo Class

UserInfo::UserInfo()
  : m_vchisq(-1.),
    m_bchisq(-1.),
    m_mchisq(-1.),
    m_id(0.),
    m_K_angl(0.),
    m_Pi_angl(0.),
    m_vmass(0.),
    m_mass(0.),
    m_vx(0.),  
    m_vy(0.),  
    m_vz(0.),  
    m_channel({{"chanel", 0}})
{
}

UserInfo::UserInfo(map <string, int> ch)
  : m_vchisq(-1.),
    m_bchisq(-1.),
    m_mchisq(-1.),
    m_id(0.),
    m_K_angl(0.),
    m_Pi_angl(0.),
    m_vmass(0.),
  
    m_mass(0.),
    m_vx(0.),
    m_vy(0.),
    m_vz(0.), 
    m_channel(ch)
{
}

UserInfo::UserInfo(const UserInfo &x)
  : m_vchisq(x.m_vchisq),
    m_bchisq(x.m_bchisq),
    m_mchisq(x.m_mchisq),
    m_id(x.m_id), 
    m_K_angl(x.m_K_angl),
    m_Pi_angl(x.m_Pi_angl),
    m_vmass(x.m_vmass),
   
    m_mass(x.m_mass),
    m_vx(x.m_vx),
    m_vy(x.m_vy),  
    m_vz(x.m_vz), 
    m_channel(x.m_channel)
{
}
UserInfo::~UserInfo()
{
}

UserInfo*
UserInfo::clone(void) const
{
  UserInfo *x = new UserInfo(*this);
  return x;
}

UserInfo &
UserInfo::operator = (const UserInfo &x)
{
  m_vchisq = x.m_vchisq;
  m_bchisq = x.m_bchisq;
  m_mchisq = x.m_mchisq;
  m_id = x.m_id;
  m_K_angl = x.m_K_angl;
  m_Pi_angl = x.m_Pi_angl;
  m_vmass  = x.m_vmass; 

  m_mass  = x.m_mass;
  m_vx  = x.m_vx; 
  m_vy  = x.m_vy; 
  m_vz  = x.m_vz;
  m_channel= x.m_channel;
}









}





