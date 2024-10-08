#include "belle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include <iostream>
#include MDST_H
#include BELLETDF_H
#include HEPEVT_H
#include "particle/utility.h"
#include "particle/combination.h"
#include "kid/atc_pid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "eid/eid.h"
#include "ip/IpProfile.h"

using namespace std;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class User_reco : public Module {
public:
  User_reco ( void );
  ~User_reco ( void ){};
  void init ( int* ){};
  void term ( void ){};
  void disp_stat ( const char* ){};
  void hist_def ( void );
  void event ( BelleEvent*, int* ); 
  void begin_run ( BelleEvent*, int* );
  void end_run   ( BelleEvent*, int* ){};
  void other ( int*, BelleEvent*, int* ){};

  /*******************************************
    Define things here that you want to use 
    inside and outside of the event function
   *******************************************/
  BelleTuple *t1;
  BelleTuple *t2;
  BelleTuple *t3;
  BelleTuple *t4;
  BelleTuple *t5;
  BelleTuple *t6;
  BelleTuple *t7;
  BelleTuple *t8;
  BelleTuple *t9;
  BelleTuple *t10;
  BelleTuple *t11;
  BelleTuple *t12;
  BelleTuple *t13;
  BelleTuple *t14;
  BelleTuple *t15;
  BelleTuple *t16;
  BelleTuple *t17;
  BelleTuple *t18;
  BelleTuple *t19;

private:
  Ptype m_ptypeGAMM;
  Ptype m_ptypeRHO0;
  Ptype m_ptypePI0;
  Ptype m_ptypePHI;
  Ptype m_ptypeETAC;
  Ptype m_ptypeCHI0;
  Ptype m_ptypePSI;
  Ptype m_ptypePSI2;
  Ptype m_ptypeLAMC;
  Ptype m_ptypeALAMC;
  Ptype m_ptypeSIGC0;
  Ptype m_ptypeDSP;
  Ptype m_ptypeDSM;
  Ptype m_ptypeDSstarP;
  Ptype m_ptypeDSstarM;
  Ptype m_ptypeDP;
  Ptype m_ptypeDM;
  Ptype m_ptypeD0;
  Ptype m_ptypeD0B;
  Ptype m_ptypeB0;
  Ptype m_ptypeB0B;
  Ptype m_ptypeBP;
  Ptype m_ptypeBM;
  Ptype m_ptypeDstarP;
  Ptype m_ptypeDstarM;
  Ptype m_ptypeKstar0;
  Ptype m_ptypeUPS4;
  Ptype m_ptypeDstar0;
  Ptype m_ptypeDstarB;
};

extern "C" Module_descr *mdcl_User_reco ()
{ /* main */
  User_reco *module = new User_reco;
  Module_descr *dscr = new Module_descr ( "User_reco", module );
  return dscr;
};

/* Constructor  */
User_reco::User_reco ( void ) 
  :  m_ptypeDSP("DS+"),
     m_ptypeDSM("DS-"),
     m_ptypeDP("D+"),
     m_ptypeDM("D-"),
     m_ptypeD0("D0"),
     m_ptypeD0B("D0B"),
     m_ptypeLAMC("LAMC"),
     m_ptypeALAMC("ALAMC"),
     m_ptypeSIGC0("SIGC0"),
     m_ptypeKstar0("K*0"),
     m_ptypeDstarP("D*+"),
     m_ptypeDstarM("D*-"),
     m_ptypeB0("B0"),
     m_ptypeB0B("B0B"),
     m_ptypeBP("B+"),
     m_ptypeBM("B-"),
     m_ptypePSI("PSI"),
     m_ptypePSI2("PSI2"),
     m_ptypePHI("PHI"),
     m_ptypeRHO0("RHO0"),
     m_ptypePI0("PI0"),
     m_ptypeETAC("ETAC"),
     m_ptypeCHI0("CHI0"),
     m_ptypeGAMM("GAMM"),
     m_ptypeUPS4("UPS4"),
     m_ptypeDstar0("D*0"),
     m_ptypeDstarB("D*B")

{
  /*******************************************
   Define things here that you want to use 
   inside and outside of the event function
   ******************************************/
};

void User_reco::begin_run ( BelleEvent*, int* )
{
  eid::init_data();
  IpProfile::begin_run();
  IpProfile::dump();
}

#if defined(BELLE_NAMESPACE)
}
#endif
