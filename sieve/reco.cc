#include "my_belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
void User_reco::hist_def( void )
{ extern BelleTupleManager* BASF_Histogram;    
  t1 = BASF_Histogram->ntuple ("lmbda_lept",
    "ecm rm2l dsm dm dsp dp ntr chxc");
};

void User_reco::event ( BelleEvent* evptr, int* status ) {
    
  *status=0;
  
  static int nevent=0;
  static int nwritt=0;
  if(++nevent<2 || !(nevent%1000)) cout << "Event number " << nevent
          << " selected " << nwritt << endl;

  
  Mdst_sim_ecl_Manager& mng_secl = Mdst_sim_ecl_Manager::get_manager();
  Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
  Gen_hepevt_Manager& hepmgr = Gen_hepevt_Manager::get_manager();
  Belle_runhead_Manager::const_iterator rhd = rhdmgr.begin();
  Evtcls_hadron_info_Manager&  ehimgr =
    Evtcls_hadron_info_Manager::get_manager();

  Evtcls_hadronic_flag_Manager&  ehadfl =
    Evtcls_hadronic_flag_Manager::get_manager();

  HepPoint3D ip_position = IpProfile::position();
  const HepSymMatrix& runIp_err = IpProfile::position_err();
  Mdst_vee2_Manager &vee2_mgr = Mdst_vee2_Manager::get_manager();
  Evtcls_hadron_info_Manager::iterator iti = ehimgr.begin();
  Evtcls_hadronic_flag_Manager::iterator eti = ehadfl.begin();
  const HepPoint3D position = IpProfile::position();
  double r2=0;
  int ntrk=0;
  double evis=0;
  double Pz=0;
  double hjmass=0;

  if (iti!=ehimgr.end()){
    r2 = (*iti).R2();
    //  ntrk = (*iti).Ntrk();
    //      evis = (*iti).Evis();
    //      Pz   = (*iti).Pz();
    //      hjmass = (*iti).HeavyJetMass();
  }

  double ecm = BeamEnergy::Ecm();
  double elec = BeamEnergy::E_HER();
  double posi = BeamEnergy::E_LER();
  VectorL beam =VectorL(elec*sin(0.022), 0.,elec*cos(0.022)-posi, elec+posi);
  
  /*************** Make particle lists ********************************/

  //Base particles

  

  std::vector<Particle> p, ap, k_p, k_m, k0, pi_p, pi_m, pi0, gamma, gam, all, e_p, e_m, mu_m, mu_p, k_s, lam, alam;

  makeKPi(k_p, k_m, pi_p, pi_m, 1);

  //filter vectors k_s, k_p, k_m, pi_p, pi_m 
  withDrDzCut(pi_p, 1., 2.);
  withDrDzCut(pi_m, 1., 2.);

  withDrDzCut(k_m, 1., 2.);
  withDrDzCut(k_p, 1., 2.);
  withKaonIdCut(k_p, k_m, 0.6);

  //fill vectors k_s, k_p, k_m, pi_p, pi_m
  makeKs(k_s);

  if(k_m.size() + k_p.size() +  k_s.size() < 0.5) return;

  for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    if (abs(l->mass()-0.4977)>0.015 || V.perp()<0.1 ||
    cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1.) {
      k_s.erase(l); --l; continue;
    }
  }

  if(k_m.size() + k_p.size() +  k_s.size() < 0.5) return;
  
  //fill vectors lam, alam, p, ap
  makeProton(p, ap, 1);
  
  if(p.size() + ap.size()  < 0.5) return;

  //filter vectors p, ap
  withProtonIdCut(p, ap, 0.6);
  withDrDzCut(p, 1., 2.);
  withDrDzCut(ap, 1., 2.);

  if(p.size() + ap.size()  < 0.5) return;

  //fill vectors pi0 gamma
  makePi0(pi0);
  //makeK0(k0, ak0);
  withEminCutPi0(pi0, 0.05);

  //Undetected particles

  //D_pm

  std::vector<Particle> D_p, D_m;

  combination(D_p, m_ptypeDP, k_m, pi_p, pi_p, 0.05);
  combination(D_m, m_ptypeDM, k_p, pi_m, pi_m, 0.05);

  combination(D_p, m_ptypeDP, k_s, pi_p, 0.05);
  combination(D_m, m_ptypeDM, k_s, pi_m, 0.05);

  combination(D_p, m_ptypeDP, k_s, pi_p, pi_p, pi_m, 0.05);
  combination(D_m, m_ptypeDM, k_s, pi_m, pi_m, pi_p, 0.05);

  combination(D_p, m_ptypeDP, k_p, k_m, pi_p, 0.05);
  combination(D_m, m_ptypeDM, k_p, k_m, pi_m, 0.05);

  //D0

  std::vector<Particle> D0, aD0, D0_to_ds, aD0_to_ds;

  combination(D0, m_ptypeD0, k_m, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi_m, 0.05);

  combination(D0, m_ptypeD0, k_s, pi_p, pi_m, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi_p, pi_m, 0.05);

  combination(D0, m_ptypeD0, k_m, pi0, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi0, pi_m, 0.05);

  combination(D0, m_ptypeD0, k_s, pi_m, pi_p, pi0, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi_m, pi_p, pi0, 0.05);


  combination(D0_to_ds, m_ptypeD0, pi0, pi0, k_m, pi_p, 0.05);
  combination(aD0_to_ds, m_ptypeD0B, pi0, pi0, k_p, pi_m, 0.05);

  combination(D0, m_ptypeD0, k_m, k_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, k_m, 0.05);

  combination(D0, m_ptypeD0, k_s, pi0, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi0, 0.05);

  combination(D0, m_ptypeD0, k_m, k_p, 0.05);
  combination(aD0, m_ptypeD0B, k_m, k_m, 0.05);

  //D star

  std::vector<Particle> D0_st, aD0_st;

  combination(D0_st, m_ptypeDstar0, D0, pi0, 0.03);
  combination(aD0_st, m_ptypeDstarB, aD0, pi0, 0.03);

  std::vector<Particle> D_p_st, D_m_st;

  combination(D_p_st, m_ptypeDstarP, D0, pi_p, 0.03);
  combination(D_m_st, m_ptypeDstarM, aD0, pi_m, 0.03);

  combination(D_p_st, m_ptypeDstarP, D0_to_ds, pi_p, 0.03);
  combination(D_m_st, m_ptypeDstarM, aD0_to_ds, pi_m, 0.03);

  combination(D_p_st, m_ptypeDstarP, D_p, pi0, 0.03);
  combination(D_m_st, m_ptypeDstarM, D_m, pi0, 0.03);
  
  for(std::vector<Particle>::iterator D = D0_st.begin(); D!=D0_st.end(); ++D) {
    if (abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015 or D->mass() - D->child(0).mass() > 0.155) {
      D0_st.erase(D); --D;
    }
  }

  for(std::vector<Particle>::iterator D = aD0_st.begin(); D!=aD0_st.end(); ++D) {
    if (abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015 or D->mass() - D->child(0).mass() > 0.155){
      aD0_st.erase(D); --D;
    }
  }

  for(std::vector<Particle>::iterator D = D_p_st.begin(); D!=D_p_st.end(); ++D) {
    if (((abs(D->child(0).mass() - D->child(0).pType().mass())) > 0.015 or D->mass() - D->child(0).mass() > 0.155)) {
      D_p_st.erase(D); --D;
    }
  }

  for(std::vector<Particle>::iterator D = D_m_st.begin(); D!=D_m_st.end(); ++D) {
    if (((abs(D->child(0).mass() - D->child(0).pType().mass())) > 0.015 or D->mass() - D->child(0).mass() > 0.155)) {
      D_m_st.erase(D); --D;
    }
  }
  
  //X_c

  std::vector<Particle> X_c, aX_c;

  combination(X_c, m_ptypeUPS4, D0, p);
  combination(aX_c, m_ptypeUPS4, aD0, ap);
  setUserInfo(X_c, {{"chanel", 1}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c, {{"chanel", 1}, {"charg", -1},  {"baryon_num", -1}});

  combination(X_c, m_ptypeUPS4, D_p, p, pi_m);
  combination(aX_c, m_ptypeUPS4, D_m, ap, pi_p);
  setUserInfo(X_c, {{"chanel", 2}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c, {{"chanel", 2}, {"charg", -1},  {"baryon_num", -1}});

  combination(X_c, m_ptypeUPS4, D0_st, p);
  combination(aX_c, m_ptypeUPS4, aD0_st, ap);
  setUserInfo(X_c,  {{"chanel", 3}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c,  {{"chanel", 3}, {"charg", -1}, {"baryon_num", -1}});

  combination(X_c, m_ptypeUPS4, D_p_st, p, pi_m);
  combination(aX_c, m_ptypeUPS4, D_m_st, ap, pi_p);
  setUserInfo(X_c,  {{"chanel", 4}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c,  {{"chanel", 4}, {"charg", -1}, {"baryon_num", -1}});


for(int j=0; j<X_c.size(); ++j){
    Particle x_c=X_c[j];

    Particle ach = x_c.child(0);
    UserInfo chxc = static_cast<UserInfo&>(x_c.userInfo());
    
    if ((beam - (x_c.p())).m2() > 3 * 3) continue;
    int ntr=0;
    
    for(int jj=0; jj<pi_p.size(); ++jj)if (!checkSame(pi_p[jj], x_c)) ntr++;
    for(int jj=0; jj<pi_m.size(); ++jj)if (!checkSame(pi_m[jj], x_c)) ntr++;
    
    if((abs(ach.pType().lund()) == 423) or (abs(ach.pType().lund()) == 413)){
      t1->column("dsm", ach.p().m());
      t1->column("ch", abs(ach.pType().lund())-400);
      t1->column("dsp", ach.p().mag());
      t1->column("dm", ach.child(0).p().m());
      t1->column("dp", ach.child(0).p().mag());
    }
    else{
      t1->column("dm", ach.p().m());
      t1->column("dp", ach.p().mag());
    }

    t1->column("rm2l", (beam - (x_c.p())).m());    
    t1->column("ecm", ecm);    
    t1->column("ntr", ntr);
    t1->column("chxc", chxc.channel().find("chanel")->second);

    t1->dumpData();
    *status = 1; 

}

for(int j=0; j<aX_c.size(); ++j){
    Particle x_c=aX_c[j];

    Particle ach = x_c.child(0);
    UserInfo chxc = static_cast<UserInfo&>(x_c.userInfo());
    
    if ((beam - (x_c.p())).m2() > 3 * 3) continue;
    int ntr=0;
    
    for(int jj=0; jj<pi_p.size(); ++jj)if (!checkSame(pi_p[jj], x_c)) ntr++;
    for(int jj=0; jj<pi_m.size(); ++jj)if (!checkSame(pi_m[jj], x_c)) ntr++;
    
    if(chxc.channel().find("chanel")->second >= 3){
      t1->column("dsm", ach.p().m());
      t1->column("dsp", ach.p().mag());
      t1->column("dm", ach.child(0).p().m());
      t1->column("dp", ach.child(0).p().mag());
    }
    else{
      t1->column("dm", ach.p().m());
      t1->column("dp", ach.p().mag());
    }
    

    t1->column("rm2l", (beam - (x_c.p())).m());    
    t1->column("ecm", ecm);    
    t1->column("ntr", ntr);
    t1->column("chxc", chxc.channel().find("chanel")->second);

    t1->dumpData();
    *status = 1; 
}

if (*status==1) {nwritt++;
  cout << "Chac " << pi_p.size() << " " <<  pi_m.size() << " ntrack " << pi_p.size() + pi_m.size() << " ks_size " << k_s.size() << " lam size " << lam.size() + alam.size() <<endl;
}
}

#if defined(BELLE_NAMESPACE)
}
#endif
