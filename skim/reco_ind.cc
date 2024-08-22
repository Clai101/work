#include "my_belle.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
void User_reco::hist_def( void )
{ extern BelleTupleManager* BASF_Histogram;    
  t1 = BASF_Histogram->ntuple ("lmbda_lept",
    "ml mach p chxc chl pl pxc cos_lam_xc ch_dt_ach chach en dm_dst ecm ntr rm2n fox rm2l rm2nu chrgl chrgach en_gam enc_gam pl");
  /*
  t2 = BASF_Histogram->ntuple ("lmbdat",
    "en ecm p ntr chu chrgach chach mach rm2lc");
  */
};

void fitRM(const Particle &p1, const Particle &p2, const VectorL UPS, const double M_fit, VectorL &P_fit_1, VectorL &P_fit_2, double &chi2) 
  {  
    int ierr;
    double lambda=0;

    // cout << "*********************" << M_fit << endl;
    // first itteration : P_fit   <-- particle.p() 
    
    P_fit_1=p1.p();
    P_fit_2=p2.p();
    
    const double M_1=p1.pType().mass();
    const double M_2=p2.pType().mass();
    
    // fill momentum vector
    
    HepVector v_p1(3,0), v_p2(3,0);
    HepVector v_p10(3,0), v_p20(3,0);
    v_p1[0]=P_fit_1.x();
    v_p1[1]=P_fit_1.y();
    v_p1[2]=P_fit_1.z();
  
    v_p2[0]=P_fit_2.x();
    v_p2[1]=P_fit_2.y();
    v_p2[2]=P_fit_2.z();
  
    v_p10=v_p1;
    v_p20=v_p2;
  
    
    // fill momenutm error matrix
    
    HepSymMatrix err_p1(3,0), err_p2(3,0), unit(3,0);
    for(int im=0; im<3; ++im) {
      unit[im][im]=1;
      for(int in=0; in<3; ++in) { 
        err_p1[im][in]=p1.momentum().dp()[im][in];
        err_p2[im][in]=p2.momentum().dp()[im][in];
        if (abs(err_p1[im][in])<1.e-10) err_p1[im][in]=1.e-10;
        if (abs(err_p2[im][in])<1.e-10) err_p2[im][in]=1.e-10;
      }
    }
    
    // base a(0,1,2) 3-momentum of p1 ; a(4,5,6) - 3 momentum of p2 ; a(6) - lambda 
    HepSymMatrix A(7,0);
    HepVector a(7,0), diff_a(7,0);
    
    HepSymMatrix diff2_p1(3,0), diff2_p2(3,0);
    HepSymMatrix diff2_p1_p2(3,0), diff2_p2_p1(3,0);
    
    double diff;
    HepVector diff_p1(3,0), diff_p2(3,0);
    
    HepVector diff_chi_p1(3,0), diff_chi_p2(3,0);
    HepVector diff_m_p1(3,0), diff_m_p2(3,0);
    HepVector diff_lam_p1(3,0), diff_lam_p2(3,0);
    int n_fit=0;
    
    for(int itt=0; itt<100; ++itt) {
      
      n_fit++;
      double E_recoil=(UPS-P_fit_1-P_fit_2).e();
      double M_recoil=(UPS-P_fit_1-P_fit_2).m();

      double E_p1=P_fit_1.e();
      double E_p2=P_fit_2.e();

      diff=M_recoil* M_recoil - M_fit* M_fit;

      double diff_msdt= M_recoil - M_fit;
      if (abs(diff_msdt)<0.0001) continue;
      
      diff_chi_p1= 2. * err_p1.inverse(ierr) *  (v_p1-v_p10) ; 
      diff_m_p1= -2. * (( E_recoil / E_p1 + 1. ) * v_p1 + v_p2);
      
      diff_chi_p2= 2. * err_p2.inverse(ierr) *  (v_p2-v_p20) ; 
      diff_m_p2= -2. * (( E_recoil / E_p2 + 1. ) * v_p2 + v_p1);
      
      diff2_p1= 2. * err_p1.inverse(ierr) + 2. * lambda * ( E_recoil / E_p1 + 1. ) * unit;
    
      diff2_p2= 2. * err_p2.inverse(ierr) + 2. * lambda * ( E_recoil / E_p2 + 1. ) * unit;
    
      diff2_p1_p2 = + 2. * lambda * unit; 
      diff2_p2_p1 = + 2. * lambda * unit; 
      
      diff_lam_p1= -1. * diff_m_p1;
      diff_lam_p2= -1. * diff_m_p2;
      
      
      for(int im=0; im<3; ++im) {
        a[im]=diff_chi_p1[im]-lambda*diff_m_p1[im];
        a[im+3]=diff_chi_p2[im]-lambda*diff_m_p2[im];
        A[im][6]=diff_lam_p1[im];
        A[6][im]=diff_lam_p1[im];
        for(int in=0; in<3; ++in) { 
          A[im][in]=diff2_p1[im][in];
          A[im+3][in+3]=diff2_p2[im][in];
          A[im+3][in]=diff2_p1_p2[im][in];
          A[im][in+3]=diff2_p2_p1[im][in];
        }
      }
      A[6][6]=0;
      a[6]=-diff;
      
      // use scale factor 0.3 for provide convergency (slow)
      diff_a = -0.3 * A.inverse(ierr) * a;
      
      for(int im=0; im<3; ++im) {
	      diff_p1[im]=diff_a[im];
	      diff_p2[im]=diff_a[im+3];
      }
      
      lambda+=diff_a[6];
      v_p1+=diff_p1;
      v_p2+=diff_p2;
      
      chi2 = ((v_p1-v_p10).T() * err_p1.inverse(ierr) * (v_p1-v_p10))[0] + ((v_p2-v_p20).T() * err_p2.inverse(ierr) * (v_p2-v_p20))[0];

      //cout << itt << " diff=" << diff_msdt << " M_rec=" << M_recoil << " M_fit=" << M_fit << "; chi2=" << chi2 << endl;
      
      P_fit_1=VectorL(v_p1[0], v_p1[1], v_p1[2], 0);
      P_fit_2=VectorL(v_p2[0], v_p2[1], v_p2[2], 0);
      
      P_fit_1.setE(sqrt(P_fit_1.vect().mag2()+M_1*M_1));
      P_fit_2.setE(sqrt(P_fit_2.vect().mag2()+M_2*M_2));
      
    }

    // double M_recoil=(UPS-P_fit_1-P_fit_2).m();
    // cout<<"in fitrm: rmyf="<<M_recoil<< " n itter " << n_fit <<endl;

  };

int fill_tup(Particle lamc, /*vector<Particle> all,*/ double elec, double posi, double ecm, double r2, BelleTuple *t)
{    
  //int chb = dynamic_cast<UserInfo&>(B.userInfo()).channel();
  
  return 1;
};


void User_reco::event ( BelleEvent* evptr, int* status ) {
  




  *status=0;
  
  static int nevent=0;
  static int nwritt=0;
  if(++nevent<2 || !(nevent%1000)) cout << "Event number " << nevent
          << " selected " << nwritt << endl;
  
  Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
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

  //fill vectors
  makeLambda(lam,alam);
  makeProton(p, ap, 1);

  for(std::vector<Particle>::iterator l = lam.begin(); l!=lam.end(); ++l) {
      HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    double p_id;
    if (l->child(0).pType().mass()>l->child(1).pType().mass())
      p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(0).mdstCharged()));
    else p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(1).mdstCharged()));
    if (abs(l->mass()-1.1157)>0.01 || l->mdstVee2().z_dist()>1. || p_id<0.6 ) {
      lam.erase(l); --l;
    }
  }


  for(std::vector<Particle>::iterator l = alam.begin(); l!=alam.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    double p_id;
    if (l->child(0).pType().mass()>l->child(1).pType().mass()) 
    p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(0).mdstCharged()));
    else p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(1).mdstCharged()));                               
    if (abs(l->mass()-1.1157)>0.01 || l->mdstVee2().z_dist()>1. || p_id<0.6 ) {
      alam.erase(l); --l;
    }
  }
  withProtonIdCut(p, ap, 0.6);
  withDrDzCut(p, 1., 2.);
  withDrDzCut(ap, 1., 2.);

  makeKs(k_s);
  makeKPi(k_p, k_m, pi_p, pi_m, 1);


  //Cuts
  withDrDzCut(e_m, 1., 2.);
  withDrDzCut(e_p, 1., 2.);
  withDrDzCut(mu_m, 1., 2.);
  withDrDzCut(mu_p, 1., 2.);
  withDrDzCut(k_m, 1., 2.);
  withDrDzCut(k_p, 1., 2.);
  withDrDzCut(pi_p, 1., 2.);
  withDrDzCut(pi_m, 1., 2.);

  withKaonIdCut(k_p, k_m, 0.6);

  


  for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    if (abs(l->mass()-0.4977)>0.03 || V.perp()<0.1 ||
    cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1.) {
      k_s.erase(l); --l; continue;
    }
  }


  //static int n_trash_track = 0, n_trash_str = 0, n_trash_bar = 0, n_trash_str_bar = 0;
  int n_str = k_s.size() + k_m.size() + k_p.size() + lam.size() + alam.size();
  if (n_str < 1.5) return;
  if (pi_p.size() + pi_m.size() > 12.5) return;
  



  makeLepton(e_p, e_m, mu_p, mu_m, 1);
  withLeptonIdCut(e_p, e_m, mu_p, mu_m, 0.01, 0.1);
  
  makePi0(pi0);
  //makeK0(k0, ak0);
  withEminCutPi0(pi0, 0.05);
  //withEminCutK0(k0, 0.05);
  makeGamma(gamma);


  float f;

  setPi0Error(pi0);
  doKmvFit(pi0, f);

  doKmvFit(k_s, f);
  doKmvFit(lam, f);
  doKmvFit(alam, f);


  //Lambdac  

  std::vector<Particle> lamc_p, lamc_m;

  combination(lamc_p, m_ptypeLAMC, lam, e_p);
  combination(lamc_m, m_ptypeLAMC, alam, e_m);
  setUserInfo(lamc_p,  {{"chanel", 1}, {"charg", 1}, {"barion_num", -1}});
  setUserInfo(lamc_m,  {{"chanel", 1}, {"charg", -1}, {"barion_num", 1}});

  combination(lamc_p, m_ptypeLAMC, lam, mu_p);
  combination(lamc_m, m_ptypeLAMC, alam, mu_m);
  setUserInfo(lamc_p,  {{"chanel", 2}, {"charg", 1}, {"barion_num", -1}});
  setUserInfo(lamc_m,  {{"chanel", 2}, {"charg", -1}, {"barion_num", 1}});

  combination(lamc_m, m_ptypeLAMC, alam, pi_m, 0.05);
  combination(lamc_p, m_ptypeLAMC, lam, pi_p, 0.05);
  setUserInfo(lamc_m,  {{"chanel", 3}, {"charg", -1}, {"barion_num", 1}});
  setUserInfo(lamc_p,  {{"chanel", 3}, {"charg", 1}, {"barion_num", -1}});

  combination(lamc_m, m_ptypeLAMC, ap, k_p, pi_m, 0.05);
  combination(lamc_p, m_ptypeLAMC, p, k_m, pi_p, 0.05);
  setUserInfo(lamc_m, {{"chanel", 5}, {"charg", -1}, {"barion_num", -1}});
  setUserInfo(lamc_p, {{"chanel", 5}, {"charg", 1}, {"barion_num", 1}});
  
  if (lamc_p.size()+lamc_m.size()==0) return;

  doKvFit(lamc_m);
  doKvFit(lamc_p);

  doKmvFit(lamc_m, f);
  doKmvFit(lamc_p, f);

  //Rho

  std::vector<Particle> ups, rho_2m, rho_2p, rho_ppm, rho_mmp, rho;
  
  combination(rho_mmp, m_ptypeRHO0, pi_m, pi_m, pi_p);
  setUserInfo(rho_mmp, {{"chanel", 1}, {"charg", -1}, {"barion_num", 0}});
  combination(rho_ppm, m_ptypeRHO0, pi_p, pi_p, pi_m);
  setUserInfo(rho_ppm, {{"chanel", 1}, {"charg", 1}, {"barion_num", 0}});

  doKvFit(rho_mmp);
  doKvFit(rho_ppm);

  doKmvFit(rho_mmp, f);
  doKmvFit(rho_ppm, f);

  //D_pm

  std::vector<Particle> D_p, D_m;

  combination(D_p, m_ptypeDP, k_m, pi_p, pi_p, 0.05);
  setUserInfo(D_p,  {{"chanel", 1}, {"charg", 1}, {"barion_num", 0}});
  combination(D_m, m_ptypeDM, k_p, pi_m, pi_m, 0.05);
  setUserInfo(D_m,  {{"chanel", 1}, {"charg", -1}, {"barion_num", 0}});

  combination(D_p, m_ptypeDP, k_s, pi_p, 0.05);
  setUserInfo(D_p, {{"chanel", 3}, {"charg", 1}, {"barion_num", 0}});
  combination(D_m, m_ptypeDM, k_s, pi_m, 0.05);
  setUserInfo(D_m, {{"chanel", 3}, {"charg", -1}, {"barion_num", 0}});

  combination(D_p, m_ptypeDP, k_s, pi_p, pi_p, pi_m, 0.05);
  setUserInfo(D_p, {{"chanel", 5}, {"charg", 1}, {"barion_num", 0}});
  combination(D_m, m_ptypeDM, k_s, pi_m, pi_m, pi_p, 0.05);
  setUserInfo(D_m, {{"chanel", 5}, {"charg", -1}, {"barion_num", 0}});

  combination(D_p, m_ptypeDP, k_p, k_m, pi_p, 0.05);
  setUserInfo(D_p, {{"chanel", 6}, {"charg", 1}, {"barion_num", 0}});
  combination(D_m, m_ptypeDM, k_p, k_m, pi_m, 0.05);
  setUserInfo(D_m, {{"chanel", 6}, {"charg", -1}, {"barion_num", 0}});

  doKvFit(D_p);
  doKvFit(D_m);

  doKmvFit(D_p, f);
  doKmvFit(D_m, f);

  //D0

  std::vector<Particle> D0, aD0;

  combination(D0, m_ptypeD0, k_m, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi_m, 0.05);
  setUserInfo(D0, {{"chanel", 1}, {"charg", 0}, {"barion_num", 0}});
  setUserInfo(aD0, {{"chanel", 1}, {"charg", 0}, {"barion_num", 0}});

  combination(D0, m_ptypeD0, k_m, pi_p, pi_p, pi_m, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi_m, pi_m, pi_p, 0.05);
  setUserInfo(D0,  {{"chanel", 2}, {"charg", 0}, {"barion_num", 0}});
  setUserInfo(aD0,  {{"chanel", 2}, {"charg", 0}, {"barion_num", 0}});

  combination(D0, m_ptypeD0, k_m, pi0, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi0, pi_m, 0.05);
  setUserInfo(D0, {{"chanel", 3}, {"charg", 0}, {"barion_num", 0}});
  setUserInfo(aD0, {{"chanel", 3}, {"charg", 0}, {"barion_num", 0}});

  combination(D0, m_ptypeD0, k_m, k_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, k_m, 0.05);
  setUserInfo(D0, {{"chanel", 4}, {"charg", 0}, {"barion_num", 0}});
  setUserInfo(aD0, {{"chanel", 4}, {"charg", 0}, {"barion_num", 0}});

  combination(D0, m_ptypeD0, k_s, pi_p, pi_m, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi_p, pi_m, 0.05);
  setUserInfo(D0, {{"chanel", 5}, {"charg", 0}, {"barion_num", 0}});
  setUserInfo(aD0, {{"chanel", 5}, {"charg", 0}, {"barion_num", 0}});

  combination(D0, m_ptypeD0, k_s, pi0, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi0, 0.05);
  setUserInfo(D0, {{"chanel", 6}, {"charg", 0}, {"barion_num", 0}});
  setUserInfo(aD0, {{"chanel", 6}, {"charg", 0}, {"barion_num", 0}});

  doKvFit(D0);
  doKvFit(aD0);

  doKmvFit(D0, f);
  doKmvFit(aD0, f);

  //D star

  std::vector<Particle> D0_st, aD0_st;

  combination(D0_st, m_ptypeDstar0, D0, pi0, 0.02);
  combination(aD0_st, m_ptypeDstarB, aD0, pi0, 0.02);
  setUserInfo(D0_st, {{"chanel", 1}, {"charg", 0}, {"barion_num", 0}});
  setUserInfo(aD0_st, {{"chanel", 1}, {"charg", 0}, {"barion_num", 0}});

  //
  std::vector<Particle> D_p_st, D_m_st;

  combination(D_p_st, m_ptypeDstar0, D0, pi_p, 0.02);
  combination(D_m_st, m_ptypeDstarB, aD0, pi_m, 0.02);
  setUserInfo(D_p_st, {{"chanel", 1}, {"charg", 1}, {"barion_num", 0}});
  setUserInfo(D_m_st, {{"chanel", 1}, {"charg", -1}, {"barion_num", 0}});

  combination(D_p_st, m_ptypeDstar0, D_p, pi0, 0.02);
  combination(D_m_st, m_ptypeDstarB, D_m, pi0, 0.02);
  setUserInfo(D_p_st, {{"chanel", 2}, {"charg", 1}, {"barion_num", 0}});
  setUserInfo(D_m_st, {{"chanel", 2}, {"charg", -1}, {"barion_num", 0}});

  for(std::vector<Particle>::iterator D = D0_st.begin(); D!=D0_st.end(); ++D) {
    if ((abs(dynamic_cast<UserInfo&>(D->child(0).userInfo()).vmass() - D->child(0).pType().mass())) > 0.015) {
      D0_st.erase(D); --D; continue;
    }
  }

  for(std::vector<Particle>::iterator D = D_p_st.begin(); D!=D_p_st.end(); ++D) {
    if (((abs(dynamic_cast<UserInfo&>(D->child(0).userInfo()).vmass() - D->child(0).pType().mass())) > 0.015)) {
      D_p_st.erase(D); --D; continue;
    }
  }

  for(std::vector<Particle>::iterator D = D_m_st.begin(); D!=D_m_st.end(); ++D) {
    if (((abs(dynamic_cast<UserInfo&>(D->child(0).userInfo()).vmass() - D->child(0).pType().mass())) > 0.015)) {
      D_m_st.erase(D); --D; continue;
    }
  }

  for(std::vector<Particle>::iterator D = aD0_st.begin(); D!=aD0_st.end(); ++D) {
    if (((abs(dynamic_cast<UserInfo&>(D->child(0).userInfo()).vmass() - D->child(0).pType().mass())) > 0.015)) {
      aD0_st.erase(D); --D; continue;
    }
  }

  doKvFit(D_p_st, f);
  doKvFit(D_m_st,f);

  doKvFit(aD0_st, f);
  doKvFit(aD0_st,f);


  doKmvFit(D_p_st, f);
  doKmvFit(D_m_st,f);

  doKmvFit(aD0_st, f);
  doKmvFit(aD0_st,f);

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

  combination(X_c, m_ptypeUPS4, D_p_st, p, pi_m);
  combination(aX_c, m_ptypeUPS4, D_m_st, ap, pi_p);
  setUserInfo(X_c,  {{"chanel", 3}, {"charg", 1}, {"barion_num", 1}});
  setUserInfo(aX_c,  {{"chanel", 3}, {"charg", -1}, {"barion_num", -1}});

  combination(X_c, m_ptypeUPS4, D0_st, p);
  combination(aX_c, m_ptypeUPS4, D0_st, ap);
  setUserInfo(X_c,  {{"chanel", 4}, {"charg", 1}, {"barion_num", 1}});
  setUserInfo(aX_c,  {{"chanel", 4}, {"charg", -1}, {"barion_num", -1}});


  combination(ups, m_ptypeUPS4, lamc_m, X_c);
  combination(ups, m_ptypeUPS4, lamc_p, aX_c);




  for(int j=0; j<ups.size(); ++j){
    Particle u=ups[j];

    short ntr=0;
    for(int jj=0; jj<pi_p.size(); ++jj) 
    if (!checkSame(pi_p[jj],u)) ntr++;

    for(int jj=0; jj<pi_m.size(); ++jj) 
    if (!checkSame(pi_m[jj],u)) ntr++;



    float en_gam=0, enc_gam=0;
    for(int jj=0; jj<gamma.size(); ++jj){
    if (!checkSame(gamma[jj],u)){
      enc_gam = enc_gam + pStar(gamma[jj], elec, posi).e();
      en_gam = en_gam + (gamma[jj]).p().e();
    }} 

    Particle lam = u.child(0);
    Particle x_c = u.child(1);
    Particle ach = x_c.child(0);
    Particle dt_ach = ach.child(0);

    map <string, int> chl = dynamic_cast<UserInfo&>(lam.userInfo());
    map <string, int> chach = dynamic_cast<UserInfo&>(ach.userInfo());
    map <string, int> ch_dt_ach = dynamic_cast<UserInfo&>(dt_ach.userInfo());
    map <string, int> chxc = dynamic_cast<UserInfo&>(x_c.userInfo());
    int chargU = 0;

    //VectorL new_lam_p, new_xc_p;
    //double chi;

    //fitRM(&lam, &x_c, u.p(), 2.28646, &new_lam_p, &new_xc_p, &chi); 


    if ((beam - (x_c.p())).m2() > 4*4 or (ntr >= 1)) continue;    

    float delta_dst = -400;
    if (chxc.channel()["chanel"] >= 3) delta_dst =  ach.mass() - ach.child(0).mass();


    t1->column("en", pStar(u, elec, posi).e());
    t1->column("ecm", ecm);
    t1->column("p", pStar(u, elec, posi).vect().mag());
    t1->column("ntr", ntr);

    t1->column("en_gam", en_gam);
    t1->column("enc_gam", enc_gam);

    t1->column("chxc", chxc.channel()["chanel"]);     
    t1->column("pxc", pStar(x_c, elec, posi).vect().mag());     

    t1->column("chach", chach.channel()["chanel"]);
    t1->column("chrgach", chach.channel()["charg"]);
    t1->column("mach", chach.vmass());

    t1->column("m_dt_ach", ch_dt_ach.vmass());
    t1->column("ch_dt_ach", ch_dt_ach.channel()["chanel"]);

    t1->column("chrgl", chl["charg"]);
    t1->column("pl", pStar(lam, elec, posi).vect().mag());     
    t1->column("chl", chl["chanel"]);
    t1->column("ml", chl.vmass());


    t1->column("cos_lam_xc", pStar(lam, elec, posi).vect() * pStar(x_c, elec, posi).vect()/(pStar(x_c, elec, posi).vect().mag()*pStar(lam, elec, posi).vect().mag()));
    t1->column("rm2n", (beam - u.p()).m2());
    t1->column("rm2l", (beam - (x_c.p())).m2());
    t1->column("rm2nu", (beam - (x_c.p() + lam.child(0).p() + lam.child(1).p())).m2());
    t1->column("new_rm2nu", (beam - (x_c.p() + lam.child(0).p() + lam.child(1).p())).m2());

    t1->column("fox", r2);




    if(chl["chanel"] <= 2){
      t1->column("pl", pStar(lam.child(1), elec, posi).vect().mag());}
    else{
      t1->column("pl", 0);}
    

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