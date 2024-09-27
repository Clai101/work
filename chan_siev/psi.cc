#include "belle.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include MDST_H
#include BELLETDF_H
#include HEPEVT_H
#include "particle/utility.h"
#include "particle/combination.h"
#include "kid/atc_pid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "eid/eid.h"
#include "psi.h"
#include "userinfo.h"
#include "myutils.h"
#include "ip/IpProfile.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

Ptype m_ptypeGAMM("GAMM");
Ptype m_ptypeRHO0("RHO0");
Ptype m_ptypePI0("Pi0");
Ptype m_ptypePHI("PHI");
Ptype m_ptypeETAC("ETAC");
Ptype m_ptypeCHI0("CHI0");
Ptype m_ptypePSI("PSI");
Ptype m_ptypePSI2("PSI2");
Ptype m_ptypeDP("D+");
Ptype m_ptypeDM("D-");
Ptype m_ptypeD0("D0");
Ptype m_ptypeD0B("D0B");
Ptype m_ptypeB0("B0");
Ptype m_ptypeB0B("B0B");
Ptype m_ptypeBP("B+");
Ptype m_ptypeBM("B-");
Ptype m_ptypeDstarP("D*+");
Ptype m_ptypeDstarM("D*-");
Ptype m_ptypeDSP("DS+");
Ptype m_ptypeDSM("DS-");
Ptype m_ptypeDSstarP("DS*+");
Ptype m_ptypeDSstarM("DS*-");
Ptype m_ptypeLAMC("LAMC");
Ptype m_ptypeALAMC("ALAMC");
Ptype m_ptypeUPS4("UPS4");

//using namespace Belle;



void
withDrDzCut(std::vector<Particle> &list1, double dr, double dz){
  for(std::vector<Particle>::iterator l = list1.begin(); l!=list1.end(); ++l){
    int mhyp=2;
    if (abs(l->lund())==2212) mhyp=4;
    if (abs(l->lund())==321) mhyp=3;
    if (abs(l->lund())==11) mhyp=0;
    if (abs(l->lund())==12) mhyp=1;
    if (l->mdstCharged() && 
	abs(l->mdstCharged().trk().mhyp(mhyp).helix(0))>dr ||
	abs(l->mdstCharged().trk().mhyp(mhyp).helix(3)
	    -IpProfile::position().z())>dz )
      {list1.erase(l); --l; }
  }
}

void
withEminCut(std::vector<Particle> &gam, 
		double Emin) {
  for(std::vector<Particle>::iterator l = gam.begin(); l!=gam.end(); ++l) 
    if (l->p().e()<Emin)
      {gam.erase(l); --l; }
  
}

void
withEminCutPi0(std::vector<Particle> &pi0, 
		double Emin) {
  for(std::vector<Particle>::iterator l = pi0.begin(); l!=pi0.end(); ++l)
    if ((l->mdstPi0().gamma(0).ecl().energy()<Emin)||
	(l->mdstPi0().gamma(1).ecl().energy()<Emin))
      {pi0.erase(l); --l; }
}

void
withPStarMinCut(std::vector<Particle> &in, 
		double Psmin) {

  for(std::vector<Particle>::iterator l = in.begin(); l!=in.end(); ++l) 
    if (pStar(l->p()).vect().mag()<Psmin)
      {in.erase(l); --l; }
}

void
withProtonIdCut(std::vector<Particle> &pr, 
		std::vector<Particle> &apr, 
		double lhpcut) {
  for(std::vector<Particle>::iterator l = pr.begin(); l!=pr.end(); ++l) 
    if (atc_pid(3,1,5,4,3).prob(&(l->mdstCharged()))<lhpcut ||
	atc_pid(3,1,5,4,2).prob(&(l->mdstCharged()))<lhpcut)
      {pr.erase(l); --l; }
  
  for(std::vector<Particle>::iterator l = apr.begin(); l!=apr.end(); ++l) 
    if (atc_pid(3,1,5,4,3).prob(&(l->mdstCharged()))<lhpcut ||
	atc_pid(3,1,5,4,2).prob(&(l->mdstCharged()))<lhpcut)
      {apr.erase(l); --l; }
}

void
withKaonIdCut(std::vector<Particle> &k_p, 
	      std::vector<Particle> &k_m, 
	      double lhkcut) {

  for(std::vector<Particle>::iterator l = k_p.begin(); l!=k_p.end(); ++l) 
    if (atc_pid(3,1,5,3,2).prob(&(l->mdstCharged()))<lhkcut)
      {k_p.erase(l); --l; }
  
  for(std::vector<Particle>::iterator l = k_m.begin(); l!=k_m.end(); ++l) 
    if (atc_pid(3,1,5,3,2).prob(&(l->mdstCharged()))<lhkcut)
      {k_m.erase(l); --l; }
}

void
withPionIdCut(std::vector<Particle> &pi_p, 
	      std::vector<Particle> &pi_m, 
	      double lhpcut) {

  for(std::vector<Particle>::iterator l = pi_p.begin(); l!=pi_p.end(); ++l) 
    if (atc_pid(3,1,5,2,3).prob(&(l->mdstCharged()))<lhpcut)
      {pi_p.erase(l); --l; }
  
  for(std::vector<Particle>::iterator l = pi_m.begin(); l!=pi_m.end(); ++l) 
    if (atc_pid(3,1,5,2,3).prob(&(l->mdstCharged()))<lhpcut)
      {pi_m.erase(l); --l; }
}

void
withLeptonIdCut(std::vector<Particle> &e_p, 
		std::vector<Particle> &e_m, 
		std::vector<Particle> &mu_p, 
		std::vector<Particle> &mu_m, 
		double lhecut, double lhmcut) {

  for(std::vector<Particle>::iterator l = mu_p.begin(); l!=mu_p.end(); ++l)
    if(Muid_mdst(l->mdstCharged()).Muon_likelihood()<=lhmcut)
      { mu_p.erase(l); --l; }
  
  for(std::vector<Particle>::iterator l = mu_m.begin(); l!=mu_m.end(); ++l)
    if(Muid_mdst(l->mdstCharged()).Muon_likelihood()<=lhmcut)
      { mu_m.erase(l); --l; }
  
  for(std::vector<Particle>::iterator l = e_p.begin(); l!=e_p.end(); ++l)
    if(eid(l->mdstCharged()).prob(3, -1, 5)<=lhecut)
      { e_p.erase(l); --l; }
  
  for(std::vector<Particle>::iterator l = e_m.begin(); l!=e_m.end(); ++l)
    if(eid(l->mdstCharged()).prob(3, -1, 5)<=lhecut)
      { e_m.erase(l); --l; }
}

void
makeBrem(std::vector<Particle> &ecl) {

  Mdst_ecl_Manager&   EclMgr    = Mdst_ecl_Manager::get_manager();
  vector<Mdst_ecl>::iterator g;
  for( g = EclMgr.begin(); g != EclMgr.end(); ++g ) {
    if ((*g).electron()==1) continue;
    double energy = (*g).energy();
    double px = energy * sin( (*g).theta() ) * cos( (*g).phi() ) ;
    double py = energy * sin( (*g).theta() ) * sin( (*g).phi() ) ;
    double pz = energy * cos( (*g).theta() )                     ;
    HepLorentzVector p4( px, py, pz, energy );
    Particle p_ecl( p4, m_ptypeGAMM );
    if( (*g).quality() != 0 ) continue;

    HepSymMatrix  errEcl( 3, 0 );
    errEcl[ 0 ][ 0 ] = (*g).error( 0 );
    errEcl[ 1 ][ 0 ] = (*g).error( 1 );
    errEcl[ 1 ][ 1 ] = (*g).error( 2 ); // Phi
    errEcl[ 2 ][ 0 ] = (*g).error( 3 );
    errEcl[ 2 ][ 1 ] = (*g).error( 4 );
    errEcl[ 2 ][ 2 ] = (*g).error( 5 ); // Theta
    
    HepSymMatrix   errCart( 4, 0 );
    HepMatrix  jacobian( 4, 3, 0 );
    
    double  cp = cos( (*g).phi() );
    double  sp = sin( (*g).phi() );
    double  ct = cos( (*g).theta() );
    double  st = sin( (*g).theta() );
    double   E = (*g).energy();
    jacobian[ 0 ][ 0 ] =       cp * st;
    jacobian[ 0 ][ 1 ] =  -E * sp * st;
    jacobian[ 0 ][ 2 ] =   E * cp * ct;
    jacobian[ 1 ][ 0 ] =       sp * st;
    jacobian[ 1 ][ 1 ] =   E * cp * st;
    jacobian[ 1 ][ 2 ] =   E * sp * ct;
    jacobian[ 2 ][ 0 ] =            ct;
    jacobian[ 2 ][ 1 ] =           0.0;
    jacobian[ 2 ][ 2 ] =  -E      * st;
    jacobian[ 3 ][ 0 ] =           1.0;
    jacobian[ 3 ][ 1 ] =           0.0;
    jacobian[ 3 ][ 2 ] =           0.0;
    errCart = errEcl.similarity( jacobian );
    p_ecl.momentum().momentum( p_ecl.p(), errCart );
    
    // Add position error matrix..(huge diagonal terms)
    HepSymMatrix dX( 3, 1 );
    Hep3Vector X;
    p_ecl.momentum().position( X, dX );
    ecl.push_back( p_ecl );
  }
}  

void
setError(Particle &gam) {

  if (!gam.mdstGamma()) return;
  VectorL Pgam=gam.p();
  Mdst_ecl g = gam.mdstGamma().ecl();

  HepSymMatrix  errEcl( 3, 0 );
  
  errEcl[ 0 ][ 0 ] = g.error( 0 );
  errEcl[ 1 ][ 0 ] = g.error( 1 );
  errEcl[ 1 ][ 1 ] = g.error( 2 ); // Phi
  errEcl[ 2 ][ 0 ] = g.error( 3 );
  errEcl[ 2 ][ 1 ] = g.error( 4 );
  errEcl[ 2 ][ 2 ] = g.error( 5 ); // Theta
    
  HepSymMatrix   errCart( 4, 0 );
  HepMatrix  jacobian( 4, 3, 0 );
  
  double  cp = cos( g.phi() );
  double  sp = sin( g.phi() );
  double  ct = cos( g.theta() );
  double  st = sin( g.theta() );
  double   E = g.energy();
  jacobian[ 0 ][ 0 ] =       cp * st;
  jacobian[ 0 ][ 1 ] =  -E * sp * st;
  jacobian[ 0 ][ 2 ] =   E * cp * ct;
  jacobian[ 1 ][ 0 ] =       sp * st;
  jacobian[ 1 ][ 1 ] =   E * cp * st;
  jacobian[ 1 ][ 2 ] =   E * sp * ct;
  jacobian[ 2 ][ 0 ] =            ct;
  jacobian[ 2 ][ 1 ] =           0.0;
  jacobian[ 2 ][ 2 ] =  -E      * st;
  jacobian[ 3 ][ 0 ] =           1.0;
  jacobian[ 3 ][ 1 ] =           0.0;
  jacobian[ 3 ][ 2 ] =           0.0;
  errCart = errEcl.similarity( jacobian );
  gam.momentum().momentum( Pgam, errCart );
  HepSymMatrix dX( 3, 1 );
  Hep3Vector   X;
  gam.momentum().position( X, dX );
}

void
makePsi(std::vector<Particle> &psi,
	std::vector<Particle> &e_p,
	std::vector<Particle> &e_m,
	std::vector<Particle> &mu_p,
	std::vector<Particle> &mu_m, 
	std::vector<Particle> &ecl, 
	const double lmcut, const double umcut) {

  for(vector<Particle>::iterator i = mu_m.begin(); i != mu_m.end(); i++){
    for(vector<Particle>::iterator j = mu_p.begin(); j != mu_p.end(); j++){
      if ((i->p()+j->p()).m()<lmcut-0.5) continue;
      kvertexfitter kvf;
      addTrack2fit(kvf,*i);
      addTrack2fit(kvf,*j);
      if (!kvf.fit()) {
	VectorL P_fit=kvf.momentum(0) + kvf.momentum(1);
	Particle p(P_fit, m_ptypePSI);
	double mass = p.mass();
	setUserInfo(p);
	p.relation().append(*i);
	p.relation().append(*j);
	p.momentum().decayVertex(kvf.get_vertex());
	dynamic_cast<UserInfo&>(p.userInfo()).vmass(mass);
	dynamic_cast<UserInfo&>(p.userInfo()).vchisq(kvf.chisq()/kvf.dgf());
	if ( mass > lmcut && mass < umcut ) psi.push_back(p);
      }
    }
  } 

  for(vector<Particle>::iterator i=e_m.begin(); i!=e_m.end(); i++){
    for(vector<Particle>::iterator j=e_p.begin(); j!=e_p.end(); j++){
      VectorL P_jpsi = ((*i).momentum().p() + (*j).momentum().p());
      
      VectorL P=i->p()+ j->p();
      VectorL Pecl;
      for(vector<Particle>::iterator k = ecl.begin();
	  k != ecl.end(); k++){
	if((*k).p().vect().mag() > 5. ) continue;
	if((*k).p().vect().angle(i->p().vect()) < 0.05) Pecl+=(*k).p();
	if((*k).p().vect().angle(j->p().vect()) < 0.05) Pecl+=(*k).p();
      }
      
      if ((P+Pecl).m()<lmcut-0.5) continue;
      kvertexfitter kvf;
      addTrack2fit(kvf,*i);
      addTrack2fit(kvf,*j);
      if (!kvf.fit()) {
	VectorL P_fit=kvf.momentum(0) + kvf.momentum(1)+Pecl;
	Particle p(P_fit, m_ptypePSI);
	double mass = p.mass();
	setUserInfo(p);
	p.relation().append(*i);
	p.relation().append(*j);
	for(vector<Particle>::iterator k = ecl.begin();
	    k != ecl.end(); k++){
	  if((*k).p().vect().mag() >  5. ) continue;
	  if((*k).p().vect().angle(i->p().vect()) < 0.05 || 
	     (*k).p().vect().angle(j->p().vect()) < 0.05) 
	    p.relation().append(*k);
	}
	p.momentum().decayVertex(kvf.get_vertex());
	dynamic_cast<UserInfo&>(p.userInfo()).vmass(mass);
	dynamic_cast<UserInfo&>(p.userInfo()).vchisq(kvf.chisq()/kvf.dgf());
	if ( mass > lmcut && mass < umcut ) psi.push_back(p);
      }
    } 
  }
}

VectorL 
reFitPsi2(Particle &p, std::vector<Particle> &ecl) {
  
  VectorL PP;
  if (p.nChildren()<2) return PP;

  kmassvertexfitter kmv;
  kmv.invariantMass(3.686);
  Particle c0=p.child(0);
  Particle c1=p.child(1);

  kmv.add_track(c0.p(), c0.x(), c0.momentum().dpx(),
		c0.charge(), c0.mass());
  kmv.add_track(c1.p(), c1.x(), c1.momentum().dpx(), 
		c1.charge(), c1.mass());
  int nfit=2;
  if (abs(c0.lund())==11)  {
    for(vector<Particle>::iterator k = ecl.begin();
	k != ecl.end(); k++){
      if (nfit>4) continue;
      if((*k).p().vect().mag() >  5.000 ) continue;
      if((*k).p().vect().angle(c0.p().vect()) < 0.05 ||
	 (*k).p().vect().angle(c1.p().vect()) < 0.05) {
	addTrack2fit(kmv,*k);
	nfit++; 
      }
    }
  }
 
  if(!kmv.fit()) {
    for(unsigned i=0; i<nfit; i++)
      PP += kmv.momentum(i); 
    return PP;
  }
  return PP;
}

VectorL 
reFitPsi2(Particle &p, Particle &p0, Particle &p1, double M_fit) {
  
  VectorL PP(0,0,0,0);

  kmassvertexfitter kmv;
  kmv.invariantMass(M_fit);
  addTrack2fit(kmv,p);
  addTrack2fit(kmv,p0);
  addTrack2fit(kmv,p1);
  if(!kmv.fit()) {
    for(unsigned i=0;i<3;++i)
      PP+=kmv.momentum(i);
  }

  return PP;
}

void
reFitPsi_make(Particle &p, std::vector<Particle> &ecl, float M_fit) {

  kmassvertexfitter kmv;
  kmv.invariantMass(M_fit);
  Particle c0=p.child(0);
  Particle c1=p.child(1);
  
  double nfit=2;
  
  addTrack2fit(kmv,c0);
  addTrack2fit(kmv,c1);
  
  if (abs(c0.lund())==11)  {
    for(vector<Particle>::iterator k = ecl.begin();
	k != ecl.end(); k++){
      if (nfit>4) continue;
      if((*k).p().vect().mag() >  5.000 ) continue;
      if((*k).p().vect().angle(c0.p().vect()) < 0.05) {
	addTrack2fit(kmv,(*k));
	nfit++; }
      if((*k).p().vect().angle(c1.p().vect()) < 0.05) {
	addTrack2fit(kmv,(*k));
	nfit++;
      }
    }
  }
  
  if(!kmv.fit()) 
    makeMother(kmv, p);
  return;
}

void
reFitPsi2(std::vector<Particle> &psi, std::vector<Particle> &ecl) {
  for(int id=0; id<psi.size(); ++id) 
    reFitPsi2(psi[id], ecl);
    return;
}

void
reFitPsi(std::vector<Particle> &psi, std::vector<Particle> &ecl) {
  
  for(int id=0; id<psi.size(); ++id) {
    float mass_fit=3.097;
    reFitPsi_make(psi[id], ecl, mass_fit);
  }
  return;
}
  
void
reFitPsiSb(std::vector<Particle> &psi, std::vector<Particle> &ecl) {

  for(int id=0; id<psi.size(); ++id) {
    float Mj =dynamic_cast<UserInfo&>(psi[id].userInfo()).vmass();
    float dwid=0.03;
    float dmass = (Mj-3.097)/2./dwid;
    int ibin = int(floor(dmass+0.5));
    float mass_fit=3.097+ibin*2.*dwid;
    reFitPsi_make(psi[id], ecl, mass_fit);
  }
}

void
reFitDSb(Particle & p,  float M_fit) {

  kmassvertexfitter kmv;
  kmv.invariantMass(M_fit);
  for(unsigned i=0; i<p.nChildren(); i++)
    addTrack2fit(kmv,p.child(i));
  int err=kmv.fit();
  if(!kmv.fit()) 
    makeMother(kmv, p);
}

#if defined(BELLE_NAMESPACE)
}
#endif
