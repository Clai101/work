#include "belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

/* Ptype m_ptypeGAMM("GAMM"); */
/* Ptype m_ptypeRHO0("RHO0"); */
/* Ptype m_ptypePHI("PHI"); */
/* Ptype m_ptypeETAC("ETAC"); */
/* Ptype m_ptypeCHI0("CHI0"); */
/* Ptype m_ptypePSI("PSI"); */
/* Ptype m_ptypePSI2("PSI2"); */
/* Ptype m_ptypeDP("D+"); */
/* Ptype m_ptypeDM("D-"); */
/* Ptype m_ptypeD0("D0"); */
/* Ptype m_ptypeD0B("D0B"); */
/* Ptype m_ptypeB0("B0"); */
/* Ptype m_ptypeB0B("B0B"); */
/* Ptype m_ptypeDstarP("D*+"); */
/* Ptype m_ptypeDstarM("D*-"); */
/* Ptype m_ptypeDSP("DS+"); */
/* Ptype m_ptypeDSM("DS-"); */
/* Ptype m_ptypeDSstarP("DS*+"); */
/* Ptype m_ptypeDSstarM("DS*-"); */

/* Ptype m_ptypeLAMC("LAMC"); */
/* Ptype m_ptypeALAMC("ALAMC"); */
/* Ptype m_ptypeUPS4("UPS4"); */

void withDrDzCut(std::vector<Particle> &list, double dr, double dz);
void withEminCut(std::vector<Particle> &gam, double Emin);
void withEminCutPi0(std::vector<Particle> &pi0, double Emin);
void withPStarMinCut(std::vector<Particle> &in, double Psmin);

void withProtonIdCut(std::vector<Particle> &pr, 
		     std::vector<Particle> &apr, 
		     double lhpcut);
void withKaonIdCut(std::vector<Particle> &k_p, 
		   std::vector<Particle> &k_m, 
		   double lhkcut);
void withPionIdCut(std::vector<Particle> &pi_p, 
		   std::vector<Particle> &pi_m, 
		   double lhpcut);

void withLeptonIdCut(std::vector<Particle> &e_p, 
		     std::vector<Particle> &e_m, 
		     std::vector<Particle> &mu_p, 
		     std::vector<Particle> &mu_m, 
		     double lhecut, double lhmcut); 

void makeBrem(std::vector<Particle> &gam);

void makePsi(std::vector<Particle> &psi, 
	     std::vector<Particle> &e_e,
	     std::vector<Particle> &e_p,
	     std::vector<Particle> &mu_m,
	     std::vector<Particle> &mu_p,
	     std::vector<Particle> &gam, 
	     const double lmcut, const double umcut);

void goodLepton( std::vector<Particle> &e_m,
		 std::vector<Particle> &e_p,
		 std::vector<Particle> &mu_m,
		 std::vector<Particle> &mu_p);

void reFitD(std::vector<Particle> &Dsp);

int checkCut(Particle &Dst, HepPoint3D ip_position);

void reFitPsi_make(Particle &p, std::vector<Particle> &ecl);
void reFitPsi_make(Particle &p, std::vector<Particle> &ecl, double M_fit);

void reFitPsi(std::vector<Particle> &psi, std::vector<Particle> &ecl);
void reFitPsi2(std::vector<Particle> &psi, std::vector<Particle> &ecl);
void reFitPsiSb(std::vector<Particle> &psi, std::vector<Particle> &ecl);

VectorL reFitPsi(Particle &p, std::vector<Particle> &ecl);

VectorL reFitPsi2(Particle &p, std::vector<Particle> &ecl);
VectorL reFitPsi2(Particle &p, Particle &p0, Particle &p1, double M_fit);

void reFitDSb(Particle & p,  float M_fit);

#if defined(BELLE_NAMESPACE)
}
#endif
