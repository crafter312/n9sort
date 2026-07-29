#include "correl2.h"

#include "wood.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

correl2::correl2() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void correl2::reset() {
	for (size_t i = 0; i < Nparticles; i++) particle[i]->mult = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void correl2::zeroMask() {
	for (size_t i = 0; i < Nparticles; i++) particle[i]->zeroMask();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void correl2::initWood(wood* w) {
	N = 0;
	
#ifdef ENABLE_DEBUG
	cout << "correl2::initWood Nparticles " << Nparticles << endl;
#endif
	
	for (size_t i = 0; i < Nparticles; i++) {
	
#ifdef ENABLE_DEBUG
		cout << "correl2::initWood i " << i << endl;
#endif
	
		for (int j = 0; j < particle[i]->mult; j++) {
		
#ifdef ENABLE_DEBUG
			cout << "correl2::initWood particle[i]->mult " << particle[i]->mult << endl;
#endif
		
			if (particle[i]->mask[j]) {
			
#ifdef ENABLE_DEBUG
				cout << "correl2::initWood adding particle" << endl;
#endif
			
				w->initFrag(N, particle[i]->detector);
				N++;
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void correl2::makeArray(bool flagMask, wood& w) {
	N = 0;
	for (size_t i = 0; i < Nparticles; i++) {
		for (int j = 0; j < particle[i]->mult; j++) {
			if (!flagMask || particle[i]->mask[j]) {
				frag[N] = particle[i]->Sol[j];
				w.loadFrag(N, particle[i]->detector, particle[i]->Sol[j]);
				N++;
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Legacy function form without wood class for use in external macros
// where one might wish to perform calculations without saving an output
void correl2::makeArray(bool flagMask) {
	N = 0;
	for (size_t i = 0; i < Nparticles; i++) {
		for (int j = 0; j < particle[i]->mult; j++) {
			if (!flagMask || particle[i]->mask[j]) {
				frag[N] = particle[i]->Sol[j];
				N++;
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void correl2::load(solution* fragment) {
	for (int i = 0; i < Nparticles; i++) {
		if (fragment->iZ == particle[i]->Z && fragment->iA == particle[i]->A)
		{
			if (particle[i]->mult < 6)
			{
				particle[i]->Sol[particle[i]->mult] = fragment;
				particle[i]->mult++;
			}
			break;
		}
	}
}

//**************************************************************
/**
 * Finds the total kinetic energy of the fragments
 * in their center-of-mass frame.
 */
float correl2::findErel()
{
	//first find total momentum
	for (int i=0;i<3;i++) Mtot[i] = 0.;
	float energyTot = 0.; // total energy for relativity, total mass for newton

	for (int i=0;i<N;i++) 
	{
		energyTot += frag[i]->energyTot;
		//cout << "frag " << i << " Z = " << frag[i]->iZ;
		//cout <<" energyTot = " << frag[i]->energyTot;
		//cout << " Ekin = " << frag[i]->Ekin;
		//cout << " mass = " << frag[i]->mass << endl;

		if(frag[i]->mass > 1000000) abort();
		
		for (int j=0;j<3;j++)
		{
			Mtot[j] += frag[i]->Mvect[j];
			// cout << "Frag = " << i << " Momentum " << j << " = " << Mtot[j] << endl;
		}
	}

	momentumCM = 0.;
	for (int j=0;j<3;j++) momentumCM += pow(Mtot[j],2);
	momentumCM = sqrt(momentumCM);


	//float mmc = 0.;
	//for (int i=0;i<3;i++) mmc += pow(momC[i],2);
	//mmc = sqrt(mmc);
	//cosThetaC = momC[2]/mmc;

	//PperpC = sqrt(pow(momC[0],2)+pow(momC[1],2));
	//PparaC = momC[2];
	//PtotC = sqrt(pow(momC[0],2)+pow(momC[1],2)+pow(momC[2],2));

	//velocity of cemter of mass
	velocityCM = momentumCM*Kinematics.c/energyTot;

	// cout << "energyTot = " << energyTot << endl;
	
	double velCM[3]={0.};
	for (int j=0;j<3;j++) velCM[j] = velocityCM/momentumCM*Mtot[j];
	thetaCM = acos(velCM[2]/velocityCM);
	phiCM = atan2(velCM[1],velCM[0]);

	float totalKE = 0.;
	for (int i=0;i<N;i++)
	{

		float eKinNew = Kinematics.transformMomentum(frag[i]->Mvect,velCM,frag[i]->energyTot,frag[i]->MomCM);
		frag[i]->energyCM = eKinNew - Kinematics.scale*frag[i]->mass;
		totalKE += eKinNew - Kinematics.scale*frag[i]->mass;
		check_ke = eKinNew - Kinematics.scale*frag[i]->mass;
		check_mass = Kinematics.scale*frag[i]->mass;
	}


	float mv = 0.;
	for (int i=0;i<3;i++)
	{
		mv += pow(frag[N-1]->MomCM[i],2);
	}
	mv = sqrt(mv);
	cos_thetaH = frag[N-1]->MomCM[2]/mv;

	//In case of 2p decay find angle between heavy fragment's momentum and CM momentum?
	if (N == 3)
	{
		float dot = 0.;
		float mm = 0.;
		for (int j=0;j<3;j++)
		{
			dot += frag[2]->MomCM[j]*momC[j];
			mm += pow(frag[2]->MomCM[j],2);
		}
		mm = sqrt(mm);
		cosAlphaQ = dot/mm/PtotC;
		// cos(x) = A*B/(|A| |B|) where A is heavy fragment's momentum and B is CM momentum
	}

	// cout << "Erel = " << totalKE << endl;

	return totalKE;
}

// reconstruct events based on Qvalue, make sure it is a 6He(d,n) reaction and not 6He(p,p)
float correl2::missingmass()
{

	float Brho = 1.09768;
	float TKE = 38.505;

	float Pcbeam_z = sqrt(pow(TKE + Mass_6He, 2) - pow(Mass_6He,2));
	float Mvect_beam[3] = {0,0,Pcbeam_z};
	float Mvect_miss[3] = {0,0,Pcbeam_z};

	for (int i=0;i<N;i++) 
	{
		for (int j=0;j<3;j++)
		{
			Mvect_miss[j] -= frag[i]->Mvect[j];
		}
	}

	float Emiss = TKE + Mass_6He + Mass_d - frag[0]->energyTot - frag[1]->energyTot;
	float Mmiss = pow(Emiss,2);
	for (int j=0;j<3;j++)
	{
		Mmiss -= Mvect_miss[j]*Mvect_miss[j];
	}
	Mmiss = sqrt(Mmiss);

	return Mmiss;
}


// reconstruct events based on Qvalue, make sure it is a 6He(d,n) reaction and not 6He(p,p)
float correl2::Qvalue()//Jack's version for 6He(d,n)6He+p
{

	float Brho = 1.09768;
	float TKE = 38.6;//From Brian's printout in logbook

	float Pcbeam_z = sqrt(pow(TKE + Mass_6He, 2) - pow(Mass_6He,2));
	float Mvect_beam[3] = {0,0,Pcbeam_z};
	float Mvect_miss[3] = {0,0,Pcbeam_z};

	for (int i=0;i<N;i++) 
	{
		for (int j=0;j<3;j++)
		{
			Mvect_miss[j] -= frag[i]->Mvect[j];
		}
	}

	float ptmiss = 0.;
	for (int j=0;j<3;j++)
	{
		ptmiss += Mvect_miss[j]*Mvect_miss[j];
	}
	float Enmiss = ptmiss/(2.*Mass_n);//E=p^2/2m
	float qval = TKE - Enmiss - frag[0]->energyTot - frag[1]->energyTot + Mass_6He + Mass_p;
	return qval;
}


// reconstruct events based on Qvalue, test to see if it is 6He(p,p)
float correl2::Qvalue2()//Jack's version for 6He(p,p)
{

	float Brho = 1.09768;
	float TKE = 38.6;

	float qval = TKE - frag[0]->energyTot - frag[1]->energyTot + Mass_6He + Mass_p;
	return qval;//q-value qval
}

// reconstruct events based on Qvalue, test to see if it is 6He(p,p)
float correl2::TargetEx()//test for 7Li(C,C)
{
	float Brho = 0.8331;
	float TKE = 42.8213;

	float Pc_proj = sqrt(pow(TKE+Mass_7Li,2) - pow(Mass_7Li,2));

	//cout << "Pc_proj " << Pc_proj << endl;

	float Pc_tar = Pc_proj - Kinematics.getMomentum(frag[0]->Ekin, Mass_t) - Kinematics.getMomentum(frag[1]->Ekin, Mass_alpha);

	//cout << "Pc_tar " << Pc_tar << endl;

	float KE_tar = Kinematics.getKE(Pc_tar, Mass_12C);

	//TODO
	//calc two excitation energies, one off of C12 one off of deuteron
	//make vectors




	//cout << "KE_tar " << KE_tar << endl;

	//float vel_tar = Pc_tar*Kinematics.c/E_tar;

	//cout << "vel_tar " << vel_tar << endl;

	//float KE_tar = 0.5*Mass_12C*pow(vel_tar,2);

	//cout << "KE_tar " << KE_tar << endl;

	float Ex_tar = TKE + Mass_7Li - Mass_alpha - Mass_t - frag[0]->Ekin - frag[1]->Ekin - KE_tar;
	
	//cout << "Ex_tar " << Ex_tar << endl;

	//float qval = TKE - frag[0]->energyTot - frag[1]->energyTot + Mass_alpha + Mass_t;
	return Ex_tar;
}




//***********************************************************
void correl2::getJacobi()
{

	for (int i=0;i<3;i++)
		{
			frag[i]->momentumCM = 0.;
			for (int k=0;k<3;k++) frag[i]->momentumCM += pow(frag[i]->MomCM[k],2);
			frag[i]->momentumCM = sqrt(frag[i]->momentumCM);
		}


	//alpha is the third fragment
	//first JacobiT
	float dot = 0.;
	float pp[3] = {0.};
	float PP = 0.;
	for (int k=0;k<3;k++)
		{
			pp[k] = frag[0]->MomCM[k] - frag[1]->MomCM[k];
			PP += pow(pp[k],2);
			dot += pp[k]*frag[2]->MomCM[k];
		}
	PP = sqrt(PP);
	cosThetaT = dot/PP/frag[2]->momentumCM;


	dot = 0;
	double PP1 = 0;
	double pp1[3]={0.};
	for (int k=0;k<3;k++) {
		pp1[k] = frag[0]->MomCM[k]/frag[0]->mass - frag[2]->MomCM[k]/frag[2]->mass;
		PP1 += pow(pp1[k],2);
		dot += pp1[k]*frag[1]->MomCM[k];
	}
	PP1 = sqrt(PP1);
	cosThetaY[0] = -dot/PP1/frag[1]->momentumCM;


	dot = 0;
	double PP2 = 0;
	double pp2[3]={0.};
	for (int k=0;k<3;k++) {
		pp2[k] = frag[1]->MomCM[k]/frag[1]->mass - frag[2]->MomCM[k]/frag[2]->mass;
		PP2 += pow(pp2[k],2);
		dot += pp2[k]*frag[0]->MomCM[k];
	}
	PP2 = sqrt(PP2);
	cosThetaY[1] = -dot/PP2/frag[0]->momentumCM;


	cosThetaV = (pp1[0]*pp2[0] + pp1[1]*pp2[1] + pp1[2]*pp2[2])/PP1/PP2;
}
