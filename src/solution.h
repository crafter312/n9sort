#ifndef _solution
#define _solution

#ifdef rel
#include "einstein.h"
#else
#include "newton.h"
#endif

class solution {

public:
#ifdef rel
	CEinstein Kinematics;
#else
	CNewton Kinematics;
#endif

	double distTarget;

	// Values loaded into solution class in Silicon when solution is found (contains Si-Si and Si-CsI variables)
	double energy;
	double energyR;
	double energylow;
	double energylowR;
	double benergy;
	double benergyR;
	double benergylow;
	double benergylowR;
	double denergy;
	double denergyR;
	double denergylow;
	double denergylowR;
	double qdc;
	double time;
	double timeR;
	double btime;
	double btimeR;
	double dtime;
	double dtimeR; 
	double CsITime;
	double CsITimeR;
	int ifront;
	int iback;
	int ide;
	int iCsI;
	int itele;
	double timediff;
	int Nbefore;
	int Norder;
	bool isSiCsI;

	// Variables filled from S800
	double theta_s800;
	double phi_s800;
	
	// Variables filled after getPID() from Silicon.cpp
	size_t ipid;
	size_t iZ;
	size_t iA;
	double mass;
	
	// Variables filled after position() and calcEloss() from Silicon.cpp
	double Xpos;
	double Ypos;
	double Zpos; // added for TexNeut
	double theta;
	double phi;
	double energyTot;
	double Ekin;
	double velocity;

	double Etot_cm;
	double theta_cm;
	double velocity_cm;
	//double Vvect[3];
	double Vcmvect[3];
	double Mvect[3];
	double MomCM[3];
	double momentum;
	double momentumCM;
	double MomRot[3];
	double MomRot2[3];
	double energyCM;

	double KE;

	double Vlab;

	void reset();
	void SetTargetDistance(double dist0);
	int setPID(int);
	int setEnergy(double deltaE, double E);
	double angle();
	void getMomentum();

	// assignment or copy operator
	solution& operator=(const solution& other) = default;
};

#endif
