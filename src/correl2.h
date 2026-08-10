#ifndef correl2_H
#define correl2_H

#include <vector>

#include "constants.h"
#include "parType.h"

class wood;

class correl2 {

public:
	correl2();

#if defined(rel) && rel == 1
	CEinstein Kinematics;
#else
	CNewton Kinematics;
#endif

	void zeroMask();
	void reset();
	void initWood(wood* w);
	void makeArrayAndOutput(bool flagMask, wood& w);
	void makeArray(bool flagMask);
	void makeOutput(bool flagMask, wood& w);
	void load(solution* fragment);
	void getJacobi();
	float findErel();
	float missingmass();
	float Qvalue();
	float Qvalue2();
	float TargetEx();
	int N;
	solution* frag[7];

	// Particles and particle list
	parType neutron{0, 1, parType::detType::TexNeut};
	parType proton{1, 1, parType::detType::Gobbi};
	parType H2{1, 2, parType::detType::Gobbi};
	parType H3{1, 3, parType::detType::Gobbi};
	parType H3_fake{1, 3, parType::detType::Gobbi};
	parType He3{2, 3, parType::detType::Gobbi};
	parType alpha{2, 4, parType::detType::Gobbi};
	parType He6{2, 6, parType::detType::Gobbi};
	parType Li6{3, 6, parType::detType::Gobbi};
	parType Li7{3, 7, parType::detType::Gobbi};
	parType Li8{3, 8, parType::detType::Gobbi};
	parType Li9{3, 9, parType::detType::Gobbi};
	parType Be7{4, 7, parType::detType::Gobbi};
	parType Be9{4, 9, parType::detType::Gobbi};
	const std::vector<parType*> particle{&neutron, &proton, &H2, &H3, &H3_fake, &He3, &alpha, &He6, &Li6, &Li7, &Li8, &Li9, &Be7, &Be9};
	const size_t Nparticles{particle.size()};

	// Various kinematic values
	float Vcmframevector[3];
	float Vcmframe;
	float VcmframePerp;
	float VcmframePara;
	float TKEL;
	float Erel;
	float KEToT;
	float thetaCM;
	float phiCM;
	float thetaInReactCM;

	int ifront;
	int istrip;
	int itele;
	float Vlab;

	float cos_thetaH;

	float thetaAlpha;
	float phiAlpha;
	float Thetaopen;
	float ppThetaRel;
	float cosThetapp;
	float pTheta0;
	float pTheta2;
	float pTheta3;
	float pTheta1;
	float aTheta0;
	float pThetaMax;
	float aTheta1;
	float thetadip;
	float Epp;
	float CosAngle;

	float CosAngleY1;
	float CosAngleY2;
	float Ey1;
	float Ey2;

	float cosBeamCMangle; //-ND

	float Ey1_10;
	float Ey2_10;
	float CosAngleAlphaY1_10;
	float CosAngleAlphaY2_10;

	float Epp_10;
	float CosAngle_10;

	float CosAngleAlpha;
	float CosAngleAlphaY1;
	float CosAngleAlphaY2;
	float Eaa;
	float ThetaP;
	float ThetaA1;
	float ThetaA2;
	float PhiP;
	float PhiA1;
	float PhiA2;
	float cosAA;

	float aarelbeamangle;
	float prelbeamangle;

	float Eap_min, Eap_max;

	float ThetaA1_P, ThetaA2_P;

	float cosAlpha[2];
	float cosProton0[2];
	float cosProton1[2];
	float energy3;

	float velocityCM;
	float momentumCM;

	// Jacobi

	float cosThetaT;
	float cosThetaY[2];
	float cosThetaV;
	float getAlphaMom();
	float x_T;
	float x_Y[2];
	void rotate();

	float momC[3];
	float Mtot[3];
	float cosThetaC;
	float PperpC;
	float PparaC;
	float PtotC;
	float cosAlphaQ;

};

#endif
