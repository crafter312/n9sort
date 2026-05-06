#ifndef _solution
#define _solution

#include <cmath>

#ifdef rel
#include "einstein.h"
#else
#include "newton.h"
#endif

/**
 * ~\brief find detected particles from hira data
 * 
 * This class finds the incident particles from the 
 * detected Si and CsI info
 */
class solution
{
  public:
#ifdef rel
  CEinstein Kinematics;
#else
  CNewton Kinematics;
#endif

  double distTarget;

  //values loaded into solution class in Silicon when solution is found
  double energy;
  double energyR;
  double benergy;
  double benergyR;
  double denergy;
  double denergyR;
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

  //variables filled from S800
  double theta_s800;
  double phi_s800;
  
  //variables filled after getPID() from Silicon.cpp
  size_t ipid;
  size_t iZ;
  size_t iA;
  double mass;
  
  //variables filled after position() and calcEloss() from Silicon.cpp
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
