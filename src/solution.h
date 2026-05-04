#ifndef _solution
#define _solution

#include <cmath>

#define rel

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

  float distTarget;

  //values loaded into solution class in Silicon when solution is found
  float energy;
  float energyR;
  float benergy;
  float benergyR;
  float denergy;
  float denergyR;
	float time;
  float timeR;
  float btime;
  float btimeR;
  float dtime;
  float dtimeR; 
	float CsITime;
	float CsITimeR;
  int ifront;
  int iback;
  int ide;
	int iCsI;
  int itele;
  float timediff;
  int Nbefore;
  int Norder;

  //variables filled from S800
  float theta_s800;
  float phi_s800;
  
  //variables filled after getPID() from Silicon.cpp
  int ipid;
  int iZ;
  int iA;
  float mass;
  
  //variables filled after position() and calcEloss() from Silicon.cpp
  float Xpos;
  float Ypos;
  float Zpos; // added for TexNeut
  float theta;
  float phi;
  float energyTot;
  float Ekin;
  float velocity;

  float Etot_cm;
  float theta_cm;
  float velocity_cm;
  //float Vvect[3];
  float Vcmvect[3];
  float Mvect[3];
  float MomCM[3];
  float momentum;
  float momentumCM;
  float MomRot[3];
  float MomRot2[3];
  float energyCM;

  float KE;

  float Vlab;

  void reset();
  void SetTargetDistance(double dist0);
  int setPID(int);
  int setEnergy(float deltaE, float E);
  float angle();
  void getMomentum();
};


#endif
