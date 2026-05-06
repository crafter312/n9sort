#ifndef kinematics_
#define kinematics_
#include <cmath>
#include <iostream>

/**
 * Base class for kinematics
 */

class CKinematics 
{
 public:

  CKinematics(){};
  virtual ~CKinematics(){};
  //virtual void AddVelocities(float*, float*, float, float*);
  //virtual void FindCenterOfMass(float* , float, float*, float);
  //virtual float getVelocity(float,float)=0;
  //virtual float getMomentum(float,float)=0;
  //virtual float transformMomentum(float* mom,float* vreference,float energyTot,
  //  float*momNew)=0;
  float theta;
  float phi;
  float energy;
  float velocity;
  float vcm[3];
  float velocitycm;
  float Erel;

	// Default assignment operator
	// Classes that store instances of this class can now use a default assignment operator
	CKinematics& operator=(const CKinematics&) = default;
};
#endif 
