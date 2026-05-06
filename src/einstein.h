#ifndef einstein_
#define einstein_
#include "kinematics.h"
#include <cmath>
#include <iostream>

/**
 * This class performs relativistics kinematics
 */

class CEinstein : public CKinematics
{
 public:
  double const c;
  double const nMass;
  double const scale;
  CEinstein();
  //void AddVelocities(double*, double*, double, double*);
  //void FindCenterOfMass(double* , double, double*, double);
  //double getVelocity(double,double);
  double getMomentum(double eKin, double mass);
  double getKE(double pc, double mass);
  double transformMomentum(double* mom, double* vreference, double energyTot, double* momNew);
  double gamma(double vel);

	// Simple assignment operator
	// Classes that store instances of this class can now use a default assignment operator
	CEinstein& operator=(const CEinstein& other) { 
		if (this == &other) return *this;
		CKinematics::operator=(other);
		return *this;
	}
};
#endif 
