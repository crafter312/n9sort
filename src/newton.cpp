#include "newton.h"

#include <iostream>

using namespace std;

CNewton::CNewton():c(.9784),nMass(1.),scale(0.),CKinematics(){} 

//*********************************************************
float CNewton::getMomentum(float eKin, float mass)
{
  float mom = sqrt(2.*mass*eKin);
  return mom;
}
//********************************************************
  //transform a momentum vector to new frame
float CNewton::transformMomentum(float* mom, float *Vreference, 
				  float mass, float* momNew)
{

  float Ekin = 0.;
  for (int j=0;j<3;j++)
    {
      momNew[j] = mom[j] - mass*Vreference[j]/c;
      Ekin += pow(momNew[j],2);
    }

  Ekin /= 2.*mass;
  return Ekin;
}
//****************************************
float CNewton::gamma(float vel)
{
  return 1.;
}

