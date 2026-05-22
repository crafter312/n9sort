#include "einstein.h"

#include <iostream>

using namespace std;

CEinstein::CEinstein():c(30.),nMass(931.478),scale(1.),CKinematics() {} 

//*********************************************************
double CEinstein::getMomentum(double eKin, double mass)
{
  double pc = sqrt(pow(eKin+mass,2) - pow(mass,2));
  return pc;
}

//*********************************************************
double CEinstein::getKE(double pc, double mass)
{
  double ek = sqrt(pow(pc,2) + pow(mass,2))-mass;
  return ek;
}
//********************************************************
  /**
   *transform a momentum vector to new frame
   * and returns the new kinetic energy in MeV
   */
double CEinstein::transformMomentum(double* mom, double *Vreference, 
				  double energyTot, double* momNew)
{

  //find momentum parallel and perpendicular to transfrom velocity
  double dot = 0.;
  double VVreference = 0.;
  double perp[3];
  double para[3];
  double paraOld = 0.;
  for (int i=0;i<3;i++) 
  {
    dot += mom[i]*Vreference[i];
    VVreference += pow(Vreference[i],2);
  }
  //take projections on vreference vector (z-axis down beam)
  for (int i=0;i<3;i++) 
  {
    para[i] = dot/VVreference*Vreference[i];
    perp[i] = mom[i] - para[i];
    paraOld += pow(para[i],2);
  }

  
  paraOld = sqrt(paraOld); // magnitude of parallel mometum
  VVreference = sqrt(VVreference); // magnidtiude of velocity shift

  //transform parallel component
  double gamma = 1./sqrt(1-pow(VVreference/c,2));
  double paraNew = (paraOld - energyTot*VVreference/c)*gamma;

  // add perpendicular and new parallel components
  for (int j=0;j<3;j++)
  {
    momNew[j] = perp[j] + paraNew/paraOld*para[j];
  }


  double energyTotNew = gamma*(energyTot - VVreference*paraOld/c);
  //cout << gamma << " " << energyTot << " " << VVreference << " " << paraOld << endl;
  //cout << mom[0] << " " << mom[1] << " " << mom[2] << endl;
  //cout << Vreference[0] << " " << Vreference[1] << " " << Vreference[1] << endl;
  //cout << "a"<< endl;
  return energyTotNew;
}
//**************************************************
double CEinstein::gamma(double vel)
{
  return 1./sqrt(1-pow(vel/c,2));
}
