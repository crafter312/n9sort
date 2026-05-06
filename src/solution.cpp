#include "solution.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


using namespace std;

void solution::reset()
{
  energy = -1;
  energyR = -1;
  benergy = -1;
  benergyR = -1;
  denergy = -1;
	denergyR = -1;
	time = -1;
  timeR = -1;
  btime = -1;
  btimeR = -1;
  dtime = -1;
	dtimeR = -1;
  ifront = -1;
  iback = -1;
  ide = -1;
  itele = -1;
  timediff = -100000.;
  
  //variables filled after getPID() from Silicon.cpp
  ipid = 0;
  iZ = 0;
  iA = 0;
  mass = 0;
  
  //variables filled after position() and calcEloss() from Silicon.cpp
  Xpos = -1;
  Ypos = -1;
  Zpos = -1;
  theta = -1;
  phi = -1;
  energyTot = -1;
  Ekin = -1;
  velocity = -1;
}

void solution::SetTargetDistance(double dist0)
{
  distTarget = dist0;
}

double solution::angle()
{
  double XYZ2 = pow(Xpos , 2) + pow(Ypos, 2) + pow(distTarget,2);
  theta = acos(distTarget/sqrt(XYZ2));
  phi = atan2(Ypos , Xpos);

  return theta;
}
//********************************************************
void solution::getMomentum()
{
  momentum = Kinematics.getMomentum(Ekin,mass);
  Mvect[0] = momentum*sin(theta)*cos(phi);
  Mvect[1] = momentum*sin(theta)*sin(phi);
  Mvect[2] = momentum*cos(theta);
  

  //scale = 1 einstein, 0 for newton
  energyTot = Ekin*Kinematics.scale + mass;

  velocity = momentum/energyTot;
}
