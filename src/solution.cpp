#include "solution.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void solution::reset() {
	energy = -1.;
	energyR = -1.;
	energylow = -1.;
	energylowR = -1.;
	benergy = -1.;
	benergyR = -1.;
	benergylow = -1.;
	benergylowR = -1.;
	denergy = -1.;
	denergyR = -1.;
	denergylow = -1.;
	denergylowR = -1.;
	qdc = -1.;
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
	isSiCsI = false;
	
	// Variables filled after getPID() from Silicon.cpp
	ipid = 0;
	iZ = 0;
	iA = 0;
	mass = 0;
	
	// Variables filled after position() and calcEloss() from Silicon.cpp
	Xpos = -1;
	Ypos = -1;
	Zpos = -1;
	theta = -1;
	phi = -1;
	energyTot = -1;
	Ekin = -1;
	velocity = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void solution::SetTargetDistance(double dist0) {
	distTarget = dist0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double solution::angle() {
	double XYZ2 = (Xpos*Xpos) + (Ypos*Ypos) + (distTarget*distTarget);
	theta = acos(distTarget / sqrt(XYZ2));
	phi = atan2(Ypos , Xpos);

	return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void solution::getMomentum() {
	momentum = Kinematics.getMomentum(Ekin, mass);
	Mvect[0] = momentum * sin(theta) * cos(phi);
	Mvect[1] = momentum * sin(theta) * sin(phi);
	Mvect[2] = momentum * cos(theta);

	// scale = 1 einstein, 0 for newton
	energyTot = (Ekin * Kinematics.scale) + mass;

	velocity = momentum / energyTot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



