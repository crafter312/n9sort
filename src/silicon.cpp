/* Class dealing with a single Hira Si telescope
 * Modified by Henry Webb (h.s.webb@wustl.edu) May 2026
 * to look more pretty (among other small modifications)
 */

#include "silicon.h"

#include <utility>

#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

silicon::silicon(double thick0, SortConfig& config) {
	TargetThickness = thick0;
	SiWidth = 6.45;
	losses = new CLosses(3, config);
	Ran = new TRandom();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

silicon::~silicon() {
	delete losses;
	delete Ran;
	delete Pid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void silicon::init(int id0, SortConfig& config) {
	id = id0;
	// -ND checked 5/12/2022 these distances are correct compared to the simulation
	double const XcenterA[4] = {4.419,2.819,-4.419,-2.819};
	double const YcenterA[4] = {2.819,-4.419,-2.819,4.419};
	Xcenter = XcenterA[id];
	Ycenter = YcenterA[id];

	ostringstream outstring;
	outstring << "pid_quad" << id+1;

	Pid = new pid(outstring.str(), config);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void silicon::SetTargetDistance(double dist)	{
	for (size_t i = 0; i < 20; i++) Solution[i].SetTargetDistance(dist);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void silicon::reset() {
	maxFront = 0.;
	multFront = 0;
	maxBack = 0.;
	multBack = 0;
	multDelta = 0;
	maxDelta = 0.;

	Front.reset();
	Back.reset();
	Delta.reset();
	if (Nsolution > 100) cout << "here post F,B,D reset, need to reset " << Nsolution << " solutions" << endl;
	for (size_t i = 0; i < Nsolution; i++) Solution[i].reset();
	Nsolution = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void silicon::Reduce() {
	multFront = Front.Reduce("F");
	multBack = Back.Reduce("B");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Subroutine to identify a single particle from strip data
int silicon::simpleFront() {
	int dstrip = abs(Front.Order[0].strip - Delta.Order[0].strip);
	if (dstrip < -1 && dstrip > 3) {
		Nsolution = 0;
		return 0;
	}

	if (fabs(Front.Order[0].energy - Back.Order[0].energy) > 2.) {
		Nsolution = 0;
		return 0;
	}

	double timediff = Front.Order[0].time - Delta.Order[0].time;
	//cout << "Front.Order[0].time " << Front.Order[0].time << " - Delta.Order[0].time " << Delta.Order[0].time << " = " << timediff << endl;
	//if ( timediff < -500. || timediff > 100) 
	//{
	//	Nsolution = 0;
	//	return 0;
	//}

	Solution[0].energy = Front.Order[0].energy;
	Solution[0].energyR = Front.Order[0].energyR;
	Solution[0].benergy = Back.Order[0].energy;
	Solution[0].benergyR = Back.Order[0].energyR;
	Solution[0].denergy = Delta.Order[0].energy;
	Solution[0].denergyR = Delta.Order[0].energyR;
	Solution[0].time = Front.Order[0].time;
	Solution[0].btime = Back.Order[0].time;
	Solution[0].dtime = Delta.Order[0].time;
	Solution[0].ifront = Front.Order[0].strip;
	Solution[0].iback = Back.Order[0].strip;
	Solution[0].ide = Delta.Order[0].strip;
	Solution[0].itele = id; 
	Solution[0].timediff = timediff;
	//Solution[0].Nbefore = Front.Order[i].Nbefore;
	//Solution[0].Norder = Front.Order[i].Norder;

	//cout << "ifront " << Solution[0].ifront << ", strip " << Front.Order[0].strip << endl;

	Nsolution = 1;
	return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Finds particle identification - checks to see if particle is inside of z - bananas
size_t silicon::getPID() {
	size_t pidmulti = 0;

	for (size_t isol = 0; isol < Nsolution; isol++) {
		Solution[isol].ipid = 0;

		double energy = Solution[isol].energy;
		double denergy = Solution[isol].denergy * cos(Solution[isol].theta);

		bool FoundPid = Pid->getPID(energy, denergy);

		// No particle id is found
		if (!FoundPid) continue;
		else pidmulti++;

		Solution[isol].iZ = Pid->Z;
		Solution[isol].iA = Pid->A;
		Solution[isol].mass = Pid->mass;

		size_t pidnum;
		pair<size_t, size_t> ZA(Pid->Z, Pid->A);
		if (ZA == sz_pair(1, 1)) pidnum = 1;      // proton
		else if (ZA == sz_pair(1, 2)) pidnum = 2; // deuteron
		else if (ZA == sz_pair(1, 3)) pidnum = 3; // triton
		else if (ZA == sz_pair(2, 3)) pidnum = 4; // 3He
		else if (ZA == sz_pair(2, 4)) pidnum = 5; // alpha
		else if (ZA == sz_pair(2, 6)) pidnum = 6; // 6He
		else if (ZA == sz_pair(3, 6)) pidnum = 7; // 6Li
		else if (ZA == sz_pair(3, 7)) pidnum = 8; // 7Li
		else pidnum = 10;                  // default case to be thorough, should only happen if there is a Z-line that doesn't have a case in this switch block

		Solution[isol].ipid = pidnum;
	}

	return pidmulti;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int silicon::calcEloss() {
	for (int isol=0; isol<Nsolution; isol++) {
		//need PID to calculate energy loss
		if (!Solution[isol].ipid)
		{
			Solution[isol].Ekin = 0;
			return 0;
		}

		//kinetics calc, add Delta and energy for total energy
		double sumEnergy = Solution[isol].denergy + Solution[isol].energy;
		double pc_before = sqrt(pow(sumEnergy+Solution[isol].mass,2) - pow(Solution[isol].mass,2));
		double velocity_before = pc_before/(sumEnergy+Solution[isol].mass);

		double thick = TargetThickness/2/cos(Solution[isol].theta);

		double ein = losses->getEin(sumEnergy,thick,Solution[isol].iZ,Solution[isol].mass/m0);

		//out << "loss correction " << ein - sumEnergy << endl;

		Solution[isol].Ekin = ein;
		//calc momentum vector, energyTot, and velocity
		Solution[isol].getMomentum();


		//protons can punch through at high energies
		if (Solution[isol].iA == 1 && Solution[isol].iZ == 1)
		{
			if (Solution[isol].Ekin > 15.5)
			{
				Solution[isol].iA = 0;
				Solution[isol].iZ = 0;
				Solution[isol].Ekin = 0;
				return 0;
			}
		}

		//deuterons can punch through
		if (Solution[isol].iA == 2 && Solution[isol].iZ == 1)
		{
			if (Solution[isol].Ekin > 20.5)
			{
				Solution[isol].iA = 0;
				Solution[isol].iZ = 0;
				Solution[isol].Ekin = 0;
				return 0;
			}
		}
		//Tritons can punch through
		if (Solution[isol].iA == 3 && Solution[isol].iZ == 1)
		{
			if (Solution[isol].Ekin > 24)
			{
				Solution[isol].iA = 0;
				Solution[isol].iZ = 0;
				Solution[isol].Ekin = 0;
				return 0;
			}
		}


	}




	return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Recursive subroutine used for multihit subroutine
void silicon::loop(int depth)
{

	if (depth == NestDim) //depth starts at 0
	{
		int dstrip = 0;
		double de = 0.;
		for (int i=0;i<NestDim;i++)
		{
			dstrip += abs(Delta.Order[NestArray[i]].strip - Front.Order[i].strip);
			de += abs(Back.Order[NestArray[i]].energy - Front.Order[i].energy);
		}

		//cout << " 1 " << zline[0].n << endl;

		if (dstrip < dstripMin)
		{
			dstripMin = dstrip;
			for (int i=0;i<NestDim;i++) {arrayD[i] = NestArray[i];}
		}
		//cout << " 2 " << zline[0].n << endl;

		if (de < deMin)
		{
			deMin = de;
			for (int i=0;i<NestDim;i++) {arrayB[i] = NestArray[i];}
		}
		//cout << "leave" << " " << zline[0].n << endl;
		return;
	}

	//cout << "recurse " << zline[0].n << endl;
	for (int i=0;i<NestDim;i++)
	{
		NestArray[depth] = i;
		int leave = 0;
		for (int j=0;j<depth;j++)
		{
			if (NestArray[j] == i)
			{
				leave =1;
				break; 
			}
		}
		if (leave) continue;
		loop(depth+1);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Extracts multiple particle from strip data 
int silicon::multiHit() {
	int Ntries = min(Front.Nstore,Back.Nstore);
	Ntries = min(Ntries,Delta.Nstore);

	if (Ntries > 4) Ntries =4;
	Nsolution = 0;
	if (Ntries <= 0) return 0;

	for (NestDim = Ntries;NestDim>0;NestDim--)
	{
		dstripMin = 1000;
		deMin = 10000.;

		//look for best solution
		loop(0);

		//check to see if best possible solution is reasonable
		int leave = 0;
		for (int i=0;i<NestDim;i++)
		{
			if (abs(Delta.Order[arrayD[i]].strip - Front.Order[i].strip) > 2) 
			{
				leave = 1;
				break;
			}
			if (fabs(Back.Order[arrayB[i]].energy - Front.Order[i].energy) > 2.) 
			{
				leave = 1;
				break;
			}
			double timediff = Front.Order[i].time - Delta.Order[arrayD[i]].time;
			//if ( timediff < -500. || timediff > 100) 
			//{
			//	Nsolution = 0;
			//	return 0;
			//}
		}

		if (leave) continue;
		// now load solution
		for (int i=0;i<NestDim;i++)
		{
			double timediff = Front.Order[i].time - Delta.Order[arrayD[i]].time;
			Solution[i].energy = Front.Order[i].energy;
			Solution[i].energyR = Front.Order[i].energyR;
			Solution[i].time = Front.Order[i].time;
			Solution[i].denergy = Delta.Order[arrayD[i]].energy;
			Solution[i].ifront = Front.Order[i].strip;
			Solution[i].iback = Back.Order[arrayB[i]].strip;
			Solution[i].ide = Delta.Order[arrayD[i]].strip;
			Solution[i].itele = id;
			Solution[i].timediff = timediff;
			//Solution[i].Nbefore = Front.Order[i].Nbefore;
			//Solution[i].Norder = Front.Order[i].Norder;
		}

		Nsolution = NestDim;

		break;
	}
	return Nsolution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Calculates the x-y position and angles in the array in cm
void silicon::position(int isol) {
	double Xpos,Ypos;

	if (id == 0) 
	{
		Xpos = Xcenter + (((double)Solution[isol].iback+Ran->Rndm())/32.-0.5)*SiWidth;
		Ypos = Ycenter + (((double)Solution[isol].ifront+Ran->Rndm())/32.-0.5)*SiWidth;
	}
	else if (id == 1)
	{
		Xpos = Xcenter + (((double)Solution[isol].ifront+Ran->Rndm())/32.-0.5)*SiWidth;
		Ypos = Ycenter + (0.5-((double)Solution[isol].iback+Ran->Rndm())/32.)*SiWidth;
	}
	else if (id == 2)
	{
		Xpos = Xcenter + (0.5-((double)Solution[isol].iback+Ran->Rndm())/32.)*SiWidth;
		Ypos = Ycenter + (0.5-((double)Solution[isol].ifront+Ran->Rndm())/32.)*SiWidth;
	}
	else if (id == 3)
	{
		Xpos = Xcenter + (0.5-((double)Solution[isol].ifront+Ran->Rndm())/32.)*SiWidth;
		Ypos = Ycenter + (((double)Solution[isol].iback+Ran->Rndm())/32.-0.5)*SiWidth;
	}

	//	Xpos += .3;

	Solution[isol].Xpos = Xpos;
	Solution[isol].Ypos = Ypos;
	double theta = Solution[isol].angle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Calculates the x-y position in the array in cm
void silicon::positionC(int isol) {
	double Xpos,Ypos;

	if (id == 0) 
	{
		Xpos = Xcenter + (((double)Solution[isol].iback+.5)/32.-0.5)*SiWidth;
		Ypos = Ycenter + (((double)Solution[isol].ifront+.5)/32.-0.5)*SiWidth;
	}
	else if (id == 1)
	{
		Xpos = Xcenter + (((double)Solution[isol].ifront+.5)/32.-0.5)*SiWidth;
		Ypos = Ycenter + (0.5-((double)Solution[isol].iback+.5)/32.)*SiWidth;
	}
	else if (id == 2)
	{
		Xpos = Xcenter + (0.5-((double)Solution[isol].iback+.5)/32.)*SiWidth;
		Ypos = Ycenter + (0.5-((double)Solution[isol].ifront+.5)/32.)*SiWidth;
	}
	else if (id == 3)
	{
		Xpos = Xcenter + (0.5-((double)Solution[isol].ifront+.5)/32.)*SiWidth;
		Ypos = Ycenter + (((double)Solution[isol].iback+.5)/32.-0.5)*SiWidth;
	}
	Solution[isol].Xpos = Xpos;
	Solution[isol].Ypos = Ypos;
	double theta = Solution[isol].angle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



