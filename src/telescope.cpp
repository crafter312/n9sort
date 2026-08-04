/* Class dealing with a single Hira Si telescope
 * Modified by Henry Webb (h.s.webb@wustl.edu) May 2026
 * to look more pretty (among other small modifications).
 * CsI option is designed to be generalized, should work
 * for both New Gobbi with 4 CsI per telescope and
 * Gobbi 28 with 7 CsI per telescope
 */

//       .-.      _______                             .  '  *   .  . '
//      {}``; |==|_______D                              .  * * -+-  
//      / ('        /|\                             .    * .    '  *
//  (  /  |        / | \                                * .  ' .  . 
//   \(_)_]]      /  |  \                            *   *  .   .
//                                                     '   *

// Upstream view of Gobbi
//        ____
//       |    |____                  
//       | 4  |    |              
//      _|____| 1  |
//     |    |_|____|
//     | 3  |    |
//     |____| 2  |
//          |____|

#include "telescope.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "algorithms.h"
#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

telescope::telescope(double thick0, SortConfig& config, histo& h, bool csi) : alThick(config.GetAlThick()), Histo(h), hasCsI(csi), Front(h), Back(h), Delta(h), CsI(h) {
	TargetThickness = thick0;
	SiWidth = 6.45;  // cm
	SiFrame = 7.237; // cm
	holeSize = config.GetGobbiHoleSize() * .1; // convert from mm to cm
	losses = new CLosses(6, config.GetLossDir(), config.GetTargetSuffix());
	Allosses = new CLosses(6, config.GetLossDir(), "Al");
	
	// Seed TRandom with current system clock
	//auto now = chrono::system_clock::now();
	//UInt_t tstamp = chrono::duration_cast<chrono::seconds>(now.time_since_epoch()).count();
	//Ran = new TRandom(tstamp);
	
	// Use default seed to be reproducible
	Ran = new TRandom();
	
	// Read in front/back CsI strip extents from file
	if (hasCsI) {
		string inextentsfile = config.GetConfigDir() + config.GetCsIStripExtentsFile();
		ifstream inextents(inextentsfile);
		if (inextents.fail()) throw invalid_argument(string(BOLDRED) + string("CsI silicon strip extents file ") + inextentsfile + string(" does not exist or failed to open") + string(RESET));
		
#ifdef ENABLE_DEBUG
		else cout << GREEN << "CsI silicon strip extents file " << inextentsfile << " opened" << RESET << endl;
#endif
		
		size_t Fmin, Fmax, Bmin, Bmax;
		while (inextents >> Fmin >> Fmax >> Bmin >> Bmax) {
			if (Fmin >= Fmax) throw invalid_argument(string(BOLDRED) + string("ERROR: Fmin >= Fmax from " + inextentsfile + string(RESET)));
			if (Bmin >= Bmax) throw invalid_argument(string(BOLDRED) + string("ERROR: Bmin >= Bmax from " + inextentsfile + string(RESET)));
			CsIFextents.emplace_back(Fmin, Fmax);
			CsIBextents.emplace_back(Bmin, Bmax);
			CsIFmids.emplace_back(((CsIFextents[NCsI].second - CsIFextents[NCsI].first) / (size_t)2) + CsIFextents[NCsI].first);
			CsIBmids.emplace_back(((CsIBextents[NCsI].second - CsIBextents[NCsI].first) / (size_t)2) + CsIBextents[NCsI].first);
			NCsI++;
		}
	}
	
	string calDir = config.GetCalDir();
	if (hasCsI) {
		calCsI_d     = new calibrate(4, NCsI, calDir + config.GetCsIEdcalFile(), 1, false);
		calCsI_t     = new calibrate(4, NCsI, calDir + config.GetCsIEtcalFile(), 1, false);
  		calCsI_Alpha = new calibrate(4, NCsI, calDir + config.GetCsIEalphacalFile(), 1, false);
	}
	else {
		calCsI_d     = nullptr;
		calCsI_t     = nullptr;
		calCsI_Alpha = nullptr;
	}
	
	PidCsI.resize(NCsI);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

telescope::~telescope() {
	delete losses;
	delete Allosses;
	delete Ran;
	
	if (hasCsI) for (size_t i = 0; i < NCsI; i++) if (PidCsI[i] != nullptr) delete PidCsI[i];
	else if (Pid != nullptr) delete Pid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Assume that is OldGobbi if no CsI, adjustible Gobbi 28 if yes CsI
// see above for Gobbi upstream view diagram
void telescope::init(int id0, SortConfig& config) {

#ifdef ENABLE_DEBUG
	cout << "id0 " << id0 << endl;
#endif

	if (id0 < 0 || id0 > 3) throw invalid_argument(string(BOLDRED) + string("ERROR: telescope id0 must be 0, 1, 2, or 3") + string(RESET));
	id = id0;
	
	// -ND checked 5/12/2022 these distances are correct compared to the simulation
	double XcenterA[4] = {4.419,2.819,-4.419,-2.819};
	double YcenterA[4] = {2.819,-4.419,-2.819,4.419};
	
#ifdef ENABLE_DEBUG
	cout << "telescope::init 1" << endl;
#endif
	
	ostringstream outstring;
	if (hasCsI) {
		double halfHole = holeSize * .5;
		double halfFrame = SiFrame * .5;
		XcenterA[0] = halfFrame + halfHole;
		XcenterA[1] = halfFrame - halfHole;
		XcenterA[2] = -halfFrame - halfHole;
		XcenterA[3] = -halfFrame + halfHole;
		YcenterA[0] = halfFrame - halfHole;
		YcenterA[1] = -halfFrame - halfHole;
		YcenterA[2] = -halfFrame + halfHole;
		YcenterA[3] = halfFrame + halfHole;
		
		for (size_t i = 0; i < NCsI; i++) {
			outstring.str("");
			outstring << "pid_" << id << "_" << i;
			try {
				PidCsI[i] = new pid(outstring.str(), config);
			}
			catch (...) {
				PidCsI[i] = nullptr;
				cout << NCsI << endl;
				cout << BOLDRED << "zline file " << outstring.str() << ".zline failed to load" << RESET << endl;
			}
		}
		Pid = nullptr;
	}
	else {
		outstring << "pid_quad" << id;
		try {
			Pid = new pid(outstring.str(), config);
		}
		catch (...) {
			Pid = nullptr;
		}
		for (size_t i = 0; i < NCsI; i++) PidCsI[i] = nullptr;
	}
	
#ifdef ENABLE_DEBUG
	cout << "telescope::init 2" << endl;
#endif
	
	Xcenter = XcenterA[id];
	Ycenter = YcenterA[id];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void telescope::SetTargetDistance(double dist)	{
	for (size_t i = 0; i < 20; i++) Solution[i].SetTargetDistance(dist);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void telescope::reset() {
	multFront = 0;
	multBack = 0;
	multDelta = 0;
	multCsI = 0;

	Front.reset();
	Back.reset();
	Delta.reset();
	CsI.reset();
	if (Nsolution > 100) cout << "here post F,B,(D|CsI) reset, need to reset " << Nsolution << " solutions" << endl;
	for (size_t i = 0; i < Nsolution; i++) Solution[i].reset();
	tempSol.reset();
	Nsolution = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void telescope::Reduce() {
	multFront = Front.Reduce("F");
	multBack = Back.Reduce("B");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Subroutine to identify a single particle from strip data
int telescope::simpleFront() {
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
	//if ( timediff < -500. || timediff > 100) 
	//{
	//	Nsolution = 0;
	//	return 0;
	//}

#ifdef ENABLE_DEBUG
	cout << "Front.Order[0].time " << Front.Order[0].time << " - Delta.Order[0].time " << Delta.Order[0].time << " = " << timediff << endl;
#endif

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
	//Solution[0].Nbefore = Front.Order[0].Nbefore;
	//Solution[0].Norder = Front.Order[0].Norder;
	Nsolution = 1;

#ifdef ENABLE_DEBUG
	cout << "ifront " << Solution[0].ifront << ", strip " << Front.Order[0].strip << endl;
#endif

	return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Subroutine to identify a position map from alpha calibrations;
// it stores the answer in a dedicated temporary solution.
int telescope::testingHitE() {
	tempSol.energy = Front.Order[0].energy;
	tempSol.energylow = Front.Order[0].energylow;
	tempSol.energylowR = Front.Order[0].energyRlow;
	tempSol.energyR = Front.Order[0].energyR;
	tempSol.benergy = Back.Order[0].energy;
	tempSol.benergyR = Back.Order[0].energyR;
	tempSol.denergy = -1;
	tempSol.denergyR = -1;;
	tempSol.ifront = Front.Order[0].strip;
	tempSol.iback = Back.Order[0].strip;
	tempSol.ide = -1;
	tempSol.iCsI = -1;
	tempSol.itele = id;
	tempSol.timediff = -1000000.0;
	tempSol.isSiCsI = false;
	
	double Xpos,Ypos;
	if (id == 0) {
		Xpos = Xcenter + (((double)tempSol.iback+Ran->Rndm())/32.-0.5)*SiWidth;
		Ypos = Ycenter + (((double)tempSol.ifront+Ran->Rndm())/32.-0.5)*SiWidth;
	}
	else if (id == 1) {
		Xpos = Xcenter + (((double)tempSol.ifront+Ran->Rndm())/32.-0.5)*SiWidth;
		Ypos = Ycenter + (0.5-((double)tempSol.iback+Ran->Rndm())/32.)*SiWidth;
	}
	else if (id == 2) {
		Xpos = Xcenter + (0.5-((double)tempSol.iback+Ran->Rndm())/32.)*SiWidth;
		Ypos = Ycenter + (0.5-((double)tempSol.ifront+Ran->Rndm())/32.)*SiWidth;
	}
	else if (id == 3) {
		Xpos = Xcenter + (0.5-((double)tempSol.ifront+Ran->Rndm())/32.)*SiWidth;
		Ypos = Ycenter + (((double)tempSol.iback+Ran->Rndm())/32.-0.5)*SiWidth;
	}

	tempSol.Xpos = Xpos;
	tempSol.Ypos = Ypos;
	tempSol.angle();
	
	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Subroutine for simple E-CsI matching (1 each F, B, and CsI)
int telescope::simpleECsI() {
	if (!hasCsI) return 0;

	// Here we check the CsI is behind the x-y potions of the E silicon
	// If it doesn't match, then we return early and do not create a solution
	int Fstrip = Front.Order[0].strip;
	int Bstrip = Back.Order[0].strip;
	if (CsI.Order[0].strip < 0) throw invalid_argument(string(BOLDRED) + string("ERROR: CsI strip # cannot be negative") + string(RESET));
	size_t CsIstrip = CsI.Order[0].strip;
	if ((Fstrip < CsIFextents[CsIstrip].first) || (Fstrip > CsIFextents[CsIstrip].second) || (Bstrip < CsIBextents[CsIstrip].first) || (Bstrip > CsIBextents[CsIstrip].second)) {
		Nsolution = 0;
		return 0;
	}

	// This won't work until your calibrations are good. but it should be turned on
	// Essentially, we want to ensure that the front and back hits are also good
	if (fabs(Front.Order[0].energy - Back.Order[0].energy) > 10.) {
		Nsolution = 0;
		return 0;
	}
	
	// CsI time gate (calibration shifts peaks to zero)
	//if (CsI.Order[0].time < -100 || CsI.Order[0].time > 100) {
	//	Nsolution = 0;
	//	return 0;
	//}

	double timediff = CsI.Order[0].time - Front.Order[0].time;
	Solution[0].energy = CsI.Order[0].energy;
	Solution[0].energyR = CsI.Order[0].energyR;
	Solution[0].energylow = CsI.Order[0].energy; //The later code uses high gain and low gain for gobbi
	Solution[0].energylowR = CsI.Order[0].energyR; //To make the code simpler, just clone CsI high into low gain
	Solution[0].benergy = Back.Order[0].energy;
	Solution[0].benergyR = Back.Order[0].energyR;
	Solution[0].denergy = Front.Order[0].energy;
	Solution[0].denergylow = Front.Order[0].energylow;
	Solution[0].denergyR = Front.Order[0].energyR;
	Solution[0].denergylowR = Front.Order[0].energyRlow;
	Solution[0].qdc = CsI.Order[0].qdc;

	Solution[0].ifront = Front.Order[0].strip;
	Solution[0].iback = Back.Order[0].strip;
	Solution[0].ide = -1;
	Solution[0].iCsI= CsI.Order[0].strip;
	Solution[0].itele = id; 
	Solution[0].CsITime = CsI.Order[0].time;
	Solution[0].timediff = timediff;
	Solution[0].time = Front.Order[0].time;
	Solution[0].btime = Back.Order[0].time;
	Solution[0].isSiCsI = true;
	Solution[0].itele = id;
	Nsolution = 1;
	
#ifdef ENABLE_DEBUG
	cout << "telescope::simple iCsI = " <<  Solution[0].iCsI << " itele = " << id << "solution = 0" << endl;
#endif
  
	return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Finds particle identification - checks to see if particle is inside of z - bananas
size_t telescope::getPID() {
	size_t pidmulti = 0;
	for (size_t isol = 0; isol < Nsolution; isol++) {
		Solution[isol].ipid = 0;
		double denergy = Solution[isol].denergy * cos(Solution[isol].theta);

		double energy;
		bool FoundPid;
		pid* pidtemp;
		bool isSiCsI = Solution[isol].isSiCsI;
		if (hasCsI && isSiCsI) {
			if (PidCsI[Solution[isol].iCsI] == nullptr) continue;
			energy = Solution[isol].energyR;
			pidtemp = PidCsI[Solution[isol].iCsI];
			FoundPid = pidtemp->getPID(energy, denergy);
		}
		else {
			if (Pid == nullptr) return 0;
			energy = Solution[isol].energy;
			pidtemp = Pid;
			FoundPid = pidtemp->getPID(energy, denergy);
		}

		// No particle id is found
		if (!FoundPid) continue;
		else pidmulti++;

		Solution[isol].iZ = pidtemp->Z;
		Solution[isol].iA = pidtemp->A;
		Solution[isol].mass = pidtemp->mass;

		size_t pidnum;
		pair<size_t, size_t> ZA(pidtemp->Z, pidtemp->A);
		if (ZA == sz_pair(1, 1)) pidnum = 1;      // proton
		else if (ZA == sz_pair(1, 2)) pidnum = 2; // deuteron
		else if (ZA == sz_pair(1, 3)) pidnum = 3; // triton
		else if (ZA == sz_pair(2, 3)) pidnum = 4; // 3He
		else if (ZA == sz_pair(2, 4)) pidnum = 5; // alpha
		else if (ZA == sz_pair(2, 6)) pidnum = 6; // 6He
		else if (ZA == sz_pair(3, 6)) pidnum = 7; // 6Li
		else if (ZA == sz_pair(3, 7)) pidnum = 8; // 7Li
		else pidnum = 10;                         // default case to be thorough, should only happen if there is a Z-line that doesn't have a case in this switch block

		Solution[isol].ipid = pidnum;
		
		// At this point, also apply CsI specific calibrations as function of PID
		if (hasCsI && isSiCsI)
			Solution[isol].energy = light2energy(ZA.first, ZA.second, Solution[isol].itele, Solution[isol].iCsI, Solution[isol].energy);
	}

	return pidmulti;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int telescope::calcEloss() {
	pidSkipped = 0;
	for (size_t isol = 0; isol < Nsolution; isol++) {

		// Need PID to calculate energy loss
		if (Solution[isol].ipid < 1) {
			Solution[isol].Ekin = 0;
			continue;
		}

		// Kinetics calc, add Delta and energy for total energy
		double sumEnergy, alThickTh, alEin, thick, ein;
		stringstream ss;
		try {
			sumEnergy = Solution[isol].denergy + Solution[isol].energy;
			alThickTh = alThick / cos(Solution[isol].theta);
			alEin = Allosses->getEin(sumEnergy, alThickTh, Solution[isol].iZ, Solution[isol].mass / m0);
			thick = (.5 * TargetThickness) / cos(Solution[isol].theta);
			ein = losses->getEin(alEin, thick, Solution[isol].iZ, Solution[isol].mass / m0);
/*
			if (Ran->Rndm() < 0.009) {
				cout << "\n===================================================="
				     << "\nZ " << Solution[isol].iZ << ", A " << Solution[isol].iA
				     << "\nCsI + Si front energy " << sumEnergy << " MeV"
				     << "\nAl absorber thickness " << alThickTh << " mg/cm^2"
				     << "\nPre Al absorber " << alEin << " MeV"
				     << "\nBe target thickness " << thick << " mg/cm^2"
				     << "\nFinal energy " << ein << " MeV"
				     << "\n====================================================" << endl;
				abort();
			}
*/
		}
		
		// If Eloss fails, for now invalidate solution and increment counter
		catch (const exception& e) {
			
#ifdef ENABLE_DEBUG
			cout << e.what() << endl;
#endif
			
			Solution[isol].ipid = 0;
			Solution[isol].Ekin = 0;
			pidSkipped++;
			continue;
		}

#ifdef ENABLE_DEBUG
		cout << "loss correction " << ein - sumEnergy << endl;
#endif

		Solution[isol].Ekin = ein;

		// Calc momentum vector, energyTot, and velocity
		Solution[isol].getMomentum();

		if (hasCsI) continue;
		
		// TODO: Double check punch through energies with Lise++; Should angle correct sumEnergy?

		// Protons can punch through at high energies
		if (Solution[isol].iA == 1 && Solution[isol].iZ == 1) {
			if (sumEnergy > 15.5) {
				Solution[isol].ipid = 0;
				Solution[isol].iA = 0;
				Solution[isol].iZ = 0;
				Solution[isol].Ekin = 0;
				pidSkipped++;
				continue;
			}
		}

		// Deuterons can punch through
		if (Solution[isol].iA == 2 && Solution[isol].iZ == 1) {
			if (sumEnergy > 20.5) {
				Solution[isol].ipid = 0;
				Solution[isol].iA = 0;
				Solution[isol].iZ = 0;
				Solution[isol].Ekin = 0;
				pidSkipped++;
				continue;
			}
		}
		
		// Tritons can punch through
		if (Solution[isol].iA == 3 && Solution[isol].iZ == 1) {
			if (sumEnergy > 24) {
				Solution[isol].ipid = 0;
				Solution[isol].iA = 0;
				Solution[isol].iZ = 0;
				Solution[isol].Ekin = 0;
				pidSkipped++;
				continue;
			}
		}
	}

	return pidSkipped;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Recursive subroutine used for multihit subroutine for fronts, backs, and deltas
void telescope::loopDEE(int depth) {

	if (depth == NestDim) { // depth starts at 0
		int dstrip = 0;
		double de = 0.;
		for (size_t i = 0; i < NestDim; i++) {
			dstrip += abs(Delta.Order[NestArray[i]].strip - Front.Order[i].strip);
			de += abs(Back.Order[NestArray[i]].energy - Front.Order[i].energy);
		}

		if (dstrip < dstripMin){
			dstripMin = dstrip;
			for (size_t i = 0; i < NestDim; i++) arrayD[i] = NestArray[i];
		}

		if (de < deMin) {
			deMin = de;
			for (size_t i = 0; i < NestDim; i++) arrayB[i] = NestArray[i];
		}

		return;
	}

	for (size_t i = 0; i < NestDim; i++) {
		NestArray[depth] = i;
		bool leave = false;
		for (size_t j = 0; j < depth; j++) {
			if (NestArray[j] == i) {
				leave = true;
				break; 
			}
		}
		if (leave) continue;
		loopDEE(depth + 1);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Extracts multiple particle from strip data for fronts, backs, and deltas
int telescope::multiHit() {
	int Ntries = min(Front.Nstore, Back.Nstore);
	Ntries = min(Ntries, Delta.Nstore);

	if (Ntries > 4) Ntries = 4;
	Nsolution = 0;
	if (Ntries <= 0) return 0;

	for (NestDim = Ntries; NestDim > 0; NestDim--) {
		dstripMin = 1000;
		deMin = 10000.;

		// Look for best solution
		loopDEE(0);

		// Check to see if best possible solution is reasonable
		bool leave = false;
		for (size_t i = 0; i < NestDim; i++) {
			if (abs(Delta.Order[arrayD[i]].strip - Front.Order[i].strip) > 2) {
				leave = true;
				break;
			}
			if (fabs(Back.Order[arrayB[i]].energy - Front.Order[i].energy) > 2.) {
				leave = true;
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

		// Now load solution
		for (size_t i = 0; i < NestDim; i++) {
			double timediff = Front.Order[i].time - Delta.Order[arrayD[i]].time;
			Solution[i].energy = Front.Order[i].energy;
			Solution[i].energyR = Front.Order[i].energyR;
			Solution[i].time = Front.Order[i].time;
			Solution[i].btime = Back.Order[arrayB[i]].time;
			Solution[i].dtime = Delta.Order[arrayD[i]].time;
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

// Recursive subroutine used for multihit subroutine for just fronts and backs (used for CsI matching).
// This version is from Johnathan Phillips' (j.s.phillips@wustl.edu) 22Si code and has deltas removed,
// along with a whole bunch of extra comments which I have left in. Originally used in below in the 
// multiHitECsI function in combination with the below getCombinationMasks function to match Si and
// CsI hits, but now replaced with Hungarian algorithm. Kept just because.
void telescope::loopE(int depth) {

	// Only work in this section if we are at the max level of recursion
	if (depth == NestDim) {

		// As an example for NestDim = 2
		// For first time through on loop(2) nestarray = {0,1} so we check
		// if highest energy delta matches with highest energy Front. Then 
		// check if second highest delta matches highest energy Front.

		double de = 0.;
		for (size_t i = 0; i < NestDim; i++) {

			// Difference in Front and back energy is how to match FrontE and BackE - need two types for Gobbi
			// NOTE: this needs high and low calibrations to work
			//TODO THIS WILL NOT WORK WITHOUT LOW GAIN CALIBRATIONS
			//cout << "Checking strips ifront=" << Front.Order[indToOrdIndF[i]].strip << " and iback=" << Back.Order[indToOrdIndB[NestArray[i]]].strip << endl;
			de += abs(Back.Order[indToOrdIndB[NestArray[i]]].energy - Front.Order[indToOrdIndF[i]].energy);
			//else de += abs(Back.Order[NestArray[i]].energylow - Front.Order[i].energylow);
		}

		// Here if it is the lowest total difference in strip# or energy, it is saved in
		// arrayB (for matching front-back)
		if (de < deMin) {
			deMin = de;
			for (size_t i = 0; i < NestDim; i++) {
				arrayF[i] = indToOrdIndF[i];
				arrayB[i] = indToOrdIndB[NestArray[i]];
			}
		}
		return;
	}

	// This section handles how deep we go into the recursion loop.
	// The key to this section if figuring out what NestArrays to check.
	// For NestDim = 2, we want to check {0,1} and {1,0}
	for (size_t i = 0; i < NestDim; i++) {
		NestArray[depth] = i;
		bool leave = false;

		// When matching we want to skip items already matched 
		for (size_t j = 0; j < depth; j++) {
			if (NestArray[j] == i) {
				leave = true;
				break; 
			}
		}
		if (leave) continue;
		loopE(depth + 1);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Originally used in below in the multiHitECsI function in combination with the
// above loopE function to match Si and CsI hits, but now replaced with Hungarian
// algorithm. Kept just because.
vector<vector<bool>> telescope::getCombinationMasks(size_t totalElements, size_t nestDim) {
	vector<vector<bool>> masks;
	if (nestDim > totalElements) return masks;
	
	// Create initial mask with totalElements # of elements, nestDim of which are true
	vector<bool> mask(totalElements, false);
	for (size_t i = 0; i < nestDim; i++)
		mask[i] = true;
		
	// Loop through all combinations and add to list of masks.
	// Despite the name `permutation` in the function below, this still works because
	// of the presence of booleans in the mask instead of other types.
	do {
		masks.push_back(mask);
	} while (prev_permutation(mask.begin(), mask.end()));
	
	return masks;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Modification of multiHitdEE() to match E-CsI events, stolen from Johnathan Phillips'
// (j.s.phillips@wustl.edu) 22Si sort code. Now modified to use the Hungarian algorithm
// to solve what is apparently just a biased assignment problem. This is mostly the
// same as the recursive function above + combinatorics, but offers a massive speed
// boost for events with large Si strip multiplicities.
int telescope::multiHitECsI(stringstream& ss) {

	// Step 1: Calculate poison value
	double S = Front.Nstore + Back.Nstore;     // padded cost array size
	double MaxEDiff = 10;                      // maximum front-back Si strip energy difference
	double POISON = (S * S * MaxEDiff) + 1000; // make sure poison value is sufficiently large

	// Step 2: Prepare cost array
	// This is done by creating an array of size S and padding the edge region around the main
	// Front.Nstore x Back.Nstore region with the cutoff value for the abs(back - front) scoring
	// calculation. Strip pairs in the main region that are either above the cutoff or not in
	// front of an activated CsI crystal have their costs "poisoned" with some large number.
	vector<vector<double>> C(S, vector<double>(S, MaxEDiff));
	bool inCsI;
	int ifront, iback, icsi;
	double dc;
	for (size_t w = 0; w < Front.Nstore; w++) {
		for (size_t j = 0; j < Back.Nstore; j++) {
			inCsI = false;
			ifront = Front.Order[w].strip;
			iback = Back.Order[j].strip;
			for (size_t i = 0; i < CsI.Nstore; i++) {
				icsi = CsI.Order[i].strip;
				if ((ifront < CsIFextents[icsi].first) || (ifront > CsIFextents[icsi].second) || (iback < CsIBextents[icsi].first) || (iback > CsIBextents[icsi].second)) continue;
				inCsI = true;
				break;
			}
			dc = abs(Back.Order[j].energy - Front.Order[w].energy);
			if (dc >= MaxEDiff || !inCsI) C[j][w] = POISON;
			else C[j][w] = dc;
		}
	}

	// Step 4: Run Hungarian algorithm on above cost matrix
	// In the language of the Hungarian algorithm found in `algorithms.h`, the front strips
	// are the "workers" and the back strips are the "jobs". This function returns a vector
	// mapping worker indices to jobs. In other words, if you plug an index for Front.Order
	// into jobs[], you get out an index for Back.Order. If jobs[w] = -1, then the worker
	// of interest was never assigned a job. However, since I am feeding in a square matrix,
	// this should never be the case. The two possibilities are that either jobs[w] points
	// to a real back strip, or it points to a padded region with the cutoff value MaxEDiff.
	vector<int> f2b = hungarian(C);
/*
	// Debug step: Print out cost array and final front-back pairings
	ss << "Cost array (x = front, y = back):" << endl;
	for (const auto& v : C) {
		for (const auto& cost : v)
			ss << cost << "\t";
		ss << endl;
	}
	ss << "Hungarian algorithm results (front-back pairs + dummy elements):" << endl;
	for (size_t i = 0; i < S; i++)
		ss << f2b[i] << " ";
	ss << endl;
*/		
	// Now assign each of these solutions a CsI detector location
	vector<vector<pair<size_t, size_t>>> sil(NCsI, vector<pair<size_t, size_t>>()); // contains a lits of (front, back) silicon solutions for each Csi

	// Look at all the front/back solutions and see how many are on each CsI
	for (size_t w = 0; w < Front.Nstore; w++) {
		if (f2b[w] >= Back.Nstore || f2b[w] == -1) continue; // not a valid pair
		int ifront = Front.Order[w].strip;
		int iback = Back.Order[f2b[w]].strip;
		for (size_t icsi = 0; icsi < NCsI; icsi++) {
			if ((ifront < CsIFextents[icsi].first) || (ifront > CsIFextents[icsi].second) || (iback < CsIBextents[icsi].first) || (iback > CsIBextents[icsi].second)) continue;

			pair<size_t, size_t> p(w, f2b[w]);
			sil[icsi].push_back(p);
			break;
		}
	}

	// Make array of detect csi energies
	vector<double> energy(NCsI, -1.);
	vector<double> energyR(NCsI, -1.);
	vector<short> order(NCsI, -1);

	// Store the CsI raw energy info in an array that corresponds to the position it is in
	for (size_t i = 0; i < CsI.Nstore; i++) {
		if (CsI.Order[i].strip >= NCsI) {
			stringstream ss;
			ss << "[ERROR] CsI ID " << CsI.Order[i].strip
			   << " greater than NCsI-1=" << NCsI-1 << endl;
			cerr << ss.str();
			abort();
		}
		//else if (CsI.Order[i].time < -100 || CsI.Order[i].time > 100) continue;
		energy[CsI.Order[i].strip] = CsI.Order[i].energy;
		order[CsI.Order[i].strip] = i;
	}

	// Define the remove constants as some index that won't happen
	int removeFirst = -1;
	int removeSecond = -1;

	// Loop over CsI location
	for (size_t icsi = 0; icsi < NCsI; icsi++) {
		size_t multSi = sil[icsi].size();

		// No solution for this location, ignore
		if (multSi == 0 || order[icsi] < 0) continue;

		// FIXME DEE matching won't work until you have good zlines. Turn off at start
		// If more than 1 si solution for a single CsI, check if it falls in a zline
		// CsI can only fire once within readout
		// Needed for events with mixed E and CsI in the same quad
		// Can only accept one solution, don't allow it to accept both
		else if (multSi > 1) {

			if (PidCsI[icsi] == nullptr) continue; // ignore if zline files are absent

			for (size_t i = 0; i < multSi; i++) {
				int iif = sil[icsi][i].first;
				int iib = sil[icsi][i].second;

				// Do zline check
				int zCheck = PidCsI[icsi]->getPID(CsI.Order[order[icsi]].energyR, Front.Order[iif].energy);
				if (zCheck == 0) continue;
				else { // need to fill stuff here using the correct iif and iib Order indices
					if (order[icsi] < 0 || order[icsi] >= CsI.Nstore) {
						stringstream ss;
						ss << "[ERROR] Invalid CsI index order[icsi]=" << order[icsi] << endl;
						cerr << ss.str();
						abort();
					}
					Solution[Nsolution].energy = energy[icsi];
					Solution[Nsolution].energyR = CsI.Order[order[icsi]].energyR;
					Solution[Nsolution].energylow = energy[icsi];
					Solution[Nsolution].energylowR = CsI.Order[order[icsi]].energyR;
					Solution[Nsolution].denergy = Front.Order[iif].energy;
					Solution[Nsolution].denergylow = Front.Order[iif].energylow;
					Solution[Nsolution].denergyR = Front.Order[iif].energyR;
					Solution[Nsolution].benergy = Back.Order[iib].energy;
					Solution[Nsolution].benergylow = Back.Order[iib].energylow;
					Solution[Nsolution].benergyR = Back.Order[iib].energyR;
					Solution[Nsolution].qdc = CsI.Order[order[icsi]].qdc;

					Solution[Nsolution].ifront = Front.Order[iif].strip;
					Solution[Nsolution].iback = Back.Order[iib].strip;
					Solution[Nsolution].ide = -1;
					Solution[Nsolution].iCsI = icsi;
					Solution[Nsolution].itele = id;
					Solution[Nsolution].isSiCsI = true;
					Solution[Nsolution].CsITime = CsI.Order[order[icsi]].time;
					float timediff = CsI.Order[order[icsi]].time - Front.Order[iif].time;
					Solution[Nsolution].timediff = timediff;
					Solution[Nsolution].time = Front.Order[iif].time;
					Solution[Nsolution].btime = Back.Order[iib].time;
					Nsolution++;

					Front.Order[iif].CsIFlag = true;
					Back.Order[iib].CsIFlag = true;
					
					/*if (CsI.Nstore > 2) {
						cout << "Matched CsI, F, and B id 's : " << icsi << " " << Front.Order[ii].strip << " " << Back.Order[arrayB[ii]].strip << endl;
					}*/

					break; //break out of loop, accept only one solution per quad. Allows for small chance of losing dE solution
					//I could find a more elegant solution using strip matching for dE and base it on a best score
				}
			}
		}

		// CsI energy < 0 should not happen, but ignore just in case
		else if (energy[icsi] <= 0.) continue;
		else {
			int iif = sil[icsi][0].first;
			int iib = sil[icsi][0].second;
			if (order[icsi] < 0 || order[icsi] >= CsI.Nstore) {
				stringstream ss;
				ss << "[ERROR] Invalid CsI index order[icsi]=" << order[icsi] << endl;
				cerr << ss.str();
				abort();
			}
			Solution[Nsolution].energy = energy[icsi];
			Solution[Nsolution].energyR = CsI.Order[order[icsi]].energyR;
			Solution[Nsolution].energylow = energy[icsi];
			Solution[Nsolution].energylowR = CsI.Order[order[icsi]].energyR;
			Solution[Nsolution].denergy = Front.Order[iif].energy;
			Solution[Nsolution].denergylow = Front.Order[iif].energylow;
			Solution[Nsolution].denergyR = Front.Order[iif].energyR;
			Solution[Nsolution].benergy = Back.Order[iib].energy;
			Solution[Nsolution].benergylow = Back.Order[iib].energylow;
			Solution[Nsolution].benergyR = Back.Order[iib].energyR;
			Solution[Nsolution].qdc = CsI.Order[order[icsi]].qdc;

			Solution[Nsolution].ifront = Front.Order[iif].strip;
			Solution[Nsolution].iback = Back.Order[iib].strip;
			Solution[Nsolution].ide = -1;
			Solution[Nsolution].iCsI = icsi;
			Solution[Nsolution].itele = id;
			Solution[Nsolution].isSiCsI = true;
			Solution[Nsolution].CsITime = CsI.Order[order[icsi]].time;
			float timediff = CsI.Order[order[icsi]].time - Front.Order[iif].time;
			Solution[Nsolution].timediff = timediff;
			Solution[Nsolution].time = Front.Order[iif].time;
			Solution[Nsolution].btime = Back.Order[iib].time;
			Nsolution++;

			Front.Order[iif].CsIFlag = true;
			Back.Order[iib].CsIFlag = true;
			
			/*if (CsI.Nstore > 2) {
				cout << "Matched CsI, F, and B id 's : " << icsi << " " << Front.Order[ii].strip << " " << Back.Order[arrayB[ii]].strip << endl;
			}*/
		}
	}
	//if (CsI.Nstore > 2 && Nsolution > 0) cout << "Nsolution " << Nsolution << endl;
	return Nsolution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Calculates the x-y position and angles in the array in cm
void telescope::position(int isol) {
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

// Calculates the x-y position and angles in the array in cm
pair<double, double> telescope::simplePosition() {
	double Xpos,Ypos;

	int ifront = Front.Order[0].strip;
	int iback = Back.Order[0].strip;
	if (id == 0) 
	{
		Xpos = Xcenter + (((double)iback+Ran->Rndm())/32.-0.5)*SiWidth;
		Ypos = Ycenter + (((double)ifront+Ran->Rndm())/32.-0.5)*SiWidth;
	}
	else if (id == 1)
	{
		Xpos = Xcenter + (((double)ifront+Ran->Rndm())/32.-0.5)*SiWidth;
		Ypos = Ycenter + (0.5-((double)iback+Ran->Rndm())/32.)*SiWidth;
	}
	else if (id == 2)
	{
		Xpos = Xcenter + (0.5-((double)iback+Ran->Rndm())/32.)*SiWidth;
		Ypos = Ycenter + (0.5-((double)ifront+Ran->Rndm())/32.)*SiWidth;
	}
	else if (id == 3)
	{
		Xpos = Xcenter + (0.5-((double)ifront+Ran->Rndm())/32.)*SiWidth;
		Ypos = Ycenter + (((double)iback+Ran->Rndm())/32.-0.5)*SiWidth;
	}

	//	Xpos += .3;

	pair<double, double> p(Xpos, Ypos);
	return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Calculates the x-y position in the array in cm
void telescope::positionC(int isol) {
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

bool telescope::isCenter(size_t ifront, size_t iback) {
	for (size_t i = 0; i < NCsI; i++) {
		if ((ifront == CsIFmids[i] || ifront == (CsIFmids[i] + 1)) && (iback == CsIBmids[i] || iback == (CsIBmids[i] + 1))) return true;
	}
	return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Converts equilivant proton energy to energy for a given isotope
// i.e. Z and A dependence of CsI light output
double telescope::light2energy(size_t Z, size_t A, size_t tel, size_t id, double energy) {
	if (!hasCsI) throw invalid_argument(string(BOLDRED) + string("ERROR: cannot invoke per-particle CsI calibrations for telescope with no CsI") + string(RESET));

	pair<size_t, size_t> ZA(Z, A);
	if (ZA == sz_pair(1, 1))
		return energy;
	else if (ZA == sz_pair(1, 2)) // deuterons
		energy = calCsI_d->getEnergy(tel, id, energy);
	else if (ZA == sz_pair(1, 3)) // tritons
		energy = calCsI_t->getEnergy(tel, id, energy);
	else if (ZA == sz_pair(2, 3)) // 3He, same calibration as alphas
		energy = calCsI_Alpha->getEnergy(tel, id, energy);
	else if (ZA == sz_pair(2, 4)) // alphas
		energy = calCsI_Alpha->getEnergy(tel, id, energy);
	else throw invalid_argument(string(BOLDRED) + string("ERROR: found no calib for Z = ") + to_string(Z) + string(" A = ") + to_string(A) + string(" tel = ") + to_string(tel) + string(" Csi ID = ") + to_string(id) + string(RESET));
	return energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



