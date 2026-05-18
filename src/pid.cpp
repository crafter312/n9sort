/* Class to detemine PID (Particle IDentification) from E-DE map via
 * banana gates stored for each particle type.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu), September 2025.
 * Mass values no longer hardcoded, now retrieved from std::unordered_map. 
 */

#include "pid.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

pid::pid(string file, SortConfig& config) : Z(-1), A(-1) {

	// Open zline file
	string name = config.GetPIDDir() + file + ".zline";
	ifstream ifile(name.c_str());
	if (!ifile.is_open()) throw invalid_argument(string(BOLDRED) + string("ERROR: Could not open zline file ") + name + string(RESET));

	// Read zline file
	ifile >> nlines;
	par = new ZApar*[nlines];
	for (size_t i = 0; i < nlines; i++)
		par[i] = new ZApar(ifile);
	ifile.close();
	ifile.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

pid::~pid() {
	for (size_t i = 0; i < nlines; i++) delete par[i];
	delete[] par;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/* Returns true if particle is in a banana gate, false otherwise. The parameters
 * Z and A are loaded with the detected particle's values.
 * \param x energy of particle
 * \param y energy loss of particle
 */
bool pid::getPID(float x, float y) {
	Z = 0;
	A = 0;
	for (size_t i = 0; i < nlines; i++) {
		if (par[i]->inBanana(x,y)) {
			Z = par[i]->Z;
			A = par[i]->A;
			mass = getMass(Z,A);
			return true;
		}
	}
	return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double pid::getMass(size_t iZ, size_t iA) {
	auto mapEntry = Mass_lookup.find({iZ, iA});
	if (mapEntry == Mass_lookup.end()) {
		cout << "No mass info for Z = "<< iZ << " A = " << iA << endl;
		abort();
		return 0;
	}
	return mapEntry->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



