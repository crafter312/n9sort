#include "losses.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * constructor
 * \param Zmax0 is the maximum Z nucleus whose loss table should be stored in this class
 * \param path is the directory relative to the program build location containing files with energy loss tables
 * \param mat is a string denoting which material the energy loss is required for, and therefore which files get loaded
 */
CLosses::CLosses(int Zmax0, string path, string mat) : material(mat) {
	Zmax = Zmax0;
#if defined(eloss_mode) && eloss_mode == 0
	loss = new CLoss3*[Zmax];
	path += "SRIM/";
#elif defined(eloss_mode) && eloss_mode == 1
	loss = new CLoss4*[Zmax];
	path += "ATIMA/";
#else
#error "compile-time variable 'eloss_mode' is undefined or set to an unsupported value!"
#endif
	const string suffix = "_" + mat + ".loss";
	string filename;
	for (size_t iZ = 1; iZ <= Zmax; iZ++) {
		switch (iZ) {
			case 1:
				filename = path + "Hydrogen" + suffix;
				break;
			case 2:
				filename = path + "Helium" + suffix;
				break;
			case 3:
				filename = path + "Lithium" + suffix;
				break;
			case 4:
				filename = path + "Beryllium" + suffix;
				break;
			case 5:
				filename = path + "Boron" + suffix;
				break;
			case 6:
				filename = path + "Carbon" + suffix;
				break;
			case 7:
				filename = path + "Nitrogen" + suffix;
				break;
			case 8:
				filename = path + "Oxygen" + suffix;
				break;
			case 9:
				filename = path + "Fluorine" + suffix;
				break;
			case 10:
				filename = path + "Neon" + suffix;
				break;
			default:
				throw_error_arg("Loss case for Z = " + to_string(iZ) + " not supported!");
		}
		
#ifdef ENABLE_DEBUG
		cout << iZ << endl;
		cout << filename << endl;
#endif

#if defined(eloss_mode) && eloss_mode == 0
		loss[iZ - 1] = new CLoss3(filename);
#elif defined(eloss_mode) && eloss_mode == 1
		loss[iZ - 1] = new CLoss4(filename);
#else
#error "compile-time variable 'eloss_mode' is undefined or set to an unsupported value!"
#endif
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CLosses::~CLosses() {
	for (size_t i = 0; i < Zmax; i++) delete loss[i];
	delete loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double CLosses::getEin(double energy, double thick, size_t Z, double A) {
	if (Z > Zmax || Z == 0) throw invalid_argument(string(BOLDRED) + string("No loss info for Z = ") + to_string(Z) + string(RESET));
	return loss[Z - 1]->getEin(energy, thick, A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double CLosses::getEout(double energy, double thick, size_t Z, double A) {
	if (Z > Zmax || Z == 0) throw invalid_argument(string(BOLDRED) + string("No loss info for Z = ") + to_string(Z) + string(RESET));
	
	return loss[Z - 1]->getEout(energy, thick, A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



