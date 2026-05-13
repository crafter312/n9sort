#include "losses.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * constructor
 * \param Zmax0 is the maximum Z nucleus whose loss table should be stored in this class
 * \param config is a reference to the sort code's config class, used to retrieve paths of files containing energy loss tables
 */
CLosses::CLosses(int Zmax0, SortConfig& config) {
	Zmax = Zmax0;
	string path = config.GetLossDir();
	string suffix = "_" + config.GetTargetSuffix() + ".loss";
	string filename;
	loss = new CLoss2*[Zmax];
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
				throw invalid_argument(string(BOLDRED) + string("loss case for Z = ") + to_string(iZ) + string( " not implemented") + string(RESET));
		}
		
#ifdef ENABLE_DEBUG
		cout << iZ << endl;
		cout << filename << endl;
#endif

		loss[iZ - 1] = new CLoss2(filename);
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



