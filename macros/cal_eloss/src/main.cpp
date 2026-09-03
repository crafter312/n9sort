#include <iostream>
#include <vector>
#include <iomanip>
#include <catima/catima.h>

#include "constants.h"

using namespace std;

int main() {

	/******** PROJECTILE ********/
	// Here, I use my mass excess lookup table to automatically calculate the isotopic
	// mass in AMU as a function of A and Z, so that all one needs to provide is A, Z,
	// and the total kinetic energy in MeV in order to define the projectile.

	const size_t A = 4;
	const size_t Z = 2;
	catima::Projectile p(A, Z);
	double T_MeV   = 280.; // total kinetic energy in MeV
	double mass_u  = (double)A + (mass_lookup.at({Z, A}) / catima::atomic_mass_unit);
	double T_per_u = T_MeV / mass_u; // total kinetic energy in MeV/u
	p.T = T_per_u;

	/******** MATERIALS ********/

	// 1/4 inch degrader (average measured real thickness in mm)
	catima::Material degrader_thin = catima::get_material(13); // Al with elemental atomic weight and standard density
	degrader_thin.thickness_cm(6.3345 * 0.1); 

	// 1/2 inch degrader (average measured real thickness in mm)
	catima::Material degrader_thick = catima::get_material(13); // Al with elemental atomic weight and standard density
	degrader_thick.thickness_cm(12.708 * 0.1);

	// 3.8 mm aluminum covers (per telescope with average measured real thickness in mm)
	const vector<double> cover_thicknesses = {3.761, 3.7505, 3.7755, 3.755};
	vector<catima::Material> covers;
	for (size_t i = 0; i < 4; i++) {
		catima::Material cover = catima::get_material(13); // Al with elemental atomic weight and standard density
		cover.thickness_cm(cover_thicknesses[i] * 0.1);
		covers.push_back(cover);
	}

	// 1500 um silicon detectors (per telescope with real thicknesses in um from Micron Semiconductors)
	const vector<double> si_thicknesses = {1536., 1500., 1537., 1460.};
	vector<catima::Material> silicons;
	for (size_t i = 0; i < 4; i++) {
		catima::Material silicon = catima::get_material(14); // Si with elemental atomic weight and standard density
		silicon.thickness_cm(si_thicknesses[i] * 0.0001);
		silicons.push_back(silicon);
	}

	/******** LAYERS ********/
	// CATima facilitates energy loss calculations for multiple materials in succession
	// by allowing you to define a `Layers` object of multiple materials and pass them
	// all into the calculation at once. It then returns a set of results for each layer
	// in the stack, as well as for the whole stack.
	
	// For now, I'll just pick one set of simple materials, ignoring effects of angle
	// and per-telescope differences.
	catima::Layers testMaterialSequence;
	testMaterialSequence.add(degrader_thin);
	testMaterialSequence.add(covers[0]);
	testMaterialSequence.add(silicons[0]);

	/******** CALCULATE ********/
	catima::MultiResult multiResult = catima::calculate(p, testMaterialSequence);

	/******** OUTPUT ********/

	// Output projectile information
	cout << fixed << setprecision(3);
	cout << "=== Energy loss calculation projectile ===" << endl;
	cout << "  A = " << A << endl;
	cout << "  Z = " << Z << endl;
	cout << "  Energy = " << T_MeV << " MeV" << endl;
	cout << "         = " << T_per_u << " MeV/u" << endl;

	// Output results per-layer
	cout << endl;
	cout << "=== Per-layer energy loss results ===" << endl;
	for (size_t i = 0; i < multiResult.results.size(); i++) {
		cout << "Layer " << i << ":" << endl;
		for (size_t j = 0; j < testMaterialSequence[i].ncomponents(); j++) {
			catima::Target atom = testMaterialSequence[i].get_element(j);
			cout << "  Atom " << j << ": A = " << atom.A << ", Z = " << atom.Z << ", weight fraction = " << atom.stn << endl;
		}
		cout << "  Thickness: " << testMaterialSequence[i].thickness_cm() << " cm" << endl;

		double Ein = multiResult.results[i].Ein*mass_u;
		double Eout = multiResult.results[i].Eout*mass_u;
		cout << "  Ein: " << Ein << " MeV" << endl;
		cout << "  Eout: " << Eout << " MeV" << endl;
		cout << "  Eloss: " << Ein - Eout << " MeV" << endl;
	}

	return 0;
}
