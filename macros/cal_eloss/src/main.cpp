#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

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

	// 1/4 inch degrader (average measured real thickness in mm), half
	catima::Material degrader_thin = catima::get_material(13); // Al with elemental atomic weight and standard density
	const double degrader_thin_half_cm = 6.3345 * 0.1 * 0.5;
	degrader_thin.thickness_cm(degrader_thin_half_cm); 

	// 1/2 inch degrader (average measured real thickness in mm), half
	catima::Material degrader_thick = catima::get_material(13); // Al with elemental atomic weight and standard density
	const double degrader_thick_half_cm = 12.708 * 0.1 * 0.5;
	degrader_thick.thickness_cm(degrader_thick_half_cm);

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

	/******** ANGLES ********/
	// Here, I calculate the polar angle corresponding to the center of each CsI crystal
	// in Gobbi. These values are rotationally symmetric around the beam axis, so I only
	// need to calculate seven for one telescope. These are used for elastic scattering
	// and angle-dependent effective thickness calculations.

	// Gobbi geometric constants (all distances in cm)
	const double silicon_frame = 7.237;    // exterior silicon frame dimension
	const double inner_hole    = 2.5;      // interior Gobbi hole dimension
	const double tele_center_x = (.5 * silicon_frame) + (.5 * inner_hole);
	const double tele_center_y = (.5 * silicon_frame) - (.5 * inner_hole);
	const double large_csi     = 3.42011;  // width of large CsI crystal
	const double small_csi     = 1.642827; // width of small CsI crystal
	const double gobbi_dist    = 54.7;     // distance between center of target and front face of silicons
	const double silicon_thick = 0.65;     // approximate distance from front face of silicons to front face of CsI crystals
	const double total_dist    = gobbi_dist + silicon_thick;

	// Calculate per-crystal values
	vector<double> front_center_x = {
		tele_center_x - (1.5 * small_csi), tele_center_x - (.5 * small_csi), tele_center_x - (.5 * small_csi), tele_center_x - (1.5 * small_csi),
		tele_center_x + (.5 * large_csi), tele_center_x + (.5 * large_csi), tele_center_x - (.5 * large_csi) 
	};
	vector<double> front_center_y = {
		tele_center_y - (1.5 * small_csi), tele_center_y - (1.5 * small_csi), tele_center_y - (.5 * small_csi), tele_center_y - (.5 * small_csi),
		tele_center_y - (.5 * large_csi), tele_center_y + (.5 * large_csi), tele_center_y + (.5 * large_csi) 
	};
	vector<double> polar_angle_rad(7, 0.);
	vector<double> polar_angle_deg(7, 0.);
	for (size_t i = 0; i < 7; i++) {
		polar_angle_rad[i] = atan(sqrt((front_center_x[i]*front_center_x[i]) + (front_center_y[i]*front_center_y[i])) / total_dist);
		polar_angle_deg[i] = polar_angle_rad[i] * rad_to_deg;
	}

	/******** ENERGY LOSS ********/

	// First, calculate loss through first half of target. This is the same for all
	// CsI crystals in all telescopes.
	catima::Result degraderResult = catima::calculate(p, degrader_thin);
	double Eout = degraderResult.Eout*mass_u;

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
		Eout = multiResult.results[i].Eout*mass_u;
		cout << "  Ein: " << Ein << " MeV" << endl;
		cout << "  Eout: " << Eout << " MeV" << endl;
		cout << "  Eloss: " << Ein - Eout << " MeV" << endl;
	}

	return 0;
}
