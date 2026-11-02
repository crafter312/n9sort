#include "loss4.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Constructor for CLoss4 class.
 * @param filename absolute or relative path of ATIMA spline file to be read in
 */
CLoss4::CLoss4(const std::string& filename) {
	ifstream file(filename);
	if (!file.is_open()) throw_error("[CLoss4] ERROR: Unable to open ATIMA file: " + filename);

	string line;

	// Skip header (first four lines)
	for (int i = 0; i < 4; ++i) {
		if (!getline(file, line)) throw_error("[CLoss4] ERROR: Failed skipping header for " + filename);
	}

	/******** READ SPLINE ********/
	// NOTE: First spline is range. Third spline is probably dE/dx. Second and fourth
	// splines are likely range straggling and angular straggling, but I'm not sure.
	// Here, we only read the first spline, since we're only interested in basic energy
	// loss calculations for now.

	// Read n (number of basis splines) and k (order)
	if (!std::getline(file, line)) throw_error("[CLoss4] ERROR: Failed to read n and k for range spline from " + filename);
	stringstream ss_nk(line);
	if (!(ss_nk >> m_n >> m_k)) throw_error("[CLoss4] ERROR: Invalid format for range spline n and k in " + filename);

	size_t num_knots = m_n + m_k;

	// Read next line of n + k knots (control points)
	if (!getline(file, line)) throw_error("[CLoss4] ERROR: Failed to read range spline knots from " + filename);
	
	gsl_vector* knots_tmp = gsl_vector_alloc(num_knots);
	stringstream ss_knots(line);
	double knot_val = 0.0;
	size_t knot_idx = 0;

	while (ss_knots >> knot_val) {
		if (knot_idx >= num_knots) {
			gsl_vector_free(knots_tmp);
			throw_error("[CLoss4] ERROR: More range spline knots than expected in " + filename);
		}
		gsl_vector_set(knots_tmp, knot_idx++, knot_val); // postfix ++ operator returns value before incrementing
	}

	// Read next line of n bspline coefficients
	if (!getline(file, line)) {
		gsl_vector_free(knots_tmp);
		throw_error("[CLoss4] ERROR: Failed to read range spline coefficients from " + filename);
	}

	m_coeffs = gsl_vector_alloc(m_n);
	stringstream ss_coeffs(line);
	double coeff_val = 0.0;
	size_t coeff_idx = 0;

	while (ss_coeffs >> coeff_val) {
		if (coeff_idx >= m_n) {
			gsl_vector_free(knots_tmp);
			gsl_vector_free(m_coeffs);
			m_coeffs = nullptr;
			throw_error("[CLoss4] ERROR: More range spline coefficients than expected in " + filename);
		}
		gsl_vector_set(m_coeffs, coeff_idx++, coeff_val); // postfix ++ operator returns value before incrementing
	}

	// Initialize GSL bspline workspace using the temporary knot vector
	m_bw = gsl_bspline_alloc_ncontrol(m_k, m_n);
	gsl_bspline_init(knots_tmp, m_bw);

	// Free temporary knot vector memory
	gsl_vector_free(knots_tmp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CLoss4::~CLoss4() {
	if (m_bw) {
		gsl_bspline_free(m_bw);
		m_bw = nullptr;
	}
	if (m_coeffs) {
		gsl_vector_free(m_coeffs);
		m_coeffs = nullptr;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Calculates the residual energy of a particle after passing through a material,
 * with the particle and material being fixed for any particular instance of
 * CLoss4. The A is only to convert from total kinetic energy to energy per AMU.
 * @param energy the initial total kinetic energy of the particle in MeV
 * @param thick the thickness of the material in mg/cm^2
 * @param A the mass of the particle in AMU, only used to convert the energy
 * @return the total kinetic energy of the particle after passing through the material
 * @throws a runtime_error if the Newton-Raphson iteration fails to converge
 *
 * NOTE: The header of the ATIMA spline file does contain the mass of the particle
 *       in MeV/u, but it's probably best to pass A into these functions to make
 *       sure that it is the same as in the rest of the code. I do not have control
 *       over or certainty of where ATIMA's mass values come from.
 */
double CLoss4::getEout(double energy, double thick, double A) {
	double epa = energy / A;
	size_t counter = 0;
	double range;
	double dedx;
	double r, rtemp;

	gsl_bspline_calc(epa, m_coeffs, &range, m_bw);
	gsl_bspline_calc_deriv(epa, m_coeffs, 1, &dedx, m_bw);
	dedx = 1. / dedx;
	if (range <= thick) return 0.;

	epa -= thick * dedx;
	while (true) {
		gsl_bspline_calc(epa, m_coeffs, &rtemp, m_bw);
		r = range - rtemp - thick;
		if (fabs(r) < Eout_th_epsilon) break;
		epa += r * dedx;
		if (epa < Ezero) return 0.;
		gsl_bspline_calc_deriv(epa, m_coeffs, 1, &dedx, m_bw);
		dedx = 1. / dedx;
		counter++;
		if (counter > Eloss_max_iter) {
			throw_error("[CLoss4] ERROR: Newton-Raphson iteration failed to converge for energy = "
				+ to_string(energy) + " MeV, thick = "
				+ to_string(thick) + " mg/cm^2, A = "
				+ to_string(A) + " u"
			);
		}
	}
	return epa * A;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Calculates the initial energy of a particle before passing through a material,
 * with the particle and material being fixed for any particular instance of
 * CLoss4. The A is only to convert from total kinetic energy to energy per AMU.
 * @param energy the residual total kinetic energy of the particle in MeV
 * @param thick the thickness of the material in mg/cm^2
 * @param A the mass of the particle in AMU, only used to convert the energy
 * @return the total kinetic energy of the particle before passing through the material
 *
 * NOTE: The header of the ATIMA spline file does contain the mass of the particle
 *       in MeV/u, but it's probably best to pass A into these functions to make
 *       sure that it is the same as in the rest of the code. I do not have control
 *       over or certainty of where ATIMA's mass values come from.
 */
double CLoss4::getEin(double energy, double thick, double A) {
	double epa = energy / A;
	size_t counter = 0;
	double range;
	double dedx;
	double r, rtemp;

	gsl_bspline_calc(epa, m_coeffs, &range, m_bw);
	gsl_bspline_calc_deriv(epa, m_coeffs, 1, &dedx, m_bw);
	dedx = 1. / dedx;
	if (range <= thick) return 0.;

	epa += thick * dedx;
	while (true) {
		gsl_bspline_calc(epa, m_coeffs, &rtemp, m_bw);
		r = range - rtemp - thick;
		if (fabs(r) < Eout_th_epsilon) break;
		epa -= r * dedx;
		if (epa < Ezero) return 0.;
		gsl_bspline_calc_deriv(epa, m_coeffs, 1, &dedx, m_bw);
		dedx = 1. / dedx;
		counter++;
		if (counter > Eloss_max_iter) {
			throw_error("[CLoss4] ERROR: Newton-Raphson iteration failed to converge for energy = "
				+ to_string(energy) + " MeV, thick = "
				+ to_string(thick) + " mg/cm^2, A = "
				+ to_string(A) + " u"
			);
		}
	}
	return epa * A;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



