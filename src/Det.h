#ifndef Det_H
#define Det_H

/* Nicolas Dronchi 2022_04_04
 * Class written to handle all specifics of the Gobbi array
 * such as communicating with HINP, calibrations, checking
 * for charge sharing in neighbor calculating geometry.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu), August 2025
 * Replaces use of HINP class for unpacking with Input class
 * for reading values from a SpecTcl-generated ROOT file.
 * This essentially offloads the work of unpacking to SpecTcl,
 * leaving this code to only do the analysis work.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu) and Johnathan
 * Phillips (j.s.phillips@wustl.edu) March 2026 for experiment
 * at TAMU Cyclotron Institute.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu) April 2026 for
 * 9N FRIB experiment. Renamed to `Det` to match other sort
 * code versions, Gobbi functionality moved to dedicated
 * Gobbi class(s).
 */

#include "calibrate.h"
#include "correl2.h"
#include "OldGobbi.h"
#include "histo.h"
#include "Input.h"
#include "silicon.h"
#include "solution.h"
#include "SortConfig.h"
#include "wood.h"

#include <iostream>
#include <memory>
#include <string>

class Det {

public:
	Det(Input& in, histo& hist, SortConfig& config, size_t run);
	~Det();

	void analyze();

	float getEnergy(int board, int chan, int Ehigh);

	void corr_4He();
	void corr_5He();
	void corr_6He();
	void corr_5Li();
	void corr_6Li();
	void corr_7Li();
	void corr_6Be();
	void corr_7Be();
	void corr_8Be();
	void corr_9B();
	void corr_8C();
	void corr_9N();

	// Each decay channel should have its own wood object instance in
	// the form of a unique pointer so that initialization can be done
	// in the constructor.
	std::unique_ptr<wood> He4_pt;
	std::unique_ptr<wood> He4_dd;
	std::unique_ptr<wood> He5_dt;
	std::unique_ptr<wood> He6_tt;
	std::unique_ptr<wood> Li5_pa;
	std::unique_ptr<wood> Li5_dHe3;
	std::unique_ptr<wood> Li6_da;
	std::unique_ptr<wood> Li7_ta_bad;
	std::unique_ptr<wood> Li7_pHe6;
	std::unique_ptr<wood> Li7_ta;
	std::unique_ptr<wood> Be6_ppa;
	std::unique_ptr<wood> Be7_He3a;
	std::unique_ptr<wood> Be7_pLi6;
	std::unique_ptr<wood> Be8_aa;
	std::unique_ptr<wood> Be8_pLi7;
	std::unique_ptr<wood> Be8_pta;
	std::unique_ptr<wood> Li7_ta_fake; // missing the p
	std::unique_ptr<wood> B9_paa;
	std::unique_ptr<wood> C8_4pa;
	std::unique_ptr<wood> N9_5pa;

	histo& Histo;
	OldGobbi gobbi;
	correl2 Correl;

	double Targetdist;      // cm
	double TargetThickness; // mg/cm^2
	size_t hinpboards;      // total number of HINP boards used
	size_t hinpchans;       // number of channels per HINP board (should always be 32)
	size_t runnum;
	
	// Record particle combinations, start with most important
	size_t a_p{0};

private:
	const Input::GobbiInput& input;
	const Input::QDCInput& input_qdc;
	const Input::TDCInput& input_tdc;

};


#endif
