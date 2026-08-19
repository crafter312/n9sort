#ifndef Gobbi28_H
#define Gobbi28_H

/* New Gobbi class for the Gobbi 28 configuration by Henry Webb (h.s.webb@wustl.edu).
 * Created 7 May 2026 as mix of `OldGobbi` class from this code and `gobbi` class
 * used by Johnathan Phillips (j.s.phillips@wustl.edu) for 22Si FRIB experiment.
 */

#include <fstream>
#include <unordered_map>

#include <TRandom.h>

#include "calibrate.h"
#include "correl2.h"
#include "CsIPSDGates.h"
#include "histo.h"
#include "Input.h"
#include "solution.h"
#include "telescope.h"

class Gobbi28 {

public:
	Gobbi28(Input& in, histo& hist, SortConfig& config);
	~Gobbi28();

	void analyze();
	size_t loadSolutions(correl2& Correl);
	solution* getNextEmptySolution(solution* sol);
	float getEin(float, float, int, float);
	
	// Counter getters
	size_t GetNsimpleECsI() const { return NsimpleECsI; }
	size_t GetNmultiECsI() const { return NmultiECsI; }
	size_t GetPidSkipped() const { return pidSkipped; }

private:
	double Targetdist;
	double TargetThickness;
	size_t hinpboards;    // total number of HINP boards used
	size_t hinpchans;     // number of channels per HINP board (should always be 32)
	size_t NCsI;          // number of CsI crystals per telescope, determined by telescope class when reading in CsI extents file
	size_t maxadcchan{0}; // the largest adc channel that the CsI crystals go up to
	size_t tdcstart{16};  // defines which tdc channel the CsI crystals start at
	size_t nTelCsIs;      // number of CsI crystals per telescope, determined by Gobbb28 class when reading in CsI mapping file

	histo& Histo;

	const Input::GobbiInput& input;
	const Input::ADCQDCInput& input_adc;
	const Input::ADCQDCInput& input_qdc;
	const Input::TDCInput& input_tdc;

	telescope* Telescope[4];

	calibrate* FrontEcal;
	calibrate* BackEcal;
	calibrate* CsIEcal;
	calibrate* FrontTimecal;
	calibrate* BackTimecal;
	calibrate* CsITimecal;
	
	// Map of CsI ADC/QDC/TDC channel to telescope and id (per-telescope CsI channel)
	std::unordered_map<size_t, size_t> telCsImap;
	std::unordered_map<size_t, size_t> idCsImap;

	// Class for managing CsI PSD gates for general data cleanup
	CsIPSDGates csiGates;
	
	// Counters
	size_t NsimpleECsI{0}; // # simple events with one each of front, back, and CsI in a telescope
	size_t NmultiECsI{0};  // # more complex events with more than one of one of front, back, and CsI in a telescope
	size_t pidSkipped{0};  // # solutions invalidated during Eloss calculations due to errors thrown by Eloss calculations
	
	TRandom* Ran;

	//ofstream testOut{"../RootFiles/TestDataOut_20Entries.txt"};

	// Inputs are telescope number, channel number, and hit index for below functions
	void addFrontHit(size_t tel, size_t ch, size_t i);
	void addBackHit(size_t tel, size_t ch, size_t i);
	void addCsIHits(); // matches CsI data from ADC, QDC, and TDC by channel number
};

#endif
