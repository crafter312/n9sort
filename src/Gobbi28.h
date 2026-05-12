#ifndef Gobbi28_H
#define Gobbi28_H

/* New Gobbi class for the Gobbi 28 configuration by Henry Webb (h.s.webb@wustl.edu).
 * Created 7 May 2026 as mix of `OldGobbi` class from this code and `gobbi` class
 * used by Johnathan Phillips (j.s.phillips@wustl.edu) for 22Si FRIB experiment.
 */

#include "calibrate.h"
#include "correl2.h"
#include "histo.h"
#include "Input.h"
#include "silicon.h"
#include "solution.h"

class Gobbi28 {

public:
	Gobbi28(Input& in, histo& hist, SortConfig& config);
	~Gobbi28();

	void analyze();
	size_t loadSolutions(correl2& Correl);
	solution* getNextEmptySolution(solution* sol);
	float getEin(float, float, int, float);

private:
	double Targetdist;
	double TargetThickness;
	size_t hinpboards; // total number of HINP boards used
	size_t hinpchans;  // number of channels per HINP board (should always be 32)

	histo& Histo;

	const Input::GobbiInput& input;
	const Input::ADCInput& input_adc;
	const Input::QDCInput& input_qdc;
	const Input::TDCInput& input_tdc;

	silicon* Silicon[4];

	calibrate* FrontEcal;
	calibrate* BackEcal;
	calibrate* CsIEcal;
	calibrate* FrontTimecal;
	calibrate* BackTimecal;
	calibrate* CsITimecal;

	// Inputs are telescope number, channel number, and hit index for below functions
	void addFrontHit(size_t tel, size_t ch, size_t i);
	void addBackHit(size_t tel, size_t ch, size_t i);
};

#endif
