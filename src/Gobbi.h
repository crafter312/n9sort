#ifndef Gobbi_H
#define Gobbi_H

/* Split off from `Det` (previously `Gobbi`) by Henry Webb (h.s.webb@wustl.edu)
 * 4 May 2026 to match other, newer versions of the Wash U radiochemistry sort
 * code than the code I started with.
 */

#include "calibrate.h"
#include "correl2.h"
#include "histo.h"
#include "Input.h"
#include "silicon.h"
#include "solution.h"

class Gobbi {

public:
	Gobbi(Input& in, histo& hist, SortConfig& config);
	~Gobbi();

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
	const Input::TDCInput& input_tdc;

	silicon* Silicon[4];

	calibrate* FrontEcal;
	calibrate* BackEcal;
	calibrate* DeltaEcal;
	calibrate* FrontTimecal;
	calibrate* BackTimecal;
	calibrate* DeltaTimecal;

	// Inputs are telescope number, channel number, and hit index for below functions
	void addFrontHit(size_t tel, size_t ch, size_t i);
	void addBackHit(size_t tel, size_t ch, size_t i);
	void addDeltaHit(size_t tel, size_t ch, size_t i);
};

#endif
