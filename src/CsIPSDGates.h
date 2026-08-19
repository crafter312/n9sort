#ifndef CsIPSDGates_H
#define CsIPSDGates_H

/**
 * Class written by Henry Webb (h.s.webb@wustl.edu) 19 August 2026
 * to handle the input, storage, and application of graphical 2D
 * cuts on CsI PSD plots (QDC vs. ADC) for general cleanup. This
 * class does not store gates for individual particle PSD using the
 * CsI crystals. It only stores general gates for data cleanup (i.e.
 * one gate per crystal).
 */

#include <unordered_map>
#include <utility>

#include <TCutG.h>

#include "constants.h"
#include "SortConfig.h"

class CsIPSDGates {

public:
	CsIPSDGates(SortConfig& config);
	~CsIPSDGates();

	bool IsParticle(size_t tele, size_t iCsI, double ER, double Q);

private:
	std::unordered_map<std::pair<size_t, size_t>, TCutG, pair_hash> crystal_gate_map; // stores one TCutG per CsI crystal, with (telescope, crystal ID) being the key

};

#endif
