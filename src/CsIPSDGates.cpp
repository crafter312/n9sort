/**
 * Class written by Henry Webb (h.s.webb@wustl.edu) 19 August 2026
 * to handle the input, storage, and application of graphical 2D
 * cuts on CsI PSD plots (QDC vs. ADC) for general cleanup. This
 * class does not store gates for individual particle PSD using the
 * CsI crystals. It only stores general gates for data cleanup (i.e.
 * one gate per crystal).
 */

#include "CsIPSDGates.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CsIPSDGates::CsIPSDGates(SortConfig& config) {

	// Open file with CsI PSD gates
	string ifname = config.GetConfigDir() + config.GetCsIPSDGateFile();
	ifstream ifile(ifname);
	if (!ifile.is_open()) throw invalid_argument(string(BOLDRED) + string("ERROR: Could not open CsI PSD gate file ") + ifname + string(RESET));

	// Read file with CsI PSD gates
	size_t ngates, tele, iCsI, npoints;
	ifile >> ngates;
	for (size_t i = 0; i < ngates; i++) {

		// Initialize empty TCutG with correct key
		ifile >> tele >> iCsI;
		pair<size_t, size_t> key{tele, iCsI};
		crystal_gate_map.try_emplace(key);

		// Now retrieve empty TCutG and add points from input list
		ifile >> npoints;
		double ER, Q;
		TCutG& gate = crystal_gate_map[key];
		for (size_t j = 0; j < npoints; j++) {
			ifile >> ER >> Q;
			gate.AddPoint(ER, Q);
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CsIPSDGates::~CsIPSDGates() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool CsIPSDGates::IsParticle(size_t tele, size_t iCsI, double ER, double Q) {
	auto mapEntry = crystal_gate_map.find({tele, iCsI});
	if (mapEntry == crystal_gate_map.end()) {
		cout << "WARNING: No PSD cleanup gate for telescope " << tele << ", CsI crystal " << iCsI << endl;
		return true;
	}
	return mapEntry->second.IsInside(ER, Q);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



