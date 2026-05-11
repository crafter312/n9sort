// Code to analyze data from Gobbi 28 Si-CsI array, including CsI PSD
// Originally written by Nicolas Dronchi, 2020
// Heavily modified by Henry Webb (h.s.webb@wustl.edu), August 2025
//   Now skips unpacking, reads values from SpecTcl-generated ROOT file
//   (i.e. SpecTcl now does the unpacking). Uses TNLIB TexNeut analysis
//   library written by Alex Alafa.
// Heavily modified by Henry Webb again, April 2026
//   All TexNeut analysis pieces removed. Gobbi class modified to handle
//   Si-CsI array with subdivided inner quadrants, along with PSD
//   analysis for the CsI and extra bits for 9N. `wood` class added
//   for per-reconstruction ROOT TTree output (significantly rewritten).

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

#include "constants.h"
#include "Det.h"
#include "histo.h"
#include "Input.h"
#include "SortConfig.h"

using namespace std;

int main() {

	// Capture the start time
	auto start = chrono::high_resolution_clock::now();

#ifndef rel
	cout << "WARNING: compile time variable `rel` not defined, using Newtonian mechanics" << endl;
#endif
	
	// Load config file for sort code
	SortConfig sortConfig("../config/sort.config");
	
	// Setup for multi-threaded progress bar
	const size_t updateRate = sortConfig.GetUpdateRate();
	long long globalProcessed{0};

	// Create output file
	string ofname = sortConfig.GetOutputDir() + sortConfig.GetOfileName();
	TFile ofile(ofname.c_str(), "RECREATE");
	cout << GREEN << "Output file: " << ofname << RESET << endl;
	
	// Initialize some variables up here so that they are accessible inside the lambda function
	size_t runnum;
	size_t numentries = 0;
	
	// Counters for certain particle combinations
	size_t count_ap;
	
	/******** RUN NUMBER LOOP ********/

	string runNumbersFile = sortConfig.GetRunNumbersFile();
	ifstream runFile(runNumbersFile);
	if (runFile.fail()) throw invalid_argument(string(BOLDRED) + string("Run numbers file ") + runNumbersFile + std::string(" does not exist or failed to open") + std::string(RESET));

	// First, loop through runs and find the total number of entries
	string itname = sortConfig.GetItreeName();
	ostringstream datastring;
	for (;;) {
		runFile >> runnum;
		if (runFile.eof() || runFile.bad()) break;

		datastring.str("");
		datastring << sortConfig.GetDataDir() << "run-" << runnum << ".root";

		// Check status of input run data file
		TFile *file = TFile::Open(datastring.str().c_str());
		if (!file || file->IsZombie()) continue;

		// Check if tree exists in the file
		TTree *tree = (TTree*)file->Get(itname.c_str());
		if (!tree) {
			file->Close();
			continue;
		}
		numentries += tree->GetEntries();
		file->Close();
	}

	// Then, loop through run numbers from numbers.beam and perform analysis on each
	runFile.clear();
	runFile.seekg(0);
	for (;;) {
		runFile >> runnum;
		if (runFile.eof() || runFile.bad()) break;

		datastring.str("");
		datastring << sortConfig.GetDataDir() << "run-" << runnum << ".root";

		// Check status of input run data file
		size_t numentries_singlefile;
		{
			TFile *file = TFile::Open(datastring.str().c_str());
			if (!file || file->IsZombie()) {
				cerr << "Error opening file for run " << runnum << "!" << endl;
				continue;
			}

			// Check if tree exists in the file
			TTree *tree = (TTree*)file->Get(itname.c_str());
			if (!tree) {
				cerr << "Tree '" << itname << "' not found in file for run " << runnum << "!" << endl;
				file->Close();
				continue;
			}
			numentries_singlefile = tree->GetEntries();
			file->Close();
		}

		cout << "Processing TTree in file: " << datastring.str() << " (" << numentries_singlefile << ")" << endl;

		// Create a TTreeProcessorMT: this class orchestrates the parallel processing of an input tree
		TTreeReader ttr(datastring.str().c_str(), itname.c_str());
		
		Input input(ttr);

		const char* otname = sortConfig.GetOtreeName().c_str();

		// Initialize analysis classes
		histo Histo(&ofile, sortConfig);
		Det det(input, Histo, sortConfig, runnum);
		
		// Thread-local event loop
		size_t localCounter = 0;
		while (reader.Next()) {

			// First, take input file from SpecTcl and refactor into usable hit list format
			input.ReadAndRefactor();
			
			// Perform analysis
			det.analyze();
			
			// Finalize per-event output, if any
			Histo.Fill();
			
			// Handle progress bar
			localCounter++;
			if (localCounter >= updateRate) {
				globalProcessed += localCounter;
				long double percentage = (long double)globalProcessed / numentries * 100.0;
				cout << "\r[ " << setw(7) << fixed << setprecision(4) 
				     << percentage << "% ] Processing entries..." << setw(10) << " " << flush;

				localCounter = 0;
			}
		}

		// Adding counters here that will tick up for different particle combinations
		// All counters should be of type atomic<> for thread safety
		count_ap += det.a_p;
		
		cout << endl;
	}

	// Output program duration
	auto end = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = end - start;
	int total_seconds = static_cast<int>(elapsed.count());
	int hours = total_seconds / 3600;
	int minutes = (total_seconds % 3600) / 60;
	int seconds = total_seconds % 60;
	cout << "Run time: "
	     << setfill('0') << setw(2) << hours << ":"
	     << setfill('0') << setw(2) << minutes << ":"
	     << setfill('0') << setw(2) << seconds << endl;

	cout << "************************************************************************" << endl;
	cout << "EVENT COUNTERS                                                          " << endl;
	cout << "1p + 1a: " << count_ap << endl;

	return 0;
}



