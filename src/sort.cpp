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

#include <ROOT/TBufferMerger.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TFile.h>
#include <TTree.h>

#include "constants.h"
#include "Det.h"
#include "histo.h"
#include "Input.h"
#include "SortConfig.h"

using namespace std;

int main() {

	// Capture the start time
	auto start = chrono::high_resolution_clock::now();
	
	// Load config file for sort code
	SortConfig sortConfig("../config/sort.config");
	
	// Setup for multi-threaded progress bar
	const size_t updateRate = sortConfig.GetUpdateRate();
	std::atomic<long long> globalProcessed{0};
	mutex consoleMutex; // To prevent text scrambling

	// Create the TBufferMerger: this class orchestrates the parallel writing to an output ROOT file
	string ofname = sortConfig.GetOutputDir() + sortConfig.GetOfileName();
	ROOT::TBufferMerger merger(ofname.c_str(), "RECREATE");
	cout << GREEN << "Output file: " << ofname << RESET << endl;

	// Enable implicit multi-threading
	int nthreads = 4;
	ROOT::EnableImplicitMT(nthreads);
	
	// Initialize some variables up here so that they are accessible inside the lambda function
	int runnum;
	size_t numentries = 0;
	
	// Counters for certain particle combinations, using atomic to be thread-safe
	// Start with 6Li -> npa
	atomic<size_t> count_ap;
	
	/******** EVENT PROCESSING LAMBDA FUNCTION ********/
	
	// Define the function that will process a subrange of the tree.
	// The function must receive only one parameter, a TTreeReader,
	// and it must be thread safe. To enforce the latter requirement,
	// TBufferMerger::GetFile will be used for the output file.
	auto f = [&](TTreeReader &reader) {
		Input input(reader);

		// Output using thread safe file
		auto f = merger.GetFile();

		const char* otname = sortConfig.GetOtreeName().c_str();

		// Initialize analysis classes
		histo Histo(f);
		Gobbi gobbi(input, Histo, sortConfig, runnum);
		
		// Thread-local event loop
		size_t localCounter = 0;
		while (reader.Next()) {

			// First, take input file from SpecTcl and refactor into usable hit list format
			input.ReadAndRefactor();
			
			// Gobbi analysis
			gobbi.analyze();
			
			// Output
			Histo.Fill();
			
			// Handle progress bar
			localCounter++;
			if (localCounter >= updateRate) {
				long long total = globalProcessed.fetch_add(localCounter);
				lock_guard<mutex> lock(consoleMutex);
				long double percentage = (long double)total / numentries * 100.0;
				cout << "\r[ " << setw(7) << fixed << setprecision(4) 
				     << percentage << "% ] Processing entries..." << setw(10) << " " << flush;

				localCounter = 0;
			}
		}

		// Adding counters here that will tick up for different particle combinations
		// All counters should be of type atomic<> for thread safety
		count_ap += gobbi.a_p;
	};
	
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
		ROOT::TTreeProcessorMT tp(datastring.str().c_str(), itname.c_str());

		// Execute multi-threaded tree processing
		tp.Process(f);
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



