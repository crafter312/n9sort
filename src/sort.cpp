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
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
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
using Clock = chrono::steady_clock;

int main() {

	// Capture the start time
	auto start = Clock::now();

	// Variables for accumulating time for profiling
	// Make the ones inside the multi-threaded lambda function atomic
	// Actually I messed this up by changing from ints to doubles,
	// and apparently `atomic<double>` does not have a += operator
	// like `atomic<int>` does, which seems stupid to me
	double timeReadFile{0};
	double timeProcessing{0};
	mutex timeCounterMutex;
	double timeInitialize{0};
	double timeReadRefactor{0};
	double timeAnalysis{0};

	// This one is weird because I'm trying to determine time across
	// the end point of the lambda function, so I store the start
	// times in a map by thread ID for thread safety and then loop
	// over the map and calculate the total elapsed time.
	mutex endTimesMutex;
	map<thread::id, Clock::time_point> lambdaEndStartTimes;
	double timeDestruct{0};

#ifndef rel
	cout << BOLDCYAN << "WARNING: compile time variable `rel` not defined, using Newtonian mechanics" << RESET << endl;
#endif
	
	// Load config file for sort code
	SortConfig sortConfig("../config/sort.config");
	if (sortConfig.GetSkipHighEStrips())
		cout << BOLDCYAN << "WARNING: skipHighEStrips is set to true, all events at least 1 silicon strip E > 90 MeV WILL BE SKIPPED" << RESET << endl;
	
	// Setup for multi-threaded progress bar
	const size_t updateRate = sortConfig.GetUpdateRate();
	std::atomic<long long> globalProcessed{0};
	mutex consoleMutex; // To prevent text scrambling

	// Setup for debug/testing output file
	mutex ofileMutex;   // To protect output file for debug information
	ofstream testOut("../RootFiles/" + sortConfig.GetDebugFileName());

	// Create the TBufferMerger: this class orchestrates the parallel writing to an output ROOT file
	string ofname = sortConfig.GetOutputDir() + sortConfig.GetOfileName();
	ROOT::TBufferMerger merger(ofname.c_str(), "RECREATE");
	cout << GREEN << "Output file: " << ofname << RESET << endl;

	// Enable implicit multi-threading
	int nthreads = 8;
	ROOT::EnableImplicitMT(nthreads);
	
	// Initialize some variables up here so that they are accessible inside the lambda function
	size_t runnum;
	size_t numentries = 0;
	
	// Counters for certain particle combinations, using atomic to be thread-safe
	atomic<size_t> count_ap     = 0;
	atomic<size_t> count_app    = 0;
	atomic<size_t> count_appp   = 0;
	atomic<size_t> count_apppp  = 0;
	atomic<size_t> count_appppp = 0;
	atomic<size_t> count_aapppp = 0;
	atomic<size_t> count_aa3He  = 0;
	
	// Other counters
	atomic<size_t> pidSkipped   = 0;
	atomic<size_t> nEvtSelected = 0;
	
	/******** EVENT PROCESSING LAMBDA FUNCTION ********/
	
	// Define the function that will process a subrange of the tree.
	// The function must receive only one parameter, a TTreeReader,
	// and it must be thread safe. To enforce the latter requirement,
	// TBufferMerger::GetFile will be used for the output file.
	auto f = [&](TTreeReader &reader) {
		auto startInit = Clock::now();
		Input input(reader, sortConfig);

		// Output using thread safe file
		auto f = merger.GetFile();

		const char* otname = sortConfig.GetOtreeName().c_str();

		// Initialize analysis classes
		histo Histo(f, sortConfig);
		Det det(input, Histo, sortConfig, runnum, nEvtSelected, ofileMutex, testOut);
		auto endInit = Clock::now();
		chrono::duration<double> elapsedInit = endInit - startInit;
		
		// Thread-local event loop
		size_t localCounter = 0;
		double localTimeReadRefactor = 0.;
		double localTimeAnalysis = 0.;
		auto startReadRefactor = Clock::now();
		while (reader.Next()) {
			auto endReadRefactor = Clock::now();
			chrono::duration<double> elapsedReadRefactor = endReadRefactor - startReadRefactor;
			localTimeReadRefactor += (double)(elapsedReadRefactor.count());

			// Reset global tree output variables
			Histo.reset();

			// First, take input file from SpecTcl and refactor into usable hit list format
			startReadRefactor = Clock::now();
			input.ReadAndRefactor();
			endReadRefactor = Clock::now();
			elapsedReadRefactor = endReadRefactor - startReadRefactor;
			localTimeReadRefactor += (double)(elapsedReadRefactor.count());
			
			// Perform analysis
			auto startAnalysis = Clock::now();
			det.analyze();
			auto endAnalysis = Clock::now();
			chrono::duration<double> elapsedAnalysis = endAnalysis - startAnalysis;
			localTimeAnalysis += (double)(elapsedAnalysis.count());
			
			// Finalize per-event output, if any
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

			startReadRefactor = Clock::now();
		}

		// Add local double time counters to global integer time counters
		{
			lock_guard<mutex> lock(timeCounterMutex);
			timeInitialize   += (double)(elapsedInit.count());
			timeReadRefactor += localTimeReadRefactor;
			timeAnalysis     += localTimeAnalysis;
		}

		// Adding counters here that will tick up for different particle combinations
		// All counters should be of type atomic<> for thread safety
		count_ap     += det.a_p;
		count_app    += det.a_pp;
		count_appp   += det.a_ppp;
		count_apppp  += det.a_pppp;
		count_appppp += det.a_ppppp;
		count_aapppp += det.aa_pppp;
		count_aa3He  += det.aa_3He;
		
		// Other counters
		pidSkipped   += det.pidSkipped;

		{
			lock_guard<mutex> lock(endTimesMutex);
			lambdaEndStartTimes[this_thread::get_id()] = Clock::now();
		}
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

		auto startRead = Clock::now();
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
		auto endRead = Clock::now();
		chrono::duration<double> elapsedRead = endRead - startRead;
		timeReadFile += (double)(elapsedRead.count());

		cout << "Processing TTree in file: " << datastring.str() << " (" << numentries_singlefile << ")" << endl;

		// Create a TTreeProcessorMT: this class orchestrates the parallel processing of an input tree
		ROOT::TTreeProcessorMT tp(datastring.str().c_str(), itname.c_str());

		// Execute multi-threaded tree processing
		lambdaEndStartTimes.clear();
		auto startProcess = Clock::now();
		tp.Process(f);
		auto endProcess = Clock::now();

		for (const auto& [tid, endStartTime] : lambdaEndStartTimes) {
			chrono::duration<double> elapsedLambdaEnd = endProcess - endStartTime;
			timeDestruct += (double)(elapsedLambdaEnd.count());
		}

		chrono::duration<double> elapsedProcess = endProcess - startProcess;
		timeProcessing += (double)(elapsedProcess.count());
		cout << endl;
	}

	// Output program duration
	auto end = Clock::now();
	chrono::duration<double> elapsed = end - start;
	double total_seconds = (double)(elapsed.count());
	int hours = (int)total_seconds / 3600;
	int minutes = ((int)total_seconds % 3600) / 60;
	int seconds = (int)total_seconds % 60;
	cout << "Run time: "
	     << setfill('0') << setw(2) << hours << ":"
	     << setfill('0') << setw(2) << minutes << ":"
	     << setfill('0') << setw(2) << seconds << endl;
	cout << endl;
	cout << "Time profiling percentage breakdown:" << endl;
	cout << "Time reading file and tree: " << timeReadFile / total_seconds * 100. << "%" << endl;
	cout << "Time from TTreeProcessorMT::Process: " << timeProcessing / total_seconds * 100. << "%" << endl;
	cout << endl;
	cout << "TTreeProcessorMT::Process further time breakdown:" << endl;
	cout << "(percent of processing time, not percent of total code run time)" << endl;
	double timeProcThreadMult = timeProcessing * nthreads;
	cout << "Time from initialization: " << timeInitialize / timeProcThreadMult * 100. << "%" << endl;
	cout << "Time from SpecTcl input tree reading and refactoring: " << timeReadRefactor / timeProcThreadMult * 100. << "%" << endl;
	cout << "Time from analysis: " << timeAnalysis / timeProcThreadMult * 100. << "%" << endl;
	cout << "Time from things going out of scope at the end of the lambda function: " << timeDestruct / timeProcThreadMult * 100. << "%" << endl;

	// Output analysis counters
	cout << "************************************************************************" << endl;
	cout << "PARTICLE COUNTERS                                                       " << endl;
	cout << "1a + 1p: "  << count_ap     << endl;
	cout << "1a + 2p: "  << count_app    << endl;
	cout << "1a + 3p: "  << count_appp   << endl;
	cout << "1a + 4p: "  << count_apppp  << endl;
	cout << "1a + 5p: "  << count_appppp << endl;
	cout << "2a + 4p: "  << count_aapppp << endl;
	cout << "2a + 3He: " << count_aa3He  << endl;
	cout << "OTHER COUNTERS                                                          " << endl;
	cout << "# solutions invalidated due to Eloss errors: " << pidSkipped << endl;

	testOut.close();
	return 0;
}



