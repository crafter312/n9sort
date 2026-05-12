/* Written by Henry Webb (h.s.webb@wustl.edu), August 2025
 * This class handles reading information from a SpecTcl-
 * generated ROOT file and reformatting it into a more
 * easily usable form. The resulting variables are stored
 * in this class on a per-event basis and can be accessed
 * via getter functions. Make sure that the static functions
 * match the format of the branch and leaf names in the
 * SpecTcl output tree. It is also important to make sure
 * that you loop through your TTreeReaderValue objects in
 * the same order in which you create them and their strings.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu) and Johnathan
 * Phillips (j.s.phillips@wustl.edu) March 2026 for experiment
 * at TAMU Cyclotron Institute
 */

#include "Input.h"

#include <cmath>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/******** STATIC FUNCTIONS ********/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

vector<string> Input::GenerateColumnNamesHINP(const string& parname) {
	vector<string> columns;
	string b, c;
	for (size_t board = 1; board <= HINP_BOARD_COUNT; board++) {
		for (size_t chan = 0; chan < HINP_CHAN_COUNT; chan++) {
			b = (board < 10 ? "0" : "") + to_string(board);
			c = (chan < 10 ? "0" : "") + to_string(chan);
			columns.push_back("SpecTcl_hinp1_mb1_" + parname + "_" + b + "." + c);
		}
	}
	return columns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

vector<string> Input::GenerateColumnNamesPSD(const string& parname) {
	vector<string> columns;
	string chi, cha;
	for (size_t chip = 1; chip <= PSD_CHIP_COUNT; chip++) {
		for (size_t chan = 0; chan < PSD_CHAN_COUNT; chan++) {
			chi = (chip < 10 ? "0" : "") + to_string(chip);
			cha = (chan < 10 ? "0" : "") + to_string(chan);
			columns.push_back("SpecTcl_psd1_" + parname + "_" + chi + "." + cha);
		}
	}
	return columns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

vector<string> Input::GenerateColumnNamesADCQDC(const string& branchname, const size_t& chancount) {
	vector<string> columns;
	string cha;
	for (size_t chan = 0; chan < chancount; chan++) {
		cha = (chan < 10 ? "0" : "") + to_string(chan);
		columns.push_back(branchname + "." + cha);
	}
	return columns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

vector<string> Input::GenerateColumnNamesQDCV965(const string& parname) {
	vector<string> columns;
	string cha;
	for (size_t chan = 0; chan < QDCV965_CHAN_COUNT; chan++) {
		cha = (chan < 10 ? "0" : "") + to_string(chan);
		columns.push_back("SpecTcl_qdc1_" + cha + "." + parname);
	}
	return columns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

vector<string> Input::GenerateColumnNamesTDC() {
	vector<string> columns;
	string cha;
	for (size_t chan = 0; chan < TDC_CHAN_COUNT; chan++) {
		for (int hit = 0; hit < TDC_HIT_COUNT; hit++) {
			cha = (chan < 10 ? "0" : "") + to_string(chan);
			columns.push_back("SpecTcl_tdc1_" + cha + "." + to_string(hit));
		}
	}
	return columns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/******** NON-STATIC FUNCTIONS ********/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Input::Input(TTreeReader& r, SortConfig& config) : reader(r) {

	// Pre-alocate required vector memory
	gobbi.eRVs.reserve(HINP_NCOLUMNS);
	gobbi.eLoRVs.reserve(HINP_NCOLUMNS);
	gobbi.tRVs.reserve(HINP_NCOLUMNS);
	adc.aqRVs.reserve(ADC_CHAN_COUNT);
	qdc.aqRVs.reserve(QDC_CHAN_COUNT);
	tdc.tRVs.reserve(TDC_NCOLUMNS);

	// Generate column names for reading from input tree
	vector<string> e_columns     = GenerateColumnNamesHINP("e");
	vector<string> eLo_columns   = GenerateColumnNamesHINP("eLo");
	vector<string> hinpt_columns = GenerateColumnNamesHINP("t");
	vector<string> adc_columns   = GenerateColumnNamesADCQDC(config.GetAdcBranchName(), ADC_CHAN_COUNT);
	vector<string> qdc_columns   = GenerateColumnNamesADCQDC(config.GetQdcBranchName(), QDC_CHAN_COUNT);
	vector<string> tdct_columns  = GenerateColumnNamesTDC();
	
	//// Create reader values for all columns, iteratively

	// HINP
	for (size_t i = 0; i < HINP_NCOLUMNS; i++) {
		gobbi.eRVs.emplace_back(reader, e_columns[i].c_str());
		gobbi.eLoRVs.emplace_back(reader, eLo_columns[i].c_str());
		gobbi.tRVs.emplace_back(reader, hinpt_columns[i].c_str());
	}
	
	// ADC
	for (size_t i = 0; i < ADC_CHAN_COUNT; i++)
		adc.aqRVs.emplace_back(reader, adc_columns[i].c_str());

	// QDC
	for (size_t i = 0; i < QDC_CHAN_COUNT; i++)
		qdc.aqRVs.emplace_back(reader, qdc_columns[i].c_str());

	// TDC
	for (size_t i = 0; i < TDC_NCOLUMNS; i++)
		tdc.tRVs.emplace_back(reader, tdct_columns[i].c_str());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Input::~Input() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Input::ReadAndRefactor() {
	gobbi.clear();
	adc.clear();
	qdc.clear();
	tdc.clear();
	
	// Loop through TDC channels to retrieve time hit information
	size_t chan;
	double tdc_t;
	for (size_t i = 0; i < TDC_NCOLUMNS; i++) {
		tdc_t = *(tdc.tRVs[i]);
		
#ifdef ENABLE_DEBUG
		cout << i / (size_t)TDC_HIT_COUNT << " " << tdc_t << endl;
#endif

		if (i == 0 && tdc_t != 0) {
		
#ifdef ENABLE_DEBUG
			cout << "bad event " << badevt << ", skip! tdc_t is " << tdc_t << endl;
#endif
			
			badevt++;
			return;
		}
		if (isnan(tdc_t) || (abs(tdc_t) >= 10000)) continue; // the nan portion only works if `tdc_t` is the same type as its column in the input tree
		chan = i / (size_t)TDC_HIT_COUNT;
		tdc.Nhits[chan]++;
		tdc.t[chan].push_back(tdc_t);
	}
	
	// Loop through HINP boards and channels and retrieve hit information
	double e;
	for (size_t i = 0; i < HINP_NCOLUMNS; i++) {
		e = *(gobbi.eRVs[i]);
		if (isnan(e) || (e == 0) || (e >= 16384)) continue; // the nan portion only works if `e` is the same type as its column in the input tree
		gobbi.Nhits++;
		gobbi.board.push_back((i / (size_t)HINP_CHAN_COUNT) + 1);
		gobbi.chan.push_back(i % (size_t)HINP_CHAN_COUNT);
		gobbi.e.push_back((size_t)e);
		gobbi.eLo.push_back((size_t)(*(gobbi.eLoRVs[i])));
		gobbi.t.push_back((size_t)(*(gobbi.tRVs[i])));
	}
	
	// Loop through ADC channels to retrieve hit information
	double a;
	for (size_t i = 0; i < ADC_CHAN_COUNT; i++) {
		a = *(adc.aqRVs[i]);
		if (isnan(a) || (a == 0) || (a >= 16384)) continue; // the nan portion only works if `a` is the same type as its column in the input tree
		adc.Nhits++;
		adc.chan.push_back(i);
		adc.aq.push_back((size_t)a);
	}
	
	// Loop through QDC channels to retrieve hit information
	double q;
	for (size_t i = 0; i < QDC_CHAN_COUNT; i++) {
		q = *(qdc.aqRVs[i]);
		if (isnan(q) || (q == 0) || (q >= 16384)) continue; // the nan portion only works if `q` is the same type as its column in the input tree
		qdc.Nhits++;
		qdc.chan.push_back(i);
		qdc.aq.push_back((size_t)q);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



