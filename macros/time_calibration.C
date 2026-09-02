
#include <TFile.h>
#include <TH1.h>
#include <TSpectrum.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/ioctl.h>
#include <unistd.h>

using namespace std;

struct Progress {
	struct winsize w; // holds console size information
	float mi{0.};     // progress bar minimum in parameter space
	float ma{100.};   // progress bar maximum in parameter space
	float pos{0.};    // progress bar position in parameter space
	size_t width{0};  // progress bar width in console space (excluding fixed/constant characters)
	size_t posi{0};   // progress bar position in console space
	
	void SetRange(float a, float b) {
		mi = min(a, b);
		ma = max(a, b);
	}

	void SetPos(float p) {
		pos = max(mi, min(ma, p));
	}

	void Increment(float i) {
		pos += i;
		pos = max(mi, min(ma, pos));
	}

	void Print() {
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w); // get console size in case of change
		if (w.ws_col < 8) return; // return early if console not wide enough to print maximum minimum width progress bar (i.e. "[] 100 %")

		// Check for visual progress
		size_t width_tmp = w.ws_col - 8;
		float progress = (pos - mi) / (ma - mi);
		size_t posi_tmp = width_tmp * progress;
		if ((posi_tmp == posi) && (width_tmp == width)) return;
		posi = posi_tmp;
		width = width_tmp;

		// Output progress bar if changed
		cout << "[";
		for (size_t i = 0; i < width; i++) {
			if (i < posi) cout << "=";
			else if (i == posi) cout << ">";
			else cout << " ";
		}
		cout << "] " << int(progress * 100.) << " %";
		cout.flush();
		cout << "\33[2K\r"; // clear the current line and carriage return
	}
};

void time_calibration() {

	// Read file and TTree
	string path = "/data4/N9/mnt/analysis/e25001/rootout/";
	string ifname = "sort_run238";
	size_t numentries;
	TFile *file = TFile::Open((path + ifname + ".root").c_str());
	if (!file || file->IsZombie()) {
		std::cerr << "Error opening file!" << std::endl;
		return;
	}

	// Macro parameters
	size_t max_channels = pow(2, 14); // 16384, max time channels from HINP
	size_t fit_bins = 600; // number of bins in histogram to fit
	double sigma = 10.; // starting sigma for fit
	double cutoff = 15000; // manually set cutoff for erroneous peaks at high time values (handfull of channels have one peak like this)

	// Default minimizer
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Migrad");

	// Set ROOT output verbosity
	gErrorIgnoreLevel = kWarning;

	// Fit setup
	TF1 fit("doubleGaus", "gaus(0) + gaus(3)");
	fit.SetParName(0, "Constant 1");
	fit.SetParName(1, "Mean 1");
	fit.SetParName(2, "Sigma 1");
	fit.SetParName(3, "Constant 2");
	fit.SetParName(4, "Mean 2");
	fit.SetParName(5, "Sigma 2");

	// Canvas, if need to save histograms
	TCanvas* c = new TCanvas("c", "Canvas", 800, 600);

	// Files to save results
	string frontName = "FrontTimecalPulser.txt";
	ofstream oFrontFile(frontName);
  if (!oFrontFile.is_open()) {
    cerr << "*WARNING* Unable to create E front time calibration file " << frontName << endl;
    abort();
  }

	string backName = "BackTimecalPulser.txt";
	ofstream oBackFile(backName);
  if (!oBackFile.is_open()) {
    cerr << "*WARNING* Unable to create E back time calibration file " << backName << endl;
    abort();
  }

	string deltaName = "DeltaTimecalPulser.txt";
	ofstream oDeltaFile(deltaName);
  if (!oDeltaFile.is_open()) {
    cerr << "*WARNING* Unable to create delta E time calibration file " << deltaName << endl;
    abort();
  }

	// File for error messages
	string ename = "timeCal_ConsoleOutput.txt";
	ofstream eofile(ename);
  if (!eofile.is_open()) {
    cerr << "*WARNING* Unable to create console output file " << ename << endl;
    abort();
  }

	// Progress bar setup
	size_t boards = 8; // # boards to loop through (12 for full HINP crate)
	size_t channels = 32; // # channels per board, should be 32
	Progress pbar;
	pbar.SetRange(0, boards * channels);

	// Loop through boards/channels to fit time peaks
	TSpectrum peakfinder; // for finding peaks, max two peaks can be found
	Int_t npeaks;
	size_t i1, i2;
	double a1, a2;
	bool hasErrors = false;
	for (size_t b = 1; b < boards + 1; b++) {
		for (size_t ch = 0; ch < channels; ch++) {
			size_t tele = (b - 1) / 2;
			bool isFront = (b % 2) == 1;
			string sumName = "FrontTime_R/Front";
			if (!isFront) sumName = "BackTime_R/Back";

			// Extract histogram from file
			TString graph_name = Form("Summary/1d%sTime_R%zu_%zu", sumName.c_str(), tele, ch);
			TH1I* hist = nullptr;
			file->GetObject(graph_name.Data(), hist);
			if (!hist) {
				cout << "WARNING: Skipping board " << b << " channel " << ch << endl;
				cout << "Histogram " << graph_name.Data() << " not found" << endl;
				continue;
			}

			// Find peaks
			peakfinder.Search(hist);
			npeaks = peakfinder.GetNPeaks(); 
			Double_t* peakx = peakfinder.GetPositionX();
			Double_t* peaky = peakfinder.GetPositionY();

			// Make sure the peaks are good
			if (npeaks != 2) {
				eofile << "Irregular time spectrum for board " << b << " channel " << ch << ", saving histogram..." << endl;
				eofile << to_string(npeaks) << " time peaks found:" << endl;
				for (int i = 0; i < npeaks; i++) eofile << "\t" << peakx[i] << " " << peaky[i] << endl;
				hist->Draw();
				hist->SetStats(0);
				
				if (npeaks < 2) {
					eofile << "Less than two time spectrum peaks for board " << b << " channel " << ch << ", skipping..." << endl;
					c->SaveAs(("histogram_" + to_string(b) + "-" + to_string(ch) + ".png").c_str());
					hasErrors = true;
					pbar.Increment(1);
					continue; // need at least two peaks for fit
				}
			}

			// Search for the two correct peaks
			a1 = 0.;
			a2 = 0.;
			for (int i = 0; i < npeaks; i++) {
				if (peakx[i] > cutoff) continue;
				if (peaky[i] > a1) {
					a2 = a1;
					i2 = i1;
					a1 = peaky[i];
					i1 = i;
				}
				else if (peaky[i] > a2) {
					a2 = peaky[i];
					i2 = i;
				}
			}
			if ((a1 == 0.) || (a2 == 0.)) {
				eofile << "No two valid time spectrum peaks found for board " << b << " channel " << ch << ", skipping..." << endl;
				c->SaveAs(("histogram_" + to_string(b) + "-" + to_string(ch) + ".png").c_str());
				hasErrors = true;
				pbar.Increment(1);
				continue;
			}

			// Fit peaks
			double avg = (peakx[i1] + peakx[i2]) / 2.;
			double min = avg - (fit_bins / 2);
			double max = avg + (fit_bins / 2);
			fit.SetParameter(0, peaky[i1]);
			fit.SetParameter(1, peakx[i1]);
			fit.SetParameter(2, sigma);
			fit.SetParameter(3, peaky[i2]);
			fit.SetParameter(4, peakx[i2]);
			fit.SetParameter(5, sigma);
			fit.SetParLimits(0, peaky[i1]*0.5, peaky[i1]*1.5);
			fit.SetParLimits(1, min, max);
			fit.SetParLimits(2, 0, 100);
			fit.SetParLimits(3, peaky[i2]*0.5, peaky[i2]*1.5);
			fit.SetParLimits(4, min, max);
			fit.SetParLimits(5, 0, 100);
			fit.SetRange(min, max);
			Int_t fitStatus = hist->Fit(&fit, "NRQ", "", min, max);

			// Check if fit failed
			if (fitStatus > 0) {
				eofile << "Time spectrum fit failed for board " << b << " channel " << ch << ", skipping and saving histogram..." << endl;
				hist->Draw();
				hist->GetXaxis()->SetRangeUser(min, max);
				fit.Draw("same");
				c->SaveAs(("histogram_" + to_string(b) + "-" + to_string(ch) + ".png").c_str());
				pbar.Increment(1);
				hasErrors = true;
				continue;
			}
			else if (fitStatus < 0) {
				eofile << "Time spectrum fit non-minimizer error for board " << b << " channel " << ch << ", skipping..." << endl;
				pbar.Increment(1);
				hasErrors = true;
				continue;
			}

			// Output to files
			double slope = abs(fit.GetParameter(1) - fit.GetParameter(4)) / 20.; // divide by 20 ns, the length of cable used for the time calibration run
			size_t quad;
			if ((b == 1) || (b == 3) || (b == 5) || (b == 7)) {
				quad = (b - 1) / 2;
				oFrontFile << quad << " " << ch << " " << slope << " " << 0 << endl;
			}
			else if ((b == 2) || (b == 4) || (b == 6) || (b == 8)) {
				quad = (b / 2) - 1;
				oBackFile << quad << " " << ch << " " << slope << " " << 0 << endl;
			}
			else if ((b == 9) || (b == 10) || (b == 11) || (b == 12)) {
				quad = b - 9;
				oDeltaFile << quad << " " << ch << " " << slope << " " << 0 << endl;
			}

			// Increment progress bar
			pbar.Increment(1);
			pbar.Print();
		}
	}

	if (hasErrors)
		cout << "Time calibration errors present, channels skipped. See " + ename + " for detailed output." << endl;

	eofile.close();
	oDeltaFile.close();
	oBackFile.close();
	oFrontFile.close();
	file->Close();
}
