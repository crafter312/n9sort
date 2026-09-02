#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <algorithm>
#include <vector>

// Run with something like `root -l -q sum_tcal_data.C`

void sum_tcal_data() {

	// Set default style attributes
	double xtitleoff = 1.;
	double ytitleoff = 1.2;
	TStyle* Sty = new TStyle("MyStyle", "MyStyle");
	Sty->SetOptTitle(0);
	Sty->SetOptStat(0);
	Sty->SetLineWidth(3);
	Sty->SetPalette(kBird);
	Sty->SetCanvasColor(10);
	Sty->SetCanvasBorderMode(0);
	Sty->SetFrameLineWidth(0);
	Sty->SetFrameFillColor(10);
	Sty->SetPadColor(10);
	Sty->SetPadTickX(1);
	Sty->SetPadTickY(1);
	Sty->SetPadBottomMargin(0.15);
	Sty->SetPadTopMargin(0.06);
	Sty->SetPadLeftMargin(0.14);
	Sty->SetPadRightMargin(0.06);
	Sty->SetHistLineWidth(3);
	Sty->SetFuncWidth(3);
	Sty->SetFuncColor(kGreen);
	Sty->SetLineWidth(3);
	Sty->SetLabelSize(0.05, "xyz");
	Sty->SetLabelOffset(0.01, "y");
	Sty->SetLabelOffset(0.01, "x");
	Sty->SetLabelColor(kBlack, "xyz");
	Sty->SetTitleSize(0.06, "y");
	Sty->SetTitleSize(0.06, "x");
	Sty->SetTitleOffset(ytitleoff, "y");
	Sty->SetTitleOffset(xtitleoff, "x");
	Sty->SetTitleFillColor(10);
	Sty->SetTitleTextColor(kBlack);
	Sty->SetTickLength(.05, "xz");
	Sty->SetTickLength(.025, "y");
	Sty->SetNdivisions(10, "y");
	Sty->SetNdivisions(510, "x");
	Sty->SetPadGridX(true);
	Sty->SetPadGridY(false);
	Sty->SetGridStyle(3); // 1=solid, 2=dashed, 3=dotted
	Sty->SetGridWidth(1);
	Sty->SetGridColor(kGray);
	Sty->SetEndErrorSize(0);
	Sty->SetTextFont(42);
	gROOT->Reset();
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();

	const double sigma_to_fwhm = 2. * sqrt(2. * log(2.));

	// 1. Open the ROOT file
	const char* filename = "/data4/N9/mnt/analysis/e25001/rootout/e21006_HINP_calibration.root";
	TFile* file = TFile::Open(filename, "READ");
	if (!file || file->IsZombie()) {
		std::cerr << "Error: Cannot open file " << filename << std::endl;
		return;
	}

	// 2. Define flattening bounds (8 chips * 32 channels = 256 total channels)
	// Index mapping: index = 32 * (i - 1) + j, ranging from 0 to 255
	const int nboards = 8;
	const int nchannels = 32;
	const int total_channels = nboards * nchannels;
	
	// Create the summary histogram
	TH1D* h_first_y = new TH1D(
		"h_first_y", 
		"First Y-Value vs Flattened Channel Index;Channel Index (32*(i-1)+j);First Y Value", 
		total_channels, 0, total_channels
	);

	// 3. Loop over i [1..8] and j [0..31]
	int found_count = 0;
	for (int i = 1; i <= nboards; ++i) {
		for (int j = 0; j < nchannels; ++j) {
			int flat_index = nchannels * (i - 1) + j;
			
			// Format object names expected in the file
			TString graph_name = Form("chip_%d_%d", i, j);
			TString fit_name   = Form("fit_%d_%d", i, j);

			// Fetch TGraphErrors
			TGraphErrors* gr = nullptr;
			file->GetObject(graph_name.Data(), gr);

			// Fetch TF1 (retrieved as requested, available for evaluation if needed)
			TF1* fit = nullptr;
			file->GetObject(fit_name.Data(), fit);

			if (gr && fit && gr->GetN() > 0) {
				// Find the point corresponding to the smallest X value
				int n_points = gr->GetN();
				double* x_vals = gr->GetX();
				double* y_vals = gr->GetY();
				double* ey_vals = gr->GetEY();

				double abs_slope = std::abs(fit->GetParameter(1));
				double sum_err_ns = 0.0;
				std::vector<double> errs_ns;
				errs_ns.reserve(n_points);

				if (abs_slope > 0.) {
					double err_ns = ey_vals[0] / abs_slope;
					errs_ns.push_back(err_ns);
					sum_err_ns += err_ns;
				}
				else {
					std::cout << "Warning: Fit slope is 0 for i=" << i << ", j=" << j << std::endl;
				}

				int min_x_idx = 0;
				double min_x = x_vals[0];

				for (int p = 1; p < n_points; ++p) {
					if (x_vals[p] < min_x) {
						min_x = x_vals[p];
						min_x_idx = p;
					}
					if (abs_slope <= 0.) continue;
					double err_ns = (ey_vals[p] / abs_slope) * sigma_to_fwhm;
					errs_ns.push_back(err_ns);
					sum_err_ns += err_ns;
				}

				double first_y = y_vals[min_x_idx];

				double mean_err_ns = sum_err_ns / n_points;
				double sum_sq_diff = 0.0;
				for (double err_ns : errs_ns) {
					double diff = err_ns - mean_err_ns;
					sum_sq_diff += diff * diff;
				}

				double std_dev_ns = 0.0;
				if (n_points > 1) {
					std_dev_ns = std::sqrt(sum_sq_diff / (n_points - 1)); // Sample std dev
				}

				std::cout << "Board (i=" << i << "), Channel (j=" << j << ")\t| "
				          << "Avg FWHM: " << mean_err_ns << " ns\t| "
				          << "Std Dev: " << std_dev_ns << " ns" << std::endl;

				// Fill histogram bin (ROOT bins are 1-indexed, flat_index + 1 maps index 0 to bin 1)
				h_first_y->SetBinContent(flat_index + 1, first_y);
				
				// Set bin error from TGraphErrors if error bars exist
				if (gr->GetEY()) {
					h_first_y->SetBinError(flat_index + 1, gr->GetEY()[min_x_idx]);
				}
				
				found_count++;
			} else {
				std::cout << "Warning: Could not find graph or empty graph for " 
						  << graph_name << std::endl;
			}
		}
	}

	std::cout << "Successfully processed " << found_count << " / " << total_channels << " channels." << std::endl;

	// 4. Style and Plotting
	TCanvas* canvas = new TCanvas("canvas", "Chip Summary", 900, 600);

	h_first_y->SetLineColor(kBlue + 2);
	h_first_y->SetLineWidth(2);
	h_first_y->SetMarkerStyle(20);
	h_first_y->SetMarkerSize(0.6);
	h_first_y->SetMarkerColor(kBlue + 2);
	h_first_y->GetXaxis()->CenterTitle();
	h_first_y->GetYaxis()->CenterTitle();

	// 1. Force exact axis range to match the 8-chip x 32-channel boundary (0 to 256)
	h_first_y->GetXaxis()->SetLimits(0, 256);

	// 2. Set Ndivisions = -(N2 * 100 + N1)
	// - 8 primary divisions (256 / 32 = 8 major blocks)
	// - 4 secondary subdivisions (32 / 8 = 4 minor spaces per major block)
	// The negative sign (-) forces ROOT to use EXACTLY 8 major intervals instead of optimizing them automatically.
	h_first_y->GetXaxis()->SetNdivisions(-408);

	h_first_y->Draw("E1"); // Draw histogram with error bars

	// 5. Export canvas to .eps vector format
	canvas->SaveAs("chip_first_y_summary.eps");

	// Clean up file handle (canvas owns the histogram duplicate upon drawing)
	file->Close();
}
