#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TMath.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TStyle.h>
#include <TView.h>

using namespace std;

// Custom 2D function for peak + diagonal background within a strip cut
Double_t peakPlusDiagonalBkg(Double_t *x, Double_t *par) {
	Double_t X = x[0];
	Double_t Y = x[1];
	
	// Fit Parameters:
	// par[0]: Peak Amplitude
	// par[1]: Peak Mean X (x0)
	// par[2]: Peak Mean Y (y0)
	// par[3]: Symmetric Gaussian Sigma (sigma)
	// par[4]: Background Offset Parameter ([3] in 1D)
	// par[5]: Background Slope Parameter ([4] in 1D)
	// par[6]: Background Non-linear Parameter ([5] in 1D)
	// par[7]: Diagonal Strip Half-Width (perpendicular distance limit from y=x)
	// par[8]: Diagonal Strip Minimum Limit (distance along y=x)
	// par[9]: Diagonal Strip Maximum Limit (distance along y=x)

	// Calculate perpendicular distance from y=x line: |y - x| / sqrt(2)
	Double_t perpDist = TMath::Abs(Y - X) / TMath::Sqrt(2.0);
	Double_t u = (X + Y) / TMath::Sqrt(2.0);
	
	// Restrict fit region to the diagonal strip
	if ((perpDist > par[7]) || (u < par[8]) || (u > par[9])) {
		TF1::RejectPoint(); // Reject bins outside the diagonal strip during fit
		return 0.;
	}

	// 1. Symmetric 2D Gaussian Peak
	Double_t r2 = ((X - par[1])*(X - par[1])) + ((Y - par[2])*(Y - par[2]));
	Double_t peak = par[0] * TMath::Exp(-0.5 * r2 / (par[3]*par[3]));

	// 2. Background evaluated along the diagonal axis u = (x + y) / sqrt(2)
	Double_t argExp = par[4] + par[5] * u;
	
	// Guard against numerical overflow in TMath::Exp
	Double_t bkg = 0.0;
	if (par[6] != 0.0) {
		bkg = TMath::Log(1.0 + TMath::Exp(argExp) * par[6]) / par[6];
	}

	return peak + bkg;
}

Double_t symGauss2D(Double_t *x, Double_t *par) {
	Double_t X = x[0];
	Double_t Y = x[1];
	
	// Fit Parameters:
	// par[0]: Peak Amplitude
	// par[1]: Peak Mean X (x0)
	// par[2]: Peak Mean Y (y0)
	// par[3]: Symmetric Gaussian Sigma (sigma)

	// Calculate perpendicular distance from y=x line: |y - x| / sqrt(2)
	Double_t perpDist = TMath::Abs(Y - X) / TMath::Sqrt(2.0);
	Double_t u = (X + Y) / TMath::Sqrt(2.0);

	Double_t r2 = ((X - par[1])*(X - par[1])) + ((Y - par[2])*(Y - par[2]));
	Double_t peak = par[0] * TMath::Exp(-0.5 * r2 / (par[3]*par[3]));

	return peak;
}

Double_t diagBackground(Double_t *x, Double_t *par) {
	Double_t X = x[0];
	Double_t Y = x[1];
	
	// Fit Parameters:
	// par[0]: Background Offset Parameter ([3] in 1D)
	// par[1]: Background Slope Parameter ([4] in 1D)
	// par[2]: Background Non-linear Parameter ([5] in 1D)

	// Calculate perpendicular distance from y=x line: |y - x| / sqrt(2)
	Double_t perpDist = TMath::Abs(Y - X) / TMath::Sqrt(2.0);
	Double_t u = (X + Y) / TMath::Sqrt(2.0);

	// Background evaluated along the diagonal axis u = (x + y) / sqrt(2)
	Double_t argExp = par[0] + par[1] * u;
	
	// Guard against numerical overflow in TMath::Exp
	Double_t bkg = 0.0;
	if (par[2] != 0.0) {
		bkg = TMath::Log(1.0 + TMath::Exp(argExp) * par[2]) / par[2];
	}

	return bkg;
}

void fit2DHistogram(TH2* h2, TF2*& fitFunc, Double_t stripHalfWidth, Double_t stripMin, Double_t stripMax) {
	if (!h2) return;

	// Create TF2 over the full histogram range
	Double_t xMin = h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetFirst());
	Double_t xMax = h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetLast());
	Double_t yMin = h2->GetYaxis()->GetBinLowEdge(h2->GetYaxis()->GetFirst());
	Double_t yMax = h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetLast());

	fitFunc = new TF2("fit2D", peakPlusDiagonalBkg, xMin, xMax, yMin, yMax, 10);

	// Set parameter names
	fitFunc->SetParNames("Amp", "MeanX", "MeanY", "Sigma", "Bkg_Off", "Bkg_Slope", "Bkg_NonLin", "StripHalfWidth", "StripMin", "StripMax");

	double Q6Be    = 1.372183; // known Q value for 6Be g.s., calculated using AME mass excesses
	double peakVal = h2->GetBinContent(h2->GetBin(Q6Be, Q6Be));

	// Initial guesses
	fitFunc->SetParameter(0, peakVal*0.8);              // Amplitude
	fitFunc->SetParLimits(0, peakVal*0.2, peakVal*1.2);
	fitFunc->FixParameter(1, Q6Be);	                    // Mean X
	fitFunc->FixParameter(2, Q6Be);	                    // Mean Y
	fitFunc->SetParameter(3, 0.1);	                    // Sigma
	fitFunc->SetParLimits(3, 0.05, 0.3);
	fitFunc->SetParameter(4, -3.3);                     // Bkg Offset
	fitFunc->SetParLimits(4, -5, -1);
	fitFunc->SetParameter(5, 4);                        // Bkg Slope
	fitFunc->SetParLimits(5, 0.1, 5);
	fitFunc->SetParameter(6, 0.021);                  	// Bkg Non-linear term
	fitFunc->FixParameter(7, stripHalfWidth);           // Diagonal cut half width
	fitFunc->FixParameter(8, stripMin);                 // Diagonal cut minimum limit
	fitFunc->FixParameter(9, stripMax);                 // Diagonal cut maximum limit

	// Perform 2D binned fit (use "R" to respect function range/rejected points, "M" for Minuit improve)
	h2->Fit(fitFunc, "RM");
}

void createClippedFitGraph(TH2* h2, TF2* func, TGraph2D*& g3d, Double_t stripHalfWidth, Double_t uMin, Double_t uMax, int nSteps = 100) {
	g3d = new TGraph2D();
	int pointCount = 0;

	Double_t xMin = h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetFirst());
	Double_t xMax = h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetLast());
	Double_t yMin = h2->GetYaxis()->GetBinLowEdge(h2->GetYaxis()->GetFirst());
	Double_t yMax = h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetLast());
	Double_t dx = (xMax - xMin) / nSteps;
	Double_t dy = (yMax - yMin) / nSteps;

	// Loop over the 2D plane within the zoomed window
	for (int i = 0; i <= nSteps; ++i) {
		Double_t x = xMin + i * dx;
		for (int j = 0; j <= nSteps; ++j) {
			Double_t y = yMin + j * dy;

			// 1. Calculate rotated coordinates
			Double_t u = (x + y) / sqrt(2.);           // Along the y=x line
			Double_t perpDist = abs(y - x) / sqrt(2.); // Transverse distance

			// 2. Apply both the strip thickness cut AND the u-coordinate bounds
			if (perpDist <= stripHalfWidth && u >= uMin && u <= uMax) {
				Double_t z = func->Eval(x, y);
				g3d->SetPoint(pointCount, x, y, z);
				pointCount++;
			}
		}
	}
}

// Function to create a clean, non-bleeding wireframe using 3D polylines
TList* createCleanWireframeLines(TH2* h2, TF2* func, TList*& lineList, Double_t stripHalfWidth, Double_t uMin, Double_t uMax, Double_t scale = 1, int nSteps = 80) {
	lineList = new TList();
	lineList->SetOwner(kTRUE); // Clean memory management

	Double_t xMin = h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetFirst());
	Double_t xMax = h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetLast());
	Double_t yMin = h2->GetYaxis()->GetBinLowEdge(h2->GetYaxis()->GetFirst());
	Double_t yMax = h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetLast());
	Double_t dx = (xMax - xMin) / nSteps;
	Double_t dy = (yMax - yMin) / nSteps;

	int nx = nSteps + 1;
	int ny = nSteps + 1;
	int totalSize = nx * ny;
	auto idx = [ny](int i, int j) { return i * ny + j; };
  
	vector<Double_t> Z(totalSize, 0.);
	vector<bool> valid(totalSize, false);
	vector<size_t> rowSizes(nx, 0.);
	vector<size_t> colSizes(ny, 0.);

	for (int i = 0; i <= nSteps; ++i) {
		Double_t x = xMin + i * dx;
		for (int j = 0; j <= nSteps; ++j) {
			Double_t y = yMin + j * dy;

			Double_t u = (x + y) / TMath::Sqrt(2.0);
			Double_t perpDist = TMath::Abs(y - x) / TMath::Sqrt(2.0);

			if (perpDist <= stripHalfWidth && u >= uMin && u <= uMax) {
				Z[idx(i, j)] = func->Eval(x, y);
				valid[idx(i, j)] = true;
				rowSizes[i]++;
				colSizes[j]++;
			}
		}
	}

	// Draw lines along rows (constant i, varying j)
	for (int i = 0; i <= nSteps; ++i) {
		if (rowSizes[i] <= 1) continue;

		Double_t x = xMin + i * dx;
		TPolyLine3D *line = new TPolyLine3D(rowSizes[i]);
		line->SetLineColor(kRed);
		line->SetLineWidth(1);

		// Loop through points in row and add to line
		int pIdx = 0;
		for (int j = 0; j <= nSteps; ++j) {
			if (!valid[idx(i, j)]) continue;

			Double_t y = yMin + j * dy;
			line->SetPoint(pIdx++, x, y, Z[idx(i, j)]*scale);
		}

		// The row size check means at this point all lines are valid
		lineList->Add(line);
	}

	// Draw lines along columns (constant j, varying i)
	for (int j = 0; j <= nSteps; ++j) {
		if (colSizes[j] <= 1) continue;

		Double_t y = yMin + j * dy;
		TPolyLine3D *line = new TPolyLine3D(colSizes[j]);
		line->SetLineColor(kRed);
		line->SetLineWidth(1);

		// Loop through points in column and add to line
		int pIdx = 0;
		for (int i = 0; i <= nSteps; ++i) {
			if (!valid[idx(i, j)]) continue;

			Double_t x = xMin + i * dx;
			line->SetPoint(pIdx++, x, y, Z[idx(i, j)]*scale);
		}

		// The column size check means at this point all lines are valid
		lineList->Add(line);
	}

	return lineList;
}

void Plot6Be6BeSimple() {
	string ifdir  = "/data4/N9/mnt/analysis/e25001/rootout/";
	string ifname = "O12_4p2a_processed_all_noneighbors_hasCsITDCQDC_SiFBGates_CsIPSD_CsIrecal_p+3MeV.root";

	double xtitleoff = 1.;
	double ytitleoff = .9;

	// Set default style attributes
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
	Sty->SetPadLeftMargin(0.12);
	Sty->SetPadRightMargin(0.14);
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
	Sty->SetNdivisions(10, "x");
	Sty->SetEndErrorSize(0);
	Sty->SetTextFont(42);
	gROOT->Reset();
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();

	double min1D = 0.;
	double max1D = 2.7;

	// Open TTree for desired invariant mass channel
	TFile* ifile = TFile::Open((ifdir + ifname).c_str());
	if (!ifile || ifile->IsZombie()) {
		cerr << "Input file " << ifdir + ifname << " does not exist or cannot be opened" << endl;
		return;
	}
	string hname = "12O/Be6_subevents_8BeVeto_x_proj";
	TH1* hist = (TH1*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Find or make canvas
	const size_t canw = 1000;
	const size_t canh = 610;
	TCanvas* mycan = (TCanvas*)gROOT->FindObjectAny("mycan");
	if (!mycan) mycan = new TCanvas("mycan", "", canw, canh);

	string xtitle = "E_{rel} ^{6}Be_{proj} (MeV)";
	string ytitle = "Counts per bin";

	hist->SetLineColor(kBlack); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	hist->SetMarkerStyle(20);
	hist->GetXaxis()->SetRangeUser(0.5, max1D);
	hist->GetXaxis()->SetTitle(xtitle.c_str());
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->SetTitle(ytitle.c_str());
	hist->GetYaxis()->CenterTitle();
	hist->Draw("EX0P");

	TF1* fit = new TF1("6Be_tot","gaus(0)+(TMath::Log(1+(TMath::Exp([3]+([4]*x))*[5]))/[5])");
	fit->SetParameter(0, 80);		// initial value estimated from plot
	fit->SetParLimits(0, 20, 120);
	fit->FixParameter(1, 1.372183);  // known Q value for 6Be g.s. using AME mass excesses
	fit->SetParameter(2, 0.1);	   // initial value estimated from plot
	fit->SetParLimits(2, 0.05, 0.3); // limits estimated from plot

	// The initial parameters for the background I got from putting the background function
	// into Desmos and playing around with it until I got a good shape. This is not an ideal
	// background function to use because to change the curvature where the linear portion
	// nears the X axis, you have to change some combination of parameters. I suspect this
	// makes it easy for the fit to fail if you don't give it a good set of starting parameters.
	fit->SetParameter(3, -3.3);	 // this is the only parameter with an obvious purpose (X-axis position)
	fit->SetParLimits(3, -5, -1);   // set range to some reasonable value just by looking at the X-axis
	fit->SetParameter(4, 4);		
	fit->SetParLimits(4, 1, 5);	 // for some reason I thought this was position and set limits (it isn't, but the limits worked anyways)
	fit->SetParameter(5, 0.021);

	fit->SetRange(min1D, max1D);
	fit->SetLineColor(kRed);
	hist->Fit(fit, "R");

	TF1* peak = new TF1("6Be_peak","gaus(0)");
	peak->SetParameter(0, fit->GetParameter(0));
	peak->SetParameter(1, fit->GetParameter(1));
	peak->SetParameter(2, fit->GetParameter(2));
	peak->SetRange(min1D, max1D);
	peak->SetLineColor(kGreen);
	peak->SetLineStyle(9);
	peak->Draw("same");

	TF1* background = new TF1("6Be_background","TMath::Log(1+(TMath::Exp([0]+([1]*x))*[2]))/[2]");
	background->SetParameter(0, fit->GetParameter(3));
	background->SetParameter(1, fit->GetParameter(4));
	background->SetParameter(2, fit->GetParameter(5));
	background->SetRange(min1D, max1D);
	background->SetLineColor(kBlue+2);
	background->SetLineStyle(9);
	background->Draw("same");

	TLegend* legend = new TLegend(0.12, 0.7, 0.32, 0.94);
	legend->AddEntry(hist, "Data", "ep");
	legend->AddEntry(fit, "Total fit", "l");
	legend->AddEntry(peak, "^{6}Be g.s.", "l");
	legend->AddEntry(background, "Background", "l");
	legend->Draw();

	const double Be6Rad = 0.292; // must be same as in O12_4p2a_process.C
	double Ipeak = peak->Integral(0.5, 1.6);
	double Iback = background->Integral(0.5, 1.6);
	double Itot  = Ipeak + Iback;
	cout << "% signal     = " << Ipeak / Itot * 100. << "%" << endl;
	cout << "% background = " << Iback / Itot * 100. << "%" << endl;

	mycan->Print("Canvas_1.eps", "eps");

	// Experiment with fitting 2D function directly to 2D histogram instead
	// of projecting first

	double min2D = 0.5;
	double max2D = 4.;

	// Get histogram
	hname = "12O/Be6_subevents_8BeVeto";
	TH2* hist1a = (TH2*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist1a) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Plot
	TCanvas* mycan1a = (TCanvas*)gROOT->FindObjectAny("mycan1a");
	if (!mycan1a) mycan1a = new TCanvas("mycan1a", "", canw, canh);

	xtitle = "E_{rel} ^{6}Be_{1} (MeV)";
	ytitle = "E_{rel} ^{6}Be_{2} (MeV)";
	string ztitle = "Counts per bin";
	hist1a->GetXaxis()->SetTitle("X Axis");
	hist1a->GetXaxis()->CenterTitle();
	hist1a->GetYaxis()->SetTitle("Y Axis");
	hist1a->GetYaxis()->CenterTitle();
	hist1a->GetZaxis()->SetTitle("Counts");
	hist1a->GetZaxis()->CenterTitle();
	hist1a->SetLineWidth(1);

	// Fit
	TF2* fitTot2D;
	fit2DHistogram(hist1a, fitTot2D, Be6Rad, min2D, max2D);

	fitTot2D->SetNpx(80);         // Increase sampling resolution along X
	fitTot2D->SetNpy(80);         // Increase sampling resolution along Y
	fitTot2D->SetLineColor(kGreen);
	fitTot2D->SetLineWidth(1);
	fitTot2D->SetFillStyle(0);    // 0 = Hollow/Transparent fill (no opaque backing)

	// "SURF1" draws a colored wireframe/surface mesh
	hist1a->GetListOfFunctions()->Clear();
	hist1a->Draw("LEGO2");
	//fitTot2D->Draw("SURF SAME");

	mycan1a->Modified();
	mycan1a->Update();

	// Get the Z scaling by first examining the Z range of
	// the full TH2 and then zooming in and getting the Z
	// range again. This Z scaling is used for the poly line
	// surface because somehow the poly lines obey the XY
	// transforms when you zoom in but not the Z transform.
	double zScale;
	TView *view = mycan1a->GetView();
	if (view) {
		Double_t min[3], max[3];
		view->GetRange(min, max);
		zScale = max[2];
	}

	hist1a->GetXaxis()->SetRangeUser(0, 3);
	hist1a->GetYaxis()->SetRangeUser(0, 3);

	mycan1a->Modified();
	mycan1a->Update();

	if (view) {
		Double_t min[3], max[3];
		view->GetRange(min, max);
		zScale /= max[2];
	}
	cout << "Z scaling value: " << zScale << endl;

	// Make custom wireframe of TF2 only inside fit region
	//TGraph2D* fitFrame;
	//createClippedFitGraph(hist1a, fitTot2D, fitFrame, Be6Rad, min2D, max2D);
	//fitFrame->Draw("P0 SAME");

	// Make a poly line surface representation of the fit and draw it
	TList* fitFrameLines;
	createCleanWireframeLines(hist1a, fitTot2D, fitFrameLines, Be6Rad, min2D, max2D, zScale);
	TIter next(fitFrameLines);
	while (TObject *obj = next()) {
		obj->Draw("SAME");
	}

	// To draw only one poly line and output its points
/*
	if (!fitFrameLines) {
		cout << "ERROR: fitFrameLines is null" << endl;
		ifile->Close();
		return;
	}
	size_t targetIndex = 20;
	int totalLines = fitFrameLines->GetEntries();
	cout << "Total # lines in list: " << totalLines << endl;
	if (targetIndex>=totalLines) {
		cout << "ERROR: targetIndex must be < " << totalLines << endl;
		ifile->Close();
		return;
	}
	TObject *obj = fitFrameLines->At(targetIndex);
	TPolyLine3D *line = dynamic_cast<TPolyLine3D*>(obj);
	if (!line) {
		cout << "ERROR: object at index " << targetIndex << " is null" << endl;
		ifile->Close();
		return;
	}
	int nPoints = line->GetN();
	Float_t *points = line->GetP();
	for (int p = 0; p < nPoints; ++p) {
		Double_t x = points[3 * p];
		Double_t y = points[3 * p + 1];
		Double_t z = points[3 * p + 2];
		cout << "  Point " << p << ": (X = " << x << ", Y = " << y << ", Z = " << z << ")" << endl;
	}
	line->Draw("SAME");
*/

	mycan1a->Modified();
	mycan1a->Update();

	mycan1a->Print("Canvas_1a.eps", "eps");

	// Now try making a fit to data comparison for this 2D fit by projecting along y = x
	hname = "12O/Be6_subevents_8BeVeto_xy_proj";
	TH1* hist1b = (TH1*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist1b) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	TCanvas* mycan1b = (TCanvas*)gROOT->FindObjectAny("mycan1b");
	if (!mycan1b) mycan1b = new TCanvas("mycan1b", "", canw, canh);

	// Clone event-based 1D projection histogram to inherit its exact binning
	TH1D* hFitProj = (TH1D*)hist1b->Clone("hFitProj");
	hFitProj->Reset(); // Clear counts, keep binning

	// Loop over every bin in your 1D projection histogram
	for (int i = 1; i <= hFitProj->GetNbinsX(); ++i) {
		Double_t u = hFitProj->GetBinCenter(i);
		if ((u < min2D) || (u > max2D)) continue;
		Double_t du = hFitProj->GetBinWidth(i);

		// We want to integrate the 2D fit function across the orthogonal strip coordinate 'v' 
		// at this specific u value, from -stripHalfWidth to +stripHalfWidth.
		// Coordinate transformation:
		// x = (u - v) / sqrt(2)
		// y = (u + v) / sqrt(2)

		Double_t stripIntegral = 0.0;
		int nIntegrationSteps = 30; // Fine steps across the strip width
		Double_t dv = (2.0 * Be6Rad) / nIntegrationSteps;

		for (int k = 0; k <= nIntegrationSteps; ++k) {
			Double_t v = -Be6Rad + k * dv;

			Double_t x = (u - v) / TMath::Sqrt(2.0);
			Double_t y = (u + v) / TMath::Sqrt(2.0);

			// Evaluate the fit function at this (x, y) point inside the strip
			Double_t val = fitTot2D->Eval(x, y);

			if (!TMath::IsNaN(val)) {
				stripIntegral += val;
			}
		}

		// Fill the model 1D histogram
		hFitProj->SetBinContent(i, stripIntegral);
	}

	// Scale the integrated/sampled projected fit histogram 
	// to match the integral with the data
	double integral_fit  = hFitProj->Integral(hFitProj->GetXaxis()->FindBin(min2D), hFitProj->GetXaxis()->FindBin(max2D));
	double integral_data = hist1b->Integral(hist1b->GetXaxis()->FindBin(min2D), hist1b->GetXaxis()->FindBin(max2D));
	hFitProj->Scale(integral_data / integral_fit);

	// Draw both on the same canvas to compare
	hist1b->SetMarkerStyle(kFullCircle);
	hist1b->SetMarkerColor(kBlack);
	hist1b->Draw("E1"); // Your event-based data projection

	hFitProj->SetLineColor(kRed);
	hFitProj->SetLineWidth(2);
	hFitProj->Draw("HIST SAME"); // Your fit model projection

	mycan1b->Modified();
	mycan1b->Update();

	mycan1b->Print("Canvas_1b.eps", "eps");

	// Peak/background % calculation from 2D integral

	// First extract peak and background components from total fit
	Double_t xMin = hist1a->GetXaxis()->GetBinLowEdge(hist1a->GetXaxis()->GetFirst());
	Double_t xMax = hist1a->GetXaxis()->GetBinUpEdge(hist1a->GetXaxis()->GetLast());
	Double_t yMin = hist1a->GetYaxis()->GetBinLowEdge(hist1a->GetYaxis()->GetFirst());
	Double_t yMax = hist1a->GetYaxis()->GetBinUpEdge(hist1a->GetYaxis()->GetLast());
	TF2* fitPeak2D = new TF2("fitPeak2D", symGauss2D, xMin, xMax, yMin, yMax, 4);
	TF2* fitBkgd2D = new TF2("fitBkgd2D", diagBackground, xMin, xMax, yMin, yMax, 3);
	double x0 = fitTot2D->GetParameter(1);
	double y0 = fitTot2D->GetParameter(2);
	fitPeak2D->SetParameter(0, fitTot2D->GetParameter(0));
	fitPeak2D->SetParameter(1, x0);
	fitPeak2D->SetParameter(2, y0);
	fitPeak2D->SetParameter(3, fitTot2D->GetParameter(3));
	fitBkgd2D->SetParameter(0, fitTot2D->GetParameter(4));
	fitBkgd2D->SetParameter(1, fitTot2D->GetParameter(5));
	fitBkgd2D->SetParameter(2, fitTot2D->GetParameter(6));

	Double_t peakIntegral = 0.0;
	Double_t bgIntegral = 0.0;

	int nStepsGrid = 100;
	Double_t dx = (2.0 * Be6Rad) / nStepsGrid;
	Double_t dy = (2.0 * Be6Rad) / nStepsGrid;
	Double_t cellArea = dx * dy;

	for (int i = 0; i <= nStepsGrid; ++i) {
		Double_t x = (x0 - Be6Rad) + i * dx;
		for (int j = 0; j <= nStepsGrid; ++j) {
			Double_t y = (y0 - Be6Rad) + j * dy;

			// Check if the point falls inside the circular region: (x - x0)^2 + (y - y0)^2 <= R^2
			Double_t r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
			if (r2 < (Be6Rad * Be6Rad)) {
				// Evaluate your peak component and background component separately
				Double_t peakVal = fitPeak2D->Eval(x, y); 
				Double_t bgVal   = fitBkgd2D->Eval(x, y);

				if (!TMath::IsNaN(peakVal)) peakIntegral += peakVal * cellArea;
				if (!TMath::IsNaN(bgVal))   bgIntegral += bgVal * cellArea;
			}
		}
	}

	// Calculate totals and percentages
	Double_t totalIntegral = peakIntegral + bgIntegral;
	Double_t percentSignal = (totalIntegral > 0) ? (peakIntegral / totalIntegral) * 100.0 : 0.0;
	Double_t percentBg     = (totalIntegral > 0) ? (bgIntegral / totalIntegral) * 100.0 : 0.0;

	// 4. Print results
	cout << "--- Circular Region Integration (Center: " << x0 << ", " << y0 << " | Radius: " << Be6Rad << ") ---" << endl;
	cout << "Peak Integral:       " << peakIntegral << endl;
	cout << "Background Integral: " << bgIntegral << endl;
	cout << "Total Integral:      " << totalIntegral << endl;
	cout << "--------------------------------------------------" << endl;
	cout << "Signal Component:    " << percentSignal << " %" << endl;
	cout << "Background Component:" << percentBg << " %" << endl;

	/******** 12O BACKGROUND ********/
	// Next, what I can do is create a sideband gate in the 6Be spectrum
	// to see what the corresponding 12O background looks like. Basically,
	// I make a 12O spectrum for the sideband gate. I then scale this
	// sideband spectrum by the ratio between the background integral in
	// the signal region and the background integral in the sideband
	// region. This gives me a background for the 12O spectrum.

	// Get 6Be gated histogram
	hname = "12O/O12_Erel_6Be6Begate_8BeVeto";
	TH1* hist2 = (TH1*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist2) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Get sideband gated histogram
	hname = "12O/O12_Erel_6Be6Besideband_8BeVeto";
	TH1F* hist3 = (TH1F*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist3) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Find or make second canvas
	TCanvas* mycan2 = (TCanvas*)gROOT->FindObjectAny("mycan2");
	if (!mycan2) mycan2 = new TCanvas("mycan2", "", canw, canh);

	// Plot data
	xtitle = "E_{rel} 4p+2#alpha (MeV, ^{8}Be veto)";
	hist2->Rebin(4);
	hist2->SetLineColor(kBlack);
	hist2->SetMarkerStyle(20);
	hist2->GetXaxis()->SetRangeUser(4, 15);
	hist2->GetXaxis()->SetTitle(xtitle.c_str());
	hist2->GetXaxis()->CenterTitle();
	hist2->GetYaxis()->SetTitle(ytitle.c_str());
	hist2->GetYaxis()->CenterTitle();
	hist2->Draw("EX0P");

	Double_t bgSideIntegral = 0.0;

	nStepsGrid = 200;
	dx = (4.0 * Be6Rad) / nStepsGrid;
	dy = (4.0 * Be6Rad) / nStepsGrid;
	cellArea = dx * dy;

	for (int i = 0; i <= nStepsGrid; ++i) {
		Double_t x = (x0 - (2.*Be6Rad)) + i * dx;
		for (int j = 0; j <= nStepsGrid; ++j) {
			Double_t y = (y0 - (2.*Be6Rad)) + j * dy;

			// Check if the point falls inside the disk region
			Double_t r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
			if ((r2 >= (Be6Rad * Be6Rad)) && (r2 < ((2.*Be6Rad)*(2.*Be6Rad)))) {
				Double_t bgVal   = fitBkgd2D->Eval(x, y);
				if (!TMath::IsNaN(bgVal))   bgSideIntegral += bgVal * cellArea;
			}
		}
	}

	// Plot scaled sideband background
	hist3->Rebin(4);
	hist3->Scale(bgIntegral / bgSideIntegral);
	hist3->SetLineColor(kBlue+2); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	hist3->SetLineStyle(9);
	hist3->Draw("HIST SAME");

	// Draw legend
	TLegend* legend2 = new TLegend(0.54, 0.76, 0.86, 0.94);
	legend2->AddEntry(hist2, "Data, ^{6}Be+^{6}Be gate", "ep");
	legend2->AddEntry(hist3, "Sideband bkgd. (scaled)", "l");
	legend2->Draw();

	mycan2->Print("Canvas_2.eps", "eps");

	// Find or make third canvas
	TCanvas* mycan3 = (TCanvas*)gROOT->FindObjectAny("mycan3");
	if (!mycan3) mycan3 = new TCanvas("mycan3", "", canw, canh);

	xtitle = "E_{rel} 4p+2#alpha (MeV, ^{6}Be gate, bkgd. sub.)";
	hist2->Add(hist3, -1);
	hist2->GetXaxis()->SetTitle(xtitle.c_str());
	hist2->GetYaxis()->SetTitle("");
	hist2->Draw("EX0P");

	mycan3->Print("Canvas_3.eps", "eps");
}



