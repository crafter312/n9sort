#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TList.h>
#include <TPad.h>
#include <TStyle.h>

using namespace std;

void Plot8CSimple() {
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

	// Open TTree for desired invariant mass channel
	TFile* ifile = TFile::Open((ifdir + ifname).c_str());
	if (!ifile || ifile->IsZombie()) {
		cerr << "Input file " << ifdir + ifname << " does not exist or cannot be opened" << endl;
		return;
	}
	string hname = "12O/C8_subevents_8BeVeto";
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

	string xtitle = "E_{rel} 4p+#alpha (MeV, 4p+2#alpha subevents)";
	string ytitle = "Counts per bin";

	hist->SetLineColor(kBlack); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	hist->SetMarkerStyle(20);
	hist->GetXaxis()->SetRangeUser(2, 6);
	hist->GetXaxis()->SetTitle(xtitle.c_str());
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->SetTitle(ytitle.c_str());
	hist->GetYaxis()->CenterTitle();
	hist->Draw("EX0P");

	TF1* fit = new TF1("8C_tot","gaus(0)+(TMath::Log(1+(TMath::Exp([3]+([4]*x))*[5]))/[5])");
	fit->SetParameter(0, 147.417);  // initial value from previous fit attempt
	fit->SetParLimits(0, 50, 200);
	fit->FixParameter(1, 3.529);    // known Q value for 8C g.s. is 3.48321 MeV, but this is from a fit to a simulated peak (from Bob)
	fit->SetParameter(2, 0.2373);   // this is also from the same fit to a simulated 8C g.s. peak,
	fit->SetParLimits(2, 0.1, 0.5); // but I will let this wander a little bit (differing detector systems)

	// The initial parameters for the background I got from putting the background function
	// into Desmos and playing around with it until I got a good shape. This is not an ideal
	// background function to use because to change the curvature where the linear portion
	// nears the X axis, you have to change some combination of parameters. I suspect this
	// makes it easy for the fit to fail if you don't give it a good set of starting parameters.
	fit->SetParameter(3, -12.7);    
	fit->SetParameter(4, 4);        // this is the only parameter with an obvious purpose (X-axis position)
	fit->SetParLimits(4, 2, 5);     // can set range just by looking at the histogram for a sensible range
	fit->SetParameter(5, 0.021);

	fit->SetRange(2, 6);
	fit->SetLineColor(kRed);
	hist->Fit(fit, "R");

	TF1* peak = new TF1("8C_peak","gaus(0)");
	peak->SetParameter(0, fit->GetParameter(0));
	peak->SetParameter(1, fit->GetParameter(1));
	peak->SetParameter(2, fit->GetParameter(2));
	peak->SetRange(2, 6);
	peak->SetLineColor(kGreen);
	peak->SetLineStyle(9);
	peak->Draw("same");

	TF1* background = new TF1("8C_background","TMath::Log(1+(TMath::Exp([0]+([1]*x))*[2]))/[2]");
	background->SetParameter(0, fit->GetParameter(3));
	background->SetParameter(1, fit->GetParameter(4));
	background->SetParameter(2, fit->GetParameter(5));
	background->SetRange(2, 6);
	background->SetLineColor(kBlue+2);
	background->SetLineStyle(9);
	background->Draw("same");

	TLegend* legend = new TLegend(0.12, 0.7, 0.32, 0.94);
	legend->AddEntry(hist, "Data", "ep");
	legend->AddEntry(fit, "Total fit", "l");
	legend->AddEntry(peak, "^{8}C g.s.", "l");
	legend->AddEntry(background, "Background", "l");
	legend->Draw();

	double Ipeak = peak->Integral(2, 4);
	double Iback = background->Integral(2, 4);
	double Itot  = Ipeak + Iback;
	cout << "% signal     = " << Ipeak / Itot * 100. << "%" << endl;
	cout << "% background = " << Iback / Itot * 100. << "%" << endl;

	mycan->Print("Canvas_1.eps", "eps");

	/******** 12O BACKGROUND ********/
	// Next, what I can do is create a sideband gate in the 8C spectrum
	// to see what the corresponding 12O background looks like. Basically,
	// I make a 12O spectrum for the sideband gate. I then scale this
	// sideband spectrum by the ratio between the background integral in
	// the signal region and the background integral in the sideband
	// region. This gives me a background for the 12O spectrum.

	// Get 8C gated histogram
	hname = "12O/O12_Erel_8Cgate_8BeVeto";
	TH1* hist2 = (TH1*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist2) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Get sideband gated histogram
	hname = "12O/O12_Erel_8Csideband_8BeVeto";
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
	xtitle = "E_{rel} 4p+2#alpha (MeV)";
	hist2->Rebin(4);
	hist2->SetLineColor(kBlack);
	hist2->SetMarkerStyle(20);
	hist2->GetXaxis()->SetRangeUser(4, 15);
	hist2->GetXaxis()->SetTitle(xtitle.c_str());
	hist2->GetXaxis()->CenterTitle();
	hist2->GetYaxis()->SetTitle(ytitle.c_str());
	hist2->GetYaxis()->CenterTitle();
	hist2->Draw("EX0P");

	// Plot scaled sideband background
	double Isideband = background->Integral(4, 5);
	hist3->Rebin(4);
	hist3->Scale(Iback / Isideband);
	hist3->SetLineColor(kBlue+2); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	hist3->SetLineStyle(9);
	hist3->Draw("HIST SAME");

	// Draw legend
	TLegend* legend2 = new TLegend(0.54, 0.76, 0.86, 0.94);
	legend2->AddEntry(hist2, "Data", "ep");
	legend2->AddEntry(hist3, "Sideband bkgd. (scaled)", "l");
	legend2->Draw();

	mycan2->Print("Canvas_2.eps", "eps");

	// Find or make third canvas
	TCanvas* mycan3 = (TCanvas*)gROOT->FindObjectAny("mycan3");
	if (!mycan3) mycan3 = new TCanvas("mycan3", "", canw, canh);

	xtitle = "E_{rel} 4p+2#alpha (MeV, ^{8}C gate, bkgd. sub.)";
	hist2->Add(hist3, -1);
	hist2->GetXaxis()->SetTitle(xtitle.c_str());
	hist2->GetYaxis()->SetTitle("");
	hist2->Draw("EX0P");

	mycan3->Print("Canvas_3.eps", "eps");
}



