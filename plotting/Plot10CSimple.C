#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TList.h>
#include <TPad.h>
#include <TStyle.h>

using namespace std;

void Plot10CSimple() {
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
	string hname = "12O/C10_subevents_8Be9Bgate";
	TH1* hist = (TH1*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}
	TH1* histCopy = (TH1*)hist->Clone();

	// Find or make canvas
	const size_t canw = 1000;
	const size_t canh = 610;
	TCanvas* mycan = (TCanvas*)gROOT->FindObjectAny("mycan");
	if (!mycan) mycan = new TCanvas("mycan", "", canw, canh);

	string xtitle = "E_{rel} 2p+2#alpha (MeV, ^{8}Be && ^{9}B g.s. gates)";
	string ytitle = "Counts per bin";

	hist->SetLineColor(kBlack); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	hist->SetMarkerStyle(20);
	hist->GetXaxis()->SetTitle(xtitle.c_str());
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->SetTitle(ytitle.c_str());
	hist->GetYaxis()->CenterTitle();
	hist->Draw("EX0P");

	// Fit Parameters:
	// par[0]: Background Amplitude
	// par[1]: Background X Position
	// par[2]: Background Tail (exponential) Parameter
	// par[3]: Background Asymmetry (gaussian sigma)
	// par[4]: Peak 1 Amplitude
	// par[5]: Peak 1 Mean
	// par[6]: Peak 1 Sigma
	// par[7]: Peak 2 Amplitude
	// par[8]: Peak 2 Mean
	// par[9]: Peak 2 Sigma

	TF1* fit = new TF1("10C_tot","((0.5*[0]/[2])*TMath::Exp((([1]-x)/[2])+(0.5*[3]*[3]/([2]*[2])))*(1-TMath::Erf((1/sqrt(2))*((([1]-x)/[3])+([3]/[2])))))+gaus(4)+gaus(7)");

	// The initial parameters for the background I got from putting the background function
	// into Desmos and playing around with it until I got a good shape, scaled down slightly
	// in amplitude because I'm now fitting a slightly different histogram
	fit->SetParameter(0, 742);
	fit->SetParameter(1, 2.8);
	fit->SetParLimits(1, 0, 10);
	fit->SetParameter(2, 2.7);
	fit->SetParameter(3, 0.967);

	// First peak
	fit->SetParameter(4, 160);       // initial visually estimated value
	fit->SetParLimits(4, 80, 200);
	fit->FixParameter(5, 1.4509120); // calculated from 5.18 MeV 10C state minus the Q value for 2p+2a (see 10.1103/PhysRevC.75.051304)
	fit->SetParameter(6, 0.676);     // also from playing around in Desmos
	fit->SetParLimits(6, 0.1, 1.5);

	// Second peak (rough estimates)
	fit->SetParameter(7, 60);
	fit->SetParLimits(7, 10, 100);
	fit->FixParameter(8, 5.960912); // calculated from 9.69 MeV 10C state minus the Q value for 2p+2a (see 10.1103/PhysRevC.84.014320)
	fit->SetParameter(9, 0.8);
	fit->SetParLimits(9, 0.2, 2);

	fit->SetLineColor(kRed);
	hist->Fit(fit);

	double xMin = hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetFirst());
	double xMax = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetLast());

	TF1* background = new TF1("10C_background","(0.5*[0]/[2])*TMath::Exp((([1]-x)/[2])+(0.5*[3]*[3]/([2]*[2])))*(1-TMath::Erf((1/sqrt(2))*((([1]-x)/[3])+([3]/[2]))))", xMin, xMax);
	background->SetParameter(0, fit->GetParameter(0));
	background->SetParameter(1, fit->GetParameter(1));
	background->SetParameter(2, fit->GetParameter(2));
	background->SetParameter(3, fit->GetParameter(3));
	background->SetLineColor(kBlue+2);
	background->SetLineStyle(9);
	background->Draw("same");

	TF1* peak1 = new TF1("10C_peak1","gaus(0)", xMin, xMax);
	peak1->SetParameter(0, fit->GetParameter(4));
	peak1->SetParameter(1, fit->GetParameter(5));
	peak1->SetParameter(2, fit->GetParameter(6));
	peak1->SetLineColor(kGreen);
	peak1->SetLineStyle(9);
	peak1->Draw("same");

	TF1* peak2 = new TF1("10C_peak2","gaus(0)", xMin, xMax);
	peak2->SetParameter(0, fit->GetParameter(7));
	peak2->SetParameter(1, fit->GetParameter(8));
	peak2->SetParameter(2, fit->GetParameter(9));
	peak2->SetLineColor(kGreen);
	peak2->SetLineStyle(9);
	peak2->Draw("same");

	

	TLegend* legend = new TLegend(0.54, 0.76, 0.86, 0.94);
	legend->AddEntry(hist, "Data", "ep");
	legend->AddEntry(fit, "Total fit", "l");
	legend->AddEntry(peak1, "^{10}C E*=5.18 MeV", "l");
	legend->AddEntry(peak2, "^{10}C E*=9.69 MeV", "l");
	legend->AddEntry(background, "Background", "l");
	legend->Draw();

	double peakMin = 0.;
	double peakMax = 2.;
	double Ipeak = peak1->Integral(peakMin, peakMax);
	double Iback = background->Integral(peakMin, peakMax);
	const double IbackSAVE = Iback;
	double Itot  = Ipeak + Iback;
	cout << "% signal     = " << Ipeak / Itot * 100. << "%" << endl;
	cout << "% background = " << Iback / Itot * 100. << "%" << endl;

	mycan->Print("Canvas_1.eps", "eps");

	//// Fit with Fermi background instead

	// Find or make canvas
	TCanvas* mycan2 = (TCanvas*)gROOT->FindObjectAny("mycan2");
	if (!mycan2) mycan2 = new TCanvas("mycan2", "", canw, canh);

	histCopy->SetLineColor(kBlack); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	histCopy->SetMarkerStyle(20);
	histCopy->GetXaxis()->SetTitle(xtitle.c_str());
	histCopy->GetXaxis()->CenterTitle();
	histCopy->GetYaxis()->SetTitle(ytitle.c_str());
	histCopy->GetYaxis()->CenterTitle();
	histCopy->Draw("EX0P");

	// Fit Parameters:
	// par[0]: Background Amplitude
	// par[1]: Background X Position
	// par[2]: Background Exponential Parameter
	// par[4]: Peak Amplitude
	// par[5]: Peak Mean
	// par[6]: Peak Sigma

	TF1* fit2 = new TF1("10C_tot2","([0]/(TMath::Exp([2]*([1]-x))+1.))+gaus(3)", 0, 4.5);
	fit2->SetRange(0, 4.5);

	// Background
	fit2->SetParameter(0, 140);       // visually estimated
	fit2->SetParLimits(0, 100, 180);
	fit2->SetParameter(1, 2.1);       // visually estimated
	fit2->SetParLimits(1, 1, 3);
	fit2->SetParameter(2, 0.1);       // idk what this one looks like, just picked a value and crossed my fingers

	// First peak
	fit2->SetParameter(3, 160);       // initial visually estimated value
	fit2->SetParLimits(3, 80, 200);
	fit2->FixParameter(4, 1.4509120); // calculated from 5.18 MeV 10C state minus the Q value for 2p+2a (see 10.1103/PhysRevC.75.051304)
	fit2->SetParameter(5, 0.676);     // also from playing around in Desmos
	fit2->SetParLimits(5, 0.1, 1.5);

	fit2->SetLineColor(kRed);
	histCopy->Fit(fit2, "R");

	xMin = histCopy->GetXaxis()->GetBinLowEdge(histCopy->GetXaxis()->GetFirst());

	TF1* background2 = new TF1("10C_background2","[0]/(TMath::Exp([2]*([1]-x))+1.)", xMin, 4.5);
	background2->SetParameter(0, fit2->GetParameter(0));
	background2->SetParameter(1, fit2->GetParameter(1));
	background2->SetParameter(2, fit2->GetParameter(2));
	background2->SetLineColor(kBlue+2);
	background2->SetLineStyle(9);
	background2->Draw("same");

	TF1* peak = new TF1("10C_peak","gaus(0)", xMin, 4.5);
	peak->SetParameter(0, fit2->GetParameter(3));
	peak->SetParameter(1, fit2->GetParameter(4));
	peak->SetParameter(2, fit2->GetParameter(5));
	peak->SetLineColor(kGreen);
	peak->SetLineStyle(9);
	peak->Draw("same");

	TLegend* legend2 = new TLegend(0.54, 0.76, 0.86, 0.94);
	legend2->AddEntry(histCopy, "Data", "ep");
	legend2->AddEntry(fit2, "Total fit", "l");
	legend2->AddEntry(peak, "^{10}C E*=5.18 MeV", "l");
	legend2->AddEntry(background2, "Background", "l");
	legend2->Draw();

	Ipeak = peak->Integral(peakMin, peakMax);
	Iback = background2->Integral(peakMin, peakMax);
	Itot  = Ipeak + Iback;
	cout << "% signal     = " << Ipeak / Itot * 100. << "%" << endl;
	cout << "% background = " << Iback / Itot * 100. << "%" << endl;

	mycan2->Print("Canvas_2.eps", "eps");




	//// Now fit 2p channel (8Be gate + 9B veto)

	hname = "12O/C10_subevents_8Begate_9Bveto";
	TH1* hist3 = (TH1*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist3) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Find or make canvas
	TCanvas* mycan3 = (TCanvas*)gROOT->FindObjectAny("mycan3");
	if (!mycan3) mycan3 = new TCanvas("mycan3", "", canw, canh);

	xtitle = "E_{rel} 2p+2#alpha (MeV, ^{8}Be && !^{9}B g.s. gates)";

	hist3->SetLineColor(kBlack); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	hist3->SetMarkerStyle(20);
	hist3->GetXaxis()->SetTitle(xtitle.c_str());
	hist3->GetXaxis()->CenterTitle();
	hist3->GetYaxis()->SetTitle(ytitle.c_str());
	hist3->GetYaxis()->CenterTitle();
	hist3->Draw("EX0P");

	// Fit Parameters:
	// par[0]: Background Amplitude
	// par[1]: Background X Position
	// par[2]: Background Tail (exponential) Parameter
	// par[3]: Background Asymmetry (gaussian sigma)
	// par[4]: Peak Amplitude
	// par[5]: Peak Mean
	// par[6]: Peak Sigma

	TF1* fit3 = new TF1("10C_tot3","((0.5*[0]/[2])*TMath::Exp((([1]-x)/[2])+(0.5*[3]*[3]/([2]*[2])))*(1-TMath::Erf((1/sqrt(2))*((([1]-x)/[3])+([3]/[2])))))+gaus(4)");
	fit3->SetParameter(0, 350);
	fit3->SetParameter(1, 2.8);
	fit3->SetParLimits(1, 0, 10);
	fit3->SetParameter(2, 2.7);
	fit3->SetParameter(3, 0.967);
	fit3->SetParameter(4, 20);        // initial visually estimated value
	fit3->SetParLimits(4, 5, 30);
	fit3->FixParameter(5, 1.5709120); // calculated from 5.30 MeV 10C state minus the Q value for 2p+2a (see 10.1103/PhysRevC.75.051304)
	fit3->SetParameter(6, 0.676);     // also from playing around in Desmos
	fit3->SetParLimits(6, 0.1, 1.5);
	fit3->SetLineColor(kRed);
	hist3->Fit(fit3);

	xMin = hist3->GetXaxis()->GetBinLowEdge(hist3->GetXaxis()->GetFirst());
	xMax = hist3->GetXaxis()->GetBinUpEdge(hist3->GetXaxis()->GetLast());

	TF1* background3 = new TF1("10C_background3","(0.5*[0]/[2])*TMath::Exp((([1]-x)/[2])+(0.5*[3]*[3]/([2]*[2])))*(1-TMath::Erf((1/sqrt(2))*((([1]-x)/[3])+([3]/[2]))))", xMin, xMax);
	background3->SetParameter(0, fit3->GetParameter(0));
	background3->SetParameter(1, fit3->GetParameter(1));
	background3->SetParameter(2, fit3->GetParameter(2));
	background3->SetParameter(3, fit3->GetParameter(3));
	background3->SetLineColor(kBlue+2);
	background3->SetLineStyle(9);
	background3->Draw("same");

	TF1* peak3 = new TF1("10C_peak3","gaus(0)", xMin, xMax);
	peak3->SetParameter(0, fit3->GetParameter(4));
	peak3->SetParameter(1, fit3->GetParameter(5));
	peak3->SetParameter(2, fit3->GetParameter(6));
	peak3->SetLineColor(kGreen);
	peak3->SetLineStyle(9);
	peak3->Draw("same");

	TLegend* legend3 = new TLegend(0.54, 0.76, 0.86, 0.94);
	legend3->AddEntry(hist3, "Data", "ep");
	legend3->AddEntry(fit3, "Total fit", "l");
	legend3->AddEntry(peak3, "^{10}C E*=5.30 MeV", "l");
	legend3->AddEntry(background3, "Background", "l");
	legend3->Draw();

	Ipeak = peak3->Integral(peakMin, peakMax);
	Iback = background3->Integral(peakMin, peakMax);
	Itot  = Ipeak + Iback;
	cout << "% signal     = " << Ipeak / Itot * 100. << "%" << endl;
	cout << "% background = " << Iback / Itot * 100. << "%" << endl;

	mycan3->Print("Canvas_3.eps", "eps");

	/******** 12O BACKGROUND ********/
	// Next, what I can do is create a sideband gate in the 10C spectrum
	// to see what the corresponding 12O background looks like. Basically,
	// I make a 12O spectrum for the sideband gate. I then scale this
	// sideband spectrum by the ratio between the background integral in
	// the signal region and the background integral in the sideband
	// region. This gives me a background for the 12O spectrum.

	// Get 10C gated histogram
	hname = "12O/O12_Erel_8Be9B10Cgates";
	TH1* hist4a = (TH1*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist4a) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Get sideband gated histogram
	hname = "12O/O12_Erel_8Be9Bgated10Csideband";
	TH1F* hist4b = (TH1F*)ifile->Get((string("InvMass/") + hname).c_str());
	if (!hist4b) {
		cerr << "Histogram " << hname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Find or make fourth canvas
	TCanvas* mycan4 = (TCanvas*)gROOT->FindObjectAny("mycan4");
	if (!mycan4) mycan4 = new TCanvas("mycan4", "", canw, canh);

	// Plot data
	xtitle = "E_{rel} 4p+2#alpha (MeV, ^{8}Be and ^{9}B gates)";
	hist4a->Rebin(4);
	hist4a->SetLineColor(kBlack);
	hist4a->SetMarkerStyle(20);
	hist4a->GetXaxis()->SetRangeUser(4, 15);
	hist4a->GetXaxis()->SetTitle(xtitle.c_str());
	hist4a->GetXaxis()->CenterTitle();
	hist4a->GetYaxis()->SetTitle(ytitle.c_str());
	hist4a->GetYaxis()->CenterTitle();
	hist4a->Draw("EX0P");

	// Plot scaled sideband background
	double Isideband = background->Integral(2, 3);
	hist4b->Rebin(4);
	hist4b->Scale(IbackSAVE / Isideband);
	hist4b->SetLineColor(kBlue+2); // `kBlue+2` is the default blue shade used when plotting histograms via `TTree::Draw`
	hist4b->SetLineStyle(9);
	hist4b->Draw("HIST SAME");

	// Draw legend
	TLegend* legend4 = new TLegend(0.5, 0.76, 0.86, 0.94);
	legend4->AddEntry(hist4a, "Data, ^{10}C E*=5.18 MeV gate", "ep");
	legend4->AddEntry(hist4b, "Sideband bkgd. (scaled)", "l");
	legend4->Draw();

	mycan4->Print("Canvas_4.eps", "eps");

	// Find or make fifth canvas
	TCanvas* mycan5 = (TCanvas*)gROOT->FindObjectAny("mycan5");
	if (!mycan5) mycan5 = new TCanvas("mycan5", "", canw, canh);

	xtitle = "E_{rel} 4p+2#alpha (MeV, ^{10}C E*=5.18 MeV, bkgd. sub.)";
	hist4a->Add(hist4b, -1);
	hist4a->GetXaxis()->SetTitle(xtitle.c_str());
	hist4a->GetYaxis()->SetTitle("");
	hist4a->Draw("EX0P");

	mycan5->Print("Canvas_5.eps", "eps");
}



