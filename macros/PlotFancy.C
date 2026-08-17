#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TH1.h>
#include <TList.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>

using namespace std;

/**
 * Helper function to retrieve an arbitrary histogram created or
 * updated as a result of a call to `TTree::Draw` on the currently
 * active pad.
 */
TH1* GetLastDrawnHistogram() {
	if (!gPad) {
		cerr << "Error: No active canvas/pad found." << endl;
		return nullptr;
	}

	// The most recently created/drawn histogram will be near the
	// end of the primitive list, according to Google Gemini
	TList* primitives = gPad->GetListOfPrimitives();
	if (!primitives) return nullptr;

	for (int i = primitives->GetSize() - 1; i >= 0; --i) {
		TObject* prim = primitives->At(i);
		if (prim && prim->InheritsFrom(TH1::Class()))
			return (TH1*)prim;
	}

	return nullptr;
}

/**
 * Function to plot from a wood tree for a particular invariant mass
 * channel, but with pre-set fancy visual attributes. Note that the
 * input file is fixed in the code, but the tree name and the draw
 * input arguments must be provided. Specifically, the tree name
 * must include as a prefix the subdirectory in the tree it is found
 * within. For example, for 9N, one would specify "9N/t_N9_5pa" as
 * input. View the desired input file in a TBrowser to see what all
 * of the various subdirectories and trees are called.
 * @param tname The subdirectory and name of the tree to be loaded
 * @param params The first argument of the `TTree::Draw` call
 * @param gate The second argument of the `TTree::Draw` call
 * @param options The third argument of the `TTree::Draw` call
 * @param xtitle The title for the X axis
 * @param ytitle The title for the Y axis
 */
void PlotFancyHist(string tname, string params, string gate, string options, string xtitle, string ytitle, double xtitleoff=1., double ytitleoff=.9) {

	// Macro constants
	string ifdir  = "/data4/N9/mnt/analysis/e25001/rootout/";
	string ifname = "sort_run16-54_noneighbors_hasCsITDC_SiFBGates_CsIrecal_p+3MeV.root";

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
	TTree* tpar = (TTree*)ifile->Get((string("InvMass/") + tname).c_str());
	if (!tpar) {
		cerr << "TTree " << tname << " not found in file or could not be opened" << endl;
		ifile->Close();
		return;
	}

	// Calculate fractional margins (default margin is 0.1 fraction of width/height)
	const size_t canw = 1000;
	const size_t canh = 610;

	TCanvas* mycan = (TCanvas*)gROOT->FindObjectAny("mycan");
	if (!mycan) mycan = new TCanvas("mycan", "", canw, canh);
	mycan->Draw();

	tpar->Draw(params.c_str(), gate.c_str(), options.c_str());

	// Get histogram and set axis titles, other visual attributes
	TH1* hist = GetLastDrawnHistogram();
	if (hist == nullptr) {
		cerr << "Problem retrieving supposedly drawn histogram, exiting..." << endl;
	}

	hist->GetXaxis()->SetTitle(xtitle.c_str());
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->SetTitle(ytitle.c_str());
	hist->GetYaxis()->CenterTitle();

	mycan->Print("Canvas_1.eps", "eps");
}
