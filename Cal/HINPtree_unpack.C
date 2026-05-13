//**********************
//Reformat E and dE calibration data to be read by asicPeakSearch.C
//Just want to output histograms for every strip
//There are a few things that need be be specified, i.e. DEs vs Es
//Change at commented locations
//**********************
#include <TFile.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#define BOARD_COUNT 12
#define CHAN_COUNT 32
#define PSD_CHIP_COUNT 24
#define PSD_CHAN_COUNT 8 

using namespace std;

vector<string> generate_column_names_hinp(const string& parname, string silname) {
	vector<string> columns;
	string b, c;

  int startboard = 1;
  int endboard = 12;

  if (silname == "DE") {
    startboard = 9;
  }
  else if (silname == "E") {
    endboard = 8;
  }
  else {
    cout << "Specify DE or E silicons" << endl;
    abort();
  }

	for (int board = startboard; board <= endboard; board++) {
		for (int chan = 0; chan < CHAN_COUNT; chan++) {
			b = (board < 10 ? "0" : "") + to_string(board);
			c = (chan < 10 ? "0" : "") + to_string(chan);
			columns.push_back("SpecTcl_hinpa_mb1_" + parname + "_" + b + "." + c);
		}
	}
	return columns;
}

void HINPtree_unpack() {

  // Get number of entries from input file
	string iprefix = "run-5";
	string path = "~/DAQ/SpecTclHira/";
	size_t numentries;

	TFile *file = TFile::Open((path + iprefix + ".root").c_str());
	if (!file || file->IsZombie()) {
		std::cerr << "Error opening file!" << std::endl;
		return;
	}

	// Check if tree exists in the file
	TTree *tree = (TTree*)file->Get("t");
	if (!tree) {
		std::cerr << "Tree 't' not found in the file!" << std::endl;
		file->Close();
		return;
	}
	numentries = tree->GetEntries();

  string ifname = path + iprefix + ".root";	
	cout << "Processing input file " << ifname << "..." << endl;
	cout << "Total entries in tree: " << numentries << endl;
  
  //Set second variable to "E" or "DE"
	vector<string> e_columns = generate_column_names_hinp("e","E");

  //Set filename depending on DE or E
  TFile * outfile = new TFile("HINPEs.root", "RECREATE");
  outfile->cd();

  for (int i=0;i<e_columns.size();i++) {
    cout << e_columns[i] << endl;
    string histparams = e_columns[i] + ">>h(4096,0,8192)";
    tree->Draw(histparams.c_str());
    auto h1 = (TH1F*)gPad->GetPrimitive("h"); 
    h1->Write(e_columns[i].c_str());
  }

  outfile->Close();
}
