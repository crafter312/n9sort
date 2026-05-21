#include <deque>
#include <vector>

#include <TFile.h>
#include <TH2I.h>
#include <TTree.h>
#include <TTreeReader.h>

#include "../src/wood.h"

using namespace std;

void O12_4p2a_process() {

	// Read in file
	TFile *file = TFile::Open("/mnt/analysis/e25001/rootout/sort_run16-54.root");
	if (!file || file->IsZombie()) return;

	// Get TTree from file
	TTree *tree = (TTree*)file->Get("InvMass/12O/t_O12_4p2a");
	if (!tree) {
		file->Close();
		return;
	}
	
	// Per-telescope counters
	vector<size_t> multsP(4, 0);                // per-event per-telescope proton multiplicity
	vector<size_t> multsA(4, 0);                // per-event per-telescope alpha multiplicity
	vector<vector<int>> pIDs(4, vector<int>()); // track the CsI crystal IDs for protons in each telescope
	vector<vector<int>> aIDs(4, vector<int>()); // track the CsI crystal IDs for alphas in each telescope
	
	// Total counters
	vector<size_t> singleMults(4, 0); // per-telescope total count of events with only one particle in this telescope
	vector<size_t> multiMults(4, 0);  // per-telescope total count of events with multiple particles in this telescope
	
	// ROOT output
	TFile* ofile = new TFile("/mnt/analysis/e25001/rootout/O12_4p2a_processed.root", "RECREATE");
	ofile->cd();
	TH2I* p2_csicombos = new TH2I("p2_csicombos", "p2_csicombos", 7, 0, 7, 7, 0, 7);
	TH2I* a2_csicombos = new TH2I("a2_csicombos", "a2_csicombos", 7, 0, 7, 7, 0, 7);
	
	// TTreeReader loop
	TTreeReader reader(tree);
	TTreeReaderValue<deque<wood::GobbiOut>> fragsRV(reader, "gobbiFrags");
	deque<wood::GobbiOut> frags;
	while (reader.Next()) {
		for (size_t i = 0; i < 4; i++) {
			multsP[i] = 0;
			multsA[i] = 0;
			pIDs[i].clear();
			aIDs[i].clear();
		}
		
		// Determine the proton and alpha multiplicities in each telescope
		frags = *fragsRV;
		for (size_t i = 0; i < frags.size(); i++) {
		
			// Particle types are ordered, so first four are protons and last two are alphas
			if (i < 4) {
				multsP[frags[i].itele]++;
				pIDs[frags[i].itele].push_back(frags[i].id);
			}
			else {
				multsA[frags[i].itele]++;
				aIDs[frags[i].itele].push_back(frags[i].id);
			}
		}
		
		// Add the per-telescope single and multi particle events to the total counters
		for (size_t i = 0; i < 4; i++) {
			if ((multsP[i] + multsA[i]) > 1) multiMults[i]++;
			else singleMults[i]++;
			
			if (multsP[i] == 2) p2_csicombos->Fill(pIDs[i][0], pIDs[i][1]);
			if (multsA[i] == 2) a2_csicombos->Fill(aIDs[i][0], aIDs[i][1]);
		}
	}
	
	cout << "=======================================================" << endl;
	for (size_t i = 0; i < 4; i++) {
		cout << "Telescope " << i << " statistics:" << endl;
		cout << "Single particle events: " << singleMults[i] << endl;
		cout << "Multi-particle events: " << multiMults[i] << endl;
		cout << "% Multi-particle events: " << (double)(multiMults[i]) / (double)(singleMults[i] + multiMults[i]) * 100. << endl;
	}
	cout << "=======================================================" << endl;
	
	ofile->Write();
	ofile->Close();
}
