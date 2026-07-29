#include <deque>
#include <stdexcept>
#include <vector>

#include <TFile.h>
#include <TH2I.h>
#include <TTree.h>
#include <TTreeReader.h>

#define rel 1
R__LOAD_LIBRARY(libn9sort.so)

#include "../src/constants.h"
#include "../src/correl2.h"
#include "../src/solution.h"
#include "../src/wood.h"

/**
 * Run this from inside the main project's build directory with something
 * like `root -l -q "../macros/N9_5pa_process.C+"`. This way, the various
 * library and ROOT dictionary files required for this to work are in the
 * current directory and can be loaded automatically. Compile using the +
 * at the end of the file name for ACLiC, otherwise it will error if used. 
 */

using namespace std;

void N9_5pa_process() {

	// Read in file
	TFile *file = TFile::Open("/data4/N9/mnt/analysis/e25001/rootout/sort_all_noneighbors.root");
	if (!file || file->IsZombie()) return;

	// Get TTree from file
	TTree *tree = (TTree*)file->Get("InvMass/9N/t_N9_5pa");
	if (!tree) {
		file->Close();
		return;
	}

	// Classes for re-processing input data
	solution frag_sols[6]; // 5 p 1 a
	correl2 correl;
	
	// Per-telescope counters
	vector<size_t> multsP(4, 0);                // per-event per-telescope proton multiplicity
	vector<size_t> multsA(4, 0);                // per-event per-telescope alpha multiplicity
	vector<vector<int>> pIDs(4, vector<int>()); // track the CsI crystal IDs for protons in each telescope
	vector<vector<int>> aIDs(4, vector<int>()); // track the CsI crystal IDs for alphas in each telescope
	
	// Total counters
	vector<size_t> singleMults(4, 0); // per-telescope total count of events with only one particle in this telescope
	vector<size_t> multiMults(4, 0);  // per-telescope total count of events with multiple particles in this telescope
	
	// ROOT output
	TFile* ofile = new TFile("/data4/N9/mnt/analysis/e25001/rootout/N9_5pa_processed_all_noneighbors.root", "RECREATE");
	ofile->cd();
	TH2I* p2_csicombos = new TH2I("p2_csicombos", "p2_csicombos", 7, 0, 7, 7, 0, 7);
	TDirectoryFile* dirInvMass = new TDirectoryFile("InvMass", "InvMass");
	TDirectory* dir9N = dirInvMass->mkdir("9N", "9N");
	correl.zeroMask();
	correl.proton.mask[0] = 1;
	correl.proton.mask[1] = 1;
	correl.proton.mask[2] = 1;
	correl.proton.mask[3] = 1;
	correl.proton.mask[4] = 1;
	correl.alpha.mask[0]  = 1;
	wood N9_5pa(correl, "t_N9_5pa", dir9N, false);
	
	// TTreeReader loop
	TTreeReader reader(tree);
	TTreeReaderValue<deque<wood::GobbiOut>> fragsRV(reader, "gobbiFrags");
	TTreeReaderValue<int> runnumRV(reader, "runnum");
	deque<wood::GobbiOut> frags;
	int runnum;
	while (reader.Next()) {
		for (size_t i = 0; i < 4; i++) {
			multsP[i] = 0;
			multsA[i] = 0;
			pIDs[i].clear();
			aIDs[i].clear();
		}
		correl.reset();
		
		// Determine the proton and alpha multiplicities in each telescope
		frags = *fragsRV;
		if (frags.size() != 6)
			throw runtime_error("9N->5p+a event must have only 6 fragments!");
		for (size_t i = 0; i < frags.size(); i++) {
			wood::GobbiOut& frag = frags[i];
		
			// Particle types are ordered, so first four are protons and last two are alphas
			if (i < 4) {
				multsP[frag.itele]++;
				pIDs[frag.itele].push_back(frag.id);
			}
			else {
				multsA[frag.itele]++;
				aIDs[frag.itele].push_back(frag.id);
			}

			// Transfer input data back into solutions for re-processing
			solution* sol = &frag_sols[i];
			sol->reset();
			sol->Mvect[0]  = frag.M[0];
			sol->Mvect[1]  = frag.M[1];
			sol->Mvect[2]  = frag.M[2];
			sol->energyTot = frag.et;
			sol->Ekin      = frag.energy_p;
			sol->theta     = frag.theta;
			sol->phi       = frag.phi;
			sol->itele     = frag.itele;
			sol->iCsI      = frag.id;
			sol->ifront    = frag.ifront;
			sol->iback     = frag.iback;
			sol->CsITime   = frag.time;
			sol->denergyR  = frag.denergy_R;
			sol->energyR   = frag.energy_p_R;

			// Set mass value (total mass in MeV)
			// 6th and last fragment should be alpha (heaviest),
			// and everything else should be protons
			sol->iZ = 1;
			sol->iA = 1;
			sol->mass = Mass_p;
			if (i == 5) sol->mass = Mass_alpha;

			// Now to load solutions into correl
			correl.load(sol);
		}
		
		// Add the per-telescope single and multi particle events to the total counters
		for (size_t i = 0; i < 4; i++) {
			if ((multsP[i] + multsA[i]) > 1) multiMults[i]++;
			else singleMults[i]++;
			
			if (multsP[i] == 2) p2_csicombos->Fill(pIDs[i][0], pIDs[i][1]);
		}

		/******** RE-PROCESSING ********/
		// With all the solutions filled and loaded into correl,
		// can now perform calculations
		
		// First, set baseline mask for 9N
		correl.zeroMask();
		correl.proton.mask[0] = 1;
		correl.proton.mask[1] = 1;
		correl.proton.mask[2] = 1;
		correl.proton.mask[3] = 1;
		correl.proton.mask[4] = 1;
		correl.alpha.mask[0]  = 1;

		// Loop through all combinations of 4p + a and perform 8C reconstruction
		bool has8C = false;
		for (size_t i = 0; i < 5; i++) {
			size_t idx = (i > 0) ? i - 1 : 0;
			correl.proton.mask[idx] = 1;
			correl.proton.mask[i] = 0;
			correl.makeArray(1);
			float Erel_8C = correl.findErel();
			has8C = (Erel_8C > 2.) && (Erel_8C < 4.);
		}
		if (!has8C) continue;

		// At this point, at least one combination of 4p + a is 8C, perform the
		// 9N reconstruction
		correl.proton.mask[4] = 1;
		correl.makeArray(1, N9_5pa);

		float Erel_9N = correl.findErel();
		float Vcm = correl.velocityCM;
		float thetaCM = correl.thetaCM*rad_to_deg;
		float cos_thetaH = correl.cos_thetaH;

		// No mass excess for 9N, so no Q value and no excitation energy
		runnum = *runnumRV;
		N9_5pa.Fill(Erel_9N, -1, Vcm, thetaCM, cos_thetaH, runnum, 8);
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
