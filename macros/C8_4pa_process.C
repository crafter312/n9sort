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
 * like `root -l -q "../macros/C8_4pa_process.C+"`. This way, the various
 * library and ROOT dictionary files required for this to work are in the
 * current directory and can be loaded automatically. Compile using the +
 * at the end of the file name for ACLiC, otherwise it will error if used. 
 */

using namespace std;

void C8_4pa_process() {

	// Read in file
	TFile *file = TFile::Open("/data4/N9/mnt/analysis/e25001/rootout/sort_run16-54_noneighbors_hasCsITDC_SiFGate.root");
	if (!file || file->IsZombie()) return;

	// Get TTree from file
	TTree *tree = (TTree*)file->Get("InvMass/8C/t_C8_4pa");
	if (!tree) {
		file->Close();
		return;
	}

	// Classes for re-processing input data
	solution frag_sols[5]; // 4 p 1 a
	correl2 correl;
	
	// Per-telescope counters
	vector<size_t> multsP(4, 0);                // per-event per-telescope proton multiplicity
	vector<size_t> multsA(4, 0);                // per-event per-telescope alpha multiplicity
	vector<vector<int>> pIDs(4, vector<int>()); // track the CsI crystal IDs for protons in each telescope
	vector<vector<int>> aIDs(4, vector<int>()); // track the CsI crystal IDs for alphas in each telescope
	
	// Total counters
	vector<size_t> singleMults(4, 0); // per-telescope total count of events with only one particle in this telescope
	vector<size_t> multiMults(4, 0);  // per-telescope total count of events with multiple particles in this telescope

	float const Q8C = mass_8C - (4*mass_p) - mass_alpha;
	
	// ROOT output
	TFile* ofile = new TFile("/data4/N9/mnt/analysis/e25001/rootout/C8_4pa_processed_run16-54_noneighbors_hasCsITDC_SiFGate.root", "RECREATE");
	ofile->cd();
	TH2I* p2_csicombos = new TH2I("p2_csicombos", "p2_csicombos", 7, 0, 7, 7, 0, 7);
	TDirectoryFile* dirInvMass = new TDirectoryFile("InvMass", "InvMass");
	TDirectory* dir8C = dirInvMass->mkdir("8C", "8C");
	dir8C->cd();
	TH1I* Be6_subevents = new TH1I("Be6_subevents", "Be6_subevents", 200, 0, 16);
	correl.zeroMask();
	correl.proton.mask[0] = 1;
	correl.proton.mask[1] = 1;
	correl.proton.mask[2] = 1;
	correl.proton.mask[3] = 1;
	correl.alpha.mask[0]  = 1;
	wood C8_4pa(correl, "t_C8_4pa", dir8C, false);
	
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
		if (frags.size() != 5)
			throw runtime_error("8C->4p+a event must have only 5 fragments!");
		for (size_t i = 0; i < frags.size(); i++) {
			wood::GobbiOut& frag = frags[i];
		
			// Particle types are ordered, so first four are protons and last is an alpha
			if (i < 4) {
				multsP[frag.itele]++;
				pIDs[frag.itele].push_back(frag.iCsI);
			}
			else {
				multsA[frag.itele]++;
				aIDs[frag.itele].push_back(frag.iCsI);
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
			sol->iCsI      = frag.iCsI;
			sol->ifront    = frag.ifront;
			sol->iback     = frag.iback;
			sol->CsITime   = frag.time;
			sol->denergyR  = frag.denergy_R;
			sol->energyR   = frag.energy_p_R;
			sol->qdc       = frag.qdc;

			// Set mass value (total mass in MeV)
			// 5th and last fragment should be alpha (heaviest),
			// and everything else should be protons
			sol->iZ = 1;
			sol->iA = 1;
			sol->mass = Mass_p;
			if (i == 4) {
				sol->iZ = 2;
				sol->iA = 4;
				sol->mass = Mass_alpha;
			}

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

		// Loop through all combinations of 2p + a and perform 6Be reconstruction
		double num6Be = 0.;
		for (size_t i = 0; i < 3; i++) {
			for (size_t j = i + 1; j < 4; j++) {
				correl.zeroMask();
				correl.alpha.mask[0]  = 1;
				correl.proton.mask[i] = 1;
				correl.proton.mask[j] = 1;
				correl.makeArray(1);
				float Erel_6Be = correl.findErel();
				Be6_subevents->Fill(Erel_6Be);
				num6Be += (double)((Erel_6Be > 0.) && (Erel_6Be < 1.8));
			}
		}
		if (num6Be == 0.) continue;

		// At this point, at least one combination of 2p + a is 6Be, perform the
		// 8C reconstruction (although there is more background inside the gate
		// here than in other cases)
		correl.zeroMask();
		correl.proton.mask[0] = 1;
		correl.proton.mask[1] = 1;
		correl.proton.mask[2] = 1;
		correl.proton.mask[3] = 1;
		correl.alpha.mask[0]  = 1;
		correl.makeArray(1);

		float Erel_8C = correl.findErel();
		float Ex = Erel_8C - Q8C;
		float Vcm = correl.velocityCM;
		float thetaCM = correl.thetaCM*rad_to_deg;

		// No mass excess for 9N, so no Q value and no excitation energy
		runnum = *runnumRV;
		correl.makeOutput(1, C8_4pa);
		C8_4pa.Fill(Erel_8C, Ex, Vcm, thetaCM, runnum, 8);
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
