#include <deque>
#include <stdexcept>
#include <vector>

#include <TCutG.h>
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
 * like `root -l -q "../macros/O12_4p2a_process.C+"`. This way, the various
 * library and ROOT dictionary files required for this to work are in the
 * current directory and can be loaded automatically.
 */

using namespace std;

void O12_4p2a_process() {

	// ROOT macro output for 6Be gate
	vector<double> cutg_vect0{ 1.175373077336973, 0.6156715764594571, 0.5923506805895606, 0.9888059103778015, 1.781716369954282, 2.13152980800273, 2.061567120393041, 1.548507411255317, 1.222014869076766, 1.175373077336973 };
	vector<double> cutg_vect1{ 2.333333293596901, 2.249999959021805, 1.333333278695738, 0.6249999348074162, 0.3749999310821259, 0.9999999403953517, 1.666666616996126, 2.208333291734256, 2.37499996088445, 2.333333293596901 };
	TCutG *cutg = new TCutG("O12_6BeGate", 10, cutg_vect0.data(), cutg_vect1.data());

	// Read in file
	TFile *file = TFile::Open("/data4/N9/mnt/analysis/e25001/rootout/sort_run16-54_noneighbors_hasCsITDC_SiFBGates.root");
	if (!file || file->IsZombie()) return;

	// Get TTree from file
	TTree *tree = (TTree*)file->Get("InvMass/12O/t_O12_4p2a");
	if (!tree) {
		file->Close();
		return;
	}

	// Classes for re-processing input data
	solution frag_sols[6]; // 4 p 2 a
	correl2 correl;
	
	// Per-telescope counters
	vector<size_t> multsP(4, 0);                // per-event per-telescope proton multiplicity
	vector<size_t> multsA(4, 0);                // per-event per-telescope alpha multiplicity
	vector<vector<int>> pIDs(4, vector<int>()); // track the CsI crystal IDs for protons in each telescope
	vector<vector<int>> aIDs(4, vector<int>()); // track the CsI crystal IDs for alphas in each telescope
	
	// Total counters
	vector<size_t> singleMults(4, 0); // per-telescope total count of events with only one particle in this telescope
	vector<size_t> multiMults(4, 0);  // per-telescope total count of events with multiple particles in this telescope

	size_t num12Os{0};
	double avgSubeventsInGate{0.};
	float const Q6Be = mass_6Be - (mass_alpha + 2.*mass_p);
	
	// ROOT output
	TFile* ofile = new TFile("/data4/N9/mnt/analysis/e25001/rootout/O12_4p2a_processed_run16-54_noneighbors_hasCsITDC_SiFBGates.root", "RECREATE");
	ofile->cd();
	TH2I* p2_csicombos = new TH2I("p2_csicombos", "p2_csicombos", 7, 0, 7, 7, 0, 7);
	TH2I* a2_csicombos = new TH2I("a2_csicombos", "a2_csicombos", 7, 0, 7, 7, 0, 7);
	TDirectoryFile* dirInvMass = new TDirectoryFile("InvMass", "InvMass");
	TDirectory* dir12O = dirInvMass->mkdir("12O", "12O");
	dir12O->cd();
	TH2I* Be6_subevents = new TH2I("Be6_subevents", "Be6_subevents", 200, 0, 10, 200, 0, 10);
	TH2I* Be6_subevents_zoomed = new TH2I("Be6_subevents_zoomed", "Be6_subevents_zoomed", 25, 0, 2.5, 25, 0, 2.5);
	TH1I* C8_subevents = new TH1I("C8_subevents", "C8_subevents", 200, 0, 16);
	correl.zeroMask();
	correl.proton.mask[0] = 1;
	correl.proton.mask[1] = 1;
	correl.proton.mask[2] = 1;
	correl.proton.mask[3] = 1;
	correl.alpha.mask[0]  = 1;
	correl.alpha.mask[1]  = 1;
	wood O12_4p2a(correl, "t_O12_4p2a", dir12O, false);

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
			throw runtime_error("12O->4p+2a event must have only 6 fragments!");
		for (size_t i = 0; i < frags.size(); i++) {
			wood::GobbiOut& frag = frags[i];
		
			// Particle types are ordered, so first four are protons and last two are alphas
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
			sol->Ekin      = frag.Ekin;
			sol->theta     = frag.theta;
			sol->phi       = frag.phi;
			sol->itele     = frag.itele;
			sol->iCsI      = frag.iCsI;
			sol->ifront    = frag.ifront;
			sol->iback     = frag.iback;
			sol->CsITime   = frag.CsITime;
			sol->time      = frag.time;
			sol->btime     = frag.btime;
			sol->timediff  = frag.timediff;
			sol->energyR   = frag.energyR;
			sol->denergyR  = frag.denergyR;
			sol->energy    = frag.energy;
			sol->denergy   = frag.denergy;
			sol->qdc       = frag.qdc;

			// Set mass value (total mass in MeV)
			// 5th and 6th fragments (last two) should be alpha (heaviest),
			// and everything else should be protons
			sol->iZ = 1;
			sol->iA = 1;
			sol->mass = Mass_p;
			if ((i == 4) || (i == 5)) {
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
			if (multsA[i] == 2) a2_csicombos->Fill(aIDs[i][0], aIDs[i][1]);
		}

		/******** RE-PROCESSING ********/
		// With all the solutions filled and loaded into correl,
		// can now perform calculations

		// Loop through all combinations of (2p + a) + (2p + a) and perform 6Be
		// reconstruction for each. To do this, only have to pick one alpha and
		// look at all pairs of protons that could go with that alpha. For each
		// case, remaining three fragments form the other 6Be.
		double num12Os{0.};
		for (size_t i = 0; i < 3; i++) {
			for (size_t j = i + 1; j < 4; j++) {
				correl.zeroMask();
				correl.alpha.mask[0]  = 1;
				correl.proton.mask[i] = 1;
				correl.proton.mask[j] = 1;
				correl.makeArray(1);
				float Erel_6Be_1 = correl.findErel();
				float Ex_6Be_1 = Erel_6Be_1 - Q6Be;

				// Calculate second 6Be fragment from remaining fragments
				correl.zeroMask();
				size_t mask = 0xF ^ ((1u << i) | (1u << j));
				size_t k = __builtin_ctz(mask);
				size_t l = __builtin_ctz(mask & (mask - 1));
				correl.alpha.mask[1]  = 1;
				correl.proton.mask[k] = 1;
				correl.proton.mask[l] = 1;
				correl.makeArray(1);
				float Erel_6Be_2 = correl.findErel();
				float Ex_6Be_2 = Erel_6Be_2 - Q6Be;

				Be6_subevents->Fill(Erel_6Be_1, Erel_6Be_2);
				Be6_subevents_zoomed->Fill(Erel_6Be_1, Erel_6Be_2);
				float radius = sqrt((Ex_6Be_1 * Ex_6Be_1) + (Ex_6Be_2 * Ex_6Be_2));
				//num12Os += (double)cutg->IsInside(Erel_6Be_1, Erel_6Be_2);
				num12Os += (double)(radius < 0.292);
			}
		}

		// Then check all combinations of 8C + a (there are only two)
		double num8Cs{0.};
		correl.zeroMask();
		correl.proton.mask[0] = 1;
		correl.proton.mask[1] = 1;
		correl.proton.mask[2] = 1;
		correl.proton.mask[3] = 1;
		correl.alpha.mask[0]  = 1;
		correl.makeArray(1);
		float Erel_8C = correl.findErel();
		C8_subevents->Fill(Erel_8C);
		num8Cs += (double)((Erel_8C > 2.) && (Erel_8C < 4.));

		correl.alpha.mask[0]  = 0;
		correl.alpha.mask[1]  = 1;
		correl.makeArray(1);
		Erel_8C = correl.findErel();
		C8_subevents->Fill(Erel_8C);
		num8Cs += (double)((Erel_8C > 2.) && (Erel_8C < 4.));

		if ((num12Os == 0.) && (num8Cs == 0.)) continue;

		// At this point, at least one combination of 4p + 2a is either
		// 6Be+6Be or 8C+a, perform the 12O reconstruction
		float const Q12O = mass_12O - (4*mass_p) - (2*mass_alpha);
		correl.alpha.mask[0]  = 1;
		correl.alpha.mask[1]  = 1;
		correl.makeArray(1);

		float Erel_12O = correl.findErel();
		float Ex = Erel_12O - Q12O;
		float Vcm = correl.velocityCM;
		float thetaCM = correl.thetaCM*rad_to_deg;

		runnum = *runnumRV;
		correl.makeOutput(1, O12_4p2a);
		O12_4p2a.Fill(Erel_12O, Ex, Vcm, thetaCM, runnum, 8);
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
