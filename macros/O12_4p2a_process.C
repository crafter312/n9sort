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

	// Radius for 6Be+6Be 2D gate
	const double Be6Rad = 0.292;

	// Read in file
	TFile *file = TFile::Open("/data4/N9/mnt/analysis/e25001/rootout/sort_all_noneighbors_hasCsITDCQDC_SiFBGates_CsIPSD_CsIrecal_p+3MeV.root");
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

	double num6Be6Be{0.};
	double num6Be6BeSideband{0.};
	double num8C{0.};
	double num8Csideband{0.};
	double num8BeGated9B{0.};
	double num8Be9BGated10C{0.};
	double num8Be9BGated10Csideband{0.};
	bool has8Be{false};
	double avgSubeventsInGate{0.};
	float const Q6Be = mass_6Be - (mass_alpha + 2.*mass_p);
	
	//// ROOT output
	TFile* ofile = new TFile("/data4/N9/mnt/analysis/e25001/rootout/O12_4p2a_processed_all_noneighbors_hasCsITDCQDC_SiFBGates_CsIPSD_CsIrecal_p+3MeV.root", "RECREATE");
	ofile->cd();

	// Make trees
	TH2I* p2_csicombos = new TH2I("p2_csicombos", "p2_csicombos", 7, 0, 7, 7, 0, 7);
	TH2I* a2_csicombos = new TH2I("a2_csicombos", "a2_csicombos", 7, 0, 7, 7, 0, 7);
	TDirectoryFile* dirInvMass = new TDirectoryFile("InvMass", "InvMass");
	TDirectory* dir12O = dirInvMass->mkdir("12O", "12O");
	dir12O->cd();
	TH2I* Be6_subevents = new TH2I("Be6_subevents_8BeVeto", "Be6_subevents_8BeVeto", 200, 0, 10, 200, 0, 10);
	TH2I* Be6_subevents_zoomed = new TH2I("Be6_subevents_8BeVeto_zoomed", "Be6_subevents_8BeVeto_zoomed", 50, 0, 2.5, 50, 0, 2.5);
	TH1I* Be6_subevents_gate_x_proj = new TH1I("Be6_subevents_8BeVeto_x_proj", "Be6_subevents_8BeVeto_proj", 100, 0, 5);
	TH1I* Be6_subevents_ExRad = new TH1I("Be6_subevents_8BeVeto_ExRad", "Be6_subevents_8BeVeto_ExRad", 200, 0, 10);
	TH1I* Be6_subevents_gate_yx_proj = new TH1I("Be6_subevents_8BeVeto_xy_proj", "Be6_subevents_8BeVeto_xy_proj", 200, 0, 10);
	TH1I* C8_subevents = new TH1I("C8_subevents_8BeVeto", "C8_subevents_8BeVeto", 200, 0, 16);
	TH1I* Be8_subevents = new TH1I("Be8_subevents", "Be8_subevents", 1000, 0, 10);
	TH1I* B9_subevents = new TH1I("B9_subevents_8Begate", "B9_subevents_8Begate", 200, 0, 12);
	TH1I* C10_subevents_8Be9Bgate = new TH1I("C10_subevents_8Be9Bgate", "C10_subevents_8Be9Bgate", 200, 0, 12);
	TH1I* C10_subevents_8Begate_9Bveto = new TH1I("C10_subevents_8Begate_9Bveto", "C10_subevents_8Begate_9Bveto", 200, 0, 12);
	TH1I* O12_Erel_8Begate = new TH1I("O12_Erel_8Begate", "O12_Erel_8Begate", 200, 4, 15);
	TH1I* O12_Erel_8Be9B10Cgates = new TH1I("O12_Erel_8Be9B10Cgates", "O12_Erel_8Be9B10Cgates", 200, 4, 15);
	TH1I* O12_Erel_8Be9Bgated10Csideband = new TH1I("O12_Erel_8Be9Bgated10Csideband", "O12_Erel_8Be9Bgated10Csideband", 200, 4, 15);
	TH1I* O12_Erel_6Be6Begate = new TH1I("O12_Erel_6Be6Begate_8BeVeto", "O12_Erel_6Be6Begate_8BeVeto", 200, 4, 15);
	TH1I* O12_Erel_6Be6Besideband = new TH1I("O12_Erel_6Be6Besideband_8BeVeto", "O12_Erel_6Be6Besideband_8BeVeto", 200, 4, 15);
	TH1I* O12_Erel_8Cgate = new TH1I("O12_Erel_8Cgate_8BeVeto", "O12_Erel_8Cgate_8BeVeto", 200, 4, 15);
	TH1D* O12_Erel_8Csideband = new TH1D("O12_Erel_8Csideband_8BeVeto", "O12_Erel_8Csideband_8BeVeto", 200, 4, 15);

	// Wood class for tree output
	correl.zeroMask();
	correl.proton.mask[0] = 1;
	correl.proton.mask[1] = 1;
	correl.proton.mask[2] = 1;
	correl.proton.mask[3] = 1;
	correl.alpha.mask[0]  = 1;
	correl.alpha.mask[1]  = 1;
	wood O12_4p2a(correl, "t_O12_4p2a", dir12O, false);
	{
		TTree* tout = O12_4p2a.getTree();
		tout->Branch("num6Be6Be", &num6Be6Be);
		tout->Branch("num6Be6BeSideband", &num6Be6BeSideband);
		tout->Branch("num8C", &num8C);
		tout->Branch("num8Csideband", &num8Csideband);
		tout->Branch("num8BeGated9B", &num8BeGated9B);
		tout->Branch("num8Be9BGated10C", &num8Be9BGated10C);
		tout->Branch("num8Be9BGated10Csideband", &num8Be9BGated10Csideband);
		tout->Branch("has8Be", &has8Be);
	}

	// TTreeReader loop
	TTreeReader reader(tree);
	TTreeReaderValue<deque<wood::GobbiOut>> fragsRV(reader, "gobbiFrags");
	TTreeReaderValue<double> Erel_12ORV(reader, "Erel");
	TTreeReaderValue<double> ExRV(reader, "Ex");
	TTreeReaderValue<double> VcmRV(reader, "Vcm");
	TTreeReaderValue<double> thetaCMRV(reader, "thetaCM");
	TTreeReaderValue<int> runnumRV(reader, "runnum");
	TTreeReaderValue<int> beamZRV(reader, "beamZ");
	deque<wood::GobbiOut> frags;
	double Erel_12O, Ex, Vcm, thetaCM;
	int runnum, beamZ;
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

		// First, reset the tree variables
		num6Be6Be                = 0.;
		num6Be6BeSideband        = 0.;
		num8C                    = 0.;
		num8Csideband            = 0.;
		num8BeGated9B            = 0.;
		num8Be9BGated10C         = 0.;
		num8Be9BGated10Csideband = 0.;
		has8Be                   = false;

		// Then, retrieve other inputs
		Erel_12O = *Erel_12ORV;
		Ex       = *ExRV;
		Vcm      = *VcmRV;
		thetaCM  = *thetaCMRV;
		runnum   = *runnumRV;
		beamZ    = *beamZRV;

		// Check for 8Be as marker of 10C decay pathway (there is only one, should be very clear with no background)
		// 12O -> 10C + 2p -> 8Be + 4p -> 2a + 4p
		// 10C may also have a 1p decay path, I'm not sure (ends up in 8Be either way)
		correl.zeroMask();
		correl.alpha.mask[0] = 1;
		correl.alpha.mask[1] = 1;
		correl.makeArray(1);
		float Erel_8Be = correl.findErel();
		Be8_subevents->Fill(Erel_8Be);
		has8Be = (Erel_8Be > 0.f) && (Erel_8Be < .2f);

		// If contains 8Be subevent, reconstruct all 9B and 10C subevents.
		// Then reconstruct 12O, fill dedicated histogram, and continue the
		// event loop
		if (has8Be) {

			// First 9B subevents
			for (size_t i = 0; i < 4; i++) {
				correl.zeroMask();
				correl.proton.mask[i] = 1;
				correl.alpha.mask[0]  = 1;
				correl.alpha.mask[1]  = 1;
				correl.makeArray(1);
				float Erel_9B = correl.findErel();
				B9_subevents->Fill(Erel_9B);
				if ((Erel_9B > 0.f) && (Erel_9B < .5f))
					num8BeGated9B += 1.;
			}

			// Then reconstruct 10C subevents
			correl.zeroMask();
			correl.alpha.mask[0] = 1;
			correl.alpha.mask[1] = 1;
			for (size_t i = 0; i < 4; i++) {
				for (size_t j = i + 1; j < 4; j++) {
					correl.proton.mask[i] = 1;
					correl.proton.mask[j] = 1;

					// Zero the complementary indices
					// i.e. if (i, j) = (0, 3), then (k, l) = (1, 2)
					size_t mask = 0xF ^ ((1u << i) | (1u << j));
					size_t k = __builtin_ctz(mask);
					size_t l = __builtin_ctz(mask & (mask - 1));
					correl.proton.mask[k] = 0;
					correl.proton.mask[l] = 0;

					correl.makeArray(1);
					float Erel_10C = correl.findErel();
					if (num8BeGated9B > 0.) {
						num8Be9BGated10C += (double)((Erel_10C > 0.f) && (Erel_10C < 2.f));
						num8Be9BGated10Csideband += (double)((Erel_10C >= 2.f) && (Erel_10C < 3.f));
						C10_subevents_8Be9Bgate->Fill(Erel_10C);
					}
					else if (num8BeGated9B == 0.)
						C10_subevents_8Begate_9Bveto->Fill(Erel_10C); 
				}
			}

			O12_Erel_8Begate->Fill(Erel_12O);
		}
		if (num8Be9BGated10C > 0.)
			O12_Erel_8Be9B10Cgates->Fill(Erel_12O);
		else if (num8Be9BGated10Csideband > 0.)
			O12_Erel_8Be9Bgated10Csideband->Fill(Erel_12O);

		// Loop through all combinations of (2p + a) + (2p + a) and perform 6Be
		// reconstruction for each. To do this, only have to pick one alpha and
		// look at all pairs of protons that could go with that alpha. For each
		// case, remaining three fragments form the other 6Be.
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
				Be6_subevents_ExRad->Fill(radius);
				//num6Be6Be += (double)cutg->IsInside(Erel_6Be_1, Erel_6Be_2);
				num6Be6Be += (double)(radius < Be6Rad);
				num6Be6BeSideband += (double)((radius >= Be6Rad) && (radius < (2.*Be6Rad)));
				bool insideDiagonalBand = (abs(Erel_6Be_2 - Erel_6Be_1) / sqrt(2.)) <= Be6Rad;
				double u = (Erel_6Be_1 + Erel_6Be_2) / sqrt(2.); // Along the y=x line
				if (insideDiagonalBand) {
					Be6_subevents_gate_x_proj->Fill(Erel_6Be_1);
					Be6_subevents_gate_yx_proj->Fill(u);
				}
			}
		}
		if (num6Be6Be > 0.)
			O12_Erel_6Be6Begate->Fill(Erel_12O);
		else if (num6Be6BeSideband > 0.)
			O12_Erel_6Be6Besideband->Fill(Erel_12O);

		// Then check all combinations of 8C + a (there are only two)
		int current_alpha_status[2] = {0, 0}; // 0: neither, 1: peak, 2: sideband

		correl.zeroMask();
		correl.proton.mask[0] = 1;
		correl.proton.mask[1] = 1;
		correl.proton.mask[2] = 1;
		correl.proton.mask[3] = 1;
		correl.alpha.mask[0]  = 1;
		correl.makeArray(1);
		float Erel_8C = correl.findErel();
		C8_subevents->Fill(Erel_8C);
		if ((Erel_8C > 2.) && (Erel_8C < 4.)) {
			num8C += 1.;
			current_alpha_status[0] = 1;
		}
		else if ((Erel_8C >= 4.) && (Erel_8C < 5.)) {
			num8Csideband += 1.0;
			current_alpha_status[0] = 2;
		}

		correl.alpha.mask[0]  = 0;
		correl.alpha.mask[1]  = 1;
		correl.makeArray(1);
		Erel_8C = correl.findErel();
		C8_subevents->Fill(Erel_8C);
		if ((Erel_8C > 2.) && (Erel_8C < 4.)) {
			num8C += 1.;
			current_alpha_status[1] = 1;
		}
		else if ((Erel_8C >= 4.) && (Erel_8C < 5.)) {
			num8Csideband += 1.0;
			current_alpha_status[1] = 2;
		}

		if (num8C > 0.) O12_Erel_8Cgate->Fill(Erel_12O);
		if ((num8C <= 0.) && (num8Csideband > 0.)) O12_Erel_8Csideband->Fill(Erel_12O);

		// Handle final output
		correl.proton.mask[0] = 1;
		correl.proton.mask[1] = 1;
		correl.proton.mask[2] = 1;
		correl.proton.mask[3] = 1;
		correl.alpha.mask[0]  = 1;
		correl.alpha.mask[1]  = 1;
		correl.makeArray(1);
		correl.makeOutput(1, O12_4p2a);
		O12_4p2a.Fill(Erel_12O, Ex, Vcm, thetaCM, runnum, beamZ);
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
