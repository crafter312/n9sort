#include <deque>
#include <iostream>
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
#include "../src/wood.h"

/**
 * Run this from inside the main project's build directory with something
 * like `root -l -q "../macros/decay_constant.C+"`. This way, the various
 * library and ROOT dictionary files required for this to work are in the
 * current directory and can be loaded automatically. Compile using the +
 * at the end of the file name for ACLiC, otherwise it will error if used. 
 */

using namespace std;

void decay_constant() {

	// Read in file
	TFile *file = TFile::Open("/data4/N9/mnt/analysis/e25001/rootout/sort_all_noneighbors_hasCsITDCQDC_SiFBGates_CsIPSD_CsIrecal.root");
	if (!file || file->IsZombie()) return;

	// Get TTree from file (6Be just because I had to pick)
	TTree *tree = (TTree*)file->Get("InvMass/6Be/t_Be6_ppa");
	if (!tree) {
		file->Close();
		return;
	}

	// Histogram ranges

	// For ADC gate 975 to 1025
	const double adc_min   = 975.;
	const double adc_max   = 1025.;
	const double adc_p_min = 500.;
	const double adc_p_max = 1000.;
	const double adc_a_min = 250.;
	const double adc_a_max = 750.;

	// For ADC gate 390 to 410
	//const double adc_min   = 390.;
	//const double adc_max   = 410.;
	//const double adc_p_min = 0.;
	//const double adc_p_max = 500.;
	//const double adc_a_min = 0.;
	//const double adc_a_max = 250.;

	// Extracted decay constants
	const double lambda_p = 0.5 * (0.00031785 + 0.000290268);
	const double lambda_a = 0.000407085;
	
	// ROOT output
	TFile* ofile = new TFile("/data4/N9/mnt/analysis/e25001/rootout/decay_constant_tele0csi0_pa.root", "RECREATE");
	ofile->cd();
	TH2I* p_tail_shape = new TH2I("p_tail_shape", "p_tail_shape", 250, -500, 500, 125, adc_p_min, adc_p_max);
	TH2I* a_tail_shape = new TH2I("a_tail_shape", "a_tail_shape", 250, -500, 500, 125, adc_a_min, adc_a_max);
	TH2I* p_PSD        = new TH2I("p_PSD", "", 1024, 0, 4096, 1024, 0, 4096);
	TH2I* a_PSD        = new TH2I("a_PSD", "", 1024, 0, 4096, 1024, 0, 4096);
	TH2I* p_corr_PSD   = new TH2I("p_corr_PSD", "", 1024, 0, 4096, 1024, 0, 4096);
	TH2I* a_corr_PSD   = new TH2I("a_corr_PSD", "", 1024, 0, 4096, 1024, 0, 4096);

	// Set PSD gates
	vector<double> gate_tele0csi0_pdt_vect0{ 189.8357211034637, 904.0386302966413, 1581.140089661602, 2316.985511254573, 2623.072472337363, 2570.512085080722, 2044.908212514315, 1544.038639833385, 894.7632678395871, 347.5168828733861, 115.6328214470296, 189.8357211034637 };
	vector<double> gate_tele0csi0_pdt_vect1{ 43.81095353244564, 586.5440395252908, 1161.971889734572, 1835.484032593163, 2116.659004854517, 2345.522354369572, 1796.250315533439, 1325.445710816754, 795.7905305104841, 311.9080201072247, 69.96676490559457, 43.81095353244564 };
	TCutG* gate_tele0csi0_pdt = new TCutG("gate_tele0csi0_pdt", 12, gate_tele0csi0_pdt_vect0.data(), gate_tele0csi0_pdt_vect1.data());
	vector<double> gate_tele0csi0_He3a_vect0{ 925.6811426964346, 1683.169076689199, 2431.381648224909, 3646.454130099016, 3686.647367412918, 2743.652184279068, 2202.589374284237, 1547.130427319069, 1064.811579552248, 316.5990080165385, 233.1207459030502, 369.1593952731793, 873.1207554397938, 910.2222052680108, 925.6811426964346 };
	vector<double> gate_tele0csi0_He3a_vect1{ 560.3882281521414, 1161.971889734572, 1770.09450416029, 2783.632194869821, 2672.469996533936, 1842.022985436451, 1397.374192092915, 867.7190117866439, 527.693463935705, 50.34990637573264, 50.34990637573264, 154.9731518683293, 514.6155582491306, 547.310322465567, 560.3882281521414 };
	TCutG* gate_tele0csi0_He3a = new TCutG("gate_tele0csi0_He3a", 15, gate_tele0csi0_He3a_vect0.data(), gate_tele0csi0_He3a_vect1.data());
	
	// TTreeReader loop
	TTreeReader reader(tree);
	TTreeReaderValue<deque<wood::GobbiOut>> fragsRV(reader, "gobbiFrags");
	TTreeReaderValue<int> runnumRV(reader, "runnum");
	deque<wood::GobbiOut> frags;
	int runnum;
	while (reader.Next()) {
		frags = *fragsRV;
		for (size_t i = 0; i < frags.size(); i++) {
			wood::GobbiOut& frag = frags[i];
		
			// First apply shared gates
			if ((frag.iCsI != 0) || (frag.itele != 0)) continue;
			if (gate_tele0csi0_pdt->IsInside(frag.energyR, frag.qdc)) {
				p_PSD->Fill(frag.energyR, frag.qdc);
				p_corr_PSD->Fill(frag.energyR, frag.qdc*exp(lambda_p*frag.CsITime));
				if ((frag.energyR > adc_min) & (frag.energyR < adc_max))
					p_tail_shape->Fill(frag.CsITime, frag.qdc);
			}
			if (gate_tele0csi0_He3a->IsInside(frag.energyR, frag.qdc))
				a_PSD->Fill(frag.energyR, frag.qdc);
				a_corr_PSD->Fill(frag.energyR, frag.qdc*exp(lambda_a*frag.CsITime));
				if ((frag.energyR > adc_min) && (frag.energyR < adc_max))
					a_tail_shape->Fill(frag.CsITime, frag.qdc);
		}
	}
	
	ofile->Write();
	ofile->Close();
}
