/* New Gobbi class for the Gobbi 28 configuration by Henry Webb (h.s.webb@wustl.edu).
 * Created 7 May 2026 as mix of `OldGobbi` class from this code and `gobbi` class
 * used by Johnathan Phillips (j.s.phillips@wustl.edu) for 22Si FRIB experiment.
 */

#include "Gobbi28.h"

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <string>
#include <utility>

#include <TH1I.h>
#include <TH2I.h>

#include "constants.h"
#include "solution.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gobbi28::Gobbi28(Input& in, histo& hist, SortConfig& config) : Targetdist(config.GetTargDist()), TargetThickness(config.GetTargThick()), hinpboards(config.GetHinpboards()), hinpchans(config.GetHinpchans()), Histo(hist), input(in.GetGobbi()), input_qdc(in.GetQDC()), input_adc(in.GetADC()), input_tdc(in.GetTDC()) {

	// Seed TRandom with current system clock
	auto now = chrono::system_clock::now();
	UInt_t tstamp = chrono::duration_cast<chrono::seconds>(now.time_since_epoch()).count();
	Ran = new TRandom(tstamp);

	for (size_t id = 0; id < 4; id++) {
		Telescope[id] = new telescope(TargetThickness, config, true);
		Telescope[id]->init(id, config); // tells Telescope what position it is in
		Telescope[id]->SetTargetDistance(Targetdist);
	}
	NCsI = Telescope[0]->GetNCsI(); // same for every telescope

#ifdef ENABLE_DEBUG
	cout << "Gobbi28::Gobbi28 1" << endl;
#endif

	string calDir = config.GetCalDir();
	FrontEcal = new calibrate(4, hinpchans, calDir + config.GetFrontEcalFile(), 1, false);
	BackEcal = new calibrate(4, hinpchans, calDir + config.GetBackEcalFile(), 1, false);
	CsIEcal = new calibrate(4, NCsI, calDir + config.GetCsIEcalFile(), 1, false);
	FrontTimecal = new calibrate(4, hinpchans, calDir + config.GetFrontTimecalFile(), 1, false);
	BackTimecal = new calibrate(4, hinpchans, calDir + config.GetBackTimecalFile(), 1, false);
	CsITimecal = new calibrate(4, NCsI, calDir + config.GetCsITimecalFile(), 1, false);

#ifdef ENABLE_DEBUG
	cout << "Gobbi28::Gobbi28 2" << endl;
#endif

	// Read in CsI mappings from file
	string incsimapfile = config.GetConfigDir() + config.GetCsIChannelMapFile();
	ifstream incsimap(incsimapfile);
	if (incsimap.fail()) throw invalid_argument(string(BOLDRED) + string("CsI channel map file ") + incsimapfile + string(" does not exist or failed to open") + string(RESET));
	
#ifdef ENABLE_DEBUG
	else cout << GREEN << "CsI channel map file " << incsimapfile << " opened" << RESET << endl;
#endif

#ifdef ENABLE_DEBUG
	cout << "Gobbi28::Gobbi28 3" << endl;
#endif

	size_t chan, tel, id;
	size_t maxid = 0;
	string name;
	while (incsimap.good()) {
		incsimap >> chan >> tel >> id;
		
		if ((tel != 0) && (tel != 1) && (tel != 2) && (tel != 3)) throw invalid_argument(string(BOLDRED) + string("CsI telescope # must be 0, 1, 2, or 3") + string(RESET));
		if (telCsImap.find(chan) != telCsImap.end()) continue; // avoid cases of duplicate ADC channel
		telCsImap[chan] = tel;
		idCsImap[chan] = id;
		
#ifdef ENABLE_DEBUG
		cout << "chan " << chan << " tel " << tel << " id " << id << endl;
#endif
		
		if (chan > maxadcchan) maxadcchan = chan;
		if (id > maxid) maxid = id;
		//cout << "MAXID " << maxid << " id " << id << endl;
		//// Make per-CsI crystal histograms
		
		Histo.dir1dCsI_Energy->cd();
		name = "CsI_Energy_" + to_string(chan) + "_R_unmatched";
		Histo.CsI_Energy_R_um[chan] = new TH1I(name.c_str(), "", 1024, 0, 4096);
		name = "CsI_Energy_" + to_string(chan) + "_cal_unmatched";
		Histo.CsI_Energy_cal_um[chan] = new TH1I(name.c_str(), "", 512, 0, 200);
		
		Histo.dir1dCsI_QDC->cd();
		name = "CsI_QDC_" + to_string(chan) + "_um";
		Histo.CsI_QDC_um[chan] = new TH1I(name.c_str(), "", 1024, 0, 4096);
		name = "CsI_QDC_" + to_string(chan) + "_matched";
		Histo.CsI_QDC_matched[chan] = new TH1I(name.c_str(), "", 1024, 0, 4096);
		
		Histo.dir1dCsI_Time->cd();
		name = "CsI_Time_" + to_string(chan) + "_um";
		Histo.CsI_Time_um[tdcstart + chan] = new TH1I(name.c_str(), "", 1000, -500, 500);
		
		Histo.dirPSD->cd();
		name = "CsIonly_PSD_" + to_string(chan);
		Histo.CsIonly_PSD[chan] = new TH2I(name.c_str(), "", 1024, 0, 4096, 1024, 0, 4096);
	}
	
#ifdef ENABLE_DEBUG
	cout << "Gobbi28::Gobbi28 4" << endl;
#endif
	
	// Make CsI histograms indexed by telescope and per-telescope ID, instead of ADC channel
	nTelCsIs = maxid + 1;
	for (size_t i = 0; i < 4; i++) {
		Histo.CsI_Energy_R[i].resize(nTelCsIs);
		Histo.CsI_Energy_R_center[i].resize(nTelCsIs);
		Histo.CsI_Energy_pcal[i].resize(nTelCsIs);
		Histo.CsI_Energy_pcal_center[i].resize(nTelCsIs);
		Histo.DEE_CsI[i].resize(nTelCsIs);
		Histo.DEE_CsI_sitgate[i].resize(nTelCsIs);
		Histo.DEE_CsI_csitgate[i].resize(nTelCsIs);
		Histo.DEE_CsI_BackE[i].resize(nTelCsIs);
		Histo.DEE_CsI_fronteven[i].resize(nTelCsIs);
		Histo.DEE_CsI_frontodd[i].resize(nTelCsIs);
		Histo.CsI_Time_matched[i].resize(nTelCsIs);
		Histo.CsI_xyhitmap[i].resize(nTelCsIs);
		
#ifdef ENABLE_DEBUG
		cout << "Gobbi28::Gobbi28 4a" << endl;
#endif
		
		for (size_t j = 0; j < nTelCsIs; j++) {
			Histo.dir1dCsI_Energy->cd();
			name = "CsI_Energy_" + to_string(i) + "_" + to_string(j); // i is telescope, j is CsI ID
			Histo.CsI_Energy_R[i][j] = new TH1I(name.c_str(), "", 4096, 0, 4096);
			name = "CsI_Energy_R_center_" + to_string(i) + "_" + to_string(j);
			Histo.CsI_Energy_R_center[i][j] = new TH1I(name.c_str(), "", 4096, 0, 4096);
			name = "CsI_Energy_pcal_" + to_string(i) + "_" + to_string(j); // i is telescope, j is CsI ID
			Histo.CsI_Energy_pcal[i][j] = new TH1I(name.c_str(), "", 4096, 0, 4096);
			name = "CsI_Energy_pcal_center_" + to_string(i) + "_" + to_string(j);
			Histo.CsI_Energy_pcal_center[i][j] = new TH1I(name.c_str(), "", 4096, 0, 4096);
			Histo.dirDEEplots->cd();
			name = "DEE_CsI_" + to_string(i) + "_" + to_string(j);
			Histo.DEE_CsI[i][j] = new TH2I(name.c_str(), "", 512, 0, 4096, 500, 0, 50); // E is x, DE is y
			name = "DEE_CsI_sitgate_" + to_string(i) + "_" + to_string(j);
	  		Histo.DEE_CsI_sitgate[i][j] = new TH2I(name.c_str(), "", 1024, 0, 4096, 500, 0, 80);
	  		name = "DEE_CsI_csitgate_" + to_string(i) + "_" + to_string(j);
	  		Histo.DEE_CsI_csitgate[i][j] = new TH2I(name.c_str(), "", 1024, 0, 4096, 500, 0, 80); 
	  		name = "DEE_CsI_BackE_" + to_string(i) + "_" + to_string(j);
	  		Histo.DEE_CsI_BackE[i][j] = new TH2I(name.c_str(), "", 1024, 0, 4096, 500, 0, 80);
			name = "DEE_CsI_fronteven_" + to_string(i) + "_" + to_string(j);
			Histo.DEE_CsI_fronteven[i][j] = new TH2I(name.c_str(), "", 1024, 0, 4096, 500, 0, 80);
			name = "DEE_CsI_frontodd_" + to_string(i) + "_" + to_string(j);
			Histo.DEE_CsI_frontodd[i][j] = new TH2I(name.c_str(), "", 1024, 0, 4096, 500, 0, 80);
			Histo.dir1dCsI_Time->cd();
			name = "CsI_Time_matched_" + to_string(i) + "_" + to_string(j);
			Histo.CsI_Time_matched[i][j] = new TH1I(name.c_str(), "", 1000, -500, 500);
			Histo.dirhitmaps->cd();
			name = "CsI_xyhitmap_" + to_string(i) + "_" + to_string(j);
			Histo.CsI_xyhitmap[i][j] = new TH2I(name.c_str(), "", 200, -10, 10, 200, -10, 10);
		}
	}

#ifdef ENABLE_DEBUG
	cout << "Gobbi28::Gobbi28 5" << endl;
#endif

	// Make CsI summary plots
	Histo.dirSummary->cd();
	Histo.sumCsIE_R_um       = new TH2I("sumCsIE_R_um", "", maxadcchan, 0, maxadcchan, 1024, 0, 4096);
	Histo.sumCsIE_cal_um     = new TH2I("sumCsIE_cal_um", "", maxadcchan, 0, maxadcchan, 512, 0, 200);
	Histo.sumCsITime_um      = new TH2I("sumCsITime_um", "", maxadcchan, 0, maxadcchan, 1000, -500, 500);
	Histo.sumCsITime_matched = new TH2I("sumCsITime_matched", "", nTelCsIs * 4, 0, nTelCsIs * 4, 1000, -500, 500);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gobbi28::~Gobbi28() {
	delete FrontEcal;
	delete BackEcal;
	delete CsIEcal;
	delete FrontTimecal;
	delete BackTimecal;
	for (int i = 0; i < 4; i++) delete Telescope[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi28::analyze() {

	// Reset member variables
	for (size_t i = 0; i < 4; i++) Telescope[i]->reset();

	// Loop through HINP hits, store in Telescope array and output to Histo
	size_t nhits = input.GetNhits();
	for (size_t i = 0; i < nhits; i++) {
		size_t board = input.GetBoard(i);
		size_t chan = input.GetChan(i);
		if (board > hinpboards || chan >= hinpchans) {
			cout << "Nhits " << nhits << endl;
			cout << "i " << i << endl;
			cout << "Board " << board << " and chan " << chan;
			cout << "Hit skipped!" << endl;
			return;
		}
		
		// Use calibration to get Energy and fill elist class in telescope
		if (board == 1 || board == 3 || board == 5 || board == 7)
			addFrontHit((board - 1) / 2, chan, i);
		else if (board == 2 || board == 4 || board == 6 || board == 8)
			addBackHit((board / 2) - 1, chan, i);
		else {
			cout << "Invalid board number " << board << endl;
			cout << "Hit skipped!" << endl;
			return;
		}
	}
	
	// Loop through CsI ADC/QDC/TDC hits, match and store in Telescope array
	addCsIHits();

	// Loop through telescopes and perform additional analysis
	size_t sumchan = 0;
	size_t Pidmulti = 0;
	int NCsI_all = 0; //number of CsI hits across all telescopes
	int NCsI_matched = 0; //number of CsI hits with good F and B
	for (size_t id = 0; id < 4; id++) {

#ifdef ENABLE_DEBUG
		cout << "id " << id << ", Nsol " << Telescope[id]->Nsolution << endl;
#endif

		// This is the spot where we run Telescope->Neighbours() for addback
		Telescope[id]->Front.Neighbours(id);
		Telescope[id]->Back.Neighbours(id);
		
		int FrontN = Telescope[id]->Front.Nstore;
		int BackN  = Telescope[id]->Back.Nstore;
		int CsIN   = Telescope[id]->CsI.Nstore;
		
		NCsI_all += CsIN;

		// Then fill summary histograms after addback
		for (size_t n = 0; n < FrontN; n++) {
			sumchan = id*hinpchans + Telescope[id]->Front.Order[n].strip;
			Histo.sumFrontE_addback->Fill(sumchan, Telescope[id]->Front.Order[n].energy);
		}
		for (size_t n = 0; n < BackN; n++) {
			sumchan = id*hinpchans + Telescope[id]->Back.Order[n].strip;
			Histo.sumBackE_addback->Fill(sumchan, Telescope[id]->Back.Order[n].energy);
		}
		
		// No point if not at least one front and one back
		if (FrontN < 1 || BackN < 1) continue;

		// Simple case of single strip multiplicity in front/back for CsI-less position map
		// Useful for alpha calibrations or other data without CsI information
		if (FrontN == 1 && BackN == 1) {
			Telescope[id]->testingHitE();
			Histo.xyhitmap_sionly->Fill(Telescope[id]->tempSol.Xpos, Telescope[id]->tempSol.Ypos);
			Histo.sumFrontTimeMult1_cal->Fill(id*hinpchans + Telescope[id]->Front.Order[0].strip, Telescope[id]->Front.Order[0].time);
			Histo.FrontvsBack_sionly[id]->Fill(Telescope[id]->tempSol.energy, Telescope[id]->tempSol.benergy);
		}
		
		// For everything else, no point if there is not at least one CsI hit
		if (CsIN < 1) continue;
		
		// Handle simple case of one hit each in front, back, and CsI
		size_t numMultCsI = 0;
		if (FrontN == 1 && BackN == 1 && CsIN == 1) NsimpleECsI = Telescope[id]->simpleECsI();
		else {
			numMultCsI = Telescope[id]->multiHitECsI();
			NmultiECsI += numMultCsI;	
		}
		
		NCsI_matched += Telescope[id]->Nsolution;

		// Next, fill E vs. dE and other plots
		for (size_t isol = 0; isol < Telescope[id]->Nsolution; isol++) {
			Telescope[id]->position(isol); // calculates (x, y) pos and lab angle
			solution& sol = Telescope[id]->Solution[isol];

#ifdef ENABLE_DEBUG
			cout << "isol " << isol << endl; 
			cout << "cos " << cos(sol.theta) << endl;
			cout << "the " << sol.theta * rad_to_deg << endl;
			cout << "front strip " << sol.ifront << ", back strip " << sol.iback << endl;
			cout << "x " << sol.Xpos << ", y " << sol.Ypos << endl;
			cout << "E " << sol.energy << ", dE " << sol.denergy << endl;
#endif

			double th = sol.theta * rad_to_deg;
			size_t icsi = sol.iCsI;
			double denergy_corr = sol.denergy * cos(sol.theta);
			Histo.DEE_CsI[id][icsi]->Fill(sol.energyR, denergy_corr);
			if (numMultCsI > 1) Histo.xyhitmap_multiCsI->Fill(sol.Xpos, sol.Ypos);
			if (sol.time >= 6000. && sol.time <= 11000.) Histo.DEE_CsI_sitgate[id][icsi]->Fill(sol.energyR, denergy_corr);
			if (sol.CsITime >= -100. && sol.CsITime <= -80.) Histo.DEE_CsI_csitgate[id][icsi]->Fill(sol.energyR, denergy_corr);
			if (sol.ifront % 2 == 0) Histo.DEE_CsI_fronteven[id][icsi]->Fill(sol.energyR, denergy_corr);
			else Histo.DEE_CsI_frontodd[id][icsi]->Fill(sol.energyR, denergy_corr);
			Histo.DEE_CsI_BackE[id][icsi]->Fill(sol.energyR, sol.benergy * cos(sol.theta));
			Histo.FrontvsBack[id]->Fill(sol.benergy, sol.denergy);
			Histo.xyhitmap->Fill(sol.Xpos, sol.Ypos);
			Histo.tphitmap->Fill(th * cos(sol.phi), th * sin(sol.phi));
			Histo.Evstheta[id]->Fill(th, sol.energy);
			Histo.Evstheta_all->Fill(th, sol.energy);
			Histo.Theta->Fill(th);
			Histo.CsI_Time_matched[id][icsi]->Fill(sol.CsITime);
			Histo.CsI_xyhitmap[id][icsi]->Fill(sol.Xpos, sol.Ypos);
			Histo.sumCsITime_matched->Fill(icsi + (nTelCsIs * id), sol.CsITime);
			
			// At this point in time, the CsI crystals should be calibrated to proton equivalent energy
			// The call to `telescope::getPID` applies PID-dependent stage 2 calibrations to all relevant non-proton particles
			// This transforms proton equivalent energy to normal energy
			Histo.CsI_Energy_R[id][icsi]->Fill(sol.energyR);
			Histo.CsI_Energy_pcal[id][icsi]->Fill(sol.energy);
			if (Telescope[id]->isCenter(sol.ifront, sol.iback)) {
				Histo.CsI_Energy_R_center[id][icsi]->Fill(sol.energyR);
				Histo.CsI_Energy_pcal_center[id][icsi]->Fill(sol.energy);
			}

			/******** ANGLE DEPENDENT CALIBRATION CORRECTION FOR HIGHER ENERGY POINTS ********/
			// This was added my Nicholas Dronchi, and can generally be ignored or bypassed
			// unless one finds that they need an angle dependent calibration. I (Henry Webb)
			// have left this in for now since all it does is fill some extra histograms.

			// NOTE: SpecTcl should already map channels to strips, and this code has been modified to skip the unpacking and directly read in a SpecTcl output tree
			int chandE = sol.ifront;

			// Make a correction to the dE silicon energy based on angle
			double angle_dEcorr = 1.0277e-5*(th*th*th) + 1.6125e-3*(th*th) + 8.3097e-4*th - 1.0227e-3;
			double dEcorr = sol.denergy + angle_dEcorr;
			double dEcorr_R = FrontEcal->reverseCal(id, chandE, dEcorr);
			Histo.AngleCorrDeltaE[id][chandE]->Fill(dEcorr);
			Histo.AngleCorrDeltaE_noCorr[id][chandE]->Fill(sol.denergy);
			Histo.AngleCorrDeltaE_R[id][chandE]->Fill(dEcorr_R);
			Histo.AngleCorrDeltaE_cal->Fill(id*hinpchans + chandE, dEcorr);
			double Etot = sol.energy + sol.denergy;
			Histo.sumEtot_cal->Fill(id*hinpchans + sol.ifront, Etot);
			Histo.AngleCorrSum_cal->Fill(id*hinpchans + sol.ifront, sol.energy + dEcorr);

#ifdef ENABLE_DEBUG
			cout << "th " << th << ", angle_dEcorr " << angle_dEcorr << " MeV" << endl;
			if (dEcorr > 5) cout << "dEnergyR " << sol.denergyR << ", dEcorr " << dEcorr << ", dEcorr_R " << dEcorr_R << endl;
#endif

		}

		// Calculate and determine particle identification (PID) in the telescope for each solution
		Pidmulti += Telescope[id]->getPID();

		// Hit maps and other plots based on Pid
		for (size_t isol = 0; isol < Telescope[id]->Nsolution; isol++) {
			solution& sol = Telescope[id]->Solution[isol];
			if (sol.ipid == 0) continue; // enforce that all histograms filled here contain hits with valid PIDs

			double xpos = sol.Xpos;
			double ypos = sol.Ypos;

			// Gate on specific (Z, A) PID results
			double tdiff = sol.timediff;
			pair<size_t, size_t> ZA(sol.iZ, sol.iA);
			if (ZA == sz_pair(1, 1)) { // protons
				Histo.protonhitmap->Fill(xpos, ypos);
				Histo.dTime_proton->Fill(tdiff);
			}
			else if (ZA == sz_pair(1, 2)) { // deuterons
				Histo.deuteronhitmap->Fill(xpos, ypos);
				Histo.dTime_deuteron->Fill(tdiff);
			}
			else if (ZA == sz_pair(1, 3)) { // tritons
				Histo.tritonhitmap->Fill(xpos, ypos);
				Histo.dTime_triton->Fill(tdiff);
			}
			else if (ZA == sz_pair(2, 4)) { // alphas
				Histo.alphahitmap->Fill(xpos, ypos);
				Histo.dTime_alpha->Fill(tdiff);
			}
			else if (ZA == sz_pair(2, 6)) { // 6He
				Histo.He6hitmap->Fill(xpos, ypos);
				Histo.dTime_He6->Fill(tdiff);
			}

			//// Other PID results

			// Any lithium
			if (sol.iZ == 3) {
				Histo.Lihitmap->Fill(xpos, ypos);
				Histo.dTime_Li->Fill(tdiff);
			}
			// Lithium veto
			if (sol.iZ != 3 && sol.iA != 7) {
				Histo.LiVETOhitmap->Fill(xpos, ypos);
			}
		}

		// Calculate sumEnergy, then account for Eloss in target, then set Ekin and momentum of solutions
		// Eloss files are loaded in telescope
		pidSkipped += Telescope[id]->calcEloss();

		// Post Eloss calculation plots and output
		for (size_t isol = 0; isol < Telescope[id]->Nsolution; isol++) {
			solution& sol = Telescope[id]->Solution[isol];
			Histo.solutions.push_back(Telescope[id]->Solution[isol]);
			if (sol.iZ == 1 && sol.iA == 1)
				Histo.ProtonEnergy->Fill(sol.Ekin, sol.theta * rad_to_deg);
			else if (sol.iZ == 2 && sol.iA == 3)
				Histo.He3Energy->Fill(sol.Ekin, sol.theta * rad_to_deg);
			else if (sol.iZ == 2 && sol.iA == 4)
				Histo.AlphaEnergy->Fill(sol.Ekin, sol.theta * rad_to_deg);
		}
	}
	
	//Multiplicity of CsI crystals with good ADC,QDC,TDC
	Histo.CsI_mult_AQT->Fill(NCsI_all);
	// " with good F and B strips
	Histo.CsI_mult_M->Fill(NCsI_matched);
	Histo.CsI_mult_M_v_R->Fill(NCsI_all, NCsI_matched);
	// " with good PID
	Histo.CsI_mult_PID->Fill(Pidmulti);
	Histo.CsI_mult_PID_v_R->Fill(NCsI_all, Pidmulti);
/*
	// Debug output for tracking events:
	if (Ran->Rndm() > 0.1) return;
	if (NCsI_all == 2 && (Pidmulti == 0 || Pidmulti == 1)) {
		stringstream ss;
		ss << GREEN << "=======================================================" << endl;
		bool foundNeqFB = false;
		for (size_t i = 0; i < 4; i++) {
			if (Telescope[i]->Front.Nstore != Telescope[i]->Back.Nstore) foundNeqFB = true;
			ss << "telescope: " << i << endl;
			ss << "Front Si hits: " << endl;
			for (size_t j = 0; j < Telescope[i]->Front.Nstore; j++) {
				order* o = &Telescope[i]->Front.Order[j];
				ss << "strip: " << o->strip << ", energy: " << o->energy << ", time: " << o->time << endl;
			}
			ss << "Back Si hits: " << endl;
			for (size_t j = 0; j < Telescope[i]->Back.Nstore; j++) {
				order* o = &Telescope[i]->Back.Order[j];
				ss << "strip: " << o->strip << ", energy: " << o->energy << ", time: " << o->time << endl;
			}
			ss << "CsI hits: " << endl;
			for (size_t j = 0; j < Telescope[i]->CsI.Nstore; j++) {
				order* o = &Telescope[i]->CsI.Order[j];
				ss << "CsI ID: " << o->strip << ", energy: " << o->energy << ", TDC: " << o->time << ", QDC: " << o->qdc << endl;
			}
		}
		if (!foundNeqFB) return;
		ss << "=======================================================" << endl;
		for (size_t i = 0; i < 4; i++) {
			ss << "telescope: " << i << endl;
			for (size_t j = 0; j < Telescope[i]->Nsolution; j++) {
				solution& sol = Telescope[i]->Solution[j];
				ss << "ifront: " << sol.ifront << ", iback: " << sol.iback << ", iCsI: " << ", denergy: " << sol.denergy << ", benergy: " << sol.benergy << ", energy: " << sol.energy << ", energyR: " << sol.energyR << ", Z: " << sol.iZ << ", A: " << sol.iA << endl;
			}
		}
		ss << "=======================================================" << RESET << endl;
		cout << ss.str();
		abort();
	}
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t Gobbi28::loadSolutions(correl2& Correl) {
	size_t goodMult = 0;
	for (size_t id = 0; id < 4; id++) {
		for (size_t isol = 0; isol < Telescope[id]->Nsolution; isol++) {
			solution& sol = Telescope[id]->Solution[isol];

			// Only keep solutions that have a valid pid
			if (sol.ipid > 0) {
				Correl.load(&sol);
				goodMult++;
			}
		}
	}

	return goodMult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

solution* Gobbi28::getNextEmptySolution(solution* sol) {
	bool foundEmptySolution = false;
	for (size_t id = 0; id < 4; id++) {
		telescope* tel = Telescope[id];
		for (size_t isol = 0; isol < tel->Nsolution; isol++) {
			if (tel->Solution[isol].ipid != sol->ipid) continue;
			else if (tel->Nsolution >= 20)
				throw invalid_argument("ERROR: max solutions reached, cannot get next empty solution");
			return &(tel->Solution[tel->Nsolution]);
		}
	}
	return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Pass through function for use with extra calculations external to this
// Gobbi class. All telescope losses should be the same, so just use the first
// one.
float Gobbi28::getEin(float energy, float thick, int Z, float A) {
	return Telescope[0]->losses->getEin(energy, thick, Z, A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi28::addFrontHit(size_t tel, size_t ch, size_t i) {
	size_t ER = input.GetE(i);
	size_t tR = input.GetT(i);
	size_t ELoR = input.GetELo(i);
	double Energy = FrontEcal->getEnergy(tel, ch, ER);
	double time = FrontTimecal->getTime(tel, ch, tR);

	Histo.sumFrontE_R->Fill(tel*hinpchans + ch, ER);
	Histo.sumFrontTime_R->Fill(tel*hinpchans + ch, tR);
	Histo.sumFrontE_cal->Fill(tel*hinpchans + ch, Energy);
	Histo.sumFrontTime_cal->Fill(tel*hinpchans + ch, time);

	Histo.FrontE_R[tel][ch]->Fill(ER);
	Histo.FrontElow_R[tel][ch]->Fill(ELoR);
	Histo.FrontTime_R[tel][ch]->Fill(tR);
	Histo.FrontE_cal[tel][ch]->Fill(Energy);

	//if (Energy > .5 && tR > 3420 && tR < 6380 && (tel != 1 || Energy > 1.8))
	//need to set thresholds just above noise
	//if (tel == 1 && (

	if (Energy > .5) { //(tel != 1 || Energy > 2)
		Telescope[tel]->Front.Add(ch, Energy, ELoR, ER, time);
		Telescope[tel]->multFront++;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi28::addBackHit(size_t tel, size_t ch, size_t i) {
	size_t ER = input.GetE(i);
	size_t tR = input.GetT(i);
	size_t ELoR = input.GetELo(i);
	double Energy = BackEcal->getEnergy(tel, ch, ER);
	double time = BackTimecal->getTime(tel, ch, tR);

	Histo.sumBackE_R->Fill(tel*hinpchans + ch, ER);
	Histo.sumBackTime_R->Fill(tel*hinpchans + ch, tR);
	Histo.sumBackE_cal->Fill(tel*hinpchans + ch, Energy);
	Histo.sumBackTime_cal->Fill(tel*hinpchans + ch, time);

	Histo.BackE_R[tel][ch]->Fill(ER);
	Histo.BackElow_R[tel][ch]->Fill(ELoR);
	Histo.BackTime_R[tel][ch]->Fill(tR);
	Histo.BackE_cal[tel][ch]->Fill(Energy);

	if (Energy > .5) {
		Telescope[tel]->Back.Add(ch, Energy, ELoR, ER, time);
		Telescope[tel]->multBack++;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi28::addCsIHits() {
	int Nfound = 0;
	int Nnotfound = 0;
	
	// Fill unmatched QDC histograms
	size_t nQhits = input_qdc.GetNhits();
	for (size_t iq = 0; iq < nQhits; iq++) {
	
		// Test for valid QDC channel
		size_t qdcchan = input_qdc.GetChan(iq);
		if (Histo.CsI_QDC_um.find(qdcchan) == Histo.CsI_QDC_um.end()) {
			//cout << string(BOLDRED) << "WARNING: QDC channel " << to_string(qdcchan) << " not found in CsI map" << string(RESET) << endl;
			
#ifdef ENABLE_DEBUG
			for (const auto& [key, value] : Histo.CsI_QDC_um)
				cout << key << endl;
#endif
			
			continue;
		}
		
		// Assume that qdc and adc channels are the same for each CsI
		Histo.CsI_QDC_um[qdcchan]->Fill(input_qdc.GetAQ(iq));
	}
	
	// Fill unmatched TDC histograms
	// This works differently from the other input classes, see `Input.h` for details
	for (size_t it = tdcstart; it <= tdcstart + maxadcchan; it++) {
		
		// Skip if no output histogram exists for this tdc channel
		if (Histo.CsI_Time_um.find(it) == Histo.CsI_Time_um.end()) continue;
		
		// Test to see if this tdc channel contains at least one hit
		optional<double> T = input_tdc.GetT(it, 0);
		if (T == nullopt) continue;
		
		Histo.CsI_Time_um[it]->Fill(T.value());
		Histo.sumCsITime_um->Fill(it - tdcstart, T.value());
	}
	
	// Fill unmatched ADC histograms
	// Also perform matching and fill matched histograms
	size_t nEhits = input_adc.GetNhits();
	
	//ADC multiplicity
	Histo.CsI_mult_R->Fill(nEhits);
	
	for (size_t ie = 0; ie < nEhits; ie++) {
		size_t adcchan = input_adc.GetChan(ie);
		size_t ER = input_adc.GetAQ(ie);
		size_t tel = telCsImap[adcchan];
		size_t id = idCsImap[adcchan];
		double Ecal = CsIEcal->getEnergy(tel, id, ER);
		
		// Test for valid ADC channel
		if (Histo.CsI_Energy_R_um.find(adcchan) == Histo.CsI_Energy_R_um.end()) {
			//cout << string(BOLDRED) << "WARNING: ADC channel " << to_string(adcchan) << " not found in CsI map" << string(RESET) << endl;
			
#ifdef ENABLE_DEBUG
			for (const auto& [key, value] : Histo.CsI_Energy_R_um)
				cout << key << endl;
#endif
			
			continue;
		}
		
		Histo.CsI_Energy_R_um[adcchan]->Fill(ER);
		Histo.CsI_Energy_cal_um[adcchan]->Fill(Ecal);
		Histo.sumCsIE_R_um->Fill(adcchan, ER);
		Histo.sumCsIE_cal_um->Fill(adcchan, Ecal);
		
		bool hasTDC = false;
		bool hasQDC = false;
		optional<double> T;
		size_t tdcchan, qdcchan, Q;
		for (size_t it = tdcstart; it <= tdcstart + maxadcchan; it++) {
			
			// Skip if no output histogram exists for this tdc channel
			if (Histo.CsI_Time_um.find(it) == Histo.CsI_Time_um.end()) continue;
			
			// Test to see if this tdc channel contains at least one hit
			T = input_tdc.GetT(it, 0);
			if (T == nullopt) continue;

			tdcchan = it - 16;
				
#ifdef ENABLE_DEBUG
			cout << "adcchan " << adcchan << " it - 16 " << tdcchan << endl;
#endif
			
			if (tdcchan != adcchan) continue;
			
			hasTDC = true;
			break;
		}
		
		for (size_t iq = 0; iq < nQhits; iq++) {
			qdcchan = input_qdc.GetChan(iq);
			Q = input_qdc.GetAQ(iq);
			if (qdcchan == qdcchan) Histo.CsI_QDC_matched[adcchan]->Fill(Q);
		
#ifdef ENABLE_DEBUG
			cout << "Gobbi28::addCsIHits adcchan " << adcchan << " qdcchan " << qdcchan << endl;
#endif
			
			if (qdcchan != adcchan) continue;
			
			hasQDC = true;
			Histo.CsIonly_PSD[adcchan]->Fill(ER, Q);
			break;
		}
		
		double t = NAN;
		if (hasTDC) t = CsITimecal->getTime(tel, id, T.value());
		if (!hasQDC) Q = -1;
		Telescope[tel]->CsI.Add(id, Ecal, 0., 0, ER, t, Q, true);
		Telescope[tel]->multCsI++;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



