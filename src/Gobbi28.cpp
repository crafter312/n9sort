/* New Gobbi class for the Gobbi 28 configuration by Henry Webb (h.s.webb@wustl.edu).
 * Created 7 May 2026 as mix of `OldGobbi` class from this code and `gobbi` class
 * used by Johnathan Phillips (j.s.phillips@wustl.edu) for 22Si FRIB experiment.
 */

#include "Gobbi28.h"

#include <stdexcept>
#include <string>
#include <utility>

#include <TH1I.h>

#include "constants.h"
#include "solution.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gobbi28::Gobbi28(Input& in, histo& hist, SortConfig& config) : Targetdist(config.GetTargDist()), TargetThickness(config.GetTargThick()), hinpboards(config.GetHinpboards()), hinpchans(config.GetHinpchans()), Histo(hist), input(in.GetGobbi()), input_qdc(in.GetQDC()), input_adc(in.GetADC()), input_tdc(in.GetTDC()) {
	for (size_t id = 0; id < 4; id++) {
		Telescope[id] = new telescope(TargetThickness, config);
		Telescope[id]->init(id, config); // tells Telescope what position it is in
		Telescope[id]->SetTargetDistance(Targetdist);
	}
	NCsI = Telescope[0]->GetNCsI(); // same for every telescope

	string calDir = config.GetCalDir();
	FrontEcal = new calibrate(4, hinpchans, calDir + config.GetFrontEcalFile(), 1, false);
	BackEcal = new calibrate(4, hinpchans, calDir + config.GetBackEcalFile(), 1, false);
	CsIEcal = new calibrate(4, NCsI, calDir + config.GetCsIEcalFile(), 1, false);
	FrontTimecal = new calibrate(4, hinpchans, calDir + config.GetFrontTimecalFile(), 1, false);
	BackTimecal = new calibrate(4, hinpchans, calDir + config.GetBackTimecalFile(), 1, false);
	CsITimecal = new calibrate(4, NCsI, calDir + config.GetCsITimecalFile(), 1, false);
	
	// Read in CsI mappings from file
	string incsimapfile = config.GetConfigDir() + config.GetCsIChannelMapFile();
	ifstream incsimap(incsimapfile);
	if (incsimap.fail()) throw invalid_argument(string(BOLDRED) + string("CsI channel map file ") + incsimapfile + string(" does not exist or failed to open") + string(RESET));
	else cout << GREEN << "CsI channel map file " << inextentsfile << " opened" << RESET << endl;
	
	size_t chan, tel, id;
	string name;
	while (incsimap.good()) {
		incsimap >> chan >> tel >> id;
		if ((tel != 0) && (tel != 1) && (tel != 2) && (tel != 3)) throw invalid_argument(string(BOLDRED) + string("CsI telescope # must be 0, 1, 2, or 3") + string(RESET));
		telCsImap[chan] = tel;
		idCsImap[chan] = id;
		
		if (chan > maxadcchan) maxadcchan = chan;
		
		//// Make per-CsI crystal histograms
		
		Histo->dir1dCsI_Energy->cd();
		name = "CsI_Energy_" + string(chan) + "_R_unmatched";
		Histo->CsI_Energy_R_um[chan] = new TH1I(name.c_str(), "", 1024, 0, 4096);
		name = "CsI_Energy_" + string(chan) + "_cal_unmatched";
		Histo->CsI_Energy_cal_um[chan] = new TH1I(name.c_str(), "", 512, 0, 200);
		
		Histo->dir1dCsI_QDC->cd();
		name = "CsI_QDC_" + string(chan) + "_um";
		Histo->CsI_QDC_um[chan] = new TH1I(name.c_str(), "", 1024, 0, 4096);
		name = "CsI_QDC_" + string(chan) + "_matched";
		Histo->CsI_QDC_matched[chan] = new TH1I(name.c_str(), "", 1024, 0, 4096);
		
		Histo->dir1dCsI_Time->cd();
		name = "CsI_Time_" + string(chan) + "_um";
		Histo->CsI_Time_um[tdcstart + chan] = new TH1I(name.c_str(), "", 1000, -500, 500);
	}
	
	// Make CsI summary plots
	Histo->dirSummary->cd();
	Histo->sumCsIE_R_um   = new TH2I("sumCsIE_R_um", "", maxadcchan, 0, maxadcchan, 1024, 0, 4096);
	Histo->sumCsIE_cal_um = new TH2I("sumCsIE_cal_um", "", maxadcchan, 0, maxadcchan, 512, 0, 200);
	Histo->sumCsITime_um  = new TH2I("sumCsITime_um", "", maxadcchan, 0, maxadcchan, 1000, -500, 500);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gobbi28::~Gobbi28() {
	delete FrontEcal;
	delete BackEcal;
	delete CsIEcal;
	delete FrontTimecal;
	delete BackTimecal;
	delete CsITimecal;
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
	

	// Loop through telescopes and perform additional analysis
	size_t sumchan = 0;
	size_t totMulti = 0;
	size_t Pidmulti = 0;
	for (size_t id = 0; id < 4; id++) {

		// This is the spot where we run Telescope->Neighbours() for addback
		Telescope[id]->Front.Neighbours(id);
		Telescope[id]->Back.Neighbours(id);

		// Then fill summary histograms after addback
		for (size_t n = 0; n < Telescope[id]->Front.Nstore; n++) {
			sumchan = id*hinpchans + Telescope[id]->Front.Order[n].strip;
			Histo.sumFrontE_addback->Fill(sumchan, Telescope[id]->Front.Order[n].energy);
		}
		for (size_t n = 0; n < Telescope[id]->Back.Nstore; n++) {
			sumchan = id*hinpchans + Telescope[id]->Back.Order[n].strip;
			Histo.sumBackE_addback->Fill(sumchan, Telescope[id]->Back.Order[n].energy);
		}

		// Handle simple case of single strip multiplicity
		bool isSimple = false;
		if (Telescope[id]->Front.Nstore == 1 && Telescope[id]->Delta.Nstore == 1) {
			Histo.frontdeltastripnum[id]->Fill(Telescope[id]->Front.Order[0].strip, Telescope[id]->Delta.Order[0].strip);
			Histo.timediff[id]->Fill(Telescope[id]->Front.Order[0].time - Telescope[id]->Delta.Order[0].time);

			if (Telescope[id]->Back.Nstore == 1) {
				totMulti += Telescope[id]->simpleFront();
				Histo.sumFrontTimeMult1_cal->Fill(id*hinpchans + Telescope[id]->Front.Order[0].strip, Telescope[id]->Front.Order[0].time);
				isSimple = true;
			}
		}

		// If higher multiplicity then worry about picking the right one
		// This also handles the case where Nstore = 0 for any of the chanels
		if (!isSimple) totMulti += Telescope[id]->multiHit();

		// Next, fill E vs. dE and other plots

#ifdef ENABLE_DEBUG
		cout << "id " << id << ", Nsol " << Telescope[id]->Nsolution << endl;
#endif

		for (size_t isol = 0; isol < Telescope[id]->Nsolution; isol++) {
			Telescope[id]->position(isol); // calculates x,y pos, and lab angle

#ifdef ENABLE_DEBUG
			cout << "isol " << isol << endl; 
			cout << "cos " << cos(Telescope[id]->Solution[isol].theta) << endl;
			cout << "the " << Telescope[id]->Solution[isol].theta * rad_to_deg << endl;
			cout << "front strip " << Telescope[id]->Solution[isol].ifront << ", back strip " << Telescope[id]->Solution[isol].iback << endl;
			cout << "x " << Telescope[id]->Solution[isol].Xpos << ", y " << Telescope[id]->Solution[isol].Ypos << endl;
			cout << "E " << Telescope[id]->Solution[isol].energy << ", dE " << Telescope[id]->Solution[isol].denergy << endl;
#endif

			Histo.FrontvsBack[id]->Fill(Telescope[id]->Solution[isol].energy, Telescope[id]->Solution[isol].benergy);
			double Ener = Telescope[id]->Solution[isol].energy + Telescope[id]->Solution[isol].denergy * (1 - cos(Telescope[id]->Solution[isol].theta));
			Histo.DEE[id]->Fill(Ener, Telescope[id]->Solution[isol].denergy * cos(Telescope[id]->Solution[isol].theta));
			Histo.xyhitmap->Fill(Telescope[id]->Solution[isol].Xpos, Telescope[id]->Solution[isol].Ypos);
			double th = Telescope[id]->Solution[isol].theta * rad_to_deg;
			Histo.Evstheta[id]->Fill(th, Telescope[id]->Solution[isol].energy);
			Histo.Evstheta_all->Fill(th, Telescope[id]->Solution[isol].energy);
			Histo.Theta->Fill(th);

			/******** ANGLE DEPENDENT CALIBRATION CORRECTION FOR HIGHER ENERGY POINTS ********/
			// This was added my Nicholas Dronchi, and can generally be ignored or bypassed
			// unless one finds that they need an angle dependent calibration. I (Henry Webb)
			// have left this in for now since all it does is fill some extra histograms.

			// NOTE: SpecTcl should already map channels to strips, and this code has been modified to skip the unpacking and directly read in a SpecTcl output tree
			int chan = Telescope[id]->Solution[isol].ifront;

			// Make a correction to the E silicon energy based on angle
			double angle_Ecorr = 1.0277e-5*(th*th*th) + 1.6125e-3*(th*th) + 8.3097e-4*th - 1.0227e-3;
			double Ecorr = Telescope[id]->Solution[isol].energy + angle_Ecorr;
			double Ecorr_R = FrontEcal->reverseCal(id, Telescope[id]->Solution[isol].ifront, Ecorr);

			Histo.AngleCorrE[id][chan]->Fill(Ecorr);
			Histo.AngleCorr_noCorr[id][chan]->Fill(Telescope[id]->Solution[isol].energy);
			Histo.AngleCorrE_R[id][chan]->Fill(Ecorr_R);
			Histo.AngleCorrFrontE_cal->Fill(id*hinpchans + Telescope[id]->Solution[isol].ifront, Ecorr);

			// NOTE: SpecTcl should already map channels to strips, and this code has been modified to skip the unpacking and directly read in a SpecTcl output tree
			int chandE = Telescope[id]->Solution[isol].ide;

			// Make a correction to the dE silicon energy based on angle
			double angle_dEcorr = -1.0971e-5*(th*th*th) - 1.1446e-3*(th*th) - 8.9371e-4*th + 1.0879e-3;
			//double angle_dEcorr = 2.0527e-6*(th*th*th) - 1.4281e-3*(th*th) + 1.5589e-4*th - 7.3389e-4; // no Au foilloss
			double dEcorr = Telescope[id]->Solution[isol].denergy + angle_dEcorr;
			double dEcorr_R = DeltaEcal->reverseCal(id,Telescope[id]->Solution[isol].ide, dEcorr);

			Histo.AngleCorrDeltaE[id][chandE]->Fill(dEcorr);
			Histo.AngleCorrDeltaE_noCorr[id][chandE]->Fill(Telescope[id]->Solution[isol].denergy);
			Histo.AngleCorrDeltaE_R[id][chandE]->Fill(dEcorr_R);
			Histo.AngleCorrDeltaE_cal->Fill(id*hinpchans + Telescope[id]->Solution[isol].ide, dEcorr);

			double Etot = Telescope[id]->Solution[isol].energy + Telescope[id]->Solution[isol].denergy;
			Histo.sumEtot_cal->Fill(id*hinpchans + Telescope[id]->Solution[isol].ifront, Etot);
			Histo.AngleCorrSum_cal->Fill(id * hinpchans + Telescope[id]->Solution[isol].ifront, Ecorr + dEcorr);

#ifdef ENABLE_DEBUG
			cout << "th " << th << ", angle_Ecorr " << angle_Ecorr << " MeV, angle_dEcorr " << angle_dEcorr << " MeV" << endl;
			if (Ecorr > 5) cout << "EnergyR " << Telescope[id]->Solution[isol].energyR << ", Ecorr " << Ecorr << ", Ecorr_R " << Ecorr_R << endl;
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

			// Gated on or A time
			if (input_tdc.t[1][0] >= -80 && input_tdc.t[1][0] <= -50)
				Histo.xyhitmap_tgate_orA->Fill(xpos, ypos);

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
		Telescope[id]->calcEloss();

		// Post Eloss calculation plots, I guess?
		for (size_t isol = 0; isol < Telescope[id]->Nsolution; isol++) {
			solution& sol = Telescope[id]->Solution[isol];
			if (sol.iZ == 1 && sol.iA == 1) {
				Histo.ProtonEnergy->Fill(sol.Ekin, sol.theta * rad_to_deg);
			}
		}
	}
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
		
		// Assume that qdc and adc channels are the same for each CsI
		Histo->CsI_QDC_um[input_qdc.GetChan(iq)]->Fill(input_qdc.GetAQ(iq));
	}
	
	// Fill unmatched TDC histograms
	// This works differently from the other input classes, see `Input.h` for details
	for (size_t it = tdcstart; it <= tdcstart + maxadcchan; it++) {
		
		// Skip if no output histogram exists for this tdc channel
		if (Histo->CsI_Time_um.find(it) == Histo->CsI_Time_um.end()) continue;
		
		// Test to see if this tdc channel contains at least one hit
		optional<double> T = input_tdc.GetT(it, 0);
		if (T == nullopt) continue;
		
		Histo->CsI_Time_um[it]->Fill(T.value());
		Histo->sumCsITime_um->Fill(it - tdcstart, T.value());
	}
	
	// Fill unmatched ADC histograms
	// Also perform matching and fill matched histograms
	size_t nEhits = input_adc.GetNhits();
	for (size_t ie = 0; ie < nEhits; ie++) {
		size_t adcchan = input_adc.GetChan(ie);
		size_t ER = input_adc.GetAQ(ie);
		size_t tel = telCsImap[adcchan];
		size_t id = idCsImap[adcchan];
		double Ecal = CsIEcal->getEnergy(tel, id, ER);
		
		Histo->CsI_Energy_R_um[adcchan]->Fill(ER);
		Histo->CsI_Energy_cal_um[adcchan]->Fill(Ecal);
		Histo->sumCsIE_R_um->Fill(adcchan, ER);
		Histo->sumCsIE_cal_um->Fill(adcchan, Ecal);
		
		for (size_t iq = 0; iq < nQhits; iq++) {
			size_t qdcchan = input_qdc.GetChan(iq);
			size_t Q = input_qdc.GetAQ(iq);
			
			if (qdcchan == qdcchan) Histo->CsI_QDC_matched[adcchan]->Fill(Q);
			
			for (size_t it = tdcstart; it <= tdcstart + maxadcchan; it++) {
				
				// Skip if no output histogram exists for this tdc channel
				if (Histo->CsI_Time_um.find(it) == Histo->CsI_Time_um.end()) continue;
				
				// Test to see if this tdc channel contains at least one hit
				optional<double> T = input_tdc.GetT(it, 0);
				if (T == nullopt) continue;
				
				if ((qdcchan != adcchan) || (it != adcchan)) continue;
				
				Telescope[tel]->CsI.Add(id, Ecal, 0., 0, ER, T.value(), Q, true);
				Telescope[tel]->multCsI++;
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



