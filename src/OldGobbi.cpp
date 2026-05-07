/* Split off from `Det` (previously `Gobbi`) by Henry Webb (h.s.webb@wustl.edu)
 * 4 May 2026 to match other, newer versions of the Wash U radiochemistry sort
 * code than the code I started with.
 */

#include "OldGobbi.h"

#include <stdexcept>
#include <string>
#include <utility>

#include "constants.h"
#include "solution.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OldGobbi::OldGobbi(Input& in, histo& hist, SortConfig& config) : Targetdist(config.GetTargDist()), TargetThickness(config.GetTargThick()), hinpboards(config.GetHinpboards()), hinpchans(config.GetHinpchans()), Histo(hist), input(in.GetGobbi()), input_tdc(in.GetTDC()) {
	for (size_t id = 0; id < 4; id++) {
		Silicon[id] = new silicon(TargetThickness, config);
		Silicon[id]->init(id, config); // tells Silicon what position it is in
		Silicon[id]->SetTargetDistance(Targetdist);
	}

	string calDir = config.GetCalDir();
	FrontEcal = new calibrate(4, hinpchans, calDir + config.GetFrontEcalFile(), 1, false);
	BackEcal = new calibrate(4, hinpchans, calDir + config.GetBackEcalFile(), 1, false);
	DeltaEcal = new calibrate(4, hinpchans, calDir + config.GetDeltaEcalFile(), 1, false);
	FrontTimecal = new calibrate(4, hinpchans, calDir + config.GetFrontTimecalFile(), 1, false);
	BackTimecal = new calibrate(4, hinpchans, calDir + config.GetBackTimecalFile(), 1, false);
	DeltaTimecal = new calibrate(4, hinpchans, calDir + config.GetDeltaTimecalFile(), 1, false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OldGobbi::~OldGobbi() {
	delete FrontEcal;
	delete BackEcal;
	delete DeltaEcal;
	delete FrontTimecal;
	delete BackTimecal;
	delete DeltaTimecal;
	for (int i = 0; i < 4; i++) delete Silicon[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OldGobbi::analyze() {

	// Reset member variables
	for (size_t i = 0; i < 4; i++) Silicon[i]->reset();

	// Loop through hits, store in Silicon array and output to Histo
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
		
		//Use calibration to get Energy and fill elist class in silicon
		if (board == 1 || board == 3 || board == 5 || board == 7)
			addFrontHit((board - 1) / 2, chan, i);
		else if (board == 2 || board == 4 || board == 6 || board == 8)
			addBackHit((board / 2) - 1, chan, i);
		else if (board == 9 || board == 10 || board == 11 || board == 12)
			addDeltaHit(board - 9, chan, i);
		else {
			cout << "Invalid board number " << board << endl;
			cout << "Hit skipped!" << endl;
			return;
		}
	}

	// Loop through telescopes and perform additional analysis
	size_t sumchan = 0;
	size_t totMulti = 0;
	size_t Pidmulti = 0;
	for (size_t id = 0; id < 4; id++) {

		// This is the spot where we run Silicon->Neighbours() for addback
		Silicon[id]->Front.Neighbours(id);
		Silicon[id]->Back.Neighbours(id);
		Silicon[id]->Delta.Neighbours(id);

		// Then fill summary histograms after addback
		for (size_t n = 0; n < Silicon[id]->Front.Nstore; n++) {
			sumchan = id*hinpchans + Silicon[id]->Front.Order[n].strip;
			Histo.sumFrontE_addback->Fill(sumchan, Silicon[id]->Front.Order[n].energy);
		}
		for (size_t n = 0; n < Silicon[id]->Back.Nstore; n++) {
			sumchan = id*hinpchans + Silicon[id]->Back.Order[n].strip;
			Histo.sumBackE_addback->Fill(sumchan, Silicon[id]->Back.Order[n].energy);
		}
		for (size_t n = 0; n < Silicon[id]->Delta.Nstore; n++) {
			sumchan = id*hinpchans + Silicon[id]->Delta.Order[n].strip;
			Histo.sumDeltaE_addback->Fill(sumchan, Silicon[id]->Delta.Order[n].energy);
		}

		// Handle simple case of single strip multiplicity
		bool isSimple = false;
		if (Silicon[id]->Front.Nstore == 1 && Silicon[id]->Delta.Nstore == 1) {
			Histo.frontdeltastripnum[id]->Fill(Silicon[id]->Front.Order[0].strip, Silicon[id]->Delta.Order[0].strip);
			Histo.timediff[id]->Fill(Silicon[id]->Front.Order[0].time - Silicon[id]->Delta.Order[0].time);

			if (Silicon[id]->Back.Nstore == 1) {
				totMulti += Silicon[id]->simpleFront();
				Histo.sumFrontTimeMult1_cal->Fill(id*hinpchans + Silicon[id]->Front.Order[0].strip, Silicon[id]->Front.Order[0].time);
				isSimple = true;
			}
		}

		// If higher multiplicity then worry about picking the right one
		// This also handles the case where Nstore = 0 for any of the chanels
		if (!isSimple) totMulti += Silicon[id]->multiHit();

		// Next, fill E vs. dE and other plots

#ifdef ENABLE_DEBUG
		cout << "id " << id << ", Nsol " << Silicon[id]->Nsolution << endl;
#endif

		for (size_t isol = 0; isol < Silicon[id]->Nsolution; isol++) {
			Silicon[id]->position(isol); // calculates x,y pos, and lab angle

#ifdef ENABLE_DEBUG
			cout << "isol " << isol << endl; 
			cout << "cos " << cos(Silicon[id]->Solution[isol].theta) << endl;
			cout << "the " << Silicon[id]->Solution[isol].theta * rad_to_deg << endl;
			cout << "front strip " << Silicon[id]->Solution[isol].ifront << ", back strip " << Silicon[id]->Solution[isol].iback << endl;
			cout << "x " << Silicon[id]->Solution[isol].Xpos << ", y " << Silicon[id]->Solution[isol].Ypos << endl;
			cout << "E " << Silicon[id]->Solution[isol].energy << ", dE " << Silicon[id]->Solution[isol].denergy << endl;
#endif

			Histo.FrontvsBack[id]->Fill(Silicon[id]->Solution[isol].energy, Silicon[id]->Solution[isol].benergy);
			double Ener = Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy * (1 - cos(Silicon[id]->Solution[isol].theta));
			Histo.DEE[id]->Fill(Ener, Silicon[id]->Solution[isol].denergy * cos(Silicon[id]->Solution[isol].theta));
			Histo.xyhitmap->Fill(Silicon[id]->Solution[isol].Xpos, Silicon[id]->Solution[isol].Ypos);
			double th = Silicon[id]->Solution[isol].theta * rad_to_deg;
			Histo.Evstheta[id]->Fill(th, Silicon[id]->Solution[isol].energy);
			Histo.Evstheta_all->Fill(th, Silicon[id]->Solution[isol].energy);
			Histo.Theta->Fill(th);

			/******** ANGLE DEPENDENT CALIBRATION CORRECTION FOR HIGHER ENERGY POINTS ********/
			// This was added my Nicholas Dronchi, and can generally be ignored or bypassed
			// unless one finds that they need an angle dependent calibration. I (Henry Webb)
			// have left this in for now since all it does is fill some extra histograms.

			// NOTE: SpecTcl should already map channels to strips, and this code has been modified to skip the unpacking and directly read in a SpecTcl output tree
			int chan = Silicon[id]->Solution[isol].ifront;

			// Make a correction to the E silicon energy based on angle
			double angle_Ecorr = 1.0277e-5*(th*th*th) + 1.6125e-3*(th*th) + 8.3097e-4*th - 1.0227e-3;
			double Ecorr = Silicon[id]->Solution[isol].energy + angle_Ecorr;
			double Ecorr_R = FrontEcal->reverseCal(id, Silicon[id]->Solution[isol].ifront, Ecorr);

			Histo.AngleCorrE[id][chan]->Fill(Ecorr);
			Histo.AngleCorr_noCorr[id][chan]->Fill(Silicon[id]->Solution[isol].energy);
			Histo.AngleCorrE_R[id][chan]->Fill(Ecorr_R);
			Histo.AngleCorrFrontE_cal->Fill(id*hinpchans + Silicon[id]->Solution[isol].ifront, Ecorr);

			// NOTE: SpecTcl should already map channels to strips, and this code has been modified to skip the unpacking and directly read in a SpecTcl output tree
			int chandE = Silicon[id]->Solution[isol].ide;

			// Make a correction to the dE silicon energy based on angle
			double angle_dEcorr = -1.0971e-5*(th*th*th) - 1.1446e-3*(th*th) - 8.9371e-4*th + 1.0879e-3;
			//double angle_dEcorr = 2.0527e-6*(th*th*th) - 1.4281e-3*(th*th) + 1.5589e-4*th - 7.3389e-4; // no Au foilloss
			double dEcorr = Silicon[id]->Solution[isol].denergy + angle_dEcorr;
			double dEcorr_R = DeltaEcal->reverseCal(id,Silicon[id]->Solution[isol].ide, dEcorr);

			Histo.AngleCorrDeltaE[id][chandE]->Fill(dEcorr);
			Histo.AngleCorrDeltaE_noCorr[id][chandE]->Fill(Silicon[id]->Solution[isol].denergy);
			Histo.AngleCorrDeltaE_R[id][chandE]->Fill(dEcorr_R);
			Histo.AngleCorrDeltaE_cal->Fill(id*hinpchans + Silicon[id]->Solution[isol].ide, dEcorr);

			double Etot = Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy;
			Histo.sumEtot_cal->Fill(id*hinpchans + Silicon[id]->Solution[isol].ifront, Etot);
			Histo.AngleCorrSum_cal->Fill(id * hinpchans + Silicon[id]->Solution[isol].ifront, Ecorr + dEcorr);

#ifdef ENABLE_DEBUG
			cout << "th " << th << ", angle_Ecorr " << angle_Ecorr << " MeV, angle_dEcorr " << angle_dEcorr << " MeV" << endl;
			if (Ecorr > 5) cout << "EnergyR " << Silicon[id]->Solution[isol].energyR << ", Ecorr " << Ecorr << ", Ecorr_R " << Ecorr_R << endl;
#endif

		}

		// Calculate and determine particle identification (PID) in the silicon for each solution
		Pidmulti += Silicon[id]->getPID();

		// Hit maps and other plots based on Pid
		for (size_t isol = 0; isol < Silicon[id]->Nsolution; isol++) {
			solution& sol = Silicon[id]->Solution[isol];
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
		// Eloss files are loaded in silicon
		Silicon[id]->calcEloss();

		// Post Eloss calculation plots, I guess?
		for (size_t isol = 0; isol < Silicon[id]->Nsolution; isol++) {
			solution& sol = Silicon[id]->Solution[isol];
			if (sol.iZ == 1 && sol.iA == 1) {
				Histo.ProtonEnergy->Fill(sol.Ekin, sol.theta * rad_to_deg);
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t OldGobbi::loadSolutions(correl2& Correl) {
	size_t goodMult = 0;
	for (size_t id = 0; id < 4; id++) {
		for (size_t isol = 0; isol < Silicon[id]->Nsolution; isol++) {
			solution& sol = Silicon[id]->Solution[isol];

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

solution* OldGobbi::getNextEmptySolution(solution* sol) {
	bool foundEmptySolution = false;
	for (size_t id = 0; id < 4; id++) {
		silicon* s = Silicon[id];
		for (size_t isol = 0; isol < s->Nsolution; isol++) {
			if (s->Solution[isol].ipid != sol->ipid) continue;
			else if (s->Nsolution >= 20)
				throw invalid_argument("ERROR: max solutions reached, cannot get next empty solution");
			return &(s->Solution[s->Nsolution]);
		}
	}
	return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Pass through function for use with extra calculations external to this
// Gobbi class. All silicon losses should be the same, so just use the first
// one.
float OldGobbi::getEin(float energy, float thick, int Z, float A) {
	return Silicon[0]->losses->getEin(energy, thick, Z, A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OldGobbi::addFrontHit(size_t tel, size_t ch, size_t i) {
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
		Silicon[tel]->Front.Add(ch, Energy, ELoR, ER, time);
		Silicon[tel]->multFront++;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OldGobbi::addBackHit(size_t tel, size_t ch, size_t i) {
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
		Silicon[tel]->Back.Add(ch, Energy, ELoR, ER, time);
		Silicon[tel]->multBack++;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OldGobbi::addDeltaHit(size_t tel, size_t ch, size_t i) {
	size_t ER = input.GetE(i);
	size_t tR = input.GetT(i);
	size_t ELoR = input.GetELo(i);
	double Energy = DeltaEcal->getEnergy(tel, ch, ER);
	double time = DeltaTimecal->getTime(tel, ch, tR);

	Histo.sumDeltaE_R->Fill(tel*hinpchans + ch, ER);
	Histo.sumDeltaTime_R->Fill(tel*hinpchans + ch, tR);
	Histo.sumDeltaE_cal->Fill(tel*hinpchans + ch, Energy);
	Histo.sumDeltaTime_cal->Fill(tel*hinpchans + ch, time);

	Histo.DeltaE_R[tel][ch]->Fill(ER);
	Histo.DeltaElow_R[tel][ch]->Fill(ELoR);
	Histo.DeltaTime_R[tel][ch]->Fill(tR);
	Histo.DeltaE_cal[tel][ch]->Fill(Energy);

	//if(Energy > .2 && tR > 1765 && tR < 8600)
	if (Energy > .2) {
		//if (tel == 0 && ch == 0) cout << "EE " << Energy << endl;

		Silicon[tel]->Delta.Add(ch, Energy, ELoR, ER, time);
		Silicon[tel]->multDelta++;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



