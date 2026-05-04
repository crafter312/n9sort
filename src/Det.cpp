/* Nicolas Dronchi 2022_04_04
 * Class written to handle all specifics of the Gobbi array
 * such as communicating with HINP, calibrations, checking
 * for charge sharing in neighbor calculating geometry.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu), August 2025
 * Replaces use of HINP class for unpacking with Input class
 * for reading values from a SpecTcl-generated ROOT file.
 * This essentially offloads the work of unpacking to SpecTcl,
 * leaving this code to only do the analysis work.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu) and Johnathan
 * Phillips (j.s.phillips@wustl.edu) March 2026 for experiment
 * at TAMU Cyclotron Institute
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu) April 2026 for
 * 9N FRIB experiment
 */

#include "Det.h"

#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gobbi::Gobbi(Input& in, histo& hist, SortConfig& config, int run) : input(in.GetGobbi()), Histo(hist), input_qdc(in.GetQDC()),input_tdc(in.GetTDC()), runnum(run) {

	// Initialize wood class instances for ROOT TTree output
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.H3.mask[0] = 1;
	He4_pt = make_unique<wood>(Correl, "t_He4_pt", Histo.dir4He, false);
	Correl.zeroMask();
	Correl.H2.mask[0] = 1;
	Correl.H2.mask[1] = 1;
	He4_dd = make_unique<wood>(Correl, "t_He4_dd", Histo.dir4He, false);
	Correl.zeroMask();
	Correl.H2.mask[0] = 1;
	Correl.H3.mask[0] = 1;
	He5_dt = make_unique<wood>(Correl, "t_He5_dt", Histo.dir5He, false);
	Correl.zeroMask();
	Correl.H3.mask[0] = 1;
	Correl.H3.mask[1] = 1;
	He6_tt = make_unique<wood>(Correl, "t_He6_tt", Histo.dir6He, false);
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.alpha.mask[0] = 1;
	Li5_pa = make_unique<wood>(Correl, "t_Li5_pa", Histo.dir5Li, false);
	Correl.zeroMask();
	Correl.H2.mask[0] = 1;
	Correl.He3.mask[0] = 1;
	Li5_dHe3 = make_unique<wood>(Correl, "t_Li5_dHe3", Histo.dir5Li, false);
	Correl.zeroMask();
	Correl.H2.mask[0] = 1;
	Correl.He3.mask[0] = 1;
	Li6_da = make_unique<wood>(Correl, "t_Li6_da", Histo.dir6Li, false);
	Correl.zeroMask();
	Correl.H2.mask[0] = 0;
	Correl.H3_fake.mask[0] = 1;
	Correl.alpha.mask[0] = 1;
	Li7_ta_bad = make_unique<wood>(Correl, "t_Li7_ta_bad", Histo.dir7Li, false);
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.He6.mask[0] = 1;
	Li7_pHe6 = make_unique<wood>(Correl, "t_Li7_pHe6", Histo.dir7Li, false);
	Correl.zeroMask();
	Correl.H3.mask[0] = 1;
	Correl.alpha.mask[0] = 1;  
	Li7_ta = make_unique<wood>(Correl, "t_Li7_ta", Histo.dir7Li, false);
	Correl.zeroMask();
	Correl.alpha.mask[0] = 1;
	Correl.proton.mask[0] = 1;
	Correl.proton.mask[1] = 1;
	Be6_ppa = make_unique<wood>(Correl, "t_Be6_ppa", Histo.dir6Be, false);
	Correl.zeroMask();
	Correl.He3.mask[0] = 1;
	Correl.alpha.mask[0] = 1;
	Be7_He3a = make_unique<wood>(Correl, "t_Be7_He3a", Histo.dir7Be, false);
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.Li6.mask[0] = 1;
	Be7_pLi6 = make_unique<wood>(Correl, "t_Be7_pLi6", Histo.dir7Be, false);
	Correl.zeroMask();
	Correl.alpha.mask[0] = 1;
	Correl.alpha.mask[1] = 1;
	Be8_aa = make_unique<wood>(Correl, "t_Be8_aa", Histo.dir8Be, false);
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.Li7.mask[0] = 1;
	Be8_pLi7 = make_unique<wood>(Correl, "t_Be8_pLi7", Histo.dir8Be, false);
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.H3.mask[0] = 1;
	Correl.alpha.mask[0] = 1;
	Be8_pta = make_unique<wood>(Correl, "t_Be8_pta", Histo.dir8Be, false);
	Correl.zeroMask();
	Correl.H3.mask[0] = 1;
	Correl.alpha.mask[0] = 1;
	Li7_ta_fake = make_unique<wood>(Correl, "t_Li7_ta_fake", Histo.dir8Be, false); // missing the p
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.alpha.mask[0] = 1;
	Correl.alpha.mask[1] = 1;
	B9_paa = make_unique<wood>(Correl, "t_B9_paa", Histo.dir9B, false);

	Correl.zeroMask();

  Targetdist = config.GetTargDist();//23.95;//23.95;//24.1;//23.5; //cm //TODO is this correct? Shoud target dist be taken from input?
  TargetThickness = config.GetTargThick();;//3.2;//2.65; //mg/cm^2 for CD2 tar1 //TODO same as targ dist but for thickness
  //TargetThickness = 3.8; //mg/cm^2

  for (int id = 0; id < 4; id++) {
    Silicon[id] = new silicon(TargetThickness, config);
    Silicon[id]->init(id,config); //tells Silicon what position it is in
    Silicon[id]->SetTargetDistance(Targetdist);
  }

  string calDir = config.GetCalDir();
  FrontEcal = new calibrate(4, Histo.channum, calDir + config.GetFrontEcalFile(), 1, false);
  BackEcal = new calibrate(4, Histo.channum, calDir + config.GetBackEcalFile(), 1, false);
  DeltaEcal = new calibrate(4, Histo.channum, calDir + config.GetDeltaEcalFile(), 1, false);
  FrontTimecal = new calibrate(4, Histo.channum, calDir + config.GetFrontTimecalFile(),1, false);
  BackTimecal = new calibrate(4, Histo.channum, calDir + config.GetBackTimecalFile(),1, false);  
  DeltaTimecal = new calibrate(4, Histo.channum, calDir + config.GetDeltaTimecalFile(),1, false);
  
  //DiamondEcal = new calibrate(1, 4, calDir + config.GetDiamondEcalFile(), 1, false);
  
  // Diamond detector not use for this experiment, but this is where you could
	// perform a run-dependent diamond calibration selection if needed.
  diamond_calch = -1;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Gobbi::~Gobbi() {
  delete FrontEcal;
  delete BackEcal;
  delete DeltaEcal;
  delete FrontTimecal;
  delete BackTimecal;
  delete DeltaTimecal;
	for (int i = 0; i < 4; i++) delete Silicon[i]; // might not need this, but will try
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Gobbi::analyze() {
	
	//Vector holding calibrated diamond energies, same indices as input_qdc.qh
	//Clear each loop
	diamond_Ecal.clear();

	//Look at diamond detector
	for (int i=0;i<input_qdc.chan.size();i++) {
	
		//Get calibrated diamond energies
		diamond_Ecal.push_back(0);
		if ((diamond_calch > -1) && (DiamondEcal != nullptr)) diamond_Ecal[i] = (DiamondEcal->getEnergy(0, diamond_calch, input_qdc.qh[i]));
	
		if (input_qdc.chan[i] == 0) {
			Histo.DiamondQDC0->Fill(input_qdc.qh[i]);
			Histo.DiamondQDC0_cal->Fill(diamond_Ecal[i]);
			
			//Gated on or A tdc
			if (input_tdc.Nhits[1] != 0) {
				if (input_tdc.t[1][0] >= -80 && input_tdc.t[1][0] <= -50) {
					Histo.DiamondQDC0_tgate_orA->Fill(input_qdc.qh[i]);
					Histo.DiamondQDC0_tgate_orA_cal->Fill(diamond_Ecal[i]);
				}
				Histo.DiamondQDC0_vs_torA->Fill(input_qdc.qh[i],input_tdc.t[1][0]);
				Histo.DiamondQDC0_vs_torA_cal->Fill(diamond_Ecal[i],input_tdc.t[1][0]);
			}
		}
		if (input_qdc.chan[i] == 1) Histo.DiamondQDC1->Fill(input_qdc.qh[i]);
		if (input_qdc.chan[i] == 1) Histo.DiamondQDC1_cal->Fill(diamond_Ecal[i]);
	}

  // Reset the Silicon class
  //cout << "here pre Si reset" << endl;
  for (int i = 0; i < 4; i++) Silicon[i]->reset();
	//cout << "here post Si reset, have " << input.GetNhits() << " hits" << endl;
	size_t nhits = input.GetNhits();
  for (int i = 0; i < nhits; i++) {
    //check if there is an array for determined chip# and chan#
    if (input.GetBoard(i) > 12 || input.GetChan(i) >= Histo.channum)
    {
      cout << "Nhits " << input.GetNhits() << endl;
      cout << "i " << i << endl;
      cout << "Board " << input.GetBoard(i) << " and chan " << input.GetChan(i);
      cout << " unpacked but not saved" << endl;
      return true;
    }
	//cout << "here pre Si storing " << i << endl;
    float Energy = 0;
    float time = 0; //can be calibrated or shifted later
		
    //Use calibration to get Energy and fill elist class in silicon
    if (input.GetBoard(i) == 1 || input.GetBoard(i) == 3 || input.GetBoard(i) == 5 || input.GetBoard(i) == 7)
    {
      int quad = (input.GetBoard(i) - 1)/2;
      Energy = FrontEcal->getEnergy(quad, input.GetChan(i), input.GetE(i));
      time = FrontTimecal->getTime(quad, input.GetChan(i), input.GetT(i));

      Histo.sumFrontE_R->Fill(quad*Histo.channum + input.GetChan(i), input.GetE(i));
      Histo.sumFrontTime_R->Fill(quad*Histo.channum + input.GetChan(i), input.GetT(i));
      Histo.sumFrontE_cal->Fill(quad*Histo.channum + input.GetChan(i), Energy);
      Histo.sumFrontTime_cal->Fill(quad*Histo.channum + input.GetChan(i), time);

      Histo.FrontE_R[quad][input.GetChan(i)]->Fill(input.GetE(i));
      Histo.FrontElow_R[quad][input.GetChan(i)]->Fill(input.GetELo(i));
      Histo.FrontTime_R[quad][input.GetChan(i)]->Fill(input.GetT(i));
      Histo.FrontE_cal[quad][input.GetChan(i)]->Fill(Energy);      

      //if (Energy > .5 && input.GetT(i) > 3420 && input.GetT(i) < 6380 && (quad != 1 || Energy > 1.8))
      //need to set thresholds just above noise
      //if (quad == 1 && (

      if (Energy > .5) //(quad != 1 || Energy > 2)
      {
          Silicon[quad]->Front.Add(input.GetChan(i), Energy, input.GetELo(i), input.GetE(i), time);
          Silicon[quad]->multFront++;
      }
    }
    if (input.GetBoard(i) == 2 || input.GetBoard(i) == 4 || input.GetBoard(i) == 6 || input.GetBoard(i) == 8)
    {
      int quad = (input.GetBoard(i)/2)-1;
      Energy = BackEcal->getEnergy(quad, input.GetChan(i), input.GetE(i));
      time = BackTimecal->getTime(quad, input.GetChan(i), input.GetT(i));

      Histo.sumBackE_R->Fill(quad*Histo.channum + input.GetChan(i), input.GetE(i));
      Histo.sumBackTime_R->Fill(quad*Histo.channum + input.GetChan(i), input.GetT(i));
      Histo.sumBackE_cal->Fill(quad*Histo.channum + input.GetChan(i), Energy);
      Histo.sumBackTime_cal->Fill(quad*Histo.channum + input.GetChan(i), time);

      Histo.BackE_R[quad][input.GetChan(i)]->Fill(input.GetE(i));
      Histo.BackElow_R[quad][input.GetChan(i)]->Fill(input.GetELo(i));
      Histo.BackTime_R[quad][input.GetChan(i)]->Fill(input.GetT(i));
      Histo.BackE_cal[quad][input.GetChan(i)]->Fill(Energy);  

      if (Energy > .5)
      {      
          Silicon[quad]->Back.Add(input.GetChan(i), Energy, input.GetELo(i), input.GetE(i), time);
          Silicon[quad]->multBack++;
      }
    }
    if (input.GetBoard(i) == 9 || input.GetBoard(i) == 10 || input.GetBoard(i) == 11 || input.GetBoard(i) == 12)
    {
      int quad = (input.GetBoard(i)-9);
      Energy = DeltaEcal->getEnergy(quad, input.GetChan(i), input.GetE(i));
      time = DeltaTimecal->getTime(quad, input.GetChan(i), input.GetT(i));

      Histo.sumDeltaE_R->Fill(quad*Histo.channum + input.GetChan(i), input.GetE(i));
      Histo.sumDeltaTime_R->Fill(quad*Histo.channum + input.GetChan(i), input.GetT(i));
      Histo.sumDeltaE_cal->Fill(quad*Histo.channum + input.GetChan(i), Energy);
      Histo.sumDeltaTime_cal->Fill(quad*Histo.channum + input.GetChan(i), time);

      Histo.DeltaE_R[quad][input.GetChan(i)]->Fill(input.GetE(i));
      Histo.DeltaElow_R[quad][input.GetChan(i)]->Fill(input.GetELo(i));
      Histo.DeltaTime_R[quad][input.GetChan(i)]->Fill(input.GetT(i));
      Histo.DeltaE_cal[quad][input.GetChan(i)]->Fill(Energy);  

      //if(Energy > .2 && input.GetT(i) > 1765 && input.GetT(i) < 8600)
      if(Energy > .2)
      {
        //if (quad == 0 && input.GetChan(i) == 0) cout << "EE " << Energy << endl;
        
        Silicon[quad]->Delta.Add(input.GetChan(i), Energy, input.GetELo(i), input.GetE(i), time);
        Silicon[quad]->multDelta++;
      }
    }
  }
  //data is unpacked and stored into Silicon class at this point

	//cout << "here post Si storing" << endl;
  //This is the spot if we run Silicon->Neighbours()
  for (int id=0;id<4;id++) 
  {
    Silicon[id]->Front.Neighbours(id);
    Silicon[id]->Back.Neighbours(id);
    Silicon[id]->Delta.Neighbours(id);
  }

  //Fill summary after addback
  int sumchan = 0;
  for (int id=0;id<4;id++) 
  {
    for (int n=0; n<Silicon[id]->Front.Nstore; n++)
    {
      sumchan = id*Histo.channum + Silicon[id]->Front.Order[n].strip;
      Histo.sumFrontE_addback->Fill(sumchan, Silicon[id]->Front.Order[n].energy);
    }
    for (int n=0; n<Silicon[id]->Back.Nstore; n++)
    {
      sumchan = id*Histo.channum + Silicon[id]->Back.Order[n].strip;
      Histo.sumBackE_addback->Fill(sumchan, Silicon[id]->Back.Order[n].energy);
    }
    for (int n=0; n<Silicon[id]->Delta.Nstore; n++)
    {
      sumchan = id*Histo.channum + Silicon[id]->Delta.Order[n].strip;
      Histo.sumDeltaE_addback->Fill(sumchan, Silicon[id]->Delta.Order[n].energy);
    }
  }

  
  for (int id=0;id<4;id++) 
  {
    if (Silicon[id]->Front.Nstore == 1 && Silicon[id]->Delta.Nstore == 1)
    {
      Histo.frontdeltastripnum[id]->Fill(Silicon[id]->Front.Order[0].strip, Silicon[id]->Delta.Order[0].strip);
      Histo.timediff[id]->Fill(Silicon[id]->Front.Order[0].time - Silicon[id]->Delta.Order[0].time);
    }
  }

  //look at all the information/multiplicities and determine how to match up the information
  int totMulti = 0;
  for (int id=0;id<4;id++) 
  {
    //save comp time if the event is obvious
    if (Silicon[id]->Front.Nstore ==1 && Silicon[id]->Back.Nstore ==1 && Silicon[id]->Delta.Nstore ==1)
    {
      totMulti += Silicon[id]->simpleFront();
			Histo.sumFrontTimeMult1_cal->Fill(id*Histo.channum + Silicon[id]->Front.Order[0].strip, Silicon[id]->Front.Order[0].time);
      //cout << "simple " << totMulti << endl;
    }
    else //if higher multiplicity then worry about picking the right one
    {    //this also handles the case where Nstore=0 for any of the chanels
      totMulti += Silicon[id]->multiHit();
      //cout << "multi " << totMulti << endl;
    }
  }

	//Make crude time gates for the different telescopes (Same for each quad's OR A right now)
	//tGates_orA_upper[4] = {13723,};
	//tGates_orA_lower[4] = {};
	
	

  //plot E vs dE bananas and hitmap of paired dE,E events
  for (int id=0;id<4;id++) 
  {
    //cout << "id " << id << "   Nsol " << Silicon[id]->Nsolution << endl;
    for (int isol=0;isol<Silicon[id]->Nsolution; isol++)
    {
      //fill in hitmap of gobbi
      Silicon[id]->position(isol); //calculates x,y pos, and lab angle
      //cout << "isol " << isol << endl; 
      //cout << "cos " << cos(Silicon[id]->Solution[isol].theta) << endl;
      //cout << "the " << Silicon[id]->Solution[isol].theta*180./3.1415 << endl;
      //cout << "front strip " << Silicon[id]->Solution[isol].ifront << "  back strip " << Silicon[id]->Solution[isol].iback << endl;
      //cout << "x " << Silicon[id]->Solution[isol].Xpos << "  y " << Silicon[id]->Solution[isol].Ypos << endl;
      //cout << "E " << Silicon[id]->Solution[isol].energy << "  dE " << Silicon[id]->Solution[isol].denergy << endl;
      
      //Fill QDC vs Gobbi Esum plots
			for (int i=0;i<input_qdc.chan.size();i++) {
      	if (input_qdc.chan[i] == 0) {
      		Histo.Diamond_vs_GobbiEsum[id]->Fill(Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy,input_qdc.qh[i]);
      		Histo.Diamond_vs_GobbiEsum_cal[id]->Fill(Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy,diamond_Ecal[i]);
      		
      		
					//Gated on or A tdc
					if (input_tdc.Nhits[1] != 0) {
						if (input_tdc.t[1][0] >= -80 && input_tdc.t[1][0] <= -50) {
      				Histo.Diamond_vs_GobbiEsum_torA[id]->Fill(Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy,input_qdc.qh[i]);
      				Histo.Diamond_vs_GobbiEsum_torA_cal[id]->Fill(Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy,diamond_Ecal[i]);
						}
					}	
      		
      	}
      }
      
      //Front vs back plots
      Histo.FrontvsBack[id]->Fill(Silicon[id]->Solution[isol].energy,Silicon[id]->Solution[isol].benergy);

      //fill in dE-E plots to select particle type
      float Ener = Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy*(1-cos(Silicon[id]->Solution[isol].theta));

      Histo.DEE[id]->Fill(Ener, Silicon[id]->Solution[isol].denergy*cos(Silicon[id]->Solution[isol].theta));

      Histo.xyhitmap->Fill(Silicon[id]->Solution[isol].Xpos, Silicon[id]->Solution[isol].Ypos);
      //fill hist on theta/phi angles
      double th = Silicon[id]->Solution[isol].theta*180./pi; //in deg, not rad
      Histo.Evstheta[id]->Fill(th, Silicon[id]->Solution[isol].energy);
      Histo.Evstheta_all->Fill(th, Silicon[id]->Solution[isol].energy);
      Histo.Theta->Fill(th);

      //used to get high energy calibration points
      //need to convert interlaced strip numbers to calibration strip numbers
      // TODO already converted to strip #, make sure this isn't changed again
      int chan = Silicon[id]->Solution[isol].ifront;
      /*if (chan%2 ==0) //evens
      {
        chan = chan/2 + 16;
      }
      else //odds
      {
        chan = (chan -1)/2;
      }*/

      //make a correction to the energy based on angle
      float angle_Ecorr =  1.0277e-5*pow(th,3) + 1.6125e-3*pow(th,2) + 8.3097e-4*th - 1.0227e-3;
      //cout << "th " << th << " Angle corr " << angle_Ecorr << " MeV" << endl;
      float Ecorr = Silicon[id]->Solution[isol].energy + angle_Ecorr;
      float Ecorr_R = FrontEcal->reverseCal(id,Silicon[id]->Solution[isol].ifront, Ecorr);
      //if (Ecorr > 5)      
        //cout << "EnergyR " << Silicon[id]->Solution[isol].energyR << " Ecorr " << Ecorr << " Ecorr_R " << Ecorr_R << endl;

      Histo.AngleCorrE[id][chan]->Fill(Ecorr);
      Histo.AngleCorr_noCorr[id][chan]->Fill(Silicon[id]->Solution[isol].energy);
      Histo.AngleCorrE_R[id][chan]->Fill(Ecorr_R);
      
      Histo.AngleCorrFrontE_cal->Fill(id*Histo.channum + Silicon[id]->Solution[isol].ifront, Ecorr);

      int chandE = Silicon[id]->Solution[isol].ide;
      /*if (chan%2 ==0) //evens
      {
        chandE = chandE/2 + 16;
      }
      else //odds
      {
        chandE = (chandE -1)/2;
      }*/

      //make a correction to the energy based on angle
      float angle_dEcorr =  -1.0971e-5*pow(th,3) - 1.1446e-3*pow(th,2) - 8.9371e-4*th + 1.0879e-3;
      //float angle_dEcorr =  2.0527e-6*pow(th,3) - 1.4281e-3*pow(th,2) + 1.5589e-4*th - 7.3389e-4;//no Au foilloss
      //cout << "th " << th << " Angle corr " << angle_dEcorr << " MeV" << endl;
      float dEcorr = Silicon[id]->Solution[isol].denergy + angle_dEcorr;
      float dEcorr_R = DeltaEcal->reverseCal(id,Silicon[id]->Solution[isol].ide, dEcorr);

      Histo.AngleCorrDeltaE[id][chandE]->Fill(dEcorr);
      Histo.AngleCorrDeltaE_noCorr[id][chandE]->Fill(Silicon[id]->Solution[isol].denergy);
      Histo.AngleCorrDeltaE_R[id][chandE]->Fill(dEcorr_R);
      
      Histo.AngleCorrDeltaE_cal->Fill(id*Histo.channum + Silicon[id]->Solution[isol].ide, dEcorr);

      float Etot = Silicon[id]->Solution[isol].energy + Silicon[id]->Solution[isol].denergy;
      Histo.sumEtot_cal->Fill(id*Histo.channum + Silicon[id]->Solution[isol].ifront, Etot);
      Histo.AngleCorrSum_cal->Fill(id*Histo.channum + Silicon[id]->Solution[isol].ifront, Ecorr+dEcorr);

    }
  }

  //calculate and determine particle identification PID in the silicon
  int Pidmulti = 0;
  for (int id=0;id<4;id++) 
  {
    Pidmulti += Silicon[id]->getPID();
  }

  //hitmaps and other plots based on Pid
  for (int id=0;id<4;id++) 
  {
    for (int isol=0; isol<Silicon[id]->Nsolution; isol++)
    {
      float xpos = Silicon[id]->Solution[isol].Xpos;
      float ypos = Silicon[id]->Solution[isol].Ypos;
      
      // Gated on or A time
      if (input_tdc.t[1][0] >= -80 && input_tdc.t[1][0] <= -50) {
      	Histo.xyhitmap_tgate_orA->Fill(xpos,ypos);
      }
      
      //Gated on E-DE blobs
      if (Silicon[id]->Solution[isol].energy > 46 && Silicon[id]->Solution[isol].energy < 53 && Silicon[id]->Solution[isol].denergy > 3.6 && Silicon[id]->Solution[isol].denergy < 6.7) {
      	Histo.xyhitmap_EdEgate_1stEL->Fill(xpos,ypos);
      }
      if (Silicon[id]->Solution[isol].energy > 37 && Silicon[id]->Solution[isol].energy < 43 && Silicon[id]->Solution[isol].denergy > 5.5 && Silicon[id]->Solution[isol].denergy < 8.5) {
      	Histo.xyhitmap_EdEgate_2ndEL->Fill(xpos,ypos);
      }
      
      //Gated on QDC energy
      if (input_qdc.chan.size() > 0) {
      	for (int i=0;i<input_qdc.chan.size();i++) {
      		if (input_qdc.chan[i] == 0 && input_qdc.qh[i] < 2000) {
      		  Histo.xyhitmap_DiamondELlow->Fill(xpos,ypos);
      		}
      		if (input_qdc.chan[i] == 0 && input_qdc.qh[i] >= 2000 && input_qdc.qh[i] <= 2300) {
      		  Histo.xyhitmap_DiamondELpeak->Fill(xpos,ypos);
      		}
      		if (input_qdc.chan[i] == 0 && input_qdc.qh[i] > 2300) {
      		  Histo.xyhitmap_DiamondELhigh->Fill(xpos,ypos);
      		}
      	}
      }

      
      //protons
      if (Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 1)
      {
        Histo.protonhitmap->Fill(xpos, ypos);
        Histo.dTime_proton->Fill(Silicon[id]->Solution[isol].timediff);
      }
      //deuterons
      if (Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 2)
      {
        Histo.deuteronhitmap->Fill(xpos, ypos);
        Histo.dTime_deuteron->Fill(Silicon[id]->Solution[isol].timediff);
      }
      //tritons
      if (Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 3)
      {
        Histo.tritonhitmap->Fill(xpos, ypos);
        Histo.dTime_triton->Fill(Silicon[id]->Solution[isol].timediff);
      }
      //alphas
      if (Silicon[id]->Solution[isol].iZ == 2 && Silicon[id]->Solution[isol].iA == 4)
      {
        Histo.alphahitmap->Fill(xpos, ypos);
        Histo.dTime_alpha->Fill(Silicon[id]->Solution[isol].timediff);
      }
      //He6
      if (Silicon[id]->Solution[isol].iZ == 2 && Silicon[id]->Solution[isol].iA == 6)
      {
        Histo.He6hitmap->Fill(xpos, ypos);
        Histo.dTime_He6->Fill(Silicon[id]->Solution[isol].timediff);
      }
      if (Silicon[id]->Solution[isol].iZ == 3)
      {
        Histo.Lihitmap->Fill(xpos, ypos);
        Histo.dTime_Li->Fill(Silicon[id]->Solution[isol].timediff);
      }
      //Li veto
      if (Silicon[id]->Solution[isol].iZ != 3 && Silicon[id]->Solution[isol].iA != 7)
      {
        Histo.LiVETOhitmap->Fill(xpos, ypos);
      }
    }
  }

  //calc sumEnergy,then account for Eloss in target, then set Ekin and momentum of solutions
  //Eloss files are loaded in silicon
  for (int id=0;id<4;id++) 
  {
    Silicon[id]->calcEloss();
  }



  for (int id=0;id<4;id++) 
  {
    for (int isol=0; isol<Silicon[id]->Nsolution; isol++)
    {
      if (Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 1)
      {
        Histo.ProtonEnergy->Fill(Silicon[id]->Solution[isol].Ekin, Silicon[id]->Solution[isol].theta*180./pi);
      }
    }
  }






  //transfer Solution class to Correl
  Correl.reset();
  int goodMult = 0;
  for (int id=0;id<4;id++) 
  {
    for (int isol=0; isol<Silicon[id]->Nsolution; isol++)
    {
      //only keep solutions that have a pid
      if(Silicon[id]->Solution[isol].ipid)
      {
        Correl.load(&Silicon[id]->Solution[isol]);
        goodMult++;
      }
    }
  }

  if (goodMult >= 2)
  {

    //list all functions to look for correlations here
    corr_4He();
    corr_5He();
    corr_6He();
    corr_5Li();
    corr_6Li();
    corr_7Li();
    corr_6Be();
    corr_7Be();
    corr_8Be();
    corr_9B();
    
		//TODO Keep this or remove it? We should see the neutrons but maybe this is useful for just 6Li + p
    //lots of Li6 but they don't come in with anything. Could be Li6 + n    Keep correlation table?
    if (goodMult == 2)
    {
      //cout << "goodMult=2" << endl;
      int pos = 0;
      int particlenum[2] = {0,0};
      for (int id=0;id<4;id++)
      {
        //cout << "id" << id << endl;
        for (int isol=0; isol<Silicon[id]->Nsolution; isol++)
        {
          //cout << "  isol" << isol << endl;
          int pidnum = 0;
          if (Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 1)
            pidnum = 1;
          if (Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 2)
            pidnum = 2;
          if (Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 3)
            pidnum = 3;
          if (Silicon[id]->Solution[isol].iZ == 2 && Silicon[id]->Solution[isol].iA == 3)
            pidnum = 4;
          if (Silicon[id]->Solution[isol].iZ == 2 && Silicon[id]->Solution[isol].iA == 4)
            pidnum = 5;
          if (Silicon[id]->Solution[isol].iZ == 2 && Silicon[id]->Solution[isol].iA == 6)
            pidnum = 6;
          if (Silicon[id]->Solution[isol].iZ == 3 && Silicon[id]->Solution[isol].iA == 6)
            pidnum = 7;
          if (Silicon[id]->Solution[isol].iZ == 3 && Silicon[id]->Solution[isol].iA == 7)
            pidnum = 8;
          
          //cout << "    pidnum " << pidnum << endl;
          if (pidnum > 0)
          {
            //cout << "    ????pidnum " << pidnum << endl;
            particlenum[pos] = pidnum;
            pos++;
          }

        }
      }
      //cout << "After    1=" << particlenum[0] << "   2=" << particlenum[1] << endl;
      Histo.CorrelationTable->Fill(particlenum[0],particlenum[1]);
    }
    
  	//Count certain particle combinations
  	if (Correl.proton.mult == 1 && Correl.alpha.mult == 1) {
  		a_p++;
  	}
    
    
  }

  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Gobbi::match() {
  //match dE, Efront, Eback as one hit
  int temp = 0;
  return temp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

float Gobbi::getEnergy(int board, int chan, int Ehigh) {
  return Ehigh*(float)slopes[board][chan] + (float)intercepts[board][chan];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_4He() {
  // p+t
  if(Correl.proton.mult == 1 && Correl.H3.mult == 1) {
    float const Q4He = mass_alpha - (mass_p + mass_t);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.H3.mask[0]=1;
    Correl.makeArray(1, *He4_pt);

    float Erel_4He = Correl.findErel();
    float Ex = Erel_4He - Q4He;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

    Histo.Erel_4He_pt->Fill(Erel_4He);
    Histo.Ex_4He_pt->Fill(Ex);
    Histo.ThetaCM_4He_pt->Fill(thetaCM);
    Histo.VCM_4He_pt->Fill(Vcm);

    if (Erel_4He > 4.0)
    {
      float xpos = Correl.frag[0]->Xpos;
      float ypos = Correl.frag[0]->Ypos;
      Histo.He4_p_hitmap->Fill(xpos, ypos);
      xpos = Correl.frag[1]->Xpos;
      ypos = Correl.frag[1]->Ypos;
      Histo.He4_t_hitmap->Fill(xpos, ypos);

      //fill in dE-E plots to select particle type
      float Ener = Correl.frag[0]->energy + Correl.frag[0]->denergy*(1-cos(Correl.frag[0]->theta));

      Histo.DEE_He4[Correl.frag[0]->itele]->Fill(Ener, Correl.frag[0]->denergy*cos(Correl.frag[0]->theta));

      Ener = Correl.frag[1]->energy + Correl.frag[1]->denergy*(1-cos(Correl.frag[1]->theta));
  
      Histo.DEE_He4[Correl.frag[1]->itele]->Fill(Ener, Correl.frag[1]->denergy*cos(Correl.frag[1]->theta));
    }

    Histo.Erel_pt_costhetaH->Fill(Erel_4He, cos_thetaH);

		He4_pt->Fill(Erel_4He, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8); // Beam Z usually comes from S800, but S800 is not used for this experiment and we have a fixed beam type
  }

  // d+d
  if(Correl.H2.mult == 2) {
    float const Q4He = mass_alpha - (mass_d + mass_d);
    Correl.zeroMask();
    Correl.H2.mask[0]=1;
    Correl.H2.mask[1]=1;
    Correl.makeArray(1, *He4_dd);

    float Erel_4He = Correl.findErel();
    float Ex = Erel_4He - Q4He;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

    Histo.Erel_4He_dd->Fill(Erel_4He);
    Histo.Ex_4He_dd->Fill(Ex);
    Histo.ThetaCM_4He_dd->Fill(thetaCM);
    Histo.VCM_4He_dd->Fill(Vcm);

    Histo.Erel_dd_costhetaH->Fill(Erel_4He, cos_thetaH);

		He4_dd->Fill(Erel_4He, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_5He() {
  // d+t
  if(Correl.H2.mult == 1 && Correl.H3.mult == 1) {
    float const Q5He = mass_5He - (mass_d + mass_t);
    Correl.zeroMask();
    Correl.H2.mask[0]=1;
    Correl.H3.mask[0]=1;
    Correl.makeArray(1, *He5_dt);

    float Erel_5He = Correl.findErel();
    float Ex = Erel_5He - Q5He;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

    Histo.Erel_5He_dt->Fill(Erel_5He);
    Histo.Ex_5He_dt->Fill(Ex);
    Histo.ThetaCM_5He_dt->Fill(thetaCM);
    Histo.VCM_5He_dt->Fill(Vcm);

		Histo.Erel_dt_costhetaH->Fill(Erel_5He, cos_thetaH);

		He5_dt->Fill(Erel_5He, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_6He() {
  // t+t
  if(Correl.H3.mult == 2) {
    float const Q6He = mass_6He - (2*mass_t);
    Correl.zeroMask();
    Correl.H3.mask[0]=1;
    Correl.H3.mask[1]=1;
    Correl.makeArray(1, *He6_tt);

    float Erel_6He = Correl.findErel();
		float Ex = Erel_6He - Q6He;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

    Histo.Erel_6He_tt->Fill(Erel_6He);
    Histo.Ex_6He_tt->Fill(Ex);
    Histo.ThetaCM_6He_tt->Fill(thetaCM);
    Histo.VCM_6He_tt->Fill(Vcm);

		Histo.Erel_tt_costhetaH->Fill(Erel_6He, cos_thetaH);

		He6_tt->Fill(Erel_6He, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_5Li() {
  // p+He4
  if(Correl.proton.mult == 1 && Correl.alpha.mult == 1) {
 
    float const Q5Li = mass_5Li - (mass_p + mass_alpha);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.alpha.mask[0]=1;
    Correl.makeArray(1, *Li5_pa);

    float Erel_5Li = Correl.findErel();
		float Ex = Erel_5Li - Q5Li;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

    Histo.Erel_5Li_pa->Fill(Erel_5Li);
    Histo.Ex_5Li_pa->Fill(Ex);
    Histo.ThetaCM_5Li_pa->Fill(thetaCM);
    Histo.VCM_5Li_pa->Fill(Vcm);

		Histo.Erel_pa_costhetaH->Fill(Erel_5Li, cos_thetaH);

		Li5_pa->Fill(Erel_5Li, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);
  }

  // H2+He3
  if(Correl.H2.mult == 1 && Correl.He3.mult == 1) {

    float const Q5Li = mass_5Li - (mass_d + mass_3He);
    Correl.zeroMask();
    Correl.H2.mask[0]=1;
    Correl.He3.mask[0]=1;
    Correl.makeArray(1, *Li5_dHe3);

    float Erel_5Li = Correl.findErel();
		float Ex = Erel_5Li - Q5Li;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

    Histo.Erel_5Li_d3He->Fill(Erel_5Li);
    Histo.Ex_5Li_d3He->Fill(Ex);
    Histo.ThetaCM_5Li_d3He->Fill(thetaCM);
    Histo.VCM_5Li_d3He->Fill(Vcm);

		Histo.Erel_d3He_costhetaH->Fill(Erel_5Li, cos_thetaH);

		Li5_dHe3->Fill(Erel_5Li, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_6Li() {

  // D+alpha
  if(Correl.H2.mult == 1 && Correl.alpha.mult == 1) {
    //cout << "inside corr_6Li_da()" << endl;
 
    float const Q6Li = mass_6Li - (mass_d + mass_alpha);
    Correl.zeroMask();
    Correl.H2.mask[0]=1;
    Correl.alpha.mask[0]=1;
    Correl.makeArray(1, *Li6_da);

    float Erel_6Li = Correl.findErel();
		float Ex = Erel_6Li - Q6Li;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

    Histo.Erel_6Li_da->Fill(Erel_6Li);
    Histo.Ex_6Li_da->Fill(Ex);
    Histo.cos_da_thetaH->Fill(cos_thetaH);
    Histo.ThetaCM_6Li_da->Fill(thetaCM);
    Histo.VCM_6Li_da->Fill(Vcm);
    Histo.VCM_vs_ThetaCM->Fill(thetaCM, Vcm);
    Histo.Erel_da_cosThetaH->Fill(Erel_6Li, cos_thetaH);

		Li6_da->Fill(Erel_6Li, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

		// ToF calculations, added by Henry Webb (h.s.webb@wustl.edu)
		// This is used for quantifying neutron time resolution when
		// using Gobbi as the reference for neutron ToF calculations
		float deut_dist  = sqrt((Correl.frag[0]->Xpos*Correl.frag[0]->Xpos) + (Correl.frag[0]->Ypos*Correl.frag[0]->Ypos) + (Targetdist*Targetdist));
		float alpha_dist = sqrt((Correl.frag[1]->Xpos*Correl.frag[1]->Xpos) + (Correl.frag[1]->Ypos*Correl.frag[1]->Ypos) + (Targetdist*Targetdist));
		float deut_ToF   = deut_dist / Correl.frag[0]->velocity;
		float alpha_ToF  = alpha_dist / Correl.frag[1]->velocity;
		float time_deut  = deut_ToF + Correl.frag[0]->time; // Correl.frag[0].time is the E front time
		float time_alpha = alpha_ToF + Correl.frag[1]->time;
		float tdiff      = time_alpha - time_deut; // order shouldn't matter
		Histo.react_origin_tdiff->Fill(tdiff);

		//OR A gate
		if (input_tdc.t[1][0] >= -80 && input_tdc.t[1][0] <= -50) {
			Histo.Erel_6Li_da_tgate_orA->Fill(Erel_6Li);
			Histo.Erel_6Li_da_vsDiamond_tgate_orA->Fill(Erel_6Li,input_qdc.qh[0]);
		}

		for (int i=0;i<input_qdc.chan.size();i++) {
			if (input_qdc.chan[i] == 0) {
				//Sum energies
				float partsum = Correl.frag[0]->energy+Correl.frag[0]->denergy + Correl.frag[1]->energy+Correl.frag[1]->denergy;
			
				//vs diamond Q
				Histo.Erel_6Li_da_vsDiamond->Fill(Erel_6Li,input_qdc.qh[0]);
				Histo.Diamond_vs_GobbiEsum_cal_6Li[Correl.frag[0]->itele]->Fill(partsum,diamond_Ecal[i]);
				Histo.Diamond_vs_GobbiEsum_cal_6Li[Correl.frag[1]->itele]->Fill(partsum,diamond_Ecal[i]);

				//E*?
				float Ex_13C = 56 - diamond_Ecal[i] - partsum; // If 13C is not at all quenched, then this gives its excitation energy
				Histo.Diamond_Ex_6Li->Fill(Ex_13C);

				//OR A gate
				if (input_tdc.t[1][0] >= -80 && input_tdc.t[1][0] <= -50) {
					Histo.Diamond_vs_GobbiEsum_cal_6Li_torA[Correl.frag[0]->itele]->Fill(partsum,diamond_Ecal[i]);
					Histo.Diamond_vs_GobbiEsum_cal_6Li_torA[Correl.frag[1]->itele]->Fill(partsum,diamond_Ecal[i]);
					Histo.Diamond_Ex_6Li_torA->Fill(Ex_13C);
				}
				
				if (Ex > 2 && Ex < 2.4) {
				
					Histo.sumDiamond_vs_GobbiEsum_cal_6Li_3plus->Fill(partsum,diamond_Ecal[i]);
					Histo.Diamond_Ex_6Li_3plus->Fill(Ex_13C);
					
					//Gated on time or A
					if (input_tdc.t[1][0] >= -80 && input_tdc.t[1][0] <= -50) {
						Histo.sumDiamond_vs_GobbiEsum_cal_6Li_3plus_torA->Fill(partsum,diamond_Ecal[i]);
						Histo.Diamond_Ex_6Li_3plus_torA->Fill(Ex_13C);
					}
				}
				
			}
		}
    
    if(fabs(Correl.cos_thetaH) < .3)
      Histo.Ex_6Li_da_trans->Fill(Ex);
    if(fabs(Correl.cos_thetaH) > .7)
      Histo.Ex_6Li_da_long->Fill(Ex);

		//Hitmap gated on 3+ state
		if (Ex > 2 && Ex < 2.4) {
			//Light then heavy
			Histo.xyhitmap_6Li_plus->Fill(Correl.frag[0]->Xpos,Correl.frag[0]->Ypos);
			Histo.xyhitmap_6Li_plus->Fill(Correl.frag[1]->Xpos,Correl.frag[1]->Ypos);
		}

    if (Ex > 2 && Ex < 2.5)
    {
      //cout << "deut Ekin " << Correl.frag[0]->Ekin << endl;
      //cout << "alpha Ekin " << Correl.frag[1]->Ekin << endl;
      Histo.deutE_gate->Fill(Correl.frag[0]->Ekin);
      Histo.alphaE_gated->Fill(Correl.frag[1]->Ekin);

      Histo.deutE_gate_cosThetaH->Fill(Correl.frag[0]->Ekin,Correl.cos_thetaH);
      Histo.alphaE_gate_cosThetaH->Fill(Correl.frag[1]->Ekin,Correl.cos_thetaH);
    }

		/******** BAD 7LI->a+t ********/
		// Reconstructing a+d as a+t, maybe to see if a triton was missidentified
		// as a deuteron or something?

    for (int id=0;id<4;id++) 
    {
      for (int isol=0; isol<Silicon[id]->Nsolution; isol++)
      {
        //pick out
        if(Silicon[id]->Solution[isol].iZ == 1 && Silicon[id]->Solution[isol].iA == 2)
        {
          Silicon[id]->Solution[isol].iZ = 1;
          Silicon[id]->Solution[isol].iA = 3; //whoops now it is a triton
          Silicon[id]->Solution[isol].mass = Mass_t;

          //now we should redo energy loss calc          

          float sumEnergy = Silicon[id]->Solution[isol].denergy + Silicon[id]->Solution[isol].energy;
          float pc_before = sqrt(pow(sumEnergy+Silicon[id]->Solution[isol].mass,2) - pow(Silicon[id]->Solution[isol].mass,2));
          float velocity_before = pc_before/(sumEnergy+Silicon[id]->Solution[isol].mass);

          float thick = TargetThickness/2/cos(Silicon[id]->Solution[isol].theta);

          float ein = Silicon[id]->losses->getEin(sumEnergy,thick,Silicon[id]->Solution[isol].iZ,Silicon[id]->Solution[isol].mass/m0);

          Silicon[id]->Solution[isol].Ekin = ein;
          //calc momentum vector, energyTot, and velocity
          Silicon[id]->Solution[isol].getMomentum();

          //reload the solution in and it will be a missidentified triton
          //particle[3] is the H3_fake particle
          Correl.particle[3]->Sol[0] = &Silicon[id]->Solution[isol];
          Correl.particle[3]->mult++;
        }
      }
    }
    
    float const Q7Li = mass_7Li - (mass_t + mass_alpha);
    Correl.zeroMask();
    Correl.H2.mask[0]=0;
    Correl.H3_fake.mask[0]=1;
    Correl.alpha.mask[0]=1;   
    Correl.makeArray(1, *Li7_ta_bad);

		float Erel_7Li = Correl.findErel();
		float Ex_7Li = Erel_7Li - Q7Li;
		float Vcm_7Li = Correl.velocityCM;
    float thetaCM_7Li = Correl.thetaCM*180./acos(-1);
		float cos_thetaH_7Li = Correl.cos_thetaH;

    Histo.Ex_7Li_ta_bad->Fill(Ex_7Li);

		Li7_ta_bad->Fill(Erel_7Li, Ex_7Li, Vcm_7Li, thetaCM_7Li, cos_thetaH_7Li, runnum, 8);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_7Li() {
  // p + 6He
  if(Correl.proton.mult == 1 && Correl.He6.mult == 1) {
    //cout << "inside corr_7Li_he6p()" << endl;
 
    float const Q7Li = mass_7Li - (mass_p + mass_6He);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.He6.mask[0]=1;   
    Correl.makeArray(1, *Li7_pHe6);

    float Erel_7Li = Correl.findErel();
    float Ex = Erel_7Li - Q7Li;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

		Li7_pHe6->Fill(Erel_7Li, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    float getqvalue = Correl.Qvalue();

    //angle between beam axis and the momentum vector between the center
    // of mass and the heavy fragment
    //float mag2 = 0;
    //for (int j=0;j<3;j++) mag2 += pow(Correl.frag[1]->MomCM[j],2);
    //float cosbeamCMtoHF = Correl.frag[1]->MomCM[2]/sqrt(mag2);
    //cout << "cosbeamCMtoHF " << cosbeamCMtoHF << endl;

    Histo.Erel_7Li_p6He->Fill(Erel_7Li);
    Histo.Erel_7Li_p6He_Q->Fill(Erel_7Li, getqvalue);
    Histo.Erel_7Li_p6He_lowres->Fill(Erel_7Li);

    if (abs(cos_thetaH) < 0.7) {
      Histo.Ex_7Li_p6He_transverse->Fill(Ex);
    }
    if (cos_thetaH > -0.7) {
      Histo.Ex_7Li_p6He_transverse2->Fill(Ex);
    }
    else  Histo.Erel_7Li_p6He_pFor->Fill(Erel_7Li);
  }

  // H3 + He4
  if(Correl.H3.mult == 1 && Correl.alpha.mult == 1) {
 
    float const Q7Li = mass_7Li - (mass_t + mass_alpha);
    Correl.zeroMask();
    Correl.H3.mask[0]=1;
    Correl.alpha.mask[0]=1;   
    Correl.makeArray(1, *Li7_ta);

    float Erel_7Li = Correl.findErel();
    float Ex = Erel_7Li - Q7Li;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

		Li7_ta->Fill(Erel_7Li, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_7Li_ta->Fill(Erel_7Li);
    Histo.Ex_7Li_ta->Fill(Ex);
    Histo.ThetaCM_7Li_ta->Fill(thetaCM);
    Histo.VCM_7Li_ta->Fill(Vcm);

    if(fabs(cos_thetaH) < .5)
      Histo.Ex_7Li_ta_trans->Fill(Ex);
    if(fabs(cos_thetaH) > .7)
      Histo.Ex_7Li_ta_long->Fill(Ex);

		float Ex_tar = Correl.TargetEx();

    Histo.cos_ta_thetaH->Fill(cos_thetaH);
    Histo.Erel_ta_cosThetaH->Fill(Erel_7Li, cos_thetaH);
    Histo.Ex_tar->Fill(Ex_tar);
    Histo.Erel_vs_Extar->Fill(Erel_7Li, Ex_tar);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_6Be() {
  // He4+p+p
  if(Correl.alpha.mult == 1 && Correl.proton.mult == 2) {
		float const Q6Be = mass_6Be - (mass_alpha + 2.*mass_p);
    Correl.zeroMask();
    Correl.alpha.mask[0]=1;
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.makeArray(1, *Be6_ppa);

    float Erel_6Be = Correl.findErel();
		float Ex = Erel_6Be - Q6Be;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

		Be6_ppa->Fill(Erel_6Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_6Be_2pa->Fill(Erel_6Be);
    Histo.ThetaCM_6Be_2pa->Fill(thetaCM);
    Histo.VCM_6Be_2pa->Fill(Vcm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_7Be() {
  // He3+He4
  if(Correl.He3.mult == 1 && Correl.alpha.mult == 1) {
    float const Q7Be = mass_7Be - (mass_3He + mass_alpha);
    Correl.zeroMask();
    Correl.He3.mask[0]=1;
    Correl.alpha.mask[0]=1;
    Correl.makeArray(1, *Be7_He3a);

		float Erel_7Be = Correl.findErel();
		float Ex = Erel_7Be - Q7Be;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

		Be7_He3a->Fill(Erel_7Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_7Be_a3He->Fill(Erel_7Be);
    Histo.Ex_7Be_a3He->Fill(Ex);
    Histo.ThetaCM_7Be_a3He->Fill(thetaCM);
    Histo.VCM_7Be_a3He->Fill(Vcm);
  }
  // p+Li6
  if(Correl.proton.mult == 1 && Correl.Li6.mult == 1) {
    float const Q7Be = mass_7Be - (mass_p + mass_6Li);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Li6.mask[0]=1;
    Correl.makeArray(1, *Be7_pLi6);

    float Erel_7Be = Correl.findErel();
    float Ex = Erel_7Be - Q7Be;
		float Vcm = Correl.velocityCM;
		float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

		Be7_pLi6->Fill(Erel_7Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_7Be_p6Li->Fill(Erel_7Be);
    Histo.Ex_7Be_p6Li->Fill(Ex);
    Histo.ThetaCM_7Be_p6Li->Fill(thetaCM);
    Histo.VCM_7Be_p6Li->Fill(Vcm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_8Be() {
  // He4+He4
  if(Correl.alpha.mult == 2) {
    //cout << "inside corr_8Be_aa()" << endl;
 
    float const Q8Be = mass_8Be - (2*mass_alpha);
    Correl.zeroMask();
    Correl.alpha.mask[0]=1;
    Correl.alpha.mask[1]=1;
    Correl.makeArray(1, *Be8_aa);

    float Erel_8Be = Correl.findErel();
		float Ex = Erel_8Be - Q8Be;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
    float cos_thetaH = Correl.cos_thetaH;

		Be7_pLi6->Fill(Erel_8Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_8Be_aa->Fill(Erel_8Be);
    Histo.Ex_8Be_aa->Fill(Ex);

    Histo.Erel_aa_cosThetaH->Fill(Erel_8Be, cos_thetaH);
    if (Ex < 0.1) {
      Histo.ThetaCM_8Be_aa->Fill(thetaCM);
      Histo.VCM_8Be_aa->Fill(Vcm);
    }

  }

  // p + Li7
  if(Correl.proton.mult == 1 && Correl.Li7.mult == 1) {
    //cout << "inside corr_6Li_da()" << endl;
 
    float const Q8Be = mass_8Be - (mass_p + mass_7Li);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Li7.mask[0]=1;
    Correl.makeArray(1, *Be8_pLi7);

    float Erel_8Be = Correl.findErel();
		float Ex = Erel_8Be - Q8Be;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
    float cos_thetaH = Correl.cos_thetaH;

		Be8_pLi7->Fill(Erel_8Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_8Be_p7Li->Fill(Erel_8Be);
    Histo.Ex_8Be_p7Li->Fill(Ex);
    if(fabs(cos_thetaH) < .5)
      Histo.Ex_8Be_p7Li_trans->Fill(Ex);

    Histo.ThetaCM_8Be_p7Li->Fill(thetaCM);
    Histo.VCM_8Be_p7Li->Fill(Vcm);
    Histo.cos_p7Li_thetaH->Fill(cos_thetaH);
    Histo.Erel_p7Li_cosThetaH->Fill(Erel_8Be, cos_thetaH);
  }

  //p+t+a
  if(Correl.proton.mult == 1 && Correl.H3.mult == 1 && Correl.alpha.mult == 1) {
    float const Q8Be = mass_8Be - (mass_p + mass_t + mass_alpha);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.H3.mask[0]=1;
    Correl.alpha.mask[0]=1;
    Correl.makeArray(1, *Be8_pta);

    float Erel_8Be = Correl.findErel();
		float Ex = Erel_8Be - Q8Be;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
		float cos_thetaH = Correl.cos_thetaH;

		Be8_pta->Fill(Erel_8Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_8Be_pta->Fill(Erel_8Be);
    Histo.Ex_8Be_pta->Fill(Ex);
    if(fabs(cos_thetaH) < .5)
      Histo.Ex_8Be_pta_trans->Fill(Ex);
    Histo.ThetaCM_8Be_pta->Fill(thetaCM);
    Histo.VCM_8Be_pta->Fill(Vcm);

    Histo.cos_pta_thetaH->Fill(cos_thetaH);
    Histo.Erel_pta_cosThetaH->Fill(Erel_8Be, cos_thetaH);
  }

  //t+a (missing p)
  if(Correl.proton.mult == 1 && Correl.H3.mult == 1 && Correl.alpha.mult == 1) {
    float const Q7Li = mass_7Li - (mass_t + mass_alpha);
    Correl.zeroMask();
    Correl.H3.mask[0]=1;
    Correl.alpha.mask[0]=1;   
    Correl.makeArray(1, *Li7_ta_fake);

    float Erel_7Li = Correl.findErel();
		float Ex = Erel_7Li - Q7Li;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
    float cos_thetaH = Correl.cos_thetaH;

		Li7_ta_fake->Fill(Erel_7Li, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_7Li_ta_fake->Fill(Erel_7Li);
    Histo.Ex_7Li_ta_fake->Fill(Ex);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Gobbi::corr_9B() {
  // He4+He4+p
  if(Correl.alpha.mult == 2 && Correl.proton.mult == 1) {
    float const Q9B = mass_9B - (2*mass_alpha) - mass_p;
    Correl.zeroMask();
    Correl.alpha.mask[0]=1;
    Correl.alpha.mask[1]=1;
    Correl.proton.mask[0]=1;
    Correl.makeArray(1, *B9_paa);

    float Erel_9B = Correl.findErel();
		float Ex = Erel_9B - Q9B;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*180./acos(-1);
    float cos_thetaH = Correl.cos_thetaH;

		B9_paa->Fill(Erel_9B, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_9B_paa->Fill(Erel_9B);
    Histo.Ex_9B_paa->Fill(Ex);
    Histo.ThetaCM_9B_paa->Fill(thetaCM);
    Histo.VCM_9B_paa->Fill(Vcm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



