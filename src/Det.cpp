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

Det::Det(Input& in, histo& hist, SortConfig& config, size_t run) : input(in.GetGobbi()), Targetdist(config.GetTargDist()), TargetThickness(config.GetTargThick()), hinpboards(config.GetHinpboards()), hinpchans(config.GetHinpchans()), Histo(hist), gobbi(in, hist, config), input_tdc(in.GetTDC()), runnum(run) {

#ifdef ENABLE_DEBUG
	cout << "Det::Det 1" << endl;
#endif

	// Initialize wood class instances for ROOT TTree output
	Correl.zeroMask();
	Correl.proton.mask[0] = 1;
	Correl.H3.mask[0] = 1;
	He4_pt = make_unique<wood>(Correl, "t_He4_pt", Histo.dir4He, false);
	
#ifdef ENABLE_DEBUG
	cout << "Det::Det 2" << endl;
#endif
	
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
	
#ifdef ENABLE_DEBUG
	cout << "Det::Det 3" << endl;
#endif
	
	Correl.zeroMask();
	Correl.proton.mask[0]=1;
	Correl.proton.mask[1]=1;
	Correl.proton.mask[2]=1;
	Correl.proton.mask[3]=1;
	Correl.alpha.mask[0]=1;
	C8_4pa = make_unique<wood>(Correl, "t_C8_4pa", Histo.dir8C, false);
	Correl.zeroMask();
	Correl.proton.mask[0]=1;
	Correl.proton.mask[1]=1;
	Correl.proton.mask[2]=1;
	Correl.proton.mask[3]=1;
	Correl.proton.mask[4]=1;
	Correl.alpha.mask[0]=1;
	N9_5pa = make_unique<wood>(Correl, "t_N9_5pa", Histo.dir9N, false);
	
#ifdef ENABLE_DEBUG
	cout << "Det::Det 4" << endl;
#endif

	Correl.zeroMask();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Det::~Det() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::analyze() {

	// Perform Gobbi specific analysis
	gobbi.analyze();

	//transfer Solution classes from detectors to Correl
	Correl.reset();
	size_t goodMult = gobbi.loadSolutions(Correl);
	if (goodMult < 2) return;

	// List all functions to look for correlations here
	//corr_4He();
	//corr_5He();
	//corr_6He();
	corr_5Li();
	//corr_6Li();
	//corr_7Li();
	corr_6Be();
	//corr_7Be();
	corr_8Be();
	//corr_9B();
	corr_8C();
	corr_9N();

	if (goodMult == 2) {
		size_t pos = 0;
		size_t particlenum[2] = { 0, 0 };
		for (size_t i = 0; i < Correl.Nparticles; i++) {
			if (i == 4) continue; // hard coded to skip H3_fake particles created in corr_6Li (TODO this should be done for all fake particles)
			for (size_t j = 0; j < Correl.particle[i]->mult; j++) {
				if (pos > 1) {
					cout << "WARNING: more than two solutions for goodMult == 2, using first two" << endl;
					break;
				}
				particlenum[pos] = Correl.particle[i]->Sol[j]->ipid;
				pos++;
			}
		}

		Histo.CorrelationTable->Fill(particlenum[0], particlenum[1]);

#ifdef ENABLE_DEBUG
		cout << "After\t1 = " << particlenum[0] << ", 2 = " << particlenum[1] << endl;
#endif

	}
	
	// Count certain particle combinations
	if ((Correl.alpha.mult == 1) && (goodMult - Correl.proton.mult - 1 == 0)) {
		switch (Correl.proton.mult) {
			case 1:
				a_p++;
				break;
			case 2:
				a_pp++;
				break;
			case 3:
				a_ppp++;
				break;
			case 4:
				a_pppp++;
				break;
			case 5:
				a_ppppp++;
				break;
		}
	}
	else if (Correl.alpha.mult == 2) {
		if ((goodMult == 6) && (Correl.proton.mult == 4))
			aa_pppp++;
		else if ((goodMult == 3) && (Correl.He3.mult == 1))
			aa_3He++;
	}
	
	// Other counters
	pidSkipped += gobbi.GetPidSkipped();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::corr_4He() {
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
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

void Det::corr_5He() {
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
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

void Det::corr_6He() {
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
    float thetaCM = Correl.thetaCM*rad_to_deg;
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

void Det::corr_5Li() {
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
    float thetaCM = Correl.thetaCM*rad_to_deg;
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
    float thetaCM = Correl.thetaCM*rad_to_deg;
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

void Det::corr_6Li() {

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
    float thetaCM = Correl.thetaCM*rad_to_deg;
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
/*
		solution* fakesol = gobbi.getNextEmptySolution(Correl.H2.Sol[0]);                  // this gets a solution we know is not being used
		if (fakesol == nullptr) {                                                          // this should never trigger, but handle nullptr case to be safe
			cout << "WARNING: next empty solution not found! Skipping bad 7Li..." << endl;
			return;
		}
		*fakesol = *(Correl.H2.Sol[0]);                                                    // make sure to copy deuteron's information to this new solution
    fakesol->iA = 3;                                                                   // whoops the copied deuteron is now a triton
    fakesol->mass = Mass_t;

    // Now we should redo energy loss calc          
    double sumEnergy = fakesol->denergy + fakesol->energy;
    double pc_before = sqrt(pow(sumEnergy + fakesol->mass, 2.) - (fakesol->mass*fakesol->mass));
    double velocity_before = pc_before / (sumEnergy + fakesol->mass);
    double thick = (.5 * TargetThickness) / cos(fakesol->theta);
    double ein = gobbi.getEin(sumEnergy, thick, fakesol->iZ, fakesol->mass / m0);

		// Calculate momentum vector, energyTot, and velocity
    fakesol->Ekin = ein;
    fakesol->getMomentum();

    // Load the new solution in and it will be a missidentified triton
    Correl.H3_fake.Sol[0] = fakesol;
    Correl.H3_fake.mult++;
    
    float const Q7Li = mass_7Li - (mass_t + mass_alpha);
    Correl.zeroMask();
    Correl.H2.mask[0]=0;
    Correl.H3_fake.mask[0]=1;
    Correl.alpha.mask[0]=1;   
    Correl.makeArray(1, *Li7_ta_bad);

		float Erel_7Li = Correl.findErel();
		float Ex_7Li = Erel_7Li - Q7Li;
		float Vcm_7Li = Correl.velocityCM;
    float thetaCM_7Li = Correl.thetaCM*rad_to_deg;
		float cos_thetaH_7Li = Correl.cos_thetaH;

    Histo.Ex_7Li_ta_bad->Fill(Ex_7Li);

		Li7_ta_bad->Fill(Erel_7Li, Ex_7Li, Vcm_7Li, thetaCM_7Li, cos_thetaH_7Li, runnum, 8);
*/
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::corr_7Li() {
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
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

void Det::corr_6Be() {
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
		float cos_thetaH = Correl.cos_thetaH;

		Be6_ppa->Fill(Erel_6Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_6Be_2pa->Fill(Erel_6Be);
    Histo.ThetaCM_6Be_2pa->Fill(thetaCM);
    Histo.VCM_6Be_2pa->Fill(Vcm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::corr_7Be() {
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
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
		float thetaCM = Correl.thetaCM*rad_to_deg;
		float cos_thetaH = Correl.cos_thetaH;

		Be7_pLi6->Fill(Erel_7Be, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_7Be_p6Li->Fill(Erel_7Be);
    Histo.Ex_7Be_p6Li->Fill(Ex);
    Histo.ThetaCM_7Be_p6Li->Fill(thetaCM);
    Histo.VCM_7Be_p6Li->Fill(Vcm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::corr_8Be() {
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
    float thetaCM = Correl.thetaCM*rad_to_deg;
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
    float thetaCM = Correl.thetaCM*rad_to_deg;
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
    float thetaCM = Correl.thetaCM*rad_to_deg;
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
    float thetaCM = Correl.thetaCM*rad_to_deg;
    float cos_thetaH = Correl.cos_thetaH;

		Li7_ta_fake->Fill(Erel_7Li, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_7Li_ta_fake->Fill(Erel_7Li);
    Histo.Ex_7Li_ta_fake->Fill(Ex);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::corr_9B() {
  // He4+He4+p
  if (Correl.alpha.mult == 2 && Correl.proton.mult == 1) {
    float const Q9B = mass_9B - (2*mass_alpha) - mass_p;
    Correl.zeroMask();
    Correl.alpha.mask[0]=1;
    Correl.alpha.mask[1]=1;
    Correl.proton.mask[0]=1;
    Correl.makeArray(1, *B9_paa);

    float Erel_9B = Correl.findErel();
		float Ex = Erel_9B - Q9B;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*rad_to_deg;
    float cos_thetaH = Correl.cos_thetaH;

		B9_paa->Fill(Erel_9B, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_9B_paa->Fill(Erel_9B);
    Histo.Ex_9B_paa->Fill(Ex);
    Histo.ThetaCM_9B_paa->Fill(thetaCM);
    Histo.VCM_9B_paa->Fill(Vcm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::corr_8C() {
	// p+p+p+p+a
	if (Correl.proton.mult == 4 && Correl.alpha.mult == 1) {
		float const Q8C = mass_8C - (4*mass_p) - mass_alpha;
		Correl.zeroMask();
		Correl.proton.mask[0]=1;
		Correl.proton.mask[1]=1;
		Correl.proton.mask[2]=1;
		Correl.proton.mask[3]=1;
		Correl.alpha.mask[0]=1;
		Correl.makeArray(1, *C8_4pa);

		float Erel_8C = Correl.findErel();
		float Ex = Erel_8C - Q8C;
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*rad_to_deg;
    float cos_thetaH = Correl.cos_thetaH;

		C8_4pa->Fill(Erel_8C, Ex, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_8C_4pa->Fill(Erel_8C);
    Histo.Ex_8C_4pa->Fill(Ex);
    Histo.ThetaCM_8C_4pa->Fill(thetaCM);
    Histo.VCM_8C_4pa->Fill(Vcm);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Det::corr_9N() {
	// p+p+p+p+p+a
	if (Correl.proton.mult == 5 && Correl.alpha.mult == 1) {
		Correl.zeroMask();
		Correl.proton.mask[0]=1;
		Correl.proton.mask[1]=1;
		Correl.proton.mask[2]=1;
		Correl.proton.mask[3]=1;
		Correl.proton.mask[4]=1;
		Correl.alpha.mask[0]=1;
		Correl.makeArray(1, *N9_5pa);

		float Erel_9N = Correl.findErel();
		float Vcm = Correl.velocityCM;
    float thetaCM = Correl.thetaCM*rad_to_deg;
    float cos_thetaH = Correl.cos_thetaH;

		// No mass excess for 9N, so no Q value and no excitation energy
		C8_4pa->Fill(Erel_9N, -1, Vcm, thetaCM, cos_thetaH, runnum, 8);

    Histo.Erel_9N_5pa->Fill(Erel_9N);
    Histo.ThetaCM_9N_5pa->Fill(thetaCM);
    Histo.VCM_9N_5pa->Fill(Vcm);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



