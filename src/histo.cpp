#include "histo.h"

#include <vector>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

histo::histo(shared_ptr<ROOT::TBufferMergerFile> f, SortConfig& config) : hinpboards(config.GetHinpboards()), hinpchans(config.GetHinpchans()) {
  file_read = f;
  file_read->cd();

  // Create global tree for storing pre-solution variables
  const char* otreename = config.GetOtreeName().c_str();
  //tpar = new TTree(otreename, otreename);
  //tpar->Branch("solutions", &solutions);

  //// Create subdirectories to store arrays of spectra

  // Directory for summary plots
  dirSummary = new TDirectoryFile("Summary", "Summary"); // name, title

  // Energies, Raw+Calibrated
  dir1dFrontE_R = dirSummary->mkdir("1dFrontE_R","1dFrontE_R");
  dir1dBackE_R  = dirSummary->mkdir("1dBackE_R","1dBackE_R");
  dir1dDeltaE_R = dirSummary->mkdir("1dDeltaE_R","1dDeltaE_R");

  dir1dFrontlowE_R = dirSummary->mkdir("1dFrontlowE_R","1dFrontlowE_R");
  dir1dBacklowE_R  = dirSummary->mkdir("1dBacklowE_R","1dBacklowE_R");
  dir1dDeltalowE_R = dirSummary->mkdir("1dDeltalowE_R","1dDeltalowE_R");

  dir1dFrontE_cal = dirSummary->mkdir("1dFrontE_cal","1dFrontE_cal");
  dir1dBackE_cal  = dirSummary->mkdir("1dBackE_cal","1dBackE_cal");
  dir1dDeltaE_cal = dirSummary->mkdir("1dDeltaE_cal","1dDeltaE_cal");

  dirAngleCorrFrontE = dirSummary->mkdir("AngleCorrFrontE","AngleCorrFrontE");
  dirAngleCorrDeltaE = dirSummary->mkdir("AngleCorrDeltaE","AngleCorrDeltaE");

  // Delta Energies, Raw+Calibrated
  dir1dFrontTime_R = dirSummary->mkdir("1dFrontTime_R","1dFrontTime_R");
  dir1dBackTime_R  = dirSummary->mkdir("1dBackTime_R","1dBackTime_R");
  dir1dDeltaTime_R = dirSummary->mkdir("1dDeltaTime_R","1dDeltaTime_R");
  
  // CsI directories
  dir1dCsI_Energy = dirSummary->mkdir("1dCsI_Energy","1dCsI_Energy");
  dir1dCsI_Time   = dirSummary->mkdir("1dCsI_Time","1dCsI_Time");
  dir1dCsI_QDC    = dirSummary->mkdir("1dCsI_QDC","1dCsI_QDC");

  // Directory for TDC detector
  dirTDC = new TDirectoryFile("dirTDC", "dirTDC"); // name, title

  // Directory for DeltaE-E plots
  dirDEEplots = new TDirectoryFile("DEEplots","DEEplots");
  dirPSD = dirDEEplots->mkdir("PSDplots","PSDplots");
  dirhitmaps  = new TDirectoryFile("hitmaps","hitmaps");

  // Directory for all correlations and inv-mass
  dirInvMass = new TDirectoryFile("InvMass","InvMass");
	singleFrags = dirInvMass->mkdir("singleFrags", "singleFrags");
  dir4He = dirInvMass->mkdir("4He","4He");
  dir5He = dirInvMass->mkdir("5He","5He");
  dir6He = dirInvMass->mkdir("6He","6He");
  dir5Li = dirInvMass->mkdir("5Li","5Li");
  dir6Li = dirInvMass->mkdir("6Li","6Li");
  dir7Li = dirInvMass->mkdir("7Li","7Li");
  dir6Be = dirInvMass->mkdir("6Be","6Be");
  dir7Be = dirInvMass->mkdir("7Be","7Be");
  dir8Be = dirInvMass->mkdir("8Be","8Be");
  dir9B  = dirInvMass->mkdir("9B","9B");
  dir8C  = dirInvMass->mkdir("8C","8C");
  dir10C = dirInvMass->mkdir("10C","10C");
	dir12C = dirInvMass->mkdir("12C","12C");
  dir9N  = dirInvMass->mkdir("9N","9N");
  dir12O = dirInvMass->mkdir("12O","12O");

  dirSummary->cd();

  int Nbin         = 5000;
  float Ecal_Emax  = 50.0;
  float Delta_Emax = 16.0;

  //// Create full summaries

  // Energies, Raw+Calibrated
  sumFrontE_R = new TH2I("sumFrontE_R","",4*hinpchans,0,4*hinpchans,1024,0,8192);
  sumFrontE_R->SetOption("colz");
  sumBackE_R = new TH2I("sumBackE_R","",4*hinpchans,0,4*hinpchans,1024,0,8192);
  sumBackE_R->SetOption("colz");
  sumDeltaE_R = new TH2I("sumDeltaE_R","",4*hinpchans,0,4*hinpchans,1024,0,8192);
  sumDeltaE_R->SetOption("colz");
  sumFrontE_cal = new TH2I("sumFrontE_cal","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  sumFrontE_cal->SetOption("colz");
  sumFrontE_addback = new TH2I("sumFrontE_addback","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  sumFrontE_addback->SetOption("colz");

  sumBackE_cal = new TH2I("sumBackE_cal","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  sumBackE_cal->SetOption("colz");
  sumBackE_addback = new TH2I("sumBackE_addback","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  sumBackE_addback->SetOption("colz");
  sumDeltaE_cal = new TH2I("sumDeltaE_cal","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  sumDeltaE_cal->SetOption("colz");
  sumDeltaE_addback = new TH2I("sumDeltaE_addback","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  sumDeltaE_addback->SetOption("colz");

  sumEtot_cal = new TH2I("sumEtot_cal","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  sumEtot_cal->SetOption("colz");
  AngleCorrSum_cal = new TH2I("AngleCorrSum_cal","",4*hinpchans,0,4*hinpchans,Nbin,0,Ecal_Emax);
  AngleCorrSum_cal->SetOption("colz");

  AngleCorrFrontE_cal = new TH2I("AngleCorrFrontE_cal","",4*hinpchans,0,4*hinpchans,Nbin/4,0,Ecal_Emax);
  AngleCorrFrontE_cal->SetOption("colz");
  AngleCorrDeltaE_cal = new TH2I("AngleCorrDeltaE_cal","",4*hinpchans,0,4*hinpchans,Nbin/4,0,Delta_Emax);
  AngleCorrDeltaE_cal->SetOption("colz");

  // Times
  sumFrontTime_R = new TH2I("sumFrontTime_R","",4*hinpchans,0,4*hinpchans,512,0,16383);
  sumFrontTime_R->SetOption("colz");
  sumFrontTime_cal = new TH2I("sumFrontTime_cal","",4*hinpchans,0,4*hinpchans,512,0,16383);
  sumFrontTime_cal->SetOption("colz");
  sumBackTime_R = new TH2I("sumBackTime_R","",4*hinpchans,0,4*hinpchans,512,0,16383);
  sumBackTime_R->SetOption("colz");
  sumBackTime_cal = new TH2I("sumBackTime_cal","",4*hinpchans,0,4*hinpchans,512,0,16383);
  sumBackTime_cal->SetOption("colz");
  sumDeltaTime_R = new TH2I("sumDeltaTime_R","",4*hinpchans,0,4*hinpchans,512,0,16383);
  sumDeltaTime_R->SetOption("colz");
  sumDeltaTime_cal = new TH2I("sumDeltaTime_cal","",4*hinpchans,0,4*hinpchans,512,0,16383);
  sumDeltaTime_cal->SetOption("colz");

  sumFrontTimeMult1_cal = new TH2I("sumFrontTimeMult1_cal","",4*hinpchans,0,4*hinpchans,512,0,16383);
  sumFrontTimeMult1_cal->SetOption("colz");
  
	//CsI multiplicity plots. A few variants
	CsI_mult_R = new TH1I("CsI_mult_R","CsI_mult_R",10,-0.5,9.5); //Raw CsI multiplicity "crystals in data stream"
	CsI_mult_AQT = new TH1I("CsI_mult_AQT","CsI_mult_AQT",10,-0.5,9.5);; //Same as raw but matched to QDC and TDC
	CsI_mult_M = new TH1I("CsI_mult_M","CsI_mult_M",10,-0.5,9.5);; //Now requiring FB matching
	CsI_mult_M_v_R = new TH2I("CsI_mult_M_v_R","CsI_mult_M_v_R",10,-0.5,9.5,10,-0.5,9.5);
	CsI_mult_PID = new TH1I("CsI_mult_PID","CsI_mult_PID",10,-0.5,9.5);; //Requires PID
	CsI_mult_PID_v_R = new TH2I("CsI_mult_PID_v_R","CsI_mult_PID_v_R",10,-0.5,9.5,10,-0.5,9.5);

  ostringstream name;
  for (int i=0;i<4;i++) {
  	name.str("");
  	name << "FrontvsBack_" << i;
  	FrontvsBack[i] = new TH2I(name.str().c_str(),"",500,0,80,500,0,80);
  	name.str("");
  	name << "FrontvsBack_sionly_" << i;
  	FrontvsBack_sionly[i] = new TH2I(name.str().c_str(),"",500,0,80,500,0,80);
  }


  // Create all 1d Front spectra
	// # of Gobbi telescopes is always 4
	FrontE_R.resize(4);
	FrontElow_R.resize(4);
	FrontTime_R.resize(4);
	FrontE_cal.resize(4);
	BackE_R.resize(4);
	BackElow_R.resize(4);
	BackTime_R.resize(4);
	BackE_cal.resize(4);
	AngleCorrE.resize(4);
	AngleCorr_noCorr.resize(4);
	AngleCorrE_R.resize(4);
	DeltaE_R.resize(4);
	DeltaElow_R.resize(4);
	DeltaTime_R.resize(4);
	DeltaE_cal.resize(4);
	AngleCorrDeltaE.resize(4);
	AngleCorrDeltaE_noCorr.resize(4);
	AngleCorrDeltaE_R.resize(4);
  for (size_t tele_i = 0; tele_i < 4; tele_i++) {
		FrontE_R[tele_i].resize(hinpchans);
		FrontElow_R[tele_i].resize(hinpchans);
		FrontTime_R[tele_i].resize(hinpchans);
		FrontE_cal[tele_i].resize(hinpchans);
		BackE_R[tele_i].resize(hinpchans);
		BackElow_R[tele_i].resize(hinpchans);
		BackTime_R[tele_i].resize(hinpchans);
		BackE_cal[tele_i].resize(hinpchans);
		AngleCorrE[tele_i].resize(hinpchans);
		AngleCorr_noCorr[tele_i].resize(hinpchans);
		AngleCorrE_R[tele_i].resize(hinpchans);
		DeltaE_R[tele_i].resize(hinpchans);
		DeltaElow_R[tele_i].resize(hinpchans);
		DeltaTime_R[tele_i].resize(hinpchans);
		DeltaE_cal[tele_i].resize(hinpchans);
		AngleCorrDeltaE[tele_i].resize(hinpchans);
		AngleCorrDeltaE_noCorr[tele_i].resize(hinpchans);
		AngleCorrDeltaE_R[tele_i].resize(hinpchans);

    for (size_t chan_i = 0; chan_i < hinpchans; chan_i++) {

      // Individual Front Energy
      dir1dFrontE_R->cd();
      name.str("");
      name << "FrontE_R" << tele_i << "_" << chan_i;
      FrontE_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",2048,0,8192);

      dir1dFrontlowE_R->cd();
      name.str("");
      name << "FrontElow_R" << tele_i << "_" << chan_i;
      FrontElow_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,4095);

      dir1dFrontTime_R->cd();
      name.str("");
      name << "FrontTime_R" << tele_i << "_" << chan_i;
      FrontTime_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,16383);

      dir1dFrontE_cal->cd();
      name.str("");
      name << "FrontE_cal" << tele_i << "_" << chan_i;
      FrontE_cal[tele_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin,5,Ecal_Emax);

      // Individual Back Energy
      dir1dBackE_R->cd();
      name.str("");
      name << "BackE_R" << tele_i << "_" << chan_i;
      BackE_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",2048,0,8192);

      dir1dBacklowE_R->cd();
      name.str("");
      name << "BackElow_R" << tele_i << "_" << chan_i;
      BackElow_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,4095);

      dir1dBackTime_R->cd();
      name.str("");
      name << "BackTime_R" << tele_i << "_" << chan_i;
      BackTime_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,16383);

      dir1dBackE_cal->cd();
      name.str("");
      name << "BackE_cal" << tele_i << "_" << chan_i;
      BackE_cal[tele_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin,0,Ecal_Emax);

      // Individual DeltaE
      dir1dDeltaE_R->cd();
      name.str("");
      name << "DeltaE_R" << tele_i << "_" << chan_i;
      DeltaE_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,4095);

      dir1dDeltalowE_R->cd();
      name.str("");
      name << "DeltaElow_R" << tele_i << "_" << chan_i;
      DeltaElow_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,4095);

      dir1dDeltaTime_R->cd();
      name.str("");
      name << "DeltaTime_R" << tele_i << "_" << chan_i;
      DeltaTime_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,16383);

      dir1dDeltaE_cal->cd();
      name.str("");
      name << "DeltaE_cal" << tele_i << "_" << chan_i;
      DeltaE_cal[tele_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin,0,Delta_Emax);

      dirAngleCorrFrontE->cd();
      name.str("");
      name << "AngleCorrFrontE" << tele_i << "_" << chan_i;
      AngleCorrE[tele_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin/3,5,Ecal_Emax);

      dirAngleCorrDeltaE->cd();
      name.str("");
      name << "AngleCorrDeltaE" << tele_i << "_" << chan_i;
      AngleCorrDeltaE[tele_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin/3,0,Delta_Emax);

			// Angle corrected energies I think
			dirAngleCorrFrontE->cd();
      name.str("");
      name << "AngleCorrE_noCorr" << tele_i << "_" << chan_i;
      AngleCorr_noCorr[tele_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin/3,5,Ecal_Emax);

      dirAngleCorrDeltaE->cd();
      name.str("");
      name << "AngleCorrDeltaE_noCorr" << tele_i << "_" << chan_i;
      AngleCorrDeltaE_noCorr[tele_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin/3,0,Delta_Emax);

			dirAngleCorrFrontE->cd();
      name.str("");
      name << "AngleCorrE_R" << tele_i << "_" << chan_i;
      AngleCorrE_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",2048,0,8192);

      dirAngleCorrDeltaE->cd();
      name.str("");
      name << "AngleCorrDeltaE_R" << tele_i << "_" << chan_i;
      AngleCorrDeltaE_R[tele_i][chan_i] = new TH1I(name.str().c_str(),"",2048,0,8192);
    }
  }
	
	//TDC plots
	dirTDC->cd();
	for (int i=0;i<16;i++) {
		name.str("");
		name << "TDCspect_" << i;
		TDC_Plot[i] = new TH1I(name.str().c_str(),"",1000,-500,500);

		if (i > 3) {
			name.str("");
			name << "TDCspect_TN_shift" << i-4;
			TDC_Plot_TN_shift[i-4] = new TH1I(name.str().c_str(),"",1000,-500,500);
		}
	}
	
	TDC_sum = new TH2I("TDC_sum","",16,-0.5,15.5,1000,-500,500);
	TDC_sum_TN = new TH2I("TDC_sum_TN","",12,-0.5,11.5,1000,-500,500);
	TDC_sum_TN_shift = new TH2I("TDC_sum_TN_shift","",12,-0.5,11.5,1000,-500,500);
  // Create all spectra based on quadrants
  
  
  dirDEEplots->cd();
  for (int quad = 0; quad < 4; quad++) {
    name.str("");
    name << "DEE_simple" << quad;
    DEE_simple[quad] = new TH2I(name.str().c_str(),"",500,0,80,500,0,16); //E is x, DE is y

    name.str("");
    name << "frontdeltastripnum_" << quad;
    frontdeltastripnum[quad] = new TH2I(name.str().c_str(),"",32,-0.5,31.5,32,-0.5,31.5); //E is x, DE is y

    name.str("");
    name << "DEE" << quad;
    DEE[quad] = new TH2I(name.str().c_str(),"",500,0,80,800,0,22); //E is x, DE is y

    name.str("");
    name << "timediff" << quad;   
    timediff[quad] = new TH1I(name.str().c_str(),"",1000,-2000,2000);
  }
  
  //Create DEE CsI plots
  //Sorry Henry this will be hard coded for now
  //Apology not accepted Johnathan (JK, I fixed the problems. CsI histograms are once again vectorized and initialized inside `Gobbi28.cpp`)

  dirhitmaps->cd();
  xyhitmap_allE = new TH2I("xyhitmap_allE","", 100,-10,10,100,-10,10);
  xyhitmap = new TH2I("xyhitmap","", 100,-10,10,100,-10,10);
  tphitmap = new TH2I("tphitmap","",1000,-40,40,1000,-40,40);
  xyhitmap_sionly = new TH2I("xyhitmap_sionly","", 100,-10,10,100,-10,10);
  xyhitmap_tgate_orA = new TH2I("xyhitmap_tgate_orA","", 100,-10,10,100,-10,10);
  xyhitmap_multiCsI = new TH2I("xyhitmap_multiCsI","", 100,-10,10,100,-10,10);
  protonhitmap = new TH2I("protonhitmap","", 100,-10,10,100,-10,10);
  deuteronhitmap = new TH2I("deuteronhitmap","", 100,-10,10,100,-10,10);
  tritonhitmap = new TH2I("tritonhitmap","", 100,-10,10,100,-10,10);
  alphahitmap = new TH2I("alphahitmap","", 100,-10,10,100,-10,10);
  He6hitmap = new TH2I("He6hitmap","", 100,-10,10,100,-10,10);
  Lihitmap = new TH2I("Lihitmap","", 100,-10,10,100,-10,10);
  LiVETOhitmap = new TH2I("LiVETOhitmap","", 100,-10,10,100,-10,10);

  hitmapof_p = new TH2I("hitmapof_p","", 100,-10,10,100,-10,10);
  hitmapof_6He = new TH2I("hitmapof_6He","", 100,-10,10,100,-10,10);

  Evstheta[0] = new TH2I("Evstheta0","",50,0,25,Nbin,0,Ecal_Emax);
  Evstheta[1] = new TH2I("Evstheta1","",50,0,25,Nbin,0,Ecal_Emax);
  Evstheta[2] = new TH2I("Evstheta2","",50,0,25,Nbin,0,Ecal_Emax);
  Evstheta[3] = new TH2I("Evstheta3","",50,0,25,Nbin,0,Ecal_Emax);
  Evstheta_all = new TH2I("Evstheta_all","",50,0,25,Nbin,0,Ecal_Emax);
  Theta = new TH1I("Theta","",50,0,25);

  ProtonEnergy = new TH2I("ProtonEnergy","",80,0,80,50,0,25);
  He3Energy = new TH2I("He3Energy","",200,0,200,50,0,25);
  AlphaEnergy = new TH2I("AlphaEnergy","",250,0,250,50,0,25);

  dTime_proton = new TH1I("dTime_proton","",1500,-4000,2000);
  dTime_deuteron = new TH1I("dTime_deuteron","",1500,-4000,2000);
  dTime_triton = new TH1I("dTime_triton","",1500,-4000,2000);
  dTime_alpha = new TH1I("dTime_alpha","",1500,-4000,2000);
  dTime_He6 = new TH1I("dTime_He6","",1500,-4000,2000);
  dTime_Li = new TH1I("dTime_Li","",1500,-4000,2000);

  CorrelationTable = new TH2I("CorrelationTable","",9,-0.5,8.5,9,0,8.5);

	NeighborStripECorrelations = new TH2I("NeighborStripECorrelations","",500,0,20,500,0,20); // Order is energy ordered, so x is larger energy and y is smaller energy


  // He4
  dir4He->cd();
  Erel_4He_pt = new TH1I("Erel_4He_pt","",800,0,30);
  Ex_4He_pt = new TH1I("Ex_4He_pt","",350,19,26);
  ThetaCM_4He_pt = new TH1I("ThetaCM_4He_pt","",200,0,10);
  VCM_4He_pt = new TH1I("VCM_4He_pt","",100,0,14);

  He4_p_hitmap = new TH2I("He4_p_hitmap","", 100,-10,10,100,-10,10);
  He4_t_hitmap = new TH2I("He4_t_hitmap","", 100,-10,10,100,-10,10);
  DEE_He4[0] = new TH2I("DEE_He4_quad0","",500,0,80,800,0,22); //E is x, DE is y
  DEE_He4[1] = new TH2I("DEE_He4_quad1","",500,0,80,800,0,22); //E is x, DE is y
  DEE_He4[2] = new TH2I("DEE_He4_quad2","",500,0,80,800,0,22); //E is x, DE is y
  DEE_He4[3] = new TH2I("DEE_He4_quad3","",500,0,80,800,0,22); //E is x, DE is y
  Erel_pt_costhetaH = new TH2I("Erel_pt_costhetaH","",200,0,8,25,-1,1);

  Erel_4He_dd = new TH1I("Erel_4He_dd","",800,0,30);
  Ex_4He_dd = new TH1I("Ex_4He_dd","",350,19,26);
  ThetaCM_4He_dd = new TH1I("ThetaCM_4He_dd","",200,0,10);
  VCM_4He_dd = new TH1I("VCM_4He_dd","",100,0,14);
  Erel_dd_costhetaH = new TH2I("Erel_dd_costhetaH","",50,0,2,25,-1,1);

  // He5
  dir5He->cd();
  Erel_5He_dt = new TH1I("Erel_5He_dt","",800,0,30);
  Ex_5He_dt = new TH1I("Ex_5He_dt","",800,-2,30);
  ThetaCM_5He_dt = new TH1I("ThetaCM_5He_dt","",200,0,10);
  VCM_5He_dt = new TH1I("VCM_5He_dt","",100,0,14);
	Erel_dt_costhetaH = new TH2I("Erel_dt_costhetaH","",750,0,30,25,-1,1);

  // He6
  dir6He->cd();
  Erel_6He_tt = new TH1I("Erel_6He_tt","",800,0,30);
  Ex_6He_tt = new TH1I("Ex_6He_tt","",800,-2,30);
  ThetaCM_6He_tt = new TH1I("ThetaCM_6He_tt","",200,0,10);
  VCM_6He_tt = new TH1I("VCM_6He_tt","",100,0,14);
	Erel_tt_costhetaH = new TH2I("Erel_tt_costhetaH","",750,0,30,25,-1,1);

  // Li5
  dir5Li->cd();
  Erel_5Li_pa = new TH1I("Erel_5Li_pa","",800,0,30);
  Ex_5Li_pa = new TH1I("Ex_5Li_pa","",800,-2,30);
  ThetaCM_5Li_pa = new TH1I("ThetaCM_5Li_pa","",200,0,10);
  VCM_5Li_pa = new TH1I("VCM_5Li_pa","",100,0,14);
	Erel_pa_costhetaH = new TH2I("Erel_pa_costhetaH","",750,0,30,25,-1,1);

  Erel_5Li_d3He = new TH1I("Erel_5Li_d3He","",800,0,30);
  Ex_5Li_d3He = new TH1I("Ex_5Li_d3He","",800,-2,30);
  ThetaCM_5Li_d3He = new TH1I("ThetaCM_5Li_d3He","",200,0,10);
  VCM_5Li_d3He = new TH1I("VCM_5Li_d3He","",100,0,14);
	Erel_d3He_costhetaH = new TH2I("Erel_d3He_costhetaH","",750,0,30,25,-1,1);

  // Li6
  dir6Li->cd();
  // -> p + n + a
  Erel_6Li_npa = new TH1I("Erel_6Li_npa","",500,0,5);
  Ex_6Li_npa_trans = new TH1I("Ex_6Li_npa_trans","",500,0,5);
  Ex_6Li_npa_long = new TH1I("Ex_6Li_npa_long","",500,0,5);
  Ex_6Li_npa = new TH1I("Ex_6Li_npa","",500,0,5);
  cos_thetaH_npa = new TH1I("cos_thetaH_npa","",100,-1.1,1.1);
  ThetaCM_6Li_npa = new TH1I("ThetaCM_6Li_npa","",200,0,25);
  VCM_6Li_npa = new TH1I("VCM_6Li_npa","",100,1.5,4.5);

  cos_npa_thetaH = new TH1I("cos_npa_thetaH","",100,-1.1,1.1);
  Erel_npa_cosThetaH = new TH2I("Erel_npa_cosThetaH","",200,0,3,25,-1,1);
  
  Erel_6Li_da = new TH1I("Erel_6Li_da","",2000,0,15);
  Erel_6Li_da_tgate_orA = new TH1I("Erel_6Li_da_tgate_orA","",2000,0,15);
  Ex_6Li_da_trans = new TH1I("Ex_6Li_da_trans","",2000,-5,10);
  Ex_6Li_da_long = new TH1I("Ex_6Li_da_long","",2000,-5,10);
  Ex_6Li_da = new TH1I("Ex_6Li_da","",2000,-5,10);
  cos_thetaH_da = new TH1I("cos_thetaH_da","",100,-1.1,1.1);
  ThetaCM_6Li_da = new TH1I("ThetaCM_6Li_da","",200,0,25);
  VCM_6Li_da = new TH1I("VCM_6Li_da","",100,1.5,4.5);
  VCM_vs_ThetaCM = new TH2I("VCM_vs_ThetaCM","", 200,0,25,100,1.5,4.5);

  cos_da_thetaH = new TH1I("cos_da_thetaH","",100,-1.1,1.1);
  Erel_da_cosThetaH = new TH2I("Erel_da_cosThetaH","",200,0,3,25,-1,1);

  deutE_gate = new TH1I("deutE_gate","",200,0,35);
  alphaE_gated = new TH1I("alphaE_gated","",200,0,35);
  deutE_gate_cosThetaH = new TH2I("deutE_gate_cosThetaH","",200,0,35,100,-1,1);
  alphaE_gate_cosThetaH = new TH2I("alphaE_gate_cosThetaH","",200,0,35,100,-1,1);

	react_origin_tdiff = new TH1I("react_origin_tdiff","",2000,-100,100);
	
	xyhitmap_6Li_plus = new TH2I("xyhitmap_6Li_plus","",100,-10,10,100,-10,10);

  // Li7
	// p + 6He
  dir7Li->cd();
  Erel_7Li_p6He = new TH1I("Erel_7Li_p6He","",200,0,8);
  Erel_7Li_p6He_Q = new TH2I("Erel_7Li_p6He_Q","",200,0,8,200,-5,15);
  Erel_7Li_cosThetaH = new TH2I("Erel_7Li_cosThetaH","",200,0,8,25,-1,1);
  Erel_7Li_p6He_lowres = new TH1I("Erel_7Li_p6He_lowres","",100,0,8);
  Ex_7Li_p6He_transverse = new TH1I("Ex_7Li_p6He_transverse","",400,10,18);
  Ex_7Li_p6He_transverse2 = new TH1I("Ex_7Li_p6He_transverse2","",400,10,18);
  Erel_7Li_p6He_pFor = new TH1I("Erel_7Li_p6He_pFor","",200,0,17);
  Ex_7Li_p6He_timegate = new TH1I("Ex_7Li_p6He_timegate","",400,10,18);

  cos_thetaH = new TH1I("cos_thetaH","",100,-1.1,1.1);
  cos_thetaH_lowErel = new TH1I("cos_thetaH_lowErel","",100,-1.1,1.1);
  missingmass = new TH1I("missingmass","",500,-150,100);
  Erel_missingmass = new TH2I("Erel_missingmass","",200,0,8,200,-10,10);
  Qvalue = new TH1I("Qvalue","",500,-5,15);
  Qvalue2 = new TH1I("Qvalue2","",500,-5,15);
  Ex_7Li_p6He_clean = new TH1I("Ex_7Li_p6He_clean","",400,10,18);

  Ex_7Li_p6He = new TH1I("Ex_7Li_p6He","",400,10,18);
  ThetaCM_7Li_p6He = new TH1I("ThetaCM_7Li_p6He","",200,0,20);
  VCM_7Li_p6He = new TH1I("VCM_7Li_p6He","",100,1.5,4);
  VCM_7Li_p6He_lowErel = new TH1I("VCM_7Li_p6He_lowErel","",100,1.5,4);
  //Vlab_LF_p6He = new TH1I("Vlab_LF_p6He","",100,0,1); //units of v/c
  //Vlab_HF_p6He = new TH1I("Vlab_HF_p6He","",100,0,1); //units of v/c
  //cosbeamCMtoHF_Ex_p6He = new TH2I("cosbeamCMtoHF_Ex_p6He","",100,-2,15,100,-1,1); //2d

  dTime_7Li_proton = new TH1I("dTime_7Li_proton","",1500,-4000,2000);
  dTime_7Li_He6 = new TH1I("dTime_7Li_He6","",1500,-4000,2000);

  Ex_7Li_p6He_ExvsEp = new TH2I("Ex_7Li_p6He_ExvsEp","", 400,10,18, 100,0,16);

	// t + alpha
  Erel_7Li_ta = new TH1I("Erel_7Li_ta","",400,0,8);
  Ex_7Li_ta = new TH1I("Ex_7Li_ta","",400,2,10);
  Ex_7Li_ta_trans = new TH1I("Ex_7Li_ta_trans","",400,2,10);
  Ex_7Li_ta_long = new TH1I("Ex_7Li_ta_long","",400,2,10);
  Ex_7Li_ta_bad = new TH1I("Ex_7Li_ta_bad","",400,2,10);
  ThetaCM_7Li_ta = new TH1I("ThetaCM_7Li_ta","",200,0,25);
  VCM_7Li_ta = new TH1I("VCM_7Li_ta","",100,1.5,4.5);
  cos_ta_thetaH = new TH1I("cos_ta_thetaH","",100,-1.1,1.1);
  Erel_ta_cosThetaH = new TH2I("Erel_ta_cosThetaH","",200,0,8,25,-1,1);
  Ex_tar = new TH1I("Ex_tar","",400,-20,20);
  Erel_vs_Extar = new TH2I("Erel_vs_Extar","",100,0,8,100,-20,20);

  Ex_7Li_ta_timegate = new TH1I("Ex_7Li_ta_timegate","",400,2,10);
  hitmapcheck1 = new TH2I("hitmapcheck1","", 100,-10,10,100,-10,10);
  hitmapcheck2 = new TH2I("hitmapcheck2","", 100,-10,10,100,-10,10);
  dTime_7Li_triton = new TH1I("dTime_7Li_triton","",1500,-4000,2000);
  dTime_7Li_alpha = new TH1I("dTime_7Li_alpha","",1500,-4000,2000);

  seperate_quad_Ex_7Li_ta = new TH1I("seperate_quad_Ex_7Li_ta","",400,2,10);

  DEE_shoulderevents = new TH2I("DEE_shoulderevents","",500,0,80,800,0,22);

  // Be6
  dir6Be->cd();
  Erel_6Be_2pa = new TH1I("Erel_6Be_2pa","",800,0,30);
  ThetaCM_6Be_2pa = new TH1I("ThetaCM_6Be_2pa","",200,0,10);
  VCM_6Be_2pa = new TH1I("VCM_6Be_2pa","",100,0,14);

  // Be7
  dir7Be->cd();
  Erel_7Be_a3He = new TH1I("Erel_7Be_a3He","",800,0,30);
  Ex_7Be_a3He = new TH1I("Ex_7Be_a3He","",800,-5,30);
  ThetaCM_7Be_a3He = new TH1I("ThetaCM_7Be_a3He","",200,0,10);
  VCM_7Be_a3He = new TH1I("VCM_7Be_a3He","",100,0,14);

  Erel_7Be_p6Li = new TH1I("Erel_7Be_p6Li","",800,0,30);
  Ex_7Be_p6Li = new TH1I("Ex_7Be_p6Li","",800,-5,30);
  ThetaCM_7Be_p6Li = new TH1I("ThetaCM_7Be_p6Li","",200,0,10);
  VCM_7Be_p6Li = new TH1I("VCM_7Be_p6Li","",100,0,14);

  // Be8
  dir8Be->cd();
  Erel_8Be_aa = new TH1I("Erel_8Be_aa","",1000,0,20);
  Ex_8Be_aa = new TH1I("Ex_8Be_aa","",1600,-1,7);
  ThetaCM_8Be_aa = new TH1I("ThetaCM_8Be_aa","",200,0,25);
  VCM_8Be_aa = new TH1I("VCM_8Be_aa","",100,0,14);
  Erel_aa_cosThetaH = new TH2I("Erel_aa_cosThetaH","",100,0,8,25,-1,1);

  Erel_8Be_p7Li = new TH1I("Erel_8Be_p7Li","",800,0,17);
  Ex_8Be_p7Li = new TH1I("Ex_8Be_p7Li","",400,17,25);
  Ex_8Be_p7Li_trans = new TH1I("Ex_8Be_p7Li_trans","",400,17,25);
  ThetaCM_8Be_p7Li = new TH1I("ThetaCM_8Be_p7Li","",100,0,15);
  VCM_8Be_p7Li = new TH1I("VCM_8Be_p7Li","",50,2,4);
  cos_p7Li_thetaH = new TH1I("cos_p7Li_thetaH","",100,-1.1,1.1);
  Erel_p7Li_cosThetaH = new TH2I("Erel_p7Li_cosThetaH","",200,0,8,25,-1,1);

  ProtonEnergies_p7Li = new TH1I("ProtonEnergies_p7Li","",200,0,40);
  LithiumEnergies_p7Li = new TH1I("LithiumEnergies_p7Li","",200,0,40);

  dTime_8Be_proton = new TH1I("dTime_8Be_proton","",1500,-4000,2000);
  dTime_8Be_Li7 = new TH1I("dTime_8Be_Li7","",1500,-4000,2000);

  Ex_8Be_p7Li_timegate = new TH1I("Ex_8Be_p7Li_timegate","",400,17,25);

  Erel_8Be_pta = new TH1I("Erel_8Be_pta","",800,0,30);
  Ex_8Be_pta = new TH1I("Ex_8Be_pta","",250,20,25);
  Ex_8Be_pta_trans = new TH1I("Ex_8Be_pta_trans","",250,20,25);
  ThetaCM_8Be_pta = new TH1I("ThetaCM_8Be_pta","",200,0,10);
  VCM_8Be_pta = new TH1I("VCM_8Be_pta","",100,0,14);
  cos_pta_thetaH = new TH1I("cos_pta_thetaH","",100,-1.1,1.1);
  Erel_pta_cosThetaH = new TH2I("Erel_pta_cosThetaH","",200,0,8,25,-1,1);

  Erel_7Li_ta_fake = new TH1I("Erel_7Li_ta_fake","",400,0,8);
  Ex_7Li_ta_fake = new TH1I("Ex_7Li_ta_fake","",400,2,10);

  Ex_8Be_7LiGate = new TH1I("Ex_8Be_pta_7LiGate","",250,20,25);

  // B9
  dir9B->cd();
  Erel_9B_paa = new TH1I("Erel_9B_paa","",800,0,17);
  Ex_9B_paa = new TH1I("Ex_9B_paa","",800,-2,15);
  Ex_9B_p8Be = new TH1I("Ex_9B_p8Be","",800,-2,15);
  Ex_9B_aa = new TH1I("Ex_8Be_in_9B_aa","",800,-2,15);
  ThetaCM_9B_paa = new TH1I("ThetaCM_9B_paa","",200,0,10);
  VCM_9B_paa = new TH1I("VCM_9B_paa","",100,0,14);

	// C8
	dir8C->cd();
	Erel_8C_4pa = new TH1I("Erel_8C_4pa","",800,0,17);
	Ex_8C_4pa = new TH1I("Ex_8C_4pa","",800,-2,15);
	ThetaCM_8C_4pa = new TH1I("ThetaCM_8C_4pa","",200,0,10);
	VCM_8C_4pa = new TH1I("VCM_8C_4pa","",100,0,14);
	
	// C10
	dir10C->cd();
	Erel_10C_2p2a = new TH1I("Erel_10C_2p2a","",800,0,17);
	Ex_10C_2p2a = new TH1I("Ex_10C_2p2a","",800,-2,15);
	ThetaCM_10C_2p2a = new TH1I("ThetaCM_10C_2p2a","",200,0,10);
	VCM_10C_2p2a = new TH1I("VCM_10C_2p2a","",100,0,14);

	// C12
	dir12C->cd();
	Erel_12C_3a = new TH1I("Erel_12C_3a","",800,0,17);
	Ex_12C_3a = new TH1I("Ex_12C_3a","",800,-2,15);
	ThetaCM_12C_3a = new TH1I("ThetaCM_12C_3a","",200,0,10);
	VCM_12C_3a = new TH1I("VCM_12C_3a","",100,0,14);

	// N9
	dir9N->cd();
	Erel_9N_5pa = new TH1I("Erel_9N_5pa","",800,0,17);
	ThetaCM_9N_5pa = new TH1I("ThetaCM_9N_5pa","",200,0,10);
	VCM_9N_5pa = new TH1I("VCM_9N_5pa","",100,0,14);
	
	// O12
	dir12O->cd();
	Erel_12O_4p2a = new TH1I("Erel_12O_4p2a","",800,0,17);
	Ex_12O_4p2a = new TH1I("Ex_12O_4p2a","",800,-2,15);
	ThetaCM_12O_4p2a = new TH1I("ThetaCM_12O_4p2a","",200,0,10);
	VCM_12O_4p2a = new TH1I("VCM_12O_4p2a","",100,0,14);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

histo::~histo() {
  file_read->Write();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void histo::reset() {
	//solutions.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Eventually, this function is where global tree should be filled with pre-solution Gobbi information
void histo::Fill() {
	//tpar->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



