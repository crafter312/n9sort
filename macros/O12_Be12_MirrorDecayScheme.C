{
 gROOT->Reset();
  TStyle * Sty = new TStyle("MyStyle","MyStyle");
  Sty->SetOptTitle(0);
  Sty->SetOptStat(0);
  Sty->SetLineWidth(4);
  Sty->SetPalette(55);
  Sty->SetCanvasColor(10);
  Sty->SetCanvasBorderMode(0);
  Sty->SetFrameLineWidth(0);
  Sty->SetFrameFillColor(10);
  Sty->SetPadColor(10);
  Sty->SetPadTickX(1);
  Sty->SetPadTickY(1);
  Sty->SetPadBottomMargin(.03);
  Sty->SetPadTopMargin(.03);
  Sty->SetPadLeftMargin(.11);
  Sty->SetPadRightMargin(.05);
  Sty->SetHistLineWidth(3);
  Sty->SetFuncWidth(3);
  Sty->SetFuncColor(kGreen);
  Sty->SetLineWidth(3);
  Sty->SetLabelSize(0.06,"xyz");
  Sty->SetLabelOffset(0.005,"y");
  Sty->SetLabelOffset(0.02,"x");
  Sty->SetLabelColor(kBlack,"xyz");
  Sty->SetTitleSize(0.06,"y");
  Sty->SetTitleSize(0.07,"x");
  Sty->SetTitleOffset(0.95,"y");
  Sty->SetTitleOffset(1.2,"x");
  Sty->SetTitleFillColor(10);
  Sty->SetTitleTextColor(kBlack);
  Sty->SetTickLength(.05,"xz");
  Sty->SetTickLength(.025,"y");
  Sty->SetNdivisions(110,"y");
  Sty->SetNdivisions(210,"x");
  Sty->SetEndErrorSize(0);
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();

  TCanvas *c1 = new TCanvas("c1","",1400,768);

  TLatex tt;
  tt.SetTextSize(.06);

  TArrow arrow;
  arrow.SetFillColor(1);
  arrow.SetAngle(35);

  tt.SetTextSize(0.045);

  Sty->SetNdivisions(18,"y");
  Sty->SetTitleOffset(0.7,"y");
  Sty->SetTitleSize(0.07,"y");

  TH2S frame2("frame","",10,0,2,10,-2,16);
  frame2.GetXaxis()->SetLabelSize(0.);
  frame2.GetXaxis()->SetLabelColor(0.);
  frame2.GetXaxis()->SetTickSize(0.);
  frame2.GetXaxis()->SetAxisColor(0);
  //frame2.GetXaxis()->SetAxisStyle(0);

  frame2.GetXaxis()->SetLabelSize(0);
  frame2.GetXaxis()->SetTickLength(0);
  frame2.SetMarkerStyle(20);
  frame2.GetYaxis()->SetTitle("Energy (MeV)");
  frame2.GetYaxis()->CenterTitle();
  frame2.Draw();

  //set energy level lines
  TLine ll;
  ll.SetLineWidth(4);
  double indentX = 0.15;
  double LevelXwide = 0.2;
  double gap = 0.09;

  double Mg21_gamma = 1.598;

  //Get masses
	const double mass_n = 8.07132;
  const double mass_p = 7.28897;
  const double mass_OPha = 2.42491;
	const double mass_12Be = 25.0778; // AME2020
	const double mass_12O = 32.013; // AME 2020

  //Get ground states
	double Be12_gs = -(mass_12Be - (4.*mass_n) - (2.*mass_OPha));
	double O12_gs  = -(mass_12O - (4.*mass_p) - (2.*mass_OPha));

  cout << "Calculated ground states" << endl;
  cout << "4n + 2a: " << Be12_gs << endl;
  cout << "4p + 2a: " << O12_gs << endl;

  //Read in level files from /levels/

  //Read in 12Be states 
  string filename_exp = "Be12_exp.dat";
  ifstream file_exp(filename_exp.c_str());
  if (file_exp.is_open() != 1)
  {
    cout << " could not open file " << filename_exp << endl;
    abort();
  }

  //Level files have a 3 line header, skip all of this
  for (int i=0;i<3;i++) file_exp.ignore(200,'\n');

  string readvals[7];
  vector<string>  BeJ;
  vector<int> BeP;
  vector<int> BeN;
  vector<double> BeE;
  vector<double> BeW;
  vector<string> BeDM;
  vector<int> BeK; //ENSDF or our experiment? changes level coloring

  //Read in values
  string fline;
  while (getline(file_exp,fline)) {
    stringstream ss(fline);
    ss >> readvals[0] >> readvals[1] >> readvals[2] >> readvals[3] >> readvals[4] >> readvals[5] >> readvals[6];

    BeJ.push_back(readvals[0]);
    BeP.push_back(stoi(readvals[1]));
    BeN.push_back(stoi(readvals[2]));
    BeE.push_back(stod(readvals[3])/1000.);
    BeW.push_back(stod(readvals[4])/1000.);
    BeDM.push_back(readvals[5]);
    BeK.push_back(stoi(readvals[6]));
  }

  file_exp.close();
  file_exp.clear();

  //Read in 12O states 
  filename_exp = "O12_exp.dat";
  file_exp.open(filename_exp.c_str());
  if (file_exp.is_open() != 1)
  {
    cout << " could not open file " << filename_exp << endl;
    abort();
  }

  //Level files have a 3 line header, skip all of this
  for (int i=0;i<3;i++) file_exp.ignore(200,'\n');

  vector<string> OJ;
  vector<int> OP;
  vector<int> ON;
  vector<double> OE;
  vector<double> OW;
  vector<string> ODM;
  vector<int> OK; //ENSDF or our experiment? changes level coloring

  //Read in values
  while (getline(file_exp,fline)) {
    stringstream ss(fline);
    ss >> readvals[0] >> readvals[1] >> readvals[2] >> readvals[3] >> readvals[4] >> readvals[5] >> readvals[6];

    OJ.push_back(readvals[0]);
    OP.push_back(stoi(readvals[1]));
    ON.push_back(stoi(readvals[2]));
    OE.push_back(stod(readvals[3])/1000.);
    OW.push_back(stod(readvals[4])/1000.);
    ODM.push_back(readvals[5]);
    OK.push_back(stoi(readvals[6]));
  }

  file_exp.close();
  file_exp.clear();






  //Plot experimental levels and widths

	double BeThreshLstart = indentX;
  double BeThreshLstop = indentX+LevelXwide;
	double BeThreshHe8a = 8.957; // from TUNL
	{
		int levelcol = 1;
		ll.SetLineColor(levelcol);
		ll.DrawLine(BeThreshLstart,BeThreshHe8a,BeThreshLstop,BeThreshHe8a);
	}
	

  //Be12
  double Bestart = indentX+LevelXwide+gap;
  double Bestop = indentX+2.*LevelXwide+gap;
  for (int i=0;i< BeJ.size();i++) {

    int levelcol;
    if (BeK[i] == 1) levelcol = 1;
    else levelcol = 4;   

    ll.SetLineColor(levelcol);  

    ll.DrawLine(Bestart,BeE[i],Bestop,BeE[i]);

    TBox *expbox = new TBox(Bestart,BeE[i]-BeW[i]/2.,Bestop,BeE[i]+BeW[i]/2.);
    expbox->SetBit(kCanDelete);
    expbox->SetFillColor(levelcol);
    expbox->SetFillStyle(3004);
    expbox->Draw("same"); 
  }

	double BeThreshRstart = indentX+2.*LevelXwide+2.*gap;
  double BeThreshRstop = indentX+3.*LevelXwide+2.*gap;
	double BeThreshBe102n = 3.672; // from TUNL
	double BeThreshHe6He6 = 10.106; // from TUNL
	{
		int levelcol = 1;
		ll.SetLineColor(levelcol);
		ll.DrawLine(BeThreshRstart,BeThreshBe102n,BeThreshRstop,BeThreshBe102n);
		ll.DrawLine(BeThreshRstart,BeThreshHe6He6,BeThreshRstop,BeThreshHe6He6);
	}

	double OThreshLstart = indentX+3.*LevelXwide+3.*gap;
  double OThreshLstop = indentX+4.*LevelXwide+3.*gap;
	double OThreshC102p = -1.638; // from TUNL
	double OThreshC8a = 5.476; // from TUNL
	{
		int levelcol = 1;
		ll.SetLineColor(levelcol);
		ll.DrawLine(OThreshLstart,OThreshC102p,OThreshLstop,OThreshC102p);
		ll.DrawLine(OThreshLstart,OThreshC8a,OThreshLstop,OThreshC8a);
	}

  //O12
  double Ostart = indentX+4.*LevelXwide+4.*gap;
  double Ostop = indentX+5.*LevelXwide+4.*gap;
  for (int i=0;i<OJ.size();i++) {

    int levelcol;
    if (OK[i] == 1) levelcol = 1;
    else levelcol = 4;   

    ll.SetLineColor(levelcol);  

    ll.DrawLine(Ostart,OE[i],Ostop,OE[i]);

    TBox *expbox = new TBox(Ostart,OE[i]-OW[i]/2.,Ostop,OE[i]+OW[i]/2.);
    expbox->SetBit(kCanDelete);
    expbox->SetFillColor(levelcol);
    expbox->SetFillStyle(3004);
    expbox->Draw("same"); 
  }

	double OThreshRstart = indentX+5.*LevelXwide+5.*gap;
  double OThreshRstop = indentX+6.*LevelXwide+5.*gap;
	double OThresh6Be6Be = 4.737; // from TUNL
	{
		int levelcol = 1;
		ll.SetLineColor(levelcol);
		ll.DrawLine(OThreshRstart,OThresh6Be6Be,OThreshRstop,OThresh6Be6Be);
	}

 


  tt.SetTextSize(0.06);
  Float_t textsize = tt.GetTextSize();

  //labels
  tt.DrawLatex(Bestop - 0.16 ,BeE[0]-0.9,"^{12}Be");
  tt.DrawLatex(Ostop - 0.2025 ,OE[0]-0.9,"^{12}O");

  tt.SetTextSize(0.045);
  //Draw state J pi labels and threshold labels
  //Be12
	tt.DrawLatex(BeThreshLstop-0.2,BeThreshHe8a+0.1,"^{8}He+#alpha");
  for (int i=0;i< BeJ.size();i++) {

    //if no spin ("UN"), skip
    if ( BeJ[i] == "UN") continue;

    string tempJ =  BeJ[i];
    if (BeP[i] == 1) tempJ += "^{+}";
    else tempJ += "^{-}";

    tt.DrawLatex(Bestop-0.08,BeE[i]+0.1,tempJ.c_str());
  }
	tt.DrawLatex(BeThreshRstop-0.2,BeThreshBe102n+0.1,"^{10}Be+2n");
	tt.DrawLatex(BeThreshRstop-0.2,BeThreshHe6He6+0.1,"^{6}He+^{6}He");


  //O12
	tt.DrawLatex(OThreshLstop-0.2,OThreshC102p+0.1,"^{10}C+2p");
	tt.DrawLatex(OThreshLstop-0.2,OThreshC8a+0.1,"^{8}C+#alpha");
  for (int i=0;i<OJ.size();i++) {

    //if no spin ("UN"), skip
    if (OJ[i] == "UN") continue;

    string tempJ = OJ[i];
    if (OP[i] == 1) tempJ += "^{+}";
    else tempJ += "^{-}"; 

    tt.DrawLatex(Ostop-0.08,OE[i]+0.1,tempJ.c_str());
  }
	tt.DrawLatex(OThreshRstop-0.2,OThresh6Be6Be+0.1,"^{6}Be+^{6}Be");

/*
  arrow.SetAngle(30);
  arrow.SetLineWidth(3);
  arrow.SetLineColor(1);
  arrow.SetFillColor(1); 

  arrow.SetLineColor(kRed);
  arrow.SetFillColor(kRed); 

  //Si23 1p decays to Al22 g.s.
  float multiarrow_sub = 0.04;
  arrow.SetAngle(45);
  arrow.Drawarrow(Bestop,BeE[1],Ostart-multiarrow_sub,Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[2],Ostart-multiarrow_sub,Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[3],Ostart-multiarrow_sub,Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[4],Ostart,Al22_gs+0.08,.02);
  arrow.Drawarrow(Bestop,BeE[5],Ostart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.Drawarrow(Bestop,BeE[6],Ostart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.Drawarrow(Bestop,BeE[7],Ostart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.Drawarrow(Bestop,BeE[8],Ostart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.SetAngle(30);

  //Si23 2p decays to Al22 g.s.
  //arrow.SetLineStyle(2);
  arrow.SetAngle(45);
  multiarrow_sub = 0.02;
  arrow.Drawarrow(Bestop,BeE[9],Ostart-multiarrow_sub,OE[3]+Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[10],Ostart-multiarrow_sub,OE[3]+Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[11],Ostart-multiarrow_sub,OE[3]+Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[12],Ostart,OE[3]+Al22_gs+0.08,.02);
  arrow.Drawarrow(Bestop,BeE[13],Ostart-multiarrow_sub,OE[3]+Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[14],Ostart-multiarrow_sub,OE[3]+Al22_gs+0.3,.001);
  arrow.Drawarrow(Bestop,BeE[15],Ostart-multiarrow_sub,OE[3]+Al22_gs+0.3,.001);
  arrow.SetAngle(30);


  //Al22 1p decays to Mg21 g.s.
  multiarrow_sub = 0.04;
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.Drawarrow(Ostop,OE[2]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001); //Small branch to g.s.
  arrow.Drawarrow(Ostop,OE[3]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001);
  arrow.Drawarrow(Ostop,OE[4]+Al22_gs,Mgstart-0.005,MgE[0]+Mg21_gs+0.08,.02);
  arrow.Drawarrow(Ostop,OE[5]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001);
  arrow.Drawarrow(Ostop,OE[7]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001);
  arrow.SetAngle(30);

  //Al22 1p decays to Mg21 1/2+_1.
  multiarrow_sub = 0.035;
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.Drawarrow(Ostop,OE[1]+Al22_gs,Mgstart-0.005,MgE[1]+Mg21_gs+0.08,.02);
  arrow.Drawarrow(Ostop,OE[2]+Al22_gs,Mgstart-multiarrow_sub,MgE[1]+Mg21_gs+0.21,.001);
  arrow.SetAngle(30);

  //Al22 1p decays to Mg21 3/2+_1.
  multiarrow_sub = 0.035;
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.Drawarrow(Ostop,OE[6]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[7]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[8]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[9]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[10]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[11]+Al22_gs,Mgstart-0.005,MgE[3]+Mg21_gs+0.08,.02);
  arrow.Drawarrow(Ostop,OE[12]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[13]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[14]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[15]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.Drawarrow(Ostop,OE[16]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.SetAngle(30);


  //Al22 2p decays through Mg21.
  multiarrow_sub = 0.035;
  //arrow.SetLineStyle(2);
  arrow.SetAngle(45);
  arrow.Drawarrow(Ostop,OE[18]+Al22_gs,Mgstart-0.005,MgE[5]+Mg21_gs+0.08,.02);
  arrow.Drawarrow(Ostop,OE[17]+Al22_gs,Mgstart-multiarrow_sub,MgE[5]+Mg21_gs+0.35,.001);
  arrow.SetAngle(30);

  //Mg21 1p to Na20 4+_1.
  arrow.SetLineStyle(2);
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.Drawarrow(MgstBeP,MgE[5]+Mg21_gs,Nastart-0.005,NaE[2]+Na20_gs+0.06,.02);
  arrow.SetAngle(30);


  //Mg21 gamma decays
  //gamma decays
  TCurlyLine *gammadecay_Mg = new TCurlyLine(MgstBeP-0.15, MgE[3]+Mg21_gs, MgstBeP-0.15, Mg21_gs+0.3,.2, 0.015);
  gammadecay_Mg->SetLineColor(417);
  gammadecay_Mg->SetLineWidth(3);
  gammadecay_Mg->SetWavy();
  gammadecay_Mg->Draw();
  arrow.SetLineColor(417);
  arrow.SetFillColor(417); 
  arrow.SetAngle(45);
  arrow.Drawarrow(MgstBeP-0.15,Mg21_gs+0.15, MgstBeP-0.15,Mg21_gs+0.035,.013);

  //Na20 gamma decays
  //gamma decays
  //4+ to g.s.
  TCurlyLine *gammadecay_Na4p = new TCurlyLine(NastBeP-0.15, NaE[2]+Na20_gs, NastBeP-0.15, Na20_gs+0.3,.2, 0.015);
  gammadecay_Na4p->SetLineColor(417);
  gammadecay_Na4p->SetLineWidth(3);
  gammadecay_Na4p->SetWavy();
  gammadecay_Na4p->Draw();
  arrow.SetLineColor(417);
  arrow.SetFillColor(417); 
  arrow.SetAngle(45);
  arrow.Drawarrow(NastBeP-0.15,Na20_gs+0.15, NastBeP-0.15,Na20_gs+0.035,.013);
*/
  gStyle->SetLineScalePS(2);
  c1->Print("O12_Be12_MirrorDecayScheme.png", "png");
  c1->Print("O12_Be12_MirrorDecayScheme.pdf", "pdf");
  c1->Print("O12_Be12_MirrorDecayScheme.eps", "eps");
}
