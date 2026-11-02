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

  Sty->SetNdivisions(12,"y");
  Sty->SetTitleOffset(0.7,"y");
  Sty->SetTitleSize(0.07,"y");

  TH2S frame2("frame","",10,0,2,10,-1,11);
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
  double gap = 0.3;

  double Mg21_gamma = 1.598;

  //Get masses
  double mass_p = 7.28897;
  double mass_20Na = 6.8505;
  double mass_21Mg = 10.9039;
  double mass_22Al = 18.092;
  double mass_23Si = 23.3629;

  //Get ground states
  double Al22_gs = -(mass_23Si - mass_p - mass_22Al);
  double Mg21_gs = -(mass_23Si - 2.*mass_p - mass_21Mg);
  double Na20_gs = -(mass_23Si - 3.*mass_p - mass_20Na);

  cout << "Calculated ground states" << endl;
  cout << "Al22 + p: " << Al22_gs << endl;
  cout << "Mg21 + 2p: " << Mg21_gs << endl;
  cout << "Na20 + 3p: " << Na20_gs << endl;

  //Read in level files from /levels/

  //Read in 23Si states 
  string filename_exp = "levels/Si23_exp.dat";
  ifstream file_exp(filename_exp.c_str());
  if (file_exp.is_open() != 1)
  {
    cout << " could not open file " << filename_exp << endl;
    abort();
  }

  //Level files have a 3 line header, skip all of this
  for (int i=0;i<3;i++) file_exp.ignore(200,'\n');

  string readvals[7];
  vector<string> SiJ;
  vector<int> SiP;
  vector<int> SiN;
  vector<double> SiE;
  vector<double> SiW;
  vector<string> SiDM;
  vector<int> SiK; //ENSDF or our experiment? changes level coloring

  //Read in values
  string fline;
  while (getline(file_exp,fline)) {
    stringstream ss(fline);
    ss >> readvals[0] >> readvals[1] >> readvals[2] >> readvals[3] >> readvals[4] >> readvals[5] >> readvals[6];

    SiJ.push_back(readvals[0]);
    SiP.push_back(stoi(readvals[1]));
    SiN.push_back(stoi(readvals[2]));
    SiE.push_back(stod(readvals[3])/1000.);
    SiW.push_back(stod(readvals[4])/1000.);
    SiDM.push_back(readvals[5]);
    SiK.push_back(stoi(readvals[6]));
  }

  file_exp.close();
  file_exp.clear();

  //Read in 22Al states 
  filename_exp = "levels/Al22_exp.dat";
  file_exp.open(filename_exp.c_str());
  if (file_exp.is_open() != 1)
  {
    cout << " could not open file " << filename_exp << endl;
    abort();
  }

  //Level files have a 3 line header, skip all of this
  for (int i=0;i<3;i++) file_exp.ignore(200,'\n');

  vector<string> AlJ;
  vector<int> AlP;
  vector<int> AlN;
  vector<double> AlE;
  vector<double> AlW;
  vector<string> AlDM;
  vector<int> AlK; //ENSDF or our experiment? changes level coloring

  //Read in values
  while (getline(file_exp,fline)) {
    stringstream ss(fline);
    ss >> readvals[0] >> readvals[1] >> readvals[2] >> readvals[3] >> readvals[4] >> readvals[5] >> readvals[6];

    AlJ.push_back(readvals[0]);
    AlP.push_back(stoi(readvals[1]));
    AlN.push_back(stoi(readvals[2]));
    AlE.push_back(stod(readvals[3])/1000.);
    AlW.push_back(stod(readvals[4])/1000.);
    AlDM.push_back(readvals[5]);
    AlK.push_back(stoi(readvals[6]));
  }

  file_exp.close();
  file_exp.clear();



  //Read in 21Mg states 
  filename_exp = "levels/Mg21_exp.dat";
  file_exp.open(filename_exp.c_str());
  if (file_exp.is_open() != 1)
  {
    cout << " could not open file " << filename_exp << endl;
    abort();
  }

  //Level files have a 3 line header, skip all of this
  for (int i=0;i<3;i++) file_exp.ignore(200,'\n');

  vector<string> MgJ;
  vector<int> MgP;
  vector<int> MgN;
  vector<double> MgE;
  vector<double> MgW;
  vector<string> MgDM;
  vector<int> MgK; //ENSDF or our experiment? changes level coloring

  //Read in values
  while (getline(file_exp,fline)) {
    stringstream ss(fline);
    ss >> readvals[0] >> readvals[1] >> readvals[2] >> readvals[3] >> readvals[4] >> readvals[5] >> readvals[6];

    MgJ.push_back(readvals[0]);
    MgP.push_back(stoi(readvals[1]));
    MgN.push_back(stoi(readvals[2]));
    MgE.push_back(stod(readvals[3])/1000.);
    MgW.push_back(stod(readvals[4])/1000.);
    MgDM.push_back(readvals[5]);
    MgK.push_back(stoi(readvals[6]));
  }

  file_exp.close();
  file_exp.clear();



  //Read in 20Na states 
  filename_exp = "levels/Na20_exp.dat";
  file_exp.open(filename_exp.c_str());
  if (file_exp.is_open() != 1)
  {
    cout << " could not open file " << filename_exp << endl;
    abort();
  }

  //Level files have a 3 line header, skip all of this
  for (int i=0;i<3;i++) file_exp.ignore(200,'\n');

  vector<string> NaJ;
  vector<int> NaP;
  vector<int> NaN;
  vector<double> NaE;
  vector<double> NaW;
  vector<string> NaDM;
  vector<int> NaK; //ENSDF or our experiment? changes level coloring

  //Read in values
  while (getline(file_exp,fline)) {
    stringstream ss(fline);
    ss >> readvals[0] >> readvals[1] >> readvals[2] >> readvals[3] >> readvals[4] >> readvals[5] >> readvals[6];

    NaJ.push_back(readvals[0]);
    NaP.push_back(stoi(readvals[1]));
    NaN.push_back(stoi(readvals[2]));
    NaE.push_back(stod(readvals[3])/1000.);
    NaW.push_back(stod(readvals[4])/1000.);
    NaDM.push_back(readvals[5]);
    NaK.push_back(stoi(readvals[6]));
  }

  file_exp.close();
  file_exp.clear();


  //Plot experimental levels and widths
  //Add g.s. separation energies

  //Si23
  double Sistart = indentX;
  double Sistop = indentX+LevelXwide;
  for (int i=0;i<SiJ.size();i++) {

    int levelcol;
    if (SiK[i] == 1) levelcol = 1;
    else levelcol = 4;   

    ll.SetLineColor(levelcol);  

    ll.DrawLine(Sistart,SiE[i],Sistop,SiE[i]);

    TBox *expbox = new TBox(Sistart,SiE[i]-SiW[i]/2.,Sistop,SiE[i]+SiW[i]/2.);
    expbox->SetBit(kCanDelete);
    expbox->SetFillColor(levelcol);
    expbox->SetFillStyle(3004);
    expbox->Draw("same"); 
  }

  //Al22
  double Alstart = indentX+LevelXwide+gap;
  double Alstop = indentX+2.*LevelXwide+gap;
  for (int i=0;i<AlJ.size();i++) {

    int levelcol;
    if (AlK[i] == 1) levelcol = 1;
    else levelcol = 4;   

    ll.SetLineColor(levelcol);  

    ll.DrawLine(Alstart,AlE[i]+Al22_gs,Alstop,AlE[i]+Al22_gs);

    TBox *expbox = new TBox(Alstart,AlE[i]+Al22_gs-AlW[i]/2.,Alstop,AlE[i]+Al22_gs+AlW[i]/2.);
    expbox->SetBit(kCanDelete);
    expbox->SetFillColor(levelcol);
    expbox->SetFillStyle(3004);
    expbox->Draw("same"); 
  }

  //Mg21
  double Mgstart = indentX+2.*LevelXwide+2.*gap;
  double Mgstop = indentX+3.*LevelXwide+2.*gap;
  for (int i=0;i<MgJ.size();i++) {

    int levelcol;
    if (MgK[i] == 1) levelcol = 1;
    else levelcol = 4;   

    ll.SetLineColor(levelcol);  

    ll.DrawLine(Mgstart,MgE[i]+Mg21_gs,Mgstop,MgE[i]+Mg21_gs);

    TBox *expbox = new TBox(Mgstart,MgE[i]+Mg21_gs-MgW[i]/2.,Mgstop,MgE[i]+Mg21_gs+MgW[i]/2.);
    expbox->SetBit(kCanDelete);
    expbox->SetFillColor(levelcol);
    expbox->SetFillStyle(3004);
    expbox->Draw("same"); 
  }

  //Na20
  double Nastart = indentX+3.*LevelXwide+3.*gap;
  double Nastop = indentX+4.*LevelXwide+3.*gap;
  for (int i=0;i<NaJ.size();i++) {

    int levelcol;
    if (NaK[i] == 1) levelcol = 1;
    else levelcol = 4;   

    ll.SetLineColor(levelcol);  

    ll.DrawLine(Nastart,NaE[i]+Na20_gs,Nastop,NaE[i]+Na20_gs);

    TBox *expbox = new TBox(Nastart,NaE[i]+Na20_gs-NaW[i]/2.,Nastop,NaE[i]+Na20_gs+NaW[i]/2.);
    expbox->SetBit(kCanDelete);
    expbox->SetFillColor(levelcol);
    expbox->SetFillStyle(3004);
    expbox->Draw("same"); 
  }

 


  tt.SetTextSize(0.06);
  Float_t textsize = tt.GetTextSize();

  //labels
  tt.DrawLatex(Sistop - 0.16 ,SiE[0]-0.9,"^{23}Si");
  tt.DrawLatex(Alstop - 0.2025 ,AlE[0]+Al22_gs-0.9,"#it{p}+#kern[-0.01]{^{22}Al}");
  tt.DrawLatex(Mgstop - 0.225 ,MgE[0]+Mg21_gs-0.9,"2#it{p}+#kern[-0.01]{^{21}Mg}");
  tt.DrawLatex(Nastop - 0.23 ,NaE[0]+Na20_gs-0.9,"3#it{p}+#kern[-0.01]{^{20}Na}");

  tt.SetTextSize(0.045);
  //Draw state J pi labels
  //Si23
  for (int i=0;i<SiJ.size();i++) {

    //if no spin ("UN"), skip
    if (SiJ[i] == "UN") continue;

    string tempJ = SiJ[i];
    if (SiP[i] == 1) tempJ += "^{+}";
    else tempJ += "^{-}";

    tt.DrawLatex(Sistop-0.08,SiE[i]+0.1,tempJ.c_str());
  }


  //Al22
  for (int i=0;i<AlJ.size();i++) {

    //if no spin ("UN"), skip
    if (AlJ[i] == "UN") continue;

    string tempJ = AlJ[i];
    if (AlP[i] == 1) tempJ += "^{+}";
    else tempJ += "^{-}"; 

    tt.DrawLatex(Alstop-0.08,AlE[i]+Al22_gs+0.1,tempJ.c_str());
  }


  //Mg21
  for (int i=0;i<MgJ.size();i++) {

    //if no spin ("UN"), skip
    if (MgJ[i] == "UN") continue;

    string tempJ = MgJ[i];
    if (MgP[i] == 1) tempJ += "^{+}";
    else tempJ += "^{-}"; 

    tt.DrawLatex(Mgstop-0.08,MgE[i]+Mg21_gs+0.1,tempJ.c_str());
  }


  //Na20
  for (int i=0;i<NaJ.size();i++) {

    //if no spin ("UN"), skip
    if (NaJ[i] == "UN") continue;

    string tempJ = NaJ[i];
    if (NaP[i] == 1) tempJ += "^{+}";
    else tempJ += "^{-}"; 

    tt.DrawLatex(Nastop-0.08,NaE[i]+Na20_gs+0.1,tempJ.c_str());
  }


  arrow.SetAngle(30);
  arrow.SetLineWidth(3);
  arrow.SetLineColor(1);
  arrow.SetFillColor(1); 

  arrow.SetLineColor(kRed);
  arrow.SetFillColor(kRed); 

  //Si23 1p decays to Al22 g.s.
  float multiarrow_sub = 0.04;
  arrow.SetAngle(45);
  arrow.DrawArrow(Sistop,SiE[1],Alstart-multiarrow_sub,Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[2],Alstart-multiarrow_sub,Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[3],Alstart-multiarrow_sub,Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[4],Alstart,Al22_gs+0.08,.02);
  arrow.DrawArrow(Sistop,SiE[5],Alstart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.DrawArrow(Sistop,SiE[6],Alstart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.DrawArrow(Sistop,SiE[7],Alstart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.DrawArrow(Sistop,SiE[8],Alstart-multiarrow_sub,Al22_gs+0.26,.001);
  arrow.SetAngle(30);

  //Si23 2p decays to Al22 g.s.
  //arrow.SetLineStyle(2);
  arrow.SetAngle(45);
  multiarrow_sub = 0.02;
  arrow.DrawArrow(Sistop,SiE[9],Alstart-multiarrow_sub,AlE[3]+Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[10],Alstart-multiarrow_sub,AlE[3]+Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[11],Alstart-multiarrow_sub,AlE[3]+Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[12],Alstart,AlE[3]+Al22_gs+0.08,.02);
  arrow.DrawArrow(Sistop,SiE[13],Alstart-multiarrow_sub,AlE[3]+Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[14],Alstart-multiarrow_sub,AlE[3]+Al22_gs+0.3,.001);
  arrow.DrawArrow(Sistop,SiE[15],Alstart-multiarrow_sub,AlE[3]+Al22_gs+0.3,.001);
  arrow.SetAngle(30);


  //Al22 1p decays to Mg21 g.s.
  multiarrow_sub = 0.04;
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.DrawArrow(Alstop,AlE[2]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001); //Small branch to g.s.
  arrow.DrawArrow(Alstop,AlE[3]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001);
  arrow.DrawArrow(Alstop,AlE[4]+Al22_gs,Mgstart-0.005,MgE[0]+Mg21_gs+0.08,.02);
  arrow.DrawArrow(Alstop,AlE[5]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001);
  arrow.DrawArrow(Alstop,AlE[7]+Al22_gs,Mgstart-multiarrow_sub,MgE[0]+Mg21_gs+0.2,.001);
  arrow.SetAngle(30);

  //Al22 1p decays to Mg21 1/2+_1.
  multiarrow_sub = 0.035;
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.DrawArrow(Alstop,AlE[1]+Al22_gs,Mgstart-0.005,MgE[1]+Mg21_gs+0.08,.02);
  arrow.DrawArrow(Alstop,AlE[2]+Al22_gs,Mgstart-multiarrow_sub,MgE[1]+Mg21_gs+0.21,.001);
  arrow.SetAngle(30);

  //Al22 1p decays to Mg21 3/2+_1.
  multiarrow_sub = 0.035;
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.DrawArrow(Alstop,AlE[6]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[7]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[8]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[9]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[10]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[11]+Al22_gs,Mgstart-0.005,MgE[3]+Mg21_gs+0.08,.02);
  arrow.DrawArrow(Alstop,AlE[12]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[13]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[14]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[15]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.DrawArrow(Alstop,AlE[16]+Al22_gs,Mgstart-multiarrow_sub,MgE[3]+Mg21_gs+0.21,.001);
  arrow.SetAngle(30);


  //Al22 2p decays through Mg21.
  multiarrow_sub = 0.035;
  //arrow.SetLineStyle(2);
  arrow.SetAngle(45);
  arrow.DrawArrow(Alstop,AlE[18]+Al22_gs,Mgstart-0.005,MgE[5]+Mg21_gs+0.08,.02);
  arrow.DrawArrow(Alstop,AlE[17]+Al22_gs,Mgstart-multiarrow_sub,MgE[5]+Mg21_gs+0.35,.001);
  arrow.SetAngle(30);

  //Mg21 1p to Na20 4+_1.
  arrow.SetLineStyle(2);
  arrow.SetLineStyle(0);
  arrow.SetAngle(45);
  arrow.DrawArrow(Mgstop,MgE[5]+Mg21_gs,Nastart-0.005,NaE[2]+Na20_gs+0.06,.02);
  arrow.SetAngle(30);


  //Mg21 gamma decays
  //gamma decays
  TCurlyLine *gammadecay_Mg = new TCurlyLine(Mgstop-0.15, MgE[3]+Mg21_gs, Mgstop-0.15, Mg21_gs+0.3,.2, 0.015);
  gammadecay_Mg->SetLineColor(417);
  gammadecay_Mg->SetLineWidth(3);
  gammadecay_Mg->SetWavy();
  gammadecay_Mg->Draw();
  arrow.SetLineColor(417);
  arrow.SetFillColor(417); 
  arrow.SetAngle(45);
  arrow.DrawArrow(Mgstop-0.15,Mg21_gs+0.15, Mgstop-0.15,Mg21_gs+0.035,.013);

  //Na20 gamma decays
  //gamma decays
  //4+ to g.s.
  TCurlyLine *gammadecay_Na4p = new TCurlyLine(Nastop-0.15, NaE[2]+Na20_gs, Nastop-0.15, Na20_gs+0.3,.2, 0.015);
  gammadecay_Na4p->SetLineColor(417);
  gammadecay_Na4p->SetLineWidth(3);
  gammadecay_Na4p->SetWavy();
  gammadecay_Na4p->Draw();
  arrow.SetLineColor(417);
  arrow.SetFillColor(417); 
  arrow.SetAngle(45);
  arrow.DrawArrow(Nastop-0.15,Na20_gs+0.15, Nastop-0.15,Na20_gs+0.035,.013);

  gStyle->SetLineScalePS(2);
  c1->Print("Si23_Al22_DecayScheme.png", "png");
  c1->Print("Si23_Al22_DecayScheme.pdf", "pdf");
  c1->Print("Si23_Al22_DecayScheme.eps", "eps");
}
