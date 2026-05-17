void CsICal()
{
  gROOT->Reset();
  
  //Proton cal runs
  //TFile* file = new TFile("/mnt/analysis/e25001/rootout/sort_cal_p_nothing.root");
  //TFile* file = new TFile("/mnt/analysis/e25001/rootout/sort_cal_p_Al0-25.root");
  TFile* file = new TFile("/mnt/analysis/e25001/rootout/sort_cal_p_Al0-50.root");
  
  //Alpha cal runs
  //TFile* file = new TFile("/mnt/analysis/e25001/rootout/sort_cal_a_nothing.root");
  
  //Proton out files
  //ofstream ofile("/user/e25001/analysis/n9sort/CsI_p_nothing.dat");
  //ofstream ofile("/user/e25001/analysis/n9sort/CsI_p_Be.dat");
  //ofstream ofile("/user/e25001/analysis/n9sort/CsI_p_Al0-25.dat");
  ofstream ofile("/user/e25001/analysis/n9sort/CsI_p_Al0-50.dat");
  
  //Alpha out files
  //ofstream ofile("/user/e25001/analysis/n9sort/CsI_a_nothing.dat");
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1);
  
  ostringstream outstring;
  ostringstream histname;
  string name;
  
  int p1= 30, p2=50; //+- fit limits up to 2 peaks. May be different.
  int const num_par = 5; //number of peaks times 2(pol1)+3(gaus).
  
  int ent = 0;
  
  gROOT->cd();
  TCanvas *mycan =(TCanvas*)gROOT->FindObjectAny("mycan");
  if(!mycan)
    {
      mycan = new TCanvas("mycan","mycan");
    }
  
  for(int it=0;it<4;it++)
    {
    for (int ic=0;ic<7;ic++) {
      outstring.str("");
      histname.str("");
      //Use the center of the CsI for some crystals. Will vary depending on proton spread
      //As Al thickness increases, outer crystals get more hits
      if (ic > 3)
        outstring << "Summary/1dCsI_Energy/CsI_Energy_" << it << "_" << ic;
      else
        outstring << "Summary/1dCsI_Energy/CsI_Energy_R_center_" << it << "_" << ic;
      name = outstring.str();
      histname << "ECsIRaw_" << it << "_" << ic << endl;
      TH1I * hist = (TH1I*)file->Get(name.c_str())->Clone(Form(histname.str().c_str()));
      hist->Draw("");


//      hist->Rebin(2);
      hist->GetXaxis()->SetRangeUser(0,1000);

      
      
      
      mycan->Modified();
      mycan->Update();
      cout << "tele = "<< it << " crystal = " << ic << endl;
      
      
      TMarker * mark;
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Get the Background limits
      int bkg_lo = mark->GetX();
      delete mark;  
      mark=(TMarker*)mycan->WaitPrimitive("TMarker");
      int bkg_hi = mark->GetX();
      delete mark;
      cout << "bkg_lo " << bkg_lo << ", bkg_hi " << bkg_hi << endl;
      
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      int peak1 = mark->GetX();
      delete mark;
      
      
      double par[num_par] = {0.};
      double out[num_par] = {0.}; 
      int peak1_lo = peak1 - p1, peak1_hi = peak1 + p1; // Peak center and limits
      cout << "peak1 " << peak1 << ", peak1_lo " << peak1_lo << ", peak1_hi " << peak1_hi << endl;
      
      
      // TF1 *l1 = new TF1("l1", "pol1", bkg_lo, bkg_hi);
      // TF1 *g1 = new TF1("g1", "gaus", peak1_lo,peak1_hi);
      
      TF1 *total = new TF1("total", "([0]*x)+[1]+([2]*TMath::Gaus(x, [3], [4]))", bkg_lo,bkg_hi); // [0] is background slope, [1] is background offset, [2] is Gaus scaling, [3] is Gaus mean, [4] is Gaus sigma
      //TF1 *total = new TF1("total", "gaus(0)", peak1_lo,peak1_hi);
      
      // hist->Fit(l1,"R");
      // hist->Fit(g1,"R+");
      
      // l1->GetParameters(&par[0]);
      // g1->GetParameters(&par[2]);
      
      total->SetParameters(0,0);
      total->SetParameters(1,0);
      total->SetParameters(2,.8*hist->GetMaximum());
      total->SetParameters(3,peak1);
      total->SetParameters(4,8);
      total->SetParLimits(2,.3*hist->GetMaximum(),hist->GetMaximum());
      total->SetParLimits(3,peak1_lo,peak1_hi);
      total->SetParLimits(4,2,12);
      hist->Fit(total,"R+");
      total->GetParameters(out);
      
      
       //   total->Write();
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      delete mark;      


      ofile << it << " " << ic << " " << out[3] << endl;
      cout << it << " " << ic << " " << out[3] << endl;
      }
    }

  ofile.close();
  file->Close();
  
  return;
}
