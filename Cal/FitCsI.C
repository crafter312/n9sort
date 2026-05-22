void FitCsI()
{

  //Proton calibration files
  ifstream FCsICh("ProtonChannels.dat");
  if(!FCsICh.is_open())
    {
      cout << "could not open csi channel text file" << endl;
      return;
    }

  ifstream FCsIE("CsIProtonE.dat");
  if(!FCsIE.is_open())
    {
      cout <<"could not open csi energy text file " << endl;
      return;
    }

  ofstream cal_out("CsI_pcal_out.dat");

  //Alpha calibration files
  
  string bogus;
  getline(FCsIE,bogus);

  float CsIchannels[56][4] = {{0.}};
  float CsIenergy[56][4] ={{0.}};

  TCanvas *fitcan1 = (TCanvas*)gROOT->FindObject("fitcan1");
  if(!fitcan1)
    {
      fitcan1 = new TCanvas("fitcan1","fitcan1",1200,800);
      fitcan1->Divide(6,6);
    }
  TCanvas *fitcan2 = (TCanvas*)gROOT->FindObject("fitcan2");
  if(!fitcan2)
    {
      fitcan2 = new TCanvas("fitcan2","fitcan2",1200,800);
      fitcan2->Divide(6,6);
    }

  int NCsI =56;
  float channel1,channel2,channel3,channel4;
  float energy1,energy2,energy3,energy4;
  int csi1,csi2;

  
  
  for(int icsi = 0;icsi<NCsI;icsi++)
    {
      //read in values for centroids of peaks and energies of protons
      FCsICh>>csi1>>channel1>>channel2>>channel3>>channel4;
      FCsIE>>csi2>>energy1>>energy2>>energy3>>energy4;
      if(csi1 !=icsi || csi2 !=icsi)
	{
	  cout <<"didn't read the input files correctly" <<endl;
	  return;
	}
      
      //store these values in arrays
      CsIchannels[icsi][0] = channel1;
      CsIchannels[icsi][1] = channel2;
      CsIchannels[icsi][2] = channel3;
      CsIchannels[icsi][3] = channel4;

      CsIenergy[icsi][0] = energy4;
      CsIenergy[icsi][1] = energy3;
      CsIenergy[icsi][2] = energy2;
      CsIenergy[icsi][3] = energy1;

      
      TGraph *mygraph = new TGraph(4,CsIchannels[icsi],CsIenergy[icsi]);
      TF1 *myfit = new TF1("myfit","pol1",0,500);
      mygraph->SetMarkerStyle(20);
      mygraph->GetXaxis()->SetTitle("CsI [channel]");
      mygraph->GetYaxis()->SetTitle("CsI [MeV]");
      mygraph->GetXaxis()->SetLimits(0,500);
      mygraph->SetMinimum(0);
      mygraph->SetMaximum(80);
      if(icsi <36)
	{
	  fitcan1->cd(icsi+1);
	  mygraph->Draw("AP");
	  mygraph->Fit(myfit,"R");
	 
	}
      else
	{
	  fitcan2->cd(icsi-36+1);
	  mygraph->Draw("AP");
	  mygraph->Fit(myfit,"R");
	}

      cal_out<<0<<" "<<icsi<<" "<<myfit->GetParameter(1)<<" "<<myfit->GetParameter(0)<<endl;
      //      break;
      
    }

  FCsICh.close();
  FCsIE.close();
  cal_out.close();
  
  return;
}
