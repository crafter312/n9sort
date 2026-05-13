 /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * asicPeakSearch.C
 * A ROOT script to do peak search for the updated TECSA+HINP setup.
 *
 * Author: Antti J. Saastamoinen, Cyclotron Insititute, Texas A&M University (ajsaasta@comp.tamu.edu)
 * Uses ROOT so GNU LGPL applies. See more https://root.cern/about/license/
 *
 * Based on various old scripts from who knows where.
 *
 * First version: 2021-10-27
 * This version 2021-11-11
 *
 * Usage: root[] .L asicPeakSearch.C++                         <-- load this stuff
 * Possible uses:
 *        root[] findPeaksFromHistogramFile("/path/to/file/with/histograms.root","Up","Sector","5peak",1,kFALSE)  <-- process histograms in given file
 *        root[] findPeaksFromTree("/path/to/file/with/tree.root","Up","Sector","5peak",1,kFALSE)                 <-- process tree in given file
 *        root[] findPeaksFromChain(nameOfTheChain,"Up","Sector","5peak",1,kFALSE)                                <-- process tree in given chain of files
 *        arguments required:
 *          inputFile = file with the data tree, should work for any file with a TTree object. In case multiple trees per file, the first is taken.
 *          detector = Up/Down (defaults to Up)
 *          detectorSide = Sector/Ring (defaults to Sector)
 *          source = 5peak/Ra226/pulser (defaults to 5peak)
 *          rebinFactor = 1 (defaults to 1)
 *          testing = kFALSE/kTRUE (not really used after debugging, defaults to kFALSE)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TObject.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TChain.h"
#include <string>


TObject *obj;
TKey *key;
TFile *treeFile;
TFile *histoFile;
TChain *chainOfFiles;
TStopwatch *clocka;
TCanvas *canvas;


Double_t gausPeak(Double_t *x,Double_t *par){
  return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2));
}



void findPeaksFromHistogramFile(TString inputFile="HINPEs.root", TString detector = "Back", TString source="U232", Int_t rebinFactor=2, Bool_t testing=kTRUE){
  
  histoFile = TFile::Open(inputFile,"READ");
  cout << "Opening file: " << histoFile->GetName() << endl;
  TList *list = new TList;  


  Int_t numPeaks = 0;
  Int_t numChips = 0;
  Int_t numChans = 32; // fix this to default per chip
  Int_t numDetectorElements = 0;

  TString histoprefix = "SpecTcl_hinpa_mb1_e_";
  

  //Not using histopath for the spectcl tree files unpacked by HINPtree_unpack.C
  TString histoPath = ""; // histogram path in the TFile, hardcoded to match the offline sort

  if (detector.Contains("Front") == kTRUE) {
    numChips = 4;
    numDetectorElements = 32*4;
    histoPath = "Summary/1dFrontE_R/";
  }
   if (detector.Contains("Back") == kTRUE) {
    numChips = 4;
    numDetectorElements = 32*4;
    histoPath = "Summary/1dBackE_R/";
  }
  if (detector.Contains("Delta") == kTRUE) {
    numChips = 4;
    numDetectorElements = 32*4;
    histoPath = "Summary/1dDeltaE_R/";
  }
  
  if (source.Contains("Ra226") == kTRUE) {
    numPeaks = 5;
  }
  
  if (source.Contains("U232") == kTRUE) {
    numPeaks = 3;
  }

  if (source.Contains("5peak") == kTRUE) {
    numPeaks = 5;
  }

  if (source.Contains("LiAu") == kTRUE) {
    numPeaks = 1;
    //histoPath = "Summary/AngleCorrDeltaE/";
    histoPath = "Summary/AngleCorrFrontE/";
  }
  
  if (source.Contains("pulser") == kTRUE) {
    numPeaks = 21;
  }

  if (source.Contains("time") == kTRUE) {
    numPeaks = 1;
    if (detector.Contains("Front") == kTRUE)
      histoPath = "Summary/1dFrontTime_R/FrontTime_R";
    if (detector.Contains("Back") == kTRUE)
      histoPath = "Summary/1dBackTime_R/BackTime_R";
    if (detector.Contains("Delta") == kTRUE)
      histoPath = "Summary/1dDeltaTime_R/DeltaTime_R";
  }

  
  if(testing==kTRUE){
    cout << "numPeaks: " << numPeaks << " numChips: " << numChips << " numChans: " << numChans << endl;
  }
  
  
  if(numPeaks == 0 || numChips == 0 ){
    cout << "INVALID CONFIGURATION, CHECK YOUR ARGUMENTS!" << endl;
    return;
  }
  
  TString outFileName = "peakpositions_";
  
  outFileName.Append(detector);outFileName.Append("_");
  outFileName.Append(source);outFileName.Append(".txt");
  
  cout << "Opening peak positions file: " << outFileName << endl;
  ofstream fout;
  fout.open(outFileName);
  histoFile = TFile::Open(inputFile,"READ");
  
  TString histoname;
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Draw();
  c2->cd();
  //c2->SetLogy();
  c2->ToggleToolBar();
  c2->ToggleEventStatus();
  c2->ToggleEditor();
  
  TDirectory *current_sourcedir = gDirectory;
  
  for(int ib=0;ib<4;ib++) //we are using 4 boards
  {
    for(int ic=0;ic<numChans;ic++) //each has 32 chan
    {

      int nDetElement = (ib*32+1)+ic;

      //if (testing==kTRUE) {
        //cout << "histoPath: " << histoPath.Data() << " detector: " << detector.Data() << " nDetElemet: " << nDetElement << " board " << ib << " chan " << ic << endl;
      //}
    
      TString histoName;
      /*histoName.Append(histoPath);
      //histoName.Append(detector);
      //histoName.Append("E_R");
      //histoName.Append("AngleCorrE_R");
      //histoName.Append("FrontTime_R");
      histoName.Append(TString::Itoa(ib,10));
      histoName.Append("_");
      histoName.Append(TString::Itoa(ic,10));*/


      histoName.Append(histoprefix);
      if (detector.Contains("Front") == kTRUE) histoName.Append("0");
      if (detector.Contains("Front") == kTRUE) histoName.Append(std::to_string((ib+1)*2 -1));
      if (detector.Contains("Front") == kTRUE) histoName.Append(".");
      if (detector.Contains("Front") == kTRUE && ic < 10) histoName.Append("0");
      if (detector.Contains("Front") == kTRUE) histoName.Append(std::to_string(ic));

      if (detector.Contains("Back") == kTRUE) histoName.Append("0");
      if (detector.Contains("Back") == kTRUE) histoName.Append(std::to_string((ib+1)*2));
      if (detector.Contains("Back") == kTRUE) histoName.Append(".");
      if (detector.Contains("Back") == kTRUE && ic < 10) histoName.Append("0");
      if (detector.Contains("Back") == kTRUE) histoName.Append(std::to_string(ic));

      if (detector.Contains("Delta") == kTRUE && ib < 1) histoName.Append("0");
      if (detector.Contains("Delta") == kTRUE) histoName.Append(std::to_string(ib+9));
      if (detector.Contains("Delta") == kTRUE) histoName.Append(".");
      if (detector.Contains("Delta") == kTRUE && ic < 10) histoName.Append("0");
      if (detector.Contains("Delta") == kTRUE) histoName.Append(std::to_string(ic));

      if (testing==kTRUE) {
        cout << "Histogram name: " << histoName.Data() << endl;
      }

      TH1F *h1 = (TH1F*)gDirectory->Get(histoName.Data());
      
      h1->Rebin(rebinFactor);
      
      if(source.Contains("Ra226") == kTRUE && detector.Contains("Front") == kTRUE)
        h1->GetXaxis()->SetRangeUser(300,1500);
      if(source.Contains("Ra226") == kTRUE && detector.Contains("Back") == kTRUE)
        h1->GetXaxis()->SetRangeUser(400,1500);
      if(source.Contains("Ra226") == kTRUE && detector.Contains("Delta") == kTRUE)
        h1->GetXaxis()->SetRangeUser(300,1300);
      if(source.Contains("U232") == kTRUE && detector.Contains("Front") == kTRUE)
        h1->GetXaxis()->SetRangeUser(1000,1600);
      if(source.Contains("U232") == kTRUE && detector.Contains("Back") == kTRUE)
        h1->GetXaxis()->SetRangeUser(950,1600);
      if(source.Contains("U232") == kTRUE && detector.Contains("Delta") == kTRUE)
        h1->GetXaxis()->SetRangeUser(1000,1600);
      if(source.Contains("5peak") == kTRUE && detector.Contains("Front"))
        h1->GetXaxis()->SetRangeUser(200,1300);
      if(source.Contains("5peak") == kTRUE && detector.Contains("Back"))
        h1->GetXaxis()->SetRangeUser(400,1300);
      if(source.Contains("5peak") == kTRUE && detector.Contains("Delta"))
        h1->GetXaxis()->SetRangeUser(400,1050);
      if(source.Contains("LiAu") == kTRUE && detector.Contains("Front"))
        h1->GetXaxis()->SetRangeUser(4000,6000);
      if(source.Contains("LiAu") == kTRUE && detector.Contains("Delta"))
        h1->GetXaxis()->SetRangeUser(500,2000);
      if(source.Contains("time") == kTRUE && detector.Contains("Front"))
        h1->GetXaxis()->SetRangeUser(4000,7000);
      if(source.Contains("time") == kTRUE && detector.Contains("Back"))
        h1->GetXaxis()->SetRangeUser(3000,7000);
      if(source.Contains("time") == kTRUE && detector.Contains("Delta"))
        h1->GetXaxis()->SetRangeUser(4000,7000);
      if(source.Contains("pulser") == kTRUE)
        h1->GetXaxis()->SetRangeUser(500,16300);
      
      cout << "Processing " << h1->GetName() << " nbins=" << h1->GetNbinsX() << endl;
      c2->cd();
      h1->Draw();
      c2->Update();
      
      //Use TSpectrum to find the peak candidates
      cout << "Start searching for " << numPeaks << " peaks!" << endl;
      
      TSpectrum *s = new TSpectrum(numPeaks); // here TSpectrum has argument of max no. peaks. Also resolution can be set here as second arg, if not set default = 1.
      //s->SetResolution(.000000000001);
      // Int_t nfound = s->Search(h1,10,"new");

                             //hin, sigma, option, threshold
      Int_t nfound = s->Search(h1,1,"new,nobackground",0.05);

      list->Add(h1);

      //Int_t nfound = s->Search(h1,1,"new,nobackground",0.2);
      // printf("Found %d candidate peaks to fitn",nfound);
      Double_t *xpeaks = s->GetPositionX();
    
      const Int_t nparams = 3*numPeaks;
      Double_t par[3], parerr[3];
      Double_t parFitted[nparams], parFittedErrors[nparams];
    
      cout<<"number of peaks expected: "<< numPeaks << ", found: " <<nfound<< ", parameters to fit: " << nparams << endl;
      
      for(int i=nfound-1;i>=0;i--){
        Float_t xp = xpeaks[i];
        Int_t bin = h1->GetXaxis()->FindBin(xp);
        Float_t yp = h1->GetBinContent(bin);
        printf("xp = %f. yp = %f\n",xp,yp);
        
        par[0] = yp;
        par[1] = xp;
        par[2] = 1.;
        
        TF1 *peak = new TF1("peak",gausPeak,xp-10.,xp+10.,3);
        //TF1 *peak = new TF1("peak",gausPeak,xp-140.,xp+140.,3); //need wider fit range for Li beam
        //TF1 *peak = new TF1("peak",gausPeak,xp-50.,xp+50.,3); //setting used for timing

        peak->SetParameters(&par[0]);
        
        h1->Fit("peak","RQ+","",xp-10,xp+10);
        //h1->Fit("peak","RQ+","",xp-140,xp+140); //need wider fit range for Li beam
        //h1->Fit("peak","RQ+","",xp-50,xp+50); //setting used for timing

        // h1->Fit("peak","R");
        //gaus->GetParameters(par  histoFile = TFile::Open(inputFile,"READ"););
        peak->GetParameters(par);
        parerr[0] =  peak->GetParError(0);
        parerr[1] =  peak->GetParError(1);
        parerr[2] =  peak->GetParError(2);

        parFitted[i*3]=par[0];
        parFitted[i*3+1]=par[1];
        parFitted[i*3+2]=par[2];

        parFittedErrors[i*3]=parerr[0];
        parFittedErrors[i*3+1]=parerr[1];
        parFittedErrors[i*3+2]=parerr[2];

        
        // peak->Draw("same");
      }
      
      cout << endl;
      // cout<<parFitted[0]<<parFitted[8];
      
      TF1 *pk1 = new TF1("pk1",gausPeak,100,4000,3);
      pk1->SetParameters(&parFitted[0]);
      
      // pk1->Draw("same");
      
      Double_t XposUS[nfound];//unsorted x position array
      Double_t XposUSerr[nfound];
    
      for (int i=0;i<nfound;i++){
        XposUS[i]=parFitted[i*3+1];
        XposUSerr[i]=parFittedErrors[i*3+1];
        //        XposUS[i]=xpeaks[i];
      }
      
      Int_t indicies[nfound];
      for(int i=0;i<nfound;i++){indicies[i]=0;}
      TMath::BubbleHigh(nfound,XposUS,indicies);
      Double_t Xpos[nfound];
      Double_t XposErr[nfound];
    
      cout <<"fitted centroids: ";
      cout << " chan: " << nDetElement << " ";
      fout << nDetElement << " ";
    
      if (nfound != numPeaks) {
        cout << " nfound NOT number of peaks! Found : " << nfound << " expected: " << numPeaks << " Setting all values to zero for safety! " << endl;
        for (Int_t i=0; i<numPeaks; i++) {
          Xpos[i]=0.0;
          XposErr[i]=0.0;
          cout<<Xpos[i]<<"+/-" << XposErr[i] << " , " ;
          fout<<Xpos[i]<<" "<<XposErr[i]<<" ";
        }
      } else {
        for (int i=0;i<nfound;i++){
          Xpos[i]=XposUS[indicies[i]];
          XposErr[i]=XposUSerr[indicies[i]];
          cout << Xpos[i] << "+/-" << XposErr[i] << " , " ;
          fout << Xpos[i] << " " << XposErr[i] << " ";
        }
      }

      cout << endl << endl;
      fout << endl;

    }
  } // end loop over channels per chip
  
  
  fout.close();

  TFile * outfile = new TFile("peaksearch.root", "RECREATE");
  outfile->cd();
  list->Write();
  outfile->Close();
  
} // end findPeaksFromHistogramFile
