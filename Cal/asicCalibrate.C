/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * asicCalibrate.C
 * A ROOT script to do peak search for the updated TECSA+HINP setup.
 *
 * Author: Antti J. Saastamoinen, Cyclotron Insititute, Texas A&M University (ajsaasta@comp.tamu.edu)
 * Uses ROOT so GNU LGPL applies. See more https://root.cern/about/license/
 *
 * Based on various old scripts from who knows where.
 *
 * First version: 2021-11-04
 * This version 2021-11-11
 *
 * Usage: root[] .L asicCalibrate.C++                         <-- load this stuff
 * Possible uses:
 *        root[] calibrate("peak_positions_from_peaksearch.txt","Up","Sector","Ra226",1)    <-- calibrate whatever detector / elements with sources
 *
 *  arguments required:
 *          inputFile = file produced by asicPeakSearch.C
 *          detector = Up/Down (defaults to Up)
 *          detectorSide = Sector/Ring (defaults to Sector)
 *          source = 4peak/Ra226/pulser (defaults to 4peak)
 *          polDeg = degree of the fitted polynomial (defaults to 1, i.e. linear)
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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2.h"



void calibrate(TString detector="Back"){

  TCanvas *cChi;
  cChi = new TCanvas("cChi","StripChiSqds");
  TCanvas *cRes;
  cRes = new TCanvas("cRes","Residuals");
  TH2F *hResiduals = new TH2F("hResiduals","Residuals",128,-0.5,127.5,4096,0,16383);

  TList *list = new TList;

  int polDeg = 1;

  TString inputFile1;
  TString inputFile2;
  TString inputFile3;

  if (detector.Contains("Front"))
  {
    inputFile1 = "peakpositions_Front_U232_save.txt";
    //inputFile2 = "peakpositions_Front_Ra226.txt";
    //inputFile3 = "peakpositions_Front_LiAu.txt";
  }
  if (detector.Contains("Back"))
  {
    inputFile1 = "peakpositions_Back_U232_saveMay12.txt";
    //inputFile2 = "peakpositions_Back_Ra226.txt";
  }
  if (detector.Contains("Delta"))
  {
    inputFile1 = "peakpositions_Delta_5peak_saveMarch16.txt";
    //inputFile2 = "peakpositions_Delta_Ra226.txt";
    //inputFile3 = "peakpositions_Delta_LiAu.txt";
  }

  cout << "Read in: inputFile = " << inputFile1 << " and File " << inputFile2 << ", detector = " << detector << ", polynomial degree = " << polDeg << endl;
  
  //FrontEcal.dat
  TString outFileName = "";
  outFileName.Append(detector);
  outFileName.Append("Ecal.dat");
  
  cout << "Opening calibration results file: " << outFileName << endl;
  ofstream fout;
  fout.open(outFileName);
  
  int numDetectorElements = 4*32;
  int numPeaks = 3;

  //if (detector.Contains("Front") || detector.Contains("Delta"))
  //{
    //numPeaks++;
  //}
  
  Double_t chisqdarray[numDetectorElements], stripno[numDetectorElements], residualarray[numDetectorElements][numPeaks], peaks[numPeaks], alphastemp[numPeaks], alphas[numDetectorElements][numPeaks], peakerr[numPeaks], alphastemperr[numPeaks], alphaserr[numDetectorElements][numPeaks];
  Double_t tempcalcval[numPeaks];
  Double_t tenparray[numPeaks][numDetectorElements+1];
  Double_t tenparrayerr[numPeaks][numDetectorElements+1];
//  double alphasDL[4][24];
  // 239Pu value seems wrong in this (found this from old SAMUTAMU script)
//  alphas[0] = 5795.,alphas[1] = 5479.3,alphas[2] = 5130.4,alphas[3] = 3182.7; // lit. values
//  alphas[0] = 5719.3,alphas[1] = 5400.3,alphas[2] = 5047.4,alphas[3] = 3073.8; // after 0.5 um Al window
  
  // Lit. values as per NNDC. In case of Cm, Am and Pu energies are weighted averages of major lines as those do not resolve to separate lines with TTT.
  //  alphas[0] = 5795.,alphas[1] = 5479.3,alphas[2] = 5148.4,alphas[3] = 3182.7;
//  alphas[0] = 5719.3,alphas[1] = 5400.3,alphas[2] = 5066.4,alphas[3] = 3073.8; // after 0.5 um Al window (LISE++ phys. calc. for loss calculation)

  //TGraph *g1 = new TGraph(numPeaks,peaks,alphastemp);
  TGraph *g2 = new TGraph(numDetectorElements,stripno,chisqdarray);
  //TGraphErrors *g1err = new TGraphErrors(numPeaks,peaks,alphastemp,peakerr,alphastemperr);
  
  for (Int_t i=0; i<numDetectorElements; i++) {
   
    //alphas[i][0] = 5804.,alphas[i][1]=5762.,alphas[i][2] = 5479.3, alphas[i][3] = 5153.6, alphas[i][4] = 3182.7; // Lit. values as per NNDC. Am and Pu energies are weighted averages      alphas[i][0] = 5795.,alphas[i][1] = 5479.3, alphas[i][2] = 5148.4, alphas[i][3] = 3182.7; // Lit. values as permajor lines as those do not resolve to separate lines
    //alphas[i][0] = 5804.,alphas[i][1]=5804.,alphas[i][2] = 5479.3, alphas[i][3] = 5153.6, alphas[i][4] = 3182.7; // Lit. values as per NNDC. Am and Pu energies are weighted averages      alphas[i][0]
    //alphas[i][0] = 7686,alphas[i][1]=6002.,alphas[i][2] = 5489., alphas[i][3] = 5304., alphas[i][4] = 4784.;
    //alphaserr[i][0] = 1.0, alphaserr[i][1] = 1.0,alphaserr[i][2] = 1.0,alphaserr[i][3] = 1.0,alphaserr[i][4] = 1.0;

    //alphas[i][4] = 7686.82, alphas[i][5] = 6002.55, alphas[i][6] = 5489.48, alphas[i][7] = 5304.33, alphas[i][8] = 4784.34;
    //alphaserr[i][4] = 1.0, alphaserr[i][5] = 1.0,alphaserr[i][6] = 1.0,alphaserr[i][7] = 1.0,alphaserr[i][8] = 1.0;

    //alphas[i][9] = 34369.04*0.98;
    
    //if (detector.Contains("Delta"))
      //alphas[i][9] = 7665.0;

    //alphaserr[i][9] = 1.0;

    alphas[i][0] = 8767.,alphas[i][1]=6758.,alphas[i][2] = 6267.; //Po-212, Po-216, and Rn-220 in U-232 chain
    alphaserr[i][0] = 1.0, alphaserr[i][1] = 1.0,alphaserr[i][2] = 1.0;
  }


  ifstream fin1;
  fin1.open(inputFile1);
  //ifstream fin2;
  //fin2.open(inputFile2);
  //ifstream fin3;
  //if (detector.Contains("Front") || detector.Contains("Delta"))
  //{
    //fin3.open(inputFile3);
  //}

  for(Int_t i=0; i<numDetectorElements; i++){
    cout << "\n\nDoing detector element " << i+1 << endl;
    stripno[i] = i+1;
    Int_t ch;
    fin1 >> ch;
    //fin2 >> ch;
    //fin3 >> ch;

    for (Int_t j=0; j<numPeaks; j++) {
      fin1 >> peaks[j] >> peakerr[j];
      peakerr[j] = 1;
      tenparray[j][i+1] = peaks[j];
      tenparrayerr[j][i+1] = peakerr[j];
    }
    //for (Int_t j=0; j<5; j++) {
      
      //fin2 >> peaks[j+4] >> peakerr[j+4];
      //peakerr[j+4] = 1;
      //tenparray[j+4][i+1] = peaks[j+4];
      //tenparrayerr[j+4][i+1] = peakerr[j+4];
    //}

    //only keep if you are doing fronts with LiAu file
    //fin3 >> peaks[9] >> peakerr[9];
    //peakerr[9] = 1;
    //tenparray[9][i+1] = peaks[9];
    //tenparrayerr[9][i+1] = peakerr[9];

    bool zeros = false;
    for (Int_t j=0; j<numPeaks; j++) {
      if (peaks[j] == 0)
        zeros = true;
    }
    if (zeros) {
      int board0 = i/32;
      int chip0 = i%32;
      fout << board0 << " " << chip0 << " " << 0 << " " << 0 << endl;
      continue;
    } 

    cout << "read peaks: ";
    for (Int_t j=0; j<numPeaks; j++){ cout << peaks[j] << " " << peakerr[j] << "  " <<  tenparray[j][i+1] << " " << tenparrayerr[j][i+1] << "  " ;}
   // for (Int_t j=0; j<5; j++){ cout << peaks[j+4] << " " << peakerr[j+4] << "  " <<  tenparray[j+4][i+1] << " " << tenparrayerr[j+4][i+1] << "  " ;}
  //  cout << peaks[9] << " " << peakerr[9] << "  " <<  tenparray[9][i+1] << " " << tenparrayerr[9][i+1] << "  " ;

    cout << "\n";
      
    // hack here to get the energies to smaller array expected by TGraph...
    
    for (Int_t j=0; j<numPeaks; j++) {
      alphastemp[j] = alphas[i][j];
      tenparray[j][0] = alphastemp[j];
      alphastemperr[j] = alphaserr[i][j];
      tenparrayerr[j][0] = alphastemperr[j];
    }
    /*for (Int_t j=0; j<5; j++) {
      alphastemp[j+4] = alphas[i][j+4];
      tenparray[j+4][0] = alphastemp[j+4];
      alphastemperr[j+4] = alphaserr[i][j+4];
      tenparrayerr[j+4][0] = alphastemperr[j+4];
    }
    alphastemp[9] = alphas[i][9];
    tenparray[9][0] = alphastemp[9];
    alphastemperr[9] = alphaserr[i][9];
    tenparrayerr[9][0] = alphastemperr[9];*/


    cout << "alphas temp peaks: ";
    for (Int_t j=0; j<numPeaks; j++){ cout << alphastemp[j] << " " << alphastemperr[j] << "  ";}
    cout << "\n";
      
   // TGraph *g1 = new TGraph(numPeaks,peaks,alphastemp);
    //TGraphErrors *g1err = new TGraphErrors(numPeaks,peaks,alphastemp,peakerr,alphastemperr);
    TGraphErrors *g1err = new TGraphErrors();
      //TGraph *g2 = new TGraph(24,stripno,chisqdarray);

    
    for (Int_t j=0; j<numPeaks; j++) {
      cout << "peakno: " << j << "  peakch: " << peaks[j] <<  "  peakcherr: " << peakerr[j] << "  peakamplitude: " << alphastemp[j] << "  peakamperr: " << alphastemperr[j] << endl;
      //g1->SetPoint(j, peaks[j], alphastemp[j]);
      g1err->SetPoint(j, peaks[j], alphastemp[j]);
      g1err->SetPointError(j, peakerr[j], alphastemperr[j]);
    }
    

    //g1->SetPoint(0, peaks[3], alphastemp[3]);
    //g1err->SetPoint(0, peaks[3], alphastemp[3]);
    //g1err->SetPointError(0, peakerr[3], alphastemperr[3]);

    //g1->SetPoint(0, peaks[4], alphastemp[4]);
    //g1err->SetPoint(0, peaks[4], alphastemp[4]);
    //g1err->SetPointError(0, peakerr[4], alphastemperr[4]);

    //g1->SetPoint(1, peaks[9], alphastemp[9]);
    //g1err->SetPoint(1, peaks[9], alphastemp[9]);
    //g1err->SetPointError(1, peakerr[9], alphastemperr[9]);

    //g1err->SetPoint(0,peaks[1],alphastemp[1]);
    //g1err->SetPointError(0, peakerr[1], alphastemperr[1]);
    //g1err->SetPoint(1,peaks[2],alphastemp[2]);
    //g1err->SetPointError(2, peakerr[2], alphastemperr[2]);
    //g1err->SetPoint(2,peaks[3],alphastemp[3]);
    //g1err->SetPointError(3, peakerr[3], alphastemperr[3]);
    //g1err->SetPoint(3,peaks[4],alphastemp[4]);
    //g1err->SetPointError(4, peakerr[4], alphastemperr[4]);


    g1err->SetMarkerStyle(20);
    g1err->SetMarkerSize(1.5);
    
    if (polDeg == 1 ) {
      //if (i == 4 ) { TF1 *f1 = new TF1("f1","pol1",100.,7675.); g1err->Fit("f1","R");}
      //if (i == 14 ) { TF1 *f1 = new TF1("f1","pol1",20.,2000.); g1err->Fit("f1","R");}
      //else {TF1 *f1 = new TF1("f1","pol1",100.,7675.); g1err->Fit("f1","R");}
      list->Add(g1err);
      TF1 *f1 = new TF1("f1","pol1",100.,16000.);
      g1err->Fit("f1","R+");
      
      TF1* myfit;
      myfit = (TF1*)g1err->GetFunction("f1");
      
      int board = i/32;
      int chip = i%32;

      fout << board << " " << chip << " " << myfit->GetParameter(1)/1000 << " " << myfit->GetParameter(0)/1000 << endl;
     
      
      chisqdarray[i] = myfit->GetChisquare()/myfit->GetNDF();
      cout << "Chisqd: " << myfit->GetChisquare() << " ndf: " << myfit->GetNDF() << " \nChiSqd/NDF = " << myfit->GetChisquare()/myfit->GetNDF() <<  " \nProbability:  " << myfit->GetProb() << endl;
      
      for (Int_t j=0; j<numPeaks; j++) {
        residualarray[i][j] = (( myfit->GetParameter(1) * peaks[j] + myfit->GetParameter(0) ) - alphastemp[j]);
      }
    }
    
    
    g2->SetPoint(i, peaks[i], chisqdarray[i]);
    
  }// end loop over strips
  

  cChi->cd();
  
  g2->Draw("AB");
  
  cRes->cd();

  
  for (Int_t i=0; i<numPeaks; i++) {
    for (Int_t j=0; j<numDetectorElements; j++) {
      hResiduals->Fill(j+1,peaks[i],residualarray[j][i]);
      hResiduals->GetXaxis()->SetTitle("Strip");
      hResiduals->GetYaxis()->SetTitle("Ch");
      hResiduals->Draw("colz");
    }
  }
  
  
  fout.close();
  fin1.close();
  //fin2.close();

  TFile * outfile = new TFile("linearfit.root", "RECREATE");
  outfile->cd();
  list->Write();
  outfile->Close();


}
