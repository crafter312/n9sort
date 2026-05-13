#ifndef _telescope
#define _telescope

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "TMath.h"
#include <cmath>
#include "CRandom.h"
#include "elist.h"
#include "solution.h"
#include "pid.h"
#include "targthick.h"
#include "losses.h"
#include "calibrate.h"
#include "constants.h"
using namespace std;

/*
//structure for storing zlines
struct lines
{
  int n; //number of points
  float *x; //pointer to x array
  float *y; //pointer to y array
};
*/

class telescope
{
 public:
  telescope(bool S800); //bool for S800 as the 5th "telescope"
  ~telescope();
  void init(int);
  void SetTarget(double, float);
  void reset();
  void Reduce();

  int testingHitE();
  int simpleECsI();
  int multiHitECsI(); //+loop is private
  int multiHitdEECsI();
  void load(int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int);

  int getPID();
  int calcEloss();

  CLosses * Targlosses;
  CLosses * Allosses;
  float TargetThickness;

  calibrate * calCsi_d;
  calibrate * calCsi_t;
  calibrate * calCsi_Alpha;
  float light2energy(int,int,int,float);

  int id;
  float maxFront;
  float maxBack;
  int imaxFront;
  int imaxBack;

  int multFront;
  int multBack;
  int multDelta;
  int multCsI;

  elist Front;
  elist Back;
  elist CsI;

  elist tempFront;
  elist tempBack;

  solution Solution[10];
  int NSisolution;
  int Nsolution;

  //Want to save CsI data for later
  //Can look for gammas in coincidence with Si-CsI particle solutions
  solution CsISolution[10];
  int CsINsolution;


  pid * Pid;
  pid * EPid;
  pid * PidECsI[4];

  int simpleFrontBack();
  void position(int);

 private:
  calibrate* calCsI_Alpha;
  calibrate* calCsI_d;
  calibrate* calCsI_t;
  calibrate* calCsI_3He;


  int FrontLow[4];
  int FrontHigh[4];
  int BackLow[4];
  int BackHigh[4];

  //position
  float Xcenter; // center of detector in cm along x axis
  float Ycenter; // center of detector in cm along y axis
  float SiWidth;
  float SiFrame;
  CRandom *Ran;
  targthick *Tthick;

  //for nested loops
  int NestDim;
  void loop(int);
  int NestArray[50];
  int arrayB[50]; //matches Front to Back
  float deMin;
  int dstripMin;

  //for testing gobbi analyze() subroutine, set to desired telescope
  int idtest = 4;
  int prtswtch = 0;

  //Switch to turn pulse-height defect correction on
  bool PHDSw = false;
};
#endif
