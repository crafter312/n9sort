#ifndef _telescope
#define _telescope

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "TMath.h"
#include <cmath>
#include "TRandom.h"
#include "elist.h"
#include "solution.h"
#include "pid.h"
#include "losses.h"
#include "SortConfig.h"
using namespace std;

/*
//structure for storing zlines
struct lines
{
  int n; //number of points
  double *x; //pointer to x array
  double *y; //pointer to y array
};
*/

class telescope
{
 public:
  telescope(double, SortConfig&);
  ~telescope();
  void reset();
  void init(int, SortConfig&);
  void Reduce();
  int simpleFront();
  int multiHit();
  void SetTargetDistance(double);
  size_t getPID();
  int calcEloss();

  CLosses * losses;
  double TargetThickness;

  int id;
  double maxFront;
  double maxBack;
  double maxDelta;
  int imaxFront;
  int imaxBack;
  int imaxDelta;
  int multFront;
  int multBack;
  int multDelta;

  elist Front;
  elist Back;
  elist Delta;

  solution Solution[20];
  int Nsolution = 0;

  pid * Pid;

  int simpleFrontBack();
  void position(int);
  void positionC(int);

 private:
  int FrontLow[4];
  int FrontHigh[4];
  int BackLow[4];
  int BackHigh[4];

  //position
  double Xcenter; // center of detector in cm along x axis
  double Ycenter; // center of detector in cm along y axis
  double SiWidth;
  TRandom *Ran;

  //for nested loops
  int NestDim;
  void loop(int);
  int NestArray[50];
  int arrayD[50];
  int arrayB[50];
  double deMin;
  int dstripMin;

};
#endif
