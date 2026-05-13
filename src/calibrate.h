#ifndef calibrate_
#define calibrate_
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

/**
 * storge of calibration coefficients
 */
struct coeff
{
  double slope; //!< slope for calibration
  double intercept; //!< intercept for calibration
  double a2;  //!< quadratic coeff if needed
  double a3; //!< cubic coeff if needed
};

class calibrate
{
 public:
  calibrate(int Ntele,int Nstrip,string file,int order,bool weave);
  ~calibrate();
  double getEnergy(int itele,int istrip,double channel);
  double getTime(int itele,int istrip,double channel);
  double reverseCal(int itele, int istrip, double energy);
  int order;
  int Nstrip;  //!< number of strips
  int Ntele;   //!<number of telescopes
  coeff ** Coeff;  //!< array with calibration coefficients for each strip

};
#endif
