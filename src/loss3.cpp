#include "loss3.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "constants.h"

using namespace std;

/**
 * constructor
\param filename is name of file containing energy loss tables of a particulat particle
*/
CLoss3::CLoss3(string filename) {
  //cout << "opening loss file: " << filename.c_str() << endl;
  ifstream File(filename.c_str());
  if (File.is_open() != 1)
  {
    cout << " could not open loss file " << filename;
    return;
  }

  char line[100];
  File.getline(line,100);
  //cout << line << endl;

  File >> N;
  if (N <= 1) throw invalid_argument(string(BOLDRED) + string("[ERROR]: N from ") + filename + string(" cannot be one or zero") + string(RESET));
  Ein.resize(N);
  dedx.resize(N);
  slope.resize(N - 1);

  File >> Ein[0] >> dedx[0];
  for (size_t i = 1; i < N; i++) {
    File >> Ein[i] >> dedx[i];
    slope[i - 1] = (dedx[i] - dedx[i - 1]) / (Ein[i] - Ein[i - 1]);
    //cout << "index " << i - 1 << " slope " << slope[i - 1] << endl;
  }
  
  if (!is_sorted(Ein.begin(), Ein.end()))
    throw invalid_argument(string(BOLDRED) + string("[ERROR]: loss file ") + filename + string(" is not energy ordered") + string(RESET));

  Emax = Ein[N-1];
}

//****************************************************************
  /**
   * destructor
   */
CLoss3::~CLoss3() {}

//*****************************************************************
  /*
   * returns a pair of first the absolute value of the slope of the
   * DeDx table between the i and i + 1 bins and second the value
   * of DeDx interpolated from table between the i and i + 1 bins,
   * where i is the largest table element <= the provided energy
   \param energy is energy of particle in MeV
   \param A is the atomic mass of the nucleus in question
   */
pair<float, float> CLoss3::getAbsSlopeDedx(float energy, float A) {
  float epa = energy / A;
  
  // Throw error if energy is outside bounds of loss table
  if (epa < Ein.front() || epa > Ein.back()) {
    stringstream ss;
    ss << BOLDRED << "[ERROR] Closs2 input energy " << epa << " per A outside bounds of energy loss table of min "
       << Ein[0] << " and max " << Ein[N - 1] << " for A=" << A << RESET << endl;
    throw runtime_error(ss.str());
  }
  
  // Binary search to find correct energy bin
  auto it = lower_bound(Ein.begin(), Ein.end(), energy);
  size_t index = distance(Ein.begin(), it);
  if (*it != energy) index--;

  // linear interpolation
  float s = slope[index];
  float de = (s * (epa - Ein[index])) + dedx[index];
  //cout << "index " << index << endl;
  return make_pair(s * A, de);
}
//********************************************************************
  /**
   * returns the residual energy of particle after passage through absorber
\param energy is initial energy of particle in MeV
\param thick is the thickness of the absorber in mg/cm2
  */
float CLoss3::getEout(float energy, float thick, float A) {
  pair<float, float> p;
  float testStep, dthick, thickness, de;
  float Eout = energy;
  for(;;) {
    p = getAbsSlopeDedx(Eout, A);
    testStep = sqrt(2. * tol / abs(p.first * p.second)); // TODO: check that this is the right formula for adaptive step size
    dthick = max(testStep, 0.1f);
    thickness = min(thick, dthick);
    de = p.second;
    Eout -= de * thickness;
    if (thickness == thick) break;
    thick -= dthick;
  }
  return Eout;
}
//********************************************************************
  /**
   * Determined the inital energy of a particle before entering an 
   * absorber given its residueal
\param energy is the residual energy of the particle
\param thick is the thickness of absorber through which the particle passed.
  */
float CLoss3::getEin(float energy, float thick, float A) {
  pair<float, float> p;
  float testStep, dthick, thickness, de;
  float Einput = energy;
  //stringstream ss;
  for(;;) {
    //ss << "thick " << thick << " dthick " << dthick << " Einput " << Einput << endl;
    p = getAbsSlopeDedx(Einput, A);
    testStep = sqrt(2. * tol / abs(p.first * p.second)); // TODO: check that this is the right formula for adaptive step size
    dthick = max(testStep, 0.1f);
    thickness = min(thick, dthick);
    de = p.second;
    Einput += de * thickness;
    if (thickness == thick) break;
    thick -= dthick;
  }
  return Einput;
}


