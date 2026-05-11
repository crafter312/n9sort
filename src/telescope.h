#ifndef _telescope
#define _telescope

/* Class dealing with a single Hira Si telescope
 * Modified by Henry Webb (h.s.webb@wustl.edu) May 2026
 * to look more pretty (among other small modifications)
 */

#include <TRandom.h>
#include <TMath.h>

#include "elist.h"
#include "losses.h"
#include "pid.h"
#include "solution.h"
#include "SortConfig.h"

class telescope {

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

	CLosses* losses;
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
	int Nsolution{0};

	pid* Pid;

	int simpleFrontBack();
	void position(int);
	void positionC(int);

private:
	int FrontLow[4];
	int FrontHigh[4];
	int BackLow[4];
	int BackHigh[4];

	// Position
	double Xcenter; // center of detector in cm along x axis
	double Ycenter; // center of detector in cm along y axis
	double SiWidth;
	TRandom *Ran;

	// For nested loops
	int NestDim;
	void loop(int);
	int NestArray[50];
	int arrayD[50];
	int arrayB[50];
	double deMin;
	int dstripMin;

};

#endif
