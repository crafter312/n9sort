#ifndef _telescope
#define _telescope

/* Class dealing with a single Hira Si telescope
 * Modified by Henry Webb (h.s.webb@wustl.edu) May 2026
 * to look more pretty (among other small modifications).
 * CsI option is designed to be generalized, should work
 * for both New Gobbi with 4 CsI per telescope and
 * Gobbi 28 with 7 CsI per telescope
 */

//       .-.      _______                             .  '  *   .  . '
//      {}``; |==|_______D                              .  * * -+-  
//      / ('        /|\                             .    * .    '  *
//  (  /  |        / | \                                * .  ' .  . 
//   \(_)_]]      /  |  \                            *   *  .   .
//                                                     '   *

// Upstream view of Gobbi
//        ____
//       |    |____                  
//       | 4  |    |              
//      _|____| 1  |
//     |    |_|____|
//     | 3  |    |
//     |____| 2  |
//          |____|

#include <utility>
#include <vector>

#include <TRandom.h>

#include "elist.h"
#include "solution.h"
#include "pid.h"
#include "losses.h"
#include "SortConfig.h"

class telescope {

public:
	telescope(double, SortConfig&, bool csi=false);
	~telescope();
	void reset();
	void init(int, SortConfig&);
	void Reduce();
	int simpleFront();
	int multiHit();
	int testingHitE();
	int simpleECsI();
	int multiHitECsI();
	void SetTargetDistance(double);
	size_t getPID();
	int calcEloss();
	bool isCenter(size_t ifront, size_t iback);

	CLosses* losses;
	CLosses* Allosses;
	double TargetThickness;

	int id;
	
	int multFront;
	int multBack;
	int multDelta;
	int multCsI;

	elist Front;
	elist Back;
	elist Delta;
	elist CsI;

	solution Solution[20];
	solution tempSol;
	int Nsolution = 0;

	pid* Pid;
	vector<pid*> PidCsI;

	int simpleFrontBack();
	void position(int);
	void positionC(int);
	
	// Getters
	size_t GetNCsI() const { return NCsI; }

 private:
	int FrontLow[4];
	int FrontHigh[4];
	int BackLow[4];
	int BackHigh[4];

	// Position
	double Xcenter; // center of detector in cm along x axis
	double Ycenter; // center of detector in cm along y axis
	double SiWidth;
	double SiFrame;
	double holeSize;
	TRandom* Ran;
	
	// FB silicon strip extents for CsI crystal matching (min, max, inclusive)
	std::vector<std::pair<size_t, size_t>> CsIFextents;
	std::vector<std::pair<size_t, size_t>> CsIBextents;
	std::vector<size_t> CsIFmids;
	std::vector<size_t> CsIBmids;
	size_t NCsI{0};

	//for nested loops
	size_t NestDim;
	void loopDEE(int);
	void loopE(int);
	int NestArray[50];
	int arrayD[50];
	int arrayB[50];
	double deMin;
	int dstripMin;
	size_t NSisolution;
	
	double alThick; // thickness of aluminum absorbers, in mg/cm^2
	bool hasCsI;

};
#endif
