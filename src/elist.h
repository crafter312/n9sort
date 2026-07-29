#ifndef elist_
#define elist_

#include "histo.h"

using namespace std;

struct order {
	double energy;
	double energyR;    // high gain channels
	double energyRlow; // low channels
	double energylow;
	double energyMax;
	int strip;
	int neighbours;    // I smell australian here
	double time;
	double qdc;        // used for PSD with CsIs
	bool qdcflag;      // denotes CsI events with qdc present
	bool CsIFlag;      // denotes CsI events? Set but never used
};

int const nnn = 60;

/**
 * !\brief Energy ordered list
 *
 * This class creates an energy ordered list of the strips
 * read out from a strip detector, keeping track of the strip 
 * numbers that fired.
 */

class elist {

public:
	int Nstore = 0; //number stored in list
	order Order[nnn];
	int mult;

  histo& Histo;

  elist(histo&);
	void Add(int, double, double, int, int, double, double, int);
	void Add(int, double, int, int, double);
	void Add(int, double, int, int);
	void Remove(int);
	int Reduce(const char*);
	void reset();
	void Neighbours(int);
	void Threshold(double);
	double threshold0;
};
#endif
