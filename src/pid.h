#ifndef pid_
#define pid_

/* Class to detemine PID (Particle IDentification) from E-DE map via
 * banana gates stored for each particle type.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu), September 2025.
 * Mass values no longer hardcoded, now retrieved from std::unordered_map. 
 */

#include <string>

#include "constants.h"
#include "ZApar.h"
#include "SortConfig.h"

class pid {

public:
	pid(std::string file,SortConfig&); 
	~pid();

	bool getPID(float x, float y);
	double getMass(int iZ,int iA);

	ZApar** par; // individual banana gates
	int nlines;  // number of banana gated stored 	
	size_t Z;    // Z of particle in gate
	size_t A;    // A of particle in gate
	double mass; // mass of particle in amu

};

#endif
