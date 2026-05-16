#ifndef parType_H
#define parType_H

#define MAX_SOLS 6

#include "constants.h"
#include "solution.h"

class parType {

public:

	//DEFINE_ENUM_WITH_STRING_CONVERSIONS(detType, (Gobbi)(TexNeut)(S800))
	enum detType {Gobbi,TexNeut,S800};

	parType(int Z, int A, detType det = Gobbi);

	void zeroMask();
	void setMask();

	int Z, A;
	solution *Sol[MAX_SOLS];
	int mult{MAX_SOLS};
	bool mask[MAX_SOLS];
	detType detector;

};

#endif
