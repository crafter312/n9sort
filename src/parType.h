#ifndef parType_H
#define parType_H

#include "constants.h"
#include "solution.h"

class parType {

public:

	DEFINE_ENUM_WITH_STRING_CONVERSIONS(detType, (Gobbi)(TexNeut)(S800))

	parType(int Z, int A, detType det = Gobbi);

	void zeroMask();
	void setMask();

	int Z, A;
	solution *Sol[6];
	int mult;
	bool mask[6];
	detType detector;

};

#endif
