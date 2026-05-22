#ifndef loss2_H
#define loss2_H

#include <string>
#include <utility>
#include <vector>

/**
 * !\brief energy loss of particles in an absorber
 */

class CLoss3 {

public:
	CLoss3(std::string);
	~CLoss3();
	
	float getEout(float, float, float);
	float getEin(float, float, float);
	
private:
	std::pair<float, float> getAbsSlopeDedx(float, float);
	
	size_t N;
	std::vector<float> Ein;
	std::vector<float> dedx;
	std::vector<float> slope; // size N - 1
	
	float Emax;
	
	// Energy tolerance in MeV for adaptive step calculation.
	// To set this, in `loss3.h`, uncomment cout statements
	// to print out the dthick values for each loop iteration
	// and then abort at the end. This will allow you to
	// inspect the dthick evolution for one energy loss or
	// gain function call and determine how it is affected
	// by the tolerance value below.
	const float tol{0.001};
	
};
#endif
