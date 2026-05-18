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
	const float tol{0.000001}; // energy tolerance in MeV for adaptive step calculation
	
};
#endif
