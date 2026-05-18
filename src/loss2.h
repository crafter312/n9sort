#ifndef loss2_H
#define loss2_H

#include <string>
#include <utility>
#include <vector>

/**
 * !\brief energy loss of particles in an absorber
 */

class CLoss2 {

public:
	CLoss2(std::string);
	~CLoss2();
	
	float getEout(float, float, float);
	float getEin(float, float, float);
	
private:
	std::pair<float, float> getAbsSlopeDedx(float, float);
	
	size_t N;
	std::vector<float> Ein;
	std::vector<float> dedx;
	std::vector<float> slope; // size N - 1
	
	float Emax;
	float tol{0.00001}; // energy tolerance in MeV for adaptive step calculation
	
};
#endif
