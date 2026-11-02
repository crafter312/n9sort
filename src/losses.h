#ifndef losses_
#define losses_
#if defined(eloss_mode) && eloss_mode == 0
#include "loss3.h"
#elif defined(eloss_mode) && eloss_mode == 1
#include "loss4.h"
#else
#error "compile-time variable 'eloss_mode' is undefined or set to an unsupported value!"
#endif

#include <string>

class CLosses {
private:
#if defined(eloss_mode) && eloss_mode == 0
	CLoss3** loss;
#elif defined(eloss_mode) && eloss_mode == 1
	CLoss4** loss;
#else
#error "compile-time variable 'eloss_mode' is undefined or set to an unsupported value!"
#endif
   int Zmax;
   std::string material;
 public:
   CLosses(int,std::string,std::string);
   ~CLosses();
   double getEin(double,double,size_t,double);
   double getEout(double,double,size_t,double);

};

#endif
