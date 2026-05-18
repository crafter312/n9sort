#ifndef losses_
#define losses_
#include "loss2.h"

#include <string>

class CLosses
{
 private:
   CLoss2 ** loss;
   int Zmax;
   std::string material;
 public:
   CLosses(int,std::string,std::string);
   ~CLosses();
   double getEin(double,double,size_t,double);
   double getEout(double,double,size_t,double);

};

#endif
