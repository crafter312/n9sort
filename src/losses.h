#ifndef losses_
#define losses_
#include "loss2.h"
#include "SortConfig.h"

class CLosses
{
 private:
   CLoss2 ** loss;
   int Zmax;
 public:
   CLosses(int,SortConfig&);
   ~CLosses();
   double getEin(double,double,size_t,double);
   double getEout(double,double,size_t,double);

};

#endif
