#include "stdheader.h"
    double minimum(double x, double y,double z)
                    {
                  
                   double smallest;
                   if((x < y)&&(x < z))
                  smallest = x;
                  else
                 {
                 if(y < z)
                 smallest = y;
                 else
                 smallest = z;
                 }
    return smallest;
    }
