#ifndef triangle_h
#define triangle_h
#include "vector.h"
class Triangle
{
    public:
    int p1,p2,p3;
    Vector normal;
    Triangle();
    Triangle(int a,int b,int c);
};
#endif
