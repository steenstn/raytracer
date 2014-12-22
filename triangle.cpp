#include "triangle.h"

class Triangle;

Triangle::Triangle()
{
    p1=p2=p3=-1;
}

Triangle::Triangle(int a,int b,int c)
{
    p1=a;
    p2=b;
    p3=c;
}

