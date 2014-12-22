#ifndef Vector_H
#define Vector_H
#include <iostream>
#include <math.h>
#include "vertex.h"
class Vector
{
    public:
    double x,y,z;

    Vector();
    Vector(const Vector& v);
    Vector(double theX,double theY,double theZ);
    Vector(const Vertex& v);
    double dot(Vector b);
    Vector cross(Vector v);
    void normalize(void);

    double length(void);

    Vector operator+ (const Vector& v);

    Vector operator* (const Vector& v);
    Vector operator- (const Vector& v);

    Vector operator* (const double &d);
    Vector operator/ (const double &d);
    Vector operator- (const double &d);
    double sum(void);
    bool operator== (Vector& d);

    const Vector& operator= (const Vector& v);

    friend std::ostream& operator<<(std::ostream& o,Vector& v);

};
#endif
