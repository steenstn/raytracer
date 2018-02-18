#ifndef Vector_H
#define Vector_H
#include <iostream>
#include <math.h>
#include "vertex.h"
class Vector
{
    public:
    float x,y,z;

    Vector();
    Vector(const Vector& v);
    Vector(float theX,float theY,float theZ);
    Vector(const Vertex& v);
    float dot(Vector b);
    Vector cross(Vector v);
    void normalize(void);

    float length(void);

    Vector operator+ (const Vector& v);

    Vector operator* (const Vector& v);
    Vector operator- (const Vector& v);

    Vector operator* (const float &d);
    Vector operator/ (const float &d);
    Vector operator- (const float &d);
    float sum(void);
    bool operator== (Vector& d);

    const Vector& operator= (const Vector& v);

    friend std::ostream& operator<<(std::ostream& o,Vector& v);

};
#endif
