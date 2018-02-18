#ifndef Plane_H
#define Plane_H
#include "mesh.h"
class Plane : public Mesh
{
    public:
    Vector normal;
    Vector p;
    Plane(float tx,float ty,float tz);
    Plane(float tx,float ty,float tz,Material theMaterial);
    Plane(float tx,float ty,float tz,Vector tp,Material theMaterial);
    bool checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos);
    Vector getNormal(Vector posx);
};
#endif
