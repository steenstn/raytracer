#ifndef Plane_H
#define Plane_H
#include "mesh.h"
class Plane : public Mesh
{
    public:
    Vector normal;
    Vector p;
    Plane(double tx,double ty,double tz);
    Plane(double tx,double ty,double tz,Material theMaterial);
    Plane(double tx,double ty,double tz,Vector tp,Material theMaterial);
    bool checkIntersection(Vector &s,Vector &d, double &theDistance,Vector &thePos);
    Vector getNormal(Vector posx);
};
#endif
