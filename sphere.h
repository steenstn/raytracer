#ifndef Sphere_H
#define Sphere_H

#include "mesh.h"
#include "material.h"
class Sphere : public Mesh
{
    public:
    double radius;
    Sphere();
    Vector emittance;
    Sphere(double theX,double theY,double theZ, double theRadius);
    Sphere(double theX,double theY,double theZ, double theRadius,Material theMaterial);
	Sphere(double theX,double theY,double theZ, double theRadius,Material theMaterial, bool isReflective);
    bool checkIntersection(Vector &s,Vector &d, double &theDistance,Vector &thePos);
    Vector getNormal(Vector posx);
    Vector getRandomPoint();
    Vector getRefractionPoint(Vector s,Vector &d);
};

#endif
