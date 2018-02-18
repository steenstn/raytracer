#ifndef Sphere_H
#define Sphere_H

#include "mesh.h"
#include "material.h"
class Sphere : public Mesh
{
    public:
    float radius;
    Sphere();
    Vector emittance;
    Sphere(float theX,float theY,float theZ, float theRadius);
    Sphere(float theX,float theY,float theZ, float theRadius,Material theMaterial);
	Sphere(float theX,float theY,float theZ, float theRadius,Material theMaterial, bool isReflective);
    bool checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos);
    Vector getNormal(Vector posx);
    Vector getRandomPoint();
    Vector getRefractionPoint(Vector s,Vector &d);
};

#endif
