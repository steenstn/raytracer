#ifndef Mesh_H
#define Mesh_H

#include "material.h"
#include "vector.h"
#include "matrix4.h"
class Mesh
{
    public:
    Material material;
    Mesh();
    virtual bool checkIntersection(Vector &s,Vector &d, double &theDistance,Vector &thePos);
    virtual Vector getNormal(Vector posx);
    int sort;
    int numConnections;
    int numTriangles;

    double radius;
    Vector position;

    virtual void transformVertices(Matrix4 m);
    virtual void calculateNormals();
    virtual void calculateBounds();
    virtual Vector getRandomPoint();
    virtual Vector getRefractionPoint(Vector s,Vector &d);
};

#endif
