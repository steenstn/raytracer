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
    virtual bool checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos);
    virtual Vector getNormal(Vector posx);
    int sort;
    int numConnections;
    int numTriangles;

    float radius;
    Vector position;

    virtual void transformVertices(Matrix4 m);
    virtual void calculateNormals();
    virtual void calculateBounds();
    virtual Vector getRandomPoint();
    virtual Vector getRefractionPoint(Vector s,Vector &d);
};

#endif
