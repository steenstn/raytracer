#ifndef cube_h
#define cube_h
#include "mesh.h"
#include "vertex.h"
#include "triangle.h"
#include "vector.h"
#include "matrix4.h"
class Cube : public Mesh
{
    public:
    Vertex* vertices;
    Triangle* triangles;
    int numTriangles,numVertices;
    int lasthit; // Triangleindex for the triangle that was last hit by a ray
    Cube();

    bool checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos);
    Vector getNormal(Vector posx);
    void transformVertices(Matrix4 m);
    void calculateNormals();
    Vector getRandomPoint();
    Vector getRefractionPoint(Vector s,Vector &d);
    bool pnpoly(unsigned char npol, float *xp, float *yp, float x, float y);


};
#endif
