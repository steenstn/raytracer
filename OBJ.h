#ifndef OBJ
#define OBJ_h
#include "vertex.h"
#include "triangle.h"
#include "vector.h"
#include "matrix4.h"
class OBJ : public Mesh
{
    public:
    Vertex* vertices;
    Triangle* triangles;
    int numTriangles,numVertices;
    int lasthit; // Triangleindex for the triangle that was last hit by a ray
    OBJ(char* filename);
    Vector boundingPos;
    float boundingRadius;
    void calculateBounds();
    bool checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos);
    Vector getNormal(Vector posx);
    void transformVertices(Matrix4 m);
    void calculateNormals();

};
#endif
