#include "mesh.h"
class Mesh;

Mesh::Mesh()
{
    sort=-1;
}

bool  Mesh::checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos)
{
    return false;
}
Vector Mesh::getNormal(Vector posx)
{
	return Vector();
}

void Mesh::transformVertices(Matrix4 m)
{
    std::cout << "Mesh::transformVertices!";

}
void Mesh::calculateNormals()
{}
void Mesh::calculateBounds()
{}

Vector Mesh::getRandomPoint()
{
	return Vector();
}

Vector Mesh::getRefractionPoint(Vector s, Vector &d)
{
	return Vector();
}


