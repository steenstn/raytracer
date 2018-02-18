#include "plane.h"
class Plane;

Plane::Plane(float tx,float ty,float tz)
{
    normal.x=tx;
    normal.y=ty;
    normal.z=tz;
}

Plane::Plane(float tx,float ty,float tz,Material theMaterial)
{
    normal.x=tx;
    normal.y=ty;
    normal.z=tz;
    material=theMaterial;
}

Plane::Plane(float tx,float ty,float tz,Vector tp,Material theMaterial)
{
    normal.x=tx;
    normal.y=ty;
    normal.z=tz;
    p=tp;
    material=theMaterial;
}

bool Plane::checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos)
{
    Vector endMovement;
    float t[1];
    t[0]=((p-s).dot(normal))/d.dot(normal);
    if(t[0]>0.000001)
    {

        endMovement=d*t[0];
        if(endMovement.length()<theDistance) // Find closest of all the spheres
        {
            theDistance=endMovement.length();
            thePos=s+endMovement;
            return true;
        }
    }
    return false;
}

Vector Plane::getNormal(Vector posx)
{
    normal.normalize();
    return normal;
}
