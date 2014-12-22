#include "plane.h"
class Plane;

Plane::Plane(double tx,double ty,double tz)
{
    normal.x=tx;
    normal.y=ty;
    normal.z=tz;
}

Plane::Plane(double tx,double ty,double tz,Material theMaterial)
{
    normal.x=tx;
    normal.y=ty;
    normal.z=tz;
    material=theMaterial;
}

Plane::Plane(double tx,double ty,double tz,Vector tp,Material theMaterial)
{
    normal.x=tx;
    normal.y=ty;
    normal.z=tz;
    p=tp;
    material=theMaterial;
}

bool Plane::checkIntersection(Vector &s,Vector &d, double &theDistance,Vector &thePos)
{
    Vector endMovement;
    double t[1];
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
