#include "sphere.h"
class Sphere;


Sphere::Sphere()
{
}

Sphere::Sphere(float theX,float theY,float theZ, float theRadius)
{
    sort=0;
    position=Vector(theX,theY,theZ);
    radius=theRadius;

}
Sphere::Sphere(float theX,float theY,float theZ, float theRadius,Material theMaterial)
{
    sort=0;
    numConnections=0;
    position=Vector(theX,theY,theZ);
    radius=theRadius;
    material=theMaterial;
}

Sphere::Sphere(float theX,float theY,float theZ, float theRadius,Material theMaterial, bool isReflective)
{
    sort=0;
    numConnections=0;
    position=Vector(theX,theY,theZ);
    radius=theRadius;
    material=theMaterial;
	material.reflective = isReflective;
}

Vector Sphere::getNormal(Vector posx)
{
    Vector n=((posx-position)/(posx-position).length());
    n.normalize();
    return n;
}

bool Sphere::checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos)
{
   // Check for intersection
    Vector c=position;
    Vector v=(s-c);
    float wee=(v.dot(d))*(v.dot(d))-(v.x*v.x+v.y*v.y+v.z*v.z-radius*radius);
    float t[2];
    Vector endMovement;
    if(wee>0) // Intersection!
    {

        t[0]=v.dot(d)*-1+sqrt(wee); // Find the two intersection points
        t[1]=v.dot(d)*-1-sqrt(wee);

      //  if(t[0]>0.000001 && t[1]>0.000001)
        {
            if(t[0]<t[1] && t[0]>0.000001) // Find closest
            {
                endMovement=d*t[0];
                    if(endMovement.length()<theDistance) // Find closest of all the spheres
                    {
                        theDistance=endMovement.length();
                        thePos=s+endMovement;
                        return true;
                    }
            }
            else if(t[1]<t[0] && t[1]>0.000001)
            {
                endMovement=d*t[1];
                    if(endMovement.length()<theDistance) // Find closest of all the spheres
                    {
                        theDistance=endMovement.length();
                        thePos=s+endMovement;
                        return true;
                    }
            }
//            if(endMovement.length()<theDistance) // Find closest of all the spheres
//            {
//                theDistance=endMovement.length();
//                thePos=s+endMovement;
//                return true;
//            }

        }


    }
    return false;
}

Vector Sphere::getRandomPoint()
{
    Vector randVector = Vector((float)rand()/RAND_MAX-(float)rand()/RAND_MAX,(float)rand()/RAND_MAX-(float)rand()/RAND_MAX,(float)rand()/RAND_MAX-(float)rand()/RAND_MAX);
    randVector.normalize();
    randVector = randVector*radius;
    return randVector;

}
Vector Sphere::getRefractionPoint(Vector s,Vector &d)
{
    d.normalize();
    float n1 = 1.0;
    float n2 = 1.4;

    Vector theNormal=getNormal(s);


    float n = n1/n2;
    float cosI = theNormal.dot(d);
    float sinT2 = n * n * (1.0 - cosI*cosI);

    Vector refracted = d*n + theNormal*(n*cosI-sqrt(1-sinT2));
    refracted.normalize();

    d=refracted;



     // Check for intersections
    Vector c=position;
    Vector v=(s-c);
    float wee=(v.dot(d))*(v.dot(d))-(v.x*v.x+v.y*v.y+v.z*v.z-radius*radius);
    float t[2];
    Vector endMovement,thePos;


    t[0]=v.dot(d)*-1+sqrt(wee); // Find the two intersection points
    t[1]=v.dot(d)*-1-sqrt(wee);


    if(t[0]>t[1]) // Find closest
        endMovement=d*t[0];
    else
        endMovement=d*t[1];

    thePos=s+endMovement;


    // Get out of the sphere

    n = n2/n1;



    theNormal=getNormal(thePos)*-1;
    cosI = theNormal.dot(d);
    sinT2 = n * n * (1.0 - cosI*cosI);

    refracted = d*n + theNormal*(n*cosI-sqrt(1.0-sinT2));
    refracted.normalize();
    d=refracted;
    if(sinT2 > 1)
        return Vector(0,0,0);
    else
        return thePos;


}
