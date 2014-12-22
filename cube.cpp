#include "material.h"
#include "cube.h"
#include "mesh.h"
#include <iostream>
class Cube;

Cube::Cube()
{
    material=Material(Vector(0.4,1,0));
    // Kub!
    vertices = new Vertex[8];
    //Botten
    vertices[0]=Vertex(0,0,0);
    vertices[1]=Vertex(1,0,0);
    vertices[2]=Vertex(0,0,1);
    vertices[3]=Vertex(1,0,1);
    //Bakre väggen
    vertices[4]=Vertex(0,1,0);
    vertices[5]=Vertex(1,1,0);
    vertices[6]=Vertex(1,1,1);
    vertices[7]=Vertex(0,1,1);

    numVertices=8;

    triangles= new Triangle[12];
    triangles[0]=Triangle(0,1,2);
    triangles[0].normal=Vector(0,-1,0);

    triangles[1]=Triangle(1,3,2);
    triangles[1].normal=Vector(0,-1,0);

    triangles[2]=Triangle(0,4,1);
    triangles[2].normal=Vector(0,0,-1);

    triangles[3]=Triangle(1,4,5);
    triangles[3].normal=Vector(0,0,-1);

    triangles[4]=Triangle(2,4,0);
    triangles[4].normal=Vector(-1,0,0);

    triangles[5]=Triangle(2,7,4);
    triangles[5].normal=Vector(-1,0,0);

    triangles[6]=Triangle(1,5,3);
    triangles[6].normal=Vector(1,0,0);

    triangles[7]=Triangle(6,3,5);
    triangles[7].normal=Vector(1,0,0);

    triangles[8]=Triangle(6,5,4);
    triangles[8].normal=Vector(0,1,0);

    triangles[9]=Triangle(7,6,4);
    triangles[9].normal=Vector(0,1,0);

    triangles[10]=Triangle(7,3,6);
    triangles[10].normal=Vector(0,0,1);

    triangles[11]=Triangle(7,2,3);
    triangles[11].normal=Vector(0,0,1);

    numTriangles=12;
// Slutkub



}
bool Cube::pnpoly(unsigned char npol, double *xp, double *yp, double x, double y)
{
	int i, j;
	bool c = false;
	for (i = 0, j = npol-1; i < npol; j = i++)
	{
		if ((((yp[i] <= y) && (y < yp[j])) ||
			((yp[j] <= y) && (y < yp[i]))) &&
				(x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
				c = !c;
	}
	return c;
}
bool Cube::checkIntersection(Vector &s,Vector &d, double &theDistance,Vector &thePos)
{

    Vector endMovement;
    bool hit=false;
    for(int i=0;i<numTriangles;i++)
    {

        // Check collision with the triangles plane
        Vector l0=s;
        Vector tempVector=vertices[triangles[i].p1];
        tempVector=tempVector-l0;
        double tempDot=tempVector.dot(triangles[i].normal);
        double res=tempDot/d.dot(triangles[i].normal);

        if(res>0.0000001) // The triangles plane hit, check collision with triangle parametically
        {
            endMovement=d*res;
            Vector p=s+endMovement;

            Vector v0=vertices[triangles[i].p2]-vertices[triangles[i].p1];
            Vector v1=vertices[triangles[i].p3]-vertices[triangles[i].p1];
            Vector v2=p-vertices[triangles[i].p1];

            double dot00=v0.dot(v0);
            double dot01=v0.dot(v1);
            double dot02=v0.dot(v2);
            double dot11=v1.dot(v1);
            double dot12=v1.dot(v2);
            double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
            double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

            // Check if point is in triangle
            if((u > -0.00001) && (v > -0.00001) && (u + v < 1.0001))
            {

                if(endMovement.length()<theDistance)
                {
                    hit=true;
                    lasthit=i;

                    thePos=s+endMovement;
                    theDistance=endMovement.length();
                }
            }
        }
    }
    if(hit)
    {
        return true;
    }
    else
        return false;
}

Vector Cube::getNormal(Vector posx)
{
    Vector n = triangles[lasthit].normal;
    n.normalize();
    return n;
}

Vector Cube::getRandomPoint()
{
    return Vector(vertices[0]);
}

void Cube::calculateNormals()
{
    for(int i=0;i<numTriangles;i++)
    {
        Vector temp=vertices[triangles[i].p3]-vertices[triangles[i].p2];
        Vector temp2=vertices[triangles[i].p1]-vertices[triangles[i].p2];
        triangles[i].normal=temp.cross(temp2);
        triangles[i].normal.normalize();
    }
}

void Cube::transformVertices(Matrix4 m)
{
    for(int i=0;i<numVertices;i++)
    {
        vertices[i]=m.mult(vertices[i]);
    }


}

Vector Cube::getRefractionPoint(Vector s,Vector &d)
{
    d.normalize();
//
    double n1=1.0;
    double n2=1.4;

    Vector theNormal = getNormal(s);


    double n = n1/n2;
    double cosI = theNormal.dot(d);
    double sinT2 = n * n * (1.0 - cosI*cosI);

    Vector refracted = d*n + theNormal*(n*cosI-sqrt(1-sinT2));
    refracted.normalize();

    d=refracted;

    Vector thePos;
    double distance=9999999;
    //s=s+d*0.0001;
    checkIntersection(s,d,distance,thePos);

    n = n2/n1;



    theNormal=getNormal(thePos)*-1;
    cosI = theNormal.dot(d);
    sinT2 = n * n * (1.0 - cosI*cosI);

    refracted = d*n + theNormal*(n*cosI-sqrt(1.0-sinT2));
    refracted.normalize();
    d=refracted;
   // thePos=thePos+d*0.0001;
    if(sinT2 > 1)
        return Vector(0,0,0);
    else
        return thePos;


}
