#include "mesh.h"
#include "OBJ.h"
#include "vector.h"
#include <vector>
#include <string>
#include <fstream>
class OBJ;
OBJ::OBJ(char* filename)
{
    material=Material(Vector(0.7,0.7,0.7));
    std::vector<Vertex*> vertexVector;
    float tempVertex[3],tempNormal[3];
    int tempTriangle[3];
    std::vector<Triangle*> triangleVector;
    std::vector<Vector*> normalVector;
    std::string temp;
    std::string temp2="";
    std::ifstream infile;
    infile.open(filename);
    std::cout << "Reading file "<< filename <<std::endl;
    while(!infile.eof())
    {
        temp=infile.get();
		std::cout << temp << std::endl;
        if(temp=="v") // Vertex or VT
        {

                infile >> tempVertex[0];
                infile >> tempVertex[1];
                infile >> tempVertex[2];
                vertexVector.push_back(new Vertex(tempVertex[0],tempVertex[1],tempVertex[2]));
                getline(infile,temp);

        }
        else if(temp=="f")
        {

                infile >> tempTriangle[0];
                infile >> tempTriangle[1];
                infile >> tempTriangle[2];
                triangleVector.push_back(new Triangle(tempTriangle[0]-1,tempTriangle[1]-1,tempTriangle[2]-1));
                getline(infile,temp);

        }


    }
    infile.close();
    numVertices = vertexVector.size();
    vertices = new Vertex[numVertices];
    for(int i=0;i<numVertices;i++)
        vertices[i]=*vertexVector.at(i);

    numTriangles = triangleVector.size();
    triangles = new Triangle[numTriangles];
    for(int i=0;i<numTriangles;i++)
        triangles[i]=*triangleVector.at(i);

    std::cout << "Triangles: "<< numTriangles << std::endl;
    std::cout << "Vertices: "<<numVertices << std::endl;
   // calculateNormals();


}

void OBJ::calculateBounds()
{

    for(int i=0;i<numVertices;i++)
    {
        boundingPos.x=boundingPos.x+vertices[0].x;
        boundingPos.y=boundingPos.y+vertices[1].y;
        boundingPos.z=boundingPos.z+vertices[2].z;
    }
    boundingPos=boundingPos/numVertices;
    float length;
    Vector temp;
    for(int i =0;i<numVertices;i++)
    {
        temp=vertices[0];
        length=(temp-boundingPos).length();
            if(length>boundingRadius)
                boundingRadius=length;
    }



}
bool OBJ::checkIntersection(Vector &s,Vector &d, float &theDistance,Vector &thePos)
{

    Vector endMovement;
    bool hit=false;

  /*  Vector c=boundingPos;
    Vector v=(s-c);
    float wee=(v.dot(d))*(v.dot(d))-(v.x*v.x+v.y*v.y+v.z*v.z-boundingRadius*boundingRadius);

    if(wee<0) // No intersection with the boundinsphere!
    {
        return false;
    }
*/


    for(int i=0;i<numTriangles;i++)
    {

        // Check collision with the triangles plane
        Vector l0=s;
        Vector tempVector=vertices[triangles[i].p1];
        tempVector=tempVector-l0;
        float tempDot=tempVector.dot(triangles[i].normal);
        float res=tempDot/d.dot(triangles[i].normal);

        if(res>0.0001) // The triangles plane hit, check collision with triangle parametically
        {

            endMovement=d*res;
            Vector p=s+endMovement;

            Vector v0=vertices[triangles[i].p2]-vertices[triangles[i].p1];
            Vector v1=vertices[triangles[i].p3]-vertices[triangles[i].p1];
            Vector v2=p-vertices[triangles[i].p1];

            float dot00=v0.dot(v0);
            float dot01=v0.dot(v1);
            float dot02=v0.dot(v2);
            float dot11=v1.dot(v1);
            float dot12=v1.dot(v2);
            float invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
            float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

            // Check if point is in triangle
            if((u > 0.00001) && (v > 0.00001) && (u + v < 1))
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

Vector OBJ::getNormal(Vector posx)
{
    return triangles[lasthit].normal;
}

void OBJ::calculateNormals()
{
    for(int i=0;i<numTriangles;i++)
    {
        Vector temp=vertices[triangles[i].p3]-vertices[triangles[i].p2];
        Vector temp2=vertices[triangles[i].p1]-vertices[triangles[i].p2];
        triangles[i].normal=temp.cross(temp2);
        triangles[i].normal.normalize();
    }
}

void OBJ::transformVertices(Matrix4 m)
{
    for(int i=0;i<numVertices;i++)
    {
        vertices[i]=m.mult(vertices[i]);
    }


}
