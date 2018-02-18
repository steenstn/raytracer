/*
Monte Carlo path tracer

Här är jag!
*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <typeinfo>
#include <omp.h>

#include "vector.h"
#include "sphere.h"
#include "mesh.h"
#include "material.h"
#include "plane.h"
#include "vertex.h"
#include "triangle.h"
#include "cube.h"
#include "matrix4.h"
#include "OBJ.h"

const int SCREENWIDTH = 1024;
const int SCREENHEIGHT = 768;

const double aaFactor=1; // Antialias factor

int numFrames=1;
int numRays=10; // Number of rays/iteration
float DoF=0.9;
const int WIDTH  = SCREENWIDTH*aaFactor;
const int HEIGHT = SCREENHEIGHT*aaFactor;
using namespace std;
unsigned char img[WIDTH*HEIGHT*3]; // Slutbilden sparad som en lång jävla sträng


short picture[SCREENWIDTH][SCREENHEIGHT][3];
double pictureAA[WIDTH][HEIGHT][3];

int bounces=0; // Number of reflected rays

int maxBounces=50; // Maximum number of bounces allowed


vector<Mesh*> theMeshes; // All of the meshes


void savebmp(const char *filename, int w,int h);
void superSample(float numPasses);
Vector shootRay(Vector s,Vector d,int index);
Vector explicitRay(Vector s,int index);
Vector shootRefractedRay(Vector s,Vector d,int index,float n1);

float step(float a,float x) {
    return (float)(x>=a);
}
float clamp(float x,float a,float b) { // Clamp x between a and b(THE CLAMPS!)
    return (x < a ? a: (x > b ? b : x));
}
float random(void) {
    return (float)rand()/(float)RAND_MAX;
}

Vector randomMixed(Vector mix) {
	float r = random();
	float g = random();
	float b = random();

	r = (r + mix.x) / 2;
	g = (g + mix.y) / 2;
	b = (b + mix.z) / 2;

	return Vector(r,g,b);
}

Vector mixWhite(Vector v) {
	return Vector((v.x + 1) / 2, (v.y + 1) / 2, (v.z + 1) / 2);
}

int main(void) {
    srand((unsigned)time(0));

/// Main scene
   // theMeshes.push_back(new Sphere(0,-4,3,1.0,Material(Vector(8,8,8),Vector(0.8,0.45,0.8))));


    //theMeshes.push_back(new Sphere(-1,-.3,0,1.8,Material(Vector(0,0,0),Vector(0.8,0.45,0.8))));
   // theMeshes.push_back(new Sphere(0,0,0,1.8,Material(Vector(2,2,2),Vector(0.8,0.45,0.8))));
	/*

    theMeshes.push_back(new OBJ("C:\\cpp\\Monte Carlo\\buddha.obj"));

    theMeshes.push_back(new Sphere(1.6,0.0,4,0.6,Material(Vector(0,0,0),Vector(0.9,0.9,0.9))));


	
    theMeshes.push_back(new Sphere(0.5,0.5,4,0.6,Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));


    int j=1;

        theMeshes.push_back(new Sphere(0,0.5,5,0.5,Material(Vector(0,0,0),Vector(0.8,0.8,0.8 ))));
        theMeshes.at(j++)->material.transparent=true;


    theMeshes.push_back(new Plane(0,-1,0,Vector(0,1,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
    theMeshes.push_back(new Plane(0,1,0,Vector(0,-4,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(0,0,1,Vector(0,0,-5),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(0,0,-1,Vector(0,0,16),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(1,0,0,Vector(-3,0,0),Material(Vector(0,0,0),Vector(0.75,0.25,0.25))));
    theMeshes.push_back(new Plane(-1,0,0,Vector(3,0,0),Material(Vector(0,0,0),Vector(0.25,0.75,0.25))));


    Matrix4 m,m2;
    m=m.rotateX(0.2);
    theMeshes.at(2)->transformVertices(m);

    m=m.scale(2,0.2,4);
    theMeshes.at(0)->transformVertices(m);
    m=m.translate(-1,-4,0);
    theMeshes.at(0)->transformVertices(m);

     m=m.translate(-1,0.2,4);
    theMeshes.at(2)->transformVertices(m);

    theMeshes.at(2)->calculateNormals();
    theMeshes.at(2)->calculateBounds();
	*/
/// End main scene

/// Korridor
	/*
    for(int i=0;i<6;i++)
        theMeshes.push_back(new Sphere(3,-3,8-10*i,0.7,Material(Vector(8,8,8),Vector(0.8,0.45,0.8))));

    for(int i=0;i<10;i++)
    {
        theMeshes.push_back(new Sphere(-0.5+2*sin(i),0.1,7-4*i,0.6,Material(Vector(0,0,0),Vector((double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX))));
        theMeshes.at(6+i)->material.transparent=true;
    }
    theMeshes.push_back(new Plane(0,-1,0,Vector(0,1,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
    theMeshes.push_back(new Plane(0,1,0,Vector(0,-4,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(0,0,1,Vector(0,0,-50),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(0,0,-1,Vector(0,0,16),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(1,0,0,Vector(-3,0,0),Material(Vector(0,0,0),Vector(0.75,0.75,0.75))));
    theMeshes.at(20)->material.reflective=true;
    theMeshes.push_back(new Plane(-1,0,0,Vector(3,0,0),Material(Vector(0,0,0),Vector(0.75,0.75,0.75))));
    theMeshes.at(21)->material.reflective=true;*/
/// End korridor

/// Icosahedrons
//
//    theMeshes.push_back(new Sphere(0,-4,3,1,Material(Vector(7,7,7),Vector(1,1,1))));
//Matrix4 m;
//
//    theMeshes.push_back(new OBJ("icosahedron.obj"));
//
//    m=m.translate(-1.2,0.2,2);
//    theMeshes.at(1)->transformVertices(m);
//
//    theMeshes.at(1)->calculateNormals();
//    theMeshes.at(1)->calculateBounds();
//theMeshes.at(1)->material.reflectance=Vector(0.8,0.8,0.1);
//    theMeshes.push_back(new OBJ("icosahedron.obj"));
//
//    m=m.translate(1.2,0.2,2);
//    theMeshes.at(2)->transformVertices(m);
//
//    theMeshes.at(2)->calculateNormals();
//    theMeshes.at(2)->calculateBounds();
//
//theMeshes.at(2)->material.reflectance=Vector(0.8,0.8,0.1);
//    theMeshes.push_back(new OBJ("icosahedron.obj"));
//
//    m=m.translate(1.5,0.2,4.7);
//    theMeshes.at(3)->transformVertices(m);
//
//    theMeshes.at(3)->calculateNormals();
//    theMeshes.at(3)->calculateBounds();
//
//theMeshes.at(3)->material.reflectance=Vector(0.8,0.8,0.1);
//    theMeshes.push_back(new OBJ("icosahedron.obj"));
//
//    m=m.translate(-1.5,0.2,4.7);
//    theMeshes.at(4)->transformVertices(m);
//
//    theMeshes.at(4)->calculateNormals();
//    theMeshes.at(4)->calculateBounds();
//
//theMeshes.at(4)->material.reflectance=Vector(0.8,0.8,0.1);
//
//
//    theMeshes.push_back(new Plane(0,-1,0,Vector(0,1,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
//    theMeshes.push_back(new Plane(0,1,0,Vector(0,-4,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
//
//    theMeshes.push_back(new Plane(0,0,1,Vector(0,0,-5),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
//
//    theMeshes.push_back(new Plane(0,0,-1,Vector(0,0,16),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
//
//    theMeshes.push_back(new Plane(1,0,0,Vector(-3,0,0),Material(Vector(0,0,0),Vector(0.75,0.25,0.25))));
//
//    theMeshes.push_back(new Plane(-1,0,0,Vector(3,0,0),Material(Vector(0,0,0),Vector(0.25,0.75,0.25))));


  //  for(int i=0;i<40;i++)
  //      theMeshes.push_back(new Sphere(20*random()-10,20*random()-10,20*random()-10,1,Material(Vector(7,7,7),Vector(1,1,1))));
  //   theMeshes.push_back(new Sphere(0,0,0,1,Material(Vector(0,0,0),Vector(0.75,0.75,0.75))));
  //   theMeshes.at(0)->material.transparent=true;
//
 //    theMeshes.push_back(new Sphere(2,-3,2,1,Material(Vector(5,0.5,0.5),Vector(1,1,1))));
 //    theMeshes.push_back(new Sphere(-2,-3,2,1,Material(Vector(0.5,5,0.5),Vector(1,1,1))));
//     theMeshes.push_back(new Sphere(0,-3,-2,1,Material(Vector(0.5,0.5,5),Vector(1,1,1))));
//
//

/// icosahedron yo
//    theMeshes.push_back(new OBJ("icosahedron.obj"));
//    theMeshes.at(0)->material.reflectance=Vector(0.5,0.5,1);
//    Matrix4 m;
//
//    m=m.translate(0,0.2,0);
//      theMeshes.at(0)->transformVertices(m);
//
//    m=m.rotateY(0.4);
// theMeshes.at(0)->transformVertices(m);
//
//    theMeshes.at(0)->calculateNormals();
//    theMeshes.at(0)->calculateBounds();
//    theMeshes.push_back(new Sphere(0,-4,0,1,Material(Vector(9,9,9),Vector(1,1,1))));
  /*      theMeshes.push_back(new Plane(0,-1,0,Vector(0,1,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
    theMeshes.push_back(new Plane(0,1,0,Vector(0,-4,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(0,0,1,Vector(0,0,-10),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(0,0,-1,Vector(0,0,16),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

    theMeshes.push_back(new Plane(1,0,0,Vector(-3,0,0),Material(Vector(0,0,0),Vector(0.75,0.25,0.25))));

    theMeshes.push_back(new Plane(-1,0,0,Vector(3,0,0),Material(Vector(0,0,0),Vector(0.25,0.75,0.25))));*/
/// End icosahedron

//

/// End Random

//////

//    Matrix4 m,m2;
//   // m=m.scale(0.4);
//  //  theMeshes.at(2)->transformVertices(m);
//    m=m.scale(2,0.2,4);
//    theMeshes.at(0)->transformVertices(m);
//    m=m.translate(-1,-4,0);
//    theMeshes.at(0)->transformVertices(m);
////
//    m=m.translate(0,-0.5,4);
//    theMeshes.at(1)->transformVertices(m);

    //for(int i=0;i<80;i++)
     //   theMeshes.push_back(new Sphere(20*random()-10,20*random()-10,50*random()-50,1,Material(Vector(0,0,0),Vector(random(),random(),random()))));
//    theMeshes.push_back(new Sphere(0,0,0,1,Material(Vector(14,14,14),Vector(1,1,1))));
//    for(float i=-0.4;i<5.6;i+=0.8)
//    {
//        theMeshes.push_back(new Sphere(2.5*cos(i),0.5,2.5*sin(i),0.5,Material(Vector(0,0,0),Vector(random(),random(),random()))));
//    }
//    theMeshes.push_back(new Plane(0,0,1,Vector(0,0,-5),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
//    theMeshes.push_back(new Plane(0,-1,0,Vector(0,1,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
//

//    m=m.rotateY(0.4);
//
//    theMeshes.at(2)->transformVertices(m);
//    m=m.translate(-1.4,0.2,4);
//    theMeshes.at(2)->transformVertices(m);
//    theMeshes.at(2)->calculateNormals();
//    theMeshes.at(2)->calculateBounds();
//
//

//Balls on the wall
/*
theMeshes.push_back(new Plane(0,-1,0,Vector(0,1,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

for(int i = 0; i < 100; i++) {
		Vector color = randomMixed(Vector(1,1,1));
		double r = random()*5+1;
		double x = random()*100-50;
		double y = -r+1;
		double z = random()*200-220;
		theMeshes.push_back(new Sphere(x,y,z,r,Material(color)));
	
}*/

/*
theMeshes.push_back(new Plane(0,-1,0,Vector(0,1,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));

theMeshes.push_back(new Plane(0,0,-1,Vector(0,0,32),Material(Vector(0,0,0),Vector(1,1,1))));

theMeshes.push_back(new Sphere(0,-80.2,-50,20,Material(Vector(10,10,10), Vector(1,1,1))));

Vector green = Vector(0.2,0.7,0.2);
int index = 4;
for(int i = 0; i < 1; i++) {
	double x = 0;//random()*50-25;
	double z = 0;//random()*50-25;
	
theMeshes.push_back(new Sphere(x,0,z,1.2,Material(Vector(0.35,0.26,0.11))));
	theMeshes.push_back(new Sphere(x,-1.2,z,1.8,Material(green)));
	theMeshes.push_back(new Sphere(x,-2.4,z,1.5,Material(green)));
	theMeshes.push_back(new Sphere(x,-3.6,z,1,Material(green)));
	theMeshes.push_back(new Sphere(x,-4.6,z,0.5,Material(green)));


		theMeshes.push_back(new Sphere(-0.1,-0.2,1.8,0.2,Material(Vector(1,0,0)), true));
		theMeshes.push_back(new Sphere(1.2,-0.6,1.2,0.2,Material(Vector(1,0,0)), true));
		theMeshes.push_back(new Sphere(-1.5,-1,1,0.2,Material(Vector(1,0,0)), true));
		theMeshes.push_back(new Sphere(0.4,-1.5,1.8,0.2,Material(Vector(1,0,0)), true));

		theMeshes.push_back(new Sphere(0.8,-2.5,1.5,0.2,Material(Vector(1,0,0)), true));
		theMeshes.push_back(new Sphere(-0.8,-2.8,1.6,0.2,Material(Vector(1,0,0)), true));
		
		theMeshes.push_back(new Sphere(0.7,-3.8,1.1,0.2,Material(Vector(1,0,0)), true));

	
	theMeshes.push_back(new Sphere(x,-5.1,z,0.3,Material(Vector(5,5,5),Vector(1,1,1))));

}*/

	
/*
for(int i = 0; i < 100; i++) {
		Vector color = randomMixed(Vector(1,1,1));
		double r = random()*5+1;
		double x = random()*100-50;
		double y = -r+1;
		double z = random()*200-220;
		theMeshes.push_back(new Sphere(x,y,z,r,Material(color)));
	
}*/


/// RANDOM BALLS
/*
double temp=random();
Vector theColor=Vector(0.8*step(0.6,temp),0.78,0.9);
   // for(float i=0;i<6;i+=2.9)
float i = 0;
        theMeshes.push_back(new Sphere(5*cos(i),4*i*sin(i),0,1,Material(Vector(0,0,0),theColor)));
    Vector pos;
    double index=0;
    Vector newDir;
	theMeshes.push_back(new Sphere(0,0,0,1,Material(Vector(0,0,0),Vector(0.8*step(0.5,random()),0.7,0.8*step(0.5,random())))));
    for(int i=0;i<1600;i++)
    {

        while(true)
        {
            double size=theMeshes.size()-1;

             index=(int)(random()*size);

            if(index<=size)
                if(theMeshes.at(index)->numConnections<1)
                    break;
        }

        theMeshes.at(index)->numConnections++;
        pos=theMeshes.at(index)->position;
        newDir=Vector((random()-random()),random()-random(),(random()-random()));
        newDir.normalize();
        pos=pos+newDir;
        temp=random();
        theColor=Vector(0.7*step(0.7,temp),0.78,0.57);
        theMeshes.push_back(new Sphere(pos.x,pos.y,pos.z,1,Material(Vector(0,0,0),theColor)));

    }
	theMeshes.push_back(new Plane(0,-1,0,Vector(0,3,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
		*/

	
 // Lonely balls
/*
Material mirror = Material(Vector(0.91,0.91,0.91));
mirror.reflective = true;
theMeshes.push_back(new Sphere(-8,-5,0,8,mirror));
theMeshes.push_back(new Sphere(4,0,3,3,Material(Vector(0.12,0.24,0.86))));

theMeshes.push_back(new Sphere(6.2,1,10,2,Material(Vector(0.2,0.64,0.16))));

	    theMeshes.push_back(new Plane(0,1,0,Vector(0,-40,0),Material(Vector(1,1,1),Vector(0.8,0.8,0.8))));
		theMeshes.push_back(new Plane(0,-1,0,Vector(0,3,0),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
		
		Vector s(0,-2.0,35); // Starting point
		*/
Material mirror = Material(mixWhite(Vector(1, 0.47, 0.47)));
mirror.reflective = true;

Material greenMirror = Material(mixWhite(Vector(0.6, 1.47, 0.47)));
mirror.reflective = true;

theMeshes.push_back(new Sphere(30, -22, -40, 25, Material(mixWhite(Vector(0.64, 0.94, 0.94)))));
theMeshes.push_back(new Sphere(0, -7, 0, 10, mirror));

theMeshes.push_back(new Sphere(-5, 1, 20, 2, Material(mixWhite(Vector(1.5,1.5,0.2)))));
theMeshes.push_back(new Sphere(3, 2.5, 20, 1, Material(mixWhite(Vector(0.23, 0.6, 2)))));

theMeshes.push_back(new Sphere(-18, -4, 2, 7, greenMirror));




theMeshes.push_back(new Sphere(-20, -30, 6, 4, Material(Vector(30,30,30), Vector(1,1,1))));

theMeshes.push_back(new Plane(0, 1, 0, Vector(0, -40, 0), Material(Vector(0.9,0.9,1), Vector(1,1, 1))));
theMeshes.push_back(new Plane(0, -1, 0, Vector(0, 3, 0), Material(Vector(0, 0, 0), Vector(0.8, 0.8, 0.8))));

Vector s(0, -3, 40); // Starting point


		
//
   // theMeshes.push_back(new Plane(0,0,1,Vector(0,0,-10),Material(Vector(0,0,0),Vector(0.8,0.8,0.5))));
//
    //theMeshes.push_back(new Plane(0,0,-1,Vector(0,0,16),Material(Vector(0,0,0),Vector(0.8,0.8,0.8))));
//
   // theMeshes.push_back(new Plane(1,0,0,Vector(-60,0,0),Material(Vector(0,0,0),Vector(0.75,0.25,0.25))));
//
  // theMeshes.push_back(new Plane(-1,0,0,Vector(60,0,0),Material(Vector(0,0,0),Vector(0.25,0.75,0.25))));
	cout << "Number of meshes: "<<theMeshes.size()<<endl<<endl;
    cout << "Number of rays/pixel: "<<numRays<<endl;
	//Vector s(0,0,15); // Starting point
	Vector s2;
	float focusLength=35;
int totalRays=0;

int numberPasses=0;
for(;;)
{
	//cout << "Samples: " << totalRays << endl;
   // s=s+Vector(0,0,-0.1);

	float xmax=5,ymax=5;
	#pragma omp parallel for
    for(int screenY=0;screenY<HEIGHT;screenY++)
    {
        for(int screenX=0;screenX<WIDTH;screenX++)
        {
            Vector endColor,dir;
            float x,y;
			x=(float)(screenX*6)/(float)WIDTH-3.0;
			y=(float)(screenY*6)*(float)HEIGHT/(float)WIDTH/(float)HEIGHT-3.0*(float)HEIGHT/(float)WIDTH;

			dir=Vector(x/xmax,y/ymax,-1); // Direction
            dir.normalize();
            for(int i=0;i<numRays;i++)
            {
                bounces=0;
                s2=s+Vector(2.0*random()-1,2.0*random()-1,2.0*random()-1)*DoF;
                Vector dir2;
                Vector position2 = s + dir*focusLength;

                dir2 = position2-s2;
                dir2.normalize();

                endColor=endColor+shootRay(s2,dir2,-1); // Fire it up
            }
            endColor=endColor/numRays;
            /*if(endColor.x>1)
            	endColor.x=1;
            if(endColor.y>1)
            	endColor.y=1;
            if(endColor.z>1)
            	endColor.z=1;*/
			

            pictureAA[screenX][screenY][0]+=endColor.x;
            pictureAA[screenX][screenY][1]+=endColor.y;
            pictureAA[screenX][screenY][2]+=endColor.z;


        }
        if((int)screenY%20==0)
            cout << screenY << " / " << HEIGHT<<endl;
    }
    totalRays+=numRays;
    cout << "Total rays: " << totalRays << endl;
    numberPasses++;
    superSample(numberPasses);
	//stringstream ss;
   // ss<<totalRays;
	//string out="bilder/aa"+ss.str()+".bmp";
	string out="C:\\a\\steen2.bmp";
	char outchar[20];
	for(int a=0;a<out.length();a++)
	{outchar[a]=out[a];}
	outchar[out.length()]='\0';

    savebmp(outchar,SCREENWIDTH,SCREENHEIGHT);




}
}
/*
shootRay()
s = starting point
d = direction
index = the index of the sphere the ray comes from


*/



Vector shootRay(Vector s,Vector d,int index)
{
    if(bounces>maxBounces){
        return Vector(0,0,0);}
    bounces++;
	d.normalize();
	Vector n,c,v,endMovement,pos;
	float wee,spec,dot=0;
	float t[2];
	Vector endColor;
    Vector refractionDir=d;
    Vector refractionPos;
	float distances=9999999;
    bool hit=false;

    for(int i=0;i<theMeshes.size();i++) // Find closest intersection
    {
        if(i!=index)
        if(theMeshes.at(i)->checkIntersection(s,d,distances,pos)==true)
        {
            index=i;
            hit=true;

        }
    }

    if(hit)
    {
        Material m=theMeshes.at(index)->material;
        if(m.emittance.x>0 || m.emittance.y>0 || m.emittance.z>0)
            return m.emittance;


        n=theMeshes.at(index)->getNormal(pos);
       /* if(index==2)
        {
                m.reflectance=Vector(1,1,1)*step(0,sin(4*pos.x)*cos(4*pos.z))+Vector(0.05,0.05,0.05); // Checkers floor
        }*/
        if(m.reflectance.sum()>((float)rand()/(float)RAND_MAX)*3.0 || m.transparent==true || m.reflective==true)
        {
            if(!m.reflective && !m.transparent)
                m.reflectance = m.reflectance*(3.0/m.reflectance.sum());

            Vector emittance=m.emittance;
            Vector reflected,refracted,explicitLight,empty;

    //   // pick a random direction from here and keep going
            Vector newDir=Vector(2*random()-1,2*random()-1,2*random()-1);

            newDir=newDir.cross(theMeshes.at(index)->getNormal(pos));
            newDir.normalize();

            float eps1 = random()*3.14159*2.0f;
            float eps2 = sqrtf(random());

            float x = cosf(eps1)*eps2;
            float y = sinf(eps1)*eps2;
            float z = sqrtf(1.0f - eps2*eps2);
            Vector tempnormal=theMeshes.at(index)->getNormal(pos);
            Vector ssx= newDir * x + tempnormal.cross(newDir) * y + tempnormal * z;
            ssx.normalize();

            if(m.transparent==true)
            {

                refractionPos=theMeshes.at(index)->getRefractionPoint(pos,refractionDir);
                if((refractionPos==empty)==false)
                    refracted=refracted+shootRay(refractionPos,refractionDir,-1);

             //   refracted=refracted+shootRefractedRay(pos,d,index,1.0); // Fixa denna sen
            }
            else if(m.reflective==true && random()>0.6)
            {
                reflected=reflected+shootRay(pos,d-n*2*(n.dot(d)),-1)*m.spec;
            }
			else
			{
				reflected = reflected + shootRay(pos, ssx, -1);
			}
           return ( reflected +refracted)* m.reflectance;
        }
        else
        {
            return Vector(0,0,0);
        }
    }
 //   float nohitdot=d.dot(Vector(0,-1,0));
   //     nohitdot=sqrt(clamp(nohitdot,0,1));
	return Vector(0,0,0);
}

Vector shootRefractedRay(Vector s,Vector d,int index,float n1)
{
    d.normalize();

    float n2 = theMeshes.at(index)->material.refractionIndex;
    float distances=9999999;
    bool hit=false;
    int index2=0;
    Vector theNormal=theMeshes.at(index)->getNormal(s);
    Vector pos=s;

    float n = n1/n2;
    float cosI = theNormal.dot(d);
    float sinT2 = n * n * (1.0 - cosI*cosI);

    Vector refracted = d*n + theNormal*(n*cosI-sqrt(1-sinT2));
    refracted.normalize();

    d=refracted;
    d.normalize();
    s=s+d*0.0001;

    for(int i=0;i<theMeshes.size();i++) // Find closest intersection
    {
        if(theMeshes.at(i)->checkIntersection(s,d,distances,pos)==true)
        {
            index2=i;
            hit=true;
        }
    }



    theNormal=theMeshes.at(index2)->getNormal(pos);
    n1 = theMeshes.at(index)->material.refractionIndex;
    if(index==index2) // No other object was inside
    {
        n2 = 1.0; // No object inside, exit to air
        theNormal = theNormal*-1;

    }

    n = n1/n2;

    cosI = theNormal.dot(d);
    sinT2 = n * n * (1.0 - cosI*cosI);

    refracted = d*n + theNormal*(n*cosI-sqrt(1.0-sinT2));
    refracted.normalize();
    d=refracted;
   // pos=pos+d*0.0001;
    if(sinT2 > 1)
        return Vector(0,0,0);
    else
        return shootRay(pos,d,-1);

}



void superSample(float numPasses)
{
    double endR,endG,endB;
    int endX=0,endY=0;
    for(int y=0;y<HEIGHT;y+=aaFactor)
    {

        for(int x=0;x<WIDTH;x+=aaFactor)
        {
            endR=endG=endB=0;

            for(int a=0;a<aaFactor;a++) //x-movement
            {
                for(int b=0;b<aaFactor;b++) // y-movement
                {
					double red = pictureAA[x + a][y + b][0];
					double green = pictureAA[x + a][y + b][1];
					double blue = pictureAA[x + a][y + b][2];

                    endR+=red*255/numPasses;
                    endG+=green*255/numPasses;
                    endB+=blue*255/numPasses;

                }
            }

        endR/=(aaFactor*aaFactor);
        endG/=(aaFactor*aaFactor);
        endB/=(aaFactor*aaFactor);
		endR = clamp(endR, 0, 255);
		endG = clamp(endG, 0, 255);
		endB = clamp(endB, 0, 255);
        picture[endX][endY][0]=endR;
        picture[endX][endY][1]=endG;
        picture[endX][endY][2]=endB;

        endX++;
        }
        endX=0;
        endY++;
    }

}


void savebmp( const char *filename, int w, int h )
{
	/*std::ofstream theFile;
	theFile.open(filename);
	// Write header
	theFile << "P6 " << w << " " << h << " " << 255 << "\n";

	for(int i=0;i<SCREENHEIGHT;i++)
    {
        for(int j=0;j<SCREENWIDTH;j++)
        {
            theFile << (unsigned char)picture[j][i][2] << " ";
            theFile << (unsigned char)picture[j][i][1] << " ";
            theFile << (unsigned char)picture[j][i][0] << " ";
            
        }
		theFile << "\n";
    }
	theFile.close();*/
	
    int i;
    FILE *f;
    int filesize = 54 + 3*w*h;

    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize    );
    bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize>>24);

    bmpinfoheader[ 4] = (unsigned char)(       w    );
    bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       w>>16);
    bmpinfoheader[ 7] = (unsigned char)(       w>>24);
    bmpinfoheader[ 8] = (unsigned char)(       h    );
    bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
    bmpinfoheader[10] = (unsigned char)(       h>>16);
    bmpinfoheader[11] = (unsigned char)(       h>>24);

    //f = fopen("img.raw","wb");
    //fwrite(img,3,w*h,f);
    //fclose(f);

    // Make the int-array into a char-array for the bmp
    int counter=0;
    for(int i=0;i<SCREENHEIGHT;i++)
    {
        for(int j=0;j<SCREENWIDTH;j++)
        {
            img[counter]=(unsigned char)picture[j][i][2];
            img[counter+1]=(unsigned char)picture[j][i][1];
            img[counter+2]=(unsigned char)picture[j][i][0];
            counter=counter+3;
        }
    }

    f = fopen(filename,"wb+");
    fwrite(bmpfileheader,1,14,f);
    fwrite(bmpinfoheader,1,40,f);

    for(i=0; i<h; i++)
    {
        fwrite(img+(w*(h-i-1)*3),3,w,f);
        fwrite(bmppad,1,(4-(w*3)%4)%4,f);
    }
    fclose(f);
}
