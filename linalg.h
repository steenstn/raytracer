#ifndef linalg_h
#define linalg_h
#include <math.h>
#include <iostream>
struct Vector
{
    float x,y,z;

    Vector()
    {
        x=0.0;
        y=0.0;
        z=0.0;
    }
    Vector(const Vector& v)
    {
        x=v.x;
        y=v.y;
        z=v.z;
    }
    Vector(float theX,float theY,float theZ)
    {
        x=theX;
        y=theY;
        z=theZ;
    }
    float dot(Vector b)
    {
        return x*b.x+y*b.y+z*b.z;
    }
    Vector cross(Vector b)
    {
        Vector3 cv ;
        cv.x=v1.y*v2.z-v1.z*v2.y;
        cv.y=v1.z*v2.x-v1.x*v2.z;
        cv.z=v1.x*v2.y-v1.y*v2.x;
        return cv;
    }

    void normalize(void)
    {
        float abs=sqrt(x*x+y*y+z*z);
        x/=abs;
        y/=abs;
        z/=abs;
    }

    float length(void)
    {
        return sqrt(x*x+y*y+z*z);

    }

    Vector operator+ (const Vector& v)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x+=v.x;
        temp.y+=v.y;
        temp.z+=v.z;
        return temp;
    }

    Vector operator* (const Vector& v)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x*=v.x;
        temp.y*=v.y;
        temp.z*=v.z;
        return temp;
    }

    Vector operator- (const Vector& v)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x-=v.x;
        temp.y-=v.y;
        temp.z-=v.z;
        return temp;
    }
    Vector operator* (const float &d)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x*=d;
        temp.y*=d;
        temp.z*=d;
        return temp;
    }
    Vector operator/ (const float &d)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x/=d;
        temp.y/=d;
        temp.z/=d;
        return temp;
    }
    Vector operator- (const float &d)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x-=d;
        temp.y-=d;
        temp.z-=d;
        return temp;
    }

    bool operator== (Vector& d)
    {
        if(x==d.x && y==d.y && z==d.z)
			return true;
		return false;
    }


    const Vector& operator= (const Vector& v)
    {

            (*this).x=v.x;
            (*this).y=v.y;
            (*this).z=v.z;


    }

    friend std::ostream& operator<<(std::ostream& o,Vector& v)
    {
        o<<"["<< v.x << " " << v.y << " " << v.z << "]";
        return o;
    }

};

struct Light
{
    float x,y,z,intensity;
    struct Vector position;
    struct Vector direction;
    struct Vector color;

    Light()
    {
    }

    Light(float theX,float theY,float theZ,Vector theColor)
    {
        position.x=theX;
        position.y=theY;
        position.z=theZ;
        direction=Vector();
        color=theColor;
        intensity=1;
    }

    Light(float theX,float theY,float theZ,float theIntensity)
    {
        position.x=theX;
        position.y=theY;
        position.z=theZ;
        direction=Vector();
        color=Vector(1,1,1);
        intensity=theIntensity;
    }

    Light(float theX,float theY,float theZ,float theIntensity,Vector theColor)
    {
        position.x=theX;
        position.y=theY;
        position.z=theZ;
        direction=Vector();
        color=theColor;
        intensity=theIntensity;
    }

    void calculateDirection(const Vector &pos)
    {
        direction.x=pos.x-position.x;
        direction.y=pos.y-position.y;
        direction.z=pos.z-position.z;
		direction.normalize();
    }
};

struct Material
{
    Vector diffuseColor;
    Vector specularColor;
    bool reflective;
    int specularPower;
    Material()
    {
    	diffuseColor=Vector(1.0,0.5,0.9);
    	specularColor=Vector(1,1,1);
    	reflective=false;
    	specularPower=400;
    }
    Material(Vector theDiff)
    {
    	diffuseColor=theDiff;
    	specularColor=Vector(1,1,1);
    	reflective=false;
    	specularPower=400;
    }
    Material(Vector theDiff,int theSpec)
    {
    	diffuseColor=theDiff;
    	specularColor=Vector(1,1,1);
    	reflective=false;
    	specularPower=theSpec;
    }
    Material(Vector theDiff,int theSpec,bool theRefl)
    {
    	diffuseColor=theDiff;
    	specularColor=Vector(1,1,1);
    	reflective=theRefl;
    	specularPower=theSpec;
    }
};

struct Sphere
{
    float radius;
    Vector position;
	Material material;
	// 0 = Diffuse
    Sphere()
    {

    }

    Sphere(float theX,float theY,float theZ, float theRadius)
    {
        position=Vector(theX,theY,theZ);
        radius=theRadius;

    }
    Sphere(float theX,float theY,float theZ, float theRadius,Material theMaterial)
    {
        position=Vector(theX,theY,theZ);
        radius=theRadius;
        material=theMaterial;
    }




};

struct Vertex4
{
    float x,y,z,w;
    Vertex4()
    {
        x=y=z=0;
        w=1;
    }
    Vertex4(float a,float b,float c)
    {
        x=a;
        y=b;
        z=c;
    }
    Vertex4(float a,float b,float c,float d)
    {
        x=a;
        y=b;
        z=c;
        w=d;
    }
    friend std::ostream& operator<<(std::ostream& o,Vertex4& v)
    {
        o<<"["<< v.x << " " << v.y << " " << v.z << "]";
        return o;
    }

};

struct Triangle
{
    int p1,p2,p3;

    Triangle()
    {
        p1=p2=p3=-1;
    }

    Triangle(int a,int b,int c)
    {
        p1=a;
        p2=b;
        p3=c;
    }
};
#endif