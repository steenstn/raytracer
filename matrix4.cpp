#include "matrix4.h"
#include "vertex.h"
#include <math.h>
class Matrix4;

Matrix4::Matrix4()
{
    a11=1;   a12=0;   a13=0;   a14=0;
    a21=0;   a22=1;   a23=0;   a24=0;
    a31=0;   a32=0;   a33=1;   a34=0;
    a41=0;   a42=0;   a43=0;   a44=1;
}

Matrix4::Matrix4(const Matrix4& m)
{
    a11=m.a11; a12=m.a12; a13=m.a13; a14=m.a14;
    a21=m.a21; a22=m.a22; a23=m.a23; a24=m.a24;
    a31=m.a31; a32=m.a32; a33=m.a33; a34=m.a34;
    a41=m.a41; a42=m.a42; a43=m.a43; a44=m.a44;
}

Matrix4 Matrix4::rotateX(double theta)
{
    Matrix4 m;
    m.a11=1; m.a12=0;               m.a13=0; 				  m.a14=0;
    m.a21=0; m.a22=cos(theta);      m.a23=-1*sin(theta);      m.a24=0;
    m.a31=0; m.a32=sin(theta);      m.a33=   cos(theta);      m.a34=0;
    m.a41=0; m.a42=0;               m.a43=0;                  m.a44=1;

    return m;
}

Matrix4 Matrix4::rotateY(double theta)
{
    Matrix4 m;
    // Set appropriate elements here
    m.a11=cos(theta);	  m.a12=0;		m.a13=sin(theta); 	  m.a14=0;
    m.a21=0; 		      m.a22=1; 	    m.a23=0; 			  m.a24=0;
    m.a31=-1*sin(theta);  m.a32=0; 		m.a33=cos(theta); 	  m.a34=0;
    m.a41=0; 			  m.a42=0;      m.a43=0;              m.a44=1;
    return m;
}
Matrix4 Matrix4::translate(double tx, double ty, double tz) {
    Matrix4 m;
    m.a14=tx; m.a24=ty; m.a34=tz;
    return m;
}

Matrix4 Matrix4::scale(double s)
{
    Matrix4 m;
    m.a11=s;
    m.a22=s;
    m.a33=s;
    return m;
}

Matrix4 Matrix4::scale(double sx,double sy,double sz)
{
    Matrix4 m;
    m.a11=sx;
    m.a22=sy;
    m.a33=sz;
    return m;
}


Vertex Matrix4::mult(Vertex v)
{
    Vertex vout;
    // Calculate x,y,z,w values of vout here
    vout.x=a11*v.x+a12*v.y+a13*v.z+a14;
    vout.y=a21*v.x+a22*v.y+a23*v.z+a24;
    vout.z=a31*v.x+a32*v.y+a33*v.z+a34;
    vout.w=a41*v.x+a42*v.y+a43*v.z+a44;

    return vout;
}



Matrix4& Matrix4::operator= (const Matrix4& m)
{

    a11=m.a11; a12=m.a12; a13=m.a13; a14=m.a14;
    a21=m.a21; a22=m.a22; a23=m.a23; a24=m.a24;
    a31=m.a31; a32=m.a32; a33=m.a33; a34=m.a34;
    a41=m.a41; a42=m.a42; a43=m.a43; a44=m.a44;

	return *this;
}
