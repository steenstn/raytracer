#ifndef Matrix4_h
#define Matrix4_h
#include "vertex.h"
 class Matrix4 {

  public:
   double a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44;

    Matrix4();
    Matrix4(const Matrix4& m);
    static Matrix4 rotateX(double theta);
    static Matrix4 rotateY(double theta);
    static Matrix4 translate(double tx, double ty, double tz);
    static Matrix4 scale(double s);
    static Matrix4 scale(double sx,double sy,double sz);


    Vertex mult(Vertex v);

    Matrix4& operator= (const Matrix4& m);
};

#endif
