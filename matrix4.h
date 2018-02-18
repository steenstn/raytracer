#ifndef Matrix4_h
#define Matrix4_h
#include "vertex.h"
 class Matrix4 {

  public:
   float a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44;

    Matrix4();
    Matrix4(const Matrix4& m);
    static Matrix4 rotateX(float theta);
    static Matrix4 rotateY(float theta);
    static Matrix4 translate(float tx, float ty, float tz);
    static Matrix4 scale(float s);
    static Matrix4 scale(float sx,float sy,float sz);


    Vertex mult(Vertex v);

    Matrix4& operator= (const Matrix4& m);
};

#endif
