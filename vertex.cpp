#include "vertex.h"
class Vertex;


    Vertex::Vertex()
    {
        w=1;
        x=y=z=0;
    }
    Vertex::Vertex(double theX,double theY,double theZ)
    {
        x=theX;
        y=theY;
        z=theZ;
        w=1;
    }
    Vertex Vertex::operator- (const Vertex& v)
    {
        Vertex temp((*this).x,(*this).y,(*this).z);
        temp.x-=v.x;
        temp.y-=v.y;
        temp.z-=v.z;
        return temp;
    }
    void Vertex::normalizeW()
    {
        x/=w;
        y/=w;
        z/=w;
        w/=w;
    }

