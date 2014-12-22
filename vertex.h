#ifndef vertex_h
#define vertex_h
class Vertex
{
    public:
    double x,y,z,w;
    Vertex();
    Vertex(double theX,double theY,double theZ);
    Vertex operator- (const Vertex& v);
    void normalizeW();

};
#endif
