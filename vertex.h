#ifndef vertex_h
#define vertex_h
class Vertex
{
    public:
    float x,y,z,w;
    Vertex();
    Vertex(float theX,float theY,float theZ);
    Vertex operator- (const Vertex& v);
    void normalizeW();

};
#endif
