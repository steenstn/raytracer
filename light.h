#ifndef Light_H
#define Light_H
#include "vector.h"
class Light
{
    public:
    float x,y,z,intensity;
    struct Vector position;
    struct Vector direction;
    struct Vector color;

    Light();

    Light(float theX,float theY,float theZ,Vector theColor);
    Light(float theX,float theY,float theZ,float theIntensity);
    Light(float theX,float theY,float theZ,float theIntensity,Vector theColor);
    void calculateDirection(const Vector &pos);
};
#endif
