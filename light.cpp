#include "light.h"
class Light;

    Light::Light()
    {
    }

    Light::Light(float theX,float theY,float theZ,Vector theColor)
    {
        position.x=theX;
        position.y=theY;
        position.z=theZ;
        direction=Vector();
        color=theColor;
        intensity=1;
    }

    Light::Light(float theX,float theY,float theZ,float theIntensity)
    {
        position.x=theX;
        position.y=theY;
        position.z=theZ;
        direction=Vector();
        color=Vector(1,1,1);
        intensity=theIntensity;
    }

    Light::Light(float theX,float theY,float theZ,float theIntensity,Vector theColor)
    {
        position.x=theX;
        position.y=theY;
        position.z=theZ;
        direction=Vector();
        color=theColor;
        intensity=theIntensity;
    }

    void Light::calculateDirection(const Vector &pos)
    {
        direction.x=pos.x-position.x;
        direction.y=pos.y-position.y;
        direction.z=pos.z-position.z;
		direction.normalize();
    }

