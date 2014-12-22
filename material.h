#ifndef Material_H
#define Material_H
#include "vector.h"
class Material
{
    public:
    Vector emittance;
    Vector spec;
    Vector reflectance;
    bool transparent;
    double refractionIndex;
    bool reflective;
    Material();
    Material(Vector theReflectance);
    Material(Vector theEmittance,Vector theReflectance);

};
#endif
