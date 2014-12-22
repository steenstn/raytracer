#include "material.h"
class Material;


Material::Material()
{
    emittance=Vector();
    reflectance=Vector(0,0.4,0.6);
    transparent=false;
    reflective=false;
    refractionIndex=0;
}
Material::Material(Vector theReflectance)
{
    transparent=false;
    emittance=Vector();
    reflective=false;
    reflectance=theReflectance;
    refractionIndex=0;
    spec=Vector(1,1,1);
}
Material::Material(Vector theEmittance,Vector theReflectance)
{
    transparent=false;
    reflectance=theReflectance;
    emittance=theEmittance;
    refractionIndex=0;
    reflective=false;
    spec=Vector(1,1,1);
}


