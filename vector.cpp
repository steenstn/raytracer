#include "vector.h"
class Vector;

    Vector::Vector()
    {
        x=0.0;
        y=0.0;
        z=0.0;
    }
    Vector::Vector(const Vector& v)
    {
        x=v.x;
        y=v.y;
        z=v.z;
    }
    Vector::Vector(const Vertex& v)
    {
        x=v.x;
        y=v.y;
        z=v.z;
    }
    Vector::Vector(double theX,double theY,double theZ)
    {
        x=theX;
        y=theY;
        z=theZ;
    }
    double Vector::dot(Vector b)
    {
        return x*b.x+y*b.y+z*b.z;
    }
    Vector Vector::cross(Vector b)
    {
        return Vector(y*b.z-z*b.y,-1*(x*b.z-z*b.x),x*b.y-y*b.x);
    }
    void Vector::normalize(void)
    {
        double abs=sqrt(x*x+y*y+z*z);
        x/=abs;
        y/=abs;
        z/=abs;
    }

    double Vector::length(void)
    {
        return sqrt(x*x+y*y+z*z);

    }

    double Vector::sum(void)
    {
        return x+y+z;
    }

    Vector Vector::operator+ (const Vector& v)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x+=v.x;
        temp.y+=v.y;
        temp.z+=v.z;
        return temp;
    }

    Vector Vector::operator* (const Vector& v)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x*=v.x;
        temp.y*=v.y;
        temp.z*=v.z;
        return temp;
    }

    Vector Vector::operator- (const Vector& v)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x-=v.x;
        temp.y-=v.y;
        temp.z-=v.z;
        return temp;
    }



    Vector Vector::operator* (const double &d)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x*=d;
        temp.y*=d;
        temp.z*=d;
        return temp;
    }
    Vector Vector::operator/ (const double &d)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x/=d;
        temp.y/=d;
        temp.z/=d;
        return temp;
    }
    Vector Vector::operator- (const double &d)
    {
        Vector temp((*this).x,(*this).y,(*this).z);
        temp.x-=d;
        temp.y-=d;
        temp.z-=d;
        return temp;
    }

    bool Vector::operator== (Vector& d)
    {
        if(x==d.x && y==d.y && z==d.z)
			return true;
		return false;
    }


    const Vector& Vector::operator= (const Vector& v)
    {

            (*this).x=v.x;
            (*this).y=v.y;
            (*this).z=v.z;

			return *this;
    }

    std::ostream& operator<<(std::ostream& o,Vector& v)
    {
        o<<"["<< v.x << " " << v.y << " " << v.z << "]";
        return o;
    }


