/**
  * Collection of routines and datatypes used for the GPU-Calculations.
  * Normally without any Qt-Classes.
  */

#ifndef SC_MATH_H
#define SC_MATH_H

#include <math.h>
#include <stdio.h>
//#ifndef __CUDACC__
//#include <QString>      // used for debugging
//#endif

#ifdef Q_OS_LINUX
#define EOL "\n"
#else
#define EOL "\r\n"
#endif


typedef bool (*progressAndAbort)(int);  //!< Callback function type to update the progress bar and get the abort information
                                        //!< ... only used without GPU


typedef bool (*progressLogging)(char*);  //!< Callback function type to put informations in the logging list
                                         //!< only used during the Fit2D function


// Welcher Quadrant wird berechnet? Wird nur bei der Bestimmung von zzmin und iimin verwendet.
typedef enum { radQ1,       //!< only the first quadrant is used (x,y >= 0)
               radQ2,       //!< two quadrants are used (y >= 0)
               radQ4        //!< all four quadrants are used
             } _radQX;      //!< Type definition of the used quadrants


#define sqr(x) ((x)*(x))            //!< helper definition


const double eps3 = 0.001;
const double eps4 = 0.0001;
const double eps5 = 0.00001;
const double eps6 = 0.000001;
const double eps7 = 0.0000001;
const double eps8 = 0.00000001;
const double eps9 = 0.000000001;

const double eps1 = 0.025;          //!< used in trapezchi() and integralchi() (FCC)


#define odd(n) (((n)&1) == 1)
#define even(n) (((n)&1) == 0)


// Optimierung von cos(M_PI*(h+k+l)):
//#define cospi(v) cos(M_PI*(v))
#define cospi(v) (1-2*((v)&1))



/**
 * @brief The Double3 class
 * Holds a vector of three double values and perform the calculations.
 */
class Double3
{
public:
#ifdef __CUDACC__
    __host__ __device__
#endif
    Double3( double x=0.0, double y=0.0, double z=0.0 ) { _x=x; _y=y; _z=z; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline double x() const { return _x; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline double y() const { return _y; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline double z() const { return _z; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline void setX( double x ) {	 _x=x; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline void setY( double y ) {	 _y=y; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline void setZ( double z ) {	 _z=z; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline double length() const { return sqrt(sqr(_x)+sqr(_y)+sqr(_z)); }

//#ifndef __CUDACC__
//    QString toString() const { return QString("Double3(%1,%2,%3)").arg(_x).arg(_y).arg(_z); }
//#endif
    const char *toString()
    {
        snprintf ( buf, sizeof(buf), "(%g,%g,%g)", _x, _y, _z );
        // Der Compiler auf dem Bürorechner hat ermittelt, dass zwischen 8 und 44 Zeichen
        // ausgegeben werden können. Daher war die ursprüngliche Variable mit 40 Zeichen zu klein...
        return buf;
    }

#ifdef __CUDACC__
    __host__ __device__
#endif
    Double3& operator=(const Double3 &v)
    {
        _x = v.x();
        _y = v.y();
        _z = v.z();
        return *this;
    }

#ifdef __CUDACC__
    __host__ __device__
#endif
    Double3 operator+(const Double3 &other) const
    {
        return Double3( _x+other.x(), _y+other.y(), _z+other.z() );
    }

#ifdef __CUDACC__
    __host__ __device__
#endif
    Double3 operator-(const Double3 &other) const
    {
        return Double3( _x-other.x(), _y-other.y(), _z-other.z() );
    }

#ifdef __CUDACC__
    __host__ __device__
#endif
    Double3 operator*(const double v) const
    {
        return Double3( _x*v, _y*v, _z*v );
    }

#ifdef __CUDACC__
    __host__ __device__
#endif
    Double3 operator/(const double v) const
    {
        if ( v != 0 )
            return Double3( _x/v, _y/v, _z/v );
        return Double3(); // Durch 0 ergibt erstmal nur 0 ...
    }

#ifdef __CUDACC__
    __host__ __device__
#endif
    friend Double3 operator/(double v, Double3 &d)
    {
        Double3 erg;
        if ( d.x() != 0 ) erg.setX( v / d.x() );
        if ( d.y() != 0 ) erg.setY( v / d.y() );
        if ( d.z() != 0 ) erg.setZ( v / d.z() );
        return erg;
    }

private:
    double _x, _y, _z;
    char buf[50]; // for toString()
};

#endif // SC_MATH_H
