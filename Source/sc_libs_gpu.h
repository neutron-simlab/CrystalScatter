
// THIS IS ONLY AN INCLUDE FILE FOR ANOTHER CLASS

#ifndef __CUDACC__
#include <QDebug>
//#include <QElapsedTimer>
#define DD(x)  //x   // Allgemeiner Output
#define D8(x)  //x   // D8 im Main
#define D8L(x) //x   // D8 in der Library
#define DCL(x) //x   // Library-Counter
#define DTL(x) //x   // Timing
#define DSM(x) //x   // Ausgaben bei den kleinen Zahlen
#define D1P(x) //x  // carr1p
#define D4P(x) //x  // carr4p
//#define DBGLIMIT 50000
#include <time.h>
#endif

#ifdef __CUDACC__
#define DD(x)
#define D8(x)
#define D8L(x)
#define DCL(x)
#define DTL(x)
#define DSM(x)
#define D1P(x)
#define D4P(x)
#endif


#ifndef USE_CLASSLIB
#include "sc_math.h"
class CLASSLIB
{
#endif

public:
    void endThread();

#ifndef USE_CLASSLIB
private:
#endif

    typedef struct
    {
        int cols;
        int rows;
        double centerx;
        double centery;
        double wavelen;
        double distance;
        double pixx;
        double pixy;
    } _latticeForFit;

#define  coeffarray_len  (150+1)
#define  imax2d_len  (130+1)
    typedef double CoeffArrayType[coeffarray_len]; // =array[0..150] of extended;
    typedef double ArrayImax2D[imax2d_len][imax2d_len];
    typedef double ArrayImax1D[20000+1];     // TODO ist das so groß nötig?

    // 20220311: Die vielen kleinen Arrays fasse ich zu einem Record zusammen
    typedef struct
    {
        double carr1p[coeffarray_len], carr2p[coeffarray_len], carr3p[coeffarray_len], carr4p[coeffarray_len],
            carr5p[coeffarray_len], carr6p[coeffarray_len], carr7p[coeffarray_len], carr8p[coeffarray_len],
            carr9p[coeffarray_len], carr1f[coeffarray_len], carr2f[coeffarray_len], carr3f[coeffarray_len],
            carr4f[coeffarray_len], carr5f[coeffarray_len], carr6f[coeffarray_len], carr7f[coeffarray_len],
            carr8f[coeffarray_len], carr9f[coeffarray_len];
        double myarray[20]; // TODO: ist 20 ausreichend?
        // Array-Parameter für coefficients(), formpq(), formfq()
        ArrayImax2D carr11pm, carr22pm;
    } _carrXX;

    // to reduce the number of parameters in each call,
    // some global variables are stored here for all classes:
    typedef struct
    {
        double rho;     // wurde bislang noch nicht gesetzt (neu 16.11.2021)
        double p1;
        double radius;
        double radiusi;
        double sigma;
        double sigmal;
        double length; // cylinder length!
        double wavelength;
        double shellno;
        double alphash;
        double ceff;        // often equal with "TwRatio" -> TODO
        double reff;        // allways const 10.0 (no input in GUI) -> TODO?
        double uca; // latt_a;      // a
        double ucb; // latt_b;      // b
        double ucc; // latt_c;      // c
        double domainsize;
        double width_zuf; // Zufalls-Faktor für domainsize/width
        double aziwidth;
        double dbeta; //{NV} NEU !!! TODO
        double psphere_r, psphere_z;  // helper var for psphere
        // Globale Parameter für ButtonHKLClick()
        int    hklmax;
        double alpha_deg, beta_deg, gamma_deg;
        double amax, bmax, cmax;

        _carrXX  *CR;

        _latticeForFit latFit;

    } _localParams;

    _localParams params;
    size_t arrCRSize;


    // Der QtCreator zeigt hier viele unbekannte Variablen an. Beim Compilieren klappt aber alles,
    //  da dieses File in die richtige Quelle eingebunden wird.
    inline void performImageExpansion()
    {
        if ( bExpandImage && calcQuadrants != radQ4 )
        {   // Punktsymetrie ausnutzen
            //zzmin = calcQuadrants == radQ4 ? -zmax : 0;     //{NV}-10771/11952 im Batch-Loop
            //zzmax = zmax;
            //iimin = calcQuadrants == radQ1 ? 0 : -zmax;
            //iimax = zmax;

            if ( calcQuadrants == radQ1 )
            {   // Quadrant von 0,0 - zmax,zmax in alle anderen kopieren
                for ( int ihex=0; ihex<zmax; ihex++ )
                    for ( int i=0; i<zmax; i++ )
                    {
                        setXYIntensity( -ihex,  i  , xyIntensity(ihex,i) );
                        setXYIntensity( -ihex, -i-1, xyIntensity(ihex,i) );
                        setXYIntensity(  ihex, -i-1, xyIntensity(ihex,i) );
                    }
            }
            else
            {   // Bildhälfte von 0,-zmax - zmax,zmax in die andere Hälfte kopieren
                for ( int ihex=0; ihex<zmax; ihex++ )
                    for ( int i=-zmax+1; i<zmax; i++ )
                        setXYIntensity( -ihex, -i, xyIntensity(ihex,i) );
            }
        }
    }

    static const int latparlen = 5000;
    #define latparIDX(a,b,l) (a*l+b)

    static const int np=5;
    static const int jmaxp=21;
    typedef double RealArrayNP[np+2]; // =array[1..np] of extended;
    typedef double RealArrayJMAXP[jmaxp+2]; // =array[1..jmaxp] of extended;
    // Die Arrays zur Sicherheit eins länger, da in den Routinen of vor dem
    // Ende der Schleife auf das nächste Element geschrieben wird.

#ifdef __CUDACC__
__host__ __device__
#endif
    inline double frac(double x) const
    {
        double intpart;
        return modf(x, &intpart);
    }


#ifdef __CUDACC__
    __host__ __device__
#endif
    double szave( int c, double a, double x, int n, double z ) const;

#ifdef __CUDACC__
__host__
#endif
    void psphere_init();

#ifdef __CUDACC__
__host__ __device__
#endif
    double psphere( double q ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double psphered(double r, double sigma, int dim, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double pspheredf(double r, double sigma, double d, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double f2dschulz( int d, double r, double sigma, double q ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double f2dhyper(double d, double alf, double r, double sigma, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double gamma(double z) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double gammaratio(double a, double b, double x) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double cosav(double di, double dj, double ei, double ej, double y, double z) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double pqcoreshell( double alf2, double d, double q ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double pqcoreshellf( double rho1, double rho2, double p1, double p2,
                         double alf1, double alf2, double rn, double d,
                         double sigma, double q ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double f2dcoreshell( double alf, double q ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double polyvesicle(double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double f2dpolyvesicle(double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double polyliposome(int dim, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double polycube(bool pf, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double burger(double del, double v, double x2) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double angleuv(double ux, double uy, double uz, double vx,double vy, double vz) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double lorentznorm3(double a) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double gaussnorm3(double a) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double pearsonnorm3(double a, double b) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void   pearsonintegral3(double at, double bt, double a, double b,
                            double &sx, int nx, int &trapzditx) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double cosavm(double d1, double d2, double e1, double e2, double p1, double p2, double y, double z) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double pgensphere1(double a, double r, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double hypcoreshellm(double alf1, double alf2, double p1, double p2, double r, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void trapcube(double ac, double bc, double a, bool pf, double q,
                  double &sc, int nc, int &trapzditc) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double cubealf(double a, double alfa, bool pf, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void trapcubebet(double ac, double bc, double a, double alfa, bool pf, double q,
                     double &sc, int nc, int &trapzditc) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double cube(double a, double alfa, double beta, double sigma, bool pf, double q) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void qrombpq( double l, double r, double p1, double rho, double alfa, double sigmal, double q,
                  int maxit, int i0, int i1, int i3, int i4,
                  double &pq ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void qrombpi2( double epsi, double r, double p1, double rho, double alfa, double sigmal, double sigma, double q,
                   int maxit, int i0, int i1, int i3, int i4, double &pq ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void polint( double *xa/*RealArrayNP*/, double *ya/*RealArrayNP*/,
                 int n, double x, double &y, double &dy, const char *dbg ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void trapezpqr( double at, double bt, double l, double r, double p1, double rho, double alfa, double sigmal, double q, int,
                    int i0, int i3, int i4, double &pq, double &nn, int n, int &trapzdit ) const; // in qrombpq()

#ifdef __CUDACC__
__host__ __device__
#endif
    void midpntchi( double a, double b, double l, double r, double p1, double rho, double alfa, double sigmal,
                    double sigma, double cc0, double cc1, double cc2, double cc3, double cc4, double cc5, double cc6,
                    double cc7, double cc8, double delta, double phi, double theta, double qxn, double qyn, double qzn,
                    double q, int i0, int i1, int i3, int i4,  double &sp, double &sf, int n, int &midpntit ) const; // in qrombpq()

#ifdef __CUDACC__
__host__ __device__
#endif
    double bessel10( double x ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void orientcylchi( double l, double r, double p1, double rho, double alfa, double sigmal, double sigma, double cc0,
                       double cc1, double cc2, double cc3, double cc4, double cc5, double cc6, double cc7, double cc8,
                       double delta, double alpha, double phi, double theta, double qxn, double qyn, double qzn,
                       double q, int i0, int i1, int i3,  double &pa, double &fa ) const;


    void ButtonHKLClick(int ltype , int *latpar1, int *latpar2) const;

    void fhkl_c( int lat, int h, int k, int l,
                 double uca, double ucb, double ucc, double ucalpha_deg, double ucbeta_deg, double ucgamma_deg,
                 double &sphno, double &fhkl, double &qhkl, double &qhkl0 ) const;

    void extinction( int lat, int h, int k, int l, int aniso,
                     int &mhkl, double &fhkl ) const;



#ifdef __CUDACC__
__host__ __device__
#endif
//    void qrombdeltac(double r, double theta, double phi, double qx, double qy, double qz,
//                     double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
//                     double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
//                     double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
//                     int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
//                     double &pq, int threadid ) const;
void qrombdeltac( double l, double r, /*p1*/ /*sigma*/ /*dbeta*/ double theta, double phi, double qx, double qy, double qz,
                  double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
                  double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
                  double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
                  int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                  double *carr1,
                  double &pq ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void trapzddeltac( double a, double b, double l, double r, double dbeta, double theta, double phi,
                       double qx, double qy, double qz, double p11, double p12, double p13, double p21,
                       double p22, double p23, double p31, double p32, double p33,
                       double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
                       double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
                       double ax3y, double ax3z, double sigx, double sigy, double sigz,
                       int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                       double *carr1,
                       double &pq, int n, int &trapzdit ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void trapzdchid( double a, double b, double l, double r, double sigma, double dbeta, double delta, double phi, double,
                     double qx, double qy, double qz, double p11, double p12, double p13, double p21, double p22, double p23,
                     double p31, double, double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n,
                     double ax2n, double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
                     double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz, int dim, int i0, int, int i1, int, int i3, int i4, double *carr1, double &pq, int n,
                     int &trapzdit ) const;
void trapzdchid_OLDVERSION( double a, double b, double r, double sigma, double dbeta, double delta, double theta, double phi,
                 double qx, double qy, double qz, double p11, double p12, double p13, double p21, double p22, double p23,
                 double p31, double p32, double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n,
                 double ax2n, double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
                 double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
                 int ordis, int dim, int i0, int i1, int i2, int i3, int i4, double *carr1, double &pq, int n,
                 int &trapzdit ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void qrombchid( double l, double r, double sigma, double dbeta, double delta, double theta, double phi,
                    double qx, double qy, double qz, double p11, double p12, double p13, double p21,
                    double p22, double p23, double p31, double p32, double p33,
                    double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
                    double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
                    double ax3y, double ax3z, double sigx, double sigy, double sigz,
                    int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                    double *carr1,
                    double &pq ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
//    double formfq(double radius, double zz, double limq1, double qx, double qy,
//                  int part, int cs,
//                  double *carr3p /*: CoeffArrayType*/ ) const;
double formfq( double length, double radius, double sigmal, double sigmar, double p1, double rho,
               double alfa, double theta, double phi, double limql, double limq1, double limq2,
               double limq3, double limq4, double limq5, double limq6, double qx, double qy,
               double qxs, double qys, double q, double norm,
               int part, int cs, int ordis, int orcase,
               const double *myarray, // CoeffArrayType
               double *carr1p, double *carr2p, double *carr3p, double *carr4p, double *carr5p,
               double *carr6p, double *carr7p //: CoeffArrayType;   /*Z0311=17704*/,
               /*ArrayImax2D &carr11pm, ArrayImax2D &carr22pm*/ ) const; //: ArrayImax2D);   /*Z0311=17705*/


#ifdef __CUDACC__
__host__ __device__
#endif
//    double formpq(double length, double radius, double p1, double rho, double zz, double qlim, double limq1, double limq2, double limq3, double limq4, double qx, double qy, double ql, double q, double norm,
//                  int part, int cs, int ordis, int orcase,
//                  double *carr1p, double *carr2p, double *carr3p, double *carr4p, /*: CoeffArrayType;*/
//                  double *carr5p /*: ArrayImax1D*/ ) const;
    double formpq(double length, double radius, double sigmal, double sigmar, double p1,
                        double rho, double alfa, double theta, double phi, double limql, double limq1,
                        double limq2, double limq3, double limq4, double limq5, double limq6,
                        double limq7, double limq8, double limq9, double qx, double qy, double qxs,
                        double qys, double q, double norm, double por,
                        int part, int cs, int ordis, int orcase,
                        const double *myarray, // CoeffArrayType
                        double *carr1p, double *carr2p, double *carr3p, // CoeffArrayType
                        double *carr4p, double *carr5p, double *carr6p, // CoeffArrayType
                        double *carr7p, double *carr8p, double *carr9p) const;   /*Z0311=14583*/


#ifdef __CUDACC__
__host__ __device__
#endif
    double qrombpi2_funcasy1( double beta, double epsi, double r, double sigma, double d, double q ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double qrombpi2_funcasy2( double beta, double epsi, double rn, double p1, double rho2, double sigma,
                              double d, int i3, double q ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void qrombpi2_trapezpqpi2( double at, double bt, double epsi, double r, double p1, double rho, double alfa,
                               double sigmal, double sigma, double q, int i0, int i1, int i3, int i4,
                               double &pq, double &nn, int n, int &trapzdit ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    void qrombpi2_midpntchi( double a, double b, double l, double r, double p1, double rho, double alfa, double sigmal,
                             double sigma, double cc0, double cc1, double cc2, double cc3, double cc4, double cc5,
                             double cc6, double cc7, double cc8, double delta, double phi, double theta, double qxn,
                             double qyn, double qzn, double q, int i0, int i1, int i3, int i4,
                             double &sp, double &sf, int n, int &midpntit ) const;


#ifdef __CUDACC__
    __host__ __device__
#endif
    double gammln( double xx ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double polycscube(double rho1, double rho2, double p1, double p2, double alf1, double alf2, double rn, double pf,
                      double sigma, double q) const;


#ifndef USE_CLASSLIB
};
#endif
