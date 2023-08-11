
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


public:
    void endThread();

#define  coeffarray_len  (150+1)
#define  imax2d_len  (130+1)
    typedef double CoeffArrayType[coeffarray_len]; // =array[0..150] of extended;
    typedef double ArrayImax2D[imax2d_len][imax2d_len];
    //typedef double ArrayImax1D[20000+1];     // TODO ist das so groß nötig?

    // 20220311: Die vielen kleinen Arrays fasse ich zu einem Record zusammen
    typedef struct
    {
        double carr1p[coeffarray_len], carr2p[coeffarray_len], carr3p[coeffarray_len],
               carr4p[coeffarray_len], carr5p[coeffarray_len], carr6p[coeffarray_len],
               carr7p[coeffarray_len], carr8p[coeffarray_len], carr9p[coeffarray_len],
               carr1f[coeffarray_len], carr2f[coeffarray_len], carr3f[coeffarray_len],
               carr4f[coeffarray_len], carr5f[coeffarray_len], carr6f[coeffarray_len],
               carr7f[coeffarray_len], carr8f[coeffarray_len], carr9f[coeffarray_len];

        double myarray[20]; // TODO: Werte besser in globale Variablen, der größte Index ist 17
        // aus prepareCalculation:
        //params.CR->myarray[ 0] = params.length;     /*  axon length  */  //Z=20162
        //params.CR->myarray[ 1] = params.radius;     /*  axon radius  */  //Z=20163
        //params.CR->myarray[ 2] = params.sigma;      //Z=20164
        //params.CR->myarray[ 3] = params.sigmal;     //Z=20165
        //params.CR->myarray[ 4] = params.radiusi;    /*  no. of layers  */  //Z=20166
        //params.CR->myarray[ 5] = params.alphash;           /*  l_ex  */  //Z=20167
        //params.CR->myarray[ 6] = params.rho;        /*  l_in  */  //Z=20168
        //params.CR->myarray[ 7] = acpl;              /*  l_lip_head  */  //Z=20169
        //params.CR->myarray[ 8] = bcpl;              /*  l_lip_tail  */  //Z=20170
        //params.CR->myarray[ 9] = params.uca;        /*  phi_axon  */  //Z=20171
        //params.CR->myarray[10] = params.ucb;        /*  phi_intra  */  //Z=20172
        //params.CR->myarray[11] = params.ucc;        /*  phi_extra  */  //Z=20173
        //params.CR->myarray[12] = params.domainsize; /*  phi_head  */  //Z=20174
        //params.CR->myarray[13] = aziwidth;          /*  phi_tail  */  //Z=20175
        //params.CR->myarray[14] = 1;                 /*  inmax  */  //Z=20176
        //params.CR->myarray[15] = 1;                 /*  vv  */  //Z=20177
        //params.CR->myarray[16] = 1;                 /*  rmax  */  //Z=20178
        //params.CR->myarray[17] = iso;               /*  test parameter  */  //Z=20179

        // Array-Parameter für coefficients(), formpq(), formfq()
        ArrayImax2D carr11pm, carr22pm;

        int latparMaxCheckCount[3];

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
        double alphash1;    // 'Alpha'
        double alpha;       // 'RotAlpha'
        double ceff;        // often equal with "TwRatio" -> TODO
        double reff;        // allways const 10.0 (no input in GUI) -> TODO?
        double uca; // latt_a;      // a
        double ucb; // latt_b;      // b
        double ucc; // latt_c;      // c
        double domainsize;
        double width_zuf; // Zufalls-Faktor für domainsize/width
        double aziwidth;
        double dbeta;
        double psphere_r, psphere_z;  // helper var for psphere

        double norm, por; // por ist ein Output von Coefficients()
        int    part, cs, orcase;
        double polTheta, polPhi;
        double limq1, limq2, limq3, limq4, limq5, limq6, limq7, limq8, limq9;
        double limq1f, limq2f, limq3f, limq4f, limq5f, limq6f, limq7f, limq8f, limq9f;

        // Globale Parameter für ButtonHKLClick()
        int    hklmax;
        double ucalpha_deg, ucbeta_deg, ucgamma_deg;
        //double amax, bmax, cmax;

        double p11,p12,p13,p21,p22,p23,p31,p32,p33;

        _carrXX  *CR;

        Double3 ax1, ax2, ax3, sig;

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

    static const int latparlen = 2500;  // war 5000, latpar Max Check: 342, 2196, 1098
    #define latparIDX(a,b,l) (a*l+b)
    // Hilfsdaten zur Ermittlung der maximal genutzten Arraylänge, um latparlen zu minimieren
    inline void latparMaxCheckInit() { params.CR->latparMaxCheckCount[0]=0; params.CR->latparMaxCheckCount[1]=0; params.CR->latparMaxCheckCount[2]=0; }
    inline void latparMaxCheck( const int i, const int a ) const
    {   if ( params.CR->latparMaxCheckCount[i]<a ) params.CR->latparMaxCheckCount[i]=a;
#ifndef __CUDACC__
        if ( a >= latparlen ) qDebug() << "##################################################" << a << latparlen;
#endif
    }

#define latpar1(a,b) latpar1ptr[latparIDX(a,b, 6)]      // [5000][6], genutzt: 0,1,2,3,4,5
#define latpar2(a,b) latpar2ptr[latparIDX(a,b, 6)]      // [5000][6], genutzt: 0,1,2,3,4,5
#define latpar3(a,b) latpar3ptr[latparIDX(a,b,17)]      // [5000][15], genutzt: 1 bis 16
    //#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]      // [5000][15], genutzt: nichts

#define setLatpar1(a,b,v) { latparMaxCheck(0,a); latpar1ptr[latparIDX(a,b, 6)] = v; }     // [5000][6], genutzt: 0,1,2,3,4,5
#define setLatpar2(a,b,v) { latparMaxCheck(1,a); latpar2ptr[latparIDX(a,b, 6)] = v; }     // [5000][6], genutzt: 0,1,2,3,4,5
#define setLatpar3(a,b,v) { latparMaxCheck(2,a); latpar3ptr[latparIDX(a,b,17)] = v; }     // [5000][15], genutzt: 1 bis 16
    //#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]        // [5000][15], genutzt: nichts


    static const int np=5;
    static const int jmaxp=21;
    typedef double RealArrayNP[np+2]; // =array[1..np] of extended;
    typedef double RealArrayJMAXP[jmaxp+2]; // =array[1..jmaxp] of extended;
    // Die Arrays zur Sicherheit eins länger, da in den Routinen oft vor dem
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
    double gamma(double z) const;

#ifdef __CUDACC__
__host__ __device__
#endif
    double polyvesicle(double q) const;

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
    void polint( double *xa/*RealArrayNP*/, double *ya/*RealArrayNP*/,
                 int n, double x, double &y, double &dy, const char *dbg ) const;

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


    void ButtonHKLClick(int ltype) const;

    void fhkl_c( int lat, int h, int k, int l,
                 double &sphno, double &fhkl, double &qhkl, double &qhkl0 ) const;

    void extinction( int lat, int h, int k, int l, int aniso,
                     int &mhkl, double &fhkl ) const;



#ifdef __CUDACC__
__host__ __device__
#endif
void qrombdeltac( double p1, double sigma, double alfa,
                  double theta, double phi, double qx, double qy, double qz,
                  double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n,
                  double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y,
                  double ax2z, double ax3x, double ax3y, double ax3z, double sigx, double sigy,
                  double sigz, int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                  double *carr1, double &pq ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
//void trapzddeltac( double a, double b, double l, double r, double dbeta, double theta, double phi,
//                   double qx, double qy, double qz, double p11, double p12, double p13, double p21,
//                   double p22, double p23, double p31, double p32, double p33,
//                   double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
//                   double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
//                   double ax3y, double ax3z, double sigx, double sigy, double sigz,
//                   int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
//                   double *carr1,
//                   double &pq, int n, int &trapzdit ) const;
void trapzddeltac( double a, double b, double l, double r, double p1, double sigma, double alfa,
             double dbeta, double theta, double phi, double qx, double qy, double qz,
             double p11, double p12, double p13, double p21, double p22, double p23,
             double p31, double p32, double p33, double qxn, double qyn, double qzn,
             double qhkl, double ax1n, double ax2n, double ax3n, double ax1x, double ax1y,
             double ax1z, double ax2x, double ax2y, double ax2z, double ax3x, double ax3y,
             double ax3z, double sigx, double sigy, double sigz,
             int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
             double *carr1, double &pq, int n, int &trapzddeltac_cnt ) const;

#ifdef __CUDACC__
__host__ __device__
#endif
//    void trapzdchid( double a, double b, double l, double r, double sigma, double dbeta, double delta, double phi, double,
//                     double qx, double qy, double qz, double p11, double p12, double p13, double p21, double p22, double p23,
//                     double p31, double, double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n,
//                     double ax2n, double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
//                     double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz, int dim, int i0, int, int i1, int, int i3, int i4, double *carr1, double &pq, int n,
//                     int &trapzdit ) const;
//void trapzdchid_OLDVERSION( double a, double b, double r, double sigma, double dbeta, double delta, double theta, double phi,
//                 double qx, double qy, double qz, double p11, double p12, double p13, double p21, double p22, double p23,
//                 double p31, double p32, double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n,
//                 double ax2n, double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
//                 double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
//                 int ordis, int dim, int i0, int i1, int i2, int i3, int i4, double *carr1, double &pq, int n,
//                 int &trapzdit ) const;
void trapzdchid( double a, double b, double l, double r, double p1, double sigma, double alfa,
                   double delta, double theta, double phi, double qx, double qy, double qz, double p11,
                   double p12, double p13, double p21, double p22, double p23, double p31, double p32,
                   double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n,
                   double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
                   double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
                   int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                   double *carr1, double &pq, int n,
                   int &trapzdchid_cnt ) const;


#ifdef __CUDACC__
__host__ __device__
#endif
    //void qrombchid( double l, double r, double sigma, double dbeta, double delta, double theta, double phi,
    //                double qx, double qy, double qz, double p11, double p12, double p13, double p21,
    //                double p22, double p23, double p31, double p32, double p33,
    //                double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
    //                double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
    //                double ax3y, double ax3z, double sigx, double sigy, double sigz,
    //                int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
    //                double *carr1,
    //                double &pq ) const;
    void qrombchid( double l, double r, double p1, double sigma, double alfa, double delta,
              double theta, double phi, double qx, double qy, double qz, double p11, double p12, double p13,
              double p21, double p22, double p23, double p31, double p32, double p33, double qxn, double qyn,
              double qzn, double qhkl, double ax1n, double ax2n, double ax3n, double ax1x, double ax1y,
              double ax1z, double ax2x, double ax2y, double ax2z, double ax3x, double ax3y, double ax3z,
              double sigx, double sigy, double sigz,
              int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
              double *carr1, double &pq ) const;


#ifdef __CUDACC__
__host__ __device__
#endif
//    double formfq(double radius, double zz, double limq1, double qx, double qy,
//                  int part, int cs,
//                  double *carr3p /*: CoeffArrayType*/ ) const;
double formfq( double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const;


#ifdef __CUDACC__
__host__ __device__
#endif
//    double formpq(double length, double radius, double p1, double rho, double zz, double qlim, double limq1, double limq2, double limq3, double limq4, double qx, double qy, double ql, double q, double norm,
//                  int part, int cs, int ordis, int orcase,
//                  double *carr1p, double *carr2p, double *carr3p, double *carr4p, /*: CoeffArrayType;*/
//                  double *carr5p /*: ArrayImax1D*/ ) const;
    double formpq( double sigmal, double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double gammln( double xx ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double polycscube(double rho1, double rho2, double p1, double p2, double alf1, double alf2, double rn, double pf,
                      double sigma, double q) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double spy( double q ) const;
