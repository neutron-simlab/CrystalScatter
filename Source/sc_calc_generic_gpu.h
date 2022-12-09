#ifndef SASCALC_GENERIC_GPU_H
#define SASCALC_GENERIC_GPU_H

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <string>
#include <map>
#include <list>

#include "sc_math.h"
#include "sc_globalConfig.h"

#ifndef __CUDACC__
#include <QDebug>
#define DCNT(x) x
#else
#define DCNT(x)
#endif



#ifdef __CUDACC__
// GPU Kernel definitions (these functions will run on the GPU)
class SasCalc_GENERIC_calculation;

__global__ void doIntCalc_GENERIC_kernel(SasCalc_GENERIC_calculation);  // Normal calculations

__global__ void doIntFitCalc_GENERIC_kernel(SasCalc_GENERIC_calculation);  // Special 2d-Fit Calculations
#endif


/**
 * @brief The SasCalc_GENERIC_calculation class
 * Class for all FCC calculation routines
 */
class SasCalc_GENERIC_calculation
{
public:
    SasCalc_GENERIC_calculation();

    void cleanup();

    //xx// Markierung für Interface-Settings, die hier nicht genutzt werden.

    // Setting / getting / defaults of general parameter for all routines
    void setRadQ1( bool f ) { if ( f ) calcQuadrants = radQ1; }
    void setRadQ2( bool f ) { if ( f ) calcQuadrants = radQ2; }
    void setRadQ4( bool f ) { if ( f ) calcQuadrants = radQ4; }
    void setExpandImage( bool f ) { bExpandImage = f && calcQuadrants != radQ4; }
    void setGridPoints( int val ) { zmax = val; }
    //xx//void setUvec( Double3 par ) { u=par; }
    //xx//void setVvec( Double3 par ) { v=par; }
    //xx//void setNvec( Double3 par ) { n=par; }
    void setHKLmax( int val ) { params.hklmax = val; }
    void setAx1( Double3 par ) { ax1=par; }
    void setAx2( Double3 par ) { ax2=par; }
    void setAx3( Double3 par ) { ax3=par; }
    void setSigXYZ( Double3 par ) { sig = 4.0 / par; }
    // Rotation parameter
    //xx//void setRotation( Double3 par ) { rot=par; }
    //xx//void setRotAngle( double par ) { rota=par; }
    //xx//void setTiltAng( double par ) { adet=par; }
    //xx//void setxycur( double x, double y ) { xcur=x; ycur=y; xycur=sqrt(xcur*xcur+ycur*ycur); }
    //xx//void setanglecur( double val ) { anglecur = val * M_PI / 180.0; }
    // Special FCC Parameter
    void setQMax( double val ) { qmax = val; }
    void setStretch( double val ) { stretch = val; }
    void setRelDis( double val ) { reldis = val; }
    void setDist( double val ) { dist = val; }
    void setDomainsize( double val ) { params.domainsize = val; }
    void setCheckBoxTwinned( bool f ) { CheckBoxTwinned = f; }
    void setDisplacement( double val ) { displacement=val; dwfactor=val*val/3.0; }  /*Z0311=19436*/
    //xx//void setMaxIter( int i ) { maxIter = i; }
    typedef enum { cbintHomogeneous,    // homogeneous:=true;
                   cbintCoreShell,      // coreshell:=true;
                   cbintCoreInShell,    // coreinshell:=true;
                   cbintMultiShell,     // lipo:=true;
                   cbintMyelin          // myelin:=true;
                 } _cbInterior;         /*Z0311=19828 aus den Kommentaren den Typ entnommen*/
    void setComboBoxInterior( int val ) { ComboBoxInterior = static_cast<_cbInterior>(val); }
    typedef enum { cbpartSphere,        /*Z0311=19787 aus den Kommentaren den Typ entnommen*/
                   cbpartCylinder,
                   cbpartDisk,
                   cbpartVesicle,
                   cbpartCube,
                   cbpartEllipsoide,
                   cbpartTriaxEllips,   // triaxial ellipsoid
                   cbpartSuperEllips    // super ellipsoid, barrel
                 } _cbParticle;
    void setComboBoxParticle( int val ) { ComboBoxParticle = static_cast<_cbParticle>(val); }

    //{NV}-10186
    // Es wird von <ComboBoxPeak.ItemIndex> die Variable <shp> gesetzt
    typedef enum { cbpeakLorentzian=1,      // shp:=1;                          Lorentzian:=true;
                   cbpeakGaussian,          // shp:=2;                          Gaussian:=true;
                   cbpeakMod1Lorentzian,    // shp:=3;                          Lorentzian1:=true;
                   cbpeakMod2Lorentzian,    // shp:=4;                          Lorentzian2:=true;
                   cbpeakPseudoVoigt,       // shp:=5; shf:=EditPeakPar.Text;   Voigt:=true;
                   cbpeakPearsonVII,        // shp:=6; shf:=EditPeakPar.Text;   Pearson:=true;
                   cbpeakGamma,             // shp:=7; shf:=EditPeakPar.Text;   BurgerG:=true;
                   cbpeakAnisotropicGaussian// shp:=8;                          -kein Flag-
                 } _cbPeak;
    void setComboBoxPeak( int v ) { shp = static_cast<_cbPeak>(v+1); }            //{NV} (*** Peak Shapes ***)
    void setPeakPar( double p ) { /*shf=p;*/ eta=p/100.; beta=p; /*bnu=p;*/ } //{NV} TODO - unklar
    // ComboBoxLattice=13,19 -> ComboBoxPeak=7 -> shp=8

    void setAzi( double v ) { aziwidth=v; }             //{NV} TODO
    typedef enum { /* 0*/ordis_Gaussian,             /* 1*/ordis_Exponential,      /* 2*/ordis_Onsager,
                   /* 3*/ordis_Maier_Saupe,          /* 4*/ordis_CutOff,           /* 5*/ordis_Laguerre,
                   /* 6*/ordis_ZDir,                 /* 7*/ordis_Isotropic,        /* 8*/ordis_mirrored_Gaussian,
                   /* 9*/ordis_mirrored_Exponential, /*10*/ordis_mirrored_Onsager, /*11*/ordis_mirrored_Maier_Saupe,
                   /*12*/ordis_mirrored_CutOff,      /*13*/ordis_FiberPattern
                 } _ordistype;
    void setOrdis( int o ) { ordis = static_cast<_ordistype>(o); }
    //{NV} - unit cell definiton
    void setUCA( double v ) { params.uca = v; }
    void setUCB( double v ) { params.ucb = v; }
    void setUCC( double v ) { params.ucc = v; }
    void setUCalpha( double v ) { params.alpha_deg = v; }
    void setUCbeta( double v ) { params.beta_deg = v; }
    void setUCgamma( double v ) { params.gamma_deg = v; }
    void setUCpsi( double v ) { ucpsi = v; }
    void setUCn1( double v ) { ucn1 = v; }
    void setUCn2( double v ) { ucn2 = v; }
    void setUCn3( double v ) { ucn3 = v; }
    void setCheckBoxWAXS( bool v ) { CheckBoxWAXS = v; }
    void setTheta( double v ) { theta = v; }
    void setPhi( double v ) { phi = v; }    /*TPV*/
    //xx//void setxrdalf( double v ) { xrdalf=v*M_PI/180.; } // StrToFloat(EditWAXSangle.Text);    (* [deg] *)
    void setwave( double val ) { params.wavelength=val;  wave=val; /*wave kann geändert werden*/ }

    void setCheckBox10( bool v ) { CheckBox10 = v; }

    void setpixnox( double v ) { pixnox=v; }
    void setpixnoy( double v ) { pixnoy=v; }
    void setpixx( double v ) { pixx=v; }
    void setpixy( double v ) { pixy=v; }
    void setdet( double val ) { det = val; }

    //xx//void setCheckBoxf2q( bool f ) { CheckBoxf2q = f; }
    void setRadioButtonDebyeScherrer( bool f ) { RadioButtonDebyeScherrer = f; }
    void setRadioButtonPara( bool f ) { RadioButtonPara = f; }
    //xx//void setRadioButtonSolid( bool f ) { RadioButtonSolid = f; }
    //xx//void setRadioButtonCS( bool f ) { RadioButtonCS = f; }
    //xx//void setRadioButtonCHS( bool f ) { RadioButtonCHS = f; }

    //void setLatticeF( double val ) { params.latt_a = val; }
    //void setLatticeCF( double val ) { params.latt_c = val; }
    void setCeffF( double val ) { params.ceff = val; }
    void setSigmaF( double val ) { params.sigma = val; }
    void setBFactorF( double val ) { bfactor = val; }
    void setRadiusF( double val ) { params.radius = val; }  /*TPV*/
    void setRadiusI( double val ) { params.radiusi = val; }
    void setDBetaF( double db ) { params.dbeta = db /* * M_PI / 180.*/ ; }
    void setRho( double val ) { params.rho = val; } // neu

    void setIso( double v ) { iso=v; }
    void setIZero( double v ) { izero=v; }  /*TPV*/
    void setBase( double v ) { base=v; }    /*TPV*/
    void setIFluc( double v ) { ifluc=v; }
    void setRFluc( double v ) { rfluc=v; }

    void setP1( double v ) { params.p1=v; }
    void setSigmaL( double v ) { params.sigmal=v; } // Wird in polyliposome() von reff überschrieben!
    void setLength( double v ) { params.length=v; }     /*TPV*/
    void setShellNo( double v ) { params.shellno=v; }
    void setReff( double v ) { params.reff=v; }
    void setAcpl( double v ) { acpl=v; }
    void setBcpl( double v ) { bcpl=v; }
    void setPor( double v ) { por=v; }

    void setLType( int v ) { ltype=v; }

    void setBeamStop( double x, double y ) { beamX0=x; beamY0=y; }

    void setLattPar( int cols, int rows, double centerx, double centery,
                    double wavelen, double distance, double pixx, double pixy )
    {
        params.latFit.cols = cols;
        params.latFit.rows = rows;
        params.latFit.centerx = centerx;
        params.latFit.centery = centery;
        params.latFit.wavelen = wavelen;
        params.latFit.distance = distance;
        params.latFit.pixx = pixx;
        params.latFit.pixy = pixy;
    }

    // Run the calculation once
    void doCalculation( int numThreads, progressAndAbort pa );
    double doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt);

    //QString tpvPerformRandom(QStringList ids);
    std::string tpvPerformRandom(std::list<std::string> ids);

    void setNoFitRect( int id, int x0, int y0, int x1, int y1 )
    {
        if ( id<0 || id>=4 ) return;
        noFitX0[id] = x0;
        noFitX1[id] = x1;
        noFitY0[id] = y0;
        noFitY1[id] = y1;
    }

    double higResTimerElapsedCalc, higResTimerElapsedPrep;

private:
    progressAndAbort progAndAbort;

    static SasCalc_GENERIC_calculation *inst;

    typedef struct
    {
        double *ptr;
        double oldval;
        enum { norm, uca, phi, azi, dw } flag;
    } tpvRandomHelper;
    //QHash< QString/*name*/, tpvRandomHelper* > tpvRandomOldValues;
    std::map< std::string, tpvRandomHelper* > tpvRandomOldValues;

    _radQX calcQuadrants;
    bool bExpandImage;
    bool CheckBoxTwinned;
    bool CheckBoxWAXS;
    bool CheckBox10;
    int zmax;
    double dwfactor, qmax, stretch, displacement, bfactor; //, dbeta; -> params.dbeta
    //xx//Double3 u, v, n;
    Double3 ax1, ax2, ax3, sig;

    //xx//Double3 rot;            //{NV} (in) rotation axis
    //xx//double rota, rota0;     //{NV} (in) rotation angle
    //xx//double rotnorm;         //{NV} (calc)
    //xx//double adet;            //{NV} (in) tilt angle
    //xx//double xcur, ycur;      //{NV} (in) qx,qy of cursor
    //xx//double xycur;           //{NV} (calc)
    //xx//double anglecur;        //{NV} (in) angle of cursor

    double beamX0, beamY0;  // neu - Beamstop

    _cbInterior ComboBoxInterior;   //{NV}
    _cbParticle ComboBoxParticle;   //{NV}
    int  cs;                        //{NV}
    bool RadioButtonDebyeScherrer, RadioButtonPara;
         //xx//CheckBoxf2q, RadioButtonSolid, RadioButtonCS, RadioButtonCHS;

    _ordistype ordis;
    //xx//int maxIter;

    int  ltype; // Zur Auswahl der Operation in buttonHKLclick und auch später für die Berechnungsmethoden

    bool prepareCalculation();

    // Variables for Lattice3D
    double /*sigmad,*/ c0, c1, c2, c3, c4; //, p1, approx, ifl, rfl, bf;
    //TODO: diese Werte müssen noch übergeben/vorbesetzt werden
    double azidom, eta, dom1, dom2, dom3;
    double beta, aziwidth, /*critdist,*/ /*shf,*/ alphash, zz;  //{NV}
    double theta, phi, /*bnu,*/ /*xrdalf,*/ wave, iso, izero, base, ifluc, rfluc;
    _cbPeak shp;                                        //{NV}
    int cdim, orcase, partdim, part;                          //{NV}
    double order, norm, limq1, limq2, limq3;            //{NV}
    double limq4, limq5, limq6;
    //20220311: TODO
    double limq7, limq8, limq9, limq1f, limq2f, limq3f, limq4f, limq5f, limq6f, limq7f, limq8f, limq9f;

    double F121,F122,F123,pqr0, sinphic, cosphic;
    double ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10;
    double ccc1,ccc2,ccc3,vv,vv3;
    double cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10;

    // TODO: Neue Parameter
    double acpl, bcpl, por;
    bool lat1d, lat2d, lat3d;

    double fhkl, qhkl, qhkl0;
    double qxhkl, qyhkl, qzhkl, qxyhkl;
    double qxhklt, qyhklt, qzhklt;
    int    mhkl, h, k, l;
    //double mildist, milq, angles, qxs, qys, norms, milqs, hs, ks, ls;
    //   werden gesetzt/berechnet, um in die GUI geschrieben zu werden, was hier nicht genutzt wird
    double peaknorm1, peaknorm2, g3, x2phihkl;

    //{NV} - unit cell definiton (TODO: Double3?)
    double /*cosa, cosb, cosg,*/ ucpsi,
           ucn1, ucn2, ucn3, /*ucn, ucn1n, ucn2n, ucn3n,*/
           ucphi, uctheta, cosphi, sinphi, costheta, sintheta,
           ucl1,ucl2,ucl3, pixnox,pixx, pixnoy,pixy, det;
    double ri11,ri12,ri13,ri21,ri22,ri23,ri31,ri32,ri33,
           rt11,rt12,rt13,rt21,rt22,rt23,rt31,rt32,rt33,
           nuvwx,nuvwy,nuvwz,uuvwx,uuvwy,uuvwz,vuvwx,vuvwy,vuvwz,
           nhklx,nhkly,nhklz,uhklx,uhkly,uhklz,vhklx,vhkly,vhklz;
    int   *latpar1ptr/* 6*/, *latpar2ptr/* 6*/;
    float *latpar3ptr/*14*/; //, *latpar4ptr/* 2*/; // : Array[1..5000,0..15] of real;
    //  latpar4 wird (noch) nicht wirklich gebraucht
    double latpar[4];
    int peakmax1, peakmax2;
    size_t latpar1Size, latpar2Size, latpar3Size;


    // Weitere Methodendefinitionen ausgelagert
#define USE_CLASSLIB
#include "sc_libs_gpu.h"
#include "sc_memory_gpu.h"


    int   latpar1( int a, int b ) const { return latpar1ptr[latparIDX(a,b, 6)]; }      // [5000][6], genutzt: 0,1,2,3,4,5
    int   latpar2( int a, int b ) const { return latpar2ptr[latparIDX(a,b, 6)]; }      // [5000][6], genutzt: 0,1,2,3,4,5
    float latpar3( int a, int b ) const { return latpar3ptr[latparIDX(a,b,14)]; }      // [5000][15], genutzt: 1 bis 12
    //#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]      // [5000][15], genutzt: nichts

    void setLatpar1( int a, int b, int v ) { latpar1ptr[latparIDX(a,b, 6)] = v; }      // [5000][6], genutzt: 0,1,2,3,4,5
    void setLatpar2( int a, int b, int v ) { latpar2ptr[latparIDX(a,b, 6)] = v; }      // [5000][6], genutzt: 0,1,2,3,4,5
    void setLatpar3( int a, int b, float v ) { latpar3ptr[latparIDX(a,b,14)] = v; }    // [5000][15], genutzt: 1 bis 13
    //#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]      // [5000][15], genutzt: nichts

    // 20220311: TODO: neue Arrays
    //ArrayImax2D carr11pm, carr22pm;  -> params

    // Internal variables outside all loops
    int zzmin, zzmax, iimin, iimax;
    double ax1n_sigx, ax2n_sigy, ax3n_sigz, cubevol, phiwidth, reldis, dist;

    int noFitX0[4], noFitX1[4], noFitY0[4], noFitY1[4];
    int fitBStopPixel, fitBorderPixel;

    // Thread parts without GPU
    static void *doThreadCalculation(void *arg);
    static void *doThreadFitCalculation(void *arg);
    pthread_t *threads;
    int *thread_args;
    int numberOfThreads;

    bool _endThread;

    void doIntCalc_GENERIC(int ihex);  /* [ 5] 11600 FCC Spheres = Face centered cubic */
    void doIntFitCalc_GENERIC(int x);

#ifdef __CUDACC__
    friend __global__ void doIntCalc_GENERIC_kernel(SasCalc_GENERIC_calculation);
    friend __global__ void doIntFitCalc_GENERIC_kernel(SasCalc_GENERIC_calculation);
#endif

#ifdef __CUDACC__
    __host__ __device__
#endif
    static void doIntCalc_GENERIC_F(const SasCalc_GENERIC_calculation& FCC, int ihex, int i);

#ifdef __CUDACC__
    __host__ __device__
#endif
    static void doIntFitCalc_GENERIC_F(const SasCalc_GENERIC_calculation& FCC, int x, int y);

#ifdef __CUDACC__
    __host__ __device__
#endif
    static double doIntCalc_GENERIC_q_xyz(const SasCalc_GENERIC_calculation& FCC, double qx, double qy, double qz);


    // Die folgenden Routinen werden nur beim Host (prepareCalculations) benötigt

    void corotations( double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg,
                      double u, double v, double w, double ephi,
                      bool lat1d, bool lat2d, bool lat3d, // Neu 20220311
                      double &m11, double &m12, double &m13, double &m21, double &m22, double &m23, double &m31, double &m32, double &m33,
                      double &mtw11, double &mtw12, double &mtw13, double &mtw21, double &mtw22, double &mtw23, double &mtw31, double &mtw32, double &mtw33, double &vvol,
                      double &nuvwx, double &nuvwy, double &nuvwz, double &uuvwx, double &uuvwy, double &uuvwz, double &vuvwx, double &vuvwy, double &vuvwz,
                      double &nhklx, double &nhkly, double &nhklz, double &uhklx, double &uhkly, double &uhklz, double &vhklx, double &vhkly, double &vhklz );

    void coefficients( double l, double r, double rm, double sigmal, double sigma,
                      double epsi, double alfa, double dbeta, double theta, double phi,
                      int part, int dim, int nmax, int ordis, int cs,
                      double *myarray /*Z0311=9310*/,
                      int &cho1,
                      double &por, double &order, double &norm, double &lim1, double &lim2,
                      double &lim3, double &lim4, double &lim5, double &lim6, double &lim7,
                      double &lim8, double &lim9,
                      double &lim1f, double &lim2f, double &lim3f, double &lim4f, double &lim5f,
                      double &lim6f, double &lim7f, double &lim8f, double &lim9f,
                      _carrXX *CR );
    //double *carr1p, double *carr2p, double *carr3p,
    //double *carr4p, double *carr5p, double *carr6p,
    //double *carr7p, double *carr8p, double *carr9p, //: CoeffArrayType  /*Z0311=9314*/,
    //double *carr1f, double *carr2f, double *carr3f,
    //double *carr4f, double *carr5f, double *carr6f,
    //double *carr7f, double *carr8f, double *carr9f //: CoeffArrayType  /*Z0311=9315*/
    // /* var carr1pm,carr2pm: ArrayImax1D; */  /*Z0311=9316*/
    // /*ArrayImax2D &carr11pm, ArrayImax2D &carr22pm*/ ); //: ArrayImax2D )

    static long dbgCount;  // Zum Test...

};

#endif // SASCALC_GENERIC_GPU_H
