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
 * Class for all calculation routines
 */
class SasCalc_GENERIC_calculation
{
public:
    SasCalc_GENERIC_calculation();

    void cleanup();

    // Neues Interface (April 2023), Parameterwerte gruppiert
    void setLatticeParams(int ltype, double uca, double ucb, double ucc,
                          double ucalpha, double ucbeta, double ucgamma,
                          double editceff, bool twinned, double ceffcyl,
                          double reff, double acpl, double bcpl,
                          double ifluc, double rfluc)
    {
        this->ltype = ltype;
        params.uca = uca;
        params.ucb = ucb;
        params.ucc = ucc;
        params.alpha_deg = ucalpha;
        params.beta_deg = ucbeta;
        params.gamma_deg = ucgamma;
        params.ceff = editceff;
        CheckBoxTwinned = twinned;
        conc = epsilon = ceffcyl;  //StrToFloat(EditCeffcyl.Text);
        params.reff = reff;
        this->acpl = acpl;
        this->bcpl = bcpl;
        this->ifluc = ifluc;
        this->rfluc = rfluc;
    }

    void setCalculationParams(int gridp, int hklmax, double qmax, bool radq1,
                              bool radq2, bool radq4, bool expand)
    {
        zmax = gridp;
        params.hklmax = hmax = kmax = lmax = hklmax;
        this->qmax = qmax;
        if ( radq1 ) calcQuadrants = radQ1;
        if ( radq2 ) calcQuadrants = radQ2;
        if ( radq4 ) calcQuadrants = radQ4;
        bExpandImage = expand && (calcQuadrants != radQ4);
    }

    void setParticleParams(int cbparticle, int cbordis, int cbinterior,
                           double radius, double radiusi, double sigma,
                           double dbeta, double length, double sigmal,
                           double alpha, double rho)
    {
        ComboBoxParticle = static_cast<_cbParticle>(cbparticle);
        ordis = static_cast<_ordistype>(cbordis);
        ComboBoxInterior = static_cast<_cbInterior>(cbinterior);
        params.radius = radius;
        params.radiusi = radiusi;
        params.sigma = sigma;
        params.dbeta = dbeta;
        params.length = length;
        params.sigmal = sigmal; // Wird in polyliposome() von reff überschrieben!
        this->alphash = alpha;
        params.rho = rho;
    }

    void setPeakShapeParams(int cbpeak, double peakpar, double debyewaller,
                            double azi, double domainsize, bool debyescherrer,
                            bool para)
    {
        shp = static_cast<_cbPeak>(cbpeak+1);            // (*** Peak Shapes ***)
        /*shf=p;*/ eta=peakpar/100.; beta=peakpar; /*bnu=p;*/  // TODO - unklar
        displacement=debyewaller; dwfactor=debyewaller*debyewaller/3.0;   /*Z0311=19436*/
        aziwidth = azi;
        params.domainsize = domainsize;
        RadioButtonDebyeScherrer = debyescherrer;
        RadioButtonPara = para;
    }

    void setPeakShapeMatrix(double ax1, double ax2, double ax3,
                            double ay1, double ay2, double ay3,
                            double az1, double az2, double az3,
                            double sigx, double sigy, double sigz)
    {
        params.ax1 = Double3(ax1,ay1,az1);
        params.ax2 = Double3(ax2,ay2,az2);
        params.ax3 = Double3(ax3,ay3,az3);
        params.sig = Double3(4.0/sigx,4.0/sigy,4.0/sigz);
    }

    void setOrientationParams(double ucpsi, double ucn1, double ucn2, double ucn3,
                              double theta, double phi, double rottheta, double rotphi)
    {
        this->ucpsi = ucpsi;
        this->ucn1 = ucn1;
        this->ucn2 = ucn2;
        this->ucn3 = ucn3;
        polTheta = theta;
        polPhi = phi;
        rotaxTheta = rottheta;
        rotaxPhi = rotphi;
    }

    void setControlParams(double bfactor, /*double rotx, double roty,*/ bool waxs, double p1)
    {
        this->bfactor = bfactor;
        CheckBoxWAXS = waxs;
        params.p1 = p1;
    }

    void setPixelManipulationParams(double iso, double i0, double base)
    {
        this->iso = iso;
        izero = i0;
        this->base = base;
    }

    void setExperimentParams(double pixelnox, double pixelnoy, double pixelsizex,
                             double pixelsizey, double detdist, double wavelength,
                             double beamcx, double beamcy)
    {
        this->pixnox = pixelnox;
        this->pixnoy = pixelnoy;
        this->pixx = pixelsizex;
        this->pixy = pixelsizey;
        this->det = detdist;
        params.wavelength=wavelength;  this->wave=wavelength; /*wave kann geändert werden*/
        this->beamX0 = beamcx; this->beamY0 = beamcy;
    }


    // Altes Interface (jeder Wert einzeln)

    // Setting / getting / defaults of general parameter for all routines
    void setRadQ1( bool f ) { if ( f ) calcQuadrants = radQ1; }//
    void setRadQ2( bool f ) { if ( f ) calcQuadrants = radQ2; }//
    void setRadQ4( bool f ) { if ( f ) calcQuadrants = radQ4; }//
    void setExpandImage( bool f ) { bExpandImage = f && calcQuadrants != radQ4; }//
    void setGridPoints( int val ) { zmax = val; }//
    void setHKLmax( int val ) { params.hklmax = hmax = kmax = lmax = val; }//
    void setAx1( Double3 par ) { params.ax1=par; }//
    void setAx2( Double3 par ) { params.ax2=par; }//
    void setAx3( Double3 par ) { params.ax3=par; }//
    void setSigXYZ( Double3 par ) { params.sig = 4.0 / par; }//
    void setQMax( double val ) { qmax = val; }//
    //void setRelDis( double val ) { reldis = val; }
    //void setDist( double val ) { dist = val; }
    void setDomainsize( double val ) { params.domainsize = val; }//
    void setCheckBoxTwinned( bool f ) { CheckBoxTwinned = f; }//
    void setDisplacement( double val ) { displacement=val; dwfactor=val*val/3.0; }//  /*Z0311=19436*/
    typedef enum { cbintHomogeneous,    // homogeneous:=true;
                   cbintCoreShell,      // coreshell:=true;
                   cbintCoreInShell,    // coreinshell:=true;
                   cbintMultiShell,     // lipo:=true;
                   cbintMyelin          // myelin:=true;
                 } _cbInterior;         /*Z0311=19828 aus den Kommentaren den Typ entnommen*/
    void setComboBoxInterior( int val ) { ComboBoxInterior = static_cast<_cbInterior>(val); }//
    typedef enum { cbpartSphere,        // sphere = 0,partsphere
                   cbpartCylinder,      // cylinder = 1,partcylinder
                   cbpartDisk,          // disk = 2,partdisk
                   cbpartVesicle,       // vesicle = 3,partvesicle
                   cbpartCube,          // cube = 4,partcube
                   cbpartEllipsoide,    // ellipsoid = 5,partellipsoid
                   cbpartTriaxEllips,   // triaxial ellipsoid = 6,parttriellipsoid
                   cbpartSuperEllips,   // super ellipsoid, barrel = 7,partbarrel
                   cbpartSuperball,     // superball = 8,partball
                   cbpartChain,         // excluded volume chain = 9,partchain
                   cbpartkpchain        // Kratky Porod chain = 10,partkpchain
                 } _cbParticle;
    void setComboBoxParticle( int val ) { ComboBoxParticle = static_cast<_cbParticle>(val); }//

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
    void setComboBoxPeak( int v ) { shp = static_cast<_cbPeak>(v+1); }//            // (*** Peak Shapes ***)
    void setPeakPar( double p ) { /*shf=p;*/ eta=p/100.; beta=p; /*bnu=p;*/ }//  // TODO - unklar
    // ComboBoxLattice=13,19 -> ComboBoxPeak=7 -> shp=8

    void setAzi( double v ) { aziwidth=v; }//             // TODO
    typedef enum { /* 0*/ordis_Gaussian,             /* 1*/ordis_Exponential,      /* 2*/ordis_Onsager,
                   /* 3*/ordis_Maier_Saupe,          /* 4*/ordis_CutOff,           /* 5*/ordis_Laguerre,
                   /* 6*/ordis_ZDir,                 /* 7*/ordis_Isotropic,        /* 8*/ordis_mirrored_Gaussian,
                   /* 9*/ordis_mirrored_Exponential, /*10*/ordis_mirrored_Onsager, /*11*/ordis_mirrored_Maier_Saupe,
                   /*12*/ordis_mirrored_CutOff,      /*13*/ordis_FiberPattern
                 } _ordistype;
    void setOrdis( int o ) { ordis = static_cast<_ordistype>(o); }//
    // unit cell definiton
    void setUCA( double v ) { params.uca = v; }//
    void setUCB( double v ) { params.ucb = v; }//
    void setUCC( double v ) { params.ucc = v; }//
    void setUCalpha( double v ) { params.alpha_deg = v; }//
    void setUCbeta( double v ) { params.beta_deg = v; }//
    void setUCgamma( double v ) { params.gamma_deg = v; }//
    void setUCpsi( double v ) { ucpsi = v; }//
    void setUCn1( double v ) { ucn1 = v; }//
    void setUCn2( double v ) { ucn2 = v; }//
    void setUCn3( double v ) { ucn3 = v; }//
    void setCheckBoxWAXS( bool v ) { CheckBoxWAXS = v; }//
    void setPolTheta( double v ) { polTheta = v; }//
    void setPolPhi( double v ) { polPhi = v; }//    /*TPV*/
    void setRotTheta( double v ) { rotaxTheta = v; }//
    void setRotPhi( double v ) { rotaxPhi = v; }//  /*TPV?*/
    void setwave( double val ) { params.wavelength=val;  wave=val; /*wave kann geändert werden*/ }//

    void setAlphash( double v ) { alphash = v; }//

    void setpixnox( double v ) { pixnox=v; }//
    void setpixnoy( double v ) { pixnoy=v; }//
    void setpixx( double v ) { pixx=v; }//
    void setpixy( double v ) { pixy=v; }//
    void setdet( double val ) { det = val; }//

    void setRadioButtonDebyeScherrer( bool f ) { RadioButtonDebyeScherrer = f; }//
    void setRadioButtonPara( bool f ) { RadioButtonPara = f; }//

    //void setLatticeF( double val ) { params.latt_a = val; }
    //void setLatticeCF( double val ) { params.latt_c = val; }
    void setCeffF( double val ) { params.ceff = val; }//
    void setSigmaF( double val ) { params.sigma = val; }//
    void setBFactorF( double val ) { bfactor = val; }//
    void setRadiusF( double val ) { params.radius = val; }//  /*TPV*/
    void setRadiusI( double val ) { params.radiusi = val; }//
    void setDBetaF( double db ) { params.dbeta = db /* * M_PI / 180.*/ ; }//
    void setRho( double val ) { params.rho = val; }//  // neu

    void setCeffCyl( double val ) { conc=epsilon=val; }//  //StrToFloat(EditCeffcyl.Text);

    void setIso( double v ) { iso=v; }//
    void setIZero( double v ) { izero=v; }//  /*TPV*/
    void setBase( double v ) { base=v; }//    /*TPV*/
    void setIFluc( double v ) { ifluc=v; }//
    void setRFluc( double v ) { rfluc=v; }//

    void setP1( double v ) { params.p1=v; }//
    void setSigmaL( double v ) { params.sigmal=v; }// // Wird in polyliposome() von reff überschrieben!
    void setLength( double v ) { params.length=v; }//     /*TPV*/
    void setReff( double v ) { params.reff=v; }//
    void setAcpl( double v ) { acpl=v; }//
    void setBcpl( double v ) { bcpl=v; }//

    void setLType( int v ) { ltype=v; }//

    void setBeamStop( double x, double y ) { beamX0=x; beamY0=y; }//

    // Run the calculation once
    void doCalculation( int numThreads );
    double doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt);

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
    double dwfactor, qmax, displacement, bfactor;
    double beamX0, beamY0;  // Beamstop
    double c0, c1, c2, c3, c4;
    double azidom, eta, dom1, dom2, dom3;
    double beta, aziwidth, critdist, alphash, zz;
    double polTheta, polPhi, rotaxTheta, rotaxPhi, wave, iso, izero, base, ifluc, rfluc;
    double order, norm, limq1, limq2, limq3, limq4, limq5, limq6;
    double limq7, limq8, limq9, limq1f, limq2f, limq3f, limq4f, limq5f, limq6f, limq7f, limq8f, limq9f;
    double F121,F122,F123,pqr0, sinphic, cosphic;
    double ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10;
    double ccc1,ccc2,ccc3,vv,vv3;
    double cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10;
    double acpl, bcpl, por; // por ist ein Output von Coefficients()
    double conc, epsilon;
    double fhkl, qhkl, qhkl0, qxhkl, qyhkl, qzhkl, qxyhkl, qxhklt, qyhklt, qzhklt;
    double peaknorm1, peaknorm2, g3, x2phihkl;
    double cosa, cosb, cosg, ucpsi, ucvol,
           ucn1, ucn2, ucn3, ucn, ucn1n, ucn2n, ucn3n,
           ucphi, uctheta, cosphi, sinphi, costheta, sintheta,
           ucl1,ucl2,ucl3, pixnox,pixx, pixnoy,pixy, det;
    double ri11,ri12,ri13,ri21,ri22,ri23,ri31,ri32,ri33,
           rt11,rt12,rt13,rt21,rt22,rt23,rt31,rt32,rt33,
           nuvwx,nuvwy,nuvwz,uuvwx,uuvwy,uuvwz,vuvwx,vuvwy,vuvwz,
           nhklx,nhkly,nhklz,uhklx,uhkly,uhklz,vhklx,vhkly,vhklz;
    int   *latpar1ptr/* 6*/, *latpar2ptr/* 6*/;
    float *latpar3ptr/*14*/; //, *latpar4ptr/* 2*/; // : Array[1..5000,0..15] of real;
    //  latpar4 wird (noch) nicht wirklich gebraucht
    size_t latpar1Size, latpar2Size, latpar3Size;




    // old order ...
    rvButtonHKLClick retvalHKLClick;//public

private:




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
    int zmax;

    _cbInterior ComboBoxInterior;
    _cbParticle ComboBoxParticle;
    int  cs;
    bool RadioButtonDebyeScherrer, RadioButtonPara;

    _ordistype ordis;

    int  ltype; // Zur Auswahl der Operation in buttonHKLclick und auch später für die Berechnungsmethoden

    bool prepareCalculation();

    // Variables for Lattice3D
    _cbPeak shp;
    int cdim, orcase, partdim, part;

    // internal flags
    bool generic, twin, debyescherrer, paracrystal, quad1, quad2, quad4,
         lattice, lat1d, lat2d, lat3d, Lorentzian, Gaussian, Lorentzian1,
         Lorentzian2, Voigt, Pearson, BurgerG,
         partsphere, partcylinder, partdisk, partvesicle, partcube,
         partellipsoid, parttriellipsoid, partbarrel, partball, partchain,
         partkpchain, homogeneous, coreshell, coreinshell, lipo, myelin;

    int    mhkl, h, k, l, hmax, kmax, lmax;
    //double mildist, milq, angles, qxs, qys, norms, milqs, hs, ks, ls;
    //   werden gesetzt/berechnet, um in die GUI geschrieben zu werden, was hier nicht genutzt wird
    double latpar[4];
    int peakmax1, peakmax2;


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
    double ax1n_sigx, ax2n_sigy, ax3n_sigz, cubevol, phiwidth; //, reldis, dist; Ausgaben

    int noFitX0[4], noFitX1[4], noFitY0[4], noFitY1[4];
    int fitBStopPixel, fitBorderPixel;

    // Thread parts without GPU
    static void *doThreadCalculation(void *arg);
    static void *doThreadFitCalculation(void *arg);
    static pthread_t *threads;
    int *thread_args;
    int numberOfThreads;

    bool _endThread;

#ifdef __CUDACC__
    friend __global__ void doIntCalc_GENERIC_kernel(SasCalc_GENERIC_calculation);
    friend __global__ void doIntFitCalc_GENERIC_kernel(SasCalc_GENERIC_calculation);
#endif

#ifdef __CUDACC__
    __host__ __device__
#endif
    static void doIntCalc_GENERIC_F(const SasCalc_GENERIC_calculation& CALC, int ihex, int i);

#ifdef __CUDACC__
    __host__ __device__
#endif
    static void doIntFitCalc_GENERIC_F(const SasCalc_GENERIC_calculation& CALC, int x, int y);

#ifdef __CUDACC__
    __host__ __device__
#endif
    static double doIntCalc_GENERIC_q_xyz(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz, bool idxCenter, bool beamCenter );


    // Die folgenden Routinen werden nur beim Host (prepareCalculations) benötigt

    void corotations( double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg,
                      double u, double v, double w, double ephi,
                      bool lat1d, bool lat2d, bool lat3d, // Neu 20220311
                      double &m11, double &m12, double &m13, double &m21, double &m22, double &m23, double &m31, double &m32, double &m33,
                      double &mtw11, double &mtw12, double &mtw13, double &mtw21, double &mtw22, double &mtw23, double &mtw31, double &mtw32, double &mtw33, double &vvol,
                      double &nuvwx, double &nuvwy, double &nuvwz, double &uuvwx, double &uuvwy, double &uuvwz, double &vuvwx, double &vuvwy, double &vuvwz,
                      double &nhklx, double &nhkly, double &nhklz, double &uhklx, double &uhkly, double &uhklz, double &vhklx, double &vhkly, double &vhklz );

    void coefficients( double l, double r, double rm, double sigmal, double sigma,
                      double epsi, double alfa, double dbeta, double polTheta, double polPhi,
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
