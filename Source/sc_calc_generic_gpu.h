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
#endif

//#ifdef __OPENCL__
//#include <CL/opencl.hpp>
//#endif


#ifdef __CUDACC__
#define DBGNAN(x)
#define DBG(x)
#else
#define DBGNAN(x) //x
#define DBG(x)    //x
#endif

#ifndef __CUDACC__
#include <QDebug>
#define DD(x)  //x   // Allgemeiner Output
#define D8(x)  //x   // D8 im Main
#define D8L(x) //x   // D8 in der Library
#define DTL(x) //x   // Timing
#define DSM(x) //x   // Ausgaben bei den kleinen Zahlen
#define D1P(x) //x  // carr1p
#define D4P(x) //x  // carr4p
#include <time.h>
#endif

#ifdef __CUDACC__
#define DD(x)
#define D8(x)
#define D8L(x)
#define DTL(x)
#define DSM(x)
#define D1P(x)
#define D4P(x)
#endif


#define CHECKENDTHREAD_VAL if(_endThread)return 0;
#define CHECKENDTHREAD_RET if(_endThread)return;


#ifdef __CUDACC__
// GPU Kernel definitions (these functions will run on the GPU)
class SasCalc_GENERIC_calculation;

__global__ void kernel_GENERIC(SasCalc_GENERIC_calculation,bool);
__global__ void kernel_partSphere_peakGauss_lattFCC_ordisGauss(SasCalc_GENERIC_calculation,bool);
__global__ void kernel_partSphere_lattNone_ordisIsotropic(SasCalc_GENERIC_calculation,bool);
__global__ void kernel_partCylinder_lattNone_ordisGauss(SasCalc_GENERIC_calculation,bool);
__global__ void kernel_partCylinder_lattNone_ordisZDir(SasCalc_GENERIC_calculation,bool);
__global__ void kernel_partCylinder_peakGauss_lattBCT_ordisZDir(SasCalc_GENERIC_calculation,bool);
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

    // PARAMETER: Lattice group
    void setLType( int v ) { ltype=v; }
    void setUCA( double v ) { params.uca = v; }
    void setUCB( double v ) { params.ucb = v; }
    void setUCC( double v ) { params.ucc = v; }
    void setUCalpha( double v ) { params.ucalpha_deg = v; }
    void setUCbeta( double v ) { params.ucbeta_deg = v; }
    void setUCgamma( double v ) { params.ucgamma_deg = v; }
    void setCeffF( double val ) { params.ceff = val; }      //ZN=20305, EditTwRatio
    void setCheckBoxTwinned( bool f ) { CheckBoxTwinned = f; }
    void setCeffCyl( double val ) { epsilon=val; }  //StrToFloat(EditCeffcyl.Text);, ='conc' wird aber nicht verwendet
    void setReff( double v ) { params.reff=v; }
    void setAcpl( double v ) { acpl=v; }
    void setBcpl( double v ) { bcpl=v; }
    void setIFluc( double v ) { ifluc=v; }
    void setRFluc( double v ) { rfluc=v; }
    void setLatticeParams(int ltype, double uca, double ucb, double ucc,
                          double ucalpha, double ucbeta, double ucgamma,
                          double editceff, bool twinned, double ceffcyl,
                          double reff, double acpl, double bcpl,
                          double ifluc, double rfluc)
    {
        this->ltype=ltype;
        params.uca = uca;
        params.ucb = ucb;
        params.ucc = ucc;
        params.ucalpha_deg = ucalpha;
        params.ucbeta_deg = ucbeta;
        params.ucgamma_deg = ucgamma;
        params.ceff = editceff;
        CheckBoxTwinned = twinned;
        epsilon = ceffcyl;
        params.reff = reff;
        this->acpl = acpl;
        this->bcpl = bcpl;
        this->ifluc = ifluc;
        this->rfluc = rfluc;
    }

    // PARAMETER: Particle group
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
    void setComboBoxParticle( int val ) { ComboBoxParticle = static_cast<_cbParticle>(val); }
    typedef enum { /* 0*/ordis_Gaussian,             /* 1*/ordis_Exponential,      /* 2*/ordis_Onsager,
                   /* 3*/ordis_Maier_Saupe,          /* 4*/ordis_CutOff,           /* 5*/ordis_Laguerre,
                   /* 6*/ordis_ZDir,                 /* 7*/ordis_Isotropic,        /* 8*/ordis_mirrored_Gaussian,
                   /* 9*/ordis_mirrored_Exponential, /*10*/ordis_mirrored_Onsager, /*11*/ordis_mirrored_Maier_Saupe,
                   /*12*/ordis_mirrored_CutOff,      /*13*/ordis_FiberPattern
    } _ordistype;
    void setOrdis( int o ) { params.ordis = static_cast<_ordistype>(o); }
    typedef enum { cbintHomogeneous,    // homogeneous:=true;
                   cbintCoreShell,      // coreshell:=true;
                   cbintCoreInShell,    // coreinshell:=true;
                   cbintMultiShell,     // lipo:=true;
                   cbintMyelin          // myelin:=true;
    } _cbInterior;         /*Z0311=19828 aus den Kommentaren den Typ entnommen*/
    void setComboBoxInterior( int val ) { ComboBoxInterior = static_cast<_cbInterior>(val); }
    void setRadiusF( double val ) { params.radius = val; }  /*TPV*/
    void setRadiusI( double val ) { params.radiusi = val; }
    void setSigmaF( double val ) { params.sigma = val; }
    void setDBetaF( double db ) { params.dbeta = db /* * M_PI / 180.*/ ; }
    void setLength( double v ) { params.length=v; }     /*TPV*/
    void setSigmaL( double v ) { params.sigmal=v; } // Wird in polyliposome() von reff überschrieben!
    void setAlpha( double v ) { params.alphash1=v; }
    void setRho( double val ) { params.rho = val; } // neu
    void setRotAlpha( double v ) { params.alpha=v; }
    // bfactor = StrToFloat(EditRotAlpha.Text);  //ZN=25425
    double getRelDis() { return reldis; }
    double getDist() { return dist; }
    void setParticleParams(int cbparticle, int cbordis, int cbinterior,
                           double radius, double radiusi, double sigma,
                           double dbeta, double length, double sigmal,
                           double alphash, double rho, double alpha)
    {
        ComboBoxParticle = static_cast<_cbParticle>(cbparticle);
        params.ordis = static_cast<_ordistype>(cbordis);
        ComboBoxInterior = static_cast<_cbInterior>(cbinterior);
        params.radius = radius;
        params.radiusi = radiusi;
        params.sigma = sigma;
        params.dbeta = dbeta;
        params.length = length;
        params.sigmal = sigmal;
        params.alphash1 = alphash;
        params.rho = rho;
        params.alpha = alpha;
    }

    // PARAMETER: Peak Shape group
    typedef enum { cbpeakLorentzian=1,      // shp:=1;                          Lorentzian:=true;
                   cbpeakGaussian,          // shp:=2;                          Gaussian:=true;
                   cbpeakMod1Lorentzian,    // shp:=3;                          Lorentzian1:=true;
                   cbpeakMod2Lorentzian,    // shp:=4;                          Lorentzian2:=true;
                   cbpeakPseudoVoigt,       // shp:=5; shf:=EditPeakPar.Text;   Voigt:=true;
                   cbpeakPearsonVII,        // shp:=6; shf:=EditPeakPar.Text;   Pearson:=true;
                   cbpeakGamma,             // shp:=7; shf:=EditPeakPar.Text;   BurgerG:=true;
                   cbpeakAnisotropicGaussian// shp:=8;                          -kein Flag-
    } _cbPeak;
    void setComboBoxPeak( int v ) { shp = static_cast<_cbPeak>(v+1); }
    void setPeakPar( double p ) { /*shf=p;*/ eta=p/100.; beta=p; /*bnu=p;*/ } // TODO - unklar
    // ComboBoxLattice=13,19 -> ComboBoxPeak=7 -> shp=8
    void setDisplacement( double val ) { displacement=val; dwfactor=val*val/3.0; }
    void setAzi( double v ) { aziwidth=v; }
    void setDomainsize( double val ) { params.domainsize = val; }
    void setRadioButtonDebyeScherrer( bool f ) { RadioButtonDebyeScherrer = f; }
    void setRadioButtonPara( bool f ) { RadioButtonPara = f; }
    void setAx1( Double3 par ) { params.ax1=par; }
    void setAx2( Double3 par ) { params.ax2=par; }
    void setAx3( Double3 par ) { params.ax3=par; }
    void setSigXYZ( Double3 par ) { params.sig = 4.0 / par; }
    void setPeakShapeParams(int cbpeak, double peakpar, double debyewaller,
                            double azi, double domainsize, bool debyescherrer,
                            bool para)
    {
        shp = static_cast<_cbPeak>(cbpeak+1);
        /*shf=peakpar;*/ eta=peakpar/100.; beta=peakpar; /*bnu=peakpar;*/  // TODO - unklar
        displacement=debyewaller; dwfactor=debyewaller*debyewaller/3.0;
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
        params.ax1.setX(ax1); params.ax1.setY(ay1); params.ax1.setZ(az1);
        params.ax2.setX(ax2); params.ax2.setY(ay2); params.ax2.setZ(az2);
        params.ax3.setX(ax3); params.ax3.setY(ay3); params.ax3.setZ(az3);
        params.sig.setX(4.0/sigx); params.sig.setY(4.0/sigy); params.sig.setZ(4.0/sigz);
    }
    void setPeakShapeMatrix(Double3 ax1, Double3 ax2, Double3 ax3, Double3 sig)
    {
        params.ax1 = ax1;
        params.ax2 = ax2;
        params.ax3 = ax3;
        params.sig = 4.0 / sig;
    }

    // PARAMETER: Orientation group
    void setUCpsi( double v ) { ucpsi = v; }
    void setUCn1( double v ) { ucn1 = v; }
    void setUCn2( double v ) { ucn2 = v; }
    void setUCn3( double v ) { ucn3 = v; }
    void setPolTheta( double v ) { params.polTheta = v; }
    void setPolPhi( double v ) { params.polPhi = v; }    /*TPV*/
    void setRotTheta( double v ) { rotaxTheta = v; }
    void setRotPhi( double v ) { rotaxPhi = v; }  /*TPV?*/
    void setOrientationParams(double ucpsi, double ucn1, double ucn2, double ucn3,
                              double theta, double phi, double rottheta, double rotphi)
    {
        this->ucpsi = ucpsi;
        this->ucn1 = ucn1;
        this->ucn2 = ucn2;
        this->ucn3 = ucn3;
        params.polTheta = theta;
        params.polPhi = phi;
        rotaxTheta = rottheta;
        rotaxPhi = rotphi;
    }

    // PARAMETER: Calculation group
    void setGridPoints( int val, bool d1 ) { zmax=val;  use1d=d1; }
    void setHKLmax( int val ) { params.hklmax = hmax = kmax = lmax = val; }
    void setQMax( double val ) { qmax = val; }
    void setQMin( double val ) { qmin = val; } // valid for 1D
    void setRadQ1( bool f ) { if ( f ) calcQuadrants = radQ1; }
    void setRadQ2( bool f ) { if ( f ) calcQuadrants = radQ2; }
    void setRadQ4( bool f ) { if ( f ) calcQuadrants = radQ4; }
    void setExpandImage( bool f ) { bExpandImage = f && calcQuadrants != radQ4; }
    void setCalculationParams(int gridp, int hklmax, double qmax, bool radq1,
                              bool radq2, bool radq4, bool expand, bool d1)
    {
        zmax = gridp;
        use1d = d1;
        params.hklmax = hmax = kmax = lmax = hklmax;
        this->qmax = qmax;
        if ( radq1 ) calcQuadrants = radQ1;
        else if ( radq2 ) calcQuadrants = radQ2;
        else if ( radq4 ) calcQuadrants = radQ4;
        bExpandImage = expand && calcQuadrants != radQ4;
    }

    // PARAMETER: Experiment group
    void setpixnox( double v ) { pixnox=v; }
    void setpixnoy( double v ) { pixnoy=v; }
    void setpixx( double v_mm ) { pixx_m=v_mm/1000.; } // Übergabe aus Eingabefeld in Millimeter, intern weiter in Meter
    void setpixy( double v_mm ) { pixy_m=v_mm/1000.; } // -"-
    void setdet( double val ) { det = val; }
    void setwave( double val ) { params.wavelength=val;  wave=val; /*wave kann geändert werden*/ }
    void setBeamStop( double x, double y ) { beamX0=x; beamY0=y; }
    void setUseBeamStop( bool f ) { useBeamStop=f; }
    void setExperimentParams(double pixelnox, double pixelnoy, double pixelsizex,
                             double pixelsizey, double detdist, double wavelength,
                             double beamcx, double beamcy, bool usebs)
    {
        pixnox = pixelnox;
        pixnoy = pixelnoy;
        pixx_m = pixelsizex/1000.;  // Übergabe aus Eingabefeld in Millimeter, intern weiter in Meter
        pixy_m = pixelsizey/1000.;  // -"-
        det = detdist;
        params.wavelength=wavelength;  wave=wavelength; /*wave kann geändert werden*/
        beamX0=beamcx; beamY0=beamcy;
        useBeamStop = usebs;
    }

    // PARAMETER: Controls group
    void setBFactorF( double val ) { bfactor = val; }
    void setCheckBoxWAXS( bool v ) { CheckBoxWAXS = v; }
    void setP1( double v ) { params.p1=v; }
    // bfactor = StrToFloat(EditRotAlpha.Text);  //ZN=25425
    void setControlParams(double bfactor, bool waxs, double p1)
    {
        this->bfactor = bfactor;
        CheckBoxWAXS = waxs;
        params.p1 = p1;
    }

    // PARAMETER: Pixel Manipulation group
    void setIso( double v ) { iso=v; }
    void setIZero( double v ) { izero=v; }  /*TPV*/
    void setBase( double v ) { base=v; }    /*TPV*/
    void setPixelManipulationParams(double iso, double i0, double base)
    {
        this->iso   = iso;
        this->izero = i0;
        this->base  = base;
    }

    // Run the calculation once
    void doCalculation( int numThreads, bool bIgnoreNewSwitch );
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

    void memcleanup( void *arr );
    void endThread();

    bool gpuAvailable() { return !noGPUavailable; }

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

    typedef enum
    {
        orcGeneral,     // 1
        orcXaxis,       // 2
        orcYaxis,       // 3
        orcZaxis        // 4
    } _orcaseType;
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
        int    part, cs;
        _orcaseType orcase;
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

        _ordistype ordis;

    } _localParams;

    _localParams params;
    size_t arrCRSize;

    static const int latparlen = 2500;  // war 5000, latpar Max Check: 342, 2196, 1098
#define latparIDX(a,b,l) (a*l+b)
    // Hilfsdaten zur Ermittlung der maximal genutzten Arraylänge, um latparlen zu minimieren
    inline void latparMaxCheckInit() { params.CR->latparMaxCheckCount[0]=0; params.CR->latparMaxCheckCount[1]=0; params.CR->latparMaxCheckCount[2]=0; }
    inline void latparMaxCheck( const int i, const int a ) const
    {   if ( params.CR->latparMaxCheckCount[i]<a ) params.CR->latparMaxCheckCount[i]=a;
#ifndef __CUDACC__
        if ( a >= latparlen ) qDebug() << "################################################## latpar" << i << a << latparlen;
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
    inline size_t IDX( int x, int y ) const
    {
        return (-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y));
    }
    int minX() { return _xmin; }
    int maxX() { return _xmax; }
    int minY() { return _ymin; }
    int maxY() { return _ymax; }
    //double xyIntensity( int x, int y ) { return arrXYIntensity[3 + (-_xmin + (x)) + (_xmax-_xmin)*(-_ymin + (y))]; }
    double xyIntensity( int x, int y ) { return arrXYIntensity[3 + IDX(x,y)]; }
    double *data() { return arrXYIntensity+3; }

    int lastX() { return arrXYIntensity==nullptr ? -10000 : static_cast<int>(arrXYIntensity[0]); }
    int lastY() { return arrXYIntensity==nullptr ? -10000 : static_cast<int>(arrXYIntensity[1]); }
    bool newSwitchUsed() const { return arrXYIntensity==nullptr ? false : (static_cast<int>(arrXYIntensity[2]) != 0); }

    void scaleIntensity( bool linlog );

#ifdef FITDATA_IN_GPU  // real func decl
    bool setFitData( int sx, int sy, const double *data );
#endif

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
    bool bIgnoreNewSwitch;
    bool CheckBoxTwinned;
    bool CheckBoxWAXS;
    int    zmax;
    double dwfactor, qmax, displacement, bfactor;
    double qmin; // for 1D
    bool   use1d;

    //??? bool twin;  // wird bei LType=4 oder 5 gesetzt und bei corotations verwendet

    double beamX0, beamY0;  // Beamstop
    bool   useBeamStop;

    _cbInterior ComboBoxInterior;
    _cbParticle ComboBoxParticle;
    bool RadioButtonDebyeScherrer, RadioButtonPara;

    int  ltype; // Zur Auswahl der Operation in buttonHKLclick und auch später für die Berechnungsmethoden

    bool prepareCalculation();

    double c0, c1, c2, c3, c4; //, p1, approx, ifl, rfl, bf;
    //TODO: diese Werte müssen noch übergeben/vorbesetzt werden
    double azidom, eta, dom1, dom2, dom3;
    double beta, aziwidth, critdist, /*shf,*/ zz;
    double rotaxTheta, rotaxPhi, /*bnu,*/ /*xrdalf,*/ wave, iso, izero, base, ifluc, rfluc;
    _cbPeak shp;
    int cdim, partdim;
    double order;

    double F121,F122,F123,pqr0, sinphic, cosphic;
    double ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10;
    double ccc1,ccc2,ccc3,vv,vv3;
    double cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10;

    double acpl, bcpl;

    // internal flags
    bool /*generic=true,*/ twin, debyescherrer, paracrystal, quad1, quad2, quad4,
         lattice, lat1d, lat2d, lat3d, Lorentzian, Gaussian, Lorentzian1,
         Lorentzian2, Voigt, Pearson, BurgerG,
         /*partsphere, partcylinder, partdisk, partvesicle, partcube,
         partellipsoid, parttriellipsoid, partbarrel, partball, partchain,
         partkpchain,*/ homogeneous, coreshell, coreinshell, lipo, myelin;

    double epsilon;        //StrToFloat(EditCeffcyl.Text);

    double fhkl, qhkl; //, qhkl0;
    double qxhkl, qyhkl, qzhkl, qxyhkl;
    double qxhklt, qyhklt, qzhklt, qhklt;
    int    mhkl, h, k, l, hmax, kmax, lmax;
    double peaknorm1, peaknorm1t, peaknorm2, /*g3, g3t,*/ x2phihkl, x2phihklt; //, x2psihkl;

    //{NV} - unit cell definiton (TODO: Double3?)
    double cosa, cosb, cosg, ucpsi, ucvol,
           ucn1, ucn2, ucn3, ucn, ucn1n, ucn2n, ucn3n,
           ucphi, uctheta, cosphi, sinphi, costheta, sintheta,
           ucl1,ucl2,ucl3, pixnox,pixx_m, pixnoy,pixy_m, det;
    double ri11,ri12,ri13,ri21,ri22,ri23,ri31,ri32,ri33,
           rt11,rt12,rt13,rt21,rt22,rt23,rt31,rt32,rt33,
           nuvwx,nuvwy,nuvwz,uuvwx,uuvwy,uuvwz,vuvwx,vuvwy,vuvwz,
           nhklx,nhkly,nhklz,uhklx,uhkly,uhklz,vhklx,vhkly,vhklz;
    int   *latpar1ptr/* 6*/, *latpar2ptr/* 6*/;
    float *latpar3ptr/*17*/; //, *latpar4ptr/* 2*/; // : Array[1..5000,0..15] of real;
    //  latpar4 wird (noch) nicht wirklich gebraucht
    double latpar[4];
    int peakmax1, peakmax2;
    size_t latpar1Size, latpar2Size, latpar3Size;


    void initMemory();

    int dbgHelper[5]; // für latpar1

    void createMemory( void **ptr, size_t lensoll, size_t &lenist, bool gpuonly, const char *dbgInfo );

    void checkArrays( int minx, int maxx, int miny, int maxy );

#ifdef __CUDACC__
    __host__ __device__
#endif
    void setXYIntensity( int x/*ihex*/, int y/*i*/, double val ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline void setNewSwitchUsed( bool val ) const
    {
        if ( arrXYIntensity == nullptr ) return;
        arrXYIntensity[2] = val ? 1 : 0;
    }

    bool   noGPUavailable;

    // Intensity-Array
    // xyintensity: array[-imax..imax,-imax..imax] of real;
    double *arrXYIntensity;     //!< array for the resulatant image
    int _xmin, _xmax,           //!< array limits for the z-values (display horizontal)
        _ymin, _ymax;           //!< array limits for the i-values (display vertical)
    size_t _arrCount;           //!< number of double values in the allocated memory
    size_t arrXYsize;

#ifdef FITDATA_IN_GPU
    double *arrFitData;
    double *arrFitFqs;
    int    _fitWidth,
           _fitHeight;
    size_t _arrFitSize;
    bool   _fitEnabled;
    size_t arrFitDSize, arrFitFSize;
#endif

//    inline int   latpar1( int a, int b ) const { latparMaxCheck(0,a); return latpar1ptr[latparIDX(a,b, 6)]; }      // [5000][6], genutzt: 0,1,2,3,4,5
//    inline int   latpar2( int a, int b ) const { latparMaxCheck(1,a); return latpar2ptr[latparIDX(a,b, 6)]; }      // [5000][6], genutzt: 0,1,2,3,4,5
//    inline float latpar3( int a, int b ) const { latparMaxCheck(2,a); return latpar3ptr[latparIDX(a,b,17)]; }      // [5000][15], genutzt: 1 bis 16
//    //#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]      // [5000][15], genutzt: nichts

//    inline void setLatpar1( int a, int b, int v ) { latparMaxCheck(0,a); latpar1ptr[latparIDX(a,b, 6)] = v; }      // [5000][6], genutzt: 0,1,2,3,4,5
//    inline void setLatpar2( int a, int b, int v ) { latparMaxCheck(1,a); latpar2ptr[latparIDX(a,b, 6)] = v; }      // [5000][6], genutzt: 0,1,2,3,4,5
//    inline void setLatpar3( int a, int b, float v ) { latparMaxCheck(2,a); latpar3ptr[latparIDX(a,b,17)] = v; }    // [5000][15], genutzt: 1 bis 16
//    //#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]      // [5000][15], genutzt: nichts

    // Internal variables outside all loops
    int zzmin, zzmax, iimin, iimax;
    double ax1n_sigx, ax2n_sigy, ax3n_sigz, cubevol, phiwidth, reldis, dist;

    int noFitX0[4], noFitX1[4], noFitY0[4], noFitY1[4];
    int fitBStopPixel, fitBorderPixel;

    void corotations(double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg,
                     double u, double v, double w, double ephi,
                     bool lat1d, bool lat2d, bool lat3d, // Neu 20220311
                     double &m11, double &m12, double &m13, double &m21, double &m22, double &m23, double &m31, double &m32, double &m33,
                     double &mtw11, double &mtw12, double &mtw13, double &mtw21, double &mtw22, double &mtw23, double &mtw31, double &mtw32, double &mtw33, double &vvol,
                     double &nuvwx, double &nuvwy, double &nuvwz, double &uuvwx, double &uuvwy, double &uuvwz, double &vuvwx, double &vuvwy, double &vuvwz,
                     double &nhklx, double &nhkly, double &nhklz, double &uhklx, double &uhkly, double &uhklz, double &vhklx, double &vhkly, double &vhklz );
    void coefficients(int dim, int nmax, double &order);
    void ButtonHKLClick(int ltype) const;
    void fhkl_c(int lat, int h, int k, int l,
                double &sphno, double &fhkl, double &qhkl, double &qhkl0 ) const;
    void extinction(int lat, int h, int k, int l, int aniso,
                    int &mhkl, double &fhkl ) const;
    double gammln( double xx ) const;


    // Thread parts without GPU
    static void *doThreadCalculation(void *arg);
    static void *doThreadFitCalculation(void *arg);
    static pthread_t *threads;
    int *thread_args;
    int numberOfThreads;

    /*static*/ bool _endThread;

    void calcCPU_selection(const SasCalc_GENERIC_calculation& CALC, bool dofit, int ihex_x, int i_y);
    void kernel_selection( int Nx, int Ny, bool dofit );

#ifdef __CUDACC__
    __host__ __device__
#endif
    static void convert_idx2q( const SasCalc_GENERIC_calculation& CALC,
                      int ihex, int i, double &qx, double &qy, double &qz );

#ifdef __CUDACC__
    __host__ __device__
#endif
    static bool checkFitRanges( const SasCalc_GENERIC_calculation& CALC, int &ihex_x, int &i_y );


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
    double gammaratio( double a, double b, double x) const;

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


#ifdef __CUDACC__
    __host__ __device__
#endif
    void qrombdeltac( double p1, double sigma, double alfa,
                    double theta, double phi, double qx/*1*/, double qy/*1*/, double qz/*1*/,
                    double qxn/*9*/, double qyn/*9*/, double qzn/*9*/, double qhkl/*9*/,
                    //double ax1n, double ax2n, double ax3n,
                    //double ax1x, double ax1y, double ax1z,
                    //double ax2x, double ax2y, double ax2z,
                    //double ax3x, double ax3y, double ax3z,
                    //double sigx, double sigy, double sigz,
                    int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                    double *carr1, double &pq ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    void trapzddeltac( double a, double b, double l, double r, double p1, double sigma, double alfa,
                     double dbeta, double theta, double phi, double qx, double qy, double qz,
                     double p11, double p12, double p13, double p21, double p22, double p23,
                     double p31, double p32, double p33, double qxn, double qyn, double qzn, double qhkl,
                     //double ax1n, double ax2n, double ax3n,
                     //double ax1x, double ax1y, double ax1z,
                     //double ax2x, double ax2y, double ax2z,
                     //double ax3x, double ax3y, double ax3z,
                     //double sigx, double sigy, double sigz,
                     int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                     double *carr1, double &pq, int n, int &trapzddeltac_cnt ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    void trapzdchid( double a, double b, double l, double r, double p1, double sigma, double alfa,
                   double delta, double theta, double phi, double qx, double qy, double qz, double p11,
                   double p12, double p13, double p21, double p22, double p23, double p31, double /*p32*/,
                   double p33, double qxn, double qyn, double qzn, double qhkl,
                   //double ax1n, double ax2n, double ax3n,
                   //double ax1x, double ax1y, double ax1z,
                   //double ax2x, double ax2y, double ax2z,
                   //double ax3x, double ax3y, double ax3z,
                   //double sigx, double sigy, double sigz,
                   int /*ordis*/, int /*dim*/, int /*i0*/, int i1, int /*i2*/, int i3, int /*i4*/,
                   double *carr1, double &pq, int n,
                   int &trapzdchid_cnt ) const;


#ifdef __CUDACC__
    __host__ __device__
#endif
    void qrombchid( double l, double r, double p1, double sigma, double alfa, double delta,
                  double theta, double phi, double qx, double qy, double qz,
                  double p11, double p12, double p13, double p21, double p22, double p23,
                  double p31, double p32, double p33, double qxn, double qyn, double qzn, double qhkl,
                  //double ax1n, double ax2n, double ax3n,
                  //double ax1x, double ax1y, double ax1z,
                  //double ax2x, double ax2y, double ax2z,
                  //double ax3x, double ax3y, double ax3z,
                  //double sigx, double sigy, double sigz,
                  int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                  double *carr1, double &pq ) const;


/*
#ifdef __CUDACC__
    __host__ __device__
#endif
    double formfq( double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const;
*/

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formfq_partSphere( double q ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formfq_partCylinder( double limql, double qxs, double qys, double q ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formfq_partDisk( /*double limql,*/ double qx, double qy, double qxs, double qys, double q/*, int ordis*/ ) const;

//#ifdef __CUDACC__
//    __host__ __device__
//#endif
//    double formfq_partCube( /*double limql,*/ double qx, double qy, /*double qxs, double qys,*/ double q, int ordis ) const;



/*
#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq( double sigmal, double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const;
*/

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partSphere( double q ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partCylinder( double qxs, double qys, double q ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partDisk( double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partCube( double qxs, double qys, double q, int ordis ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partEllips( double sigmal, double qx, double qy, double qxs, double qys, double q ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partTriaxEllips( double qxs, double qys, double q ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partSuperEllips( double qxs, double qys, double q, int ordis ) const;

#ifdef __CUDACC__
    __host__ __device__
#endif
    double formpq_partSuperball( double qxs, double qys, double q, int ordis ) const;


#ifdef procnotused
#ifdef __CUDACC__
    __host__ __device__
#endif
    double polycscube(double rho1, double rho2, double p1, double p2, double alf1, double alf2, double rn, double pf,
                   double sigma, double q) const;
#endif

#ifdef undef
#ifdef __CUDACC__
    __host__ __device__
#endif
    double spy( double q ) const;
#endif


#ifdef __CUDACC__
    friend __global__ void kernel_GENERIC(SasCalc_GENERIC_calculation,bool);
    friend __global__ void kernel_partSphere_peakGauss_lattFCC_ordisGauss(SasCalc_GENERIC_calculation,bool);
    friend __global__ void kernel_partSphere_lattNone_ordisIsotropic(SasCalc_GENERIC_calculation,bool);
    friend __global__ void kernel_partCylinder_lattNone_ordisGauss(SasCalc_GENERIC_calculation,bool);
    friend __global__ void kernel_partCube_lattNone_ordisIsotropic(SasCalc_GENERIC_calculation,bool);
    friend __global__ void kernel_partCylinder_lattNone_ordisZDir(SasCalc_GENERIC_calculation,bool);
    friend __global__ void kernel_partCylinder_peakGauss_lattBCT_ordisZDir(SasCalc_GENERIC_calculation,bool);
#endif

#ifdef __CUDACC__
    __host__ __device__
#endif
    static double calc_GENERIC(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz );

#ifdef __CUDACC__
    __host__ __device__
#endif
    static double calc_partSphere_peakGauss_lattFCC_ordisGauss(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz );

#ifdef __CUDACC__
    __host__ __device__
#endif
    static double calc_partSphere_lattNone_ordisIsotropic(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz );

#ifdef __CUDACC__
    __host__ __device__
#endif
    static double calc_partCylinder_lattNone_ordisGauss(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz );

#ifdef __CUDACC__
    __host__ __device__
#endif
    static double calc_partCube_lattNone_ordisIsotropic(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz );

#ifdef __CUDACC__
    __host__ __device__
#endif
        static double calc_partCylinder_lattNone_ordisZDir(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz );

#ifdef __CUDACC__
    __host__ __device__
#endif
        static double calc_partCylinder_peakGauss_lattBCT_ordisZDir(const SasCalc_GENERIC_calculation& CALC, double qx, double qy, double qz );

};


#endif // SASCALC_GENERIC_GPU_H
