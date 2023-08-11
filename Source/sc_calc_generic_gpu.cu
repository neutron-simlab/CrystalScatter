/** Crystal calculation.
  *
  * A class with in interface for all parameters an the result and a calculation function
  * to calculate the image.
  * This was originally written in Pascal and translated into C++ by m.wagener@fz-juelich.de
  * During the translation, the variable and type names are preserved.
  *
  * Calculation functions:
  *     09. Jul. 2020:  Only for FCC Spheres with GPU
  *     FCC Spheres (Face centered cubic)
  */

#include "sc_calc_generic_gpu.h"
#include <stdlib.h>
#include <string.h>

#ifndef __CUDACC__
#define DBGFILE(x) //x
#include <QFile>
#else
#define DBGFILE(x)
#endif

#define DBGSPEC(x) //x
// Diese Debugausgaben sind innerhalb der Schleifen und produzieren viel Output.
//  Werden aber sowohl für die CPU als auch für die GPU genutzt.
//  Vorwiegend in sc_libs_gpu.cu

#include <iostream>
#include <chrono>
#include <unistd.h>
#include <signal.h>

DBGFILE( static QFile *fdbg; )

#ifdef __CUDACC__
#include <cuda.h>
#endif

/*static*/ long SasCalc_GENERIC_calculation::dbgCount;  // Zum Test...


#define myacos(x) ((x)<-1)?0.0:(((x)>+1)?M_PI:acos((x)))


#include "sc_libs_gpu.cu"
#include "sc_memory_gpu.cu"


#define GPU_2D
// If defined, use block.x and block.y etc. (new version)
// Otherwise use only block.x and iterate in each thread over y (old version)



//#############################################################################################################
// Attention for CUDA version 12 on Linux:
//
//  The Lazy Loading feature (introduced in CUDA 11.7) is now enabled by default on Linux with the 535 driver.
//  To disable this feature on Linux, set the environment variable CUDA_MODULE_LOADING=EAGER before launch.
//  Default enablement for Windows will happen in a future CUDA driver release. To enable this feature on
//  Windows, set the environment variable CUDA_MODULE_LOADING=LAZY before launch.
//
// With this Lazy Loading activated, this program will not run (August 2023).
//#############################################################################################################



SasCalc_GENERIC_calculation *SasCalc_GENERIC_calculation::inst;              //!< class instance pointer for the threads
pthread_t *SasCalc_GENERIC_calculation::threads;


SasCalc_GENERIC_calculation::SasCalc_GENERIC_calculation()
{
    inst = this;
    params.CR = nullptr;
    latpar1ptr = nullptr;
    latpar2ptr = nullptr;
    latpar3ptr = nullptr;
    //latpar4ptr = nullptr;
    threads = nullptr;
    initMemory();
    params.width_zuf = 1; // für die domainsize / width TPV Berechnung (old value)
}


#ifdef __CUDACC__
__global__ void doIntCalc_GENERIC_kernel( SasCalc_GENERIC_calculation CALC )
{
    int ihex = blockIdx.x * blockDim.x + threadIdx.x + CALC.zzmin;
    if ( ihex >= CALC.zzmax ) return;
#ifdef GPU_2D
    int i    = blockIdx.y * blockDim.y + threadIdx.y + CALC.iimin;
    if ( i >= CALC.iimax ) return;
    SasCalc_GENERIC_calculation::doIntCalc_GENERIC_F( CALC, ihex, i );
#else
    for ( int i=threadIdx.x+CALC.iimin; i<CALC.iimax; i+=blockDim.x )
        SasCalc_GENERIC_calculation::doIntCalc_GENERIC_F( CALC, ihex, i );
#endif
}

__global__ void doIntFitCalc_GENERIC_kernel( SasCalc_GENERIC_calculation CALC )
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if ( x >= CALC._fitWidth ) return;
#ifdef GPU_2D
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if ( y >= CALC._fitHeight ) return;
    SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC_F( CALC, x, y );
#else
    for ( int y=threadIdx.x; y<CALC._fitHeight; y+=blockDim.x )
        SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC_F( CALC, x, y );
#endif
}
#endif


bool SasCalc_GENERIC_calculation::prepareCalculation()
{
    dbgCount=0;  // Zum Test...

    /* ******* main program for calculation **** */  //Z=20012

    // Clear all internal flags
    //generic = true;  //Z=20055 ==> in der GUI ist kein Toggle dafür, ist aber immer gesetzt
    twin = false;  //Z=20056
    debyescherrer = false;  //Z=20057
    paracrystal = false;  //Z=20058
    quad1 = false;  //Z=20059
    quad2 = false;  //Z=20060
    quad4 = false;  //Z=20061
    lattice = true;  //Z=20062
    lat1d = false;  //Z=20063
    lat2d = false;  //Z=20064
    lat3d = true;  //Z=20065

    Lorentzian = false;  //Z=20067
    Gaussian = false;  //Z=20068
    Lorentzian1 = false;  //Z=20069
    Lorentzian2 = false;  //Z=20070
    Voigt = false;  //Z=20071
    Pearson = false;  //Z=20072
    BurgerG = false;  //Z=20073

    partsphere = false;  //Z=20075
    partcylinder = false;  //Z=20076
    partdisk = false;  //Z=20077
    partvesicle = false;  //Z=20078
    partcube = false;  //Z=20079
    partellipsoid = false;  //Z=20080
    parttriellipsoid = false;  //Z=20081
    partbarrel = false;  //Z=20082
    partball = false;  //Z=20083
    partchain = false;  //Z=20084
    partkpchain = false;  //Z=20085

    homogeneous = false;  //Z=20087
    coreshell = false;  //Z=20088
    coreinshell = false;  //Z=20089
    lipo = false;  //Z=20090
    myelin = false;  //Z=20091

    if ( ltype == 0/*Lamellae*/ )
    {
        lat1d = true;  //Z=20100
        lat3d = false;  //Z=20101
    }
    if ( (ltype==1/*hex cyl*/) || (ltype==2/*sq cyl*/) || (ltype==3/*rec cyl*/) )
    {
        lat2d = true;  //Z=20104
        lat3d = false;  //Z=20105
    }
    if ( ltype==14/*2DHex (GISAXS)*/ )
    {   /*  GiSAXS-pattern  */  //ZN=21813
        calcQuadrants = radQ4;  //ZN=21818
    }   //ZN=21819
    if ( ltype==15/*2DSquare (GISAXS)*/ )
    {   /*  GiSAXS-pattern  */  //ZN=21820
        calcQuadrants = radQ4;  //ZN=21825
    }   //ZN=21826

    zzmin = calcQuadrants == radQ4 ? -zmax : 0;     /* im Batch-Loop */
    zzmax = zmax;
    iimin = calcQuadrants == radQ1 ? 0 : -zmax;
    iimax = zmax;

    ax1n_sigx = params.ax1.length() * params.sig.x();
    ax2n_sigy = params.ax2.length() * params.sig.y();
    ax3n_sigz = params.ax3.length() * params.sig.z();

    /* ********************** */  //Z=20108
    /* ** input parameters ** */  //Z=20109
    /* ********************** */  //Z=20110
    //? lattice0 = uca;  //Z=20112
    //? p = radius/radiusi;  //Z=20119
    //? alphash = StrToFloat(EditAlpha.Text);  //Z=20122
    //? maxit = StrToInt(EditOrincr.Text);  //Z=20131
    //? dwfactoriso = ucb;  //Z=20135
    //? caille = 8*displacement*displacement/(uca*uca);  //Z=20138
    //? theocaille = 4*(1-2*radius/a)*(1-2*radius/a)/3.0;  //Z=20139
    //? delta = 2*radius;  //Z=20143
    //?später width = 4/domainsize;  //Z=20144
    //? widthiso = 1/a;  //Z=20145
    phiwidth = 4./aziwidth;  //Z=20147
    //? sigmad = displacement/(radiusi-radius+eps);  //Z=20148
    //? bnu = StrToFloat(EditPeakPar.Text);  //Z=20149
    params.p1 = params.radius/params.radiusi;  //Z=20151
    //? if ( CheckBoxHighq.Checked==false ) highq = 0  //Z=20154
    //?         else highq = 1;  //Z=20155
    //? maxit = highq;        /*  check  */  //Z=20156
    //? ella = radius;        /*  ellipsoid  */  //Z=20157
    //? ellb = length;        /*  ellipsoid  */  //Z=20158
    //? ellc = radiusi;       /*  ellipsoid  */  //Z=20159

    // --> SetLength(carr11pm^, (130+1));  //Z=20454ff.  die Arrays hier werden anders generiert
    if ( latpar1ptr == nullptr )
    {
        createMemory( (void **)(&latpar1ptr), sizeof(int) * latparlen * 6, latpar1Size, true, "latpar1" );
        createMemory( (void **)(&latpar2ptr), sizeof(int) * latparlen * 6, latpar2Size, true, "latpar2" );
        createMemory( (void **)(&latpar3ptr), sizeof(float) * latparlen * 17, latpar3Size, true, "latpar3" );
    }
    if ( params.CR == nullptr )
    {
        //        carr1p,carr2p,carr3p,carr4p,carr5p,carr6p,carr7p,carr8p,carr9p,
        //        carr1f,carr2f,carr3f,carr4f,carr5f,carr6f,carr7f,carr8f,carr9f: ^CoeffArrayType;
        //        coeffarraytype=array[0..150] of extended;
        // besser einen Speicher als Record anlegen
        createMemory( (void **)(&params.CR), sizeof(_carrXX), arrCRSize, true, "CR" );
        latparMaxCheckInit();
    }

    /*  Parameter array for liposome and myelin  */  //Z=20161
    params.CR->myarray[ 0] = params.length;     /*  axon length  */  //Z=20162
    params.CR->myarray[ 1] = params.radius;     /*  axon radius  */  //Z=20163
    params.CR->myarray[ 2] = params.sigma;      //Z=20164
    params.CR->myarray[ 3] = params.sigmal;     //Z=20165
    params.CR->myarray[ 4] = params.radiusi;    /*  no. of layers  */  //Z=20166
    params.CR->myarray[ 5] = params.alphash1;           /*  l_ex  */  //Z=20167
    params.CR->myarray[ 6] = params.rho;        /*  l_in  */  //Z=20168
    params.CR->myarray[ 7] = acpl;              /*  l_lip_head  */  //Z=20169
    params.CR->myarray[ 8] = bcpl;              /*  l_lip_tail  */  //Z=20170
    params.CR->myarray[ 9] = params.uca;        /*  phi_axon  */  //Z=20171
    params.CR->myarray[10] = params.ucb;        /*  phi_intra  */  //Z=20172
    params.CR->myarray[11] = params.ucc;        /*  phi_extra  */  //Z=20173
    params.CR->myarray[12] = params.domainsize; /*  phi_head  */  //Z=20174
    params.CR->myarray[13] = aziwidth;          /*  phi_tail  */  //Z=20175
    params.CR->myarray[14] = 1;                 /*  inmax  */  //Z=20176
    params.CR->myarray[15] = 1;                 /*  vv  */  //Z=20177
    params.CR->myarray[16] = 1;                 /*  rmax  */  //Z=20178
    params.CR->myarray[17] = iso;               /*  test parameter  */  //Z=20179

    /*  polar angles  */  //Z=20182
    ucphi = params.polPhi;  //Z=20185
    uctheta = params.polTheta;  //Z=20186
    cosphi = cos(params.polPhi*M_PI/180.0);  //Z=20189
    sinphi = sin(params.polPhi*M_PI/180.0);  //Z=20190
    costheta = cos(params.polTheta*M_PI/180.0);  //Z=20191
    sintheta = sin(params.polTheta*M_PI/180.0);  //Z=20192

    /* *************************** */  //Z=20200
    /* ** unit cell definition *** */  //Z=20201
    /* *************************** */  //Z=20202
    cosa = cos(params.ucalpha_deg*M_PI/180.0);  //Z=20209
    cosb = cos(params.ucbeta_deg*M_PI/180.0);  //Z=20210
    cosg = cos(params.ucgamma_deg*M_PI/180.0);  //Z=20211
    ucvol = params.uca*params.ucb*params.ucc*sqrt(1.0-sqr(cosa)-sqr(cosb)-sqr(cosg)+2*cosa*cosb*cosg);  //Z=20212

    /*  crystal orientation  */  //Z=20215
    ucn = sqrt(ucn1*ucn1+ucn2*ucn2+ucn3*ucn3+eps9);  //Z=20219
    ucn1n = ucn1/ucn;  //Z=20220
    ucn2n = ucn2/ucn;  //Z=20221
    ucn3n = ucn3/ucn;  //Z=20222

    if ( lat1d )
    {/*2*/ /*  for 1D-lattices  */  //Z=20225
        params.polPhi = (180/M_PI)*atan2(ucn2n,(ucn1n+eps9));  //Z=20226
        params.polTheta = (180/M_PI)*atan2(sqrt(ucn1n*ucn1n+ucn2n*ucn2n),(ucn3n+eps9));  //Z=20227
        cosphi = cos(params.polPhi*M_PI/180.0);  //Z=20230
        sinphi = sin(params.polPhi*M_PI/180.0);  //Z=20231
        costheta = cos(params.polTheta*M_PI/180.0);  //Z=20232
        sintheta = sin(params.polTheta*M_PI/180.0);  //Z=20233
    }/*2*/  //Z=20234

    if ( lat2d )
    {/*2*/ /*  for 2D-lattices  */  //Z=20236
        params.polPhi = (180/M_PI)*atan2(ucn2n,(ucn1n+eps9));  //Z=20237
        params.polTheta = (180/M_PI)*atan2(sqrt(ucn1n*ucn1n+ucn2n*ucn2n),(ucn3n+eps9));  //Z=20238
        cosphi = cos(params.polPhi*M_PI/180.0);  //Z=20241
        sinphi = sin(params.polPhi*M_PI/180.0);  //Z=20242
        costheta = cos(params.polTheta*M_PI/180.0);  //Z=20243
        sintheta = sin(params.polTheta*M_PI/180.0);  //Z=20244
        ucpsi = 90-params.polPhi;  //Z=20245
    }/*2*/  //Z=20247

    /* ************************* */  //Z=20249
    /* ** detector definition ** */  //Z=20250
    /* ************************* */  //Z=20251
    //? xrdalf = StrToFloat(EditWAXSangle.Text);    /*  [deg]  */  //Z=20253
    //? xrdalf = xrdalf*M_PI/180.0;  //Z=20254
    //? tilt = StrToFloat(EditTilt.Text);  //Z=20263
    //? tilt = -tilt*M_PI/180.0;  //Z=20264

    /* ********************* */  //Z=20266
    /* ** beam definition ** */  //Z=20267
    /* ********************* */  //Z=20268
    //? lambda = StrToFloat(EditWavelength.Text);  //Z=20271
    //? beamx0 = StrToFloat(EditCenterX.Text);  //Z=20277
    //? beamy0 = StrToFloat(EditCenterY.Text);  //Z=20278
    //? nbeam = StrToInt(EditBeamN.Text);  //Z=20279
    //? mbeam = StrToInt(EditBeamM.Text);  //Z=20280
    //? sigbeamx = StrToFloat(EditSigN.Text);  //Z=20281
    //? sigbeamy = StrToFloat(EditSigM.Text);  //Z=20282
    //? midn = round((nbeam+1)/2.0);  //Z=20283
    //? midm = round((mbeam+1)/2.0);  //Z=20284

    /* *********************** */  //Z=20287
    /* ** GISAXS-parameters ** */  //Z=20288
    /* *********************** */  //Z=20289
    /*  TODO: hierfür gibt es noch keine Eingabefelder ...
    alphai = StrToFloat(EditAlphai.Text);  //Z=20290
    betag = StrToFloat(EditAbsorb.Text);  //Z=20291
    betag2 = StrToFloat(EditAbsorb2.Text);  //Z=20292
    alphacrit = StrToFloat(EditAcrit.Text);  //Z=20293
    alphacrit2 = StrToFloat(EditAcrit2.Text);  //Z=20294
    layerh = StrToFloat(EditLayer.Text);  //Z=20295
    qhoriz = StrToFloat(EditQhoriz.Text);  //Z=20296
    amppeak = StrToFloat(EditPeakAmp.Text);  //Z=20297
    cutoff = StrToFloat(EditCutOff.Text);  //Z=20298
    alphacrit = alphacrit*M_PI/180.0;  //Z=20299
    alphacrit2 = alphacrit2*M_PI/180.0;  //Z=20300
    alphai = alphai*M_PI/180.0;  //Z=20301
    k00 = qhoriz/sin(alphai);  //Z=20302

    gisthetamax = 5.0;          //  default maximum angle in °    //Z=20304
    gisthetamax = gisthetamax*M_PI/180.0;  //Z=20305
    gisqmax = 2.0;              //  default maximum q    //Z=20306
    gislambda = 2*M_PI*sin(gisthetamax)/gisqmax;  //Z=20307
    gisk0 = 2*M_PI/gislambda;     //  default k0    //Z=20308
    gisn0 = 1.0;  //Z=20309
    gisn1 = gisn0*cos(alphacrit);  //Z=20310
    gisn2 = gisn0*cos(alphacrit2);  //Z=20311
    Editn0.Text = FloatToStr(gisn0);  //Z=20312
    Editn1.Text = FloatToStr(gisn1);  //Z=20313
    Editn2.Text = FloatToStr(gisn2);  //Z=20314
    alphai1 = acos(cos(alphai)/gisn1);  //Z=20315

    qcrit = qhoriz+k00*sin(alphacrit);  //Z=20317
    EditQcrit.Text = FloatToStr(qcrit);  //Z=20318
    rough = StrToFloat(EditVertRough.Text);  //Z=20319
    arough = StrToFloat(EditRoughAmp.Text);  //Z=20320
    if ( CheckBoxGISAXS.Checked==true ) gisaxs = true;  //Z=20321
    */

    /* ******************************* */  //Z=20324
    /* ** settings for calculations ** */  //Z=20325
    /* ******************************* */  //Z=20326
    hmax = params.hklmax;  //Z=20330
    kmax = params.hklmax;  //Z=20331
    lmax = params.hklmax;  //Z=20332
    //? zmax = StrToInt(EditGridPoints.Text);  //Z=20333  sonst zzmax...
    /*  Color Maps  */  //Z=20334
    //? ColMap  =  ComboBoxMap.ItemIndex + 1;  //Z=20335
    /* 1: Color, 2: Color Wheel, 3: Temperature, 4: Geo, 5: DESY, 6: grey, */  //Z=20336
    /* 7: inverse grey, 8: jet */  //Z=20337

    /* ******************************** */  //Z=20343
    /*  training parameters and update  */  //Z=20344
    /* ******************************** */  //Z=20345
    // --> werden an anderer Stelle gesetzt, daher hier weglassen
    // --> if ( CheckBoxTraining.checked==true )  wird auch ausgelassen

    /* *************************** */  //Z=20517
    /* ** peak shape parameters ** */  //Z=20518
    /* *************************** */  //Z=20519

    // ax1n = sqrt(ax1x*ax1x+ax1y*ax1y+ax1z*ax1z);  //Z=20533 = ax1.length()
    // ax2n = sqrt(ax2x*ax2x+ax2y*ax2y+ax2z*ax2z);  //Z=20534 = ax2.length()
    // ax3n = sqrt(ax3x*ax3x+ax3y*ax3y+ax3z*ax3z);  //Z=20535 = ax3.length()

    //? axmp^[1][1] = ax1x/ax1n;     axmp^[1][2] = ax1y/ax1n;     axmp^[1][3] = ax1z/ax1n;  //Z=20537
    //? axmp^[2][1] = ax2x/ax2n;     axmp^[2][2] = ax2y/ax2n;     axmp^[2][3] = ax2z/ax2n;  //Z=20538
    //? axmp^[3][1] = ax3x/ax3n;     axmp^[3][2] = ax3y/ax3n;     axmp^[3][3] = ax3z/ax3n;  //Z=20539

    c0 = 1;  //Z=20545
    c1 = 4;  //Z=20546
    c2 = 4*(sqrt(2.0)-1);  //Z=20547
    c3 = 4*(exp(2*log(2.0)/3.0)-1);  //Z=20548
    if ( beta == 0 )
        c4 = 1; // da bei "case cbpeakPearsonVII /*Pearson*/:" damit multipliziert wird
    else
        c4=exp(log(2.)/beta)-1;  //Z=20549

    /*  critdist:=StrToFloat(EditCritDist.Text);  //Z=20551 */
    critdist = 0.5;      /*  for allowed/non-allowed reflection overview  */  //Z=20552

    /* ******************* */  //Z=20554
    /* ** lattice tpyes ** */  //Z=20555
    /* ******************* */  //Z=20556
    //if ( CheckBoxTest.Checked==true ) generic = true; /*  generic  */  //Z=20580 ==> lassen wir jetzt mal fest

    //if ( (ltype==4/*BCC*/) || (ltype==5/*FCC*/) ) twin = true;     /*  twinning  */  //Z=20582
    //230728 Lt. Hr. Förster könnte das Twinning hier noch Probleme bereiten, daher nicht automatisch an machen
    twin = CheckBoxTwinned;  // von außen gesteuert

    if ( (ltype==12/*None*/) || (ltype==20/*Percus-Yevick*/) || (ltype==21/*Teubner-Strey*/) ) lattice = false;  //Z=20583

    if ( ltype == 14 /*2D-Hex, GiSAXS*/ ) calcQuadrants = radQ4;  /*Z=21005*/
    if ( ltype == 15 /*2D-quare, GiSAXS*/ ) calcQuadrants = radQ4;
    // RadioButtonSlice wird hier ignoriert

    //? if ( RadioButtonDebyeScherrer.Checked==true ) debyescherrer = true;  //Z=20585
    //? if ( RadioButtonPara.Checked==true ) paracrystal = true;  //Z=20586

    //? if ( ComboBoxLattice.ItemIndex==9 ) ButtonGyrVesClick(Sender);  //Z=20588
    //? if ( ComboBoxLattice.ItemIndex==10 ) ButtonPn3mVesClick(Sender);  //Z=20589
    //? if ( ComboBoxLattice.ItemIndex==11 ) ButtonIm3mVesClick(Sender);  //Z=20590

    /* ** this generates reflection tables latpar1, latpar2, latpar1a,latpar2a ** */  //Z=20634
    // ==> ebenfalls ab Zeile 21045, daher hier weglassen
    /*  end of lattice types  */  //Z=20640

    /* ***************** */  //Z=20643
    /* ** Peak Shapes ** */  //Z=20644
    /* ***************** */  //Z=20645
    switch ( shp )
    {
    case 1: Lorentzian = true; break;  //Z=20663
    case 2: Gaussian = true; break;  //Z=20664
    case 3: Lorentzian1 = true; break;  //Z=20665
    case 4: Lorentzian2 = true; break;  //Z=20666
    case 5: Voigt = true; break;  //Z=20667
    case 6: Pearson = true; break;  //Z=20668
    case 7: BurgerG = true; break;  //Z=20669
    case 8: break; // ????
    }

    /*  Particle type  */  //Z=20672
    switch ( ComboBoxParticle )
    {
    case cbpartSphere/*0*/:      //(* sphere *)  //Z=20673
        params.part = 0;      /*  sphere  */  //Z=20674
        partsphere = true;  //Z=20675
        partdim = 3;  //Z=20676
        cdim = 3;    /*  sphere  */  //Z=21000
        break;
    case cbpartCylinder/*1*/:    //(* cylinder *)  //Z=20678
        params.part = 1;      /*  cylinder  */  //Z=20679
        partcylinder = true;  //Z=20680
        partdim = 2;  //Z=20681
        cdim = 1;    /*  cylinder  */  //Z=21001
        break;
    case cbpartDisk/*2*/:        //(* disk *)  //Z=20683
        params.part = 2;      /*  disk  */  //Z=20684
        partdisk = true;  //Z=20685
        partdim = 1;  //Z=20686
        cdim = 2;    /*  disk  */  //Z=21002
        break;
    case cbpartVesicle/*3*/:     //(* vesicle *)  //Z=20688
        params.part = 3;      /*  vesicle  */  //Z=20689
        partvesicle = true;  //Z=20690
        partdim = 3;  //Z=20691
        cdim = 3;    /*  vesicle  */  //Z=21003
        break;
    case cbpartCube/*4*/:        //(* cube *)  //Z=20693
        params.part = 4;      /*  cube  */  //Z=20694
        partcube = true;  //Z=20695
        partdim = 4;  //Z=20696
        cdim = 4;    /*  cube  */  //Z=21004
        break;
    case cbpartEllipsoide/*5*/:  //(* ellipsoid *)  //Z=20698
        params.part = 5;      /*  ellipsoid  */  //Z=20699
        partellipsoid = true;  //Z=20700
        partdim = 5;  //Z=20701
        cdim = 5;    /*  ellipsoid  */  //Z=21005
        break;
    case cbpartTriaxEllips/*6*/: //(* triaxial ellipsoid *)  //Z=20703
        params.part = 6;      /*  triaxial ellipsoid  */  //Z=20704
        parttriellipsoid = true;  //Z=20705
        partdim = 6;  //Z=20706
        cdim = 6;    /*  triaxial ellipsoid  */  //Z=21006
        break;
    case cbpartSuperEllips/*7*/: //(* super ellipsoid, barrel *)  //Z=20708
        params.part = 7;      /*  super ellipsoid, barrel  */  //Z=20709
        partbarrel = true;  //Z=20710
        partdim = 7;  //Z=20711
        cdim = 7;    /*  super ellipsoid, barrel  */  //Z=21007
        break;
    case cbpartSuperball/*8*/: //  //Z=20713
        params.part = 8;      /*  superball  */  //Z=20714
        partball = true;  //Z=20715
        partdim = 8;  //Z=20716
        cdim = 8;    /*  superball  */  //Z=21008
        break;
    case cbpartChain/*9*/:  //Z=20718
        params.part = 9;      /*  excluded volume chain  */  //Z=20719
        partchain = true;  //Z=20720
        partdim = 9;  //Z=20721
        cdim = 9;    /*  excluded volume chain  */  //Z=21009
        break;
    case cbpartkpchain/*10*/:  //Z=20724
        params.part = 10;      /*  Kratky Porod chain  */  //Z=20725
        partkpchain = true;  //Z=20726
        partdim = 10;  //Z=20727
        cdim = 10;  /*  Kratky Porod chain  */  //Z=21010
        break;
    }

    /*  Particle Interior  */  //Z=20730
    switch ( ComboBoxInterior )
    {
    case cbintHomogeneous/*0*/:
        params.cs = 0;        /*  homogeneous  */  //Z=20731
        homogeneous = true;  //Z=20738
        break;
    case cbintCoreShell/*1*/:
        params.cs = 1;        /*  core/homogeneous shell  */  //Z=20732
        params.alphash1 = 0.0001;  //Z=20736
        coreshell = true;  //Z=20739
        break;
    case cbintCoreInShell/*2*/:
        params.cs = 2;        /*  core/inhomogeneous shell  */  //Z=20733
        coreinshell = true;  //Z=20740
        break;
    case cbintMultiShell/*3*/:
        params.cs = 3;        /*  multi-shell  */  //Z=20734
        lipo = true;  //Z=20741
        break;
    case cbintMyelin/*4*/:
        params.cs = 4;        /*  myelin  */  //Z=20735
        myelin = true;  //Z=20742
        break;
    }

    // Setzen der Variablen 'dist' und 'reldis' abhängig vom LatticeType
    switch ( ltype )        // ZZ= Zeilennummern aus dem File C:\SimLab\sas-crystal\Scatter Jun 2023\crystal3d1.pas
    {
    case  0:    // Lamellae             ZZ=21879
    case  1:    // hex cyl              ZZ=22779
    case  2:    // sq cyl               ZZ=23077
    case  3:    // rec cyl              ZZ=23369
    case  6:    // HCP                  ZZ=26788
    case  7:    // SC                   ZZ=26466
    case  8:    // BCT                  ZZ=27224
    case 10:    // Pn3m                 ZZ=28997
    case 12:    // None                 ZZ=?
    case 14:    // 2D-Hex, GISAXS       ZZ=29808
    case 15:    // 2D-square, GISAXS    ZZ=30038
    case 16:    // 1D-lam, GISAXS       ZZ=30235
    case 18:    // orthorombic spheres  ZZ=27550
    case 20:    // Percus-Yevick        ZZ=?
    case 21:    // Teubner-Strey        ZZ=?
    case 22:    // Pm3n, A15            ZZ=?
    case 23:    // P42/mnm (sigma)      ZZ=? (neu)
    case 24:    // Fddd-network         ZZ=? (neu)
        dist = params.uca;
        break;
    case  5:    // FCC                  ZZ=24730
    case 13:    // CP-Layers            ZZ=23688
    case 17:    // Fd3m, diamond        ZZ=27889
    case 19:    // QC                   ZZ=23992, auskommentiert
        dist = sqrt(2.0) * params.uca / 2.0;
        break;
    case  4:    // BCC                  ZZ=24309
    case  9:    // Ia3d                 ZZ=28370
    case 11:    // Im3m                 ZZ=29400
        dist = sqrt(3.0) * params.uca / 2.0;
        break;
    }
    reldis = displacement / dist;
    //qDebug() << ltype << "reldis=displacement/dist" << reldis << displacement << dist;

    //? rotnorm = sqrt(rotx*rotx+roty*roty+rotz*rotz);  //Z=20765
    //? xycur = sqrt(xcur*xcur+ycur*ycur);  //Z=20771
    //? anglecur = 2*M_PI*StrToFloat(EditAnglexy.Text)/360.0;    /*  angle of cursor  */  //Z=20772

    if ( (cdim!=9) && (cdim!=10) )
    {/*4*/  //Z=21015
        coefficients( cdim, 120, ordis, order );   //Z=21016
    }/*4*/  //Z=21020

    zz = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=21028

    //? nnnorm = norm;  //Z=21031

    // }  //Z=21036, ende von if(generic) aus Zeile 20999

    if ( /*generic=true &&*/ lattice )
    {/*3*/  //Z=21042

        /* ** this generates reflection tables latpar1, latpar2, latpar1a,latpar2a ** */  //Z=21044
        ButtonHKLClick( ltype );  //Z=21045,  ltype=0..11 und nicht mehr ....
        peakmax1 = latpar1(1,0);  //Z=21046
        peakmax2 = latpar2(1,0);  //Z=21047

        if ( peakmax1 >= latparlen )
        {
#ifndef __CUDACC__
            qDebug() << "PrepareCalculation ERROR peakmax1" << peakmax1 << " (max" << latparlen << ")";
#endif
            peakmax1 = 0;
        }
        if ( peakmax2 >= latparlen )
        {
#ifndef __CUDACC__
            qDebug() << "PrepareCalculation ERROR peakmax2" << peakmax2 << " (max" << latparlen << ")";
#endif
            peakmax2 = 0;
        }

        corotations(params.uca,params.ucb,params.ucc,params.ucalpha_deg,params.ucbeta_deg,params.ucgamma_deg,
                    ucn1,ucn2,ucn3,ucpsi,lat1d,lat2d,lat3d,ri11,ri12,ri13,ri21,ri22,ri23,ri31,ri32,ri33,  //Z=21049
                    rt11,rt12,rt13,rt21,rt22,rt23,rt31,rt32,rt33,ucvol,  //Z=21050
                    nuvwx,nuvwy,nuvwz,uuvwx,uuvwy,uuvwz,vuvwx,vuvwy,vuvwz,  //Z=21051
                    nhklx,nhkly,nhklz,uhklx,uhkly,uhklz,vhklx,vhkly,vhklz);  //Z=21052

        //? rimp^[1][1] = ri11;    rimp^[1][2] = ri12;    rimp^[1][3] = ri13;  //Z=21106
        //? rimp^[2][1] = ri21;    rimp^[2][2] = ri22;    rimp^[2][3] = ri23;  //Z=21107
        //? rimp^[3][1] = ri31;    rimp^[1][2] = ri32;    rimp^[1][3] = ri33;  //Z=21108
        //? rtmp^[1][1] = rt11;    rtmp^[1][2] = rt12;    rtmp^[1][3] = rt13;  //Z=21109
        //? rtmp^[2][1] = rt21;    rtmp^[2][2] = rt22;    rtmp^[2][3] = rt23;  //Z=21110
        //? rtmp^[3][1] = rt31;    rtmp^[1][2] = rt32;    rtmp^[1][3] = rt33;  //Z=21111

    }/*3*/  //Z=21113

    /* ** reciprocal space vector list ** */  //Z=21115
    if ( /*generic=true &&*/ lattice )
    {/*3*/  //Z=21116
        /* ************************************************* */  //Z=21117
        /* ** isotropic peak list for allowed reflections ** */  //Z=21118
        /* ************************************************* */  //Z=21119

        double sphno, qhkl0;
        for ( int peakct1=1; peakct1<=peakmax1; peakct1++ )
        {/*4*/  //Z=21121
            h = latpar1(peakct1,1);  //Z=21122
            k = latpar1(peakct1,2);  //Z=21123
            l = latpar1(peakct1,3);  //Z=21124
            mhkl = latpar1(peakct1,4);    /*  check  */  //Z=21125
            fhkl_c(ltype,h,k,l,sphno,fhkl,qhkl,qhkl0);  //Z=21126

            setLatpar1(peakct1,5, round(fhkl) ); //Z=21128, wird danach aber nie wieder gebraucht ...
            latpar[1] = ucvol;  //Z=21129
            latpar[2] = sphno;  //Z=21130
            latpar[3] = dwfactor;  //Z=21131
        }/*4*/  //Z=21140
        if ( peakmax1 == 0 )
        {   // Damit auf jeden Fall gültige Werte in das latpar[] Array geschrieben werden...
            // Sonst wird jeder Pixel zu nan weil durch cubevol geteilt wird ...
            latpar[1]=ucvol;
            latpar[2]=1;
            latpar[3]=dwfactor;
        }

        /*  anisotropic peak list  */  //Z=21143

        double g3, g3t;
        for ( int peakct2=1; peakct2<=peakmax2; peakct2++ )
        {/*4*/  //Z=21150
            h = latpar2(peakct2,1);  //Z=21151
            k = latpar2(peakct2,2);  //Z=21152
            l = latpar2(peakct2,3);  //Z=21153
            mhkl = 1;  //Z=21154
            fhkl_c(ltype,h,k,l,sphno,fhkl,qhkl,qhkl0);  //Z=21155
            setLatpar2(peakct2, 5, round(fhkl) );  //Z=21156
            setLatpar2(peakct2, 4, mhkl);  //Z=21161

            /* latpar3[peakct2,1]:=2*pi/uca;  //Z=21164 */
            setLatpar3(peakct2,1, 2*M_PI/params.uca);  //Z=21165
            qxhkl = (2*M_PI)*(ri11*h+ri21*k+ri31*l);  //Z=21166
            qyhkl = (2*M_PI)*(ri12*h+ri22*k+ri32*l);  //Z=21167
            qzhkl = (2*M_PI)*(ri13*h+ri23*k+ri33*l);  //Z=21168
            qxyhkl = (2*M_PI)*sqrt(h*h*(ri11*ri11+ri12*ri12+ri13*ri13)+k*k*(ri21*ri21+ri22*ri22+ri23*ri23));  //Z=21169
            qhkl = sqrt(qxhkl*qxhkl+qyhkl*qyhkl+qzhkl*qzhkl);  //Z=21170

            /* latpar3[peakct2,2]:=qxhkl;  //Z=21172 */
            setLatpar3(peakct2,2, qxhkl);  //Z=21173
            /* latpar3[peakct2,3]:=qyhkl;  //Z=21174 */
            setLatpar3(peakct2,3, qyhkl);  //Z=21175
            /* latpar3[peakct2,4]:=qzhkl;  //Z=21176 */
            setLatpar3(peakct2,4, qzhkl);  //Z=21177
            /* latpar3[peakct2,5]:=qhkl;  //Z=21178 */
            setLatpar3(peakct2,5, qhkl);  //Z=21179
            /* latpar3[peakct2,6]:=qhkl/(2*pi);  //Z=21180 */
            setLatpar3(peakct2,6, qhkl/(2.0*M_PI));  //Z=21181
            qxhklt = (2*M_PI)*(rt11*h+rt21*k+rt31*l);  //Z=21182
            qyhklt = (2*M_PI)*(rt12*h+rt22*k+rt32*l);  //Z=21183
            qzhklt = (2*M_PI)*(rt13*h+rt23*k+rt33*l);  //Z=21184
            qhklt = sqrt(qxhklt*qxhklt+qyhklt*qyhklt+qzhklt*qzhklt);  //ZN=21380

            setLatpar3(peakct2,7, qxhklt);  //Z=21187
            setLatpar3(peakct2,8, qyhklt);  //Z=21189
            setLatpar3(peakct2,9, qzhklt);  //Z=21191

            g3 = 4*M_PI*M_PI/(2.0*M_PI*qhkl*qhkl);  //Z=21192
            if ( twin ) g3t = 4*M_PI*M_PI/(2.0*M_PI*qhklt*qhklt);  //ZN=21390
            x2phihkl = 4*qhkl*qhkl/(M_PI*phiwidth*phiwidth);  //Z=21193
            x2phihklt = 4*qhklt*qhklt/(M_PI*phiwidth*phiwidth);  //ZN=21392
            peaknorm2=0; // zur Sicherheit

            switch ( shp )      /*Z=20614*/
            {
            case cbpeakLorentzian: // 1
                peaknorm1 = lorentznorm3(x2phihkl);  //Z=21194
                break;
            case cbpeakGaussian: // 2
                peaknorm1 = gaussnorm3(x2phihkl);  //Z=21195
                peaknorm1t = gaussnorm3(x2phihklt);  //ZN=21396
                break;
            case cbpeakMod1Lorentzian: // 3
                peaknorm1 = pearsonnorm3(x2phihkl,2);  //Z=21196
                break;
            case cbpeakMod2Lorentzian: // 4
                peaknorm1 = pearsonnorm3(x2phihkl,3/2.0);  //Z=21197
                break;
            case cbpeakPseudoVoigt: // 5
                peaknorm1 = lorentznorm3(x2phihkl);  //Z=21199
                //peaknorm2 = gaussnorm3(x2psihkl);  //ZN=21402  x2psihkl wird erst in der Pixelberechnung bestimmt
                break;
            case cbpeakPearsonVII: // 6
                peaknorm1 = pearsonnorm3(x2phihkl,beta);  //Z=21202
                break;
            case cbpeakGamma: // 7
                peaknorm1 = gaussnorm3(x2phihkl);  //Z=21203
                break;
            case cbpeakAnisotropicGaussian: // 8
                peaknorm1 = 0;  // wird lt. Prof. Förster hier nicht benötigt
                peaknorm2 = 0;
                D8( if ( peakct2 == 1 ) qDebug() << "prepCalc cbpeakAnisotropicGaussian" );
                break;
            }
#ifndef __CUDACC__
            if ( peaknorm1 == 0 ) qDebug() << "prepareCalc:" << peakct2 << shp << x2phihkl;
#endif
            setLatpar3(peakct2,10, g3);  //Z=21205
            setLatpar3(peakct2,11, peaknorm1);  //Z=21207
            setLatpar3(peakct2,12, peaknorm2);  //Z=21209
            setLatpar3(peakct2,13, qxyhkl);  //Z=21211
            setLatpar3(peakct2,14, qhklt);  //ZN=21415
            setLatpar3(peakct2,15, g3t);  //ZN=21416
            setLatpar3(peakct2,16, peaknorm1t);  //ZN=21417

        }/*4*/  //Z=21294
    }/*3*/  /*  of test or ltype=30  */  //Z=21295

    /*Z=24607 Der folgende Code ist direkt vor den ihex/i Schleifen */
    sinphic=sin(-params.polPhi*M_PI/180.);  //ZN=25418
    cosphic=cos(-params.polPhi*M_PI/180.);

    ucl1 = sintheta*cosphi;  //ZN=25421
    ucl2 = sintheta*sinphi;
    ucl3 = -costheta;

    //bfactor:=StrToFloat(EditRotAlpha.Text);  //ZN=25425

    if ( shp==8 && lattice )
    {
        D8( qDebug() << "prepCalc start qrombdeltac, theta phi:" << params.polTheta << params.polPhi );
        qrombdeltac( params.p1, params.sigma, params.alpha, params.polTheta, params.polPhi, 1,1,1,
                     9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,
                     2, 0, 0, 0, 0, params.CR->carr1p, params.norm );  //Z=25221
        D8( qDebug() << "prepCalc qrombdeltac fertig" );
    }

    // global init to avoid same calculations every time
    params.psphere_r = /*(1. + 6.*sqr(params.sigma)) * */ params.radius;
    params.psphere_z = (1. - sqr(params.sigma)) / sqr(params.sigma);

    if ( bExpandImage )
    {   // Es wird immer Speicher für das ganze Image angelegt
        checkArrays( -zmax, zmax, -zmax, zmax );
    }
    else
    {   // Es wird ein kleinerer Speicher genutzt
        checkArrays( zzmin, zzmax, iimin, iimax );
    }
#ifndef __CUDACC__
    setDebugIntensity( true );   // Damit nur beim ersten Durchlauf eine Kontrollausgabe kommt.
#endif

    return true;
}


void SasCalc_GENERIC_calculation::cleanup()
{
    std::cerr << "latpar Max Check: " << params.CR->latparMaxCheckCount[0] << ", "
              << params.CR->latparMaxCheckCount[1] << ", " << params.CR->latparMaxCheckCount[2] << std::endl;

    std::cerr << "GENERIC::cleanup:";
#ifdef COPY_FITDATA_TO_GPU  // FITDATA_IN_GPU ok, cleanup
    if ( arrDataForFitUsed || arrDataForFit != nullptr )
    {
        memcleanup( arrDataForFit );   arrDataForFit = nullptr;   std::cerr << " arrDataForFit";
    }
#endif
#ifdef FITDATA_IN_GPU  // cleanup
    if ( arrFitFqs != nullptr )
    {
        memcleanup( arrFitFqs );   arrFitFqs = nullptr;   std::cerr << " arrFitFqs";
    }
    if ( arrFitData != nullptr )
    {
        memcleanup( arrFitData );   arrFitData = nullptr;   std::cerr << " arrFitData";
    }
#endif
    memcleanup( arrXYIntensity );   arrXYIntensity = nullptr;   std::cerr << " xyIntensity";
    memcleanup( params.CR );        params.CR = nullptr;        std::cerr << " CR";
    memcleanup( latpar3ptr );       latpar3ptr = nullptr;       std::cerr << " latpar3";
    memcleanup( latpar2ptr );       latpar2ptr = nullptr;       std::cerr << " latpar2";
    memcleanup( latpar1ptr );       latpar1ptr = nullptr;       std::cerr << " latpar1";
    std::cerr << std::endl;

}


/**
 * @brief SasCalc_GENERIC_calculation::tpvPerformRandom
 * @param ids - Liste der zufällig zu ändernden Variablen mit den Faktoren
 *              in der Art: "<name>=<faktor>" pro Eintrag.
 * Mit dieser Routine werden die angegebenen Variablen per Zufall geändert. Das kann
 * leider nicht weiter oben gemacht werden, da einige Variable aus den übergebenen
 * Parametern berechnet und dann erst geändert werden.
 */
std::string SasCalc_GENERIC_calculation::tpvPerformRandom(std::list<std::string> ids)
{
    //std::cerr << "     ----- tpvPerformRandom ----- " << ids.size() << std::endl;
    //std::cerr << "tpvPerformRandom() " << ids. << std::endl;

    if ( ids.empty() )
    {   // Leeres Argument, die internen Daten löschen
        for (std::map<std::string, tpvRandomHelper*>::const_iterator it = tpvRandomOldValues.begin(); it != tpvRandomOldValues.end(); it++)
        //      std::cout << it->first << " = " << it->second << "; ";
        {
            delete it->second;
        }
        //foreach ( tpvRandomHelper *hlp, tpvRandomOldValues )
        //{
        //    delete hlp;
        //}
        tpvRandomOldValues.clear();
        params.width_zuf = 1.0;
        //std::cerr << "     ----- tpvPerformRandom ----- CLEAR" << std::endl;
        return "TPV: clear";
    }

//    typedef struct
//    {
//        double *ptr;
//        double oldval;
//        enum { norm, uca, phi, azi, dw } flag;
//    } tpvRandomHelper;
//    QHash< QString/*name*/, tpvRandomHelper* > tpvRandomOldValues;
    if ( tpvRandomOldValues.empty() ) // tpvPerformRandom()
    {   // Erster Aufruf, füllen der Daten
        //std::cerr << "     ----- tpvPerformRandom ----- INIT" << std::endl;
        tpvRandomHelper *hlp;
#define TPVINIT(nam,var,flg) hlp = new tpvRandomHelper; \
                             hlp->ptr = &var; \
                             hlp->oldval = var; \
                             hlp->flag = flg; \
                             tpvRandomOldValues.insert({nam,hlp});
        TPVINIT( "I0",              izero,             tpvRandomHelper::norm );
        TPVINIT( "EditRadius",      params.radius,     tpvRandomHelper::norm );
        TPVINIT( "EditSigma",       params.sigma,      tpvRandomHelper::norm );
        TPVINIT( "Length",          params.length,     tpvRandomHelper::norm );
        TPVINIT( "SigmaL",          params.sigmal,     tpvRandomHelper::norm );
        TPVINIT( "uca",             params.uca,        tpvRandomHelper::norm );
        TPVINIT( "ucb",             params.ucb,        tpvRandomHelper::uca  );
        TPVINIT( "ucc",             params.ucc,        tpvRandomHelper::uca  );
        TPVINIT( "ucpsi",           ucpsi,             tpvRandomHelper::norm );
        TPVINIT( "EditDbeta",       params.dbeta,      tpvRandomHelper::norm );
        TPVINIT( "EditDomainSize",  params.width_zuf,  tpvRandomHelper::norm ); // width ist nur lokal
        TPVINIT( "EditAzi",         phiwidth,          tpvRandomHelper::azi  ); hlp->oldval = aziwidth;
        TPVINIT( "EditDebyeWaller", dwfactor,          tpvRandomHelper::dw   ); hlp->oldval = displacement;
        TPVINIT( "phi",             params.polPhi,               tpvRandomHelper::phi  );
        // TODO: noch nicht in der GUI
        TPVINIT( "Base",            base,              tpvRandomHelper::norm );
        TPVINIT( "EditRho",         params.rho,        tpvRandomHelper::norm );
        TPVINIT( "EditDbeta",       params.dbeta,      tpvRandomHelper::norm );
    }

#define zuf() (((static_cast<double>(rand())/RAND_MAX * 200.0) - 100.0) / 100.0) // -1 .. 1

    std::string retval="";
    for ( const std::string &id : ids )
    {
        int pos = id.find("=");
        std::string key = id.substr(0,pos);
        std::string val = id.substr(pos+1);
        double dval = std::stod(val);
        if ( dval < 0.00001 && val.length() > 2 )
        {   // Typische Werte sind hier 0.1 bis 0.9 und keine negativen Werte!
            std::string::size_type pk = val.find(".");
            if ( pk != std::string::npos )
                val[pk] = ',';
            else
            {
                pk = val.find(",");
                val[pk] = '.';
            }
            dval = std::stod(val);
        }
        //std::cerr << "     ----- tpvPerformRandom ----- Key=" << key << ", val=" << dval << std::endl;
        if ( tpvRandomOldValues.count(key) == 0 ) continue; // Dann fehlt oben etwas ...
        tpvRandomHelper *hlp = tpvRandomOldValues[key];
        switch ( hlp->flag )
        {
        case tpvRandomHelper::norm:
            *hlp->ptr = hlp->oldval * ( 1 + dval * zuf() );
            break;
        case tpvRandomHelper::uca:
            *hlp->ptr = params.uca;
            break;
        case tpvRandomHelper::phi:
            *hlp->ptr = 2.0 + (zuf() + 1) * 43.0;   // [2, 88]
            break;
        case tpvRandomHelper::azi:
            aziwidth = hlp->oldval;
            *hlp->ptr = 4.0 / aziwidth;
            *hlp->ptr = *hlp->ptr * ( 1 + dval * zuf() );
            break;
        case tpvRandomHelper::dw:
            displacement = hlp->oldval;
            *hlp->ptr = displacement*displacement/3.0;
            *hlp->ptr = *hlp->ptr * ( 1 + dval * zuf() );
            break;
        }
        retval += key + " = " + std::to_string(*(hlp->ptr)) + EOL;
        //std::cerr << "     ----- tpvPerformRandom ----- NewVal=" << *(hlp->ptr) << std::endl;
    }
    // Sicherheitsprüfungen
    if ( (ComboBoxInterior == cbintCoreShell || ComboBoxParticle == cbpartVesicle) &&
         (params.length <= params.radius) )
        params.length = params.radius + params.length;
    if ( ComboBoxParticle == cbpartDisk && params.length <= 2.*params.radius )
        params.length = 2.*params.radius + params.length;

    return retval;
}


/**
 * @brief SasCalc_gpu::doCalculation - main calculation procedure.
 * @param numThreads - number of threads to use
 * All other parameters must be stored in global variables before calling this procedure.
 */
void SasCalc_GENERIC_calculation::doCalculation( int numThreads )
{
    numberOfThreads = numThreads;

    auto start1 = std::chrono::high_resolution_clock::now();
    if ( !prepareCalculation() ) return;
    // Der Speicher wird erst am Ende von prepareCalculation angelegt
    setXYIntensity( zzmin, iimin, 0.0 );    // reset internal indices
/*
    if ( tpvRandomOldValues.empty() ) // doCalculation()
    {   // Jetzt ist keine TPV-Aktion, daher hier die width bestimmen,
        // vorher war dies eine lokale Variable in dem GPU-Code.
        std::cerr << "TPV DO CALCULATION width_zuf = " << params.width_zuf << std::endl;
    }
*/
    DBGFILE(if ( numThreads > 1 ) fdbg=nullptr;
            else { fdbg = new QFile("dbgfile.txt");
            if ( !fdbg->open(QIODevice::Append) )
                fdbg->open(QIODevice::WriteOnly);
            else
                fdbg->write("------\nihex; i;   qx; qy; qz;   pixval;\n");
            qDebug() << "DBGFILE:" << fdbg->fileName();
            } )

    auto end1 = std::chrono::high_resolution_clock::now();
    auto calc_time1 = std::chrono::duration_cast<std::chrono::duration<float>>(end1-start1);
    higResTimerElapsedPrep = calc_time1.count()*1000.0;

    // Use a high resolution clock to get the calculation time of the GPUs
    auto start2 = std::chrono::high_resolution_clock::now();

#ifdef __CUDACC__
    if ( gpuAvailable() && numThreads == 0 )    // GPU abschaltbar
    {
        _endThread = false;

#ifdef GPU_2D
        // GPU Programming with CUDA @ JSC (Kurs Apr.2022)
        // GPU 01 Introduction.pdf - Seite 52
        int Nx = zzmax - zzmin;
        int Ny = iimax - iimin;
        dim3 blockDim(16, 16);
        int gx = (Nx % blockDim.x == 0) ? (Nx / blockDim.x) : (Nx / blockDim.x + 1);
        int gy = (Ny % blockDim.y == 0) ? (Ny / blockDim.y) : (Ny / blockDim.y + 1);
        dim3 gridDim(gx, gy);

        //std::cerr << "Kernel-Call: "
        //          << "ImgDim=" << Nx<<"*"<<Ny<<"="<<Nx*Ny
        //          << ", blockDim=" << blockDim.x<<"*"<<blockDim.y
        //          << " = " << blockDim.x*blockDim.y << " threads"
        //          << ", gridDim=" << gx<<"*"<<gy << std::endl;

        doIntCalc_GENERIC_kernel<<<gridDim,blockDim>>>(*this);
#else
        doIntCalc_GENERIC_kernel<<<zzmax-zzmin,256>>>(*this);
#endif
        cudaError_t err = cudaGetLastError();
        if ( err != cudaSuccess )
            std::cerr << "#####   CUDA ERROR: " << cudaGetErrorString(err) << std::endl;

        cudaDeviceSynchronize();
    }
    else
    {
#endif

        D8( numThreads=1; numberOfThreads=1; )

        _endThread = false;

        if ( numThreads <= 1 )
        {   // Without threading it is simple ...
            for ( int ihex=zzmin; ihex<zzmax; /*ihex++*/ )  // (***   z-loop  ***)
            {
                if ( _endThread ) break;
                doIntCalc_GENERIC( ihex++ );
            } // for ihex
        } // numThreads <= 1

        else
        {   // Special work for Multi-Threads...
            if ( threads != nullptr ) { delete threads; delete thread_args; }
            threads     = new pthread_t[numThreads];    // Weil die Anzahl der Threads geändert werden kann
            thread_args = new int[numThreads];
            memset( threads, 0, sizeof(pthread_t)*numThreads ); // for ( int t=0; t<numThreads; t++ ) threads[t] = 0;
            int ihex=zzmin;
            while ( ihex<zzmax )  // (***   z-loop  ***)
            {
                if ( _endThread ) break;

                // If the number of threads is greater than the number of cores, this loop must be optimized
                // so that it will not wait for all threads. This might speed up a little more.
                for ( int t=0; t<numThreads && ihex<=zzmax; t++ )
                {
                    thread_args[t] = ihex++;  // threadid > 0
                    pthread_create( &threads[t], nullptr, doThreadCalculation, &thread_args[t] );
                    if ( _endThread ) break;
                }
                for ( int t=0; t<numThreads; t++ )
                {
                    if ( _endThread ) break;
                    if ( threads[t] )
                    {
                        pthread_join(threads[t], nullptr);
                        threads[t] = 0;
                    }
                }
            } // for ihex
            for ( int t=0; t<numThreads; t++ )
            {
                if ( _endThread ) break;
                if ( threads[t] )
                {
                    pthread_join(threads[t], nullptr);
                    threads[t] = 0;
                }
            }
            //delete threads;
            //delete thread_args;
            //threads = nullptr;
        } // numThreads > 1

#ifdef __CUDACC__
    } // if ( ! noGPUavailable ) else {
#endif

    performImageExpansion();

    DBGFILE( if ( fdbg ) fdbg->close(); )

    // Use a high resolution clock to get the calculation time of the GPUs
    auto end2 = std::chrono::high_resolution_clock::now();
    auto calc_time2 = std::chrono::duration_cast<std::chrono::duration<float>>(end2-start2);
    higResTimerElapsedCalc = calc_time2.count()*1000.0;
}



/**
 * @brief SasCalc_gpu::doCalculation - main calculation procedure.
 * @param numThreads - number of threads to use
 * All other parameters must be stored in global variables before calling this procedure.
 */
double SasCalc_GENERIC_calculation::doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt )
{
    numberOfThreads = numThreads;
    fitBStopPixel  = bstop;
    fitBorderPixel = border;

    auto start1 = std::chrono::high_resolution_clock::now();
    if ( !prepareCalculation() ) return 1e20; // Fehlermeldung beim Fit
    // Der Speicher wird erst am Ende von prepareCalculation angelegt
    setXYIntensity( zzmin, iimin, 0.0 );    // reset internal indices

    auto end1 = std::chrono::high_resolution_clock::now();
    auto calc_time1 = std::chrono::duration_cast<std::chrono::duration<float>>(end1-start1);
    higResTimerElapsedPrep = calc_time1.count()*1000.0;

    DBGFILE(if ( numThreads > 1 ) fdbg=nullptr;
            else { fdbg = new QFile("dbgfitfile.txt");
            if ( !fdbg->open(QIODevice::Append) )
                fdbg->open(QIODevice::WriteOnly);
            else
                fdbg->write("------\nx; y;   xdet; ydet;   qx; qy; qz;   pixval; fitdata; FQS\n");
            qDebug() << "DBGFILE:" << fdbg->fileName();
            } )

    // Use a high resolution clock to get the calculation time of the GPUs
    auto start2 = std::chrono::high_resolution_clock::now();

#ifdef __CUDACC__
    if ( gpuAvailable() && numThreads == 0 )    // GPU abschaltbar
    {
        _endThread = false;

#ifdef GPU_2D
        // GPU Programming with CUDA @ JSC (Kurs Apr.2022)
        // GPU 01 Introduction.pdf - Seite 52
        int Nx = _fitWidth;
        int Ny = _fitHeight;
        dim3 blockDim(16, 16);
        int gx = (Nx % blockDim.x == 0) ? (Nx / blockDim.x) : (Nx / blockDim.x + 1);
        int gy = (Ny % blockDim.y == 0) ? (Ny / blockDim.y) : (Ny / blockDim.y + 1);
        dim3 gridDim(gx, gy);

        doIntFitCalc_GENERIC_kernel<<<gridDim,blockDim>>>(*this);
#else
        doIntFitCalc_GENERIC_kernel<<<_fitWidth,256>>>(*this);
#endif
        cudaError_t err = cudaGetLastError();
        if ( err != cudaSuccess )
            std::cerr << cudaGetErrorString(err) << std::endl;

        cudaDeviceSynchronize();
    }
    else
    {
#endif

        D8( numThreads=1; numberOfThreads=1; )

        _endThread = false;

        if ( numThreads <= 1 )
        {   // Without threading it is simple ...
            for ( int x=0; x<_fitWidth; /*x++*/ )
            {
                if ( _endThread ) break;
                doIntFitCalc_GENERIC( x++ );
            } // for x
        } // numThreads <= 1

        else
        {   // Special work for Multi-Threads...
            if ( threads != nullptr ) { delete threads; delete thread_args; }
            threads     = new pthread_t[numThreads];    // Weil die Anzahl der Threads geändert werden kann
            thread_args = new int[numThreads];
            memset( threads, 0, sizeof(pthread_t)*numThreads ); // for ( int t=0; t<numThreads; t++ ) threads[t] = 0;
            int x=0;
            while ( x<_fitWidth /*&& x<10/ *TEST*/ )
            {
                if ( _endThread ) break;

                // If the number of threads is greater than the number of cores, this loop must be optimized
                // so that it will not wait for all threads. This might speed up a little more.
                for ( int t=0; t<numThreads && x<_fitWidth; t++ )
                {
                    thread_args[t] = x++;  // threadid > 0
                    pthread_create( &threads[t], nullptr, doThreadFitCalculation, &thread_args[t] );
                    if ( _endThread ) break;
                }
                for ( int t=0; t<numThreads; t++ )
                {
                    if ( _endThread ) break;
                    if ( threads[t] )
                    {
                        pthread_join(threads[t], nullptr);
                        threads[t] = 0;
                    }
                }
            } // for ihex
            for ( int t=0; t<numThreads; t++ )
            {
                if ( _endThread ) break;
                if ( threads[t] )
                {
                    pthread_join(threads[t], nullptr);
                    threads[t] = 0;
                }
            }
            //delete threads;
            //delete thread_args;
            threads = nullptr;
        } // numThreads > 1

#ifdef __CUDACC__
    } // if ( ! noGPUavailable ) else {
#endif

    // Ergebnis bestimmen
    double erg = 0.0;
    double *p  = arrFitFqs;
    cnt = 0;
    nancnt = 0;
    for ( size_t i=0; i<_arrFitSize; i++, p++ )
    {
        if ( isnan(*p) )  // kommt leider (noch) vor, daher hier ignorieren, damit der Fit klappt
            nancnt++;
        else
        {
            if ( *p != 0 ) cnt++;   // Echte NUll in der FQS ist (fast) nur bei noFitRect Pixeln
            erg += *p;
        }
    }

    // Use a high resolution clock to get the calculation time of the GPUs
    auto end2 = std::chrono::high_resolution_clock::now();
    auto calc_time2 = std::chrono::duration_cast<std::chrono::duration<float>>(end2-start2);
    higResTimerElapsedCalc = calc_time2.count()*1000.0;

    DBGFILE( if ( fdbg ) { fdbg->write(qPrintable(QString("Erg=%1  cnt=%2  nan=%3\n").arg(erg).arg(cnt).arg(nancnt))); fdbg->close(); } )

    if ( _endThread ) return 1e20;
    return erg;
}


/**
 * @brief SasCalc_gpu::doThreadCalculation
 * @param arg - parameter for the calculation function
 * @return allways nullptr
 * This is a static function for the threads used to calculate.
 */
void *SasCalc_GENERIC_calculation::doThreadCalculation(void *arg)
{
    pthread_setcancelstate( PTHREAD_CANCEL_ENABLE, nullptr );
    pthread_setcanceltype( PTHREAD_CANCEL_ASYNCHRONOUS, nullptr );
    int ihex = *(static_cast<int*>(arg));
    inst->doIntCalc_GENERIC( ihex );
    return nullptr;
}

/**
 * @brief SasCalc_gpu::doThreadCalculation
 * @param arg - parameter for the calculation function
 * @return allways nullptr
 * This is a static function for the threads used to calculate.
 */
void *SasCalc_GENERIC_calculation::doThreadFitCalculation(void *arg)
{
    pthread_setcancelstate( PTHREAD_CANCEL_ENABLE, nullptr );
    pthread_setcanceltype( PTHREAD_CANCEL_ASYNCHRONOUS, nullptr );
    int x = *(static_cast<int*>(arg));
    inst->doIntFitCalc_GENERIC( x );
    return nullptr;
}


/**
 * @brief SasCalc_gpu::doIntCalc_GENERIC
 * @param ihex - vertical pixel index
 */
void SasCalc_GENERIC_calculation::doIntCalc_GENERIC(int ihex)
{
    for ( int i=iimin; i<iimax; i++ )
    {
        doIntCalc_GENERIC_F( *this, ihex, i );
        if ( _endThread ) break;
    }
}

/**
 * @brief SasCalc_gpu::doIntCalc_GENERIC
 * @param ihex - vertical pixel index
 */
void SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC(int x)
{
//#ifndef __CUDACC__
//    qDebug() << "Fit-x:" << x << _fitWidth << "y-bis" << _fitHeight;
//#endif
    for ( int y=0; y<_fitHeight; y++ )
    {
        doIntFitCalc_GENERIC_F( *this, x, y );
        if ( _endThread ) break;
    }
}




/**
 * @brief SasCalc_gpu::doIntCalc_GENERIC_F
 * @param CALC  - reference to class with parameters and subfunctions
 * @param ihex - vertical pixel index
 * @param i    - horizontal pixel index
 * @param lamv - running vertical axis value
 * Calculation from Pascalprogram (20210818-crystal3d1.pas + updates) - only part Generic
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline void SasCalc_GENERIC_calculation::doIntCalc_GENERIC_F( const SasCalc_GENERIC_calculation& CALC,
                                                              int ihex, int i )
{
#ifdef COPY_FITDATA_TO_GPU  // Eckpixeltest
    if ( CALC.arrDataForFitUsed )
    {   // Spezielle Aktionen (Maskieren und FQS) für das Simplex-2D-Fit
        if ( CALC.arrDataForFit[1+CALC.IDX(ihex,i)] <= CALC.arrDataForFit[1] )
        {   // Datenwert <= Eckpixel --> Ergebnis immer 0.
            CALC.setXYIntensity( ihex, i, 0.0 );
            return;
        }
    }
#endif

    // Bestimmen von qx,qy,qz
    double qx, qy, qz;
    if ( CALC.useBeamStop )
    {   // Einrechnen des Beamstops (d.h. Verschiebung des Zentrums)
        // Hier wird die Pixelsize aber nicht der qmax Wert beachtet.
        double mdet = (ihex/*+CALC.zmax+1*/)*CALC.pixnoy/(2.0*CALC.zmax);      /* mth pixel */
        double ndet = (i   /*+CALC.zmax+1*/)*CALC.pixnox/(2.0*CALC.zmax);      /* nth pixel */
        /*Z=24635*/
        // Im Pascal-Programm ist der Beamstop für Null in der Ecke.
        // Hier ist der Beamstop für Null in der Mitte gerechnet.
        double xdet = CALC.pixx_m * (ndet - CALC.beamX0);
        double ydet = CALC.pixy_m * (mdet - CALC.beamY0);
        double rdet = sqrt(xdet*xdet+ydet*ydet);
        double phidet = atan2(ydet,xdet);
        double thetadet = atan2(rdet,CALC.det);
        qx = 2*M_PI*cos(phidet)*sin(thetadet)/CALC.wave;
        qy = 2*M_PI*sin(phidet)*sin(thetadet)/CALC.wave;
        qz = 2*M_PI*(1-cos(thetadet))/CALC.wave;
        //if ( ihex == 10 && i >= 10 && i <= 20 )
        //    qDebug() << "ihex/i" << ihex << i << "BSx/y" << CALC.beamX0 << CALC.beamY0 << "zmax" << CALC.zmax << "m/ndet" << mdet << ndet << "qxy" << qx << qy;
    }
    else
    {   // Den Mittelpunkt nutzen (ihex=0,i=0)
        // hier wird nur der qmax Wert beachtet und die Pixelsize außer Acht gelassen.
        qx = CALC.qmax * i    / (double)(CALC.zmax); // lamu
        qy = CALC.qmax * ihex / (double)(CALC.zmax); // lamv
        qz = 1e-20;
        //if ( ihex == 10 && i >= 10 && i <= 20 )
        //    qDebug() << "ihex/i" << ihex << i << "BSx/y" << CALC.beamX0 << CALC.beamY0 << "zmax" << CALC.zmax << "qxy" << qx << qy;
    }

    // alles mit "RadioButtonSlice.Checked=true" lasse ich weg, da ich hier noch keine Schnitte behandele.

    // Call q(x,y,z) Routine (also from 2dFit)
    double pixval = doIntCalc_GENERIC_q_xyz( CALC, qx, qy, qz,
#ifndef __CUDACC__
                                            (ihex==0 && i==0), fabs(ihex-CALC.beamX0)<1 && fabs(i-CALC.beamY0)<1 );
#else
                                            false, false );
#endif
    if ( CALC._endThread ) return;

    CALC.setXYIntensity( ihex, i, pixval );

#ifndef __CUDACC__
    //if ( pixval < 0.2 || isnan(pixval) ) qDebug() << pixval << "NULL" << ihex << i << qx << qy << qz;
    CALC.setDebugIntensity( false );   // Damit nur beim ersten Durchlauf eine Kontrollausgabe kommt.
#endif

    DBGFILE( if ( fdbg ) fdbg->write(qPrintable(QString("%1; %2;   %3; %4; %5;   %6;\n")
                                       .arg(ihex).arg(i)
                                       .arg(qx).arg(qy).arg(qz)
                                       .arg(pixval)
                                   )); )


} /* doIntCalc_GENERIC_F() */


/**
 * @brief SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC_F
 * @param CALC - Reference to the class to access the variables
 * @param x    - Pixelindex hor.
 * @param y    - Pixelindex vert.
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline void SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC_F( const SasCalc_GENERIC_calculation& CALC,
                                                        int x, int y )
{
    // There are up to 4 rectangle regions to be ignored during the fit.
    for ( int i=0; i<4; i++ )
    {
        if ( CALC.noFitX0[i] < 0 ) continue;
        if ( x < CALC.noFitX0[i] ) continue;
        if ( CALC.noFitX1[i] < x ) continue;
        if ( y < CALC.noFitY0[i] ) continue;
        if ( CALC.noFitY1[i] < y ) continue;
        size_t idx = x + (CALC._fitWidth * y);
        CALC.arrFitFqs[idx] = 0.0;
        DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; FitRect %4\n")
                                           .arg(x).arg(y).arg(CALC.arrFitData[idx]).arg(i)
                                       )); )
        return;
    }

    if ( CALC.fitBorderPixel > 0 || CALC.fitBStopPixel > 0 )
    {   // Ausblendungen über Pixelangaben an Rand und Mitte
        if ( x < CALC._xmin + CALC.fitBorderPixel || x >= CALC._xmax - CALC.fitBorderPixel ||
             y < CALC._ymin + CALC.fitBorderPixel || y >= CALC._ymax - CALC.fitBorderPixel )
        {
            size_t idx = x + (CALC._fitWidth * y);
            CALC.arrFitFqs[idx] = 0.0;
            DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                    fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; BorderPixel\n")
                                               .arg(x).arg(y).arg(CALC.arrFitData[idx])
                                           )); )
            return;
        }
        if ( x >= CALC.beamX0 - CALC.fitBStopPixel &&
             x <  CALC.beamX0 + CALC.fitBStopPixel &&
             y >= CALC.beamY0 - CALC.fitBStopPixel &&
             y <  CALC.beamY0 + CALC.fitBStopPixel )
        {
            size_t idx = x + (CALC._fitWidth * y);
            CALC.arrFitFqs[idx] = 0.0;
            DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                    fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; BSpixel\n")
                                               .arg(x).arg(y).arg(CALC.arrFitData[idx])
                                           )); )
            return;
        }
    }
    else if ( CALC.fitBStopPixel == -1 )
    {   // Ausblendungen per Eck-Pixel-Wert
        size_t idx = x + (CALC._fitWidth * y);
        if ( CALC.arrFitData[idx] <= CALC.arrFitData[0] )
        {
            CALC.arrFitFqs[idx] = 0.0;
            DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                    fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; Ecke %4 <= %5\n")
                                               .arg(x).arg(y).arg(CALC.arrFitData[idx]).arg(CALC.arrFitData[idx]).arg(CALC.arrFitData[0])
                                           )); )
            return;
        }
    }
    //else
    //{   // Keine Ausblenungen
    //}

    int xx = x - CALC.beamX0;
    int yy = y - CALC.beamY0;
    // Einrechnen des Beamstops (d.h. Verschiebung des Zentrums)
    double xdet = CALC.pixx_m * xx;
    double ydet = CALC.pixy_m * yy;
    double rdet = sqrt(xdet*xdet+ydet*ydet);
    double phidet = atan2(ydet,xdet);
    double thetadet = atan2(rdet,CALC.det);
    double qx = 2*M_PI*cos(phidet)*sin(thetadet)/CALC.params.wavelength;
    double qy = 2*M_PI*sin(phidet)*sin(thetadet)/CALC.params.wavelength;
    double qz = 2*M_PI*(1-cos(thetadet))/CALC.params.wavelength;

    //DBGFILE( if ( !( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 ) ) return; )

    // Call q(x,y,z) Routine (also from 2dFit)
    double pixval = doIntCalc_GENERIC_q_xyz( CALC, qx, qy, qz,
#ifndef __CUDACC__
                                            (x==0 && y==0), xx<1 && yy<1 );
#else
                                            false, false );
#endif
    if ( CALC._endThread ) return;

    size_t idx = x + (CALC._fitWidth * y);
    if ( pixval > 0 && CALC.arrFitData[idx] > 0 )
        CALC.arrFitFqs[idx] = FQSVERGL( pixval, CALC.arrFitData[idx] );
    else
        CALC.arrFitFqs[idx] = 0.0;

    DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                fdbg->write(qPrintable(QString("%1; %2;   %3; %4;   %5; %6; %7;   %8; %9; %10\n")
                                           .arg(x).arg(y)
                                           .arg(xdet).arg(ydet)
                                           .arg(qx).arg(qy).arg(qz)
                                           .arg(pixval).arg(CALC.arrFitData[idx]).arg(CALC.arrFitFqs[idx])
                                       )); )

//    if ( (x < 4 /*&& y < 5*/) || (x >= CALC._fitWidth-5 && y >= CALC._fitHeight) )
//        qDebug() << "Fit-v:" << x << y << idx
//                 << "val" << pixval << CALC.arrFitData[idx] << CALC.arrFitFqs[idx]
//                 << "q" << qx << qy << qz;

} /* doIntFitCalc_GENERIC_F() */



/**
 * @brief SasCalc_GENERIC_calculation::doIntCalc_GENERIC_q_xyz
 * @param CALC - reference to class with parameters and subfunctions
 * @param qx   - Coordinates in the q dimension
 * @param qy   - "
 * @param qz   - "
 * @param idxCenter  - true if current pixel is in the image center (pixel indices==0)
 * @param beamCenter - true if current pixel is in the beam center
 *                     Used specially for some debug or control output if not running in CUDA environment
 * @return     - calculated value for this coordinate
 * Calculation from Pascalprogram (20210818-crystal3d1.pas + updates) - only part Generic
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::doIntCalc_GENERIC_q_xyz(const SasCalc_GENERIC_calculation& CALC,
                                                                   double qx, double qy, double qz,
                                                                   bool idxCenter, bool /*beamCenter*/)
{
    double /*lamu,*/ /*qabs,*/ pq, fq, intensity, radintensity;
    double /*shkl,*/ /*fhkl,*/ x2;
    double dqs1, dqs2, dqs3, sq;
    //Double3 q3hkl; // qxhkl, qyhkl, qzhkl
    //Double3 dq3;   // dqx, dqy, dqz

    double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254
    //s = q/(2*M_PI);

    //Z=25257 if gisaxs ... erstmal weglassen und das auch bei allen weiteren Abfragen!!!

    double delta = 2.0 * CALC.params.radius;
    double pqiso=1;
    double limql=0;  // Wird als Parameter bei formpq() verwendet und z.T. kurz vorher berechnet
    double limqlf=0; // Wird nur als Parameter bei formfq() verwendet, sonst nicht
    double szqiso, szq;

    fq = 1; // Verhindern einer Compiler-Warnung
    // <fq> wird bei der Berechnung von szq und somit für die Intensität verwendet (=fq/pq).
    // Aber nur dann, wenn lattice gesetzt ist. Und nur dann wird es auch per formfq berechnet.


    /* ************* */  //Z=25269  //ZN=25476
    /* ** spheres ** */  //Z=25270
    /* ************* */  //Z=25271
    if ( CALC.partsphere )
    {/*6*/  //Z=25272

        limql = 1;
        pq = CALC.formpq( CALC.params.sigmal, limql, qx, qy, qx, qy, q, CALC.ordis );  //Z=25276

        /* fq = pq;  //Z=25277 */
        if ( CALC.lattice )  //Z=25283
            fq=CALC.formfq( limqlf, qx, qy, qx, qy, q, CALC.ordis );  //Z=25285
        pqiso = pq;  //Z=25297

    }/*6*/  /*  of part=0  */  //Z=25299


    /* *************** */  //Z=25302
    /* ** cylinders ** */  //Z=25303
    /* *************** */  //Z=25304
    if ( CALC.partcylinder )
    {/*6*/     /*  cylinders  */  //Z=25305

        /*  isotropic cases  */  //Z=25307
        if ( CALC.ordis==7 )
        {/*7*/  //Z=25308
            pq = CALC.formpq(CALC.params.sigmal, q, qx, qy, qx, qy, q, CALC.ordis);
            /* fq = pq;  //Z=25311 */
            if ( CALC.lattice )  //Z=25312
                fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );
        }/*7*/  //Z=25315

        /*  perfect orientation  */               /*  abs((cosphi*qx+sinphi*qy)  //Z=25318 */
        else if ( CALC.ordis==ordis_ZDir /*6*/ )
        {/*7*/  //Z=25319
            switch ( CALC.params.orcase )
            {
            case 1:     // General: phi!=0 && phi!=90 && theta!=0 && theta!=90
                //Z=25320
                limql = sqrt(sqr((qx*CALC.cosphi-qy*CALC.sinphi)*CALC.costheta*CALC.sintheta)
                             +sqr(qx+CALC.sinphi*(-qx*CALC.sinphi+qy*CALC.cosphi)*CALC.sintheta*CALC.sintheta)
                             +sqr(qy-CALC.cosphi*(-qx*CALC.sinphi+qy*CALC.cosphi)*CALC.sintheta*CALC.sintheta));  //Z=25321
                pq = CALC.formpq(CALC.params.sigmal,limql,qx,qy,qx*CALC.cosphi*CALC.sintheta,
                                qy*CALC.sinphi*CALC.sintheta,q,CALC.ordis);  //Z=25323
                /* fq = pq;  //Z=25324 */
                if ( CALC.lattice )  //Z=25325
                    fq = CALC.formfq( sqrt(CALC.cosphi*qx*qx+CALC.sinphi*qy*qy+eps9), qx, qy,
                                      qx*CALC.cosphi*CALC.sintheta, qy*CALC.sinphi*CALC.sintheta,q, CALC.ordis );  //Z=25327
                break;   //Z=25328
            case 2:     // X-Axis phi==0 && theta==90
                //Z=25329
                pq = CALC.formpq(CALC.params.sigmal, fabs(qx), qx, qy, qx, 0, q, CALC.ordis);   //Z=25331
                /* fq = pq;  //Z=25332 */
                if ( CALC.lattice )
                    fq = CALC.formfq( fabs(qx), qx, qy, qx, 0, q, CALC.ordis );   //Z=25335
                break;   //Z=25336
            case 3:     // Y-Axis phi==90 && theta==90
                /*Z=24733*/
                pq = CALC.formpq(CALC.params.sigmal, fabs(qy), qx, qy, 0, qy, q, CALC.ordis);   //Z=25339
                /* fq = pq;  //Z=25340 */
                if ( CALC.lattice )
                    fq = CALC.formfq( fabs(qy), qx, qy, 0, qy, q, CALC.ordis );   //Z=25343
                break;   //Z=25344
            case 4:     // Z-Axis (phi==0 || phi==90) && theta==0
                /*Z=24741*/
                pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx, qy, q, CALC.ordis);   //Z=25347
                /* fq = pq;  //Z=25348 */
                if ( CALC.lattice )
                    fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );   //Z=25351
                break;   //Z=25352
            } // switch orcase
        } // if ( CALC.ordis==6 ) //Z=25353

        /*  general orientation  */  //Z=25357
        else if ( CALC.ordis==ordis_Gaussian /*0*/ )
        {   //Z=25358
            switch ( CALC.params.orcase )
            {
            case 1:    /*  general orientation  */  //Z=25359
                pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic,
                                qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis);   //Z=25361
                /* fq = pq;  //Z=25362 */
                if ( CALC.lattice ) /* fq:=pq; */   /*Z=24758*/
                    fq = CALC.formfq( q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic,
                                      qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis );   //Z=25365
                break;   //Z=25366
            case 2:   /*  x-axis  */  //Z=25367
                pq = CALC.formpq(CALC.params.sigmal, sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx, qy, q, CALC.ordis);   //Z=25369
                /* fq = pq;  //Z=25370 */
                //ffq = pq;  //Z=25371
                if ( CALC.lattice )
                    fq = CALC.formfq( sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx, qy, q, CALC.ordis );   //Z=25374
                szq = pq; // ffq;  //Z=25375  TODO das macht hier keinen Sinn
                break;   //Z=25376
            case 3:  /*  y-axis  */  //Z=25377
                pq = CALC.formpq(CALC.params.sigmal, sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx, qy, q, CALC.ordis);   //Z=25379
                /* fq = pq;  //Z=25380 */
                if ( CALC.lattice ) /* fq:=pq; */   /*Z=24776*/
                    fq = CALC.formfq( sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx, qy, q, CALC.ordis );   //Z=25383
                break;   //Z=25384
            case 4:  /*  z-axis  */  //Z=25385
                pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx, qy, q, CALC.ordis);   //Z=25387
                /* fq = pq;  //Z=25388 */
                if ( CALC.lattice )
                    fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );   //Z=25391
                break;   /*Z=24787*/
            } // switch orcase

        } // if ( CALC.ordis==0 ) //Z=25393

    }  /*  of part=1 (Cylinder)  */   //Z=25395


    /* *********** */  //Z=25398
    /* ** disks ** */  //Z=25399
    /* *********** */  //Z=25400
    if ( CALC.partdisk )
    {/*6*/     /*  disks  */  //Z=25401

        /*  isotropic cases  */  //Z=25403
        if ( CALC.ordis == 7 )
        {
            pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);   //Z=25406
            fq = pq;  //Z=25407
            //if ( CALC.lattice )
            //    fq=CALC.formfq(CALC.params.length,CALC.params.radius,CALC.params.sigmal,CALC.params.sigma,CALC.params.p1,
            //                    CALC.params.rho,CALC.params.alphash,CALC.polTheta,CALC.polPhi,q,CALC.params.limq1f,CALC.params.limq2f,
            //                    CALC.params.limq3f,CALC.params.limq4f,CALC.params.limq5f,CALC.params.limq6f,qx,qy,qx,qy,q,CALC.norm,
            //                    CALC.part,CALC.cs,CALC.ordis,CALC.orcase,CALC.params.CR->myarray,CALC.params.CR->carr1f,
            //                    CALC.params.CR->carr2f,CALC.params.CR->carr3f,CALC.params.CR->carr4f,CALC.params.CR->carr5f,
            //                    CALC.params.CR->carr6f,CALC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );
        }   //Z=25411

        /*  perfect orientation  */    //Z=25412
        if ( CALC.ordis==6 ) //Z=25413
        {
            switch ( CALC.params.orcase )
            {
            case 1:
                limql = sqrt(  sqr((qx*CALC.sinphi+qy*CALC.cosphi)*CALC.costheta*CALC.sintheta)
                             + sqr(qx-CALC.cosphi*(qx*CALC.cosphi+qy*CALC.sinphi)*CALC.sintheta*CALC.sintheta)
                             + sqr(qy-CALC.sinphi*(qx*CALC.cosphi+qy*CALC.sinphi)*CALC.sintheta*CALC.sintheta));  //Z=25415
                pq=CALC.formpq(CALC.params.sigmal, limql,
                                qx,qy,qx*CALC.cosphi*CALC.sintheta/q,qy*CALC.sinphi*CALC.sintheta/q,q,CALC.ordis);   //Z=25417
                fq = pq;  //Z=25418
                //if ( CALC.lattice )
                //    fq=CALC.formfq(CALC.params.length,CALC.params.radius,CALC.params.sigmal,CALC.params.sigma,CALC.params.p1,
                //                    CALC.params.rho,CALC.params.alphash,CALC.polTheta,CALC.polPhi,
                //                    sqrt(CALC.sinphi*qx*qx+CALC.cosphi*qy*qy+eps9),
                //                    CALC.params.limq1f,CALC.params.limq2f,CALC.params.limq3f,CALC.params.limq4f,CALC.params.limq5f,CALC.params.limq6f,
                //                    qx,qy,qx*CALC.cosphi*CALC.sintheta/q,qy*CALC.sinphi*CALC.sintheta/q,q,CALC.norm,
                //                    CALC.part,CALC.cs,CALC.ordis,CALC.orcase,CALC.params.CR->myarray,CALC.params.CR->carr1f,
                //                    CALC.params.CR->carr2f,CALC.params.CR->carr3f,CALC.params.CR->carr4f,CALC.params.CR->carr5f,
                //                    CALC.params.CR->carr6f,CALC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );  //Z=25421
                break;   //Z=25422
            case 2:
                pq=CALC.formpq(CALC.params.sigmal, fabs(qy), qx,qy,qx/q,0,q,CALC.ordis);   //Z=25425
                fq = pq;  //Z=25426
                //if ( CALC.lattice )
                //    fq=CALC.formfq(CALC.params.length,CALC.params.radius,CALC.params.sigmal,CALC.params.sigma,CALC.params.p1,
                //                    CALC.params.rho,CALC.params.alphash,CALC.polTheta,CALC.polPhi, fabs(qy),
                //                    CALC.params.limq1f,CALC.params.limq2f,CALC.params.limq3f,CALC.params.limq4f,CALC.params.limq5f,CALC.params.limq6f,
                //                    qx,qy,qx/q,0,q,CALC.norm,
                //                    CALC.part,CALC.cs,CALC.ordis,CALC.orcase,CALC.params.CR->myarray,CALC.params.CR->carr1f,
                //                    CALC.params.CR->carr2f,CALC.params.CR->carr3f,CALC.params.CR->carr4f,CALC.params.CR->carr5f,
                //                    CALC.params.CR->carr6f,CALC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );   //Z=25429
                break;  //Z=25430
            case 3:
                pq=CALC.formpq(CALC.params.sigmal, fabs(qx), qx,qy,0,qy/q,q,CALC.ordis);  //Z=25433
                fq = pq;  //Z=25434
                //if ( CALC.lattice )
                //    fq=CALC.formfq(CALC.params.length,CALC.params.radius,CALC.params.sigmal,CALC.params.sigma,CALC.params.p1,
                //                    CALC.params.rho,CALC.params.alphash,CALC.polTheta,CALC.polPhi, fabs(qx),
                //                    CALC.params.limq1f,CALC.params.limq2f,CALC.params.limq3f,CALC.params.limq4f,CALC.params.limq5f,CALC.params.limq6f,
                //                    qx,qy,0,qy/q,q,CALC.norm,
                //                    CALC.part,CALC.cs,CALC.ordis,CALC.orcase,CALC.params.CR->myarray,CALC.params.CR->carr1f, /*TODO: hier stand carr1p*/
                //                    CALC.params.CR->carr2f,CALC.params.CR->carr3f,CALC.params.CR->carr4f,CALC.params.CR->carr5f,
                //                    CALC.params.CR->carr6f,CALC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );  //Z=25437
                break;  //Z=25438
            case 4:
                pq=CALC.formpq(CALC.params.sigmal, q, qx,qy,qx,qy,q,CALC.ordis);   //Z=25441
                fq = pq;  //Z=25442
                //if ( CALC.lattice )
                //    fq=CALC.formfq(CALC.params.length,CALC.params.radius,CALC.params.sigmal,CALC.params.sigma,CALC.params.p1,
                //                    CALC.params.rho,CALC.params.alphash,CALC.polTheta,CALC.polPhi, q,
                //                    CALC.params.limq1f,CALC.params.limq2f,CALC.params.limq3f,CALC.params.limq4f,CALC.params.limq5f,CALC.params.limq6f,
                //                    qx,qy,qx,qy,q,CALC.norm,
                //                    CALC.part,CALC.cs,CALC.ordis,CALC.orcase,CALC.params.CR->myarray,CALC.params.CR->carr1f, /*TODO ...*/
                //                    CALC.params.CR->carr2f,CALC.params.CR->carr3f,CALC.params.CR->carr4f,CALC.params.CR->carr5f,
                //                    CALC.params.CR->carr6f,CALC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );   //Z=25445
                break;  //Z=25446
            } // switch orcase
        } /* if ( CALC.ordis==6 ) */  //Z=25447

        /*  isotropic fraction  */  //Z=25449
        if ( CALC.iso>0 )
            pqiso = CALC.formpq(CALC.params.sigmal, q, qx,qy,qx,qy,q, 7/*ordis*/);  //Z=25451
        else
            pqiso = 0.0;  //Z=25452

        /*  general orientation  */  //Z=25455
        if ( CALC.ordis==0 )
        {
            // TODO: bei den folgenden 4 Aufrufen von formfq() stehen die Arrays carr1p anstatt carr1f im Aufruf.
            //       Nach einem Gespräch mit Prof. Förster: das ist wohl ein Copy/Paste-Fehler
#ifndef __CUDACC__
            //if ( idxCenter ) qDebug() << "TODO formfq, orcase:"<<CALC.params.orcase
            //             << "ordis:"<<CALC.ordis << "part:"<<CALC.params.part << "cs:"<<CALC.params.cs;
            //Debug: TODO formfq, orcase: 1 ordis: 0 part: 2 cs: 0
#endif
            switch ( CALC.params.orcase )
            {
            case 1:
                pq=CALC.formpq(CALC.params.sigmal,
                                q,qx,qy,qx*CALC.cosphic/q-qy*CALC.sinphic/q,qx*CALC.sinphic/q+qy*CALC.cosphic/q,q,CALC.ordis);   //Z=25459
                /* fq = pq;  //Z=25460 */
                if ( CALC.lattice )
                    fq=CALC.formfq( q, qx, qy, qx*CALC.cosphic/q-qy*CALC.sinphic/q, qx*CALC.sinphic/q+qy*CALC.cosphic/q, q, CALC.ordis );   //Z=25463
                break;  //Z=25464
            case 2:  //Z=25465
                pq=CALC.formpq(CALC.params.sigmal, sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx,qy,qx/q,qy/q,q,CALC.ordis);
                /* fq = pq;  //Z=25468 */
                if ( CALC.lattice ) //fq:=pq;
                    fq=CALC.formfq( sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx/q, qy/q, q, CALC.ordis );  //Z=25471
                break;  //Z=25472
            case 3:  //Z=25473
                pq=CALC.formpq(CALC.params.sigmal,sqrt(qx*qx+(1-CALC.order)*qy*qy),qx,qy,qx/q,qy/q,q,CALC.ordis );  //Z=25475
                /* fq = pq;  //Z=25476 */
                if ( CALC.lattice )
                    fq=CALC.formfq( sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx/q, qy/q, q, CALC.ordis );  //Z=25479
                break;  //Z=25480
            case 4:  //Z=25481
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25483
                /* fq = pq;  //Z=25484 */
                if ( CALC.lattice ) //fq:=pq;
                    fq=CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );  //Z=25487
                break;  //Z=25488
            } // switch orcase
        } // if ( CALC.ordis==0 )  //Z=25489

    } /*  of part=2, disk  */  //Z=25490


    if ( CALC.partvesicle )
    {  //Z=25493
        pq = CALC.polyvesicle(/*length,radius,sigma,sigmal,*/q);  //Z=25494
        //if ( CheckBoxf2q.Checked==true )
        //    fq = f2dpolyvesicle(length,radius,sigma,sigmal,q);  //Z=25495
        //else
        fq = pq;  //Z=25496
    }  //Z=25497


    // Wird aktuell nicht verwendet
    //if ( CALC.partliposome )
    //{  //Z=25499
    //    /* pq:=polyliposome(length,radius,sigma,sigmal,shellno,alphash,ceff,reff,a,b,c,domainsize,aziwidth,3,q);  //Z=25500 */
    //    /* pq:=liposome(a,rho,radiusi,p1,shellno,q);  //Z=25501 */
    //    pq = liposome(5,0,10,0.8,3,q);  //Z=25502
    //    fq = pq;  //Z=25503
    //}  //Z=25504


    if ( CALC.partcube )
    {/*6*/  //Z=25506
        /* pq:=polycube(radius,sigma,0,q);  //Z=25507 */
        /* fq:=polycube(radius,sigma,1,q);  //Z=25508 */

        /*  isotropic cases  */  //Z=25511
        if ( (CALC.ordis==7) && CALC.homogeneous )
        {/*7*/  //Z=25512
            pq = CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25514
            fq = pq;  //Z=25515
            /* if lattice then  //Z=25516 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25517 */
            /*      part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25518 */
        }/*7*/  //Z=25519

        /*  perfect orientation  */  //Z=25521
        if ( (CALC.ordis==6) && CALC.homogeneous )
        {/*7*/  //Z=25522
            pq = CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25524
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25525 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25526 */
            fq = pq;  //Z=25527
        }/*7*/  //Z=25528
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^);  //Z=25529 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^);  //Z=25530 */

    }/*6*/  //Z=25532 cbpartCube


    if ( CALC.partellipsoid )
    {   //Z=25534

        /*  isotropic cases  */  //Z=25535
        if ( CALC.ordis==7 && CALC.homogeneous )
        {   /*Z221118=25536*/
            pq=CALC.formpq(CALC.epsilon,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25538
            fq = pq;  //Z=25539
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25540 */
            /*      part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25541 */
        }  //Z=25542

        /* perfect orientation */
        if ( CALC.ordis==6 && CALC.homogeneous )
        {   //Z=25545
            switch ( CALC.params.orcase )
            {
            case 1:
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25548
                break;
            case 2:
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,0,q,CALC.ordis);  //Z=25551
                break;
            case 3:
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qy,qx,q,CALC.ordis);  //Z=25554
                break;
            case 4:
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25557
                break;
            }
            fq = pq;  //Z=25559
        }

        /*  general orientation  */  //Z=25562
        if ( CALC.ordis==0 )
        {   //Z=25563
            switch ( CALC.params.orcase )
            {
            case 1:   /*  general orientation  */  //Z=25564
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx*CALC.cosphic-qy*CALC.sinphic,qx*CALC.sinphic+qy*CALC.cosphic,q,CALC.ordis);  //Z=25566
                /* fq = pq;  //Z=25567 */
                if ( CALC.lattice )
                    fq=CALC.formfq( q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic, qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis );  //Z=25570
                break;  //Z=25571
            case 2:  /*  x-axis  */  //Z=25572
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25574
                fq = pq;  //Z=25575
                /* ffq:=pq;  //Z=25576 */
                /* if lattice then //fq:=pq;  //Z=25577 */
                /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1f,limq2f,limq3f,limq4f,limq5f,limq6f,qx,qy,qx,qy,q,norm,  //Z=25578 */
                /*    part,cs,ordis,orcase,myarray,carr1p^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25579 */
                /* szq:=ffq;  //Z=25580 */
                break;  //Z=25581
            case 3:  /*  y-axis  */  //Z=25582
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25584
                /* fq = pq;  //Z=25585 */
                if ( CALC.lattice )
                    fq=CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );  //Z=25588
                break;  //Z=25589
            case 4:  /*  z-axis  */  //Z=25590
                pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25592
                /* fq = pq;  //Z=25593 */
                if ( CALC.lattice )
                    fq=CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );  //Z=25596
                break;  //Z=25597
            }

            //fq = pq;  //Z=25600  TODO  Wenn if ( CALC.lattice ) greift, dann macht die dortige Rechnung keinen Sinn mehr

        } // if ordis==0  //Z=25598
    } // partellipsoid  //Z=25601


    if ( CALC.parttriellipsoid )
    {   //Z=25603

        /*  isotropic cases  */  //Z=25604
        if ( (CALC.ordis==7) && CALC.homogeneous )
        {   //Z=25605
            pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25607
            //fq = pq;  //Z=25608
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25609 */
            /*      part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25610 */
        }  //Z=25611

        /*  perfect orientation  */  //Z=25613
        if ( (CALC.ordis==6) && CALC.homogeneous )
        {   //Z=25614
            pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);
            //fq = pq;  //Z=25617
            /* if lattice then  //Z=25618 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25619 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25620 */
        }  //Z=25621
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^,carr11pm^,carr22pm^);  //Z=25622 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^,carr11pm^,carr22pm^);  //Z=25623 */

        fq = pq;  //Z=25625
    }  //Z=25626


    if ( CALC.partbarrel )
    {   //Z=25629

        /*  isotropic cases  */  //Z=25630
        if ( (CALC.ordis==7) && CALC.homogeneous )
        {   //Z=25631
            pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25633
            //fq = pq;  //Z=25634
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25635 */
            /*      part,cs,ordis,orcase,myarray,carr1f^ ,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25636 */
        }  //Z=25637

        /*  perfect orientation  */  //Z=25639
        if ( (CALC.ordis==6) && CALC.homogeneous )
        {   //Z=25640
            pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25642
            //fq = pq;  //Z=25643
            /* if lattice then  //Z=25644 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25645 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25646 */
        }  //Z=25647
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^,carr11pm^,carr22pm^);  //Z=25648 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^,carr11pm^,carr22pm^);  //Z=25649 */

        fq = pq;  //Z=25651
    }  //Z=25652


    if ( CALC.partball )
    {   //Z=25654

        /*  isotropic cases  */  //Z=25655
        if ( (CALC.ordis==7) && CALC.homogeneous )
        {   //Z=25656
            pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25658
            //fq = pq;  //Z=25659
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25660 */
            /*      part,cs,ordis,orcase,myarray,carr1f^ ,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25661 */
        }  //Z=25662

        /*  perfect orientation  */  //Z=25664
        if ( (CALC.ordis==6) && CALC.homogeneous )
        {   //Z=25665
            pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25667
            //fq = pq;  //Z=25668
            /* if lattice then  //Z=25669 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25670 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25671 */
        }  //Z=25672
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^,carr11pm^,carr22pm^);  //Z=25673 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^,carr11pm^,carr22pm^);  //Z=25674 */

        fq = pq;  //Z=25676
    }  //Z=25677


#ifdef nochnichteingebaut
    /*  excluded volume chain  */  //Z=25679
    /*  descloi4(rg,nu,gammas,q)  */  //Z=25680
    if ( partchain )
    {  //Z=25681
        pq = descloi4(radius,sigma,length,q);  //Z=25682
        fq = pq;  //Z=25683
    }  //Z=25684

    if ( partkpchain )
    {  //Z=25686
        pq = kpchain(length,sigmal,q)*psphered(radius,sigma,2,q);  ;      /*  sigmal = persistence length  */  //Z=25687
        fq = pq;  //Z=25688
    }  //Z=25689


    if ( Checkbox10 )
    {  //Z=25692
        if ( partsphere ) pq = schulz(length,radius,radiusi,1,1,sigma,q,1);  //Z=25693
        if ( partcylinder ) pq = schulztheta(length,radius,radiusi,1,sigma,q,2);  //Z=25702
        if ( partellipsoid ) pq = schulztheta(length,radius,radiusi,1,sigma,q,4);  //Z=25703
        if ( parttriellipsoid ) pq = schulzthetaphi(length,radius,radiusi,sigma,q,5);  //Z=25718
        if ( partcube ) pq = schulzthetaphi(length,radius,radiusi,sigma,q,7);  //Z=25719
    }  //Z=25720
#endif


    double width = 1;
    if ( CALC.RadioButtonDebyeScherrer )
        width = 4.0/CALC.params.domainsize;  //Z=25722
    if ( CALC.RadioButtonPara )
        width = (4.0/CALC.params.domainsize)+sqr(CALC.reldis)*CALC.dist*sqr(q);  //Z=25723
    if ( CALC.lattice && CALC.shp==cbpeakAnisotropicGaussian/*8*/ && CALC.ordis==7/*isotropic*/ )
        width = CALC.params.sig.length() /*sqrt(sigx*sigx+sigy*sigy+sigz*sigz)*/ /3.0;  //Z=26118
#ifndef __CUDACC__
    if ( idxCenter && width != 1 && CALC.params.width_zuf != 1 )
        qDebug() << "WIDTH" << width << CALC.params.width_zuf;
#endif

    //if ( ! CALC.tpvRandomOldValues.empty() )   // doIntCalc... - nicht auf der GPU zulässig
    {
        //if ( fabs(qx) < 0.1 && fabs(qy) < 0.1 && fabs(qz) < 0.1 )
        //    std::cerr << "TPV CALC " << CALC.params.width_zuf << std::endl << std::flush;
        width = width * CALC.params.width_zuf;
    }

    /*  Percus-Yevick  */  //Z=25725
    if ( CALC.ltype==20 )
    {/*6*/  //Z=25726
        pq = pq * CALC.spy( q );  //Z=25727
    }/*6*/  //Z=25728

#ifdef nochnichteingebaut
    /*  Teubner-Strey  */  //Z=25730
    if ( ltype==21 )
    {/*6*/  //Z=25731
        tsrr = sqr(2*M_PI*domainsize/uca);  //Z=25732
        tsc0 = sqr(1+tsrr);  //Z=25733
        tsc1 = 2*sqr(domainsize)*(1-tsrr);  //Z=25734
        tsc2 = pow(domainsize,4);  //Z=25735
        pq = tsc0/(tsc0+tsc1*q*q+tsc2*q*q*q*q);  //Z=25736
    }/*6*/  //Z=25737
#endif

/*#ifndef __CUDACC__
    if ( idxCenter )
    {
        if ( CALC.RadioButtonDebyeScherrer )
            qDebug() << "WIDTH" << width << "Scherrer" << CALC.params.domainsize << "zuf" << CALC.params.width_zuf;
        if ( CALC.RadioButtonPara )
            qDebug() << "WIDTH" << width << "Para" << CALC.params.domainsize << CALC.reldis << CALC.dist << q << "zuf" << CALC.params.width_zuf;
        if ( ! CALC.RadioButtonDebyeScherrer && ! CALC.RadioButtonPara )
            qDebug() << "WIDTH" << width << "zuf" << CALC.params.width_zuf;
    }
#endif*/

    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

    double /*qhkl0,*/ qxhkl, qyhkl, qzhkl, qhkl, qxhklt, qyhklt, qzhklt, qhklt, g3, g3t, peaknorm1; //, peaknorm2;
    double dqx, dqy, dqz, yphi, psiord, phiord, x2phi;

    double widthiso = 1.0 / CALC.params.uca;  //ZN=20321
    double sphno=0, cubevol=0; // fürs Debuggen besser

    radintensity = 0.0;     // Immer nur für ein Pixel berechnet, somit kein Array nötig
    intensity = 0.0;

    /* ** lattice hkl-factor calcuation  */  //Z=25744
    if ( CALC.lattice )
    {/*6*/  //Z=25745

        /*  PC-version  */  //Z=25748
        sphno = CALC.latpar[1];  //Z=25749
        cubevol = CALC.latpar[2];  //Z=25750
        //global: dwfactor = CALC.latpar[3];  //Z=25751

        /* isotropic peaks */
        for ( int ii=1; ii<=CALC.peakmax1; ii++ )
        {   //Z=25759  //ZN=25967
            if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat
            //int h    = CALC.latpar1(ii,1);  //Z=25760
            //int k    = CALC.latpar1(ii,2);
            //int l    = CALC.latpar1(ii,3);
            int mhkl    = CALC.latpar1(ii,4);  //Z=25763
            int fhkl    = CALC.latpar1(ii,5);
            qhkl        = CALC.latpar3(ii,5);
            double x2   = sqr(q-qhkl)/sqr(widthiso);  //Z=25766
            double sq   = exp(-4*x2/M_PI)/(M_PI*widthiso/2.0);  //Z=25767
            radintensity += sq*mhkl*fhkl;  //Z=25768
        } /* of peak loop */  //Z=25769

#ifndef __CUDACC__
        int ignHklData = 0;
#endif
        for ( int ii=1; ii<=CALC.peakmax2; ii++ )
        {   //Z=25779  //ZN=25987
            if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat
            int h = CALC.latpar2(ii,1);  //Z=25780  //ZN=25988
            int k = CALC.latpar2(ii,2);
            int l = CALC.latpar2(ii,3);
            int mhkl = CALC.latpar2(ii,4);
            int fhkl = CALC.latpar2(ii,5);

            //qhkl0 = CALC.latpar3(ii,1);  //Z=25820  //ZN=26028
            qxhkl = CALC.latpar3(ii,2);
            qyhkl = CALC.latpar3(ii,3);
            qzhkl = CALC.latpar3(ii,4);
            qhkl  = CALC.latpar3(ii,5);
            qxhklt = CALC.latpar3(ii,7);
            qyhklt = CALC.latpar3(ii,8);
            qzhklt = CALC.latpar3(ii,9);
            g3     = CALC.latpar3(ii,10);
            qhklt  = CALC.latpar3(ii,14);  //ZN=26037
            g3t    = CALC.latpar3(ii,15);  //ZN=26038

            switch ( CALC.shp )
            {
            case cbpeakAnisotropicGaussian /*8*/:  //Z=25853
                /*  perfect orientation  */  //Z=25854
                if ( CALC.ordis==6/*z-dir*/ || fabs(CALC.params.dbeta) < eps9/*dbeta==0*/ )
                {   //Z=25855
                    D8( qDebug() << "CPU: shp=8(A): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 );
                    dqx = qx-qxhkl;
                    dqy = qy-qyhkl;
                    dqz = qz-qzhkl;
                    dqs1 = (dqx*CALC.params.ax1.x()+dqy*CALC.params.ax1.y()+dqz*CALC.params.ax1.z())/(CALC.ax1n_sigx);
                    dqs2 = (dqx*CALC.params.ax2.x()+dqy*CALC.params.ax2.y()+dqz*CALC.params.ax2.z())/(CALC.ax2n_sigy);
                    dqs3 = (dqx*CALC.params.ax3.x()+dqy*CALC.params.ax3.y()+dqz*CALC.params.ax3.z())/(CALC.ax3n_sigz);
                    x2 = dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                    sq = exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8.0);
                }
                else if ( CALC.ordis==7/*isotropic*/ )
                {   /* isotropic orientation */   //Z=25890
                    D8( qDebug() << "CPU: shp=8(B): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 );
                    x2 = (q-qhkl)*(q-qhkl)/(width*width);  //Z=25893
                    sq = g3*exp(-4*x2/M_PI)/(M_PI*width/2.0);  //Z=25894
                }
                else if ( CALC.ordis==13/*fiber pattern*/ )
                {   /* fiber pattern */
                    D8( qDebug() << "CPU: shp=8(C): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 );
                    // Update 20220608
                    // rotaxphi   = CALC.polPhi
                    // rotaxtheta = CALC.polTheta
#define sign(x) ((x<0)?-1:1)
                    /* rotated around y-axis */
                    if ( fabs(CALC.params.polPhi)<eps9 /*(rotaxphi=0)*/ &&
                         fabs(CALC.params.polTheta-90)<eps9 /*(rotaxtheta=90)*/ )
                    {  //Z=25899
                        double signq=sign(qyhkl);
                        qyhkl=signq*sqrt(qyhkl*qyhkl+qzhkl*qzhkl);  //Z=25901
                        qzhkl=1E-20;
                        dqx=qx-qxhkl;
                        dqy=qy-qyhkl;
                        dqz=qz-qzhkl;
                        dqs1=(dqx*CALC.params.ax1.x()+dqy*CALC.params.ax1.y()+dqz*CALC.params.ax1.z())/(CALC.ax1n_sigx);  //Z=25906
                        dqs2=(dqx*CALC.params.ax2.x()+dqy*CALC.params.ax2.y()+dqz*CALC.params.ax2.z())/(CALC.ax2n_sigy);
                        dqs3=(dqx*CALC.params.ax3.x()+dqy*CALC.params.ax3.y()+dqz*CALC.params.ax3.z())/(CALC.ax3n_sigz);
                        x2=dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                        sq=exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8.);
                    }  //Z=25911
                    /* rotated around x-axis */
                    else if ( fabs(CALC.params.polPhi-90)<eps9 /*(rotaxphi=90)*/ &&
                              fabs(CALC.params.polTheta-90)<eps9 /*(rotaxtheta=90)*/ )
                    {  //Z=25913
                        double signq=sign(qxhkl);
                        qxhkl=signq*sqrt(qxhkl*qxhkl+qzhkl*qzhkl);
                        qzhkl=1E-20;
                        dqx=qx-qxhkl;
                        dqy=qy-qyhkl;
                        dqz=qz-qzhkl;
                        dqs1=(dqx*CALC.params.ax1.x()+dqy*CALC.params.ax1.y()+dqz*CALC.params.ax1.z())/(CALC.ax1n_sigx);  //Z=25920
                        dqs2=(dqx*CALC.params.ax2.x()+dqy*CALC.params.ax2.y()+dqz*CALC.params.ax2.z())/(CALC.ax2n_sigy);
                        dqs3=(dqx*CALC.params.ax3.x()+dqy*CALC.params.ax3.y()+dqz*CALC.params.ax3.z())/(CALC.ax3n_sigz);
                        x2=dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                        sq=exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8.);
                    }  //Z=25925
                    /* rotated round an oblique axis */
                    else // if ((rotaxphi<>0) and (rotaxphi<>90)) then begin
                    {/*10*/  //Z=25927
                        CALC.qrombchid(CALC.params.length,CALC.ucl1,CALC.params.p1,CALC.ucl2,CALC.params.alpha,/*CALC.ucl3, (war dbeta)*/
                                      delta,CALC.params.polTheta*M_PI/180.0,CALC.params.polPhi*M_PI/180.0,qx,qy,qz,CALC.ri11,CALC.ri12,CALC.ri13,
                                      CALC.ri21,CALC.ri22,CALC.ri23,CALC.ri31,CALC.ri32,CALC.ri33,  //Z=25928
                                      qxhkl,qyhkl,qzhkl,qhkl,
                                      CALC.params.ax1.length(),CALC.params.ax2.length(),CALC.params.ax3.length(),
                                      CALC.params.ax1.x(),CALC.params.ax1.y(),CALC.params.ax1.z(),
                                      CALC.params.ax2.x(),CALC.params.ax2.y(),CALC.params.ax2.z(),
                                      CALC.params.ax3.x(),CALC.params.ax3.y(),CALC.params.ax3.z(),
                                      CALC.params.sig.x(),CALC.params.sig.y(),CALC.params.sig.z(),  //Z=25929
                                      CALC.ordis,3,5,7,h,k,l, CALC.params.CR->carr1p, sq);  //Z=25930
                        sq = sq*2*M_PI*qhkl/(2.0*M_PI*sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8.0);  //Z=25931
                    }/*10*/  //Z=25932
                }
                else // if ( (CALC.ordis!=6) && (CALC.ordis!=7) && (CALC.ordis!=13) && (CALC.dbeta!=0) )      //Z=25985
                {   /* other anisotropic cases, rotated around the qhkl-axis direction */                                                             //20210812-D
                    D8( qDebug() << "CPU: shp=8(D): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 );
                    double phi = atan2(qyhkl,(qxhkl+eps6));
                    double theta = atan2(sqrt(qxhkl*qxhkl+qyhkl*qyhkl),(qzhkl+eps6));  //Z=25987
                    phi = phi*180/M_PI;
                    theta = theta*180/M_PI;
                    CALC.qrombdeltac(CALC.params.p1, CALC.params.sigma, CALC.params.alpha,
                                    theta, phi, qx, qy, qz,
                                    qxhkl, qyhkl, qzhkl, qhkl,
                                    CALC.params.ax1.length(), CALC.params.ax2.length(), CALC.params.ax3.length(),
                                    CALC.params.ax1.x(), CALC.params.ax1.y(), CALC.params.ax1.z(),
                                    CALC.params.ax2.x(), CALC.params.ax2.y(), CALC.params.ax2.z(),
                                    CALC.params.ax3.x(), CALC.params.ax3.y(), CALC.params.ax3.z(),
                                    CALC.params.sig.x(), CALC.params.sig.y(), CALC.params.sig.z(),
                                    CALC.ordis,3, /*i0=*/5, /*i1=*/6,0,0,0, CALC.params.CR->carr1p, sq );
                    sq = sq*2*M_PI*qhkl/CALC.params.norm;  //Z=25991
                }

                psiord = 1;  //Z=25994
#ifdef undef
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {   //ZN=26205
                    TODO: im Prinzip nochmal die gleichen Verzweigungen wie oben
                }   //ZN=26270
#endif
                break;  // shp == cbpeakAnisotropicGaussian  //Z=26062

            case cbpeakLorentzian /*1*/:
                peaknorm1 = CALC.latpar3(ii,11);  //Z=26066
                yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = sqrt(CALC.c1)/(M_PI*width*(1+CALC.c1*x2));
                x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** b-factor ***/
                psiord = g3/(peaknorm1*(1+x2phi*yphi*yphi));  //Z=26071
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {
                    yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*(1+x2phi*yphi*yphi));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26076
                break;

            case cbpeakGaussian /*2*/:                                                          //20210812-E
                peaknorm1 = CALC.latpar3(ii,11);  //Z=26081
                yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = exp(-4*x2/M_PI)/(M_PI*width/2.0);
                x2phi = 4*q*q/(M_PI*sqr(CALC.phiwidth));             /*** a-factor ***/
                psiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
#ifndef __CUDACC__
                if ( isinf(psiord) ) qDebug() << "psiord1" << g3 << x2phi << yphi << peaknorm1;
#endif
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {
                    double peaknorm1t = CALC.latpar3(ii,16);   //ZN=26305
                    yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhklt));
                    phiord = g3t*exp(-x2phi*yphi*yphi)/peaknorm1t;
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;   //ZN=26313
#ifndef __CUDACC__
                    if ( isinf(psiord) )
                        qDebug() << "psiord2" << CALC.params.ceff << phiord << "=" << g3t << x2phi << yphi << peaknorm1t;

                    if ( idxCenter && (ii < 10 || ii > CALC.peakmax2-10) )
                        qDebug() << ii << "q?hkl" << qxhkl << qyhkl << qzhkl << "q?hklt" << qxhklt << qyhklt << qzhklt;
#endif
                }
                // TODO Hier tauchen Probleme bim BCC / BCT auf:
                // Wenn TwRatio(ceff)=0 und CheckBoxTwinned=True dann wird falsch gerechnet. Bei allen anderen
                // Kombinationen stimmen die Ergebnisse.
                // 230804: CheckBoxTwinned soll lt. Hr. Förster langfristig rausfliegen. Die Parts mit "if twin then" sind
                //          schon angepasst, die anderen müssen noch überarbeitet werden.
                break;

            case cbpeakMod1Lorentzian /*Lorentzian1*/:
                peaknorm1 = CALC.latpar3(ii,11);  //Z=26100
                yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = 2*sqrt(CALC.c2)/(M_PI*width*(1+CALC.c1*x2)*(1+CALC.c1*x2));
                x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*(1+x2phi*yphi*yphi)*(1+x2phi*yphi*yphi));
                if ( CALC.CheckBoxTwinned )
                {  //Z=26106
                    yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*(1+x2phi*yphi*yphi)*(1+x2phi*yphi*yphi));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26110
                break;

            case cbpeakMod2Lorentzian /*Lorentzian2*/:
                peaknorm1 = CALC.latpar3(ii,11);  //Z=26115
                yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = sqrt(CALC.c3)/(2*width*exp(3*log(1+CALC.c1*x2)/2.0));
                x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*exp(3*log(1+x2phi*yphi*yphi)/2.0));
                if ( CALC.CheckBoxTwinned )
                {  //Z=26121
                    yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*exp(3*log(1+x2phi*yphi*yphi)/2.0));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26125
                break;

            case cbpeakPseudoVoigt /*Voigt*/:
            {
                peaknorm1 = CALC.latpar3(ii,11);  //Z=26131
                double peaknorm2 = CALC.latpar3(ii,12);
                yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = CALC.eta*sqrt(CALC.c1)/(M_PI*width*(1+CALC.c1*x2))+(1-CALC.eta)*sqrt(CALC.c0)*exp(-CALC.c0*x2)/(sqrt(M_PI)*width);  //Z=26135
                double x2psi = 4*q*q/(M_PI*sqr(CALC.phiwidth));             /*** a-factor ***/
                //double x2psihkl = 4*qhkl*qhkl/(M_PI*sqr(CALC.phiwidth));  //Z=26137 wird nicht weiter verwendet
                x2phi = 4*q*q/sqr(CALC.phiwidth);                /*** b-factor ***/
                psiord = g3*(CALC.eta*(1/(1+x2phi*yphi*yphi))+(1-CALC.eta)*exp(-x2psi*yphi*yphi))/(CALC.eta*peaknorm1+(1-CALC.eta)*peaknorm2);  //Z=26139
                if ( CALC.CheckBoxTwinned )
                {  //Z=26140
                    yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3*(CALC.eta*(1/(1+x2phi*yphi*yphi))+(1-CALC.eta)*exp(-x2psi*yphi*yphi))/(CALC.eta*peaknorm1+(1-CALC.eta)*peaknorm2);
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26144
                break;
            }

            case cbpeakPearsonVII /*Pearson*/:
                peaknorm1 = CALC.latpar3(ii,11);  //Z=26149
                yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = CALC.gamma(CALC.beta)*CALC.c4*2/(CALC.gamma(CALC.beta-0.5)*M_PI*width*exp(CALC.beta*log(1+4*CALC.c4*x2)));
                x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*exp(CALC.beta*log(1+x2phi*yphi*yphi)));
                if ( CALC.CheckBoxTwinned )
                {  //Z=26155
                    yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*exp(CALC.beta*log(1+x2phi*yphi*yphi)));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26159
                break;

#ifdef undef
                if ( BurgerG )
                {/*8*/  //Z=26161
                    /* x2phihkl:=4*qhkl*qhkl/(pi*phiwidth*phiwidth);  //Z=26162 */
                    /* peaknorm1:=gaussnorm3(x2phihkl);  //Z=26163 */
                    peaknorm1 = latpar3p^[ii][11];  //Z=26164
                    yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));  //Z=26165
                    x2 = (q-qhkl)*(q-qhkl)/(width*width);  //Z=26166
                    sq = burger(width,bnu,x2);  //Z=26167
                    x2phi = 4*q*q/(M_PI*phiwidth*phiwidth);             /* ** a-factor ** */  //Z=26168
                    psiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;  //Z=26169
                    if ( CheckBoxTwinned.Checked==true )
                    {/*9*/  //Z=26170
                        yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));  //Z=26171
                        phiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;  //Z=26172
                        psiord = ceff*psiord+(1-ceff)*phiord;  //Z=26173
                    }/*9*/  //Z=26174
                }/*8*/  //Z=26175
#endif

            default:
                return 0;
            } // switch shp

            if ( qhkl > 0 )
            {
                intensity += sq*mhkl*fhkl*psiord;  //Z=26176
                if ( isinf(intensity) )
                {
#ifndef __CUDACC__
                    qDebug() << "INTENS=INF" << sq << mhkl << fhkl << psiord << "ii" << ii << CALC.peakmax2;
#endif
                    break;
                }
            }
#ifndef __CUDACC__
            else ignHklData++;
#endif

        }/*2*/  /* of peak-loop */  //Z=26182

#ifndef __CUDACC__
        if ( ignHklData > 0 && idxCenter )
            qDebug() << "IgnHklData" << ignHklData << "max" << CALC.peakmax2;
#endif

        szqiso = (1+(1)*(8*M_PI*M_PI*M_PI*radintensity/(4*M_PI*sphno*cubevol)-1)*exp(-CALC.params.ucb/*dwfactoriso*/*q*q));   /*Z0311=24804*/
        szq = (1+(fq/pq)*(2*M_PI*intensity/(sphno*cubevol)-1)*exp(-CALC.dwfactor*q*q));   /*Z0311=24805*/
    }
    else // (if lattice)
    {  //Z=26200
        szq = 1.0;
        szqiso = 1.0;   // zur Sicherheit
    }

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    //Z=26207: szq*pq + CALC.iso*szqiso*pqiso
    //Z=30277: xyintensity^[ihex+zmax][i+zmax] = base+izero*xyintensity^[ihex+zmax][i+zmax]+ifluc/(1.0+q*q*rfluc*rfluc);  //Z=30277
    // retval ist der Pixelwert bei [ihex+zmax][i+zmax]

#ifndef __CUDACC__
    if ( retval < -1e8 || isnan(retval) || isinf(retval) )
        qDebug() << "szq"<<szq << "pq"<<pq <<"fq"<<fq << "szqiso"<<szqiso << "pqiso"<<pqiso << "q"<<q
                 << "intens"<<intensity << "radint"<<radintensity << CALC.peakmax1
                 << "="<<retval;
    // szq inf pq 2.16215e-05 fq 2.08036e-07 szqiso 1 pqiso 2.16215e-05 q 2.74004 intens inf radint 2.45523 107 = inf

#else
    //if ( retval < -1e6 )
    //    printf( "szq=%lf pq=%lf fq=%lf szqiso=%lf pqiso=%lf q=%lf intens=%lf radint=%lf erg=%lf\n", szq, pq, fq, szqiso, pqiso, q, intensity, radintensity, retval );
#endif

    return retval;
} /* doIntCalc_GENERIC_q_xyz() */






void SasCalc_GENERIC_calculation::corotations(double a, double b, double c, double alpha, double beta, double gamma,
                                              double u, double v, double w, double ephi, /*int ltype,*/
                                              bool lat1d, bool lat2d, bool /*lat3d*/, /*bool twin,*/
                                              double &m11, double &m12, double &m13, double &m21, double &m22, double &m23, double &m31, double &m32, double &m33,
                                              double &mtw11, double &mtw12, double &mtw13, double &mtw21, double &mtw22, double &mtw23, double &mtw31, double &mtw32, double &mtw33, double &vvol,
                                              double &nuvwx, double &nuvwy, double &nuvwz, double &uuvwx, double &uuvwy, double &uuvwz, double &vuvwx, double &vuvwy, double &vuvwz,
                                              double &nhklx, double &nhkly, double &nhklz, double &uhklx, double &uhkly, double &uhklz, double &vhklx, double &vhkly, double &vhklz )
{
    // int i;  //Z=9224
    double ca, cb, cg, /*sa,*/ sb, sg, /*vol,*/ n1, n2, n3, l1, l2, l3, m1, m2, m3;  //Z=9225
    double msi11, msi12, msi13, msi21, msi22, msi23, msi31, msi32, msi33; //, detmsi;  //Z=9226
    double /*ms11, ms12, ms13, ms21, ms22, ms23, ms31, ms32, ms33,*/ lenn;  //Z=9227
    //double msi11n, msi12n, msi13n, msi21n, msi22n, msi23n, msi31n, msi32n, msi33n, detmsin;  //Z=9228
    double m11n, m12n, m13n, m21n, m22n, m23n, m31n, m32n, m33n;  //Z=9229
    //double mst11, mst12, mst13, mst21, mst22, mst23, mst31, mst32, mst33, detmst;  //Z=9230
    //double mt11, mt12, mt13, mt21, mt22, mt23, mt31, mt32, mt33;  //Z=9231
    float g11, g12, g13, g21, g22, g23, g31, g32, g33, gi11, gi12, gi13, gi21, gi22, gi23, gi31, gi32, gi33, detg;  //Z=9232
    float g11t, g12t, g13t, g21t, g22t, g23t, g31t, g32t, g33t, gi11t, gi12t, gi13t, gi21t, gi22t, gi23t, gi31t, gi32t, gi33t, detgt;  //Z=9233
    //float mi11, mi12, mi13, mi21, mi22, mi23, mi31, mi32, mi33, detm;  //Z=9234
    double ri11, ri12, ri13, ri21, ri22, ri23, ri31, ri32, ri33; //, detr;  //Z=9235
    float aax, aay, aaz, bax, bay, baz, cax, cay, caz, nnhkl, unhkl, vnhkl;  //Z=9236
    float aaxt, aayt, aazt, baxt, bayt, bazt, caxt, cayt, cazt;  //Z=9237
    //float aex, aey, aez, bex, bey, bez, cex, cey, cez, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;  //Z=9238
    float aexn, aeyn, aezn, bexn, beyn, bezn, cexn, ceyn, cezn, asxn, asyn, aszn, bsxn, bsyn, bszn, csxn, csyn, cszn;  //Z=9239
    float aext, aeyt, aezt, bext, beyt, bezt, cext, ceyt, cezt, asxt, asyt, aszt, bsxt, bsyt, bszt, csxt, csyt, cszt;  //Z=9240
    //float asxd, asyd, aszd, bsxd, bsyd, bszd, csxd, csyd, cszd;  //Z=9241
    float epsi, etheta, len, r11, r12, r13, r21, r22, r23, r31, r32, r33;  //Z=9242
    //float aedx, aedy, aedz, bedx, bedy, bedz, cedx, cedy, cedz, vols;  //Z=9243
    float nxyzx, nxyzy, nxyzz, uxyzx, uxyzy, uxyzz, vxyzx, vxyzy, vxyzz, nxyznx, /*nxyzny,*/ nxyznz;  //Z=9244
    //double asax, asay, asaz, bsax, bsay, bsaz, csax, csay, csaz, asex, asey, asez, bsex, bsey, bsez, csex, csey, csez;  //Z=9245
    //double asedx, asedy, asedz, bsedx, bsedy, bsedz, csedx, csedy, csedz;  //Z=9246
    double asedxn, asedyn, asedzn, bsedxn, bsedyn, bsedzn, csedxn, csedyn, csedzn;  //Z=9247
    double asedxt, asedyt, asedzt, bsedxt, bsedyt, bsedzt, csedxt, csedyt, csedzt;  //Z=9248
    //double aedxt, aedyt, aedzt, bedxt, bedyt, bedzt, cedxt, cedyt, cedzt;  //Z=9249
    double m11DD, m12DD, m13DD, m21DD, m22DD, m23DD; //, m31DD, m32DD, m33DD;  //Z=9250
    //double mi11DD, mi12DD, mi13DD, mi21DD, mi22DD, mi23DD, mi31DD, mi32DD, mi33DD;  //Z=9251
    double aexDD, aeyDD, /*aezDD,*/ bexDD, beyDD; //, bezDD, cexDD, ceyDD, cezDD;  //Z=9252
    double asxDD, asyDD, aszDD, bsxDD, bsyDD, bszDD, csxDD, csyDD, cszDD, area;  //Z=9253
    //double aedxDD, aedyDD, aedzDD, bedxDD, bedyDD, bedzDD, cedxDD, cedyDD, cedzDD;  //Z=9254
    double asedxDD, asedyDD, asedzDD, bsedxDD, bsedyDD, bsedzDD, csedxDD, csedyDD, csedzDD;  //Z=9255
    //double m11DL, m12DL, m13DL, m21DL, m22DL, m23DL, m31DL, m32DL, m33DL;  //Z=9256
    //double mi11DL, mi12DL, mi13DL, mi21DL, mi22DL, mi23DL, mi31DL, mi32DL, mi33DL;  //Z=9257
    //double aexDL, aeyDL, aezDL, bexDL, beyDL, bezDL, cexDL, ceyDL, cezDL;  //Z=9258
    double asxDL, asyDL, aszDL, bsxDL, bsyDL, bszDL, csxDL, csyDL, cszDL, area1;  //Z=9259
    //double aedxDL, aedyDL, aedzDL, bedxDL, bedyDL, bedzDL, cedxDL, cedyDL, cedzDL;  //Z=9260
    double asedxDL, asedyDL, asedzDL, bsedxDL, bsedyDL, bsedzDL, csedxDL, csedyDL, csedzDL;  //Z=9261
    double thetarot, ct, st, ux, uy, uz, ul;  //Z=9262
    double rott[3+1][3+1];  //Z=9263

    //qDebug() << "corotations: twin =" << twin;

    /*begin*/  //Z=9265
    /* ***************** */  //Z=9266
    /* ** unit cell **** */  //Z=9267
    alpha = alpha*M_PI/180.0;  //Z=9268
    beta = beta*M_PI/180.0;  //Z=9269
    gamma = gamma*M_PI/180.0;  //Z=9270
    ephi = ephi*M_PI/180.0;  //Z=9271
    ca = cos(alpha);  //Z=9272
    cb = cos(beta);  //Z=9273
    cg = cos(gamma);  //Z=9274
    //sa = sin(alpha);  //Z=9275
    sb = sin(beta);  //Z=9276
    sg = sin(gamma);  //Z=9277
    //vol = a*b*c*sqrt(1.0-ca*ca-cb*cb-cg*cg+2*ca*cb*cg);  //Z=9278

    n1 = 0;  //Z=9280
    n2 = 0;  //Z=9281
    n3 = 1;  //Z=9282
    l1 = sb;  //Z=9283
    l2 = 0;  //Z=9284
    l3 = cb;  //Z=9285
    m1 = (cg-cb*ca)/sb;  //Z=9286
    m2 = sqrt(1.0-m1*m1-ca*ca);  //Z=9287
    m3 = ca;  //Z=9288

    /* ***************** */  //Z=9290
    /* ** Msi matrix *** */  //Z=9291
    msi11 = a*l1;         msi12 = a*l2;       msi13 = a*l3;  //Z=9292
    msi21 = b*m1;         msi22 = b*m2;       msi23 = b*m3;  //Z=9293
    msi31 = c*n1;         msi32 = c*n2;       msi33 = c*n3;  //Z=9294

    /*  new: transposed Msi matrix  */  //Z=9296
    /* msi11n:=a*l1;         msi21n:=a*l2;       msi31n:=a*l3;  //Z=9297
       msi12n:=b*m1;         msi22n:=b*m2;       msi32n:=b*m3;  //Z=9298
       msi13n:=c*n1;         msi23n:=c*n2;       msi33n:=c*n3;  */  //Z=9299

    /* ********************* */  //Z=9301
    /* ** new: = M matrix ** */  //Z=9302
    m11n = msi11;          m12n = msi12;       m13n = msi13;  //Z=9303
    m21n = msi21;          m22n = msi22;       m32n = msi23;  //Z=9304
    m31n = msi31;          m23n = msi32;       m33n = msi33;  //Z=9305

    m11DD = a;     m12DD = b*cg;    m13DD = 0;  //Z=9307
    m21DD = 0;     m22DD = b*sg;    m23DD = 0;  //Z=9308
    //m31DD = 0;     m32DD = 0;       m33DD = 1;  //Z=9309

    //m11DL = a;     m12DL = 0;    m13DL = 0;  //Z=9311
    //m21DL = 0;     m22DL = 1;    m23DL = 0;  //Z=9312
    //m31DL = 0;     m32DL = 0;    m33DL = 1;  //Z=9313

    /*  old: Ms matrix  */  //Z=9315
    /* detmsi:=msi11*msi22*msi33+msi12*msi23*msi31+msi13*msi21*msi32-msi31*msi22*msi13-msi32*msi23*msi11-msi33*msi21*msi12;  //Z=9316
       ms11:=(msi22*msi33-msi23*msi32)/detmsi;  //Z=9317
       ms12:=(msi13*msi32-msi12*msi33)/detmsi;  //Z=9318
       ms13:=(msi12*msi23-msi13*msi22)/detmsi;  //Z=9319
       ms21:=(msi23*msi31-msi21*msi33)/detmsi;  //Z=9320
       ms22:=(msi11*msi33-msi13*msi31)/detmsi;  //Z=9321
       ms23:=(msi13*msi21-msi11*msi23)/detmsi;  //Z=9322
       ms31:=(msi21*msi32-msi22*msi31)/detmsi;  //Z=9323
       ms32:=(msi12*msi31-msi11*msi32)/detmsi;  //Z=9324
       ms33:=(msi11*msi22-msi12*msi21)/detmsi;    */  //Z=9325

    /*  old: Ms matrix, transposed  */  //Z=9327
    /* mst11:=ms11;    mst12:=ms21;    mst13:=ms31;  //Z=9328
       mst21:=ms12;    mst22:=ms22;    mst23:=ms32;  //Z=9329
       mst31:=ms13;    mst32:=ms23;    mst33:=ms33;   */  //Z=9330

    /*  old: M matrix, transposed  */  //Z=9332
    /* detmst:=mst11*mst22*mst33+mst12*mst23*mst31+mst13*mst21*mst32-mst31*mst22*mst13-mst32*mst23*mst11-mst33*mst21*mst12;  //Z=9333
       mt11:=(mst22*mst33-mst23*mst32)/detmst;  //Z=9334
       mt12:=(mst13*mst32-mst12*mst33)/detmst;  //Z=9335
       mt13:=(mst12*mst23-mst13*mst22)/detmst;  //Z=9336
       mt21:=(mst23*mst31-mst21*mst33)/detmst;  //Z=9337
       mt22:=(mst11*mst33-mst13*mst31)/detmst;  //Z=9338
       mt23:=(mst13*mst21-mst11*mst23)/detmst;  //Z=9339
       mt31:=(mst21*mst32-mst22*mst31)/detmst;  //Z=9340
       mt32:=(mst12*mst31-mst11*mst32)/detmst;  //Z=9341
       mt33:=(mst11*mst22-mst12*mst21)/detmst;   */  //Z=9342

    /*  old: M matrix  */  //Z=9344
    /*  rxyz=M ruvw  */  //Z=9345
    /* m11:=mt11;       m12:=mt21;        m13:=mt31;  //Z=9346
       m21:=mt12;       m22:=mt22;        m23:=mt32;  //Z=9347
       m31:=mt13;       m32:=mt23;        m33:=mt33;      */  //Z=9348


    /*  Mi inverse matrix  */  //Z=9351
    /*  ruvw=Mi rxyz  */  //Z=9352
    /* detm:=m11*m22*m33+m12*m23*m31+m13*m21*m32-m31*m22*m13-m32*m23*m11-m33*m21*m12;  //Z=9353
       mi11:=(m22*m33-m23*m32)/detm;  //Z=9354
       mi12:=(m13*m32-m12*m33)/detm;  //Z=9355
       mi13:=(m12*m23-m13*m22)/detm;  //Z=9356
       mi21:=(m23*m31-m21*m33)/detm;  //Z=9357
       mi22:=(m11*m33-m13*m31)/detm;  //Z=9358
       mi23:=(m13*m21-m11*m23)/detm;  //Z=9359
       mi31:=(m21*m32-m22*m31)/detm;  //Z=9360
       mi32:=(m12*m31-m11*m32)/detm;  //Z=9361
       mi33:=(m11*m22-m12*m21)/detm;  //Z=9362

       mi11DL:=1/a;      mi12DL:=0;        mi13DL:=0;  //Z=9364
       mi21DL:=0;        mi22DL:=1;        mi23DL:=0;  //Z=9365
       mi31DL:=0;        mi32DL:=0;        mi33DL:=1;   */  //Z=9366

    /*  new Mi inverse matrix  */  //Z=9368
    //detm = m11n*m22n*m33n+m12n*m23n*m31n+m13n*m21n*m32n-m31n*m22n*m13n-m32n*m23n*m11n-m33n*m21n*m12n;  //Z=9369
    //mi11 = (m22n*m33n-m23n*m32n)/detm;  //Z=9370
    //mi12 = (m13n*m32n-m12n*m33n)/detm;  //Z=9371
    //mi13 = (m12n*m23n-m13n*m22n)/detm;  //Z=9372
    //mi21 = (m23n*m31n-m21n*m33n)/detm;  //Z=9373
    //mi22 = (m11n*m33n-m13n*m31n)/detm;  //Z=9374
    //mi23 = (m13n*m21n-m11n*m23n)/detm;  //Z=9375
    //mi31 = (m21n*m32n-m22n*m31n)/detm;  //Z=9376
    //mi32 = (m12n*m31n-m11n*m32n)/detm;  //Z=9377
    //mi33 = (m11n*m22n-m12n*m21n)/detm;  //Z=9378

    /* *************************** */  //Z=9380
    /*  base vectors of unit cell  */  //Z=9381
    /*  aa, ba, ca in uvw-system   */  //Z=9382
    aax = 1;         aay = 0;      aaz = 0;  //Z=9383
    bax = 0;         bay = 1;      baz = 0;  //Z=9384
    cax = 0;         cay = 0;      caz = 1;  //Z=9385


    /*  generate fcc- and bcc-twin from rotation matrix  */  //Z=9388
    if ( twin )
    {/*2*/  //Z=9389
        /*  bcc  */  //Z=9390
        if ( ltype==4 )
        {/*3*/  //Z=9391
            ux = 1;  //Z=9392
            uy = 1;  //Z=9393
            uz = 1;  //Z=9394
            thetarot = 60;  //Z=9395
        }/*3*/  //Z=9396
        /*  fcc  */  //Z=9397
        if ( ltype==5 )
        {/*3*/  //Z=9398
            ux = -1;  //Z=9399
            uy = 1;  //Z=9400
            uz = 1;  //Z=9401
            thetarot = 180;  //Z=9402
        }/*3*/  //Z=9403
        ul = sqrt(ux*ux+uy*uy+uz*uz);  //Z=9404
        ux = ux/ul;  //Z=9405
        uy = uy/ul;  //Z=9406
        uz = uz/ul;  //Z=9407
        ct = cos(thetarot*M_PI/180.0);  //Z=9408
        st = sin(thetarot*M_PI/180.0);  //Z=9409
        rott[1][1] = ct+ux*ux*(1-ct);  //Z=9410
        rott[1][2] = ux*uy*(1-ct)-uz*st;  //Z=9411
        rott[1][3] = ux*uz*(1-ct)+uy*st;  //Z=9412
        rott[2][1] = uy*ux*(1-ct)+uz*st;  //Z=9413
        rott[2][2] = ct+uy*uy*(1-ct);  //Z=9414
        rott[2][3] = uy*uz*(1-ct)-ux*st;  //Z=9415
        rott[3][1] = uz*ux*(1-ct)-uy*st;  //Z=9416
        rott[3][2] = uz*uy*(1-ct)+ux*st;  //Z=9417
        rott[3][3] = ct+uz*uz*(1-ct);  //Z=9418
        aaxt = rott[1][1]*aax+rott[1][2]*aay+rott[1][3]*aaz;  //Z=9419
        aayt = rott[2][1]*aax+rott[2][2]*aay+rott[2][3]*aaz;  //Z=9420
        aazt = rott[3][1]*aax+rott[3][2]*aay+rott[3][3]*aaz;  //Z=9421
        baxt = rott[1][1]*bax+rott[1][2]*bay+rott[1][3]*baz;  //Z=9422
        bayt = rott[2][1]*bax+rott[2][2]*bay+rott[2][3]*baz;  //Z=9423
        bazt = rott[3][1]*bax+rott[3][2]*bay+rott[3][3]*baz;  //Z=9424
        caxt = rott[1][1]*cax+rott[1][2]*cay+rott[1][3]*caz;  //Z=9425
        cayt = rott[2][1]*cax+rott[2][2]*cay+rott[2][3]*caz;  //Z=9426
        cazt = rott[3][1]*cax+rott[3][2]*cay+rott[3][3]*caz;  //Z=9427
    }/*2*/  //Z=9428



    /*  unit vectors in carthesian coordinate system  */  //Z=9432
    /*  ae=M aa, be=M ba, ce=M ca in xyz-system  */  //Z=9433
    /* aex:=m11*aax+m12*aay+m13*aaz;  //Z=9434
       aey:=m21*aax+m22*aay+m23*aaz;  //Z=9435
       aez:=m31*aax+m32*aay+m33*aaz;  //Z=9436
       bex:=m11*bax+m12*bay+m13*baz;  //Z=9437
       bey:=m21*bax+m22*bay+m23*baz;  //Z=9438
       bez:=m31*bax+m32*bay+m33*baz;  //Z=9439
       cex:=m11*cax+m12*cay+m13*caz;  //Z=9440
       cey:=m21*cax+m22*cay+m23*caz;  //Z=9441
       cez:=m31*cax+m32*cay+m33*caz;   */  //Z=9442

    /* ********************************************** */  //Z=9444
    /*  new unit vectors in carthesian coordinate system  */  //Z=9445
    aexn = m11n*aax+m12n*aay+m13n*aaz;  //Z=9446
    aeyn = m21n*aax+m22n*aay+m23n*aaz;  //Z=9447
    aezn = m31n*aax+m32n*aay+m33n*aaz;  //Z=9448
    bexn = m11n*bax+m12n*bay+m13n*baz;  //Z=9449
    beyn = m21n*bax+m22n*bay+m23n*baz;  //Z=9450
    bezn = m31n*bax+m32n*bay+m33n*baz;  //Z=9451
    cexn = m11n*cax+m12n*cay+m13n*caz;  //Z=9452
    ceyn = m21n*cax+m22n*cay+m23n*caz;  //Z=9453
    cezn = m31n*cax+m32n*cay+m33n*caz;  //Z=9454

    if ( twin )
    {/*2*/  //Z=9456
        aext = m11n*aaxt+m12n*aayt+m13n*aazt;  //Z=9457
        aeyt = m21n*aaxt+m22n*aayt+m23n*aazt;  //Z=9458
        aezt = m31n*aaxt+m32n*aayt+m33n*aazt;  //Z=9459
        bext = m11n*baxt+m12n*bayt+m13n*bazt;  //Z=9460
        beyt = m21n*baxt+m22n*bayt+m23n*bazt;  //Z=9461
        bezt = m31n*baxt+m32n*bayt+m33n*bazt;  //Z=9462
        cext = m11n*caxt+m12n*cayt+m13n*cazt;  //Z=9463
        ceyt = m21n*caxt+m22n*cayt+m23n*cazt;  //Z=9464
        cezt = m31n*caxt+m32n*cayt+m33n*cazt;  //Z=9465
    }/*2*/  //Z=9466

    aexDD = m11DD*aax+m12DD*aay+m13DD*aaz;  //Z=9468
    aeyDD = m21DD*aax+m22DD*aay+m23DD*aaz;  //Z=9469
    //aezDD = m31DD*aax+m32DD*aay+m33DD*aaz;  //Z=9470
    bexDD = m11DD*bax+m12DD*bay+m13DD*baz;  //Z=9471
    beyDD = m21DD*bax+m22DD*bay+m23DD*baz;  //Z=9472
    //bezDD = m31DD*bax+m32DD*bay+m33DD*baz;  //Z=9473
    //cexDD = m11DD*cax+m12DD*cay+m13DD*caz;  //Z=9474
    //ceyDD = m21DD*cax+m22DD*cay+m23DD*caz;  //Z=9475
    //cezDD = m31DD*cax+m32DD*cay+m33DD*caz;  //Z=9476

    //aexDL = m11DL*aax+m12DL*aay+m13DL*aaz;  //Z=9478
    //aeyDL = m21DL*aax+m22DL*aay+m23DL*aaz;  //Z=9479
    //aezDL = m31DL*aax+m32DL*aay+m33DL*aaz;  //Z=9480
    //bexDL = m11DL*bax+m12DL*bay+m13DL*baz;  //Z=9481
    //beyDL = m21DL*bax+m22DL*bay+m23DL*baz;  //Z=9482
    //bezDL = m31DL*bax+m32DL*bay+m33DL*baz;  //Z=9483
    //cexDL = m11DL*cax+m12DL*cay+m13DL*caz;  //Z=9484
    //ceyDL = m21DL*cax+m22DL*cay+m23DL*caz;  //Z=9485
    //cezDL = m31DL*cax+m32DL*cay+m33DL*caz;  //Z=9486

    /*  generate fcc-twin from rotation matrix  */  //Z=9488
    /* ux:=1;  //Z=9489
       uy:=1;  //Z=9490
       uz:=1;  //Z=9491
       ul:=sqrt(ux*ux+uy*uy+uz*uz);  //Z=9492
       ux:=ux/ul;  //Z=9493
       uy:=uy/ul;  //Z=9494
       uz:=uz/ul;  //Z=9495
       thetarot:=60;  //Z=9496
       ct:=cos(thetarot*pi/180);  //Z=9497
       st:=sin(thetarot*pi/180);  //Z=9498
       rott[1,1]:=ct+ux*ux*(1-ct);  //Z=9499
       rott[1,2]:=ux*uy*(1-ct)-uz*st;  //Z=9500
       rott[1,3]:=ux*uz*(1-ct)+uy*st;  //Z=9501
       rott[2,1]:=uy*ux*(1-ct)+uz*st;  //Z=9502
       rott[2,2]:=ct+uy*uy*(1-ct);  //Z=9503
       rott[2,3]:=uy*uz*(1-ct)-ux*st;  //Z=9504
       rott[3,1]:=uz*ux*(1-ct)-uy*st;  //Z=9505
       rott[3,2]:=uz*uy*(1-ct)+ux*st;  //Z=9506
       rott[3,3]:=ct+uz*uz*(1-ct);  //Z=9507
       aext:=rott[1,1]*aexn+rott[1,2]*aeyn+rott[1,3]*aezn;  //Z=9508
       aeyt:=rott[2,1]*aexn+rott[2,2]*aeyn+rott[2,3]*aezn;  //Z=9509
       aezt:=rott[3,1]*aexn+rott[3,2]*aeyn+rott[3,3]*aezn;  //Z=9510
       bext:=rott[1,1]*bexn+rott[1,2]*beyn+rott[1,3]*bezn;  //Z=9511
       beyt:=rott[2,1]*bexn+rott[2,2]*beyn+rott[2,3]*bezn;  //Z=9512
       bezt:=rott[3,1]*bexn+rott[3,2]*beyn+rott[3,3]*bezn;  //Z=9513
       cext:=rott[1,1]*cexn+rott[1,2]*ceyn+rott[1,3]*cezn;  //Z=9514
       ceyt:=rott[2,1]*cexn+rott[2,2]*ceyn+rott[2,3]*cezn;  //Z=9515
       cezt:=rott[3,1]*cexn+rott[3,2]*ceyn+rott[3,3]*cezn;  */  //Z=9516


    /*  old: reciprocal space vector in carthesian coordinate system  */  //Z=9519
    /*  ase, bse, cse in xyz-system  */  //Z=9520
    /*  Mok  */  //Z=9521
    /* vvol:=aex*(bey*cez-bez*cey)+aey*(bez*cex-bex*cez)+aez*(bex*cey-bey*cex);  //Z=9522
       asx:=(bey*cez-bez*cey)/vvol;  //Z=9523
       asy:=(bez*cex-bex*cez)/vvol;  //Z=9524
       asz:=(bex*cey-bey*cex)/vvol;  //Z=9525
       bsx:=(aez*cey-aey*cez)/vvol;  //Z=9526
       bsy:=(aex*cez-aez*cex)/vvol;  //Z=9527
       bsz:=(aey*cex-aex*cey)/vvol;  //Z=9528
       csx:=(aey*bez-aez*bey)/vvol;  //Z=9529
       csy:=(aez*bex-aex*bez)/vvol;  //Z=9530
       csz:=(aex*bey-aey*bex)/vvol;  //Z=9531

       asxt:=(beyt*cezt-bezt*ceyt)/vvol;  //Z=9533
       asyt:=(bezt*cext-bext*cezt)/vvol;  //Z=9534
       aszt:=(bext*ceyt-beyt*cext)/vvol;  //Z=9535
       bsxt:=(aezt*ceyt-aeyt*cezt)/vvol;  //Z=9536
       bsyt:=(aext*cezt-aezt*cext)/vvol;  //Z=9537
       bszt:=(aeyt*cext-aext*ceyt)/vvol;  //Z=9538
       csxt:=(aeyt*bezt-aezt*beyt)/vvol;  //Z=9539
       csyt:=(aezt*bext-aext*bezt)/vvol;  //Z=9540
       cszt:=(aext*beyt-aeyt*bext)/vvol;   */  //Z=9541

    /* ***************************************** */  //Z=9543
    /*  new: G metric matrix in xyz-coordinates  */  //Z=9544
    g11 = aexn;     g12 = bexn;     g13 = cexn;  //Z=9545
    g21 = aeyn;     g22 = beyn;     g23 = ceyn;  //Z=9546
    g31 = aezn;     g32 = bezn;     g33 = cezn;  //Z=9547

    g11t = aext;     g12t = bext;     g13t = cext;  //Z=9549
    g21t = aeyt;     g22t = beyt;     g23t = ceyt;  //Z=9550
    g31t = aezt;     g32t = bezt;     g33t = cezt;  //Z=9551

    /*  old G matrix  */  //Z=9553
    /* g11:=aex*aex+aey*aey+aez*aez;  //Z=9554
       g12:=aex*bex+aey*bey+aez*bez;  //Z=9555
       g13:=aex*cex+aey*cey+aez*cez;  //Z=9556
       g21:=bex*aex+bey*aey+bez*aez;  //Z=9557
       g22:=bex*bex+bey*bey+bez*bez;  //Z=9558
       g23:=bex*cex+bey*cey+bez*cez;  //Z=9559
       g31:=cex*aex+cey*aey+cez*aez;  //Z=9560
       g32:=cex*bex+cey*bey+cez*bez;  //Z=9561
       g33:=cex*cex+cey*cey+cez*cez;  */  //Z=9562

    /*  Gs inverse metric matrix  */  //Z=9564
    detg = g11*g22*g33+g12*g23*g31+g13*g21*g32-g31*g22*g13-g32*g23*g11-g33*g21*g12;  //Z=9565
    gi11 = (g22*g33-g23*g32)/detg;  //Z=9566
    gi12 = (g13*g32-g12*g33)/detg;  //Z=9567
    gi13 = (g12*g23-g13*g22)/detg;  //Z=9568
    gi21 = (g23*g31-g21*g33)/detg;  //Z=9569
    gi22 = (g11*g33-g13*g31)/detg;  //Z=9570
    gi23 = (g13*g21-g11*g23)/detg;  //Z=9571
    gi31 = (g21*g32-g22*g31)/detg;  //Z=9572
    gi32 = (g12*g31-g11*g32)/detg;  //Z=9573
    gi33 = (g11*g22-g12*g21)/detg;  //Z=9574

    if ( twin )
    {/*2*/  //Z=9576
        detgt = g11t*g22t*g33t+g12t*g23t*g31t+g13t*g21t*g32t-g31t*g22t*g13t-g32t*g23t*g11t-g33t*g21t*g12t;  //Z=9577
        gi11t = (g22t*g33t-g23t*g32t)/detgt;  //Z=9578
        gi12t = (g13t*g32t-g12t*g33t)/detgt;  //Z=9579
        gi13t = (g12t*g23t-g13t*g22t)/detgt;  //Z=9580
        gi21t = (g23t*g31t-g21t*g33t)/detgt;  //Z=9581
        gi22t = (g11t*g33t-g13t*g31t)/detgt;  //Z=9582
        gi23t = (g13t*g21t-g11t*g23t)/detgt;  //Z=9583
        gi31t = (g21t*g32t-g22t*g31t)/detgt;  //Z=9584
        gi32t = (g12t*g31t-g11t*g32t)/detgt;  //Z=9585
        gi33t = (g11t*g22t-g12t*g21t)/detgt;  //Z=9586
    }/*2*/  //Z=9587

    /*  new: reciprocal space vector in carthesian coordinate system  */  //Z=9589
    /*  ase, bse, cse in xyz-system  */  //Z=9590
    asxn = gi11;   asyn = gi12;   aszn = gi13;  //Z=9591
    bsxn = gi21;   bsyn = gi22;   bszn = gi23;  //Z=9592
    csxn = gi31;   csyn = gi32;   cszn = gi33;  //Z=9593

    if ( twin )
    {/*2*/  //Z=9595
        asxt = gi11t;   asyt = gi12t;   aszt = gi13t;  //Z=9596
        bsxt = gi21t;   bsyt = gi22t;   bszt = gi23t;  //Z=9597
        csxt = gi31t;   csyt = gi32t;   cszt = gi33t;  //Z=9598
    }/*2*/  //Z=9599

    area = a*b*sg;  //Z=9601
    asxDD = beyDD/area;    asyDD = -bexDD/area;     aszDD = 0;  //Z=9602
    bsxDD = -aeyDD/area;   bsyDD = aexDD/area;      bszDD = 0;  //Z=9603
    csxDD = 0;             csyDD = 0;               cszDD = 1;  //Z=9604

    area1 = a;  //Z=9606
    asxDL = 1/area1;       asyDL = 0;               aszDL = 0;  //Z=9607
    bsxDL = 0;             bsyDL = 1;               bszDL = 0;  //Z=9608
    csxDL = 0;             csyDL = 0;               cszDL = 1;  //Z=9609


    /*  a*,b*,c* reciprocal space vectors  */  //Z=9612
    /*  in uvw-space  */  //Z=9613
    /* asax:=g11*aax+g12*aay+g13*aaz;  //Z=9614
       asay:=g21*aax+g22*aay+g23*aaz;  //Z=9615
       asaz:=g31*aax+g32*aay+g33*aaz;  //Z=9616
       bsax:=g11*bax+g12*bay+g13*baz;  //Z=9617
       bsay:=g21*bax+g22*bay+g23*baz;  //Z=9618
       bsaz:=g31*bax+g32*bay+g33*baz;  //Z=9619
       csax:=g11*cax+g12*cay+g13*caz;  //Z=9620
       csay:=g21*cax+g22*cay+g23*caz;  //Z=9621
       csaz:=g31*cax+g32*cay+g33*caz; */  //Z=9622

    /*  a*, b*, c* reciprocal space vectors  */  //Z=9624
    /*  in xyz-space  */  //Z=9625
    /* asex:=m11*asax+m12*asay+m13*asaz;  //Z=9626
       asey:=m21*asax+m22*asay+m23*asaz;  //Z=9627
       asez:=m31*asax+m32*asay+m33*asaz;  //Z=9628
       bsex:=m11*bsax+m12*bsay+m13*bsaz;  //Z=9629
       bsey:=m21*bsax+m22*bsay+m23*bsaz;  //Z=9630
       bsez:=m31*bsax+m32*bsay+m33*bsaz;  //Z=9631
       csex:=m11*csax+m12*csay+m13*csaz;  //Z=9632
       csey:=m21*csax+m22*csay+m23*csaz;  //Z=9633
       csez:=m31*csax+m32*csay+m33*csaz; */  //Z=9634


    /*  nuvw-vector || beam  */  //Z=9637
    /*  nuvw in unit cell uvw-system  */  //Z=9638
    nuvwx = u;  //Z=9639
    nuvwy = v;  //Z=9640
    nuvwz = w;  //Z=9641
    lenn = sqrt(u*u+v*v+w*w);  //Z=9642

    /*  nhkl-vector  */  //Z=9644
    /*  nhkl=G nuvw in unit cell uvw-system  */  //Z=9645
    /* nhklx:=g11*nuvwx+g12*nuvwy+g13*nuvwz;  //Z=9646
       nhkly:=g21*nuvwx+g22*nuvwy+g23*nuvwz;  //Z=9647
       nhklz:=g31*nuvwx+g32*nuvwy+g33*nuvwz; */  //Z=9648

    /*  nhkl-vector  */  //Z=9650
    /*  nxyz=M nuvw in xyz-system  */  //Z=9651
    /*  with (a,b,c)=(2,1,1) and (u,v,w)=(2,1,1) this is a (ua,vb,wc)=(4,1,1)-vector  */  //Z=9652
    nxyzx = m11n*nuvwx+m12n*nuvwy+m13n*nuvwz;  //Z=9653
    nxyzy = m21n*nuvwx+m22n*nuvwy+m23n*nuvwz;  //Z=9654
    nxyzz = m31n*nuvwx+m32n*nuvwy+m33n*nuvwz;  //Z=9655

    /*  unit nxyz = director  */  //Z=9657
    /*  in xyz-system  */  //Z=9658
    len = sqrt(nxyzx*nxyzx+nxyzy*nxyzy+nxyzz*nxyzz);  //Z=9659
    nxyznx = nxyzx/len;  //Z=9660
    //nxyzny = nxyzy/len;  //Z=9661
    nxyznz = nxyzz/len;  //Z=9662

    /*  R rotation matrix  */  //Z=9664
    /*  in xyz-system  */  //Z=9665
    etheta = myacos(-nxyznz);  //Z=9666
    epsi = asin(nxyznx/sin(etheta));  //Z=9667

    r11 = cos(epsi)*cos(ephi)-sin(epsi)*cos(etheta)*sin(ephi);  //Z=9669
    r12 = -cos(epsi)*sin(ephi)-sin(epsi)*cos(etheta)*cos(ephi);  //Z=9670
    r13 = sin(epsi)*sin(etheta);  //Z=9671
    r21 = sin(epsi)*cos(ephi)+cos(epsi)*cos(etheta)*sin(ephi);  //Z=9672
    r22 = -sin(epsi)*sin(ephi)+cos(epsi)*cos(etheta)*cos(ephi);  //Z=9673
    r23 = -cos(epsi)*sin(etheta);  //Z=9674
    r31 = sin(etheta)*sin(ephi);  //Z=9675
    r32 = sin(etheta)*cos(ephi);  //Z=9676
    r33 = cos(etheta);  //Z=9677

    /*  Ri inverse rotation matrix  */  //Z=9679
    /*  in xyz-system  */  //Z=9680
    /* detr:=r11*r22*r33+r12*r23*r31+r13*r21*r32-r31*r22*r13-r32*r23*r11-r33*r21*r12;  //Z=9681
       ri11:=(r22*r33-r23*r32)/detr;  //Z=9682
       ri12:=(r13*r32-r12*r33)/detr;  //Z=9683
       ri13:=(r12*r23-r13*r22)/detr;  //Z=9684
       ri21:=(r23*r31-r21*r33)/detr;  //Z=9685
       ri22:=(r11*r33-r13*r31)/detr;  //Z=9686
       ri23:=(r13*r21-r11*r23)/detr;  //Z=9687
       ri31:=(r21*r32-r22*r31)/detr;  //Z=9688
       ri32:=(r12*r31-r11*r32)/detr;  //Z=9689
       ri33:=(r11*r22-r12*r21)/detr;     */  //Z=9690

    ri11 = r11;  //Z=9692
    ri12 = r21;  //Z=9693
    ri13 = r31;  //Z=9694
    ri21 = r12;  //Z=9695
    ri22 = r22;  //Z=9696
    ri23 = r32;  //Z=9697
    ri31 = r13;  //Z=9698
    ri32 = r23;  //Z=9699
    ri33 = r33;  //Z=9700

    /*  rotated base vectors a,b,c in carthesian coordinate system  */  //Z=9702
    /*  aed=Ri ae, bed=Ri be, ced=Ri ce in xyz-system  */  //Z=9703
    /*  needed for calculation of fiber pattern  */  //Z=9704
    /*  Mok  */  //Z=9705
    /* aedx:=ri11*aex+ri12*aey+ri13*aez;  //Z=9706
       aedy:=ri21*aex+ri22*aey+ri23*aez;  //Z=9707
       aedz:=ri31*aex+ri32*aey+ri33*aez;  //Z=9708
       bedx:=ri11*bex+ri12*bey+ri13*bez;  //Z=9709
       bedy:=ri21*bex+ri22*bey+ri23*bez;  //Z=9710
       bedz:=ri31*bex+ri32*bey+ri33*bez;  //Z=9711
       cedx:=ri11*cex+ri12*cey+ri13*cez;  //Z=9712
       cedy:=ri21*cex+ri22*cey+ri23*cez;  //Z=9713
       cedz:=ri31*cex+ri32*cey+ri33*cez;  //Z=9714

       aedxDD:=ri11*aexDD+ri12*aeyDD+ri13*aezDD;  //Z=9716
       aedyDD:=ri21*aexDD+ri22*aeyDD+ri23*aezDD;  //Z=9717
       aedzDD:=ri31*aexDD+ri32*aeyDD+ri33*aezDD;  //Z=9718
       bedxDD:=ri11*bexDD+ri12*beyDD+ri13*bezDD;  //Z=9719
       bedyDD:=ri21*bexDD+ri22*beyDD+ri23*bezDD;  //Z=9720
       bedzDD:=ri31*bexDD+ri32*beyDD+ri33*bezDD;  //Z=9721
       cedxDD:=ri11*cexDD+ri12*ceyDD+ri13*cezDD;  //Z=9722
       cedyDD:=ri21*cexDD+ri22*ceyDD+ri23*cezDD;  //Z=9723
       cedzDD:=ri31*cexDD+ri32*ceyDD+ri33*cezDD;  //Z=9724

       aedxDL:=ri11*aexDL+ri12*aeyDL+ri13*aezDL;  //Z=9726
       aedyDL:=ri21*aexDL+ri22*aeyDL+ri23*aezDL;  //Z=9727
       aedzDL:=ri31*aexDL+ri32*aeyDL+ri33*aezDL;  //Z=9728
       bedxDL:=ri11*bexDL+ri12*beyDL+ri13*bezDL;  //Z=9729
       bedyDL:=ri21*bexDL+ri22*beyDL+ri23*bezDL;  //Z=9730
       bedzDL:=ri31*bexDL+ri32*beyDL+ri33*bezDL;  //Z=9731
       cedxDL:=ri11*cexDL+ri12*ceyDL+ri13*cezDL;  //Z=9732
       cedyDL:=ri21*cexDL+ri22*ceyDL+ri23*cezDL;  //Z=9733
       cedzDL:=ri31*cexDL+ri32*ceyDL+ri33*cezDL;  //Z=9734

       aedxt:=ri11*aext+ri12*aeyt+ri13*aezt;  //Z=9736
       aedyt:=ri21*aext+ri22*aeyt+ri23*aezt;  //Z=9737
       aedzt:=ri31*aext+ri32*aeyt+ri33*aezt;  //Z=9738
       bedxt:=ri11*bext+ri12*beyt+ri13*bezt;  //Z=9739
       bedyt:=ri21*bext+ri22*beyt+ri23*bezt;  //Z=9740
       bedzt:=ri31*bext+ri32*beyt+ri33*bezt;  //Z=9741
       cedxt:=ri11*cext+ri12*ceyt+ri13*cezt;  //Z=9742
       cedyt:=ri21*cext+ri22*ceyt+ri23*cezt;  //Z=9743
       cedzt:=ri31*cext+ri32*ceyt+ri33*cezt;      */  //Z=9744

    /*  rotated reciprocal space vectors a*,b*,c* in carthesian coordinate system  */  //Z=9746
    /*  calculated from cross-product of aed,bed,ced  */  //Z=9747
    /*  into output matrix  */  //Z=9748
    /*  Mok  */  //Z=9749
    /* asedx:=(bedy*cedz-bedz*cedy)/vvol;  //Z=9750
       asedy:=(bedz*cedx-bedx*cedz)/vvol;  //Z=9751
       asedz:=(bedx*cedy-bedy*cedx)/vvol;  //Z=9752
       bsedx:=(aedz*cedy-aedy*cedz)/vvol;  //Z=9753
       bsedy:=(aedx*cedz-aedz*cedx)/vvol;  //Z=9754
       bsedz:=(aedy*cedx-aedx*cedy)/vvol;  //Z=9755
       csedx:=(aedy*bedz-aedz*bedy)/vvol;  //Z=9756
       csedy:=(aedz*bedx-aedx*bedz)/vvol;  //Z=9757
       csedz:=(aedx*bedy-aedy*bedx)/vvol;  */  //Z=9758

    /*  new: a*_r=R a*  */  //Z=9760
    /*  output  */  //Z=9761
    asedxn = ri11*asxn+ri12*asyn+ri13*aszn;  //Z=9762
    asedyn = ri21*asxn+ri22*asyn+ri23*aszn;  //Z=9763
    asedzn = ri31*asxn+ri32*asyn+ri33*aszn;  //Z=9764
    bsedxn = ri11*bsxn+ri12*bsyn+ri13*bszn;  //Z=9765
    bsedyn = ri21*bsxn+ri22*bsyn+ri23*bszn;  //Z=9766
    bsedzn = ri31*bsxn+ri32*bsyn+ri33*bszn;  //Z=9767
    csedxn = ri11*csxn+ri12*csyn+ri13*cszn;  //Z=9768
    csedyn = ri21*csxn+ri22*csyn+ri23*cszn;  //Z=9769
    csedzn = ri31*csxn+ri32*csyn+ri33*cszn;  //Z=9770

    if ( twin )
    {/*2*/  //Z=9772
        asedxt = ri11*asxt+ri12*asyt+ri13*aszt;  //Z=9773
        asedyt = ri21*asxt+ri22*asyt+ri23*aszt;  //Z=9774
        asedzt = ri31*asxt+ri32*asyt+ri33*aszt;  //Z=9775
        bsedxt = ri11*bsxt+ri12*bsyt+ri13*bszt;  //Z=9776
        bsedyt = ri21*bsxt+ri22*bsyt+ri23*bszt;  //Z=9777
        bsedzt = ri31*bsxt+ri32*bsyt+ri33*bszt;  //Z=9778
        csedxt = ri11*csxt+ri12*csyt+ri13*cszt;  //Z=9779
        csedyt = ri21*csxt+ri22*csyt+ri23*cszt;  //Z=9780
        csedzt = ri31*csxt+ri32*csyt+ri33*cszt;  //Z=9781
    }/*2*/  //Z=9782

    asedxDD = ri11*asxDD+ri12*asyDD+ri13*aszDD;  //Z=9784
    asedyDD = ri21*asxDD+ri22*asyDD+ri23*aszDD;  //Z=9785
    asedzDD = ri31*asxDD+ri32*asyDD+ri33*aszDD;  //Z=9786
    bsedxDD = ri11*bsxDD+ri12*bsyDD+ri13*bszDD;  //Z=9787
    bsedyDD = ri21*bsxDD+ri22*bsyDD+ri23*bszDD;  //Z=9788
    bsedzDD = ri31*bsxDD+ri32*bsyDD+ri33*bszDD;  //Z=9789
    csedxDD = ri11*csxDD+ri12*csyDD+ri13*cszDD;  //Z=9790
    csedyDD = ri21*csxDD+ri22*csyDD+ri23*cszDD;  //Z=9791
    csedzDD = ri31*csxDD+ri32*csyDD+ri33*cszDD;  //Z=9792

    asedxDL = ri11*asxDL+ri12*asyDL+ri13*aszDL;  //Z=9794
    asedyDL = ri21*asxDL+ri22*asyDL+ri23*aszDL;  //Z=9795
    asedzDL = ri31*asxDL+ri32*asyDL+ri33*aszDL;  //Z=9796
    bsedxDL = ri11*bsxDL+ri12*bsyDL+ri13*bszDL;  //Z=9797
    bsedyDL = ri21*bsxDL+ri22*bsyDL+ri23*bszDL;  //Z=9798
    bsedzDL = ri31*bsxDL+ri32*bsyDL+ri33*bszDL;  //Z=9799
    csedxDL = ri11*csxDL+ri12*csyDL+ri13*cszDL;  //Z=9800
    csedyDL = ri21*csxDL+ri22*csyDL+ri23*cszDL;  //Z=9801
    csedzDL = ri31*csxDL+ri32*csyDL+ri33*cszDL;  //Z=9802

    /* asedxt:=(bedyt*cedzt-bedzt*cedyt)/vvol;  //Z=9804
       asedyt:=(bedzt*cedxt-bedxt*cedzt)/vvol;  //Z=9805
       asedzt:=(bedxt*cedyt-bedyt*cedxt)/vvol;  //Z=9806
       bsedxt:=(aedzt*cedyt-aedyt*cedzt)/vvol;  //Z=9807
       bsedyt:=(aedxt*cedzt-aedzt*cedxt)/vvol;  //Z=9808
       bsedzt:=(aedyt*cedxt-aedxt*cedyt)/vvol;  //Z=9809
       csedxt:=(aedyt*bedzt-aedzt*bedyt)/vvol;  //Z=9810
       csedyt:=(aedzt*bedxt-aedxt*bedzt)/vvol;  //Z=9811
       csedzt:=(aedxt*bedyt-aedyt*bedxt)/vvol;  */  //Z=9812

    /*  a*, b*, c* rotated  */  //Z=9814
    /*  in xyz-space  */  //Z=9815
    /* asedx:=r11*asex+r12*asey+r13*asez;  //Z=9816
       asedy:=r21*asex+r22*asey+r23*asez;  //Z=9817
       asedz:=r31*asex+r32*asey+r33*asez;  //Z=9818
       bsedx:=r11*bsex+r12*bsey+r13*bsez;  //Z=9819
       bsedy:=r21*bsex+r22*bsey+r23*bsez;  //Z=9820
       bsedz:=r31*bsex+r32*bsey+r33*bsez;  //Z=9821
       csedx:=r11*csex+r12*csey+r13*csez;  //Z=9822
       csedy:=r21*csex+r22*csey+r23*csez;  //Z=9823
       csedz:=r31*csex+r32*csey+r33*csez; */  //Z=9824


    /*  n,u,v axis vectors in Carthesian coordinate system  */  //Z=9827
    /*  in xyz-system  */  //Z=9828
    nxyzx = 0;     nxyzy = 0;     nxyzz = -1;  //Z=9829
    uxyzx = 1;     uxyzy = 0;     uxyzz = 0;  //Z=9830
    vxyzx = 0;     vxyzy = 1;     vxyzz = 0;  //Z=9831

    /*  back-rotated n,u,v axis vectors  */  //Z=9833
    /*  nu=R n, uu=R u, vu=R v  in xyz-system  */  //Z=9834
    nuvwx = r11*nxyzx+r12*nxyzy+r13*nxyzz;  //Z=9835
    nuvwy = r21*nxyzx+r22*nxyzy+r23*nxyzz;  //Z=9836
    nuvwz = r31*nxyzx+r32*nxyzy+r33*nxyzz;  //Z=9837
    uuvwx = r11*uxyzx+r12*uxyzy+r13*uxyzz;  //Z=9838
    uuvwy = r21*uxyzx+r22*uxyzy+r23*uxyzz;  //Z=9839
    uuvwz = r31*uxyzx+r32*uxyzy+r33*uxyzz;  //Z=9840
    vuvwx = r11*vxyzx+r12*vxyzy+r13*vxyzz;  //Z=9841
    vuvwy = r21*vxyzx+r22*vxyzy+r23*vxyzz;  //Z=9842
    vuvwz = r31*vxyzx+r32*vxyzy+r33*vxyzz;  //Z=9843

    /*  reciprocal space n,u,v axis vectors  */  //Z=9845
    /*  nua=G nu, uua=G uu, vua=G vu  */  //Z=9846
    nhklx = g11*nuvwx+g12*nuvwy+g13*nuvwz;  //Z=9847
    nhkly = g21*nuvwx+g22*nuvwy+g23*nuvwz;  //Z=9848
    nhklz = g31*nuvwx+g32*nuvwy+g33*nuvwz;  //Z=9849
    uhklx = g11*uuvwx+g12*uuvwy+g13*uuvwz;  //Z=9850
    uhkly = g21*uuvwx+g22*uuvwy+g23*uuvwz;  //Z=9851
    uhklz = g31*uuvwx+g32*uuvwy+g33*uuvwz;  //Z=9852
    vhklx = g11*vuvwx+g12*vuvwy+g13*vuvwz;  //Z=9853
    vhkly = g21*vuvwx+g22*vuvwy+g23*vuvwz;  //Z=9854
    vhklz = g31*vuvwx+g32*vuvwy+g33*vuvwz;  //Z=9855

    nnhkl = sqrt(nhklx*nhklx+nhkly*nhkly+nhklz*nhklz+1e-10);  //Z=9857
    unhkl = sqrt(uhklx*uhklx+uhkly*uhkly+uhklz*uhklz+1e-10);  //Z=9858
    vnhkl = sqrt(vhklx*vhklx+vhkly*vhkly+vhklz*vhklz+1e-10);  //Z=9859

    nhklx = nhklx/nnhkl;  //Z=9861
    nhkly = nhkly/nnhkl;  //Z=9862
    nhklz = nhklz/nnhkl;  //Z=9863
    uhklx = uhklx/unhkl;  //Z=9864
    uhkly = uhkly/unhkl;  //Z=9865
    uhklz = uhklz/unhkl;  //Z=9866
    vhklx = vhklx/vnhkl;  //Z=9867
    vhkly = vhkly/vnhkl;  //Z=9868
    vhklz = vhklz/vnhkl;  //Z=9869


    /*  rotated n,u,v axis vectors  */  //Z=9872
    /*  nhkl=G nuvw, uhkl=G uuvw, vhkl=G vuvw  in unit cell uvw-system  */  //Z=9873
    /*  needed to indicate reciprocal space directions <hkl> in 2D-pattern  */  //Z=9874
    /* vols:=nuvwx*(uuvwy*vuvwz-uuvwz*vuvwy)+nuvwy*(uuvwz*vuvwx-uuvwx*vuvwz)+nuvwz*(uuvwx*vuvwy-uuvwy*vuvwx);  //Z=9875
       nhklx:=(uuvwy*vuvwz-uuvwz*vuvwy)/vols;  //Z=9876
       nhkly:=(uuvwz*vuvwx-uuvwx*vuvwz)/vols;  //Z=9877
       nhklz:=(uuvwx*vuvwy-uuvwy*vuvwx)/vols;  //Z=9878
       uhklx:=(nuvwz*vuvwy-nuvwy*vuvwz)/vols;  //Z=9879
       uhkly:=(nuvwx*vuvwz-nuvwz*vuvwx)/vols;  //Z=9880
       uhklz:=(nuvwy*vuvwx-nuvwx*vuvwy)/vols;  //Z=9881
       vhklx:=(nuvwy*uuvwz-nuvwz*uuvwy)/vols;  //Z=9882
       vhkly:=(nuvwz*uuvwx-nuvwx*uuvwz)/vols;  //Z=9883
       vhklz:=(nuvwx*uuvwy-nuvwy*uuvwx)/vols;   */  //Z=9884

    /* nhklx:=g11*nuvwx+g12*nuvwy+g13*nuvwz;  //Z=9886
       nhkly:=g21*nuvwx+g22*nuvwy+g23*nuvwz;  //Z=9887
       nhklz:=g31*nuvwx+g32*nuvwy+g33*nuvwz;  //Z=9888
       uhklx:=g11*uuvwx+g12*uuvwy+g13*uuvwz;  //Z=9889
       uhkly:=g21*uuvwx+g22*uuvwy+g23*uuvwz;  //Z=9890
       uhklz:=g31*uuvwx+g32*uuvwy+g33*uuvwz;  //Z=9891
       vhklx:=g11*vuvwx+g12*vuvwy+g13*vuvwz;  //Z=9892
       vhkly:=g21*vuvwx+g22*vuvwy+g23*vuvwz;  //Z=9893
       vhklz:=g31*vuvwx+g32*vuvwy+g33*vuvwz; */  //Z=9894

    /*  rotated reciprocal space base vectors  */  //Z=9896
    /*  aed=R ae, bed=R be, ced=R ce in xyz-system  */  //Z=9897
    /* asxd:=r11*asx+r12*asy+r13*asz;  //Z=9898
       asyd:=r21*asx+r22*asy+r23*asz;  //Z=9899
       aszd:=r31*asx+r32*asy+r33*asz;  //Z=9900
       bsxd:=r11*bsx+r12*bsy+r13*bsz;  //Z=9901
       bsyd:=r21*bsx+r22*bsy+r23*bsz;  //Z=9902
       bszd:=r31*bsx+r32*bsy+r33*bsz;  //Z=9903
       csxd:=r11*csx+r12*csy+r13*csz;  //Z=9904
       csyd:=r21*csx+r22*csy+r23*csz;  //Z=9905
       cszd:=r31*csx+r32*csy+r33*csz;   */  //Z=9906

    /* *************** */  //Z=9908
    /*  output matrix  */  //Z=9909
    m11 = asedxn;   m21 = bsedxn;   m31 = csedxn;  //Z=9910
    m12 = asedyn;   m22 = bsedyn;   m32 = csedyn;  //Z=9911
    m13 = asedzn;   m23 = bsedzn;   m33 = csedzn;  //Z=9912

    if ( twin )
    {/*2*/  //Z=9914
        mtw11 = asedxt;   mtw21 = bsedxt;   mtw31 = csedxt;  //Z=9915
        mtw12 = asedyt;   mtw22 = bsedyt;   mtw32 = csedyt;  //Z=9916
        mtw13 = asedzt;   mtw23 = bsedzt;   mtw33 = csedzt;  //Z=9917
    }/*2*/  //Z=9918

    /* (* test output *)  //Z=9920
       m11:=asxn;   m21:=bsxn;   m31:=csxn;  //Z=9921
       m12:=asyn;   m22:=bsyn;   m32:=csyn;  //Z=9922
       m13:=aszn;   m23:=bszn;   m33:=cszn;  //Z=9923

       mtw11:=asxt;   mtw21:=bsxt;   mtw31:=csxt;  //Z=9925
       mtw12:=asyt;   mtw22:=bsyt;   mtw32:=csyt;  //Z=9926
       mtw13:=aszt;   mtw23:=bszt;   mtw33:=cszt; */  //Z=9927

    if ( lat2d )
    {/*2*/  //Z=9929
        m11 = asedxDD;   m21 = bsedxDD;   m31 = csedxDD;  //Z=9930
        m12 = asedyDD;   m22 = bsedyDD;   m32 = csedyDD;  //Z=9931
        m13 = asedzDD;   m23 = bsedzDD;   m33 = csedzDD;  //Z=9932
        vvol = area;  //Z=9933
    }/*2*/  //Z=9934

    if ( lat1d )
    {/*2*/  //Z=9936
        m11 = asedxDL;   m21 = bsedxDL;   m31 = csedxDL;  //Z=9937
        m12 = asedyDL;   m22 = bsedyDL;   m32 = csedyDL;  //Z=9938
        m13 = asedzDL;   m23 = bsedzDL;   m33 = csedzDL;  //Z=9939

        m11 = u/(lenn*area1);  m21 = 0;   m31 = 0;  //Z=9941
        m12 = v/(lenn*area1);  m22 = 0;   m32 = 0;  //Z=9942
        m13 = w/(lenn*area1);  m23 = 0;   m33 = 0;  //Z=9943
        vvol = area1;  //Z=9944
    }/*2*/  //Z=9945

    /*10:*/  //Z=9947
}



// Bei den Parametern ist in /* */ der Variablenname im rufenden Programm angegeben.
void SasCalc_GENERIC_calculation::coefficients(int dim/*cdim*/, int nmax/*120*/, int ordis/*ordis*/, double &order/*order*/)
{/*1*/  //Z=9861

    static const int  nnmax = 120; // 500;  //Z=9865
    static const double min =  1E-40;  //Z=9866

    int    i,j,k,ll,m,n,ns,ms,n1,n2,n3,n4,n5,n6,n7,n8,n9,n1f,n2f,n3f,n4f,n5f,n6f,n7f,n8f,n9f,/*nmax2,*/cnt,inmax;
    int    /*ii,jj,kk,pp,qq,rr,iis,iie,jjs,jje,pps,ppe,qqs,qqe,*/ncell;
    double z,zz,zl,zzl,xlz,xl2z,xrz,xr2z,b1s,/*v,ee0,ee1,gb1s,preg1,pz2v,pz2v1,pz2v2,com,*/binsum/*,binsum1*/;
    double a1,intl,p,vvm,eps,area,vol; // p11,p12,p13,p21,p22,p23,p31,p32,p33 jetzt in params.*
    long double /*u1z,u1n,u2,u3,u4,u5,*/xrmz,xrm2z,sump,sump1,sump2,sumi,sumf,radv,aell,bell,cell;
    double carr2i[nnmax+1],fsum[nnmax+1]; //: Array[0..2*nnmax] of extended;  /*Z0311=9329*/;
    double lqm[10/*3*nnmax+1*/],qqm[10/*3*nnmax+1*/],phim[3*nnmax+1],llm[3*nnmax+1],rrm[3*nnmax+1],ppm[3*nnmax+1],
           a1m[3*nnmax+1]/*,a2m[3*nnmax+1]*/; //: Array[0..3*nnmax] of extended;  /*Z=9378*/
    double /*uell[3*nnmax+1],*/u1ell[3*nnmax+1],u2ell[3*nnmax+1],u3ell[3*nnmax+1],/*vell[3*nnmax+1],*/v1ell[3*nnmax+1],
           v2ell[3*nnmax+1],v3ell[3*nnmax+1],gell[3*nnmax+1],g3ell[3*nnmax+1]; //: Array[0..3*nnmax] of extended;  /*Z=9378*/
    double intlar[2*nnmax+1][2*nnmax+1]; //: Array[0..2*nnmax,0..2*nnmax] of extended;  /*Z0311=9331*/;
    //bool   search1,/*search2,*/search3,search4,search5,search6/*,search7,search8,search9*/;
    float  philiph,philipt,phiax,phiin,phiout,rad,lliph,llipt,lin,lout,len,rmax;  /*Z0311=9333*/;
    float  /*xradm,xrzm,*/x1zm,x12zm;  /*Z0311=9334*/;

    double b1sv_n = 1;

    // Arrays von oben mit erhöhter Genauigkeit
    long double gam3[nnmax+1], fkv[nnmax+1], fk2v[nnmax+1], fkvm[nnmax+1], pn[nnmax+1],xrn[nnmax+1], xln[nnmax+1], z12v[nnmax+1], z12vl[nnmax+1];

    // Arrays von oben, von denen aber nur das letzte Element genutzt wird.
    // Können also als normale Variablen geschrieben werden um Speicher zu sparen
    long double xrmn_n; // , fkvm_n;

    // TODO: Vielleicht sollten einige Arrays nur lokal angelegt werden, sobald diese auch verwendet werden.
    //       Es bleibt zu untersuchen, ob man diese Routine nicht auch auf die GPU portieren kann. Dann sollten
    //       einige Arrays durch die passende Mathematik ersetzt werden.

    // Sicherheitsabfrage: nmax (Parameter) wird für Schleifen genutzt (i.d.R. 120),
    //                     nnmax (const hier) wird für die Arrays genutzt (z.Zt. 120)
    if ( nnmax+1 < nmax )
    {
#ifndef __CUDACC__
        qDebug() << "coefficients() FATAL error: nmax" << nmax << "> nnmax" << nnmax;
#endif
        exit(255);
    }

    // globale Variablen, die jetzt hier definiert werden
    // TODO
    const double a = params.uca; // latt_a  oder  "EditLattice"

    // Variablen vorbesetzen um Compiler-Warnungen zu unterdrücken
    params.p11=1; params.p12=1; params.p13=1;
    params.p21=1; params.p22=1; params.p23=1;
    params.p31=1; params.p32=1; params.p33=1;

    z = (1-sqr(params.sigma))/sqr(params.sigma);  //Z=9891
    zl = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=9892
    zz = z;  //Z=9893
    zzl = zl;  //Z=9894
    DBG( qDebug() << "coefficients()" << "z=" << z << "zl=" << zl; )
    /* z:=z+2*dim;          (* z-average *)  //Z=9895 */
    /* zl:=zl+2*dim;  //Z=9896 */
    //nmax2 = ceil(nmax/2.0);  //Z=9897
    p = params.radius/params.radiusi;  //Z=9898
    eps = params.length/params.radius;  //Z=9899

    xlz = params.length/(2.0*(zzl+1));  //Z=9901
    xl2z = xlz*xlz;  //Z=9902
    xrz = params.radius/(2.0*(zz+1));  //Z=9903
    xr2z = xrz*xrz;  //Z=9904
    xrmz = params.radiusi/(2.0*(zz+1));  //Z=9905
    xrm2z = xrmz*xrmz;  //Z=9906

    /*  ellipsoid semiaxis  */  //Z=9908
    aell = params.radius;  //Z=9909
    bell = params.length;  //Z=9910
    cell = params.radiusi;  //Z=9911

    xln[0] = 1;  //Z=9913
    xrn[0] = 1;  //Z=9914
    xrmn_n = 1;  //Z=9915
    pn[0] = 1;  //Z=9916
    z12v[0] = 1;  //Z=9917
    z12vl[0] = 1;  //Z=9918
    fsum[0] = 1;  //Z=9919

    /*  factor for cross-sectional formfactors  */  //Z=9922
    if ( dim==1 ) b1s = 2;         /*  cylinder, cross-section is disk  */  //Z=9923
    if ( dim==2 ) b1s = 3/2.0;       /*  disk, cross-section is cylinder  */  //Z=9924
    if ( dim==3 ) b1s = 5/2.0;       /*  sphere  */  //Z=9925
    if ( dim==4 ) b1s = 3/2.0;       /*  cube  */  //Z=9926

    /*  start values for recursive parameters  */  //Z=9928
    //b1sv[0] = 1;  //Z=9929
    fkv[0] = 1;  //Z=9930
    fk2v[0] = 1;  //Z=9931
    //e1[0] = 1;  //Z=9932
    //gam1[0] = sqrt(M_PI);  //Z=9933
    //gam2[0] = 1.0;  //Z=9934
    gam3[0] = sqrt(M_PI)/2.0;  //Z=9935

    /*  initialize  */  //Z=9937
    /* for n:=0 to 10000 do begin  //Z=9938 */
    /*    carr1pm[n]:=1;  //Z=9939 */
    /*    carr2pm[n]:=1;  //Z=9940 */
    /* end;  //Z=9941 */
    for ( n=0; n<=150; n++ )        // Diese Initialisierungen sind in der neuen Version raus.
    {   //Z=9942
        params.CR->carr1p[n] = 1;      params.CR->carr1f[n] = 1;  //Z=9943
        params.CR->carr2p[n] = 1;      params.CR->carr2f[n] = 1;  //Z=9944
        params.CR->carr3p[n] = 1;      params.CR->carr3f[n] = 1;  //Z=9945
        params.CR->carr4p[n] = 1;      params.CR->carr4f[n] = 1;  //Z=9946
        params.CR->carr5p[n] = 1;      params.CR->carr5f[n] = 1;  //Z=9947
        params.CR->carr6p[n] = 1;      params.CR->carr6f[n] = 1;  //Z=9948
        params.CR->carr7p[n] = 1;      params.CR->carr7f[n] = 1;  //Z=9949
        params.CR->carr8p[n] = 1;      params.CR->carr8f[n] = 1;  //Z=9950
        params.CR->carr9p[n] = 1;      params.CR->carr9f[n] = 1;  //Z=9951
        //CR->carr2i[n] = 1;  //Z=9952
    }
    for ( n=0; n<=130; n++ )        // Diese Initialisierungen sind in der neuen Version raus.
    {   //Z=9954
        for ( m=0; m<=130; m++ )
        {   //Z=9955
            params.CR->carr11pm[n][m] = 1;  //Z=9956
            params.CR->carr22pm[n][m] = 1;  //Z=9957
        }
    }

    /*  multi-shell or liposome structure parameters  */  //Z=9961
    if ( params.cs==3 )
    {/*2*/  //Z=9962
        philiph = params.CR->myarray[12];     /*  water  */  //Z=9963
        philipt = params.CR->myarray[13];     /*  bilayer  */  //Z=9964

        rad = params.CR->myarray[1];          /*  vesicle inner radius  */  //Z=9966
        lliph = params.CR->myarray[7];        /*  water  */  //Z=9967
        llipt = params.CR->myarray[8];        /*  bilayer  */  //Z=9968

        len = lliph+llipt;  //Z=9970
        ncell = round(params.CR->myarray[4]);  //Z=9971
        rmax = rad+ncell*len;  //Z=9972

        lqm[1] = lliph;  //Z=9974
        lqm[2] = lqm[1]+llipt;  //Z=9975

        for ( i=1; i<=2; i++ ) qqm[i] = lqm[i]/len;  //Z=9977

        phim[1] = philiph;       /*  vesicle interior  */  //Z=9979
        llm[1] = 0.0;  //Z=9980
        rrm[1] = rad+llm[1];  //Z=9981
        ppm[1] = 1.0;  //Z=9982

        radv = rrm[1];        /*  vesicle radius  */  //Z=9984
        cnt = 1;  //Z=9985
        for ( i=1; i<=ncell; i++ )
        {/*3*/  //Z=9986
            phim[cnt+1] = philipt;         /*  bilayer  */  //Z=9987
            phim[cnt+2] = philiph;         /*  water  */  //Z=9988
            llm[cnt+1] = (i-1+qqm[1])*len;  //Z=9989
            llm[cnt+2] = (i-1+qqm[2])*len;  //Z=9990
            rrm[cnt+1] = radv+llm[cnt+1];  //Z=9991
            rrm[cnt+2] = radv+llm[cnt+2];  //Z=9992
            ppm[cnt+1] = rrm[cnt+1]/rad;  //Z=9993
            ppm[cnt+2] = rrm[cnt+2]/rad;  //Z=9994
            cnt = cnt+2;  //Z=9995
        }/*3*/  //Z=9996
        inmax = cnt;  //Z=9997
        phim[cnt+1] = 0.0;  //Z=9998

        //xradm = rad;  //Z=10000
        //xrzm = xradm/(z+1);  //Z=10001
        x1zm = rad/(2.0*(z+1));  //Z=10002
        x12zm = x1zm*x1zm;  //Z=10003
        //xmax = q*rmax;  //Z=10004

        for ( i=1; i<=inmax; i++ )
        {/*3*/  //Z=10006
            if ( params.part==0 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],3); /*  spheres  */  //Z=10007
            if ( params.part==1 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],2); /*  cylinders  */  //Z=10008
            if ( params.part==2 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],1); /*  disks  */  //Z=10009
            params.CR->carr7p[i] = ppm[i];  //Z=10010
            params.CR->carr3p[i] = llm[i];  //Z=10011
            params.CR->carr5p[i] = a1m[i];  //Z=10012
        }/*3*/  //Z=10013

        fkvm[0] = 1;  //Z=10015
        gam3[0] = sqrt(M_PI)/2.0;  //Z=10016
        for ( n=1; n<=nmax; n++ )
        {/*3*/  //Z=10017
            if ( params.part==0 )
            {/*4*/       /*  spheres  */  //Z=10018
                fkvm[n] = fkvm[n-1]*n;  //Z=10019
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10020
                params.CR->carr6p[n] = (n+3/2.0)*gam3[n]*fkvm[n]*4/(3.0*sqrt(M_PI));  //Z=10021
            }/*4*/  //Z=10022
            if ( params.part==1 )
            {/*4*/       /*  cylinders  */  //Z=10023
                fkvm[n] = fkvm[n-1]*n;  //Z=10024
                params.CR->carr6p[n] = (n+1)*fkvm[n]*fkvm[n];  //Z=10025
            }/*4*/  //Z=10026
            if ( params.part==2 )
            {/*4*/       /*  disks  */  //Z=10027
                fkvm[n] = fkvm[n-1]*n;  //Z=10028
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10029
                params.CR->carr6p[n] = gam3[n]*fkvm[n]*2/sqrt(M_PI);  //Z=10030
            }/*4*/  //Z=10031
        }/*3*/  //Z=10032

        vvm = 0;  //Z=10034
        for ( i=1; i<=inmax; i++ )
        {/*3*/  //Z=10035
            for ( j=1; j<=inmax; j++ ) vvm = vvm+a1m[i]*a1m[j];  //Z=10036
        }/*3*/  //Z=10037

        params.CR->myarray[14] = inmax;  //Z=10039
        params.CR->myarray[15] = vvm;  //Z=10040
        params.CR->myarray[16] = rmax;  //Z=10041
    }/*2*/  //Z=10042


    /*  myelin structure parameters  */  //Z=10045
    if ( params.cs==4 )
    {/*2*/  //Z=10046
        philiph = params.CR->myarray[12];     /*  head group  */  //Z=10047
        philipt = params.CR->myarray[13];     /*  tail  */  //Z=10048
        phiax = params.CR->myarray[9];        /*  axon  */  //Z=10049
        phiin = params.CR->myarray[10];       /*  intra cell  */  //Z=10050
        phiout = params.CR->myarray[11];      /*  extra cell  */  //Z=10051

        rad = params.CR->myarray[1];          /*  vesicle inner radius  */  //Z=10053
        lliph = params.CR->myarray[7];        /*  head group  */  //Z=10054
        llipt = params.CR->myarray[8];        /*  tail  */  //Z=10055
        lin = params.CR->myarray[6];          /*  intra cell  */  //Z=10056
        lout = params.CR->myarray[5];         /*  extra cell  */  //Z=10057

        len = lout+2*(2*lliph+llipt)+lin;  //Z=10059
        ncell = round(params.CR->myarray[4]);  //Z=10060
        rmax = rad+ncell*len;  //Z=10061

        lqm[1] = lout;  //Z=10063
        lqm[2] = lqm[1]+lliph;  //Z=10064
        lqm[3] = lqm[2]+llipt;  //Z=10065
        lqm[4] = lqm[3]+lliph;  //Z=10066
        lqm[5] = lqm[4]+lin;  //Z=10067
        lqm[6] = lqm[5]+lliph;  //Z=10068
        lqm[7] = lqm[6]+llipt;  //Z=10069
        lqm[8] = lqm[7]+lliph;  //Z=10070

        for ( i=1; i<=8; i++ ) qqm[i] = lqm[i]/len;  //Z=10072

        phim[1] = phiax;       /*  vesicle interior  */  //Z=10074
        llm[1] = 0.0;  //Z=10075
        rrm[1] = rad+llm[1];  //Z=10076
        ppm[1] = 1.0;  //Z=10077
        phim[2] = philiph;     /*  vesicle bilayer: head group  */  //Z=10078
        llm[2] = lliph;  //Z=10079
        rrm[2] = rad+llm[2];  //Z=10080
        ppm[2] = rrm[2]/rad;  //Z=10081
        phim[3] = philipt;     /*  vesicle bilayer: tail group  */  //Z=10082
        llm[3] = llm[2]+llipt;  //Z=10083
        rrm[3] = rad+llm[3];  //Z=10084
        ppm[3] = rrm[3]/rad;  //Z=10085
        phim[4] = philiph;     /*  vesicle bilayer: head group  */  //Z=10086
        llm[4] = llm[3]+lliph;  //Z=10087
        rrm[4] = rad+llm[4];  //Z=10088
        ppm[4] = rrm[4]/rad;  //Z=10089

        radv = rrm[4];        /*  vesicle radius + bilayer  */  //Z=10091
        cnt = 4;  //Z=10092
        for ( i=1; i<=ncell; i++ )
        {/*3*/  //Z=10093
            phim[cnt+1] = phiout;          /*  extra cell  */  //Z=10094
            phim[cnt+2] = philiph;         /*  head group  */  //Z=10095
            phim[cnt+3] = philipt;         /*  tail group  */  //Z=10096
            phim[cnt+4] = philiph;         /*  head group  */  //Z=10097
            phim[cnt+5] = phiin;           /*  intra cell  */  //Z=10098
            phim[cnt+6] = philiph;         /*  head group  */  //Z=10099
            phim[cnt+7] = philipt;         /*  tail group  */  //Z=10100
            phim[cnt+8] = philiph;         /*  head group  */  //Z=10101
            llm[cnt+1] = (i-1+qqm[1])*len;  //Z=10102
            llm[cnt+2] = (i-1+qqm[2])*len;  //Z=10103
            llm[cnt+3] = (i-1+qqm[3])*len;  //Z=10104
            llm[cnt+4] = (i-1+qqm[4])*len;  //Z=10105
            llm[cnt+5] = (i-1+qqm[5])*len;  //Z=10106
            llm[cnt+6] = (i-1+qqm[6])*len;  //Z=10107
            llm[cnt+7] = (i-1+qqm[7])*len;  //Z=10108
            llm[cnt+8] = (i-1+qqm[8])*len;  //Z=10109
            rrm[cnt+1] = radv+llm[cnt+1];  //Z=10110
            rrm[cnt+2] = radv+llm[cnt+2];  //Z=10111
            rrm[cnt+3] = radv+llm[cnt+3];  //Z=10112
            rrm[cnt+4] = radv+llm[cnt+4];  //Z=10113
            rrm[cnt+5] = radv+llm[cnt+5];  //Z=10114
            rrm[cnt+6] = radv+llm[cnt+6];  //Z=10115
            rrm[cnt+7] = radv+llm[cnt+7];  //Z=10116
            rrm[cnt+8] = radv+llm[cnt+8];  //Z=10117
            ppm[cnt+1] = rrm[cnt+1]/rad;  //Z=10118
            ppm[cnt+2] = rrm[cnt+2]/rad;  //Z=10119
            ppm[cnt+3] = rrm[cnt+3]/rad;  //Z=10120
            ppm[cnt+4] = rrm[cnt+4]/rad;  //Z=10121
            ppm[cnt+5] = rrm[cnt+5]/rad;  //Z=10122
            ppm[cnt+6] = rrm[cnt+6]/rad;  //Z=10123
            ppm[cnt+7] = rrm[cnt+7]/rad;  //Z=10124
            ppm[cnt+8] = rrm[cnt+8]/rad;  //Z=10125
            cnt = cnt+8;  //Z=10126
        }/*3*/  //Z=10127
        inmax = cnt;  //Z=10128
        phim[cnt+1] = 0.0;  //Z=10129

        //xradm = rad;  //Z=10131
        //xrzm = xradm/(z+1);  //Z=10132
        x1zm = rad/(2.0*(z+1));  //Z=10133
        x12zm = x1zm*x1zm;  //Z=10134
        //xmax = q*rmax;  //Z=10135

        for ( i=1; i<=inmax; i++ )
        {/*3*/  //Z=10137
            if ( params.part==0 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],3); /*  spheres  */  //Z=10138
            if ( params.part==1 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],2); /*  cylinders  */  //Z=10139
            if ( params.part==2 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],1); /*  disks  */  //Z=10140
            params.CR->carr7p[i] = ppm[i];  //Z=10141
            params.CR->carr3p[i] = llm[i];  //Z=10142
            params.CR->carr5p[i] = a1m[i];  //Z=10143
        }/*3*/  //Z=10144

        fkvm[0] = 1;  //Z=10146
        gam3[0] = sqrt(M_PI)/2.0;  //Z=10147
        for ( n=1; n<=nmax; n++ )
        {/*3*/  //Z=10148
            if ( params.part==0 )
            {/*4*/       /*  spheres  */  //Z=10149
                fkvm[n] = fkvm[n-1]*n;  //Z=10150
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10151
                params.CR->carr6p[n] = (n+3/2.0)*gam3[n]*fkvm[n]*4/(3.0*sqrt(M_PI));  //Z=10152
            }/*4*/  //Z=10153
            if ( params.part==1 )
            {/*4*/       /*  cylinders  */  //Z=10154
                fkvm[n] = fkvm[n-1]*n;  //Z=10155
                params.CR->carr6p[n] = (n+1)*fkvm[n]*fkvm[n];  //Z=10156
            }/*4*/  //Z=10157
            if ( params.part==2 )
            {/*4*/       /*  disks  */  //Z=10158
                fkvm[n] = fkvm[n-1]*n;  //Z=10159
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10160
                params.CR->carr6p[n] = gam3[n]*fkvm[n]*2/sqrt(M_PI);  //Z=10161
            }/*4*/  //Z=10162
        }/*3*/  //Z=10163

        vvm = 0;  //Z=10165
        for ( i=1; i<=inmax; i++ )
        {/*3*/  //Z=10166
            for ( j=1; j<=inmax; j++ ) vvm = vvm+a1m[i]*a1m[j];  //Z=10167
        }/*3*/  //Z=10168

        params.CR->myarray[14] = inmax;  //Z=10170
        params.CR->myarray[15] = vvm;  //Z=10171
        params.CR->myarray[16] = rmax;  //Z=10172
    }/*2*/  //Z=10173



    //search1 = true;  //Z=10176
    //search2 = true;  //Z=10177
    //search3 = true;  //Z=10178
    //search4 = true;  //Z=10179
    //search5 = true;  //Z=10180
    //search6 = true;  //Z=10181
    //search7 = true;  //Z=10182
    //search8 = true;  //Z=10183
    //search9 = true;  //Z=10184
    n1 = nmax;      n1f = nmax;  //Z=10185
    n2 = nmax;      n2f = nmax;  //Z=10186
    n3 = nmax;      n3f = nmax;  //Z=10187
    n4 = nmax;      n4f = nmax;  //Z=10188
    n5 = nmax;      n5f = nmax;  //Z=10189
    n6 = nmax;      n6f = nmax;  //Z=10190
    n7 = nmax;      n7f = nmax;  //Z=10191
    n8 = nmax;      n8f = nmax;  //Z=10192
    n9 = nmax;      n9f = nmax;  //Z=10193

    /*  cho1 = orientation case  */  //Z=10195
    params.orcase = 1;                                                         /*  general  */  //Z=10196
    if ( (params.polPhi== 0) && (params.polTheta==90) ) params.orcase = 2;     /*  x-axis  */  //Z=10197
    if ( (params.polPhi==90) && (params.polTheta==90) ) params.orcase = 3;     /*  y-axis  */  //Z=10198
    if ( (params.polPhi== 0) && (params.polTheta== 0) ) params.orcase = 4;     /*  z-axis  */  //Z=10199
    if ( (params.polPhi==90) && (params.polTheta== 0) ) params.orcase = 4;     /*  z-axis  */  //Z=10200

#ifndef __CUDACC__
    qDebug() << "coefficients()" << "ordis"<<ordis << "dim"<<dim << "part"<<params.part << "cs"<<params.cs
             << "cho1|orcase"<<params.orcase << params.polPhi << params.polTheta;
#endif
    //Debug: coefficients ordis 0="CBOrdis.Gaussian"
    //                    dim 1=im Vorfeld aus CBPart.Cylinder bestimmt
    //                    part 1="CBPart.Cylinder"
    //                    nmax 120="120" festgeschieben
    //                    cs 0="CBInt.homogeneous" - geht von 0 bis 2, hier wird aber 3 und mehr genutzt TODO


    /* ** isotropic case for spheres ** */  //Z=10203
    if ( dim==3 )
    {/*2*/  //Z=10204
        params.norm = 1;  //Z=10205
        order = 0;  //Z=10206
        /*  homogeneous  */  //Z=10207
        if ( params.cs==0 )
        {/*3*/  //Z=10208

            //[10310]         (* monodisperse sphere test *)
            //[10311]         {xrn[0]:=1;
            //[10312]         carr4p[0]:=1;
            //[10313]         for n:=1 to nmax do begin
            //[10314]            xrn[n]:=xrn[n-1]*(-r*r)/((n+3/2)*n);
            //[10315]            carr4p[n]:=6*xrn[n]/((n+3)*(n+2));
            //[10316]            if (abs(carr4p[n])<min) then begin
            //[10317]               if n<n4 then n4:=n;
            //[10318]            end;
            //[10319]         end;  }

            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10209
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10210
                fkv[n] = fkv[n-1]*n;  //Z=10211
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10212
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10213 */
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10214
                /*  P(q)-coefficient  */  //Z=10215
                params.CR->carr4p[n] = 9.*sqrt(M_PI)*pow(4.0,n)*(z12v[n]*xrn[n])/(2.0*(n+3.)*(n+2.)*(n+3./2.0)*gam3[n]*fkv[n]);  //Z=10216
                /*  F(q)-coefficient  */  //Z=10217
                binsum = 0.0;  //Z=10218
                for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+3./2.0)*gam3[m]*(n-m+3./2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10219
                params.CR->carr4f[n] = 9*M_PI*xrn[n]*binsum/16.0;  //Z=10220
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10221
                    if ( n<n4 ) n4 = n;  //Z=10222
                }/*5*/  //Z=10223
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10224
                    if ( n<n4f ) n4f = n;  //Z=10225
                }/*5*/  //Z=10226
            }/*4*/  //Z=10227
//#ifndef __CUDACC__
//            qDebug() << "              " << "n4"<<n4<<params.CR->carr4p[n4] << "n4f"<<n4f<<params.CR->carr4f[n4f] << "#" << params.CR->carr4f[n4];
//#endif
            goto Label99;  //Z=10228
        }/*3*/   /*  of homogeneous  */  //Z=10229

        /*  core/shell  */  //Z=10231
        if ( params.cs==1 )
        {/*3*/  //Z=10232
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10233
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10234
                fkv[n] = fkv[n-1]*n;  //Z=10235
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10236
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10237 */
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10238
                xrmn_n = -xrmn_n*xrm2z;  //Z=10239
                pn[n] = pn[n-1]*p*p;  //Z=10240
                /*  P(q)-coefficients  */  //Z=10241
                sump = 0.0;  //Z=10242
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10243
                    sump = sump+pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10244
                }/*5*/  //Z=10245
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn[n]/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10246
                params.CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;  //Z=10247
                params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=10248
                /*  F(q)-coefficients  */  //Z=10249
                sump = 0.0;  //Z=10250
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10251
                    sump = sump+z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10252
                }/*5*/  //Z=10253
                params.CR->carr4f[n] = (9*M_PI/16.0)*xrn[n]*sump;  //Z=10254
                sump = 0.0;  //Z=10255
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10256
                    sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10257
                }/*5*/  //Z=10258
                params.CR->carr5f[n] = (9*M_PI/16.0)*xrmn_n*sump;  //Z=10259
                params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=10260
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10261
                    if ( n<n4 ) n4 = n;  //Z=10262
                }/*5*/  //Z=10263
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=10264
                    if ( n<n5 ) n5 = n;  //Z=10265
                }/*5*/  //Z=10266
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=10267
                    if ( n<n6 ) n6 = n;  //Z=10268
                }/*5*/  //Z=10269
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10270
                    if ( n<n4f ) n4f = n;  //Z=10271
                }/*5*/  //Z=10272
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=10273
                    if ( n<n5f ) n5f = n;  //Z=10274
                }/*5*/  //Z=10275
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=10276
                    if ( n<n6f ) n6f = n;  //Z=10277
                }/*5*/  //Z=10278
            }/*4*/  //Z=10279
            goto Label99;  //Z=10280
        }/*3*/   /*  of core/shell  */  //Z=10281

        /*  inhomogeneous core/shell  */  //Z=10283
        if ( params.cs==2 )
        {/*3*/  //Z=10284
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10285
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10286
                fkv[n] = fkv[n-1]*n;  //Z=10287
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10288
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10289 */
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10290
                xrmn_n = -xrmn_n*xrm2z;  //Z=10291
                pn[n] = pn[n-1]*p*p;  //Z=10292
                /*  P(q)-coefficients  */  //Z=10293
                params.CR->carr1p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn[n]/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10294
                sump = 0.0;  //Z=10295
                sump1 = 0.0;  //Z=10296
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10297
                    sumi = 1/((n-m+3/2.0)*gam3[n-m]*(m+3/2.0-params.alphash1/2.0)*gam3[m]*fkv[m]*fkv[n-m]);  //Z=10298
                    sump = sump+pn[n-m]*sumi;  //Z=10299
                    sump1 = sump1+sumi;  //Z=10300
                }/*5*/  //Z=10301
                params.CR->carr2p[n] = (3*M_PI*(3-params.alphash1)/16.0)*z12v[n]*xrmn_n*sump;  //Z=10302
                params.CR->carr3p[n] = (3*M_PI*(3-params.alphash1)/16.0)*z12v[n]*xrn[n]*sump1;  //Z=10303
                sump = 0.0;  //Z=10304
                sump1 = 0.0;  //Z=10305
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10306
                    sumi = 1/((n-m+3/2.0-params.alphash1/2.0)*(m+3/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=10307
                    sump = sump+sumi;  //Z=10308
                    sump1 = sump1+pn[n-m]*sumi;  //Z=10309
                }/*5*/  //Z=10310
                params.CR->carr4p[n] = (sqr(3-params.alphash1)*M_PI/16.0)*z12v[n]*xrmn_n*sump;  //Z=10311
                params.CR->carr5p[n] = (sqr(3-params.alphash1)*M_PI/16.0)*z12v[n]*xrmn_n*sump1;  //Z=10312
                params.CR->carr6p[n] = (sqr(3-params.alphash1)*M_PI/16.0)*z12v[n]*xrn[n]*sump;  //Z=10313

                /*  F(q)-coefficients  */  //Z=10315
                params.CR->carr4f[n] = (3*sqrt(M_PI)/4.0)*z12v[n]*xrn[n]/((n+3/2.0)*gam3[n]*fkv[n]);  //Z=10316
                params.CR->carr5f[n] = (sqrt(M_PI)*(3-params.alphash1)/4.0)*z12v[n]*xrmn_n/((n+3/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=10317
                params.CR->carr6f[n] = (sqrt(M_PI)*(3-params.alphash1)/4.0)*z12v[n]*xrn[n]/((n+3/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=10318
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=10319
                    if ( n<n1 ) n1 = n;  //Z=10320
                }/*5*/  //Z=10321
                if ( fabs(params.CR->carr2p[n])<min )
                {/*5*/  //Z=10322
                    if ( n<n2 ) n2 = n;  //Z=10323
                }/*5*/  //Z=10324
                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=10325
                    if ( n<n3 ) n3 = n;  //Z=10326
                }/*5*/  //Z=10327
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10328
                    if ( n<n4 ) n4 = n;  //Z=10329
                }/*5*/  //Z=10330
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=10331
                    if ( n<n5 ) n5 = n;  //Z=10332
                }/*5*/  //Z=10333
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=10334
                    if ( n<n6 ) n6 = n;  //Z=10335
                }/*5*/  //Z=10336
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10337
                    if ( n<n4f ) n4f = n;  //Z=10338
                }/*5*/  //Z=10339
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=10340
                    if ( n<n5f ) n5f = n;  //Z=10341
                }/*5*/  //Z=10342
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=10343
                    if ( n<n6f ) n6f = n;  //Z=10344
                }/*5*/  //Z=10345
            }/*4*/  //Z=10346
            goto Label99;  //Z=10347
        }/*3*/   /*  of inhomogeneous core/shell  */  //Z=10348

        /*  myelin  */  //Z=10350
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=10351
            i = 2;  //Z=10352
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10353
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10354
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10355
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=10356
                fkv[n] = fkv[n-1]*n;  //Z=10357
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10358
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10359
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10360 */
                xln[n] = -xln[n-1]*xl2z;  //Z=10361
                /* xrn[n]:=-xrn[n-1]*xr2z;  //Z=10362 */
                xrn[n] = -xrn[n-1]*x12zm;         /*  myelin radius  */  //Z=10363

                /*  P(q)  */  //Z=10365
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10366
                    /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=10367 */
                    /* i:=i+1;  //Z=10368 */
                    params.CR->carr11pm[n][m] = (9*M_PI/16.0)*(1/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]));  //Z=10369
                }/*5*/  //Z=10370
                params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=10371
                /* carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=10372 */

                /*  F(q)  */  //Z=10375
                binsum = 0.0;  //Z=10376
                for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10377
                params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=10378
                binsum = 0.0;  //Z=10379
                for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=10380
                params.CR->carr4f[n] = xrn[n]*binsum;  //Z=10381

                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=10384
                    if ( n<n1 ) n1 = n;  //Z=10385
                }/*5*/  //Z=10386
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=10387
                    if ( n<n1f ) n1f = n;  //Z=10388
                }/*5*/  //Z=10389
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10390
                    if ( n<n4 ) n4 = n;  //Z=10391
                }/*5*/  //Z=10392
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10393
                    if ( n<n4f ) n4f = n;  //Z=10394
                }/*5*/  //Z=10395
            }/*4*/  //Z=10396
            goto Label99;
        }/*3*/ /*  of myelin  */  //Z=10397

    }/*2*/  /*  of dim=3, spheres  */  //Z=10401


    /*  isotropic case for cubes  */  //Z=10403
    if ( (ordis==7) && (dim==4) )
    {/*2*/    /*  cubes  */  //Z=10404
        params.norm = 1;  //Z=10405
        order = 0;  //Z=10406
        /*  homogeneous  */  //Z=10407
        if ( params.cs==0 )
        {/*3*/  //Z=10408
            area = 6*4*sqr(params.radius);  //Z=10409
            vol = 8*params.radius*sqr(params.radius);  //Z=10410
            params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10411

#ifndef __CUDACC__
            qDebug() << "coeff() ordis=7, dim=4, cs=0, por:" << params.por;
#endif

            u1ell[0] = 2;  //Z=10413
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10414
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10415
                fkv[n] = fkv[n-1]*n;  //Z=10416
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10417
                /* u1ell[n]:=z12v[n]/((n+1/2)*(n+1)*fkv[n]);  //Z=10418 */
                u1ell[n] = 1/((n+1/2.0)*(n+1)*fkv[n]);  //Z=10419
            }/*4*/  //Z=10420

            for ( n=0; n<=nmax; n++ )
            {/*4*/  //Z=10422
                sump1 = 0.0;  //Z=10423
                for ( m=0; m<=n; m++ ) sump1 = sump1+u1ell[n-m]*u1ell[m];  //Z=10424
                v1ell[n] = sump1;  //Z=10425
            }/*4*/  //Z=10426

            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10428
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10429
                sumf = 0.0;  //Z=10430
                for ( m=0; m<=n; m++ ) sumf = sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10431
                fsum[n] = sumf;  //Z=10432
                /*  P(q)-coefficient  */  //Z=10433
                sump = 0.0;  //Z=10434
                /* for m:=0 to n do begin  //Z=10435
                   sump1:=0.0;  //Z=10436
                   for k:=0 to m do sump1:=sump1+z12v[m-k]*z12v[k]/((k+1/2)*(m-k+1/2)*(k+1)*fkv[k]*(m-k+1)*fkv[m-k]);  //Z=10437
                   sump:=sump+z12v[n-m]*sump1/((n-m+1/2)*(n-m+1)*fkv[n-m]);  //Z=10438
                   end; */  //Z=10439

                /* for m:=0 to n do begin  //Z=10441
                   sump1:=0.0;  //Z=10442
                   for k:=0 to m do sump1:=sump1+u1ell[m-k]*u1ell[k];  //Z=10443
                   sump:=sump+u1ell[n-m]*sump1;  //Z=10444
                   end; */  //Z=10445

                for ( m=0; m<=n; m++ ) sump = sump+u1ell[n-m]*v1ell[m];  //Z=10447

                /* carr4p[n]:=sqrt(pi)*power(4,n)*xrn[n]*sump/(16*gam3[n]);  //Z=10449 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn[n]*sump/(16.0*gam3[n]);  //Z=10450

                /*  F(q)-coefficient  */  //Z=10452
                sump = 0.0;  //Z=10453
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10454
                    sump1 = 0.0;  //Z=10455
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10456
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10457
                }/*5*/  //Z=10458
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]*sump/(128.0*gam3[n]);  //Z=10459
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10460
                    if ( n<n4 ) n4 = n;  //Z=10461
                }/*5*/  //Z=10462
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10463
                    if ( n<n4f ) n4f = n;  //Z=10464
                }/*5*/  //Z=10465

//#ifndef __CUDACC__
//                qDebug() << "coeff() carr4p[]:" << n << params.CR->carr4p[n] << n4;
//#endif

            }/*4*/  //Z=10466
            goto Label99;  //Z=10467
        }/*3*/   /*  of homogeneous  */  //Z=10468

        /*  core/shell  */          /*  not yet ready  */  //Z=10470
        if ( params.cs==1 )
        {/*3*/  //Z=10471
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10472
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10473
                fkv[n] = fkv[n-1]*n;  //Z=10474
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10475
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10476 */
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10477
                xrmn_n = -xrmn_n*xrm2z;  //Z=10478
                pn[n] = pn[n-1]*p*p;  //Z=10479
                /*  P(q)-coefficient  */  //Z=10480
                sump = 0.0;  //Z=10481
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10482
                    sump = sump+pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10483
                }/*5*/  //Z=10484
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn[n]/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10485
                params.CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;  //Z=10486
                params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=10487
                /*  F(q)-coefficient  */  //Z=10488
                /* carr3[n]:=3*sqrt(pi)*z12v[n]*xrn[n]/(4*(n+3/2)*gam3[n]*fkv[n]);  //Z=10489 */
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10490
                    if ( n<n4 ) n4 = n;  //Z=10491
                }/*5*/  //Z=10492
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=10493
                    if ( n<n5 ) n5 = n;  //Z=10494
                }/*5*/  //Z=10495
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=10496
                    if ( n<n6 ) n6 = n;  //Z=10497
                }/*5*/  //Z=10498
            }/*4*/  //Z=10499
            goto Label99;  //Z=10500
        }/*3*/   /*  of core/shell  */  //Z=10501
    }/*2*/  /*  of dim=4, cubes  */  //Z=10502

    /*  perfect orientation case for cubes  */  //Z=10504
    if ( (ordis==6) && (dim==4) )
    {/*2*/    /*  cubes  */  //Z=10505
        params.norm = 1;  //Z=10506
        order = 1;  //Z=10507
        /*  homogeneous  */  //Z=10508
        if ( params.cs==0 )
        {/*3*/  //Z=10509
            i = 2;  //Z=10510
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10511
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10512
                fkv[n] = fkv[n-1]*n;  //Z=10513
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10514
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10515 */
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10516
                /*  P(q)-coefficient  */  //Z=10517
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10518
                    /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10519 */
                    /* carr1fm[i]:=carr1pm[i];  //Z=10520 */
                    params.CR->carr11pm[n][m] = (M_PI/4.0)*pow(4.0,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10521
                    i = i+1;  //Z=10522
                }/*5*/  //Z=10523
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn[n]/(16.0*gam3[n]);  //Z=10524
                /*  F(q)-coefficient  */  //Z=10525
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=10526
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10527
                    if ( n<n4 ) n4 = n;  //Z=10528
                }/*5*/  //Z=10529
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10530
                    if ( n<n4f ) n4f = n;  //Z=10531
                }/*5*/  //Z=10532
            }/*4*/  //Z=10533
            goto Label99;  //Z=10534
        }/*3*/   /*  of homogeneous  */  //Z=10535

        /*  core/shell  */          /*  not yet ready  */  //Z=10537
        if ( params.cs==1 )
        {/*3*/  //Z=10538
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10539
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10540
                fkv[n] = fkv[n-1]*n;  //Z=10541
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10542
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10543 */
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10544
                xrmn_n = -xrmn_n*xrm2z;  //Z=10545
                pn[n] = pn[n-1]*p*p;  //Z=10546
                /*  P(q)-coefficient  */  //Z=10547
                sump = 0.0;  //Z=10548
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10549
                    sump = sump+pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10550
                }/*5*/  //Z=10551
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn[n]/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10552
                params.CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;  //Z=10553
                params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=10554
                /*  F(q)-coefficient  */  //Z=10555
                /* carr3[n]:=3*sqrt(pi)*z12v[n]*xrn[n]/(4*(n+3/2)*gam3[n]*fkv[n]);  //Z=10556 */
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10557
                    if ( n<n4 ) n4 = n;  //Z=10558
                }/*5*/  //Z=10559
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=10560
                    if ( n<n5 ) n5 = n;  //Z=10561
                }/*5*/  //Z=10562
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=10563
                    if ( n<n6 ) n6 = n;  //Z=10564
                }/*5*/  //Z=10565
            }/*4*/  //Z=10566
            goto Label99;  //Z=10567
        }/*3*/   /*  of core/shell  */  //Z=10568
    }/*2*/  /*  of dim=4, cubes  */  //Z=10569


    /*  isotropic case for ellipsoids  */  //Z=10571
    if ( (ordis==7) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10572
        params.norm = 1;  //Z=10573
        order = 0;  //Z=10574
        // Die Abfrage eps==1 sollte nicht so vorkommen
        if ( eps>1 )
            area = 2*M_PI*params.radius*(params.radius+(sqr(params.length)/sqrt(sqr(params.length)-sqr(params.radius)))*asin(sqrt(sqr(params.length)-sqr(params.radius))/params.length));  //Z=10576
        else if ( eps<1 )
            area = 2*M_PI*params.radius*(params.radius+(sqr(params.length)/sqrt(sqr(params.radius)-sqr(params.length)))*asinh(sqrt(sqr(params.radius)-sqr(params.length))/params.length));  //Z=10577
        else //if ( eps==1 )
            area = 4*M_PI*sqr(params.radius);  //Z=10575
        vol = (4*M_PI/3.0)*sqr(params.radius)*params.length;  //Z=10578
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10579

        /*  homogeneous  */  //Z=10581
        if ( params.cs==0 )
        {/*3*/  //Z=10582

            //[10695]         (* monodisperse ellipsoid test *)
            //[10696]         {e1[0]:=0;
            //[10697]         xrn[0]:=1;
            //[10698]         carr4p[0]:=1;
            //[10699]         carr11pm[0,0]:=1.0;
            //[10700]         for n:=1 to nmax do begin
            //[10701]            sump:=1.0;
            //[10702]            carr11pm[n,0]:=1.0;
            //[10703]            for m:=1 to n do begin
            //[10704]               carr11pm[n,m]:=carr11pm[n,m-1]*((n+1-m)/m)*(eps*eps-1
            //[10705]               sump:=sump+carr11pm[n,m]/(2*m+1);
            //[10706]            end;
            //[10707]            e1[n]:=sump;
            //[10708]            xrn[n]:=xrn[n-1]*(-r*r)/((n+3/2)*n);
            //[10709]            carr4p[n]:=6*e1[n]*xrn[n]/((n+3)*(n+2));
            //[10710]            if (abs(carr4p[n])<min) then begin
            //[10711]               if n<n4 then n4:=n;
            //[10712]            end;
            //[10713]         end; }

            double e1[nnmax+1];  // lokal angelegt, da nur hier verwendet
            e1[0] = 1;
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10583
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10584
                fkv[n] = fkv[n-1]*n;  //Z=10585
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10586
                e1[n] = e1[n-1]*(eps*eps-1);  //Z=10587
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10588
                sump = 0.0;  //Z=10589
                for ( m=0; m<=n; m++ ) sump = sump+e1[n-m]/(fkv[n-m]*fkv[m]*(2*(n-m)+1));  //Z=10590
                /* sumf:=0.0;  //Z=10591 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10592 */
                /* fsum[n]:=sumf;  //Z=10593 */
                /*  P(q)-coefficient  */  //Z=10594
                params.CR->carr4p[n] = (9*sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);  //Z=10595
                /*  F(q)-coefficient  */  //Z=10596
                sump = 0.0;  //Z=10597
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10598
                    sump1 = 0.0;  //Z=10599
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10600
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10601
                }/*5*/  //Z=10602
                params.CR->carr4f[n]  = M_PI*M_PI*xrn[n]*sump/(128.0*gam3[n]);  //Z=10603
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10604
                    if ( n<n4 ) n4 = n;  //Z=10605
                }/*5*/  //Z=10606
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10607
                    if ( n<n4f ) n4f = n;  //Z=10608
                }/*5*/  //Z=10609
            }/*4*/  //Z=10610
            goto Label99;  //Z=10611
        }/*3*/   /*  of homogeneous  */  //Z=10612

        /*  core/shell  */          /*  not yet ready  */  //Z=10614
        /* if cs=1 then begin  //Z=10615
           end;   (* of core/shell *) */  //Z=10645

    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=10646


    /*  perfect orientation case for ellipsoids  */  //Z=10648
    if ( (ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10649
        params.norm = 1;  //Z=10650
        order = 1;  //Z=10651
        /*  homogeneous  */  //Z=10652
        if ( params.cs==0 )
        {/*3*/  //Z=10653
            for ( n=1; n<=2*nmax+2; n++ )
            {/*4*/  //Z=10654
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10655
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10656
                fkv[n] = fkv[n-1]*n;  //Z=10657
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10658
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=10659
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10660
                xln[n] = -xln[n-1]*xl2z;  //Z=10661
            }/*4*/  //Z=10662
            for ( n=0; n<=nmax; n++ )
            {/*4*/  //Z=10663
                a1 = sqr(3/4.0);  //Z=10664
                for ( m=0; m<=nmax; m++ )
                {/*5*/  //Z=10665
                    sump1 = 0.0;  //Z=10666
                    for ( k=0; k<=n; k++ )
                    {/*6*/  //Z=10667
                        sump2 = 0.0;  //Z=10668
                        for ( ll=0; ll<=m; ll++ )  //Z=10669
                            sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10670
                        sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);  //Z=10671
                    }/*6*/  //Z=10672
                    params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10673
                }/*5*/  //Z=10674

                params.CR->carr4p[n] = a1*z12vl[n]*xln[n];  //Z=10676
                params.CR->carr5p[n] = z12v[n]*xrn[n];  //Z=10677

                /*  P(q)-coefficient  */  //Z=10679
                /* for m:=0 to n do begin  //Z=10680 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10681 */
                /* carr1fm[i]:=carr1pm[i];  //Z=10682 */
                /* i:=i+1;  //Z=10683 */
                /* end;  //Z=10684 */

                /*  F(q)-coefficient  */  //Z=10686
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=10687
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10688
                    if ( n<n4 ) n4 = n;  //Z=10689
                }/*5*/  //Z=10690
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10691
                    if ( n<n4f ) n4f = n;  //Z=10692
                }/*5*/  //Z=10693
            }/*4*/  //Z=10694
            goto Label99;  //Z=10695
        }/*3*/   /*  of homogeneous  */  //Z=10696

        /*  core/shell  */          /*  not yet ready  */  //Z=10698
        /* if cs=1 then begin  //Z=10699
            goto 99;  //Z=10728
           end;   (* of core/shell *)   */  //Z=10729
    }/*2*/  /*  of dim=5, ellipsoid  */  //Z=10730


    /* ** ellipsoid orientational distribution  */  //Z=10732
    if ( (ordis==0) && (params.orcase==2) && (dim==5) )
    {/*2*/  //Z=10733

        if ( params.cs==0 )
        {/*3*/  //Z=10735
            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,2,0,0,0,0,params.CR->carr1p,params.norm);  //Z=10736
            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,3,0,0,0,0,params.CR->carr1p,order);  //Z=10737
            order = order/params.norm;  //Z=10738

            for ( n=1; n<=2*nmax+2; n++ )
            {/*4*/  //Z=10740
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10741
                z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10742
                fkv[n] = fkv[n-1]*n;  //Z=10743
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10744
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10745
                xln[n] = -xln[n-1]*xl2z;  //Z=10746
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10747
            }/*4*/  //Z=10748

            for ( n=0; n<=nmax; n++ )
            {/*4*/  //Z=10750
                for ( m=0; m<=nmax; m++ )
                {/*5*/  //Z=10751
                    sump1 = 0.0;  //Z=10752
                    for ( k=0; k<=n; k++ )
                    {/*6*/  //Z=10753
                        sump2 = 0.0;  //Z=10754
                        for ( ll=0; ll<=m; ll++ )  //Z=10755
                            sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10756
                        sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);  //Z=10757
                    }/*6*/  //Z=10758
                    params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10759
                }/*5*/  //Z=10760
                sump1 = 0.0;  //Z=10761
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10762
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=10763
                    sump1 = sump1+pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=10764
                    params.CR->carr22pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=10765
                }/*5*/  //Z=10766

                /*  all coefficients: Mok  */  //Z=10768
                /*  integral part  */  //Z=10769
                params.CR->carr4p[n] = sump1;  //Z=10770
                /*  qR-part  */  //Z=10771
                params.CR->carr5p[n] = z12v[n]*xrn[n];  //Z=10772
                /*  qL-part  */  //Z=10773
                /* carr3p[n]:=sqr(3/4)*fk2v[n]*z12v[n]*xln[n]/power(4,n);  //Z=10774 */
                params.CR->carr6p[n] = sqr(3/4.0)*fk2v[n]*z12vl[n]*xln[n];  //Z=10775

                /* (* P(q)-coefficient *)  //Z=10778
                //carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=10779
                (* F(q)-coefficient *)  //Z=10780
                sump:=0.0;  //Z=10781
                for m:=0 to n do sump:=sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10782
                carr1f[n]:=fk2v[n]*xln[n]*sump/(power(4,n));  //Z=10783

                (* cross-section *)  //Z=10785
                carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=10786
                sump:=0.0;  //Z=10787
                for m:=0 to n do sump:=sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=10788
                carr4f[n]:=xrn[n]*sump;  //Z=10789
                */  //Z=10790

                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=10792
                    if ( n<n3 ) n3 = n;  //Z=10793
                }/*5*/  //Z=10794
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10795
                    if ( n<n4 ) n4 = n;  //Z=10796
                }/*5*/  //Z=10797
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=10798
                    if ( n<n5 ) n5 = n;  //Z=10799
                }/*5*/  //Z=10800
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=10801
                    if ( n<n6 ) n6 = n;  //Z=10802
                }/*5*/  //Z=10803
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=10804
                    if ( n<n1f ) n1f = n;  //Z=10805
                }/*5*/  //Z=10806
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10807
                    if ( n<n4f ) n4f = n;  //Z=10808
                }/*5*/  //Z=10809
            }/*4*/  /*  of n-loop  */  //Z=10810
            goto Label99;  //Z=10812
        }/*3*/  /*  of cs-loop  */  //Z=10813
    }/*2*/  /*  of ordis-loop  if ( (ordis==0) && (cho1==2) && (dim==5) ) */  //Z=10814


    /*  isotropic case for triaxial ellipsoids  */  //Z=10817
    if ( (ordis==7) && (dim==6) )
    {/*2*/    /*  triaxial ellipsoids  */  //Z=10818
        params.norm = 1;  //Z=10819
        order = 0;  //Z=10820
        if ( (aell==bell) && (bell==cell) )
            area = 4*M_PI*aell*aell;  //Z=10821
        else
            area = 4*M_PI*pow((pow(aell*bell,8/5.0)+pow(bell*cell,8/5.0)+pow(aell*cell,8/5.0))/3.0,5/8.0);  //Z=10822
        vol = (4*M_PI/3.0)*aell*bell*cell;  //Z=10823
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10824

        /*  homogeneous  */  //Z=10826
        if ( params.cs==0 )
        {/*3*/  //Z=10827
            u1ell[0] = 1;  //Z=10829
            u2ell[0] = 1;  //Z=10830
            v1ell[0] = 1;  //Z=10831
            v2ell[0] = 1;  //Z=10832
            gell[0] = sqrt(M_PI);  //Z=10833
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10834
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10835
                fkv[n] = fkv[n-1]*n;  //Z=10836
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10837
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10838
                u1ell[n] = (aell*aell-bell*bell)*u1ell[n-1]/n;  //Z=10839
                u2ell[n] = bell*bell*u2ell[n-1]/n;  //Z=10840
                v1ell[n] = (bell*bell-aell*aell)*v1ell[n-1]/n;  //Z=10841
                v2ell[n] = (cell*cell-bell*bell)*v2ell[n-1]/n;  //Z=10842
                gell[n] = gam3[n]/((n+1/2.0)*fkv[n]);  //Z=10843
            }/*4*/  //Z=10844

            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10846
                sump = 0.0;  //Z=10859
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10860
                    sump1 = 0.0;  //Z=10861
                    for ( k=0; k<=m; k++ )
                    {/*6*/  //Z=10862
                        sump2 = 0.0;  //Z=10863
                        for ( ll=0; ll<=n-m; ll++ ) sump2 = sump2+gell[k+ll]*v1ell[ll]*v2ell[n-m-ll];  //Z=10864
                        sump1 = sump1+u1ell[k]*u2ell[m-k]*sump2;  //Z=10865
                    }/*6*/  //Z=10866
                    sump = sump+sump1/((2.0*(n-m))+1);  //Z=10867
                }/*5*/  //Z=10868

                /* fsum[n]:=sumf;  //Z=10870 */
                sumf = 0.0;  //Z=10871
                for ( m=0; m<=n; m++ ) sumf = sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10872
                fsum[n] = sumf;  //Z=10873
                /*  P(q)-coefficient  */  //Z=10874
                params.CR->carr4p[n] = (9*M_PI)*pow(4.0,n-1)*z12v[n]*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);  //Z=10875
                /*  F(q)-coefficient  */  //Z=10876
                sump = 0.0;  //Z=10877
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10878
                    sump1 = 0.0;  //Z=10879
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10880
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10881
                }/*5*/  //Z=10882
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]*sump/(128.0*gam3[n]);  //Z=10883
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10884
                    if ( n<n4 ) n4 = n;  //Z=10885
                }/*5*/  //Z=10886
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10887
                    if ( n<n4f ) n4f = n;  //Z=10888
                }/*5*/  //Z=10889
            }/*4*/  //Z=10890
            goto Label99;  //Z=10891
        }/*3*/   /*  of homogeneous  */  //Z=10892

        /*  core/shell  */          /*  not yet ready  */  //Z=10894
        /* if cs=1 then begin  //Z=10895
             goto 99;  //Z=10924
           end;   (* of core/shell *) */  //Z=10925
    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=10926


    /*  perfect orientation case for ellipsoids  */  //Z=10928  TODO: siehe unten | cs==0
    if ( (ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10929
        params.norm = 1;  //Z=10930
        order = 1;  //Z=10931
        /*  homogeneous  */  //Z=10932
        if ( params.cs==0 )
        {/*3*/  //Z=10933
            i = 2;  //Z=10934
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10935
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10936
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10937
                fkv[n] = fkv[n-1]*n;  //Z=10938
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10939
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=10940
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10941
                xln[n] = -xln[n-1]*xl2z;  //Z=10942
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=10943
                    sump = 0.0;  //Z=10944
                    for ( ns=0; ns<=n; ns++ )
                    {/*6*/  //Z=10945
                        sump1 = 0.0;  //Z=10946
                        for ( ms=0; ms<=m; ms++ ) sump1 = sump1+1/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);  //Z=10947
                        sump = sump+sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);  //Z=10948
                    }/*6*/  //Z=10949
                    params.CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;  //Z=10950
                }/*5*/  //Z=10951
                /*  P(q)-coefficient  */  //Z=10952
                /* for m:=0 to n do begin  //Z=10953 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10954 */
                /* carr1fm[i]:=carr1pm[i];  //Z=10955 */
                /* i:=i+1;  //Z=10956 */
                /* end;  //Z=10957 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn[n]/(16.0*gam3[n]);  //Z=10958
                /*  F(q)-coefficient  */  //Z=10959
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=10960
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10961
                    if ( n<n4 ) n4 = n;  //Z=10962
                }/*5*/  //Z=10963
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10964
                    if ( n<n4f ) n4f = n;  //Z=10965
                }/*5*/  //Z=10966
            }/*4*/  //Z=10967
            goto Label99;  //Z=10968
        }/*3*/   /*  of homogeneous  */  //Z=10969

        /*  core/shell  */          /*  not yet ready  */  //Z=10971
        /* if cs=1 then begin  //Z=10972
             goto 99;  //Z=11001
           end;   (* of core/shell *)   */  //Z=11002
    }/*2*/  /*  of dim=6, triellipsoid  */  //Z=11003


    /*  isotropic case for super ellipsoids, barrel  */  //Z=11006
    if ( (ordis==7) && (dim==7) )
    {/*2*/    /*  super ellipsoids  */  //Z=11007

        // TODO
        double qx=1, qy=1, qz=1, qhkl=1;
        // Diese Variablen sind im Orginal-Pascalprogramm an der Stelle des coefficients() Aufrufes nicht bekannt.

        params.norm = 1;  //Z=11008
        order = 0;  //Z=11009
        qrombchid(params.length,params.radius,params.alphash1,params.sigma,params.alphash1,params.polPhi,params.polTheta,params.polPhi,qx,qy,qz,params.p11,params.p12,params.p13,params.p21,params.p22,params.p23,params.p31,params.p32,params.p33,qx,qy,0,qhkl,
                  params.ax1.length(),params.ax2.length(),params.ax3.length(),
                  params.ax1.x(),params.ax1.y(),params.ax1.z(),
                  params.ax2.x(),params.ax2.y(),params.ax2.z(),
                  params.ax3.x(),params.ax3.y(),params.ax3.z(),
                  params.sig.x(),params.sig.y(),params.sig.z(),
                  ordis,3,8,15,7,0,0,params.CR->carr1p,area);   //Z=11010
        area = 2*M_PI*area;  //Z=11011
        vol = 2*M_PI*sqr(params.radius)*params.length*gamma((2+params.alphash1)/params.alphash1)*gamma(1/params.alphash1)/(params.alphash1*gamma((3+params.alphash1)/params.alphash1));  //Z=11012
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=11013

        /*  homogeneous  */  //Z=11015
        u1ell[0] = 1;             /*  for barrel  */  //Z=11016
        u2ell[0] = 1;             /*  for barrel  */  //Z=11017
        for ( n=1; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11018
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11019
            fkv[n] = fkv[n-1]*n;  //Z=11020
            gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11021
            u1ell[n] = u1ell[n-1]/((n-1/2.0)*n);   /*  for barrel  */  //Z=11022
            u2ell[n] = u2ell[n-1]/((n+1)*n);     /*  for barrel  */  //Z=11023
        }/*3*/  //Z=11024

        for ( n=0; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11026
            v1ell[n] = gamma((2*n+1)/params.alphash1)*u1ell[n];         /*  for barrel  */  //Z=11027
            v2ell[n] = gamma((2*n+2+params.alphash1)/params.alphash1)*u2ell[n];    /*  for barrel  */  //Z=11028
            gell[n] = gamma((2*n+3+params.alphash1)/params.alphash1);              /*  for barrel  */  //Z=11029
        }/*3*/  //Z=11030

        if ( params.cs==0 )
        {/*3*/  //Z=11032
            for ( n=0; n<=nmax; n++ )
            {/*4*/  //Z=11033

                if ( params.alphash1==2 )
                {/*5*/  //Z=11047
                    a1 = sqr(3/4.0);  //Z=11048
                    for ( m=0; m<=nmax; m++ )
                    {/*6*/  //Z=11049
                        sump1 = 0.0;  //Z=11050
                        for ( k=0; k<=n; k++ )
                        {/*7*/  //Z=11051
                            sump2 = 0.0;  //Z=11052
                            for ( ll=0; ll<=m; ll++ )  //Z=11053
                                sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=11054
                            sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);  //Z=11055
                        }/*7*/  //Z=11056
                        params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=11057
                    }/*6*/  //Z=11058
                }/*5*/  //Z=11059
                else
                {/*5*/  //Z=11060
                    a1 = gamma((params.alphash1+3)/params.alphash1)/(gamma((params.alphash1+2)/params.alphash1)*gamma(1/params.alphash1));  //Z=11061
                    a1 = a1*a1;  //Z=11062
                    for ( m=0; m<=nmax; m++ )
                    {/*6*/  //Z=11063
                        sump1 = 0.0;  //Z=11064
                        for ( ns=0; ns<=n; ns++ )
                        {/*7*/  //Z=11065
                            sump2 = 0.0;  //Z=11066
                            for ( ms=0; ms<=m; ms++ )  //Z=11067
                                sump2 = sump2+v2ell[m-ms]*v2ell[ms]/(gell[ns+ms]*gell[n-ns+m-ms]);  //Z=11068
                            sump1 = sump1+v1ell[n-ns]*v1ell[ns]*sump2;  //Z=11069
                        }/*7*/  //Z=11070
                        params.CR->carr11pm[n][m] = sump1;  //Z=11071
                    }/*6*/  //Z=11072
                }/*5*/  //Z=11073

                /*  orientational average  */  //Z=11075
                sump = 0.0;  //Z=11076
                for ( m=0; m<=n; m++ )
                    sump += gam3[n-m]*z12v[n-m]*z12v[m]*pow(sqr(params.length),n-m)*pow(sqr(params.radius),m)*fkv[m]*params.CR->carr11pm[n-m][m]/(n-m+1/2.0);  //Z=11077
                params.CR->carr4p[n] = a1*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/(2.0*gam3[n]);  //Z=11078

                /* fsum[n]:=sumf;  //Z=11081 */
                /* sumf:=0.0;  //Z=11082 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11083 */
                /* fsum[n]:=sumf;  //Z=11084 */
                /*  P(q)-coefficient  */  //Z=11085
                /*  F(q)-coefficient  */  //Z=11086
                sump = 0.0;  //Z=11087
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=11088
                    sump1 = 0.0;  //Z=11089
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=11090
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=11091
                }/*5*/  //Z=11092
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]*sump/(128.0*gam3[n]);  //Z=11093
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11094
                    if ( n<n4 ) n4 = n;  //Z=11095
                }/*5*/  //Z=11096
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11097
                    if ( n<n4f ) n4f = n;  //Z=11098
                }/*5*/  //Z=11099
            }/*4*/  //Z=11100
            goto Label99;  //Z=11101
        }/*3*/   /*  of homogeneous  */  //Z=11102

        /*  core/shell  */          /*  not yet ready  */  //Z=11104
        /* if cs=1 then begin  //Z=11105
             goto 99;  //Z=11134
           end;   (* of core/shell *) */  //Z=11135
    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=11136


    /*  perfect orientation case for ellipsoids  */  //Z=11138  TODO: gleich wie oben | cs==0  -  in neuer Version auskommentiert
    if ( (ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=11139
#ifdef undef
        params.norm = 1;  //Z=11140
        order = 1;  //Z=11141
        /*  homogeneous  */  //Z=11142
        if ( params.cs==0 )
        {/*3*/  //Z=11143
            i = 2;  //Z=11144
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11145
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11146
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11147
                fkv[n] = fkv[n-1]*n;  //Z=11148
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11149
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=11150
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11151
                xln[n] = -xln[n-1]*xl2z;  //Z=11152
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=11153
                    sump = 0.0;  //Z=11154
                    for ( ns=0; ns<=n; ns++ )
                    {/*6*/  //Z=11155
                        sump1 = 0.0;  //Z=11156
                        for ( ms=0; ms<=m; ms++ ) sump1 = sump1+1/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);  //Z=11157
                        sump = sump+sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);  //Z=11158
                    }/*6*/  //Z=11159
                    params.CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;  //Z=11160
                }/*5*/  //Z=11161
                /*  P(q)-coefficient  */  //Z=11162
                /* for m:=0 to n do begin  //Z=11163 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=11164 */
                /* carr1fm[i]:=carr1pm[i];  //Z=11165 */
                /* i:=i+1;  //Z=11166 */
                /* end;  //Z=11167 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn[n]/(16.0*gam3[n]);  //Z=11168
                /*  F(q)-coefficient  */  //Z=11169
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=11170
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11171
                    if ( n<n4 ) n4 = n;  //Z=11172
                }/*5*/  //Z=11173
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11174
                    if ( n<n4f ) n4f = n;  //Z=11175
                }/*5*/  //Z=11176
            }/*4*/  //Z=11177
            goto Label99;  //Z=11178
        }/*3*/   /*  of homogeneous  */  //Z=11179

        /*  core/shell  */          /*  not yet ready  */  //Z=11181
        /* if cs=1 then begin  //Z=11182
             goto 99;  //Z=11211
           end;   (* of core/shell *)   */  //Z=11212
#endif
    }/*2*/  /*  of dim=7, super ellipsoid  */  //Z=11213


    /*  isotropic case for superballs  */  //Z=11217
    if ( (ordis==7) && (dim==8) )
    {/*2*/    /*  superball  */  //Z=11218

        //type ArrayImax3D=array[0..50,0..50,0..50] of real;
        //carr111pm: ^ArrayImax3D;
        //new(carr111pm);   /*Z=10728*/
        static const int ardim = 50;  //Z=11221
        //SetLength(carr111pm^, (ardim+1));  //Z=11222
        //for ( i=0; i<=ardim; i++ ) SetLength(carr111pm^[i], (ardim+1));  //Z=11223
        //for ( i=0; i<=ardim; i++ )
        //{/*3*/  //Z=11224
        //    for ( j=0; j<=ardim; j++ ) SetLength(carr111pm^[i][j], (ardim+1));  //Z=11225
        //}/*3*/  //Z=11226

        float carr111pm[ardim+1][ardim+1][ardim+1];

        nmax = 40;  //Z=11231
        params.norm = 1;  //Z=11232
        order = 0;  //Z=11233
        /* l:=r;  //Z=11234 */

        /*  radius=a, rm=b, length=c  */  //Z=11236
        qrombdeltac(params.radiusi,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,8,17,0,0,0,params.CR->carr1p,area);  //Z=11237
        area = 8*area;  //Z=11238
        vol = 8*params.radius*params.radiusi*params.length*pow(gamma(1/params.alphash1),3)/(pow(params.alphash1,3)*gamma(1+3/params.alphash1));  //Z=11239
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=11240

        /*  homogeneous  */  //Z=11243
        u3ell[0] = 1;             /*  for superball  */  //Z=11244
        for ( n=1; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11245
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11246
            fkv[n] = fkv[n-1]*n;  //Z=11247
            gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11248
            u3ell[n] = u3ell[n-1]/((n-1/2.0)*n);   /*  for superball  */  //Z=11249
        }/*3*/  //Z=11250

        /* for n:=0 to 3*nmax+1 do begin  //Z=11252 */
        /*    v3ell[n]:=gamma((2*n+1)/alfa)*u3ell[n];         (* for superball *)  //Z=11253 */
        /*    g3ell[n]:=gamma(((2*n+3)/alfa)+1);              (* for superball *)  //Z=11254 */
        /* end;  //Z=11255 */

        for ( n=0; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11257
            v3ell[n] = exp(gammln((2*n+1)/params.alphash1))*u3ell[n];         /*  for superball  */  //Z=11258
            g3ell[n] = exp(gammln(((2*n+3)/params.alphash1)+1));              /*  for superball  */  //Z=11259
        }/*3*/  //Z=11260

        if ( params.cs==0 )
        {/*3*/  //Z=11262
            for ( n=0; n<=nmax; n++ )
            {/*4*/  //Z=11263
                if ( params.alphash1==200000 )
                {/*5*/  //Z=11264
                    a1 = sqr(3/4.0);  //Z=11265
                    for ( m=0; m<=nmax; m++ )
                    {/*6*/  //Z=11266
                        sump1 = 0.0;  //Z=11267
                        for ( k=0; k<=n; k++ )
                        {/*7*/  //Z=11268
                            sump2 = 0.0;  //Z=11269
                            for ( ll=0; ll<=m; ll++ )  //Z=11270
                                sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=11271
                            sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);  //Z=11272
                        }/*7*/  //Z=11273
                        params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=11274
                    }/*6*/  //Z=11275
                }/*5*/  //Z=11276
                else
                {/*5*/  //Z=11277
                    a1 = gamma(1+3/params.alphash1)/pow(gamma(1/params.alphash1),3);  //Z=11278
                    a1 = a1*a1;  //Z=11279
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11280
                        for ( k=0; k<=nmax; k++ )
                        {/*7*/  //Z=11281
                            sump1 = 0.0;  //Z=11282
                            for ( ns=0; ns<=n; ns++ )
                            {/*8*/  //Z=11283
                                sump2 = 0.0;  //Z=11284
                                for ( ms=0; ms<=m; ms++ )
                                {/*9*/  //Z=11285
                                    double sump3 = 0.0;  //Z=11286
                                    for ( int ks=0; ks<=k; ks++ )  //Z=11287
                                        sump3 = sump3+v3ell[k-ks]*v3ell[ks]/(g3ell[ns+ms+ks]*g3ell[n-ns+m-ms+k-ks]);  //Z=11288
                                    sump2 = sump2+v3ell[m-ms]*v3ell[ms]*sump3;  //Z=11289
                                }/*9*/  //Z=11290
                                sump1 = sump1+v3ell[n-ns]*v3ell[ns]*sump2;  //Z=11291
                            }/*8*/  //Z=11292
                            carr111pm[n][m][k] = sump1;  //Z=11293
                            carr111pm[m][n][k] = sump1;  //Z=11294
                        }/*7*/  //Z=11295
                    }/*6*/  //Z=11296
                }/*5*/  //Z=11297

                /*  orientational average for superball  */  //Z=11299
                sump = 0.0;  //Z=11300
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=11301
                    sump1 = 0.0;  //Z=11302
                    for ( k=0; k<=m; k++ )  //Z=11303
                        sump1 += z12v[m-k]*z12v[k]*gam3[m-k]*gam3[k]*pow(sqr(params.radiusi),m-k)*pow(sqr(params.length),k)*carr111pm[n-m][m-k][k]/((m-k+1/2.0)*(k+1/2.0));   //Z=11304
                    sump = sump+z12v[n-m]*gam3[n-m]*pow(sqr(params.radius),n-m)*sump1/(n-m+1/2.0);  //Z=11305
                }/*5*/  //Z=11306
                params.CR->carr4p[n] = (a1/(2.0*M_PI))*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/gam3[n];  //Z=11307


                /* fsum[n]:=sumf;  //Z=11310 */
                /* sumf:=0.0;  //Z=11311 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11312 */
                /* fsum[n]:=sumf;  //Z=11313 */
                /*  P(q)-coefficient  */  //Z=11314
                /*  F(q)-coefficient  */  //Z=11315
                sump = 0.0;  //Z=11316
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=11317
                    sump1 = 0.0;  //Z=11318
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=11319
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=11320
                }/*5*/  //Z=11321
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]*sump/(128.0*gam3[n]);  //Z=11322
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11323
                    if ( n<n4 ) n4 = n;  //Z=11324
                }/*5*/  //Z=11325
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11326
                    if ( n<n4f ) n4f = n;  //Z=11327
                }/*5*/  //Z=11328
            }/*4*/  //Z=11329
            goto Label99;  //Z=11330
        }/*3*/   /*  of homogeneous  */  //Z=11331

        /*  core/shell  */          /*  not yet ready  */  //Z=11333
        /* if cs=1 then begin  //Z=11334
             goto 99;  //Z=11363
           end;   (* of core/shell *) */  //Z=11364
        //dispose(carr111pm);  //Z=11365
    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=11366


    /*  perfect orientation case for ellipsoids  */  //Z=11368  TODO: gleich wie oben | cs==0  -  in neuer Version auskommentiert
    if ( (ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=11369
#ifdef undef
        params.norm = 1;  //Z=11370
        order = 1;  //Z=11371
        /*  homogeneous  */  //Z=11372
        if ( params.cs==0 )
        {/*3*/  //Z=11373
            i = 2;  //Z=11374
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11375
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11376
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11377
                fkv[n] = fkv[n-1]*n;  //Z=11378
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11379
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=11380
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11381
                xln[n] = -xln[n-1]*xl2z;  //Z=11382
                for ( m=0; m<=n; m++ )
                {/*5*/  //Z=11383
                    sump = 0.0;  //Z=11384
                    for ( ns=0; ns<=n; ns++ )
                    {/*6*/  //Z=11385
                        sump1 = 0.0;  //Z=11386
                        for ( ms=0; ms<=m; ms++ ) sump1 = sump1+1/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);  //Z=11387
                        sump = sump+sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);  //Z=11388
                    }/*6*/  //Z=11389
                    params.CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;  //Z=11390
                }/*5*/  //Z=11391
                /*  P(q)-coefficient  */  //Z=11392
                /* for m:=0 to n do begin  //Z=11393 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=11394 */
                /* carr1fm[i]:=carr1pm[i];  //Z=11395 */
                /* i:=i+1;  //Z=11396 */
                /* end;  //Z=11397 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn[n]/(16.0*gam3[n]);  //Z=11398
                /*  F(q)-coefficient  */  //Z=11399
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=11400
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11401
                    if ( n<n4 ) n4 = n;  //Z=11402
                }/*5*/  //Z=11403
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11404
                    if ( n<n4f ) n4f = n;  //Z=11405
                }/*5*/  //Z=11406
            }/*4*/  //Z=11407
            goto Label99;  //Z=11408
        }/*3*/   /*  of homogeneous  */  //Z=11409

        /*  core/shell  */          /*  not yet ready  */  //Z=11411
        /* if cs=1 then begin  //Z=11412
             goto 99;  //Z=11441
           end;   (* of core/shell *)   */  //Z=11442
#endif
    }/*2*/  /*  of dim=8, superball  */  //Z=11443


    /* ** isotropic case for cylinders and disks ** */  //Z=11447
    if ( (ordis==7) && (dim!=3) )
    {/*2*/  //Z=11448
        params.norm = 1;  //Z=11449
        order = 0;  //Z=11450
        /*  homogeneous  */  //Z=11451
        if ( params.cs==0 )
        {/*3*/  //Z=11452

            /*  small axial ratios  */  //Z=11454
            if ( (params.length/params.radius)<2 )
            {/*4*/  //Z=11455
                area = 2*M_PI*sqr(params.radius)+2*M_PI*params.radius*(2*params.length);  //Z=11456
                vol = M_PI*sqr(params.radius)*(2*params.length);  //Z=11457
                params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=11458

                for ( n=1; n<=nmax; n++ )
                {/*5*/  //Z=11460
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11461
                    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11462
                    b1sv_n = b1sv_n*(b1s-1+n);  //Z=11463
                    fkv[n] = fkv[n-1]*n;  //Z=11464
                    fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11465
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11466
                    /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11467 */
                    xln[n] = -xln[n-1]*xl2z;  //Z=11468
                    xrn[n] = -xrn[n-1]*xr2z;  //Z=11469
                }/*5*/  //Z=11470

                for ( n=0; n<=nmax; n++ )
                {/*5*/  //Z=11472
                    binsum = 0.0;  //Z=11473
                    /* for m:=0 to n do    (* Cauchy sum *)  //Z=11474 */
                    /*       binsum:=binsum+gam3[m]*z12vl[n-m]*z12v[m]*xln[n-m]*xrn[m]/((n-m+1/2)*(n-m+1)*fkv[n-m]*(m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]);  //Z=11475 */
                    params.CR->carr11pm[n][0] = sqrt(M_PI)/gam3[n];  //Z=11476
                    a1 = sqrt(M_PI)/(2.0*gam3[n]);  //Z=11477
                    for ( m=1; m<=nmax; m++ )
                    {/*6*/    /*  double sum  */  //Z=11478
                        a1 = a1*(m+1/2.0)/(n+m+1/2.0);  //Z=11479
                        /* carr11pm[n,m]:=power(4,m+1)*gam3[m]*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]*gam3[n+m]);   (* ok *)  //Z=11480 */
                        params.CR->carr11pm[n][m] = pow(4.0,m+1)*a1*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]);       /*  Mok  */  //Z=11481
                    }/*6*/  //Z=11482
                    params.CR->carr2p[n] = pow(4.0,n-1)*z12vl[n]*xln[n]/((n+1/2.0)*(n+1)*fkv[n]);     /*  double sum  */  //Z=11483
                    params.CR->carr3p[n] = pow(4.0,n)*binsum/gam3[n];      /*  Cauchy sum  */  //Z=11484
                }/*5*/  //Z=11485

            }/*4*/  //Z=11487
            /*  large axial ratios  */  //Z=11488
            else
            {/*4*/  //Z=11489

                for ( n=1; n<=nmax; n++ )
                {/*5*/  //Z=11491
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11492
                    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11493
                    b1sv_n = b1sv_n*(b1s-1+n);  //Z=11494
                    fkv[n] = fkv[n-1]*n;  //Z=11495
                    fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11496
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11497
                    /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11498 */
                    xln[n] = -xln[n-1]*xl2z;  //Z=11499
                    xrn[n] = -xrn[n-1]*xr2z;  //Z=11500
                    /*  cylinder, ok */  //Z=11501
                    if ( dim==1 )
                    {/*6*/  //Z=11502

                        /*  P(q) factorization  */  //Z=11504
                        params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);      /*  P||iso(q)  */  //Z=11505
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*  P-(q)  */  //Z=11506
                        /*  F(q)  */  //Z=11507
                        binsum = 0.0;  //Z=11508
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11509
                        params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11510
                        binsum = 0.0;  //Z=11511
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11512
                        params.CR->carr4f[n] = xrn[n]*binsum;  //Z=11513
                    }/*6*/  //Z=11514
                    /*  disk, ok  */  //Z=11515
                    if ( dim==2 )
                    {/*6*/  //Z=11516
                        /*  P(q)  */  //Z=11517
                        params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11518
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11519
                        /*  F(q)  */  //Z=11520
                        binsum = 0.0;  //Z=11521
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11522
                        params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11523
                        binsum = 0.0;  //Z=11524
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11525
                        params.CR->carr4f[n] = M_PI*xrn[n]*binsum/4.0;  //Z=11526
                    }/*6*/  //Z=11527
                }/*5*/  /*  of large axial ratios  */  //Z=11528

                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11530
                    if ( n<n1 ) n1 = n;  //Z=11531
                }/*5*/  //Z=11532
                if ( fabs(params.CR->carr2p[n])<min )
                {/*5*/  //Z=11533
                    if ( n<n2 ) n2 = n;  //Z=11534
                }/*5*/  //Z=11535
                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=11536
                    if ( n<n3 ) n3 = n;  //Z=11537
                }/*5*/  //Z=11538
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11539
                    if ( n<n1f ) n1f = n;  //Z=11540
                }/*5*/  //Z=11541
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11542
                    if ( n<n4 ) n4 = n;  //Z=11543
                }/*5*/  //Z=11544
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11545
                    if ( n<n4f ) n4f = n;  //Z=11546
                }/*5*/  //Z=11547
            }/*4*/  //Z=11548
        }/*3*/ /*  of homogeneous  */  //Z=11549

        /*  core/shell  */  //Z=11551
        if ( params.cs==1 )
        {/*3*/  //Z=11552
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11553
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11554
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11555
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=11556
                fkv[n] = fkv[n-1]*n;  //Z=11557
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11558
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11559
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11560 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11561
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11562
                xrmn_n = -xrmn_n*xrm2z;  //Z=11563
                pn[n] = pn[n-1]*p*p;  //Z=11564
                /* ** cylinder ** */  //Z=11565
                if ( dim==1 )
                {/*5*/  //Z=11566
                    /*  longitudinal P(q)  */  //Z=11567
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11568
                    /*  cross-sectional P(q)  */  //Z=11569
                    /*  F121  */  //Z=11570
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11571
                    /*  F122  */  //Z=11572
                    sump = 0.0;  //Z=11573
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11574
                        sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11575
                    }/*6*/  //Z=11576
                    params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=11577
                    /*  F123  */  //Z=11578
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11579

                    /*  longitudinal F(q)  */  //Z=11581
                    binsum = 0.0;  //Z=11582
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11583
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11584
                    /*  cross-sectional F(q)  */  //Z=11585
                    /*  F121  */  //Z=11586
                    sump = 0.0;  //Z=11587
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11588
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11589
                    }/*6*/  //Z=11590
                    params.CR->carr4f[n] = xrn[n]*sump;  //Z=11591
                    /*  F122  */  //Z=11592
                    sump = 0.0;  //Z=11593
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11594
                        sump = sump+z12v[m]*z12v[n-m]*pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11595
                    }/*6*/  //Z=11596
                    params.CR->carr5f[n] = xrmn_n*sump;  //Z=11597
                    /*  F123  */  //Z=11598
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=11599
                }/*5*/  //Z=11600

                /* ** disk ** */  //Z=11602
                if ( dim==2 )
                {/*5*/  //Z=11603
                    /*  longitudinal  */  //Z=11604
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11605
                    /*  cross-sectional P(q)  */  //Z=11606
                    /*  F121  */  //Z=11607
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11608
                    /*  F122  */  //Z=11609
                    sump = 0.0;  //Z=11610
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11611
                        sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11612
                    }/*6*/  //Z=11613
                    params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=11614
                    /*  F123  */  //Z=11615
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11616

                    /*  longitudinal F(q)  */  //Z=11618
                    binsum = 0.0;  //Z=11619
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11620
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11621
                    /*  cross-sectional F(q)  */  //Z=11622
                    /*  F121  */  //Z=11623
                    sump = 0.0;  //Z=11624
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11625
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11626
                    }/*6*/  //Z=11627
                    params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=11628
                    /*  F122  */  //Z=11629
                    sump = 0.0;  //Z=11630
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11631
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11632
                    }/*6*/  //Z=11633
                    params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=11634
                    /*  F123  */  //Z=11635
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=11636
                }/*5*/  //Z=11637
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11638
                    if ( n<n1 ) n1 = n;  //Z=11639
                }/*5*/  //Z=11640
                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=11641
                    if ( n<n3 ) n3 = n;  //Z=11642
                }/*5*/  //Z=11643
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11644
                    if ( n<n4 ) n4 = n;  //Z=11645
                }/*5*/  //Z=11646
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=11647
                    if ( n<n5 ) n5 = n;  //Z=11648
                }/*5*/  //Z=11649
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=11650
                    if ( n<n6 ) n6 = n;  //Z=11651
                }/*5*/  //Z=11652
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11653
                    if ( n<n1f ) n1f = n;  //Z=11654
                }/*5*/  //Z=11655
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11656
                    if ( n<n4f ) n4f = n;  //Z=11657
                }/*5*/  //Z=11658
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=11659
                    if ( n<n5f ) n5f = n;  //Z=11660
                }/*5*/  //Z=11661
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=11662
                    if ( n<n6f ) n6f = n;  //Z=11663
                }/*5*/  //Z=11664
            }/*4*/  //Z=11665
        }/*3*/ /*  of core/shell  */  //Z=11666

        /*  inhomogeneous core/shell  */  //Z=11668
        if ( params.cs==2 )
        {/*3*/  //Z=11669
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11670
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11671
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11672
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=11673
                fkv[n] = fkv[n-1]*n;  //Z=11674
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11675
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11676
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11677 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11678
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11679
                xrmn_n = -xrmn_n*xrm2z;  //Z=11680
                pn[n] = pn[n-1]*p*p;  //Z=11681
                /* ** cylinder ** */  //Z=11682
                if ( dim==1 )
                {/*5*/  //Z=11683
                    /*  longitudinal P(q)  */  //Z=11684
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11685

                    /*  cross-sectional P(q)  */  //Z=11687
                    params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11688
                    sump = 0.0;  //Z=11689
                    sump1 = 0.0;  //Z=11690
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11691
                        sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=11692
                        sump = sump+pn[n-m]*sumi;  //Z=11693
                        sump1 = sump1+sumi;  //Z=11694
                    }/*6*/  //Z=11695
                    params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=11696
                    params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n]*sump1;  //Z=11697
                    sump = 0.0;  //Z=11698
                    sump1 = 0.0;  //Z=11699
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11700
                        sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=11701
                        sump = sump+sumi;  //Z=11702
                        sump1 = sump1+pn[n-m]*sumi;  //Z=11703
                    }/*6*/  //Z=11704
                    params.CR->carr7p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=11705
                    params.CR->carr8p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=11706
                    params.CR->carr9p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=11707


                    /* (* cross-sectional P(q) *)  //Z=11710
                       (* F121 *)  //Z=11711
                       carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11712
                       (* F122 *)  //Z=11713
                       sump:=0.0;  //Z=11714
                          for m:=0 to n do begin  //Z=11715
                          sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11716
                       end;  //Z=11717
                       carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=11718
                       (* F123 *)  //Z=11719
                       carr6p[n]:=carr4p[n]/pn[n];  */  //Z=11720

                    /*  longitudinal F(q)  */  //Z=11722
                    binsum = 0.0;  //Z=11723
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11724
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11725
                    /*  cross-sectional F(q)  */  //Z=11726
                    params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=11727
                    params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11728
                    params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11729

                }/*5*/  //Z=11731

                /* ** disk ** */  //Z=11733
                if ( dim==2 )
                {/*5*/  //Z=11734
                    /*  longitudinal  */  //Z=11735
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11736


                    /*  cross-sectional P(q)  */  //Z=11739
                    params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=11740
                    sump = 0.0;  //Z=11741
                    sump1 = 0.0;  //Z=11742
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11743
                        sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11744
                        sump = sump+pn[n-m]*sumi;  //Z=11745
                        sump1 = sump1+sumi;  //Z=11746
                    }/*6*/  //Z=11747
                    params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=11748
                    params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=11749
                    sump = 0.0;  //Z=11750
                    sump1 = 0.0;  //Z=11751
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11752
                        sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=11753
                        sump = sump+sumi;  //Z=11754
                        sump1 = sump1+pn[n-m]*sumi;  //Z=11755
                    }/*6*/  //Z=11756
                    params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;   //Z=11757
                    params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=11758
                    params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;   //Z=11759


                    /*  cross-sectional P(q)  */  //Z=11762
                    /* (* F121 *)  //Z=11763
                       carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=11764
                       (* F122 *)  //Z=11765
                       sump:=0.0;  //Z=11766
                          for m:=0 to n do begin  //Z=11767
                          sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11768
                       end;  //Z=11769
                       carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=11770
                       (* F123 *)  //Z=11771
                       carr6p[n]:=carr4p[n]/pn[n];   */  //Z=11772

                    /*  longitudinal F(q)  */  //Z=11774
                    binsum = 0.0;  //Z=11775
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11776
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11777
                    /*  cross-sectional F(q)  */  //Z=11778
                    params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=11779
                    params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11780
                    params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11781
                }/*5*/  //Z=11782
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11783
                    if ( n<n1 ) n1 = n;  //Z=11784
                }/*5*/  //Z=11785
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11786
                    if ( n<n4 ) n4 = n;  //Z=11787
                }/*5*/  //Z=11788
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=11789
                    if ( n<n5 ) n5 = n;  //Z=11790
                }/*5*/  //Z=11791
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=11792
                    if ( n<n6 ) n6 = n;  //Z=11793
                }/*5*/  //Z=11794
                if ( fabs(params.CR->carr7p[n])<min )
                {/*5*/  //Z=11795
                    if ( n<n7 ) n7 = n;  //Z=11796
                }/*5*/  //Z=11797
                if ( fabs(params.CR->carr8p[n])<min )
                {/*5*/  //Z=11798
                    if ( n<n8 ) n8 = n;  //Z=11799
                }/*5*/  //Z=11800
                if ( fabs(params.CR->carr9p[n])<min )
                {/*5*/  //Z=11801
                    if ( n<n9 ) n9 = n;  //Z=11802
                }/*5*/  //Z=11803
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11804
                    if ( n<n1f ) n1f = n;  //Z=11805
                }/*5*/  //Z=11806
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11807
                    if ( n<n4f ) n4f = n;  //Z=11808
                }/*5*/  //Z=11809
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=11810
                    if ( n<n5f ) n5f = n;  //Z=11811
                }/*5*/  //Z=11812
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=11813
                    if ( n<n6f ) n6f = n;  //Z=11814
                }/*5*/  //Z=11815
                if ( fabs(params.CR->carr7f[n])<min )
                {/*5*/  //Z=11816
                    if ( n<n7f ) n7f = n;  //Z=11817
                }/*5*/  //Z=11818
                if ( fabs(params.CR->carr8f[n])<min )
                {/*5*/  //Z=11819
                    if ( n<n8f ) n8f = n;  //Z=11820
                }/*5*/  //Z=11821
                if ( fabs(params.CR->carr9f[n])<min )
                {/*5*/  //Z=11822
                    if ( n<n9f ) n9f = n;  //Z=11823
                }/*5*/  //Z=11824
            }/*4*/  //Z=11825
        }/*3*/ /*  of inhomogeneous core/shell  */  //Z=11826


        /*  myelin  */  //Z=11829
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=11830
            i = 2;  //Z=11831
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11832
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11833
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11834
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=11835
                fkv[n] = fkv[n-1]*n;  //Z=11836
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11837
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11838
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11839 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11840
                /* xrn[n]:=-xrn[n-1]*xr2z;  //Z=11841 */
                xrn[n] = -xrn[n-1]*x12zm;         /*  myelin radius  */  //Z=11842
                /*  cylinder, ok */  //Z=11843
                if ( dim==1 )
                {/*5*/  //Z=11844
                    /*  P(q)  */  //Z=11845
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11846

                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11848
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11849 */
                        /* i:=i+1;  //Z=11850 */
                        params.CR->carr11pm[n][m] = 1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11851
                    }/*6*/  //Z=11852
                    params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=11853
                    /* carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11854 */


                    /*  F(q)  */  //Z=11857
                    binsum = 0.0;  //Z=11858
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11859
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11860
                    binsum = 0.0;  //Z=11861
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11862
                    params.CR->carr4f[n] = xrn[n]*binsum;  //Z=11863
                }/*5*/  //Z=11864
                /*  disk, ok  */  //Z=11865
                if ( dim==2 )
                {/*5*/  //Z=11866
                    /*  P(q)  */  //Z=11867
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11868
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11869
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11870 */
                        /* i:=i+1;  //Z=11871 */
                        params.CR->carr11pm[n][m] = (M_PI/4.0)*(1/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]));  //Z=11872
                    }/*6*/  //Z=11873
                    params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=11874

                    /* carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=11876 */
                    /*  F(q)  */  //Z=11877
                    binsum = 0.0;  //Z=11878
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11879
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11880
                    binsum = 0.0;  //Z=11881
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11882
                    params.CR->carr4f[n] = M_PI*xrn[n]*binsum/4.0;  //Z=11883
                }/*5*/  //Z=11884
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11885
                    if ( n<n1 ) n1 = n;  //Z=11886
                }/*5*/  //Z=11887
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11888
                    if ( n<n1f ) n1f = n;  //Z=11889
                }/*5*/  //Z=11890
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11891
                    if ( n<n4 ) n4 = n;  //Z=11892
                }/*5*/  //Z=11893
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11894
                    if ( n<n4f ) n4f = n;  //Z=11895
                }/*5*/  //Z=11896
            }/*4*/  //Z=11897
        }/*3*/ /*  of myelin  */  //Z=11898

    }/*2*/  //Z=11900


    /* ** perfect orientation for cylinders and disks ** */  //Z=11902
    if ( (ordis==6) && (dim!=3) )
    {/*2*/  //Z=11903
        params.norm = 1;  //Z=11904
        order = 1;  //Z=11905
        /*  homogeneous  */  //Z=11906
        if ( params.cs==0 )
        {/*3*/  //Z=11907
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11908
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11909
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11910
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=11911
                fkv[n] = fkv[n-1]*n;  //Z=11912
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11913
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11914
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11915 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11916
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11917
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11918
                    /*  P(q)-coefficients  */  //Z=11919
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);       /*  P||(q)  */  //Z=11920
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  /*  P-(q)  */  //Z=11921
                    /*  F(q)-coefficients  */  //Z=11922
                    sump = 0.0;  //Z=11923
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11924
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11925
                    }/*6*/  //Z=11926
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=11927
                    sump = 0.0;  //Z=11928
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11929
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11930
                    }/*6*/  //Z=11931
                    params.CR->carr4f[n] = xrn[n]*sump;  //Z=11932
                }/*5*/  //Z=11933
                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=11934
                    /*  P(q)-coefficients  */  //Z=11935
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);    /*  P-(q)  */  //Z=11936
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);           /*  P||(q)  */  //Z=11937
                    /*  F(q)-coefficients  */  //Z=11938
                    sump = 0.0;  //Z=11939
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11940
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11941
                    }/*6*/  //Z=11942
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=11943
                    sump = 0.0;  //Z=11944
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11945
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11946
                    }/*6*/  //Z=11947
                    params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=11948
                }/*5*/  //Z=11949
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11950
                    if ( n<n1 ) n1 = n;  //Z=11951
                }/*5*/  //Z=11952
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11953
                    if ( n<n4 ) n4 = n;  //Z=11954
                }/*5*/  //Z=11955
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11956
                    if ( n<n1f ) n1f = n;  //Z=11957
                }/*5*/  //Z=11958
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11959
                    if ( n<n4f ) n4f = n;  //Z=11960
                }/*5*/  //Z=11961
            }/*4*/  //Z=11962
        }/*3*/  //Z=11963

        /*  core/shell  */  //Z=11965
        if ( params.cs==1 )
        {/*3*/  //Z=11966
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11967
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11968
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11969
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=11970
                fkv[n] = fkv[n-1]*n;  //Z=11971
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11972
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11973
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11974 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11975
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11976
                xrmn_n = -xrmn_n*xrm2z;  //Z=11977
                pn[n] = pn[n-1]*p*p;  //Z=11978
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11979
                    /*  P(q)-coefficients  */  //Z=11980
                    /*  longitudinal  */  //Z=11981
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11982
                    /*  cross-sectional  */  //Z=11983
                    /*  F121  */  //Z=11984
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11985
                    /*  F122  */  //Z=11986
                    sump = 0.0;  //Z=11987
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11988
                        sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11989
                    }/*6*/  //Z=11990
                    params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=11991
                    /*  F123  */  //Z=11992
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11993

                    /*  F(q)-coefficients  */  //Z=11995
                    /*  longitudinal  */  //Z=11996
                    sump = 0.0;  //Z=11997
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11998
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11999
                    }/*6*/  //Z=12000
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=12001
                    /*  cross-sectional  */  //Z=12002
                    /*  F121  */  //Z=12003
                    sump = 0.0;  //Z=12004
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12005
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12006
                    }/*6*/  //Z=12007
                    params.CR->carr4f[n] = xrn[n]*sump;  //Z=12008
                    /*  F122  */  //Z=12009
                    sump = 0.0;  //Z=12010
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12011
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12012
                    }/*6*/  //Z=12013
                    params.CR->carr5f[n] = xrmn_n*sump;  //Z=12014
                    /*  F123  */  //Z=12015
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=12016
                }/*5*/  //Z=12017

                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=12019
                    /*  P(q)-coefficients  */  //Z=12020
                    /*  longitudinal  */  //Z=12021
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12022
                    /*  cross-sectional  */  //Z=12023
                    /*  F121  */  //Z=12024
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12025
                    /*  F122  */  //Z=12026
                    /* sump:=0.0;  //Z=12027 */
                    /*    for m:=0 to n do begin  //Z=12028 */
                    /*    sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12029 */
                    /* end;  //Z=12030 */
                    /* carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=12031 */

                    /*  F122  */  //Z=12033
                    sump = 0.0;  //Z=12034
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12035
                        sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12036
                    }/*6*/  //Z=12037
                    params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=12038

                    /*  F123  */  //Z=12041
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=12042
                    /*  F(q)-coefficients  */  //Z=12043
                    /*  longitudinal  */  //Z=12044
                    sump = 0.0;  //Z=12045
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12046
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12047
                    }/*6*/  //Z=12048
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12049
                    /*  cross-sectional  */  //Z=12050
                    /*  F121  */  //Z=12051
                    sump = 0.0;  //Z=12052
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12053
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12054
                    }/*6*/  //Z=12055
                    params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=12056
                    /*  F122  */  //Z=12057
                    sump = 0.0;  //Z=12058
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12059
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12060
                    }/*6*/  //Z=12061
                    params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=12062
                    /*  F123  */  //Z=12063
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=12064
                }/*5*/  //Z=12065
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=12066
                    if ( n<n1 ) n1 = n;  //Z=12067
                }/*5*/  //Z=12068
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=12069
                    if ( n<n4 ) n4 = n;  //Z=12070
                }/*5*/  //Z=12071
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=12072
                    if ( n<n5 ) n5 = n;  //Z=12073
                }/*5*/  //Z=12074
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=12075
                    if ( n<n6 ) n6 = n;  //Z=12076
                }/*5*/  //Z=12077
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=12078
                    if ( n<n1f ) n1f = n;  //Z=12079
                }/*5*/  //Z=12080
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=12081
                    if ( n<n4f ) n4f = n;  //Z=12082
                }/*5*/  //Z=12083
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=12084
                    if ( n<n5f ) n5f = n;  //Z=12085
                }/*5*/  //Z=12086
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=12087
                    if ( n<n6f ) n6f = n;  //Z=12088
                }/*5*/  //Z=12089
            }/*4*/  //Z=12090
        }/*3*/  /*  of core/shell  */  //Z=12091

        /*  inhomogeneous core/shell  */  //Z=12093
        if ( params.cs==2 )
        {/*3*/  //Z=12094
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=12095
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12096
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12097
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=12098
                fkv[n] = fkv[n-1]*n;  //Z=12099
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12100
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12101
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12102 */
                xln[n] = -xln[n-1]*xl2z;  //Z=12103
                xrn[n] = -xrn[n-1]*xr2z;  //Z=12104
                xrmn_n = -xrmn_n*xrm2z;  //Z=12105
                pn[n] = pn[n-1]*p*p;  //Z=12106
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=12107
                    /*  P(q)-coefficients  */  //Z=12108
                    /*  longitudinal  */  //Z=12109
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12110

                    /*  cross-sectional P(q)  */  //Z=12112
                    params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12113
                    sump = 0.0;  //Z=12114
                    sump1 = 0.0;  //Z=12115
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12116
                        sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=12117
                        sump = sump+pn[n-m]*sumi;  //Z=12118
                        sump1 = sump1+sumi;  //Z=12119
                    }/*6*/  //Z=12120
                    params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=12121
                    params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n]*sump1;  //Z=12122
                    sump = 0.0;  //Z=12123
                    sump1 = 0.0;  //Z=12124
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12125
                        sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=12126
                        sump = sump+sumi;  //Z=12127
                        sump1 = sump1+pn[n-m]*sumi;  //Z=12128
                    }/*6*/  //Z=12129
                    params.CR->carr7p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=12130
                    params.CR->carr8p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=12131
                    params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=12132

                    /*  cross-sectional  */  //Z=12134
                    /*  (* F121 *)  //Z=12135
                   carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12136
                   (* F122 *)  //Z=12137
                   sump:=0.0;  //Z=12138
                   for m:=0 to n do begin  //Z=12139
                      sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12140
                   end;  //Z=12141
                   carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=12142
                   (* F123 *)  //Z=12143
                   carr6p[n]:=carr4p[n]/pn[n];    */  //Z=12144

                    /*  F(q)-coefficients  */  //Z=12146
                    /*  longitudinal  */  //Z=12147
                    sump = 0.0;  //Z=12148
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12149
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12150
                    }/*6*/  //Z=12151
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=12152
                    /*  cross-sectional  */  //Z=12153
                    params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=12154
                    params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12155
                    params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12156
                }/*5*/  //Z=12157

                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=12159
                    /*  P(q)-coefficients  */  //Z=12160
                    /*  longitudinal  */  //Z=12161
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12162

                    /*  cross-sectional P(q)  */  //Z=12164
                    params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=12165
                    sump = 0.0;  //Z=12166
                    sump1 = 0.0;  //Z=12167
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12168
                        sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=12169
                        sump = sump+pn[n-m]*sumi;  //Z=12170
                        sump1 = sump1+sumi;  //Z=12171
                    }/*6*/  //Z=12172
                    params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=12173
                    params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=12174
                    sump = 0.0;  //Z=12175
                    sump1 = 0.0;  //Z=12176
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12177
                        sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=12178
                        sump = sump+sumi;  //Z=12179
                        sump1 = sump1+pn[n-m]*sumi;  //Z=12180
                    }/*6*/  //Z=12181
                    params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=12182
                    params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=12183
                    params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=12184

                    /*  cross-sectional P(q)  */  //Z=12187
                    /*  F121  */  //Z=12188
                    /*  carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=12189
                        (* F122 *)  //Z=12190
                    //sump:=0.0;  //Z=12191
                    //   for m:=0 to n do begin  //Z=12192
                    //   sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12193
                    //end;  //Z=12194
                    //carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=12195

                    (* F122 *)  //Z=12197
                    sump:=0.0;  //Z=12198
                    for m:=0 to n do begin  //Z=12199
                       sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12200
                    end;  //Z=12201
                    carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=12202
                    (* F123 *)  //Z=12203
                    carr6p[n]:=carr4p[n]/pn[n];        */  //Z=12204

                    /*  F(q)-coefficients  */  //Z=12206
                    /*  longitudinal  */  //Z=12207
                    sump = 0.0;  //Z=12208
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12209
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12210
                    }/*6*/  //Z=12211
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12212
                    /*  cross-sectional  */  //Z=12213
                    params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=12214
                    params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12215
                    params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12216

                }/*5*/  //Z=12218
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=12219
                    if ( n<n1 ) n1 = n;  //Z=12220
                }/*5*/  //Z=12221
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=12222
                    if ( n<n4 ) n4 = n;  //Z=12223
                }/*5*/  //Z=12224
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=12225
                    if ( n<n5 ) n5 = n;  //Z=12226
                }/*5*/  //Z=12227
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=12228
                    if ( n<n6 ) n6 = n;  //Z=12229
                }/*5*/  //Z=12230
                if ( fabs(params.CR->carr7p[n])<min )
                {/*5*/  //Z=12231
                    if ( n<n7 ) n7 = n;  //Z=12232
                }/*5*/  //Z=12233
                if ( fabs(params.CR->carr8p[n])<min )
                {/*5*/  //Z=12234
                    if ( n<n8 ) n8 = n;  //Z=12235
                }/*5*/  //Z=12236
                if ( fabs(params.CR->carr9p[n])<min )
                {/*5*/  //Z=12237
                    if ( n<n9 ) n9 = n;  //Z=12238
                }/*5*/  //Z=12239
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=12240
                    if ( n<n1f ) n1f = n;  //Z=12241
                }/*5*/  //Z=12242
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=12243
                    if ( n<n4f ) n4f = n;  //Z=12244
                }/*5*/  //Z=12245
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=12246
                    if ( n<n5f ) n5f = n;  //Z=12247
                }/*5*/  //Z=12248
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=12249
                    if ( n<n6f ) n6f = n;  //Z=12250
                }/*5*/  //Z=12251
                if ( fabs(params.CR->carr7f[n])<min )
                {/*5*/  //Z=12252
                    if ( n<n7f ) n7f = n;  //Z=12253
                }/*5*/  //Z=12254
                if ( fabs(params.CR->carr8f[n])<min )
                {/*5*/  //Z=12255
                    if ( n<n8f ) n8f = n;  //Z=12256
                }/*5*/  //Z=12257
                if ( fabs(params.CR->carr9f[n])<min )
                {/*5*/  //Z=12258
                    if ( n<n9f ) n9f = n;  //Z=12259
                }/*5*/  //Z=12260
            }/*4*/  //Z=12261
        }/*3*/  /*  of inhomogeneous core/shell  */  //Z=12262


        /*  myelin  */  //Z=12265
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=12266
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=12267
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12268
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12269
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=12270
                fkv[n] = fkv[n-1]*n;  //Z=12271
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12272
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12273
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12274 */
                xln[n] = -xln[n-1]*xl2z;  //Z=12275
                xrn[n] = -xrn[n-1]*x12zm;  //Z=12276
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=12277
                    /*  P(q)-coefficients  */  //Z=12278
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12279

                    params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=12281

                    /*  F(q)-coefficients  */  //Z=12283
                    sump = 0.0;  //Z=12284
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12285
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12286
                    }/*6*/  //Z=12287
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=12288
                    sump = 0.0;  //Z=12289
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12290
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12291
                    }/*6*/  //Z=12292
                    params.CR->carr4f[n] = xrn[n]*sump;  //Z=12293
                }/*5*/  //Z=12294
                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=12295
                    /*  P(q)-coefficients  */  //Z=12296
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12297
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12298
                    /*  F(q)-coefficients  */  //Z=12299
                    sump = 0.0;  //Z=12300
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12301
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12302
                    }/*6*/  //Z=12303
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12304
                    sump = 0.0;  //Z=12305
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12306
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12307
                    }/*6*/  //Z=12308
                    params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=12309
                }/*5*/  //Z=12310
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=12311
                    if ( n<n1 ) n1 = n;  //Z=12312
                }/*5*/  //Z=12313
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=12314
                    if ( n<n4 ) n4 = n;  //Z=12315
                }/*5*/  //Z=12316
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=12317
                    if ( n<n1f ) n1f = n;  //Z=12318
                }/*5*/  //Z=12319
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=12320
                    if ( n<n4f ) n4f = n;  //Z=12321
                }/*5*/  //Z=12322
            }/*4*/  //Z=12323
        }/*3*/  /*  of myelin  */  //Z=12324

    }/*2*/  //Z=12328


    /* ** orientational distribution for cylinders and disks ** */  //Z=12330
    if ( (ordis==0) && (dim!=3) )
    {/*2*/  //Z=12331
        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,2,0,0,0,0,params.CR->carr1p,params.norm);  //Z=12332
        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,3,0,0,0,0,params.CR->carr1p,order);  //Z=12333
        order = order/params.norm;  //Z=12334

        /*  use phi=0 and rotate qx/qy-axis  */  //Z=12336
        if ( params.orcase==1 )
        {/*3*/  //Z=12337
            params.polPhi = 0;  //Z=12338
            params.p11 = -cos(params.polPhi*M_PI/180.0)*cos(params.polTheta*M_PI/180.0);       /*  = -cos(theta*pi/180);  //Z=12339 */
            params.p12 = sin(params.polPhi*M_PI/180.0);                                        /*  = 0;  //Z=12340 */
            params.p13 = cos(params.polPhi*M_PI/180.0)*sin(params.polTheta*M_PI/180.0);        /*  =  sin(theta*pi/180);  //Z=12341 */
            params.p21 = -cos(params.polPhi*M_PI/180.0);                                       /*  = -1;  //Z=12342 */
            params.p22 = -sin(params.polPhi*M_PI/180.0)*cos(params.polTheta*M_PI/180.0);       /*  = 0;  //Z=12343 */
            params.p23 = sin(params.polPhi*M_PI/180.0)*sin(params.polTheta*M_PI/180.0);        /*  = 0;  //Z=12344 */
            params.p31 = -sin(params.polTheta*M_PI/180.0);                                     /*  = 0;  //Z=12345 */
            params.p32 = 0;  //Z=12346
            params.p33 = -cos(params.polTheta*M_PI/180.0);                                      /*  = -cos(theta*pi/180);  //Z=12347 */

            for ( n=0; n<=nmax; n++ )
            {/*4*/  //Z=12349
                for ( m=0; m<=nmax; m++ )
                {/*5*/  //Z=12350
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,0,0,2*n,2*m,params.CR->carr1p,intl);  //Z=12351
                    intlar[n][m] = intl;  //Z=12352
                }/*5*/  //Z=12353
            }/*4*/  //Z=12354

            /* for n:=1 to nmax do begin  //Z=12356
                fkv[n]:=fkv[n-1]*n;  //Z=12357
                fk2v[n]:=fk2v[n-1]*(2*n-1)*(2*n);  //Z=12358
                if odd(n) then gam1[n]:=fkv[round((n-1)/2)]  //Z=12359
                   else gam1[n]:=fkv[n]*sqrt(pi)/(fkv[round(n/2)]*power(2,n));  //Z=12360
                gam2[n]:=(n/2)*fkv[n-1]*sqrt(pi)/(gam1[n]*power(2,n-1));  //Z=12361
                gam3[n]:=gam3[n-1]*(2*n+1)/2;  //Z=12362
             end;  */  //Z=12363

            /*  cylinder form factor coefficients  */  //Z=12365
            if ( dim==1 )
            {/*4*/  //Z=12366
                /*  homogeneous  */  //Z=12367
                params.p11 = sin(params.polTheta*M_PI/180.0);  //Z=12368
                params.p13 = cos(params.polTheta*M_PI/180.0);  //Z=12369
                if ( params.cs==0 )
                {/*5*/  //Z=12370
                    i = 2;  //Z=12371
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12372
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12373
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12374
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=12375
                        fkv[n] = fkv[n-1]*n;  //Z=12376
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12377
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12378
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12379 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12380
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12381
                        /*  longitudinal  */  //Z=12382
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12383
                            sump = 0.0;  //Z=12384
                            for ( ll=0; ll<=m; ll++ )  //Z=12385
                                sump = sump+pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12386
                            /* carr1pm[i]:=power(4,m)*power(p21*p21,n-m)*sump/(fkv[n-m]);  //Z=12387 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12388 */
                            /* i:=i+1;  //Z=12389 */
                            params.CR->carr11pm[n][m] = sump/(fkv[n-m]);  //Z=12390
                        }/*7*/  //Z=12391
                        /*  P(q)-coefficients  */  //Z=12392
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=12393 */

                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=12395
                        /* carr1p[n]:=(1/2)*power(4,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=12396 */


                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12399
                        /*  F(q)-coefficients  */  //Z=12400
                        sump = 0.0;  //Z=12401
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12402
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12403
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=12404
                        /*  cross-sectional  */  //Z=12405
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=12406
                        sump = 0.0;  //Z=12407
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12408
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=12409
                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=12410
                            if ( n<n1 ) n1 = n;  //Z=12411
                        }/*7*/  //Z=12412
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=12413
                            if ( n<n4 ) n4 = n;  //Z=12414
                        }/*7*/  //Z=12415
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12416
                            if ( n<n1f ) n1f = n;  //Z=12417
                        }/*7*/  //Z=12418
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12419
                            if ( n<n4f ) n4f = n;  //Z=12420
                        }/*7*/  //Z=12421
                    }/*6*/  //Z=12422
                }/*5*/  /*  of homogeneous cylinder  */  //Z=12423

                /*  core/shell  */  //Z=12425
                if ( params.cs==1 )
                {/*5*/  //Z=12426
                    i = 2;  //Z=12427
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12428
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12429
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12430
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=12431
                        fkv[n] = fkv[n-1]*n;  //Z=12432
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12433
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12434
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12435 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12436
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12437
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12438
                        pn[n] = pn[n-1]*p*p;  //Z=12439
                        /*  longitudinal  */  //Z=12440
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12441
                            sump = 0.0;  //Z=12442
                            for ( ll=0; ll<=m; ll++ )  //Z=12443
                                sump = sump+pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12444
                            /* carr1pm[i]:=power(4,m)*power(p21*p21,n-m)*sump/(fkv[n-m]);  //Z=12445 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12446 */
                            /* i:=i+1;  //Z=12447 */
                            params.CR->carr11pm[n][m] = sump/(fkv[n-m]);  //Z=12448
                        }/*7*/  //Z=12449
                        /*  P(q)-coefficients  */  //Z=12450
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=12451 */

                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=12453
                        /* carr1p[n]:=(1/2)*power(4,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=12454 */

                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12456
                        /*  F(q)-coefficients  */  //Z=12457
                        sump = 0.0;  //Z=12458
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12459
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12460
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=12461

                        /*  P(q)-coefficients  */  //Z=12463
                        /*  cross-sectional  */  //Z=12464
                        /*  F121  */  //Z=12465
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12466
                        /*  F122  */  //Z=12467
                        sump = 0.0;  //Z=12468
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12469
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12470
                        }/*7*/  //Z=12471
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=12472
                        /*  F123  */  //Z=12473
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=12474

                        /*  F(q)-coefficients  */  //Z=12476
                        /*  cross-sectional  */  //Z=12477
                        /*  F121  */  //Z=12478
                        sump = 0.0;  //Z=12479
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12480
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12481
                        }/*7*/  //Z=12482
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=12483
                        /*  F122  */  //Z=12484
                        sump = 0.0;  //Z=12485
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12486
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12487
                        }/*7*/  //Z=12488
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=12489
                        /*  F123  */  //Z=12490
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=12491

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=12493
                            if ( n<n1 ) n1 = n;  //Z=12494
                        }/*7*/  //Z=12495
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=12496
                            if ( n<n4 ) n4 = n;  //Z=12497
                        }/*7*/  //Z=12498
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=12499
                            if ( n<n5 ) n5 = n;  //Z=12500
                        }/*7*/  //Z=12501
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=12502
                            if ( n<n6 ) n6 = n;  //Z=12503
                        }/*7*/  //Z=12504
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12505
                            if ( n<n1f ) n1f = n;  //Z=12506
                        }/*7*/  //Z=12507
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12508
                            if ( n<n4f ) n4f = n;  //Z=12509
                        }/*7*/  //Z=12510
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12511
                            if ( n<n5f ) n5f = n;  //Z=12512
                        }/*7*/  //Z=12513
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12514
                            if ( n<n6f ) n6f = n;  //Z=12515
                        }/*7*/  //Z=12516
                    }/*6*/  /*  n-loop  */  //Z=12517
                }/*5*/  /*  homogeneous loop  */  //Z=12518

                /*  inhomogeneous core/shell  */  //Z=12522
                if ( params.cs==2 )
                {/*5*/  //Z=12523
                    i = 2;  //Z=12524
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12525
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12526
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12527
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=12528
                        fkv[n] = fkv[n-1]*n;  //Z=12529
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12530
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12531
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12532 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12533
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12534
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12535
                        pn[n] = pn[n-1]*p*p;  //Z=12536
                        /*  longitudinal  */  //Z=12537
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12538
                            sump = 0.0;  //Z=12539
                            for ( ll=0; ll<=m; ll++ )  //Z=12540
                                sump = sump+pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12541
                            /* carr1pm[i]:=power(4,m)*power(p21*p21,n-m)*sump/(fkv[n-m]);  //Z=12542 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12543 */
                            /* i:=i+1;  //Z=12544 */
                            params.CR->carr11pm[n][m] = sump/(fkv[n-m]);  //Z=12545
                        }/*7*/  //Z=12546
                        /*  P(q)-coefficients  */  //Z=12547
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=12548 */

                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=12550
                        /* carr1p[n]:=(1/2)*power(4,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=12551 */

                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12553
                        /*  F(q)-coefficients  */  //Z=12554
                        sump = 0.0;  //Z=12555
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12556
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12557
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=12558


                        /*  cross-sectional P(q)  */  //Z=12561
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12562
                        sump = 0.0;  //Z=12563
                        sump1 = 0.0;  //Z=12564
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12565
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=12566
                            sump = sump+pn[n-m]*sumi;  //Z=12567
                            sump1 = sump1+sumi;  //Z=12568
                        }/*7*/  //Z=12569
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=12570
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n]*sump1;  //Z=12571
                        sump = 0.0;  //Z=12572
                        sump1 = 0.0;  //Z=12573
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12574
                            sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=12575
                            sump = sump+sumi;  //Z=12576
                            sump1 = sump1+pn[n-m]*sumi;  //Z=12577
                        }/*7*/  //Z=12578
                        params.CR->carr7p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=12579
                        params.CR->carr8p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=12580
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=12581

                        /*  (* P(q)-coefficients *)  //Z=12585
                            (* cross-sectional *)  //Z=12586
                            (* F121 *)  //Z=12587
                            carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12588
                            (* F122 *)  //Z=12589
                            sump:=0.0;  //Z=12590
                               for m:=0 to n do begin  //Z=12591
                               sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12592
                            end;  //Z=12593
                            carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=12594
                            (* F123 *)  //Z=12595
                            carr6p[n]:=carr4p[n]/pn[n];   */  //Z=12596

                        /*  F(q)-coefficients  */  //Z=12598
                        params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=12599
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12600
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12601

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=12603
                            if ( n<n1 ) n1 = n;  //Z=12604
                        }/*7*/  //Z=12605
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=12606
                            if ( n<n4 ) n4 = n;  //Z=12607
                        }/*7*/  //Z=12608
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=12609
                            if ( n<n5 ) n5 = n;  //Z=12610
                        }/*7*/  //Z=12611
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=12612
                            if ( n<n6 ) n6 = n;  //Z=12613
                        }/*7*/  //Z=12614
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=12615
                            if ( n<n7 ) n7 = n;  //Z=12616
                        }/*7*/  //Z=12617
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=12618
                            if ( n<n8 ) n8 = n;  //Z=12619
                        }/*7*/  //Z=12620
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=12621
                            if ( n<n9 ) n9 = n;  //Z=12622
                        }/*7*/  //Z=12623
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12624
                            if ( n<n1f ) n1f = n;  //Z=12625
                        }/*7*/  //Z=12626
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12627
                            if ( n<n4f ) n4f = n;  //Z=12628
                        }/*7*/  //Z=12629
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12630
                            if ( n<n5f ) n5f = n;  //Z=12631
                        }/*7*/  //Z=12632
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12633
                            if ( n<n6f ) n6f = n;  //Z=12634
                        }/*7*/  //Z=12635
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=12636
                            if ( n<n7f ) n7f = n;  //Z=12637
                        }/*7*/  //Z=12638
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=12639
                            if ( n<n8f ) n8f = n;  //Z=12640
                        }/*7*/  //Z=12641
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=12642
                            if ( n<n9f ) n9f = n;  //Z=12643
                        }/*7*/  //Z=12644
                    }/*6*/   /*  of n-loop  */  //Z=12645
                }/*5*/  /*  of inhomogeneous core/shell  */  //Z=12646

                /*  myelin  */  //Z=12648
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=12649
                    i = 2;  //Z=12650
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12651
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12652
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12653
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=12654
                        fkv[n] = fkv[n-1]*n;  //Z=12655
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12656
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12657
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12658 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12659
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=12660
                        /*  longitudinal  */  //Z=12661
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12662
                            sump = 0.0;  //Z=12663
                            for ( ll=0; ll<=m; ll++ )  //Z=12664
                                sump = sump+pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12665
                            /* carr1pm[i]:=power(4,m)*power(p21*p21,n-m)*sump/(fkv[n-m]);  //Z=12666 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12667 */
                            /* i:=i+1;  //Z=12668 */
                            params.CR->carr11pm[n][m] = sump/(fkv[n-m]);  //Z=12669
                        }/*7*/  //Z=12670
                        /*  P(q)-coefficients  */  //Z=12671
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=12672 */

                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=12674
                        /* carr1p[n]:=(1/2)*power(4,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=12675 */

                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12677
                        /*  F(q)-coefficients  */  //Z=12678
                        sump = 0.0;  //Z=12679
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12680
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12681
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=12682
                        /*  cross-sectional  */  //Z=12683
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=12684

                        sump = 0.0;  //Z=12686
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12687
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=12688
                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=12689
                            if ( n<n1 ) n1 = n;  //Z=12690
                        }/*7*/  //Z=12691
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=12692
                            if ( n<n4 ) n4 = n;  //Z=12693
                        }/*7*/  //Z=12694
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12695
                            if ( n<n1f ) n1f = n;  //Z=12696
                        }/*7*/  //Z=12697
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12698
                            if ( n<n4f ) n4f = n;  //Z=12699
                        }/*7*/  //Z=12700
                    }/*6*/  //Z=12701
                }/*5*/        /*  of myelin  */  //Z=12702
            }/*4*/   /*  of cylinder  */  //Z=12703

            /* if (dim=1) then begin  //Z=12705
                for n:=1 to nmax do begin  //Z=12706
                   z12v[n]:=z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12707
                   b1sv[n]:=b1sv[n-1]*(b1s-1+n);  //Z=12708
                   xln[n]:=-xln[n-1]*xl2z;  //Z=12709
                   xrn[n]:=-xrn[n-1]*xr2z;  //Z=12710
                   carr1[n]:=power(4,2*n)*z12v[n]*xln[n]/((2*n+1)*(n+1));  //Z=12711
                   carr1p[n]:=carr1[n];  //Z=12712
                   carr3[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*b1sv[n]*fkv[n]);  //Z=12713
                   carr3p[n]:=carr3[n];  //Z=12714
                   if search1 then begin  //Z=12715
                      if carr1[n]<1e-50 then begin  //Z=12716
                         n1:=n;  //Z=12717
                         search1:=false;  //Z=12718
                      end;  //Z=12719
                   end;  //Z=12720
                   if search3 then begin  //Z=12721
                      if carr3[n]<1e-50 then begin  //Z=12722
                         n3:=n;  //Z=12723
                         search3:=false;  //Z=12724
                      end;  //Z=12725
                   end;  //Z=12726
                end;  //Z=12727
             end;      */  //Z=12728

            /*  disk form factor coefficients  */  //Z=12730
            if ( dim==2 )
            {/*4*/  //Z=12731
                /*  homogeneous  */  //Z=12732
                params.p11 = sin(params.polTheta*M_PI/180.0);  //Z=12733
                params.p13 = cos(params.polTheta*M_PI/180.0);  //Z=12734
                if ( params.cs==0 )
                {/*5*/  //Z=12735
                    i = 2;  //Z=12736
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12737
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12738
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12739
                        fkv[n] = fkv[n-1]*n;  //Z=12740
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12741
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12742
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12743 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12744
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12745
                        /*  longitudinal  */  //Z=12746
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12747
                            sump = 0.0;  //Z=12748
                            for ( ll=0; ll<=m; ll++ )  //Z=12749
                                sump = sump+pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12750
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12751 */
                            params.CR->carr22pm[n][m] = sump/(fkv[n-m]);  //Z=12752
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=12753 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=12754
                            /* carr2fm[i]:=carr2pm[i];  //Z=12755 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12756 */
                            i = i+1;  //Z=12757
                        }/*7*/  //Z=12758
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=12759

                        sump1 = 0.0;  //Z=12761
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12762
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12763
                        }/*7*/  //Z=12764
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12765

                        /*  cross-sectional  */  //Z=12767
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12768
                        sump = 0.0;  //Z=12769
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12770
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12771
                        }/*7*/  //Z=12772
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=12773

                        /*  series for <...> integration  */  //Z=12775
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12776
                        sump = 0.0;  //Z=12777
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12778
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12779
                        }/*7*/  //Z=12780
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12781

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=12783
                            if ( n<n1 ) n1 = n;  //Z=12784
                        }/*7*/  //Z=12785
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=12786
                            if ( n<n4 ) n4 = n;  //Z=12787
                        }/*7*/  //Z=12788
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12789
                            if ( n<n1f ) n1f = n;  //Z=12790
                        }/*7*/  //Z=12791
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12792
                            if ( n<n4f ) n4f = n;  //Z=12793
                        }/*7*/  //Z=12794
                    }/*6*/  //Z=12795
                }/*5*/  //Z=12796

                /*  core/shell  */  //Z=12798
                if ( params.cs==1 )
                {/*5*/  //Z=12799
                    i = 2;  //Z=12800
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12801
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12802
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12803
                        fkv[n] = fkv[n-1]*n;  //Z=12804
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12805
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12806
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12807 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12808
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12809
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12810
                        pn[n] = pn[n-1]*p*p;  //Z=12811
                        /*  longitudinal  */  //Z=12812
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12813
                            sump = 0.0;  //Z=12814
                            for ( ll=0; ll<=m; ll++ )  //Z=12815
                                sump = sump+pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12816
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12817 */
                            params.CR->carr22pm[n][m] = sump/(fkv[n-m]);  //Z=12818
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=12819 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=12820
                            /* carr2fm[i]:=carr2pm[i];  //Z=12821 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12822 */
                            i = i+1;  //Z=12823
                        }/*7*/  //Z=12824
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=12825

                        sump1 = 0.0;  //Z=12827
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12828
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12829
                        }/*7*/  //Z=12830
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12831

                        /*  P(q)-coefficients  */  //Z=12833
                        /*  cross-sectional  */  //Z=12834
                        /*  F121  */  //Z=12835
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12836
                        /*  F122  */  //Z=12837
                        sump = 0.0;  //Z=12838
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12839
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12840
                        }/*7*/  //Z=12841
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=12842

                        /*  F122  */  //Z=12845
                        sump = 0.0;  //Z=12846
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12847
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12848
                        }/*7*/  //Z=12849
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=12850

                        /*  F123  */  //Z=12853
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=12854
                        /*  F121  */  //Z=12855
                        sump = 0.0;  //Z=12856
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12857
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12858
                        }/*7*/  //Z=12859
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=12860
                        /*  F122  */  //Z=12861
                        sump = 0.0;  //Z=12862
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12863
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12864
                        }/*7*/  //Z=12865
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=12866
                        /*  F123  */  //Z=12867
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=12868

                        /*  series for <...> integration  */  //Z=12870
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12871
                        sump = 0.0;  //Z=12872
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12873
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12874
                        }/*7*/  //Z=12875
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12876

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=12878
                            if ( n<n1 ) n1 = n;  //Z=12879
                        }/*7*/  //Z=12880
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=12881
                            if ( n<n4 ) n4 = n;  //Z=12882
                        }/*7*/  //Z=12883
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=12884
                            if ( n<n5 ) n5 = n;  //Z=12885
                        }/*7*/  //Z=12886
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=12887
                            if ( n<n6 ) n6 = n;  //Z=12888
                        }/*7*/  //Z=12889
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12890
                            if ( n<n1f ) n1f = n;  //Z=12891
                        }/*7*/  //Z=12892
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12893
                            if ( n<n4f ) n4f = n;  //Z=12894
                        }/*7*/  //Z=12895
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12896
                            if ( n<n5f ) n5f = n;  //Z=12897
                        }/*7*/  //Z=12898
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12899
                            if ( n<n6f ) n6f = n;  //Z=12900
                        }/*7*/  //Z=12901
                    }/*6*/ /*  of n-loop  */  //Z=12902
                }/*5*/  /*  of core/shell  */  //Z=12903

                /*  inhomogeneous core/shell  */  //Z=12905
                if ( params.cs==2 )
                {/*5*/  //Z=12906
                    i = 2;  //Z=12907
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12908
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12909
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12910
                        fkv[n] = fkv[n-1]*n;  //Z=12911
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12912
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12913
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12914 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12915
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12916
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12917
                        pn[n] = pn[n-1]*p*p;  //Z=12918
                        /*  longitudinal  */  //Z=12919
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12920
                            sump = 0.0;  //Z=12921
                            for ( ll=0; ll<=m; ll++ )  //Z=12922
                                sump = sump+pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12923
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12924 */
                            params.CR->carr22pm[n][m] = sump/(fkv[n-m]);  //Z=12925
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=12926 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=12927
                            /* carr2fm[i]:=carr2pm[i];  //Z=12928 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12929 */
                            i = i+1;  //Z=12930
                        }/*7*/  //Z=12931
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=12932

                        sump1 = 0.0;  //Z=12934
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12935
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12936
                        }/*7*/  //Z=12937
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12938

                        /*  cross-sectional P(q)  */  //Z=12940
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=12941
                        sump = 0.0;  //Z=12942
                        sump1 = 0.0;  //Z=12943
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12944
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=12945
                            sump = sump+pn[n-m]*sumi;  //Z=12946
                            sump1 = sump1+sumi;  //Z=12947
                        }/*7*/  //Z=12948
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=12949
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=12950
                        sump = 0.0;  //Z=12951
                        sump1 = 0.0;  //Z=12952
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12953
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=12954
                            sump = sump+sumi;  //Z=12955
                            sump1 = sump1+pn[n-m]*sumi;  //Z=12956
                        }/*7*/  //Z=12957
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=12958
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=12959
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=12960

                        /*  (* P(q)-coefficients *)  //Z=12962
                            (* cross-sectional *)  //Z=12963
                            (* F121 *)  //Z=12964
                            carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=12965
                            (* F122 *)  //Z=12966
                            sump:=0.0;  //Z=12967
                               for m:=0 to n do begin  //Z=12968
                               sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12969
                            end;  //Z=12970
                            carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=12971
                            (* F122 *)  //Z=12972
                            sump:=0.0;  //Z=12973
                            for m:=0 to n do begin  //Z=12974
                               sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12975
                            end;  //Z=12976
                            carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=12977
                            (* F123 *)  //Z=12978
                            carr6p[n]:=carr4p[n]/pn[n];      */  //Z=12979

                        /*  F(q)  */  //Z=12981
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=12982
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12983
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12984

                        /*  series for <...> integration  */  //Z=12986
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12987
                        sump = 0.0;  //Z=12988
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12989
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12990
                        }/*7*/  //Z=12991
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12992

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=12994
                            if ( n<n1 ) n1 = n;  //Z=12995
                        }/*7*/  //Z=12996
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=12997
                            if ( n<n4 ) n4 = n;  //Z=12998
                        }/*7*/  //Z=12999
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=13000
                            if ( n<n5 ) n5 = n;  //Z=13001
                        }/*7*/  //Z=13002
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=13003
                            if ( n<n6 ) n6 = n;  //Z=13004
                        }/*7*/  //Z=13005
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=13006
                            if ( n<n7 ) n7 = n;  //Z=13007
                        }/*7*/  //Z=13008
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=13009
                            if ( n<n8 ) n8 = n;  //Z=13010
                        }/*7*/  //Z=13011
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=13012
                            if ( n<n9 ) n9 = n;  //Z=13013
                        }/*7*/  //Z=13014
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13015
                            if ( n<n1f ) n1f = n;  //Z=13016
                        }/*7*/  //Z=13017
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13018
                            if ( n<n4f ) n4f = n;  //Z=13019
                        }/*7*/  //Z=13020
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13021
                            if ( n<n5f ) n5f = n;  //Z=13022
                        }/*7*/  //Z=13023
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13024
                            if ( n<n6f ) n6f = n;  //Z=13025
                        }/*7*/  //Z=13026
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=13027
                            if ( n<n7f ) n7f = n;  //Z=13028
                        }/*7*/  //Z=13029
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=13030
                            if ( n<n8f ) n8f = n;  //Z=13031
                        }/*7*/  //Z=13032
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=13033
                            if ( n<n9f ) n9f = n;  //Z=13034
                        }/*7*/  //Z=13035
                    }/*6*/ /*  of n-loop  */  //Z=13036
                }/*5*/  /*  of inhomogeneous core/shell  */  //Z=13037

            }/*4*/   /*  of disk  */  //Z=13041

        }/*3*/   /*  of cho1=1  */  //Z=13158

        /* ** general orientation case ** */  //Z=13160
        /*  for phi<>0, too slow, only for cylinders  */  //Z=13161
        if ( params.orcase==5 )
        {/*3*/  //Z=13162

            params.p11 = -cos(params.polPhi*M_PI/180.0)*cos(params.polTheta*M_PI/180.0);  //Z=13164
            params.p12 = sin(params.polPhi*M_PI/180.0);  //Z=13165
            params.p13 = cos(params.polPhi*M_PI/180.0)*sin(params.polTheta*M_PI/180.0);  //Z=13166
            params.p21 = -cos(params.polPhi*M_PI/180.0);  //Z=13167
            params.p22 = -sin(params.polPhi*M_PI/180.0)*cos(params.polTheta*M_PI/180.0);  //Z=13168
            params.p23 = sin(params.polPhi*M_PI/180.0)*sin(params.polTheta*M_PI/180.0);  //Z=13169
            params.p31 = -sin(params.polTheta*M_PI/180.0);  //Z=13170
            params.p32 = 0;  //Z=13171
            params.p33 = -cos(params.polTheta*M_PI/180.0);  //Z=13172

            for ( n=1; n<=2*nmax; n++ ) fkv[n] = fkv[n-1]*n;  //Z=13174

            i = 2;  //Z=13176
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=13177
                for ( m=0; m<=2*n; m++ )
                {/*5*/  //Z=13178
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,5,0,m,2*n-m,params.CR->carr1p,intl);  //Z=13179
                    /* carr1pm[i]:=intl/(fkv[m]*fkv[2*n-m]*norm);  //Z=13180 */
                    i = i+1;  //Z=13181
                }/*5*/  //Z=13182
            }/*4*/  //Z=13183

            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=13185
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13186
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=13187
                xln[n] = -xln[n-1]*xl2z;  //Z=13188
                xrn[n] = -xrn[n-1]*xr2z;  //Z=13189
                params.CR->carr1p[n] = pow(4.0,2*n)*z12v[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=13190
                params.CR->carr3p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*b1sv_n*fkv[n]);  //Z=13191

                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=13193
                    if ( n<n1 ) n1 = n;  //Z=13194
                }/*5*/  //Z=13195
                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=13196
                    if ( n<n3 ) n3 = n;  //Z=13197
                }/*5*/  //Z=13198
            }/*4*/  //Z=13199
        }/*3*/   /*  of cho1=5  */  //Z=13200

        /* ** x-axis, orientational distribution centered around x-axis ** */  //Z=13205
        if ( (params.orcase==2) && (dim!=3) )
        {/*3*/  //Z=13206
            /* ** cylinders ** */  //Z=13207
            if ( dim==1 )
            {/*4*/  //Z=13208
                /*  homogeneous  */  //Z=13209
                if ( params.cs==0 )
                {/*5*/  //Z=13210
                    i = 2;  //Z=13211
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13212
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13213
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13214
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=13215
                        fkv[n] = fkv[n-1]*n;  //Z=13216
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13217
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13218
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13219 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13220
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13221
                        /*  longitudinal  */  //Z=13222
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13223
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13224
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13225
                            /* carr1fm[i]:=carr1pm[i];  //Z=13226 */
                            /* i:=i+1;  //Z=13227 */
                        }/*7*/  //Z=13228
                        /*  P(q)-coefficient  */  //Z=13229
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=13230 */
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13231 */
                        params.CR->carr1p[n] = (1/2.0)*pow(4.0,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2.0));  //Z=13232

                        /*  F(q)-coefficient  */  //Z=13235
                        sump = 0.0;  //Z=13236
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13237
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=13238

                        /*  cross-section  */  //Z=13240
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=13241
                        sump = 0.0;  //Z=13242
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13243
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13244

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13246
                            if ( n<n1 ) n1 = n;  //Z=13247
                        }/*7*/  //Z=13248
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13249
                            if ( n<n4 ) n4 = n;  //Z=13250
                        }/*7*/  //Z=13251
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13252
                            if ( n<n1f ) n1f = n;  //Z=13253
                        }/*7*/  //Z=13254
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13255
                            if ( n<n4f ) n4f = n;  //Z=13256
                        }/*7*/  //Z=13257
                    }/*6*/  /*  of n-loop  */  //Z=13258
                }/*5*/  /*  of cs=0  */  //Z=13259

                /*  core/shell  */  //Z=13261
                if ( params.cs==1 )
                {/*5*/  //Z=13262
                    i = 2;  //Z=13263
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13264
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13265
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13266
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=13267
                        fkv[n] = fkv[n-1]*n;  //Z=13268
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13269
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13270
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13271 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13272
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13273
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13274
                        pn[n] = pn[n-1]*p*p;  //Z=13275
                        /*  longitudinal  */  //Z=13276
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13277
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13278
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13279
                            /* carr1fm[i]:=carr1pm[i];  //Z=13280 */
                            /* i:=i+1;  //Z=13281 */
                        }/*7*/  //Z=13282
                        /*  P(q)-coefficient  */  //Z=13283
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=13284 */
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13285 */
                        params.CR->carr1p[n] = (1/2.0)*pow(4.0,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2.0));  //Z=13286

                        /*  F(q)-coefficient  */  //Z=13288
                        sump = 0.0;  //Z=13289
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13290
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=13291

                        /*  P(q)-coefficients  */  //Z=13293
                        /*  cross-sectional  */  //Z=13294
                        /*  F121  */  //Z=13295
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13296
                        /*  F122  */  //Z=13297
                        sump = 0.0;  //Z=13298
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13299
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13300
                        }/*7*/  //Z=13301
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=13302
                        /*  F123  */  //Z=13303
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=13304

                        /*  F(q)-coefficients  */  //Z=13306
                        /*  cross-sectional  */  //Z=13307
                        /*  F121  */  //Z=13308
                        sump = 0.0;  //Z=13309
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13310
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13311
                        }/*7*/  //Z=13312
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13313
                        /*  F122  */  //Z=13314
                        sump = 0.0;  //Z=13315
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13316
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13317
                        }/*7*/  //Z=13318
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=13319
                        /*  F123  */  //Z=13320
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=13321

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13323
                            if ( n<n1 ) n1 = n;  //Z=13324
                        }/*7*/  //Z=13325
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13326
                            if ( n<n4 ) n4 = n;  //Z=13327
                        }/*7*/  //Z=13328
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=13329
                            if ( n<n5 ) n5 = n;  //Z=13330
                        }/*7*/  //Z=13331
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=13332
                            if ( n<n6 ) n6 = n;  //Z=13333
                        }/*7*/  //Z=13334
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13335
                            if ( n<n1f ) n1f = n;  //Z=13336
                        }/*7*/  //Z=13337
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13338
                            if ( n<n4f ) n4f = n;  //Z=13339
                        }/*7*/  //Z=13340
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13341
                            if ( n<n5f ) n5f = n;  //Z=13342
                        }/*7*/  //Z=13343
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13344
                            if ( n<n6f ) n6f = n;  //Z=13345
                        }/*7*/  //Z=13346
                    }/*6*/  /*  of n-loop  */  //Z=13347
                }/*5*/  /*  of cs=1  */  //Z=13348

                /*  inhomogeneous core/shell  */  //Z=13350
                if ( params.cs==2 )
                {/*5*/  //Z=13351
                    i = 2;  //Z=13352
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13353
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13354
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13355
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=13356
                        fkv[n] = fkv[n-1]*n;  //Z=13357
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13358
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13359
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13360 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13361
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13362
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13363
                        pn[n] = pn[n-1]*p*p;  //Z=13364
                        /*  longitudinal  */  //Z=13365
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13366
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13367
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13368
                            /* carr1fm[i]:=carr1pm[i];  //Z=13369 */
                            /* i:=i+1;  //Z=13370 */
                        }/*7*/  //Z=13371
                        /*  P(q)-coefficient  */  //Z=13372
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=13373 */
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13374 */
                        params.CR->carr1p[n] = (1/2.0)*pow(4.0,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2.0));  //Z=13375

                        /*  F(q)-coefficient  */  //Z=13377
                        sump = 0.0;  //Z=13378
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13379
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=13380

                        /*  cross-sectional P(q)  */  //Z=13383
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13384
                        sump = 0.0;  //Z=13385
                        sump1 = 0.0;  //Z=13386
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13387
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=13388
                            sump = sump+pn[n-m]*sumi;  //Z=13389
                            sump1 = sump1+sumi;  //Z=13390
                        }/*7*/  //Z=13391
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=13392
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n]*sump1;  //Z=13393
                        sump = 0.0;  //Z=13394
                        sump1 = 0.0;  //Z=13395
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13396
                            sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=13397
                            sump = sump+sumi;  //Z=13398
                            sump1 = sump1+pn[n-m]*sumi;  //Z=13399
                        }/*7*/  //Z=13400
                        params.CR->carr7p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=13401
                        params.CR->carr8p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=13402
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=13403

                        /*  F(q)-coefficients  */  //Z=13418
                        /*  cross-sectional  */  //Z=13419
                        params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=13420
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=13421
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=13422

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13424
                            if ( n<n1 ) n1 = n;  //Z=13425
                        }/*7*/  //Z=13426
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13427
                            if ( n<n4 ) n4 = n;  //Z=13428
                        }/*7*/  //Z=13429
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=13430
                            if ( n<n5 ) n5 = n;  //Z=13431
                        }/*7*/  //Z=13432
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=13433
                            if ( n<n6 ) n6 = n;  //Z=13434
                        }/*7*/  //Z=13435
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=13436
                            if ( n<n7 ) n7 = n;  //Z=13437
                        }/*7*/  //Z=13438
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=13439
                            if ( n<n8 ) n8 = n;  //Z=13440
                        }/*7*/  //Z=13441
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=13442
                            if ( n<n9 ) n9 = n;  //Z=13443
                        }/*7*/  //Z=13444
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13445
                            if ( n<n1f ) n1f = n;  //Z=13446
                        }/*7*/  //Z=13447
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13448
                            if ( n<n4f ) n4f = n;  //Z=13449
                        }/*7*/  //Z=13450
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13451
                            if ( n<n5f ) n5f = n;  //Z=13452
                        }/*7*/  //Z=13453
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13454
                            if ( n<n6f ) n6f = n;  //Z=13455
                        }/*7*/  //Z=13456
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=13457
                            if ( n<n7f ) n7f = n;  //Z=13458
                        }/*7*/  //Z=13459
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=13460
                            if ( n<n8f ) n8f = n;  //Z=13461
                        }/*7*/  //Z=13462
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=13463
                            if ( n<n9f ) n9f = n;  //Z=13464
                        }/*7*/  //Z=13465
                    }/*6*/  /*  of n-loop  */  //Z=13466
                }/*5*/  /*  of cs=1  */  //Z=13467

                /*  myelin  */  //Z=13469
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=13470
                    i = 2;  //Z=13471
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13472
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13473
                        z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13474
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=13475
                        fkv[n] = fkv[n-1]*n;  //Z=13476
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13477
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13478
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13479 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13480
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=13481
                        /*  longitudinal  */  //Z=13482
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13483
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13484
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13485
                            /* carr1fm[i]:=carr1pm[i];  //Z=13486 */
                            /* i:=i+1;  //Z=13487 */
                        }/*7*/  //Z=13488
                        /*  P(q)-coefficient  */  //Z=13489
                        /* carr1p[n]:=power(4,n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=13490 */
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13491 */
                        params.CR->carr1p[n] = (1/2.0)*pow(4.0,n)*z12vl[n]*xln[n]/((n+1)*(n+1/2.0));  //Z=13492

                        /*  F(q)-coefficient  */  //Z=13494
                        sump = 0.0;  //Z=13495
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13496
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=13497

                        /*  cross-section  */  //Z=13499
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=13500

                        sump = 0.0;  //Z=13502
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13503
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13504

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13506
                            if ( n<n1 ) n1 = n;  //Z=13507
                        }/*7*/  //Z=13508
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13509
                            if ( n<n4 ) n4 = n;  //Z=13510
                        }/*7*/  //Z=13511
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13512
                            if ( n<n1f ) n1f = n;  //Z=13513
                        }/*7*/  //Z=13514
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13515
                            if ( n<n4f ) n4f = n;  //Z=13516
                        }/*7*/  //Z=13517
                    }/*6*/  /*  of n-loop  */  //Z=13518
                }/*5*/  /*  of cs=3  */  //Z=13519
            }/*4*/  /*  of cylinders  */  //Z=13520

            /* ** disks ** */  //Z=13522
            if ( dim==2 )
            {/*4*/  //Z=13523
                /*  homogeneous  */  //Z=13524
                if ( params.cs==0 )
                {/*5*/  //Z=13525
                    i = 2;  //Z=13526
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13527
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13528
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13529
                        fkv[n] = fkv[n-1]*n;  //Z=13530
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13531
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13532
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13533 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13534
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13535
                        /*  longitudinal  */  //Z=13536
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13537
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13538
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=13539 */
                            params.CR->carr22pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13540
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=13541 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=13542
                            /* carr2fm[i]:=carr2pm[i];  //Z=13543 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13544 */
                            /* i:=i+1;  //Z=13545 */
                        }/*7*/  //Z=13546
                        /*  P(q)-coefficient  */  //Z=13547
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13548

                        /*  F(q)-coefficient  */  //Z=13550
                        sump = 0.0;  //Z=13551
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13552
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13553
                        }/*7*/  //Z=13554
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13555

                        /*  cross-section  */  //Z=13557
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13558
                        sump = 0.0;  //Z=13559
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13560
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13561
                        }/*7*/  //Z=13562
                        params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=13563

                        /*  series for <...> integration  */  //Z=13565
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13566
                        sump = 0.0;  //Z=13567
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13568
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13569
                        }/*7*/  //Z=13570
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13571

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13573
                            if ( n<n1 ) n1 = n;  //Z=13574
                        }/*7*/  //Z=13575
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13576
                            if ( n<n4 ) n4 = n;  //Z=13577
                        }/*7*/  //Z=13578
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13579
                            if ( n<n1f ) n1f = n;  //Z=13580
                        }/*7*/  //Z=13581
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13582
                            if ( n<n4f ) n4f = n;  //Z=13583
                        }/*7*/  //Z=13584
                    }/*6*/ /*  n-loop  */  //Z=13585
                }/*5*/  /*  of cs=0  */  //Z=13586

                /*  core/shell  */  //Z=13588
                if ( params.cs==1 )
                {/*5*/  //Z=13589
                    i = 2;  //Z=13590
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13591
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13592
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13593
                        fkv[n] = fkv[n-1]*n;  //Z=13594
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13595
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13596
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13597 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13598
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13599
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13600
                        pn[n] = pn[n-1]*p*p;  //Z=13601
                        /*  longitudinal  */  //Z=13602
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13603
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13604
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=13605 */
                            params.CR->carr22pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13606
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=13607 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=13608
                            /* carr2fm[i]:=carr2pm[i];  //Z=13609 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13610 */
                            /* i:=i+1;  //Z=13611 */
                        }/*7*/  //Z=13612
                        /*  P(q)-coefficient  */  //Z=13613
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13614

                        /*  F(q)-coefficient  */  //Z=13616
                        sump = 0.0;  //Z=13617
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13618
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13619
                        }/*7*/  //Z=13620
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13621

                        /*  F121  */  //Z=13623
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13624
                        /*  F122  */  //Z=13625
                        sump = 0.0;  //Z=13626
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13627
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13628
                        }/*7*/  //Z=13629
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=13630

                        /*  F122  */  //Z=13632
                        sump = 0.0;  //Z=13633
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13634
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13635
                        }/*7*/  //Z=13636
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=13637

                        /*  F123  */  //Z=13640
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=13641
                        /*  F121  */  //Z=13642
                        sump = 0.0;  //Z=13643
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13644
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13645
                        }/*7*/  //Z=13646
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=13647
                        /*  F122  */  //Z=13648
                        sump = 0.0;  //Z=13649
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13650
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13651
                        }/*7*/  //Z=13652
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=13653
                        /*  F123  */  //Z=13654
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=13655

                        /*  series for <...> integration  */  //Z=13657
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13658
                        sump = 0.0;  //Z=13659
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13660
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13661
                        }/*7*/  //Z=13662
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13663

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13665
                            if ( n<n1 ) n1 = n;  //Z=13666
                        }/*7*/  //Z=13667
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13668
                            if ( n<n4 ) n4 = n;  //Z=13669
                        }/*7*/  //Z=13670
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=13671
                            if ( n<n5 ) n5 = n;  //Z=13672
                        }/*7*/  //Z=13673
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=13674
                            if ( n<n6 ) n6 = n;  //Z=13675
                        }/*7*/  //Z=13676
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13677
                            if ( n<n1f ) n1f = n;  //Z=13678
                        }/*7*/  //Z=13679
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13680
                            if ( n<n4f ) n4f = n;  //Z=13681
                        }/*7*/  //Z=13682
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13683
                            if ( n<n5f ) n5f = n;  //Z=13684
                        }/*7*/  //Z=13685
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13686
                            if ( n<n6f ) n6f = n;  //Z=13687
                        }/*7*/  //Z=13688
                    }/*6*/ /*  n-loop  */  //Z=13689
                }/*5*/  /*  of cs=1  */  //Z=13690

                /*  inhomogeneous core/shell  */  //Z=13692
                if ( params.cs==2 )
                {/*5*/  //Z=13693
                    i = 2;  //Z=13694
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13695
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13696
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13697
                        fkv[n] = fkv[n-1]*n;  //Z=13698
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13699
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13700
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13701 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13702
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13703
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13704
                        pn[n] = pn[n-1]*p*p;  //Z=13705
                        /*  longitudinal  */  //Z=13706
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13707
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13708
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=13709 */
                            params.CR->carr22pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13710
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=13711 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=13712
                            /* carr2fm[i]:=carr2pm[i];  //Z=13713 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13714 */
                            /* i:=i+1;  //Z=13715 */
                        }/*7*/  //Z=13716
                        /*  P(q)-coefficient  */  //Z=13717
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13718

                        /*  F(q)-coefficient  */  //Z=13720
                        sump = 0.0;  //Z=13721
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13722
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13723
                        }/*7*/  //Z=13724
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13725

                        /*  cross-sectional P(q)  */  //Z=13727
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13728
                        sump = 0.0;  //Z=13729
                        sump1 = 0.0;  //Z=13730
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13731
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=13732
                            sump = sump+pn[n-m]*sumi;  //Z=13733
                            sump1 = sump1+sumi;  //Z=13734
                        }/*7*/  //Z=13735
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=13736
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=13737
                        sump = 0.0;  //Z=13738
                        sump1 = 0.0;  //Z=13739
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13740
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=13741
                            sump = sump+sumi;  //Z=13742
                            sump1 = sump1+pn[n-m]*sumi;  //Z=13743
                        }/*7*/  //Z=13744
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=13745
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=13746
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=13747

                        /*  F(q)   */  //Z=13767
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=13768
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=13769
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=13770

                        /*  series for <...> integration  */  //Z=13772
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13773
                        sump = 0.0;  //Z=13774
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13775
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13776
                        }/*7*/  //Z=13777
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13778

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13780
                            if ( n<n1 ) n1 = n;  //Z=13781
                        }/*7*/  //Z=13782
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13783
                            if ( n<n4 ) n4 = n;  //Z=13784
                        }/*7*/  //Z=13785
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=13786
                            if ( n<n5 ) n5 = n;  //Z=13787
                        }/*7*/  //Z=13788
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=13789
                            if ( n<n6 ) n6 = n;  //Z=13790
                        }/*7*/  //Z=13791
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=13792
                            if ( n<n7 ) n7 = n;  //Z=13793
                        }/*7*/  //Z=13794
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=13795
                            if ( n<n8 ) n8 = n;  //Z=13796
                        }/*7*/  //Z=13797
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=13798
                            if ( n<n9 ) n9 = n;  //Z=13799
                        }/*7*/  //Z=13800
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13801
                            if ( n<n1f ) n1f = n;  //Z=13802
                        }/*7*/  //Z=13803
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13804
                            if ( n<n4f ) n4f = n;  //Z=13805
                        }/*7*/  //Z=13806
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13807
                            if ( n<n5f ) n5f = n;  //Z=13808
                        }/*7*/  //Z=13809
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13810
                            if ( n<n6f ) n6f = n;  //Z=13811
                        }/*7*/  //Z=13812
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=13813
                            if ( n<n7f ) n7f = n;  //Z=13814
                        }/*7*/  //Z=13815
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=13816
                            if ( n<n8f ) n8f = n;  //Z=13817
                        }/*7*/  //Z=13818
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=13819
                            if ( n<n9f ) n9f = n;  //Z=13820
                        }/*7*/  //Z=13821
                    }/*6*/ /*  n-loop  */  //Z=13822
                }/*5*/  /*  of inhomogeneous core/shell  */  //Z=13823

            }/*4*/  /*  of disk  */  //Z=13825
        }/*3*/  //Z=13826

        /* ** y-axis, orientational distribution centered around y-axis ** */  //Z=13830
        if ( (params.orcase==3) && (dim!=3) )
        {/*3*/  //Z=13831

            /* ** cylinder ** */  //Z=13833
            if ( dim==1 )
            {/*4*/  //Z=13834
                /*  homogeneous  */  //Z=13835
                if ( params.cs==0 )
                {/*5*/  //Z=13836
                    i = 2;  //Z=13837
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13838
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13839
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13840
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=13841
                        fkv[n] = fkv[n-1]*n;  //Z=13842
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13843
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13844
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13845 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13846
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13847
                        /*  longitudinal  */  //Z=13848
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13849
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13850
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13851

                            /* carr1fm[i]:=carr1pm[i];  //Z=13853 */
                            /* i:=i+1;  //Z=13854 */
                        }/*7*/  //Z=13855
                        /*  P(q)-coefficient  */  //Z=13856
                        /* carr1p[n]:=power(4,2*n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=13857 */
                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13858
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13859 */
                        /* carr1p[n]:=(1/2)*power(4,2*n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=13860 */

                        /*  F(q)-coefficient  */  //Z=13862
                        sump = 0.0;  //Z=13863
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13864
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=13865

                        /*  cross-section  */  //Z=13867
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=13868

                        /* carr4p[n]:=power(4,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(pi)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13870 */

                        sump = 0.0;  //Z=13872
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13873
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13874

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13876
                            if ( n<n1 ) n1 = n;  //Z=13877
                        }/*7*/  //Z=13878
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13879
                            if ( n<n4 ) n4 = n;  //Z=13880
                        }/*7*/  //Z=13881
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13882
                            if ( n<n1f ) n1f = n;  //Z=13883
                        }/*7*/  //Z=13884
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13885
                            if ( n<n4f ) n4f = n;  //Z=13886
                        }/*7*/  //Z=13887
                    }/*6*/  /*  of n-loop  */  //Z=13888
                }/*5*/  /*  cs=0  */  //Z=13889

                /*  core/shell  */  //Z=13891
                if ( params.cs==1 )
                {/*5*/  //Z=13892
                    i = 2;  //Z=13893
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13894
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13895
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13896
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=13897
                        fkv[n] = fkv[n-1]*n;  //Z=13898
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13899
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13900
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13901 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13902
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13903
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13904
                        pn[n] = pn[n-1]*p*p;  //Z=13905
                        /*  longitudinal  */  //Z=13906
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13907
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13908
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13909

                            /* carr1fm[i]:=carr1pm[i];  //Z=13911 */
                            /* i:=i+1;  //Z=13912 */
                        }/*7*/  //Z=13913
                        /*  P(q)-coefficient  */  //Z=13914
                        /* carr1p[n]:=power(4,2*n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=13915 */
                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13916
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13917 */
                        /* carr1p[n]:=(1/2)*power(4,2*n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=13918 */

                        /*  F(q)-coefficient  */  //Z=13920
                        sump = 0.0;  //Z=13921
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13922
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=13923

                        /*  P(q)-coefficients  */  //Z=13925
                        /*  cross-sectional  */  //Z=13926
                        /*  F121  */  //Z=13927
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13928
                        /*  F122  */  //Z=13929
                        sump = 0.0;  //Z=13930
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13931
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13932
                        }/*7*/  //Z=13933
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=13934
                        /*  F123  */  //Z=13935
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=13936

                        /*  F(q)-coefficients  */  //Z=13938
                        /*  cross-sectional  */  //Z=13939
                        /*  F121  */  //Z=13940
                        sump = 0.0;  //Z=13941
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13942
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13943
                        }/*7*/  //Z=13944
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13945
                        /*  F122  */  //Z=13946
                        sump = 0.0;  //Z=13947
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13948
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13949
                        }/*7*/  //Z=13950
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=13951
                        /*  F123  */  //Z=13952
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=13953

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=13955
                            if ( n<n1 ) n1 = n;  //Z=13956
                        }/*7*/  //Z=13957
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=13958
                            if ( n<n4 ) n4 = n;  //Z=13959
                        }/*7*/  //Z=13960
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=13961
                            if ( n<n5 ) n5 = n;  //Z=13962
                        }/*7*/  //Z=13963
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=13964
                            if ( n<n6 ) n6 = n;  //Z=13965
                        }/*7*/  //Z=13966
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13967
                            if ( n<n1f ) n1f = n;  //Z=13968
                        }/*7*/  //Z=13969
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13970
                            if ( n<n4f ) n4f = n;  //Z=13971
                        }/*7*/  //Z=13972
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13973
                            if ( n<n5f ) n5f = n;  //Z=13974
                        }/*7*/  //Z=13975
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13976
                            if ( n<n6f ) n6f = n;  //Z=13977
                        }/*7*/  //Z=13978
                    }/*6*/  /*  of n-loop  */  //Z=13979
                }/*5*/  /*  cs=1  */  //Z=13980

                /*  inhomogeneous core/shell  */  //Z=13982
                if ( params.cs==2 )
                {/*5*/  //Z=13983
                    i = 2;  //Z=13984
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13985
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13986
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13987
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=13988
                        fkv[n] = fkv[n-1]*n;  //Z=13989
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13990
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13991
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13992 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13993
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13994
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13995
                        pn[n] = pn[n-1]*p*p;  //Z=13996
                        /*  longitudinal  */  //Z=13997
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13998
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13999
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=14000

                            /* carr1fm[i]:=carr1pm[i];  //Z=14002 */
                            /* i:=i+1;  //Z=14003 */
                        }/*7*/  //Z=14004
                        /*  P(q)-coefficient  */  //Z=14005
                        /* carr1p[n]:=power(4,2*n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=14006 */
                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14007
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14008 */
                        /* carr1p[n]:=(1/2)*power(4,2*n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=14009 */

                        /*  F(q)-coefficient  */  //Z=14011
                        sump = 0.0;  //Z=14012
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14013
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=14014

                        /*  cross-sectional P(q)  */  //Z=14016
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14017
                        sump = 0.0;  //Z=14018
                        sump1 = 0.0;  //Z=14019
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14020
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=14021
                            sump = sump+pn[n-m]*sumi;  //Z=14022
                            sump1 = sump1+sumi;  //Z=14023
                        }/*7*/  //Z=14024
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=14025
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n]*sump1;  //Z=14026
                        sump = 0.0;  //Z=14027
                        sump1 = 0.0;  //Z=14028
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14029
                            sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=14030
                            sump = sump+sumi;  //Z=14031
                            sump1 = sump1+pn[n-m]*sumi;  //Z=14032
                        }/*7*/  //Z=14033
                        params.CR->carr7p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=14034
                        params.CR->carr8p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=14035
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=14036

                        /*  F(q)-coefficients  */  //Z=14051
                        /*  cross-sectional  */  //Z=14052
                        params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=14053
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14054
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14055

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14057
                            if ( n<n1 ) n1 = n;  //Z=14058
                        }/*7*/  //Z=14059
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14060
                            if ( n<n4 ) n4 = n;  //Z=14061
                        }/*7*/  //Z=14062
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=14063
                            if ( n<n5 ) n5 = n;  //Z=14064
                        }/*7*/  //Z=14065
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=14066
                            if ( n<n6 ) n6 = n;  //Z=14067
                        }/*7*/  //Z=14068
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=14069
                            if ( n<n7 ) n7 = n;  //Z=14070
                        }/*7*/  //Z=14071
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=14072
                            if ( n<n8 ) n8 = n;  //Z=14073
                        }/*7*/  //Z=14074
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=14075
                            if ( n<n9 ) n9 = n;  //Z=14076
                        }/*7*/  //Z=14077
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14078
                            if ( n<n1f ) n1f = n;  //Z=14079
                        }/*7*/  //Z=14080
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14081
                            if ( n<n4f ) n4f = n;  //Z=14082
                        }/*7*/  //Z=14083
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14084
                            if ( n<n5f ) n5f = n;  //Z=14085
                        }/*7*/  //Z=14086
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14087
                            if ( n<n6f ) n6f = n;  //Z=14088
                        }/*7*/  //Z=14089
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=14090
                            if ( n<n7f ) n7f = n;  //Z=14091
                        }/*7*/  //Z=14092
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=14093
                            if ( n<n8f ) n8f = n;  //Z=14094
                        }/*7*/  //Z=14095
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=14096
                            if ( n<n9f ) n9f = n;  //Z=14097
                        }/*7*/  //Z=14098
                    }/*6*/  /*  of n-loop  */  //Z=14099
                }/*5*/  /*  cs=2  */  //Z=14100

                /*  myelin  */  //Z=14103
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=14104
                    i = 2;  //Z=14105
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14106
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14107
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14108
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=14109
                        fkv[n] = fkv[n-1]*n;  //Z=14110
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14111
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14112
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14113 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14114
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=14115
                        /*  longitudinal  */  //Z=14116
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14117
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14118
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=14119

                            /* carr1fm[i]:=carr1pm[i];  //Z=14121 */
                            /* i:=i+1;  //Z=14122 */
                        }/*7*/  //Z=14123
                        /*  P(q)-coefficient  */  //Z=14124
                        /* carr1p[n]:=power(4,2*n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));  //Z=14125 */
                        params.CR->carr1p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14126
                        /* carr1p[n]:=(sqrt(pi)/2)*fk2v[n]*z12vl[n]*xln[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14127 */
                        /* carr1p[n]:=(1/2)*power(4,2*n)*z12vl[n]*xln[n]/((n+1)*(n+1/2));  //Z=14128 */

                        /*  F(q)-coefficient  */  //Z=14130
                        sump = 0.0;  //Z=14131
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14132
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=14133

                        /*  cross-section  */  //Z=14135
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=14136

                        sump = 0.0;  //Z=14138
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14139
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=14140

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14142
                            if ( n<n1 ) n1 = n;  //Z=14143
                        }/*7*/  //Z=14144
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14145
                            if ( n<n4 ) n4 = n;  //Z=14146
                        }/*7*/  //Z=14147
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14148
                            if ( n<n1f ) n1f = n;  //Z=14149
                        }/*7*/  //Z=14150
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14151
                            if ( n<n4f ) n4f = n;  //Z=14152
                        }/*7*/  //Z=14153
                    }/*6*/  /*  of n-loop  */  //Z=14154
                }/*5*/  /*  cs=3  */  //Z=14155

            }/*4*/  /*  of cylinders  */  //Z=14157


            /* ** disks ** */  //Z=14160
            if ( dim==2 )
            {/*4*/  //Z=14161
                /*  homogeneous  */  //Z=14162
                if ( params.cs==0 )
                {/*5*/  //Z=14163
                    i = 2;  //Z=14164
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14165
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14166
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14167
                        fkv[n] = fkv[n-1]*n;  //Z=14168
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14169
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14170
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14171 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14172
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14173
                        /*  longitudinal  */  //Z=14174
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14175
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14176
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=14177 */
                            params.CR->carr22pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=14178
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=14179 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=14180
                            /* carr2fm[i]:=carr2pm[i];  //Z=14181 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=14182 */
                            /* i:=i+1;  //Z=14183 */
                        }/*7*/  //Z=14184
                        /*  P(q)-coefficient  */  //Z=14185
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14186
                        /*  F(q)-coefficient  */  //Z=14187
                        sump = 0.0;  //Z=14188
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14189
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14190
                        }/*7*/  //Z=14191
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=14192

                        /*  cross-section  */  //Z=14194
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14195
                        sump = 0.0;  //Z=14196
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14197
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14198
                        }/*7*/  //Z=14199
                        params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=14200

                        /*  series for <...> integration  */  //Z=14202
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14203
                        sump = 0.0;  //Z=14204
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14205
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14206
                        }/*7*/  //Z=14207
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14208

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14210
                            if ( n<n1 ) n1 = n;  //Z=14211
                        }/*7*/  //Z=14212
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14213
                            if ( n<n4 ) n4 = n;  //Z=14214
                        }/*7*/  //Z=14215
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14216
                            if ( n<n1f ) n1f = n;  //Z=14217
                        }/*7*/  //Z=14218
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14219
                            if ( n<n4f ) n4f = n;  //Z=14220
                        }/*7*/  //Z=14221
                    }/*6*/  /*  of n-loop  */  //Z=14222
                }/*5*/  /*  of cs=0  */  //Z=14223

                /*  core/shell  */  //Z=14225
                if ( params.cs==1 )
                {/*5*/  //Z=14226
                    i = 2;  //Z=14227
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14228
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14229
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14230
                        fkv[n] = fkv[n-1]*n;  //Z=14231
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14232
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14233
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14234 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14235
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14236
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14237
                        pn[n] = pn[n-1]*p*p;  //Z=14238
                        /*  longitudinal  */  //Z=14239
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14240
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14241
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=14242 */
                            params.CR->carr22pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=14243
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=14244 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=14245
                            /* carr2fm[i]:=carr2pm[i];  //Z=14246 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=14247 */
                            /* i:=i+1;  //Z=14248 */
                        }/*7*/  //Z=14249
                        /*  P(q)-coefficient  */  //Z=14250
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14251

                        /*  F(q)-coefficient  */  //Z=14253
                        sump = 0.0;  //Z=14254
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14255
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14256
                        }/*7*/  //Z=14257
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=14258

                        /*  cross-sectional  */  //Z=14260
                        /*  F121  */  //Z=14261
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14262
                        /*  F122  */  //Z=14263
                        sump = 0.0;  //Z=14264
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14265
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14266
                        }/*7*/  //Z=14267
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=14268
                        /*  F122  */  //Z=14269
                        sump = 0.0;  //Z=14270
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14271
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14272
                        }/*7*/  //Z=14273
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=14274
                        /*  F123  */  //Z=14275
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=14276
                        /*  F121  */  //Z=14277
                        sump = 0.0;  //Z=14278
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14279
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14280
                        }/*7*/  //Z=14281
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=14282
                        /*  F122  */  //Z=14283
                        sump = 0.0;  //Z=14284
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14285
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14286
                        }/*7*/  //Z=14287
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=14288
                        /*  F123  */  //Z=14289
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=14290

                        /*  series for <...> integration  */  //Z=14292
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14293
                        sump = 0.0;  //Z=14294
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14295
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14296
                        }/*7*/  //Z=14297
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14298

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14300
                            if ( n<n1 ) n1 = n;  //Z=14301
                        }/*7*/  //Z=14302
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14303
                            if ( n<n4 ) n4 = n;  //Z=14304
                        }/*7*/  //Z=14305
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=14306
                            if ( n<n5 ) n5 = n;  //Z=14307
                        }/*7*/  //Z=14308
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=14309
                            if ( n<n6 ) n6 = n;  //Z=14310
                        }/*7*/  //Z=14311
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14312
                            if ( n<n1f ) n1f = n;  //Z=14313
                        }/*7*/  //Z=14314
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14315
                            if ( n<n4f ) n4f = n;  //Z=14316
                        }/*7*/  //Z=14317
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14318
                            if ( n<n5f ) n5f = n;  //Z=14319
                        }/*7*/  //Z=14320
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14321
                            if ( n<n6f ) n6f = n;  //Z=14322
                        }/*7*/  //Z=14323
                    }/*6*/  /*  of n-loop  */  //Z=14324
                }/*5*/  /*  of cs=1  */  //Z=14325

                /*  inhomogeneous core/shell  */  //Z=14327
                if ( params.cs==2 )
                {/*5*/  //Z=14328
                    i = 2;  //Z=14329
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14330
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14331
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14332
                        fkv[n] = fkv[n-1]*n;  //Z=14333
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14334
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14335
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14336 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14337
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14338
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14339
                        pn[n] = pn[n-1]*p*p;  //Z=14340
                        /*  longitudinal  */  //Z=14341
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14342
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14343
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=14344 */
                            params.CR->carr22pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=14345
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=14346 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=14347
                            /* carr2fm[i]:=carr2pm[i];  //Z=14348 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=14349 */
                            /* i:=i+1;  //Z=14350 */
                        }/*7*/  //Z=14351
                        /*  P(q)-coefficient  */  //Z=14352
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14353

                        /*  F(q)-coefficient  */  //Z=14355
                        sump = 0.0;  //Z=14356
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14357
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14358
                        }/*7*/  //Z=14359
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=14360

                        /*  cross-sectional P(q)  */  //Z=14362
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14363
                        sump = 0.0;  //Z=14364
                        sump1 = 0.0;  //Z=14365
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14366
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=14367
                            sump = sump+pn[n-m]*sumi;  //Z=14368
                            sump1 = sump1+sumi;  //Z=14369
                        }/*7*/  //Z=14370
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14371
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=14372
                        sump = 0.0;  //Z=14373
                        sump1 = 0.0;  //Z=14374
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14375
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=14376
                            sump = sump+sumi;  //Z=14377
                            sump1 = sump1+pn[n-m]*sumi;  //Z=14378
                        }/*7*/  //Z=14379
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14380
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=14381
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=14382

                        /*  F(q)  */  //Z=14401
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=14402
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14403
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14404

                        /*  series for <...> integration  */  //Z=14406
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14407
                        sump = 0.0;  //Z=14408
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14409
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14410
                        }/*7*/  //Z=14411
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14412

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14414
                            if ( n<n1 ) n1 = n;  //Z=14415
                        }/*7*/  //Z=14416
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14417
                            if ( n<n4 ) n4 = n;  //Z=14418
                        }/*7*/  //Z=14419
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=14420
                            if ( n<n5 ) n5 = n;  //Z=14421
                        }/*7*/  //Z=14422
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=14423
                            if ( n<n6 ) n6 = n;  //Z=14424
                        }/*7*/  //Z=14425
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=14426
                            if ( n<n7 ) n7 = n;  //Z=14427
                        }/*7*/  //Z=14428
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=14429
                            if ( n<n8 ) n8 = n;  //Z=14430
                        }/*7*/  //Z=14431
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=14432
                            if ( n<n9 ) n9 = n;  //Z=14433
                        }/*7*/  //Z=14434
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14435
                            if ( n<n1f ) n1f = n;  //Z=14436
                        }/*7*/  //Z=14437
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14438
                            if ( n<n4f ) n4f = n;  //Z=14439
                        }/*7*/  //Z=14440
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14441
                            if ( n<n5f ) n5f = n;  //Z=14442
                        }/*7*/  //Z=14443
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14444
                            if ( n<n6f ) n6f = n;  //Z=14445
                        }/*7*/  //Z=14446
                    }/*6*/  /*  of n-loop  */  //Z=14447
                }/*5*/  /*  of cs=2  */  //Z=14448

            }/*4*/  /*  of disks  */  //Z=14450
        }/*3*/  //Z=14451

        /* ** z-axis ** */  //Z=14453
        if ( (params.orcase==4) && (dim!=3) )
        {/*3*/  //Z=14454
            /* ** cylinders ** */  //Z=14455
            if ( dim==1 )
            {/*4*/  //Z=14456
                /*  homogeneous  */  //Z=14457
                if ( params.cs==0 )
                {/*5*/  //Z=14458
                    i = 2;  //Z=14459
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14460
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14461
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14462
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=14463
                        fkv[n] = fkv[n-1]*n;  //Z=14464
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14465
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14466
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14467 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14468
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14469
                        /*  longitudinal  */  //Z=14470
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14471
                        /*  P(q)-coefficient  */  //Z=14472
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14473
                        /*  F(q)-coefficient  */  //Z=14474
                        sump = 0.0;  //Z=14475
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14476
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n)*fkv[n]*fkv[n]*params.norm);  //Z=14477
                        /*  cross-sectional  */  //Z=14478
                        /*  P(q)-coefficient  */  //Z=14479
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=14480
                        /*  F(q)-coefficient  */  //Z=14481
                        sump = 0.0;  //Z=14482
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14483
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=14484

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14486
                            if ( n<n1 ) n1 = n;  //Z=14487
                        }/*7*/  //Z=14488
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14489
                            if ( n<n4 ) n4 = n;  //Z=14490
                        }/*7*/  //Z=14491
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14492
                            if ( n<n1f ) n1f = n;  //Z=14493
                        }/*7*/  //Z=14494
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14495
                            if ( n<n4f ) n4f = n;  //Z=14496
                        }/*7*/  //Z=14497
                    }/*6*/  /*  of n-loop  */  //Z=14498
                }/*5*/  /*  of cs=0  */  //Z=14499

                /*  core/shell  */  //Z=14501
                if ( params.cs==1 )
                {/*5*/  //Z=14502
                    i = 2;  //Z=14503
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14504
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14505
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14506
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=14507
                        fkv[n] = fkv[n-1]*n;  //Z=14508
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14509
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14510
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14511 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14512
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14513
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14514
                        pn[n] = pn[n-1]*p*p;  //Z=14515
                        /*  longitudinal  */  //Z=14516
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14517
                        /*  P(q)-coefficient  */  //Z=14518
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14519
                        /*  F(q)-coefficient  */  //Z=14520
                        sump = 0.0;  //Z=14521
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14522
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14523
                        /*  P(q)-coefficients  */  //Z=14524
                        /*  cross-sectional  */  //Z=14525
                        /*  F121  */  //Z=14526
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14527
                        /*  F122  */  //Z=14528
                        sump = 0.0;  //Z=14529
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14530
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14531
                        }/*7*/  //Z=14532
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=14533
                        /*  F123  */  //Z=14534
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=14535

                        /*  F(q)-coefficients  */  //Z=14537
                        /*  cross-sectional  */  //Z=14538
                        /*  F121  */  //Z=14539
                        sump = 0.0;  //Z=14540
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14541
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14542
                        }/*7*/  //Z=14543
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=14544
                        /*  F122  */  //Z=14545
                        sump = 0.0;  //Z=14546
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14547
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14548
                        }/*7*/  //Z=14549
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=14550
                        /*  F123  */  //Z=14551
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=14552

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14554
                            if ( n<n1 ) n1 = n;  //Z=14555
                        }/*7*/  //Z=14556
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14557
                            if ( n<n4 ) n4 = n;  //Z=14558
                        }/*7*/  //Z=14559
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=14560
                            if ( n<n5 ) n5 = n;  //Z=14561
                        }/*7*/  //Z=14562
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=14563
                            if ( n<n6 ) n6 = n;  //Z=14564
                        }/*7*/  //Z=14565
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14566
                            if ( n<n1f ) n1f = n;  //Z=14567
                        }/*7*/  //Z=14568
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14569
                            if ( n<n4f ) n4f = n;  //Z=14570
                        }/*7*/  //Z=14571
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14572
                            if ( n<n5f ) n5f = n;  //Z=14573
                        }/*7*/  //Z=14574
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14575
                            if ( n<n6f ) n6f = n;  //Z=14576
                        }/*7*/  //Z=14577
                    }/*6*/  /*  of n-loop  */  //Z=14578
                }/*5*/  /*  of cs=1  */  //Z=14579

                /*  inhomogeneous core/shell  */  //Z=14582
                if ( params.cs==2 )
                {/*5*/  //Z=14583
                    i = 2;  //Z=14584
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14585
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14586
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14587
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=14588
                        fkv[n] = fkv[n-1]*n;  //Z=14589
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14590
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14591
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14592 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14593
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14594
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14595
                        pn[n] = pn[n-1]*p*p;  //Z=14596
                        /*  longitudinal  */  //Z=14597
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14598
                        /*  P(q)-coefficient  */  //Z=14599
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14600
                        /*  F(q)-coefficient  */  //Z=14601
                        sump = 0.0;  //Z=14602
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14603
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14604

                        /*  cross-sectional P(q)  */  //Z=14606
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14607
                        sump = 0.0;  //Z=14608
                        sump1 = 0.0;  //Z=14609
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14610
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=14611
                            sump = sump+pn[n-m]*sumi;  //Z=14612
                            sump1 = sump1+sumi;  //Z=14613
                        }/*7*/  //Z=14614
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=14615
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n]*sump1;  //Z=14616
                        sump = 0.0;  //Z=14617
                        sump1 = 0.0;  //Z=14618
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14619
                            sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=14620
                            sump = sump+sumi;  //Z=14621
                            sump1 = sump1+pn[n-m]*sumi;  //Z=14622
                        }/*7*/  //Z=14623
                        params.CR->carr7p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=14624
                        params.CR->carr8p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=14625
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=14626

                        /*  F(q)-coefficients  */  //Z=14641
                        /*  cross-sectional  */  //Z=14642
                        params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=14643
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14644
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14645

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14647
                            if ( n<n1 ) n1 = n;  //Z=14648
                        }/*7*/  //Z=14649
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14650
                            if ( n<n4 ) n4 = n;  //Z=14651
                        }/*7*/  //Z=14652
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=14653
                            if ( n<n5 ) n5 = n;  //Z=14654
                        }/*7*/  //Z=14655
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=14656
                            if ( n<n6 ) n6 = n;  //Z=14657
                        }/*7*/  //Z=14658
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=14659
                            if ( n<n7 ) n7 = n;  //Z=14660
                        }/*7*/  //Z=14661
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=14662
                            if ( n<n8 ) n8 = n;  //Z=14663
                        }/*7*/  //Z=14664
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=14665
                            if ( n<n9 ) n9 = n;  //Z=14666
                        }/*7*/  //Z=14667
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14668
                            if ( n<n1f ) n1f = n;  //Z=14669
                        }/*7*/  //Z=14670
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14671
                            if ( n<n4f ) n4f = n;  //Z=14672
                        }/*7*/  //Z=14673
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14674
                            if ( n<n5f ) n5f = n;  //Z=14675
                        }/*7*/  //Z=14676
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14677
                            if ( n<n6f ) n6f = n;  //Z=14678
                        }/*7*/  //Z=14679
                    }/*6*/  /*  of n-loop  */  //Z=14680
                }/*5*/  /*  of cs=2  */  //Z=14681

                /*  myelin  */  //Z=14683
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=14684
                    i = 2;  //Z=14685
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14686
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14687
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14688
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=14689
                        fkv[n] = fkv[n-1]*n;  //Z=14690
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14691
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14692
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14693 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14694
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=14695
                        /*  longitudinal  */  //Z=14696
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14697
                        /*  P(q)-coefficient  */  //Z=14698
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14699
                        /*  F(q)-coefficient  */  //Z=14700
                        sump = 0.0;  //Z=14701
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14702
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n)*fkv[n]*fkv[n]*params.norm);  //Z=14703
                        /*  cross-sectional  */  //Z=14704
                        /*  P(q)-coefficient  */  //Z=14705
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=14706

                        /*  F(q)-coefficient  */  //Z=14708
                        sump = 0.0;  //Z=14709
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14710
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=14711

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14713
                            if ( n<n1 ) n1 = n;  //Z=14714
                        }/*7*/  //Z=14715
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14716
                            if ( n<n4 ) n4 = n;  //Z=14717
                        }/*7*/  //Z=14718

                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14720
                            if ( n<n1f ) n1f = n;  //Z=14721
                        }/*7*/  //Z=14722
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14723
                            if ( n<n4f ) n4f = n;  //Z=14724
                        }/*7*/  //Z=14725
                    }/*6*/  /*  of n-loop  */  //Z=14726
                }/*5*/  /*  of cs=3  */  //Z=14727
            }/*4*/  /*  of cylinders  */  //Z=14728

            /* ** disks ** */  //Z=14730
            if ( dim==2 )
            {/*4*/  //Z=14731
                /*  homogeneous  */  //Z=14732
                if ( params.cs==0 )
                {/*5*/  //Z=14733
                    carr2i[0] = 1;  //Z=14734
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14735
                        fkv[n] = fkv[n-1]*n;  //Z=14736
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14737
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14738
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14739
                    }/*6*/  //Z=14740
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14741
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14742
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14743
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14744
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14745 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14746
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14747
                        /*  longitudinal  */  //Z=14748
                        sump = 0.0;  //Z=14749
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14750
                            sump = sump+carr2i[m]/fkv[n-m];  //Z=14751
                        }/*7*/  //Z=14752
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14753
                        sump1 = 0.0;  //Z=14754
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14755
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14756
                        }/*7*/  //Z=14757
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;  //Z=14758

                        /*  cross-sectional  */  //Z=14760
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14761
                        sump = 0.0;  //Z=14762
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14763
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14764
                        }/*7*/  //Z=14765
                        params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=14766

                        /*  series for <...> integration  */  //Z=14768
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14769
                        sump = 0.0;  //Z=14770
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14771
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14772
                        }/*7*/  //Z=14773
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14774

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14776
                            if ( n<n1 ) n1 = n;  //Z=14777
                        }/*7*/  //Z=14778
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14779
                            if ( n<n4 ) n4 = n;  //Z=14780
                        }/*7*/  //Z=14781
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14782
                            if ( n<n1f ) n1f = n;  //Z=14783
                        }/*7*/  //Z=14784
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14785
                            if ( n<n4f ) n4f = n;  //Z=14786
                        }/*7*/  //Z=14787
                    }/*6*/  /*  of n-loop  */  //Z=14788
                }/*5*/  /*  of cs=0  */  //Z=14789


                /*  core/shell  */  //Z=14792
                if ( params.cs==1 )
                {/*5*/  //Z=14793
                    carr2i[0] = 1;  //Z=14794
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14795
                        fkv[n] = fkv[n-1]*n;  //Z=14796
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14797
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14798
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14799
                    }/*6*/  //Z=14800
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14801
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14802
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14803
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14804
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14805 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14806
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14807
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14808
                        pn[n] = pn[n-1]*p*p;  //Z=14809
                        /*  longitudinal  */  //Z=14810
                        sump = 0.0;  //Z=14811
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14812
                            sump = sump+carr2i[m]/fkv[n-m];  //Z=14813
                        }/*7*/  //Z=14814
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14815
                        sump1 = 0.0;  //Z=14816
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14817
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14818
                        }/*7*/  //Z=14819
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;  //Z=14820

                        /*  cross-sectional  */  //Z=14822
                        /*  F121  */  //Z=14823
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14824
                        /*  F122  */  //Z=14825
                        sump = 0.0;  //Z=14826
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14827
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14828
                        }/*7*/  //Z=14829
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=14830
                        /*  F122  */  //Z=14831
                        sump = 0.0;  //Z=14832
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14833
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14834
                        }/*7*/  //Z=14835
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=14836
                        /*  F123  */  //Z=14837
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=14838
                        /*  F121  */  //Z=14839
                        sump = 0.0;  //Z=14840
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14841
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14842
                        }/*7*/  //Z=14843
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=14844
                        /*  F122  */  //Z=14845
                        sump = 0.0;  //Z=14846
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14847
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14848
                        }/*7*/  //Z=14849
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=14850
                        /*  F123  */  //Z=14851
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=14852

                        /*  series for <...> integration  */  //Z=14854
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14855
                        sump = 0.0;  //Z=14856
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14857
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14858
                        }/*7*/  //Z=14859
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14860

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14862
                            if ( n<n1 ) n1 = n;  //Z=14863
                        }/*7*/  //Z=14864
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14865
                            if ( n<n4 ) n4 = n;  //Z=14866
                        }/*7*/  //Z=14867
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=14868
                            if ( n<n5 ) n5 = n;  //Z=14869
                        }/*7*/  //Z=14870
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=14871
                            if ( n<n6 ) n6 = n;  //Z=14872
                        }/*7*/  //Z=14873
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14874
                            if ( n<n1f ) n1f = n;  //Z=14875
                        }/*7*/  //Z=14876
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14877
                            if ( n<n4f ) n4f = n;  //Z=14878
                        }/*7*/  //Z=14879
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14880
                            if ( n<n5f ) n5f = n;  //Z=14881
                        }/*7*/  //Z=14882
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14883
                            if ( n<n6f ) n6f = n;  //Z=14884
                        }/*7*/  //Z=14885
                    }/*6*/  /*  of n-loop  */  //Z=14886
                }/*5*/  /*  of cs=1  */  //Z=14887

                /*  inhomogeneous core/shell  */  //Z=14889
                if ( params.cs==2 )
                {/*5*/  //Z=14890
                    carr2i[0] = 1;  //Z=14891
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14892
                        fkv[n] = fkv[n-1]*n;  //Z=14893
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14894
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14895
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14896
                    }/*6*/  //Z=14897
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14898
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14899
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14900
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14901
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14902 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14903
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14904
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14905
                        pn[n] = pn[n-1]*p*p;  //Z=14906
                        /*  longitudinal  */  //Z=14907
                        sump = 0.0;  //Z=14908
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14909
                            sump = sump+carr2i[m]/fkv[n-m];  //Z=14910
                        }/*7*/  //Z=14911
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14912
                        sump1 = 0.0;  //Z=14913
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14914
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14915
                        }/*7*/  //Z=14916
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;  //Z=14917

                        /*  cross-sectional P(q)  */  //Z=14919
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14920
                        sump = 0.0;  //Z=14921
                        sump1 = 0.0;  //Z=14922
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14923
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=14924
                            sump = sump+pn[n-m]*sumi;  //Z=14925
                            sump1 = sump1+sumi;  //Z=14926
                        }/*7*/  //Z=14927
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14928
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=14929
                        sump = 0.0;  //Z=14930
                        sump1 = 0.0;  //Z=14931
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14932
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=14933
                            sump = sump+sumi;  //Z=14934
                            sump1 = sump1+pn[n-m]*sumi;  //Z=14935
                        }/*7*/  //Z=14936
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14937
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=14938
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=14939

                        /*  F(q)  */  //Z=14959
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=14960
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14961
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14962

                        /*  series for <...> integration  */  //Z=14964
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14965
                        sump = 0.0;  //Z=14966
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14967
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14968
                        }/*7*/  //Z=14969
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14970

                        if ( fabs(params.CR->carr1p[n])<min )
                        {/*7*/  //Z=14972
                            if ( n<n1 ) n1 = n;  //Z=14973
                        }/*7*/  //Z=14974
                        if ( fabs(params.CR->carr4p[n])<min )
                        {/*7*/  //Z=14975
                            if ( n<n4 ) n4 = n;  //Z=14976
                        }/*7*/  //Z=14977
                        if ( fabs(params.CR->carr5p[n])<min )
                        {/*7*/  //Z=14978
                            if ( n<n5 ) n5 = n;  //Z=14979
                        }/*7*/  //Z=14980
                        if ( fabs(params.CR->carr6p[n])<min )
                        {/*7*/  //Z=14981
                            if ( n<n6 ) n6 = n;  //Z=14982
                        }/*7*/  //Z=14983
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=14984
                            if ( n<n7 ) n7 = n;  //Z=14985
                        }/*7*/  //Z=14986
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=14987
                            if ( n<n8 ) n8 = n;  //Z=14988
                        }/*7*/  //Z=14989
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=14990
                            if ( n<n9 ) n9 = n;  //Z=14991
                        }/*7*/  //Z=14992
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14993
                            if ( n<n1f ) n1f = n;  //Z=14994
                        }/*7*/  //Z=14995
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14996
                            if ( n<n4f ) n4f = n;  //Z=14997
                        }/*7*/  //Z=14998
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14999
                            if ( n<n5f ) n5f = n;  //Z=15000
                        }/*7*/  //Z=15001
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=15002
                            if ( n<n6f ) n6f = n;  //Z=15003
                        }/*7*/  //Z=15004
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=15005
                            if ( n<n7f ) n7f = n;  //Z=15006
                        }/*7*/  //Z=15007
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=15008
                            if ( n<n8f ) n8f = n;  //Z=15009
                        }/*7*/  //Z=15010
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=15011
                            if ( n<n9f ) n9f = n;  //Z=15012
                        }/*7*/  //Z=15013
                    }/*6*/  /*  of n-loop  */  //Z=15014
                }/*5*/  /*  of cs=2  */  //Z=15015


            }/*4*/  /*  of disk  */  //Z=15018
        }/*3*/  /*  of z-axis  */  //Z=15019
    }/*2*/  /*  of ordis=0  */  //Z=15020


Label99:  //Z=15068
    // Es kommt vor, dass die carr??[] beim genutzten Index (n?) den Wert 'inf' haben. Damit
    // lässt sich aber schlecht weiterrechnen. Daher hier die passenden Indizes verringern auf
    // den letzten gültigen Wert.
#ifndef __CUDACC__
    //int n1sav=n1, n2sav=n2, n3sav=n3, n4sav=n4, n5sav=n5, n6sav=n6, n7sav=n7, n8sav=n8, n9sav=n9;
    //int n1fsav=n1f, n2fsav=n2f, n3fsav=n3f, n4fsav=n4f, n5fsav=n5f, n6fsav=n6f, n7fsav=n7f, n8fsav=n8f, n9fsav=n9f;
#endif
    while ( n1  > 1 && (isinf(params.CR->carr1p[n1 ]) || fabs(params.CR->carr1p[n1 ])<1e-99) ) n1--;
    while ( n2  > 1 && (isinf(params.CR->carr2p[n2 ]) || fabs(params.CR->carr2p[n2 ])<1e-99) ) n2--;
    while ( n3  > 1 && (isinf(params.CR->carr3p[n3 ]) || fabs(params.CR->carr3p[n3 ])<1e-99) ) n3--;
    while ( n4  > 1 && (isinf(params.CR->carr4p[n4 ]) || fabs(params.CR->carr4p[n4 ])<1e-99) ) n4--;
    while ( n5  > 1 && (isinf(params.CR->carr5p[n5 ]) || fabs(params.CR->carr5p[n5 ])<1e-99) ) n5--;
    while ( n6  > 1 && (isinf(params.CR->carr6p[n6 ]) || fabs(params.CR->carr6p[n6 ])<1e-99) ) n6--;
    while ( n7  > 1 && (isinf(params.CR->carr7p[n7 ]) || fabs(params.CR->carr7p[n7 ])<1e-99) ) n7--;
    while ( n8  > 1 && (isinf(params.CR->carr8p[n8 ]) || fabs(params.CR->carr8p[n8 ])<1e-99) ) n8--;
    while ( n9  > 1 && (isinf(params.CR->carr9p[n9 ]) || fabs(params.CR->carr9p[n9 ])<1e-99) ) n9--;
    while ( n1f > 1 && (isinf(params.CR->carr1f[n1f]) || fabs(params.CR->carr1f[n1f])<1e-99) ) n1f--;
    while ( n2f > 1 && (isinf(params.CR->carr2f[n2f]) || fabs(params.CR->carr2f[n2f])<1e-99) ) n2f--;
    while ( n3f > 1 && (isinf(params.CR->carr3f[n3f]) || fabs(params.CR->carr3f[n3f])<1e-99) ) n3f--;
    while ( n4f > 1 && (isinf(params.CR->carr4f[n4f]) || fabs(params.CR->carr4f[n4f])<1e-99) ) n4f--;
    while ( n5f > 1 && (isinf(params.CR->carr5f[n5f]) || fabs(params.CR->carr5f[n5f])<1e-99) ) n5f--;
    while ( n6f > 1 && (isinf(params.CR->carr6f[n6f]) || fabs(params.CR->carr6f[n6f])<1e-99) ) n6f--;
    while ( n7f > 1 && (isinf(params.CR->carr7f[n7f]) || fabs(params.CR->carr7f[n7f])<1e-99) ) n7f--;
    while ( n8f > 1 && (isinf(params.CR->carr8f[n8f]) || fabs(params.CR->carr8f[n8f])<1e-99) ) n8f--;
    while ( n9f > 1 && (isinf(params.CR->carr9f[n9f]) || fabs(params.CR->carr9f[n9f])<1e-99) ) n9f--;
    // Jetzt folgt erst die lim? Berechnung
    params.limq1 = pow(fabs(params.CR->carr1p[n1]),-1/(2.0*n1));  //Z=15070
    params.limq2 = pow(fabs(params.CR->carr2p[n2]),-1/(2.0*n2));  //Z=15071
    params.limq3 = pow(fabs(params.CR->carr3p[n3]),-1/(2.0*n3));  //Z=15072
    params.limq4 = pow(fabs(params.CR->carr4p[n4]),-1/(2.0*n4));  //Z=15073
    params.limq5 = pow(fabs(params.CR->carr5p[n5]),-1/(2.0*n5));  //Z=15074
    params.limq6 = pow(fabs(params.CR->carr6p[n6]),-1/(2.0*n6));  //Z=15075
    params.limq7 = pow(fabs(params.CR->carr7p[n7]),-1/(2.0*n7));  //Z=15076
    params.limq8 = pow(fabs(params.CR->carr8p[n8]),-1/(2.0*n8));  //Z=15077
    params.limq9 = pow(fabs(params.CR->carr9p[n9]),-1/(2.0*n9));  //Z=15078
    params.limq1f = pow(fabs(params.CR->carr1f[n1f]),-1/(2.0*n1f));  //Z=15079 Im Orginal-Pascalprogramm wurden hier
    params.limq2f = pow(fabs(params.CR->carr2f[n2f]),-1/(2.0*n2f));  //Z=15080  auch die Variablen n1 bis n9 genutzt
    params.limq3f = pow(fabs(params.CR->carr3f[n3f]),-1/(2.0*n3f));  //Z=15081  und nicht die jetzt hier stehenden
    params.limq4f = pow(fabs(params.CR->carr4f[n4f]),-1/(2.0*n4f));  //Z=15082  n1f bis n9f, obwohl diese oben passend
    params.limq5f = pow(fabs(params.CR->carr5f[n5f]),-1/(2.0*n5f));  //Z=15083  bestimmt worden sind...
    params.limq6f = pow(fabs(params.CR->carr6f[n6f]),-1/(2.0*n6f));  //Z=15084
    params.limq7f = pow(fabs(params.CR->carr7f[n7f]),-1/(2.0*n7f));  //Z=15085  TODO ?
    params.limq8f = pow(fabs(params.CR->carr8f[n8f]),-1/(2.0*n8f));  //Z=15086
    params.limq9f = pow(fabs(params.CR->carr9f[n9f]),-1/(2.0*n9f));  //Z=15087
#ifndef __CUDACC__
    //qDebug() << "              " << "limq4"<<params.limq4 << "limq4f"<<params.limq4f;
#endif

#ifndef __CUDACC__
/*
#define CI(a1,n1,a2,n2,l) qPrintable(QString("%1:%2  %3:%4  %5").arg(n1,3).arg(a1[n1],12,'e',5).arg(n2,3).arg(a2[n2],12,'e',5).arg(l,8,'f',6))
    qDebug() << "Label99: ?p:  n?:carr?p[n?]  n?sav:carr?p[n?sav]  lim?    | ?f: n?f:carr?f[n?f] n?fsav:carr?f[n?fsav]  lim?f";
    qDebug() << "        " << "1p:" << CI(params.CR->carr1p,n1, params.CR->carr1p,n1sav, params.limq1) << "| 1f:" << CI(params.CR->carr1f,n1f, params.CR->carr1f,n1fsav, params.limq1f);
    qDebug() << "        " << "2p:" << CI(params.CR->carr2p,n2, params.CR->carr2p,n2sav, params.limq2) << "| 2f:" << CI(params.CR->carr2f,n2f, params.CR->carr2f,n2fsav, params.limq2f);
    qDebug() << "        " << "3p:" << CI(params.CR->carr3p,n3, params.CR->carr3p,n3sav, params.limq3) << "| 3f:" << CI(params.CR->carr3f,n3f, params.CR->carr3f,n3fsav, params.limq3f);
    qDebug() << "        " << "4p:" << CI(params.CR->carr4p,n4, params.CR->carr4p,n4sav, params.limq4) << "| 4f:" << CI(params.CR->carr4f,n4f, params.CR->carr4f,n4fsav, params.limq4f);
    qDebug() << "        " << "5p:" << CI(params.CR->carr5p,n5, params.CR->carr5p,n5sav, params.limq5) << "| 5f:" << CI(params.CR->carr5f,n5f, params.CR->carr5f,n5fsav, params.limq5f);
    qDebug() << "        " << "6p:" << CI(params.CR->carr6p,n6, params.CR->carr6p,n6sav, params.limq6) << "| 6f:" << CI(params.CR->carr6f,n6f, params.CR->carr6f,n6fsav, params.limq6f);
    qDebug() << "        " << "7p:" << CI(params.CR->carr7p,n7, params.CR->carr7p,n7sav, params.limq7) << "| 7f:" << CI(params.CR->carr7f,n7f, params.CR->carr7f,n7fsav, params.limq7f);
    qDebug() << "        " << "8p:" << CI(params.CR->carr8p,n8, params.CR->carr8p,n8sav, params.limq8) << "| 8f:" << CI(params.CR->carr8f,n8f, params.CR->carr8f,n8fsav, params.limq8f);
    qDebug() << "        " << "9p:" << CI(params.CR->carr9p,n9, params.CR->carr9p,n9sav, params.limq9) << "| 9f:" << CI(params.CR->carr9f,n9f, params.CR->carr9f,n9fsav, params.limq9f);
#undef CI
*/
/*
    QFile fout(QString("carr11pm_ordis%1_dim%2_cs%3.csv").arg(ordis).arg(dim).arg(cs));
    if ( fout.open(QIODevice::WriteOnly) )
    {
        QString str = "carr11pm[n][m]";
        for ( int m=0; m<imax2d_len; m++ )
            str += QString("; m=%1").arg(m,3);
        fout.write(qPrintable(str+EOL));
        for ( int n=0; n<imax2d_len && n<=120; n++ ) // imax2d_len=130+1
        {   str = QString("n=%1").arg(n,3);
            for ( int m=0; m<imax2d_len; m++ )
                str += QString("; %1").arg(params.CR->carr11pm[n][m]);
            fout.write(qPrintable(str+EOL));
        }
        //---
        str = "Array";
        for ( int n=0; n<coeffarray_len && n<=120; n++ )
            str += QString("; n=%1").arg(n,3);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr1p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr1p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr2p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr2p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr3p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr3p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr4p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr4p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr5p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr5p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr6p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr6p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr7p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr7p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr8p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr8p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr9p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(params.CR->carr9p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        fout.close();
    }
    else
        qDebug() << fout.errorString();
*/
#endif
}
