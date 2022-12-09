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


#define CLASSLIB SasCalc_GENERIC_calculation
//#define USE_szave
#define USE_psphere // nur init
//#define USE_psphered
//#define USE_pspheredf
//#define USE_f2dschulz
//#define USE_f2dhyper
#define USE_gamma
//#define USE_cosav
//#define USE_gammaratio
//#define USE_pqcoreshell
//#define USE_f2dcoreshell
//#define USE_polyvesicle
//#define USE_f2dpolyvesicle
//#define USE_polyliposome
//#define USE_polycube
//#define USE_burger
//#define USE_angleuv
#define USE_lorentznorm3
#define USE_gaussnorm3
#define USE_pearsonnorm3
#define USE_pearsonintegral3 // TODO raus?
//#define USE_cosavm
//#define USE_pgensphere1
//#define USE_hypcoreshellm
//#define USE_trapcube
//#define USE_cubealf
//#define USE_trapcubebet
//#define USE_cube
//#define USE_trapezpqr
//#define USE_midpntchi
//#define USE_qrombpq
//#define USE_qrombpi2
//#define USE_bessel10
#define USE_ButtonHKLClick
#define USE_fhkl_c
#define USE_extinction
#define USE_trapzddeltac
#define USE_polint
#define USE_qrombdeltac
#define USE_qrombchid
#define USE_trapzdchid

#define USE_CLASSLIB
#include "sc_libs_gpu.cu"
#include "sc_memory_gpu.cu"


#define GPU_2D
// If defined, use block.x and block.y etc. (new version)
// Otherwise use only block.x and iterate in each thread over y (old version)


SasCalc_GENERIC_calculation *SasCalc_GENERIC_calculation::inst;              //!< class instance pointer for the threads


SasCalc_GENERIC_calculation::SasCalc_GENERIC_calculation()
{
    inst = this;
    progAndAbort = nullptr;
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
//__launch_bounds_(256, 1)
__global__ void doIntCalc_GENERIC_kernel( SasCalc_GENERIC_calculation FCC )
{
    int ihex = blockIdx.x * blockDim.x + threadIdx.x + FCC.zzmin;
    if ( ihex >= FCC.zzmax ) return;
#ifdef GPU_2D
    int i    = blockIdx.y * blockDim.y + threadIdx.y + FCC.iimin;
    if ( i >= FCC.iimax ) return;
    SasCalc_GENERIC_calculation::doIntCalc_GENERIC_F( FCC, ihex, i );
#else
    for ( int i=threadIdx.x+FCC.iimin; i<FCC.iimax; i+=blockDim.x )
        SasCalc_GENERIC_calculation::doIntCalc_GENERIC_F( FCC, ihex, i );
#endif
}

__global__ void doIntFitCalc_GENERIC_kernel( SasCalc_GENERIC_calculation FCC )
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if ( x >= FCC._fitWidth ) return;
#ifdef GPU_2D
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if ( y >= FCC._fitHeight ) return;
    SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC_F( FCC, x, y );
#else
    for ( int y=threadIdx.x; y<FCC._fitHeight; i+=blockDim.x )
        SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC_F( FCC, x, y );
#endif
}
#endif


bool SasCalc_GENERIC_calculation::prepareCalculation()
{
    dbgCount=0;  // Zum Test...

    /*Z=19733 (******** main program for calculation *****) */

    if ( ltype == 14 /* GiSAXS-pattern */ ) calcQuadrants = radQ4;  /*Z=21005*/
    if ( ltype == 15 /* GiSAXS-pattern */ ) calcQuadrants = radQ4;
    // RadioButtonSlice wird hier noch ignoriert

    zzmin = calcQuadrants == radQ4 ? -zmax : 0;     /*Z=20987 im Batch-Loop */
    zzmax = zmax;
    iimin = calcQuadrants == radQ1 ? 0 : -zmax;
    iimax = zmax;

    //if ( ! tpvRandomOldValues.empty() )   // prepareCalculation()
    //{
    //    std::cerr << "TPV PREPARE CALCULATION init " << params.width_zuf << std::endl;
    //}
    //else
    //    std::cerr << "TPV prepare width = " << params.width_zuf << std::endl;

    //qDebug() << "PrepCalc" << zzmin << zzmax << iimin << iimax;

    ax1n_sigx = ax1.length() * sig.x();
    ax2n_sigy = ax2.length() * sig.y();
    ax3n_sigz = ax3.length() * sig.z();

    lat1d = false;  /*Z=19783*/
    lat2d = false;
    lat3d = true;
    if ( ltype == 0 ) { lat1d=true; lat3d=false; }     /*Z=19819*/
    if ( ltype == 1 || ltype == 2 || ltype == 3 ) { lat2d=true; lat3d=false; }     /*Z=19823*/

    if ( CheckBoxWAXS )  /*Z=20000  && RadioButtonPreset */
    {
        double thetmax = atan2( det, pixnox * pixx / 2.0 );
        wave = 2. * M_PI * sin(thetmax) / qmax;
    }

    c0=1;       /*Z=20097*/
    c1=4;
    c2=4*(sqrt(2.)-1);
    c3=4*(exp(2*log(2.)/3.)-1);
    if ( beta == 0 )
        c4 = 1; // da bei "case cbpeakPearsonVII /*Pearson*/:" damit multipliziert wird
    else
        c4=exp(log(2.)/beta)-1;

    //aziwidth:=StrToFloat(EditAzi.Text);
    phiwidth=4./aziwidth;   /*Z=19867*/
    // critdist:=StrToFloat(EditCritDist.Text);
    //critdist=0.5;  Echter Input      //(* for allowed/non-allowed reflection overview *)

    //(*** this generates reflection tables latpar1, latpar2, latpar1a,latpar2a ***)
    params.amax = 10; // TODO, noch kein GUI-Element vorhanden
    params.bmax = 10;
    params.cmax = 10;

    /*Z=19757*/
    if ( latpar1ptr == nullptr )
    {
        createMemory( (void **)(&latpar1ptr), sizeof(int) * latparlen * 6, latpar1Size, true, "latpar1" );
        createMemory( (void **)(&latpar2ptr), sizeof(int) * latparlen * 6, latpar2Size, true, "latpar2" );
        createMemory( (void **)(&latpar3ptr), sizeof(float) * latparlen * 14, latpar3Size, true, "latpar3" );
        //memory->createMemory( static_cast<void **>(&latpar4ptr), sizeof(float) * SC_Libs::latparlen * 2 );
    }
    if ( params.CR == nullptr )
    {
//        carr1p,carr2p,carr3p,carr4p,carr5p,carr6p,carr7p,carr8p,carr9p,
//        carr1f,carr2f,carr3f,carr4f,carr5f,carr6f,carr7f,carr8f,carr9f: ^CoeffArrayType;
//        coeffarraytype=array[0..150] of extended;
        // besser einen Speicher als Record anlegen
        createMemory( (void **)(&params.CR), sizeof(_carrXX), arrCRSize, true, "CR" );
    }

    ucphi=phi;      /*Z=19904*/
    uctheta=theta;
    cosphi=cos(phi*M_PI/180.);  /*Z=19908*/     /*TPV*/
    sinphi=sin(phi*M_PI/180.);                  /*TPV*/
    costheta=cos(theta*M_PI/180.);
    sintheta=sin(theta*M_PI/180.);

    if ( ltype != 12 )
    {   /*Z=20185 aber auch in Z=20483 ... */
        ButtonHKLClick( ltype, latpar1ptr, latpar2ptr );    // ltype=0..11 und nicht mehr ....
        peakmax1 = latpar1(1,0);
        peakmax2 = latpar2(1,0);    /*Z0311=19753*/
        //qDebug() << "nach ButtonHKLClick" << peakmax1 << peakmax2;
    }
    else
    {   // Zur Sicherheit, nicht im Pasccal-Programm
        peakmax1=0;
        peakmax2=0;
    }

    //(*** Peak Shapes ***)
    // Die Werte von <shp> und <shf> werden von der GUI übernommen.

    //(* Particle Interior *)  /*Z0311=19828*/
    cs = ComboBoxInterior;
    //if ( ComboBoxInterior == cbintHomogeneous/*0*/ ) cs=0;
    //if ( ComboBoxInterior == cbintCoreInShell/*2*/ ) cs=2;
    //if ( ComboBoxInterior == cbintMyelin/*3*/      ) cs=3;  //Neu
    if ( ComboBoxInterior == cbintCoreShell/*1*/   ) { cs=1; alphash=0.0001; }

    //(* Orientational distribution functions *) /*Z0311=19843*/
    // Der Wert von <ordis> wird von der GUI übernommen.

    /*Z0311=19860*/
    //rotx:=StrToFloat(EditRotx.Text);     (* rotation axis *)
    //roty:=StrToFloat(EditRoty.Text);
    //rotz:=StrToFloat(EditRotz.Text);
    //xx//rotnorm=sqrt( sqr(rot.x()) + sqr(rot.y()) + sqr(rot.z()) );
    //rota:=StrToFloat(EditRotAlpha.Text);    (* rotation angle *)
    //xx//rota0=rota;
    //adet:=StrToFloat(EditDetector.Text);    (* tilt angle *)
    //xcur:=StrToFloat(Editx.Text);           (* qx of cursor *)
    //ycur:=StrToFloat(Edity.Text);           (* qy of cursor *)
    //xycur=sqrt(xcur*xcur+ycur*ycur); - sofort beim setzen von xcur,ycur
    //anglecur=2*M_PI*StrToFloat(EditAnglexy.Text)/360;    (* angle of cursor *)


    /*Z=20221 Particle type */
    /*Z=20435 für die cdim=?, im Batch-Loop */
    switch ( ComboBoxParticle )
    {
    case cbpartSphere/*0*/:      //(* sphere *)
        cdim=3;
        part=0;
        //partsphere=true;
        partdim=3;
        break;
    case cbpartCylinder/*1*/:    //(* cylinder *)
        cdim=1;
        part=1;
        //partcylinder=true;
        partdim=2;
        break;
    case cbpartDisk/*2*/:        //(* disk *)
        cdim=2;
        part=2;
        //partdisk=true;
        partdim=1;
        break;
    case cbpartVesicle/*3*/:     //(* vesicle *)
        cdim=3;
        part=3;
        //partvesicle:=true;
        partdim=3;
        break;
    case cbpartCube/*4*/:        //(* cube *)
        cdim=4;
        part=4;
        //partcube:=true;
        partdim=4;
        break;
    case cbpartEllipsoide/*5*/:  //(* ellipsoid *)
        cdim=5;
        part=5;
        //partellipsoid:=true;
        partdim=5;
        break;
    case cbpartTriaxEllips/*6*/: //(* triaxial ellipsoid *)
        cdim=6;
        part=6;
        //parttriellipsoid:=true;
        partdim=6;
        break;
    case cbpartSuperEllips/*7*/: //(* super ellipsoid, barrel *)
        cdim=7;
        part=7;
        //partbarrel:=true;
        partdim=7;
        break;
        /*Z=20443 TODO neue Werte bei cdim, sicher auch vorher für die anderen */
    } // switch ( ComboBoxParticle )

    /*Z=19880 Parameter array for liposome and myelin */
    params.CR->myarray[0]=params.length;     /* axon length */
    params.CR->myarray[1]=params.radius;      /* axon radius */
    params.CR->myarray[2]=params.sigma;
    params.CR->myarray[3]=params.sigmal;
    params.CR->myarray[4]=params.radiusi;     /* no. of layers */
    params.CR->myarray[5]=alphash;            /* l_ex */
    params.CR->myarray[6]=params.rho;         /* l_in */
    params.CR->myarray[7]=0; // TODO: acpl;        /* l_lip_head */
    params.CR->myarray[8]=0; // TODO: bcpl;        /* l_lip_tail */
    params.CR->myarray[9]=params.uca;         /* phi_axon */
    params.CR->myarray[10]=params.ucb;        /* phi_intra */
    params.CR->myarray[11]=params.ucc;        /* phi_extra */
    params.CR->myarray[12]=params.domainsize; /* phi_head */
    params.CR->myarray[13]=aziwidth;          /* phi_tail */
    params.CR->myarray[14]=1;                 /* inmax */
    params.CR->myarray[15]=1;                 /* vv */
    params.CR->myarray[16]=1;                 /* rmax */
    params.CR->myarray[17]=iso;               /* test parameter */

    if ( cdim != 9 && cdim != 10 )  // Neue Werte für cdim (s.o.)
    {   /*Z=20452*/
        coefficients( params.length, params.radius, params.radiusi, params.sigmal, params.sigma, phi, alphash, params.dbeta,
                     theta, phi, part, cdim, 120, ordis, cs, params.CR->myarray,
                     /* ab hier outputs */
                     orcase, por, order, norm,
                     limq1,  limq2,  limq3,  limq4,  limq5,  limq6,  limq7,  limq8,  limq9,
                     limq1f, limq2f, limq3f, limq4f, limq5f, limq6f, limq7f, limq8f, limq9f, // TODO-F
                     params.CR );  // Für GPU ist ein Zeiger besser als viele!
    }

/*#ifndef __CUDACC__
    for ( int i=0; i<coeffarray_len; i++ )
    {
        //qDebug() << "carr1p[]" << i << params.CR->carr1p[i];
        if ( isnan(params.CR->carr1p[i]) )
        {
            if ( i > 0 )
                qDebug() << "carr1p["<<i-1<<"] =" << params.CR->carr1p[i-1];
            qDebug() << "carr1p["<<i<<"] = NaN";
            break;
        }
        if ( isinf(params.CR->carr1p[i]) )
        {
            if ( i > 0 )
                qDebug() << "carr1p["<<i-1<<"] =" << params.CR->carr1p[i-1];
            qDebug() << "carr1p["<<i<<"] =" << params.CR->carr1p[i];
            break;
        }
        if ( i > 1 && params.CR->carr1p[i] == 1 ) break;
    }
#endif*/
    /* Ändert nichts am Ergebnis ...
    for ( int i=0; i<coeffarray_len; i++ )
    {
        if ( isnan(params.CR->carr1p[i]) || isinf(params.CR->carr1p[i]) ) params.CR->carr1p[i] = 1;
        if ( isnan(params.CR->carr2p[i]) || isinf(params.CR->carr2p[i]) ) params.CR->carr2p[i] = 1;
        if ( isnan(params.CR->carr3p[i]) || isinf(params.CR->carr3p[i]) ) params.CR->carr3p[i] = 1;
        if ( isnan(params.CR->carr4p[i]) || isinf(params.CR->carr4p[i]) ) params.CR->carr4p[i] = 1;
        if ( isnan(params.CR->carr5p[i]) || isinf(params.CR->carr5p[i]) ) params.CR->carr5p[i] = 1;
        if ( isnan(params.CR->carr6p[i]) || isinf(params.CR->carr6p[i]) ) params.CR->carr6p[i] = 1;
        if ( isnan(params.CR->carr7p[i]) || isinf(params.CR->carr7p[i]) ) params.CR->carr7p[i] = 1;
        if ( isnan(params.CR->carr8p[i]) || isinf(params.CR->carr8p[i]) ) params.CR->carr8p[i] = 1;
        if ( isnan(params.CR->carr9p[i]) || isinf(params.CR->carr9p[i]) ) params.CR->carr9p[i] = 1;

        if ( isnan(params.CR->carr1f[i]) || isinf(params.CR->carr1f[i]) ) params.CR->carr1f[i] = 1;
        if ( isnan(params.CR->carr2f[i]) || isinf(params.CR->carr2f[i]) ) params.CR->carr2f[i] = 1;
        if ( isnan(params.CR->carr3f[i]) || isinf(params.CR->carr3f[i]) ) params.CR->carr3f[i] = 1;
        if ( isnan(params.CR->carr4f[i]) || isinf(params.CR->carr4f[i]) ) params.CR->carr4f[i] = 1;
        if ( isnan(params.CR->carr5f[i]) || isinf(params.CR->carr5f[i]) ) params.CR->carr5f[i] = 1;
        if ( isnan(params.CR->carr6f[i]) || isinf(params.CR->carr6f[i]) ) params.CR->carr6f[i] = 1;
        if ( isnan(params.CR->carr7f[i]) || isinf(params.CR->carr7f[i]) ) params.CR->carr7f[i] = 1;
        if ( isnan(params.CR->carr8f[i]) || isinf(params.CR->carr8f[i]) ) params.CR->carr8f[i] = 1;
        if ( isnan(params.CR->carr9f[i]) || isinf(params.CR->carr9f[i]) ) params.CR->carr9f[i] = 1;
    }
    */

    zz=(1-sqr(params.sigma))/sqr(params.sigma);     /*Z=20464*/

    //for ii:=0 to 100 do MemoTesti.Lines[ii]:=FloatToStr(carr4[ii]);   /*Z=20471*/
    //for ii:=0 to 100 do MemoTestj.Lines[ii]:=FloatToStr(carr5[ii]);
    //for ii:=0 to 100 do MemoTestk.Lines[ii]:=FloatToStr(carr6[ii]);

    //(***************************************)
    //(*** calculate unit cell orientation ***)
    //(***************************************)
    double ucvol;
    if ( ltype != 12 )
    {   /*Z=20483 aber auch in Z=20185 ... */
        ButtonHKLClick( ltype, latpar1ptr, latpar2ptr );    // ltype=0..11 und nicht mehr ....
        peakmax1 = latpar1(1,0);
        peakmax2 = latpar2(1,0);    /*Z0311=19753*/
//#ifndef __CUDACC__
//        qDebug() << "nach ButtonHKLClick(2)" << peakmax1 << peakmax2;
//#endif

        /*Z=20487*/
        corotations(params.uca,params.ucb,params.ucc,params.alpha_deg,params.beta_deg,params.gamma_deg,          // in
                    ucn1,ucn2,ucn3,ucpsi, lat1d,lat2d,lat3d,                 // in
                    ri11,ri12,ri13,ri21,ri22,ri23,ri31,ri32,ri33,            // out
                    rt11,rt12,rt13,rt21,rt22,rt23,rt31,rt32,rt33,ucvol,      // out
                    nuvwx,nuvwy,nuvwz,uuvwx,uuvwy,uuvwz,vuvwx,vuvwy,vuvwz,   // out
                    nhklx,nhkly,nhklz,uhklx,uhkly,uhklz,vhklx,vhkly,vhklz);  // out

        //        rimp^[1,1]:=ri11;    rimp^[1,2]:=ri12;    rimp^[1,3]:=ri13;  /*Z=20534*/
        //        rimp^[2,1]:=ri21;    rimp^[2,2]:=ri22;    rimp^[2,3]:=ri23;
        //        rimp^[3,1]:=ri31;    rimp^[1,2]:=ri32;    rimp^[1,3]:=ri33;
        //        rtmp^[1,1]:=rt11;    rtmp^[1,2]:=rt12;    rtmp^[1,3]:=rt13;
        //        rtmp^[2,1]:=rt21;    rtmp^[2,2]:=rt22;    rtmp^[2,3]:=rt23;
        //        rtmp^[3,1]:=rt31;    rtmp^[1,2]:=rt32;    rtmp^[1,3]:=rt33;
    }
    peaknorm2=0;

    //(***************************************************)
    //(*** isotropic peak list for allowed reflections ***)
    //(***************************************************)

    if ( ltype != 12 )
    {   /*Z=20545*/
        double sphno;
        for ( int peakct1=1; peakct1<=peakmax1; peakct1++ )     /*Z=20549*/
        {
            h=latpar1(peakct1,1);
            k=latpar1(peakct1,2);
            l=latpar1(peakct1,3);
            mhkl=latpar1(peakct1,4);    //(* check *)
            fhkl_c(ltype,h,k,l,params.uca,params.ucb,params.ucc,params.alpha_deg,params.beta_deg,params.gamma_deg,
                   sphno,fhkl,qhkl,qhkl0);

            //if ( peakct1 < 5 ) qDebug() << "PEAK1" << peakct1 << sphno << fhkl << qhkl << qhkl0;

            setLatpar1(peakct1,5, round(fhkl) ); // wird danach aber nie wieder gebraucht ...
            latpar[1]=ucvol;         /*Z=20557*/
            latpar[2]=sphno;         // Zwischenspeicher für 'die letzten' berechneten Werte ...
            latpar[3]=dwfactor;
            // [1] und [2] wurden bei NewCPU anders ausgelesen. Macht zwar nichts aus, da diese beiden Werte
            // beim Ergebnis miteinander multipliziert und sonst nicht verwendet werden.
        }
        if ( peakmax1 == 0 )
        {   // Damit auf jeden Fall gültige Werte in das latpar[] Array geschrieben werden...
            // Sonst wird jeder Pixel zu nan weil durch cubevol geteilt wird ...
            latpar[1]=1;
            latpar[2]=ucvol;
        }

        //(* anisotropic peak list *)
        //int hkli=0;
        for ( int peakct2=1; peakct2<=peakmax2; peakct2++ )     /*Z=20572*/
        {
            h=latpar2(peakct2,1);
            k=latpar2(peakct2,2);
            l=latpar2(peakct2,3);
            //int mhkl=1;
            fhkl_c(ltype,h,k,l,params.uca,params.ucb,params.ucc,params.alpha_deg,params.beta_deg,params.gamma_deg,
                   sphno,fhkl,qhkl,qhkl0);

            //if ( peakct2 < 5 ) qDebug() << "PEAK2" << peakct2 << sphno << fhkl << qhkl << qhkl0;

            setLatpar2(peakct2,5, round(fhkl) );
            setLatpar2(peakct2,4, 1 );      /*Z=20583  =mhkl welches vorher auf 1 gesetzt wurde*/

            setLatpar3(peakct2,1, 2*M_PI/params.uca );
            qxhkl=(2*M_PI)*(ri11*h+ri21*k+ri31*l);
            qyhkl=(2*M_PI)*(ri12*h+ri22*k+ri32*l);
            qzhkl=(2*M_PI)*(ri13*h+ri23*k+ri33*l);
            qxyhkl=(2*M_PI)*sqrt(h*h*(ri11*ri11+ri12*ri12+ri13*ri13)+k*k*(ri21*ri21+ri22*ri22+ri23*ri23));
            qhkl=sqrt(qxhkl*qxhkl+qyhkl*qyhkl+qzhkl*qzhkl);
            setLatpar3(peakct2,2, qxhkl );
            setLatpar3(peakct2,3, qyhkl );
            setLatpar3(peakct2,4, qzhkl );
            setLatpar3(peakct2,5, qhkl );
            setLatpar3(peakct2,6, qhkl/(2*M_PI) );
            qxhklt=(2*M_PI)*(rt11*h+rt21*k+rt31*l);
            qyhklt=(2*M_PI)*(rt12*h+rt22*k+rt32*l);
            qzhklt=(2*M_PI)*(rt13*h+rt23*k+rt33*l);
            setLatpar3(peakct2,7, qxhklt );
            setLatpar3(peakct2,8, qyhklt );
            setLatpar3(peakct2,9, qzhklt );
            g3=4*M_PI*M_PI/(2*M_PI*qhkl*qhkl);
            x2phihkl=4*qhkl*qhkl/(M_PI*phiwidth*phiwidth);
            switch ( shp )      /*Z=20614*/
            {
            case cbpeakLorentzian: // 1
                peaknorm1=lorentznorm3(x2phihkl);
                break;
            case cbpeakGaussian: // 2
                peaknorm1=gaussnorm3(x2phihkl);
                break;
            case cbpeakMod1Lorentzian: // 3
                peaknorm1=pearsonnorm3(x2phihkl,2);
                break;
            case cbpeakMod2Lorentzian: // 4
                peaknorm1=pearsonnorm3(x2phihkl,3/2.);
                break;
            case cbpeakPseudoVoigt: // 5
                peaknorm1=lorentznorm3(x2phihkl);
                // UNKLAR peaknorm2=gaussnorm3(x2psihkl);
                peaknorm2=gaussnorm3(x2phihkl);
                // x2psihkl ist auch im Quellcode nicht gesetzt
                // wird erst bei späteren Berechnungen gesetzt
                break;
            case cbpeakPearsonVII: // 6
                peaknorm1=pearsonnorm3(x2phihkl,beta);
                break;
            case cbpeakGamma: // 7
                peaknorm1=gaussnorm3(x2phihkl);
                break;
            case cbpeakAnisotropicGaussian: // 8
                peaknorm1 = 0;  // wird lt. Prof. Förster hier nicht benötigt
                peaknorm2 = 0;
                D8( if ( peakct2 == 1 ) qDebug() << "prepCalc cbpeakAnisotropicGaussian" );
                break;
            }
            setLatpar3(peakct2,10, g3 );
            setLatpar3(peakct2,11, peaknorm1 );
            setLatpar3(peakct2,12, peaknorm2 );
            setLatpar3(peakct2,13, qxyhkl );    /*Z=20631*/

            //(* Peak reference table entries *)
            /*Z=20633 Der folgende Code wird nur für Debugausgaben verwendet und sonst nicht. Daher hier weglassen. */

        } // for ( int peakct2=1; peakct2<=peakmax2; peakct2++ )
    } /*Z=20694  if ltype != 12  <=>  if generic and lattice */




    /*Z=24607 Der folgende Code ist direkt vor den ihex/i Schleifen */
    sinphic=sin(-phi*M_PI/180.);
    cosphic=cos(-phi*M_PI/180.);

    ucl1 = sintheta*cosphi;
    ucl2 = sintheta*sinphi;
    ucl3 = -costheta;

    //bfactor:=StrToFloat(EditRotAlpha.Text);

    if ( shp==8 )
    {
        D8( qDebug() << "prepCalc start qrombdeltac, theta phi:" << theta << phi );
        qrombdeltac(params.length, params.radius, /*p1, sigma, dbeta,*/ theta, phi, 1,1,1,
                    9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,
                    /*i0=*/2, /*i1=*/0, /*i2=*/0, /*i3=*/0, /*i4=*/0, params.CR->carr1p, /*Erg=*/norm );
        D8( qDebug() << "prepCalc qrombdeltac fertig" );
    }

    psphere_init();

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
        TPVINIT( "phi",             phi,               tpvRandomHelper::phi  );
        // TODO: noch nicht in der GUI
        TPVINIT( "Base",            base,              tpvRandomHelper::norm );
        TPVINIT( "EditRho",         params.rho,        tpvRandomHelper::norm );
        TPVINIT( "EditDbeta",       params.dbeta,      tpvRandomHelper::norm );
    }

#define zuf() (((static_cast<double>(rand())/RAND_MAX * 200.0) - 100.0) / 100.0) // -1 .. 1

    std::string retval="";
    for ( std::string id : ids )
    {
        int pos = id.find("=");
        std::string key = id.substr(0,pos);
        std::string val = id.substr(pos+1);
        double dval = std::stod(val);
        if ( dval < 0.00001 && val.length() > 2 )
        {   // Typische Werte sind hier 0.1 bis 0.9 und keine negativen Werte!
            int pk = val.find(".");
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
void SasCalc_GENERIC_calculation::doCalculation( int numThreads, progressAndAbort pa )
{
    if ( pa != nullptr ) progAndAbort = pa;
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
        if ( progAndAbort )
            progAndAbort( -1 );     // Signal the usage of GPU and disables the progress bar

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
                if ( progAndAbort )
                {   // --> wird nur noch bei Tests genutzt, im Normalfall nicht mehr!
                    // Displays the current progess level in percent
                    // and checks if the cancel button was clicked
                    if ( progAndAbort( (100*(ihex-zzmin))/(zzmax-zzmin) ) ) break;
                }
                if ( _endThread ) break;
                doIntCalc_GENERIC( ihex++ );
            } // for ihex
        } // numThreads <= 1

        else
        {   // Special work for Multi-Threads...
            threads     = new pthread_t[numThreads];
            thread_args = new int[numThreads];
            memset( threads, 0, sizeof(pthread_t)*numThreads ); // for ( int t=0; t<numThreads; t++ ) threads[t] = 0;
            int ihex=zzmin;
            while ( ihex<zzmax )  // (***   z-loop  ***)
            {
                if ( progAndAbort )
                {   // --> wird nur noch bei Tests genutzt, im Normalfall nicht mehr!
                    // Displays the current progess level in percent
                    // and checks if the cancel button was clicked
                    if ( progAndAbort( (100*(ihex-zzmin))/(zzmax-zzmin) ) ) break; // Wert in %
                }
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
                    if ( threads[t] )
                    {
                        pthread_join(threads[t], nullptr);
                        threads[t] = 0;
                    }
                }
            } // for ihex
            for ( int t=0; t<numThreads; t++ )
            {
                if ( threads[t] )
                {
                    pthread_join(threads[t], nullptr);
                    threads[t] = 0;
                }
            }
            delete threads;
            delete thread_args;
            threads = nullptr;
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
    //if ( pa != nullptr ) progAndAbort = pa;
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

    //std::cerr << "FitParams: [X] pix=" << params.latFit.pixx << ", cols=" << params.latFit.cols << ", center=" << params.latFit.centerx << std::endl;
    //std::cerr << "           [Y] pix=" << params.latFit.pixy << ", rows=" << params.latFit.rows << ", center=" << params.latFit.centery << std::endl;

    // Use a high resolution clock to get the calculation time of the GPUs
    auto start2 = std::chrono::high_resolution_clock::now();

#ifdef __CUDACC__
    if ( gpuAvailable() && numThreads == 0 )    // GPU abschaltbar
    {
        _endThread = false;
        if ( progAndAbort )
            progAndAbort( -1 );     // Signal the usage of GPU and disables the progress bar

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
                if ( progAndAbort )
                {   // --> wird nur noch bei Tests genutzt, im Normalfall nicht mehr!
                    // Displays the current progess level in percent
                    // and checks if the cancel button was clicked
                    if ( progAndAbort( (100*x)/_fitWidth ) ) break;
                }
                if ( _endThread ) break;
                doIntFitCalc_GENERIC( x++ );
            } // for x
        } // numThreads <= 1

        else
        {   // Special work for Multi-Threads...
            threads     = new pthread_t[numThreads];
            thread_args = new int[numThreads];
            memset( threads, 0, sizeof(pthread_t)*numThreads ); // for ( int t=0; t<numThreads; t++ ) threads[t] = 0;
            int x=0;
            while ( x<_fitWidth /*&& x<10/ *TEST*/ )
            {
                if ( progAndAbort )
                {   // --> wird nur noch bei Tests genutzt, im Normalfall nicht mehr!
                    // Displays the current progess level in percent
                    // and checks if the cancel button was clicked
                    if ( progAndAbort( (100*x)/_fitWidth ) ) break; // Wert in %
                }
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
                    if ( threads[t] )
                    {
                        pthread_join(threads[t], nullptr);
                        threads[t] = 0;
                    }
                }
            } // for ihex
            for ( int t=0; t<numThreads; t++ )
            {
                if ( threads[t] )
                {
                    pthread_join(threads[t], nullptr);
                    threads[t] = 0;
                }
            }
            delete threads;
            delete thread_args;
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
 * @param FCC  - reference to class with parameters and subfunctions
 * @param ihex - vertical pixel index
 * @param i    - horizontal pixel index
 * @param lamv - running vertical axis value
 * Calculation from Pascalprogram (20210818-crystal3d1.pas + updates) - only part Generic
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline void SasCalc_GENERIC_calculation::doIntCalc_GENERIC_F( const SasCalc_GENERIC_calculation& FCC,
                                                              int ihex, int i )
{
#ifdef COPY_FITDATA_TO_GPU  // Eckpixeltest
    if ( FCC.arrDataForFitUsed )
    {   // Spezielle Aktionen (Maskieren und FQS) für das Simplex-2D-Fit
        if ( FCC.arrDataForFit[1+FCC.IDX(ihex,i)] <= FCC.arrDataForFit[1] )
        {   // Datenwert <= Eckpixel --> Ergebnis immer 0.
            FCC.setXYIntensity( ihex, i, 0.0 );
            return;
        }
    }
#endif

    //lamu   = (FCC.qmax * i) / FCC.iimax;    // (*** running index for x-axis ***)
    //{NV} alles mit "RadioButtonSlice.Checked=true" lasse ich weg, da ich hier noch keine
    //     Schnitte behandele.

    // Einrechnen des Beamstops (d.h. Verschiebung des Zentrums)
    //int mdet = ihex + BCT.zmax + 1;     // (* mth pixel *)
    //int ndet = i + BCT.zmax + 1;        // (* nth pixel *)
    /*Z=24635*/
    // Im Pascal-Programm ist der Beamstop für Null in der Ecke.
    // Hier ist der Beamstop für Null in der Mitte gerechnet.
    double xdet = FCC.pixx * (i   /*+FCC.zmax+1*/ - FCC.beamX0);
    double ydet = FCC.pixy * (ihex/*+FCC.zmax+1*/ - FCC.beamY0);
    double rdet = sqrt(xdet*xdet+ydet*ydet);
    double phidet = atan2(ydet,xdet);
    double thetadet = atan2(rdet,FCC.det);
    double qx = 2*M_PI*cos(phidet)*sin(thetadet)/FCC.wave;
    double qy = 2*M_PI*sin(phidet)*sin(thetadet)/FCC.wave;
    double qz = 2*M_PI*(1-cos(thetadet))/FCC.wave;

    // Call q(x,y,z) Routine (also from 2dFit)
    double pixval = doIntCalc_GENERIC_q_xyz( FCC, qx, qy, qz );
    if ( FCC._endThread ) return;

    FCC.setXYIntensity( ihex, i, pixval );

#ifndef __CUDACC__
    FCC.setDebugIntensity( false );   // Damit nur beim ersten Durchlauf eine Kontrollausgabe kommt.
#endif

    DBGFILE( if ( fdbg ) fdbg->write(qPrintable(QString("%1; %2;   %3; %4; %5;   %6;\n")
                                       .arg(ihex).arg(i)
                                       .arg(qx).arg(qy).arg(qz)
                                       .arg(pixval)
                                   )); )


} /* doIntCalc_GENERIC_F() */


#ifdef __CUDACC__
__host__ __device__
#endif
inline void SasCalc_GENERIC_calculation::doIntFitCalc_GENERIC_F( const SasCalc_GENERIC_calculation& FCC,
                                                        int x, int y )
{
    for ( int i=0; i<4; i++ )
    {
        if ( FCC.noFitX0[i] < 0 ) continue;
        if ( x < FCC.noFitX0[i] ) continue;
        if ( FCC.noFitX1[i] < x ) continue;
        if ( y < FCC.noFitY0[i] ) continue;
        if ( FCC.noFitY1[i] < y ) continue;
        size_t idx = x + (FCC._fitWidth * y);
        FCC.arrFitFqs[idx] = 0.0;
        DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; FitRect %4\n")
                                           .arg(x).arg(y).arg(FCC.arrFitData[idx]).arg(i)
                                       )); )
        return;
    }

    int xx = x - FCC.params.latFit.centerx;
    int yy = y - FCC.params.latFit.centery;
    if ( FCC.fitBorderPixel > 0 || FCC.fitBStopPixel > 0 )
    {   // Ausblendungen über Pixelangaben an Rand und Mitte
        if ( x < FCC._xmin + FCC.fitBorderPixel || x >= FCC._xmax - FCC.fitBorderPixel ||
             y < FCC._ymin + FCC.fitBorderPixel || y >= FCC._ymax - FCC.fitBorderPixel )
        {
            size_t idx = x + (FCC._fitWidth * y);
            FCC.arrFitFqs[idx] = 0.0;
            DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                    fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; BorderPixel\n")
                                               .arg(x).arg(y).arg(FCC.arrFitData[idx])
                                           )); )
            return;
        }
        if ( x >= FCC.params.latFit.centerx - FCC.fitBStopPixel &&
             x <  FCC.params.latFit.centerx + FCC.fitBStopPixel &&
             y >= FCC.params.latFit.centery - FCC.fitBStopPixel &&
             y <  FCC.params.latFit.centery + FCC.fitBStopPixel )
        {
            size_t idx = x + (FCC._fitWidth * y);
            FCC.arrFitFqs[idx] = 0.0;
            DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                    fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; BSpixel\n")
                                               .arg(x).arg(y).arg(FCC.arrFitData[idx])
                                           )); )
            return;
        }
    }
    else if ( FCC.fitBStopPixel == -1 )
    {   // Ausblendungen per Eck-Pixel-Wert
        size_t idx = x + (FCC._fitWidth * y);
        if ( FCC.arrFitData[idx] <= FCC.arrFitData[0] )
        {
            FCC.arrFitFqs[idx] = 0.0;
            DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                    fdbg->write(qPrintable(QString("%1; %2;   -/-; -/-;   -/-; -/-; -/-;   -/-; %3; 0; Ecke %4 <= %5\n")
                                               .arg(x).arg(y).arg(FCC.arrFitData[idx]).arg(FCC.arrFitData[idx]).arg(FCC.arrFitData[0])
                                           )); )
            return;
        }
    }
    //else
    //{   // Keine Ausblenungen
    //}

    // Einrechnen des Beamstops (d.h. Verschiebung des Zentrums)
    double xdet = FCC.params.latFit.pixx * xx;
    double ydet = FCC.params.latFit.pixy * yy;
    double rdet = sqrt(xdet*xdet+ydet*ydet);
    double phidet = atan2(ydet,xdet);
    double thetadet = atan2(rdet,FCC.params.latFit.distance);
    double qx = 2*M_PI*cos(phidet)*sin(thetadet)/FCC.params.latFit.wavelen;
    double qy = 2*M_PI*sin(phidet)*sin(thetadet)/FCC.params.latFit.wavelen;
    double qz = 2*M_PI*(1-cos(thetadet))/FCC.params.latFit.wavelen;

    //DBGFILE( if ( !( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 ) ) return; )

    // Call q(x,y,z) Routine (also from 2dFit)
    double pixval = doIntCalc_GENERIC_q_xyz( FCC, qx, qy, qz );
    if ( FCC._endThread ) return;

    size_t idx = x + (FCC._fitWidth * y);
    if ( pixval > 0 && FCC.arrFitData[idx] > 0 )
        FCC.arrFitFqs[idx] = FQSVERGL( pixval, FCC.arrFitData[idx] );
    else
        FCC.arrFitFqs[idx] = 0.0;

    DBGFILE( if ( fdbg ) //if ( (x<3 && y<3) || fabs(qy-(-1.28525))<0.001 )
                fdbg->write(qPrintable(QString("%1; %2;   %3; %4;   %5; %6; %7;   %8; %9; %10\n")
                                           .arg(x).arg(y)
                                           .arg(xdet).arg(ydet)
                                           .arg(qx).arg(qy).arg(qz)
                                           .arg(pixval).arg(FCC.arrFitData[idx]).arg(FCC.arrFitFqs[idx])
                                       )); )

//    if ( (x < 4 /*&& y < 5*/) || (x >= FCC._fitWidth-5 && y >= FCC._fitHeight) )
//        qDebug() << "Fit-v:" << x << y << idx
//                 << "val" << pixval << FCC.arrFitData[idx] << FCC.arrFitFqs[idx]
//                 << "q" << qx << qy << qz;

} /* doIntFitCalc_GENERIC_F() */



/**
 * @brief SasCalc_GENERIC_calculation::doIntCalc_GENERIC_q_xyz
 * @param FCC  - reference to class with parameters and subfunctions
 * @param qx   - Coordinates in the q dimension
 * @param qy   - "
 * @param qz   - "
 * @return     - calculated value for this coordinate
 * Calculation from Pascalprogram (20210818-crystal3d1.pas + updates) - only part Generic
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::doIntCalc_GENERIC_q_xyz(const SasCalc_GENERIC_calculation& FCC,
                                                         double qx, double qy, double qz)
{
    double /*lamu,*/ /*qabs,*/ pq, fq, intensity, radintensity;
    double /*shkl,*/ /*fhkl,*/ x2;
    double dqs1, dqs2, dqs3, sq;
    //Double3 q3hkl; // qxhkl, qyhkl, qzhkl
    //Double3 dq3;   // dqx, dqy, dqz

    double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;   /*Z=24650*/
    //s = q/(2*M_PI);

    /*Z=24654 if gisaxs ... erstmal weglassen und das auch bei allen weiteren Abfragen!!!  */

    double delta = 2.0 * FCC.params.radius;
    double pqiso=1;
    double limql=0;  // Achtung: L statt 1 !!!
    double limqlf=0; // Wird nur als Parameter bei formfq() verwendet, sonst nicht
    double szqiso, szq;

#define lattice (FCC.ltype != 12)  // 12=None,  Im Pascalprogramm als globale boolsche Variable in /*Z=20133*/

    fq = 1; //TODO
    // fq wird bei der Berechnung von szq und somit für die Intensität verwendet und nur bei 'lattice' gesetzt

    if ( FCC.ComboBoxParticle==cbpartSphere /*0*/ )
    {
        /*Z=24671*/
        pq=FCC.formpq( FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                       FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, limql, FCC.limq1, FCC.limq2, FCC.limq3,
                       FCC.limq4, FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx, qy, q, FCC.norm,
                       FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p, FCC.params.CR->carr3p,
                       FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p, FCC.params.CR->carr9p /*, FCC.carr11pm, FCC.carr22pm*/ );
        //fq:=pq;
        if ( lattice )  /*Z=24679*/
            fq=FCC.formfq( FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                           FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, limqlf, FCC.limq1f, FCC.limq2f, FCC.limq3f,
                           FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx, qy, q, FCC.norm,
                           FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1f, FCC.params.CR->carr2f, FCC.params.CR->carr3f,
                           FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f /*, FCC.carr11pm, FCC.carr22pm*/ );
        pqiso = pq;    /*Z=24693*/
    } // if sphere

    if ( FCC.ComboBoxParticle==cbpartCylinder /*1*/ )     /*Z=24701*/
    {   /* cylinders */   /*Z0311=24157*/

        /* isotropic cases */   /*Z=24703*/
        if ( FCC.ordis==ordis_Isotropic /*7*/ )
        {   /*Z=24705*/
            pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                            FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1, FCC.limq2, FCC.limq3,
                            FCC.limq4, FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx, qy, q,
                            FCC.norm,
                            FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                            FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                            FCC.params.CR->carr9p /*,carr11pm^,carr22pm^*/ );
            /* fq:=pq; */
            if ( lattice )   /*Z=24708*/
                fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx, qy, q, FCC.norm,
                                FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1f, FCC.params.CR->carr2f,
                                FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                /*,carr11pm^,carr22pm^*/ );
        }   /*Z=24711*/

        /* perfect orientation */   /*Z=24714*/
        else if ( FCC.ordis==ordis_ZDir /*6*/ )
        {   /*Z=24716*/
            switch ( FCC.orcase )
            {
            case 1:     // General: phi!=0 && phi!=90 && theta!=0 && theta!=90
                /*Z=24717*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, sqrt(FCC.cosphi*qx*qx+FCC.sinphi*qy*qy+eps9),
                                FCC.limq1, FCC.limq2, FCC.limq3, FCC.limq4, FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8,
                                FCC.limq9, qx, qy, qx*FCC.cosphi*FCC.sintheta, qy*FCC.sinphi*FCC.sintheta, q, FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p, FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p,
                                FCC.params.CR->carr8p, FCC.params.CR->carr9p /*,carr11pm^,carr22pm^*/ );   /*Z0311=24172*/
                /* fq:=pq; */   /*Z=24719*/
                if ( lattice )
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, sqrt(FCC.cosphi*qx*qx+FCC.sinphi*qy*qy+eps9),
                                    FCC.limq1f, FCC.limq2f, FCC.limq3f, FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy,
                                    qx*FCC.cosphi*FCC.sintheta, qy*FCC.sinphi*FCC.sintheta, q, FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2f,  // TODO: wirklich carr1p und nicht carr1f?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f/*,
                                    carr11pm^,carr22pm^*/ );   /*Z0311=24176*/
                break;   /*Z=24723*/
            case 2:     // X-Axis phi==0 && theta==90
                /*Z=24725*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, fabs(qx), FCC.limq1, FCC.limq2, FCC.limq3,
                                FCC.limq4, FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx, 0, q,FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                                FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                                FCC.params.CR->carr9p /*, carr11pm^,carr22pm^*/ );   /*Z0311=24180*/
                /* fq:=pq; */   /*Z=24727*/
                if ( lattice )
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, fabs(qx), FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                    FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx, 0, q,FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2f,  // TODO: carr1p?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                    /*,carr11pm^,carr22pm^*/ );   /*Z0311=24184*/
                break;   /*Z=24731*/
            case 3:     // Y-Axis phi==90 && theta==90
                /*Z=24733*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, fabs(qy), FCC.limq1, FCC.limq2, FCC.limq3,
                                FCC.limq4, FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, 0, qy, q, FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                                FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                                FCC.params.CR->carr9p /*, carr11pm^,carr22pm^*/ );   /*Z0311=24188*/
                /* fq:=pq; */   /*Z=24735*/
                if ( lattice )
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, fabs(qy), FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                    FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, 0, qy, q, FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2f,  // TODO: carr1p?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                    /*,carr11pm^,carr22pm^*/ );   /*Z0311=24192*/
                break;   /*Z=24739*/
            case 4:     // Z-Axis (phi==0 || phi==90) && theta==0
                /*Z=24741*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1, FCC.limq2, FCC.limq3, FCC.limq4,
                                FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx, qy, q, FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                                FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                                FCC.params.CR->carr9p /*,carr11pm^,carr22pm^*/ );   /*Z0311=24196*/
                /* fq:=pq; */   /*Z=24743*/
                if ( lattice )
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                    FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx, qy, q, FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2f,  // TODO carr1p?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                    /*,carr11pm^,carr22pm^*/ );   /*Z0311=24200*/
                break;   /*Z=24747*/
            } // switch orcase   /*Z=24748*/
        } // if ( FCC.ordis==6 )

        /* general orientation */   /*Z0311=24204*/
        else if ( FCC.ordis==ordis_Gaussian /*0*/ )
        {   /*Z=24754*/
            switch ( FCC.orcase )
            {
            case 1:
                /* general orientation */   /*Z=24755*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1, FCC.limq2, FCC.limq3, FCC.limq4,
                                FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx*FCC.cosphic-qy*FCC.sinphic,
                                qx*FCC.sinphic+qy*FCC.cosphic, q, FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                                FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                                FCC.params.CR->carr9p /*,carr11pm^,carr22pm^*/ );   /*Z0311=24208*/
                //fq=pq;    /*Z=24757, neu weg*/
                if ( lattice ) /* fq:=pq; */   /*Z=24758*/
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                    FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx*FCC.cosphic-qy*FCC.sinphic,
                                    qx*FCC.sinphic+qy*FCC.cosphic, q, FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2f,  // TODO: carr1p?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                    /*,carr11pm^,carr22pm^*/ );   /*Z0311=24212*/
                break;   /*Z=24761*/
            case 2:
                /* x-axis */   /*Z=24763*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1, FCC.limq2, FCC.limq3, FCC.limq4,
                                FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx, qy, q, FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                                FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                                FCC.params.CR->carr9p /*,carr11pm^,carr22pm^ */);   /*Z0311=24216*/
                //fq=pq;    /*Z=24765, neu weg*/
                //ffq = pq;   /*Z=24766, neu, wird aber nur hier verwendet*/
                if ( lattice ) /* fq:=pq; */   /*Z0311=24219*/
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                    FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx, qy, q, FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1f, FCC.params.CR->carr2f,  // TODO: carr1p?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                    /*,carr11pm^,carr22pm^*/ );   /*Z0311=24221*/
                szq = /*ffq*/ pq;   /*Z=24770 neu*/
                break;   /*Z=24771*/
            case 3:
                /* y-axis */   /*Z0311=24224*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1, FCC.limq2, FCC.limq3, FCC.limq4,
                                FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx, qy, q, FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                                FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                                FCC.params.CR->carr9p /*,carr11pm^,carr22pm^*/ );   /*Z0311=24226*/
                //fq=pq;    /*Z=24775 neu weg*/
                if ( lattice ) /* fq:=pq; */   /*Z=24776*/
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                    FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx, qy, q, FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1f, FCC.params.CR->carr2f,  // TODO: carr1p?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                    /*,carr11pm^,carr22pm^*/ );   /*Z0311=24230*/
                break;   /*Z=24779*/
            case 4:
                /* z-axis */   /*Z=24780*/
                pq = FCC.formpq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1, FCC.limq2, FCC.limq3,
                                FCC.limq4, FCC.limq5, FCC.limq6, FCC.limq7, FCC.limq8, FCC.limq9, qx, qy, qx, qy, q, FCC.norm,
                                FCC.por, FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1p, FCC.params.CR->carr2p,
                                FCC.params.CR->carr3p, FCC.params.CR->carr4p, FCC.params.CR->carr5p, FCC.params.CR->carr6p, FCC.params.CR->carr7p, FCC.params.CR->carr8p,
                                FCC.params.CR->carr9p /*,carr11pm^,carr22pm^*/ );   /*Z0311=24234*/
                fq=pq;    /*Z0311=24235*/
                if ( lattice ) /* fq:=pq; */   /*Z=24784*/
                    fq = FCC.formfq(FCC.params.length, FCC.params.radius, FCC.params.sigmal, FCC.params.sigma, FCC.params.p1,
                                    FCC.params.rho, FCC.alphash, FCC.theta, FCC.phi, q, FCC.limq1f, FCC.limq2f, FCC.limq3f,
                                    FCC.limq4f, FCC.limq5f, FCC.limq6f, qx, qy, qx, qy, q, FCC.norm,
                                    FCC.part, FCC.cs, FCC.ordis, FCC.orcase, FCC.params.CR->myarray, FCC.params.CR->carr1f, FCC.params.CR->carr2f,  // TODO: carr1p?
                                    FCC.params.CR->carr3f, FCC.params.CR->carr4f, FCC.params.CR->carr5f, FCC.params.CR->carr6f, FCC.params.CR->carr7f
                                    /*,carr11pm^,carr22pm^*/ );   /*Z0311=24238*/
                break;   /*Z=24787*/
            } // switch orcase   /*Z=24788*/

        } // if ( FCC.ordis==0 )

    }  /* of part=1 (Cylinder) */   /*Z=24790*/

    if ( FCC.ComboBoxParticle==cbpartDisk /*2*/ )     /*Z0311=24252*/  /*Z=24796*/
    {
        /* isotropic cases */
        if ( FCC.ordis == 7 )
        {
            pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                          FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,q,FCC.limq1,FCC.limq2,FCC.limq3,
                          FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,qx,qy,qx,qy,q,FCC.norm,
                          FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,
                          FCC.params.CR->carr1p,FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,
                          FCC.params.CR->carr5p,FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,
                          FCC.params.CR->carr9p /*,carr11pm^,carr22pm^*/ );
            //fq:=pq;
            if ( lattice )
                fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,q,FCC.limq1f,FCC.limq2f,
                                FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,qx,qy,qx,qy,q,FCC.norm,
                                FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1f,
                                FCC.params.CR->carr2f,FCC.params.CR->carr3f,FCC.params.CR->carr4f,FCC.params.CR->carr5f,
                                FCC.params.CR->carr6f,FCC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );
        }

        /* perfect orientation */
        if ( FCC.ordis==6 ) /*Z=24808*/
        {
            switch ( FCC.orcase )
            {
            case 1:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                sqrt(FCC.sinphi*qx*qx+FCC.cosphi*qy*qy+eps9),
                                FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,
                                qx,qy,qx*FCC.cosphi*FCC.sintheta/q,qy*FCC.sinphi*FCC.sintheta/q,q,FCC.norm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice )
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    sqrt(FCC.sinphi*qx*qx+FCC.cosphi*qy*qy+eps9),
                                    FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,qx*FCC.cosphi*FCC.sintheta/q,qy*FCC.sinphi*FCC.sintheta/q,q,FCC.norm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1f, /*TODO: hier stand carr1p*/
                                    FCC.params.CR->carr2f,FCC.params.CR->carr3f,FCC.params.CR->carr4f,FCC.params.CR->carr5f,
                                    FCC.params.CR->carr6f,FCC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );
                break;
            case 2:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                fabs(qy),
                                FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,
                                qx,qy,qx/q,0,q,FCC.norm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice )
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    fabs(qy),
                                    FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,qx/q,0,q,FCC.norm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1f, /*TODO: hier stand carr1p*/
                                    FCC.params.CR->carr2f,FCC.params.CR->carr3f,FCC.params.CR->carr4f,FCC.params.CR->carr5f,
                                    FCC.params.CR->carr6f,FCC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );
                break;
            case 3:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                fabs(qx),FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,
                                FCC.limq8,FCC.limq9,
                                qx,qy,0,qy/q,q,FCC.norm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice )
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    fabs(qx),
                                    FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,0,qy/q,q,FCC.norm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1f, /*TODO: hier stand carr1p*/
                                    FCC.params.CR->carr2f,FCC.params.CR->carr3f,FCC.params.CR->carr4f,FCC.params.CR->carr5f,
                                    FCC.params.CR->carr6f,FCC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );
                break;
            case 4:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                q,
                                FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,
                                qx,qy,qx,qy,q,FCC.norm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice )
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    q,
                                    FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,qx,qy,q,FCC.norm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1f, /*TODO ...*/
                                    FCC.params.CR->carr2f,FCC.params.CR->carr3f,FCC.params.CR->carr4f,FCC.params.CR->carr5f,
                                    FCC.params.CR->carr6f,FCC.params.CR->carr7f /*,carr11pm^,carr22pm^*/ );
                break;
            } // switch orcase
        } /* if ( FCC.ordis==6 ) */  /*Z=24841*/

        //(* isotropic fraction *)
        //if (iso>0) then pqiso:=formpq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,
        //       por,part,cs,7,orcase,myarray,carr1p^,carr2p^,carr3p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr11pm^,carr22pm^)
        //else pqiso:=0.0;

        /* general orientation */
        if ( FCC.ordis==0 ) /*Z=24850*/
        {
#define nnnorm FCC.norm  /*Z=20467*/
            switch ( FCC.orcase )
            {
            case 1:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                q,FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,
                                qx,qy,qx*FCC.cosphic/q-qy*FCC.sinphic/q,qx*FCC.sinphic/q+qy*FCC.cosphic/q,q,nnnorm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice )
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    q,FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,qx*FCC.cosphic/q-qy*FCC.sinphic/q,qx*FCC.sinphic/q+qy*FCC.cosphic/q,q,nnnorm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                    FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                    FCC.params.CR->carr6p,FCC.params.CR->carr7p /*,carr11pm^,carr22pm^*/ );
                break;
            case 2:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                q,FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,
                                qx,qy,qx/q,qy/q,q,nnnorm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice ) //fq:=pq;
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    q,FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,qx/q,qy/q,q,nnnorm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                    FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                    FCC.params.CR->carr6p,FCC.params.CR->carr7p /*,carr11pm^,carr22pm^*/ );
                break;
            case 3:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                q,FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,
                                qx,qy,qx/q,qy/q,q,nnnorm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice ) //fq:=pq;
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    q,FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,qx/q,qy/q,q,nnnorm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                    FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                    FCC.params.CR->carr6p,FCC.params.CR->carr7p /*,carr11pm^,carr22pm^*/ );
                break;
            case 4:
                pq=FCC.formpq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                q,FCC.limq1,FCC.limq2,FCC.limq3,FCC.limq4,FCC.limq5,FCC.limq6,FCC.limq7,FCC.limq8,FCC.limq9,
                                qx,qy,qx,qy,q,nnnorm,
                                FCC.por,FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                FCC.params.CR->carr6p,FCC.params.CR->carr7p,FCC.params.CR->carr8p,FCC.params.CR->carr9p
                                /*,carr11pm^,carr22pm^*/ );
                //fq:=pq;
                if ( lattice ) //fq:=pq;
                    fq=FCC.formfq(FCC.params.length,FCC.params.radius,FCC.params.sigmal,FCC.params.sigma,FCC.params.p1,
                                    FCC.params.rho,FCC.params.alphash,FCC.theta,FCC.phi,
                                    q,FCC.limq1f,FCC.limq2f,FCC.limq3f,FCC.limq4f,FCC.limq5f,FCC.limq6f,
                                    qx,qy,qx,qy,q,nnnorm,
                                    FCC.part,FCC.cs,FCC.ordis,FCC.orcase,FCC.params.CR->myarray,FCC.params.CR->carr1p,
                                    FCC.params.CR->carr2p,FCC.params.CR->carr3p,FCC.params.CR->carr4p,FCC.params.CR->carr5p,
                                    FCC.params.CR->carr6p,FCC.params.CR->carr7p /*,carr11pm^,carr22pm^*/ );
                break;
            } // switch orcase
#undef nnnorm
        } // if ( FCC.ordis==0 )

    } // if ( FCC.ComboBoxParticle==cbpartDisk ) /*Z=24884*/



#ifdef nur_sphere   // Als erster Ansatz


    if ( FCC.ComboBoxParticle==cbpartVesicle )      /*Z0311=24343*/
    {   //{NV}-X/16446
        //pq:=polyvesicle(length,radius,sigma,sigmal,q);
        //if CheckBoxf2q.Checked=true then fq:=f2dpolyvesicle(length,radius,sigma,sigmal,q)
        //else fq:=pq;
    }


    if ( FCC.ComboBoxParticle==cbpartLiposome )     /*Z0311=24349*/
    {   //{NV}-X/16452
        //pq=polyliposome(length,radius,sigma,sigmal,shellno,alphash,ceff,reff,a,b,c,domainsize,aziwidth,3,q);
        //fq=pq;
    }

    if ( FCC.ComboBoxParticle==cbpartCube )     /*Z0311=24356*/
    {   //{NV}-X/16458
        /* isotropic cases */
        //   if ( (ordis=7) and homogeneous) then pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,FCC.params.CR->carr1,FCC.params.CR->carr3,FCC.params.CR->carr5,FCC.params.CR->carr6,FCC.params.CR->carr2);
        //fq:=pq;
    }

        cbpartEllipsoide        /*Z0311=24384*/
        cbpartTriaxEllips       /*Z0311=24453*/
        cbpartSuperEllips       /*Z0311=24479*/

        if gisax  /*Z0311=24504*/

#endif // nur_sphere (erster Ansatz)

/*Z=25084
    if ( Checkbox10 )
    {
    if partsphere then pq:=schulz(length,radius,radiusi,1,1,sigma,q,1);
    if partcylinder then pq:=schulztheta(length,radius,radiusi,1,sigma,q,2);
    if partellipsoid then pq:=schulztheta(length,radius,radiusi,1,sigma,q,4);
    if parttriellipsoid then pq:=schulzthetaphi(length,radius,radiusi,sigma,q,5);
    if partcube then pq:=schulzthetaphi(length,radius,radiusi,sigma,q,7);
    }
*/

    double width = 1;
    if ( FCC.RadioButtonDebyeScherrer )    /*Z=25114*/
        width = 4.0/FCC.params.domainsize;
    if ( FCC.RadioButtonPara )
        width = (4.0/FCC.params.domainsize)+FCC.reldis*FCC.reldis*FCC.dist*q*q;
    //if ( ! FCC.tpvRandomOldValues.empty() )   // doIntCalc... - nicht auf der GPU zulässig
    {
        //if ( fabs(qx) < 0.1 && fabs(qy) < 0.1 && fabs(qz) < 0.1 )
        //    std::cerr << "TPV CALC " << FCC.params.width_zuf << std::endl << std::flush;
        width = width * FCC.params.width_zuf;
    }

/*#ifndef __CUDACC__
    if ( ihex == 0 && i == 0 )
    {
        if ( FCC.RadioButtonDebyeScherrer )
            qDebug() << "WIDTH" << width << "Scherrer" << FCC.params.domainsize;
        if ( FCC.RadioButtonPara )
            qDebug() << "WIDTH" << width << "Para" << FCC.params.domainsize << FCC.reldis << FCC.dist << q;
        if ( ! FCC.RadioButtonDebyeScherrer && ! FCC.RadioButtonPara )
            qDebug() << "WIDTH" << width;
    }
#endif*/

    double /*qhkl0,*/ qxhkl, qyhkl, qzhkl, qhkl, qxhklt, qyhklt, qzhklt, g3, peaknorm1; //, peaknorm2;
    double dqx, dqy, dqz, yphi, psiord, phiord, x2phi;

    double widthiso = 1.0 / FCC.params.uca;  /*Z0311=19448*/
    double sphno=0, cubevol=0; // fürs Debuggen besser

    radintensity = 0.0;     // Immer nur für ein Pixel, somit kein Array
    intensity = 0.0;

    /*** lattice hkl-factor calcuation */
    if ( lattice )  /*Z=25123*/
    {
        sphno    = FCC.latpar[1];   /*Z=25127, nur PC-Version*/
        cubevol  = FCC.latpar[2];
        // dwfactor ist global

        /* isotropic peaks */
        for ( int ii=1; ii<=FCC.peakmax1; ii++ )
        {
            //int h    = FCC.latpar1(ii,1);
            //int k    = FCC.latpar1(ii,2);
            //int l    = FCC.latpar1(ii,3);
            int mhkl    = FCC.latpar1(ii,4);
            int fhkl    = FCC.latpar1(ii,5);
            double qhkl = FCC.latpar3(ii,5);
            double x2   = sqr(q-qhkl)/sqr(widthiso);  /*Z=25144*/
            double sq   = exp(-4*x2/M_PI)/(M_PI*widthiso/2.0);
            radintensity += sq*mhkl*fhkl;
        } /* of peak loop */

        for ( int ii=1; ii<=FCC.peakmax2; ii++ )
        {
            if ( FCC._endThread ) return 0;  // Falls Anwender abgebrochen hat
            int h = FCC.latpar2(ii,1);  /*Z=25158*/
            int k = FCC.latpar2(ii,2);
            int l = FCC.latpar2(ii,3);
            int mhkl = FCC.latpar2(ii,4);
            int fhkl = FCC.latpar2(ii,5);

            //qhkl0 = FCC.latpar3(ii,1);    /*Z=25198*/
            qxhkl = FCC.latpar3(ii,2);
            qyhkl = FCC.latpar3(ii,3);
            qzhkl = FCC.latpar3(ii,4);
            qhkl  = FCC.latpar3(ii,5);
            qxhklt = FCC.latpar3(ii,7);
            qyhklt = FCC.latpar3(ii,8);
            qzhklt = FCC.latpar3(ii,9);
            g3     = FCC.latpar3(ii,10);

            switch ( FCC.shp )  /*Z=25231*/
            {
            case cbpeakAnisotropicGaussian /*8*/:
                if ( FCC.ordis==6/*z-dir*/ || fabs(FCC.params.dbeta) < eps9 )   /*Z=25233*/
                {   /* perfect orientation */
                    D8( qDebug() << "CPU: shp=8(A): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                    dqx = qx-qxhkl;
                    dqy = qy-qyhkl;
                    dqz = qz-qzhkl;
                    dqs1 = (dqx*FCC.ax1.x()+dqy*FCC.ax1.y()+dqz*FCC.ax1.z())/(FCC.ax1n_sigx);
                    dqs2 = (dqx*FCC.ax2.x()+dqy*FCC.ax2.y()+dqz*FCC.ax2.z())/(FCC.ax2n_sigy);
                    dqs3 = (dqx*FCC.ax3.x()+dqy*FCC.ax3.y()+dqz*FCC.ax3.z())/(FCC.ax3n_sigz);
                    x2 = dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                    sq = exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*FCC.sig.x()*FCC.sig.y()*FCC.sig.z()/8.0);
                }
                else if ( FCC.ordis==7/*isotropic*/ )
                {   /* isotropic orientation */
                    D8( qDebug() << "CPU: shp=8(B): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                    width = FCC.sig.length() /*sqrt(sigx*sigx+sigy*sigy+sigz*sigz)*/ /3.0;      /*Z0311=24626*/
/*#ifndef __CUDACC__
                    if ( ihex == 0 && i == 0 )
                        qDebug() << "WIDTH" << width << "Ordis=7, cbpeakAnisotropicGaussian";
#endif*/
                    x2 = (q-qhkl)*(q-qhkl)/(width*width);
                    sq = g3*exp(-4*x2/M_PI)/(M_PI*width/2.0);
                }
                else if ( FCC.ordis==13/*fiber pattern*/ )
                {   /* fiber pattern */
                    D8( qDebug() << "CPU: shp=8(C): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                    // Update 20220608
                    // rotaxphi   = FCC.phi
                    // rotaxtheta = FCC.theta
// TODO: Dieser Code (die ersten beiden if) ist in dem Programm vom 30.07.22 nicht mehr dabei ???
#define sign(x) ((x<0)?-1:1)
                    /* rotated around y-axis */
                    if ( fabs(FCC.phi)<eps9 /*(rotaxphi=0)*/ &&
                        fabs(FCC.theta-90)<eps9 /*(rotaxtheta=90)*/ )
                    {
                        double signq=sign(qyhkl);
                        qyhkl=signq*sqrt(qyhkl*qyhkl+qzhkl*qzhkl);
                        qzhkl=1E-20;
                        dqx=qx-qxhkl;
                        dqy=qy-qyhkl;
                        dqz=qz-qzhkl;
                        dqs1=(dqx*FCC.ax1.x()+dqy*FCC.ax1.y()+dqz*FCC.ax1.z())/(FCC.ax1n_sigx);
                        dqs2=(dqx*FCC.ax2.x()+dqy*FCC.ax2.y()+dqz*FCC.ax2.z())/(FCC.ax2n_sigy);
                        dqs3=(dqx*FCC.ax3.x()+dqy*FCC.ax3.y()+dqz*FCC.ax3.z())/(FCC.ax3n_sigz);
                        x2=dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                        sq=exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*FCC.sig.x()*FCC.sig.y()*FCC.sig.z()/8.);
                    }
                    /* rotated around x-axis */
                    else if ( fabs(FCC.phi-90)<eps9 /*(rotaxphi=90)*/ &&
                             fabs(FCC.theta-90)<eps9 /*(rotaxtheta=90)*/ )
                    {
                        double signq=sign(qxhkl);
                        qxhkl=signq*sqrt(qxhkl*qxhkl+qzhkl*qzhkl);
                        qzhkl=1E-20;
                        dqx=qx-qxhkl;
                        dqy=qy-qyhkl;
                        dqz=qz-qzhkl;
                        dqs1=(dqx*FCC.ax1.x()+dqy*FCC.ax1.y()+dqz*FCC.ax1.z())/(FCC.ax1n_sigx);
                        dqs2=(dqx*FCC.ax2.x()+dqy*FCC.ax2.y()+dqz*FCC.ax2.z())/(FCC.ax2n_sigy);
                        dqs3=(dqx*FCC.ax3.x()+dqy*FCC.ax3.y()+dqz*FCC.ax3.z())/(FCC.ax3n_sigz);
                        x2=dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                        sq=exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*FCC.sig.x()*FCC.sig.y()*FCC.sig.z()/8.);
                    }
                    /* rotated round an oblique axis */
                    else // if ((rotaxphi<>0) and (rotaxphi<>90)) then begin
                    {
                        FCC.qrombchid(FCC.params.length,FCC.ucl1,FCC.ucl2,FCC.ucl3,delta,FCC.uctheta*M_PI/180.0,FCC.ucphi*M_PI/180.0,qx,qy,qz,
                                      FCC.ri11,FCC.ri12,FCC.ri13,FCC.ri21,FCC.ri22,FCC.ri23,FCC.ri31,FCC.ri32,FCC.ri33,
                                      qxhkl,qyhkl,qzhkl,qhkl,FCC.ax1.length(),FCC.ax2.length(),FCC.ax3.length(),
                                      FCC.ax1.x(),FCC.ax1.y(),FCC.ax1.z(),FCC.ax2.x(),FCC.ax2.y(),FCC.ax2.z(),
                                      FCC.ax3.x(),FCC.ax3.y(),FCC.ax3.z(),FCC.sig.x(),FCC.sig.y(),FCC.sig.z(),
                                      FCC.ordis,3,5,7,h,k,l, FCC.params.CR->carr1p, sq );                         /*Z0311=24644*/
                        sq = sq*2*M_PI*qhkl/(2*M_PI*sqrt(M_PI*M_PI*M_PI)*FCC.sig.x()*FCC.sig.y()*FCC.sig.z()/8.0);
                    }
                }
                else if ( FCC.ordis == 14 /* TEST:(x²*y³) */ )
                {
                    D8( qDebug() << "CPU: shp=8(T): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                    double phi = atan2(qyhkl,(qxhkl+eps9));
                    double theta = atan2(sqrt(qxhkl*qxhkl+qyhkl*qyhkl),(qzhkl+eps9));
                    phi = phi*180/M_PI;
                    theta = theta*180/M_PI;
                    FCC.qrombdeltac(FCC.params.length, FCC.params.radius, theta, phi, qx, qy, qz,
                                    qxhkl, qyhkl, qzhkl, qhkl, FCC.ax1.length(), FCC.ax2.length(), FCC.ax3.length(),
                                    FCC.ax1.x(), FCC.ax1.y(), FCC.ax1.z(), FCC.ax2.x(), FCC.ax2.y(), FCC.ax2.z(),
                                    FCC.ax3.x(), FCC.ax3.y(), FCC.ax3.z(), FCC.sig.x(), FCC.sig.y(), FCC.sig.z(),
                                    FCC.ordis,3, /*i0=*/99, /*i1=*/99, 0,0,0, FCC.params.CR->carr1p, sq );
                }
                else // if ( (FCC.ordis!=6) && (FCC.ordis!=7) && (FCC.ordis!=13) && (FCC.dbeta!=0) )    /*Z0311=24676*/
                {   /* other anisotropic cases, rotated around the qhkl-axis direction */                                                             //20210812-D
                    D8( qDebug() << "CPU: shp=8(D): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                    double phi = atan2(qyhkl,(qxhkl+eps6));
                    double theta = atan2(sqrt(qxhkl*qxhkl+qyhkl*qyhkl),(qzhkl+eps6));
                    phi = phi*180/M_PI;
                    theta = theta*180/M_PI;
                    FCC.qrombdeltac(FCC.params.length, FCC.params.radius, theta, phi, qx, qy, qz,
                                    qxhkl, qyhkl, qzhkl, qhkl, FCC.ax1.length(), FCC.ax2.length(), FCC.ax3.length(),
                                    FCC.ax1.x(), FCC.ax1.y(), FCC.ax1.z(), FCC.ax2.x(), FCC.ax2.y(), FCC.ax2.z(),
                                    FCC.ax3.x(), FCC.ax3.y(), FCC.ax3.z(), FCC.sig.x(), FCC.sig.y(), FCC.sig.z(),
                                    FCC.ordis,3, /*i0=*/5, /*i1=*/6,0,0,0, FCC.params.CR->carr1p, sq );
                    sq = sq*2*M_PI*qhkl/FCC.norm;
                }

                psiord = 1;     //{NV}-X/16581
#ifdef undef
                if ( FCC.CheckBoxTwinned )
                {
                    double sqtwin;
                    if ( FCC.ordis==6 || fabs(FCC.params.dbeta) < eps9 )
                    {   /* perfect orientation */
                        D8( qDebug() << "CPU: shp=8(TA): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                        dqx = qx-qxhklt;
                        dqy = qy-qyhklt;
                        dqz = qz-qzhklt;
                        dqs1 = (dqx*FCC.ax1.x()+dqy*FCC.ax1.y()+dqz*FCC.ax1.z())/(FCC.ax1n_sigx);
                        dqs2 = (dqx*FCC.ax2.x()+dqy*FCC.ax2.y()+dqz*FCC.ax2.z())/(FCC.ax2n_sigy);
                        dqs3 = (dqx*FCC.ax3.x()+dqy*FCC.ax3.y()+dqz*FCC.ax3.z())/(FCC.ax3n_sigz);
                        x2 = dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                        sqtwin = exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*FCC.sig.x()*FCC.sig.y()*FCC.sig.z()/8.0);
                        /* sqtwin:=sqtwin*2*pi*qhkl; */
                    }
                    else if ( FCC.ordis==7 )
                    {   /* isotropic orientation */
                        D8( qDebug() << "CPU: shp=8(TB): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                        width = FCC.sig.length() /*sqrt(sigx*sigx+sigy*sigy+sigz*sigz)*/ /3.0;
                        x2 = (q-qhkl)*(q-qhkl)/(width*width);
                        sqtwin = g3*exp(-4*x2/M_PI)/(M_PI*width/2.0);
                    }
                    else if ( FCC.ordis==13 )
                    {   /* fiber pattern */
                        D8( qDebug() << "CPU: shp=8(TC): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                        FCC.qrombchid(FCC.ucl1,FCC.ucl2,FCC.ucl3,delta,FCC.uctheta*M_PI/180.0,FCC.ucphi*M_PI/180.0,qx,qy,qz,
                                      FCC.ri11,FCC.ri12,FCC.ri13,FCC.ri21,FCC.ri22,FCC.ri23,FCC.ri31,FCC.ri32,FCC.ri33,
                                      qxhkl,qyhkl,qzhkl,qhkl,FCC.ax1.length(),FCC.ax2.length(),FCC.ax3.length(),
                                      FCC.ax1.x(),FCC.ax1.y(),FCC.ax1.z(),FCC.ax2.x(),FCC.ax2.y(),FCC.ax2.z(),
                                      FCC.ax3.x(),FCC.ax3.y(),FCC.ax3.z(),FCC.sig.x(),FCC.sig.y(),FCC.sig.z(),
                                      FCC.ordis,3,5,7,h,k,l,sqtwin, threadid );
                        sqtwin = sqtwin*2*M_PI*qhkl/(2*M_PI*sqrt(M_PI*M_PI*M_PI)*FCC.sig.x()*FCC.sig.y()*FCC.sig.z()/8.0);
                    }
                    else // if ( (FCC.ordis!=6) && (FCC.ordis!=7) && (FCC.ordis!=13) && (FCC.dbeta!=0) )
                    {   /* other anistropic cases */
                        D8( qDebug() << "CPU: shp=8(TD): ihex,i" << ihex << i << "ii" << ii << "/" << FCC.peakmax2 );
                        double phi = atan2(qyhklt,(qxhklt+eps9));
                        double theta = atan2(sqrt(qxhklt*qxhklt+qyhklt*qyhklt),(qzhklt+eps9));
                        phi = phi*180/M_PI;
                        theta = theta*180/M_PI;
                        FCC.qrombdeltac( FCC.params.radius, theta, phi, qx, qy, qz,
                                        qxhkl, qyhkl, qzhkl, qhkl, FCC.ax1.length(), FCC.ax2.length(), FCC.ax3.length(),
                                        FCC.ax1.x(), FCC.ax1.y(), FCC.ax1.z(), FCC.ax2.x(), FCC.ax2.y(), FCC.ax2.z(),
                                        FCC.ax3.x(), FCC.ax3.y(), FCC.ax3.z(), FCC.sig.x(), FCC.sig.y(), FCC.sig.z(),
                                        FCC.ordis,3,5,6,0,0,0,sqtwin, threadid );
                        sqtwin = sqtwin*2*M_PI*qhkl/FCC.norm;
                        //qrombdeltac(radius,sigma,dbeta,theta,phi,qx,qy,qz,qxhklt,qyhklt,qzhklt,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,3,5,6,0,0,0,sqtwin);
                        //sqtwin = sqtwin*2*M_PI*qhkl/norm;
                    }
                    sq = FCC.params.ceff*sq+(1-FCC.params.ceff)*sqtwin;
                    //sq = ceff*sq+(1-ceff)*sqtwin;
                }
#endif
                break;  // shp == cbpeakAnisotropicGaussian

            case cbpeakLorentzian /*1*/:
                peaknorm1 = FCC.latpar3(ii,11);				/*Z0311=24690*/
                yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = sqrt(FCC.c1)/(M_PI*width*(1+FCC.c1*x2));
                x2phi = 4*q*q/sqr(FCC.phiwidth);             /*** b-factor ***/
                psiord = g3/(peaknorm1*(1+x2phi*yphi*yphi));
                if ( FCC.CheckBoxTwinned )
                {
                    yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*(1+x2phi*yphi*yphi));
                    psiord = FCC.params.ceff*psiord+(1-FCC.params.ceff)*phiord;
                }
                break;

            case cbpeakGaussian /*2*/:                                                          //20210812-E
                peaknorm1 = FCC.latpar3(ii,11);
                yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));   /*Z0311=24706*/
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = exp(-4*x2/M_PI)/(M_PI*width/2.0);
                x2phi = 4*q*q/(M_PI*sqr(FCC.phiwidth));             /*** a-factor ***/
                psiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
                if ( FCC.CheckBoxTwinned )
                {
                    yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
                    psiord = FCC.params.ceff*psiord+(1-FCC.params.ceff)*phiord;
                }
                break;

            case cbpeakMod1Lorentzian /*Lorentzian1*/:
                peaknorm1 = FCC.latpar3(ii,11);
                yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = 2*sqrt(FCC.c2)/(M_PI*width*(1+FCC.c1*x2)*(1+FCC.c1*x2));
                x2phi = 4*q*q/sqr(FCC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*(1+x2phi*yphi*yphi)*(1+x2phi*yphi*yphi));
                if ( FCC.CheckBoxTwinned )
                {
                    yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*(1+x2phi*yphi*yphi)*(1+x2phi*yphi*yphi));
                    psiord = FCC.params.ceff*psiord+(1-FCC.params.ceff)*phiord;
                }
                break;

            case cbpeakMod2Lorentzian /*Lorentzian2*/:
                peaknorm1 = FCC.latpar3(ii,11);
                yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = sqrt(FCC.c3)/(2*width*exp(3*log(1+FCC.c1*x2)/2.0));
                x2phi = 4*q*q/sqr(FCC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*exp(3*log(1+x2phi*yphi*yphi)/2.0));
                if ( FCC.CheckBoxTwinned )
                {
                    yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*exp(3*log(1+x2phi*yphi*yphi)/2.0));
                    psiord = FCC.params.ceff*psiord+(1-FCC.params.ceff)*phiord;
                }
                break;

            case cbpeakPseudoVoigt /*Voigt*/:
            {
                peaknorm1 = FCC.latpar3(ii,11);
                double peaknorm2 = FCC.latpar3(ii,12);
                yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = FCC.eta*sqrt(FCC.c1)/(M_PI*width*(1+FCC.c1*x2))+(1-FCC.eta)*sqrt(FCC.c0)*exp(-FCC.c0*x2)/(sqrt(M_PI)*width);
                double x2psi = 4*q*q/(M_PI*sqr(FCC.phiwidth));             /*** a-factor ***/
                //double x2psihkl = 4*qhkl*qhkl/(M_PI*sqr(FCC.phiwidth));
                x2phi = 4*q*q/sqr(FCC.phiwidth);                /*** b-factor ***/
                psiord = g3*(FCC.eta*(1/(1+x2phi*yphi*yphi))+(1-FCC.eta)*exp(-x2psi*yphi*yphi))/(FCC.eta*peaknorm1+(1-FCC.eta)*peaknorm2);
                if ( FCC.CheckBoxTwinned )
                {
                    yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3*(FCC.eta*(1/(1+x2phi*yphi*yphi))+(1-FCC.eta)*exp(-x2psi*yphi*yphi))/(FCC.eta*peaknorm1+(1-FCC.eta)*peaknorm2);
                    psiord = FCC.params.ceff*psiord+(1-FCC.params.ceff)*phiord;
                }
                break;
            }

            case cbpeakPearsonVII /*Pearson*/:
                peaknorm1 = FCC.latpar3(ii,11);
                yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = FCC.gamma(FCC.beta)*FCC.c4*2/(FCC.gamma(FCC.beta-0.5)*M_PI*width*exp(FCC.beta*log(1+4*FCC.c4*x2)));
                x2phi = 4*q*q/sqr(FCC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*exp(FCC.beta*log(1+x2phi*yphi*yphi)));
                if ( FCC.CheckBoxTwinned )
                {
                    yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    phiord = g3/(peaknorm1*exp(FCC.beta*log(1+x2phi*yphi*yphi)));
                    psiord = FCC.params.ceff*psiord+(1-FCC.params.ceff)*phiord;
                }
                break;

#ifdef undef
                if ( BurgerG )
                {
                    peaknorm1 = latpar3p^[ii][11];
                    yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                    x2 = (q-qhkl)*(q-qhkl)/(width*width);
                    sq = burger(width,bnu,x2);
                    x2phi = 4*q*q/(M_PI*phiwidth*phiwidth);             /*** a-factor ***/
                    psiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
                    if ( FCC.CheckBoxTwinned )
                    {
                        yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                        phiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
                        psiord = ceff*psiord+(1-ceff)*phiord;
                    }
                }
#endif

            default:
                return 0;
            } // switch shp

            if ( qhkl > 0 ) intensity += sq*mhkl*fhkl*psiord; /*Z0311=24796*/
        }/*2*/  /* of peak-loop */

        /*Z0311=19416:  b:=StrToFloat(EditLatticeb.Text); */
        /*Z0311=19438:  dwfactoriso:=b; */

        /*Z0311=28841 Abschlussberechnungen (izero,base) machen nur hier Sinn, später werden intensity und radintensity nicht mehr verwendet */
        //for i:=1 to imax do
        //intensity = FCC.base + FCC.izero*intensity + conc*radintensity + ifluc/(1+q*q*rfluc*rfluc);
        //for i:=iimin to iimax do
        //radintensity = FCC.base + FCC.izero*radintensity + ifluc/(1+q*q*rfluc*rfluc);

        szqiso = (1+(1)*(8*M_PI*M_PI*M_PI*radintensity/(4*M_PI*sphno*cubevol)-1)*exp(-FCC.params.ucb/*dwfactoriso*/*q*q));   /*Z0311=24804*/
        szq = (1+(fq/pq)*(2*M_PI*intensity/(sphno*cubevol)-1)*exp(-FCC.dwfactor*q*q));   /*Z0311=24805*/
    }
    else // (if ! lattice)
    {
        szq = 1.0;           /*Z=25483*/
        szqiso = 1.0;   // zur Sicherheit
    }

    double retval = FCC.base + FCC.izero*(szq*pq + FCC.iso*szqiso*pqiso) + FCC.ifluc/(1+q*q*FCC.rfluc*FCC.rfluc); /*Z=29549*/
#ifndef __CUDACC__
    if ( retval < -1e6 )
        qDebug() << "szq"<<szq << "pq"<<pq <<"fq"<<fq << "szqiso"<<szqiso << "pqiso"<<pqiso << "q"<<q
                 << "intens"<<intensity << "radint"<<radintensity
                 << "="<<retval;
#else
    //if ( retval < -1e6 )
    //    printf( "szq=%lf pq=%lf fq=%lf szqiso=%lf pqiso=%lf q=%lf intens=%lf radint=%lf erg=%lf\n", szq, pq, fq, szqiso, pqiso, q, intensity, radintensity, retval );
#endif

    return retval;
} /* doIntCalc_GENERIC_q_xyz() */






/*Z0311=8711
    procedure corotations(a,b,c,alpha,beta,gamma: extended;
                          u,v,w,ephi: extended;
                          lat1d,lat2d,lat3d: boolean;
                       var m11,m12,m13,m21,m22,m23,m31,m32,m33: extended;
                       var mtw11,mtw12,mtw13,mtw21,mtw22,mtw23,mtw31,mtw32,mtw33,vvol: extended;
                       var nuvwx,nuvwy,nuvwz,uuvwx,uuvwy,uuvwz,vuvwx,vuvwy,vuvwz: extended;
                       var nhklx,nhkly,nhklz,uhklx,uhkly,uhklz,vhklx,vhkly,vhklz: extended);
    */
void SasCalc_GENERIC_calculation::corotations(double a, double b, double c, double alpha, double beta, double gamma,
                                              double u, double v, double w, double ephi,
                                              bool lat1d, bool lat2d, bool /*lat3d*/,
                                              double &m11, double &m12, double &m13, double &m21, double &m22, double &m23, double &m31, double &m32, double &m33,
                                              double &mtw11, double &mtw12, double &mtw13, double &mtw21, double &mtw22, double &mtw23, double &mtw31, double &mtw32, double &mtw33, double &vvol,
                                              double &nuvwx, double &nuvwy, double &nuvwz, double &uuvwx, double &uuvwy, double &uuvwz, double &vuvwx, double &vuvwy, double &vuvwz,
                                              double &nhklx, double &nhkly, double &nhklz, double &uhklx, double &uhkly, double &uhklz, double &vhklx, double &vhkly, double &vhklz)
{   /*Z0311=8717*/
    //int    i;
    double ca,cb,cg,/*sa,*/sb,sg,/*vol,*/n1,n2,n3,l1,l2,l3,m1,m2,m3;
    double msi11,msi12,msi13,msi21,msi22,msi23,msi31,msi32,msi33,detmsi;
    double ms11,ms12,ms13,ms21,ms22,ms23,ms31,ms32,ms33;
    //double msi11n,msi12n,msi13n,msi21n,msi22n,msi23n,msi31n,msi32n,msi33n,detmsin;
    //double m11n,m12n,m13n,m21n,m22n,m23n,m31n,m32n,m33n;
    double mst11,mst12,mst13,mst21,mst22,mst23,mst31,mst32,mst33,detmst;
    double mt11,mt12,mt13,mt21,mt22,mt23,mt31,mt32,mt33;
    float g11,g12,g13,g21,g22,g23,g31,g32,g33,gi11,gi12,gi13,gi21,gi22,gi23,gi31,gi32,gi33,detg;
    float mi11,mi12,mi13,mi21,mi22,mi23,mi31,mi32,mi33,detm;
    double ri11,ri12,ri13,ri21,ri22,ri23,ri31,ri32,ri33,detr;
    float aax,aay,aaz,bax,bay,baz,cax,cay,caz;
    float aaxt,aayt,aazt,baxt,bayt,bazt,caxt,cayt,cazt;
    float aex,aey,aez,bex,bey,bez,cex,cey,cez; // ,asx,asy,asz,bsx,bsy,bsz,csx,csy,csz;
    float /*aexn,aeyn,aezn,bexn,beyn,bezn,cexn,ceyn,cezn,*/asxn,asyn,aszn,bsxn,bsyn,bszn,csxn,csyn,cszn;
    float aext,aeyt,aezt,bext,beyt,bezt,cext,ceyt,cezt; // ,asxt,asyt,aszt,bsxt,bsyt,bszt,csxt,csyt,cszt;
    //float asxd,asyd,aszd,bsxd,bsyd,bszd,csxd,csyd,cszd;
    float epsi,etheta,len,r11,r12,r13,r21,r22,r23,r31,r32,r33;
    float aedx,aedy,aedz,bedx,bedy,bedz,cedx,cedy,cedz,vols;
    float nxyzx,nxyzy,nxyzz,uxyzx,uxyzy,uxyzz,vxyzx,vxyzy,vxyzz,nxyznx,/*nxyzny,*/nxyznz;
    //double asax,asay,asaz,bsax,bsay,bsaz,csax,csay,csaz,asex,asey,asez,bsex,bsey,bsez,csex,csey,csez;
    double asedx,asedy,asedz,bsedx,bsedy,bsedz,csedx,csedy,csedz;
    double asedxn,asedyn,asedzn,bsedxn,bsedyn,bsedzn,csedxn,csedyn,csedzn;
    double asedxt,asedyt,asedzt,bsedxt,bsedyt,bsedzt,csedxt,csedyt,csedzt;
    double aedxt,aedyt,aedzt,bedxt,bedyt,bedzt,cedxt,cedyt,cedzt;
    double m11DD,m12DD,m13DD,m21DD,m22DD,m23DD,m31DD,m32DD,m33DD;
    //double mi11DD,mi12DD,mi13DD,mi21DD,mi22DD,mi23DD,mi31DD,mi32DD,mi33DD;
    double aexDD,aeyDD,aezDD,bexDD,beyDD,bezDD,cexDD,ceyDD,cezDD;
    double asxDD,asyDD,aszDD,bsxDD,bsyDD,bszDD,csxDD,csyDD,cszDD,area;
    //double aedxDD,aedyDD,aedzDD,bedxDD,bedyDD,bedzDD,cedxDD,cedyDD,cedzDD;
    double asedxDD,asedyDD,asedzDD,bsedxDD,bsedyDD,bsedzDD,csedxDD,csedyDD,csedzDD;
    double m11DL,m12DL,m13DL,m21DL,m22DL,m23DL,m31DL,m32DL,m33DL;
    //double mi11DL,mi12DL,mi13DL,mi21DL,mi22DL,mi23DL,mi31DL,mi32DL,mi33DL;
    //double aexDL,aeyDL,aezDL,bexDL,beyDL,bezDL,cexDL,ceyDL,cezDL;
    double asxDL,asyDL,aszDL,bsxDL,bsyDL,bszDL,csxDL,csyDL,cszDL,area1;
    //double aedxDL,aedyDL,aedzDL,bedxDL,bedyDL,bedzDL,cedxDL,cedyDL,cedzDL;
    double asedxDL,asedyDL,asedzDL,bsedxDL,bsedyDL,bsedzDL,csedxDL,csedyDL,csedzDL;

    //    if ( a == 0 ) a = 1;        // Kantenlängen der Einheitszelle nicht 0
    //    if ( b == 0 ) b = 1;
    //    if ( c == 0 ) c = 1;

    /* unit cell */   /*Z0311=8758*/
    alpha = alpha*M_PI/180.0;   /*Z0311=8759*/
    beta = beta*M_PI/180.0;   /*Z0311=8760*/
    gamma = gamma*M_PI/180.0;   /*Z0311=8761*/
    ephi = ephi*M_PI/180.0;   /*Z0311=8762*/
    ca = cos(alpha);   /*Z0311=8763*/
    cb = cos(beta);   /*Z0311=8764*/
    cg = cos(gamma);   /*Z0311=8765*/
    //sa = sin(alpha);   /*Z0311=8766*/
    sb = sin(beta);   /*Z0311=8767*/
    sg = sin(gamma);   /*Z0311=8768*/
    //vol = a*b*c*sqrt(1.0-ca*ca-cb*cb-cg*cg+2*ca*cb*cg);   /*Z0311=8769*/

    n1 = 0;   /*Z0311=8771*/
    n2 = 0;   /*Z0311=8772*/
    n3 = 1;   /*Z0311=8773*/
    l1 = sb;   /*Z0311=8774*/
    l2 = 0;   /*Z0311=8775*/
    l3 = cb;   /*Z0311=8776*/
    m1 = (cg-cb*ca)/sb;   /*Z0311=8777*/
    m2 = sqrt(1.0-m1*m1-ca*ca);   /*Z0311=8778*/
    m3 = ca;   /*Z0311=8779*/

    /* Msi matrix */   /*Z0311=8781*/
    msi11 = a*l1;         msi12 = a*l2;       msi13 = a*l3;   /*Z0311=8782*/
    msi21 = b*m1;         msi22 = b*m2;       msi23 = b*m3;   /*Z0311=8783*/
    msi31 = c*n1;         msi32 = c*n2;       msi33 = c*n3;   /*Z0311=8784*/

    /* new: transposed Msi matrix */   /*Z0311=8786*/
    //msi11n = a*l1;         msi21n = a*l2;       msi31n = a*l3;   /*Z0311=8787*/
    //msi12n = b*m1;         msi22n = b*m2;       msi32n = b*m3;   /*Z0311=8788*/
    //msi13n = c*n1;         msi23n = c*n2;       msi33n = c*n3;   /*Z0311=8789*/

    /* new: = M matrix */   /*Z0311=8791*/
    //m11n = msi11;          m12n = msi12;       m13n = msi13;   /*Z0311=8792*/
    //m21n = msi21;          m22n = msi22;       m32n = msi23;   /*Z0311=8793*/
    //m31n = msi31;          m23n = msi32;       m33n = msi33;   /*Z0311=8794*/

    /* old: Ms matrix */   /*Z0311=8796*/
    detmsi = msi11*msi22*msi33+msi12*msi23*msi31+msi13*msi21*msi32-msi31*msi22*msi13-msi32*msi23*msi11-msi33*msi21*msi12;   /*Z0311=8797*/
    ms11 = (msi22*msi33-msi23*msi32)/detmsi;   /*Z0311=8798*/
    ms12 = (msi13*msi32-msi12*msi33)/detmsi;   /*Z0311=8799*/
    ms13 = (msi12*msi23-msi13*msi22)/detmsi;   /*Z0311=8800*/
    ms21 = (msi23*msi31-msi21*msi33)/detmsi;   /*Z0311=8801*/
    ms22 = (msi11*msi33-msi13*msi31)/detmsi;   /*Z0311=8802*/
    ms23 = (msi13*msi21-msi11*msi23)/detmsi;   /*Z0311=8803*/
    ms31 = (msi21*msi32-msi22*msi31)/detmsi;   /*Z0311=8804*/
    ms32 = (msi12*msi31-msi11*msi32)/detmsi;   /*Z0311=8805*/
    ms33 = (msi11*msi22-msi12*msi21)/detmsi;   /*Z0311=8806*/

    /* old: Ms matrix, transposed */   /*Z0311=8808*/
    mst11 = ms11;    mst12 = ms21;    mst13 = ms31;   /*Z0311=8809*/
    mst21 = ms12;    mst22 = ms22;    mst23 = ms32;   /*Z0311=8810*/
    mst31 = ms13;    mst32 = ms23;    mst33 = ms33;   /*Z0311=8811*/

    /* old: M matrix, transposed */   /*Z0311=8813*/
    detmst = mst11*mst22*mst33+mst12*mst23*mst31+mst13*mst21*mst32-mst31*mst22*mst13-mst32*mst23*mst11-mst33*mst21*mst12;   /*Z0311=8814*/
    mt11 = (mst22*mst33-mst23*mst32)/detmst;   /*Z0311=8815*/
    mt12 = (mst13*mst32-mst12*mst33)/detmst;   /*Z0311=8816*/
    mt13 = (mst12*mst23-mst13*mst22)/detmst;   /*Z0311=8817*/
    mt21 = (mst23*mst31-mst21*mst33)/detmst;   /*Z0311=8818*/
    mt22 = (mst11*mst33-mst13*mst31)/detmst;   /*Z0311=8819*/
    mt23 = (mst13*mst21-mst11*mst23)/detmst;   /*Z0311=8820*/
    mt31 = (mst21*mst32-mst22*mst31)/detmst;   /*Z0311=8821*/
    mt32 = (mst12*mst31-mst11*mst32)/detmst;   /*Z0311=8822*/
    mt33 = (mst11*mst22-mst12*mst21)/detmst;   /*Z0311=8823*/

    /* old: M matrix */   /*Z0311=8825*/
    /* rxyz=M ruvw */   /*Z0311=8826*/
    m11 = mt11;       m12 = mt21;        m13 = mt31;   /*Z0311=8827*/
    m21 = mt12;       m22 = mt22;        m23 = mt32;   /*Z0311=8828*/
    m31 = mt13;       m32 = mt23;        m33 = mt33;   /*Z0311=8829*/

    m11DD = a;     m12DD = b*cg;    m13DD = 0;   /*Z0311=8831*/
    m21DD = 0;     m22DD = b*sg;    m23DD = 0;   /*Z0311=8832*/
    m31DD = 0;     m32DD = 0;       m33DD = 1;   /*Z0311=8833*/

    m11DL = a;     m12DL = 0;    m13DL = 0;   /*Z0311=8835*/
    m21DL = 0;     m22DL = 1;    m23DL = 0;   /*Z0311=8836*/
    m31DL = 0;     m32DL = 0;    m33DL = 1;   /*Z0311=8837*/

    /* Mi inverse matrix */   /*Z0311=8840*/
    /* ruvw=Mi rxyz */   /*Z0311=8841*/
    detm = m11*m22*m33+m12*m23*m31+m13*m21*m32-m31*m22*m13-m32*m23*m11-m33*m21*m12;   /*Z0311=8842*/
    mi11 = (m22*m33-m23*m32)/detm;   /*Z0311=8843*/
    mi12 = (m13*m32-m12*m33)/detm;   /*Z0311=8844*/
    mi13 = (m12*m23-m13*m22)/detm;   /*Z0311=8845*/
    mi21 = (m23*m31-m21*m33)/detm;   /*Z0311=8846*/
    mi22 = (m11*m33-m13*m31)/detm;   /*Z0311=8847*/
    mi23 = (m13*m21-m11*m23)/detm;   /*Z0311=8848*/
    mi31 = (m21*m32-m22*m31)/detm;   /*Z0311=8849*/
    mi32 = (m12*m31-m11*m32)/detm;   /*Z0311=8850*/
    mi33 = (m11*m22-m12*m21)/detm;   /*Z0311=8851*/

    //mi11DL = 1/a;      mi12DL = 0;        mi13DL = 0;   /*Z0311=8853*/
    //mi21DL = 0;        mi22DL = 1;        mi23DL = 0;   /*Z0311=8854*/
    //mi31DL = 0;        mi32DL = 0;        mi33DL = 1;   /*Z0311=8855*/

    /* base vectors of unit cell */   /*Z0311=8858*/
    /* aa, ba, ca in uvw-system */   /*Z0311=8859*/
    aax = 1;         aay = 0;      aaz = 0;   /*Z0311=8860*/
    bax = 0;         bay = 1;      baz = 0;   /*Z0311=8861*/
    cax = 0;         cay = 0;      caz = 1;   /*Z0311=8862*/

    /* base vectors of twinned unit cell */   /*Z0311=8864*/
    aaxt = 2/3.0;      aayt = 2/3.0;     aazt = -1/3.0;   /*Z0311=8865*/
    baxt = -1/3.0;     bayt = 2/3.0;     bazt = 2/3.0;   /*Z0311=8866*/
    caxt = 2/3.0;      cayt = -1/3.0;    cazt = 2/3.0;   /*Z0311=8867*/

    /* unit vectors in carthesian coordinate system */   /*Z0311=8869*/
    /* ae=M aa, be=M ba, ce=M ca in xyz-system */   /*Z0311=8870*/
    /* Mok */   /*Z0311=8871*/
    /* old */   /*Z0311=8872*/
    aex = m11*aax+m12*aay+m13*aaz;   /*Z0311=8873*/
    aey = m21*aax+m22*aay+m23*aaz;   /*Z0311=8874*/
    aez = m31*aax+m32*aay+m33*aaz;   /*Z0311=8875*/
    bex = m11*bax+m12*bay+m13*baz;   /*Z0311=8876*/
    bey = m21*bax+m22*bay+m23*baz;   /*Z0311=8877*/
    bez = m31*bax+m32*bay+m33*baz;   /*Z0311=8878*/
    cex = m11*cax+m12*cay+m13*caz;   /*Z0311=8879*/
    cey = m21*cax+m22*cay+m23*caz;   /*Z0311=8880*/
    cez = m31*cax+m32*cay+m33*caz;   /*Z0311=8881*/

    /* new */   /*Z0311=8883*/
    //aexn = m11n*aax+m12n*aay+m13n*aaz;   /*Z0311=8884*/
    //aeyn = m21n*aax+m22n*aay+m23n*aaz;   /*Z0311=8885*/
    //aezn = m31n*aax+m32n*aay+m33n*aaz;   /*Z0311=8886*/
    //bexn = m11n*bax+m12n*bay+m13n*baz;   /*Z0311=8887*/
    //beyn = m21n*bax+m22n*bay+m23n*baz;   /*Z0311=8888*/
    //bezn = m31n*bax+m32n*bay+m33n*baz;   /*Z0311=8889*/
    //cexn = m11n*cax+m12n*cay+m13n*caz;   /*Z0311=8890*/
    //ceyn = m21n*cax+m22n*cay+m23n*caz;   /*Z0311=8891*/
    //cezn = m31n*cax+m32n*cay+m33n*caz;   /*Z0311=8892*/

    aexDD = m11DD*aax+m12DD*aay+m13DD*aaz;   /*Z0311=8894*/
    aeyDD = m21DD*aax+m22DD*aay+m23DD*aaz;   /*Z0311=8895*/
    aezDD = m31DD*aax+m32DD*aay+m33DD*aaz;   /*Z0311=8896*/
    bexDD = m11DD*bax+m12DD*bay+m13DD*baz;   /*Z0311=8897*/
    beyDD = m21DD*bax+m22DD*bay+m23DD*baz;   /*Z0311=8898*/
    bezDD = m31DD*bax+m32DD*bay+m33DD*baz;   /*Z0311=8899*/
    cexDD = m11DD*cax+m12DD*cay+m13DD*caz;   /*Z0311=8900*/
    ceyDD = m21DD*cax+m22DD*cay+m23DD*caz;   /*Z0311=8901*/
    cezDD = m31DD*cax+m32DD*cay+m33DD*caz;   /*Z0311=8902*/

    //aexDL = m11DL*aax+m12DL*aay+m13DL*aaz;   /*Z0311=8904*/
    //aeyDL = m21DL*aax+m22DL*aay+m23DL*aaz;   /*Z0311=8905*/
    //aezDL = m31DL*aax+m32DL*aay+m33DL*aaz;   /*Z0311=8906*/
    //bexDL = m11DL*bax+m12DL*bay+m13DL*baz;   /*Z0311=8907*/
    //beyDL = m21DL*bax+m22DL*bay+m23DL*baz;   /*Z0311=8908*/
    //bezDL = m31DL*bax+m32DL*bay+m33DL*baz;   /*Z0311=8909*/
    //cexDL = m11DL*cax+m12DL*cay+m13DL*caz;   /*Z0311=8910*/
    //ceyDL = m21DL*cax+m22DL*cay+m23DL*caz;   /*Z0311=8911*/
    //cezDL = m31DL*cax+m32DL*cay+m33DL*caz;   /*Z0311=8912*/

    aext = m11*aaxt+m12*aayt+m13*aazt;   /*Z0311=8914*/
    aeyt = m21*aaxt+m22*aayt+m23*aazt;   /*Z0311=8915*/
    aezt = m31*aaxt+m32*aayt+m33*aazt;   /*Z0311=8916*/
    bext = m11*baxt+m12*bayt+m13*bazt;   /*Z0311=8917*/
    beyt = m21*baxt+m22*bayt+m23*bazt;   /*Z0311=8918*/
    bezt = m31*baxt+m32*bayt+m33*bazt;   /*Z0311=8919*/
    cext = m11*caxt+m12*cayt+m13*cazt;   /*Z0311=8920*/
    ceyt = m21*caxt+m22*cayt+m23*cazt;   /*Z0311=8921*/
    cezt = m31*caxt+m32*cayt+m33*cazt;   /*Z0311=8922*/

    /* old: reciprocal space vector in carthesian coordinate system */   /*Z0311=8924*/
    /* ase, bse, cse in xyz-system */   /*Z0311=8925*/
    /* Mok */   /*Z0311=8926*/
    vvol = aex*(bey*cez-bez*cey)+aey*(bez*cex-bex*cez)+aez*(bex*cey-bey*cex);   /*Z0311=8927*/
    //asx = (bey*cez-bez*cey)/vvol;   /*Z0311=8928*/
    //asy = (bez*cex-bex*cez)/vvol;   /*Z0311=8929*/
    //asz = (bex*cey-bey*cex)/vvol;   /*Z0311=8930*/
    //bsx = (aez*cey-aey*cez)/vvol;   /*Z0311=8931*/
    //bsy = (aex*cez-aez*cex)/vvol;   /*Z0311=8932*/
    //bsz = (aey*cex-aex*cey)/vvol;   /*Z0311=8933*/
    //csx = (aey*bez-aez*bey)/vvol;   /*Z0311=8934*/
    //csy = (aez*bex-aex*bez)/vvol;   /*Z0311=8935*/
    //csz = (aex*bey-aey*bex)/vvol;   /*Z0311=8936*/

    area = a*b*sg;   /*Z0311=8938*/
    asxDD = beyDD/area;   /*Z0311=8939*/
    asyDD = -bexDD/area;   /*Z0311=8940*/
    aszDD = 0;   /*Z0311=8941*/
    bsxDD = -aeyDD/area;   /*Z0311=8942*/
    bsyDD = aexDD/area;   /*Z0311=8943*/
    bszDD = 0;   /*Z0311=8944*/
    csxDD = 0;   /*Z0311=8945*/
    csyDD = 0;   /*Z0311=8946*/
    cszDD = 1;   /*Z0311=8947*/

    area1 = a;   /*Z0311=8949*/
    asxDL = 1/area1;   /*Z0311=8950*/
    asyDL = 0;   /*Z0311=8951*/
    aszDL = 0;   /*Z0311=8952*/
    bsxDL = 0;   /*Z0311=8953*/
    bsyDL = 1;   /*Z0311=8954*/
    bszDL = 0;   /*Z0311=8955*/
    csxDL = 0;   /*Z0311=8956*/
    csyDL = 0;   /*Z0311=8957*/
    cszDL = 1;   /*Z0311=8958*/

    //asxt = (beyt*cezt-bezt*ceyt)/vvol;   /*Z0311=8960*/
    //asyt = (bezt*cext-bext*cezt)/vvol;   /*Z0311=8961*/
    //aszt = (bext*ceyt-beyt*cext)/vvol;   /*Z0311=8962*/
    //bsxt = (aezt*ceyt-aeyt*cezt)/vvol;   /*Z0311=8963*/
    //bsyt = (aext*cezt-aezt*cext)/vvol;   /*Z0311=8964*/
    //bszt = (aeyt*cext-aext*ceyt)/vvol;   /*Z0311=8965*/
    //csxt = (aeyt*bezt-aezt*beyt)/vvol;   /*Z0311=8966*/
    //csyt = (aezt*bext-aext*bezt)/vvol;   /*Z0311=8967*/
    //cszt = (aext*beyt-aeyt*bext)/vvol;   /*Z0311=8968*/

    /* new: G metric matrix in xyz-coordinates */   /*Z0311=8970*/
    g11 = aex;     g12 = bex;     g13 = cex;   /*Z0311=8971*/
    g21 = aey;     g22 = bey;     g23 = cey;   /*Z0311=8972*/
    g31 = aez;     g32 = bez;     g33 = cez;   /*Z0311=8973*/
    /*Z0311=8974*/
    /* old G matrix */   /*Z0311=8975*/

#ifdef PascalComment
    g11 = aex*aex+aey*aey+aez*aez;   /*Z0311=8976*/
    g12 = aex*bex+aey*bey+aez*bez;   /*Z0311=8977*/
    g13 = aex*cex+aey*cey+aez*cez;   /*Z0311=8978*/
    g21 = bex*aex+bey*aey+bez*aez;   /*Z0311=8979*/
    g22 = bex*bex+bey*bey+bez*bez;   /*Z0311=8980*/
    g23 = bex*cex+bey*cey+bez*cez;   /*Z0311=8981*/
    g31 = cex*aex+cey*aey+cez*aez;   /*Z0311=8982*/
    g32 = cex*bex+cey*bey+cez*bez;   /*Z0311=8983*/
    g33 = cex*cex+cey*cey+cez*cez;
#endif

    /* Gs inverse metric matrix */   /*Z0311=8986*/
    detg = g11*g22*g33+g12*g23*g31+g13*g21*g32-g31*g22*g13-g32*g23*g11-g33*g21*g12;   /*Z0311=8987*/
    gi11 = (g22*g33-g23*g32)/detg;   /*Z0311=8988*/
    gi12 = (g13*g32-g12*g33)/detg;   /*Z0311=8989*/
    gi13 = (g12*g23-g13*g22)/detg;   /*Z0311=8990*/
    gi21 = (g23*g31-g21*g33)/detg;   /*Z0311=8991*/
    gi22 = (g11*g33-g13*g31)/detg;   /*Z0311=8992*/
    gi23 = (g13*g21-g11*g23)/detg;   /*Z0311=8993*/
    gi31 = (g21*g32-g22*g31)/detg;   /*Z0311=8994*/
    gi32 = (g12*g31-g11*g32)/detg;   /*Z0311=8995*/
    gi33 = (g11*g22-g12*g21)/detg;   /*Z0311=8996*/

    /* new: reciprocal space vector in carthesian coordinate system */   /*Z0311=8998*/
    /* ase, bse, cse in xyz-system */   /*Z0311=8999*/
    asxn = gi11;   asyn = gi12;   aszn = gi13;   /*Z0311=9000*/
    bsxn = gi21;   bsyn = gi22;   bszn = gi23;   /*Z0311=9001*/
    csxn = gi31;   csyn = gi32;   cszn = gi33;   /*Z0311=9002*/

    /* a*,b*,c* reciprocal space vectors */   /*Z0311=9004*/
    /* in uvw-space */   /*Z0311=9005*/

#ifdef PascalComment
    asax = g11*aax+g12*aay+g13*aaz;   /*Z0311=9006*/
    asay = g21*aax+g22*aay+g23*aaz;   /*Z0311=9007*/
    asaz = g31*aax+g32*aay+g33*aaz;   /*Z0311=9008*/
    bsax = g11*bax+g12*bay+g13*baz;   /*Z0311=9009*/
    bsay = g21*bax+g22*bay+g23*baz;   /*Z0311=9010*/
    bsaz = g31*bax+g32*bay+g33*baz;   /*Z0311=9011*/
    csax = g11*cax+g12*cay+g13*caz;   /*Z0311=9012*/
    csay = g21*cax+g22*cay+g23*caz;   /*Z0311=9013*/
    csaz = g31*cax+g32*cay+g33*caz;
#endif

    /* a*, b*, c* reciprocal space vectors */   /*Z0311=9016*/
    /* in xyz-space */   /*Z0311=9017*/

#ifdef PascalComment
    asex = m11*asax+m12*asay+m13*asaz;   /*Z0311=9018*/
    asey = m21*asax+m22*asay+m23*asaz;   /*Z0311=9019*/
    asez = m31*asax+m32*asay+m33*asaz;   /*Z0311=9020*/
    bsex = m11*bsax+m12*bsay+m13*bsaz;   /*Z0311=9021*/
    bsey = m21*bsax+m22*bsay+m23*bsaz;   /*Z0311=9022*/
    bsez = m31*bsax+m32*bsay+m33*bsaz;   /*Z0311=9023*/
    csex = m11*csax+m12*csay+m13*csaz;   /*Z0311=9024*/
    csey = m21*csax+m22*csay+m23*csaz;   /*Z0311=9025*/
    csez = m31*csax+m32*csay+m33*csaz;
#endif

    /* nuvw-vector || beam */   /*Z0311=9029*/
    /* nuvw in unit cell uvw-system */   /*Z0311=9030*/
    nuvwx = u;   /*Z0311=9031*/
    nuvwy = v;   /*Z0311=9032*/
    nuvwz = w;   /*Z0311=9033*/

    /* nhkl-vector */   /*Z0311=9035*/
    /* nhkl=G nuvw in unit cell uvw-system */   /*Z0311=9036*/

#ifdef PascalComment
    nhklx = g11*nuvwx+g12*nuvwy+g13*nuvwz;   /*Z0311=9037*/
    nhkly = g21*nuvwx+g22*nuvwy+g23*nuvwz;   /*Z0311=9038*/
    nhklz = g31*nuvwx+g32*nuvwy+g33*nuvwz;
#endif

    /* nhkl-vector */   /*Z0311=9041*/
    /* nxyz=M nuvw in xyz-system */   /*Z0311=9042*/
    /* with (a,b,c)=(2,1,1) and (u,v,w)=(2,1,1) this is a (ua,vb,wc)=(4,1,1)-vector */   /*Z0311=9043*/
    nxyzx = m11*nuvwx+m12*nuvwy+m13*nuvwz;   /*Z0311=9044*/
    nxyzy = m21*nuvwx+m22*nuvwy+m23*nuvwz;   /*Z0311=9045*/
    nxyzz = m31*nuvwx+m32*nuvwy+m33*nuvwz;   /*Z0311=9046*/

    /* unit nxyz = director */   /*Z0311=9048*/
    /* in xyz-system */   /*Z0311=9049*/
    len = sqrt(1.0*nxyzx*nxyzx+nxyzy*nxyzy+nxyzz*nxyzz);   /*Z0311=9050*/ //220908
    nxyznx = nxyzx/len;   /*Z0311=9051*/
    //nxyzny = nxyzy/len;   /*Z0311=9052*/
    nxyznz = nxyzz/len;   /*Z0311=9053*/

    /* R rotation matrix */   /*Z0311=9055*/
    /* in xyz-system */   /*Z0311=9056*/
    etheta = acos(-nxyznz);   /*Z0311=9057*/
    epsi = asin(nxyznx/sin(etheta));   /*Z0311=9058*/

    r11 = cos(epsi)*cos(ephi)-sin(epsi)*cos(etheta)*sin(ephi);   /*Z0311=9060*/
    r12 = -cos(epsi)*sin(ephi)-sin(epsi)*cos(etheta)*cos(ephi);   /*Z0311=9061*/
    r13 = sin(epsi)*sin(etheta);   /*Z0311=9062*/
    r21 = sin(epsi)*cos(ephi)+cos(epsi)*cos(etheta)*sin(ephi);   /*Z0311=9063*/
    r22 = -sin(epsi)*sin(ephi)+cos(epsi)*cos(etheta)*cos(ephi);   /*Z0311=9064*/
    r23 = -cos(epsi)*sin(etheta);   /*Z0311=9065*/
    r31 = sin(etheta)*sin(ephi);   /*Z0311=9066*/
    r32 = sin(etheta)*cos(ephi);   /*Z0311=9067*/
    r33 = cos(etheta);   /*Z0311=9068*/

    /* Ri inverse rotation matrix */   /*Z0311=9070*/
    /* in xyz-system */   /*Z0311=9071*/
    detr = r11*r22*r33+r12*r23*r31+r13*r21*r32-r31*r22*r13-r32*r23*r11-r33*r21*r12;   /*Z0311=9072*/
    ri11 = (r22*r33-r23*r32)/detr;   /*Z0311=9073*/
    ri12 = (r13*r32-r12*r33)/detr;   /*Z0311=9074*/
    ri13 = (r12*r23-r13*r22)/detr;   /*Z0311=9075*/
    ri21 = (r23*r31-r21*r33)/detr;   /*Z0311=9076*/
    ri22 = (r11*r33-r13*r31)/detr;   /*Z0311=9077*/
    ri23 = (r13*r21-r11*r23)/detr;   /*Z0311=9078*/
    ri31 = (r21*r32-r22*r31)/detr;   /*Z0311=9079*/
    ri32 = (r12*r31-r11*r32)/detr;   /*Z0311=9080*/
    ri33 = (r11*r22-r12*r21)/detr;   /*Z0311=9081*/

    /* rotated base vectors a,b,c in carthesian coordinate system */   /*Z0311=9083*/
    /* aed=Ri ae, bed=Ri be, ced=Ri ce in xyz-system */   /*Z0311=9084*/
    /* needed for calculation of fiber pattern */   /*Z0311=9085*/
    /* Mok */   /*Z0311=9086*/
    aedx = ri11*aex+ri12*aey+ri13*aez;   /*Z0311=9087*/
    aedy = ri21*aex+ri22*aey+ri23*aez;   /*Z0311=9088*/
    aedz = ri31*aex+ri32*aey+ri33*aez;   /*Z0311=9089*/
    bedx = ri11*bex+ri12*bey+ri13*bez;   /*Z0311=9090*/
    bedy = ri21*bex+ri22*bey+ri23*bez;   /*Z0311=9091*/
    bedz = ri31*bex+ri32*bey+ri33*bez;   /*Z0311=9092*/
    cedx = ri11*cex+ri12*cey+ri13*cez;   /*Z0311=9093*/
    cedy = ri21*cex+ri22*cey+ri23*cez;   /*Z0311=9094*/
    cedz = ri31*cex+ri32*cey+ri33*cez;   /*Z0311=9095*/

    //aedxDD = ri11*aexDD+ri12*aeyDD+ri13*aezDD;   /*Z0311=9097*/
    //aedyDD = ri21*aexDD+ri22*aeyDD+ri23*aezDD;   /*Z0311=9098*/
    //aedzDD = ri31*aexDD+ri32*aeyDD+ri33*aezDD;   /*Z0311=9099*/
    //bedxDD = ri11*bexDD+ri12*beyDD+ri13*bezDD;   /*Z0311=9100*/
    //bedyDD = ri21*bexDD+ri22*beyDD+ri23*bezDD;   /*Z0311=9101*/
    //bedzDD = ri31*bexDD+ri32*beyDD+ri33*bezDD;   /*Z0311=9102*/
    //cedxDD = ri11*cexDD+ri12*ceyDD+ri13*cezDD;   /*Z0311=9103*/
    //cedyDD = ri21*cexDD+ri22*ceyDD+ri23*cezDD;   /*Z0311=9104*/
    //cedzDD = ri31*cexDD+ri32*ceyDD+ri33*cezDD;   /*Z0311=9105*/

    //aedxDL = ri11*aexDL+ri12*aeyDL+ri13*aezDL;   /*Z0311=9107*/
    //aedyDL = ri21*aexDL+ri22*aeyDL+ri23*aezDL;   /*Z0311=9108*/
    //aedzDL = ri31*aexDL+ri32*aeyDL+ri33*aezDL;   /*Z0311=9109*/
    //bedxDL = ri11*bexDL+ri12*beyDL+ri13*bezDL;   /*Z0311=9110*/
    //bedyDL = ri21*bexDL+ri22*beyDL+ri23*bezDL;   /*Z0311=9111*/
    //bedzDL = ri31*bexDL+ri32*beyDL+ri33*bezDL;   /*Z0311=9112*/
    //cedxDL = ri11*cexDL+ri12*ceyDL+ri13*cezDL;   /*Z0311=9113*/
    //cedyDL = ri21*cexDL+ri22*ceyDL+ri23*cezDL;   /*Z0311=9114*/
    //cedzDL = ri31*cexDL+ri32*ceyDL+ri33*cezDL;   /*Z0311=9115*/

    aedxt = ri11*aext+ri12*aeyt+ri13*aezt;   /*Z0311=9117*/
    aedyt = ri21*aext+ri22*aeyt+ri23*aezt;   /*Z0311=9118*/
    aedzt = ri31*aext+ri32*aeyt+ri33*aezt;   /*Z0311=9119*/
    bedxt = ri11*bext+ri12*beyt+ri13*bezt;   /*Z0311=9120*/
    bedyt = ri21*bext+ri22*beyt+ri23*bezt;   /*Z0311=9121*/
    bedzt = ri31*bext+ri32*beyt+ri33*bezt;   /*Z0311=9122*/
    cedxt = ri11*cext+ri12*ceyt+ri13*cezt;   /*Z0311=9123*/
    cedyt = ri21*cext+ri22*ceyt+ri23*cezt;   /*Z0311=9124*/
    cedzt = ri31*cext+ri32*ceyt+ri33*cezt;   /*Z0311=9125*/

    /* rotated reciprocal space vectors a*,b*,c* in carthesian coordinate system */   /*Z0311=9127*/
    /* calculated from cross-product of aed,bed,ced */   /*Z0311=9128*/
    /* into output matrix */   /*Z0311=9129*/
    /* Mok */   /*Z0311=9130*/
    asedx = (bedy*cedz-bedz*cedy)/vvol;   /*Z0311=9131*/
    asedy = (bedz*cedx-bedx*cedz)/vvol;   /*Z0311=9132*/
    asedz = (bedx*cedy-bedy*cedx)/vvol;   /*Z0311=9133*/
    bsedx = (aedz*cedy-aedy*cedz)/vvol;   /*Z0311=9134*/
    bsedy = (aedx*cedz-aedz*cedx)/vvol;   /*Z0311=9135*/
    bsedz = (aedy*cedx-aedx*cedy)/vvol;   /*Z0311=9136*/
    csedx = (aedy*bedz-aedz*bedy)/vvol;   /*Z0311=9137*/
    csedy = (aedz*bedx-aedx*bedz)/vvol;   /*Z0311=9138*/
    csedz = (aedx*bedy-aedy*bedx)/vvol;   /*Z0311=9139*/

    /* new: a*_r=R a* */   /*Z0311=9141*/
    /* output */   /*Z0311=9142*/
    asedxn = ri11*asxn+ri12*asyn+ri13*aszn;   /*Z0311=9143*/
    asedyn = ri21*asxn+ri22*asyn+ri23*aszn;   /*Z0311=9144*/
    asedzn = ri31*asxn+ri32*asyn+ri33*aszn;   /*Z0311=9145*/
    bsedxn = ri11*bsxn+ri12*bsyn+ri13*bszn;   /*Z0311=9146*/
    bsedyn = ri21*bsxn+ri22*bsyn+ri23*bszn;   /*Z0311=9147*/
    bsedzn = ri31*bsxn+ri32*bsyn+ri33*bszn;   /*Z0311=9148*/
    csedxn = ri11*csxn+ri12*csyn+ri13*cszn;   /*Z0311=9149*/
    csedyn = ri21*csxn+ri22*csyn+ri23*cszn;   /*Z0311=9150*/
    csedzn = ri31*csxn+ri32*csyn+ri33*cszn;   /*Z0311=9151*/

    asedxDD = ri11*asxDD+ri12*asyDD+ri13*aszDD;   /*Z0311=9153*/
    asedyDD = ri21*asxDD+ri22*asyDD+ri23*aszDD;   /*Z0311=9154*/
    asedzDD = ri31*asxDD+ri32*asyDD+ri33*aszDD;   /*Z0311=9155*/
    bsedxDD = ri11*bsxDD+ri12*bsyDD+ri13*bszDD;   /*Z0311=9156*/
    bsedyDD = ri21*bsxDD+ri22*bsyDD+ri23*bszDD;   /*Z0311=9157*/
    bsedzDD = ri31*bsxDD+ri32*bsyDD+ri33*bszDD;   /*Z0311=9158*/
    csedxDD = ri11*csxDD+ri12*csyDD+ri13*cszDD;   /*Z0311=9159*/
    csedyDD = ri21*csxDD+ri22*csyDD+ri23*cszDD;   /*Z0311=9160*/
    csedzDD = ri31*csxDD+ri32*csyDD+ri33*cszDD;   /*Z0311=9161*/

    asedxDL = ri11*asxDL+ri12*asyDL+ri13*aszDL;   /*Z0311=9163*/
    asedyDL = ri21*asxDL+ri22*asyDL+ri23*aszDL;   /*Z0311=9164*/
    asedzDL = ri31*asxDL+ri32*asyDL+ri33*aszDL;   /*Z0311=9165*/
    bsedxDL = ri11*bsxDL+ri12*bsyDL+ri13*bszDL;   /*Z0311=9166*/
    bsedyDL = ri21*bsxDL+ri22*bsyDL+ri23*bszDL;   /*Z0311=9167*/
    bsedzDL = ri31*bsxDL+ri32*bsyDL+ri33*bszDL;   /*Z0311=9168*/
    csedxDL = ri11*csxDL+ri12*csyDL+ri13*cszDL;   /*Z0311=9169*/
    csedyDL = ri21*csxDL+ri22*csyDL+ri23*cszDL;   /*Z0311=9170*/
    csedzDL = ri31*csxDL+ri32*csyDL+ri33*cszDL;   /*Z0311=9171*/

    asedxt = (bedyt*cedzt-bedzt*cedyt)/vvol;   /*Z0311=9173*/
    asedyt = (bedzt*cedxt-bedxt*cedzt)/vvol;   /*Z0311=9174*/
    asedzt = (bedxt*cedyt-bedyt*cedxt)/vvol;   /*Z0311=9175*/
    bsedxt = (aedzt*cedyt-aedyt*cedzt)/vvol;   /*Z0311=9176*/
    bsedyt = (aedxt*cedzt-aedzt*cedxt)/vvol;   /*Z0311=9177*/
    bsedzt = (aedyt*cedxt-aedxt*cedyt)/vvol;   /*Z0311=9178*/
    csedxt = (aedyt*bedzt-aedzt*bedyt)/vvol;   /*Z0311=9179*/
    csedyt = (aedzt*bedxt-aedxt*bedzt)/vvol;   /*Z0311=9180*/
    csedzt = (aedxt*bedyt-aedyt*bedxt)/vvol;   /*Z0311=9181*/

    /* a*, b*, c* rotated */   /*Z0311=9183*/
    /* in xyz-space */   /*Z0311=9184*/

#ifdef PascalComment
    asedx = r11*asex+r12*asey+r13*asez;   /*Z0311=9185*/
    asedy = r21*asex+r22*asey+r23*asez;   /*Z0311=9186*/
    asedz = r31*asex+r32*asey+r33*asez;   /*Z0311=9187*/
    bsedx = r11*bsex+r12*bsey+r13*bsez;   /*Z0311=9188*/
    bsedy = r21*bsex+r22*bsey+r23*bsez;   /*Z0311=9189*/
    bsedz = r31*bsex+r32*bsey+r33*bsez;   /*Z0311=9190*/
    csedx = r11*csex+r12*csey+r13*csez;   /*Z0311=9191*/
    csedy = r21*csex+r22*csey+r23*csez;   /*Z0311=9192*/
    csedz = r31*csex+r32*csey+r33*csez;
#endif

    /* n,u,v axis vectors in Carthesian coordinate system */   /*Z0311=9196*/
    /* in xyz-system */   /*Z0311=9197*/
    nxyzx = 0;     nxyzy = 0;     nxyzz = -1;   /*Z0311=9198*/
    uxyzx = 1;     uxyzy = 0;     uxyzz = 0;   /*Z0311=9199*/
    vxyzx = 0;     vxyzy = 1;     vxyzz = 0;   /*Z0311=9200*/

    /* rotated n,u,v axis vectors */   /*Z0311=9202*/
    /* nd=R n, ud=R u, vd=R v  in xyz-system */   /*Z0311=9203*/
    double nxyzdx = r11*nxyzx+r12*nxyzy+r13*nxyzz;   /*Z0311=9204*/
    double nxyzdy = r21*nxyzx+r22*nxyzy+r23*nxyzz;   /*Z0311=9205*/
    double nxyzdz = r31*nxyzx+r32*nxyzy+r33*nxyzz;   /*Z0311=9206*/
    double uxyzdx = r11*uxyzx+r12*uxyzy+r13*uxyzz;   /*Z0311=9207*/
    double uxyzdy = r21*uxyzx+r22*uxyzy+r23*uxyzz;   /*Z0311=9208*/
    double uxyzdz = r31*uxyzx+r32*uxyzy+r33*uxyzz;   /*Z0311=9209*/
    double vxyzdx = r11*vxyzx+r12*vxyzy+r13*vxyzz;   /*Z0311=9210*/
    double vxyzdy = r21*vxyzx+r22*vxyzy+r23*vxyzz;   /*Z0311=9211*/
    double vxyzdz = r31*vxyzx+r32*vxyzy+r33*vxyzz;   /*Z0311=9212*/

    /* rotated n,u,v axis vectors */   /*Z0311=9214*/
    /* nuvw=Mi nd, uuvw=Mi ud, vuvw=Mi vd  in unit cell uvw-system */   /*Z0311=9215*/
    /* needed to indicate unit cell directions <uvw> in 2D-pattern */   /*Z0311=9216*/
    /* output */   /*Z0311=9217*/
    nuvwx = mi11*nxyzdx+mi12*nxyzdy+mi13*nxyzdz;   /*Z0311=9218*/
    nuvwy = mi21*nxyzdx+mi22*nxyzdy+mi23*nxyzdz;   /*Z0311=9219*/
    nuvwz = mi31*nxyzdx+mi32*nxyzdy+mi33*nxyzdz;   /*Z0311=9220*/
    uuvwx = mi11*uxyzdx+mi12*uxyzdy+mi13*uxyzdz;   /*Z0311=9221*/
    uuvwy = mi21*uxyzdx+mi22*uxyzdy+mi23*uxyzdz;   /*Z0311=9222*/
    uuvwz = mi31*uxyzdx+mi32*uxyzdy+mi33*uxyzdz;   /*Z0311=9223*/
    vuvwx = mi11*vxyzdx+mi12*vxyzdy+mi13*vxyzdz;   /*Z0311=9224*/
    vuvwy = mi21*vxyzdx+mi22*vxyzdy+mi23*vxyzdz;   /*Z0311=9225*/
    vuvwz = mi31*vxyzdx+mi32*vxyzdy+mi33*vxyzdz;   /*Z0311=9226*/

    /* rotated n,u,v axis vectors */   /*Z0311=9228*/
    /* nhkl=G nuvw, uhkl=G uuvw, vhkl=G vuvw  in unit cell uvw-system */   /*Z0311=9229*/
    /* needed to indicate reciprocal space directions <hkl> in 2D-pattern */   /*Z0311=9230*/
    vols = nuvwx*(uuvwy*vuvwz-uuvwz*vuvwy)+nuvwy*(uuvwz*vuvwx-uuvwx*vuvwz)+nuvwz*(uuvwx*vuvwy-uuvwy*vuvwx);   /*Z0311=9231*/
    nhklx = (uuvwy*vuvwz-uuvwz*vuvwy)/vols;   /*Z0311=9232*/
    nhkly = (uuvwz*vuvwx-uuvwx*vuvwz)/vols;   /*Z0311=9233*/
    nhklz = (uuvwx*vuvwy-uuvwy*vuvwx)/vols;   /*Z0311=9234*/
    uhklx = (nuvwz*vuvwy-nuvwy*vuvwz)/vols;   /*Z0311=9235*/
    uhkly = (nuvwx*vuvwz-nuvwz*vuvwx)/vols;   /*Z0311=9236*/
    uhklz = (nuvwy*vuvwx-nuvwx*vuvwy)/vols;   /*Z0311=9237*/
    vhklx = (nuvwy*uuvwz-nuvwz*uuvwy)/vols;   /*Z0311=9238*/
    vhkly = (nuvwz*uuvwx-nuvwx*uuvwz)/vols;   /*Z0311=9239*/
    vhklz = (nuvwx*uuvwy-nuvwy*uuvwx)/vols;   /*Z0311=9240*/

#ifdef PascalComment
    nhklx = g11*nuvwx+g12*nuvwy+g13*nuvwz;   /*Z0311=9242*/
    nhkly = g21*nuvwx+g22*nuvwy+g23*nuvwz;   /*Z0311=9243*/
    nhklz = g31*nuvwx+g32*nuvwy+g33*nuvwz;   /*Z0311=9244*/
    uhklx = g11*uuvwx+g12*uuvwy+g13*uuvwz;   /*Z0311=9245*/
    uhkly = g21*uuvwx+g22*uuvwy+g23*uuvwz;   /*Z0311=9246*/
    uhklz = g31*uuvwx+g32*uuvwy+g33*uuvwz;   /*Z0311=9247*/
    vhklx = g11*vuvwx+g12*vuvwy+g13*vuvwz;   /*Z0311=9248*/
    vhkly = g21*vuvwx+g22*vuvwy+g23*vuvwz;   /*Z0311=9249*/
    vhklz = g31*vuvwx+g32*vuvwy+g33*vuvwz;
#endif

    /* rotated reciprocal space base vectors */   /*Z0311=9252*/
    /* aed=R ae, bed=R be, ced=R ce in xyz-system */   /*Z0311=9253*/
    //asxd = r11*asx+r12*asy+r13*asz;   /*Z0311=9254*/
    //asyd = r21*asx+r22*asy+r23*asz;   /*Z0311=9255*/
    //aszd = r31*asx+r32*asy+r33*asz;   /*Z0311=9256*/
    //bsxd = r11*bsx+r12*bsy+r13*bsz;   /*Z0311=9257*/
    //bsyd = r21*bsx+r22*bsy+r23*bsz;   /*Z0311=9258*/
    //bszd = r31*bsx+r32*bsy+r33*bsz;   /*Z0311=9259*/
    //csxd = r11*csx+r12*csy+r13*csz;   /*Z0311=9260*/
    //csyd = r21*csx+r22*csy+r23*csz;   /*Z0311=9261*/
    //cszd = r31*csx+r32*csy+r33*csz;   /*Z0311=9262*/

    /* output matrix */   /*Z0311=9264*/
    /* m11:=ri11;   m12:=ri12;   m13:=ri13; */   /*Z0311=9265*/
    /* m21:=ri21;   m22:=ri22;   m23:=ri23; */   /*Z0311=9266*/
    /* m31:=ri31;   m32:=ri32;   m33:=ri33; */   /*Z0311=9267*/

    /* m11:=asxd;   m12:=bsxd;   m13:=csxd; */   /*Z0311=9269*/
    /* m21:=asyd;   m22:=bsyd;   m23:=csyd; */   /*Z0311=9270*/
    /* m31:=aszd;   m32:=bszd;   m33:=cszd; */   /*Z0311=9271*/

    m11 = asedx;   m21 = bsedx;   m31 = csedx;   /*Z0311=9273*/
    m12 = asedy;   m22 = bsedy;   m32 = csedy;   /*Z0311=9274*/
    m13 = asedz;   m23 = bsedz;   m33 = csedz;   /*Z0311=9275*/

    if ( lat2d )
    {   /*Z0311=9277*/
        m11 = asedxDD;   m21 = bsedxDD;   m31 = csedxDD;   /*Z0311=9278*/
        m12 = asedyDD;   m22 = bsedyDD;   m32 = csedyDD;   /*Z0311=9279*/
        m13 = asedzDD;   m23 = bsedzDD;   m33 = csedzDD;   /*Z0311=9280*/
        vvol = area;   /*Z0311=9281*/
    }   /*Z0311=9282*/

    if ( lat1d )
    {   /*Z0311=9284*/
        m11 = asedxDL;   m21 = bsedxDL;   m31 = csedxDL;   /*Z0311=9285*/
        m12 = asedyDL;   m22 = bsedyDL;   m32 = csedyDL;   /*Z0311=9286*/
        m13 = asedzDL;   m23 = bsedzDL;   m33 = csedzDL;   /*Z0311=9287*/
        vvol = area1;   /*Z0311=9288*/
    }   /*Z0311=9289*/

    mtw11 = asedxt;   mtw21 = bsedxt;   mtw31 = csedxt;   /*Z0311=9291*/
    mtw12 = asedyt;   mtw22 = bsedyt;   mtw32 = csedyt;   /*Z0311=9292*/
    mtw13 = asedzt;   mtw23 = bsedzt;   mtw33 = csedzt;   /*Z0311=9293*/

    /* test output */   /*Z0311=9295*/
    nhklx = asedxn;   /*Z0311=9296*/
    nhkly = asedyn;   /*Z0311=9297*/
    nhklz = asedzn;   /*Z0311=9298*/
    uhklx = bsedxn;   /*Z0311=9299*/
    uhkly = bsedyn;   /*Z0311=9300*/
    uhklz = bsedzn;   /*Z0311=9301*/
    vhklx = csedxn;   /*Z0311=9302*/
    vhkly = csedyn;   /*Z0311=9303*/
    vhklz = csedzn;   /*Z0311=9304*/
}   /*Z0311=9305*/



// myarray wird genutzt, um diverse Parameter als "Struktur" zu übergeben!
void CLASSLIB::coefficients( double l, double r, double rm, double sigmal, double sigma,
                             double epsi, double alfa, double dbeta, double theta, double phi,
                             int part, int dim, int nmax, int ordis, int cs,
                             double *myarray /*Z0311=9310*/,
                             int &cho1/*orcase*/,
                             double &por, double &order, double &norm, double &lim1, double &lim2,
                             double &lim3, double &lim4, double &lim5, double &lim6, double &lim7,
                             double &lim8, double &lim9,
                             double &lim1f, double &lim2f, double &lim3f, double &lim4f, double &lim5f,
                             double &lim6f, double &lim7f, double &lim8f, double &lim9f,
                             _carrXX *CR )  // Für GPU ist ein Zeiger besser als viele!
                             //double *carr1p, double *carr2p, double *carr3p,
                             //double *carr4p, double *carr5p, double *carr6p,
                             //double *carr7p, double *carr8p, double *carr9p, //: CoeffArrayType  /*Z0311=9314*/,
                             //double *carr1f, double *carr2f, double *carr3f,
                             //double *carr4f, double *carr5f, double *carr6f,
                             //double *carr7f, double *carr8f, double *carr9f)  //: ArrayImax2D )
{   /*Z0311=9317*/ /*Z=9365*/

    static const int nnmax = 120; // 500;  /*Z0311=9321*/
    const double min = 1E-30;  /*Z0311=9322*/

    int    i,j,k,ll,m,n,ns,ms,n1,n2,n3,n4,n5,n6,n7,n8,n9,n1f,/*n2f,n3f,*/n4f,n5f,n6f,n7f,n8f,n9f,/*nmax2,*/cnt,inmax;
    int    /*ii,jj,kk,pp,qq,rr,iis,iie,jjs,jje,pps,ppe,qqs,qqe,*/ncell;
    double z,zz,zl,zzl,xlz,xl2z,xrz,xr2z,b1s,/*v,ee0,ee1,gb1s,preg1,pz2v,pz2v1,pz2v2,com,*/binsum/*,binsum1*/;
    double a1,intl,p11,p12,p13,p21,p22,p23,p31,p32,p33,p,vvm,eps,area,vol;
    long double /*u1z,u1n,u2,u3,u4,u5,*/xrmz,xrm2z,sump,sump1,sump2,sumi,sumf,radv,aell,bell,cell;
    double carr2i[nnmax+1],fsum[nnmax+1]; //: Array[0..2*nnmax] of extended;  /*Z0311=9329*/;
    double lqm[10/*3*nnmax+1*/],qqm[10/*3*nnmax+1*/],phim[3*nnmax+1],llm[3*nnmax+1],rrm[3*nnmax+1],ppm[3*nnmax+1],
           a1m[3*nnmax+1]/*,a2m[3*nnmax+1]*/; //: Array[0..3*nnmax] of extended;  /*Z=9378*/
    double /*uell[3*nnmax+1],*/u1ell[3*nnmax+1],u2ell[3*nnmax+1],u3ell[3*nnmax+1],/*vell[3*nnmax+1],*/v1ell[3*nnmax+1],
           v2ell[3*nnmax+1],v3ell[3*nnmax+1],gell[3*nnmax+1],g3ell[3*nnmax+1]; //: Array[0..3*nnmax] of extended;  /*Z=9378*/
    double intlar[2*nnmax+1][2*nnmax+1]; //: Array[0..2*nnmax,0..2*nnmax] of extended;  /*Z0311=9331*/;
    bool   search1,/*search2,*/search3,search4,search5,search6/*,search7,search8,search9*/;
    float  philiph,philipt,phiax,phiin,phiout,rad,lliph,llipt,lin,lout,len,rmax;  /*Z0311=9333*/;
    float  /*xradm,xrzm,*/x1zm,x12zm;  /*Z0311=9334*/;

    // Arrays von oben mit erhöhter Genauigkeit
    long double gam3[nnmax+1], fkv[nnmax+1], fk2v[nnmax+1], pn[nnmax+1],xrn[nnmax+1], xln[nnmax+1], z12v[nnmax+1], z12vl[nnmax+1];

    // Arrays von oben, von denen aber nur das letzte Element genutzt wird.
    // Können also als normale Variablen geschrieben werden um Speicher zu sparen
    long double xrmn_n, fkvm_n;

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
    p11=1; p12=1; p13=1; p21=1; p22=1; p23=1; p31=1; p32=1; p33=1;

    /*begin*/  /*Z0311=9336*/
    z = (1-sqr(sigma))/sqr(sigma);  /*Z0311=9337*/ /*Z=9386*/
    zl = (1-sqr(sigmal))/sqr(sigmal);  /*Z0311=9338*/
    zz = z;  /*Z0311=9339*/
    zzl = zl;  /*Z0311=9340*/
    DBG( qDebug() << "coefficients()" << "z=" << z << "zl=" << zl; )
    // z:=z+2*dim;          /* z-average */ */  /*Z0311=9341*/
    /* zl:=zl+2*dim; */  /*Z0311=9342*/
    //nmax2 = ceil(nmax/2.0);  /*Z0311=9343*/
    p = r/rm;  /*Z0311=9344*/
    eps = l/r;  /*Z0311=9345*/

    xlz = l/(2.*(zzl+1));  /*Z0311=9347*/
    xl2z = xlz*xlz;  /*Z0311=9348*/
    xrz = r/(2.*(zz+1));  /*Z0311=9349*/
    xr2z = xrz*xrz;  /*Z0311=9350*/
    xrmz = rm/(2.*(zz+1));  /*Z0311=9351*/
    xrm2z = xrmz*xrmz;  /*Z0311=9352*/

    /* ellipsoid semiaxis */  /*Z0311=9354*/
    aell = r;  /*Z0311=9355*/
    bell = l;  /*Z0311=9356*/
    cell = rm;  /*Z0311=9357*/

    xln[0] = 1;  /*Z0311=9359*/
    xrn[0] = 1;  /*Z0311=9360*/
    xrmn_n = 1;  /*Z0311=9361*/
    pn[0] = 1;  /*Z0311=9362*/
    z12v[0] = 1;  /*Z0311=9363*/
    z12vl[0] = 1;  /*Z0311=9364*/
    fsum[0] = 1;  /*Z0311=9365*/

    /* factor for cross-sectional formfactors */  /*Z0311=9368*/ /*Z=9417*/
    if ( dim==1 ) b1s = 2;         /* cylinder, cross-section is disk */  /*Z0311=9369*/
    if ( dim==2 ) b1s = 3/2.0;       /* disk, cross-section is cylinder */  /*Z0311=9370*/
    if ( dim==3 ) b1s = 5/2.0;       /* sphere */  /*Z0311=9371*/
    if ( dim==4 ) b1s = 3/2.0;       /* cube */  /*Z0311=9372*/

    /* start values for recursive parameters */  /*Z0311=9374*/
    //b1sv[0] = 1;  /*Z0311=9375*/
    fkv[0] = 1;  /*Z0311=9376*/
    fk2v[0] = 1;  /*Z0311=9377*/
    //e1[0] = 1;  /*Z0311=9378*/
    //gam1[0] = sqrt(M_PI);  /*Z0311=9379*/
    //gam2[0] = 1.0;  /*Z0311=9380*/
    gam3[0] = sqrt(M_PI)/2.0;  /*Z0311=9381*/

    /* initialize */  /*Z0311=9383*/
    /* for n:=0 to 10000 do begin */  /*Z0311=9384*/
    /*    carr1pm[n]:=1; */  /*Z0311=9385*/
    /*    carr2pm[n]:=1; */  /*Z0311=9386*/
    /* end; */  /*Z0311=9387*/
    for ( n=0; n<coeffarray_len; n++ )
    {  /*Z0311=9388*/
        CR->carr1p[n] = 1;      CR->carr1f[n] = 1;  /*Z0311=9389*/
        CR->carr2p[n] = 1;      CR->carr2f[n] = 1;  /*Z0311=9390*/
        CR->carr3p[n] = 1;      CR->carr3f[n] = 1;  /*Z0311=9391*/
        CR->carr4p[n] = 1;      CR->carr4f[n] = 1;  /*Z0311=9392*/
        CR->carr5p[n] = 1;      CR->carr5f[n] = 1;  /*Z0311=9393*/
        CR->carr6p[n] = 1;      CR->carr6f[n] = 1;  /*Z0311=9394*/
        CR->carr7p[n] = 1;      CR->carr7f[n] = 1;  /*Z0311=9395*/
        CR->carr8p[n] = 1;      CR->carr8f[n] = 1;  /*Z0311=9396*/
        CR->carr9p[n] = 1;      CR->carr9f[n] = 1;  /*Z0311=9397*/
        carr2i[n] = 1;  /*Z0311=9398*/
    }  /*Z0311=9399*/
    for ( n=0; n<imax2d_len; n++ ) // imax2d_len=130+1
    {  /*Z0311=9400*/
        for ( m=0; m<imax2d_len; m++ )
        {  /*Z0311=9401*/
            CR->carr11pm[n][m] = 1;  /*Z0311=9402, init*/
            CR->carr22pm[n][m] = 1;  /*Z0311=9403*/
        }  /*Z0311=9404*/
    }  /*Z0311=9405*/

    /* multi-shell or liposome structure parameters */  /*Z=9456, nicht angepasst*/
    if ( cs==3 )
    {   /*Z0311=9408*/
        philiph = myarray[12];     /* water */  /*Z0311=9409*/
        philipt = myarray[13];     /* bilayer */  /*Z0311=9410*/

        rad = myarray[1];          /* vesicle inner radius */  /*Z0311=9412*/
        lliph = myarray[7];        /* water */  /*Z0311=9413*/
        llipt = myarray[8];        /* bilayer */  /*Z0311=9414*/

        len = lliph+llipt;  /*Z0311=9416*/
        ncell = round(myarray[4]);  /*Z0311=9417*/
        rmax = rad+ncell*len;  /*Z0311=9418*/

        lqm[1] = lliph;  /*Z0311=9420*/
        lqm[2] = lqm[1]+llipt;  /*Z0311=9421*/

        for ( i=1; i<=2; i++ ) qqm[i] = lqm[i]/len;  /*Z0311=9423*/

        phim[1] = philiph;       /* vesicle interior */  /*Z0311=9425*/
        llm[1] = 0.0;  /*Z0311=9426*/
        rrm[1] = rad+llm[1];  /*Z0311=9427*/
        ppm[1] = 1.0;  /*Z0311=9428*/

        radv = rrm[1];        /* vesicle radius */  /*Z0311=9430*/
        cnt = 1;  /*Z0311=9431*/
        for ( i=1; i<=ncell; i++ )
        {  /*Z0311=9432*/
            phim[cnt+1] = philipt;         /* bilayer */  /*Z0311=9433*/
            phim[cnt+2] = philiph;         /* water */  /*Z0311=9434*/
            llm[cnt+1] = (i-1+qqm[1])*len;  /*Z0311=9435*/
            llm[cnt+2] = (i-1+qqm[2])*len;  /*Z0311=9436*/
            rrm[cnt+1] = radv+llm[cnt+1];  /*Z0311=9437*/
            rrm[cnt+2] = radv+llm[cnt+2];  /*Z0311=9438*/
            ppm[cnt+1] = rrm[cnt+1]/rad;  /*Z0311=9439*/
            ppm[cnt+2] = rrm[cnt+2]/rad;  /*Z0311=9440*/
            cnt = cnt+2;  /*Z0311=9441*/
        }  /*Z0311=9442*/
        inmax = cnt;  /*Z0311=9443*/
        phim[cnt+1] = 0.0;  /*Z0311=9444*/

        //xradm = rad;  /*Z0311=9446*/
        //xrzm = xradm/(z+1);  /*Z0311=9447*/
        x1zm = rad/(2.*(z+1));  /*Z0311=9448*/
        x12zm = x1zm*x1zm;  /*Z0311=9449*/
        // nicht verwendet: xmax = q*rmax;  /*Z0311=9450*/

        for ( i=1; i<=inmax; i++ )
        {  /*Z0311=9452*/           //220908 - alle pow() auf double,double anpassen
            if ( part==0 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],3.); /* spheres */  /*Z0311=9453*/
            if ( part==1 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],2.); /* cylinders */  /*Z0311=9454*/
            if ( part==2 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],1.); /* disks */  /*Z0311=9455*/
            CR->carr7p[i] = ppm[i];  /*Z0311=9456*/
            CR->carr3p[i] = llm[i];  /*Z0311=9457*/
            CR->carr5p[i] = a1m[i];  /*Z0311=9458*/
        }  /*Z0311=9459*/

        fkvm_n = 1;  /*Z0311=9461*/
        gam3[0] = sqrt(M_PI)/2.0;  /*Z0311=9462*/
        for ( n=1; n<=nmax; n++ )
        {  /*Z0311=9463*/
            if ( part==0 )
            {   /* spheres */  /*Z0311=9464*/
                fkvm_n = fkvm_n*n;  /*Z0311=9465*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9466*/
                CR->carr6p[n] = (n+3/2.0)*gam3[n]*fkvm_n*4/(3.0*sqrt(M_PI));  /*Z0311=9467*/
            }  /*Z0311=9468*/
            if ( part==1 )
            {   /* cylinders */  /*Z0311=9469*/
                fkvm_n = fkvm_n*n;  /*Z0311=9470*/
                CR->carr6p[n] = (n+1)*fkvm_n*fkvm_n;  /*Z0311=9471*/
            }  /*Z0311=9472*/
            if ( part==2 )
            {   /* disks */  /*Z0311=9473*/
                fkvm_n = fkvm_n*n;  /*Z0311=9474*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9475*/
                CR->carr6p[n] = gam3[n]*fkvm_n*2.0/sqrt(M_PI);  /*Z0311=9476*/
            }  /*Z0311=9477*/
        }  /*Z0311=9478*/

        vvm = 0;  /*Z0311=9480*/
        for ( i=1; i<=inmax; i++ )
        {  /*Z0311=9481*/
            for ( j=1; j<=inmax; j++ ) vvm = vvm+a1m[i]*a1m[j];  /*Z0311=9482*/
        }  /*Z0311=9483*/

        myarray[14] = inmax;  /*Z0311=9485*/
        myarray[15] = vvm;  /*Z0311=9486*/
        myarray[16] = rmax;  /*Z0311=9487*/
    }  /*Z0311=9488*/


    /* myelin structure parameters */  /*Z=9540, nicht angepasst*/
    if ( cs==4 )
    {   /*Z0311=9492*/
        philiph = myarray[12];     /* head group */  /*Z0311=9493*/
        philipt = myarray[13];     /* tail */  /*Z0311=9494*/
        phiax = myarray[9];        /* axon */  /*Z0311=9495*/
        phiin = myarray[10];       /* intra cell */  /*Z0311=9496*/
        phiout = myarray[11];      /* extra cell */  /*Z0311=9497*/

        rad = myarray[1];          /* vesicle inner radius */  /*Z0311=9499*/
        lliph = myarray[7];        /* head group */  /*Z0311=9500*/
        llipt = myarray[8];        /* tail */  /*Z0311=9501*/
        lin = myarray[6];          /* intra cell */  /*Z0311=9502*/
        lout = myarray[5];         /* extra cell */  /*Z0311=9503*/

        len = lout+2*(2*lliph+llipt)+lin;  /*Z0311=9505*/
        ncell = round(myarray[4]);  /*Z0311=9506*/
        rmax = rad+ncell*len;  /*Z0311=9507*/

        lqm[1] = lout;  /*Z0311=9509*/
        lqm[2] = lqm[1]+lliph;  /*Z0311=9510*/
        lqm[3] = lqm[2]+llipt;  /*Z0311=9511*/
        lqm[4] = lqm[3]+lliph;  /*Z0311=9512*/
        lqm[5] = lqm[4]+lin;  /*Z0311=9513*/
        lqm[6] = lqm[5]+lliph;  /*Z0311=9514*/
        lqm[7] = lqm[6]+llipt;  /*Z0311=9515*/
        lqm[8] = lqm[7]+lliph;  /*Z0311=9516*/

        for ( i=1; i<=8; i++ ) qqm[i] = lqm[i]/len;  /*Z0311=9518*/

        phim[1] = phiax;       /* vesicle interior */  /*Z0311=9520*/
        llm[1] = 0.0;  /*Z0311=9521*/
        rrm[1] = rad+llm[1];  /*Z0311=9522*/
        ppm[1] = 1.0;  /*Z0311=9523*/
        phim[2] = philiph;     /* vesicle bilayer: head group */  /*Z0311=9524*/
        llm[2] = lliph;  /*Z0311=9525*/
        rrm[2] = rad+llm[2];  /*Z0311=9526*/
        ppm[2] = rrm[2]/rad;  /*Z0311=9527*/
        phim[3] = philipt;     /* vesicle bilayer: tail group */  /*Z0311=9528*/
        llm[3] = llm[2]+llipt;  /*Z0311=9529*/
        rrm[3] = rad+llm[3];  /*Z0311=9530*/
        ppm[3] = rrm[3]/rad;  /*Z0311=9531*/
        phim[4] = philiph;     /* vesicle bilayer: head group */  /*Z0311=9532*/
        llm[4] = llm[3]+lliph;  /*Z0311=9533*/
        rrm[4] = rad+llm[4];  /*Z0311=9534*/
        ppm[4] = rrm[4]/rad;  /*Z0311=9535*/

        radv = rrm[4];        /* vesicle radius + bilayer */  /*Z0311=9537*/
        cnt = 4;  /*Z0311=9538*/
        for ( i=1; i<=ncell; i++ )
        {   /*Z0311=9539*/
            phim[cnt+1] = phiout;          /* extra cell */  /*Z0311=9540*/
            phim[cnt+2] = philiph;         /* head group */  /*Z0311=9541*/
            phim[cnt+3] = philipt;         /* tail group */  /*Z0311=9542*/
            phim[cnt+4] = philiph;         /* head group */  /*Z0311=9543*/
            phim[cnt+5] = phiin;           /* intra cell */  /*Z0311=9544*/
            phim[cnt+6] = philiph;         /* head group */  /*Z0311=9545*/
            phim[cnt+7] = philipt;         /* tail group */  /*Z0311=9546*/
            phim[cnt+8] = philiph;         /* head group */  /*Z0311=9547*/
            llm[cnt+1] = (i-1+qqm[1])*len;  /*Z0311=9548*/
            llm[cnt+2] = (i-1+qqm[2])*len;  /*Z0311=9549*/
            llm[cnt+3] = (i-1+qqm[3])*len;  /*Z0311=9550*/
            llm[cnt+4] = (i-1+qqm[4])*len;  /*Z0311=9551*/
            llm[cnt+5] = (i-1+qqm[5])*len;  /*Z0311=9552*/
            llm[cnt+6] = (i-1+qqm[6])*len;  /*Z0311=9553*/
            llm[cnt+7] = (i-1+qqm[7])*len;  /*Z0311=9554*/
            llm[cnt+8] = (i-1+qqm[8])*len;  /*Z0311=9555*/
            rrm[cnt+1] = radv+llm[cnt+1];  /*Z0311=9556*/
            rrm[cnt+2] = radv+llm[cnt+2];  /*Z0311=9557*/
            rrm[cnt+3] = radv+llm[cnt+3];  /*Z0311=9558*/
            rrm[cnt+4] = radv+llm[cnt+4];  /*Z0311=9559*/
            rrm[cnt+5] = radv+llm[cnt+5];  /*Z0311=9560*/
            rrm[cnt+6] = radv+llm[cnt+6];  /*Z0311=9561*/
            rrm[cnt+7] = radv+llm[cnt+7];  /*Z0311=9562*/
            rrm[cnt+8] = radv+llm[cnt+8];  /*Z0311=9563*/
            ppm[cnt+1] = rrm[cnt+1]/rad;  /*Z0311=9564*/
            ppm[cnt+2] = rrm[cnt+2]/rad;  /*Z0311=9565*/
            ppm[cnt+3] = rrm[cnt+3]/rad;  /*Z0311=9566*/
            ppm[cnt+4] = rrm[cnt+4]/rad;  /*Z0311=9567*/
            ppm[cnt+5] = rrm[cnt+5]/rad;  /*Z0311=9568*/
            ppm[cnt+6] = rrm[cnt+6]/rad;  /*Z0311=9569*/
            ppm[cnt+7] = rrm[cnt+7]/rad;  /*Z0311=9570*/
            ppm[cnt+8] = rrm[cnt+8]/rad;  /*Z0311=9571*/
            cnt = cnt+8;  /*Z0311=9572*/
        }   /*Z0311=9573*/
        inmax = cnt;  /*Z0311=9574*/
        phim[cnt+1] = 0.0;  /*Z0311=9575*/

        //xradm = rad;  /*Z0311=9577*/
        //xrzm = xradm/(z+1);  /*Z0311=9578*/
        x1zm = rad/(2.*(z+1));  /*Z0311=9579*/
        x12zm = x1zm*x1zm;  /*Z0311=9580*/
        // nicht verwendet: xmax = q*rmax;  /*Z0311=9581*/

        for ( i=1; i<=inmax; i++ )
        {   /*Z0311=9583*/
            if ( part==0 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],3.); /* spheres */  /*Z0311=9584*/
            if ( part==1 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],2.); /* cylinders */  /*Z0311=9585*/
            if ( part==2 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],1.); /* disks */  /*Z0311=9586*/
            CR->carr7p[i] = ppm[i];  /*Z0311=9587*/
            CR->carr3p[i] = llm[i];  /*Z0311=9588*/
            CR->carr5p[i] = a1m[i];  /*Z0311=9589*/
        }   /*Z0311=9590*/

        fkvm_n = 1;  /*Z0311=9592*/
        gam3[0] = sqrt(M_PI)/2.0;  /*Z0311=9593*/
        for ( n=1; n<=nmax; n++ )
        {   /*Z0311=9594*/
            if ( part==0 )
            {   /* spheres */  /*Z0311=9595*/
                fkvm_n = fkvm_n*n;  /*Z0311=9596*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9597*/
                CR->carr6p[n] = (n+3/2.0)*gam3[n]*fkvm_n*4/(3*sqrt(M_PI));  /*Z0311=9598*/
            }   /*Z0311=9599*/
            if ( part==1 )
            {   /* cylinders */  /*Z0311=9600*/
                fkvm_n = fkvm_n*n;  /*Z0311=9601*/
                CR->carr6p[n] = (n+1)*fkvm_n*fkvm_n;  /*Z0311=9602*/
            }   /*Z0311=9603*/
            if ( part==2 )
            {   /* disks */  /*Z0311=9604*/
                fkvm_n = fkvm_n*n;  /*Z0311=9605*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9606*/
                CR->carr6p[n] = gam3[n]*fkvm_n*2/sqrt(M_PI);  /*Z0311=9607*/
            }   /*Z0311=9608*/
        }   /*Z0311=9609*/

        vvm = 0;  /*Z0311=9611*/
        for ( i=1; i<=inmax; i++ )
        {   /*Z0311=9612*/
            for ( j=1; j<=inmax; j++ ) vvm = vvm+a1m[i]*a1m[j];  /*Z0311=9613*/
        }   /*Z0311=9614*/

        myarray[14] = inmax;  /*Z0311=9616*/
        myarray[15] = vvm;  /*Z0311=9617*/
        myarray[16] = rmax;  /*Z0311=9618*/
    }   /*Z0311=9619*/


    search1 = true;  /*Z0311=9624*/ /*Z=9673*/
    //search2 = true;  /*Z0311=9625*/
    search3 = true;  /*Z0311=9626*/
    search4 = true;  /*Z0311=9627*/
    search5 = true;  /*Z0311=9628*/
    search6 = true;  /*Z0311=9629*/
    //search7 = true;  /*Z0311=9630*/
    //search8 = true;  /*Z0311=9631*/
    //search9 = true;  /*Z0311=9632*/
    n1 = nmax;      n1f = nmax;  /*Z0311=9633*/
    n2 = nmax;      //n2f = nmax;  /*Z0311=9634*/
    n3 = nmax;      //n3f = nmax;  /*Z0311=9635*/
    n4 = nmax;      n4f = nmax;  /*Z0311=9636*/
    n5 = nmax;      n5f = nmax;  /*Z0311=9637*/
    n6 = nmax;      n6f = nmax;  /*Z0311=9638*/
    n7 = nmax;      n7f = nmax;  /*Z0311=9639*/
    n8 = nmax;      n8f = nmax;  /*Z0311=9640*/
    n9 = nmax;      n9f = nmax;  /*Z0311=9641*/

    /* orientation case */  /*Z0311=9643*/ /*Z=9692*/
    cho1 = 1;                                     /* general */  /*Z0311=9644*/
    if ( (phi== 0) && (theta==90) ) cho1 = 2;     /* x-axis */   /*Z0311=9645*/
    if ( (phi==90) && (theta==90) ) cho1 = 3;     /* y-axis */   /*Z0311=9646*/
    if ( (phi== 0) && (theta== 0) ) cho1 = 4;     /* z-axis */   /*Z0311=9647*/
    if ( (phi==90) && (theta== 0) ) cho1 = 4;     /* z-axis */   /*Z0311=9648*/

#ifndef __CUDACC__
    //DBG( qDebug() << "coefficients()" << "ordis"<<ordis << "dim"<<dim << "part"<<part << "nmax"<<nmax << "cs"<<cs << "orcase" << cho1; )
#endif
    //Debug: coefficients ordis 0="CBOrdis.Gaussian"
    //                    dim 1=im Vorfeld aus CBPart.Cylinder bestimmt
    //                    part 1="CBPart.Cylinder"
    //                    nmax 120="120" festgeschieben
    //                    cs 0="CBInt.homogeneous" - geht von 0 bis 2, hier wird aber 3 und mehr genutzt TODO


    /*** isotropic case for spheres ***/  /*Z=9700, nicht angepasst*/
    if ( dim==3 )
    {   /*Z0311=9652*/
        norm = 1;  /*Z0311=9653*/
        order = 0;  /*Z0311=9654*/
        /* homogeneous */  /*Z0311=9655*/
        if ( cs==0 )
        {   /*Z0311=9656*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z0311=9657*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  /*Z0311=9658*/
                fkv[n] = fkv[n-1]*n;  /*Z0311=9659*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9660*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */  /*Z0311=9661*/
                xrn[n] = -xrn[n-1]*xr2z;  /*Z0311=9662*/
                /* P(q)-coefficient */  /*Z0311=9663*/
                CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.,1.0*n)*z12v[n]*xrn[n]/(2.*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  /*Z0311=9664*/
                /* F(q)-coefficient */  /*Z0311=9665*/
                binsum = 0.0;  /*Z0311=9666*/
                for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  /*Z0311=9667*/
                CR->carr4f[n] = 9*M_PI*xrn[n-1]*binsum/16.0;  /*Z0311=9668*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z0311=9669*/
                    if ( n<n4 ) n4 = n;  /*Z0311=9670*/
                }   /*Z0311=9671*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z0311=9672*/
                    if ( n<n4f ) n4f = n;  /*Z0311=9673*/
                }   /*Z0311=9674*/
            }   /*Z0311=9675*/
            goto Label99;  /*Z0311=9676*/
        }   /* of homogeneous */  /*Z0311=9677*/

        /* core/shell */  /*Z0311=9679*/
        if ( cs==1 )
        {   /*Z0311=9680*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z0311=9681*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  /*Z0311=9682*/
                fkv[n] = fkv[n-1]*n;  /*Z0311=9683*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9684*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */  /*Z0311=9685*/
                xrn[n] = -xrn[n-1]*xr2z;  /*Z0311=9686*/
                xrmn_n = -xrmn_n*xrm2z;  /*Z0311=9687*/
                pn[n] = pn[n-1]*p*p;  /*Z0311=9688*/
                /* P(q)-coefficients */  /*Z0311=9689*/
                sump = 0.0;  /*Z0311=9690*/
                for ( m=0; m<=n; m++ )
                {   /*Z0311=9691*/
                    sump = sump+pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  /*Z0311=9692*/
                }   /*Z0311=9693*/
                CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.,1.0*n)*z12v[n]*xrn[n]/(2.*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  /*Z0311=9694*/
                CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;  /*Z0311=9695*/
                CR->carr6p[n] = CR->carr4p[n]/pn[n];  /*Z0311=9696*/
                /* F(q)-coefficients */  /*Z0311=9697*/
                sump = 0.0;  /*Z0311=9698*/
                for ( m=0; m<=n; m++ )
                {   /*Z0311=9699*/
                    sump = sump+z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  /*Z0311=9700*/
                }   /*Z0311=9701*/
                CR->carr4f[n] = (9*M_PI/16.0)*xrn[n-1]*sump;  /*Z0311=9702*/
                sump = 0.0;  /*Z0311=9703*/
                for ( m=0; m<=n; m++ )
                {   /*Z0311=9704*/
                    sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  /*Z0311=9705*/
                }   /*Z0311=9706*/
                CR->carr5f[n] = (9*M_PI/16.0)*xrmn_n*sump;  /*Z0311=9707*/
                CR->carr6f[n] = CR->carr4f[n]/pn[n];  /*Z0311=9708*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z0311=9709*/
                    if ( n<n4 ) n4 = n;  /*Z0311=9710*/
                }   /*Z0311=9711*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z0311=9712*/
                    if ( n<n5 ) n5 = n;  /*Z0311=9713*/
                }   /*Z0311=9714*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z0311=9715*/
                    if ( n<n6 ) n6 = n;  /*Z0311=9716*/
                }   /*Z0311=9717*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z0311=9718*/
                    if ( n<n4f ) n4f = n;  /*Z0311=9719*/
                }   /*Z0311=9720*/
                if ( fabs(CR->carr5f[n])<min )
                {   /*Z0311=9721*/
                    if ( n<n5f ) n5f = n;  /*Z0311=9722*/
                }   /*Z0311=9723*/
                if ( fabs(CR->carr6f[n])<min )
                {   /*Z0311=9724*/
                    if ( n<n6f ) n6f = n;  /*Z0311=9725*/
                }   /*Z0311=9726*/
            }   /*Z0311=9727*/
            goto Label99;  /*Z0311=9728*/
        }   /* of core/shell */  /*Z0311=9729*/

        /* inhomogeneous core/shell */  /*Z0311=9731*/
        if ( cs==2 )
        {   /*Z0311=9732*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z0311=9733*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  /*Z0311=9734*/
                fkv[n] = fkv[n-1]*n;  /*Z0311=9735*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9736*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */  /*Z0311=9737*/
                xrn[n] = -xrn[n-1]*xr2z;  /*Z0311=9738*/
                xrmn_n = -xrmn_n*xrm2z;  /*Z0311=9739*/
                pn[n] = pn[n-1]*p*p;  /*Z0311=9740*/
                /* P(q)-coefficients */  /*Z0311=9741*/
                CR->carr1p[n] = 9*sqrt(M_PI)*pow(4.,1.0*n)*z12v[n]*xrn[n]/(2.*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  /*Z0311=9742*/
                sump = 0.0;  /*Z0311=9743*/
                sump1 = 0.0;  /*Z0311=9744*/
                for ( m=0; m<=n; m++ )
                {   /*Z0311=9745*/
                    sumi = 1/((n-m+3/2.0)*gam3[n-m]*(m+3/2.0-alfa/2.0)*gam3[m]*fkv[m]*fkv[n-m]);  /*Z0311=9746*/
                    sump = sump+pn[n-m]*sumi;  /*Z0311=9747*/
                    sump1 = sump1+sumi;  /*Z0311=9748*/
                }   /*Z0311=9749*/
                CR->carr2p[n] = (3*M_PI*(3-alfa)/16.0)*z12v[n]*xrmn_n*sump;  /*Z0311=9750*/
                CR->carr3p[n] = (3*M_PI*(3-alfa)/16.0)*z12v[n]*xrn[n-1]*sump1;  /*Z0311=9751*/
                sump = 0.0;  /*Z0311=9752*/
                sump1 = 0.0;  /*Z0311=9753*/
                for ( m=0; m<=n; m++ )
                {   /*Z0311=9754*/
                    sumi = 1/((n-m+3/2.0-alfa/2.0)*(m+3/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  /*Z0311=9755*/
                    sump = sump+sumi;  /*Z0311=9756*/
                    sump1 = sump1+pn[n-m]*sumi;  /*Z0311=9757*/
                }   /*Z0311=9758*/
                CR->carr4p[n] = ((3-alfa)*(3-alfa)*M_PI/16.0)*z12v[n]*xrmn_n*sump;  /*Z0311=9759*/
                CR->carr5p[n] = ((3-alfa)*(3-alfa)*M_PI/16.0)*z12v[n]*xrmn_n*sump1;  /*Z0311=9760*/
                CR->carr6p[n] = ((3-alfa)*(3-alfa)*M_PI/16.0)*z12v[n]*xrn[n-1]*sump;  /*Z0311=9761*/

                /* F(q)-coefficients */  /*Z0311=9763*/
                CR->carr4f[n] = (3*sqrt(M_PI)/4.0)*z12v[n]*xrn[n]/((n+3/2.0)*gam3[n]*fkv[n]);  /*Z0311=9764*/
                CR->carr5f[n] = (sqrt(M_PI)*(3-alfa)/4.0)*z12v[n]*xrmn_n/((n+3/2.0-alfa/2.0)*gam3[n]*fkv[n]);  /*Z0311=9765*/
                CR->carr6f[n] = (sqrt(M_PI)*(3-alfa)/4.0)*z12v[n]*xrn[n]/((n+3/2.0-alfa/2.0)*gam3[n]*fkv[n]);  /*Z0311=9766*/
                if ( fabs(CR->carr1p[n])<min )
                {   /*Z0311=9767*/
                    if ( n<n1 ) n1 = n;  /*Z0311=9768*/
                }   /*Z0311=9769*/
                if ( fabs(CR->carr2p[n])<min )
                {   /*Z0311=9770*/
                    if ( n<n2 ) n2 = n;  /*Z0311=9771*/
                }   /*Z0311=9772*/
                if ( fabs(CR->carr3p[n])<min )
                {   /*Z0311=9773*/
                    if ( n<n3 ) n3 = n;  /*Z0311=9774*/
                }   /*Z0311=9775*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z0311=9776*/
                    if ( n<n4 ) n4 = n;  /*Z0311=9777*/
                }   /*Z0311=9778*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z0311=9779*/
                    if ( n<n5 ) n5 = n;  /*Z0311=9780*/
                }   /*Z0311=9781*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z0311=9782*/
                    if ( n<n6 ) n6 = n;  /*Z0311=9783*/
                }   /*Z0311=9784*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z0311=9785*/
                    if ( n<n4f ) n4f = n;  /*Z0311=9786*/
                }   /*Z0311=9787*/
                if ( fabs(CR->carr5f[n])<min )
                {   /*Z0311=9788*/
                    if ( n<n5f ) n5f = n;  /*Z0311=9789*/
                }   /*Z0311=9790*/
                if ( fabs(CR->carr6f[n])<min )
                {   /*Z0311=9791*/
                    if ( n<n6f ) n6f = n;  /*Z0311=9792*/
                }   /*Z0311=9793*/
            }   /*Z0311=9794*/
            goto Label99;  /*Z0311=9795*/
        }   /* of inhomogeneous core/shell */  /*Z0311=9796*/

        /* myelin */  /*Z0311=9798*/
        if ( (cs==3) || (cs==4) )
        {   /*Z0311=9799*/
            i = 2;  /*Z0311=9800*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z0311=9801*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  /*Z0311=9802*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  /*Z0311=9803*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  /*Z0311=9804*/
                fkv[n] = fkv[n-1]*n;  /*Z0311=9805*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  /*Z0311=9806*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  /*Z0311=9807*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */  /*Z0311=9808*/
                xln[n] = -xln[n-1]*xl2z;  /*Z0311=9809*/
                /* xrn[n]:=-xrn[n-1]*xr2z; */  /*Z0311=9810*/
                xrn[n] = -xrn[n-1]*x12zm;         /* myelin radius */  /*Z0311=9811*/

                /* P(q) */  /*Z0311=9813*/
                for ( m=0; m<=n; m++ )
                {   /*Z0311=9814*/
                    /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]); */  /*Z0311=9815*/
                    /* i:=i+1; */  /*Z0311=9816*/
                    CR->carr11pm[n][m] = (9*M_PI/16.0)*(1/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]));  /*Z0311=9817, cs=3|4,dim=3*/
                }   /*Z0311=9818*/
                CR->carr4p[n] = z12v[n]*xrn[n];  /*Z0311=9819*/
                /* CR->carr4[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]); */  /*Z0311=9820*/

                /* F(q) */  /*Z0311=9823*/
                binsum = 0.0;  /*Z0311=9824*/
                for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  /*Z0311=9825*/
                CR->carr1f[n] = M_PI*xln[n]*binsum/(4.*(2*n+1));  /*Z0311=9826*/
                binsum = 0.0;  /*Z0311=9827*/
                for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  /*Z0311=9828*/
                CR->carr4f[n] = xrn[n-1]*binsum;  /*Z0311=9829*/

                if ( fabs(CR->carr1p[n])<min )
                {   /*Z0311=9832*/
                    if ( n<n1 ) n1 = n;  /*Z0311=9833*/
                }   /*Z0311=9834*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z0311=9835*/
                    if ( n<n1f ) n1f = n;  /*Z0311=9836*/
                }   /*Z0311=9837*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z0311=9838*/
                    if ( n<n4 ) n4 = n;  /*Z0311=9839*/
                }   /*Z0311=9840*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z0311=9841*/
                    if ( n<n4f ) n4f = n;  /*Z0311=9842*/
                }   /*Z0311=9843*/
            }   /*Z0311=9844*/
        }  /* of myelin */  /*Z0311=9845*/

    }   /* of dim=3, spheres */  /*Z0311=9849*/

    /* isotropic case for cubes */  /*Z=9900*/
    if ( (ordis==7) && (dim==4) )
    {    /* cubes */
        norm = 1;   /*Z=9902*/
        order = 0;   /*Z=9903*/
        /* homogeneous */   /*Z=9904*/
        if ( cs==0 )
        {   /*Z=9905*/
            area = 6*4*r*r;   /*Z=9906*/
            vol = 8*r*r*r;   /*Z=9907*/
            por = 2*M_PI*pow(z+1,4.)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);   /*Z=9908*/

            u1ell[0] = 2;   /*Z=9910*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=9911*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=9912*/
                fkv[n]  = fkv[n-1]*n;   /*Z=9913*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=9914*/
                /* u1ell[n]:=z12v[n]/((n+1/2)*(n+1)*fkv[n]); */   /*Z=9915*/
                u1ell[n] = 1/((n+1/2.0)*(n+1)*fkv[n]);   /*Z=9916*/
            }   /*Z=9917*/

            for ( n=0; n<=nmax; n++ )
            {   /*Z=9919*/
                sump1 = 0.0;   /*Z=9920*/
                for ( m=0; m<=n; m++ ) sump1 = sump1+u1ell[n-m]*u1ell[m];   /*Z=9921*/
                v1ell[n] = sump1;   /*Z=9922*/
            }   /*Z=9923*/

            for ( n=1; n<=nmax; n++ )
            {   /*Z=9925*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=9926*/
                sumf = 0.0;   /*Z=9927*/
                for ( m=0; m<=n; m++ )
                    sumf = sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=9928*/
                fsum[n] = sumf;   /*Z=9929*/
                /* P(q)-coefficient */   /*Z=9930*/
                sump = 0.0;   /*Z=9931*/
                for ( m=0; m<=n; m++ ) sump = sump+u1ell[n-m]*v1ell[m];   /*Z=9944*/

                /* carr4p[n]:=sqrt(pi)*power(4,n)*xrn[n]*sump/(16*gam3[n]); */   /*Z=9946*/
                CR->carr4p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12v[n]*xrn[n-1]*sump/(16*gam3[n]);   /*Z=9947*/

                /* F(q)-coefficient */   /*Z=9949*/
                sump = 0.0;   /*Z=9950*/
                for ( m=0; m<=n; m++ )
                {   /*Z=9951*/
                    sump1 = 0.0;   /*Z=9952*/
                    for ( k=0; k<=m; k++ )
                        sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));   /*Z=9953*/
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);   /*Z=9954*/
                }   /*Z=9955*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n-1]*sump/(128*gam3[n]);   /*Z=9956*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=9957*/
                    if ( n<n4 ) n4 = n;   /*Z=9958*/
                }   /*Z=9959*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=9960*/
                    if ( n<n4f ) n4f = n;   /*Z=9961*/
                }   /*Z=9962*/
            }   /*Z=9963*/
            goto Label99;   /*Z=9964*/
        }   /* of homogeneous */   /*Z=9965*/

        /* core/shell */          /* not yet ready */   /*Z=9967*/
        if ( cs==1 )
        {   /*Z=9968*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=9969*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=9970*/
                fkv[n] = fkv[n-1]*n;   /*Z=9971*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=9972*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=9973*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=9974*/
                xrmn_n = -xrmn_n*xrm2z;   /*Z=9975*/
                pn[n] = pn[n-1]*p*p;   /*Z=9976*/
                /* P(q)-coefficient */   /*Z=9977*/
                sump = 0.0;   /*Z=9978*/
                for ( m=0; m<=n; m++ )
                {   /*Z=9979*/
                    sump = sump+pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=9980*/
                }   /*Z=9981*/
                CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.,1.0*n)*z12v[n]*xrn[n]/(2.*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);   /*Z=9982*/
                CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;   /*Z=9983*/
                CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=9984*/
                /* F(q)-coefficient */   /*Z=9985*/
                /* carr3[n]:=3*sqrt(pi)*z12v[n]*xrn[n]/(4*(n+3/2)*gam3[n]*fkv[n]); */   /*Z=9986*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=9987*/
                    if ( n<n4 ) n4 = n;   /*Z=9988*/
                }   /*Z=9989*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z=9990*/
                    if ( n<n5 ) n5 = n;   /*Z=9991*/
                }   /*Z=9992*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z=9993*/
                    if ( n<n6 ) n6 = n;   /*Z=9994*/
                }   /*Z=9995*/
            }   /*Z=9996*/
            goto Label99;   /*Z=9997*/
        }   /* of core/shell */   /*Z=9998*/
    }  /* of dim=4, cubes */  /*Z=9999*/

    /* perfect orientation case for cubes */  /*Z=10001*/
    if ( (ordis==6) && (dim==4) )
    {   /* cubes */
        norm = 1;   /*Z=10003*/
        order = 1;   /*Z=10004*/
        /* homogeneous */   /*Z=10005*/
        if ( cs==0 )
        {   /*Z=10006*/
            i = 2;   /*Z=10007*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=10008*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10009*/
                fkv[n] = fkv[n-1]*n;   /*Z=10010*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10011*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=10012*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10013*/
                /* P(q)-coefficient */   /*Z=10014*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10015*/
                    /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]); */   /*Z=10016*/
                    /* carr1fm[i]:=carr1pm[i]; */   /*Z=10017*/
                    CR->carr11pm[n][m] = (M_PI/4.0)*pow(4.,1.0*n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);   /*Z=10018, cs=0,ordis=6,dim=4*/
                    i = i+1;   /*Z=10019*/
                }   /*Z=10020*/
                CR->carr4p[n] = sqrt(M_PI)*pow(4.,1.0*n)*xrn[n]/(16*gam3[n]);   /*Z=10021*/
                /* F(q)-coefficient */   /*Z=10022*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128*gam3[n]);   /*Z=10023*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10024*/
                    if ( n<n4 ) n4 = n;   /*Z=10025*/
                }   /*Z=10026*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10027*/
                    if ( n<n4f ) n4f = n;   /*Z=10028*/
                }   /*Z=10029*/
            }   /*Z=10030*/
            goto Label99;   /*Z=10031*/
        }   /* of homogeneous */   /*Z=10032*/

        /* core/shell */          /* not yet ready */   /*Z=10034*/
        if ( cs==1 )
        {   /*Z=10035*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=10036*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10037*/
                fkv[n] = fkv[n-1]*n;   /*Z=10038*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10039*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=10040*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10041*/
                xrmn_n = -xrmn_n*xrm2z;   /*Z=10042*/
                pn[n] = pn[n-1]*p*p;   /*Z=10043*/
                /* P(q)-coefficient */   /*Z=10044*/
                sump = 0.0;   /*Z=10045*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10046*/
                    sump = sump+pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=10047*/
                }   /*Z=10048*/
                CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.,1.0*n)*z12v[n]*xrn[n]/(2.*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);   /*Z=10049*/
                CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;   /*Z=10050*/
                CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=10051*/
                /* F(q)-coefficient */   /*Z=10052*/
                /* carr3[n]:=3*sqrt(pi)*z12v[n]*xrn[n]/(4*(n+3/2)*gam3[n]*fkv[n]); */   /*Z=10053*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10054*/
                    if ( n<n4 ) n4 = n;   /*Z=10055*/
                }   /*Z=10056*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z=10057*/
                    if ( n<n5 ) n5 = n;   /*Z=10058*/
                }   /*Z=10059*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z=10060*/
                    if ( n<n6 ) n6 = n;   /*Z=10061*/
                }   /*Z=10062*/
            }   /*Z=10063*/
            goto Label99;   /*Z=10064*/
        }   /* of core/shell */   /*Z=10065*/
    }   /* of dim=4, cubes */  /*Z=10066*/

    /* isotropic case for ellipsoids */  /*Z=10068*/
    if ( (ordis==7) && (dim==5) )
    {   /* ellipsoids */
        norm = 1;   /*Z=10070*/
        order = 0;   /*Z=10071*/
        if ( eps==1 ) area = 4*M_PI*r*r;   /*Z=10072*/
        if ( eps> 1 ) area = 2*M_PI*r*(r+(l*l/sqrt(l*l-r*r))*asin(sqrt(l*l-r*r)/l));   /*Z=10073*/
        if ( eps< 1 ) area = 2*M_PI*r*(r+(l*l/sqrt(r*r-l*l))*asinh(sqrt(r*r-l*l)/l));   /*Z=10074*/
        vol = (4*M_PI/3.0)*r*r*l;   /*Z=10075*/
        por = 2*M_PI*pow(z+1,4.)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);   /*Z=10076*/

        /* homogeneous */   /*Z=10078*/
        if ( cs==0 )
        {   /*Z=10079*/
            double e1[nnmax+1];
            e1[0] = 1;
            for ( n=1; n<=nmax; n++ )
            {   /*Z=10080*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10081*/
                fkv[n] = fkv[n-1]*n;   /*Z=10082*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10083*/
                e1[n] = e1[n-1]*(eps*eps-1);   /*Z=10084*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10085*/
                sump = 0.0;   /*Z=10086*/
                for ( m=0; m<=n; m++ ) sump = sump+e1[n-m]/(fkv[n-m]*fkv[m]*(2*(n-m)+1));   /*Z=10087*/
                /* sumf:=0.0; */   /*Z=10088*/
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]); */   /*Z=10089*/
                /* fsum[n]:=sumf; */   /*Z=10090*/
                /* P(q)-coefficient */   /*Z=10091*/
                CR->carr4p[n] = (9*sqrt(M_PI)/2.0)*pow(4.,1.0*n)*z12v[n]*xrn[n-1]*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);   /*Z=10092*/
                /* F(q)-coefficient */   /*Z=10093*/
                sump = 0.0;   /*Z=10094*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10095*/
                    sump1 = 0.0;   /*Z=10096*/
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));   /*Z=10097*/
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);   /*Z=10098*/
                }   /*Z=10099*/
                CR->carr4f[n]  = M_PI*M_PI*xrn[n-1]*sump/(128*gam3[n]);   /*Z=10100*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10101*/
                    if ( n<n4 ) n4 = n;   /*Z=10102*/
                }   /*Z=10103*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10104*/
                    if ( n<n4f ) n4f = n;   /*Z=10105*/
                }   /*Z=10106*/
            }   /*Z=10107*/
            goto Label99;   /*Z=10108*/
        }   /* of homogeneous */   /*Z=10109*/

        /* core/shell */          /* not yet ready */   /*Z=10111*/
        //if cs=1 then begin   /*Z=10112*/

    }   /* of dim=5, ellipsoids */  /*Z=10143*/

    /* perfect orientation case for ellipsoids */  /*Z=10145*/
    if ( (ordis==6) && (dim==5) )
    {   /* ellipsoids */
        norm = 1;   /*Z=10147*/
        order = 1;   /*Z=10148*/
        /* homogeneous */   /*Z=10149*/
        if ( cs==0 )
        {   /*Z=10150*/
            for ( n=1; n<=2*nmax+2; n++ )
            {   /*Z=10151*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10152*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=10153*/
                fkv[n] = fkv[n-1]*n;   /*Z=10154*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10155*/
                //e1[n] = e1[n-1]*(epsi*epsi-1);   /*Z=10156*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10157*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=10158*/
            }   /*Z=10159*/
            for ( n=0; n<=nmax; n++ )
            {   /*Z=10160*/
                a1 = sqr(3/4.0);   /*Z=10161*/
                for ( m=0; m<=nmax; m++ )
                {   /*Z=10162*/
                    sump1 = 0.0;   /*Z=10163*/
                    for ( k=0; k<=n; k++ )
                    {   /*Z=10164*/
                        sump2 = 0.0;   /*Z=10165*/
                        for ( ll=0; ll<=m; ll++ )   /*Z=10166*/
                            sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);   /*Z=10167*/
                        sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);   /*Z=10168*/
                    }   /*Z=10169*/
                    CR->carr11pm[n][m] = M_PI*sump1;   /*Z=10170, cs=0,ordis=6,dim=5*/
                }   /*Z=10171*/

                CR->carr4p[n] = a1*z12vl[n]*xln[n];   /*Z=10173*/
                CR->carr5p[n] = z12v[n]*xrn[n];   /*Z=10174*/

                /* P(q)-coefficient */   /*Z=10176*/
                /* for m:=0 to n do begin */   /*Z=10177*/
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]); */   /*Z=10178*/
                /* carr1fm[i]:=carr1pm[i]; */   /*Z=10179*/
                /* i:=i+1; */   /*Z=10180*/
                /* end; */   /*Z=10181*/

                /* F(q)-coefficient */   /*Z=10183*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128*gam3[n]);   /*Z=10184*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10185*/
                    if ( n<n4 ) n4 = n;   /*Z=10186*/
                }   /*Z=10187*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10188*/
                    if ( n<n4f ) n4f = n;   /*Z=10189*/
                }   /*Z=10190*/
            }   /*Z=10191*/
            goto Label99;   /*Z=10192*/
        }   /* of homogeneous */   /*Z=10193*/

        /* core/shell */          /* not yet ready */   /*Z=10195*/
        // if cs=1 then begin   /*Z=10196*/

    }  /* of dim=5, ellipsoid */  /*Z=10227*/

    /*** ellipsoid orientational distribution */  /*Z=10229*/
    if ( (ordis==0) && (cho1==2) && (dim==5) )
    {
        if ( cs==0 )
        {   /*Z=10232*/
            //qrombdeltac(l,r,/*params.p1,sigma,dbeta,*/theta,phi,1,1,1, 9,9,9, 9,9,9,9, 9,9,9, 9,9,9, 9,9,9, 9,9,9,
            //            ordis, 3, 2, 0, 0, 0, 0, CR->carr1, norm );  /*Z0311=10162*/
            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,2,0,0,0,0,CR->carr1p,norm);   /*Z=10233*/
            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,3,0,0,0,0,CR->carr1p,order);   /*Z=10234*/
            //TODO 30.07.2022: Parameter alpha ist hier neu, ich warte noch auf die passende Routine als Source
            order = order/norm;   /*Z=10235*/

            for ( n=1; n<=2*nmax+2; n++ )
            {   /*Z=10237*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10238*/
                z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10239*/
                fkv[n] = fkv[n-1]*n;   /*Z=10240*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=10241*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10242*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=10243*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10244*/
            }   /*Z=10245*/

            for ( n=0; n<=nmax; n++ )
            {   /*Z=10247*/
                for ( m=0; m<=nmax; m++ )
                {   /*Z=10248*/
                    sump1 = 0.0;   /*Z=10249*/
                    for ( k=0; k<=n; k++ )
                    {   /*Z=10250*/
                        sump2 = 0.0;   /*Z=10251*/
                        for ( ll=0; ll<=m; ll++ )   /*Z=10252*/
                            sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);   /*Z=10253*/
                        sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);   /*Z=10254*/
                    }   /*Z=10255*/
                    CR->carr11pm[n][m] = M_PI*sump1;   /*Z=10256, cs=0,ordis=0,cho1=2,dim=5*/
                }   /*Z=10257*/
                sump1 = 0.0;   /*Z=10258*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10259*/
                    qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=10260*/
                    // TOD: alpha ist hier neu (s.o.)
                    sump1 = sump1+pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=10261*/
                    CR->carr22pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=10262*/
                }   /*Z=10263*/

                /* all coefficients: Mok */   /*Z=10265*/
                /* integral part */   /*Z=10266*/
                CR->carr4p[n] = sump1;   /*Z=10267*/
                /* qR-part */   /*Z=10268*/
                CR->carr5p[n] = z12v[n]*xrn[n];   /*Z=10269*/
                /* qL-part */   /*Z=10270*/
                /* carr3p[n]:=sqr(3/4)*fk2v[n]*z12v[n]*xln[n]/power(4,n); */   /*Z=10271*/
                CR->carr6p[n] = sqr(3/4.0)*fk2v[n]*z12vl[n]*xln[n];   /*Z=10272*/

                if ( search5 )
                {   /*Z=10289*/
                    if ( fabs(CR->carr5p[n])<1e-50 )
                    {   /*Z=10290*/
                        n5 = n;   /*Z=10291*/
                        search1 = false;   /*Z=10292*/
                    }   /*Z=10293*/
                }   /*Z=10294*/
                if ( search6 )
                {   /*Z=10295*/
                    if ( fabs(CR->carr6p[n])<1e-50 )
                    {   /*Z=10296*/
                        n6 = n;   /*Z=10297*/
                        search6 = false;   /*Z=10298*/
                    }   /*Z=10299*/
                }   /*Z=10300*/
                if ( search3 )
                {   /*Z=10301*/
                    if ( fabs(CR->carr3p[n])<1e-50 )
                    {   /*Z=10302*/
                        n3 = n;   /*Z=10303*/
                        search3 = false;   /*Z=10304*/
                    }   /*Z=10305*/
                }   /*Z=10306*/
                if ( search4 )
                {   /*Z=10307*/
                    if ( fabs(CR->carr4p[n])<1e-50 )
                    {   /*Z=10308*/
                        n4 = n;   /*Z=10309*/
                        search4 = false;   /*Z=10310*/
                    }   /*Z=10311*/
                }   /*Z=10312*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=10313*/
                    if ( n<n1f ) n1f = n;   /*Z=10314*/
                }   /*Z=10315*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10316*/
                    if ( n<n4f ) n4f = n;   /*Z=10317*/
                }   /*Z=10318*/
            }  /* of n-loop */   /*Z=10319*/

            goto Label99;   /*Z=10321*/
        }  /* of cs-loop */   /*Z=10322*/
    }   /* if ( (ordis==0) && (cho1==2) && (dim==5) ) */  /*Z=10323*/

    /* isotropic case for triaxial ellipsoids */  /*Z=10326*/
    if ( (ordis==7) && (dim==6) )     /* triaxial ellipsoids */
    {
        norm = 1;   /*Z=10328*/
        order = 0;   /*Z=10329*/
        if ( (aell==bell) && (bell==cell) )
            area = 4*M_PI*aell*aell;   /*Z=10330*/
        else
        {   double ab = aell * bell;    //220908
            double bc = bell * cell;
            double ac = aell * cell;
            area = 4*M_PI*pow((pow(ab,8/5.0)+pow(bc,8/5.0)+pow(ac,8/5.0))/3.0,5/8.0);   /*Z=10331*/
        }
        vol = (4*M_PI/3.0)*aell*bell*cell;   /*Z=10332*/
        por = 2*M_PI*pow(z+1,4.)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);   /*Z=10333*/

        /* homogeneous */   /*Z=10335*/
        if ( cs==0 )
        {   /*Z=10336*/
            u1ell[0] = 1;   /*Z=10338*/
            u2ell[0] = 1;   /*Z=10339*/
            v1ell[0] = 1;   /*Z=10340*/
            v2ell[0] = 1;   /*Z=10341*/
            gell[0] = sqrt(M_PI);   /*Z=10342*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=10343*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10344*/
                fkv[n] = fkv[n-1]*n;   /*Z=10345*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10346*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10347*/
                u1ell[n] = (aell*aell-bell*bell)*u1ell[n-1]/double(n);   /*Z=10348*/
                u2ell[n] = bell*bell*u2ell[n-1]/n;   /*Z=10349*/
                v1ell[n] = (bell*bell-aell*aell)*v1ell[n-1]/double(n);   /*Z=10350*/
                v2ell[n] = (cell*cell-bell*bell)*v2ell[n-1]/double(n);   /*Z=10351*/
                gell[n] = gam3[n]/((n+1/2.0)*fkv[n]);   /*Z=10352*/
            }   /*Z=10353*/

            for ( n=1; n<=nmax; n++ )
            {   /*Z=10355*/
                sump = 0.0;   /*Z=10368*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10369*/
                    sump1 = 0.0;   /*Z=10370*/
                    for ( k=0; k<=m; k++ )
                    {   /*Z=10371*/
                        sump2 = 0.0;   /*Z=10372*/
                        for ( ll=0; ll<=n-m; ll++ ) sump2 = sump2+gell[k+ll]*v1ell[ll]*v2ell[n-m-ll];   /*Z=10373*/
                        sump1 = sump1+u1ell[k]*u2ell[m-k]*sump2;   /*Z=10374*/
                    }   /*Z=10375*/
                    sump = sump+sump1/((2.*(n-m))+1);   /*Z=10376*/
                }   /*Z=10377*/

                /* fsum[n]:=sumf; */   /*Z=10379*/
                sumf = 0.0;   /*Z=10380*/
                for ( m=0; m<=n; m++ ) sumf = sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=10381*/
                fsum[n] = sumf;   /*Z=10382*/
                /* P(q)-coefficient */   /*Z=10383*/
                CR->carr4p[n] = (9*M_PI)*pow(4.,n-1.)*z12v[n]*pow(-1./(4.*(z+1)*(z+1)),1.0*n)*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);   /*Z=10384*/
                /* F(q)-coefficient */   /*Z=10385*/
                sump = 0.0;   /*Z=10386*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10387*/
                    sump1 = 0.0;   /*Z=10388*/
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));   /*Z=10389*/
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);   /*Z=10390*/
                }   /*Z=10391*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n-1]*sump/(128*gam3[n]);   /*Z=10392*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10393*/
                    if ( n<n4 ) n4 = n;   /*Z=10394*/
                }   /*Z=10395*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10396*/
                    if ( n<n4f ) n4f = n;   /*Z=10397*/
                }   /*Z=10398*/
            }   /*Z=10399*/
            goto Label99;   /*Z=10400*/
        }   /* of homogeneous */   /*Z=10401*/

        /* core/shell */          /* not yet ready */   /*Z=10403*/
        //if cs=1 then begin   /*Z=10404*/

    }   /* of dim=6, ellipsoids */  /*Z=10435*/

    /* perfect orientation case for ellipsoids */  /*Z=10437*/
    if ( (ordis==6) && (dim==5) )     /* ellipsoids */
    {
        norm = 1;   /*Z=10439*/
        order = 1;   /*Z=10440*/
        /* homogeneous */   /*Z=10441*/
        if ( cs==0 )
        {   /*Z=10442*/
            i = 2;   /*Z=10443*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=10444*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10445*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=10446*/
                fkv[n] = fkv[n-1]*n;   /*Z=10447*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10448*/
                //e1[n] = e1[n-1]*(epsi*epsi-1);   /*Z=10449*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10450*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=10451*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10452*/
                    sump = 0.0;   /*Z=10453*/
                    for ( ns=0; ns<=n; ns++ )
                    {   /*Z=10454*/
                        sump1 = 0.0;   /*Z=10455*/
                        for ( ms=0; ms<=m; ms++ ) sump1 = sump1+1/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);   /*Z=10456*/
                        sump = sump+sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);   /*Z=10457*/
                    }   /*Z=10458*/
                    CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;   /*Z=10459, cs=0,ordis=6,dim=5*/
                }   /*Z=10460*/
                /* P(q)-coefficient */   /*Z=10461*/
                /* for m:=0 to n do begin */   /*Z=10462*/
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]); */   /*Z=10463*/
                /* carr1fm[i]:=carr1pm[i]; */   /*Z=10464*/
                /* i:=i+1; */   /*Z=10465*/
                /* end; */   /*Z=10466*/
                CR->carr4p[n] = sqrt(M_PI)*pow(4.,1.0*n)*xrn[n]/(16.*gam3[n]);   /*Z=10467*/
                /* F(q)-coefficient */   /*Z=10468*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128*gam3[n]);   /*Z=10469*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10470*/
                    if ( n<n4 ) n4 = n;   /*Z=10471*/
                }   /*Z=10472*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10473*/
                    if ( n<n4f ) n4f = n;   /*Z=10474*/
                }   /*Z=10475*/
            }   /*Z=10476*/
            goto Label99;   /*Z=10477*/
        }   /* of homogeneous */   /*Z=10478*/

        /* core/shell */          /* not yet ready */   /*Z=10480*/
        //if cs=1 then begin   /*Z=10481*/

    }  /* of dim=6, triellipsoid */  /*Z=10512*/

    /* isotropic case for super ellipsoids, barrel */  /*Z=10515*/
    if ( (ordis==7) && (dim==7) )    /* super ellipsoids */
    {
        // TODO
        double qx=1, qy=1, qz=1, qhkl=1;
        // Diese Variablen sind im Orginal-Pascalprogramm an der Stelle des coefficients() Aufrufes nicht bekannt.

        // TODO: die Variablen p11 bis p33 werden im Vorfeld zu diesem Aufruf nicht gesetzt
//        p11=p12=p13=p21=p22=p23=p31=p32=p33 = 1; // da später damit multipliziert wird macht eine 0 weniger Sinn

        norm = 1;   /*Z=10517*/
        order = 0;   /*Z=10518*/
        //    void qrombchid( double l, double r, double sigma, double dbeta, double delta, double theta, double phi,
        qrombchid(params.length,params.radius,/*alfa,*/sigma,/*alfa,*/dbeta,epsi,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qx,qy,0,qhkl,
                  ax1.length(),ax2.length(),ax3.length(),ax1.x(),ax1.y(),ax1.z(),ax2.x(),ax2.y(),ax2.z(),ax3.x(),ax3.y(),ax3.z(),
                  sig.x(),sig.y(),sig.z(),ordis,3,8,15,7,0,0,CR->carr1p,area);   /*Z=10519*/
        area = 2*M_PI*area;   /*Z=10520*/
        vol = 2*M_PI*r*r*l*gamma((2+alfa)/alfa)*gamma(1/alfa)/(alfa*gamma((3+alfa)/alfa));   /*Z=10521*/
        por = 2*M_PI*pow(z+1,4.)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);   /*Z=10522*/

        /* homogeneous */   /*Z=10524*/
        u1ell[0] = 1;             /* for barrel */   /*Z=10525*/
        u2ell[0] = 1;             /* for barrel */   /*Z=10526*/
        for ( n=1; n<=3*nmax+1; n++ )
        {   /*Z=10527*/
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10528*/
            fkv[n] = fkv[n-1]*n;   /*Z=10529*/
            gam3[n] = gam3[n-1]*(2*n+1.)/2.0;   /*Z=10530*/
            u1ell[n] = u1ell[n-1]/((n-1./2.0)*n);   /* for barrel */   /*Z=10531*/
            u2ell[n] = u2ell[n-1]/((n+1.)*n);     /* for barrel */   /*Z=10532*/
        }   /*Z=10533*/

        for ( n=0; n<=3*nmax+1; n++ )
        {   /*Z=10535*/
            v1ell[n] = gamma((2*n+1)/alfa)*u1ell[n];         /* for barrel */   /*Z=10536*/
            v2ell[n] = gamma((2*n+2+alfa)/alfa)*u2ell[n];    /* for barrel */   /*Z=10537*/
            gell[n] = gamma((2*n+3+alfa)/alfa);              /* for barrel */   /*Z=10538*/
        }   /*Z=10539*/

        if ( cs==0 )
        {   /*Z=10541*/
            for ( n=0; n<=nmax; n++ )
            {   /*Z=10542*/
                /* for k variable, very slow */   /*Z=10543*/
                if ( alfa==2 )
                {   /*Z=10556*/
                    a1 = sqr(3/4.0);   /*Z=10557*/
                    for ( m=0; m<=nmax; m++ )
                    {   /*Z=10558*/
                        sump1 = 0.0;   /*Z=10559*/
                        for ( k=0; k<=n; k++ )
                        {   /*Z=10560*/
                            sump2 = 0.0;   /*Z=10561*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=10562*/
                                sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);   /*Z=10563*/
                            sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);   /*Z=10564*/
                        }   /*Z=10565*/
                        CR->carr11pm[n][m] = M_PI*sump1;   /*Z=10566, cs=0,ordis=7,dim=7,alfa=2*/
                    }   /*Z=10567*/
                }   /*Z=10568*/
                else
                {   /*Z=10569*/
                    a1 = gamma((alfa+3)/alfa)/(gamma((alfa+2)/alfa)*gamma(1/alfa));   /*Z=10570*/
                    a1 = a1*a1;   /*Z=10571*/
                    for ( m=0; m<=nmax; m++ )
                    {   /*Z=10572*/
                        sump1 = 0.0;   /*Z=10573*/
                        for ( ns=0; ns<=n; ns++ )
                        {   /*Z=10574*/
                            sump2 = 0.0;   /*Z=10575*/
                            for ( ms=0; ms<=m; ms++ )   /*Z=10576*/
                                sump2 = sump2+v2ell[m-ms]*v2ell[ms]/(gell[ns+ms]*gell[n-ns+m-ms]);   /*Z=10577*/
                            sump1 = sump1+v1ell[n-ns]*v1ell[ns]*sump2;   /*Z=10578*/
                        }   /*Z=10579*/
                        CR->carr11pm[n][m] = sump1;   /*Z=10580, cs=0,ordis=7,dim=7,alfa!=2*/
                    }   /*Z=10581*/
                }   /*Z=10582*/

                /* orientational average */   /*Z=10584*/
                sump = 0.0;   /*Z=10585*/
                for ( m=0; m<=n; m++ ) sump = sump+gam3[n-m]*z12v[n-m]*z12v[m]*pow(l*l,1.0*n-m)*pow(r*r,1.0*m)*fkv[m]*CR->carr11pm[n-m][m]/(n-m+1/2.0);   /*Z=10586*/
                CR->carr4p[n] = a1*pow(-1./(4.*(z+1)*(z+1)),1.0*n)*sump/(2*gam3[n]);   /*Z=10587*/

                /* fsum[n]:=sumf; */   /*Z=10590*/
                /* sumf:=0.0; */   /*Z=10591*/
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]); */   /*Z=10592*/
                /* fsum[n]:=sumf; */   /*Z=10593*/
                /* P(q)-coefficient */   /*Z=10594*/
                /* F(q)-coefficient */   /*Z=10595*/
                sump = 0.0;   /*Z=10596*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10597*/
                    sump1 = 0.0;   /*Z=10598*/
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));   /*Z=10599*/
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);   /*Z=10600*/
                }   /*Z=10601*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n-1]*sump/(128*gam3[n]);   /*Z=10602*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10603*/
                    if ( n<n4 ) n4 = n;   /*Z=10604*/
                }   /*Z=10605*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10606*/
                    if ( n<n4f ) n4f = n;   /*Z=10607*/
                }   /*Z=10608*/
            }   /*Z=10609*/
            goto Label99;   /*Z=10610*/
        }   /* of homogeneous */   /*Z=10611*/

        /* core/shell */          /* not yet ready */   /*Z=10613*/
        //if cs=1 then begin   /*Z=10614*/

    }  /* of dim=5, ellipsoids */  /*Z=10645*/

    /* perfect orientation case for ellipsoids */  /*Z=10647*/
    if ( (ordis==6) && (dim==5) )    /* ellipsoids */
    {
        norm = 1;   /*Z=10649*/
        order = 1;   /*Z=10650*/
        /* homogeneous */   /*Z=10651*/
        if ( cs==0 )
        {   /*Z=10652*/
            i = 2;   /*Z=10653*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=10654*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10655*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=10656*/
                fkv[n] = fkv[n-1]*n;   /*Z=10657*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10658*/
                //e1[n] = e1[n-1]*(epsi*epsi-1);   /*Z=10659*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10660*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=10661*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10662*/
                    sump = 0.0;   /*Z=10663*/
                    for ( ns=0; ns<=n; ns++ )
                    {   /*Z=10664*/
                        sump1 = 0.0;   /*Z=10665*/
                        for ( ms=0; ms<=m; ms++ ) sump1 = sump1+1/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);   /*Z=10666*/
                        sump = sump+sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);   /*Z=10667*/
                    }   /*Z=10668*/
                    CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;   /*Z=10669, cs=0,ordis=6,dim=5*/
                }   /*Z=10670*/
                /* P(q)-coefficient */   /*Z=10671*/
                /* for m:=0 to n do begin */   /*Z=10672*/
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]); */   /*Z=10673*/
                /* carr1fm[i]:=carr1pm[i]; */   /*Z=10674*/
                /* i:=i+1; */   /*Z=10675*/
                /* end; */   /*Z=10676*/
                CR->carr4p[n] = sqrt(M_PI)*pow(4.,1.0*n)*xrn[n]/(16*gam3[n]);   /*Z=10677*/
                /* F(q)-coefficient */   /*Z=10678*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128*gam3[n]);   /*Z=10679*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10680*/
                    if ( n<n4 ) n4 = n;   /*Z=10681*/
                }   /*Z=10682*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10683*/
                    if ( n<n4f ) n4f = n;   /*Z=10684*/
                }   /*Z=10685*/
            }   /*Z=10686*/
            goto Label99;   /*Z=10687*/
        }   /* of homogeneous */   /*Z=10688*/

        /* core/shell */          /* not yet ready */   /*Z=10690*/
        //if cs=1 then begin   /*Z=10691*/

    }  /* of dim=7, super ellipsoid */  /*Z=10722*/

    /* isotropic case for superballs */   /*Z=10726*/
    if ( (ordis==7) && (dim==8) ) // NEU Aug.2022 - TODO: gammlog() unbekannt
    {    /* superball */   /*Z=10727*/
        //type ArrayImax3D=array[0..50,0..50,0..50] of real;
        //carr111pm: ^ArrayImax3D;
        //new(carr111pm);   /*Z=10728*/
        float carr111pm[51][51][51];
        nmax = 40;   /*Z=10729*/
        norm = 1;   /*Z=10730*/
        order = 0;   /*Z=10731*/
        /* l:=r; */   /*Z=10732*/
        /*Z=10733*/
        /* radius=a, rm=b, length=c */   /*Z=10734*/
        qrombdeltac(l,r,/*rm,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,8,17,0,0,0,CR->carr1p,area);   /*Z=10735*/
        area = 8*area;   /*Z=10736*/
        vol = 8*r*rm*l*pow(gamma(1/alfa),3.)/(pow(alfa,3.)*gamma(1+3/alfa));   /*Z=10737*/
        por = 2*M_PI*pow(z+1,4.)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);   /*Z=10738*/

        /* homogeneous */   /*Z=10741*/
        u3ell[0] = 1;             /* for superball */   /*Z=10742*/
        for ( n=1; n<=3*nmax+1; n++ )
        {   /*Z=10743*/
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10744*/
            fkv[n] = fkv[n-1]*n;   /*Z=10745*/
            gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10746*/
            u3ell[n] = u3ell[n-1]/((n-1/2.0)*n);   /* for superball */   /*Z=10747*/
        }   /*Z=10748*/

        /* for n:=0 to 3*nmax+1 do begin */   /*Z=10750*/
        /*    v3ell[n]:=gamma((2*n+1)/alfa)*u3ell[n];  */       /* for superball */    /*Z=10751*/
        /*    g3ell[n]:=gamma(((2*n+3)/alfa)+1);       */       /* for superball */    /*Z=10752*/
        /* end; */   /*Z=10753*/

        for ( n=0; n<=3*nmax+1; n++ )
        {   /*Z=10755*/
#define gammlog(x) (log10(gamma(x)))
            v3ell[n] = exp(gammlog((2*n+1)/alfa))*u3ell[n];         /* for superball */   /*Z=10756*/
            g3ell[n] = exp(gammlog(((2*n+3)/alfa)+1));              /* for superball */   /*Z=10757*/
        }   /*Z=10758*/

        if ( cs==0 )
        {   /*Z=10760*/
            for ( n=0; n<=nmax; n++ )
            {   /*Z=10761*/
                if ( alfa==200000 )
                {   /*Z=10762*/
                    a1 = sqr(3/4.0);   /*Z=10763*/
                    for ( m=0; m<=nmax; m++ )
                    {   /*Z=10764*/
                        sump1 = 0.0;   /*Z=10765*/
                        for ( k=0; k<=n; k++ )
                        {   /*Z=10766*/
                            sump2 = 0.0;   /*Z=10767*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=10768*/
                                sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);   /*Z=10769*/
                            sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);   /*Z=10770*/
                        }   /*Z=10771*/
                        CR->carr11pm[n][m] = M_PI*sump1;   /*Z=10772, cs=0,ordis=7,dim=8,alfa=200000*/
                    }   /*Z=10773*/
                }   /*Z=10774*/
                else
                {   /*Z=10775*/
                    a1 = gamma(1+3/alfa)/pow(gamma(1/alfa),3.);   /*Z=10776*/
                    a1 = a1*a1;   /*Z=10777*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=10778*/
                        for ( k=0; k<=nmax; k++ )
                        {   /*Z=10779*/
                            sump1 = 0.0;   /*Z=10780*/
                            for ( ns=0; ns<=n; ns++ )
                            {   /*Z=10781*/
                                sump2 = 0.0;   /*Z=10782*/
                                for ( ms=0; ms<=m; ms++ )
                                {   /*Z=10783*/
                                    double sump3 = 0.0;   /*Z=10784*/
                                    for ( int ks=0; ks<=k; ks++ )   /*Z=10785*/
                                        sump3 = sump3+v3ell[k-ks]*v3ell[ks]/(g3ell[ns+ms+ks]*g3ell[n-ns+m-ms+k-ks]);   /*Z=10786*/
                                    sump2 = sump2+v3ell[m-ms]*v3ell[ms]*sump3;   /*Z=10787*/
                                }   /*Z=10788*/
                                sump1 = sump1+v3ell[n-ns]*v3ell[ns]*sump2;   /*Z=10789*/
                            }   /*Z=10790*/
                            carr111pm[n][m][k] = sump1;   /*Z=10791*/
                            carr111pm[m][n][k] = sump1;   /*Z=10792*/
                        }   /*Z=10793*/
                    }   /*Z=10794*/
                }   /*Z=10795*/

                /* orientational average for superball */   /*Z=10797*/
                sump = 0.0;   /*Z=10798*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10799*/
                    sump1 = 0.0;   /*Z=10800*/
                    for ( k=0; k<=m; k++ )   /*Z=10801*/
                        sump1 = sump1+z12v[m-k]*z12v[k]*gam3[m-k]*gam3[k]*pow(rm*rm,1.0*m-k)*pow(l*l,1.0*k)*carr111pm[n-m][m-k][k]/((m-k+1/2.0)*(k+1/2.0));     /*  */   /*Z=10802*/
                    sump = sump+z12v[n-m]*gam3[n-m]*pow(r*r,1.0*n-m)*sump1/(n-m+1/2.0);   /*Z=10803*/
                }   /*Z=10804*/
                CR->carr4p[n] = (a1/(2*M_PI))*pow(-1./(4.*(z+1)*(z+1)),1.0*n)*sump/gam3[n];   /*Z=10805*/

                /* fsum[n]:=sumf; */   /*Z=10808*/
                /* sumf:=0.0; */   /*Z=10809*/
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]); */   /*Z=10810*/
                /* fsum[n]:=sumf; */   /*Z=10811*/
                /* P(q)-coefficient */   /*Z=10812*/
                /* F(q)-coefficient */   /*Z=10813*/
                sump = 0.0;   /*Z=10814*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10815*/
                    sump1 = 0.0;   /*Z=10816*/
                    for ( k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));   /*Z=10817*/
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);   /*Z=10818*/
                }   /*Z=10819*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n-1]*sump/(128*gam3[n]);   /*Z=10820*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10821*/
                    if ( n<n4 ) n4 = n;   /*Z=10822*/
                }   /*Z=10823*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=10824*/
                    if ( n<n4f ) n4f = n;   /*Z=10825*/
                }   /*Z=10826*/
            }   /*Z=10827*/
            goto Label99;   /*Z=10828*/
        }   /* of homogeneous */   /*Z=10829*/

        /* core/shell */          /* not yet ready */   /*Z=10831*/
        //if cs=1 then begin   /*Z=10832*/

        //dispose(carr111pm);   /*Z=10863*/
    }  /* of isotropic case for superballs */   /*Z=10864*/

    /* perfect orientation case for ellipsoids */   /*Z=10866*/
    if ( (ordis==6) && (dim==5) )  // NEU Aug.2022
    {    /* ellipsoids */   /*Z=10867*/
        norm = 1;   /*Z=10868*/
        order = 1;   /*Z=10869*/
        /* homogeneous */   /*Z=10870*/
        if ( cs==0 )
        {   /*Z=10871*/
            i = 2;   /*Z=10872*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=10873*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10874*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=10875*/
                fkv[n] = fkv[n-1]*n;   /*Z=10876*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10877*/
                //e1[n] = e1[n-1]*(epsi*epsi-1);   /*Z=10878*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=10879*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=10880*/
                for ( m=0; m<=n; m++ )
                {   /*Z=10881*/
                    sump = 0.0;   /*Z=10882*/
                    for ( ns=0; ns<=n; ns++ )
                    {   /*Z=10883*/
                        sump1 = 0.0;   /*Z=10884*/
                        for ( ms=0; ms<=m; ms++ ) sump1 = sump1+1/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);   /*Z=10885*/
                        sump = sump+sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);   /*Z=10886*/
                    }   /*Z=10887*/
                    CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;   /*Z=10888, cs=0,ordis=6,dim=5*/
                }   /*Z=10889*/
                /* P(q)-coefficient */   /*Z=10890*/
                /* for m:=0 to n do begin */   /*Z=10891*/
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]); */   /*Z=10892*/
                /* carr1fm[i]:=carr1pm[i]; */   /*Z=10893*/
                /* i:=i+1; */   /*Z=10894*/
                /* end; */   /*Z=10895*/
                CR->carr4p[n] = sqrt(M_PI)*pow(4.,1.0*n)*xrn[n]/(16*gam3[n]);   /*Z=10896*/
                /* F(q)-coefficient */   /*Z=10897*/
                CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128*gam3[n]);   /*Z=10898*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=10899*/
                    if ( n<n4 ) n4 = n;   /*Z=10900*/
                }   /*Z=10901*/
                if ( abs(CR->carr4f[n])<min )
                {   /*Z=10902*/
                    if ( n<n4f ) n4f = n;   /*Z=10903*/
                }   /*Z=10904*/
            }   /*Z=10905*/
            goto Label99;   /*Z=10906*/
        }   /* of homogeneous */   /*Z=10907*/

        /* core/shell */          /* not yet ready */   /*Z=10909*/
        //if cs=1 then begin   /*Z=10910*/

    }  /* of dim=8, superball */   /*Z=10941*/

    /*** isotropic case for cylinders and disks ***/  /*Z=10945*/
    if ( (ordis==7) && (dim!=3) )
    {
        norm = 1;   /*Z=10947*/
        order = 0;   /*Z=10948*/
        /* homogeneous */   /*Z=10949*/
        if ( cs==0 )
        {   /*Z=10950*/

            /* small axial ratios */   /*Z=10952*/
            if ( (l/r)<2 )
            {   /*Z=10953*/
                area = 2*M_PI*r*r+2*M_PI*r*(2*l);   /*Z=10954*/
                vol = M_PI*r*r*(2*l);   /*Z=10955*/
                por = 2*M_PI*pow(z+1,4.)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);   /*Z=10956*/

                for ( n=1; n<=nmax; n++ )
                {   /*Z=10958*/
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10959*/
                    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=10960*/
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=10961*/
                    fkv[n] = fkv[n-1]*n;   /*Z=10962*/
                    fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=10963*/
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10964*/
                    /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=10965*/
                    xln[n] = -xln[n-1]*xl2z;   /*Z=10966*/
                    xrn[n] = -xrn[n-1]*xr2z;   /*Z=10967*/
                }   /*Z=10968*/

                for ( n=0; n<=nmax; n++ )
                {   /*Z=10970*/
                    binsum = 0.0;   /*Z=10971*/
                    /* for m:=0 to n do */   /* Cauchy sum */    /*Z=10972*/
                    /*       binsum:=binsum+gam3[m]*z12vl[n-m]*z12v[m]*xln[n-m]*xrn[m]/((n-m+1/2)*(n-m+1)*fkv[n-m]*(m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]); */   /*Z=10973*/
                    CR->carr11pm[n][0] = sqrt(M_PI)/gam3[n];   /*Z=10974, cs=0,ordis=7,dim!=3*/
                    a1 = sqrt(M_PI)/(2*gam3[n]);   /*Z=10975*/
                    for ( m=1; m<=nmax; m++ )
                    {    /* double sum */   /*Z=10976*/
                        a1 = a1*(m+1/2.0)/(n+m+1/2.0);   /*Z=10977*/
                        /* carr11pm[n,m]:=power(4,m+1)*gam3[m]*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]*gam3[n+m]); */  /* ok */   /*Z=10978*/
                        CR->carr11pm[n][m] = pow(4.,m+1.)*a1*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]);       /* Mok */   /*Z=10979, cs=0,ordis=7,dim!=3*/
                    }   /*Z=10980*/
                    CR->carr2p[n] = pow(4.,n-1.)*z12vl[n]*xln[n]/((n+1/2.0)*(n+1)*fkv[n]);     /* double sum */   /*Z=10981*/
                    CR->carr3p[n] = pow(4.,1.0*n)*binsum/gam3[n];      /* Cauchy sum */   /*Z=10982*/
                }   /*Z=10983*/

            }   /*Z=10985*/
            /* large axial ratios */   /*Z=10986*/
            else
            {   /*Z=10987*/
                for ( n=1; n<=nmax; n++ )
                {   /*Z=10989*/
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=10990*/
                    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=10991*/
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=10992*/
                    fkv[n] = fkv[n-1]*n;   /*Z=10993*/
                    fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=10994*/
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=10995*/
                    /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=10996*/
                    xln[n] = -xln[n-1]*xl2z;   /*Z=10997*/
                    xrn[n] = -xrn[n-1]*xr2z;   /*Z=10998*/
                    /* cylinder, ok*/   /*Z=10999*/
                    if ( dim==1 )
                    {   /*Z=11000*/
                        /* P(q) factorization */   /*Z=11002*/
                        CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(2*n+1)*(n+1)*gam3[n]*fkv[n]);      /* P||iso(q) */   /*Z=11003*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /* P-(q) */   /*Z=11004*/
                        /* F(q) */   /*Z=11005*/
                        binsum = 0.0;   /*Z=11006*/
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11007*/
                        CR->carr1f[n] = M_PI*xln[n]*binsum/(4.*(2*n+1));   /*Z=11008*/
                        binsum = 0.0;   /*Z=11009*/
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11010*/
                        CR->carr4f[n] = xrn[n-1]*binsum;   /*Z=11011*/
                    }   /*Z=11012*/
                    /* disk, ok */   /*Z=11013*/
                    if ( dim==2 )
                    {   /*Z=11014*/
                        /* P(q) */   /*Z=11015*/
                        CR->carr1p[n] = 2*pow(4.,1.0*n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=11016*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=11017*/
                        /* F(q) */   /*Z=11018*/
                        binsum = 0.0;   /*Z=11019*/
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11020*/
                        CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2*gam3[n]);   /*Z=11021*/
                        binsum = 0.0;   /*Z=11022*/
                        for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11023*/
                        CR->carr4f[n] = M_PI*xrn[n-1]*binsum/4.0;   /*Z=11024*/
                    }   /*Z=11025*/
                }  /* of large axial ratios */   /*Z=11026*/

                if ( fabs(CR->carr1p[n])<min )
                {   /*Z=11028*/
                    if ( n<n1 ) n1 = n;   /*Z=11029*/
                }   /*Z=11030*/
                if ( fabs(CR->carr2p[n])<min )
                {   /*Z=11031*/
                    if ( n<n2 ) n2 = n;   /*Z=11032*/
                }   /*Z=11033*/
                if ( fabs(CR->carr3p[n])<min )
                {   /*Z=11034*/
                    if ( n<n3 ) n3 = n;   /*Z=11035*/
                }   /*Z=11036*/
                if ( fabs(CR->carr1f[n])<min )
                    {   /*Z=11037*/
                    if ( n<n1f ) n1f = n;   /*Z=11038*/
                }   /*Z=11039*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=11040*/
                    if ( n<n4 ) n4 = n;   /*Z=11041*/
                }   /*Z=11042*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=11043*/
                    if ( n<n4f ) n4f = n;   /*Z=11044*/
                }   /*Z=11045*/
            }   /*Z=11046*/
        } /* of homogeneous */   /*Z=11047*/

        /* core/shell */   /*Z=11049*/
        if ( cs==1 )
        {   /*Z=11050*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=11051*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11052*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11053*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11054*/
                fkv[n] = fkv[n-1]*n;   /*Z=11055*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11056*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11057*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11058*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=11059*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=11060*/
                xrmn_n = -xrmn_n*xrm2z;   /*Z=11061*/
                pn[n] = pn[n-1]*p*p;   /*Z=11062*/
                /*** cylinder ***/   /*Z=11063*/
                if ( dim==1 )
                {   /*Z=11064*/
                    /* longitudinal P(q) */   /*Z=11065*/
                    CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(2*n+1)*(n+1)*gam3[n]*fkv[n]);   /*Z=11066*/
                    /* cross-sectional P(q) */   /*Z=11067*/
                    /* F121 */   /*Z=11068*/
                    CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=11069*/
                    /* F122 */   /*Z=11070*/
                    sump = 0.0;   /*Z=11071*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11072*/
                        sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11073*/
                    }   /*Z=11074*/
                    CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=11075*/
                    /* F123 */   /*Z=11076*/
                    CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=11077*/

                    /* longitudinal F(q) */   /*Z=11079*/
                    binsum = 0.0;   /*Z=11080*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11081*/
                    CR->carr1f[n] = M_PI*xln[n]*binsum/(4.*(2*n+1));   /*Z=11082*/
                    /* cross-sectional F(q) */   /*Z=11083*/
                    /* F121 */   /*Z=11084*/
                    sump = 0.0;   /*Z=11085*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11086*/
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11087*/
                    }   /*Z=11088*/
                    CR->carr4f[n] = xrn[n-1]*sump;   /*Z=11089*/
                    /* F122 */   /*Z=11090*/
                    sump = 0.0;   /*Z=11091*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11092*/
                        sump = sump+z12v[m]*z12v[n-m]*pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11093*/
                    }   /*Z=11094*/
                    CR->carr5f[n] = xrmn_n*sump;   /*Z=11095*/
                    /* F123 */   /*Z=11096*/
                    CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=11097*/
                }   /*Z=11098*/

                /*** disk ***/   /*Z=11100*/
                if ( dim==2 )
                {   /*Z=11101*/
                    /* longitudinal */   /*Z=11102*/
                    CR->carr1p[n] = 2*pow(4.,1.0*n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=11103*/
                    /* cross-sectional P(q) */   /*Z=11104*/
                    /* F121 */   /*Z=11105*/
                    CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=11106*/
                    /* F122 */   /*Z=11107*/
                    sump = 0.0;   /*Z=11108*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11109*/
                        sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11110*/
                    }   /*Z=11111*/
                    CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;   /*Z=11112*/
                    /* F123 */   /*Z=11113*/
                    CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=11114*/

                    /* longitudinal F(q) */   /*Z=11116*/
                    binsum = 0.0;   /*Z=11117*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11118*/
                    CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2*gam3[n]);   /*Z=11119*/
                    /* cross-sectional F(q) */   /*Z=11120*/
                    /* F121 */   /*Z=11121*/
                    sump = 0.0;   /*Z=11122*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11123*/
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11124*/
                    }   /*Z=11125*/
                    CR->carr4f[n] = M_PI*xrn[n-1]*sump/4.0;   /*Z=11126*/
                    /* F122 */   /*Z=11127*/
                    sump = 0.0;   /*Z=11128*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11129*/
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11130*/
                    }   /*Z=11131*/
                    CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;   /*Z=11132*/
                    /* F123 */   /*Z=11133*/
                    CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=11134*/
                }   /*Z=11135*/
                if ( fabs(CR->carr1p[n])<min )
                {   /*Z=11136*/
                    if ( n<n1 ) n1 = n;   /*Z=11137*/
                }   /*Z=11138*/
                if ( fabs(CR->carr3p[n])<min )
                {   /*Z=11139*/
                    if ( n<n3 ) n3 = n;   /*Z=11140*/
                }   /*Z=11141*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=11142*/
                    if ( n<n4 ) n4 = n;   /*Z=11143*/
                }   /*Z=11144*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z=11145*/
                    if ( n<n5 ) n5 = n;   /*Z=11146*/
                }   /*Z=11147*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z=11148*/
                    if ( n<n6 ) n6 = n;   /*Z=11149*/
                }   /*Z=11150*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=11151*/
                    if ( n<n1f ) n1f = n;   /*Z=11152*/
                }   /*Z=11153*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=11154*/
                    if ( n<n4f ) n4f = n;   /*Z=11155*/
                }   /*Z=11156*/
                if ( fabs(CR->carr5f[n])<min )
                {   /*Z=11157*/
                    if ( n<n5f ) n5f = n;   /*Z=11158*/
                }   /*Z=11159*/
                if ( fabs(CR->carr6f[n])<min )
                {   /*Z=11160*/
                    if ( n<n6f ) n6f = n;   /*Z=11161*/
                }   /*Z=11162*/
            }   /*Z=11163*/
        } /* of core/shell */   /*Z=11164*/

        /* inhomogeneous core/shell */   /*Z=11166*/
        if ( cs==2 )
        {   /*Z=11167*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=11168*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11169*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11170*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11171*/
                fkv[n] = fkv[n-1]*n;   /*Z=11172*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11173*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11174*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11175*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=11176*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=11177*/
                xrmn_n = -xrmn_n*xrm2z;   /*Z=11178*/
                pn[n] = pn[n-1]*p*p;   /*Z=11179*/
                /*** cylinder ***/   /*Z=11180*/
                if ( dim==1 )
                {   /*Z=11181*/
                    /* longitudinal P(q) */   /*Z=11182*/
                    CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(2*n+1)*(n+1)*gam3[n]*fkv[n]);   /*Z=11183*/

                    /* cross-sectional P(q) */   /*Z=11185*/
                    CR->carr4p[n] = pow(4.,n+1.)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=11186*/
                    sump = 0.0;   /*Z=11187*/
                    sump1 = 0.0;   /*Z=11188*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11189*/
                        sumi = 1/((m+1-alfa/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);   /*Z=11190*/
                        sump = sump+pn[n-m]*sumi;   /*Z=11191*/
                        sump1 = sump1+sumi;   /*Z=11192*/
                    }   /*Z=11193*/
                    CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;   /*Z=11194*/
                    CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n-1]*sump1;   /*Z=11195*/
                    sump = 0.0;   /*Z=11196*/
                    sump1 = 0.0;   /*Z=11197*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11198*/
                        sumi = 1/((n-m+1-alfa/2.0)*(m+1-alfa/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);   /*Z=11199*/
                        sump = sump+sumi;   /*Z=11200*/
                        sump1 = sump1+pn[n-m]*sumi;   /*Z=11201*/
                    }   /*Z=11202*/
                    CR->carr7p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump;   /*Z=11203*/
                    CR->carr8p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump1;   /*Z=11204*/
                    CR->carr9p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrn[n-1]*sump;   /*Z=11205*/

                    /* longitudinal F(q) */   /*Z=11220*/
                    binsum = 0.0;   /*Z=11221*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11222*/
                    CR->carr1f[n] = M_PI*xln[n]*binsum/(4.*(2*n+1));   /*Z=11223*/
                    /* cross-sectional F(q) */   /*Z=11224*/
                    CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);   /*Z=11225*/
                    CR->carr5f[n] = (1-alfa/2.0)*z12v[n]*xrmn_n/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=11226*/
                    CR->carr6f[n] = (1-alfa/2.0)*z12v[n]*xrn[n]/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=11227*/
                }   /*Z=11229*/

                /*** disk ***/   /*Z=11231*/
                if ( dim==2 )
                {   /*Z=11232*/
                    /* longitudinal */   /*Z=11233*/
                    CR->carr1p[n] = 2*pow(4.,1.0*n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=11234*/

                    /* cross-sectional P(q) */   /*Z=11237*/
                    CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.,1.0*n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);   /*Z=11238*/
                    sump = 0.0;   /*Z=11239*/
                    sump1 = 0.0;   /*Z=11240*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11241*/
                        sumi = (m+1/2.0)/((m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=11242*/
                        sump = sump+pn[n-m]*sumi;   /*Z=11243*/
                        sump1 = sump1+sumi;   /*Z=11244*/
                    }   /*Z=11245*/
                    CR->carr5p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=11246*/
                    CR->carr6p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrn[n-1]*sump1;   /*Z=11247*/
                    sump = 0.0;   /*Z=11248*/
                    sump1 = 0.0;   /*Z=11249*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11250*/
                        sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-alfa/2.0)*(m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);   /*Z=11251*/
                        sump = sump+sumi;   /*Z=11252*/
                        sump1 = sump1+pn[n-m]*sumi;   /*Z=11253*/
                    }   /*Z=11254*/
                    CR->carr7p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=11255*/
                    CR->carr8p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump1;   /*Z=11256*/
                    CR->carr9p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrn[n-1]*sump;   /*Z=11257*/

                    /* cross-sectional P(q) */   /*Z=11260*/

                    /* longitudinal F(q) */   /*Z=11272*/
                    binsum = 0.0;   /*Z=11273*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11274*/
                    CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2*gam3[n]);   /*Z=11275*/
                    /* cross-sectional F(q) */   /*Z=11276*/
                    CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);   /*Z=11277*/
                    CR->carr5f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=11278*/
                    CR->carr6f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrn[n-1]*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=11279*/
                }   /*Z=11280*/
                if ( fabs(CR->carr1p[n])<min )
                {   /*Z=11281*/
                    if ( n<n1 ) n1 = n;   /*Z=11282*/
                }   /*Z=11283*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=11284*/
                    if ( n<n4 ) n4 = n;   /*Z=11285*/
                }   /*Z=11286*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z=11287*/
                    if ( n<n5 ) n5 = n;   /*Z=11288*/
                }   /*Z=11289*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z=11290*/
                    if ( n<n6 ) n6 = n;   /*Z=11291*/
                }   /*Z=11292*/
                if ( fabs(CR->carr7p[n])<min )
                {   /*Z=11293*/
                    if ( n<n7 ) n7 = n;   /*Z=11294*/
                }   /*Z=11295*/
                if ( fabs(CR->carr8p[n])<min )
                {   /*Z=11296*/
                    if ( n<n8 ) n8 = n;   /*Z=11297*/
                }   /*Z=11298*/
                if ( fabs(CR->carr9p[n])<min )
                {   /*Z=11299*/
                    if ( n<n9 ) n9 = n;   /*Z=11300*/
                }   /*Z=11301*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=11302*/
                    if ( n<n1f ) n1f = n;   /*Z=11303*/
                }   /*Z=11304*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=11305*/
                    if ( n<n4f ) n4f = n;   /*Z=11306*/
                }   /*Z=11307*/
                if ( fabs(CR->carr5f[n])<min )
                {   /*Z=11308*/
                    if ( n<n5f ) n5f = n;   /*Z=11309*/
                }   /*Z=11310*/
                if ( fabs(CR->carr6f[n])<min )
                {   /*Z=11311*/
                    if ( n<n6f ) n6f = n;   /*Z=11312*/
                }   /*Z=11313*/
                if ( fabs(CR->carr7f[n])<min )
                {   /*Z=11314*/
                    if ( n<n7f ) n7f = n;   /*Z=11315*/
                }   /*Z=11316*/
                if ( fabs(CR->carr8f[n])<min )
                {   /*Z=11317*/
                    if ( n<n8f ) n8f = n;   /*Z=11318*/
                }   /*Z=11319*/
                if ( fabs(CR->carr9f[n])<min )
                {   /*Z=11320*/
                    if ( n<n9f ) n9f = n;   /*Z=11321*/
                }   /*Z=11322*/
            }   /*Z=11323*/
        } /* of inhomogeneous core/shell */   /*Z=11324*/

        /* myelin */   /*Z=11327*/
        if ( (cs==3) || (cs==4) )
        {   /*Z=11328*/
            i = 2;   /*Z=11329*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=11330*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11331*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11332*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11333*/
                fkv[n] = fkv[n-1]*n;   /*Z=11334*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11335*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11336*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11337*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=11338*/
                /* xrn[n]:=-xrn[n-1]*xr2z; */   /*Z=11339*/
                xrn[n] = -xrn[n-1]*x12zm;         /* myelin radius */   /*Z=11340*/
                /* cylinder, ok*/   /*Z=11341*/
                if ( dim==1 )
                {   /*Z=11342*/
                    /* P(q) */   /*Z=11343*/
                    CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(2*n+1)*(n+1)*gam3[n]*fkv[n]);   /*Z=11344*/
                    /*Z=11345*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11346*/
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]); */   /*Z=11347*/
                        /* i:=i+1; */   /*Z=11348*/
                        CR->carr11pm[n][m] = 1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11349, dim=1,cs=3|4,ordis=7*/
                    }   /*Z=11350*/
                    CR->carr4p[n] = z12v[n]*xrn[n];   /*Z=11351*/
                    /* carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]); */   /*Z=11352*/

                    /* F(q) */   /*Z=11355*/
                    binsum = 0.0;   /*Z=11356*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11357*/
                    CR->carr1f[n] = M_PI*xln[n]*binsum/(4*(2*n+1));   /*Z=11358*/
                    binsum = 0.0;   /*Z=11359*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11360*/
                    CR->carr4f[n] = xrn[n-1]*binsum;   /*Z=11361*/
                }   /*Z=11362*/
                /* disk, ok */   /*Z=11363*/
                if ( dim==2 )
                {   /*Z=11364*/
                    /* P(q) */   /*Z=11365*/
                    CR->carr1p[n] = 2*pow(4.,1.0*n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=11366*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11367*/
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]); */   /*Z=11368*/
                        /* i:=i+1; */   /*Z=11369*/
                        CR->carr11pm[n][m] = (M_PI/4.0)*(1/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]));   /*Z=11370, dim=2,cs=3|4,ordis=7*/
                    }   /*Z=11371*/
                    CR->carr4p[n] = z12v[n]*xrn[n];   /*Z=11372*/

                    /* carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]); */   /*Z=11374*/
                    /* F(q) */   /*Z=11375*/
                    binsum = 0.0;   /*Z=11376*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11377*/
                    CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2*gam3[n]);   /*Z=11378*/
                    binsum = 0.0;   /*Z=11379*/
                    for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11380*/
                    CR->carr4f[n] = M_PI*xrn[n-1]*binsum/4.0;   /*Z=11381*/
                }   /*Z=11382*/
                if ( fabs(CR->carr1p[n])<min )
                {   /*Z=11383*/
                    if ( n<n1 ) n1 = n;   /*Z=11384*/
                }   /*Z=11385*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=11386*/
                    if ( n<n1f ) n1f = n;   /*Z=11387*/
                }   /*Z=11388*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=11389*/
                    if ( n<n4 ) n4 = n;   /*Z=11390*/
                }   /*Z=11391*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=11392*/
                    if ( n<n4f ) n4f = n;   /*Z=11393*/
                }   /*Z=11394*/
            }   /*Z=11395*/
        } /* of myelin */   /*Z=11396*/
    }   /*Z=11398*/

    /*** perfect orientation for cylinders and disks ***/  /*Z=11400*/
    if ( (ordis==6) && (dim!=3) )
    {
        norm = 1;   /*Z=11402*/
        order = 1;   /*Z=11403*/
        /* homogeneous */   /*Z=11404*/
        if ( cs==0 )
        {   /*Z=11405*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=11406*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11407*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11408*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11409*/
                fkv[n] = fkv[n-1]*n;   /*Z=11410*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11411*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11412*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11413*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=11414*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=11415*/
                if ( dim==1 )
                { /* cylinder */   /*Z=11416*/
                    /* P(q)-coefficients */   /*Z=11417*/
                    CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(n+1)*gam3[n]*fkv[n]);       /* P||(q) */   /*Z=11418*/
                    CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  /* P-(q) */   /*Z=11419*/
                    /* F(q)-coefficients */   /*Z=11420*/
                    sump = 0.0;   /*Z=11421*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11422*/
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11423*/
                    }   /*Z=11424*/
                    CR->carr1f[n] = M_PI*xln[n]*sump/4.0;   /*Z=11425*/
                    sump = 0.0;   /*Z=11426*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11427*/
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11428*/
                    }   /*Z=11429*/
                    CR->carr4f[n] = xrn[n-1]*sump;   /*Z=11430*/
                }   /*Z=11431*/
                if ( dim==2 )
                { /* disk */   /*Z=11432*/
                    /* P(q)-coefficients */   /*Z=11433*/
                    CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);    /* P-(q) */   /*Z=11434*/
                    CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);           /* P||(q) */   /*Z=11435*/
                    /* F(q)-coefficients */   /*Z=11436*/
                    sump = 0.0;   /*Z=11437*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11438*/
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11439*/
                    }   /*Z=11440*/
                    CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=11441*/
                    sump = 0.0;   /*Z=11442*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11443*/
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11444*/
                    }   /*Z=11445*/
                    CR->carr4f[n] = M_PI*xln[n]*sump/4.0;   /*Z=11446*/
                }   /*Z=11447*/
                if ( search1 )
                {   /*Z=11448*/
                    if ( fabs(CR->carr1p[n])<1e-50 )
                    {   /*Z=11449*/
                        n1 = n;   /*Z=11450*/
                        search1 = false;   /*Z=11451*/
                    }   /*Z=11452*/
                }   /*Z=11453*/
                if ( search4 )
                {   /*Z=11454*/
                    if ( fabs(CR->carr4p[n])<1e-50 )
                    {   /*Z=11455*/
                        n4 = n;   /*Z=11456*/
                        search4 = false;   /*Z=11457*/
                    }   /*Z=11458*/
                }   /*Z=11459*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=11460*/
                    if ( n<n1f ) n1f = n;   /*Z=11461*/
                }   /*Z=11462*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=11463*/
                    if ( n<n4f ) n4f = n;   /*Z=11464*/
                }   /*Z=11465*/
            }   /*Z=11466*/
        }   /*Z=11467*/

        /* core/shell */   /*Z=11469*/
        if ( cs==1 )
        {   /*Z=11470*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=11471*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11472*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11473*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11474*/
                fkv[n] = fkv[n-1]*n;   /*Z=11475*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11476*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11477*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11478*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=11479*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=11480*/
                xrmn_n = -xrmn_n*xrm2z;   /*Z=11481*/
                pn[n] = pn[n-1]*p*p;   /*Z=11482*/
                if ( dim==1 )
                { /* cylinder */   /*Z=11483*/
                    /* P(q)-coefficients */   /*Z=11484*/
                    /* longitudinal */   /*Z=11485*/
                    CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=11486*/
                    /* cross-sectional */   /*Z=11487*/
                    /* F121 */   /*Z=11488*/
                    CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=11489*/
                    /* F122 */   /*Z=11490*/
                    sump = 0.0;   /*Z=11491*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11492*/
                        sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11493*/
                    }   /*Z=11494*/
                    CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=11495*/
                    /* F123 */   /*Z=11496*/
                    CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=11497*/

                    /* F(q)-coefficients */   /*Z=11499*/
                    /* longitudinal */   /*Z=11500*/
                    sump = 0.0;   /*Z=11501*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11502*/
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11503*/
                    }   /*Z=11504*/
                    CR->carr1f[n] = M_PI*xln[n]*sump/4.0;   /*Z=11505*/
                    /* cross-sectional */   /*Z=11506*/
                    /* F121 */   /*Z=11507*/
                    sump = 0.0;   /*Z=11508*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11509*/
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11510*/
                    }   /*Z=11511*/
                    CR->carr4f[n] = xrn[n-1]*sump;   /*Z=11512*/
                    /* F122 */   /*Z=11513*/
                    sump = 0.0;   /*Z=11514*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11515*/
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11516*/
                    }   /*Z=11517*/
                    CR->carr5f[n] = xrmn_n*sump;   /*Z=11518*/
                    /* F123 */   /*Z=11519*/
                    CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=11520*/
                }   /*Z=11521*/

                if ( dim==2 )
                { /* disk */   /*Z=11523*/
                    /* P(q)-coefficients */   /*Z=11524*/
                    /* longitudinal */   /*Z=11525*/
                    CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=11526*/
                    /* cross-sectional */   /*Z=11527*/
                    /* F121 */   /*Z=11528*/
                    CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=11529*/
                    /* F122 */   /*Z=11530*/
                    /* sump:=0.0; */   /*Z=11531*/
                    /*    for m:=0 to n do begin */   /*Z=11532*/
                    /*    sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]); */   /*Z=11533*/
                    /* end; */   /*Z=11534*/
                    /* carr5p[n]:=pi*z12v[n]*xrmn_n*sump/4; */   /*Z=11535*/

                    /* F122 */   /*Z=11537*/
                    sump = 0.0;   /*Z=11538*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11539*/
                        sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11540*/
                    }   /*Z=11541*/
                    CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;   /*Z=11542*/

                    /* F123 */   /*Z=11545*/
                    CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=11546*/
                    /* F(q)-coefficients */   /*Z=11547*/
                    /* longitudinal */   /*Z=11548*/
                    sump = 0.0;   /*Z=11549*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11550*/
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11551*/
                    }   /*Z=11552*/
                    CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=11553*/
                    /* cross-sectional */   /*Z=11554*/
                    /* F121 */   /*Z=11555*/
                    sump = 0.0;   /*Z=11556*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11557*/
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11558*/
                    }   /*Z=11559*/
                    CR->carr4f[n] = M_PI*xrn[n-1]*sump/4.0;   /*Z=11560*/
                    /* F122 */   /*Z=11561*/
                    sump = 0.0;   /*Z=11562*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11563*/
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11564*/
                    }   /*Z=11565*/
                    CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;   /*Z=11566*/
                    /* F123 */   /*Z=11567*/
                    CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=11568*/
                }   /*Z=11569*/
                if ( fabs(CR->carr1p[n])<min )
                {   /*Z=11570*/
                    if ( n<n1 ) n1 = n;   /*Z=11571*/
                }   /*Z=11572*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=11573*/
                    if ( n<n4 ) n4 = n;   /*Z=11574*/
                }   /*Z=11575*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z=11576*/
                    if ( n<n5 ) n5 = n;   /*Z=11577*/
                }   /*Z=11578*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z=11579*/
                    if ( n<n6 ) n6 = n;   /*Z=11580*/
                }   /*Z=11581*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=11582*/
                    if ( n<n1f ) n1f = n;   /*Z=11583*/
                }   /*Z=11584*/
                if ( abs(CR->carr4f[n])<min )
                {   /*Z=11585*/
                    if ( n<n4f ) n4f = n;   /*Z=11586*/
                }   /*Z=11587*/
                if ( fabs(CR->carr5f[n])<min )
                {   /*Z=11588*/
                    if ( n<n5f ) n5f = n;   /*Z=11589*/
                }   /*Z=11590*/
                if ( fabs(CR->carr6f[n])<min )
                {   /*Z=11591*/
                    if ( n<n6f ) n6f = n;   /*Z=11592*/
                }   /*Z=11593*/
            }   /*Z=11594*/
        }  /* of core/shell */   /*Z=11595*/

        /* inhomogeneous core/shell */   /*Z=11597*/
        if ( cs==2 )
        {   /*Z=11598*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=11599*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11600*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11601*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11602*/
                fkv[n] = fkv[n-1]*n;   /*Z=11603*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11604*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11605*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11606*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=11607*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=11608*/
                xrmn_n = -xrmn_n*xrm2z;   /*Z=11609*/
                pn[n] = pn[n-1]*p*p;   /*Z=11610*/
                if ( dim==1 )
                { /* cylinder */   /*Z=11611*/
                    /* P(q)-coefficients */   /*Z=11612*/
                    /* longitudinal */   /*Z=11613*/
                    CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=11614*/

                    /* cross-sectional P(q) */   /*Z=11616*/
                    CR->carr4p[n] = pow(4.,n+1.)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=11617*/
                    sump = 0.0;   /*Z=11618*/
                    sump1 = 0.0;   /*Z=11619*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11620*/
                        sumi = 1/((m+1-alfa/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);   /*Z=11621*/
                        sump = sump+pn[n-m]*sumi;   /*Z=11622*/
                        sump1 = sump1+sumi;   /*Z=11623*/
                    }   /*Z=11624*/
                    CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;   /*Z=11625*/
                    CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n-1]*sump1;   /*Z=11626*/
                    sump = 0.0;   /*Z=11627*/
                    sump1 = 0.0;   /*Z=11628*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11629*/
                        sumi = 1/((n-m+1-alfa/2.0)*(m+1-alfa/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);   /*Z=11630*/
                        sump = sump+sumi;   /*Z=11631*/
                        sump1 = sump1+pn[n-m]*sumi;   /*Z=11632*/
                    }   /*Z=11633*/
                    CR->carr7p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump;   /*Z=11634*/
                    CR->carr8p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump1;   /*Z=11635*/
                    CR->carr9p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrn[n-1]*sump;   /*Z=11636*/

                    /* cross-sectional */   /*Z=11638*/

                    /* F(q)-coefficients */   /*Z=11650*/
                    /* longitudinal */   /*Z=11651*/
                    sump = 0.0;   /*Z=11652*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11653*/
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11654*/
                    }   /*Z=11655*/
                    CR->carr1f[n] = M_PI*xln[n]*sump/4.0;   /*Z=11656*/
                    /* cross-sectional */   /*Z=11657*/
                    CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);   /*Z=11658*/
                    CR->carr5f[n] = (1-alfa/2.0)*z12v[n]*xrmn_n/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=11659*/
                    CR->carr6f[n] = (1-alfa/2.0)*z12v[n]*xrn[n]/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=11660*/
                }   /*Z=11661*/

                if ( dim==2 )
                { /* disk */   /*Z=11663*/
                    /* P(q)-coefficients */   /*Z=11664*/
                    /* longitudinal */   /*Z=11665*/
                    CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=11666*/

                    /* cross-sectional P(q) */   /*Z=11668*/
                    CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.,1.0*n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);   /*Z=11669*/
                    sump = 0.0;   /*Z=11670*/
                    sump1 = 0.0;   /*Z=11671*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11672*/
                        sumi = (m+1/2.0)/((m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=11673*/
                        sump = sump+pn[n-m]*sumi;   /*Z=11674*/
                        sump1 = sump1+sumi;   /*Z=11675*/
                    }   /*Z=11676*/
                    CR->carr5p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=11677*/
                    CR->carr6p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrn[n-1]*sump1;   /*Z=11678*/
                    sump = 0.0;   /*Z=11679*/
                    sump1 = 0.0;   /*Z=11680*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11681*/
                        sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-alfa/2.0)*(m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);   /*Z=11682*/
                        sump = sump+sumi;   /*Z=11683*/
                        sump1 = sump1+pn[n-m]*sumi;   /*Z=11684*/
                    }   /*Z=11685*/
                    CR->carr7p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=11686*/
                    CR->carr8p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump1;   /*Z=11687*/
                    CR->carr9p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrn[n-1]*sump;   /*Z=11688*/

                    /* cross-sectional P(q) */   /*Z=11691*/
                    /* F121 */   /*Z=11692*/

                    /* F(q)-coefficients */   /*Z=11710*/
                    /* longitudinal */   /*Z=11711*/
                    sump = 0.0;   /*Z=11712*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11713*/
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11714*/
                    }   /*Z=11715*/
                    CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=11716*/
                    /* cross-sectional */   /*Z=11717*/
                    CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);   /*Z=11718*/
                    CR->carr5f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=11719*/
                    CR->carr6f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrn[n-1]*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=11720*/
                }   /*Z=11722*/
                if ( fabs(CR->carr1p[n])<min )
                {   /*Z=11723*/
                    if ( n<n1 ) n1 = n;   /*Z=11724*/
                }   /*Z=11725*/
                if ( fabs(CR->carr4p[n])<min )
                {   /*Z=11726*/
                    if ( n<n4 ) n4 = n;   /*Z=11727*/
                }   /*Z=11728*/
                if ( fabs(CR->carr5p[n])<min )
                {   /*Z=11729*/
                    if ( n<n5 ) n5 = n;   /*Z=11730*/
                }   /*Z=11731*/
                if ( fabs(CR->carr6p[n])<min )
                {   /*Z=11732*/
                    if ( n<n6 ) n6 = n;   /*Z=11733*/
                }   /*Z=11734*/
                if ( fabs(CR->carr7p[n])<min )
                {   /*Z=11735*/
                    if ( n<n7 ) n7 = n;   /*Z=11736*/
                }   /*Z=11737*/
                if ( abs(CR->carr8p[n])<min )
                {   /*Z=11738*/
                    if ( n<n8 ) n8 = n;   /*Z=11739*/
                }   /*Z=11740*/
                if ( fabs(CR->carr9p[n])<min )
                {   /*Z=11741*/
                    if ( n<n9 ) n9 = n;   /*Z=11742*/
                }   /*Z=11743*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=11744*/
                    if ( n<n1f ) n1f = n;   /*Z=11745*/
                }   /*Z=11746*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=11747*/
                    if ( n<n4f ) n4f = n;   /*Z=11748*/
                }   /*Z=11749*/
                if ( fabs(CR->carr5f[n])<min )
                {   /*Z=11750*/
                    if ( n<n5f ) n5f = n;   /*Z=11751*/
                }   /*Z=11752*/
                if ( fabs(CR->carr6f[n])<min )
                {   /*Z=11753*/
                    if ( n<n6f ) n6f = n;   /*Z=11754*/
                }   /*Z=11755*/
                if ( fabs(CR->carr7f[n])<min )
                {   /*Z=11756*/
                    if ( n<n7f ) n7f = n;   /*Z=11757*/
                }   /*Z=11758*/
                if ( fabs(CR->carr8f[n])<min )
                {   /*Z=11759*/
                    if ( n<n8f ) n8f = n;   /*Z=11760*/
                }   /*Z=11761*/
                if ( fabs(CR->carr9f[n])<min )
                {   /*Z=11762*/
                    if ( n<n9f ) n9f = n;   /*Z=11763*/
                }   /*Z=11764*/
            }   /*Z=11765*/
        }  /* of inhomogeneous core/shell */   /*Z=11766*/

        /* myelin */   /*Z=11769*/
        if ( (cs==3) || (cs==4) )
        {   /*Z=11770*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=11771*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11772*/
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11773*/
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11774*/
                fkv[n] = fkv[n-1]*n;   /*Z=11775*/
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11776*/
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11777*/
                /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11778*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=11779*/
                xrn[n] = -xrn[n-1]*x12zm;   /*Z=11780*/
                if ( dim==1 )
                { /* cylinder */   /*Z=11781*/
                    /* P(q)-coefficients */   /*Z=11782*/
                    CR->carr1p[n] = sqrt(M_PI)*pow(4.,1.0*n)*z12vl[n]*xln[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=11783*/

                    CR->carr4p[n] = z12v[n]*xrn[n];   /*Z=11785*/

                    /* F(q)-coefficients */   /*Z=11787*/
                    sump = 0.0;   /*Z=11788*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11789*/
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11790*/
                    }   /*Z=11791*/
                    CR->carr1f[n] = M_PI*xln[n]*sump/4.0;   /*Z=11792*/
                    sump = 0.0;   /*Z=11793*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11794*/
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11795*/
                    }   /*Z=11796*/
                    CR->carr4f[n] = xrn[n-1]*sump;   /*Z=11797*/
                }   /*Z=11798*/
                if ( dim==2 )
                { /* disk */   /*Z=11799*/
                    /* P(q)-coefficients */   /*Z=11800*/
                    CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=11801*/
                    CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=11802*/
                    /* F(q)-coefficients */   /*Z=11803*/
                    sump = 0.0;   /*Z=11804*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11805*/
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11806*/
                    }   /*Z=11807*/
                    CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=11808*/
                    sump = 0.0;   /*Z=11809*/
                    for ( m=0; m<=n; m++ )
                    {   /*Z=11810*/
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11811*/
                    }   /*Z=11812*/
                    CR->carr4f[n] = M_PI*xln[n]*sump/4.0;   /*Z=11813*/
                }   /*Z=11814*/
                if ( search1 )
                {   /*Z=11815*/
                    if ( fabs(CR->carr1p[n])<1e-50 )
                    {   /*Z=11816*/
                        n1 = n;   /*Z=11817*/
                        search1 = false;   /*Z=11818*/
                    }   /*Z=11819*/
                }   /*Z=11820*/
                if ( search4 )
                {   /*Z=11821*/
                    if ( fabs(CR->carr4p[n])<1e-50 )
                    {   /*Z=11822*/
                        n4 = n;   /*Z=11823*/
                        search4 = false;   /*Z=11824*/
                    }   /*Z=11825*/
                }   /*Z=11826*/
                if ( fabs(CR->carr1f[n])<min )
                {   /*Z=11827*/
                    if ( n<n1f ) n1f = n;   /*Z=11828*/
                }   /*Z=11829*/
                if ( fabs(CR->carr4f[n])<min )
                {   /*Z=11830*/
                    if ( n<n4f ) n4f = n;   /*Z=11831*/
                }   /*Z=11832*/
            }   /*Z=11833*/
        }  /* of myelin */   /*Z=11834*/
    }   /*Z=11838*/

    /*** orientational distribution for cylinders and disks ***/  /*Z0311=11509*/ /*Z=11840*/
    if ( (ordis==0) && (dim!=3) )   /*Z0311=11510*/
    {
        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,2,0,0,0,0,CR->carr1p,norm);  /*Z=11842*/
        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,3,0,0,0,0,CR->carr1p,order);  /*Z=11843*/
        order = order/norm;  /*Z=11844*/

        DBG( qDebug() << "coefficients()" << "order" << order << "norm" << norm << "cho1"<<cho1; )

        /* use phi=0 and rotate qx/qy-axis */  /*Z=11846*/
        if ( cho1==1 )
        {
            phi = 0;   /*Z=11848*/
            p11 = -cos(phi*M_PI/180.0)*cos(theta*M_PI/180.0);       /*  = -cos(theta*pi/180); */   /*Z=11849*/
            p12 = sin(phi*M_PI/180.0);                          /*  = 0; */   /*Z=11850*/
            p13 = cos(phi*M_PI/180.0)*sin(theta*M_PI/180.0);        /*  =  sin(theta*pi/180); */   /*Z=11851*/
            p21 = -cos(phi*M_PI/180.0);                         /*  = -1; */   /*Z=11852*/
            p22 = -sin(phi*M_PI/180.0)*cos(theta*M_PI/180.0);       /*  = 0; */   /*Z=11853*/
            p23 = sin(phi*M_PI/180.0)*sin(theta*M_PI/180.0);        /*  = 0; */   /*Z=11854*/
            p31 = -sin(theta*M_PI/180.0);                       /*  = 0; */   /*Z=11855*/
            p32 = 0;   /*Z=11856*/
            p33 = -cos(theta*M_PI/180.0);                       /*  = -cos(theta*pi/180); */   /*Z=11857*/

            for ( n=0; n<=nmax; n++ )
            {   /*Z=11859*/
                for ( m=0; m<=nmax; m++ )
                {   /*Z=11860*/
                    qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,0,0,2*n,2*m,CR->carr1p,intl);   /*Z=11861*/
                    intlar[n][m] = intl;   /*Z=11862*/
                }   /*Z=11863*/
            }   /*Z=11864*/

            /* cylinder form factor coefficients */   /*Z=11875*/
            if ( dim==1 )
            {   /*Z=11876*/
                /* homogeneous */   /*Z=11877*/
                if ( cs==0 )
                {   /*Z=11878*/
                    i = 2;   /*Z=11879*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=11880*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11881*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11882*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11883*/
                        fkv[n] = fkv[n-1]*n;   /*Z=11884*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11885*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11886*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11887*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=11888*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=11889*/
                        /* longitudinal */   /*Z=11890*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=11891*/
                            sump = 0.0;   /*Z=11892*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=11893*/
                                sump = sump+pow(p11*p11,1.0*ll)*pow(p13*p13,1.0*m-ll)*intlar[m-ll][n-m+ll]/(pow(4.,1.0*ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*norm);   /*Z=11894*/
                            /* carr1pm[i]:=power(4,m)*power(p21*p21,n-m)*sump/(fkv[n-m]); */   /*Z=11895*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=11896*/
                            /* i:=i+1; */   /*Z=11897*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*pow(p21*p21,1.0*n-m)*sump/(fkv[n-m]);   /*Z=11898, cs=0,dim=1,cho1=1*/
                        }   /*Z=11899*/
                        /* P(q)-coefficients */   /*Z=11900*/
                        CR->carr1p[n] = pow(4.,1.0*n)*((z12vl[n]*xln[n])/((2.*n+1)*(n+1)));   /*Z=11901*/
                        //D1P( std::cerr << "CARR1P(a) " << n << " " << CR->carr1p[n] << " (" << search1 << ") " << n1 << " | "
                        //               << z12vl[n] << " " << xln[n] << " " << pow(4.,1.0*n) << " " << (2.*n+1)*(n+1) << std::endl; )
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=11902*/
                        /* F(q)-coefficients */   /*Z=11903*/
                        sump = 0.0;   /*Z=11904*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11905*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=11906*/
                        CR->carr2f[n] = CR->carr2p[n];   /*Z=11907*/
                        /* cross-sectional */   /*Z=11908*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);   /*Z=11909*/
                        //  carr4p[n]:= 4*(n+1/2  )*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);
                        D4P( qDebug() << "CARR4P(a)" << n << CR->carr4p[n]; )
                        sump = 0.0;   /*Z=11910*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=11911*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=11912*/
                        if ( search1 )
                        {   /*Z=11913*/
                            if ( fabs(CR->carr1p[n])<1e-50 || fabs(CR->carr1p[n])>1e+50 ) //NEU
                            {   /*Z=11914*/
                                n1 = n;   /*Z=11915*/
                                search1 = false;   /*Z=11916*/
                            }   /*Z=11917*/
                        }   /*Z=11918*/
                        if ( search4 )
                        {   /*Z=11919*/
                            if ( fabs(CR->carr4p[n])<1e-50 || fabs(CR->carr4p[n])>1e+50 ) //NEU
                            {   /*Z=11920*/
                                n4 = n;   /*Z=11921*/
                                search4 = false;   /*Z=11922*/
                            }   /*Z=11923*/
                        }   /*Z=11924*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=11925*/
                            if ( n<n1f ) n1f = n;   /*Z=11926*/
                        }   /*Z=11927*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=11928*/
                            if ( n<n4f ) n4f = n;   /*Z=11929*/
                        }   /*Z=11930*/
                    }   /*Z=11931*/
                }  /* of homogeneous cylinder */   /*Z=11932*/

                /* core/shell */   /*Z=11934*/
                if ( cs==1 )
                {   /*Z=11935*/
                    i = 2;   /*Z=11936*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=11937*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=11938*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=11939*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=11940*/
                        fkv[n] = fkv[n-1]*n;   /*Z=11941*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=11942*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=11943*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=11944*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=11945*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=11946*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=11947*/
                        pn[n] = pn[n-1]*p*p;   /*Z=11948*/
                        /* longitudinal */   /*Z=11949*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=11950*/
                            sump = 0.0;   /*Z=11951*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=11952*/
                                sump = sump+pow(p11*p11,1.0*ll)*pow(p13*p13,1.0*m-ll)*intlar[m-ll][n-m+ll]/(pow(4.,1.0*ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*norm);   /*Z=11953*/
                            /* carr1pm[i]:=power(4,m)*sump/(fkv[n-m]); */   /*Z=11954*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=11955*/
                            /* i:=i+1; */   /*Z=11956*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*sump/(fkv[n-m]);   /*Z=11957, cs=1,dim=1,cho1=1*/
                        }   /*Z=11958*/
                        /* P(q)-coefficients */   /*Z=11959*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=11960*/
                        D1P( qDebug() << "CARR1P(b)" << n << CR->carr1p[n]; )
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=11961*/
                        /* F(q)-coefficients */   /*Z=11962*/
                        sump = 0.0;   /*Z=11963*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=11964*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=11965*/
                        CR->carr2f[n] = CR->carr2p[n];   /*Z=11966*/

                        /* P(q)-coefficients */   /*Z=11968*/
                        /* cross-sectional */   /*Z=11969*/
                        /* F121 */   /*Z=11970*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=11971*/
                        //  carr4p[n]:= 4*(n+1/2  )*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);
                        //D4P( qDebug() << "CARR4P(b)" << n << CR->carr4p[n]; )
                        /* F122 */   /*Z=11972*/
                        sump = 0.0;   /*Z=11973*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=11974*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11975*/
                        }   /*Z=11976*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=11977*/
                        /* F123 */   /*Z=11978*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=11979*/

                        /* F(q)-coefficients */   /*Z=11981*/
                        /* cross-sectional */   /*Z=11982*/
                        /* F121 */   /*Z=11983*/
                        sump = 0.0;   /*Z=11984*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=11985*/
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11986*/
                        }   /*Z=11987*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=11988*/
                        /* F122 */   /*Z=11989*/
                        sump = 0.0;   /*Z=11990*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=11991*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=11992*/
                        }   /*Z=11993*/
                        CR->carr5f[n] = xrmn_n*sump;   /*Z=11994*/
                        /* F123 */   /*Z=11995*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=11996*/

                        if ( search1 )
                        {   /*Z=11998*/
                            if ( fabs(CR->carr1p[n])<1e-50 || fabs(CR->carr1p[n])>1e+50 ) //NEU
                            {   /*Z=11999*/
                                n1 = n;   /*Z=12000*/
                                search1 = false;   /*Z=12001*/
                            }   /*Z=12002*/
                        }   /*Z=12003*/
                        if ( search4 )
                        {   /*Z=12004*/
                            if ( fabs(CR->carr4p[n])<1e-50 || fabs(CR->carr4p[n])>1e+50 ) //NEU
                            {   /*Z=12005*/
                                n4 = n;   /*Z=12006*/
                                search4 = false;   /*Z=12007*/
                            }   /*Z=12008*/
                        }   /*Z=12009*/
                        if ( search5 )
                        {   /*Z=12010*/
                            if ( fabs(CR->carr5p[n])<1e-50 || fabs(CR->carr5p[n])>1e+50 ) //NEU
                            {   /*Z=12011*/
                                n5 = n;   /*Z=12012*/
                                search5 = false;   /*Z=12013*/
                            }   /*Z=12014*/
                        }   /*Z=12015*/
                        if ( search6 )
                        {   /*Z=12016*/
                            if ( fabs(CR->carr6p[n])<1e-50 || fabs(CR->carr6p[n])>1e+50 ) //NEU
                            {   /*Z=12017*/
                                n6 = n;   /*Z=12018*/
                                search6 = false;   /*Z=12019*/
                            }   /*Z=12020*/
                        }   /*Z=12021*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12022*/
                            if ( n<n1f ) n1f = n;   /*Z=12023*/
                        }   /*Z=12024*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12025*/
                            if ( n<n4f ) n4f = n;   /*Z=12026*/
                        }   /*Z=12027*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=12028*/
                            if ( n<n5f ) n5f = n;   /*Z=12029*/
                        }   /*Z=12030*/
                        if ( abs(CR->carr6f[n])<min )
                        {   /*Z=12031*/
                            if ( n<n6f ) n6f = n;   /*Z=12032*/
                        }   /*Z=12033*/
                    }  /* n-loop */   /*Z=12034*/
                }  /* homogeneous loop */   /*Z=12035*/

                /* inhomogeneous core/shell */   /*Z=12039*/
                if ( cs==2 )
                {   /*Z=12040*/
                    i = 2;   /*Z=12041*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12042*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12043*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=12044*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=12045*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12046*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12047*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12048*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12049*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12050*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=12051*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=12052*/
                        pn[n] = pn[n-1]*p*p;   /*Z=12053*/
                        /* longitudinal */   /*Z=12054*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12055*/
                            sump = 0.0;   /*Z=12056*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=12057*/
                                sump = sump+pow(p11*p11,1.0*ll)*pow(p13*p13,1.0*m-ll)*intlar[m-ll][n-m+ll]/(pow(4.,1.0*ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*norm);   /*Z=12058*/
                            /* carr1pm[i]:=power(4,m)*sump/(fkv[n-m]); */   /*Z=12059*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12060*/
                            /* i:=i+1; */   /*Z=12061*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*sump/(fkv[n-m]);   /*Z=12062, cs=2,dim=1,cho1=1*/
                        }   /*Z=12063*/
                        /* P(q)-coefficients */   /*Z=12064*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=12065*/
                        D1P( qDebug() << "CARR1P(c)" << n << CR->carr1p[n]; )
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=12066*/
                        /* F(q)-coefficients */   /*Z=12067*/
                        sump = 0.0;   /*Z=12068*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12069*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=12070*/
                        CR->carr2f[n] = CR->carr2p[n];   /*Z=12071*/

                        /* cross-sectional P(q) */   /*Z=12074*/
                        CR->carr4p[n] = pow(4.,n+1.)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=12075*/
                        //  carr4p[n]:=power(4,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(pi)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);
                        //D4P( qDebug() << "CARR4P(c)" << n << CR->carr4p[n]; )
                        sump = 0.0;   /*Z=12076*/
                        sump1 = 0.0;   /*Z=12077*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12078*/
                            sumi = 1/((m+1-alfa/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);   /*Z=12079*/
                            sump = sump+pn[n-m]*sumi;   /*Z=12080*/
                            sump1 = sump1+sumi;   /*Z=12081*/
                        }   /*Z=12082*/
                        CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;   /*Z=12083*/
                        CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n-1]*sump1;   /*Z=12084*/
                        sump = 0.0;   /*Z=12085*/
                        sump1 = 0.0;   /*Z=12086*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12087*/
                            sumi = 1/((n-m+1-alfa/2.0)*(m+1-alfa/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);   /*Z=12088*/
                            sump = sump+sumi;   /*Z=12089*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=12090*/
                        }   /*Z=12091*/
                        CR->carr7p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump;   /*Z=12092*/
                        CR->carr8p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump1;   /*Z=12093*/
                        CR->carr9p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrn[n-1]*sump;   /*Z=12094*/

                        /* F(q)-coefficients */   /*Z=12111*/
                        CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);   /*Z=12112*/
                        CR->carr5f[n] = (1-alfa/2.0)*z12v[n]*xrmn_n/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=12113*/
                        CR->carr6f[n] = (1-alfa/2.0)*z12v[n]*xrn[n]/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=12114*/

                        if ( search1 )
                        {   /*Z=12116*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=12117*/
                                n1 = n;   /*Z=12118*/
                                search1 = false;   /*Z=12119*/
                            }   /*Z=12120*/
                        }   /*Z=12121*/
                        if ( search4 )
                        {   /*Z=12122*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=12123*/
                                n4 = n;   /*Z=12124*/
                                search4 = false;   /*Z=12125*/
                            }   /*Z=12126*/
                        }   /*Z=12127*/
                        if ( search5 )
                        {   /*Z=12128*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=12129*/
                                n5 = n;   /*Z=12130*/
                                search5 = false;   /*Z=12131*/
                            }   /*Z=12132*/
                        }   /*Z=12133*/
                        if ( search6 )
                        {   /*Z=12134*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=12135*/
                                n6 = n;   /*Z=12136*/
                                search6 = false;   /*Z=12137*/
                            }   /*Z=12138*/
                        }   /*Z=12139*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=12140*/
                            if ( n<n7 ) n7 = n;   /*Z=12141*/
                        }   /*Z=12142*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=12143*/
                            if ( n<n8 ) n8 = n;   /*Z=12144*/
                        }   /*Z=12145*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=12146*/
                            if ( n<n9 ) n9 = n;   /*Z=12147*/
                        }   /*Z=12148*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12149*/
                            if ( n<n1f ) n1f = n;   /*Z=12150*/
                        }   /*Z=12151*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12152*/
                            if ( n<n4f ) n4f = n;   /*Z=12153*/
                        }   /*Z=12154*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=12155*/
                            if ( n<n5f ) n5f = n;   /*Z=12156*/
                        }   /*Z=12157*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=12158*/
                            if ( n<n6f ) n6f = n;   /*Z=12159*/
                        }   /*Z=12160*/
                        if ( fabs(CR->carr7f[n])<min )
                        {   /*Z=12161*/
                            if ( n<n7f ) n7f = n;   /*Z=12162*/
                        }   /*Z=12163*/
                        if ( fabs(CR->carr8f[n])<min )
                        {   /*Z=12164*/
                            if ( n<n8f ) n8f = n;   /*Z=12165*/
                        }   /*Z=12166*/
                        if ( fabs(CR->carr9f[n])<min )
                        {   /*Z=12167*/
                            if ( n<n9f ) n9f = n;   /*Z=12168*/
                        }   /*Z=12169*/
                    }   /* of n-loop */   /*Z=12170*/
                }  /* of inhomogeneous core/shell */   /*Z=12171*/

                /* myelin */   /*Z=12173*/
                if ( (cs==3) || (cs==4) )
                {   /*Z=12174*/
                    i = 2;   /*Z=12175*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12176*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12177*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=12178*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=12179*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12180*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12181*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12182*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12183*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12184*/
                        xrn[n] = -xrn[n-1]*x12zm;   /*Z=12185*/
                        /* longitudinal */   /*Z=12186*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12187*/
                            sump = 0.0;   /*Z=12188*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=12189*/
                                sump = sump+pow(p11*p11,1.0*ll)*pow(p13*p13,1.0*m-ll)*intlar[m-ll][n-m+ll]/(pow(4.,1.0*ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*norm);   /*Z=12190*/
                            /* carr1pm[i]:=power(4,m)*sump/(fkv[n-m]); */   /*Z=12191*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12192*/
                            /* i:=i+1; */   /*Z=12193*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*sump/(fkv[n-m]);   /*Z=12194, cs=3|4,dim=1,cho1=1*/
                        }   /*Z=12195*/
                        /* P(q)-coefficients */   /*Z=12196*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=12197*/
                        D1P( qDebug() << "CARR1P(d)" << n << CR->carr1p[n]; )
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=12198*/
                        /* F(q)-coefficients */   /*Z=12199*/
                        sump = 0.0;   /*Z=12200*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12201*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=12202*/
                        CR->carr2f[n] = CR->carr2p[n];   /*Z=12203*/
                        /* cross-sectional */   /*Z=12204*/
                        CR->carr4p[n] = z12v[n]*xrn[n];   /*Z=12205*/
                        //D4P( qDebug() << "CARR4P(d)" << n << CR->carr4p[n]; )

                        sump = 0.0;   /*Z=12207*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=12208*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=12209*/
                        if ( search1 )
                        {   /*Z=12210*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=12211*/
                                n1 = n;   /*Z=12212*/
                                search1 = false;   /*Z=12213*/
                            }   /*Z=12214*/
                        }   /*Z=12215*/
                        if ( search4 )
                        {   /*Z=12216*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=12217*/
                                n4 = n;   /*Z=12218*/
                                search4 = false;   /*Z=12219*/
                            }   /*Z=12220*/
                        }   /*Z=12221*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12222*/
                            if ( n<n1f ) n1f = n;   /*Z=12223*/
                        }   /*Z=12224*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12225*/
                            if ( n<n4f ) n4f = n;   /*Z=12226*/
                        }   /*Z=12227*/
                    }   /*Z=12228*/
                }        /* of myelin */   /*Z=12229*/
            }   /* of cylinder */   /*Z=12230*/

            //if (dim=1) then begin   /*Z=12232*/

            /* disk form factor coefficients */   /*Z=12257*/
            if ( dim==2 )
            {   /*Z=12258*/
                /* homogeneous */   /*Z=12259*/
                if ( cs==0 )
                {   /*Z=12260*/
                    i = 2;   /*Z=12261*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12262*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12263*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=12264*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12265*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12266*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12267*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12268*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12269*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=12270*/
                        /* longitudinal */   /*Z=12271*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12272*/
                            sump = 0.0;   /*Z=12273*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=12274*/
                                sump = sump+pow(p11*p11,1.0*ll)*pow(p13*p13,1.0*m-ll)*intlar[m-ll][n-m+ll]/(pow(4.,1.0*ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*norm);   /*Z=12275*/
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]); */   /*Z=12276*/
                            CR->carr22pm[n][m] = pow(4.,1.0*m)*sump/(fkv[n-m]);   /*Z=12277*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]); */   /*Z=12278*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(pow(4.,1.0*m)*fkv[m]*fkv[n-m]);   /*Z=12279, cs=0,dim=2,cho1=1*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=12280*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12281*/
                            i = i+1;   /*Z=12282*/
                        }   /*Z=12283*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=12284*/
                        sump1 = 0.0;   /*Z=12285*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12286*/
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=12287*/
                        }   /*Z=12288*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump1;   /*Z=12289*/
                        /*Z=12290*/
                        /* cross-sectional */   /*Z=12291*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=12292*/
                        sump = 0.0;   /*Z=12293*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12294*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12295*/
                        }   /*Z=12296*/
                        CR->carr4f[n] = M_PI*xrn[n-1]*sump/4.0;   /*Z=12297*/
                        /*Z=12298*/
                        /* series for <...> integration */   /*Z=12299*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=12300*/
                        sump = 0.0;   /*Z=12301*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12302*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=12303*/
                        }   /*Z=12304*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=12305*/
                        /*Z=12306*/
                        if ( search1 )
                        {   /*Z=12307*/
                            if ( fabs(CR->carr1p[n])<1e-50 || fabs(CR->carr1p[n])>1e+50 ) //NEU
                            {   /*Z=12308*/
                                n1 = n;   /*Z=12309*/
                                search1 = false;   /*Z=12310*/
                            }   /*Z=12311*/
                        }   /*Z=12312*/
                        if ( search4 )
                        {   /*Z=12313*/
                            if ( fabs(CR->carr4p[n])<1e-50 || fabs(CR->carr4p[n])>1e+50 ) //NEU
                            {   /*Z=12314*/
                                n4 = n;   /*Z=12315*/
                                search4 = false;   /*Z=12316*/
                            }   /*Z=12317*/
                        }   /*Z=12318*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12319*/
                            if ( n<n1f ) n1f = n;   /*Z=12320*/
                        }   /*Z=12321*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12322*/
                            if ( n<n4f ) n4f = n;   /*Z=12323*/
                        }   /*Z=12324*/
                    }   /*Z=12325*/
                }   /*Z=12326*/

                /* core/shell */   /*Z=12328*/
                if ( cs==1 )
                {   /*Z=12329*/
                    i = 2;   /*Z=12330*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12331*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12332*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=12333*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12334*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12335*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12336*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12337*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12338*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=12339*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=12340*/
                        pn[n] = pn[n-1]*p*p;   /*Z=12341*/
                        /* longitudinal */   /*Z=12342*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12343*/
                            sump = 0.0;   /*Z=12344*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=12345*/
                                sump = sump+pow(p11*p11,1.0*ll)*pow(p13*p13,1.0*m-ll)*intlar[m-ll][n-m+ll]/(pow(4.,1.0*ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*norm);   /*Z=12346*/
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]); */   /*Z=12347*/
                            CR->carr22pm[n][m] = pow(4.,1.0*m)*sump/(fkv[n-m]);   /*Z=12348*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]); */   /*Z=12349*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(pow(4.,1.0*m)*fkv[m]*fkv[n-m]);   /*Z=12350, cs=1,dim=2,cho1=1*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=12351*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12352*/
                            /* i:=i+1; */   /*Z=12353*/
                        }   /*Z=12354*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=12355*/
                        sump1 = 0.0;   /*Z=12356*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12357*/
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=12358*/
                        }   /*Z=12359*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump1;   /*Z=12360*/

                        /* P(q)-coefficients */   /*Z=12362*/
                        /* cross-sectional */   /*Z=12363*/
                        /* F121 */   /*Z=12364*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=12365*/
                        /* F122 */   /*Z=12366*/
                        sump = 0.0;   /*Z=12367*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12368*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=12369*/
                        }   /*Z=12370*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=12371*/

                        /* F122 */   /*Z=12374*/
                        sump = 0.0;   /*Z=12375*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12376*/
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12377*/
                        }   /*Z=12378*/
                        CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;   /*Z=12379*/

                        /* F123 */   /*Z=12382*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=12383*/
                        /* F121 */   /*Z=12384*/
                        sump = 0.0;   /*Z=12385*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12386*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12387*/
                        }   /*Z=12388*/
                        CR->carr4f[n] = M_PI*xrn[n-1]*sump/4.0;   /*Z=12389*/
                        /* F122 */   /*Z=12390*/
                        sump = 0.0;   /*Z=12391*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12392*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12393*/
                        }   /*Z=12394*/
                        CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;   /*Z=12395*/
                        /* F123 */   /*Z=12396*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=12397*/

                        /* series for <...> integration */   /*Z=12399*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=12400*/
                        sump = 0.0;   /*Z=12401*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12402*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=12403*/
                        }   /*Z=12404*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=12405*/

                        if ( search1 )
                        {   /*Z=12407*/
                            if ( fabs(CR->carr1p[n])<1e-50 || fabs(CR->carr1p[n])>1e+50 )
                            {   /*Z=12408*/
                                n1 = n;   /*Z=12409*/
                                search1 = false;   /*Z=12410*/
                            }   /*Z=12411*/
                        }   /*Z=12412*/
                        if ( search4 )
                        {   /*Z=12413*/
                            if ( fabs(CR->carr4p[n])<1e-50 || fabs(CR->carr4p[n])>1e+50 ) //NEU
                            {   /*Z=12414*/
                                n4 = n;   /*Z=12415*/
                                search4 = false;   /*Z=12416*/
                            }   /*Z=12417*/
                        }   /*Z=12418*/
                        if ( search5 )
                        {   /*Z=12419*/
                            if ( fabs(CR->carr5p[n])<1e-50 || fabs(CR->carr5p[n])>1e+50 ) //NEU
                            {   /*Z=12420*/
                                n5 = n;   /*Z=12421*/
                                search5 = false;   /*Z=12422*/
                            }   /*Z=12423*/
                        }   /*Z=12424*/
                        if ( search6 )
                        {   /*Z=12425*/
                            if ( fabs(CR->carr6p[n])<1e-50 || fabs(CR->carr6p[n])>1e+50 ) //NEU
                            {   /*Z=12426*/
                                n6 = n;   /*Z=12427*/
                                search6 = false;   /*Z=12428*/
                            }   /*Z=12429*/
                        }   /*Z=12430*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12431*/
                            if ( n<n1f ) n1f = n;   /*Z=12432*/
                        }   /*Z=12433*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12434*/
                            if ( n<n4f ) n4f = n;   /*Z=12435*/
                        }   /*Z=12436*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=12437*/
                            if ( n<n5f ) n5f = n;   /*Z=12438*/
                        }   /*Z=12439*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=12440*/
                            if ( n<n6f ) n6f = n;   /*Z=12441*/
                        }   /*Z=12442*/
                    } /* of n-loop */   /*Z=12443*/
                }  /* of core/shell */   /*Z=12444*/

                /* inhomogeneous core/shell */   /*Z=12446*/
                if ( cs==2 )
                {   /*Z=12447*/
                    i = 2;   /*Z=12448*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12449*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12450*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=12451*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12452*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12453*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12454*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12455*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12456*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=12457*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=12458*/
                        pn[n] = pn[n-1]*p*p;   /*Z=12459*/
                        /* longitudinal */   /*Z=12460*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12461*/
                            sump = 0.0;   /*Z=12462*/
                            for ( ll=0; ll<=m; ll++ )   /*Z=12463*/
                                sump = sump+pow(p11*p11,1.0*ll)*pow(p13*p13,1.0*m-ll)*intlar[m-ll][n-m+ll]/(pow(4.,1.0*ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*norm);   /*Z=12464*/
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]); */   /*Z=12465*/
                            CR->carr22pm[n][m] = pow(4.,1.0*m)*sump/(fkv[n-m]);   /*Z=12466*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]); */   /*Z=12467*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(pow(4.,1.0*m)*fkv[m]*fkv[n-m]);   /*Z=12468, cs=2,dim=2,cho1=1*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=12469*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12470*/
                            /* i:=i+1; */   /*Z=12471*/
                        }   /*Z=12472*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=12473*/
                        sump1 = 0.0;   /*Z=12474*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12475*/
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=12476*/
                        }   /*Z=12477*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump1;   /*Z=12478*/

                        /* cross-sectional P(q) */   /*Z=12480*/
                        CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.,1.0*n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);   /*Z=12481*/
                        sump = 0.0;   /*Z=12482*/
                        sump1 = 0.0;   /*Z=12483*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12484*/
                            sumi = (m+1/2.0)/((m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=12485*/
                            sump = sump+pn[n-m]*sumi;   /*Z=12486*/
                            sump1 = sump1+sumi;   /*Z=12487*/
                        }   /*Z=12488*/
                        CR->carr5p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=12489*/
                        CR->carr6p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrn[n-1]*sump1;   /*Z=12490*/
                        sump = 0.0;   /*Z=12491*/
                        sump1 = 0.0;   /*Z=12492*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12493*/
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-alfa/2.0)*(m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);   /*Z=12494*/
                            sump = sump+sumi;   /*Z=12495*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=12496*/
                        }   /*Z=12497*/
                        CR->carr7p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=12498*/
                        CR->carr8p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump1;   /*Z=12499*/
                        CR->carr9p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrn[n-1]*sump;   /*Z=12500*/

                        /* F(q) */   /*Z=12521*/
                        CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);   /*Z=12522*/
                        CR->carr5f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=12523*/
                        CR->carr6f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrn[n-1]*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=12524*/

                        /* series for <...> integration */   /*Z=12526*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=12527*/
                        sump = 0.0;   /*Z=12528*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12529*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=12530*/
                        }   /*Z=12531*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=12532*/

                        if ( search1 )
                        {   /*Z=12534*/
                            if ( fabs(CR->carr1p[n])<1e-50 || fabs(CR->carr1p[n])>1e+50 ) //NEU
                            {   /*Z=12535*/
                                n1 = n;   /*Z=12536*/
                                search1 = false;   /*Z=12537*/
                            }   /*Z=12538*/
                        }   /*Z=12539*/
                        if ( search4 )
                        {   /*Z=12540*/
                            if ( fabs(CR->carr4p[n])<1e-50 || fabs(CR->carr4p[n])>1e+50 ) //NEU
                            {   /*Z=12541*/
                                n4 = n;   /*Z=12542*/
                                search4 = false;   /*Z=12543*/
                            }   /*Z=12544*/
                        }   /*Z=12545*/
                        if ( search5 )
                        {   /*Z=12546*/
                            if ( fabs(CR->carr5p[n])<1e-50 || fabs(CR->carr5p[n])>1e+50 ) //NEU
                            {   /*Z=12547*/
                                n5 = n;   /*Z=12548*/
                                search5 = false;   /*Z=12549*/
                            }   /*Z=12550*/
                        }   /*Z=12551*/
                        if ( search6 )
                        {   /*Z=12552*/
                            if ( fabs(CR->carr6p[n])<1e-50 || fabs(CR->carr6p[n])>1e+50 )
                            {   /*Z=12553*/
                                n6 = n;   /*Z=12554*/
                                search6 = false;   /*Z=12555*/
                            }   /*Z=12556*/
                        }   /*Z=12557*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=12558*/
                            if ( n<n7 ) n7 = n;   /*Z=12559*/
                        }   /*Z=12560*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=12561*/
                            if ( n<n8 ) n8 = n;   /*Z=12562*/
                        }   /*Z=12563*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=12564*/
                            if ( n<n9 ) n9 = n;   /*Z=12565*/
                        }   /*Z=12566*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12567*/
                            if ( n<n1f ) n1f = n;   /*Z=12568*/
                        }   /*Z=12569*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12570*/
                            if ( n<n4f ) n4f = n;   /*Z=12571*/
                        }   /*Z=12572*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=12573*/
                            if ( n<n5f ) n5f = n;   /*Z=12574*/
                        }   /*Z=12575*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=12576*/
                            if ( n<n6f ) n6f = n;   /*Z=12577*/
                        }   /*Z=12578*/
                        if ( fabs(CR->carr7f[n])<min )
                        {   /*Z=12579*/
                            if ( n<n7f ) n7f = n;   /*Z=12580*/
                        }   /*Z=12581*/
                        if ( fabs(CR->carr8f[n])<min )
                        {   /*Z=12582*/
                            if ( n<n8f ) n8f = n;   /*Z=12583*/
                        }   /*Z=12584*/
                        if ( fabs(CR->carr9f[n])<min )
                        {   /*Z=12585*/
                            if ( n<n9f ) n9f = n;   /*Z=12586*/
                        }   /*Z=12587*/
                    } /* of n-loop */   /*Z=12588*/
                }  /* of inhomogeneous core/shell */   /*Z=12589*/
            }   /* of disk */   /*Z=12593*/

        }   /* of cho1=1 */  /*Z=12710*/


        /*** general orientation case ***/   /*Z=12712*/
        /* for phi<>0, too slow, only for cylinders */   /*Z=12713*/
        if ( cho1==5 )
        {   /*Z=12714*/
            /*Z=12715*/
            p11 = -cos(phi*M_PI/180.0)*cos(theta*M_PI/180.0);   /*Z=12716*/
            p12 = sin(phi*M_PI/180.0);   /*Z=12717*/
            p13 = cos(phi*M_PI/180.0)*sin(theta*M_PI/180.0);   /*Z=12718*/
            p21 = -cos(phi*M_PI/180.0);   /*Z=12719*/
            p22 = -sin(phi*M_PI/180.0)*cos(theta*M_PI/180.0);   /*Z=12720*/
            p23 = sin(phi*M_PI/180.0)*sin(theta*M_PI/180.0);   /*Z=12721*/
            p31 = -sin(theta*M_PI/180.0);   /*Z=12722*/
            p32 = 0;   /*Z=12723*/
            p33 = -cos(theta*M_PI/180.0);   /*Z=12724*/
            /*Z=12725*/
            for ( n=1; n<=2*nmax; n++ ) fkv[n] = fkv[n-1]*n;   /*Z=12726*/
            /*Z=12727*/
            i = 2;   /*Z=12728*/
            for ( n=1; n<=nmax; n++ )
            {   /*Z=12729*/
                for ( m=0; m<=2*n; m++ )
                {   /*Z=12730*/
                    qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,5,0,m,2*n-m,CR->carr1p,intl);   /*Z=12731*/
                    /* carr1pm[i]:=intl/(fkv[m]*fkv[2*n-m]*norm); */   /*Z=12732*/
                    i = i+1;   /*Z=12733*/
                }   /*Z=12734*/
            }   /*Z=12735*/

            double b1sv_n = 1;
            for ( n=1; n<=nmax; n++ )
            {   /*Z=12737*/
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12738*/
                b1sv_n = b1sv_n*(b1s-1+n);   /*Z=12739*/
                xln[n] = -xln[n-1]*xl2z;   /*Z=12740*/
                xrn[n] = -xrn[n-1]*xr2z;   /*Z=12741*/
                CR->carr1p[n] = pow(4.,2.*n)*z12v[n]*xln[n]/((2*n+1)*(n+1));   /*Z=12742*/
                CR->carr3p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*b1sv_n*fkv[n]);   /*Z=12743*/
                if ( search1 )
                {   /*Z=12744*/
                    if ( fabs(CR->carr1p[n])<1e-50 )
                    {   /*Z=12745*/
                        n1 = n;   /*Z=12746*/
                        search1 = false;   /*Z=12747*/
                    }   /*Z=12748*/
                }   /*Z=12749*/
                if ( search3 )
                {   /*Z=12750*/
                    if ( fabs(CR->carr3p[n])<1e-50 )
                    {   /*Z=12751*/
                        n3 = n;   /*Z=12752*/
                        search3 = false;   /*Z=12753*/
                    }   /*Z=12754*/
                }   /*Z=12755*/
            }   /*Z=12756*/
        }   /* of cho1=5 */   /*Z=12757*/

        /*** x-axis ***/   /*Z=12762*/
        if ( (cho1==2) && (dim!=3) )
        {   /*Z=12763*/
            /*** cylinders ***/   /*Z=12764*/
            if ( dim==1 )
            {   /*Z=12765*/
                /* homogeneous */   /*Z=12766*/
                if ( cs==0 )
                {   /*Z=12767*/
                    i = 2;   /*Z=12768*/

                    DBG( qDebug() << "sigma" << sigma << "sigmal" << sigmal << "b1s" << b1s << "length" << l
                                 << "radius" << r << "xlz" << xlz << "xrz" << xrz; )

                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12769*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12770*/
                        z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12771*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=12772*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12773*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12774*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12775*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12776*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12777*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=12778*/
                        /* longitudinal */   /*Z=12779*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12780*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=12781*/
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm); */   /*Z=12782*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=12783, cs=0,dim=1,cho1=2*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12784*/
                            /* i:=i+1; */   /*Z=12785*/
                        }   /*Z=12786*/
                        /* P(q)-coefficient */   /*Z=12787*/
                        //CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=12788*/
                        CR->carr1p[n] = (1./2.)*pow(4.,1.0*n)*z12vl[n]*xln[n]/((n+1)*(n+1./2.));   /*Mail 20220809*/
                        D1P( qDebug() << "CARR1P(A)" << n << CR->carr1p[n]; )

                        /* F(q)-coefficient */   /*Z=12789*/
                        sump = 0.0;   /*Z=12790*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12791*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=12792*/

                        /* cross-section */   /*Z=12794*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);   /*Z=12795*/
                        //D4P( qDebug() << "CARR4P(A)" << n << CR->carr4p[n]; )
                        sump = 0.0;   /*Z=12796*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=12797*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=12798*/

                        if ( search1 )
                        {   /*Z=12801*/
                            if ( fabs(CR->carr1p[n])<1e-50 || fabs(CR->carr1p[n])>1e+50 )
                            {   /*Z=12802*/
                                n1 = n;   /*Z=12803*/
                                search1 = false;   /*Z=12804*/
                            }   /*Z=12805*/
                        }   /*Z=12806*/
                        if ( search4 )
                        {   /*Z=12807*/
                            if ( fabs(CR->carr4p[n])<1e-50 || fabs(CR->carr4p[n])>1e+50 )
                            {   /*Z=12808*/
                                n4 = n;   /*Z=12809*/
                                search4 = false;   /*Z=12810*/
                            }   /*Z=12811*/
                        }   /*Z=12812*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12813*/
                            if ( n<n1f ) n1f = n;   /*Z=12814*/
                        }   /*Z=12815*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12816*/
                            if ( n<n4f ) n4f = n;   /*Z=12817*/
                        }   /*Z=12818*/
                    }  /* of n-loop */   /*Z=12819*/
                }  /* of cs=0 */   /*Z=12820*/

                /* core/shell */   /*Z=12822*/
                if ( cs==1 )
                {   /*Z=12823*/
                    i = 2;   /*Z=12824*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12825*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12826*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=12827*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=12828*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12829*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12830*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12831*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12832*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12833*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=12834*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=12835*/
                        pn[n] = pn[n-1]*p*p;   /*Z=12836*/
                        /* longitudinal */   /*Z=12837*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12838*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=12839*/
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm); */   /*Z=12840*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=12841, cs=1,dim=1,cho1=2*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12842*/
                            /* i:=i+1; */   /*Z=12843*/
                        }   /*Z=12844*/
                        /* P(q)-coefficient */   /*Z=12845*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]/((2*n+1)*(n+1));   /*Z=12846*/
                        D1P( qDebug() << "CARR1P(B)" << n << CR->carr1p[n]; )
                        /* F(q)-coefficient */   /*Z=12847*/
                        sump = 0.0;   /*Z=12848*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12849*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=12850*/

                        /* P(q)-coefficients */   /*Z=12852*/
                        /* cross-sectional */   /*Z=12853*/
                        /* F121 */   /*Z=12854*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=12855*/
                        //D4P( qDebug() << "CARR4P(B)" << n << CR->carr4p[n]; )
                        /* F122 */   /*Z=12856*/
                        sump = 0.0;   /*Z=12857*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12858*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=12859*/
                        }   /*Z=12860*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=12861*/
                        /* F123 */   /*Z=12862*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=12863*/

                        /* F(q)-coefficients */   /*Z=12865*/
                        /* cross-sectional */   /*Z=12866*/
                        /* F121 */   /*Z=12867*/
                        sump = 0.0;   /*Z=12868*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12869*/
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=12870*/
                        }   /*Z=12871*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=12872*/
                        /* F122 */   /*Z=12873*/
                        sump = 0.0;   /*Z=12874*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12875*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=12876*/
                        }   /*Z=12877*/
                        CR->carr5f[n] = xrmn_n*sump;   /*Z=12878*/
                        /* F123 */   /*Z=12879*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=12880*/

                        if ( search1 )
                        {   /*Z=12883*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=12884*/
                                n1 = n;   /*Z=12885*/
                                search1 = false;   /*Z=12886*/
                            }   /*Z=12887*/
                        }   /*Z=12888*/
                        if ( search4 )
                        {   /*Z=12889*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=12890*/
                                n4 = n;   /*Z=12891*/
                                search4 = false;   /*Z=12892*/
                            }   /*Z=12893*/
                        }   /*Z=12894*/
                        if ( search5 )
                        {   /*Z=12895*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=12896*/
                                n5 = n;   /*Z=12897*/
                                search5 = false;   /*Z=12898*/
                            }   /*Z=12899*/
                        }   /*Z=12900*/
                        if ( search6 )
                        {   /*Z=12901*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=12902*/
                                n6 = n;   /*Z=12903*/
                                search6 = false;   /*Z=12904*/
                            }   /*Z=12905*/
                        }   /*Z=12906*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=12907*/
                            if ( n<n1f ) n1f = n;   /*Z=12908*/
                        }   /*Z=12909*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=12910*/
                            if ( n<n4f ) n4f = n;   /*Z=12911*/
                        }   /*Z=12912*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=12913*/
                            if ( n<n5f ) n5f = n;   /*Z=12914*/
                        }   /*Z=12915*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=12916*/
                            if ( n<n6f ) n6f = n;   /*Z=12917*/
                        }   /*Z=12918*/
                    }  /* of n-loop */   /*Z=12919*/
                }  /* of cs=1 */   /*Z=12920*/

                /* inhomogeneous core/shell */   /*Z=12922*/
                if ( cs==2 )
                {   /*Z=12923*/
                    i = 2;   /*Z=12924*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=12925*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=12926*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=12927*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=12928*/
                        fkv[n] = fkv[n-1]*n;   /*Z=12929*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=12930*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=12931*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=12932*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=12933*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=12934*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=12935*/
                        pn[n] = pn[n-1]*p*p;   /*Z=12936*/
                        /* longitudinal */   /*Z=12937*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12938*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=12939*/
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm); */   /*Z=12940*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=12941, cs=2,dim=1,cho1=2*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=12942*/
                            /* i:=i+1; */   /*Z=12943*/
                        }   /*Z=12944*/
                        /* P(q)-coefficient */   /*Z=12945*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=12946*/
                        D1P( qDebug() << "CARR1P(C)" << n << CR->carr1p[n]; )
                        /* F(q)-coefficient */   /*Z=12947*/
                        sump = 0.0;   /*Z=12948*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=12949*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=12950*/

                        /* cross-sectional P(q) */   /*Z=12953*/
                        CR->carr4p[n] = pow(4.,n+1.)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=12954*/
                        //D4P( qDebug() << "CARR4P(C)" << n << CR->carr4p[n]; )
                        sump = 0.0;   /*Z=12955*/
                        sump1 = 0.0;   /*Z=12956*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12957*/
                            sumi = 1/((m+1-alfa/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);   /*Z=12958*/
                            sump = sump+pn[n-m]*sumi;   /*Z=12959*/
                            sump1 = sump1+sumi;   /*Z=12960*/
                        }   /*Z=12961*/
                        CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;   /*Z=12962*/
                        CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n-1]*sump1;   /*Z=12963*/
                        sump = 0.0;   /*Z=12964*/
                        sump1 = 0.0;   /*Z=12965*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=12966*/
                            sumi = 1/((n-m+1-alfa/2.0)*(m+1-alfa/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);   /*Z=12967*/
                            sump = sump+sumi;   /*Z=12968*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=12969*/
                        }   /*Z=12970*/
                        CR->carr7p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump;   /*Z=12971*/
                        CR->carr8p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump1;   /*Z=12972*/
                        CR->carr9p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrn[n-1]*sump;   /*Z=12973*/

                        /* F(q)-coefficients */   /*Z=12988*/
                        /* cross-sectional */   /*Z=12989*/
                        CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);   /*Z=12990*/
                        CR->carr5f[n] = (1-alfa/2.0)*z12v[n]*xrmn_n/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=12991*/
                        CR->carr6f[n] = (1-alfa/2.0)*z12v[n]*xrn[n]/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=12992*/

                        if ( search1 )
                        {   /*Z=12995*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=12996*/
                                n1 = n;   /*Z=12997*/
                                search1 = false;   /*Z=12998*/
                            }   /*Z=12999*/
                        }   /*Z=13000*/
                        if ( search4 )
                        {   /*Z=13001*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13002*/
                                n4 = n;   /*Z=13003*/
                                search4 = false;   /*Z=13004*/
                            }   /*Z=13005*/
                        }   /*Z=13006*/
                        if ( search5 )
                        {   /*Z=13007*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=13008*/
                                n5 = n;   /*Z=13009*/
                                search5 = false;   /*Z=13010*/
                            }   /*Z=13011*/
                        }   /*Z=13012*/
                        if ( search6 )
                        {   /*Z=13013*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=13014*/
                                n6 = n;   /*Z=13015*/
                                search6 = false;   /*Z=13016*/
                            }   /*Z=13017*/
                        }   /*Z=13018*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=13019*/
                            if ( n<n7 ) n7 = n;   /*Z=13020*/
                        }   /*Z=13021*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=13022*/
                            if ( n<n8 ) n8 = n;   /*Z=13023*/
                        }   /*Z=13024*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=13025*/
                            if ( n<n9 ) n9 = n;   /*Z=13026*/
                        }   /*Z=13027*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13028*/
                            if ( n<n1f ) n1f = n;   /*Z=13029*/
                        }   /*Z=13030*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13031*/
                            if ( n<n4f ) n4f = n;   /*Z=13032*/
                        }   /*Z=13033*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=13034*/
                            if ( n<n5f ) n5f = n;   /*Z=13035*/
                        }   /*Z=13036*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=13037*/
                            if ( n<n6f ) n6f = n;   /*Z=13038*/
                        }   /*Z=13039*/
                        if ( fabs(CR->carr7f[n])<min )
                        {   /*Z=13040*/
                            if ( n<n7f ) n7f = n;   /*Z=13041*/
                        }   /*Z=13042*/
                        if ( fabs(CR->carr8f[n])<min )
                        {   /*Z=13043*/
                            if ( n<n8f ) n8f = n;   /*Z=13044*/
                        }   /*Z=13045*/
                        if ( fabs(CR->carr9f[n])<min )
                        {   /*Z=13046*/
                            if ( n<n9f ) n9f = n;   /*Z=13047*/
                        }   /*Z=13048*/
                    }  /* of n-loop */   /*Z=13049*/
                }  /* of cs=1 */   /*Z=13050*/

                /* myelin */   /*Z=13052*/
                if ( (cs==3) || (cs==4) )
                {   /*Z=13053*/
                    i = 2;   /*Z=13054*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13055*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13056*/
                        z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13057*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=13058*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13059*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13060*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13061*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13062*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13063*/
                        xrn[n] = -xrn[n-1]*x12zm;   /*Z=13064*/
                        /* longitudinal */   /*Z=13065*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13066*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=13067*/
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm); */   /*Z=13068*/
                            CR->carr11pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=13069, cs=3|4,dim=1,cho1=2*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13070*/
                            /* i:=i+1; */   /*Z=13071*/
                        }   /*Z=13072*/
                        /* P(q)-coefficient */   /*Z=13073*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=13074*/
                        /* F(q)-coefficient */   /*Z=13075*/
                        sump = 0.0;   /*Z=13076*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13077*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.,1.0*n));   /*Z=13078*/

                        /* cross-section */   /*Z=13080*/
                        CR->carr4p[n] = z12v[n]*xrn[n];   /*Z=13081*/
                        /*Z=13082*/
                        sump = 0.0;   /*Z=13083*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13084*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=13085*/

                        if ( search1 )
                        {   /*Z=13088*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13089*/
                                n1 = n;   /*Z=13090*/
                                search1 = false;   /*Z=13091*/
                            }   /*Z=13092*/
                        }   /*Z=13093*/
                        if ( search4 )
                        {   /*Z=13094*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13095*/
                                n4 = n;   /*Z=13096*/
                                search4 = false;   /*Z=13097*/
                            }   /*Z=13098*/
                        }   /*Z=13099*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13100*/
                            if ( n<n1f ) n1f = n;   /*Z=13101*/
                        }   /*Z=13102*/
                        if ( abs(CR->carr4f[n])<min )
                        {   /*Z=13103*/
                            if ( n<n4f ) n4f = n;   /*Z=13104*/
                        }   /*Z=13105*/
                    }  /* of n-loop */   /*Z=13106*/
                }  /* of cs=3 */   /*Z=13107*/
            }  /* of cylinders */   /*Z=13108*/

            /*** disks ***/   /*Z=13110*/
            if ( dim==2 )
            {   /*Z=13111*/
                /* homogeneous */   /*Z=13112*/
                if ( cs==0 )
                {   /*Z=13113*/
                    i = 2;   /*Z=13114*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13115*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13116*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13117*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13118*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13119*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13120*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13121*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13122*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13123*/
                        /* longitudinal */   /*Z=13124*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13125*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=13126*/
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm); */   /*Z=13127*/
                            CR->carr22pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=13128*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]); */   /*Z=13129*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(pow(4.,1.0*m)*fkv[m]*fkv[n-m]);   /*Z=13130, cs=0,dim=2,cho1=2*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=13131*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13132*/
                            /* i:=i+1; */   /*Z=13133*/
                        }   /*Z=13134*/
                        /* P(q)-coefficient */   /*Z=13135*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=13136*/

                        /* F(q)-coefficient */   /*Z=13138*/
                        sump = 0.0;   /*Z=13139*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13140*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13141*/
                        }   /*Z=13142*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump;   /*Z=13143*/

                        /* cross-section */   /*Z=13145*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=13146*/
                        sump = 0.0;   /*Z=13147*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13148*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13149*/
                        }   /*Z=13150*/
                        CR->carr4f[n] = M_PI*xln[n]*sump/4.0;   /*Z=13151*/

                        /* series for <...> integration */   /*Z=13153*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=13154*/
                        sump = 0.0;   /*Z=13155*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13156*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13157*/
                        }   /*Z=13158*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=13159*/

                        if ( search1 )
                        {   /*Z=13161*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13162*/
                                n1 = n;   /*Z=13163*/
                                search1 = false;   /*Z=13164*/
                            }   /*Z=13165*/
                        }   /*Z=13166*/
                        if ( search4 )
                        {   /*Z=13167*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13168*/
                                n4 = n;   /*Z=13169*/
                                search4 = false;   /*Z=13170*/
                            }   /*Z=13171*/
                        }   /*Z=13172*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13173*/
                            if ( n<n1f ) n1f = n;   /*Z=13174*/
                        }   /*Z=13175*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13176*/
                            if ( n<n4f ) n4f = n;   /*Z=13177*/
                        }   /*Z=13178*/
                    } /* n-loop */   /*Z=13179*/
                }  /* of cs=0 */   /*Z=13180*/

                /* core/shell */   /*Z=13182*/
                if ( cs==1 )
                {   /*Z=13183*/
                    i = 2;   /*Z=13184*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13185*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13186*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13187*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13188*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13189*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13190*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13191*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13192*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13193*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=13194*/
                        pn[n] = pn[n-1]*p*p;   /*Z=13195*/
                        /* longitudinal */   /*Z=13196*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13197*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=13198*/
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm); */   /*Z=13199*/
                            CR->carr22pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=13200*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]); */   /*Z=13201*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(pow(4.,1.0*m)*fkv[m]*fkv[n-m]);   /*Z=13202, cs=1,dim=2,cho1=2*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=13203*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13204*/
                            /* i:=i+1; */   /*Z=13205*/
                        }   /*Z=13206*/
                        /* P(q)-coefficient */   /*Z=13207*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=13208*/
                        /* F(q)-coefficient */   /*Z=13209*/
                        sump = 0.0;   /*Z=13210*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13211*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13212*/
                        }   /*Z=13213*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump;   /*Z=13214*/

                        /* F121 */   /*Z=13216*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=13217*/
                        /* F122 */   /*Z=13218*/
                        sump = 0.0;   /*Z=13219*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13220*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13221*/
                        }   /*Z=13222*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=13223*/

                        /* F122 */   /*Z=13225*/
                        sump = 0.0;   /*Z=13226*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13227*/
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13228*/
                        }   /*Z=13229*/
                        CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;   /*Z=13230*/

                        /* F123 */   /*Z=13233*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=13234*/
                        /* F121 */   /*Z=13235*/
                        sump = 0.0;   /*Z=13236*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13237*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13238*/
                        }   /*Z=13239*/
                        CR->carr4f[n] = M_PI*xrn[n-1]*sump/4.0;   /*Z=13240*/
                        /* F122 */   /*Z=13241*/
                        sump = 0.0;   /*Z=13242*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13243*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13244*/
                        }   /*Z=13245*/
                        CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;   /*Z=13246*/
                        /* F123 */   /*Z=13247*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=13248*/

                        /* series for <...> integration */   /*Z=13250*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=13251*/
                        sump = 0.0;   /*Z=13252*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13253*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13254*/
                        }   /*Z=13255*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=13256*/

                        if ( search1 )
                        {   /*Z=13258*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13259*/
                                n1 = n;   /*Z=13260*/
                                search1 = false;   /*Z=13261*/
                            }   /*Z=13262*/
                        }   /*Z=13263*/
                        if ( search4 )
                        {   /*Z=13264*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13265*/
                                n4 = n;   /*Z=13266*/
                                search4 = false;   /*Z=13267*/
                            }   /*Z=13268*/
                        }   /*Z=13269*/
                        if ( search5 )
                        {   /*Z=13270*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=13271*/
                                n5 = n;   /*Z=13272*/
                                search5 = false;   /*Z=13273*/
                            }   /*Z=13274*/
                        }   /*Z=13275*/
                        if ( search6 )
                        {   /*Z=13276*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=13277*/
                                n6 = n;   /*Z=13278*/
                                search6 = false;   /*Z=13279*/
                            }   /*Z=13280*/
                        }   /*Z=13281*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13282*/
                            if ( n<n1f ) n1f = n;   /*Z=13283*/
                        }   /*Z=13284*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13285*/
                            if ( n<n4f ) n4f = n;   /*Z=13286*/
                        }   /*Z=13287*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=13288*/
                            if ( n<n5f ) n5f = n;   /*Z=13289*/
                        }   /*Z=13290*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=13291*/
                            if ( n<n6f ) n6f = n;   /*Z=13292*/
                        }   /*Z=13293*/
                    } /* n-loop */   /*Z=13294*/
                }  /* of cs=1 */   /*Z=13295*/

                /* inhomogeneous core/shell */   /*Z=13297*/
                if ( cs==2 )
                {   /*Z=13298*/
                    i = 2;   /*Z=13299*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13300*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13301*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13302*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13303*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13304*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13305*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13306*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13307*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13308*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=13309*/
                        pn[n] = pn[n-1]*p*p;   /*Z=13310*/
                        /* longitudinal */   /*Z=13311*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13312*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,CR->carr1p,intl);   /*Z=13313*/
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm); */   /*Z=13314*/
                            CR->carr22pm[n][m] = pow(4.,1.0*m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);   /*Z=13315*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]); */   /*Z=13316*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(pow(4.,1.0*m)*fkv[m]*fkv[n-m]);   /*Z=13317, cs=2,dim=2,cho1=2*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=13318*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13319*/
                            /* i:=i+1; */   /*Z=13320*/
                        }   /*Z=13321*/
                        /* P(q)-coefficient */   /*Z=13322*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=13323*/
                        /* F(q)-coefficient */   /*Z=13324*/
                        sump = 0.0;   /*Z=13325*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13326*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13327*/
                        }   /*Z=13328*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump;   /*Z=13329*/

                        /* cross-sectional P(q) */   /*Z=13331*/
                        CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.,1.0*n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);   /*Z=13332*/
                        sump = 0.0;   /*Z=13333*/
                        sump1 = 0.0;   /*Z=13334*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13335*/
                            sumi = (m+1/2.0)/((m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=13336*/
                            sump = sump+pn[n-m]*sumi;   /*Z=13337*/
                            sump1 = sump1+sumi;   /*Z=13338*/
                        }   /*Z=13339*/
                        CR->carr5p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=13340*/
                        CR->carr6p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrn[n-1]*sump1;   /*Z=13341*/
                        sump = 0.0;   /*Z=13342*/
                        sump1 = 0.0;   /*Z=13343*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13344*/
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-alfa/2.0)*(m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);   /*Z=13345*/
                            sump = sump+sumi;   /*Z=13346*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=13347*/
                        }   /*Z=13348*/
                        CR->carr7p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=13349*/
                        CR->carr8p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump1;   /*Z=13350*/
                        CR->carr9p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrn[n-1]*sump;   /*Z=13351*/

                        /* F(q)  */   /*Z=13371*/
                        CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);   /*Z=13372*/
                        CR->carr5f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=13373*/
                        CR->carr6f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrn[n-1]*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=13374*/

                        /* series for <...> integration */   /*Z=13376*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=13377*/
                        sump = 0.0;   /*Z=13378*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13379*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13380*/
                        }   /*Z=13381*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=13382*/

                        if ( search1 )
                        {   /*Z=13384*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13385*/
                                n1 = n;   /*Z=13386*/
                                search1 = false;   /*Z=13387*/
                            }   /*Z=13388*/
                        }   /*Z=13389*/
                        if ( search4 )
                        {   /*Z=13390*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13391*/
                                n4 = n;   /*Z=13392*/
                                search4 = false;   /*Z=13393*/
                            }   /*Z=13394*/
                        }   /*Z=13395*/
                        if ( search5 )
                        {   /*Z=13396*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=13397*/
                                n5 = n;   /*Z=13398*/
                                search5 = false;   /*Z=13399*/
                            }   /*Z=13400*/
                        }   /*Z=13401*/
                        if ( search6 )
                        {   /*Z=13402*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=13403*/
                                n6 = n;   /*Z=13404*/
                                search6 = false;   /*Z=13405*/
                            }   /*Z=13406*/
                        }   /*Z=13407*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=13408*/
                            if ( n<n7 ) n7 = n;   /*Z=13409*/
                        }   /*Z=13410*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=13411*/
                            if ( n<n8 ) n8 = n;   /*Z=13412*/
                        }   /*Z=13413*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=13414*/
                            if ( n<n9 ) n9 = n;   /*Z=13415*/
                        }   /*Z=13416*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13417*/
                            if ( n<n1f ) n1f = n;   /*Z=13418*/
                        }   /*Z=13419*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13420*/
                            if ( n<n4f ) n4f = n;   /*Z=13421*/
                        }   /*Z=13422*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=13423*/
                            if ( n<n5f ) n5f = n;   /*Z=13424*/
                        }   /*Z=13425*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=13426*/
                            if ( n<n6f ) n6f = n;   /*Z=13427*/
                        }   /*Z=13428*/
                        if ( fabs(CR->carr7f[n])<min )
                        {   /*Z=13429*/
                            if ( n<n7f ) n7f = n;   /*Z=13430*/
                        }   /*Z=13431*/
                        if ( fabs(CR->carr8f[n])<min )
                        {   /*Z=13432*/
                            if ( n<n8f ) n8f = n;   /*Z=13433*/
                        }   /*Z=13434*/
                        if ( fabs(CR->carr9f[n])<min )
                        {   /*Z=13435*/
                            if ( n<n9f ) n9f = n;   /*Z=13436*/
                        }   /*Z=13437*/
                    } /* n-loop */   /*Z=13438*/
                }  /* of inhomogeneous core/shell */   /*Z=13439*/
            }  /* of disk */   /*Z=13441*/
        } // if ( (cho1==2) && (dim!=3) )    /*Z=13442*/



        /*** y-axis ***/   /*Z=13446*/
        if ( (cho1==3) && (dim!=3) )
        {   /*Z=13447*/

            /*** cylinder ***/   /*Z=13449*/
            if ( dim==1 )
            {   /*Z=13450*/
                /* homogeneous */   /*Z=13451*/
                if ( cs==0 )
                {   /*Z=13452*/
                    i = 2;   /*Z=13453*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13454*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13455*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13456*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=13457*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13458*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13459*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13460*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13461*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13462*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13463*/
                        /* longitudinal */   /*Z=13464*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13465*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,CR->carr1p,intl);   /*Z=13466*/
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm); */   /*Z=13467*/
                            CR->carr11pm[n][m] = intl/(pow(4.,1.0*m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);   /*Z=13468, cs=0,dim=1,cho1=3*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13469*/
                            /* i:=i+1; */   /*Z=13470*/
                        }   /*Z=13471*/
                        /* P(q)-coefficient */   /*Z=13472*/
                        CR->carr1p[n] = pow(4.,2.*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=13473*/
                        D1P( qDebug() << "CARR1P(D)" << n << CR->carr1p[n]; )
                        /* F(q)-coefficient */   /*Z=13474*/
                        sump = 0.0;   /*Z=13475*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13476*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump;   /*Z=13477*/
                        //if ( n<2 ) std::cerr << "CARR11PM - C " << norm << std::endl;

                        /* cross-section */   /*Z=13479*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);   /*Z=13480*/
                        //CR->carr4p[n] = 4*(n+1/2.0)*(fk2v[n]/fkv[n])*(z12v[n]*xrn[n])/((n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]);   /*Z=13480*/
                        //D4P( std::cerr << "CARR4P(D) " << n<<" " << CR->carr4p[n]<<" " << fk2v[n]<<" " << z12v[n]<<" " << xrn[n]<<" " << fkv[n] << std::endl; )
                        sump = 0.0;   /*Z=13481*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13482*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=13483*/

                        if ( search1 )
                        {   /*Z=13485*/
                            if ( fabs(CR->carr1p[n])<1e-50 )    // fabs
                            {   /*Z=13486*/
                                n1 = n;   /*Z=13487*/
                                search1 = false;   /*Z=13488*/
                            }   /*Z=13489*/
                        }   /*Z=13490*/
                        if ( search4 )
                        {   /*Z=13491*/
                            if ( fabs(CR->carr4p[n])<1e-50 )    // fabs
                            {   /*Z=13492*/
                                n4 = n;   /*Z=13493*/
                                search4 = false;   /*Z=13494*/
                            }   /*Z=13495*/
                        }   /*Z=13496*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13497*/
                            if ( n<n1f ) n1f = n;   /*Z=13498*/
                        }   /*Z=13499*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13500*/
                            if ( n<n4f ) n4f = n;   /*Z=13501*/
                        }   /*Z=13502*/
                    }  /* of n-loop */   /*Z=13503*/
                }  /* cs=0 */   /*Z=13504*/

                /* core/shell */   /*Z=13506*/
                if ( cs==1 )
                {   /*Z=13507*/
                    i = 2;   /*Z=13508*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13509*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13510*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13511*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=13512*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13513*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13514*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13515*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13516*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13517*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13518*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=13519*/
                        pn[n] = pn[n-1]*p*p;   /*Z=13520*/
                        /* longitudinal */   /*Z=13521*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13522*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,CR->carr1p,intl);   /*Z=13523*/
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm); */   /*Z=13524*/
                            CR->carr11pm[n][m] = intl/(pow(4.,1.0*m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);   /*Z=13525, cs=1,dim=1,cho1=3*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13526*/
                            /* i:=i+1; */   /*Z=13527*/
                        }   /*Z=13528*/
                        /* P(q)-coefficient */   /*Z=13529*/
                        CR->carr1p[n] = pow(4.,2.*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=13530*/
                        /* F(q)-coefficient */   /*Z=13531*/
                        sump = 0.0;   /*Z=13532*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13533*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump;   /*Z=13534*/

                        /* P(q)-coefficients */   /*Z=13536*/
                        /* cross-sectional */   /*Z=13537*/
                        /* F121 */   /*Z=13538*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=13539*/
                        /* F122 */   /*Z=13540*/
                        sump = 0.0;   /*Z=13541*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13542*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13543*/
                        }   /*Z=13544*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=13545*/
                        /* F123 */   /*Z=13546*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=13547*/

                        /* F(q)-coefficients */   /*Z=13549*/
                        /* cross-sectional */   /*Z=13550*/
                        /* F121 */   /*Z=13551*/
                        sump = 0.0;   /*Z=13552*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13553*/
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13554*/
                        }   /*Z=13555*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=13556*/
                        /* F122 */   /*Z=13557*/
                        sump = 0.0;   /*Z=13558*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13559*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13560*/
                        }   /*Z=13561*/
                        CR->carr5f[n] = xrmn_n*sump;   /*Z=13562*/
                        /* F123 */   /*Z=13563*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=13564*/

                        if ( search1 )
                        {   /*Z=13566*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13567*/
                                n1 = n;   /*Z=13568*/
                                search1 = false;   /*Z=13569*/
                            }   /*Z=13570*/
                        }   /*Z=13571*/
                        if ( search4 )
                        {   /*Z=13572*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13573*/
                                n4 = n;   /*Z=13574*/
                                search4 = false;   /*Z=13575*/
                            }   /*Z=13576*/
                        }   /*Z=13577*/
                        if ( search5 )
                        {   /*Z=13578*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=13579*/
                                n5 = n;   /*Z=13580*/
                                search5 = false;   /*Z=13581*/
                            }   /*Z=13582*/
                        }   /*Z=13583*/
                        if ( search6 )
                        {   /*Z=13584*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=13585*/
                                n6 = n;   /*Z=13586*/
                                search6 = false;   /*Z=13587*/
                            }   /*Z=13588*/
                        }   /*Z=13589*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13590*/
                            if ( n<n1f ) n1f = n;   /*Z=13591*/
                        }   /*Z=13592*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13593*/
                            if ( n<n4f ) n4f = n;   /*Z=13594*/
                        }   /*Z=13595*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=13596*/
                            if ( n<n5f ) n5f = n;   /*Z=13597*/
                        }   /*Z=13598*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=13599*/
                            if ( n<n6f ) n6f = n;   /*Z=13600*/
                        }   /*Z=13601*/
                    }  /* of n-loop */   /*Z=13602*/
                }  /* cs=1 */   /*Z=13603*/

                /* inhomogeneous core/shell */   /*Z=13605*/
                if ( cs==2 )
                {   /*Z=13606*/
                    i = 2;   /*Z=13607*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13608*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13609*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13610*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=13611*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13612*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13613*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13614*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13615*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13616*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13617*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=13618*/
                        pn[n] = pn[n-1]*p*p;   /*Z=13619*/
                        /* longitudinal */   /*Z=13620*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13621*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,CR->carr1p,intl);   /*Z=13622*/
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm); */   /*Z=13623*/
                            CR->carr11pm[n][m] = intl/(pow(4.,1.0*m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);   /*Z=13624, cs=2,dim=1,cho1=3*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13625*/
                            /* i:=i+1; */   /*Z=13626*/
                        }   /*Z=13627*/
                        /* P(q)-coefficient */   /*Z=13628*/
                        CR->carr1p[n] = pow(4.,2.*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=13629*/
                        /* F(q)-coefficient */   /*Z=13630*/
                        sump = 0.0;   /*Z=13631*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13632*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump;   /*Z=13633*/

                        /* cross-sectional P(q) */   /*Z=13635*/
                        CR->carr4p[n] = pow(4.,n+1.)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=13636*/
                        sump = 0.0;   /*Z=13637*/
                        sump1 = 0.0;   /*Z=13638*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13639*/
                            sumi = 1/((m+1-alfa/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);   /*Z=13640*/
                            sump = sump+pn[n-m]*sumi;   /*Z=13641*/
                            sump1 = sump1+sumi;   /*Z=13642*/
                        }   /*Z=13643*/
                        CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;   /*Z=13644*/
                        CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n-1]*sump1;   /*Z=13645*/
                        sump = 0.0;   /*Z=13646*/
                        sump1 = 0.0;   /*Z=13647*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13648*/
                            sumi = 1/((n-m+1-alfa/2.0)*(m+1-alfa/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);   /*Z=13649*/
                            sump = sump+sumi;   /*Z=13650*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=13651*/
                        }   /*Z=13652*/
                        CR->carr7p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump;   /*Z=13653*/
                        CR->carr8p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump1;   /*Z=13654*/
                        CR->carr9p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrn[n-1]*sump;   /*Z=13655*/

                        /* F(q)-coefficients */   /*Z=13670*/
                        /* cross-sectional */   /*Z=13671*/
                        CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);   /*Z=13672*/
                        CR->carr5f[n] = (1-alfa/2.0)*z12v[n]*xrmn_n/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=13673*/
                        CR->carr6f[n] = (1-alfa/2.0)*z12v[n]*xrn[n]/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=13674*/

                        if ( search1 )
                        {   /*Z=13676*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13677*/
                                n1 = n;   /*Z=13678*/
                                search1 = false;   /*Z=13679*/
                            }   /*Z=13680*/
                        }   /*Z=13681*/
                        if ( search4 )
                        {   /*Z=13682*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13683*/
                                n4 = n;   /*Z=13684*/
                                search4 = false;   /*Z=13685*/
                            }   /*Z=13686*/
                        }   /*Z=13687*/
                        if ( search5 )
                        {   /*Z=13688*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=13689*/
                                n5 = n;   /*Z=13690*/
                                search5 = false;   /*Z=13691*/
                            }   /*Z=13692*/
                        }   /*Z=13693*/
                        if ( search6 )
                        {   /*Z=13694*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=13695*/
                                n6 = n;   /*Z=13696*/
                                search6 = false;   /*Z=13697*/
                            }   /*Z=13698*/
                        }   /*Z=13699*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=13700*/
                            if ( n<n7 ) n7 = n;   /*Z=13701*/
                        }   /*Z=13702*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=13703*/
                            if ( n<n8 ) n8 = n;   /*Z=13704*/
                        }   /*Z=13705*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=13706*/
                            if ( n<n9 ) n9 = n;   /*Z=13707*/
                        }   /*Z=13708*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13709*/
                            if ( n<n1f ) n1f = n;   /*Z=13710*/
                        }   /*Z=13711*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13712*/
                            if ( n<n4f ) n4f = n;   /*Z=13713*/
                        }   /*Z=13714*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=13715*/
                            if ( n<n5f ) n5f = n;   /*Z=13716*/
                        }   /*Z=13717*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=13718*/
                            if ( n<n6f ) n6f = n;   /*Z=13719*/
                        }   /*Z=13720*/
                        if ( fabs(CR->carr7f[n])<min )
                        {   /*Z=13721*/
                            if ( n<n7f ) n7f = n;   /*Z=13722*/
                        }   /*Z=13723*/
                        if ( fabs(CR->carr8f[n])<min )
                        {   /*Z=13724*/
                            if ( n<n8f ) n8f = n;   /*Z=13725*/
                        }   /*Z=13726*/
                        if ( fabs(CR->carr9f[n])<min )
                        {   /*Z=13727*/
                            if ( n<n9f ) n9f = n;   /*Z=13728*/
                        }   /*Z=13729*/
                    }  /* of n-loop */   /*Z=13730*/
                }  /* cs=2 */   /*Z=13731*/

                /* myelin */   /*Z=13734*/
                if ( (cs==3) || (cs==4) )
                {   /*Z=13735*/
                    i = 2;   /*Z=13736*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13737*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13738*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13739*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=13740*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13741*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13742*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13743*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13744*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13745*/
                        xrn[n] = -xrn[n-1]*x12zm;   /*Z=13746*/
                        /* longitudinal */   /*Z=13747*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13748*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,CR->carr1p,intl);   /*Z=13749*/
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm); */   /*Z=13750*/
                            CR->carr11pm[n][m] = intl/(pow(4.,1.0*m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);   /*Z=13751, cs=3|4,dim=1,cho1=3*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13752*/
                            /* i:=i+1; */   /*Z=13753*/
                        }   /*Z=13754*/
                        /* P(q)-coefficient */   /*Z=13755*/
                        CR->carr1p[n] = pow(4.,2.*n)*z12vl[n]*xln[n]/((2.*n+1)*(n+1));   /*Z=13756*/
                        /* F(q)-coefficient */   /*Z=13757*/
                        sump = 0.0;   /*Z=13758*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13759*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump;   /*Z=13760*/

                        /* cross-section */   /*Z=13762*/
                        CR->carr4p[n] = z12v[n]*xrn[n];   /*Z=13763*/

                        sump = 0.0;   /*Z=13765*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13766*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=13767*/

                        if ( search1 )
                        {   /*Z=13769*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13770*/
                                n1 = n;   /*Z=13771*/
                                search1 = false;   /*Z=13772*/
                            }   /*Z=13773*/
                        }   /*Z=13774*/
                        if ( search4 )
                        {   /*Z=13775*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13776*/
                                n4 = n;   /*Z=13777*/
                                search4 = false;   /*Z=13778*/
                            }   /*Z=13779*/
                        }   /*Z=13780*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13781*/
                            if ( n<n1f ) n1f = n;   /*Z=13782*/
                        }   /*Z=13783*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13784*/
                            if ( n<n4f ) n4f = n;   /*Z=13785*/
                        }   /*Z=13786*/
                    }  /* of n-loop */   /*Z=13787*/
                }  /* cs=3 */   /*Z=13788*/

            }  /* of cylinders */   /*Z=13790*/

            /*** disks ***/   /*Z=13793*/
            if ( dim==2 )
            {   /*Z=13794*/
                /* homogeneous */   /*Z=13795*/
                if ( cs==0 )
                {   /*Z=13796*/
                    i = 2;   /*Z=13797*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13798*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13799*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13800*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13801*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13802*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13803*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13804*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13805*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13806*/
                        /* longitudinal */   /*Z=13807*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13808*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,CR->carr1p,intl);   /*Z=13809*/
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm); */   /*Z=13810*/
                            CR->carr22pm[n][m] = intl/(pow(4.,1.0*m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);   /*Z=13811*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]); */   /*Z=13812*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(fkv[m]*fkv[n-m]);   /*Z=13813, cs=0,dim=2,cho1=3*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=13814*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13815*/
                            /* i:=i+1; */   /*Z=13816*/
                        }   /*Z=13817*/
                        /* P(q)-coefficient */   /*Z=13818*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=13819*/
                        /* F(q)-coefficient */   /*Z=13820*/
                        sump = 0.0;   /*Z=13821*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13822*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13823*/
                        }   /*Z=13824*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump;   /*Z=13825*/

                        /* cross-section */   /*Z=13827*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=13828*/
                        sump = 0.0;   /*Z=13829*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13830*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13831*/
                        }   /*Z=13832*/
                        CR->carr4f[n] = M_PI*xln[n]*sump/4.0;   /*Z=13833*/

                        /* series for <...> integration */   /*Z=13835*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=13836*/
                        sump = 0.0;   /*Z=13837*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13838*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13839*/
                        }   /*Z=13840*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=13841*/

                        if ( search1 )
                        {   /*Z=13843*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13844*/
                                n1 = n;   /*Z=13845*/
                                search1 = false;   /*Z=13846*/
                            }   /*Z=13847*/
                        }   /*Z=13848*/
                        if ( search4 )
                        {   /*Z=13849*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13850*/
                                n4 = n;   /*Z=13851*/
                                search4 = false;   /*Z=13852*/
                            }   /*Z=13853*/
                        }   /*Z=13854*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13855*/
                            if ( n<n1f ) n1f = n;   /*Z=13856*/
                        }   /*Z=13857*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13858*/
                            if ( n<n4f ) n4f = n;   /*Z=13859*/
                        }   /*Z=13860*/
                    }  /* of n-loop */   /*Z=13861*/
                }  /* of cs=0 */   /*Z=13862*/

                /* core/shell */   /*Z=13864*/
                if ( cs==1 )
                {   /*Z=13865*/
                    i = 2;   /*Z=13866*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13867*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13868*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13869*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13870*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13871*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13872*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13873*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13874*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13875*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=13876*/
                        pn[n] = pn[n-1]*p*p;   /*Z=13877*/
                        /* longitudinal */   /*Z=13878*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13879*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,CR->carr1p,intl);   /*Z=13880*/
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm); */   /*Z=13881*/
                            CR->carr22pm[n][m] = intl/(pow(4.,1.0*m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);   /*Z=13882*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]); */   /*Z=13883*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(fkv[m]*fkv[n-m]);   /*Z=13884, cs=1,dim=2,cho1=3*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=13885*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13886*/
                            /* i:=i+1; */   /*Z=13887*/
                        }   /*Z=13888*/
                        /* P(q)-coefficient */   /*Z=13889*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=13890*/
                        /* F(q)-coefficient */   /*Z=13891*/
                        sump = 0.0;   /*Z=13892*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13893*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=13894*/
                        }   /*Z=13895*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump;   /*Z=13896*/

                        /* cross-sectional */   /*Z=13898*/
                        /* F121 */   /*Z=13899*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=13900*/
                        /* F122 */   /*Z=13901*/
                        sump = 0.0;   /*Z=13902*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13903*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13904*/
                        }   /*Z=13905*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=13906*/
                        /* F122 */   /*Z=13907*/
                        sump = 0.0;   /*Z=13908*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13909*/
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13910*/
                        }   /*Z=13911*/
                        CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;   /*Z=13912*/
                        /* F123 */   /*Z=13913*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=13914*/
                        /* F121 */   /*Z=13915*/
                        sump = 0.0;   /*Z=13916*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13917*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13918*/
                        }   /*Z=13919*/
                        CR->carr4f[n] = M_PI*xrn[n-1]*sump/4.0;   /*Z=13920*/
                        /* F122 */   /*Z=13921*/
                        sump = 0.0;   /*Z=13922*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13923*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=13924*/
                        }   /*Z=13925*/
                        CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;   /*Z=13926*/
                        /* F123 */   /*Z=13927*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=13928*/

                        /* series for <...> integration */   /*Z=13930*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=13931*/
                        sump = 0.0;   /*Z=13932*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13933*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=13934*/
                        }   /*Z=13935*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=13936*/

                        if ( search1 )
                        {   /*Z=13938*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=13939*/
                                n1 = n;   /*Z=13940*/
                                search1 = false;   /*Z=13941*/
                            }   /*Z=13942*/
                        }   /*Z=13943*/
                        if ( search4 )
                        {   /*Z=13944*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=13945*/
                                n4 = n;   /*Z=13946*/
                                search4 = false;   /*Z=13947*/
                            }   /*Z=13948*/
                        }   /*Z=13949*/
                        if ( search5 )
                        {   /*Z=13950*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=13951*/
                                n5 = n;   /*Z=13952*/
                                search5 = false;   /*Z=13953*/
                            }   /*Z=13954*/
                        }   /*Z=13955*/
                        if ( search6 )
                        {   /*Z=13956*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=13957*/
                                n6 = n;   /*Z=13958*/
                                search6 = false;   /*Z=13959*/
                            }   /*Z=13960*/
                        }   /*Z=13961*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=13962*/
                            if ( n<n1f ) n1f = n;   /*Z=13963*/
                        }   /*Z=13964*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=13965*/
                            if ( n<n4f ) n4f = n;   /*Z=13966*/
                        }   /*Z=13967*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=13968*/
                            if ( n<n5f ) n5f = n;   /*Z=13969*/
                        }   /*Z=13970*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=13971*/
                            if ( n<n6f ) n6f = n;   /*Z=13972*/
                        }   /*Z=13973*/
                    }  /* of n-loop */   /*Z=13974*/
                }  /* of cs=1 */   /*Z=13975*/

                /* inhomogeneous core/shell */   /*Z=13977*/
                if ( cs==2 )
                {   /*Z=13978*/
                    i = 2;   /*Z=13979*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=13980*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=13981*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=13982*/
                        fkv[n] = fkv[n-1]*n;   /*Z=13983*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=13984*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=13985*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=13986*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=13987*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=13988*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=13989*/
                        pn[n] = pn[n-1]*p*p;   /*Z=13990*/
                        /* longitudinal */   /*Z=13991*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=13992*/
                            qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,CR->carr1p,intl);   /*Z=13993*/
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm); */   /*Z=13994*/
                            CR->carr22pm[n][m] = intl/(pow(4.,1.0*m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);   /*Z=13995*/
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]); */   /*Z=13996*/
                            CR->carr11pm[n][m] = pow(-1.,1.0*m)*fk2v[m]/(fkv[m]*fkv[n-m]);   /*Z=13997, cs=2,dim=2,cho1=3*/
                            /* carr2fm[i]:=carr2pm[i]; */   /*Z=13998*/
                            /* carr1fm[i]:=carr1pm[i]; */   /*Z=13999*/
                            /* i:=i+1; */   /*Z=14000*/
                        }   /*Z=14001*/
                        /* P(q)-coefficient */   /*Z=14002*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=14003*/
                        /* F(q)-coefficient */   /*Z=14004*/
                        sump = 0.0;   /*Z=14005*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14006*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=14007*/
                        }   /*Z=14008*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump;   /*Z=14009*/

                        /* cross-sectional P(q) */   /*Z=14011*/
                        CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.,1.0*n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);   /*Z=14012*/
                        sump = 0.0;   /*Z=14013*/
                        sump1 = 0.0;   /*Z=14014*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14015*/
                            sumi = (m+1/2.0)/((m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=14016*/
                            sump = sump+pn[n-m]*sumi;   /*Z=14017*/
                            sump1 = sump1+sumi;   /*Z=14018*/
                        }   /*Z=14019*/
                        CR->carr5p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=14020*/
                        CR->carr6p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrn[n-1]*sump1;   /*Z=14021*/
                        sump = 0.0;   /*Z=14022*/
                        sump1 = 0.0;   /*Z=14023*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14024*/
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-alfa/2.0)*(m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);   /*Z=14025*/
                            sump = sump+sumi;   /*Z=14026*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=14027*/
                        }   /*Z=14028*/
                        CR->carr7p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=14029*/
                        CR->carr8p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump1;   /*Z=14030*/
                        CR->carr9p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrn[n-1]*sump;   /*Z=14031*/

                        /* F(q) */   /*Z=14050*/
                        CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);   /*Z=14051*/
                        CR->carr5f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=14052*/
                        CR->carr6f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrn[n-1]*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=14053*/

                        /* series for <...> integration */   /*Z=14055*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=14056*/
                        sump = 0.0;   /*Z=14057*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14058*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=14059*/
                        }   /*Z=14060*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=14061*/

                        if ( search1 )
                        {   /*Z=14063*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14064*/
                                n1 = n;   /*Z=14065*/
                                search1 = false;   /*Z=14066*/
                            }   /*Z=14067*/
                        }   /*Z=14068*/
                        if ( search4 )
                        {   /*Z=14069*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14070*/
                                n4 = n;   /*Z=14071*/
                                search4 = false;   /*Z=14072*/
                            }   /*Z=14073*/
                        }   /*Z=14074*/
                        if ( search5 )
                        {   /*Z=14075*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=14076*/
                                n5 = n;   /*Z=14077*/
                                search5 = false;   /*Z=14078*/
                            }   /*Z=14079*/
                        }   /*Z=14080*/
                        if ( search6 )
                        {   /*Z=14081*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=14082*/
                                n6 = n;   /*Z=14083*/
                                search6 = false;   /*Z=14084*/
                            }   /*Z=14085*/
                        }   /*Z=14086*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=14087*/
                            if ( n<n7 ) n7 = n;   /*Z=14088*/
                        }   /*Z=14089*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=14090*/
                            if ( n<n8 ) n8 = n;   /*Z=14091*/
                        }   /*Z=14092*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=14093*/
                            if ( n<n9 ) n9 = n;   /*Z=14094*/
                        }   /*Z=14095*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14096*/
                            if ( n<n1f ) n1f = n;   /*Z=14097*/
                        }   /*Z=14098*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14099*/
                            if ( n<n4f ) n4f = n;   /*Z=14100*/
                        }   /*Z=14101*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=14102*/
                            if ( n<n5f ) n5f = n;   /*Z=14103*/
                        }   /*Z=14104*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=14105*/
                            if ( n<n6f ) n6f = n;   /*Z=14106*/
                        }   /*Z=14107*/
                    }  /* of n-loop */   /*Z=14108*/
                }  /* of cs=2 */   /*Z=14109*/

            }  /* of disks */   /*Z=14111*/
        }   /*Z=14112*/

        /*** z-axis ***/   /*Z=14114*/
        if ( (cho1==4) && (dim!=3) )
        {   /*Z=14115*/
            /*** cylinders ***/   /*Z=14116*/
            if ( dim==1 )
            {   /*Z=14117*/
                /* homogeneous */   /*Z=14118*/
                if ( cs==0 )
                {   /*Z=14119*/
                    i = 2;   /*Z=14120*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14121*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=14122*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=14123*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=14124*/
                        fkv[n] = fkv[n-1]*n;   /*Z=14125*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=14126*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=14127*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=14128*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=14129*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=14130*/
                        /* longitudinal */   /*Z=14131*/
                        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,CR->carr1p,intl);   /*Z=14132*/
                        /* P(q)-coefficient */   /*Z=14133*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]*intl/((2*n+1)*(n+1)*fkv[n]*fkv[n]*norm);   /*Z=14134*/
                        /* F(q)-coefficient */   /*Z=14135*/
                        sump = 0.0;   /*Z=14136*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14137*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.,1.0*n)*fkv[n]*fkv[n]*norm);   /*Z=14138*/
                        /* cross-sectional */   /*Z=14139*/
                        /* P(q)-coefficient */   /*Z=14140*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);   /*Z=14141*/
                        /* F(q)-coefficient */   /*Z=14142*/
                        sump = 0.0;   /*Z=14143*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=14144*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=14145*/

                        if ( search1 )
                        {   /*Z=14147*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14148*/
                                n1 = n;   /*Z=14149*/
                                search1 = false;   /*Z=14150*/
                            }   /*Z=14151*/
                        }   /*Z=14152*/
                        if ( search4 )
                        {   /*Z=14153*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14154*/
                                n4 = n;   /*Z=14155*/
                                search4 = false;   /*Z=14156*/
                            }   /*Z=14157*/
                        }   /*Z=14158*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14159*/
                            if ( n<n1f ) n1f = n;   /*Z=14160*/
                        }   /*Z=14161*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14162*/
                            if ( n<n4f ) n4f = n;   /*Z=14163*/
                        }   /*Z=14164*/
                    }  /* of n-loop */   /*Z=14165*/
                }  /* of cs=0 */   /*Z=14166*/

                /* core/shell */   /*Z=14168*/
                if ( cs==1 )
                {   /*Z=14169*/
                    i = 2;   /*Z=14170*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14171*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=14172*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=14173*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=14174*/
                        fkv[n] = fkv[n-1]*n;   /*Z=14175*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=14176*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=14177*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=14178*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=14179*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=14180*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=14181*/
                        pn[n] = pn[n-1]*p*p;   /*Z=14182*/
                        /* longitudinal */   /*Z=14183*/
                        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,CR->carr1p,intl);   /*Z=14184*/
                        /* P(q)-coefficient */   /*Z=14185*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]*intl/((2*n+1)*(n+1)*fkv[n]*fkv[n]*norm);   /*Z=14186*/
                        /* F(q)-coefficient */   /*Z=14187*/
                        sump = 0.0;   /*Z=14188*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14189*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.,n+1.)*fkv[n]*fkv[n]*norm);   /*Z=14190*/
                        /* P(q)-coefficients */   /*Z=14191*/
                        /* cross-sectional */   /*Z=14192*/
                        /* F121 */   /*Z=14193*/
                        CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=14194*/
                        /* F122 */   /*Z=14195*/
                        sump = 0.0;   /*Z=14196*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14197*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=14198*/
                        }   /*Z=14199*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=14200*/
                        /* F123 */   /*Z=14201*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=14202*/

                        /* F(q)-coefficients */   /*Z=14204*/
                        /* cross-sectional */   /*Z=14205*/
                        /* F121 */   /*Z=14206*/
                        sump = 0.0;   /*Z=14207*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14208*/
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=14209*/
                        }   /*Z=14210*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=14211*/
                        /* F122 */   /*Z=14212*/
                        sump = 0.0;   /*Z=14213*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14214*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=14215*/
                        }   /*Z=14216*/
                        CR->carr5f[n] = xrmn_n*sump;   /*Z=14217*/
                        /* F123 */   /*Z=14218*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=14219*/

                        if ( search1 )
                        {   /*Z=14221*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14222*/
                                n1 = n;   /*Z=14223*/
                                search1 = false;   /*Z=14224*/
                            }   /*Z=14225*/
                        }   /*Z=14226*/
                        if ( search4 )
                        {   /*Z=14227*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14228*/
                                n4 = n;   /*Z=14229*/
                                search4 = false;   /*Z=14230*/
                            }   /*Z=14231*/
                        }   /*Z=14232*/
                        if ( search5 )
                        {   /*Z=14233*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=14234*/
                                n5 = n;   /*Z=14235*/
                                search5 = false;   /*Z=14236*/
                            }   /*Z=14237*/
                        }   /*Z=14238*/
                        if ( search6 )
                        {   /*Z=14239*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=14240*/
                                n6 = n;   /*Z=14241*/
                                search6 = false;   /*Z=14242*/
                            }   /*Z=14243*/
                        }   /*Z=14244*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14245*/
                            if ( n<n1f ) n1f = n;   /*Z=14246*/
                        }   /*Z=14247*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14248*/
                            if ( n<n4f ) n4f = n;   /*Z=14249*/
                        }   /*Z=14250*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=14251*/
                            if ( n<n5f ) n5f = n;   /*Z=14252*/
                        }   /*Z=14253*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=14254*/
                            if ( n<n6f ) n6f = n;   /*Z=14255*/
                        }   /*Z=14256*/
                    }  /* of n-loop */   /*Z=14257*/
                }  /* of cs=1 */   /*Z=14258*/

                /* inhomogeneous core/shell */   /*Z=14261*/
                if ( cs==2 )
                {   /*Z=14262*/
                    i = 2;   /*Z=14263*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14264*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=14265*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=14266*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=14267*/
                        fkv[n] = fkv[n-1]*n;   /*Z=14268*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=14269*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=14270*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=14271*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=14272*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=14273*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=14274*/
                        pn[n] = pn[n-1]*p*p;   /*Z=14275*/
                        /* longitudinal */   /*Z=14276*/
                        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,CR->carr1p,intl);   /*Z=14277*/
                        /* P(q)-coefficient */   /*Z=14278*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]*intl/((2*n+1)*(n+1)*fkv[n]*fkv[n]*norm);   /*Z=14279*/
                        /* F(q)-coefficient */   /*Z=14280*/
                        sump = 0.0;   /*Z=14281*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14282*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.,n+1.)*fkv[n]*fkv[n]*norm);   /*Z=14283*/

                        /* cross-sectional P(q) */   /*Z=14285*/
                        CR->carr4p[n] = pow(4.,n+1.)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*Z=14286*/
                        sump = 0.0;   /*Z=14287*/
                        sump1 = 0.0;   /*Z=14288*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14289*/
                            sumi = 1/((m+1-alfa/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);   /*Z=14290*/
                            sump = sump+pn[n-m]*sumi;   /*Z=14291*/
                            sump1 = sump1+sumi;   /*Z=14292*/
                        }   /*Z=14293*/
                        CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;   /*Z=14294*/
                        CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn[n-1]*sump1;   /*Z=14295*/
                        sump = 0.0;   /*Z=14296*/
                        sump1 = 0.0;   /*Z=14297*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14298*/
                            sumi = 1/((n-m+1-alfa/2.0)*(m+1-alfa/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);   /*Z=14299*/
                            sump = sump+sumi;   /*Z=14300*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=14301*/
                        }   /*Z=14302*/
                        CR->carr7p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump;   /*Z=14303*/
                        CR->carr8p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrmn_n*sump1;   /*Z=14304*/
                        CR->carr9p[n] = (1-alfa/2.0)*(1-alfa/2.0)*z12v[n]*xrn[n-1]*sump;   /*Z=14305*/

                        /* F(q)-coefficients */   /*Z=14320*/
                        /* cross-sectional */   /*Z=14321*/
                        CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);   /*Z=14322*/
                        CR->carr5f[n] = (1-alfa/2.0)*z12v[n]*xrmn_n/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=14323*/
                        CR->carr6f[n] = (1-alfa/2.0)*z12v[n]*xrn[n]/((n+1-alfa/2.0)*fkv[n]*fkv[n]);   /*Z=14324*/

                        if ( search1 )
                        {   /*Z=14326*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14327*/
                                n1 = n;   /*Z=14328*/
                                search1 = false;   /*Z=14329*/
                            }   /*Z=14330*/
                        }   /*Z=14331*/
                        if ( search4 )
                        {   /*Z=14332*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14333*/
                                n4 = n;   /*Z=14334*/
                                search4 = false;   /*Z=14335*/
                            }   /*Z=14336*/
                        }   /*Z=14337*/
                        if ( search5 )
                        {   /*Z=14338*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=14339*/
                                n5 = n;   /*Z=14340*/
                                search5 = false;   /*Z=14341*/
                            }   /*Z=14342*/
                        }   /*Z=14343*/
                        if ( search6 )
                        {   /*Z=14344*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=14345*/
                                n6 = n;   /*Z=14346*/
                                search6 = false;   /*Z=14347*/
                            }   /*Z=14348*/
                        }   /*Z=14349*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=14350*/
                            if ( n<n7 ) n7 = n;   /*Z=14351*/
                        }   /*Z=14352*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=14353*/
                            if ( n<n8 ) n8 = n;   /*Z=14354*/
                        }   /*Z=14355*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=14356*/
                            if ( n<n9 ) n9 = n;   /*Z=14357*/
                        }   /*Z=14358*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14359*/
                            if ( n<n1f ) n1f = n;   /*Z=14360*/
                        }   /*Z=14361*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14362*/
                            if ( n<n4f ) n4f = n;   /*Z=14363*/
                        }   /*Z=14364*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=14365*/
                            if ( n<n5f ) n5f = n;   /*Z=14366*/
                        }   /*Z=14367*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=14368*/
                            if ( n<n6f ) n6f = n;   /*Z=14369*/
                        }   /*Z=14370*/
                    }  /* of n-loop */   /*Z=14371*/
                }  /* of cs=2 */   /*Z=14372*/

                /* myelin */   /*Z=14374*/
                if ( (cs==3) || (cs==4) )
                {   /*Z=14375*/
                    i = 2;   /*Z=14376*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14377*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=14378*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=14379*/
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=14380*/
                        fkv[n] = fkv[n-1]*n;   /*Z=14381*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=14382*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=14383*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=14384*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=14385*/
                        xrn[n] = -xrn[n-1]*x12zm;   /*Z=14386*/
                        /* longitudinal */   /*Z=14387*/
                        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,CR->carr1p,intl);   /*Z=14388*/
                        /* P(q)-coefficient */   /*Z=14389*/
                        CR->carr1p[n] = pow(4.,1.0*n)*z12vl[n]*xln[n]*intl/((2*n+1)*(n+1)*fkv[n]*fkv[n]*norm);   /*Z=14390*/
                        /* F(q)-coefficient */   /*Z=14391*/
                        sump = 0.0;   /*Z=14392*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14393*/
                        CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.,1.0*n)*fkv[n]*fkv[n]*norm);   /*Z=14394*/
                        /* cross-sectional */   /*Z=14395*/
                        /* P(q)-coefficient */   /*Z=14396*/
                        CR->carr4p[n] = z12v[n]*xrn[n];   /*Z=14397*/

                        /* F(q)-coefficient */   /*Z=14399*/
                        sump = 0.0;   /*Z=14400*/
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=14401*/
                        CR->carr4f[n] = xrn[n-1]*sump;   /*Z=14402*/

                        if ( search1 )
                        {   /*Z=14404*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14405*/
                                n1 = n;   /*Z=14406*/
                                search1 = false;   /*Z=14407*/
                            }   /*Z=14408*/
                        }   /*Z=14409*/
                        if ( search4 )
                        {   /*Z=14410*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14411*/
                                n4 = n;   /*Z=14412*/
                                search4 = false;   /*Z=14413*/
                            }   /*Z=14414*/
                        }   /*Z=14415*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14416*/
                            if ( n<n1f ) n1f = n;   /*Z=14417*/
                        }   /*Z=14418*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14419*/
                            if ( n<n4f ) n4f = n;   /*Z=14420*/
                        }   /*Z=14421*/
                    }  /* of n-loop */   /*Z=14422*/
                }  /* of cs=3 */   /*Z=14423*/
            }  /* of cylinders */   /*Z=14424*/

            /*** disks ***/   /*Z=14426*/
            if ( dim==2 )
            {   /*Z=14427*/
                /* homogeneous */   /*Z=14428*/
                if ( cs==0 )
                {   /*Z=14429*/
                    carr2i[0] = 1;   /*Z=14430*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14431*/
                        fkv[n] = fkv[n-1]*n;   /*Z=14432*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=14433*/
                        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,CR->carr1p,intl);   /*Z=14434*/
                        carr2i[n] = pow(-1.,1.0*n)*fk2v[n]*intl/(pow(4.,1.0*n)*fkv[n]*fkv[n]*fkv[n]*norm);   /*Z=14435*/
                    }   /*Z=14436*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14437*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=14438*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=14439*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=14440*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=14441*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=14442*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=14443*/
                        /* longitudinal */   /*Z=14444*/
                        sump = 0.0;   /*Z=14445*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14446*/
                            sump = sump+carr2i[m]/fkv[n-m];   /*Z=14447*/
                        }   /*Z=14448*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=14449*/
                        sump1 = 0.0;   /*Z=14450*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14451*/
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=14452*/
                        }   /*Z=14453*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;   /*Z=14454*/

                        /* cross-sectional */   /*Z=14456*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=14457*/
                        sump = 0.0;   /*Z=14458*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14459*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14460*/
                        }   /*Z=14461*/
                        CR->carr4f[n] = M_PI*xln[n]*sump/4.0;   /*Z=14462*/

                        /* series for <...> integration */   /*Z=14464*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=14465*/
                        sump = 0.0;   /*Z=14466*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14467*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=14468*/
                        }   /*Z=14469*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=14470*/

                        if ( search1 )
                        {   /*Z=14472*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14473*/
                                n1 = n;   /*Z=14474*/
                                search1 = false;   /*Z=14475*/
                            }   /*Z=14476*/
                        }   /*Z=14477*/
                        if ( search4 )
                        {   /*Z=14478*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14479*/
                                n4 = n;   /*Z=14480*/
                                search4 = false;   /*Z=14481*/
                            }   /*Z=14482*/
                        }   /*Z=14483*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14484*/
                            if ( n<n1f ) n1f = n;   /*Z=14485*/
                        }   /*Z=14486*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14487*/
                            if ( n<n4f ) n4f = n;   /*Z=14488*/
                        }   /*Z=14489*/
                    }  /* of n-loop */   /*Z=14490*/
                }  /* of cs=0 */   /*Z=14491*/

                /* core/shell */   /*Z=14494*/
                if ( cs==1 )
                {   /*Z=14495*/
                    carr2i[0] = 1;   /*Z=14496*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14497*/
                        fkv[n] = fkv[n-1]*n;   /*Z=14498*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=14499*/
                        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,CR->carr1p,intl);   /*Z=14500*/
                        carr2i[n] = pow(-1.,1.0*n)*fk2v[n]*intl/(pow(4.,1.0*n)*fkv[n]*fkv[n]*fkv[n]*norm);   /*Z=14501*/
                    }   /*Z=14502*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14503*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=14504*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=14505*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=14506*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=14507*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=14508*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=14509*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=14510*/
                        pn[n] = pn[n-1]*p*p;   /*Z=14511*/
                        /* longitudinal */   /*Z=14512*/
                        sump = 0.0;   /*Z=14513*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14514*/
                            sump = sump+carr2i[m]/fkv[n-m];   /*Z=14515*/
                        }   /*Z=14516*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=14517*/
                        sump1 = 0.0;   /*Z=14518*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14519*/
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=14520*/
                        }   /*Z=14521*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;   /*Z=14522*/

                        /* cross-sectional */   /*Z=14524*/
                        /* F121 */   /*Z=14525*/
                        CR->carr4p[n] = pow(4.,1.0*n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);   /*Z=14526*/
                        /* F122 */   /*Z=14527*/
                        sump = 0.0;   /*Z=14528*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14529*/
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=14530*/
                        }   /*Z=14531*/
                        CR->carr5p[n] = z12v[n]*xrmn_n*sump;   /*Z=14532*/
                        /* F122 */   /*Z=14533*/
                        sump = 0.0;   /*Z=14534*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14535*/
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14536*/
                        }   /*Z=14537*/
                        CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;   /*Z=14538*/
                        /* F123 */   /*Z=14539*/
                        CR->carr6p[n] = CR->carr4p[n]/pn[n];   /*Z=14540*/
                        /* F121 */   /*Z=14541*/
                        sump = 0.0;   /*Z=14542*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14543*/
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14544*/
                        }   /*Z=14545*/
                        CR->carr4f[n] = M_PI*xrn[n-1]*sump/4.0;   /*Z=14546*/
                        /* F122 */   /*Z=14547*/
                        sump = 0.0;   /*Z=14548*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14549*/
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);   /*Z=14550*/
                        }   /*Z=14551*/
                        CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;   /*Z=14552*/
                        /* F123 */   /*Z=14553*/
                        CR->carr6f[n] = CR->carr4f[n]/pn[n];   /*Z=14554*/

                        /* series for <...> integration */   /*Z=14556*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=14557*/
                        sump = 0.0;   /*Z=14558*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14559*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=14560*/
                        }   /*Z=14561*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=14562*/

                        if ( search1 )
                        {   /*Z=14564*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14565*/
                                n1 = n;   /*Z=14566*/
                                search1 = false;   /*Z=14567*/
                            }   /*Z=14568*/
                        }   /*Z=14569*/
                        if ( search4 )
                        {   /*Z=14570*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14571*/
                                n4 = n;   /*Z=14572*/
                                search4 = false;   /*Z=14573*/
                            }   /*Z=14574*/
                        }   /*Z=14575*/
                        if ( search5 )
                        {   /*Z=14576*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=14577*/
                                n5 = n;   /*Z=14578*/
                                search5 = false;   /*Z=14579*/
                            }   /*Z=14580*/
                        }   /*Z=14581*/
                        if ( search6 )
                        {   /*Z=14582*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=14583*/
                                n6 = n;   /*Z=14584*/
                                search6 = false;   /*Z=14585*/
                            }   /*Z=14586*/
                        }   /*Z=14587*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14588*/
                            if ( n<n1f ) n1f = n;   /*Z=14589*/
                        }   /*Z=14590*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14591*/
                            if ( n<n4f ) n4f = n;   /*Z=14592*/
                        }   /*Z=14593*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=14594*/
                            if ( n<n5f ) n5f = n;   /*Z=14595*/
                        }   /*Z=14596*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=14597*/
                            if ( n<n6f ) n6f = n;   /*Z=14598*/
                        }   /*Z=14599*/
                    }  /* of n-loop */   /*Z=14600*/
                }  /* of cs=1 */   /*Z=14601*/

                /* inhomogeneous core/shell */   /*Z=14603*/
                if ( cs==2 )
                {   /*Z=14604*/
                    carr2i[0] = 1;   /*Z=14605*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14606*/
                        fkv[n] = fkv[n-1]*n;   /*Z=14607*/
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);   /*Z=14608*/
                        qrombdeltac(l,r,/*p1,sigma,alfa,dbeta,*/theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,CR->carr1p,intl);   /*Z=14609*/
                        carr2i[n] = pow(-1.,1.0*n)*fk2v[n]*intl/(pow(4.,1.0*n)*fkv[n]*fkv[n]*fkv[n]*norm);   /*Z=14610*/
                    }   /*Z=14611*/
                    for ( n=1; n<=nmax; n++ )
                    {   /*Z=14612*/
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);   /*Z=14613*/
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);   /*Z=14614*/
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z=14615*/
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1); */   /*Z=14616*/
                        xln[n] = -xln[n-1]*xl2z;   /*Z=14617*/
                        xrn[n] = -xrn[n-1]*xr2z;   /*Z=14618*/
                        xrmn_n = -xrmn_n*xrm2z;   /*Z=14619*/
                        pn[n] = pn[n-1]*p*p;   /*Z=14620*/
                        /* longitudinal */   /*Z=14621*/
                        sump = 0.0;   /*Z=14622*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14623*/
                            sump = sump+carr2i[m]/fkv[n-m];   /*Z=14624*/
                        }   /*Z=14625*/
                        CR->carr1p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);   /*Z=14626*/
                        sump1 = 0.0;   /*Z=14627*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14628*/
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);   /*Z=14629*/
                        }   /*Z=14630*/
                        CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;   /*Z=14631*/

                        /* cross-sectional P(q) */   /*Z=14633*/
                        CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.,1.0*n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);   /*Z=14634*/
                        sump = 0.0;   /*Z=14635*/
                        sump1 = 0.0;   /*Z=14636*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14637*/
                            sumi = (m+1/2.0)/((m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);   /*Z=14638*/
                            sump = sump+pn[n-m]*sumi;   /*Z=14639*/
                            sump1 = sump1+sumi;   /*Z=14640*/
                        }   /*Z=14641*/
                        CR->carr5p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=14642*/
                        CR->carr6p[n] = (M_PI/4.0)*(1-alfa)*z12v[n]*xrn[n-1]*sump1;   /*Z=14643*/
                        sump = 0.0;   /*Z=14644*/
                        sump1 = 0.0;   /*Z=14645*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14646*/
                            sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-alfa/2.0)*(m+1/2.0-alfa/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);   /*Z=14647*/
                            sump = sump+sumi;   /*Z=14648*/
                            sump1 = sump1+pn[n-m]*sumi;   /*Z=14649*/
                        }   /*Z=14650*/
                        CR->carr7p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump;   /*Z=14651*/
                        CR->carr8p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrmn_n*sump1;   /*Z=14652*/
                        CR->carr9p[n] = (M_PI/4.0)*(1-alfa)*(1-alfa)*z12v[n]*xrn[n-1]*sump;   /*Z=14653*/

                        /* F(q) */   /*Z=14673*/
                        CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);   /*Z=14674*/
                        CR->carr5f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=14675*/
                        CR->carr6f[n] = (sqrt(M_PI)*(1-alfa)/2.0)*z12v[n]*xrn[n-1]*(n+1)/((n+1/2.0-alfa/2.0)*gam3[n]*fkv[n]);   /*Z=14676*/

                        /* series for <...> integration */   /*Z=14678*/
                        CR->carr2p[n] = pow(4.,n+1.)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);   /*Z=14679*/
                        sump = 0.0;   /*Z=14680*/
                        for ( m=0; m<=n; m++ )
                        {   /*Z=14681*/
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);   /*Z=14682*/
                        }   /*Z=14683*/
                        CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2*gam3[n]);   /*Z=14684*/

                        if ( search1 )
                        {   /*Z=14686*/
                            if ( fabs(CR->carr1p[n])<1e-50 )
                            {   /*Z=14687*/
                                n1 = n;   /*Z=14688*/
                                search1 = false;   /*Z=14689*/
                            }   /*Z=14690*/
                        }   /*Z=14691*/
                        if ( search4 )
                        {   /*Z=14692*/
                            if ( fabs(CR->carr4p[n])<1e-50 )
                            {   /*Z=14693*/
                                n4 = n;   /*Z=14694*/
                                search4 = false;   /*Z=14695*/
                            }   /*Z=14696*/
                        }   /*Z=14697*/
                        if ( search5 )
                        {   /*Z=14698*/
                            if ( fabs(CR->carr5p[n])<1e-50 )
                            {   /*Z=14699*/
                                n5 = n;   /*Z=14700*/
                                search5 = false;   /*Z=14701*/
                            }   /*Z=14702*/
                        }   /*Z=14703*/
                        if ( search6 )
                        {   /*Z=14704*/
                            if ( fabs(CR->carr6p[n])<1e-50 )
                            {   /*Z=14705*/
                                n6 = n;   /*Z=14706*/
                                search6 = false;   /*Z=14707*/
                            }   /*Z=14708*/
                        }   /*Z=14709*/
                        if ( fabs(CR->carr7p[n])<min )
                        {   /*Z=14710*/
                            if ( n<n7 ) n7 = n;   /*Z=14711*/
                        }   /*Z=14712*/
                        if ( fabs(CR->carr8p[n])<min )
                        {   /*Z=14713*/
                            if ( n<n8 ) n8 = n;   /*Z=14714*/
                        }   /*Z=14715*/
                        if ( fabs(CR->carr9p[n])<min )
                        {   /*Z=14716*/
                            if ( n<n9 ) n9 = n;   /*Z=14717*/
                        }   /*Z=14718*/
                        if ( fabs(CR->carr1f[n])<min )
                        {   /*Z=14719*/
                            if ( n<n1f ) n1f = n;   /*Z=14720*/
                        }   /*Z=14721*/
                        if ( fabs(CR->carr4f[n])<min )
                        {   /*Z=14722*/
                            if ( n<n4f ) n4f = n;   /*Z=14723*/
                        }   /*Z=14724*/
                        if ( fabs(CR->carr5f[n])<min )
                        {   /*Z=14725*/
                            if ( n<n5f ) n5f = n;   /*Z=14726*/
                        }   /*Z=14727*/
                        if ( fabs(CR->carr6f[n])<min )
                        {   /*Z=14728*/
                            if ( n<n6f ) n6f = n;   /*Z=14729*/
                        }   /*Z=14730*/
                        if ( fabs(CR->carr7f[n])<min )
                        {   /*Z=14731*/
                            if ( n<n7f ) n7f = n;   /*Z=14732*/
                        }   /*Z=14733*/
                        if ( fabs(CR->carr8f[n])<min )
                        {   /*Z=14734*/
                            if ( n<n8f ) n8f = n;   /*Z=14735*/
                        }   /*Z=14736*/
                        if ( fabs(CR->carr9f[n])<min )
                        {   /*Z=14737*/
                            if ( n<n9f ) n9f = n;   /*Z=14738*/
                        }   /*Z=14739*/
                    }  /* of n-loop */   /*Z=14740*/
                }  /* of cs=2 */   /*Z=14741*/
            }  /* of disk */   /*Z=14744*/
        }  /* of z-axis */   /*Z=14745*/
    }  /* of ordis=0 */   /*Z=14746*/

Label99:  /*Z=14794*/
    // Es kommt vor, dass die carr??[] beim genutzten Index (n?) den Wert 'inf' haben. Damit
    // lässt sich aber schlecht weiterrechnen. Daher hier die passenden Indizes verringern auf
    // den letzten gültigen Wert.
#ifndef __CUDACC__
    //int n1sav=n1, n2sav=n2, n3sav=n3, n4sav=n4, n5sav=n5, n6sav=n6, n7sav=n7, n8sav=n8, n9sav=n9;
#endif
    while ( n1 > 1 && isinf(CR->carr1p[n1]) ) n1--;
    while ( n2 > 1 && isinf(CR->carr2p[n2]) ) n2--;
    while ( n3 > 1 && isinf(CR->carr3p[n3]) ) n3--;
    while ( n4 > 1 && isinf(CR->carr4p[n4]) ) n4--;
    while ( n5 > 1 && isinf(CR->carr5p[n5]) ) n5--;
    while ( n6 > 1 && isinf(CR->carr6p[n6]) ) n6--;
    while ( n7 > 1 && isinf(CR->carr7p[n7]) ) n7--;
    while ( n8 > 1 && isinf(CR->carr8p[n8]) ) n8--;
    while ( n9 > 1 && isinf(CR->carr9p[n9]) ) n9--;
    while ( n1 > 1 && isinf(CR->carr1f[n1]) ) n1--;
    while ( n2 > 1 && isinf(CR->carr2f[n2]) ) n2--;
    while ( n3 > 1 && isinf(CR->carr3f[n3]) ) n3--;
    while ( n4 > 1 && isinf(CR->carr4f[n4]) ) n4--;
    while ( n5 > 1 && isinf(CR->carr5f[n5]) ) n5--;
    while ( n6 > 1 && isinf(CR->carr6f[n6]) ) n6--;
    while ( n7 > 1 && isinf(CR->carr7f[n7]) ) n7--;
    while ( n8 > 1 && isinf(CR->carr8f[n8]) ) n8--;
    while ( n9 > 1 && isinf(CR->carr9f[n9]) ) n9--;
    // Jetzt folgt erst die lim? Berechnung
    lim1 = pow(fabs(CR->carr1p[n1]),-1/(2.*n1));
    lim2 = pow(fabs(CR->carr2p[n2]),-1/(2.*n2));
    lim3 = pow(fabs(CR->carr3p[n3]),-1/(2.*n3));
    lim4 = pow(fabs(CR->carr4p[n4]),-1/(2.*n4));
    lim5 = pow(fabs(CR->carr5p[n5]),-1/(2.*n5));
    lim6 = pow(fabs(CR->carr6p[n6]),-1/(2.*n6));
    lim7 = pow(fabs(CR->carr7p[n7]),-1/(2.*n7));
    lim8 = pow(fabs(CR->carr8p[n8]),-1/(2.*n8));
    lim9 = pow(fabs(CR->carr9p[n9]),-1/(2.*n9));
    lim1f = pow(fabs(CR->carr1f[n1]),-1/(2.*n1));
    lim2f = pow(fabs(CR->carr2f[n2]),-1/(2.*n2));
    lim3f = pow(fabs(CR->carr3f[n3]),-1/(2.*n3));
    lim4f = pow(fabs(CR->carr4f[n4]),-1/(2.*n4));
    lim5f = pow(fabs(CR->carr5f[n5]),-1/(2.*n5));
    lim6f = pow(fabs(CR->carr6f[n6]),-1/(2.*n6));
    lim7f = pow(fabs(CR->carr7f[n7]),-1/(2.*n7));
    lim8f = pow(fabs(CR->carr8f[n8]),-1/(2.*n8));
    lim9f = pow(fabs(CR->carr9f[n9]),-1/(2.*n9));

#ifndef __CUDACC__
/*
    qDebug() << "Label99:" << n1 << n1sav << "1p" << CR->carr1p[n1] << lim1 << "1f" << CR->carr1f[n1] << lim1f;
    qDebug() << "        " << n2 << n2sav << "2p" << CR->carr2p[n2] << lim2 << "2f" << CR->carr2f[n2] << lim2f;
    qDebug() << "        " << n3 << n3sav << "3p" << CR->carr3p[n3] << lim3 << "3f" << CR->carr3f[n3] << lim3f;
    qDebug() << "        " << n4 << n4sav << "4p" << CR->carr4p[n4] << lim4 << "4f" << CR->carr4f[n4] << lim4f;
    qDebug() << "        " << n5 << n5sav << "5p" << CR->carr5p[n5] << lim5 << "5f" << CR->carr5f[n5] << lim5f;
    qDebug() << "        " << n6 << n6sav << "6p" << CR->carr6p[n6] << lim6 << "6f" << CR->carr6f[n6] << lim6f;
    qDebug() << "        " << n7 << n7sav << "7p" << CR->carr7p[n7] << lim7 << "7f" << CR->carr7f[n7] << lim7f;
    qDebug() << "        " << n8 << n8sav << "8p" << CR->carr8p[n8] << lim8 << "8f" << CR->carr8f[n8] << lim8f;
    qDebug() << "        " << n9 << n9sav << "9p" << CR->carr9p[n9] << lim9 << "9f" << CR->carr9f[n9] << lim9f;
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
                str += QString("; %1").arg(CR->carr11pm[n][m]);
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
            str += QString("; %1").arg(CR->carr1p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr2p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr2p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr3p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr3p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr4p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr4p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr5p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr5p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr6p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr6p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr7p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr7p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr8p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr8p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        str = "carr9p[n]";
        for ( int n=0; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
            str += QString("; %1").arg(CR->carr9p[n]);
        fout.write(qPrintable(str+EOL));
        //---
        fout.close();
    }
    else
        qDebug() << fout.errorString();
*/
#endif
}
