/** Crystal calculation.
  *
  * A class with in interface for all parameters an the result and a calculation function
  * to calculate the image.
  * This was originally written in Pascal and translated into C++ by m.wagener@fz-juelich.de
  * During the translation, the variable and type names are preserved.
  *
  * ONLY CPU routines here!
  */

// Dez.2023 ff.
// Versuch, ob man die vielen Switch-Statements etwas aufbrechen kann.
// Hintergrund: durch einen kleineren Kernel könnte die Berechnungszeit in der GPU gesenkt werden,
//              da u.U. mehr Kernel gleichzeitig laufen können.
// Zum Test werden die Routinen natürlich auch in der CPU-Variante genutzt.
// Zum Vergleich können die neuen Funktionen per Flag zur Laufzeit deaktiviert werden.
//
// Die Switch-Statements zum Auflösen und deren Kennungen in den Routinennamen sind:
//      ComboBoxParticle        partXXXX        * 11
//                               partSphere, partCylinder, partDisk, partVesicle, partCube, partEllips
//                               partTriaxEllips, partSuperEllips, partSuperball, partChain, partKPChain
//      shp                     peakXXXX        *  8
//                               peakLorentzian, peakGauss, peakMod1Lor, peakMod2Lor, peakPseudoVoigt
//                               peakPearsonVII, peakGamma, peakAnisoGauss
//      ltype (ltype)           lattXXXX        * 25 - unklar, ob es Sinn macht
//                               lattLamellae, lattHexPackCyl, lattSquarePackCyl, lattRectCentCyl, lattBCC,
//                               lattFCC, lattHCP, lattSC, lattBCT, lattGyroid, lattOBDD, lattPlumNightmare,
//                               lattNone, lattCPLayers, latt2DHexG, latt2DSquareG, latt1DLamG, lattFd3m,
//                               lattOrthorombic, lattLDQ12, lattPercYevick, lattTeubnerStrey, lattPm3n,
//                               lattP42mnm, lattFdddNetwork
//      ordis                   ordisXXXX       * 14 - unklar, ob es Sinn macht
//                               ordisGauss, ordisExponent, ordisOnsager, ordisMaierSaupe, ordisCutOff,
//                               ordisLaguerre, ordisZDir, ordisIsotropic, ordisMirrorGauss, ordisMirrorExponent,
//                               ordisMirrorOnsager, ordisMirrorMaierSaupe, ordisMirrorCutOff, ordisFiberPattern
//  XXXX = die textuelle Kurz-Bezeichnung.
// Sollte eine noch nicht implementierte Kombination vorkommen, dann wird immer die alte Version
// verwendet (calc_GENERIC / kernel_GENERIC) und eine entsprechende Logmeldung generiert.
//
// Neue Berechnungs-Routinen und die passenden Kernel sind:
//      calc_partSphere_peakGauss_lattFCC_ordisGauss,  kernel_partSphere_peakGauss_lattFCC_ordisGauss
//      calc_partSphere_lattNone_ordisIsotropic,       kernel_partSphere_lattNone_ordisIsotropic
//      calc_partCylinder_lattNone_ordisGauss,         kernel_partCylinder_lattNone_ordisGauss
//      calc_partCylinder_lattNone_ordisZDir,          kernel_partCylinder_lattNone_ordisZDir
//      calc_partCylinder_peakGauss_lattBCT_ordisZDir, kernel_partCylinder_peakGauss_lattBCT_ordisZDir
//      calc_partCube_lattNone_ordisIsotropic,         kernel_partCube_lattNone_ordisIsotropic
//
//      partDisk_lattNone_ordisIsotropic
//      partDisk_peakAnisoGauss_lattFCC_ordisGauss
//      partEllips_lattNone_ordisIsotropic
//      partSphere_lattNone_ordisGauss
//      partSphere_lattNone_ordisZDir
//      partSphere_peakAnisoGauss_lattBCC_ordisZDir
//      partSphere_peakAnisoGauss_lattBCT_ordisFiberPattern
//      partSphere_peakAnisoGauss_lattBCT_ordisIsotropic
//      partSphere_peakAnisoGauss_lattBCT_ordisZDir
//      partSphere_peakAnisoGauss_lattFCC_ordisFiberPattern
//      partSphere_peakAnisoGauss_lattFCC_ordisGauss
//      partSphere_peakAnisoGauss_lattFCC_ordisIsotropic
//      partSphere_peakAnisoGauss_lattFCC_ordisZDir
//      partSphere_peakAnisoGauss_lattPercYevick_ordisIsotropic
//      partSphere_peakGauss_lattBCC_ordisIsotropic
//      partSphere_peakGauss_lattBCT_ordisIsotropic
//      partSphere_peakGauss_lattFCC_ordisIsotropic
//      partSphere_peakGauss_lattSC_ordisIsotropic
//      partTriaxEllips_lattNone_ordisIsotropic
//      partVesicle_lattNone_ordisIsotropic
//
// Auch bei den Library-Funktionen kann optimiert werden:
//      formpq_partSphere,   formfq_partSphere
//      formpq_partCylinder, formfq_partCylinder
//      formpq_partDisk,     formfq_partDisk
//      formpq_partCube,     formfq_partCube(nicht verwendet)
//      formpq_partEllips
//      formpq_partTriaxEllips
//      formpq_partSuperEllips
//      formpq_partSuperball



#include "sc_calc_generic_gpu.h"
#include <stdlib.h>
#include <string.h>

#define DBGSPEC(x) //x
// Diese Debugausgaben sind innerhalb der Schleifen und produzieren viel Output.
//  Werden aber sowohl für die CPU als auch für die GPU genutzt.

#include <iostream>
#include <chrono>
#include <unistd.h>
#include <signal.h>


#define myacos(x) ((x)<-1)?0.0:(((x)>+1)?M_PI:acos((x)))



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
    _endThread = false;
}


bool SasCalc_GENERIC_calculation::prepareCalculation()
{
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

    //partsphere = false;  //Z=20075
    //partcylinder = false;  //Z=20076
    //partdisk = false;  //Z=20077
    //partvesicle = false;  //Z=20078
    //partcube = false;  //Z=20079
    //partellipsoid = false;  //Z=20080
    //parttriellipsoid = false;  //Z=20081
    //partbarrel = false;  //Z=20082
    //partball = false;  //Z=20083
    //partchain = false;  //Z=20084
    //partkpchain = false;  //Z=20085

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

#ifdef undef
    //20240301 - ab hier auskommentiert (im 20240301 ist es auch raus)
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
    //20240301 - bis hier auskommentiert (TODO?)
#endif

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
    critdist = params.alpha; // :=StrToFloat(EditRotAlpha.Text); //20240301 neu

    /* ******************* */  //Z=20554
    /* ** lattice tpyes ** */  //Z=20555
    /* ******************* */  //Z=20556
    //if ( CheckBoxTest.Checked==true ) generic = true; /*  generic  */  //Z=20580 ==> lassen wir jetzt mal fest

    //if ( (ltype==4/*BCC*/) || (ltype==5/*FCC*/) ) twin = true;     /*  twinning  */  //Z=20582
    //230728 Lt. Hr. Förster könnte das Twinning hier noch Probleme bereiten, daher nicht automatisch an machen
    twin = CheckBoxTwinned;  // von außen gesteuert

    if ( ltype==12/*None*/ ) lattice = false;  //Z=20583

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
        //partsphere = true;  //Z=20675
        partdim = 3;  //Z=20676
        cdim = 3;    /*  sphere  */  //Z=21000
        break;
    case cbpartCylinder/*1*/:    //(* cylinder *)  //Z=20678
        params.part = 1;      /*  cylinder  */  //Z=20679
        //partcylinder = true;  //Z=20680
        partdim = 2;  //Z=20681
        cdim = 1;    /*  cylinder  */  //Z=21001
        break;
    case cbpartDisk/*2*/:        //(* disk *)  //Z=20683
        params.part = 2;      /*  disk  */  //Z=20684
        //partdisk = true;  //Z=20685
        partdim = 1;  //Z=20686
        cdim = 2;    /*  disk  */  //Z=21002
        break;
    case cbpartVesicle/*3*/:     //(* vesicle *)  //Z=20688
        params.part = 3;      /*  vesicle  */  //Z=20689
        //partvesicle = true;  //Z=20690
        partdim = 3;  //Z=20691
        cdim = 3;    /*  vesicle  */  //Z=21003
        break;
    case cbpartCube/*4*/:        //(* cube *)  //Z=20693
        params.part = 4;      /*  cube  */  //Z=20694
        //partcube = true;  //Z=20695
        partdim = 4;  //Z=20696
        cdim = 4;    /*  cube  */  //Z=21004
        break;
    case cbpartEllipsoide/*5*/:  //(* ellipsoid *)  //Z=20698
        params.part = 5;      /*  ellipsoid  */  //Z=20699
        //partellipsoid = true;  //Z=20700
        partdim = 5;  //Z=20701
        cdim = 5;    /*  ellipsoid  */  //Z=21005
        break;
    case cbpartTriaxEllips/*6*/: //(* triaxial ellipsoid *)  //Z=20703
        params.part = 6;      /*  triaxial ellipsoid  */  //Z=20704
        //parttriellipsoid = true;  //Z=20705
        partdim = 6;  //Z=20706
        cdim = 6;    /*  triaxial ellipsoid  */  //Z=21006
        break;
    case cbpartSuperEllips/*7*/: //(* super ellipsoid, barrel *)  //Z=20708
        params.part = 7;      /*  super ellipsoid, barrel  */  //Z=20709
        //partbarrel = true;  //Z=20710
        partdim = 7;  //Z=20711
        cdim = 7;    /*  super ellipsoid, barrel  */  //Z=21007
        break;
    case cbpartSuperball/*8*/: //  //Z=20713
        params.part = 8;      /*  superball  */  //Z=20714
        //partball = true;  //Z=20715
        partdim = 8;  //Z=20716
        cdim = 8;    /*  superball  */  //Z=21008
        break;
    case cbpartChain/*9*/:  //Z=20718
        params.part = 9;      /*  excluded volume chain  */  //Z=20719
        //partchain = true;  //Z=20720
        partdim = 9;  //Z=20721
        cdim = 9;    /*  excluded volume chain  */  //Z=21009
        break;
    case cbpartkpchain/*10*/:  //Z=20724
        params.part = 10;      /*  Kratky Porod chain  */  //Z=20725
        //partkpchain = true;  //Z=20726
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
    //case 20:    // Percus-Yevick        ZZ=?  //20240301 die ltype 20 bis 24 raus
    //case 21:    // Teubner-Strey        ZZ=?
    //case 22:    // Pm3n, A15            ZZ=?
    //case 23:    // P42/mnm (sigma)      ZZ=? (neu)
    //case 24:    // Fddd-network         ZZ=? (neu)
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
    {   //Z=21015
        coefficients( cdim, 120, order );   //Z=21016
    }   //Z=21020

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
/*#ifndef __CUDACC__
        qDebug() << "ri1*" << ri11 << ri12 << ri13;
        qDebug() << "ri2*" << ri21 << ri22 << ri23;
        qDebug() << "ri3*" << ri31 << ri32 << ri33;
#endif*/
    }/*3*/  //Z=21113

    /* ** reciprocal space vector list ** */  //Z=21115
    if ( /*generic=true &&*/ lattice )
    {/*3*/  //Z=21116
        /* ************************************************* */  //Z=21117
        /* ** isotropic peak list for allowed reflections ** */  //Z=21118
        /* ************************************************* */  //Z=21119

        double sphno=1, qhkl0;
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

        double g3; //, g3t=0;
        for ( int peakct2=1; peakct2<=peakmax2; peakct2++ )
        {/*4*/  //Z=21150
            h = latpar2(peakct2,1);  //20240301 - Z=20573
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
            //qhklt = sqrt(qxhklt*qxhklt+qyhklt*qyhklt+qzhklt*qzhklt);  //ZN=21380

            setLatpar3(peakct2,7, qxhklt);  //Z=21187
            setLatpar3(peakct2,8, qyhklt);  //Z=21189
            setLatpar3(peakct2,9, qzhklt);  //Z=21191

            g3 = 4*M_PI*M_PI/(2.0*M_PI*qhkl*qhkl);  //Z=21192
            //20240301 raus: if ( twin ) g3t = 4*M_PI*M_PI/(2.0*M_PI*qhklt*qhklt);  //ZN=21390
            x2phihkl = 4*qhkl*qhkl/(M_PI*phiwidth*phiwidth);  //Z=21193
            //20240301 raus: x2phihklt = 4*qhklt*qhklt/(M_PI*phiwidth*phiwidth);  //ZN=21392
            peaknorm2=0; // zur Sicherheit

            switch ( shp )      /*Z=20614*/
            {
            case cbpeakLorentzian: // 1
                peaknorm1 = lorentznorm3(x2phihkl);  //Z=21194
                break;
            case cbpeakGaussian: // 2
                peaknorm1 = gaussnorm3(x2phihkl);  //Z=21195
                //20240301 raus: peaknorm1t = gaussnorm3(x2phihklt);  //ZN=21396
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

            //20240301 raus: setLatpar3(peakct2,14, qhklt);  //ZN=21415
            //20240301 raus: setLatpar3(peakct2,15, g3t);  //ZN=21416
            //20240301 raus: setLatpar3(peakct2,16, peaknorm1t);  //ZN=21417

        }/*4*/  //Z=21294
    }/*3*/  /*  of test or ltype=30  */  //20240301-Z=20692

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
        qrombdeltac( params.p1, params.sigma, params.alpha, params.polTheta, params.polPhi, 1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,
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
    return true;
}


void SasCalc_GENERIC_calculation::cleanup()
{
    if ( params.CR != nullptr )
        std::cerr << "latpar Max Check: " << params.CR->latparMaxCheckCount[0] << ", "
                  << params.CR->latparMaxCheckCount[1] << ", " << params.CR->latparMaxCheckCount[2] << std::endl;

    std::cerr << "GENERIC::cleanup:";
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
 * Parametern in prepareCalculation() berechnet werden und dann erst geändert werden
 * dürfen. Danach geht es in die doCalculation().
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
void SasCalc_GENERIC_calculation::doCalculation( int numThreads, bool bIgnNewSwitch )
{
    numberOfThreads = numThreads;
    bIgnoreNewSwitch = bIgnNewSwitch;

    auto start1 = std::chrono::high_resolution_clock::now();
    if ( !prepareCalculation() ) return;
    // Der Speicher wird erst am Ende von prepareCalculation angelegt
    setXYIntensity( zzmin, iimin, 0.0 );    // reset internal indices

    auto end1 = std::chrono::high_resolution_clock::now();
    auto calc_time1 = std::chrono::duration_cast<std::chrono::duration<float>>(end1-start1);
    higResTimerElapsedPrep = calc_time1.count()*1000.0;

    // Use a high resolution clock to get the calculation time of the GPUs
    auto start2 = std::chrono::high_resolution_clock::now();

    if ( gpuAvailable() && numThreads == 0 )    // GPU abschaltbar
    {
        kernel_selection( zzmax - zzmin, iimax - iimin, false );
    }
    else
    {
        D8( numThreads=1; numberOfThreads=1; )

        _endThread = false;

        if ( numThreads <= 1 )
        {   // Without threading it is simple ...
            for ( int ihex=zzmin; ihex<zzmax; ihex++ )  // (***   z-loop  ***)
            {
                if ( _endThread ) break;
                for ( int i=iimin; i<iimax; i++ )
                {
                    calcCPU_selection( *this, false, ihex, i );
                    if ( _endThread ) break;
                }
                //doIntCalc_GENERIC( ihex++ );
            } // for ihex
        } // numThreads <= 1

        else
        {   // Special work for Multi-Threads...
            if ( threads != nullptr ) { delete threads; delete thread_args; }
            threads     = new pthread_t[numThreads];    // Weil die Anzahl der Threads geändert werden kann
            thread_args = new int[numThreads];
            memset( threads, 0, sizeof(pthread_t)*numThreads ); // for ( int t=0; t<numThreads; t++ ) threads[t] = 0;
            int ihex=zzmin;
            // Es macht keinen Sinn, hier eine doppelte Schleife (ihex,i) zu machen und für
            //  jeden Pixel einen eigenen Thread zu machen. Da ist der Overhead der Threads
            //  größer und macht es ineffizient.
            while ( ihex<zzmax )
            {
                if ( _endThread ) break;

                // If the number of threads is greater than the number of cores, this loop must be optimized
                // so that it will not wait for all threads. This might speed up a little more.
                for ( int t=0; t<numThreads && ihex<zzmax; t++ )
                {
                    thread_args[t] = ihex++;
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
            } // while ( ihex<=zzmax && i<=iimax )
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
    } // if ( ! noGPUavailable ) else {

    //performImageExpansion();
    //inline void performImageExpansion()
    //{
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
    //}

    // Use a high resolution clock to get the calculation time of the GPUs
    auto end2 = std::chrono::high_resolution_clock::now();
    auto calc_time2 = std::chrono::duration_cast<std::chrono::duration<float>>(end2-start2);
    higResTimerElapsedCalc = calc_time2.count()*1000.0;
}


/**
 * @brief endThread
 * Helper function to cancel all calculation threads.
 * The QtCreator marks many undefined symbols because this is an include. It compiles without errors.
 */
void SasCalc_GENERIC_calculation::endThread()
{
    if ( _endThread ) return;   // endThread()
    // Jetzt wirklich die Threads killen
    _endThread = true; // Damit die Threads nicht weitermachen
#ifdef __CUDACC__
    if ( gpuAvailable() && numberOfThreads == 0 )
    {
        //cudaReset(); // TODO
    }
    else
#endif
        if ( threads != nullptr && numberOfThreads > 1 )
        {   // Special work for Multi-Threads...
            std::cerr << "endThread() Canceling ..." << std::endl;
            int s;
            bool doJoin=true;
            for ( int i=0; i<numberOfThreads; i++ )
            {
                std::cerr << "    Thread #" << i << ": ";
                if ( threads != nullptr && threads[i] != 0 )
                {
                    s = pthread_cancel( threads[i] );
                    //pthread_kill( threads[i], SIGTERM );
                    if ( s == 0 )
                    {
                        std::cerr << "Canceled." << std::endl << std::flush;
                        if ( doJoin )
                        {   // Nur beim ersten wird auf das Ende gewartet, die anderen sind dann auch beendet
                            pthread_join( threads[i], nullptr );
                            std::cerr << "               finished." << std::endl << std::flush;
                            doJoin = false;
                        }
                    }
                    else
                        std::cerr << "Cancelerror! " << s << std::endl << std::flush;
                    if ( threads != nullptr ) threads[i] = 0;
                }
                else
                    std::cerr << "not running" << std::endl << std::flush;
            } // for i
        } // if threads != 0
} /* endThread() */



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

    // Use a high resolution clock to get the calculation time of the GPUs
    auto start2 = std::chrono::high_resolution_clock::now();

    if ( gpuAvailable() && numThreads == 0 )    // GPU abschaltbar
    {
        kernel_selection( _fitWidth, _fitHeight, true );
    }
    else
    {
        D8( numThreads=1; numberOfThreads=1; )

        _endThread = false;

        if ( numThreads <= 1 )
        {   // Without threading it is simple ...
            for ( int x=0; x<_fitWidth; x++ )
            {
                if ( _endThread ) break;
                for ( int y=0; y<_fitHeight; y++ )
                {
                    calcCPU_selection( *this, true, x, y );
                    if ( _endThread ) break;
                }
                //doIntFitCalc_GENERIC( x++ );
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
    } // if ( ! noGPUavailable ) else {

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
    for ( int i=inst->iimin; i<inst->iimax; i++ )
    {
        inst->calcCPU_selection( *inst, false, ihex, i );
        if ( inst->_endThread ) break;
    }
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
    for ( int y=0; y<inst->_fitHeight; y++ )
    {
        inst->calcCPU_selection( *inst, true, x, y );
        if ( inst->_endThread ) break;
    }
    return nullptr;
}



void SasCalc_GENERIC_calculation::corotations(double a, double b, double c, double alpha, double beta, double gamma,
                                              double u, double v, double w, double ephi,
                                              bool lat1d, bool lat2d, bool /*lat3d*/,
                                              double &m11, double &m12, double &m13, double &m21, double &m22, double &m23, double &m31, double &m32, double &m33,
                                              double &mtw11, double &mtw12, double &mtw13, double &mtw21, double &mtw22, double &mtw23, double &mtw31, double &mtw32, double &mtw33, double &vvol,
                                              double &nuvwx, double &nuvwy, double &nuvwz, double &uuvwx, double &uuvwy, double &uuvwz, double &vuvwx, double &vuvwy, double &vuvwz,
                                              double &nhklx, double &nhkly, double &nhklz, double &uhklx, double &uhkly, double &uhklz, double &vhklx, double &vhkly, double &vhklz )
{/*1*/

    double ca, cb, cg, /*sa,*/ sb, sg, /*vol,*/ n1, n2, n3, l1, l2, l3, m1, m2, m3;
    double msi11, msi12, msi13, msi21, msi22, msi23, msi31, msi32, msi33, detmsi;
    double ms11, ms12, ms13, ms21, ms22, ms23, ms31, ms32, ms33;
    //double msi11n, msi12n, msi13n, msi21n, msi22n, msi23n, msi31n, msi32n, msi33n; //, detmsin;
    //double m11n, m12n, m13n, m21n, m22n, m23n, m31n, m32n, m33n;
    double mst11, mst12, mst13, mst21, mst22, mst23, mst31, mst32, mst33, detmst;
    double mt11, mt12, mt13, mt21, mt22, mt23, mt31, mt32, mt33;
    float g11, g12, g13, g21, g22, g23, g31, g32, g33, gi11, gi12, gi13, gi21, gi22, gi23, gi31, gi32, gi33, detg;
    float mi11, mi12, mi13, mi21, mi22, mi23, mi31, mi32, mi33, detm;
    double ri11, ri12, ri13, ri21, ri22, ri23, ri31, ri32, ri33, detr;
    float aax, aay, aaz, bax, bay, baz, cax, cay, caz;
    float aaxt, aayt, aazt, baxt, bayt, bazt, caxt, cayt, cazt;
    float aex, aey, aez, bex, bey, bez, cex, cey, cez; //, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
    float /*aexn, aeyn, aezn, bexn, beyn, bezn, cexn, ceyn, cezn,*/ asxn, asyn, aszn, bsxn, bsyn, bszn, csxn, csyn, cszn;
    float aext, aeyt, aezt, bext, beyt, bezt, cext, ceyt, cezt; //, asxt, asyt, aszt, bsxt, bsyt, bszt, csxt, csyt, cszt;
    //float asxd, asyd, aszd, bsxd, bsyd, bszd, csxd, csyd, cszd;
    float epsi, etheta, len, r11, r12, r13, r21, r22, r23, r31, r32, r33;
    float aedx, aedy, aedz, bedx, bedy, bedz, cedx, cedy, cedz, vols;
    float nxyzx, nxyzy, nxyzz, uxyzx, uxyzy, uxyzz, vxyzx, vxyzy, vxyzz, nxyznx, /*nxyzny,*/ nxyznz;
    //double asax, asay, asaz, bsax, bsay, bsaz, csax, csay, csaz, asex, asey, asez, bsex, bsey, bsez, csex, csey, csez;
    double asedx, asedy, asedz, bsedx, bsedy, bsedz, csedx, csedy, csedz;
    double asedxn, asedyn, asedzn, bsedxn, bsedyn, bsedzn, csedxn, csedyn, csedzn;
    double asedxt, asedyt, asedzt, bsedxt, bsedyt, bsedzt, csedxt, csedyt, csedzt;
    double aedxt, aedyt, aedzt, bedxt, bedyt, bedzt, cedxt, cedyt, cedzt;
    double m11DD, m12DD, m13DD, m21DD, m22DD, m23DD; //, m31DD, m32DD, m33DD;
    //double mi11DD, mi12DD, mi13DD, mi21DD, mi22DD, mi23DD, mi31DD, mi32DD, mi33DD;
    double aexDD, aeyDD, /*aezDD,*/ bexDD, beyDD; //, bezDD, cexDD, ceyDD, cezDD;
    double asxDD, asyDD, aszDD, bsxDD, bsyDD, bszDD, csxDD, csyDD, cszDD, area;
    //double aedxDD, aedyDD, aedzDD, bedxDD, bedyDD, bedzDD, cedxDD, cedyDD, cedzDD;
    double asedxDD, asedyDD, asedzDD, bsedxDD, bsedyDD, bsedzDD, csedxDD, csedyDD, csedzDD;
    //double m11DL, m12DL, m13DL, m21DL, m22DL, m23DL, m31DL, m32DL, m33DL;
    //double mi11DL, mi12DL, mi13DL, mi21DL, mi22DL, mi23DL, mi31DL, mi32DL, mi33DL;
    //double aexDL, aeyDL, aezDL, bexDL, beyDL, bezDL, cexDL, ceyDL, cezDL;
    double asxDL, asyDL, aszDL, bsxDL, bsyDL, bszDL, csxDL, csyDL, cszDL, area1;
    //double aedxDL, aedyDL, aedzDL, bedxDL, bedyDL, bedzDL, cedxDL, cedyDL, cedzDL;
    double asedxDL, asedyDL, asedzDL, bsedxDL, bsedyDL, bsedzDL, csedxDL, csedyDL, csedzDL;

    /*  unit cell  */
    alpha = alpha*M_PI/180.0;
    beta = beta*M_PI/180.0;
    gamma = gamma*M_PI/180.0;
    ephi = ephi*M_PI/180.0;
    ca = cos(alpha);
    cb = cos(beta);
    cg = cos(gamma);
    //sa = sin(alpha);
    sb = sin(beta);
    sg = sin(gamma);
    //vol = a*b*c*sqrt(1.0-ca*ca-cb*cb-cg*cg+2*ca*cb*cg);

    n1 = 0;
    n2 = 0;
    n3 = 1;
    l1 = sb;
    l2 = 0;
    l3 = cb;
    m1 = (cg-cb*ca)/sb;
    m2 = sqrt(1.0-m1*m1-ca*ca);
    m3 = ca;

    /*  Msi matrix  */
    msi11 = a*l1;         msi12 = a*l2;       msi13 = a*l3;
    msi21 = b*m1;         msi22 = b*m2;       msi23 = b*m3;
    msi31 = c*n1;         msi32 = c*n2;       msi33 = c*n3;

    /*  new: transposed Msi matrix  */
    //msi11n = a*l1;         msi21n = a*l2;       msi31n = a*l3;
    //msi12n = b*m1;         msi22n = b*m2;       msi32n = b*m3;
    //msi13n = c*n1;         msi23n = c*n2;       msi33n = c*n3;

    /*  new: = M matrix  */
    //m11n = msi11;          m12n = msi12;       m13n = msi13;
    //m21n = msi21;          m22n = msi22;       m32n = msi23;
    //m31n = msi31;          m23n = msi32;       m33n = msi33;

    /*  old: Ms matrix  */
    detmsi = msi11*msi22*msi33+msi12*msi23*msi31+msi13*msi21*msi32-msi31*msi22*msi13-msi32*msi23*msi11-msi33*msi21*msi12;
    ms11 = (msi22*msi33-msi23*msi32)/detmsi;
    ms12 = (msi13*msi32-msi12*msi33)/detmsi;
    ms13 = (msi12*msi23-msi13*msi22)/detmsi;
    ms21 = (msi23*msi31-msi21*msi33)/detmsi;
    ms22 = (msi11*msi33-msi13*msi31)/detmsi;
    ms23 = (msi13*msi21-msi11*msi23)/detmsi;
    ms31 = (msi21*msi32-msi22*msi31)/detmsi;
    ms32 = (msi12*msi31-msi11*msi32)/detmsi;
    ms33 = (msi11*msi22-msi12*msi21)/detmsi;

    /*  old: Ms matrix, transposed  */
    mst11 = ms11;    mst12 = ms21;    mst13 = ms31;
    mst21 = ms12;    mst22 = ms22;    mst23 = ms32;
    mst31 = ms13;    mst32 = ms23;    mst33 = ms33;

    /*  old: M matrix, transposed  */
    detmst = mst11*mst22*mst33+mst12*mst23*mst31+mst13*mst21*mst32-mst31*mst22*mst13-mst32*mst23*mst11-mst33*mst21*mst12;
    mt11 = (mst22*mst33-mst23*mst32)/detmst;
    mt12 = (mst13*mst32-mst12*mst33)/detmst;
    mt13 = (mst12*mst23-mst13*mst22)/detmst;
    mt21 = (mst23*mst31-mst21*mst33)/detmst;
    mt22 = (mst11*mst33-mst13*mst31)/detmst;
    mt23 = (mst13*mst21-mst11*mst23)/detmst;
    mt31 = (mst21*mst32-mst22*mst31)/detmst;
    mt32 = (mst12*mst31-mst11*mst32)/detmst;
    mt33 = (mst11*mst22-mst12*mst21)/detmst;

    /*  old: M matrix  */
    /*  rxyz=M ruvw  */
    m11 = mt11;       m12 = mt21;        m13 = mt31;
    m21 = mt12;       m22 = mt22;        m23 = mt32;
    m31 = mt13;       m32 = mt23;        m33 = mt33;

    m11DD = a;     m12DD = b*cg;    m13DD = 0;
    m21DD = 0;     m22DD = b*sg;    m23DD = 0;
    //m31DD = 0;     m32DD = 0;       m33DD = 1;

    //m11DL = a;     m12DL = 0;    m13DL = 0;
    //m21DL = 0;     m22DL = 1;    m23DL = 0;
    //m31DL = 0;     m32DL = 0;    m33DL = 1;


    /*  Mi inverse matrix  */
    /*  ruvw=Mi rxyz  */
    detm = m11*m22*m33+m12*m23*m31+m13*m21*m32-m31*m22*m13-m32*m23*m11-m33*m21*m12;
    mi11 = (m22*m33-m23*m32)/detm;
    mi12 = (m13*m32-m12*m33)/detm;
    mi13 = (m12*m23-m13*m22)/detm;
    mi21 = (m23*m31-m21*m33)/detm;
    mi22 = (m11*m33-m13*m31)/detm;
    mi23 = (m13*m21-m11*m23)/detm;
    mi31 = (m21*m32-m22*m31)/detm;
    mi32 = (m12*m31-m11*m32)/detm;
    mi33 = (m11*m22-m12*m21)/detm;

    //mi11DL = 1/a;      mi12DL = 0;        mi13DL = 0;
    //mi21DL = 0;        mi22DL = 1;        mi23DL = 0;
    //mi31DL = 0;        mi32DL = 0;        mi33DL = 1;


    /*  base vectors of unit cell  */
    /*  aa, ba, ca in uvw-system  */
    aax = 1;         aay = 0;      aaz = 0;
    bax = 0;         bay = 1;      baz = 0;
    cax = 0;         cay = 0;      caz = 1;

    /*  base vectors of twinned unit cell  */
    aaxt = 2/3.0;      aayt = 2/3.0;     aazt = -1/3.0;
    baxt = -1/3.0;     bayt = 2/3.0;     bazt = 2/3.0;
    caxt = 2/3.0;      cayt = -1/3.0;    cazt = 2/3.0;

    /*  unit vectors in carthesian coordinate system  */
    /*  ae=M aa, be=M ba, ce=M ca in xyz-system  */
    /*  Mok  */
    /*  old  */
    aex = m11*aax+m12*aay+m13*aaz;
    aey = m21*aax+m22*aay+m23*aaz;
    aez = m31*aax+m32*aay+m33*aaz;
    bex = m11*bax+m12*bay+m13*baz;
    bey = m21*bax+m22*bay+m23*baz;
    bez = m31*bax+m32*bay+m33*baz;
    cex = m11*cax+m12*cay+m13*caz;
    cey = m21*cax+m22*cay+m23*caz;
    cez = m31*cax+m32*cay+m33*caz;

    /*  new  */
    //aexn = m11n*aax+m12n*aay+m13n*aaz;
    //aeyn = m21n*aax+m22n*aay+m23n*aaz;
    //aezn = m31n*aax+m32n*aay+m33n*aaz;
    //bexn = m11n*bax+m12n*bay+m13n*baz;
    //beyn = m21n*bax+m22n*bay+m23n*baz;
    //bezn = m31n*bax+m32n*bay+m33n*baz;
    //cexn = m11n*cax+m12n*cay+m13n*caz;
    //ceyn = m21n*cax+m22n*cay+m23n*caz;
    //cezn = m31n*cax+m32n*cay+m33n*caz;

    aexDD = m11DD*aax+m12DD*aay+m13DD*aaz;
    aeyDD = m21DD*aax+m22DD*aay+m23DD*aaz;
    //aezDD = m31DD*aax+m32DD*aay+m33DD*aaz;
    bexDD = m11DD*bax+m12DD*bay+m13DD*baz;
    beyDD = m21DD*bax+m22DD*bay+m23DD*baz;
    //bezDD = m31DD*bax+m32DD*bay+m33DD*baz;
    //cexDD = m11DD*cax+m12DD*cay+m13DD*caz;
    //ceyDD = m21DD*cax+m22DD*cay+m23DD*caz;
    //cezDD = m31DD*cax+m32DD*cay+m33DD*caz;

    //aexDL = m11DL*aax+m12DL*aay+m13DL*aaz;
    //aeyDL = m21DL*aax+m22DL*aay+m23DL*aaz;
    //aezDL = m31DL*aax+m32DL*aay+m33DL*aaz;
    //bexDL = m11DL*bax+m12DL*bay+m13DL*baz;
    //beyDL = m21DL*bax+m22DL*bay+m23DL*baz;
    //bezDL = m31DL*bax+m32DL*bay+m33DL*baz;
    //cexDL = m11DL*cax+m12DL*cay+m13DL*caz;
    //ceyDL = m21DL*cax+m22DL*cay+m23DL*caz;
    //cezDL = m31DL*cax+m32DL*cay+m33DL*caz;

    aext = m11*aaxt+m12*aayt+m13*aazt;
    aeyt = m21*aaxt+m22*aayt+m23*aazt;
    aezt = m31*aaxt+m32*aayt+m33*aazt;
    bext = m11*baxt+m12*bayt+m13*bazt;
    beyt = m21*baxt+m22*bayt+m23*bazt;
    bezt = m31*baxt+m32*bayt+m33*bazt;
    cext = m11*caxt+m12*cayt+m13*cazt;
    ceyt = m21*caxt+m22*cayt+m23*cazt;
    cezt = m31*caxt+m32*cayt+m33*cazt;

    /*  old: reciprocal space vector in carthesian coordinate system  */
    /*  ase, bse, cse in xyz-system  */
    /*  Mok  */
    vvol = aex*(bey*cez-bez*cey)+aey*(bez*cex-bex*cez)+aez*(bex*cey-bey*cex);
    //asx = (bey*cez-bez*cey)/vvol;
    //asy = (bez*cex-bex*cez)/vvol;
    //asz = (bex*cey-bey*cex)/vvol;
    //bsx = (aez*cey-aey*cez)/vvol;
    //bsy = (aex*cez-aez*cex)/vvol;
    //bsz = (aey*cex-aex*cey)/vvol;
    //csx = (aey*bez-aez*bey)/vvol;
    //csy = (aez*bex-aex*bez)/vvol;
    //csz = (aex*bey-aey*bex)/vvol;

    area = a*b*sg;
    asxDD = beyDD/area;
    asyDD = -bexDD/area;
    aszDD = 0;
    bsxDD = -aeyDD/area;
    bsyDD = aexDD/area;
    bszDD = 0;
    csxDD = 0;
    csyDD = 0;
    cszDD = 1;

    area1 = a;
    asxDL = 1/area1;
    asyDL = 0;
    aszDL = 0;
    bsxDL = 0;
    bsyDL = 1;
    bszDL = 0;
    csxDL = 0;
    csyDL = 0;
    cszDL = 1;

    //asxt = (beyt*cezt-bezt*ceyt)/vvol;
    //asyt = (bezt*cext-bext*cezt)/vvol;
    //aszt = (bext*ceyt-beyt*cext)/vvol;
    //bsxt = (aezt*ceyt-aeyt*cezt)/vvol;
    //bsyt = (aext*cezt-aezt*cext)/vvol;
    //bszt = (aeyt*cext-aext*ceyt)/vvol;
    //csxt = (aeyt*bezt-aezt*beyt)/vvol;
    //csyt = (aezt*bext-aext*bezt)/vvol;
    //cszt = (aext*beyt-aeyt*bext)/vvol;

    /*  new: G metric matrix in xyz-coordinates  */
    g11 = aex;     g12 = bex;     g13 = cex;
    g21 = aey;     g22 = bey;     g23 = cey;
    g31 = aez;     g32 = bez;     g33 = cez;

    /*  old G matrix  */
    /* g11:=aex*aex+aey*aey+aez*aez;
       g12:=aex*bex+aey*bey+aez*bez;
       g13:=aex*cex+aey*cey+aez*cez;
       g21:=bex*aex+bey*aey+bez*aez;
       g22:=bex*bex+bey*bey+bez*bez;
       g23:=bex*cex+bey*cey+bez*cez;
       g31:=cex*aex+cey*aey+cez*aez;
       g32:=cex*bex+cey*bey+cez*bez;
       g33:=cex*cex+cey*cey+cez*cez;  */

    /*  Gs inverse metric matrix  */
    detg = g11*g22*g33+g12*g23*g31+g13*g21*g32-g31*g22*g13-g32*g23*g11-g33*g21*g12;
    gi11 = (g22*g33-g23*g32)/detg;
    gi12 = (g13*g32-g12*g33)/detg;
    gi13 = (g12*g23-g13*g22)/detg;
    gi21 = (g23*g31-g21*g33)/detg;
    gi22 = (g11*g33-g13*g31)/detg;
    gi23 = (g13*g21-g11*g23)/detg;
    gi31 = (g21*g32-g22*g31)/detg;
    gi32 = (g12*g31-g11*g32)/detg;
    gi33 = (g11*g22-g12*g21)/detg;

    /*  new: reciprocal space vector in carthesian coordinate system  */
    /*  ase, bse, cse in xyz-system  */
    asxn = gi11;   asyn = gi12;   aszn = gi13;
    bsxn = gi21;   bsyn = gi22;   bszn = gi23;
    csxn = gi31;   csyn = gi32;   cszn = gi33;

    /*  a*,b*,c* reciprocal space vectors  */
    /*  in uvw-space  */
    /* asax:=g11*aax+g12*aay+g13*aaz;
       asay:=g21*aax+g22*aay+g23*aaz;
       asaz:=g31*aax+g32*aay+g33*aaz;
       bsax:=g11*bax+g12*bay+g13*baz;
       bsay:=g21*bax+g22*bay+g23*baz;
       bsaz:=g31*bax+g32*bay+g33*baz;
       csax:=g11*cax+g12*cay+g13*caz;
       csay:=g21*cax+g22*cay+g23*caz;
       csaz:=g31*cax+g32*cay+g33*caz; */

    /*  a*, b*, c* reciprocal space vectors  */
    /*  in xyz-space  */
    /* asex:=m11*asax+m12*asay+m13*asaz;
       asey:=m21*asax+m22*asay+m23*asaz;
       asez:=m31*asax+m32*asay+m33*asaz;
       bsex:=m11*bsax+m12*bsay+m13*bsaz;
       bsey:=m21*bsax+m22*bsay+m23*bsaz;
       bsez:=m31*bsax+m32*bsay+m33*bsaz;
       csex:=m11*csax+m12*csay+m13*csaz;
       csey:=m21*csax+m22*csay+m23*csaz;
       csez:=m31*csax+m32*csay+m33*csaz; */


    /*  nuvw-vector || beam  */
    /*  nuvw in unit cell uvw-system  */
    nuvwx = u;
    nuvwy = v;
    nuvwz = w;

    /*  nhkl-vector  */
    /*  nhkl=G nuvw in unit cell uvw-system  */
    /* nhklx:=g11*nuvwx+g12*nuvwy+g13*nuvwz;
       nhkly:=g21*nuvwx+g22*nuvwy+g23*nuvwz;
       nhklz:=g31*nuvwx+g32*nuvwy+g33*nuvwz; */

    /*  nhkl-vector  */
    /*  nxyz=M nuvw in xyz-system  */
    /*  with (a,b,c)=(2,1,1) and (u,v,w)=(2,1,1) this is a (ua,vb,wc)=(4,1,1)-vector  */
    nxyzx = m11*nuvwx+m12*nuvwy+m13*nuvwz;
    nxyzy = m21*nuvwx+m22*nuvwy+m23*nuvwz;
    nxyzz = m31*nuvwx+m32*nuvwy+m33*nuvwz;

    /*  unit nxyz = director  */
    /*  in xyz-system  */
    len = sqrt(nxyzx*nxyzx+nxyzy*nxyzy+nxyzz*nxyzz);
    nxyznx = nxyzx/len;
    //nxyzny = nxyzy/len;
    nxyznz = nxyzz/len;

    /*  R rotation matrix  */
    /*  in xyz-system  */
    etheta = acos(-nxyznz);
    epsi = asin(nxyznx/sin(etheta));

    r11 = cos(epsi)*cos(ephi)-sin(epsi)*cos(etheta)*sin(ephi);
    r12 = -cos(epsi)*sin(ephi)-sin(epsi)*cos(etheta)*cos(ephi);
    r13 = sin(epsi)*sin(etheta);
    r21 = sin(epsi)*cos(ephi)+cos(epsi)*cos(etheta)*sin(ephi);
    r22 = -sin(epsi)*sin(ephi)+cos(epsi)*cos(etheta)*cos(ephi);
    r23 = -cos(epsi)*sin(etheta);
    r31 = sin(etheta)*sin(ephi);
    r32 = sin(etheta)*cos(ephi);
    r33 = cos(etheta);

    /*  Ri inverse rotation matrix  */
    /*  in xyz-system  */
    detr = r11*r22*r33+r12*r23*r31+r13*r21*r32-r31*r22*r13-r32*r23*r11-r33*r21*r12;
    ri11 = (r22*r33-r23*r32)/detr;
    ri12 = (r13*r32-r12*r33)/detr;
    ri13 = (r12*r23-r13*r22)/detr;
    ri21 = (r23*r31-r21*r33)/detr;
    ri22 = (r11*r33-r13*r31)/detr;
    ri23 = (r13*r21-r11*r23)/detr;
    ri31 = (r21*r32-r22*r31)/detr;
    ri32 = (r12*r31-r11*r32)/detr;
    ri33 = (r11*r22-r12*r21)/detr;

    /*  rotated base vectors a,b,c in carthesian coordinate system  */
    /*  aed=Ri ae, bed=Ri be, ced=Ri ce in xyz-system  */
    /*  needed for calculation of fiber pattern  */
    /*  Mok  */
    aedx = ri11*aex+ri12*aey+ri13*aez;
    aedy = ri21*aex+ri22*aey+ri23*aez;
    aedz = ri31*aex+ri32*aey+ri33*aez;
    bedx = ri11*bex+ri12*bey+ri13*bez;
    bedy = ri21*bex+ri22*bey+ri23*bez;
    bedz = ri31*bex+ri32*bey+ri33*bez;
    cedx = ri11*cex+ri12*cey+ri13*cez;
    cedy = ri21*cex+ri22*cey+ri23*cez;
    cedz = ri31*cex+ri32*cey+ri33*cez;

    //aedxDD = ri11*aexDD+ri12*aeyDD+ri13*aezDD;
    //aedyDD = ri21*aexDD+ri22*aeyDD+ri23*aezDD;
    //aedzDD = ri31*aexDD+ri32*aeyDD+ri33*aezDD;
    //bedxDD = ri11*bexDD+ri12*beyDD+ri13*bezDD;
    //bedyDD = ri21*bexDD+ri22*beyDD+ri23*bezDD;
    //bedzDD = ri31*bexDD+ri32*beyDD+ri33*bezDD;
    //cedxDD = ri11*cexDD+ri12*ceyDD+ri13*cezDD;
    //cedyDD = ri21*cexDD+ri22*ceyDD+ri23*cezDD;
    //cedzDD = ri31*cexDD+ri32*ceyDD+ri33*cezDD;

    //aedxDL = ri11*aexDL+ri12*aeyDL+ri13*aezDL;
    //aedyDL = ri21*aexDL+ri22*aeyDL+ri23*aezDL;
    //aedzDL = ri31*aexDL+ri32*aeyDL+ri33*aezDL;
    //bedxDL = ri11*bexDL+ri12*beyDL+ri13*bezDL;
    //bedyDL = ri21*bexDL+ri22*beyDL+ri23*bezDL;
    //bedzDL = ri31*bexDL+ri32*beyDL+ri33*bezDL;
    //cedxDL = ri11*cexDL+ri12*ceyDL+ri13*cezDL;
    //cedyDL = ri21*cexDL+ri22*ceyDL+ri23*cezDL;
    //cedzDL = ri31*cexDL+ri32*ceyDL+ri33*cezDL;

    aedxt = ri11*aext+ri12*aeyt+ri13*aezt;
    aedyt = ri21*aext+ri22*aeyt+ri23*aezt;
    aedzt = ri31*aext+ri32*aeyt+ri33*aezt;
    bedxt = ri11*bext+ri12*beyt+ri13*bezt;
    bedyt = ri21*bext+ri22*beyt+ri23*bezt;
    bedzt = ri31*bext+ri32*beyt+ri33*bezt;
    cedxt = ri11*cext+ri12*ceyt+ri13*cezt;
    cedyt = ri21*cext+ri22*ceyt+ri23*cezt;
    cedzt = ri31*cext+ri32*ceyt+ri33*cezt;

    /*  rotated reciprocal space vectors a*,b*,c* in carthesian coordinate system  */
    /*  calculated from cross-product of aed,bed,ced  */
    /*  into output matrix  */
    /*  Mok  */
    asedx = (bedy*cedz-bedz*cedy)/vvol;
    asedy = (bedz*cedx-bedx*cedz)/vvol;
    asedz = (bedx*cedy-bedy*cedx)/vvol;
    bsedx = (aedz*cedy-aedy*cedz)/vvol;
    bsedy = (aedx*cedz-aedz*cedx)/vvol;
    bsedz = (aedy*cedx-aedx*cedy)/vvol;
    csedx = (aedy*bedz-aedz*bedy)/vvol;
    csedy = (aedz*bedx-aedx*bedz)/vvol;
    csedz = (aedx*bedy-aedy*bedx)/vvol;

    /*  new: a*_r=R a*  */
    /*  output  */
    asedxn = ri11*asxn+ri12*asyn+ri13*aszn;
    asedyn = ri21*asxn+ri22*asyn+ri23*aszn;
    asedzn = ri31*asxn+ri32*asyn+ri33*aszn;
    bsedxn = ri11*bsxn+ri12*bsyn+ri13*bszn;
    bsedyn = ri21*bsxn+ri22*bsyn+ri23*bszn;
    bsedzn = ri31*bsxn+ri32*bsyn+ri33*bszn;
    csedxn = ri11*csxn+ri12*csyn+ri13*cszn;
    csedyn = ri21*csxn+ri22*csyn+ri23*cszn;
    csedzn = ri31*csxn+ri32*csyn+ri33*cszn;

    asedxDD = ri11*asxDD+ri12*asyDD+ri13*aszDD;
    asedyDD = ri21*asxDD+ri22*asyDD+ri23*aszDD;
    asedzDD = ri31*asxDD+ri32*asyDD+ri33*aszDD;
    bsedxDD = ri11*bsxDD+ri12*bsyDD+ri13*bszDD;
    bsedyDD = ri21*bsxDD+ri22*bsyDD+ri23*bszDD;
    bsedzDD = ri31*bsxDD+ri32*bsyDD+ri33*bszDD;
    csedxDD = ri11*csxDD+ri12*csyDD+ri13*cszDD;
    csedyDD = ri21*csxDD+ri22*csyDD+ri23*cszDD;
    csedzDD = ri31*csxDD+ri32*csyDD+ri33*cszDD;

    asedxDL = ri11*asxDL+ri12*asyDL+ri13*aszDL;
    asedyDL = ri21*asxDL+ri22*asyDL+ri23*aszDL;
    asedzDL = ri31*asxDL+ri32*asyDL+ri33*aszDL;
    bsedxDL = ri11*bsxDL+ri12*bsyDL+ri13*bszDL;
    bsedyDL = ri21*bsxDL+ri22*bsyDL+ri23*bszDL;
    bsedzDL = ri31*bsxDL+ri32*bsyDL+ri33*bszDL;
    csedxDL = ri11*csxDL+ri12*csyDL+ri13*cszDL;
    csedyDL = ri21*csxDL+ri22*csyDL+ri23*cszDL;
    csedzDL = ri31*csxDL+ri32*csyDL+ri33*cszDL;

    asedxt = (bedyt*cedzt-bedzt*cedyt)/vvol;
    asedyt = (bedzt*cedxt-bedxt*cedzt)/vvol;
    asedzt = (bedxt*cedyt-bedyt*cedxt)/vvol;
    bsedxt = (aedzt*cedyt-aedyt*cedzt)/vvol;
    bsedyt = (aedxt*cedzt-aedzt*cedxt)/vvol;
    bsedzt = (aedyt*cedxt-aedxt*cedyt)/vvol;
    csedxt = (aedyt*bedzt-aedzt*bedyt)/vvol;
    csedyt = (aedzt*bedxt-aedxt*bedzt)/vvol;
    csedzt = (aedxt*bedyt-aedyt*bedxt)/vvol;

    /*  a*, b*, c* rotated  */
    /*  in xyz-space  */
    /* asedx:=r11*asex+r12*asey+r13*asez;
       asedy:=r21*asex+r22*asey+r23*asez;
       asedz:=r31*asex+r32*asey+r33*asez;
       bsedx:=r11*bsex+r12*bsey+r13*bsez;
       bsedy:=r21*bsex+r22*bsey+r23*bsez;
       bsedz:=r31*bsex+r32*bsey+r33*bsez;
       csedx:=r11*csex+r12*csey+r13*csez;
       csedy:=r21*csex+r22*csey+r23*csez;
       csedz:=r31*csex+r32*csey+r33*csez; */


    /*  n,u,v axis vectors in Carthesian coordinate system  */
    /*  in xyz-system  */
    nxyzx = 0;     nxyzy = 0;     nxyzz = -1;
    uxyzx = 1;     uxyzy = 0;     uxyzz = 0;
    vxyzx = 0;     vxyzy = 1;     vxyzz = 0;

    /*  rotated n,u,v axis vectors  */
    /*  nd=R n, ud=R u, vd=R v  in xyz-system  */
    double nxyzdx = r11*nxyzx+r12*nxyzy+r13*nxyzz;
    double nxyzdy = r21*nxyzx+r22*nxyzy+r23*nxyzz;
    double nxyzdz = r31*nxyzx+r32*nxyzy+r33*nxyzz;
    double uxyzdx = r11*uxyzx+r12*uxyzy+r13*uxyzz;
    double uxyzdy = r21*uxyzx+r22*uxyzy+r23*uxyzz;
    double uxyzdz = r31*uxyzx+r32*uxyzy+r33*uxyzz;
    double vxyzdx = r11*vxyzx+r12*vxyzy+r13*vxyzz;
    double vxyzdy = r21*vxyzx+r22*vxyzy+r23*vxyzz;
    double vxyzdz = r31*vxyzx+r32*vxyzy+r33*vxyzz;

    /*  rotated n,u,v axis vectors  */
    /*  nuvw=Mi nd, uuvw=Mi ud, vuvw=Mi vd  in unit cell uvw-system  */
    /*  needed to indicate unit cell directions <uvw> in 2D-pattern  */
    /*  output  */
    nuvwx = mi11*nxyzdx+mi12*nxyzdy+mi13*nxyzdz;
    nuvwy = mi21*nxyzdx+mi22*nxyzdy+mi23*nxyzdz;
    nuvwz = mi31*nxyzdx+mi32*nxyzdy+mi33*nxyzdz;
    uuvwx = mi11*uxyzdx+mi12*uxyzdy+mi13*uxyzdz;
    uuvwy = mi21*uxyzdx+mi22*uxyzdy+mi23*uxyzdz;
    uuvwz = mi31*uxyzdx+mi32*uxyzdy+mi33*uxyzdz;
    vuvwx = mi11*vxyzdx+mi12*vxyzdy+mi13*vxyzdz;
    vuvwy = mi21*vxyzdx+mi22*vxyzdy+mi23*vxyzdz;
    vuvwz = mi31*vxyzdx+mi32*vxyzdy+mi33*vxyzdz;

    /*  rotated n,u,v axis vectors  */
    /*  nhkl=G nuvw, uhkl=G uuvw, vhkl=G vuvw  in unit cell uvw-system  */
    /*  needed to indicate reciprocal space directions <hkl> in 2D-pattern  */
    vols = nuvwx*(uuvwy*vuvwz-uuvwz*vuvwy)+nuvwy*(uuvwz*vuvwx-uuvwx*vuvwz)+nuvwz*(uuvwx*vuvwy-uuvwy*vuvwx);
    nhklx = (uuvwy*vuvwz-uuvwz*vuvwy)/vols;
    nhkly = (uuvwz*vuvwx-uuvwx*vuvwz)/vols;
    nhklz = (uuvwx*vuvwy-uuvwy*vuvwx)/vols;
    uhklx = (nuvwz*vuvwy-nuvwy*vuvwz)/vols;
    uhkly = (nuvwx*vuvwz-nuvwz*vuvwx)/vols;
    uhklz = (nuvwy*vuvwx-nuvwx*vuvwy)/vols;
    vhklx = (nuvwy*uuvwz-nuvwz*uuvwy)/vols;
    vhkly = (nuvwz*uuvwx-nuvwx*uuvwz)/vols;
    vhklz = (nuvwx*uuvwy-nuvwy*uuvwx)/vols;

    /* nhklx:=g11*nuvwx+g12*nuvwy+g13*nuvwz;
       nhkly:=g21*nuvwx+g22*nuvwy+g23*nuvwz;
       nhklz:=g31*nuvwx+g32*nuvwy+g33*nuvwz;
       uhklx:=g11*uuvwx+g12*uuvwy+g13*uuvwz;
       uhkly:=g21*uuvwx+g22*uuvwy+g23*uuvwz;
       uhklz:=g31*uuvwx+g32*uuvwy+g33*uuvwz;
       vhklx:=g11*vuvwx+g12*vuvwy+g13*vuvwz;
       vhkly:=g21*vuvwx+g22*vuvwy+g23*vuvwz;
       vhklz:=g31*vuvwx+g32*vuvwy+g33*vuvwz; */

    /*  rotated reciprocal space base vectors  */
    /*  aed=R ae, bed=R be, ced=R ce in xyz-system  */
    //asxd = r11*asx+r12*asy+r13*asz;
    //asyd = r21*asx+r22*asy+r23*asz;
    //aszd = r31*asx+r32*asy+r33*asz;
    //bsxd = r11*bsx+r12*bsy+r13*bsz;
    //bsyd = r21*bsx+r22*bsy+r23*bsz;
    //bszd = r31*bsx+r32*bsy+r33*bsz;
    //csxd = r11*csx+r12*csy+r13*csz;
    //csyd = r21*csx+r22*csy+r23*csz;
    //cszd = r31*csx+r32*csy+r33*csz;

    /*  output matrix  */
    /* m11:=ri11;   m12:=ri12;   m13:=ri13; */
    /* m21:=ri21;   m22:=ri22;   m23:=ri23; */
    /* m31:=ri31;   m32:=ri32;   m33:=ri33; */

    /* m11:=asxd;   m12:=bsxd;   m13:=csxd; */
    /* m21:=asyd;   m22:=bsyd;   m23:=csyd; */
    /* m31:=aszd;   m32:=bszd;   m33:=cszd; */

    m11 = asedx;   m21 = bsedx;   m31 = csedx;
    m12 = asedy;   m22 = bsedy;   m32 = csedy;
    m13 = asedz;   m23 = bsedz;   m33 = csedz;

    if ( lat2d )
    {/*2*/
        m11 = asedxDD;   m21 = bsedxDD;   m31 = csedxDD;
        m12 = asedyDD;   m22 = bsedyDD;   m32 = csedyDD;
        m13 = asedzDD;   m23 = bsedzDD;   m33 = csedzDD;
        vvol = area;
    }/*2*/

    if ( lat1d )
    {/*2*/
        m11 = asedxDL;   m21 = bsedxDL;   m31 = csedxDL;
        m12 = asedyDL;   m22 = bsedyDL;   m32 = csedyDL;
        m13 = asedzDL;   m23 = bsedzDL;   m33 = csedzDL;
        vvol = area1;
    }/*2*/

    mtw11 = asedxt;   mtw21 = bsedxt;   mtw31 = csedxt;
    mtw12 = asedyt;   mtw22 = bsedyt;   mtw32 = csedyt;
    mtw13 = asedzt;   mtw23 = bsedzt;   mtw33 = csedzt;

    /*  test output  */
    nhklx = asedxn;
    nhkly = asedyn;
    nhklz = asedzn;
    uhklx = bsedxn;
    uhkly = bsedyn;
    uhklz = bsedzn;
    vhklx = csedxn;
    vhkly = csedyn;
    vhklz = csedzn;
}



// Bei den Parametern ist in /* */ der Variablenname im rufenden Programm angegeben.
void SasCalc_GENERIC_calculation::coefficients_old(int dim/*cdim*/, int nmax/*120*/, double &order/*order*/)
{/*1*/  //Z=9861

    static const int  nnmax = 120; // 500;  //Z=9865
    static const double min =  1E-40;  //Z=9866

    int    j,ll,m,n,n1,n2,n3,n4,n5,n6,n7,n8,n9,n1f,n2f,n3f,n4f,n5f,n6f,n7f,n8f,n9f,/*nmax2,*/cnt,inmax;
    int    /*ii,jj,kk,pp,qq,rr,iis,iie,jjs,jje,pps,ppe,qqs,qqe,*/ncell;
    double z,zz,zl,zzl,xlz,xl2z,xrz,xr2z,b1s,/*v,ee0,ee1,gb1s,preg1,pz2v,pz2v1,pz2v2,com,*/binsum/*,binsum1*/;
    double a1,intl,p,vvm,eps,area,vol; // p11,p12,p13,p21,p22,p23,p31,p32,p33 jetzt in params.*
    long double /*u1z,u1n,u2,u3,u4,u5,*/xrmz,xrm2z,sump,sump1,sump2,sumi,radv,aell,bell,cell;
    double carr2i[nnmax+1],fsum[nnmax+1]; //: Array[0..2*nnmax] of extended;  /*Z0311=9329*/;
    double lqm[10/*3*nnmax+1*/],qqm[10/*3*nnmax+1*/],phim[3*nnmax+1],llm[3*nnmax+1],rrm[3*nnmax+1],ppm[3*nnmax+1],
           a1m[3*nnmax+1]/*,a2m[3*nnmax+1]*/; //: Array[0..3*nnmax] of extended;  /*Z=9378*/
    double /*uell[3*nnmax+1],*/u1ell[3*nnmax+1],u2ell[3*nnmax+1],u3ell[3*nnmax+1],/*vell[3*nnmax+1],*/v1ell[3*nnmax+1+1],
           v2ell[3*nnmax+1+1],v3ell[3*nnmax+1+1],gell[3*nnmax+1+1],g3ell[3*nnmax+1+1]; //: Array[0..3*nnmax] of extended;  /*Z=9378*/
    double intlar[2*nnmax+1][2*nnmax+1]; //: Array[0..2*nnmax,0..2*nnmax] of extended;  /*Z0311=9331*/;
    //bool   search1,/*search2,*/search3,search4,search5,search6/*,search7,search8,search9*/;
    float  philiph,philipt,phiax,phiin,phiout,rad,lliph,llipt,lin,lout,len,rmax;  /*Z0311=9333*/;
    float  /*xradm,xrzm,*/x1zm,x12zm=1.0;  /*Z0311=9334*/;

    double b1sv_n = 1;

    // Arrays von oben mit erhöhter Genauigkeit
    long double gam3[3*nmax+1+1], fkv[3*nmax+1+1], fk2v[nnmax+1], fkvm[nnmax+1], pn[nnmax+1], xln[nnmax+1], z12v[3*nmax+1+1], z12vl[nnmax+1];

    // Arrays von oben, von denen aber nur das letzte Element genutzt wird.
    // Können also als normale Variablen geschrieben werden um Speicher zu sparen
    long double xrmn_n = 1.0; // , fkvm_n;

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

    params.norm = 1;  // steht (fast) überall, damit es nicht unbesetzt bleibt, hier zur Sicherheit setzen.

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
    //xrn[0] = 1;  //Z=9914
    //xrmn[0] = 1;  //Z=9915
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
        //params.CR->carr2i[n] = 1;  //Z=9952
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

        for ( int i=1; i<=2; i++ ) qqm[i] = lqm[i]/len;  //Z=9977

        phim[1] = philiph;       /*  vesicle interior  */  //Z=9979
        llm[1] = 0.0;  //Z=9980
        rrm[1] = rad+llm[1];  //Z=9981
        ppm[1] = 1.0;  //Z=9982

        radv = rrm[1];        /*  vesicle radius  */  //Z=9984
        cnt = 1;  //Z=9985
        for ( int i=1; i<=ncell; i++ )
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

        for ( int i=1; i<=inmax; i++ )
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
        for ( int i=1; i<=inmax; i++ )
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

        for ( int i=1; i<=8; i++ ) qqm[i] = lqm[i]/len;  //Z=10072

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
        for ( int i=1; i<=ncell; i++ )
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

        for ( int i=1; i<=inmax; i++ )
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
        for ( int i=1; i<=inmax; i++ )
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

//#ifndef __CUDACC__
//    qDebug() << "coefficients()" << "ordis"<<ordis << "dim"<<dim << "part"<<params.part << "cs"<<params.cs
//             << "cho1|orcase"<<params.orcase << "phi,theta" << params.polPhi << params.polTheta;
//#endif
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
            long double xrn_n = 1.0;
            //qDebug() << "TESTI" << z << xr2z << nmax;
            //n=0; std::cerr<<"TEST0 "<<n<<", "<<z12v[n]<<", "<<fkv[n]<<", "<<gam3[n]<<", "<<xrn[n]<<", "<<params.CR->carr4p[n]<<std::endl;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10209
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10210
                fkv[n] = fkv[n-1]*n;  //Z=10211
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10212
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10213 */
                xrn_n = -xrn_n*xr2z;  //Z=10214
                /*  P(q)-coefficient  */  //Z=10215
                params.CR->carr4p[n] = 9.*sqrt(M_PI)*pow(4.0,n)*(z12v[n]*xrn_n)/(2.0*(n+3.)*(n+2.)*(n+3./2.0)*gam3[n]*fkv[n]);  //Z=10216
                /*  F(q)-coefficient  */  //Z=10217
                double binsum = 0.0;  //Z=10218
                for ( m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/((m+3./2.0)*gam3[m]*(n-m+3./2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10219
                params.CR->carr4f[n] = 9*M_PI*xrn_n*binsum/16.0;  //Z=10220
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10221
                    if ( n<n4 ) n4 = n;  //Z=10222
                }/*5*/  //Z=10223
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10224
                    if ( n<n4f ) n4f = n;  //Z=10225
                }/*5*/  //Z=10226
                //if ( n < 5 ) std::cerr<<"TESTV "<<n<<", "<<z12v[n]<<", "<<fkv[n]<<", "<<gam3[n]<<", "<<xrn[n]<<", "<<params.CR->carr4p[n]<<std::endl;
            }/*4*/  //Z=10227
//#ifndef __CUDACC__
//            qDebug() << "              " << "n4"<<n4<<params.CR->carr4p[n4] << "n4f"<<n4f<<params.CR->carr4f[n4f] << "#" << params.CR->carr4f[n4];
//#endif
            goto Label99;  //Z=10228
        }/*3*/   /*  of homogeneous  */  //Z=10229

        /*  core/shell  */  //Z=10231
        if ( params.cs==1 )
        {/*3*/  //Z=10232
            long double xrmn_n = 1.0;
            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10233
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10234
                fkv[n] = fkv[n-1]*n;  //Z=10235
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10236
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10237 */
                xrn_n = -xrn_n*xr2z;  //Z=10238
                xrmn_n = -xrmn_n*xrm2z;  //Z=10239
                pn[n] = pn[n-1]*p*p;  //Z=10240
                /*  P(q)-coefficients  */  //Z=10241
                double sump = 0.0;  //Z=10242
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10243
                    sump += pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10244
                }/*5*/  //Z=10245
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10246
                params.CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;  //Z=10247
                params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=10248
                /*  F(q)-coefficients  */  //Z=10249
                sump = 0.0;  //Z=10250
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10251
                    sump += z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10252
                }/*5*/  //Z=10253
                params.CR->carr4f[n] = (9*M_PI/16.0)*xrn_n*sump;  //Z=10254
                sump = 0.0;  //Z=10255
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10256
                    sump += pn[m]*z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10257
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
            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10285
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10286
                fkv[n] = fkv[n-1]*n;  //Z=10287
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10288
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10289 */
                xrn_n = -xrn_n*xr2z;  //Z=10290
                xrmn_n = -xrmn_n*xrm2z;  //Z=10291
                pn[n] = pn[n-1]*p*p;  //Z=10292
                /*  P(q)-coefficients  */  //Z=10293
                params.CR->carr1p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10294
                double sump = 0.0;  //Z=10295
                double sump1 = 0.0;  //Z=10296
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10297
                    const double sumi = 1/((n-m+3/2.0)*gam3[n-m]*(m+3/2.0-params.alphash1/2.0)*gam3[m]*fkv[m]*fkv[n-m]);  //Z=10298
                    sump += pn[n-m]*sumi;  //Z=10299
                    sump1 += sumi;  //Z=10300
                }/*5*/  //Z=10301
                params.CR->carr2p[n] = (3*M_PI*(3-params.alphash1)/16.0)*z12v[n]*xrmn_n*sump;  //Z=10302
                params.CR->carr3p[n] = (3*M_PI*(3-params.alphash1)/16.0)*z12v[n]*xrn_n*sump1;  //Z=10303
                sump = 0.0;  //Z=10304
                sump1 = 0.0;  //Z=10305
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10306
                    const double sumi = 1/((n-m+3/2.0-params.alphash1/2.0)*(m+3/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=10307
                    sump += sumi;  //Z=10308
                    sump1 += pn[n-m]*sumi;  //Z=10309
                }/*5*/  //Z=10310
                params.CR->carr4p[n] = (sqr(3-params.alphash1)*M_PI/16.0)*z12v[n]*xrmn_n*sump;  //Z=10311
                params.CR->carr5p[n] = (sqr(3-params.alphash1)*M_PI/16.0)*z12v[n]*xrmn_n*sump1;  //Z=10312
                params.CR->carr6p[n] = (sqr(3-params.alphash1)*M_PI/16.0)*z12v[n]*xrn_n*sump;  //Z=10313

                /*  F(q)-coefficients  */  //Z=10315
                params.CR->carr4f[n] = (3*sqrt(M_PI)/4.0)*z12v[n]*xrn_n/((n+3/2.0)*gam3[n]*fkv[n]);  //Z=10316
                params.CR->carr5f[n] = (sqrt(M_PI)*(3-params.alphash1)/4.0)*z12v[n]*xrmn_n/((n+3/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=10317
                params.CR->carr6f[n] = (sqrt(M_PI)*(3-params.alphash1)/4.0)*z12v[n]*xrn_n/((n+3/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=10318
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
            //i = 2;  //Z=10352
            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
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
                xrn_n = -xrn_n*x12zm;         /*  myelin radius  */  //Z=10363

                /*  P(q)  */  //Z=10365
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10366
                    /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=10367 */
                    /* i:=i+1;  //Z=10368 */
                    params.CR->carr11pm[n][m] = (9*M_PI/16.0)*(1/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]));  //Z=10369
                }/*5*/  //Z=10370
                params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=10371
                /* carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=10372 */

                /*  F(q)  */  //Z=10375
                binsum = 0.0;  //Z=10376
                for ( m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10377
                params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=10378
                binsum = 0.0;  //Z=10379
                for ( m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=10380
                params.CR->carr4f[n] = xrn_n*binsum;  //Z=10381

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
    if ( (params.ordis==7) && (dim==4) )
    {/*2*/    /*  cubes  */  //Z=10404
        params.norm = 1;  //Z=10405
        order = 0;  //Z=10406
        /*  homogeneous  */  //Z=10407
        if ( params.cs==0 )
        {/*3*/  //Z=10408
            const double area = 6*4*sqr(params.radius);  //Z=10409
            const double vol = 8*params.radius*sqr(params.radius);  //Z=10410
            params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10411
            u1ell[0] = 2;  //Z=10413
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10414
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10415
                fkv[n] = fkv[n-1]*n;  //Z=10416
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10417
                /* u1ell[n]:=z12v[n]/((n+1/2)*(n+1)*fkv[n]);  //Z=10418 */
                u1ell[n] = 1/((n+1/2.0)*(n+1)*fkv[n]);  //Z=10419
            }/*4*/  //Z=10420

            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=10422
                double sump1 = 0.0;  //Z=10423
                for ( int m=0; m<=n; m++ ) sump1 += u1ell[n-m]*u1ell[m];  //Z=10424
                v1ell[n] = sump1;  //Z=10425
            }/*4*/  //Z=10426

            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10428
                xrn_n = -xrn_n*xr2z;  //Z=10429
                double sumf = 0.0;  //Z=10430
                for ( int m=0; m<=n; m++ ) sumf += z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10431
                fsum[n] = sumf;  //Z=10432
                /*  P(q)-coefficient  */  //Z=10433
                double sump = 0.0;  //Z=10434
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

                for ( int m=0; m<=n; m++ ) sump += u1ell[n-m]*v1ell[m];  //Z=10447

                /* carr4p[n]:=sqrt(pi)*power(4,n)*xrn[n]*sump/(16*gam3[n]);  //Z=10449 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n*sump/(16.0*gam3[n]);  //Z=10450

                /*  F(q)-coefficient  */  //Z=10452
                sump = 0.0;  //Z=10453
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10454
                    double sump1 = 0.0;  //Z=10455
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10456
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10457
                }/*5*/  //Z=10458
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=10459
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10460
                    if ( n<n4 ) n4 = n;  //Z=10461
                }/*5*/  //Z=10462
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10463
                    if ( n<n4f ) n4f = n;  //Z=10464
                }/*5*/  //Z=10465
            }/*4*/  //Z=10466
            goto Label99;  //Z=10467
        }/*3*/   /*  of homogeneous  */  //Z=10468

        /*  core/shell  */          /*  not yet ready  */  //Z=10470
        if ( params.cs==1 )
        {/*3*/  //Z=10471
            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10472
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10473
                fkv[n] = fkv[n-1]*n;  //Z=10474
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10475
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10476 */
                xrn_n = -xrn_n*xr2z;  //Z=10477
                xrmn_n = -xrmn_n*xrm2z;  //Z=10478
                pn[n] = pn[n-1]*p*p;  //Z=10479
                /*  P(q)-coefficient  */  //Z=10480
                double sump = 0.0;  //Z=10481
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10482
                    sump += pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10483
                }/*5*/  //Z=10484
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10485
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
    if ( (params.ordis==6) && (dim==4) )
    {/*2*/    /*  cubes  */  //Z=10505
        params.norm = 1;  //Z=10506
        order = 1;  //Z=10507
        long double xrn_n = 1.0;
        /*  homogeneous  */  //Z=10508
        if ( params.cs==0 )
        {/*3*/  //Z=10509
            //i = 2;  //Z=10510
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10511
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10512
                fkv[n] = fkv[n-1]*n;  //Z=10513
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10514
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10515 */
                xrn_n = -xrn_n*xr2z;  //Z=10516
                /*  P(q)-coefficient  */  //Z=10517
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10518
                    /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10519 */
                    /* carr1fm[i]:=carr1pm[i];  //Z=10520 */
                    params.CR->carr11pm[n][m] = (M_PI/4.0)*pow(4.0,n)*z12v[n-m]*z12v[m]*xrn_n/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10521
                    //i++;  //Z=10522
                }/*5*/  //Z=10523
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn_n/(16.0*gam3[n]);  //Z=10524
                /*  F(q)-coefficient  */  //Z=10525
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n/(128.0*gam3[n]);  //Z=10526
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
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10539
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10540
                fkv[n] = fkv[n-1]*n;  //Z=10541
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10542
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10543 */
                xrn_n = -xrn_n*xr2z;  //Z=10544
                xrmn_n = -xrmn_n*xrm2z;  //Z=10545
                pn[n] = pn[n-1]*p*p;  //Z=10546
                /*  P(q)-coefficient  */  //Z=10547
                double sump = 0.0;  //Z=10548
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10549
                    sump += pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10550
                }/*5*/  //Z=10551
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10552
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
    if ( (params.ordis==7) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10572
        params.norm = 1;  //Z=10573
        order = 0;  //Z=10574
        double area;
        // Die Abfrage eps==1 sollte nicht so vorkommen
        if ( eps>1 )
            area = 2*M_PI*params.radius*(params.radius+(sqr(params.length)/sqrt(sqr(params.length)-sqr(params.radius)))*asin(sqrt(sqr(params.length)-sqr(params.radius))/params.length));  //Z=10576
        else if ( eps<1 )
            area = 2*M_PI*params.radius*(params.radius+(sqr(params.length)/sqrt(sqr(params.radius)-sqr(params.length)))*asinh(sqrt(sqr(params.radius)-sqr(params.length))/params.length));  //Z=10577
        else //if ( eps==1 )
            area = 4*M_PI*sqr(params.radius);  //Z=10575
        const double vol = (4*M_PI/3.0)*sqr(params.radius)*params.length;  //Z=10578
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10579

        /*  homogeneous  */  //Z=10581
        if ( params.cs==0 )
        {/*3*/  //Z=10582
            double e1[nnmax+1];  // lokal angelegt, da nur hier verwendet
            e1[0] = 1;
            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10583
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10584
                fkv[n] = fkv[n-1]*n;  //Z=10585
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10586
                e1[n] = e1[n-1]*(eps*eps-1);  //Z=10587
                xrn_n = -xrn_n*xr2z;  //Z=10588
                double sump = 0.0;  //Z=10589
                for ( int m=0; m<=n; m++ ) sump += e1[n-m]/(fkv[n-m]*fkv[m]*(2*(n-m)+1));  //Z=10590
                /* sumf:=0.0;  //Z=10591 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10592 */
                /* fsum[n]:=sumf;  //Z=10593 */
                /*  P(q)-coefficient  */  //Z=10594
                params.CR->carr4p[n] = (9*sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn_n*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);  //Z=10595
                /*  F(q)-coefficient  */  //Z=10596
                sump = 0.0;  //Z=10597
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10598
                    double sump1 = 0.0;  //Z=10599
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10600
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10601
                }/*5*/  //Z=10602
                params.CR->carr4f[n]  = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=10603
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
    if ( (params.ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10649
        params.norm = 1;  //Z=10650
        order = 1;  //Z=10651
        /*  homogeneous  */  //Z=10652
        if ( params.cs==0 )
        {/*3*/  //Z=10653
            long double xrn_n = 1.0;
            //for ( int n=1; n<=2*nmax+2; n++ ) -> unten gehen alle Schleifen nur bis nmax
            //{/*4*/  //Z=10654
            //    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10655
            //    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10656
            //    fkv[n] = fkv[n-1]*n;  //Z=10657
            //    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10658
            //    //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=10659
            //    xrn[n] = -xrn[n-1]*xr2z;  //Z=10660
            //    xln[n] = -xln[n-1]*xl2z;  //Z=10661
            //}/*4*/  //Z=10662
            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=10663
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10655
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10656
                fkv[n] = fkv[n-1]*n;  //Z=10657
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10658
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=10659
                xrn_n = -xrn_n*xr2z;  //Z=10660
                xln[n] = -xln[n-1]*xl2z;  //Z=10661

                const double a1 = sqr(3/4.0);  //Z=10664
                for ( int m=0; m<=nmax; m++ )
                {/*5*/  //Z=10665
                    double sump1 = 0.0;  //Z=10666
                    for ( int k=0; k<=n; k++ )
                    {/*6*/  //Z=10667
                        double sump2 = 0.0;  //Z=10668
                        for ( int ll=0; ll<=m; ll++ )  //Z=10669
                            sump2 += 1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10670
                        sump1 += sump2/(fkv[n-k]*fkv[k]);  //Z=10671
                    }/*6*/  //Z=10672
                    params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10673
                }/*5*/  //Z=10674

                params.CR->carr4p[n] = a1*z12vl[n]*xln[n];  //Z=10676
                params.CR->carr5p[n] = z12v[n]*xrn_n;  //Z=10677

                /*  P(q)-coefficient  */  //Z=10679
                /* for m:=0 to n do begin  //Z=10680 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10681 */
                /* carr1fm[i]:=carr1pm[i];  //Z=10682 */
                /* i:=i+1;  //Z=10683 */
                /* end;  //Z=10684 */

                /*  F(q)-coefficient  */  //Z=10686
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n/(128.0*gam3[n]);  //Z=10687
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
    if ( (params.ordis==0) && (params.orcase==2) && (dim==5) )
    {/*2*/  //Z=10733

        if ( params.cs==0 )
        {/*3*/  //Z=10735
            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,2,0,0,0,0,params.CR->carr1p,params.norm);  //Z=10736
            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,3,0,0,0,0,params.CR->carr1p,order);  //Z=10737
            order = order/params.norm;  //Z=10738

            long double xrn_n = 1.0;
            //for ( n=1; n<=2*nmax+2; n++ )
            //{/*4*/  //Z=10740
            //    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10741
            //    z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10742
            //    fkv[n] = fkv[n-1]*n;  //Z=10743
            //    fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10744
            //    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10745
            //    xln[n] = -xln[n-1]*xl2z;  //Z=10746
            //    xrn[n] = -xrn[n-1]*xr2z;  //Z=10747
            //}/*4*/  //Z=10748

            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=10750
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10741
                z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10742
                fkv[n] = fkv[n-1]*n;  //Z=10743
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10744
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10745
                xln[n] = -xln[n-1]*xl2z;  //Z=10746
                xrn_n = -xrn_n*xr2z;  //Z=10747

                for ( int m=0; m<=nmax; m++ )
                {/*5*/  //Z=10751
                    double sump1 = 0.0;  //Z=10752
                    for ( int k=0; k<=n; k++ )
                    {/*6*/  //Z=10753
                        double sump2 = 0.0;  //Z=10754
                        for ( int ll=0; ll<=m; ll++ )  //Z=10755
                            sump2 += 1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10756
                        sump1 += sump2/(fkv[n-k]*fkv[k]);  //Z=10757
                    }/*6*/  //Z=10758
                    params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10759
                }/*5*/  //Z=10760
                double sump1 = 0.0;  //Z=10761
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10762
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=10763
                    sump1 += pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=10764
                    params.CR->carr22pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=10765
                }/*5*/  //Z=10766

                /*  all coefficients: Mok  */  //Z=10768
                /*  integral part  */  //Z=10769
                params.CR->carr4p[n] = sump1;  //Z=10770
                /*  qR-part  */  //Z=10771
                params.CR->carr5p[n] = z12v[n]*xrn_n;  //Z=10772
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
    if ( (params.ordis==7) && (dim==6) )
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
            //for ( n=1; n<=nmax; n++ )
            //{/*4*/  //Z=10834
            //    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10835
            //    fkv[n] = fkv[n-1]*n;  //Z=10836
            //    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10837
            //    xrn[n] = -xrn[n-1]*xr2z;  //Z=10838
            //    u1ell[n] = (aell*aell-bell*bell)*u1ell[n-1]/n;  //Z=10839
            //    u2ell[n] = bell*bell*u2ell[n-1]/n;  //Z=10840
            //    v1ell[n] = (bell*bell-aell*aell)*v1ell[n-1]/n;  //Z=10841
            //    v2ell[n] = (cell*cell-bell*bell)*v2ell[n-1]/n;  //Z=10842
            //    gell[n] = gam3[n]/((n+1/2.0)*fkv[n]);  //Z=10843
            //}/*4*/  //Z=10844

            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10846
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10835
                fkv[n] = fkv[n-1]*n;  //Z=10836
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10837
                xrn_n = -xrn_n*xr2z;  //Z=10838
                u1ell[n] = (aell*aell-bell*bell)*u1ell[n-1]/n;  //Z=10839
                u2ell[n] = bell*bell*u2ell[n-1]/n;  //Z=10840
                v1ell[n] = (bell*bell-aell*aell)*v1ell[n-1]/n;  //Z=10841
                v2ell[n] = (cell*cell-bell*bell)*v2ell[n-1]/n;  //Z=10842
                gell[n] = gam3[n]/((n+1/2.0)*fkv[n]);  //Z=10843

                double sump = 0.0;  //Z=10859
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10860
                    double sump1 = 0.0;  //Z=10861
                    for ( int k=0; k<=m; k++ )
                    {/*6*/  //Z=10862
                        double sump2 = 0.0;  //Z=10863
                        for ( int ll=0; ll<=n-m; ll++ ) sump2 += gell[k+ll]*v1ell[ll]*v2ell[n-m-ll];  //Z=10864
                        sump1 += u1ell[k]*u2ell[m-k]*sump2;  //Z=10865
                    }/*6*/  //Z=10866
                    sump += sump1/((2.0*(n-m))+1);  //Z=10867
                }/*5*/  //Z=10868

                /* fsum[n]:=sumf;  //Z=10870 */
                double sumf = 0.0;  //Z=10871
                for ( int m=0; m<=n; m++ ) sumf += z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10872
                fsum[n] = sumf;  //Z=10873
                /*  P(q)-coefficient  */  //Z=10874
                params.CR->carr4p[n] = (9*M_PI)*pow(4.0,n-1)*z12v[n]*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);  //Z=10875
                /*  F(q)-coefficient  */  //Z=10876
                sump = 0.0;  //Z=10877
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10878
                    double sump1 = 0.0;  //Z=10879
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10880
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10881
                }/*5*/  //Z=10882
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=10883
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
    if ( (params.ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10929
        params.norm = 1;  //Z=10930
        order = 1;  //Z=10931
        /*  homogeneous  */  //Z=10932
        if ( params.cs==0 )
        {/*3*/  //Z=10933
            //i = 2;  //Z=10934
            long double xrn[nmax+2];
            xrn[0] = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10935
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10936
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10937
                fkv[n] = fkv[n-1]*n;  //Z=10938
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10939
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=10940
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10941
                xln[n] = -xln[n-1]*xl2z;  //Z=10942
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10943
                    double sump = 0.0;  //Z=10944
                    for ( int ns=0; ns<=n; ns++ )
                    {/*6*/  //Z=10945
                        double sump1 = 0.0;  //Z=10946
                        for ( int ms=0; ms<=m; ms++ ) sump1 += 1.0/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);  //Z=10947
                        sump += sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);  //Z=10948
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
    if ( (params.ordis==7) && (dim==7) )
    {/*2*/    /*  super ellipsoids  */  //Z=11007

        //double qx=1, qy=1, qz=1, qhkl=1;
        // TODO Diese Variablen sind im Orginal-Pascalprogramm an der Stelle des coefficients() Aufrufes nicht bekannt.

        double xrn_n = 1.0;
        // TODO Dieses Array xrn[] wird unten verwendet, aber nicht berechnet....

        params.norm = 1;  //Z=11008
        order = 0;  //Z=11009
        double area, a1;
        qrombchid(params.length,params.radius,params.alphash1,params.sigma,params.alphash1,params.polPhi,params.polTheta,params.polPhi,
                  1/*qx*/,1/*qy*/,1/*qz*/,
                  params.p11,params.p12,params.p13,params.p21,params.p22,params.p23,params.p31,params.p32,params.p33,
                  9/*qx*/,9/*qy*/,9/*0*/,9/*qhkl*/,
                  //params.ax1.length(),params.ax2.length(),params.ax3.length(),
                  //params.ax1.x(),params.ax1.y(),params.ax1.z(),
                  //params.ax2.x(),params.ax2.y(),params.ax2.z(),
                  //params.ax3.x(),params.ax3.y(),params.ax3.z(),
                  //params.sig.x(),params.sig.y(),params.sig.z(),
                  params.ordis,3,8,15,7,0,0,params.CR->carr1p,area);   //Z=11010
        area = 2*M_PI*area;  //Z=11011
        const double vol = 2*M_PI*sqr(params.radius)*params.length*gamma((2+params.alphash1)/params.alphash1)*gamma(1/params.alphash1)/(params.alphash1*gamma((3+params.alphash1)/params.alphash1));  //Z=11012
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=11013

        /*  homogeneous  */  //Z=11015
        u1ell[0] = 1;             /*  for barrel  */  //Z=11016
        u2ell[0] = 1;             /*  for barrel  */  //Z=11017
        for ( int n=1; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11018
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11019
            fkv[n] = fkv[n-1]*n;  //Z=11020
            gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11021
            u1ell[n] = u1ell[n-1]/((n-1/2.0)*n);   /*  for barrel  */  //Z=11022
            u2ell[n] = u2ell[n-1]/((n+1)*n);     /*  for barrel  */  //Z=11023
        }/*3*/  //Z=11024

        for ( int n=0; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11026
            v1ell[n] = gamma((2*n+1)/params.alphash1)*u1ell[n];         /*  for barrel  */  //Z=11027
            v2ell[n] = gamma((2*n+2+params.alphash1)/params.alphash1)*u2ell[n];    /*  for barrel  */  //Z=11028
            gell[n] = gamma((2*n+3+params.alphash1)/params.alphash1);              /*  for barrel  */  //Z=11029
        }/*3*/  //Z=11030

        if ( params.cs==0 )
        {/*3*/  //Z=11032
            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=11033

                if ( params.alphash1==2 )
                {/*5*/  //Z=11047
                    a1 = sqr(3/4.0);  //Z=11048
                    for ( int m=0; m<=nmax; m++ )
                    {/*6*/  //Z=11049
                        double sump1 = 0.0;  //Z=11050
                        for ( int k=0; k<=n; k++ )
                        {/*7*/  //Z=11051
                            double sump2 = 0.0;  //Z=11052
                            for ( int ll=0; ll<=m; ll++ )  //Z=11053
                                sump2 += 1.0/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=11054
                            sump1 += sump2/(fkv[n-k]*fkv[k]);  //Z=11055
                        }/*7*/  //Z=11056
                        params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=11057
                    }/*6*/  //Z=11058
                }/*5*/  //Z=11059
                else
                {/*5*/  //Z=11060
                    a1 = gamma((params.alphash1+3)/params.alphash1)/(gamma((params.alphash1+2)/params.alphash1)*gamma(1/params.alphash1));  //Z=11061
                    a1 = a1*a1;  //Z=11062
                    for ( int m=0; m<=nmax; m++ )
                    {/*6*/  //Z=11063
                        double sump1 = 0.0;  //Z=11064
                        for ( int ns=0; ns<=n; ns++ )
                        {/*7*/  //Z=11065
                            double sump2 = 0.0;  //Z=11066
                            for ( int ms=0; ms<=m; ms++ )  //Z=11067
                                sump2 += v2ell[m-ms]*v2ell[ms]/(gell[ns+ms]*gell[n-ns+m-ms]);  //Z=11068
                            sump1 += v1ell[n-ns]*v1ell[ns]*sump2;  //Z=11069
                        }/*7*/  //Z=11070
                        params.CR->carr11pm[n][m] = sump1;  //Z=11071
                    }/*6*/  //Z=11072
                }/*5*/  //Z=11073

                /*  orientational average  */  //Z=11075
                double sump = 0.0;  //Z=11076
                for ( int m=0; m<=n; m++ )
                    sump += gam3[n-m]*z12v[n-m]*z12v[m]*pow(sqr(params.length),n-m)*pow(sqr(params.radius),m)*fkv[m]*params.CR->carr11pm[n-m][m]/(n-m+1/2.0);  //Z=11077
                params.CR->carr4p[n] = a1*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/(2.0*gam3[n]);  //Z=11078

                /* fsum[n]:=sumf;  //Z=11081 */
                /* sumf:=0.0;  //Z=11082 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11083 */
                /* fsum[n]:=sumf;  //Z=11084 */
                /*  P(q)-coefficient  */  //Z=11085
                /*  F(q)-coefficient  */  //Z=11086
                sump = 0.0;  //Z=11087
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=11088
                    double sump1 = 0.0;  //Z=11089
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=11090
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=11091
                }/*5*/  //Z=11092
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=11093
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
    if ( (params.ordis==6) && (dim==5) )
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
    if ( (params.ordis==7) && (dim==8) )
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

        double xrn_n = 1.0;
        // TODO Dieses Array xrn[] wird unten verwendet, aber nicht berechnet....

        float carr111pm[ardim+1][ardim+1][ardim+1];

        nmax = 40;  //Z=11231 - war Parameter (=120)
        params.norm = 1;  //Z=11232
        order = 0;  //Z=11233
        /* l:=r;  //Z=11234 */

        /*  radius=a, rm=b, length=c  */  //Z=11236
        qrombdeltac(params.radiusi,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,8,17,0,0,0,params.CR->carr1p,area);  //Z=11237
        area = 8*area;  //Z=11238
        vol = 8*params.radius*params.radiusi*params.length*pow(gamma(1/params.alphash1),3)/(pow(params.alphash1,3)*gamma(1+3/params.alphash1));  //Z=11239
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=11240

        /*  homogeneous  */  //Z=11243
        u3ell[0] = 1;             /*  for superball  */  //Z=11244
        for ( int n=1; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11245
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11246     ERROR: Iteration 120 invokes undefined behaviour
            fkv[n] = fkv[n-1]*n;  //Z=11247
            gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11248
            u3ell[n] = u3ell[n-1]/((n-1/2.0)*n);   /*  for superball  */  //Z=11249
        }/*3*/  //Z=11250

        /* for n:=0 to 3*nmax+1 do begin  //Z=11252 */
        /*    v3ell[n]:=gamma((2*n+1)/alfa)*u3ell[n];         (* for superball *)  //Z=11253 */
        /*    g3ell[n]:=gamma(((2*n+3)/alfa)+1);              (* for superball *)  //Z=11254 */
        /* end;  //Z=11255 */

        for ( int n=0; n<=3*nmax+1; n++ )
        {/*3*/  //Z=11257
            v3ell[n] = exp(gammln((2*n+1)/params.alphash1))*u3ell[n];         /*  for superball  */  //Z=11258
            g3ell[n] = exp(gammln(((2*n+3)/params.alphash1)+1));              /*  for superball  */  //Z=11259
        }/*3*/  //Z=11260

        if ( params.cs==0 )
        {/*3*/  //Z=11262
            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=11263
                if ( params.alphash1==200000 )
                {/*5*/  //Z=11264
                    a1 = sqr(3/4.0);  //Z=11265
                    for ( int m=0; m<=nmax; m++ )
                    {/*6*/  //Z=11266
                        sump1 = 0.0;  //Z=11267
                        for ( int k=0; k<=n; k++ )
                        {/*7*/  //Z=11268
                            sump2 = 0.0;  //Z=11269
                            for ( int ll=0; ll<=m; ll++ )  //Z=11270
                                sump2 += 1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=11271
                            sump1 += sump2/(fkv[n-k]*fkv[k]);  //Z=11272
                        }/*7*/  //Z=11273
                        params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=11274
                    }/*6*/  //Z=11275
                }/*5*/  //Z=11276
                else
                {/*5*/  //Z=11277
                    a1 = gamma(1+3/params.alphash1)/pow(gamma(1/params.alphash1),3);  //Z=11278
                    a1 = a1*a1;  //Z=11279
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11280
                        for ( int k=0; k<=nmax; k++ )
                        {/*7*/  //Z=11281
                            double sump1 = 0.0;  //Z=11282
                            for ( int ns=0; ns<=n; ns++ )
                            {/*8*/  //Z=11283
                                double sump2 = 0.0;  //Z=11284
                                for ( int ms=0; ms<=m; ms++ )
                                {/*9*/  //Z=11285
                                    double sump3 = 0.0;  //Z=11286
                                    for ( int ks=0; ks<=k; ks++ )  //Z=11287
                                        sump3 += v3ell[k-ks]*v3ell[ks]/(g3ell[ns+ms+ks]*g3ell[n-ns+m-ms+k-ks]);  //Z=11288
                                    sump2 += v3ell[m-ms]*v3ell[ms]*sump3;  //Z=11289
                                }/*9*/  //Z=11290
                                sump1 += v3ell[n-ns]*v3ell[ns]*sump2;  //Z=11291
                            }/*8*/  //Z=11292
                            carr111pm[n][m][k] = sump1;  //Z=11293
                            carr111pm[m][n][k] = sump1;  //Z=11294
                        }/*7*/  //Z=11295
                    }/*6*/  //Z=11296
                }/*5*/  //Z=11297

                /*  orientational average for superball  */  //Z=11299
                double sump = 0.0;  //Z=11300
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=11301
                    double sump1 = 0.0;  //Z=11302
                    for ( int k=0; k<=m; k++ )  //Z=11303
                        sump1 += z12v[m-k]*z12v[k]*gam3[m-k]*gam3[k]*pow(sqr(params.radiusi),m-k)*pow(sqr(params.length),k)*carr111pm[n-m][m-k][k]/((m-k+1/2.0)*(k+1/2.0));   //Z=11304
                    sump += z12v[n-m]*gam3[n-m]*pow(sqr(params.radius),n-m)*sump1/(n-m+1/2.0);  //Z=11305
                }/*5*/  //Z=11306
                params.CR->carr4p[n] = (a1/(2.0*M_PI))*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/gam3[n];  //Z=11307


                /* fsum[n]:=sumf;  //Z=11310 */
                /* sumf:=0.0;  //Z=11311 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11312 */
                /* fsum[n]:=sumf;  //Z=11313 */
                /*  P(q)-coefficient  */  //Z=11314
                /*  F(q)-coefficient  */  //Z=11315
                sump = 0.0;  //Z=11316
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=11317
                    double sump1 = 0.0;  //Z=11318
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=11319
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=11320
                }/*5*/  //Z=11321
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=11322
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
    if ( (params.ordis==6) && (dim==5) )
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
    if ( (params.ordis==7) && (dim!=3) )
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

                long double xrn[nmax+1];
                xrn[0] = 1.0;

                for ( int n=1; n<=nmax; n++ )
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

                for ( int n=0; n<=nmax; n++ )
                {/*5*/  //Z=11472
                    double binsum = 0.0;  //Z=11473
                    /* for m:=0 to n do    (* Cauchy sum *)  //Z=11474 */
                    /*       binsum:=binsum+gam3[m]*z12vl[n-m]*z12v[m]*xln[n-m]*xrn[m]/((n-m+1/2)*(n-m+1)*fkv[n-m]*(m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]);  //Z=11475 */
                    params.CR->carr11pm[n][0] = sqrt(M_PI)/gam3[n];  //Z=11476
                    double a1 = sqrt(M_PI)/(2.0*gam3[n]);  //Z=11477
                    for ( int m=1; m<=nmax; m++ )
                    {/*6*/    /*  double sum  */  //Z=11478
                        a1 = a1*(m+1/2.0)/(n+m+1/2.0);  //Z=11479
                        /* carr11pm[n,m]:=power(4,m+1)*gam3[m]*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]*gam3[n+m]);   (* ok *)  //Z=11480 */
                        params.CR->carr11pm[n][m] = pow(4.0,m+1)*a1*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]);       /*  Mok  */  //Z=11481
                    }/*6*/  //Z=11482
                    params.CR->carr2p[n] = pow(4.0,n-1)*z12vl[n]*xln[n]/((n+1/2.0)*(n+1)*fkv[n]);     /*  double sum  */  //Z=11483
                    params.CR->carr3p[n] = pow(4.0,n)*binsum/gam3[n]; // TODO binsum wird nicht berechnet     /*  Cauchy sum  */  //Z=11484
                }/*5*/  //Z=11485

            }/*4*/  //Z=11487
            /*  large axial ratios  */  //Z=11488
            else
            {/*4*/  //Z=11489

                long double xrn_n = 1.0;

                for ( int n=1; n<=nmax; n++ )
                {/*5*/  //Z=11491
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11492
                    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11493
                    b1sv_n = b1sv_n*(b1s-1+n);  //Z=11494
                    fkv[n] = fkv[n-1]*n;  //Z=11495
                    fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11496
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11497
                    /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11498 */
                    xln[n] = -xln[n-1]*xl2z;  //Z=11499
                    xrn_n = -xrn_n*xr2z;  //Z=11500
                    /*  cylinder, ok */  //Z=11501
                    if ( dim==1 )
                    {/*6*/  //Z=11502

                        /*  P(q) factorization  */  //Z=11504
                        params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);      /*  P||iso(q)  */  //Z=11505
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*  P-(q)  */  //Z=11506
                        /*  F(q)  */  //Z=11507
                        double binsum = 0.0;  //Z=11508
                        for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11509
                        params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11510
                        binsum = 0.0;  //Z=11511
                        for ( int m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11512
                        params.CR->carr4f[n] = xrn_n*binsum;  //Z=11513
                    }/*6*/  //Z=11514
                    /*  disk, ok  */  //Z=11515
                    if ( dim==2 )
                    {/*6*/  //Z=11516
                        /*  P(q)  */  //Z=11517
                        params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11518
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11519
                        /*  F(q)  */  //Z=11520
                        double binsum = 0.0;  //Z=11521
                        for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11522
                        params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11523
                        binsum = 0.0;  //Z=11524
                        for ( int m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11525
                        params.CR->carr4f[n] = M_PI*xrn_n*binsum/4.0;  //Z=11526
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
            long double xrn_n = 1.0;
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
                xrn_n = -xrn_n*xr2z;  //Z=11562
                xrmn_n = -xrmn_n*xrm2z;  //Z=11563
                pn[n] = pn[n-1]*p*p;  //Z=11564
                /* ** cylinder ** */  //Z=11565
                if ( dim==1 )
                {/*5*/  //Z=11566
                    /*  longitudinal P(q)  */  //Z=11567
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11568
                    /*  cross-sectional P(q)  */  //Z=11569
                    /*  F121  */  //Z=11570
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11571
                    /*  F122  */  //Z=11572
                    double sump = 0.0;  //Z=11573
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11574
                        sump += pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11575
                    }/*6*/  //Z=11576
                    params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=11577
                    /*  F123  */  //Z=11578
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11579

                    /*  longitudinal F(q)  */  //Z=11581
                    double binsum = 0.0;  //Z=11582
                    for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11583
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11584
                    /*  cross-sectional F(q)  */  //Z=11585
                    /*  F121  */  //Z=11586
                    sump = 0.0;  //Z=11587
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11588
                        sump += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11589
                    }/*6*/  //Z=11590
                    params.CR->carr4f[n] = xrn_n*sump;  //Z=11591
                    /*  F122  */  //Z=11592
                    sump = 0.0;  //Z=11593
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11594
                        sump += z12v[m]*z12v[n-m]*pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11595
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
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11608
                    /*  F122  */  //Z=11609
                    double sump = 0.0;  //Z=11610
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11611
                        sump += pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11612
                    }/*6*/  //Z=11613
                    params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=11614
                    /*  F123  */  //Z=11615
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11616

                    /*  longitudinal F(q)  */  //Z=11618
                    double binsum = 0.0;  //Z=11619
                    for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11620
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11621
                    /*  cross-sectional F(q)  */  //Z=11622
                    /*  F121  */  //Z=11623
                    sump = 0.0;  //Z=11624
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=11625
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11626
                    }/*6*/  //Z=11627
                    params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=11628
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
            long double xrn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=11670
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11671
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11672
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=11673
                fkv[n] = fkv[n-1]*n;  //Z=11674
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11675
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11676
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11677 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11678
                xrn_n = -xrn_n*xr2z;  //Z=11679
                xrmn_n = -xrmn_n*xrm2z;  //Z=11680
                pn[n] = pn[n-1]*p*p;  //Z=11681
                /* ** cylinder ** */  //Z=11682
                if ( dim==1 )
                {/*5*/  //Z=11683
                    /*  longitudinal P(q)  */  //Z=11684
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11685

                    /*  cross-sectional P(q)  */  //Z=11687
                    params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn_n/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11688
                    double sump = 0.0;  //Z=11689
                    double sump1 = 0.0;  //Z=11690
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11691
                        const double sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=11692
                        sump += pn[n-m]*sumi;  //Z=11693
                        sump1 += sumi;  //Z=11694
                    }/*6*/  //Z=11695
                    params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=11696
                    params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn_n*sump1;  //Z=11697
                    sump = 0.0;  //Z=11698
                    sump1 = 0.0;  //Z=11699
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11700
                        const double sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=11701
                        sump += sumi;  //Z=11702
                        sump1 += pn[n-m]*sumi;  //Z=11703
                    }/*6*/  //Z=11704
                    params.CR->carr7p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=11705
                    params.CR->carr8p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=11706
                    params.CR->carr9p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrn_n*sump;  //Z=11707


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
                    params.CR->carr4f[n] = z12v[n]*xrn_n/((n+1)*fkv[n]*fkv[n]);  //Z=11727
                    params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11728
                    params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11729

                }/*5*/  //Z=11731

                /* ** disk ** */  //Z=11733
                if ( dim==2 )
                {/*5*/  //Z=11734
                    /*  longitudinal  */  //Z=11735
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11736


                    /*  cross-sectional P(q)  */  //Z=11739
                    params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn_n/((n+1)*gam3[n]*fkv[n]);  //Z=11740
                    double sump = 0.0;  //Z=11741
                    double sump1 = 0.0;  //Z=11742
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11743
                        const double sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11744
                        sump += pn[n-m]*sumi;  //Z=11745
                        sump1 += sumi;  //Z=11746
                    }/*6*/  //Z=11747
                    params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=11748
                    params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn_n*sump1;  //Z=11749
                    sump = 0.0;  //Z=11750
                    sump1 = 0.0;  //Z=11751
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11752
                        const double sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=11753
                        sump += sumi;  //Z=11754
                        sump1 += pn[n-m]*sumi;  //Z=11755
                    }/*6*/  //Z=11756
                    params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;   //Z=11757
                    params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=11758
                    params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn_n*sump;   //Z=11759


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
                    double binsum = 0.0;  //Z=11775
                    for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11776
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11777
                    /*  cross-sectional F(q)  */  //Z=11778
                    params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn_n/(gam3[n]*fkv[n]);  //Z=11779
                    params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11780
                    params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11781
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
            //i = 2;  //Z=11831
            long double xrn_n = 1.0;
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
                xrn_n = -xrn_n*x12zm;         /*  myelin radius  */  //Z=11842
                /*  cylinder, ok */  //Z=11843
                if ( dim==1 )
                {/*5*/  //Z=11844
                    /*  P(q)  */  //Z=11845
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11846

                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11848
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11849 */
                        /* i:=i+1;  //Z=11850 */
                        params.CR->carr11pm[n][m] = 1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11851
                    }/*6*/  //Z=11852
                    params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=11853
                    /* carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11854 */


                    /*  F(q)  */  //Z=11857
                    double binsum = 0.0;  //Z=11858
                    for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11859
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11860
                    binsum = 0.0;  //Z=11861
                    for ( int m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11862
                    params.CR->carr4f[n] = xrn_n*binsum;  //Z=11863
                }/*5*/  //Z=11864
                /*  disk, ok  */  //Z=11865
                if ( dim==2 )
                {/*5*/  //Z=11866
                    /*  P(q)  */  //Z=11867
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11868
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11869
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11870 */
                        /* i:=i+1;  //Z=11871 */
                        params.CR->carr11pm[n][m] = (M_PI/4.0)*(1/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]));  //Z=11872
                    }/*6*/  //Z=11873
                    params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=11874

                    /* carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=11876 */
                    /*  F(q)  */  //Z=11877
                    double binsum = 0.0;  //Z=11878
                    for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11879
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11880
                    binsum = 0.0;  //Z=11881
                    for ( int m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11882
                    params.CR->carr4f[n] = M_PI*xrn_n*binsum/4.0;  //Z=11883
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
    if ( (params.ordis==6) && (dim!=3) )
    {/*2*/  //Z=11903
        params.norm = 1;  //Z=11904
        order = 1;  //Z=11905
        long double xrn_n = 1.0;
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
                xrn_n = -xrn_n*xr2z;  //Z=11917
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11918
                    /*  P(q)-coefficients  */  //Z=11919
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);       /*  P||(q)  */  //Z=11920
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  /*  P-(q)  */  //Z=11921
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
                    params.CR->carr4f[n] = xrn_n*sump;  //Z=11932
                }/*5*/  //Z=11933
                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=11934
                    /*  P(q)-coefficients  */  //Z=11935
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);    /*  P-(q)  */  //Z=11936
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);           /*  P||(q)  */  //Z=11937
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
                xrn_n = -xrn_n*xr2z;  //Z=11976
                xrmn_n = -xrmn_n*xrm2z;  //Z=11977
                pn[n] = pn[n-1]*p*p;  //Z=11978
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11979
                    /*  P(q)-coefficients  */  //Z=11980
                    /*  longitudinal  */  //Z=11981
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11982
                    /*  cross-sectional  */  //Z=11983
                    /*  F121  */  //Z=11984
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11985
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
                    params.CR->carr4f[n] = xrn_n*sump;  //Z=12008
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
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12025
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
                    params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=12056
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
                xrn_n = -xrn_n*xr2z;  //Z=12104
                xrmn_n = -xrmn_n*xrm2z;  //Z=12105
                pn[n] = pn[n-1]*p*p;  //Z=12106
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=12107
                    /*  P(q)-coefficients  */  //Z=12108
                    /*  longitudinal  */  //Z=12109
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12110

                    /*  cross-sectional P(q)  */  //Z=12112
                    params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn_n/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12113
                    sump = 0.0;  //Z=12114
                    sump1 = 0.0;  //Z=12115
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12116
                        sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=12117
                        sump = sump+pn[n-m]*sumi;  //Z=12118
                        sump1 = sump1+sumi;  //Z=12119
                    }/*6*/  //Z=12120
                    params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=12121
                    params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn_n*sump1;  //Z=12122
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
                    params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn_n*sump;  //Z=12132

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
                    params.CR->carr4f[n] = z12v[n]*xrn_n/((n+1)*fkv[n]*fkv[n]);  //Z=12154
                    params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12155
                    params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12156
                }/*5*/  //Z=12157

                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=12159
                    /*  P(q)-coefficients  */  //Z=12160
                    /*  longitudinal  */  //Z=12161
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12162

                    /*  cross-sectional P(q)  */  //Z=12164
                    params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn_n/((n+1)*gam3[n]*fkv[n]);  //Z=12165
                    sump = 0.0;  //Z=12166
                    sump1 = 0.0;  //Z=12167
                    for ( m=0; m<=n; m++ )
                    {/*6*/  //Z=12168
                        sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=12169
                        sump = sump+pn[n-m]*sumi;  //Z=12170
                        sump1 = sump1+sumi;  //Z=12171
                    }/*6*/  //Z=12172
                    params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=12173
                    params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn_n*sump1;  //Z=12174
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
                    params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn_n*sump;  //Z=12184

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
                    params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn_n/(gam3[n]*fkv[n]);  //Z=12214
                    params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12215
                    params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12216

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
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=12267
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12268
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12269
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=12270
                fkv[n] = fkv[n-1]*n;  //Z=12271
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12272
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12273
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12274 */
                xln[n] = -xln[n-1]*xl2z;  //Z=12275
                xrn_n = -xrn_n*x12zm;  //Z=12276
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=12277
                    /*  P(q)-coefficients  */  //Z=12278
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12279

                    params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=12281

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
                    params.CR->carr4f[n] = xrn_n*sump;  //Z=12293
                }/*5*/  //Z=12294
                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=12295
                    /*  P(q)-coefficients  */  //Z=12296
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12297
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12298
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
    if ( (params.ordis==0) && (dim!=3) )
    {/*2*/  //Z=12331
        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,2,0,0,0,0,params.CR->carr1p,params.norm);  //Z=12332
        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,3,0,0,0,0,params.CR->carr1p,order);  //Z=12333
        order = order/params.norm;  //Z=12334

        long double xrn_n = 1.0;

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

            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=12349
                for ( int m=0; m<=nmax; m++ )
                {/*5*/  //Z=12350
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,0,0,2*n,2*m,params.CR->carr1p,intl);  //Z=12351
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
                    //i = 2;  //Z=12371
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12372
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12373
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12374
                        b1sv_n = b1sv_n*(b1s-1+n);  //Z=12375
                        fkv[n] = fkv[n-1]*n;  //Z=12376
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12377
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12378
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12379 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12380
                        xrn_n = -xrn_n*xr2z;  //Z=12381
                        /*  longitudinal  */  //Z=12382
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12383
                            double sump = 0.0;  //Z=12384
                            for ( int ll=0; ll<=m; ll++ )  //Z=12385
                                sump += pow(4.0,ll)*pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[ll][n-ll]/(fk2v[ll]*fkv[m-ll]*fkv[n-ll]*params.norm);  //Z=12386
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
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=12406
                        sump = 0.0;  //Z=12407
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12408
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=12409
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
                    //i = 2;  //Z=12427
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
                        xrn_n = -xrn_n*xr2z;  //Z=12437
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
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12466
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
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=12483
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
                    //i = 2;  //Z=12524
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
                        xrn_n = -xrn_n*xr2z;  //Z=12534
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
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn_n/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12562
                        sump = 0.0;  //Z=12563
                        sump1 = 0.0;  //Z=12564
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12565
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=12566
                            sump = sump+pn[n-m]*sumi;  //Z=12567
                            sump1 = sump1+sumi;  //Z=12568
                        }/*7*/  //Z=12569
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=12570
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn_n*sump1;  //Z=12571
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
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn_n*sump;  //Z=12581

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
                        params.CR->carr4f[n] = z12v[n]*xrn_n/((n+1)*fkv[n]*fkv[n]);  //Z=12599
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12600
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12601

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
                    //i = 2;  //Z=12650
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
                        xrn_n = -xrn_n*x12zm;  //Z=12660
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
                        params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=12684

                        sump = 0.0;  //Z=12686
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12687
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=12688
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

            /*  disk form factor coefficients  */  //Z=12730
            if ( dim==2 )
            {/*4*/  //Z=12731
                /*  homogeneous  */  //Z=12732
                params.p11 = sin(params.polTheta*M_PI/180.0);  //Z=12733
                params.p13 = cos(params.polTheta*M_PI/180.0);  //Z=12734
                if ( params.cs==0 )
                {/*5*/  //Z=12735
                    //i = 2;  //Z=12736
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12737
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12738
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12739
                        fkv[n] = fkv[n-1]*n;  //Z=12740
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12741
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12742
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12743 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12744
                        xrn_n = -xrn_n*xr2z;  //Z=12745
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
                            //i = i+1;  //Z=12757
                        }/*7*/  //Z=12758
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=12759

                        sump1 = 0.0;  //Z=12761
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12762
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12763
                        }/*7*/  //Z=12764
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12765

                        /*  cross-sectional  */  //Z=12767
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12768
                        sump = 0.0;  //Z=12769
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12770
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12771
                        }/*7*/  //Z=12772
                        params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=12773

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
                    //i = 2;  //Z=12800
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12801
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12802
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12803
                        fkv[n] = fkv[n-1]*n;  //Z=12804
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12805
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12806
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12807 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12808
                        xrn_n = -xrn_n*xr2z;  //Z=12809
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
                            //i = i+1;  //Z=12823
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
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12836
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
                        params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=12860
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
                    //i = 2;  //Z=12907
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12908
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12909
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12910
                        fkv[n] = fkv[n-1]*n;  //Z=12911
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12912
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12913
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12914 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12915
                        xrn_n = -xrn_n*xr2z;  //Z=12916
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
                            //i = i+1;  //Z=12930
                        }/*7*/  //Z=12931
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=12932

                        sump1 = 0.0;  //Z=12934
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12935
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12936
                        }/*7*/  //Z=12937
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12938

                        /*  cross-sectional P(q)  */  //Z=12940
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn_n/((n+1)*gam3[n]*fkv[n]);  //Z=12941
                        sump = 0.0;  //Z=12942
                        sump1 = 0.0;  //Z=12943
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=12944
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=12945
                            sump = sump+pn[n-m]*sumi;  //Z=12946
                            sump1 = sump1+sumi;  //Z=12947
                        }/*7*/  //Z=12948
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=12949
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn_n*sump1;  //Z=12950
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
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn_n*sump;  //Z=12960

                        /*  F(q)  */  //Z=12981
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn_n/(gam3[n]*fkv[n]);  //Z=12982
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12983
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12984

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

            //i = 2;  //Z=13176
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=13177
                for ( m=0; m<=2*n; m++ )
                {/*5*/  //Z=13178
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,5,0,m,2*n-m,params.CR->carr1p,intl);  //Z=13179
                    /* carr1pm[i]:=intl/(fkv[m]*fkv[2*n-m]*norm);  //Z=13180 */
                    //i = i+1;  //Z=13181
                }/*5*/  //Z=13182
            }/*4*/  //Z=13183

            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=13185
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13186
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=13187
                xln[n] = -xln[n-1]*xl2z;  //Z=13188
                xrn_n = -xrn_n*xr2z;  //Z=13189
                params.CR->carr1p[n] = pow(4.0,2*n)*z12v[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=13190
                params.CR->carr3p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*b1sv_n*fkv[n]);  //Z=13191

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
                    //i = 2;  //Z=13211
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
                        xrn_n = -xrn_n*xr2z;  //Z=13221
                        /*  longitudinal  */  //Z=13222
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13223
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13224
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
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=13241
                        sump = 0.0;  //Z=13242
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13243
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=13244

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
                    //i = 2;  //Z=13263
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
                        xrn_n = -xrn_n*xr2z;  //Z=13273
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13274
                        pn[n] = pn[n-1]*p*p;  //Z=13275
                        /*  longitudinal  */  //Z=13276
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13277
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13278
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
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13296
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
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=13313
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
                    //i = 2;  //Z=13352
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
                        xrn_n = -xrn_n*xr2z;  //Z=13362
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13363
                        pn[n] = pn[n-1]*p*p;  //Z=13364
                        /*  longitudinal  */  //Z=13365
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13366
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13367
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
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn_n/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13384
                        sump = 0.0;  //Z=13385
                        sump1 = 0.0;  //Z=13386
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13387
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=13388
                            sump = sump+pn[n-m]*sumi;  //Z=13389
                            sump1 = sump1+sumi;  //Z=13390
                        }/*7*/  //Z=13391
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=13392
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn_n*sump1;  //Z=13393
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
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn_n*sump;  //Z=13403

                        /*  F(q)-coefficients  */  //Z=13418
                        /*  cross-sectional  */  //Z=13419
                        params.CR->carr4f[n] = z12v[n]*xrn_n/((n+1)*fkv[n]*fkv[n]);  //Z=13420
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=13421
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=13422

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
                    //i = 2;  //Z=13471
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
                        xrn_n = -xrn_n*x12zm;  //Z=13481
                        /*  longitudinal  */  //Z=13482
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13483
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13484
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
                        params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=13500

                        sump = 0.0;  //Z=13502
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13503
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=13504

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
                    //i = 2;  //Z=13526
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13527
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13528
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13529
                        fkv[n] = fkv[n-1]*n;  //Z=13530
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13531
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13532
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13533 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13534
                        xrn_n = -xrn_n*xr2z;  //Z=13535
                        /*  longitudinal  */  //Z=13536
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13537
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13538
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
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13558
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
                    //i = 2;  //Z=13590
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13591
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13592
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13593
                        fkv[n] = fkv[n-1]*n;  //Z=13594
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13595
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13596
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13597 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13598
                        xrn_n = -xrn_n*xr2z;  //Z=13599
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13600
                        pn[n] = pn[n-1]*p*p;  //Z=13601
                        /*  longitudinal  */  //Z=13602
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13603
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13604
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
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13624
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
                        params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=13647
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
                    //i = 2;  //Z=13694
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13695
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13696
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13697
                        fkv[n] = fkv[n-1]*n;  //Z=13698
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13699
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13700
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13701 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13702
                        xrn_n = -xrn_n*xr2z;  //Z=13703
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13704
                        pn[n] = pn[n-1]*p*p;  //Z=13705
                        /*  longitudinal  */  //Z=13706
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13707
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13708
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
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn_n/((n+1)*gam3[n]*fkv[n]);  //Z=13728
                        sump = 0.0;  //Z=13729
                        sump1 = 0.0;  //Z=13730
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13731
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=13732
                            sump = sump+pn[n-m]*sumi;  //Z=13733
                            sump1 = sump1+sumi;  //Z=13734
                        }/*7*/  //Z=13735
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=13736
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn_n*sump1;  //Z=13737
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
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn_n*sump;  //Z=13747

                        /*  F(q)   */  //Z=13767
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn_n/(gam3[n]*fkv[n]);  //Z=13768
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=13769
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=13770

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
                    //i = 2;  //Z=13837
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
                        xrn_n = -xrn_n*xr2z;  //Z=13847
                        /*  longitudinal  */  //Z=13848
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13849
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13850
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
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=13868

                        /* carr4p[n]:=power(4,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(pi)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13870 */

                        sump = 0.0;  //Z=13872
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13873
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=13874

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
                    //i = 2;  //Z=13893
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
                        xrn_n = -xrn_n*xr2z;  //Z=13903
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13904
                        pn[n] = pn[n-1]*p*p;  //Z=13905
                        /*  longitudinal  */  //Z=13906
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13907
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13908
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
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13928
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
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=13945
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
                    //i = 2;  //Z=13984
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
                        xrn_n = -xrn_n*xr2z;  //Z=13994
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13995
                        pn[n] = pn[n-1]*p*p;  //Z=13996
                        /*  longitudinal  */  //Z=13997
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=13998
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13999
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
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn_n/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14017
                        sump = 0.0;  //Z=14018
                        sump1 = 0.0;  //Z=14019
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14020
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=14021
                            sump = sump+pn[n-m]*sumi;  //Z=14022
                            sump1 = sump1+sumi;  //Z=14023
                        }/*7*/  //Z=14024
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=14025
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn_n*sump1;  //Z=14026
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
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn_n*sump;  //Z=14036

                        /*  F(q)-coefficients  */  //Z=14051
                        /*  cross-sectional  */  //Z=14052
                        params.CR->carr4f[n] = z12v[n]*xrn_n/((n+1)*fkv[n]*fkv[n]);  //Z=14053
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14054
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14055

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
                    //i = 2;  //Z=14105
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
                        xrn_n = -xrn_n*x12zm;  //Z=14115
                        /*  longitudinal  */  //Z=14116
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14117
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14118
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
                        params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=14136

                        sump = 0.0;  //Z=14138
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14139
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=14140

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
                    //i = 2;  //Z=14164
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14165
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14166
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14167
                        fkv[n] = fkv[n-1]*n;  //Z=14168
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14169
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14170
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14171 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14172
                        xrn_n = -xrn_n*xr2z;  //Z=14173
                        /*  longitudinal  */  //Z=14174
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14175
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14176
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
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14195
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
                    //i = 2;  //Z=14227
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14228
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14229
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14230
                        fkv[n] = fkv[n-1]*n;  //Z=14231
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14232
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14233
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14234 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14235
                        xrn_n = -xrn_n*xr2z;  //Z=14236
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14237
                        pn[n] = pn[n-1]*p*p;  //Z=14238
                        /*  longitudinal  */  //Z=14239
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14240
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14241
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
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14262
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
                        params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=14282
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
                    //i = 2;  //Z=14329
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14330
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14331
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14332
                        fkv[n] = fkv[n-1]*n;  //Z=14333
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14334
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14335
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14336 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14337
                        xrn_n = -xrn_n*xr2z;  //Z=14338
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14339
                        pn[n] = pn[n-1]*p*p;  //Z=14340
                        /*  longitudinal  */  //Z=14341
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14342
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=14343
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
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn_n/((n+1)*gam3[n]*fkv[n]);  //Z=14363
                        sump = 0.0;  //Z=14364
                        sump1 = 0.0;  //Z=14365
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14366
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=14367
                            sump = sump+pn[n-m]*sumi;  //Z=14368
                            sump1 = sump1+sumi;  //Z=14369
                        }/*7*/  //Z=14370
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14371
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn_n*sump1;  //Z=14372
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
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn_n*sump;  //Z=14382

                        /*  F(q)  */  //Z=14401
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn_n/(gam3[n]*fkv[n]);  //Z=14402
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14403
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14404

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
                    //i = 2;  //Z=14459
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
                        xrn_n = -xrn_n*xr2z;  //Z=14469
                        /*  longitudinal  */  //Z=14470
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14471
                        /*  P(q)-coefficient  */  //Z=14472
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14473
                        /*  F(q)-coefficient  */  //Z=14474
                        sump = 0.0;  //Z=14475
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14476
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n)*fkv[n]*fkv[n]*params.norm);  //Z=14477
                        /*  cross-sectional  */  //Z=14478
                        /*  P(q)-coefficient  */  //Z=14479
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=14480
                        /*  F(q)-coefficient  */  //Z=14481
                        sump = 0.0;  //Z=14482
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14483
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=14484

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
                    //i = 2;  //Z=14503
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
                        xrn_n = -xrn_n*xr2z;  //Z=14513
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14514
                        pn[n] = pn[n-1]*p*p;  //Z=14515
                        /*  longitudinal  */  //Z=14516
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14517
                        /*  P(q)-coefficient  */  //Z=14518
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14519
                        /*  F(q)-coefficient  */  //Z=14520
                        sump = 0.0;  //Z=14521
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14522
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14523
                        /*  P(q)-coefficients  */  //Z=14524
                        /*  cross-sectional  */  //Z=14525
                        /*  F121  */  //Z=14526
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14527
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
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=14544
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
                    //i = 2;  //Z=14584
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
                        xrn_n = -xrn_n*xr2z;  //Z=14594
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14595
                        pn[n] = pn[n-1]*p*p;  //Z=14596
                        /*  longitudinal  */  //Z=14597
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14598
                        /*  P(q)-coefficient  */  //Z=14599
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14600
                        /*  F(q)-coefficient  */  //Z=14601
                        sump = 0.0;  //Z=14602
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14603
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14604

                        /*  cross-sectional P(q)  */  //Z=14606
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn_n/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14607
                        sump = 0.0;  //Z=14608
                        sump1 = 0.0;  //Z=14609
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14610
                            sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=14611
                            sump = sump+pn[n-m]*sumi;  //Z=14612
                            sump1 = sump1+sumi;  //Z=14613
                        }/*7*/  //Z=14614
                        params.CR->carr5p[n] = (1-a/2.0)*z12v[n]*xrmn_n*sump;  //Z=14615
                        params.CR->carr6p[n] = (1-a/2.0)*z12v[n]*xrn_n*sump1;  //Z=14616
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
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v[n]*xrn_n*sump;  //Z=14626

                        /*  F(q)-coefficients  */  //Z=14641
                        /*  cross-sectional  */  //Z=14642
                        params.CR->carr4f[n] = z12v[n]*xrn_n/((n+1)*fkv[n]*fkv[n]);  //Z=14643
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14644
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14645

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
                    //i = 2;  //Z=14685
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
                        xrn_n = -xrn_n*x12zm;  //Z=14695
                        /*  longitudinal  */  //Z=14696
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14697
                        /*  P(q)-coefficient  */  //Z=14698
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14699
                        /*  F(q)-coefficient  */  //Z=14700
                        sump = 0.0;  //Z=14701
                        for ( m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14702
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n)*fkv[n]*fkv[n]*params.norm);  //Z=14703
                        /*  cross-sectional  */  //Z=14704
                        /*  P(q)-coefficient  */  //Z=14705
                        params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=14706

                        /*  F(q)-coefficient  */  //Z=14708
                        sump = 0.0;  //Z=14709
                        for ( m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14710
                        params.CR->carr4f[n] = xrn_n*sump;  //Z=14711

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
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14738
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14739
                    }/*6*/  //Z=14740
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14741
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14742
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14743
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14744
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14745 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14746
                        xrn_n = -xrn_n*xr2z;  //Z=14747
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
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14761
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
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14798
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14799
                    }/*6*/  //Z=14800
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14801
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14802
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14803
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14804
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14805 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14806
                        xrn_n = -xrn_n*xr2z;  //Z=14807
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
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14824
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
                        params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=14844
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
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14895
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14896
                    }/*6*/  //Z=14897
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14898
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14899
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14900
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14901
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14902 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14903
                        xrn_n = -xrn_n*xr2z;  //Z=14904
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
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn_n/((n+1)*gam3[n]*fkv[n]);  //Z=14920
                        sump = 0.0;  //Z=14921
                        sump1 = 0.0;  //Z=14922
                        for ( m=0; m<=n; m++ )
                        {/*7*/  //Z=14923
                            sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=14924
                            sump = sump+pn[n-m]*sumi;  //Z=14925
                            sump1 = sump1+sumi;  //Z=14926
                        }/*7*/  //Z=14927
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14928
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn_n*sump1;  //Z=14929
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
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn_n*sump;  //Z=14939

                        /*  F(q)  */  //Z=14959
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn_n/(gam3[n]*fkv[n]);  //Z=14960
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14961
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14962

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
    // Für Python-Tests hier mal nur das carr4p[] ausgeben, aber komplett
    QString str = "carr4p[]";
    str += QString("%1").arg(params.CR->carr4p[0]);
    for ( int n=1; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
        str += QString(", %1").arg(params.CR->carr4p[n]);
    qDebug() << str;
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



void SasCalc_GENERIC_calculation::coefficients(int dim/*cdim*/, int nmax/*120*/, double &order/*order*/)
{/*1*/  //Z=9363
    // aus dem File 20240301 - crystal3d1.pas mit den dortigen Zeilennummern !!!

    const int nnmax = 120; // 500;  //Z=9367
    const double min =  1E-30;  //Z=9368

    int  n, n1, n2, n3, n4, n5, n6, n7, n8, n9, n1f, n2f, n3f, n4f, n5f, n6f, n7f, n8f, n9f, cnt, inmax;  //Z=9372
    double z, zz, zl, zzl, xlz, xl2z, xrz, xr2z, b1s, binsum;  //Z=9374
    double a1, intl, p, vvm, eps, area, vol;  //Z=9375
    double xrmz, xrm2z, sump, sump1, sump2, sump3, radv, aell, bell, cell;  //Z=9376
    // TODO: ab sump noch weiter optimieren (und auch die s=s+...)
    double xrn[nnmax+1], e1[nnmax+1], carr2i[coeffarray_len], fsum[nnmax+1];  //Z=9377
    double lqm[10], qqm[10], phim[3*nnmax+1], llm[3*nnmax+1], rrm[3*nnmax+1], ppm[3*nnmax+1], a1m[3*nnmax+1], u1ell[nnmax+1],
           u2ell[nnmax+1], v1ell[3*nnmax+2], v2ell[3*nnmax+2], v3ell[3*nnmax+2], gell[3*nnmax+2], g3ell[3*nnmax+2];  //Z=9378
    double intlar[nnmax+1][nnmax+1];  //Z=9379
    bool search1, /*search2,*/ search3, search4, search5, search6; //, search7, search8, search9;  //Z=9380
    float x1zm, x12zm=1.0;  //Z=9382

    // Arrays von oben mit erhöhter Genauigkeit
    long double gam3[2*nnmax+2], fkv[2*nnmax+2], fk2v[nnmax+1], fkvm[nnmax+1], pn[nnmax+1], xln[nnmax+1], z12vl[nnmax+1];

    // Das Array gam3[] scheint größer sein zu müssen. Im Pascal-Code ist es bis 500 deklariert (nmax=nnmax=120)
    typedef struct { int max; QString dbg; } gam3s_struct;
#define GAM3S_INIT(m,s)  gam3s_struct m; m.max=0; m.dbg=s;
#define GAM3S_CHECK(m,i) if ( (i) > m.max ) m.max = (i);
#define GAM3S_PRINT(m)   if ( m.max > 0 ) qDebug() << "coeeficients() gam3[" << m.dbg << "]" << m.max << "nmax" << nmax;

    // TODO: Vielleicht sollten einige Arrays nur lokal angelegt werden, sobald diese auch verwendet werden.

    // Sicherheitsabfrage: nmax (Parameter) wird für Schleifen genutzt (i.d.R. 120),
    //                     nnmax (const hier) wird für die Arrays genutzt (z.Zt. 120)
    if ( nnmax+1 < nmax )
    {
#ifndef __CUDACC__
        qDebug() << "coefficients() FATAL error: nmax" << nmax << "> nnmax" << nnmax;
#endif
        exit(255);
    }

    /*begin*/  //Z=9385
    z = (1.0-sqr(params.sigma))/sqr(params.sigma);  //Z=9386
    zl = (1.0-sqr(params.sigmal))/sqr(params.sigmal);  //Z=9387
    zz = z;  //Z=9388
    zzl = zl;  //Z=9389
    /* z:=z+2*dim;          (* z-average *)  //Z=9390 */
    /* zl:=zl+2*dim;  //Z=9391 */
    //nmax2 = ceil(nmax/2.0);  //Z=9392
    p = params.radius/params.radiusi;  //Z=9393  r/rm
    eps = params.length/params.radius;  //Z=9394  l/r

    xlz = params.length/(2.0*(zzl+1));  //Z=9396
    xl2z = xlz*xlz;  //Z=9397
    xrz = params.radius/(2.0*(zz+1));  //Z=9398
    xr2z = xrz*xrz;  //Z=9399
    xrmz = params.radiusi/(2.0*(zz+1));  //Z=9400
    xrm2z = xrmz*xrmz;  //Z=9401

    /*  ellipsoid semiaxis  */  //Z=9403
    aell = params.radius;  //Z=9404
    bell = params.length;  //Z=9405
    cell = params.radiusi;  //Z=9406

    xln[0] = 1;  //Z=9408
    xrn[0] = 1;  //Z=9409
    //xrmn[0] = 1;  //Z=9410
    pn[0] = 1;  //Z=9411
    //z12v[0] = 1;  //Z=9412
    z12vl[0] = 1;  //Z=9413
    fsum[0] = 1;  //Z=9414


    /*  factor for cross-sectional formfactors  */  //Z=9417
    if ( dim==1 ) b1s = 2;         /*  cylinder, cross-section is disk  */  //Z=9418
    if ( dim==2 ) b1s = 3/2.0;       /*  disk, cross-section is cylinder  */  //Z=9419
    if ( dim==3 ) b1s = 5/2.0;       /*  sphere  */  //Z=9420
    if ( dim==4 ) b1s = 3/2.0;       /*  cube  */  //Z=9421

    /*  start values for recursive parameters  */  //Z=9423
    //b1sv[0] = 1;  //Z=9424
    fkv[0] = 1;  //Z=9425
    fk2v[0] = 1;  //Z=9426
    e1[0] = 1;  //Z=9427
    //gam1[0] = sqrt(M_PI);  //Z=9428
    //gam2[0] = 1.0;  //Z=9429
    gam3[0] = sqrt(M_PI)/2.0;  //Z=9430

    /*  initialize  */  //Z=9432
    /* for n:=0 to 10000 do begin  //Z=9433 */
    /*    carr1pm[n]:=1;  //Z=9434 */
    /*    carr2pm[n]:=1;  //Z=9435 */
    /* end;  //Z=9436 */
    for ( int n=0; n<coeffarray_len; n++ )
    {   //Z=9437
        params.CR->carr1p[n] = 1;      params.CR->carr1f[n] = 1;  //Z=9438
        params.CR->carr2p[n] = 1;      params.CR->carr2f[n] = 1;  //Z=9439
        params.CR->carr3p[n] = 1;      params.CR->carr3f[n] = 1;  //Z=9440
        params.CR->carr4p[n] = 1;      params.CR->carr4f[n] = 1;  //Z=9441
        params.CR->carr5p[n] = 1;      params.CR->carr5f[n] = 1;  //Z=9442
        params.CR->carr6p[n] = 1;      params.CR->carr6f[n] = 1;  //Z=9443
        params.CR->carr7p[n] = 1;      params.CR->carr7f[n] = 1;  //Z=9444
        params.CR->carr8p[n] = 1;      params.CR->carr8f[n] = 1;  //Z=9445
        params.CR->carr9p[n] = 1;      params.CR->carr9f[n] = 1;  //Z=9446
        carr2i[n] = 1;  //Z=9447
    }   //Z=9448
    for ( int n=0; n<imax2d_len; n++ )
    {   //Z=9449
        for ( int m=0; m<imax2d_len; m++ )
        {   //Z=9450
            params.CR->carr11pm[n][m] = 1;  //Z=9451
            params.CR->carr11pm[n][m] = 1;  //Z=9452
        }   //Z=9453
    }   //Z=9454

    /*  multi-shell or liposome structure parameters  */  //Z=9456
    if ( params.cs==3 )
    {/*2*/  //Z=9457
        const double philiph = params.CR->myarray[12];     /*  water  */  //Z=9458
        const double philipt = params.CR->myarray[13];     /*  bilayer  */  //Z=9459

        const double rad = params.CR->myarray[1];          /*  vesicle inner radius  */  //Z=9461
        const double lliph = params.CR->myarray[7];        /*  water  */  //Z=9462
        const double llipt = params.CR->myarray[8];        /*  bilayer  */  //Z=9463

        const double len = lliph+llipt;  //Z=9465
        const int    ncell = round(params.CR->myarray[4]);  //Z=9466
        const double rmax = rad+ncell*len;  //Z=9467

        lqm[1] = lliph;  //Z=9469
        lqm[2] = lqm[1]+llipt;  //Z=9470

        for ( int i=1; i<=2; i++ ) qqm[i] = lqm[i]/len;  //Z=9472

        phim[1] = philiph;       /*  vesicle interior  */  //Z=9474
        llm[1] = 0.0;  //Z=9475
        rrm[1] = rad+llm[1];  //Z=9476
        ppm[1] = 1.0;  //Z=9477

        radv = rrm[1];        /*  vesicle radius  */  //Z=9479
        cnt = 1;  //Z=9480
        for ( int i=1; i<=ncell; i++ )
        {/*3*/  //Z=9481
            phim[cnt+1] = philipt;         /*  bilayer  */  //Z=9482
            phim[cnt+2] = philiph;         /*  water  */  //Z=9483
            llm[cnt+1] = (i-1+qqm[1])*len;  //Z=9484
            llm[cnt+2] = (i-1+qqm[2])*len;  //Z=9485
            rrm[cnt+1] = radv+llm[cnt+1];  //Z=9486
            rrm[cnt+2] = radv+llm[cnt+2];  //Z=9487
            ppm[cnt+1] = rrm[cnt+1]/rad;  //Z=9488
            ppm[cnt+2] = rrm[cnt+2]/rad;  //Z=9489
            cnt = cnt+2;  //Z=9490
        }/*3*/  //Z=9491
        inmax = cnt;  //Z=9492
        phim[cnt+1] = 0.0;  //Z=9493

        //xradm = rad;  //Z=9495
        //xrzm = xradm/(z+1);  //Z=9496
        x1zm = rad/(2.0*(z+1));  //Z=9497
        x12zm = x1zm*x1zm;  //Z=9498
        //xmax = q*rmax;  //Z=9499

        for ( int i=1; i<=inmax; i++ )
        {/*3*/  //Z=9501
            if ( params.part==0 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],3); /*  spheres  */  //Z=9502
            if ( params.part==1 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],2); /*  cylinders  */  //Z=9503
            if ( params.part==2 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],1); /*  disks  */  //Z=9504
            params.CR->carr7p[i] = ppm[i];  //Z=9505
            params.CR->carr3p[i] = llm[i];  //Z=9506
            params.CR->carr5p[i] = a1m[i];  //Z=9507
        }/*3*/  //Z=9508

        fkvm[0] = 1;  //Z=9510
        //for ( int n=1; n<=nmax; n++ )
        //{   //Z=9512
        if ( params.part==0 )
        {       /*  spheres  */  //Z=9513
            double gam3_n = sqrt(M_PI)/2.0;  //Z=9511
            for ( int n=1; n<=nmax; n++ )
            {
                fkvm[n] = fkvm[n-1]*n;  //Z=9514
                gam3_n = gam3_n*(2*n+1)/2.0;  //Z=9515
                params.CR->carr6p[n] = (n+3/2.0)*gam3_n*fkvm[n]*4/(3.0*sqrt(M_PI));  //Z=9516
            }
        }   //Z=9517
        else if ( params.part==1 )
        {       /*  cylinders  */  //Z=9518
            for ( int n=1; n<=nmax; n++ )
            {
                fkvm[n] = fkvm[n-1]*n;  //Z=9519
                params.CR->carr6p[n] = (n+1)*fkvm[n]*fkvm[n];  //Z=9520
            }
        }   //Z=9521
        else if ( params.part==2 )
        {       /*  disks  */  //Z=9522
            double gam3_n = sqrt(M_PI)/2.0;  //Z=9511
            for ( int n=1; n<=nmax; n++ )
            {
                fkvm[n] = fkvm[n-1]*n;  //Z=9523
                gam3_n = gam3_n*(2*n+1)/2.0;  //Z=9524
                params.CR->carr6p[n] = gam3_n*fkvm[n]*2/sqrt(M_PI);  //Z=9525
            }
        }   //Z=9526
        //}   //Z=9527

        vvm = 0;  //Z=9529
        for ( int i=1; i<=inmax; i++ )
        {   //Z=9530
            for ( int j=1; j<=inmax; j++ ) vvm = vvm+a1m[i]*a1m[j];  //Z=9531
        }   //Z=9532

        params.CR->myarray[14] = inmax;  //Z=9534
        params.CR->myarray[15] = vvm;  //Z=9535
        params.CR->myarray[16] = rmax;  //Z=9536
    }/*2*/  //Z=9537


    /*  myelin structure parameters  */  //Z=9540
    if ( params.cs==4 )
    {/*2*/  //Z=9541
        const double philiph = params.CR->myarray[12];     /*  head group  */  //Z=9542
        const double philipt = params.CR->myarray[13];     /*  tail  */  //Z=9543
        const double phiax = params.CR->myarray[9];        /*  axon  */  //Z=9544
        const double phiin = params.CR->myarray[10];       /*  intra cell  */  //Z=9545
        const double phiout = params.CR->myarray[11];      /*  extra cell  */  //Z=9546

        const double rad = params.CR->myarray[1];          /*  vesicle inner radius  */  //Z=9548
        const double lliph = params.CR->myarray[7];        /*  head group  */  //Z=9549
        const double llipt = params.CR->myarray[8];        /*  tail  */  //Z=9550
        const double lin = params.CR->myarray[6];          /*  intra cell  */  //Z=9551
        const double lout = params.CR->myarray[5];         /*  extra cell  */  //Z=9552

        const double len = lout+2*(2*lliph+llipt)+lin;  //Z=9554
        const int    ncell = round(params.CR->myarray[4]);  //Z=9555
        const double rmax = rad+ncell*len;  //Z=9556

        lqm[1] = lout;  //Z=9558
        lqm[2] = lqm[1]+lliph;  //Z=9559
        lqm[3] = lqm[2]+llipt;  //Z=9560
        lqm[4] = lqm[3]+lliph;  //Z=9561
        lqm[5] = lqm[4]+lin;  //Z=9562
        lqm[6] = lqm[5]+lliph;  //Z=9563
        lqm[7] = lqm[6]+llipt;  //Z=9564
        lqm[8] = lqm[7]+lliph;  //Z=9565

        for ( int i=1; i<=8; i++ ) qqm[i] = lqm[i]/len;  //Z=9567

        phim[1] = phiax;       /*  vesicle interior  */  //Z=9569
        llm[1] = 0.0;  //Z=9570
        rrm[1] = rad+llm[1];  //Z=9571
        ppm[1] = 1.0;  //Z=9572
        phim[2] = philiph;     /*  vesicle bilayer: head group  */  //Z=9573
        llm[2] = lliph;  //Z=9574
        rrm[2] = rad+llm[2];  //Z=9575
        ppm[2] = rrm[2]/rad;  //Z=9576
        phim[3] = philipt;     /*  vesicle bilayer: tail group  */  //Z=9577
        llm[3] = llm[2]+llipt;  //Z=9578
        rrm[3] = rad+llm[3];  //Z=9579
        ppm[3] = rrm[3]/rad;  //Z=9580
        phim[4] = philiph;     /*  vesicle bilayer: head group  */  //Z=9581
        llm[4] = llm[3]+lliph;  //Z=9582
        rrm[4] = rad+llm[4];  //Z=9583
        ppm[4] = rrm[4]/rad;  //Z=9584

        radv = rrm[4];        /*  vesicle radius + bilayer  */  //Z=9586
        cnt = 4;  //Z=9587
        for ( int i=1; i<=ncell; i++ )
        {/*3*/  //Z=9588
            phim[cnt+1] = phiout;          /*  extra cell  */  //Z=9589
            phim[cnt+2] = philiph;         /*  head group  */  //Z=9590
            phim[cnt+3] = philipt;         /*  tail group  */  //Z=9591
            phim[cnt+4] = philiph;         /*  head group  */  //Z=9592
            phim[cnt+5] = phiin;           /*  intra cell  */  //Z=9593
            phim[cnt+6] = philiph;         /*  head group  */  //Z=9594
            phim[cnt+7] = philipt;         /*  tail group  */  //Z=9595
            phim[cnt+8] = philiph;         /*  head group  */  //Z=9596
            llm[cnt+1] = (i-1+qqm[1])*len;  //Z=9597
            llm[cnt+2] = (i-1+qqm[2])*len;  //Z=9598
            llm[cnt+3] = (i-1+qqm[3])*len;  //Z=9599
            llm[cnt+4] = (i-1+qqm[4])*len;  //Z=9600
            llm[cnt+5] = (i-1+qqm[5])*len;  //Z=9601
            llm[cnt+6] = (i-1+qqm[6])*len;  //Z=9602
            llm[cnt+7] = (i-1+qqm[7])*len;  //Z=9603
            llm[cnt+8] = (i-1+qqm[8])*len;  //Z=9604
            rrm[cnt+1] = radv+llm[cnt+1];  //Z=9605
            rrm[cnt+2] = radv+llm[cnt+2];  //Z=9606
            rrm[cnt+3] = radv+llm[cnt+3];  //Z=9607
            rrm[cnt+4] = radv+llm[cnt+4];  //Z=9608
            rrm[cnt+5] = radv+llm[cnt+5];  //Z=9609
            rrm[cnt+6] = radv+llm[cnt+6];  //Z=9610
            rrm[cnt+7] = radv+llm[cnt+7];  //Z=9611
            rrm[cnt+8] = radv+llm[cnt+8];  //Z=9612
            ppm[cnt+1] = rrm[cnt+1]/rad;  //Z=9613
            ppm[cnt+2] = rrm[cnt+2]/rad;  //Z=9614
            ppm[cnt+3] = rrm[cnt+3]/rad;  //Z=9615
            ppm[cnt+4] = rrm[cnt+4]/rad;  //Z=9616
            ppm[cnt+5] = rrm[cnt+5]/rad;  //Z=9617
            ppm[cnt+6] = rrm[cnt+6]/rad;  //Z=9618
            ppm[cnt+7] = rrm[cnt+7]/rad;  //Z=9619
            ppm[cnt+8] = rrm[cnt+8]/rad;  //Z=9620
            cnt = cnt+8;  //Z=9621
        }/*3*/  //Z=9622
        inmax = cnt;  //Z=9623
        phim[cnt+1] = 0.0;  //Z=9624

        //xradm = rad;  //Z=9626
        //xrzm = xradm/(z+1);  //Z=9627
        x1zm = rad/(2.0*(z+1));  //Z=9628
        x12zm = x1zm*x1zm;  //Z=9629
        //xmax = q*rmax;  //Z=9630

        for ( int i=1; i<=inmax; i++ )
        {/*3*/  //Z=9632
            if ( params.part==0 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],3); /*  spheres  */  //Z=9633
            if ( params.part==1 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],2); /*  cylinders  */  //Z=9634
            if ( params.part==2 ) a1m[i] = (phim[i]-phim[i+1])*pow(rrm[i],1); /*  disks  */  //Z=9635
            params.CR->carr7p[i] = ppm[i];  //Z=9636
            params.CR->carr3p[i] = llm[i];  //Z=9637
            params.CR->carr5p[i] = a1m[i];  //Z=9638
        }/*3*/  //Z=9639

        fkvm[0] = 1;  //Z=9641
        //for ( int n=1; n<=nmax; n++ )
        //{   //Z=9643
        if ( params.part==0 )
        {/*4*/       /*  spheres  */  //Z=9644
            double gam3_n = sqrt(M_PI)/2.0;  //Z=9511
            for ( int n=1; n<=nmax; n++ )
            {
                fkvm[n] = fkvm[n-1]*n;  //Z=9645
                gam3_n = gam3_n*(2*n+1)/2.0;  //Z=9646
                params.CR->carr6p[n] = (n+3/2.0)*gam3_n*fkvm[n]*4/(3.0*sqrt(M_PI));  //Z=9647
            }
        }/*4*/  //Z=9648
        else if ( params.part==1 )
        {/*4*/       /*  cylinders  */  //Z=9649
            for ( int n=1; n<=nmax; n++ )
            {
                fkvm[n] = fkvm[n-1]*n;  //Z=9650
                params.CR->carr6p[n] = (n+1)*fkvm[n]*fkvm[n];  //Z=9651
            }
        }/*4*/  //Z=9652
        else if ( params.part==2 )
        {/*4*/       /*  disks  */  //Z=9653
            double gam3_n = sqrt(M_PI)/2.0;  //Z=9511
            for ( int n=1; n<=nmax; n++ )
            {
                fkvm[n] = fkvm[n-1]*n;  //Z=9654
                gam3_n = gam3_n*(2*n+1)/2.0;  //Z=9655
                params.CR->carr6p[n] = gam3_n*fkvm[n]*2/sqrt(M_PI);  //Z=9656
            }
        }/*4*/  //Z=9657
        //}   //Z=9658

        vvm = 0;  //Z=9660
        for ( int i=1; i<=inmax; i++ )
        {   //Z=9661
            for ( int j=1; j<=inmax; j++ ) vvm = vvm+a1m[i]*a1m[j];  //Z=9662
        }   //Z=9663

        params.CR->myarray[14] = inmax;  //Z=9665
        params.CR->myarray[15] = vvm;  //Z=9666
        params.CR->myarray[16] = rmax;  //Z=9667
    }/*2*/  //Z=9668




    search1 = true;  //Z=9673
    //search2 = true;  //Z=9674
    search3 = true;  //Z=9675
    search4 = true;  //Z=9676
    search5 = true;  //Z=9677
    search6 = true;  //Z=9678
    //search7 = true;  //Z=9679
    //search8 = true;  //Z=9680
    //search9 = true;  //Z=9681
    n1 = nmax;      n1f = nmax;  //Z=9682
    n2 = nmax;      n2f = nmax;  //Z=9683
    n3 = nmax;      n3f = nmax;  //Z=9684
    n4 = nmax;      n4f = nmax;  //Z=9685
    n5 = nmax;      n5f = nmax;  //Z=9686
    n6 = nmax;      n6f = nmax;  //Z=9687
    n7 = nmax;      n7f = nmax;  //Z=9688
    n8 = nmax;      n8f = nmax;  //Z=9689
    n9 = nmax;      n9f = nmax;  //Z=9690

    /*  cho1 = orientation case  */  //Z=9692
    params.orcase = 1;                                       /*  general  */  //Z=9693
    if ( (params.polPhi== 0) && (params.polTheta==90) ) params.orcase = 2;      /*  x-axis  */  //Z=9694
    if ( (params.polPhi==90) && (params.polTheta==90) ) params.orcase = 3;     /*  y-axis  */  //Z=9695
    if ( (params.polPhi== 0) && (params.polTheta== 0) ) params.orcase = 4;       /*  z-axis  */  //Z=9696
    if ( (params.polPhi==90) && (params.polTheta== 0) ) params.orcase = 4;      /*  z-axis  */  //Z=9697


    /* ** isotropic case for spheres ** */  //Z=9700
    if ( dim==3 )
    {/*2*/  //Z=9701
        params.norm = 1;  //Z=9702
        order = 0;  //Z=9703
        /*  homogeneous  */  //Z=9704
        if ( params.cs==0 )
        {/*3*/  //Z=9705
            double xrn_n = 1.0;
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=9706
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=9707
                fkv[n] = fkv[n-1]*n;  //Z=9708
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=9709
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=9710 */
                xrn_n = -xrn_n*xr2z;  //Z=9711
                /*  P(q)-coefficient  */  //Z=9712
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=9713
                /*  F(q)-coefficient  */  //Z=9714
                double binsum = 0.0;  //Z=9715
                for ( int m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=9716
                params.CR->carr4f[n] = 9*M_PI*xrn_n*binsum/16.0;  //Z=9717
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=9718
                    if ( n<n4 ) n4 = n;  //Z=9719
                }/*5*/  //Z=9720
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=9721
                    if ( n<n4f ) n4f = n;  //Z=9722
                }/*5*/  //Z=9723
            }/*4*/  //Z=9724
            goto Label99;  //Z=9725
        }/*3*/   /*  of homogeneous  */  //Z=9726

        /*  core/shell  */  //Z=9728
        if ( params.cs==1 )
        {/*3*/  //Z=9729
            double xrn_n = 1.0;
            double xrmn_n = 1.0;
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=9730
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=9731
                fkv[n] = fkv[n-1]*n;  //Z=9732
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=9733
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=9734 */
                xrn_n = -xrn_n*xr2z;  //Z=9735
                xrmn_n = -xrmn_n*xrm2z;  //Z=9736
                pn[n] = pn[n-1]*p*p;  //Z=9737
                /*  P(q)-coefficients  */  //Z=9738
                double sump = 0.0;  //Z=9739
                for ( int m=0; m<=n; m++ )
                {   //Z=9740
                    sump += pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=9741
                }   //Z=9742
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=9743
                params.CR->carr5p[n] = (9*M_PI/16.0)*z12v[n]*xrmn_n*sump;  //Z=9744
                params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=9745
                /*  F(q)-coefficients  */  //Z=9746
                sump = 0.0;  //Z=9747
                for ( int m=0; m<=n; m++ )
                {   //Z=9748
                    sump += z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=9749
                }   //Z=9750
                params.CR->carr4f[n] = (9*M_PI/16.0)*xrn_n*sump;  //Z=9751
                sump = 0.0;  //Z=9752
                for ( int m=0; m<=n; m++ )
                {   //Z=9753
                    sump += pn[m]*z12v[m]*z12v[n-m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=9754
                }   //Z=9755
                params.CR->carr5f[n] = (9*M_PI/16.0)*xrmn_n*sump;  //Z=9756
                params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=9757
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=9758
                    if ( n<n4 ) n4 = n;  //Z=9759
                }/*5*/  //Z=9760
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=9761
                    if ( n<n5 ) n5 = n;  //Z=9762
                }/*5*/  //Z=9763
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=9764
                    if ( n<n6 ) n6 = n;  //Z=9765
                }/*5*/  //Z=9766
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=9767
                    if ( n<n4f ) n4f = n;  //Z=9768
                }/*5*/  //Z=9769
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=9770
                    if ( n<n5f ) n5f = n;  //Z=9771
                }/*5*/  //Z=9772
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=9773
                    if ( n<n6f ) n6f = n;  //Z=9774
                }/*5*/  //Z=9775
            }/*4*/  //Z=9776
            goto Label99;  //Z=9777
        }/*3*/   /*  of core/shell  */  //Z=9778

        /*  inhomogeneous core/shell  */  //Z=9780
        if ( params.cs==2 )
        {/*3*/  //Z=9781
            double xrmn_n = 1.0;
            double z12v_n = 1.0;
            double xrn_n  = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=9782
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=9783
                fkv[n] = fkv[n-1]*n;  //Z=9784
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=9785
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=9786 */
                xrn_n = -xrn_n*xr2z;  //Z=9787
                xrmn_n = -xrmn_n*xrm2z;  //Z=9788
                pn[n] = pn[n-1]*p*p;  //Z=9789
                /*  P(q)-coefficients  */  //Z=9790
                params.CR->carr1p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v_n*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=9791
                double sump = 0.0;  //Z=9792
                double sump1 = 0.0;  //Z=9793
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=9794
                    const double sumi = 1/((n-m+3/2.0)*gam3[n-m]*(m+3/2.0-params.alphash1/2.0)*gam3[m]*fkv[m]*fkv[n-m]);  //Z=9795
                    sump += pn[n-m]*sumi;  //Z=9796
                    sump1 += sumi;  //Z=9797
                }/*5*/  //Z=9798
                params.CR->carr2p[n] = (3*M_PI*(3-params.alphash1)/16.0)*z12v_n*xrmn_n*sump;  //Z=9799
                params.CR->carr3p[n] = (3*M_PI*(3-params.alphash1)/16.0)*z12v_n*xrn_n*sump1;  //Z=9800
                sump = 0.0;  //Z=9801
                sump1 = 0.0;  //Z=9802
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=9803
                    const double sumi = 1/((n-m+3/2.0-params.alphash1/2.0)*(m+3/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=9804
                    sump += sumi;  //Z=9805
                    sump1 += pn[n-m]*sumi;  //Z=9806
                }/*5*/  //Z=9807
                params.CR->carr4p[n] = ((3-params.alphash1)*(3-params.alphash1)*M_PI/16.0)*z12v_n*xrmn_n*sump;  //Z=9808
                params.CR->carr5p[n] = ((3-params.alphash1)*(3-params.alphash1)*M_PI/16.0)*z12v_n*xrmn_n*sump1;  //Z=9809
                params.CR->carr6p[n] = ((3-params.alphash1)*(3-params.alphash1)*M_PI/16.0)*z12v_n*xrn_n*sump;  //Z=9810

                /*  F(q)-coefficients  */  //Z=9812
                params.CR->carr4f[n] = (3*sqrt(M_PI)/4.0)*z12v_n*xrn_n/((n+3/2.0)*gam3[n]*fkv[n]);  //Z=9813
                params.CR->carr5f[n] = (sqrt(M_PI)*(3-params.alphash1)/4.0)*z12v_n*xrmn_n/((n+3/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=9814
                params.CR->carr6f[n] = (sqrt(M_PI)*(3-params.alphash1)/4.0)*z12v_n*xrn_n/((n+3/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=9815
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=9816
                    if ( n<n1 ) n1 = n;  //Z=9817
                }/*5*/  //Z=9818
                if ( fabs(params.CR->carr2p[n])<min )
                {/*5*/  //Z=9819
                    if ( n<n2 ) n2 = n;  //Z=9820
                }/*5*/  //Z=9821
                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=9822
                    if ( n<n3 ) n3 = n;  //Z=9823
                }/*5*/  //Z=9824
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=9825
                    if ( n<n4 ) n4 = n;  //Z=9826
                }/*5*/  //Z=9827
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=9828
                    if ( n<n5 ) n5 = n;  //Z=9829
                }/*5*/  //Z=9830
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=9831
                    if ( n<n6 ) n6 = n;  //Z=9832
                }/*5*/  //Z=9833
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=9834
                    if ( n<n4f ) n4f = n;  //Z=9835
                }/*5*/  //Z=9836
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=9837
                    if ( n<n5f ) n5f = n;  //Z=9838
                }/*5*/  //Z=9839
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=9840
                    if ( n<n6f ) n6f = n;  //Z=9841
                }/*5*/  //Z=9842
            }/*4*/  //Z=9843
            goto Label99;  //Z=9844
        }/*3*/   /*  of inhomogeneous core/shell  */  //Z=9845

        /*  myelin  */  //Z=9847
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=9848
            //i = 2;  //Z=9849
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            double xrn_n = 1.0;
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=9850
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=9851
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=9852
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=9853
                fkv[n] = fkv[n-1]*n;  //Z=9854
                //fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=9855
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=9856
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=9857 */
                xln[n] = -xln[n-1]*xl2z;  //Z=9858
                /* xrn[n]:=-xrn[n-1]*xr2z;  //Z=9859 */
                xrn_n = -xrn_n*x12zm;         /*  myelin radius  */  //Z=9860

                /*  P(q)  */  //Z=9862
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=9863
                    /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=9864 */
                    /* i:=i+1;  //Z=9865 */
                    params.CR->carr11pm[n][m] = (9*M_PI/16.0)*(1/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]));  //Z=9866
                }/*5*/  //Z=9867
                params.CR->carr4p[n] = z12v[n]*xrn_n;  //Z=9868
                /* carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=9869 */


                /*  F(q)  */  //Z=9872
                binsum = 0.0;  //Z=9873
                for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=9874
                params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=9875
                binsum = 0.0;  //Z=9876
                for ( int m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=9877
                params.CR->carr4f[n] = xrn_n*binsum;  //Z=9878


                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=9881
                    if ( n<n1 ) n1 = n;  //Z=9882
                }/*5*/  //Z=9883
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=9884
                    if ( n<n1f ) n1f = n;  //Z=9885
                }/*5*/  //Z=9886
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=9887
                    if ( n<n4 ) n4 = n;  //Z=9888
                }/*5*/  //Z=9889
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=9890
                    if ( n<n4f ) n4f = n;  //Z=9891
                }/*5*/  //Z=9892
            }/*4*/  //Z=9893
        }/*3*/ /*  of myelin  */  //Z=9894



    }/*2*/  /*  of dim=3, spheres  */  //Z=9898

    /*  isotropic case for cubes  */  //Z=9900
    if ( (params.ordis==7) && (dim==4) )
    {/*2*/    /*  cubes  */  //Z=9901
        params.norm = 1;  //Z=9902
        order = 0;  //Z=9903
        /*  homogeneous  */  //Z=9904
        if ( params.cs==0 )
        {/*3*/  //Z=9905
            area = 6*4*params.radius*params.radius;  //Z=9906
            vol = 8*params.radius*params.radius*params.radius;  //Z=9907
            params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=9908

            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            double xrn_n = 1.0;

            u1ell[0] = 2;  //Z=9910
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=9911
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=9912
                fkv[n] = fkv[n-1]*n;  //Z=9913
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=9914
                /* u1ell[n]:=z12v[n]/((n+1/2)*(n+1)*fkv[n]);  //Z=9915 */
                u1ell[n] = 1/((n+1/2.0)*(n+1)*fkv[n]);  //Z=9916
            }/*4*/  //Z=9917

            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=9919
                double sump1 = 0.0;  //Z=9920
                for ( int m=0; m<=n; m++ ) sump1 += u1ell[n-m]*u1ell[m];  //Z=9921
                v1ell[n] = sump1;  //Z=9922
            }/*4*/  //Z=9923

            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=9925
                xrn_n = -xrn_n*xr2z;  //Z=9926
                double sumf = 0.0;  //Z=9927
                for ( int m=0; m<=n; m++ ) sumf += z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=9928
                fsum[n] = sumf;  //Z=9929
                /*  P(q)-coefficient  */  //Z=9930
                double sump = 0.0;  //Z=9931
                /* for m:=0 to n do begin  //Z=9932
                   sump1:=0.0;  //Z=9933
                   for k:=0 to m do sump1:=sump1+z12v[m-k]*z12v[k]/((k+1/2)*(m-k+1/2)*(k+1)*fkv[k]*(m-k+1)*fkv[m-k]);  //Z=9934
                   sump:=sump+z12v[n-m]*sump1/((n-m+1/2)*(n-m+1)*fkv[n-m]);  //Z=9935
                end; */  //Z=9936

                /* for m:=0 to n do begin  //Z=9938
                   sump1:=0.0;  //Z=9939
                   for k:=0 to m do sump1:=sump1+u1ell[m-k]*u1ell[k];  //Z=9940
                   sump:=sump+u1ell[n-m]*sump1;  //Z=9941
                end; */  //Z=9942

                for ( int m=0; m<=n; m++ ) sump += u1ell[n-m]*v1ell[m];  //Z=9944

                /* carr4p[n]:=sqrt(pi)*power(4,n)*xrn[n]*sump/(16*gam3[n]);  //Z=9946 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*z12v[n]*xrn_n*sump/(16.0*gam3[n]);  //Z=9947

                /*  F(q)-coefficient  */  //Z=9949
                sump = 0.0;  //Z=9950
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=9951
                    double sump1 = 0.0;  //Z=9952
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=9953
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=9954
                }/*5*/  //Z=9955
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=9956
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=9957
                    if ( n<n4 ) n4 = n;  //Z=9958
                }/*5*/  //Z=9959
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=9960
                    if ( n<n4f ) n4f = n;  //Z=9961
                }/*5*/  //Z=9962
            }/*4*/  //Z=9963
            goto Label99;  //Z=9964
        }/*3*/   /*  of homogeneous  */  //Z=9965

        /*  core/shell  */          /*  not yet ready  */  //Z=9967
        if ( params.cs==1 )
        {/*3*/  //Z=9968
            double xrmn_n = 1.0;
            double z12v_n = 1.0;
            double xrn_n  = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=9969
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=9970
                fkv[n] = fkv[n-1]*n;  //Z=9971
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=9972
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=9973 */
                xrn_n = -xrn_n*xr2z;  //Z=9974
                xrmn_n = -xrmn_n*xrm2z;  //Z=9975
                pn[n] = pn[n-1]*p*p;  //Z=9976
                /*  P(q)-coefficient  */  //Z=9977
                double sump = 0.0;  //Z=9978
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=9979
                    sump += pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=9980
                }/*5*/  //Z=9981
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v_n*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=9982
                params.CR->carr5p[n] = (9*M_PI/16.0)*z12v_n*xrmn_n*sump;  //Z=9983
                params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=9984
                /*  F(q)-coefficient  */  //Z=9985
                /* carr3[n]:=3*sqrt(pi)*z12v[n]*xrn[n]/(4*(n+3/2)*gam3[n]*fkv[n]);  //Z=9986 */
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=9987
                    if ( n<n4 ) n4 = n;  //Z=9988
                }/*5*/  //Z=9989
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=9990
                    if ( n<n5 ) n5 = n;  //Z=9991
                }/*5*/  //Z=9992
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=9993
                    if ( n<n6 ) n6 = n;  //Z=9994
                }/*5*/  //Z=9995
            }/*4*/  //Z=9996
            goto Label99;  //Z=9997
        }/*3*/   /*  of core/shell  */  //Z=9998
    }/*2*/  /*  of dim=4, cubes  */  //Z=9999

    /*  perfect orientation case for cubes  */  //Z=10001
    if ( (params.ordis==6) && (dim==4) )
    {/*2*/    /*  cubes  */  //Z=10002
        params.norm = 1;  //Z=10003
        order = 1;  //Z=10004
        /*  homogeneous  */  //Z=10005
        if ( params.cs==0 )
        {/*3*/  //Z=10006
            //i = 2;  //Z=10007
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            double xrn_n = 1.0;
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=10008
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10009
                fkv[n] = fkv[n-1]*n;  //Z=10010
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10011
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10012 */
                xrn_n = -xrn_n*xr2z;  //Z=10013
                /*  P(q)-coefficient  */  //Z=10014
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10015
                    /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10016 */
                    /* carr1fm[i]:=carr1pm[i];  //Z=10017 */
                    params.CR->carr11pm[n][m] = (M_PI/4.0)*pow(4.0,n)*z12v[n-m]*z12v[m]*xrn_n/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10018
                    //i = i+1;  //Z=10019
                }/*5*/  //Z=10020
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn_n/(16.0*gam3[n]);  //Z=10021
                /*  F(q)-coefficient  */  //Z=10022
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n/(128.0*gam3[n]);  //Z=10023
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10024
                    if ( n<n4 ) n4 = n;  //Z=10025
                }/*5*/  //Z=10026
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10027
                    if ( n<n4f ) n4f = n;  //Z=10028
                }/*5*/  //Z=10029
            }/*4*/  //Z=10030
            goto Label99;  //Z=10031
        }/*3*/   /*  of homogeneous  */  //Z=10032

        /*  core/shell  */          /*  not yet ready  */  //Z=10034
        if ( params.cs==1 )
        {/*3*/  //Z=10035
            double xrmn_n = 1.0;
            double z12v_n = 1.0;
            double xrn_n  = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10036
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10037
                fkv[n] = fkv[n-1]*n;  //Z=10038
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10039
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10040 */
                xrn_n = -xrn_n*xr2z;  //Z=10041
                xrmn_n = -xrmn_n*xrm2z;  //Z=10042
                pn[n] = pn[n-1]*p*p;  //Z=10043
                /*  P(q)-coefficient  */  //Z=10044
                double sump = 0.0;  //Z=10045
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10046
                    sump += pn[m]/((m+3/2.0)*gam3[m]*(n-m+3/2.0)*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=10047
                }/*5*/  //Z=10048
                params.CR->carr4p[n] = 9*sqrt(M_PI)*pow(4.0,n)*z12v_n*xrn_n/(2.0*(n+3)*(n+2)*(n+3/2.0)*gam3[n]*fkv[n]);  //Z=10049
                params.CR->carr5p[n] = (9*M_PI/16.0)*z12v_n*xrmn_n*sump;  //Z=10050
                params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=10051
                /*  F(q)-coefficient  */  //Z=10052
                /* carr3[n]:=3*sqrt(pi)*z12v[n]*xrn[n]/(4*(n+3/2)*gam3[n]*fkv[n]);  //Z=10053 */
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10054
                    if ( n<n4 ) n4 = n;  //Z=10055
                }/*5*/  //Z=10056
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=10057
                    if ( n<n5 ) n5 = n;  //Z=10058
                }/*5*/  //Z=10059
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=10060
                    if ( n<n6 ) n6 = n;  //Z=10061
                }/*5*/  //Z=10062
            }/*4*/  //Z=10063
            goto Label99;  //Z=10064
        }/*3*/   /*  of core/shell  */  //Z=10065
    }/*2*/  /*  of dim=4, cubes  */  //Z=10066

    /*  isotropic case for ellipsoids  */  //Z=10068
    if ( (params.ordis==7) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10069
        params.norm = 1;  //Z=10070
        order = 0;  //Z=10071
        if ( eps==1 ) area = 4*M_PI*params.radius*params.radius;  //Z=10072
        if ( eps>1 ) area = 2*M_PI*params.radius*(params.radius+(sqr(params.length)/sqrt(sqr(params.length)-sqr(params.radius)))*asin(sqrt(sqr(params.length)-sqr(params.radius))/params.length));  //Z=10073
        if ( eps<1 ) area = 2*M_PI*params.radius*(params.radius+(sqr(params.length)/sqrt(sqr(params.radius)-sqr(params.length)))*asinh(sqrt(sqr(params.radius)-sqr(params.length))/params.length));  //Z=10074
        vol = (4*M_PI/3.0)*sqr(params.radius)*params.length;  //Z=10075
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10076

        /*  homogeneous  */  //Z=10078
        if ( params.cs==0 )
        {/*3*/  //Z=10079
            double z12v_n = 1.0;
            double xrn_n  = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10080
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10081
                fkv[n] = fkv[n-1]*n;  //Z=10082
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10083
                e1[n] = e1[n-1]*(eps*eps-1);  //Z=10084
                xrn_n = -xrn_n*xr2z;  //Z=10085
                double sump = 0.0;  //Z=10086
                for ( int m=0; m<=n; m++ ) sump += e1[n-m]/(fkv[n-m]*fkv[m]*(2*(n-m)+1));  //Z=10087
                /* sumf:=0.0;  //Z=10088 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10089 */
                /* fsum[n]:=sumf;  //Z=10090 */
                /*  P(q)-coefficient  */  //Z=10091
                params.CR->carr4p[n] = (9*sqrt(M_PI)/2.0)*pow(4.0,n)*z12v_n*xrn_n*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);  //Z=10092
                /*  F(q)-coefficient  */  //Z=10093
                sump = 0.0;  //Z=10094
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10095
                    sump1 = 0.0;  //Z=10096
                    for ( int k=0; k<=m; k++ ) sump1 = sump1+fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10097
                    sump = sump+sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10098
                }/*5*/  //Z=10099
                params.CR->carr4f[n]  = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=10100
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10101
                    if ( n<n4 ) n4 = n;  //Z=10102
                }/*5*/  //Z=10103
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10104
                    if ( n<n4f ) n4f = n;  //Z=10105
                }/*5*/  //Z=10106
            }/*4*/  //Z=10107
            goto Label99;  //Z=10108
        }/*3*/   /*  of homogeneous  */  //Z=10109

        /*  core/shell  */          /*  not yet ready  */  //Z=10111
        /* if cs=1 then begin  //Z=10112
             goto 99;  //Z=10141
           end;   (* of core/shell *) */  //Z=10142
    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=10143

    /*  perfect orientation case for ellipsoids  */  //Z=10145
    if ( (params.ordis==6) && (dim==5) )
    {   /*  ellipsoids  */  //Z=10146
        params.norm = 1;  //Z=10147
        order = 1;  //Z=10148
        /*  homogeneous  */  //Z=10149
        if ( params.cs==0 )
        {   //Z=10150
            double z12v_n = 1.0;
            double xrn_n  = 1.0;
            double z12vl_n = 1.0;
            double xln_n = 1.0;
            for ( n=0; n<=2*nmax+1; n++ )
            {/*4*/  //Z=10160
                fkv[n]  = fkv[n-1]*n;  //Z=10154
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10155  gam3[SIZE
            }
            GAM3S_INIT(max1,"ordis=6, cs=0: m-ll+n-k")
            GAM3S_INIT(max2,"ordis=6, cs=0: ll+k")
            for ( n=0; n<=nmax; n++ )
            {/*4*/  //Z=10160
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10152
                z12vl_n = z12vl_n*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10153
                xrn_n = -xrn_n*xr2z;  //Z=10157
                xln_n = -xln_n*xl2z;  //Z=10158

                a1 = sqr(3/4.0);  //Z=10161
                for ( int m=0; m<=nmax; m++ )
                {/*5*/  //Z=10162
                    double sump1 = 0.0;  //Z=10163
                    for ( int k=0; k<=n; k++ )
                    {/*6*/  //Z=10164
                        double sump2 = 0.0;  //Z=10165
                        for ( int ll=0; ll<=m; ll++ )  //Z=10166
                        {
                            GAM3S_CHECK( max1, m-ll+n-k )
                            GAM3S_CHECK( max2, ll+k )
                            sump2 += 1./(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3./2.0)*(ll+k+3./2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10167
                        }
                        sump1 += sump2/(fkv[n-k]*fkv[k]);  //Z=10168
                    }/*6*/  //Z=10169
                    params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10170
                }/*5*/  //Z=10171

                params.CR->carr4p[n] = a1*z12vl_n*xln_n;  //Z=10173
                params.CR->carr5p[n] = z12v_n*xrn_n;  //Z=10174

                /*  P(q)-coefficient  */  //Z=10176
                /* for m:=0 to n do begin  //Z=10177 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10178 */
                /* carr1fm[i]:=carr1pm[i];  //Z=10179 */
                /* i:=i+1;  //Z=10180 */
                /* end;  //Z=10181 */

                /*  F(q)-coefficient  */  //Z=10183
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n/(128.0*gam3[n]);  //Z=10184
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10185
                    if ( n<n4 ) n4 = n;  //Z=10186
                }/*5*/  //Z=10187
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10188
                    if ( n<n4f ) n4f = n;  //Z=10189
                }/*5*/  //Z=10190
            }/*4*/  //Z=10191
            GAM3S_PRINT(max1)
            GAM3S_PRINT(max2)
            goto Label99;  //Z=10192
        }   /*  of homogeneous  */  //Z=10193

        /*  core/shell  */          /*  not yet ready  */  //Z=10195
        /* if cs=1 then begin  //Z=10196
             goto 99;  //Z=10225
           end;   (* of core/shell *)   */  //Z=10226
    } /*  of dim=5, ellipsoid  */  //Z=10227

    /* ** ellipsoid orientational distribution  */  //Z=10229
    if ( (params.ordis==0) && (params.orcase==2) && (dim==5) )
    {/*2*/  //Z=10230

        if ( params.cs==0 )
        {/*3*/  //Z=10232
            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,2,0,0,0,0,params.CR->carr1p,params.norm);
            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,3,0,0,0,0,params.CR->carr1p,order);
            order = order/params.norm;  //Z=10235

            for ( int n=0; n<=2*nmax+1; n++ )
            {/*4*/  //Z=10247
                //z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10239
                //fkv[n] = fkv[n-1]*n;  //Z=10240
                //fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10241
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10242
                //xln[n] = -xln[n-1]*xl2z;  //Z=10243
            }
            double z12v_n = 1.0;
            double xrn_n  = 1.0;
            GAM3S_INIT(max1,"ordis=0, cs=0: m-ll+n-k")
            GAM3S_INIT(max2,"ordis=0, cs=0: ll+k")
            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=10247
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10238
                z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10239
                fkv[n] = fkv[n-1]*n;  //Z=10240
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10241
                xln[n] = -xln[n-1]*xl2z;  //Z=10243
                xrn_n = -xrn_n*xr2z;  //Z=10244

                for ( int m=0; m<=nmax; m++ )
                {/*5*/  //Z=10248
                    double sump1 = 0.0;  //Z=10249
                    for ( int k=0; k<=n; k++ )
                    {/*6*/  //Z=10250
                        double sump2 = 0.0;  //Z=10251
                        for ( int ll=0; ll<=m; ll++ )  //Z=10252
                        {
                            GAM3S_CHECK(max1,m-ll+n-k)
                            GAM3S_CHECK(max2,ll+k)
                            sump2 += 1.0/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3.0/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10253
                        }
                        sump1 += sump2/(fkv[n-k]*fkv[k]);  //Z=10254
                    }/*6*/  //Z=10255
                    params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10256
                }/*5*/  //Z=10257
                sump1 = 0.0;  //Z=10258
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10259
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                    sump1 += pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=10261
                    params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=10262
                }/*5*/  //Z=10263

                /*  all coefficients: Mok  */  //Z=10265
                /*  integral part  */  //Z=10266
                params.CR->carr4p[n] = sump1;  //Z=10267
                /*  qR-part  */  //Z=10268
                params.CR->carr5p[n] = z12v_n*xrn_n;  //Z=10269
                /*  qL-part  */  //Z=10270
                /* carr3p[n]:=sqr(3/4)*fk2v[n]*z12v[n]*xln[n]/power(4,n);  //Z=10271 */
                params.CR->carr6p[n] = sqr(3/4.0)*fk2v[n]*z12vl[n]*xln[n];  //Z=10272

                if ( search5 )
                {/*5*/  //Z=10289
                    if ( params.CR->carr5p[n]<1e-50 )
                    {/*6*/  //Z=10290
                        n5 = n;  //Z=10291
                        search5 = false;  //Z=10292
                    }/*6*/  //Z=10293
                }/*5*/  //Z=10294
                if ( search6 )
                {/*5*/  //Z=10295
                    if ( params.CR->carr6p[n]<1e-50 )
                    {/*6*/  //Z=10296
                        n6 = n;  //Z=10297
                        search6 = false;  //Z=10298
                    }/*6*/  //Z=10299
                }/*5*/  //Z=10300
                if ( search3 )
                {/*5*/  //Z=10301
                    if ( params.CR->carr3p[n]<1e-50 )
                    {/*6*/  //Z=10302
                        n3 = n;  //Z=10303
                        search3 = false;  //Z=10304
                    }/*6*/  //Z=10305
                }/*5*/  //Z=10306
                if ( search4 )
                {/*5*/  //Z=10307
                    if ( params.CR->carr4p[n]<1e-50 )
                    {/*6*/  //Z=10308
                        n4 = n;  //Z=10309
                        search4 = false;  //Z=10310
                    }/*6*/  //Z=10311
                }/*5*/  //Z=10312
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=10313
                    if ( n<n1f ) n1f = n;  //Z=10314
                }/*5*/  //Z=10315
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10316
                    if ( n<n4f ) n4f = n;  //Z=10317
                }/*5*/  //Z=10318
            }/*4*/  /*  of n-loop  */  //Z=10319
            GAM3S_PRINT(max1)
            GAM3S_PRINT(max2)
            goto Label99;  //Z=10321
        }/*3*/  /*  of cs-loop  */  //Z=10322
    }/*2*/  /*  of ordis-loop  */  //Z=10323

    /*  isotropic case for triaxial ellipsoids  */  //Z=10326
    if ( (params.ordis==7) && (dim==6) )
    {/*2*/    /*  triaxial ellipsoids  */  //Z=10327
        params.norm = 1;  //Z=10328
        order = 0;  //Z=10329
        if ( (aell==bell) && (bell==cell) ) area = 4*M_PI*aell*aell;  //Z=10330
                else area = 4*M_PI*pow((pow(aell*bell,8/5.0)+pow(bell*cell,8/5.0)+pow(aell*cell,8/5.0))/3.0,5/8.0);  //Z=10331
        vol = (4*M_PI/3.0)*aell*bell*cell;  //Z=10332
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10333

        /*  homogeneous  */  //Z=10335
        if ( params.cs==0 )
        {/*3*/  //Z=10336

            u1ell[0] = 1;  //Z=10338
            u2ell[0] = 1;  //Z=10339
            v1ell[0] = 1;  //Z=10340
            v2ell[0] = 1;  //Z=10341
            gell[0] = sqrt(M_PI);  //Z=10342

            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            double xrn_n = 1.0;

            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10355
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10344
                fkv[n] = fkv[n-1]*n;  //Z=10345
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10346
                xrn_n = -xrn_n*xr2z;  //Z=10347
                u1ell[n] = (aell*aell-bell*bell)*u1ell[n-1]/n;  //Z=10348
                u2ell[n] = bell*bell*u2ell[n-1]/n;  //Z=10349
                v1ell[n] = (bell*bell-aell*aell)*v1ell[n-1]/n;  //Z=10350
                v2ell[n] = (cell*cell-bell*bell)*v2ell[n-1]/n;  //Z=10351
                gell[n] = gam3[n]/((n+1/2.0)*fkv[n]);  //Z=10352

                double sump = 0.0;  //Z=10368
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10369
                    double sump1 = 0.0;  //Z=10370
                    for ( int k=0; k<=m; k++ )
                    {/*6*/  //Z=10371
                        double sump2 = 0.0;  //Z=10372
                        for ( int ll=0; ll<=n-m; ll++ )
                            sump2 += gell[k+ll]*v1ell[ll]*v2ell[n-m-ll];  //Z=10373
                        sump1 += u1ell[k]*u2ell[m-k]*sump2;  //Z=10374
                    }/*6*/  //Z=10375
                    sump += sump1/((2.0*(n-m))+1);  //Z=10376
                }/*5*/  //Z=10377

                /* fsum[n]:=sumf;  //Z=10379 */
                double sumf = 0.0;  //Z=10380
                for ( int m=0; m<=n; m++ ) sumf += z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10381
                fsum[n] = sumf;  //Z=10382
                /*  P(q)-coefficient  */  //Z=10383
                params.CR->carr4p[n] = (9*M_PI)*pow(4.0,n-1)*z12v[n]*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/((n+3)*(n+2)*(n+3/2.0)*gam3[n]);  //Z=10384
                /*  F(q)-coefficient  */  //Z=10385
                sump = 0.0;  //Z=10386
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10387
                    double sump1 = 0.0;  //Z=10388
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10389
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10390
                }/*5*/  //Z=10391
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=10392
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10393
                    if ( n<n4 ) n4 = n;  //Z=10394
                }/*5*/  //Z=10395
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10396
                    if ( n<n4f ) n4f = n;  //Z=10397
                }/*5*/  //Z=10398
            }/*4*/  //Z=10399
            goto Label99;  //Z=10400
        }/*3*/   /*  of homogeneous  */  //Z=10401

        /*  core/shell  */          /*  not yet ready  */  //Z=10403
        /* if cs=1 then begin  //Z=10404
             goto 99;  //Z=10433
           end;   (* of core/shell *) */  //Z=10434
    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=10435

#ifdef doppelteAbfrage
    /*  perfect orientation case for ellipsoids  */  //Z=10437
    if ( (params.ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10438
        params.norm = 1;  //Z=10439
        order = 1;  //Z=10440
        /*  homogeneous  */  //Z=10441
        if ( params.cs==0 )
        {/*3*/  //Z=10442
            //i = 2;  //Z=10443
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            for ( int n=1; n<=2*nmax+1; n++ )
            {/*4*/  //Z=10444
                //z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10445
                //z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10446
                //fkv[n] = fkv[n-1]*n;  //Z=10447
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10448  gam3[SIZE
                //xrn[n] = -xrn[n-1]*xr2z;  //Z=10450
                //xln[n] = -xln[n-1]*xl2z;  //Z=10451
            }
            GAM3S_INIT(max1,"ordis=6, cs=0: ms+ns")
            GAM3S_INIT(max2,"ordis=6, cs=0: m-ms+n-ns")
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10444
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10445
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10446
                fkv[n] = fkv[n-1]*n;  //Z=10447
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10450
                xln[n] = -xln[n-1]*xl2z;  //Z=10451
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10452
                    double sump = 0.0;  //Z=10453
                    for ( int ns=0; ns<=n; ns++ )
                    {/*6*/  //Z=10454
                        double sump1 = 0.0;  //Z=10455
                        for ( int ms=0; ms<=m; ms++ )
                        {
                            GAM3S_CHECK(max1,ms+ns)
                            GAM3S_CHECK(max2,m-ms+n-ns)
                            sump1 += 1.0/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);  //Z=10456
                        }
                        sump += sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);  //Z=10457
                    }/*6*/  //Z=10458
                    params.CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;  //Z=10459
                }/*5*/  //Z=10460
                /*  P(q)-coefficient  */  //Z=10461
                /* for m:=0 to n do begin  //Z=10462 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10463 */
                /* carr1fm[i]:=carr1pm[i];  //Z=10464 */
                /* i:=i+1;  //Z=10465 */
                /* end;  //Z=10466 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn[n]/(16.0*gam3[n]);  //Z=10467
                /*  F(q)-coefficient  */  //Z=10468
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=10469
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10470
                    if ( n<n4 ) n4 = n;  //Z=10471
                }/*5*/  //Z=10472
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10473
                    if ( n<n4f ) n4f = n;  //Z=10474
                }/*5*/  //Z=10475
            }/*4*/  //Z=10476
            GAM3S_PRINT(max1)
            GAM3S_PRINT(max2)
            goto Label99;  //Z=10477
        }/*3*/   /*  of homogeneous  */  //Z=10478

        /*  core/shell  */          /*  not yet ready  */  //Z=10480
        /* if cs=1 then begin  //Z=10481
             goto 99;  //Z=10510
           end;   (* of core/shell *)   */  //Z=10511
    }/*2*/  /*  of dim=6, triellipsoid  */  //Z=10512
#endif

    /*  isotropic case for super ellipsoids, barrel  */  //Z=10515
    if ( (params.ordis==7) && (dim==7) )
    {/*2*/    /*  super ellipsoids  */  //Z=10516
        params.norm = 1;  //Z=10517
        order = 0;  //Z=10518
        double area, a1;
        qrombchid(params.length,params.radius,params.alphash1,params.sigma,params.alphash1,params.polPhi,params.polTheta,params.polPhi,
                  1/*qx*/,1/*qy*/,1/*qz*/,
                  params.p11,params.p12,params.p13,params.p21,params.p22,params.p23,params.p31,params.p32,params.p33,
                  9/*qx*/,9/*qy*/,9/*0*/,9/*qhkl*/,
                  //params.ax1.length(),params.ax2.length(),params.ax3.length(),
                  //params.ax1.x(),params.ax1.y(),params.ax1.z(),
                  //params.ax2.x(),params.ax2.y(),params.ax2.z(),
                  //params.ax3.x(),params.ax3.y(),params.ax3.z(),
                  //params.sig.x(),params.sig.y(),params.sig.z(),
                  params.ordis,3,8,15,7,0,0,params.CR->carr1p,area);
        area = 2*M_PI*area;  //Z=10520
        vol = 2*M_PI*sqr(params.radius)*params.length*gamma((2+params.alphash1)/params.alphash1)*gamma(1/params.alphash1)/(params.alphash1*gamma((3+params.alphash1)/params.alphash1));  //Z=10521
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10522

        /*  homogeneous  */  //Z=10524
        double u1ell_n = 1.0;             /*  for barrel  */  //Z=10525
        double u2ell_n = 1.0;             /*  for barrel  */  //Z=10526
        v1ell[0] = gamma((2*0+1)/params.alphash1)*u1ell_n;         /*  for barrel  */  //Z=10536
        v2ell[0] = gamma((2*0+2+params.alphash1)/params.alphash1)*u2ell_n;    /*  for barrel  */  //Z=10537
        gell[0]  = gamma((2*0+3+params.alphash1)/params.alphash1);              /*  for barrel  */  //Z=10538

        double z12v[nnmax+1];   // local
        z12v[0] = 1.0;

        //for ( int n=1; n<=3*nmax+1; n++ )
        for ( int n=1; n<=2*nmax+1; n++ )
        {/*3*/  //Z=10527
            //z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10528
            //fkv[n] = fkv[n-1]*n;  //Z=10529
            gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10530  gam3[SIZE
            //v1ell[n] = gamma((2*n+1)/params.alphash1)*u1ell_n;         /*  for barrel  */  //Z=10536
            //v2ell[n] = gamma((2*n+2+params.alphash1)/params.alphash1)*u2ell_n;    /*  for barrel  */  //Z=10537
            //gell[n] = gamma((2*n+3+params.alphash1)/params.alphash1);              /*  for barrel  */  //Z=10538
        }/*3*/  //Z=10533

        if ( params.cs==0 )
        {/*3*/  //Z=10541
            GAM3S_INIT(max1,"ordis=7, cs=0: m-ll+n-k")
            GAM3S_INIT(max2,"ordis=7, cs=0: ll+k")
            double xrn_n = 1.0;
            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=10542
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10528
                fkv[n] = fkv[n-1]*n;  //Z=10529
                u1ell_n = u1ell_n/((n-1/2.0)*n);   /*  for barrel  */  //Z=10531
                u2ell_n = u2ell_n/((n+1)*n);     /*  for barrel  */  //Z=10532
                v1ell[n] = gamma((2*n+1)/params.alphash1)*u1ell_n;         /*  for barrel  */  //Z=10536
                v2ell[n] = gamma((2*n+2+params.alphash1)/params.alphash1)*u2ell_n;    /*  for barrel  */  //Z=10537
                gell[n] = gamma((2*n+3+params.alphash1)/params.alphash1);              /*  for barrel  */  //Z=10538
                xrn_n = -xrn_n*xr2z;
                if ( params.alphash1==2 )
                {/*5*/  //Z=10556
                    a1 = sqr(3/4.0);  //Z=10557
                    for ( int m=0; m<=nmax; m++ )
                    {/*6*/  //Z=10558
                        double sump1 = 0.0;  //Z=10559
                        for ( int k=0; k<=n; k++ )
                        {/*7*/  //Z=10560
                            double sump2 = 0.0;  //Z=10561
                            for ( int ll=0; ll<=m; ll++ )  //Z=10562
                            {
                                GAM3S_CHECK(max1,m-ll+n-k)
                                GAM3S_CHECK(max2,ll+k)
                                sump2 += 1.0/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10563
                            }
                            sump1 += sump2/(fkv[n-k]*fkv[k]);  //Z=10564
                        }/*7*/  //Z=10565
                        params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10566
                    }/*6*/  //Z=10567
                }/*5*/  //Z=10568
                else
                {/*5*/  //Z=10569
                    a1 = gamma((params.alphash1+3)/params.alphash1)/(gamma((params.alphash1+2)/params.alphash1)*gamma(1/params.alphash1));  //Z=10570
                    a1 = a1*a1;  //Z=10571
                    for ( int m=0; m<=nmax; m++ )
                    {/*6*/  //Z=10572
                        double sump1 = 0.0;  //Z=10573
                        for ( int ns=0; ns<=n; ns++ )
                        {/*7*/  //Z=10574
                            double sump2 = 0.0;  //Z=10575
                            for ( int ms=0; ms<=m; ms++ )  //Z=10576
                                sump2 = sump2+v2ell[m-ms]*v2ell[ms]/(gell[ns+ms]*gell[n-ns+m-ms]);  //Z=10577
                            sump1 = sump1+v1ell[n-ns]*v1ell[ns]*sump2;  //Z=10578
                        }/*7*/  //Z=10579
                        params.CR->carr11pm[n][m] = sump1;  //Z=10580
                    }/*6*/  //Z=10581
                }/*5*/  //Z=10582

                /*  orientational average  */  //Z=10584
                double sump = 0.0;  //Z=10585
                for ( int m=0; m<=n; m++ ) sump += gam3[n-m]*z12v[n-m]*z12v[m]*pow(sqr(params.length),n-m)*pow(sqr(params.radius),m)*fkv[m]*params.CR->carr11pm[n-m][m]/(n-m+1/2.0);  //Z=10586
                params.CR->carr4p[n] = a1*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/(2.0*gam3[n]);  //Z=10587


                /* fsum[n]:=sumf;  //Z=10590 */
                /* sumf:=0.0;  //Z=10591 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10592 */
                /* fsum[n]:=sumf;  //Z=10593 */
                /*  P(q)-coefficient  */  //Z=10594
                /*  F(q)-coefficient  */  //Z=10595
                sump = 0.0;  //Z=10596
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10597
                    double sump1 = 0.0;  //Z=10598
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10599
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10600
                }/*5*/  //Z=10601
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=10602
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10603
                    if ( n<n4 ) n4 = n;  //Z=10604
                }/*5*/  //Z=10605
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10606
                    if ( n<n4f ) n4f = n;  //Z=10607
                }/*5*/  //Z=10608
            }/*4*/  //Z=10609
            GAM3S_PRINT(max1)
            GAM3S_PRINT(max2)
            goto Label99;  //Z=10610
        }/*3*/   /*  of homogeneous  */  //Z=10611

        /*  core/shell  */          /*  not yet ready  */  //Z=10613
        /* if cs=1 then begin  //Z=10614
             goto 99;  //Z=10643
           end;   (* of core/shell *) */  //Z=10644
    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=10645

#ifdef doppelteAbfrage
    /*  perfect orientation case for ellipsoids  */  //Z=10647
    if ( (params.ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10648
        params.norm = 1;  //Z=10649
        order = 1;  //Z=10650
        /*  homogeneous  */  //Z=10651
        if ( params.cs==0 )
        {/*3*/  //Z=10652
            //i = 2;  //Z=10653
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10654
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10655
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10656
                fkv[n] = fkv[n-1]*n;  //Z=10657
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10658
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=10659
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10660
                xln[n] = -xln[n-1]*xl2z;  //Z=10661
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10662
                    double sump = 0.0;  //Z=10663
                    for ( int ns=0; ns<=n; ns++ )
                    {/*6*/  //Z=10664
                        double sump1 = 0.0;  //Z=10665
                        for ( int ms=0; ms<=m; ms++ )
                            sump1 += 1.0/(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);  //Z=10666
                        sump += sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);  //Z=10667
                    }/*6*/  //Z=10668
                    params.CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;  //Z=10669
                }/*5*/  //Z=10670
                /*  P(q)-coefficient  */  //Z=10671
                /* for m:=0 to n do begin  //Z=10672 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10673 */
                /* carr1fm[i]:=carr1pm[i];  //Z=10674 */
                /* i:=i+1;  //Z=10675 */
                /* end;  //Z=10676 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn[n]/(16.0*gam3[n]);  //Z=10677
                /*  F(q)-coefficient  */  //Z=10678
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=10679
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10680
                    if ( n<n4 ) n4 = n;  //Z=10681
                }/*5*/  //Z=10682
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10683
                    if ( n<n4f ) n4f = n;  //Z=10684
                }/*5*/  //Z=10685
            }/*4*/  //Z=10686
            goto Label99;  //Z=10687
        }/*3*/   /*  of homogeneous  */  //Z=10688

        /*  core/shell  */          /*  not yet ready  */  //Z=10690
        /* if cs=1 then begin  //Z=10691
             goto 99;  //Z=10720
           end;   (* of core/shell *)   */  //Z=10721
    }/*2*/  /*  of dim=7, super ellipsoid  */  //Z=10722
#endif

    /*  isotropic case for superballs  */  //Z=10726
    if ( (params.ordis==7) && (dim==8) )
    {/*2*/    /*  superball  */  //Z=10727
        //new(carr111pm);  //Z=10728
        static const int ardim = 50;  //Z=11221
        float carr111pm[ardim+1][ardim+1][ardim+1];

        nmax = 40;  //Z=10729
        params.norm = 1;  //Z=10730
        order = 0;  //Z=10731
        /* l:=r;  //Z=10732 */

        /*  radius=a, rm=b, length=c  */  //Z=10734
        qrombdeltac(params.radiusi,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,8,17,0,0,0,params.CR->carr1p,area);  //Z=10735
        area = 8*area;  //Z=10736
        vol = 8*params.radius*params.radiusi*params.length*pow(gamma(1/params.alphash1),3)/(pow(params.alphash1,3)*gamma(1+3/params.alphash1));  //Z=10737
        params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10738


        /*  homogeneous  */  //Z=10741
        double u3ell_n = 1.0;             /*  for superball  */  //Z=10742
        v3ell[0] = exp(gammln((2*0+1)/params.alphash1))*u3ell_n;         /*  for superball  */  //Z=10756
        g3ell[0] = exp(gammln(((2*0+3)/params.alphash1)+1));              /*  for superball  */  //Z=10757
        double z12v[3*nmax+1+1];   // local
        z12v[0] = 1.0;
        for ( n=1; n<=3*nmax+1; n++ )
        {/*3*/  //Z=10743
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10744
            fkv[n] = fkv[n-1]*n;  //Z=10745
            gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10746
            u3ell_n = u3ell_n/((n-1/2.0)*n);   /*  for superball  */  //Z=10747
            v3ell[n] = exp(gammln((2*n+1)/params.alphash1))*u3ell_n;          /*  for superball  */  //Z=10756
            g3ell[n] = exp(gammln(((2*n+3)/params.alphash1)+1));              /*  for superball  */  //Z=10757
        }/*3*/  //Z=10748

        if ( params.cs==0 )
        {/*3*/  //Z=10760
            double xrn_n = 1.0;
            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=10761
                xrn_n = -xrn_n*xr2z;
                if ( params.alphash1==200000 )
                {/*5*/  //Z=10762
                    a1 = sqr(3/4.0);  //Z=10763
                    for ( int m=0; m<=nmax; m++ )
                    {/*6*/  //Z=10764
                        sump1 = 0.0;  //Z=10765
                        for ( int k=0; k<=n; k++ )
                        {/*7*/  //Z=10766
                            sump2 = 0.0;  //Z=10767
                            for ( int ll=0; ll<=m; ll++ )  //Z=10768
                                sump2 = sump2+1/(fkv[m-ll]*fkv[ll]*(m-ll+n-k+3/2.0)*(ll+k+3/2.0)*gam3[m-ll+n-k]*gam3[ll+k]);  //Z=10769
                            sump1 = sump1+sump2/(fkv[n-k]*fkv[k]);  //Z=10770
                        }/*7*/  //Z=10771
                        params.CR->carr11pm[n][m] = M_PI*sump1;  //Z=10772
                    }/*6*/  //Z=10773
                }/*5*/  //Z=10774
                else
                {/*5*/  //Z=10775
                    a1 = gamma(1+3/params.alphash1)/pow(gamma(1/params.alphash1),3);  //Z=10776
                    a1 = a1*a1;  //Z=10777
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=10778
                        for ( int k=0; k<=nmax; k++ )
                        {/*7*/  //Z=10779
                            sump1 = 0.0;  //Z=10780
                            for ( int ns=0; ns<=n; ns++ )
                            {/*8*/  //Z=10781
                                sump2 = 0.0;  //Z=10782
                                for ( int ms=0; ms<=m; ms++ )
                                {/*9*/  //Z=10783
                                    sump3 = 0.0;  //Z=10784
                                    for ( int ks=0; ks<=k; ks++ )  //Z=10785
                                        sump3 = sump3+v3ell[k-ks]*v3ell[ks]/(g3ell[ns+ms+ks]*g3ell[n-ns+m-ms+k-ks]);  //Z=10786
                                    sump2 = sump2+v3ell[m-ms]*v3ell[ms]*sump3;  //Z=10787
                                }/*9*/  //Z=10788
                                sump1 = sump1+v3ell[n-ns]*v3ell[ns]*sump2;  //Z=10789
                            }/*8*/  //Z=10790
                            carr111pm[n][m][k] = sump1;  //Z=10791
                            carr111pm[m][n][k] = sump1;  //Z=10792
                        }/*7*/  //Z=10793
                    }/*6*/  //Z=10794
                }/*5*/  //Z=10795

                /*  orientational average for superball  */  //Z=10797
                double sump = 0.0;  //Z=10798
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10799
                    double sump1 = 0.0;  //Z=10800
                    for ( int k=0; k<=m; k++ )  //Z=10801
                        sump1 += z12v[m-k]*z12v[k]*gam3[m-k]*gam3[k]*pow(params.radiusi*params.radiusi,m-k)*pow(sqr(params.length),k)*carr111pm[n-m][m-k][k]/((m-k+1/2.0)*(k+1/2.0));     /*   //Z=10802 */
                    sump += z12v[n-m]*gam3[n-m]*pow(params.radius*params.radius,n-m)*sump1/(n-m+1/2.0);  //Z=10803
                }/*5*/  //Z=10804
                params.CR->carr4p[n] = (a1/(2.0*M_PI))*pow(-1/(4.0*(z+1)*(z+1)),n)*sump/gam3[n];  //Z=10805


                /* fsum[n]:=sumf;  //Z=10808 */
                /* sumf:=0.0;  //Z=10809 */
                /* for m:=0 to n do sumf:=sumf+z12v[n-m]*z12v[m]/(gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=10810 */
                /* fsum[n]:=sumf;  //Z=10811 */
                /*  P(q)-coefficient  */  //Z=10812
                /*  F(q)-coefficient  */  //Z=10813
                sump = 0.0;  //Z=10814
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10815
                    double sump1 = 0.0;  //Z=10816
                    for ( int k=0; k<=m; k++ ) sump1 += fsum[m-k]*fsum[k]*gam3[m-k]*gam3[k]/((k+1/2.0)*(m-k+1/2.0));  //Z=10817
                    sump += sump1*fsum[n-m]*gam3[n-m]/(n-m+1/2.0);  //Z=10818
                }/*5*/  //Z=10819
                params.CR->carr4f[n] = M_PI*M_PI*xrn_n*sump/(128.0*gam3[n]);  //Z=10820
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10821
                    if ( n<n4 ) n4 = n;  //Z=10822
                }/*5*/  //Z=10823
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10824
                    if ( n<n4f ) n4f = n;  //Z=10825
                }/*5*/  //Z=10826
            }/*4*/  //Z=10827
            goto Label99;  //Z=10828
        }/*3*/   /*  of homogeneous  */  //Z=10829

        /*  core/shell  */          /*  not yet ready  */  //Z=10831
        /* if cs=1 then begin  //Z=10832
             goto 99;  //Z=10861
           end;   (* of core/shell *) */  //Z=10862
        //dispose(carr111pm);  //Z=10863
    }/*2*/  /*  of dim=5, ellipsoids  */  //Z=10864

#ifdef doppelteAbfrage
    /*  perfect orientation case for ellipsoids  */  //Z=10866
    if ( (params.ordis==6) && (dim==5) )
    {/*2*/    /*  ellipsoids  */  //Z=10867
        params.norm = 1;  //Z=10868
        order = 1;  //Z=10869
        /*  homogeneous  */  //Z=10870
        if ( params.cs==0 )
        {/*3*/  //Z=10871
            //i = 2;  //Z=10872
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            GAM3S_INIT(max1,"ordis=6, cs=0: ms+ns")
            GAM3S_INIT(max2,"ordis=6, cs=0: m-ms+n-ns")
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=10873
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10874
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10875
                fkv[n] = fkv[n-1]*n;  //Z=10876
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10877     gam3[SIZE
                //e1[n] = e1[n-1]*(epsi*epsi-1);  //Z=10878
                xrn[n] = -xrn[n-1]*xr2z;  //Z=10879
                xln[n] = -xln[n-1]*xl2z;  //Z=10880
                for ( int m=0; m<=n; m++ )
                {/*5*/  //Z=10881
                    double sump = 0.0;  //Z=10882
                    for ( ns=0; ns<=n; ns++ )
                    {/*6*/  //Z=10883
                        double sump1 = 0.0;  //Z=10884
                        for ( int ms=0; ms<=m; ms++ )
                        {
                            GAM3S_CHECK(max1,ms+ns)
                            GAM3S_CHECK(max2,m-ms+n-ns)
                            sump1 += 1./(fkv[ms]*fkv[m-ms]*(ms+ns+3/2.0)*gam3[ms+ns]*(m-ms+n-ns+3/2.0)*gam3[m-ms+n-ns]);  //Z=10885
                        }
                        sump = sump+sump1*fsum[n-m]/(fkv[ns]*fkv[n-ns]);  //Z=10886
                    }/*6*/  //Z=10887
                    params.CR->carr11pm[n][m] = (9/16.0)*z12vl[n]*z12v[m]*xln[n]*xrn[m]*sump;  //Z=10888
                }/*5*/  //Z=10889
                /*  P(q)-coefficient  */  //Z=10890
                /* for m:=0 to n do begin  //Z=10891 */
                /* carr1pm[i]:=(pi/4)*power(4,n)*z12v[n-m]*z12v[m]*xrn[n]/(gam3[n-m]*gam3[m]*(n-m+1)*fkv[n-m]*(m+1)*fkv[m]);  //Z=10892 */
                /* carr1fm[i]:=carr1pm[i];  //Z=10893 */
                /* i:=i+1;  //Z=10894 */
                /* end;  //Z=10895 */
                params.CR->carr4p[n] = sqrt(M_PI)*pow(4.0,n)*xrn[n]/(16.0*gam3[n]);  //Z=10896
                /*  F(q)-coefficient  */  //Z=10897
                params.CR->carr4f[n] = M_PI*M_PI*xrn[n]/(128.0*gam3[n]);  //Z=10898
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=10899
                    if ( n<n4 ) n4 = n;  //Z=10900
                }/*5*/  //Z=10901
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=10902
                    if ( n<n4f ) n4f = n;  //Z=10903
                }/*5*/  //Z=10904
            }/*4*/  //Z=10905
            GAM3S_PRINT(max1)
            GAM3S_PRINT(max2)
            goto Label99;  //Z=10906
        }/*3*/   /*  of homogeneous  */  //Z=10907

        /*  core/shell  */          /*  not yet ready  */  //Z=10909
        /* if cs=1 then begin  //Z=10910
             goto 99;  //Z=10939
           end;   (* of core/shell *)   */  //Z=10940
    }/*2*/  /*  of dim=8, superball  */  //Z=10941
#endif

    /* ** isotropic case for cylinders and disks ** */  //Z=10945
    if ( (params.ordis==7) && (dim!=3) )
    {/*2*/  //Z=10946
        params.norm = 1;  //Z=10947
        order = 0;  //Z=10948
        /*  homogeneous  */  //Z=10949
        if ( params.cs==0 )
        {/*3*/  //Z=10950
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;

            /*  small axial ratios  */  //Z=10952
            if ( (params.length/params.radius)<2 )
            {/*4*/  //Z=10953
                area = 2*M_PI*sqr(params.radius)+2*M_PI*params.radius*(2.*params.length);  //Z=10954
                vol = M_PI*sqr(params.radius)*(2.*params.length);  //Z=10955
                params.por = 2*M_PI*pow(z+1,4)*area/(z*(z-1)*(z-2)*(z-3)*vol*vol);  //Z=10956

                for ( int n=1; n<=nmax; n++ )
                {/*5*/  //Z=10958
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10959
                    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10960
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=10961
                    fkv[n] = fkv[n-1]*n;  //Z=10962
                    //fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10963
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10964
                    /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10965 */
                    xln[n] = -xln[n-1]*xl2z;  //Z=10966
                    xrn[n] = -xrn[n-1]*xr2z;  //Z=10967
                }/*5*/  //Z=10968

                for ( n=0; n<=nmax; n++ )
                {/*5*/  //Z=10970
                    binsum = 0.0;  //Z=10971
                    /* for m:=0 to n do    (* Cauchy sum *)  //Z=10972 */
                    /*       binsum:=binsum+gam3[m]*z12vl[n-m]*z12v[m]*xln[n-m]*xrn[m]/((n-m+1/2)*(n-m+1)*fkv[n-m]*(m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]);  //Z=10973 */
                    params.CR->carr11pm[n][0] = sqrt(M_PI)/gam3[n];  //Z=10974
                    a1 = sqrt(M_PI)/(2.0*gam3[n]);  //Z=10975
                    for ( int m=1; m<=nmax; m++ )
                    {/*6*/    /*  double sum  */  //Z=10976
                        a1 = a1*(m+1/2.0)/(n+m+1/2.0);  //Z=10977
                        /* carr11pm[n,m]:=power(4,m+1)*gam3[m]*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]*gam3[n+m]);   (* ok *)  //Z=10978 */
                        params.CR->carr11pm[n][m] = pow(4.0,m+1)*a1*z12v[m]*xrn[m]/((m+2)*(m+1)*fkv[m]*(m+1)*fkv[m]);       /*  Mok  */  //Z=10979
                    }/*6*/  //Z=10980
                    params.CR->carr2p[n] = pow(4.0,n-1)*z12vl[n]*xln[n]/((n+1/2.0)*(n+1)*fkv[n]);     /*  double sum  */  //Z=10981
                    params.CR->carr3p[n] = pow(4.0,n)*binsum/gam3[n];      /*  Cauchy sum  */  //Z=10982
                }/*5*/  //Z=10983

            }/*4*/  //Z=10985
            /*  large axial ratios  */  //Z=10986
            else
            {/*4*/  //Z=10987

                double xrn_n = 1.0;
                for ( n=1; n<=nmax; n++ )
                {/*5*/  //Z=10989
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=10990
                    z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=10991
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=10992
                    fkv[n] = fkv[n-1]*n;  //Z=10993
                    fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=10994
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=10995
                    /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=10996 */
                    xln[n] = -xln[n-1]*xl2z;  //Z=10997
                    xrn_n = -xrn_n*xr2z;  //Z=10998
                    /*  cylinder, ok */  //Z=10999
                    if ( dim==1 )
                    {/*6*/  //Z=11000

                        /*  P(q) factorization  */  //Z=11002
                        params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);      /*  P||iso(q)  */  //Z=11003
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);   /*  P-(q)  */  //Z=11004
                        /*  F(q)  */  //Z=11005
                        binsum = 0.0;  //Z=11006
                        for ( int m=0; m<=n; m++ ) binsum += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11007
                        params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11008
                        binsum = 0.0;  //Z=11009
                        for ( int m=0; m<=n; m++ ) binsum += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11010
                        params.CR->carr4f[n] = xrn_n*binsum;  //Z=11011
                    }/*6*/  //Z=11012
                    /*  disk, ok  */  //Z=11013
                    if ( dim==2 )
                    {/*6*/  //Z=11014
                        /*  P(q)  */  //Z=11015
                        params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11016
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11017
                        /*  F(q)  */  //Z=11018
                        binsum = 0.0;  //Z=11019
                        for ( int m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11020
                        params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11021
                        binsum = 0.0;  //Z=11022
                        for ( int m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11023
                        params.CR->carr4f[n] = M_PI*xrn_n*binsum/4.0;  //Z=11024
                    }/*6*/  //Z=11025
                }/*5*/  /*  of large axial ratios  */  //Z=11026

                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11028
                    if ( n<n1 ) n1 = n;  //Z=11029
                }/*5*/  //Z=11030
                if ( fabs(params.CR->carr2p[n])<min )
                {/*5*/  //Z=11031
                    if ( n<n2 ) n2 = n;  //Z=11032
                }/*5*/  //Z=11033
                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=11034
                    if ( n<n3 ) n3 = n;  //Z=11035
                }/*5*/  //Z=11036
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11037
                    if ( n<n1f ) n1f = n;  //Z=11038
                }/*5*/  //Z=11039
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11040
                    if ( n<n4 ) n4 = n;  //Z=11041
                }/*5*/  //Z=11042
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11043
                    if ( n<n4f ) n4f = n;  //Z=11044
                }/*5*/  //Z=11045
            }/*4*/  //Z=11046
        }/*3*/ /*  of homogeneous  */  //Z=11047

        /*  core/shell  */  //Z=11049
        if ( params.cs==1 )
        {/*3*/  //Z=11050
            double xrmn_n = 1.0;
            double xrn_n  = 1.0;
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=11051
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11052
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11053
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11054
                fkv[n] = fkv[n-1]*n;  //Z=11055
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11056
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11057
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11058 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11059
                xrn_n = -xrn_n*xr2z;  //Z=11060
                xrmn_n = -xrmn_n*xrm2z;  //Z=11061
                pn[n] = pn[n-1]*p*p;  //Z=11062
                /* ** cylinder ** */  //Z=11063
                if ( dim==1 )
                {/*5*/  //Z=11064
                    /*  longitudinal P(q)  */  //Z=11065
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11066
                    /*  cross-sectional P(q)  */  //Z=11067
                    /*  F121  */  //Z=11068
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn_n/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11069
                    /*  F122  */  //Z=11070
                    double sump = 0.0;  //Z=11071
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11072
                        sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11073
                    }/*6*/  //Z=11074
                    params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=11075
                    /*  F123  */  //Z=11076
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11077

                    /*  longitudinal F(q)  */  //Z=11079
                    binsum = 0.0;  //Z=11080
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11081
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11082
                    /*  cross-sectional F(q)  */  //Z=11083
                    /*  F121  */  //Z=11084
                    sump = 0.0;  //Z=11085
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11086
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11087
                    }/*6*/  //Z=11088
                    params.CR->carr4f[n] = xrn_n*sump;  //Z=11089
                    /*  F122  */  //Z=11090
                    sump = 0.0;  //Z=11091
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11092
                        sump = sump+z12v[m]*z12v[n-m]*pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11093
                    }/*6*/  //Z=11094
                    params.CR->carr5f[n] = xrmn_n*sump;  //Z=11095
                    /*  F123  */  //Z=11096
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=11097
                }/*5*/  //Z=11098

                /* ** disk ** */  //Z=11100
                if ( dim==2 )
                {/*5*/  //Z=11101
                    /*  longitudinal  */  //Z=11102
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11103
                    /*  cross-sectional P(q)  */  //Z=11104
                    /*  F121  */  //Z=11105
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn_n/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11106
                    /*  F122  */  //Z=11107
                    double sump = 0.0;  //Z=11108
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11109
                        sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11110
                    }/*6*/  //Z=11111
                    params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=11112
                    /*  F123  */  //Z=11113
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11114

                    /*  longitudinal F(q)  */  //Z=11116
                    binsum = 0.0;  //Z=11117
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11118
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11119
                    /*  cross-sectional F(q)  */  //Z=11120
                    /*  F121  */  //Z=11121
                    sump = 0.0;  //Z=11122
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11123
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11124
                    }/*6*/  //Z=11125
                    params.CR->carr4f[n] = M_PI*xrn_n*sump/4.0;  //Z=11126
                    /*  F122  */  //Z=11127
                    sump = 0.0;  //Z=11128
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11129
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11130
                    }/*6*/  //Z=11131
                    params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=11132
                    /*  F123  */  //Z=11133
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=11134
                }/*5*/  //Z=11135
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11136
                    if ( n<n1 ) n1 = n;  //Z=11137
                }/*5*/  //Z=11138
                if ( fabs(params.CR->carr3p[n])<min )
                {/*5*/  //Z=11139
                    if ( n<n3 ) n3 = n;  //Z=11140
                }/*5*/  //Z=11141
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11142
                    if ( n<n4 ) n4 = n;  //Z=11143
                }/*5*/  //Z=11144
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=11145
                    if ( n<n5 ) n5 = n;  //Z=11146
                }/*5*/  //Z=11147
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=11148
                    if ( n<n6 ) n6 = n;  //Z=11149
                }/*5*/  //Z=11150
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11151
                    if ( n<n1f ) n1f = n;  //Z=11152
                }/*5*/  //Z=11153
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11154
                    if ( n<n4f ) n4f = n;  //Z=11155
                }/*5*/  //Z=11156
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=11157
                    if ( n<n5f ) n5f = n;  //Z=11158
                }/*5*/  //Z=11159
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=11160
                    if ( n<n6f ) n6f = n;  //Z=11161
                }/*5*/  //Z=11162
            }/*4*/  //Z=11163
        }/*3*/ /*  of core/shell  */  //Z=11164

        /*  inhomogeneous core/shell  */  //Z=11166
        if ( params.cs==2 )
        {/*3*/  //Z=11167
            double xrmn_n = 1.0;
            double z12v_n = 1.0;
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11168
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11169
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11170
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11171
                fkv[n] = fkv[n-1]*n;  //Z=11172
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11173
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11174
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11175 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11176
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11177
                xrmn_n = -xrmn_n*xrm2z;  //Z=11178
                pn[n] = pn[n-1]*p*p;  //Z=11179
                /* ** cylinder ** */  //Z=11180
                if ( dim==1 )
                {/*5*/  //Z=11181
                    /*  longitudinal P(q)  */  //Z=11182
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11183

                    /*  cross-sectional P(q)  */  //Z=11185
                    params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v_n*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11186
                    double sump = 0.0;  //Z=11187
                    double sump1 = 0.0;  //Z=11188
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11189
                        const double sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=11190
                        sump += pn[n-m]*sumi;  //Z=11191
                        sump1 += sumi;  //Z=11192
                    }/*6*/  //Z=11193
                    params.CR->carr5p[n] = (1-params.uca/2.0)*z12v_n*xrmn_n*sump;  //Z=11194
                    params.CR->carr6p[n] = (1-params.uca/2.0)*z12v_n*xrn[n]*sump1;  //Z=11195
                    sump = 0.0;  //Z=11196
                    sump1 = 0.0;  //Z=11197
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11198
                        const double sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=11199
                        sump += sumi;  //Z=11200
                        sump1 += pn[n-m]*sumi;  //Z=11201
                    }/*6*/  //Z=11202
                    params.CR->carr7p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump;  //Z=11203
                    params.CR->carr8p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump1;  //Z=11204
                    params.CR->carr9p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v_n*xrn[n]*sump;  //Z=11205


                    /* (* cross-sectional P(q) *)  //Z=11208
                       (* F121 *)  //Z=11209
                       carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11210
                       (* F122 *)  //Z=11211
                       sump:=0.0;  //Z=11212
                          for m:=0 to n do begin  //Z=11213
                          sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11214
                       end;  //Z=11215
                       carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=11216
                       (* F123 *)  //Z=11217
                       carr6p[n]:=carr4p[n]/pn[n];  */  //Z=11218

                    /*  longitudinal F(q)  */  //Z=11220
                    binsum = 0.0;  //Z=11221
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11222
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11223
                    /*  cross-sectional F(q)  */  //Z=11224
                    params.CR->carr4f[n] = z12v_n*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=11225
                    params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v_n*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11226
                    params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v_n*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11227

                }/*5*/  //Z=11229

                /* ** disk ** */  //Z=11231
                if ( dim==2 )
                {/*5*/  //Z=11232
                    /*  longitudinal  */  //Z=11233
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11234


                    /*  cross-sectional P(q)  */  //Z=11237
                    params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v_n*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=11238
                    double sump = 0.0;  //Z=11239
                    double sump1 = 0.0;  //Z=11240
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11241
                        const double sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11242
                        sump += pn[n-m]*sumi;  //Z=11243
                        sump1 += sumi;  //Z=11244
                    }/*6*/  //Z=11245
                    params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v_n*xrmn_n*sump;  //Z=11246
                    params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v_n*xrn[n]*sump1;  //Z=11247
                    sump = 0.0;  //Z=11248
                    sump1 = 0.0;  //Z=11249
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11250
                        const double sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=11251
                        sump += sumi;  //Z=11252
                        sump1 += pn[n-m]*sumi;  //Z=11253
                    }/*6*/  //Z=11254
                    params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrmn_n*sump;  //Z=11255
                    params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrmn_n*sump1;  //Z=11256
                    params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrn[n]*sump;  //Z=11257


                    /*  cross-sectional P(q)  */  //Z=11260
                    /* (* F121 *)  //Z=11261
                       carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=11262
                       (* F122 *)  //Z=11263
                       sump:=0.0;  //Z=11264
                          for m:=0 to n do begin  //Z=11265
                          sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11266
                       end;  //Z=11267
                       carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=11268
                       (* F123 *)  //Z=11269
                       carr6p[n]:=carr4p[n]/pn[n];   */  //Z=11270

                    /*  longitudinal F(q)  */  //Z=11272
                    binsum = 0.0;  //Z=11273
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11274
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11275
                    /*  cross-sectional F(q)  */  //Z=11276
                    params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v_n*xrn[n]/(gam3[n]*fkv[n]);  //Z=11277
                    params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v_n*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11278
                    params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v_n*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11279
                }/*5*/  //Z=11280
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11281
                    if ( n<n1 ) n1 = n;  //Z=11282
                }/*5*/  //Z=11283
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11284
                    if ( n<n4 ) n4 = n;  //Z=11285
                }/*5*/  //Z=11286
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=11287
                    if ( n<n5 ) n5 = n;  //Z=11288
                }/*5*/  //Z=11289
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=11290
                    if ( n<n6 ) n6 = n;  //Z=11291
                }/*5*/  //Z=11292
                if ( fabs(params.CR->carr7p[n])<min )
                {/*5*/  //Z=11293
                    if ( n<n7 ) n7 = n;  //Z=11294
                }/*5*/  //Z=11295
                if ( fabs(params.CR->carr8p[n])<min )
                {/*5*/  //Z=11296
                    if ( n<n8 ) n8 = n;  //Z=11297
                }/*5*/  //Z=11298
                if ( fabs(params.CR->carr9p[n])<min )
                {/*5*/  //Z=11299
                    if ( n<n9 ) n9 = n;  //Z=11300
                }/*5*/  //Z=11301
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11302
                    if ( n<n1f ) n1f = n;  //Z=11303
                }/*5*/  //Z=11304
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11305
                    if ( n<n4f ) n4f = n;  //Z=11306
                }/*5*/  //Z=11307
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=11308
                    if ( n<n5f ) n5f = n;  //Z=11309
                }/*5*/  //Z=11310
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=11311
                    if ( n<n6f ) n6f = n;  //Z=11312
                }/*5*/  //Z=11313
                if ( fabs(params.CR->carr7f[n])<min )
                {/*5*/  //Z=11314
                    if ( n<n7f ) n7f = n;  //Z=11315
                }/*5*/  //Z=11316
                if ( fabs(params.CR->carr8f[n])<min )
                {/*5*/  //Z=11317
                    if ( n<n8f ) n8f = n;  //Z=11318
                }/*5*/  //Z=11319
                if ( fabs(params.CR->carr9f[n])<min )
                {/*5*/  //Z=11320
                    if ( n<n9f ) n9f = n;  //Z=11321
                }/*5*/  //Z=11322
            }/*4*/  //Z=11323
        }/*3*/ /*  of inhomogeneous core/shell  */  //Z=11324


        /*  myelin  */  //Z=11327
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=11328
            //i = 2;  //Z=11329
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11330
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11331
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11332
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11333
                fkv[n] = fkv[n-1]*n;  //Z=11334
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11335
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11336
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11337 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11338
                /* xrn[n]:=-xrn[n-1]*xr2z;  //Z=11339 */
                xrn[n] = -xrn[n-1]*x12zm;         /*  myelin radius  */  //Z=11340
                /*  cylinder, ok */  //Z=11341
                if ( dim==1 )
                {/*5*/  //Z=11342
                    /*  P(q)  */  //Z=11343
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(2*n+1)*(n+1)*gam3[n]*fkv[n]);  //Z=11344

                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11346
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11347 */
                        /* i:=i+1;  //Z=11348 */
                        params.CR->carr11pm[n][m] = 1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11349
                    }/*6*/  //Z=11350
                    params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=11351
                    /* carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11352 */


                    /*  F(q)  */  //Z=11355
                    binsum = 0.0;  //Z=11356
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11357
                    params.CR->carr1f[n] = M_PI*xln[n]*binsum/(4.0*(2*n+1));  //Z=11358
                    binsum = 0.0;  //Z=11359
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11360
                    params.CR->carr4f[n] = xrn[n]*binsum;  //Z=11361
                }/*5*/  //Z=11362
                /*  disk, ok  */  //Z=11363
                if ( dim==2 )
                {/*5*/  //Z=11364
                    /*  P(q)  */  //Z=11365
                    params.CR->carr1p[n] = 2*pow(4.0,n)*z12vl[n]*xln[n]/((n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=11366
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11367
                        /* carr1pm[i]:=1/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11368 */
                        /* i:=i+1;  //Z=11369 */
                        params.CR->carr11pm[n][m] = (M_PI/4.0)*(1/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]));  //Z=11370
                    }/*6*/  //Z=11371
                    params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=11372

                    /* carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=11374 */
                    /*  F(q)  */  //Z=11375
                    binsum = 0.0;  //Z=11376
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11377
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*binsum/(2.0*gam3[n]);  //Z=11378
                    binsum = 0.0;  //Z=11379
                    for ( int m=0; m<=n; m++ ) binsum = binsum+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11380
                    params.CR->carr4f[n] = M_PI*xrn[n]*binsum/4.0;  //Z=11381
                }/*5*/  //Z=11382
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11383
                    if ( n<n1 ) n1 = n;  //Z=11384
                }/*5*/  //Z=11385
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11386
                    if ( n<n1f ) n1f = n;  //Z=11387
                }/*5*/  //Z=11388
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11389
                    if ( n<n4 ) n4 = n;  //Z=11390
                }/*5*/  //Z=11391
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11392
                    if ( n<n4f ) n4f = n;  //Z=11393
                }/*5*/  //Z=11394
            }/*4*/  //Z=11395
        }/*3*/ /*  of myelin  */  //Z=11396

    }/*2*/  //Z=11398

    /* ** perfect orientation for cylinders and disks ** */  //Z=11400
    if ( (params.ordis==6) && (dim!=3) )
    {/*2*/  //Z=11401
        params.norm = 1;  //Z=11402
        order = 1;  //Z=11403
        double z12v[nnmax+1];   // local
        z12v[0] = 1.0;
        /*  homogeneous  */  //Z=11404
        if ( params.cs==0 )
        {/*3*/  //Z=11405
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11406
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11407
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11408
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11409
                fkv[n] = fkv[n-1]*n;  //Z=11410
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11411
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11412
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11413 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11414
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11415
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11416
                    /*  P(q)-coefficients  */  //Z=11417
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);       /*  P||(q)  */  //Z=11418
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  /*  P-(q)  */  //Z=11419
                    /*  F(q)-coefficients  */  //Z=11420
                    double sump = 0.0;  //Z=11421
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11422
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11423
                    }/*6*/  //Z=11424
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=11425
                    sump = 0.0;  //Z=11426
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11427
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11428
                    }/*6*/  //Z=11429
                    params.CR->carr4f[n] = xrn[n]*sump;  //Z=11430
                }/*5*/  //Z=11431
                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=11432
                    /*  P(q)-coefficients  */  //Z=11433
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);    /*  P-(q)  */  //Z=11434
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);           /*  P||(q)  */  //Z=11435
                    /*  F(q)-coefficients  */  //Z=11436
                    double sump = 0.0;  //Z=11437
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11438
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11439
                    }/*6*/  //Z=11440
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=11441
                    sump = 0.0;  //Z=11442
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11443
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11444
                    }/*6*/  //Z=11445
                    params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=11446
                }/*5*/  //Z=11447
                if ( search1 )
                {/*5*/  //Z=11448
                    if ( params.CR->carr1p[n]<1e-50 )
                    {/*6*/  //Z=11449
                        n1 = n;  //Z=11450
                        search1 = false;  //Z=11451
                    }/*6*/  //Z=11452
                }/*5*/  //Z=11453
                if ( search4 )
                {/*5*/  //Z=11454
                    if ( params.CR->carr4p[n]<1e-50 )
                    {/*6*/  //Z=11455
                        n4 = n;  //Z=11456
                        search4 = false;  //Z=11457
                    }/*6*/  //Z=11458
                }/*5*/  //Z=11459
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11460
                    if ( n<n1f ) n1f = n;  //Z=11461
                }/*5*/  //Z=11462
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11463
                    if ( n<n4f ) n4f = n;  //Z=11464
                }/*5*/  //Z=11465
            }/*4*/  //Z=11466
        }/*3*/  //Z=11467

        /*  core/shell  */  //Z=11469
        if ( params.cs==1 )
        {/*3*/  //Z=11470
            double xrmn_n = 1.0;
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=11471
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11472
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11473
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11474
                fkv[n] = fkv[n-1]*n;  //Z=11475
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11476
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11477
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11478 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11479
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11480
                xrmn_n = -xrmn_n*xrm2z;  //Z=11481
                pn[n] = pn[n-1]*p*p;  //Z=11482
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11483
                    /*  P(q)-coefficients  */  //Z=11484
                    /*  longitudinal  */  //Z=11485
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11486
                    /*  cross-sectional  */  //Z=11487
                    /*  F121  */  //Z=11488
                    params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11489
                    /*  F122  */  //Z=11490
                    double sump = 0.0;  //Z=11491
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11492
                        sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11493
                    }/*6*/  //Z=11494
                    params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=11495
                    /*  F123  */  //Z=11496
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11497

                    /*  F(q)-coefficients  */  //Z=11499
                    /*  longitudinal  */  //Z=11500
                    sump = 0.0;  //Z=11501
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11502
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11503
                    }/*6*/  //Z=11504
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=11505
                    /*  cross-sectional  */  //Z=11506
                    /*  F121  */  //Z=11507
                    sump = 0.0;  //Z=11508
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11509
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11510
                    }/*6*/  //Z=11511
                    params.CR->carr4f[n] = xrn[n]*sump;  //Z=11512
                    /*  F122  */  //Z=11513
                    sump = 0.0;  //Z=11514
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11515
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11516
                    }/*6*/  //Z=11517
                    params.CR->carr5f[n] = xrmn_n*sump;  //Z=11518
                    /*  F123  */  //Z=11519
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=11520
                }/*5*/  //Z=11521

                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=11523
                    /*  P(q)-coefficients  */  //Z=11524
                    /*  longitudinal  */  //Z=11525
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=11526
                    /*  cross-sectional  */  //Z=11527
                    /*  F121  */  //Z=11528
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11529
                    /*  F122  */  //Z=11530
                    /* sump:=0.0;  //Z=11531 */
                    /*    for m:=0 to n do begin  //Z=11532 */
                    /*    sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11533 */
                    /* end;  //Z=11534 */
                    /* carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=11535 */

                    /*  F122  */  //Z=11537
                    double sump = 0.0;  //Z=11538
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11539
                        sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11540
                    }/*6*/  //Z=11541
                    params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=11542


                    /*  F123  */  //Z=11545
                    params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11546
                    /*  F(q)-coefficients  */  //Z=11547
                    /*  longitudinal  */  //Z=11548
                    sump = 0.0;  //Z=11549
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11550
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11551
                    }/*6*/  //Z=11552
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=11553
                    /*  cross-sectional  */  //Z=11554
                    /*  F121  */  //Z=11555
                    sump = 0.0;  //Z=11556
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11557
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11558
                    }/*6*/  //Z=11559
                    params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=11560
                    /*  F122  */  //Z=11561
                    sump = 0.0;  //Z=11562
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11563
                        sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11564
                    }/*6*/  //Z=11565
                    params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=11566
                    /*  F123  */  //Z=11567
                    params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=11568
                }/*5*/  //Z=11569
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11570
                    if ( n<n1 ) n1 = n;  //Z=11571
                }/*5*/  //Z=11572
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11573
                    if ( n<n4 ) n4 = n;  //Z=11574
                }/*5*/  //Z=11575
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=11576
                    if ( n<n5 ) n5 = n;  //Z=11577
                }/*5*/  //Z=11578
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=11579
                    if ( n<n6 ) n6 = n;  //Z=11580
                }/*5*/  //Z=11581
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11582
                    if ( n<n1f ) n1f = n;  //Z=11583
                }/*5*/  //Z=11584
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11585
                    if ( n<n4f ) n4f = n;  //Z=11586
                }/*5*/  //Z=11587
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=11588
                    if ( n<n5f ) n5f = n;  //Z=11589
                }/*5*/  //Z=11590
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=11591
                    if ( n<n6f ) n6f = n;  //Z=11592
                }/*5*/  //Z=11593
            }/*4*/  //Z=11594
        }/*3*/  /*  of core/shell  */  //Z=11595

        /*  inhomogeneous core/shell  */  //Z=11597
        if ( params.cs==2 )
        {/*3*/  //Z=11598
            double xrmn_n = 1.0;
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=11599
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11600
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11601
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11602
                fkv[n] = fkv[n-1]*n;  //Z=11603
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11604
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11605
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11606 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11607
                xrn[n] = -xrn[n-1]*xr2z;  //Z=11608
                xrmn_n = -xrmn_n*xrm2z;  //Z=11609
                pn[n] = pn[n-1]*p*p;  //Z=11610
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11611
                    /*  P(q)-coefficients  */  //Z=11612
                    /*  longitudinal  */  //Z=11613
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11614

                    /*  cross-sectional P(q)  */  //Z=11616
                    params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11617
                    double sump = 0.0;  //Z=11618
                    double sump1 = 0.0;  //Z=11619
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11620
                        const double sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=11621
                        sump += pn[n-m]*sumi;  //Z=11622
                        sump1 += sumi;  //Z=11623
                    }/*6*/  //Z=11624
                    params.CR->carr5p[n] = (1-params.uca/2.0)*z12v[n]*xrmn_n*sump;  //Z=11625
                    params.CR->carr6p[n] = (1-params.uca/2.0)*z12v[n]*xrn[n]*sump1;  //Z=11626
                    sump = 0.0;  //Z=11627
                    sump1 = 0.0;  //Z=11628
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11629
                        const double sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=11630
                        sump += sumi;  //Z=11631
                        sump1 += pn[n-m]*sumi;  //Z=11632
                    }/*6*/  //Z=11633
                    params.CR->carr7p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=11634
                    params.CR->carr8p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=11635
                    params.CR->carr9p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=11636

                    /*  cross-sectional  */  //Z=11638
                    /*  (* F121 *)  //Z=11639
                       carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11640
                       (* F122 *)  //Z=11641
                       sump:=0.0;  //Z=11642
                       for m:=0 to n do begin  //Z=11643
                          sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11644
                       end;  //Z=11645
                       carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=11646
                       (* F123 *)  //Z=11647
                       carr6p[n]:=carr4p[n]/pn[n];    */  //Z=11648

                    /*  F(q)-coefficients  */  //Z=11650
                    /*  longitudinal  */  //Z=11651
                    sump = 0.0;  //Z=11652
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11653
                        sump += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11654
                    }/*6*/  //Z=11655
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=11656
                    /*  cross-sectional  */  //Z=11657
                    params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=11658
                    params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11659
                    params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=11660
                }/*5*/  //Z=11661

                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=11663
                    /*  P(q)-coefficients  */  //Z=11664
                    /*  longitudinal  */  //Z=11665
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=11666

                    /*  cross-sectional P(q)  */  //Z=11668
                    params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=11669
                    double sump = 0.0;  //Z=11670
                    double sump1 = 0.0;  //Z=11671
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11672
                        const double sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=11673
                        sump += pn[n-m]*sumi;  //Z=11674
                        sump1 += sumi;  //Z=11675
                    }/*6*/  //Z=11676
                    params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=11677
                    params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=11678
                    sump = 0.0;  //Z=11679
                    sump1 = 0.0;  //Z=11680
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11681
                        const double sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=11682
                        sump += sumi;  //Z=11683
                        sump1 += pn[n-m]*sumi;  //Z=11684
                    }/*6*/  //Z=11685
                    params.CR->carr7p[n] = (M_PI/4.0)*(1-params.alphash1)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=11686
                    params.CR->carr8p[n] = (M_PI/4.0)*(1-params.alphash1)*(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=11687
                    params.CR->carr9p[n] = (M_PI/4.0)*(1-params.alphash1)*(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=11688


                    /*  cross-sectional P(q)  */  //Z=11691
                    /*  F121  */  //Z=11692
                    /*  carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=11693
                        (* F122 *)  //Z=11694
                        //sump:=0.0;  //Z=11695
                        //   for m:=0 to n do begin  //Z=11696
                        //   sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11697
                        //end;  //Z=11698
                        //carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=11699

                        (* F122 *)  //Z=11701
                        sump:=0.0;  //Z=11702
                        for m:=0 to n do begin  //Z=11703
                           sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11704
                        end;  //Z=11705
                        carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=11706
                        (* F123 *)  //Z=11707
                        carr6p[n]:=carr4p[n]/pn[n];        */  //Z=11708

                    /*  F(q)-coefficients  */  //Z=11710
                    /*  longitudinal  */  //Z=11711
                    sump = 0.0;  //Z=11712
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11713
                        sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11714
                    }/*6*/  //Z=11715
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=11716
                    /*  cross-sectional  */  //Z=11717
                    params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=11718
                    params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11719
                    params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=11720

                }/*5*/  //Z=11722
                if ( fabs(params.CR->carr1p[n])<min )
                {/*5*/  //Z=11723
                    if ( n<n1 ) n1 = n;  //Z=11724
                }/*5*/  //Z=11725
                if ( fabs(params.CR->carr4p[n])<min )
                {/*5*/  //Z=11726
                    if ( n<n4 ) n4 = n;  //Z=11727
                }/*5*/  //Z=11728
                if ( fabs(params.CR->carr5p[n])<min )
                {/*5*/  //Z=11729
                    if ( n<n5 ) n5 = n;  //Z=11730
                }/*5*/  //Z=11731
                if ( fabs(params.CR->carr6p[n])<min )
                {/*5*/  //Z=11732
                    if ( n<n6 ) n6 = n;  //Z=11733
                }/*5*/  //Z=11734
                if ( fabs(params.CR->carr7p[n])<min )
                {/*5*/  //Z=11735
                    if ( n<n7 ) n7 = n;  //Z=11736
                }/*5*/  //Z=11737
                if ( fabs(params.CR->carr8p[n])<min )
                {/*5*/  //Z=11738
                    if ( n<n8 ) n8 = n;  //Z=11739
                }/*5*/  //Z=11740
                if ( fabs(params.CR->carr9p[n])<min )
                {/*5*/  //Z=11741
                    if ( n<n9 ) n9 = n;  //Z=11742
                }/*5*/  //Z=11743
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11744
                    if ( n<n1f ) n1f = n;  //Z=11745
                }/*5*/  //Z=11746
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11747
                    if ( n<n4f ) n4f = n;  //Z=11748
                }/*5*/  //Z=11749
                if ( fabs(params.CR->carr5f[n])<min )
                {/*5*/  //Z=11750
                    if ( n<n5f ) n5f = n;  //Z=11751
                }/*5*/  //Z=11752
                if ( fabs(params.CR->carr6f[n])<min )
                {/*5*/  //Z=11753
                    if ( n<n6f ) n6f = n;  //Z=11754
                }/*5*/  //Z=11755
                if ( fabs(params.CR->carr7f[n])<min )
                {/*5*/  //Z=11756
                    if ( n<n7f ) n7f = n;  //Z=11757
                }/*5*/  //Z=11758
                if ( fabs(params.CR->carr8f[n])<min )
                {/*5*/  //Z=11759
                    if ( n<n8f ) n8f = n;  //Z=11760
                }/*5*/  //Z=11761
                if ( fabs(params.CR->carr9f[n])<min )
                {/*5*/  //Z=11762
                    if ( n<n9f ) n9f = n;  //Z=11763
                }/*5*/  //Z=11764
            }/*4*/  //Z=11765
        }/*3*/  /*  of inhomogeneous core/shell  */  //Z=11766


        /*  myelin  */  //Z=11769
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=11770
            for ( int n=1; n<=nmax; n++ )
            {/*4*/  //Z=11771
                z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11772
                z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11773
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11774
                fkv[n] = fkv[n-1]*n;  //Z=11775
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11776
                gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11777
                /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11778 */
                xln[n] = -xln[n-1]*xl2z;  //Z=11779
                xrn[n] = -xrn[n-1]*x12zm;  //Z=11780
                if ( dim==1 )
                {/*5*/ /*  cylinder  */  //Z=11781
                    /*  P(q)-coefficients  */  //Z=11782
                    params.CR->carr1p[n] = sqrt(M_PI)*pow(4.0,n)*z12vl[n]*xln[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11783

                    params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=11785

                    /*  F(q)-coefficients  */  //Z=11787
                    double sump = 0.0;  //Z=11788
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11789
                        sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11790
                    }/*6*/  //Z=11791
                    params.CR->carr1f[n] = M_PI*xln[n]*sump/4.0;  //Z=11792
                    sump = 0.0;  //Z=11793
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11794
                        sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11795
                    }/*6*/  //Z=11796
                    params.CR->carr4f[n] = xrn[n]*sump;  //Z=11797
                }/*5*/  //Z=11798
                if ( dim==2 )
                {/*5*/ /*  disk  */  //Z=11799
                    /*  P(q)-coefficients  */  //Z=11800
                    params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=11801
                    params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=11802
                    /*  F(q)-coefficients  */  //Z=11803
                    double sump = 0.0;  //Z=11804
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11805
                        sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11806
                    }/*6*/  //Z=11807
                    params.CR->carr1f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=11808
                    sump = 0.0;  //Z=11809
                    for ( int m=0; m<=n; m++ )
                    {/*6*/  //Z=11810
                        sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11811
                    }/*6*/  //Z=11812
                    params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=11813
                }/*5*/  //Z=11814
                if ( search1 )
                {/*5*/  //Z=11815
                    if ( params.CR->carr1p[n]<1e-50 )
                    {/*6*/  //Z=11816
                        n1 = n;  //Z=11817
                        search1 = false;  //Z=11818
                    }/*6*/  //Z=11819
                }/*5*/  //Z=11820
                if ( search4 )
                {/*5*/  //Z=11821
                    if ( params.CR->carr4p[n]<1e-50 )
                    {/*6*/  //Z=11822
                        n4 = n;  //Z=11823
                        search4 = false;  //Z=11824
                    }/*6*/  //Z=11825
                }/*5*/  //Z=11826
                if ( fabs(params.CR->carr1f[n])<min )
                {/*5*/  //Z=11827
                    if ( n<n1f ) n1f = n;  //Z=11828
                }/*5*/  //Z=11829
                if ( fabs(params.CR->carr4f[n])<min )
                {/*5*/  //Z=11830
                    if ( n<n4f ) n4f = n;  //Z=11831
                }/*5*/  //Z=11832
            }/*4*/  //Z=11833
        }/*3*/  /*  of myelin  */  //Z=11834



    }/*2*/  //Z=11838

    /* ** orientational distribution for cylinders and disks ** */  //Z=11840
    if ( (params.ordis==0) && (dim!=3) )
    {/*2*/  //Z=11841
        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,2,0,0,0,0,params.CR->carr1p,params.norm);  //Z=11842
        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,3,0,0,0,0,params.CR->carr1p,order);  //Z=11843
        order = order/params.norm;  //Z=11844

        /*  use phi=0 and rotate qx/qy-axis  */  //Z=11846
        if ( params.orcase==1 )
        {/*3*/  //Z=11847
            params.polPhi = 0;  //Z=11848
            params.p11 = -cos(params.polPhi*M_PI/180.0)*cos(params.polTheta*M_PI/180.0);       /*  = -cos(theta*pi/180);  //Z=11849 */
            params.p12 = sin(params.polPhi*M_PI/180.0);                                        /*  = 0;  //Z=11850 */
            params.p13 = cos(params.polPhi*M_PI/180.0)*sin(params.polTheta*M_PI/180.0);        /*  =  sin(theta*pi/180);  //Z=11851 */
            params.p21 = -cos(params.polPhi*M_PI/180.0);                                       /*  = -1;  //Z=11852 */
            params.p22 = -sin(params.polPhi*M_PI/180.0)*cos(params.polTheta*M_PI/180.0);       /*  = 0;  //Z=11853 */
            params.p23 = sin(params.polPhi*M_PI/180.0)*sin(params.polTheta*M_PI/180.0);        /*  = 0;  //Z=11854 */
            params.p31 = -sin(params.polTheta*M_PI/180.0);                                     /*  = 0;  //Z=11855 */
            params.p32 = 0;  //Z=11856
            params.p33 = -cos(params.polTheta*M_PI/180.0);                                     /*  = -cos(theta*pi/180);  //Z=11857 */

            for ( int n=0; n<=nmax; n++ )
            {/*4*/  //Z=11859
                for ( int m=0; m<=nmax; m++ )
                {/*5*/  //Z=11860
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,0,0,2*n,2*m,params.CR->carr1p,intl);  //Z=11861
                    intlar[n][m] = intl;  //Z=11862
                }/*5*/  //Z=11863
            }/*4*/  //Z=11864

            /* for n:=1 to nmax do begin  //Z=11866
                fkv[n]:=fkv[n-1]*n;  //Z=11867
                fk2v[n]:=fk2v[n-1]*(2*n-1)*(2*n);  //Z=11868
                if odd(n) then gam1[n]:=fkv[round((n-1)/2)]  //Z=11869
                   else gam1[n]:=fkv[n]*sqrt(pi)/(fkv[round(n/2)]*power(2,n));  //Z=11870
                gam2[n]:=(n/2)*fkv[n-1]*sqrt(pi)/(gam1[n]*power(2,n-1));  //Z=11871
                gam3[n]:=gam3[n-1]*(2*n+1)/2;  //Z=11872
             end;  */  //Z=11873

            /*  cylinder form factor coefficients  */  //Z=11875
            if ( dim==1 )
            {/*4*/  //Z=11876
                /*  homogeneous  */  //Z=11877
                if ( params.cs==0 )
                {/*5*/  //Z=11878
                    //i = 2;  //Z=11879
                    double z12v[nnmax+1];   // local
                    z12v[0] = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=11880
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11881
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11882
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11883
                        fkv[n] = fkv[n-1]*n;  //Z=11884
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11885
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11886
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11887 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=11888
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=11889
                        /*  longitudinal  */  //Z=11890
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=11891
                            double sump = 0.0;  //Z=11892
                            for ( int ll=0; ll<=m; ll++ )  //Z=11893
                                sump += pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[m-ll][n-m+ll]/(pow(4.0,ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*params.norm);  //Z=11894
                            /* carr1pm[i]:=power(4,m)*power(p21*p21,n-m)*sump/(fkv[n-m]);  //Z=11895 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=11896 */
                            /* i:=i+1;  //Z=11897 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*pow(sqr(params.p21),n-m)*sump/(fkv[n-m]);  //Z=11898
                        }/*7*/  //Z=11899
                        /*  P(q)-coefficients  */  //Z=11900
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=11901
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=11902
                        /*  F(q)-coefficients  */  //Z=11903
                        double sump = 0.0;  //Z=11904
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11905
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=11906
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=11907
                        /*  cross-sectional  */  //Z=11908
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=11909
                        sump = 0.0;  //Z=11910
                        for ( int m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=11911
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=11912
                        if ( search1 )
                        {/*7*/  //Z=11913
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=11914
                                n1 = n;  //Z=11915
                                search1 = false;  //Z=11916
                            }/*8*/  //Z=11917
                        }/*7*/  //Z=11918
                        if ( search4 )
                        {/*7*/  //Z=11919
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=11920
                                n4 = n;  //Z=11921
                                search4 = false;  //Z=11922
                            }/*8*/  //Z=11923
                        }/*7*/  //Z=11924
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=11925
                            if ( n<n1f ) n1f = n;  //Z=11926
                        }/*7*/  //Z=11927
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=11928
                            if ( n<n4f ) n4f = n;  //Z=11929
                        }/*7*/  //Z=11930
                    }/*6*/  //Z=11931
                }/*5*/  /*  of homogeneous cylinder  */  //Z=11932

                /*  core/shell  */  //Z=11934
                if ( params.cs==1 )
                {/*5*/  //Z=11935
                    //i = 2;  //Z=11936
                    double xrmn_n = 1.0;
                    double z12v[nnmax+1];   // local
                    z12v[0] = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=11937
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=11938
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=11939
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=11940
                        fkv[n] = fkv[n-1]*n;  //Z=11941
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=11942
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=11943
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=11944 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=11945
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=11946
                        xrmn_n = -xrmn_n*xrm2z;  //Z=11947
                        pn[n] = pn[n-1]*p*p;  //Z=11948
                        /*  longitudinal  */  //Z=11949
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=11950
                            double sump = 0.0;  //Z=11951
                            for ( int ll=0; ll<=m; ll++ )  //Z=11952
                                sump = sump+pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[m-ll][n-m+ll]/(pow(4.0,ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*params.norm);  //Z=11953
                            /* carr1pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=11954 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=11955 */
                            /* i:=i+1;  //Z=11956 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*sump/(fkv[n-m]);  //Z=11957
                        }/*7*/  //Z=11958
                        /*  P(q)-coefficients  */  //Z=11959
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=11960
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=11961
                        /*  F(q)-coefficients  */  //Z=11962
                        double sump = 0.0;  //Z=11963
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=11964
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=11965
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=11966

                        /*  P(q)-coefficients  */  //Z=11968
                        /*  cross-sectional  */  //Z=11969
                        /*  F121  */  //Z=11970
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=11971
                        /*  F122  */  //Z=11972
                        sump = 0.0;  //Z=11973
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=11974
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11975
                        }/*7*/  //Z=11976
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=11977
                        /*  F123  */  //Z=11978
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=11979

                        /*  F(q)-coefficients  */  //Z=11981
                        /*  cross-sectional  */  //Z=11982
                        /*  F121  */  //Z=11983
                        sump = 0.0;  //Z=11984
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=11985
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11986
                        }/*7*/  //Z=11987
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=11988
                        /*  F122  */  //Z=11989
                        sump = 0.0;  //Z=11990
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=11991
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=11992
                        }/*7*/  //Z=11993
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=11994
                        /*  F123  */  //Z=11995
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=11996

                        if ( search1 )
                        {/*7*/  //Z=11998
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=11999
                                n1 = n;  //Z=12000
                                search1 = false;  //Z=12001
                            }/*8*/  //Z=12002
                        }/*7*/  //Z=12003
                        if ( search4 )
                        {/*7*/  //Z=12004
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12005
                                n4 = n;  //Z=12006
                                search4 = false;  //Z=12007
                            }/*8*/  //Z=12008
                        }/*7*/  //Z=12009
                        if ( search5 )
                        {/*7*/  //Z=12010
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=12011
                                n5 = n;  //Z=12012
                                search5 = false;  //Z=12013
                            }/*8*/  //Z=12014
                        }/*7*/  //Z=12015
                        if ( search6 )
                        {/*7*/  //Z=12016
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=12017
                                n6 = n;  //Z=12018
                                search6 = false;  //Z=12019
                            }/*8*/  //Z=12020
                        }/*7*/  //Z=12021
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12022
                            if ( n<n1f ) n1f = n;  //Z=12023
                        }/*7*/  //Z=12024
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12025
                            if ( n<n4f ) n4f = n;  //Z=12026
                        }/*7*/  //Z=12027
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12028
                            if ( n<n5f ) n5f = n;  //Z=12029
                        }/*7*/  //Z=12030
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12031
                            if ( n<n6f ) n6f = n;  //Z=12032
                        }/*7*/  //Z=12033
                    }/*6*/  /*  n-loop  */  //Z=12034
                }/*5*/  /*  homogeneous loop  */  //Z=12035



                /*  inhomogeneous core/shell  */  //Z=12039
                if ( params.cs==2 )
                {/*5*/  //Z=12040
                    //i = 2;  //Z=12041
                    double xrmn_n = 1.0;
                    double z12v_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12042
                        z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12043
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12044
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=12045
                        fkv[n] = fkv[n-1]*n;  //Z=12046
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12047
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12048
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12049 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12050
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12051
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12052
                        pn[n] = pn[n-1]*p*p;  //Z=12053
                        /*  longitudinal  */  //Z=12054
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12055
                            double sump = 0.0;  //Z=12056
                            for ( int ll=0; ll<=m; ll++ )  //Z=12057
                                sump += pow(sqr(params.p11),ll)*pow(sqr(params.p13),m-ll)*intlar[m-ll][n-m+ll]/(pow(4.0,ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*params.norm);  //Z=12058
                            /* carr1pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12059 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12060 */
                            /* i:=i+1;  //Z=12061 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*sump/(fkv[n-m]);  //Z=12062
                        }/*7*/  //Z=12063
                        /*  P(q)-coefficients  */  //Z=12064
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=12065
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12066
                        /*  F(q)-coefficients  */  //Z=12067
                        double sump = 0.0;  //Z=12068
                        for ( int m=0; m<=n; m++ ) sump += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12069
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12070
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=12071


                        /*  cross-sectional P(q)  */  //Z=12074
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v_n*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12075
                        sump = 0.0;  //Z=12076
                        double sump1 = 0.0;  //Z=12077
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12078
                            const double sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=12079
                            sump += pn[n-m]*sumi;  //Z=12080
                            sump1 += sumi;  //Z=12081
                        }/*7*/  //Z=12082
                        params.CR->carr5p[n] = (1-params.uca/2.0)*z12v_n*xrmn_n*sump;  //Z=12083
                        params.CR->carr6p[n] = (1-params.uca/2.0)*z12v_n*xrn[n]*sump1;  //Z=12084
                        sump = 0.0;  //Z=12085
                        sump1 = 0.0;  //Z=12086
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12087
                            const double sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=12088
                            sump += sumi;  //Z=12089
                            sump1 += pn[n-m]*sumi;  //Z=12090
                        }/*7*/  //Z=12091
                        params.CR->carr7p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump;  //Z=12092
                        params.CR->carr8p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump1;  //Z=12093
                        params.CR->carr9p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v_n*xrn[n]*sump;  //Z=12094

                        /*  F(q)-coefficients  */  //Z=12111
                        params.CR->carr4f[n] = z12v_n*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=12112
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v_n*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12113
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v_n*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12114

                        if ( search1 )
                        {/*7*/  //Z=12116
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12117
                                n1 = n;  //Z=12118
                                search1 = false;  //Z=12119
                            }/*8*/  //Z=12120
                        }/*7*/  //Z=12121
                        if ( search4 )
                        {/*7*/  //Z=12122
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12123
                                n4 = n;  //Z=12124
                                search4 = false;  //Z=12125
                            }/*8*/  //Z=12126
                        }/*7*/  //Z=12127
                        if ( search5 )
                        {/*7*/  //Z=12128
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=12129
                                n5 = n;  //Z=12130
                                search5 = false;  //Z=12131
                            }/*8*/  //Z=12132
                        }/*7*/  //Z=12133
                        if ( search6 )
                        {/*7*/  //Z=12134
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=12135
                                n6 = n;  //Z=12136
                                search6 = false;  //Z=12137
                            }/*8*/  //Z=12138
                        }/*7*/  //Z=12139
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=12140
                            if ( n<n7 ) n7 = n;  //Z=12141
                        }/*7*/  //Z=12142
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=12143
                            if ( n<n8 ) n8 = n;  //Z=12144
                        }/*7*/  //Z=12145
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=12146
                            if ( n<n9 ) n9 = n;  //Z=12147
                        }/*7*/  //Z=12148
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12149
                            if ( n<n1f ) n1f = n;  //Z=12150
                        }/*7*/  //Z=12151
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12152
                            if ( n<n4f ) n4f = n;  //Z=12153
                        }/*7*/  //Z=12154
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12155
                            if ( n<n5f ) n5f = n;  //Z=12156
                        }/*7*/  //Z=12157
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12158
                            if ( n<n6f ) n6f = n;  //Z=12159
                        }/*7*/  //Z=12160
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=12161
                            if ( n<n7f ) n7f = n;  //Z=12162
                        }/*7*/  //Z=12163
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=12164
                            if ( n<n8f ) n8f = n;  //Z=12165
                        }/*7*/  //Z=12166
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=12167
                            if ( n<n9f ) n9f = n;  //Z=12168
                        }/*7*/  //Z=12169
                    }/*6*/   /*  of n-loop  */  //Z=12170
                }/*5*/  /*  of inhomogeneous core/shell  */  //Z=12171

                /*  myelin  */  //Z=12173
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=12174
                    //i = 2;  //Z=12175
                    double z12v[nnmax+1];   // local
                    z12v[0] = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12176
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12177
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12178
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=12179
                        fkv[n] = fkv[n-1]*n;  //Z=12180
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12181
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12182
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12183 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12184
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=12185
                        /*  longitudinal  */  //Z=12186
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12187
                            double sump = 0.0;  //Z=12188
                            for ( int ll=0; ll<=m; ll++ )  //Z=12189
                                sump = sump+pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[m-ll][n-m+ll]/(pow(4.0,ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*params.norm);  //Z=12190
                            /* carr1pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12191 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12192 */
                            /* i:=i+1;  //Z=12193 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*sump/(fkv[n-m]);  //Z=12194
                        }/*7*/  //Z=12195
                        /*  P(q)-coefficients  */  //Z=12196
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=12197
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12198
                        /*  F(q)-coefficients  */  //Z=12199
                        double sump = 0.0;  //Z=12200
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12201
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12202
                        params.CR->carr2f[n] = params.CR->carr2p[n];  //Z=12203
                        /*  cross-sectional  */  //Z=12204
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=12205

                        sump = 0.0;  //Z=12207
                        for ( int m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12208
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=12209
                        if ( search1 )
                        {/*7*/  //Z=12210
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12211
                                n1 = n;  //Z=12212
                                search1 = false;  //Z=12213
                            }/*8*/  //Z=12214
                        }/*7*/  //Z=12215
                        if ( search4 )
                        {/*7*/  //Z=12216
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12217
                                n4 = n;  //Z=12218
                                search4 = false;  //Z=12219
                            }/*8*/  //Z=12220
                        }/*7*/  //Z=12221
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12222
                            if ( n<n1f ) n1f = n;  //Z=12223
                        }/*7*/  //Z=12224
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12225
                            if ( n<n4f ) n4f = n;  //Z=12226
                        }/*7*/  //Z=12227
                    }/*6*/  //Z=12228
                }/*5*/        /*  of myelin  */  //Z=12229
            }/*4*/   /*  of cylinder  */  //Z=12230

            /* if (dim=1) then begin  //Z=12232
                for n:=1 to nmax do begin  //Z=12233
                   z12v[n]:=z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12234
                   b1sv[n]:=b1sv[n-1]*(b1s-1+n);  //Z=12235
                   xln[n]:=-xln[n-1]*xl2z;  //Z=12236
                   xrn[n]:=-xrn[n-1]*xr2z;  //Z=12237
                   carr1[n]:=power(4,2*n)*z12v[n]*xln[n]/((2*n+1)*(n+1));  //Z=12238
                   carr1p[n]:=carr1[n];  //Z=12239
                   carr3[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*b1sv[n]*fkv[n]);  //Z=12240
                   carr3p[n]:=carr3[n];  //Z=12241
                   if search1 then begin  //Z=12242
                      if carr1[n]<1e-50 then begin  //Z=12243
                         n1:=n;  //Z=12244
                         search1:=false;  //Z=12245
                      end;  //Z=12246
                   end;  //Z=12247
                   if search3 then begin  //Z=12248
                      if carr3[n]<1e-50 then begin  //Z=12249
                         n3:=n;  //Z=12250
                         search3:=false;  //Z=12251
                      end;  //Z=12252
                   end;  //Z=12253
                end;  //Z=12254
             end;      */  //Z=12255

            /*  disk form factor coefficients  */  //Z=12257
            if ( dim==2 )
            {/*4*/  //Z=12258
                /*  homogeneous  */  //Z=12259
                if ( params.cs==0 )
                {/*5*/  //Z=12260
                    //int i = 2;  //Z=12261
                    double z12v[nnmax+1];   // local
                    z12v[0] = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12262
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12263
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12264
                        fkv[n] = fkv[n-1]*n;  //Z=12265
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12266
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12267
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12268 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12269
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12270
                        /*  longitudinal  */  //Z=12271
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12272
                            double sump = 0.0;  //Z=12273
                            for ( int ll=0; ll<=m; ll++ )  //Z=12274
                                sump = sump+pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[m-ll][n-m+ll]/(pow(4.0,ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*params.norm);  //Z=12275
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12276 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*sump/(fkv[n-m]);  //Z=12277
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=12278 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=12279
                            /* carr2fm[i]:=carr2pm[i];  //Z=12280 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12281 */
                            //i = i+1;  //Z=12282
                        }/*7*/  //Z=12283
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=12284
                        sump1 = 0.0;  //Z=12285
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12286
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12287
                        }/*7*/  //Z=12288
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12289

                        /*  cross-sectional  */  //Z=12291
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12292
                        double sump = 0.0;  //Z=12293
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12294
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12295
                        }/*7*/  //Z=12296
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=12297

                        /*  series for <...> integration  */  //Z=12299
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12300
                        sump = 0.0;  //Z=12301
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12302
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12303
                        }/*7*/  //Z=12304
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12305

                        if ( search1 )
                        {/*7*/  //Z=12307
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12308
                                n1 = n;  //Z=12309
                                search1 = false;  //Z=12310
                            }/*8*/  //Z=12311
                        }/*7*/  //Z=12312
                        if ( search4 )
                        {/*7*/  //Z=12313
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12314
                                n4 = n;  //Z=12315
                                search4 = false;  //Z=12316
                            }/*8*/  //Z=12317
                        }/*7*/  //Z=12318
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12319
                            if ( n<n1f ) n1f = n;  //Z=12320
                        }/*7*/  //Z=12321
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12322
                            if ( n<n4f ) n4f = n;  //Z=12323
                        }/*7*/  //Z=12324
                    }/*6*/  //Z=12325
                }/*5*/  //Z=12326

                /*  core/shell  */  //Z=12328
                if ( params.cs==1 )
                {/*5*/  //Z=12329
                    //i = 2;  //Z=12330
                    double xrmn_n = 1.0;
                    double z12v[nnmax+1];   // local
                    z12v[0] = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12331
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12332
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12333
                        fkv[n] = fkv[n-1]*n;  //Z=12334
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12335
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12336
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12337 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12338
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12339
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12340
                        pn[n] = pn[n-1]*p*p;  //Z=12341
                        /*  longitudinal  */  //Z=12342
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12343
                            double sump = 0.0;  //Z=12344
                            for ( int ll=0; ll<=m; ll++ )  //Z=12345
                                sump = sump+pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[m-ll][n-m+ll]/(pow(4.0,ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*params.norm);  //Z=12346
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12347 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*sump/(fkv[n-m]);  //Z=12348
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=12349 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=12350
                            /* carr2fm[i]:=carr2pm[i];  //Z=12351 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12352 */
                            /* i:=i+1;  //Z=12353 */
                        }/*7*/  //Z=12354
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12355
                        sump1 = 0.0;  //Z=12356
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12357
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12358
                        }/*7*/  //Z=12359
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12360

                        /*  P(q)-coefficients  */  //Z=12362
                        /*  cross-sectional  */  //Z=12363
                        /*  F121  */  //Z=12364
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=12365
                        /*  F122  */  //Z=12366
                        double sump = 0.0;  //Z=12367
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12368
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12369
                        }/*7*/  //Z=12370
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=12371


                        /*  F122  */  //Z=12374
                        sump = 0.0;  //Z=12375
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12376
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12377
                        }/*7*/  //Z=12378
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=12379


                        /*  F123  */  //Z=12382
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=12383
                        /*  F121  */  //Z=12384
                        sump = 0.0;  //Z=12385
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12386
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12387
                        }/*7*/  //Z=12388
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=12389
                        /*  F122  */  //Z=12390
                        sump = 0.0;  //Z=12391
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12392
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12393
                        }/*7*/  //Z=12394
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=12395
                        /*  F123  */  //Z=12396
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=12397

                        /*  series for <...> integration  */  //Z=12399
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12400
                        sump = 0.0;  //Z=12401
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12402
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12403
                        }/*7*/  //Z=12404
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12405

                        if ( search1 )
                        {/*7*/  //Z=12407
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12408
                                n1 = n;  //Z=12409
                                search1 = false;  //Z=12410
                            }/*8*/  //Z=12411
                        }/*7*/  //Z=12412
                        if ( search4 )
                        {/*7*/  //Z=12413
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12414
                                n4 = n;  //Z=12415
                                search4 = false;  //Z=12416
                            }/*8*/  //Z=12417
                        }/*7*/  //Z=12418
                        if ( search5 )
                        {/*7*/  //Z=12419
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=12420
                                n5 = n;  //Z=12421
                                search5 = false;  //Z=12422
                            }/*8*/  //Z=12423
                        }/*7*/  //Z=12424
                        if ( search6 )
                        {/*7*/  //Z=12425
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=12426
                                n6 = n;  //Z=12427
                                search6 = false;  //Z=12428
                            }/*8*/  //Z=12429
                        }/*7*/  //Z=12430
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12431
                            if ( n<n1f ) n1f = n;  //Z=12432
                        }/*7*/  //Z=12433
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12434
                            if ( n<n4f ) n4f = n;  //Z=12435
                        }/*7*/  //Z=12436
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12437
                            if ( n<n5f ) n5f = n;  //Z=12438
                        }/*7*/  //Z=12439
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12440
                            if ( n<n6f ) n6f = n;  //Z=12441
                        }/*7*/  //Z=12442
                    }/*6*/ /*  of n-loop  */  //Z=12443
                }/*5*/  /*  of core/shell  */  //Z=12444

                /*  inhomogeneous core/shell  */  //Z=12446
                if ( params.cs==2 )
                {/*5*/  //Z=12447
                    //i = 2;  //Z=12448
                    double xrmn_n = 1.0;
                    double z12v_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12449
                        z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12450
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12451
                        fkv[n] = fkv[n-1]*n;  //Z=12452
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12453
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12454
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12455 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12456
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12457
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12458
                        pn[n] = pn[n-1]*p*p;  //Z=12459
                        /*  longitudinal  */  //Z=12460
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12461
                            double sump = 0.0;  //Z=12462
                            for ( int ll=0; ll<=m; ll++ )  //Z=12463
                                sump += pow(params.p11*params.p11,ll)*pow(params.p13*params.p13,m-ll)*intlar[m-ll][n-m+ll]/(pow(4.0,ll)*fk2v[m-ll]*fkv[n-m+ll]*fkv[ll]*params.norm);  //Z=12464
                            /* carr2pm[i]:=power(4,m)*sump/(fkv[n-m]);  //Z=12465 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*sump/(fkv[n-m]);  //Z=12466
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=12467 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=12468
                            /* carr2fm[i]:=carr2pm[i];  //Z=12469 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=12470 */
                            /* i:=i+1;  //Z=12471 */
                        }/*7*/  //Z=12472
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12473
                        double sump1 = 0.0;  //Z=12474
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12475
                            sump1 += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12476
                        }/*7*/  //Z=12477
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1;  //Z=12478

                        /*  cross-sectional P(q)  */  //Z=12480
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v_n*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=12481
                        double sump = 0.0;  //Z=12482
                        sump1 = 0.0;  //Z=12483
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12484
                            const double sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=12485
                            sump += pn[n-m]*sumi;  //Z=12486
                            sump1 += sumi;  //Z=12487
                        }/*7*/  //Z=12488
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v_n*xrmn_n*sump;  //Z=12489
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v_n*xrn[n]*sump1;  //Z=12490
                        sump = 0.0;  //Z=12491
                        sump1 = 0.0;  //Z=12492
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12493
                            const double sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=12494
                            sump += sumi;  //Z=12495
                            sump1 += pn[n-m]*sumi;  //Z=12496
                        }/*7*/  //Z=12497
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrmn_n*sump;  //Z=12498
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrmn_n*sump1;  //Z=12499
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrn[n]*sump;  //Z=12500

                        /*  (* P(q)-coefficients *)  //Z=12502
                            (* cross-sectional *)  //Z=12503
                            (* F121 *)  //Z=12504
                            carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=12505
                            (* F122 *)  //Z=12506
                            sump:=0.0;  //Z=12507
                               for m:=0 to n do begin  //Z=12508
                               sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12509
                            end;  //Z=12510
                            carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=12511
                            (* F122 *)  //Z=12512
                            sump:=0.0;  //Z=12513
                            for m:=0 to n do begin  //Z=12514
                               sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12515
                            end;  //Z=12516
                            carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=12517
                            (* F123 *)  //Z=12518
                            carr6p[n]:=carr4p[n]/pn[n];      */  //Z=12519

                        /*  F(q)  */  //Z=12521
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v_n*xrn[n]/(gam3[n]*fkv[n]);  //Z=12522
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v_n*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12523
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v_n*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=12524

                        /*  series for <...> integration  */  //Z=12526
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=12527
                        sump = 0.0;  //Z=12528
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12529
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12530
                        }/*7*/  //Z=12531
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=12532

                        if ( search1 )
                        {/*7*/  //Z=12534
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12535
                                n1 = n;  //Z=12536
                                search1 = false;  //Z=12537
                            }/*8*/  //Z=12538
                        }/*7*/  //Z=12539
                        if ( search4 )
                        {/*7*/  //Z=12540
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12541
                                n4 = n;  //Z=12542
                                search4 = false;  //Z=12543
                            }/*8*/  //Z=12544
                        }/*7*/  //Z=12545
                        if ( search5 )
                        {/*7*/  //Z=12546
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=12547
                                n5 = n;  //Z=12548
                                search5 = false;  //Z=12549
                            }/*8*/  //Z=12550
                        }/*7*/  //Z=12551
                        if ( search6 )
                        {/*7*/  //Z=12552
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=12553
                                n6 = n;  //Z=12554
                                search6 = false;  //Z=12555
                            }/*8*/  //Z=12556
                        }/*7*/  //Z=12557
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=12558
                            if ( n<n7 ) n7 = n;  //Z=12559
                        }/*7*/  //Z=12560
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=12561
                            if ( n<n8 ) n8 = n;  //Z=12562
                        }/*7*/  //Z=12563
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=12564
                            if ( n<n9 ) n9 = n;  //Z=12565
                        }/*7*/  //Z=12566
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12567
                            if ( n<n1f ) n1f = n;  //Z=12568
                        }/*7*/  //Z=12569
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12570
                            if ( n<n4f ) n4f = n;  //Z=12571
                        }/*7*/  //Z=12572
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12573
                            if ( n<n5f ) n5f = n;  //Z=12574
                        }/*7*/  //Z=12575
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12576
                            if ( n<n6f ) n6f = n;  //Z=12577
                        }/*7*/  //Z=12578
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=12579
                            if ( n<n7f ) n7f = n;  //Z=12580
                        }/*7*/  //Z=12581
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=12582
                            if ( n<n8f ) n8f = n;  //Z=12583
                        }/*7*/  //Z=12584
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=12585
                            if ( n<n9f ) n9f = n;  //Z=12586
                        }/*7*/  //Z=12587
                    }/*6*/ /*  of n-loop  */  //Z=12588
                }/*5*/  /*  of inhomogeneous core/shell  */  //Z=12589



            }/*4*/   /*  of disk  */  //Z=12593

           /* i:=2;   (* first element: carr2[1]:=1 for n=0, m=0; *)  //Z=12598
               end;  (* of n-loop *)    */  //Z=12709
        }/*3*/   /*  of cho1=1  */  //Z=12710

        /* ** general orientation case ** */  //Z=12712
        /*  for phi<>0, too slow, only for cylinders  */  //Z=12713
        if ( params.orcase==5 )
        {/*3*/  //Z=12714

            const double phi = params.polPhi;
            const double theta = params.polTheta;

            params.p11 = -cos(phi*M_PI/180.0)*cos(theta*M_PI/180.0);  //Z=12716
            params.p12 = sin(phi*M_PI/180.0);  //Z=12717
            params.p13 = cos(phi*M_PI/180.0)*sin(theta*M_PI/180.0);  //Z=12718
            params.p21 = -cos(phi*M_PI/180.0);  //Z=12719
            params.p22 = -sin(phi*M_PI/180.0)*cos(theta*M_PI/180.0);  //Z=12720
            params.p23 = sin(phi*M_PI/180.0)*sin(theta*M_PI/180.0);  //Z=12721
            params.p31 = -sin(theta*M_PI/180.0);  //Z=12722
            params.p32 = 0;  //Z=12723
            params.p33 = -cos(theta*M_PI/180.0);  //Z=12724

            for ( n=1; n<=2*nmax; n++ ) fkv[n] = fkv[n-1]*n;  //Z=12726

            //i = 2;  //Z=12728
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=12729
                for ( int m=0; m<=2*n; m++ )
                {/*5*/  //Z=12730
                    qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,5,0,m,2*n-m,params.CR->carr1p,intl);
                    /* carr1pm[i]:=intl/(fkv[m]*fkv[2*n-m]*norm);  //Z=12732 */
                    //i = i+1;  //Z=12733
                }/*5*/  //Z=12734
            }/*4*/  //Z=12735

            double b1sv_n = 1.0;
            double z12v_n = 1.0;
            for ( n=1; n<=nmax; n++ )
            {/*4*/  //Z=12737
                z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12738
                b1sv_n = b1sv_n*(b1s-1+n);  //Z=12739
                xln[n] = -xln[n-1]*xl2z;  //Z=12740
                xrn[n] = -xrn[n-1]*xr2z;  //Z=12741
                params.CR->carr1p[n] = pow(4.0,2*n)*z12v_n*xln[n]/((2.0*n+1)*(n+1));  //Z=12742
                params.CR->carr3p[n] = 4*(n+1/2.0)*fk2v[n]*z12v_n*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*b1sv_n*fkv[n]);  //Z=12743
                if ( search1 )
                {/*5*/  //Z=12744
                    if ( params.CR->carr1p[n]<1e-50 )
                    {/*6*/  //Z=12745
                        n1 = n;  //Z=12746
                        search1 = false;  //Z=12747
                    }/*6*/  //Z=12748
                }/*5*/  //Z=12749
                if ( search3 )
                {/*5*/  //Z=12750
                    if ( params.CR->carr3p[n]<1e-50 )
                    {/*6*/  //Z=12751
                        n3 = n;  //Z=12752
                        search3 = false;  //Z=12753
                    }/*6*/  //Z=12754
                }/*5*/  //Z=12755
            }/*4*/  //Z=12756
        }/*3*/   /*  of cho1=5  */  //Z=12757




        /* ** x-axis ** */  //Z=12762
        if ( (params.orcase==2) && (dim!=3) )
        {/*3*/  //Z=12763
            /* ** cylinders ** */  //Z=12764
            if ( dim==1 )
            {/*4*/  //Z=12765
                double z12v[nnmax+1];   // local
                z12v[0] = 1.0;
                /*  homogeneous  */  //Z=12766
                if ( params.cs==0 )
                {/*5*/  //Z=12767
                    //i = 2;  //Z=12768
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12769
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12770
                        z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12771
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=12772
                        fkv[n] = fkv[n-1]*n;  //Z=12773
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12774
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12775
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12776 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12777
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12778
                        /*  longitudinal  */  //Z=12779
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12780
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=12781
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=12782 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=12783
                            /* carr1fm[i]:=carr1pm[i];  //Z=12784 */
                            /* i:=i+1;  //Z=12785 */
                        }/*7*/  //Z=12786
                        /*  P(q)-coefficient  */  //Z=12787
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=12788
                        /*  F(q)-coefficient  */  //Z=12789
                        double sump = 0.0;  //Z=12790
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12791
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12792

                        /*  cross-section  */  //Z=12794
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=12795
                        sump = 0.0;  //Z=12796
                        for ( int m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=12797
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=12798


                        if ( search1 )
                        {/*7*/  //Z=12801
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12802
                                n1 = n;  //Z=12803
                                search1 = false;  //Z=12804
                            }/*8*/  //Z=12805
                        }/*7*/  //Z=12806
                        if ( search4 )
                        {/*7*/  //Z=12807
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12808
                                n4 = n;  //Z=12809
                                search4 = false;  //Z=12810
                            }/*8*/  //Z=12811
                        }/*7*/  //Z=12812
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12813
                            if ( n<n1f ) n1f = n;  //Z=12814
                        }/*7*/  //Z=12815
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12816
                            if ( n<n4f ) n4f = n;  //Z=12817
                        }/*7*/  //Z=12818
                    }/*6*/  /*  of n-loop  */  //Z=12819
                }/*5*/  /*  of cs=0  */  //Z=12820

                /*  core/shell  */  //Z=12822
                if ( params.cs==1 )
                {/*5*/  //Z=12823
                    //i = 2;  //Z=12824
                    double xrmn_n = 1.0;
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12825
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12826
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12827
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=12828
                        fkv[n] = fkv[n-1]*n;  //Z=12829
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12830
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12831
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12832 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12833
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12834
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12835
                        pn[n] = pn[n-1]*p*p;  //Z=12836
                        /*  longitudinal  */  //Z=12837
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12838
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=12839
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=12840 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=12841
                            /* carr1fm[i]:=carr1pm[i];  //Z=12842 */
                            /* i:=i+1;  //Z=12843 */
                        }/*7*/  //Z=12844
                        /*  P(q)-coefficient  */  //Z=12845
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=12846
                        /*  F(q)-coefficient  */  //Z=12847
                        double sump = 0.0;  //Z=12848
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12849
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12850

                        /*  P(q)-coefficients  */  //Z=12852
                        /*  cross-sectional  */  //Z=12853
                        /*  F121  */  //Z=12854
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12855
                        /*  F122  */  //Z=12856
                        sump = 0.0;  //Z=12857
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12858
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12859
                        }/*7*/  //Z=12860
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=12861
                        /*  F123  */  //Z=12862
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=12863

                        /*  F(q)-coefficients  */  //Z=12865
                        /*  cross-sectional  */  //Z=12866
                        /*  F121  */  //Z=12867
                        sump = 0.0;  //Z=12868
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12869
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12870
                        }/*7*/  //Z=12871
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=12872
                        /*  F122  */  //Z=12873
                        sump = 0.0;  //Z=12874
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12875
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12876
                        }/*7*/  //Z=12877
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=12878
                        /*  F123  */  //Z=12879
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=12880


                        if ( search1 )
                        {/*7*/  //Z=12883
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12884
                                n1 = n;  //Z=12885
                                search1 = false;  //Z=12886
                            }/*8*/  //Z=12887
                        }/*7*/  //Z=12888
                        if ( search4 )
                        {/*7*/  //Z=12889
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=12890
                                n4 = n;  //Z=12891
                                search4 = false;  //Z=12892
                            }/*8*/  //Z=12893
                        }/*7*/  //Z=12894
                        if ( search5 )
                        {/*7*/  //Z=12895
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=12896
                                n5 = n;  //Z=12897
                                search5 = false;  //Z=12898
                            }/*8*/  //Z=12899
                        }/*7*/  //Z=12900
                        if ( search6 )
                        {/*7*/  //Z=12901
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=12902
                                n6 = n;  //Z=12903
                                search6 = false;  //Z=12904
                            }/*8*/  //Z=12905
                        }/*7*/  //Z=12906
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=12907
                            if ( n<n1f ) n1f = n;  //Z=12908
                        }/*7*/  //Z=12909
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=12910
                            if ( n<n4f ) n4f = n;  //Z=12911
                        }/*7*/  //Z=12912
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=12913
                            if ( n<n5f ) n5f = n;  //Z=12914
                        }/*7*/  //Z=12915
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=12916
                            if ( n<n6f ) n6f = n;  //Z=12917
                        }/*7*/  //Z=12918
                    }/*6*/  /*  of n-loop  */  //Z=12919
                }/*5*/  /*  of cs=1  */  //Z=12920

                /*  inhomogeneous core/shell  */  //Z=12922
                if ( params.cs==2 )
                {/*5*/  //Z=12923
                    //i = 2;  //Z=12924
                    double xrmn_n = 1.0;
                    double z12v_n = 1.0;
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=12925
                        z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=12926
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=12927
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=12928
                        fkv[n] = fkv[n-1]*n;  //Z=12929
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=12930
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=12931
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=12932 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=12933
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=12934
                        xrmn_n = -xrmn_n*xrm2z;  //Z=12935
                        pn[n] = pn[n-1]*p*p;  //Z=12936
                        /*  longitudinal  */  //Z=12937
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12938
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=12939
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=12940 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=12941
                            /* carr1fm[i]:=carr1pm[i];  //Z=12942 */
                            /* i:=i+1;  //Z=12943 */
                        }/*7*/  //Z=12944
                        /*  P(q)-coefficient  */  //Z=12945
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=12946
                        /*  F(q)-coefficient  */  //Z=12947
                        double sump = 0.0;  //Z=12948
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=12949
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=12950


                        /*  cross-sectional P(q)  */  //Z=12953
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v_n*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12954
                        sump = 0.0;  //Z=12955
                        sump1 = 0.0;  //Z=12956
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12957
                            const double sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=12958
                            sump += pn[n-m]*sumi;  //Z=12959
                            sump1 += sumi;  //Z=12960
                        }/*7*/  //Z=12961
                        params.CR->carr5p[n] = (1-params.uca/2.0)*z12v_n*xrmn_n*sump;  //Z=12962
                        params.CR->carr6p[n] = (1-params.uca/2.0)*z12v_n*xrn[n]*sump1;  //Z=12963
                        sump = 0.0;  //Z=12964
                        sump1 = 0.0;  //Z=12965
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=12966
                            const double sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=12967
                            sump += sumi;  //Z=12968
                            sump1 += pn[n-m]*sumi;  //Z=12969
                        }/*7*/  //Z=12970
                        params.CR->carr7p[n] = sqr(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump;  //Z=12971
                        params.CR->carr8p[n] = sqr(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump1;  //Z=12972
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v_n*xrn[n]*sump;  //Z=12973

                        /* (* P(q)-coefficients *)  //Z=12975
                        (* cross-sectional *)  //Z=12976
                        (* F121 *)  //Z=12977
                        carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=12978
                        (* F122 *)  //Z=12979
                        sump:=0.0;  //Z=12980
                           for m:=0 to n do begin  //Z=12981
                           sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=12982
                        end;  //Z=12983
                        carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=12984
                        (* F123 *)  //Z=12985
                        carr6p[n]:=carr4p[n]/pn[n];   */  //Z=12986

                        /*  F(q)-coefficients  */  //Z=12988
                        /*  cross-sectional  */  //Z=12989
                        params.CR->carr4f[n] = z12v_n*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=12990
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v_n*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12991
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v_n*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=12992


                        if ( search1 )
                        {/*7*/  //Z=12995
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=12996
                                n1 = n;  //Z=12997
                                search1 = false;  //Z=12998
                            }/*8*/  //Z=12999
                        }/*7*/  //Z=13000
                        if ( search4 )
                        {/*7*/  //Z=13001
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13002
                                n4 = n;  //Z=13003
                                search4 = false;  //Z=13004
                            }/*8*/  //Z=13005
                        }/*7*/  //Z=13006
                        if ( search5 )
                        {/*7*/  //Z=13007
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=13008
                                n5 = n;  //Z=13009
                                search5 = false;  //Z=13010
                            }/*8*/  //Z=13011
                        }/*7*/  //Z=13012
                        if ( search6 )
                        {/*7*/  //Z=13013
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=13014
                                n6 = n;  //Z=13015
                                search6 = false;  //Z=13016
                            }/*8*/  //Z=13017
                        }/*7*/  //Z=13018
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=13019
                            if ( n<n7 ) n7 = n;  //Z=13020
                        }/*7*/  //Z=13021
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=13022
                            if ( n<n8 ) n8 = n;  //Z=13023
                        }/*7*/  //Z=13024
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=13025
                            if ( n<n9 ) n9 = n;  //Z=13026
                        }/*7*/  //Z=13027
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13028
                            if ( n<n1f ) n1f = n;  //Z=13029
                        }/*7*/  //Z=13030
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13031
                            if ( n<n4f ) n4f = n;  //Z=13032
                        }/*7*/  //Z=13033
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13034
                            if ( n<n5f ) n5f = n;  //Z=13035
                        }/*7*/  //Z=13036
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13037
                            if ( n<n6f ) n6f = n;  //Z=13038
                        }/*7*/  //Z=13039
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=13040
                            if ( n<n7f ) n7f = n;  //Z=13041
                        }/*7*/  //Z=13042
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=13043
                            if ( n<n8f ) n8f = n;  //Z=13044
                        }/*7*/  //Z=13045
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=13046
                            if ( n<n9f ) n9f = n;  //Z=13047
                        }/*7*/  //Z=13048
                    }/*6*/  /*  of n-loop  */  //Z=13049
                }/*5*/  /*  of cs=1  */  //Z=13050

                /*  myelin  */  //Z=13052
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=13053
                    //i = 2;  //Z=13054
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13055
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13056
                        z12vl[n] = z12vl[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13057
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=13058
                        fkv[n] = fkv[n-1]*n;  //Z=13059
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13060
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13061
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13062 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13063
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=13064
                        /*  longitudinal  */  //Z=13065
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13066
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13067
                            /* carr1pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=13068 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13069
                            /* carr1fm[i]:=carr1pm[i];  //Z=13070 */
                            /* i:=i+1;  //Z=13071 */
                        }/*7*/  //Z=13072
                        /*  P(q)-coefficient  */  //Z=13073
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=13074
                        /*  F(q)-coefficient  */  //Z=13075
                        double sump = 0.0;  //Z=13076
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13077
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump/(pow(4.0,n));  //Z=13078

                        /*  cross-section  */  //Z=13080
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=13081

                        sump = 0.0;  //Z=13083
                        for ( int m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13084
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13085


                        if ( search1 )
                        {/*7*/  //Z=13088
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13089
                                n1 = n;  //Z=13090
                                search1 = false;  //Z=13091
                            }/*8*/  //Z=13092
                        }/*7*/  //Z=13093
                        if ( search4 )
                        {/*7*/  //Z=13094
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13095
                                n4 = n;  //Z=13096
                                search4 = false;  //Z=13097
                            }/*8*/  //Z=13098
                        }/*7*/  //Z=13099
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13100
                            if ( n<n1f ) n1f = n;  //Z=13101
                        }/*7*/  //Z=13102
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13103
                            if ( n<n4f ) n4f = n;  //Z=13104
                        }/*7*/  //Z=13105
                    }/*6*/  /*  of n-loop  */  //Z=13106
                }/*5*/  /*  of cs=3  */  //Z=13107
            }/*4*/  /*  of cylinders  */  //Z=13108

            /* ** disks ** */  //Z=13110
            if ( dim==2 )
            {/*4*/  //Z=13111
                double z12v[nnmax+1];   // local
                z12v[0] = 1.0;
                /*  homogeneous  */  //Z=13112
                if ( params.cs==0 )
                {/*5*/  //Z=13113
                    //i = 2;  //Z=13114
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13115
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13116
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13117
                        fkv[n] = fkv[n-1]*n;  //Z=13118
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13119
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13120
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13121 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13122
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13123
                        /*  longitudinal  */  //Z=13124
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13125
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13126
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=13127 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13128
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=13129 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=13130
                            /* carr2fm[i]:=carr2pm[i];  //Z=13131 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13132 */
                            /* i:=i+1;  //Z=13133 */
                        }/*7*/  //Z=13134
                        /*  P(q)-coefficient  */  //Z=13135
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13136

                        /*  F(q)-coefficient  */  //Z=13138
                        double sump = 0.0;  //Z=13139
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13140
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13141
                        }/*7*/  //Z=13142
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13143

                        /*  cross-section  */  //Z=13145
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13146
                        sump = 0.0;  //Z=13147
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13148
                            sump += z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13149
                        }/*7*/  //Z=13150
                        params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=13151

                        /*  series for <...> integration  */  //Z=13153
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13154
                        sump = 0.0;  //Z=13155
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13156
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13157
                        }/*7*/  //Z=13158
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13159

                        if ( search1 )
                        {/*7*/  //Z=13161
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13162
                                n1 = n;  //Z=13163
                                search1 = false;  //Z=13164
                            }/*8*/  //Z=13165
                        }/*7*/  //Z=13166
                        if ( search4 )
                        {/*7*/  //Z=13167
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13168
                                n4 = n;  //Z=13169
                                search4 = false;  //Z=13170
                            }/*8*/  //Z=13171
                        }/*7*/  //Z=13172
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13173
                            if ( n<n1f ) n1f = n;  //Z=13174
                        }/*7*/  //Z=13175
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13176
                            if ( n<n4f ) n4f = n;  //Z=13177
                        }/*7*/  //Z=13178
                    }/*6*/ /*  n-loop  */  //Z=13179
                }/*5*/  /*  of cs=0  */  //Z=13180

                /*  core/shell  */  //Z=13182
                if ( params.cs==1 )
                {/*5*/  //Z=13183
                    //i = 2;  //Z=13184
                    double xrmn_n = 1.0;
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13185
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13186
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13187
                        fkv[n] = fkv[n-1]*n;  //Z=13188
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13189
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13190
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13191 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13192
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13193
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13194
                        pn[n] = pn[n-1]*p*p;  //Z=13195
                        /*  longitudinal  */  //Z=13196
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13197
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13198
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=13199 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13200
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=13201 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=13202
                            /* carr2fm[i]:=carr2pm[i];  //Z=13203 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13204 */
                            /* i:=i+1;  //Z=13205 */
                        }/*7*/  //Z=13206
                        /*  P(q)-coefficient  */  //Z=13207
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13208
                        /*  F(q)-coefficient  */  //Z=13209
                        double sump = 0.0;  //Z=13210
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13211
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13212
                        }/*7*/  //Z=13213
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13214

                        /*  F121  */  //Z=13216
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13217
                        /*  F122  */  //Z=13218
                        sump = 0.0;  //Z=13219
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13220
                            sump += pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13221
                        }/*7*/  //Z=13222
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=13223

                        /*  F122  */  //Z=13225
                        sump = 0.0;  //Z=13226
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13227
                            sump += pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13228
                        }/*7*/  //Z=13229
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=13230


                        /*  F123  */  //Z=13233
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=13234
                        /*  F121  */  //Z=13235
                        sump = 0.0;  //Z=13236
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13237
                            sump += z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13238
                        }/*7*/  //Z=13239
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=13240
                        /*  F122  */  //Z=13241
                        sump = 0.0;  //Z=13242
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13243
                            sump += pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13244
                        }/*7*/  //Z=13245
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=13246
                        /*  F123  */  //Z=13247
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=13248

                        /*  series for <...> integration  */  //Z=13250
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13251
                        sump = 0.0;  //Z=13252
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13253
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13254
                        }/*7*/  //Z=13255
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13256

                        if ( search1 )
                        {/*7*/  //Z=13258
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13259
                                n1 = n;  //Z=13260
                                search1 = false;  //Z=13261
                            }/*8*/  //Z=13262
                        }/*7*/  //Z=13263
                        if ( search4 )
                        {/*7*/  //Z=13264
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13265
                                n4 = n;  //Z=13266
                                search4 = false;  //Z=13267
                            }/*8*/  //Z=13268
                        }/*7*/  //Z=13269
                        if ( search5 )
                        {/*7*/  //Z=13270
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=13271
                                n5 = n;  //Z=13272
                                search5 = false;  //Z=13273
                            }/*8*/  //Z=13274
                        }/*7*/  //Z=13275
                        if ( search6 )
                        {/*7*/  //Z=13276
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=13277
                                n6 = n;  //Z=13278
                                search6 = false;  //Z=13279
                            }/*8*/  //Z=13280
                        }/*7*/  //Z=13281
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13282
                            if ( n<n1f ) n1f = n;  //Z=13283
                        }/*7*/  //Z=13284
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13285
                            if ( n<n4f ) n4f = n;  //Z=13286
                        }/*7*/  //Z=13287
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13288
                            if ( n<n5f ) n5f = n;  //Z=13289
                        }/*7*/  //Z=13290
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13291
                            if ( n<n6f ) n6f = n;  //Z=13292
                        }/*7*/  //Z=13293
                    }/*6*/ /*  n-loop  */  //Z=13294
                }/*5*/  /*  of cs=1  */  //Z=13295

                /*  inhomogeneous core/shell  */  //Z=13297
                if ( params.cs==2 )
                {/*5*/  //Z=13298
                    //i = 2;  //Z=13299
                    double xrmn_n = 1.0;
                    double z12v_n = 1.0;
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13300
                        z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13301
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13302
                        fkv[n] = fkv[n-1]*n;  //Z=13303
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13304
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13305
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13306 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13307
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13308
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13309
                        pn[n] = pn[n-1]*p*p;  //Z=13310
                        /*  longitudinal  */  //Z=13311
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13312
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,2,0,2*m,2*n-2*m,params.CR->carr1p,intl);  //Z=13313
                            /* carr2pm[i]:=power(4,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*norm);  //Z=13314 */
                            params.CR->carr11pm[n][m] = pow(4.0,m)*intl/(fk2v[m]*fkv[n-m]*fkv[n-m]*params.norm);  //Z=13315
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(power(4,m)*fkv[m]*fkv[n-m]);  //Z=13316 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(pow(4.0,m)*fkv[m]*fkv[n-m]);  //Z=13317
                            /* carr2fm[i]:=carr2pm[i];  //Z=13318 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13319 */
                            /* i:=i+1;  //Z=13320 */
                        }/*7*/  //Z=13321
                        /*  P(q)-coefficient  */  //Z=13322
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13323
                        /*  F(q)-coefficient  */  //Z=13324
                        double sump = 0.0;  //Z=13325
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13326
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13327
                        }/*7*/  //Z=13328
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13329

                        /*  cross-sectional P(q)  */  //Z=13331
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v_n*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=13332
                        sump = 0.0;  //Z=13333
                        sump1 = 0.0;  //Z=13334
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13335
                            const double sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=13336
                            sump += pn[n-m]*sumi;  //Z=13337
                            sump1 += sumi;  //Z=13338
                        }/*7*/  //Z=13339
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v_n*xrmn_n*sump;  //Z=13340
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v_n*xrn[n]*sump1;  //Z=13341
                        sump = 0.0;  //Z=13342
                        sump1 = 0.0;  //Z=13343
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13344
                            const double sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=13345
                            sump += sumi;  //Z=13346
                            sump1 += pn[n-m]*sumi;  //Z=13347
                        }/*7*/  //Z=13348
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrmn_n*sump;  //Z=13349
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrmn_n*sump1;  //Z=13350
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v_n*xrn[n]*sump;  //Z=13351


                        /* (* F121 *)  //Z=13354
                         carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=13355
                         (* F122 *)  //Z=13356
                         sump:=0.0;  //Z=13357
                            for m:=0 to n do begin  //Z=13358
                            sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13359
                         end;  //Z=13360
                         carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=13361
                         (* F122 *)  //Z=13362
                         sump:=0.0;  //Z=13363
                            for m:=0 to n do begin  //Z=13364
                            sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13365
                         end;  //Z=13366
                         carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=13367
                         (* F123 *)  //Z=13368
                         carr6p[n]:=carr4p[n]/pn[n];   */  //Z=13369

                        /*  F(q)   */  //Z=13371
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v_n*xrn[n]/(gam3[n]*fkv[n]);  //Z=13372
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v_n*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=13373
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v_n*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=13374

                        /*  series for <...> integration  */  //Z=13376
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13377
                        sump = 0.0;  //Z=13378
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13379
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13380
                        }/*7*/  //Z=13381
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13382

                        if ( search1 )
                        {/*7*/  //Z=13384
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13385
                                n1 = n;  //Z=13386
                                search1 = false;  //Z=13387
                            }/*8*/  //Z=13388
                        }/*7*/  //Z=13389
                        if ( search4 )
                        {/*7*/  //Z=13390
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13391
                                n4 = n;  //Z=13392
                                search4 = false;  //Z=13393
                            }/*8*/  //Z=13394
                        }/*7*/  //Z=13395
                        if ( search5 )
                        {/*7*/  //Z=13396
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=13397
                                n5 = n;  //Z=13398
                                search5 = false;  //Z=13399
                            }/*8*/  //Z=13400
                        }/*7*/  //Z=13401
                        if ( search6 )
                        {/*7*/  //Z=13402
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=13403
                                n6 = n;  //Z=13404
                                search6 = false;  //Z=13405
                            }/*8*/  //Z=13406
                        }/*7*/  //Z=13407
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=13408
                            if ( n<n7 ) n7 = n;  //Z=13409
                        }/*7*/  //Z=13410
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=13411
                            if ( n<n8 ) n8 = n;  //Z=13412
                        }/*7*/  //Z=13413
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=13414
                            if ( n<n9 ) n9 = n;  //Z=13415
                        }/*7*/  //Z=13416
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13417
                            if ( n<n1f ) n1f = n;  //Z=13418
                        }/*7*/  //Z=13419
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13420
                            if ( n<n4f ) n4f = n;  //Z=13421
                        }/*7*/  //Z=13422
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13423
                            if ( n<n5f ) n5f = n;  //Z=13424
                        }/*7*/  //Z=13425
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13426
                            if ( n<n6f ) n6f = n;  //Z=13427
                        }/*7*/  //Z=13428
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=13429
                            if ( n<n7f ) n7f = n;  //Z=13430
                        }/*7*/  //Z=13431
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=13432
                            if ( n<n8f ) n8f = n;  //Z=13433
                        }/*7*/  //Z=13434
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=13435
                            if ( n<n9f ) n9f = n;  //Z=13436
                        }/*7*/  //Z=13437
                    }/*6*/ /*  n-loop  */  //Z=13438
                }/*5*/  /*  of inhomogeneous core/shell  */  //Z=13439

            }/*4*/  /*  of disk  */  //Z=13441
        }/*3*/  //Z=13442



        /* ** y-axis ** */  //Z=13446
        if ( (params.orcase==3) && (dim!=3) )
        {/*3*/  //Z=13447

            /* ** cylinder ** */  //Z=13449
            if ( dim==1 )
            {/*4*/  //Z=13450
                double z12v[nnmax+1];   // local
                z12v[0] = 1.0;
                /*  homogeneous  */  //Z=13451
                if ( params.cs==0 )
                {/*5*/  //Z=13452
                    //i = 2;  //Z=13453
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13454
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13455
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13456
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=13457
                        fkv[n] = fkv[n-1]*n;  //Z=13458
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13459
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13460
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13461 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13462
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13463
                        /*  longitudinal  */  //Z=13464
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13465
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13466
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=13467 */
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13468
                            /* carr1fm[i]:=carr1pm[i];  //Z=13469 */
                            /* i:=i+1;  //Z=13470 */
                        }/*7*/  //Z=13471
                        /*  P(q)-coefficient  */  //Z=13472
                        params.CR->carr1p[n] = pow(4.0,2*n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=13473
                        /*  F(q)-coefficient  */  //Z=13474
                        double sump = 0.0;  //Z=13475
                        for ( int m=0; m<=n; m++ ) sump += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13476
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=13477

                        /*  cross-section  */  //Z=13479
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=13480
                        sump = 0.0;  //Z=13481
                        for ( int m=0; m<=n; m++ ) sump += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13482
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13483

                        if ( search1 )
                        {/*7*/  //Z=13485
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13486
                                n1 = n;  //Z=13487
                                search1 = false;  //Z=13488
                            }/*8*/  //Z=13489
                        }/*7*/  //Z=13490
                        if ( search4 )
                        {/*7*/  //Z=13491
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13492
                                n4 = n;  //Z=13493
                                search4 = false;  //Z=13494
                            }/*8*/  //Z=13495
                        }/*7*/  //Z=13496
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13497
                            if ( n<n1f ) n1f = n;  //Z=13498
                        }/*7*/  //Z=13499
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13500
                            if ( n<n4f ) n4f = n;  //Z=13501
                        }/*7*/  //Z=13502
                    }/*6*/  /*  of n-loop  */  //Z=13503
                }/*5*/  /*  cs=0  */  //Z=13504

                /*  core/shell  */  //Z=13506
                if ( params.cs==1 )
                {/*5*/  //Z=13507
                    //i = 2;  //Z=13508
                    double xrmn_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13509
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13510
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13511
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=13512
                        fkv[n] = fkv[n-1]*n;  //Z=13513
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13514
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13515
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13516 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13517
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13518
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13519
                        pn[n] = pn[n-1]*p*p;  //Z=13520
                        /*  longitudinal  */  //Z=13521
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13522
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13523
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=13524 */
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13525
                            /* carr1fm[i]:=carr1pm[i];  //Z=13526 */
                            /* i:=i+1;  //Z=13527 */
                        }/*7*/  //Z=13528
                        /*  P(q)-coefficient  */  //Z=13529
                        params.CR->carr1p[n] = pow(4.0,2*n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=13530
                        /*  F(q)-coefficient  */  //Z=13531
                        double sump = 0.0;  //Z=13532
                        for ( int m=0; m<=n; m++ ) sump += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13533
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=13534

                        /*  P(q)-coefficients  */  //Z=13536
                        /*  cross-sectional  */  //Z=13537
                        /*  F121  */  //Z=13538
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13539
                        /*  F122  */  //Z=13540
                        sump = 0.0;  //Z=13541
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13542
                            sump += pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13543
                        }/*7*/  //Z=13544
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=13545
                        /*  F123  */  //Z=13546
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=13547

                        /*  F(q)-coefficients  */  //Z=13549
                        /*  cross-sectional  */  //Z=13550
                        /*  F121  */  //Z=13551
                        sump = 0.0;  //Z=13552
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13553
                            sump += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13554
                        }/*7*/  //Z=13555
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13556
                        /*  F122  */  //Z=13557
                        sump = 0.0;  //Z=13558
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13559
                            sump += pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13560
                        }/*7*/  //Z=13561
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=13562
                        /*  F123  */  //Z=13563
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=13564

                        if ( search1 )
                        {/*7*/  //Z=13566
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13567
                                n1 = n;  //Z=13568
                                search1 = false;  //Z=13569
                            }/*8*/  //Z=13570
                        }/*7*/  //Z=13571
                        if ( search4 )
                        {/*7*/  //Z=13572
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13573
                                n4 = n;  //Z=13574
                                search4 = false;  //Z=13575
                            }/*8*/  //Z=13576
                        }/*7*/  //Z=13577
                        if ( search5 )
                        {/*7*/  //Z=13578
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=13579
                                n5 = n;  //Z=13580
                                search5 = false;  //Z=13581
                            }/*8*/  //Z=13582
                        }/*7*/  //Z=13583
                        if ( search6 )
                        {/*7*/  //Z=13584
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=13585
                                n6 = n;  //Z=13586
                                search6 = false;  //Z=13587
                            }/*8*/  //Z=13588
                        }/*7*/  //Z=13589
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13590
                            if ( n<n1f ) n1f = n;  //Z=13591
                        }/*7*/  //Z=13592
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13593
                            if ( n<n4f ) n4f = n;  //Z=13594
                        }/*7*/  //Z=13595
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13596
                            if ( n<n5f ) n5f = n;  //Z=13597
                        }/*7*/  //Z=13598
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13599
                            if ( n<n6f ) n6f = n;  //Z=13600
                        }/*7*/  //Z=13601
                    }/*6*/  /*  of n-loop  */  //Z=13602
                }/*5*/  /*  cs=1  */  //Z=13603

                /*  inhomogeneous core/shell  */  //Z=13605
                if ( params.cs==2 )
                {/*5*/  //Z=13606
                    //i = 2;  //Z=13607
                    double xrmn_n = 1.0;
                    double z12v_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13608
                        z12v_n = z12v_n*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13609
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13610
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=13611
                        fkv[n] = fkv[n-1]*n;  //Z=13612
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13613
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13614
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13615 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13616
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13617
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13618
                        pn[n] = pn[n-1]*p*p;  //Z=13619
                        /*  longitudinal  */  //Z=13620
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13621
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13622
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=13623 */
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13624
                            /* carr1fm[i]:=carr1pm[i];  //Z=13625 */
                            /* i:=i+1;  //Z=13626 */
                        }/*7*/  //Z=13627
                        /*  P(q)-coefficient  */  //Z=13628
                        params.CR->carr1p[n] = pow(4.0,2*n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=13629
                        /*  F(q)-coefficient  */  //Z=13630
                        double sump = 0.0;  //Z=13631
                        for ( int m=0; m<=n; m++ ) sump += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13632
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=13633

                        /*  cross-sectional P(q)  */  //Z=13635
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v_n*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13636
                        sump = 0.0;  //Z=13637
                        sump1 = 0.0;  //Z=13638
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13639
                            const double sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=13640
                            sump += pn[n-m]*sumi;  //Z=13641
                            sump1 += sumi;  //Z=13642
                        }/*7*/  //Z=13643
                        params.CR->carr5p[n] = (1-params.uca/2.0)*z12v_n*xrmn_n*sump;  //Z=13644
                        params.CR->carr6p[n] = (1-params.uca/2.0)*z12v_n*xrn[n]*sump1;  //Z=13645
                        sump = 0.0;  //Z=13646
                        sump1 = 0.0;  //Z=13647
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13648
                            const double sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=13649
                            sump += sumi;  //Z=13650
                            sump1 += pn[n-m]*sumi;  //Z=13651
                        }/*7*/  //Z=13652
                        params.CR->carr7p[n] = sqr(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump;  //Z=13653
                        params.CR->carr8p[n] = sqr(1-params.alphash1/2.0)*z12v_n*xrmn_n*sump1;  //Z=13654
                        params.CR->carr9p[n] = sqr(1-params.alphash1/2.0)*z12v_n*xrn[n]*sump;  //Z=13655

                        /* (* P(q)-coefficients *)  //Z=13657
                        (* cross-sectional *)  //Z=13658
                        (* F121 *)  //Z=13659
                        carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=13660
                        (* F122 *)  //Z=13661
                        sump:=0.0;  //Z=13662
                           for m:=0 to n do begin  //Z=13663
                           sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13664
                        end;  //Z=13665
                        carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=13666
                        (* F123 *)  //Z=13667
                        carr6p[n]:=carr4p[n]/pn[n];  */  //Z=13668

                        /*  F(q)-coefficients  */  //Z=13670
                        /*  cross-sectional  */  //Z=13671
                        params.CR->carr4f[n] = z12v_n*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=13672
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v_n*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=13673
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v_n*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=13674

                        if ( search1 )
                        {/*7*/  //Z=13676
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13677
                                n1 = n;  //Z=13678
                                search1 = false;  //Z=13679
                            }/*8*/  //Z=13680
                        }/*7*/  //Z=13681
                        if ( search4 )
                        {/*7*/  //Z=13682
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13683
                                n4 = n;  //Z=13684
                                search4 = false;  //Z=13685
                            }/*8*/  //Z=13686
                        }/*7*/  //Z=13687
                        if ( search5 )
                        {/*7*/  //Z=13688
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=13689
                                n5 = n;  //Z=13690
                                search5 = false;  //Z=13691
                            }/*8*/  //Z=13692
                        }/*7*/  //Z=13693
                        if ( search6 )
                        {/*7*/  //Z=13694
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=13695
                                n6 = n;  //Z=13696
                                search6 = false;  //Z=13697
                            }/*8*/  //Z=13698
                        }/*7*/  //Z=13699
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=13700
                            if ( n<n7 ) n7 = n;  //Z=13701
                        }/*7*/  //Z=13702
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=13703
                            if ( n<n8 ) n8 = n;  //Z=13704
                        }/*7*/  //Z=13705
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=13706
                            if ( n<n9 ) n9 = n;  //Z=13707
                        }/*7*/  //Z=13708
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13709
                            if ( n<n1f ) n1f = n;  //Z=13710
                        }/*7*/  //Z=13711
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13712
                            if ( n<n4f ) n4f = n;  //Z=13713
                        }/*7*/  //Z=13714
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13715
                            if ( n<n5f ) n5f = n;  //Z=13716
                        }/*7*/  //Z=13717
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13718
                            if ( n<n6f ) n6f = n;  //Z=13719
                        }/*7*/  //Z=13720
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=13721
                            if ( n<n7f ) n7f = n;  //Z=13722
                        }/*7*/  //Z=13723
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=13724
                            if ( n<n8f ) n8f = n;  //Z=13725
                        }/*7*/  //Z=13726
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=13727
                            if ( n<n9f ) n9f = n;  //Z=13728
                        }/*7*/  //Z=13729
                    }/*6*/  /*  of n-loop  */  //Z=13730
                }/*5*/  /*  cs=2  */  //Z=13731


                /*  myelin  */  //Z=13734
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=13735
                    //i = 2;  //Z=13736
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13737
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13738
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13739
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=13740
                        fkv[n] = fkv[n-1]*n;  //Z=13741
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13742
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13743
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13744 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13745
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=13746
                        /*  longitudinal  */  //Z=13747
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13748
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13749
                            /* carr1pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=13750 */
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13751
                            /* carr1fm[i]:=carr1pm[i];  //Z=13752 */
                            /* i:=i+1;  //Z=13753 */
                        }/*7*/  //Z=13754
                        /*  P(q)-coefficient  */  //Z=13755
                        params.CR->carr1p[n] = pow(4.0,2*n)*z12vl[n]*xln[n]/((2.0*n+1)*(n+1));  //Z=13756
                        /*  F(q)-coefficient  */  //Z=13757
                        double sump = 0.0;  //Z=13758
                        for ( int m=0; m<=n; m++ ) sump += z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13759
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump;  //Z=13760

                        /*  cross-section  */  //Z=13762
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=13763

                        sump = 0.0;  //Z=13765
                        for ( int m=0; m<=n; m++ ) sump += z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13766
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=13767

                        if ( search1 )
                        {/*7*/  //Z=13769
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13770
                                n1 = n;  //Z=13771
                                search1 = false;  //Z=13772
                            }/*8*/  //Z=13773
                        }/*7*/  //Z=13774
                        if ( search4 )
                        {/*7*/  //Z=13775
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13776
                                n4 = n;  //Z=13777
                                search4 = false;  //Z=13778
                            }/*8*/  //Z=13779
                        }/*7*/  //Z=13780
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13781
                            if ( n<n1f ) n1f = n;  //Z=13782
                        }/*7*/  //Z=13783
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13784
                            if ( n<n4f ) n4f = n;  //Z=13785
                        }/*7*/  //Z=13786
                    }/*6*/  /*  of n-loop  */  //Z=13787
                }/*5*/  /*  cs=3  */  //Z=13788

            }/*4*/  /*  of cylinders  */  //Z=13790


            /* ** disks ** */  //Z=13793
            if ( dim==2 )
            {/*4*/  //Z=13794
                double z12v[nnmax+1];   // local
                z12v[0] = 1.0;
                /*  homogeneous  */  //Z=13795
                if ( params.cs==0 )
                {/*5*/  //Z=13796
                    //i = 2;  //Z=13797
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13798
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13799
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13800
                        fkv[n] = fkv[n-1]*n;  //Z=13801
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13802
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13803
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13804 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13805
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13806
                        /*  longitudinal  */  //Z=13807
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13808
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13809
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=13810 */
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13811
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=13812 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=13813
                            /* carr2fm[i]:=carr2pm[i];  //Z=13814 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13815 */
                            /* i:=i+1;  //Z=13816 */
                        }/*7*/  //Z=13817
                        /*  P(q)-coefficient  */  //Z=13818
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13819
                        /*  F(q)-coefficient  */  //Z=13820
                        double sump = 0.0;  //Z=13821
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13822
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13823
                        }/*7*/  //Z=13824
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13825

                        /*  cross-section  */  //Z=13827
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13828
                        sump = 0.0;  //Z=13829
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13830
                            sump += z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13831
                        }/*7*/  //Z=13832
                        params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=13833

                        /*  series for <...> integration  */  //Z=13835
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13836
                        sump = 0.0;  //Z=13837
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13838
                            sump += z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13839
                        }/*7*/  //Z=13840
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13841

                        if ( search1 )
                        {/*7*/  //Z=13843
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13844
                                n1 = n;  //Z=13845
                                search1 = false;  //Z=13846
                            }/*8*/  //Z=13847
                        }/*7*/  //Z=13848
                        if ( search4 )
                        {/*7*/  //Z=13849
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13850
                                n4 = n;  //Z=13851
                                search4 = false;  //Z=13852
                            }/*8*/  //Z=13853
                        }/*7*/  //Z=13854
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13855
                            if ( n<n1f ) n1f = n;  //Z=13856
                        }/*7*/  //Z=13857
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13858
                            if ( n<n4f ) n4f = n;  //Z=13859
                        }/*7*/  //Z=13860
                    }/*6*/  /*  of n-loop  */  //Z=13861
                }/*5*/  /*  of cs=0  */  //Z=13862

                /*  core/shell  */  //Z=13864
                if ( params.cs==1 )
                {/*5*/  //Z=13865
                    //i = 2;  //Z=13866
                    double xrmn_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13867
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13868
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13869
                        fkv[n] = fkv[n-1]*n;  //Z=13870
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13871
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13872
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13873 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13874
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13875
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13876
                        pn[n] = pn[n-1]*p*p;  //Z=13877
                        /*  longitudinal  */  //Z=13878
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13879
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13880
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=13881 */
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13882
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=13883 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=13884
                            /* carr2fm[i]:=carr2pm[i];  //Z=13885 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13886 */
                            /* i:=i+1;  //Z=13887 */
                        }/*7*/  //Z=13888
                        /*  P(q)-coefficient  */  //Z=13889
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=13890
                        /*  F(q)-coefficient  */  //Z=13891
                        sump = 0.0;  //Z=13892
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13893
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=13894
                        }/*7*/  //Z=13895
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=13896

                        /*  cross-sectional  */  //Z=13898
                        /*  F121  */  //Z=13899
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=13900
                        /*  F122  */  //Z=13901
                        sump = 0.0;  //Z=13902
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13903
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13904
                        }/*7*/  //Z=13905
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=13906
                        /*  F122  */  //Z=13907
                        sump = 0.0;  //Z=13908
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13909
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13910
                        }/*7*/  //Z=13911
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=13912
                        /*  F123  */  //Z=13913
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=13914
                        /*  F121  */  //Z=13915
                        sump = 0.0;  //Z=13916
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13917
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13918
                        }/*7*/  //Z=13919
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=13920
                        /*  F122  */  //Z=13921
                        sump = 0.0;  //Z=13922
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13923
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=13924
                        }/*7*/  //Z=13925
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=13926
                        /*  F123  */  //Z=13927
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=13928

                        /*  series for <...> integration  */  //Z=13930
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=13931
                        sump = 0.0;  //Z=13932
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13933
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=13934
                        }/*7*/  //Z=13935
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=13936

                        if ( search1 )
                        {/*7*/  //Z=13938
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=13939
                                n1 = n;  //Z=13940
                                search1 = false;  //Z=13941
                            }/*8*/  //Z=13942
                        }/*7*/  //Z=13943
                        if ( search4 )
                        {/*7*/  //Z=13944
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=13945
                                n4 = n;  //Z=13946
                                search4 = false;  //Z=13947
                            }/*8*/  //Z=13948
                        }/*7*/  //Z=13949
                        if ( search5 )
                        {/*7*/  //Z=13950
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=13951
                                n5 = n;  //Z=13952
                                search5 = false;  //Z=13953
                            }/*8*/  //Z=13954
                        }/*7*/  //Z=13955
                        if ( search6 )
                        {/*7*/  //Z=13956
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=13957
                                n6 = n;  //Z=13958
                                search6 = false;  //Z=13959
                            }/*8*/  //Z=13960
                        }/*7*/  //Z=13961
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=13962
                            if ( n<n1f ) n1f = n;  //Z=13963
                        }/*7*/  //Z=13964
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=13965
                            if ( n<n4f ) n4f = n;  //Z=13966
                        }/*7*/  //Z=13967
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=13968
                            if ( n<n5f ) n5f = n;  //Z=13969
                        }/*7*/  //Z=13970
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=13971
                            if ( n<n6f ) n6f = n;  //Z=13972
                        }/*7*/  //Z=13973
                    }/*6*/  /*  of n-loop  */  //Z=13974
                }/*5*/  /*  of cs=1  */  //Z=13975

                /*  inhomogeneous core/shell  */  //Z=13977
                if ( params.cs==2 )
                {/*5*/  //Z=13978
                    //i = 2;  //Z=13979
                    double xrmn_n = 1.0;
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=13980
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=13981
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=13982
                        fkv[n] = fkv[n-1]*n;  //Z=13983
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=13984
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=13985
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=13986 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=13987
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=13988
                        xrmn_n = -xrmn_n*xrm2z;  //Z=13989
                        pn[n] = pn[n-1]*p*p;  //Z=13990
                        /*  longitudinal  */  //Z=13991
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=13992
                            qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);
                            //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,3,0,2*n-2*m,2*m,params.CR->carr1p,intl);  //Z=13993
                            /* carr2pm[i]:=intl/(power(4,m)*fk2v[n-m]*fkv[m]*fkv[m]*norm);  //Z=13994 */
                            params.CR->carr11pm[n][m] = intl/(pow(4.0,m)*fk2v[n-m]*fkv[m]*fkv[m]*params.norm);  //Z=13995
                            /* carr1pm[i]:=power(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=13996 */
                            params.CR->carr11pm[n][m] = pow(-1,m)*fk2v[m]/(fkv[m]*fkv[n-m]);  //Z=13997
                            /* carr2fm[i]:=carr2pm[i];  //Z=13998 */
                            /* carr1fm[i]:=carr1pm[i];  //Z=13999 */
                            /* i:=i+1;  //Z=14000 */
                        }/*7*/  //Z=14001
                        /*  P(q)-coefficient  */  //Z=14002
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14003
                        /*  F(q)-coefficient  */  //Z=14004
                        sump = 0.0;  //Z=14005
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14006
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14007
                        }/*7*/  //Z=14008
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump;  //Z=14009

                        /*  cross-sectional P(q)  */  //Z=14011
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14012
                        sump = 0.0;  //Z=14013
                        sump1 = 0.0;  //Z=14014
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14015
                            const double sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=14016
                            sump = sump+pn[n-m]*sumi;  //Z=14017
                            sump1 = sump1+sumi;  //Z=14018
                        }/*7*/  //Z=14019
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14020
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=14021
                        sump = 0.0;  //Z=14022
                        sump1 = 0.0;  //Z=14023
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14024
                            const double sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=14025
                            sump = sump+sumi;  //Z=14026
                            sump1 = sump1+pn[n-m]*sumi;  //Z=14027
                        }/*7*/  //Z=14028
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14029
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=14030
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=14031

                        /*  (* F121 *)  //Z=14033
                         carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=14034
                         (* F122 *)  //Z=14035
                         sump:=0.0;  //Z=14036
                            for m:=0 to n do begin  //Z=14037
                            sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14038
                         end;  //Z=14039
                         carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=14040
                         (* F122 *)  //Z=14041
                         sump:=0.0;  //Z=14042
                            for m:=0 to n do begin  //Z=14043
                            sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14044
                         end;  //Z=14045
                         carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=14046
                         (* F123 *)  //Z=14047
                         carr6p[n]:=carr4p[n]/pn[n];       */  //Z=14048

                        /*  F(q)  */  //Z=14050
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=14051
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14052
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14053

                        /*  series for <...> integration  */  //Z=14055
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14056
                        sump = 0.0;  //Z=14057
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14058
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14059
                        }/*7*/  //Z=14060
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14061

                        if ( search1 )
                        {/*7*/  //Z=14063
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14064
                                n1 = n;  //Z=14065
                                search1 = false;  //Z=14066
                            }/*8*/  //Z=14067
                        }/*7*/  //Z=14068
                        if ( search4 )
                        {/*7*/  //Z=14069
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14070
                                n4 = n;  //Z=14071
                                search4 = false;  //Z=14072
                            }/*8*/  //Z=14073
                        }/*7*/  //Z=14074
                        if ( search5 )
                        {/*7*/  //Z=14075
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=14076
                                n5 = n;  //Z=14077
                                search5 = false;  //Z=14078
                            }/*8*/  //Z=14079
                        }/*7*/  //Z=14080
                        if ( search6 )
                        {/*7*/  //Z=14081
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=14082
                                n6 = n;  //Z=14083
                                search6 = false;  //Z=14084
                            }/*8*/  //Z=14085
                        }/*7*/  //Z=14086
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=14087
                            if ( n<n7 ) n7 = n;  //Z=14088
                        }/*7*/  //Z=14089
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=14090
                            if ( n<n8 ) n8 = n;  //Z=14091
                        }/*7*/  //Z=14092
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=14093
                            if ( n<n9 ) n9 = n;  //Z=14094
                        }/*7*/  //Z=14095
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14096
                            if ( n<n1f ) n1f = n;  //Z=14097
                        }/*7*/  //Z=14098
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14099
                            if ( n<n4f ) n4f = n;  //Z=14100
                        }/*7*/  //Z=14101
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14102
                            if ( n<n5f ) n5f = n;  //Z=14103
                        }/*7*/  //Z=14104
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14105
                            if ( n<n6f ) n6f = n;  //Z=14106
                        }/*7*/  //Z=14107
                    }/*6*/  /*  of n-loop  */  //Z=14108
                }/*5*/  /*  of cs=2  */  //Z=14109

            }/*4*/  /*  of disks  */  //Z=14111
        }/*3*/  //Z=14112

        /* ** z-axis ** */  //Z=14114
        if ( (params.orcase==4) && (dim!=3) )
        {/*3*/  //Z=14115
            double z12v[nnmax+1];   // local
            z12v[0] = 1.0;
            /* ** cylinders ** */  //Z=14116
            if ( dim==1 )
            {/*4*/  //Z=14117
                /*  homogeneous  */  //Z=14118
                if ( params.cs==0 )
                {/*5*/  //Z=14119
                    //i = 2;  //Z=14120
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14121
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14122
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14123
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=14124
                        fkv[n] = fkv[n-1]*n;  //Z=14125
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14126
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14127
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14128 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14129
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14130
                        /*  longitudinal  */  //Z=14131
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);
                        //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14132
                        /*  P(q)-coefficient  */  //Z=14133
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14134
                        /*  F(q)-coefficient  */  //Z=14135
                        sump = 0.0;  //Z=14136
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14137
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n)*fkv[n]*fkv[n]*params.norm);  //Z=14138
                        /*  cross-sectional  */  //Z=14139
                        /*  P(q)-coefficient  */  //Z=14140
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]);  //Z=14141
                        /*  F(q)-coefficient  */  //Z=14142
                        sump = 0.0;  //Z=14143
                        for ( int m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14144
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=14145

                        if ( search1 )
                        {/*7*/  //Z=14147
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14148
                                n1 = n;  //Z=14149
                                search1 = false;  //Z=14150
                            }/*8*/  //Z=14151
                        }/*7*/  //Z=14152
                        if ( search4 )
                        {/*7*/  //Z=14153
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14154
                                n4 = n;  //Z=14155
                                search4 = false;  //Z=14156
                            }/*8*/  //Z=14157
                        }/*7*/  //Z=14158
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14159
                            if ( n<n1f ) n1f = n;  //Z=14160
                        }/*7*/  //Z=14161
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14162
                            if ( n<n4f ) n4f = n;  //Z=14163
                        }/*7*/  //Z=14164
                    }/*6*/  /*  of n-loop  */  //Z=14165
                }/*5*/  /*  of cs=0  */  //Z=14166

                /*  core/shell  */  //Z=14168
                if ( params.cs==1 )
                {/*5*/  //Z=14169
                    //i = 2;  //Z=14170
                    double xrmn_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14171
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14172
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14173
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=14174
                        fkv[n] = fkv[n-1]*n;  //Z=14175
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14176
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14177
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14178 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14179
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14180
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14181
                        pn[n] = pn[n-1]*p*p;  //Z=14182
                        /*  longitudinal  */  //Z=14183
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);
                        //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14184
                        /*  P(q)-coefficient  */  //Z=14185
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14186
                        /*  F(q)-coefficient  */  //Z=14187
                        sump = 0.0;  //Z=14188
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14189
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14190
                        /*  P(q)-coefficients  */  //Z=14191
                        /*  cross-sectional  */  //Z=14192
                        /*  F121  */  //Z=14193
                        params.CR->carr4p[n] = 4*(n+1/2.0)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14194
                        /*  F122  */  //Z=14195
                        sump = 0.0;  //Z=14196
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14197
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14198
                        }/*7*/  //Z=14199
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=14200
                        /*  F123  */  //Z=14201
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=14202

                        /*  F(q)-coefficients  */  //Z=14204
                        /*  cross-sectional  */  //Z=14205
                        /*  F121  */  //Z=14206
                        sump = 0.0;  //Z=14207
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14208
                            sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14209
                        }/*7*/  //Z=14210
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=14211
                        /*  F122  */  //Z=14212
                        sump = 0.0;  //Z=14213
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14214
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14215
                        }/*7*/  //Z=14216
                        params.CR->carr5f[n] = xrmn_n*sump;  //Z=14217
                        /*  F123  */  //Z=14218
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=14219

                        if ( search1 )
                        {/*7*/  //Z=14221
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14222
                                n1 = n;  //Z=14223
                                search1 = false;  //Z=14224
                            }/*8*/  //Z=14225
                        }/*7*/  //Z=14226
                        if ( search4 )
                        {/*7*/  //Z=14227
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14228
                                n4 = n;  //Z=14229
                                search4 = false;  //Z=14230
                            }/*8*/  //Z=14231
                        }/*7*/  //Z=14232
                        if ( search5 )
                        {/*7*/  //Z=14233
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=14234
                                n5 = n;  //Z=14235
                                search5 = false;  //Z=14236
                            }/*8*/  //Z=14237
                        }/*7*/  //Z=14238
                        if ( search6 )
                        {/*7*/  //Z=14239
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=14240
                                n6 = n;  //Z=14241
                                search6 = false;  //Z=14242
                            }/*8*/  //Z=14243
                        }/*7*/  //Z=14244
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14245
                            if ( n<n1f ) n1f = n;  //Z=14246
                        }/*7*/  //Z=14247
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14248
                            if ( n<n4f ) n4f = n;  //Z=14249
                        }/*7*/  //Z=14250
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14251
                            if ( n<n5f ) n5f = n;  //Z=14252
                        }/*7*/  //Z=14253
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14254
                            if ( n<n6f ) n6f = n;  //Z=14255
                        }/*7*/  //Z=14256
                    }/*6*/  /*  of n-loop  */  //Z=14257
                }/*5*/  /*  of cs=1  */  //Z=14258


                /*  inhomogeneous core/shell  */  //Z=14261
                if ( params.cs==2 )
                {/*5*/  //Z=14262
                    //i = 2;  //Z=14263
                    double xrmn_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14264
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14265
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14266
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=14267
                        fkv[n] = fkv[n-1]*n;  //Z=14268
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14269
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14270
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14271 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14272
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14273
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14274
                        pn[n] = pn[n-1]*p*p;  //Z=14275
                        /*  longitudinal  */  //Z=14276
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);
                        //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14277
                        /*  P(q)-coefficient  */  //Z=14278
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14279
                        /*  F(q)-coefficient  */  //Z=14280
                        sump = 0.0;  //Z=14281
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14282
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14283

                        /*  cross-sectional P(q)  */  //Z=14285
                        params.CR->carr4p[n] = pow(4.0,n+1)*gam3[n]*z12v[n]*xrn[n]/(sqrt(M_PI)*(n+2)*(n+1)*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14286
                        sump = 0.0;  //Z=14287
                        sump1 = 0.0;  //Z=14288
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14289
                            const double sumi = 1/((m+1-params.alphash1/2.0)*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]*fkv[m]);  //Z=14290
                            sump = sump+pn[n-m]*sumi;  //Z=14291
                            sump1 = sump1+sumi;  //Z=14292
                        }/*7*/  //Z=14293
                        params.CR->carr5p[n] = (1-params.uca/2.0)*z12v[n]*xrmn_n*sump;  //Z=14294
                        params.CR->carr6p[n] = (1-params.uca/2.0)*z12v[n]*xrn[n]*sump1;  //Z=14295
                        sump = 0.0;  //Z=14296
                        sump1 = 0.0;  //Z=14297
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14298
                            const double sumi = 1/((n-m+1-params.alphash1/2.0)*(m+1-params.alphash1/2.0)*fkv[n-m]*fkv[m]*fkv[m]*fkv[n-m]);  //Z=14299
                            sump = sump+sumi;  //Z=14300
                            sump1 = sump1+pn[n-m]*sumi;  //Z=14301
                        }/*7*/  //Z=14302
                        params.CR->carr7p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump;  //Z=14303
                        params.CR->carr8p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrmn_n*sump1;  //Z=14304
                        params.CR->carr9p[n] = (1-params.alphash1/2.0)*(1-params.alphash1/2.0)*z12v[n]*xrn[n]*sump;  //Z=14305

                        /*  (* P(q)-coefficients *)  //Z=14307
                        (* cross-sectional *)  //Z=14308
                        (* F121 *)  //Z=14309
                        carr4p[n]:=4*(n+1/2)*fk2v[n]*z12v[n]*xrn[n]/((n+2)*(n+1)*fkv[n]*fkv[n]*(n+1)*fkv[n]*fkv[n]);  //Z=14310
                        (* F122 *)  //Z=14311
                        sump:=0.0;  //Z=14312
                           for m:=0 to n do begin  //Z=14313
                           sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14314
                        end;  //Z=14315
                        carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=14316
                        (* F123 *)  //Z=14317
                        carr6p[n]:=carr4p[n]/pn[n];    */  //Z=14318

                        /*  F(q)-coefficients  */  //Z=14320
                        /*  cross-sectional  */  //Z=14321
                        params.CR->carr4f[n] = z12v[n]*xrn[n]/((n+1)*fkv[n]*fkv[n]);  //Z=14322
                        params.CR->carr5f[n] = (1-params.alphash1/2.0)*z12v[n]*xrmn_n/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14323
                        params.CR->carr6f[n] = (1-params.alphash1/2.0)*z12v[n]*xrn[n]/((n+1-params.alphash1/2.0)*fkv[n]*fkv[n]);  //Z=14324

                        if ( search1 )
                        {/*7*/  //Z=14326
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14327
                                n1 = n;  //Z=14328
                                search1 = false;  //Z=14329
                            }/*8*/  //Z=14330
                        }/*7*/  //Z=14331
                        if ( search4 )
                        {/*7*/  //Z=14332
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14333
                                n4 = n;  //Z=14334
                                search4 = false;  //Z=14335
                            }/*8*/  //Z=14336
                        }/*7*/  //Z=14337
                        if ( search5 )
                        {/*7*/  //Z=14338
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=14339
                                n5 = n;  //Z=14340
                                search5 = false;  //Z=14341
                            }/*8*/  //Z=14342
                        }/*7*/  //Z=14343
                        if ( search6 )
                        {/*7*/  //Z=14344
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=14345
                                n6 = n;  //Z=14346
                                search6 = false;  //Z=14347
                            }/*8*/  //Z=14348
                        }/*7*/  //Z=14349
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=14350
                            if ( n<n7 ) n7 = n;  //Z=14351
                        }/*7*/  //Z=14352
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=14353
                            if ( n<n8 ) n8 = n;  //Z=14354
                        }/*7*/  //Z=14355
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=14356
                            if ( n<n9 ) n9 = n;  //Z=14357
                        }/*7*/  //Z=14358
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14359
                            if ( n<n1f ) n1f = n;  //Z=14360
                        }/*7*/  //Z=14361
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14362
                            if ( n<n4f ) n4f = n;  //Z=14363
                        }/*7*/  //Z=14364
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14365
                            if ( n<n5f ) n5f = n;  //Z=14366
                        }/*7*/  //Z=14367
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14368
                            if ( n<n6f ) n6f = n;  //Z=14369
                        }/*7*/  //Z=14370
                    }/*6*/  /*  of n-loop  */  //Z=14371
                }/*5*/  /*  of cs=2  */  //Z=14372

                /*  myelin  */  //Z=14374
                if ( (params.cs==3) || (params.cs==4) )
                {/*5*/  //Z=14375
                    //i = 2;  //Z=14376
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14377
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14378
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14379
                        //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=14380
                        fkv[n] = fkv[n-1]*n;  //Z=14381
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14382
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14383
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14384 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14385
                        xrn[n] = -xrn[n-1]*x12zm;  //Z=14386
                        /*  longitudinal  */  //Z=14387
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);
                        //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14388
                        /*  P(q)-coefficient  */  //Z=14389
                        params.CR->carr1p[n] = pow(4.0,n)*z12vl[n]*xln[n]*intl/((2.0*n+1)*(n+1)*fkv[n]*fkv[n]*params.norm);  //Z=14390
                        /*  F(q)-coefficient  */  //Z=14391
                        sump = 0.0;  //Z=14392
                        for ( int m=0; m<=n; m++ ) sump = sump+z12vl[m]*z12vl[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14393
                        params.CR->carr1f[n] = fk2v[n]*xln[n]*sump*intl/(pow(4.0,n)*fkv[n]*fkv[n]*params.norm);  //Z=14394
                        /*  cross-sectional  */  //Z=14395
                        /*  P(q)-coefficient  */  //Z=14396
                        params.CR->carr4p[n] = z12v[n]*xrn[n];  //Z=14397

                        /*  F(q)-coefficient  */  //Z=14399
                        sump = 0.0;  //Z=14400
                        for ( int m=0; m<=n; m++ ) sump = sump+z12v[m]*z12v[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14401
                        params.CR->carr4f[n] = xrn[n]*sump;  //Z=14402

                        if ( search1 )
                        {/*7*/  //Z=14404
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14405
                                n1 = n;  //Z=14406
                                search1 = false;  //Z=14407
                            }/*8*/  //Z=14408
                        }/*7*/  //Z=14409
                        if ( search4 )
                        {/*7*/  //Z=14410
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14411
                                n4 = n;  //Z=14412
                                search4 = false;  //Z=14413
                            }/*8*/  //Z=14414
                        }/*7*/  //Z=14415
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14416
                            if ( n<n1f ) n1f = n;  //Z=14417
                        }/*7*/  //Z=14418
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14419
                            if ( n<n4f ) n4f = n;  //Z=14420
                        }/*7*/  //Z=14421
                    }/*6*/  /*  of n-loop  */  //Z=14422
                }/*5*/  /*  of cs=3  */  //Z=14423
            }/*4*/  /*  of cylinders  */  //Z=14424

            /* ** disks ** */  //Z=14426
            if ( dim==2 )
            {/*4*/  //Z=14427
                /*  homogeneous  */  //Z=14428
                if ( params.cs==0 )
                {/*5*/  //Z=14429
                    carr2i[0] = 1;  //Z=14430
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14431
                        fkv[n] = fkv[n-1]*n;  //Z=14432
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14433
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);
                        //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14434
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14435
                    }/*6*/  //Z=14436
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14437
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14438
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14439
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14440
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14441 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14442
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14443
                        /*  longitudinal  */  //Z=14444
                        sump = 0.0;  //Z=14445
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14446
                            sump = sump+carr2i[m]/fkv[n-m];  //Z=14447
                        }/*7*/  //Z=14448
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14449
                        sump1 = 0.0;  //Z=14450
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14451
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14452
                        }/*7*/  //Z=14453
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;  //Z=14454

                        /*  cross-sectional  */  //Z=14456
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14457
                        sump = 0.0;  //Z=14458
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14459
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14460
                        }/*7*/  //Z=14461
                        params.CR->carr4f[n] = M_PI*xln[n]*sump/4.0;  //Z=14462

                        /*  series for <...> integration  */  //Z=14464
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14465
                        sump = 0.0;  //Z=14466
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14467
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14468
                        }/*7*/  //Z=14469
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14470

                        if ( search1 )
                        {/*7*/  //Z=14472
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14473
                                n1 = n;  //Z=14474
                                search1 = false;  //Z=14475
                            }/*8*/  //Z=14476
                        }/*7*/  //Z=14477
                        if ( search4 )
                        {/*7*/  //Z=14478
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14479
                                n4 = n;  //Z=14480
                                search4 = false;  //Z=14481
                            }/*8*/  //Z=14482
                        }/*7*/  //Z=14483
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14484
                            if ( n<n1f ) n1f = n;  //Z=14485
                        }/*7*/  //Z=14486
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14487
                            if ( n<n4f ) n4f = n;  //Z=14488
                        }/*7*/  //Z=14489
                    }/*6*/  /*  of n-loop  */  //Z=14490
                }/*5*/  /*  of cs=0  */  //Z=14491


                /*  core/shell  */  //Z=14494
                if ( params.cs==1 )
                {/*5*/  //Z=14495
                    carr2i[0] = 1;  //Z=14496
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14497
                        fkv[n] = fkv[n-1]*n;  //Z=14498
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14499
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);
                        //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14500
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14501
                    }/*6*/  //Z=14502
                    double xrmn_n = 1.0;
                    for ( n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14503
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14504
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14505
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14506
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14507 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14508
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14509
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14510
                        pn[n] = pn[n-1]*p*p;  //Z=14511
                        /*  longitudinal  */  //Z=14512
                        sump = 0.0;  //Z=14513
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14514
                            sump = sump+carr2i[m]/fkv[n-m];  //Z=14515
                        }/*7*/  //Z=14516
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14517
                        sump1 = 0.0;  //Z=14518
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14519
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14520
                        }/*7*/  //Z=14521
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;  //Z=14522

                        /*  cross-sectional  */  //Z=14524
                        /*  F121  */  //Z=14525
                        params.CR->carr4p[n] = pow(4.0,n)*sqrt(M_PI)*z12v[n]*xrn[n]/(2.0*(n+1)*gam3[n]*fkv[n]);  //Z=14526
                        /*  F122  */  //Z=14527
                        sump = 0.0;  //Z=14528
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14529
                            sump = sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14530
                        }/*7*/  //Z=14531
                        params.CR->carr5p[n] = z12v[n]*xrmn_n*sump;  //Z=14532
                        /*  F122  */  //Z=14533
                        sump = 0.0;  //Z=14534
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14535
                            sump = sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14536
                        }/*7*/  //Z=14537
                        params.CR->carr5p[n] = M_PI*z12v[n]*xrmn_n*sump/4.0;  //Z=14538
                        /*  F123  */  //Z=14539
                        params.CR->carr6p[n] = params.CR->carr4p[n]/pn[n];  //Z=14540
                        /*  F121  */  //Z=14541
                        sump = 0.0;  //Z=14542
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14543
                            sump = sump+z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14544
                        }/*7*/  //Z=14545
                        params.CR->carr4f[n] = M_PI*xrn[n]*sump/4.0;  //Z=14546
                        /*  F122  */  //Z=14547
                        sump = 0.0;  //Z=14548
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14549
                            sump = sump+pn[m]*z12v[m]*z12v[n-m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14550
                        }/*7*/  //Z=14551
                        params.CR->carr5f[n] = M_PI*xrmn_n*sump/4.0;  //Z=14552
                        /*  F123  */  //Z=14553
                        params.CR->carr6f[n] = params.CR->carr4f[n]/pn[n];  //Z=14554

                        /*  series for <...> integration  */  //Z=14556
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14557
                        sump = 0.0;  //Z=14558
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14559
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14560
                        }/*7*/  //Z=14561
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14562

                        if ( search1 )
                        {/*7*/  //Z=14564
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14565
                                n1 = n;  //Z=14566
                                search1 = false;  //Z=14567
                            }/*8*/  //Z=14568
                        }/*7*/  //Z=14569
                        if ( search4 )
                        {/*7*/  //Z=14570
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14571
                                n4 = n;  //Z=14572
                                search4 = false;  //Z=14573
                            }/*8*/  //Z=14574
                        }/*7*/  //Z=14575
                        if ( search5 )
                        {/*7*/  //Z=14576
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=14577
                                n5 = n;  //Z=14578
                                search5 = false;  //Z=14579
                            }/*8*/  //Z=14580
                        }/*7*/  //Z=14581
                        if ( search6 )
                        {/*7*/  //Z=14582
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=14583
                                n6 = n;  //Z=14584
                                search6 = false;  //Z=14585
                            }/*8*/  //Z=14586
                        }/*7*/  //Z=14587
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14588
                            if ( n<n1f ) n1f = n;  //Z=14589
                        }/*7*/  //Z=14590
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14591
                            if ( n<n4f ) n4f = n;  //Z=14592
                        }/*7*/  //Z=14593
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14594
                            if ( n<n5f ) n5f = n;  //Z=14595
                        }/*7*/  //Z=14596
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14597
                            if ( n<n6f ) n6f = n;  //Z=14598
                        }/*7*/  //Z=14599
                    }/*6*/  /*  of n-loop  */  //Z=14600
                }/*5*/  /*  of cs=1  */  //Z=14601

                /*  inhomogeneous core/shell  */  //Z=14603
                if ( params.cs==2 )
                {/*5*/  //Z=14604
                    carr2i[0] = 1;  //Z=14605
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14606
                        fkv[n] = fkv[n-1]*n;  //Z=14607
                        fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);  //Z=14608
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,1,1,1,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);
                        //qrombdeltac(params.length,r,p1,sigma,alfa,dbeta,theta,phi,1,1,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,dim,1,4,0,0,2*n,params.CR->carr1p,intl);  //Z=14609
                        carr2i[n] = pow(-1,n)*fk2v[n]*intl/(pow(4.0,n)*fkv[n]*fkv[n]*fkv[n]*params.norm);  //Z=14610
                    }/*6*/  //Z=14611
                    double xrmn_n = 1.0;
                    for ( int n=1; n<=nmax; n++ )
                    {/*6*/  //Z=14612
                        z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);  //Z=14613
                        z12vl[n] = z12vl[n-1]*((zl+1)-2+2*n)*((zl+1)-1+2*n);  //Z=14614
                        gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=14615
                        /* e1[n]:=e1[n-1]*(epsi*epsi-1);  //Z=14616 */
                        xln[n] = -xln[n-1]*xl2z;  //Z=14617
                        xrn[n] = -xrn[n-1]*xr2z;  //Z=14618
                        xrmn_n = -xrmn_n*xrm2z;  //Z=14619
                        pn[n] = pn[n-1]*p*p;  //Z=14620
                        /*  longitudinal  */  //Z=14621
                        sump = 0.0;  //Z=14622
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14623
                            sump = sump+carr2i[m]/fkv[n-m];  //Z=14624
                        }/*7*/  //Z=14625
                        params.CR->carr1p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]*sump/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]);  //Z=14626
                        sump1 = 0.0;  //Z=14627
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14628
                            sump1 = sump1+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]);  //Z=14629
                        }/*7*/  //Z=14630
                        params.CR->carr1f[n] = xln[n]*fkv[n]*sump1*sump;  //Z=14631

                        /*  cross-sectional P(q)  */  //Z=14633
                        params.CR->carr4p[n] = (sqrt(M_PI)/2.0)*pow(4.0,n)*z12v[n]*xrn[n]/((n+1)*gam3[n]*fkv[n]);  //Z=14634
                        sump = 0.0;  //Z=14635
                        sump1 = 0.0;  //Z=14636
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14637
                            const double sumi = (m+1/2.0)/((m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[n-m]*fkv[m]);  //Z=14638
                            sump = sump+pn[n-m]*sumi;  //Z=14639
                            sump1 = sump1+sumi;  //Z=14640
                        }/*7*/  //Z=14641
                        params.CR->carr5p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14642
                        params.CR->carr6p[n] = (M_PI/4.0)*(1-params.alphash1)*z12v[n]*xrn[n]*sump1;  //Z=14643
                        sump = 0.0;  //Z=14644
                        sump1 = 0.0;  //Z=14645
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14646
                            const double sumi = (n-m+1/2.0)*(m+1/2.0)/((n-m+1/2.0-params.alphash1/2.0)*(m+1/2.0-params.alphash1/2.0)*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=14647
                            sump = sump+sumi;  //Z=14648
                            sump1 = sump1+pn[n-m]*sumi;  //Z=14649
                        }/*7*/  //Z=14650
                        params.CR->carr7p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump;  //Z=14651
                        params.CR->carr8p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrmn_n*sump1;  //Z=14652
                        params.CR->carr9p[n] = (M_PI/4.0)*sqr(1-params.alphash1)*z12v[n]*xrn[n]*sump;  //Z=14653

                        /* (* cross-sectional *)  //Z=14655
                        (* F121 *)  //Z=14656
                         carr4p[n]:=power(4,n)*sqrt(pi)*z12v[n]*xrn[n]/(2*(n+1)*gam3[n]*fkv[n]);  //Z=14657
                         (* F122 *)  //Z=14658
                         sump:=0.0;  //Z=14659
                            for m:=0 to n do begin  //Z=14660
                            sump:=sump+pn[m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14661
                         end;  //Z=14662
                         carr5p[n]:=z12v[n]*xrmn[n]*sump;  //Z=14663
                         (* F122 *)  //Z=14664
                         sump:=0.0;  //Z=14665
                            for m:=0 to n do begin  //Z=14666
                            sump:=sump+pn[m]/(gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=14667
                         end;  //Z=14668
                         carr5p[n]:=pi*z12v[n]*xrmn[n]*sump/4;  //Z=14669
                         (* F123 *)  //Z=14670
                         carr6p[n]:=carr4p[n]/pn[n];    */  //Z=14671

                        /*  F(q)  */  //Z=14673
                        params.CR->carr4f[n] = (sqrt(M_PI)/2.0)*z12v[n]*xrn[n]/(gam3[n]*fkv[n]);  //Z=14674
                        params.CR->carr5f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrmn_n*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14675
                        params.CR->carr6f[n] = (sqrt(M_PI)*(1-params.alphash1)/2.0)*z12v[n]*xrn[n]*(n+1)/((n+1/2.0-params.alphash1/2.0)*gam3[n]*fkv[n]);  //Z=14676

                        /*  series for <...> integration  */  //Z=14678
                        params.CR->carr2p[n] = pow(4.0,n+1)*gam3[n]*z12vl[n]*xln[n]/(sqrt(M_PI)*(n+2)*(n+1)*(n+1)*fkv[n]*fkv[n]*fkv[n]);  //Z=14679
                        sump = 0.0;  //Z=14680
                        for ( int m=0; m<=n; m++ )
                        {/*7*/  //Z=14681
                            sump = sump+z12vl[m]*z12vl[n-m]/((m+1)*fkv[m]*(n-m+1)*fkv[n-m]*fkv[m]*fkv[n-m]);  //Z=14682
                        }/*7*/  //Z=14683
                        params.CR->carr2f[n] = sqrt(M_PI)*fkv[n]*xln[n]*sump/(2.0*gam3[n]);  //Z=14684

                        if ( search1 )
                        {/*7*/  //Z=14686
                            if ( params.CR->carr1p[n]<1e-50 )
                            {/*8*/  //Z=14687
                                n1 = n;  //Z=14688
                                search1 = false;  //Z=14689
                            }/*8*/  //Z=14690
                        }/*7*/  //Z=14691
                        if ( search4 )
                        {/*7*/  //Z=14692
                            if ( params.CR->carr4p[n]<1e-50 )
                            {/*8*/  //Z=14693
                                n4 = n;  //Z=14694
                                search4 = false;  //Z=14695
                            }/*8*/  //Z=14696
                        }/*7*/  //Z=14697
                        if ( search5 )
                        {/*7*/  //Z=14698
                            if ( params.CR->carr5p[n]<1e-50 )
                            {/*8*/  //Z=14699
                                n5 = n;  //Z=14700
                                search5 = false;  //Z=14701
                            }/*8*/  //Z=14702
                        }/*7*/  //Z=14703
                        if ( search6 )
                        {/*7*/  //Z=14704
                            if ( params.CR->carr6p[n]<1e-50 )
                            {/*8*/  //Z=14705
                                n6 = n;  //Z=14706
                                search6 = false;  //Z=14707
                            }/*8*/  //Z=14708
                        }/*7*/  //Z=14709
                        if ( fabs(params.CR->carr7p[n])<min )
                        {/*7*/  //Z=14710
                            if ( n<n7 ) n7 = n;  //Z=14711
                        }/*7*/  //Z=14712
                        if ( fabs(params.CR->carr8p[n])<min )
                        {/*7*/  //Z=14713
                            if ( n<n8 ) n8 = n;  //Z=14714
                        }/*7*/  //Z=14715
                        if ( fabs(params.CR->carr9p[n])<min )
                        {/*7*/  //Z=14716
                            if ( n<n9 ) n9 = n;  //Z=14717
                        }/*7*/  //Z=14718
                        if ( fabs(params.CR->carr1f[n])<min )
                        {/*7*/  //Z=14719
                            if ( n<n1f ) n1f = n;  //Z=14720
                        }/*7*/  //Z=14721
                        if ( fabs(params.CR->carr4f[n])<min )
                        {/*7*/  //Z=14722
                            if ( n<n4f ) n4f = n;  //Z=14723
                        }/*7*/  //Z=14724
                        if ( fabs(params.CR->carr5f[n])<min )
                        {/*7*/  //Z=14725
                            if ( n<n5f ) n5f = n;  //Z=14726
                        }/*7*/  //Z=14727
                        if ( fabs(params.CR->carr6f[n])<min )
                        {/*7*/  //Z=14728
                            if ( n<n6f ) n6f = n;  //Z=14729
                        }/*7*/  //Z=14730
                        if ( fabs(params.CR->carr7f[n])<min )
                        {/*7*/  //Z=14731
                            if ( n<n7f ) n7f = n;  //Z=14732
                        }/*7*/  //Z=14733
                        if ( fabs(params.CR->carr8f[n])<min )
                        {/*7*/  //Z=14734
                            if ( n<n8f ) n8f = n;  //Z=14735
                        }/*7*/  //Z=14736
                        if ( fabs(params.CR->carr9f[n])<min )
                        {/*7*/  //Z=14737
                            if ( n<n9f ) n9f = n;  //Z=14738
                        }/*7*/  //Z=14739
                    }/*6*/  /*  of n-loop  */  //Z=14740
                }/*5*/  /*  of cs=2  */  //Z=14741


            }/*4*/  /*  of disk  */  //Z=14744
        }/*3*/  /*  of z-axis  */  //Z=14745
    }/*2*/  /*  of ordis=0  */  //Z=14746


Label99:
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
    params.limq1 = pow(fabs(params.CR->carr1p[n1]),-1/(2.0*n1));
    params.limq2 = pow(fabs(params.CR->carr2p[n2]),-1/(2.0*n2));
    params.limq3 = pow(fabs(params.CR->carr3p[n3]),-1/(2.0*n3));
    params.limq4 = pow(fabs(params.CR->carr4p[n4]),-1/(2.0*n4));
    params.limq5 = pow(fabs(params.CR->carr5p[n5]),-1/(2.0*n5));
    params.limq6 = pow(fabs(params.CR->carr6p[n6]),-1/(2.0*n6));
    params.limq7 = pow(fabs(params.CR->carr7p[n7]),-1/(2.0*n7));
    params.limq8 = pow(fabs(params.CR->carr8p[n8]),-1/(2.0*n8));
    params.limq9 = pow(fabs(params.CR->carr9p[n9]),-1/(2.0*n9));
    params.limq1f = pow(fabs(params.CR->carr1f[n1f]),-1/(2.0*n1f));  // Im Orginal-Pascalprogramm wurden hier
    params.limq2f = pow(fabs(params.CR->carr2f[n2f]),-1/(2.0*n2f));  //  auch die Variablen n1 bis n9 genutzt
    params.limq3f = pow(fabs(params.CR->carr3f[n3f]),-1/(2.0*n3f));  //  und nicht die jetzt hier stehenden
    params.limq4f = pow(fabs(params.CR->carr4f[n4f]),-1/(2.0*n4f));  //  n1f bis n9f, obwohl diese oben passend
    params.limq5f = pow(fabs(params.CR->carr5f[n5f]),-1/(2.0*n5f));  //  bestimmt worden sind...
    params.limq6f = pow(fabs(params.CR->carr6f[n6f]),-1/(2.0*n6f));  //
    params.limq7f = pow(fabs(params.CR->carr7f[n7f]),-1/(2.0*n7f));  //
    params.limq8f = pow(fabs(params.CR->carr8f[n8f]),-1/(2.0*n8f));  //
    params.limq9f = pow(fabs(params.CR->carr9f[n9f]),-1/(2.0*n9f));  //
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
    // Für Python-Tests hier mal nur das carr4p[] ausgeben, aber komplett
    QString str = "carr4p[]";
    str += QString("%1").arg(params.CR->carr4p[0]);
    for ( int n=1; n<coeffarray_len && n<=120; n++ ) // coeffarray_len=150+1
        str += QString(", %1").arg(params.CR->carr4p[n]);
    qDebug() << str;
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


/**
 * @brief SC_GpuMemory::checkXYarray - helper function to initialize the dynamic result array
 * @param minx - minimal horizontal index (incl.)
 * @param maxx - maximal horizontal index (incl.)
 * @param miny - minimal vertical index (incl.)
 * @param maxy - maximal vertical index (incl.)
 */
void SasCalc_GENERIC_calculation::checkArrays( int minx, int maxx, int miny, int maxy )
{
    _arrCount = (maxx-minx) * (maxy-miny) + 1;
    size_t fullsize = sizeof(double) * (_arrCount + 3); // Es werden die letzten (x,y) Indizes und ein Debugflag mit gespeichert

    if ( arrXYIntensity == nullptr || (_xmax-_xmin) != (maxx-minx) || (_ymax-_ymin) != (maxy-miny) )
    {
        createMemory( (void **)(&arrXYIntensity), fullsize, arrXYsize, false, "xyIntensity" );
        //qDebug() << "checkArrays(" << minx << maxx << miny << maxy << ")" << _arrCount << fullsize;
    }
    // save geometry
    _xmin = minx;
    _xmax = maxx;
    _ymin = miny;
    _ymax = maxy;
}


#ifdef FITDATA_IN_GPU  // real func prog
bool SasCalc_GENERIC_calculation::setFitData( int sx, int sy, const double *data )
{
    _fitWidth  = sx;
    _fitHeight = sy;
    _arrFitSize = sx * sy;
    if ( data == nullptr || _arrFitSize == 0 )
    {
        if ( arrFitData != nullptr ) { memcleanup( arrFitData );   arrFitData = nullptr; }
        if ( arrFitFqs  != nullptr ) { memcleanup( arrFitFqs  );   arrFitFqs  = nullptr; }
        _fitEnabled = false;
        return false;
    }
    createMemory( (void **)(&arrFitData), _arrFitSize * sizeof(double), arrFitDSize, false, "arrFitData" );
    createMemory( (void **)(&arrFitFqs),  _arrFitSize * sizeof(double), arrFitFSize, false, "arrFitFqs" );
    _fitEnabled = arrFitData != nullptr && arrFitFqs != nullptr;
    if ( _fitEnabled )
        memcpy( arrFitData, data, arrFitDSize );

    std::cerr << "Fit-Craete: " << sx << "x" << sy << "=" << _arrFitSize << std::endl;

    return _fitEnabled;
}
#endif


// Nur in prepareCalculation
void SasCalc_GENERIC_calculation::ButtonHKLClick( int ltype ) const
{   //Z=43355

    //const int np=20000;

    int /*i,j,*/ii/*,jj*/,h,k,l,hmax,kmax,lmax;//,index,xpos,ypos,image1width,image1height;
    int /*c0,*/c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,mult,/*multmax,hlev,klev,llev,zlev,zvar,*/ct1;
    // Unbenutzte Variablen: c3, c5, c9, c11, c12, c14 - sie werden aber hier gelassen, da ich sie nicht aus den Zuweisungs-case rausnehmen will
    Q_UNUSED(c3); Q_UNUSED(c5); Q_UNUSED(c9); Q_UNUSED(c11); Q_UNUSED(c12); Q_UNUSED(c14);
    int ccase=0,c1case,c2case/*,c3case*/,c4case/*,c5case*/,c6case,c7case,c8case/*,c9case*/,c10case;//,c11case,c12case;
    float a,b,c,alf,gam; //,xmin,xmax,wave,ttheta,sinarg;
#ifdef UseStringGrid12
    int   ct;
    float q, invd, bet;
#endif
    //float xcord1,ycord1,xcord2,ycord2,xcord3,ycord3,xcord4,ycord4,boxscale,xdir1,ydir1,xdir2,ydir2;
#ifdef UseStringGrid12
    float s2a,c1a,c2a,c3a,c1b,c2b,s2b,c1g,c2g,s2g,scale,amax,bmax,cmax,mx,dhkl;
#endif
    int /*xpos1,ypos1,xpos2,ypos2,*/hmin,kmin,lmin;
    //bool fq,na,nah;
    //std::string trigP;
    char trigP;
    //float nnl,nnal,nnbl,uu;
    //float nn[4],pp[4],nna[4],nnb[4]; // Array[1..3] of real;
    //float qq[31][4],qp[31][4]/*,qs[31][4]*/; // Array[1..30,1..3] of real;

    enum {
        _000,
        _0kl,       //(* one index = zero AND two indices not equal *)
        _h0l,
        _hk0,
        _h00,       //(* two indices = zero *)
        _0k0,
        _00l,
        _hhl,       //(* two indices are equal AND the other is not zero *)
        _hlh,       //(* 'hlh' *)
        _lhh,       //(* 'lhh' *)
        _hh0,       //(* one index = zero   AND   the other two are equal *)
        _h0h,       //(* 'h0l' *)
        _0hh,       //(* '0hh' *)
        _hhh,       //(* all indices are equal *)
        _hkl        //(* most general *)
    } hklcase;

#ifdef UseStringGrid12
    float wave,ttheta,sinarg; //,dhkl;
    QStringList StringGrid12;   // Wird für Debug-Ausgaben verwendet
    QString hklcase2str[] = { "000", "0kl", "h0l", "hk0", "h00", "0k0", "00l", "hhl", "hlh", "lhh", "hh0", "h0h", "0hh", "hhh", "hkl" };
#endif

    // type RealArrayNP = Array[1..np] of real;
    // index1,index2: RealArrayNP; --> werden nicht verwendet, nur an einer Stelle gesetzt...

    // Lokale Prozedur sort() wird nicht verwendet.

    bool    RadioButtonSysCubic=false,
        RadioButtonSysTetragonal=false,
        RadioButtonSys1D=false,
        RadioButtonSys2D=false,
        RadioButtonSysHex=false,
        RadioButtonSysTrigonal=false,   // unused
        RadioButtonSysOrtho=false,      // unused
        RadioButtonSysMonoclinic=false, // unused
        RadioButtonSysTriclinic=false;  // unused
    int     ComboBoxCubic,
        ComboBoxTetragonal,
        ComboBox1D,
        ComboBox2D,
        ComboBoxHexagonal,
        ComboBoxTrigonal=0,     // unused
        ComboBoxOrthorhombic=0, // unused
        ComboBoxMonoclinic=0,   // unused
        ComboBoxTriclinic=0;    // unused

    setLatpar1(  1,0,0);
    setLatpar2(  1,0,0);

    switch ( ltype )
    {
    // Die jeweiligen "RadioButtonSys*Click(Sender);" setzen ComboBox-Auswahlen und Visibility von GUI-Elementen
    // und können somit hier ignoriert werden.
    case  0: /*  Lam  */  //Z=43424
        RadioButtonSys1D = true;  //Z=43425
        //RadioButtonSys1DClick(Sender);  //Z=43426
        ComboBox1D = 0;  //Z=43427
        break;
    case  1: /*  hex cyl  */  //Z=43429
        RadioButtonSys2D = true;  //Z=43430
        //RadioButtonSys2DClick(Sender);  //Z=43431
        ComboBox2D = 12;  //Z=43432
        break;
    case  2: /*  sq Cyl  */  //Z=43434
        RadioButtonSys2D = true;  //Z=43435
        //RadioButtonSys2DClick(Sender);  //Z=43436
        ComboBox2D = 10;  //Z=43437
        break;
    case  3: /*  rec cyl  */  //Z=43439
        RadioButtonSys2D = true;  //Z=43440
        //RadioButtonSys2DClick(Sender);  //Z=43441
        ComboBox2D = 8;  //Z=43442
        break;
    case  4: /*  BCC  */  //Z=43444
        RadioButtonSysCubic = true;  //Z=43445
        //RadioButtonSysCubicClick(Sender);  //Z=43446
        ComboBoxCubic = 7;  //Z=43447
        break;
    case 30: //  Generic
    case  5: /*  FCC  */  //Z=43449
        RadioButtonSysCubic = true;  //Z=43450
        //RadioButtonSysCubicClick(Sender);  //Z=43451
        ComboBoxCubic = 12;  //Z=43452
        break;
    case  6: /*  HCP  */  //Z=43454
        RadioButtonSysHex = true;  //Z=43455
        //RadioButtonSysHexClick(Sender);  //Z=43456
        ComboBoxHexagonal = 0;  //Z=43457
        break;
    case  7: /*  SC  */  //Z=43459
        RadioButtonSysCubic = true;  //Z=43460
        //RadioButtonSysCubicClick(Sender);  //Z=43461
        ComboBoxCubic = 0;  //Z=43462
        break;
    case  8: /*  BCT  */  //Z=43464
        RadioButtonSysTetragonal = true;  //Z=43465
        //RadioButtonSysTetragonalClick(Sender);  //Z=43466
        ComboBoxTetragonal = 23;  //Z=43467
        break;
    case  9: /*  Ia3d  */  //Z=43469
        RadioButtonSysCubic = true;  //Z=43470
        //RadioButtonSysCubicClick(Sender);  //Z=43471
        ComboBoxCubic = 11;  //Z=43472
        break;
    case 10: /*  Pn3m  */  //Z=43474
        RadioButtonSysCubic = true;  //Z=43475
        //RadioButtonSysCubicClick(Sender);  //Z=43476
        ComboBoxCubic = 5;  //Z=43477
        break;
    case 11: /*  Im3m  */  //Z=43479
        RadioButtonSysCubic = true;  //Z=43480
        //RadioButtonSysCubicClick(Sender);  //Z=43481
        ComboBoxCubic = 7;  //Z=43482
        break;
    case 17: /*  Fd3m  */  //Z=43484
        RadioButtonSysCubic = true;  //Z=43485
        //RadioButtonSysCubicClick(Sender);  //Z=43486
        ComboBoxCubic = 15;  //Z=43487
        break;
    case 18: /*  Orthorhombic  */  //Z=43489
        RadioButtonSysOrtho = true;  //Z=43490
        //RadioButtonSysOrthoClick(Sender);  //Z=43491
        ComboBoxCubic = 7;               /*  P121212  */  //Z=43492
        break;
    case 22: /*  Pm3n, A15  */  //Z=43495
        RadioButtonSysOrtho = true;  //Z=43496
        //RadioButtonSysOrthoClick(Sender);  //Z=43497
        ComboBoxCubic = 3;               /*  P__n  */  //Z=43498
        break;

    default:
        return; // TODO Fehlermeldung?
        //break;
    } // switch ltype

// Image?? werden hier komplett ausgeblendet

// StringGrid12 wird als "Debug"-Tabelle angezeigt, daher hier ein paar Infos dazu erhalten
#ifdef UseStringGrid12
    StringGrid12.clear();
    StringGrid12        << "h" << "k" << "l" << "hklcase" << "ccase" << "dhkl/nm" << "mult" << "q" << "2 theta";
    // StringGrid12.Cells  [0,0]  [1,0]  [2,0]  [3,0]        [4,0]      [5,0]        [6,0]     [7,0]  [8,0]
    qDebug() << "SG12:" << StringGrid12.join(" ; ");
#endif

    //(* default hkl-range *)   //{NV}-30759
    hmin=0;
    kmin=0;
    lmin=0;
    hmax=params.hklmax; // StrToInt(Edithklmax.Text);
    kmax=hmax;
    lmax=hmax;

    trigP='h';  //(* initial value *)
    //c0=-1;      //(* R centered points *)
    //c13=-1;     //(* unique axis *)
    //c14=-1;     //(* 2D cells *)

    //(**************** Reflection conditions **************)

    /* ***** 1D systems **** */  //Z=43550
    if ( RadioButtonSys1D )
    {/*2*/  //Z=43551

        /*  lamellar structure  */  //Z=43553
        if ( ComboBox1D==0 )
        {/*3*/ /*  P1 (1)  */  //Z=43554
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c14 = 0;  //Z=43555
            //EditSG.Text = 'LAM';
        }/*3*/  //Z=43556
    }/*2*/  //Z=43557


    /* ***** 2D systems **** */  //Z=43560
    if ( RadioButtonSys2D )
    {
        switch ( ComboBox2D )
        {
            /*  parallelogram (monoclinic)  */  //Z=43563
        case  0: /*  P1 (1)  */  //Z=43564
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c14 = 0;  //Z=43565
            //EditSG.Text = 'P1(1)'; }/*3*/  //Z=43568
            break;
        case  1: /*  P2 (3) */  //Z=43569
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2; c14 = 0;  //Z=43570
            //EditSG.Text = 'P112(3);unique c-axis'; }/*3*/  //Z=43573
            break;

            /*  rectangular (monoclinic, b-axis unique  */  //Z=43575
        case  2: /*  Pm (6)  */  //Z=43576
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1; c14 = 0;  //Z=43577
            //EditSG.Text = 'P1m1(6);unique b-axis'; }/*3*/  //Z=43581
            break;
        case  3: /*  Pg (7) h0l:h, h00:h  */  //Z=43582
            c1 = -1; c2 = -1; c3 = 1; c4 = -1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1; c14 = 1;  //Z=43583
            //EditSG.Text = 'P1a1(7);unique b-axis'; }/*3*/  //Z=43587
            break;
        case  4: /*  Cm (8) hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k  */  //Z=43588
            c1 = 1; c2 = 1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1; c14 = 2;  //Z=43589
            //EditSG.Text = 'C1m1(8);unique b-axis'; }/*3*/  //Z=43593
            break;

            /*  rectangular (orthorhombic)  */  //Z=43595
        case  5: /*  Pmm2 (25)  */  //Z=43596
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c14 = 0;  //Z=43597
            //EditSG.Text = 'Pmm2(25)'; }/*3*/  //Z=43601
            break;
        case  6: /*  Pmg2 (28) 0kl:k, 0k0:k  */  //Z=43602
            c1 = -1; c2 = 1; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c14 = 3;  //Z=43603
            //EditSG.Text = 'Pbm2(28)'; }/*3*/  //Z=43607
            break;
        case  7: /*  Pgg2 (32) 0kl:k, h0l:h, h00:h, 0k0:k  */  //Z=43608
            c1 = -1; c2 = 1; c3 = 1; c4 = -1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c14 = 4;  //Z=43609
            //EditSG.Text = 'Pba2(32)'; }/*3*/  //Z=43613
            break;
        case  8: /*  Cmm2 (35) hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k  */  //Z=43614
            c1 = 1; c2 = 1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c14 = 5;  //Z=43615
            //EditSG.Text = 'Cmm2(35)'; }/*3*/  //Z=43619
            break;

            /*  trigonal  */  //Z=43621
        case  9: /*  P3-- (143,156,157)  */  //Z=43622
            /*c0=0;*/ c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1; c14 = 0;  //Z=43623
            //EditSG.Text = 'P3(143),P3m1(156),P31m(157)'; trigP = 'h'; }/*3*/  //Z=43628
            break;

            /*  tetragonal  */  //Z=43630
        case 10: /*  P4-- (75,99)  */  //Z=43631
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0; c14 = 0;  //Z=43632
            //EditSG.Text = 'P4(75),P4mm(99)'; }/*3*/  //Z=43637
            break;
        case 11: /*  P4gm (100) 0kl:k, h0l:h, 0k0:k  */  //Z=43638
            c1 = -1; c2 = 1; c3 = 1; c4 = -1; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0; c14 = 6;  //Z=43639
            //EditSG.Text = 'P4bm(100)'; }/*3*/  //Z=43644
            break;

            /*  hexagonal  */  //Z=43646
        case 12: /*  P6-- (168,183)  */  //Z=43647
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1; c14 = 0;  //Z=43648
            //EditSG.Text = 'P6(168),P6mm(183)'; }/*3*/  //Z=43653
            break;
        } // switch ComboBox2D
    } // if RadioButtonSys2D


    /* *************** Cubic crystal system ************** */  //Z=43657
    if ( RadioButtonSysCubic )
    {
        switch ( ComboBoxCubic )
        {
        /*  P___  */  //Z=43660
        case  0: /*  P--- |  */  //Z=43661
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 2;  //Z=43662
            //EditSG.Text = 'P23(195),Pm-3(200),P432(207),P-43m(215),Pm-3m(221)'; }/*3*/  //Z=43663
            break;
        case  1: /*  P21--, P42-- | 00l:l  */  //Z=43664
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0;  c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 2;  //Z=43665
            //EditSG.Text = 'P213(198),P4232(208)'; }/*3*/  //Z=43666
            break;
        case  2: /*  P41-- | 00l:4n  */  //Z=43667
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1;  c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 2;  //Z=43668
            //EditSG.Text = 'P4132(213),P4332(212)'; }/*3*/  //Z=43669
            break;
        case  3: /*  P--n | hhl:l, 00l:l  */  //Z=43670
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0;  c8 = -1; c9 = -1; c10 = 0;  c11 = -1; c12 = 2;  //Z=43671
            //EditSG.Text = 'P-43n(218),Pm-3n(223)'; }/*3*/  //Z=43672
            break;
        case  4: /*  Pa-- | cyclic permutations: 0kl:k, h0l:l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=43673
            c1 = -1; c2 = 5;  c3 = -1;  c4 = -1;  c5 = 0;  c6 = 0;  c7 = 0;  c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 2;  //Z=43674
            //EditSG.Text = 'Pa-3(205)'; }/*3*/  //Z=43675
            break;
        case  5: /*  Pn-- | 0kl:k+l, 00l:l  */  //Z=43676
            c1 = -1; c2 = 0;  c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0;  c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 2;  //Z=43677
            //EditSG.Text = 'Pn-3(201),Pn-3m(224)'; }/*3*/  //Z=43678
            break;
        case  6: /*  Pn-n | 0kl:k+l, hhl:l, 00l:l  */  //Z=43679
            c1 = -1; c2 = 0;  c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0;  c8 = -1; c9 = -1; c10 = 0;  c11 = -1; c12 = 2;  //Z=43680
            //EditSG.Text = 'Pn-3n(222)'; }/*3*/  //Z=43681
            break;

        /*  I___  */  //Z=43683
        case  7: /*  I--- | hkl:h+k+l, 0kl:k+l, hhl:l, 00l:l  */  //Z=43684
            c1 = 0; c2 = 0;  c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0;  c8 = -1; c9 = -1; c10 = 0;  c11 = -1; c12 = 2;  //Z=43685
            //EditSG.Text = 'I23(197),I213(199),Im-3(204),I432(211),I-43m(217),Im-3m(229)'; }/*3*/  //Z=43686
            break;
        case  8: /*  I41-- | hkl:h+k+l, 0kl:k+l, hhl:l, 00l:4n  */  //Z=43687
            c1 = 0;  c2 = 0; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1;  c8 = -1; c9 = -1; c10 = 0;  c11 = -1; c12 = 2;  //Z=43688
            //EditSG.Text = 'I4132(214)'; }/*3*/  //Z=43689
            break;
        case  9: /*  I--d | hkl:h+k+l, 0kl:k+l, hhl:2h+l=4n,l, 00l:4n  */  //Z=43690
            c1 = 0; c2 = 0;  c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1;  c8 = -1; c9 = -1; c10 = 4;  c11 = -1; c12 = 2;  //Z=43691
            //EditSG.Text = 'I-43d(220)'; }/*3*/  //Z=43692
            break;
        case 10: /*  Ia-- | hkl:h+k+l, 0kl:k,l, hhl:l, 00l:l  */  //Z=43693
            c1 = 0;  c2 = 4; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0;  c8 = -1; c9 = -1; c10 = 0;  c11 = -1; c12 = 2;  //Z=43694
            //EditSG.Text = 'Ia-3(206)'; }/*3*/  //Z=43695
            break;
        case 11: /*  Ia-d | hkl:h+k+l, 0kl:k,l, hhl:2h+l=4n,l, 00l:4n  */  //Z=43696
            c1 = 0;  c2 = 4; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1; c8 = -1; c9 = -1;  c10 = 4;  c11 = -1; c12 = 2;  //Z=43697
            //EditSG.Text = 'Ia-3d(230)'; }/*3*/  //Z=43698
            break;

        /*  F___  */  //Z=43700
        case 12: /*  F--- | hkl:h+k,h+l,k+l, 0kl:k,l, hhl:h+l, 00l:l  */  //Z=43701
            c1 = 6;  c2 = 4; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = 2;  c11 = -1; c12 = 2;  //Z=43702
            //EditSG.Text = 'F23(196),Fm-3(202),F432(209),F-43m(216),Fm-3m(225)'; }/*3*/  //Z=43703
            break;
        case 13: /*  F41-- | hkl:h+k,h+l,k+l, 0kl: k,l, hhl:h+l, 00l:4n  */  //Z=43704
            c1 = 6; c2 = 4; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1;  c8 = -1; c9 = -1; c10 = 2; c11 = -1; c12 = 2;  //Z=43705
            //EditSG.Text = 'F4132(210)'; }/*3*/  //Z=43706
            break;
        case 14: /*  F--c | hkl:h+k,h+l,k+l, 0kl:k,l, hhl:h,l, 00l:l  */  //Z=43707
            c1 = 6; c2 = 4; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = 5; c11 = -1; c12 = 2;  //Z=43708
            //EditSG.Text = 'F-43c(219),Fm-3c(226)'; }/*3*/  //Z=43709
            break;
        case 15: /*  Fd-- | hkl:h+k,h+l,k+l, 0kl:k+l=4n,k,l, hhl:h+l, 00l:4n  */  //Z=43710
            c1 = 6; c2 = 3; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1; c8 = -1; c9 = -1; c10 = 2; c11 = -1; c12 = 2;  //Z=43711
            //EditSG.Text = 'Fd-3(203),Fd-3m(227)'; }/*3*/  //Z=43712
            break;
        case 16: /*  Fd-c | hkl:h+k,h+l,k+l, 0kl:k+l=4n,k,l, hhl:h,l, 00l:4n  */  //Z=43713
            c1 = 6; c2 = 3; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1; c8 = -1; c9 = -1; c10 = 5; c11 = -1; c12 = 2;  //Z=43714
            //EditSG.Text = 'Fd-3c(228)'; }/*3*/  //Z=43715
            break;
        } // switch ComboBoxCubic
    } // if RadioButtonSysCubic


    /* *************** Hexagonal crystal system ************** */  //Z=43718
    if ( RadioButtonSysHex )
    {
        hmin = -hmax;  //Z=43720
        kmin = -kmax;  //Z=43721
        lmin = -lmax;  //Z=43722
        switch ( ComboBoxHexagonal )
        {
        /*  P___  */  //Z=43723
        case  0: /*  P--- |  */  //Z=43724
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1;  //Z=43725
            //EditSG.Text = 'P6(168),P-6(174),P6/m(175),P622(177),P6mm(183),P-62m(189),P-6m2(187),P6/mmm(191)'; }/*3*/  //Z=43726
            break;
        case  1: /*  P63-- | 00l:l  */  //Z=43727
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1;  //Z=43728
            //EditSG.Text = 'P63(173),P63/m(176),P6322(182)'; }/*3*/  //Z=43729
            break;
        case  2: /*  P62-- | 00l:3n  */  //Z=43730
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 2; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1;  //Z=43731
            //EditSG.Text = 'P62(171),P64(172),P6222(180),P6422(181)'; }/*3*/  //Z=43732
            break;
        case  3: /*  P61-- | 00l:6n  */  //Z=43733
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 3; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1;  //Z=43734
            //EditSG.Text = 'P61(169),P63(170),P6122(178),P6522(179)'; }/*3*/  //Z=43735
            break;
        case  4: /*  P--c | hhl:l, 00l:l  */  //Z=43736
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 1;  //Z=43737
            //EditSG.Text = 'P63mc(186),P-62c(190),P63/mmc(194)'; }/*3*/  //Z=43738
            break;
        case  5: /*  P-c- | h-hl:l, 00l:l  */  //Z=43739
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = 0; c10 = -1; c11 = -1; c12 = 1;  //Z=43740
            //EditSG.Text = 'P63cm(185),P-6c2(188),P63/mcm(193)'; }/*3*/  //Z=43741
            break;
        case  6: /*  P-cc | h-hl:l, hhl:l, 00l:l  */  //Z=43742
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = 0; c10 = 0; c11 = -1; c12 = 1;  //Z=43743
            //EditSG.Text = 'P6cc(184),P6/mcc(192)'; }/*3*/  //Z=43744
            break;
        } // switch ComboBoxHexagonal
    } // if RadioButtonSysHex


    /* *************** Trigonal crystal system ************** */  //Z=43747
    if ( RadioButtonSysTrigonal )
    {
        hmin = -hmax;  //Z=43749
        kmin = -kmax;  //Z=43750
        lmin = -lmax;  //Z=43751
        switch ( ComboBoxTrigonal )
        {
            /*  hexagonal axes, P___  */  //Z=43752
        case  0: /*  P--- |  */  //Z=43753
            /*c0=0;*/ c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1;  //Z=43754
            //EditSG.Text = 'P3(143),P-3(147),P321(150),P3m1(156),P-3m1(164),P312(149),P31m(157),P-31m(162)'; }/*3*/  //Z=43755
            trigP = 'h';
            break;
        case  1: /*  P31-- | 00l:3n  */  //Z=43756
            /*c0=0;*/ c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 2; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1;  //Z=43757
            //EditSG.Text = 'P31(144),P32(145),P3121(152),P3221(154),P3112(151),P3212(153)'; }/*3*/  //Z=43758
            trigP = 'h';
            break;
        case  2: /*  P--c | hhl:l, 00l:l  */  //Z=43759
            /*c0=0;*/ c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 1;  //Z=43760
            //EditSG.Text = 'P31c(159),P-31c(163)'; }/*3*/  //Z=43761
            trigP = 'h';
            break;
        case  3: /*  P-c- | h-hl:l, 00l:l  */  //Z=43762
            /*c0=0;*/ c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = 0; c10 = -1; c11 = -1; c12 = 1;  //Z=43763
            //EditSG.Text = 'P3c1(158),P-3c1(165)'; }/*3*/  //Z=43764
            trigP = 'h';
            break;

            /*  hexagonal axes, R___  */  //Z=43766
        case  4: /*  R(obv)-- | hkl:-h+k+l=3n, h-hl:h+l=3n, hhl:l=3n, 00l:3n  */  //Z=43767
            /*c0=1;*/ c1 = 4; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 2; c8 = -1; c9 = 1; c10 = 1; c11 = -1; c12 = 1;  //Z=43768
            //EditSG.Text = 'R3(146),R-3(148),R32(155),R3m(160),R-3m(166)'; }/*3*/  //Z=43769
            trigP = 'h';
            break;
        case  5: /*  R(obv)-c | hkl:-h+k+l=3n, h-hl:h+l=3n,l, hhl:l=3n, 00l:6n  */  //Z=43770
            /*c0=1;*/ c1 = 4; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 3; c8 = -1; c9 = 3; c10 = 1; c11 = -1; c12 = 1;  //Z=43771
            //EditSG.Text = 'R3c(161),R-3c(148)'; trigP = 'h'; }/*3*/  //Z=43772
            break;
        case  6: /*  R(rev)-- | hkl:h-k+l=3n, h-hl:-h+l=3n, hhl:l=3n, 00l:3n  */  //Z=43773
            /*c0=2;*/ c1 = 5; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 2; c8 = -1; c9 = 2; c10 = 1; c11 = -1; c12 = 1;  //Z=43774
            //EditSG.Text = 'R3(146),R-3(148),R32(155),R3m(160),R-3m(166)'; }/*3*/  //Z=43775
            trigP = 'h';
            break;
        case  7: /*  R(rev)-c | hkl:h-k+l=3n, h-hl:-h+l=3n,l, hhl:l=3n, 00l:6n  */  //Z=43776
            /*c0=2;*/ c1 = 5; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 3; c8 = -1; c9 = 3; c10 = 1; c11 = 0; c12 = 1;  //Z=43777
            //EditSG.Text = 'R3c(161),R-3c(167)'; }/*3*/  //Z=43778
            trigP = 'h';
            break;

            /*  rhombohedral axes, R___  */  //Z=43780
        case  8: /*  R-- |   */  //Z=43781
            /*c0=-1;*/ c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 1;  //Z=43782
            //EditSG.Text = 'R3(146),R-3(148),R32(155),R3m(160),R-3m(166); rhombohedral'; }/*3*/  //Z=43783
            trigP='r';
            break;
        case  9: /*  R-c | hhl:l, hhh:h  */  //Z=43784
            /*c0=-1;*/ c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = 0; c11 = 0; c12 = 1;  //Z=43785
            //EditSG.Text = 'R3c(161),R-3c(167); rhombohedral'; }/*3*/  //Z=43786
            trigP = 'r';
            break;
        } // switch ComboBoxTrigonal
    } // if RadioButtonSysTrigonal


    /* *************** Tetragonal crystal system ************** */  //Z=43790
    if ( RadioButtonSysTetragonal )
    {
        /* hmin:=0;  //Z=43792 */
        /* kmin:=0;  //Z=43793 */
        /* lmin:=0;  //Z=43794 */
        /*  P___ */  //Z=43795
        switch ( ComboBoxTetragonal )
        {
        case  0: /*  P--- |  */  //Z=43796
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43797
            //EditSG.Text = 'P4(75),P-4(81),P4/m(83),P422(89),P4mm(99),P-42m(111),P-4m2(115),P4/mmm(123)'; }/*3*/  //Z=43798
            break;
        case  1: /*  P-21- | 0k0:k &permu: h00:h  */  //Z=43799
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43800
            //EditSG.Text = 'P4212(90),P-421m(113)'; }/*3*/  //Z=43801
            break;
        case  2: /*  P42-- | 00l:l  */  //Z=43802
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43803
            //EditSG.Text = 'P42(77),P42/m(84),P4222(93)'; }/*3*/  //Z=43804
            break;
        case  3: /*  P4221- | 00l:l, 0k0:k &permu: h00:h  */  //Z=43805
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43806
            //EditSG.Text = 'P42212(94)'; }/*3*/  //Z=43807
            break;
        case  4: /*  P41-- | 00l:4n  */  //Z=43808
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43809
            //EditSG.Text = 'P41(76),P43(78),P4122(91),P4322(95)'; }/*3*/  //Z=43810
            break;
        case  5: /*  P4121- | 00l:4n, 0k0:k &permu: h00:h  */  //Z=43811
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = 1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43812
            //EditSG.Text = 'P41212(92),P43212(96)'; }/*3*/  //Z=43813
            break;
        case  6: /*  P--c | hhl:l, 00l:l  */  //Z=43814
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43815
            //EditSG.Text = 'P42mc(105),P-42c(112),P42/mmc(131)'; }/*3*/  //Z=43816
            break;
        case  7: /*  P-21c | hhl:l, 00l:l, 0k0:k &permu: h00:h  */  //Z=43817
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43818
            //EditSG.Text = 'P-421c(114)'; }/*3*/  //Z=43819
            break;
        case  8: /*  P-b- | 0kl:k, 0k0:k &permu: h0l:h  */  //Z=43820
            c1 = -1; c2 = 1; c3 = 1; c4 = -1; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43821
            //EditSG.Text = 'P4bm(100),P-4b2(117),P4/mbm(127)'; }/*3*/  //Z=43822
            break;
        case  9: /*  P-bc | 0kl:k, hhl:l, 00l:l, 0k0:k &permu: h0l:h  */  //Z=43823
            c1 = -1; c2 = 1; c3 = 1; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43824
            //EditSG.Text = 'P42bc(106),P42/mbc(135)'; }/*3*/  //Z=43825
            break;
        case 10: /*  P-c- | 0kl:l, 00l:l &permu: h0l:l */  //Z=43826
            c1 = -1; c2 = 2; c3 = 2; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43827
            //EditSG.Text = 'P42cm(101),P-4c2(116),P42mcm(132)'; }/*3*/  //Z=43828
            break;
        case 11: /*  P-cc | 0kl:l, hhl:l, 00l:l &permu: h0l:l  */  //Z=43829
            c1 = -1; c2 = 2; c3 = 2; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43830
            //EditSG.Text = 'P4cc(103),P4/mcc(124)'; }/*3*/  //Z=43831
            break;
        case 12: /*  P-n- | 0kl:k+l, 00l:l, 0k0:k &permu: h0l:h+l  */  //Z=43832
            c1 = -1; c2 = 0; c3 = 0; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43833
            //EditSG.Text = 'P42nm(102),P-4n2(118),P42/mnm(136)'; }/*3*/  //Z=43834
            break;
        case 13: /*  P-nc | 0kl:k+l, hhl:l, 00l:l, 0k0:k &permu: h0l:h+l  */  //Z=43835
            c1 = -1; c2 = 0; c3 = 0; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43836
            //EditSG.Text = 'P4nc(104),P4/mnc(128)'; }/*3*/  //Z=43837
            break;
        case 14: /*  Pn-- | hk0:h+k, 0k0:k  */  //Z=43838
            c1 = -1; c2 = -1; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43839
            //EditSG.Text = 'P4/n(85),P4/nmm(129)'; }/*3*/  //Z=43840
            break;
        case 15: /*  P42/n-- | hk0:h+k, 00l:l, 0k0:k  */  //Z=43841
            c1 = -1; c2 = -1; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43842
            //EditSG.Text = 'P42/n(86)'; }/*3*/  //Z=43843
            break;
        case 16: /*  Pn-c | hk0:h+k, hhl:l, 00l:l, 0k0:k  */  //Z=43844
            c1 = -1; c2 = -1; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43845
            //EditSG.Text = 'P42/nmc(137)'; }/*3*/  //Z=43846
            break;
        case 17: /*  Pnb- | hk0:h+k, 0kl:k, 0k0:k &permu: h0l:h  */  //Z=43847
            c1 = -1; c2 = 1; c3 = 1; c4 = 0; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43848
            //EditSG.Text = 'P4/nbm(125)'; }/*3*/  //Z=43849
            break;
        case 18: /*  Pnbc | hk0:h+k, 0kl:k, hhl:l, 00l:l, 0k0:k &permu: h0l:h */  //Z=43850
            c1 = -1; c2 = 1; c3 = 1; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43851
            //EditSG.Text = 'P42/nbc(133)'; }/*3*/  //Z=43852
            break;
        case 19: /*  Pnc- | hk0:h+k, 0kl:l, 00l:l, 0k0:k &permu: h0l:l  */  //Z=43853
            c1 = -1; c2 = 2; c3 = 2; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43854
            //EditSG.Text = 'P42/ncm(138)'; }/*3*/  //Z=43855
            break;
        case 20: /*  Pncc | hk0:h+k, 0kl:l, hhl:l, 00l:l, 0k0:k &permu: h0l:l  */  //Z=43856
            c1 = -1; c2 = 2; c3 = 2; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43857
            //EditSG.Text = 'P4/ncc(130)'; }/*3*/  //Z=43858
            break;
        case 21: /*  Pnn- | hk0:h+k, 0kl:k+l, 00l:l, 0k0:k  */  //Z=43859
            c1 = -1; c2 = 0; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = 0;  //Z=43860
            //EditSG.Text = 'P42/nnm(134)'; }/*3*/  //Z=43861
            break;
        case 22: /*  Pnnc | hk0:h+k, 0kl:k+l, hhl:l, 00l:l, 0k0:k  */  //Z=43862
            c1 = -1; c2 = 0; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43863
            //EditSG.Text = 'P4/ncc(126)'; }/*3*/  //Z=43864
            break;

            /*  I___  */  //Z=43866
        case 23: /*  I--- | hkl:h+k+l, hk0:h+k, 0kl:k+l, hhl:l, 00l:l, 0k0:k  */  //Z=43867
            c1 = 0; c2 = 0; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43868
            //EditSG.Text = 'I4(79),I-4(82),I4/m(87),I422(97),I4mm(107),I-42m(121),I-4m2(119),I4/mmm(13)'; }/*3*/  //Z=43869
            break;
        case 24: /*  I41-- | hkl:h+k+l, hk0:h+k, 0kl:k+l, hhl:l, 00l:4n, 0k0:k  */  //Z=43870
            c1 = 0; c2 = 0; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 1; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43871
            //EditSG.Text = 'I41(80),I4122(98)'; }/*3*/  //Z=43872
            break;
        case 25: /*  I--d | hkl:h+k+l, hk0:h+k, 0kl:k+l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h  */  //Z=43873
            c1 = 0; c2 = 0; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 1; c8 = 0; c9 = -1; c10 = 4; c11 = -1; c12 = 0;  //Z=43874
            //EditSG.Text = 'I41md(109),I-42d(122)'; }/*3*/  //Z=43875
            break;
        case 26: /*  I-c- | hkl:h+k+l, hk0:h+k, 0kl:k,l, hhl:l, 00l:l, 0k0:k  */  //Z=43876
            c1 = 0; c2 = 4; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43877
            //EditSG.Text = 'I4cm(108),I-4c2(120),I4/mcm(140)'; }/*3*/  //Z=43878
            break;
        case 27: /*  I-cd | hkl:h+k+l, hk0:h+k, 0kl:k,l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h  */  //Z=43879
            c1 = 0; c2 = 4; c3 = -1; c4 = 0; c5 = -1; c6 = 0; c7 = 1; c8 = 0; c9 = -1; c10 = 4; c11 = -1; c12 = 0;  //Z=43880
            //EditSG.Text = 'I41cd(110)'; }/*3*/  //Z=43881
            break;
        case 28: /*  I41/a-- | hkl:h+k+l, hk0:h,k, 0kl:k+l, hhl:l, 00l:4n, 0k0:k  */  //Z=43882
            c1 = 0; c2 = 0; c3 = -1; c4 = 4; c5 = -1; c6 = 0; c7 = 1; c8 = -1; c9 = -1; c10 = 0; c11 = -1; c12 = 0;  //Z=43883
            //EditSG.Text = 'I41/a(88)'; }/*3*/  //Z=43884
            break;
        case 29: /*  Ia-d | hkl:h+k+l, hk0:h,k, 0kl:k+l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h  */  //Z=43885
            c1 = 0; c2 = 0; c3 = -1; c4 = 4; c5 = -1; c6 = 0; c7 = 1; c8 = 0; c9 = -1; c10 = 4; c11 = -1; c12 = 0;  //Z=43886
            //EditSG.Text = 'I41/amd(141)'; }/*3*/  //Z=43887
            break;
        case 30: /*  Iacd | hkl:h+k+l, hk0:h,k, 0kl:k,l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h  */  //Z=43888
            c1 = 0; c2 = 4; c3 = -1; c4 = 4; c5 = -1; c6 = 0; c7 = 1; c8 = 0; c9 = -1; c10 = 4; c11 = -1; c12 = 0;  //Z=43889
            //EditSG.Text = 'I41/acd(142)'; }/*3*/  //Z=43890
            break;
        } // switch ( ComboBoxTetragonal )
    } // if ( RadioButtonSysTetragonal )


    /* *************** Orthorhombic crystal system ************** */  //Z=43893
    if ( RadioButtonSysOrtho )
    {
        /* hmin:=0;  //Z=43895 */
        /* kmin:=0;  //Z=43896 */
        /* lmin:=0;  //Z=43897 */

        switch ( ComboBoxOrthorhombic ) // 0 .. 110
        {
            /*  primitive  */  //Z=43899
        case   0: /*  P--- |  */  //Z=43900
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43901
            //EditSG.Text = 'P222(16),Pmm2(25),Pm2m(25),P2mm(25),Pmmm(47)'; }/*3*/  //Z=43902
            break;
        case   1: /*  P--21 | 00l:l  */  //Z=43903
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43904
            //EditSG.Text = 'P2221(17)'; }/*3*/  //Z=43905
            break;
        case   2: /*  P-21- | 0k0:k  */  //Z=43906
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43907
            //EditSG.Text = 'P2212(17)'; }/*3*/  //Z=43908
            break;
        case   3: /*  P-2121 | 0k0:k, 00l:l  */  //Z=43909
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43910
            //EditSG.Text = 'P22121(18)'; }/*3*/  //Z=43911
            break;
        case   4: /*  P21-- | h00:h  */  //Z=43912
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43913
            //EditSG.Text = 'P2122(17)'; }/*3*/  //Z=43914
            break;
        case   5: /*  P21-21 | h00:h, 00l:l  */  //Z=43915
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43916
            //EditSG.Text = 'P21221(18)'; }/*3*/  //Z=43917
            break;
        case   6: /*  P2121- | h00:h, 0k0:k  */  //Z=43918
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43919
            //EditSG.Text = 'P21212(18)'; }/*3*/  //Z=43920
            break;
        case   7: /*  P212121 | h00:h, 0k0:k, 00l:l  */  //Z=43921
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43922
            //EditSG.Text = 'P212121(19)'; }/*3*/  //Z=43923
            break;
        case   8: /*  P--a | hk0:h, h00:h  */  //Z=43924
            c1 = -1; c2 = -1; c3 = -1; c4 = 1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43925
            //EditSG.Text = 'Pm2a(28),P21ma(26),Pmma(51)'; }/*3*/  //Z=43926
            break;
        case   9: /*  P--b | hk0:k, 0k0:k  */  //Z=43927
            c1 = -1; c2 = -1; c3 = -1; c4 = 2; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43928
            //EditSG.Text = 'Pm21b(26),P2mb(28),Pmmb(51)'; }/*3*/  //Z=43929
            break;
        case  10: /*  P--n | hk0:h+k, h00:h, 0k0:k  */  //Z=43930
            c1 = -1; c2 = -1; c3 = -1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43931
            //EditSG.Text = 'Pm21n(31),P21mn(31),Pmmn(59)'; }/*3*/  //Z=43932
            break;
        case  11: /*  P-a- | h0l:h, h00:h  */  //Z=43933
            c1 = -1; c2 = -1; c3 = 1; c4 = -1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43934
            //EditSG.Text = 'Pma2(28),P21am(26),Pmam(51)'; }/*3*/  //Z=43935
            break;
        case  12: /*  P-aa | h0l:h, hk0:h, h00:h  */  //Z=43936
            c1 = -1; c2 = -1; c3 = 1; c4 = 1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43937
            //EditSG.Text = 'P2aa(27),Pmaa(49)'; }/*3*/  //Z=43938
            break;
        case  13: /*  P-ab | h0l:h, hk0:k, h00:h, 0k0:k  */  //Z=43939
            c1 = -1; c2 = -1; c3 = 1; c4 = 1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43940
            //EditSG.Text = 'P21ab(29),Pmab(57)'; }/*3*/  //Z=43941
            break;
        case  14: /*  P-an | h0l:h, hk0:h+k, h00:h, 0k0:k  */  //Z=43942
            c1 = -1; c2 = -1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43943
            //EditSG.Text = 'P2an(30),Pman(53)'; }/*3*/  //Z=43944
            break;
        case  15: /*  P-c- | h0l:l, 00l:l  */  //Z=43945
            c1 = -1; c2 = -1; c3 = 2; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43946
            //EditSG.Text = 'Pmc21(26),P2cm(28),Pmcm(51)'; }/*3*/  //Z=43947
            break;
        case  16: /*  P-ca | h0l:l, hk0:h, h00:h, 00l:l  */  //Z=43948
            c1 = -1; c2 = -1; c3 = 2; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43949
            //EditSG.Text = 'P21ca(29),Pmca(57)'; }/*3*/  //Z=43950
            break;
        case  17: /*  P-cb | h0l:l, hk0:k, 0k0:k, 00l:l  */  //Z=43951
            c1 = -1; c2 = -1; c3 = 2; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43952
            //EditSG.Text = 'P2cb(32),Pmcb(55)'; }/*3*/  //Z=43953
            break;
        case  18: /*  P-cn | h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=43954
            c1 = -1; c2 = -1; c3 = 2; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43955
            //EditSG.Text = 'P21cn(33),Pmcn(62)'; }/*3*/  //Z=43956
            break;
        case  19: /*  P-n- | h0l:h+l, h00:h, 00l:l  */  //Z=43957
            c1 = -1; c2 = -1; c3 = 0; c4 = -1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43958
            //EditSG.Text = 'Pmn21(31),P21nm(31),Pmnm(59)'; }/*3*/  //Z=43959
            break;
        case  20: /*  P-na | h0l:h+l, hk0:h, h00:h, 00l:l  */  //Z=43960
            c1 = -1; c2 = -1; c3 = 0; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43961
            //EditSG.Text = 'P2na(30),Pmna(53)'; }/*3*/  //Z=43962
            break;
        case  21: /*  P-nb | h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=43963
            c1 = -1; c2 = -1; c3 = 0; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43964
            //EditSG.Text = 'P21nb(33),Pmnb(62)'; }/*3*/  //Z=43965
            break;
        case  22: /*  P-nn | h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=43966
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43967
            //EditSG.Text = 'P2nn(34),Pmnn(58)'; }/*3*/  //Z=43968
            break;
        case  23: /*  Pb-- | 0kl:k, 0k0:k  */  //Z=43969
            c1 = -1; c2 = 1; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43970
            //EditSG.Text = 'Pbm2(28),Pb21m(26),Pbmm(51)'; }/*3*/  //Z=43971
            break;
        case  24: /*  Pb-a | 0kl:k, hk0:h, h00:h, 0k0:k  */  //Z=43972
            c1 = -1; c2 = 1; c3 = -1; c4 = 1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43973
            //EditSG.Text = 'Pb21a(29),Pbma(57)'; }/*3*/  //Z=43974
            break;
        case  25: /*  Pb-b | 0kl:k, hk0:k, 0k0:k  */  //Z=43975
            c1 = -1; c2 = 1; c3 = -1; c4 = 2; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43976
            //EditSG.Text = 'Pb2b(27),Pbmb(49)'; }/*3*/  //Z=43977
            break;
        case 26: /*  Pb-n | 0kl:k, hk0:h+k, h00:h, 0k0:k  */  //Z=43978
            c1 = -1; c2 = 1; c3 = -1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43979
            //EditSG.Text = 'Pb2n(30),Pbmn(53)'; }/*3*/  //Z=43980
            break;
        case 27: /*  Pba- | 0kl:k, h0l:h, h00:h, 0k0:k  */  //Z=43981
            c1 = -1; c2 = 1; c3 = 1; c4 = -1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43982
            //EditSG.Text = 'Pba2(32),Pbam(55)'; }/*3*/  //Z=43983
            break;
        case 28: /*  Pbaa | 0kl:k, h0l:h, hk0:h, h00:h, 0k0:k  */  //Z=43984
            c1 = -1; c2 = 1; c3 = 1; c4 = 1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43985
            //EditSG.Text = 'Pbaa(54)'; }/*3*/  //Z=43986
            break;
        case 29: /*  Pbab | 0kl:k, h0l:h, hk0:k, h00:h, 0k0:k  */  //Z=43987
            c1 = -1; c2 = 1; c3 = 1; c4 = 2; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43988
            //EditSG.Text = 'Pbab(54)'; }/*3*/  //Z=43989
            break;
        case 30: /*  Pban | 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k  */  //Z=43990
            c1 = -1; c2 = 1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43991
            //EditSG.Text = 'Pban(50)'; }/*3*/  //Z=43992
            break;
        case 31: /*  Pbc- | 0kl:k, h0l:l, 0k0:k, 00l:l  */  //Z=43993
            c1 = -1; c2 = 1; c3 = 2; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43994
            //EditSG.Text = 'Pbc21(29),Pbcm(57)'; }/*3*/  //Z=43995
            break;
        case 32: /*  Pbca | 0kl:k, h0l:l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=43996
            c1 = -1; c2 = 1; c3 = 2; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=43997
            //EditSG.Text = 'Pbca(61)'; }/*3*/  //Z=43998
            break;
        case 33: /*  Pbcb | 0kl:k, h0l:l, hk0:k, 0k0:k, 00l:l  */  //Z=43999
            c1 = -1; c2 = 1; c3 = 2; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44000
            //EditSG.Text = 'Pbcb(54)'; }/*3*/  //Z=44001
            break;
        case 34: /*  Pbcn | 0kl:k, h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44002
            c1 = -1; c2 = 1; c3 = 2; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44003
            //EditSG.Text = 'Pbcn(60)'; }/*3*/  //Z=44004
            break;
        case 35: /*  Pbn- | 0kl:k, h0l:h+l, h00:h, 0k0:k, 00l:l  */  //Z=44005
            c1 = -1; c2 = 1; c3 = 0; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44006
            //EditSG.Text = 'Pbn21(33),Pbnm(62)'; }/*3*/  //Z=44007
            break;
        case 36: /*  Pbna | 0kl:k, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44008
            c1 = -1; c2 = 1; c3 = 0; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44009
            //EditSG.Text = 'Pbna(60)'; }/*3*/  //Z=44010
            break;
        case 37: /*  Pbnb | 0kl:k, h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44011
            c1 = -1; c2 = 1; c3 = 0; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44012
            //EditSG.Text = 'Pbnb(56)'; }/*3*/  //Z=44013
            break;
        case 38: /*  Pbnn | 0kl:k, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44014
            c1 = -1; c2 = 1; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44015
            //EditSG.Text = 'Pbnn(52)'; }/*3*/  //Z=44016
            break;
        case 39: /*  Pc-- | 0kl:l, 00l:l  */  //Z=44017
            c1 = -1; c2 = 2; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44018
            //EditSG.Text = 'Pcm21(26),Pc2m(28),Pcmm(51)'; }/*3*/  //Z=44019
            break;
        case 40: /*  Pc-a | 0kl:l, hk0:h, h00:h, 00l:l  */  //Z=44020
            c1 = -1; c2 = 2; c3 = -1; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44021
            //EditSG.Text = 'Pc2a(32),Pcma(55)'; }/*3*/  //Z=44022
            break;
        case 41: /*  Pc-b | 0kl:l, hk0:k, 0k0:k, 00l:l  */  //Z=44023
            c1 = -1; c2 = 2; c3 = -1; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44024
            //EditSG.Text = 'Pc21b(29),Pcmb(57)'; }/*3*/  //Z=44025
            break;
        case 42: /*  Pc-n | 0kl:l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44026
            c1 = -1; c2 = 2; c3 = -1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44027
            //EditSG.Text = 'Pc21n(33),Pcmn(62)'; }/*3*/  //Z=44028
            break;
        case 43: /*  Pca- | 0kl:l, h0l:h, h00:h, 00l:l  */  //Z=44029
            c1 = -1; c2 = 2; c3 = 1; c4 = -1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44030
            //EditSG.Text = 'Pca21(29),Pcam(57)'; }/*3*/  //Z=44031
            break;
        case 44: /*  Pcaa | 0kl:l, h0l:h, hk0:h, h00:h, 00l:l  */  //Z=44032
            c1 = -1; c2 = 2; c3 = 1; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44033
            //EditSG.Text = 'Pcaa(54)'; }/*3*/  //Z=44034
            break;
        case 45: /*  Pcab | 0kl:l, h0l:h, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44035
            c1 = -1; c2 = 2; c3 = 1; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44036
            //EditSG.Text = 'Pcab(61)'; }/*3*/  //Z=44037
            break;
        case 46: /*  Pcan | 0kl:l, h0l:h, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44038
            c1 = -1; c2 = 2; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44039
            //EditSG.Text = 'Pcan(60)'; }/*3*/  //Z=44040
            break;
        case 47: /*  Pcc- | 0kl:l, h0l:l, 00l:l  */  //Z=44041
            c1 = -1; c2 = 2; c3 = 2; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44042
            //EditSG.Text = 'Pcc2(27),Pccm(49)'; }/*3*/  //Z=44043
            break;
        case 48: /*  Pcca | 0kl:l, h0l:l, hk0:h, h00:h, 00l:l  */  //Z=44044
            c1 = -1; c2 = 2; c3 = 2; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44045
            //EditSG.Text = 'Pcca(54)'; }/*3*/  //Z=44046
            break;
        case 49: /*  Pccb | 0kl:l, h0l:l, hk0:k, 0k0:k, 00l:l  */  //Z=44047
            c1 = -1; c2 = 2; c3 = 2; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44048
            //EditSG.Text = 'Pccb(54)'; }/*3*/  //Z=44049
            break;
        case 50: /*  Pccn | 0kl:l, h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44050
            c1 = -1; c2 = 2; c3 = 2; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44051
            //EditSG.Text = 'Pccn(56)'; }/*3*/  //Z=44052
            break;
        case 51: /*  Pcn- | 0kl:l, h0l:h+l, h00:h, 00l:l  */  //Z=44053
            c1 = -1; c2 = 2; c3 = 0; c4 = -1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44054
            //EditSG.Text = 'Pcn2(30),Pcnm(53)'; }/*3*/  //Z=44055
            break;
        case 52: /*  Pcna | 0kl:l, h0l:h+l, hk0:h, h00:h, 00l:l  */  //Z=44056
            c1 = -1; c2 = 2; c3 = 0; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44057
            //EditSG.Text = 'Pcna(50)'; }/*3*/  //Z=44058
            break;
        case 53: /*  Pcnb | 0kl:l, h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44059
            c1 = -1; c2 = 2; c3 = 0; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44060
            //EditSG.Text = 'Pcnb(60)'; }/*3*/  //Z=44061
            break;
        case 54: /*  Pcnn | 0kl:l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44062
            c1 = -1; c2 = 2; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44063
            //EditSG.Text = 'Pcnn(52)'; }/*3*/  //Z=44064
            break;
        case 55: /*  Pn-- | 0kl:k+l, 0k0:k, 00l:l  */  //Z=44065
            c1 = -1; c2 = 0; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44066
            //EditSG.Text = 'Pnm21(31),Pn21m(31),Pnmm(59)'; }/*3*/  //Z=44067
            break;
        case 56: /*  Pn-a | 0kl:k+l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44068
            c1 = -1; c2 = 0; c3 = -1; c4 = 1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44069
            //EditSG.Text = 'Pn21a(33),Pnma(62)'; }/*3*/  //Z=44070
            break;
        case 57: /*  Pn-b | 0kl:k+l, hk0:k, 0k0:k, 00l:l  */  //Z=44071
            c1 = -1; c2 = 0; c3 = -1; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44072
            //EditSG.Text = 'Pn2b(30),Pnmb(53)'; }/*3*/  //Z=44073
            break;
        case 58: /*  Pn-n | 0kl:k+l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44074
            c1 = -1; c2 = 0; c3 = -1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44075
            //EditSG.Text = 'Pn2n(34),Pnmn(58)'; }/*3*/  //Z=44076
            break;
        case 59: /*  Pna- | 0kl:k+l, h0l:h, h00:h, 0k0:k, 00l:l  */  //Z=44077
            c1 = -1; c2 = 0; c3 = 1; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44078
            //EditSG.Text = 'Pna21(33),Pnam(62)'; }/*3*/  //Z=44079
            break;
        case 60: /*  Pnaa | 0kl:k+l, h0l:h, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44080
            c1 = -1; c2 = 0; c3 = 1; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44081
            //EditSG.Text = 'Pnaa(56)'; }/*3*/  //Z=44082
            break;
        case 61: /*  Pnab | 0kl:k+l, h0l:h, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44083
            c1 = -1; c2 = 0; c3 = 1; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44084
            //EditSG.Text = 'Pnab(60)'; }/*3*/  //Z=44085
            break;
        case 62: /*  Pnan | 0kl:k+l, h0l:h, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44086
            c1 = -1; c2 = 0; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44087
            //EditSG.Text = 'Pnan(52)'; }/*3*/  //Z=44088
            break;
        case 63: /*  Pnc- | 0kl:k+l, h0l:l, 0k0:k, 00l:l  */  //Z=44089
            c1 = -1; c2 = 0; c3 = 2; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44090
            //EditSG.Text = 'Pnc2(30),Pncm(53)'; }/*3*/  //Z=44091
            break;
        case 64: /*  Pnca | 0kl:k+l, h0l:l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44092
            c1 = -1; c2 = 0; c3 = 2; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44093
            //EditSG.Text = 'Pnca(60)'; }/*3*/  //Z=44094
            break;
        case 65: /*  Pncb | 0kl:k+l, h0l:l, hk0:k, 0k0:k, 00l:l |   */  //Z=44095
            c1 = -1; c2 = 0; c3 = 2; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44096
            //EditSG.Text = 'Pncb(50)'; }/*3*/  //Z=44097
            break;
        case 66: /*  Pncn | 0kl:k+l, h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44098
            c1 = -1; c2 = 0; c3 = 2; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44099
            //EditSG.Text = 'Pncn(52)'; }/*3*/  //Z=44100
            break;
        case 67: /*  Pnn- | 0kl:k+l, h0l:h+l, h00:h, 0k0:k, 00l:l  */  //Z=44101
            c1 = -1; c2 = 0; c3 = 0; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44102
            //EditSG.Text = 'Pnn2(34),Pnnm(58)'; }/*3*/  //Z=44103
            break;
        case 68: /*  Pnna | 0kl:k+l, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44104
            c1 = -1; c2 = 0; c3 = 0; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44105
            //EditSG.Text = 'Pnna(52)'; }/*3*/  //Z=44106
            break;
        case 69: /*  Pnnb | 0kl:k+l, h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44107
            c1 = -1; c2 = 0; c3 = 0; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44108
            //EditSG.Text = 'Pnnb(52)'; }/*3*/  //Z=44109
            break;
        case 70: /*  Pnnn | 0kl:k+l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44110
            c1 = -1; c2 = 0; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44111
            //EditSG.Text = 'Pnnn(48)'; }/*3*/  //Z=44112
            break;

            /*  centered C___  */  //Z=44115
        case 71: /*  C--- | hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k  */  //Z=44116
            c1 = 1; c2 = 1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44117
            //EditSG.Text = 'C222(21),Cmm2(35),Cm2m(38),C2mm(38),Cmmm(65)'; }/*3*/  //Z=44118
            break;
        case 72: /*  C--21 | hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k, 00l:.  */  //Z=44119
            c1 = 1; c2 = 1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44120
            //EditSG.Text = 'C2221(20)'; }/*3*/  //Z=44121
            break;
        case 73: /*  C--(ab) | hkl:h+k, 0kl:k, h0l:h, hk0:h,k, h00:h, 0k0:k  */  //Z=44122
            c1 = 1; c2 = 1; c3 = 1; c4 = 4; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44123
            //EditSG.Text = 'Cm2e(39),C2me(39),Cmme(67)'; }/*3*/  //Z=44124
            break;
        case 74: /*  C-c- | hkl:h+k, 0kl:k, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44125
            c1 = 1; c2 = 1; c3 = 4; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44126
            //EditSG.Text = 'Ccm21(36),C2cm(40),Cmcm(63)'; }/*3*/  //Z=44127
            break;
        case 75: /*  C-c(ab) | hkl:h+k, 0kl:k, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44128
            c1 = 1; c2 = 1; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44129
            //EditSG.Text = 'C2ce(41),Cmce(64)'; }/*3*/  //Z=44130
            break;
        case 76: /*  Cc-- | hkl:h+k, 0kl:k,l, h0l:h, hk0:h+k, h00:h, 0k0:k, 001:l  */  //Z=44131
            c1 = 1; c2 = 4; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44132
            //EditSG.Text = 'Ccm21(36),Cc2m(40),Ccmm(63)'; }/*3*/  //Z=44133
            break;
        case 77: /*  Cc-(ab) | hkl:h+k, 0kl:k,l, h0l:h, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44134
            c1 = 1; c2 = 4; c3 = 1; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44135
            //EditSG.Text = 'Cc2e(41),Ccme(64)'; }/*3*/  //Z=44136
            break;
        case 78: /*  Ccc- | hkl:h+k, 0kl:k,l, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44137
            c1 = 1; c2 = 4; c3 = 4; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44138
            //EditSG.Text = 'Ccc2(37),Cccm(66)'; }/*3*/  //Z=44139
            break;
        case 79: /*  Ccc(ab) | hkl:h+k, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44140
            c1 = 1; c2 = 4; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44141
            //EditSG.Text = 'Ccce(68)'; }/*3*/  //Z=44142
            break;

            /*  centered B___  */  //Z=44145
        case 80: /*  B--- | hkl:h+l, 0kl:l, h0l:h+l, hk0:h, h00:h, 00l:l  */  //Z=44146
            c1 = 3; c2 = 2; c3 = 0; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44147
            //EditSG.Text = 'B222(21),Bmm2(38),Bm2m(35),B2mm(38),Bmmm(65)'; }/*3*/  //Z=44148
            break;
        case 81: /*  B-21- | hkl:h+l, 0kl:l, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44149
            c1 = 3; c2 = 2; c3 = 0; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44150
            //EditSG.Text = 'B2212(20)'; }/*3*/  //Z=44151
            break;
        case 82: /*  B--b | hkl:h+l, 0kl:l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44152
            c1 = 3; c2 = 2; c3 = 0; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44153
            //EditSG.Text = 'Bm21b(36),B2mb(40),Bmmb(63)'; }/*3*/  //Z=44154
            break;
        case 83: /*  B-(ac)- | hkl:h+l, 0kl:l, h0l;h,l, hk0:h, h00:h, 00l:l  */  //Z=44155
            c1 = 3; c2 = 2; c3 = 4; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44156
            //EditSG.Text = 'Bme2(39),B2em(39),Bmem(67)'; }/*3*/  //Z=44157
            break;
        case 84: /*  B-(ac)b | hkl:h+l, 0kl:l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44158
            c1 = 3; c2 = 2; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44159
            //EditSG.Text = 'B2eb(41),Bmeb(64)'; }/*3*/  //Z=44160
            break;
        case 85: /*  Bb-- | hkl:h+l, 0kl:k,l, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44161
            c1 = 3; c2 = 4; c3 = 0; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44162
            //EditSG.Text = 'Bbm2(40),Bb21m(36),Bbmm(63)'; }/*3*/  //Z=44163
            break;
        case 86: /*  Bb-b | hkl:h+l, 0kl:k,l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44164
            c1 = 3; c2 = 4; c3 = 0; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44165
            //EditSG.Text = 'Bb2b(37),Bbmb(66)'; }/*3*/  //Z=44166
            break;
        case 87: /*  Bb(ac)- | hkl:h+l, 0kl:k,l, h0l:h,l, hk0:h, h00:h, 0k0:k, 00l:l  */  //Z=44167
            c1 = 3; c2 = 4; c3 = 4; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44168
            //EditSG.Text = 'Bbe2(41),Bbem(64)'; }/*3*/  //Z=44169
            break;
        case 88: /*  Bb(ac)b | hkl:h+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44170
            c1 = 3; c2 = 4; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44171
            //EditSG.Text = 'Bbeb(68)'; }/*3*/  //Z=44172
            break;

            /*  centered A___  */  //Z=44175
        case 89: /*  A--- | hkl:k+l, 0kl:k+l, h0l:l, hk0:k, 0k0:k, 00l:l  */  //Z=44176
            c1 = 2; c2 = 0; c3 = 2; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44177
            //EditSG.Text = 'A222(21),Amm2(38),Am2m(38),A2mm(35),Ammm(65)'; }/*3*/  //Z=44178
            break;
        case 90: /*  A21-- | hkl:k+l, 0kl:k+l, h0l:l, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44179
            c1 = 2; c2 = 0; c3 = 2; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44180
            //EditSG.Text = 'A2122(20)'; }/*3*/  //Z=44181
            break;
        case 91: /*  A--a | hkl:k+l, 0kl:k+l, h0l:l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44182
            c1 = 2; c2 = 0; c3 = 2; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44183
            //EditSG.Text = 'Am2a(40),A21ma(36),Amma(63)'; }/*3*/  //Z=44184
            break;
        case 92: /*  A-a-| hkl:k+l, 0kl:k+l, h0l:h,l, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44185
            c1 = 2; c2 = 0; c3 = 4; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44186
            //EditSG.Text = 'Ama2(40),A21am(36),Amam(63)'; }/*3*/  //Z=44187
            break;
        case 93: /*  A-aa | hkl:k+l, 0kl:k+l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44188
            c1 = 2; c2 = 0; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44189
            //EditSG.Text = 'A2aa(37),Amaa(66)'; }/*3*/  //Z=44190
            break;
        case 94: /*  A(bc)-- | hkl:k+l, 0kl:k,l, h0l:l, hk0:k, 0k0:k, 00l:l  */  //Z=44191
            c1 = 2; c2 = 4; c3 = 2; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44192
            //EditSG.Text = 'Aem2(39),Ae2m(39),Aemm(67)'; }/*3*/  //Z=44193
            break;
        case 95: /*  A(bc)-a | hkl:k+l, 0kl:k,l, h0l:l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44194
            c1 = 2; c2 = 4; c3 = 2; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44195
            //EditSG.Text = 'Ae2a(41),Aema(64)'; }/*3*/  //Z=44196
            break;
        case 96: /*  A(bc)a- | hkl:k+l, 0kl:k,l, h0l:h,l, hk0:k, h00:h, 0k0:k, 00l:l  */  //Z=44197
            c1 = 2; c2 = 4; c3 = 4; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44198
            //EditSG.Text = 'Aea2(41),Aeam(64)'; }/*3*/  //Z=44199
            break;
        case 97: /*  A(bc)aa | hkl:k+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44200
            c1 = 2; c2 = 4; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44201
            //EditSG.Text = 'Aeaa(68)'; }/*3*/  //Z=44202
            break;

            /*  centered I___  */  //Z=44205
        case 98: /*  I--- | hkl:h+k+l, 0kl:k+l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44206
            c1 = 0; c2 = 0; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44207
            //EditSG.Text = 'I222(23),I212121(24),Imm2(44),Im2m(44),I2mm(44),Immm(71)'; }/*3*/  //Z=44208
            break;
        case 99: /*  I--(ab) | hkl:h+k+l, 0kl:k+l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44209
            c1 = 0; c2 = 0; c3 = 0; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44210
            //EditSG.Text = 'Im2a(46),I2mb(46),Imma(74),Immb(74)'; }/*3*/  //Z=44211
            break;
        case 100: /*  I-(ac)- | hkl:h+k+l, 0kl:k+l, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44212
            c1 = 0; c2 = 0; c3 = 4; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44213
            //EditSG.Text = 'Ima2(46),I2cm(46),Imam(74),Imcm(74)'; }/*3*/  //Z=44214
            break;
        case 101: /*  I-cb | hkl:h+k+l, 0kl:k+l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44215
            c1 = 0; c2 = 0; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44216
            //EditSG.Text = 'I2cb(45),Imcb(72)'; }/*3*/  //Z=44217
            break;
        case 102: /*  I(bc)-- | hkl:h+k+l, 0kl:k,l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44218
            c1 = 0; c2 = 4; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44219
            //EditSG.Text = 'Iem2(46),Ie2m(46),Iemm(74)'; }/*3*/  //Z=44220
            break;
        case 103: /*  Ic-a | hkl:h+k+l, 0kl:k,l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44221
            c1 = 0; c2 = 4; c3 = 0; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44222
            //EditSG.Text = 'Ic2a(45),Icma(72)'; }/*3*/  //Z=44223
            break;
        case 104: /*  Iba- | hkl:h+k+l, 0kl:k,l, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l  */  //Z=44224
            c1 = 0; c2 = 4; c3 = 4; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44225
            //EditSG.Text = 'Iba2(45),Ibam(72)'; }/*3*/  //Z=44226
            break;
        case 105: /*  Ibca | hkl:h+k+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44227
            c1 = 0; c2 = 4; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44228
            //EditSG.Text = 'Ibca(73),Icab(73)'; }/*3*/  //Z=44229
            break;

            /*  centered F___  */  //Z=44232
        case 106: /*  F--- | hkl:h+k,h+l,k+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l  */  //Z=44233
            c1 = 6; c2 = 4; c3 = 4; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44234
            //EditSG.Text = 'F222(22),Fmm2(42),Fm2m(42),F2mm(42),Fmmm(69)'; }/*3*/  //Z=44235
            break;
        case 107: /*  F-dd | hkl:h+k,h+l,k+l, 0kl:k,l, h0l:4n,h,l, hk0:4n,h,k, h00:4n, 0k0:4n, 00l:4n  */  //Z=44236
            c1 = 6; c2 = 4; c3 = 3; c4 = 3; c5 = 1; c6 = 1; c7 = 1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44237
            //EditSG.Text = 'F2dd(43)'; }/*3*/  //Z=44238
            break;
        case 108: /*  Fd-d | hkl:h+k,h+l,k+l, 0kl:4n,k,l, h0l:h,l, hk0:4n,h,k, h00:4n, 0k0:4n, 00l:4n  */  //Z=44239
            c1 = 6; c2 = 3; c3 = 4; c4 = 3; c5 = 1; c6 = 1; c7 = 1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44240
            //EditSG.Text = 'Fd2d(43)'; }/*3*/  //Z=44241
            break;
        case 109: /*  Fdd- | hkl:h+k,h+l,k+l, 0kl:4n,k,l, h0l:4n,h,l, hk0:h,k, h00:4n, 0k0:4n, 00l:4n  */  //Z=44242
            c1 = 6; c2 = 3; c3 = 3; c4 = 4; c5 = 1; c6 = 1; c7 = 1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44243
            //EditSG.Text = 'Fdd2(43)'; }/*3*/  //Z=44244
            break;
        case 110: /*  Fddd | hkl:h+k,h+l,k+l, 0kl:4n,k,l, h0l:4n,h,l, hk0:4n,h,k, h00:4n, 0k0:4n, 00l:4n  */  //Z=44245
            c1 = 6; c2 = 3; c3 = 3; c4 = 3; c5 = 1; c6 = 1; c7 = 1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44246
            //EditSG.Text = 'Fddd(70)'; }/*3*/  //Z=44247
            break;
        } // switch ComboBoxOrthorhombic
    } // if RadioButtonSysOrtho


    /* *************** Monoclinic crystal system ************** */  //Z=44252
    if ( RadioButtonSysMonoclinic )
    {
        /* hmin:=0;  //Z=44254 */
        /* kmin:=0;  //Z=44255 */
        /* lmin:=-lmax;  //Z=44256 */
        switch ( ComboBoxMonoclinic )
        {
            /*  P___ unique axis b  */  //Z=44258
        case 0: /*  P1-1  */  //Z=44259
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44260
            //EditSG.Text = 'P121(3),P1m1(6),P12/m1(10);unique b-axis'; }/*3*/  //Z=44261
            //>>>>>> Das hier ist im Screenshot von Prof. Förster zu lesen ....
            break;
        case 1: /*  P1211 | 0k0:k  */  //Z=44262
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44263
            //EditSG.Text = 'P1211(4),P21/m1(11);unique b-axis'; }/*3*/  //Z=44264
            break;
        case 2: /*  P1a1 | h0l,h00:h  */  //Z=44265
            c1 = -1; c2 = -1; c3 = 1; c4 = -1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44266
            //EditSG.Text = 'P1a1(7),P12/a1(13);unique b-axis'; }/*3*/  //Z=44267
            break;
        case 3: /*  P121/a1 | h0l,h00:h, 0k0:k  */  //Z=44268
            c1 = -1; c2 = -1; c3 = 1; c4 = -1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44269
            //EditSG.Text = 'P121/a1(14);unique b-axis'; }/*3*/  //Z=44270
            break;
        case 4: /*  P1c1 | h0l,00l:l  */  //Z=44271
            c1 = -1; c2 = -1; c3 = 2; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44272
            //EditSG.Text = 'P1c1(7),P12/c1(13);unique b-axis'; }/*3*/  //Z=44273
            break;
        case 5: /*  P121/c1 | h0l,00l:l, 0k0:k  */  //Z=44274
            c1 = -1; c2 = -1; c3 = 2; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44275
            //EditSG.Text = 'P121/c1(14);unique b-axis'; }/*3*/  //Z=44276
            break;
        case 6: /*  P1n1 | h0l,h00,00l:h+l  */  //Z=44277
            c1 = -1; c2 = -1; c3 = 0; c4 = -1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44278
            //EditSG.Text = 'P1n1(7),P12/n1(13);unique b-axis'; }/*3*/  //Z=44279
            break;
        case 7: /*  P121/n1 | h0l,h00,00l:h+l, 0k0:k  */  //Z=44280
            c1 = -1; c2 = -1; c3 = 0; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44281
            //EditSG.Text = 'P121/n1(14);unique b-axis'; }/*3*/  //Z=44282
            break;

            /*  C___  */  //Z=44284
        case 8: /*  C1-1 | hkl,0kl,hk0:h+k, h0l,h00,00l:h, 0k0:k  */  //Z=44285
            c1 = 1; c2 = 1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44286
            //EditSG.Text = 'C121(5),C1m1(8),C12/m1(12);unique b-axis'; }/*3*/  //Z=44287
            break;
        case 9: /*  C1c1 | hkl,0kl,hk0:h+k, h0l,h00,00l:h,l, 0k0:k  */  //Z=44288
            c1 = 1; c2 = 1; c3 = 4; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44289
            //EditSG.Text = 'C1c1(9),C12/c1(15);unique b-axis'; }/*3*/  //Z=44290
            break;

            /*  A___  */  //Z=44292
        case 10: /*  A1-1 | hkl,0kl,hk0:k+l, h0l,h00,00l:l, 0k0:k  */  //Z=44293
            c1 = 2; c2 = 0; c3 = 2; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44294
            //EditSG.Text = 'A121(15),A1m1(8),A12/m1(12);unique b-axis'; }/*3*/  //Z=44295
            break;
        case 11: /*  A1n1 | hkl,0kl,hk0:k+l, h0l,h00,00l:h,l, 0k0:k  */  //Z=44296
            c1 = 2; c2 = 0; c3 = 4; c4 = 2; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44297
            //EditSG.Text = 'A1n1(9),A12/n1(15);unique b-axis'; }/*3*/  //Z=44298
            break;

            /*  I___  */  //Z=44300
        case 12: /*  I1-1 | hkl,0kl,hk0:h+k+l, h0l,h00,00l:h+l, 0k0:k  */  //Z=44301
            c1 = 0; c2 = 0; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44302
            //EditSG.Text = 'I121(5),I1m1(8),I12/m1(12)'; }/*3*/  //Z=44303
            break;
        case 13: /*  I1a1 | hkl,0kl,hk0:h+k+l, h0l,h00,00l:h,l, 0k0:k  */  //Z=44304
            c1 = 0; c2 = 0; c3 = 4; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 1;  //Z=44305
            //EditSG.Text = 'I1a1(9),I12/a1(15);unique b-axis'; }/*3*/  //Z=44306
            break;

            /*  P___ unique axis c  */  //Z=44311
        case 14: /*  P11-  */  //Z=44312
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44313
            //EditSG.Text = 'P112(3),P11m(6),P112/m(10);unique c-axis'; }/*3*/  //Z=44314
            break;
        case 15: /*  P1121 | 00l:l  */  //Z=44315
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44316
            //EditSG.Text = 'P1121(4),P1121/m(11);unique c-axis'; }/*3*/  //Z=44317
            break;
        case 16: /*  P11a | hk0,h00,0k0:h  */  //Z=44318
            c1 = -1; c2 = -1; c3 = -1; c4 = 1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44319
            //EditSG.Text = 'P11a(7),P112/a(13);unique c-axis'; }/*3*/  //Z=44320
            break;
        case 17: /*  P1121/a | hk0,h00,0k0:h, 00l:l  */  //Z=44321
            c1 = -1; c2 = -1; c3 = -1; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44322
            //EditSG.Text = 'P1121/a(14);unique c-axis'; }/*3*/  //Z=44323
            break;
        case 18: /*  P11b | hk0,h00,0k0:k  */  //Z=44324
            c1 = -1; c2 = -1; c3 = -1; c4 = 2; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44325
            //EditSG.Text = 'P11b(7),P112/b(13);unique c-axis'; }/*3*/  //Z=44326
            break;
        case 19: /*  P1121/b | hk0,h00,0k0:k, 00l:l  */  //Z=44327
            c1 = -1; c2 = -1; c3 = -1; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44328
            //EditSG.Text = 'P1121/b(14);unique c-axis'; }/*3*/  //Z=44329
            break;
        case 20: /*  P11n | hk0,h00,0k0:h+k  */  //Z=44330
            c1 = -1; c2 = -1; c3 = -1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44331
            //EditSG.Text = 'P11n(7),P112/n(13);unique c-axis'; }/*3*/  //Z=44332
            break;
        case 21: /*  P1121/n | hk0,h00,0k0:h+k, 00l:l  */  //Z=44333
            c1 = -1; c2 = -1; c3 = -1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44334
            //EditSG.Text = 'P1121/n(14);unique c-axis'; }/*3*/  //Z=44335
            break;

            /*  B___  */  //Z=44337
        case 22: /*  B11- | hkl,0kl,h0l:h+l, hk0,h00,0k0:h, 00l:l  */  //Z=44338
            c1 = 3; c2 = 2; c3 = 0; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  c13 = 2;  //Z=44339
            //EditSG.Text = 'B112(5),B11m(8),B112/m(12);unique c-axis'; }/*3*/  //Z=44340
            break;
        case 23: /*  B11n | hkl,0kl,h0l:h+l, hk0,h00,0k0:h,k, 00l:l  */  //Z=44341
            c1 = 3; c2 = 2; c3 = 0; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  c13 = 2;  //Z=44342
            //EditSG.Text = 'B11n(9),B112/n(15);unique c-axis'; }/*3*/  //Z=44343
            break;

            /*  A___  */  //Z=44345
        case 24: /*  A11- | hkl,0kl,h0l:k+l, hk0,h00,0k0:k, 00l:l  */  //Z=44346
            c1 = 2; c2 = 0; c3 = 2; c4 = 2; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44347
            //EditSG.Text = 'A112(15),A11m(8),A112/m(12);unique c-axis'; }/*3*/  //Z=44348
            break;
        case 25: /*  A11a | hkl,0kl,h0l:k+l, hk0,h00,0k0:h,k, 00l:l  */  //Z=44349
            c1 = 2; c2 = 0; c3 = 2; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44350
            //EditSG.Text = 'A11a(9),A112/a(15);unique c-axis'; }/*3*/  //Z=44351
            break;

            /*  I___  */  //Z=44353
        case 26: /*  I11- | hkl,0kl,h0l:h+k+l, hk0,h00,0k0:h+k, 00l:l  */  //Z=44354
            c1 = 0; c2 = 0; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44355
            //EditSG.Text = 'I112(5),I11m(8),I112/m(12);unique c-axis'; }/*3*/  //Z=44356
            break;
        case 27: /*  I11b | hkl,0kl,h0l:h+k+l, hk0,h00,0k0:h,k, 00l:l  */  //Z=44357
            c1 = 0; c2 = 0; c3 = 0; c4 = 4; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 2;  //Z=44358
            //EditSG.Text = 'I11b(9),I112/b(15);unique c-axis'; }/*3*/  //Z=44359
            break;

            /*  P___ unique axis a  */  //Z=44363
        case 28: /*  P-11  */  //Z=44364
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44365
            //EditSG.Text = 'P211(3),Pm11(6),P2/m11(10);unique a-axis'; }/*3*/  //Z=44366
            break;
        case 29: /*  P2111 | h00:h  */  //Z=44367
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = 0; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44368
            //EditSG.Text = 'P2111(4),P21/m11(11);unique a-axis'; }/*3*/  //Z=44369
            break;
        case 30: /*  Pb11 | 0kl,0k0,00l:k  */  //Z=44370
            c1 = -1; c2 = 1; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44371
            //EditSG.Text = 'Pb11(7),P2/b11(13);unique a-axis'; }/*3*/  //Z=44372
            break;
        case 31: /*  P21/b11 | 0kl,0k0,00l:k, h00:h  */  //Z=44373
            c1 = -1; c2 = 1; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44374
            //EditSG.Text = 'P21/b11(14);unique q-axis'; }/*3*/  //Z=44375
            break;
        case 32: /*  Pc11 | 0kl,0k0,00l:l  */  //Z=44376
            c1 = -1; c2 = 2; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44377
            //EditSG.Text = 'Pc11(7),P2/c11(13);unique a-axis'; }/*3*/  //Z=44378
            break;
        case 33: /*  P21/c11 | 0kl,0k0,00l:l, h00:h  */  //Z=44379
            c1 = -1; c2 = 2; c3 = -1; c4 = -1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44380
            //EditSG.Text = 'P21/c11(14);unique a-axis'; }/*3*/  //Z=44381
            break;
        case 34: /*  Pn11 | 0kl,0k0,00l:k+l  */  //Z=44382
            c1 = -1; c2 = 0; c3 = -1; c4 = -1; c5 = -1; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44383
            //EditSG.Text = 'Pn11(7),P2/n11(13);unique a-axis'; }/*3*/  //Z=44384
            break;
        case 35: /*  P21/n11 | 0kl,0k0,00l:k+l, h00:h  */  //Z=44385
            c1 = -1; c2 = 0; c3 = -1; c4 = -1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44386
            //EditSG.Text = 'P21/n11(14);unique c-axis'; }/*3*/  //Z=44387
            break;

            /*  C___  */  //Z=44389
        case 36: /*  C-11 | hkl,h0l,hk0:h+k, 0kl,0k0,00l:k, h00:h  */  //Z=44390
            c1 = 1; c2 = 1; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44391
            //EditSG.Text = 'C211(5),Cm11(8),C2/m11(12);unique a-axis'; }/*3*/  //Z=44392
            break;
        case 37: /*  Cn11 | hkl,h0l,hk0:h+k, 0kl,0k0,00l:k,l, h00:h  */  //Z=44393
            c1 = 1; c2 = 4; c3 = 1; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44394
            //EditSG.Text = 'Cn11(9),C2/n11(15);unique a-axis'; }/*3*/  //Z=44395
            break;

            /*  B___  */  //Z=44397
        case 38: /*  B-11 | hkl,h0l,hk0:h+l, 0kl,0k0,00l:l, h00:h  */  //Z=44398
            c1 = 3; c2 = 2; c3 = 0; c4 = 1; c5 = 0; c6 = -1; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44399
            //EditSG.Text = 'B211(15),Bm11(8),B2/m11(12);unique a-axis'; }/*3*/  //Z=44400
            break;
        case 40: /*  Bb11 | hkl,h0l,hk0:h+l, 0kl,0k0,00l:k,l, h00:h  */  //Z=44401
            c1 = 3; c2 = 4; c3 = 0; c4 = 1; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44402
            //EditSG.Text = 'Bb11(9),B2/b11(15);unique a-axis'; }/*3*/  //Z=44403
            break;

            /*  I___  */  //Z=44405
        case 41: /*  I-11 | hkl,h0l,hk0:h+k+l, 0kl,0k0,00l:k+l, h00:h  */  //Z=44406
            c1 = 0; c2 = 0; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44407
            //EditSG.Text = 'I211(5),Im11(8),I2/m11(12);unique a-axis'; }/*3*/  //Z=44408
            break;
        case 42: /*  Ic11 | hkl,h0l,hk0:h+k+l, 0kl,0k0,00l:k,l, h00:h  */  //Z=44409
            c1 = 0; c2 = 0; c3 = 0; c4 = 0; c5 = 0; c6 = 0; c7 = 0; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1; c13 = 0;  //Z=44410
            //EditSG.Text = 'Ic11(9),I2/c11(15);unique a-axis'; }/*3*/  //Z=44411
            break;
        } // switch ComboBoxMonoclinic
    } // if RadioButtonSysMonoclinic


    /* *************** Triclinic crystal system ************** done  */  //Z=44415
    if ( RadioButtonSysTriclinic )
    {
        hmin = -hmax;  //Z=44417
        kmin = -kmax;  //Z=44418
        lmin = 0;  //Z=44419
        switch ( ComboBoxTriclinic )
        {
        /*  P___  */  //Z=44420
        case 0: /*  P- |  */  //Z=44421
            c1 = -1; c2 = -1; c3 = -1; c4 = -1; c5 = -1; c6 = -1; c7 = -1; c8 = -1; c9 = -1; c10 = -1; c11 = -1; c12 = -1;  //Z=44422
            //EditSG.Text = 'P1(1),P-1(2)'; }/*3*/  //Z=44423
            break;
        } // switch ComboBoxTriclinic
    } // if RadioButtonSysTriclinic


    //(* update unit cell values *)         //{NV}-31652
    a=params.uca; // StrToFloat(EditCellA.Text);
    b=params.ucb; // StrToFloat(EditCellB.Text);
    c=params.ucc; // StrToFloat(EditCellC.Text);
    if ( a<=0 ) a=1e-10;
    if ( b<=0 ) b=1e-10;
    if ( c<=0 ) c=1e-10;
    alf=params.ucalpha_deg * M_PI/180.;  // StrToFloat(EditCellAlpha.Text);
    gam=params.ucgamma_deg * M_PI/180.;  // StrToFloat(EditCellGamma.Text);

#ifdef UseStringGrid12
    bet=params.ucbeta_deg  * M_PI/180.;  // StrToFloat(EditCellBeta.Text);
    amax=params.amax; // StrToFloat(EditAmax.Text);
    bmax=params.bmax; // StrToFloat(EditBmax.Text);
    cmax=params.cmax; // StrToFloat(EditCmax.Text);

    wave = params.length;       //{NV}-30753, ist aber dort auskommentiert
    // Ansonsten wird <wave> nirgendwo gesetzt, aber genutzt....
#endif


    //(* unit cell coordinates *) -31671
    if ( RadioButtonSys1D )
    {
        c=a*100;
        alf=M_PI/2.;
#ifdef UseStringGrid12
        bet=M_PI/2.;
#endif
        kmin=0;
        kmax=0;
        lmin=0;
        lmax=0;
    }
    if ( RadioButtonSys2D )
    {
        c=a*100;
        alf=M_PI/2.;
#ifdef UseStringGrid12
        bet=M_PI/2.;
#endif
        lmin=0;
        lmax=0;
    }
    if ( RadioButtonSysCubic ) // -31686
    {
        b=a;  c=a;    alf=M_PI/2.;  gam=M_PI/2.;    //{NV}-31687
#ifdef UseStringGrid12
        bet=M_PI/2.;  bmax=amax;    cmax=amax;
#endif \
    //EditCellb.Text:=FloatToStr(b);  EditCellc.Text:=FloatToStr(c);
    }
    if ( RadioButtonSysTetragonal )
    {
        b=a;  alf=M_PI/2.0;  gam=M_PI/2.0;
#ifdef UseStringGrid12
        bet=M_PI/2.0;  bmax=amax;
#endif \
    //EditCellb.Text:=FloatToStr(b);
    }
    if ( RadioButtonSysOrtho )
    {
        alf=M_PI/2.;  gam=M_PI/2.;
#ifdef UseStringGrid12
        bet=M_PI/2.;
#endif
    }
    if ( RadioButtonSysHex )
    {
        b=a; alf=M_PI/2.;  gam=2*M_PI/3.;
#ifdef UseStringGrid12
        bet=M_PI/2.;  bmax=amax;
#endif \
    //EditCellb.Text:=FloatToStr(b);
    }
    if ( RadioButtonSysTrigonal )   //(* reads angle alfa *)
    {
        if ( trigP=='h' )  //(* hexagonal axes *)
        {
            b=a;
#ifdef UseStringGrid12
            bmax=amax;
            bet=M_PI/2.;
#endif
            alf=M_PI/2.;
            gam=2*M_PI/3.;
            //EditCellA.Text:=FloatToStr(a);
            //EditCellB.Text:=FloatToStr(b);
            //EditCellC.Text:=FloatToStr(c);
            //EditCellAlpha.Text:=FloatToStr(alf*180/pi);
            //EditCellBeta.Text:=FloatToStr(bet*180/pi);
            //EditCellGamma.Text:=FloatToStr(gam*180/pi);
            //TrackBarCellC.Visible:=true;
            //TrackBarCellAlpha.Visible:=false;
        }
        if ( trigP=='r' )  //(* rhomobohedral axes *)
        {
            b=a;
            c=a;
#ifdef UseStringGrid12
            bmax=amax;
            cmax=amax;
            bet=alf;
#endif
            gam=alf;
            //EditCellA.Text:=FloatToStr(a);
            //EditCellB.Text:=FloatToStr(b);
            //EditCellC.Text:=FloatToStr(c);
            //EditCellAlpha.Text:=FloatToStr(alf*180/pi);
            //EditCellBeta.Text:=FloatToStr(bet*180/pi);
            //EditCellGamma.Text:=FloatToStr(gam*180/pi);
            //TrackBarCellC.Visible:=false;
            //TrackBarCellAlpha.Visible:=true;
        }
    }
    if ( RadioButtonSysMonoclinic )   //(* reads angle gamma *)
    {
        switch ( c13 )
        {
        case 0:  //(* unique axis a *)
            alf=gam;      gam=M_PI/2.;
#ifdef UseStringGrid12
            bet=M_PI/2.;
#endif
            break;
        case 1:  //(* unique axis b *)
            alf=M_PI/2.;    gam=M_PI/2.;
#ifdef UseStringGrid12
            bet=gam;
#endif
            break;
        case 2:  //(* unique axis c *)
            alf=M_PI/2.;
#ifdef UseStringGrid12
            bet=M_PI/2.;
#endif
            break;
        } // switch c13
    }

#ifdef UseStringGrid12
    s2a=sin(alf)*sin(alf);
    c1a=cos(alf);
    c2a=c1a*c1a;
    c3a=c1a*c1a*c1a;
    c1b=cos(bet);
    c2b=c1b*c1b;
    s2b=sin(bet)*sin(bet);
    c1g=cos(gam);
    c2g=c1g*c1g;
    s2g=sin(gam)*sin(gam);

//if ((amax>=bmax) && (amax>=cmax)) mx=amax;      //{NV}-31760
//if ((bmax>=amax) && (bmax>=cmax)) mx=bmax;
//if ((cmax>=amax) && (cmax>=bmax)) mx=cmax;
#endif

    //{NV}-31907 bis 32410 ist nur für Image38

    //(************ generate hkl-Table *************)       //{NV}-32417
    //(* loop 1 for isotropic, loop 2 for anisotropic *)
    for ( ii=1; ii<=2; ii++ )
    {
#ifdef UseStringGrid12
        ct=1;   // Zähler der Schleife für StrinGrid12[*][ct]
#endif
        ct1=1;  // laufender Index in latpar1/2[ct1][*]
        if ( ii==2 )
        {
            hmin=-hmax;
            kmin=-hmax;
            lmin=-hmax;
        }

        for ( h=hmin; h<=hmax; h++ )
        {
            for ( k=kmin; k<=kmax; k++ )
            {
                for ( l=lmin; l<=lmax; l++ )
                {
#ifdef UseStringGrid12
                    StringGrid12[0] = QString::number(h);   //StringGrid12.Cells[0,ct]:=IntToStr(h);
                    StringGrid12[1] = QString::number(k);   //StringGrid12.Cells[1,ct]:=IntToStr(k);
                    StringGrid12[2] = QString::number(l);   //StringGrid12.Cells[2,ct]:=IntToStr(l);

                    //(* calculate q-values *)
                    if ( RadioButtonSys1D )
                        invd=(h*h/(a*a));
                    else if ( RadioButtonSys2D )
                        invd=(h*h*b*b+2.*h*k*a*b*cos(gam)+k*k*a*a)/(a*a*b*b*sin(gam)*sin(gam));
                    else if ( RadioButtonSysCubic )
                        invd=(h*h+k*k+l*l)/(a*a);               //{NV}-32442
                    else if ( RadioButtonSysTetragonal )
                        invd=(h*h+k*k)/(a*a)+l*l/(c*c);
                    else if ( RadioButtonSysOrtho )
                        invd=h*h/(a*a)+k*k/(b*b)+l*l/(c*c);
                    else if ( RadioButtonSysHex )
                        invd=4.*(h*h+k*k+h*k)/(3.*a*a)+l*l/(c*c);
                    else if ( RadioButtonSysTrigonal && trigP=='h' )    //(* hexagonal axes *)
                        invd=4.*(h*h+k*k+h*k)/(3.*a*a)+l*l/(c*c);
                    else if ( RadioButtonSysTrigonal && trigP=='r' )    //(* rhombohedral axes *)
                        invd=(1./(a*a))*((h*h+k*k+l*l)*s2a+2.*(h*k+h*l+k*l)*(c2a-c1a))/(1.+2*c3a-3*c2a);
                    else if ( RadioButtonSysMonoclinic )
                        invd=h*h/(a*a*s2b)+k*k/(b*b)+l*l/(c*c*s2b)-2.*h*l*c1b/(a*c*s2b);
                    else if ( RadioButtonSysTriclinic )
                        invd=(h*h*s2a/(a*a)+k*k*s2b/(b*b)+l*l*s2g/(c*c)+2.*k*l*(c1b*c1g-c1a)/(b*c)+2.*l*h*(c1g*c1a-c1b)/(c*a)+2.*h*k*(c1a*c1b-c1g)/(a*b))/(1.-c2a-c2b-c2g+2*c1a*c1b*c1g);
                    else
                        invd = 0;

                    if ( invd>0 )
                    {
                        q=2*M_PI*sqrt(invd);
                        dhkl=2*M_PI/q;
                    }
                    else
                    {
                        q=0.0;
                        dhkl=0.0;
                    }
                    //index1[ct]=q;  diese Arrays werden sonst nirgendwo verwendet ...
                    //index2[ct]=ct;
                    StringGrid12[7] = QString::number(q)+" "+QString::number(invd);
                    //StringGrid12.Cells[7,ct]:=FloatToStrF(q,FFFixed,5,3);
                    sinarg=q*wave/(4.*M_PI);
                    if ( sinarg>1 ) sinarg=0.0;
                    ttheta=(180./M_PI)*2*asin(sinarg);
                    StringGrid12[8] = QString::number(ttheta);   //StringGrid12.Cells[8,ct]:=FloatToStrF(ttheta,FFFixed,5,3);
#endif \
    //(* sort into categories for multiplicity *)
                    hklcase=_hkl; // Default

                    if ((h==0) && (k==0) && (l==0)) hklcase=_000;

                    if ((h==0) && (k!=0) && (l!=0) && (k!=l)) hklcase=_0kl;       //(* one index = zero AND two indices not equal *)
                    if ((k==0) && (h!=0) && (l!=0) && (h!=l)) hklcase=_h0l;
                    if ((l==0) && (h!=0) && (k!=0) && (h!=k)) hklcase=_hk0;

                    if ((h!=0) && (k==0) && (l==0)) hklcase=_h00;        //(* two indices = zero *)
                    if ((k!=0) && (h==0) && (l==0)) hklcase=_0k0;
                    if ((l!=0) && (h==0) && (k==0)) hklcase=_00l;

                    if ((h!=0) && (h==k) && (l!=0) && (l!=h)) hklcase=_hhl;       //(* two indices are equal AND the other is not zero *)
                    if ((h!=0) && (h==l) && (k!=0) && (k!=h)) hklcase=_hlh;       //(* 'hlh' *)
                    if ((k!=0) && (k==l) && (h!=0) && (h!=k)) hklcase=_lhh;       //(* 'lhh' *)

                    if ((h!=0) && (h==k) && (l==0)) hklcase=_hh0;        //(* one index = zero   AND   the other two are equal *)
                    if ((h!=0) && (h==l) && (k==0)) hklcase=_h0h;        //(* 'h0l' *)
                    if ((k!=0) && (k==l) && (h==0)) hklcase=_0hh;        //(* '0hh' *)

                    if ((h!=0) && (h==k) && (h==l)) hklcase=_hhh;        //(* all indices are equal *)

//if ( hklcase=="" ) hklcase="hkl";      //(* most general *)
//{NV}-32497

#ifdef UseStringGrid12
                    StringGrid12[3] = hklcase2str[hklcase];      //StringGrid12.Cells[3,ct]:=hklcase;
                    StringGrid12[5] = QString::number(dhkl);     //StringGrid12.Cells[5,ct]:=FloatToStrF(dhkl,FFFixed,5,3);
#endif
                    mult = 0;   // Default, wird aber dann ausgeblendet

                    //(* Multiplicities *)
                    if ( RadioButtonSys1D || RadioButtonSys2D )
                        mult=1;

                    else if ( RadioButtonSysCubic )    //{NV}-32507
                    {
                        if ( hklcase==_hkl ) mult=48;
                        if ( hklcase==_hhl ) mult=24;
                        if ( hklcase==_hlh ) mult=24;
                        if ( hklcase==_lhh ) mult=24;
                        if ( hklcase==_hh0 ) mult=12;
                        if ( hklcase==_h0h ) mult=12;
                        if ( hklcase==_0hh ) mult=12;
                        if ( hklcase==_hhh ) mult=8;
                        if ( hklcase==_hk0 ) mult=24;
                        if ( hklcase==_h0l ) mult=24;
                        if ( hklcase==_0kl ) mult=24;
                        if ( hklcase==_h00 ) mult=6;
                        if ( hklcase==_0k0 ) mult=6;
                        if ( hklcase==_00l ) mult=6;
                        //multmax=48;
                    }

                    else if ( RadioButtonSysTetragonal )
                    {
                        if ( hklcase==_hkl ) mult=16;
                        if ( hklcase==_hhl ) mult=8;
                        if ( hklcase==_hlh ) mult=8;
                        if ( hklcase==_lhh ) mult=8;
                        if ( hklcase==_hh0 ) mult=4;
                        if ( hklcase==_h0h ) mult=8;
                        if ( hklcase==_0hh ) mult=8;
                        if ( hklcase==_hhh ) mult=8;
                        if ( hklcase==_hk0 ) mult=8;
                        if ( hklcase==_h0l ) mult=8;
                        if ( hklcase==_0kl ) mult=8;
                        if ( hklcase==_h00 ) mult=4;
                        if ( hklcase==_0k0 ) mult=4;
                        if ( hklcase==_00l ) mult=2;
                        //multmax=16;
                    }

                    else if ( RadioButtonSysHex )
                    {
                        if ( hklcase==_hkl ) mult=24;
                        if ( hklcase==_hhl ) mult=12;
                        if ( hklcase==_hlh ) mult=12;
                        if ( hklcase==_lhh ) mult=12;
                        if ( hklcase==_hh0 ) mult=6;
                        if ( hklcase==_h0h ) mult=12;
                        if ( hklcase==_0hh ) mult=12;
                        if ( hklcase==_hhh ) mult=12;
                        if ( hklcase==_hk0 ) mult=12;
                        if ( hklcase==_h0l ) mult=12;
                        if ( hklcase==_0kl ) mult=12;
                        if ( hklcase==_h00 ) mult=6;
                        if ( hklcase==_0k0 ) mult=6;
                        if ( hklcase==_00l ) mult=2;
                        //multmax=24;
                    }

                    else if ( RadioButtonSysOrtho )
                    {
                        if ( hklcase==_hkl ) mult=8;
                        if ( hklcase==_hhl ) mult=8;
                        if ( hklcase==_hlh ) mult=8;
                        if ( hklcase==_lhh ) mult=8;
                        if ( hklcase==_hh0 ) mult=8;
                        if ( hklcase==_h0h ) mult=8;
                        if ( hklcase==_0hh ) mult=8;
                        if ( hklcase==_hhh ) mult=8;
                        if ( hklcase==_hk0 ) mult=4;
                        if ( hklcase==_h0l ) mult=4;
                        if ( hklcase==_0kl ) mult=4;
                        if ( hklcase==_h00 ) mult=2;
                        if ( hklcase==_0k0 ) mult=2;
                        if ( hklcase==_00l ) mult=2;
                        //multmax:=8;
                    }

                    else if ( RadioButtonSysMonoclinic )
                    {
                        if ( hklcase==_hkl ) mult=4;
                        if ( hklcase==_hhl ) mult=4;
                        if ( hklcase==_hlh ) mult=4;
                        if ( hklcase==_lhh ) mult=4;
                        if ( hklcase==_hh0 ) mult=4;
                        if ( hklcase==_h0h ) mult=4;
                        if ( hklcase==_0hh ) mult=4;
                        if ( hklcase==_hhh ) mult=4;
                        if ( hklcase==_hk0 ) mult=4;
                        if ( hklcase==_h0l ) mult=2;
                        if ( hklcase==_0kl ) mult=4;
                        if ( hklcase==_h00 ) mult=2;
                        if ( hklcase==_0k0 ) mult=2;
                        if ( hklcase==_00l ) mult=2;
                        //multmax:=4;
                    }

                    else if ( RadioButtonSysTriclinic || RadioButtonSysTrigonal )
                    {
                        if ( hklcase==_hkl ) mult=2;
                        if ( hklcase==_hhl ) mult=2;
                        if ( hklcase==_hlh ) mult=2;
                        if ( hklcase==_lhh ) mult=2;
                        if ( hklcase==_hh0 ) mult=2;
                        if ( hklcase==_h0h ) mult=2;
                        if ( hklcase==_0hh ) mult=2;
                        if ( hklcase==_hhh ) mult=2;
                        if ( hklcase==_hk0 ) mult=2;
                        if ( hklcase==_h0l ) mult=2;
                        if ( hklcase==_0kl ) mult=2;
                        if ( hklcase==_h00 ) mult=2;
                        if ( hklcase==_0k0 ) mult=2;
                        if ( hklcase==_00l ) mult=2;
                        //multmax:=2;
                    }

#ifdef UseStringGrid12
                    StringGrid12[6] = QString::number(mult);    //StringGrid12.Cells[6,ct]:=IntToStr(mult);
#endif

                    //(* go through extinction conditions *)    //{NV}-32618
                    //(* 1D and 2D systems *)
                    if ( RadioButtonSys1D || RadioButtonSys2D )
                    {
                        ccase=1;
                    }

                    //(* cubic system *)
                    if ( RadioButtonSysCubic )
                    {
                        c1case=0;      //(* hkl *)
                        c2case=1;      //(* 0kl *)
                        c7case=1;      //(* 00l *)
                        c10case=1;     //(* hhl *)

                        //(* integral condition *)
                        if ( c1==-1 ) c1case=1;
                        if ( c1== 0 ) if (((h+k+l) % 2)==0) c1case=1;   //(* F *)
                        if ( c1== 6 ) if ((((h+k) % 2)==0) && (((h+l) % 2)==0) && (((k+l) % 2)==0)) c1case=1;  //(* I *)

                        //(* zonal conditions *)
                        //(* 0kl:k+l *)
                        if ((c2==0) && (h==0) && ! (((k+l) % 2)==0)) c2case=0;  //(* 0kl:k+l *)
                        if ((c2==0) && (k==0) && ! (((h+l) % 2)==0)) c2case=0;  //(* h0l:h+l *)
                        if ((c2==0) && (l==0) && ! (((h+k) % 2)==0)) c2case=0;  //(* hk0:h+k *)
                        //(* 0kl:k+l=4n,k,l *)
                        if ((c2==3) && (h==0) && ! ((((k+l) % 4)==0) && ((k % 2)==0) && ((l % 2)==0))) c2case=0; //(* 0kl *)
                        if ((c2==3) && (k==0) && ! ((((h+l) % 4)==0) && ((h % 2)==0) && ((l % 2)==0))) c2case=0; //(* h0l *)
                        if ((c2==3) && (l==0) && ! ((((h+k) % 4)==0) && ((h % 2)==0) && ((k % 2)==0))) c2case=0; //(* hk0 *)
                        //(* 0kl:k,l *)
                        if ((c2==4) && (h==0) && ! (((k % 2)==0) && ((l % 2)==0))) c2case=0; //(* 0kl:k,l *)
                        if ((c2==4) && (k==0) && ! (((h % 2)==0) && ((l % 2)==0))) c2case=0; //(* h0l:h,l *)
                        if ((c2==4) && (l==0) && ! (((h % 2)==0) && ((k % 2)==0))) c2case=0; //(* hk0:h,k *)
                        //(* 0kl:k and cyclic permutations *)
                        if ((c2==5) && (h==0) && ! ((k % 2)==0)) c2case=0; //(* 0kl:k *)
                        if ((c2==5) && (k==0) && ! ((l % 2)==0)) c2case=0; //(* h0l:l *)
                        if ((c2==5) && (l==0) && ! ((h % 2)==0)) c2case=0; //(* hk0:h *)

                        //(* hhl:l *)
                        if ((c10==0) && (h==k) && (l!=0) && ! ((l % 2)==0)) c10case=0;  //(* hhl:l *)
                        if ((c10==0) && (h==l) && (k!=0) && ! ((k % 2)==0)) c10case=0;  //(* hkh:k *)
                        if ((c10==0) && (k==l) && (h!=0) && ! ((h % 2)==0)) c10case=0;  //(* hkk:h *)
                        //(* hhl:h+l *)
                        if ((c10==2) && (h==k) && (l!=0) && ! (((h+l) % 2)==0)) c10case=0;  //(* hhl:h+l *)
                        if ((c10==2) && (h==l) && (k!=0) && ! (((h+k) % 2)==0)) c10case=0;  //(* hkh:h+k *)
                        if ((c10==2) && (k==l) && (h!=0) && ! (((h+k) % 2)==0)) c10case=0;  //(* hkk:h+k *)
                        //(* hhl:2h+l=4n,l *)
                        if ((c10==4) && (h==k) && (l!=0) && ! ((((2*h+l) % 4)==0) && ((l % 2)==0))) c10case=0; //(* hhl *)
                        if ((c10==4) && (h==l) && (k!=0) && ! ((((2*h+k) % 4)==0) && ((k % 2)==0))) c10case=0; //(* hkh *)
                        if ((c10==4) && (k==l) && (h!=0) && ! ((((2*k+h) % 4)==0) && ((h % 2)==0))) c10case=0; //(* hkk *)
                        //(* hhl:h,l *)
                        if ((c10==5) && (h==k) && (l!=0) && ! (((h % 2)==0) && ((l % 2)==0))) c10case=0; //(* hhl:h,l *)
                        if ((c10==5) && (h==l) && (k!=0) && ! (((h % 2)==0) && ((k % 2)==0))) c10case=0; //(* hkh:h,k *)
                        if ((c10==5) && (k==l) && (h!=0) && ! (((k % 2)==0) && ((h % 2)==0))) c10case=0; //(* hkk:h,k *)

                        //(* serial conditions *)
                        //(* 00l:l *)
                        if ((c7==0) && (h==0) && (k==0) && ! ((l % 2)==0)) c7case=0; //(* 00l:0 *)
                        if ((c7==0) && (h==0) && (l==0) && ! ((k % 2)==0)) c7case=0; //(* 0k0:k *)
                        if ((c7==0) && (k==0) && (l==0) && ! ((h % 2)==0)) c7case=0; //(* h00:h *)
                        //(* 00l:4n *)
                        if ((c7==1) && (h==0) && (k==0) && ! ((l % 4)==0)) c7case=0; //(* 00l:4n *)
                        if ((c7==1) && (h==0) && (l==0) && ! ((k % 4)==0)) c7case=0; //(* 0k0:4n *)
                        if ((c7==1) && (k==0) && (l==0) && ! ((h % 4)==0)) c7case=0; //(* h00:4n *)

                        ccase=c1case*c2case*c7case*c10case;
                    }

                    //(* tetragonal system *)
                    if ( RadioButtonSysTetragonal )
                    {
                        c1case = 0;      /* hkl */
                        c2case = 1;      /* 0kl, h0l */
                        c4case = 1;      /* hk0 */
                        c6case = 1;      /* h00, 0k0 */
                        c7case = 1;      /* 00l */
                        c8case = 1;      /* hh0 */
                        c10case = 1;     /* hhl */

                        /* integral condition */
                        if ( c1==-1 ) c1case = 1;
                        if ( c1==0 ) if ( (h+k+l) % 2==0 ) c1case = 1;     /* I */

                        /* zonal conditions */
                        /* 0kl:k+l */
                        if ( (c2==0) && (h==0) && ! ((k+l) % 2==0) ) c2case = 0; /* 0kl:k+l */
                        if ( (c2==0) && (k==0) && ! ((h+l) % 2==0) ) c2case = 0; /* h0l:h+l */
                        /* 0kl:k */
                        if ( (c2==1) && (h==0) && ! (k % 2==0) ) c2case = 0; /* 0kl:k */
                        if ( (c2==1) && (k==0) && ! (h % 2==0) ) c2case = 0; /* h0l:h */
                        /* 0kl:l */
                        if ( (c2==2) && (h==0) && ! (l % 2==0) ) c2case = 0; /* 0kl:l */
                        if ( (c2==2) && (k==0) && ! (l % 2==0) ) c2case = 0; /* h0l:l */
                        /* 0kl:k,l */
                        if ( (c2==4) && (h==0) && ! ((k % 2==0) && (l % 2==0)) ) c2case = 0; /* 0kl:k,l */
                        if ( (c2==4) && (k==0) && ! ((h % 2==0) && (l % 2==0)) ) c2case = 0; /* h0l:h,l */

                        /* hk0:h+k */
                        if ( (c4==0) && (l==0) && ! ((h+k) % 2==0) ) c4case = 0; /* hk0:h+k */
                        /* hk0:h,k */
                        if ( (c4==4) && (l==0) && ! ((h % 2==0) && (k % 2==0)) ) c4case = 0;  /* hk0:h,k */

                        /* hhl:l */
                        if ( (c10==0) && (h==k) && (l!=0) && ! (l % 2==0) ) c10case = 0; /* hhl:l */

                        /* hhl:2h+l=4n,l */
                        if ( (c10==4) && (h==k) && (l!=0) && ! (((2*h+l) % 4==0) && (l % 2==0)) ) c10case = 0; /* hhl:2h+l=4n,l */

                        /* serial conditions */
                        /* 0k0:k */
                        if ( (c6==0) && (h==0) && (l==0) && ! (k % 2==0) ) c6case = 0; /* 0k0:k */
                        if ( (c6==0) && (k==0) && (l==0) && ! (h % 2==0) ) c6case = 0; /* h00:h */

                        /* 00l:l */
                        if ( (c7==0) && (h==0) && (k==0) && ! (l % 2==0) ) c7case = 0; /* 00l:l */
                        /* 00l:4n */
                        if ( (c7==1) && (h==0) && (k==0) && ! (l % 4==0) ) c7case = 0; /* 00l:4n */

                        /* hh0:h */
                        if ( (c8==0) && (h==k) && (l==0) && ! (h % 2==0) ) c8case = 1; /* hh0:h */

                        ccase = c1case*c2case*c4case*c6case*c7case*c8case*c10case;
                    }

#ifndef __CUDACC__
                    if ( RadioButtonSysHex ) qDebug() << "RadioButtonSysHex";
#endif
                        /*(*** Hexagonal system ***)
                    if RadioButtonSysHex.Checked=true then begin
                       ccase:=0;
                       c7case:=1;      (* 00l *)
                       c9case:=1;      (* h-hl *)
                       c10case:=1;     (* hhl *)

                       (* zonal conditions *)
                       (* h-hl:l, only for h0l-case *)
                       if ((c9=0) and (h<>0) and (k=-h) and not (l mod 2=0)) then c9case:=0; (* h-hl:l *)

                       (* hhl:l *)
                       if ((c10=0) and (h=k) and (l<>0) and not (l mod 2=0)) then c10case:=0; (* hhl:l *)

                       (* serial conditions *)
                       (* 00l:l *)
                       if ((c7=0) and (h=0) and (k=0) and not (l mod 2=0)) then c7case:=0; (* 00l:l *)
                       (* 00l:3n *)
                       if ((c7=2) and (h=0) and (k=0) and not (l mod 3=0)) then c7case:=0; (* 00l:3n *)
                       (* 00l:6n *)
                       if ((c7=3) and (h=0) and (k=0) and not (l mod 6=0)) then c7case:=0; (* 00l:6n *)

                       ccase:=c7case*c9case*c10case;
                    end;*/

#ifndef __CUDACC__
                    if ( RadioButtonSysOrtho ) qDebug() << "RadioButtonSysOrtho";
#endif
                        /*(*** Orthorhombic system ***)
                    if ((RadioButtonSysOrtho.Checked=true) or (RadioButtonSys2D.Checked=true)) then begin
                       ccase:=0;
                       c1case:=0;      (* hkl *)
                       c2case:=1;      (* 0kl *)
                       c3case:=1;      (* h0l *)
                       c4case:=1;      (* hk0 *)
                       c5case:=1;      (* h00 *)
                       c6case:=1;      (* 0k0 *)
                       c7case:=1;      (* 00l *)

                       (* integral condition *)
                       if c1=-1 then c1case:=1; (* hkl *)
                       if c1=0 then if ((h+k+l) mod 2=0) then c1case:=1; (* h+k+l *)   (* I *)
                       if c1=1 then if ((h+k) mod 2=0) then c1case:=1;   (* h+k *)     (* C *)
                       if c1=2 then if ((k+l) mod 2=0) then c1case:=1;   (* k+l *)     (* A *)
                       if c1=3 then if ((h+l) mod 2=0) then c1case:=1;   (* h+l *)     (* B *)
                       if c1=6 then if (   ((h+k) mod 2=0) and ((h+l) mod 2=0) and ((k+l) mod 2=0)) then c1case:=1;  (* F *)

                       (* zonal conditions *)
                       (* 0kl *)
                       if ((c2=0) and (h=0) and not ((k+l) mod 2=0)) then c2case:=0;  (* 0kl:k+l *) (* glide plane n-- *)
                       if ((c2=1) and (h=0) and not (k mod 2=0)) then c2case:=0;  (* 0kl:k *)  (* glide plane b-- *)
                       if ((c2=2) and (h=0) and not (l mod 2=0)) then c2case:=0;  (* 0kl:l *)  (* glide plane c-- *)
                       if ((c2=3) and (h=0) and not (((k+l) mod 4=0) and (k mod 2=0) and (l mod 2=0))) then c2case:=0; (* 0kl:k+l=4n,k,l *)  (* glide plane d-- *)
                       if ((c2=4) and (h=0) and not ((k mod 2=0) and (l mod 2=0))) then c2case:=0;

                       (* h0l *)
                       if ((c3=0) and (k=0) and not ((h+l) mod 2=0)) then c3case:=0;  (* h0l:h+l *) (* glide plane -n- *)
                       if ((c3=1) and (k=0) and not (h mod 2=0)) then c3case:=0;  (* h0l:h *)  (* glide plane -a- *)
                       if ((c3=2) and (k=0) and not (l mod 2=0)) then c3case:=0;  (* h0l:l *)  (* glide plane -c- *)
                       if ((c3=3) and (k=0) and not (((h+l) mod 4=0) and (h mod 2=0) and (l mod 2=0))) then c3case:=0; (* h0l:h+l=4n,h,l *)  (* glide plane -d- *)
                       if ((c3=4) and (k=0) and not ((h mod 2=0) and (l mod 2=0))) then c3case:=0;

                       (* hk0 *)
                       if ((c4=0) and (l=0) and not ((h+k) mod 2=0)) then c4case:=0;  (* hk0:h+k *) (* glide plane --n *)
                       if ((c4=1) and (l=0) and not (h mod 2=0)) then c4case:=0;  (* hk0:h *)  (* glide plane --a *)
                       if ((c4=2) and (l=0) and not (k mod 2=0)) then c4case:=0;  (* hk0:k *)  (* glide plane --b *)
                       if ((c4=3) and (l=0) and not (((h+k) mod 4=0) and (h mod 2=0) and (k mod 2=0))) then c4case:=0; (* hk0:h+k=4n,h,k *)  (* glide plane --d *)
                       if ((c4=4) and (l=0) and not ((h mod 2=0) and (k mod 2=0))) then c4case:=0;


                       (* serial conditions *)
                       (* h00 *)
                       if ((c5=0) and (k=0) and (l=0) and not (h mod 2=0)) then c5case:=0;  (* h00:h *) (* screw axis 21-- ||a *)
                       if ((c5=1) and (k=0) and (l=0) and not (h mod 4=0)) then c5case:=0;  (* h00:4n *) (* screw axis 41-- ||a *)

                       (* 0k0 *)
                       if c6=0 then if ((c6=0) and (h=0) and (l=0) and not (k mod 2=0)) then c6case:=0;  (* 0k0:k *) (* screw axis -21- ||b *)
                       if c6=1 then if ((c6=1) and (h=0) and (l=0) and not (k mod 4=0)) then c6case:=0;  (* 0k0:4n *) (* screw axis -41- ||b *)

                       (* 00l *)
                       if ((c7=0) and (h=0) and (k=0) and not (l mod 2=0)) then c7case:=0;  (* 00l:h *) (* screw axis --21 ||c *)
                       if ((c7=1) and (h=0) and (k=0) and not (l mod 4=0)) then c7case:=0;  (* 00l:4n *) (* screw axis --41 ||c *)

                       ccase:=c1case*c2case*c3case*c4case*c5case*c6case*c7case;
                    end;*/

#ifndef __CUDACC__
                    if ( RadioButtonSysMonoclinic ) qDebug() << "RadioButtonSysMonoclinic";
#endif
                        /*(*** Monoclinic system ***)
                    if RadioButtonSysMonoclinic.Checked=true then begin
                       ccase:=0;
                       c1case:=0;      (* hkl *)
                       c2case:=1;      (* 0kl *)
                       c3case:=1;      (* h0l *)
                       c4case:=1;      (* hk0 *)
                       c5case:=1;      (* h00 *)
                       c6case:=1;      (* 0k0 *)
                       c7case:=1;      (* 00l *)

                       (* integral condition *)
                       if c1=-1 then c1case:=1; (* hkl *)
                       if c1=0 then if ((h+k+l) mod 2=0) then c1case:=1; (* h+k+l *)   (* I *)
                       if c1=1 then if ((h+k) mod 2=0) then c1case:=1;   (* h+k *)     (* C *)
                       if c1=2 then if ((k+l) mod 2=0) then c1case:=1;   (* k+l *)     (* A *)
                       if c1=3 then if ((h+l) mod 2=0) then c1case:=1;   (* h+l *)     (* B *)

                       (* zonal conditions *)
                       (* 0kl *)
                       if ((c2=0) and (h=0) and not ((k+l) mod 2=0)) then c2case:=0;  (* 0kl:k+l *) (* glide plane n *)
                       if ((c2=1) and (h=0) and not (k mod 2=0)) then c2case:=0;  (* 0kl:k *)  (* glide plane b *)
                       if ((c2=2) and (h=0) and not (l mod 2=0)) then c2case:=0;  (* 0kl:l *)  (* glide plane c *)
                       if ((c2=4) and (h=0) and not ((k mod 2=0) and (l mod 2=0))) then c2case:=0; (* 0kl:k,l *)

                       (* h0l *)
                       if ((c3=0) and (k=0) and not ((h+l) mod 2=0)) then c3case:=0;  (* h0l:h+l *) (* glide plane n *)
                       if ((c3=1) and (k=0) and not (h mod 2=0)) then c3case:=0;  (* h0l:h *)  (* glide plane -a- *)
                       if ((c3=2) and (k=0) and not (l mod 2=0)) then c3case:=0;  (* h0l:l *)  (* glide plane -c- *)
                       if ((c3=4) and (k=0) and not ((h mod 2=0) and (l mod 2=0))) then c3case:=0; (* h0l:h,l *)

                       (* hk0 *)
                       if ((c4=0) and (l=0) and not ((h+k) mod 2=0)) then c4case:=0;  (* hk0:h+k *) (* glide plane --n *)
                       if ((c4=1) and (l=0) and not (h mod 2=0)) then c4case:=0;  (* hk0:h *)  (* glide plane --a *)
                       if ((c4=2) and (l=0) and not (k mod 2=0)) then c4case:=0;  (* hk0:k *)  (* glide plane --b *)
                       if ((c4=4) and (l=0) and not ((h mod 2=0) and (k mod 2=0))) then c4case:=0;  (* hk0:h,k *)


                       (* serial conditions *)
                       (* h00 *)
                       if ((c5=0) and (k=0) and (l=0) and not (h mod 2=0)) then c5case:=0;  (* h00:h *) (* screw axis 21-- ||a *)

                       (* 0k0 *)
                       if ((c6=0) and (h=0) and (l=0) and not (k mod 2=0)) then c6case:=0;  (* 0k0:k *) (* screw axis -21- ||b *)

                       (* 00l *)
                       if c7=0 then if ((c7=0) and (h=0) and (k=0) and not (l mod 2=0)) then c7case:=0;  (* 00l:h *) (* screw axis --21 ||c *)

                       ccase:=c1case*c2case*c3case*c4case*c5case*c6case*c7case;
                    end;*/

#ifndef __CUDACC__
                    if ( RadioButtonSysTrigonal ) qDebug() << "RadioButtonSysTrigonal";
#endif
                    /*(*** Trigonal system ***)
                    if RadioButtonSysTrigonal.checked then begin
                       ccase:=0;
                       c1case:=0;      (* hkl *)
                       c7case:=1;      (* 00l *)
                       c9case:=1;      (* h-hl *)
                       c10case:=1;     (* hhl *)
                       c11case:=1;

                       (* integral condition *)
                       if c1=-1 then c1case:=1; (* hkl *)
                       if c1=4 then if ((-h+k+l) mod 3=0) then c1case:=1;   (* -h+k+l *)
                       if c1=5 then if ((h-k+l) mod 3=0) then c1case:=1;   (* h-k+l *)


                       (* zonal conditions *)
                       (* h-hl *)
                       if ((c9=0) and (h<>0) and (k=-h) and not (l mod 2=0)) then c9case:=0;      (* h-hl:l *)
                       if ((c9=1) and (h<>0) and (k=-h) and not ((h+l) mod 3=0)) then c9case:=0;  (* h-hl:h+l=3n *)
                       if ((c9=2) and (h<>0) and (k=-h) and not ((-h+l) mod 3=0)) then c9case:=0;  (* h-hl:-h+l=3n *)
                       if ((c9=3) and (h<>0) and (k=-h) and not (((h+l) mod 3=0) and (l mod 2=0))) then c9case:=0;  (* h-hl:h+l=3n, l *)

                       (* hhl *)
                       if ((c10=0) and (h=k) and (l<>0) and not (l mod 2=0)) then c10case:=0;     (* hhl:l *)
                       if ((c10=1) and (h=k) and (l<>0) and not (l mod 3=0)) then c10case:=0;     (* hhl:l=3n *)

                       (* hhh *)
                       if ((c11=0) and (h=k) and (k=l) and not (h mod 2=0)) then c11case:=0;


                       (* serial conditions *)

                       (* 00l *)
                       if ((c7=0) and (h=0) and (k=0) and not (l mod 2=0)) then c7case:=0;  (* 00l:l *)
                       if ((c7=2) and (h=0) and (k=0) and not (l mod 3=0)) then c7case:=0;  (* 00l:l=3n *)
                       if ((c7=3) and (h=0) and (k=0) and not (l mod 6=0)) then c7case:=0;  (* 00l:l=6n *)

                       ccase:=c1case*c7case*c9case*c10case*c11case;
                    end;*/

                    if ( RadioButtonSysTriclinic )
                        ccase=1;

                    if ((h==0) && (k==0) && (l==0)) ccase=0;

                    if ( ccase==1 ) //{NV}-32926
                    {
                        //fq=true;
#ifdef UseStringGrid12
                        StringGrid12[4] = "TRUE";   //StringGrid12.Cells[4,ct]:='true';
#endif
                        if ( ii==1 )
                        {
                            setLatpar1(  1,0,ct1);
                            setLatpar1(ct1,1,h);
                            setLatpar1(ct1,2,k);
                            setLatpar1(ct1,3,l);
                            setLatpar1(ct1,4,mult);
                            //qDebug() << "latpar1:" << ct1 << h << k << l << mult;
                            ct1++;
                        }
                        else //if ( ii==2 )
                        {
                            setLatpar2(  1,0,ct1);
                            setLatpar2(ct1,1,h);
                            setLatpar2(ct1,2,k);
                            setLatpar2(ct1,3,l);
                            setLatpar2(ct1,4,mult);
                            //qDebug() << "latpar2:" << ct1 << h << k << l << mult;
                            ct1++;
                        }
                    }
                    else
                    {
                        //fq=false;
#ifdef UseStringGrid12
                        StringGrid12[4] = "false";  //StringGrid12.Cells[4,ct]:='false';
#endif
                    }

#ifdef UseStringGrid12
                    ct++;
                    qDebug() << "SG12:" << StringGrid12.join(" ; ") << "ct="<<ct << "ct1="<<ct1 << "ii="<<ii;
#endif
                }
            }
        }  //(* of hkl-loop *)
    }  //(* of ii=1,2-case *)

} /* ButtonHKLClick() */



// Nur in prepareCalculation
void SasCalc_GENERIC_calculation::fhkl_c( int lat, int h, int k, int l,
                                         double &sphno, double &fhkl, double &qhkl,
                                         double &qhkl0 ) const
{
    // 'fhkl' wird im rufenden Programm zwar noch in einem Array (latpar?) gespeichert, aber dann
    // nirgendwo weiter verwendet. Daher wird es hier auch nicht gesetzt und kann langfristig
    // als Parameter entfernt werden - sofern Herr Förster keine andere Idee hat.

    int i,nn;
    double a,b,c,/*alf,bet,*/gam,sumc,sums,arg,invd;
    //double /*s2a,*/c1a,/*c2a,c3a,*/c1b,/*c2b,s2b,*/c1g; //,c2g,s2g;
    float r[13][4]; // array[1..12,1..3] of real;

    a=params.uca;
    b=params.ucb;
    c=params.ucc;
    //alf=params.ucalpha_deg*M_PI/180.;
    //bet=params.ucbeta_deg*M_PI/180.;
    gam=params.ucgamma_deg*M_PI/180.;

    //s2a=sin(alf)*sin(alf);
    //c1a=cos(alf);
    //c2a=c1a*c1a;
    //c3a=c1a*c1a*c1a;
    //c1b=cos(bet);
    //c2b=c1b*c1b;
    //s2b=sin(bet)*sin(bet);
    //c1g=cos(gam);
    //c2g=c1g*c1g;
    //s2g=sin(gam)*sin(gam);

    switch ( lat )
    {
    case 0:  //(* LAM *)
        sphno=1;
        fhkl=1;
        qhkl=2*M_PI*h/a;
        qhkl0=2*M_PI/a;
        break;

    case 1:  //(* hex cyl *)
        sphno=1;
        fhkl=1;
        invd=(h*h*b*b+2*h*k*a*b*cos(gam)+k*k*a*a)/(a*a*b*b*sin(gam)*sin(gam));
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI/a;
        break;

    case 2:  //(* sq cyl *)
        sphno=1;
        fhkl=1;
        invd=(h*h*b*b+2*h*k*a*b*cos(gam)+k*k*a*a)/(a*a*b*b*sin(gam)*sin(gam));
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI/a;
        break;

    case 3:  //(* c rect cyl *)
        sphno=2;
        nn=2;
        r[1][1]=0;    r[1][2]=0;
        r[2][1]=1/2.; r[2][2]=1/2.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);
        invd=(h*h*b*b+2*h*k*a*b*cos(gam)+k*k*a*a)/(a*a*b*b*sin(gam)*sin(gam));
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI*sqrt(2.)/a;
        break;

    case 4:  //(* BCC *)
        sphno=2;
        nn=2;
        r[1][1]=0;    r[1][2]=0;     r[1][3]=0;
        r[2][1]=1/2.; r[2][2]=1/2.;  r[2][3]=1/2.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]+l*r[i][3]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);

        //if odd(h+k+l) then fhkl:=0
        //   else fhkl:=4;

        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI*sqrt(2.)/a;
        break;

    case 5:  //(* FCC *)
    case 30: // Generic / Testmode
        sphno=4;
        nn=4;
        r[1][1]=0;    r[1][2]=0;     r[1][3]=0;
        r[2][1]=1/2.; r[2][2]=1/2.;  r[2][3]=0;
        r[3][1]=0;    r[3][2]=1/2.;  r[3][3]=1/2.;
        r[4][1]=1/2.; r[4][2]=0;     r[4][3]=1/2.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]+l*r[i][3]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);

        //if ((odd(h) and odd(k) and odd(l)) or ((not odd(h)) and (not odd(k)) and (not odd(l)))) then fhkl:=16
        //   else fhkl:=0;

        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI*sqrt(3.)/a;
        break;

    case 6:  //(* HCP *)
        sphno=2;
        nn=2;
        r[1][1]=0;    r[1][2]=0;     r[1][3]=0;
        r[2][1]=1/3.; r[2][2]=2/3.;  r[2][3]=1/2.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]+l*r[i][3]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);

        //fhkl:=0;
        //if ((((h+2*k) mod 3)=0) and odd(l)) then fhkl:=0;
        //if ((((h+2*k) mod 3)=0) and not(odd(l))) then fhkl:=4;
        //if ((fabs((h+2*k) mod 3)=1) and odd(l)) then fhkl:=3;
        //if ((fabs((h+2*k) mod 3)=1) and not(odd(l))) then fhkl:=1;

        invd=4*(h*h+k*k+h*k)/(3*a*a)+l*l/(c*c);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI/a;
        break;

    case 7:  //(* SC *)
        sphno=1;
        fhkl=1;
        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI/a;
        break;

    case 8:  //(* BCT *)
        sphno=2;
        nn=2;
        r[1][1]=0;    r[1][2]=0;     r[1][3]=0;
        r[2][1]=1/2.; r[2][2]=1/2.;  r[2][3]=1/2.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]+l*r[i][3]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);

        //if odd(h+k+l) then fhkl:=0
        //   else fhkl:=4;

        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);  // qhkl:=2*pi*sqrt(invd);
        qhkl0=2*M_PI*sqrt(2.)/a; // qhkl:=2*pi*sqrt(2)/a;
        break;

    case 9:  //(* Ia3d *)
        sphno=2;
        nn=2;
        r[1][1]=0;    r[1][2]=0;     r[1][3]=0;
        r[2][1]=1/2.; r[2][2]=1/2.;  r[2][3]=1/2.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]+l*r[i][3]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);

        //if odd(h+k+l) then fhkl:=0
        //   else fhkl:=4;

        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);  // qhkl:=2*pi*sqrt(invd);
        qhkl0=2*M_PI*sqrt(2.)/a; // qhkl:=2*pi*sqrt(2)/a;
        break;

    case 10: //(* Pn3m *)
        sphno=1;
        fhkl=1;
        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI/a;
        break;

    case 11: //(* Im3m *)
        sphno=2;
        nn=2;
        r[1][1]=0;    r[1][2]=0;     r[1][3]=0;
        r[2][1]=1/2.; r[2][2]=1/2.;  r[2][3]=1/2.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]+l*r[i][3]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);

        //if odd(h+k+l) then fhkl:=0
        //   else fhkl:=4;

        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI*sqrt(2.)/2.;
        break;

    case 17: //(* Fd3m *)
        sphno=8;
        nn=8;
        r[1][1]=0;    r[1][2]=0;     r[1][3]=0;
        r[2][1]=1/2.; r[2][2]=1/2.;  r[2][3]=0;
        r[3][1]=0;    r[3][2]=1/2.;  r[3][3]=1/2.;
        r[4][1]=1/2.; r[4][2]=0;     r[4][3]=1/2.;
        r[5][1]=1/4.; r[5][2]=1/4.;  r[5][3]=1/4.;
        r[6][1]=3/4.; r[6][2]=3/4.;  r[6][3]=1/4.;
        r[7][1]=1/4.; r[7][2]=3/4.;  r[7][3]=3/4.;
        r[8][1]=3/4.; r[8][2]=1/4.;  r[8][3]=3/4.;
        sumc=0.0;
        sums=0.0;
        for ( i=1; i<=nn; i++ )
        {
            arg=2*M_PI*(h*r[i][1]+k*r[i][2]+l*r[i][3]);
            sumc=sumc+cos(arg);
            sums=sums+sin(arg);
        }
        fhkl=round(sumc*sumc+sums*sums);

        //fhkl:=0;
        //if ((h+k+l) mod 4)=0 then fhkl:=64;
        //if ((h+k+l) mod 2)=1 then fhkl:=32;
        //if ((h+k+l) mod 4)=2 then fhkl:=0;

        invd=(h*h+k*k+l*l)/(a*a);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI*sqrt(3.)/a;
        break;

    case 18: //(* orthorhombic *)
        sphno=1;
        fhkl=1;
        invd=h*h/(a*a)+k*k/(b*b)+l*l/(c*c);
        qhkl=2*M_PI*sqrt(invd);
        qhkl0=2*M_PI/a;
        break;
    } // switch lat

} /* fhkl_c() */



// Nur in prepareCalculation verwendet
void SasCalc_GENERIC_calculation::extinction( int lat, int h, int k, int l, int aniso,
                                             int &mhkl, double &fhkl ) const
{
    //int ah,ak,al;

    switch ( lat )
    {
    case 0:     //(* LAM *)
        if (aniso==0)
        {
            mhkl=1;
            fhkl=1;
        }
        else if (aniso==1)
        {
            mhkl=1;
            fhkl=1;
        }
        break;

    case 1:     //(* HEX *)
        if (aniso==0)
        {
            //(*** Multiplicity ***)
            if (k==0 || k==h) mhkl=6;
            else mhkl=12;
            fhkl=1;
        }
        else if (aniso==1)
        {
            if (h==0 && k==0) mhkl=0;
            else mhkl=1;
            fhkl=1;
        }
        break;

    case 2:     //(* sq. Cyl *)
        if (aniso==0)
        {
            //(*** Multiplicity ***)
            if (k==0 || k==h) mhkl=4;
            else mhkl=8;
            fhkl=1;
        }
        else if (aniso==1)
        {
            if (h==0 && k==0) mhkl=0;
            else mhkl=1;
            fhkl=1;
        }
        break;

    case 3:     //(* centered rect. Cyl *)
        if (aniso==0)
        {
            //(*** Multiplicity ***)
            if (k==0 && k==0) mhkl=0;
            else mhkl=1;
            fhkl=1;
            if ( (h+k) % 2 !=0 ) fhkl=0;     //(* I *)
            if ( (h==0) && (k % 2 != 0) ) fhkl=0; //(* 0kl:k *)
            if ( (k==0) && (h % 2 != 0) ) fhkl=0; //(* h0l:h *)
            if ( (l==0) && ((h+k) % 2 !=0) ) fhkl=0; //(* hk0:h+k *)
            if ( (k==0) && (l==0) && (h % 2 != 0) ) fhkl=0; //(* h00:h *)
            if ( (h==0) && (l==0) && (k % 2 != 0) ) fhkl=0; //(* 0k0:k *)
        }
        else if (aniso==1)
        {
            if ( (h==0) && (k==0) ) mhkl=0;
            else mhkl=1;
            //(*** Extinction Rules ***)
            fhkl=1;
            if ( (h+k) % 2 != 0 ) fhkl=0;     //(* I *)
            if ( (h==0) && (k % 2 != 0) ) fhkl=0; //(* 0kl:k *)
            if ( (k==0) && (h % 2 != 0) ) fhkl=0; //(* h0l:h *)
            if ( (l==0) && ((h+k) % 2 != 0) ) fhkl=0; //(* hk0:h+k *)
            if ( (k==0) && (l==0) && (h % 2 != 0) ) fhkl=0; //(* h00:h *)
            if ( (h==0) && (l==0) && (k % 2 != 0) ) fhkl=0; //(* 0k0:k *)
        }
        break;

    case 4:     //(* BCC lattice *)
        if (aniso==0)
        {
            //(*** Multiplicity ***)
            /*h00*/if ( ((h==0) && (k==0)) || ((h==0) && (l==0)) || ((k==0) && (l==0)) ) mhkl=6;
            /*hh0*/if ( ((h==0) && (k==l) && (k!=0)) || ((k==0) && (h==l) && (h!=0))
                        || ((l==0) && (h==k) && (h!=0)) ) mhkl=12;
            /*hhh*/if ( (h==k) && (k==l) ) mhkl=8;
            /*hk0*/if ( ((h!=k) && (l==0) && (h!=0) && (k!=0)) || ((h!=l) && (k==0) && (h!=0) && (l!=0))
                        || ((k!=l) && (h==0) && (k!=0) && (l!=0)) ) mhkl=24;
            /*hhk*/if ( ((h==k) && (l!=0) && (l!=h) && (h!=0)) || ((h==l) && (k!=0) && (k!=h) && (h!=0))
                        || ((k==l) && (h!=0) && (h!=k) && (k!=0)) ) mhkl=24;
            /*hkl*/if ( (h!=k) && (k!=l) && (l!=h) && (h!=0) && (k!=0) && (l!=0) ) mhkl=48;
            //(*** Extinction Rules ***)
            fhkl=1+cospi(h+k+l); // 1+cos(M_PI*(h+k+l));
        }
        else if (aniso==1)
        {
            if ( (h==0) && (k==0) && (l==0) ) mhkl=0;
            else mhkl=1;
            //(* if ((ah+ak+al) mod 2) = 0 then fhkl:=2
            //    else fhkl:=0; *)
            fhkl=1+cospi(h+k+l); // 1+cos(M_PI*(h+k+l));
        }
        break;

    case 5:     //(* FCC lattice *)
        if (aniso==0)
        {
            //(*** Multiplicity ***)
            /*h00*/if ( ((h==0) && (k==0)) || ((h==0) && (l==0)) || ((k==0) && (l==0)) ) mhkl=6;
            /*hh0*/if ( ((h==0) && (k==l) && (k!=0)) || ((k==0) && (h==l) && (h!=0))
                        || ((l==0) && (h==k) && (h!=0)) ) mhkl=12;
            /*hhh*/if ( (h==k) && (k==l) ) mhkl=8;
            /*hk0*/if ( ((h!=k) && (l!=0) && (h!=0) && (k!=0)) || ((h!=l) && (k==0) && (h!=0) && (l!=0))
                        || ((k!=l) && (h==0) && (k!=0) && (l!=0)) ) mhkl=24;
            /*hhk*/if ( ((h==k) && (l!=0) && (l!=h) && (h!=0)) || ((h==l) && (k!=0) && (k!=h) && (h!=0))
                        || ((k==l) && (h!=0) && (h!=k) && (k!=0)) ) mhkl=24;
            /*hkl*/if ( (h!=k) && (k!=l) && (l!=h) && (h!=0) && (k!=0) && (l!=0) ) mhkl=48;

            //(*** Extinction Rules ***)
            //(* if (((h mod 2)=0) and ((k mod 2)=0) and ((l mod 2)=0)) or
            //    (((h mod 2)=1) and ((k mod 2)=1) and ((l mod 2)=1)) then fhkl:=4
            //    else fhkl:=0; *)
            //fhkl=1+cos(M_PI*(h+k))+cos(M_PI*(h+l))+cos(M_PI*(k+l));
            fhkl=1+cospi(h+k)+cospi(h+l)+cospi(k+l);
        }
        else if (aniso==1)
        {
            if ((h==0) && (k==0) && (l==0)) mhkl=0;
            else mhkl=1;
            if ( ((odd(h) && odd(k) && odd(l)) ||
                 ((even(h)) && (even(k)) && (even(l)))) ) fhkl=4;
            else fhkl=0;
        }
        break;

    case 6:     // (* HCP lattice *)
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            if ( l!=0 )
            {
                /*00l*/if ( (h==0) && (k==0) ) mhkl = 2;
                /*h0l*/if ( (((h==0) && (k!=0)) || ((h!=0) && (k==0))) ) mhkl = 8;
                /*hhl*/if ( ((h!=0) && (k==h)) ) mhkl = 8;
                /*hkl*/if ( ((h!=0) && (k!=0) && (k!=h)) ) mhkl = 16;
            }
            else // if ( (l==0) )
            {
                /*h00*/if ( (((h!=0) && (k==0)) || ((h==0) && (k!=0))) ) mhkl = 4;
                /*hh0*/if ( ((h!=0) && (h==k)) ) mhkl = 4;
                /*hk0*/if ( ((h!=0) && (k!=0) && (h!=k)) ) mhkl = 8;
                /*000*/if ( ((h==0) && (k==0)) ) mhkl = 0;
            }

            /* Exctinction rules */
            fhkl = 1;
            /* hhl:l */
            if ( (h==k) && (l!=0) && ! (l % 2==0) ) fhkl = 0; /* hhl:l */
            /* 00l:l */
            if ( (h==0) && (k==0) && ! (l % 2==0) ) fhkl = 0; /* 00l:l */

            /*fhkl:=2*cos(2*pi*(((h+2*k)/3)+(l/4)));
              fhkl2:=fhkl*fhkl;
              if (fhkl2<0.1) then fhkl2:=0;  */
        }
        else if ( aniso==1 )
        {
            if ( (h==0) && (k==0) && (l==0) ) mhkl = 0;
            else mhkl = 1;

            fhkl = 2*cos(2*M_PI*(((h+2*k)/3.0)+(l/4.0)));
            fhkl = fhkl*fhkl;
        }
        break;

    case 7:     // (* SC lattice *)
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            /*h00*/if ( ((h==0) && (k==0)) || ((h==0) && (l==0)) || ((k==0) && (l==0)) ) mhkl = 6;
            /*hh0*/if ( ((h==0) && (k==l) && (k!=0)) || ((k==0) && (h==l) && (h!=0))
                        || ((l==0) && (h==k) && (h!=0)) ) mhkl = 12;
            /*hhh*/if ( ((h==k) && (k==l)) ) mhkl = 8;
            /*hk0*/if ( ((h!=k) && (l==0) && (h!=0) && (k!=0)) || ((h!=l) && (k==0) && (h!=0) && (l!=0))
                        || ((k!=l) && (h==0) && (k!=0) && (l!=0)) ) mhkl = 24;
            /*hhk*/if ( ((h==k) && (l!=0) && (l!=h) && (h!=0)) || ((h==l) && (k!=0) && (k!=h) && (h!=0))
                        || ((k==l) && (h!=0) && (h!=k) && (k!=0)) ) mhkl = 24;
            /*hkl*/if ( ((h!=k) && (k!=l) && (l!=h) && (h!=0) && (k!=0) && (l!=0)) ) mhkl = 48;

            fhkl = 1;
        }
        if ( aniso==1 )
        {
            if ( (h==0) && (k==0) && (l==0) ) mhkl = 0;
            else mhkl = 1;
            fhkl = 1;
        }
        break;

    case 8:     // (* BCT lattice *)
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            if ( l!=0 )
            {
                /*00l*/if ( (h==0) && (k==0) ) mhkl = 2;
                /*h0l*/if ( ((h==0) && (k!=0)) || ((h!=0) && (k==0)) ) mhkl = 8;
                /*hhl*/if ( (h!=0) && (k==h) ) mhkl = 8;
                /*hkl*/if ( (h!=0) && (k!=0) && (k!=h) ) mhkl = 16;
            }
            else // if ( l==0 )
            {
                /*h00*/if ( ((h!=0) && (k==0)) || ((h==0) && (k!=0)) ) mhkl = 4;
                /*hh0*/if ( (h!=0) && (h==k) ) mhkl = 4;
                /*hk0*/if ( (h!=0) && (k!=0) && (h!=k) ) mhkl = 8;
                /*000*/if ( (h==0) && (k==0) ) mhkl = 0;
            }
            /*** Exctinction Rules ***/
            fhkl = 1+cospi(h+k+l);
        }
        if ( aniso==1 )
        {
            if ( (h==0) && (k==0) && (l==0) ) mhkl = 0;
            else mhkl = 1;

            fhkl = 1+cospi(h+k+l);
        }
        break;

    case 9: /* Ia3d lattice */
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            /*h00*/if ( ((h==0) && (k==0)) || ((h==0) && (l==0)) || ((k==0) && (l==0)) ) mhkl = 6;
            /*hh0*/if ( ((h==0) && (k==l) && (k!=0)) || ((k==0) && (h==l) && (h!=0))
                        || ((l==0) && (h==k) && (h!=0)) ) mhkl = 12;
            /*hhh*/if ( ((h==k) && (k==l)) ) mhkl = 8;
            /*hk0*/if ( ((h!=k) && (l==0) && (h!=0) && (k!=0)) || ((h!=l) && (k==0) && (h!=0) && (l!=0))
                        || ((k!=l) && (h==0) && (k!=0) && (l!=0)) ) mhkl = 24;
            /*hhk*/if ( ((h==k) && (l!=0) && (l!=h) && (h!=0)) || ((h==l) && (k!=0) && (k!=h) && (h!=0))
                        || ((k==l) && (h!=0) && (h!=k) && (k!=0)) ) mhkl = 24;
            /*hkl*/if ( ((h!=k) && (k!=l) && (l!=h) && (h!=0) && (k!=0) && (l!=0)) ) mhkl = 48;

            fhkl = 0;
            /**** hkl: h+k+l=2n condition ****/
            if ( (h!=0) && (k!=0) && (l!=0) &&
                (h!=k) && (k!=l) && (h!=l) &&
                (((h+k+l) % 2) == 0) ) fhkl = 2;

            /**** hhl: 2*h+l=4n, l=2n condition ****/
            if ( (h==k) && (h!=0) && (k!=0) && (l!=0) &&
                (l % 2 == 0) && (((2*h+l) % 4) == 0) ) fhkl = 2;
            if ( (h==l) && (h!=0) && (k!=0) && (l!=0) &&
                (k % 2 == 0) && (((2*h+k) % 4) == 0) ) fhkl = 2;
            if ( (k==l) && (h!=0) && (k!=0) && (l!=0) &&
                (h % 2 == 0) && (((2*k+h) % 4) == 0) ) fhkl = 2;

            /**** hk0: h=2n, k=2n, h+k=2n condition ****/
            if ( (h==0) && (k!=0) && (l!=0) &&
                (k % 2 == 0) && (l % 2 == 0) ) fhkl = 2;
            if ( (k==0) && (h!=0) && (l!=0) &&
                (h % 2 == 0) && (l % 2 == 0) ) fhkl = 2;
            if ( (l==0) && (h!=0) && (k!=0) &&
                (h % 2 == 0) && (k % 2 == 0) ) fhkl = 2;

            /**** h00: h=4n condition ****/
            if ( (h==0) && (k==0) && (l!=0) && ((l % 4) == 0) ) fhkl = 2;
            if ( (h==0) && (l==0) && (k!=0) && ((k % 4) == 0) ) fhkl = 2;
            if ( (k==0) && (l==0) && (h!=0) && ((h % 4) == 0) ) fhkl = 2;

        }
        else if ( aniso==1 )
        {
            int ah = fabs(h);
            int ak = fabs(k);
            int al = fabs(l);

            mhkl = 0;
            /*h00*/if ( ((ah==0) && (ak==0)) || ((ah==0) && (al==0)) || ((ak==0) && (al==0)) ) mhkl = 1;
            /*hh0*/if ( ((ah==0) && (ak==al) && (ak!=0)) || ((ak==0) && (ah==al) && (ah!=0))
                        || ((al==0) && (ah==ak) && (ah!=0)) ) mhkl = 1;
            /*hhh*/if ( ((ah==ak) && (ak==al)) ) mhkl = 1;
            /*hk0*/if ( ((ah!=ak) && (al==0) && (ah!=0) && (ak!=0)) || ((ah!=al) && (ak==0) && (ah!=0) && (al!=0))
                        || ((ak!=al) && (ah==0) && (ak!=0) && (al!=0)) ) mhkl = 1;
            /*hhk*/if ( ((ah==ak) && (al!=0) && (al!=ah) && (ah!=0)) || ((ah==al) && (ak!=0) && (ak!=ah) && (ah!=0))
                        || ((ak==al) && (ah!=0) && (ah!=ak) && (ak!=0)) ) mhkl = 1;
            /*hkl*/if ( ((ah!=ak) && (ak!=al) && (al!=ah) && (ah!=0) && (ak!=0) && (al!=0)) ) mhkl = 1;

            fhkl = 0;
            /**** hkl: h+k+l=2n condition ****/
            if ( (ah!=0) && (ak!=0) && (al!=0) &&
                (ah!=ak) && (ak!=al) && (ah!=al) &&
                (((ah+ak+al) % 2) == 0) ) fhkl = 2;

            /**** hhl: 2*h+l=4n, l=2n condition ****/
            if ( (ah==ak) && (ah!=0) && (ak!=0) && (al!=0) &&
                (al % 2 == 0) && (((2*ah+al) % 4) == 0) ) fhkl = 2;
            if ( (ah==al) && (ah!=0) && (ak!=0) && (al!=0) &&
                (ak % 2 == 0) && (((2*ah+ak) % 4) == 0) ) fhkl = 2;
            if ( (ak==al) && (ah!=0) && (ak!=0) && (al!=0) &&
                (ah % 2 == 0) && (((2*ak+ah) % 4) == 0) ) fhkl = 2;

            /**** hk0: h=2n, k=2n, h+k=2n condition ****/
            if ( (ah==0) && (ak!=0) && (al!=0) &&
                (ak % 2 == 0) && (al % 2 == 0) ) fhkl = 2;
            if ( (ak==0) && (ah!=0) && (al!=0) &&
                (ah % 2 == 0) && (al % 2 == 0) ) fhkl = 2;
            if ( (al==0) && (ah!=0) && (ak!=0) &&
                (ah % 2 == 0) && (ak % 2 == 0) ) fhkl = 2;

            /**** h00: h=4n condition ****/
            if ( (ah==0) && (ak==0) && (al!=0) && ((al % 4) == 0) ) fhkl = 2;
            if ( (ah==0) && (al==0) && (ak!=0) && ((ak % 4) == 0) ) fhkl = 2;
            if ( (ak==0) && (al==0) && (ah!=0) && ((ah % 4) == 0) ) fhkl = 2;
        }
        break;

    case 10:    /* Pn3m lattice */
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            /*h00*/if ( ((h==0) && (k==0)) || ((h==0) && (l==0)) || ((k==0) && (l==0)) ) mhkl = 6;
            /*hh0*/if ( ((h==0) && (k==l) && (k!=0)) || ((k==0) && (h==l) && (h!=0))
                        || ((l==0) && (h==k) && (h!=0)) ) mhkl = 12;
            /*hhh*/if ( ((h==k) && (k==l)) ) mhkl = 8;
            /*hk0*/if ( ((h!=k) && (l==0) && (h!=0) && (k!=0)) || ((h!=l) && (k==0) && (h!=0) && (l!=0))
                        || ((k!=l) && (h==0) && (k!=0) && (l!=0)) ) mhkl = 24;
            /*hhk*/if ( ((h==k) && (l!=0) && (l!=h) && (h!=0)) || ((h==l) && (k!=0) && (k!=h) && (h!=0))
                        || ((k==l) && (h!=0) && (h!=k) && (k!=0)) ) mhkl = 24;
            /*hkl*/if ( ((h!=k) && (k!=l) && (l!=h) && (h!=0) && (k!=0) && (l!=0)) ) mhkl = 48;

            fhkl = 0;
            /**** hkl: general reflection condition ****/
            if ( (h!=0) && (k!=0) && (l!=0) ) fhkl = 2;

            /**** hk0: h+k=2n condition ****/
            if ( (l==0) && (((h+k) % 2) == 0) ) fhkl = 2;
            if ( (k==0) && (((h+l) % 2) == 0) ) fhkl = 2;
            if ( (h==0) && (((k+l) % 2) == 0) ) fhkl = 2;

            /**** h00: h=2n condition ****/
            if ( (k==0) && (l==0) && ((h % 2) == 0) ) fhkl = 2;
            if ( (h==0) && (l==0) && ((k % 2) == 0) ) fhkl = 2;
            if ( (h==0) && (k==0) && ((l % 2) == 0) ) fhkl = 2;
        }
        else if ( aniso==1 )
        {
            int ah = fabs(h);
            int ak = fabs(k);
            int al = fabs(l);

            mhkl = 0;
            /*h00*/if ( ((ah==0) && (ak==0)) || ((ah==0) && (al==0)) || ((ak==0) && (al==0)) ) mhkl = 1;
            /*hh0*/if ( ((ah==0) && (ak==al) && (ak!=0)) || ((ak==0) && (ah==al) && (ah!=0))
                        || ((al==0) && (ah==ak) && (ah!=0)) ) mhkl = 1;
            /*hhh*/if ( ((ah==ak) && (ak==al)) ) mhkl = 1;
            /*hk0*/if ( ((ah!=ak) && (al==0) && (ah!=0) && (ak!=0)) || ((ah!=al) && (ak==0) && (ah!=0) && (al!=0))
                        || ((ak!=al) && (ah==0) && (ak!=0) && (al!=0)) ) mhkl = 1;
            /*hhk*/if ( ((ah==ak) && (al!=0) && (al!=ah) && (ah!=0)) || ((ah==al) && (ak!=0) && (ak!=ah) && (ah!=0))
                        || ((ak==al) && (ah!=0) && (ah!=ak) && (ak!=0)) ) mhkl = 1;
            /*hkl*/if ( ((ah!=ak) && (ak!=al) && (al!=ah) && (ah!=0) && (ak!=0) && (al!=0)) ) mhkl = 1;

            fhkl = 0;
            /**** hkl: general reflection condition ****/
            if ( (ah!=0) && (ak!=0) && (al!=0) ) fhkl = 2;

            /**** hk0: h+k=2n condition ****/
            if ( (al==0) && (((ah+ak) % 2) == 0) ) fhkl = 2;
            if ( (ak==0) && (((ah+al) % 2) == 0) ) fhkl = 2;
            if ( (ah==0) && (((ak+al) % 2) == 0) ) fhkl = 2;

            /**** h00: h=2n condition ****/
            if ( (ak==0) && (al==0) && ((ah % 2) == 0) ) fhkl = 2;
            if ( (ah==0) && (al==0) && ((ak % 2) == 0) ) fhkl = 2;
            if ( (ah==0) && (ak==0) && ((al % 2) == 0) ) fhkl = 2;

            if ( ((h==0) && (k==0) && (l==0)) ) fhkl = 0;
        }
        break;

    case 11:    /* Im3m lattice */
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            /*h00*/if ( ((h==0) && (k==0)) || ((h==0) && (l==0)) || ((k==0) && (l==0)) ) mhkl = 6;
            /*hh0*/if ( ((h==0) && (k==l) && (k!=0)) || ((k==0) && (h==l) && (h!=0))
                        || ((l==0) && (h==k) && (h!=0)) ) mhkl = 12;
            /*hhh*/if ( ((h==k) && (k==l)) ) mhkl = 8;
            /*hk0*/if ( ((h!=k) && (l==0) && (h!=0) && (k!=0)) || ((h!=l) && (k==0) && (h!=0) && (l!=0))
                        || ((k!=l) && (h==0) && (k!=0) && (l!=0)) ) mhkl = 24;
            /*hhk*/if ( ((h==k) && (l!=0) && (l!=h) && (h!=0)) || ((h==l) && (k!=0) && (k!=h) && (h!=0))
                        || ((k==l) && (h!=0) && (h!=k) && (k!=0)) ) mhkl = 24;
            /*hkl*/if ( ((h!=k) && (k!=l) && (l!=h) && (h!=0) && (k!=0) && (l!=0)) ) mhkl = 48;

            fhkl = 0;
            /**** hkl: h+k+l=2n condition (Bravais centering) ****/
            /**** 321, 431, 521, ... ****/
            if ( (h!=0) && (k!=0) && (l!=0) &&
                (h!=k) && (k!=l) && (h!=l) &&
                (((h+k+l) % 2) == 0) ) fhkl = 2;

            /**** hk0: h+k=2n condition ****/
            /**** 110, 220, 310, ... ****/
            if ( (l==0) && (h!=0) && (k!=0) && (((h+k) % 2) == 0) ) fhkl = 2;
            if ( (k==0) && (h!=0) && (l!=0) && (((h+l) % 2) == 0) ) fhkl = 2;
            if ( (h==0) && (k!=0) && (l!=0) && (((k+l) % 2) == 0) ) fhkl = 2;

            /**** hhl: l=2n condition ****/
            /**** 110, 220, 330, 440, ... , 112, 222, 332, ... ****/

            if ( (h!=0) && (k!=0) &&
                (h==k) && ((l % 2) == 0) ) fhkl = 2;
            if ( (h!=0) && (l!=0) &&
                (h==l) && ((k % 2) == 0) ) fhkl = 2;
            if ( (k!=0) && (l!=0) &&
                (k==l) && ((h % 2) == 0) ) fhkl = 2;

            /**** h00: h=2n condition ****/
            /**** 200, 400, 600, 800, ... ****/

            if ( (h==0) && (k==0) && (l!=0) && ((l % 2) == 0) ) fhkl = 2;
            if ( (h==0) && (l==0) && (k!=0) && ((k % 2) == 0) ) fhkl = 2;
            if ( (k==0) && (l==0) && (h!=0) && ((h % 2) == 0) ) fhkl = 2;
        }
        if ( aniso==1 )
        {
            int ah = fabs(h);
            int ak = fabs(k);
            int al = fabs(l);
            mhkl = 0;
            /*h00*/if ( ((ah==0) && (ak==0)) || ((ah==0) && (al==0)) || ((ak==0) && (al==0)) ) mhkl = 1;
            /*hh0*/if ( ((ah==0) && (ak==al) && (ak!=0)) || ((ak==0) && (ah==al) && (ah!=0))
                        || ((al==0) && (ah==ak) && (ah!=0)) ) mhkl = 1;
            /*hhh*/if ( ((ah==ak) && (ak==al)) ) mhkl = 1;
            /*hk0*/if ( ((ah!=ak) && (al==0) && (ah!=0) && (ak!=0)) || ((ah!=al) && (ak==0) && (ah!=0) && (al!=0))
                        || ((ak!=al) && (ah==0) && (ak!=0) && (al!=0)) ) mhkl = 1;
            /*hhk*/if ( ((ah==ak) && (al!=0) && (al!=ah) && (ah!=0)) || ((ah==al) && (ak!=0) && (ak!=ah) && (ah!=0))
                        || ((ak==al) && (ah!=0) && (ah!=ak) && (ak!=0)) ) mhkl = 1;
            /*hkl*/if ( ((ah!=ak) && (ak!=al) && (al!=ah) && (ah!=0) && (ak!=0) && (al!=0)) ) mhkl = 1;

            fhkl = 0;
            /**** hkl: h+k+l=2n condition (Bravais centering) ****/
            /**** 321, 431, 521, ... ****/

            if ( (ah!=0) && (ak!=0) && (al!=0) &&
                (ah!=ak) && (ak!=al) && (ah!=al) &&
                (((ah+ak+al) % 2) == 0) ) fhkl = 2;

            /**** hk0: h+k=2n condition ****/
            /**** 110, 220, 310, ... ****/

            if ( (al==0) && (ah!=0) && (ak!=0) && (((ah+ak) % 2) == 0) ) fhkl = 2;
            if ( (ak==0) && (ah!=0) && (al!=0) && (((ah+al) % 2) == 0) ) fhkl = 2;
            if ( (ah==0) && (ak!=0) && (al!=0) && (((ak+al) % 2) == 0) ) fhkl = 2;

            /**** hhl: l=2n condition ****/
            /**** 110, 220, 330, 440, ... , 112, 222, 332, ... ****/

            if ( (ah!=0) && (ak!=0) &&
                (ah==ak) && ((al % 2) == 0) ) fhkl = 2;
            if ( (ah!=0) && (al!=0) &&
                (ah==al) && ((ak % 2) == 0) ) fhkl = 2;
            if ( (ak!=0) && (al!=0) &&
                (ak==al) && ((ah % 2) == 0) ) fhkl = 2;

            /**** h00: h=2n condition ****/
            /**** 200, 400, 600, 800, ... ****/

            if ( (ah==0) && (ak==0) && (al!=0) && ((al % 2) == 0) ) fhkl = 2;
            if ( (ah==0) && (al==0) && (ak!=0) && ((ak % 2) == 0) ) fhkl = 2;
            if ( (ak==0) && (al==0) && (ah!=0) && ((ah % 2) == 0) ) fhkl = 2;
        }
        break;

    case 13:    /* CP-Layers */
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            if ( (k==0) || (k==h) ) mhkl = 6;
            else mhkl = 12;
            fhkl = 1;
        }
        else if ( aniso==1 )
        {
            if ( (h==0) && (k==0) ) mhkl = 0;
            else mhkl = 1;
            fhkl = 1;
        }
        break;

    case 17:    /* Fd3m lattice */
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            /*h00*/if ( ((h==0) && (k==0)) || ((h==0) && (l==0)) || ((k==0) && (l==0)) ) mhkl = 6;
            /*hh0*/if ( ((h==0) && (k==l) && (k!=0)) || ((k==0) && (h==l) && (h!=0))
                        || ((l==0) && (h==k) && (h!=0)) ) mhkl = 12;
            /*hhh*/if ( ((h==k) && (k==l)) ) mhkl = 8;
            /*hk0*/if ( ((h!=k) && (l==0) && (h!=0) && (k!=0)) || ((h!=l) && (k==0) && (h!=0) && (l!=0))
                        || ((k!=l) && (h==0) && (k!=0) && (l!=0)) ) mhkl = 24;
            /*hhk*/if ( ((h==k) && (l!=0) && (l!=h) && (h!=0)) || ((h==l) && (k!=0) && (k!=h) && (h!=0))
                        || ((k==l) && (h!=0) && (h!=k) && (k!=0)) ) mhkl = 24;
            /*hkl*/if ( ((h!=k) && (k!=l) && (l!=h) && (h!=0) && (k!=0) && (l!=0)) ) mhkl = 48;

            /* exctinction rules */
            fhkl = 1;
            /* F centering */
            if ( ! (((h+k) % 2==0) && ((h+l) % 2==0) && ((k+l) % 2==0)) ) fhkl = 0;  /* F */

            /* 0kl:k+l=4n,k,l */
            if ( ((h==0) && ! (((k+l) % 4==0) && (k % 2==0) && (l % 2==0))) ) fhkl = 0; /* 0kl */
            if ( ((k==0) && ! (((h+l) % 4==0) && (h % 2==0) && (l % 2==0))) ) fhkl = 0; /* h0l */
            if ( ((l==0) && ! (((h+k) % 4==0) && (h % 2==0) && (k % 2==0))) ) fhkl = 0; /* hk0 */

            /* hhl:h+l */
            if ( ((h==k) && (l!=0) && ! ((h+l) % 2==0)) ) fhkl = 0;  /* hhl:h+l */
            if ( ((h==l) && (k!=0) && ! ((h+k) % 2==0)) ) fhkl = 0;  /* hkh:h+k */
            if ( ((k==l) && (h!=0) && ! ((h+k) % 2==0)) ) fhkl = 0;  /* hkk:h+k */

            /* 00l:4n */
            if ( ((h==0) && (k==0) && ! (l % 4==0)) ) fhkl = 0; /* 00l:4n */
            if ( ((h==0) && (l==0) && ! (k % 4==0)) ) fhkl = 0; /* 0k0:4n */
            if ( ((k==0) && (l==0) && ! (h % 4==0)) ) fhkl = 0; /* h00:4n */
        }
        else if ( aniso==1 )
        {
            if ( ((h==0) && (k==0) && (l==0)) ) mhkl = 0;
            else mhkl = 1;

            /* exctinction rules */
            fhkl = 1;
            /* F centering */
            if ( ! (((h+k) % 2==0) && ((h+l) % 2==0) && ((k+l) % 2==0)) ) fhkl = 0;  /* F */

            /* 0kl:k+l=4n,k,l */
            if ( ((h==0) && ! (((k+l) % 4==0) && (k % 2==0) && (l % 2==0))) ) fhkl = 0; /* 0kl */
            if ( ((k==0) && ! (((h+l) % 4==0) && (h % 2==0) && (l % 2==0))) ) fhkl = 0; /* h0l */
            if ( ((l==0) && ! (((h+k) % 4==0) && (h % 2==0) && (k % 2==0))) ) fhkl = 0; /* hk0 */

            /* hhl:h+l */
            if ( ((h==k) && (l!=0) && ! ((h+l) % 2==0)) ) fhkl = 0;  /* hhl:h+l */
            if ( ((h==l) && (k!=0) && ! ((h+k) % 2==0)) ) fhkl = 0;  /* hkh:h+k */
            if ( ((k==l) && (h!=0) && ! ((h+k) % 2==0)) ) fhkl = 0;  /* hkk:h+k */

            /* 00l:4n */
            if ( ((h==0) && (k==0) && ! (l % 4==0)) ) fhkl = 0; /* 00l:4n */
            if ( ((h==0) && (l==0) && ! (k % 4==0)) ) fhkl = 0; /* 0k0:4n */
            if ( ((k==0) && (l==0) && ! (h % 4==0)) ) fhkl = 0; /* h00:4n */
        }
        break;

    case 18:    /* Orthorhombic lattice */
        if ( aniso==0 )
        {
            /*** Multiplicity ***/
            if ( l!=0 )
            {
                /*00l*/if ( ((h==0) && (k==0)) ) mhkl = 2;
                /*h0l*/if ( (((h==0) && (k!=0)) || ((h!=0) && (k==0))) ) mhkl = 8;
                /*hhl*/if ( ((h!=0) && (k==h)) ) mhkl = 8;
                /*hkl*/if ( ((h!=0) && (k!=0) && (k!=h)) ) mhkl = 16;
            }
            else // if ( l==0 )
            {
                /*h00*/if ( (((h!=0) && (k==0)) || ((h==0) && (k!=0))) ) mhkl = 4;
                /*hh0*/if ( ((h!=0) && (h==k)) ) mhkl = 4;
                /*hk0*/if ( ((h!=0) && (k!=0) && (h!=k)) ) mhkl = 8;
                /*000*/if ( ((h==0) && (k==0)) ) mhkl = 0;
            }
            /*** Exctinction Rules ***/
            fhkl = 1+cos(M_PI*(h+k+l));
        }
        else if ( aniso==1 )
        {
            mhkl = 1;
            if ( ((h==0) && (k==0) && (l==0)) ) mhkl = 0;

            /* P-type */
            fhkl = 1;
            if ( (h==0) && (k==0) && odd(l) ) fhkl = 0;
            if ( (h==0) && (l==0) && odd(k) ) fhkl = 0;
            if ( (k==0) && (l==0) && odd(h) ) fhkl = 0;

            /*  fhkl:=1+cos(pi*(h+k+l));                                  (* for C-type *) */
            /*  fhkl:=1+cos(pi*(h+k))+cos(pi*(h+l))+cos(pi*(k+l));        (* for F-type *) */
            /*  ...                                                       (* for E-type *) */
        }
        break;

    default:
        DD( qDebug() << "***extinction*** lat unknown" << lat );
        break;

    } // switch lat

} /* extinction() */


double SasCalc_GENERIC_calculation::gammln( double xx ) const
{
    const double stp = 2.50662827465;

    double x,tmp,ser;

    x=xx-1.0;
    tmp=x+5.5;
    tmp=(x+0.5)*log(tmp)-tmp;
    ser=1.0+76.18009173/(x+1.0)-86.50532033/(x+2.0)+24.01409822/(x+3.0)
          -1.231739516/(x+4.0)+0.120858003e-2/(x+5.0)-0.536382e-5/(x+6.0);
    return tmp+log(stp*ser);
}


void SasCalc_GENERIC_calculation::scaleIntensity( bool linlog )
{
    double vmin=1e20, vmax=1e-20;
    double *p=(arrXYIntensity+3);
    if ( linlog )
    {   // Lin
        for ( size_t i=0; i<_arrCount; i++,p++ )
            if ( *p < vmin ) vmin = *p;
            else if ( *p > vmax ) vmax = *p;
    }
    else
    {   // Log
        for ( size_t i=0; i<_arrCount; i++,p++ )
            if ( *p > 0 )
            {
                if ( log10(*p) < vmin )
                    vmin = log10(*p);
                else if ( log10(*p) > vmax )
                    vmax = log10(*p);
            }
    }
    if ( fabs(vmin-vmax) < 1e-6 )
    {   // Bei lin und log ...
        vmin *= 0.998;
        vmax *= 1.002;
        if ( vmax == 0 ) vmax = 0.2;
    }
    p = (arrXYIntensity+3);
    if ( linlog )
    {   // Lin
        for ( size_t i=0; i<_arrCount; i++,p++ )
            *p = (*p - vmin) / (vmax-vmin);
    }
    else
    {   // Log
        for ( size_t i=0; i<_arrCount; i++,p++ )
            if ( *p > 0 )
                *p = (log10(*p) - vmin) / (vmax-vmin);
            else
                *p = 0;
    }
} // scaleIntensity()
