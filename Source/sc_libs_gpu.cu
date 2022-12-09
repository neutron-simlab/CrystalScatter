/*
 * Collection of common calculatoion routines.
 * Advantage:    avoid code duplicates and save program memory
 * Disadvantage: no code optimization due to partly used functions possible
 */

// THIS IS ONLY AN INCLUDE FILE FOR ANOTHER CLASS



// Um den Code in den jeweiligen Berechnungsklassen doch etwas
// kleiner zu halten, habe ich jede Routine mittels #ifdef ...
// ausgeblendet und nur die verwendeten können eingeblendet werden.
// Damit der QtCreator damit keine Probleme bekommt, und es leichter
// zu schreiben ist, werden hier alle eingeblendet.
/* */
#ifndef USE_CLASSLIB
#include "sc_libs_gpu.h"

#define USE_szave
#define USE_psphere // nur init genutzt
#define USE_psphered
#define USE_pspheredf
#define USE_f2dschulz
#define USE_f2dhyper
#define USE_gamma
#define USE_cosav
#define USE_gammaratio
#define USE_pqcoreshell
#define USE_f2dcoreshell
#define USE_polyvesicle
#define USE_f2dpolyvesicle
#define USE_polyliposome
#define USE_polycube
#define USE_burger
#define USE_angleuv
#define USE_lorentznorm3
#define USE_gaussnorm3
#define USE_pearsonnorm3
#define USE_pearsonintegral3
#define USE_cosavm
#define USE_pgensphere1
#define USE_hypcoreshellm
#define USE_trapcube
#define USE_cubealf
#define USE_trapcubebet
#define USE_cube
#define USE_trapezpqr
#define USE_midpntchi
#define USE_qrombpq
#define USE_qrombpi2
#define USE_bessel10
#define USE_ButtonHKLClick
#define USE_fhkl_c
#define USE_extinction
#define USE_trapzddeltac
#define USE_polint
#define USE_qrombdeltac
#define USE_qrombchid
#define USE_trapzdchid
#endif // USE_CLASSLIB
/* */


#ifdef __CUDACC__
#define DBGNAN(x)
#define DBG(x)
#else
#define DBGNAN(x) //x
#define DBG(x)    //x
#endif


/**
 * @brief endThread
 * Helper function to cancel all calculation threads.
 * The QtCreator marks many undefined symbols because this is an include. It compiles without errors.
 */
void CLASSLIB::endThread()
{
    if ( !_endThread )
    {   // Jetzt wirklich die Threads killen
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
                for ( int i=0; i<numberOfThreads; i++ )
                    if ( threads != nullptr && threads[i] != 0 )
                    {
                        std::cerr << "endThread() Cancel Thread #" << i << std::endl;
                        pthread_cancel( threads[i] );
                        //pthread_kill( threads[i], SIGTERM );
                        threads[i] = 0;
                    }
            }
    }
}



#ifdef USE_szave
/**
 * @brief CLASSLIB::szave - calculates SZ-averages of trigonometric functions
 * @param c - select the x^n function to use
 * @param a - factor in the trigonometric functions
 * @param x - q sphere radius
 * @param n - select the gamma functions to use
 * @param z - (1 - sigma^2)/sigma^2 with sigma read from GUI
 * @return  the calculated average
 */
#ifdef __CUDACC__
    __host__ __device__
#endif
double CLASSLIB::szave( int c, double a, double x, int n, double z ) const  //{NV}
{
    double f10,f1b,yz;

    if ( fabs(static_cast<int>(z)-z) < eps9 ) z = z + eps9; //  (* avoid division by zero for integer values of z *)

    yz = x / (z+1);

    // (*** gamma(z+n+1)*yz^n/gamma(z+1) ***)
    switch ( n )
    {
    case  10: f10 = (z+10)*(z+9)*(z+8)*(z+7)*(z+6)*(z+5)*(z+4)*(z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   9: f10 = (z+9)*(z+8)*(z+7)*(z+6)*(z+5)*(z+4)*(z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   8: f10 = (z+8)*(z+7)*(z+6)*(z+5)*(z+4)*(z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   7: f10 = (z+7)*(z+6)*(z+5)*(z+4)*(z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   6: f10 = (z+6)*(z+5)*(z+4)*(z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   5: f10 = (z+5)*(z+4)*(z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   4: f10 = (z+4)*(z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   3: f10 = (z+3)*(z+2)*(z+1)*pow(yz,n); break;
    case   2: f10 = (z+2)*(z+1)*pow(yz,n); break;
    case   1: f10 = (z+1)*pow(yz,n); break;
    case   0: f10 = pow(yz,n); break;
    case  -1: f10 = pow(yz,n)/z; break;
    case  -2: f10 = pow(yz,n)/((z-1)*z); break;
    case  -3: f10 = pow(yz,n)/((z-2)*(z-1)*z); break;
    case  -4: f10 = pow(yz,n)/((z-3)*(z-2)*(z-1)*z); break;
    case  -5: f10 = pow(yz,n)/((z-4)*(z-3)*(z-2)*(z-1)*z); break;
    case  -6: f10 = pow(yz,n)/((z-5)*(z-4)*(z-3)*(z-2)*(z-1)*z); break;
    case  -7: f10 = pow(yz,n)/((z-6)*(z-5)*(z-4)*(z-3)*(z-2)*(z-1)*z); break;
    case  -8: f10 = pow(yz,n)/((z-7)*(z-6)*(z-5)*(z-4)*(z-3)*(z-2)*(z-1)*z); break;
    case  -9: f10 = pow(yz,n)/((z-8)*(z-7)*(z-6)*(z-5)*(z-4)*(z-3)*(z-2)*(z-1)*z); break;
    case -10: f10 = pow(yz,n)/((z-9)*(z-8)*(z-7)*(z-6)*(z-5)*(z-4)*(z-3)*(z-2)*(z-1)*z); break;
    default:
        if ( z > 1000 )                                     // Neue Version
            f10 = gammaratio(n+1,1,z) * pow(yz,n);          // Neue Version
        else                                                // Neue Version
            f10 = gamma(z+n+1) * pow(yz,n) / (gamma(z+1));  // Neue Version
        break;
    } // switch n

    // (* **** <x^n> **** *)
    switch ( c )
    {
    case 1:
        return f10;

    case 2: // (* **** <cos(ax) x^n> **** *)
        f1b = cos((z+n+1)*atan(a*yz))*pow(1+a*a*yz*yz,-(z+n+1)/2.);
        return f10*f1b;

    case 3: // (* **** <sin(ax) x^n> **** *)
        f1b = sin((z+n+1)*atan(a*yz))*pow(1+a*a*yz*yz,-(z+n+1)/2.);
        return f10*f1b;

    case 4: // (* **** <cos(ax)^2 x^n> **** *)
        f1b = 1+cos((z+n+1)*atan(2*a*yz))*pow(1+4*a*a*yz*yz,-(z+n+1)/2.);
        return f10*f1b/2.;

    case 5: // (* **** <sin(ax)^2 x^n> **** *)
        f1b = 1-cos((z+n+1)*atan(2*a*yz))*pow(1+4*a*a*yz*yz,-(z+n+1)/2.);
        return f10*f1b/2.;

    case 6: // (* **** <cos(ax) sin(ax) x^n> **** *)
        f1b = sin((z+n+1)*atan(2*a*yz))*pow(1+4*a*a*yz*yz,-(z+n+1)/2.);
        return f10*f1b/2.;
    }
    return 0.0;
}
#endif // USE_szave



#ifdef USE_psphere
#ifdef __CUDACC__
__host__
#endif
void CLASSLIB::psphere_init()       //{NV} zur Optimierung
{
    // global init to avoid same calculations every time
    params.psphere_r = /*(1. + 6.*sqr(params.sigma)) * */ params.radius; // {NV}
    params.psphere_z = (1. - sqr(params.sigma)) / sqr(params.sigma);
}

/**
 * @brief CLASSLIB::psphere - Polysphere calculation
 * @param radius = CLASSLIB::radius
 * @param sigma  = CLASSLIB::sigma
 * @param q      = Local var
 * @return value
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::psphere( double q ) const      // {NV} NEU function psphere(r,sigma,q: extended): extended;
{
    int dim = 3;    // Auswahl: 1=disk, 2=cylinder, 3=sphere

    if ( params.sigma < 0.02 )
        return 9. * sqr(sin(q*params.psphere_r)-q*params.psphere_r*cos(q*params.psphere_r))*exp(-6.*log(q*params.psphere_r));

    double x1 = q * params.psphere_r;
    double x1z = x1 / ( 2. * (params.psphere_z+1.) );
    double x12z = x1z * x1z;

    double lim = 18. * exp(-5. * params.sigma);

    const int nmax = 300;
    const double delc = 0.0001;

    double F12sez,oldF12sez,F12as1z,F12as2z,F12as3z;
    double b1s,v,ee0,ee1,lim1,arg11,nen11,arg12,nen12,arg13,nen13;
    double pz2v,pz2v1,pz2v2,preg1,del;
    double z12v[nmax+1], b1sv[nmax+1], fkv[nmax+1], sum12[nmax+1]; // Array[0..nmax] of extended;

    switch ( dim )
    {
    case 1:     //(* disk *)
        b1s = 3. / 2.;
        v = -1.;
        ee0 = 1.;
        ee1 = 0.;
        //gb1s = sqrt(M_PI) / 2.;
        preg1 = 1. / 2.;
        pz2v = 1./(params.psphere_z * (params.psphere_z - 1.));
        pz2v1 = pz2v / (params.psphere_z - 2.);
        pz2v2 = pz2v1 / (params.psphere_z - 3.);
        lim1 = lim * 1.2;           //(* for term 1 *)
        break;
    case 2:     //(* cylinder *)
        b1s = 2.;
        v = -3. / 2.;
        ee0 = 1.;
        ee1 = -9. / 16.;
        //gb1s = 1.;
        preg1 = 1./sqrt(M_PI);
        pz2v = 1. / (params.psphere_z * (params.psphere_z-1.) * (params.psphere_z-2.));
        pz2v1 = pz2v / (params.psphere_z-3.);
        pz2v2 = pz2v1 / (params.psphere_z-4.);
        lim1 = lim * 1.2;           //(* for term 1 *)
        break;
    case 3:     //(* sphere *)
        b1s = 5. / 2.;
        v = -2.;
        ee0 = 1.;
        ee1 = -1.;
        //gb1s = 3.*sqrt(M_PI)/4.;
        preg1 = 3./4.;
        pz2v = 1./(params.psphere_z*(params.psphere_z-1.)*(params.psphere_z-2.)*(params.psphere_z-3.));
        pz2v1 = pz2v/(params.psphere_z-4.);
        pz2v2 = pz2v1/(params.psphere_z-5.);
        lim1 = lim * 1.2;           //(* for term 1 *)
        break;
    } // switch dim

    if ( x1 < lim1 )
    {
        z12v[0] = 1.;
        b1sv[0] = 1.;
        fkv[0] = 1.;
        F12sez = 1.0;
        oldF12sez = 1.0;
        for ( int n=1; n<=nmax; n++ )
        {
            z12v[n] = z12v[n-1] * ((params.psphere_z+1)-2+2*n)*((params.psphere_z+1)-1+2*n);
            b1sv[n] = b1sv[n-1] * (b1s-1+n);
            fkv[n]  = fkv[n-1] * n;
            sum12[n] = 0;
            for ( int m=0; m<=n; m++ )
                sum12[n] += 1./(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);
            F12sez += pow(-x12z,n)*z12v[n]*sum12[n];
            del = fabs((F12sez-oldF12sez)/F12sez);
            if ( del < delc ) break;
            oldF12sez = F12sez;
        }
        return F12sez;
    } // if ( x1 < lim1 )

    //(*** term #1 ***)
    // if ( x1 >= lim1 ) then begin
    arg11 = (params.psphere_z+2.*v+1.)*atan(4.*x1z);
    nen11 = pow(1.+16.*x1z*x1z,(params.psphere_z+2*v+1)/2.);
    arg12 = (params.psphere_z+2*v)*atan(4.*x1z);
    nen12 = pow(1.+16.*x1z*x1z,(params.psphere_z+2*v)/2.);
    arg13 = (params.psphere_z+2.*v-1.)*atan(4.*x1z);
    nen13 = pow(1.+16.*x1z*x1z,(params.psphere_z+2.*v-1)/2.);

    F12as1z = ee0*ee0*pz2v*(1.+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);
    F12as2z = 2.*ee0*ee1*(1./(2.*x1z))*pz2v1*(cos(M_PI*(2.*v-1.)/2.)*cos(arg12)/nen12-sin(M_PI*(2.*v-1.)/2.)*sin(arg12)/nen12);
    F12as3z = ee1*ee1*(1./(4.*x1z*x1z))*pz2v2*(1.+cos(M_PI*(v-1.))*cos(arg13)/nen13-sin(M_PI*(v-1.))*sin(arg13)/nen13);
    //F12asz = preg1*preg1*pow(x1z,2.*v)*(1./2.)*(F12as1z+F12as2z+F12as3z);
    return /*F12asz0=*/ preg1*preg1*pow(x1z,2.*v)*(1./2.)*(F12as1z+F12as2z+F12as3z);
}
#endif // USE_psphere



#ifdef USE_psphered
//{NV}-Upq.pas-1475
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::psphered(double r, double sigma, int dim, double q) const
{
    const int nmax = 300;
    const double delc = 0.0001;

    double z;
    double F12as1z,F12as2z,F12as3z,zz;
    double b1s,v,ee0,ee1,x1,x1z,x12z,lim1,arg11,nen11,arg12,nen12,arg13,nen13;
    double pz2v,pz2v1,pz2v2,preg1,del,com;

    z = (1-sqr(sigma))/sqr(sigma);
    zz = z;
    z = z+2*dim;          /* z-average */

    x1 = q*r;
    x1z = x1/(2*(zz+1));
    x12z = x1z*x1z;

    lim1 = 18*exp(-5*sigma);

    switch ( dim )
    {
    case 1: /* disk */
        b1s = 3/2.0;
        v = -1;
        ee0 = 1;
        ee1 = 0;
        //gb1s = sqrt(M_PI)/2.0;
        preg1 = 1/2.0;
        pz2v = 1/(z*(z-1));
        pz2v1 = pz2v/(z-2);
        pz2v2 = pz2v1/(z-3);
        lim1 = lim1*1.2;           /* for term 1 */
        break;
    case 2: /* cylinder */
        b1s = 2;
        v = -3/2.0;
        ee0 = 1;
        ee1 = -9/16.0;
        //gb1s = 1;
        preg1 = 1/sqrt(M_PI);
        pz2v = 1/(z*(z-1)*(z-2));
        pz2v1 = pz2v/(z-3);
        pz2v2 = pz2v1/(z-4);
        lim1 = lim1*1.2;           /* for term 1 */
        break;
    case 3: /* sphere */
        b1s = 5/2.0;
        v = -2;
        ee0 = 1;
        ee1 = -1;
        //gb1s = 3*sqrt(M_PI)/4.0;
        preg1 = 3/4.0;
        pz2v = 1/(z*(z-1)*(z-2)*(z-3));
        pz2v1 = pz2v/(z-4);
        pz2v2 = pz2v1/(z-5);
        lim1 = lim1*1.1;           /* for term 1 */
        break;
    } // switch dim

    if ( x1<=lim1 )
    {
        double z12v = 1;
        double b1sv = 1;
        double fkv = 1;
        double fk2v = 1;
        double F12sez = 1.0;
        double oldF12sez = 1.0;
        for ( int n=1; n<=nmax; n++ )
        {
            z12v = z12v*((z+1)-2+2*n)*((z+1)-1+2*n);
            b1sv = b1sv*(b1s-1+n);
            fkv = fkv*n;
            switch ( dim )
            {
            case 1:
                com = pow(4,n)/(n+1.0);
                break;
            case 2:
                fk2v = fk2v*(2*n-1)*(2*n);
                com = 4*(n+1/2.0)*fk2v/((n+2)*(n+1)*fkv*fkv);
                break;
            case 3:
                com = pow(4.0,n)*6.0/((n+2.0)*(n+3.0));
                break;
            } // switch dim
            F12sez = F12sez+pow(-x12z,n)*z12v*com/(b1sv*fkv);
            del = fabs((F12sez-oldF12sez)/F12sez);
            if ( del<delc ) break; // goto 101;
            oldF12sez = F12sez;
        }
        //101:
        return F12sez;
    }

    /*** term #1 ***/
    else // if ( x1>lim1 )
    {
        arg11 = (z+2*v+1)*atan(4*x1z);
        nen11 = pow(1+16*x1z*x1z,(z+2*v+1)/2.0);
        arg12 = (z+2*v)*atan(4*x1z);
        nen12 = pow(1+16*x1z*x1z,(z+2*v)/2.0);
        arg13 = (z+2*v-1)*atan(4*x1z);
        nen13 = pow(1+16*x1z*x1z,(z+2*v-1)/2.0);

        F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);
        F12as2z = 2*ee0*ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);
        F12as3z = ee1*ee1*(1/(4*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);
        return preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);
    }
}
#endif // USE_psphered



#ifdef USE_pspheredf
/* ************************** Polysphere ****************************/
//{NV}-Upq.pas-1610
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::pspheredf(double r, double sigma, double d, double q) const
{

    const int nmax = 300;
    const double delc = 0.0001;

    double /*x,*/z; //,pq,p1,p2,p3,fct;
    int    n; //,m;
    double F12sez,oldF12sez,F12as1z,F12as2z,/*F12as3z,*/F12asz0,/*F12asz,*/dim,zz;
    double b1s,v,ee0,ee1,x1,x1z,x12z,lim,lim1,arg11,nen11,arg12,nen12; //,arg13,nen13;
    double /*gb1s,*/pz2v,pz2v1,/*pz2v2,*/preg1,del;
    double z12v[nmax+1],b1sv[nmax+1],fkv[nmax+1]; //,sum12[nmax+1]; // : Array[0..nmax] of extended;;

    /*begin*/
    z = (1-sqr(sigma))/sqr(sigma);
    dim = d;
    zz = z;
    z = z+2*dim;        /* z-average */

    /* conversion from input n-average to calculation w-average */
    /* r:=(1+6*sigma*sigma)*r; */
    /*x:=q*r;
       p1:=szave(5,1,x,-6,z);
       p2:=szave(6,1,x,-5,z);
       p3:=szave(4,1,x,-4,z);
       fct:=9*(p1-2*p2+p3);
       psphere:=fct;  */

    x1 = q*r;
    x1z = x1/(2*(zz+1));
    x12z = x1z*x1z;

    lim = 18*exp(-5*sigma);

    if ( dim==1 )
    {     /* disk */
        b1s = 3/2.0;
        v = -1;
        ee0 = 1;
        ee1 = 0;
        //gb1s = sqrt(M_PI)/2.0;
        preg1 = 1/2.0;
        pz2v = 1/z;
        pz2v1 = pz2v/(z-1);
        //pz2v2 = pz2v1/(z-2);
        lim1 = lim*1.2;           /* for term 1 */
    }
    if ( dim==2 )
    {     /* cylinder */
        b1s = 2;
        v = -3/2.0;
        ee0 = 1;
        ee1 = -9/16.0;
        //gb1s = 1;
        preg1 = 1/sqrt(M_PI);
        pz2v = gamma(z-1/2.0)/gamma(z+1);
        pz2v1 = pz2v1/(z-3/2.0);
        //pz2v2 = pz2v1/(z-5/2.0);
        lim1 = lim*1.2;           /* for term 1 */
    }
    if ( dim==3 )
    {     /* sphere */
        b1s = 5/2.0;
        v = -2;
        ee0 = 1;
        ee1 = -1;
        //gb1s = 3*sqrt(M_PI)/4.0;
        preg1 = 3/4.0;
        pz2v = 1/(z*(z-1));
        pz2v1 = pz2v/(z-2);
        //pz2v2 = pz2v1/(z-3);
        lim1 = lim*1.2;           /* for term 1 */
    }

    if ( x1<lim1 )
    {
        z12v[0] = 1;
        b1sv[0] = 1;
        fkv[0] = 1;
        F12sez = 1.0;
        oldF12sez = 1.0;
        for ( n=1; n<=nmax; n++ )
        {
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);
            b1sv[n] = b1sv[n-1]*(b1s-1+n);
            fkv[n] = fkv[n-1]*n;
            F12sez = F12sez+pow(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]);
            del = fabs((F12sez-oldF12sez)/F12sez);
            if ( del<delc ) break; // goto 101;
            oldF12sez = F12sez;
        }
        //101:
        return F12sez*F12sez;
    }

    /*** term #1 ***/
    else // if ( x1>=lim1 )
    {
        arg11 = (z+v+1)*atan(2*x1z);
        nen11 = pow(1+4*x1z*x1z,(z+v+1)/2.0);
        arg12 = (z+v)*atan(2*x1z);
        nen12 = pow(1+4*x1z*x1z,(z+v)/2.0);

        F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);
        F12as2z = ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);
        F12asz0 = preg1*pow(x1z,v)*(F12as1z+F12as2z);
        return F12asz0*F12asz0;
    }
}
#endif // USE_pspheredf



#ifdef USE_f2dschulz
/**
 * @brief SasCalculation::f2dschulz
 * @param d
 * @param r
 * @param sigma
 * @param q
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::f2dschulz(int d, double r, double sigma, double q) const       //{NV}
{
    double z;
    double hgeo2;

    z  = (1-sqr(sigma))/sqr(sigma);

    if ( d == 1 )
        hgeo2 = szave(3,1,q*r,-1,z);
    else if ( d == 2 )
        hgeo2 = f2dhyper(2,0.001,r,sigma,q);
    else //if ( d == 3 )
        hgeo2 = 3*(szave(3,1,q*r,-3,z)-szave(2,1,q*r,-2,z));
    return hgeo2*hgeo2;
}
#endif // USE_f2dschulz



#ifdef USE_f2dhyper
//(* ************************ d-dimensional <F(q)> average *************** *)
/**
 * @brief SasCalculation::f2dhyper
 * @param d
 * @param alf
 * @param r
 * @param sigma
 * @param q
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::f2dhyper(double d, double alf, double r, double sigma, double q) const     //{NV}
{
    const int maxsum = 250;         // (* number of summations *)

    int n,sig;
    double cc[6], dd[6], ee[6]; // array[0..5] of real;
    double a1,a2,a3,b1,b2,gam,z,x,y,
           hgeo,oldhgeo,
           s1,s1a,s2,s3,t1,t2,t3a,t3,f2;

    z=(1.-sqr(sigma))/sqr(sigma);
    a1=(d-alf)/2.;
    a2=(z+1)/2.;
    a3=(z+2)/2.;
    b1=(d+2-alf)/2.;
    b2=d/2.;
    x=sqr(q*r/(z+1));
    gam=gamma(b1)*gamma(b2)/(gamma(a1)*gamma(a2)*gamma(a3));
    y=q*r/(z+1);

    //(* **** Asymptotic expansion **** *)
    if ( q*r > 10 )
    {
        s1=a1-b1-b2;
        s1a=(1./2.)+s1;
        s2=a1*a1-b1*b1-b2*b2;
        s3=a1*a1*a1-b1*b1*b1-b2*b2*b2;
        t1=(1./3.)+s1-s1a*s1a-2*(-(1./3.)-s1+s2);
        t2=t1*t1/8.+(1./6.)*((3./2.)*s1a*s1a-s1a*s1a*s1a-(1./2.)*s1a+4.*((1./2.)*s1-(3./2.)*s2+s3));
        t3a=gamma(a1);
        t3=t3a/(gamma(b1-a1)*gamma(b2-a1));
        cc[1]=exp(-s1a*log(4.)/2.)*exp(-log(M_PI)/2.);
        cc[2]=t1*cc[1]/2.;
        cc[3]=t2*cc[1];
        cc[4]=exp(a1*log(4.))*t3;
        cc[5]=exp((a1+1)*log(4.))*a1*(1+a1-b1)*(1+a1-b2)*t3;
        dd[1]=((1./2.)+s1)*M_PI/2.;
        dd[2]=((3./2.)+s1)*M_PI/2.;
        dd[3]=((5./2.)+s1)*M_PI/2.;
        dd[4]=-1000;
        dd[5]=-1000;
        ee[1]=s1+1./2.;
        ee[2]=s1-1./2.;
        ee[3]=s1-3./2.;
        ee[4]=-2*a1;
        ee[5]=-2*a1-2;

        f2=0;
        for ( n=1; n<=5; n++ )
            f2 += cc[n]*cosav(dd[n],-1000,ee[n],0,y,z);
        return f2*gamma(b1)*gamma(b2)/gamma(a1);
    }

    //(* **** sum over hypergeometric functions **** *)
    hgeo=0;
    for ( n=0; n<=maxsum; n++ )
    {
        oldhgeo=hgeo;
        if ( (n % 2) == 0 ) sig=1;
        else sig=-1;
        hgeo += gamma(a1+n)*gamma(a2+n)*gamma(a3+n)*sig*exp(n*log(x))/(gamma(b1+n)*gamma(b2+n)*gamma(n+1));
        if ( fabs(hgeo-oldhgeo) < eps6) break;
    }
    return gam*hgeo;
}
#endif // USE_f2dhyper



#ifdef USE_gamma
//(************************** Gamma function ***************************)
/**
 * @brief CLASSLIB::gamma
 * @param z
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::gamma(double z) const          //{NV}
{
    // Die einzelnen Zahlen werden in Arrays umkopiert. Daher hier direkt die Arrys fertigstellen.
    // Dabei ist der erste Wert immer 0, da die Pascal-Arrays bei 1 beginnen und ich nicht alle Indizes umsetzen will.
    const double bm[] = { 0.0,
                          -0.577191652, 0.988205891, -0.897056937, 0.918206857, -0.756704078, 0.482199394, -0.193527818, 0.035868343 };
    const double cm[] = {  0.0,
                           1.0000000000000000,  0.5772156649015329, -0.6558780715202538,  -0.0420026350340952,  0.1665386113822915,
                          -0.0421977345555443, -0.0096219715278770,  0.00721894322466630, -0.0011651675918591, -0.0002152416741149,
                           0.0001280502823882, -0.0000201348547807, -0.0000012504934821,   0.0000011330272320, -0.0000002056338417,
                           0.0000000061160950,  0.0000000050020075, -0.0000000011812746,   0.0000000001043427,  0.0000000000077823,
                          -0.0000000000036968,  0.0000000000005100, -0.0000000000000206,  -0.0000000000000054,  0.0000000000000014,
                           0.0000000000000001 };

    int i, di;
    double x, fak, y, f, fct, absz;
    bool negative;

    if ( fabs(z) < eps9 ) z = z + eps9;
    if ( z > 0 )
        negative = false;
    else //if ( z < 0 )
    {
        negative = true;
        if ( fabs(frac(z)) < eps9 ) z = z+eps9;
    }
    absz = fabs(z);

    if ( absz > 0.0 && absz <= 1 )
    {
        fct = 0;
        for ( i=1; i<=26; i++ ) fct += cm[i]*pow(absz,i);
        if ( negative ) return -M_PI*fct/(absz*sin(M_PI*absz));
        return 1./fct;
    }
    if ( absz > 1 && absz <= 1000 )
    {
        if ( fabs(frac(absz)) < eps9 )  // (* integer z *)
        {
            fct = 1;
            for ( i=1; i<=trunc(absz)-1; i++ ) fct = fct*i;
            if ( negative ) return -M_PI/(absz*fct*sin(M_PI*absz));
            return fct;
        }
        y = absz-1;
        di = trunc(y);
        x = absz-trunc(absz);
        fak = y;
        if ( di == 0 ) fak = 1;
        else for ( i=1; i<di; i++ ) fak = fak*(y-i);
        f = 1;
        for ( i=1; i<=8; i++ ) f = f+bm[i]*pow(x,i);
        if ( negative ) return -M_PI/(absz*fak*f*sin(M_PI*absz));
        return fak*f;
    }
    //if ( absz > 1000 )
    //{
        fct = exp(-absz)*exp((absz-0.5)*log(absz))*sqrt(2*M_PI);
        if ( negative ) fct = -M_PI/(absz*fct*sin(M_PI*absz));
        return fct;
    //}
    //return 0.0;
}
#endif // USE_gamma



#ifdef USE_cosav
//(************* function <cos(di+y)cos(dj+y)y^(ei+ej)> **************)
/**
 * @brief CLASSLIB::cosav
 * @param di
 * @param dj
 * @param ei
 * @param ej
 * @param y
 * @param z
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::cosav(double di, double dj, double ei, double ej, double y, double z) const        //{NV}
{
    double nav,p2a,p2b,p2c,ava,avb,avc,f2a,f2b,f2c,yz,pyz;

    if ( fabs(int(z)-z) < eps9 ) z=z+eps9;   //(* avoid division by zero for integer values of z *)

    nav=ei+ej;
    yz=y;
    pyz=pow(yz,nav);

    // Vergleiche, ob di oder dj == -1000 oder nicht sind. Ich stelle diese Gruppen etwas um,
    //  so dass die double-Vergleiche sicherer sind.
    if ( di > -1000 )
    {
        if ( dj > -1000 )
        {   // if ( di != -1000 && dj != -1000 )
            if ( z > 1000 )
            {
                p2a=gammaratio(nav+1,1,z)*pyz/2.;
                p2b=gammaratio(nav+2,1,z)*pyz/((z+nav+1)*2);
                p2c=gammaratio(nav+3,1,z)*pyz/((z+nav+1)*(z+nav+2)*2);
            }
            else
            {
                p2a=gamma(z+nav+1)*pyz/(gamma(z+1)*2);
                p2b=gamma(z+nav+2)*pyz/(gamma(z+1)*(z+nav+1)*2);
                p2c=gamma(z+nav+3)*pyz/(gamma(z+1)*(z+nav+1)*(z+nav+2)*2);
            }
            //ava:=(1+cos((z+nav+1)*arctan(2*yz))/exp((z+nav+1)*log(1+4*yz*yz)/2));
            ava=(1+cos((z+nav+1)*atan(2*yz))*pow(1+4*yz*yz,-(z+nav+1)/2.));
            //avb:=   sin((z+nav+1)*arctan(2*yz))/exp((z+nav+1)*log(1+4*yz*yz)/2);
            avb=   sin((z+nav+1)*atan(2*yz))*pow(1+4*yz*yz,-(z+nav+1)/2.);
            //avc:=(1-cos((z+nav+1)*arctan(2*yz))/exp((z+nav+1)*log(1+4*yz*yz)/2));
            avc=(1-cos((z+nav+1)*atan(2*yz))*pow(1+4*yz*yz,-(z+nav+1)/2.));
            f2a=cos(di)*cos(dj)*p2a*ava;
            f2b=sin(di+dj)*p2b*avb;
            f2c=sin(di)*sin(dj)*p2c*avc;
            return f2a-f2b+f2c;
        }
        else
        {   // if ( di != -1000 && dj == -1000 )
            if ( z > 1000 )
            {
                p2a=gammaratio(nav+1,1,z)*pyz;
                p2b=gammaratio(nav+2,1,z)*pyz/((z+nav+1));
            }
            else
            {
                p2a=gamma(z+nav+1)*pyz/(gamma(z+1));
                p2b=gamma(z+nav+2)*pyz/(gamma(z+1)*(z+nav+1));
            }
            //ava:=cos((z+nav+1)*arctan(yz))/exp((z+nav+1)*log(1+yz*yz)/2);
            ava=cos((z+nav+1)*atan(yz))*pow(1+yz*yz,-(z+nav+1)/2.);
            //avb:=sin((z+nav+1)*arctan(yz))/exp((z+nav+1)*log(1+yz*yz)/2);
            avb=sin((z+nav+1)*atan(yz))*pow(1+yz*yz,-(z+nav+1)/2.);
            return cos(di)*p2a*ava-sin(di)*p2b*avb;
        }
    }
    else
    {
        if ( dj > -1000 )
        {   // if ( di == -1000 && dj != -1000 )
            if ( z > 1000 )
            {
                p2a=gammaratio(nav+1,1,z)*pyz;
                p2b=gammaratio(nav+2,1,z)*pyz/((z+nav+1));
            }
            else
            {
                p2a=gamma(z+nav+1)*pyz/(gamma(z+1));
                p2b=gamma(z+nav+2)*pyz/(gamma(z+1)*(z+nav+1));
            }
            //ava:=cos((z+nav+1)*arctan(yz))/exp((z+nav+1)*log(1+yz*yz)/2);
            ava=cos((z+nav+1)*atan(yz))*pow(1+yz*yz,-(z+nav+1)/2.);
            //avb:=sin((z+nav+1)*arctan(yz))/exp((z+nav+1)*log(1+yz*yz)/2);
            avb=sin((z+nav+1)*atan(yz))*pow(1+yz*yz,-(z+nav+1)/2.);
            return cos(dj)*p2a*ava-sin(dj)*p2b*avb;
        }
        else
        {   // if ( di == -1000 && dj == -1000 )
            if ( z > 1000 ) return gammaratio(nav+1,1,z)*pyz;
            return gamma(z+nav+1)*pyz/gamma(z+1);
        }
    }
    return 0;
}
#endif // USE_cosav



#ifdef USE_gammaratio
//(********************* gamma(z+a)/gamma(z+b) for large z *******************)
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::gammaratio(double a, double b, double x) const         //{NV}
{
    double c1 = pow(1+a/x,x+a-0.5);
    double c2 = pow(1+b/x,x+b-0.5);
    return exp(-a+b)*exp((a-b)*log(x))*c1/c2;
}
#endif // USE_gammaratio



#ifdef USE_pqcoreshell
/**
 * @brief SasCalculation::pqcoreshell
 * @param rho1  = 1.0
 * @param rho2  = global var rho => params.rho
 * @param p1    = global var p1 => params.p1
 * @param p2    = 1.0
 * @param alf1  = 0.001
 * @param alf2  = local var
 * @param rn    = global var radiusi => params.radiusi
 * @param d     = local var
 * @param sigma = global var sigma => params.sigma
 * @param q     = local var
 * @return
 *
// BCC:  pqcoreshell(1.0,rho,p1,1.0,0.001,0.0001 ,radiusi,3,sigma,q);   TODO: dieses überprüfen!
// BCC:  pqcoreshell(1.0,rho,p1,1.0,0.001,alphash,radiusi,3,sigma,q);
// FCC?: pqcoreshell(1.0,rho,p1,1.0,0.001,0.0001 ,rm     ,3,sig  ,q);
// FCC?: pqcoreshell(1.0,rho,p1,1.0,0.001,alf    ,rm     ,3,sig  ,q);
// FD3M: pqcoreshell(1.0,rho,p1,1.0,0.001,0.0001 ,radiusi,3,sigma,q);
// FD3M: pqcoreshell(1.0,rho,p1,1.0,0.001,alphash,radiusi,3,sigma,q);
// HCP : pqcoreshell(1.0,rho,p1,1.0,0.001,0.0001 ,radiusi,3,sigma,q);
// HCP : pqcoreshell(1.0,rho,p1,1.0,0.001,alphash,radiusi,3,sigma,q);
// SC  : pqcoreshell(1.0,rho,p1,1.0,0.001,0.0001 ,radiusi,3,sigma,q);
// SC  : pqcoreshell(1.0,rho,p1,1.0,0.001,alphash,radiusi,3,sigma,q);

TODO:
   function pqcoreshell(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended): extended;  <== hier, Rest ist neu ...
   function pqcoreshellf(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended): extended;
   function pqcoreshellin(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended): extended;
   function pqcoreshellinf(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended): extended;
   function pqcoreshelld3(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q,slit: extended): extended;
   function pqcoreshelld2(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended): extended;
   function pqcoreshelld1(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended): extended;
   function pqcoreshell1(alf,rn,d,sigma,q: extended): extended;
   function pqcoreshell2(rn,sigma,q: extended): extended;
   function pqcoreshellhves(rhoi,rho2,p,pp,alf1,alf2,rn,rv,d,sigma,q: extended): extended;
   function polypqcoreshellhves(rhoi,rho2,p,sigmav,alf1,alf2,rn,rv,d,sigma,q: extended): extended;
   function pqcoreshell_c(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended;
                      car1,car2,car3: array of extended): extended;
   function pqcoreshellin_c(rho1,rho2,p1,p2,alf1,alf2,rn,d,sigma,q: extended;
                      car1,car2,car3,car4,car5,car6: array of extended): extended;

 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::pqcoreshell( double /*alf2*/, double d, double q ) const         //{NV} NEU plus Param d
{
    // Constant parameters in all calls
    //- const double rho1 = 1.0;
    //- const double p2   = 1.0;
    //- const double alf1 = 0.001; //not used

    // Other Parameters
    const double rn   = params.radiusi;
    const double rho2 = params.rho;

    const int nmax = 300;
    const double delc = 0.0001;

    //- double cpq[7], fpq[7]; // array[1..6] of extended;
    //- double nenner,zaehler;
    double c1,c2,c3/*,c4,f1,f2,f3,fq,x,y*/,z/*,yz,da*/;
    int /*i,j,*/n,m,dim;
    double /*zl,x1,xz,*/x1z,x12z,b1s,v,e0,e1,/*gb1s,*/pz2v,pz2v1,pz2v2,preg1,/*lim,*/lim1,lim2,lim3,vv,del;
    double xrad,xradp,x2z,x22z,com,zz,epsi;
    double F12as1z,F12as2z,F12as3z,F12asz,F32as1z,F32as2z,F32as3z,F32asz;
    double arg11,nen11,arg12,nen12,arg13,nen13,arg31,nen31,arg32,nen32,arg33,nen33;
    double F12sez1,oldF12sez1,F121,F12sez2,oldF12sez2,F122/*,F122a,F122b*/,F12sez3,oldF12sez3,F123;
    double /*F121as1z,F121as2z,F121as3z,F121as4z,F121asz,*/F122as1z,F122as2z,F122as3z,F122as4z,F122asz;
    //- double F123as1z,F123as2z,F123as3z,F123as4z,F123asz,F124as1z,F124as2z,F124as3z,F124as4z,F124asz;
    double /*arg0,nen0,arg,nen,arg1,nen1,arg2,nen2,*/xijm,arglmz,nenlmz,xijp,arglpz,nenlpz;
    double z12v[nmax+1], b1sv[nmax+1], fkv[nmax+1], fk2v[nmax+1], /*sum121[nmax+1],*/ sum122[nmax+1], /*sum123[nmax+1],*/ eps1[nmax+1]; // Array[0..nmax] of extended;
    double cc[11], ac[11]; // Array[1..10] of extended;
    //- int bin[nmax+1][nmax+1]; // Array[0..nmax,0..nmax] of integer;
    bool ex;
    double x1zz,x2zz/*,ffz*/,vv3,pz4,pz5,pz6/*,F12*/,lima,limb/*,ipq*/,binsum;
    double arga,nena,argbm,nenbm,argbp,nenbp,argc,nenc,argd,nend,argep,nenep,argem,nenem;
    double argf,nenf,arggp,nengp,arggm,nengm,argh,nenh,/*argim,nenim,argip,nenip,*/argj,nenj;

    dim = round(d);
    z   = (1-sqr(params.sigma))/sqr(params.sigma);
    zz  = z;
    z   = z+2.*dim;           //(* z-average *)
    xrad  = q*rn;
    xradp = q*params.p1*rn;
    x1z   = q*params.p1*rn/(2.*(zz+1));
    x1zz  = 2.*x1z;
    x12z  = x1z*x1z;
    x2z   = q*rn/(2.*(zz+1));
    x2zz  = 2.*x2z;
    x22z  = x2z*x2z;

    epsi = 1.0;

    //lim  = 18.*exp(-5.*params.sigma);
    lima = 18.*exp(-5.*params.sigma);
    limb = 18.*exp(-2.5*params.sigma);

    switch ( dim )
    {
    case 1:     //(* disk *)
        ex = false;
        b1s = 3./2.;
        v = -1;
        e0 = 1;
        e1 = 0;
        //gb1s = sqrt(M_PI)/2.;
        preg1 = 1./2.;
        pz2v = 1./(z*(z-1));
        pz2v1 = pz2v/(z-2);
        pz2v2 = pz2v1/(z-3);
        lim1 = lima*1.2;      //(* ok *)
        lim2 = lima*1.7;
        lim3 = lima*1.2;
        limb = limb*1.1;
        break;
    case 2:     //(* cylinder *)
        ex = false;
        b1s = 2;
        v = -3./2.;
        e0 = 1;
        e1 = -9./16.;
        //gb1s = 1;
        preg1 = 1./sqrt(M_PI);
        pz2v = 1./(z*(z-1)*(z-2));
        pz2v1 = pz2v/(z-3);
        pz2v2 = pz2v1/(z-4);
        lim1 = lima*1.2;            //(* ok *)
        lim2 = lima*1.5;
        lim3 = lima*1.2;
        limb = limb*1.1;
        break;
    case 3:     //(* sphere *)
        ex = true;
        b1s = 5./2.;
        v = -2;
        e0 = 1;
        e1 = -1;
        //gb1s = 3.*sqrt(M_PI)/4.;
        preg1 = 3./4.;
        pz2v = 1./(z*(z-1)*(z-2)*(z-3));
        pz2v1 = pz2v/(z-4);
        pz2v2 = pz2v1/(z-5);
        lim1 = lima*1.2;             //(* ok *)
        lim2 = lima*1.7;
        lim3 = lima*1.2;
        limb = limb*1.2;
        break;
    } // switch dim

    pz4 = 1./(z*(z-1)*(z-2)*(z-3));
    pz5 = pz4/(z-4);
    pz6 = pz5/(z-5);

    c1 = sqr(1-rho2)*pow(params.p1,2*dim);
    c2 = 2.*rho2*(1-rho2)*pow(params.p1,dim);
    c3 = rho2*rho2;
    vv = c1+c2+c3;
    vv3 = sqr((1-rho2)*pow(params.p1,3)+rho2);

    //(* series expansion: provides F121, F122, F123 *)
    if ( xradp <= lim1 )
    {
        z12v[0] = 1;
        b1sv[0] = 1;
        fkv[0] = 1;
        fk2v[0] = 1;
        eps1[0] = 1;
        F12sez1 = 1.0;
        oldF12sez1 = 1.0;
        for ( n=1; n<=nmax; n++ )
        {
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);
            b1sv[n] = b1sv[n-1]*(b1s-1+n);
            fkv[n]  = fkv[n-1]*n;
            eps1[n] = eps1[n-1]*(epsi*epsi-1);
            switch ( dim )
            {
            case 1:
                com = pow(4.,n)/(n+1.);
                break;
            case 2:
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);
                com = 4.*(n+1./2.)*fk2v[n]/((n+2.)*(n+1.)*fkv[n]*fkv[n]);
                break;
            case 3:
                com = pow(4.,n)*6/((n+2.)*(n+3.));
                break;
            }
            binsum = 0.0;
            for ( m=0; m<=n; m++ ) binsum += (fkv[n]/(fkv[m]*fkv[n-m]))*eps1[n-m]/(2.*(n-m)+1);
            F12sez1 += pow(-x12z,n)*z12v[n]*com*binsum/(b1sv[n]*fkv[n]);
            del = fabs((F12sez1-oldF12sez1)/F12sez1);
            if ( del < delc ) break;
            oldF12sez1 = F12sez1;
        }
        F121 = c1*F12sez1;
    } // if ( xradp <= lim1 )

    if ( xrad <= lim2)
    {
        z12v[0] = 1;
        b1sv[0] = 1;
        fkv[0] = 1;
        eps1[0] = 1;
        F12sez2 = 1.0;
        oldF12sez2 = 1.0;
        for ( n=1; n<=nmax; n++ )
        {
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);
            b1sv[n] = b1sv[n-1]*(b1s-1+n);
            fkv[n] = fkv[n-1]*n;
            eps1[n] = eps1[n-1]*(epsi*epsi-1);
            sum122[n] = 0;
            for ( m=0; m<=n; m++ ) sum122[n] += pow(params.p1*params.p1,m)/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);
            binsum = 0.0;
            for ( m=0; m<=n; m++ ) binsum += (fkv[n]/(fkv[m]*fkv[n-m]))*eps1[n-m]/(2.*(n-m)+1.);
            F12sez2 += pow(-x22z,n)*z12v[n]*sum122[n]*binsum;
            del = fabs((F12sez2-oldF12sez2)/F12sez2);
            if ( del < delc ) break;
            oldF12sez2 = F12sez2;
        }
        F122 = c2*F12sez2;
    } // if ( xrad <= lim2)

    if ( xrad <= lim3 )
    {
        z12v[0] = 1;
        b1sv[0] = 1;
        fkv[0] = 1;
        eps1[0] = 1;
        F12sez3 = 1.0;
        oldF12sez3 = 1.0;
        for ( n=1; n<=nmax; n++)
        {
            z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);
            b1sv[n] = b1sv[n-1]*(b1s-1+n);
            fkv[n] = fkv[n-1]*n;
            eps1[n] = eps1[n-1]*(epsi*epsi-1);
            switch ( dim )
            {
            case 1:
                com = pow(4.,n)/(n+1.);
                break;
            case 2:
                fk2v[n] = fk2v[n-1]*(2*n-1)*(2*n);
                com = 4.*(n+1./2.)*fk2v[n]/((n+2.)*(n+1.)*fkv[n]*fkv[n]);
                break;
            case 3:
                com = pow(4.,n)*6./((n+2.)*(n+3.));
                break;
            }
            binsum = 0.0;
            for ( m=0; m<=n; m++ ) binsum += (fkv[n]/(fkv[m]*fkv[n-m]))*eps1[n-m]/(2.*(n-m)+1.);
            F12sez3 += pow(-x22z,n)*z12v[n]*com*binsum/(b1sv[n]*fkv[n]);
            del = fabs((F12sez3-oldF12sez3)/F12sez3);
            if ( del < delc ) break;
            oldF12sez3 = F12sez3;
        }
        F123 = c3*F12sez3;
    } // if ( xrad <= lim3 )

    //(*** asymptotic expansion for d=1,2 ***)
    if ( !ex )
    {
        //(*** F121 integration and asymptote ***)
        if ( epsi == 1 )   //(* sphere *)
        {
            if ( (xradp > lim1) && (xradp < limb) ) //(* numerical integration for d=1,2 *)
            {
                if ( dim == 1 ) qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,12,1,F121);
                if ( dim == 2 ) qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,9,2,F121);
            }
            if ( xradp >= limb )
            {
                arg11 = (z+2.*v+1.)*atan(4.*x1z);
                nen11 = pow(1.+16.*x1z*x1z,(z+2.*v+1)/2.);
                arg12 = (z+2.*v)*atan(4.*x1z);
                nen12 = pow(1.+16.*x1z*x1z,(z+2.*v)/2.);
                arg13 = (z+2.*v-1.)*atan(4.*x1z);
                nen13 = pow(1.+16.*x1z*x1z,(z+2.*v-1.)/2.);

                F12as1z = e0*e0*pz2v*(1.+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);
                F12as2z = 2.*e0*e1*(1./(2.*x1z))*pz2v1*(cos(M_PI*(2.*v-1.)/2.)*cos(arg12)/nen12-sin(M_PI*(2.*v-1.)/2.)*sin(arg12)/nen12);
                F12as3z = e1*e1*(1./(4.*x1z*x1z))*pz2v2*(1.+cos(M_PI*(v-1.))*cos(arg13)/nen13-sin(M_PI*(v-1.))*sin(arg13)/nen13);
                F12asz = preg1*preg1*pow(x1z,2.*v)*(1./2.)*(F12as1z+F12as2z+F12as3z);
                F121 = c1*F12asz;
            }
        }
        else //(* for ellpisoids *)
        {
            if ( xradp > lim1 )
                qrombpi2(epsi,rn,params.p1,rho2,1,1,params.sigma,q,1,1,11,11,dim,F121);
        }

        //(*** F122 integration and asymptote ***)
        if ( epsi == 1 )   //(* sphere *)
        {
            if ( (xrad > lim2) && (xrad < limb) ) //(* numerical integration for d=2 *)
            {
                if ( dim == 1 ) qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,13,1,F122);
                if ( dim == 2 ) qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,10,2,F122);
            }

            if ( xrad > limb )
            {
                xijm = x1z-x2z;
                arglmz = (z+1.)*atan(2.*xijm);
                nenlmz = pow(1.+4.*xijm*xijm,(z+1)/2.);
                xijp = x1z+x2z;
                arglpz = (z+1)*atan(2.*xijp);
                nenlpz = pow(1.+4.*xijp*xijp,(z+1)/2.);
                F122as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+cos(M_PI*v)*cos(arglpz)/nenlpz-sin(M_PI*v)*sin(arglpz)/nenlpz);
                F122as2z = e0*e1*(1./(2.*x2z))*pz2v1*(-sin(arglmz)/nenlmz+cos(M_PI*(2.*v-1.)/2.)*cos(arglpz)/nenlpz-sin(M_PI*(2.*v-1.)/2.)*sin(arglpz)/nenlpz);
                F122as3z = e1*e0*(1./(2.*x1z))*pz2v1*(sin(arglmz)/nenlmz+cos(M_PI*(2.*v-1.)/2.)*cos(arglpz)/nenlpz-sin(M_PI*(2.*v-1.)/2.)*sin(arglpz)/nenlpz);
                F122as4z = e1*e1*(1./(4.*x1z*x2z))*pz2v2*(cos(arglmz)/nenlmz+cos(M_PI*(v-1.))*cos(arglpz)/nenlpz-sin(M_PI*(v-1.))*sin(arglpz)/nenlpz);
                F122asz = (F122as1z+F122as2z+F122as3z+F122as4z);
                F122 = c2*preg1*preg1*pow(x1z,v)*pow(x2z,v)*(1./2.)*F122asz;
            }
        }
        else //(* for ellpisoids *)
        {
            if ( xrad > lim2 )
                qrombpi2(epsi,rn,params.p1,rho2,1,1,params.sigma,q,1,1,11,12,dim,F122);
        }

        //(*** F123 integration and asymptote ***)
        if ( epsi == 1 )   //(* sphere *)
        {
            if ( (xrad > lim3) && (xrad < limb) )  //(* numerical integration for d=2 *)
            {
                if ( dim == 1 ) qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,14,1,F123);
                if ( dim == 2 ) qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,11,2,F123);
            }

            if ( xrad > limb )
            {
                arg31 = (z+2.*v+1)*atan(4.*x2z);
                nen31 = pow(1.+16.*x2z*x2z,(z+2*v+1)/2.);
                arg32 = (z+2.*v)*atan(4.*x2z);
                nen32 = pow(1.+16.*x2z*x2z,(z+2.*v)/2.);
                arg33 = (z+2.*v-1)*atan(4.*x2z);
                nen33 = pow(1.+16.*x2z*x2z,(z+2*v-1)/2.);

                F32as1z = e0*e0*pz2v*(1.+cos(M_PI*v)*cos(arg31)/nen31-sin(M_PI*v)*sin(arg31)/nen31);
                F32as2z = 2.*e0*e1*(1./(2.*x2z))*pz2v1*(cos(M_PI*(2.*v-1)/2.)*cos(arg32)/nen32-sin(M_PI*(2.*v-1)/2.)*sin(arg32)/nen32);
                F32as3z = e1*e1*(1./(4.*x2z*x2z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg33)/nen33-sin(M_PI*(v-1))*sin(arg33)/nen33);
                F32asz = preg1*preg1*pow(x2z,2*v)*(1./2.)*(F32as1z+F32as2z+F32as3z);
                F123 = c3*F32asz;
            }
        }
        else  //(* for ellpisoids *)
        {
            if ( xrad > lim3 )
                qrombpi2(epsi,rn,params.p1,rho2,1,1,params.sigma,q,1,1,11,13,dim,F123);
        }

        return (F121+F122+F123)/vv;
    }               // (* of ex = false = asymptotic expansion *)

    //(*** ex = true, analytical continuation for d=3 ***)
    else          //(* for d=3 *)
    {
        //(*** F121 integration and continuation ***)
        if ( epsi == 1 )
        {
            if ( (xradp > lim1) && (xradp < limb) )  //(* numerical integration for d=3 *)
                qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,6,3,F121);

            if ( xradp >= limb )  //(* analytical continuation for d=3 *)
            {
                cc[3] = sqr(1.-rho2)*sqr(params.p1);
                argc = (z-4.+1)*atan(2.*x1zz);
                nenc = pow(1.+sqr(2.*x1zz),(z-4+1)/2.);
                ac[3] = pz4*pow(x2zz,-4)*(1./2.)*(1.+cos(argc)/nenc);
                cc[8] = -sqr(1.-rho2)*2.*params.p1;
                argh = (z-5.+1)*atan(2.*x1zz);
                nenh = pow(1.+sqr(2.*x1z),(z-5+1)/2.);
                ac[8] = pz5*pow(x2zz,-5)*(1./2.)*sin(argh)/nenh;
                cc[10] = sqr(1.-rho2);
                argj = (z-6+1.)*atan(2.*x1zz);
                nenj = pow(1.+sqr(2.*x1zz),(z-6+1)/2.);
                ac[10] = pz6*pow(x2zz,-6)*(1./2.)*(1-cos(argj)/nenj);
                F121 = 9.*(cc[3]*ac[3]+cc[8]*ac[8]+cc[10]*ac[10]);
            }
        }
        else  //(* for ellpisoids *)
        {
            if ( xradp > lim1 )
                qrombpi2(epsi,rn,params.p1,rho2,1,1,params.sigma,q,1,1,11,14,dim,F121);
        }

        //(*** F122 integration and continuation ***)
        if ( epsi == 1 )
        {
            if ( (xrad > lim2) && (xrad < limb) )   //(* numerical integration for d=3 *)
                qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,7,3,F122);

            if ( xrad >= limb )    //(* analytical continuation for d=3 *)
            {
                cc[2] = 2.*params.p1*rho2*(1-rho2);
                argbm = (z-4.+1)*atan(x1zz-x2zz);
                nenbm = pow(1.+sqr(x1zz-x2zz),(z-4+1)/2.);
                argbp = (z-4.+1)*atan(x1zz+x2zz);
                nenbp = pow(1.+sqr(x1zz+x2zz),(z-4+1)/2.);
                ac[2] = pz4*pow(x2zz,-4)*(1./2.)*(cos(argbm)/nenbm+cos(argbp)/nenbp);
                cc[5] = -2.*params.p1*rho2*(1.-rho2);
                argep = (z-5.+1)*atan(x1zz+x2zz);
                nenep = pow(1.+sqr(x1zz+x2zz),(z-5+1)/2.);
                argem = (z-5.+1)*atan(x1zz-x2zz);
                nenem = pow(1.+sqr(x1zz-x2zz),(z-5+1)/2.);
                ac[5] = pz5*pow(x2zz,-5)*(1./2.)*(sin(argep)/nenep-sin(argem)/nenem);
                cc[7] = -2.*rho2*(1.-rho2);
                arggp = (z-5.+1)*atan(x2zz+x1zz);
                nengp = pow(1.+sqr(x2zz+x1zz),(z-5+1)/2.);
                arggm = (z-5.+1)*atan(x2zz-x1zz);
                nengm = pow(1.+sqr(x2zz-x1zz),(z-5+1)/2.);
                ac[7] = pz5*pow(x2zz,-5)*(1./2.)*(sin(arggp)/nengp-sin(arggm)/nengm);
                cc[9] = 2.*rho2*(1.-rho2);
                //argim = (z-6.+1)*atan(x1zz-x2zz);
                //nenim = pow(1.+sqr(x1zz-x2zz),(z-6+1)/2.);
                //argip = (z-6.+1)*atan(x1zz+x2zz);
                //nenip = pow(1.+sqr(x1zz+x2z),(z-6+ 1)/2.);
                ac[9] = pz6*pow(x2zz,-6)*(1./2.)*(cos(argbm)/nenbm-cos(argbp)/nenbp);
                F122 = 9.*(cc[2]*ac[2]+cc[5]*ac[5]+cc[7]*ac[7]+cc[9]*ac[9]);
            }
        }
        else  //(* for ellpisoids *)
        {
            if ( xrad > lim2 )
                qrombpi2(epsi,rn,params.p1,rho2,1,1,params.sigma,q,1,1,11,15,dim,F122);
        }

        //(*** F123 integration and analytical continuation ***)
        if ( epsi == 1 )  //(* sphere *)
        {
            if ( (xrad > lim3) && (xrad < limb) )     //(* numercial integration for d=3 *)
                qrombpq(1,rn,params.p1,rho2,1,0.1,q,10,0,0,8,3,F123);

            if ( xrad >= lim3 )      //(* analytical continuation for d=3 *)
            {
                cc[1] = sqr(rho2);
                arga = (z-4.+1)*atan(2.*x2zz);
                nena = pow(1.+sqr(2*x2zz),(z-4+1)/2.);
                ac[1] = pz4*pow(x2zz,-4)*(1./2.)*(1+cos(arga)/nena);
                cc[4] = -2.*sqr(rho2);
                argd = (z-5.+1)*atan(2.*x2zz);
                nend = pow(1.+sqr(2*x2zz),(z-5+1)/2.);
                ac[4] = pz5*pow(x2zz,-5)*(1./2.)*sin(argd)/nend;
                cc[6] = sqr(rho2);
                argf = (z-6.+1)*atan(2.*x2zz);
                nenf = pow(1.+sqr(2*x2zz),(z-6+1)/2.);
                ac[6] = pz6*pow(x2zz,-6)*(1./2.)*(1-cos(argf)/nenf);
                F123 = 9.*(cc[1]*ac[1]+cc[4]*ac[4]+cc[6]*ac[6]);
            }
        }
        else  //(* for ellpisoids *)
        {
            if ( xrad > lim3 )
                qrombpi2(epsi,rn,params.p1,rho2,1,1,params.sigma,q,1,1,11,16,dim,F123);
        }
        return (F121+F122+F123)/vv3;
    }             //(* of d=3 *)
}



/* ***************** for hyperbolic r^(-a) profiles, 3D **************** */
double CLASSLIB::pqcoreshellf( double /*rho1*/, double rho2, double p1, double /*p2*/,
                               double /*alf1*/, double /*alf2*/, double rn, double d,
                               double sigma, double q ) const
{
      //const double eps  = 0.0000000001;
      const int    nmax = 300;
      const double delc = 0.0001;

      //double cpq[7], fpq[7]; // : array[1..6] of extended;;
      //double nenner,zaehler;
      double c1/*,c2*/,c3/*,c4,f1,f2,f3,fq,x,y*/,z/*,yz,da*/,zz;
      int    /*i,j,*/n/*,m*/;
      double /*zl,x1,xz,*/x1z,x12z,b1s,v,e0,e1/*,gb1s*/,pz2v,pz2v1/*,pz2v2*/,preg1,lim,lim1,lim2,lim3,vv,del,dim;
      double xrad,xradp,x2z,x22z;
      double F12as1z,F12as2z/*,F12as3z*/,F12asz,F32as1z,F32as2z/*,F32as3z*/,F32asz,FF1;
      double arg11,nen11,arg12,nen12/*,arg13,nen13*/,arg31,nen31,arg32,nen32/*,arg33,nen33*/;
      double F12sez1,oldF12sez1,F121,F12sez2,oldF12sez2/*,F122,F122a,F122b,F12sez3,oldF12sez3*/,F123;
      //double F121as1z,F121as2z,F121as3z,F121as4z,F121asz,F122as1z,F122as2z,F122as3z,F122as4z,F122asz;
      //double F123as1z,F123as2z,F123as3z,F123as4z,F123asz,F124as1z,F124as2z,F124as3z,F124as4z,F124asz;
      //double arg0,nen0,arg,nen,arg1,nen1,arg2,nen2,xijm,arglmz,nenlmz,xijp,arglpz,nenlpz;
      //double z12v[nmax+1],b1sv[nmax+1],fkv[nmax+1],sum121[nmax+1],sum122[nmax+1],sum123[nmax+1]; // : Array[0..nmax] of extended;;

      // p1:=fabs(p1);        /* ignore negative radii */

      /* convert from input n-average to calculation w-average */
      /* rn:=(1+6*sigma*sigma)*rn; */


      /* ################################### */
      dim = d;
      z = (1-sqr(sigma))/sqr(sigma);
      zz = z;
      z = z+2*dim;
      xrad = q*rn;
      xradp = q*p1*rn;
      x1z = q*p1*rn/(2*(zz+1));
      x12z = x1z*x1z;
      x2z = q*rn/(2*(zz+1));
      x22z = x2z*x2z;

      lim = 18*exp(-5*sigma);

      if ( dim==1 )
      {     /* disk */
          b1s = 3/2.0;
          v = -1;
          e0 = 1;
          e1 = 0;
          //gb1s = sqrt(M_PI)/2.0;
          preg1 = 1/2.0;
          pz2v = 1/z;
          pz2v1 = pz2v/(z-1);
          //pz2v2 = pz2v1/(z-2);
          lim1 = lim*1.2;           /* for term 1 */
          lim2 = lim*1.4;           /* for term 2 */
          lim3 = lim*1.0;           /* for term 3 */
      }
      if ( dim==2 )
      {     /* cylinder */
          b1s = 2;
          v = -3/2.0;
          e0 = 1;
          e1 = -9/16.0;
          //gb1s = 1;
          preg1 = 1/sqrt(M_PI);
          pz2v = gamma(z-1/2.0)/gamma(z+1);
          pz2v1 = pz2v/(z-3/2.0);
          //pz2v2 = pz2v1/(z-5/2.0);
          lim1 = lim*1.2;           /* for term 1 */
          lim2 = lim*1.4;           /* for term 2 */
          lim3 = lim*1.0;           /* for term 3 */
      }
      if ( dim==3 )
      {     /* sphere */
          b1s = 5/2.0;
          v = -2;
          e0 = 1;
          e1 = -1;
          //gb1s = 3*sqrt(M_PI)/4.0;
          preg1 = 3/4.0;
          pz2v = 1/(z*(z-1));
          pz2v1 = pz2v/(z-2);
          //pz2v2 = pz2v1/(z-3);
          lim1 = lim*1.2;           /* for term 1 */
          lim2 = lim*1.4;           /* for term 2 */
          lim3 = lim*1.0;           /* for term 3 */
      }

      c1 = (1-rho2)*pow(p1,dim);
      c3 = rho2;
      vv = c1+c3;

      if ( xradp<lim1 )
      {
          double z12v = 1; // Kein Array nötig
          double b1sv = 1;
          double fkv = 1;
          F12sez1 = 1.0;
          oldF12sez1 = 1.0;
          for ( n=1; n<=nmax; n++ )
          {
              z12v = z12v*((z+1)-2+2*n)*((z+1)-1+2*n);
              b1sv = b1sv*(b1s-1+n);
              fkv  = fkv*n;
              F12sez1 = F12sez1+pow(-x12z,n)*z12v/(b1sv*fkv);
              del = fabs((F12sez1-oldF12sez1)/F12sez1);
              if ( del<delc ) break;
              oldF12sez1 = F12sez1;
          }
          F121 = F12sez1;
      }

      if ( xrad<lim2 )
      {
          double z12v = 1;
          double b1sv = 1;
          double fkv = 1;
          F12sez2 = 1.0;
          oldF12sez2 = 1.0;
          for ( n=1; n<=nmax; n++ )
          {
              z12v = z12v*((z+1)-2+2*n)*((z+1)-1+2*n);
              b1sv = b1sv*(b1s-1+n);
              fkv = fkv*n;
              F12sez2 = F12sez2+pow(-x22z,n)*z12v/(b1sv*fkv);
              del = fabs((F12sez2-oldF12sez2)/F12sez2);
              if ( del<delc ) break;
              oldF12sez2 = F12sez2;
          }
          F123 = F12sez2;
      }

      if ( xradp>=lim1 )
      {
          arg11 = (z+v+1)*atan(2*x1z);
          nen11 = pow(1+4*x1z*x1z,(z+v+1)/2.0);
          arg12 = (z+v)*atan(2*x1z);
          nen12 = pow(1+4*x1z*x1z,(z+v)/2.0);
          //arg13 = (z+v-1)*atan(2*x1z);
          //nen13 = pow(1+4*x1z*x1z,(z+v-1)/2.0);

          F12as1z = e0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);
          F12as2z = e1*(1/(2*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);
          F12asz = preg1*pow(x1z,v)*(F12as1z+F12as2z);
          F121 = F12asz;
      }

      if ( xrad>=lim3 )
      {
          arg31 = (z+v+1)*atan(2*x2z);
          nen31 = pow(1+4*x2z*x2z,(z+v+1)/2.0);
          arg32 = (z+v)*atan(2*x2z);
          nen32 = pow(1+4*x2z*x2z,(z+v)/2.0);
          //arg33 = (z+v-1)*atan(2*x2z);
          //nen33 = pow(1+4*x2z*x2z,(z+v-1)/2.0);

          F32as1z = e0*pz2v*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);
          F32as2z = e1*(1/(2*x2z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);
          F32asz = preg1*pow(x2z,v)*(F32as1z+F32as2z);
          F123 = F32asz;
      }

      FF1 = (c1*F121+c3*F123)/vv;
      return FF1*FF1;
}
#endif // USE_pqcoreshell



#ifdef USE_f2dcoreshell
/**
 * @brief CLASSLIB::f2dcoreshell - SZ-averaged <F(q)>^2 for core/shell-systems
 * @param d     = 3
 * @param rho   = CLASSLIB::rho
 * @param p     = CLASSLIB::p1
 * @param alf   = local var
 * @param r     = CLASSLIB::radiusi
 * @param sigma = CLASSLIB::sigma
 * @param q     = local var
 * @return
 *
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::f2dcoreshell( double alf, double q ) const     //{NV}
{
    const double d = 3;

    double a1, a2, a3, b1, b2, b3;
    double nenner, hgeo;

    nenner = (1/d)-(params.rho/(d-alf))+(params.rho*exp((alf-d)*log(params.p1))/(d-alf));
    a1 = (1/d);
    a2 = -(params.rho/(d-alf));
    a3 = (params.rho*exp((alf-d)*log(params.p1))/(d-alf));
    b1 = a1*f2dhyper(d,0.001,params.p1*params.radiusi,params.sigma,q);
    b2 = a2*f2dhyper(d,alf,params.p1*params.radiusi,params.sigma,q);
    b3 = a3*f2dhyper(d,alf,params.radiusi,params.sigma,q);
    hgeo = (b1+b2+b3)/nenner;
    return hgeo*hgeo;
}
#endif // USE_f2dcoreshell



#ifdef USE_polyvesicle
//(* *************************** vesicles **************************** *)
/**
 * @brief CLASSLIB::polyvesicle
 * @param ro     = CLASSLIB::length
 * @param ri     = CLASSLIB::radius
 * @param sigmar = CLASSLIB::sigma
 * @param sigmal = CLASSLIB::sigmal
 * @param q      = local var from calling routine
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::polyvesicle(double q) const        //{NV}
{
    double r, d, yr, yd, zr, zd, c,
           gr1, gr2, gr3, gd1, gd2, gd3,
           hr1, hr2, hr3, hd1, hd2, hd3;

    r = (params.length+params.radius)/2;
    d = (params.length-params.radius)/2;
    zr = (1-sqr(params.sigma))/sqr(params.sigma);
    zd = (1-sqr(params.sigmal))/sqr(params.sigmal);
    yr = q*r/(zr+1);
    yd = q*d/(zd+1);
    c = 1+d*d/(12*r*r);

    if ( q*r < 0.001 )
        return 1;
    gr1 = 1./(2*zr*(zr-1)*yr*yr);
    gr2 = gr1/((zr-2)*yr);
    gr3 = gr2/((zr-3)*yr);
    gd1 = 1./(2*zd*(zd-1)*yd*yd);
    gd2 = 1./(2*zd*yd);
    gd3 = 1./2.;

    hr1 = gr1*(1-cos((zr-2+1)*atan(2*yr))/exp((zr-2+1)*log(1+4*yr*yr)/2.));    //(* <sin^2(yr)yr^-2> *)
    hr2 = gr2*(sin((zr-3+1)*atan(2*yr))/exp((zr-3+1)*log(1+4*yr*yr)/2.));      //(* <sin(yr)cos(yr)yr^-3> *)
    hr3 = gr3*(1+cos((zr-4+1)*atan(2*yr))/exp((zr-4+1)*log(1+4*yr*yr)/2.));    //(* <cos^2(yr)yr^-4> *)
    hd1 = gd1*(1-cos((zd-2+1)*atan(2*yd))/exp((zd-2+1)*log(1+4*yd*yd)/2.));    //(* <sin^2(yd)yd^-2> *)
    hd2 = gd2*(sin((zd-1+1)*atan(2*yd))/exp((zd-1+1)*log(1+4*yd*yd)/2.));      //(* <sin(yd)cos(yd)yd^-1> *)
    hd3 = gd3*(1+cos((zd+0+1)*atan(2*yd))/exp((zd+0+1)*log(1+4*yd*yd)/2.));    //(* <cos^2(yd)yd^0> *)

    return (hr1*hd1-2*hr2*hd2+2*hr2*hd1+hr3*hd3-2*hr3*hd2+hr3*hd1)/(c*c);
}
#endif // USE_polyvesicle



#ifdef USE_f2dpolyvesicle
//(* *************************** vesicles **************************** *)
/**
 * @brief CLASSLIB::f2dpolyvesicle
 * @param ro     = CLASSLIB::length
 * @param ri     = CLASSLIB::radius
 * @param sigmar = CLASSLIB::sigma
 * @param sigmad = CLASSLIB::sigmal
 * @param q      = local var from calling routine
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::f2dpolyvesicle(double q) const         //{NV}
{
    double r, d, yr, yd, zr, zd, c, pq,
           gr1, gr2, gd1, gd2,
           hr1, hr2, hr3, hd1, hd2, hd3;

    r = (params.length+params.radius)/2.;
    d = (params.length-params.radius)/2.;
    zr = (1-sqr(params.sigma))/sqr(params.sigma);
    zd = (1-sqr(params.sigmal))/sqr(params.sigmal);
    yr = q*r/(zr+1.);
    yd = q*d/(zd+1.);
    c = 1+d*d/(12.*r*r);

    if ( q*r < 0.001 )
        return 1;
    if ( params.sigma < 0.05 )
    {
        gr1 = gammaratio(-1+1,1,zr)*exp(-1*log(yr));
        gr2 = gammaratio(-2+1,1,zr)*exp(-2*log(yr));
    }
    else
    {
        gr1 = gamma(zr-1+1)*exp(-1*log(yr))/gamma(zr+1);            // (* <sin(yr)yr^-1> *)
        gr2 = gamma(zr-2+1)*exp(-2*log(yr))/gamma(zr+1);            // (* <cos(yr)yr^-2> *)
    }
    if ( params.sigmal < 0.05 )
    {
        gd1 = gammaratio(-1+1,1,zd)*exp(-1*log(yd));
        gd2 = gammaratio(-0+1,1,zd);
    }
    else
    {
        gd1 = gamma(zd-1+1)*exp(-1*log(yd))/gamma(zd+1);            // (* <sin(yd)yd^-1> *)
        gd2 = gamma(zd-0+1)/gamma(zd+1);                            // (* <cos(yd)yd^0> *)
    }

    hr1 = gr1*sin((zr-1+1)*atan(yr))/exp((zr-1+1)*log(1+yr*yr)/2.);  // (* <sin(yr)yr^-1> *)
    hr2 = gr2*cos((zr-2+1)*atan(yr))/exp((zr-2+1)*log(1+yr*yr)/2.);  // (* <cos(yr)yr^-2> *)
    hr3 = hr2;                                                       // (* <cos(yr)yr^-2> *)
    hd1 = gd1*sin((zd-1+1)*atan(yd))/exp((zd-1+1)*log(1+yd*yd)/2.);  // (* <sin(yd)yd^-1> *)
    hd2 = gd2*cos((zd-0+1)*atan(yd))/exp((zd-0+1)*log(1+yd*yd)/2.);  // (* <cos(yd)yd^0> *)
    hd3 = hd1;                                                       // (* <sin(yd)yd^-1> *)

    pq = (hr1*hd1-hr2*hd2+hr3*hd3)/c;

    return pq*pq;
}
#endif



#ifdef USE_polyliposome
//(*** ********* polydisperse liposomes ********* *)
/**
 * @brief CLASSLIB::polyliposome
 * @param lt     = CLASSLIB::length
 * @param rc     = CLASSLIB::radius
 * @param lh     = CLASSLIB::sigma
 * @param lin    = CLASSLIB::sigmal
 * @param lex    = CLASSLIB::shellno
 * @param nl     = CLASSLIB::alphash
 * @param sigmar = CLASSLIB::ceff
 * @param sigmal = CLASSLIB::reff
 * @param cc     = CLASSLIB::latt_a
 * @param ch     = CLASSLIB::latt_b
 * @param ct     = CLASSLIB::latt_c
 * @param cin    = CLASSLIB::domainsize
 * @param cex    = CLASSLIB::aziwidth
 * @param dim    = NV NEU
 * @param q      = local var from calling routine
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::polyliposome(int dim, double q) const      //{NV} NEU
{
    const double lt     = params.length; // TODO
    const double rc     = params.radius;
    const double lh     = params.sigma;
    const double lin    = params.sigmal;
    const double lex    = params.shellno;
    const double nl     = params.alphash;
    const double sigmar = params.ceff;
    const double sigmal = params.reff;
    const double cc     = params.uca;
    const double ch     = params.ucb;
    const double ct     = params.ucc;
    const double cin    = params.domainsize;
    const double cex    = params.aziwidth;

    const int nmax = 300;
    const double delc = 0.0001;

    int i,j,ii,jj,n,m,ncell,cnt,inmax;
    double philiph,philipt,phiax,phiin,phiout,lliph,llipt,lout,rad,radv,len,vv;
    double z,zl,lim,lim1,b1s,v,e0,e1,/*gb1s,*/preg1,pz2v,pz2v1,pz2v2,xrad,x1z,x12z,xrz;
    double lq[nmax+1],phi[nmax+1],ll[nmax+1],qq[nmax+1],rr[nmax+1],pp[nmax+1],z12v[nmax+1],b1sv[nmax+1],fkv[nmax+1],sum12[nmax+1]; // Array[0..nmax] of extended;
    double F12sum,a1,a2,F12sez,oldF12sez,del,/*F12ser,F12,*/rmax,xmax;
    double /*arg0,nen0,*/arg,nen,arg1,nen1,arg2,nen2,F12asz,xijm,arglmz,nenlmz,xijp,arglpz,nenlpz;
    double F12as1z,F12as2z,F12as3z,F12as4z/*,F12asy*/;

    philiph = ch;      // (* head group *)
    philipt = ct;      // (* tail *)
    phiax = cc;        // (* axon *)
    phiin = cin;       // (* intra cell *)
    phiout = cex;      // (* extra cell *)

    rad = rc;          // (* vesicle inner radius *)
    lliph = lh;        // (* head group *)
    llipt = lt;        // (* tail *)
    //?lin = lin;         // (* intra cell *)
    lout = lex;        // (* extra cell *)

    len = lout+2*(2*lliph+llipt)+lin;
    ncell = static_cast<int>(round(nl));
    rmax = rad+ncell*len;

    lq[1] = lout;
    lq[2] = lq[1]+lliph;
    lq[3] = lq[2]+llipt;
    lq[4] = lq[3]+lliph;
    lq[5] = lq[4]+lin;
    lq[6] = lq[5]+lliph;
    lq[7] = lq[6]+llipt;
    lq[8] = lq[7]+lliph;

    for ( i=1; i<=8; i++ ) qq[i] = lq[i]/len;

    phi[1] = phiax;       //(* vesicle interior *)
    ll[1]  = 0.0;
    rr[1]  = rad+ll[1];
    pp[1]  = 1.0;
    phi[2] = philiph;     //(* vesicle bilayer: head group *)
    ll[2]  = lliph;
    rr[2]  = rad+ll[2];
    pp[2]  = rr[2]/rad;
    phi[3] = philipt;     //(* vesicle bilayer: tail group *)
    ll[3]  = ll[2]+llipt;
    rr[3]  = rad+ll[3];
    pp[3]  = rr[3]/rad;
    phi[4] = philiph;     //(* vesicle bilayer: head group *)
    ll[4]  = ll[3]+lliph;
    rr[4]  = rad+ll[4];
    pp[4]  = rr[4]/rad;

    radv = rr[4];        //(* vesicle radius + bilayer *)
    cnt = 4;
    for ( i=1; i<=ncell; i++ )
    {
        phi[cnt+1] = phiout;          //(* extra cell *)
        phi[cnt+2] = philiph;         //(* head group *)
        phi[cnt+3] = philipt;         //(* tail group *)
        phi[cnt+4] = philiph;         //(* head group *)
        phi[cnt+5] = phiin;           //(* intra cell *)
        phi[cnt+6] = philiph;         //(* head group *)
        phi[cnt+7] = philipt;         //(* tail group *)
        phi[cnt+8] = philiph;         //(* head group *)
        ll[cnt+1] = (i-1+qq[1])*len;
        ll[cnt+2] = (i-1+qq[2])*len;
        ll[cnt+3] = (i-1+qq[3])*len;
        ll[cnt+4] = (i-1+qq[4])*len;
        ll[cnt+5] = (i-1+qq[5])*len;
        ll[cnt+6] = (i-1+qq[6])*len;
        ll[cnt+7] = (i-1+qq[7])*len;
        ll[cnt+8] = (i-1+qq[8])*len;
        rr[cnt+1] = radv+ll[cnt+1];
        rr[cnt+2] = radv+ll[cnt+2];
        rr[cnt+3] = radv+ll[cnt+3];
        rr[cnt+4] = radv+ll[cnt+4];
        rr[cnt+5] = radv+ll[cnt+5];
        rr[cnt+6] = radv+ll[cnt+6];
        rr[cnt+7] = radv+ll[cnt+7];
        rr[cnt+8] = radv+ll[cnt+8];
        pp[cnt+1] = rr[cnt+1]/rad;
        pp[cnt+2] = rr[cnt+2]/rad;
        pp[cnt+3] = rr[cnt+3]/rad;
        pp[cnt+4] = rr[cnt+4]/rad;
        pp[cnt+5] = rr[cnt+5]/rad;
        pp[cnt+6] = rr[cnt+6]/rad;
        pp[cnt+7] = rr[cnt+7]/rad;
        pp[cnt+8] = rr[cnt+8]/rad;
        cnt += 8;
    }
    inmax = cnt;
    phi[cnt+1] = 0.0;

    vv = 0;
    for ( i=1; i<=inmax; i++ )
    {
        for ( j=1; j<=inmax; i++ )
            vv += (phi[i]-phi[i+1])*(phi[j]-phi[j+1])*pow(rr[i],dim)*pow(rr[j],dim);
    }

    z = (1-sigmar*sigmar)/(sigmar*sigmar);
    zl = (1-sigmal*sigmal)/(sigmal*sigmal);

    lim = 18 * exp(-5*sigmar);

    switch ( dim )
    {
    case 1:     //(* disk *)
        b1s = 3./2.;
        v = -1;
        e0 = 1;
        e1 = 0;
        //gb1s = sqrt(M_PI)/2.;
        preg1 = 1./2.;
        pz2v = 1./(z*(z-1));
        pz2v1 = pz2v/(z-2);
        pz2v2 = pz2v1/(z-3);
        lim1 = lim*1.2;           //(* for term 1 *)
        break;
    case 2:     //(* cylinder *)
        b1s = 2;
        v = -3./2.;
        e0 = 1;
        e1 = -9./16.;
        //gb1s = 1;
        preg1 = 1./sqrt(M_PI);
        pz2v = 1./(z*(z-1)*(z-2));
        pz2v1 = pz2v/(z-3);
        pz2v2 = pz2v1/(z-4);
        lim1 = lim*1.2;           //(* for term 1 *)
        break;
    case 3:     //(* sphere *)
        b1s = 5./2.;
        v = -2;
        e0 = 1;
        e1 = -1;
        //gb1s = 3.*sqrt(M_PI)/4.;
        preg1 = 3./4.;
        pz2v = 1./(z*(z-1)*(z-2)*(z-3));
        pz2v1 = pz2v/(z-4);
        pz2v2 = pz2v1/(z-5);
        lim1 = lim*1.4;           //(* for term 1 *)
        break;
    }

    xrad = q*rad;
    xrz = xrad/(z+1);
    x1z = q*rad/(2*(z+1));
    x12z = x1z*x1z;
    xmax = q*rmax;

    if ( xmax < lim1 )
    {
        F12sum = 0.0;
        for ( ii=1; ii<=inmax; ii++ )
        {
            a1 = (phi[ii]-phi[ii+1])*pow(rr[ii],dim);
            for ( jj=1; jj<=inmax; jj++ )
            {
                a2 = (phi[jj]-phi[jj+1])*pow(rr[jj],dim);
                z12v[0] = 1;
                b1sv[0] = 1;
                fkv[0] = 1;
                F12sez = 1.0;
                oldF12sez = 1.0;
                for ( n=1; n<=nmax; n++ )
                {
                    z12v[n] = z12v[n-1]*((z+1)-2+2*n)*((z+1)-1+2*n);
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);
                    fkv[n] = fkv[n-1]*n;
                    sum12[n] = 0;
                    for ( m=0; m<=n; m++ )
                        sum12[n] += pow(pp[ii],2*m)*pow(pp[jj],2*(n-m))/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);
                    F12sez += pow(-x12z,n)*z12v[n]*sum12[n];
                    del = fabs((F12sez-oldF12sez)/F12sez);
                    if ( del < delc ) break;
                    oldF12sez = F12sez;
                } // for n
                F12sum += a1*a2*F12sez;
            } // for jj
        } // for ii
        return F12sum/vv;
    }

    //(*** term #1 ***)
    else //if ( xmax >= lim1 )
    {
        //arg0 = 2*q*rad;
        //nen0 = 1.0;
        arg = (z+2*v+1)*atan(2*xrz);
        nen = pow(1+4*xrz*xrz,(z+2*v+1)/2.);
        arg1 = (z+2*v)*atan(2*xrz);
        nen1 = pow(1+4*xrz*xrz,(z+2*v)/2.);
        arg2 = (z+2*v-1)*atan(2*xrz);
        nen2 = pow(1+4*xrz*xrz,(z+2*v-1)/2.);

        F12asz = 0.0;
        for ( ii=1; ii<=inmax; ii++ )
        {
            a1 = (phi[ii]-phi[ii+1])*pow(rr[ii],dim)*pow(pp[ii],v);
            for ( jj=1; jj<=inmax; jj++ )
            {
                a2 = (phi[jj]-phi[jj+1])*pow(rr[jj],dim)*pow(pp[jj],v);
                xijm = (ll[ii]-ll[jj])*q/(zl+1);
                arglmz = (zl+1)*atan(xijm);
                nenlmz = pow(1+xijm*xijm,(zl+1)/2.);
                xijp = (ll[ii]+ll[jj])*q/(zl+1);
                arglpz = (zl+1)*atan(xijp);
                nenlpz = pow(1+xijp*xijp,(zl+1)/2.);
                F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(M_PI*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(M_PI*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));
                F12as2z = e0*e1*(1./(pp[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));
                F12as3z = e1*e0*(1./(pp[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));
                F12as4z = e1*e1*(1./(pp[ii]*pp[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(M_PI*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(M_PI*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));
                F12asz = F12asz+a1*a2*(F12as1z+F12as2z+F12as3z+F12as4z);
            } // for jj
        } // for ii
        return preg1*preg1*pow(xrz/2,2*v)*(1./2.)*F12asz/vv;
    }
}
#endif // USE_polyliposome



#ifdef USE_polycube
//(*** ********* polydisperse cubes (isotropic average) ********* *)
/**
 * @brief SasCalculation::polycube
 * @param a     = radius
 * @param sigma = sigma
 * @param pf    = Flag 0/1
 * @param q     = local var
 * @return
 *             pq=CLASSLIB::polycube(CLASSLIB::radius,CLASSLIB::sigma,0,q);
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::polycube(bool pf, double q) const      //{NV}
{
    const double a = params.radius;

    const int maxit = 7; // 10;

    double intrc, oldintrc;
    int itrc, trapzditrc;

    intrc = 0;
    //(*** integral over alpha = 0 .. pi ***)
    //(*** integrate over alpha = 0 .. pi/2 due to 4-fold symmetry ***)
    trapcube(0+eps3,M_PI/2-eps3,a,pf,q,intrc,1,trapzditrc);
    oldintrc = intrc;
    for ( itrc=2; itrc<=maxit; itrc++ )
    {
        trapcube(0+eps3,M_PI/2-eps3,a,pf,q,intrc,itrc,trapzditrc);
        if ( oldintrc != 0.0 && fabs(1-intrc/oldintrc) < eps3 ) break;
        oldintrc = intrc;
    }
    //(* normalization: division by  integral 0..pi/2  sin(theta) d(theta) = 1 *)
    return intrc;
}
#endif // USE_polycube



#ifdef USE_burger
//(********************** Burger Peak-Shape ***************************)
/**
 * @brief SasCalculation::burger
 * @param del
 * @param v
 * @param x2
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::burger(double del, double v, double x2) const      //{NV}
{
    int i, imax;
    double gv, w;

    gv = sqrt(M_PI)*gamma((v+1)/2.)/gamma(v/2.);
    w = 2./(M_PI*del);
    imax = round(10*v+50);
    for ( i=0; i<=imax; i++ )
        w = w/(1+gv*gv*4.*x2/(M_PI*M_PI*(v/2.+i)*(v/2.+i)));
    return w;
}
#endif // USE_burger



#ifdef USE_angleuv
//(**************** enclosed angle between vectors u and v *******************)
/**
 * @brief SasCalculation::angleuv
 * @param ux
 * @param uy
 * @param uz
 * @param vx
 * @param vy
 * @param vz
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::angleuv(double ux, double uy, double uz, double vx,double vy, double vz) const     //{NV}
{
    double un, vn, uv;

    un = sqrt(ux*ux+uy*uy+uz*uz);
    vn = sqrt(vx*vx+vy*vy+vz*vz);
    if ( fabs(un) < eps9 || fabs(vn) < eps9 )
        return 0;
    uv = (ux*vx+uy*vy+uz*vz)/(un*vn);
    if ( fabs(uv-1) < eps9 ) // uv == +1
        return 0;
    if ( fabs(uv+1) < eps9 ) // uv == -1
        return M_PI;
    return M_PI/2.-atan(uv/sqrt(1.000000001-uv*uv));
}
#endif // USE_angleuv



#ifdef USE_lorentznorm3
//(* ************************** 3d-integral over lorentz(x)*sin(x) ********************************* *)
/**
 * @brief SasCalculation::lorentznorm3
 * @param a
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::lorentznorm3(double a) const           //{NV}
{
    double a0, a1, a2, a3, b0, b1, b2, b3;

    if ( a > 0.01 )
    {
        b0 = log(1+a*M_PI*M_PI)/2.;
        b1 = (a*M_PI*M_PI-log(1+a*M_PI*M_PI))/12.;
        b2 = (-2*a*M_PI*M_PI+a*a*M_PI*M_PI*M_PI*M_PI+2*log(1+a*M_PI*M_PI))/480.;
        b3 = (a*M_PI*M_PI*(6-3*a*M_PI*M_PI+2*a*a*M_PI*M_PI*M_PI*M_PI)-6*log(1+a*M_PI*M_PI))/60480.;
        return b0/a-b1/(a*a)+b2/(a*a*a)-b3/(a*a*a*a);
    }
    //if a<=0.01 then begin
    a0 = 2;
    a1 = (M_PI*M_PI-4);
    a2 = (18-12*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI)/2.;
    a3 = (-1440+360*M_PI*M_PI-30*M_PI*M_PI*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI*M_PI*M_PI)/6.;
    return a0-a1*a+a2*a*a-a3*a*a*a;
}
#endif // USE_lorentznorm3



#ifdef USE_gaussnorm3
//(* ************************** 3d-integral over gauss(x)*sin(x) ********************************* *)
/**
 * @brief SasCalculation::gaussnorm3
 * @param a
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::gaussnorm3(double a) const         //{NV}
{
    double a0,a1,a2,a3,b0,b1,b2,b3;

    if ( a > 0.01 )
    {
        b0 = (1-exp(-a*M_PI*M_PI))/2.;
        b1 = (1-exp(-a*M_PI*M_PI)*(1+a*M_PI*M_PI))/12.;
        b2 = (2-exp(-a*M_PI*M_PI)*(2+2*a*M_PI*M_PI+a*a*M_PI*M_PI*M_PI*M_PI))/240.;
        b3 = (6-exp(-a*M_PI*M_PI)*(6+6*a*M_PI*M_PI+3*a*a*M_PI*M_PI*M_PI*M_PI+a*a*a*M_PI*M_PI*M_PI*M_PI*M_PI*M_PI))/10080.;
        return b0/a-b1/(a*a)+b2/(a*a*a)-b3/(a*a*a*a);
    }
    //if a<=0.01 then begin
    a0 = 2;
    a1 = (M_PI*M_PI-4);
    a2 = (18-12*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI)/2.;
    a3 = (-1440+360*M_PI*M_PI-30*M_PI*M_PI*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI*M_PI*M_PI)/6.;
    return a0-a1*a+a2*a*a-a3*a*a*a;
}
#endif // USE_gaussnorm3



#ifdef USE_pearsonnorm3
//(* ************************** 3d-integral over pearson(x)*sin(x) ********************************* *)
/**
 * @brief SasCalculation::pearsonnorm3
 * @param a
 * @param b
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::pearsonnorm3(double a, double b) const         //{NV}
{
    const int    maxit = 10;

    double inttf,oldinttf,a0,a1,a2,a3,b0; //,b1,b2,b3;
    int    itf,trapzditx;

    b = b+0.000001;
    if ( a > 20 )
    {
        b0 = (1-exp((1-b)*log(1+a*M_PI*M_PI)))/(2*(b-1));
        //b1 = (-1+a*a*M_PI*M_PI-b*a*M_PI*M_PI*(1+a*M_PI*M_PI)+exp(b*log(1+a*M_PI*M_PI)))/(12*(b-2)*(b-1)*exp(b*log(1+a*M_PI*M_PI)));
        //b2 = (-2*a*M_PI*M_PI+a*a*M_PI*M_PI*M_PI*M_PI+2*log(1+a*M_PI*M_PI))/480;
        //b3 = (a*M_PI*M_PI*(6-3*a*M_PI*M_PI+2*a*a*M_PI*M_PI*M_PI*M_PI)-6*log(1+a*M_PI*M_PI))/60480;
        return b0/a;
    }
    if ( a <= 0.01 )
    {
        a0 = 2;
        a1 = (M_PI*M_PI-4)*b;
        a2 = (18-12*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI)*b*(b+1)/2.;
        a3 = (-1440+360*M_PI*M_PI-30*M_PI*M_PI*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI*M_PI*M_PI)*b*(b+1)*(b+2)/6.;
        return a0-a1*a+a2*a*a-a3*a*a*a;
    }
    inttf=0;
    pearsonintegral3(0,M_PI,a,b,inttf,1,trapzditx);
    oldinttf=inttf;
    for ( itf=2; itf<=maxit; itf++ )
    {
        pearsonintegral3(0,M_PI,a,b,inttf,itf,trapzditx);
        if ( fabs(1-inttf/oldinttf) < eps3 ) break;
        oldinttf = inttf;
    }
    return inttf;
}
#endif /// USE_pearsonnorm3



#ifdef USE_pearsonintegral3
//(* ********************* integration procedure for pearson(x)*sin(x) ***************************** *)
/**
 * @brief SasCalculation::pearsonintegral3
 * @param at
 * @param bt
 * @param a
 * @param b
 * @param sx
 * @param nx
 * @param trapzditx
 */
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::pearsonintegral3(double at, double bt, double a, double b,
                                      double &sx, int nx, int &trapzditx) const         //{NV}
{
    double xt,fa,fb,fx,tnmx,sumx,delx;

    if ( nx == 1 )
    {
        fa = sin(at)/exp(b*log(1+a*at*at));
        fb = sin(bt)/exp(b*log(1+a*bt*bt));
        sx = 0.5*(bt-at)*(fa+fb);
        trapzditx = 1;
    }
    else
    {
        tnmx = trapzditx;
        delx = (bt-at)/tnmx;
        xt = at+0.5*delx;
        sumx = 0.0;
        for ( int jx=1; jx<=trapzditx; jx++ )
        {
            fx = sin(xt)/exp(b*log(1+a*xt*xt));
            sumx += fx;
            xt = xt+delx;
        }
        sx = 0.5*(sx+(bt-at)*sumx/tnmx);
        trapzditx = 2*trapzditx;
    }
}
#endif // USE_pearsonintegral3



#ifdef USE_cosavm
//(************* function <cos(di+piy)cos(dj+pjy)y^(ei+ej)> **************)
/**
 * @brief SasCalculation::cosavm
 * @param d1 - (-pi/2) or 0
 * @param d2 - (-pi/2) or 0
 * @param e1
 * @param e2
 * @param p1
 * @param p2
 * @param y
 * @param z
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::cosavm(double d1, double d2, double e1, double e2, double p1, double p2, double y, double z) const //{NV}
{
    double nav,p2a,p2b,p2c,p2d,ava,avb,avc,avd,f2a,f2b,f2c,f2d,yz,pyz;

    if ( fabs(int(z)-z) < eps9 ) z=z+eps9;   //(* avoid division by zero for integer values of z *)

    nav=e1+e2;
    yz=y;
    pyz=pow(yz,nav);

    // Vergleiche, ob di oder dj == -1000 oder nicht sind. Ich stelle diese Gruppen etwas um,
    //  so dass die double-Vergleiche sicherer sind.
    if ( d1 > -1000 )
    {
        if ( d2 > -1000 )
        {   // if ( d1 != -1000 && d2 != -1000 )
            //(* (p1^e1)(p2^e2) <cos(d1+p1y) cos(d2+p2y) y^nav> *)
            if ( z > 1000 )
            {
                p2a=gammaratio(nav+1,1,z)*pyz/2.;
                p2b=gammaratio(nav+1,1,z)*pyz/2.;
                p2c=gammaratio(nav+2,1,z)*pyz/((z+nav+1)*2);
                p2d=gammaratio(nav+2,1,z)*pyz/((z+nav+1)*2);
            }
            else
            {
                p2a=gamma(z+nav+1)*pyz/(gamma(z+1)*2);
                p2b=gamma(z+nav+1)*pyz/(gamma(z+1)*2);
                p2c=gamma(z+nav+2)*pyz/(gamma(z+1)*(z+nav+1)*2);
                p2d=gamma(z+nav+2)*pyz/(gamma(z+1)*(z+nav+1)*2);
            }
            //ava:=cos((z+nav+1)*arctan((p1+p2)*yz))/exp((z+nav+1)*log(1+(p1+p2)*(p1+p2)*yz*yz)/2);
            ava=cos((z+nav+1)*atan((p1+p2)*yz))*pow(1+(p1+p2)*(p1+p2)*yz*yz,-(z+nav+1)/2.);
            //avb:=cos((z+nav+1)*arctan((p1-p2)*yz))/exp((z+nav+1)*log(1+(p1-p2)*(p1-p2)*yz*yz)/2);
            avb=cos((z+nav+1)*atan((p1-p2)*yz))*pow(1+(p1-p2)*(p1-p2)*yz*yz,-(z+nav+1)/2.);
            //avc:=sin((z+nav+1)*arctan((p1+p2)*yz))/exp((z+nav+1)*log(1+(p1+p2)*(p1+p2)*yz*yz)/2);
            avc=sin((z+nav+1)*atan((p1+p2)*yz))*pow(1+(p1+p2)*(p1+p2)*yz*yz,-(z+nav+1)/2.);
            //avd:=sin((z+nav+1)*arctan((p1-p2)*yz))/exp((z+nav+1)*log(1+(p1-p2)*(p1-p2)*yz*yz)/2);
            avd=sin((z+nav+1)*atan((p1-p2)*yz))*pow(1+(p1-p2)*(p1-p2)*yz*yz,-(z+nav+1)/2.);
            f2a=cos(d1+d2)*p2a*ava;
            f2b=cos(d1-d2)*p2b*avb;
            f2c=sin(d1+d2)*p2c*avc;
            f2d=sin(d1-d2)*p2d*avd;
            return pow(p1,e1)*pow(p2,e2)*(f2a+f2b-f2c-f2d);
        }
        else
        {   // if ( d1 != -1000 && d2 == -1000 )
            //(* (p1^e1)(p2^e2) <cos(d1+p1y) y^nav> *)
            if ( z > 1000 )
            {
                p2a=gammaratio(nav+1,1,z)*pyz;
                p2b=gammaratio(nav+2,1,z)*pyz/(z+nav+1);
            }
            else
            {
                p2a=gamma(z+nav+1)*pyz/(gamma(z+1));
                p2b=gamma(z+nav+2)*pyz/(gamma(z+1)*(z+nav+1));
            }
            //ava:=cos((z+nav+1)*arctan(p1*yz))/exp((z+nav+1)*log(1+p1*p1*yz*yz)/2);
            ava=cos((z+nav+1)*atan(p1*yz))*pow(1+p1*p1*yz*yz,-(z+nav+1)/2.);
            //avb:=sin((z+nav+1)*arctan(p1*yz))/exp((z+nav+1)*log(1+p1*p1*yz*yz)/2);
            avb=sin((z+nav+1)*atan(p1*yz))*pow(1+p1*p1*yz*yz,-(z+nav+1)/2.);
            return pow(p1,e1)*pow(p2,e2)*(cos(d1)*p2a*ava-sin(d1)*p2b*avb);
        }
    }
    else
    {
        if ( d2 > -1000 )
        {   // if ( d1 == -1000 && d2 != -1000 )
            //(* (p1^e1)(p2^e2) <cos(d2+p2y) y^nav> *)
            if ( z > 1000 )
            {
                p2a=gammaratio(nav+1,1,z)*pyz;
                p2b=gammaratio(nav+2,1,z)*pyz/(z+nav+1);
            }
            else
            {
                p2a=gamma(z+nav+1)*pyz/(gamma(z+1));
                p2b=gamma(z+nav+2)*pyz/(gamma(z+1)*(z+nav+1));
            }
            //ava:=cos((z+nav+1)*arctan(p2*yz))/exp((z+nav+1)*log(1+p2*p2*yz*yz)/2);
            ava=cos((z+nav+1)*atan(p2*yz))*pow(1+p2*p2*yz*yz,-(z+nav+1)/2.);
            //avb:=sin((z+nav+1)*arctan(p2*yz))/exp((z+nav+1)*log(1+p2*p2*yz*yz)/2);
            avb=sin((z+nav+1)*atan(p2*yz))*pow(1+p2*p2*yz*yz,-(z+nav+1)/2.);
            return pow(p1,e1)*pow(p2,e2)*(cos(d2)*p2a*ava-sin(d2)*p2b*avb);
        }
        else
        {   // if ( d1 == -1000 && d2 == -1000 )
            //(* yz:=y/(z+1); *)
            if ( z > 1000 ) p2a=gammaratio(nav+1,1,z)*pyz;
            else p2a=gamma(z+nav+1)*pyz/gamma(z+1);
            return pow(p1,e1)*pow(p2,e2)*p2a;
        }
    }


    return 0;
}
#endif // USE_cosavm



#ifdef USE_pgensphere1
//(* *********** P(q) fr harte Kugeln mit Burger-I-Asymptote *********** *)
//(* **** do not use for d<3 ! **** *)
/**
 * @brief SasCalculation::pgensphere1
 * @param a
 * @param r
 * @param sigma
 * @param q
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::pgensphere1(double a, double r, double q) const        //{NV}
{
    const int dim    = 3;            //(* dimension *)
    const int maxsum = 250;          //(* number of summations *)

    double pa1[maxsum+1], pb1[maxsum+1], pb2[maxsum+1], pz1[maxsum+1], fak[maxsum+1], pn[maxsum+1], xn[maxsum+1]; // array[0..maxsum] of extended;
    int n, m, sig;
    double a1,b1,b2,z,y,x,msum,p,
           zaehler,nenner,hgeo,oldhgeo,
           rt,sa1,sa2,sa3,sb1,sb2,sb3,sb4,sb5,sb6,sc1,sc2,f1,f2,f3,f4,f5;

    a1 = (dim-a)/2.;
    b1 = dim/2.;
    b2 = (dim+2-a)/2.;
    z = (1-sqr(params.sigma))/sqr(params.sigma);
    p = 1;
    y = q*r/(z+1);
    x = sqr(y/2.);
    pa1[0] = 1;
    pb1[0] = 1;
    pb2[0] = 1;
    pz1[0] = 1;
    fak[0] = 1;
    pn[0] = 1;
    xn[0] = 1;

    //(*** simplification ***)
    pa1[1] = a1;
    pb1[1] = b1;
    pb2[1] = b2;
    pz1[1] = (z+2)*(z+1);
    fak[1] = 1;
    pn[1] = p;
    xn[1] = x;

    //(*** now uses recurrence relation ***)
    for ( n=2; n<=maxsum; n++ )
    {
        pa1[n] = pa1[n-1]*(a1+n-1);
        pb1[n] = pb1[n-1]*(b1+n-1);
        pb2[n] = pb2[n-1]*(b2+n-1);
        pz1[n] = (z+2*n)*(z+2*n-1)*pz1[n-1];
        fak[n] = fak[n-1]*n;
        pn[n]  = pn[n-1]*p;
        xn[n]  = xn[n-1]*x;
    }

    //(* **** Asymptotic expansion **** *)
    if ( q*r > 5 )
    {
        sa1 = sqr(3-a);
        if ( fabs(a - 0) < eps ) rt = 0;
        else if ( fabs(a - 1) < eps ) rt = 0.5641895835;
        else if ( fabs(a - 4/3.0) < eps ) rt = 0.8335958021;
        else rt = gamma((3-a)/2.)/(gamma(a/2.));
        sc1 = sqrt(M_PI)*exp(((1-a)/2.)*log(4.))*rt;
        sc2 = ((37/6.)+2*(-(13/12.)+(sqr(3-a)/4.)+((5-a)/2.)-(sqr(5-a)/4.)+((a-3)/2.)))/2.;
        sa2 = sc1/gamma(z+1);
        sb1 = sc1*gamma(z-5+2*a)*exp((2*a-6)*log(y));
        sb2 = 2*gamma(z-4+a)*exp((a-5)*log(y));
        sb3 = 2*sc2*gamma(z-4+a)*exp((a-6)*log(y))/(z-5+a);
        sa3 = gamma(z-3)/(sqr(y)*sqr(y)*sqr(y)*2*gamma(z+1));
        sb4 = sqr(y);
        sb5 = 2*sc2*y/(z-4);
        sb6 = sc2*sc2/((z-5)*(z-4));
        f1 = exp(-(z-4+a)*log(1+sqr(y))/2.)*cos((z-4+a)*atan(y));
        f2 = exp(-(z-5+a)*log(1+sqr(y))/2.)*sin((z-5+a)*atan(y));
        f3 = 1+exp(-(z-3)*log(1+4*sqr(y))/2.)*cos((z-3)*atan(2*y));
        f4 =   exp(-(z-4)*log(1+4*sqr(y))/2.)*sin((z-4)*atan(2*y));
        f5 = 1-exp(-(z-5)*log(1+4*sqr(y))/2.)*cos((z-5)*atan(2*y));
        return sa1*(sa2*(sb1-sb2*f1-sb3*f2)+sa3*(sb4*f3+sb5*f4+sb6*f5));
    }

    //(* **** Hypergeometric function **** *)
    hgeo = 0;
    for ( n=0; n<=maxsum; n++ )
    {
        oldhgeo = hgeo;
        if ( (n % 2) == 0 ) sig = 1;
        else sig = -1;
        msum = 0;
        for ( m=0; m<=n; m++ )
        {
            zaehler = pa1[m]*pa1[n-m]*pn[m];
            nenner = pb1[m]*pb1[n-m]*pb2[m]*pb2[n-m]*fak[m]*fak[n-m];
            msum += zaehler/nenner;
        }
        hgeo += msum*pz1[n]*sig*xn[n];
        if (fabs(hgeo-oldhgeo) < eps) break;
    }
    return hgeo;
}
#endif // USE_pgensphere1



#ifdef USE_hypcoreshellm
/**
 * @brief SasCalculation::hypcoreshellm
 * @param alf1
 * @param alf2
 * @param p1
 * @param p2
 * @param r
 * @param sigma
 * @param q
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::hypcoreshellm(double alf1, double alf2, double p1, double p2, double r, double q) const    //{NV}
{
    const int dim    = 3;            //(* dimension *)
    const int maxsum = 250;          //(* number of summations *)

    double pa11[maxsum+1], pa12[maxsum+1], pb11[maxsum+1], pb12[maxsum+1], pb21[maxsum+1], pb22[maxsum+1],
            pz1[maxsum+1], fak[maxsum+1], ex1[maxsum+1], ex2[maxsum+1], xn[maxsum+1]; // array[0..maxsum] of extended;
    double c1[6], c2[6], d1[6], d2[6], e1[6], e2[6]; // array[0..5] of extended;
    int n,nn,m,mm,sig;
    double z,y,x,msum,
           a11,a12,b11,b12,b21,b22,ga11,ga12,gb11,gb12,gb21,gb22,
           zaehler,nenner,hgeo,oldhgeo,
           s11,s12,s1a1,s1a2,s21,s22,s31,s32,
           t11,t12,t21,t22,t31,t32,t3a1,t3a2,
           f2;

    a11 = (dim-alf1)/2.;
    a12 = (dim-alf2)/2.;
    b11 = dim/2.;
    b12 = dim/2.;
    b21 = (dim+2-alf1)/2.;
    b22 = (dim+2-alf2)/2.;
    z = (1-sqr(params.sigma))/sqr(params.sigma);
    y = q*r/(z+1);
    x = sqr(y/2.);
    ga11 = gamma(a11);
    ga12 = gamma(a12);
    gb11 = gamma(b11);
    gb12 = gamma(b12);
    gb21 = gamma(b21);
    gb22 = gamma(b22);
    //gz1 = gamma(z+1);
    pa11[0] = 1;
    pa12[0] = 1;
    pb11[0] = 1;
    pb12[0] = 1;
    pb21[0] = 1;
    pb22[0] = 1;
    pz1[0] = 1;
    fak[0] = 1;
    ex1[0] = 1;
    ex2[0] = 1;
    xn[0] = 1;

    //(*** simplification ***)
    pa11[1] = a11;
    pa12[1] = a12;
    pb11[1] = b11;
    pb12[1] = b12;
    pb21[1] = b21;
    pb22[1] = b22;
    pz1[1] = (z+2)*(z+1);
    fak[1] = 1;
    ex1[1] = p1*p1/(p2*p2);
    ex2[1] = p2*p2;
    xn[1] = x;

    for ( n=2; n<=maxsum; n++ )
    {
        //pa11[n] = gamma(a11+n)/ga11;
        pa11[n] = pa11[n-1]*(a11+n-1);
        //pa12[n] = gamma(a12+n)/ga12;
        pa12[n] = pa12[n-1]*(a12+n-1);
        //pb11[n] = gamma(b11+n)/gb11;
        pb11[n] = pb11[n-1]*(b11+n-1);
        //pb12[n] = gamma(b12+n)/gb12;
        pb12[n] = pb12[n-1]*(b12+n-1);
        //pb21[n] = gamma(b21+n)/gb21;
        pb21[n] = pb21[n-1]*(b21+n-1);
        //pb22[n] = gamma(b22+n)/gb22;
        pb22[n] = pb22[n-1]*(b22+n-1);
        //pz1[n] = gamma(z+1+2*n)/gz1;
        pz1[n] = (z+2*n)*(z+2*n-1)*pz1[n-1];
        fak[n] = fak[n-1]*n;
        ex1[n] = ex1[n-1]*ex1[1];
        ex2[n] = ex2[n-1]*ex2[1];
        xn[n] = xn[n-1]*x;
    }

    //(* **** asymptotic expansion **** *)
    if ( q*r > 5 )
    {
        s11 = a11-b11-b21;
        s12 = a12-b12-b22;
        s1a1 = (1/2.)+s11;
        s1a2 = (1/2.)+s12;
        s21 = a11*a11-b11*b11-b21*b21;
        s22 = a12*a12-b12*b12-b22*b22;
        s31 = a11*a11*a11-b11*b11*b11-b21*b21*b21;
        s32 = a12*a12*a12-b12*b12*b12-b22*b22*b22;
        t11 = (1/3.)+s11-s1a1*s1a1-2*(-(1/3.)-s11+s21);
        t12 = (1/3.)+s12-s1a2*s1a2-2*(-(1/3.)-s12+s22);
        t21 = t11*t11/8.+(1/6.)*((3/2.)*s1a1*s1a1-s1a1*s1a1*s1a1-(1/2.)*s1a1+4*((1/2.)*s11-(3/2.)*s21+s31));
        t22 = t12*t12/8.+(1/6.)*((3/2.)*s1a2*s1a2-s1a2*s1a2*s1a2-(1/2.)*s1a2+4*((1/2.)*s12-(3/2.)*s22+s32));
        //t3a1 = gamma(a11);
        t3a1 = ga11;
        //t3a2 = gamma(a12);
        t3a2 = ga12;
        t31 = t3a1/(gamma(b11-a11)*gamma(b21-a11));
        t32 = t3a2/(gamma(b12-a12)*gamma(b22-a12));
        c1[1] = exp(-s1a1*log(4.)/2.)*exp(-log(M_PI)/2.);
        c2[1] = exp(-s1a2*log(4.)/2.)*exp(-log(M_PI)/2.);
        c1[2] = t11*c1[1]/2.;
        c2[2] = t12*c2[1]/2.;
        c1[3] = t21*c1[1];
        c2[3] = t22*c2[1];
        c1[4] = exp(a11*log(4.))*t31;
        c2[4] = exp(a12*log(4.))*t32;
        c1[5] = exp((a11+1)*log(4.))*a11*(1+a11-b11)*(1+a11-b21)*t31;
        c2[5] = exp((a12+1)*log(4.))*a12*(1+a12-b12)*(1+a12-b22)*t32;
        d1[1] = ((1/2)+s11)*M_PI/2.;
        d2[1] = ((1/2)+s12)*M_PI/2.;
        d1[2] = ((3/2)+s11)*M_PI/2.;
        d2[2] = ((3/2)+s12)*M_PI/2.;
        d1[3] = ((5/2)+s11)*M_PI/2.;
        d2[3] = ((5/2)+s12)*M_PI/2.;
        d1[4] = -1000;
        d2[4] = -1000;
        d1[5] = -1000;
        d2[5] = -1000;
        e1[1] = s11+1/2.;
        e2[1] = s12+1/2.;
        e1[2] = s11-1/2.;
        e2[2] = s12-1/2.;
        e1[3] = s11-3/2.;
        e2[3] = s12-3/2.;
        e1[4] = -2*a11;
        e2[4] = -2*a12;
        e1[5] = -2*a11-2;
        e2[5] = -2*a12-2;
        f2 = 0;
        for ( nn=1; nn<=5; nn++ )
            for ( mm=1; mm<=5; mm++ )
                f2 += c1[nn]*c2[mm]*cosavm(d1[nn],d2[mm],e1[nn],e2[mm],p1,p2,y,z);
        //f2 = f2*gamma(b11)*gamma(b21)*gamma(b12)*gamma(b22)/(gamma(a11)*gamma(a12));
        return f2*gb11*gb21*gb12*gb22/(ga11*ga12);
    }

    //(* **** hypergeometric function **** *)
    hgeo = 0;
    for ( n=0; n<=maxsum; n++ )
    {
        oldhgeo = hgeo;
        if ( (n % 2) == 0 ) sig = 1;
        else sig = -1;
        msum = 0;
        for ( m=0; m<=n; m++ )
        {
            //zaehler = pa11[m]*pa12[n-m]*exp(2*m*log(p1))*exp(2*(n-m)*log(p2));
            zaehler = pa11[m]*pa12[n-m]*ex1[m]*ex2[n];
            nenner = pb11[m]*pb12[n-m]*pb21[m]*pb22[n-m]*fak[m]*fak[n-m];
            msum += zaehler/nenner;
        }
        //hgeo = hgeo+msum*pz1[n]*sig*exp(n*log(x));
        hgeo += msum*pz1[n]*sig*xn[n];
        if (fabs(hgeo-oldhgeo) < eps) break;
    }
    return hgeo;
}
#endif // USE_hypcoreshellm



#ifdef USE_trapcube
//(* *********************** Cube- Integral ****************************** *)
/**
 * @brief SasCalculation::trapcube
 * @param ac
 * @param bc
 * @param a
 * @param sigma
 * @param pf
 * @param q
 * @param sc
 * @param nc
 * @param trapzditc
 */
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::trapcube(double ac, double bc, double a, bool pf, double q,
                              double &sc, int nc,int &trapzditc) const          //{NV}
{
    long jc;
    double fa,fb,fx,xc,tnmc,sumc,delc;

    if ( nc == 1 )
    {
        fa=cubealf(a,ac,pf,q)*sin(ac);
        fb=cubealf(a,bc,pf,q)*sin(bc);
        sc=0.5*(bc-ac)*(fa+fb);
        trapzditc=1;
    }
    else
    {
        tnmc=trapzditc;
        delc=(bc-ac)/tnmc;
        xc=ac+0.5*delc;
        sumc=0.0;
        for ( jc=1; jc<trapzditc; jc++ )
        {
            fx=cubealf(a,xc,pf,q)*sin(xc);
            sumc=sumc+fx;
            xc=xc+delc;
        }
        sc=0.5*(sc+(bc-ac)*sumc/tnmc);
        trapzditc=2*trapzditc;
    }
}
#endif // USE_trapcube



#ifdef USE_cubealf
//(*** ********* polydisperse cubes ********* *)
/**
 * @brief SasCalculation::cubealf
 * @param a
 * @param alfa
 * @param sigma
 * @param pf
 * @param q
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::cubealf(double a, double alfa, bool pf, double q) const        //{NV}
{
    const int maxit = 7; // früher  10;

    double intrc,oldintrc;
    int    itrc,trapzditrc;

    intrc=0;
    //(*** integral over beta = 0 .. 2*pi ***)
    //(*** integrate over beta = 0 .. pi/2 due to 4-fold symmetry ***)
    trapcubebet(0+eps3,M_PI/2.-eps3,a,alfa,pf,q,intrc,1,trapzditrc);
    oldintrc=intrc;
    for ( itrc=2; itrc<=maxit; itrc++ )
    {
        trapcubebet(0+eps3,M_PI/2.-eps3,a,alfa,pf,q,intrc,itrc,trapzditrc);
        if (oldintrc != 0.0 && fabs(1-intrc/oldintrc) < eps3) break;
        oldintrc=intrc;
    }
    return intrc*2/M_PI;   //(* normalization *)
}
#endif // USE_cubealf



#ifdef USE_trapcubebet
//(* *********************** Cube- Integral ****************************** *)
/**
 * @brief SasCalculation::trapcubebet
 * @param ac
 * @param bc
 * @param a
 * @param alfa
 * @param pf
 * @param q
 * @param sc
 * @param nc
 * @param trapzditc
 */
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::trapcubebet(double ac, double bc, double a, double alfa, bool pf, double q,
                                 double &sc, int nc, int &trapzditc) const                                  //{NV}
{
    long jc;
    double fa,fb,fx,xc,tnmc,sumc,delc;

    if ( nc == 1 )
    {
        fa=cube(a,alfa,ac,params.sigma,pf,q);
        fb=cube(a,alfa,bc,params.sigma,pf,q);
        sc=0.5*(bc-ac)*(fa+fb);
        trapzditc=1;
    }
    else
    {
        tnmc=trapzditc;
        delc=(bc-ac)/tnmc;
        xc=ac+0.5*delc;
        sumc=0.0;
        for ( jc=1; jc<=trapzditc; jc++ )
        {
            fx=cube(a,alfa,xc,params.sigma,pf,q);
            sumc += fx;
            xc += delc;
        }
        sc = 0.5*(sc+(bc-ac)*sumc/tnmc);
        trapzditc = 2*trapzditc;
    }
}
#endif // USE_trapcubebet



#ifdef USE_cube
//(* *************** polydisperse cube (alpha, beta ) ****************************)
//(* input in polar coordinates: q, alfa, beta *)
//(* alfa = polar angle 0..pi, beta = azimuthal angle -pi .. pi *)
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::cube(double a, double alfa, double beta, double sigma, bool pf, double q) const        //{NV}
{
    double z,b,c,xa,xb,xc,fave,i1,i2,i3;

    a=a/2.;       //(*** a in calculation is half base length ***)
    z=(1-sqr(sigma))/sqr(sigma);    // TODO: sigma -> params ???
    b=a;
    c=a;
    xa=q*a*sin(alfa)*cos(beta)+eps9;
    xb=q*b*sin(alfa)*sin(beta)+eps9;
    xc=q*c*cos(alfa)+eps9;
    //ia:=sin(xa)*sin(xb)*sin(xc)/(xa*xb*xc);
    if ( pf /*==1*/ )
    {
        if (xa < 0.001) i1=1; else i1=szave(3,1,xa,-1,z);
        if (xb < 0.001) i2=1; else i2=szave(3,1,xb,-1,z);
        if (xc < 0.001) i3=1; else i3=szave(3,1,xc,-1,z);
        fave=i1*i2*i3;
        //fave:=szave(3,1,xa,-1,z)*szave(3,1,xb,-1,z)*szave(3,1,xc,-1,z);
        return fave*fave;
    }
    // if pf == 0
    if (xa < 0.001) i1=1; else i1=szave(5,1,xa,-2,z);
    if (xb < 0.001) i2=1; else i2=szave(5,1,xb,-2,z);
    if (xc < 0.001) i3=1; else i3=szave(5,1,xc,-2,z);
    return i1*i2*i3;
    // iave:=szave(5,1,xa,-2,z)*szave(5,1,xb,-2,z)*szave(5,1,xc,-2,z);
}
#endif // USE_cube




//{NV} - neue Routinen ------------------------------------------------------------------------


//(* *********************** Romberg integration ****************************** *)

// C++ kann leider nicht so einfach lokale Funktionen nutzen. Somit müssen die lokalen Pascal-Routinen hier auch global sein.
// Das bedeutet, dass es auch Variablen / Konstanten gibt, die hier global definiert werden. Diese werden später richtig in
// die Klasse integriert (und dabei passend benamt).

//static const int jmax=20;
//static const int k=5;


#ifdef USE_trapezpqr
//(*** integration routine use trapezoidal rule ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::trapezpqr( double at, double bt, double /*l*/, double r, double p1, double rho, double /*alfa*/, double /*sigmal*/, double q,
                          int /*i0*/, int /*i1*/, int i3, int i4, double &pq, double &nn, int n, int &trapzdit ) const
{
    int j, dim, m;
    double xt, tnm, sump, sumn, del;
    double pa, pb, px, z, pa1, pb1, px1, na, nb, nx;
    double pa1c, pa1s, pa1cs, pb1c, pb1s, pb1cs, px1c, px1s, px1cs;

    z = (1-sqr(params.sigma))/sqr(params.sigma);
    dim = i4;
    m = 2*dim;

    if ( n == 1 )
    {
        switch ( i3 )
        {
        case 0:   //(* sphere *)
            pa1 = 3*(sin(q*at)-q*at*cos(q*at))/pow(q*at,3);
            pa = sqr(pa1)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1 = 3*(sin(q*bt)-q*bt*cos(q*bt))/pow(q*bt,3);
            pb = sqr(pb1)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 1:   //(* cylinder *)
            pa1 = 2*bessel10(q*at)/(q*at);
            pa = sqr(pa1)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1 = 2*bessel10(q*bt)/(q*bt);
            pb = sqr(pb1)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 2:   //(* disk *)
            pa1 = sin(q*at/2.)/(q*at/2.);
            pa = sqr(pa1)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1 = sin(q*bt/2.)/(q*bt/2.);
            pb = sqr(pb1)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 3:   //(* core/shell sphere *)
            pa1c = 3*(sin(q*p1*at)-q*p1*at*cos(q*p1*at))/pow(q*p1*at,3);
            pa1s = 3*(sin(q*at)-q*at*cos(q*at))/pow(q*at,3);
            pa1cs = (pow(p1,dim)*(1-rho)*pa1c+rho*pa1s)/(pow(p1,dim)*(1-rho)+rho);
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 3*(sin(q*p1*bt)-q*p1*bt*cos(q*p1*bt))/pow(q*p1*bt,3);
            pb1s = 3*(sin(q*bt)-q*bt*cos(q*bt))/pow(q*bt,3);
            pb1cs = (pow(p1,dim)*(1-rho)*pb1c+rho*pb1s)/(pow(p1,dim)*(1-rho)+rho);
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 4:   //(* core/shell cylinder *)
            pa1c = 2*bessel10(q*p1*at)/(q*p1*at);
            pa1s = 2*bessel10(q*at)/(q*at);
            pa1cs = (pow(p1,dim)*(1-rho)*pa1c+rho*pa1s)/(pow(p1,dim)*(1-rho)+rho);
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 2*bessel10(q*p1*bt)/(q*p1*bt);
            pb1s = 2*bessel10(q*bt)/(q*bt);
            pb1cs = (pow(p1,dim)*(1-rho)*pb1c+rho*pb1s)/(pow(p1,dim)*(1-rho)+rho);
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 5:   //(* core/shell disk *)
            pa1c = sin(q*p1*at/2.)/(q*p1*at/2.);
            pa1s = sin(q*at/2.)/(q*at/2.);
            pa1cs = (pow(p1,dim)*(1-rho)*pa1c+rho*pa1s)/(pow(p1,dim)*(1-rho)+rho);
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = sin(q*p1*bt/2.)/(q*p1*bt/2.);
            pb1s = sin(q*bt/2.)/(q*bt/2.);
            pb1cs = (pow(p1,dim)*(1-rho)*pb1c+rho*pb1s)/(pow(p1,dim)*(1-rho)+rho);
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 6:   //(* core/shell sphere term #1 *)
            pa1c = 3*(sin(q*p1*at)-q*p1*at*cos(q*p1*at))/pow(q*p1*at,3);
            pa1cs = pow(p1,dim)*(1-rho)*pa1c;
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 3*(sin(q*p1*bt)-q*p1*bt*cos(q*p1*bt))/pow(q*p1*bt,3);
            pb1cs = pow(p1,dim)*(1-rho)*pb1c;
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 7:   //(* core/shell sphere term #2 *)
            pa1c = 9*(sin(q*p1*at)-q*p1*at*cos(q*p1*at))*(sin(q*at)-q*at*cos(q*at))/(pow(q*p1*at,3)*pow(q*at,3));
            pa1cs = 2*(pow(p1,dim)*(1-rho)*rho)*pa1c;
            pa = pa1cs*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 9*(sin(q*p1*bt)-q*p1*bt*cos(q*p1*bt))*(sin(q*bt)-q*bt*cos(q*bt))/(pow(q*p1*bt,3)*pow(q*bt,3));
            pb1cs = 2*(pow(p1,dim)*(1-rho)*rho)*pb1c;
            pb = pb1cs*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 8:   //(* core/shell sphere term #3 *)
            pa1c = 3*(sin(q*at)-q*at*cos(q*at))/pow(q*at,3);
            pa1cs = rho*pa1c;
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 3*(sin(q*bt)-q*bt*cos(q*bt))/pow(q*bt,3);
            pb1cs = rho*pb1c;
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 9:   //(* core/shell cylinder term #1 *)
            pa1c = 2*bessel10(q*p1*at)/(q*p1*at);
            pa1cs = pow(p1,dim)*(1-rho)*pa1c;
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 2*bessel10(q*p1*bt)/(q*p1*bt);
            pb1cs = pow(p1,dim)*(1-rho)*pb1c;
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 10:   //(* core/shell cylinder term #2 *)
            pa1c = 4*bessel10(q*p1*at)*bessel10(q*at)/(q*at*q*p1*at);
            pa1cs = 2*(pow(p1,dim)*(1-rho)*rho)*pa1c;
            pa = pa1cs*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 4*bessel10(q*p1*bt)*bessel10(q*bt)/(q*bt*q*p1*bt);
            pb1cs = 2*(pow(p1,dim)*(1-rho)*rho)*pb1c;
            pb = pb1cs*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 11:   //(* core/shell cylinder term #3 *)
            pa1c = 2*bessel10(q*at)/(q*at);
            pa1cs = rho*pa1c;
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = 2*bessel10(q*bt)/(q*bt);
            pb1cs = rho*pb1c;
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 12:   //(* core/shell disk term #1 *)
            pa1c = sin(q*p1*at)/(q*p1*at);
            pa1cs = pow(p1,dim)*(1-rho)*pa1c;
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = sin(q*p1*bt)/(q*p1*bt);
            pb1cs = pow(p1,dim)*(1-rho)*pb1c;
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 13:   //(* core/shell disk term #2 *)
            pa1c = sin(q*p1*at)*sin(q*at)/(q*at*q*p1*at);
            pa1cs = 2*(pow(p1,dim)*(1-rho)*rho)*pa1c;
            pa = pa1cs*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = sin(q*p1*bt)*sin(q*bt)/(q*bt*q*p1*bt);
            pb1cs = 2*(pow(p1,dim)*(1-rho)*rho)*pb1c;
            pb = pb1cs*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        case 14:   //(* core/shell disk term #3 *)
            pa1c = sin(q*at)/(q*at);
            pa1cs = rho*pa1c;
            pa = sqr(pa1cs)*pow(at,z+m)*exp(-(z+1)*at/r);
            pb1c = sin(q*bt)/(q*bt);
            pb1cs = rho*pb1c;
            pb = sqr(pb1cs)*pow(bt,z+m)*exp(-(z+1)*bt/r);
            na = pow(at,z+m)*exp(-(z+1)*at/r);
            nb = pow(bt,z+m)*exp(-(z+1)*bt/r);
            break;
        } // switch ( i3 )
        pq = 0.5*(bt-at)*(pa+pb);
        nn = 0.5*(bt-at)*(na+nb);
        trapzdit = 1;
    } // if ( n == 1 )
    else
    {
        tnm = trapzdit;
        del = (bt-at)/tnm;
        xt = at+0.5*del;
        sump = 0.0;
        sumn = 0.0;
        for ( j=1; j<=trapzdit; j++ )
        {
            switch ( i3 )
            {
            case 0: //(* sphere *)
                px1 = 3*(sin(q*xt)-q*xt*cos(q*xt))/pow(q*xt,3);
                px = sqr(px1)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 1: //(* cylinder *)
                px1 = 2*bessel10(q*xt)/(q*xt);
                px = sqr(px1)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 2: //(* disk *)
                px1 = sin(q*xt/2.)/(q*xt/2.);
                px = sqr(px1)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 3:   //(* core/shell sphere *)
                px1c = 3*(sin(q*p1*xt)-q*p1*xt*cos(q*p1*xt))/pow(q*p1*xt,3);
                px1s = 3*(sin(q*xt)-q*xt*cos(q*xt))/pow(q*xt,3);
                px1cs = (pow(p1,dim)*(1-rho)*px1c+rho*px1s)/(pow(p1,dim)*(1-rho)+rho);
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 4:   //(* core/shell cylinder *)
                px1c = 2*bessel10(q*p1*xt)/(q*p1*xt);
                px1s = 2*bessel10(q*xt)/(q*xt);
                px1cs = (pow(p1,dim)*(1-rho)*px1c+rho*px1s)/(pow(p1,dim)*(1-rho)+rho);
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 5:   //(* core/shell disk *)
                px1c = sin(q*p1*xt/2.)/(q*p1*xt/2.);
                px1s = sin(q*xt/2.)/(q*xt/2.);
                px1cs = (pow(p1,dim)*(1-rho)*px1c+rho*px1s)/(pow(p1,dim)*(1-rho)+rho);
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 6:   //(* core/shell sphere term #1 *)
                px1c = 3*(sin(q*p1*xt)-q*p1*xt*cos(q*p1*xt))/pow(q*p1*xt,3);
                px1cs = pow(p1,dim)*(1-rho)*px1c;
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 7:   //(* core/shell sphere term #2 *)
                px1c = 9*(sin(q*p1*xt)-q*p1*xt*cos(q*p1*xt))*(sin(q*xt)-q*xt*cos(q*xt))/(pow(q*p1*xt,3)*pow(q*xt,3));
                px1cs = 2*(pow(p1,dim)*(1-rho)*rho)*px1c;
                px = px1cs*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 8:   //(* core/shell sphere term #3 *)
                px1c = 3*(sin(q*xt)-q*xt*cos(q*xt))/pow(q*xt,3);
                px1cs = rho*px1c;
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 9:   //(* core/shell cylinder term #1 *)
                px1c = 2*bessel10(q*p1*xt)/(q*p1*xt);
                px1cs = pow(p1,dim)*(1-rho)*px1c;
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 10:   //(* core/shell cylinder term #2 *)
                px1c = 4*bessel10(q*p1*xt)*bessel10(q*xt)/(q*xt*q*p1*xt);
                px1cs = 2*(pow(p1,dim)*(1-rho)*rho)*px1c;
                px = px1cs*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 11:   //(* core/shell cylinder term #3 *)
                px1c = 2*bessel10(q*xt)/(q*xt);
                px1cs = rho*px1c;
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 12:   //(* core/shell disk term #1 *)
                px1c = sin(q*p1*xt)/(q*p1*xt);
                px1cs = pow(p1,dim)*(1-rho)*px1c;
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 13:   //(* core/shell disk term #2 *)
                px1c = sin(q*p1*xt)*sin(q*xt)/(q*p1*xt*q*xt);
                px1cs = 2*(pow(p1,dim)*(1-rho)*rho)*px1c;
                px = px1cs*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            case 14:   //(* core/shell disk term #3 *)
                px1c = sin(q*xt)/(q*xt);
                px1cs = rho*px1c;
                px = sqr(px1cs)*pow(xt,z+m)*exp(-(z+1)*xt/r);
                nx = pow(xt,z+m)*exp(-(z+1)*xt/r);
                break;
            } // switch ( i3 )
            sump += px;
            sumn += nx;
            xt   += del;
        } // for j
        pq = 0.5*(pq+(bt-at)*sump/tnm);
        nn = 0.5*(nn+(bt-at)*sumn/tnm);
        trapzdit = 2*trapzdit;
    } // if ( n==1 ) else ...
} /* trapezpqr() */
#endif // USE_trapezpqr



#ifdef USE_midpntchi
//(*** integration routine use trapezoidal rule with midpoints ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::midpntchi( double a, double b, double l, double r, double p1, double rho, double alfa, double sigmal,
                          double sigma, double cc0, double cc1, double cc2, double cc3, double cc4, double cc5, double cc6,
                          double cc7, double cc8, double delta, double phi, double theta, double qxn, double qyn, double qzn,
                          double q, int i0, int i1, int i3, int i4,  double &sp, double &sf, int n, int &midpntit ) const
{
    int j;
    double x, x1, tnm, sump, sumf, del, ddel, spp=0, sff=0; // =0 to avoid compiler warnings
    DD( qDebug() << "midpntchi nyi" );
    if ( n == 1 )
    {
        x1 = 0.5*(a+b);
        switch ( i4 )
        {
        case 0:
            //TBD- orientcylchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x1,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
            break;
        case 1:
            //TBD- orientdiskchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x1,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
            break;
        case 2:
            //TBD- orientcubechi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x1,qxn,qyn,qzn,q,i0,i1,spp,sff);
            break;
        case 3:
            //TBD- orientlatchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x1,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
            break;
        }
        sp = (b-a)*spp; // Werte werden in den Routinen oben berechnet ...
        sf = (b-a)*sff;
        midpntit = 1;
    } // if n == 1
    else
    {
        tnm = midpntit;
        del = (b-a)/(3.0*tnm);
        ddel = del+del;
        x = a+0.5*del;
        sump = 0.0;
        sumf = 0.0;
        for ( j=1; j<=midpntit; j++ )
        {
            switch ( i4 )
            {
            case 0:
                //TBD- orientcylchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            case 1:
                //TBD- orientdiskchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            case 2:
                //TBD- orientcubechi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,qxn,qyn,qzn,q,i0,i1,spp,sff);
                break;
            case 3:
                //TBD- orientlatchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            }
            sump += spp; // Werte werden in den Routinen oben berechnet ...
            sumf += sff;
            x    += ddel;
            switch ( i4 )
            {
            case 0:
                //TBD- orientcylchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            case 1:
                //TBD- orientdiskchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            case 2:
                //TBD- orientcubechi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,qxn,qyn,qzn,q,i0,i1,spp,sff);
                break;
            case 3:
                //TBD- orientlatchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            }
            sump += spp; // Werte werden in den Routinen oben berechnet ...
            sumf += sff;
            x    += del;
        } // for j
        sp = (sp+(b-a)*sump/tnm)/3.0;
        sf = (sf+(b-a)*sumf/tnm)/3.0;
        midpntit = 3*midpntit;
    } // if ( n==1 ) else ...
} /* midpntchi() */
#endif // USE_midpntchi



#ifdef USE_qrombpq
//(* returns integral in the limits a and b *)
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::qrombpq( double l, double r, double p1, double rho, double alfa, double sigmal, double q,
                        int /*maxit*/, int i0, int i1, int i3, int i4,
                        double &pq ) const
{
    const int jmax=20;
    //const int jmaxp=21;  --> im h File
    const int k=5;
    //const int np=5;  --> im h File

    //type RealArrayJMAXP=array[1..jmaxp] of extended;
    //     RealArrayNP=array[1..np] of extended;

    int i, j; //, m;
    double dss; //, delmax;
    double hp[jmaxp+1], hn[jmaxp+1], sp[jmaxp+1], sn[jmaxp+1]; // ^RealArrayJMAXP;
    double cp[np+1], cn[np+1], dp[np+1], dn[np+1]; // ^RealArrayNP;
    double alim, blim;     //(* integration limits *)
    double /*z, norm,*/ nn;
    int trapzdit;

    //(* integration limits and prefactors for SZ-Distribution *)
    alim = r*1e-20;
    blim = r*(1+50*params.sigma);
    //z = (1-sigma*sigma)/(sigma*sigma);
    //m = 2*i4;
    //norm = pow(z+1,z+m+1)/(pow(r,z+m+1)*gamma(z+m+1));

    //new(hp); new(hn);
    //new(sp); new(sn);
    //new(cp); new(cn);
    //new(dp); new(dn);
    hp[1] = 1.0;
    hn[1] = 1.0;
    for ( j=1; j<=jmax; j++ )
    {
        trapezpqr(alim,blim,l,r,p1,rho,alfa,sigmal,q,i0,i1,i3,i4,sp[j],sn[j],j,trapzdit);
        //trapzd(a,b,rr,z,q,s^[j],j);
        //midpnt(a,b,rr,z,q,sp^[j],sf^[j],j);
        //midpntchi(alim,blim,l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,i4,sp^[j],sf^[j],j);
        if ( j >= k )
        {
            for ( i=1; i<=k; i++ )
            {
                cp[i] = hp[j-k+i];
                cn[i] = hn[j-k+i];
                dp[i] = sp[j-k+i];
                dn[i] = sn[j-k+i];
            }
            polint(cp,dp,k,0.0,pq,dss,"qrombpq(a)");
            polint(cn,dn,k,0.0,nn,dss,"qrombpq(b)");
            if ( fabs(dss) < eps3*fabs(pq) ) break;
        } // if j >= k
        sp[j+1] = sp[j];
        sn[j+1] = sn[j];
        hp[j+1] = 0.25*hp[j];
        hn[j+1] = 0.25*hn[j];
    } // for j

    pq = pq/nn;
    //dispose(dp); dispose(dn);
    //dispose(cp); dispose(cn);
    //dispose(sp); dispose(sn);
    //dispose(hp); dispose(hn);
}
#endif // USE_qrombpq



#ifdef USE_qrombpi2


//(*** function, over which is integrated ***)
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::qrombpi2_funcasy1( double beta, double epsi, double r, double sigma, double d, double q ) const
{

    double arg11,nen11,arg12,nen12,arg13,nen13,z,zz,x1,x1z;//,x12z;
    double v,ee0,ee1,pz2v,pz2v1,pz2v2,preg1,F12as1z,F12as2z,F12as3z;//,F12asz,F12asz0;
    int    dim;

    z=(1-sqr(sigma))/sqr(sigma);
    zz=z;
    dim=trunc(d);
    z=z+2*dim;          //(* z-average *)

    //r:=r*sqrt(sin(beta)*sin(beta)+epsi*epsi*cos(beta)*cos(beta));
    r=r*sqrt(1+(epsi*epsi-1)*beta*beta);
    x1=q*r;
    x1z=x1/(2*(zz+1));
    //x12z=x1z*x1z;

    switch ( dim )
    {
    case 1:
        v=-1;
        ee0=1;
        ee1=0;
        pz2v=1/(z*(z-1));
        pz2v1=pz2v/(z-2);
        pz2v2=pz2v1/(z-3);
        preg1=1/2.0;
        break;
    case 2:
        v=-3/2.0;
        ee0=1;
        ee1=-9/16.0;
        pz2v=1/(z*(z-1)*(z-2));
        pz2v1=pz2v/(z-3);
        pz2v2=pz2v1/(z-4);
        preg1=1/sqrt(M_PI);
        break;
    case 3:
        v=-2;
        ee0=1;
        ee1=-1;
        pz2v=1/(z*(z-1)*(z-2)*(z-3));
        pz2v1=pz2v/(z-4);
        pz2v2=pz2v1/(z-5);
        preg1=3/4.0;
        break;
    } // switch dim

   arg11=(z+2*v+1)*atan(4*x1z);
   nen11=pow(1+16*x1z*x1z,(z+2*v+1)/2.0);
   arg12=(z+2*v)*atan(4*x1z);
   nen12=pow(1+16*x1z*x1z,(z+2*v)/2.0);
   arg13=(z+2*v-1)*atan(4*x1z);
   nen13=pow(1+16*x1z*x1z,(z+2*v-1)/2.0);

   F12as1z=ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);
   F12as2z=2*ee0*ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);
   F12as3z=ee1*ee1*(1/(4*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);
   return preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);
}



#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::qrombpi2_funcasy2( double beta, double epsi, double rn, double p1, double rho2, double sigma,
                                    double d, int i3, double q ) const   //(* rho2 *)
{
    int dim;
    double z,zz/*,xrad,xradp*/,x1z,x1zz/*,x12z*/,x2z,x2zz/*,x22z*/,c1,c2,c3;//,ee0,ee1;
    double /*b1s,*/v,e0,e1,/*gb1s,*/preg1,pz2v,pz2v1,pz2v2,pz4,pz5,pz6;
    double arg11,nen11,arg12,nen12,arg13,nen13,F12as1z,F12as2z,F12as3z,F12asz;//,F121;
    double xijm,arglmz,nenlmz,xijp,arglpz,nenlpz,F122as1z,F122as2z,F122as3z,F122as4z,F122asz;//,F122;
    double arg31,nen31,arg32,nen32,arg33,nen33,F32as1z,F32as2z,F32as3z,F32asz;//,F123;
    double argbm,nenbm,argbp,nenbp,argep,nenep,argem,nenem,arggp,nengp,arggm,nengm;//,argim,nenim,argip,nenip;
    double argc,nenc,argh,nenh,argj,nenj;//,nen;
    double arga,nena,argd,nend,argf,nenf;
    double cc[11],ac[11]; //: Array[1..10] of extended;

    //(* ################################### *)
    dim=round(d);
    z=(1-sqr(sigma))/sqr(sigma);
    zz=z;
    z=z+2*dim;           //(* z-average *)

    //rn:=rn*sqrt(sin(beta)*sin(beta)+epsi*epsi*cos(beta)*cos(beta));
    rn=rn*sqrt(1+(epsi*epsi-1)*beta*beta);

    //xrad=q*rn;
    //xradp=q*p1*rn;
    x1z=q*p1*rn/(2*(zz+1));
    x1zz=2*x1z;
    //x12z=x1z*x1z;
    x2z=q*rn/(2*(zz+1));
    x2zz=2*x2z;
    //x22z=x2z*x2z;

    switch ( dim )
    {
    case 1: //(* disk *)
        //b1s=3/2.0;
        v=-1;
        e0=1;
        e1=0;
        //ee0=1;
        //ee1=0;
        //gb1s=sqrt(M_PI)/2.0;
        preg1=1/2.0;
        pz2v=1/(z*(z-1));
        pz2v1=pz2v/(z-2);
        pz2v2=pz2v1/(z-3);
        break;
    case 2 : //(* cylinder *)
        //b1s=2;
        v=-3/2.0;
        e0=1;
        e1=-9/16.0;
        //ee0=1;
        //ee1=-9/16.0;
        //gb1s=1;
        preg1=1/sqrt(M_PI);
        pz2v=1/(z*(z-1)*(z-2));
        pz2v1=pz2v/(z-3);
        pz2v2=pz2v1/(z-4);
        break;
    case 3: //(* sphere *)
        //b1s=5/2.0;
        v=-2;
        e0=1;
        e1=-1;
        //ee0=1;
        //ee1=-1;
        //gb1s=3*sqrt(M_PI)/4.0;
        preg1=3/4.0;
        pz2v=1/(z*(z-1)*(z-2)*(z-3));
        pz2v1=pz2v/(z-4);
        pz2v2=pz2v1/(z-5);
        break;
    } // switch dim

    c1=sqr(1-rho2)*pow(p1,2*dim);
    c2=2*rho2*(1-rho2)*pow(p1,dim);
    c3=rho2*rho2;

    pz4=1/(z*(z-1)*(z-2)*(z-3));
    pz5=pz4/(z-4);
    pz6=pz5/(z-5);

    switch ( i3 )
    {
    case 11: // (* dim 1,2; term #1 *)
        arg11=(z+2*v+1)*atan(4*x1z);
        nen11=pow(1+16*x1z*x1z,(z+2*v+1)/2.0);
        arg12=(z+2*v)*atan(4*x1z);
        nen12=pow(1+16*x1z*x1z,(z+2*v)/2.0);
        arg13=(z+2*v-1)*atan(4*x1z);
        nen13=pow(1+16*x1z*x1z,(z+2*v-1)/2.0);
        F12as1z=e0*e0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);
        F12as2z=2*e0*e1*(1/(2*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);
        F12as3z=e1*e1*(1/(4*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);
        F12asz=preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);
        return c1*F12asz;

    case 12: // (* dim 1,2; term #2 *)
        xijm=x1z-x2z;
        arglmz=(z+1)*atan(2*xijm);
        nenlmz=pow(1+4*xijm*xijm,(z+1)/2.0);
        xijp=x1z+x2z;
        arglpz=(z+1)*atan(2*xijp);
        nenlpz=pow(1+4*xijp*xijp,(z+1)/2.0);
        F122as1z=e0*e0*pz2v*(cos(arglmz)/nenlmz+cos(M_PI*v)*cos(arglpz)/nenlpz-sin(M_PI*v)*sin(arglpz)/nenlpz);
        F122as2z=e0*e1*(1/(2*x2z))*pz2v1*(-sin(arglmz)/nenlmz+cos(M_PI*(2*v-1)/2.0)*cos(arglpz)/nenlpz-sin(M_PI*(2*v-1)/2.0)*sin(arglpz)/nenlpz);
        F122as3z=e1*e0*(1/(2*x1z))*pz2v1*(sin(arglmz)/nenlmz+cos(M_PI*(2*v-1)/2.0)*cos(arglpz)/nenlpz-sin(M_PI*(2*v-1)/2.0)*sin(arglpz)/nenlpz);
        F122as4z=e1*e1*(1/(4*x1z*x2z))*pz2v2*(cos(arglmz)/nenlmz+cos(M_PI*(v-1))*cos(arglpz)/nenlpz-sin(M_PI*(v-1))*sin(arglpz)/nenlpz);
        F122asz=(F122as1z+F122as2z+F122as3z+F122as4z);
        return c2*preg1*preg1*pow(x1z,v)*pow(x2z,v)*(1/2.0)*F122asz;

    case 13: // (* dim 1,2; term #3 *)
        arg31=(z+2*v+1)*atan(4*x2z);
        nen31=pow(1+16*x2z*x2z,(z+2*v+1)/2.0);
        arg32=(z+2*v)*atan(4*x2z);
        nen32=pow(1+16*x2z*x2z,(z+2*v)/2.0);
        arg33=(z+2*v-1)*atan(4*x2z);
        nen33=pow(1+16*x2z*x2z,(z+2*v-1)/2.0);
        F32as1z=e0*e0*pz2v*(1+cos(M_PI*v)*cos(arg31)/nen31-sin(M_PI*v)*sin(arg31)/nen31);
        F32as2z=2*e0*e1*(1/(2*x2z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(2*v-1)/2.0)*sin(arg32)/nen32);
        F32as3z=e1*e1*(1/(4*x2z*x2z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg33)/nen33-sin(M_PI*(v-1))*sin(arg33)/nen33);
        F32asz=preg1*preg1*pow(x2z,2*v)*(1/2.0)*(F32as1z+F32as2z+F32as3z);
        return c3*F32asz;

    case 14: // (* dim 3; term #1 *)
        cc[3] = sqr(1-rho2)*sqr(p1);
        argc = (z-4+1)*atan(2*x1zz);
        nenc = pow(1+sqr(2*x1zz),(z-4+1)/2.0);
        ac[3] = pz4*pow(x2zz,-4)*(1/2.0)*(1+cos(argc)/nenc);
        cc[8] = -sqr(1-rho2)*2*p1;
        argh = (z-5+1)*atan(2*x1zz);
        nenh = pow(1+sqr(2*x1z),(z-5+1)/2.0);
        ac[8] = pz5*pow(x2zz,-5)*(1/2.0)*sin(argh)/nenh;
        cc[10] = sqr(1-rho2);
        argj = (z-6+1)*atan(2*x1zz);
        nenj = pow(1+sqr(2*x1zz),(z-6+1)/2.0);
        ac[10] = pz6*pow(x2zz,-6)*(1/2.0)*(1-cos(argj)/nenj);
        return 9*(cc[3]*ac[3]+cc[8]*ac[8]+cc[10]*ac[10]);

    case 15: // (* dim 3; term #2 *)
        cc[2] = 2*p1*rho2*(1-rho2);
        argbm = (z-4+1)*atan(x1zz-x2zz);
        nenbm = pow(1+sqr(x1zz-x2zz),(z-4+1)/2.0);
        argbp = (z-4+1)*atan(x1zz+x2zz);
        nenbp = pow(1+sqr(x1zz+x2zz),(z-4+1)/2.0);
        ac[2] = pz4*pow(x2zz,-4)*(1/2.0)*(cos(argbm)/nenbm+cos(argbp)/nenbp);
        cc[5] = -2*p1*rho2*(1-rho2);
        argep = (z-5+1)*atan(x1zz+x2zz);
        nenep = pow(1+sqr(x1zz+x2zz),(z-5+1)/2.0);
        argem = (z-5+1)*atan(x1zz-x2zz);
        nenem = pow(1+sqr(x1zz-x2zz),(z-5+1)/2.0);
        ac[5] = pz5*pow(x2zz,-5)*(1/2.0)*(sin(argep)/nenep-sin(argem)/nenem);
        cc[7] = -2*rho2*(1-rho2) ;
        arggp = (z-5+1)*atan(x2zz+x1zz);
        nengp = pow(1+sqr(x2zz+x1zz),(z-5+1)/2.0);
        arggm = (z-5+1)*atan(x2zz-x1zz);
        nengm = pow(1+sqr(x2zz-x1zz),(z-5+1)/2.0);
        ac[7] = pz5*pow(x2zz,-5)*(1/2.0)*(sin(arggp)/nengp-sin(arggm)/nengm);
        cc[9] = 2*rho2*(1-rho2);
        //argim = (z-6+1)*atan(x1zz-x2zz);
        //nenim = pow(1+sqr(x1zz-x2zz),(z-6+1)/2.0);
        //argip = (z-6+1)*atan(x1zz+x2zz);
        //nenip = pow(1+sqr(x1zz+x2z),(z-6+ 1)/2.0);
        ac[9] = pz6*pow(x2zz,-6)*(1/2.0)*(cos(argbm)/nenbm-cos(argbp)/nenbp);
        return 9*(cc[2]*ac[2]+cc[5]*ac[5]+cc[7]*ac[7]+cc[9]*ac[9]);

    case 16: // (* dim = 3; term #3 *)
        cc[1] = sqr(rho2);
        arga = (z-4+1)*atan(2*x2zz);
        nena = pow(1+sqr(2*x2zz),(z-4+1)/2.0);
        ac[1] = pz4*pow(x2zz,-4)*(1/2.0)*(1+cos(arga)/nena);
        cc[4] = -2*sqr(rho2);
        argd = (z-5+1)*atan(2*x2zz);
        nend = pow(1+sqr(2*x2zz),(z-5+1)/2.0);
        ac[4] = pz5*pow(x2zz,-5)*(1/2.0)*sin(argd)/nend;
        cc[6] = sqr(rho2);
        argf = (z-6+1)*atan(2*x2zz);
        nenf = pow(1+sqr(2*x2zz),(z-6+1)/2.0);
        ac[6] = pz6*pow(x2zz,-6)*(1/2.0)*(1-cos(argf)/nenf);
        return 9*(cc[1]*ac[1]+cc[4]*ac[4]+cc[6]*ac[6]);
    } // switch i3
    return 0.0;
}

/*** integration routine use trapezoidal rule ***/
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::qrombpi2_trapezpqpi2( double at, double bt, double epsi, double r, double p1, double rho, double /*alfa*/,
                                     double /*sigmal*/, double sigma, double q, int /*i0*/, int i1, int i3, int i4,
                                     double &pq, double &nn, int n, int &trapzdit ) const
{
    int    j,dim;//,m;
    double xt,tnm,sump,sumn,del;
    double /*pa,pb,*/px/*,z,pa1,pb1,px1,na,nb*/,nx;
    //double pa1c,pa1s,pa1cs,pb1c,pb1s,pb1cs,px1c,px1s,px1cs;

    //z = (1-sigma*sigma)/(sigma*sigma);
    dim = i4;
    //m = 2*dim;

    if ( n==1 )
    {
        if ( i1==10 )
        {
            /* pq:=0.5*(bt-at)*(sin(at)*funcasy1(at,epsi,r,sigma,dim,q)-sin(bt)*funcasy1(bt,epsi,r,sigma,dim,q)); */
            /* nn:=0.5*(bt-at)*(sin(at)-sin(bt)); */
            pq = 0.5*(bt-at)*(qrombpi2_funcasy1(at,epsi,r,sigma,dim,q)-qrombpi2_funcasy1(bt,epsi,r,sigma,dim,q));
            nn = 0.5*(bt-at)*(bt-at);
        }
        if ( i1==11 )
        {
            pq = 0.5*(bt-at)*(qrombpi2_funcasy2(at,epsi,r,p1,rho,sigma,dim,i3,q)-qrombpi2_funcasy2(bt,epsi,r,p1,rho,sigma,i3,dim,q));
            nn = 0.5*(bt-at)*(bt-at);
        }
        trapzdit = 1;
    }
    else
    {
        tnm = trapzdit;
        del = (bt-at)/tnm;
        xt = at+0.5*del;
        sump = 0.0;
        sumn = 0.0;
        for ( j=1; j<=trapzdit; j++ )
        {
            if ( i1==10 )
            {
                /* px:=sin(xt)*funcasy1(xt,epsi,r,sigma,dim,q); */
                /* nx:=sin(xt); */
                px = qrombpi2_funcasy1(xt,epsi,r,sigma,dim,q);
                nx = 1;
            }
            if ( i1==11 )
            {
                px = qrombpi2_funcasy2(xt,epsi,r,p1,rho,sigma,dim,i3,q);
                nx = 1;
            }
            sump = sump+px;
            sumn = sumn+nx;
            xt = xt+del;
        }
        pq = 0.5*(pq+(bt-at)*sump/tnm);
        nn = 0.5*(nn+(bt-at)*sumn/tnm);
        trapzdit = 2*trapzdit;
    }
}


/*** integration routine use trapezoidal rule with midpoints ***/
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::qrombpi2_midpntchi( double a, double b, double l, double r, double p1, double rho, double alfa, double sigmal,
                                   double sigma, double cc0, double cc1, double cc2, double cc3, double cc4, double cc5,
                                   double cc6, double cc7, double cc8, double delta, double phi, double theta, double qxn,
                                   double qyn, double qzn, double q, int i0, int i1, int i3, int i4,
                                   double &sp, double &sf, int n, int &midpntit ) const
{
#ifndef __CUDACC__
    qDebug() << "qrombpi2_midpntchi  i4" << i4 << "unbekannte Routinen"; // TODO
    exit(1);
#endif
#ifdef undef
    int    j;
    double x,x1,tnm,sump,sumf,del,ddel,spp,sff;

    if ( n==1 )
    {
        x1 = 0.5*(a+b);
        switch ( i4 )
        {
        case 0:
            orientcylchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x1,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
            break;
        case 1:
            //orientdiskchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x1,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
            break;
        case 2:
            //orientcubechi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x1,qxn,qyn,qzn,q,i0,i1,spp,sff);
            break;
        case 3:
            //orientlatchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x1,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
            break;
        }
        sp = (b-a)*spp;
        sf = (b-a)*sff;
        midpntit = 1;
    }
    else
    {
        tnm = midpntit;
        del = (b-a)/(3.0*tnm);
        ddel = del+del;
        x = a+0.5*del;
        sump = 0.0;
        sumf = 0.0;
        for ( j=1; j<=midpntit; j++ )
        {
            // Die lokalen Additionen sind effektiver von der Laufzeit als ein zweites switch()
            switch ( i4 )
            {
            case 0:
                orientcylchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                sump += spp;
                sumf += sff;
                x    += ddel;     // Doppelschritt
                orientcylchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            case 1:
                //orientdiskchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                sump += spp;
                sumf += sff;
                x    += ddel;     // Doppelschritt
                //orientdiskchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            case 2:
                //orientcubechi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,qxn,qyn,qzn,q,i0,i1,spp,sff);
                sump += spp;
                sumf += sff;
                x    += ddel;     // Doppelschritt
                //orientcubechi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,delta,x,qxn,qyn,qzn,q,i0,i1,spp,sff);
                break;
            case 3:
                //orientlatchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                sump += spp;
                sumf += sff;
                x    += ddel;     // Doppelschritt
                //orientlatchi(l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,x,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,spp,sff);
                break;
            }
            sump += spp;
            sumf += sff;
            x    += del;      // Einfachschritt
        }
        sp = (sp+(b-a)*sump/tnm)/3.0;
        sf = (sf+(b-a)*sumf/tnm)/3.0;
        midpntit = 3*midpntit;
    }
#endif
}


//(* *********************** Romberg integration ****************************** *)
//(* returns integral in the limits a and b *)
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::qrombpi2( double epsi, double r, double p1, double rho, double alfa, double sigmal, double sigma, double q,
                         int /*maxit*/, int i0, int i1, int i3, int i4, double &pq ) const
{
    const double eps=1.0e-5;
    const int jmax=20;
    //const int jmaxp=21;  --> im h File
    const int k=5;
    //const int np=5;  --> im h File

    int    i,j;//,m;
    int    /*trapzdit,*/midpntit;
    double dss;//,delmax;
    RealArrayJMAXP hp,hn,sp,sn; // ^RealArrayJMAXP;  =array[1..jmaxp=21] of extended;
    RealArrayNP    cp,cn,dp,dn; // ^RealArrayNP;     =array[1..np=5] of extended;

    double alim,blim;     /* integration limits */;
    double /*z,norm,*/nn;

    /* integration limits and prefactors for SZ-Distribution */
    alim = 0;
    blim = M_PI/2.0;

    /* alternative integration limits */
    alim = 0;
    blim = 1;

    hp[1] = 1.0;
    hn[1] = 1.0;
    for ( j=1; j<=jmax; j++ )
    {
        qrombpi2_trapezpqpi2(alim,blim,epsi,r,p1,rho,alfa,sigmal,sigma,q,i0,i1,i3,i4,sp[j],sn[j],j,midpntit);
        /* trapzd(a,b,rr,z,q,s^[j],j); */
        /* midpnt(a,b,rr,z,q,sp^[j],sf^[j],j); */
        /* midpntchi(alim,blim,l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,delta,phi,theta,qxn,qyn,qzn,q,i0,i1,i3,i4,sp^[j],sf^[j],j); */
        if ( j>=k )
        {
            for ( i=1; i<=k; i++ )
            {
                cp[i] = hp[j-k+i];
                cn[i] = hn[j-k+i];
                dp[i] = sp[j-k+i];
                dn[i] = sn[j-k+i];
            }
            polint(cn,dn,k,0.0,nn,dss,"qrombpi2-1");
            polint(cp,dp,k,0.0,pq,dss,"qrombpi2-2");
            if ( fabs(dss)<eps*fabs(pq) ) break;
        }
        sp[j+1] = sp[j];
        sn[j+1] = sn[j];
        hp[j+1] = 0.25*hp[j];
        hn[j+1] = 0.25*hn[j];
    }
}
#endif // USE_qrombpi2



#ifdef USE_bessel10
//(* ******************* Function BesselJ1 ************************ *)
#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::bessel10( double x ) const
{
    //const double a1 = 0.56249985;
    //const double a2 = 0.21093573;
    //const double a3 = 0.03954289;
    //const double a4 = 0.00443319;
    //const double a5 = 0.00031761;
    //const double a6 = 0.00001109;

    const double  b1 = 2;
    const double  b2 = 16;
    const double  b3 = 384;
    const double  b4 = 18432;
    const double  b5 = 1474560;
    const double  b6 = 176947200;
    const double  b7 = 29727129600;
    const double  b8 = 6658877030400;
    const double  b9 = 1917756584755200;

    double /*term,*/ c1, /*f1, sig,*/ z;
    //int i;

    z = fabs(x);
    //if ( z == 0 ) return 0;

    if ( z > 0 && z <= 5.5 )
        return z/b1-pow(z,3)/b2+pow(z,5)/b3-pow(z,7)/b4+pow(z,9)/b5-pow(z,11)/b6+pow(z,13)/b7-pow(z,15)/b8+pow(z,17)/b9;
    if ( z > 5.5 )
    {
        c1 = sqrt(2/(M_PI*z))*(cos(z-3*M_PI/4.));
        if ( x < 0 ) return -c1;
        if ( x > 0 ) return  c1;
    }
    return 0.0;
}
#endif // USE_bessel10


#ifdef undef
//(*** provides pa(chi), fa(chi); chi (paper) = alpha (routine); for fixed delta ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::orientcylchi( double l, double r, double p1, double rho, double alfa, double sigmal, double sigma, double /*cc0*/,
                             double /*cc1*/, double /*cc2*/, double /*cc3*/, double /*cc4*/, double /*cc5*/, double /*cc6*/, double /*cc7*/, double /*cc8*/,
                             double delta, double alpha, double phi, double theta, double qxn, double qyn, double qzn,
                             double q, int i0, int i1, int i3,  double &pa, double &fa ) const
{
    const int nmax = 300;
    const double delc = 0.0001;

    double radiusi,cosarg,sinarg,/*phienc,beta,xl,*/xlq,xr/*,xrq,pla,pra,fla,fra*/,lim,lim1,zzlen;
    double cylvec[4],cylvecrot[4]; // array[1..3] of extended;
    double conerotarray[4][4]; // array[1..3,1..3] of extended;
    int /*i,*/ii,jj,n; //,m,ncell;
    //float phiax,philiph,philipt,phiin,phiout,lliph,llipt,lin,lout,lcell,frsum,vrsum,flq,frq,vrq,qr,ldel;
    //float phic[101],rad[101],frs[101],vrs[101]; // array[0..100] of real;
    double zlen,xlz,xl2z,b1sl,vl,ee0,ee1,preg1l,pz2vl,pz2v1l,pz2v2l,del,Flsezor,oldFlsezor;
    double z12vl[nmax+1],b1svl[nmax+1],fkv[nmax+1],sum12l[nmax+1]; // Array[0..nmax] of extended;
    double arg11c,nen11c,arg12c,nen12c,arg13c,nen13c,F12as1zc,F12as2zc,F12as3zc,F12aszc;

    if ( i3 == 1 )  //(* crystal lattice with (001)-orientation *)
    {
        cylvecrot[1]=sin(alpha)*sin(delta);
        cylvecrot[2]=cos(alpha)*sin(delta);
        cylvecrot[3]=cos(delta);
        cosarg=qxn*cylvecrot[1]+qyn*cylvecrot[2]+qzn*cylvecrot[3];
    }
    else if ( i3 == 0 )  //(* free rotation *)
    {
        cylvec[1]=sin(theta-delta)*cos(phi);
        cylvec[2]=cos(theta-delta);
        cylvec[3]=-sin(theta-delta)*sin(phi);
        conerotarray[1][1]=cos(alpha)+(1-cos(alpha))*sin(theta)*sin(theta)*cos(phi)*cos(phi);
        conerotarray[1][2]=-sin(theta)*sin(phi)*sin(alpha)+(1-cos(alpha))*cos(theta)*sin(theta)*cos(phi);
        conerotarray[1][3]=-cos(theta)*sin(alpha)-(1-cos(alpha))*sin(theta)*sin(theta)*cos(phi)*sin(phi);
        conerotarray[2][1]=sin(theta)*sin(phi)*sin(alpha)+(1-cos(alpha))*cos(theta)*sin(theta)*cos(phi);
        conerotarray[2][2]=cos(alpha)+(1-cos(alpha))*cos(theta)*cos(theta);
        conerotarray[2][3]=sin(theta)*cos(phi)*sin(alpha)-(1-cos(alpha))*cos(theta)*sin(theta)*sin(phi);
        conerotarray[3][1]=cos(theta)*sin(alpha)-(1-cos(alpha))*sin(theta)*sin(theta)*cos(phi)*sin(phi);
        conerotarray[3][2]=-sin(theta)*cos(phi)*sin(alpha)-(1-cos(alpha))*cos(theta)*sin(theta)*sin(phi);
        conerotarray[3][3]=cos(alpha)+(1-cos(alpha))*sin(theta)*sin(theta)*sin(phi)*sin(phi);
        for ( ii=1; ii<=3; ii++ )
        {
            cylvecrot[ii]=0;
            for ( jj=1; jj<=3; jj++ )
            {
                cylvecrot[ii]=cylvecrot[ii]+conerotarray[ii][jj]*cylvec[jj];
            }
        }
        cosarg=qxn*cylvecrot[1]+qyn*cylvecrot[2]+qzn*cylvecrot[3];
    } // if i3==0
    sinarg=sqrt(1-cosarg*cosarg);

    radiusi=r/p1;

    zlen=(1-sigmal*sigmal)/(sigmal*sigmal);
    zzlen=zlen;
    zlen=zlen+2;   //(* z-aveage for dim=1 *)
    xlq=q*l*cosarg;
    xlz=xlq/(zzlen+1);
    xl2z=xlz*xlz;

    b1sl=3./2.;
    vl=-1;
    ee0=1;
    ee1=0;
    preg1l=1./2.;

    lim=18*exp(-5*sigmal);
    lim1=lim*1.5;

    if ( i0==0 || i0==1 )
    {
        //(*** P(q) ***)
        pz2vl=1/(zlen*(zlen-1));
        pz2v1l=pz2vl/(zlen-2);
        pz2v2l=pz2v1l/(zlen-3);
        //(*** oriented thin rod ***)
        if ( fabs(xlq) < lim1 )
        {
            z12vl[0]=1;
            b1svl[0]=1;
            fkv[0]=1;
            Flsezor=1.0;
            oldFlsezor=1.0;
            for ( n=1; n<=nmax; n++ )
            {
                z12vl[n]=z12vl[n-1]*((zlen+1)-2+2*n)*((zlen+1)-1+2*n);
                b1svl[n]=b1svl[n-1]*(b1sl-1+n);
                fkv[n]=fkv[n-1]*n;
                //sum12l[n]:=0;
                //for m:=0 to n do sum12l[n]:=sum12l[n]+1/(b1svl[m]*b1svl[n-m]*fkv[m]*fkv[n-m]);
                //Flsezor:=Flsezor+power(-xl2z/4,n)*z12vl[n]*sum12l[n];    (******)
                Flsezor=Flsezor+pow(-xl2z,n)*z12vl[n]/((n+1)*b1svl[n]*fkv[n]);
                del=fabs((Flsezor-oldFlsezor)/Flsezor);
                if ( del < delc ) break;
                oldFlsezor=Flsezor;
            } // for n
            pa=Flsezor;
        }

        //(* thin rod oriented asymptote *)
        //(*** term #1 for thin rod ***)
        else // if ( fabs(xlq) >= lim1 )
        {
            arg11c=(zlen+2*vl+1)*atan(xlz);
            nen11c=pow(1+xl2z,(zlen+2*vl+1)/2.);
            arg12c=(zlen+2*vl)*atan(xlz);
            nen12c=pow(1+xl2z,(zlen+2*vl)/2.);
            arg13c=(zlen+2*vl-1)*atan(xlz);
            nen13c=pow(1+xl2z,(zlen+2*vl-1)/2.);

            F12as1zc=ee0*ee0*pz2vl*(1+cos(M_PI*vl)*cos(arg11c)/nen11c-sin(M_PI*vl)*sin(arg11c)/nen11c);
            F12as2zc=2*ee0*ee1*(1/(xlz))*pz2v1l*(cos(M_PI*(2*vl-1)/2)*cos(arg12c)/nen12c-sin(M_PI*(2*vl-1)/2)*sin(arg12c)/nen12c);
            F12as3zc=ee1*ee1*(1/(xl2z))*pz2v2l*(1+cos(M_PI*(vl-1))*cos(arg13c)/nen13c-sin(M_PI*(vl-1))*sin(arg13c)/nen13c);
            F12aszc=preg1l*preg1l*pow(xlz/2.,2*vl)*(1./2.)*(F12as1zc+F12as2zc+F12as3zc);
            pa=F12aszc;
        }

        //(*** F(q) ***)
        if ( i3 == 1 )   //(* crystal lattice, needs F(q) *)
        {
            pz2vl=1/zlen;
            pz2v1l=pz2vl/(zlen-1);
            pz2v2l=pz2v1l/(zlen-2);
            //(*** oriented thin rod ***)
            if ( fabs(xlq) < lim1 )
            {
                z12vl[0]=1;
                b1svl[0]=1;
                fkv[0]=1;
                Flsezor=1.0;
                oldFlsezor=1.0;
                for ( n=1; n<=nmax; n++ )
                {
                    z12vl[n]=z12vl[n-1]*((zlen+1)-2+2*n)*((zlen+1)-1+2*n);
                    b1svl[n]=b1svl[n-1]*(b1sl-1+n);
                    fkv[n]=fkv[n-1]*n;
                    sum12l[n]=0;
                    Flsezor=Flsezor+pow(-xl2z/4.,n)*z12vl[n]/(b1svl[n]*fkv[n]);    //(******)
                    del=fabs((Flsezor-oldFlsezor)/Flsezor);
                    if ( del < delc ) break;
                    oldFlsezor=Flsezor;
                }
                fa=Flsezor*Flsezor;
            }

            //(* thin rod oriented asymptote *)
            //(*** term #1 for thin rod ***)
            else // if ( fabs(xlq) >= lim1 )
            {
                arg11c=(zlen+vl+1)*atan(xlz);
                nen11c=pow(1+xl2z,(zlen+vl+1)/2.);
                arg12c=(zlen+vl)*atan(xlz);
                nen12c=pow(1+xl2z,(zlen+vl)/2.);
                arg13c=(zlen+vl-1)*atan(xlz);
                nen13c=pow(1+xl2z,(zlen+vl-1)/2.);

                F12as1zc=ee0*pz2vl*(cos(M_PI*vl/2.)*cos(arg11c)/nen11c-sin(M_PI*vl/2.)*sin(arg11c)/nen11c);
                F12as2zc=ee1*(1/(xlz))*pz2v1l*(cos(M_PI*(vl-1)/2.)*cos(arg12c)/nen12c-sin(M_PI*(vl-1)/2.)*sin(arg12c)/nen12c);
                F12aszc=preg1l*pow(xlz/2.,vl)*(F12as1zc+F12as2zc);
                fa=F12aszc*F12aszc;
            }
        }
        if ( i3 == 0 ) fa=pa;   //(* no lattice *)
    } //     if ( i0==0 || i0==1 )


    if ( i0 == 0 )     //(* exact calculation *)
    {
        switch ( i1 )
        {
        case 0:
            xr=r*sinarg+eps;
            //TBD- pa=pa*psphered(xr,sigma,2,q);
            //TBD- fa=fa*pspheredf(xr,sigma,2,q);
            break;
        case 1:
            xr=radiusi*sinarg+eps;
            /* @brief SasCalculation::pqcoreshell
            * @param rho1  = 1.0
            * @param rho2  = global var rho => params.rho
            * @param p1    = global var p1 => params.p1
            * @param p2    = 1.0
            * @param alf1  = 0.001
            * @param alf2  = local var
            * @param rn    = global var radiusi => params.radiusi
            * @param d     = 3
            * @param sigma = global var sigma => params.sigma
            * @param q     = local var
            * @return
            */
            //pa=pa*pqcoreshell(1.0,rho,p1,1.0,0.001,0.0001,xr,2,sigma,q);
            params.rho = rho;  // Hier Parameter
            params.p1  = p1;
            params.radiusi = xr;
            params.sigma = sigma;
            pa=pa*pqcoreshell(0.0001,2,q);
            //TBD- fa=fa*pqcoreshellf(1.0,rho,p1,1.0,0.001,0.0001,xr,2,sigma,q);
            break;
        case 2:
            xr=radiusi*sinarg+eps;
            //TBD- pa=pa*pqcoreshellin(1.0,rho,p1,1.0,0.001,alfa,xr,2,sigma,q);
            //TBD- fa=fa*pqcoreshellin(1.0,rho,p1,1.0,0.001,alfa,xr,2,sigma,q);
            break;
        }
    }
}
#endif // undef


#ifndef __CUDACC__
//#define UseStringGrid12
#endif

#define latpar1(a,b) latpar1ptr[latparIDX(a,b, 6)]      // [5000][6], genutzt: 0,1,2,3,4,5
#define latpar2(a,b) latpar2ptr[latparIDX(a,b, 6)]      // [5000][6], genutzt: 0,1,2,3,4,5
#define latpar3(a,b) latpar3ptr[latparIDX(a,b,14)]      // [5000][15], genutzt: 1 bis 12
//#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]      // [5000][15], genutzt: nichts

#ifdef USE_ButtonHKLClick
// Nur in prepareCalculation
void CLASSLIB::ButtonHKLClick( int ltype, int *latpar1ptr, int *latpar2ptr ) const
{
    //const int np=20000;

    int /*i,j,*/ii/*,jj*/,h,k,l,hmax,kmax,lmax;//,index,xpos,ypos,image1width,image1height;
    int /*c0,*/c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,mult,/*multmax,hlev,klev,llev,zlev,zvar,*/ct1;
    // Unbenutzte Variablen: c3, c5, c9, c11, c12, c14 - sie werden aber hier gelassen, da ich sie nicht aus den Zuweisungs-case rausnehmen will
    int ccase,c1case,c2case/*,c3case*/,c4case/*,c5case*/,c6case,c7case,c8case/*,c9case*/,c10case;//,c11case,c12case;
    float a,b,c,alf,bet,gam; //,xmin,xmax,wave,ttheta,sinarg;
    // bet wird gebraucht, wenn UseStringGrid12 definiert ist.
#ifdef UseStringGrid12
    // bet gehört eigendlich hier hin, ist aber an vielen Stellen gesetzt
    int   ct;
    float q, invd;
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
            //ComboBox1D,
            ComboBox2D,
            ComboBoxHexagonal,
            ComboBoxTrigonal=0,     // unused
            ComboBoxOrthorhombic=0, // unused
            ComboBoxMonoclinic=0,   // unused
            ComboBoxTriclinic=0;    // unused

    latpar1(  1,0)=0;
    latpar2(  1,0)=0;

    switch ( ltype )
    {
    // Die jeweiligen "RadioButtonSys*Click(Sender);" setzen ComboBox-Auswahlen und Visibility von GUI-Elementen
    // und können somit hier ignoriert werden.
    case  0: // (* Lam *)
        RadioButtonSys1D=true;
        //RadioButtonSys1DClick(Sender);
        //ComboBox1D=0;
        break;
    case  1: //   (* hex cyl *)
        RadioButtonSys2D=true;
        //RadioButtonSys2DClick(Sender);
        ComboBox2D=12;
        break;
    case  2: //   (* sq Cyl *)
        RadioButtonSys2D=true;
        //RadioButtonSys2DClick(Sender);
        ComboBox2D=10;
        break;
    case  3: //   (* rec cyl *)
        RadioButtonSys2D=true;
        //RadioButtonSys2DClick(Sender);
        ComboBox2D=8;
        break;
    case  4: //   (* BCC *)
        RadioButtonSysCubic=true;
        //RadioButtonSysCubicClick(Sender);
        ComboBoxCubic=7;
        break;
    case  5: //  (* FCC *)
    case 30: //  Generic
        RadioButtonSysCubic=true;
        //RadioButtonSysCubicClick(Sender);
        ComboBoxCubic=12;
        break;
    case  6: //   (* HCP *)
        RadioButtonSysHex=true;
        //RadioButtonSysHexClick(Sender);
        ComboBoxHexagonal=0;
        break;
    case  7: //   (* SC *)
        RadioButtonSysCubic=true;
        //RadioButtonSysCubicClick(Sender);
        ComboBoxCubic=0;
        break;
    case  8: //   (* BCT *)
        RadioButtonSysTetragonal=true;
        //RadioButtonSysTetragonalClick(Sender);
        ComboBoxTetragonal=23;
        break;
    case  9: //   (* Ia3d *)
        RadioButtonSysCubic=true;
        //RadioButtonSysCubicClick(Sender);
        ComboBoxCubic=11;
        break;
    case 10: //  (* Pn3m *)
        RadioButtonSysCubic=true;
        //RadioButtonSysCubicClick(Sender);
        ComboBoxCubic=5;
        break;
    case 11: //  (* Im3m *)
        RadioButtonSysCubic=true;
        //RadioButtonSysCubicClick(Sender);
        ComboBoxCubic=7;
        break;
    default:
        //return; // TODO Fehlermeldung
        break;
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

    //(****** 1D systems *****)
    if ( RadioButtonSys1D )
    {
        //(* lamellar structure *)
        //if ( ComboBox1D == 0 ) // (* P1 (1) *) - geht nicht anders (für die Zukunft?)
        //{
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c14=0;
            //EditSG.Text:='LAM'; end;
        //}
    } // if RadioButtonSys1D


    //(****** 2D systems *****)
    if ( RadioButtonSys2D )
    {
        switch ( ComboBox2D )
        {
        //(* parallelogram (monoclinic) *)
        case  0: // (* P1 (1) *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c14=0;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=true;
            //EditSG.Text:='P1(1)'; end;
            break;
        case  1: // (* P2 (3)*)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2; c14=0;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=true;
            //EditSG.Text:='P112(3);unique c-axis'; end;
            break;
        //(* rectangular (monoclinic, b-axis unique *)
        case  2: // (* Pm (6) *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1; c14=0;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='P1m1(6);unique b-axis'; end;
            break;
        case  3: // (* Pg (7) h0l:h, h00:h *)
            c1=-1; c2=-1; c3=1; c4=-1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1; c14=1;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='P1a1(7);unique b-axis'; end;
            break;
        case  4: // (* Cm (8) hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k *)
            c1=1; c2=1; c3=1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1; c14=2;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='C1m1(8);unique b-axis'; end;
            break;
        //(* rectangular (orthorhombic) *)
        case  5: // (* Pmm2 (25) *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c14=0;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='Pmm2(25)'; end;
            break;
        case  6: // (* Pmg2 (28) 0kl:k, 0k0:k *)
            c1=-1; c2=1; c3=-1; c4=-1; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c14=3;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='Pbm2(28)'; end;
            break;
        case  7: // (* Pgg2 (32) 0kl:k, h0l:h, h00:h, 0k0:k *)
            c1=-1; c2=1; c3=1; c4=-1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c14=4;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='Pba2(32)'; end;
            break;
        case  8: // (* Cmm2 (35) hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k *)
            c1=1; c2=1; c3=1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c14=5;
            //TrackBarCellB.Visible:=true;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='Cmm2(35)'; end;
            break;
        //(* trigonal *)
        case  9: // (* P3-- (143,156,157) *)
            /*c0=0;*/ c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=1; c14=0;
            //TrackBarCellB.Visible:=false;
            //EditCellB.Text:=EditCellA.Text;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='120';
            //EditSG.Text:='P3(143),P3m1(156),P31m(157)'; trigP:='h'; end;
            break;
        //(* tetragonal *)
        case 10: // (* P4-- (75,99) *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0; c14=0;
            //TrackBarCellB.Visible:=false;
            //EditCellB.Text:=EditCellA.Text;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='P4(75),P4mm(99)'; end;
            break;
        case 11: // (* P4gm (100) 0kl:k, h0l:h, 0k0:k *)
            c1=-1; c2=1; c3=1; c4=-1; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0; c14=6;
            //TrackBarCellB.Visible:=false;
            //EditCellB.Text:=EditCellA.Text;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='90';
            //EditSG.Text:='P4bm(100)'; end;
            break;
        //(* hexagonal *)
        case 12: // (* P6-- (168,183) *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=1; c14=0;
            //TrackBarCellB.Visible:=false;
            //EditCellB.Text:=EditCellA.Text;
            //TrackBarCellGamma.Visible:=false;
            //EditCellGamma.Text:='120';
            //EditSG.Text:='P6(168),P6mm(183)'; end;
            break;
        } // switch ComboBox2D
    } // if RadioButtonSys2D


    //(**************** Cubic crystal system ***************)
    if ( RadioButtonSysCubic )
    {
        switch ( ComboBoxCubic )
        {
        case  0: //(* P--- | *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=2;
            //EditSG.Text:='P23(195),Pm-3(200),P432(207),P-43m(215),Pm-3m(221)'; end;
            break;
        case  1: // (* P21--, P42-- | 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0;  c8=-1; c9=-1; c10=-1; c11=-1; c12=2;
            //EditSG.Text:='P213(198),P4232(208)'; end;
            break;
        case  2: // (* P41-- | 00l:4n *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=1;  c8=-1; c9=-1; c10=-1; c11=-1; c12=2;
            //EditSG.Text:='P4132(213),P4332(212)'; end;
            break;
        case  3: // (* P--n | hhl:l, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0;  c8=-1; c9=-1; c10=0;  c11=-1; c12=2;
            //EditSG.Text:='P-43n(218),Pm-3n(223)'; end;
            break;
        case  4: // (* Pa-- | cyclic permutations: 0kl:k, h0l:l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=5;  c3=-1;  c4=-1;  c5=0;  c6=0;  c7=0;  c8=-1; c9=-1; c10=-1; c11=-1; c12=2;
            //EditSG.Text:='Pa-3(205)'; end;
            break;
        case  5: // (* Pn-- | 0kl:k+l, 00l:l *)
            c1=-1; c2=0;  c3=-1; c4=-1; c5=-1; c6=-1; c7=0;  c8=-1; c9=-1; c10=-1; c11=-1; c12=2;
            //EditSG.Text:='Pn-3(201),Pn-3m(224)'; end;
            break;
        case  6: // (* Pn-n | 0kl:k+l, hhl:l, 00l:l *)
            c1=-1; c2=0;  c3=-1; c4=-1; c5=-1; c6=-1; c7=0;  c8=-1; c9=-1; c10=0;  c11=-1; c12=2;
            //EditSG.Text:='Pn-3n(222)'; end;
            break;
        case  7: // (* I--- | hkl:h+k+l, 0kl:k+l, hhl:l, 00l:l *)
            c1=0; c2=0;  c3=-1; c4=-1; c5=-1; c6=-1; c7=0;  c8=-1; c9=-1; c10=0;  c11=-1; c12=2;
            //EditSG.Text:='I23(197),I213(199),Im-3(204),I432(211),I-43m(217),Im-3m(229)'; end;
            break;
        case  8: // (* I41-- | hkl:h+k+l, 0kl:k+l, hhl:l, 00l:4n *)
            c1=0;  c2=0; c3=-1; c4=-1; c5=-1; c6=-1; c7=1;  c8=-1; c9=-1; c10=0;  c11=-1; c12=2;
            //EditSG.Text:='I4132(214)'; end;
            break;
        case  9: // (* I--d | hkl:h+k+l, 0kl:k+l, hhl:2h+l=4n,l, 00l:4n *)
            c1=0; c2=0;  c3=-1; c4=-1; c5=-1; c6=-1; c7=1;  c8=-1; c9=-1; c10=4;  c11=-1; c12=2;
            //EditSG.Text:='I-43d(220)'; end;
            break;
        case 10: // (* Ia-- | hkl:h+k+l, 0kl:k,l, hhl:l, 00l:l *)
            c1=0;  c2=4; c3=-1; c4=-1; c5=-1; c6=-1; c7=0;  c8=-1; c9=-1; c10=0;  c11=-1; c12=2;
            //EditSG.Text:='Ia-3(206)'; end;
            break;
        case 11: // (* Ia-d | hkl:h+k+l, 0kl:k,l, hhl:2h+l=4n,l, 00l:4n *)
            c1=0;  c2=4; c3=-1; c4=-1; c5=-1; c6=-1; c7=1; c8=-1; c9=-1;  c10=4;  c11=-1; c12=2;
            //EditSG.Text:='Ia-3d(230)'; end;
            break;
        //(* F___ *)    //{NV}-30924
        case 12: // (* F--- | hkl:h+k,h+l,k+l, 0kl:k,l, hhl:h+l, 00l:l *)
            c1=6;  c2=4; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=2;  c11=-1; c12=2;
            //EditSG.Text:='F23(196),Fm-3(202),F432(209),F-43m(216),Fm-3m(225)'; end;
            break;
        case 13: // (* F41-- | hkl:h+k,h+l,k+l, 0kl: k,l, hhl:h+l, 00l:4n *)
            c1=6; c2=4; c3=-1; c4=-1; c5=-1; c6=-1; c7=1;  c8=-1; c9=-1; c10=2; c11=-1; c12=2;
            //EditSG.Text:='F4132(210)'; end;
            break;
        case 14: // (* F--c | hkl:h+k,h+l,k+l, 0kl:k,l, hhl:h,l, 00l:l *)
            c1=6; c2=4; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=5; c11=-1; c12=2;
            //EditSG.Text:='F-43c(219),Fm-3c(226)'; end;
            break;
        case 15: // (* Fd-- | hkl:h+k,h+l,k+l, 0kl:k+l=4n,k,l, hhl:h+l, 00l:4n *)
            c1=6; c2=3; c3=-1; c4=-1; c5=-1; c6=-1; c7=1; c8=-1; c9=-1; c10=2; c11=-1; c12=2;
            //EditSG.Text:='Fd-3(203),Fd-3m(227)'; end;
            break;
        case 16: // (* Fd-c | hkl:h+k,h+l,k+l, 0kl:k+l=4n,k,l, hhl:h,l, 00l:4n *)
            c1=6; c2=3; c3=-1; c4=-1; c5=-1; c6=-1; c7=1; c8=-1; c9=-1; c10=5; c11=-1; c12=2;
            //EditSG.Text:='Fd-3c(228)'; end;
            break;
        } // switch ComboBoxCubic
    } // if RadioButtonSysCubic


    //(**************** Hexagonal crystal system ***************)
    if ( RadioButtonSysHex )
    {
        hmin=-hmax;
        kmin=-kmax;
        lmin=-lmax;
        switch ( ComboBoxHexagonal )
        {
        //(* P___ *)
        case  0: // (* P--- | *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P6(168),P-6(174),P6/m(175),P622(177),P6mm(183),P-62m(189),P-6m2(187),P6/mmm(191)'; end;
            break;
        case  1: // (* P63-- | 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P63(173),P63/m(176),P6322(182)'; end;
            break;
        case  2: // (* P62-- | 00l:3n *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=2; c8=-1; c9=-1; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P62(171),P64(172),P6222(180),P6422(181)'; end;
            break;
        case  3: // (* P61-- | 00l:6n *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=3; c8=-1; c9=-1; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P61(169),P63(170),P6122(178),P6522(179)'; end;
            break;
        case  4: // (* P--c | hhl:l, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=1;
            //EditSG.Text:='P63mc(186),P-62c(190),P63/mmc(194)'; end;
            break;
        case  5: // (* P-c- | h-hl:l, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=0; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P63cm(185),P-6c2(188),P63/mcm(193)'; end;
            break;
        case  6: // (* P-cc | h-hl:l, hhl:l, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=0; c10=0; c11=-1; c12=1;
            //EditSG.Text:='P6cc(184),P6/mcc(192)'; end;
            break;
        } // switch ComboBoxHexagonal
    } // if RadioButtonSysHex


    //(**************** Trigonal crystal system ***************)
    if ( RadioButtonSysTrigonal )
    {
        hmin=-hmax;
        kmin=-kmax;
        lmin=-lmax;
        switch ( ComboBoxTrigonal )
        {
        //(* hexagonal axes, P___ *)
        case  0: // (* P--- | *)
            /*c0=0;*/ c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P3(143),P-3(147),P321(150),P3m1(156),P-3m1(164),P312(149),P31m(157),P-31m(162)'; trigP:='h'; end;
            break;
        case  1: //  (* P31-- | 00l:3n *)
            /*c0=0;*/ c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=2; c8=-1; c9=-1; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P31(144),P32(145),P3121(152),P3221(154),P3112(151),P3212(153)'; trigP:='h'; end;
            break;
        case  2: //  (* P--c | hhl:l, 00l:l *)
            /*c0=0;*/ c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=1;
            //EditSG.Text:='P31c(159),P-31c(163)'; trigP:='h'; end;
            break;
        case  3: // (* P-c- | h-hl:l, 00l:l *)
            /*c0=0;*/ c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=0; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='P3c1(158),P-3c1(165)'; trigP:='h'; end;
            break;
        //(* hexagonal axes, R___ *)
        case  4: // (* R(obv)-- | hkl:-h+k+l=3n, h-hl:h+l=3n, hhl:l=3n, 00l:3n *)
            /*c0=1;*/ c1=4; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=2; c8=-1; c9=1; c10=1; c11=-1; c12=1;
            //EditSG.Text:='R3(146),R-3(148),R32(155),R3m(160),R-3m(166)'; trigP:='h'; end;
            break;
        case  5: // (* R(obv)-c | hkl:-h+k+l=3n, h-hl:h+l=3n,l, hhl:l=3n, 00l:6n *)
            /*c0=1;*/ c1=4; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=3; c8=-1; c9=3; c10=1; c11=-1; c12=1;
            //EditSG.Text:='R3c(161),R-3c(148)'; trigP:='h'; end;
            break;
        case  6: // (* R(rev)-- | hkl:h-k+l=3n, h-hl:-h+l=3n, hhl:l=3n, 00l:3n *)
            /*c0=2;*/ c1=5; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=2; c8=-1; c9=2; c10=1; c11=-1; c12=1;
            //EditSG.Text:='R3(146),R-3(148),R32(155),R3m(160),R-3m(166)'; trigP:='h'; end;
            break;
        case  7: // (* R(rev)-c | hkl:h-k+l=3n, h-hl:-h+l=3n,l, hhl:l=3n, 00l:6n *)
            /*c0=2;*/ c1=5; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=3; c8=-1; c9=3; c10=1; c11=0; c12=1;
            //EditSG.Text:='R3c(161),R-3c(167)'; trigP:='h'; end;
            break;
        //(* rhombohedral axes, R___ *)
        case  8: // (* R-- |  *)
            /*c0=-1;*/ c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=1;
            //EditSG.Text:='R3(146),R-3(148),R32(155),R3m(160),R-3m(166); rhombohedral'; trigP:='r'; end;
            trigP='r';
            break;
        case  9: // (* R-c | hhl:l, hhh:h *)
            /*c0=-1;*/ c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=0; c11=0; c12=1;
            //EditSG.Text:='R3c(161),R-3c(167); rhombohedral'; trigP:='r'; end;
            trigP='r';
            break;
        } // switch ComboBoxTrigonal
    } // if RadioButtonSysTrigonal


    //(**************** Tetragonal crystal system ***************)
    if ( RadioButtonSysTetragonal )
    {
        //hmin=0;
        //kmin=0;
        //lmin=0;
        switch ( ComboBoxTetragonal )
        {
        case  0: // (* P--- | *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P4(75),P-4(81),P4/m(83),P422(89),P4mm(99),P-42m(111),P-4m2(115),P4/mmm(123)'; end;
            break;
        case  1: // (* P-21- | 0k0:k &permu: h00:h *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P4212(90),P-421m(113)'; end;
            break;
        case  2: // (* P42-- | 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P42(77),P42/m(84),P4222(93)'; end;
            break;
        case  3: // (* P4221- | 00l:l, 0k0:k &permu: h00:h *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P42212(94)'; end;
            break;
        case  4: // (* P41-- | 00l:4n *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P41(76),P43(78),P4122(91),P4322(95)'; end;
            break;
        case  5: // (* P4121- | 00l:4n, 0k0:k &permu: h00:h *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=0; c7=1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P41212(92),P43212(96)'; end;
            break;
        case  6: // (* P--c | hhl:l, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P42mc(105),P-42c(112),P42/mmc(131)'; end;
            break;
        case  7: // (* P-21c | hhl:l, 00l:l, 0k0:k &permu: h00:h *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P-421c(114)'; end;
            break;
        case  8: // (* P-b- | 0kl:k, 0k0:k &permu: h0l:h *)
            c1=-1; c2=1; c3=1; c4=-1; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P4bm(100),P-4b2(117),P4/mbm(127)'; end;
            break;
        case  9: // (* P-bc | 0kl:k, hhl:l, 00l:l, 0k0:k &permu: h0l:h *)
            c1=-1; c2=1; c3=1; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P42bc(106),P42/mbc(135)'; end;
            break;
        case 10: // (* P-c- | 0kl:l, 00l:l &permu: h0l:l*)
            c1=-1; c2=2; c3=2; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P42cm(101),P-4c2(116),P42mcm(132)'; end;
            break;
        case 11: // (* P-cc | 0kl:l, hhl:l, 00l:l &permu: h0l:l *)
            c1=-1; c2=2; c3=2; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P4cc(103),P4/mcc(124)'; end;
            break;
        case 12: // (* P-n- | 0kl:k+l, 00l:l, 0k0:k &permu: h0l:h+l *)
            c1=-1; c2=0; c3=0; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P42nm(102),P-4n2(118),P42/mnm(136)'; end;
            break;
        case 13: // (* P-nc | 0kl:k+l, hhl:l, 00l:l, 0k0:k &permu: h0l:h+l *)
            c1=-1; c2=0; c3=0; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P4nc(104),P4/mnc(128)'; end;
            break;
        case 14: // (* Pn-- | hk0:h+k, 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=0; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P4/n(85),P4/nmm(129)'; end;
            break;
        case 15: // (* P42/n-- | hk0:h+k, 00l:l, 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P42/n(86)'; end;
            break;
        case 16: // (* Pn-c | hk0:h+k, hhl:l, 00l:l, 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P42/nmc(137)'; end;
            break;
        case 17: // (* Pnb- | hk0:h+k, 0kl:k, 0k0:k &permu: h0l:h *)
            c1=-1; c2=1; c3=1; c4=0; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P4/nbm(125)'; end;
            break;
        case 18: // (* Pnbc | hk0:h+k, 0kl:k, hhl:l, 00l:l, 0k0:k &permu: h0l:h*)
            c1=-1; c2=1; c3=1; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P42/nbc(133)'; end;
            break;
        case 19: // (* Pnc- | hk0:h+k, 0kl:l, 00l:l, 0k0:k &permu: h0l:l *)
            c1=-1; c2=2; c3=2; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P42/ncm(138)'; end;
            break;
        case 20: // (* Pncc | hk0:h+k, 0kl:l, hhl:l, 00l:l, 0k0:k &permu: h0l:l *)
            c1=-1; c2=2; c3=2; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P4/ncc(130)'; end;
            break;
        case 21: // (* Pnn- | hk0:h+k, 0kl:k+l, 00l:l, 0k0:k *)
            c1=-1; c2=0; c3=-1; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=0;
            //EditSG.Text='P42/nnm(134)'; end;
            break;
        case 22: // (* Pnnc | hk0:h+k, 0kl:k+l, hhl:l, 00l:l, 0k0:k *)
            c1=-1; c2=0; c3=-1; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='P4/ncc(126)'; end;
            break;

        case 23: // (* I--- | hkl:h+k+l, hk0:h+k, 0kl:k+l, hhl:l, 00l:l, 0k0:k *)
            c1=0; c2=0; c3=-1; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='I4(79),I-4(82),I4/m(87),I422(97),I4mm(107),I-42m(121),I-4m2(119),I4/mmm(13)'; end;
            break;
        case 24: // (* I41-- | hkl:h+k+l, hk0:h+k, 0kl:k+l, hhl:l, 00l:4n, 0k0:k *)
            c1=0; c2=0; c3=-1; c4=0; c5=-1; c6=0; c7=1; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='I41(80),I4122(98)'; end;
            break;
        case 25: // (* I--d | hkl:h+k+l, hk0:h+k, 0kl:k+l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h *)
            c1=0; c2=0; c3=-1; c4=0; c5=-1; c6=0; c7=1; c8=0; c9=-1; c10=4; c11=-1; c12=0;
            //EditSG.Text='I41md(109),I-42d(122)'; end;
            break;
        case 26: // (* I-c- | hkl:h+k+l, hk0:h+k, 0kl:k,l, hhl:l, 00l:l, 0k0:k *)
            c1=0; c2=4; c3=-1; c4=0; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='I4cm(108),I-4c2(120),I4/mcm(140)'; end;
            break;
        case 27: // (* I-cd | hkl:h+k+l, hk0:h+k, 0kl:k,l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h *)
            c1=0; c2=4; c3=-1; c4=0; c5=-1; c6=0; c7=1; c8=0; c9=-1; c10=4; c11=-1; c12=0;
            //EditSG.Text='I41cd(110)'; end;
            break;
        case 28: // (* I41/a-- | hkl:h+k+l, hk0:h,k, 0kl:k+l, hhl:l, 00l:4n, 0k0:k *)
            c1=0; c2=0; c3=-1; c4=4; c5=-1; c6=0; c7=1; c8=-1; c9=-1; c10=0; c11=-1; c12=0;
            //EditSG.Text='I41/a(88)'; end;
            break;
        case 29: // (* Ia-d | hkl:h+k+l, hk0:h,k, 0kl:k+l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h *)
            c1=0; c2=0; c3=-1; c4=4; c5=-1; c6=0; c7=1; c8=0; c9=-1; c10=4; c11=-1; c12=0;
            //EditSG.Text='I41/amd(141)'; end;
            break;
        case 30: // (* Iacd | hkl:h+k+l, hk0:h,k, 0kl:k,l, hhl:2h+l=4n,l, 00l:4n, 0k0:k, hh0:h *)
            c1=0; c2=4; c3=-1; c4=4; c5=-1; c6=0; c7=1; c8=0; c9=-1; c10=4; c11=-1; c12=0;
            //EditSG.Text='I41/acd(142)'; end;
            break;
        } // switch ( ComboBoxTetragonal )
    } // if ( RadioButtonSysTetragonal )


    //(**************** Orthorhombic crystal system ***************)
    if ( RadioButtonSysOrtho )
    {
        //hmin=0;
        //kmin=0;
        //lmin=0;
        switch ( ComboBoxOrthorhombic ) // 0 .. 110
        {
        //(* primitive *)
        case   0: // (* P--- | *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P222(16),Pmm2(25),Pm2m(25),P2mm(25),Pmmm(47)'; end;
            break;
        case   1: // (* P--21 | 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2221(17)'; end;
            break;
        case   2: // (* P-21- | 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2212(17)'; end;
            break;
        case   3: // (* P-2121 | 0k0:k, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P22121(18)'; end;
            break;
        case   4: // (* P21-- | h00:h *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2122(17)'; end;
            break;
        case   5: // (* P21-21 | h00:h, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P21221(18)'; end;
            break;
        case   6: // (* P2121- | h00:h, 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P21212(18)'; end;
            break;
        case   7: // (* P212121 | h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P212121(19)'; end;
            break;
        case   8: // (* P--a | hk0:h, h00:h *)
            c1=-1; c2=-1; c3=-1; c4=1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pm2a(28),P21ma(26),Pmma(51)'; end;
            break;
        case   9: // (* P--b | hk0:k, 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=2; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pm21b(26),P2mb(28),Pmmb(51)'; end;
            break;
        case  10: // (* P--n | hk0:h+k, h00:h, 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pm21n(31),P21mn(31),Pmmn(59)'; end;
            break;
        case  11: // (* P-a- | h0l:h, h00:h *)
            c1=-1; c2=-1; c3=1; c4=-1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pma2(28),P21am(26),Pmam(51)'; end;
            break;
        case  12: // (* P-aa | h0l:h, hk0:h, h00:h *)
            c1=-1; c2=-1; c3=1; c4=1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2aa(27),Pmaa(49)'; end;
            break;
        case  13: // (* P-ab | h0l:h, hk0:k, h00:h, 0k0:k *)
            c1=-1; c2=-1; c3=1; c4=1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P21ab(29),Pmab(57)'; end;
            break;
        case  14: // (* P-an | h0l:h, hk0:h+k, h00:h, 0k0:k *)
            c1=-1; c2=-1; c3=1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2an(30),Pman(53)'; end;
            break;
        case  15: // (* P-c- | h0l:l, 00l:l *)
            c1=-1; c2=-1; c3=2; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pmc21(26),P2cm(28),Pmcm(51)'; end;
            break;
        case  16: // (* P-ca | h0l:l, hk0:h, h00:h, 00l:l *)
            c1=-1; c2=-1; c3=2; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P21ca(29),Pmca(57)'; end;
            break;
        case  17: // (* P-cb | h0l:l, hk0:k, 0k0:k, 00l:l *)
            c1=-1; c2=-1; c3=2; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2cb(32),Pmcb(55)'; end;
            break;
        case  18: // (* P-cn | h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=-1; c3=2; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P21cn(33),Pmcn(62)'; end;
            break;
        case  19: // (* P-n- | h0l:h+l, h00:h, 00l:l *)
            c1=-1; c2=-1; c3=0; c4=-1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pmn21(31),P21nm(31),Pmnm(59)'; end;
            break;
        case 20: // (* P-na | h0l:h+l, hk0:h, h00:h, 00l:l *)
            c1=-1; c2=-1; c3=0; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2na(30),Pmna(53)'; end;
            break;
        case 21: // (* P-nb | h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=-1; c3=0; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P21nb(33),Pmnb(62)'; end;
            break;
        case 22: // (* P-nn | h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P2nn(34),Pmnn(58)'; end;
            break;
        case 23: // (* Pb-- | 0kl:k, 0k0:k *)
            c1=-1; c2=1; c3=-1; c4=-1; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbm2(28),Pb21m(26),Pbmm(51)'; end;
            break;
        case 24: // (* Pb-a | 0kl:k, hk0:h, h00:h, 0k0:k *)
            c1=-1; c2=1; c3=-1; c4=1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pb21a(29),Pbma(57)'; end;
            break;
        case 25: // (* Pb-b | 0kl:k, hk0:k, 0k0:k *)
            c1=-1; c2=1; c3=-1; c4=2; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pb2b(27),Pbmb(49)'; end;
            break;
        case 26: // (* Pb-n | 0kl:k, hk0:h+k, h00:h, 0k0:k *)
            c1=-1; c2=1; c3=-1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pb2n(30),Pbmn(53)'; end;
            break;
        case 27: // (* Pba- | 0kl:k, h0l:h, h00:h, 0k0:k *)
            c1=-1; c2=1; c3=1; c4=-1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pba2(32),Pbam(55)'; end;
            break;
        case 28: // (* Pbaa | 0kl:k, h0l:h, hk0:h, h00:h, 0k0:k *)
            c1=-1; c2=1; c3=1; c4=1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbaa(54)'; end;
            break;
        case 29: // (* Pbab | 0kl:k, h0l:h, hk0:k, h00:h, 0k0:k *)
            c1=-1; c2=1; c3=1; c4=2; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbab(54)'; end;
            break;
        case 30: // (* Pban | 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k *)
            c1=-1; c2=1; c3=1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pban(50)'; end;
            break;
        case 31: // (* Pbc- | 0kl:k, h0l:l, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=2; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbc21(29),Pbcm(57)'; end;
            break;
        case 32: // (* Pbca | 0kl:k, h0l:l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=2; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbca(61)'; end;
            break;
        case 33: // (* Pbcb | 0kl:k, h0l:l, hk0:k, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=2; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbcb(54)'; end;
            break;
        case 34: // (* Pbcn | 0kl:k, h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=2; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbcn(60)'; end;
            break;
        case 35: // (* Pbn- | 0kl:k, h0l:h+l, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=0; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbn21(33),Pbnm(62)'; end;
            break;
        case 36: // (* Pbna | 0kl:k, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=0; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbna(60)'; end;
            break;
        case 37: // (* Pbnb | 0kl:k, h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=0; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbnb(56)'; end;
            break;
        case 38: // (* Pbnn | 0kl:k, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=1; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pbnn(52)'; end;
            break;
        case 39: // (* Pc-- | 0kl:l, 00l:l *)
            c1=-1; c2=2; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcm21(26),Pc2m(28),Pcmm(51)'; end;
            break;
        case 40: // (* Pc-a | 0kl:l, hk0:h, h00:h, 00l:l *)
            c1=-1; c2=2; c3=-1; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pc2a(32),Pcma(55)'; end;
            break;
        case 41: // (* Pc-b | 0kl:l, hk0:k, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=-1; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pc21b(29),Pcmb(57)'; end;
            break;
        case 42: // (* Pc-n | 0kl:l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=-1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pc21n(33),Pcmn(62)'; end;
            break;
        case 43: // (* Pca- | 0kl:l, h0l:h, h00:h, 00l:l *)
            c1=-1; c2=2; c3=1; c4=-1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pca21(29),Pcam(57)'; end;
            break;
        case 44: // (* Pcaa | 0kl:l, h0l:h, hk0:h, h00:h, 00l:l *)
            c1=-1; c2=2; c3=1; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcaa(54)'; end;
            break;
        case 45: // (* Pcab | 0kl:l, h0l:h, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=1; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcab(61)'; end;
            break;
        case 46: // (* Pcan | 0kl:l, h0l:h, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcan(60)'; end;
            break;
        case 47: // (* Pcc- | 0kl:l, h0l:l, 00l:l *)
            c1=-1; c2=2; c3=2; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcc2(27),Pccm(49)'; end;
            break;
        case 48: // (* Pcca | 0kl:l, h0l:l, hk0:h, h00:h, 00l:l *)
            c1=-1; c2=2; c3=2; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcca(54)'; end;
            break;
        case 49: // (* Pccb | 0kl:l, h0l:l, hk0:k, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=2; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pccb(54)'; end;
            break;
        case 50: // (* Pccn | 0kl:l, h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=2; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pccn(56)'; end;
            break;
        case 51: // (* Pcn- | 0kl:l, h0l:h+l, h00:h, 00l:l *)
            c1=-1; c2=2; c3=0; c4=-1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcn2(30),Pcnm(53)'; end;
            break;
        case 52: // (* Pcna | 0kl:l, h0l:h+l, hk0:h, h00:h, 00l:l *)
            c1=-1; c2=2; c3=0; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcna(50)'; end;
            break;
        case 53: // (* Pcnb | 0kl:l, h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=0; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcnb(60)'; end;
            break;
        case 54: // (* Pcnn | 0kl:l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=2; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pcnn(52)'; end;
            break;
        case 55: // (* Pn-- | 0kl:k+l, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=-1; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnm21(31),Pn21m(31),Pnmm(59)'; end;
            break;
        case 56: // (* Pn-a | 0kl:k+l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=-1; c4=1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pn21a(33),Pnma(62)'; end;
            break;
        case 57: // (* Pn-b | 0kl:k+l, hk0:k, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=-1; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pn2b(30),Pnmb(53)'; end;
            break;
        case 58: // (* Pn-n | 0kl:k+l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=-1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pn2n(34),Pnmn(58)'; end;
            break;
        case 59: // (* Pna- | 0kl:k+l, h0l:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=1; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pna21(33),Pnam(62)'; end;
            break;
        case 60: // (* Pnaa | 0kl:k+l, h0l:h, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=1; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnaa(56)'; end;
            break;
        case 61: // (* Pnab | 0kl:k+l, h0l:h, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=1; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnab(60)'; end;
            break;
        case 62: // (* Pnan | 0kl:k+l, h0l:h, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnan(52)'; end;
            break;
        case 63: // (* Pnc- | 0kl:k+l, h0l:l, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=2; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnc2(30),Pncm(53)'; end;
            break;
        case 64: // (* Pnca | 0kl:k+l, h0l:l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=2; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnca(60)'; end;
            break;
        case 65: // (* Pncb | 0kl:k+l, h0l:l, hk0:k, 0k0:k, 00l:l |  *)
            c1=-1; c2=0; c3=2; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pncb(50)'; end;
            break;
        case 66: // (* Pncn | 0kl:k+l, h0l:l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=2; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pncn(52)'; end;
            break;
        case 67: // (* Pnn- | 0kl:k+l, h0l:h+l, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=0; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnn2(34),Pnnm(58)'; end;
            break;
        case 68: // (* Pnna | 0kl:k+l, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=0; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnna(52)'; end;
            break;
        case 69: // (* Pnnb | 0kl:k+l, h0l:h+l, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=0; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnnb(52)'; end;
            break;
        case 70: // (* Pnnn | 0kl:k+l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=-1; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Pnnn(48)'; end;
            break;

        //(* centered C___ *)
        case 71: // (* C--- | hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k *)
            c1=1; c2=1; c3=1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='C222(21),Cmm2(35),Cm2m(38),C2mm(38),Cmmm(65)'; end;
            break;
        case 72: // (* C--21 | hkl:h+k, 0kl:k, h0l:h, hk0:h+k, h00:h, 0k0:k, 00l:. *)
            c1=1; c2=1; c3=1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='C2221(20)'; end;
            break;
        case 73: // (* C--(ab) | hkl:h+k, 0kl:k, h0l:h, hk0:h,k, h00:h, 0k0:k *)
            c1=1; c2=1; c3=1; c4=4; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Cm2e(39),C2me(39),Cmme(67)'; end;
            break;
        case 74: // (* C-c- | hkl:h+k, 0kl:k, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=1; c2=1; c3=4; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ccm21(36),C2cm(40),Cmcm(63)'; end;
            break;
        case 75: // (* C-c(ab) | hkl:h+k, 0kl:k, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=1; c2=1; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='C2ce(41),Cmce(64)'; end;
            break;
        case 76: // (* Cc-- | hkl:h+k, 0kl:k,l, h0l:h, hk0:h+k, h00:h, 0k0:k, 001:l *)
            c1=1; c2=4; c3=1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ccm21(36),Cc2m(40),Ccmm(63)'; end;
            break;
        case 77: // (* Cc-(ab) | hkl:h+k, 0kl:k,l, h0l:h, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=1; c2=4; c3=1; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Cc2e(41),Ccme(64)'; end;
            break;
        case 78: // (* Ccc- | hkl:h+k, 0kl:k,l, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=1; c2=4; c3=4; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ccc2(37),Cccm(66)'; end;
            break;
        case 79: // (* Ccc(ab) | hkl:h+k, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=1; c2=4; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ccce(68)'; end;
            break;

        //(* centered B___ *)
        case 80: // (* B--- | hkl:h+l, 0kl:l, h0l:h+l, hk0:h, h00:h, 00l:l *)
            c1=3; c2=2; c3=0; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='B222(21),Bmm2(38),Bm2m(35),B2mm(38),Bmmm(65)'; end;
            break;
        case 81: // (* B-21- | hkl:h+l, 0kl:l, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=3; c2=2; c3=0; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='B2212(20)'; end;
            break;
        case 82: // (* B--b | hkl:h+l, 0kl:l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=3; c2=2; c3=0; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Bm21b(36),B2mb(40),Bmmb(63)'; end;
            break;
        case 83: // (* B-(ac)- | hkl:h+l, 0kl:l, h0l;h,l, hk0:h, h00:h, 00l:l *)
            c1=3; c2=2; c3=4; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Bme2(39),B2em(39),Bmem(67)'; end;
            break;
        case 84: // (* B-(ac)b | hkl:h+l, 0kl:l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=3; c2=2; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='B2eb(41),Bmeb(64)'; end;
            break;
        case 85: // (* Bb-- | hkl:h+l, 0kl:k,l, h0l:h+l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=3; c2=4; c3=0; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Bbm2(40),Bb21m(36),Bbmm(63)'; end;
            break;
        case 86: // (* Bb-b | hkl:h+l, 0kl:k,l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=3; c2=4; c3=0; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Bb2b(37),Bbmb(66)'; end;
            break;
        case 87: // (* Bb(ac)- | hkl:h+l, 0kl:k,l, h0l:h,l, hk0:h, h00:h, 0k0:k, 00l:l *)
            c1=3; c2=4; c3=4; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Bbe2(41),Bbem(64)'; end;
            break;
        case 88: // (* Bb(ac)b | hkl:h+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=3; c2=4; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Bbeb(68)'; end;
            break;

        //(* centered A___ *)
        case 89: // (* A--- | hkl:k+l, 0kl:k+l, h0l:l, hk0:k, 0k0:k, 00l:l *)
            c1=2; c2=0; c3=2; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='A222(21),Amm2(38),Am2m(38),A2mm(35),Ammm(65)'; end;
            break;
        case 90: // (* A21-- | hkl:k+l, 0kl:k+l, h0l:l, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=2; c2=0; c3=2; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='A2122(20)'; end;
            break;
        case 91: // (* A--a | hkl:k+l, 0kl:k+l, h0l:l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=2; c2=0; c3=2; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Am2a(40),A21ma(36),Amma(63)'; end;
            break;
        case 92: // (* A-a-| hkl:k+l, 0kl:k+l, h0l:h,l, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=2; c2=0; c3=4; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ama2(40),A21am(36),Amam(63)'; end;
            break;
        case 93: // (* A-aa | hkl:k+l, 0kl:k+l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=2; c2=0; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='A2aa(37),Amaa(66)'; end;
            break;
        case 94: // (* A(bc)-- | hkl:k+l, 0kl:k,l, h0l:l, hk0:k, 0k0:k, 00l:l *)
            c1=2; c2=4; c3=2; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Aem2(39),Ae2m(39),Aemm(67)'; end;
            break;
        case 95: // (* A(bc)-a | hkl:k+l, 0kl:k,l, h0l:l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=2; c2=4; c3=2; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ae2a(41),Aema(64)'; end;
            break;
        case 96: // (* A(bc)a- | hkl:k+l, 0kl:k,l, h0l:h,l, hk0:k, h00:h, 0k0:k, 00l:l *)
            c1=2; c2=4; c3=4; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Aea2(41),Aeam(64)'; end;
            break;
        case 97: // (* A(bc)aa | hkl:k+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=2; c2=4; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Aeaa(68)'; end;
            break;

        //(* centered I___ *)
        case 98: // (* I--- | hkl:h+k+l, 0kl:k+l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='I222(23),I212121(24),Imm2(44),Im2m(44),I2mm(44),Immm(71)'; end;
            break;
        case 99: // (* I--(ab) | hkl:h+k+l, 0kl:k+l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=0; c3=0; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Im2a(46),I2mb(46),Imma(74),Immb(74)'; end;
            break;
        case 100: // (* I-(ac)- | hkl:h+k+l, 0kl:k+l, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=0; c3=4; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ima2(46),I2cm(46),Imam(74),Imcm(74)'; end;
            break;
        case 101: // (* I-cb | hkl:h+k+l, 0kl:k+l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=0; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='I2cb(45),Imcb(72)'; end;
            break;
        case 102: // (* I(bc)-- | hkl:h+k+l, 0kl:k,l, h0l:h+l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=4; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Iem2(46),Ie2m(46),Iemm(74)'; end;
            break;
        case 103: // (* Ic-a | hkl:h+k+l, 0kl:k,l, h0l:h+l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=4; c3=0; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ic2a(45),Icma(72)'; end;
            break;
        case 104: // (* Iba- | hkl:h+k+l, 0kl:k,l, h0l:h,l, hk0:h+k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=4; c3=4; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Iba2(45),Ibam(72)'; end;
            break;
        case 105: // (* Ibca | hkl:h+k+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=0; c2=4; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Ibca(73),Icab(73)'; end;
            break;

        //(* centered F___ *)
        case 106: // (* F--- | hkl:h+k,h+l,k+l, 0kl:k,l, h0l:h,l, hk0:h,k, h00:h, 0k0:k, 00l:l *)
            c1=6; c2=4; c3=4; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='F222(22),Fmm2(42),Fm2m(42),F2mm(42),Fmmm(69)'; end;
            break;
        case 107: // (* F-dd | hkl:h+k,h+l,k+l, 0kl:k,l, h0l:4n,h,l, hk0:4n,h,k, h00:4n, 0k0:4n, 00l:4n *)
            c1=6; c2=4; c3=3; c4=3; c5=1; c6=1; c7=1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='F2dd(43)'; end;
            break;
        case 108: // (* Fd-d | hkl:h+k,h+l,k+l, 0kl:4n,k,l, h0l:h,l, hk0:4n,h,k, h00:4n, 0k0:4n, 00l:4n *)
            c1=6; c2=3; c3=4; c4=3; c5=1; c6=1; c7=1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Fd2d(43)'; end;
            break;
        case 109: // (* Fdd- | hkl:h+k,h+l,k+l, 0kl:4n,k,l, h0l:4n,h,l, hk0:h,k, h00:4n, 0k0:4n, 00l:4n *)
            c1=6; c2=3; c3=3; c4=4; c5=1; c6=1; c7=1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Fdd2(43)'; end;
            break;
        case 110: // (* Fddd | hkl:h+k,h+l,k+l, 0kl:4n,k,l, h0l:4n,h,l, hk0:4n,h,k, h00:4n, 0k0:4n, 00l:4n *)
            c1=6; c2=3; c3=3; c4=3; c5=1; c6=1; c7=1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='Fddd(70)'; end;
            break;
        } // switch ComboBoxOrthorhombic
    } // if RadioButtonSysOrtho


    //(**************** Monoclinic crystal system ***************)
    if ( RadioButtonSysMonoclinic )
    {
        //hmin=0;
        //kmin=0;
        //lmin=-lmax;
        switch ( ComboBoxMonoclinic )
        {
        //(* P___ unique axis b *)
        case 0: // (* P1-1 *)  <<< das hier ist auch zu sehen!
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;

            //EditSG.Text:='P121(3),P1m1(6),P12/m1(10);unique b-axis'; end;
            //>>>>>> Das hier ist im Screenshot von Prof. Förster zu lesen .... //{NV}-31484

            break;
        case 1: // (* P1211 | 0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='P1211(4),P21/m1(11);unique b-axis'; end;
            break;
        case 2: // (* P1a1 | h0l,h00:h *)
            c1=-1; c2=-1; c3=1; c4=-1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='P1a1(7),P12/a1(13);unique b-axis'; end;
            break;
        case 3: // (* P121/a1 | h0l,h00:h, 0k0:k *)
            c1=-1; c2=-1; c3=1; c4=-1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='P121/a1(14);unique b-axis'; end;
            break;
        case 4: // (* P1c1 | h0l,00l:l *)
            c1=-1; c2=-1; c3=2; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='P1c1(7),P12/c1(13);unique b-axis'; end;
            break;
        case 5: // (* P121/c1 | h0l,00l:l, 0k0:k *)
            c1=-1; c2=-1; c3=2; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='P121/c1(14);unique b-axis'; end;
            break;
        case 6: // (* P1n1 | h0l,h00,00l:h+l *)
            c1=-1; c2=-1; c3=0; c4=-1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='P1n1(7),P12/n1(13);unique b-axis'; end;
            break;
        case 7: // (* P121/n1 | h0l,h00,00l:h+l, 0k0:k *)
            c1=-1; c2=-1; c3=0; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='P121/n1(14);unique b-axis'; end;
            break;

            //(* C___ *)
        case 8: // (* C1-1 | hkl,0kl,hk0:h+k, h0l,h00,00l:h, 0k0:k *)
            c1=1; c2=1; c3=1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='C121(5),C1m1(8),C12/m1(12);unique b-axis'; end;
            break;
        case 9: // (* C1c1 | hkl,0kl,hk0:h+k, h0l,h00,00l:h,l, 0k0:k *)
            c1=1; c2=1; c3=4; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='C1c1(9),C12/c1(15);unique b-axis'; end;
            break;

            //(* A___ *)
        case 10: // (* A1-1 | hkl,0kl,hk0:k+l, h0l,h00,00l:l, 0k0:k *)
            c1=2; c2=0; c3=2; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='A121(15),A1m1(8),A12/m1(12);unique b-axis'; end;
            break;
        case 11: // (* A1n1 | hkl,0kl,hk0:k+l, h0l,h00,00l:h,l, 0k0:k *)
            c1=2; c2=0; c3=4; c4=2; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='A1n1(9),A12/n1(15);unique b-axis'; end;
            break;

            //(* I___ *)
        case 12: // (* I1-1 | hkl,0kl,hk0:h+k+l, h0l,h00,00l:h+l, 0k0:k *)
            c1=0; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='I121(5),I1m1(8),I12/m1(12)'; end;
            break;
        case 13: //  (* I1a1 | hkl,0kl,hk0:h+k+l, h0l,h00,00l:h,l, 0k0:k *)
            c1=0; c2=0; c3=4; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=1;
            //EditSG.Text:='I1a1(9),I12/a1(15);unique b-axis'; end;
            break;

            //(* P___ unique axis c *)
        case 14: // (* P11- *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P112(3),P11m(6),P112/m(10);unique c-axis'; end;
            break;
        case 15: // (* P1121 | 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P1121(4),P1121/m(11);unique c-axis'; end;
            break;
        case 16: // (* P11a | hk0,h00,0k0:h *)
            c1=-1; c2=-1; c3=-1; c4=1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P11a(7),P112/a(13);unique c-axis'; end;
            break;
        case 17: // (* P1121/a | hk0,h00,0k0:h, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P1121/a(14);unique c-axis'; end;
            break;
        case 18: // (* P11b | hk0,h00,0k0:k *)
            c1=-1; c2=-1; c3=-1; c4=2; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P11b(7),P112/b(13);unique c-axis'; end;
            break;
        case 19: // (* P1121/b | hk0,h00,0k0:k, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P1121/b(14);unique c-axis'; end;
            break;
        case 20: // (* P11n | hk0,h00,0k0:h+k *)
            c1=-1; c2=-1; c3=-1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P11n(7),P112/n(13);unique c-axis'; end;
            break;
        case 21: // (* P1121/n | hk0,h00,0k0:h+k, 00l:l *)
            c1=-1; c2=-1; c3=-1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='P1121/n(14);unique c-axis'; end;
            break;

            //(* B___ *)
        case 22: // (* B11- | hkl,0kl,h0l:h+l, hk0,h00,0k0:h, 00l:l *)
            c1=3; c2=2; c3=0; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;  c13=2;
            //EditSG.Text:='B112(5),B11m(8),B112/m(12);unique c-axis'; end;
            break;
        case 23: // (* B11n | hkl,0kl,h0l:h+l, hk0,h00,0k0:h,k, 00l:l *)
            c1=3; c2=2; c3=0; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;  c13=2;
            //EditSG.Text:='B11n(9),B112/n(15);unique c-axis'; end;
            break;

            //(* A___ *)
        case 24: // (* A11- | hkl,0kl,h0l:k+l, hk0,h00,0k0:k, 00l:l *)
            c1=2; c2=0; c3=2; c4=2; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='A112(15),A11m(8),A112/m(12);unique c-axis'; end;
            break;
        case 25: // (* A11a | hkl,0kl,h0l:k+l, hk0,h00,0k0:h,k, 00l:l *)
            c1=2; c2=0; c3=2; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='A11a(9),A112/a(15);unique c-axis'; end;
            break;

            //(* I___ *)
        case 26: // (* I11- | hkl,0kl,h0l:h+k+l, hk0,h00,0k0:h+k, 00l:l *)
            c1=0; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='I112(5),I11m(8),I112/m(12);unique c-axis'; end;
            break;
        case 27: //  (* I11b | hkl,0kl,h0l:h+k+l, hk0,h00,0k0:h,k, 00l:l *)
            c1=0; c2=0; c3=0; c4=4; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=2;
            //EditSG.Text:='I11b(9),I112/b(15);unique c-axis'; end;
            break;

            //(* P___ unique axis a *)
        case 28: // (* P-11 *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='P211(3),Pm11(6),P2/m11(10);unique a-axis'; end;
            break;
        case 29: // (* P2111 | h00:h *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=0; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='P2111(4),P21/m11(11);unique a-axis'; end;
            break;
        case 30: // (* Pb11 | 0kl,0k0,00l:k *)
            c1=-1; c2=1; c3=-1; c4=-1; c5=-1; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='Pb11(7),P2/b11(13);unique a-axis'; end;
            break;
        case 31: // (* P21/b11 | 0kl,0k0,00l:k, h00:h *)
            c1=-1; c2=1; c3=-1; c4=-1; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='P21/b11(14);unique q-axis'; end;
            break;
        case 32: // (* Pc11 | 0kl,0k0,00l:l *)
            c1=-1; c2=2; c3=-1; c4=-1; c5=-1; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='Pc11(7),P2/c11(13);unique a-axis'; end;
            break;
        case 33: // (* P21/c11 | 0kl,0k0,00l:l, h00:h *)
            c1=-1; c2=2; c3=-1; c4=-1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='P21/c11(14);unique a-axis'; end;
            break;
        case 34: // (* Pn11 | 0kl,0k0,00l:k+l *)
            c1=-1; c2=0; c3=-1; c4=-1; c5=-1; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='Pn11(7),P2/n11(13);unique a-axis'; end;
            break;
        case 35: // (* P21/n11 | 0kl,0k0,00l:k+l, h00:h *)
            c1=-1; c2=0; c3=-1; c4=-1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='P21/n11(14);unique c-axis'; end;
            break;

            //(* C___ *)
        case 36: // (* C-11 | hkl,h0l,hk0:h+k, 0kl,0k0,00l:k, h00:h *)
            c1=1; c2=1; c3=1; c4=0; c5=0; c6=0; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='C211(5),Cm11(8),C2/m11(12);unique a-axis'; end;
            break;
        case 37: // (* Cn11 | hkl,h0l,hk0:h+k, 0kl,0k0,00l:k,l, h00:h *)
            c1=1; c2=4; c3=1; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='Cn11(9),C2/n11(15);unique a-axis'; end;
            break;

            //(* B___ *)
        case 38: // (* B-11 | hkl,h0l,hk0:h+l, 0kl,0k0,00l:l, h00:h *)
            c1=3; c2=2; c3=0; c4=1; c5=0; c6=-1; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='B211(15),Bm11(8),B2/m11(12);unique a-axis'; end;
            break;
        case 40: // (* Bb11 | hkl,h0l,hk0:h+l, 0kl,0k0,00l:k,l, h00:h *)
            c1=3; c2=4; c3=0; c4=1; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='Bb11(9),B2/b11(15);unique a-axis'; end;
            break;

            //(* I___ *)
        case 41: // (* I-11 | hkl,h0l,hk0:h+k+l, 0kl,0k0,00l:k+l, h00:h *)
            c1=0; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='I211(5),Im11(8),I2/m11(12);unique a-axis'; end;
            break;
        case 42: //  (* Ic11 | hkl,h0l,hk0:h+k+l, 0kl,0k0,00l:k,l, h00:h *)
            c1=0; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1; c13=0;
            //EditSG.Text:='Ic11(9),I2/c11(15);unique a-axis'; end;
            break;
        } // switch ComboBoxMonoclinic
    } // if RadioButtonSysMonoclinic


    //(**************** Triclinic crystal system ************** done *)
    if ( RadioButtonSysTriclinic )
    {
        hmin=-hmax;
        kmin=-kmax;
        lmin=0;
        switch ( ComboBoxTriclinic )
        {
        //(* P___ *)
        case 0: // (* P- | *)
            c1=-1; c2=-1; c3=-1; c4=-1; c5=-1; c6=-1; c7=-1; c8=-1; c9=-1; c10=-1; c11=-1; c12=-1;
            //EditSG.Text:='P1(1),P-1(2)'; end;
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
    alf=params.alpha_deg * M_PI/180.;  // StrToFloat(EditCellAlpha.Text);
    bet=params.beta_deg  * M_PI/180.;  // StrToFloat(EditCellBeta.Text);
    gam=params.gamma_deg * M_PI/180.;  // StrToFloat(EditCellGamma.Text);

#ifdef UseStringGrid12
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
        bet=M_PI/2.;
        kmin=0;
        kmax=0;
        lmin=0;
        lmax=0;
    }
    if ( RadioButtonSys2D )
    {
        c=a*100;
        alf=M_PI/2.;
        bet=M_PI/2.;
        lmin=0;
        lmax=0;
    }
    if ( RadioButtonSysCubic ) // -31686
    {
        b=a;  c=a;  alf=M_PI/2.;    bet=M_PI/2.;    gam=M_PI/2.;    //{NV}-31687
#ifdef UseStringGrid12
        bmax=amax;  cmax=amax;
#endif
        //EditCellb.Text:=FloatToStr(b);  EditCellc.Text:=FloatToStr(c);
    }
    if ( RadioButtonSysTetragonal )
    {
        b=a;  alf=M_PI/2.0;  bet=M_PI/2.0;   gam=M_PI/2.0;
#ifdef UseStringGrid12
        bmax=amax;
#endif
        //EditCellb.Text:=FloatToStr(b);
    }
    if ( RadioButtonSysOrtho )
    {
        alf=M_PI/2.;  bet=M_PI/2.;  gam=M_PI/2.;
    }
    if ( RadioButtonSysHex )
    {
        b=a; alf=M_PI/2.;  bet=M_PI/2.;  gam=2*M_PI/3.;
#ifdef UseStringGrid12
        bmax=amax;
#endif
        //EditCellb.Text:=FloatToStr(b);
    }
    if ( RadioButtonSysTrigonal )   //(* reads angle alfa *)
    {
        if ( trigP=='h' )  //(* hexagonal axes *)
        {
            b=a;
#ifdef UseStringGrid12
            bmax=amax;
#endif
            alf=M_PI/2.;
            bet=M_PI/2.;
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
#endif
            bet=alf;
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
            alf=gam;      bet=M_PI/2.;    gam=M_PI/2.;
            break;
        case 1:  //(* unique axis b *)
            bet=gam;      alf=M_PI/2.;    gam=M_PI/2.;
            break;
        case 2:  //(* unique axis c *)
            alf=M_PI/2.;     bet=M_PI/2.;
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
#endif
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
                        ccase=0;
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
                        ccase = 0;
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
                            latpar1(  1,0)=ct1;
                            latpar1(ct1,1)=h;
                            latpar1(ct1,2)=k;
                            latpar1(ct1,3)=l;
                            latpar1(ct1,4)=mult;
                            //qDebug() << "latpar1:" << ct1 << h << k << l << mult;
                            ct1++;
                        }
                        else //if ( ii==2 )
                        {
                            latpar2(  1,0)=ct1;
                            latpar2(ct1,1)=h;
                            latpar2(ct1,2)=k;
                            latpar2(ct1,3)=l;
                            latpar2(ct1,4)=mult;
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
#endif // USE_ButtonHKLClick



#ifdef USE_fhkl_c
// Nur in prepareCalculation
void CLASSLIB::fhkl_c( int lat, int h, int k, int l,
                       double uca, double ucb, double ucc, double /*ucalpha_deg*/, double /*ucbeta_deg*/,
                       double ucgamma_deg, double &sphno, double &fhkl, double &qhkl, double &qhkl0 ) const
{
    // 'fhkl' wird im rufenden Programm zwar noch in einem Array (latpar?) gespeichert, aber dann
    // nirgendwo weiter verwendet. Daher wird es hier auch nicht gesetzt und kann langfristig
    // als Parameter entfernt werden - sofern Herr Förster keine andere Idee hat.

    int i,nn;
    double a,b,c,/*alf,bet,*/gam,sumc,sums,arg,invd;
    //double /*s2a,*/c1a,/*c2a,c3a,*/c1b,/*c2b,s2b,*/c1g; //,c2g,s2g;
    float r[13][4]; // array[1..12,1..3] of real;

    a=uca;
    b=ucb;
    c=ucc;
    //alf=ucalpha_deg*M_PI/180.;
    //bet=ucbeta_deg*M_PI/180.;
    gam=ucgamma_deg*M_PI/180.;

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
#endif // USE_fhkl_c



#ifdef USE_extinction
// Nur in prepareCalculation verwendet
void CLASSLIB::extinction( int lat, int h, int k, int l, int aniso,
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
#endif // USE_extinction




#ifdef USE_trapzddeltac
//(* *********************** Romberg integration ****************************** *)
//(* returns integral in the limits a and b *)

// LOKALE ROUTINE !!!
//(*** integration routine use trapezoidal rule ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::trapzddeltac( double a, double b, double l, double r, double dbeta, double theta, double phi,
                             double qx, double qy, double qz, double p11, double p12, double p13, double p21,
                             double p22, double p23, double p31, double p32, double p33,
                             double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
                             double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
                             double ax3y, double ax3z, double sigx, double sigy, double sigz,
                             int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                             double *carr1,
                             double &pq, int n, int &trapzddeltac_cnt ) const
{
    int j;
    double x,sump/*,sumf,sumn,sums*/,del;
    double /*fa,fb,fx,*/pa,pb,px/*,na,nb,nx,sa,sb,sx*/;

    if ( n==1 )     /*Z-UPQ1:1878*/
    {
        /*D8L(
        static int ii0=-1, ii1=-1;
        if ( ii0!=i0 || ii1!=i1 )
        {   ii0=i0; ii1=i1;
            qDebug() << "trapzddeltac(0) a,b" << a << b << "i0,i1" << i0 << i1 << "dbeta" << dbeta;
        }
        )*/
        switch ( i0 )
        {
        case 1:   //(* delta-integration *)     /*Z-UPQ1:1879*/
            if ( i1==5 )  //(* general case *)
            {
                qrombchid(l,r,params.sigma,dbeta,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                          carr1, pa );
                qrombchid(l,r,params.sigma,dbeta,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                          carr1, pb );
                pa=pa/(2*M_PI);
                pb=pb/(2*M_PI);
            }
            else
            {
                pa=pow(cos(a),i3)*pow(sin(a),i4);
                pb=pow(cos(b),i3)*pow(sin(b),i4);
            }
            pa=sin(a)*exp(-a*a/(dbeta*dbeta))*pa;
            pb=sin(b)*exp(-b*b/(dbeta*dbeta))*pb;
            break;
        case 2:   //(* norm *)
            //pa:=sin(a)*orfunc(dbeta,a,i2);
            //pb:=sin(b)*orfunc(dbeta,b,i2);
            pa=sin(a)*exp(-a*a/(dbeta*dbeta));
            pb=sin(b)*exp(-b*b/(dbeta*dbeta));
            //pa:=sin(a)*exp(-(a-pi/2)*(a-pi/2)/(dbeta*dbeta));
            //pb:=sin(b)*exp(-(b-pi/2)*(b-pi/2)/(dbeta*dbeta));
            break;
        case 3:   //(* order parameter *)
            //pa:=sin(a)*orfunc(dbeta,a,i2)*(3*cos(a)*cos(a)-1)/2;
            //pb:=sin(b)*orfunc(dbeta,b,i2)*(3*cos(b)*cos(b)-1)/2;
            pa=sin(a)*exp(-a*a/(dbeta*dbeta))*(3*cos(a)*cos(a)-1)/2.;
            pb=sin(b)*exp(-b*b/(dbeta*dbeta))*(3*cos(b)*cos(b)-1)/2.;
            //pa:=sin(a)*exp(-(a-pi/2)*(a-pi/2)/(dbeta*dbeta))*(3*cos(a)*cos(a)-1)/2;
            //pb:=sin(b)*exp(-(b-pi/2)*(b-pi/2)/(dbeta*dbeta))*(3*cos(b)*cos(b)-1)/2;
            break;
        case 4:   //(* cylinder formfactor *)
            qrombchid(l,r,params.sigma,dbeta,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pa );
            qrombchid(l,r,params.sigma,dbeta,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pb );
            pa=sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2*M_PI);
            pb=sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2*M_PI);
            break;
        case 5:   //(* unit cell rotation *)
            // Bei den Anisotropic Gaussian Tests genutzt
            qrombchid(l,r,params.sigma,dbeta,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pa );
            qrombchid(l,r,params.sigma,dbeta,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pb );
            pa=8*sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2*M_PI*pow(M_PI,3/2.)*sigx*sigy*sigz);
            pb=8*sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2*M_PI*pow(M_PI,3/2.)*sigx*sigy*sigz);
            break;
        case 6:   //(* disk formfactor *)
            qrombchid(l,r,params.sigma,dbeta,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pa );
            qrombchid(l,r,params.sigma,dbeta,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pb );
            pa=sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2*M_PI);
            pb=sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2*M_PI);
            //pa:=sin(a)*exp(-(a-pi/2)*(a-pi/2)/(dbeta*dbeta))*pa/(2*pi);
            //pb:=sin(b)*exp(-(b-pi/2)*(b-pi/2)/(dbeta*dbeta))*pb/(2*pi);
            break;
        case 99:   // TEST x²*y³
            qrombchid(l,r,params.sigma,dbeta,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pa );
            qrombchid(l,r,params.sigma,dbeta,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                      carr1, pb );
            pa = pa * a * a * a;
            pb = pb * b * b * b;
            //qDebug() << "TEST trapzddeltac (1)" << pa << pb;
            break;
        } // switch i0
        if ( i0 == 99 )
            pq = pa + pb;
        else
            pq=0.5*(b-a)*(pa+pb);
        DSM(
        //if ( fabs(pq) < 1.0e-170 )
        //    qDebug() << "trapzddeltac(1)" << pa << pb << pq << "i?" << i0 << i1 << i2 << i3 << i4;
        )
        trapzddeltac_cnt=1;
        DCNT( dbgCount++; )
    }
    else
    {
        del=(b-a)/trapzddeltac_cnt;
        DCNT(
        #ifdef DBGLIMIT
                    if ( ++dbgCount > DBGLIMIT ) return;
        #else
                    dbgCount++;
        #endif
                )
        x=a+0.5*del;
        sump=0.0;
        for ( j=1; j<=trapzddeltac_cnt; j++ )
        {
            //if i0=1 then px:=sin(x)*orfunc(dbeta,x,i2)*power(cos(x),i3)*power(sin(x),i4);
            //if i0=2 then px:=sin(x)*orfunc(dbeta,x,i2);
            //if i0=3 then px:=sin(x)*orfunc(dbeta,x,i2)*(3*cos(x)*cos(x)-1)/2;
            //if i0=2 then px:=sin(x)*exp(-(x-pi/2)*(x-pi/2)/(dbeta*dbeta));
            //if i0=3 then px:=sin(x)*exp(-(x-pi/2)*(x-pi/2)/(dbeta*dbeta))*(3*cos(x)*cos(x)-1)/2;
            switch ( i0 )
            {
            case 1:
                px=sin(x)*exp(-x*x/(dbeta*dbeta))*pow(cos(x),i3)*pow(sin(x),i4);
                break;
            case 2:
                px=sin(x)*exp(-x*x/(dbeta*dbeta));
                break;
            case 3:
                px=sin(x)*exp(-x*x/(dbeta*dbeta))*(3*cos(x)*cos(x)-1)/2.;
                break;
            case 4:   //(* cylinder formfactor *)
                qrombchid(l,r,params.sigma,dbeta,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                          carr1, px );
                px=sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2*M_PI);
                break;
            case 5:   //(* unit cell rotation *)
                // Bei den Anisotropic Gaussian Tests genutzt
                qrombchid(l,r,params.sigma,dbeta,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                          carr1, px );
                px=8*sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2*M_PI*pow(M_PI,3/2.)*sigx*sigy*sigz);
                break;
            case 6:   //(* disk formfactor *)
                qrombchid(l,r,params.sigma,dbeta,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                          carr1, px );
                px=sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2*M_PI);
                //px:=sin(x)*exp(-(x-pi/2)*(x-pi/2)/(dbeta*dbeta))*px/(2*pi);
                break;
            case 99:   // TEST x²*y³
                // Bei den Anisotropic Gaussian Tests genutzt
                qrombchid(l,r,params.sigma,dbeta,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,
                          carr1, px );
                px = px * x * x * x;
                //if ( j == 1 ) qDebug() << "TEST trapzddeltac (n)" << px;
                break;
            } // switch i0
            sump += px;
            x += del;
        } // for j
        if ( i0 == 99 )
            pq = sump / trapzddeltac_cnt;
        else
            pq = 0.5*(pq+(b-a)*sump/trapzddeltac_cnt);
        D8L( qDebug() << "trapzddeltac(2)" << sump << trapzddeltac_cnt << del << pq );
        DSM(
        if ( fabs(pq) < 1.0e-170 )
            qDebug() << "trapzddeltac(3)" << sump << trapzddeltac_cnt << del << pq << "i?" << i0 << i1 << i2 << i3 << i4;
        )
        trapzddeltac_cnt = 2*trapzddeltac_cnt;
    }
}
#endif // USE_trapzddeltac



#ifdef USE_polint
//(*** interpolation routine ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::polint( double *xa/*RealArrayNP*/, double *ya/*RealArrayNP*/,
                       int n, double /*x*/, double &y, double &dy, const char * /* DSM(dbg) */ ) const
{
    // Parameter x ist bei allen Aufrufen immer 0.0 - daher kann der auch weggelassen werden!
    // Parameter xa= 1 0.25 0.0625 0.015625 0.00390625 bei allen Aufrufen in qrombdeltac() und qrombchid()
    int ns,m,i;
    double w,hp,ho,dift,dif,den;
    RealArrayNP c, d; // ^RealArrayNP;     =array[1..np=5] of extended;

    //D8L( qDebug() << "   polint Anf xa=" <<xa[1]<<xa[2]<<xa[3]<<xa[4]<<xa[5] << "ya=" <<ya[1]<<ya[2]<<ya[3]<<ya[4]<<ya[5] << dbg );
    ns=1;
    dif=fabs(xa[1]);           // fabs(x-xa[1])
    for ( i=1/*?? 2*/; i<=n; i++ )
    {
        dift=fabs(xa[i]);      // fabs(x-xa[i])
        if ( dift < dif )
        {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    //D8L( qDebug() << "   polint Inner" << dif << ns << ya[ns] );
    y=ya[ns];
    ns=ns-1;
    for ( m=1; m<=n-1; m++ )
    {
        for ( i=1; i<=n-m; i++ )
        {
            ho=xa[i];       // -x
            hp=xa[i+m];     // -x
            w=c[i+1]-d[i];
            den=ho-hp;
            //if ( fabs(den) > eps5 ) den=w/den;
            if ( den != 0 ) den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        } // for i
        if ( 2*ns < n-m )
            dy=c[ns+1];
        else
        {
            dy=d[ns];
            ns=ns-1;
        }
        y=y+dy;
        //qDebug() << m << "c=" << arr2str(c,n) << "d=" << arr2str(d,n) << y << dy;
        //D8L( qDebug() << "   polint m" << m << y << dy );
    } // for m

    DSM( if ( fabs(dy) < 1.0e-200 && fabs(dy) > 0 )
            qDebug() << "polint ya=" <<ya[1]<<ya[2]<<ya[3]<<ya[4]<<ya[5] << dbg
                     << y << dy; )

}
#endif // USE_polint



#ifdef USE_qrombdeltac

// Old:
//procedure qrombdeltac(r,sigma,dbeta,theta,phi,qx,qy,qz: extended;
//                   qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz: extended;
//                   ordis,dim,i0,i1,i2,i3,i4: integer;
//                   var pq: extended);
// 20220311:
//procedure qrombdeltac(l,r,p1,sigma,dbeta,theta,phi,qx,qy,qz: extended;
//                   qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz: extended;
//                   ordis,dim,i0,i1,i2,i3,i4: integer;
//                   var carr1: Array of extended;
//                   var pq: extended);


#ifdef __CUDACC__
__host__ __device__
#endif
// Parameter r darf nicht wegfallen, da er meistens <params.radius> ist, aber auch mal <params.length> oder lokale Variablen übergeben werden
void CLASSLIB::qrombdeltac( double l, double r, /*p1*/ /*sigma*/ /*dbeta*/ double theta, double phi, double qx, double qy, double qz,
                            double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
                            double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
                            double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
                            int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                            double *carr1,
                            double &pq ) const
{
    const int jmax=20;
    const int k=5;
    const double prmin = 0.00001;       //(* minimum pr(delta)-value *)

    DCL( static long dbgCountMin=0;
         static long dbgCountMax=0;
         static double dbgCountSum = 0 );
    DTL( static time_t startzeit = 0;
         static long   dbgCntDiff = 0 );

    int i,j;
    //int midpntit;
    int trapzddeltac_cnt;
    double dssp,delmax/*,norm*/;
    RealArrayJMAXP hp, sp; // ^RealArrayJMAXP;  =array[1..jmaxp=21] of extended;
    RealArrayNP    cp, dp; // ^RealArrayNP;     =array[1..np=5] of extended;
    float  alim,blim,p11,p12,p13,p21,p22,p23,p31,p32,p33;     //(* integration limits *)  x² bringt nichts

    DTL( if ( startzeit == 0 ) startzeit = time(nullptr) );

    //D8L( qDebug() << "qrombdeltac(" << r << theta << phi << "qx/y/z" << qx << qy << qz << "q?n" << qxn << qyn << qzn << qhkl << "ax?n" << ax1n << ax2n << ax3n
    //                               << "ax1?" << ax1x << ax1y << ax1z << "ax2?" << ax2x << ax2y << ax2z << "ax3?" << ax3x << ax3y << ax3z
    //                               << "sig?" << sigx << sigy << sigz << "ordis" << ordis << dim << "i?" << i0 << i1 << i2 << i3 << i4 << "&pq )" << params.dbeta );

    //memset( hp, 0, sizeof(hp) );
    //memset( sp, 0, sizeof(sp) );
    //memset( cp, 0, sizeof(cp) );
    //memset( dp, 0, sizeof(dp) );

    // params.dbeta ist die globale Variable in Grad.
    double dbeta = params.dbeta * M_PI / 180.;
    theta=theta*M_PI/180.;
    phi=phi*M_PI/180.;

/*
Debug:    qrombdeltac( 10 1.23096 -2.61799 qx/y/z 0.625 0 1e-20 ... i? 5 6 0 0 0  )
Debug: ** qrombchid 0 5 3.13487e-224 3.80432e-209 16 249592
Debug: ** qrombchid 0 5 2.38098e-215 3.80432e-209 16 249613
Debug: ** qrombchid 0 5 1.17959e-215 3.80432e-209 16 249637
*/
    /*if ( i0 < 99 )
    {
        if ( qx >= 0.6 || qy >= 0.6 ) { pq=1.0; return; } // hilft etwas, ist aber keine Lösung!
    }*/

    DCL( qDebug() << "   qrombdeltac(" << r << theta << phi << "qx/y/z" << qx << qy << qz << "..." // nur Kurzform
                                       << "i?" << i0 << i1 << i2 << i3 << i4 << " ) sigma:" << params.sigma );

    // aus fcc*.cu kommen: ,FCC.ordis,  3, 5, 6, 0, 0, 0,
    //                     ,int ordis,dim,i0,i1,i2,i3,i4,

    //(* search for maximum integration angle *)
    delmax = 1; // to avoid compiler warning
    // D8L: i2=0 -> delmax = dbeta * sqrt( log( 1.0/0.00001 ) )
    //                               2.2360679775
    if ( (i2==0) || (i2==2) || (i2==3) || (i2==5) ) delmax=dbeta*sqrt(log(1./prmin));  //(* Gaussian, Onsager, Maier-Saupe, Laguerre *)
    if ( i2==1 ) delmax=dbeta*log(1./prmin);        //(* Exponential *)
    //if i2=2 then delmax:=arcsin(dbeta*ln(1/prmin));
    if ( i2==4 ) delmax=dbeta;                    //(* cut-off *)
    if ( (i2==7) || (i2==8) || (i2==9) || (i2==10) || (i2==11) || (i2==12) ) delmax=M_PI/2.; //(* isotropic, mirrored distributions *)
    // i2 = 0..12, 6=s.u.

    if ( delmax > M_PI/2. ) delmax=M_PI/2.;
    alim=0.0;
    blim=delmax;
    if ( i0 == 99 )
    {   // TEST x²*y³
        blim = 4.0;
    }

    DCL( qDebug() << "               dbeta" << dbeta << "delmax" << delmax << "a/blim" << alim << blim; )

    if ( i1==12 ) //   (* for cubes *)
    {
        alim=eps4;
        blim=M_PI/2.0-eps4;
    }

    //(* for disks: mirrored Gaussian *)
    /*if ( qzn==2 || i0==6 )    - im Pascalcode auskommentiert
    {
       alim=M_PI/2.-delmax;
       blim=M_PI/2.;
    }*/

    p11=-cos(phi)*cos(theta);   /*Z-UPQ1:2017*/
    p12=sin(phi);
    p13=cos(phi)*sin(theta);
    p21=-cos(phi);
    p22=-sin(phi)*cos(theta);
    p23=sin(phi)*sin(theta);
    p31=-sin(theta);
    p32=0;
    p33=-cos(theta);

    hp[1]=1.0;
    for ( j=1; j<=jmax; j++ )
    {
        DCNT(
        #ifdef DBGLIMIT
                    if ( ++dbgCount > DBGLIMIT ) break;
        #else
                    dbgCount++;
        #endif
                )

        trapzddeltac(alim,blim, l, r,/*p1*/ /*sigma*/ dbeta,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                     ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,
                     ordis,dim,i0,i1,i2,i3,i4, carr1/*neu*/, sp[j],j, trapzddeltac_cnt );
        //D8L(
        //if ( j == 1 ) qDebug() << "sp[j]" << j << sp[j] << alim << blim << dbgCount;
        //)
        /*if ( j == 1 && dbgCount < 1000 )
        {
            std::cerr.precision(10);
            std::cerr << "sp[1] = " << sp[j] << ", DbgCount = " << dbgCount << std::endl;
        }*/
        if ( j>=k )     /*Z-UPQ1:2034*/
        {
            for ( i=1; i<=k; i++ )      // k=5
            {
                cp[i]=hp[j-k+i];
                dp[i]=sp[j-k+i];
            }
            polint(cp,dp,k,0.0,pq,dssp,"qrombdeltac");
            D8L( qDebug() << "qrombdeltac inner" << j << dssp << pq );
            if ( fabs(dssp) < eps4*fabs(pq) ) break;
            if ( fabs(pq) < 1.0e-100 ) break;  // nicht im Pascal-Programm
        }
        sp[j+1]=sp[j];
        hp[j+1]=0.25*hp[j];
    }
    DCL( if ( dbgCount > 20 )    // Init weglassen
         {
           dbgCountSum += dbgCount;
           if ( dbgCount > dbgCountMax ) dbgCountMax = dbgCount;
           if ( dbgCount < dbgCountMin || dbgCountMin == 0 ) dbgCountMin = dbgCount;
         } );
    D8L( qDebug() << "qrombdeltac: trapzddeltac_cnt =" << trapzddeltac_cnt );
    DCL(
    qDebug() << "DBGCNT ****************************   qrombdeltac" << dbgCount
                << "*********" << dbgCountMin << dbgCountMax << dbgCountSum
                DTL(<< " t:" << time(nullptr) - startzeit);
    //if ( dbgCount > 500000 || dbgCount == dbgCountMin/*zum Wertevergleich*/ )
    //    qDebug() << "qrombdeltac(" << r << theta << phi << "qx/y/z" << qx << qy << qz << "q?n" << qxn << qyn << qzn << qhkl << "ax?n" << ax1n << ax2n << ax3n
    //             << "ax1?" << ax1x << ax1y << ax1z << "ax2?" << ax2x << ax2y << ax2z << "ax3?" << ax3x << ax3y << ax3z
    //             << "sig?" << sigx << sigy << sigz << "ordis" << ordis << dim << "i?" << i0 << i1 << i2 << i3 << i4 << pq << ")";
    )

}
#endif // USE_qrombdeltac



#ifdef USE_qrombchid
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::qrombchid( double l, double r, double sigma, double dbeta, double delta, double theta, double phi,
                          double qx, double qy, double qz, double p11, double p12, double p13, double p21,
                          double p22, double p23, double p31, double p32, double p33,
                          double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
                          double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
                          double ax3y, double ax3z, double sigx, double sigy, double sigz,
                          int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                          double *carr1, // NEU
                          double &pq ) const
{
    //const double eps=1.0e-4;
    const int jmax=20;
    const int k=5;

    int i,j, trapzdchid_cnt;
    double dssp/*,delmax*/;
    RealArrayJMAXP hp, sp; // ^RealArrayJMAXP;  =array[1..jmaxp=21] of extended;
    RealArrayNP    cp, dp; // ^RealArrayNP;     =array[1..np=5] of extended;
    double /*float */ alim,blim;     //(* integration limits *)      x² bringt nichts

    //qDebug() << "qrombchid(" << r << sigma << dbeta << delta << theta << phi << qx << qy << qz << p11 << p12 << p13 << p21
    //         << p22 << p23 << p31 << p32 << p33 << qxn << qyn << qzn << qhkl << ax1n << ax2n << ax3n << ax1x << ax1y
    //         << ax1z << ax2x << ax2y << ax2z << ax3x << ax3y << ax3z << sigx << sigy << sigz << "ordis" << ordis << dim
    //         << i0 << i1 << i2 << i3 << i4 << "&pq" << " );";

    alim=0.0;
    blim=2*M_PI;
    if ( i0 == 99 )
    {   // TEST x²*y³
        blim = 2.0;
    }

    hp[1]=1.0;
    for ( j=1; j<=jmax; j++ )
    {
        DCL(
#ifdef DBGLIMIT
        if ( ++dbgCount > DBGLIMIT ) return;
#else
        dbgCount++;
#endif
                )
        //trapzdchid(alim,blim,l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,dbeta,phi,theta,qxn,qyn,qzn,q,maxit,i0,i1,i2,i3,i4,sp^[j],j);
        trapzdchid(alim,blim, l, r,sigma,dbeta,delta,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                   ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,
                   ordis,dim,i0,i1,i2,i3,i4,carr1/*NEU*/, sp[j],j, trapzdchid_cnt );
        if ( j >= k )
        {
            for ( i=1; i<=k; i++ )
            {
                cp[i]=hp[j-k+i];
                dp[i]=sp[j-k+i];
            }
            polint(cp,dp,k,0.0,pq,dssp,"qrombchid");
            if ( fabs(dssp) < eps4*fabs(pq) ) break;
            if ( fabs(pq) < 1.0e-300 )
            {
                //qDebug() << " ... qrombchid" << dssp << pq;
                break; // TODO ?
            }
        }
        sp[j+1]=sp[j];
        hp[j+1]=0.25*hp[j];
    } // for j
} /* qrombchid() */
#endif // USE_qrombchid



#ifdef USE_trapzdchid
#ifdef __CUDACC__
__host__ __device__
#endif
void CLASSLIB::trapzdchid( double a, double b, double l, double r, double sigma, double dbeta, double delta, double /*theta*/, double /*phi*/,
                           double qx, double qy, double qz, double p11, double p12, double p13, double p21, double p22, double p23,
                           double p31, double /*p32*/, double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n,
                           double ax2n, double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
                           double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
                           int /*ordis*/, int /*dim*/, int /*i0*/, int i1, int /*i2*/, int i3, int i4,
                           double *carr1,  // NEU
                           double &pq, int n, int &trapzdchid_cnt ) const
{
    /*Z0311= gilt für "20220311 - upq1.pas" */
    //const double eps = 0.0000000001;   /*Z0311=126  wird zu eps9 (global)*/
    //const double eps1 = 0.00001;   /*Z0311=127*/
    // label 700,701;   /*Z0311=129*/

    int    i,j;
    double x,tnm,sump,/*sumf,sumn,sums,*/del/*,pqsum,oldpqsum,delser*/;
    double /*fa,fb,fx,*/pa,pb,px,/*na,nb,nx,sa,sb,sx,*/epsi;
    double z,a1,/*x1z,*/arga,argb,argx,l1,l2,l3;//,qxhkl,qyhkl,qzhkl;
    double mla1,mla2,/*mla3,*/mlb1,mlb2,/*mlb3,*/mlx1,mlx2; //,mlx3;
    //double ma11,ma12,ma13,ma21,ma22,ma23,ma31,ma32,ma33;
    //double mb11,mb12,mb13,mb21,mb22,mb23,mb31,mb32,mb33;
    //double mx11,mx12,mx13,mx21,mx22,mx23,mx31,mx32,mx33;
    double dqxa,dqya,dqza,dqxb,dqyb,dqzb,dqs1a,dqs2a,dqs3a,dqs1b,dqs2b,dqs3b;
    double dqxx,dqyx,dqzx,dqs1x,dqs2x,dqs3x; //,vola,volb,volx;
    double r11a,r12a,r13a,r21a,r22a,r23a,r31a,r32a,r33a;
    double r11b,r12b,r13b,r21b,r22b,r23b,r31b,r32b,r33b;
    double r11x,r12x,r13x,r21x,r22x,r23x,r31x,r32x,r33x;
    double ax1xa,ax1ya,ax1za,ax1xb,ax1yb,ax1zb,ax1xx,ax1yx,ax1zx;
    double ax2xa,ax2ya,ax2za,ax2xb,ax2yb,ax2zb,ax2xx,ax2yx,ax2zx;
    double ax3xa,ax3ya,ax3za,ax3xb,ax3yb,ax3zb,ax3xx,ax3yx,ax3zx;
    double qn,qxhkla,qyhkla,qzhkla,qxhklb,qyhklb,qzhklb,qxhklx,qyhklx,qzhklx;
    //double aexa,aeya,aeza,bexa,beya,beza,cexa,ceya,ceza,aexb,aeyb,aezb,bexb,beyb,bezb,cexb,ceyb,cezb;
    //double aexx,aeyx,aezx,bexx,beyx,bezx,cexx,ceyx,cezx;
    double pa1,pa2,pa3,pb1,pb2,pb3,px1,px2,px3,ella,ellb,ellc;
    double qxl,qyl,qnna,qnnb,qnnx; //,argas,argbs,z1,z2,z3,z4,z5,z6;
    double arga1,arga2,arga3,argb1,argb2,argb3,argx1,argx2,argx3;

    arga = 1; // to avoid compiler warnings
    argb = 1;

    if ( n==1 )
    {   /*Z0311=156*/
        if ( i1==1 )
        {    /* cylinder, general case */   /*Z0311=157*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=158*/
            a1 = 1/(2*z*(z-1));   /*Z0311=159*/
            mla1 = p11*cos(a)*sin(delta)+p12*sin(a)*sin(delta)+p13*cos(delta);   /*Z0311=160*/
            mla2 = p21*sin(a)*sin(delta)+p22*cos(a)*sin(delta)+p23*cos(delta);   /*Z0311=161*/
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta); */   /*Z0311=162*/
            mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);   /*Z0311=163*/
            mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);   /*Z0311=164*/
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta); */   /*Z0311=165*/
            arga = (qx*mla1+qy*mla2+eps9)*l/(z+1);   /*Z0311=166*/     // TODO len-wave
            argb = (qx*mlb1+qy*mlb2+eps9)*l/(z+1);   /*Z0311=167*/     // TODO len-wave
            if ( i3==0 )
            { /* P(q) */   /*Z0311=168*/
                a1 = 1/(2*z*(z-1));   /*Z0311=169*/
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.0));   /*Z0311=170*/
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.0));   /*Z0311=171*/
            }   /*Z0311=172*/
            if ( i3==1 )
            { /* F(q) */   /*Z0311=173*/
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1+arga*arga,z/2.0);   /*Z0311=174*/
                pa = pa1*pa1;   /*Z0311=175*/
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1+argb*argb,z/2.0);   /*Z0311=176*/
                pb = pb1*pb1;   /*Z0311=177*/
            }   /*Z0311=178*/
        }   /*Z0311=179*/
        if ( i1==2 )
        {   /* cylinder, x-axis */   /*Z0311=180*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=181*/
            arga = (qx*cos(delta)-qy*sin(a)*sin(delta)+eps9)*l/(z+1);   /*Z0311=182*/
            argb = (qx*cos(delta)-qy*sin(b)*sin(delta)+eps9)*l/(z+1);   /*Z0311=183*/
            if ( i3==0 )
            { /* P(q) */   /*Z0311=184*/
                a1 = 1/(2*z*(z-1));   /*Z0311=185*/
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.0));   /*Z0311=186*/
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.0));   /*Z0311=187*/
            }   /*Z0311=188*/
            if ( i3==1 )
            { /* F(q) */   /*Z0311=189*/
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1+arga*arga,z/2.0);   /*Z0311=190*/
                pa = pa1*pa1;   /*Z0311=191*/
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1+argb*argb,z/2.0);   /*Z0311=192*/
                pb = pb1*pb1;   /*Z0311=193*/
            }   /*Z0311=194*/
        }   /*Z0311=195*/
        if ( i1==3 )
        {   /* cylinder, y-axis */   /*Z0311=196*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=197*/
            arga = (qx*sin(a)*sin(delta)+qy*cos(delta)+eps9)*l/(z+1);   /*Z0311=198*/
            argb = (qx*sin(b)*sin(delta)+qy*cos(delta)+eps9)*l/(z+1);   /*Z0311=199*/
            if ( i3==0 )
            { /* P(q) */   /*Z0311=200*/
                a1 = 1/(2*z*(z-1));   /*Z0311=201*/
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.0));   /*Z0311=202*/
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.0));   /*Z0311=203*/
            }   /*Z0311=204*/
            if ( i3==1 )
            { /* F(q) */   /*Z0311=205*/
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1+arga*arga,z/2.0);   /*Z0311=206*/
                pa = pa1*pa1;   /*Z0311=207*/
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1+argb*argb,z/2.0);   /*Z0311=208*/
                pb = pb1*pb1;   /*Z0311=209*/
            }   /*Z0311=210*/
        }   /*Z0311=211*/
        if ( i1==4 )
        {   /* cylinder, -z-axis */   /*Z0311=212*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=213*/
            arga = (qx*sin(a)*sin(delta)-qy*cos(a)*sin(delta)+eps9)*l/(z+1);   /*Z0311=214*/
            argb = (qx*sin(b)*sin(delta)-qy*cos(b)*sin(delta)+eps9)*l/(z+1);   /*Z0311=215*/
            if ( (i3==0) )
            { /* P(q) */   /*Z0311=216*/
                a1 = 1/(2*z*(z-1));   /*Z0311=217*/
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.0));   /*Z0311=218*/
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.0));   /*Z0311=219*/
            }   /*Z0311=220*/
            if ( (i3==1) )
            { /* F(q) */   /*Z0311=221*/
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1+arga*arga,z/2.0);   /*Z0311=222*/
                pa = pa1*pa1;   /*Z0311=223*/
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1+argb*argb,z/2.0);   /*Z0311=224*/
                pb = pb1*pb1;   /*Z0311=225*/
            }   /*Z0311=226*/
        }   /*Z0311=227*/
        if ( (i1==5) )
        {    /* general series expansion */   /*Z0311=228*/
            mla1 = p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);   /*Z0311=229*/
            mla2 = p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);   /*Z0311=230*/
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta); */   /*Z0311=231*/
            mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);   /*Z0311=232*/
            mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);   /*Z0311=233*/
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta); */   /*Z0311=234*/
            pa = pow(mla1,i3)*pow(mla2,i4);   /*Z0311=235*/
            pb = pow(mlb1,i3)*pow(mlb2,i4);   /*Z0311=236*/
        }   /*Z0311=237*/
        if ( i1==6 )
        {    /* unit cell rotation */   /*Z0311=238*/
            dqxa = qx-qhkl*(p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));   /*Z0311=239*/
            dqya = qy-qhkl*(p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));   /*Z0311=240*/
            dqza = qz-qhkl*(p31*cos(a)*sin(delta)+p33*cos(delta));   /*Z0311=241*/
            dqxb = qx-qhkl*(p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));   /*Z0311=242*/
            dqyb = qy-qhkl*(p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));   /*Z0311=243*/
            dqzb = qz-qhkl*(p31*cos(b)*sin(delta)+p33*cos(delta));   /*Z0311=244*/
            dqs1a = (dqxa*ax1x+dqya*ax1y+dqza*ax1z)/(ax1n*sigx);;   /*Z0311=245*/
            dqs2a = (dqxa*ax2x+dqya*ax2y+dqza*ax2z)/(ax2n*sigy);   /*Z0311=246*/
            dqs3a = (dqxa*ax3x+dqya*ax3y+dqza*ax3z)/(ax3n*sigz);   /*Z0311=247*/
            dqs1b = (dqxb*ax1x+dqyb*ax1y+dqzb*ax1z)/(ax1n*sigx);   /*Z0311=248*/
            dqs2b = (dqxb*ax2x+dqyb*ax2y+dqzb*ax2z)/(ax2n*sigy);   /*Z0311=249*/
            dqs3b = (dqxb*ax3x+dqyb*ax3y+dqzb*ax3z)/(ax3n*sigz);   /*Z0311=250*/
            pa = exp(-4*(dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a)/M_PI);   /*Z0311=251*/
            pb = exp(-4*(dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b)/M_PI);   /*Z0311=252*/
        }   /*Z0311=253*/
        if ( i1==7 )
        {    /* fiber rotation */   /*Z0311=254*/
            /* rotation axis, director */   /*Z0311=255*/
            /* l1:=sin(theta)*cos(phi); */   /*Z0311=256*/
            /* l2:=sin(theta)*sin(phi); */   /*Z0311=257*/
            /* l3:=-cos(theta); */   /*Z0311=258*/
            l1 = r;   /*Z0311=259*/
            l2 = sigma;   /*Z0311=260*/
            l3 = dbeta;   /*Z0311=261*/

            /* rotation matrix Ra */   /*Z0311=263*/
            r11a = cos(a)+(1-cos(a))*l1*l1;   /*Z0311=264*/
            r12a = -l3*sin(a)+(1-cos(a))*l1*l2;   /*Z0311=265*/
            r13a = -l2*sin(a)+(1-cos(a))*l1*l3;   /*Z0311=266*/
            r21a = l3*sin(a)+(1-cos(a))*l1*l2;   /*Z0311=267*/
            r22a = cos(a)+(1-cos(a))*l2*l2;   /*Z0311=268*/
            r23a = l1*sin(a)+(1-cos(a))*l2*l3;   /*Z0311=269*/
            r31a = l2*sin(a)+(1-cos(a))*l1*l3;   /*Z0311=270*/
            r32a = -l1*sin(a)+(1-cos(a))*l2*l3;   /*Z0311=271*/
            r33a = cos(a)+(1-cos(a))*l3*l3;   /*Z0311=272*/

            /* rotation matrix Rb */   /*Z0311=274*/
            r11b = cos(b)+(1-cos(b))*l1*l1;   /*Z0311=275*/
            r12b = -l3*sin(b)+(1-cos(b))*l1*l2;   /*Z0311=276*/
            r13b = -l2*sin(b)+(1-cos(b))*l1*l3;   /*Z0311=277*/
            r21b = l3*sin(b)+(1-cos(b))*l1*l2;   /*Z0311=278*/
            r22b = cos(b)+(1-cos(b))*l2*l2;   /*Z0311=279*/
            r23b = l1*sin(b)+(1-cos(b))*l2*l3;   /*Z0311=280*/
            r31b = l2*sin(b)+(1-cos(b))*l1*l3;   /*Z0311=281*/
            r32b = -l1*sin(b)+(1-cos(b))*l2*l3;   /*Z0311=282*/
            r33b = cos(b)+(1-cos(b))*l3*l3;   /*Z0311=283*/

            /* scattering vector hkl */   /*Z0311=285*/
#ifdef PascalComment
            qxhkl = (2*pi/1)*(p11*i2+p21*i3+p31*i4);   /*Z0311=286*/
            qyhkl = (2*M_PI/1.0)*(p12*i2+p22*i3+p32*i4);   /*Z0311=287*/
            qzhkl = (2*M_PI/1.0)*(p13*i2+p23*i3+p33*i4);
#endif

            /* rotate scattering vector */   /*Z0311=290*/
            qxhkla = r11a*qxn+r12a*qyn+r13a*qzn;   /*Z0311=291*/
            qyhkla = r21a*qxn+r22a*qyn+r23a*qzn;   /*Z0311=292*/
            qzhkla = r31a*qxn+r32a*qyn+r33a*qzn;   /*Z0311=293*/
            qxhklb = r11b*qxn+r12b*qyn+r13b*qzn;   /*Z0311=294*/
            qyhklb = r21b*qxn+r22b*qyn+r23b*qzn;   /*Z0311=295*/
            qzhklb = r31b*qxn+r32b*qyn+r33b*qzn;   /*Z0311=296*/

#ifdef PascalComment
            /* rotated a */   /*Z0311=298*/
            aexa = r11a*p11+r12a*p12+r13a*p13;   /*Z0311=299*/
            aeya = r21a*p11+r22a*p12+r23a*p13;   /*Z0311=300*/
            aeza = r31a*p11+r32a*p12+r33a*p13;   /*Z0311=301*/
            aexb = r11b*p11+r12b*p12+r13b*p13;   /*Z0311=302*/
            aeyb = r21b*p11+r22b*p12+r23b*p13;   /*Z0311=303*/
            aezb = r31b*p11+r32b*p12+r33b*p13;   /*Z0311=304*/

            /* rotated b */   /*Z0311=306*/
            bexa = r11a*p21+r12a*p22+r13a*p23;   /*Z0311=307*/
            beya = r21a*p21+r22a*p22+r23a*p23;   /*Z0311=308*/
            beza = r31a*p21+r32a*p22+r33a*p23;   /*Z0311=309*/
            bexb = r11b*p21+r12b*p22+r13b*p23;   /*Z0311=310*/
            beyb = r21b*p21+r22b*p22+r23b*p23;   /*Z0311=311*/
            bezb = r31b*p21+r32b*p22+r33b*p23;   /*Z0311=312*/

            /* rotated c */   /*Z0311=314*/
            cexa = r11a*p31+r12a*p32+r13a*p33;   /*Z0311=315*/
            ceya = r21a*p31+r22a*p32+r23a*p33;   /*Z0311=316*/
            ceza = r31a*p31+r32a*p32+r33a*p33;   /*Z0311=317*/
            cexb = r11b*p31+r12b*p32+r13b*p33;   /*Z0311=318*/
            ceyb = r21b*p31+r22b*p32+r23b*p33;   /*Z0311=319*/
            cezb = r31b*p31+r32b*p32+r33b*p33;   /*Z0311=320*/

            ma11 = aexa;   /*Z0311=322*/
            ma12 = aeya;   /*Z0311=323*/
            ma13 = aeza;   /*Z0311=324*/
            ma21 = bexa;   /*Z0311=325*/
            ma22 = beya;   /*Z0311=326*/
            ma23 = beza;   /*Z0311=327*/
            ma31 = cexa;   /*Z0311=328*/
            ma32 = ceya;   /*Z0311=329*/
            ma33 = ceza;   /*Z0311=330*/

            mb11 = aexb;   /*Z0311=332*/
            mb12 = aeyb;   /*Z0311=333*/
            mb13 = aezb;   /*Z0311=334*/
            mb21 = bexb;   /*Z0311=335*/
            mb22 = beyb;   /*Z0311=336*/
            mb23 = bezb;   /*Z0311=337*/
            mb31 = cexb;   /*Z0311=338*/
            mb32 = ceyb;   /*Z0311=339*/
            mb33 = cezb;   /*Z0311=340*/

            /* sigma is unit cell volume */   /*Z0311=342*/
            /* reciprocal space base vectors */   /*Z0311=343*/
            vola = aexa*(beya*ceza-beza*ceya)+aeya*(beza*cexa-bexa*ceza)+aeza*(bexa*ceya-beya*cexa);   /*Z0311=344*/
            ma11 = (beya*ceza-beza*ceya)/vola;   /*Z0311=345*/
            ma12 = (beza*cexa-bexa*ceza)/vola;   /*Z0311=346*/
            ma13 = (bexa*ceya-beya*cexa)/vola;   /*Z0311=347*/
            ma21 = (aeza*ceya-aeya*ceza)/vola;   /*Z0311=348*/
            ma22 = (aexa*ceza-aeza*cexa)/vola;   /*Z0311=349*/
            ma23 = (aeya*cexa-aexa*ceya)/vola;   /*Z0311=350*/
            ma31 = (aeya*beza-aeza*beya)/vola;   /*Z0311=351*/
            ma32 = (aeza*bexa-aexa*beza)/vola;   /*Z0311=352*/
            ma33 = (aexa*beya-aeya*bexa)/vola;   /*Z0311=353*/

            volb = aexb*(beyb*cezb-bezb*ceyb)+aeyb*(bezb*cexb-bexb*cezb)+aezb*(bexb*ceyb-beyb*cexb);   /*Z0311=355*/
            mb11 = (beyb*cezb-bezb*ceyb)/volb;   /*Z0311=356*/
            mb12 = (bezb*cexb-bexb*cezb)/volb;   /*Z0311=357*/
            mb13 = (bexb*ceyb-beyb*cexb)/volb;   /*Z0311=358*/
            mb21 = (aezb*ceyb-aeyb*cezb)/volb;   /*Z0311=359*/
            mb22 = (aexb*cezb-aezb*cexb)/volb;   /*Z0311=360*/
            mb23 = (aeyb*cexb-aexb*ceyb)/volb;   /*Z0311=361*/
            mb31 = (aeyb*bezb-aezb*beyb)/volb;   /*Z0311=362*/
            mb32 = (aezb*bexb-aexb*bezb)/volb;   /*Z0311=363*/
            mb33 = (aexb*beyb-aeyb*bexb)/volb;   /*Z0311=364*/

            qxhkla = (2*M_PI/1.0)*(ma11*i2+ma21*i3+ma31*i4);   /*Z0311=366*/
            qyhkla = (2*M_PI/1.0)*(ma12*i2+ma22*i3+ma32*i4);   /*Z0311=367*/
            qzhkla = (2*M_PI/1.0)*(ma13*i2+ma23*i3+ma33*i4);   /*Z0311=368*/

            qxhklb = (2*M_PI/1.0)*(mb11*i2+mb21*i3+mb31*i4);   /*Z0311=370*/
            qyhklb = (2*M_PI/1.0)*(mb12*i2+mb22*i3+mb32*i4);   /*Z0311=371*/
            qzhklb = (2*M_PI/1.0)*(mb13*i2+mb23*i3+mb33*i4);
#endif

            dqxa = qx-qxhkla;   /*Z0311=374*/
            dqya = qy-qyhkla;   /*Z0311=375*/
            dqza = qz-qzhkla;   /*Z0311=376*/

            dqxb = qx-qxhklb;   /*Z0311=378*/
            dqyb = qy-qyhklb;   /*Z0311=379*/
            dqzb = qz-qzhklb;   /*Z0311=380*/

#ifdef PascalComment
            dqxa = qx-(r11a*qxn+r12a*qyn+r13a*qzn);   /*Z0311=382*/
            dqya = qy-(r21a*qxn+r22a*qyn+r23a*qzn);   /*Z0311=383*/
            dqza = qz-(r31a*qxn+r32a*qyn+r33a*qzn);   /*Z0311=384*/
            dqxb = qx-(r11b*qxn+r12b*qyn+r13b*qzn);   /*Z0311=385*/
            dqyb = qy-(r21b*qxn+r22b*qyn+r23b*qzn);   /*Z0311=386*/
            dqzb = qz-(r31b*qxn+r32b*qyn+r33b*qzn);
#endif

            ax1xa = (r11a*ax1x+r12a*ax1y+r13a*ax1z);   /*Z0311=389*/
            ax1ya = (r21a*ax1x+r22a*ax1y+r23a*ax1z);   /*Z0311=390*/
            ax1za = (r31a*ax1x+r32a*ax1y+r33a*ax1z);   /*Z0311=391*/
            ax1xb = (r11b*ax1x+r12b*ax1y+r13b*ax1z);   /*Z0311=392*/
            ax1yb = (r21b*ax1x+r22b*ax1y+r23b*ax1z);   /*Z0311=393*/
            ax1zb = (r31b*ax1x+r32b*ax1y+r33b*ax1z);   /*Z0311=394*/
            ax2xa = (r11a*ax2x+r12a*ax2y+r13a*ax2z);   /*Z0311=395*/
            ax2ya = (r21a*ax2x+r22a*ax2y+r23a*ax2z);   /*Z0311=396*/
            ax2za = (r31a*ax2x+r32a*ax2y+r33a*ax2z);   /*Z0311=397*/
            ax2xb = (r11b*ax2x+r12b*ax2y+r13b*ax2z);   /*Z0311=398*/
            ax2yb = (r21b*ax2x+r22b*ax2y+r23b*ax2z);   /*Z0311=399*/
            ax2zb = (r31b*ax2x+r32b*ax2y+r33b*ax2z);   /*Z0311=400*/
            ax3xa = (r11a*ax3x+r12a*ax3y+r13a*ax3z);   /*Z0311=401*/
            ax3ya = (r21a*ax3x+r22a*ax3y+r23a*ax3z);   /*Z0311=402*/
            ax3za = (r31a*ax3x+r32a*ax3y+r33a*ax3z);   /*Z0311=403*/
            ax3xb = (r11b*ax3x+r12b*ax3y+r13b*ax3z);   /*Z0311=404*/
            ax3yb = (r21b*ax3x+r22b*ax3y+r23b*ax3z);   /*Z0311=405*/
            ax3zb = (r31b*ax3x+r32b*ax3y+r33b*ax3z);   /*Z0311=406*/

            dqs1a = (dqxa*ax1xa+dqya*ax1ya+dqza*ax1za)/(ax1n*sigx);   /*Z0311=408*/
            dqs2a = (dqxa*ax2xa+dqya*ax2ya+dqza*ax2za)/(ax2n*sigy);   /*Z0311=409*/
            dqs3a = (dqxa*ax3xa+dqya*ax3ya+dqza*ax3za)/(ax3n*sigz);   /*Z0311=410*/
            dqs1b = (dqxb*ax1xb+dqyb*ax1yb+dqzb*ax1zb)/(ax1n*sigx);   /*Z0311=411*/
            dqs2b = (dqxb*ax2xb+dqyb*ax2yb+dqzb*ax2zb)/(ax2n*sigy);   /*Z0311=412*/
            dqs3b = (dqxb*ax3xb+dqyb*ax3yb+dqzb*ax3zb)/(ax3n*sigz);   /*Z0311=413*/

#ifdef PascalComment
            dqs1a = (dqxa*ax1x+dqya*ax1y+dqza*ax1z)/(ax1n*sigx);   /*Z0311=415*/
            dqs2a = (dqxa*ax2x+dqya*ax2y+dqza*ax2z)/(ax2n*sigy);   /*Z0311=416*/
            dqs3a = (dqxa*ax3x+dqya*ax3y+dqza*ax3z)/(ax3n*sigz);   /*Z0311=417*/
            dqs1b = (dqxb*ax1x+dqyb*ax1y+dqzb*ax1z)/(ax1n*sigx);   /*Z0311=418*/
            dqs2b = (dqxb*ax2x+dqyb*ax2y+dqzb*ax2z)/(ax2n*sigy);   /*Z0311=419*/
            dqs3b = (dqxb*ax3x+dqyb*ax3y+dqzb*ax3z)/(ax3n*sigz);
#endif

            arga = dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a;   /*Z0311=422*/
            argb = dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b;   /*Z0311=423*/
            pa = exp(-4*arga/M_PI);   /*Z0311=424*/
            pb = exp(-4*argb/M_PI);   /*Z0311=425*/
            /* if (arga>2) then pa:=eps */   /*Z0311=426*/
            /*    else pa:=exp(-4*arga/pi); */   /*Z0311=427*/
            /* if (argb>2) then pb:=eps */   /*Z0311=428*/
            /*    else pb:=exp(-4*argb/pi); */   /*Z0311=429*/
        }   /*Z0311=430*/
        if ( i1==8 )
        {    /* disk, general case */   /*Z0311=431*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=432*/
            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=433*/
            qxl = qx/qn;   /*Z0311=434*/
            qyl = qy/qn;   /*Z0311=435*/
            mla1 = p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);   /*Z0311=436*/
            mla2 = p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);   /*Z0311=437*/
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta); */   /*Z0311=438*/
            mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);   /*Z0311=439*/
            mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);   /*Z0311=440*/
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta); */   /*Z0311=441*/
            qnna = (qxl*mla1+qyl*mla2);   /*Z0311=442*/
            qnnb = (qxl*mlb1+qyl*mlb2);   /*Z0311=443*/
            arga = sqrt(1.0-qnna*qnna+eps9)*l*qn/(z+1);   /*Z0311=444*/
            argb = sqrt(1.0-qnnb*qnnb+eps9)*l*qn/(z+1);   /*Z0311=445*/

            if ( sigma<0.15 )
            {  /* series expansion/asymptote */   /*Z0311=447*/
                if ( arga<0.015 )
                {   /*Z0311=448*/
                    pa = 1;   /*Z0311=449*/
                    for ( i=1; i<=30; i++ ) pa = pa+carr1[i]*pow(arga/2.0,2*i);   /*Z0311=450*/
                    if ( (i3==1) ) pa = pa*pa;   /*Z0311=451*/
                }   /*Z0311=452*/
                else
                {   /*Z0311=453*/
                    if ( i3==0 )
                    { /* P(q) */   /*Z0311=454*/
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=455*/
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2*arga))/pow(1+4*arga*arga,(z-2)/2.0);   /*Z0311=456*/
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2*arga))/pow(1+4*arga*arga,(z-3)/2.0);   /*Z0311=457*/
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);   /*Z0311=458*/
                    }   /*Z0311=459*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=460*/
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1+arga*arga,(z-1/2.0)/2.0);   /*Z0311=461*/
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1+arga*arga,(z-3/2.0)/2.0);   /*Z0311=462*/
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);   /*Z0311=463*/
                        pa = pa3*pa3;   /*Z0311=464*/
                    }   /*Z0311=465*/
                }   /*Z0311=466*/
                if ( argb<0.015 )
                {   /*Z0311=467*/
                    pb = 1;   /*Z0311=468*/
                    for ( i=1; i<=30; i++ ) pb = pb+carr1[i]*pow(argb/2.0,2*i);   /*Z0311=469*/
                    if ( (i3==1) ) pb = pb*pb;   /*Z0311=470*/
                }   /*Z0311=471*/
                else
                {   /*Z0311=472*/
                    if ( i3==0 )
                    {  /* P(q) */   /*Z0311=473*/
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=474*/
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2*argb))/pow(1+4*argb*argb,(z-2)/2.0);   /*Z0311=475*/
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2*argb))/pow(1+4*argb*argb,(z-3)/2.0);   /*Z0311=476*/
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);   /*Z0311=477*/
                    }   /*Z0311=478*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=479*/
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1+argb*argb,(z-1/2.0)/2.0);   /*Z0311=480*/
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1+argb*argb,(z-3/2.0)/2.0);   /*Z0311=481*/
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);   /*Z0311=482*/
                        pb = pb3*pb3;   /*Z0311=483*/
                    }   /*Z0311=484*/
                }   /*Z0311=485*/
            }   /*Z0311=486*/
            else
            {  /* OZ-type */   /*Z0311=487*/
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=488*/
                pa = 1/(1+M_PI/(4*pa1));   /*Z0311=489*/
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=490*/
                pb = 1/(1+M_PI/(4*pb1));   /*Z0311=491*/
            }   /*Z0311=492*/
        }   /*Z0311=493*/
        /*Z0311=494*/
        if ( i1==9 )
        {   /* disk, x-axis */   /*Z0311=495*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=496*/
            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=497*/
            qxl = qx/qn;   /*Z0311=498*/
            qyl = qy/qn;   /*Z0311=499*/
            qnna = qxl*cos(delta)-qyl*sin(a)*sin(delta);   /*Z0311=500*/
            qnnb = qxl*cos(delta)-qyl*sin(b)*sin(delta);   /*Z0311=501*/
            arga = sqrt(1.0-qnna*qnna+eps9)*l*qn/(z+1);   /*Z0311=502*/
            argb = sqrt(1.0-qnnb*qnnb+eps9)*l*qn/(z+1);   /*Z0311=503*/
            /*Z0311=504*/
            if ( sigma<0.15 )
            {  /* series expansion/asymptote */   /*Z0311=505*/
                if ( arga<0.015 )
                {   /*Z0311=506*/
                    pa = 1;   /*Z0311=507*/
                    for ( i=1; i<=30; i++ ) pa = pa+carr1[i]*pow(arga/2.0,2*i);   /*Z0311=508*/
                    if ( i3==1 ) pa = pa*pa;   /*Z0311=509*/
                }   /*Z0311=510*/
                else
                {   /*Z0311=511*/
                    if ( i3==0 )
                    { /* P(q) */   /*Z0311=512*/
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=513*/
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2*arga))/pow(1+4*arga*arga,(z-2)/2.0);   /*Z0311=514*/
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2*arga))/pow(1+4*arga*arga,(z-3)/2.0);   /*Z0311=515*/
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);   /*Z0311=516*/
                    }   /*Z0311=517*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=518*/
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1+arga*arga,(z-1/2.0)/2.0);   /*Z0311=519*/
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1+arga*arga,(z-3/2.0)/2.0);   /*Z0311=520*/
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);   /*Z0311=521*/
                        pa = pa3*pa3;   /*Z0311=522*/
                    }   /*Z0311=523*/
                }   /*Z0311=524*/
                if ( argb<0.015 )
                {   /*Z0311=525*/
                    pb = 1;   /*Z0311=526*/
                    for ( i=1; i<=30; i++ ) pb = pb+carr1[i]*pow(argb/2.0,2*i);   /*Z0311=527*/
                    if ( i3==1 ) pb = pb*pb;   /*Z0311=528*/
                }   /*Z0311=529*/
                else
                {   /*Z0311=530*/
                    if ( i3==0 )
                    {  /* P(q) */   /*Z0311=531*/
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=532*/
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2*argb))/pow(1+4*argb*argb,(z-2)/2.0);   /*Z0311=533*/
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2*argb))/pow(1+4*argb*argb,(z-3)/2.0);   /*Z0311=534*/
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);   /*Z0311=535*/
                    }   /*Z0311=536*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=537*/
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1+argb*argb,(z-1/2.0)/2.0);   /*Z0311=538*/
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1+argb*argb,(z-3/2.0)/2.0);   /*Z0311=539*/
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);   /*Z0311=540*/
                        pb = pb3*pb3;   /*Z0311=541*/
                    }   /*Z0311=542*/
                }   /*Z0311=543*/
            }   /*Z0311=544*/
            else
            {  /* OZ-type */   /*Z0311=545*/
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=546*/
                pa = 1/(1+M_PI/(4*pa1));   /*Z0311=547*/
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=548*/
                pb = 1/(1+M_PI/(4*pb1));   /*Z0311=549*/
            }   /*Z0311=550*/
        }   /*Z0311=551*/
        /*Z0311=552*/
        if ( i1==10 )
        {   /* disk, y-axis */   /*Z0311=553*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=554*/
            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=555*/
            qxl = qx/qn;   /*Z0311=556*/
            qyl = qy/qn;   /*Z0311=557*/
            qnna = qxl*sin(a)*sin(delta)+qyl*cos(delta);   /*Z0311=558*/
            qnnb = qxl*sin(b)*sin(delta)+qyl*cos(delta);   /*Z0311=559*/
            arga = sqrt(1.0-qnna*qnna+eps9)*l*qn/(z+1);   /*Z0311=560*/
            argb = sqrt(1.0-qnnb*qnnb+eps9)*l*qn/(z+1);   /*Z0311=561*/
            /*Z0311=562*/
            if ( sigma<0.15 )
            {  /* series expansion/asymptote */   /*Z0311=563*/
                if ( arga<0.015 )
                {   /*Z0311=564*/
                    pa = 1;   /*Z0311=565*/
                    for ( i=1; i<=30; i++ ) pa = pa+carr1[i]*pow(arga/2.0,2*i);   /*Z0311=566*/
                    if ( i3==1 ) pa = pa*pa;   /*Z0311=567*/
                }   /*Z0311=568*/
                else
                {   /*Z0311=569*/
                    if ( i3==0 )
                    { /* P(q) */   /*Z0311=570*/
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=571*/
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2*arga))/pow(1+4*arga*arga,(z-2)/2.0);   /*Z0311=572*/
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2*arga))/pow(1+4*arga*arga,(z-3)/2.0);   /*Z0311=573*/
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);   /*Z0311=574*/
                    }   /*Z0311=575*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=576*/
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1+arga*arga,(z-1/2.0)/2.0);   /*Z0311=577*/
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1+arga*arga,(z-3/2.0)/2.0);   /*Z0311=578*/
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);   /*Z0311=579*/
                        pa = pa3*pa3;   /*Z0311=580*/
                    }   /*Z0311=581*/
                }   /*Z0311=582*/
                if ( argb<0.015 )
                {   /*Z0311=583*/
                    pb = 1;   /*Z0311=584*/
                    for ( i=1; i<=30; i++ ) pb = pb+carr1[i]*pow(argb/2.0,2*i);   /*Z0311=585*/
                    if ( i3==1 ) pb = pb*pb;   /*Z0311=586*/
                }   /*Z0311=587*/
                else
                {   /*Z0311=588*/
                    if ( i3==0 )
                    {  /* P(q) */   /*Z0311=589*/
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=590*/
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2*argb))/pow(1+4*argb*argb,(z-2)/2.0);   /*Z0311=591*/
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2*argb))/pow(1+4*argb*argb,(z-3)/2.0);   /*Z0311=592*/
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);   /*Z0311=593*/
                    }   /*Z0311=594*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=595*/
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1+argb*argb,(z-1/2.0)/2.0);   /*Z0311=596*/
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1+argb*argb,(z-3/2.0)/2.0);   /*Z0311=597*/
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);   /*Z0311=598*/
                        pb = pb3*pb3;   /*Z0311=599*/
                    }   /*Z0311=600*/
                }   /*Z0311=601*/
            }   /*Z0311=602*/
            else
            {  /* OZ-type */   /*Z0311=603*/
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=604*/
                pa = 1/(1+M_PI/(4*pa1));   /*Z0311=605*/
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=606*/
                pb = 1/(1+M_PI/(4*pb1));   /*Z0311=607*/
            }   /*Z0311=608*/
        }   /*Z0311=609*/
        /*Z0311=610*/
        if ( i1==11 )
        {   /* disk, -z-axis */   /*Z0311=611*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=612*/
            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=613*/
            qxl = qx/qn;   /*Z0311=614*/
            qyl = qy/qn;   /*Z0311=615*/
            qnna = qxl*sin(a)*sin(delta)-qyl*cos(a)*sin(delta);   /*Z0311=616*/
            qnnb = qxl*sin(b)*sin(delta)-qyl*cos(b)*sin(delta);   /*Z0311=617*/
            arga = sqrt(1.0-qnna*qnna+eps9)*l*qn/(z+1);   /*Z0311=618*/
            argb = sqrt(1.0-qnnb*qnnb+eps9)*l*qn/(z+1);   /*Z0311=619*/
            /*Z0311=620*/
            if ( sigma<0.15 )
            {  /* series expansion/asymptote */   /*Z0311=621*/
                if ( arga<0.015 )
                {   /*Z0311=622*/
                    pa = 1;   /*Z0311=623*/
                    for ( i=1; i<=30; i++ ) pa = pa+carr1[i]*pow(arga/2.0,2*i);   /*Z0311=624*/
                    if ( i3==1 ) pa = pa*pa;   /*Z0311=625*/
                }   /*Z0311=626*/
                else
                {   /*Z0311=627*/
                    if ( i3==0 )
                    { /* P(q) */   /*Z0311=628*/
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=629*/
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2*arga))/pow(1+4*arga*arga,(z-2)/2.0);   /*Z0311=630*/
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2*arga))/pow(1+4*arga*arga,(z-3)/2.0);   /*Z0311=631*/
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);   /*Z0311=632*/
                    }   /*Z0311=633*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=634*/
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1+arga*arga,(z-1/2.0)/2.0);   /*Z0311=635*/
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1+arga*arga,(z-3/2.0)/2.0);   /*Z0311=636*/
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);   /*Z0311=637*/
                        pa = pa3*pa3;   /*Z0311=638*/
                    }   /*Z0311=639*/
                }   /*Z0311=640*/
                if ( argb<0.015 )
                {   /*Z0311=641*/
                    pb = 1;   /*Z0311=642*/
                    for ( i=1; i<=30; i++ ) pb = pb+carr1[i]*pow(argb/2.0,2*i);   /*Z0311=643*/
                    if ( i3==1 ) pb = pb*pb;   /*Z0311=644*/
                }   /*Z0311=645*/
                else
                {   /*Z0311=646*/
                    if ( i3==0 )
                    {  /* P(q) */   /*Z0311=647*/
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=648*/
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2*argb))/pow(1+4*argb*argb,(z-2)/2.0);   /*Z0311=649*/
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2*argb))/pow(1+4*argb*argb,(z-3)/2.0);   /*Z0311=650*/
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);   /*Z0311=651*/
                    }   /*Z0311=652*/
                    if ( i3==1 )
                    {  /* F(q) */   /*Z0311=653*/
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1+argb*argb,(z-1/2.0)/2.0);   /*Z0311=654*/
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1+argb*argb,(z-3/2.0)/2.0);   /*Z0311=655*/
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);   /*Z0311=656*/
                        pb = pb3*pb3;   /*Z0311=657*/
                    }   /*Z0311=658*/
                }   /*Z0311=659*/
            }   /*Z0311=660*/
            else
            {  /* OZ-type */   /*Z0311=661*/
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);   /*Z0311=662*/
                pa = 1/(1+M_PI/(4*pa1));   /*Z0311=663*/
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);   /*Z0311=664*/
                pb = 1/(1+M_PI/(4*pb1));   /*Z0311=665*/
            }   /*Z0311=666*/
        }   /*Z0311=667*/
        /*Z0311=668*/
        if ( i1==12 )
        {   /* isotropic cube */   /*Z0311=669*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=670*/

#ifdef PascalComment
            arga1 = qn*sin(delta)*cos(a)*r+eps;   /*Z0311=672*/
            arga2 = qn*sin(delta)*sin(a)*r+eps;   /*Z0311=673*/
            arga3 = qn*cos(delta)*r+eps;   /*Z0311=674*/
            if ( (i3==0) )
            { /* P(q) */   /*Z0311=675*/
                if ( (arga1 < 0.001) ) pa1 = 1 else pa1 = szave(5,1,arga1,-2,z);   /*Z0311=676*/
                if ( (arga2 < 0.001) ) pa2 = 1 else pa2 = szave(5,1,arga2,-2,z);   /*Z0311=677*/
                if ( (arga3 < 0.001) ) pa3 = 1 else pa3 = szave(5,1,arga3,-2,z);   /*Z0311=678*/
                pa = pa1*pa2*pa3;   /*Z0311=679*/
            }   /*Z0311=680*/
            if ( (i3==1) )
            { /* F(q) */   /*Z0311=681*/
                if ( (arga1 < 0.001) ) pa1 = 1 else pa1 = szave(3,1,arga1,-1,z);   /*Z0311=682*/
                if ( (arga2 < 0.001) ) pa2 = 1 else pa2 = szave(3,1,arga2,-1,z);   /*Z0311=683*/
                if ( (arga3 < 0.001) ) pa3 = 1 else pa3 = szave(3,1,arga3,-1,z);   /*Z0311=684*/
                pa = pa1*pa1*pa2*pa2*pa3*pa3;   /*Z0311=685*/
            }   /*Z0311=686*/
            argb1 = qn*sin(delta)*cos(b)*r+eps;   /*Z0311=687*/
            argb2 = qn*sin(delta)*sin(b)*r+eps;   /*Z0311=688*/
            argb3 = qn*cos(delta)*r+eps;   /*Z0311=689*/
            if ( (i3==0) )
            { /* P(q) */   /*Z0311=690*/
                if ( (argb1 < 0.001) ) pb1 = 1 else pb1 = szave(5,1,argb1,-2,z);   /*Z0311=691*/
                if ( (argb2 < 0.001) ) pb2 = 1 else pb2 = szave(5,1,argb2,-2,z);   /*Z0311=692*/
                if ( (argb3 < 0.001) ) pb3 = 1 else pb3 = szave(5,1,argb3,-2,z);   /*Z0311=693*/
                pb = pb1*pb2*pb3;   /*Z0311=694*/
            }   /*Z0311=695*/
            if ( (i3==1) )
            { /* F(q) */   /*Z0311=696*/
                if ( (argb1 < 0.001) ) pb1 = 1 else pb1 = szave(3,1,argb1,-1,z);   /*Z0311=697*/
                if ( (argb2 < 0.001) ) pb2 = 1 else pb2 = szave(3,1,argb2,-1,z);   /*Z0311=698*/
                if ( (argb3 < 0.001) ) pb3 = 1 else pb3 = szave(3,1,argb3,-1,z);   /*Z0311=699*/
                pb = pb1*pb1*pb2*pb2*pb3*pb3;   /*Z0311=700*/
            }
#endif

            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=703*/
            arga1 = qn*sin(delta)*cos(a)*r/(z+1)+eps9;   /*Z0311=704*/
            arga2 = qn*sin(delta)*sin(a)*r/(z+1)+eps9;   /*Z0311=705*/
            arga3 = qn*cos(delta)*r/(z+1)+eps9;   /*Z0311=706*/
            if ( i3==0 )
            {  /* P(q) */   /*Z0311=707*/
                a1 = 1/(2*z*(z-1));   /*Z0311=708*/
                if ( arga1 < 0.001 )
                    pa1 = 1;
                else   /*Z0311=709*/
                    pa1 = (a1/(arga1*arga1+eps9))*(1-cos((z-1)*atan(2*arga1))/pow(1+4*arga1*arga1,(z-1)/2.0));   /*Z0311=710*/
                if ( arga2 < 0.001 )
                    pa2 = 1;
                else   /*Z0311=711*/
                    pa2 = (a1/(arga2*arga2+eps9))*(1-cos((z-1)*atan(2*arga2))/pow(1+4*arga2*arga2,(z-1)/2.0));   /*Z0311=712*/
                if ( arga3 < 0.001 )
                    pa3 = 1;
                else   /*Z0311=713*/
                    pa3 = (a1/(arga3*arga3+eps9))*(1-cos((z-1)*atan(2*arga3))/pow(1+4*arga3*arga3,(z-1)/2.0));   /*Z0311=714*/
                pa = pa1*pa2*pa3;   /*Z0311=715*/
            }   /*Z0311=716*/
            if ( i3==1 )
            { /* F(q) */   /*Z0311=717*/
                a1 = 1/z;   /*Z0311=718*/
                if ( arga1 < 0.001 )
                    pa1 = 1;
                else   /*Z0311=719*/
                    pa1 = (a1/(arga1+eps9))*sin(z*atan(arga1))/pow(1+arga1*arga1,z/2.0);   /*Z0311=720*/
                if ( arga2 < 0.001 )
                    pa2 = 1;
                else   /*Z0311=721*/
                    pa2 = (a1/(arga2+eps9))*sin(z*atan(arga2))/pow(1+arga2*arga2,z/2.0);   /*Z0311=722*/
                if ( arga3 < 0.001 )
                    pa3 = 1;
                else   /*Z0311=723*/
                    pa3 = (a1/(arga3+eps9))*sin(z*atan(arga3))/pow(1+arga3*arga3,z/2.0);   /*Z0311=724*/
                pa = pa1*pa1*pa2*pa2*pa3*pa3;   /*Z0311=725*/
            }   /*Z0311=726*/
            argb1 = qn*sin(delta)*cos(b)*r/(z+1);   /*Z0311=727*/
            argb2 = qn*sin(delta)*sin(b)*r/(z+1);   /*Z0311=728*/
            argb3 = qn*cos(delta)*r/(z+1);   /*Z0311=729*/
            if ( i3==0 )
            {   /* P(q) */   /*Z0311=730*/
                a1 = 1/(2*z*(z-1));   /*Z0311=731*/
                if ( argb1 < 0.001 )
                    pb1 = 1;
                else   /*Z0311=732*/
                    pb1 = (a1/(argb1*argb1+eps9))*(1-cos((z-1)*atan(2*argb1))/pow(1+4*argb1*argb1,(z-1)/2.0));   /*Z0311=733*/
                if ( argb2 < 0.001 )
                    pb2 = 1;
                else   /*Z0311=734*/
                    pb2 = (a1/(argb2*argb2+eps9))*(1-cos((z-1)*atan(2*argb2))/pow(1+4*argb2*argb2,(z-1)/2.0));   /*Z0311=735*/
                if ( argb3 < 0.001 )
                    pb3 = 1;
                else   /*Z0311=736*/
                    pb3 = (a1/(argb3*argb3+eps9))*(1-cos((z-1)*atan(2*argb3))/pow(1+4*argb3*argb3,(z-1)/2.0));   /*Z0311=737*/
                pb = pb1*pb2*pb3;   /*Z0311=738*/
            }   /*Z0311=739*/
            if ( i3==1 )
            { /* F(q) */   /*Z0311=740*/
                a1 = 1/z;   /*Z0311=741*/
                if ( argb1 < 0.001 )
                    pb1 = 1;
                else   /*Z0311=742*/
                    pb1 = (a1/(argb1+eps9))*sin(z*atan(argb1))/pow(1+argb1*argb1,z/2.0);   /*Z0311=743*/
                if ( argb2 < 0.001 )
                    pb2 = 1;
                else   /*Z0311=744*/
                    pb2 = (a1/(argb2+eps9))*sin(z*atan(argb2))/pow(1+argb2*argb2,z/2.0);   /*Z0311=745*/
                if ( arga3 < 0.001 )
                    pb3 = 1;
                else   /*Z0311=746*/
                    pb3 = (a1/(argb3+eps9))*sin(z*atan(argb3))/pow(1+argb3*argb3,z/2.0);   /*Z0311=747*/
                pb = pb1*pb1*pb2*pb2*pb3*pb3;   /*Z0311=748*/
            }   /*Z0311=749*/
        }   /*Z0311=750*/
        if ( i1==13 )
        {   /* biaxial ellipsoid, isotropic */   /*Z0311=751*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=752*/
            epsi = l/r;   /*Z0311=753*/
            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=754*/
            arga = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(a)*cos(a))/(z+1);   /*Z0311=755*/
            a1 = (1/(2*z*(z-1)*(z-2)*(z-3)));   /*Z0311=756*/
            pa1 = a1*pow(arga,-4)*(1+cos((z-3)*atan(2*arga))/pow(1+4*arga*arga,(z-3)/2.0));   /*Z0311=757*/
            pa2 = (a1/(z-4))*pow(arga,-5)*sin((z-4)*atan(2*arga))/pow(1+4*arga*arga,(z-4)/2.0);   /*Z0311=758*/
            pa3 = (a1/((z-4)*(z-5)))*pow(arga,-6)*(1-cos((z-5)*atan(2*arga))/pow(1+4*arga*arga,(z-5)/2.0));   /*Z0311=759*/
            pa = 9*(pa1-2*pa2+pa3)*sin(a);   /*Z0311=760*/
            argb = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(b)*cos(b))/(z+1);   /*Z0311=761*/
            pb1 = a1*pow(argb,-4)*(1+cos((z-3)*atan(2*argb))/pow(1+4*argb*argb,(z-3)/2.0));   /*Z0311=762*/
            pb2 = (a1/(z-4))*pow(argb,-5)*sin((z-4)*atan(2*argb))/pow(1+4*argb*argb,(z-4)/2.0);   /*Z0311=763*/
            pb3 = (a1/((z-4)*(z-5)))*pow(argb,-6)*(1-cos((z-5)*atan(2*argb))/pow(1+4*argb*argb,(z-5)/2.0));   /*Z0311=764*/
            pb = 9*(pb1-2*pb2+pb3)*sin(b);   /*Z0311=765*/
        }  /* of biaxial ellipsoid */   /*Z0311=766*/

        if ( i1==14 )
        {   /* triaxial ellipsoid, isotropic */   /*Z0311=768*/
            ella = r;   /*Z0311=769*/
            ellb = l;   /*Z0311=770*/
            ellc = r/params.p1;   /*Z0311=771*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=772*/
            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=773*/
            arga = qn*sqrt(pow(ella*cos(a)*sin(delta),2)+pow(ellb*sin(a)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);   /*Z0311=774*/
            a1 = (1/(2*z*(z-1)*(z-2)*(z-3)));   /*Z0311=775*/
            pa1 = a1*pow(arga,-4)*(1+cos((z-3)*atan(2*arga))/pow(1+4*arga*arga,(z-3)/2.0));   /*Z0311=776*/
            pa2 = (a1/(z-4))*pow(arga,-5)*sin((z-4)*atan(2*arga))/pow(1+4*arga*arga,(z-4)/2.0);   /*Z0311=777*/
            pa3 = (a1/((z-4)*(z-5)))*pow(arga,-6)*(1-cos((z-5)*atan(2*arga))/pow(1+4*arga*arga,(z-5)/2.0));   /*Z0311=778*/
            pa = 9*(pa1-2*pa2+pa3);   /*Z0311=779*/
            /* pa:=power(3*(sin(qn*r)-qn*r*cos(qn*r))/(qn*qn*qn*r*r*r+eps),2); */   /*Z0311=780*/
            /* pa:=1.05; */   /*Z0311=781*/
            argb = qn*sqrt(pow(ella*cos(b)*sin(delta),2)+pow(ellb*sin(b)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);   /*Z0311=782*/
            pb1 = a1*pow(argb,-4)*(1+cos((z-3)*atan(2*argb))/pow(1+4*argb*argb,(z-3)/2.0));   /*Z0311=783*/
            pb2 = (a1/(z-4))*pow(argb,-5)*sin((z-4)*atan(2*argb))/pow(1+4*argb*argb,(z-4)/2.0);   /*Z0311=784*/
            pb3 = (a1/((z-4)*(z-5)))*pow(argb,-6)*(1-cos((z-5)*atan(2*argb))/pow(1+4*argb*argb,(z-5)/2.0));   /*Z0311=785*/
            pb = 9*(pb1-2*pb2+pb3);   /*Z0311=786*/
            /* pb:=power(3*(sin(qn*r*1.04)-qn*r*cos(qn*r))/(qn*qn*qn*r*r*r+eps),2); */   /*Z0311=787*/
            /* pb:=1.04; */   /*Z0311=788*/
        }  /* of ellipsoid */   /*Z0311=789*/

        if ( i1==15 )
        {   /* barrel area integration */   /*Z0311=791*/
            pa1 = r*pow(1-pow(a/l,params.p1),1/params.p1);   /*Z0311=792*/
            arga = -r*pow(a/l,params.p1-1)*pow(1-pow(a/l,params.p1),(1/params.p1)-1)/l;   /*Z0311=793*/
            pa = pa1*sqrt(1.0+arga*arga);   /*Z0311=794*/
            pb1 = r*pow(1-pow(b/l,params.p1),1/params.p1);   /*Z0311=795*/
            argb = -r*pow(b/l,params.p1-1)*pow(1-pow(b/l,params.p1),(1/params.p1)-1)/l;   /*Z0311=796*/
            pb = pb1*sqrt(1.0+argb*argb);   /*Z0311=797*/
        }   /*Z0311=798*/

        if ( i1==16 )
        {   /* barrel, x-axis */   /*Z0311=800*/
            z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=801*/
            qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=802*/
            arga = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(a)*sin(delta)+eps9));   /*Z0311=803*/
            argb = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(b)*sin(delta)+eps9));   /*Z0311=804*/
            if ( i3==0 )
            { /* P(q) */   /*Z0311=805*/
                a1 = 9*pow(z+1,4)/(2*z*(z-1)*(z-2)*(z-3));   /*Z0311=806*/
                pa = (a1/(arga*arga));   /*Z0311=807*/
                pb = (a1/(argb*argb));   /*Z0311=808*/
            }   /*Z0311=809*/
            if ( i3==1 )
            { /* F(q) */   /*Z0311=810*/
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1+arga*arga,z/2.0);   /*Z0311=811*/
                pa = pa1*pa1;   /*Z0311=812*/
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1+argb*argb,z/2.0);   /*Z0311=813*/
                pb = pb1*pb1;   /*Z0311=814*/
            }   /*Z0311=815*/
        }   /*Z0311=816*/
        /*Z0311=817*/
        pq = 0.5*(b-a)*(pa+pb);   /*Z0311=818*/
        trapzdchid_cnt = 1;   /*Z0311=819*/
    }   /*Z0311=820*/
    else
    {   /*Z0311=821*/
        tnm = trapzdchid_cnt;   /*Z0311=822*/
        del = (b-a)/tnm;   /*Z0311=823*/
        x = a+0.5*del;   /*Z0311=824*/
        sump = 0.0;   /*Z0311=825*/
        for ( j=1; j<=trapzdchid_cnt; j++ )
        {   /*Z0311=826*/
            if ( i1==1 )
            {  /* cylinder, general case */   /*Z0311=827*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=828*/
                mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);   /*Z0311=829*/
                mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);   /*Z0311=830*/
                /* mlx3:=p31*cos(x)*sin(delta)+p33*cos(delta); */   /*Z0311=831*/
                argx = (qx*mlx1+qy*mlx2+eps9)*l/(z+1);   /*Z0311=832*/
                if ( i3==0 )
                { /* P(q) */   /*Z0311=833*/
                    a1 = 1/(2*z*(z-1));   /*Z0311=834*/
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.0));   /*Z0311=835*/
                }   /*Z0311=836*/
                if ( i3==1 )
                { /* F(q) */   /*Z0311=837*/
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1+argx*argx,z/2.0);   /*Z0311=838*/
                    px = px1*px1;   /*Z0311=839*/
                }   /*Z0311=840*/
            }   /*Z0311=841*/
            if ( i1==2 )
            {   /* cylinder, x-axis */   /*Z0311=842*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=843*/
                argx = (qx*cos(delta)-qy*sin(x)*sin(delta)+eps9)*l/(z+1);   /*Z0311=844*/
                if ( i3==0 )
                { /* P(q) */   /*Z0311=845*/
                    a1 = 1/(2*z*(z-1));   /*Z0311=846*/
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.0));   /*Z0311=847*/
                }   /*Z0311=848*/
                if ( i3==1 )
                { /* F(q) */   /*Z0311=849*/
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1+argx*argx,z/2.0);   /*Z0311=850*/
                    px = px1*px1;   /*Z0311=851*/
                }   /*Z0311=852*/
            }   /*Z0311=853*/
            if ( i1==3 )
            {   /* cylinder, y-axis */   /*Z0311=854*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=855*/
                a1 = 1/(2*z*(z-1));   /*Z0311=856*/
                argx = (qx*sin(x)*sin(delta)+qy*cos(delta)+eps9)*l/(z+1);   /*Z0311=857*/
                if ( i3==0 )
                { /* P(q) */   /*Z0311=858*/
                    a1 = 1/(2*z*(z-1));   /*Z0311=859*/
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.0));   /*Z0311=860*/
                }   /*Z0311=861*/
                if ( i3==1 )
                { /* F(q) */   /*Z0311=862*/
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1+argx*argx,z/2.0);   /*Z0311=863*/
                    px = px1*px1;   /*Z0311=864*/
                }   /*Z0311=865*/
            }   /*Z0311=866*/
            if ( i1==4 )
            {   /* cylinder, -z-axis */   /*Z0311=867*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=868*/
                a1 = 1/(2*z*(z-1));   /*Z0311=869*/
                argx = (qx*sin(x)*sin(delta)-qy*cos(x)*sin(delta)+eps9)*l/(z+1);   /*Z0311=870*/
                if ( i3==0 )
                { /* P(q) */   /*Z0311=871*/
                    a1 = 1/(2*z*(z-1));   /*Z0311=872*/
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.0));   /*Z0311=873*/
                }   /*Z0311=874*/
                if ( i3==1 )
                { /* F(q) */   /*Z0311=875*/
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1+argx*argx,z/2.0);   /*Z0311=876*/
                    px = px1*px1;   /*Z0311=877*/
                }   /*Z0311=878*/
            }   /*Z0311=879*/
            if ( i1==5 )
            {   /* general series expansion */   /*Z0311=880*/
                mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);   /*Z0311=881*/
                mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);   /*Z0311=882*/
                /* mlx3 = p31*cos(x)*sin(delta)+p33*cos(delta); */   /*Z0311=883*/
                px = pow(mlx1,i3)*pow(mlx2,i4);   /*Z0311=884*/
            }   /*Z0311=885*/
            if ( i1==6 )
            {    /* unit cell rotation */   /*Z0311=886*/
                dqxx = qx-qhkl*(p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta));   /*Z0311=887*/
                dqyx = qy-qhkl*(p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta));   /*Z0311=888*/
                dqzx = qz-qhkl*(p31*cos(x)*sin(delta)+p33*cos(delta));   /*Z0311=889*/
                dqs1x = (dqxx*ax1x+dqyx*ax1y+dqzx*ax1z)/(ax1n*sigx);   /*Z0311=890*/
                dqs2x = (dqxx*ax2x+dqyx*ax2y+dqzx*ax2z)/(ax2n*sigy);   /*Z0311=891*/
                dqs3x = (dqxx*ax3x+dqyx*ax3y+dqzx*ax3z)/(ax3n*sigz);   /*Z0311=892*/
                px = exp(-4*(dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x)/M_PI);   /*Z0311=893*/
            }   /*Z0311=894*/
            if ( i1==7 )
            {    /* fiber unit cell rotation */   /*Z0311=895*/
                /* l1:=sin(theta)*cos(phi); */   /*Z0311=896*/
                /* l2:=sin(theta)*sin(phi); */   /*Z0311=897*/
                /* l3:=-cos(theta); */   /*Z0311=898*/
                l1 = r;   /*Z0311=899*/
                l2 = sigma;   /*Z0311=900*/
                l3 = dbeta;   /*Z0311=901*/

                r11x = cos(x)+(1-cos(x))*l1*l1;   /*Z0311=903*/
                r12x = -l3*sin(x)+(1-cos(x))*l1*l2;   /*Z0311=904*/
                r13x = -l2*sin(x)+(1-cos(x))*l1*l3;   /*Z0311=905*/
                r21x = l3*sin(x)+(1-cos(x))*l1*l2;   /*Z0311=906*/
                r22x = cos(x)+(1-cos(x))*l2*l2;   /*Z0311=907*/
                r23x = l1*sin(x)+(1-cos(x))*l2*l3;   /*Z0311=908*/
                r31x = l2*sin(x)+(1-cos(x))*l1*l3;   /*Z0311=909*/
                r32x = -l1*sin(x)+(1-cos(x))*l2*l3;   /*Z0311=910*/
                r33x = cos(x)+(1-cos(x))*l3*l3;   /*Z0311=911*/

#ifdef PascalComment
                /* scattering vector hkl */   /*Z0311=913*/
                qxhkl = (2*M_PI/1.0)*(p11*i2+p21*i3+p31*i4);   /*Z0311=914*/
                qyhkl = (2*M_PI/1.0)*(p12*i2+p22*i3+p32*i4);   /*Z0311=915*/
                qzhkl = (2*M_PI/1.0)*(p13*i2+p23*i3+p33*i4);
#endif

                /* rotate this scattering vector */   /*Z0311=918*/
                qxhklx = r11x*qxn+r12x*qyn+r13x*qzn;   /*Z0311=919*/
                qyhklx = r21x*qxn+r22x*qyn+r23x*qzn;   /*Z0311=920*/
                qzhklx = r31x*qxn+r32x*qyn+r33x*qzn;   /*Z0311=921*/

#ifdef PascalComment
                11:=0.707;   /*Z0311=923*/
                p12 = 0.707;   /*Z0311=924*/
                p13 = 0;   /*Z0311=925*/
                p21 = 0.408;   /*Z0311=926*/
                p22 = -0.408;   /*Z0311=927*/
                p23 = 0.816;   /*Z0311=928*/
                p31 = 0.577;   /*Z0311=929*/
                p32 = -0.577;   /*Z0311=930*/
                p33 = -0.577;   /*Z0311=931*/

                /* rotated a */   /*Z0311=933*/
                aexx = r11x*p11+r12x*p12+r13x*p13;   /*Z0311=934*/
                aeyx = r21x*p11+r22x*p12+r23x*p13;   /*Z0311=935*/
                aezx = r31x*p11+r32x*p12+r33x*p13;   /*Z0311=936*/

                /* rotated b */   /*Z0311=938*/
                bexx = r11x*p21+r12x*p22+r13x*p23;   /*Z0311=939*/
                beyx = r21x*p21+r22x*p22+r23x*p23;   /*Z0311=940*/
                bezx = r31x*p21+r32x*p22+r33x*p23;   /*Z0311=941*/

                /* rotated c */   /*Z0311=943*/
                cexx = r11x*p31+r12x*p32+r13x*p33;   /*Z0311=944*/
                ceyx = r21x*p31+r22x*p32+r23x*p33;   /*Z0311=945*/
                cezx = r31x*p31+r32x*p32+r33x*p33;   /*Z0311=946*/

                mx11 = aexx;   /*Z0311=948*/
                mx12 = aeyx;   /*Z0311=949*/
                mx13 = aezx;   /*Z0311=950*/
                mx21 = bexx;   /*Z0311=951*/
                mx22 = beyx;   /*Z0311=952*/
                mx23 = bezx;   /*Z0311=953*/
                mx31 = cexx;   /*Z0311=954*/
                mx32 = ceyx;   /*Z0311=955*/
                mx33 = cezx;   /*Z0311=956*/

                /* reciprocal space base vectors */   /*Z0311=958*/
                volx = aexx*(beyx*cezx-bezx*ceyx)+aeyx*(bezx*cexx-bexx*cezx)+aezx*(bexx*ceyx-beyx*cexx);   /*Z0311=959*/
                mx11 = (beyx*cezx-bezx*ceyx)/volx;   /*Z0311=960*/
                mx12 = (bezx*cexx-bexx*cezx)/volx;   /*Z0311=961*/
                mx13 = (bexx*ceyx-beyx*cexx)/volx;   /*Z0311=962*/
                mx21 = (aezx*ceyx-aeyx*cezx)/volx;   /*Z0311=963*/
                mx22 = (aexx*cezx-aezx*cexx)/volx;   /*Z0311=964*/
                mx23 = (aeyx*cexx-aexx*ceyx)/volx;   /*Z0311=965*/
                mx31 = (aeyx*bezx-aezx*beyx)/volx;   /*Z0311=966*/
                mx32 = (aezx*bexx-aexx*bezx)/volx;   /*Z0311=967*/
                mx33 = (aexx*beyx-aeyx*bexx)/volx;   /*Z0311=968*/

                qxhklx = (2*M_PI/1.0)*(mx11*i2+mx21*i3+mx31*i4);   /*Z0311=970*/
                qyhklx = (2*M_PI/1.0)*(mx12*i2+mx22*i3+mx32*i4);   /*Z0311=971*/
                qzhklx = (2*M_PI/1.0)*(mx13*i2+mx23*i3+mx33*i4);   /*Z0311=972*/

                r11x = cos(x);   /*Z0311=974*/
                r12x = -sin(x);   /*Z0311=975*/
                r13x = 0;   /*Z0311=976*/
                r21x = sin(x);   /*Z0311=977*/
                r22x = cos(x);   /*Z0311=978*/
                r23x = 0;   /*Z0311=979*/
                r31x = 0;   /*Z0311=980*/
                r32x = 0;   /*Z0311=981*/
                r33x = 1;   /*Z0311=982*/

                dqxx = qx-(r11x*qxn+r12x*qyn+r13x*qzn);   /*Z0311=984*/
                dqyx = qy-(r21x*qxn+r22x*qyn+r23x*qzn);   /*Z0311=985*/
                dqzx = qz-(r31x*qxn+r32x*qyn+r33x*qzn);
#endif

                dqxx = qx-qxhklx;   /*Z0311=988*/
                dqyx = qy-qyhklx;   /*Z0311=989*/
                dqzx = qz-qzhklx;   /*Z0311=990*/

                ax1xx = (r11x*ax1x+r12x*ax1y+r13x*ax1z);   /*Z0311=992*/
                ax1yx = (r21x*ax1x+r22x*ax1y+r23x*ax1z);   /*Z0311=993*/
                ax1zx = (r31x*ax1x+r32x*ax1y+r33x*ax1z);   /*Z0311=994*/
                ax2xx = (r11x*ax2x+r12x*ax2y+r13x*ax2z);   /*Z0311=995*/
                ax2yx = (r21x*ax2x+r22x*ax2y+r23x*ax2z);   /*Z0311=996*/
                ax2zx = (r31x*ax2x+r32x*ax2y+r33x*ax2z);   /*Z0311=997*/
                ax3xx = (r11x*ax3x+r12x*ax3y+r13x*ax3z);   /*Z0311=998*/
                ax3yx = (r21x*ax3x+r22x*ax3y+r23x*ax3z);   /*Z0311=999*/
                ax3zx = (r31x*ax3x+r32x*ax3y+r33x*ax3z);   /*Z0311=1000*/

                dqs1x = (dqxx*ax1xx+dqyx*ax1yx+dqzx*ax1zx)/(ax1n*sigx);   /*Z0311=1002*/
                dqs2x = (dqxx*ax2xx+dqyx*ax2yx+dqzx*ax2zx)/(ax2n*sigy);   /*Z0311=1003*/
                dqs3x = (dqxx*ax3xx+dqyx*ax3yx+dqzx*ax3zx)/(ax3n*sigz);   /*Z0311=1004*/

#ifdef PascalComment
                dqs1x = (dqxx*ax1x+dqyx*ax1y+dqzx*ax1z)/(ax1n*sigx);   /*Z0311=1006*/
                dqs2x = (dqxx*ax2x+dqyx*ax2y+dqzx*ax2z)/(ax2n*sigy);   /*Z0311=1007*/
                dqs3x = (dqxx*ax3x+dqyx*ax3y+dqzx*ax3z)/(ax3n*sigz);
#endif

                argx = dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x;   /*Z0311=1010*/
                px = exp(-4*argx/M_PI);   /*Z0311=1011*/
                /* if (argx>2) then px:=eps */   /*Z0311=1012*/
                /*    else px:=exp(-4*argx/pi); */   /*Z0311=1013*/
            }   /*Z0311=1014*/
            if ( i1==8 )
            {  /* disk, general case */   /*Z0311=1015*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1016*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1017*/
                qxl = qx/qn;   /*Z0311=1018*/
                qyl = qy/qn;   /*Z0311=1019*/
                mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);   /*Z0311=1020*/
                mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);   /*Z0311=1021*/
                /* mlx3:=p31*cos(x)*sin(delta)+p33*cos(delta); */   /*Z0311=1022*/
                qnnx = qxl*mlx1+qyl*mlx2;   /*Z0311=1023*/
                argx = sqrt(1.0-qnnx*qnnx+eps9)*l*qn/(z+1);   /*Z0311=1024*/

                if ( sigma<0.15 )
                {   /*Z0311=1026*/
                    if ( argx<0.015 )
                    {   /*Z0311=1027*/
                        px = 1;   /*Z0311=1028*/
                        for ( i=1; i<=30; i++ ) px = px+carr1[i]*pow(argx/2.0,2*i);   /*Z0311=1029*/
                        if ( i3==1 ) px = px*px;   /*Z0311=1030*/
                    }   /*Z0311=1031*/
                    else
                    {   /*Z0311=1032*/
                        if ( i3==0 )
                        {  /* P(q) */   /*Z0311=1033*/
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1034*/
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2*argx))/pow(1+4*argx*argx,(z-2)/2.0);   /*Z0311=1035*/
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2*argx))/pow(1+4*argx*argx,(z-3)/2.0);   /*Z0311=1036*/
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);   /*Z0311=1037*/
                        }   /*Z0311=1038*/
                        if ( i3==1 )
                        {  /* F(q) */   /*Z0311=1039*/
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1+argx*argx,(z-1/2.0)/2.0);   /*Z0311=1040*/
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1+argx*argx,(z-3/2.0)/2.0);   /*Z0311=1041*/
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);   /*Z0311=1042*/
                            px = px3*px3;   /*Z0311=1043*/
                        }   /*Z0311=1044*/
                    }   /*Z0311=1045*/
                }   /*Z0311=1046*/
                else
                {   /*Z0311=1047*/
                    px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1048*/
                    px = 1/(1+M_PI/(4*px1));   /*Z0311=1049*/
                }   /*Z0311=1050*/
            }   /*Z0311=1051*/

            if ( i1==9 )
            {   /* disk, x-axis */   /*Z0311=1053*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1054*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1055*/
                qxl = qx/qn;   /*Z0311=1056*/
                qyl = qy/qn;   /*Z0311=1057*/
                qnnx = qxl*cos(delta)-qyl*sin(x)*sin(delta);   /*Z0311=1058*/
                argx = sqrt(1.0-qnnx*qnnx+eps9)*l*qn/(z+1);   /*Z0311=1059*/

                if ( sigma<0.15 )
                {   /*Z0311=1061*/
                    if ( argx<0.015 )
                    {   /*Z0311=1062*/
                        px = 1;   /*Z0311=1063*/
                        for ( i=1; i<=30; i++ ) px = px+carr1[i]*pow(argx/2.0,2*i);   /*Z0311=1064*/
                        if ( i3==1 ) px = px*px;   /*Z0311=1065*/
                    }   /*Z0311=1066*/
                    else
                    {   /*Z0311=1067*/
                        if ( i3==0 )
                        {  /* P(q) */   /*Z0311=1068*/
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1069*/
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2*argx))/pow(1+4*argx*argx,(z-2)/2.0);   /*Z0311=1070*/
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2*argx))/pow(1+4*argx*argx,(z-3)/2.0);   /*Z0311=1071*/
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);   /*Z0311=1072*/
                        }   /*Z0311=1073*/
                        if ( i3==1 )
                        {  /* F(q) */   /*Z0311=1074*/
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1+argx*argx,(z-1/2.0)/2.0);   /*Z0311=1075*/
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1+argx*argx,(z-3/2.0)/2.0);   /*Z0311=1076*/
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);   /*Z0311=1077*/
                            px = px3*px3;   /*Z0311=1078*/
                        }   /*Z0311=1079*/
                    }   /*Z0311=1080*/
                }   /*Z0311=1081*/
                else
                {   /*Z0311=1082*/
                    px = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1083*/
                    px = 1/(1+M_PI/(4*px1));   /*Z0311=1084*/
                }   /*Z0311=1085*/
            }   /*Z0311=1086*/

            if ( i1==10 )
            {   /* disk, y-axis */   /*Z0311=1088*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1089*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1090*/
                qxl = qx/qn;   /*Z0311=1091*/
                qyl = qy/qn;   /*Z0311=1092*/
                qnnx = qxl*sin(x)*sin(delta)+qyl*cos(delta);   /*Z0311=1093*/
                argx = sqrt(1.0-qnnx*qnnx+eps9)*l*qn/(z+1);   /*Z0311=1094*/
                /*Z0311=1095*/
                if ( sigma<0.15 )
                {   /*Z0311=1096*/
                    if ( argx<0.015 )
                    {   /*Z0311=1097*/
                        px = 1;   /*Z0311=1098*/
                        for ( i=1; i<=30; i++ ) px = px+carr1[i]*pow(argx/2.0,2*i);   /*Z0311=1099*/
                        if ( i3==1 ) px = px*px;   /*Z0311=1100*/
                    }   /*Z0311=1101*/
                    else
                    {   /*Z0311=1102*/
                        if ( i3==0 )
                        {  /* P(q) */   /*Z0311=1103*/
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1104*/
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2*argx))/pow(1+4*argx*argx,(z-2)/2.0);   /*Z0311=1105*/
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2*argx))/pow(1+4*argx*argx,(z-3)/2.0);   /*Z0311=1106*/
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);   /*Z0311=1107*/
                        }   /*Z0311=1108*/
                        if ( i3==1 )
                        {  /* F(q) */   /*Z0311=1109*/
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1+argx*argx,(z-1/2.0)/2.0);   /*Z0311=1110*/
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1+argx*argx,(z-3/2.0)/2.0);   /*Z0311=1111*/
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);   /*Z0311=1112*/
                            px = px3*px3;   /*Z0311=1113*/
                        }   /*Z0311=1114*/
                    }   /*Z0311=1115*/
                }   /*Z0311=1116*/
                else
                {   /*Z0311=1117*/
                    px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1118*/
                    px = 1/(1+M_PI/(4*px1));   /*Z0311=1119*/
                }   /*Z0311=1120*/
            }   /*Z0311=1121*/

            if ( i1==11 )
            {   /* disk, z-axis */   /*Z0311=1123*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1124*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1125*/
                qxl = qx/qn;   /*Z0311=1126*/
                qyl = qy/qn;   /*Z0311=1127*/
                qnnx = qxl*sin(x)*sin(delta)-qyl*cos(x)*sin(delta);   /*Z0311=1128*/
                argx = sqrt(1.0-qnnx*qnnx+eps9)*l*qn/(z+1);   /*Z0311=1129*/
                /*Z0311=1130*/
                if ( sigma<0.15 )
                {   /*Z0311=1131*/
                    if ( argx<0.015 )
                    {   /*Z0311=1132*/
                        px = 1;   /*Z0311=1133*/
                        for ( i=1; i<=30; i++ ) px = px+carr1[i]*pow(argx/2.0,2*i);   /*Z0311=1134*/
                        if ( i3==1 ) px = px*px;   /*Z0311=1135*/
                    }   /*Z0311=1136*/
                    else
                    {   /*Z0311=1137*/
                        if ( i3==0 )
                        {  /* P(q) */   /*Z0311=1138*/
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1139*/
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2*argx))/pow(1+4*argx*argx,(z-2)/2.0);   /*Z0311=1140*/
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2*argx))/pow(1+4*argx*argx,(z-3)/2.0);   /*Z0311=1141*/
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);   /*Z0311=1142*/
                        }   /*Z0311=1143*/
                        if ( i3==1 )
                        {  /* F(q) */   /*Z0311=1144*/
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1+argx*argx,(z-1/2.0)/2.0);   /*Z0311=1145*/
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1+argx*argx,(z-3/2.0)/2.0);   /*Z0311=1146*/
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);   /*Z0311=1147*/
                            px = px3*px3;   /*Z0311=1148*/
                        }   /*Z0311=1149*/
                    }   /*Z0311=1150*/
                }   /*Z0311=1151*/
                else
                {   /*Z0311=1152*/
                    px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);   /*Z0311=1153*/
                    px = 1/(1+M_PI/(4*px1));   /*Z0311=1154*/
                }   /*Z0311=1155*/
            }   /*Z0311=1156*/

            if ( i1==12 )
            {   /* isotropic cube */   /*Z0311=1158*/

#ifdef PascalComment
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1159*/
                qn = sqrt(qx*qx+qy*qy+eps);   /*Z0311=1160*/
                argx1 = qn*sin(delta)*cos(x)*r+eps;   /*Z0311=1161*/
                argx2 = qn*sin(delta)*sin(x)*r+eps;   /*Z0311=1162*/
                argx3 = qn*cos(delta)*r+eps;   /*Z0311=1163*/
                if ( (i3==0) )
                {   /* P(q) */   /*Z0311=1164*/
                    if ( (argx1 < 0.001) ) px1 = 1 else px1 = szave(5,1,argx1,-2,z);   /*Z0311=1165*/
                    if ( (arga2 < 0.001) ) px2 = 1 else px2 = szave(5,1,argx2,-2,z);   /*Z0311=1166*/
                    if ( (argx3 < 0.001) ) px3 = 1 else px3 = szave(5,1,argx3,-2,z);   /*Z0311=1167*/
                    px = px1*px2*px3;   /*Z0311=1168*/
                }   /*Z0311=1169*/
                if ( (i3==1) )
                {   /* F(q) */   /*Z0311=1170*/
                    if ( (argx1 < 0.001) ) px1 = 1 else px1 = szave(3,1,argx1,-1,z);   /*Z0311=1171*/
                    if ( (arga2 < 0.001) ) px2 = 1 else px2 = szave(3,1,argx2,-1,z);   /*Z0311=1172*/
                    if ( (argx3 < 0.001) ) px3 = 1 else px3 = szave(3,1,argx3,-1,z);   /*Z0311=1173*/
                    px = px1*px1*px2*px2*px3*px3;   /*Z0311=1174*/
                }
#endif

                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1177*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1178*/
                argx1 = qn*sin(delta)*cos(x)*r/(z+1);   /*Z0311=1179*/
                argx2 = qn*sin(delta)*sin(x)*r/(z+1);   /*Z0311=1180*/
                argx3 = qn*cos(delta)*l/(z+1);   /*Z0311=1181*/
                if ( i3==0 )
                {   /* P(q) */   /*Z0311=1182*/
                    a1 = 1/(2*z*(z-1));   /*Z0311=1183*/
                    if ( (argx1 < 0.001) )
                        px1 = 1;
                    else   /*Z0311=1184*/
                        px1 = (a1/(argx1*argx1+eps9))*(1-cos((z-1)*atan(2*argx1))/pow(1+4*argx1*argx1,(z-1)/2.0));   /*Z0311=1185*/
                    if ( (argx2 < 0.001) )
                        px2 = 1;
                    else   /*Z0311=1186*/
                        px2 = (a1/(argx2*argx2+eps9))*(1-cos((z-1)*atan(2*argx2))/pow(1+4*argx2*argx2,(z-1)/2.0));   /*Z0311=1187*/
                    if ( (argx3 < 0.001) )
                        px3 = 1;
                    else   /*Z0311=1188*/
                        px3 = (a1/(argx3*argx3+eps9))*(1-cos((z-1)*atan(2*argx3))/pow(1+4*argx3*argx3,(z-1)/2.0));   /*Z0311=1189*/
                    px = px1*px2*px3;   /*Z0311=1190*/
                }   /*Z0311=1191*/
                if ( i3==1 )
                { /* F(q) */   /*Z0311=1192*/
                    a1 = 1/z;   /*Z0311=1193*/
                    if ( (argx1 < 0.001) )
                        px1 = 1;
                    else   /*Z0311=1194*/
                        px1 = (a1/(argx1+eps9))*sin(z*atan(argx1))/pow(1+argx1*argx1,z/2.0);   /*Z0311=1195*/
                    if ( (argx2 < 0.001) )
                        px2 = 1;
                    else   /*Z0311=1196*/
                        px2 = (a1/(argx2+eps9))*sin(z*atan(argx2))/pow(1+argx2*argx2,z/2.0);   /*Z0311=1197*/
                    if ( (arga3 < 0.001) )
                        px3 = 1;
                    else   /*Z0311=1198*/
                        px3 = (a1/(argx3+eps9))*sin(z*atan(argx3))/pow(1+argx3*argx3,z/2.0);   /*Z0311=1199*/
                    px = px1*px1*px2*px2*px3*px3;   /*Z0311=1200*/
                }   /*Z0311=1201*/
            }   /*Z0311=1202*/

            if ( i1==13 )
            {   /* biaxial ellipsoid, isotropic */   /*Z0311=1204*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1205*/
                epsi = l/r;   /*Z0311=1206*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1207*/
                argx = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(x)*cos(x))/(z+1);   /*Z0311=1208*/
                a1 = (1/(2*z*(z-1)*(z-2)*(z-3)));   /*Z0311=1209*/
                px1 = a1*pow(argx,-4)*(1+cos((z-3)*atan(2*argx))/pow(1+4*argx*argx,(z-3)/2.0));   /*Z0311=1210*/
                px2 = (a1/(z-4))*pow(argx,-5)*sin((z-4)*atan(2*argx))/pow(1+4*argx*argx,(z-4)/2.0);   /*Z0311=1211*/
                px3 = (a1/((z-4)*(z-5)))*pow(argx,-6)*(1-cos((z-5)*atan(2*argx))/pow(1+4*argx*argx,(z-5)/2.0));   /*Z0311=1212*/
                px = 9*(px1-2*px2+px3)*sin(x);   /*Z0311=1213*/
            }  /* of biaxial ellipsoid */   /*Z0311=1214*/

            if ( i1==14 )
            {   /* triaxial ellipsoid, isotropic */   /*Z0311=1216*/
                ella = r;   /*Z0311=1217*/
                ellb = l;   /*Z0311=1218*/
                ellc = r/params.p1;   /*Z0311=1219*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1220*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1221*/
                argx = qn*sqrt(pow(ella*cos(x)*sin(delta),2)+pow(ellb*sin(x)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);   /*Z0311=1222*/
                a1 = (1/(2*z*(z-1)*(z-2)*(z-3)));   /*Z0311=1223*/
                px1 = a1*pow(argx,-4)*(1+cos((z-3)*atan(2*argx))/pow(1+4*argx*argx,(z-3)/2.0));   /*Z0311=1224*/
                px2 = (a1/(z-4))*pow(argx,-5)*sin((z-4)*atan(2*argx))/pow(1+4*argx*argx,(z-4)/2.0);   /*Z0311=1225*/
                px3 = (a1/((z-4)*(z-5)))*pow(argx,-6)*(1-cos((z-5)*atan(2*argx))/pow(1+4*argx*argx,(z-5)/2.0));   /*Z0311=1226*/
                px = 9*(px1-2*px2+px3);   /*Z0311=1227*/
            }  /* of triaxial ellipsoid */   /*Z0311=1228*/

            if ( i1==15 )
            {   /* barrel area integration */   /*Z0311=1230*/
                px1 = r*pow(1-pow(x/l,params.p1),1/params.p1);   /*Z0311=1231*/
                argx = -r*pow(x/l,params.p1-1)*pow(1-pow(x/l,params.p1),(1/params.p1)-1)/l;   /*Z0311=1232*/
                px = px1*sqrt(1.0+argx*argx);   /*Z0311=1233*/
            }   /*Z0311=1234*/

            if ( i1==16 )
            {   /* barrel, x-axis */   /*Z0311=1236*/
                z = (1-sigma*sigma)/(sigma*sigma);   /*Z0311=1237*/
                qn = sqrt(qx*qx+qy*qy+eps9);   /*Z0311=1238*/
                argx = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(x)*sin(delta)+eps9));   /*Z0311=1239*/
                if ( i3==0 )
                { /* P(q) */   /*Z0311=1240*/
                    a1 = 9*pow(z+1,4)/(2*z*(z-1)*(z-2)*(z-3));   /*Z0311=1241*/
                    px = (a1/(argx*argx));   /*Z0311=1242*/
                }   /*Z0311=1243*/
                if ( i3==1 )
                { /* F(q) */   /*Z0311=1244*/
                    // 'arga' wird weiter oben gesetzt. Ist dann aber dummerweise der Erstaufruf ...
                    pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1+arga*arga,z/2.0);   /*Z0311=1245*/
                    pa = pa1*pa1;   /*Z0311=1246*/
                    pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1+argb*argb,z/2.0);   /*Z0311=1247*/
                    pb = pb1*pb1;   /*Z0311=1248*/
                }   /*Z0311=1249*/
            }   /*Z0311=1250*/

            sump = sump+px;   /*Z0311=1253*/
            x = x+del;   /*Z0311=1254*/
        }   /*Z0311=1255*/
        pq = 0.5*(pq+(b-a)*sump/tnm);   /*Z0311=1256*/
        trapzdchid_cnt = 2*trapzdchid_cnt;   /*Z0311=1257*/
    }   /*Z0311=1258*/
}   /*Z0311=1259*/


#ifdef USE_OLDVERSION
void CLASSLIB::trapzdchid_OLDVERSION( double a, double b, double r, double sigma, double dbeta, double delta, double /*theta*/, double /*phi*/,
                           double qx, double qy, double qz, double p11, double p12, double p13, double p21, double p22, double p23,
                           double p31, double /*p32*/, double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n,
                           double ax2n, double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
                           double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
                           int /*ordis*/, int /*dim*/, int /*i0*/, int i1, int /*i2*/, int i3, int i4,
                           double *carr1,  // NEU
                           double &pq, int n, int &trapzdchid_cnt ) const
{
    int j;
    double x,sump,/*sumf,sumn,sums,*/del;
    double /*fa,fb,fx,*/pa=0,pb=0,px=0;//,na,nb,nx,sa,sb,sx;
    double z,a1,/*x1z,*/arga,argb,argx,l1,l2,l3;//,qxhkl,qyhkl,qzhkl;
    double mla1,mla2,/*mla3,*/mlb1,mlb2,/*mlb3,*/mlx1,mlx2;//,mlx3;
    //double ma11,ma12,ma13,ma21,ma22,ma23,ma31,ma32,ma33;
    //double mb11,mb12,mb13,mb21,mb22,mb23,mb31,mb32,mb33;
    //double mx11,mx12,mx13,mx21,mx22,mx23,mx31,mx32,mx33;
    double dqxa,dqya,dqza,dqxb,dqyb,dqzb,dqs1a,dqs2a,dqs3a,dqs1b,dqs2b,dqs3b;
    double dqxx,dqyx,dqzx,dqs1x,dqs2x,dqs3x;//,vola,volb,volx;
    double r11a,r12a,r13a,r21a,r22a,r23a,r31a,r32a,r33a;
    double r11b,r12b,r13b,r21b,r22b,r23b,r31b,r32b,r33b;
    double r11x,r12x,r13x,r21x,r22x,r23x,r31x,r32x,r33x;
    double ax1xa,ax1ya,ax1za,ax1xb,ax1yb,ax1zb,ax1xx,ax1yx,ax1zx;
    double ax2xa,ax2ya,ax2za,ax2xb,ax2yb,ax2zb,ax2xx,ax2yx,ax2zx;
    double ax3xa,ax3ya,ax3za,ax3xb,ax3yb,ax3zb,ax3xx,ax3yx,ax3zx;
    double qn,qxhkla,qyhkla,qzhkla,qxhklb,qyhklb,qzhklb,qxhklx,qyhklx,qzhklx;
    //double aexa,aeya,aeza,bexa,beya,beza,cexa,ceya,ceza,aexb,aeyb,aezb,bexb,beyb,bezb,cexb,ceyb,cezb;
    //double aexx,aeyx,aezx,bexx,beyx,bezx,cexx,ceyx,cezx;
    //double pa1,pa2,pa3,pb1,pb2,pb3,px1,px2,px3;
    double qxl,qyl,qnna,qnnb,qnnx;

    //qDebug() << "trapzdchid( ... " << i1 << i3 << i4 << "&pq" << n << trapzdchid_cnt << ")" << dbgCount;

    if ( n==1 )
    {
        switch ( i1 )
        {
        case 1:    //(* cylinder, general case *)
            z=(1-sigma*sigma)/(sigma*sigma);
            a1=1/(2*z*(z-1));
            mla1=p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);
            mla2=p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);
            mlb1=p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);
            mlb2=p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);
            arga=(qx*mla1+qy*mla2+eps7)*r/(z+1);
            argb=(qx*mlb1+qy*mlb2+eps7)*r/(z+1);
            pa=(a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.));
            pb=(a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.));
            break;
        case 2:   //(* cylinder, x-axis *)
            z=(1-sigma*sigma)/(sigma*sigma);
            a1=1/(2*z*(z-1));
            arga=(qx*cos(delta)-qy*sin(a)*sin(delta)+eps7)*r/(z+1);
            argb=(qx*cos(delta)-qy*sin(b)*sin(delta)+eps7)*r/(z+1);
            pa=(a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.));
            pb=(a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.));
            break;
        case 3:   //(* cylinder, y-axis *)
            z=(1-sigma*sigma)/(sigma*sigma);
            a1=1/(2*z*(z-1));
            arga=(qx*sin(a)*sin(delta)+qy*cos(delta)+eps7)*r/(z+1);
            argb=(qx*sin(b)*sin(delta)+qy*cos(delta)+eps7)*r/(z+1);
            pa=(a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.));
            pb=(a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.));
            break;
        case 4:   //(* cylinder, -z-axis *)
            z=(1-sigma*sigma)/(sigma*sigma);
            a1=1/(2*z*(z-1));
            arga=(qx*sin(a)*sin(delta)-qy*cos(a)*sin(delta)+eps7)*r/(z+1);
            argb=(qx*sin(b)*sin(delta)-qy*cos(b)*sin(delta)+eps7)*r/(z+1);
            pa=(a1/(arga*arga))*(1-cos((z-1)*atan(2*arga))/pow(1+4*arga*arga,(z-1)/2.));
            pb=(a1/(argb*argb))*(1-cos((z-1)*atan(2*argb))/pow(1+4*argb*argb,(z-1)/2.));
            break;
        case 5:    //(* general series expansion *)
            mla1=p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);
            mla2=p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);
            mlb1=p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);
            mlb2=p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);
            pa=pow(mla1,i3)*pow(mla2,i4);
            pb=pow(mlb1,i3)*pow(mlb2,i4);
            break;
        case 6:    //(* unit cell rotation *)
            // Bei den Anisotropic Gaussian Tests genutzt
            dqxa=qx-qhkl*(p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));
            dqya=qy-qhkl*(p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));
            dqza=qz-qhkl*(p31*cos(a)*sin(delta)+p33*cos(delta));
            dqxb=qx-qhkl*(p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));
            dqyb=qy-qhkl*(p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));
            dqzb=qz-qhkl*(p31*cos(b)*sin(delta)+p33*cos(delta));
            dqs1a=(dqxa*ax1x+dqya*ax1y+dqza*ax1z)/(ax1n*sigx);;
            dqs2a=(dqxa*ax2x+dqya*ax2y+dqza*ax2z)/(ax2n*sigy);
            dqs3a=(dqxa*ax3x+dqya*ax3y+dqza*ax3z)/(ax3n*sigz);
            dqs1b=(dqxb*ax1x+dqyb*ax1y+dqzb*ax1z)/(ax1n*sigx);;
            dqs2b=(dqxb*ax2x+dqyb*ax2y+dqzb*ax2z)/(ax2n*sigy);
            dqs3b=(dqxb*ax3x+dqyb*ax3y+dqzb*ax3z)/(ax3n*sigz);
            pa=exp(-4*(dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a)/M_PI);
            pb=exp(-4*(dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b)/M_PI);
            break;
        case 7:    //(* fiber rotation *)
            //(* rotation axis, director *)
            l1=r;
            l2=sigma;
            l3=dbeta;

            //(* rotation matrix Ra *)
            r11a=cos(a)+(1-cos(a))*l1*l1;
            r12a=-l3*sin(a)+(1-cos(a))*l1*l2;
            r13a=-l2*sin(a)+(1-cos(a))*l1*l3;
            r21a=l3*sin(a)+(1-cos(a))*l1*l2;
            r22a=cos(a)+(1-cos(a))*l2*l2;
            r23a=l1*sin(a)+(1-cos(a))*l2*l3;
            r31a=l2*sin(a)+(1-cos(a))*l1*l3;
            r32a=-l1*sin(a)+(1-cos(a))*l2*l3;
            r33a=cos(a)+(1-cos(a))*l3*l3;

            //(* rotation matrix Rb *)
            r11b=cos(b)+(1-cos(b))*l1*l1;
            r12b=-l3*sin(b)+(1-cos(b))*l1*l2;
            r13b=-l2*sin(b)+(1-cos(b))*l1*l3;
            r21b=l3*sin(b)+(1-cos(b))*l1*l2;
            r22b=cos(b)+(1-cos(b))*l2*l2;
            r23b=l1*sin(b)+(1-cos(b))*l2*l3;
            r31b=l2*sin(b)+(1-cos(b))*l1*l3;
            r32b=-l1*sin(b)+(1-cos(b))*l2*l3;
            r33b=cos(b)+(1-cos(b))*l3*l3;

            //(* rotate scattering vector *)
            qxhkla=r11a*qxn+r12a*qyn+r13a*qzn;
            qyhkla=r21a*qxn+r22a*qyn+r23a*qzn;
            qzhkla=r31a*qxn+r32a*qyn+r33a*qzn;
            qxhklb=r11b*qxn+r12b*qyn+r13b*qzn;
            qyhklb=r21b*qxn+r22b*qyn+r23b*qzn;
            qzhklb=r31b*qxn+r32b*qyn+r33b*qzn;

            dqxa=qx-qxhkla;
            dqya=qy-qyhkla;
            dqza=qz-qzhkla;

            dqxb=qx-qxhklb;
            dqyb=qy-qyhklb;
            dqzb=qz-qzhklb;

            ax1xa=(r11a*ax1x+r12a*ax1y+r13a*ax1z);
            ax1ya=(r21a*ax1x+r22a*ax1y+r23a*ax1z);
            ax1za=(r31a*ax1x+r32a*ax1y+r33a*ax1z);
            ax1xb=(r11b*ax1x+r12b*ax1y+r13b*ax1z);
            ax1yb=(r21b*ax1x+r22b*ax1y+r23b*ax1z);
            ax1zb=(r31b*ax1x+r32b*ax1y+r33b*ax1z);
            ax2xa=(r11a*ax2x+r12a*ax2y+r13a*ax2z);
            ax2ya=(r21a*ax2x+r22a*ax2y+r23a*ax2z);
            ax2za=(r31a*ax2x+r32a*ax2y+r33a*ax2z);
            ax2xb=(r11b*ax2x+r12b*ax2y+r13b*ax2z);
            ax2yb=(r21b*ax2x+r22b*ax2y+r23b*ax2z);
            ax2zb=(r31b*ax2x+r32b*ax2y+r33b*ax2z);
            ax3xa=(r11a*ax3x+r12a*ax3y+r13a*ax3z);
            ax3ya=(r21a*ax3x+r22a*ax3y+r23a*ax3z);
            ax3za=(r31a*ax3x+r32a*ax3y+r33a*ax3z);
            ax3xb=(r11b*ax3x+r12b*ax3y+r13b*ax3z);
            ax3yb=(r21b*ax3x+r22b*ax3y+r23b*ax3z);
            ax3zb=(r31b*ax3x+r32b*ax3y+r33b*ax3z);

            dqs1a=(dqxa*ax1xa+dqya*ax1ya+dqza*ax1za)/(ax1n*sigx);
            dqs2a=(dqxa*ax2xa+dqya*ax2ya+dqza*ax2za)/(ax2n*sigy);
            dqs3a=(dqxa*ax3xa+dqya*ax3ya+dqza*ax3za)/(ax3n*sigz);
            dqs1b=(dqxb*ax1xb+dqyb*ax1yb+dqzb*ax1zb)/(ax1n*sigx);
            dqs2b=(dqxb*ax2xb+dqyb*ax2yb+dqzb*ax2zb)/(ax2n*sigy);
            dqs3b=(dqxb*ax3xb+dqyb*ax3yb+dqzb*ax3zb)/(ax3n*sigz);

            arga=dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a;
            argb=dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b;
            pa=exp(-4*arga/M_PI);
            pb=exp(-4*argb/M_PI);
            break;
        case 8:    //(* disk, general case *)
            z=(1-sigma*sigma)/(sigma*sigma);
            qn=sqrt(qx*qx+qy*qy+eps7);
            qxl=qx/qn;
            qyl=qy/qn;
            mla1=p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);
            mla2=p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);
            mlb1=p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);
            mlb2=p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);
            qnna=(qxl*mla1+qyl*mla2);
            qnnb=(qxl*mlb1+qyl*mlb2);
            arga=sqrt(1-qnna*qnna+eps7)*r*qn;
            argb=sqrt(1-qnnb*qnnb+eps7)*r*qn;
            pa=1/(1+(M_PI/4.)*pow(arga,3));
            pb=1/(1+(M_PI/4.)*pow(argb,3));
            break;
        case 9:   //(* disk, x-axis *)
            z=(1-sigma*sigma)/(sigma*sigma);
            qn=sqrt(qx*qx+qy*qy+eps7);
            qxl=qx/qn;
            qyl=qy/qn;
            qnna=qxl*cos(delta)-qyl*sin(a)*sin(delta);
            qnnb=qxl*cos(delta)-qyl*sin(b)*sin(delta);
            break;
        case 10:   //(* disk, y-axis *)
            z=(1-sigma*sigma)/(sigma*sigma);
            qn=sqrt(qx*qx+qy*qy+eps7);
            qxl=qx/qn;
            qyl=qy/qn;
            qnna=qxl*sin(a)*sin(delta)+qyl*cos(delta);
            qnnb=qxl*sin(b)*sin(delta)+qyl*cos(delta);
            arga=sqrt(1-qnna*qnna+eps7)*r*qn;
            argb=sqrt(1-qnnb*qnnb+eps7)*r*qn;
            pa=1/(1+(M_PI/4.)*pow(arga,3));
            pb=1/(1+(M_PI/4.)*pow(argb,3));
            break;
        case 11:   //(* disk, -z-axis *)
            z=(1-sigma*sigma)/(sigma*sigma);
            qn=sqrt(qx*qx+qy*qy+eps7);
            qxl=qx/qn;
            qyl=qy/qn;
            qnna=qxl*sin(a)*sin(delta)-qyl*cos(a)*sin(delta);
            qnnb=qxl*sin(b)*sin(delta)-qyl*cos(b)*sin(delta);
            arga=sqrt(1-qnna*qnna+eps7)*r*qn;
            argb=sqrt(1-qnnb*qnnb+eps7)*r*qn;
            pa=1/(1+(M_PI/4.)*pow(arga,3));
            pb=1/(1+(M_PI/4.)*pow(argb,3));
            break;
        case 99:    // TEST x²*y³
            pa = a * a;
            pb = b * b;
            //pq=pa+pb;   // 0 + 4
            //qDebug() << "TEST trapzdchid (1)" << pa << pb << pq;
            //trapzdchid_cnt=1;
            //dbgCount++;
            //return;
            break;
        } // switch ( i1 )
        if ( i1 == 99 )
            pq = pa + pb;
        else
            pq=0.5*(b-a)*(pa+pb);
        trapzdchid_cnt=1;
        DCL( dbgCount++; )
    } // if n == 1
    else
    {
        int tnm = 0; // trapzdchid_cnt;
        del=(b-a)/trapzdchid_cnt;
        x=a+0.5*del;
        //D8L( qDebug() << "trapzdchid" << trapzdit << del << x << dbgCount );
        sump=0.0;
        for ( j=1; j<=trapzdchid_cnt; j++ )
        {
            switch ( i1 )
            {
            case 1:  //(* cylinder, general case *)
                z=(1-sigma*sigma)/(sigma*sigma);
                a1=1/(2*z*(z-1));
                mlx1=p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);
                mlx2=p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);
                argx=(qx*mlx1+qy*mlx2+eps7)*r/(z+1);
                px=(a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.));
                break;
            case 2:   //(* cylinder, x-axis *)
                z=(1-sigma*sigma)/(sigma*sigma);
                a1=1/(2*z*(z-1));
                argx=(qx*cos(delta)-qy*sin(x)*sin(delta)+eps7)*r/(z+1);
                px=(a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.));
                break;
            case 3:   //(* cylinder, y-axis *)
                z=(1-sigma*sigma)/(sigma*sigma);
                a1=1/(2*z*(z-1));
                argx=(qx*sin(x)*sin(delta)+qy*cos(delta)+eps7)*r/(z+1);
                px=(a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.));
                break;
            case 4:   //(* cylinder, -z-axis *)
                z=(1-sigma*sigma)/(sigma*sigma);
                a1=1/(2*z*(z-1));
                argx=(qx*sin(x)*sin(delta)-qy*cos(x)*sin(delta)+eps7)*r/(z+1);
                px=(a1/(argx*argx))*(1-cos((z-1)*atan(2*argx))/pow(1+4*argx*argx,(z-1)/2.));
                break;
            case 5:   //(* general series expansion *)
                mlx1=p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);
                mlx2=p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);
                //mlx3=p31*cos(x)*sin(delta)+p33*cos(delta);
                px=pow(mlx1,i3)*pow(mlx2,i4);
                break;
            case 6:    //(* unit cell rotation *)
                // Bei den Anisotropic Gaussian Tests genutzt
                dqxx=qx-qhkl*(p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta));
                dqyx=qy-qhkl*(p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta));
                dqzx=qz-qhkl*(p31*cos(x)*sin(delta)+p33*cos(delta));
                dqs1x=(dqxx*ax1x+dqyx*ax1y+dqzx*ax1z)/(ax1n*sigx);
                dqs2x=(dqxx*ax2x+dqyx*ax2y+dqzx*ax2z)/(ax2n*sigy);
                dqs3x=(dqxx*ax3x+dqyx*ax3y+dqzx*ax3z)/(ax3n*sigz);
                px=exp(-4*(dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x)/M_PI);

                DSM( if ( fabs(px) < 1.0e-200 )
                {
                    qDebug() << "trapzdchid px" << px << "*** q?" << qx << qy << qz;
                    qDebug() << "dqxx" << dqxx << "=" << qx << "-" << qhkl << "*(" << p11 << "*" << cos(x); // *sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta));
                    qDebug() << "dqyx" << dqyx; //=qy-qhkl*(p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta));
                    qDebug() << "dqzx" << dqzx; //=qz-qhkl*(p31*cos(x)*sin(delta)+p33*cos(delta));
                    qDebug() << "dqs1x" << dqs1x; //=(dqxx*ax1x+dqyx*ax1y+dqzx*ax1z)/(ax1n*sigx);
                    qDebug() << "dqs2x" << dqs2x; //=(dqxx*ax2x+dqyx*ax2y+dqzx*ax2z)/(ax2n*sigy);
                    qDebug() << "dqs3x" << dqs3x; //=(dqxx*ax3x+dqyx*ax3y+dqzx*ax3z)/(ax3n*sigz);
                    qDebug() << "px=exp(" << -4*(dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x)/M_PI << ")";
                    //exit(0);
                }
                //else
                //    qDebug() << "trapzdchid px=exp(" << -4*(dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x)/M_PI << ")";
                )

                break;
            case 7:    //(* fiber unit cell rotation *)
                l1=r;
                l2=sigma;
                l3=dbeta;

                r11x=cos(x)+(1-cos(x))*l1*l1;
                r12x=-l3*sin(x)+(1-cos(x))*l1*l2;
                r13x=-l2*sin(x)+(1-cos(x))*l1*l3;
                r21x=l3*sin(x)+(1-cos(x))*l1*l2;
                r22x=cos(x)+(1-cos(x))*l2*l2;
                r23x=l1*sin(x)+(1-cos(x))*l2*l3;
                r31x=l2*sin(x)+(1-cos(x))*l1*l3;
                r32x=-l1*sin(x)+(1-cos(x))*l2*l3;
                r33x=cos(x)+(1-cos(x))*l3*l3;

                //(* rotate this scattering vector *)
                qxhklx=r11x*qxn+r12x*qyn+r13x*qzn;
                qyhklx=r21x*qxn+r22x*qyn+r23x*qzn;
                qzhklx=r31x*qxn+r32x*qyn+r33x*qzn;

                dqxx=qx-qxhklx;
                dqyx=qy-qyhklx;
                dqzx=qz-qzhklx;

                ax1xx=(r11x*ax1x+r12x*ax1y+r13x*ax1z);
                ax1yx=(r21x*ax1x+r22x*ax1y+r23x*ax1z);
                ax1zx=(r31x*ax1x+r32x*ax1y+r33x*ax1z);
                ax2xx=(r11x*ax2x+r12x*ax2y+r13x*ax2z);
                ax2yx=(r21x*ax2x+r22x*ax2y+r23x*ax2z);
                ax2zx=(r31x*ax2x+r32x*ax2y+r33x*ax2z);
                ax3xx=(r11x*ax3x+r12x*ax3y+r13x*ax3z);
                ax3yx=(r21x*ax3x+r22x*ax3y+r23x*ax3z);
                ax3zx=(r31x*ax3x+r32x*ax3y+r33x*ax3z);

                dqs1x=(dqxx*ax1xx+dqyx*ax1yx+dqzx*ax1zx)/(ax1n*sigx);
                dqs2x=(dqxx*ax2xx+dqyx*ax2yx+dqzx*ax2zx)/(ax2n*sigy);
                dqs3x=(dqxx*ax3xx+dqyx*ax3yx+dqzx*ax3zx)/(ax3n*sigz);

                argx=dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x;
                px=exp(-4*argx/M_PI);
                break;
            case 8:  //(* disk, general case *)
                z=(1-sigma*sigma)/(sigma*sigma);
                qn=sqrt(qx*qx+qy*qy+eps7);
                qxl=qx/qn;
                qyl=qy/qn;
                mlx1=p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);
                mlx2=p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);
                //mlx3:=p31*cos(x)*sin(delta)+p33*cos(delta);
                qnnx=qxl*mlx1+qyl*mlx2;
                argx=sqrt(1-qnnx*qnnx+eps7)*r*qn;
                px=1/(1+(M_PI/4.)*pow(argx,3));
                break;
            case 9:   //(* disk, x-axis *)
                z=(1-sigma*sigma)/(sigma*sigma);
                qn=sqrt(qx*qx+qy*qy+eps7);
                qxl=qx/qn;
                qyl=qy/qn;
                qnnx=qxl*cos(delta)-qyl*sin(x)*sin(delta);
                argx=sqrt(1-qnnx*qnnx+eps7)*r*qn;
                px=1/(1+(M_PI/4.)*pow(argx,3));
                break;
            case 10:   //(* disk, y-axis *)
                z=(1-sigma*sigma)/(sigma*sigma);
                qn=sqrt(qx*qx+qy*qy+eps7);
                qxl=qx/qn;
                qyl=qy/qn;
                qnnx=qxl*sin(x)*sin(delta)+qyl*cos(delta);
                argx=sqrt(1-qnnx*qnnx+eps7)*r*qn;
                px=1/(1+(M_PI/4.)*pow(argx,3));
                break;
            case 11:   //(* disk, -z-axis *)
                z=(1-sigma*sigma)/(sigma*sigma);
                qn=sqrt(qx*qx+qy*qy+eps7);
                qxl=qx/qn;
                qyl=qy/qn;
                qnnx=qxl*sin(x)*sin(delta)-qyl*cos(x)*sin(delta);
                argx=sqrt(1-qnnx*qnnx+eps7)*r*qn;
                px=1/(1+(M_PI/4.)*pow(argx,3));
                break;
            case 99:    // TEST x²*y³
                px = x * x;
                //if ( j == 1 ) qDebug() << "TEST trapzdchid (n)" << px;
                break;
            } // switch i1
            sump += px;
            x += del;
            tnm++;  // Echte Anzahl nehmen, damit die Mittelung unten passt
            DCL(
#ifdef DBGLIMIT
            if ( ++dbgCount > DBGLIMIT ) return;
#else
            dbgCount++;
#endif
                    )
        } // for j
        //double pqold = pq;
        if ( i1 == 99 )
        {
            pq = sump/tnm;
            //if ( trapzdchid_cnt < 10 )
            //qDebug() << "TEST trapzdchid (n)" << trapzdchid_cnt << pq;
        }
        else
            pq=0.5*(pq+(b-a)*sump/tnm);
        //if ( fabs(pq) < 1.0e-170 )
        //    qDebug() << "trapzdchid" << pqold << pq << "ab" << a << b << "sum" << sump << tnm;
        trapzdchid_cnt=2*trapzdchid_cnt;
    } // if n > 1
} /* trapzdchid() */
#endif  // USE_OLDVERSION
#endif // USE_trapzdchid










#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::formpq(double length, double radius, double sigmal, double sigmar, double p1,
                        double rho, double alfa, double theta, double phi, double limql, double limq1,
                        double limq2, double /*limq3*/, double limq4, double limq5, double limq6,
                        double /*limq7*/, double /*limq8*/, double /*limq9*/, double qx, double qy, double qxs,
                        double qys, double q, double norm, double por,
                        int part, int cs, int ordis, int orcase,
                        const double *myarray, // CoeffArrayType
                        double *carr1p, double *carr2p, double *carr3p, // CoeffArrayType
                        double *carr4p, double *carr5p, double *carr6p, // CoeffArrayType
                        double *carr7p, double *carr8p, double *carr9p // CoeffArrayType   /*Z0311=14582*/
                        /*ArrayImax2D carr11pm, ArrayImax2D carr22pm*/) const   /*Z=14910*/
{
    int    ii,jj,n,nser,m,mser,lser,/*indx,*/inmax;
    double pqsum,oldpqsum,binsum,delser,argq,arglq,/*argp1,*/pqr,pql,pq1,pq2,pq3,epsi;
    double ccc1,ccc2,ccc3,vv3,zl,zr,radiusm,argqx,argqy,pqrx,pqry;//,ella,ellb,ellc;
    double cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10;
    double ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10;
    double argbm,nenbm,argbp,nenbp,argep,nenep,argem,nenem,arggp,nengp;
    double arggm,nengm,/*argim,nenim,argip,nenip,*/F121,F122,F123; //,F124,F125,F126;
    double qqn[200+1],/*qqnx[200+1],qqny[200+1],*/fkv[200+1],gam3[200+1]; //: array[0..200] of extended;   /*Z0311=14596*/;
    double F12,F12ser,F12asz,F12sum,F12sez,oldF12sez,F12asy,F12as1z,F12as2z,F12as3z,F12as4z;
    double v,e0,e1,pz2v,pz2v1,pz2v2,lim,lim1,xrz,arg,nen,arg1,nen1,arg2,nen2;
    double a1m,a2m,xijm,arglmz,nenlmz,xijp,arglpz,nenlpz,vvm,rmax,del,delc;
    //double nom,lout,lin,lliph,llipt,phiax,phiin,phiout,philiph,philipt;
    double dim,xrad,xradp,lim2,lim3,lim4,lim5,lim6;
    double a1,b1,b2,b1s,d0,/*d1,*/ee0,ee1,x1z,x12z,x2z,x22z,gb1s;
    double gz1,preg1,preg3,preg4,pzvc,pzvc1,pzvc2,pzac,pzac1,pzc,pzc1,pza,pzva,pzva1; // ,pzac2,dnv0,pvav0,pvav10,pva0;
    double F22sez,oldF22sez,F22,F32sez,oldF32sez,F32;
    double F42sez,oldF42sez,F42,F52sez,oldF52sez,F52,F62sez,oldF62sez,F62;
    double arg11,nen11,arg12,nen12,arg13,nen13;
    double /*arg21,nen21,arg22,nen22,*/arg210,nen210,arg220,nen220,arg23,nen23,arg24,nen24,arg25,nen25,arg26,nen26,arg27,nen27,arg28,nen28;
    double /*F22as1sum1z,F22as1sum2z,*/F22as10z,/*F22as1z,*/F22as1sum1z0,F22as1sum2z0,F22as1z0;
    double a22as21z,F22as21z,a22as22z,F22as22z,a22as23z,F22as23z,a22as24z,F22as24z,F22as20z,F22as2z,/*F22asz,*/F22asz0;
    double /*arg31,nen31,arg32,nen32,*/arg310,nen310,arg320,nen320,arg33,nen33,arg34,nen34,arg35,nen35;
    double /*F32as1sum1z,F32as1sum2z,*/F32as10z,/*F32as1z,*/F32as1sum1z0,F32as1sum2z0,F32as1z0;
    double F32as21z,F32as22z,F32as23z,F32as24z,F32as20z,F32as2z,/*F32asz,*/F32asz0;
    double arg41,nen41,arg42,nen42,/*arg43,nen43,*/arg44,nen44,arg45,nen45;
    double F42as10z,/*F42as1sumz,F42as1z,*/F42as1z0,F42as20z,F42as21,F42as22,/*F42as23,F42as2z,*/F42as2z0;
    double F42as30z,F42as24,/*F42as25,*/F42as26,/*F42as3z,*/F42as3z0,F42as40z,F42as27,F42as28,F42as29,F42as4z,/*F42asz,*/F42asz0;
    double arg51,nen51,arg52,nen52,/*arg53,nen53,*/arg54,nen54,/*arg55,nen55,*/arg56,nen56;
    double arg57,nen57,arg58,nen58,arg59,nen59,arg510,nen510;
    double F52as10z,/*F52as1sumz,F52as1z,*/F52as1z0,F52as20z,F52as21,F52as22,/*F52as23,F52as2z,*/F52as2z0;
    double F52as30z,F52as24,/*F52as25,*/F52as26,/*F52as3z,*/F52as3z0,F52as40z,F52as27,F52as28,F52as29,F52as4z,/*F52asz,*/F52asz0;
    double arg61,nen61,arg62,nen62,/*arg63,nen63,*/arg64,nen64,arg65,nen65;
    double F62as10z,/*F62as1sumz,F62as1z,*/F62as1z0,F62as20z,F62as21,F62as22,/*F62as23,F62as2z,*/F62as2z0;
    double F62as30z,F62as24,/*F62as25,*/F62as26,/*F62as3z,*/F62as3z0,F62as40z,F62as27,F62as28,F62as29,F62as4z,/*F62asz,*/F62asz0;
    double z12v[200+1],a1v[200+1],b1v[200+1],b2v[200+1],b1sv[200+1],/*sum12[200+1],*/sum22[200+1],
           sum32[200+1]; // ,sum42[200+1],sum52[200+1],sum62[200+1]; //: Array[0..200] of extended;   /*Z0311=14623*/;

    double qxn[200+1], qyn[200+1];  // TODO fehlende (bzw. globale) Variablen
    double pqr1, pqr2, pqr3, binsum1;
    const double qz = 1; // TODO
    const double z = 1; // TODO

    //PQRTEST
    //double pqlt1, pqlt2, pqrt1, pqrt2;


    /*begin*/   /*Z=14956*/
    zl = (1-sigmal*sigmal)/(sigmal*sigmal);
    zr = (1-sigmar*sigmar)/(sigmar*sigmar);
    radiusm = radius/p1;   /* outer radius of core/shell particle */

    //ella = radius;   /*Z0311=14630*/
    //ellb = length;   /*Z0311=14631*/
    //ellc = radiusm;   /*Z0311=14632*/

    double zz = (1-sqr(params.sigma))/sqr(params.sigma);  /*Z0311=9337 in coefficients*/
    // noch undefinierte Variablen
    double argpq, c;

    /*Z0311=11518 und 12385*/
    double p11=-cos(phi*M_PI/180.)*cos(theta*M_PI/180.);      // = -cos(theta*pi/180);
    double p12=sin(phi*M_PI/180.);                            // = 0;
    double p13=cos(phi*M_PI/180.)*sin(theta*M_PI/180.);       // =  sin(theta*pi/180);
    double p21=-cos(phi*M_PI/180.);                           // = -1;
    double p22=-sin(phi*M_PI/180.)*cos(theta*M_PI/180.);      // = 0;
    double p23=sin(phi*M_PI/180.)*sin(theta*M_PI/180.);       // = 0;
    double p31=-sin(theta*M_PI/180.);                         // = 0;
    double p32=0;
    double p33=-cos(theta*M_PI/180.);                         // = -cos(theta*pi/180);


    /**************/   /*Z0311=14635*/
    /*** sphere ***/   /*Z=14967*/
    /**************/
    if ( part==0 )
    {
        //qDebug() << "FORMPQ" << cs << q << "limq4,5,6" << limq4 << limq5 << limq6;
        /*** homogeneous sphere ***/   /*Z0311=14639*/
        if ( cs==0 )
        {   /*Z0311=14640*/
            // war 0.6, machte aber bei Test-2 einen schwarzen Kringel...
            // Herr Förster meinte bei einem Test, 0.4 sei besser
            // In der Mail vom 30.05.2022 soll ich mal 0.7 versuchen
            if ( q<0.5*limq4 )  // 221025 - von 0.7 auf 0.5 im Gespräch geändert
            {   /*Z0311=14641*/
                //qDebug() << "FORMPQ" << cs << q << "limq4,5,6" << limq4 << limq5 << limq6;
                double pqsum = 1.0;   /*Z0311=14642*/
                double oldpqsum = 0.0;   /*Z0311=14643*/
                double qqn = 1.0;   /*Z0311=14644*/
                for ( nser=1; nser<=100; nser++ )
                {   /*Z0311=14645*/
                    qqn = qqn*q*q;   /*Z0311=14646*/
                    pqsum = pqsum+carr4p[nser]*qqn;   /*Z0311=14647*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=14648*/
                    if ( delser<0.0001 ) break;   /*Z0311=14649*/
                    oldpqsum = pqsum;   /*Z0311=14650*/
                }   /*Z0311=14651*/
                //50:   /*Z0311=14652*/
                return pqsum;   /*Z0311=14653*/
            }   /*Z0311=14654*/
            else
            {   /*Z0311=14655*/
                argq = q*radius/(zr+1);   /*Z0311=14656*/
                pqr = (1/(2*zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4);   /*Z0311=14657*/
                pq1 = pqr*(1+cos((zr-3)*atan(2*argq))/pow(1+4*argq*argq,(zr-3)/2.0));   /*Z0311=14658*/
                pq2 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2*argq))/pow(1+4*argq*argq,(zr-4)/2.0);   /*Z0311=14659*/
                pq3 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2*argq))/pow(1+4*argq*argq,zr-5)/2.0);   /*Z0311=14660*/
                return 9*(pq1-2*pq2+pq3);   /*Z0311=14661*/
            }   /*Z0311=14662*/
        }   /* of homogeneous sphere*/   /*Z0311=14663*/

        /*** core/shell sphere ***/   /*Z0311=14665*/
        if ( cs==1 )
        {   /*Z0311=14666*/
            /*Z0311=14667*/
            cc1 = sqr(rho);   /*Z0311=14668*/
            cc2 = 2*p1*rho*(1-rho);   /*Z0311=14669*/
            cc3 = sqr(1-rho)*sqr(p1);   /*Z0311=14670*/
            cc4 = -2*sqr(rho);   /*Z0311=14671*/
            cc5 = -2*p1*rho*(1-rho);   /*Z0311=14672*/
            cc6 = sqr(rho);   /*Z0311=14673*/
            cc7 = -2*rho*(1-rho);   /*Z0311=14674*/
            cc8 = -sqr(1-rho)*2*p1;   /*Z0311=14675*/
            cc9 = 2*rho*(1-rho);   /*Z0311=14676*/
            cc10 = sqr(1-rho);   /*Z0311=14677*/

            ccc1 = sqr(1-rho)*pow(p1,6);   /*Z0311=14679*/
            ccc2 = 2*rho*(1-rho)*pow(p1,3);   /*Z0311=14680*/
            ccc3 = rho*rho;   /*Z0311=14681*/
            vv3 = sqr((1-rho)*pow(p1,3)+rho);   /*Z0311=14682*/

            argq = q*radiusm/(zz+1);   /*Z0311=14684*/
            argpq = q*radius/(zz+1);   /*Z0311=14685*/
            pqr = (1/(2*zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4);   /*Z0311=14686*/

            /* F121 sphere */   /*Z0311=14688*/
            if ( q<(0.5*limq4) )
            {   /*Z0311=14689*/
                qqn[0] = 1.0;   /*Z0311=14690*/
                pqsum = 1.0;   /*Z0311=14691*/
                oldpqsum = 0.0;   /*Z0311=14692*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=14693*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=14694*/
                    pqsum = pqsum+qqn[nser]*carr4p[nser];   /*Z0311=14695*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=14696*/
                    if ( delser<0.0001 ) break;   /*Z0311=14697*/
                    oldpqsum = pqsum;   /*Z0311=14698*/
                }   /*Z0311=14699*/
                //51:   /*Z0311=14700*/
                F121 = ccc1*pqsum/vv3;   /*Z0311=14701*/
            }   /*Z0311=14702*/
            else
            {   /*Z0311=14703*/
                ac3 = pqr*(1+cos((zr-3)*atan(2*argpq))/pow(1+4*argpq*argpq,(zr-3)/2.0));   /*Z0311=14704*/
                ac8 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2*argpq))/pow(1+4*argpq*argpq,(zr-4)/2.0);   /*Z0311=14705*/
                ac10 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2*argpq))/pow(1+4*argpq*argpq,(zr-5)/2.0));   /*Z0311=14706*/
                F121 = 9*(cc3*ac3+cc8*ac8+cc10*ac10)/vv3;   /*Z0311=14707*/
            }   /*Z0311=14708*/

            /* F122 sphere */   /*Z0311=14710*/
            if ( q<(limq5/2.0) )
            {   /*Z0311=14711*/
                qqn[0] = 1.0;   /*Z0311=14712*/
                pqsum = 1.0;   /*Z0311=14713*/
                oldpqsum = 0.0;   /*Z0311=14714*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=14715*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=14716*/
                    pqsum = pqsum+qqn[nser]*carr5p[nser];   /*Z0311=14717*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=14718*/
                    if ( delser<0.0001 ) break;   /*Z0311=14719*/
                    oldpqsum = pqsum;   /*Z0311=14720*/
                }   /*Z0311=14721*/
                //52:   /*Z0311=14722*/
                F122 = ccc2*pqsum/vv3;   /*Z0311=14723*/
            }   /*Z0311=14724*/
            else
            {   /*Z0311=14725*/
                argbm = (zr-3)*atan(argpq-argq);   /*Z0311=14726*/
                nenbm = pow(1+sqr(argpq-argq),(zr-3)/2.0);   /*Z0311=14727*/
                argbp = (zr-3)*atan(argpq+argq);   /*Z0311=14728*/
                nenbp = pow(1+sqr(argpq+argq),(zr-3)/2.0);   /*Z0311=14729*/
                ac2 = pqr*(cos(argbm)/nenbm+cos(argbp)/nenbp);   /*Z0311=14730*/
                argep = (zr-4)*atan(argpq+argq);   /*Z0311=14731*/
                nenep = pow(1+sqr(argpq+argq),(zr-4)/2.0);   /*Z0311=14732*/
                argem = (zr-4)*atan(argpq-argq);   /*Z0311=14733*/
                nenem = pow(1+sqr(argpq-argq),(zr-4)/2.0);   /*Z0311=14734*/
                ac5 = (pqr/((zr-4)*argq))*(sin(argep)/nenep-sin(argem)/nenem);   /*Z0311=14735*/
                arggp = (zr-4)*atan(argpq+argq);   /*Z0311=14736*/
                nengp = pow(1+sqr(argpq+argq),(zr-4)/2.0);   /*Z0311=14737*/
                arggm = (zr-4)*atan(argq-argpq);   /*Z0311=14738*/
                nengm = pow(1+sqr(argq-argpq),(zr-4)/2.0);   /*Z0311=14739*/
                ac7 = (pqr/((zr-4)*argq))*(sin(arggp)/nengp-sin(arggm)/nengm);   /*Z0311=14740*/
                //argim = (zr-5)*atan(argpq-argq);   /*Z0311=14741*/
                //nenim = pow(1+sqr(argpq-argq),(zr-5)/2.0);   /*Z0311=14742*/
                //argip = (zr-5)*atan(argpq+argq);   /*Z0311=14743*/
                //nenip = pow(1+sqr(argpq+argq),(zr-5)/2.0);   /*Z0311=14744*/
                ac9 = (pqr/((zr-4)*(zr-5)*argq*argq))*(cos(argbm)/nenbm-cos(argbp)/nenbp);   /*Z0311=14745*/

                F122 = 9*(cc2*ac2+cc5*ac5+cc7*ac7+cc9*ac9)/vv3;   /*Z0311=14747*/
            }   /*Z0311=14748*/

            /* F123 sphere */   /*Z0311=14750*/
            if ( q<(limq6/2.0) )
            {   /*Z0311=14751*/
                qqn[0] = 1.0;   /*Z0311=14752*/
                pqsum = 1.0;   /*Z0311=14753*/
                oldpqsum = 0.0;   /*Z0311=14754*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=14755*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=14756*/
                    pqsum = pqsum+qqn[nser]*carr6p[nser];   /*Z0311=14757*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=14758*/
                    if ( delser<0.0001 ) break;   /*Z0311=14759*/
                    oldpqsum = pqsum;   /*Z0311=14760*/
                }   /*Z0311=14761*/
                //53:   /*Z0311=14762*/
                F123 = ccc3*pqsum/vv3;   /*Z0311=14763*/
            }   /*Z0311=14764*/
            else
            {   /*Z0311=14765*/
                ac1 = pqr*(1+cos((zr-3)*atan(2*argq))/pow(1+4*argq*argq,(zr-3)/2.0));   /*Z0311=14766*/
                ac4 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2*argq))/pow(1+4*argq*argq,(zr-4)/2.0);   /*Z0311=14767*/
                ac6 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2*argq))/pow(1+4*argq*argq,(zr-5)/2.0));   /*Z0311=14768*/
                F123 = 9*(cc1*ac1+cc4*ac4+cc6*ac6)/vv3;   /*Z0311=14769*/
            }   /*Z0311=14770*/
            return F121+F122+F123;   /*Z0311=14771*/
        }  /* of core/shell sphere */   /*Z0311=14772*/

        /*** inhomogeneous core/shell sphere ***/   /*Z0311=14774*/
        if ( cs==2 )
        {   /*Z0311=14775*/

            dim = 3;   /*Z0311=14777*/
            delc = 0.0001;   /*Z0311=14778*/
            xrad = q*radiusm;   /*Z0311=14779*/
            xradp = q*radius;   /*Z0311=14780*/
            x1z = q*radius/(2*(zr+1));   /*Z0311=14781*/
            x12z = x1z*x1z;   /*Z0311=14782*/
            x2z = q*radiusm/(2*(zr+1));   /*Z0311=14783*/
            x22z = x2z*x2z;   /*Z0311=14784*/

            lim = 18*exp(-5*sigmar);   /*Z0311=14786*/
            lim1 = lim;   /*Z0311=14787*/
            lim2 = lim*0.7;   /*Z0311=14788*/
            lim3 = lim;   /*Z0311=14789*/
            lim4 = lim;   /*Z0311=14790*/
            lim5 = lim*0.7;   /*Z0311=14791*/
            lim6 = lim*1.2;   /*Z0311=14792*/

            a1 = (dim-alfa)/2.0;   /*Z0311=14794*/
            b1 = dim/2.0;   /*Z0311=14795*/
            b2 = (dim+2-alfa)/2.0;   /*Z0311=14796*/
            b1s = (dim+2)/2.0;   /*Z0311=14797*/
            v = -b1s+1/2.0;   /*Z0311=14798*/
            c = a1-b1-b2+1/2.0;   /*Z0311=14799*/
            d0 = 1;   /*Z0311=14800*/
            //d1 = a1*(1+a1-b1)*(1+a1-b2);   /*Z0311=14801*/
            e0 = 1.0;   /*Z0311=14802*/
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;   /*Z0311=14803*/
            ee0 = 1.0;   /*Z0311=14804*/
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16*(1-b1s));   /*Z0311=14805*/

            gb1s = 3*sqrt(M_PI)/4.0;   /*Z0311=14807*/
            pz2v = 1/(zr*(zr-1)*(zr-2)*(zr-3));   /*Z0311=14808*/
            pz2v1 = pz2v/(zr-4);   /*Z0311=14809*/
            pz2v2 = pz2v1/(zr-5);   /*Z0311=14810*/

            gz1 = gamma(zr+1);   /*Z0311=14812*/
            preg1 = gb1s/sqrt(M_PI);   /*Z0311=14813*/
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));   /*Z0311=14814*/
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));   /*Z0311=14815*/
            pzvc = gamma(zr+1+v+c)/gz1;   /*Z0311=14816*/
            pzvc1 = gamma(zr+1+v+c-1)/gz1;   /*Z0311=14817*/
            pzvc2 = gamma(zr+1+v+c-2)/gz1;   /*Z0311=14818*/
            pzac = gamma(zr+1-2*a1+c)/gz1;   /*Z0311=14819*/
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;   /*Z0311=14820*/
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;   /*Z0311=14821*/
            pzc = gamma(zr+1+2*c)/gz1;   /*Z0311=14822*/
            pzc1 = gamma(zr+1+2*c-1)/gz1;   /*Z0311=14823*/
            pza = gamma(zr+1-4*a1)/gz1;   /*Z0311=14824*/
            pzva = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=14825*/
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=14826*/
            //dnv0 = 1;   /*Z0311=14827*/
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=14828*/
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=14829*/
            //pva0 = gamma(zr+1-4*a1)/gz1;   /*Z0311=14830*/

            cc1 = 1/(dim*dim);   /*Z0311=14832*/
            cc2 = 2*rho/(dim*(dim-alfa)*pow(p1,dim-alfa));   /*Z0311=14833*/
            cc3 = -2*rho/(dim*(dim-alfa));   /*Z0311=14834*/
            cc4 = rho*rho/((dim-alfa)*(dim-alfa)*pow(p1*p1,dim-alfa));   /*Z0311=14835*/
            cc5 = -2*rho*rho/((dim-alfa)*(dim-alfa)*pow(p1,dim-alfa));   /*Z0311=14836*/
            cc6 = rho*rho/((dim-alfa)*(dim-alfa));   /*Z0311=14837*/
            vv3 = cc1+cc2+cc3+cc4+cc5+cc6;   /*Z0311=14838*/

            /* term #1 series */   /*Z0311=14840*/
            if ( (xradp)<lim1 )
            {   /*Z0311=14841*/
                z12v[0] = 1;   /*Z0311=14842*/
                b1sv[0] = 1;   /*Z0311=14843*/
                fkv[0] = 1;   /*Z0311=14844*/
                gam3[0] = sqrt(M_PI)/2.0;   /*Z0311=14845*/
                qqn[0] = 1.0;   /*Z0311=14846*/
                F12sez = 1.0;   /*Z0311=14847*/
                oldF12sez = 0.0;   /*Z0311=14848*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=14849*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=14850*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=14851*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=14852*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=14853*/
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z0311=14854*/
                    //sum12[n] = 0;   /*Z0311=14855*/
                    /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=14856*/
                    /* sum12[n]:=(9*sqrt(pi)/2)*power(4,n)/((n+3)*(n+2)*(n+3/2)*gam3[n]*fkv[n]); */   /*Z0311=14857*/
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n]; */   /*Z0311=14858*/

                    F12sez = F12sez+carr1p[n]*qqn[n];   /*Z0311=14860*/

                    del = fabs((F12sez-oldF12sez)/F12sez);   /*Z0311=14862*/
                    if ( del<delc ) break;   /*Z0311=14863*/
                    oldF12sez = F12sez;   /*Z0311=14864*/
                }   /*Z0311=14865*/
                //201:   /*Z0311=14866*/
                F12 = F12sez;   /*Z0311=14867*/
            }   /*Z0311=14868*/

            /* term #2 series */   /*Z0311=14870*/
            if ( (xradp)<lim2 )
            {   /*Z0311=14871*/
                z12v[0] = 1;   /*Z0311=14872*/
                a1v[0] = 1;   /*Z0311=14873*/
                b1v[0] = 1;   /*Z0311=14874*/
                b2v[0] = 1;   /*Z0311=14875*/
                b1sv[0] = 1;   /*Z0311=14876*/
                fkv[0] = 1;   /*Z0311=14877*/
                gam3[0] = sqrt(M_PI)/2.0;   /*Z0311=14878*/
                qqn[0] = 1.0;   /*Z0311=14879*/
                F22sez = 1.0;   /*Z0311=14880*/
                oldF22sez = 0.0;   /*Z0311=14881*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=14882*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=14883*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=14884*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=14885*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=14886*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=14887*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=14888*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=14889*/
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z0311=14890*/
                    sum22[n] = 0;   /*Z0311=14891*/
                    /* for m:=0 to n do sum22[n]:=sum22[n]+a1v[n-m]*power(p1*p1,m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=14892*/
                    /* F22sez:=F22sez+power(-x22z,n)*z12v[n]*sum22[n]; */   /*Z0311=14893*/

                    /* for m:=0 to n do sum22[n]:=sum22[n]+power(p1*p1,m)/((n-m+3/2)*(m+(3/2)-(alfa/2))*gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=14895*/
                    /* F22sez:=F22sez+(3*pi*(3-alfa)/16)*power(-x22z,n)*z12v[n]*sum22[n]; */   /*Z0311=14896*/

                    F22sez = F22sez+carr2p[n]*qqn[n];   /*Z0311=14898*/

                    del = fabs((F22sez-oldF22sez)/F22sez);   /*Z0311=14900*/
                    if ( del<delc ) break;   /*Z0311=14901*/
                    oldF22sez = F22sez;   /*Z0311=14902*/
                }   /*Z0311=14903*/
                //202:   /*Z0311=14904*/
                F22 = F22sez;   /*Z0311=14905*/
            }   /*Z0311=14906*/

            /* term #3 series */   /*Z0311=14908*/
            if ( (xradp)<lim3 )
            {   /*Z0311=14909*/
                z12v[0] = 1;   /*Z0311=14910*/
                a1v[0] = 1;   /*Z0311=14911*/
                b1v[0] = 1;   /*Z0311=14912*/
                b2v[0] = 1;   /*Z0311=14913*/
                b1sv[0] = 1;   /*Z0311=14914*/
                fkv[0] = 1;   /*Z0311=14915*/
                gam3[0] = sqrt(M_PI)/2.0;   /*Z0311=14916*/
                qqn[0] = 1.0;   /*Z0311=14917*/
                F32sez = 1.0;   /*Z0311=14918*/
                oldF32sez = 0.0;   /*Z0311=14919*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=14920*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=14921*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=14922*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=14923*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=14924*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=14925*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=14926*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=14927*/
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z0311=14928*/
                    sum32[n] = 0;   /*Z0311=14929*/
                    /* for m:=0 to n do sum32[n]:=sum32[n]+a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=14930*/
                    /* F32sez:=F32sez+power(-x12z,n)*z12v[n]*sum32[n]; */   /*Z0311=14931*/
                    /*Z0311=14932*/
                    /* for m:=0 to n do sum32[n]:=sum32[n]+1/((n-m+3/2)*(m+(3/2)-(alfa/2))*gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=14933*/
                    /* F32sez:=F32sez+(3*pi*(3-alfa)/16)*power(-x12z,n)*z12v[n]*sum32[n]; */   /*Z0311=14934*/
                    /*Z0311=14935*/
                    F32sez = F32sez+carr3p[n]*qqn[n];   /*Z0311=14936*/
                    /*Z0311=14937*/
                    del = fabs((F32sez-oldF32sez)/F32sez);   /*Z0311=14938*/
                    if ( del<delc ) break;   /*Z0311=14939*/
                    oldF32sez = F32sez;   /*Z0311=14940*/
                }   /*Z0311=14941*/
                //203:   /*Z0311=14942*/
                F32 = F32sez;   /*Z0311=14943*/
            }   /*Z0311=14944*/

            /* term #4 series */   /*Z0311=14946*/
            if ( (xradp)<lim4 )
            {   /*Z0311=14947*/
                z12v[0] = 1;   /*Z0311=14948*/
                a1v[0] = 1;   /*Z0311=14949*/
                b1v[0] = 1;   /*Z0311=14950*/
                b2v[0] = 1;   /*Z0311=14951*/
                b1sv[0] = 1;   /*Z0311=14952*/
                fkv[0] = 1;   /*Z0311=14953*/
                gam3[0] = sqrt(M_PI)/2.0;   /*Z0311=14954*/
                qqn[0] = 1.0;   /*Z0311=14955*/
                F42sez = 1.0;   /*Z0311=14956*/
                oldF42sez = 0.0;   /*Z0311=14957*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=14958*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=14959*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=14960*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=14961*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=14962*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=14963*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=14964*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=14965*/
                    gam3[n] = gam3[n-1]*(2*n+1)/2.0;   /*Z0311=14966*/
                    //sum42[n] = 0;   /*Z0311=14967*/
                    /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=14968*/
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n]; */   /*Z0311=14969*/

                    /* for m:=0 to n do sum42[n]:=sum42[n]+1/((n-m+(3/2)-(alfa/2))*(m+(3/2)-(alfa/2))*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]); */   /*Z0311=14971*/
                    /* F42sez:=F42sez+((3-alfa)*(3-alfa)*pi/16)*power(-x22z,n)*z12v[n]*sum42[n]; */   /*Z0311=14972*/

                    F42sez = F42sez+carr4p[n]*qqn[n];   /*Z0311=14974*/

                    del = fabs((F42sez-oldF42sez)/F42sez);   /*Z0311=14976*/
                    if ( del<delc ) break;   /*Z0311=14977*/
                    oldF42sez = F42sez;   /*Z0311=14978*/
                }   /*Z0311=14979*/
                //204:   /*Z0311=14980*/
                F42 = F42sez;   /*Z0311=14981*/
            }   /*Z0311=14982*/

            /* term #5 series */   /*Z0311=14984*/
            if ( (xradp)<lim5 )
            {   /*Z0311=14985*/
                z12v[0] = 1;   /*Z0311=14986*/
                a1v[0] = 1;   /*Z0311=14987*/
                b1v[0] = 1;   /*Z0311=14988*/
                b2v[0] = 1;   /*Z0311=14989*/
                b1sv[0] = 1;   /*Z0311=14990*/
                fkv[0] = 1;   /*Z0311=14991*/
                qqn[0] = 1.0;   /*Z0311=14992*/
                F52sez = 1.0;   /*Z0311=14993*/
                oldF52sez = 0.0;   /*Z0311=14994*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=14995*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=14996*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=14997*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=14998*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=14999*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=15000*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=15001*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=15002*/
                    //sum52[n] = 0;   /*Z0311=15003*/
                    /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=15004*/
                    /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n]; */   /*Z0311=15005*/

                    /* for m:=0 to n do sum52[n]:=sum52[n]+power(p1*p1,n-m)/((n-m+(3/2)-(alfa/2))*(m+(3/2)-(alfa/2))*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]); */   /*Z0311=15007*/
                    /* F52sez:=F52sez+((3-alfa)*(3-alfa)*pi/16)*power(-x22z,n)*z12v[n]*sum52[n]; */   /*Z0311=15008*/

                    F52sez = F52sez+carr5p[n]*qqn[n];   /*Z0311=15010*/

                    del = fabs((F52sez-oldF52sez)/F52sez);   /*Z0311=15012*/
                    if ( del<delc ) break;   /*Z0311=15013*/
                    oldF52sez = F52sez;   /*Z0311=15014*/
                }   /*Z0311=15015*/
                //205:   /*Z0311=15016*/
                F52 = F52sez;   /*Z0311=15017*/
            }   /*Z0311=15018*/

            /* term #6 series */   /*Z0311=15020*/
            if ( (xradp)<lim6 )
            {   /*Z0311=15021*/
                z12v[0] = 1;   /*Z0311=15022*/
                a1v[0] = 1;   /*Z0311=15023*/
                b1v[0] = 1;   /*Z0311=15024*/
                b2v[0] = 1;   /*Z0311=15025*/
                b1sv[0] = 1;   /*Z0311=15026*/
                fkv[0] = 1;   /*Z0311=15027*/
                qqn[0] = 1.0;   /*Z0311=15028*/
                F62sez = 1.0;   /*Z0311=15029*/
                oldF62sez = 0.0;   /*Z0311=15030*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=15031*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=15032*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=15033*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=15034*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=15035*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=15036*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=15037*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=15038*/
                    //sum62[n] = 0;   /*Z0311=15039*/
                    /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=15040*/
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n]; */   /*Z0311=15041*/

                    /* for m:=0 to n do sum62[n]:=sum62[n]+1/((n-m+(3/2)-(alfa/2))*(m+(3/2)-(alfa/2))*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]); */   /*Z0311=15043*/
                    /* F62sez:=F62sez+((3-alfa)*(3-alfa)*pi/16)*power(-x12z,n)*z12v[n]*sum62[n]; */   /*Z0311=15044*/

                    F62sez = F62sez+carr6p[n]*qqn[n];   /*Z0311=15046*/
                    /*Z0311=15047*/
                    del = fabs((F62sez-oldF62sez)/F62sez);   /*Z0311=15048*/
                    if ( del<delc ) break;   /*Z0311=15049*/
                    oldF62sez = F62sez;   /*Z0311=15050*/
                }   /*Z0311=15051*/
                //206:   /*Z0311=15052*/
                F62 = F62sez;   /*Z0311=15053*/
            }   /*Z0311=15054*/

            /*** term #1 asymptote ***/   /*Z0311=15057*/
            if ( xradp>=lim1 )
            {   /*Z0311=15058*/
                arg11 = (zr+2*v+1)*atan(4*x1z);   /*Z0311=15059*/
                nen11 = pow(1+16*x1z*x1z,(zr+2*v+1)/2.0);   /*Z0311=15060*/
                arg12 = (zr+2*v)*atan(4*x1z);   /*Z0311=15061*/
                nen12 = pow(1+16*x1z*x1z,(zr+2*v)/2.0);   /*Z0311=15062*/
                arg13 = (zr+2*v-1)*atan(4*x1z);   /*Z0311=15063*/
                nen13 = pow(1+16*x1z*x1z,(zr+2*v-1)/2.0);   /*Z0311=15064*/

                F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);   /*Z0311=15066*/
                F12as2z = 2*ee0*ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);   /*Z0311=15067*/
                F12as3z = ee1*ee1*(1/(4*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);   /*Z0311=15068*/
                F12asz = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);   /*Z0311=15069*/
                F12 = F12asz;   /*Z0311=15070*/
            }   /*Z0311=15071*/

            /*** term #2 asymptote ***/   /*Z0311=15073*/
            if ( xradp>=lim2 )
            {   /*Z0311=15074*/
                //arg21 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=15075*/
                //nen21 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=15076*/
                //arg22 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=15077*/
                //nen22 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=15078*/
                //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);   /*Z0311=15079*/
                //F22as1sum2z = dnv0*ee1*(1/(2*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);   /*Z0311=15080*/
                F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);   /*Z0311=15081*/
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);   /*Z0311=15082*/

                arg210 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=15084*/
                nen210 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=15085*/
                arg220 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=15086*/
                nen220 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=15087*/
                F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);   /*Z0311=15088*/
                F22as1sum2z0 = ee1*(1/(2*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);   /*Z0311=15089*/
                F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);   /*Z0311=15090*/
                arg23 = (zr+v+c+1)*atan(2*(x1z-x2z));   /*Z0311=15091*/
                nen23 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);   /*Z0311=15092*/
                arg24 = (zr+v+c+1)*atan(2*(x1z+x2z));   /*Z0311=15093*/
                nen24 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);   /*Z0311=15094*/
                arg25 = (zr+v+c)*atan(2*(x1z-x2z));   /*Z0311=15095*/
                nen25 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);   /*Z0311=15096*/
                arg26 = (zr+v+c)*atan(2*(x1z+x2z));   /*Z0311=15097*/
                nen26 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);   /*Z0311=15098*/
                arg27 = (zr+v+c-1)*atan(2*(x1z-x2z));   /*Z0311=15099*/
                nen27 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);   /*Z0311=15100*/
                arg28 = (zr+v+c-1)*atan(2*(x1z+x2z));   /*Z0311=15101*/
                nen28 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);   /*Z0311=15102*/

                a22as21z = (1/2.0)*ee0*e0*pzvc;   /*Z0311=15104*/
                F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);   /*Z0311=15105*/
                a22as22z = (1/2.0)*ee0*e1*(1/(2*x2z))*pzvc1;   /*Z0311=15106*/
                F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);   /*Z0311=15107*/
                a22as23z = (1/2.0)*ee1*e0*(1/(2*x1z))*pzvc1;   /*Z0311=15108*/
                F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);   /*Z0311=15109*/
                a22as24z = (1/2.0)*ee1*e1*(1/(2*x1z))*(1/(2*x2z))*pzvc2;   /*Z0311=15110*/
                F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);   /*Z0311=15111*/
                F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);   /*Z0311=15112*/
                F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);   /*Z0311=15113*/
                //F22asz = F22as1z+F22as2z;   /*Z0311=15114*/
                F22asz0 = F22as1z0+F22as2z;   /*Z0311=15115*/
                F22 = F22asz0;   /*Z0311=15116*/
            }   /*Z0311=15117*/

            /*** term #3 asymptote ***/   /*Z0311=15119*/
            if ( xradp>=lim3 )
            {   /*Z0311=15120*/
                //arg31 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=15121*/
                //nen31 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=15122*/
                //arg32 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=15123*/
                //nen32 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=15124*/
                //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);   /*Z0311=15125*/
                //F32as1sum2z = dnv0*ee1*(1/(2*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);   /*Z0311=15126*/
                F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);   /*Z0311=15127*/
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);   /*Z0311=15128*/

                arg310 = (z+v-2*a1+1)*atan(2*x1z);   /*Z0311=15130*/
                nen310 = pow(1+4*x1z*x1z,(z+v-2*a1+1)/2.0);   /*Z0311=15131*/
                arg320 = (z+v-2*a1)*atan(2*x1z);   /*Z0311=15132*/
                nen320 = pow(1+4*x1z*x1z,(z+v-2*a1)/2.0);   /*Z0311=15133*/
                F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);   /*Z0311=15134*/
                F32as1sum2z0 = ee1*(1/(2*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);   /*Z0311=15135*/
                F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);   /*Z0311=15136*/

                arg33 = (zr+v+c+1)*atan(4*x1z);   /*Z0311=15138*/
                nen33 = pow(1+16*x1z*x1z,(zr+v+c+1)/2.0);   /*Z0311=15139*/
                arg34 = (zr+v+c)*atan(4*x1z);   /*Z0311=15140*/
                nen34 = pow(1+16*x1z*x1z,(zr+v+c)/2.0);   /*Z0311=15141*/
                arg35 = (zr+v+c-1)*atan(4*x1z);   /*Z0311=15142*/
                nen35 = pow(1+16*x1z*x1z,(zr+v+c-1)/2.0);   /*Z0311=15143*/
                F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);   /*Z0311=15144*/
                F32as22z = (1/2.0)*ee0*e1*(1/(2*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);   /*Z0311=15145*/
                F32as23z = (1/2.0)*ee1*e0*(1/(2*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);   /*Z0311=15146*/
                F32as24z = (1/2.0)*ee1*e1*(1/(4*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);   /*Z0311=15147*/
                F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);   /*Z0311=15148*/
                F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);   /*Z0311=15149*/
                //F32asz = F32as1z+F32as2z;   /*Z0311=15150*/
                F32asz0 = F32as1z0+F32as2z;   /*Z0311=15151*/
                F32 = F32asz0;   /*Z0311=15152*/
            }   /*Z0311=15153*/

            /*** term #4 asymptote ***/   /*Z0311=15156*/
            if ( xrad>=lim4 )
            {   /*Z0311=15157*/
                F42as10z = preg4*preg4*pow(x22z,-2*a1);   /*Z0311=15158*/
                //F42as1sumz = pva0;   /*Z0311=15159*/
                //F42as1z = F42as10z*F42as1sumz;   /*Z0311=15160*/
                F42as1z0 = F42as10z*pza;   /*Z0311=15161*/

                arg41 = (zr-2*a1+c+1)*atan(2*x2z);   /*Z0311=15163*/
                nen41 = pow(1+4*x2z*x2z,(zr-2*a1+c+1)/2.0);   /*Z0311=15164*/
                arg42 = (zr-2*a1+c)*atan(2*x2z);   /*Z0311=15165*/
                nen42 = pow(1+4*x2z*x2z,(zr-2*a1+c)/2.0);   /*Z0311=15166*/
                //arg43 = (zr-2*a1+c+3)*atan(2*x2z);   /*Z0311=15167*/
                //nen43 = pow(1+4*x2z*x2z,(zr-2*a1+c+3)/2.0);   /*Z0311=15168*/
                F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);   /*Z0311=15169*/
                F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);   /*Z0311=15170*/
                F42as22 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);   /*Z0311=15171*/
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);   /*Z0311=15172*/
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);   /*Z0311=15173*/
                F42as2z0 = F42as20z*(F42as21+F42as22);   /*Z0311=15174*/

                F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);   /*Z0311=15176*/
                F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);   /*Z0311=15177*/
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);   /*Z0311=15178*/
                F42as26 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);   /*Z0311=15179*/
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);   /*Z0311=15180*/
                F42as3z0 = F42as30z*(F42as24+F42as26);   /*Z0311=15181*/

                F42as40z = preg3*preg3*pow(x2z*x2z,c);   /*Z0311=15183*/
                arg44 = (zr+2*c+1)*atan(4*x2z);   /*Z0311=15184*/
                nen44 = pow(1+16*x2z*x2z,(zr+2*c+1)/2.0);   /*Z0311=15185*/
                arg45 = (zr+2*c)*atan(4*x2z);   /*Z0311=15186*/
                nen45 = pow(1+16*x2z*x2z,(zr+2*c)/2.0);   /*Z0311=15187*/
                F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);   /*Z0311=15188*/
                F42as28 = (1/2.0)*e0*e1*(1/(2*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=15189*/
                F42as29 = (1/2.0)*e1*e0*(1/(2*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=15190*/
                F42as4z = F42as40z*(F42as27+F42as28+F42as29);   /*Z0311=15191*/
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;   /*Z0311=15192*/
                F42asz0 = F42as1z0+F42as2z0+F42as3z0+F42as4z;   /*Z0311=15193*/
                F42 = F42asz0;   /*Z0311=15194*/
            }   /*Z0311=15195*/

            /*** term #5 asymptote ***/   /*Z0311=15198*/
            if ( xradp>=lim5 )
            {   /*Z0311=15199*/
                F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);   /*Z0311=15200*/
                //F52as1sumz = pva0;   /*Z0311=15201*/
                //F52as1z = F52as10z*F52as1sumz;   /*Z0311=15202*/
                F52as1z0 = F52as10z*pza;   /*Z0311=15203*/

                F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);   /*Z0311=15205*/
                arg51 = (zr-2*a1+c+1)*atan(2*x2z);   /*Z0311=15206*/
                nen51 = pow(1+4*x2z*x2z,(zr-2*a1+c+1)/2.0);   /*Z0311=15207*/
                arg52 = (zr-2*a1+c)*atan(2*x2z);   /*Z0311=15208*/
                nen52 = pow(1+4*x2z*x2z,(zr-2*a1+c)/2.0);   /*Z0311=15209*/
                //arg53 = (zr-2*a1+c+3)*atan(2*x2z);   /*Z0311=15210*/
                //nen53 = pow(1+4*x2z*x2z,(zr-2*a1+c+3)/2.0);   /*Z0311=15211*/
                F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);   /*Z0311=15212*/
                F52as22 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);   /*Z0311=15213*/
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);   /*Z0311=15214*/
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);   /*Z0311=15215*/
                F52as2z0 = F52as20z*(F52as21+F52as22);   /*Z0311=15216*/

                F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);   /*Z0311=15218*/
                arg54 = (zr-2*a1+c+1)*atan(2*x1z);   /*Z0311=15219*/
                nen54 = pow(1+4*x1z*x1z,(zr-2*a1+c+1)/2.0);   /*Z0311=15220*/
                //arg55 = (zr-2*a1+c+3)*atan(2*x1z);   /*Z0311=15221*/
                //nen55 = pow(1+4*x1z*x1z,(zr-2*a1+c+3)/2.0);   /*Z0311=15222*/
                arg56 = (zr-2*a1+c)*atan(2*x1z);   /*Z0311=15223*/
                nen56 = pow(1+4*x1z*x1z,(zr-2*a1+c)/2.0);   /*Z0311=15224*/
                F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);   /*Z0311=15225*/
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);   /*Z0311=15226*/
                F52as26 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);   /*Z0311=15227*/
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);   /*Z0311=15228*/
                F52as3z0 = F52as30z*(F52as24+F52as26);   /*Z0311=15229*/

                F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);   /*Z0311=15231*/
                arg57 = (zr+2*c+1)*atan(2*(x1z-x2z));   /*Z0311=15232*/
                nen57 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);   /*Z0311=15233*/
                arg58 = (zr+2*c+1)*atan(2*(x1z+x2z));   /*Z0311=15234*/
                nen58 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);   /*Z0311=15235*/
                arg59 = (zr+2*c)*atan(2*(x1z-x2z));   /*Z0311=15236*/
                nen59 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);   /*Z0311=15237*/
                arg510 = (zr+2*c)*atan(2*(x1z+x2z));   /*Z0311=15238*/
                nen510 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);   /*Z0311=15239*/
                F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);   /*Z0311=15240*/
                F52as28 = (1/2.0)*e0*e1*(1/(2*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);   /*Z0311=15241*/
                F52as29 = (1/2.0)*e1*e0*(1/(2*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);   /*Z0311=15242*/
                F52as4z = F52as40z*(F52as27+F52as28+F52as29);   /*Z0311=15243*/
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;   /*Z0311=15244*/
                F52asz0 = F52as1z0+F52as2z0+F52as3z0+F52as4z;   /*Z0311=15245*/
                F52 = F52asz0;   /*Z0311=15246*/
            }   /*Z0311=15247*/

            /*** term #6 asymptote ***/   /*Z0311=15249*/
            if ( xradp>=lim6 )
            {   /*Z0311=15250*/
                F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);   /*Z0311=15251*/
                //F62as1sumz = pva0;   /*Z0311=15252*/
                //F62as1z = F62as10z*F62as1sumz;   /*Z0311=15253*/
                F62as1z0 = F62as10z*pza;   /*Z0311=15254*/

                F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);   /*Z0311=15256*/
                arg61 = (zr-2*a1+c+1)*atan(2*x1z);   /*Z0311=15257*/
                nen61 = pow(1+4*x1z*x1z,(zr-2*a1+c+1)/2.0);   /*Z0311=15258*/
                arg62 = (zr-2*a1+c)*atan(2*x1z);   /*Z0311=15259*/
                nen62 = pow(1+4*x1z*x1z,(zr-2*a1+c)/2.0);   /*Z0311=15260*/
                //arg63 = (zr-2*a1+c+3)*atan(2*x1z);   /*Z0311=15261*/
                //nen63 = pow(1+4*x1z*x1z,(zr-2*a1+c+3)/2.0);   /*Z0311=15262*/
                F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);   /*Z0311=15263*/
                F62as22 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);   /*Z0311=15264*/
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);   /*Z0311=15265*/
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);   /*Z0311=15266*/
                F62as2z0 = F62as20z*(F62as21+F62as22);   /*Z0311=15267*/

                F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);   /*Z0311=15269*/
                F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);   /*Z0311=15270*/
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);   /*Z0311=15271*/
                F62as26 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);   /*Z0311=15272*/
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);   /*Z0311=15273*/
                F62as3z0 = F62as30z*(F62as24+F62as26);   /*Z0311=15274*/

                F62as40z = preg3*preg3*pow(x1z*x1z,c);   /*Z0311=15276*/
                arg64 = (zr+2*c+1)*atan(4*x1z);   /*Z0311=15277*/
                nen64 = pow(1+16*x1z*x1z,(zr+2*c+1)/2.0);   /*Z0311=15278*/
                arg65 = (zr+2*c)*atan(4*x1z);   /*Z0311=15279*/
                nen65 = pow(1+16*x1z*x1z,(zr+2*c)/2.0);   /*Z0311=15280*/
                F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);   /*Z0311=15281*/
                F62as28 = (1/2.0)*e0*e1*(1/(2*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=15282*/
                F62as29 = (1/2.0)*e1*e0*(1/(2*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=15283*/
                F62as4z = F62as40z*(F62as27+F62as28+F62as29);   /*Z0311=15284*/
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;   /*Z0311=15285*/
                F62asz0 = F62as1z0+F62as2z0+F62as3z0+F62as4z;   /*Z0311=15286*/
                F62 = F62asz0;   /*Z0311=15287*/
            }   /*Z0311=15288*/

            return (cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;   /*Z0311=15290*/

            /* formpq:=pqcoreshellin(1.0,rho,p1,1.0,0.001,alfa,radiusm,3,sigmar,q); */   /*Z0311=15293*/
        } /* of inhomogeneous core/shell sphere */   /*Z0311=15294*/

        /* myelin sphere */   /*Z0311=15298*/
        if ( (cs==3) || (cs==4) )
        {   /*Z0311=15299*/

            /* sphere parameters */   /*Z0311=15301*/
            v = -2;   /*Z0311=15302*/
            e0 = 1;   /*Z0311=15303*/
            e1 = -1;   /*Z0311=15304*/
            preg1 = 3/4.0;   /*Z0311=15305*/
            pz2v = 1/(zr*(zr-1)*(zr-2)*(zr-3));   /*Z0311=15306*/
            pz2v1 = pz2v/(zr-4);   /*Z0311=15307*/
            pz2v2 = pz2v1/(zr-5);   /*Z0311=15308*/
            lim = 18*exp(-5*sigmar);   /*Z0311=15309*/
            lim1 = lim*1.4;   /*Z0311=15310*/
            double rad = myarray[1];   /*Z0311=15311 TODO*/
            inmax = round(myarray[14]);   /*Z0311=15312*/
            vvm = myarray[15];   /*Z0311=15313*/
            rmax = myarray[16];   /*Z0311=15314*/
            double xmax = q*rmax;   /*Z0311=15315 TODO*/

            if ( xmax<(lim1) )
            {   /*Z0311=15317*/
                /* fkv[0]:=1; */   /*Z0311=15318*/
                qqn[0] = 1.0;   /*Z0311=15319*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=15320*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15321*/
                    /* fkv[nser]:=fkv[nser-1]*nser; */   /*Z0311=15322*/
                }   /*Z0311=15323*/

                F12sum = 0.0;   /*Z0311=15325*/
                for ( ii=1; ii<=inmax; ii++ )
                {   /*Z0311=15326*/
                    for ( jj=1; jj<=inmax; jj++ )
                    {   /*Z0311=15327*/
                        F12sez = 1.0;   /*Z0311=15328*/
                        oldF12sez = 1.0;   /*Z0311=15329*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=15330*/
                            pqsum = 0;   /*Z0311=15331*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=15332*/
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]); */   /*Z0311=15333*/
                                pqsum = pqsum+pow(carr7p[ii],2*mser)*pow(carr7p[jj],2*(nser-mser))/(carr6p[mser]*carr6p[nser-mser]);   /*Z0311=15334*/

                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=15336*/
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx]; */   /*Z0311=15337*/
                            }   /*Z0311=15338*/
                            F12sez = F12sez+carr4p[nser]*qqn[nser]*pqsum;   /*Z0311=15339*/
                            delser = fabs((F12sez-oldF12sez)/F12sez);   /*Z0311=15340*/
                            if ( delser<0.0001 ) break;   /*Z0311=15341*/
                            oldF12sez = F12sez;   /*Z0311=15342*/
                        }   /*Z0311=15343*/
                        //250:   /*Z0311=15344*/
                        F12sum = F12sum+carr5p[ii]*carr5p[jj]*F12sez;   /*Z0311=15345*/
                    }   /*Z0311=15346*/
                }   /*Z0311=15347*/
                F12ser = F12sum/vvm;   /*Z0311=15348*/
                F12 = F12ser;   /*Z0311=15349*/
            }   /*Z0311=15350*/
            else
            {   /*Z0311=15351*/
                xrz = q*rad/(zr+1);   /*Z0311=15352*/
                arg = (zr+2*v+1)*atan(2*xrz);   /*Z0311=15353*/
                nen = pow(1+4*xrz*xrz,(zr+2*v+1)/2.0);   /*Z0311=15354*/
                arg1 = (zr+2*v)*atan(2*xrz);   /*Z0311=15355*/
                nen1 = pow(1+4*xrz*xrz,(zr+2*v)/2.0);   /*Z0311=15356*/
                arg2 = (zr+2*v-1)*atan(2*xrz);   /*Z0311=15357*/
                nen2 = pow(1+4*xrz*xrz,(zr+2*v-1)/2.0);   /*Z0311=15358*/

                F12asz = 0.0;   /*Z0311=15360*/
                for ( ii=1; ii<=inmax; ii++ )
                {   /*Z0311=15361*/
                    a1m = carr5p[ii]*pow(carr7p[ii],v);   /*  carr7p[ii]:=pp[ii]; */   /*Z0311=15362*/
                    for ( jj=1; jj<=inmax; jj++ )
                    {   /*Z0311=15363*/
                        a2m = carr5p[jj]*pow(carr7p[jj],v);   /*Z0311=15364*/
                        xijm = (carr3p[ii]-carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii]; */   /*Z0311=15365*/
                        arglmz = (zr+1)*atan(xijm);   /*Z0311=15366*/
                        nenlmz = pow(1+xijm*xijm,(zr+1)/2.0);   /*Z0311=15367*/
                        xijp = (carr3p[ii]+carr3p[jj])*q/(zr+1);   /*Z0311=15368*/
                        arglpz = (zr+1)*atan(xijp);   /*Z0311=15369*/
                        nenlpz = pow(1+xijp*xijp,(zr+1)/2.0);   /*Z0311=15370*/
                        F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(M_PI*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(M_PI*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));   /*Z0311=15371*/
                        F12as2z = e0*e1*(1/(carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));   /*Z0311=15372*/
                        F12as3z = e1*e0*(1/(carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));   /*Z0311=15373*/
                        F12as4z = e1*e1*(1/(carr7p[ii]*carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(M_PI*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(M_PI*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));   /*Z0311=15374*/

                        F12asz = F12asz+a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);   /*Z0311=15376*/
                    }   /*Z0311=15377*/
                }   /*Z0311=15378*/
                F12asy = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;   /*Z0311=15379*/
                F12 = F12asy;   /*Z0311=15380*/
            }   /*Z0311=15381*/
            return F12;   /*Z0311=15382*/

            /* formpq:=polyliposome(llipt,radius,lliph,lin,lout,nom,sigmar,sigmal,phiax,philiph,philipt,phiin,phiout,2,q); */   /*Z0311=15384*/
            /* formpq:=polyliposome(2.0,200,1.0,3.5,3.5,1,sigmar,sigmal,0.001,-0.55,-0.7,0.001,0.001,3,q); */   /*Z0311=15385*/
            /* formpq:=pql; */   /*Z0311=15386*/
        } /* of myelin sphere */   /*Z0311=15387*/

    } /* of sphere */   /*Z0311=15389*/


    /**************/   /*Z0311=15393*/
    /*** ellipsoid ***/   /*Z0311=15394*/
    /**************/   /*Z0311=15395*/
    if ( part==6 )
    {   /*Z0311=15396*/
        epsi = sigmal;   /*Z0311=15397*/
        /*** homogeneous ellipsoid ***/   /*Z0311=15398*/
        if ( ordis==7 )
        {   /* isotropic */   /*Z0311=15399*/
            if ( cs==0 )
            {   /*Z0311=15400*/
                if ( (q<0.4*limq4) )
                {   /*Z0311=15401*/
                    pqsum = 1.0;   /*Z0311=15402*/
                    oldpqsum = 0.0;   /*Z0311=15403*/
                    qqn[0] = 1.0;   /*Z0311=15404*/
                    for ( nser=1; nser<=100; nser++ )
                    {   /*Z0311=15405*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15406*/
                        pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=15407*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15408*/
                        if ( delser<0.0001 ) break;   /*Z0311=15409*/
                        oldpqsum = pqsum;   /*Z0311=15410*/
                    }   /*Z0311=15411*/
                    //90:   /*Z0311=15412*/
                    return pqsum;   /*Z0311=15413*/
                }   /*Z0311=15414*/
                else
                {   /*Z0311=15415*/
                    /* if (q>(4*limq4)) then pq:=(3*pi/(4*zr*(zr-1)*(zr-2)*(zr-3)))*power(q*radius/(zr+1),-4) */   /*Z0311=15416*/
                    /*    else begin */   /*Z0311=15417*/

                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,orcase,0,0,0,carr1p,pq); */   /*Z0311=15419*/

                    const double qz = 1.0; // TODO
                    qrombchid(length,radius,/*p1,*/params.sigma,params.dbeta,epsi,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qx,qy,0,qhkl,ax1.length(),ax2.length(),ax3.length(),ax1.x(),ax1.y(),ax1.z(),ax2.x(),ax2.y(),ax2.z(),ax3.x(),ax3.y(),ax3.z(),sig.x(),sig.y(),sig.z(),ordis,3,7,13,7,0,0,carr1p,pql);   /*Z0311=15421*/
                    return pql;   /*Z0311=15422*/
                }   /*Z0311=15423*/
            } /* of homogeneous sphere*/   /*Z0311=15424*/
        }  /* of isotropic */   /*Z0311=15425*/

        /* perfect */   /*Z0311=15428*/
        if ( ordis==6 )
        {   /*Z0311=15429*/
            if ( orcase==4 )
                pql = 1.0;   /*Z0311=15430*/
            else
            {   /*Z0311=15431*/
                if ( q<(0.6*limq1) )
                {   /*Z0311=15432*/
                    pqsum = 1.0;   /*Z0311=15433*/
                    oldpqsum = 0.0;   /*Z0311=15434*/
                    qxn[0] = 1.0;   /*Z0311=15435*/
                    qyn[0] = 1.0;   /*Z0311=15436*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=15437*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=15438*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=15439*/
                        binsum = 0.0;   /*Z0311=15440*/
                        for ( mser=0; mser<=nser; mser++ )   /*Z0311=15441*/
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[mser]*qyn[nser-mser];   /*Z0311=15442*/
                        pqsum = pqsum+binsum;   /*Z0311=15443*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15444*/
                        if ( delser<0.0001 ) break;   /*Z0311=15445*/
                        oldpqsum = pqsum;   /*Z0311=15446*/
                    }   /*Z0311=15447*/
                    //91:   /*Z0311=15448*/
                    pql = pqsum;   /*Z0311=15449*/
                }   /*Z0311=15450*/
                else
                {   /*Z0311=15451*/
                    arglq = (qxs+qys+eps5)*length/(zl+1);   /*Z0311=15452*/
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2)); */   /*Z0311=15453*/
                    /* pql:=(pi/(2*zl))*(1/arglq); */   /*Z0311=15454*/
                    pql = (1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-1)/2.0));   /*Z0311=15455*/
                }   /*Z0311=15456*/
            }   /*Z0311=15457*/
        }   /* of perfect */   /*Z0311=15458*/

    }   /* of ellipsoid */   /*Z0311=15460*/


    /************/   /*Z0311=15464*/
    /* cylinder */   /*Z=15799*/
    /************/
    if ( part==1 )
    {

#ifndef __CUDACC__
        if ( dbgFlag() )
        {   // Kommt nur einmal pro Thread zu Beginn der Berechnungen
            qDebug() << "formpq:" << "ordis"<<ordis << "orcase"<<orcase << "q"<<q << "limq1"<<limq1 << "limq4"<<limq4
                     << "cs"<<cs << "len/rad"<<length/radius << "norm"<<norm << "limql"<<limql << "q?s"<<qxs<<qys;
        }
#endif

        /*** longitudinal part ***/   /*Z0311=15469*/
        /*** isotropic ***/   /*Z0311=15470*/
        if ( ordis==7 )
        {   /*Z0311=15471*/
            /* exact average */   /*Z0311=15472*/
            if ( (length/radius)<2 )
            {   /*Z0311=15473*/
                if ( q<(5*limq2) )
                {   /*Z0311=15474*/
                    /* Cauchy sum */   /*Z0311=15475*/

#ifdef PascalComment
                    pqsum = 1.0;   /*Z0311=15476*/
                    oldpqsum = 0.0;   /*Z0311=15477*/
                    qqn[0] = 1.0;   /*Z0311=15478*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=15479*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15480*/
                        pqsum = pqsum+carr3p[nser]*qqn[nser];   /*Z0311=15481*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15482*/
                        if ( delser<0.0001 ) goto 601;   /*Z0311=15483*/
                        oldpqsum = pqsum;   /*Z0311=15484*/
                    }
#endif

                    /* double sum */   /*Z0311=15487*/
                    qqn[0] = 1.0;   /*Z0311=15488*/
                    for ( nser=1; nser<=100; nser++ ) qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15489*/
                    pqsum = 0.0;   /*Z0311=15490*/
                    oldpqsum = -10.0;   /*Z0311=15491*/
                    for ( nser=0; nser<=100; nser++ )
                    {   /*Z0311=15492*/
                        binsum = 0.0;   /*Z0311=15493*/
                        for ( mser=0; mser<=100; mser++ ) binsum = binsum+params.CR->carr11pm[nser][mser]*qqn[mser];   /*Z0311=15494*/
                        pqsum = pqsum+carr2p[nser]*qqn[nser]*binsum;   /*Z0311=15495*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15496*/
                        if ( delser<0.0001 ) break;   /*Z0311=15497*/
                        oldpqsum = pqsum;   /*Z0311=15498*/
                    }   /*Z0311=15499*/
                    //601:   /*Z0311=15501*/
                    pql = pqsum;   /*Z0311=15502*/
                }   /*Z0311=15503*/
                else
                {   /*Z0311=15504*/
                    pql = por/pow(q,4);   /*Z0311=15505*/
                }   /*Z0311=15506*/
            }   /*Z0311=15507*/

            /* factorization */   /*Z0311=15509*/
            else
            {   /*Z0311=15510*/
                if ( q<(0.6*limq1) )
                {   /*Z0311=15511*/
                    pqsum = 1.0;   /*Z0311=15512*/
                    oldpqsum = 0.0;   /*Z0311=15513*/
                    qqn[0] = 1.0;   /*Z0311=15514*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=15515*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15516*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser];   /*Z0311=15517*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15518*/
                        if ( delser<0.0001 ) break;   /*Z0311=15519*/
                        oldpqsum = pqsum;   /*Z0311=15520*/
                    }   /*Z0311=15521*/
                    //60:   /*Z0311=15522*/
                    pql = pqsum;   /*Z0311=15523*/
                }   /*Z0311=15524*/
                else
                {   /*Z0311=15525*/
                    arglq = q*length/(zl+1);   /*Z0311=15526*/
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2)); */   /*Z0311=15527*/
                    pql = (M_PI/(2*zl))*(1/arglq);   /*Z0311=15528*/
                    pql = pql-(1/(2*zl*(zl-1)*arglq*arglq))*cos((zl-1)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-1)/2.0);   /*Z0311=15529*/
                }   /*Z0311=15530*/
            }   /*Z0311=15531*/
        }   /* of isotropic */   /*Z0311=15532*/

        /* perfect */   /*Z0311=15534*/
        if ( ordis==/* 6*/ordis_ZDir )
        {   /*Z0311=15535, updates vom 22.06.2022 (220622-cylinder.pas) */
            if ( orcase==4 )
                pql = 1.0;   /*Z0311=15536*/
            else
            {   /*Z0311=15537*/
                if ( limql<(4*limq1) )  /* 220622-cylinder.pas: if (limql<(4*limq1)) then begin */
                {   /*Z0311=15538*/
                    /* if (sqrt(qx*qx*length*length+qy*qy*radius*radius+eps)<10) then begin */   /*Z0311=15539*/
                    pqsum = 1.0;   /*Z0311=15540*/
                    oldpqsum = 0.0;   /*Z0311=15541*/
                    qqn[0] = 1.0;   /*Z0311=15542*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=15543*/
                        qqn[nser] = qqn[nser-1]*(qxs+qys)*(qxs+qys);     /* (qxs,qys)=(qx,0) for x, (0,qy) for y, (qx,qy) for z */   /*Z0311=15544*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser];   /*Z0311=15545*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15546*/
                        if ( delser<0.0001 ) break;   /*Z0311=15547*/
                        oldpqsum = pqsum;   /*Z0311=15548*/
                    }   /*Z0311=15549*/
                    pql = pqsum;   /*Z0311=15551*/
                }   /*Z0311=15552*/
                else
                {   /*Z0311=15553*/
                    arglq = (qxs+qys+eps5)*length/(zl+1);   /*Z0311=15554*/
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2)); */   /*Z0311=15555*/
                    /* pql:=(pi/(2*zl))*(1/arglq); */   /*Z0311=15556*/
                    pql = (1./(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-1)/2.0));   /*Z0311=15557*/
                }   /*Z0311=15558*/
            }   /*Z0311=15559*/
        }   /* of perfect */   /*Z0311=15560*/

        /* general */   /*Z0311=15562*/
        if ( ordis==0 )
        {   /*Z=15898*/
            if ( orcase==4/*z*/ )
                pql = 1.0;
            else
            {
                lim = myarray[17];   /*Z=15900*/

                if ( limql < (2.0 *limq1) ) /* vom 18.11.2022 */
/*PQRTEST*/     //if ( q < limq1 )        /* 20220730-crystal1.pas Z=15901  0.65*limq1 */
                //if ( q < 0.5 )              /* Mail 20220809 */
                {
#ifndef __CUDACC__
                    //qDebug() << "formpq: q<limq1" << q << limql << limq1;
#endif
                    pqsum = 1.0;   /*Z=15902*/
                    oldpqsum = 0.0;
                    qxn[0] = 1.0;
                    qyn[0] = 1.0;
                    if ( orcase==1/*general*/ )
                    {
                        for ( nser=1; nser<=120; nser++ )
                        {
                            qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=15574*/
                            qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=15575*/
                            binsum = 0.0;   /*Z0311=15576*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=15577*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser]; */
                                //binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];   /*Z=15914*/
                                binsum += params.CR->carr11pm[nser][mser]*qxn[mser]*qyn[nser-mser];   /*Mail 20220809*/
                            }   /*Z0311=15581*/
                            pqsum += carr1p[nser]*binsum;   /*Z0311=15582*/

                            if ( isnan(pqsum) || isinf(pqsum) )
                            {
                                DBG( qDebug() << "formpq orcase=1, nser=" << nser << oldpqsum; )
                                pqsum = oldpqsum;
                                break;
                            }

                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15583*/
#ifndef __CUDACC__
                            //qDebug() << "formpq1" << nser << pqsum << oldpqsum << "dif" << delser;
#endif
                            if ( delser<0.0001 ) break;   /*Z0311=15584*/
                            oldpqsum = pqsum;
                        }   /*Z=15920*/
                    }

                    else if ( orcase==2/*x*/ )
                    {  /* x-axis */
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z=15924*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=15591*/
                            qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=15592*/
                            binsum = 0.0;   /*Z0311=15593*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=15594*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=15595*/
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser]; */   /*Z0311=15596*/
                                //binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];   /*Z0311=15597*/
                                binsum += params.CR->carr11pm[nser][mser]*qxn[mser]*qyn[nser-mser];   /*Mail 20220809*/
                            }   /*Z0311=15598*/
                            pqsum += carr1p[nser]*binsum;   /*Z0311=15599*/

                            DBGNAN( if ( isnan(pqsum) ) qDebug() << "NAN1" << nser << carr1p[nser] << binsum; )
                            if ( isnan(pqsum) || isinf(pqsum) )
                            {
                                DBG( qDebug() << "formpq orcase=2, nser=" << nser << oldpqsum; )
                                pqsum = oldpqsum;
                                break;
                            }

                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15600*/
#ifndef __CUDACC__
                            //qDebug() << "formpq2" << nser << pqsum << oldpqsum << "dif" << delser;
#endif
                            if ( delser<0.0001 ) break;   /*Z0311=15601*/
                            oldpqsum = pqsum;   /*Z0311=15602*/
                        }   /*Z=15937*/
                    }

                    else if ( orcase==3/*y*/ )
                    {  /* y-axis */
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z=15942*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs;
                            qyn[nser] = qyn[nser-1]*qys*qys;
                            binsum = 0.0;
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=15612*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=15613*/
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser]; */   /*Z0311=15614*/
                                //binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qyn[mser]*qxn[nser-mser];   /*Z0311=15615*/
                                binsum += params.CR->carr11pm[nser][mser]*qxn[mser]*qyn[nser-mser];   /*Mail 20220809*/
                                // Gleich wie oben, aber carr1pm wurde anders berechnet
                            }   /*Z0311=15616*/
                            pqsum += carr1p[nser]*binsum;   /*Z0311=15617*/

                            if ( isnan(pqsum) || isinf(pqsum) )
                            {
                                DBG( qDebug() << "formpq orcase=3, nser=" << nser << oldpqsum; )
                                pqsum = oldpqsum;
                                break;
                            }

                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15618*/
                            if ( delser<0.0001 ) break;   /*Z0311=15619*/
                            oldpqsum = pqsum;   /*Z0311=15620*/
                        }   /*Z=15955*/
                    }
                    //66:   /*Z0311=15623*/
                    pql = pqsum;   /*Z=15958*/
                    //pqlt1 = pql; //PQRTEST
                }
/*PQRTEST*/     else
                {   /*Z0311=15626*/
                    qrombdeltac(length,radius,/*p1,sigmal, ?alfa,? dbeta,*/theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,orcase,0,0,0,carr1p,pql);   /*Z=15961*/

                    DBGNAN( if ( isnan(pql) ) qDebug() << "NAN2"; )
#ifndef __CUDACC__
                    //qDebug() << "formpq: qrombdeltac" << q << pql << norm << "pql" << pql/norm;
#endif
                    pql = pql/norm;   /*Z0311=15628*/
                    //pqlt2 = pql; //PQRTEST
                }   /*Z0311=15629*/

                //pql = ( q < limq1 ) ? pqlt1 : pqlt2;  //PQRTEST

                //pql = 1;  // ergibt einen homogenen Kreis
                // LOAD PARAMS C:\\SimLab\\sas-crystal\\20220608 - 2D Fits und Fiber-Pattern\\TestParamsCylinder3-foe.ini
            }   /*Z0311=15630*/
        }   /* of general, ordis==0 */   /*Z=15965*/

        /* transverse part */   /*Z=15968*/
        /* homogeneous cylinder */
        if ( cs==0 )
        {
            /* exact average */
            if ( (length/radius)<2 )
                pqr = 1;   /*Z0311=15638*/
            /* factorization */
            else
            {
                if ( q < (0.3 * limq4) )     /*Z=16251 vom 18.11.2022 */
/*PQRTEST*/     //if ( q< limq4 )     /*Z=15975, auch Mail 20220809  vorher:  0.65*limq4 */
                {
                    /* if (sqrt(qx*qx*length*length+qy*qy*radius*radius+eps)<10) then begin */   /*Z0311=15642*/
                    pqsum = 1.0;   /*Z=15977*/
                    oldpqsum = 0.0;
                    double qqn_n = 1.0;
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=15646*/
                        qqn_n = qqn_n*q*q;   /*Z0311=15647*/
                        pqsum += carr4p[nser]*qqn_n;   /*Z0311=15649*/

                        DBGNAN( if ( isnan(pqsum) ) qDebug() << "NAN3" << nser << carr4p[nser] << qqn_n; )
                        if ( isnan(pqsum) || isinf(pqsum) )
                        {
                            pqsum = oldpqsum;
                            break;
                        }

                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15650*/
                        if ( delser<0.0001 ) break;   /*Z0311=15651*/
                        oldpqsum = pqsum;   /*Z0311=15652*/
                    }   /*Z=15988*/
                    //61:   /*Z0311=15654*/
                    pqr = pqsum;
                    //pqrt1 = pqr; //PQRTEST
#ifdef __CUDACC__
                    DBGSPEC( printf( "formpq q=%lf < l4=%lf pql=%lf pqr=%lf\n", q, limq4, pql, pqr ); )
#else
                    DBGSPEC( qDebug() << "formpq q" << q << "< l4" << limq4 << "pql" << pql << "pqr" << pqr << nser << delser; )
#endif
                }
/*PQRTEST*/     else
                {   /*Z0311=15657*/
                    /* q:=sqrt(qxs*qxs+qys*qys+eps);*/       /* for perfect orientation */
                    argpq = q*radius/(zr+1.0);   /*Z=15995*/
                    pqr1 = (1.0/(zr*(zr-1)*(zr-2.)))*pow(argpq,-3.);
                    pqr2 = (1.0/(zr*(zr-1)*(zr-2.)))*pow(argpq,-3.)*sin((zr-2.)*atan(2.*argpq))/pow(1.+4.*argpq*argpq,(zr-2.)/2.0);   /*Z0311=15661*/
                    pqr3 = (1.0/(zr*(zr-1)*(zr-2.)*(zr-3.)))*pow(argpq,-4.)*cos((zr-3.)*atan(2.*argpq))/pow(1.+4.*argpq*argpq,(zr-3.)/2.0);   /*Z0311=15662*/
                    //pqr = 1; // (TEST) (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);   /*Z0311=15663*/
                    pqr = (4./M_PI)*(pqr1-pqr2-(9./8.0)*pqr3);   /*Z=15999 und Mail 20220809*/
                    //pqrt2 = pqr; //PQRTEST

#ifdef __CUDACC__
                    DBGSPEC( printf( "formpq q=%lf >= l4=%lf pql=%lf pqr=%lf\n", q, limq4, pql, pqr ); )
#else
                    DBGSPEC( qDebug() << "formpq q" << q << ">= l4" << limq4 << "pql" << pql << "pq=" << pqr; )
#endif
                    DBGNAN( if ( isnan(pqr) ) qDebug() << "NAN4" << argpq << pqr1 << pqr2 << pqr3; )

                    /* add more terms, if necessary */   /*Z0311=15664*/
                }   /*Z0311=15665*/
            }   /*Z=16002*/
            //pqr = ( q < limq4 ) ? pqrt1 : pqrt2; //PQRTEST
            //qDebug() << "//PQRTEST" << q << "l" << pqlt1 << pqlt2 << "r" << pqrt1 << pqrt2;

            //qDebug() << pql << pqr << pql*pqr;
            //pqr = 1;  streckt die 8 weiter
            //LOAD PARAMS C:\\SimLab\\sas-crystal\\20220608 - 2D Fits und Fiber-Pattern\\TestParamsCylinder3-foe.ini
            return pql * pqr;   /*Z0311=15667*/
            /* formpq:=pql; */   /*Z0311=15668*/
        } /* of homogeneous */   /*Z0311=15669*/

        /* core/shell cylinder */   /*Z0311=15671*/
        if ( cs==1 )
        {   /*Z0311=15672*/
            cc1 = sqr(rho);   /*Z0311=15673*/
            cc2 = 2*p1*rho*(1-rho);   /*Z0311=15674*/
            cc3 = sqr(1-rho)*sqr(p1);   /*Z0311=15675*/
            cc4 = -2*sqr(rho);   /*Z0311=15676*/
            cc5 = -2*p1*rho*(1-rho);   /*Z0311=15677*/
            cc6 = sqr(rho);   /*Z0311=15678*/
            cc7 = -2*rho*(1-rho);   /*Z0311=15679*/
            cc8 = -sqr(1-rho)*2*p1;   /*Z0311=15680*/
            cc9 = 2*rho*(1-rho);   /*Z0311=15681*/
            cc10 = sqr(1-rho);   /*Z0311=15682*/

            ccc1 = sqr(1-rho)*pow(p1,4);   /*Z0311=15684*/
            ccc2 = 2*rho*(1-rho)*pow(p1,2);   /*Z0311=15685*/
            ccc3 = rho*rho;   /*Z0311=15686*/
            vv3 = sqr((1-rho)*pow(p1,2)+rho);   /*Z0311=15687*/

            argq = q*radiusm/(zz+1);   /*Z0311=15689*/
            argpq = q*radius/(zz+1);   /*Z0311=15690*/

            /* F121 cylinder */   /*Z0311=15692*/
            if ( q<(0.7*limq4) )
            {   /*Z0311=15693*/
                /*** series expansion ***/   /*Z0311=15694*/
                pqsum = 1.0;   /*Z0311=15695*/
                oldpqsum = 0.0;   /*Z0311=15696*/
                qqn[0] = 1.0;   /*Z0311=15697*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=15698*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15699*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=15700*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15701*/
                    if ( delser<0.0001 ) break;   /*Z0311=15702*/
                    oldpqsum = pqsum;   /*Z0311=15703*/
                }   /*Z0311=15704*/
                //62:   /*Z0311=15705*/
                F121 = ccc1*pqsum/vv3;   /*Z0311=15706*/
            }   /*Z0311=15707*/
            else
            {   /*Z0311=15708*/
                pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3);   /*Z0311=15709*/
                pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3)*sin((zr-2)*atan(2*argpq))/pow(1+4*argpq*argpq,(zr-2)/2.0);   /*Z0311=15710*/
                pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argpq,-4)*cos((zr-3)*atan(2*argpq))/pow(1+4*argpq*argpq,(zr-3)/2.0);   /*Z0311=15711*/
                pqr = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);   /*Z0311=15712*/
                F121 = ccc1*pqr/vv3;   /*Z0311=15713*/
                /* add more terms, if necessary */   /*Z0311=15714*/
            }   /*Z0311=15715*/

            /* F122 cylinder */   /*Z0311=15717*/
            if ( q<(1.5*limq5) )
            {   /*Z0311=15718*/
                /*** series expansion ***/   /*Z0311=15719*/
                pqsum = 1.0;   /*Z0311=15720*/
                oldpqsum = 0.0;   /*Z0311=15721*/
                qqn[0] = 1.0;   /*Z0311=15722*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=15723*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15724*/
                    pqsum = pqsum+carr5p[nser]*qqn[nser];   /*Z0311=15725*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15726*/
                    if ( delser<0.0001 ) break;   /*Z0311=15727*/
                    oldpqsum = pqsum;   /*Z0311=15728*/
                }   /*Z0311=15729*/
                //63:   /*Z0311=15730*/
                F122 = ccc2*pqsum/vv3;   /*Z0311=15731*/
            }   /*Z0311=15732*/
            else
            {   /*Z0311=15733*/
                argbm = (zr-2)*atan(argpq-argq);   /*Z0311=15734*/
                nenbm = pow(1+sqr(argpq-argq),(zr-2)/2.0);   /*Z0311=15735*/
                argbp = (zr-2)*atan(argpq+argq);   /*Z0311=15736*/
                nenbp = pow(1+sqr(argpq+argq),(zr-2)/2.0);   /*Z0311=15737*/
                argem = (zr-3)*atan(argpq-argq);   /*Z0311=15738*/
                nenem = pow(1+sqr(argpq-argq),(zr-3)/2.0);   /*Z0311=15739*/
                argep = (zr-3)*atan(argpq+argq);   /*Z0311=15740*/
                nenep = pow(1+sqr(argpq+argq),(zr-3)/2.0);   /*Z0311=15741*/
                arggm = (zr-4)*atan(argpq-argq);   /*Z0311=15742*/
                nengm = pow(1+sqr(argpq-argq),(zr-4)/2.0);   /*Z0311=15743*/
                arggp = (zr-4)*atan(argpq+argq);   /*Z0311=15744*/
                nengp = pow(1+sqr(argpq+argq),(zr-4)/2.0);   /*Z0311=15745*/

                pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3)*(cos(argbm)/nenbm-sin(argbp)/nenbp);   /*Z0311=15747*/
                pqr2 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4)*((1-1/p1)*sin(argem)/nenem-(1+1/p1)*cos(argep)/nenep);   /*Z0311=15748*/
                pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)*(zr-4)))*pow(argq,-5)*(1/p1)*(cos(arggm)/nengm-sin(arggp)/nengp);   /*Z0311=15749*/
                pqr = (4/M_PI)*pow(p1,-3/2.0)*(pqr1+(9/16.0)*pqr2+(9/16.0)*(9/16.0)*pqr3);   /*Z0311=15750*/
                F122 = ccc2*pqr/vv3;   /*Z0311=15751*/
            }   /*Z0311=15752*/

            /* F123 cylinder */   /*Z0311=15754*/
            if ( q<(0.6*limq6) )
            {   /*Z0311=15755*/
                /*** series expansion ***/   /*Z0311=15756*/
                pqsum = 1.0;   /*Z0311=15757*/
                oldpqsum = 0.0;   /*Z0311=15758*/
                qqn[0] = 1.0;   /*Z0311=15759*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=15760*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=15761*/
                    pqsum = pqsum+carr6p[nser]*qqn[nser];   /*Z0311=15762*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=15763*/
                    if ( delser<0.0001 ) break;   /*Z0311=15764*/
                    oldpqsum = pqsum;   /*Z0311=15765*/
                }   /*Z0311=15766*/
                //64:   /*Z0311=15767*/
                F123 = ccc3*pqsum/vv3;   /*Z0311=15768*/
            }   /*Z0311=15769*/
            else
            {   /*Z0311=15770*/
                pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3);   /*Z0311=15771*/
                pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3)*sin((zr-2)*atan(2*argq))/pow(1+4*argq*argq,(zr-2)/2.0);   /*Z0311=15772*/
                pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4)*cos((zr-3)*atan(2*argq))/pow(1+4*argq*argq,(zr-3)/2.0);   /*Z0311=15773*/
                pqr = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);   /*Z0311=15774*/
                F123 = ccc3*pqr/vv3;   /*Z0311=15775*/
                /* add more terms, if necessary */   /*Z0311=15776*/
            }   /*Z0311=15777*/
            return pql*(F121+F122+F123);   /*Z0311=15778*/
            /* formpq:=F122; */   /*Z0311=15779*/
        } /* of core/shell */   /*Z0311=15780*/

        /* inhomogeneous core/shell cylinder */   /*Z0311=15782*/
        if ( cs==2 )
        {   /*Z0311=15783*/

            dim = 2;   /*Z0311=15785*/
            delc = 0.0001;   /*Z0311=15786*/
            xrad = q*radiusm;   /*Z0311=15787*/
            xradp = q*radius;   /*Z0311=15788*/
            x1z = q*radius/(2*(zr+1));   /*Z0311=15789*/
            x12z = x1z*x1z;   /*Z0311=15790*/
            x2z = q*radiusm/(2*(zr+1));   /*Z0311=15791*/
            x22z = x2z*x2z;   /*Z0311=15792*/
            /*Z0311=15793*/
            lim = 18*exp(-5*sigmar);   /*Z0311=15794*/
            lim1 = lim;   /*Z0311=15795*/
            lim2 = lim*0.7;   /*Z0311=15796*/
            lim3 = lim;   /*Z0311=15797*/
            lim4 = lim;   /*Z0311=15798*/
            lim5 = lim*0.7;   /*Z0311=15799*/
            lim6 = lim*1.2;   /*Z0311=15800*/

            a1 = (dim-alfa)/2.0;   /*Z0311=15802*/
            b1 = dim/2.0;   /*Z0311=15803*/
            b2 = (dim+2-alfa)/2.0;   /*Z0311=15804*/
            b1s = (dim+2)/2.0;   /*Z0311=15805*/
            v = -b1s+1/2.0;   /*Z0311=15806*/
            c = a1-b1-b2+1/2.0;   /*Z0311=15807*/
            d0 = 1;   /*Z0311=15808*/
            //d1 = a1*(1+a1-b1)*(1+a1-b2);   /*Z0311=15809*/
            e0 = 1.0;   /*Z0311=15810*/
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;   /*Z0311=15811*/
            ee0 = 1.0;   /*Z0311=15812*/
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16*(1-b1s));   /*Z0311=15813*/

            gb1s = 1;   /*Z0311=15815*/
            pz2v = 1/(zr*(zr-1)*(zr-2));   /*Z0311=15816*/
            pz2v1 = pz2v/(zr-3);   /*Z0311=15817*/
            pz2v2 = pz2v1/(zr-4);   /*Z0311=15818*/

            gz1 = gamma(zr+1);   /*Z0311=15820*/
            preg1 = gb1s/sqrt(M_PI);   /*Z0311=15821*/
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));   /*Z0311=15822*/
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));   /*Z0311=15823*/
            pzvc = gamma(zr+1+v+c)/gz1;   /*Z0311=15824*/
            pzvc1 = gamma(zr+1+v+c-1)/gz1;   /*Z0311=15825*/
            pzvc2 = gamma(zr+1+v+c-2)/gz1;   /*Z0311=15826*/
            pzac = gamma(zr+1-2*a1+c)/gz1;   /*Z0311=15827*/
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;   /*Z0311=15828*/
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;   /*Z0311=15829*/
            pzc = gamma(zr+1+2*c)/gz1;   /*Z0311=15830*/
            pzc1 = gamma(zr+1+2*c-1)/gz1;   /*Z0311=15831*/
            pza = gamma(zr+1-4*a1)/gz1;   /*Z0311=15832*/
            pzva = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=15833*/
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=15834*/
            //dnv0 = 1;   /*Z0311=15835*/
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=15836*/
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=15837*/
            //pva0 = gamma(zr+1-4*a1)/gz1;   /*Z0311=15838*/

            cc1 = 1/(dim*dim);   /*Z0311=15840*/
            cc2 = 2*rho/(dim*(dim-alfa)*pow(p1,dim-alfa));   /*Z0311=15841*/
            cc3 = -2*rho/(dim*(dim-alfa));   /*Z0311=15842*/
            cc4 = rho*rho/((dim-alfa)*(dim-alfa)*pow(p1*p1,dim-alfa));   /*Z0311=15843*/
            cc5 = -2*rho*rho/((dim-alfa)*(dim-alfa)*pow(p1,dim-alfa));   /*Z0311=15844*/
            cc6 = rho*rho/((dim-alfa)*(dim-alfa));   /*Z0311=15845*/
            vv3 = cc1+cc2+cc3+cc4+cc5+cc6;   /*Z0311=15846*/

            /* term #1 series */   /*Z=16184 NEU Aug.2022*/
            if ( xradp<lim1 )
            {   /*Z=16185*/
                z12v[0] = 1;   /*Z=16186*/
                b1sv[0] = 1;   /*Z=16187*/
                fkv[0] = 1;   /*Z=16188*/
                F12sez = 1.0;   /*Z=16189*/
                oldF12sez = 0.0;   /*Z=16190*/
                qqn[0] = 1.0;   /*Z=16191*/
                for ( n=1; n<=120; n++ )
                {   /*Z=16192*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z=16193*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z=16194*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=16195*/
                    fkv[n] = fkv[n-1]*n;   /*Z=16196*/
                    //sum12[n] = 0;   /*Z=16197*/
                    /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]); */   /*Z=16198*/
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n]; */   /*Z=16199*/

                    F12sez = F12sez+carr4p[n]*qqn[n];   /*Z=16201*/

                    del = abs((F12sez-oldF12sez)/F12sez);   /*Z=16203*/
                    if ( del<delc ) break; // goto 211;   /*Z=16204*/
                    oldF12sez = F12sez;   /*Z=16205*/
                }   /*Z=16206*/
                //211:   /*Z=16207*/
                F12 = F12sez;   /*Z=16208*/
            }   /*Z=16209*/

            /* term #2 series */   /*Z=16211*/
            if ( xradp<lim2 )
            {   /*Z=16212*/
                z12v[0] = 1;   /*Z=16213*/
                a1v[0] = 1;   /*Z=16214*/
                b1v[0] = 1;   /*Z=16215*/
                b2v[0] = 1;   /*Z=16216*/
                b1sv[0] = 1;   /*Z=16217*/
                fkv[0] = 1;   /*Z=16218*/
                qqn[0] = 1.0;   /*Z=16219*/
                F22sez = 1.0;   /*Z=16220*/
                oldF22sez = 0.0;   /*Z=16221*/
                for ( n=1; n<=120; n++ )
                {   /*Z=16222*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z=16223*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z=16224*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z=16225*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z=16226*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z=16227*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=16228*/
                    fkv[n] = fkv[n-1]*n;   /*Z=16229*/
                    sum22[n] = 0;   /*Z=16230*/
                    for ( m=0; m<=n; m++ ) sum22[n] = sum22[n]+a1v[n-m]*pow(p1*p1,m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);   /*Z=16231*/
                    F22sez = F22sez+pow(-x22z,n)*z12v[n]*sum22[n];   /*Z=16232*/

                    /* F22sez:=F22sez+carr5p[n]*qqn[n]; */   /*Z=16234*/

                    del = abs((F22sez-oldF22sez)/F22sez);   /*Z=16236*/
                    if ( del<delc ) break; //goto 212;   /*Z=16237*/
                    oldF22sez = F22sez;   /*Z=16238*/
                }   /*Z=16239*/
                //212:   /*Z=16240*/
                F22 = F22sez;   /*Z=16241*/
            }   /*Z=16242*/

            /* term #3 series */   /*Z0311=15908*/
            if ( (xradp)<lim3 )
            {   /*Z0311=15909*/
                z12v[0] = 1;   /*Z0311=15910*/
                a1v[0] = 1;   /*Z0311=15911*/
                b1v[0] = 1;   /*Z0311=15912*/
                b2v[0] = 1;   /*Z0311=15913*/
                b1sv[0] = 1;   /*Z0311=15914*/
                fkv[0] = 1;   /*Z0311=15915*/
                qqn[0] = 1.0;   /*Z0311=15916*/
                F32sez = 1.0;   /*Z0311=15917*/
                oldF32sez = 0.0;   /*Z0311=15918*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=15919*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=15920*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=15921*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=15922*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=15923*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=15924*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=15925*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=15926*/
                    sum32[n] = 0;   /*Z0311=15927*/
                    for ( m=0; m<=n; m++ ) sum32[n] = sum32[n]+a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);   /*Z0311=15928*/
                    F32sez = F32sez+pow(-x12z,n)*z12v[n]*sum32[n];   /*Z0311=15929*/

                    /* F32sez:=F32sez+carr6p[n]*qqn[n]; */   /*Z0311=15931*/

                    del = fabs((F32sez-oldF32sez)/F32sez);   /*Z0311=15933*/
                    if ( del<delc ) break;   /*Z0311=15934*/
                    oldF32sez = F32sez;   /*Z0311=15935*/
                }   /*Z0311=15936*/
                //213:   /*Z0311=15937*/
                F32 = F32sez;   /*Z0311=15938*/
            }   /*Z0311=15939*/

            /* term #4 series */   /*Z0311=15941*/
            if ( (xradp)<lim4 )
            {   /*Z0311=15942*/
                z12v[0] = 1;   /*Z0311=15943*/
                a1v[0] = 1;   /*Z0311=15944*/
                b1v[0] = 1;   /*Z0311=15945*/
                b2v[0] = 1;   /*Z0311=15946*/
                b1sv[0] = 1;   /*Z0311=15947*/
                fkv[0] = 1;   /*Z0311=15948*/
                qqn[0] = 1.0;   /*Z0311=15949*/
                F42sez = 1.0;   /*Z0311=15950*/
                oldF42sez = 0.0;   /*Z0311=15951*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=15952*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=15953*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=15954*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=15955*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=15956*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=15957*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=15958*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=15959*/
                    //sum42[n] = 0;   /*Z0311=15960*/
                    /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=15961*/
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n]; */   /*Z0311=15962*/

                    F42sez = F42sez+carr7p[n]*qqn[n];   /*Z0311=15964*/

                    del = fabs((F42sez-oldF42sez)/F42sez);   /*Z0311=15966*/
                    if ( del<delc ) break;   /*Z0311=15967*/
                    oldF42sez = F42sez;   /*Z0311=15968*/
                }   /*Z0311=15969*/
                //214:   /*Z0311=15970*/
                F42 = F42sez;   /*Z0311=15971*/
            }   /*Z0311=15972*/

            /* term #5 series */   /*Z0311=15974*/
            if ( (xradp)<lim5 )
            {   /*Z0311=15975*/
                z12v[0] = 1;   /*Z0311=15976*/
                a1v[0] = 1;   /*Z0311=15977*/
                b1v[0] = 1;   /*Z0311=15978*/
                b2v[0] = 1;   /*Z0311=15979*/
                b1sv[0] = 1;   /*Z0311=15980*/
                fkv[0] = 1;   /*Z0311=15981*/
                qqn[0] = 1.0;   /*Z0311=15982*/
                F52sez = 1.0;   /*Z0311=15983*/
                oldF52sez = 0.0;   /*Z0311=15984*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=15985*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=15986*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=15987*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=15988*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=15989*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=15990*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=15991*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=15992*/
                    //sum52[n] = 0;   /*Z0311=15993*/
                    /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=15994*/
                    /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n]; */   /*Z0311=15995*/

                    F52sez = F52sez+carr8p[n]*qqn[n];   /*Z0311=15997*/

                    del = fabs((F52sez-oldF52sez)/F52sez);   /*Z0311=15999*/
                    if ( del<delc ) break;   /*Z0311=16000*/
                    oldF52sez = F52sez;   /*Z0311=16001*/
                }   /*Z0311=16002*/
                //215:   /*Z0311=16003*/
                F52 = F52sez;   /*Z0311=16004*/
            }   /*Z0311=16005*/

            /* term #6 series */   /*Z0311=16007*/
            if ( (xradp)<lim6 )
            {   /*Z0311=16008*/
                z12v[0] = 1;   /*Z0311=16009*/
                a1v[0] = 1;   /*Z0311=16010*/
                b1v[0] = 1;   /*Z0311=16011*/
                b2v[0] = 1;   /*Z0311=16012*/
                b1sv[0] = 1;   /*Z0311=16013*/
                fkv[0] = 1;   /*Z0311=16014*/
                qqn[0] = 1.0;   /*Z0311=16015*/
                F62sez = 1.0;   /*Z0311=16016*/
                oldF62sez = 0.0;   /*Z0311=16017*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=16018*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=16019*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=16020*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=16021*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=16022*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=16023*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=16024*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=16025*/
                    //sum62[n] = 0;   /*Z0311=16026*/
                    /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=16027*/
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n]; */   /*Z0311=16028*/

                    F62sez = F62sez+carr9p[n]*qqn[n];   /*Z0311=16030*/

                    del = fabs((F62sez-oldF62sez)/F62sez);   /*Z0311=16032*/
                    if ( del<delc ) break;   /*Z0311=16033*/
                    oldF62sez = F62sez;   /*Z0311=16034*/
                }   /*Z0311=16035*/
                //216:   /*Z0311=16036*/
                F62 = F62sez;   /*Z0311=16037*/
            }   /*Z0311=16038*/

            /*** term #1 asymptote ***/   /*Z0311=16041*/
            if ( xradp>=lim1 )
            {   /*Z0311=16042*/
                arg11 = (zr+2*v+1)*atan(4*x1z);   /*Z0311=16043*/
                nen11 = pow(1+16*x1z*x1z,(zr+2*v+1)/2.0);   /*Z0311=16044*/
                arg12 = (zr+2*v)*atan(4*x1z);   /*Z0311=16045*/
                nen12 = pow(1+16*x1z*x1z,(zr+2*v)/2.0);   /*Z0311=16046*/
                arg13 = (zr+2*v-1)*atan(4*x1z);   /*Z0311=16047*/
                nen13 = pow(1+16*x1z*x1z,(zr+2*v-1)/2.0);   /*Z0311=16048*/

                F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);   /*Z0311=16050*/
                F12as2z = 2*ee0*ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);   /*Z0311=16051*/
                F12as3z = ee1*ee1*(1/(4*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);   /*Z0311=16052*/
                F12asz = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);   /*Z0311=16053*/
                F12 = F12asz;   /*Z0311=16054*/
            }   /*Z0311=16055*/

            /*** term #2 asymptote ***/   /*Z0311=16057*/
            if ( xradp>=lim2 )
            {   /*Z0311=16058*/
                //arg21 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=16059*/
                //nen21 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=16060*/
                //arg22 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=16061*/
                //nen22 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=16062*/
                //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);   /*Z0311=16063*/
                //F22as1sum2z = dnv0*ee1*(1/(2*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);   /*Z0311=16064*/
                F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);   /*Z0311=16065*/
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);   /*Z0311=16066*/

                arg210 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=16068*/
                nen210 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=16069*/
                arg220 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=16070*/
                nen220 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=16071*/
                F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);   /*Z0311=16072*/
                F22as1sum2z0 = ee1*(1/(2*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);   /*Z0311=16073*/
                F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);   /*Z0311=16074*/
                arg23 = (zr+v+c+1)*atan(2*(x1z-x2z));   /*Z0311=16075*/
                nen23 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);   /*Z0311=16076*/
                arg24 = (zr+v+c+1)*atan(2*(x1z+x2z));   /*Z0311=16077*/
                nen24 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);   /*Z0311=16078*/
                arg25 = (zr+v+c)*atan(2*(x1z-x2z));   /*Z0311=16079*/
                nen25 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);   /*Z0311=16080*/
                arg26 = (zr+v+c)*atan(2*(x1z+x2z));   /*Z0311=16081*/
                nen26 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);   /*Z0311=16082*/
                arg27 = (zr+v+c-1)*atan(2*(x1z-x2z));   /*Z0311=16083*/
                nen27 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);   /*Z0311=16084*/
                arg28 = (zr+v+c-1)*atan(2*(x1z+x2z));   /*Z0311=16085*/
                nen28 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);   /*Z0311=16086*/

                a22as21z = (1/2.0)*ee0*e0*pzvc;   /*Z0311=16088*/
                F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);   /*Z0311=16089*/
                a22as22z = (1/2.0)*ee0*e1*(1/(2*x2z))*pzvc1;   /*Z0311=16090*/
                F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);   /*Z0311=16091*/
                a22as23z = (1/2.0)*ee1*e0*(1/(2*x1z))*pzvc1;   /*Z0311=16092*/
                F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);   /*Z0311=16093*/
                a22as24z = (1/2.0)*ee1*e1*(1/(2*x1z))*(1/(2*x2z))*pzvc2;   /*Z0311=16094*/
                F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);   /*Z0311=16095*/
                F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);   /*Z0311=16096*/
                F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);   /*Z0311=16097*/
                //F22asz = F22as1z+F22as2z;   /*Z0311=16098*/
                F22asz0 = F22as1z0+F22as2z;   /*Z0311=16099*/
                F22 = F22asz0;   /*Z0311=16100*/
            }   /*Z0311=16101*/

            /*** term #3 asymptote ***/   /*Z0311=16103*/
            if ( xradp>=lim3 )
            {   /*Z0311=16104*/
                //arg31 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=16105*/
                //nen31 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=16106*/
                //arg32 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=16107*/
                //nen32 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=16108*/
                //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);   /*Z0311=16109*/
                //F32as1sum2z = dnv0*ee1*(1/(2*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);   /*Z0311=16110*/
                F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);   /*Z0311=16111*/
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);   /*Z0311=16112*/

                arg310 = (z+v-2*a1+1)*atan(2*x1z);   /*Z0311=16114*/
                nen310 = pow(1+4*x1z*x1z,(z+v-2*a1+1)/2.0);   /*Z0311=16115*/
                arg320 = (z+v-2*a1)*atan(2*x1z);   /*Z0311=16116*/
                nen320 = pow(1+4*x1z*x1z,(z+v-2*a1)/2.0);   /*Z0311=16117*/
                F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);   /*Z0311=16118*/
                F32as1sum2z0 = ee1*(1/(2*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);   /*Z0311=16119*/
                F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);   /*Z0311=16120*/

                arg33 = (zr+v+c+1)*atan(4*x1z);   /*Z0311=16122*/
                nen33 = pow(1+16*x1z*x1z,(zr+v+c+1)/2.0);   /*Z0311=16123*/
                arg34 = (zr+v+c)*atan(4*x1z);   /*Z0311=16124*/
                nen34 = pow(1+16*x1z*x1z,(zr+v+c)/2.0);   /*Z0311=16125*/
                arg35 = (zr+v+c-1)*atan(4*x1z);   /*Z0311=16126*/
                nen35 = pow(1+16*x1z*x1z,(zr+v+c-1)/2.0);   /*Z0311=16127*/
                F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);   /*Z0311=16128*/
                F32as22z = (1/2.0)*ee0*e1*(1/(2*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);   /*Z0311=16129*/
                F32as23z = (1/2.0)*ee1*e0*(1/(2*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);   /*Z0311=16130*/
                F32as24z = (1/2.0)*ee1*e1*(1/(4*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);   /*Z0311=16131*/
                F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);   /*Z0311=16132*/
                F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);   /*Z0311=16133*/
                //F32asz = F32as1z+F32as2z;   /*Z0311=16134*/
                F32asz0 = F32as1z0+F32as2z;   /*Z0311=16135*/
                F32 = F32asz0;   /*Z0311=16136*/
            }   /*Z0311=16137*/

            /*** term #4 asymptote ***/   /*Z0311=16140*/
            if ( xrad>=lim4 )
            {   /*Z0311=16141*/
                F42as10z = preg4*preg4*pow(x22z,-2*a1);   /*Z0311=16142*/
                //F42as1sumz = pva0;   /*Z0311=16143*/
                //F42as1z = F42as10z*F42as1sumz;   /*Z0311=16144*/
                F42as1z0 = F42as10z*pza;   /*Z0311=16145*/

                arg41 = (zr-2*a1+c+1)*atan(2*x2z);   /*Z0311=16147*/
                nen41 = pow(1+4*x2z*x2z,(zr-2*a1+c+1)/2.0);   /*Z0311=16148*/
                arg42 = (zr-2*a1+c)*atan(2*x2z);   /*Z0311=16149*/
                nen42 = pow(1+4*x2z*x2z,(zr-2*a1+c)/2.0);   /*Z0311=16150*/
                //arg43 = (zr-2*a1+c+3)*atan(2*x2z);   /*Z0311=16151*/
                //nen43 = pow(1+4*x2z*x2z,(zr-2*a1+c+3)/2.0);   /*Z0311=16152*/
                F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);   /*Z0311=16153*/
                F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);   /*Z0311=16154*/
                F42as22 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);   /*Z0311=16155*/
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);   /*Z0311=16156*/
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);   /*Z0311=16157*/
                F42as2z0 = F42as20z*(F42as21+F42as22);   /*Z0311=16158*/

                F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);   /*Z0311=16160*/
                F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);   /*Z0311=16161*/
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);   /*Z0311=16162*/
                F42as26 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);   /*Z0311=16163*/
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);   /*Z0311=16164*/
                F42as3z0 = F42as30z*(F42as24+F42as26);   /*Z0311=16165*/

                F42as40z = preg3*preg3*pow(x2z*x2z,c);   /*Z0311=16167*/
                arg44 = (zr+2*c+1)*atan(4*x2z);   /*Z0311=16168*/
                nen44 = pow(1+16*x2z*x2z,(zr+2*c+1)/2.0);   /*Z0311=16169*/
                arg45 = (zr+2*c)*atan(4*x2z);   /*Z0311=16170*/
                nen45 = pow(1+16*x2z*x2z,(zr+2*c)/2.0);   /*Z0311=16171*/
                F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);   /*Z0311=16172*/
                F42as28 = (1/2.0)*e0*e1*(1/(2*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=16173*/
                F42as29 = (1/2.0)*e1*e0*(1/(2*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=16174*/
                F42as4z = F42as40z*(F42as27+F42as28+F42as29);   /*Z0311=16175*/
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;   /*Z0311=16176*/
                F42asz0 = F42as1z0+F42as2z0+F42as3z0+F42as4z;   /*Z0311=16177*/
                F42 = F42asz0;   /*Z0311=16178*/
            }   /*Z0311=16179*/

            /*** term #5 asymptote ***/   /*Z0311=16182*/
            if ( xradp>=lim5 )
            {   /*Z0311=16183*/
                F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);   /*Z0311=16184*/
                //F52as1sumz = pva0;   /*Z0311=16185*/
                //F52as1z = F52as10z*F52as1sumz;   /*Z0311=16186*/
                F52as1z0 = F52as10z*pza;   /*Z0311=16187*/

                F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);   /*Z0311=16189*/
                arg51 = (zr-2*a1+c+1)*atan(2*x2z);   /*Z0311=16190*/
                nen51 = pow(1+4*x2z*x2z,(zr-2*a1+c+1)/2.0);   /*Z0311=16191*/
                arg52 = (zr-2*a1+c)*atan(2*x2z);   /*Z0311=16192*/
                nen52 = pow(1+4*x2z*x2z,(zr-2*a1+c)/2.0);   /*Z0311=16193*/
                //arg53 = (zr-2*a1+c+3)*atan(2*x2z);   /*Z0311=16194*/
                //nen53 = pow(1+4*x2z*x2z,(zr-2*a1+c+3)/2.0);   /*Z0311=16195*/
                F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);   /*Z0311=16196*/
                F52as22 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);   /*Z0311=16197*/
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);   /*Z0311=16198*/
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);   /*Z0311=16199*/
                F52as2z0 = F52as20z*(F52as21+F52as22);   /*Z0311=16200*/

                F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);   /*Z0311=16202*/
                arg54 = (zr-2*a1+c+1)*atan(2*x1z);   /*Z0311=16203*/
                nen54 = pow(1+4*x1z*x1z,(zr-2*a1+c+1)/2.0);   /*Z0311=16204*/
                //arg55 = (zr-2*a1+c+3)*atan(2*x1z);   /*Z0311=16205*/
                //nen55 = pow(1+4*x1z*x1z,(zr-2*a1+c+3)/2.0);   /*Z0311=16206*/
                arg56 = (zr-2*a1+c)*atan(2*x1z);   /*Z0311=16207*/
                nen56 = pow(1+4*x1z*x1z,(zr-2*a1+c)/2.0);   /*Z0311=16208*/
                F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);   /*Z0311=16209*/
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);   /*Z0311=16210*/
                F52as26 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);   /*Z0311=16211*/
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);   /*Z0311=16212*/
                F52as3z0 = F52as30z*(F52as24+F52as26);   /*Z0311=16213*/

                F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);   /*Z0311=16215*/
                arg57 = (zr+2*c+1)*atan(2*(x1z-x2z));   /*Z0311=16216*/
                nen57 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);   /*Z0311=16217*/
                arg58 = (zr+2*c+1)*atan(2*(x1z+x2z));   /*Z0311=16218*/
                nen58 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);   /*Z0311=16219*/
                arg59 = (zr+2*c)*atan(2*(x1z-x2z));   /*Z0311=16220*/
                nen59 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);   /*Z0311=16221*/
                arg510 = (zr+2*c)*atan(2*(x1z+x2z));   /*Z0311=16222*/
                nen510 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);   /*Z0311=16223*/
                F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);   /*Z0311=16224*/
                F52as28 = (1/2.0)*e0*e1*(1/(2*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);   /*Z0311=16225*/
                F52as29 = (1/2.0)*e1*e0*(1/(2*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);   /*Z0311=16226*/
                F52as4z = F52as40z*(F52as27+F52as28+F52as29);   /*Z0311=16227*/
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;   /*Z0311=16228*/
                F52asz0 = F52as1z0+F52as2z0+F52as3z0+F52as4z;   /*Z0311=16229*/
                F52 = F52asz0;   /*Z0311=16230*/
            }   /*Z0311=16231*/

            /*** term #6 asymptote ***/   /*Z0311=16233*/
            if ( xradp>=lim6 )
            {   /*Z0311=16234*/
                F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);   /*Z0311=16235*/
                //F62as1sumz = pva0;   /*Z0311=16236*/
                //F62as1z = F62as10z*F62as1sumz;   /*Z0311=16237*/
                F62as1z0 = F62as10z*pza;   /*Z0311=16238*/

                F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);   /*Z0311=16240*/
                arg61 = (zr-2*a1+c+1)*atan(2*x1z);   /*Z0311=16241*/
                nen61 = pow(1+4*x1z*x1z,(zr-2*a1+c+1)/2.0);   /*Z0311=16242*/
                arg62 = (zr-2*a1+c)*atan(2*x1z);   /*Z0311=16243*/
                nen62 = pow(1+4*x1z*x1z,(zr-2*a1+c)/2.0);   /*Z0311=16244*/
                //arg63 = (zr-2*a1+c+3)*atan(2*x1z);   /*Z0311=16245*/
                //nen63 = pow(1+4*x1z*x1z,(zr-2*a1+c+3)/2.0);   /*Z0311=16246*/
                F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);   /*Z0311=16247*/
                F62as22 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);   /*Z0311=16248*/
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);   /*Z0311=16249*/
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);   /*Z0311=16250*/
                F62as2z0 = F62as20z*(F62as21+F62as22);   /*Z0311=16251*/

                F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);   /*Z0311=16253*/
                F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);   /*Z0311=16254*/
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);   /*Z0311=16255*/
                F62as26 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);   /*Z0311=16256*/
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);   /*Z0311=16257*/
                F62as3z0 = F62as30z*(F62as24+F62as26);   /*Z0311=16258*/

                F62as40z = preg3*preg3*pow(x1z*x1z,c);   /*Z0311=16260*/
                arg64 = (zr+2*c+1)*atan(4*x1z);   /*Z0311=16261*/
                nen64 = pow(1+16*x1z*x1z,(zr+2*c+1)/2.0);   /*Z0311=16262*/
                arg65 = (zr+2*c)*atan(4*x1z);   /*Z0311=16263*/
                nen65 = pow(1+16*x1z*x1z,(zr+2*c)/2.0);   /*Z0311=16264*/
                F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);   /*Z0311=16265*/
                F62as28 = (1/2.0)*e0*e1*(1/(2*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=16266*/
                F62as29 = (1/2.0)*e1*e0*(1/(2*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=16267*/
                F62as4z = F62as40z*(F62as27+F62as28+F62as29);   /*Z0311=16268*/
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;   /*Z0311=16269*/
                F62asz0 = F62as1z0+F62as2z0+F62as3z0+F62as4z;   /*Z0311=16270*/
                F62 = F62asz0;   /*Z0311=16271*/
            }   /*Z0311=16272*/

            return pql*(cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;   /*Z0311=16274*/
            /* formpq:=pql*(cc2*F22)/vv3; */   /*Z0311=16275*/
            /* formpq:=pqcoreshellin(1.0,rho,p1,1.0,0.001,alfa,radiusm,2,sigmar,q); */   /*Z0311=16278*/
        } /* of inhomogeneous core/shell */   /*Z0311=16280*/

        /* myelin cylinder */   /*Z0311=16283*/
        if ( (cs==3) || (cs==4) )
        {   /*Z0311=16284*/

            /* cylinder parameters */   /*Z0311=16286*/
            v = -3/2.0;   /*Z0311=16287*/
            e0 = 1;   /*Z0311=16288*/
            e1 = -9/16.0;   /*Z0311=16289*/
            preg1 = 1/sqrt(M_PI);   /*Z0311=16290*/
            pz2v = 1/(zr*(zr-1)*(zr-2));   /*Z0311=16291*/
            pz2v1 = pz2v/(zr-3);   /*Z0311=16292*/
            pz2v2 = pz2v1/(zr-4);   /*Z0311=16293*/
            lim = 18*exp(-5*sigmar);   /*Z0311=16294*/
            lim1 = lim*1.2;   /*Z0311=16295*/
            double rad = myarray[1];   /*Z0311=16296 TODO*/
            inmax = round(myarray[14]);   /*Z0311=16297*/
            vvm = myarray[15];   /*Z0311=16298*/
            rmax = myarray[16];   /*Z0311=16299*/
            double xmax = q*rmax;   /*Z0311=16300 TODO*/

            if ( xmax<(lim1) )
            {   /*Z0311=16302*/
                /* fkv[0]:=1; */   /*Z0311=16303*/
                qqn[0] = 1.0;   /*Z0311=16304*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=16305*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=16306*/
                    /* fkv[nser]:=fkv[nser-1]*nser; */   /*Z0311=16307*/
                }   /*Z0311=16308*/

                F12sum = 0.0;   /*Z0311=16310*/
                for ( ii=1; ii<=inmax; ii++ )
                {   /*Z0311=16311*/
                    for ( jj=1; jj<=inmax; jj++ )
                    {   /*Z0311=16312*/
                        F12sez = 1.0;   /*Z0311=16313*/
                        oldF12sez = 1.0;   /*Z0311=16314*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=16315*/
                            pqsum = 0;   /*Z0311=16316*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=16317*/
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]); */   /*Z0311=16318*/
                                pqsum = pqsum+pow(carr7p[ii],2*mser)*pow(carr7p[jj],2*(nser-mser))/(carr6p[mser]*carr6p[nser-mser]);   /*Z0311=16319*/

                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=16321*/
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx]; */   /*Z0311=16322*/
                            }   /*Z0311=16323*/
                            F12sez = F12sez+carr4p[nser]*qqn[nser]*pqsum;   /*Z0311=16324*/
                            delser = fabs((F12sez-oldF12sez)/F12sez);   /*Z0311=16325*/
                            if ( delser<0.0001 ) break;   /*Z0311=16326*/
                            oldF12sez = F12sez;   /*Z0311=16327*/
                        }   /*Z0311=16328*/
                        //101:   /*Z0311=16329*/
                        F12sum = F12sum+carr5p[ii]*carr5p[jj]*F12sez;   /*Z0311=16330*/
                    }   /*Z0311=16331*/
                }   /*Z0311=16332*/
                F12ser = F12sum/vvm;   /*Z0311=16333*/
                F12 = F12ser;   /*Z0311=16334*/
            }   /*Z0311=16335*/
            else
            {   /*Z0311=16336*/
                xrz = q*rad/(zr+1);   /*Z0311=16337*/
                arg = (zr+2*v+1)*atan(2*xrz);   /*Z0311=16338*/
                nen = pow(1+4*xrz*xrz,(zr+2*v+1)/2.0);   /*Z0311=16339*/
                arg1 = (zr+2*v)*atan(2*xrz);   /*Z0311=16340*/
                nen1 = pow(1+4*xrz*xrz,(zr+2*v)/2.0);   /*Z0311=16341*/
                arg2 = (zr+2*v-1)*atan(2*xrz);   /*Z0311=16342*/
                nen2 = pow(1+4*xrz*xrz,(zr+2*v-1)/2.0);   /*Z0311=16343*/

                F12asz = 0.0;   /*Z0311=16345*/
                for ( ii=1; ii<=inmax; ii++ )
                {   /*Z0311=16346*/
                    a1m = carr5p[ii]*pow(carr7p[ii],v);   /*  carr7p[ii]:=pp[ii]; */   /*Z0311=16347*/
                    for ( jj=1; jj<=inmax; jj++ )
                    {   /*Z0311=16348*/
                        a2m = carr5p[jj]*pow(carr7p[jj],v);   /*Z0311=16349*/
                        xijm = (carr3p[ii]-carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii]; */   /*Z0311=16350*/
                        arglmz = (zr+1)*atan(xijm);   /*Z0311=16351*/
                        nenlmz = pow(1+xijm*xijm,(zr+1)/2.0);   /*Z0311=16352*/
                        xijp = (carr3p[ii]+carr3p[jj])*q/(zr+1);   /*Z0311=16353*/
                        arglpz = (zr+1)*atan(xijp);   /*Z0311=16354*/
                        nenlpz = pow(1+xijp*xijp,(zr+1)/2.0);   /*Z0311=16355*/
                        /* F12as1z:=e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(pi*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(pi*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz)); */   /*Z0311=16356*/
                        F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(0-(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));   /*Z0311=16357*/
                        /* F12as2z:=e0*e1*(1/(carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(pi*(2*v-1)/2)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(pi*(2*v-1)/2)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz)); */   /*Z0311=16358*/
                        F12as2z = e0*e1*(1/(carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(1*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-0)/(nen1*nenlpz));   /*Z0311=16359*/
                        /* F12as3z:=e1*e0*(1/(carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(pi*(2*v-1)/2)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(pi*(2*v-1)/2)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz)); */   /*Z0311=16360*/
                        F12as3z = e1*e0*(1/(carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(1*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-0)/(nen1*nenlpz));   /*Z0311=16361*/
                        /* F12as4z:=e1*e1*(1/(carr7p[ii]*carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(pi*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(pi*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz)); */   /*Z0311=16362*/
                        F12as4z = e1*e1*(1/(carr7p[ii]*carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(0+1*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));   /*Z0311=16363*/

                        F12asz = F12asz+a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);   /*Z0311=16365*/
                    }   /*Z0311=16366*/
                }   /*Z0311=16367*/
                F12asy = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;   /*Z0311=16368*/
                F12 = F12asy;   /*Z0311=16369*/
            }   /*Z0311=16370*/
            return pql*F12;   /*Z0311=16371*/
            /* formpq:=pql; */   /*Z0311=16372*/
            /* formpq:=pql*polyliposome(llipt,radius,lliph,lin,lout,nom,sigmar,sigmal,phiax,philiph,philipt,phiin,phiout,2,q); */   /*Z0311=16374*/
            /* formpq:=pql*polyliposome(2.0,200,1.0,3.5,3.5,1,sigmar,sigmal,0.001,-0.55,-0.7,0.001,0.001,2,q); */   /*Z0311=16375*/
            /* formpq:=pql; */   /*Z0311=16376*/
        } /* of myelin */   /*Z0311=16377*/

    } /* of cylinder */   /*Z0311=16379*/


    /********/   /*Z0311=16383*/
    /* disk */   /*Z0311=16384*/
    /********/   /*Z0311=16385*/
    if ( part==2 )
    {   /*Z=16772*/
#ifndef __CUDACC__
        if ( dbgFlag() )
        {   // Kommt nur einmal pro Thread zu Beginn der Berechnungen
            qDebug() << "formpq: disk" << "ordis"<<ordis << "orcase"<<orcase << "q"<<q << "limq1"<<limq1 << "limq4"<<limq4
                     << "cs"<<cs << "length"<<length << "norm"<<norm << "limql"<<limql << "zl"<<zl;
        }
#endif
        /*** longitudinal part ***/
        /*** isotropic ***/
        if ( ordis==7 ) /*Z=16726*/
        {
            if ( q<(0.5*limq1) )
            {
                pqsum = 1.0;
                oldpqsum = 0.0;
                double qqn = 1.0;
                for ( nser=1; nser<=80; nser++ )
                {
                    qqn = qqn*q*q;
                    pqsum = pqsum+carr1p[nser]*qqn;
                    delser = fabs((pqsum-oldpqsum)/pqsum);
                    if ( delser<0.0001 ) break;
                    oldpqsum = pqsum;
                }
                //70:
                pql = pqsum;
            }
            else
            {   /*Z=16742*/
                arglq = q*length/(zl+1);
                pql = (2.0/(zl*(zl-1)))*pow(arglq,-2);
            }
        }  /* of isotropic */

        /* perfect */
        if ( ordis==6 ) /*Z=16748*/
        {
            if ( (limql*length) < 5 )
            {
                pqsum = 1.0;
                oldpqsum = 0.0;
                double qqn = 1.0;
                switch ( orcase )
                {
                case 1: /*Z=16753*/
                    argq = qxs+qys;
                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn = qqn*q*q;
                        pqsum = pqsum+carr1p[nser]*qqn*pow(1-argq*argq,nser);
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    break;
                case 2: /*Z=16763*/
                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn = qqn*q*q;
                        pqsum = pqsum+carr1p[nser]*qqn*pow(1-qxs*qxs,nser);
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    break;
                case 3: /*Z=16772*/
                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn = qqn*q*q;
                        pqsum = pqsum+carr1p[nser]*qqn*pow(1-qys*qys,nser);
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    break;
                case 4: /*Z=16781*/
                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn = qqn*q*q;
                        pqsum = pqsum+carr1p[nser]*qqn;
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    break;
                }   // switch orcase
                pql = pqsum;
            }
            else
            {   /*Z=16793*/
                if ( orcase==1 ) arglq = sqrt(1.0-sqr(qxs+qys))*q*length/(zl+1)+eps4;
                if ( orcase==2 ) arglq = sqrt(1.0-sqr(qxs))*q*length/(zl+1)+eps4;
                if ( orcase==3 ) arglq = sqrt(1.0-sqr(qys))*q*length/(zl+1)+eps4;
                if ( orcase==4 ) arglq = q*length/(zl+1)+eps4;

                pqr1 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3);
                pqr2 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3)*sin((zl-2)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-2)/2.0);
                pqr3 = (1/(zl*(zl-1)*(zl-2)*(zl-3)))*pow(arglq,-4)*cos((zl-3)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-3)/2.0);
                pql = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);
            }
        }   /* of perfect */

        /* orientational distribution */
        if ( ordis==0 ) /*Z=16810*/
        {
            if ( orcase==1 )
            {
                if ( q<(6*limq1) )
                {
                    pqsum = 1.0;
                    oldpqsum = 0.0;
                    qqn[0] = 1.0;
                    qxn[0] = 1.0;
                    qyn[0] = 1.0;

                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn[nser] = qqn[nser-1]*q*q;
                        qxn[nser] = qxn[nser-1]*qxs*qxs;
                        qyn[nser] = qyn[nser-1]*qys*qys;

                        binsum = 0.0;
                        for ( mser=0; mser<=nser; mser++ )
                        {
                            binsum1 = 0.0;
                            for ( lser=0; lser<=mser; lser++ )
                            {
                                /* indx:=lser+1+round(mser*(mser+1)/2); */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser]; */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];
                            }
                            /* indx:=mser+1+round(nser*(nser+1)/2); */
                            /* binsum:=binsum+carr1pm[indx]*binsum1; */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;
                        }
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*binsum;
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    pql = pqsum;
                }
                else
                {
                    /* disk: length = disk radius */
                    /* always use Bessel function approximation */
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);
                    pql = pql/norm;
                }
            } // if ( orcase==1 ) /*Z=16850*/

            if ( orcase==2 )
            {
                if ( q<(6000*limq1) )
                {
                    pqsum = 1.0;
                    oldpqsum = 0.0;
                    qqn[0] = 1.0;
                    qxn[0] = 1.0;
                    qyn[0] = 1.0;

                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn[nser] = qqn[nser-1]*q*q;
                        qxn[nser] = qxn[nser-1]*qxs*qxs;
                        qyn[nser] = qyn[nser-1]*qys*qys;

                        binsum = 0.0;
                        for ( mser=0; mser<=nser; mser++ )
                        {
                            binsum1 = 0.0;
                            for ( lser=0; lser<=mser; lser++ )
                            {
                                /* indx:=lser+1+round(mser*(mser+1)/2); */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser]; */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];
                            }
                            /* indx:=mser+1+round(nser*(nser+1)/2); */
                            /* binsum:=binsum+carr1pm[indx]*binsum1; */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;
                        }
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*binsum;
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    pql = pqsum;
                }
                else
                {
                    /* disk: length = disk radius */
                    /* always use Bessel function approximation */
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);
                    pql = pql/norm;
                    /* pql:=pq/1e-5; */
                    /* pql:=0.5; */
                }
            } // if ( orcase==2 )

            if ( orcase==3 )
            {
                if ( q<(6000*limq1) )
                {
                    pqsum = 1.0;
                    oldpqsum = 0.0;
                    qqn[0] = 1.0;
                    qxn[0] = 1.0;
                    qyn[0] = 1.0;

                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn[nser] = qqn[nser-1]*q*q;
                        qxn[nser] = qxn[nser-1]*qxs*qxs;
                        qyn[nser] = qyn[nser-1]*qys*qys;

                        binsum = 0.0;
                        for ( mser=0; mser<=nser; mser++ )
                        {
                            binsum1 = 0.0;
                            for ( lser=0; lser<=mser; lser++ )
                            {
                                /* indx:=lser+1+round(mser*(mser+1)/2); */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser]; */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];
                            }
                            /* indx:=mser+1+round(nser*(nser+1)/2); */
                            /* binsum:=binsum+carr1pm[indx]*binsum1; */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;
                        }
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*binsum;
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    //79:
                    pql = pqsum;
                }
                else
                {
                    /* disk: length = disk radius */
                    /* always use Bessel function approximation */
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);
                    pql = pql/norm;
                }
            }

            if ( orcase==4 )
            {
                if ( q<(6000*limq1) )
                {
                    pqsum = 1.0;
                    oldpqsum = 0.0;
                    qqn[0] = 1.0;
                    for ( nser=1; nser<=120; nser++ )
                    {
                        qqn[nser] = qqn[nser-1]*q*q;
                        pqsum = pqsum+carr1p[nser]*qqn[nser];
                        delser = fabs((pqsum-oldpqsum)/pqsum);
                        if ( delser<0.0001 ) break;
                        oldpqsum = pqsum;
                    }
                    //80:
                    pql = pqsum;
                }
                else
                {
                    /* disk: length = disk radius */
                    /* always use Bessel function approximation */
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);
                    pql = pql/norm;
                }
            }
        } // if ( ordis==0 )   /* of orientational distribution */

        /* transverse part */   /*Z0311=16625*/
        /* disk: radius = disk thickness/2 */   /*Z0311=16626*/
        /* homogeneous disk */   /*Z0311=16627*/
        if ( cs==0 )
        {   /*Z0311=16628*/
            if ( q<(0.5*limq4) )
            {   /*Z0311=16629*/
                pqsum = 1.0;   /*Z0311=16630*/
                oldpqsum = 0.0;   /*Z0311=16631*/
                qqn[0] = 1.0;   /*Z0311=16632*/
                for ( nser=1; nser<=100; nser++ )
                {   /*Z0311=16633*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=16634*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=16635*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=16636*/
                    if ( delser<0.0001 ) break;   /*Z0311=16637*/
                    oldpqsum = pqsum;   /*Z0311=16638*/
                }   /*Z0311=16639*/
                //71:   /*Z0311=16640*/
                pqr = pqsum;   /*Z0311=16641*/
            }   /*Z0311=16642*/
            else
            {   /*Z0311=16643*/
                argpq = q*radius/(zr+1);   /*Z0311=16644*/
                pqr = (1/(2*zr*(zr-1)))*pow(argpq,-2)*(1-cos((zr-1)*atan(2*argpq))/pow(1+4*argpq*argpq,(zr-1)/2.0));   /*Z0311=16645*/
            }   /*Z0311=16646*/
            return pql*pqr;   /*Z0311=16647*/
        } /* of homogeneous */   /*Z0311=16648*/

        /* core/shell disk */   /*Z0311=16650*/
        if ( cs==1 )
        {   /*Z0311=16651*/
            ccc1 = sqr(1-rho)*pow(p1,2);   /*Z0311=16652*/
            ccc2 = 2*rho*(1-rho)*pow(p1,1);   /*Z0311=16653*/
            ccc3 = rho*rho;   /*Z0311=16654*/
            vv3 = sqr((1-rho)*pow(p1,1)+rho);   /*Z0311=16655*/

            argq = q*radiusm/(zz+1);   /*Z0311=16657*/
            argpq = q*radius/(zz+1);   /*Z0311=16658*/

            /* F121 disk */   /*Z0311=16660*/
            if ( q<(0.8*limq4) )
            {   /*Z0311=16661*/
                /*** series expansion ***/   /*Z0311=16662*/
                pqsum = 1.0;   /*Z0311=16663*/
                oldpqsum = 0.0;   /*Z0311=16664*/
                qqn[0] = 1.0;   /*Z0311=16665*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=16666*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=16667*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=16668*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=16669*/
                    if ( delser<0.0001 ) break;   /*Z0311=16670*/
                    oldpqsum = pqsum;   /*Z0311=16671*/
                }   /*Z0311=16672*/
                //72:   /*Z0311=16673*/
                F121 = ccc1*pqsum/vv3;   /*Z0311=16674*/
            }   /*Z0311=16675*/
            else
            {   /*Z0311=16676*/
                pqr = (1/(2*zr*(zr-1)))*pow(argpq,-2)*(1-cos((zr-1)*atan(2*argpq))/pow(1+4*argpq*argpq,(zr-1)/2.0));   /*Z0311=16677*/
                F121 = ccc1*pqr/vv3;   /*Z0311=16678*/
            }   /*Z0311=16679*/

            /* F122 disk */   /*Z0311=16681*/
            if ( q<(2.0*limq5) )
            {   /*Z0311=16682*/
                /*** series expansion ***/   /*Z0311=16683*/
                pqsum = 1.0;   /*Z0311=16684*/
                oldpqsum = 0.0;   /*Z0311=16685*/
                qqn[0] = 1.0;   /*Z0311=16686*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=16687*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=16688*/
                    pqsum = pqsum+carr5p[nser]*qqn[nser];   /*Z0311=16689*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=16690*/
                    if ( delser<0.0001 ) break;   /*Z0311=16691*/
                    oldpqsum = pqsum;   /*Z0311=16692*/
                }   /*Z0311=16693*/
                //73:   /*Z0311=16694*/
                F122 = ccc2*pqsum/vv3;   /*Z0311=16695*/
            }   /*Z0311=16696*/
            else
            {   /*Z0311=16697*/
                argbm = (zr+1)*atan(argpq-argq);   /*Z0311=16698*/
                nenbm = pow(1+4*sqr(argpq-argq),(zr+1)/2.0);   /*Z0311=16699*/
                argbp = (zr+1)*atan(argpq+argq);   /*Z0311=16700*/
                nenbp = pow(1+4*sqr(argpq+argq),(zr+1)/2.0);   /*Z0311=16701*/

                pqr = (1/(2*zr*(zr-1)*(zr-2)*argpq*argq))*(cos(argbm)/nenbm-cos(argbp)/nenbp);   /*Z0311=16703*/
                F122 = ccc2*pqr/vv3;   /*Z0311=16704*/
            }   /*Z0311=16705*/

            /* F123 disk */   /*Z0311=16707*/
            if ( q<(0.3*limq6) )
            {   /*Z0311=16708*/
                /*** series expansion ***/   /*Z0311=16709*/
                pqsum = 1.0;   /*Z0311=16710*/
                oldpqsum = 0.0;   /*Z0311=16711*/
                qqn[0] = 1.0;   /*Z0311=16712*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=16713*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=16714*/
                    pqsum = pqsum+carr6p[nser]*qqn[nser];   /*Z0311=16715*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=16716*/
                    if ( delser<0.0001 ) break;   /*Z0311=16717*/
                    oldpqsum = pqsum;   /*Z0311=16718*/
                }   /*Z0311=16719*/
                //74:   /*Z0311=16720*/
                F123 = ccc3*pqsum/vv3;   /*Z0311=16721*/
            }   /*Z0311=16722*/
            else
            {   /*Z0311=16723*/
                pqr = (1/(2*zr*(zr-1)))*pow(argq,-2)*(1-cos((zr-1)*atan(2*argq))/pow(1+4*argq*argq,(zr-1)/2.0));   /*Z0311=16724*/
                F123 = ccc3*pqr/vv3;   /*Z0311=16725*/
                /* add more terms, if necessary */   /*Z0311=16726*/
            }   /*Z0311=16727*/
            return pql*(F121+F122+F123);   /*Z0311=16728*/
            /* formpq:=F122; */   /*Z0311=16729*/
        } /* of core/shell-disk */   /*Z0311=16730*/

        /* inhomogeneous core/shell disk */   /*Z0311=16732*/
        if ( cs==2 )
        {   /*Z0311=16733*/

            dim = 1;   /*Z0311=16735*/
            delc = 0.0001;   /*Z0311=16736*/
            xrad = q*radiusm;   /*Z0311=16737*/
            xradp = q*radius;   /*Z0311=16738*/
            x1z = q*radius/(2*(zr+1));   /*Z0311=16739*/
            x12z = x1z*x1z;   /*Z0311=16740*/
            x2z = q*radiusm/(2*(zr+1));   /*Z0311=16741*/
            x22z = x2z*x2z;   /*Z0311=16742*/

            lim = 18*exp(-5*sigmar);   /*Z0311=16744*/
            lim1 = lim;   /*Z0311=16745*/
            lim2 = lim*0.7;   /*Z0311=16746*/
            lim3 = lim;   /*Z0311=16747*/
            lim4 = lim;   /*Z0311=16748*/
            lim5 = lim*0.7;   /*Z0311=16749*/
            lim6 = lim*1.2;   /*Z0311=16750*/

            a1 = (dim-alfa)/2.0;   /*Z0311=16752*/
            b1 = dim/2.0;   /*Z0311=16753*/
            b2 = (dim+2-alfa)/2.0;   /*Z0311=16754*/
            b1s = (dim+2)/2.0;   /*Z0311=16755*/
            v = -b1s+1/2.0;   /*Z0311=16756*/
            c = a1-b1-b2+1/2.0;   /*Z0311=16757*/
            d0 = 1;   /*Z0311=16758*/
            //d1 = a1*(1+a1-b1)*(1+a1-b2);   /*Z0311=16759*/
            e0 = 1.0;   /*Z0311=16760*/
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;   /*Z0311=16761*/
            ee0 = 1.0;   /*Z0311=16762*/
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16*(1-b1s));   /*Z0311=16763*/

            gb1s = sqrt(M_PI)/2.0;   /*Z=17101 Neu Aug.2022*/
            pz2v = 1/(zr*(zr-1));   /*Z=17102*/
            pz2v1 = pz2v/(zr-2);   /*Z=17103*/
            pz2v2 = pz2v1/(zr-3);   /*Z=17104*/

            gz1 = gamma(zr+1);   /*Z=17106*/
            preg1 = gb1s/sqrt(M_PI);   /*Z=17107*/
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));   /*Z=17108*/
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));   /*Z=17109*/
            pzvc = gamma(zr+1+v+c)/gz1;   /*Z=17110*/
            pzvc1 = gamma(zr+1+v+c-1)/gz1;   /*Z=17111*/
            pzvc2 = gamma(zr+1+v+c-2)/gz1;   /*Z=17112*/
            pzac = gamma(zr+1-2*a1+c)/gz1;   /*Z=17113*/
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;   /*Z=17114*/
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;   /*Z=17115*/
            pzc = gamma(zr+1+2*c)/gz1;   /*Z=17116*/
            pzc1 = gamma(zr+1+2*c-1)/gz1;   /*Z=17117*/
            pza = gamma(zr+1-4*a1)/gz1;   /*Z=17118*/
            pzva = gamma(zr+1+v-2*a1)/gz1;   /*Z=17119*/
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z=17120*/
            //dnv0 = 1;   /*Z=17121*/
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;   /*Z=17122*/
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z=17123*/
            //pva0 = gamma(zr+1-4*a1)/gz1;   /*Z=17124*/

            cc1 = 1/(dim*dim);   /*Z=17126*/
            cc2 = 2*rho/(dim*(dim-alfa)*pow(p1,dim-alfa));   /*Z=17127*/
            cc3 = -2*rho/(dim*(dim-alfa));   /*Z=17128*/
            cc4 = rho*rho/((dim-alfa)*(dim-alfa)*pow(p1*p1,dim-alfa));   /*Z=17129*/
            cc5 = -2*rho*rho/((dim-alfa)*(dim-alfa)*pow(p1,dim-alfa));   /*Z=17130*/
            cc6 = rho*rho/((dim-alfa)*(dim-alfa));   /*Z=17131*/
            vv3 = cc1+cc2+cc3+cc4+cc5+cc6;   /*Z=17132*/

            /* term #1 series */   /*Z=17134*/
            if ( xradp<lim1 )
            {   /*Z=17135*/
                z12v[0] = 1;   /*Z=17136*/
                b1sv[0] = 1;   /*Z=17137*/
                fkv[0] = 1;   /*Z=17138*/
                qqn[0] = 1.0;   /*Z=17139*/
                F12sez = 1.0;   /*Z=17140*/
                oldF12sez = 0.0;   /*Z=17141*/
                for ( n=1; n<=120; n++ )
                {   /*Z=17142*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z=17143*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z=17144*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=17145*/
                    fkv[n] = fkv[n-1]*n;   /*Z=17146*/
                    //sum12[n] = 0;   /*Z=17147*/
                    /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]); */   /*Z=17148*/
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n]; */   /*Z=17149*/

                    F12sez = F12sez+carr4p[n]*qqn[n];   /*Z=17151*/

                    del = abs((F12sez-oldF12sez)/F12sez);   /*Z=17153*/
                    if ( del<delc ) break; //goto 221;   /*Z=17154*/
                    oldF12sez = F12sez;   /*Z=17155*/
                }   /*Z=17156*/
                //221:   /*Z=17157*/
                F12 = F12sez;   /*Z=17158*/
            }   /*Z=17159*/

            /* term #2 series */   /*Z=17161*/
            if ( xradp<lim2 )
            {   /*Z=17162*/
                z12v[0] = 1;   /*Z=17163*/
                a1v[0] = 1;   /*Z=17164*/
                b1v[0] = 1;   /*Z=17165*/
                b2v[0] = 1;   /*Z=17166*/
                b1sv[0] = 1;   /*Z=17167*/
                fkv[0] = 1;   /*Z=17168*/
                qqn[0] = 1.0;   /*Z=17169*/
                F22sez = 1.0;   /*Z=17170*/
                oldF22sez = 0.0;   /*Z=17171*/
                for ( n=1; n<=120; n++ )
                {   /*Z=17172*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z=17173*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z=17174*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z=17175*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z=17176*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z=17177*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=17178*/
                    fkv[n] = fkv[n-1]*n;   /*Z=17179*/
                    sum22[n] = 0;   /*Z=17180*/
                    /* for m:=0 to n do sum22[n]:=sum22[n]+a1v[n-m]*power(p1*p1,m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z=17181*/
                    /* F22sez:=F22sez+power(-x22z,n)*z12v[n]*sum22[n]; */   /*Z=17182*/

                    F22sez = F22sez+carr5p[n]*qqn[n];   /*Z=17184*/

                    del = abs((F22sez-oldF22sez)/F22sez);   /*Z=17186*/
                    if ( del<delc ) break; //goto 222;   /*Z=17187*/
                    oldF22sez = F22sez;   /*Z=17188*/
                }   /*Z=17189*/
                //222:   /*Z=17190*/
                F22 = F22sez;   /*Z=17191*/
            }   /*Z=17192*/

            /* term #3 series */   /*Z=17194*/
            if ( xradp<lim3 )
            {   /*Z=17195*/
                z12v[0] = 1;   /*Z=17196*/
                a1v[0] = 1;   /*Z=17197*/
                b1v[0] = 1;   /*Z=17198*/
                b2v[0] = 1;   /*Z=17199*/
                b1sv[0] = 1;   /*Z=17200*/
                fkv[0] = 1;   /*Z=17201*/
                qqn[0] = 1.0;   /*Z=17202*/
                F32sez = 1.0;   /*Z=17203*/
                oldF32sez = 0.0;   /*Z=17204*/
                for ( n=1; n<=120; n++ )
                {   /*Z=17205*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z=17206*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z=17207*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z=17208*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z=17209*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z=17210*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=17211*/
                    fkv[n] = fkv[n-1]*n;   /*Z=17212*/
                    sum32[n] = 0;   /*Z=17213*/
                    /* for m:=0 to n do sum32[n]:=sum32[n]+a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z=17214*/
                    /* F32sez:=F32sez+power(-x12z,n)*z12v[n]*sum32[n]; */   /*Z=17215*/

                    F32sez = F32sez+carr6p[n]*qqn[n];   /*Z=17217*/

                    del = abs((F32sez-oldF32sez)/F32sez);   /*Z=17219*/
                    if ( del<delc ) break; //goto 223;   /*Z=17220*/
                    oldF32sez = F32sez;   /*Z=17221*/
                }   /*Z=17222*/
                //223:   /*Z=17223*/
                F32 = F32sez;   /*Z=17224*/
            }   /*Z=17225*/

            /* term #4 series */   /*Z=17227*/
            if ( xradp<lim4 )
            {   /*Z=17228*/
                z12v[0] = 1;   /*Z=17229*/
                a1v[0] = 1;   /*Z=17230*/
                b1v[0] = 1;   /*Z=17231*/
                b2v[0] = 1;   /*Z=17232*/
                b1sv[0] = 1;   /*Z=17233*/
                fkv[0] = 1;   /*Z=17234*/
                qqn[0] = 1.0;   /*Z=17235*/
                F42sez = 1.0;   /*Z=17236*/
                oldF42sez = 0.0;   /*Z=17237*/
                for ( n=1; n<=120; n++ )
                {   /*Z=17238*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z=17239*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z=17240*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z=17241*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z=17242*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z=17243*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z=17244*/
                    fkv[n] = fkv[n-1]*n;   /*Z=17245*/
                    //sum42[n] = 0;   /*Z=17246*/
                    /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z=17247*/
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n]; */   /*Z=17248*/

                    F42sez = F42sez+carr7p[n]*qqn[n];   /*Z=17250*/

                    del = abs((F42sez-oldF42sez)/F42sez);   /*Z=17252*/
                    if ( del<delc ) break; //goto 224;   /*Z=17253*/
                    oldF42sez = F42sez;   /*Z=17254*/
                }   /*Z=17255*/
                //224:   /*Z=17256*/
                F42 = F42sez;   /*Z=17257*/
            }   /*Z=17258*/

            /* term #5 series */   /*Z0311=16924*/
            if ( (xradp)<lim5 )
            {   /*Z0311=16925*/
                z12v[0] = 1;   /*Z0311=16926*/
                a1v[0] = 1;   /*Z0311=16927*/
                b1v[0] = 1;   /*Z0311=16928*/
                b2v[0] = 1;   /*Z0311=16929*/
                b1sv[0] = 1;   /*Z0311=16930*/
                fkv[0] = 1;   /*Z0311=16931*/
                qqn[0] = 1.0;   /*Z0311=16932*/
                F52sez = 1.0;   /*Z0311=16933*/
                oldF52sez = 0.0;   /*Z0311=16934*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=16935*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=16936*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=16937*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=16938*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=16939*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=16940*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=16941*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=16942*/
                    //sum52[n] = 0;   /*Z0311=16943*/
                    /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=16944*/
                    /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n]; */   /*Z0311=16945*/

                    F52sez = F52sez+carr8p[n]*qqn[n];   /*Z0311=16947*/

                    del = fabs((F52sez-oldF52sez)/F52sez);   /*Z0311=16949*/
                    if ( del<delc ) break;   /*Z0311=16950*/
                    oldF52sez = F52sez;   /*Z0311=16951*/
                }   /*Z0311=16952*/
                //225:   /*Z0311=16953*/
                F52 = F52sez;   /*Z0311=16954*/
            }   /*Z0311=16955*/

            /* term #6 series */   /*Z0311=16957*/
            if ( (xradp)<lim6 )
            {   /*Z0311=16958*/
                z12v[0] = 1;   /*Z0311=16959*/
                a1v[0] = 1;   /*Z0311=16960*/
                b1v[0] = 1;   /*Z0311=16961*/
                b2v[0] = 1;   /*Z0311=16962*/
                b1sv[0] = 1;   /*Z0311=16963*/
                fkv[0] = 1;   /*Z0311=16964*/
                qqn[0] = 1.0;   /*Z0311=16965*/
                F62sez = 1.0;   /*Z0311=16966*/
                oldF62sez = 0.0;   /*Z0311=16967*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=16968*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=16969*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=16970*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=16971*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=16972*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=16973*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=16974*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=16975*/
                    //sum62[n] = 0;   /*Z0311=16976*/
                    /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]); */   /*Z0311=16977*/
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n]; */   /*Z0311=16978*/

                    F62sez = F62sez+carr9p[n]*qqn[n];   /*Z0311=16980*/

                    del = fabs((F62sez-oldF62sez)/F62sez);   /*Z0311=16982*/
                    if ( del<delc ) break;   /*Z0311=16983*/
                    oldF62sez = F62sez;   /*Z0311=16984*/
                }   /*Z0311=16985*/
                //226:   /*Z0311=16986*/
                F62 = F62sez;   /*Z0311=16987*/
            }   /*Z0311=16988*/

            /*** term #1 asymptote ***/   /*Z0311=16991*/
            if ( xradp>=lim1 )
            {   /*Z0311=16992*/
                arg11 = (zr+2*v+1)*atan(4*x1z);   /*Z0311=16993*/
                nen11 = pow(1+16*x1z*x1z,(zr+2*v+1)/2.0);   /*Z0311=16994*/
                arg12 = (zr+2*v)*atan(4*x1z);   /*Z0311=16995*/
                nen12 = pow(1+16*x1z*x1z,(zr+2*v)/2.0);   /*Z0311=16996*/
                arg13 = (zr+2*v-1)*atan(4*x1z);   /*Z0311=16997*/
                nen13 = pow(1+16*x1z*x1z,(zr+2*v-1)/2.0);   /*Z0311=16998*/

                F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);   /*Z0311=17000*/
                F12as2z = 2*ee0*ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);   /*Z0311=17001*/
                F12as3z = ee1*ee1*(1/(4*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);   /*Z0311=17002*/
                F12asz = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);   /*Z0311=17003*/
                F12 = F12asz;   /*Z0311=17004*/
            }   /*Z0311=17005*/

            /*** term #2 asymptote ***/   /*Z0311=17007*/
            if ( xradp>=lim2 )
            {   /*Z0311=17008*/
                //arg21 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=17009*/
                //nen21 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=17010*/
                //arg22 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=17011*/
                //nen22 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=17012*/
                //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);   /*Z0311=17013*/
                //F22as1sum2z = dnv0*ee1*(1/(2*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);   /*Z0311=17014*/
                F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);   /*Z0311=17015*/
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);   /*Z0311=17016*/

                arg210 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=17018*/
                nen210 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=17019*/
                arg220 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=17020*/
                nen220 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=17021*/
                F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);   /*Z0311=17022*/
                F22as1sum2z0 = ee1*(1/(2*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);   /*Z0311=17023*/
                F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);   /*Z0311=17024*/
                arg23 = (zr+v+c+1)*atan(2*(x1z-x2z));   /*Z0311=17025*/
                nen23 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);   /*Z0311=17026*/
                arg24 = (zr+v+c+1)*atan(2*(x1z+x2z));   /*Z0311=17027*/
                nen24 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);   /*Z0311=17028*/
                arg25 = (zr+v+c)*atan(2*(x1z-x2z));   /*Z0311=17029*/
                nen25 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);   /*Z0311=17030*/
                arg26 = (zr+v+c)*atan(2*(x1z+x2z));   /*Z0311=17031*/
                nen26 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);   /*Z0311=17032*/
                arg27 = (zr+v+c-1)*atan(2*(x1z-x2z));   /*Z0311=17033*/
                nen27 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);   /*Z0311=17034*/
                arg28 = (zr+v+c-1)*atan(2*(x1z+x2z));   /*Z0311=17035*/
                nen28 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);   /*Z0311=17036*/

                a22as21z = (1/2.0)*ee0*e0*pzvc;   /*Z0311=17038*/
                F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);   /*Z0311=17039*/
                a22as22z = (1/2.0)*ee0*e1*(1/(2*x2z))*pzvc1;   /*Z0311=17040*/
                F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);   /*Z0311=17041*/
                a22as23z = (1/2.0)*ee1*e0*(1/(2*x1z))*pzvc1;   /*Z0311=17042*/
                F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);   /*Z0311=17043*/
                a22as24z = (1/2.0)*ee1*e1*(1/(2*x1z))*(1/(2*x2z))*pzvc2;   /*Z0311=17044*/
                F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);   /*Z0311=17045*/
                F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);   /*Z0311=17046*/
                F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);   /*Z0311=17047*/
                //F22asz = F22as1z+F22as2z;   /*Z0311=17048*/
                F22asz0 = F22as1z0+F22as2z;   /*Z0311=17049*/
                F22 = F22asz0;   /*Z0311=17050*/
            }   /*Z0311=17051*/

            /*** term #3 asymptote ***/   /*Z0311=17053*/
            if ( xradp>=lim3 )
            {   /*Z0311=17054*/
                //arg31 = (zr+v-2*a1+1)*atan(2*x1z);   /*Z0311=17055*/
                //nen31 = pow(1+4*x1z*x1z,(zr+v-2*a1+1)/2.0);   /*Z0311=17056*/
                //arg32 = (zr+v-2*a1)*atan(2*x1z);   /*Z0311=17057*/
                //nen32 = pow(1+4*x1z*x1z,(zr+v-2*a1)/2.0);   /*Z0311=17058*/
                //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);   /*Z0311=17059*/
                //F32as1sum2z = dnv0*ee1*(1/(2*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);   /*Z0311=17060*/
                F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);   /*Z0311=17061*/
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);   /*Z0311=17062*/

                arg310 = (z+v-2*a1+1)*atan(2*x1z);   /*Z0311=17064*/
                nen310 = pow(1+4*x1z*x1z,(z+v-2*a1+1)/2.0);   /*Z0311=17065*/
                arg320 = (z+v-2*a1)*atan(2*x1z);   /*Z0311=17066*/
                nen320 = pow(1+4*x1z*x1z,(z+v-2*a1)/2.0);   /*Z0311=17067*/
                F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);   /*Z0311=17068*/
                F32as1sum2z0 = ee1*(1/(2*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);   /*Z0311=17069*/
                F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);   /*Z0311=17070*/

                arg33 = (zr+v+c+1)*atan(4*x1z);   /*Z0311=17072*/
                nen33 = pow(1+16*x1z*x1z,(zr+v+c+1)/2.0);   /*Z0311=17073*/
                arg34 = (zr+v+c)*atan(4*x1z);   /*Z0311=17074*/
                nen34 = pow(1+16*x1z*x1z,(zr+v+c)/2.0);   /*Z0311=17075*/
                arg35 = (zr+v+c-1)*atan(4*x1z);   /*Z0311=17076*/
                nen35 = pow(1+16*x1z*x1z,(zr+v+c-1)/2.0);   /*Z0311=17077*/
                F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);   /*Z0311=17078*/
                F32as22z = (1/2.0)*ee0*e1*(1/(2*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);   /*Z0311=17079*/
                F32as23z = (1/2.0)*ee1*e0*(1/(2*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);   /*Z0311=17080*/
                F32as24z = (1/2.0)*ee1*e1*(1/(4*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);   /*Z0311=17081*/
                F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);   /*Z0311=17082*/
                F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);   /*Z0311=17083*/
                //F32asz = F32as1z+F32as2z;   /*Z0311=17084*/
                F32asz0 = F32as1z0+F32as2z;   /*Z0311=17085*/
                F32 = F32asz0;   /*Z0311=17086*/
            }   /*Z0311=17087*/

            /*** term #4 asymptote ***/   /*Z0311=17090*/
            if ( xrad>=lim4 )
            {   /*Z0311=17091*/
                F42as10z = preg4*preg4*pow(x22z,-2*a1);   /*Z0311=17092*/
                //F42as1sumz = pva0;   /*Z0311=17093*/
                //F42as1z = F42as10z*F42as1sumz;   /*Z0311=17094*/
                F42as1z0 = F42as10z*pza;   /*Z0311=17095*/

                arg41 = (zr-2*a1+c+1)*atan(2*x2z);   /*Z0311=17097*/
                nen41 = pow(1+4*x2z*x2z,(zr-2*a1+c+1)/2.0);   /*Z0311=17098*/
                arg42 = (zr-2*a1+c)*atan(2*x2z);   /*Z0311=17099*/
                nen42 = pow(1+4*x2z*x2z,(zr-2*a1+c)/2.0);   /*Z0311=17100*/
                //arg43 = (zr-2*a1+c+3)*atan(2*x2z);   /*Z0311=17101*/
                //nen43 = pow(1+4*x2z*x2z,(zr-2*a1+c+3)/2.0);   /*Z0311=17102*/
                F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);   /*Z0311=17103*/
                F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);   /*Z0311=17104*/
                F42as22 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);   /*Z0311=17105*/
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);   /*Z0311=17106*/
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);   /*Z0311=17107*/
                F42as2z0 = F42as20z*(F42as21+F42as22);   /*Z0311=17108*/

                F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);   /*Z0311=17110*/
                F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);   /*Z0311=17111*/
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);   /*Z0311=17112*/
                F42as26 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);   /*Z0311=17113*/
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);   /*Z0311=17114*/
                F42as3z0 = F42as30z*(F42as24+F42as26);   /*Z0311=17115*/

                F42as40z = preg3*preg3*pow(x2z*x2z,c);   /*Z0311=17117*/
                arg44 = (zr+2*c+1)*atan(4*x2z);   /*Z0311=17118*/
                nen44 = pow(1+16*x2z*x2z,(zr+2*c+1)/2.0);   /*Z0311=17119*/
                arg45 = (zr+2*c)*atan(4*x2z);   /*Z0311=17120*/
                nen45 = pow(1+16*x2z*x2z,(zr+2*c)/2.0);   /*Z0311=17121*/
                F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);   /*Z0311=17122*/
                F42as28 = (1/2.0)*e0*e1*(1/(2*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=17123*/
                F42as29 = (1/2.0)*e1*e0*(1/(2*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=17124*/
                F42as4z = F42as40z*(F42as27+F42as28+F42as29);   /*Z0311=17125*/
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;   /*Z0311=17126*/
                F42asz0 = F42as1z0+F42as2z0+F42as3z0+F42as4z;   /*Z0311=17127*/
                F42 = F42asz0;   /*Z0311=17128*/
            }   /*Z0311=17129*/

            /*** term #5 asymptote ***/   /*Z0311=17132*/
            if ( xradp>=lim5 )
            {   /*Z0311=17133*/
                F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);   /*Z0311=17134*/
                //F52as1sumz = pva0;   /*Z0311=17135*/
                //F52as1z = F52as10z*F52as1sumz;   /*Z0311=17136*/
                F52as1z0 = F52as10z*pza;   /*Z0311=17137*/

                F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);   /*Z0311=17139*/
                arg51 = (zr-2*a1+c+1)*atan(2*x2z);   /*Z0311=17140*/
                nen51 = pow(1+4*x2z*x2z,(zr-2*a1+c+1)/2.0);   /*Z0311=17141*/
                arg52 = (zr-2*a1+c)*atan(2*x2z);   /*Z0311=17142*/
                nen52 = pow(1+4*x2z*x2z,(zr-2*a1+c)/2.0);   /*Z0311=17143*/
                //arg53 = (zr-2*a1+c+3)*atan(2*x2z);   /*Z0311=17144*/
                //nen53 = pow(1+4*x2z*x2z,(zr-2*a1+c+3)/2.0);   /*Z0311=17145*/
                F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);   /*Z0311=17146*/
                F52as22 = d0*e1*pzac1*(1/(2*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);   /*Z0311=17147*/
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);   /*Z0311=17148*/
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);   /*Z0311=17149*/
                F52as2z0 = F52as20z*(F52as21+F52as22);   /*Z0311=17150*/

                F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);   /*Z0311=17152*/
                arg54 = (zr-2*a1+c+1)*atan(2*x1z);   /*Z0311=17153*/
                nen54 = pow(1+4*x1z*x1z,(zr-2*a1+c+1)/2.0);   /*Z0311=17154*/
                //arg55 = (zr-2*a1+c+3)*atan(2*x1z);   /*Z0311=17155*/
                //nen55 = pow(1+4*x1z*x1z,(zr-2*a1+c+3)/2.0);   /*Z0311=17156*/
                arg56 = (zr-2*a1+c)*atan(2*x1z);   /*Z0311=17157*/
                nen56 = pow(1+4*x1z*x1z,(zr-2*a1+c)/2.0);   /*Z0311=17158*/
                F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);   /*Z0311=17159*/
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);   /*Z0311=17160*/
                F52as26 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);   /*Z0311=17161*/
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);   /*Z0311=17162*/
                F52as3z0 = F52as30z*(F52as24+F52as26);   /*Z0311=17163*/

                F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);   /*Z0311=17165*/
                arg57 = (zr+2*c+1)*atan(2*(x1z-x2z));   /*Z0311=17166*/
                nen57 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);   /*Z0311=17167*/
                arg58 = (zr+2*c+1)*atan(2*(x1z+x2z));   /*Z0311=17168*/
                nen58 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);   /*Z0311=17169*/
                arg59 = (zr+2*c)*atan(2*(x1z-x2z));   /*Z0311=17170*/
                nen59 = pow(1+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);   /*Z0311=17171*/
                arg510 = (zr+2*c)*atan(2*(x1z+x2z));   /*Z0311=17172*/
                nen510 = pow(1+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);   /*Z0311=17173*/
                F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);   /*Z0311=17174*/
                F52as28 = (1/2.0)*e0*e1*(1/(2*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);   /*Z0311=17175*/
                F52as29 = (1/2.0)*e1*e0*(1/(2*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);   /*Z0311=17176*/
                F52as4z = F52as40z*(F52as27+F52as28+F52as29);   /*Z0311=17177*/
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;   /*Z0311=17178*/
                F52asz0 = F52as1z0+F52as2z0+F52as3z0+F52as4z;   /*Z0311=17179*/
                F52 = F52asz0;   /*Z0311=17180*/
            }   /*Z0311=17181*/

            /*** term #6 asymptote ***/   /*Z0311=17183*/
            if ( xradp>=lim6 )
            {   /*Z0311=17184*/
                F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);   /*Z0311=17185*/
                //F62as1sumz = pva0;   /*Z0311=17186*/
                //F62as1z = F62as10z*F62as1sumz;   /*Z0311=17187*/
                F62as1z0 = F62as10z*pza;   /*Z0311=17188*/

                F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);   /*Z0311=17190*/
                arg61 = (zr-2*a1+c+1)*atan(2*x1z);   /*Z0311=17191*/
                nen61 = pow(1+4*x1z*x1z,(zr-2*a1+c+1)/2.0);   /*Z0311=17192*/
                arg62 = (zr-2*a1+c)*atan(2*x1z);   /*Z0311=17193*/
                nen62 = pow(1+4*x1z*x1z,(zr-2*a1+c)/2.0);   /*Z0311=17194*/
                //arg63 = (zr-2*a1+c+3)*atan(2*x1z);   /*Z0311=17195*/
                //nen63 = pow(1+4*x1z*x1z,(zr-2*a1+c+3)/2.0);   /*Z0311=17196*/
                F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);   /*Z0311=17197*/
                F62as22 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);   /*Z0311=17198*/
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);   /*Z0311=17199*/
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);   /*Z0311=17200*/
                F62as2z0 = F62as20z*(F62as21+F62as22);   /*Z0311=17201*/

                F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);   /*Z0311=17203*/
                F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);   /*Z0311=17204*/
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);   /*Z0311=17205*/
                F62as26 = d0*e1*pzac1*(1/(2*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);   /*Z0311=17206*/
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);   /*Z0311=17207*/
                F62as3z0 = F62as30z*(F62as24+F62as26);   /*Z0311=17208*/

                F62as40z = preg3*preg3*pow(x1z*x1z,c);   /*Z0311=17210*/
                arg64 = (zr+2*c+1)*atan(4*x1z);   /*Z0311=17211*/
                nen64 = pow(1+16*x1z*x1z,(zr+2*c+1)/2.0);   /*Z0311=17212*/
                arg65 = (zr+2*c)*atan(4*x1z);   /*Z0311=17213*/
                nen65 = pow(1+16*x1z*x1z,(zr+2*c)/2.0);   /*Z0311=17214*/
                F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);   /*Z0311=17215*/
                F62as28 = (1/2.0)*e0*e1*(1/(2*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=17216*/
                F62as29 = (1/2.0)*e1*e0*(1/(2*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=17217*/
                F62as4z = F62as40z*(F62as27+F62as28+F62as29);   /*Z0311=17218*/
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;   /*Z0311=17219*/
                F62asz0 = F62as1z0+F62as2z0+F62as3z0+F62as4z;   /*Z0311=17220*/
                F62 = F62asz0;   /*Z0311=17221*/
            }   /*Z0311=17222*/

            return pql*(cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;   /*Z0311=17224*/
            /* formpq:=pql*(cc5*F52)/vv3; */   /*Z0311=17225*/
            /* formpq:=pql*(F121+F122+F123); */   /*Z0311=17228*/
        } /* of inhomogeneous core/shell-disk */   /*Z0311=17230*/

        /* myelin disk */   /*Z0311=17232*/
        if ( (cs==3) || (cs==4) )
        {   /*Z0311=17233*/

            /* disk parameters */   /*Z0311=17235*/
            v = -1;   /*Z0311=17236*/
            e0 = 1;   /*Z0311=17237*/
            e1 = 0;   /*Z0311=17238*/
            preg1 = 1/2.0;   /*Z0311=17239*/
            pz2v = 1/(zr*(zr-1));   /*Z0311=17240*/
            pz2v1 = pz2v/(zr-2);   /*Z0311=17241*/
            pz2v2 = pz2v1/(zr-3);   /*Z0311=17242*/
            lim = 18*exp(-5*sigmar);   /*Z0311=17243*/
            lim1 = lim*1.2;   /*Z0311=17244*/
            double rad = myarray[1];   /*Z0311=17245 TODO*/
            inmax = round(myarray[14]);   /*Z0311=17246*/
            vvm = myarray[15];   /*Z0311=17247*/
            rmax = myarray[16];   /*Z0311=17248*/
            double xmax = q*rmax;   /*Z0311=17249 TODO*/

            if ( xmax<(lim1) )
            {   /*Z0311=17251*/
                /* fkv[0]:=1; */   /*Z0311=17252*/
                qqn[0] = 1.0;   /*Z0311=17253*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=17254*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17255*/
                    /* fkv[nser]:=fkv[nser-1]*nser; */   /*Z0311=17256*/
                }   /*Z0311=17257*/

                F12sum = 0.0;   /*Z0311=17259*/
                for ( ii=1; ii<=inmax; ii++ )
                {   /*Z0311=17260*/
                    for ( jj=1; jj<=inmax; jj++ )
                    {   /*Z0311=17261*/
                        F12sez = 1.0;   /*Z0311=17262*/
                        oldF12sez = 1.0;   /*Z0311=17263*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=17264*/
                            pqsum = 0;   /*Z0311=17265*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=17266*/
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]); */   /*Z0311=17267*/
                                pqsum = pqsum+pow(carr7p[ii],2*mser)*pow(carr7p[jj],2*(nser-mser))/(carr6p[mser]*carr6p[nser-mser]);   /*Z0311=17268*/
                                /*Z0311=17269*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=17270*/
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx]; */   /*Z0311=17271*/
                            }   /*Z0311=17272*/
                            F12sez = F12sez+carr4p[nser]*qqn[nser]*pqsum;   /*Z0311=17273*/
                            delser = fabs((F12sez-oldF12sez)/F12sez);   /*Z0311=17274*/
                            if ( delser<0.0001 ) break;   /*Z0311=17275*/
                            oldF12sez = F12sez;   /*Z0311=17276*/
                        }   /*Z0311=17277*/
                        //251:   /*Z0311=17278*/
                        F12sum = F12sum+carr5p[ii]*carr5p[jj]*F12sez;   /*Z0311=17279*/
                    }   /*Z0311=17280*/
                }   /*Z0311=17281*/
                F12ser = F12sum/vvm;   /*Z0311=17282*/
                F12 = F12ser;   /*Z0311=17283*/
            }   /*Z0311=17284*/
            else
            {   /*Z0311=17285*/
                xrz = q*rad/(zr+1);   /*Z0311=17286*/
                arg = (zr+2*v+1)*atan(2*xrz);   /*Z0311=17287*/
                nen = pow(1+4*xrz*xrz,(zr+2*v+1)/2.0);   /*Z0311=17288*/
                arg1 = (zr+2*v)*atan(2*xrz);   /*Z0311=17289*/
                nen1 = pow(1+4*xrz*xrz,(zr+2*v)/2.0);   /*Z0311=17290*/
                arg2 = (zr+2*v-1)*atan(2*xrz);   /*Z0311=17291*/
                nen2 = pow(1+4*xrz*xrz,(zr+2*v-1)/2.0);   /*Z0311=17292*/

                F12asz = 0.0;   /*Z0311=17294*/
                for ( ii=1; ii<=inmax; ii++ )
                {   /*Z0311=17295*/
                    a1m = carr5p[ii]*pow(carr7p[ii],v);   /*  carr7p[ii]:=pp[ii]; */   /*Z0311=17296*/
                    for ( jj=1; jj<=inmax; jj++ )
                    {   /*Z0311=17297*/
                        a2m = carr5p[jj]*pow(carr7p[jj],v);   /*Z0311=17298*/
                        xijm = (carr3p[ii]-carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii]; */   /*Z0311=17299*/
                        arglmz = (zr+1)*atan(xijm);   /*Z0311=17300*/
                        nenlmz = pow(1+xijm*xijm,(zr+1)/2.0);   /*Z0311=17301*/
                        xijp = (carr3p[ii]+carr3p[jj])*q/(zr+1);   /*Z0311=17302*/
                        arglpz = (zr+1)*atan(xijp);   /*Z0311=17303*/
                        nenlpz = pow(1+xijp*xijp,(zr+1)/2.0);   /*Z0311=17304*/
                        F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(M_PI*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(M_PI*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));   /*Z0311=17305*/
                        F12as2z = e0*e1*(1/(carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));   /*Z0311=17306*/
                        F12as3z = e1*e0*(1/(carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));   /*Z0311=17307*/
                        F12as4z = e1*e1*(1/(carr7p[ii]*carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(M_PI*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(M_PI*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));   /*Z0311=17308*/

                        F12asz = F12asz+a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);   /*Z0311=17310*/
                    }   /*Z0311=17311*/
                }   /*Z0311=17312*/
                F12asy = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;   /*Z0311=17313*/
                F12 = F12asy;   /*Z0311=17314*/
            }   /*Z0311=17315*/
            return pql*F12;   /*Z0311=17316*/
            /* formpq:=pql; */   /*Z0311=17317*/
            /* formpq:=pql*polyliposome(llipt,radius,lliph,lin,lout,nom,sigmar,sigmal,phiax,philiph,philipt,phiin,phiout,2,q); */   /*Z0311=17319*/
            /* formpq:=pql*polyliposome(2.0,200,1.0,3.5,3.5,1,sigmar,sigmal,0.001,-0.55,-0.7,0.001,0.001,1,q); */   /*Z0311=17320*/
            /* formpq:=pql; */   /*Z0311=17321*/
        } /* of myelin disk */   /*Z0311=17322*/

    } /* of disk */   /*Z0311=17324*/


    /* cube */   /*Z0311=17327*/
    if ( part==4 )
    {   /*Z0311=17328*/
        /* homogeneous isotropic cube */   /*Z0311=17329*/
        if ( ordis==7 )
        {   /*Z0311=17330*/
            if ( cs==0 )
            {   /*Z0311=17331*/
                if ( q<0.8*limq4 )
                {   /*Z0311=17332*/
                    pqsum = 1.0;   /*Z0311=17333*/
                    oldpqsum = 0.0;   /*Z0311=17334*/
                    qqn[0] = 1.0;   /*Z0311=17335*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=17336*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17337*/
                        pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=17338*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17339*/
                        if ( delser<0.0001 ) break;   /*Z0311=17340*/
                        oldpqsum = pqsum;   /*Z0311=17341*/
                    }   /*Z0311=17342*/
                    //81:   /*Z0311=17343*/
                    return pqsum;   /*Z0311=17344*/
                }   /*Z0311=17345*/
                else
                    return por/(q*q*q*q+eps4);   /*Z0311=17346*/
            } /* of homogeneous isotropic cube*/   /*Z0311=17347*/

            /* core/shell isotropic cube */   /*Z0311=17349*/
            if ( cs==1 )
            {   /*Z0311=17350*/
                return polycscube(1.0,rho,p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);   /*Z0311=17351*/
            }   /*Z0311=17352*/
        }  /* of isotropic cube */   /*Z0311=17353*/

        /* perfectly oriented cube */   /*Z0311=17355*/
        if ( ordis==6 )
        {   /*Z0311=17356*/
            /* if (orcase=1) then begin */   /*Z0311=17357*/
            if ( 1==1 )     // TODO ???
            {   /*Z0311=17358*/

                if ( q<(1/radius) )
                {   /*Z0311=17360*/
                    pqsum = 1.0;   /*Z0311=17361*/
                    oldpqsum = 0.0;   /*Z0311=17362*/
                    qxn[0] = 1.0;   /*Z0311=17363*/
                    qyn[0] = 1.0;   /*Z0311=17364*/
                    for ( nser=1; nser<=80; nser++ )
                    {   /*Z0311=17365*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=17366*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=17367*/
                        /*Z0311=17368*/
                        binsum = 0.0;   /*Z0311=17369*/
                        for ( mser=0; mser<=nser; mser++ )
                        {   /*Z0311=17370*/
                            /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=17371*/
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser]; */   /*Z0311=17372*/
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];   /*Z0311=17373*/
                        }   /*Z0311=17374*/
                        pqsum = pqsum+binsum;   /*Z0311=17375*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17376*/
                        if ( delser<0.0001 ) break;   /*Z0311=17377*/
                        oldpqsum = pqsum;   /*Z0311=17378*/
                    }   /*Z0311=17379*/
                    //84:   /*Z0311=17380*/
                    pql = pqsum;   /*Z0311=17381*/
                }   /*Z0311=17382*/
                else
                {   /*Z0311=17383*/
                    argqx = qxs*radius/(zr+1)+eps4;   /*Z0311=17384*/
                    argqy = qys*radius/(zr+1)+eps4;   /*Z0311=17385*/
                    pqrx = (1/(2*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2*argqx))/pow(1+4*argqx*argqx,(zr-1)/2.0));   /*Z0311=17386*/
                    pqry = (1/(2*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2*argqy))/pow(1+4*argqy*argqy,(zr-1)/2.0));   /*Z0311=17387*/
                    pql = pqrx*pqry;   /*Z0311=17388*/
                }   /*Z0311=17389*/
                return pql;   /*Z0311=17390*/
            }  /* of orcase=1 */   /*Z0311=17391*/
        }  /* of perfect cube */   /*Z0311=17392*/
    }  /* of cube */   /*Z0311=17393*/


    /*** biaxial ellipsoid ***/   /*Z=17732*/
    if ( part==5 )
    {   /*Z0311=17397*/

#ifndef __CUDACC__
        if ( dbgFlag() )
        {   // Kommt nur einmal pro Thread zu Beginn der Berechnungen
            qDebug() << "formpq: part=5" << "ordis"<<ordis << "orcase"<<orcase << "q"<<q << "limq1"<<limq1 << "limq4"<<limq4
                     << "cs"<<cs << "len/rad"<<length/radius << "norm"<<norm << "limql"<<limql << "q?s"<<qxs<<qys;
        }
#endif

        /* homogeneous isotropic ellipsoid */   /*Z0311=17398*/
        if ( ordis==7 )
        {   /*Z0311=17399*/
            if ( cs==0 )
            {   /*Z0311=17400*/
                if ( q<0.8*limq4 )      // früher: 0.8, im Code vom 4.10.22 steht 0.8
                {   /*Z0311=17401*/
                    pqsum = 1.0;   /*Z0311=17402*/
                    oldpqsum = 0.0;   /*Z0311=17403*/
                    qqn[0] = 1.0;   /*Z0311=17404*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=17405*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17406*/
                        pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=17407*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17408*/
                        if ( delser<0.0001 ) break;   /*Z0311=17409*/
                        oldpqsum = pqsum;   /*Z0311=17410*/
                    }   /*Z0311=17411*/
                    //260:   /*Z0311=17412*/
                    return pqsum;   /*Z0311=17413*/
                }   /*Z0311=17414*/
                else
                {   /*Z0311=17415*/
                    if ( (q>1.5*limq4) )  // früher: 2, im Code vom 4.10.22 steht 1.5
                        return por/(q*q*q*q);   /*Z0311=17416*/
                    else
                    {   /*Z0311=17417*/
                        qrombchid(length,radius,/*p1,*/params.sigma,params.dbeta,epsi,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qx,qy,0,qhkl,ax1.length(),ax2.length(),ax3.length(),ax1.x(),ax1.y(),ax1.z(),ax2.x(),ax2.y(),ax2.z(),ax3.x(),ax3.y(),ax3.z(),sig.x(),sig.y(),sig.z(),ordis,3,8,13,7,0,0,carr1p,pql);   /*Z0311=17418*/
                        return pql;   /*Z0311=17419*/
                    }   /*Z0311=17420*/
                    //return pq;   /*Z0311=17421*/
                }   /*Z0311=17422*/
            } /* of homogeneous isotropic ellipsoid*/   /*Z0311=17423*/

            /* core/shell isotropic ellipsoid */   /*Z0311=17425*/
            if ( cs==1 )
            {   /*Z0311=17426*/
                return polycscube(1.0,rho,p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);   /*Z0311=17427*/
            }   /*Z0311=17428*/
        }  /* of isotropic ellipsoid */   /*Z0311=17429*/

        /* perfectly oriented ellipsoid */   /*Z0311=17431*/
        if ( ordis==6 )
        {   /*Z0311=17432*/
            /* if (orcase=1) then begin */   /*Z0311=17433*/
            if ( 1==1 )     // TODO ???
            {   /*Z0311=17434*/

                if ( sqrt(qx*qx*length*length+qy*qy*radius*radius+eps4)<15 )
                {   /*Z0311=17436*/
                    qxn[0] = 1.0;   /*Z0311=17437*/
                    qyn[0] = 1.0;   /*Z0311=17438*/
                    for ( nser=1; nser<=81; nser++ )
                    {   /*Z0311=17439*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=17440*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=17441*/
                    }   /*Z0311=17442*/
                    pqsum = 0.0;   /*Z0311=17443*/
                    oldpqsum = 0.0;   /*Z0311=17444*/
                    for ( nser=0; nser<=80; nser++ )
                    {   /*Z0311=17445*/
                        binsum = 0.0;   /*Z0311=17446*/
                        for ( mser=0; mser<=80; mser++ ) binsum = binsum+carr5p[mser]*params.CR->carr11pm[nser][mser]*qyn[mser];   /*Z0311=17447*/
                        pqsum = pqsum+carr4p[nser]*qxn[nser]*binsum;   /*Z0311=17448*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17449*/
                        if ( delser<0.0001 ) break;   /*Z0311=17450*/
                        oldpqsum = pqsum;   /*Z0311=17451*/
                    }   /*Z0311=17452*/
                    //261:   /*Z0311=17453*/
                    pql = pqsum;   /*Z0311=17454*/
                }   /*Z0311=17455*/
                else
                {   /*Z0311=17456*/
                    a1 = 9*pow(zr+1,4)/(2*zr*(zr-1)*(zr-2)*(zr-3));   /*Z0311=17457*/
                    pql = a1/(sqr(qy*qy*radius*radius+qx*qx*length*length)+eps4);   /*Z0311=17458*/
                }   /*Z0311=17459*/
                return pql;   /*Z0311=17460*/
            }  /* of orcase=1 */   /*Z0311=17461*/
        }  /* of perfect ellipsoid */   /*Z0311=17462*/


        /* general */   /*Z0311=17465*/
        if ( ordis==0 )
        {   /*Z0311=17466*/
            if ( orcase==4 )
                pql = 1.0;   /*Z0311=17467*/
            else
            {   /*Z0311=17468*/
                if ( sqrt(qx*qx*length*length+qy*qy*radius*radius+eps4)<10 )
                {   /*Z0311=17469*/
                    pqsum = 0.0;   /*Z0311=17470*/
                    oldpqsum = -10.0;   /*Z0311=17471*/
                    qxn[0] = 1.0;   /*Z0311=17472*/
                    qyn[0] = 1.0;   /*Z0311=17473*/
                    qqn[0] = 1.0;   /*Z0311=17474*/

                    if ( orcase==1 )
                    {   /*Z0311=17476*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=17477*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs/(q*q);   /*Z0311=17478*/
                            qyn[nser] = qyn[nser-1]*qys*qys/(q*q);   /*Z0311=17479*/
                            qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17480*/
                        }   /*Z0311=17481*/
                        for ( nser=0; nser<=120; nser++ )
                        {   /*Z0311=17482*/
                            binsum1 = 0.0;   /*Z0311=17483*/
                            for ( lser=0; lser<=nser; lser++ )   /*Z0311=17484*/
                                binsum1 = binsum1+carr1p[lser]*qxn[lser]*qyn[nser-lser];   /*Z0311=17485*/
                            binsum = 0.0;   /*Z0311=17486*/
                            for ( mser=0; mser<=120; mser++ )   /*Z0311=17487*/
                                binsum = binsum+carr2p[mser]*qqn[mser]*params.CR->carr11pm[nser][mser];   /*Z0311=17488*/
                            pqsum = pqsum+carr3p[nser]*qqn[nser]*binsum*binsum1/pow(4.0,nser);   /*Z0311=17489 TODO hier war vorher pow(4,n)*/
                            // Da n aber eine globale Schleifenvariable ist, macht das hier keinen Sinn.
                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17490*/
                            if ( delser<0.0001 ) break;   /*Z0311=17491*/
                            oldpqsum = pqsum;   /*Z0311=17492*/
                        }   /*Z0311=17493*/
                    }   /*Z0311=17494*/

                    if ( orcase==2 )
                    {  /* x-axis */   /*Z0311=17496*/
                        for ( nser=1; nser<=110; nser++ )
                        {   /*Z0311=17497*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs/(q*q);   /*Z0311=17498*/
                            qyn[nser] = qyn[nser-1]*qys*qys/(q*q);   /*Z0311=17499*/
                            qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17500*/
                        }   /*Z0311=17501*/
                        for ( nser=0; nser<=100; nser++ )
                        {   /*Z0311=17502*/
                            binsum1 = 0.0;   /*Z0311=17503*/
                            for ( lser=0; lser<=nser; lser++ )   /*Z0311=17504*/
                                /* binsum1:=binsum1+carr4p[lser]*qxn[lser]*qyn[nser-lser]; */   /*Z0311=17505*/
                                binsum1 = binsum1+params.CR->carr22pm[nser][lser]*qxn[lser]*qyn[nser-lser];   /*Z0311=17506*/
                            binsum = 0.0;   /*Z0311=17507*/
                            for ( mser=0; mser<=100; mser++ )   /*Z0311=17508*/
                                binsum = binsum+carr5p[mser]*qqn[mser]*params.CR->carr11pm[nser][mser];   /*Z0311=17509*/
                            pqsum = pqsum+carr6p[nser]*qqn[nser]*binsum*binsum1/pow(4,nser);   /*Z0311=17510*/
                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17511*/
                            if ( delser<0.0001 ) break;   /*Z0311=17512*/
                            oldpqsum = pqsum;   /*Z0311=17513*/
                        }   /*Z0311=17514*/
                    }   /*Z0311=17515*/

                    if ( orcase==3 )
                    {  /* y-axis */   /*Z0311=17517*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=17518*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=17519*/
                            qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=17520*/
                            binsum = 0.0;   /*Z0311=17521*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=17522*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=17523*/
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser]; */   /*Z0311=17524*/
                                binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qyn[mser]*qxn[nser-mser];   /*Z0311=17525*/
                            }   /*Z0311=17526*/
                            pqsum = pqsum+carr1p[nser]*binsum;   /*Z0311=17527*/
                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17528*/
                            if ( delser<0.0001 ) break;   /*Z0311=17529*/
                            oldpqsum = pqsum;   /*Z0311=17530*/
                        }   /*Z0311=17531*/
                    }   /*Z0311=17532*/
                    //273:   /*Z0311=17533*/
                    pql = pqsum;   /*Z0311=17534*/
                }   /* of q<lim */   /*Z0311=17535*/
                else
                {   /*Z0311=17536*/
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,4,16,0,0,0,carr1p,pql);   /*Z0311=17537*/
                    pql = pql/norm;   /*Z0311=17538*/
                }   /*Z0311=17539*/
            }  /* of orcase=1,2,3 */   /*Z0311=17540*/
            return pql;   /*Z0311=17541*/
        }   /* of general */   /*Z0311=17542*/
    }  /* of biaxial ellipsoid */   /*Z0311=17543*/


    /*** triaxial ellipsoid ***/   /*Z0311=17546*/
    if ( part==6 )
    {   /*Z0311=17547*/
        /* homogeneous isotropic triaxial ellipsoid */   /*Z0311=17548*/
        if ( ordis==7 )
        {   /*Z0311=17549*/
            if ( cs==0 )
            {   /*Z0311=17550*/
                if ( q<0.05*limq4 )
                {   /*Z0311=17551*/
                    pqsum = 1.0;   /*Z0311=17552*/
                    oldpqsum = 0.0;   /*Z0311=17553*/
                    qqn[0] = 1.0;   /*Z0311=17554*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=17555*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17556*/
                        pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=17557*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17558*/
                        if ( delser<0.0001 ) break;   /*Z0311=17559*/
                        oldpqsum = pqsum;   /*Z0311=17560*/
                    }   /*Z0311=17561*/
                    //263:   /*Z0311=17562*/
                    return pqsum;   /*Z0311=17563*/
                }   /*Z0311=17564*/
                else
                {   /*Z0311=17565*/
                    if ( q>2*limq4 )
                        return por/(q*q*q*q);   /*Z0311=17566*/
                    else
                    {   /*Z0311=17567*/
                        qrombdeltac(length,radius,/*p1,sigma,dbeta,*/theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,7,14,7,0,0,carr1p,pql);   /*Z0311=17568*/
                        return pql/(M_PI/2.0);   /*Z0311=17569*/
                    }   /*Z0311=17570*/
                    //return pq;   /*Z0311=17571*/
                }   /*Z0311=17572*/
            } /* of homogeneous isotropic ellipsoid*/   /*Z0311=17573*/

            /* core/shell isotropic ellipsoid */   /*Z0311=17575*/
            if ( cs==1 )
            {   /*Z0311=17576*/
                return polycscube(1.0,rho,p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);   /*Z0311=17577*/
            }   /*Z0311=17578*/
        }  /* of isotropic triaxial ellipsoid */   /*Z0311=17579*/

        /* perfectly oriented triaxial ellipsoid */   /*Z0311=17581*/
        if ( ordis==6 )
        {   /*Z0311=17582*/
            /* if (orcase=1) then begin */   /*Z0311=17583*/
            if ( 1==1 )     // TODO ???
            {   /*Z0311=17584*/

                if ( q<(1/radius) )
                {   /*Z0311=17586*/
                    pqsum = 1.0;   /*Z0311=17587*/
                    oldpqsum = 0.0;   /*Z0311=17588*/
                    qxn[0] = 1.0;   /*Z0311=17589*/
                    qyn[0] = 1.0;   /*Z0311=17590*/
                    for ( nser=1; nser<=80; nser++ )
                    {   /*Z0311=17591*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=17592*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=17593*/
                        /*Z0311=17594*/
                        binsum = 0.0;   /*Z0311=17595*/
                        for ( mser=0; mser<=nser; mser++ )
                        {   /*Z0311=17596*/
                            /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=17597*/
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser]; */   /*Z0311=17598*/
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];   /*Z0311=17599*/
                        }   /*Z0311=17600*/
                        pqsum = pqsum+binsum;   /*Z0311=17601*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17602*/
                        if ( delser<0.0001 ) break;   /*Z0311=17603*/
                        oldpqsum = pqsum;   /*Z0311=17604*/
                    }   /*Z0311=17605*/
                    //264:   /*Z0311=17606*/
                    pql = pqsum;   /*Z0311=17607*/
                }   /*Z0311=17608*/
                else
                {   /*Z0311=17609*/
                    argqx = qxs*radius/(zr+1)+eps4;   /*Z0311=17610*/
                    argqy = qys*radius/(zr+1)+eps4;   /*Z0311=17611*/
                    pqrx = (1/(2*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2*argqx))/pow(1+4*argqx*argqx,(zr-1)/2.0));   /*Z0311=17612*/
                    pqry = (1/(2*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2*argqy))/pow(1+4*argqy*argqy,(zr-1)/2.0));   /*Z0311=17613*/
                    pql = pqrx*pqry;   /*Z0311=17614*/
                }   /*Z0311=17615*/
                return pql;   /*Z0311=17616*/
            }  /* of orcase=1 */   /*Z0311=17617*/
        }  /* of perfect triaxial ellipsoid */   /*Z0311=17618*/
    }  /* of triaxial ellipsoid */   /*Z0311=17619*/


    /*** super ellipsoid, barrel ***/   /*Z0311=17622*/
    if ( part==7 )
    {   /*Z0311=17623*/
        /* homogeneous isotropic super ellipsoid */   /*Z0311=17624*/
        if ( ordis==7 )
        {   /*Z0311=17625*/
            if ( cs==0 )
            {   /*Z0311=17626*/
                if ( q<1.5*limq4 )
                {   /*Z0311=17627*/
                    pqsum = 1.0;   /*Z0311=17628*/
                    oldpqsum = 0.0;   /*Z0311=17629*/
                    qqn[0] = 1.0;   /*Z0311=17630*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=17631*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17632*/
                        pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=17633*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17634*/
                        if ( delser<0.0001 ) break;   /*Z0311=17635*/
                        oldpqsum = pqsum;   /*Z0311=17636*/
                    }   /*Z0311=17637*/
                    //270:   /*Z0311=17638*/
                    return pqsum;   /*Z0311=17639*/
                }   /*Z0311=17640*/
                else
                {   /*Z0311=17641*/
                    if ( q>=1.5*limq4 )
                        return por/(q*q*q*q);   /*Z0311=17642*/
                    else
                    {   /*Z0311=17643*/
                        qrombdeltac(length,radius,/*p1,sigma,dbeta,*/theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,7,14,7,0,0,carr1p,pql);   /*Z0311=17644*/
                        return pql/(M_PI/2.0);   /*Z0311=17645*/
                    }   /*Z0311=17646*/
                    //return pq;   /*Z0311=17647*/
                }   /*Z0311=17648*/
            } /* of homogeneous isotropic ellipsoid*/   /*Z0311=17649*/

            /* core/shell isotropic ellipsoid */   /*Z0311=17651*/
            if ( cs==1 )
            {   /*Z0311=17652*/
                return polycscube(1.0,rho,p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);   /*Z0311=17653*/
            }   /*Z0311=17654*/
        }  /* of isotropic triaxial ellipsoid */   /*Z0311=17655*/

        /* perfectly oriented super ellipsoid */   /*Z0311=17657*/
        if ( ordis==6 )
        {   /*Z0311=17658*/
            /* if (orcase=1) then begin */   /*Z0311=17659*/
            if ( 1==1 )     // TODO ???
            {   /*Z0311=17660*/
                /*Z0311=17661*/
                if ( q<(1/radius) )
                {   /*Z0311=17662*/
                    pqsum = 1.0;   /*Z0311=17663*/
                    oldpqsum = 0.0;   /*Z0311=17664*/
                    qxn[0] = 1.0;   /*Z0311=17665*/
                    qyn[0] = 1.0;   /*Z0311=17666*/
                    for ( nser=1; nser<=80; nser++ )
                    {   /*Z0311=17667*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=17668*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=17669*/

                        binsum = 0.0;   /*Z0311=17671*/
                        for ( mser=0; mser<=nser; mser++ )
                        {   /*Z0311=17672*/
                            /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=17673*/
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser]; */   /*Z0311=17674*/
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];   /*Z0311=17675*/
                        }   /*Z0311=17676*/
                        pqsum = pqsum+binsum;   /*Z0311=17677*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17678*/
                        if ( delser<0.0001 ) break;   /*Z0311=17679*/
                        oldpqsum = pqsum;   /*Z0311=17680*/
                    }   /*Z0311=17681*/
                    //271:   /*Z0311=17682*/
                    pql = pqsum;   /*Z0311=17683*/
                }   /*Z0311=17684*/
                else
                {   /*Z0311=17685*/
                    argqx = qxs*radius/(zr+1)+eps4;   /*Z0311=17686*/
                    argqy = qys*radius/(zr+1)+eps4;   /*Z0311=17687*/
                    pqrx = (1/(2*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2*argqx))/pow(1+4*argqx*argqx,(zr-1)/2.0));   /*Z0311=17688*/
                    pqry = (1/(2*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2*argqy))/pow(1+4*argqy*argqy,(zr-1)/2.0));   /*Z0311=17689*/
                    pql = pqrx*pqry;   /*Z0311=17690*/
                }   /*Z0311=17691*/
                return pql;   /*Z0311=17692*/
            }  /* of orcase=1 */   /*Z0311=17693*/
        }  /* of perfect triaxial ellipsoid */   /*Z0311=17694*/
    }  /* of super ellipsoid */   /*Z0311=17695*/
    //qDebug() << "formpq default return";
    return 0; // kommt das Programm jemals hierhin?
}   /*Z0311=17697*/


#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::formfq( double length, double radius, double sigmal, double sigmar, double p1, double rho,
                         double alfa, double theta, double phi, double /*limql*/, double limq1, double /*limq2*/,
                         double /*limq3*/, double limq4, double limq5, double limq6, double qx, double qy,
                         double qxs, double qys, double q, double norm,
                         int part, int cs, int ordis, int orcase,
                         const double */*myarray*/, // CoeffArrayType
                         double *carr1p, double *carr2p, double */*carr3p*/, double *carr4p, double *carr5p,
                         double *carr6p, double */*carr7p*/ //: CoeffArrayType;   /*Z0311=17704*/,
                         /*ArrayImax2D &carr11pm, ArrayImax2D &carr22pm*/ ) const //: ArrayImax2D);   /*Z0311=17705*/
{
    // label 50,51,52,53,60,61,62,63,64,65,66,70,71,72,73,74,75,76,77,78,79,80,81,82,83;   /*Z0311=17707*/
    // label 101,102,103,104,105,106,111,112,113,114,115,116,121,122,123,124,125,126;   /*Z0311=17708*/

    int    n,/*m,*/nser,mser,lser;//,indx;
    double pqsum,oldpqsum,binsum,delser,argq,arglq,/*argp1,*/pqr,pqr1,pqr2,pql,pq1,pq2,pq3,pq4,pq5,pq6;
    double ccc1,ccc2,ccc3,vv3,zl,zr,radiusm;
    double cc1,cc4,cc6; //,cc2,cc3,cc5,cc7,cc8,cc9,cc10;
    //double ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10;
    //double argbm,nenbm,argbp,nenbp,argep,nenep,argem,nenem,arggp,nengp;
    double /*arggm,nengm,argim,nenim,argip,nenip,*/F121,F122,F123;
    double qqn[200+1],/*qqnx[200+1],qqny[200+1],*/z12v[200+1],a1v[200+1],b1v[200+1],b2v[200+1],b1sv[200+1],fkv[200+1]; //,gam3[200+1]; //: array[0..200] of extended;   /*Z0311=17717*/;
    double dim,xrad,xradp,x1z,x12z,x2z,x22z,lim,lim1,lim4,lim6; //,lim2,lim3,lim5
    double a1,b1,b2,b1s,v,c,/*d0,d1,*/e0,e1,ee0,ee1;
    double gb1s,pz2v,pz2v1,/*pz2v2,*/gz1,preg1,preg3,preg4,/*pzvc,*/pzvc1,pzvc2,pzac,pzac1,pzac2;
    double pzc,pzc1,pza,pzva,pzva1,dnv0,pvav0,pvav10,pva0,sumc;
    double del,delc,F12,F12sez,oldF12sez,F42,F42sez,oldF42sez,F62,F62sez,oldF62sez;
    double arg11,nen11,arg12,nen12,F12as1z,F12as2z,F12asz;
    double arg44,nen44,arg45,nen45,F42as10z,F42as1sumz,F42as1z,F42as1z0,F42as40z,F42as27,F42as28,F42as4z,F42asz,F42asz0;
    double arg64,nen64,arg65,nen65,F62as10z,F62as1z,F62as1z0,F62as40z,F62as1sumz,F62as27,F62as28,F62as4z,F62asz,F62asz0,FF1;

    // Im Orginal globale Variablen
    double argpq, zz;
    double qxn[200+1], qyn[200+1];  // TODO fehlende (bzw. globale) Variablen
    double pqr3, pqr4, binsum1;
    const double qz = 1; // TODO
    //const double z = 1; // TODO

    zl = (1-sigmal*sigmal)/(sigmal*sigmal);   /*Z0311=17728*/
    zr = (1-sigmar*sigmar)/(sigmar*sigmar);   /*Z0311=17729*/
    radiusm = radius/p1;   /* outer radius of core/shell particle */   /*Z0311=17730*/
    zz = zr;   /*Z0311=17855, 'zz' wird an anderer Stelle uninitialisiert verwendet */

    /**************/   /*Z0311=17732*/
    /*** sphere ***/   /*Z0311=17733*/
    /**************/   /*Z0311=17734*/
    if ( part==0 )
    {   /*Z0311=17735*/
        /*** homogeneous sphere ***/   /*Z0311=17736*/
        if ( cs==0 )
        {   /*Z0311=17737*/
            if ( q<(0.6*limq4) )
            {   /*Z0311=17738*/
                pqsum = 1.0;   /*Z0311=17739*/
                oldpqsum = 0.0;   /*Z0311=17740*/
                qqn[0] = 1.0;   /*Z0311=17741*/
                for ( nser=1; nser<=100; nser++ )
                {   /*Z0311=17742*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17743*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=17744*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17745*/
                    if ( delser<0.0001 ) break;   /*Z0311=17746*/
                    oldpqsum = pqsum;   /*Z0311=17747*/
                }   /*Z0311=17748*/
                //50:   /*Z0311=17749*/
                return pqsum;   /*Z0311=17750*/
            }   /*Z0311=17751*/
            else
            {   /*Z0311=17752*/
                argq = q*radius/(zr+1);   /*Z0311=17753*/
                pqr = (1/(zr*(zr-1))*pow(argq,-2));   /*Z0311=17754*/
                pq1 = pqr*cos((zr-1)*atan(argq))/pow(1+argq*argq,(zr-1)/2.0);   /*Z0311=17755*/
                pq2 = (pqr/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1+argq*argq,(zr-2)/2.0);   /*Z0311=17756*/
                pq3 = 3*(pq2-pq1);   /*Z0311=17757*/
                return pq3*pq3;   /*Z0311=17758*/
            }   /*Z0311=17759*/
        } /* of homogeneous sphere*/   /*Z0311=17760*/

        /*** core/shell sphere ***/   /*Z0311=17762*/
        if ( cs==1 )
        {   /*Z0311=17763*/
            ccc1 = sqr(1-rho)*pow(p1,6);   /*Z0311=17764*/
            ccc2 = 2*rho*(1-rho)*pow(p1,3);   /*Z0311=17765*/
            ccc3 = rho*rho;   /*Z0311=17766*/
            vv3 = sqr((1-rho)*pow(p1,3)+rho);   /*Z0311=17767*/

            argq = q*radiusm/(zr+1);   /*Z0311=17769*/
            argpq = q*radius/(zr+1);   /*Z0311=17770*/

            /* F121 sphere */   /*Z0311=17772*/
            if ( q<(limq4/2) )
            {   /*Z0311=17773*/
                qqn[0] = 1.0;   /*Z0311=17774*/
                pqsum = 1.0;   /*Z0311=17775*/
                oldpqsum = 0.0;   /*Z0311=17776*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=17777*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17778*/
                    pqsum = pqsum+qqn[nser]*carr4p[nser];   /*Z0311=17779*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17780*/
                    if ( delser<0.0001 ) break;   /*Z0311=17781*/
                    oldpqsum = pqsum;   /*Z0311=17782*/
                }   /*Z0311=17783*/
                //51:   /*Z0311=17784*/
                F121 = ccc1*pqsum/vv3;   /*Z0311=17785*/
            }   /*Z0311=17786*/
            else
            {   /*Z0311=17787*/
                pqr = (1/(zr*(zr-1))*pow(argpq,-2));   /*Z0311=17788*/
                pq1 = pqr*cos((zr-1)*atan(argpq))/pow(1+argpq*argpq,(zr-1)/2.0);   /*Z0311=17789*/
                pq2 = (pqr/((zr-2)*argpq))*sin((zr-2)*atan(argpq))/pow(1+argpq*argpq,(zr-2)/2.0);   /*Z0311=17790*/
                pq3 = 3*(pq2-pq1);   /*Z0311=17791*/
                F121 = ccc1*pq3*pq3/vv3;   /*Z0311=17792*/
            }   /*Z0311=17793*/

            /* F122 sphere */   /*Z0311=17795*/
            if ( q<(0.3*limq5) )
            {   /*Z0311=17796*/
                qqn[0] = 1.0;   /*Z0311=17797*/
                pqsum = 1.0;   /*Z0311=17798*/
                oldpqsum = 0.0;   /*Z0311=17799*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=17800*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17801*/
                    pqsum = pqsum+qqn[nser]*carr5p[nser];   /*Z0311=17802*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17803*/
                    if ( delser<0.0001 ) break;   /*Z0311=17804*/
                    oldpqsum = pqsum;   /*Z0311=17805*/
                }   /*Z0311=17806*/
                //52:   /*Z0311=17807*/
                F122 = ccc2*pqsum/vv3;   /*Z0311=17808*/
            }   /*Z0311=17809*/
            else
            {   /*Z0311=17810*/
                pqr1 = (1/(zr*(zr-1))*pow(argpq,-2));   /*Z0311=17811*/
                pq1 = pqr1*cos((zr-1)*atan(argpq))/pow(1+argpq*argpq,(zr-1)/2.0);   /*Z0311=17812*/
                pq2 = (pqr1/((zr-2)*argpq))*sin((zr-2)*atan(argpq))/pow(1+argpq*argpq,(zr-2)/2.0);   /*Z0311=17813*/
                pq3 = 3*(pq2-pq1);   /*Z0311=17814*/

                pqr2 = (1/(zr*(zr-1))*pow(argq,-2));   /*Z0311=17816*/
                pq4 = pqr2*cos((zr-1)*atan(argq))/pow(1+argq*argq,(zr-1)/2.0);   /*Z0311=17817*/
                pq5 = (pqr2/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1+argq*argq,(zr-2)/2.0);   /*Z0311=17818*/
                pq6 = 3*(pq5-pq4);   /*Z0311=17819*/

                F122 = ccc2*pq3*pq6/vv3;   /*Z0311=17821*/
            }   /*Z0311=17822*/

            /* F123 sphere */   /*Z0311=17824*/
            if ( q<(limq6/2) )
            {   /*Z0311=17825*/
                qqn[0] = 1.0;   /*Z0311=17826*/
                pqsum = 1.0;   /*Z0311=17827*/
                oldpqsum = 0.0;   /*Z0311=17828*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=17829*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=17830*/
                    pqsum = pqsum+qqn[nser]*carr6p[nser];   /*Z0311=17831*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=17832*/
                    if ( delser<0.0001 ) break;   /*Z0311=17833*/
                    oldpqsum = pqsum;   /*Z0311=17834*/
                }   /*Z0311=17835*/
                //53:   /*Z0311=17836*/
                F123 = ccc3*pqsum/vv3;   /*Z0311=17837*/
            }   /*Z0311=17838*/
            else
            {   /*Z0311=17839*/
                pqr = (1/(zr*(zr-1))*pow(argq,-2));   /*Z0311=17840*/
                pq1 = pqr*cos((zr-1)*atan(argq))/pow(1+argq*argq,(zr-1)/2.0);   /*Z0311=17841*/
                pq2 = (pqr/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1+argq*argq,(zr-2)/2.0);   /*Z0311=17842*/
                pq3 = 3*(pq2-pq1);   /*Z0311=17843*/
                F123 = ccc1*pq3*pq3/vv3;   /*Z0311=17844*/
            }   /*Z0311=17845*/
            /* formfq:=F121+F122+F123; */   /*Z0311=17846*/
            return F123;   /*Z0311=17847*/
        }  /* of core/shell sphere */   /*Z0311=17848*/

        /*** inhomogeneous core/shell sphere ***/   /*Z0311=17850*/
        if ( cs==2 )
        {   /*Z0311=17851*/

            dim = 3;   /*Z0311=17853*/
            delc = 0.0001;   /*Z0311=17854*/
            zz = zr;   /*Z0311=17855*/
            xrad = q*radiusm;   /*Z0311=17856*/
            xradp = q*radius;   /*Z0311=17857*/
            x1z = q*radius/(2*(zz+1));   /*Z0311=17858*/
            x12z = x1z*x1z;   /*Z0311=17859*/
            x2z = q*radiusm/(2*(zz+1));   /*Z0311=17860*/
            x22z = x2z*x2z;   /*Z0311=17861*/

            lim = 18*exp(-5*params.sigma);   /*Z0311=17863*/
            lim1 = lim;   /*Z0311=17864*/
            //lim2 = lim*0.7;   /*Z0311=17865*/
            //lim3 = lim;   /*Z0311=17866*/
            lim4 = lim;   /*Z0311=17867*/
            //lim5 = lim*0.7;   /*Z0311=17868*/
            lim6 = lim*1.2;   /*Z0311=17869*/

            a1 = (dim-alfa)/2.0;   /*Z0311=17871*/
            b1 = dim/2.0;   /*Z0311=17872*/
            b2 = (dim+2-alfa)/2.0;   /*Z0311=17873*/
            b1s = (dim+2)/2.0;   /*Z0311=17874*/
            v = -b1s+1/2.0;   /*Z0311=17875*/
            c = a1-b1-b2+1/2.0;   /*Z0311=17876*/
            //d0 = 1;   /*Z0311=17877*/
            //d1 = a1*(1+a1-b1)*(1+a1-b2);   /*Z0311=17878*/
            e0 = 1.0;   /*Z0311=17879*/
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;   /*Z0311=17880*/
            ee0 = 1.0;   /*Z0311=17881*/
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16*(1-b1s));   /*Z0311=17882*/

            gb1s = 3*sqrt(M_PI)/4.0;   /*Z0311=17884*/
            pz2v = 1/(zr*(zr-1));   /*Z0311=17885*/
            pz2v1 = pz2v/(zr-2);   /*Z0311=17886*/
            //pz2v2 = pz2v1/(zr-3);   /*Z0311=17887*/

            gz1 = gamma(zr+1);   /*Z0311=17889*/
            preg1 = gb1s/sqrt(M_PI);   /*Z0311=17890*/
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));   /*Z0311=17891*/
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));   /*Z0311=17892*/
            //pzvc = gamma(zr+1+v+c)/gz1;   /*Z0311=17893*/
            pzvc1 = gamma(zr+1+v+c-1)/gz1;   /*Z0311=17894*/
            pzvc2 = gamma(zr+1+v+c-2)/gz1;   /*Z0311=17895*/
            pzac = gamma(zr+1-2*a1+c)/gz1;   /*Z0311=17896*/
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;   /*Z0311=17897*/
            pzac2 = gamma(zr+1-2*a1+c+2)/gz1;   /*Z0311=17898*/
            pzc = gamma(zr+1+c)/gz1;   /*Z0311=17899*/
            pzc1 = gamma(zr+1+c-1)/gz1;   /*Z0311=17900*/
            pza = gamma(zr+1-2*a1)/gz1;   /*Z0311=17901*/
            pzva = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=17902*/
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=17903*/
            dnv0 = 1;   /*Z0311=17904*/
            pvav0 = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=17905*/
            pvav10 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=17906*/
            pva0 = gamma(zr+1-2*a1)/gz1;   /*Z0311=17907*/

            cc1 = 1/dim;   /*Z0311=17909*/
            cc4 = rho/((dim-alfa)*pow(p1,dim-alfa));   /*Z0311=17910*/
            cc6 = -rho/(dim-alfa);   /*Z0311=17911*/
            sumc = cc1+cc4+cc6;   /*Z0311=17912*/

            /* term #1 series */   /*Z0311=17914*/
            if ( xradp<lim1 )
            {   /*Z0311=17915*/
                z12v[0] = 1;   /*Z0311=17916*/
                b1sv[0] = 1;   /*Z0311=17917*/
                fkv[0] = 1;   /*Z0311=17918*/
                qqn[0] = 1.0;   /*Z0311=17919*/
                F12sez = 1.0;   /*Z0311=17920*/
                oldF12sez = 1.0;   /*Z0311=17921*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=17922*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=17923 TODO-F: hier war qnn[] geschrieben*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=17924*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=17925*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=17926*/
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]); */   /*Z0311=17927*/

                    F12sez = F12sez+carr4p[n]*qqn[n];   /*Z0311=17929*/

                    del = fabs((F12sez-oldF12sez)/F12sez);   /*Z0311=17931*/
                    if ( del<delc ) break;   /*Z0311=17932*/
                    oldF12sez = F12sez;   /*Z0311=17933*/
                }   /*Z0311=17934*/
                //101:   /*Z0311=17935*/
                F12 = F12sez;   /*Z0311=17936*/
            }   /*Z0311=17937*/
            /*Z0311=17938*/
            /* term #4 series */   /*Z0311=17939*/
            if ( xradp<lim4 )
            {   /*Z0311=17940*/
                z12v[0] = 1;   /*Z0311=17941*/
                a1v[0] = 1;   /*Z0311=17942*/
                b1v[0] = 1;   /*Z0311=17943*/
                b2v[0] = 1;   /*Z0311=17944*/
                b1sv[0] = 1;   /*Z0311=17945*/
                fkv[0] = 1;   /*Z0311=17946*/
                qqn[0] = 1.0;   /*Z0311=17947*/
                F42sez = 1.0;   /*Z0311=17948*/
                oldF42sez = 1.0;   /*Z0311=17949*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=17950*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=17951*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=17952*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=17953*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=17954*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=17955*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=17956*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=17957*/
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]); */   /*Z0311=17958*/

                    F42sez = F42sez+carr5p[n]*qqn[n];   /*Z0311=17960*/

                    del = fabs((F42sez-oldF42sez)/F42sez);   /*Z0311=17962*/
                    if ( del<delc ) break;   /*Z0311=17963*/
                    oldF42sez = F42sez;   /*Z0311=17964*/
                }   /*Z0311=17965*/
                //104:   /*Z0311=17966*/
                F42 = F42sez;   /*Z0311=17967*/
            }   /*Z0311=17968*/

            /* term #6 series */   /*Z0311=17970*/
            if ( xradp<lim6 )
            {   /*Z0311=17971*/
                z12v[0] = 1;   /*Z0311=17972*/
                a1v[0] = 1;   /*Z0311=17973*/
                b1v[0] = 1;   /*Z0311=17974*/
                b2v[0] = 1;   /*Z0311=17975*/
                b1sv[0] = 1;   /*Z0311=17976*/
                fkv[0] = 1;   /*Z0311=17977*/
                qqn[0] = 1.0;   /*Z0311=17978*/
                F62sez = 1.0;   /*Z0311=17979*/
                oldF62sez = 1.0;   /*Z0311=17980*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=17981*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=17982*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=17983*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=17984*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=17985*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=17986*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=17987*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=17988*/
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]); */   /*Z0311=17989*/

                    F62sez = F62sez+carr6p[n]*qqn[n];   /*Z0311=17991*/

                    del = fabs((F62sez-oldF62sez)/F62sez);   /*Z0311=17993*/
                    if ( del<delc ) break;   /*Z0311=17994*/
                    oldF62sez = F62sez;   /*Z0311=17995*/
                }   /*Z0311=17996*/
                //106:   /*Z0311=17997*/
                F62 = F62sez;   /*Z0311=17998*/
            }   /*Z0311=17999*/

            /*** term #1 asymptote ***/   /*Z0311=18001*/
            if ( xradp>=lim1 )
            {   /*Z0311=18002*/
                arg11 = (zr+v+1)*atan(2*x1z);   /*Z0311=18003*/
                nen11 = pow(1+4*x1z*x1z,(zr+v+1)/2.0);   /*Z0311=18004*/
                arg12 = (zr+v)*atan(2*x1z);   /*Z0311=18005*/
                nen12 = pow(1+4*x1z*x1z,(zr+v)/2.0);   /*Z0311=18006*/

                F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);   /*Z0311=18008*/
                F12as2z = ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);   /*Z0311=18009*/
                F12asz = preg1*pow(x1z,v)*(F12as1z+F12as2z);   /*Z0311=18010*/
                F12 = F12asz;   /*Z0311=18011*/
            }   /*Z0311=18012*/

            /*** term #4 asymptote ***/   /*Z0311=18014*/
            if ( xrad>=lim4 )
            {   /*Z0311=18015*/
                F42as10z = preg4*pow(x22z,-a1);   /*Z0311=18016*/
                F42as1sumz = pva0;   /*Z0311=18017*/
                F42as1z = F42as10z*F42as1sumz;   /*Z0311=18018*/
                F42as1z0 = F42as10z*pza;   /***/   /*Z0311=18019*/

                F42as40z = preg3*pow(x2z,c);   /*Z0311=18021*/
                arg44 = (zr+c+1)*atan(2*x2z);   /*Z0311=18022*/
                nen44 = pow(1+4*x2z*x2z,(zr+c+1)/2.0);   /*Z0311=18023*/
                arg45 = (zr+c)*atan(2*x2z);   /*Z0311=18024*/
                nen45 = pow(1+4*x2z*x2z,(zr+c)/2.0);   /*Z0311=18025*/
                F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);   /*Z0311=18026*/
                F42as28 = e1*(1/(2*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=18027*/
                F42as4z = F42as40z*(F42as27+F42as28);   /*Z0311=18028*/
                F42asz = F42as1z+F42as4z;   /*Z0311=18029*/
                F42asz0 = F42as1z0+F42as4z;   /*Z0311=18030*/
                F42 = F42asz0;   /*Z0311=18031*/
            }   /*Z0311=18032*/

            /*** term #6 asymptote ***/   /*Z0311=18034*/
            if ( xradp>=lim6 )
            {   /*Z0311=18035*/
                F62as10z = preg4*pow(x12z,-a1);   /*Z0311=18036*/
                F62as1sumz = pva0;   /*Z0311=18037*/
                F62as1z = F62as10z*F62as1sumz;   /*Z0311=18038*/
                F62as1z0 = F62as10z*pza;     /***/   /*Z0311=18039*/

                F62as40z = preg3*pow(x1z,c);   /*Z0311=18041*/
                arg64 = (zr+c+1)*atan(2*x1z);   /*Z0311=18042*/
                nen64 = pow(1+4*x1z*x1z,(zr+c+1)/2.0);   /*Z0311=18043*/
                arg65 = (zr+c)*atan(2*x1z);   /*Z0311=18044*/
                nen65 = pow(1+4*x1z*x1z,(zr+c)/2.0);   /*Z0311=18045*/
                F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);   /*Z0311=18046*/
                F62as28 = e1*(1/(2*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=18047*/
                F62as4z = F62as40z*(F62as27+F62as28);   /*Z0311=18048*/
                F62asz = F62as1z+F62as4z;   /*Z0311=18049*/
                F62asz0 = F62as1z0+F62as4z;   /*Z0311=18050*/
                F62 = F62asz0;   /*Z0311=18051*/
            }   /*Z0311=18052*/

            /* FF1:=(cc1*F12+cc4*F42+cc6*F62)/sumc; */   /*Z0311=18054*/
            FF1 = (cc1*F12)/sumc;   /*Z0311=18055*/

            return FF1*FF1;   /*Z0311=18057*/
            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,3,sigmar,q); */   /*Z0311=18060*/

        }  /* of inhomogeneous core/shell sphere */   /*Z0311=18063*/

    } /* of sphere */   /*Z0311=18065*/


    /************/   /*Z0311=18068*/
    /* cylinder */   /*Z0311=18069*/
    /************/   /*Z0311=18070*/
    if ( part==1 )
    {   /*Z0311=18071*/
        /*** longitudinal part ***/   /*Z0311=18073*/
        /*** isotropic ***/   /*Z0311=18074*/
        if ( ordis==7 )
        {   /*Z0311=18075*/
            if ( q<(0.6*limq1) )
            {   /*Z0311=18076*/
                pqsum = 1.0;   /*Z0311=18077*/
                oldpqsum = 0.0;   /*Z0311=18078*/
                qqn[0] = 1.0;   /*Z0311=18079*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18080*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18081*/
                    pqsum = pqsum+carr1p[nser]*qqn[nser];   /*Z0311=18082*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18083*/
                    if ( delser<0.0001 ) break;   /*Z0311=18084*/
                    oldpqsum = pqsum;   /*Z0311=18085*/
                }   /*Z0311=18086*/
                //60:   /*Z0311=18087*/
                pql = pqsum;   /*Z0311=18088*/
            }   /*Z0311=18089*/
            else
            {   /* = P(q) */   /*Z0311=18090*/
                arglq = q*length/(zl+1);   /*Z0311=18091*/
                /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2)); */   /*Z0311=18092*/
                pql = (M_PI/(2*zl))*(1/arglq);   /*Z0311=18093*/
                pql = pql-(1/(2*zl*(zl-1)*arglq*arglq))*cos((zl-1)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-1)/2.0);   /*Z0311=18094*/
            }   /*Z0311=18095*/
        }   /* of isotropic */   /*Z0311=18096*/

        /* perfect */   /*Z0311=18098*/
        if ( ordis==6 )
        {   /*Z0311=18099*/
            if ( orcase==4 )
                pql = 1.0;   /*Z0311=18100*/
            else
            {   /*Z0311=18101*/
                if ( q<(0.6*limq1) )
                {   /*Z0311=18102*/
                    pqsum = 1.0;   /*Z0311=18103*/
                    oldpqsum = 0.0;   /*Z0311=18104*/
                    qqn[0] = 1.0;   /*Z0311=18105*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18106*/
                        qqn[nser] = qqn[nser-1]*(qxs+qys)*(qxs+qys);   /*Z0311=18107*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser];   /*Z0311=18108*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18109*/
                        if ( delser<0.0001 ) break;   /*Z0311=18110*/
                        oldpqsum = pqsum;   /*Z0311=18111*/
                    }   /*Z0311=18112*/
                    //65:   /*Z0311=18113*/
                    pql = pqsum;   /*Z0311=18114*/
                }   /*Z0311=18115*/
                else
                {   /* yes */   /*Z0311=18116*/
                    arglq = (qxs+qys+eps4)*length/(zl+1);   /*Z0311=18117*/
                    /* F(q) */   /*Z0311=18118*/
                    /* pql:=(1/zl)*(1/arglq)*sin(zl*arctan(arglq))/power(1+arglq*arglq,zl/2); */   /*Z0311=18119*/
                    /* pql:=pql*pql; */   /*Z0311=18120*/
                    /* P(q) */   /*Z0311=18121*/
                    pql = (1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-1)/2.0));   /*Z0311=18122*/
                }   /*Z0311=18123*/
            }   /*Z0311=18124*/
        }   /* of perfect */   /*Z0311=18125*/

        /* general */   /*Z0311=18127*/
        if ( ordis==0 )
        {   /*Z0311=18128*/
            if ( orcase==4 )
                pql = 1.0;   /*Z0311=18129*/
            else
            {   /*Z0311=18130*/
                if ( q<(0.2*limq1) )
                {   /*Z0311=18131*/
                    pqsum = 1.0;   /*Z0311=18132*/
                    oldpqsum = 0.0;   /*Z0311=18133*/
                    qxn[0] = 1.0;   /*Z0311=18134*/
                    qyn[0] = 1.0;   /*Z0311=18135*/
                    if ( orcase==1 )
                    {   /*Z0311=18136*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=18137*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=18138*/
                            qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=18139*/
                            binsum = 0.0;   /*Z0311=18140*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=18141*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=18142*/
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser]; */   /*Z0311=18143*/
                                binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qyn[mser]*qxn[nser-mser];   /*Z0311=18144*/
                            }   /*Z0311=18145*/
                            pqsum = pqsum+carr1p[nser]*binsum;   /*Z0311=18146*/
                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18147*/
                            if ( delser<0.0001 ) break; // goto 66;   /*Z0311=18148*/
                            oldpqsum = pqsum;   /*Z0311=18149*/
                        }   /*Z0311=18150*/
                    }   /*Z0311=18151*/

                    if ( orcase==2 )
                    {  /* x-axis */   /*Z0311=18153*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=18154*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=18155*/
                            qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=18156*/
                            binsum = 0.0;   /*Z0311=18157*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=18158*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=18159*/
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser]; */   /*Z0311=18160*/
                                binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];   /*Z0311=18161*/
                            }   /*Z0311=18162*/
                            pqsum = pqsum+carr1p[nser]*binsum;   /*Z0311=18163*/
                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18164*/
                            if ( delser<0.0001 ) break; // goto 66;   /*Z0311=18165*/
                            oldpqsum = pqsum;   /*Z0311=18166*/
                        }   /*Z0311=18167*/
                    }   /*Z0311=18168*/
                    //pql = pqsum;   /*Z0311=18169*/

                    if ( orcase==3 )
                    {  /* y-axis */   /*Z0311=18171*/
                        for ( nser=1; nser<=120; nser++ )
                        {   /*Z0311=18172*/
                            qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=18173*/
                            qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=18174*/
                            binsum = 0.0;   /*Z0311=18175*/
                            for ( mser=0; mser<=nser; mser++ )
                            {   /*Z0311=18176*/
                                /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=18177*/
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser]; */   /*Z0311=18178*/
                                binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[mser]*qyn[nser-mser];   /*Z0311=18179*/
                            }   /*Z0311=18180*/
                            pqsum = pqsum+carr1p[nser]*binsum;   /*Z0311=18181*/
                            delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18182*/
                            if ( delser<0.0001 ) break; // goto 66;   /*Z0311=18183*/
                            oldpqsum = pqsum;   /*Z0311=18184*/
                        }   /*Z0311=18185*/
                    }   /*Z0311=18186*/
                    //66:   /*Z0311=18187*/
                    pql = pqsum;   /*Z0311=18188*/
                }   /*Z0311=18189*/
                else
                {   /*Z0311=18190*/
                    /* F(q) */   /*Z0311=18191*/
                    /* qrombdeltac(length,sigmal,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,orcase,0,1,0,carr1p,pql); */   /*Z0311=18192*/
                    /* P(q) */   /*Z0311=18193*/
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,orcase,0,0,0,carr1p,pql);   /*Z0311=18194*/
                    pql = pql/norm;   /*Z0311=18195*/
                }   /*Z0311=18196*/
            }   /*Z0311=18197*/
        }   /* of general */   /*Z0311=18198*/

        /* transverse part */   /*Z0311=18200*/
        /* homogeneous cylinder */   /*Z0311=18201*/
        if ( cs==0 )
        {   /*Z0311=18202*/
            if ( q<(1.0*limq4) )
            {   /*Z0311=18203*/
                pqsum = 1.0;   /*Z0311=18204*/
                oldpqsum = 0.0;   /*Z0311=18205*/
                qqn[0] = 1.0;   /*Z0311=18206*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18207*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18208*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=18209*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18210*/
                    if ( delser<0.0001 ) break;   /*Z0311=18211*/
                    oldpqsum = pqsum;   /*Z0311=18212*/
                }   /*Z0311=18213*/
                //61:   /*Z0311=18214*/
                pqr = pqsum;   /*Z0311=18215*/
            }   /*Z0311=18216*/
            else
            {   /*Z0311=18217*/
                argpq = q*radius/(zr+1);   /*Z0311=18218*/
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1+argpq*argpq,(zr-1/2.0)/2.0);   /*Z0311=18219*/
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1+argpq*argpq,(zr-3/2.0)/2.0);   /*Z0311=18220*/
                pqr3 = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);   /*Z0311=18221*/
                pqr = pqr3*pqr3;   /*Z0311=18222*/
            }   /*Z0311=18223*/
            return pql*pqr;   /*Z0311=18224*/
            /* formfq:=pql; */   /*Z0311=18225*/
        } /* of homogeneous cylinder */   /*Z0311=18226*/

        /* homogeneous core/shell cylinder */   /*Z0311=18228*/
        if ( cs==1 )
        {   /*Z0311=18229*/
            ccc1 = sqr(1-rho)*pow(p1,4);   /*Z0311=18230*/
            ccc2 = 2*rho*(1-rho)*pow(p1,2);   /*Z0311=18231*/
            ccc3 = rho*rho;   /*Z0311=18232*/
            vv3 = sqr((1-rho)*pow(p1,2)+rho);   /*Z0311=18233*/

            argq = q*radiusm/(zz+1);   /*Z0311=18235*/
            argpq = q*radius/(zz+1);   /*Z0311=18236*/

            /* F121 cylinder */   /*Z0311=18238*/
            if ( q<(0.7*limq4) )
            {   /*Z0311=18239*/
                /*** series expansion ***/   /*Z0311=18240*/
                pqsum = 1.0;   /*Z0311=18241*/
                oldpqsum = 0.0;   /*Z0311=18242*/
                qqn[0] = 1.0;   /*Z0311=18243*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18244*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18245*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=18246*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18247*/
                    if ( delser<0.0001 ) break;   /*Z0311=18248*/
                    oldpqsum = pqsum;   /*Z0311=18249*/
                }   /*Z0311=18250*/
                //62:   /*Z0311=18251*/
                F121 = ccc1*pqsum/vv3;   /*Z0311=18252*/
            }   /*Z0311=18253*/
            else
            {   /*Z0311=18254*/
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1+argpq*argpq,(zr-1/2.0)/2.0);   /*Z0311=18255*/
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1+argpq*argpq,(zr-3/2.0)/2.0);   /*Z0311=18256*/
                pqr3 = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);   /*Z0311=18257*/
                pqr = pqr3;   /*Z0311=18258*/
                F121 = ccc1*pqr*pqr/vv3;   /*Z0311=18259*/
            }   /*Z0311=18260*/

            /* F122 cylinder */   /*Z0311=18262*/
            if ( q<(0.3*limq5) )
            {   /*Z0311=18263*/
                /*** series expansion ***/   /*Z0311=18264*/
                pqsum = 1.0;   /*Z0311=18265*/
                oldpqsum = 0.0;   /*Z0311=18266*/
                qqn[0] = 1.0;   /*Z0311=18267*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18268*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18269*/
                    pqsum = pqsum+carr5p[nser]*qqn[nser];   /*Z0311=18270*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18271*/
                    if ( delser<0.0001 ) break;   /*Z0311=18272*/
                    oldpqsum = pqsum;   /*Z0311=18273*/
                }   /*Z0311=18274*/
                //63:   /*Z0311=18275*/
                F122 = ccc2*pqsum/vv3;   /*Z0311=18276*/
            }   /*Z0311=18277*/
            else
            {   /*Z0311=18278*/
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1+argpq*argq,(zr-1/2.0)/2.0);   /*Z0311=18279*/
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1+argpq*argq,(zr-3/2.0)/2.0);   /*Z0311=18280*/
                pqr3 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argq,-3/2.0)*(sin((zr-1/2.0)*atan(argq))-cos((zr-1/2.0)*atan(argq)))/pow(1+argq*argq,(zr-1/2.0)/2.0);   /*Z0311=18281*/
                pqr4 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argq,-5/2.0)*(sin((zr-3/2.0)*atan(argq))+cos((zr-3/2.0)*atan(argq)))/pow(1+argq*argq,(zr-3/2.0)/2.0);   /*Z0311=18282*/
                pqr = (4/M_PI)*(pqr1+(9/16.0)*pqr2)*(pqr3+(9/16.0)*pqr4);   /*Z0311=18283*/
                F122 = ccc2*pqr/vv3;   /*Z0311=18284*/
            }   /*Z0311=18285*/

            /* F123 cylinder */   /*Z0311=18287*/
            if ( q<(0.6*limq6) )
            {   /*Z0311=18288*/
                /*** series expansion ***/   /*Z0311=18289*/
                pqsum = 1.0;   /*Z0311=18290*/
                oldpqsum = 0.0;   /*Z0311=18291*/
                qqn[0] = 1.0;   /*Z0311=18292*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18293*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18294*/
                    pqsum = pqsum+carr6p[nser]*qqn[nser];   /*Z0311=18295*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18296*/
                    if ( delser<0.0001 ) break;   /*Z0311=18297*/
                    oldpqsum = pqsum;   /*Z0311=18298*/
                }   /*Z0311=18299*/
                //64:   /*Z0311=18300*/
                F123 = ccc3*pqsum/vv3;   /*Z0311=18301*/
            }   /*Z0311=18302*/
            else
            {   /*Z0311=18303*/
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argq,-3/2.0)*(sin((zr-1/2.0)*atan(argq))-cos((zr-1/2.0)*atan(argq)))/pow(1+argq*argq,(zr-1/2.0)/2.0);   /*Z0311=18304*/
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argq,-5/2.0)*(sin((zr-3/2.0)*atan(argq))+cos((zr-3/2.0)*atan(argq)))/pow(1+argq*argq,(zr-3/2.0)/2.0);   /*Z0311=18305*/
                pqr3 = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);   /*Z0311=18306*/
                pqr = pqr3;   /*Z0311=18307*/
                F123 = ccc3*pqr*pqr/vv3;   /*Z0311=18308*/
            }   /*Z0311=18309*/
            return pql*(F121+F122+F123);   /*Z0311=18310*/
            /* formfq:=(F121+F122+F123); */   /*Z0311=18311*/
        } /* of homogeneous core/shell cylinder */   /*Z0311=18312*/

        /*** inhomogeneous core/shell cylinder ***/   /*Z0311=18314*/
        if ( cs==2 )
        {   /*Z0311=18315*/
            dim = 2;   /*Z0311=18317*/
            delc = 0.0001;   /*Z0311=18318*/
            zz = zr;   /*Z0311=18319*/
            xrad = q*radiusm;   /*Z0311=18320*/
            xradp = q*radius;   /*Z0311=18321*/
            x1z = q*radius/(2*(zz+1));   /*Z0311=18322*/
            x12z = x1z*x1z;   /*Z0311=18323*/
            x2z = q*radiusm/(2*(zz+1));   /*Z0311=18324*/
            x22z = x2z*x2z;   /*Z0311=18325*/

            lim = 18*exp(-5*params.sigma);   /*Z0311=18327*/
            lim1 = lim;   /*Z0311=18328*/
            //lim2 = lim*0.7;   /*Z0311=18329*/
            //lim3 = lim;   /*Z0311=18330*/
            lim4 = lim;   /*Z0311=18331*/
            //lim5 = lim*0.7;   /*Z0311=18332*/
            lim6 = lim*1.2;   /*Z0311=18333*/

            a1 = (dim-alfa)/2.0;   /*Z0311=18335*/
            b1 = dim/2.0;   /*Z0311=18336*/
            b2 = (dim+2-alfa)/2.0;   /*Z0311=18337*/
            b1s = (dim+2)/2.0;   /*Z0311=18338*/
            v = -b1s+1/2.0;   /*Z0311=18339*/
            c = a1-b1-b2+1/2.0;   /*Z0311=18340*/
            //d0 = 1;   /*Z0311=18341*/
            //d1 = a1*(1+a1-b1)*(1+a1-b2);   /*Z0311=18342*/
            e0 = 1.0;   /*Z0311=18343*/
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;   /*Z0311=18344*/
            ee0 = 1.0;   /*Z0311=18345*/
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16*(1-b1s));   /*Z0311=18346*/

            gb1s = 3*sqrt(M_PI)/4.0;   /*Z0311=18348*/
            pz2v = 1/(zr*(zr-1));   /*Z0311=18349*/
            pz2v1 = pz2v/(zr-2);   /*Z0311=18350*/
            //pz2v2 = pz2v1/(zr-3);   /*Z0311=18351*/

            gz1 = gamma(zr+1);   /*Z0311=18353*/
            preg1 = gb1s/sqrt(M_PI);   /*Z0311=18354*/
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));   /*Z0311=18355*/
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));   /*Z0311=18356*/
            //pzvc = gamma(zr+1+v+c)/gz1;   /*Z0311=18357*/
            pzvc1 = gamma(zr+1+v+c-1)/gz1;   /*Z0311=18358*/
            pzvc2 = gamma(zr+1+v+c-2)/gz1;   /*Z0311=18359*/
            pzac = gamma(zr+1-2*a1+c)/gz1;   /*Z0311=18360*/
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;   /*Z0311=18361*/
            pzac2 = gamma(zr+1-2*a1+c+2)/gz1;   /*Z0311=18362*/
            pzc = gamma(zr+1+c)/gz1;   /*Z0311=18363*/
            pzc1 = gamma(zr+1+c-1)/gz1;   /*Z0311=18364*/
            pza = gamma(zr+1-2*a1)/gz1;   /*Z0311=18365*/
            pzva = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=18366*/
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=18367*/
            dnv0 = 1;   /*Z0311=18368*/
            pvav0 = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=18369*/
            pvav10 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=18370*/
            pva0 = gamma(zr+1-2*a1)/gz1;   /*Z0311=18371*/

            cc1 = 1/dim;   /*Z0311=18373*/
            cc4 = rho/((dim-alfa)*pow(p1,dim-alfa));   /*Z0311=18374*/
            cc6 = -rho/(dim-alfa);   /*Z0311=18375*/
            sumc = cc1+cc4+cc6;   /*Z0311=18376*/

            /* term #1 series */   /*Z0311=18378*/
            if ( (xradp)<lim1 )
            {   /*Z0311=18379*/
                z12v[0] = 1;   /*Z0311=18380*/
                b1sv[0] = 1;   /*Z0311=18381*/
                fkv[0] = 1;   /*Z0311=18382*/
                qqn[0] = 1.0;   /*Z0311=18383*/
                F12sez = 1.0;   /*Z0311=18384*/
                oldF12sez = 1.0;   /*Z0311=18385*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=18386*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=18387*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=18388*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=18389*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=18390*/
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]); */   /*Z0311=18391*/

                    F12sez = F12sez+carr4p[n]*qqn[n];   /*Z0311=18393*/

                    del = fabs((F12sez-oldF12sez)/F12sez);   /*Z0311=18395*/
                    if ( del<delc ) break;   /*Z0311=18396*/
                    oldF12sez = F12sez;   /*Z0311=18397*/
                }   /*Z0311=18398*/
                //111:   /*Z0311=18399*/
                F12 = F12sez;   /*Z0311=18400*/
            }   /*Z0311=18401*/

            /* term #4 series */   /*Z0311=18403*/
            if ( xradp<lim4 )
            {   /*Z0311=18404*/
                z12v[0] = 1;   /*Z0311=18405*/
                a1v[0] = 1;   /*Z0311=18406*/
                b1v[0] = 1;   /*Z0311=18407*/
                b2v[0] = 1;   /*Z0311=18408*/
                b1sv[0] = 1;   /*Z0311=18409*/
                fkv[0] = 1;   /*Z0311=18410*/
                qqn[0] = 1.0;   /*Z0311=18411*/
                F42sez = 1.0;   /*Z0311=18412*/
                oldF42sez = 1.0;   /*Z0311=18413*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=18414*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=18415*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=18416*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=18417*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=18418*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=18419*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=18420*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=18421*/
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]); */   /*Z0311=18422*/

                    F42sez = F42sez+carr5p[n]*qqn[n];   /*Z0311=18424*/

                    del = fabs((F42sez-oldF42sez)/F42sez);   /*Z0311=18426*/
                    if ( del<delc ) break;   /*Z0311=18427*/
                    oldF42sez = F42sez;   /*Z0311=18428*/
                }   /*Z0311=18429*/
                //114:   /*Z0311=18430*/
                F42 = F42sez;   /*Z0311=18431*/
            }   /*Z0311=18432*/

            /* term #6 series */   /*Z0311=18434*/
            if ( xradp<lim6 )
            {   /*Z0311=18435*/
                z12v[0] = 1;   /*Z0311=18436*/
                a1v[0] = 1;   /*Z0311=18437*/
                b1v[0] = 1;   /*Z0311=18438*/
                b2v[0] = 1;   /*Z0311=18439*/
                b1sv[0] = 1;   /*Z0311=18440*/
                fkv[0] = 1;   /*Z0311=18441*/
                qqn[0] = 1.0;   /*Z0311=18442*/
                F62sez = 1.0;   /*Z0311=18443*/
                oldF62sez = 1.0;   /*Z0311=18444*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=18445*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=18446*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=18447*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=18448*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=18449*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=18450*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=18451*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=18452*/
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]); */   /*Z0311=18453*/

                    F62sez = F62sez+carr6p[n]*qqn[n];   /*Z0311=18455*/

                    del = fabs((F62sez-oldF62sez)/F62sez);   /*Z0311=18457*/
                    if ( del<delc ) break;   /*Z0311=18458*/
                    oldF62sez = F62sez;   /*Z0311=18459*/
                }   /*Z0311=18460*/
                //116:   /*Z0311=18461*/
                F62 = F62sez;   /*Z0311=18462*/
            }   /*Z0311=18463*/

            /*** term #1 asymptote ***/   /*Z0311=18465*/
            if ( xradp>=lim1 )
            {   /*Z0311=18466*/
                arg11 = (zr+v+1)*atan(2*x1z);   /*Z0311=18467*/
                nen11 = pow(1+4*x1z*x1z,(zr+v+1)/2.0);   /*Z0311=18468*/
                arg12 = (zr+v)*atan(2*x1z);   /*Z0311=18469*/
                nen12 = pow(1+4*x1z*x1z,(zr+v)/2.0);   /*Z0311=18470*/

                F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);   /*Z0311=18472*/
                F12as2z = ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);   /*Z0311=18473*/
                F12asz = preg1*pow(x1z,v)*(F12as1z+F12as2z);   /*Z0311=18474*/
                F12 = F12asz;   /*Z0311=18475*/
            }   /*Z0311=18476*/

            /*** term #4 asymptote ***/   /*Z0311=18478*/
            if ( xrad>=lim4 )
            {   /*Z0311=18479*/
                F42as10z = preg4*pow(x22z,-a1);   /*Z0311=18480*/
                F42as1sumz = pva0;   /*Z0311=18481*/
                F42as1z = F42as10z*F42as1sumz;   /*Z0311=18482*/
                F42as1z0 = F42as10z*pza;   /***/   /*Z0311=18483*/

                F42as40z = preg3*pow(x2z,c);   /*Z0311=18485*/
                arg44 = (zr+c+1)*atan(2*x2z);   /*Z0311=18486*/
                nen44 = pow(1+4*x2z*x2z,(zr+c+1)/2.0);   /*Z0311=18487*/
                arg45 = (zr+c)*atan(2*x2z);   /*Z0311=18488*/
                nen45 = pow(1+4*x2z*x2z,(zr+c)/2.0);   /*Z0311=18489*/
                F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);   /*Z0311=18490*/
                F42as28 = e1*(1/(2*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=18491*/
                F42as4z = F42as40z*(F42as27+F42as28);   /*Z0311=18492*/
                F42asz = F42as1z+F42as4z;   /*Z0311=18493*/
                F42asz0 = F42as1z0+F42as4z;   /*Z0311=18494*/
                F42 = F42asz0;   /*Z0311=18495*/
            }   /*Z0311=18496*/

            /*** term #6 asymptote ***/   /*Z0311=18498*/
            if ( xradp>=lim6 )
            {   /*Z0311=18499*/
                F62as10z = preg4*pow(x12z,-a1);   /*Z0311=18500*/
                F62as1sumz = pva0;   /*Z0311=18501*/
                F62as1z = F62as10z*F62as1sumz;   /*Z0311=18502*/
                F62as1z0 = F62as10z*pza;     /***/   /*Z0311=18503*/

                F62as40z = preg3*pow(x1z,c);   /*Z0311=18505*/
                arg64 = (zr+c+1)*atan(2*x1z);   /*Z0311=18506*/
                nen64 = pow(1+4*x1z*x1z,(zr+c+1)/2.0);   /*Z0311=18507*/
                arg65 = (zr+c)*atan(2*x1z);   /*Z0311=18508*/
                nen65 = pow(1+4*x1z*x1z,(zr+c)/2.0);   /*Z0311=18509*/
                F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);   /*Z0311=18510*/
                F62as28 = e1*(1/(2*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=18511*/
                F62as4z = F62as40z*(F62as27+F62as28);   /*Z0311=18512*/
                F62asz = F62as1z+F62as4z;   /*Z0311=18513*/
                F62asz0 = F62as1z0+F62as4z;   /*Z0311=18514*/
                F62 = F62asz0;   /*Z0311=18515*/
            }   /*Z0311=18516*/

            FF1 = (cc1*F12+cc4*F42+cc6*F62)/sumc;   /*Z0311=18518*/
            /* FF1:=(cc6*F62)/sumc; */   /*Z0311=18519*/

            return FF1*FF1;   /*Z0311=18522*/
            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,2,sigmar,q); */   /*Z0311=18524*/
        } /* of inhomogeneous core/shell cylinder*/   /*Z0311=18525*/

    } /* of cylinder */   /*Z0311=18528*/


    /********/   /*Z0311=18532*/
    /* disk */   /*Z0311=18533*/
    /********/   /*Z0311=18534*/
    if ( part==2 )
    {   /*Z0311=18535*/
        /*** longitudinal part ***/   /*Z0311=18537*/
        /*** isotropic ***/   /*Z0311=18538*/
        if ( ordis==7 )
        {   /*Z0311=18539*/
            if ( q<(0.5*limq1) )
            {   /*Z0311=18540*/
                pqsum = 1.0;   /*Z0311=18541*/
                oldpqsum = 0.0;   /*Z0311=18542*/
                qqn[0] = 1.0;   /*Z0311=18543*/
                for ( nser=1; nser<=80; nser++ )
                {   /*Z0311=18544*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18545*/
                    pqsum = pqsum+carr1p[nser]*qqn[nser];   /*Z0311=18546*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18547*/
                    if ( delser<0.0001 ) break;   /*Z0311=18548*/
                    oldpqsum = pqsum;   /*Z0311=18549*/
                }   /*Z0311=18550*/
                //70:   /*Z0311=18551*/
                pql = pqsum;   /*Z0311=18552*/
            }   /*Z0311=18553*/
            else
            {  /*  = P||(q)  */   /*Z0311=18554*/
                arglq = q*length/(zl+1);   /*Z0311=18555*/
                pql = (2/(zl*(zl-1)))*pow(arglq,-2);   /*Z0311=18556*/
            }   /*Z0311=18557*/
        }  /* of isotropic */   /*Z0311=18558*/

        /* perfect */   /*Z0311=18560*/
        if ( ordis==6 )
        {   /*Z0311=18561*/
            if ( q<(0.5*limq1) )
            {   /*Z0311=18562*/
                pqsum = 1.0;   /*Z0311=18563*/
                oldpqsum = 0.0;   /*Z0311=18564*/
                qqn[0] = 1.0;   /*Z0311=18565*/
                if ( orcase==1 )
                {   /*Z0311=18566*/
                    argq = qxs+qys;   /*Z0311=18567*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18568*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18569*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*pow(1-argq*argq,nser);   /*Z0311=18570*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18571*/
                        if ( delser<0.0001 ) break;   /*Z0311=18572*/
                        oldpqsum = pqsum;   /*Z0311=18573*/
                    }   /*Z0311=18574*/
                }   /*Z0311=18575*/
                if ( orcase==2 )
                {   /*Z0311=18576*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18577*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18578*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*pow(1-qxs*qxs,nser);   /*Z0311=18579*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18580*/
                        if ( delser<0.0001 ) break;   /*Z0311=18581*/
                        oldpqsum = pqsum;   /*Z0311=18582*/
                    }   /*Z0311=18583*/
                }   /*Z0311=18584*/
                if ( orcase==3 )
                {   /*Z0311=18585*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18586*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18587*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*pow(1-qys*qys,nser);   /*Z0311=18588*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18589*/
                        if ( delser<0.0001 ) break;   /*Z0311=18590*/
                        oldpqsum = pqsum;   /*Z0311=18591*/
                    }   /*Z0311=18592*/
                }   /*Z0311=18593*/
                if ( orcase==4 )
                {   /*Z0311=18594*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18595*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18596*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser];   /*Z0311=18597*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18598*/
                        if ( delser<0.0001 ) break;   /*Z0311=18599*/
                        oldpqsum = pqsum;   /*Z0311=18600*/
                    }   /*Z0311=18601*/
                }   /*Z0311=18602*/
                //76:   /*Z0311=18603*/
                pql = pqsum;   /*Z0311=18604*/
            }   /*Z0311=18605*/
            else
            {   /* ok */   /*Z0311=18606*/
                if ( (orcase==1) )
                {   /*Z0311=18607*/
                    double qnarg = qxs+qys;   /*Z0311=18608 TODO*/
                    arglq = sqrt(1.0-qnarg*qnarg)*q*length/(zl+1)+eps4;   /*Z0311=18609*/
                }   /*Z0311=18610*/
                if ( orcase==2 ) arglq = sqrt(1.0-qxs*qxs)*q*length/(zl+1)+eps4;   /*Z0311=18611*/
                if ( orcase==3 ) arglq = sqrt(1.0-qys*qys)*q*length/(zl+1)+eps4;   /*Z0311=18612*/
                if ( orcase==4 ) arglq = q*length/(zl+1)+eps4;   /*Z0311=18613*/

                /* F(q) */   /*Z0311=18615*/
                /* pqr1:=(gamma(zl-1/2)/gamma(zl+1))*power(arglq,-3/2)*(sin((zl-1/2)*arctan(arglq))-cos((zl-1/2)*arctan(arglq)))/power(1+arglq*arglq,(zl-1/2)/2); */   /*Z0311=18616*/
                /* pqr2:=(gamma(zl-3/2)/gamma(zl+1))*power(arglq,-5/2)*(sin((zl-3/2)*arctan(arglq))+cos((zl-3/2)*arctan(arglq)))/power(1+arglq*arglq,(zl-3/2)/2); */   /*Z0311=18617*/
                /* pqr3:=(2/sqrt(pi))*(pqr1+(9/16)*pqr2); */   /*Z0311=18618*/
                /* pql:=pqr3*pqr3; */   /*Z0311=18619*/

                /* P(q) */   /*Z0311=18621*/
                pqr1 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3);   /*Z0311=18622*/
                pqr2 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3)*sin((zl-2)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-2)/2.0);   /*Z0311=18623*/
                pqr3 = (1/(zl*(zl-1)*(zl-2)*(zl-3)))*pow(arglq,-4)*cos((zl-3)*atan(2*arglq))/pow(1+4*arglq*arglq,(zl-3)/2.0);   /*Z0311=18624*/
                pql = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);   /*Z0311=18625*/
            }   /*Z0311=18626*/
        }   /* of perfect */   /*Z0311=18627*/

        /* orientational distribution */   /*Z0311=18629*/
        if ( ordis==0 )
        {   /*Z0311=18630*/
            if ( orcase==1 )
            {   /*Z0311=18631*/
                if ( q<(1.2*limq1) )
                {   /*Z0311=18632*/
                    pqsum = 1.0;   /*Z0311=18633*/
                    oldpqsum = 0.0;   /*Z0311=18634*/
                    qqn[0] = 1.0;   /*Z0311=18635*/
                    qxn[0] = 1.0;   /*Z0311=18636*/
                    qyn[0] = 1.0;   /*Z0311=18637*/
                    /*Z0311=18638*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18639*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18640*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=18641*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=18642*/

                        binsum = 0.0;   /*Z0311=18644*/
                        for ( mser=0; mser<=nser; mser++ )
                        {   /*Z0311=18645*/
                            binsum1 = 0.0;   /*Z0311=18646*/
                            for ( lser=0; lser<=mser; lser++ )
                            {   /*Z0311=18647*/
                                /* indx:=lser+1+round(mser*(mser+1)/2); */   /*Z0311=18648*/
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser]; */   /*Z0311=18649*/
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];   /*Z0311=18650*/
                            }   /*Z0311=18651*/
                            /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=18652*/
                            /* binsum:=binsum+carr1pm[indx]*binsum1; */   /*Z0311=18653*/
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;   /*Z0311=18654*/
                        }   /*Z0311=18655*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*binsum;   /*Z0311=18656*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18657*/
                        if ( delser<0.0001 ) break;   /*Z0311=18658*/
                        oldpqsum = pqsum;   /*Z0311=18659*/
                    }   /*Z0311=18660*/
                    //77:   /*Z0311=18661*/
                    pql = pqsum;   /*Z0311=18662*/
                }   /*Z0311=18663*/
                else
                {   /*Z0311=18664*/
                    /* disk: length = disk radius */   /*Z0311=18665*/
                    /* always use Bessel function approximation */   /*Z0311=18666*/
                    /* F(q) */   /*Z0311=18667*/
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql); */   /*Z0311=18668*/
                    /* P(q) */   /*Z0311=18669*/
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);   /*Z0311=18670*/
                    pql = pql/norm;   /*Z0311=18671*/
                }   /*Z0311=18672*/
            }   /*Z0311=18673*/

            if ( orcase==2 )
            {   /*Z0311=18675*/
                if ( q<(0.9*limq1) )
                {   /*Z0311=18676*/
                    pqsum = 1.0;   /*Z0311=18677*/
                    oldpqsum = 0.0;   /*Z0311=18678*/
                    qqn[0] = 1.0;   /*Z0311=18679*/
                    qxn[0] = 1.0;   /*Z0311=18680*/
                    qyn[0] = 1.0;   /*Z0311=18681*/

                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18683*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18684*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=18685*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=18686*/

                        binsum = 0.0;   /*Z0311=18688*/
                        for ( mser=0; mser<=nser; mser++ )
                        {   /*Z0311=18689*/
                            binsum1 = 0.0;   /*Z0311=18690*/
                            for ( lser=0; lser<=mser; lser++ )
                            {   /*Z0311=18691*/
                                /* indx:=lser+1+round(mser*(mser+1)/2); */   /*Z0311=18692*/
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser]; */   /*Z0311=18693*/
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];   /*Z0311=18694*/
                            }   /*Z0311=18695*/
                            /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=18696*/
                            /* binsum:=binsum+carr1pm[indx]*binsum1; */   /*Z0311=18697*/
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;   /*Z0311=18698*/
                        }   /*Z0311=18699*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*binsum;   /*Z0311=18700*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18701*/
                        if ( delser<0.0001 ) break;   /*Z0311=18702*/
                        oldpqsum = pqsum;   /*Z0311=18703*/
                    }   /*Z0311=18704*/
                    //78:   /*Z0311=18705*/
                    pql = pqsum;   /*Z0311=18706*/
                }   /*Z0311=18707*/
                else
                {   /*Z0311=18708*/
                    /* disk: length = disk radius */   /*Z0311=18709*/
                    /* always use Bessel function approximation */   /*Z0311=18710*/
                    /* F(q) */   /*Z0311=18711*/
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql); */   /*Z0311=18712*/
                    /* P(q) */   /*Z0311=18713*/
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);   /*Z0311=18714*/
                    pql = pql/norm;   /*Z0311=18715*/
                }   /*Z0311=18716*/
            }   /*Z0311=18717*/

            if ( orcase==3 )
            {   /*Z0311=18719*/
                if ( q<(0.9*limq1) )
                {   /*Z0311=18720*/
                    pqsum = 1.0;   /*Z0311=18721*/
                    oldpqsum = 0.0;   /*Z0311=18722*/
                    qqn[0] = 1.0;   /*Z0311=18723*/
                    qxn[0] = 1.0;   /*Z0311=18724*/
                    qyn[0] = 1.0;   /*Z0311=18725*/

                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18727*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18728*/
                        qxn[nser] = qxn[nser-1]*qxs*qxs;   /*Z0311=18729*/
                        qyn[nser] = qyn[nser-1]*qys*qys;   /*Z0311=18730*/

                        binsum = 0.0;   /*Z0311=18732*/
                        for ( mser=0; mser<=nser; mser++ )
                        {   /*Z0311=18733*/
                            binsum1 = 0.0;   /*Z0311=18734*/
                            for ( lser=0; lser<=mser; lser++ )
                            {   /*Z0311=18735*/
                                /* indx:=lser+1+round(mser*(mser+1)/2); */   /*Z0311=18736*/
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser]; */   /*Z0311=18737*/
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];   /*Z0311=18738*/
                            }   /*Z0311=18739*/
                            /* indx:=mser+1+round(nser*(nser+1)/2); */   /*Z0311=18740*/
                            /* binsum:=binsum+carr1pm[indx]*binsum1; */   /*Z0311=18741*/
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;   /*Z0311=18742*/
                        }   /*Z0311=18743*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser]*binsum;   /*Z0311=18744*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18745*/
                        if ( delser<0.0001 ) break;   /*Z0311=18746*/
                        oldpqsum = pqsum;   /*Z0311=18747*/
                    }   /*Z0311=18748*/
                    //79:   /*Z0311=18749*/
                    pql = pqsum;   /*Z0311=18750*/
                }   /*Z0311=18751*/
                else
                {   /*Z0311=18752*/
                    /* disk: length = disk radius */   /*Z0311=18753*/
                    /* always use Bessel function approximation */   /*Z0311=18754*/
                    /* F(q) */   /*Z0311=18755*/
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql); */   /*Z0311=18756*/
                    /* P(q) */   /*Z0311=18757*/
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);   /*Z0311=18758*/
                    pql = pql/norm;   /*Z0311=18759*/
                }   /*Z0311=18760*/
            }   /*Z0311=18761*/

            if ( orcase==4 )
            {   /*Z0311=18763*/
                if ( q<(0.5*limq1) )
                {   /*Z0311=18764*/
                    pqsum = 1.0;   /*Z0311=18765*/
                    oldpqsum = 0.0;   /*Z0311=18766*/
                    qqn[0] = 1.0;   /*Z0311=18767*/
                    for ( nser=1; nser<=120; nser++ )
                    {   /*Z0311=18768*/
                        qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18769*/
                        pqsum = pqsum+carr1p[nser]*qqn[nser];   /*Z0311=18770*/
                        delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18771*/
                        if ( delser<0.0001 ) break;   /*Z0311=18772*/
                        oldpqsum = pqsum;   /*Z0311=18773*/
                    }   /*Z0311=18774*/
                    //80:   /*Z0311=18775*/
                    pql = pqsum;   /*Z0311=18776*/
                }   /*Z0311=18777*/
                else
                {   /*Z0311=18778*/
                    /* disk: length = disk radius */   /*Z0311=18779*/
                    /* always use Bessel function approximation */   /*Z0311=18780*/
                    /* F(q) */   /*Z0311=18781*/
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql); */   /*Z0311=18782*/
                    /* P(q) */   /*Z0311=18783*/
                    qrombdeltac(length,radius,/*p1,sigmal,dbeta,*/theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,0,0,carr2p,pql);   /*Z0311=18784*/
                    pql = pql/norm;   /*Z0311=18785*/
                }   /*Z0311=18786*/
            }   /*Z0311=18787*/
        }   /* of orientational distribution */   /*Z0311=18788*/

        /* transverse part */   /*Z0311=18790*/
        /* disk: radius = disk thickness/2 */   /*Z0311=18791*/
        /* homogeneous disk */   /*Z0311=18792*/
        if ( cs==0 )
        {   /*Z0311=18793*/
            if ( q<(0.5*limq4) )
            {   /*Z0311=18794*/
                pqsum = 1.0;   /*Z0311=18795*/
                oldpqsum = 0.0;   /*Z0311=18796*/
                qqn[0] = 1.0;   /*Z0311=18797*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18798*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18799*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=18800*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18801*/
                    if ( delser<0.0001 ) break;   /*Z0311=18802*/
                    oldpqsum = pqsum;   /*Z0311=18803*/
                }   /*Z0311=18804*/
                //71:   /*Z0311=18805*/
                pqr = pqsum;   /*Z0311=18806*/
            }   /*Z0311=18807*/
            else
            {   /*Z0311=18808*/
                argpq = q*radius/(zr+1);   /*Z0311=18809*/
                pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1+argpq*argpq,zr/2.0);   /*Z0311=18810*/
            }   /*Z0311=18811*/
            return pql*pqr*pqr;   /*Z0311=18812*/
            /* formfq:=pql;; */   /*Z0311=18813*/
        } /* of homogeneous */   /*Z0311=18814*/

        /* core/shell disk */   /*Z0311=18816*/
        if ( cs==1 )
        {   /*Z0311=18817*/
            ccc1 = sqr(1-rho)*pow(p1,2);   /*Z0311=18818*/
            ccc2 = 2*rho*(1-rho)*pow(p1,1);   /*Z0311=18819*/
            ccc3 = rho*rho;   /*Z0311=18820*/
            vv3 = sqr((1-rho)*pow(p1,1)+rho);   /*Z0311=18821*/

            argq = q*radiusm/(zz+1);   /*Z0311=18823*/
            argpq = q*radius/(zz+1);   /*Z0311=18824*/

            /* F121 disk */   /*Z0311=18826*/
            if ( q<(0.8*limq4) )
            {   /*Z0311=18827*/
                /*** series expansion ***/   /*Z0311=18828*/
                pqsum = 1.0;   /*Z0311=18829*/
                oldpqsum = 0.0;   /*Z0311=18830*/
                qqn[0] = 1.0;   /*Z0311=18831*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18832*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18833*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=18834*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18835*/
                    if ( delser<0.0001 ) break;   /*Z0311=18836*/
                    oldpqsum = pqsum;   /*Z0311=18837*/
                }   /*Z0311=18838*/
                //72:   /*Z0311=18839*/
                F121 = ccc1*pqsum/vv3;   /*Z0311=18840*/
            }   /*Z0311=18841*/
            else
            {   /*Z0311=18842*/
                pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1+argpq*argpq,zr/2.0);   /*Z0311=18843*/
                F121 = ccc1*pqr*pqr/vv3;   /*Z0311=18844*/
            }   /*Z0311=18845*/

            /* F122 disk */   /*Z0311=18847*/
            if ( q<(1.0*limq5) )
            {   /*Z0311=18848*/
                /*** series expansion ***/   /*Z0311=18849*/
                pqsum = 1.0;   /*Z0311=18850*/
                oldpqsum = 0.0;   /*Z0311=18851*/
                qqn[0] = 1.0;   /*Z0311=18852*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18853*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18854*/
                    pqsum = pqsum+carr5p[nser]*qqn[nser];   /*Z0311=18855*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18856*/
                    if ( delser<0.0001 ) break;   /*Z0311=18857*/
                    oldpqsum = pqsum;   /*Z0311=18858*/
                }   /*Z0311=18859*/
                //73:   /*Z0311=18860*/
                F122 = ccc2*pqsum/vv3;   /*Z0311=18861*/
            }   /*Z0311=18862*/
            else
            {   /*Z0311=18863*/
                pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1+argpq*argpq,zr/2.0);   /*Z0311=18864*/
                pqr1 = (1/zr)*pow(argq,-1)*sin(zr*atan(argq))/pow(1+argq*argq,zr/2.0);   /*Z0311=18865*/
                F122 = ccc2*pqr*pqr1/vv3;   /*Z0311=18866*/
            }   /*Z0311=18867*/

            /* F123 disk */   /*Z0311=18869*/
            if ( q<(0.3*limq6) )
            {   /*Z0311=18870*/
                /*** series expansion ***/   /*Z0311=18871*/
                pqsum = 1.0;   /*Z0311=18872*/
                oldpqsum = 0.0;   /*Z0311=18873*/
                qqn[0] = 1.0;   /*Z0311=18874*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=18875*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=18876*/
                    pqsum = pqsum+carr6p[nser]*qqn[nser];   /*Z0311=18877*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=18878*/
                    if ( delser<0.0001 ) break;   /*Z0311=18879*/
                    oldpqsum = pqsum;   /*Z0311=18880*/
                }   /*Z0311=18881*/
                //74:   /*Z0311=18882*/
                F123 = ccc3*pqsum/vv3;   /*Z0311=18883*/
            }   /*Z0311=18884*/
            else
            {   /*Z0311=18885*/
                pqr = (1/zr)*pow(argq,-1)*sin(zr*atan(argq))/pow(1+argq*argq,zr/2.0);   /*Z0311=18886*/
                F123 = ccc3*pqr*pqr/vv3;   /*Z0311=18887*/
                /* add more terms, if necessary */   /*Z0311=18888*/
            }   /*Z0311=18889*/
            return pql*(F121+F122+F123);   /*Z0311=18890*/
            /* formfq:=(F121+F122+F123); */   /*Z0311=18891*/
        } /* of core/shell-disk */   /*Z0311=18892*/

        /*** inhomogeneous core/shell disk ***/   /*Z0311=18894*/
        if ( cs==2 )
        {   /*Z0311=18895*/
            dim = 1;   /*Z0311=18897*/
            delc = 0.0001;   /*Z0311=18898*/
            zz = zr;   /*Z0311=18899*/
            xrad = q*radiusm;   /*Z0311=18900*/
            xradp = q*radius;   /*Z0311=18901*/
            x1z = q*radius/(2*(zz+1));   /*Z0311=18902*/
            x12z = x1z*x1z;   /*Z0311=18903*/
            x2z = q*radiusm/(2*(zz+1));   /*Z0311=18904*/
            x22z = x2z*x2z;   /*Z0311=18905*/

            lim = 18*exp(-5*params.sigma);   /*Z0311=18907*/
            lim1 = lim;   /*Z0311=18908*/
            //lim2 = lim*0.7;   /*Z0311=18909*/
            //lim3 = lim;   /*Z0311=18910*/
            lim4 = lim;   /*Z0311=18911*/
            //lim5 = lim*0.7;   /*Z0311=18912*/
            lim6 = lim*1.2;   /*Z0311=18913*/

            a1 = (dim-alfa)/2.0;   /*Z0311=18915*/
            b1 = dim/2.0;   /*Z0311=18916*/
            b2 = (dim+2-alfa)/2.0;   /*Z0311=18917*/
            b1s = (dim+2)/2.0;   /*Z0311=18918*/
            v = -b1s+1/2.0;   /*Z0311=18919*/
            c = a1-b1-b2+1/2.0;   /*Z0311=18920*/
            //d0 = 1;   /*Z0311=18921*/
            //d1 = a1*(1+a1-b1)*(1+a1-b2);   /*Z0311=18922*/
            e0 = 1.0;   /*Z0311=18923*/
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;   /*Z0311=18924*/
            ee0 = 1.0;   /*Z0311=18925*/
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16*(1-b1s));   /*Z0311=18926*/

            gb1s = 3*sqrt(M_PI)/4.0;   /*Z0311=18928*/
            pz2v = 1/(zr*(zr-1));   /*Z0311=18929*/
            pz2v1 = pz2v/(zr-2);   /*Z0311=18930*/
            //pz2v2 = pz2v1/(zr-3);   /*Z0311=18931*/

            gz1 = gamma(zr+1);   /*Z0311=18933*/
            preg1 = gb1s/sqrt(M_PI);   /*Z0311=18934*/
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));   /*Z0311=18935*/
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));   /*Z0311=18936*/
            //pzvc = gamma(zr+1+v+c)/gz1;   /*Z0311=18937*/
            pzvc1 = gamma(zr+1+v+c-1)/gz1;   /*Z0311=18938*/
            pzvc2 = gamma(zr+1+v+c-2)/gz1;   /*Z0311=18939*/
            pzac = gamma(zr+1-2*a1+c)/gz1;   /*Z0311=18940*/
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;   /*Z0311=18941*/
            pzac2 = gamma(zr+1-2*a1+c+2)/gz1;   /*Z0311=18942*/
            pzc = gamma(zr+1+c)/gz1;   /*Z0311=18943*/
            pzc1 = gamma(zr+1+c-1)/gz1;   /*Z0311=18944*/
            pza = gamma(zr+1-2*a1)/gz1;   /*Z0311=18945*/
            pzva = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=18946*/
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=18947*/
            dnv0 = 1;   /*Z0311=18948*/
            pvav0 = gamma(zr+1+v-2*a1)/gz1;   /*Z0311=18949*/
            pvav10 = gamma(zr+1+v-2*a1-1)/gz1;   /*Z0311=18950*/
            pva0 = gamma(zr+1-2*a1)/gz1;   /*Z0311=18951*/

            cc1 = 1/dim;   /*Z0311=18953*/
            cc4 = rho/((dim-alfa)*pow(p1,dim-alfa));   /*Z0311=18954*/
            cc6 = -rho/(dim-alfa);   /*Z0311=18955*/
            sumc = cc1+cc4+cc6;   /*Z0311=18956*/

            /* term #1 series */   /*Z0311=18958*/
            if ( (xradp)<lim1 )
            {   /*Z0311=18959*/
                z12v[0] = 1;   /*Z0311=18960*/
                b1sv[0] = 1;   /*Z0311=18961*/
                fkv[0] = 1;   /*Z0311=18962*/
                qqn[0] = 1.0;   /*Z0311=18963*/
                F12sez = 1.0;   /*Z0311=18964*/
                oldF12sez = 1.0;   /*Z0311=18965*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=18966*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=18967*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=18968*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=18969*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=18970*/
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]); */   /*Z0311=18971*/

                    F12sez = F12sez+carr4p[n]*qqn[n];   /*Z0311=18973*/

                    del = fabs((F12sez-oldF12sez)/F12sez);   /*Z0311=18975*/
                    if ( del<delc ) break;   /*Z0311=18976*/
                    oldF12sez = F12sez;   /*Z0311=18977*/
                }   /*Z0311=18978*/
                //121:   /*Z0311=18979*/
                F12 = F12sez;   /*Z0311=18980*/
            }   /*Z0311=18981*/

            /* term #4 series */   /*Z0311=18983*/
            if ( xradp<lim4 )
            {   /*Z0311=18984*/
                z12v[0] = 1;   /*Z0311=18985*/
                a1v[0] = 1;   /*Z0311=18986*/
                b1v[0] = 1;   /*Z0311=18987*/
                b2v[0] = 1;   /*Z0311=18988*/
                b1sv[0] = 1;   /*Z0311=18989*/
                fkv[0] = 1;   /*Z0311=18990*/
                qqn[0] = 1.0;   /*Z0311=18991*/
                F42sez = 1.0;   /*Z0311=18992*/
                oldF42sez = 1.0;   /*Z0311=18993*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=18994*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=18995*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=18996*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=18997*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=18998*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=18999*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=19000*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=19001*/
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]); */   /*Z0311=19002*/

                    F42sez = F42sez+carr5p[n]*qqn[n];   /*Z0311=19004*/

                    del = fabs((F42sez-oldF42sez)/F42sez);   /*Z0311=19006*/
                    if ( del<delc ) break;   /*Z0311=19007*/
                    oldF42sez = F42sez;   /*Z0311=19008*/
                }   /*Z0311=19009*/
                //124:   /*Z0311=19010*/
                F42 = F42sez;   /*Z0311=19011*/
            }   /*Z0311=19012*/

            /* term #6 series */   /*Z0311=19014*/
            if ( xradp<lim6 )
            {   /*Z0311=19015*/
                z12v[0] = 1;   /*Z0311=19016*/
                a1v[0] = 1;   /*Z0311=19017*/
                b1v[0] = 1;   /*Z0311=19018*/
                b2v[0] = 1;   /*Z0311=19019*/
                b1sv[0] = 1;   /*Z0311=19020*/
                fkv[0] = 1;   /*Z0311=19021*/
                qqn[0] = 1.0;   /*Z0311=19022*/
                F62sez = 1.0;   /*Z0311=19023*/
                oldF62sez = 1.0;   /*Z0311=19024*/
                for ( n=1; n<=120; n++ )
                {   /*Z0311=19025*/
                    qqn[n] = qqn[n-1]*q*q;   /*Z0311=19026*/
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);   /*Z0311=19027*/
                    a1v[n] = a1v[n-1]*(a1-1+n);   /*Z0311=19028*/
                    b1v[n] = b1v[n-1]*(b1-1+n);   /*Z0311=19029*/
                    b2v[n] = b2v[n-1]*(b2-1+n);   /*Z0311=19030*/
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);   /*Z0311=19031*/
                    fkv[n] = fkv[n-1]*n;   /*Z0311=19032*/
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]); */   /*Z0311=19033*/

                    F62sez = F62sez+carr6p[n]*qqn[n];   /*Z0311=19035*/

                    del = fabs((F62sez-oldF62sez)/F62sez);   /*Z0311=19037*/
                    if ( del<delc ) break;   /*Z0311=19038*/
                    oldF62sez = F62sez;   /*Z0311=19039*/
                }   /*Z0311=19040*/
                //126:   /*Z0311=19041*/
                F62 = F62sez;   /*Z0311=19042*/
            }   /*Z0311=19043*/

            /*** term #1 asymptote ***/   /*Z0311=19045*/
            if ( xradp>=lim1 )
            {   /*Z0311=19046*/
                arg11 = (zr+v+1)*atan(2*x1z);   /*Z0311=19047*/
                nen11 = pow(1+4*x1z*x1z,(zr+v+1)/2.0);   /*Z0311=19048*/
                arg12 = (zr+v)*atan(2*x1z);   /*Z0311=19049*/
                nen12 = pow(1+4*x1z*x1z,(zr+v)/2.0);   /*Z0311=19050*/

                F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);   /*Z0311=19052*/
                F12as2z = ee1*(1/(2*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);   /*Z0311=19053*/
                F12asz = preg1*pow(x1z,v)*(F12as1z+F12as2z);   /*Z0311=19054*/
                F12 = F12asz;   /*Z0311=19055*/
            }   /*Z0311=19056*/

            /*** term #4 asymptote ***/   /*Z0311=19058*/
            if ( xrad>=lim4 )
            {   /*Z0311=19059*/
                F42as10z = preg4*pow(x22z,-a1);   /*Z0311=19060*/
                F42as1sumz = pva0;   /*Z0311=19061*/
                F42as1z = F42as10z*F42as1sumz;   /*Z0311=19062*/
                F42as1z0 = F42as10z*pza;   /***/   /*Z0311=19063*/

                F42as40z = preg3*pow(x2z,c);   /*Z0311=19065*/
                arg44 = (zr+c+1)*atan(2*x2z);   /*Z0311=19066*/
                nen44 = pow(1+4*x2z*x2z,(zr+c+1)/2.0);   /*Z0311=19067*/
                arg45 = (zr+c)*atan(2*x2z);   /*Z0311=19068*/
                nen45 = pow(1+4*x2z*x2z,(zr+c)/2.0);   /*Z0311=19069*/
                F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);   /*Z0311=19070*/
                F42as28 = e1*(1/(2*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);   /*Z0311=19071*/
                F42as4z = F42as40z*(F42as27+F42as28);   /*Z0311=19072*/
                F42asz = F42as1z+F42as4z;   /*Z0311=19073*/
                F42asz0 = F42as1z0+F42as4z;   /*Z0311=19074*/
                F42 = F42asz0;   /*Z0311=19075*/
            }   /*Z0311=19076*/

            /*** term #6 asymptote ***/   /*Z0311=19078*/
            if ( xradp>=lim6 )
            {   /*Z0311=19079*/
                F62as10z = preg4*pow(x12z,-a1);   /*Z0311=19080*/
                F62as1sumz = pva0;   /*Z0311=19081*/
                F62as1z = F62as10z*F62as1sumz;   /*Z0311=19082*/
                F62as1z0 = F62as10z*pza;     /***/   /*Z0311=19083*/

                F62as40z = preg3*pow(x1z,c);   /*Z0311=19085*/
                arg64 = (zr+c+1)*atan(2*x1z);   /*Z0311=19086*/
                nen64 = pow(1+4*x1z*x1z,(zr+c+1)/2.0);   /*Z0311=19087*/
                arg65 = (zr+c)*atan(2*x1z);   /*Z0311=19088*/
                nen65 = pow(1+4*x1z*x1z,(zr+c)/2.0);   /*Z0311=19089*/
                F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);   /*Z0311=19090*/
                F62as28 = e1*(1/(2*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);   /*Z0311=19091*/
                F62as4z = F62as40z*(F62as27+F62as28);   /*Z0311=19092*/
                F62asz = F62as1z+F62as4z;   /*Z0311=19093*/
                F62asz0 = F62as1z0+F62as4z;   /*Z0311=19094*/
                F62 = F62asz0;   /*Z0311=19095*/
            }   /*Z0311=19096*/

            FF1 = (cc1*F12+cc4*F42+cc6*F62)/sumc;   /*Z0311=19098*/
            /* FF1:=(cc1*F12)/sumc; */   /*Z0311=19099*/

            return FF1*FF1;   /*Z0311=19101*/
            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,1,sigmar,q); */   /*Z0311=19103*/

        } /* of inhomogeneous core/shell disk */   /*Z0311=19105*/
    } /* of disk */   /*Z0311=19106*/


    /* cube */   /*Z0311=19109*/
    if ( part==5 )
    {   /*Z0311=19110*/
        /* homogeneous cube */   /*Z0311=19111*/
        if ( cs==0 )
        {   /*Z0311=19112*/
            if ( q<0.7*limq4 )
            {   /*Z0311=19113*/
                pqsum = 1.0;   /*Z0311=19114*/
                oldpqsum = 0.0;   /*Z0311=19115*/
                qqn[0] = 1.0;   /*Z0311=19116*/
                for ( nser=1; nser<=120; nser++ )
                {   /*Z0311=19117*/
                    qqn[nser] = qqn[nser-1]*q*q;   /*Z0311=19118*/
                    pqsum = pqsum+carr4p[nser]*qqn[nser];   /*Z0311=19119*/
                    delser = fabs((pqsum-oldpqsum)/pqsum);   /*Z0311=19120*/
                    if ( delser<0.0001 ) break;   /*Z0311=19121*/
                    oldpqsum = pqsum;   /*Z0311=19122*/
                }   /*Z0311=19123*/
                //81:   /*Z0311=19124*/
                return pqsum;   /*Z0311=19125*/
            }   /*Z0311=19126*/
            else
            {   /*Z0311=19127*/
                qrombdeltac(length,radius,/*p1,sigma,dbeta,*/theta,phi,qx,qy,qz,qxhklt,qyhklt,qzhklt,qhkl,ax1.length(),ax2.length(),ax3.length(),ax1.x(),ax1.y(),ax1.z(),ax2.x(),ax2.y(),ax2.z(),ax3.x(),ax3.y(),ax3.z(),sig.x(),sig.y(),sig.z(),ordis,3,7,12,7,1,0,carr1p,pql);   /*Z0311=19128*/
                return pql/(M_PI/2.0);   /*Z0311=19129*/
            }   /*Z0311=19130*/
        } /* of homogeneous cube*/   /*Z0311=19131*/

        /* core/shell cube */   /*Z0311=19133*/
        if ( cs==1 )
        {   /*Z0311=19134*/
            return polycscube(1.0,rho,p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);   /*Z0311=19135*/
        }   /*Z0311=19136*/

    }  /* of cube */   /*Z0311=19138*/

    return 0.0;
}   /*Z0311=19139*/



#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::gammln( double xx ) const
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


#ifdef __CUDACC__
__host__ __device__
#endif
double CLASSLIB::polycscube(double /*rho1*/, double /*rho2*/, double /*p1*/, double /*p2*/, double /*alf1*/, double /*alf2*/, double /*rn*/, double /*pf*/,
                            double /*sigma*/, double /*q*/) const
{
    return 0;
#ifdef undef
    const int maxit = 10;
    const double eps = 0.00001;

    double intrc,oldintrc;
    int    itrc,trapzditrc;

    intrc = 0;
    /*** integral over alpha = 0 .. pi ***/
    /*** integrate over alpha = 0 .. pi/2 due to 4-fold symmetry ***/
    trapcscube(0+eps,M_PI/2.0-eps,rho1,rho2,p1,p2,alf1,alf2,rn,pf,sigma,q,intrc,1,trapzditrc);
    oldintrc = intrc;
    for ( itrc=2; itrc<=maxit; itrc++ )
    {
        trapcscube(0+eps,M_PI/2.0-eps,rho1,rho2,p1,p2,alf1,alf2,rn,pf,sigma,q,intrc,itrc,trapzditrc);
        if ( (oldintrc != 0.0) && (fabs(1-intrc/oldintrc) < eps) ) break;
        oldintrc = intrc;
    }
    /* normalization: division by  integral 0..pi/2  sin(theta) d(theta) = 1 */
    return intrc;
#endif
}


#ifdef undef
/* *********************** Core-shell cube- Integral ****************************** */
void trapcscube(double ac, double bc, double rho1, double rho2, double p1, double p2, double alf1, double alf2, double rn, double pf, double sigma, double q,
                double &sc,
                int nc,
                int &trapzditc)
{
    long jc;
    double fa,fb,fx,xc,tnmc,sumc,delc;

    if ( nc==1 )
    {
        fa = cscubealf(ac,rho1,rho2,p1,p2,alf1,alf2,rn,pf,sigma,q)*sin(ac);
        fb = cscubealf(bc,rho1,rho2,p1,p2,alf1,alf2,rn,pf,sigma,q)*sin(bc);
        sc = 0.5*(bc-ac)*(fa+fb);
        trapzditc = 1;
    }
    else
{
        tnmc = trapzditc;
        delc = (bc-ac)/tnmc;
        xc = ac+0.5*delc;
        sumc = 0.0;
        for ( jc=1; jc<=trapzditc; jc++ )
{
            fx = cscubealf(xc,rho1,rho2,p1,p2,alf1,alf2,rn,pf,sigma,q)*sin(xc);
            sumc = sumc+fx;
            xc = xc+delc;
        }
        sc = 0.5*(sc+(bc-ac)*sumc/tnmc);
        trapzditc = 2*trapzditc;
    }
}
#endif
