/*
 * Collection of common calculatoion routines.
 * Advantage:    avoid code duplicates and save program memory
 * Disadvantage: no code optimization due to partly used functions possible
 */

// THIS IS ONLY AN INCLUDE FILE FOR ANOTHER CLASS

// TODO: raus if ( i0 == 99 )
// TODO: raus {   // TEST x²*y³


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
void SasCalc_GENERIC_calculation::endThread()
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
                            std::cerr << "Canceled." << std::endl;
                            if ( doJoin )
                            {   // Nur beim ersten wird auf das Ende gewartet, die anderen sind dann auch beendet
                                pthread_join( threads[i], nullptr );
                                //std::cerr << "               finished." << std::endl;
                                doJoin = false;
                            }
                        }
                        else
                            std::cerr << "Cancelerror!" << s << std::endl;
                        if ( threads != nullptr ) threads[i] = 0;
                    }
                    else
                        std::cerr << "not running" << std::endl;
                }
            }
    }
}



//(************************** Gamma function ***************************)
/**
 * @brief SasCalc_GENERIC_calculation::gamma
 * @param z
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::gamma(double z) const          //{NV}
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



//(* *************************** vesicles **************************** *)
/**
 * @brief SasCalc_GENERIC_calculation::polyvesicle
 * @param ro     = SasCalc_GENERIC_calculation::length
 * @param ri     = SasCalc_GENERIC_calculation::radius
 * @param sigmar = SasCalc_GENERIC_calculation::sigma
 * @param sigmal = SasCalc_GENERIC_calculation::sigmal
 * @param q      = local var from calling routine
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::polyvesicle(double q) const        //{NV}
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



//(* ************************** 3d-integral over lorentz(x)*sin(x) ********************************* *)
/**
 * @brief SasCalc_GENERIC_calculation::lorentznorm3
 * @param a
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::lorentznorm3(double a) const           //{NV}
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



//(* ************************** 3d-integral over gauss(x)*sin(x) ********************************* *)
/**
 * @brief SasCalc_GENERIC_calculation::gaussnorm3
 * @param a
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::gaussnorm3(double a) const         //{NV}
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



//(* ************************** 3d-integral over pearson(x)*sin(x) ********************************* *)
/**
 * @brief SasCalc_GENERIC_calculation::pearsonnorm3
 * @param a
 * @param b
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::pearsonnorm3(double a, double b) const         //{NV}
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



//(* ********************* integration procedure for pearson(x)*sin(x) ***************************** *)
/**
 * @brief SasCalc_GENERIC_calculation::pearsonintegral3
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
void SasCalc_GENERIC_calculation::pearsonintegral3(double at, double bt, double a, double b,
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



#ifndef __CUDACC__
//#define UseStringGrid12
#endif

#define latpar1(a,b) latpar1ptr[latparIDX(a,b, 6)]      // [5000][6], genutzt: 0,1,2,3,4,5
#define latpar2(a,b) latpar2ptr[latparIDX(a,b, 6)]      // [5000][6], genutzt: 0,1,2,3,4,5
#define latpar3(a,b) latpar3ptr[latparIDX(a,b,14)]      // [5000][15], genutzt: 1 bis 12
//#define latpar4(a,b) latpar4ptr[latparIDX(a,b, 2)]      // [5000][15], genutzt: nichts

// Nur in prepareCalculation
void SasCalc_GENERIC_calculation::ButtonHKLClick( int ltype, int *latpar1ptr, int *latpar2ptr ) const
{   //Z=43355

    //const int np=20000;

    int /*i,j,*/ii/*,jj*/,h,k,l,hmax,kmax,lmax;//,index,xpos,ypos,image1width,image1height;
    int /*c0,*/c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,mult,/*multmax,hlev,klev,llev,zlev,zvar,*/ct1;
    // Unbenutzte Variablen: c3, c5, c9, c11, c12, c14 - sie werden aber hier gelassen, da ich sie nicht aus den Zuweisungs-case rausnehmen will
    int ccase,c1case,c2case/*,c3case*/,c4case/*,c5case*/,c6case,c7case,c8case/*,c9case*/,c10case;//,c11case,c12case;
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

    latpar1(  1,0)=0;
    latpar2(  1,0)=0;

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
#endif
        //EditCellb.Text:=FloatToStr(b);  EditCellc.Text:=FloatToStr(c);
    }
    if ( RadioButtonSysTetragonal )
    {
        b=a;  alf=M_PI/2.0;  gam=M_PI/2.0;
#ifdef UseStringGrid12
        bet=M_PI/2.0;  bmax=amax;
#endif
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




//(* *********************** Romberg integration ****************************** *)
//(* returns integral in the limits a and b *)

// LOKALE ROUTINE !!!
//(*** integration routine use trapezoidal rule ***)
#ifdef __CUDACC__
__host__ __device__
#endif
//void SasCalc_GENERIC_calculation::trapzddeltac( double a, double b, double l, double r, double dbeta, double theta, double phi,
//                             double qx, double qy, double qz, double p11, double p12, double p13, double p21,
//                             double p22, double p23, double p31, double p32, double p33,
//                             double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
//                             double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
//                             double ax3y, double ax3z, double sigx, double sigy, double sigz,
//                             int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
//                             double *carr1,
//                             double &pq, int n, int &trapzddeltac_cnt ) const
void SasCalc_GENERIC_calculation::trapzddeltac( double a, double b, double l, double r, double p1, double sigma, double alfa,
                             double dbeta, double theta, double phi, double qx, double qy, double qz,
                             double p11, double p12, double p13, double p21, double p22, double p23,
                             double p31, double p32, double p33, double qxn, double qyn, double qzn,
                             double qhkl, double ax1n, double ax2n, double ax3n, double ax1x, double ax1y,
                             double ax1z, double ax2x, double ax2y, double ax2z, double ax3x, double ax3y,
                             double ax3z, double sigx, double sigy, double sigz,
                             int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                             double *carr1, double &pq, int n, int &trapzddeltac_cnt ) const
{
    int j;  //Zupq1=2045
    double x, tnm, sump, /*sumf, sumn, sums,*/ del;  //Zupq1=2046
    double /*fa, fb, fx,*/ pa, pb, px; //, na, nb, nx, sa, sb, sx;  //Zupq1=2047

    if ( n==1 )
    {/*2*/  //Zupq1=2050
        switch ( i0 )
        {
        case 1:   /*  delta/chi integration  */  //Zupq1=2051
            switch ( i2 )
            {
            case 5: /*  delta and chi integration  */  //Zupq1=2052
                qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2053
                qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pb);  //Zupq1=2054
                pa = pa/(2.0*M_PI);  //Zupq1=2055
                pb = pb/(2.0*M_PI);  //Zupq1=2056
                break;  //Zupq1=2057
            default: /*  just delta integration  */  //Zupq1=2058
                pa = pow(cos(a),i3)*pow(sin(a),i4);  //Zupq1=2059
                pb = pow(cos(b),i3)*pow(sin(b),i4);  //Zupq1=2060
                break;  //Zupq1=2061
            case 0: /*  Gauss  */  //Zupq1=2062
                pa = sin(a)*exp(-a*a/(dbeta*dbeta))*pa;  //Zupq1=2063
                pb = sin(b)*exp(-b*b/(dbeta*dbeta))*pb;  //Zupq1=2064
                break;  //Zupq1=2065
            case 1: /*  Exponential  */  //Zupq1=2066
                pa = sin(a)*exp(-a/dbeta)*pa;  //Zupq1=2067
                pb = sin(b)*exp(-b/dbeta)*pb;  //Zupq1=2068
                break;  //Zupq1=2069
            case 2: /*  Onsager  */  //Zupq1=2070
                pa = sin(a)*exp(-sin(a)/dbeta)*pa;  //Zupq1=2071
                pb = sin(b)*exp(-sin(b)/dbeta)*pb;  //Zupq1=2072
                break;  //Zupq1=2073
            }  //Zupq1=2074
            break;
        case 2:   /*  norm  */  //Zupq1=2075
            switch ( i2 )
            {
            case 0: /*  Gauss  */  //Zupq1=2078
                pa = sin(a)*exp(-a*a/(dbeta*dbeta));  //Zupq1=2079
                pb = sin(b)*exp(-b*b/(dbeta*dbeta));  //Zupq1=2080
                break;  //Zupq1=2081
            case 1: /*  Exponential  */  //Zupq1=2082
                pa = sin(a)*exp(-a/dbeta*dbeta);  //Zupq1=2083
                pb = sin(b)*exp(-b/dbeta*dbeta);  //Zupq1=2084
                break;  //Zupq1=2085
            case 2: /*  Onsager  */  //Zupq1=2086
                pa = sin(a)*exp(-sin(a)/dbeta*dbeta);  //Zupq1=2087
                pb = sin(b)*exp(-sin(b)/dbeta*dbeta);  //Zupq1=2088
                break;  //Zupq1=2089
            }  //Zupq1=2092
            break;
        case 3:   /*  order parameter  */  //Zupq1=2093
            switch ( i2 )
            {
            case 0: /*  Gauss  */  //Zupq1=2096
                pa = sin(a)*exp(-a*a/(dbeta*dbeta))*(3*cos(a)*cos(a)-1)/2.0;  //Zupq1=2097
                pb = sin(b)*exp(-b*b/(dbeta*dbeta))*(3*cos(b)*cos(b)-1)/2.0;  //Zupq1=2098
                break;;  //Zupq1=2099
            case 1: /*  Exponential  */  //Zupq1=2100
                pa = sin(a)*exp(-a/dbeta)*(3*cos(a)*cos(a)-1)/2.0;  //Zupq1=2101
                pb = sin(b)*exp(-b/dbeta)*(3*cos(b)*cos(b)-1)/2.0;  //Zupq1=2102
                break;  //Zupq1=2103
            case 2: /*  Onsager  */  //Zupq1=2104
                pa = sin(a)*exp(-sin(a)/dbeta)*(3*cos(a)*cos(a)-1)/2.0;  //Zupq1=2105
                pb = sin(b)*exp(-sin(b)/dbeta)*(3*cos(b)*cos(b)-1)/2.0;  //Zupq1=2106
                break;  //Zupq1=2107
            }  //Zupq1=2110
            break;
        case 4:   /*  cylinder formfactor  */  //Zupq1=2111
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2112
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pb);  //Zupq1=2113
            switch ( i2 )
            {
            case 0: /*  Gauss  */  //Zupq1=2114
                pa = sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2.0*M_PI);  //Zupq1=2115
                pb = sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2.0*M_PI);  //Zupq1=2116
                break;  //Zupq1=2117
            case 1: /*  Exponential  */  //Zupq1=2118
                pa = sin(a)*exp(-a/dbeta)*pa/(2.0*M_PI);  //Zupq1=2119
                pb = sin(b)*exp(-b/dbeta)*pb/(2.0*M_PI);  //Zupq1=2120
                break;  //Zupq1=2121
            case 2: /*  Onsager  */  //Zupq1=2122
                pa = sin(a)*exp(-sin(a)/dbeta)*pa/(2.0*M_PI);  //Zupq1=2123
                pb = sin(b)*exp(-sin(b)/dbeta)*pb/(2.0*M_PI);  //Zupq1=2124
                break;  //Zupq1=2125
            }  //Zupq1=2126
            break;
        case 5:   /*  unit cell rotation  */  //Zupq1=2127
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2128
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pb);  //Zupq1=2129
            switch ( i2 )
            {
            case 0: /*  Gauss  */  //Zupq1=2130
                pa = 8*sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  //Zupq1=2131
                pb = 8*sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  //Zupq1=2132
                break;  //Zupq1=2133
            case 1: /*  Exponential  */  //Zupq1=2134
                pa = 8*sin(a)*exp(-a/dbeta)*pa/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  //Zupq1=2135
                pb = 8*sin(b)*exp(-b/dbeta)*pb/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  //Zupq1=2136
                break;  //Zupq1=2137
            case 2: /*  Onsager  */  //Zupq1=2138
                pa = 8*sin(a)*exp(-sin(a)/dbeta)*pa/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  //Zupq1=2139
                pb = 8*sin(b)*exp(-sin(b)/dbeta)*pb/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  //Zupq1=2140
                break;  //Zupq1=2141
            }  //Zupq1=2142
            break;
        case 6:   /*  disk formfactor  */  //Zupq1=2143
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2144
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pb);  //Zupq1=2145
            switch ( i2 )
            {
            case 0: /*  Gauss  */  //Zupq1=2146
                pa = sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2.0*M_PI);  //Zupq1=2147
                pb = sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2.0*M_PI);  //Zupq1=2148
                break;  //Zupq1=2149
            case 1: /*  Exponential  */  //Zupq1=2150
                pa = sin(a)*exp(-a/dbeta)*pa/(2.0*M_PI);  //Zupq1=2151
                pb = sin(b)*exp(-b/dbeta)*pb/(2.0*M_PI);  //Zupq1=2152
                break;  //Zupq1=2153
            case 2: /*  Onsager  */  //Zupq1=2154
                pa = sin(a)*exp(-sin(a)/dbeta)*pa/(2.0*M_PI);  //Zupq1=2155
                pb = sin(b)*exp(-sin(b)/dbeta)*pb/(2.0*M_PI);  //Zupq1=2156
                break;  //Zupq1=2157
                /* pa:=sin(a)*exp(-(a-pi/2)*(a-pi/2)/(dbeta*dbeta))*pa/(2*pi);  //Zupq1=2158 */
                /* pb:=sin(b)*exp(-(b-pi/2)*(b-pi/2)/(dbeta*dbeta))*pb/(2*pi);  //Zupq1=2159 */
            }  //Zupq1=2160
            break;
        case 7:    /*  cube-, triaxial ellipsoid-integration  */  //Zupq1=2161
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2162
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pb);  //Zupq1=2163
            pa = pa*sin(a);  //Zupq1=2164
            pb = pb*sin(b);  //Zupq1=2165
            break;  //Zupq1=2166
        case 8:   /*  superball integration  */  //Zupq1=2167
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2168
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pb);  //Zupq1=2169
            /* pa:=1.05;  //Zupq1=2170 */
            /* pb:=1.0;  //Zupq1=2171 */
            break;  //Zupq1=2172
        } // switch i0
        if ( i0 == 99 )
            pq = pa + pb;
        else
            pq = 0.5*(b-a)*(pa+pb);  //Zupq1=2173
        //DSM(
        //if ( fabs(pq) < 1.0e-170 )
        //    qDebug() << "trapzddeltac(1)" << pa << pb << pq << "i?" << i0 << i1 << i2 << i3 << i4;
        //)
        trapzddeltac_cnt = 1;  //Zupq1=2174
        DCNT( dbgCount++; )
    }
    else
    {/*2*/  //Zupq1=2176
        tnm = trapzddeltac_cnt;  //Zupq1=2177
        del = (b-a)/tnm;  //Zupq1=2178
        DCNT(
#ifdef DBGLIMIT
            if ( ++dbgCount > DBGLIMIT ) return;
#else
            dbgCount++;
#endif
            )
        x = a+0.5*del;  //Zupq1=2179
        sump = 0.0;  //Zupq1=2180
        for ( j=1; j<=trapzddeltac_cnt; j++ )
        {/*3*/  //Zupq1=2181
            switch ( i0 )
            {
            case 1:  /*  delta integration  */  //Zupq1=2182
                switch ( i2 )
                {
                case 5:  /*  delta and chi integration  */  //Zupq1=2184
                    qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,px);  //Zupq1=2185
                    px = px/(2.0*M_PI);  //Zupq1=2186
                    break;  //Zupq1=2187
                    /*  just delta integration  */  //Zupq1=2188
                default:
                    px = pow(cos(x),i3)*pow(sin(x),i4);  //Zupq1=2189
                    break;
                case 0:
                    px = sin(x)*exp(-x*x/(dbeta*dbeta))*px;  /*  Gauss  */  //Zupq1=2190
                    break;
                case 1:
                    px = sin(x)*exp(-x/dbeta)*px;  /*  Exponential  */  //Zupq1=2191
                    break;
                case 2:
                    px = sin(x)*exp(-sin(x)/dbeta)*px;  /*  Onsager  */  //Zupq1=2192
                    break;
                }/*4*/  //Zupq1=2193
                break;

            case 2:  /*  norm  */  //Zupq1=2195
                if ( i2==0 ) px = sin(x)*exp(-x*x/(dbeta*dbeta));  /*  Gauss  */  //Zupq1=2197
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta);  /*  Exponential  */  //Zupq1=2198
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta);  /*  Onsager  */  //Zupq1=2199
                break;  //Zupq1=2200

            case 3:  /*  order parameter  */  //Zupq1=2202
                if ( i2==0 ) px = sin(x)*exp(-x*x/(dbeta*dbeta))*(3*cos(x)*cos(x)-1)/2.0;  /*  Gauss  */  //Zupq1=2204
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta)*(3*cos(x)*cos(x)-1)/2.0;  /*  Exponential  */  //Zupq1=2205
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta)*(3*cos(x)*cos(x)-1)/2.0;  /*  Onsager  */  //Zupq1=2206
                break;  //Zupq1=2207
            /* if i0=3 then px:=sin(x)*exp(-(x-pi/2)*(x-pi/2)/(dbeta*dbeta))*(3*cos(x)*cos(x)-1)/2;  //Zupq1=2208 */

            case 4:  /*  cylinder formfactor  */  //Zupq1=2210
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,px);  //Zupq1=2211
                if ( i2==0 ) px = sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2.0*M_PI);   /*  Gauss  */  //Zupq1=2212
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta)*px/(2.0*M_PI);   /*  Exponential  */  //Zupq1=2213
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta)*px/(2.0*M_PI);   /*  Onsager  */  //Zupq1=2214
                break;  //Zupq1=2215

            case 5:  /*  unit cell rotation  */  //Zupq1=2216
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,px);  //Zupq1=2217
                if ( i2==0 ) px = 8*sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  /*  Gauss  */  //Zupq1=2218
                if ( i2==1 ) px = 8*sin(x)*exp(-x/dbeta)*px/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  /*  Exponential  */  //Zupq1=2219
                if ( i2==2 ) px = 8*sin(x)*exp(-sin(x)/dbeta)*px/(2.0*M_PI*pow(M_PI,3/2.0)*sigx*sigy*sigz);  /*  Onsager  */  //Zupq1=2220
                break;  //Zupq1=2221

            case 6:  /*  disk formfactor  */  //Zupq1=2222
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,px);  //Zupq1=2223
                if ( i2==0 ) px = sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2.0*M_PI);  /*  Gauss  */  //Zupq1=2224
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta)*px/(2.0*M_PI);  /*  Exponential  */  //Zupq1=2225
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta)*px/(2.0*M_PI);  /*  Onsager  */  //Zupq1=2226
                /* px:=sin(x)*exp(-(x-pi/2)*(x-pi/2)/(dbeta*dbeta))*px/(2*pi);  //Zupq1=2227 */
                break;  //Zupq1=2228

            case 7:  /*  cube-, triaxial ellipsoid-integration  */  //Zupq1=2229
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,px);  //Zupq1=2230
                px = px*sin(x);  //Zupq1=2231
                break;  //Zupq1=2232

            case 8:  /*  superball integration  */  //Zupq1=2233
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,px);  //Zupq1=2234
                /* px:=0.99;  //Zupq1=2235 */
                break;  //Zupq1=2236
            } // switch i0
            sump = sump+px;  //Zupq1=2237
            x = x+del;  //Zupq1=2238
        }/* for j */  //Zupq1=2239
        pq = 0.5*(pq+(b-a)*sump/tnm);  //Zupq1=2240
        trapzddeltac_cnt = 2*trapzddeltac_cnt;  //Zupq1=2241
    }/* n==1 else */  //Zupq1=2242
}  //Zupq1=2243



//(*** interpolation routine ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void SasCalc_GENERIC_calculation::polint( double *xa/*RealArrayNP*/, double *ya/*RealArrayNP*/,
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
// 20220730:
//procedure qrombdeltac(l,r,p1,sigma,alfa,dbeta,theta,phi,qx,qy,qz: extended;
//                      qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz: extended;
//                      ordis,dim,i0,i1,i2,i3,i4: integer;
//                      var carr1: Array of extended;
//                      var pq: extended);

#ifdef __CUDACC__
__host__ __device__
#endif
//void SasCalc_GENERIC_calculation::qrombdeltac(
//    double l,     = params.length   -> raus
//    double r,     = params.radius   -> raus
//    double p1,    = params.p1 | params.radiusi(nur einmal)   ==> bleibt
//    double sigma, = params.sigma | sigmal(formpq-Parameter) | params.sigmal   ==> bleibt
//    double alfa,  = params.alphash | params.alpha   ==> bleibt
//    double dbeta, = params.dbeta   -> raus
//    double theta, = params.polTheta | theta(loc.Var)   ==> bleibt
//    double phi,   = params.polPhi | phi(loc.Var)   ==> bleibt
//    double qx,    = qx(loc.Var) | qxs(Param) | 1   ==> bleibt
//    double qy,    = qy(loc.Var) | qys(Param) | 1   ==> bleibt
//    double qz,    = qz(loc.Var) | 1   ==> bleibt
//    double qxn,   = 9 | qxhkl(loc.Var)   ==> bleibt
//    double qyn,   = 9 | qyhkl(loc.Var)   ==> bleibt
//    double qzn,   = 9 | qzhkl(loc.Var)   ==> bleibt
//    double qhkl,  = 9 | qhkl(loc.Var)   ==> bleibt
//    double ax1n,  = 9 | params.ax1.length()
//    double ax2n,  = 9 | params.ax2.length()
//    double ax3n,  = 9 | params.ax3.length()
//    double ax1x,  = 9 | params.ax1.x()
//    double ax1y,  = 9 | params.ax1.y()
//    double ax1z,  = 9 | params.ax1.z()
//    double ax2x,  = 9 | params.ax2.x()
//    double ax2y,  = 9 | params.ax2.y()
//    double ax2z,  = 9 | params.ax2.z()
//    double ax3x,  = 9 | params.ax3.x()
//    double ax3y,  = 9 | params.ax3.y()
//    double ax3z,  = 9 | params.ax3.z()
//    double sigx,  = 9 | params.sig.x()
//    double sigy,  = 9 | params.sig.y()
//    double sigz,  = 9 | params.sig.z()
//    int ordis,    = ordis(Param in formpq) | params.ordis   ==> bleibt
//    int dim,      = 2 | 3 u.a.   ==> bleibt
//    int i0,       = 6 | 5 u.a.   ==> bleibt
//    int i1,       = params.orcase+7 | 6 u.a.   ==> bleibt
//    int i2,       = 0 u.a.   ==> bleibt
//    int i3,       = 0 u.a.   ==> bleibt
//    int i4,       = 0 u.a.   ==> bleibt
//    double *carr1,= params.CR->carr2p | params.CR->carr1p   ==> bleibt
//    double &pq    = (Returnwert)
//   ) const

void SasCalc_GENERIC_calculation::qrombdeltac( double p1, double sigma, double alfa,
                            double theta, double phi, double qx, double qy, double qz,
                            double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n,
                            double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y,
                            double ax2z, double ax3x, double ax3y, double ax3z, double sigx, double sigy,
                            double sigz, int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                            double *carr1, double &pq ) const
{
    /* label 99; */  //Zupq1=1968
    const double eps = 1.0e-4;  //Zupq1=1970
    //const double eps1 = 1.0e-8;  //Zupq1=1971
    const int    jmax = 20;  //Zupq1=1972
    //const int    jmaxp = 21;  //Zupq1=1973 global
    const int    k = 5;  //Zupq1=1974
    //const int    np = 5;  //Zupq1=1975 global
    const double prmin =  0.00001;       /*  minimum pr(delta)-value  */  //Zupq1=1976

    DCL( static long dbgCountMin=0;
         static long dbgCountMax=0;
         static double dbgCountSum = 0 );
    DTL( static time_t startzeit = 0;
         static long   dbgCntDiff = 0 );

    //typedef double RealArrayJMAXP[jmaxp+1];  //Zupq1=1978 global
    //typedef double RealArrayNP[np+1];  //Zupq1=1979 global

    int i, j;  //Zupq1=1981
    //int trapzdit, midpntit;  //Zupq1=1982
    double dssp, delmax; //, norm;  //Zupq1=1983
    RealArrayJMAXP hp, sp;  //Zupq1=1984
    RealArrayNP    cp, dp;  //Zupq1=1985
    float alim, blim, p11, p12, p13, p21, p22, p23, p31, p32, p33;       //Zupq1=1986

    int trapzddeltac_cnt;

    DTL( if ( startzeit == 0 ) startzeit = time(nullptr) );

    //D8L( qDebug() << "qrombdeltac(" << r << theta << phi << "qx/y/z" << qx << qy << qz << "q?n" << qxn << qyn << qzn << qhkl << "ax?n" << ax1n << ax2n << ax3n
    //                               << "ax1?" << ax1x << ax1y << ax1z << "ax2?" << ax2x << ax2y << ax2z << "ax3?" << ax3x << ax3y << ax3z
    //                               << "sig?" << sigx << sigy << sigz << "ordis" << ordis << dim << "i?" << i0 << i1 << i2 << i3 << i4 << "&pq )" << params.dbeta );

    //memset( hp, 0, sizeof(hp) );
    //memset( sp, 0, sizeof(sp) );
    //memset( cp, 0, sizeof(cp) );
    //memset( dp, 0, sizeof(dp) );

    // params.dbeta ist die globale Variable in Grad.
    double dbeta = params.dbeta*M_PI/180.0;  //Zupq1=2247
    theta = theta*M_PI/180.0;  //Zupq1=2248
    phi = phi*M_PI/180.0;  //Zupq1=2249

    DCL( qDebug() << "   qrombdeltac(" << params.radius << theta << phi << "qx/y/z" << qx << qy << qz << "..." // nur Kurzform
                                       << "i?" << i0 << i1 << i2 << i3 << i4 << " ) sigma:" << params.sigma );

    /*  search for maximum integration angle  */  //Zupq1=2250
    if ( (i2==0) || (i2==2) || (i2==3) || (i2==5) ) delmax = dbeta*sqrt(log(1.0/prmin));  /*  Gaussian, Onsager, Maier-Saupe, Laguerre  */  //Zupq1=2251
    if ( i2==1 ) delmax = dbeta*log(1.0/prmin);        /*  Exponential  */  //Zupq1=2252
    /* if i2=2 then delmax:=arcsin(dbeta*ln(1/prmin));  //Zupq1=2253 */
    if ( i2==4 ) delmax = dbeta;                    /*  cut-off  */  //Zupq1=2254
    if ( (i2==7) || (i2==8) || (i2==9) || (i2==10) || (i2==11) || (i2==12) ) delmax = M_PI/2.0; /*  isotropic, mirrored distributions  */  //Zupq1=2255

    if ( delmax>M_PI/2.0 ) delmax = M_PI/2.0;  //Zupq1=2257
    alim = 0.0;  //Zupq1=2258
    blim = delmax;  //Zupq1=2259

    DCL( qDebug() << "               dbeta" << dbeta << "delmax" << delmax << "a/blim" << alim << blim; )

    if ( i1==12 )
    {/*2*/   /*  for cubes  */  //Zupq1=2261
        alim = eps;  //Zupq1=2262
        blim = M_PI/2.0-eps;  //Zupq1=2263
    }/*2*/  //Zupq1=2264

    if ( i1==17 )
    {/*2*/   /*  for superballs  */  //Zupq1=2266
        alim = 0.0001;  //Zupq1=2267
        /* blim:=r-0.0001;  //Zupq1=2268 */
        blim = M_PI/2.0-0.0001;  //Zupq1=2269
    }/*2*/  //Zupq1=2270

    /*  for disks: mirrored Gaussian  */  //Zupq1=2273
    /* if ((qzn=2) or (i0=6)) then begin  //Zupq1=2274 */
    /*    alim:=pi/2-delmax;  //Zupq1=2275 */
    /*    blim:=pi/2;  //Zupq1=2276 */
    /* end;  //Zupq1=2277 */

    p11 = -cos(phi)*cos(theta);  //Zupq1=2279
    p12 = sin(phi);  //Zupq1=2280
    p13 = cos(phi)*sin(theta);  //Zupq1=2281
    p21 = -cos(phi);  //Zupq1=2282
    p22 = -sin(phi)*cos(theta);  //Zupq1=2283
    p23 = sin(phi)*sin(theta);  //Zupq1=2284
    p31 = -sin(theta);  //Zupq1=2285
    p32 = 0;  //Zupq1=2286
    p33 = -cos(theta);  //Zupq1=2287

    hp[1] = 1.0;  //Zupq1=2293
    for ( j=1; j<=jmax; j++ )
    {/*2*/  //Zupq1=2294
        trapzddeltac(alim,blim,params.length,params.radius,p1,sigma,alfa,dbeta,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,
                     qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,
                     ordis,dim,i0,i1,i2,i3,i4,carr1,sp[j],j, trapzddeltac_cnt );  //Zupq1=2295
        if ( j>=k )
        {/*3*/  //Zupq1=2296
            for ( i=1; i<=k; i++ )
            {/*4*/  //Zupq1=2297
                cp[i] = hp[j-k+i];  //Zupq1=2298
                dp[i] = sp[j-k+i];  //Zupq1=2299
            }/*4*/  //Zupq1=2300
            polint(cp,dp,k,0.0,pq,dssp,"qrombdeltac");  //Zupq1=2301
            D8L( qDebug() << "qrombdeltac inner" << j << dssp << pq );
            if ( abs(dssp)<(eps*abs(pq)) ) break; /* goto 99 */  //Zupq1=2302
            if ( fabs(pq) < 1.0e-100 ) break;  // nicht im Pascal-Programm
        }/*3*/  //Zupq1=2303
        sp[j+1] = sp[j];  //Zupq1=2304
        hp[j+1] = 0.25*hp[j];  //Zupq1=2305
    }/*2*/  //Zupq1=2306

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



#ifdef __CUDACC__
__host__ __device__
#endif
//void SasCalc_GENERIC_calculation::qrombchid( double l, double r, double sigma, double dbeta, double delta, double theta, double phi,
//                          double qx, double qy, double qz, double p11, double p12, double p13, double p21,
//                          double p22, double p23, double p31, double p32, double p33,
//                          double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n, double ax3n,
//                          double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z, double ax3x,
//                          double ax3y, double ax3z, double sigx, double sigy, double sigz,
//                          int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
//                          double *carr1, // NEU
//                          double &pq ) const
void SasCalc_GENERIC_calculation::qrombchid( double l, double r, double p1, double sigma, double alfa, double delta,
                    double theta, double phi, double qx, double qy, double qz, double p11, double p12, double p13,
                    double p21, double p22, double p23, double p31, double p32, double p33, double qxn, double qyn,
                    double qzn, double qhkl, double ax1n, double ax2n, double ax3n, double ax1x, double ax1y,
                    double ax1z, double ax2x, double ax2y, double ax2z, double ax3x, double ax3y, double ax3z,
                    double sigx, double sigy, double sigz,
                    int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
                    double *carr1, double &pq ) const
{/*1*/  //Zupq1=53

    /* label 99; */  //Zupq1=55
    const double eps = 1.0e-4;  //Zupq1=57
    //const double eps1 = 1.0e-8;  //Zupq1=58
    //const double epschi =  0.001;  //Zupq1=59
    const int    jmax = 20;  //Zupq1=60
    //const int    jmaxp = 21;  //Zupq1=61
    const int    k = 5;  //Zupq1=62
    //const int    np = 5;  //Zupq1=63
    //const double prmin =  0.00001;       /*  minimum pr(delta)-value  */  //Zupq1=64

    int i, j;  //Zupq1=69
    //int trapzdit;  //Zupq1=70
    double dssp; //, delmax;  //Zupq1=71
    RealArrayJMAXP hp, sp;  //Zupq1=72
    RealArrayNP    cp, dp;  //Zupq1=73
    float alim, blim;       //Zupq1=74

    int trapzdchid_cnt;

    //qDebug() << "qrombchid(" << r << sigma << dbeta << delta << theta << phi << qx << qy << qz << p11 << p12 << p13 << p21
    //         << p22 << p23 << p31 << p32 << p33 << qxn << qyn << qzn << qhkl << ax1n << ax2n << ax3n << ax1x << ax1y
    //         << ax1z << ax2x << ax2y << ax2z << ax3x << ax3y << ax3z << sigx << sigy << sigz << "ordis" << ordis << dim
    //         << i0 << i1 << i2 << i3 << i4 << "&pq" << " );";

    /*begin*/  //Zupq1=1420
    alim = 0.0;  //Zupq1=1421
    blim = 2*M_PI;  //Zupq1=1422

    if ( i1==12 )
    {/*2*/       /*  for cubes  */  //Zupq1=1424
        alim = eps;  //Zupq1=1425
        blim = M_PI/2.0-eps;  //Zupq1=1426
    }/*2*/  //Zupq1=1427

    if ( (i1==13) || (i1==14) )
    {/*2*/       /*  for ellipsoids  */  //Zupq1=1429
        alim = 0;  //Zupq1=1430
        blim = M_PI/2.0;  //Zupq1=1431
    }/*2*/  //Zupq1=1432

    if ( i1==15 )
    {/*2*/           /*  for area integration  */  //Zupq1=1434
        alim = 0.0;  //Zupq1=1435
        blim = l;  //Zupq1=1436
    }/*2*/  //Zupq1=1437

    if ( i1==17 )
    {/*2*/        /*  superball integration  */  //Zupq1=1439
        alim = 0;  //Zupq1=1440
        /* blim:=r*power(1-power(delta/r,p1),1/p1)-0.0005;  //Zupq1=1441 */
        blim = M_PI/2.0;  //Zupq1=1442
    }/*2*/  //Zupq1=1443

    hp[1] = 1.0;  //Zupq1=1449
    for ( j=1; j<=jmax; j++ )
    {/*2*/  //Zupq1=1450
        /* trapzdchid(alim,blim,l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,dbeta,phi,theta,qxn,qyn,qzn,q,maxit,i0,i1,i2,i3,i4,sp^[j],j);  //Zupq1=1451 */
        trapzdchid(alim,blim,l,r,p1,sigma,alfa,delta,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,
                   p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,
                   ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,sp[j],j, trapzdchid_cnt );  //Zupq1=1452
        if ( j>=k )
        {/*3*/  //Zupq1=1453
            for ( i=1; i<=k; i++ )
            {/*4*/  //Zupq1=1454
                cp[i] = hp[j-k+i];  //Zupq1=1455
                dp[i] = sp[j-k+i];  //Zupq1=1456
            }/*4*/  //Zupq1=1457
            polint(cp,dp,k,0.0,pq,dssp,"qrombchid");  //Zupq1=1458
            if ( abs(dssp)<(eps*abs(pq)) ) break; /* goto 99*/  //Zupq1=1459
            if ( fabs(pq) < 1.0e-300 )
            {
                //qDebug() << " ... qrombchid" << dssp << pq;
                break; // TODO ?
            }
        }/*3*/  //Zupq1=1460
        sp[j+1] = sp[j];  //Zupq1=1461
        hp[j+1] = 0.25*hp[j];  //Zupq1=1462
    }/*2*/  //Zupq1=1463
    /*99:*/  //Zupq1=1464
} /* qrombchid() */



#ifdef __CUDACC__
__host__ __device__
#endif
//void SasCalc_GENERIC_calculation::trapzdchid( double a, double b, double l, double r, double p1, double sigma, double alfa, double dbeta,
//                                            double delta, double theta, double phi, double qx, double qy, double qz, double p11,
//                                            double p12, double p13, double p21, double p22, double p23, double p31, double p32,
//                                            double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n,
//                                            double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
//                                            double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
//                                            int ordis, int dim, int i0, int i1, int i2, int i3, int i4,
//                                            double *carr1, double &pq, int n,
//                                            int &trapzdchid_cnt ) const
void SasCalc_GENERIC_calculation::trapzdchid( double a, double b, double l, double r, double p1, double sigma, double alfa,
               double delta, double theta, double phi, double qx, double qy, double qz, double p11,
               double p12, double p13, double p21, double p22, double p23, double p31, double /*p32*/,
               double p33, double qxn, double qyn, double qzn, double qhkl, double ax1n, double ax2n,
               double ax3n, double ax1x, double ax1y, double ax1z, double ax2x, double ax2y, double ax2z,
               double ax3x, double ax3y, double ax3z, double sigx, double sigy, double sigz,
               int /*ordis*/, int /*dim*/, int /*i0*/, int i1, int /*i2*/, int i3, int i4,
               double *carr1, double &pq, int n,
               int &trapzdchid_cnt ) const
{/*1*/  //Zupq1=131
    const double eps =  0.0000000001;  //Zupq1=133

    /*var*/ int i, j;  //Zupq1=138
    double x, tnm, sump, /*sumf, sumn, sums,*/ del, /*pqsum, oldpqsum,*/ delser, argser;  //Zupq1=139
    double /*fa, fb, fx,*/ pa, pb, px, /*na, nb, nx, sa, sb, sx,*/ epsi, aa, bb, cc, oldpa, oldpb, oldpx;  //Zupq1=140
    double z, a1, /*x1z,*/ arga, argb, argx, l1, l2, l3; //, qxhkl, qyhkl, qzhkl;  //Zupq1=141
    double mla1, mla2, /*mla3,*/ mlb1, mlb2, /*mlb3,*/ mlx1, mlx2; //, mlx3;  //Zupq1=142
    //double ma11, ma12, ma13, ma21, ma22, ma23, ma31, ma32, ma33;  //Zupq1=143
    //double mb11, mb12, mb13, mb21, mb22, mb23, mb31, mb32, mb33;  //Zupq1=144
    //double mx11, mx12, mx13, mx21, mx22, mx23, mx31, mx32, mx33;  //Zupq1=145
    double dqxa, dqya, dqza, dqxb, dqyb, dqzb, dqs1a, dqs2a, dqs3a, dqs1b, dqs2b, dqs3b;  //Zupq1=146
    double dqxx, dqyx, dqzx, dqs1x, dqs2x, dqs3x; //, vola, volb, volx;  //Zupq1=147
    double r11a, r12a, r13a, r21a, r22a, r23a, r31a, r32a, r33a;  //Zupq1=148
    double r11b, r12b, r13b, r21b, r22b, r23b, r31b, r32b, r33b;  //Zupq1=149
    double r11x, r12x, r13x, r21x, r22x, r23x, r31x, r32x, r33x;  //Zupq1=150
    double ax1xa, ax1ya, ax1za, ax1xb, ax1yb, ax1zb, ax1xx, ax1yx, ax1zx;  //Zupq1=151
    double ax2xa, ax2ya, ax2za, ax2xb, ax2yb, ax2zb, ax2xx, ax2yx, ax2zx;  //Zupq1=152
    double ax3xa, ax3ya, ax3za, ax3xb, ax3yb, ax3zb, ax3xx, ax3yx, ax3zx;  //Zupq1=153
    double qn, qxhkla, qyhkla, qzhkla, qxhklb, qyhklb, qzhklb, qxhklx, qyhklx, qzhklx;  //Zupq1=154
    //double aexa, aeya, aeza, bexa, beya, beza, cexa, ceya, ceza, aexb, aeyb, aezb, bexb, beyb, bezb, cexb, ceyb, cezb;  //Zupq1=155
    //double aexx, aeyx, aezx, bexx, beyx, bezx, cexx, ceyx, cezx;  //Zupq1=156
    double pa1, pa2, pa3, pb1, pb2, pb3, px1, px2, px3, ella, ellb, ellc;  //Zupq1=157
    double qxl, qyl, qnna, qnnb, qnnx; //, argas, argbs, z1, z2, z3, z4, z5, z6;  //Zupq1=158
    double arga1, arga2, arga3, argb1, argb2, argb3, argx1, argx2, argx3;  //Zupq1=159

    arga = 1; // to avoid compiler warnings
    argb = 1;

    if ( n==1 )
    {/*2*/  //Zupq1=163
        if ( i1==1 )
        {/*3*/    /*  cylinder, general case  */  //Zupq1=164
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=165
            a1 = 1/(2.0*z*(z-1));  //Zupq1=166
            mla1 = p11*cos(a)*sin(delta)+p12*sin(a)*sin(delta)+p13*cos(delta);  //Zupq1=167
            mla2 = p21*sin(a)*sin(delta)+p22*cos(a)*sin(delta)+p23*cos(delta);  //Zupq1=168
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta);  //Zupq1=169 */
            mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=170
            mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=171
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta);  //Zupq1=172 */
            arga = (qx*mla1+qy*mla2+eps)*l/(z+1);  //Zupq1=173
            argb = (qx*mlb1+qy*mlb2+eps)*l/(z+1);  //Zupq1=174
            if ( i3==0 )
            {/*4*/ /*  P(q)  */  //Zupq1=175
                a1 = 1/(2.0*z*(z-1));  //Zupq1=176
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=177
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=178
            }/*4*/  //Zupq1=179
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=180
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=181
                pa = pa1*pa1;  //Zupq1=182
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=183
                pb = pb1*pb1;  //Zupq1=184
            }/*4*/  //Zupq1=185
        }/*3*/  //Zupq1=186
        if ( i1==2 )
        {/*3*/   /*  cylinder, x-axis  */  //Zupq1=187
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=188
            arga = (qx*cos(delta)-qy*sin(a)*sin(delta)+eps)*l/(z+1);  //Zupq1=189
            argb = (qx*cos(delta)-qy*sin(b)*sin(delta)+eps)*l/(z+1);  //Zupq1=190
            if ( i3==0 )
            {/*4*/ /*  P(q)  */  //Zupq1=191
                a1 = 1/(2.0*z*(z-1));  //Zupq1=192
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=193
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=194
            }/*4*/  //Zupq1=195
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=196
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=197
                pa = pa1*pa1;  //Zupq1=198
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=199
                pb = pb1*pb1;  //Zupq1=200
            }/*4*/  //Zupq1=201
        }/*3*/  //Zupq1=202
        if ( i1==3 )
        {/*3*/   /*  cylinder, y-axis  */  //Zupq1=203
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=204
            arga = (qx*sin(a)*sin(delta)+qy*cos(delta)+eps)*l/(z+1);  //Zupq1=205
            argb = (qx*sin(b)*sin(delta)+qy*cos(delta)+eps)*l/(z+1);  //Zupq1=206
            if ( i3==0 )
            {/*4*/ /*  P(q)  */  //Zupq1=207
                a1 = 1/(2.0*z*(z-1));  //Zupq1=208
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=209
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=210
            }/*4*/  //Zupq1=211
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=212
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=213
                pa = pa1*pa1;  //Zupq1=214
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=215
                pb = pb1*pb1;  //Zupq1=216
            }/*4*/  //Zupq1=217
        }/*3*/  //Zupq1=218
        if ( i1==4 )
        {/*3*/   /*  cylinder, -z-axis  */  //Zupq1=219
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=220
            arga = (qx*sin(a)*sin(delta)-qy*cos(a)*sin(delta)+eps)*l/(z+1);  //Zupq1=221
            argb = (qx*sin(b)*sin(delta)-qy*cos(b)*sin(delta)+eps)*l/(z+1);  //Zupq1=222
            if ( i3==0 )
            {/*4*/ /*  P(q)  */  //Zupq1=223
                a1 = 1/(2.0*z*(z-1));  //Zupq1=224
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=225
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=226
            }/*4*/  //Zupq1=227
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=228
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=229
                pa = pa1*pa1;  //Zupq1=230
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=231
                pb = pb1*pb1;  //Zupq1=232
            }/*4*/  //Zupq1=233
        }/*3*/  //Zupq1=234
        if ( i1==5 )
        {/*3*/    /*  general series expansion  */  //Zupq1=235
            mla1 = p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=236
            mla2 = p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=237
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta);  //Zupq1=238 */
            mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=239
            mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=240
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta);  //Zupq1=241 */
            pa = pow(mla1,i3)*pow(mla2,i4);  //Zupq1=242
            pb = pow(mlb1,i3)*pow(mlb2,i4);  //Zupq1=243
        }/*3*/  //Zupq1=244
        if ( i1==6 )
        {/*3*/    /*  unit cell rotation  */  //Zupq1=245
            dqxa = qx-qhkl*(p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));  //Zupq1=246
            dqya = qy-qhkl*(p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));  //Zupq1=247
            dqza = qz-qhkl*(p31*cos(a)*sin(delta)+p33*cos(delta));  //Zupq1=248
            dqxb = qx-qhkl*(p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));  //Zupq1=249
            dqyb = qy-qhkl*(p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));  //Zupq1=250
            dqzb = qz-qhkl*(p31*cos(b)*sin(delta)+p33*cos(delta));  //Zupq1=251
            dqs1a = (dqxa*ax1x+dqya*ax1y+dqza*ax1z)/(ax1n*sigx);;  //Zupq1=252
            dqs2a = (dqxa*ax2x+dqya*ax2y+dqza*ax2z)/(ax2n*sigy);  //Zupq1=253
            dqs3a = (dqxa*ax3x+dqya*ax3y+dqza*ax3z)/(ax3n*sigz);  //Zupq1=254
            dqs1b = (dqxb*ax1x+dqyb*ax1y+dqzb*ax1z)/(ax1n*sigx);  //Zupq1=255
            dqs2b = (dqxb*ax2x+dqyb*ax2y+dqzb*ax2z)/(ax2n*sigy);  //Zupq1=256
            dqs3b = (dqxb*ax3x+dqyb*ax3y+dqzb*ax3z)/(ax3n*sigz);  //Zupq1=257
            pa = exp(-4*(dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a)/M_PI);  //Zupq1=258
            pb = exp(-4*(dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b)/M_PI);  //Zupq1=259
        }/*3*/  //Zupq1=260
        if ( i1==7 )
        {/*3*/    /*  fiber rotation  */  //Zupq1=261
            /*  rotation axis, director  */  //Zupq1=262
            l1 = sin(theta)*cos(phi);  //Zupq1=263
            l2 = sin(theta)*sin(phi);  //Zupq1=264
            l3 = -cos(theta);  //Zupq1=265

            /*  rotation matrix Ra  */  //Zupq1=267
            r11a = cos(a)+(1-cos(a))*l1*l1;  //Zupq1=268
            r12a = -l3*sin(a)+(1-cos(a))*l1*l2;  //Zupq1=269
            r13a = -l2*sin(a)+(1-cos(a))*l1*l3;  //Zupq1=270
            r21a = l3*sin(a)+(1-cos(a))*l1*l2;  //Zupq1=271
            r22a = cos(a)+(1-cos(a))*l2*l2;  //Zupq1=272
            r23a = l1*sin(a)+(1-cos(a))*l2*l3;  //Zupq1=273
            r31a = l2*sin(a)+(1-cos(a))*l1*l3;  //Zupq1=274
            r32a = -l1*sin(a)+(1-cos(a))*l2*l3;  //Zupq1=275
            r33a = cos(a)+(1-cos(a))*l3*l3;  //Zupq1=276

            /*  rotation matrix Rb  */  //Zupq1=278
            r11b = cos(b)+(1-cos(b))*l1*l1;  //Zupq1=279
            r12b = -l3*sin(b)+(1-cos(b))*l1*l2;  //Zupq1=280
            r13b = -l2*sin(b)+(1-cos(b))*l1*l3;  //Zupq1=281
            r21b = l3*sin(b)+(1-cos(b))*l1*l2;  //Zupq1=282
            r22b = cos(b)+(1-cos(b))*l2*l2;  //Zupq1=283
            r23b = l1*sin(b)+(1-cos(b))*l2*l3;  //Zupq1=284
            r31b = l2*sin(b)+(1-cos(b))*l1*l3;  //Zupq1=285
            r32b = -l1*sin(b)+(1-cos(b))*l2*l3;  //Zupq1=286
            r33b = cos(b)+(1-cos(b))*l3*l3;  //Zupq1=287

            /*  rotate scattering vector  */  //Zupq1=294
            qxhkla = r11a*qxn+r12a*qyn+r13a*qzn;  //Zupq1=295
            qyhkla = r21a*qxn+r22a*qyn+r23a*qzn;  //Zupq1=296
            qzhkla = r31a*qxn+r32a*qyn+r33a*qzn;  //Zupq1=297
            qxhklb = r11b*qxn+r12b*qyn+r13b*qzn;  //Zupq1=298
            qyhklb = r21b*qxn+r22b*qyn+r23b*qzn;  //Zupq1=299
            qzhklb = r31b*qxn+r32b*qyn+r33b*qzn;  //Zupq1=300

            dqxa = qx-qxhkla;  //Zupq1=378
            dqya = qy-qyhkla;  //Zupq1=379
            dqza = qz-qzhkla;  //Zupq1=380

            dqxb = qx-qxhklb;  //Zupq1=382
            dqyb = qy-qyhklb;  //Zupq1=383
            dqzb = qz-qzhklb;  //Zupq1=384

            ax1xa = (r11a*ax1x+r12a*ax1y+r13a*ax1z);  //Zupq1=393
            ax1ya = (r21a*ax1x+r22a*ax1y+r23a*ax1z);  //Zupq1=394
            ax1za = (r31a*ax1x+r32a*ax1y+r33a*ax1z);  //Zupq1=395
            ax1xb = (r11b*ax1x+r12b*ax1y+r13b*ax1z);  //Zupq1=396
            ax1yb = (r21b*ax1x+r22b*ax1y+r23b*ax1z);  //Zupq1=397
            ax1zb = (r31b*ax1x+r32b*ax1y+r33b*ax1z);  //Zupq1=398
            ax2xa = (r11a*ax2x+r12a*ax2y+r13a*ax2z);  //Zupq1=399
            ax2ya = (r21a*ax2x+r22a*ax2y+r23a*ax2z);  //Zupq1=400
            ax2za = (r31a*ax2x+r32a*ax2y+r33a*ax2z);  //Zupq1=401
            ax2xb = (r11b*ax2x+r12b*ax2y+r13b*ax2z);  //Zupq1=402
            ax2yb = (r21b*ax2x+r22b*ax2y+r23b*ax2z);  //Zupq1=403
            ax2zb = (r31b*ax2x+r32b*ax2y+r33b*ax2z);  //Zupq1=404
            ax3xa = (r11a*ax3x+r12a*ax3y+r13a*ax3z);  //Zupq1=405
            ax3ya = (r21a*ax3x+r22a*ax3y+r23a*ax3z);  //Zupq1=406
            ax3za = (r31a*ax3x+r32a*ax3y+r33a*ax3z);  //Zupq1=407
            ax3xb = (r11b*ax3x+r12b*ax3y+r13b*ax3z);  //Zupq1=408
            ax3yb = (r21b*ax3x+r22b*ax3y+r23b*ax3z);  //Zupq1=409
            ax3zb = (r31b*ax3x+r32b*ax3y+r33b*ax3z);  //Zupq1=410

            dqs1a = (dqxa*ax1xa+dqya*ax1ya+dqza*ax1za)/(ax1n*sigx);  //Zupq1=412
            dqs2a = (dqxa*ax2xa+dqya*ax2ya+dqza*ax2za)/(ax2n*sigy);  //Zupq1=413
            dqs3a = (dqxa*ax3xa+dqya*ax3ya+dqza*ax3za)/(ax3n*sigz);  //Zupq1=414
            dqs1b = (dqxb*ax1xb+dqyb*ax1yb+dqzb*ax1zb)/(ax1n*sigx);  //Zupq1=415
            dqs2b = (dqxb*ax2xb+dqyb*ax2yb+dqzb*ax2zb)/(ax2n*sigy);  //Zupq1=416
            dqs3b = (dqxb*ax3xb+dqyb*ax3yb+dqzb*ax3zb)/(ax3n*sigz);  //Zupq1=417

            arga = dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a;  //Zupq1=426
            argb = dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b;  //Zupq1=427
            pa = exp(-4*arga/M_PI);  //Zupq1=428
            pb = exp(-4*argb/M_PI);  //Zupq1=429
            /* if (arga>2) then pa:=eps  //Zupq1=430 */
            /*    else pa:=exp(-4*arga/pi);  //Zupq1=431 */
            /* if (argb>2) then pb:=eps  //Zupq1=432 */
            /*    else pb:=exp(-4*argb/pi);  //Zupq1=433 */
        }/*3*/  //Zupq1=434
        if ( i1==8 )
        {/*3*/    /*  disk, general case  */  //Zupq1=435
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=436
            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=437
            qxl = qx/qn;  //Zupq1=438
            qyl = qy/qn;  //Zupq1=439
            mla1 = p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=440
            mla2 = p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=441
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta);  //Zupq1=442 */
            mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=443
            mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=444
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta);  //Zupq1=445 */
            qnna = (qxl*mla1+qyl*mla2);  //Zupq1=446
            qnnb = (qxl*mlb1+qyl*mlb2);  //Zupq1=447
            arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=448
            argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=449

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=451
                if ( arga<0.015 )
                {/*5*/  //Zupq1=452
                    pa = 1;  //Zupq1=453
                    oldpa = 0;  //Zupq1=454
                    argser = 1;  //Zupq1=455
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=456
                        argser = argser*arga*arga/4.0;  //Zupq1=457
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=458 */
                        pa = pa+carr1[i]*argser;  //Zupq1=459
                        delser = fabs((pa-oldpa)/pa);  //Zupq1=460
                        if ( delser<0.0001 ) break; /* goto 12; */  //Zupq1=461
                        oldpa = pa;  //Zupq1=462
                    }/*6*/  //Zupq1=463
                    /*12:*/  //Zupq1=464
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=465
                }/*5*/  //Zupq1=466
                else
                {/*5*/  //Zupq1=467
                    if ( i3==0 )
                    {/*6*/ /*  P(q)  */  //Zupq1=468
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=469
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=470
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=471
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=472
                    }/*6*/  //Zupq1=473
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=474
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=475
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=476
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=477
                        pa = pa3*pa3;  //Zupq1=478
                    }/*6*/  //Zupq1=479
                }/*5*/  //Zupq1=480
                if ( argb<0.015 )
                {/*5*/  //Zupq1=481
                    pb = 1;  //Zupq1=482
                    oldpb = 0;  //Zupq1=483
                    argser = 1;  //Zupq1=484
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=485
                        argser = argser*argb*argb/4.0;  //Zupq1=486
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=487 */
                        pb = pb+carr1[i]*argser;  //Zupq1=488
                        delser = fabs((pb-oldpb)/pb);  //Zupq1=489
                        if ( delser<0.0001 ) break; /* goto 13; */  //Zupq1=490
                        oldpb = pb;  //Zupq1=491
                    }/*6*/  //Zupq1=492
                    /*13:*/  //Zupq1=493
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=494
                }/*5*/  //Zupq1=495
                else
                {/*5*/  //Zupq1=496
                    if ( i3==0 )
                    {/*6*/  /*  P(q)  */  //Zupq1=497
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=498
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=499
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=500
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=501
                    }/*6*/  //Zupq1=502
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=503
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=504
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=505
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=506
                        pb = pb3*pb3;  //Zupq1=507
                    }/*6*/  //Zupq1=508
                }/*5*/  //Zupq1=509
            }/*4*/  //Zupq1=510
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=511
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=512
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=513
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=514
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=515
            }/*4*/  //Zupq1=516
        }/*3*/  //Zupq1=517

        if ( i1==9 )
        {/*3*/   /*  disk, x-axis  */  //Zupq1=519
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=520
            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=521
            qxl = qx/qn;  //Zupq1=522
            qyl = qy/qn;  //Zupq1=523
            qnna = qxl*cos(delta)-qyl*sin(a)*sin(delta);  //Zupq1=524
            qnnb = qxl*cos(delta)-qyl*sin(b)*sin(delta);  //Zupq1=525
            arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=526
            argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=527

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=529
                if ( arga<0.015 )
                {/*5*/  //Zupq1=530
                    pa = 1;  //Zupq1=531
                    oldpa = 0;  //Zupq1=532
                    argser = 1;  //Zupq1=533
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=534
                        argser = argser*arga*arga/4.0;  //Zupq1=535
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=536 */
                        pa = pa+carr1[i]*argser;  //Zupq1=537
                        delser = fabs((pa-oldpa)/pa);  //Zupq1=538
                        if ( delser<0.0001 ) break; /* goto 15; */  //Zupq1=539
                        oldpa = pa;  //Zupq1=540
                    }/*6*/  //Zupq1=541
                    /*15:*/  //Zupq1=542
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=543
                }/*5*/  //Zupq1=544
                else
                {/*5*/  //Zupq1=545
                    if ( i3==0 )
                    {/*6*/ /*  P(q)  */  //Zupq1=546
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=547
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=548
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=549
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=550
                    }/*6*/  //Zupq1=551
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=552
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=553
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=554
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=555
                        pa = pa3*pa3;  //Zupq1=556
                    }/*6*/  //Zupq1=557
                }/*5*/  //Zupq1=558
                if ( argb<0.015 )
                {/*5*/  //Zupq1=559
                    pb = 1;  //Zupq1=560
                    oldpb = 0;  //Zupq1=561
                    argser = 1;  //Zupq1=562
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=563
                        argser = argser*argb*argb/4.0;  //Zupq1=564
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=565 */
                        pb = pb+carr1[i]*argser;  //Zupq1=566
                        delser = fabs((pb-oldpb)/pb);  //Zupq1=567
                        if ( delser<0.0001 ) break; /* goto 16; */  //Zupq1=568
                        oldpb = pb;  //Zupq1=569
                    }/*6*/  //Zupq1=570
                    /*16:*/  //Zupq1=571
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=572
                }/*5*/  //Zupq1=573
                else
                {/*5*/  //Zupq1=574
                    if ( i3==0 )
                    {/*6*/  /*  P(q)  */  //Zupq1=575
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=576
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=577
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=578
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=579
                    }/*6*/  //Zupq1=580
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=581
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=582
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=583
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=584
                        pb = pb3*pb3;  //Zupq1=585
                    }/*6*/  //Zupq1=586
                }/*5*/  //Zupq1=587
            }/*4*/  //Zupq1=588
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=589
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=590
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=591
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=592
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=593
            }/*4*/  //Zupq1=594
        }/*3*/  //Zupq1=595

        if ( i1==10 )
        {/*3*/   /*  disk, y-axis  */  //Zupq1=597
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=598
            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=599
            qxl = qx/qn;  //Zupq1=600
            qyl = qy/qn;  //Zupq1=601
            qnna = qxl*sin(a)*sin(delta)+qyl*cos(delta);  //Zupq1=602
            qnnb = qxl*sin(b)*sin(delta)+qyl*cos(delta);  //Zupq1=603
            arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=604
            argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=605

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=607
                if ( arga<0.015 )
                {/*5*/  //Zupq1=608
                    pa = 1;  //Zupq1=609
                    oldpa = 0;  //Zupq1=610
                    argser = 1;  //Zupq1=611
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=612
                        argser = argser*arga*arga/4.0;  //Zupq1=613
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=614 */
                        pa = pa+carr1[i]*argser;  //Zupq1=615
                        delser = fabs((pa-oldpa)/pa);  //Zupq1=616
                        if ( delser<0.0001 ) break; /* goto 18; */  //Zupq1=617
                        oldpa = pa;  //Zupq1=618
                    }/*6*/  //Zupq1=619
                    /*18:*/  //Zupq1=620
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=621
                }/*5*/  //Zupq1=622
                else
                {/*5*/  //Zupq1=623
                    if ( i3==0 )
                    {/*6*/ /*  P(q)  */  //Zupq1=624
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=625
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=626
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=627
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=628
                    }/*6*/  //Zupq1=629
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=630
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=631
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=632
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=633
                        pa = pa3*pa3;  //Zupq1=634
                    }/*6*/  //Zupq1=635
                }/*5*/  //Zupq1=636
                if ( argb<0.015 )
                {/*5*/  //Zupq1=637
                    pb = 1;  //Zupq1=638
                    oldpb = 0;  //Zupq1=639
                    argser = 1;  //Zupq1=640
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=641
                        argser = argser*argb*argb/4.0;  //Zupq1=642
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=643 */
                        pb = pb+carr1[i]*argser;  //Zupq1=644
                        delser = fabs((pb-oldpb)/pb);  //Zupq1=645
                        if ( delser<0.0001 ) break; /* goto 19; */  //Zupq1=646
                        oldpb = pb;  //Zupq1=647
                    }/*6*/  //Zupq1=648
                    /*19:*/  //Zupq1=649
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=650
                }/*5*/  //Zupq1=651
                else
                {/*5*/  //Zupq1=652
                    if ( i3==0 )
                    {/*6*/  /*  P(q)  */  //Zupq1=653
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=654
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=655
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=656
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=657
                    }/*6*/  //Zupq1=658
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=659
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=660
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=661
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=662
                        pb = pb3*pb3;  //Zupq1=663
                    }/*6*/  //Zupq1=664
                }/*5*/  //Zupq1=665
            }/*4*/  //Zupq1=666
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=667
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=668
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=669
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=670
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=671
            }/*4*/  //Zupq1=672
        }/*3*/  //Zupq1=673

        if ( i1==11 )
        {/*3*/   /*  disk, -z-axis  */  //Zupq1=675
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=676
            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=677
            qxl = qx/qn;  //Zupq1=678
            qyl = qy/qn;  //Zupq1=679
            qnna = qxl*sin(a)*sin(delta)-qyl*cos(a)*sin(delta);  //Zupq1=680
            qnnb = qxl*sin(b)*sin(delta)-qyl*cos(b)*sin(delta);  //Zupq1=681
            arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=682
            argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=683

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=685
                if ( arga<0.015 )
                {/*5*/  //Zupq1=686
                    pa = 1;  //Zupq1=687
                    oldpa = 0;  //Zupq1=688
                    argser = 1;  //Zupq1=689
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=690
                        argser = argser*arga*arga/4.0;  //Zupq1=691
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=692 */
                        pa = pa+carr1[i]*argser;  //Zupq1=693
                        delser = fabs((pa-oldpa)/pa);  //Zupq1=694
                        if ( delser<0.0001 ) break; /* goto 21; */  //Zupq1=695
                        oldpa = pa;  //Zupq1=696
                    }/*6*/  //Zupq1=697
                    /*21:*/  //Zupq1=698
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=699
                }/*5*/  //Zupq1=700
                else
                {/*5*/  //Zupq1=701
                    if ( i3==0 )
                    {/*6*/ /*  P(q)  */  //Zupq1=702
                        pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=703
                        pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=704
                        pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=705
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=706
                    }/*6*/  //Zupq1=707
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=708
                        pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=709
                        pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=710
                        pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=711
                        pa = pa3*pa3;  //Zupq1=712
                    }/*6*/  //Zupq1=713
                }/*5*/  //Zupq1=714
                if ( argb<0.015 )
                {/*5*/  //Zupq1=715
                    pb = 1;  //Zupq1=716
                    oldpb = 0;  //Zupq1=717
                    argser = 1;  //Zupq1=718
                    for ( i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=719
                        argser = argser*argb*argb/4.0;  //Zupq1=720
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=721 */
                        pb = pb+carr1[i]*argser;  //Zupq1=722
                        delser = fabs((pb-oldpb)/pb);  //Zupq1=723
                        if ( delser<0.0001 ) break; /* goto 22; */  //Zupq1=724
                        oldpb = pb;  //Zupq1=725
                    }/*6*/  //Zupq1=726
                    /*22:*/  //Zupq1=727
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=728
                }/*5*/  //Zupq1=729
                else
                {/*5*/  //Zupq1=730
                    if ( i3==0 )
                    {/*6*/  /*  P(q)  */  //Zupq1=731
                        pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=732
                        pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=733
                        pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=734
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=735
                    }/*6*/  //Zupq1=736
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=737
                        pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=738
                        pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=739
                        pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=740
                        pb = pb3*pb3;  //Zupq1=741
                    }/*6*/  //Zupq1=742
                }/*5*/  //Zupq1=743
            }/*4*/  //Zupq1=744
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=745
                pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=746
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=747
                pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=748
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=749
            }/*4*/  //Zupq1=750
        }/*3*/  //Zupq1=751

        if ( i1==12 )
        {/*3*/   /*  isotropic cube  */  //Zupq1=753
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=754

            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=787
            arga1 = qn*sin(delta)*cos(a)*r/(z+1)+eps;  //Zupq1=788
            arga2 = qn*sin(delta)*sin(a)*r/(z+1)+eps;  //Zupq1=789
            arga3 = qn*cos(delta)*r/(z+1)+eps;  //Zupq1=790
            if ( i3==0 )
            {/*4*/  /*  P(q)  */  //Zupq1=791
                a1 = 1/(2.0*z*(z-1));  //Zupq1=792
                if ( arga1 < 0.001 )
                    pa1 = 1;
                else  //Zupq1=793
                    pa1 = (a1/(arga1*arga1+eps))*(1-cos((z-1)*atan(2.0*arga1))/pow(1.0+4*arga1*arga1,(z-1)/2.0));  //Zupq1=794
                if ( arga2 < 0.001 )
                    pa2 = 1;
                else  //Zupq1=795
                    pa2 = (a1/(arga2*arga2+eps))*(1-cos((z-1)*atan(2.0*arga2))/pow(1.0+4*arga2*arga2,(z-1)/2.0));  //Zupq1=796
                if ( arga3 < 0.001 )
                    pa3 = 1;
                else  //Zupq1=797
                    pa3 = (a1/(arga3*arga3+eps))*(1-cos((z-1)*atan(2.0*arga3))/pow(1.0+4*arga3*arga3,(z-1)/2.0));  //Zupq1=798
                pa = pa1*pa2*pa3;  //Zupq1=799
            }/*4*/  //Zupq1=800
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=801
                a1 = 1/z;  //Zupq1=802
                if ( arga1 < 0.001 ) pa1 = 1; else  //Zupq1=803
                        pa1 = (a1/(arga1+eps))*sin(z*atan(arga1))/pow(1.0+arga1*arga1,z/2.0);  //Zupq1=804
                if ( arga2 < 0.001 ) pa2 = 1; else  //Zupq1=805
                        pa2 = (a1/(arga2+eps))*sin(z*atan(arga2))/pow(1.0+arga2*arga2,z/2.0);  //Zupq1=806
                if ( arga3 < 0.001 ) pa3 = 1; else  //Zupq1=807
                        pa3 = (a1/(arga3+eps))*sin(z*atan(arga3))/pow(1.0+arga3*arga3,z/2.0);  //Zupq1=808
                pa = pa1*pa1*pa2*pa2*pa3*pa3;  //Zupq1=809
            }/*4*/  //Zupq1=810
            argb1 = qn*sin(delta)*cos(b)*r/(z+1);  //Zupq1=811
            argb2 = qn*sin(delta)*sin(b)*r/(z+1);  //Zupq1=812
            argb3 = qn*cos(delta)*r/(z+1);  //Zupq1=813
            if ( i3==0 )
            {/*4*/   /*  P(q)  */  //Zupq1=814
                a1 = 1/(2.0*z*(z-1));  //Zupq1=815
                if ( argb1 < 0.001 ) pb1 = 1; else  //Zupq1=816
                        pb1 = (a1/(argb1*argb1+eps))*(1-cos((z-1)*atan(2.0*argb1))/pow(1.0+4*argb1*argb1,(z-1)/2.0));  //Zupq1=817
                if ( argb2 < 0.001 ) pb2 = 1; else  //Zupq1=818
                        pb2 = (a1/(argb2*argb2+eps))*(1-cos((z-1)*atan(2.0*argb2))/pow(1.0+4*argb2*argb2,(z-1)/2.0));  //Zupq1=819
                if ( argb3 < 0.001 ) pb3 = 1; else  //Zupq1=820
                        pb3 = (a1/(argb3*argb3+eps))*(1-cos((z-1)*atan(2.0*argb3))/pow(1.0+4*argb3*argb3,(z-1)/2.0));  //Zupq1=821
                pb = pb1*pb2*pb3;  //Zupq1=822
            }/*4*/  //Zupq1=823
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=824
                a1 = 1/z;  //Zupq1=825
                if ( argb1 < 0.001 ) pb1 = 1; else  //Zupq1=826
                        pb1 = (a1/(argb1+eps))*sin(z*atan(argb1))/pow(1.0+argb1*argb1,z/2.0);  //Zupq1=827
                if ( argb2 < 0.001 ) pb2 = 1; else  //Zupq1=828
                        pb2 = (a1/(argb2+eps))*sin(z*atan(argb2))/pow(1.0+argb2*argb2,z/2.0);  //Zupq1=829
                if ( arga3 < 0.001 ) pb3 = 1; else  //Zupq1=830
                        pb3 = (a1/(argb3+eps))*sin(z*atan(argb3))/pow(1.0+argb3*argb3,z/2.0);  //Zupq1=831
                pb = pb1*pb1*pb2*pb2*pb3*pb3;  //Zupq1=832
            }/*4*/  //Zupq1=833
        }/*3*/  //Zupq1=834
        if ( i1==13 )
        {/*3*/   /*  biaxial ellipsoid, isotropic  */  //Zupq1=835
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=836
            epsi = l/r;  //Zupq1=837
            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=838
            arga = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(a)*cos(a))/(z+1);  //Zupq1=839
            a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=840
            pa1 = a1*pow(arga,-4)*(1+cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0));  //Zupq1=841
            pa2 = (a1/(z-4))*pow(arga,-5)*sin((z-4)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-4)/2.0);  //Zupq1=842
            pa3 = (a1/((z-4)*(z-5)))*pow(arga,-6)*(1-cos((z-5)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-5)/2.0));  //Zupq1=843
            pa = 9*(pa1-2*pa2+pa3)*sin(a);  //Zupq1=844
            argb = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(b)*cos(b))/(z+1);  //Zupq1=845
            pb1 = a1*pow(argb,-4)*(1+cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0));  //Zupq1=846
            pb2 = (a1/(z-4))*pow(argb,-5)*sin((z-4)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-4)/2.0);  //Zupq1=847
            pb3 = (a1/((z-4)*(z-5)))*pow(argb,-6)*(1-cos((z-5)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-5)/2.0));  //Zupq1=848
            pb = 9*(pb1-2*pb2+pb3)*sin(b);  //Zupq1=849
        }/*3*/  /*  of biaxial ellipsoid  */  //Zupq1=850

        if ( i1==14 )
        {/*3*/   /*  triaxial ellipsoid, isotropic  */  //Zupq1=852
            ella = r;  //Zupq1=853
            ellb = l;  //Zupq1=854
            ellc = r/p1;  //Zupq1=855
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=856
            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=857
            arga = qn*sqrt(pow(ella*cos(a)*sin(delta),2)+pow(ellb*sin(a)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);  //Zupq1=858
            a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=859
            pa1 = a1*pow(arga,-4)*(1+cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0));  //Zupq1=860
            pa2 = (a1/(z-4))*pow(arga,-5)*sin((z-4)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-4)/2.0);  //Zupq1=861
            pa3 = (a1/((z-4)*(z-5)))*pow(arga,-6)*(1-cos((z-5)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-5)/2.0));  //Zupq1=862
            pa = 9*(pa1-2*pa2+pa3);  //Zupq1=863
            /* pa:=power(3*(sin(qn*r)-qn*r*cos(qn*r))/(qn*qn*qn*r*r*r+eps),2);   }  //Zupq1=864 */
            /* pa:=1.05;  //Zupq1=865 */
            argb = qn*sqrt(pow(ella*cos(b)*sin(delta),2)+pow(ellb*sin(b)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);  //Zupq1=866
            pb1 = a1*pow(argb,-4)*(1+cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0));  //Zupq1=867
            pb2 = (a1/(z-4))*pow(argb,-5)*sin((z-4)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-4)/2.0);  //Zupq1=868
            pb3 = (a1/((z-4)*(z-5)))*pow(argb,-6)*(1-cos((z-5)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-5)/2.0));  //Zupq1=869
            pb = 9*(pb1-2*pb2+pb3);  //Zupq1=870
            /* pb:=power(3*(sin(qn*r*1.04)-qn*r*cos(qn*r))/(qn*qn*qn*r*r*r+eps),2);   }  //Zupq1=871 */
            /* pb:=1.04;  //Zupq1=872 */
        }/*3*/  /*  of ellipsoid  */  //Zupq1=873

        if ( i1==15 )
        {/*3*/   /*  barrel area integration  */  //Zupq1=875
            pa1 = r*pow(1.0-pow(a/l,p1),1/p1);  //Zupq1=876
            arga = -r*pow(a/l,p1-1)*pow(1.0-pow(a/l,p1),(1/p1)-1)/l;  //Zupq1=877
            pa = pa1*sqrt(1.0+arga*arga);  //Zupq1=878
            pb1 = r*pow(1.0-pow(b/l,p1),1/p1);  //Zupq1=879
            argb = -r*pow(b/l,p1-1)*pow(1.0-pow(b/l,p1),(1/p1)-1)/l;  //Zupq1=880
            pb = pb1*sqrt(1.0+argb*argb);  //Zupq1=881
        }/*3*/  //Zupq1=882

        if ( i1==16 )
        {/*3*/   /*  barrel, x-axis  */  //Zupq1=884
            z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=885
            qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=886
            arga = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(a)*sin(delta)+eps));  //Zupq1=887
            argb = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(b)*sin(delta)+eps));  //Zupq1=888
            if ( i3==0 )
            {/*4*/ /*  P(q)  */  //Zupq1=889
                a1 = 9*pow(z+1,4)/(2.0*z*(z-1)*(z-2)*(z-3));  //Zupq1=890
                pa = (a1/(arga*arga));  //Zupq1=891
                pb = (a1/(argb*argb));  //Zupq1=892
            }/*4*/  //Zupq1=893
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=894
                pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=895
                pa = pa1*pa1;  //Zupq1=896
                pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=897
                pb = pb1*pb1;  //Zupq1=898
            }/*4*/  //Zupq1=899
        }/*3*/  //Zupq1=900

        if ( i1==17 )
        {/*3*/   /*  superball integration  */  //Zupq1=902
            aa = r;  //Zupq1=903
            bb = p1;  //Zupq1=904
            cc = l;  //Zupq1=905
            /* arga1:=-(l/r)*power(delta/r,p1-1)*power(1-power(delta/r,p1)-power(a/r,p1),(1/p1)-1);  //Zupq1=906 */
            /* arga2:=-(l/r)*power(a/r,p1-1)*power(1-power(delta/r,p1)-power(a/r,p1),(1/p1)-1);  //Zupq1=907 */
            /* pa:=sqrt(1+arga1*arga1+arga2*arga2);  //Zupq1=908 */
            /* pa:=power(r,4)*power(l,2)*(power(power(r*r*cos(delta),p1)+(power(cos(a),p1)+power(sin(a),p1))*power(r*l*sin(delta),p1),-(2+p1)/p1))*  //Zupq1=909 */
            /*       sqrt(power(r,4*p1)*power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(a),2*p1-2)+power(sin(a),2*p1-2))*power(r*l*sin(delta),2*p1));  //Zupq1=910 */
            pa = aa*aa*bb*bb*cc*cc*(pow(pow(aa*bb*cos(delta),alfa)+(pow(bb*cc*cos(a),alfa)+pow(aa*cc*sin(a),alfa))*pow(sin(delta),alfa),-(2+alfa)/alfa))*  //Zupq1=911
                 sqrt(pow(aa*bb,2*alfa)*pow(cos(delta),2*alfa-2)*pow(sin(delta),2)+(pow(bb*cc,2*alfa)*pow(cos(a),2*alfa-2)+pow(aa*cc,2*alfa)*pow(sin(a),2*alfa-2))*pow(sin(delta),2*alfa));  //Zupq1=912
            /* pa:=r*r*(power(power(cos(delta),p1)+(power(cos(a),p1)+power(sin(a),p1))*power(sin(delta),p1),-(2+p1)/p1))*  //Zupq1=913 */
            /*       sqrt(power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(a),2*p1-2)+power(sin(a),2*p1-2))*power(sin(delta),2*p1));  //Zupq1=914 */
            /* argb1:=-(l/r)*power(delta/r,p1-1)*power(1-power(delta/r,p1)-power(b/r,p1),(1/p1)-1);  //Zupq1=915 */
            /* argb2:=-(l/r)*power(b/r,p1-1)*power(1-power(delta/r,p1)-power(b/r,p1),(1/p1)-1);  //Zupq1=916 */
            /* pb:=sqrt(1+argb1*argb1+argb2*argb2);  //Zupq1=917 */
            /* pb:=power(r,4)*power(l,2)*(power(power(r*r*cos(delta),p1)+(power(cos(b),p1)+power(sin(b),p1))*power(r*l*sin(delta),p1),-(2+p1)/p1))*  //Zupq1=918 */
            /*       sqrt(power(r,4*p1)*power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(b),2*p1-2)+power(sin(b),2*p1-2))*power(r*l*sin(delta),2*p1));  //Zupq1=919 */
            pb = aa*aa*bb*bb*cc*cc*(pow(pow(aa*bb*cos(delta),alfa)+(pow(bb*cc*cos(b),alfa)+pow(aa*cc*sin(b),alfa))*pow(sin(delta),alfa),-(2+alfa)/alfa))*  //Zupq1=920
                 sqrt(pow(aa*bb,2*alfa)*pow(cos(delta),2*alfa-2)*pow(sin(delta),2)+(pow(bb*cc,2*alfa)*pow(cos(b),2*alfa-2)+pow(aa*cc,2*alfa)*pow(sin(b),2*alfa-2))*pow(sin(delta),2*alfa));  //Zupq1=921
            /* pb:=r*r*(power(power(cos(delta),p1)+(power(cos(b),p1)+power(sin(b),p1))*power(sin(delta),p1),-(2+p1)/p1))*  //Zupq1=922 */
            /*       sqrt(power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(b),2*p1-2)+power(sin(b),2*p1-2))*power(sin(delta),2*p1));  //Zupq1=923 */
        }/*3*/  //Zupq1=924

        pq = 0.5*(b-a)*(pa+pb);  //Zupq1=926
        trapzdchid_cnt = 1;  //Zupq1=927
    }/*2*/  //Zupq1=928
    else
    {/*2*/  //Zupq1=929
        tnm = trapzdchid_cnt;  //Zupq1=930
        del = (b-a)/tnm;  //Zupq1=931
        x = a+0.5*del;  //Zupq1=932
        sump = 0.0;  //Zupq1=933
        for ( j=1; j<=trapzdchid_cnt; j++ )
        {/*3*/  //Zupq1=934
            if ( i1==1 )
            {/*4*/  /*  cylinder, general case  */  //Zupq1=935
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=936
                mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);  //Zupq1=937
                mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);  //Zupq1=938
                /* mlx3:=p31*cos(x)*sin(delta)+p33*cos(delta);  //Zupq1=939 */
                argx = (qx*mlx1+qy*mlx2+eps)*l/(z+1);  //Zupq1=940
                if ( i3==0 )
                {/*5*/ /*  P(q)  */  //Zupq1=941
                    a1 = 1/(2.0*z*(z-1));  //Zupq1=942
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=943
                }/*5*/  //Zupq1=944
                if ( i3==1 )
                {/*5*/ /*  F(q)  */  //Zupq1=945
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=946
                    px = px1*px1;  //Zupq1=947
                }/*5*/  //Zupq1=948
            }/*4*/  //Zupq1=949
            if ( i1==2 )
            {/*4*/   /*  cylinder, x-axis  */  //Zupq1=950
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=951
                argx = (qx*cos(delta)-qy*sin(x)*sin(delta)+eps)*l/(z+1);  //Zupq1=952
                if ( i3==0 )
                {/*5*/ /*  P(q)  */  //Zupq1=953
                    a1 = 1/(2.0*z*(z-1));  //Zupq1=954
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=955
                }/*5*/  //Zupq1=956
                if ( i3==1 )
                {/*5*/ /*  F(q)  */  //Zupq1=957
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=958
                    px = px1*px1;  //Zupq1=959
                }/*5*/  //Zupq1=960
            }/*4*/  //Zupq1=961
            if ( i1==3 )
            {/*4*/   /*  cylinder, y-axis  */  //Zupq1=962
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=963
                a1 = 1/(2.0*z*(z-1));  //Zupq1=964
                argx = (qx*sin(x)*sin(delta)+qy*cos(delta)+eps)*l/(z+1);  //Zupq1=965
                if ( i3==0 )
                {/*5*/ /*  P(q)  */  //Zupq1=966
                    a1 = 1/(2.0*z*(z-1));  //Zupq1=967
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=968
                }/*5*/  //Zupq1=969
                if ( i3==1 )
                {/*5*/ /*  F(q)  */  //Zupq1=970
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=971
                    px = px1*px1;  //Zupq1=972
                }/*5*/  //Zupq1=973
            }/*4*/  //Zupq1=974
            if ( i1==4 )
            {/*4*/   /*  cylinder, -z-axis  */  //Zupq1=975
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=976
                a1 = 1/(2.0*z*(z-1));  //Zupq1=977
                argx = (qx*sin(x)*sin(delta)-qy*cos(x)*sin(delta)+eps)*l/(z+1);  //Zupq1=978
                if ( i3==0 )
                {/*5*/ /*  P(q)  */  //Zupq1=979
                    a1 = 1/(2.0*z*(z-1));  //Zupq1=980
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=981
                }/*5*/  //Zupq1=982
                if ( i3==1 )
                {/*5*/ /*  F(q)  */  //Zupq1=983
                    px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=984
                    px = px1*px1;  //Zupq1=985
                }/*5*/  //Zupq1=986
            }/*4*/  //Zupq1=987
            if ( i1==5 )
            {/*4*/   /*  general series expansion  */  //Zupq1=988
                mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);  //Zupq1=989
                mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);  //Zupq1=990
                //mlx3 = p31*cos(x)*sin(delta)+p33*cos(delta);  //Zupq1=991
                px = pow(mlx1,i3)*pow(mlx2,i4);  //Zupq1=992
            }/*4*/  //Zupq1=993
            if ( i1==6 )
            {/*4*/    /*  unit cell rotation  */  //Zupq1=994
                dqxx = qx-qhkl*(p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta));  //Zupq1=995
                dqyx = qy-qhkl*(p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta));  //Zupq1=996
                dqzx = qz-qhkl*(p31*cos(x)*sin(delta)+p33*cos(delta));  //Zupq1=997
                dqs1x = (dqxx*ax1x+dqyx*ax1y+dqzx*ax1z)/(ax1n*sigx);  //Zupq1=998
                dqs2x = (dqxx*ax2x+dqyx*ax2y+dqzx*ax2z)/(ax2n*sigy);  //Zupq1=999
                dqs3x = (dqxx*ax3x+dqyx*ax3y+dqzx*ax3z)/(ax3n*sigz);  //Zupq1=1000
                px = exp(-4*(dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x)/M_PI);  //Zupq1=1001
            }/*4*/  //Zupq1=1002
            if ( i1==7 )
            {/*4*/    /*  fiber unit cell rotation  */  //Zupq1=1003
                l1 = sin(theta)*cos(phi);  //Zupq1=1004
                l2 = sin(theta)*sin(phi);  //Zupq1=1005
                l3 = -cos(theta);  //Zupq1=1006

                r11x = cos(x)+(1-cos(x))*l1*l1;  //Zupq1=1008
                r12x = -l3*sin(x)+(1-cos(x))*l1*l2;  //Zupq1=1009
                r13x = -l2*sin(x)+(1-cos(x))*l1*l3;  //Zupq1=1010
                r21x = l3*sin(x)+(1-cos(x))*l1*l2;  //Zupq1=1011
                r22x = cos(x)+(1-cos(x))*l2*l2;  //Zupq1=1012
                r23x = l1*sin(x)+(1-cos(x))*l2*l3;  //Zupq1=1013
                r31x = l2*sin(x)+(1-cos(x))*l1*l3;  //Zupq1=1014
                r32x = -l1*sin(x)+(1-cos(x))*l2*l3;  //Zupq1=1015
                r33x = cos(x)+(1-cos(x))*l3*l3;  //Zupq1=1016

                /*  rotate this scattering vector  */  //Zupq1=1023
                qxhklx = r11x*qxn+r12x*qyn+r13x*qzn;  //Zupq1=1024
                qyhklx = r21x*qxn+r22x*qyn+r23x*qzn;  //Zupq1=1025
                qzhklx = r31x*qxn+r32x*qyn+r33x*qzn;  //Zupq1=1026

                dqxx = qx-qxhklx;  //Zupq1=1093
                dqyx = qy-qyhklx;  //Zupq1=1094
                dqzx = qz-qzhklx;  //Zupq1=1095

                ax1xx = (r11x*ax1x+r12x*ax1y+r13x*ax1z);  //Zupq1=1097
                ax1yx = (r21x*ax1x+r22x*ax1y+r23x*ax1z);  //Zupq1=1098
                ax1zx = (r31x*ax1x+r32x*ax1y+r33x*ax1z);  //Zupq1=1099
                ax2xx = (r11x*ax2x+r12x*ax2y+r13x*ax2z);  //Zupq1=1100
                ax2yx = (r21x*ax2x+r22x*ax2y+r23x*ax2z);  //Zupq1=1101
                ax2zx = (r31x*ax2x+r32x*ax2y+r33x*ax2z);  //Zupq1=1102
                ax3xx = (r11x*ax3x+r12x*ax3y+r13x*ax3z);  //Zupq1=1103
                ax3yx = (r21x*ax3x+r22x*ax3y+r23x*ax3z);  //Zupq1=1104
                ax3zx = (r31x*ax3x+r32x*ax3y+r33x*ax3z);  //Zupq1=1105

                dqs1x = (dqxx*ax1xx+dqyx*ax1yx+dqzx*ax1zx)/(ax1n*sigx);  //Zupq1=1107
                dqs2x = (dqxx*ax2xx+dqyx*ax2yx+dqzx*ax2zx)/(ax2n*sigy);  //Zupq1=1108
                dqs3x = (dqxx*ax3xx+dqyx*ax3yx+dqzx*ax3zx)/(ax3n*sigz);  //Zupq1=1109

                argx = dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x;  //Zupq1=1115
                px = exp(-4*argx/M_PI);  //Zupq1=1116
                /* if (argx>2) then px:=eps  //Zupq1=1117 */
                /*    else px:=exp(-4*argx/pi);  //Zupq1=1118 */
            }/*4*/  //Zupq1=1119
            if ( i1==8 )
            {/*4*/  /*  disk, general case  */  //Zupq1=1120
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1121
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1122
                qxl = qx/qn;  //Zupq1=1123
                qyl = qy/qn;  //Zupq1=1124
                mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);  //Zupq1=1125
                mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);  //Zupq1=1126
                /* mlx3:=p31*cos(x)*sin(delta)+p33*cos(delta);  //Zupq1=1127 */
                qnnx = qxl*mlx1+qyl*mlx2;  //Zupq1=1128
                argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1129

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1131
                    if ( argx<0.015 )
                    {/*6*/  //Zupq1=1132
                        px = 1;  //Zupq1=1133
                        oldpx = 0;  //Zupq1=1134
                        argser = 1;  //Zupq1=1135
                        for ( i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1136
                            argser = argser*argx*argx/4.0;  //Zupq1=1137
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1138 */
                            px = px+carr1[i]*argser;  //Zupq1=1139
                            delser = fabs((px-oldpx)/px);  //Zupq1=1140
                            if ( delser<0.0001 ) break; /* goto 14; */  //Zupq1=1141
                            oldpx = px;  //Zupq1=1142
                        }/*7*/  //Zupq1=1143
                        /*14:*/  //Zupq1=1144
                        if ( i3==1 ) px = px*px;  //Zupq1=1145
                    }/*6*/  //Zupq1=1146
                    else
                    {/*6*/  //Zupq1=1147
                        if ( i3==0 )
                        {/*7*/  /*  P(q)  */  //Zupq1=1148
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1149
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1150
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1151
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1152
                        }/*7*/  //Zupq1=1153
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1154
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1155
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1156
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1157
                            px = px3*px3;  //Zupq1=1158
                        }/*7*/  //Zupq1=1159
                    }/*6*/  //Zupq1=1160
                }/*5*/  //Zupq1=1161
                else
                {/*5*/  //Zupq1=1162
                    px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1163
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1164
                }/*5*/  //Zupq1=1165
            }/*4*/  //Zupq1=1166

            if ( i1==9 )
            {/*4*/   /*  disk, x-axis  */  //Zupq1=1168
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1169
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1170
                qxl = qx/qn;  //Zupq1=1171
                qyl = qy/qn;  //Zupq1=1172
                qnnx = qxl*cos(delta)-qyl*sin(x)*sin(delta);  //Zupq1=1173
                argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1174

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1176
                    if ( argx<0.015 )
                    {/*6*/  /*  series expansion  */  //Zupq1=1177
                        px = 1;  //Zupq1=1178
                        oldpx = 0;  //Zupq1=1179
                        argser = 1;  //Zupq1=1180
                        for ( i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1181
                            argser = argser*argx*argx/4.0;  //Zupq1=1182
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1183 */
                            px = px+carr1[i]*argser;  //Zupq1=1184
                            delser = fabs((px-oldpx)/px);  //Zupq1=1185
                            if ( delser<0.0001 ) break; /* goto 17; */  //Zupq1=1186
                            oldpx = px;  //Zupq1=1187
                        }/*7*/  //Zupq1=1188
                        /*17:*/  //Zupq1=1189
                        if ( i3==1 ) px = px*px;  //Zupq1=1190
                    }/*6*/  //Zupq1=1191
                    else
                    {/*6*/  //Zupq1=1192
                        if ( i3==0 )
                        {/*7*/  /*  P(q)  */  //Zupq1=1193
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1194
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1195
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1196
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1197
                        }/*7*/  //Zupq1=1198
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1199
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1200
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1201
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1202
                            px = px3*px3;  //Zupq1=1203
                        }/*7*/  //Zupq1=1204
                    }/*6*/  //Zupq1=1205
                }/*5*/  //Zupq1=1206
                else
                {/*5*/  //Zupq1=1207
                    px = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1208
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1209
                }/*5*/  //Zupq1=1210
            }/*4*/  //Zupq1=1211

            if ( i1==10 )
            {/*4*/   /*  disk, y-axis  */  //Zupq1=1213
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1214
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1215
                qxl = qx/qn;  //Zupq1=1216
                qyl = qy/qn;  //Zupq1=1217
                qnnx = qxl*sin(x)*sin(delta)+qyl*cos(delta);  //Zupq1=1218
                argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1219

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1221
                    if ( argx<0.015 )
                    {/*6*/  //Zupq1=1222
                        px = 1;  //Zupq1=1223
                        oldpx = 0;  //Zupq1=1224
                        argser = 1;  //Zupq1=1225
                        for ( i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1226
                            argser = argser*argx*argx/4.0;  //Zupq1=1227
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1228 */
                            px = px+carr1[i]*argser;  //Zupq1=1229
                            delser = fabs((px-oldpx)/px);  //Zupq1=1230
                            if ( delser<0.0001 ) break; /* goto 20; */  //Zupq1=1231
                            oldpx = px;  //Zupq1=1232
                        }/*7*/  //Zupq1=1233
                        /*20:*/  //Zupq1=1234
                        if ( i3==1 ) px = px*px;  //Zupq1=1235
                    }/*6*/  //Zupq1=1236
                    else
                    {/*6*/  //Zupq1=1237
                        if ( i3==0 )
                        {/*7*/  /*  P(q)  */  //Zupq1=1238
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1239
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1240
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1241
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1242
                        }/*7*/  //Zupq1=1243
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1244
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1245
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1246
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1247
                            px = px3*px3;  //Zupq1=1248
                        }/*7*/  //Zupq1=1249
                    }/*6*/  //Zupq1=1250
                }/*5*/  //Zupq1=1251
                else
                {/*5*/  //Zupq1=1252
                    px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1253
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1254
                }/*5*/  //Zupq1=1255
            }/*4*/  //Zupq1=1256

            if ( i1==11 )
            {/*4*/   /*  disk, z-axis  */  //Zupq1=1258
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1259
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1260
                qxl = qx/qn;  //Zupq1=1261
                qyl = qy/qn;  //Zupq1=1262
                qnnx = qxl*sin(x)*sin(delta)-qyl*cos(x)*sin(delta);  //Zupq1=1263
                argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1264

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1266
                    if ( argx<0.015 )
                    {/*6*/  //Zupq1=1267
                        px = 1;  //Zupq1=1268
                        oldpx = 0;  //Zupq1=1269
                        argser = 1;  //Zupq1=1270
                        for ( i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1271
                            argser = argser*argx*argx/4.0;  //Zupq1=1272
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1273 */
                            px = px+carr1[i]*argser;  //Zupq1=1274
                            delser = fabs((px-oldpx)/px);  //Zupq1=1275
                            if ( delser<0.0001 ) break; /* goto 23; */  //Zupq1=1276
                            oldpx = px;  //Zupq1=1277
                        }/*7*/  //Zupq1=1278
                        /*23:*/  //Zupq1=1279
                        if ( i3==1 ) px = px*px;  //Zupq1=1280
                    }/*6*/  //Zupq1=1281
                    else
                    {/*6*/  //Zupq1=1282
                        if ( i3==0 )
                        {/*7*/  /*  P(q)  */  //Zupq1=1283
                            px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1284
                            px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1285
                            px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1286
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1287
                        }/*7*/  //Zupq1=1288
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1289
                            px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1290
                            px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1291
                            px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1292
                            px = px3*px3;  //Zupq1=1293
                        }/*7*/  //Zupq1=1294
                    }/*6*/  //Zupq1=1295
                }/*5*/  //Zupq1=1296
                else
                {/*5*/  //Zupq1=1297
                    px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1298
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1299
                }/*5*/  //Zupq1=1300
            }/*4*/  //Zupq1=1301

            if ( i1==12 )
            {/*4*/   /*  isotropic cube  */  //Zupq1=1303
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1322
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1323
                argx1 = qn*sin(delta)*cos(x)*r/(z+1);  //Zupq1=1324
                argx2 = qn*sin(delta)*sin(x)*r/(z+1);  //Zupq1=1325
                argx3 = qn*cos(delta)*l/(z+1);  //Zupq1=1326
                if ( i3==0 )
                {/*5*/   /*  P(q)  */  //Zupq1=1327
                    a1 = 1/(2.0*z*(z-1));  //Zupq1=1328
                    if ( argx1 < 0.001 ) px1 = 1; else  //Zupq1=1329
                            px1 = (a1/(argx1*argx1+eps))*(1-cos((z-1)*atan(2.0*argx1))/pow(1.0+4*argx1*argx1,(z-1)/2.0));  //Zupq1=1330
                    if ( argx2 < 0.001 ) px2 = 1; else  //Zupq1=1331
                            px2 = (a1/(argx2*argx2+eps))*(1-cos((z-1)*atan(2.0*argx2))/pow(1.0+4*argx2*argx2,(z-1)/2.0));  //Zupq1=1332
                    if ( argx3 < 0.001 ) px3 = 1; else  //Zupq1=1333
                            px3 = (a1/(argx3*argx3+eps))*(1-cos((z-1)*atan(2.0*argx3))/pow(1.0+4*argx3*argx3,(z-1)/2.0));  //Zupq1=1334
                    px = px1*px2*px3;  //Zupq1=1335
                }/*5*/  //Zupq1=1336
                if ( i3==1 )
                {/*5*/ /*  F(q)  */  //Zupq1=1337
                    a1 = 1/z;  //Zupq1=1338
                    if ( argx1 < 0.001 ) px1 = 1; else  //Zupq1=1339
                            px1 = (a1/(argx1+eps))*sin(z*atan(argx1))/pow(1.0+argx1*argx1,z/2.0);  //Zupq1=1340
                    if ( argx2 < 0.001 ) px2 = 1; else  //Zupq1=1341
                            px2 = (a1/(argx2+eps))*sin(z*atan(argx2))/pow(1.0+argx2*argx2,z/2.0);  //Zupq1=1342
                    if ( argx3 < 0.001 ) px3 = 1; else  //Zupq1=1343 TODO: es stand arga3 hier...
                            px3 = (a1/(argx3+eps))*sin(z*atan(argx3))/pow(1.0+argx3*argx3,z/2.0);  //Zupq1=1344
                    px = px1*px1*px2*px2*px3*px3;  //Zupq1=1345
                }/*5*/  //Zupq1=1346
            }/*4*/  //Zupq1=1347

            if ( i1==13 )
            {/*4*/   /*  biaxial ellipsoid, isotropic  */  //Zupq1=1349
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1350
                epsi = l/r;  //Zupq1=1351
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1352
                argx = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(x)*cos(x))/(z+1);  //Zupq1=1353
                a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=1354
                px1 = a1*pow(argx,-4)*(1+cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0));  //Zupq1=1355
                px2 = (a1/(z-4))*pow(argx,-5)*sin((z-4)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-4)/2.0);  //Zupq1=1356
                px3 = (a1/((z-4)*(z-5)))*pow(argx,-6)*(1-cos((z-5)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-5)/2.0));  //Zupq1=1357
                px = 9*(px1-2*px2+px3)*sin(x);  //Zupq1=1358
            }/*4*/  /*  of biaxial ellipsoid  */  //Zupq1=1359

            if ( i1==14 )
            {/*4*/   /*  triaxial ellipsoid, isotropic  */  //Zupq1=1361
                ella = r;  //Zupq1=1362
                ellb = l;  //Zupq1=1363
                ellc = r/p1;  //Zupq1=1364
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1365
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1366
                argx = qn*sqrt(pow(ella*cos(x)*sin(delta),2)+pow(ellb*sin(x)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);  //Zupq1=1367
                a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=1368
                px1 = a1*pow(argx,-4)*(1+cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0));  //Zupq1=1369
                px2 = (a1/(z-4))*pow(argx,-5)*sin((z-4)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-4)/2.0);  //Zupq1=1370
                px3 = (a1/((z-4)*(z-5)))*pow(argx,-6)*(1-cos((z-5)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-5)/2.0));  //Zupq1=1371
                px = 9*(px1-2*px2+px3);  //Zupq1=1372
            }/*4*/  /*  of triaxial ellipsoid  */  //Zupq1=1373

            if ( i1==15 )
            {/*4*/   /*  barrel area integration  */  //Zupq1=1375
                px1 = r*pow(1.0-pow(x/l,p1),1/p1);  //Zupq1=1376
                argx = -r*pow(x/l,p1-1)*pow(1.0-pow(x/l,p1),(1/p1)-1)/l;  //Zupq1=1377
                px = px1*sqrt(1.0+argx*argx);  //Zupq1=1378
            }/*4*/  //Zupq1=1379

            if ( i1==16 )
            {/*4*/   /*  barrel, x-axis  */  //Zupq1=1381
                z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1382
                qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1383
                argx = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(x)*sin(delta)+eps));  //Zupq1=1384
                if ( i3==0 )
                {/*5*/ /*  P(q)  */  //Zupq1=1385
                    a1 = 9*pow(z+1,4)/(2.0*z*(z-1)*(z-2)*(z-3));  //Zupq1=1386
                    px = (a1/(argx*argx));  //Zupq1=1387
                }/*5*/  //Zupq1=1388
                if ( i3==1 )
                {/*5*/ /*  F(q)  */  //Zupq1=1389
                    pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=1390
                    pa = pa1*pa1;  //Zupq1=1391
                    pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=1392
                    pb = pb1*pb1;  //Zupq1=1393
                }/*5*/  //Zupq1=1394
            }/*4*/  //Zupq1=1395

            if ( i1==17 )
            {/*4*/   /*  superball integration  */  //Zupq1=1397
                aa = r;     // von oben kopiert... (TODO?)
                bb = p1;
                cc = l;
                /*  argx1:=-(l/r)*power(delta/r,p1-1)*power(1-power(delta/r,p1)-power(x/r,p1),(1/p1)-1);  //Zupq1=1398 */
                /*  argx2:=-(l/r)*power(x/r,p1-1)*power(1-power(delta/r,p1)-power(x/r,p1),(1/p1)-1);  //Zupq1=1399 */
                /*  px:=sqrt(1+argx1*argx1+argx2*argx2);  //Zupq1=1400 */
                /* px:=power(r,4)*power(l,2)*(power(power(r*r*cos(delta),p1)+(power(cos(x),p1)+power(sin(x),p1))*power(r*l*sin(delta),p1),-(2+p1)/p1))*  //Zupq1=1401 */
                /*    sqrt(power(r,4*p1)*power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(x),2*p1-2)+power(sin(x),2*p1-2))*power(r*l*sin(delta),2*p1));  //Zupq1=1402 */
                px = aa*aa*bb*bb*cc*cc*(pow(pow(aa*bb*cos(delta),alfa)+(pow(bb*cc*cos(x),alfa)+pow(aa*cc*sin(x),alfa))*pow(sin(delta),alfa),-(2+alfa)/alfa))*  //Zupq1=1403
                     sqrt(pow(aa*bb,2*alfa)*pow(cos(delta),2*alfa-2)*pow(sin(delta),2)+(pow(bb*cc,2*alfa)*pow(cos(x),2*alfa-2)+pow(aa*cc,2*alfa)*pow(sin(x),2*alfa-2))*pow(sin(delta),2*alfa));  //Zupq1=1404
                /* px:=r*r*(power(power(cos(delta),p1)+(power(cos(x),p1)+power(sin(x),p1))*power(sin(delta),p1),-(2+p1)/p1))*  //Zupq1=1405 */
                /*    sqrt(power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(x),2*p1-2)+power(sin(x),2*p1-2))*power(sin(delta),2*p1));  //Zupq1=1406 */
            }/*4*/  //Zupq1=1407

            sump = sump+px;  //Zupq1=1409
            x = x+del;  //Zupq1=1410
        }/*3*/  //Zupq1=1411
        pq = 0.5*(pq+(b-a)*sump/tnm);  //Zupq1=1412
        trapzdchid_cnt = 2*trapzdchid_cnt;  //Zupq1=1413
    }/*2*/  //Zupq1=1414
}/*1*/  //Zupq1=1415






/*#ifndef __CUDACC__
    if ( dbgFlag() )
    {   // Kommt nur einmal pro Thread zu Beginn der Berechnungen
        qDebug() << "formpq:"  << "part"<<part << "ordis"<<ordis << "orcase"<<orcase << "cs"<<cs
                 << "q"<<q << "limq1"<<limq1 << "limq4"<<limq4 << "limql"<<limql;
    }
#endif*/




#ifdef __CUDACC__
__host__ __device__
#endif
//    double SasCalc_GENERIC_calculation::formpq(double length, double radius, double sigmal, double sigmar, double p1,
//                                        double rho, double alfa, double theta, double phi, double limql, double limq1,
//                                        double limq2, double /*limq3*/, double limq4, double limq5, double limq6,
//                                        double /*limq7*/, double /*limq8*/, double /*limq9*/, double qx, double qy, double qxs,
//                                        double qys, double q, double norm, double por,
//                                        int part, int cs, int ordis, int orcase,
//                                        const double *myarray, // CoeffArrayType
//                                        double *carr1p, double *carr2p, double *carr3p, // CoeffArrayType
//                                        double *carr4p, double *carr5p, double *carr6p, // CoeffArrayType
//                                        double *carr7p, double *carr8p, double *carr9p // CoeffArrayType   /*Z0311=14582*/
//                                        /*ArrayImax2D carr11pm, ArrayImax2D carr22pm*/) const   /*Z=14910*/
double SasCalc_GENERIC_calculation::formpq(double sigmal, double limql, double qx, double qy, double qxs, double qys, double q, int ordis) const   /*Z=14910*/
{/*1*/  //Z=15188
    // double sigmal ist nicht zu ersetzen, da es an einer Stelle CALC.epsilon ist, sonst nur CALC.params.sigmal
    // int ordis ist nicht zu ersetzen, da es einmal eine feste Zahl ist, sonst nur CALC.ordis

    int    ii, jj, n, nser, m, mser, lser, /*indx,*/ inmax;  //Z=15194
    double pqsum, oldpqsum, binsum, delser, argq, arglq, /*argp1,*/ pqr, pql, pq1, pq2, pq3, epsi, qq2;  //Z=15195
    double ccc1, ccc2, ccc3, vv3, zl, zr, radiusm, argqx, argqy, pqrx, pqry; //, ella, ellb, ellc;  //Z=15196
    double cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8, cc9, cc10;  //Z=15197
    double ac1, ac2, ac3, ac4, ac5, ac6, ac7, ac8, ac9, ac10;  //Z=15198
    double argbm, nenbm, argbp, nenbp, argep, nenep, argem, nenem, arggp, nengp;  //Z=15199
    double arggm, nengm, /*argim, nenim, argip, nenip,*/ F121, F122, F123; //, F124, F125, F126;  //Z=15200
    double qqn[200], /*qqnx[200], qqny[200],*/ fkv[200]; //, gam3[200];  //Z=15201
    double F12, F12ser, F12asz, F12sum, F12sez, oldF12sez, F12asy, F12as1z, F12as2z, F12as3z, F12as4z;  //Z=15202
    double v, e0, e1, pz2v, pz2v1, pz2v2, lim, lim1, xrz, arg, nen, arg1, nen1, arg2, nen2;  //Z=15203
    double a1m, a2m, xijm, arglmz, nenlmz, xijp, arglpz, nenlpz, vvm, rmax, del, delc;  //Z=15204
    //double nom, lout, lin, lliph, llipt, phiax, phiin, phiout, philiph, philipt;  //Z=15205
    double dim, xrad, xradp, lim2, lim3, lim4, lim5, lim6;  //Z=15206
    double a1, b1, b2, b1s, d0, /*d1,*/ ee0, ee1, x1z, x12z, x2z, x22z, gb1s;  //Z=15207
    double gz1, preg1, preg3, preg4, pzvc, pzvc1, pzvc2, pzac, pzac1, /*pzac2,*/ pzc, pzc1, pza, pzva, pzva1; //, dnv0,  pvav0, pvav10, pva0;  //Z=15208
    double F22sez, oldF22sez, F22, F32sez, oldF32sez, F32;  //Z=15209
    double F42sez, oldF42sez, F42, F52sez, oldF52sez, F52, F62sez, oldF62sez, F62;  //Z=15210
    double arg11, nen11, arg12, nen12, arg13, nen13;  //Z=15211
    double /*arg21, nen21, arg22, nen22,*/ arg210, nen210, arg220, nen220, arg23, nen23, arg24, nen24, arg25, nen25, arg26, nen26, arg27, nen27, arg28, nen28;  //Z=15212
    double /*F22as1sum1z, F22as1sum2z,*/ F22as10z, /*F22as1z,*/ F22as1sum1z0, F22as1sum2z0, F22as1z0;  //Z=15213
    double a22as21z, F22as21z, a22as22z, F22as22z, a22as23z, F22as23z, a22as24z, F22as24z, F22as20z, F22as2z, /*F22asz,*/ F22asz0;  //Z=15214
    double /*arg31, nen31, arg32, nen32,*/ arg310, nen310, arg320, nen320, arg33, nen33, arg34, nen34, arg35, nen35;  //Z=15215
    double /*F32as1sum1z, F32as1sum2z,*/ F32as10z, /*F32as1z,*/ F32as1sum1z0, F32as1sum2z0, F32as1z0;  //Z=15216
    double F32as21z, F32as22z, F32as23z, F32as24z, F32as20z, F32as2z, /*F32asz,*/ F32asz0;  //Z=15217
    double arg41, nen41, arg42, nen42, /*arg43, nen43,*/ arg44, nen44, arg45, nen45;  //Z=15218
    double F42as10z, /*F42as1sumz, F42as1z,*/ F42as1z0, F42as20z, F42as21, F42as22, /*F42as23, F42as2z,*/ F42as2z0;  //Z=15219
    double F42as30z, F42as24, /*F42as25,*/ F42as26, /*F42as3z,*/ F42as3z0, F42as40z, F42as27, F42as28, F42as29, F42as4z, /*F42asz,*/ F42asz0;  //Z=15220
    double arg51, nen51, arg52, nen52, /*arg53, nen53,*/ arg54, nen54, /*arg55, nen55,*/ arg56, nen56;  //Z=15221
    double arg57, nen57, arg58, nen58, arg59, nen59, arg510, nen510;  //Z=15222
    double F52as10z, /*F52as1sumz, F52as1z,*/ F52as1z0, F52as20z, F52as21, F52as22, /*F52as23, F52as2z,*/ F52as2z0;  //Z=15223
    double F52as30z, F52as24, /*F52as25,*/ F52as26, /*F52as3z,*/ F52as3z0, F52as40z, F52as27, F52as28, F52as29, F52as4z, /*F52asz,*/ F52asz0;  //Z=15224
    double arg61, nen61, arg62, nen62, /*arg63, nen63,*/ arg64, nen64, arg65, nen65;  //Z=15225
    double F62as10z, /*F62as1sumz, F62as1z,*/ F62as1z0, F62as20z, F62as21, F62as22, /*F62as23, F62as2z,*/ F62as2z0;  //Z=15226
    double F62as30z, F62as24, /*F62as25,*/ F62as26, /*F62as3z,*/ F62as3z0, F62as40z, F62as27, F62as28, F62as29, F62as4z, /*F62asz,*/ F62asz0;  //Z=15227
    double z12v[200], a1v[200], b1v[200], b2v[200], b1sv[200], /*sum12[200],*/ sum22[200], sum32[200]; //, sum42[200], sum52[200], sum62[200];  //Z=15228

    /*begin*/  //Z=15230
    zl = (1-sigmal*sigmal)/(sigmal*sigmal);  //Z=15231
    zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232
    radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=15233

    // Noch fehlende (globale) Variablen und sonstige Anpassungen:
    // An machen Stellen muss "params." eingefügt werden.
    // Die Konstante 'eps' muss durch 'eps9' ersetzt werden.
    // Bei machen if-else muss noch ein ';' eingefügt werden.
    double zz = zr;
    double z  = zl;
    double argpq, c, rad, xmax, pqr1, pqr2, pqr3, qnarg, binsum1, pq;
    double qxn[200], qyn[200];
    double qz = 1; // wird für qrombchid() verwendet (zuerst bei Z=16029)
    // Dort werden auch p11..p33, ax1[nxyz],ax2[nxyz],ax3[nxyz], sig[xyz] verwendet

    //ella = params.radius;  //Z=15235
    //ellb = params.length;  //Z=15236
    //ellc = radiusm;  //Z=15237

    /* ************ */  //Z=15240
    /* ** sphere ** */  //Z=15241
    /* ************ */  //Z=15242
    if ( params.part==0 )
    {/*2*/  //Z=15243
        /* ** homogeneous sphere ** */  //Z=15244
        if ( params.cs==0 )
        {/*3*/  //Z=15245
            //if ( q > 2.2 ) qDebug() << "  formpq" << 0.4*params.limq4 << q << "r"<<params.radius;
            if ( q<0.4*params.limq4 )
            {/*4*/  //Z=15246
                pqsum = 1.0;  //Z=15247
                oldpqsum = 0.0;  //Z=15248
                /* qqn[0]:=1.0;  //Z=15249 */
                qq2 = 1.0;  //Z=15250
                for ( nser=1; nser<=100; nser++ )
                {/*5*/  //Z=15251
                    /* qqn[nser]:=qqn[nser-1]*q*q;  //Z=15252 */
                    qq2 = qq2*q*q;  //Z=15253
                    /* pqsum:=pqsum+carr4p[nser]*qqn[nser];  //Z=15254 */
                    pqsum = pqsum+params.CR->carr4p[nser]*qq2;  //Z=15255
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15256
                    if ( delser<0.0001 ) break; /* goto 50; */  //Z=15257
                    oldpqsum = pqsum;  //Z=15258
                }/*5*/  //Z=15259
                /*50:*/  //Z=15260
                /*formpq:=*/ return pqsum;  //Z=15261
            }/*4*/  //Z=15262
            else
            {/*4*/  //Z=15263
                argq = q*params.radius/(zr+1);  //Z=15264
                pqr = (1/(2.0*zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4);  //Z=15265
                pq1 = pqr*(1+cos((zr-3)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-3)/2.0));  //Z=15266
                pq2 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-4)/2.0);  //Z=15267
                pq3 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2.0*argq))/pow(1.0+4*argq*argq,zr-5)/2.0);  //Z=15268
                /*formpq:=*/ return 9*(pq1-2*pq2+pq3);  //Z=15269
            }/*4*/  //Z=15270
        }/*3*/ /*  of homogeneous sphere */  //Z=15271

        /* ** core/shell sphere ** */  //Z=15273
        if ( params.cs==1 )
        {/*3*/  //Z=15274

            cc1 = sqr(params.rho);  //Z=15276
            cc2 = 2*params.p1*params.rho*(1-params.rho);  //Z=15277
            cc3 = sqr(1-params.rho)*sqr(params.p1);  //Z=15278
            cc4 = -2*sqr(params.rho);  //Z=15279
            cc5 = -2*params.p1*params.rho*(1-params.rho);  //Z=15280
            cc6 = sqr(params.rho);  //Z=15281
            cc7 = -2*params.rho*(1-params.rho);  //Z=15282
            cc8 = -sqr(1-params.rho)*2*params.p1;  //Z=15283
            cc9 = 2*params.rho*(1-params.rho);  //Z=15284
            cc10 = sqr(1-params.rho);  //Z=15285

            ccc1 = sqr(1-params.rho)*pow(params.p1,6);  //Z=15287
            ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,3);  //Z=15288
            ccc3 = params.rho*params.rho;  //Z=15289
            vv3 = sqr((1-params.rho)*pow(params.p1,3)+params.rho);  //Z=15290

            argq = q*radiusm/(zz+1);  //Z=15292
            argpq = q*params.radius/(zz+1);  //Z=15293
            pqr = (1/(2.0*zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4);  //Z=15294

            /*  F121 sphere  */  //Z=15296
            if ( q<(0.3*params.limq4) )
            {/*4*/  //Z=15297
                qqn[0] = 1.0;  //Z=15298
                pqsum = 1.0;  //Z=15299
                oldpqsum = 0.0;  //Z=15300
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=15301
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=15302
                    pqsum = pqsum+qqn[nser]*params.CR->carr4p[nser];  //Z=15303
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15304
                    if ( delser<0.0001 ) break; /* goto 51; */  //Z=15305
                    oldpqsum = pqsum;  //Z=15306
                }/*5*/  //Z=15307
                /*51:*/  //Z=15308
                F121 = ccc1*pqsum/vv3;  //Z=15309
            }/*4*/  //Z=15310
            else
            {/*4*/  //Z=15311
                ac3 = pqr*(1+cos((zr-3)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-3)/2.0));  //Z=15312
                ac8 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-4)/2.0);  //Z=15313
                ac10 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-5)/2.0));  //Z=15314
                F121 = 9*(cc3*ac3+cc8*ac8+cc10*ac10)/vv3;  //Z=15315
            }/*4*/  //Z=15316

            /*  F122 sphere  */  //Z=15318
            if ( q<(0.3*params.limq5) )
            {/*4*/  //Z=15319
                qqn[0] = 1.0;  //Z=15320
                pqsum = 1.0;  //Z=15321
                oldpqsum = 0.0;  //Z=15322
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=15323
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=15324
                    pqsum = pqsum+qqn[nser]*params.CR->carr5p[nser];  //Z=15325
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15326
                    if ( delser<0.0001 ) break; /* goto 52; */  //Z=15327
                    oldpqsum = pqsum;  //Z=15328
                }/*5*/  //Z=15329
                /*52:*/  //Z=15330
                F122 = ccc2*pqsum/vv3;  //Z=15331
            }/*4*/  //Z=15332
            else
            {/*4*/  //Z=15333
                argbm = (zr-3)*atan(argpq-argq);  //Z=15334
                nenbm = pow(1.0+sqr(argpq-argq),(zr-3)/2.0);  //Z=15335
                argbp = (zr-3)*atan(argpq+argq);  //Z=15336
                nenbp = pow(1.0+sqr(argpq+argq),(zr-3)/2.0);  //Z=15337
                ac2 = pqr*(cos(argbm)/nenbm+cos(argbp)/nenbp);  //Z=15338
                argep = (zr-4)*atan(argpq+argq);  //Z=15339
                nenep = pow(1.0+sqr(argpq+argq),(zr-4)/2.0);  //Z=15340
                argem = (zr-4)*atan(argpq-argq);  //Z=15341
                nenem = pow(1.0+sqr(argpq-argq),(zr-4)/2.0);  //Z=15342
                ac5 = (pqr/((zr-4)*argq))*(sin(argep)/nenep-sin(argem)/nenem);  //Z=15343
                arggp = (zr-4)*atan(argpq+argq);  //Z=15344
                nengp = pow(1.0+sqr(argpq+argq),(zr-4)/2.0);  //Z=15345
                arggm = (zr-4)*atan(argq-argpq);  //Z=15346
                nengm = pow(1.0+sqr(argq-argpq),(zr-4)/2.0);  //Z=15347
                ac7 = (pqr/((zr-4)*argq))*(sin(arggp)/nengp-sin(arggm)/nengm);  //Z=15348
                //argim = (zr-5)*atan(argpq-argq);  //Z=15349
                //nenim = pow(1.0+sqr(argpq-argq),(zr-5)/2.0);  //Z=15350
                //argip = (zr-5)*atan(argpq+argq);  //Z=15351
                //nenip = pow(1.0+sqr(argpq+argq),(zr-5)/2.0);  //Z=15352
                ac9 = (pqr/((zr-4)*(zr-5)*argq*argq))*(cos(argbm)/nenbm-cos(argbp)/nenbp);  //Z=15353

                F122 = 9*(cc2*ac2+cc5*ac5+cc7*ac7+cc9*ac9)/vv3;  //Z=15355
            }/*4*/  //Z=15356

            /*  F123 sphere  */  //Z=15358
            if ( q<(0.3*params.limq6) )
            {/*4*/  //Z=15359
                qqn[0] = 1.0;  //Z=15360
                pqsum = 1.0;  //Z=15361
                oldpqsum = 0.0;  //Z=15362
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=15363
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=15364
                    pqsum = pqsum+qqn[nser]*params.CR->carr6p[nser];  //Z=15365
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15366
                    if ( delser<0.0001 ) break; /* goto 53; */  //Z=15367
                    oldpqsum = pqsum;  //Z=15368
                }/*5*/  //Z=15369
                /*53:*/  //Z=15370
                F123 = ccc3*pqsum/vv3;  //Z=15371
            }/*4*/  //Z=15372
            else
            {/*4*/  //Z=15373
                ac1 = pqr*(1+cos((zr-3)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-3)/2.0));  //Z=15374
                ac4 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-4)/2.0);  //Z=15375
                ac6 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-5)/2.0));  //Z=15376
                F123 = 9*(cc1*ac1+cc4*ac4+cc6*ac6)/vv3;  //Z=15377
            }/*4*/  //Z=15378
            /*formpq:=*/ return F121+F122+F123;  //Z=15379
        }/*3*/  /*  of core/shell sphere  */  //Z=15380

        /* ** inhomogeneous core/shell sphere ** */  //Z=15382
        if ( params.cs==2 )
        {/*3*/  //Z=15383

            dim = 3;  //Z=15385
            delc = 0.0001;  //Z=15386
            xrad = q*radiusm;  //Z=15387
            xradp = q*params.radius;  //Z=15388
            x1z = q*params.radius/(2.0*(zr+1));  //Z=15389
            x12z = x1z*x1z;  //Z=15390
            x2z = q*radiusm/(2.0*(zr+1));  //Z=15391
            x22z = x2z*x2z;  //Z=15392

            lim = 18*exp(-5*params.sigma);  //Z=15394
            lim1 = lim;  //Z=15395
            lim2 = lim*0.7;  //Z=15396
            lim3 = lim;  //Z=15397
            lim4 = lim;  //Z=15398
            lim5 = lim*0.7;  //Z=15399
            lim6 = lim*1.2;  //Z=15400

            a1 = (dim-params.alphash1)/2.0;  //Z=15402
            b1 = dim/2.0;  //Z=15403
            b2 = (dim+2-params.alphash1)/2.0;  //Z=15404
            b1s = (dim+2)/2.0;  //Z=15405
            v = -b1s+1/2.0;  //Z=15406
            c = a1-b1-b2+1/2.0;  //Z=15407
            d0 = 1;  //Z=15408
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=15409
            e0 = 1.0;  //Z=15410
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=15411
            ee0 = 1.0;  //Z=15412
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=15413

            gb1s = 3*sqrt(M_PI)/4.0;  //Z=15415
            pz2v = 1/(zr*(zr-1)*(zr-2)*(zr-3));  //Z=15416
            pz2v1 = pz2v/(zr-4);  //Z=15417
            pz2v2 = pz2v1/(zr-5);  //Z=15418

            gz1 = gamma(zr+1);  //Z=15420
            preg1 = gb1s/sqrt(M_PI);  //Z=15421
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=15422
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=15423
            pzvc = gamma(zr+1+v+c)/gz1;  //Z=15424
            pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=15425
            pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=15426
            pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=15427
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=15428
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=15429
            pzc = gamma(zr+1+2*c)/gz1;  //Z=15430
            pzc1 = gamma(zr+1+2*c-1)/gz1;  //Z=15431
            pza = gamma(zr+1-4*a1)/gz1;  //Z=15432
            pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=15433
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=15434
            //dnv0 = 1;  //Z=15435
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=15436
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=15437
            //pva0 = gamma(zr+1-4*a1)/gz1;  //Z=15438

            cc1 = 1/(dim*dim);  //Z=15440
            cc2 = 2*params.rho/(dim*(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=15441
            cc3 = -2*params.rho/(dim*(dim-params.alphash1));  //Z=15442
            cc4 = sqr(params.rho)/(sqr(dim-params.alphash1)*pow(sqr(params.p1),dim-params.alphash1));  //Z=15443
            cc5 = -2*sqr(params.rho)/(sqr(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=15444
            cc6 = sqr(params.rho)/sqr(dim-params.alphash1);  //Z=15445
            vv3 = cc1+cc2+cc3+cc4+cc5+cc6;  //Z=15446

            /*  term #1 series  */  //Z=15448
            if ( (xradp)<lim1 )
            {/*4*/  //Z=15449
                //z12v[0] = 1;  //Z=15450
                //b1sv[0] = 1;  //Z=15451
                //fkv[0] = 1;  //Z=15452
                //gam3[0] = sqrt(M_PI)/2.0;  //Z=15453
                double qqnn = 1.0; // qqn[0] = 1.0;  //Z=15454
                F12sez = 1.0;  //Z=15455
                oldF12sez = 0.0;  //Z=15456
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=15457
                    //qqn[n] = qqn[n-1]*q*q;  //Z=15458
                    qqnn = qqnn * q * q;
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=15459
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=15460
                    //fkv[n] = fkv[n-1]*n;  //Z=15461
                    //gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=15462
                    //sum12[n] = 0;  //Z=15463
                    /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);  //Z=15464 */
                    /* sum12[n]:=(9*sqrt(pi)/2)*power(4,n)/((n+3)*(n+2)*(n+3/2)*gam3[n]*fkv[n]);  //Z=15465 */
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n];  //Z=15466 */

                    F12sez += params.CR->carr1p[n]*qqnn; //qqn[n];  //Z=15468

                    del = fabs((F12sez-oldF12sez)/F12sez);  //Z=15470
                    if ( del<delc ) break; /* goto 201; */  //Z=15471
                    oldF12sez = F12sez;  //Z=15472
                }/*5*/  //Z=15473
                /*201:*/  //Z=15474
                F12 = F12sez;  //Z=15475
            }/*4*/  //Z=15476

            /*  term #2 series  */  //Z=15478
            if ( (xradp)<lim2 )
            {/*4*/  //Z=15479
                //z12v[0] = 1;  //Z=15480
                //a1v[0] = 1;  //Z=15481
                //b1v[0] = 1;  //Z=15482
                //b2v[0] = 1;  //Z=15483
                //b1sv[0] = 1;  //Z=15484
                //fkv[0] = 1;  //Z=15485
                //gam3[0] = sqrt(M_PI)/2.0;  //Z=15486
                double qqnn = 1.0; // qqn[0] = 1.0;  //Z=15487
                F22sez = 1.0;  //Z=15488
                oldF22sez = 0.0;  //Z=15489
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=15490
                    //qqn[n] = qqn[n-1]*q*q;  //Z=15491
                    qqnn = qqnn * q * q;
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=15492
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=15493
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=15494
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=15495
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=15496
                    //fkv[n] = fkv[n-1]*n;  //Z=15497
                    //gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=15498
                    //sum22[n] = 0;  //Z=15499
                    /* for m:=0 to n do sum22[n]:=sum22[n]+a1v[n-m]*power(p1*p1,m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=15500 */
                    /* F22sez:=F22sez+power(-x22z,n)*z12v[n]*sum22[n];  //Z=15501 */

                    /* for m:=0 to n do sum22[n]:=sum22[n]+power(p1*p1,m)/((n-m+3/2)*(m+(3/2)-(alfa/2))*gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=15503 */
                    /* F22sez:=F22sez+(3*pi*(3-alfa)/16)*power(-x22z,n)*z12v[n]*sum22[n];  //Z=15504 */

                    F22sez += params.CR->carr2p[n]*qqnn; // qqn[n];  //Z=15506

                    del = fabs((F22sez-oldF22sez)/F22sez);  //Z=15508
                    if ( del<delc ) break; /* goto 202; */  //Z=15509
                    oldF22sez = F22sez;  //Z=15510
                }/*5*/  //Z=15511
                /*202:*/  //Z=15512
                F22 = F22sez;  //Z=15513
            }/*4*/  //Z=15514

            /*  term #3 series  */  //Z=15516
            if ( (xradp)<lim3 )
            {/*4*/  //Z=15517
                //z12v[0] = 1;  //Z=15518
                //a1v[0] = 1;  //Z=15519
                //b1v[0] = 1;  //Z=15520
                //b2v[0] = 1;  //Z=15521
                //b1sv[0] = 1;  //Z=15522
                //fkv[0] = 1;  //Z=15523
                //gam3[0] = sqrt(M_PI)/2.0;  //Z=15524
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=15525
                F32sez = 1.0;  //Z=15526
                oldF32sez = 0.0;  //Z=15527
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=15528
                    //qqn[n] = qqn[n-1]*q*q;  //Z=15529
                    qqnn = qqnn * q * q;
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=15530
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=15531
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=15532
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=15533
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=15534
                    //fkv[n] = fkv[n-1]*n;  //Z=15535
                    //gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=15536
                    //sum32[n] = 0;  //Z=15537
                    /* for m:=0 to n do sum32[n]:=sum32[n]+a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=15538 */
                    /* F32sez:=F32sez+power(-x12z,n)*z12v[n]*sum32[n];  //Z=15539 */

                    /* for m:=0 to n do sum32[n]:=sum32[n]+1/((n-m+3/2)*(m+(3/2)-(alfa/2))*gam3[m]*gam3[n-m]*fkv[m]*fkv[n-m]);  //Z=15541 */
                    /* F32sez:=F32sez+(3*pi*(3-alfa)/16)*power(-x12z,n)*z12v[n]*sum32[n];  //Z=15542 */

                    F32sez += params.CR->carr3p[n]*qqnn; //qqn[n];  //Z=15544

                    del = fabs((F32sez-oldF32sez)/F32sez);  //Z=15546
                    if ( del<delc ) break; /* goto 203; */  //Z=15547
                    oldF32sez = F32sez;  //Z=15548
                }/*5*/  //Z=15549
                /*203:*/  //Z=15550
                F32 = F32sez;  //Z=15551
            }/*4*/  //Z=15552

            /*  term #4 series  */  //Z=15554
            if ( (xradp)<lim4 )
            {/*4*/  //Z=15555
                //z12v[0] = 1;  //Z=15556
                //a1v[0] = 1;  //Z=15557
                //b1v[0] = 1;  //Z=15558
                //b2v[0] = 1;  //Z=15559
                //b1sv[0] = 1;  //Z=15560
                //fkv[0] = 1;  //Z=15561
                //gam3[0] = sqrt(M_PI)/2.0;  //Z=15562
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=15563
                F42sez = 1.0;  //Z=15564
                oldF42sez = 0.0;  //Z=15565
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=15566
                    //qqn[n] = qqn[n-1]*q*q;  //Z=15567
                    qqnn = qqnn * q * q;
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=15568
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=15569
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=15570
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=15571
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=15572
                    //fkv[n] = fkv[n-1]*n;  //Z=15573
                    //gam3[n] = gam3[n-1]*(2*n+1)/2.0;  //Z=15574
                    //sum42[n] = 0;  //Z=15575
                    /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=15576 */
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n];  //Z=15577 */

                    /* for m:=0 to n do sum42[n]:=sum42[n]+1/((n-m+(3/2)-(alfa/2))*(m+(3/2)-(alfa/2))*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=15579 */
                    /* F42sez:=F42sez+((3-alfa)*(3-alfa)*pi/16)*power(-x22z,n)*z12v[n]*sum42[n];  //Z=15580 */

                    F42sez += params.CR->carr4p[n]*qqnn; //qqn[n];  //Z=15582

                    del = fabs((F42sez-oldF42sez)/F42sez);  //Z=15584
                    if ( del<delc ) break; /* goto 204; */  //Z=15585
                    oldF42sez = F42sez;  //Z=15586
                }/*5*/  //Z=15587
                /*204:*/  //Z=15588
                F42 = F42sez;  //Z=15589
            }/*4*/  //Z=15590

            /*  term #5 series  */  //Z=15592
            if ( (xradp)<lim5 )
            {/*4*/  //Z=15593
                //z12v[0] = 1;  //Z=15594
                //a1v[0] = 1;  //Z=15595
                //b1v[0] = 1;  //Z=15596
                //b2v[0] = 1;  //Z=15597
                //b1sv[0] = 1;  //Z=15598
                //fkv[0] = 1;  //Z=15599
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=15600
                F52sez = 1.0;  //Z=15601
                oldF52sez = 0.0;  //Z=15602
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=15603
                    //qqn[n] = qqn[n-1]*q*q;  //Z=15604
                    qqnn = qqnn * q * q;
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=15605
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=15606
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=15607
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=15608
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=15609
                    //fkv[n] = fkv[n-1]*n;  //Z=15610
                    //sum52[n] = 0;  //Z=15611
                    /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=15612 */
                    /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n];  //Z=15613 */

                    /* for m:=0 to n do sum52[n]:=sum52[n]+power(p1*p1,n-m)/((n-m+(3/2)-(alfa/2))*(m+(3/2)-(alfa/2))*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=15615 */
                    /* F52sez:=F52sez+((3-alfa)*(3-alfa)*pi/16)*power(-x22z,n)*z12v[n]*sum52[n];  //Z=15616 */

                    F52sez += params.CR->carr5p[n]*qqnn; //qqn[n];  //Z=15618

                    del = fabs((F52sez-oldF52sez)/F52sez);  //Z=15620
                    if ( del<delc ) break; /* goto 205; */  //Z=15621
                    oldF52sez = F52sez;  //Z=15622
                }/*5*/  //Z=15623
                /*205:*/  //Z=15624
                F52 = F52sez;  //Z=15625
            }/*4*/  //Z=15626

            /*  term #6 series  */  //Z=15628
            if ( (xradp)<lim6 )
            {/*4*/  //Z=15629
                //z12v[0] = 1;  //Z=15630
                //a1v[0] = 1;  //Z=15631
                //b1v[0] = 1;  //Z=15632
                //b2v[0] = 1;  //Z=15633
                //b1sv[0] = 1;  //Z=15634
                //fkv[0] = 1;  //Z=15635
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=15636
                F62sez = 1.0;  //Z=15637
                oldF62sez = 0.0;  //Z=15638
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=15639
                    //qqn[n] = qqn[n-1]*q*q;  //Z=15640
                    qqnn = qqnn * sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=15641
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=15642
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=15643
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=15644
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=15645
                    //fkv[n] = fkv[n-1]*n;  //Z=15646
                    //sum62[n] = 0;  //Z=15647
                    /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=15648 */
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n];  //Z=15649 */

                    /* for m:=0 to n do sum62[n]:=sum62[n]+1/((n-m+(3/2)-(alfa/2))*(m+(3/2)-(alfa/2))*gam3[n-m]*gam3[m]*fkv[m]*fkv[n-m]);  //Z=15651 */
                    /* F62sez:=F62sez+((3-alfa)*(3-alfa)*pi/16)*power(-x12z,n)*z12v[n]*sum62[n];  //Z=15652 */

                    F62sez += params.CR->carr6p[n]*qqnn; //qqn[n];  //Z=15654

                    del = fabs((F62sez-oldF62sez)/F62sez);  //Z=15656
                    if ( del<delc ) break; /* goto 206; */  //Z=15657
                    oldF62sez = F62sez;  //Z=15658
                }/*5*/  //Z=15659
                /*206:*/  //Z=15660
                F62 = F62sez;  //Z=15661
            }/*4*/  //Z=15662


            /* ** term #1 asymptote ** */  //Z=15665
            if ( xradp>=lim1 )
            {/*4*/  //Z=15666
                arg11 = (zr+2*v+1)*atan(4.0*x1z);  //Z=15667
                nen11 = pow(1.0+16*x1z*x1z,(zr+2*v+1)/2.0);  //Z=15668
                arg12 = (zr+2*v)*atan(4.0*x1z);  //Z=15669
                nen12 = pow(1.0+16*x1z*x1z,(zr+2*v)/2.0);  //Z=15670
                arg13 = (zr+2*v-1)*atan(4.0*x1z);  //Z=15671
                nen13 = pow(1.0+16*x1z*x1z,(zr+2*v-1)/2.0);  //Z=15672

                F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);  //Z=15674
                F12as2z = 2*ee0*ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);  //Z=15675
                F12as3z = ee1*ee1*(1/(4.0*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);  //Z=15676
                F12asz = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);  //Z=15677
                F12 = F12asz;  //Z=15678
            }/*4*/  //Z=15679

            /* ** term #2 asymptote ** */  //Z=15681
            if ( xradp>=lim2 )
            {/*4*/  //Z=15682
                //arg21 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=15683
                //nen21 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=15684
                //arg22 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=15685
                //nen22 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=15686
                //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);  //Z=15687
                //F22as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);  //Z=15688
                F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);  //Z=15689
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);  //Z=15690

                arg210 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=15692
                nen210 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=15693
                arg220 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=15694
                nen220 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=15695
                F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);  //Z=15696
                F22as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);  //Z=15697
                F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);  //Z=15698
                arg23 = (zr+v+c+1)*atan(2.0*(x1z-x2z));  //Z=15699
                nen23 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);  //Z=15700
                arg24 = (zr+v+c+1)*atan(2.0*(x1z+x2z));  //Z=15701
                nen24 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);  //Z=15702
                arg25 = (zr+v+c)*atan(2.0*(x1z-x2z));  //Z=15703
                nen25 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);  //Z=15704
                arg26 = (zr+v+c)*atan(2.0*(x1z+x2z));  //Z=15705
                nen26 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);  //Z=15706
                arg27 = (zr+v+c-1)*atan(2.0*(x1z-x2z));  //Z=15707
                nen27 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);  //Z=15708
                arg28 = (zr+v+c-1)*atan(2.0*(x1z+x2z));  //Z=15709
                nen28 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);  //Z=15710

                a22as21z = (1/2.0)*ee0*e0*pzvc;  //Z=15712
                F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);  //Z=15713
                a22as22z = (1/2.0)*ee0*e1*(1/(2.0*x2z))*pzvc1;  //Z=15714
                F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);  //Z=15715
                a22as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1;  //Z=15716
                F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);  //Z=15717
                a22as24z = (1/2.0)*ee1*e1*(1/(2.0*x1z))*(1/(2.0*x2z))*pzvc2;  //Z=15718
                F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);  //Z=15719
                F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);  //Z=15720
                F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);  //Z=15721
                //F22asz = F22as1z+F22as2z;  //Z=15722
                F22asz0 = F22as1z0+F22as2z;  //Z=15723
                F22 = F22asz0;  //Z=15724
            }/*4*/  //Z=15725

            /* ** term #3 asymptote ** */  //Z=15727
            if ( xradp>=lim3 )
            {/*4*/  //Z=15728
                //arg31 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=15729
                //nen31 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=15730
                //arg32 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=15731
                //nen32 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=15732
                //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);  //Z=15733
                //F32as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);  //Z=15734
                F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);  //Z=15735
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);  //Z=15736

                arg310 = (z+v-2*a1+1)*atan(2.0*x1z);  //Z=15738
                nen310 = pow(1.0+4*x1z*x1z,(z+v-2*a1+1)/2.0);  //Z=15739
                arg320 = (z+v-2*a1)*atan(2.0*x1z);  //Z=15740
                nen320 = pow(1.0+4*x1z*x1z,(z+v-2*a1)/2.0);  //Z=15741
                F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);  //Z=15742
                F32as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);  //Z=15743
                F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);  //Z=15744

                arg33 = (zr+v+c+1)*atan(4.0*x1z);  //Z=15746
                nen33 = pow(1.0+16*x1z*x1z,(zr+v+c+1)/2.0);  //Z=15747
                arg34 = (zr+v+c)*atan(4.0*x1z);  //Z=15748
                nen34 = pow(1.0+16*x1z*x1z,(zr+v+c)/2.0);  //Z=15749
                arg35 = (zr+v+c-1)*atan(4.0*x1z);  //Z=15750
                nen35 = pow(1.0+16*x1z*x1z,(zr+v+c-1)/2.0);  //Z=15751
                F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);  //Z=15752
                F32as22z = (1/2.0)*ee0*e1*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);  //Z=15753
                F32as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);  //Z=15754
                F32as24z = (1/2.0)*ee1*e1*(1/(4.0*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);  //Z=15755
                F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);  //Z=15756
                F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);  //Z=15757
                //F32asz = F32as1z+F32as2z;  //Z=15758
                F32asz0 = F32as1z0+F32as2z;  //Z=15759
                F32 = F32asz0;  //Z=15760
            }/*4*/  //Z=15761


            /* ** term #4 asymptote ** */  //Z=15764
            if ( xrad>=lim4 )
            {/*4*/  //Z=15765
                F42as10z = preg4*preg4*pow(x22z,-2*a1);  //Z=15766
                //F42as1sumz = pva0;  //Z=15767
                //F42as1z = F42as10z*F42as1sumz;  //Z=15768
                F42as1z0 = F42as10z*pza;  //Z=15769

                arg41 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=15771
                nen41 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=15772
                arg42 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=15773
                nen42 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=15774
                //arg43 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=15775
                //nen43 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=15776
                F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=15777
                F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=15778
                F42as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);  //Z=15779
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);  //Z=15780
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);  //Z=15781
                F42as2z0 = F42as20z*(F42as21+F42as22);  //Z=15782

                F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=15784
                F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=15785
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);  //Z=15786
                F42as26 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);  //Z=15787
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);  //Z=15788
                F42as3z0 = F42as30z*(F42as24+F42as26);  //Z=15789

                F42as40z = preg3*preg3*pow(x2z*x2z,c);  //Z=15791
                arg44 = (zr+2*c+1)*atan(4.0*x2z);  //Z=15792
                nen44 = pow(1.0+16*x2z*x2z,(zr+2*c+1)/2.0);  //Z=15793
                arg45 = (zr+2*c)*atan(4.0*x2z);  //Z=15794
                nen45 = pow(1.0+16*x2z*x2z,(zr+2*c)/2.0);  //Z=15795
                F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);  //Z=15796
                F42as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=15797
                F42as29 = (1/2.0)*e1*e0*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=15798
                F42as4z = F42as40z*(F42as27+F42as28+F42as29);  //Z=15799
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;  //Z=15800
                F42asz0 = F42as1z0+F42as2z0+F42as3z0+F42as4z;  //Z=15801
                F42 = F42asz0;  //Z=15802
            }/*4*/  //Z=15803


            /* ** term #5 asymptote ** */  //Z=15806
            if ( xradp>=lim5 )
            {/*4*/  //Z=15807
                F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);  //Z=15808
                //F52as1sumz = pva0;  //Z=15809
                //F52as1z = F52as10z*F52as1sumz;  //Z=15810
                F52as1z0 = F52as10z*pza;  //Z=15811

                F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);  //Z=15813
                arg51 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=15814
                nen51 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=15815
                arg52 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=15816
                nen52 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=15817
                //arg53 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=15818
                //nen53 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=15819
                F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);  //Z=15820
                F52as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);  //Z=15821
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);  //Z=15822
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);  //Z=15823
                F52as2z0 = F52as20z*(F52as21+F52as22);  //Z=15824

                F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);  //Z=15826
                arg54 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=15827
                nen54 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=15828
                //arg55 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=15829
                //nen55 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=15830
                arg56 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=15831
                nen56 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=15832
                F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);  //Z=15833
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);  //Z=15834
                F52as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);  //Z=15835
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);  //Z=15836
                F52as3z0 = F52as30z*(F52as24+F52as26);  //Z=15837

                F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);  //Z=15839
                arg57 = (zr+2*c+1)*atan(2.0*(x1z-x2z));  //Z=15840
                nen57 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);  //Z=15841
                arg58 = (zr+2*c+1)*atan(2.0*(x1z+x2z));  //Z=15842
                nen58 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);  //Z=15843
                arg59 = (zr+2*c)*atan(2.0*(x1z-x2z));  //Z=15844
                nen59 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);  //Z=15845
                arg510 = (zr+2*c)*atan(2.0*(x1z+x2z));  //Z=15846
                nen510 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);  //Z=15847
                F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);  //Z=15848
                F52as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=15849
                F52as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=15850
                F52as4z = F52as40z*(F52as27+F52as28+F52as29);  //Z=15851
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;  //Z=15852
                F52asz0 = F52as1z0+F52as2z0+F52as3z0+F52as4z;  //Z=15853
                F52 = F52asz0;  //Z=15854
            }/*4*/  //Z=15855

            /* ** term #6 asymptote ** */  //Z=15857
            if ( xradp>=lim6 )
            {/*4*/  //Z=15858
                F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);  //Z=15859
                //F62as1sumz = pva0;  //Z=15860
                //F62as1z = F62as10z*F62as1sumz;  //Z=15861
                F62as1z0 = F62as10z*pza;  //Z=15862

                F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=15864
                arg61 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=15865
                nen61 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=15866
                arg62 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=15867
                nen62 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=15868
                //arg63 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=15869
                //nen63 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=15870
                F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=15871
                F62as22 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=15872
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);  //Z=15873
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);  //Z=15874
                F62as2z0 = F62as20z*(F62as21+F62as22);  //Z=15875

                F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=15877
                F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=15878
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);  //Z=15879
                F62as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=15880
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);  //Z=15881
                F62as3z0 = F62as30z*(F62as24+F62as26);  //Z=15882

                F62as40z = preg3*preg3*pow(x1z*x1z,c);  //Z=15884
                arg64 = (zr+2*c+1)*atan(4.0*x1z);  //Z=15885
                nen64 = pow(1.0+16*x1z*x1z,(zr+2*c+1)/2.0);  //Z=15886
                arg65 = (zr+2*c)*atan(4.0*x1z);  //Z=15887
                nen65 = pow(1.0+16*x1z*x1z,(zr+2*c)/2.0);  //Z=15888
                F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);  //Z=15889
                F62as28 = (1/2.0)*e0*e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=15890
                F62as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=15891
                F62as4z = F62as40z*(F62as27+F62as28+F62as29);  //Z=15892
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;  //Z=15893
                F62asz0 = F62as1z0+F62as2z0+F62as3z0+F62as4z;  //Z=15894
                F62 = F62asz0;  //Z=15895
            }/*4*/  //Z=15896

            /*formpq:=*/ return (cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;  //Z=15898


            /* formpq:=pqcoreshellin(1.0,rho,p1,1.0,0.001,alfa,radiusm,3,sigmar,q);  //Z=15901 */
        }/*3*/ /*  of inhomogeneous core/shell sphere  */  //Z=15902



        /*  myelin sphere  */  //Z=15906
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=15907

            /*  sphere parameters  */  //Z=15909
            v = -2;  //Z=15910
            e0 = 1;  //Z=15911
            e1 = -1;  //Z=15912
            preg1 = 3/4.0;  //Z=15913
            pz2v = 1/(zr*(zr-1)*(zr-2)*(zr-3));  //Z=15914
            pz2v1 = pz2v/(zr-4);  //Z=15915
            pz2v2 = pz2v1/(zr-5);  //Z=15916
            lim = 18*exp(-5*params.sigma);  //Z=15917
            lim1 = lim*1.4;  //Z=15918
            rad = params.CR->myarray[1];  //Z=15919
            inmax = round(params.CR->myarray[14]);  //Z=15920
            vvm = params.CR->myarray[15];  //Z=15921
            rmax = params.CR->myarray[16];  //Z=15922
            xmax = q*rmax;  //Z=15923

            if ( xmax<(lim1) )
            {/*4*/  //Z=15925
                /* fkv[0]:=1;  //Z=15926 */
                qqn[0] = 1.0;  //Z=15927
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=15928
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=15929
                    /* fkv[nser]:=fkv[nser-1]*nser;  //Z=15930 */
                }/*5*/  //Z=15931

                F12sum = 0.0;  //Z=15933
                for ( ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=15934
                    for ( jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=15935
                        F12sez = 1.0;  //Z=15936
                        oldF12sez = 1.0;  //Z=15937
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=15938
                            pqsum = 0;  //Z=15939
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=15940
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]);  //Z=15941 */
                                pqsum = pqsum+pow(params.CR->carr7p[ii],2*mser)*pow(params.CR->carr7p[jj],2*(nser-mser))/(params.CR->carr6p[mser]*params.CR->carr6p[nser-mser]);  //Z=15942

                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=15944 */
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx];  //Z=15945 */
                            }/*8*/  //Z=15946
                            F12sez = F12sez+params.CR->carr4p[nser]*qqn[nser]*pqsum;  //Z=15947
                            delser = fabs((F12sez-oldF12sez)/F12sez);  //Z=15948
                            if ( delser<0.0001 ) break; /* goto 250; */  //Z=15949
                            oldF12sez = F12sez;  //Z=15950
                        }/*7*/  //Z=15951
                        /*250:*/  //Z=15952
                        F12sum = F12sum+params.CR->carr5p[ii]*params.CR->carr5p[jj]*F12sez;  //Z=15953
                    }/*6*/  //Z=15954
                }/*5*/  //Z=15955
                F12ser = F12sum/vvm;  //Z=15956
                F12 = F12ser;  //Z=15957
            }/*4*/  //Z=15958
            else
            {/*4*/  //Z=15959
                xrz = q*rad/(zr+1);  //Z=15960
                arg = (zr+2*v+1)*atan(2.0*xrz);  //Z=15961
                nen = pow(1.0+4*xrz*xrz,(zr+2*v+1)/2.0);  //Z=15962
                arg1 = (zr+2*v)*atan(2.0*xrz);  //Z=15963
                nen1 = pow(1.0+4*xrz*xrz,(zr+2*v)/2.0);  //Z=15964
                arg2 = (zr+2*v-1)*atan(2.0*xrz);  //Z=15965
                nen2 = pow(1.0+4*xrz*xrz,(zr+2*v-1)/2.0);  //Z=15966

                F12asz = 0.0;  //Z=15968
                for ( ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=15969
                    a1m = params.CR->carr5p[ii]*pow(params.CR->carr7p[ii],v);   /*  carr7p[ii]:=pp[ii];  //Z=15970 */
                    for ( jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=15971
                        a2m = params.CR->carr5p[jj]*pow(params.CR->carr7p[jj],v);  //Z=15972
                        xijm = (params.CR->carr3p[ii]-params.CR->carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii];  //Z=15973 */
                        arglmz = (zr+1)*atan(xijm);  //Z=15974
                        nenlmz = pow(1.0+xijm*xijm,(zr+1)/2.0);  //Z=15975
                        xijp = (params.CR->carr3p[ii]+params.CR->carr3p[jj])*q/(zr+1);  //Z=15976
                        arglpz = (zr+1)*atan(xijp);  //Z=15977
                        nenlpz = pow(1.0+xijp*xijp,(zr+1)/2.0);  //Z=15978
                        F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(M_PI*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(M_PI*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=15979
                        F12as2z = e0*e1*(1/(params.CR->carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=15980
                        F12as3z = e1*e0*(1/(params.CR->carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=15981
                        F12as4z = e1*e1*(1/(params.CR->carr7p[ii]*params.CR->carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(M_PI*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(M_PI*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=15982

                        F12asz = F12asz+a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);  //Z=15984
                    }/*6*/  //Z=15985
                }/*5*/  //Z=15986
                F12asy = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;  //Z=15987
                F12 = F12asy;  //Z=15988
            }/*4*/  //Z=15989
            /*formpq:=*/ return F12;  //Z=15990

            /* formpq:=polyliposome(llipt,radius,lliph,lin,lout,nom,sigmar,sigmal,phiax,philiph,philipt,phiin,phiout,2,q);  //Z=15992 */
            /* formpq:=polyliposome(2.0,200,1.0,3.5,3.5,1,sigmar,sigmal,0.001,-0.55,-0.7,0.001,0.001,3,q);  //Z=15993 */
            /* formpq:=pql;  //Z=15994 */
        }/*3*/ /*  of myelin sphere  */  //Z=15995

    }/*2*/ /*  of sphere  */  //Z=15997



    /* ************ */  //Z=16001
    /* ** triaxial ellipsoid ** */  //Z=16002
    /* ************ */  //Z=16003
    if ( params.part==116 )
    {/*2*/  //Z=16004
        epsi = sigmal;  //Z=16005
        /* ** homogeneous triaxial ellipsoid ** */  //Z=16006
        if ( ordis==7 )
        {/*3*/   /*  isotropic  */  //Z=16007
            if ( params.cs==0 )
            {/*4*/  //Z=16008
                if ( q<0.4*params.limq4 )
                {/*5*/  //Z=16009
                    pqsum = 1.0;  //Z=16010
                    oldpqsum = 0.0;  //Z=16011
                    qqn[0] = 1.0;  //Z=16012
                    for ( nser=1; nser<=100; nser++ )
                    {/*6*/  //Z=16013
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=16014
                        pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=16015
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16016
                        if ( delser<0.0001 ) break; /* goto 90; */  //Z=16017
                        oldpqsum = pqsum;  //Z=16018
                    }/*6*/  //Z=16019
                    /*90:*/  //Z=16020
                    /*formpq:=*/ return pqsum;  //Z=16021
                }/*5*/  //Z=16022
                else
                {/*5*/  //Z=16023
                    /* if (q>(4*limq4)) then pq:=(3*pi/(4*zr*(zr-1)*(zr-2)*(zr-3)))*power(q*radius/(zr+1),-4)  //Z=16024 */
                    /*    else begin  //Z=16025 */

                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,orcase,0,0,0,carr1p,pq);  //Z=16027 */

                    qrombchid(params.length,params.radius,params.p1,params.sigma,params.alphash1,epsi,params.polTheta,params.polPhi,qx,qy,qz,
                              params.p11,params.p12,params.p13,params.p21,params.p22,params.p23,params.p31,params.p32,params.p33,
                              qx,qy,0,qhkl,
                              params.ax1.length(),params.ax2.length(),params.ax3.length(),
                              params.ax1.x(),params.ax1.y(),params.ax1.z(),
                              params.ax2.x(),params.ax2.y(),params.ax2.z(),
                              params.ax3.x(),params.ax3.y(),params.ax3.z(),
                              params.sig.x(),params.sig.y(),params.sig.z(),
                              ordis,3,7,13,7,0,0,params.CR->carr1p,pql);  //Z=16029
                    /*formpq:=*/ return pql;  //Z=16030
                }/*5*/  //Z=16031
            }/*4*/ /*  of homogeneous sphere */  //Z=16032
        }/*3*/  /*  of isotropic  */  //Z=16033


        /*  perfect  */  //Z=16036
        if ( ordis==6 )
        {/*3*/  //Z=16037
            if ( params.orcase==4 )
                pql = 1.0;  //Z=16038
            else
            {/*4*/  //Z=16039
                if ( q<(0.6*params.limq1) )
                {/*5*/  //Z=16040
                    pqsum = 1.0;  //Z=16041
                    oldpqsum = 0.0;  //Z=16042
                    qxn[0] = 1.0;  //Z=16043
                    qyn[0] = 1.0;  //Z=16044
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=16045
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=16046
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=16047
                        binsum = 0.0;  //Z=16048
                        for ( mser=0; mser<=nser; mser++ )  //Z=16049
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[mser]*qyn[nser-mser];  //Z=16050
                        pqsum = pqsum+binsum;  //Z=16051
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16052
                        if ( delser<0.0001 ) break; /* goto 91; */  //Z=16053
                        oldpqsum = pqsum;  //Z=16054
                    }/*6*/  //Z=16055
                    /*91:*/  //Z=16056
                    pql = pqsum;  //Z=16057
                }/*5*/  //Z=16058
                else
                {/*5*/  //Z=16059
                    arglq = (qxs+qys+eps9)*params.length/(zl+1);  //Z=16060
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2));  //Z=16061 */
                    /* pql:=(pi/(2*zl))*(1/arglq);  //Z=16062 */
                    pql = (1/(2.0*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-1)/2.0));  //Z=16063
                }/*5*/  //Z=16064
            }/*4*/  //Z=16065
        }/*3*/   /*  of perfect  */  //Z=16066

    }/*2*/   /*  of ellipsoid  */  //Z=16068



    /* ********** */  //Z=16072
    /*  cylinder  */  //Z=16073
    /* ********** */  //Z=16074
    if ( params.part==1 )
    {/*2*/  //Z=16075

        /* ** longitudinal part ** */  //Z=16077
        /* ** isotropic ** */  //Z=16078
        if ( ordis==7 )
        {/*3*/  //Z=16079
            /*  exact average  */  //Z=16080
            if ( (params.length/params.radius)<2 )
            {/*4*/  //Z=16081
                if ( q<(5*params.limq2) )
                {/*5*/  //Z=16082
                    /*  Cauchy sum  */  //Z=16083
                    /* pqsum:=1.0;  //Z=16084
                       oldpqsum:=0.0;  //Z=16085
                       qqn[0]:=1.0;  //Z=16086
                       for nser:=1 to 120 do begin  //Z=16087
                          qqn[nser]:=qqn[nser-1]*q*q;  //Z=16088
                          pqsum:=pqsum+carr3p[nser]*qqn[nser];  //Z=16089
                          delser:=abs((pqsum-oldpqsum)/pqsum);  //Z=16090
                          if delser<0.0001 then goto 601;  //Z=16091
                          oldpqsum:=pqsum;  //Z=16092
                       end;  */  //Z=16093

                    /*  double sum  */  //Z=16095
                    qqn[0] = 1.0;  //Z=16096
                    for ( nser=1; nser<=100; nser++ ) qqn[nser] = qqn[nser-1]*q*q;  //Z=16097
                    pqsum = 0.0;  //Z=16098
                    oldpqsum = -10.0;  //Z=16099
                    for ( nser=0; nser<=100; nser++ )
                    {/*6*/  //Z=16100
                        binsum = 0.0;  //Z=16101
                        for ( mser=0; mser<=100; mser++ ) binsum = binsum+params.CR->carr11pm[nser][mser]*qqn[mser];  //Z=16102
                        pqsum = pqsum+params.CR->carr2p[nser]*qqn[nser]*binsum;  //Z=16103
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16104
                        if ( delser<0.0001 ) break; /* goto 601; */  //Z=16105
                        oldpqsum = pqsum;  //Z=16106
                    }/*6*/  //Z=16107

                    /*601:*/  //Z=16109
                    pql = pqsum;  //Z=16110
                }/*5*/  //Z=16111
                else
                {/*5*/  //Z=16112
                    pql = params.por/pow(q,4);  //Z=16113
                }/*5*/  //Z=16114
            }/*4*/  //Z=16115

            /*  factorization  */  //Z=16117
            else
            {/*4*/  //Z=16118
                if ( q<(0.6*params.limq1) )
                {/*5*/  //Z=16119
                    pqsum = 1.0;  //Z=16120
                    oldpqsum = 0.0;  //Z=16121
                    qqn[0] = 1.0;  //Z=16122
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=16123
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=16124
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser];  //Z=16125
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16126
                        if ( delser<0.0001 ) break; /* goto 60; */  //Z=16127
                        oldpqsum = pqsum;  //Z=16128
                    }/*6*/  //Z=16129
                    /*60:*/  //Z=16130
                    pql = pqsum;  //Z=16131
                }/*5*/  //Z=16132
                else
                {/*5*/  //Z=16133
                    arglq = q*params.length/(zl+1);  //Z=16134
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2));  //Z=16135 */
                    pql = (M_PI/(2.0*zl))*(1/arglq);  //Z=16136
                    pql = pql-(1/(2.0*zl*(zl-1)*arglq*arglq))*cos((zl-1)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-1)/2.0);  //Z=16137
                }/*5*/  //Z=16138
            }/*4*/  //Z=16139
        }/*3*/   /*  of isotropic  */  //Z=16140

        /*  perfect  */  //Z=16142
        if ( ordis==6 )
        {/*3*/  //Z=16143
            if ( params.orcase==4 )
                pql = 1.0;  //Z=16144
            else
            {/*4*/  //Z=16145
                if ( limql<(0.5*params.limq1) )
                {/*5*/  //Z=16146
                    /* if (sqrt(qx*qx*length*length+qy*qy*radius*radius+eps)<10) then begin  //Z=16147 */
                    pqsum = 1.0;  //Z=16148
                    oldpqsum = 0.0;  //Z=16149
                    qqn[0] = 1.0;  //Z=16150
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=16151
                        qqn[nser] = qqn[nser-1]*(qxs+qys)*(qxs+qys);     /*  (qxs,qys)=(qx,0) for x, (0,qy) for y, (qx,qy) for z  */  //Z=16152
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser];  //Z=16153
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16154
                        if ( delser<0.0001 ) break; /* goto 65; */  //Z=16155
                        oldpqsum = pqsum;  //Z=16156
                    }/*6*/  //Z=16157
                    /*65:*/  //Z=16158
                    pql = pqsum;  //Z=16159
                }/*5*/  //Z=16160
                else
                {/*5*/  //Z=16161
                    arglq = (qxs+qys+eps9)*params.length/(zl+1);  //Z=16162
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2));  //Z=16163 */
                    /* pql:=(pi/(2*zl))*(1/arglq);  //Z=16164 */
                    pql = (1/(2.0*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-1)/2.0));  //Z=16165
                }/*5*/  //Z=16166
            }/*4*/  //Z=16167
        }/*3*/   /*  of perfect  */  //Z=16168

        /*  general  */  //Z=16170
        if ( ordis==0 )
        {/*3*/  //Z=16171
            if ( params.orcase==4 )
                pql = 1.0;  //Z=16172
            else
            {/*4*/  //Z=16173
                lim = params.CR->myarray[17];  //Z=16174
                if ( limql<(2*params.limq1) )
                {/*5*/  //Z=16175
                    /* if (limql<0.1) then begin  //Z=16176 */
                    /* if (q<0.05) then begin  //Z=16177 */
                    pqsum = 1.0;  //Z=16178
                    oldpqsum = 0.0;  //Z=16179
                    qxn[0] = 1.0;  //Z=16180
                    qyn[0] = 1.0;  //Z=16181
                    if ( params.orcase==1 )
                    {/*6*/  //Z=16182
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16183
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=16184
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=16185
                            binsum = 0.0;  //Z=16186
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16187
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16188 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=16189 */
                                binsum = binsum+params.CR->carr11pm[nser][mser]*qxn[mser]*qyn[nser-mser];  //Z=16190
                            }/*8*/  //Z=16191
                            pqsum = pqsum+params.CR->carr1p[nser]*binsum;  //Z=16192
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16193
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=16194
                            oldpqsum = pqsum;  //Z=16195
                        }/*7*/  //Z=16196
                    }/*6*/  //Z=16197

                    if ( params.orcase==2 )
                    {/*6*/  /*  x-axis  */  //Z=16199
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16200
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=16201
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=16202
                            binsum = 0.0;  //Z=16203
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16204
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16205 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=16206 */
                                binsum = binsum+params.CR->carr11pm[nser][mser]*qxn[mser]*qyn[nser-mser];  //Z=16207
                            }/*8*/  //Z=16208
                            pqsum = pqsum+params.CR->carr1p[nser]*binsum;  //Z=16209
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16210
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=16211
                            oldpqsum = pqsum;  //Z=16212
                        }/*7*/  //Z=16213
                    }/*6*/  //Z=16214
                    pql = pqsum;  //Z=16215

                    if ( params.orcase==3 )
                    {/*6*/  /*  y-axis  */  //Z=16217
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16218
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=16219
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=16220
                            binsum = 0.0;  //Z=16221
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16222
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16223 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=16224 */
                                binsum = binsum+params.CR->carr11pm[nser][mser]*qxn[mser]*qyn[nser-mser];  //Z=16225
                            }/*8*/  //Z=16226
                            pqsum = pqsum+params.CR->carr1p[nser]*binsum;  //Z=16227
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16228
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=16229
                            oldpqsum = pqsum;  //Z=16230
                        }/*7*/  //Z=16231
                    }/*6*/  //Z=16232
                    /*66:*/  //Z=16233
                    pql = pqsum;  //Z=16234
                }/*5*/  //Z=16235
                else
                {/*5*/  //Z=16236
                    qrombdeltac(params.p1,sigmal,params.alphash1,params.polTheta,0,qxs,qys,qz,
                                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,params.orcase,0,0,0,params.CR->carr1p,pql);  //Z=16237
                    pql = pql/params.norm;  //Z=16238
                }/*5*/  //Z=16239
            }/*4*/  //Z=16240
        }/*3*/   /*  of general  */  //Z=16241


        /*  transverse part  */  //Z=16244
        /*  homogeneous cylinder  */  //Z=16245
        if ( params.cs==0 )
        {/*3*/  //Z=16246
            /*  exact average  */  //Z=16247
            if ( (params.length/params.radius)<2 )
                pqr = 1;  //Z=16248
            /*  factorization  */  //Z=16249
            else
            {/*4*/  //Z=16250
                if ( q<(0.3*params.limq4) )
                {/*5*/  //Z=16251
                    /* if (sqrt(qx*qx*length*length+qy*qy*radius*radius+eps)<10) then begin  //Z=16252 */
                    pqsum = 1.0;  //Z=16253
                    oldpqsum = 0.0;  //Z=16254
                    qqn[0] = 1.0;  //Z=16255
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=16256
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=16257
                        /* qqn[nser]:=qqn[nser-1]*(qxs+qys)*(qxs+qys);     (* (qxs,qys)=(qx,0) for x, (0,qy) for y, (qx,qy) for z *)  //Z=16258 */
                        /* qqn[nser]:=qqn[nser-1]*qx*qx;       (* GISAXS *)  //Z=16259 */
                        pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=16260
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16261
                        if ( delser<0.0001 ) break; /* goto 61; */  //Z=16262
                        oldpqsum = pqsum;  //Z=16263
                    }/*6*/  //Z=16264
                    /*61:*/  //Z=16265
                    pqr = pqsum;  //Z=16266
                }/*5*/  //Z=16267
                else
                {/*5*/  //Z=16268
                    /* q:=sqrt(qxs*qxs+qys*qys+eps);       (* for perfect orientation *)  //Z=16269 */
                    /* q:=abs(qx+eps);                               (* for GISAXS *)  //Z=16270 */
                    argpq = q*params.radius/(zr+1);  //Z=16271
                    pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3);  //Z=16272
                    pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3)*sin((zr-2)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-2)/2.0);  //Z=16273
                    pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argpq,-4)*cos((zr-3)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-3)/2.0);  //Z=16274
                    pqr = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=16275
                    /*  add more terms, if necessary  */  //Z=16276
                }/*5*/  //Z=16277
            }/*4*/  //Z=16278
            /*formpq:=*/ return pql*pqr;  //Z=16279
            /* formpq:=pqr;  //Z=16280 */
            /* formpq:=pql;  //Z=16281 */
        }/*3*/ /*  of homogeneous  */  //Z=16282

        /*  core/shell cylinder  */  //Z=16284
        if ( params.cs==1 )
        {/*3*/  //Z=16285
            cc1 = sqr(params.rho);  //Z=16286
            cc2 = 2*params.p1*params.rho*(1-params.rho);  //Z=16287
            cc3 = sqr(1-params.rho)*sqr(params.p1);  //Z=16288
            cc4 = -2*sqr(params.rho);  //Z=16289
            cc5 = -2*params.p1*params.rho*(1-params.rho);  //Z=16290
            cc6 = sqr(params.rho);  //Z=16291
            cc7 = -2*params.rho*(1-params.rho);  //Z=16292
            cc8 = -sqr(1-params.rho)*2*params.p1;  //Z=16293
            cc9 = 2*params.rho*(1-params.rho);  //Z=16294
            cc10 = sqr(1-params.rho);  //Z=16295

            ccc1 = sqr(1-params.rho)*pow(params.p1,4);  //Z=16297
            ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,2);  //Z=16298
            ccc3 = params.rho*params.rho;  //Z=16299
            vv3 = sqr((1-params.rho)*pow(params.p1,2)+params.rho);  //Z=16300

            argq = q*radiusm/(zz+1);  //Z=16302
            argpq = q*params.radius/(zz+1);  //Z=16303

            /*  F121 cylinder  */  //Z=16305
            if ( q<(0.7*params.limq4) )
            {/*4*/  //Z=16306
                /* ** series expansion ** */  //Z=16307
                pqsum = 1.0;  //Z=16308
                oldpqsum = 0.0;  //Z=16309
                qqn[0] = 1.0;  //Z=16310
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=16311
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=16312
                    pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=16313
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16314
                    if ( delser<0.0001 ) break; /* goto 62; */  //Z=16315
                    oldpqsum = pqsum;  //Z=16316
                }/*5*/  //Z=16317
                /*62:*/  //Z=16318
                F121 = ccc1*pqsum/vv3;  //Z=16319
            }/*4*/  //Z=16320
            else
            {/*4*/  //Z=16321
                pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3);  //Z=16322
                pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3)*sin((zr-2)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-2)/2.0);  //Z=16323
                pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argpq,-4)*cos((zr-3)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-3)/2.0);  //Z=16324
                pqr = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=16325
                F121 = ccc1*pqr/vv3;  //Z=16326
                /*  add more terms, if necessary  */  //Z=16327
            }/*4*/  //Z=16328

            /*  F122 cylinder  */  //Z=16330
            if ( q<(1.5*params.limq5) )
            {/*4*/  //Z=16331
                /* ** series expansion ** */  //Z=16332
                pqsum = 1.0;  //Z=16333
                oldpqsum = 0.0;  //Z=16334
                qqn[0] = 1.0;  //Z=16335
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=16336
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=16337
                    pqsum = pqsum+params.CR->carr5p[nser]*qqn[nser];  //Z=16338
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16339
                    if ( delser<0.0001 ) break; /* goto 63; */  //Z=16340
                    oldpqsum = pqsum;  //Z=16341
                }/*5*/  //Z=16342
                /*63:*/  //Z=16343
                F122 = ccc2*pqsum/vv3;  //Z=16344
            }/*4*/  //Z=16345
            else
            {/*4*/  //Z=16346
                argbm = (zr-2)*atan(argpq-argq);  //Z=16347
                nenbm = pow(1.0+sqr(argpq-argq),(zr-2)/2.0);  //Z=16348
                argbp = (zr-2)*atan(argpq+argq);  //Z=16349
                nenbp = pow(1.0+sqr(argpq+argq),(zr-2)/2.0);  //Z=16350
                argem = (zr-3)*atan(argpq-argq);  //Z=16351
                nenem = pow(1.0+sqr(argpq-argq),(zr-3)/2.0);  //Z=16352
                argep = (zr-3)*atan(argpq+argq);  //Z=16353
                nenep = pow(1.0+sqr(argpq+argq),(zr-3)/2.0);  //Z=16354
                arggm = (zr-4)*atan(argpq-argq);  //Z=16355
                nengm = pow(1.0+sqr(argpq-argq),(zr-4)/2.0);  //Z=16356
                arggp = (zr-4)*atan(argpq+argq);  //Z=16357
                nengp = pow(1.0+sqr(argpq+argq),(zr-4)/2.0);  //Z=16358

                pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3)*(cos(argbm)/nenbm-sin(argbp)/nenbp);  //Z=16360
                pqr2 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4)*((1-1/params.p1)*sin(argem)/nenem-(1+1/params.p1)*cos(argep)/nenep);  //Z=16361
                pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)*(zr-4)))*pow(argq,-5)*(1/params.p1)*(cos(arggm)/nengm-sin(arggp)/nengp);  //Z=16362
                pqr = (4/M_PI)*pow(params.p1,-3/2.0)*(pqr1+(9/16.0)*pqr2+(9/16.0)*(9/16.0)*pqr3);  //Z=16363
                F122 = ccc2*pqr/vv3;  //Z=16364
            }/*4*/  //Z=16365

            /*  F123 cylinder  */  //Z=16367
            if ( q<(0.6*params.limq6) )
            {/*4*/  //Z=16368
                /* ** series expansion ** */  //Z=16369
                pqsum = 1.0;  //Z=16370
                oldpqsum = 0.0;  //Z=16371
                qqn[0] = 1.0;  //Z=16372
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=16373
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=16374
                    pqsum = pqsum+params.CR->carr6p[nser]*qqn[nser];  //Z=16375
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16376
                    if ( delser<0.0001 ) break; /* goto 64; */  //Z=16377
                    oldpqsum = pqsum;  //Z=16378
                }/*5*/  //Z=16379
                /*64:*/  //Z=16380
                F123 = ccc3*pqsum/vv3;  //Z=16381
            }/*4*/  //Z=16382
            else
            {/*4*/  //Z=16383
                pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3);  //Z=16384
                pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3)*sin((zr-2)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-2)/2.0);  //Z=16385
                pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4)*cos((zr-3)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-3)/2.0);  //Z=16386
                pqr = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=16387
                F123 = ccc3*pqr/vv3;  //Z=16388
                /*  add more terms, if necessary  */  //Z=16389
            }/*4*/  //Z=16390
            /*formpq:=*/ return pql*(F121+F122+F123);  //Z=16391
            /* formpq:=F122;  //Z=16392 */
        }/*3*/ /*  of core/shell  */  //Z=16393

        /*  inhomogeneous core/shell cylinder  */  //Z=16395
        if ( params.cs==2 )
        {/*3*/  //Z=16396

            dim = 2;  //Z=16398
            delc = 0.0001;  //Z=16399
            xrad = q*radiusm;  //Z=16400
            xradp = q*params.radius;  //Z=16401
            x1z = q*params.radius/(2.0*(zr+1));  //Z=16402
            x12z = x1z*x1z;  //Z=16403
            x2z = q*radiusm/(2.0*(zr+1));  //Z=16404
            x22z = x2z*x2z;  //Z=16405

            lim = 18*exp(-5*params.sigma);  //Z=16407
            lim1 = lim;  //Z=16408
            lim2 = lim*0.7;  //Z=16409
            lim3 = lim;  //Z=16410
            lim4 = lim;  //Z=16411
            lim5 = lim*0.7;  //Z=16412
            lim6 = lim*1.2;  //Z=16413

            a1 = (dim-params.alphash1)/2.0;  //Z=16415
            b1 = dim/2.0;  //Z=16416
            b2 = (dim+2-params.alphash1)/2.0;  //Z=16417
            b1s = (dim+2)/2.0;  //Z=16418
            v = -b1s+1/2.0;  //Z=16419
            c = a1-b1-b2+1/2.0;  //Z=16420
            d0 = 1;  //Z=16421
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=16422
            e0 = 1.0;  //Z=16423
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=16424
            ee0 = 1.0;  //Z=16425
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=16426

            gb1s = 1;  //Z=16428
            pz2v = 1/(zr*(zr-1)*(zr-2));  //Z=16429
            pz2v1 = pz2v/(zr-3);  //Z=16430
            pz2v2 = pz2v1/(zr-4);  //Z=16431

            gz1 = gamma(zr+1);  //Z=16433
            preg1 = gb1s/sqrt(M_PI);  //Z=16434
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=16435
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=16436
            pzvc = gamma(zr+1+v+c)/gz1;  //Z=16437
            pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=16438
            pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=16439
            pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=16440
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=16441
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=16442
            pzc = gamma(zr+1+2*c)/gz1;  //Z=16443
            pzc1 = gamma(zr+1+2*c-1)/gz1;  //Z=16444
            pza = gamma(zr+1-4*a1)/gz1;  //Z=16445
            pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=16446
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=16447
            //dnv0 = 1;  //Z=16448
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=16449
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=16450
            //pva0 = gamma(zr+1-4*a1)/gz1;  //Z=16451

            cc1 = 1/(dim*dim);  //Z=16453
            cc2 = 2*params.rho/(dim*(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=16454
            cc3 = -2*params.rho/(dim*(dim-params.alphash1));  //Z=16455
            cc4 = sqr(params.rho)/(sqr(dim-params.alphash1)*pow(sqr(params.p1),dim-params.alphash1));  //Z=16456
            cc5 = -2*sqr(params.rho)/(sqr(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=16457
            cc6 = sqr(params.rho)/sqr(dim-params.alphash1);  //Z=16458
            vv3 = cc1+cc2+cc3+cc4+cc5+cc6;  //Z=16459

            /*  term #1 series  */  //Z=16461
            if ( (xradp)<lim1 )
            {/*4*/  //Z=16462
                //z12v[0] = 1;  //Z=16463
                //b1sv[0] = 1;  //Z=16464
                //fkv[0] = 1;  //Z=16465
                F12sez = 1.0;  //Z=16466
                oldF12sez = 0.0;  //Z=16467
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16468
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=16469
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16470
                    qqnn = qqnn * sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16471
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16472
                    //fkv[n] = fkv[n-1]*n;  //Z=16473
                    //sum12[n] = 0;  //Z=16474
                    /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);  //Z=16475 */
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n];  //Z=16476 */

                    F12sez += params.CR->carr4p[n]*qqnn; //qqn[n];  //Z=16478

                    del = fabs((F12sez-oldF12sez)/F12sez);  //Z=16480
                    if ( del<delc ) break; /* goto 211; */  //Z=16481
                    oldF12sez = F12sez;  //Z=16482
                }/*5*/  //Z=16483
                /*211:*/  //Z=16484
                F12 = F12sez;  //Z=16485
            }/*4*/  //Z=16486

            /*  term #2 series  */  //Z=16488
            if ( (xradp)<lim2 )
            {/*4*/  //Z=16489
                z12v[0] = 1;  //Z=16490
                a1v[0] = 1;  //Z=16491
                b1v[0] = 1;  //Z=16492
                b2v[0] = 1;  //Z=16493
                b1sv[0] = 1;  //Z=16494
                fkv[0] = 1;  //Z=16495
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16496
                F22sez = 1.0;  //Z=16497
                oldF22sez = 0.0;  //Z=16498
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=16499
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16500
                    qqnn = qqnn * sqr(q);
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16501
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16502
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16503
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16504
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16505
                    fkv[n] = fkv[n-1]*n;  //Z=16506
                    sum22[n] = 0;  //Z=16507
                    for ( m=0; m<=n; m++ ) sum22[n] = sum22[n]+a1v[n-m]*pow(sqr(params.p1),m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16508
                    F22sez = F22sez+pow(-x22z,n)*z12v[n]*sum22[n];  //Z=16509

                    /* F22sez:=F22sez+carr5p[n]*qqn[n];  //Z=16511 */

                    del = fabs((F22sez-oldF22sez)/F22sez);  //Z=16513
                    if ( del<delc ) break; /* goto 212; */  //Z=16514
                    oldF22sez = F22sez;  //Z=16515
                }/*5*/  //Z=16516
                /*212:*/  //Z=16517
                F22 = F22sez;  //Z=16518
            }/*4*/  //Z=16519

            /*  term #3 series  */  //Z=16521
            if ( (xradp)<lim3 )
            {/*4*/  //Z=16522
                z12v[0] = 1;  //Z=16523
                a1v[0] = 1;  //Z=16524
                b1v[0] = 1;  //Z=16525
                b2v[0] = 1;  //Z=16526
                b1sv[0] = 1;  //Z=16527
                fkv[0] = 1;  //Z=16528
                qqn[0] = 1.0;  //Z=16529
                F32sez = 1.0;  //Z=16530
                oldF32sez = 0.0;  //Z=16531
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=16532
                    qqn[n] = qqn[n-1]*q*q;  //Z=16533
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16534
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16535
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16536
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16537
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16538
                    fkv[n] = fkv[n-1]*n;  //Z=16539
                    sum32[n] = 0;  //Z=16540
                    for ( m=0; m<=n; m++ ) sum32[n] = sum32[n]+a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16541
                    F32sez = F32sez+pow(-x12z,n)*z12v[n]*sum32[n];  //Z=16542

                    /* F32sez:=F32sez+carr6p[n]*qqn[n];  //Z=16544 */

                    del = fabs((F32sez-oldF32sez)/F32sez);  //Z=16546
                    if ( del<delc ) break; /* goto 213; */  //Z=16547
                    oldF32sez = F32sez;  //Z=16548
                }/*5*/  //Z=16549
                /*213:*/  //Z=16550
                F32 = F32sez;  //Z=16551
            }/*4*/  //Z=16552

            /*  term #4 series  */  //Z=16554
            if ( (xradp)<lim4 )
            {/*4*/  //Z=16555
                //z12v[0] = 1;  //Z=16556
                //a1v[0] = 1;  //Z=16557
                //b1v[0] = 1;  //Z=16558
                //b2v[0] = 1;  //Z=16559
                //b1sv[0] = 1;  //Z=16560
                //fkv[0] = 1;  //Z=16561
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16562
                F42sez = 1.0;  //Z=16563
                oldF42sez = 0.0;  //Z=16564
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=16565
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16566
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16567
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16568
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16569
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16570
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16571
                    //fkv[n] = fkv[n-1]*n;  //Z=16572
                    //sum42[n] = 0;  //Z=16573
                    /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16574 */
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n];  //Z=16575 */

                    F42sez += params.CR->carr7p[n]*qqnn; //qqn[n];  //Z=16577

                    del = fabs((F42sez-oldF42sez)/F42sez);  //Z=16579
                    if ( del<delc ) break; /* goto 214; */  //Z=16580
                    oldF42sez = F42sez;  //Z=16581
                }/*5*/  //Z=16582
                /*214:*/  //Z=16583
                F42 = F42sez;  //Z=16584
            }/*4*/  //Z=16585

            /*  term #5 series  */  //Z=16587
            if ( (xradp)<lim5 )
            {/*4*/  //Z=16588
                //z12v[0] = 1;  //Z=16589
                //a1v[0] = 1;  //Z=16590
                //b1v[0] = 1;  //Z=16591
                //b2v[0] = 1;  //Z=16592
                //b1sv[0] = 1;  //Z=16593
                //fkv[0] = 1;  //Z=16594
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16595
                F52sez = 1.0;  //Z=16596
                oldF52sez = 0.0;  //Z=16597
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=16598
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16599
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16600
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16601
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16602
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16603
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16604
                    //fkv[n] = fkv[n-1]*n;  //Z=16605
                    //sum52[n] = 0;  //Z=16606
                    /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16607 */
                    /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n];  //Z=16608 */

                    F52sez += params.CR->carr8p[n]*qqnn; //qqn[n];  //Z=16610

                    del = fabs((F52sez-oldF52sez)/F52sez);  //Z=16612
                    if ( del<delc ) break; /* goto 215; */  //Z=16613
                    oldF52sez = F52sez;  //Z=16614
                }/*5*/  //Z=16615
                /*215:*/  //Z=16616
                F52 = F52sez;  //Z=16617
            }/*4*/  //Z=16618

            /*  term #6 series  */  //Z=16620
            if ( (xradp)<lim6 )
            {/*4*/  //Z=16621
                //z12v[0] = 1;  //Z=16622
                //a1v[0] = 1;  //Z=16623
                //b1v[0] = 1;  //Z=16624
                //b2v[0] = 1;  //Z=16625
                //b1sv[0] = 1;  //Z=16626
                //fkv[0] = 1;  //Z=16627
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16628
                F62sez = 1.0;  //Z=16629
                oldF62sez = 0.0;  //Z=16630
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=16631
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16632
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16633
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16634
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16635
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16636
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16637
                    //fkv[n] = fkv[n-1]*n;  //Z=16638
                    //sum62[n] = 0;  //Z=16639
                    /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16640 */
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n];  //Z=16641 */

                    F62sez += params.CR->carr9p[n]*qqnn; //qqn[n];  //Z=16643

                    del = fabs((F62sez-oldF62sez)/F62sez);  //Z=16645
                    if ( del<delc ) break; /* goto 216; */  //Z=16646
                    oldF62sez = F62sez;  //Z=16647
                }/*5*/  //Z=16648
                /*216:*/  //Z=16649
                F62 = F62sez;  //Z=16650
            }/*4*/  //Z=16651


            /* ** term #1 asymptote ** */  //Z=16654
            if ( xradp>=lim1 )
            {/*4*/  //Z=16655
                arg11 = (zr+2*v+1)*atan(4.0*x1z);  //Z=16656
                nen11 = pow(1.0+16*x1z*x1z,(zr+2*v+1)/2.0);  //Z=16657
                arg12 = (zr+2*v)*atan(4.0*x1z);  //Z=16658
                nen12 = pow(1.0+16*x1z*x1z,(zr+2*v)/2.0);  //Z=16659
                arg13 = (zr+2*v-1)*atan(4.0*x1z);  //Z=16660
                nen13 = pow(1.0+16*x1z*x1z,(zr+2*v-1)/2.0);  //Z=16661

                F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);  //Z=16663
                F12as2z = 2*ee0*ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);  //Z=16664
                F12as3z = ee1*ee1*(1/(4.0*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);  //Z=16665
                F12asz = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);  //Z=16666
                F12 = F12asz;  //Z=16667
            }/*4*/  //Z=16668

            /* ** term #2 asymptote ** */  //Z=16670
            if ( xradp>=lim2 )
            {/*4*/  //Z=16671
                //arg21 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=16672
                //nen21 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=16673
                //arg22 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=16674
                //nen22 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=16675
                //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);  //Z=16676
                //F22as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);  //Z=16677
                F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);  //Z=16678
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);  //Z=16679

                arg210 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=16681
                nen210 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=16682
                arg220 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=16683
                nen220 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=16684
                F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);  //Z=16685
                F22as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);  //Z=16686
                F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);  //Z=16687
                arg23 = (zr+v+c+1)*atan(2.0*(x1z-x2z));  //Z=16688
                nen23 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);  //Z=16689
                arg24 = (zr+v+c+1)*atan(2.0*(x1z+x2z));  //Z=16690
                nen24 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);  //Z=16691
                arg25 = (zr+v+c)*atan(2.0*(x1z-x2z));  //Z=16692
                nen25 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);  //Z=16693
                arg26 = (zr+v+c)*atan(2.0*(x1z+x2z));  //Z=16694
                nen26 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);  //Z=16695
                arg27 = (zr+v+c-1)*atan(2.0*(x1z-x2z));  //Z=16696
                nen27 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);  //Z=16697
                arg28 = (zr+v+c-1)*atan(2.0*(x1z+x2z));  //Z=16698
                nen28 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);  //Z=16699

                a22as21z = (1/2.0)*ee0*e0*pzvc;  //Z=16701
                F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);  //Z=16702
                a22as22z = (1/2.0)*ee0*e1*(1/(2.0*x2z))*pzvc1;  //Z=16703
                F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);  //Z=16704
                a22as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1;  //Z=16705
                F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);  //Z=16706
                a22as24z = (1/2.0)*ee1*e1*(1/(2.0*x1z))*(1/(2.0*x2z))*pzvc2;  //Z=16707
                F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);  //Z=16708
                F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);  //Z=16709
                F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);  //Z=16710
                //F22asz = F22as1z+F22as2z;  //Z=16711
                F22asz0 = F22as1z0+F22as2z;  //Z=16712
                F22 = F22asz0;  //Z=16713
            }/*4*/  //Z=16714

            /* ** term #3 asymptote ** */  //Z=16716
            if ( xradp>=lim3 )
            {/*4*/  //Z=16717
                //arg31 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=16718
                //nen31 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=16719
                //arg32 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=16720
                //nen32 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=16721
                //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);  //Z=16722
                //F32as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);  //Z=16723
                F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);  //Z=16724
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);  //Z=16725

                arg310 = (z+v-2*a1+1)*atan(2.0*x1z);  //Z=16727
                nen310 = pow(1.0+4*x1z*x1z,(z+v-2*a1+1)/2.0);  //Z=16728
                arg320 = (z+v-2*a1)*atan(2.0*x1z);  //Z=16729
                nen320 = pow(1.0+4*x1z*x1z,(z+v-2*a1)/2.0);  //Z=16730
                F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);  //Z=16731
                F32as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);  //Z=16732
                F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);  //Z=16733

                arg33 = (zr+v+c+1)*atan(4.0*x1z);  //Z=16735
                nen33 = pow(1.0+16*x1z*x1z,(zr+v+c+1)/2.0);  //Z=16736
                arg34 = (zr+v+c)*atan(4.0*x1z);  //Z=16737
                nen34 = pow(1.0+16*x1z*x1z,(zr+v+c)/2.0);  //Z=16738
                arg35 = (zr+v+c-1)*atan(4.0*x1z);  //Z=16739
                nen35 = pow(1.0+16*x1z*x1z,(zr+v+c-1)/2.0);  //Z=16740
                F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);  //Z=16741
                F32as22z = (1/2.0)*ee0*e1*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);  //Z=16742
                F32as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);  //Z=16743
                F32as24z = (1/2.0)*ee1*e1*(1/(4.0*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);  //Z=16744
                F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);  //Z=16745
                F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);  //Z=16746
                //F32asz = F32as1z+F32as2z;  //Z=16747
                F32asz0 = F32as1z0+F32as2z;  //Z=16748
                F32 = F32asz0;  //Z=16749
            }/*4*/  //Z=16750


            /* ** term #4 asymptote ** */  //Z=16753
            if ( xrad>=lim4 )
            {/*4*/  //Z=16754
                F42as10z = preg4*preg4*pow(x22z,-2*a1);  //Z=16755
                //F42as1sumz = pva0;  //Z=16756
                //F42as1z = F42as10z*F42as1sumz;  //Z=16757
                F42as1z0 = F42as10z*pza;  //Z=16758

                arg41 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=16760
                nen41 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=16761
                arg42 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=16762
                nen42 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=16763
                //arg43 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=16764
                //nen43 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=16765
                F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=16766
                F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=16767
                F42as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);  //Z=16768
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);  //Z=16769
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);  //Z=16770
                F42as2z0 = F42as20z*(F42as21+F42as22);  //Z=16771

                F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=16773
                F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=16774
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);  //Z=16775
                F42as26 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);  //Z=16776
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);  //Z=16777
                F42as3z0 = F42as30z*(F42as24+F42as26);  //Z=16778

                F42as40z = preg3*preg3*pow(x2z*x2z,c);  //Z=16780
                arg44 = (zr+2*c+1)*atan(4.0*x2z);  //Z=16781
                nen44 = pow(1.0+16*x2z*x2z,(zr+2*c+1)/2.0);  //Z=16782
                arg45 = (zr+2*c)*atan(4.0*x2z);  //Z=16783
                nen45 = pow(1.0+16*x2z*x2z,(zr+2*c)/2.0);  //Z=16784
                F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);  //Z=16785
                F42as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=16786
                F42as29 = (1/2.0)*e1*e0*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=16787
                F42as4z = F42as40z*(F42as27+F42as28+F42as29);  //Z=16788
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;  //Z=16789
                F42asz0 = F42as1z0+F42as2z0+F42as3z0+F42as4z;  //Z=16790
                F42 = F42asz0;  //Z=16791
            }/*4*/  //Z=16792


            /* ** term #5 asymptote ** */  //Z=16795
            if ( xradp>=lim5 )
            {/*4*/  //Z=16796
                F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);  //Z=16797
                //F52as1sumz = pva0;  //Z=16798
                //F52as1z = F52as10z*F52as1sumz;  //Z=16799
                F52as1z0 = F52as10z*pza;  //Z=16800

                F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);  //Z=16802
                arg51 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=16803
                nen51 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=16804
                arg52 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=16805
                nen52 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=16806
                //arg53 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=16807
                //nen53 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=16808
                F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);  //Z=16809
                F52as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);  //Z=16810
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);  //Z=16811
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);  //Z=16812
                F52as2z0 = F52as20z*(F52as21+F52as22);  //Z=16813

                F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);  //Z=16815
                arg54 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=16816
                nen54 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=16817
                //arg55 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=16818
                //nen55 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=16819
                arg56 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=16820
                nen56 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=16821
                F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);  //Z=16822
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);  //Z=16823
                F52as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);  //Z=16824
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);  //Z=16825
                F52as3z0 = F52as30z*(F52as24+F52as26);  //Z=16826

                F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);  //Z=16828
                arg57 = (zr+2*c+1)*atan(2.0*(x1z-x2z));  //Z=16829
                nen57 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);  //Z=16830
                arg58 = (zr+2*c+1)*atan(2.0*(x1z+x2z));  //Z=16831
                nen58 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);  //Z=16832
                arg59 = (zr+2*c)*atan(2.0*(x1z-x2z));  //Z=16833
                nen59 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);  //Z=16834
                arg510 = (zr+2*c)*atan(2.0*(x1z+x2z));  //Z=16835
                nen510 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);  //Z=16836
                F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);  //Z=16837
                F52as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=16838
                F52as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=16839
                F52as4z = F52as40z*(F52as27+F52as28+F52as29);  //Z=16840
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;  //Z=16841
                F52asz0 = F52as1z0+F52as2z0+F52as3z0+F52as4z;  //Z=16842
                F52 = F52asz0;  //Z=16843
            }/*4*/  //Z=16844

            /* ** term #6 asymptote ** */  //Z=16846
            if ( xradp>=lim6 )
            {/*4*/  //Z=16847
                F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);  //Z=16848
                //F62as1sumz = pva0;  //Z=16849
                //F62as1z = F62as10z*F62as1sumz;  //Z=16850
                F62as1z0 = F62as10z*pza;  //Z=16851

                F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=16853
                arg61 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=16854
                nen61 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=16855
                arg62 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=16856
                nen62 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=16857
                //arg63 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=16858
                //nen63 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=16859
                F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=16860
                F62as22 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=16861
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);  //Z=16862
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);  //Z=16863
                F62as2z0 = F62as20z*(F62as21+F62as22);  //Z=16864

                F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=16866
                F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=16867
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);  //Z=16868
                F62as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=16869
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);  //Z=16870
                F62as3z0 = F62as30z*(F62as24+F62as26);  //Z=16871

                F62as40z = preg3*preg3*pow(x1z*x1z,c);  //Z=16873
                arg64 = (zr+2*c+1)*atan(4.0*x1z);  //Z=16874
                nen64 = pow(1.0+16*x1z*x1z,(zr+2*c+1)/2.0);  //Z=16875
                arg65 = (zr+2*c)*atan(4.0*x1z);  //Z=16876
                nen65 = pow(1.0+16*x1z*x1z,(zr+2*c)/2.0);  //Z=16877
                F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);  //Z=16878
                F62as28 = (1/2.0)*e0*e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=16879
                F62as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=16880
                F62as4z = F62as40z*(F62as27+F62as28+F62as29);  //Z=16881
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;  //Z=16882
                F62asz0 = F62as1z0+F62as2z0+F62as3z0+F62as4z;  //Z=16883
                F62 = F62asz0;  //Z=16884
            }/*4*/  //Z=16885

            /*formpq:=*/ return pql*(cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;  //Z=16887
            /* formpq:=pql*(cc2*F22)/vv3;  //Z=16888 */


            /* formpq:=pqcoreshellin(1.0,rho,p1,1.0,0.001,alfa,radiusm,2,sigmar,q);  //Z=16891 */

        }/*3*/ /*  of inhomogeneous core/shell  */  //Z=16893


        /*  myelin cylinder  */  //Z=16896
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=16897

            /*  cylinder parameters  */  //Z=16899
            v = -3/2.0;  //Z=16900
            e0 = 1;  //Z=16901
            e1 = -9/16.0;  //Z=16902
            preg1 = 1/sqrt(M_PI);  //Z=16903
            pz2v = 1/(zr*(zr-1)*(zr-2));  //Z=16904
            pz2v1 = pz2v/(zr-3);  //Z=16905
            pz2v2 = pz2v1/(zr-4);  //Z=16906
            lim = 18*exp(-5*params.sigma);  //Z=16907
            lim1 = lim*1.2;  //Z=16908
            rad = params.CR->myarray[1];  //Z=16909
            inmax = round(params.CR->myarray[14]);  //Z=16910
            vvm = params.CR->myarray[15];  //Z=16911
            rmax = params.CR->myarray[16];  //Z=16912
            xmax = q*rmax;  //Z=16913

            if ( xmax<(lim1) )
            {/*4*/  //Z=16915
                /* fkv[0]:=1;  //Z=16916 */
                qqn[0] = 1.0;  //Z=16917
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=16918
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=16919
                    /* fkv[nser]:=fkv[nser-1]*nser;  //Z=16920 */
                }/*5*/  //Z=16921

                F12sum = 0.0;  //Z=16923
                for ( ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=16924
                    for ( jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=16925
                        F12sez = 1.0;  //Z=16926
                        oldF12sez = 1.0;  //Z=16927
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16928
                            pqsum = 0;  //Z=16929
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16930
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]);  //Z=16931 */
                                pqsum = pqsum+pow(params.CR->carr7p[ii],2*mser)*pow(params.CR->carr7p[jj],2*(nser-mser))/(params.CR->carr6p[mser]*params.CR->carr6p[nser-mser]);  //Z=16932

                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16934 */
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx];  //Z=16935 */
                            }/*8*/  //Z=16936
                            F12sez = F12sez+params.CR->carr4p[nser]*qqn[nser]*pqsum;  //Z=16937
                            delser = fabs((F12sez-oldF12sez)/F12sez);  //Z=16938
                            if ( delser<0.0001 ) break; /* goto 101; */  //Z=16939
                            oldF12sez = F12sez;  //Z=16940
                        }/*7*/  //Z=16941
                        /*101:*/  //Z=16942
                        F12sum = F12sum+params.CR->carr5p[ii]*params.CR->carr5p[jj]*F12sez;  //Z=16943
                    }/*6*/  //Z=16944
                }/*5*/  //Z=16945
                F12ser = F12sum/vvm;  //Z=16946
                F12 = F12ser;  //Z=16947
            }/*4*/  //Z=16948
            else
            {/*4*/  //Z=16949
                xrz = q*rad/(zr+1);  //Z=16950
                arg = (zr+2*v+1)*atan(2.0*xrz);  //Z=16951
                nen = pow(1.0+4*xrz*xrz,(zr+2*v+1)/2.0);  //Z=16952
                arg1 = (zr+2*v)*atan(2.0*xrz);  //Z=16953
                nen1 = pow(1.0+4*xrz*xrz,(zr+2*v)/2.0);  //Z=16954
                arg2 = (zr+2*v-1)*atan(2.0*xrz);  //Z=16955
                nen2 = pow(1.0+4*xrz*xrz,(zr+2*v-1)/2.0);  //Z=16956

                F12asz = 0.0;  //Z=16958
                for ( ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=16959
                    a1m = params.CR->carr5p[ii]*pow(params.CR->carr7p[ii],v);   /*  carr7p[ii]:=pp[ii];  //Z=16960 */
                    for ( jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=16961
                        a2m = params.CR->carr5p[jj]*pow(params.CR->carr7p[jj],v);  //Z=16962
                        xijm = (params.CR->carr3p[ii]-params.CR->carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii];  //Z=16963 */
                        arglmz = (zr+1)*atan(xijm);  //Z=16964
                        nenlmz = pow(1.0+xijm*xijm,(zr+1)/2.0);  //Z=16965
                        xijp = (params.CR->carr3p[ii]+params.CR->carr3p[jj])*q/(zr+1);  //Z=16966
                        arglpz = (zr+1)*atan(xijp);  //Z=16967
                        nenlpz = pow(1.0+xijp*xijp,(zr+1)/2.0);  //Z=16968
                        /* F12as1z:=e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(pi*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(pi*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=16969 */
                        F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(0-(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=16970
                        /* F12as2z:=e0*e1*(1/(carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(pi*(2*v-1)/2)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(pi*(2*v-1)/2)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=16971 */
                        F12as2z = e0*e1*(1/(params.CR->carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(1*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-0)/(nen1*nenlpz));  //Z=16972
                        /* F12as3z:=e1*e0*(1/(carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(pi*(2*v-1)/2)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(pi*(2*v-1)/2)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=16973 */
                        F12as3z = e1*e0*(1/(params.CR->carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(1*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-0)/(nen1*nenlpz));  //Z=16974
                        /* F12as4z:=e1*e1*(1/(carr7p[ii]*carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(pi*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(pi*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=16975 */
                        F12as4z = e1*e1*(1/(params.CR->carr7p[ii]*params.CR->carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(0+1*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=16976

                        F12asz = F12asz+a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);  //Z=16978
                    }/*6*/  //Z=16979
                }/*5*/  //Z=16980
                F12asy = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;  //Z=16981
                F12 = F12asy;  //Z=16982
            }/*4*/  //Z=16983
            /*formpq:=*/ return pql*F12;  //Z=16984
            /* formpq:=pql;  //Z=16985 */

            /* formpq:=pql*polyliposome(llipt,radius,lliph,lin,lout,nom,sigmar,sigmal,phiax,philiph,philipt,phiin,phiout,2,q);  //Z=16987 */
            /* formpq:=pql*polyliposome(2.0,200,1.0,3.5,3.5,1,sigmar,sigmal,0.001,-0.55,-0.7,0.001,0.001,2,q);  //Z=16988 */
            /* formpq:=pql;  //Z=16989 */
        }/*3*/ /*  of myelin  */  //Z=16990

    }/*2*/ /*  of cylinder  */  //Z=16992



    /* ****** */  //Z=16996
    /*  disk  */  //Z=16997
    /* ****** */  //Z=16998
    if ( params.part==2 )
    {/*2*/  //Z=16999

        /* ** longitudinal part ** */  //Z=17001
        /* ** isotropic ** */  //Z=17002
        if ( ordis==7 )
        {/*3*/  //Z=17003
            if ( q<(0.5*params.limq1) )
            {/*4*/  //Z=17004
                pqsum = 1.0;  //Z=17005
                oldpqsum = 0.0;  //Z=17006
                qqn[0] = 1.0;  //Z=17007
                for ( nser=1; nser<=80; nser++ )
                {/*5*/  //Z=17008
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=17009
                    pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser];  //Z=17010
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17011
                    if ( delser<0.0001 ) break; /* goto 70; */  //Z=17012
                    oldpqsum = pqsum;  //Z=17013
                }/*5*/  //Z=17014
                /*70:*/  //Z=17015
                pql = pqsum;  //Z=17016
            }/*4*/  //Z=17017
            else
            {/*4*/  //Z=17018
                arglq = q*params.length/(zl+1);  //Z=17019
                pql = (2/(zl*(zl-1)))*pow(arglq,-2);  //Z=17020
            }/*4*/  //Z=17021
        }/*3*/  /*  of isotropic  */  //Z=17022

        /*  perfect  */  //Z=17024
        if ( ordis==6 )
        {/*3*/  //Z=17025
            if ( limql<0.7*params.limq1 )
            {/*4*/  //Z=17026

                pqsum = 1.0;  //Z=17028
                oldpqsum = 0.0;  //Z=17029
                qqn[0] = 1.0;  //Z=17030
                if ( params.orcase==1 )
                {/*5*/  //Z=17031
                    argq = qxs+qys;  //Z=17032
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17033
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17034
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser]*pow(1.0-argq*argq,nser);  //Z=17035
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17036
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=17037
                        oldpqsum = pqsum;  //Z=17038
                    }/*6*/  //Z=17039
                }/*5*/  //Z=17040
                if ( params.orcase==2 )
                {/*5*/  //Z=17041
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17042
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17043
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser]*pow(1.0-qxs*qxs,nser);  //Z=17044
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17045
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=17046
                        oldpqsum = pqsum;  //Z=17047
                    }/*6*/  //Z=17048
                }/*5*/  //Z=17049
                if ( params.orcase==3 )
                {/*5*/  //Z=17050
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17051
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17052
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser]*pow(1.0-qys*qys,nser);  //Z=17053
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17054
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=17055
                        oldpqsum = pqsum;  //Z=17056
                    }/*6*/  //Z=17057
                }/*5*/  //Z=17058
                if ( params.orcase==4 )
                {/*5*/  //Z=17059
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17060
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17061
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser];  //Z=17062
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17063
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=17064
                        oldpqsum = pqsum;  //Z=17065
                    }/*6*/  //Z=17066
                }/*5*/  //Z=17067
                /*76:*/  //Z=17068
                pql = pqsum;  //Z=17069
            }/*4*/  //Z=17070
            else
            {/*4*/  //Z=17071
                if ( params.orcase==1 )
                {/*5*/  //Z=17072
                    qnarg = qxs+qys;  //Z=17073
                    arglq = sqrt(1.0-qnarg*qnarg)*q*params.length/(zl+1)+eps9;  //Z=17074
                }/*5*/  //Z=17075
                if ( params.orcase==2 ) arglq = sqrt(1.0-qxs*qxs)*q*params.length/(zl+1)+eps9;  //Z=17076
                if ( params.orcase==3 ) arglq = sqrt(1.0-qys*qys)*q*params.length/(zl+1)+eps9;  //Z=17077
                if ( params.orcase==4 ) arglq = q*params.length/(zl+1)+eps9;  //Z=17078

                pqr1 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3);  //Z=17080
                pqr2 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3)*sin((zl-2)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-2)/2.0);  //Z=17081
                pqr3 = (1/(zl*(zl-1)*(zl-2)*(zl-3)))*pow(arglq,-4)*cos((zl-3)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-3)/2.0);  //Z=17082
                pql = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=17083
            }/*4*/  //Z=17084
        }/*3*/   /*  of perfect  */  //Z=17085

        /*  orientational distribution  */  //Z=17087
        if ( ordis==0 )
        {/*3*/       /*  general  */  //Z=17088
            if ( params.orcase==1 )
            {/*4*/  //Z=17089
                if ( limql<(0.3*params.limq1) )
                {/*5*/  //Z=17090
                    pqsum = 1.0;  //Z=17091
                    oldpqsum = 0.0;  //Z=17092
                    qqn[0] = 1.0;  //Z=17093
                    qxn[0] = 1.0;  //Z=17094
                    qyn[0] = 1.0;  //Z=17095

                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17097
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17098
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17099
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17100

                        binsum = 0.0;  //Z=17102
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=17103
                            binsum1 = 0.0;  //Z=17104
                            for ( lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=17105
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=17106 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=17107 */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=17108
                            }/*8*/  //Z=17109
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17110 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=17111 */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;  //Z=17112
                        }/*7*/  //Z=17113
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser]*binsum;  //Z=17114
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17115
                        if ( delser<0.0001 ) break; /* goto 77; */  //Z=17116
                        oldpqsum = pqsum;  //Z=17117
                    }/*6*/  //Z=17118
                    /*77:*/  //Z=17119
                    pql = pqsum;  //Z=17120
                }/*5*/  //Z=17121
                else
                {/*5*/  //Z=17122
                    /*  disk: length = disk radius  */  //Z=17123
                    /*  always use Bessel function approximation  */  //Z=17124
                    qrombdeltac(params.p1,sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17125
                    pql = pql/params.norm;  //Z=17126
                }/*5*/  //Z=17127
            }/*4*/  //Z=17128

            if ( params.orcase==2 )
            {/*4*/   /*  x-axis  */  //Z=17130
                if ( limql<(0.9*params.limq1) )
                {/*5*/  //Z=17131
                    pqsum = 1.0;  //Z=17132
                    oldpqsum = 0.0;  //Z=17133
                    qqn[0] = 1.0;  //Z=17134
                    qxn[0] = 1.0;  //Z=17135
                    qyn[0] = 1.0;  //Z=17136

                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17138
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17139
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17140
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17141

                        binsum = 0.0;  //Z=17143
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=17144
                            binsum1 = 0.0;  //Z=17145
                            for ( lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=17146
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=17147 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=17148 */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=17149
                            }/*8*/  //Z=17150
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17151 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=17152 */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;  //Z=17153
                        }/*7*/  //Z=17154
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser]*binsum;  //Z=17155
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17156
                        if ( delser<0.0001 ) break; /* goto 78; */  //Z=17157
                        oldpqsum = pqsum;  //Z=17158
                    }/*6*/  //Z=17159
                    /*78:*/  //Z=17160
                    pql = pqsum;  //Z=17161
                }/*5*/  //Z=17162
                else
                {/*5*/  //Z=17163
                    /*  disk: length = disk radius  */  //Z=17164
                    /*  always use Bessel function approximation  */  //Z=17165
                    qrombdeltac(params.p1,sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17166
                    pql = pql/params.norm;  //Z=17167
                    /* pql:=pq/1e-5;  //Z=17168 */
                    /* pql:=0.5;  //Z=17169 */
                }/*5*/  //Z=17170
            }/*4*/  //Z=17171

            if ( params.orcase==3 )
            {/*4*/     /*  y-axis  */  //Z=17173
                if ( limql<(0.9*params.limq1) )
                {/*5*/  //Z=17174
                    pqsum = 1.0;  //Z=17175
                    oldpqsum = 0.0;  //Z=17176
                    qqn[0] = 1.0;  //Z=17177
                    qxn[0] = 1.0;  //Z=17178
                    qyn[0] = 1.0;  //Z=17179

                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17181
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17182
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17183
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17184

                        binsum = 0.0;  //Z=17186
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=17187
                            binsum1 = 0.0;  //Z=17188
                            for ( lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=17189
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=17190 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=17191 */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=17192
                            }/*8*/  //Z=17193
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17194 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=17195 */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;  //Z=17196
                        }/*7*/  //Z=17197
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser]*binsum;  //Z=17198
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17199
                        if ( delser<0.0001 ) break; /* goto 79; */  //Z=17200
                        oldpqsum = pqsum;  //Z=17201
                    }/*6*/  //Z=17202
                    /*79:*/  //Z=17203
                    pql = pqsum;  //Z=17204
                }/*5*/  //Z=17205
                else
                {/*5*/  //Z=17206
                    /*  disk: length = disk radius  */  //Z=17207
                    /*  always use Bessel function approximation  */  //Z=17208
                    qrombdeltac(params.p1,sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17209
                    pql = pql/params.norm;  //Z=17210
                }/*5*/  //Z=17211
            }/*4*/  //Z=17212

            if ( params.orcase==4 )
            {/*4*/  //Z=17214
                if ( limql<(0.7*params.limq1) )
                {/*5*/  //Z=17215
                    pqsum = 1.0;  //Z=17216
                    oldpqsum = 0.0;  //Z=17217
                    qqn[0] = 1.0;  //Z=17218
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17219
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=17220
                        pqsum = pqsum+params.CR->carr1p[nser]*qqn[nser];  //Z=17221
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17222
                        if ( delser<0.0001 ) break; /* goto 80; */  //Z=17223
                        oldpqsum = pqsum;  //Z=17224
                    }/*6*/  //Z=17225
                    /*80:*/  //Z=17226
                    pql = pqsum;  //Z=17227
                }/*5*/  //Z=17228
                else
                {/*5*/  //Z=17229
                    /*  disk: length = disk radius  */  //Z=17230
                    /*  always use Bessel function approximation  */  //Z=17231
                    qrombdeltac(params.p1,sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17232
                    pql = pql/params.norm;  //Z=17233
                }/*5*/  //Z=17234
            }/*4*/  //Z=17235
        }/*3*/   /*  of orientational distribution  */  //Z=17236


        /*  transverse part  */  //Z=17239
        /*  disk: radius = disk thickness/2  */  //Z=17240
        /*  homogeneous disk  */  //Z=17241
        if ( params.cs==0 )
        {/*3*/  //Z=17242
            if ( q<(0.3*params.limq4) )
            {/*4*/  //Z=17243
                pqsum = 1.0;  //Z=17244
                oldpqsum = 0.0;  //Z=17245
                qqn[0] = 1.0;  //Z=17246
                for ( nser=1; nser<=100; nser++ )
                {/*5*/  //Z=17247
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=17248
                    pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=17249
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17250
                    if ( delser<0.0001 ) break; /* goto 71; */  //Z=17251
                    oldpqsum = pqsum;  //Z=17252
                }/*5*/  //Z=17253
                /*71:*/  //Z=17254
                pqr = pqsum;  //Z=17255
            }/*4*/  //Z=17256
            else
            {/*4*/  //Z=17257
                argpq = q*params.radius/(zr+1);  //Z=17258
                pqr = (1/(2.0*zr*(zr-1)))*pow(argpq,-2)*(1-cos((zr-1)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-1)/2.0));  //Z=17259
            }/*4*/  //Z=17260
            /*formpq:=*/ return pql*pqr;  //Z=17261
            /* formpq:=pql;  //Z=17262 */
        }/*3*/ /*  of homogeneous  */  //Z=17263

        /*  core/shell disk  */  //Z=17265
        if ( params.cs==1 )
        {/*3*/  //Z=17266
            ccc1 = sqr(1-params.rho)*pow(params.p1,2);  //Z=17267
            ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,1);  //Z=17268
            ccc3 = params.rho*params.rho;  //Z=17269
            vv3 = sqr((1-params.rho)*pow(params.p1,1)+params.rho);  //Z=17270

            argq = q*radiusm/(zz+1);  //Z=17272
            argpq = q*params.radius/(zz+1);  //Z=17273

            /*  F121 disk  */  //Z=17275
            if ( q<(0.8*params.limq4) )
            {/*4*/  //Z=17276
                /* ** series expansion ** */  //Z=17277
                pqsum = 1.0;  //Z=17278
                oldpqsum = 0.0;  //Z=17279
                qqn[0] = 1.0;  //Z=17280
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=17281
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=17282
                    pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=17283
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17284
                    if ( delser<0.0001 ) break; /* goto 72; */  //Z=17285
                    oldpqsum = pqsum;  //Z=17286
                }/*5*/  //Z=17287
                /*72:*/  //Z=17288
                F121 = ccc1*pqsum/vv3;  //Z=17289
            }/*4*/  //Z=17290
            else
            {/*4*/  //Z=17291
                pqr = (1/(2.0*zr*(zr-1)))*pow(argpq,-2)*(1-cos((zr-1)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-1)/2.0));  //Z=17292
                F121 = ccc1*pqr/vv3;  //Z=17293
            }/*4*/  //Z=17294

            /*  F122 disk  */  //Z=17296
            if ( q<(2.0*params.limq5) )
            {/*4*/  //Z=17297
                /* ** series expansion ** */  //Z=17298
                pqsum = 1.0;  //Z=17299
                oldpqsum = 0.0;  //Z=17300
                qqn[0] = 1.0;  //Z=17301
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=17302
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=17303
                    pqsum = pqsum+params.CR->carr5p[nser]*qqn[nser];  //Z=17304
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17305
                    if ( delser<0.0001 ) break; /* goto 73; */  //Z=17306
                    oldpqsum = pqsum;  //Z=17307
                }/*5*/  //Z=17308
                /*73:*/  //Z=17309
                F122 = ccc2*pqsum/vv3;  //Z=17310
            }/*4*/  //Z=17311
            else
            {/*4*/  //Z=17312
                argbm = (zr+1)*atan(argpq-argq);  //Z=17313
                nenbm = pow(1.0+4*sqr(argpq-argq),(zr+1)/2.0);  //Z=17314
                argbp = (zr+1)*atan(argpq+argq);  //Z=17315
                nenbp = pow(1.0+4*sqr(argpq+argq),(zr+1)/2.0);  //Z=17316

                pqr = (1/(2.0*zr*(zr-1)*(zr-2)*argpq*argq))*(cos(argbm)/nenbm-cos(argbp)/nenbp);  //Z=17318
                F122 = ccc2*pqr/vv3;  //Z=17319
            }/*4*/  //Z=17320

            /*  F123 disk  */  //Z=17322
            if ( q<(0.3*params.limq6) )
            {/*4*/  //Z=17323
                /* ** series expansion ** */  //Z=17324
                pqsum = 1.0;  //Z=17325
                oldpqsum = 0.0;  //Z=17326
                qqn[0] = 1.0;  //Z=17327
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=17328
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=17329
                    pqsum = pqsum+params.CR->carr6p[nser]*qqn[nser];  //Z=17330
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17331
                    if ( delser<0.0001 ) break; /* goto 74; */  //Z=17332
                    oldpqsum = pqsum;  //Z=17333
                }/*5*/  //Z=17334
                /*74:*/  //Z=17335
                F123 = ccc3*pqsum/vv3;  //Z=17336
            }/*4*/  //Z=17337
            else
            {/*4*/  //Z=17338
                pqr = (1/(2.0*zr*(zr-1)))*pow(argq,-2)*(1-cos((zr-1)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-1)/2.0));  //Z=17339
                F123 = ccc3*pqr/vv3;  //Z=17340
                /*  add more terms, if necessary  */  //Z=17341
            }/*4*/  //Z=17342
            /*formpq:=*/ return pql*(F121+F122+F123);  //Z=17343
            /* formpq:=F122;  //Z=17344 */
        }/*3*/ /*  of core/shell-disk  */  //Z=17345

        /*  inhomogeneous core/shell disk  */  //Z=17347
        if ( params.cs==2 )
        {/*3*/  //Z=17348

            dim = 1;  //Z=17350
            delc = 0.0001;  //Z=17351
            xrad = q*radiusm;  //Z=17352
            xradp = q*params.radius;  //Z=17353
            x1z = q*params.radius/(2.0*(zr+1));  //Z=17354
            x12z = x1z*x1z;  //Z=17355
            x2z = q*radiusm/(2.0*(zr+1));  //Z=17356
            x22z = x2z*x2z;  //Z=17357

            lim = 18*exp(-5*params.sigma);  //Z=17359
            lim1 = lim;  //Z=17360
            lim2 = lim*0.7;  //Z=17361
            lim3 = lim;  //Z=17362
            lim4 = lim;  //Z=17363
            lim5 = lim*0.7;  //Z=17364
            lim6 = lim*1.2;  //Z=17365

            a1 = (dim-params.alphash1)/2.0;  //Z=17367
            b1 = dim/2.0;  //Z=17368
            b2 = (dim+2-params.alphash1)/2.0;  //Z=17369
            b1s = (dim+2)/2.0;  //Z=17370
            v = -b1s+1/2.0;  //Z=17371
            c = a1-b1-b2+1/2.0;  //Z=17372
            d0 = 1;  //Z=17373
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=17374
            e0 = 1.0;  //Z=17375
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=17376
            ee0 = 1.0;  //Z=17377
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=17378

            gb1s = sqrt(M_PI)/2.0;  //Z=17380
            pz2v = 1/(zr*(zr-1));  //Z=17381
            pz2v1 = pz2v/(zr-2);  //Z=17382
            pz2v2 = pz2v1/(zr-3);  //Z=17383

            gz1 = gamma(zr+1);  //Z=17385
            preg1 = gb1s/sqrt(M_PI);  //Z=17386
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=17387
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=17388
            pzvc = gamma(zr+1+v+c)/gz1;  //Z=17389
            pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=17390
            pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=17391
            pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=17392
            pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=17393
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=17394
            pzc = gamma(zr+1+2*c)/gz1;  //Z=17395
            pzc1 = gamma(zr+1+2*c-1)/gz1;  //Z=17396
            pza = gamma(zr+1-4*a1)/gz1;  //Z=17397
            pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=17398
            pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=17399
            //dnv0 = 1;  //Z=17400
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=17401
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=17402
            //pva0 = gamma(zr+1-4*a1)/gz1;  //Z=17403

            cc1 = 1/(dim*dim);  //Z=17405
            cc2 = 2*params.rho/(dim*(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=17406
            cc3 = -2*params.rho/(dim*(dim-params.alphash1));  //Z=17407
            cc4 = sqr(params.rho)/(sqr(dim-params.alphash1)*pow(sqr(params.p1),dim-params.alphash1));  //Z=17408
            cc5 = -2*sqr(params.rho)/(sqr(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=17409
            cc6 = sqr(params.rho)/sqr(dim-params.alphash1);  //Z=17410
            vv3 = cc1+cc2+cc3+cc4+cc5+cc6;  //Z=17411

            /*  term #1 series  */  //Z=17413
            if ( (xradp)<lim1 )
            {/*4*/  //Z=17414
                //z12v[0] = 1;  //Z=17415
                //b1sv[0] = 1;  //Z=17416
                //fkv[0] = 1;  //Z=17417
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17418
                F12sez = 1.0;  //Z=17419
                oldF12sez = 0.0;  //Z=17420
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=17421
                    //qqn[n] = qqn[n-1]*q*q;  //Z=17422
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17423
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17424
                    //fkv[n] = fkv[n-1]*n;  //Z=17425
                    //sum12[n] = 0;  //Z=17426
                    /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);  //Z=17427 */
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n];  //Z=17428 */

                    F12sez += params.CR->carr4p[n]*qqnn; //qqn[n];  //Z=17430

                    del = fabs((F12sez-oldF12sez)/F12sez);  //Z=17432
                    if ( del<delc ) break; /* goto 221; */  //Z=17433
                    oldF12sez = F12sez;  //Z=17434
                }/*5*/  //Z=17435
                /*221:*/  //Z=17436
                F12 = F12sez;  //Z=17437
            }/*4*/  //Z=17438

            /*  term #2 series  */  //Z=17440
            if ( (xradp)<lim2 )
            {/*4*/  //Z=17441
                //z12v[0] = 1;  //Z=17442
                //a1v[0] = 1;  //Z=17443
                //b1v[0] = 1;  //Z=17444
                //b2v[0] = 1;  //Z=17445
                //b1sv[0] = 1;  //Z=17446
                //fkv[0] = 1;  //Z=17447
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17448
                F22sez = 1.0;  //Z=17449
                oldF22sez = 0.0;  //Z=17450
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=17451
                    //qqn[n] = qqn[n-1]*q*q;  //Z=17452
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17453
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17454
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17455
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17456
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17457
                    //fkv[n] = fkv[n-1]*n;  //Z=17458
                    //sum22[n] = 0;  //Z=17459
                    /* for m:=0 to n do sum22[n]:=sum22[n]+a1v[n-m]*power(p1*p1,m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17460 */
                    /* F22sez:=F22sez+power(-x22z,n)*z12v[n]*sum22[n];  //Z=17461 */

                    F22sez += params.CR->carr5p[n]*qqnn; //qqn[n];  //Z=17463

                    del = fabs((F22sez-oldF22sez)/F22sez);  //Z=17465
                    if ( del<delc ) break; /* goto 222; */  //Z=17466
                    oldF22sez = F22sez;  //Z=17467
                }/*5*/  //Z=17468
                /*222:*/  //Z=17469
                F22 = F22sez;  //Z=17470
            }/*4*/  //Z=17471

            /*  term #3 series  */  //Z=17473
            if ( (xradp)<lim3 )
            {/*4*/  //Z=17474
                //z12v[0] = 1;  //Z=17475
                //a1v[0] = 1;  //Z=17476
                //b1v[0] = 1;  //Z=17477
                //b2v[0] = 1;  //Z=17478
                //b1sv[0] = 1;  //Z=17479
                //fkv[0] = 1;  //Z=17480
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17481
                F32sez = 1.0;  //Z=17482
                oldF32sez = 0.0;  //Z=17483
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=17484
                    //qqn[n] = qqn[n-1]*q*q;  //Z=17485
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17486
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17487
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17488
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17489
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17490
                    //fkv[n] = fkv[n-1]*n;  //Z=17491
                    //sum32[n] = 0;  //Z=17492
                    /* for m:=0 to n do sum32[n]:=sum32[n]+a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17493 */
                    /* F32sez:=F32sez+power(-x12z,n)*z12v[n]*sum32[n];  //Z=17494 */

                    F32sez += params.CR->carr6p[n]*qqnn; //qqn[n];  //Z=17496

                    del = fabs((F32sez-oldF32sez)/F32sez);  //Z=17498
                    if ( del<delc ) break; /* goto 223; */  //Z=17499
                    oldF32sez = F32sez;  //Z=17500
                }/*5*/  //Z=17501
                /*223:*/  //Z=17502
                F32 = F32sez;  //Z=17503
            }/*4*/  //Z=17504

            /*  term #4 series  */  //Z=17506
            if ( (xradp)<lim4 )
            {/*4*/  //Z=17507
                //z12v[0] = 1;  //Z=17508
                //a1v[0] = 1;  //Z=17509
                //b1v[0] = 1;  //Z=17510
                //b2v[0] = 1;  //Z=17511
                //b1sv[0] = 1;  //Z=17512
                //fkv[0] = 1;  //Z=17513
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17514
                F42sez = 1.0;  //Z=17515
                oldF42sez = 0.0;  //Z=17516
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=17517
                    //qqn[n] = qqn[n-1]*q*q;  //Z=17518
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17519
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17520
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17521
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17522
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17523
                    //fkv[n] = fkv[n-1]*n;  //Z=17524
                    //sum42[n] = 0;  //Z=17525
                    /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17526 */
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n];  //Z=17527 */

                    F42sez += params.CR->carr7p[n]*qqnn; //qqn[n];  //Z=17529

                    del = fabs((F42sez-oldF42sez)/F42sez);  //Z=17531
                    if ( del<delc ) break; /* goto 224; */  //Z=17532
                    oldF42sez = F42sez;  //Z=17533
                }/*5*/  //Z=17534
                /*224:*/  //Z=17535
                F42 = F42sez;  //Z=17536
            }/*4*/  //Z=17537

            /*  term #5 series  */  //Z=17539
            if ( (xradp)<lim5 )
            {/*4*/  //Z=17540
                //z12v[0] = 1;  //Z=17541
                //a1v[0] = 1;  //Z=17542
                //b1v[0] = 1;  //Z=17543
                //b2v[0] = 1;  //Z=17544
                //b1sv[0] = 1;  //Z=17545
                //fkv[0] = 1;  //Z=17546
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17547
                F52sez = 1.0;  //Z=17548
                oldF52sez = 0.0;  //Z=17549
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=17550
                    //qqn[n] = qqn[n-1]*q*q;  //Z=17551
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17552
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17553
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17554
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17555
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17556
                    //fkv[n] = fkv[n-1]*n;  //Z=17557
                    //sum52[n] = 0;  //Z=17558
                    /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17559 */
                    /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n];  //Z=17560 */

                    F52sez += params.CR->carr8p[n]*qqnn; //qqn[n];  //Z=17562

                    del = fabs((F52sez-oldF52sez)/F52sez);  //Z=17564
                    if ( del<delc ) break; /* goto 225; */  //Z=17565
                    oldF52sez = F52sez;  //Z=17566
                }/*5*/  //Z=17567
                /*225:*/  //Z=17568
                F52 = F52sez;  //Z=17569
            }/*4*/  //Z=17570

            /*  term #6 series  */  //Z=17572
            if ( (xradp)<lim6 )
            {/*4*/  //Z=17573
                //z12v[0] = 1;  //Z=17574
                //a1v[0] = 1;  //Z=17575
                //b1v[0] = 1;  //Z=17576
                //b2v[0] = 1;  //Z=17577
                //b1sv[0] = 1;  //Z=17578
                //fkv[0] = 1;  //Z=17579
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17580
                F62sez = 1.0;  //Z=17581
                oldF62sez = 0.0;  //Z=17582
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=17583
                    //qqn[n] = qqn[n-1]*q*q;  //Z=17584
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17585
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17586
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17587
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17588
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17589
                    //fkv[n] = fkv[n-1]*n;  //Z=17590
                    //sum62[n] = 0;  //Z=17591
                    /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17592 */
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n];  //Z=17593 */

                    F62sez += params.CR->carr9p[n]*qqnn; //qqn[n];  //Z=17595

                    del = fabs((F62sez-oldF62sez)/F62sez);  //Z=17597
                    if ( del<delc ) break; /* goto 226; */  //Z=17598
                    oldF62sez = F62sez;  //Z=17599
                }/*5*/  //Z=17600
                /*226:*/  //Z=17601
                F62 = F62sez;  //Z=17602
            }/*4*/  //Z=17603


            /* ** term #1 asymptote ** */  //Z=17606
            if ( xradp>=lim1 )
            {/*4*/  //Z=17607
                arg11 = (zr+2*v+1)*atan(4.0*x1z);  //Z=17608
                nen11 = pow(1.0+16*x1z*x1z,(zr+2*v+1)/2.0);  //Z=17609
                arg12 = (zr+2*v)*atan(4.0*x1z);  //Z=17610
                nen12 = pow(1.0+16*x1z*x1z,(zr+2*v)/2.0);  //Z=17611
                arg13 = (zr+2*v-1)*atan(4.0*x1z);  //Z=17612
                nen13 = pow(1.0+16*x1z*x1z,(zr+2*v-1)/2.0);  //Z=17613

                F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);  //Z=17615
                F12as2z = 2*ee0*ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);  //Z=17616
                F12as3z = ee1*ee1*(1/(4.0*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);  //Z=17617
                F12asz = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);  //Z=17618
                F12 = F12asz;  //Z=17619
            }/*4*/  //Z=17620

            /* ** term #2 asymptote ** */  //Z=17622
            if ( xradp>=lim2 )
            {/*4*/  //Z=17623
                //arg21 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=17624
                //nen21 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=17625
                //arg22 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=17626
                //nen22 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=17627
                //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);  //Z=17628
                //F22as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);  //Z=17629
                F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);  //Z=17630
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);  //Z=17631

                arg210 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=17633
                nen210 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=17634
                arg220 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=17635
                nen220 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=17636
                F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);  //Z=17637
                F22as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);  //Z=17638
                F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);  //Z=17639
                arg23 = (zr+v+c+1)*atan(2.0*(x1z-x2z));  //Z=17640
                nen23 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);  //Z=17641
                arg24 = (zr+v+c+1)*atan(2.0*(x1z+x2z));  //Z=17642
                nen24 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);  //Z=17643
                arg25 = (zr+v+c)*atan(2.0*(x1z-x2z));  //Z=17644
                nen25 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);  //Z=17645
                arg26 = (zr+v+c)*atan(2.0*(x1z+x2z));  //Z=17646
                nen26 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);  //Z=17647
                arg27 = (zr+v+c-1)*atan(2.0*(x1z-x2z));  //Z=17648
                nen27 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);  //Z=17649
                arg28 = (zr+v+c-1)*atan(2.0*(x1z+x2z));  //Z=17650
                nen28 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);  //Z=17651

                a22as21z = (1/2.0)*ee0*e0*pzvc;  //Z=17653
                F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);  //Z=17654
                a22as22z = (1/2.0)*ee0*e1*(1/(2.0*x2z))*pzvc1;  //Z=17655
                F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);  //Z=17656
                a22as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1;  //Z=17657
                F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);  //Z=17658
                a22as24z = (1/2.0)*ee1*e1*(1/(2.0*x1z))*(1/(2.0*x2z))*pzvc2;  //Z=17659
                F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);  //Z=17660
                F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);  //Z=17661
                F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);  //Z=17662
                //F22asz = F22as1z+F22as2z;  //Z=17663
                F22asz0 = F22as1z0+F22as2z;  //Z=17664
                F22 = F22asz0;  //Z=17665
            }/*4*/  //Z=17666

            /* ** term #3 asymptote ** */  //Z=17668
            if ( xradp>=lim3 )
            {/*4*/  //Z=17669
                //arg31 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=17670
                //nen31 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=17671
                //arg32 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=17672
                //nen32 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=17673
                //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);  //Z=17674
                //F32as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);  //Z=17675
                F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);  //Z=17676
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);  //Z=17677

                arg310 = (z+v-2*a1+1)*atan(2.0*x1z);  //Z=17679
                nen310 = pow(1.0+4*x1z*x1z,(z+v-2*a1+1)/2.0);  //Z=17680
                arg320 = (z+v-2*a1)*atan(2.0*x1z);  //Z=17681
                nen320 = pow(1.0+4*x1z*x1z,(z+v-2*a1)/2.0);  //Z=17682
                F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);  //Z=17683
                F32as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);  //Z=17684
                F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);  //Z=17685

                arg33 = (zr+v+c+1)*atan(4.0*x1z);  //Z=17687
                nen33 = pow(1.0+16*x1z*x1z,(zr+v+c+1)/2.0);  //Z=17688
                arg34 = (zr+v+c)*atan(4.0*x1z);  //Z=17689
                nen34 = pow(1.0+16*x1z*x1z,(zr+v+c)/2.0);  //Z=17690
                arg35 = (zr+v+c-1)*atan(4.0*x1z);  //Z=17691
                nen35 = pow(1.0+16*x1z*x1z,(zr+v+c-1)/2.0);  //Z=17692
                F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);  //Z=17693
                F32as22z = (1/2.0)*ee0*e1*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);  //Z=17694
                F32as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);  //Z=17695
                F32as24z = (1/2.0)*ee1*e1*(1/(4.0*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);  //Z=17696
                F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);  //Z=17697
                F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);  //Z=17698
                //F32asz = F32as1z+F32as2z;  //Z=17699
                F32asz0 = F32as1z0+F32as2z;  //Z=17700
                F32 = F32asz0;  //Z=17701
            }/*4*/  //Z=17702


            /* ** term #4 asymptote ** */  //Z=17705
            if ( xrad>=lim4 )
            {/*4*/  //Z=17706
                F42as10z = preg4*preg4*pow(x22z,-2*a1);  //Z=17707
                //F42as1sumz = pva0;  //Z=17708
                //F42as1z = F42as10z*F42as1sumz;  //Z=17709
                F42as1z0 = F42as10z*pza;  //Z=17710

                arg41 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=17712
                nen41 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=17713
                arg42 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=17714
                nen42 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=17715
                //arg43 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=17716
                //nen43 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=17717
                F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=17718
                F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=17719
                F42as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);  //Z=17720
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);  //Z=17721
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);  //Z=17722
                F42as2z0 = F42as20z*(F42as21+F42as22);  //Z=17723

                F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=17725
                F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=17726
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);  //Z=17727
                F42as26 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);  //Z=17728
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);  //Z=17729
                F42as3z0 = F42as30z*(F42as24+F42as26);  //Z=17730

                F42as40z = preg3*preg3*pow(x2z*x2z,c);  //Z=17732
                arg44 = (zr+2*c+1)*atan(4.0*x2z);  //Z=17733
                nen44 = pow(1.0+16*x2z*x2z,(zr+2*c+1)/2.0);  //Z=17734
                arg45 = (zr+2*c)*atan(4.0*x2z);  //Z=17735
                nen45 = pow(1.0+16*x2z*x2z,(zr+2*c)/2.0);  //Z=17736
                F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);  //Z=17737
                F42as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=17738
                F42as29 = (1/2.0)*e1*e0*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=17739
                F42as4z = F42as40z*(F42as27+F42as28+F42as29);  //Z=17740
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;  //Z=17741
                F42asz0 = F42as1z0+F42as2z0+F42as3z0+F42as4z;  //Z=17742
                F42 = F42asz0;  //Z=17743
            }/*4*/  //Z=17744


            /* ** term #5 asymptote ** */  //Z=17747
            if ( xradp>=lim5 )
            {/*4*/  //Z=17748
                F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);  //Z=17749
                //F52as1sumz = pva0;  //Z=17750
                //F52as1z = F52as10z*F52as1sumz;  //Z=17751
                F52as1z0 = F52as10z*pza;  //Z=17752

                F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);  //Z=17754
                arg51 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=17755
                nen51 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=17756
                arg52 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=17757
                nen52 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=17758
                //arg53 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=17759
                //nen53 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=17760
                F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);  //Z=17761
                F52as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);  //Z=17762
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);  //Z=17763
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);  //Z=17764
                F52as2z0 = F52as20z*(F52as21+F52as22);  //Z=17765

                F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);  //Z=17767
                arg54 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=17768
                nen54 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=17769
                //arg55 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=17770
                //nen55 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=17771
                arg56 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=17772
                nen56 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=17773
                F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);  //Z=17774
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);  //Z=17775
                F52as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);  //Z=17776
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);  //Z=17777
                F52as3z0 = F52as30z*(F52as24+F52as26);  //Z=17778

                F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);  //Z=17780
                arg57 = (zr+2*c+1)*atan(2.0*(x1z-x2z));  //Z=17781
                nen57 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);  //Z=17782
                arg58 = (zr+2*c+1)*atan(2.0*(x1z+x2z));  //Z=17783
                nen58 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);  //Z=17784
                arg59 = (zr+2*c)*atan(2.0*(x1z-x2z));  //Z=17785
                nen59 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);  //Z=17786
                arg510 = (zr+2*c)*atan(2.0*(x1z+x2z));  //Z=17787
                nen510 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);  //Z=17788
                F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);  //Z=17789
                F52as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=17790
                F52as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=17791
                F52as4z = F52as40z*(F52as27+F52as28+F52as29);  //Z=17792
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;  //Z=17793
                F52asz0 = F52as1z0+F52as2z0+F52as3z0+F52as4z;  //Z=17794
                F52 = F52asz0;  //Z=17795
            }/*4*/  //Z=17796

            /* ** term #6 asymptote ** */  //Z=17798
            if ( xradp>=lim6 )
            {/*4*/  //Z=17799
                F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);  //Z=17800
                //F62as1sumz = pva0;  //Z=17801
                //F62as1z = F62as10z*F62as1sumz;  //Z=17802
                F62as1z0 = F62as10z*pza;  //Z=17803

                F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=17805
                arg61 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=17806
                nen61 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=17807
                arg62 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=17808
                nen62 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=17809
                //arg63 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=17810
                //nen63 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=17811
                F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=17812
                F62as22 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=17813
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);  //Z=17814
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);  //Z=17815
                F62as2z0 = F62as20z*(F62as21+F62as22);  //Z=17816

                F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=17818
                F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=17819
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);  //Z=17820
                F62as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=17821
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);  //Z=17822
                F62as3z0 = F62as30z*(F62as24+F62as26);  //Z=17823

                F62as40z = preg3*preg3*pow(x1z*x1z,c);  //Z=17825
                arg64 = (zr+2*c+1)*atan(4.0*x1z);  //Z=17826
                nen64 = pow(1.0+16*x1z*x1z,(zr+2*c+1)/2.0);  //Z=17827
                arg65 = (zr+2*c)*atan(4.0*x1z);  //Z=17828
                nen65 = pow(1.0+16*x1z*x1z,(zr+2*c)/2.0);  //Z=17829
                F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);  //Z=17830
                F62as28 = (1/2.0)*e0*e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=17831
                F62as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=17832
                F62as4z = F62as40z*(F62as27+F62as28+F62as29);  //Z=17833
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;  //Z=17834
                F62asz0 = F62as1z0+F62as2z0+F62as3z0+F62as4z;  //Z=17835
                F62 = F62asz0;  //Z=17836
            }/*4*/  //Z=17837

            /*formpq:=*/ return pql*(cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;  //Z=17839
            /* formpq:=pql*(cc5*F52)/vv3;  //Z=17840 */


            /* formpq:=pql*(F121+F122+F123);  //Z=17843 */

        }/*3*/ /*  of inhomogeneous core/shell-disk  */  //Z=17845

        /*  myelin disk  */  //Z=17847
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=17848

            /*  disk parameters  */  //Z=17850
            v = -1;  //Z=17851
            e0 = 1;  //Z=17852
            e1 = 0;  //Z=17853
            preg1 = 1/2.0;  //Z=17854
            pz2v = 1/(zr*(zr-1));  //Z=17855
            pz2v1 = pz2v/(zr-2);  //Z=17856
            pz2v2 = pz2v1/(zr-3);  //Z=17857
            lim = 18*exp(-5*params.sigma);  //Z=17858
            lim1 = lim*1.2;  //Z=17859
            rad = params.CR->myarray[1];  //Z=17860
            inmax = round(params.CR->myarray[14]);  //Z=17861
            vvm = params.CR->myarray[15];  //Z=17862
            rmax = params.CR->myarray[16];  //Z=17863
            xmax = q*rmax;  //Z=17864

            if ( xmax<(lim1) )
            {/*4*/  //Z=17866
                /* fkv[0]:=1;  //Z=17867 */
                qqn[0] = 1.0;  //Z=17868
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=17869
                    qqn[nser] = qqn[nser-1]*q*q;  //Z=17870
                    /* fkv[nser]:=fkv[nser-1]*nser;  //Z=17871 */
                }/*5*/  //Z=17872

                F12sum = 0.0;  //Z=17874
                for ( ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=17875
                    for ( jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=17876
                        F12sez = 1.0;  //Z=17877
                        oldF12sez = 1.0;  //Z=17878
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=17879
                            pqsum = 0;  //Z=17880
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=17881
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]);  //Z=17882 */
                                pqsum = pqsum+pow(params.CR->carr7p[ii],2*mser)*pow(params.CR->carr7p[jj],2*(nser-mser))/(params.CR->carr6p[mser]*params.CR->carr6p[nser-mser]);  //Z=17883

                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17885 */
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx];  //Z=17886 */
                            }/*8*/  //Z=17887
                            F12sez = F12sez+params.CR->carr4p[nser]*qqn[nser]*pqsum;  //Z=17888
                            delser = fabs((F12sez-oldF12sez)/F12sez);  //Z=17889
                            if ( delser<0.0001 ) break; /* goto 251; */  //Z=17890
                            oldF12sez = F12sez;  //Z=17891
                        }/*7*/  //Z=17892
                        /*251:*/  //Z=17893
                        F12sum = F12sum+params.CR->carr5p[ii]*params.CR->carr5p[jj]*F12sez;  //Z=17894
                    }/*6*/  //Z=17895
                }/*5*/  //Z=17896
                F12ser = F12sum/vvm;  //Z=17897
                F12 = F12ser;  //Z=17898
            }/*4*/  //Z=17899
            else
            {/*4*/  //Z=17900
                xrz = q*rad/(zr+1);  //Z=17901
                arg = (zr+2*v+1)*atan(2.0*xrz);  //Z=17902
                nen = pow(1.0+4*xrz*xrz,(zr+2*v+1)/2.0);  //Z=17903
                arg1 = (zr+2*v)*atan(2.0*xrz);  //Z=17904
                nen1 = pow(1.0+4*xrz*xrz,(zr+2*v)/2.0);  //Z=17905
                arg2 = (zr+2*v-1)*atan(2.0*xrz);  //Z=17906
                nen2 = pow(1.0+4*xrz*xrz,(zr+2*v-1)/2.0);  //Z=17907

                F12asz = 0.0;  //Z=17909
                for ( ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=17910
                    a1m = params.CR->carr5p[ii]*pow(params.CR->carr7p[ii],v);   /*  carr7p[ii]:=pp[ii];  //Z=17911 */
                    for ( jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=17912
                        a2m = params.CR->carr5p[jj]*pow(params.CR->carr7p[jj],v);  //Z=17913
                        xijm = (params.CR->carr3p[ii]-params.CR->carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii];  //Z=17914 */
                        arglmz = (zr+1)*atan(xijm);  //Z=17915
                        nenlmz = pow(1.0+xijm*xijm,(zr+1)/2.0);  //Z=17916
                        xijp = (params.CR->carr3p[ii]+params.CR->carr3p[jj])*q/(zr+1);  //Z=17917
                        arglpz = (zr+1)*atan(xijp);  //Z=17918
                        nenlpz = pow(1.0+xijp*xijp,(zr+1)/2.0);  //Z=17919
                        F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(M_PI*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(M_PI*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=17920
                        F12as2z = e0*e1*(1/(params.CR->carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=17921
                        F12as3z = e1*e0*(1/(params.CR->carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=17922
                        F12as4z = e1*e1*(1/(params.CR->carr7p[ii]*params.CR->carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(M_PI*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(M_PI*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=17923

                        F12asz = F12asz+a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);  //Z=17925
                    }/*6*/  //Z=17926
                }/*5*/  //Z=17927
                F12asy = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;  //Z=17928
                F12 = F12asy;  //Z=17929
            }/*4*/  //Z=17930
            /*formpq:=*/ return pql*F12;  //Z=17931
            /* formpq:=pql;  //Z=17932 */

            /* formpq:=pql*polyliposome(llipt,radius,lliph,lin,lout,nom,sigmar,sigmal,phiax,philiph,philipt,phiin,phiout,2,q);  //Z=17934 */
            /* formpq:=pql*polyliposome(2.0,200,1.0,3.5,3.5,1,sigmar,sigmal,0.001,-0.55,-0.7,0.001,0.001,1,q);  //Z=17935 */
            /* formpq:=pql;  //Z=17936 */
        }/*3*/ /*  of myelin disk  */  //Z=17937

    }/*2*/ /*  of disk  */  //Z=17939


    /*  cube  */  //Z=17942
    if ( params.part==4 )
    {/*2*/  //Z=17943

        //if ( q<0.02 )
        //    qDebug() << "formpq (cube): ordis"<<ordis << "cs"<<cs << "q"<<q << "limq4"<<limq4;

        /*  homogeneous isotropic cube  */  //Z=17944
        if ( ordis==7 )
        {/*3*/  //Z=17945
            if ( params.cs==0 )
            {/*4*/  //Z=17946
                if ( q<0.7*params.limq4 )
                {/*5*/  //Z=17947
                    pqsum = 1.0;  //Z=17948
                    oldpqsum = 0.0;  //Z=17949
                    double qqn = 1.0;  //Z=17950
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=17951
                        qqn = qqn*q*q;  //Z=17952
                        pqsum = pqsum+params.CR->carr4p[nser]*qqn;  //Z=17953
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17954
                        if ( delser<0.0001 ) break; /* goto 81; */  //Z=17955
                        oldpqsum = pqsum;  //Z=17956
                    }/*6*/  //Z=17957
                    /*81:*/  //Z=17958
                    //qDebug() << "formpq A" << q << 0.7*limq4 << "=" << pqsum << nser;
                    /*formpq:=*/ return pqsum;  //Z=17959
                }/*5*/  //Z=17960
                else
                {
                    //if ( q < 0.75*limq4 )
                    //    qDebug() << "formpq B" << q << 0.7*limq4 << "=" << por/(q*q*q*q) << "por" << por;
                    /*formpq:=*/ return params.por/(q*q*q*q);  //Z=17961
                }
            }/*4*/ /*  of homogeneous isotropic cube */  //Z=17962

            /*  core/shell isotropic cube  */  //Z=17964
            if ( params.cs==1 )
            {/*4*/  //Z=17965
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=17966
            }/*4*/  //Z=17967
        }/*3*/  /*  of isotropic cube  */  //Z=17968

        /*  perfectly oriented cube  */  //Z=17970
        if ( ordis==6 )
        {/*3*/  //Z=17971
            /* if (orcase=1) then begin  //Z=17972 */
            if ( 1==1 )
            {/*4*/  //Z=17973

                if ( q<(1/params.radius) )
                {/*5*/  //Z=17975
                    pqsum = 1.0;  //Z=17976
                    oldpqsum = 0.0;  //Z=17977
                    qxn[0] = 1.0;  //Z=17978
                    qyn[0] = 1.0;  //Z=17979
                    for ( nser=1; nser<=80; nser++ )
                    {/*6*/  //Z=17980
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17981
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17982

                        binsum = 0.0;  //Z=17984
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=17985
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17986 */
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=17987 */
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=17988
                        }/*7*/  //Z=17989
                        pqsum = pqsum+binsum;  //Z=17990
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17991
                        if ( delser<0.0001 ) break; /* goto 84; */  //Z=17992
                        oldpqsum = pqsum;  //Z=17993
                    }/*6*/  //Z=17994
                    /*84:*/  //Z=17995
                    pql = pqsum;  //Z=17996
                }/*5*/  //Z=17997
                else
                {/*5*/  //Z=17998
                    argqx = qxs*params.radius/(zr+1)+eps9;  //Z=17999
                    argqy = qys*params.radius/(zr+1)+eps9;  //Z=18000
                    pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18001
                    pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18002
                    pql = pqrx*pqry;  //Z=18003
                }/*5*/  //Z=18004
                /*formpq:=*/ return pql;  //Z=18005
            }/*4*/  /*  of orcase=1  */  //Z=18006
        }/*3*/  /*  of perfect cube  */  //Z=18007
    }/*2*/  /*  of cube  */  //Z=18008


    /* ** biaxial ellipsoid ** */  //Z=18011
    if ( params.part==5 )
    {/*2*/  //Z=18012
        /*  homogeneous isotropic ellipsoid  */  //Z=18013
        if ( ordis==7 )
        {/*3*/  //Z=18014
            if ( params.cs==0 )
            {/*4*/  //Z=18015
                if ( q<0.8*params.limq4 )
                {/*5*/  //Z=18016
                    pqsum = 1.0;  //Z=18017
                    oldpqsum = 0.0;  //Z=18018
                    qqn[0] = 1.0;  //Z=18019
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18020
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=18021
                        pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=18022
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18023
                        if ( delser<0.0001 ) break; /* goto 260; */  //Z=18024
                        oldpqsum = pqsum;  //Z=18025
                    }/*6*/  //Z=18026
                    /*260:*/  //Z=18027
                    /*formpq:=*/ return pqsum;  //Z=18028
                }/*5*/  //Z=18029
                else
                {/*5*/  //Z=18030
                    if ( q>=1.5*params.limq4 )
                        pq = params.por/(q*q*q*q);  //Z=18031
                    else
                    {/*6*/  //Z=18032
                        //qrombchid(l,r,p1,sigma,alfa,dbeta,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2053
                        double qhkl = 1;
                        qrombchid(params.length,params.radius,params.p1,params.sigma,params.alphash1,epsi,params.polTheta,params.polPhi,
                                  qx,qy,qz,params.p11,params.p12,params.p13,params.p21,params.p22,
                                  params.p23,params.p31,params.p32,params.p33,qx,qy,0,qhkl,
                                  params.ax1.length(),params.ax2.length(),params.ax3.length(),
                                  params.ax1.x(),params.ax1.y(),params.ax1.z(),
                                  params.ax2.x(),params.ax2.y(),params.ax2.z(),
                                  params.ax3.x(),params.ax3.y(),params.ax3.z(),
                                  params.sig.x(),params.sig.y(),params.sig.z(),
                                  ordis,3,8,13,7,0,0,params.CR->carr1p,pql);  //Z=18033
                        pq = pql;  //Z=18034
                    }/*6*/  //Z=18035
                    /*formpq:=*/ return pq;  //Z=18036
                }/*5*/  //Z=18037
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18038

            /*  core/shell isotropic ellipsoid  */  //Z=18040
            if ( params.cs==1 )
            {/*4*/  //Z=18041
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18042
            }/*4*/  //Z=18043
        }/*3*/  /*  of isotropic ellipsoid  */  //Z=18044

        /*  perfectly oriented ellipsoid  */  //Z=18046
        if ( ordis==6 )
        {/*3*/  //Z=18047
            /* if (orcase=1) then begin  //Z=18048 */
            if ( 1==1 )
            {/*4*/  //Z=18049

                if ( sqrt(qx*qx*sqr(params.length)+qy*qy*sqr(params.radius)+eps9)<15 )
                {/*5*/  //Z=18051
                    qxn[0] = 1.0;  //Z=18052
                    qyn[0] = 1.0;  //Z=18053
                    for ( nser=1; nser<=81; nser++ )
                    {/*6*/  //Z=18054
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18055
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18056
                    }/*6*/  //Z=18057
                    pqsum = 0.0;  //Z=18058
                    oldpqsum = 0.0;  //Z=18059
                    for ( nser=0; nser<=80; nser++ )
                    {/*6*/  //Z=18060
                        binsum = 0.0;  //Z=18061
                        for ( mser=0; mser<=80; mser++ ) binsum = binsum+params.CR->carr5p[mser]*params.CR->carr11pm[nser][mser]*qyn[mser];  //Z=18062
                        pqsum = pqsum+params.CR->carr4p[nser]*qxn[nser]*binsum;  //Z=18063
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18064
                        if ( delser<0.0001 ) break; /* goto 261; */  //Z=18065
                        oldpqsum = pqsum;  //Z=18066
                    }/*6*/  //Z=18067
                    /*261:*/  //Z=18068
                    pql = pqsum;  //Z=18069
                }/*5*/  //Z=18070
                else
                {/*5*/  //Z=18071
                    a1 = 9*pow(zr+1,4)/(2.0*zr*(zr-1)*(zr-2)*(zr-3));  //Z=18072
                    pql = a1/(sqr(qy*qy*sqr(params.radius)+qx*qx*sqr(params.length))+eps9);  //Z=18073
                }/*5*/  //Z=18074
                /*formpq:=*/ return pql;  //Z=18075
            }/*4*/  /*  of orcase=1  */  //Z=18076
        }/*3*/  /*  of perfect ellipsoid  */  //Z=18077


        /*  general  */  //Z=18080
        if ( ordis==0 )
        {/*3*/  //Z=18081
            if ( params.orcase==4 )
                pql = 1.0;  //Z=18082
            else
            {/*4*/  //Z=18083
                if ( sqrt(qx*qx*sqr(params.length)+qy*qy*sqr(params.radius)+eps9)<10 )
                {/*5*/  //Z=18084
                    pqsum = 0.0;  //Z=18085
                    oldpqsum = -10.0;  //Z=18086
                    qxn[0] = 1.0;  //Z=18087
                    qyn[0] = 1.0;  //Z=18088
                    qqn[0] = 1.0;  //Z=18089

                    if ( params.orcase==1 )
                    {/*6*/  //Z=18091
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18092
                            qxn[nser] = qxn[nser-1]*qxs*qxs/(q*q);  //Z=18093
                            qyn[nser] = qyn[nser-1]*qys*qys/(q*q);  //Z=18094
                            qqn[nser] = qqn[nser-1]*q*q;  //Z=18095
                        }/*7*/  //Z=18096
                        for ( nser=0; nser<=120; nser++ )
                        {/*7*/  //Z=18097
                            binsum1 = 0.0;  //Z=18098
                            for ( lser=0; lser<=nser; lser++ )  //Z=18099
                                binsum1 = binsum1+params.CR->carr1p[lser]*qxn[lser]*qyn[nser-lser];  //Z=18100
                            binsum = 0.0;  //Z=18101
                            for ( mser=0; mser<=120; mser++ )  //Z=18102
                                binsum = binsum+params.CR->carr2p[mser]*qqn[mser]*params.CR->carr11pm[nser][mser];  //Z=18103
                            pqsum = pqsum+params.CR->carr3p[nser]*qqn[nser]*binsum*binsum1/pow(4.0,nser);  //Z=18104  TODO es stand: pow(4,n)
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18105
                            if ( delser<0.0001 ) break; /* goto 273; */  //Z=18106
                            oldpqsum = pqsum;  //Z=18107
                        }/*7*/  //Z=18108
                    }/*6*/  //Z=18109

                    if ( params.orcase==2 )
                    {/*6*/  /*  x-axis  */  //Z=18111
                        for ( nser=1; nser<=110; nser++ )
                        {/*7*/  //Z=18112
                            qxn[nser] = qxn[nser-1]*qxs*qxs/(q*q);  //Z=18113
                            qyn[nser] = qyn[nser-1]*qys*qys/(q*q);  //Z=18114
                            qqn[nser] = qqn[nser-1]*q*q;  //Z=18115
                        }/*7*/  //Z=18116
                        for ( nser=0; nser<=100; nser++ )
                        {/*7*/  //Z=18117
                            binsum1 = 0.0;  //Z=18118
                            for ( lser=0; lser<=nser; lser++ )  //Z=18119
                                /* binsum1:=binsum1+carr4p[lser]*qxn[lser]*qyn[nser-lser];  //Z=18120 */
                                binsum1 = binsum1+params.CR->carr22pm[nser][lser]*qxn[lser]*qyn[nser-lser];  //Z=18121
                            binsum = 0.0;  //Z=18122
                            for ( mser=0; mser<=100; mser++ )  //Z=18123
                                binsum = binsum+params.CR->carr5p[mser]*qqn[mser]*params.CR->carr11pm[nser][mser];  //Z=18124
                            pqsum = pqsum+params.CR->carr6p[nser]*qqn[nser]*binsum*binsum1/pow(4.0,nser);  //Z=18125
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18126
                            if ( delser<0.0001 ) break; /* goto 273; */  //Z=18127
                            oldpqsum = pqsum;  //Z=18128
                        }/*7*/  //Z=18129
                    }/*6*/  //Z=18130

                    if ( params.orcase==3 )
                    {/*6*/  /*  y-axis  */  //Z=18132
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18133
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18134
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18135
                            binsum = 0.0;  //Z=18136
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18137
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18138 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18139 */
                                binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qyn[mser]*qxn[nser-mser];  //Z=18140
                            }/*8*/  //Z=18141
                            pqsum = pqsum+params.CR->carr1p[nser]*binsum;  //Z=18142
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18143
                            if ( delser<0.0001 ) break; /* goto 273; */  //Z=18144
                            oldpqsum = pqsum;  //Z=18145
                        }/*7*/  //Z=18146
                    }/*6*/  //Z=18147
                    /*273:*/  //Z=18148
                    pql = pqsum;  //Z=18149
                }/*5*/   /*  of q<lim  */  //Z=18150
                else
                {/*5*/  //Z=18151
                    qrombdeltac(params.p1,sigmal,params.alphash1,params.polTheta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,4,16,0,0,0,params.CR->carr1p,pql);  //Z=18152
                    pql = pql/params.norm;  //Z=18153
                }/*5*/  //Z=18154
            }/*4*/  /*  of orcase=1,2,3  */  //Z=18155
            /*formpq:=*/ return pql;  //Z=18156
        }/*3*/   /*  of general  */  //Z=18157
    }/*2*/  /*  of biaxial ellipsoid  */  //Z=18158


    /* ** triaxial ellipsoid ** */  //Z=18161
    if ( params.part==6 )
    {/*2*/  //Z=18162
        /*  homogeneous isotropic triaxial ellipsoid  */  //Z=18163
        if ( ordis==7 )
        {/*3*/  //Z=18164
            if ( params.cs==0 )
            {/*4*/  //Z=18165
                if ( q<0.05*params.limq4 )
                {/*5*/  //Z=18166
                    pqsum = 1.0;  //Z=18167
                    oldpqsum = 0.0;  //Z=18168
                    qqn[0] = 1.0;  //Z=18169
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18170
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=18171
                        pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=18172
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18173
                        if ( delser<0.0001 ) break; /* goto 263; */  //Z=18174
                        oldpqsum = pqsum;  //Z=18175
                    }/*6*/  //Z=18176
                    /*263:*/  //Z=18177
                    /*formpq:=*/ return pqsum;  //Z=18178
                }/*5*/  //Z=18179
                else
                {/*5*/  //Z=18180
                    if ( q>2*params.limq4 )
                        pq = params.por/(q*q*q*q);  //Z=18181
                    else
                    {/*6*/  //Z=18182
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,7,14,7,0,0,params.CR->carr1p,pql);  //Z=18183
                        pq = pql/(M_PI/2.0);  //Z=18184
                    }/*6*/  //Z=18185
                    /*formpq:=*/ return pq;  //Z=18186
                }/*5*/  //Z=18187
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18188

            /*  core/shell isotropic ellipsoid  */  //Z=18190
            if ( params.cs==1 )
            {/*4*/  //Z=18191
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18192
            }/*4*/  //Z=18193
        }/*3*/  /*  of isotropic triaxial ellipsoid  */  //Z=18194

        /*  perfectly oriented triaxial ellipsoid  */  //Z=18196
        if ( ordis==6 )
        {/*3*/  //Z=18197
            /* if (orcase=1) then begin  //Z=18198 */
            if ( 1==1 )
            {/*4*/  //Z=18199

                if ( q<(1/params.radius) )
                {/*5*/  //Z=18201
                    pqsum = 1.0;  //Z=18202
                    oldpqsum = 0.0;  //Z=18203
                    qxn[0] = 1.0;  //Z=18204
                    qyn[0] = 1.0;  //Z=18205
                    for ( nser=1; nser<=80; nser++ )
                    {/*6*/  //Z=18206
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18207
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18208

                        binsum = 0.0;  //Z=18210
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=18211
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18212 */
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=18213 */
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=18214
                        }/*7*/  //Z=18215
                        pqsum = pqsum+binsum;  //Z=18216
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18217
                        if ( delser<0.0001 ) break; /* goto 264; */  //Z=18218
                        oldpqsum = pqsum;  //Z=18219
                    }/*6*/  //Z=18220
                    /*264:*/  //Z=18221
                    pql = pqsum;  //Z=18222
                }/*5*/  //Z=18223
                else
                {/*5*/  //Z=18224
                    argqx = qxs*params.radius/(zr+1)+eps9;  //Z=18225
                    argqy = qys*params.radius/(zr+1)+eps9;  //Z=18226
                    pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18227
                    pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18228
                    pql = pqrx*pqry;  //Z=18229
                }/*5*/  //Z=18230
                /*formpq:=*/ return pql;  //Z=18231
            }/*4*/  /*  of orcase=1  */  //Z=18232
        }/*3*/  /*  of perfect triaxial ellipsoid  */  //Z=18233
    }/*2*/  /*  of triaxial ellipsoid  */  //Z=18234


    /* ** super ellipsoid, barrel ** */  //Z=18237
    if ( params.part==7 )
    {/*2*/  //Z=18238
        /*  homogeneous isotropic super ellipsoid  */  //Z=18239
        if ( ordis==7 )
        {/*3*/  //Z=18240
            if ( params.cs==0 )
            {/*4*/  //Z=18241
                if ( q<0.9*params.limq4 )
                {/*5*/  //Z=18242
                    pqsum = 1.0;  //Z=18243
                    oldpqsum = 0.0;  //Z=18244
                    qqn[0] = 1.0;  //Z=18245
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18246
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=18247
                        pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=18248
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18249
                        if ( delser<0.0001 ) break; /* goto 270; */  //Z=18250
                        oldpqsum = pqsum;  //Z=18251
                    }/*6*/  //Z=18252
                    /*270:*/  //Z=18253
                    /*formpq:=*/ return pqsum;  //Z=18254
                }/*5*/  //Z=18255
                else
                {/*5*/  //Z=18256
                    //if ( q>=0.9*limq4 ) <<< macht keinen Sinn
                    pq = params.por/(q*q*q*q);  //Z=18257
                    /* else begin  //Z=18258 */
                    /*    qrombdeltac(length,radius,p1,sigma,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,7,14,7,0,0,carr1p,pql);  //Z=18259 */
                    /*    pq:=pql/(pi/2);  //Z=18260 */
                    /* end;  //Z=18261 */
                    /*formpq:=*/ return pq;  //Z=18262
                }/*5*/  //Z=18263
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18264

            /*  core/shell isotropic ellipsoid  */  //Z=18266
            if ( params.cs==1 )
            {/*4*/  //Z=18267
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18268
            }/*4*/  //Z=18269
        }/*3*/  /*  of isotropic triaxial ellipsoid  */  //Z=18270

        /*  perfectly oriented super ellipsoid  */  //Z=18272
        if ( ordis==6 )
        {/*3*/  //Z=18273
            /* if (orcase=1) then begin  //Z=18274 */
            if ( 1==1 )
            {/*4*/  //Z=18275

                if ( q<(1/params.radius) )
                {/*5*/  //Z=18277
                    pqsum = 1.0;  //Z=18278
                    oldpqsum = 0.0;  //Z=18279
                    qxn[0] = 1.0;  //Z=18280
                    qyn[0] = 1.0;  //Z=18281
                    for ( nser=1; nser<=80; nser++ )
                    {/*6*/  //Z=18282
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18283
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18284

                        binsum = 0.0;  //Z=18286
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=18287
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18288 */
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=18289 */
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=18290
                        }/*7*/  //Z=18291
                        pqsum = pqsum+binsum;  //Z=18292
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18293
                        if ( delser<0.0001 ) break; /* goto 271; */  //Z=18294
                        oldpqsum = pqsum;  //Z=18295
                    }/*6*/  //Z=18296
                    /*271:*/  //Z=18297
                    pql = pqsum;  //Z=18298
                }/*5*/  //Z=18299
                else
                {/*5*/  //Z=18300
                    argqx = qxs*params.radius/(zr+1)+eps9;  //Z=18301
                    argqy = qys*params.radius/(zr+1)+eps9;  //Z=18302
                    pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18303
                    pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18304
                    pql = pqrx*pqry;  //Z=18305
                }/*5*/  //Z=18306
                /*formpq:=*/ return pql;  //Z=18307
            }/*4*/  /*  of orcase=1  */  //Z=18308
        }/*3*/  /*  of perfect triaxial ellipsoid  */  //Z=18309
    }/*2*/  /*  of super ellipsoid  */  //Z=18310

    /* ** superball ** */  //Z=18312
    if ( params.part==8 )
    {/*2*/  //Z=18313
        /*  homogeneous isotropic superball  */  //Z=18314
        if ( ordis==7 )
        {/*3*/  //Z=18315
            if ( params.cs==0 )
            {/*4*/  //Z=18316
                if ( q<1.1*params.limq4 )
                {/*5*/  //Z=18317
                    pqsum = 1.0;  //Z=18318
                    oldpqsum = 0.0;  //Z=18319
                    qqn[0] = 1.0;  //Z=18320
                    for ( nser=1; nser<=40; nser++ )
                    {/*6*/  //Z=18321
                        qqn[nser] = qqn[nser-1]*q*q;  //Z=18322
                        pqsum = pqsum+params.CR->carr4p[nser]*qqn[nser];  //Z=18323
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18324
                        if ( delser<0.0001 ) break; /* goto 274; */  //Z=18325
                        oldpqsum = pqsum;  //Z=18326
                    }/*6*/  //Z=18327
                    /*274:*/  //Z=18328
                    /*formpq:=*/ return pqsum;  //Z=18329
                }/*5*/  //Z=18330
                else
                {/*5*/  //Z=18331
                    //if ( q>=1.1*limq4 ) <<< macht keinen Sinn
                    pq = params.por/(q*q*q*q);  //Z=18332
                    /* else begin  //Z=18333 */
                    /*    qrombdeltac(length,radius,p1,sigma,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,7,14,7,0,0,carr1p,pql);  //Z=18334 */
                    /*    pq:=pql/(pi/2);  //Z=18335 */
                    /* end;  //Z=18336 */
                    /*formpq:=*/ return pq;  //Z=18337
                }/*5*/  //Z=18338
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18339

            /*  core/shell isotropic ellipsoid  */  //Z=18341
            if ( params.cs==1 )
            {/*4*/  //Z=18342
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18343
            }/*4*/  //Z=18344
        }/*3*/  /*  of isotropic superball  */  //Z=18345

        /*  perfectly oriented superball  */  //Z=18347
        if ( ordis==6 )
        {/*3*/  //Z=18348
            /* if (orcase=1) then begin  //Z=18349 */
            if ( 1==1 )
            {/*4*/  //Z=18350

                if ( q<(1/params.radius) )
                {/*5*/  //Z=18352
                    pqsum = 1.0;  //Z=18353
                    oldpqsum = 0.0;  //Z=18354
                    qxn[0] = 1.0;  //Z=18355
                    qyn[0] = 1.0;  //Z=18356
                    for ( nser=1; nser<=80; nser++ )
                    {/*6*/  //Z=18357
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18358
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18359

                        binsum = 0.0;  //Z=18361
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=18362
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18363 */
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=18364 */
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=18365
                        }/*7*/  //Z=18366
                        pqsum = pqsum+binsum;  //Z=18367
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18368
                        if ( delser<0.0001 ) break; /* goto 275; */  //Z=18369
                        oldpqsum = pqsum;  //Z=18370
                    }/*6*/  //Z=18371
                    /*275:*/  //Z=18372
                    pql = pqsum;  //Z=18373
                }/*5*/  //Z=18374
                else
                {/*5*/  //Z=18375
                    argqx = qxs*params.radius/(zr+1)+eps9;  //Z=18376
                    argqy = qys*params.radius/(zr+1)+eps9;  //Z=18377
                    pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18378
                    pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18379
                    pql = pqrx*pqry;  //Z=18380
                }/*5*/  //Z=18381
                /*formpq:=*/ return pql;  //Z=18382
            }/*4*/  /*  of orcase=1  */  //Z=18383
        }/*3*/  /*  of perfect triaxial ellipsoid  */  //Z=18384
    }/*2*/  /*  of superball  */  //Z=18385

#ifndef __CUDACC__
    qDebug() << "formpq: part unknown" << params.part;
#endif
    return 0.0;
}/*1*/  //Z=18387








#ifdef __CUDACC__
__host__ __device__
#endif
//double SasCalc_GENERIC_calculation::formfq( double length, double radius, double sigmal, double sigmar, double p1, double rho,
//                                        double alfa, double theta, double phi, double limql, double limq1, double /*limq2*/,
//                                        double /*limq3*/, double limq4, double limq5, double limq6, double qx, double qy,
//                                        double qxs, double qys, double q, double norm,
//                                        int part, int cs, int ordis, int orcase,
//                                        const double */*myarray*/, // CoeffArrayType
//                                        double *carr1p, double *carr2p, double */*carr3p*/, double *carr4p, double *carr5p,
//                                        double *carr6p, double */*carr7p*/ //: CoeffArrayType;   /*Z0311=17704*/,
//                                        /*ArrayImax2D &carr11pm, ArrayImax2D &carr22pm*/ ) const //: ArrayImax2D);   /*Z0311=17705*/
double SasCalc_GENERIC_calculation::formfq( double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const
{  //Z=18395

    // Tests mit LOAD PARAMS C:\SimLab\sas-crystal\Manuskript HyperSeries (Prof.Förster)\SimulTest(Disks-gaus).ini
    // Gespräch mit Prof. Förster (05.Jun.2023): Es sollten immer die mit f hinten genutzt werden...

    int    n, /*m,*/ nser, mser, lser; //, indx;  //Z=18400
    double pqsum, oldpqsum, binsum, delser, argq, arglq, /*argp1,*/ pqr, pqr1, pqr2, pql, pq1, pq2, pq3, pq4, pq5, pq6;  //Z=18401
    double ccc1, ccc2, ccc3, vv3, zl, zr, radiusm;  //Z=18402
    double cc1, /*cc2, cc3,*/ cc4, /*cc5,*/ cc6; //, cc7, cc8, cc9, cc10;  //Z=18403
    //double ac1, ac2, ac3, ac4, ac5, ac6, ac7, ac8, ac9, ac10;  //Z=18404
    //double argbm, nenbm, argbp, nenbp, argep, nenep, argem, nenem, arggp, nengp;  //Z=18405
    double /*arggm, nengm, argim, nenim, argip, nenip,*/ F121, F122, F123;  //Z=18406
    double /*qqn[200], qqnx[200], qqny[200],*/ z12v[200], a1v[200], b1v[200], b2v[200], b1sv[200], fkv[200]; //, gam3[200];  //Z=18407
    double dim, xrad, xradp, x1z, x12z, x2z, x22z, lim, lim1, /*lim2, lim3,*/ lim4, /*lim5,*/ lim6;  //Z=18408
    double a1, b1, b2, b1s, v, c, /*d0, d1,*/ e0, e1, ee0, ee1;  //Z=18409
    double gb1s, pz2v, pz2v1, /*pz2v2,*/ gz1, preg1, preg3, preg4; //, pzvc, pzvc1, pzvc2, pzac, pzac1, pzac2;  //Z=18410
    double pzc, pzc1, pza, /*pzva, pzva1, dnv0, pvav0, pvav10, pva0,*/ sumc;  //Z=18411
    double del, delc, F12, F12sez, oldF12sez, F42, F42sez, oldF42sez, F62, F62sez, oldF62sez;  //Z=18412
    double arg11, nen11, arg12, nen12, F12as1z, F12as2z, F12asz;  //Z=18413
    double arg44, nen44, arg45, nen45, F42as10z, /*F42as1sumz, F42as1z,*/ F42as1z0, F42as40z, F42as27, F42as28, F42as4z, /*F42asz,*/ F42asz0;  //Z=18414
    double arg64, nen64, arg65, nen65, F62as10z, /*F62as1z,*/ F62as1z0, F62as40z, /*F62as1sumz,*/ F62as27, F62as28, F62as4z, /*F62asz,*/ F62asz0, FF1;  //Z=18415

    zl = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=18418
    zr = (1-sqr(params.sigma))/sqr(params.sigma);  //Z=18419
    radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=18420

    // TODO Unbekannte Variablen:
    double argpq, zz, pqr3, pqr4, qnarg, binsum1;
    double qxn[121], qyn[121];
    double qz=1.0; // TODO: in qrombdeltac Aufruf verwendet, siehe auch bei formpq()
    double qxhklt=0,qyhklt=0,qzhklt=0,qhkl=0;


    /* ************ */  //Z=18422
    /* ** sphere ** */  //Z=18423
    /* ************ */  //Z=18424
    if ( params.part==0 )
    {/*2*/  //Z=18425
        /* ** homogeneous sphere ** */  //Z=18426
        if ( params.cs==0 )
        {/*3*/  //Z=18427
            if ( q<(0.4*params.limq4f) )
            {/*4*/  //Z=18428
                pqsum = 1.0;  //Z=18429
                oldpqsum = 0.0;  //Z=18430
                double qqnn = 1.0;  //Z=18431
                for ( nser=1; nser<=100; nser++ )
                {/*5*/  //Z=18432
                    qqnn = qqnn*q*q;  //Z=18433
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=18434
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18435
                    if ( delser<0.0001 ) break; /* goto 50; */  //Z=18436
                    oldpqsum = pqsum;  //Z=18437
                }/*5*/  //Z=18438
                /*50:*/  //Z=18439
                /*formfq:=*/ return pqsum;  //Z=18440
            }/*4*/  //Z=18441
            else
            {/*4*/  //Z=18442
                argq = q*params.radius/(zr+1);  //Z=18443
                pqr = (1/(zr*(zr-1))*pow(argq,-2));  //Z=18444
                pq1 = pqr*cos((zr-1)*atan(argq))/pow(1.0+argq*argq,(zr-1)/2.0);  //Z=18445
                pq2 = (pqr/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1.0+argq*argq,(zr-2)/2.0);  //Z=18446
                pq3 = 3*(pq2-pq1);  //Z=18447
                /*formfq:=*/ return pq3*pq3;  //Z=18448
            }/*4*/  //Z=18449
        }/*3*/ /*  of homogeneous sphere */  //Z=18450

        /* ** core/shell sphere ** */  //Z=18452
        if ( params.cs==1 )
        {/*3*/  //Z=18453
            ccc1 = sqr(1-params.rho)*pow(params.p1,6);  //Z=18454
            ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,3);  //Z=18455
            ccc3 = sqr(params.rho);  //Z=18456
            vv3 = sqr((1-params.rho)*pow(params.p1,3)+params.rho);  //Z=18457

            argq = q*radiusm/(zr+1);  //Z=18459
            argpq = q*params.radius/(zr+1);  //Z=18460

            /*  F121 sphere  */  //Z=18462
            if ( q<(params.limq4f/2.0) )
            {/*4*/  //Z=18463
                double qqnn = 1.0;  //Z=18464
                pqsum = 1.0;  //Z=18465
                oldpqsum = 0.0;  //Z=18466
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18467
                    qqnn = qqnn*q*q;  //Z=18468
                    pqsum = pqsum+qqnn*params.CR->carr4f[nser];  //Z=18469
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18470
                    if ( delser<0.0001 ) break; /* goto 51; */  //Z=18471
                    oldpqsum = pqsum;  //Z=18472
                }/*5*/  //Z=18473
                /*51:*/  //Z=18474
                F121 = ccc1*pqsum/vv3;  //Z=18475
            }/*4*/  //Z=18476
            else
            {/*4*/  //Z=18477
                pqr = (1/(zr*(zr-1))*pow(argpq,-2));  //Z=18478
                pq1 = pqr*cos((zr-1)*atan(argpq))/pow(1.0+argpq*argpq,(zr-1)/2.0);  //Z=18479
                pq2 = (pqr/((zr-2)*argpq))*sin((zr-2)*atan(argpq))/pow(1.0+argpq*argpq,(zr-2)/2.0);  //Z=18480
                pq3 = 3*(pq2-pq1);  //Z=18481
                F121 = ccc1*pq3*pq3/vv3;  //Z=18482
            }/*4*/  //Z=18483

            /*  F122 sphere  */  //Z=18485
            if ( q<(0.3*params.limq5f) )
            {/*4*/  //Z=18486
                double qqnn = 1.0;  //Z=18487
                pqsum = 1.0;  //Z=18488
                oldpqsum = 0.0;  //Z=18489
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18490
                    qqnn = qqnn*q*q;  //Z=18491
                    pqsum = pqsum+qqnn*params.CR->carr5f[nser];  //Z=18492
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18493
                    if ( delser<0.0001 ) break; /* goto 52; */  //Z=18494
                    oldpqsum = pqsum;  //Z=18495
                }/*5*/  //Z=18496
                /*52:*/  //Z=18497
                F122 = ccc2*pqsum/vv3;  //Z=18498
            }/*4*/  //Z=18499
            else
            {/*4*/  //Z=18500
                pqr1 = (1/(zr*(zr-1))*pow(argpq,-2));  //Z=18501
                pq1 = pqr1*cos((zr-1)*atan(argpq))/pow(1.0+argpq*argpq,(zr-1)/2.0);  //Z=18502
                pq2 = (pqr1/((zr-2)*argpq))*sin((zr-2)*atan(argpq))/pow(1.0+argpq*argpq,(zr-2)/2.0);  //Z=18503
                pq3 = 3*(pq2-pq1);  //Z=18504

                pqr2 = (1/(zr*(zr-1))*pow(argq,-2));  //Z=18506
                pq4 = pqr2*cos((zr-1)*atan(argq))/pow(1.0+argq*argq,(zr-1)/2.0);  //Z=18507
                pq5 = (pqr2/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1.0+argq*argq,(zr-2)/2.0);  //Z=18508
                pq6 = 3*(pq5-pq4);  //Z=18509

                F122 = ccc2*pq3*pq6/vv3;  //Z=18511
            }/*4*/  //Z=18512

            /*  F123 sphere  */  //Z=18514
            if ( q<(params.limq6f/2.0) )
            {/*4*/  //Z=18515
                double qqnn = 1.0;  //Z=18516
                pqsum = 1.0;  //Z=18517
                oldpqsum = 0.0;  //Z=18518
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18519
                    qqnn = qqnn*q*q;  //Z=18520
                    pqsum = pqsum+qqnn*params.CR->carr6f[nser];  //Z=18521
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18522
                    if ( delser<0.0001 ) break; /* goto 53; */  //Z=18523
                    oldpqsum = pqsum;  //Z=18524
                }/*5*/  //Z=18525
                /*53:*/  //Z=18526
                F123 = ccc3*pqsum/vv3;  //Z=18527
            }/*4*/  //Z=18528
            else
            {/*4*/  //Z=18529
                pqr = (1/(zr*(zr-1))*pow(argq,-2));  //Z=18530
                pq1 = pqr*cos((zr-1)*atan(argq))/pow(1.0+argq*argq,(zr-1)/2.0);  //Z=18531
                pq2 = (pqr/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1.0+argq*argq,(zr-2)/2.0);  //Z=18532
                pq3 = 3*(pq2-pq1);  //Z=18533
                F123 = ccc1*pq3*pq3/vv3;  //Z=18534
            }/*4*/  //Z=18535
            /* formfq:=F121+F122+F123;  //Z=18536 */
            /*formfq:=*/ return F123;  //Z=18537
        }/*3*/  /*  of core/shell sphere  */  //Z=18538

        /* ** inhomogeneous core/shell sphere ** */  //Z=18540
        if ( params.cs==2 )
        {/*3*/  //Z=18541

            dim = 3;  //Z=18543
            delc = 0.0001;  //Z=18544
            zz = zr;  //Z=18545
            xrad = q*radiusm;  //Z=18546
            xradp = q*params.radius;  //Z=18547
            x1z = q*params.radius/(2.0*(zz+1));  //Z=18548
            x12z = x1z*x1z;  //Z=18549
            x2z = q*radiusm/(2.0*(zz+1));  //Z=18550
            x22z = x2z*x2z;  //Z=18551

            lim = 18*exp(-5*params.sigma);  //Z=18553
            lim1 = lim;  //Z=18554
            //lim2 = lim*0.7;  //Z=18555
            //lim3 = lim;  //Z=18556
            lim4 = lim;  //Z=18557
            //lim5 = lim*0.7;  //Z=18558
            lim6 = lim*1.2;  //Z=18559

            a1 = (dim-params.alphash1)/2.0;  //Z=18561
            b1 = dim/2.0;  //Z=18562
            b2 = (dim+2-params.alphash1)/2.0;  //Z=18563
            b1s = (dim+2)/2.0;  //Z=18564
            v = -b1s+1/2.0;  //Z=18565
            c = a1-b1-b2+1/2.0;  //Z=18566
            //d0 = 1;  //Z=18567
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=18568
            e0 = 1.0;  //Z=18569
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=18570
            ee0 = 1.0;  //Z=18571
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=18572

            gb1s = 3*sqrt(M_PI)/4.0;  //Z=18574
            pz2v = 1/(zr*(zr-1));  //Z=18575
            pz2v1 = pz2v/(zr-2);  //Z=18576
            //pz2v2 = pz2v1/(zr-3);  //Z=18577

            gz1 = gamma(zr+1);  //Z=18579
            preg1 = gb1s/sqrt(M_PI);  //Z=18580
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=18581
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=18582
            //pzvc = gamma(zr+1+v+c)/gz1;  //Z=18583
            //pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=18584
            //pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=18585
            //pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=18586
            //pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=18587
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=18588
            pzc = gamma(zr+1+c)/gz1;  //Z=18589
            pzc1 = gamma(zr+1+c-1)/gz1;  //Z=18590
            pza = gamma(zr+1-2*a1)/gz1;  //Z=18591
            //pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=18592
            //pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=18593
            //dnv0 = 1;  //Z=18594
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=18595
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=18596
            //pva0 = gamma(zr+1-2*a1)/gz1;  //Z=18597

            cc1 = 1/dim;  //Z=18599
            cc4 = params.rho/((dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=18600
            cc6 = -params.rho/(dim-params.alphash1);  //Z=18601
            sumc = cc1+cc4+cc6;  //Z=18602

            /*  term #1 series  */  //Z=18604
            if ( (xradp)<lim1 )
            {/*4*/  //Z=18605
                z12v[0] = 1;  //Z=18606
                b1sv[0] = 1;  //Z=18607
                fkv[0] = 1;  //Z=18608
                double qqnn = 1.0;  //Z=18609
                F12sez = 1.0;  //Z=18610
                oldF12sez = 1.0;  //Z=18611
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=18612
                    qqnn = qqnn*q*q;  //Z=18613 war: qnn[..], ist zwar definiert, aber oben wurde qqn[0]=1 gesetzt... (TODO)
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=18614
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=18615
                    fkv[n] = fkv[n-1]*n;  //Z=18616
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]);  //Z=18617 */

                    F12sez = F12sez+params.CR->carr4f[n]*qqnn;  //Z=18619

                    del = fabs((F12sez-oldF12sez)/F12sez);  //Z=18621
                    if ( del<delc ) break; /* goto 101; */  //Z=18622
                    oldF12sez = F12sez;  //Z=18623
                }/*5*/  //Z=18624
                /*101:*/  //Z=18625
                F12 = F12sez;  //Z=18626
            }/*4*/  //Z=18627

            /*  term #4 series  */  //Z=18629
            if ( xradp<lim4 )
            {/*4*/  //Z=18630
                z12v[0] = 1;  //Z=18631
                a1v[0] = 1;  //Z=18632
                b1v[0] = 1;  //Z=18633
                b2v[0] = 1;  //Z=18634
                b1sv[0] = 1;  //Z=18635
                fkv[0] = 1;  //Z=18636
                double qqnn = 1.0;  //Z=18637
                F42sez = 1.0;  //Z=18638
                oldF42sez = 1.0;  //Z=18639
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=18640
                    qqnn = qqnn*q*q;  //Z=18641
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=18642
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=18643
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=18644
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=18645
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=18646
                    fkv[n] = fkv[n-1]*n;  //Z=18647
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=18648 */

                    F42sez = F42sez+params.CR->carr5f[n]*qqnn;  //Z=18650

                    del = fabs((F42sez-oldF42sez)/F42sez);  //Z=18652
                    if ( del<delc ) break; /* goto 104; */  //Z=18653
                    oldF42sez = F42sez;  //Z=18654
                }/*5*/  //Z=18655
                /*104:*/  //Z=18656
                F42 = F42sez;  //Z=18657
            }/*4*/  //Z=18658

            /*  term #6 series  */  //Z=18660
            if ( xradp<lim6 )
            {/*4*/  //Z=18661
                z12v[0] = 1;  //Z=18662
                a1v[0] = 1;  //Z=18663
                b1v[0] = 1;  //Z=18664
                b2v[0] = 1;  //Z=18665
                b1sv[0] = 1;  //Z=18666
                fkv[0] = 1;  //Z=18667
                double qqnn = 1.0;  //Z=18668
                F62sez = 1.0;  //Z=18669
                oldF62sez = 1.0;  //Z=18670
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=18671
                    qqnn = qqnn*q*q;  //Z=18672
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=18673
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=18674
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=18675
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=18676
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=18677
                    fkv[n] = fkv[n-1]*n;  //Z=18678
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=18679 */

                    F62sez = F62sez+params.CR->carr6f[n]*qqnn;  //Z=18681

                    del = fabs((F62sez-oldF62sez)/F62sez);  //Z=18683
                    if ( del<delc ) break; /* goto 106; */  //Z=18684
                    oldF62sez = F62sez;  //Z=18685
                }/*5*/  //Z=18686
                /*106:*/  //Z=18687
                F62 = F62sez;  //Z=18688
            }/*4*/  //Z=18689

            /* ** term #1 asymptote ** */  //Z=18691
            if ( xradp>=lim1 )
            {/*4*/  //Z=18692
                arg11 = (zr+v+1)*atan(2.0*x1z);  //Z=18693
                nen11 = pow(1.0+4*x1z*x1z,(zr+v+1)/2.0);  //Z=18694
                arg12 = (zr+v)*atan(2.0*x1z);  //Z=18695
                nen12 = pow(1.0+4*x1z*x1z,(zr+v)/2.0);  //Z=18696

                F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);  //Z=18698
                F12as2z = ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);  //Z=18699
                F12asz = preg1*pow(x1z,v)*(F12as1z+F12as2z);  //Z=18700
                F12 = F12asz;  //Z=18701
            }/*4*/  //Z=18702

            /* ** term #4 asymptote ** */  //Z=18704
            if ( xrad>=lim4 )
            {/*4*/  //Z=18705
                F42as10z = preg4*pow(x22z,-a1);  //Z=18706
                //F42as1sumz = pva0;  //Z=18707
                //F42as1z = F42as10z*F42as1sumz;  //Z=18708
                F42as1z0 = F42as10z*pza;   /* * */  //Z=18709

                F42as40z = preg3*pow(x2z,c);  //Z=18711
                arg44 = (zr+c+1)*atan(2.0*x2z);  //Z=18712
                nen44 = pow(1.0+4*x2z*x2z,(zr+c+1)/2.0);  //Z=18713
                arg45 = (zr+c)*atan(2.0*x2z);  //Z=18714
                nen45 = pow(1.0+4*x2z*x2z,(zr+c)/2.0);  //Z=18715
                F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);  //Z=18716
                F42as28 = e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);  //Z=18717
                F42as4z = F42as40z*(F42as27+F42as28);  //Z=18718
                //F42asz = F42as1z+F42as4z;  //Z=18719
                F42asz0 = F42as1z0+F42as4z;  //Z=18720
                F42 = F42asz0;  //Z=18721
            }/*4*/  //Z=18722

            /* ** term #6 asymptote ** */  //Z=18724
            if ( xradp>=lim6 )
            {/*4*/  //Z=18725
                F62as10z = preg4*pow(x12z,-a1);  //Z=18726
                //F62as1sumz = pva0;  //Z=18727
                //F62as1z = F62as10z*F62as1sumz;  //Z=18728
                F62as1z0 = F62as10z*pza;     /* * */  //Z=18729

                F62as40z = preg3*pow(x1z,c);  //Z=18731
                arg64 = (zr+c+1)*atan(2.0*x1z);  //Z=18732
                nen64 = pow(1.0+4*x1z*x1z,(zr+c+1)/2.0);  //Z=18733
                arg65 = (zr+c)*atan(2.0*x1z);  //Z=18734
                nen65 = pow(1.0+4*x1z*x1z,(zr+c)/2.0);  //Z=18735
                F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);  //Z=18736
                F62as28 = e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);  //Z=18737
                F62as4z = F62as40z*(F62as27+F62as28);  //Z=18738
                //F62asz = F62as1z+F62as4z;  //Z=18739
                F62asz0 = F62as1z0+F62as4z;  //Z=18740
                F62 = F62asz0;  //Z=18741
            }/*4*/  //Z=18742

            /* FF1:=(cc1*F12+cc4*F42+cc6*F62)/sumc;  //Z=18744 */
            FF1 = (cc1*F12)/sumc;  //Z=18745

            /*formfq:=*/ return FF1*FF1;  //Z=18747


            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,3,sigmar,q);  //Z=18750 */


        }/*3*/  /*  of inhomogeneous core/shell sphere  */  //Z=18753

    }/*2*/ /*  of sphere  */  //Z=18755


    /* ********** */  //Z=18758
    /*  cylinder  */  //Z=18759
    /* ********** */  //Z=18760
    if ( params.part==1 )
    {/*2*/  //Z=18761

        /* ** longitudinal part ** */  //Z=18763
        /* ** isotropic ** */  //Z=18764
        if ( ordis==7 )
        {/*3*/  //Z=18765
            if ( q<(0.6*params.limq1f) )
            {/*4*/  //Z=18766
                pqsum = 1.0;  //Z=18767
                oldpqsum = 0.0;  //Z=18768
                double qqnn = 1.0;  //Z=18769
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18770
                    qqnn = qqnn*q*q;  //Z=18771
                    pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=18772
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18773
                    if ( delser<0.0001 ) break; /* goto 60; */  //Z=18774
                    oldpqsum = pqsum;  //Z=18775
                }/*5*/  //Z=18776
                /*60:*/  //Z=18777
                pql = pqsum;  //Z=18778
            }/*4*/  //Z=18779
            else
            {/*4*/   /*  = P(q)  */  //Z=18780
                arglq = q*params.length/(zl+1);  //Z=18781
                /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2));  //Z=18782 */
                pql = (M_PI/(2.0*zl))*(1/arglq);  //Z=18783
                pql = pql-(1/(2.0*zl*(zl-1)*arglq*arglq))*cos((zl-1)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-1)/2.0);  //Z=18784
            }/*4*/  //Z=18785
        }/*3*/   /*  of isotropic  */  //Z=18786

        /*  perfect  */  //Z=18788
        if ( ordis==6 )
        {/*3*/  //Z=18789
            if ( params.orcase==4 )
                pql = 1.0;  //Z=18790
            else
            {/*4*/  //Z=18791
                if ( limql<(4*params.limq1f) )
                {/*5*/  //Z=18792
                    pqsum = 1.0;  //Z=18793
                    oldpqsum = 0.0;  //Z=18794
                    double qqnn = 1.0;  //Z=18795
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18796
                        qqnn = qqnn*(qxs+qys)*(qxs+qys);  //Z=18797
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=18798
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18799
                        if ( delser<0.0001 ) break; /* goto 65; */  //Z=18800
                        oldpqsum = pqsum;  //Z=18801
                    }/*6*/  //Z=18802
                    /*65:*/  //Z=18803
                    pql = pqsum;  //Z=18804
                }/*5*/  //Z=18805
                else
                {/*5*/  //Z=18806
                    arglq = (qxs+qys+eps9)*params.length/(zl+1);  //Z=18807
                    /*  F(q)  */  //Z=18808
                    /* pql:=(1/zl)*(1/arglq)*sin(zl*arctan(arglq))/power(1+arglq*arglq,zl/2);  //Z=18809 */
                    /* pql:=pql*pql;  //Z=18810 */
                    /*  P(q)  */  //Z=18811
                    pql = (1/(2.0*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-1)/2.0));  //Z=18812
                }/*5*/  //Z=18813
                }/*4*/  //Z=18814
        }/*3*/   /*  of perfect  */  //Z=18815

        /*  general  */  //Z=18817
        if ( ordis==0 )
        {/*3*/  //Z=18818
            if ( params.orcase==4 )
                pql = 1.0;  //Z=18819
            else
            {/*4*/  //Z=18820
                if ( limql<(0.7*params.limq1f) )
                {/*5*/  //Z=18821
                    pqsum = 1.0;  //Z=18822
                    oldpqsum = 0.0;  //Z=18823
                    qxn[0] = 1.0;  //Z=18824
                    qyn[0] = 1.0;  //Z=18825
                    if ( params.orcase==1 )
                    {/*6*/  //Z=18826
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18827
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18828
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18829
                            binsum = 0.0;  //Z=18830
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18831
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18832 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18833 */
                                binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qyn[mser]*qxn[nser-mser];  //Z=18834
                            }/*8*/  //Z=18835
                            pqsum = pqsum+params.CR->carr1f[nser]*binsum;  //Z=18836
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18837
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=18838
                            oldpqsum = pqsum;  //Z=18839
                        }/*7*/  //Z=18840
                    }/*6*/  //Z=18841

                    if ( params.orcase==2 )
                    {/*6*/  /*  x-axis  */  //Z=18843
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18844
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18845
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18846
                            binsum = 0.0;  //Z=18847
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18848
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18849 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18850 */
                                binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];  //Z=18851
                            }/*8*/  //Z=18852
                            pqsum = pqsum+params.CR->carr1f[nser]*binsum;  //Z=18853
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18854
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=18855
                            oldpqsum = pqsum;  //Z=18856
                        }/*7*/  //Z=18857
                    }/*6*/  //Z=18858
                    pql = pqsum;  //Z=18859

                    if ( params.orcase==3 )
                    {/*6*/  /*  y-axis  */  //Z=18861
                        for ( nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18862
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18863
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18864
                            binsum = 0.0;  //Z=18865
                            for ( mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18866
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18867 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18868 */
                                binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[mser]*qyn[nser-mser];  //Z=18869
                            }/*8*/  //Z=18870
                            pqsum = pqsum+params.CR->carr1f[nser]*binsum;  //Z=18871
                            delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18872
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=18873
                            oldpqsum = pqsum;  //Z=18874
                        }/*7*/  //Z=18875
                    }/*6*/  //Z=18876
                    /*66:*/  //Z=18877
                    pql = pqsum;  //Z=18878
                }/*5*/  //Z=18879
                else
                {/*5*/  //Z=18880
                    /*  F(q)  */  //Z=18881
                    /* qrombdeltac(length,sigmal,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,orcase,0,1,0,carr1p,pql);  //Z=18882 */
                    /*  P(q)  */  //Z=18883
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,0,qxs,qys,qz,
                                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,1,4,params.orcase,0,0,0,params.CR->carr1f,pql);  //Z=18884
                    pql = pql/params.norm;  //Z=18885
                }/*5*/  //Z=18886
            }/*4*/  //Z=18887
        }/*3*/   /*  of general  */  //Z=18888

        /*  transverse part  */  //Z=18890
        /*  homogeneous cylinder  */  //Z=18891
        if ( params.cs==0 )
        {/*3*/  //Z=18892
            if ( q<(1.5*params.limq4f) )
            {/*4*/  //Z=18893
                pqsum = 1.0;  //Z=18894
                oldpqsum = 0.0;  //Z=18895
                double qqnn = 1.0;  //Z=18896
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18897
                    qqnn = qqnn*q*q;  //Z=18898
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=18899
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18900
                    if ( delser<0.0001 ) break; /* goto 61; */  //Z=18901
                    oldpqsum = pqsum;  //Z=18902
                }/*5*/  //Z=18903
                /*61:*/  //Z=18904
                pqr = pqsum;  //Z=18905
            }/*4*/  //Z=18906
            else
            {/*4*/  //Z=18907
                argpq = q*params.radius/(zr+1);  //Z=18908
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-1/2.0)/2.0);  //Z=18909
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-3/2.0)/2.0);  //Z=18910
                pqr3 = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);  //Z=18911
                pqr = pqr3*pqr3;  //Z=18912
            }/*4*/  //Z=18913
            /*formfq:=*/ return pql*pqr;  //Z=18914
            /* formfq:=pql;  //Z=18915 */
        }/*3*/ /*  of homogeneous cylinder  */  //Z=18916

        /*  homogeneous core/shell cylinder  */  //Z=18918
        if ( params.cs==1 )
        {/*3*/  //Z=18919
            ccc1 = sqr(1-params.rho)*pow(params.p1,4);  //Z=18920
            ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,2);  //Z=18921
            ccc3 = sqr(params.rho);  //Z=18922
            vv3 = sqr((1-params.rho)*pow(params.p1,2)+params.rho);  //Z=18923

            zz = zr;  // TODO: zz war in diesem Zweig nicht gesetzt
            argq = q*radiusm/(zz+1);  //Z=18925
            argpq = q*params.radius/(zz+1);  //Z=18926

            /*  F121 cylinder  */  //Z=18928
            if ( q<(0.7*params.limq4f) )
            {/*4*/  //Z=18929
                /* ** series expansion ** */  //Z=18930
                pqsum = 1.0;  //Z=18931
                oldpqsum = 0.0;  //Z=18932
                double qqnn = 1.0;  //Z=18933
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18934
                    qqnn = qqnn*q*q;  //Z=18935
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=18936
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18937
                    if ( delser<0.0001 ) break; /* goto 62; */  //Z=18938
                    oldpqsum = pqsum;  //Z=18939
                }/*5*/  //Z=18940
                /*62:*/  //Z=18941
                F121 = ccc1*pqsum/vv3;  //Z=18942
            }/*4*/  //Z=18943
            else
            {/*4*/  //Z=18944
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-1/2.0)/2.0);  //Z=18945
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-3/2.0)/2.0);  //Z=18946
                pqr3 = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);  //Z=18947
                pqr = pqr3;  //Z=18948
                F121 = ccc1*pqr*pqr/vv3;  //Z=18949
            }/*4*/  //Z=18950

            /*  F122 cylinder  */  //Z=18952
            if ( q<(0.3*params.limq5f) )
            {/*4*/  //Z=18953
                /* ** series expansion ** */  //Z=18954
                pqsum = 1.0;  //Z=18955
                oldpqsum = 0.0;  //Z=18956
                double qqnn = 1.0;  //Z=18957
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18958
                    qqnn = qqnn*q*q;  //Z=18959
                    pqsum = pqsum+params.CR->carr5f[nser]*qqnn;  //Z=18960
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18961
                    if ( delser<0.0001 ) break; /* goto 63; */  //Z=18962
                    oldpqsum = pqsum;  //Z=18963
                }/*5*/  //Z=18964
                /*63:*/  //Z=18965
                F122 = ccc2*pqsum/vv3;  //Z=18966
            }/*4*/  //Z=18967
            else
            {/*4*/  //Z=18968
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1.0+argpq*argq,(zr-1/2.0)/2.0);  //Z=18969
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1.0+argpq*argq,(zr-3/2.0)/2.0);  //Z=18970
                pqr3 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argq,-3/2.0)*(sin((zr-1/2.0)*atan(argq))-cos((zr-1/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-1/2.0)/2.0);  //Z=18971
                pqr4 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argq,-5/2.0)*(sin((zr-3/2.0)*atan(argq))+cos((zr-3/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-3/2.0)/2.0);  //Z=18972
                pqr = (4/M_PI)*(pqr1+(9/16.0)*pqr2)*(pqr3+(9/16.0)*pqr4);  //Z=18973
                F122 = ccc2*pqr/vv3;  //Z=18974
            }/*4*/  //Z=18975

            /*  F123 cylinder  */  //Z=18977
            if ( q<(0.6*params.limq6f) )
            {/*4*/  //Z=18978
                /* ** series expansion ** */  //Z=18979
                pqsum = 1.0;  //Z=18980
                oldpqsum = 0.0;  //Z=18981
                double qqnn = 1.0;  //Z=18982
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18983
                    qqnn = qqnn*q*q;  //Z=18984
                    pqsum = pqsum+params.CR->carr6f[nser]*qqnn;  //Z=18985
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18986
                    if ( delser<0.0001 ) break; /* goto 64; */  //Z=18987
                    oldpqsum = pqsum;  //Z=18988
                }/*5*/  //Z=18989
                /*64:*/  //Z=18990
                F123 = ccc3*pqsum/vv3;  //Z=18991
            }/*4*/  //Z=18992
            else
            {/*4*/  //Z=18993
                pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argq,-3/2.0)*(sin((zr-1/2.0)*atan(argq))-cos((zr-1/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-1/2.0)/2.0);  //Z=18994
                pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argq,-5/2.0)*(sin((zr-3/2.0)*atan(argq))+cos((zr-3/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-3/2.0)/2.0);  //Z=18995
                pqr3 = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);  //Z=18996
                pqr = pqr3;  //Z=18997
                F123 = ccc3*pqr*pqr/vv3;  //Z=18998
            }/*4*/  //Z=18999
            /*formfq:=*/ return pql*(F121+F122+F123);  //Z=19000
            /* formfq:=(F121+F122+F123);  //Z=19001 */
        }/*3*/ /*  of homogeneous core/shell cylinder  */  //Z=19002

        /* ** inhomogeneous core/shell cylinder ** */  //Z=19004
        if ( params.cs==2 )
        {/*3*/  //Z=19005

            dim = 2;  //Z=19007
            delc = 0.0001;  //Z=19008
            zz = zr;  //Z=19009
            xrad = q*radiusm;  //Z=19010
            xradp = q*params.radius;  //Z=19011
            x1z = q*params.radius/(2.0*(zz+1));  //Z=19012
            x12z = x1z*x1z;  //Z=19013
            x2z = q*radiusm/(2.0*(zz+1));  //Z=19014
            x22z = x2z*x2z;  //Z=19015

            lim = 18*exp(-5*params.sigma);  //Z=19017
            lim1 = lim;  //Z=19018
            //lim2 = lim*0.7;  //Z=19019
            //lim3 = lim;  //Z=19020
            lim4 = lim;  //Z=19021
            //lim5 = lim*0.7;  //Z=19022
            lim6 = lim*1.2;  //Z=19023

            a1 = (dim-params.alphash1)/2.0;  //Z=19025
            b1 = dim/2.0;  //Z=19026
            b2 = (dim+2-params.alphash1)/2.0;  //Z=19027
            b1s = (dim+2)/2.0;  //Z=19028
            v = -b1s+1/2.0;  //Z=19029
            c = a1-b1-b2+1/2.0;  //Z=19030
            //d0 = 1;  //Z=19031
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=19032
            e0 = 1.0;  //Z=19033
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=19034
            ee0 = 1.0;  //Z=19035
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=19036

            gb1s = 3*sqrt(M_PI)/4.0;  //Z=19038
            pz2v = 1/(zr*(zr-1));  //Z=19039
            pz2v1 = pz2v/(zr-2);  //Z=19040
            //pz2v2 = pz2v1/(zr-3);  //Z=19041

            gz1 = gamma(zr+1);  //Z=19043
            preg1 = gb1s/sqrt(M_PI);  //Z=19044
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=19045
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=19046
            //pzvc = gamma(zr+1+v+c)/gz1;  //Z=19047
            //pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=19048
            //pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=19049
            //pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=19050
            //pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=19051
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=19052
            pzc = gamma(zr+1+c)/gz1;  //Z=19053
            pzc1 = gamma(zr+1+c-1)/gz1;  //Z=19054
            pza = gamma(zr+1-2*a1)/gz1;  //Z=19055
            //pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=19056
            //pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19057
            //dnv0 = 1;  //Z=19058
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=19059
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19060
            //pva0 = gamma(zr+1-2*a1)/gz1;  //Z=19061

            cc1 = 1/dim;  //Z=19063
            cc4 = params.rho/((dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=19064
            cc6 = -params.rho/(dim-params.alphash1);  //Z=19065
            sumc = cc1+cc4+cc6;  //Z=19066

            /*  term #1 series  */  //Z=19068
            if ( (xradp)<lim1 )
            {/*4*/  //Z=19069
                z12v[0] = 1;  //Z=19070
                b1sv[0] = 1;  //Z=19071
                fkv[0] = 1;  //Z=19072
                double qqnn = 1.0;  //Z=19073
                F12sez = 1.0;  //Z=19074
                oldF12sez = 1.0;  //Z=19075
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=19076
                    qqnn = qqnn*q*q;  //Z=19077
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19078
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19079
                    fkv[n] = fkv[n-1]*n;  //Z=19080
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]);  //Z=19081 */

                    F12sez = F12sez+params.CR->carr4f[n]*qqnn;  //Z=19083

                    del = fabs((F12sez-oldF12sez)/F12sez);  //Z=19085
                    if ( del<delc ) break; /* goto 111; */  //Z=19086
                    oldF12sez = F12sez;  //Z=19087
                }/*5*/  //Z=19088
                /*111:*/  //Z=19089
                F12 = F12sez;  //Z=19090
            }/*4*/  //Z=19091

            /*  term #4 series  */  //Z=19093
            if ( xradp<lim4 )
            {/*4*/  //Z=19094
                z12v[0] = 1;  //Z=19095
                a1v[0] = 1;  //Z=19096
                b1v[0] = 1;  //Z=19097
                b2v[0] = 1;  //Z=19098
                b1sv[0] = 1;  //Z=19099
                fkv[0] = 1;  //Z=19100
                double qqnn = 1.0;  //Z=19101
                F42sez = 1.0;  //Z=19102
                oldF42sez = 1.0;  //Z=19103
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=19104
                    qqnn = qqnn*q*q;  //Z=19105
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19106
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19107
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19108
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19109
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19110
                    fkv[n] = fkv[n-1]*n;  //Z=19111
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19112 */

                    F42sez = F42sez+params.CR->carr5f[n]*qqnn;  //Z=19114

                    del = fabs((F42sez-oldF42sez)/F42sez);  //Z=19116
                    if ( del<delc ) break; /* goto 114; */  //Z=19117
                    oldF42sez = F42sez;  //Z=19118
                }/*5*/  //Z=19119
                /*114:*/  //Z=19120
                F42 = F42sez;  //Z=19121
            }/*4*/  //Z=19122

            /*  term #6 series  */  //Z=19124
            if ( xradp<lim6 )
            {/*4*/  //Z=19125
                z12v[0] = 1;  //Z=19126
                a1v[0] = 1;  //Z=19127
                b1v[0] = 1;  //Z=19128
                b2v[0] = 1;  //Z=19129
                b1sv[0] = 1;  //Z=19130
                fkv[0] = 1;  //Z=19131
                double qqnn = 1.0;  //Z=19132
                F62sez = 1.0;  //Z=19133
                oldF62sez = 1.0;  //Z=19134
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=19135
                    qqnn = qqnn*q*q;  //Z=19136
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19137
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19138
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19139
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19140
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19141
                    fkv[n] = fkv[n-1]*n;  //Z=19142
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19143 */

                    F62sez = F62sez+params.CR->carr6f[n]*qqnn;  //Z=19145

                    del = fabs((F62sez-oldF62sez)/F62sez);  //Z=19147
                    if ( del<delc ) break; /* goto 116; */  //Z=19148
                    oldF62sez = F62sez;  //Z=19149
                }/*5*/  //Z=19150
                /*116:*/  //Z=19151
                F62 = F62sez;  //Z=19152
            }/*4*/  //Z=19153

            /* ** term #1 asymptote ** */  //Z=19155
            if ( xradp>=lim1 )
            {/*4*/  //Z=19156
                arg11 = (zr+v+1)*atan(2.0*x1z);  //Z=19157
                nen11 = pow(1.0+4*x1z*x1z,(zr+v+1)/2.0);  //Z=19158
                arg12 = (zr+v)*atan(2.0*x1z);  //Z=19159
                nen12 = pow(1.0+4*x1z*x1z,(zr+v)/2.0);  //Z=19160

                F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);  //Z=19162
                F12as2z = ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);  //Z=19163
                F12asz = preg1*pow(x1z,v)*(F12as1z+F12as2z);  //Z=19164
                F12 = F12asz;  //Z=19165
            }/*4*/  //Z=19166

            /* ** term #4 asymptote ** */  //Z=19168
            if ( xrad>=lim4 )
            {/*4*/  //Z=19169
                F42as10z = preg4*pow(x22z,-a1);  //Z=19170
                //F42as1sumz = pva0;  //Z=19171
                //F42as1z = F42as10z*F42as1sumz;  //Z=19172
                F42as1z0 = F42as10z*pza;   /* * */  //Z=19173

                F42as40z = preg3*pow(x2z,c);  //Z=19175
                arg44 = (zr+c+1)*atan(2.0*x2z);  //Z=19176
                nen44 = pow(1.0+4*x2z*x2z,(zr+c+1)/2.0);  //Z=19177
                arg45 = (zr+c)*atan(2.0*x2z);  //Z=19178
                nen45 = pow(1.0+4*x2z*x2z,(zr+c)/2.0);  //Z=19179
                F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);  //Z=19180
                F42as28 = e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);  //Z=19181
                F42as4z = F42as40z*(F42as27+F42as28);  //Z=19182
                //F42asz = F42as1z+F42as4z;  //Z=19183
                F42asz0 = F42as1z0+F42as4z;  //Z=19184
                F42 = F42asz0;  //Z=19185
            }/*4*/  //Z=19186

            /* ** term #6 asymptote ** */  //Z=19188
            if ( xradp>=lim6 )
            {/*4*/  //Z=19189
                F62as10z = preg4*pow(x12z,-a1);  //Z=19190
                //F62as1sumz = pva0;  //Z=19191
                //F62as1z = F62as10z*F62as1sumz;  //Z=19192
                F62as1z0 = F62as10z*pza;     /* * */  //Z=19193

                F62as40z = preg3*pow(x1z,c);  //Z=19195
                arg64 = (zr+c+1)*atan(2.0*x1z);  //Z=19196
                nen64 = pow(1.0+4*x1z*x1z,(zr+c+1)/2.0);  //Z=19197
                arg65 = (zr+c)*atan(2.0*x1z);  //Z=19198
                nen65 = pow(1.0+4*x1z*x1z,(zr+c)/2.0);  //Z=19199
                F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);  //Z=19200
                F62as28 = e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);  //Z=19201
                F62as4z = F62as40z*(F62as27+F62as28);  //Z=19202
                //F62asz = F62as1z+F62as4z;  //Z=19203
                F62asz0 = F62as1z0+F62as4z;  //Z=19204
                F62 = F62asz0;  //Z=19205
            }/*4*/  //Z=19206

            FF1 = (cc1*F12+cc4*F42+cc6*F62)/sumc;  //Z=19208
            /* FF1:=(cc6*F62)/sumc;  //Z=19209 */


            /*formfq:=*/ return FF1*FF1;  //Z=19212

            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,2,sigmar,q);  //Z=19214 */
        }/*3*/ /*  of inhomogeneous core/shell cylinder */  //Z=19215


    }/*2*/ /*  of cylinder  */  //Z=19218



    /* ****** */  //Z=19222
    /*  disk  */  //Z=19223
    /* ****** */  //Z=19224
    if ( params.part==2 )
    {/*2*/  //Z=19225

        /* ** longitudinal part ** */  //Z=19227
        /* ** isotropic ** */  //Z=19228
        if ( ordis==7 )
        {/*3*/  //Z=19229
            if ( q<(0.5*params.limq1f) )
            {/*4*/  //Z=19230
                pqsum = 1.0;  //Z=19231
                oldpqsum = 0.0;  //Z=19232
                double qqnn = 1.0;  //Z=19233
                for ( nser=1; nser<=80; nser++ )
                {/*5*/  //Z=19234
                    qqnn = qqnn*q*q;  //Z=19235
                    pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=19236
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19237
                    if ( delser<0.0001 ) break; /* goto 70; */  //Z=19238
                    oldpqsum = pqsum;  //Z=19239
                }/*5*/  //Z=19240
                /*70:*/  //Z=19241
                pql = pqsum;  //Z=19242
            }/*4*/  //Z=19243
            else
            {/*4*/  /*   = P||(q)   */  //Z=19244
                arglq = q*params.length/(zl+1);  //Z=19245
                pql = (2/(zl*(zl-1)))*pow(arglq,-2);  //Z=19246
            }/*4*/  //Z=19247
        }/*3*/  /*  of isotropic  */  //Z=19248

        /*  perfect  */  //Z=19250
        if ( ordis==6 )
        {/*3*/  //Z=19251
            if ( q<(0.5*params.limq1f) )
            {/*4*/  //Z=19252
                pqsum = 1.0;  //Z=19253
                oldpqsum = 0.0;  //Z=19254
                double qqnn = 1.0;  //Z=19255
                if ( params.orcase==1 )
                {/*5*/  //Z=19256
                    argq = qxs+qys;  //Z=19257
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19258
                        qqnn = qqnn*q*q;  //Z=19259
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn*pow(1.0-argq*argq,nser);  //Z=19260
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19261
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19262
                        oldpqsum = pqsum;  //Z=19263
                    }/*6*/  //Z=19264
                }/*5*/  //Z=19265
                if ( params.orcase==2 )
                {/*5*/  //Z=19266
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19267
                        qqnn = qqnn*q*q;  //Z=19268
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn*pow(1.0-qxs*qxs,nser);  //Z=19269
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19270
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19271
                        oldpqsum = pqsum;  //Z=19272
                    }/*6*/  //Z=19273
                }/*5*/  //Z=19274
                if ( params.orcase==3 )
                {/*5*/  //Z=19275
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19276
                        qqnn = qqnn*q*q;  //Z=19277
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn*pow(1.0-qys*qys,nser);  //Z=19278
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19279
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19280
                        oldpqsum = pqsum;  //Z=19281
                    }/*6*/  //Z=19282
                }/*5*/  //Z=19283
                if ( params.orcase==4 )
                {/*5*/  //Z=19284
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19285
                        qqnn = qqnn*q*q;  //Z=19286
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=19287
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19288
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19289
                        oldpqsum = pqsum;  //Z=19290
                    }/*6*/  //Z=19291
                }/*5*/  //Z=19292
                /*76:*/  //Z=19293
                pql = pqsum;  //Z=19294
            }/*4*/  //Z=19295
            else
            {/*4*/      /*  ok  */  //Z=19296
                if ( params.orcase==1 )
                {/*5*/  //Z=19297
                    qnarg = qxs+qys;  //Z=19298
                    arglq = sqrt(1.0-qnarg*qnarg)*q*params.length/(zl+1)+eps9;  //Z=19299
                }/*5*/  //Z=19300
                if ( params.orcase==2 ) arglq = sqrt(1.0-qxs*qxs)*q*params.length/(zl+1)+eps9;  //Z=19301
                if ( params.orcase==3 ) arglq = sqrt(1.0-qys*qys)*q*params.length/(zl+1)+eps9;  //Z=19302
                if ( params.orcase==4 ) arglq = q*params.length/(zl+1)+eps9;  //Z=19303

                /*  F(q)  */  //Z=19305
                /* pqr1:=(gamma(zl-1/2)/gamma(zl+1))*power(arglq,-3/2)*(sin((zl-1/2)*arctan(arglq))-cos((zl-1/2)*arctan(arglq)))/power(1+arglq*arglq,(zl-1/2)/2);  //Z=19306 */
                /* pqr2:=(gamma(zl-3/2)/gamma(zl+1))*power(arglq,-5/2)*(sin((zl-3/2)*arctan(arglq))+cos((zl-3/2)*arctan(arglq)))/power(1+arglq*arglq,(zl-3/2)/2);  //Z=19307 */
                /* pqr3:=(2/sqrt(pi))*(pqr1+(9/16)*pqr2);  //Z=19308 */
                /* pql:=pqr3*pqr3;  //Z=19309 */

                /*  P(q)  */  //Z=19311
                pqr1 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3);  //Z=19312
                pqr2 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3)*sin((zl-2)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-2)/2.0);  //Z=19313
                pqr3 = (1/(zl*(zl-1)*(zl-2)*(zl-3)))*pow(arglq,-4)*cos((zl-3)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-3)/2.0);  //Z=19314
                pql = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=19315
            }/*4*/  //Z=19316
        }/*3*/   /*  of perfect  */  //Z=19317

        /*  orientational distribution  */  //Z=19319
        if ( ordis==0 )
        {/*3*/  //Z=19320
            if ( params.orcase==1 )
            {/*4*/  //Z=19321
                if ( q<(1.2*params.limq1f) )  // war 1.2
                {/*5*/  //Z=19322
                    pqsum = 1.0;  //Z=19323
                    oldpqsum = 0.0;  //Z=19324
                    double qqnn = 1.0;  //Z=19325
                    qxn[0] = 1.0;  //Z=19326
                    qyn[0] = 1.0;  //Z=19327

                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19329
                        qqnn = qqnn*q*q;  //Z=19330
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19331
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19332

                        binsum = 0.0;  //Z=19334
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19335
                            binsum1 = 0.0;  //Z=19336
                            for ( lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=19337
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=19338 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=19339 */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=19340
                            }/*8*/  //Z=19341
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=19342 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=19343 */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;  //Z=19344
                        }/*7*/  //Z=19345
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn*binsum;  //Z=19346
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19347
                        if ( delser<0.0001 ) break; /* goto 77; */  //Z=19348
                        oldpqsum = pqsum;  //Z=19349
                    }/*6*/  //Z=19350
                    /*77:*/  //Z=19351
                    pql = pqsum;  //Z=19352
                }/*5*/  //Z=19353
                else
                {/*5*/  //Z=19354
                    /*  disk: length = disk radius  */  //Z=19355
                    /*  always use Bessel function approximation  */  //Z=19356
                    /*  F(q)  */  //Z=19357
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19358 */
                    /*  P(q)  */  //Z=19359
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,
                                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19360
                    pql = pql/params.norm;  //Z=19361
                }/*5*/  //Z=19362
            }/*4*/  //Z=19363

            if ( params.orcase==2 )
            {/*4*/  //Z=19365
                if ( q<(0.9*params.limq1f) )
                {/*5*/  //Z=19366
                    pqsum = 1.0;  //Z=19367
                    oldpqsum = 0.0;  //Z=19368
                    double qqnn = 1.0;  //Z=19369
                    qxn[0] = 1.0;  //Z=19370
                    qyn[0] = 1.0;  //Z=19371

                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19373
                        qqnn = qqnn*q*q;  //Z=19374
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19375
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19376

                        binsum = 0.0;  //Z=19378
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19379
                            binsum1 = 0.0;  //Z=19380
                            for ( lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=19381
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=19382 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=19383 */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=19384
                            }/*8*/  //Z=19385
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=19386 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=19387 */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;  //Z=19388
                        }/*7*/  //Z=19389
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn*binsum;  //Z=19390
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19391
                        if ( delser<0.0001 ) break; /* goto 78; */  //Z=19392
                        oldpqsum = pqsum;  //Z=19393
                    }/*6*/  //Z=19394
                    /*78:*/  //Z=19395
                    pql = pqsum;  //Z=19396
                }/*5*/  //Z=19397
                else
                {/*5*/  //Z=19398
                    /*  disk: length = disk radius  */  //Z=19399
                    /*  always use Bessel function approximation  */  //Z=19400
                        /*  F(q)  */  //Z=19401
                        /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19402 */
                        /*  P(q)  */  //Z=19403
                        qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,
                                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19404
                        pql = pql/params.norm;  //Z=19405
                }/*5*/  //Z=19406
            }/*4*/  //Z=19407

            if ( params.orcase==3 )
            {/*4*/  //Z=19409
                if ( q<(0.9*params.limq1f) )
                {/*5*/  //Z=19410
                    pqsum = 1.0;  //Z=19411
                    oldpqsum = 0.0;  //Z=19412
                    double qqnn = 1.0;  //Z=19413
                    qxn[0] = 1.0;  //Z=19414
                    qyn[0] = 1.0;  //Z=19415

                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19417
                        qqnn = qqnn*q*q;  //Z=19418
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19419
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19420

                        binsum = 0.0;  //Z=19422
                        for ( mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19423
                            binsum1 = 0.0;  //Z=19424
                            for ( lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=19425
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=19426 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=19427 */
                                binsum1 = binsum1+params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=19428
                            }/*8*/  //Z=19429
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=19430 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=19431 */
                            binsum = binsum+params.CR->carr11pm[nser][mser]*binsum1;  //Z=19432
                        }/*7*/  //Z=19433
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn*binsum;  //Z=19434
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19435
                        if ( delser<0.0001 ) break; /* goto 79; */  //Z=19436
                        oldpqsum = pqsum;  //Z=19437
                    }/*6*/  //Z=19438
                    /*79:*/  //Z=19439
                    pql = pqsum;  //Z=19440
                }/*5*/  //Z=19441
                else
                {/*5*/  //Z=19442
                    /*  disk: length = disk radius  */  //Z=19443
                    /*  always use Bessel function approximation  */  //Z=19444
                    /*  F(q)  */  //Z=19445
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19446 */
                    /*  P(q)  */  //Z=19447
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,
                                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19448
                    pql = pql/params.norm;  //Z=19449
                }/*5*/  //Z=19450
            }/*4*/  //Z=19451

            if ( params.orcase==4 )
            {/*4*/  //Z=19453
                if ( q<(0.5*params.limq1f) )
                {/*5*/  //Z=19454
                    pqsum = 1.0;  //Z=19455
                    oldpqsum = 0.0;  //Z=19456
                    double qqnn = 1.0;  //Z=19457
                    for ( nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19458
                        qqnn = qqnn*q*q;  //Z=19459
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=19460
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19461
                        if ( delser<0.0001 ) break; /* goto 80; */  //Z=19462
                        oldpqsum = pqsum;  //Z=19463
                    }/*6*/  //Z=19464
                    /*80:*/  //Z=19465
                    pql = pqsum;  //Z=19466
                }/*5*/  //Z=19467
                else
                {/*5*/  //Z=19468
                    /*  disk: length = disk radius  */  //Z=19469
                    /*  always use Bessel function approximation  */  //Z=19470
                    /*  F(q)  */  //Z=19471
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19472 */
                    /*  P(q)  */  //Z=19473
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,
                                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19474
                    pql = pql/params.norm;  //Z=19475
                }/*5*/  //Z=19476
            }/*4*/  //Z=19477
        }/*3*/   /*  of orientational distribution  */  //Z=19478

        /*  transverse part  */  //Z=19480
        /*  disk: radius = disk thickness/2  */  //Z=19481
        /*  homogeneous disk  */  //Z=19482
        if ( params.cs==0 )
        {/*3*/  //Z=19483
            if ( q<(0.5*params.limq4f) )
            {/*4*/  //Z=19484
                pqsum = 1.0;  //Z=19485
                oldpqsum = 0.0;  //Z=19486
                double qqnn = 1.0;  //Z=19487
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19488
                    qqnn = qqnn*q*q;  //Z=19489
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=19490
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19491
                    if ( delser<0.0001 ) break; /* goto 71; */  //Z=19492
                    oldpqsum = pqsum;  //Z=19493
                }/*5*/  //Z=19494
                /*71:*/  //Z=19495
                pqr = pqsum;  //Z=19496
            }/*4*/  //Z=19497
            else
            {/*4*/  //Z=19498
                argpq = q*params.radius/(zr+1);  //Z=19499
                pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1.0+argpq*argpq,zr/2.0);  //Z=19500
            }/*4*/  //Z=19501
            /*formfq:=*/ return pql*pqr*pqr;  //Z=19502
            /* formfq:=pql;;  //Z=19503 */
        }/*3*/ /*  of homogeneous  */  //Z=19504

        /*  core/shell disk  */  //Z=19506
        if ( params.cs==1 )
        {/*3*/  //Z=19507
            ccc1 = sqr(1-params.rho)*pow(params.p1,2);  //Z=19508
            ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,1);  //Z=19509
            ccc3 = sqr(params.rho);  //Z=19510
            vv3 = sqr((1-params.rho)*pow(params.p1,1)+params.rho);  //Z=19511

            zz = zr;  // TODO: zz war in diesem Zweig nicht gesetzt
            argq = q*radiusm/(zz+1);  //Z=19513
            argpq = q*params.radius/(zz+1);  //Z=19514

            /*  F121 disk  */  //Z=19516
            if ( q<(0.8*params.limq4f) )
            {/*4*/  //Z=19517
                /* ** series expansion ** */  //Z=19518
                pqsum = 1.0;  //Z=19519
                oldpqsum = 0.0;  //Z=19520
                double qqnn = 1.0;  //Z=19521
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19522
                    qqnn = qqnn*q*q;  //Z=19523
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=19524
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19525
                    if ( delser<0.0001 ) break; /* goto 72; */  //Z=19526
                    oldpqsum = pqsum;  //Z=19527
                }/*5*/  //Z=19528
                /*72:*/  //Z=19529
                F121 = ccc1*pqsum/vv3;  //Z=19530
            }/*4*/  //Z=19531
            else
            {/*4*/  //Z=19532
                pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1.0+argpq*argpq,zr/2.0);  //Z=19533
                F121 = ccc1*pqr*pqr/vv3;  //Z=19534
            }/*4*/  //Z=19535

            /*  F122 disk  */  //Z=19537
            if ( q<(1.0*params.limq5f) )
            {/*4*/  //Z=19538
                /* ** series expansion ** */  //Z=19539
                pqsum = 1.0;  //Z=19540
                oldpqsum = 0.0;  //Z=19541
                double qqnn = 1.0;  //Z=19542
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19543
                    qqnn = qqnn*q*q;  //Z=19544
                    pqsum = pqsum+params.CR->carr5f[nser]*qqnn;  //Z=19545
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19546
                    if ( delser<0.0001 ) break; /* goto 73; */  //Z=19547
                    oldpqsum = pqsum;  //Z=19548
                }/*5*/  //Z=19549
                /*73:*/  //Z=19550
                F122 = ccc2*pqsum/vv3;  //Z=19551
            }/*4*/  //Z=19552
            else
            {/*4*/  //Z=19553
                pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1.0+argpq*argpq,zr/2.0);  //Z=19554
                pqr1 = (1/zr)*pow(argq,-1)*sin(zr*atan(argq))/pow(1.0+argq*argq,zr/2.0);  //Z=19555
                F122 = ccc2*pqr*pqr1/vv3;  //Z=19556
            }/*4*/  //Z=19557

            /*  F123 disk  */  //Z=19559
            if ( q<(0.3*params.limq6f) )
            {/*4*/  //Z=19560
                /* ** series expansion ** */  //Z=19561
                pqsum = 1.0;  //Z=19562
                oldpqsum = 0.0;  //Z=19563
                double qqnn = 1.0;  //Z=19564
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19565
                    qqnn = qqnn*q*q;  //Z=19566
                    pqsum = pqsum+params.CR->carr6f[nser]*qqnn;  //Z=19567
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19568
                    if ( delser<0.0001 ) break; /* goto 74; */  //Z=19569
                    oldpqsum = pqsum;  //Z=19570
                }/*5*/  //Z=19571
                /*74:*/  //Z=19572
                F123 = ccc3*pqsum/vv3;  //Z=19573
            }/*4*/  //Z=19574
            else
            {/*4*/  //Z=19575
                pqr = (1/zr)*pow(argq,-1)*sin(zr*atan(argq))/pow(1.0+argq*argq,zr/2.0);  //Z=19576
                F123 = ccc3*pqr*pqr/vv3;  //Z=19577
                /*  add more terms, if necessary  */  //Z=19578
            }/*4*/  //Z=19579
            /*formfq:=*/ return pql*(F121+F122+F123);  //Z=19580
            /* formfq:=(F121+F122+F123);  //Z=19581 */
        }/*3*/ /*  of core/shell-disk  */  //Z=19582

        /* ** inhomogeneous core/shell disk ** */  //Z=19584
        if ( params.cs==2 )
        {/*3*/  //Z=19585

            dim = 1;  //Z=19587
            delc = 0.0001;  //Z=19588
            zz = zr;  //Z=19589
            xrad = q*radiusm;  //Z=19590
            xradp = q*params.radius;  //Z=19591
            x1z = q*params.radius/(2.0*(zz+1));  //Z=19592
            x12z = x1z*x1z;  //Z=19593
            x2z = q*radiusm/(2.0*(zz+1));  //Z=19594
            x22z = x2z*x2z;  //Z=19595

            lim = 18*exp(-5*params.sigma);  //Z=19597
            lim1 = lim;  //Z=19598
            //lim2 = lim*0.7;  //Z=19599
            //lim3 = lim;  //Z=19600
            lim4 = lim;  //Z=19601
            //lim5 = lim*0.7;  //Z=19602
            lim6 = lim*1.2;  //Z=19603

            a1 = (dim-params.alphash1)/2.0;  //Z=19605
            b1 = dim/2.0;  //Z=19606
            b2 = (dim+2-params.alphash1)/2.0;  //Z=19607
            b1s = (dim+2)/2.0;  //Z=19608
            v = -b1s+1/2.0;  //Z=19609
            c = a1-b1-b2+1/2.0;  //Z=19610
            //d0 = 1;  //Z=19611
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=19612
            e0 = 1.0;  //Z=19613
            e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=19614
            ee0 = 1.0;  //Z=19615
            ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=19616

            gb1s = 3*sqrt(M_PI)/4.0;  //Z=19618
            pz2v = 1/(zr*(zr-1));  //Z=19619
            pz2v1 = pz2v/(zr-2);  //Z=19620
            //pz2v2 = pz2v1/(zr-3);  //Z=19621

            gz1 = gamma(zr+1);  //Z=19623
            preg1 = gb1s/sqrt(M_PI);  //Z=19624
            preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=19625
            preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=19626
            //pzvc = gamma(zr+1+v+c)/gz1;  //Z=19627
            //pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=19628
            //pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=19629
            //pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=19630
            //pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=19631
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=19632
            pzc = gamma(zr+1+c)/gz1;  //Z=19633
            pzc1 = gamma(zr+1+c-1)/gz1;  //Z=19634
            pza = gamma(zr+1-2*a1)/gz1;  //Z=19635
            //pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=19636
            //pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19637
            //dnv0 = 1;  //Z=19638
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=19639
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19640
            //pva0 = gamma(zr+1-2*a1)/gz1;  //Z=19641

            cc1 = 1/dim;  //Z=19643
            cc4 = params.rho/((dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=19644
            cc6 = -params.rho/(dim-params.alphash1);  //Z=19645
            sumc = cc1+cc4+cc6;  //Z=19646

            /*  term #1 series  */  //Z=19648
            if ( (xradp)<lim1 )
            {/*4*/  //Z=19649
                z12v[0] = 1;  //Z=19650
                b1sv[0] = 1;  //Z=19651
                fkv[0] = 1;  //Z=19652
                double qqnn = 1.0;  //Z=19653
                F12sez = 1.0;  //Z=19654
                oldF12sez = 1.0;  //Z=19655
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=19656
                    qqnn = qqnn*q*q;  //Z=19657
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19658
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19659
                    fkv[n] = fkv[n-1]*n;  //Z=19660
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]);  //Z=19661 */

                    F12sez = F12sez+params.CR->carr4f[n]*qqnn;  //Z=19663

                    del = fabs((F12sez-oldF12sez)/F12sez);  //Z=19665
                    if ( del<delc ) break; /* goto 121; */  //Z=19666
                    oldF12sez = F12sez;  //Z=19667
                }/*5*/  //Z=19668
                /*121:*/  //Z=19669
                F12 = F12sez;  //Z=19670
            }/*4*/  //Z=19671

            /*  term #4 series  */  //Z=19673
            if ( xradp<lim4 )
            {/*4*/  //Z=19674
                z12v[0] = 1;  //Z=19675
                a1v[0] = 1;  //Z=19676
                b1v[0] = 1;  //Z=19677
                b2v[0] = 1;  //Z=19678
                b1sv[0] = 1;  //Z=19679
                fkv[0] = 1;  //Z=19680
                double qqnn = 1.0;  //Z=19681
                F42sez = 1.0;  //Z=19682
                oldF42sez = 1.0;  //Z=19683
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=19684
                    qqnn = qqnn*q*q;  //Z=19685
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19686
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19687
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19688
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19689
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19690
                    fkv[n] = fkv[n-1]*n;  //Z=19691
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19692 */

                    F42sez = F42sez+params.CR->carr5f[n]*qqnn;  //Z=19694

                    del = fabs((F42sez-oldF42sez)/F42sez);  //Z=19696
                    if ( del<delc ) break; /* goto 124; */  //Z=19697
                    oldF42sez = F42sez;  //Z=19698
                }/*5*/  //Z=19699
                /*124:*/  //Z=19700
                F42 = F42sez;  //Z=19701
            }/*4*/  //Z=19702

            /*  term #6 series  */  //Z=19704
            if ( xradp<lim6 )
            {/*4*/  //Z=19705
                z12v[0] = 1;  //Z=19706
                a1v[0] = 1;  //Z=19707
                b1v[0] = 1;  //Z=19708
                b2v[0] = 1;  //Z=19709
                b1sv[0] = 1;  //Z=19710
                fkv[0] = 1;  //Z=19711
                double qqnn = 1.0;  //Z=19712
                F62sez = 1.0;  //Z=19713
                oldF62sez = 1.0;  //Z=19714
                for ( n=1; n<=120; n++ )
                {/*5*/  //Z=19715
                    qqnn = qqnn*q*q;  //Z=19716
                    z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19717
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19718
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19719
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19720
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19721
                    fkv[n] = fkv[n-1]*n;  //Z=19722
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19723 */

                    F62sez = F62sez+params.CR->carr6f[n]*qqnn;  //Z=19725

                    del = fabs((F62sez-oldF62sez)/F62sez);  //Z=19727
                    if ( del<delc ) break; /* goto 126; */  //Z=19728
                    oldF62sez = F62sez;  //Z=19729
                }/*5*/  //Z=19730
                /*126:*/  //Z=19731
                F62 = F62sez;  //Z=19732
            }/*4*/  //Z=19733

            /* ** term #1 asymptote ** */  //Z=19735
            if ( xradp>=lim1 )
            {/*4*/  //Z=19736
                arg11 = (zr+v+1)*atan(2.0*x1z);  //Z=19737
                nen11 = pow(1.0+4*x1z*x1z,(zr+v+1)/2.0);  //Z=19738
                arg12 = (zr+v)*atan(2.0*x1z);  //Z=19739
                nen12 = pow(1.0+4*x1z*x1z,(zr+v)/2.0);  //Z=19740

                F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);  //Z=19742
                F12as2z = ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);  //Z=19743
                F12asz = preg1*pow(x1z,v)*(F12as1z+F12as2z);  //Z=19744
                F12 = F12asz;  //Z=19745
            }/*4*/  //Z=19746

            /* ** term #4 asymptote ** */  //Z=19748
            if ( xrad>=lim4 )
            {/*4*/  //Z=19749
                F42as10z = preg4*pow(x22z,-a1);  //Z=19750
                //F42as1sumz = pva0;  //Z=19751
                //F42as1z = F42as10z*F42as1sumz;  //Z=19752
                F42as1z0 = F42as10z*pza;   /* * */  //Z=19753

                F42as40z = preg3*pow(x2z,c);  //Z=19755
                arg44 = (zr+c+1)*atan(2.0*x2z);  //Z=19756
                nen44 = pow(1.0+4*x2z*x2z,(zr+c+1)/2.0);  //Z=19757
                arg45 = (zr+c)*atan(2.0*x2z);  //Z=19758
                nen45 = pow(1.0+4*x2z*x2z,(zr+c)/2.0);  //Z=19759
                F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);  //Z=19760
                F42as28 = e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);  //Z=19761
                F42as4z = F42as40z*(F42as27+F42as28);  //Z=19762
                //F42asz = F42as1z+F42as4z;  //Z=19763
                F42asz0 = F42as1z0+F42as4z;  //Z=19764
                F42 = F42asz0;  //Z=19765
            }/*4*/  //Z=19766

            /* ** term #6 asymptote ** */  //Z=19768
            if ( xradp>=lim6 )
            {/*4*/  //Z=19769
                F62as10z = preg4*pow(x12z,-a1);  //Z=19770
                //F62as1sumz = pva0;  //Z=19771
                //F62as1z = F62as10z*F62as1sumz;  //Z=19772
                F62as1z0 = F62as10z*pza;     /* * */  //Z=19773

                F62as40z = preg3*pow(x1z,c);  //Z=19775
                arg64 = (zr+c+1)*atan(2.0*x1z);  //Z=19776
                nen64 = pow(1.0+4*x1z*x1z,(zr+c+1)/2.0);  //Z=19777
                arg65 = (zr+c)*atan(2.0*x1z);  //Z=19778
                nen65 = pow(1.0+4*x1z*x1z,(zr+c)/2.0);  //Z=19779
                F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);  //Z=19780
                F62as28 = e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);  //Z=19781
                F62as4z = F62as40z*(F62as27+F62as28);  //Z=19782
                //F62asz = F62as1z+F62as4z;  //Z=19783
                F62asz0 = F62as1z0+F62as4z;  //Z=19784
                F62 = F62asz0;  //Z=19785
            }/*4*/  //Z=19786

            FF1 = (cc1*F12+cc4*F42+cc6*F62)/sumc;  //Z=19788
            /* FF1:=(cc1*F12)/sumc;  //Z=19789 */

            /*formfq:=*/ return FF1*FF1;  //Z=19791

            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,1,sigmar,q);  //Z=19793 */

        }/*3*/ /*  of inhomogeneous core/shell disk  */  //Z=19795
    }/*2*/ /*  of disk  */  //Z=19796


    /*  cube  */  //Z=19799
    if ( params.part==5 )
    {/*2*/  //Z=19800
        /*  homogeneous cube  */  //Z=19801
        if ( params.cs==0 )
        {/*3*/  //Z=19802
            if ( q<0.7*params.limq4f )
            {/*4*/  //Z=19803
                pqsum = 1.0;  //Z=19804
                oldpqsum = 0.0;  //Z=19805
                double qqnn = 1.0;  //Z=19806
                for ( nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19807
                    qqnn = qqnn*q*q;  //Z=19808
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=19809
                    delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19810
                    if ( delser<0.0001 ) break; /* goto 81; */  //Z=19811
                    oldpqsum = pqsum;  //Z=19812
                }/*5*/  //Z=19813
                /*81:*/  //Z=19814
                /*formfq:=*/ return pqsum;  //Z=19815
            }/*4*/  //Z=19816
            else
            {/*4*/  //Z=19817
                qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,qxhklt,qyhklt,qzhklt,qhkl,
                            params.ax1.length(),params.ax2.length(),params.ax3.length(),
                            params.ax1.x(),params.ax1.y(),params.ax1.z(),
                            params.ax2.x(),params.ax2.y(),params.ax2.z(),
                            params.ax3.x(),params.ax3.y(),params.ax3.z(),
                            params.sig.x(),params.sig.y(),params.sig.z(),
                            ordis,3,7,12,7,1,0,params.CR->carr1f,pql);  //Z=19818
                    /*formfq:=*/ return pql/(M_PI/2.0);  //Z=19819
                }/*4*/  //Z=19820
            }/*3*/ /*  of homogeneous cube */  //Z=19821

            /*  core/shell cube  */  //Z=19823
            if ( params.cs==1 )
            {/*3*/  //Z=19824
                /*formfq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=19825
            }/*3*/  //Z=19826

    }/*2*/  /*  of cube  */  //Z=19828

    return 0.0;
}  //Z=19829



#ifdef __CUDACC__
__host__ __device__
#endif
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


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::polycscube(double /*rho1*/, double /*rho2*/, double /*p1*/, double /*p2*/,
                                            double /*alf1*/, double /*alf2*/, double /*rn*/, double /*pf*/,
                                            double /*sigma*/, double /*q*/) const
{
    // Gespräch mit Prof. Förster (05.Jun.2023): Diese Routinen werden noch nicht verwendet, sie werden bei Erweiterungen kommen.
    return 0;
}
