
#ifdef __CUDACC__
__global__ void kernel_partSphere_peakGauss_lattFCC_ordisGauss( SasCalc_GENERIC_calculation CALC, bool dofit )
{
    KERNEL_PREWORK
    const double pixval = CALC.calc_partSphere_peakGauss_lattFCC_ordisGauss( CALC, qx, qy, qz );
    KERNEL_POSTWORK
}
#endif


#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::calc_partSphere_peakGauss_lattFCC_ordisGauss(const SasCalc_GENERIC_calculation& CALC,
                                                                     double qx, double qy, double qz)
{
    const double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254

    /* ************* */  //Z=25269  //ZN=25476
    /* ** spheres ** */  //Z=25270
    /* ************* */  //Z=25271
    //    if ( CALC.ComboBoxParticle == cbpartSphere /*CALC.partsphere*/ )
    //    {/*6*/  //Z=25272

    const double pq = CALC.formpq_partSphere(q);  //Z=25276

    /* fq = pq;  //Z=25277 */
    //if ( CALC.lattice )  //Z=25283
    const double fq=CALC.formfq_partSphere( q );  //Z=25285
    const double pqiso = pq;  //Z=25297

    //    }/*6*/  /*  of part=0  */  //Z=25299

    double width = 1;
    if ( CALC.RadioButtonDebyeScherrer )
        width = 4.0/CALC.params.domainsize;  //Z=25722
    if ( CALC.RadioButtonPara )
        width = (4.0/CALC.params.domainsize)+sqr(CALC.reldis)*CALC.dist*sqr(q);  //Z=25723
    //if ( CALC.lattice && CALC.shp==cbpeakAnisotropicGaussian/*8*/ && CALC.ordis==7/*isotropic*/ )
    //    width = CALC.params.sig.length() /*sqrt(sigx*sigx+sigy*sigy+sigz*sigz)*/ /3.0;  //Z=26118

    //if ( ! CALC.tpvRandomOldValues.empty() )   // doIntCalc... - nicht auf der GPU zulässig
    //{
    //if ( fabs(qx) < 0.1 && fabs(qy) < 0.1 && fabs(qz) < 0.1 )
    //    std::cerr << "TPV CALC " << CALC.params.width_zuf << std::endl << std::flush;
    width *= CALC.params.width_zuf;
    //}

    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

    double psiord;

    const double widthiso = 1.0 / CALC.params.uca;  //ZN=20321

    double radintensity = 0.0;     // Immer nur für ein Pixel berechnet, somit kein Array nötig
    double intensity = 0.0;

    /* ** lattice hkl-factor calcuation  */  //Z=25744
    //if ( CALC.lattice )
    //{/*6*/  //Z=25745

    //const double sphno = CALC.latpar[1];  //Z=25749
    //const double cubevol = CALC.latpar[2];  //Z=25750
    //global: dwfactor = CALC.latpar[3];  //Z=25751
    const double sphno_cubevol = CALC.latpar[1] * CALC.latpar[2]; // werden immer so genutzt

    /* isotropic peaks */
    for ( int ii=1; ii<=CALC.peakmax1; ii++ )
    {   //Z=25759  //ZN=25967
        if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat
        //int h    = CALC.latpar1(ii,1);  //Z=25760
        //int k    = CALC.latpar1(ii,2);
        //int l    = CALC.latpar1(ii,3);
        const int mhkl    = CALC.latpar1(ii,4);  //Z=25763
        const int fhkl    = CALC.latpar1(ii,5);
        if ( mhkl==0 || fhkl==0 ) continue;
        const double qhkl = CALC.latpar3(ii,5);
        const double x2   = sqr(q-qhkl)/sqr(widthiso);  //Z=25766
        const double sq   = exp(-4*x2/M_PI)/(M_PI*widthiso/2.0);  //Z=25767

        radintensity += sq*mhkl*fhkl;  //Z=25768
    } /* of peak loop */  //Z=25769

    for ( int ii=1; ii<=CALC.peakmax2; ii++ )
    {   //Z=25779  //ZN=25987
        if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat
        const double qhkl   = CALC.latpar3(ii,5);
        //const int mhkl = CALC.latpar2(ii,4);
        //const int fhkl = CALC.latpar2(ii,5);
        const int mhkl_fhkl = CALC.latpar2(ii,4) * CALC.latpar2(ii,5); // werden nur so genutzt
        if ( qhkl > 0 && mhkl_fhkl != 0 )
        {
            //int h = CALC.latpar2(ii,1);  //Z=25780  //ZN=25988
            //int k = CALC.latpar2(ii,2);
            //int l = CALC.latpar2(ii,3);

            //qhkl0 = CALC.latpar3(ii,1);  //Z=25820  //ZN=26028
            const double qxhkl  = CALC.latpar3(ii,2);
            const double qyhkl  = CALC.latpar3(ii,3);
            const double qzhkl  = CALC.latpar3(ii,4);
            const double g3     = CALC.latpar3(ii,10);

            //switch ( CALC.shp )
            //{
            //case cbpeakGaussian /*2*/:                                                          //20210812-E
            const double peaknorm1 = CALC.latpar3(ii,11);  //Z=26081
            const double yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
            const double x2 = (q-qhkl)*(q-qhkl)/(width*width);
            const double sq = exp(-4*x2/M_PI)/(M_PI*width/2.0);
            const double x2phi = 4*q*q/(M_PI*sqr(CALC.phiwidth));             /*** a-factor ***/
            psiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
            if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
            {
                const double qxhklt = CALC.latpar3(ii,7);
                const double qyhklt = CALC.latpar3(ii,8);
                const double qzhklt = CALC.latpar3(ii,9);
                /*//20240301 raus:*/ const double qhklt  = CALC.latpar3(ii,14);  //ZN=26037
                /*//20240301 raus:*/ const double g3t    = CALC.latpar3(ii,15);  //ZN=26038
                /*//20240301 raus:*/ const double peaknorm1t = CALC.latpar3(ii,16);   //ZN=26305
                const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhklt));
                const double phiord = g3t*exp(-x2phi*yphi*yphi)/peaknorm1t;
                psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;   //ZN=26313
            }
            //default:
            //    return 0;
            //} // switch shp

            //if ( qhkl > 0 )
            //{ --> kann weiter oben geprüft werden und spart u.U. viele Berechnungen
            const double savintens = intensity;
            intensity += sq*mhkl_fhkl*psiord;  //Z=26176
            if ( isinf(intensity) )
            {
                intensity = savintens; // wieder einen Schritt zurück gehen...
                //#ifndef __CUDACC__
                //    qDebug() << "INTENS=INF" << sq << mhkl << fhkl << psiord << "ii" << ii << CALC.peakmax2;
                //#endif
                break;
            }
        }

    }/*2*/  /* of peak-loop */  //Z=26182

    const double szqiso = (1+(1)*(8*M_PI*M_PI*M_PI*radintensity/(4*M_PI*sphno_cubevol)-1)*exp(-CALC.params.ucb/*dwfactoriso*/*q*q));   /*Z0311=24804*/
    const double szq = (1+(fq/pq)*(2*M_PI*intensity/sphno_cubevol-1)*exp(-CALC.dwfactor*q*q));   /*Z0311=24805*/
    //}
    //else // (if lattice)
    //{  //Z=26200
    //    szq = 1.0;
    //    szqiso = 1.0;   // zur Sicherheit
    //}

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    const double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    // retval ist der Pixelwert bei [ihex+zmax][i+zmax]

/*#ifndef __CUDACC__
    if ( retval < -1e8 || isnan(retval) || isinf(retval) )
        //qDebug() << "szq"<<szq << "pq"<<pq <<"fq"<<fq << "szqiso"<<szqiso << "pqiso"<<pqiso << "q"<<q
                 << "intens"<<intensity << "radint"<<radintensity << CALC.peakmax1
                 << "="<<retval;
        // szq inf pq 2.16215e-05 fq 2.08036e-07 szqiso 1 pqiso 2.16215e-05 q 2.74004 intens inf radint 2.45523 107 = inf

#else
    //if ( retval < -1e6 )
    //    printf( "szq=%lf pq=%lf fq=%lf szqiso=%lf pqiso=%lf q=%lf intens=%lf radint=%lf erg=%lf\n", szq, pq, fq, szqiso, pqiso, q, intensity, radintensity, retval );
#endif*/

    return retval;
} /* calc_partSphere_peakGauss_lattFCC_ordisGauss() */
