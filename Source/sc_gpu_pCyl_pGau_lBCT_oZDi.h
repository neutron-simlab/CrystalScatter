
#ifdef __CUDACC__
// Default kernel if no specific kernel is defined
__global__ void kernel_partCylinder_peakGauss_lattBCT_ordisZDir( SasCalc_GENERIC_calculation CALC, bool dofit )
{
    KERNEL_PREWORK
    const double pixval = CALC.calc_partCylinder_peakGauss_lattBCT_ordisZDir( CALC, qx, qy, qz );
    KERNEL_POSTWORK
}
#endif


#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::calc_partCylinder_peakGauss_lattBCT_ordisZDir(const SasCalc_GENERIC_calculation& CALC,
                                                                   double qx, double qy, double qz)
{
    double pq=0, fq=0, intensity, radintensity;

    const double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254

    //const double delta = 2.0 * CALC.params.radius;
    double pqiso=1;
    double szqiso, szq;

    /* *************** */  //Z=25302
    /* ** cylinders ** */  //Z=25303
    /* *************** */  //Z=25304
    //if ( CALC.partcylinder )
    //{/*6*/     /*  cylinders  */  //Z=25305

    /*  perfect orientation  */               /*  abs((cosphi*qx+sinphi*qy)  //Z=25318 */
    //if ( CALC.params.ordis==ordis_ZDir /*6*/ )
    //{/*7*/  //Z=25319
    switch ( CALC.params.orcase )
    {
    case orcGeneral:     // General: phi!=0 && phi!=90 && theta!=0 && theta!=90
    {
        //Z=25320
        //const double limql = sqrt(sqr((qx*CALC.cosphi-qy*CALC.sinphi)*CALC.costheta*CALC.sintheta)
        //                          +sqr(qx+CALC.sinphi*(-qx*CALC.sinphi+qy*CALC.cosphi)*CALC.sintheta*CALC.sintheta)
        //                          +sqr(qy-CALC.cosphi*(-qx*CALC.sinphi+qy*CALC.cosphi)*CALC.sintheta*CALC.sintheta));  //Z=25321
        //pq = CALC.formpq(CALC.params.sigmal,limql,qx,qy,qx*CALC.cosphi*CALC.sintheta,
        //                 qy*CALC.sinphi*CALC.sintheta,q,CALC.ordis);  //Z=25323
        pq = CALC.formpq_partCylinder(qx*CALC.cosphi*CALC.sintheta,
                                      qy*CALC.sinphi*CALC.sintheta,q);  //Z=25323
        /* fq = pq;  //Z=25324 */
        //if ( CALC.lattice )  //Z=25325
        //fq = CALC.formfq( sqrt(CALC.cosphi*qx*qx+CALC.sinphi*qy*qy+eps9), qx, qy,
        //                 qx*CALC.cosphi*CALC.sintheta, qy*CALC.sinphi*CALC.sintheta,q, CALC.ordis );  //Z=25327
        fq = CALC.formfq_partCylinder( sqrt(CALC.cosphi*qx*qx+CALC.sinphi*qy*qy+eps9),
                                      qx*CALC.cosphi*CALC.sintheta, qy*CALC.sinphi*CALC.sintheta,q );  //Z=25327
        break;   //Z=25328
    }
    case orcXaxis:     // X-Axis phi==0 && theta==90
        //Z=25329
        //pq = CALC.formpq(CALC.params.sigmal, fabs(qx), qx, qy, qx, 0, q, CALC.ordis);   //Z=25331
        pq = CALC.formpq_partCylinder(qx, 0, q);   //Z=25331
        /* fq = pq;  //Z=25332 */
        //if ( CALC.lattice )
        //fq = CALC.formfq( fabs(qx), qx, qy, qx, 0, q, CALC.ordis );   //Z=25335
        fq = CALC.formfq_partCylinder( fabs(qx), qx, 0, q );   //Z=25335
        break;   //Z=25336
    case orcYaxis:     // Y-Axis phi==90 && theta==90
        /*Z=24733*/
        //pq = CALC.formpq(CALC.params.sigmal, fabs(qy), qx, qy, 0, qy, q, CALC.ordis);   //Z=25339
        pq = CALC.formpq_partCylinder(0, qy, q);   //Z=25339
        /* fq = pq;  //Z=25340 */
        //if ( CALC.lattice )
        //fq = CALC.formfq( fabs(qy), qx, qy, 0, qy, q, CALC.ordis );   //Z=25343
        fq = CALC.formfq_partCylinder( fabs(qy), 0, qy, q );   //Z=25343
        break;   //Z=25344
    case orcZaxis:     // Z-Axis (phi==0 || phi==90) && theta==0
        /*Z=24741*/
        //pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx, qy, q, CALC.ordis);   //Z=25347
        pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25347
        /* fq = pq;  //Z=25348 */
        //if ( CALC.lattice )
        //fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );   //Z=25351
        fq = CALC.formfq_partCylinder( q, qx, qy, q );   //Z=25351
        break;   //Z=25352
    } // switch orcase
    //} // if ( CALC.ordis==6 ) //Z=25353

    //}  /*  of part=1 (Cylinder)  */   //Z=25395

    double width = 1;
    if ( CALC.RadioButtonDebyeScherrer )
        width = 4.0/CALC.params.domainsize;  //Z=25722
    if ( CALC.RadioButtonPara )
        width = (4.0/CALC.params.domainsize)+sqr(CALC.reldis)*CALC.dist*sqr(q);  //Z=25723

    //if ( ! CALC.tpvRandomOldValues.empty() )   // doIntCalc... - nicht auf der GPU zulässig
    //{
    //if ( fabs(qx) < 0.1 && fabs(qy) < 0.1 && fabs(qz) < 0.1 )
    //    std::cerr << "TPV CALC " << CALC.params.width_zuf << std::endl << std::flush;
    width *= CALC.params.width_zuf;
    //}


    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

    double widthiso = 1.0 / CALC.params.uca;  //ZN=20321

    radintensity = 0.0;     // Immer nur für ein Pixel berechnet, somit kein Array nötig
    intensity = 0.0;

    /* ** lattice hkl-factor calcuation  */  //Z=25744
    //if ( CALC.lattice )
    //{/*6*/  //Z=25745

    /*  PC-version  */  //Z=25748
    //sphno = CALC.latpar[1];  //Z=25749
    //cubevol = CALC.latpar[2];  //Z=25750
    //global: dwfactor = CALC.latpar[3];  //Z=25751
    const double sphno_cubevol = CALC.latpar[1] * CALC.latpar[2]; // werden nur so verwendet

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
        const double qhkl  = CALC.latpar3(ii,5);
        //const int h = CALC.latpar2(ii,1);  //Z=25780  //ZN=25988
        //const int k = CALC.latpar2(ii,2);
        //const int l = CALC.latpar2(ii,3);
        //const int mhkl = CALC.latpar2(ii,4);
        //const int fhkl = CALC.latpar2(ii,5);
        const int mhkl_fhkl = CALC.latpar2(ii,4) * CALC.latpar2(ii,5); // werden nur so genutzt
        if ( qhkl > 0 && mhkl_fhkl != 0 )
        {

            //qhkl0 = CALC.latpar3(ii,1);  //Z=25820  //ZN=26028
            double qxhkl = CALC.latpar3(ii,2);
            double qyhkl = CALC.latpar3(ii,3);
            double qzhkl = CALC.latpar3(ii,4);
            const double qxhklt = CALC.latpar3(ii,7);
            const double qyhklt = CALC.latpar3(ii,8);
            const double qzhklt = CALC.latpar3(ii,9);
            const double g3     = CALC.latpar3(ii,10);
            /*//20240301 raus const:*/ double qhklt  = CALC.latpar3(ii,14);  //ZN=26037
            /*//20240301 raus const:*/ double g3t    = CALC.latpar3(ii,15);  //ZN=26038

            double sq=1.0;
            double psiord=1.0;

            //switch ( CALC.shp )
            //{
            //case cbpeakGaussian /*2*/:                                                          //20210812-E
            //{
            const double peaknorm1 = CALC.latpar3(ii,11);  //Z=26081
            const double yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
            const double x2 = (q-qhkl)*(q-qhkl)/(width*width);
            sq = exp(-4*x2/M_PI)/(M_PI*width/2.0);
            const double x2phi = 4*q*q/(M_PI*sqr(CALC.phiwidth));             /*** a-factor ***/
            psiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
#ifndef __CUDACC__
            if ( isinf(psiord) ) qDebug() << "psiord1" << g3 << x2phi << yphi << peaknorm1;
#endif
            if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
            {
                /*//20240301 raus:*/ const double peaknorm1t = CALC.latpar3(ii,16);   //ZN=26305
                const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhklt));
                const double phiord = g3t*exp(-x2phi*yphi*yphi)/peaknorm1t;
                psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;   //ZN=26313
#ifndef __CUDACC__
                if ( isinf(psiord) )
                    qDebug() << "psiord2" << CALC.params.ceff << phiord << "=" << g3t << x2phi << yphi << peaknorm1t;
#endif
            }
            // TODO Hier tauchen Probleme bim BCC / BCT auf:
            // Wenn TwRatio(ceff)=0 und CheckBoxTwinned=True dann wird falsch gerechnet. Bei allen anderen
            // Kombinationen stimmen die Ergebnisse.
            // 230804: CheckBoxTwinned soll lt. Hr. Förster langfristig rausfliegen. Die Parts mit "if twin then" sind
            //          schon angepasst, die anderen müssen noch überarbeitet werden.
            //break;
            //}

            //default:
            //    return 0;
            //} // switch shp

            //if ( qhkl > 0 ) wird schon oben abgefragt
            //{
            const double savintens = intensity;
            intensity += sq*mhkl_fhkl*psiord;  //Z=26176
            if ( isinf(intensity) )
            {
                intensity = savintens;
#ifndef __CUDACC__
                qDebug() << "INTENS=INF" << sq << mhkl_fhkl << psiord << "ii" << ii << CALC.peakmax2;
#endif
                break;
            }
            //}
        }

    }/*2*/  /* of peak-loop */  //Z=26182

    szqiso = (1+(1)*(8*M_PI*M_PI*M_PI*radintensity/(4*M_PI*sphno_cubevol)-1)*exp(-CALC.params.ucb/*dwfactoriso*/*q*q));   /*Z0311=24804*/
    szq = (1+(fq/pq)*(2*M_PI*intensity/(sphno_cubevol)-1)*exp(-CALC.dwfactor*q*q));   /*Z0311=24805*/
    //}
    //else // (if lattice)
    //{  //Z=26200
    //    szq = 1.0;
    //    szqiso = 1.0;   // zur Sicherheit
    //}

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    const double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    //Z=26207: szq*pq + CALC.iso*szqiso*pqiso
    //Z=30277: xyintensity^[ihex+zmax][i+zmax] = base+izero*xyintensity^[ihex+zmax][i+zmax]+ifluc/(1.0+q*q*rfluc*rfluc);  //Z=30277
    // retval ist der Pixelwert bei [ihex+zmax][i+zmax]

#ifndef __CUDACC__
    if ( retval < -1e8 || isnan(retval) || isinf(retval) )
        qDebug() << "_pCyl_pGau_lBCT_oZDi.h" << "szq"<<szq << "pq"<<pq <<"fq"<<fq << "szqiso"<<szqiso << "pqiso"<<pqiso << "q"<<q
                 << "intens"<<intensity << "radint"<<radintensity << CALC.peakmax1 << "orcase"<<CALC.params.orcase
                 << "cs"<<CALC.params.cs
                 << "="<<retval;
        // szq inf pq 2.16215e-05 fq 2.08036e-07 szqiso 1 pqiso 2.16215e-05 q 2.74004 intens inf radint 2.45523 107 = inf

#else
    //if ( retval < -1e6 )
    //    printf( "szq=%lf pq=%lf fq=%lf szqiso=%lf pqiso=%lf q=%lf intens=%lf radint=%lf erg=%lf\n", szq, pq, fq, szqiso, pqiso, q, intensity, radintensity, retval );
#endif

    return retval;
} /* calc_partCylinder_peakGauss_lattBCT_ordisZDir() */
