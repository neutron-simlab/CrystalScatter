
#ifdef __CUDACC__
__global__ void kernel_partCylinder_lattNone_ordisGauss( SasCalc_GENERIC_calculation CALC, bool dofit )
{
    KERNEL_PREWORK
    const double pixval = CALC.calc_partCylinder_lattNone_ordisGauss( CALC, qx, qy, qz );
    KERNEL_POSTWORK
}
#endif


#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::calc_partCylinder_lattNone_ordisGauss(const SasCalc_GENERIC_calculation& CALC,
                                                         double qx, double qy, double qz)
{
    double pq, fq;

    const double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254

    fq = 1; // Verhindern einer Compiler-Warnung
    // <fq> wird bei der Berechnung von szq und somit für die Intensität verwendet (=fq/pq).
    // Aber nur dann, wenn lattice gesetzt ist. Und nur dann wird es auch per formfq berechnet.

    /* *************** */  //Z=25302
    /* ** cylinders ** */  //Z=25303
    /* *************** */  //Z=25304
    //if ( CALC.partcylinder )
    //{/*6*/     /*  cylinders  */  //Z=25305

    /*  general orientation  */  //Z=25357
    //else if ( CALC.ordis==ordis_Gaussian /*0*/ )
    //{   //Z=25358

    switch ( CALC.params.orcase )
    {
    case 1:    /*  general orientation  */  //Z=25359
        //pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic,
        //                 qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis);   //Z=25361
        pq = CALC.formpq_partCylinder(qx*CALC.cosphic-qy*CALC.sinphic,
                                      qx*CALC.sinphic+qy*CALC.cosphic, q);   //Z=25361
        /* fq = pq;  //Z=25362 */
        if ( CALC.lattice ) /* fq:=pq; */   /*Z=24758*/
            //fq = CALC.formfq( q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic,
            //                 qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis );   //Z=25365
            fq = CALC.formfq_partCylinder( q, qx*CALC.cosphic-qy*CALC.sinphic,
                                          qx*CALC.sinphic+qy*CALC.cosphic, q );   //Z=25365
        break;   //Z=25366
    case 2:   /*  x-axis  */  //Z=25367
        //pq = CALC.formpq(CALC.params.sigmal, sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx, qy, q, CALC.ordis);   //Z=25369
        pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25369
        /* fq = pq;  //Z=25370 */
        //ffq = pq;  //Z=25371
        if ( CALC.lattice )
            //fq = CALC.formfq( sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx, qy, q, CALC.ordis );   //Z=25374
            fq = CALC.formfq_partCylinder( sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, q );   //Z=25374
        break;   //Z=25376
    case 3:  /*  y-axis  */  //Z=25377
        //pq = CALC.formpq(CALC.params.sigmal, sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx, qy, q, CALC.ordis);   //Z=25379
        pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25379
        /* fq = pq;  //Z=25380 */
        if ( CALC.lattice ) /* fq:=pq; */   /*Z=24776*/
            //fq = CALC.formfq( sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx, qy, q, CALC.ordis );   //Z=25383
            fq = CALC.formfq_partCylinder( sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, q );   //Z=25383
        break;   //Z=25384
    case 4:  /*  z-axis  */  //Z=25385
        //pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx, qy, q, CALC.ordis);   //Z=25387
        pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25387
        /* fq = pq;  //Z=25388 */
        if ( CALC.lattice )
            //fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );   //Z=25391
            fq = CALC.formfq_partCylinder( q, qx, qy, q );   //Z=25391
        break;   /*Z=24787*/
    } // switch orcase

    //} // if ( CALC.ordis==0 ) //Z=25393

    //}  /*  of part=1 (Cylinder)  */   //Z=25395

    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

    //const double szq = 1.0;
    //const double szqiso = 1.0;
    //const double pqiso=1;

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    //double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    double retval = CALC.base + CALC.izero*(pq + CALC.iso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);

#ifndef __CUDACC__
    if ( retval < -1e8 || isnan(retval) || isinf(retval) )
        qDebug() << "_pCyl_lNon_oGau.h" << "pq"<<pq <<"fq"<<fq << "q"<<q << "="<<retval;
//#else
//    if ( retval < -1e6 )
//        printf( "pq=%lf fq=%lf q=%lf erg=%lf\n", pq, fq, q, retval );
#endif

    return retval;
} /* calc_partCylinder_lattNone_ordisGauss() */

