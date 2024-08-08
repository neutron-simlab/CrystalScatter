
#ifdef __CUDACC__
// Default kernel if no specific kernel is defined
__global__ void kernel_partCylinder_lattNone_ordisZDir( SasCalc_GENERIC_calculation CALC, bool dofit )
{
    KERNEL_PREWORK
    const double pixval = CALC.calc_partCylinder_lattNone_ordisZDir( CALC, qx, qy, qz );
    KERNEL_POSTWORK
}
#endif



#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::calc_partCylinder_lattNone_ordisZDir(const SasCalc_GENERIC_calculation& CALC,
                                                                      double qx, double qy, double qz)
{
    double pq=0;
    const double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254

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
        break;   //Z=25328
    }
    case orcXaxis:     // X-Axis phi==0 && theta==90
        //Z=25329
        //pq = CALC.formpq(CALC.params.sigmal, fabs(qx), qx, qy, qx, 0, q, CALC.ordis);   //Z=25331
        pq = CALC.formpq_partCylinder(qx, 0, q);   //Z=25331
        break;   //Z=25336
    case orcYaxis:     // Y-Axis phi==90 && theta==90
        /*Z=24733*/
        //pq = CALC.formpq(CALC.params.sigmal, fabs(qy), qx, qy, 0, qy, q, CALC.ordis);   //Z=25339
        pq = CALC.formpq_partCylinder(0, qy, q);   //Z=25339
        break;   //Z=25344
    case orcZaxis:     // Z-Axis (phi==0 || phi==90) && theta==0
        /*Z=24741*/
        //pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx, qy, q, CALC.ordis);   //Z=25347
        pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25347
        break;   //Z=25352
    } // switch orcase
    //} // if ( CALC.ordis==6 ) //Z=25353

    //}  /*  of part=1 (Cylinder)  */   //Z=25395

    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

    //const double szq = 1.0;
    //const double szqiso = 1.0;
    //const double pqiso=1;

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    //const double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    const double retval = CALC.base + CALC.izero*(pq + CALC.iso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    // retval ist der Pixelwert bei [ihex+zmax][i+zmax]

#ifndef __CUDACC__
    if ( retval < -1e8 || isnan(retval) || isinf(retval) )
        qDebug() << "pq"<<pq << "q"<<q
                 << "="<<retval;

//#else
    //if ( retval < -1e6 )
    //    printf( "pq=%lf fq=%lf q=%lf intens=%lf radint=%lf erg=%lf\n", pq, fq, q, intensity, radintensity, retval );
#endif

    return retval;
} /* calc_partCylinder_lattNone_ordisZDir() */
