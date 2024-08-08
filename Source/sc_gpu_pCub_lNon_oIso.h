
#ifdef __CUDACC__
// Default kernel if no specific kernel is defined
__global__ void kernel_partCube_lattNone_ordisIsotropic( SasCalc_GENERIC_calculation CALC, bool dofit )
{
    KERNEL_PREWORK
    const double pixval = CALC.calc_partCube_lattNone_ordisIsotropic( CALC, qx, qy, qz );
    KERNEL_POSTWORK
}
#endif


#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::calc_partCube_lattNone_ordisIsotropic(const SasCalc_GENERIC_calculation& CALC,
                                                                   double qx, double qy, double qz)
{
    double pq; //, fq;

    const double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254

    //const double pqiso=1;
    //const double szqiso=1;
    //const double szq=1;

    //if ( CALC.partcube )
    //{/*6*/  //Z=25506

    //switch ( ComboBoxInterior )
    //{
    //case cbintHomogeneous/*0*/:
    //    params.cs = 0;        /*  homogeneous  */  //Z=20731
    //    homogeneous = true;  //Z=20738

    /*  isotropic cases  */  //Z=25511
    if ( /*(CALC.ordis==7) &&*/ CALC.homogeneous )
    {/*7*/  //Z=25512
        //pq = CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25514
        pq = CALC.formpq_partCube(qx,qy,q,CALC.params.ordis);  //Z=25514
    }/*7*/  //Z=25519
    else
    {
        pq = 1.0;
    }
    //fq = pq;  //Z=25515

    //}/*6*/  //Z=25532 cbpartCube

    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

    // (if !lattice)
    //{  //Z=26200
    //szq = 1.0;
    //szqiso = 1.0;
    //}

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    //const double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    const double retval = CALC.base + CALC.izero*(pq + CALC.iso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    //Z=26207: szq*pq + CALC.iso*szqiso*pqiso
    //Z=30277: xyintensity^[ihex+zmax][i+zmax] = base+izero*xyintensity^[ihex+zmax][i+zmax]+ifluc/(1.0+q*q*rfluc*rfluc);  //Z=30277
    // retval ist der Pixelwert bei [ihex+zmax][i+zmax]

#ifndef __CUDACC__
    if ( retval < -1e8 || isnan(retval) || isinf(retval) )
        qDebug() << "_pCub_lNon_oIso.h" << "pq"<<pq << "q"<<q
                 << "="<<retval;
    // szq inf pq 2.16215e-05 fq 2.08036e-07 szqiso 1 pqiso 2.16215e-05 q 2.74004 intens inf radint 2.45523 107 = inf

#else
    //if ( retval < -1e6 )
    //    printf( "szq=%lf pq=%lf fq=%lf szqiso=%lf pqiso=%lf q=%lf intens=%lf radint=%lf erg=%lf\n", szq, pq, fq, szqiso, pqiso, q, intensity, radintensity, retval );
#endif

    return retval;
} /* calc_partCube_lattNone_ordisIsotropic() */
