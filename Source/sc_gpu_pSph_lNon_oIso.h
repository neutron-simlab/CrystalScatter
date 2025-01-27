
#ifdef __CUDACC__
__global__ void kernel_partSphere_lattNone_ordisIsotropic( SasCalc_GENERIC_calculation CALC, bool dofit )
{
    KERNEL_PREWORK
    const double pixval = CALC.calc_partSphere_lattNone_ordisIsotropic( CALC, qx, qy, qz );
    KERNEL_POSTWORK
}
#endif


#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::calc_partSphere_lattNone_ordisIsotropic(const SasCalc_GENERIC_calculation& CALC,
                                                         double qx, double qy, double qz)
{
    const double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254

    /* ************* */  //Z=25269  //ZN=25476
    /* ** spheres ** */  //Z=25270
    /* ************* */  //Z=25271
    //if ( CALC.ComboBoxParticle == cbpartSphere /*CALC.partsphere*/ )
    //{/*6*/  //Z=25272

#define limql 1.0
    double pq = CALC.formpq_partSphere( q );  //Z=25276
//#define pqiso pq

    //}/*6*/  /*  of part=0  */  //Z=25299

    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

//#define szq 1.0
//#define szqiso 1.0

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    //double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    double retval = CALC.base + CALC.izero*(pq + CALC.iso*pq) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);

    //if ( q > 2.0 )
    //    qDebug() << "pSph_lNon_oIso" << qx << qy << qz << "cs=" << CALC.params.cs << "0.4*limq4=" << (0.4*CALC.params.limq4)
    //             << "pq=" << pq << "ret=" << retval;

    return retval;

#undef limql
#undef pqiso
#undef szq
#undef szqiso
} /* calc_partSphere_lattNone_ordisIsotropic() */

