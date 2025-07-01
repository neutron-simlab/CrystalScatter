
#ifdef __CUDACC__
// Default kernel if no specific kernel is defined
__global__ void kernel_GENERIC( SasCalc_GENERIC_calculation CALC, bool dofit )
{
    KERNEL_PREWORK
    const double pixval = CALC.calc_GENERIC( CALC, qx, qy, qz );
    KERNEL_POSTWORK
}
#endif


/**
 * @brief SasCalc_GENERIC_calculation::calc_GENERIC
 * @param CALC - reference to class with parameters and subfunctions
 * @param qx   - Coordinates in the q dimension
 * @param qy   - "
 * @param qz   - "
 * @return     - calculated value for this coordinate
 * Calculation from Pascalprogram (20210818-crystal3d1.pas + updates) - only part Generic
 * This function is called from the CPU-Threads and the GPU-Threads to calculate the pixelvalue if no
 *  specialized function is defined (see top of file).
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline double SasCalc_GENERIC_calculation::calc_GENERIC(const SasCalc_GENERIC_calculation& CALC,
                                              double qx, double qy, double qz)
{
    double pq, fq, intensity, radintensity;

    const double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;  //Z=25254

    //Z=25257 if gisaxs ... erstmal weglassen und das auch bei allen weiteren Abfragen!!!

    const double delta = 2.0 * CALC.params.radius;
    double pqiso=1;
    double szqiso=1, szq=1;

    /* ************* */  //Z=25269  //ZN=25476
    /* ** spheres ** */  //Z=25270
    /* ************* */  //Z=25271
    if ( CALC.ComboBoxParticle == cbpartSphere /*CALC.partsphere*/ )
    {/*6*/  //Z=25272

        //limql = 1;
        //pq = CALC.formpq( CALC.params.sigmal, limql, qx, qy, qx, qy, q, CALC.ordis );  //Z=25276
        pq = CALC.formpq_partSphere( q );  //Z=25276

        /* fq = pq;  //Z=25277 */
        if ( CALC.lattice )  //Z=25283
            fq=CALC.formfq_partSphere( q );  //Z=25285
            //fq=CALC.formfq( limqlf, qx, qy, qx, qy, q, CALC.ordis );  //Z=25285
        pqiso = pq;  //Z=25297

    }/*6*/  /*  of part=0  */  //Z=25299


    /* *************** */  //Z=25302
    /* ** cylinders ** */  //Z=25303
    /* *************** */  //Z=25304
    else if ( CALC.ComboBoxParticle == cbpartCylinder /*CALC.partcylinder*/ )
    {/*6*/     /*  cylinders  */  //Z=25305

        /*  isotropic cases  */  //Z=25307
        if ( CALC.params.ordis==7 )
        {/*7*/  //Z=25308
            //pq = CALC.formpq(CALC.params.sigmal, q, qx, qy, qx, qy, q, CALC.ordis);
            pq = CALC.formpq_partCylinder(qx, qy, q);
            if ( CALC.lattice )  //Z=25312
                //fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );
                fq = CALC.formfq_partCylinder( q, qx, qy, q );
        }/*7*/  //Z=25315

        /*  perfect orientation  */               /*  abs((cosphi*qx+sinphi*qy)  //Z=25318 */
        else if ( CALC.params.ordis==ordis_ZDir /*6*/ )
        {/*7*/  //Z=25319
            switch ( CALC.params.orcase )
            {
            case orcGeneral:     // General: phi!=0 && phi!=90 && theta!=0 && theta!=90
            {
                //Z=25320
                //const double limql = sqrt(sqr((qx*CALC.cosphi-qy*CALC.sinphi)*CALC.costheta*CALC.sintheta)
                //             +sqr(qx+CALC.sinphi*(-qx*CALC.sinphi+qy*CALC.cosphi)*CALC.sintheta*CALC.sintheta)
                //             +sqr(qy-CALC.cosphi*(-qx*CALC.sinphi+qy*CALC.cosphi)*CALC.sintheta*CALC.sintheta));  //Z=25321
                //pq = CALC.formpq(CALC.params.sigmal,limql,qx,qy,qx*CALC.cosphi*CALC.sintheta,
                //                 qy*CALC.sinphi*CALC.sintheta,q,CALC.ordis);  //Z=25323
                pq = CALC.formpq_partCylinder( qx*CALC.cosphi*CALC.sintheta, qy*CALC.sinphi*CALC.sintheta, q );
                if ( CALC.lattice )  //Z=25325
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
                if ( CALC.lattice )
                    //fq = CALC.formfq( fabs(qx), qx, qy, qx, 0, q, CALC.ordis );   //Z=25335
                    fq = CALC.formfq_partCylinder( fabs(qx), qx, 0, q );   //Z=25335
                break;   //Z=25336
            case orcYaxis:     // Y-Axis phi==90 && theta==90
                /*Z=24733*/
                //pq = CALC.formpq(CALC.params.sigmal, fabs(qy), qx, qy, 0, qy, q, CALC.ordis);   //Z=25339
                pq = CALC.formpq_partCylinder(0, qy, q);   //Z=25339
                if ( CALC.lattice )
                    //fq = CALC.formfq( fabs(qy), qx, qy, 0, qy, q, CALC.ordis );   //Z=25343
                    fq = CALC.formfq_partCylinder( fabs(qy), 0, qy, q );   //Z=25343
                break;   //Z=25344
            case orcZaxis:     // Z-Axis (phi==0 || phi==90) && theta==0
                /*Z=24741*/
                //pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx, qy, q, CALC.ordis);   //Z=25347
                pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25347
                if ( CALC.lattice )
                    //fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );   //Z=25351
                    fq = CALC.formfq_partCylinder( q, qx, qy, q );   //Z=25351
                break;   //Z=25352
            } // switch orcase
        } // if ( CALC.ordis==6 ) //Z=25353

        /*  general orientation  */  //Z=25357
        else if ( CALC.params.ordis==ordis_Gaussian /*0*/ )
        {   //Z=25358
            switch ( CALC.params.orcase )
            {
            case orcGeneral:    /*  general orientation  */  //Z=25359
                //pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic,
                //                 qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis);   //Z=25361
                pq = CALC.formpq_partCylinder(qx*CALC.cosphic-qy*CALC.sinphic,
                                 qx*CALC.sinphic+qy*CALC.cosphic, q);   //Z=25361
                if ( CALC.lattice ) /* fq:=pq; */   /*Z=24758*/
                    //fq = CALC.formfq( q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic,
                    //                 qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis );   //Z=25365
                    fq = CALC.formfq_partCylinder( q, qx*CALC.cosphic-qy*CALC.sinphic,
                                     qx*CALC.sinphic+qy*CALC.cosphic, q );   //Z=25365
                break;   //Z=25366
            case orcXaxis:   /*  x-axis  */  //Z=25367
                //pq = CALC.formpq(CALC.params.sigmal, sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx, qy, q, CALC.ordis);   //Z=25369
                pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25369
                if ( CALC.lattice )
                    //fq = CALC.formfq( sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx, qy, q, CALC.ordis );   //Z=25374
                    fq = CALC.formfq_partCylinder( sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, q );   //Z=25374
                szq = pq; // ffq;  //Z=25375  TODO das macht hier keinen Sinn
                break;   //Z=25376
            case orcYaxis:  /*  y-axis  */  //Z=25377
                //pq = CALC.formpq(CALC.params.sigmal, sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx, qy, q, CALC.ordis);   //Z=25379
                pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25379
                /* fq = pq;  //Z=25380 */
                if ( CALC.lattice ) /* fq:=pq; */   /*Z=24776*/
                    //fq = CALC.formfq( sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx, qy, q, CALC.ordis );   //Z=25383
                    fq = CALC.formfq_partCylinder( sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, q );   //Z=25383
                break;   //Z=25384
            case orcZaxis:  /*  z-axis  */  //Z=25385
                //pq = CALC.formpq(CALC.params.sigmal,  q, qx, qy, qx, qy, q, CALC.ordis);   //Z=25387
                pq = CALC.formpq_partCylinder(qx, qy, q);   //Z=25387
                /* fq = pq;  //Z=25388 */
                if ( CALC.lattice )
                    //fq = CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );   //Z=25391
                    fq = CALC.formfq_partCylinder( q, qx, qy, q );   //Z=25391
                break;   /*Z=24787*/
            } // switch orcase

        } // if ( CALC.ordis==0 ) //Z=25393

    }  /*  of part=1 (Cylinder)  */   //Z=25395


    /* *********** */  //Z=25398
    /* ** disks ** */  //Z=25399
    /* *********** */  //Z=25400
    else if ( CALC.ComboBoxParticle == cbpartDisk /*CALC.partdisk*/ )
    {/*6*/     /*  disks  */  //Z=24796

        /*  isotropic cases  */  //Z=24798
        if ( CALC.params.ordis == 7 )
        {
            //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);
            pq=CALC.formpq_partDisk(q,qx,qy,qx,qy,q,CALC.params.ordis);
            if ( CALC.lattice )
                fq = CALC.formfq_partDisk( qx, qy, qx, qy, q );
        }

        /*  perfect orientation  */    //Z=24807
        if ( CALC.params.ordis==6 )
        {
            switch ( CALC.params.orcase )
            {
            case orcGeneral:
            {
                const double limql = sqrt(CALC.sinphi*qx*qx+CALC.cosphi*qy*qy+eps6);
                pq=CALC.formpq_partDisk( limql, qx, qy, qx*CALC.cosphi*CALC.sintheta/q, qy*CALC.sinphi*CALC.sintheta/q, q, CALC.params.ordis);
                if ( CALC.lattice )
                    fq = CALC.formfq_partDisk( qx, qy, qx*CALC.cosphi*CALC.sintheta/q, qy*CALC.sinphi*CALC.sintheta/q, q );
                break;
            }
            case orcXaxis:
                //pq=CALC.formpq(CALC.params.sigmal, fabs(qy), qx,qy,qx/q,0,q,CALC.ordis);   //Z=25425
                pq=CALC.formpq_partDisk( fabs(qy), qx, qy, qx/q, 0, q, CALC.params.ordis );
                if ( CALC.lattice )
                    fq = CALC.formfq_partDisk( qx, qy, qx/q, 0, q );
                break;
            case orcYaxis:
                //pq=CALC.formpq(CALC.params.sigmal, fabs(qx), qx,qy,0,qy/q,q,CALC.ordis);  //Z=25433
                pq=CALC.formpq_partDisk(fabs(qx), qx,qy,0,qy/q,q,CALC.params.ordis);
                if ( CALC.lattice )
                    fq = CALC.formfq_partDisk( qx, qy, 0, qy/q, q );
                break;
            case orcZaxis:
                //pq=CALC.formpq(CALC.params.sigmal, q, qx,qy,qx,qy,q,CALC.ordis);   //Z=25441
                pq=CALC.formpq_partDisk(q, qx,qy,qx,qy,q,CALC.params.ordis);
                if ( CALC.lattice )
                    fq = CALC.formfq_partDisk( qx, qy, qx, qy, q );
                break;
            } // switch orcase
        } /* if ( CALC.ordis==6 ) */  //Z=25447

        /*  isotropic fraction  */  //Z=25449
        if ( CALC.iso>0 )
            //pqiso = CALC.formpq(CALC.params.sigmal, q, qx,qy,qx,qy,q, 7/*ordis*/);  //Z=25451
            pqiso = CALC.formpq_partDisk(q, qx,qy,qx,qy,q, 7/*ordis*/);  //Z=25451
        else
            pqiso = 0.0;  //Z=25452

        /*  general orientation  */  //Z=25455
        if ( CALC.params.ordis==0 )
        {
#ifndef __CUDACC__
            //if ( idxCenter ) qDebug() << "TODO formfq, orcase:"<<CALC.params.orcase
            //             << "ordis:"<<CALC.ordis << "part:"<<CALC.params.part << "cs:"<<CALC.params.cs;
            //Debug: TODO formfq, orcase: 1 ordis: 0 part: 2 cs: 0
#endif
            switch ( CALC.params.orcase )
            {
            case orcGeneral:
                //pq=CALC.formpq(CALC.params.sigmal,
                //                 q,qx,qy,qx*CALC.cosphic/q-qy*CALC.sinphic/q,qx*CALC.sinphic/q+qy*CALC.cosphic/q,q,CALC.ordis);   //Z=25459
                pq=CALC.formpq_partDisk(q,qx,qy,qx*CALC.cosphic/q-qy*CALC.sinphic/q,qx*CALC.sinphic/q+qy*CALC.cosphic/q,q,CALC.params.ordis);   //Z=25459
                if ( CALC.lattice )
                    //fq=CALC.formfq( q, qx, qy, qx*CALC.cosphic/q-qy*CALC.sinphic/q, qx*CALC.sinphic/q+qy*CALC.cosphic/q, q, CALC.ordis );   //Z=25463
                    fq=CALC.formfq_partDisk( /*q,*/ qx, qy, qx*CALC.cosphic/q-qy*CALC.sinphic/q, qx*CALC.sinphic/q+qy*CALC.cosphic/q, q/*, CALC.params.ordis*/ );   //Z=25463
                break;  //Z=25464
            case orcXaxis:  //Z=25465
                //pq=CALC.formpq(CALC.params.sigmal, sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx,qy,qx/q,qy/q,q,CALC.ordis);
                pq=CALC.formpq_partDisk(sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx,qy,qx/q,qy/q,q,CALC.params.ordis);
                if ( CALC.lattice ) //fq:=pq;
                    //fq=CALC.formfq( sqrt((1.0-CALC.order)*qx*qx+qy*qy), qx, qy, qx/q, qy/q, q, CALC.ordis );  //Z=25471
                    fq=CALC.formfq_partDisk( /*sqrt((1.0-CALC.order)*qx*qx+qy*qy),*/ qx, qy, qx/q, qy/q, q/*, CALC.params.ordis*/ );  //Z=25471
                break;  //Z=25472
            case orcYaxis:  //Z=25473
                //pq=CALC.formpq(CALC.params.sigmal,sqrt(qx*qx+(1-CALC.order)*qy*qy),qx,qy,qx/q,qy/q,q,CALC.ordis );  //Z=25475
                pq=CALC.formpq_partDisk(sqrt(qx*qx+(1-CALC.order)*qy*qy),qx,qy,qx/q,qy/q,q,CALC.params.ordis );  //Z=25475
                if ( CALC.lattice )
                    //fq=CALC.formfq( sqrt(qx*qx+(1-CALC.order)*qy*qy), qx, qy, qx/q, qy/q, q, CALC.ordis );  //Z=25479
                    fq=CALC.formfq_partDisk( /*sqrt(qx*qx+(1-CALC.order)*qy*qy),*/ qx, qy, qx/q, qy/q, q/*, CALC.params.ordis*/ );  //Z=25479
                break;  //Z=25480
            case orcZaxis:  //Z=25481
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25483
                pq=CALC.formpq_partDisk(q,qx,qy,qx,qy,q,CALC.params.ordis);  //Z=25483
                if ( CALC.lattice ) //fq:=pq;
                    //fq=CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );  //Z=25487
                    fq=CALC.formfq_partDisk( /*q,*/ qx, qy, qx, qy, q/*, CALC.params.ordis*/ );  //Z=25487
                break;  //Z=25488
            } // switch orcase
        } // if ( CALC.ordis==0 )  //Z=25489

    } /*  of part=2, disk  */  //Z=25490


    else if ( CALC.ComboBoxParticle == cbpartVesicle /*CALC.partvesicle*/ )
    {  //Z=25493
        pq = CALC.polyvesicle(/*length,radius,sigma,sigmal,*/q);
        //if ( CheckBoxf2q.Checked==true )  TODO
        //    fq = CALC.f2dpolyvesicle(/*length,radius,sigma,sigmal,*/q);
        //else
        fq = pq;  //Z=25496
    }  //Z=25497


    // Wird aktuell nicht verwendet
    //else if ( CALC.partliposome )
    //{  //Z=25499
    //    /* pq:=polyliposome(length,radius,sigma,sigmal,shellno,alphash,ceff,reff,a,b,c,domainsize,aziwidth,3,q);  //Z=25500 */
    //    /* pq:=liposome(a,rho,radiusi,p1,shellno,q);  //Z=25501 */
    //    pq = liposome(5,0,10,0.8,3,q);  //Z=25502
    //    fq = pq;  //Z=25503
    //}  //Z=25504


    else if ( CALC.ComboBoxParticle == cbpartCube /*CALC.partcube*/ )
    {/*6*/  //Z=25506
        /* pq:=polycube(radius,sigma,0,q);  //Z=25507 */
        /* fq:=polycube(radius,sigma,1,q);  //Z=25508 */

        /*  isotropic cases  */  //Z=25511
        if ( (CALC.params.ordis==7) && CALC.homogeneous )
        {/*7*/  //Z=25512
            //pq = CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25514
            pq = CALC.formpq_partCube(qx,qy,q,CALC.params.ordis);  //Z=25514
            fq = pq;  //Z=25515
            /* if lattice then  //Z=25516 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25517 */
            /*      part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25518 */
        }/*7*/  //Z=25519

        /*  perfect orientation  */  //Z=25521
        if ( (CALC.params.ordis==6) && CALC.homogeneous )
        {/*7*/  //Z=25522
            //pq = CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25524
            pq = CALC.formpq_partCube(qx,qy,q,CALC.params.ordis);  //Z=25524
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25525 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25526 */
            fq = pq;  //Z=25527
        }/*7*/  //Z=25528
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^);  //Z=25529 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^);  //Z=25530 */

    }/*6*/  //Z=25532 cbpartCube


    else if ( CALC.ComboBoxParticle == cbpartEllipsoide /*CALC.partellipsoid*/ )
    {   //Z=25534

        /*  isotropic cases  */  //Z=25535
        if ( CALC.params.ordis==7 && CALC.homogeneous )
        {   /*Z221118=25536*/
            //pq=CALC.formpq(CALC.epsilon,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25538 TODO: ist das epsilon hier richtig?
            pq=CALC.formpq_partEllips(CALC.epsilon,qx,qy,qx,qy,q);  //Z=25538 TODO: ist das epsilon hier richtig?
            fq = pq;  //Z=25539
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25540 */
            /*      part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25541 */
        }  //Z=25542

        /* perfect orientation */
        if ( CALC.params.ordis==6 && CALC.homogeneous )
        {   //Z=25545
            switch ( CALC.params.orcase )
            {
                //double formpq_partEllips( double sigmal, double qx, double qy, double qxs, double qys, double q ) const;

            case orcGeneral:
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25548
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qx*CALC.cosphi,qy*CALC.sinphi,q);  //Z=25548
                break;
            case orcXaxis:
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,0,q,CALC.ordis);  //Z=25551
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qx,qy,q);  //Z=25551
                break;
            case orcYaxis:
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qy,qx,q,CALC.ordis);  //Z=25554
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qy,qx,q);  //Z=25554
                break;
            case orcZaxis:
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25557
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qx,qy,q);  //Z=25557
                break;
            }
            fq = pq;  //Z=25559
        }

        /*  general orientation  */  //Z=25562
        if ( CALC.params.ordis==0 )
        {   //Z=25563
            switch ( CALC.params.orcase )
            {
            case orcGeneral:   /*  general orientation  */  //Z=25564
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx*CALC.cosphic-qy*CALC.sinphic,qx*CALC.sinphic+qy*CALC.cosphic,q,CALC.ordis);  //Z=25566
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qx*CALC.cosphic-qy*CALC.sinphic,qx*CALC.sinphic+qy*CALC.cosphic,q);  //Z=25566
                /* fq = pq;  //Z=25567 */
                //if ( CALC.lattice )
                //    fq=pq; // undef CALC.formfq( q, qx, qy, qx*CALC.cosphic-qy*CALC.sinphic, qx*CALC.sinphic+qy*CALC.cosphic, q, CALC.ordis );  //Z=25570
                break;  //Z=25571
            case orcXaxis:  /*  x-axis  */  //Z=25572
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25574
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qx,qy,q);  //Z=25574
                fq = pq;  //Z=25575
                /* ffq:=pq;  //Z=25576 */
                /* if lattice then //fq:=pq;  //Z=25577 */
                /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1f,limq2f,limq3f,limq4f,limq5f,limq6f,qx,qy,qx,qy,q,norm,  //Z=25578 */
                /*    part,cs,ordis,orcase,myarray,carr1p^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25579 */
                /* szq:=ffq;  //Z=25580 */
                break;  //Z=25581
            case orcYaxis:  /*  y-axis  */  //Z=25582
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25584
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qx,qy,q);  //Z=25584
                /* fq = pq;  //Z=25585 */
                //if ( CALC.lattice )
                //    fq=pq; // undef CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );  //Z=25588
                break;  //Z=25589
            case orcZaxis:  /*  z-axis  */  //Z=25590
                //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25592
                pq=CALC.formpq_partEllips(CALC.params.sigmal,qx,qy,qx,qy,q);  //Z=25592
                /* fq = pq;  //Z=25593 */
                //if ( CALC.lattice )
                //    fq=pq; // undef CALC.formfq( q, qx, qy, qx, qy, q, CALC.ordis );  //Z=25596
                break;  //Z=25597
            }

            fq = pq;  //Z=25600  Wenn if ( CALC.lattice ) greift, dann kommt bei formfq() Null zurück (Ellipsoid undefiniert)

        } // if ordis==0  //Z=25598
    } // partellipsoid  //Z=25601


    else if ( CALC.ComboBoxParticle == cbpartTriaxEllips /*CALC.parttriellipsoid*/ )
    {   //Z=25603

        /*  isotropic cases  */  //Z=25604
        if ( (CALC.params.ordis==7) && CALC.homogeneous )
        {   //Z=25605
            //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25607
            pq=CALC.formpq_partTriaxEllips(qx,qy,q);  //Z=25607
            //fq = pq;  //Z=25608
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25609 */
            /*      part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25610 */
        }  //Z=25611

        /*  perfect orientation  */  //Z=25613
        if ( (CALC.params.ordis==6) && CALC.homogeneous )
        {   //Z=25614
            //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);
            pq=CALC.formpq_partTriaxEllips(qx,qy,q);
            //fq = pq;  //Z=25617
            /* if lattice then  //Z=25618 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25619 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25620 */
        }  //Z=25621
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^,carr11pm^,carr22pm^);  //Z=25622 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^,carr11pm^,carr22pm^);  //Z=25623 */

        fq = pq;  //Z=25625
    }  //Z=25626


    else if ( CALC.ComboBoxParticle == cbpartSuperEllips /*CALC.partbarrel*/ )
    {   //Z=25629

        /*  isotropic cases  */  //Z=25630
        if ( (CALC.params.ordis==7) && CALC.homogeneous )
        {   //Z=25631
            //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25633
            pq=CALC.formpq_partSuperEllips(qx,qy,q,CALC.params.ordis);  //Z=25633
            //fq = pq;  //Z=25634
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25635 */
            /*      part,cs,ordis,orcase,myarray,carr1f^ ,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25636 */
        }  //Z=25637

        /*  perfect orientation  */  //Z=25639
        if ( (CALC.params.ordis==6) && CALC.homogeneous )
        {   //Z=25640
            //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25642
            pq=CALC.formpq_partSuperEllips(qx,qy,q,CALC.params.ordis);  //Z=25642
            //fq = pq;  //Z=25643
            /* if lattice then  //Z=25644 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25645 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25646 */
        }  //Z=25647
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^,carr11pm^,carr22pm^);  //Z=25648 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^,carr11pm^,carr22pm^);  //Z=25649 */

        fq = pq;  //Z=25651
    }  //Z=25652


    else if ( CALC.ComboBoxParticle == cbpartSuperball /*CALC.partball*/ )
    {   //Z=25654

        /*  isotropic cases  */  //Z=25655
        if ( (CALC.params.ordis==7) && CALC.homogeneous )
        {   //Z=25656
            //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25658
            pq=CALC.formpq_partSuperball(qx,qy,q,CALC.params.ordis);  //Z=25658
            //fq = pq;  //Z=25659
            /* if lattice then fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25660 */
            /*      part,cs,ordis,orcase,myarray,carr1f^ ,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25661 */
        }  //Z=25662

        /*  perfect orientation  */  //Z=25664
        if ( (CALC.params.ordis==6) && CALC.homogeneous )
        {   //Z=25665
            //pq=CALC.formpq(CALC.params.sigmal,q,qx,qy,qx,qy,q,CALC.ordis);  //Z=25667
            pq=CALC.formpq_partSuperball(qx,qy,q,CALC.params.ordis);  //Z=25667
            //fq = pq;  //Z=25668
            /* if lattice then  //Z=25669 */
            /* fq:=formfq(length,radius,sigmal,sigma,p1,rho,alphash,theta,phi,q,limq1,limq2,limq3,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,qx,qy,q,norm,  //Z=25670 */
            /*          part,cs,ordis,orcase,myarray,carr1f^,carr2f^,carr3f^,carr4f^,carr5f^,carr6f^,carr7f^,carr11pm^,carr22pm^);  //Z=25671 */
        }  //Z=25672
        /* pq:=formpq(length,radius,1,1,zz,q,limq1,limq3,1,1,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr3p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr7p^,carr2p^,carr11pm^,carr22pm^);  //Z=25673 */
        /* if ((ordis=7) and coreshell) then pq:=formpq(length,radiusi,p1,rho,alphash,theta,phi,zz,q,limq1,limq4,limq5,limq6,limq7,limq8,limq9,qx,qy,1,q,norm,part,cs,ordis,orcase,carr1p^,carr4p^,carr5p^,carr6p^,carr7p^,carr8p^,carr9p^,carr2p^,carr11pm^,carr22pm^);  //Z=25674 */

        fq = pq;  //Z=25676
    }  //Z=25677


#ifdef nochnichteingebaut
    /*  excluded volume chain  */  //Z=25679
    /*  descloi4(rg,nu,gammas,q)  */  //Z=25680
    else if ( CALC.ComboBoxParticle == cbpartChain /*CALC.partchain*/ )
    {  //Z=25681
        pq = descloi4(radius,sigma,length,q);  //Z=25682
        fq = pq;  //Z=25683
    }  //Z=25684

    else if ( CALC.ComboBoxParticle == cbpartkpchain /*CALC.partkpchain*/ )
    {  //Z=25686
        pq = kpchain(length,sigmal,q)*psphered(radius,sigma,2,q);  ;      /*  sigmal = persistence length  */  //Z=25687
        fq = pq;  //Z=25688
    }  //Z=25689
#endif


    double width = 1;
    if ( CALC.RadioButtonDebyeScherrer )
        width = 4.0/CALC.params.domainsize;  //Z=25722
    if ( CALC.RadioButtonPara )
        width = (4.0/CALC.params.domainsize)+sqr(CALC.reldis)*CALC.dist*sqr(q);  //Z=25723
    if ( CALC.lattice && CALC.shp==cbpeakAnisotropicGaussian/*8*/ && CALC.params.ordis==7/*isotropic*/ )
        width = CALC.params.sig.length() /*sqrt(sigx*sigx+sigy*sigy+sigz*sigz)*/ /3.0;  //Z=26118

    //if ( ! CALC.tpvRandomOldValues.empty() )   // doIntCalc... - nicht auf der GPU zulässig
    //{
    //if ( fabs(qx) < 0.1 && fabs(qy) < 0.1 && fabs(qz) < 0.1 )
    //    std::cerr << "TPV CALC " << CALC.params.width_zuf << std::endl << std::flush;
    width *= CALC.params.width_zuf;
    //}


    if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat

    const double widthiso = 1.0 / CALC.params.uca;  //ZN=20321

    radintensity = 0.0;     // Immer nur für ein Pixel berechnet, somit kein Array nötig
    intensity = 0.0;

    /* ** lattice hkl-factor calcuation  */  //Z=25744
    if ( CALC.lattice /* bei LType=None(12) immer false */ )
    {/*6*/  //Z=25745

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
            const double qhkl = CALC.latpar3(ii,5);
            const double x2   = sqr(q-qhkl)/sqr(widthiso);  //Z=25766
            const double sq   = exp(-4*x2/M_PI)/(M_PI*widthiso/2.0);  //Z=25767

            const double savintens = radintensity;
            radintensity += sq*mhkl*fhkl;  //Z=25768
            if ( isinf(radintensity) || isnan(radintensity) )
            {
                radintensity = savintens;
//#ifndef __CUDACC__
//                qDebug() << "RADINTENS=INF" << sq << mhkl << fhkl << "ii" << ii << x2 << widthiso << "q" << q << qhkl;
//#endif
                break;
            }

        } /* of peak loop */  //Z=25769

        for ( int ii=1; ii<=CALC.peakmax2; ii++ )
        {   //Z=25779  //ZN=25987
            if ( CALC._endThread ) return 0;  // Falls Anwender abgebrochen hat
            const double qhkl  = CALC.latpar3(ii,5);
            if ( qhkl <= 0 || isnan(qhkl) ) continue; // Abfrage vor intensity+=...
            const int h = CALC.latpar2(ii,1);
            const int k = CALC.latpar2(ii,2);
            const int l = CALC.latpar2(ii,3);
            const int mhkl = CALC.latpar2(ii,4);
            const int fhkl = CALC.latpar2(ii,5);

            //qhkl0 = CALC.latpar3(ii,1);
            double qxhkl = CALC.latpar3(ii,2);  //Z240301=25199
            double qyhkl = CALC.latpar3(ii,3);
            double qzhkl = CALC.latpar3(ii,4);
            const double qxhklt = CALC.latpar3(ii,7);
            const double qyhklt = CALC.latpar3(ii,8);
            const double qzhklt = CALC.latpar3(ii,9);
            const double g3     = CALC.latpar3(ii,10);
            //20240301 raus: const double qhklt  = CALC.latpar3(ii,14);  //ZN=26037
            //20240301 raus: const double g3t    = CALC.latpar3(ii,15);  //ZN=26038

            double sq=1.0;
            double psiord=1.0;

            switch ( CALC.shp )
            {
            case cbpeakAnisotropicGaussian /*8*/:  //Z=25853
                if ( ! CALC.twin /*CALC.CheckBoxTwinned*/ )
                {   /*  perfect orientation  */  //Z=25854
                    if ( CALC.params.ordis==6/*z-dir*/ || fabs(CALC.params.dbeta) < eps9/*dbeta==0*/ )
                    {   //Z=25855
                        D8( qDebug() << "CPU: shp=8(A): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 );
                        const double dqx = qx-qxhkl;
                        const double dqy = qy-qyhkl;
                        const double dqz = qz-qzhkl;
                        const double dqs1 = (dqx*CALC.params.ax1.x()+dqy*CALC.params.ax1.y()+dqz*CALC.params.ax1.z())/(CALC.ax1n_sigx);
                        const double dqs2 = (dqx*CALC.params.ax2.x()+dqy*CALC.params.ax2.y()+dqz*CALC.params.ax2.z())/(CALC.ax2n_sigy);
                        const double dqs3 = (dqx*CALC.params.ax3.x()+dqy*CALC.params.ax3.y()+dqz*CALC.params.ax3.z())/(CALC.ax3n_sigz);
                        const double x2 = dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                        sq = exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8.0);
                    }
                    else if ( CALC.params.ordis==7/*isotropic*/ )
                    {   /* isotropic orientation */   //Z=25890
                        const double locwidth = CALC.params.sig.length() / 3.0;
                        D8( qDebug() << "CPU: shp=8(B): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 << locwidth );
                        const double x2 = (q-qhkl)*(q-qhkl)/(locwidth*locwidth);  //Z=25893
                        sq = g3*exp(-4*x2/M_PI)/(M_PI*locwidth/2.0);  //Z=25894
                    }
                    else if ( CALC.params.ordis==13/*fiber pattern*/ )
                    {   /* fiber pattern */
                        D8( qDebug() << "CPU: shp=8(C): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 );
                        // rotaxphi   = CALC.polPhi
                        // rotaxtheta = CALC.polTheta
                        CALC.qrombchid(CALC.params.length,CALC.ucl1,CALC.params.p1,CALC.ucl2,CALC.params.alpha,/*CALC.ucl3, (war dbeta)*/
                                       delta,CALC.params.polTheta*M_PI/180.0,CALC.params.polPhi*M_PI/180.0,qx,qy,qz,
                                       CALC.ri11,CALC.ri12,CALC.ri13, CALC.ri21,CALC.ri22,CALC.ri23,
                                       CALC.ri31,CALC.ri32,CALC.ri33,  //Z=25928
                                       qxhkl,qyhkl,qzhkl,qhkl,
                                       /*CALC.params.ordis,3,5,*/7,h,k,l, CALC.params.CR->carr1p, sq);  //Z=25930
                        sq = sq*2*M_PI*qhkl/(2.0*M_PI*sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8.0);  //Z=25931
                    }
                    else // if ( (CALC.ordis!=6) && (CALC.ordis!=7) && (CALC.ordis!=13) && (CALC.dbeta!=0) )      //Z=25985
                    {   /* other anisotropic cases, rotated around the qhkl-axis direction */                                                             //20210812-D
                        D8( qDebug() << "CPU: shp=8(D): ihex,i" << ihex << i << "ii" << ii << "/" << CALC.peakmax2 );
                        const double phi = atan2(qyhkl,(qxhkl+eps6)) * 180.0/M_PI;
                        const double theta = atan2(sqrt(qxhkl*qxhkl+qyhkl*qyhkl),(qzhkl+eps6)) * 180.0/M_PI;
                        CALC.qrombdeltac(CALC.params.p1, CALC.params.sigma, CALC.params.alpha, // CALC.params.ax*
                                         theta, phi, qx, qy, qz, qxhkl, qyhkl, qzhkl, qhkl,
                                         /*CALC.params.ordis,3,*/ /*i0=*/5, /*i1=*/6,0,0,0, CALC.params.CR->carr1p, sq );
                        sq = sq*2*M_PI*qhkl/CALC.params.norm;  //Z=25991
                    }
                } // if ! twin

                psiord = 1.0;

                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {   //Z240301=25312
                    double sqtwin=0;
                    /* perfect orientation */
                    if ( CALC.params.ordis==6/*z-dir*/ || fabs(CALC.params.dbeta) < eps9/*dbeta==0*/ )
                    {
                        const double dqx = qx-qxhklt;
                        const double dqy = qy-qyhklt;
                        const double dqz = qz-qzhklt;
                        const double dqs1 = (dqx*CALC.params.ax1.x()+dqy*CALC.params.ax1.y()+dqz*CALC.params.ax1.z())/(CALC.ax1n_sigx);
                        const double dqs2 = (dqx*CALC.params.ax2.x()+dqy*CALC.params.ax2.y()+dqz*CALC.params.ax2.z())/(CALC.ax2n_sigy);
                        const double dqs3 = (dqx*CALC.params.ax3.x()+dqy*CALC.params.ax3.y()+dqz*CALC.params.ax3.z())/(CALC.ax3n_sigz);
                        const double x2 = dqs1*dqs1+dqs2*dqs2+dqs3*dqs3;                        /*** different for twin ***/
                        sqtwin = exp(-4*x2/M_PI)/(sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8);
                        //sqtwin:=sqtwin*2*pi*qhkl;
                    }
                    else if ( CALC.params.ordis==7/*isotropic*/ )
                    {   /* isotropic orientation */
                        const double locwidth = CALC.params.sig.length() / 3.0;
                        const double x2 = (q-qhkl)*(q-qhkl)/sqr(locwidth);
                        sqtwin = g3*exp(-4*x2/M_PI)/(M_PI*locwidth/2.0);
                    }
                    else if ( CALC.params.ordis==13/*fiber pattern*/ )
                    {   /* fiber pattern */
                        CALC.qrombchid(CALC.params.length,CALC.ucl1,CALC.params.p1,CALC.ucl2,CALC.params.alpha,/*CALC.ucl3, (war dbeta)*/
                                       delta,CALC.params.polTheta*M_PI/180.0,CALC.params.polPhi*M_PI/180.0,qx,qy,qz,
                                       CALC.ri11,CALC.ri12,CALC.ri13, CALC.ri21,CALC.ri22,CALC.ri23,
                                       CALC.ri31,CALC.ri32,CALC.ri33,  //Z=25928
                                       qxhkl,qyhkl,qzhkl,qhkl,
                                       /*CALC.params.ordis,3,5,*/7,h,k,l, CALC.params.CR->carr1p, sqtwin);  //Z=25930
                        sqtwin = sqtwin*2*M_PI*qhkl/(2*M_PI*sqrt(M_PI*M_PI*M_PI)*CALC.params.sig.x()*CALC.params.sig.y()*CALC.params.sig.z()/8.0);
                    }
                    else // if ((ordis<>6) and (ordis<>7) and (ordis<>13) and (dbeta<>0)) then begin
                    {    /* other anistropic cases */
                        const double phi = atan2(qyhkl,(qxhkl+eps6)) * 180.0/M_PI;
                        const double theta = atan2(sqrt(qxhkl*qxhkl+qyhkl*qyhkl),(qzhkl+eps6)) * 180.0/M_PI;
                        CALC.qrombdeltac(CALC.params.p1, CALC.params.sigma, CALC.params.alpha, // CALC.params.ax*
                                         theta, phi, qx, qy, qz, qxhkl, qyhkl, qzhkl, qhkl,
                                         /*CALC.params.ordis,3,*/ /*i0=*/5, /*i1=*/6,0,0,0, CALC.params.CR->carr1p, sq );
                        sqtwin = sqtwin*2*M_PI*qhkl/CALC.params.norm;
                    }
                    sq = CALC.params.ceff*sq+(1-CALC.params.ceff)*sqtwin;
                } // if twin
                break;  // shp == cbpeakAnisotropicGaussian  //Z=26062

            case cbpeakLorentzian /*1*/:
            {
                const double peaknorm1 = CALC.latpar3(ii,11);  //Z=26066
                const double yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                const double x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = sqrt(CALC.c1)/(M_PI*width*(1.0+CALC.c1*x2));
                const double x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** b-factor ***/
                psiord = g3/(peaknorm1*(1+x2phi*yphi*yphi));  //Z=26071
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {
                    const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    const double phiord = g3/(peaknorm1*(1+x2phi*yphi*yphi));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26076
                break;
            }

            case cbpeakGaussian /*2*/:                                                          //20210812-E
            {
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
                    const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    const double phiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;   //ZN=26313
#ifndef __CUDACC__
                    if ( isinf(psiord) )
                    {
                        const double peaknorm1t = CALC.latpar3(ii,16);   // nur für Debug
                        qDebug() << "psiord2" << CALC.params.ceff << phiord << "=" << g3 << x2phi << yphi << peaknorm1t;
                    }
#endif
                }
                // TODO Hier tauchen Probleme bim BCC / BCT auf:
                // Wenn TwRatio(ceff)=0 und CheckBoxTwinned=True dann wird falsch gerechnet. Bei allen anderen
                // Kombinationen stimmen die Ergebnisse.
                // 230804: CheckBoxTwinned soll lt. Hr. Förster langfristig rausfliegen. Die Parts mit "if twin then" sind
                //          schon angepasst, die anderen müssen noch überarbeitet werden.
                break;
            }

            case cbpeakMod1Lorentzian /*Lorentzian1*/:
            {
                const double peaknorm1 = CALC.latpar3(ii,11);  //Z=26100
                const double yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                const double x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = 2*sqrt(CALC.c2)/(M_PI*width*(1+CALC.c1*x2)*(1+CALC.c1*x2));
                const double x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*(1+x2phi*yphi*yphi)*(1+x2phi*yphi*yphi));
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {  //Z=26106
                    const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    const double phiord = g3/(peaknorm1*(1+x2phi*yphi*yphi)*(1+x2phi*yphi*yphi));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26110
                break;
            }

            case cbpeakMod2Lorentzian /*Lorentzian2*/:
            {
                const double peaknorm1 = CALC.latpar3(ii,11);  //Z=26115
                const double yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                const double x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = sqrt(CALC.c3)/(2*width*exp(3*log(1+CALC.c1*x2)/2.0));
                const double x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*exp(3*log(1+x2phi*yphi*yphi)/2.0));
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {  //Z=26121
                    const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    const double phiord = g3/(peaknorm1*exp(3*log(1+x2phi*yphi*yphi)/2.0));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26125
                break;
            }

            case cbpeakPseudoVoigt /*Voigt*/:
            {
                const double peaknorm1 = CALC.latpar3(ii,11);  //Z=26131
                const double peaknorm2 = CALC.latpar3(ii,12);
                const double yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                const double x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = CALC.eta*sqrt(CALC.c1)/(M_PI*width*(1+CALC.c1*x2))+(1-CALC.eta)*sqrt(CALC.c0)*exp(-CALC.c0*x2)/(sqrt(M_PI)*width);  //Z=26135
                const double x2psi = 4*q*q/(M_PI*sqr(CALC.phiwidth));             /*** a-factor ***/
                //double x2psihkl = 4*qhkl*qhkl/(M_PI*sqr(CALC.phiwidth));  //Z=26137 wird nicht weiter verwendet
                const double x2phi = 4*q*q/sqr(CALC.phiwidth);                /*** b-factor ***/
                psiord = g3*(CALC.eta*(1/(1+x2phi*yphi*yphi))+(1-CALC.eta)*exp(-x2psi*yphi*yphi))/(CALC.eta*peaknorm1+(1-CALC.eta)*peaknorm2);  //Z=26139
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {  //Z=26140
                    const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    const double phiord = g3*(CALC.eta*(1/(1+x2phi*yphi*yphi))+(1-CALC.eta)*exp(-x2psi*yphi*yphi))/(CALC.eta*peaknorm1+(1-CALC.eta)*peaknorm2);
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26144
                break;
            }

            case cbpeakPearsonVII /*Pearson*/:
            {
                const double peaknorm1 = CALC.latpar3(ii,11);  //Z=26149
                const double yphi = myacos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                const double x2 = (q-qhkl)*(q-qhkl)/(width*width);
                sq = CALC.gamma(CALC.beta)*CALC.c4*2/(CALC.gamma(CALC.beta-0.5)*M_PI*width*exp(CALC.beta*log(1+4*CALC.c4*x2)));
                const double x2phi = 4*q*q/sqr(CALC.phiwidth);             /*** c-factor ***/
                psiord = g3/(peaknorm1*exp(CALC.beta*log(1+x2phi*yphi*yphi)));
                if ( CALC.twin /*CALC.CheckBoxTwinned*/ )
                {  //Z=26155
                    const double yphi = myacos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                    const double phiord = g3/(peaknorm1*exp(CALC.beta*log(1+x2phi*yphi*yphi)));
                    psiord = CALC.params.ceff*psiord+(1-CALC.params.ceff)*phiord;
                }  //Z=26159
                break;
            }

#ifdef undef
                if ( BurgerG )
                {/*8*/  //Z240301=25443
                    /* x2phihkl:=4*qhkl*qhkl/(pi*phiwidth*phiwidth);   */
                    /* peaknorm1:=gaussnorm3(x2phihkl);   */
                    peaknorm1 = latpar3p^[ii][11];
                    yphi = acos((qx*qxhkl+qy*qyhkl+qz*qzhkl)/(q*qhkl));
                    x2 = (q-qhkl)*(q-qhkl)/(width*width);
                    sq = burger(width,bnu,x2);
                    x2phi = 4*q*q/(M_PI*phiwidth*phiwidth);             /* ** a-factor ** */
                    psiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
                    if ( CheckBoxTwinned.Checked==true )
                    {/*9*/
                        yphi = acos((qx*qxhklt+qy*qyhklt+qz*qzhklt)/(q*qhkl));
                        phiord = g3*exp(-x2phi*yphi*yphi)/peaknorm1;
                        psiord = ceff*psiord+(1-ceff)*phiord;
                    }/*9*/
                }/*8*/
#endif

            default:
                return 0;
            } // switch shp

            //if ( qhkl > 0 ) wird schon oben abgefragt
            {
                const double savintens = intensity;
                intensity += sq*mhkl*fhkl*psiord;  //Z=26176
                if ( isinf(intensity) || isnan(intensity) )
                {
                    intensity = savintens;
//#ifndef __CUDACC__
//                    qDebug() << "INTENS=INF" << sq << mhkl << fhkl << psiord << "ii" << ii << CALC.peakmax2;
//#endif
                    break;
                }
            }

        }/*2*/  /* of peak-loop */  //Z=26182

        szqiso = (1+(1)*(8*M_PI*M_PI*M_PI*radintensity/(4*M_PI*sphno_cubevol)-1)*exp(-CALC.params.ucb/*dwfactoriso*/*q*q));   /*Z0311=24804*/
        szq = (1+(fq/pq)*(2*M_PI*intensity/(sphno_cubevol)-1)*exp(-CALC.dwfactor*q*q));   /*Z0311=24805*/
    }
    else // (if lattice)
    {  //Z=26200
        szq = 1.0;
        szqiso = 1.0;   // zur Sicherheit
    }

    // Abschlussberechnungen (izero,base) machen nur hier Sinn. Im Pascalprogramm wurde dies nach den kompletten Schleifen gemacht
    const double retval = CALC.base + CALC.izero*(szq*pq + CALC.iso*szqiso*pqiso) + CALC.ifluc/(1+q*q*CALC.rfluc*CALC.rfluc);
    //Z=26207: szq*pq + CALC.iso*szqiso*pqiso
    //Z=30277: xyintensity^[ihex+zmax][i+zmax] = base+izero*xyintensity^[ihex+zmax][i+zmax]+ifluc/(1.0+q*q*rfluc*rfluc);  //Z=30277
    // retval ist der Pixelwert bei [ihex+zmax][i+zmax]

//#ifndef __CUDACC__
//    if ( retval < -1e8 || isnan(retval) || isinf(retval) )
//        qDebug() << "_generic.h" << "szq"<<szq << "pq"<<pq <<"fq"<<fq << "szqiso"<<szqiso << "pqiso"<<pqiso << "q"<<q
//                 << "intens"<<intensity << "radint"<<radintensity << CALC.peakmax1
//                 << "="<<retval;
    // szq inf pq 2.16215e-05 fq 2.08036e-07 szqiso 1 pqiso 2.16215e-05 q 2.74004 intens inf radint 2.45523 107 = inf

//#else
    //if ( retval < -1e6 )
    //    printf( "szq=%lf pq=%lf fq=%lf szqiso=%lf pqiso=%lf q=%lf intens=%lf radint=%lf erg=%lf\n", szq, pq, fq, szqiso, pqiso, q, intensity, radintensity, retval );
//#endif

    return retval;
} /* calc_GENERIC() */
