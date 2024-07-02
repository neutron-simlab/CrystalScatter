#ifndef SC_LIB_FORMFQ_partDisk_H
#define SC_LIB_FORMFQ_partDisk_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formfq_partDisk( /*double limql,*/ double qx, double qy, double qxs, double qys, double q/*, int ordis*/ ) const
{  //Z=18395

    double pqr, pql;  //Z=18401

    const double zl = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=18418
    const double zr = (1-sqr(params.sigma))/sqr(params.sigma);  //Z=18419
    const double radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=18420
    const double qz=1.0; // TODO: in qrombdeltac Aufruf verwendet, siehe auch bei formpq()

    CHECKENDTHREAD_VAL;

    /* ****** */  //Z=19222
    /*  disk  */  //Z=19223
    /* ****** */  //Z=19224
    //if ( params.part==2 )
    //{/*2*/  //Z=19225

        /* ** longitudinal part ** */  //Z=19227
        /* ** isotropic ** */  //Z=19228
        if ( params.ordis==7 )
        {/*3*/  //Z=19229
            if ( q<(0.5*params.limq1f) )
            {/*4*/  //Z=19230
                pql = 1.0;  //Z=19231
                double oldpqsum = 0.0;  //Z=19232
                double qqnn = 1.0;  //Z=19233
                for ( int nser=1; nser<=80; nser++ )
                {/*5*/  //Z=19234
                    qqnn = qqnn*q*q;  //Z=19235
                    pql += params.CR->carr1f[nser]*qqnn;  //Z=19236
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=19237
                    if ( delser<0.0001 ) break; /* goto 70; */  //Z=19238
                    oldpqsum = pql;  //Z=19239
                }/*5*/  //Z=19240
                /*70:*/  //Z=19241
                //pql = pqsum;  //Z=19242
            }/*4*/  //Z=19243
            else
            {/*4*/  /*   = P||(q)   */  //Z=19244
                const double arglq = q*params.length/(zl+1);  //Z=19245
                pql = (2/(zl*(zl-1)))*pow(arglq,-2);  //Z=19246
            }/*4*/  //Z=19247
        }/*3*/  /*  of isotropic  */  //Z=19248

        /*  perfect  */  //Z=19250
        if ( params.ordis==6 )
        {/*3*/  //Z=19251
            if ( q<(0.5*params.limq1f) )
            {/*4*/  //Z=19252
                pql = 1.0;  //Z=19253
                double oldpqsum = 0.0;  //Z=19254
                double qqnn = 1.0;  //Z=19255
                if ( params.orcase==1 )
                {/*5*/  //Z=19256
                    const double argq = qxs+qys;  //Z=19257
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19258
                        qqnn = qqnn*q*q;  //Z=19259
                        pql += params.CR->carr1f[nser]*qqnn*pow(1.0-argq*argq,nser);  //Z=19260
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19261
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19262
                        oldpqsum = pql;  //Z=19263
                    }/*6*/  //Z=19264
                }/*5*/  //Z=19265
                else if ( params.orcase==2 )
                {/*5*/  //Z=19266
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19267
                        qqnn = qqnn*q*q;  //Z=19268
                        pql += params.CR->carr1f[nser]*qqnn*pow(1.0-qxs*qxs,nser);  //Z=19269
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19270
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19271
                        oldpqsum = pql;  //Z=19272
                    }/*6*/  //Z=19273
                }/*5*/  //Z=19274
                else if ( params.orcase==3 )
                {/*5*/  //Z=19275
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19276
                        qqnn = qqnn*q*q;  //Z=19277
                        pql += params.CR->carr1f[nser]*qqnn*pow(1.0-qys*qys,nser);  //Z=19278
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19279
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19280
                        oldpqsum = pql;  //Z=19281
                    }/*6*/  //Z=19282
                }/*5*/  //Z=19283
                else if ( params.orcase==4 )
                {/*5*/  //Z=19284
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19285
                        qqnn = qqnn*q*q;  //Z=19286
                        pql += params.CR->carr1f[nser]*qqnn;  //Z=19287
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19288
                        if ( delser<0.0001 ) break; /* goto 76; */  //Z=19289
                        oldpqsum = pql;  //Z=19290
                    }/*6*/  //Z=19291
                }/*5*/  //Z=19292
                /*76:*/  //Z=19293
                //pql = pqsum;  //Z=19294
            }/*4*/  //Z=19295
            else
            {/*4*/      /*  ok  */  //Z=19296
                double arglq;
                if ( params.orcase==1 )
                {/*5*/  //Z=19297
                    const double qnarg = qxs+qys;  //Z=19298
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
                const double pqr1 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3);  //Z=19312
                const double pqr2 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3)*sin((zl-2)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-2)/2.0);  //Z=19313
                const double pqr3 = (1/(zl*(zl-1)*(zl-2)*(zl-3)))*pow(arglq,-4)*cos((zl-3)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-3)/2.0);  //Z=19314
                pql = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=19315
            }/*4*/  //Z=19316
        }/*3*/   /*  of perfect  */  //Z=19317

        /*  orientational distribution  */  //Z=19319
        if ( params.ordis==0 )
        {/*3*/  //Z=19320
            double qxn[121], qyn[121];
            if ( params.orcase==1 )
            {/*4*/  //Z=19321
                if ( q<(1.2*params.limq1f) )  // war 1.2
                {/*5*/  //Z=19322
                    pql = 1.0;  //Z=19323
                    double oldpqsum = 0.0;  //Z=19324
                    double qqnn = 1.0;  //Z=19325
                    qxn[0] = 1.0;  //Z=19326
                    qyn[0] = 1.0;  //Z=19327

                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19329
                        qqnn = qqnn*q*q;  //Z=19330
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19331
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19332

                        double binsum = 0.0;  //Z=19334
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19335
                            double binsum1 = 0.0;  //Z=19336
                            for ( int lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=19337
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=19338 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=19339 */
                                binsum1 += params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=19340
                            }/*8*/  //Z=19341
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=19342 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=19343 */
                            binsum += params.CR->carr11pm[nser][mser]*binsum1;  //Z=19344
                        }/*7*/  //Z=19345
                        pql += params.CR->carr1f[nser]*qqnn*binsum;  //Z=19346
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19347
                        if ( delser<0.0001 ) break; /* goto 77; */  //Z=19348
                        oldpqsum = pql;  //Z=19349
                    }/*6*/  //Z=19350
                    /*77:*/  //Z=19351
                    //pql = pqsum;  //Z=19352
                }/*5*/  //Z=19353
                else
                {/*5*/  //Z=19354
                    /*  disk: length = disk radius  */  //Z=19355
                    /*  always use Bessel function approximation  */  //Z=19356
                    /*  F(q)  */  //Z=19357
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19358 */
                    /*  P(q)  */  //Z=19359
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz, // 9,9,
                                9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19360
                    pql = pql/params.norm;  //Z=19361
                }/*5*/  //Z=19362
            }/*4*/  //Z=19363

            else if ( params.orcase==2 )
            {/*4*/  //Z=19365
                if ( q<(0.9*params.limq1f) )
                {/*5*/  //Z=19366
                    pql = 1.0;  //Z=19367
                    double oldpqsum = 0.0;  //Z=19368
                    double qqnn = 1.0;  //Z=19369
                    qxn[0] = 1.0;  //Z=19370
                    qyn[0] = 1.0;  //Z=19371

                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19373
                        qqnn = qqnn*q*q;  //Z=19374
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19375
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19376

                        double binsum = 0.0;  //Z=19378
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19379
                            double binsum1 = 0.0;  //Z=19380
                            for ( int lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=19381
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=19382 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=19383 */
                                binsum1 += params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=19384
                            }/*8*/  //Z=19385
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=19386 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=19387 */
                            binsum += params.CR->carr11pm[nser][mser]*binsum1;  //Z=19388
                        }/*7*/  //Z=19389
                        pql += params.CR->carr1f[nser]*qqnn*binsum;  //Z=19390
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19391
                        if ( delser<0.0001 ) break; /* goto 78; */  //Z=19392
                        oldpqsum = pql;  //Z=19393
                    }/*6*/  //Z=19394
                    /*78:*/  //Z=19395
                    //pql = pqsum;  //Z=19396
                }/*5*/  //Z=19397
                else
                {/*5*/  //Z=19398
                    /*  disk: length = disk radius  */  //Z=19399
                    /*  always use Bessel function approximation  */  //Z=19400
                    /*  F(q)  */  //Z=19401
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19402 */
                    /*  P(q)  */  //Z=19403
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz, // 9,9,
                                9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19404
                    pql = pql/params.norm;  //Z=19405
                }/*5*/  //Z=19406
            }/*4*/  //Z=19407

            else if ( params.orcase==3 )
            {/*4*/  //Z=19409
                if ( q<(0.9*params.limq1f) )
                {/*5*/  //Z=19410
                    pql = 1.0;  //Z=19411
                    double oldpqsum = 0.0;  //Z=19412
                    double qqnn = 1.0;  //Z=19413
                    qxn[0] = 1.0;  //Z=19414
                    qyn[0] = 1.0;  //Z=19415

                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19417
                        qqnn = qqnn*q*q;  //Z=19418
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19419
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19420

                        double binsum = 0.0;  //Z=19422
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19423
                            double binsum1 = 0.0;  //Z=19424
                            for ( int lser=0; lser<=mser; lser++ )
                            {/*8*/  //Z=19425
                                /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=19426 */
                                /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=19427 */
                                binsum1 += params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=19428
                            }/*8*/  //Z=19429
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=19430 */
                            /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=19431 */
                            binsum += params.CR->carr11pm[nser][mser]*binsum1;  //Z=19432
                        }/*7*/  //Z=19433
                        pql += params.CR->carr1f[nser]*qqnn*binsum;  //Z=19434
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19435
                        if ( delser<0.0001 ) break; /* goto 79; */  //Z=19436
                        oldpqsum = pql;  //Z=19437
                    }/*6*/  //Z=19438
                    /*79:*/  //Z=19439
                    //pql = pqsum;  //Z=19440
                }/*5*/  //Z=19441
                else
                {/*5*/  //Z=19442
                    /*  disk: length = disk radius  */  //Z=19443
                    /*  always use Bessel function approximation  */  //Z=19444
                    /*  F(q)  */  //Z=19445
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19446 */
                    /*  P(q)  */  //Z=19447
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz, // 9,9,
                                9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19448
                    pql = pql/params.norm;  //Z=19449
                }/*5*/  //Z=19450
            }/*4*/  //Z=19451

            else if ( params.orcase==4 )
            {/*4*/  //Z=19453
                if ( q<(0.5*params.limq1f) )
                {/*5*/  //Z=19454
                    pql = 1.0;  //Z=19455
                    double oldpqsum = 0.0;  //Z=19456
                    double qqnn = 1.0;  //Z=19457
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19458
                        qqnn = qqnn*q*q;  //Z=19459
                        pql += params.CR->carr1f[nser]*qqnn;  //Z=19460
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=19461
                        if ( delser<0.0001 ) break; /* goto 80; */  //Z=19462
                        oldpqsum = pql;  //Z=19463
                    }/*6*/  //Z=19464
                    /*80:*/  //Z=19465
                    //pql = pqsum;  //Z=19466
                }/*5*/  //Z=19467
                else
                {/*5*/  //Z=19468
                    /*  disk: length = disk radius  */  //Z=19469
                    /*  always use Bessel function approximation  */  //Z=19470
                    /*  F(q)  */  //Z=19471
                    /* qrombdeltac(length,radius,p1,sigmal,dbeta,theta,phi,qx,qy,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,2,6,orcase+7,0,1,0,carr2p,pql);  //Z=19472 */
                    /*  P(q)  */  //Z=19473
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz, // 9,9,
                                9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2f,pql);  //Z=19474
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
                pqr = 1.0;  //Z=19485
                double oldpqsum = 0.0;  //Z=19486
                double qqnn = 1.0;  //Z=19487
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19488
                    qqnn = qqnn*q*q;  //Z=19489
                    pqr += params.CR->carr4f[nser]*qqnn;  //Z=19490
                    const double delser = fabs((pqr-oldpqsum)/pqr);  //Z=19491
                    if ( delser<0.0001 ) break; /* goto 71; */  //Z=19492
                    oldpqsum = pqr;  //Z=19493
                }/*5*/  //Z=19494
                /*71:*/  //Z=19495
                //pqr = pqsum;  //Z=19496
            }/*4*/  //Z=19497
            else
            {/*4*/  //Z=19498
                const double argpq = q*params.radius/(zr+1);  //Z=19499
                pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1.0+argpq*argpq,zr/2.0);  //Z=19500
            }/*4*/  //Z=19501
            /*formfq:=*/ return pql*pqr*pqr;  //Z=19502
            /* formfq:=pql;;  //Z=19503 */
        }/*3*/ /*  of homogeneous  */  //Z=19504

        /*  core/shell disk  */  //Z=19506
        if ( params.cs==1 )
        {/*3*/  //Z=19507
            const double ccc1 = sqr(1-params.rho)*pow(params.p1,2);  //Z=19508
            const double ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,1);  //Z=19509
            const double ccc3 = sqr(params.rho);  //Z=19510
            const double vv3 = sqr((1-params.rho)*pow(params.p1,1)+params.rho);  //Z=19511

            //zz = zr;  // TODO: zz war in diesem Zweig nicht gesetzt
            const double argq = q*radiusm/(zr+1);  //Z=19513
            const double argpq = q*params.radius/(zr+1);  //Z=19514

            double F121, F122, F123;

            /*  F121 disk  */  //Z=19516
            if ( q<(0.8*params.limq4f) )
            {/*4*/  //Z=19517
                /* ** series expansion ** */  //Z=19518
                double pqsum = 1.0;  //Z=19519
                double oldpqsum = 0.0;  //Z=19520
                double qqnn = 1.0;  //Z=19521
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19522
                    qqnn = qqnn*q*q;  //Z=19523
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=19524
                    const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19525
                    if ( delser<0.0001 ) break; /* goto 72; */  //Z=19526
                    oldpqsum = pqsum;  //Z=19527
                }/*5*/  //Z=19528
                /*72:*/  //Z=19529
                F121 = ccc1*pqsum/vv3;  //Z=19530
            }/*4*/  //Z=19531
            else
            {/*4*/  //Z=19532
                const double pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1.0+argpq*argpq,zr/2.0);  //Z=19533
                F121 = ccc1*pqr*pqr/vv3;  //Z=19534
            }/*4*/  //Z=19535

            /*  F122 disk  */  //Z=19537
            if ( q<(1.0*params.limq5f) )
            {/*4*/  //Z=19538
                /* ** series expansion ** */  //Z=19539
                double pqsum = 1.0;  //Z=19540
                double oldpqsum = 0.0;  //Z=19541
                double qqnn = 1.0;  //Z=19542
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19543
                    qqnn = qqnn*q*q;  //Z=19544
                    pqsum = pqsum+params.CR->carr5f[nser]*qqnn;  //Z=19545
                    const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19546
                    if ( delser<0.0001 ) break; /* goto 73; */  //Z=19547
                    oldpqsum = pqsum;  //Z=19548
                }/*5*/  //Z=19549
                /*73:*/  //Z=19550
                F122 = ccc2*pqsum/vv3;  //Z=19551
            }/*4*/  //Z=19552
            else
            {/*4*/  //Z=19553
                const double pqr = (1/zr)*pow(argpq,-1)*sin(zr*atan(argpq))/pow(1.0+argpq*argpq,zr/2.0);  //Z=19554
                const double pqr1 = (1/zr)*pow(argq,-1)*sin(zr*atan(argq))/pow(1.0+argq*argq,zr/2.0);  //Z=19555
                F122 = ccc2*pqr*pqr1/vv3;  //Z=19556
            }/*4*/  //Z=19557

            /*  F123 disk  */  //Z=19559
            if ( q<(0.3*params.limq6f) )
            {/*4*/  //Z=19560
                /* ** series expansion ** */  //Z=19561
                double pqsum = 1.0;  //Z=19562
                double oldpqsum = 0.0;  //Z=19563
                double qqnn = 1.0;  //Z=19564
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=19565
                    qqnn = qqnn*q*q;  //Z=19566
                    pqsum = pqsum+params.CR->carr6f[nser]*qqnn;  //Z=19567
                    const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19568
                    if ( delser<0.0001 ) break; /* goto 74; */  //Z=19569
                    oldpqsum = pqsum;  //Z=19570
                }/*5*/  //Z=19571
                /*74:*/  //Z=19572
                F123 = ccc3*pqsum/vv3;  //Z=19573
            }/*4*/  //Z=19574
            else
            {/*4*/  //Z=19575
                const double pqr = (1/zr)*pow(argq,-1)*sin(zr*atan(argq))/pow(1.0+argq*argq,zr/2.0);  //Z=19576
                F123 = ccc3*pqr*pqr/vv3;  //Z=19577
                /*  add more terms, if necessary  */  //Z=19578
            }/*4*/  //Z=19579
            /*formfq:=*/ return pql*(F121+F122+F123);  //Z=19580
            /* formfq:=(F121+F122+F123);  //Z=19581 */
        }/*3*/ /*  of core/shell-disk  */  //Z=19582

        /* ** inhomogeneous core/shell disk ** */  //Z=19584
        if ( params.cs==2 )
        {/*3*/  //Z=19585

            const double dim = 1;  //Z=19587
            const double delc = 0.0001;  //Z=19588
            const double zz = zr;  //Z=19589
            const double xrad = q*radiusm;  //Z=19590
            const double xradp = q*params.radius;  //Z=19591
            const double x1z = q*params.radius/(2.0*(zz+1));  //Z=19592
            const double x12z = x1z*x1z;  //Z=19593
            const double x2z = q*radiusm/(2.0*(zz+1));  //Z=19594
            const double x22z = x2z*x2z;  //Z=19595

            const double lim = 18*exp(-5*params.sigma);  //Z=19597
            const double lim1 = lim;  //Z=19598
            //lim2 = lim*0.7;  //Z=19599
            //lim3 = lim;  //Z=19600
            const double lim4 = lim;  //Z=19601
            //lim5 = lim*0.7;  //Z=19602
            const double lim6 = lim*1.2;  //Z=19603

            const double a1 = (dim-params.alphash1)/2.0;  //Z=19605
            const double b1 = dim/2.0;  //Z=19606
            const double b2 = (dim+2-params.alphash1)/2.0;  //Z=19607
            const double b1s = (dim+2)/2.0;  //Z=19608
            const double v = -b1s+1/2.0;  //Z=19609
            const double c = a1-b1-b2+1/2.0;  //Z=19610
            //d0 = 1;  //Z=19611
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=19612
            const double e0 = 1.0;  //Z=19613
            const double e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=19614
            const double ee0 = 1.0;  //Z=19615
            const double ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=19616

            const double gb1s = 3*sqrt(M_PI)/4.0;  //Z=19618
            const double pz2v = 1/(zr*(zr-1));  //Z=19619
            const double pz2v1 = pz2v/(zr-2);  //Z=19620
            //pz2v2 = pz2v1/(zr-3);  //Z=19621

            const double gz1 = gamma(zr+1);  //Z=19623
            const double preg1 = gb1s/sqrt(M_PI);  //Z=19624
            const double preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=19625
            const double preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=19626
            //pzvc = gamma(zr+1+v+c)/gz1;  //Z=19627
            //pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=19628
            //pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=19629
            //pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=19630
            //pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=19631
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=19632
            const double pzc = gamma(zr+1+c)/gz1;  //Z=19633
            const double pzc1 = gamma(zr+1+c-1)/gz1;  //Z=19634
            const double pza = gamma(zr+1-2*a1)/gz1;  //Z=19635
            //pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=19636
            //pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19637
            //dnv0 = 1;  //Z=19638
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=19639
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19640
            //pva0 = gamma(zr+1-2*a1)/gz1;  //Z=19641

            const double cc1 = 1.0/dim;  //Z=19643
            const double cc4 = params.rho/((dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=19644
            const double cc6 = -params.rho/(dim-params.alphash1);  //Z=19645
            const double sumc = cc1+cc4+cc6;  //Z=19646

            double F12, F42, F62;

            /*  term #1 series  */  //Z=19648
            if ( (xradp)<lim1 )
            {/*4*/  //Z=19649
                //z12v[0] = 1;  //Z=19650
                //b1sv[0] = 1;  //Z=19651
                //fkv[0] = 1;  //Z=19652
                double qqnn = 1.0;  //Z=19653
                F12 = 1.0;  //Z=19654
                double oldF12sez = 1.0;  //Z=19655
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=19656
                    qqnn = qqnn*q*q;  //Z=19657
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19658
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19659
                    //fkv[n] = fkv[n-1]*n;  //Z=19660
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]);  //Z=19661 */

                    F12 += params.CR->carr4f[n]*qqnn;  //Z=19663

                    const double del = fabs((F12-oldF12sez)/F12);  //Z=19665
                    if ( del<delc ) break; /* goto 121; */  //Z=19666
                    oldF12sez = F12;  //Z=19667
                }/*5*/  //Z=19668
                /*121:*/  //Z=19669
                //F12 = F12sez;  //Z=19670
            }/*4*/  //Z=19671

            /*  term #4 series  */  //Z=19673
            if ( xradp<lim4 )
            {/*4*/  //Z=19674
                //z12v[0] = 1;  //Z=19675
                //a1v[0] = 1;  //Z=19676
                //b1v[0] = 1;  //Z=19677
                //b2v[0] = 1;  //Z=19678
                //b1sv[0] = 1;  //Z=19679
                //fkv[0] = 1;  //Z=19680
                double qqnn = 1.0;  //Z=19681
                F42 = 1.0;  //Z=19682
                double oldF42sez = 1.0;  //Z=19683
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=19684
                    qqnn = qqnn*q*q;  //Z=19685
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19686
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19687
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19688
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19689
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19690
                    //fkv[n] = fkv[n-1]*n;  //Z=19691
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19692 */

                    F42 += params.CR->carr5f[n]*qqnn;  //Z=19694

                    const double del = fabs((F42-oldF42sez)/F42);  //Z=19696
                    if ( del<delc ) break; /* goto 124; */  //Z=19697
                    oldF42sez = F42;  //Z=19698
                }/*5*/  //Z=19699
                /*124:*/  //Z=19700
                //F42 = F42sez;  //Z=19701
            }/*4*/  //Z=19702

            /*  term #6 series  */  //Z=19704
            if ( xradp<lim6 )
            {/*4*/  //Z=19705
                //z12v[0] = 1;  //Z=19706
                //a1v[0] = 1;  //Z=19707
                //b1v[0] = 1;  //Z=19708
                //b2v[0] = 1;  //Z=19709
                //b1sv[0] = 1;  //Z=19710
                //fkv[0] = 1;  //Z=19711
                double qqnn = 1.0;  //Z=19712
                F62 = 1.0;  //Z=19713
                double oldF62sez = 1.0;  //Z=19714
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=19715
                    qqnn = qqnn*q*q;  //Z=19716
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19717
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19718
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19719
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19720
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19721
                    //fkv[n] = fkv[n-1]*n;  //Z=19722
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19723 */

                    F62 += params.CR->carr6f[n]*qqnn;  //Z=19725

                    const double del = fabs((F62-oldF62sez)/F62);  //Z=19727
                    if ( del<delc ) break; /* goto 126; */  //Z=19728
                    oldF62sez = F62;  //Z=19729
                }/*5*/  //Z=19730
                /*126:*/  //Z=19731
                //F62 = F62sez;  //Z=19732
            }/*4*/  //Z=19733

            /* ** term #1 asymptote ** */  //Z=19735
            if ( xradp>=lim1 )
            {/*4*/  //Z=19736
                const double arg11 = (zr+v+1)*atan(2.0*x1z);  //Z=19737
                const double nen11 = pow(1.0+4*x1z*x1z,(zr+v+1)/2.0);  //Z=19738
                const double arg12 = (zr+v)*atan(2.0*x1z);  //Z=19739
                const double nen12 = pow(1.0+4*x1z*x1z,(zr+v)/2.0);  //Z=19740

                const double F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);  //Z=19742
                const double F12as2z = ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);  //Z=19743
                F12 = preg1*pow(x1z,v)*(F12as1z+F12as2z);  //Z=19744
                //F12 = F12asz;  //Z=19745
            }/*4*/  //Z=19746

            /* ** term #4 asymptote ** */  //Z=19748
            if ( xrad>=lim4 )
            {/*4*/  //Z=19749
                const double F42as10z = preg4*pow(x22z,-a1);  //Z=19750
                //F42as1sumz = pva0;  //Z=19751
                //F42as1z = F42as10z*F42as1sumz;  //Z=19752
                const double F42as1z0 = F42as10z*pza;   /* * */  //Z=19753

                const double F42as40z = preg3*pow(x2z,c);  //Z=19755
                const double arg44 = (zr+c+1)*atan(2.0*x2z);  //Z=19756
                const double nen44 = pow(1.0+4*x2z*x2z,(zr+c+1)/2.0);  //Z=19757
                const double arg45 = (zr+c)*atan(2.0*x2z);  //Z=19758
                const double nen45 = pow(1.0+4*x2z*x2z,(zr+c)/2.0);  //Z=19759
                const double F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);  //Z=19760
                const double F42as28 = e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);  //Z=19761
                const double F42as4z = F42as40z*(F42as27+F42as28);  //Z=19762
                //F42asz = F42as1z+F42as4z;  //Z=19763
                F42 = F42as1z0+F42as4z;  //Z=19764
                //F42 = F42asz0;  //Z=19765
            }/*4*/  //Z=19766

            /* ** term #6 asymptote ** */  //Z=19768
            if ( xradp>=lim6 )
            {/*4*/  //Z=19769
                const double F62as10z = preg4*pow(x12z,-a1);  //Z=19770
                //F62as1sumz = pva0;  //Z=19771
                //F62as1z = F62as10z*F62as1sumz;  //Z=19772
                const double F62as1z0 = F62as10z*pza;     /* * */  //Z=19773

                const double F62as40z = preg3*pow(x1z,c);  //Z=19775
                const double arg64 = (zr+c+1)*atan(2.0*x1z);  //Z=19776
                const double nen64 = pow(1.0+4*x1z*x1z,(zr+c+1)/2.0);  //Z=19777
                const double arg65 = (zr+c)*atan(2.0*x1z);  //Z=19778
                const double nen65 = pow(1.0+4*x1z*x1z,(zr+c)/2.0);  //Z=19779
                const double F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);  //Z=19780
                const double F62as28 = e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);  //Z=19781
                const double F62as4z = F62as40z*(F62as27+F62as28);  //Z=19782
                //F62asz = F62as1z+F62as4z;  //Z=19783
                F62 = F62as1z0+F62as4z;  //Z=19784
                //F62 = F62asz0;  //Z=19785
            }/*4*/  //Z=19786

            const double FF1 = (cc1*F12+cc4*F42+cc6*F62)/sumc;  //Z=19788
            /* FF1:=(cc1*F12)/sumc;  //Z=19789 */

            /*formfq:=*/ return FF1*FF1;  //Z=19791

            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,1,sigmar,q);  //Z=19793 */

        }/*3*/ /*  of inhomogeneous core/shell disk  */  //Z=19795
    //}/*2*/ /*  of disk  */  //Z=19796

    return 0.0;
}  //Z=19829


#endif // SC_LIB_FORMFQ_partDisk_H
