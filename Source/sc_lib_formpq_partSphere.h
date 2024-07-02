#ifndef SC_LIB_FORMPQ_partSphere_H
#define SC_LIB_FORMPQ_partSphere_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partSphere(double q) const   /*Z=14910*/
{/*1*/  //Z=15188

    const double z  /*zl*/ = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=15231
    const double zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232
    const double radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=15233

    // Bei sigma=0.5 wird zr=3 und somit das Ergebnis hier zum Großteil inf ...

    const double zz = zr;
    //double z  = zl;

    //ella = params.radius;  //Z=15235
    //ellb = params.length;  //Z=15236
    //ellc = radiusm;  //Z=15237

    CHECKENDTHREAD_VAL;

    /* ************ */  //Z=15240
    /* ** sphere ** */  //Z=15241
    /* ************ */  //Z=15242
    //if ( params.part==0 )
    //{/*2*/  //Z=15243
        /* ** homogeneous sphere ** */  //Z=15244
        if ( params.cs==0 )
        {/*3*/  //Z=15245
            //if ( q > 2.2 ) qDebug() << "  formpq" << 0.4*params.limq4 << q << "r"<<params.radius;
            //[15386]        //if (q<1.5*limq4) then begin    (* for monodisperse sphere test *)
            if ( q<0.4*params.limq4 )   // alt: 1.5, neu:0.4
            {/*4*/  //Z=15246
                double pqsum = 1.0;  //Z=15247
                double oldpqsum = 0.0;  //Z=15248
                /* qqn[0]:=1.0;  //Z=15249 */
                double qq2 = 1.0;  //Z=15250
                for ( int nser=1; nser<=100; nser++ )
                {/*5*/  //Z=15251
                    /* qqn[nser]:=qqn[nser-1]*q*q;  //Z=15252 */
                    qq2 = qq2*q*q;  //Z=15253
                    /* pqsum:=pqsum+carr4p[nser]*qqn[nser];  //Z=15254 */
                    pqsum = pqsum+params.CR->carr4p[nser]*qq2;  //Z=15255
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15256
                    if ( delser<0.0001 ) break; /* goto 50; */  //Z=15257
                    oldpqsum = pqsum;  //Z=15258
                }/*5*/  //Z=15259
                /*50:*/  //Z=15260
                /*formpq:=*/ return pqsum;  //Z=15261
            }/*4*/  //Z=15262
            else
            {/*4*/  //Z=15263
                const double argq = q*params.radius/(zr+1);  //Z=15264
                double pqr;
                if ( zr*(zr-1)*(zr-2)*(zr-3)*(zr-4)*(zr-5) == 0 )
                    pqr = 0.1;
                else
                    pqr = (1/(2.0*zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4);  //Z=15265
                const double pq1 = pqr*(1+cos((zr-3)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-3)/2.0));  //Z=15266
                const double pq2 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-4)/2.0);  //Z=15267
                const double pq3 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2.0*argq))/pow(1.0+4*argq*argq,zr-5)/2.0);  //Z=15268
                // qDebug() << argq << pqr << pq1 << pq2 << pq3 << "zr" << zr;
                // Debug: 1.74127 inf inf inf inf zr 3  <<== bei sigma=0.5, Ergebnis = inf bei vielen q am Rand (OHNE DAS IF OBEN)
                /*formpq:=*/ return 9.0*(pq1-2.0*pq2+pq3);  //Z=15269
            }/*4*/  //Z=15270
        }/*3*/ /*  of homogeneous sphere */  //Z=15271

        /* ** core/shell sphere ** */  //Z=15273
        if ( params.cs==1 )
        {/*3*/  //Z=15274

            const double cc1 = sqr(params.rho);  //Z=15276
            const double cc2 = 2*params.p1*params.rho*(1-params.rho);  //Z=15277
            const double cc3 = sqr(1-params.rho)*sqr(params.p1);  //Z=15278
            const double cc4 = -2*sqr(params.rho);  //Z=15279
            const double cc5 = -2*params.p1*params.rho*(1-params.rho);  //Z=15280
            const double cc6 = sqr(params.rho);  //Z=15281
            const double cc7 = -2*params.rho*(1-params.rho);  //Z=15282
            const double cc8 = -sqr(1-params.rho)*2*params.p1;  //Z=15283
            const double cc9 = 2*params.rho*(1-params.rho);  //Z=15284
            const double cc10 = sqr(1-params.rho);  //Z=15285

            const double ccc1 = sqr(1-params.rho)*pow(params.p1,6);  //Z=15287
            const double ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,3);  //Z=15288
            const double ccc3 = params.rho*params.rho;  //Z=15289
            const double vv3 = sqr((1-params.rho)*pow(params.p1,3)+params.rho);  //Z=15290

            const double argq = q*radiusm/(zz+1);  //Z=15292
            const double argpq = q*params.radius/(zz+1);  //Z=15293
            const double pqr = (1/(2.0*zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4);  //Z=15294

            double F121, F122, F123;

            /*  F121 sphere  */  //Z=15296
            if ( q<(0.5*params.limq4) ) //20240301 - war 0.3
            {/*4*/  //Z=15297
                double qqnn = 1.0;  //Z=15298
                double pqsum = 1.0;  //Z=15299
                double oldpqsum = 0.0;  //Z=15300
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=15301
                    qqnn = qqnn*q*q;  //Z=15302
                    pqsum = pqsum+qqnn*params.CR->carr4p[nser];  //Z=15303
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15304
                    if ( delser<0.0001 ) break; /* goto 51; */  //Z=15305
                    oldpqsum = pqsum;  //Z=15306
                }/*5*/  //Z=15307
                /*51:*/  //Z=15308
                F121 = ccc1*pqsum/vv3;  //Z=15309
            }/*4*/  //Z=15310
            else
            {/*4*/  //Z=15311
                const double ac3 = pqr*(1+cos((zr-3)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-3)/2.0));  //Z=15312
                const double ac8 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-4)/2.0);  //Z=15313
                const double ac10 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-5)/2.0));  //Z=15314
                F121 = 9*(cc3*ac3+cc8*ac8+cc10*ac10)/vv3;  //Z=15315
            }/*4*/  //Z=15316

            /*  F122 sphere  */  //Z=15318
            if ( q<(0.5*params.limq5) ) //20240301 - war 0.3
            {/*4*/  //Z=15319
                double qqnn = 1.0;  //Z=15320
                double pqsum = 1.0;  //Z=15321
                double oldpqsum = 0.0;  //Z=15322
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=15323
                    qqnn = qqnn*q*q;  //Z=15324
                    pqsum = pqsum+qqnn*params.CR->carr5p[nser];  //Z=15325
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15326
                    if ( delser<0.0001 ) break; /* goto 52; */  //Z=15327
                    oldpqsum = pqsum;  //Z=15328
                }/*5*/  //Z=15329
                /*52:*/  //Z=15330
                F122 = ccc2*pqsum/vv3;  //Z=15331
            }/*4*/  //Z=15332
            else
            {/*4*/  //Z=15333
                const double argbm = (zr-3)*atan(argpq-argq);  //Z=15334
                const double nenbm = pow(1.0+sqr(argpq-argq),(zr-3)/2.0);  //Z=15335
                const double argbp = (zr-3)*atan(argpq+argq);  //Z=15336
                const double nenbp = pow(1.0+sqr(argpq+argq),(zr-3)/2.0);  //Z=15337
                const double ac2 = pqr*(cos(argbm)/nenbm+cos(argbp)/nenbp);  //Z=15338
                const double argep = (zr-4)*atan(argpq+argq);  //Z=15339
                const double nenep = pow(1.0+sqr(argpq+argq),(zr-4)/2.0);  //Z=15340
                const double argem = (zr-4)*atan(argpq-argq);  //Z=15341
                const double nenem = pow(1.0+sqr(argpq-argq),(zr-4)/2.0);  //Z=15342
                const double ac5 = (pqr/((zr-4)*argq))*(sin(argep)/nenep-sin(argem)/nenem);  //Z=15343
                const double arggp = (zr-4)*atan(argpq+argq);  //Z=15344
                const double nengp = pow(1.0+sqr(argpq+argq),(zr-4)/2.0);  //Z=15345
                const double arggm = (zr-4)*atan(argq-argpq);  //Z=15346
                const double nengm = pow(1.0+sqr(argq-argpq),(zr-4)/2.0);  //Z=15347
                const double ac7 = (pqr/((zr-4)*argq))*(sin(arggp)/nengp-sin(arggm)/nengm);  //Z=15348
                //argim = (zr-5)*atan(argpq-argq);  //Z=15349
                //nenim = pow(1.0+sqr(argpq-argq),(zr-5)/2.0);  //Z=15350
                //argip = (zr-5)*atan(argpq+argq);  //Z=15351
                //nenip = pow(1.0+sqr(argpq+argq),(zr-5)/2.0);  //Z=15352
                const double ac9 = (pqr/((zr-4)*(zr-5)*argq*argq))*(cos(argbm)/nenbm-cos(argbp)/nenbp);  //Z=15353

                F122 = 9*(cc2*ac2+cc5*ac5+cc7*ac7+cc9*ac9)/vv3;  //Z=15355
            }/*4*/  //Z=15356

            /*  F123 sphere  */  //Z=15358
            if ( q<(0.5*params.limq6) ) //20240301 - war 0.3
            {/*4*/  //Z=15359
                double qqnn = 1.0;  //Z=15360
                double pqsum = 1.0;  //Z=15361
                double oldpqsum = 0.0;  //Z=15362
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=15363
                    qqnn = qqnn*q*q;  //Z=15364
                    pqsum = pqsum+qqnn*params.CR->carr6p[nser];  //Z=15365
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=15366
                    if ( delser<0.0001 ) break; /* goto 53; */  //Z=15367
                    oldpqsum = pqsum;  //Z=15368
                }/*5*/  //Z=15369
                /*53:*/  //Z=15370
                F123 = ccc3*pqsum/vv3;  //Z=15371
            }/*4*/  //Z=15372
            else
            {/*4*/  //Z=15373
                const double ac1 = pqr*(1+cos((zr-3)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-3)/2.0));  //Z=15374
                const double ac4 = (pqr/((zr-4)*argq))*sin((zr-4)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-4)/2.0);  //Z=15375
                const double ac6 = (pqr/((zr-4)*(zr-5)*argq*argq))*(1-cos((zr-5)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-5)/2.0));  //Z=15376
                F123 = 9*(cc1*ac1+cc4*ac4+cc6*ac6)/vv3;  //Z=15377
            }/*4*/  //Z=15378
            /*formpq:=*/ return F121+F122+F123;  //Z=15379
        }/*3*/  /*  of core/shell sphere  */  //Z=15380

        /* ** inhomogeneous core/shell sphere ** */  //Z=15382
        if ( params.cs==2 )
        {/*3*/  //Z=15383

            const double dim = 3;  //Z=15385
            const double delc = 0.0001;  //Z=15386
            const double xrad = q*radiusm;  //Z=15387
            const double xradp = q*params.radius;  //Z=15388
            const double x1z = q*params.radius/(2.0*(zr+1));  //Z=15389
            const double x12z = x1z*x1z;  //Z=15390
            const double x2z = q*radiusm/(2.0*(zr+1));  //Z=15391
            const double x22z = x2z*x2z;  //Z=15392

            const double lim = 18*exp(-5*params.sigma);  //Z=15394
            const double lim1 = lim;  //Z=15395
            const double lim2 = lim*0.7;  //Z=15396
            const double lim3 = lim;  //Z=15397
            const double lim4 = lim;  //Z=15398
            const double lim5 = lim*0.7;  //Z=15399
            const double lim6 = lim*1.2;  //Z=15400

            const double a1 = (dim-params.alphash1)/2.0;  //Z=15402
            const double b1 = dim/2.0;  //Z=15403
            const double b2 = (dim+2-params.alphash1)/2.0;  //Z=15404
            const double b1s = (dim+2)/2.0;  //Z=15405
            const double v = -b1s+1/2.0;  //Z=15406
            const double c = a1-b1-b2+1/2.0;  //Z=15407
            const double d0 = 1;  //Z=15408
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=15409
            const double e0 = 1.0;  //Z=15410
            const double e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=15411
            const double ee0 = 1.0;  //Z=15412
            const double ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=15413

            const double gb1s = 3*sqrt(M_PI)/4.0;  //Z=15415
            const double pz2v = 1/(zr*(zr-1)*(zr-2)*(zr-3));  //Z=15416
            const double pz2v1 = pz2v/(zr-4);  //Z=15417
            const double pz2v2 = pz2v1/(zr-5);  //Z=15418

            const double gz1 = gamma(zr+1);  //Z=15420
            const double preg1 = gb1s/sqrt(M_PI);  //Z=15421
            const double preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=15422
            const double preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=15423
            const double pzvc = gamma(zr+1+v+c)/gz1;  //Z=15424
            const double pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=15425
            const double pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=15426
            const double pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=15427
            const double pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=15428
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=15429
            const double pzc = gamma(zr+1+2*c)/gz1;  //Z=15430
            const double pzc1 = gamma(zr+1+2*c-1)/gz1;  //Z=15431
            const double pza = gamma(zr+1-4*a1)/gz1;  //Z=15432
            const double pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=15433
            const double pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=15434
            //dnv0 = 1;  //Z=15435
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=15436
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=15437
            //pva0 = gamma(zr+1-4*a1)/gz1;  //Z=15438

            const double cc1 = 1/(dim*dim);  //Z=15440
            const double cc2 = 2*params.rho/(dim*(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=15441
            const double cc3 = -2*params.rho/(dim*(dim-params.alphash1));  //Z=15442
            const double cc4 = sqr(params.rho)/(sqr(dim-params.alphash1)*pow(sqr(params.p1),dim-params.alphash1));  //Z=15443
            const double cc5 = -2*sqr(params.rho)/(sqr(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=15444
            const double cc6 = sqr(params.rho)/sqr(dim-params.alphash1);  //Z=15445
            const double vv3 = cc1+cc2+cc3+cc4+cc5+cc6;  //Z=15446

            double F12, F22, F32, F42, F52, F62;

            /*  term #1 series  */  //Z=15448
            if ( xradp<lim1 )
            {/*4*/  //Z=15449
                //z12v[0] = 1;  //Z=15450
                //b1sv[0] = 1;  //Z=15451
                //fkv[0] = 1;  //Z=15452
                //gam3[0] = sqrt(M_PI)/2.0;  //Z=15453
                double qqnn = 1.0; // qqn[0] = 1.0;  //Z=15454
                F12 = 1.0;  //Z=15455
                double oldF12 = 0.0;  //Z=15456
                for ( int n=1; n<=120; n++ )
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

                    F12 += params.CR->carr1p[n]*qqnn; //qqn[n];  //Z=15468

                    double del = fabs((F12-oldF12)/F12);  //Z=15470
                    if ( del<delc ) break; /* goto 201; */  //Z=15471
                    oldF12 = F12;  //Z=15472
                }/*5*/  //Z=15473
                /*201:*/  //Z=15474
                //F12 = F12sez;  //Z=15475
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
                F22 = 1.0;  //Z=15488
                double oldF22 = 0.0;  //Z=15489
                for ( int n=1; n<=120; n++ )
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

                    F22 += params.CR->carr2p[n]*qqnn; // qqn[n];  //Z=15506

                    double del = fabs((F22-oldF22)/F22);  //Z=15508
                    if ( del<delc ) break; /* goto 202; */  //Z=15509
                    oldF22 = F22;  //Z=15510
                }/*5*/  //Z=15511
                /*202:*/  //Z=15512
                //F22 = F22sez;  //Z=15513
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
                F32 = 1.0;  //Z=15526
                double oldF32 = 0.0;  //Z=15527
                for ( int n=1; n<=120; n++ )
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

                    F32 += params.CR->carr3p[n]*qqnn; //qqn[n];  //Z=15544

                    double del = fabs((F32-oldF32)/F32);  //Z=15546
                    if ( del<delc ) break; /* goto 203; */  //Z=15547
                    oldF32 = F32;  //Z=15548
                }/*5*/  //Z=15549
                /*203:*/  //Z=15550
                //F32 = F32sez;  //Z=15551
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
                F42 = 1.0;  //Z=15564
                double oldF42 = 0.0;  //Z=15565
                for ( int n=1; n<=120; n++ )
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

                    F42 += params.CR->carr4p[n]*qqnn; //qqn[n];  //Z=15582

                    double del = fabs((F42-oldF42)/F42);  //Z=15584
                    if ( del<delc ) break; /* goto 204; */  //Z=15585
                    oldF42 = F42;  //Z=15586
                }/*5*/  //Z=15587
                /*204:*/  //Z=15588
                //F42 = F42sez;  //Z=15589
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
                F52 = 1.0;  //Z=15601
                double oldF52 = 0.0;  //Z=15602
                for ( int n=1; n<=120; n++ )
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

                    F52 += params.CR->carr5p[n]*qqnn; //qqn[n];  //Z=15618

                    double del = fabs((F52-oldF52)/F52);  //Z=15620
                    if ( del<delc ) break; /* goto 205; */  //Z=15621
                    oldF52 = F52;  //Z=15622
                }/*5*/  //Z=15623
                /*205:*/  //Z=15624
                //F52 = F52sez;  //Z=15625
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
                F62 = 1.0;  //Z=15637
                double oldF62 = 0.0;  //Z=15638
                for ( int n=1; n<=120; n++ )
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

                    F62 += params.CR->carr6p[n]*qqnn; //qqn[n];  //Z=15654

                    double del = fabs((F62-oldF62)/F62);  //Z=15656
                    if ( del<delc ) break; /* goto 206; */  //Z=15657
                    oldF62 = F62;  //Z=15658
                }/*5*/  //Z=15659
                /*206:*/  //Z=15660
                //F62 = F62sez;  //Z=15661
            }/*4*/  //Z=15662


            /* ** term #1 asymptote ** */  //Z=15665
            if ( xradp>=lim1 )
            {/*4*/  //Z=15666
                const double arg11 = (zr+2*v+1)*atan(4.0*x1z);  //Z=15667
                const double nen11 = pow(1.0+16*x1z*x1z,(zr+2*v+1)/2.0);  //Z=15668
                const double arg12 = (zr+2*v)*atan(4.0*x1z);  //Z=15669
                const double nen12 = pow(1.0+16*x1z*x1z,(zr+2*v)/2.0);  //Z=15670
                const double arg13 = (zr+2*v-1)*atan(4.0*x1z);  //Z=15671
                const double nen13 = pow(1.0+16*x1z*x1z,(zr+2*v-1)/2.0);  //Z=15672

                const double F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);  //Z=15674
                const double F12as2z = 2*ee0*ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);  //Z=15675
                const double F12as3z = ee1*ee1*(1/(4.0*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);  //Z=15676
                F12 = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);  //Z=15677
                //F12 = F12asz;  //Z=15678
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
                const double F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);  //Z=15689
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);  //Z=15690

                const double arg210 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=15692
                const double nen210 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=15693
                const double arg220 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=15694
                const double nen220 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=15695
                const double F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);  //Z=15696
                const double F22as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);  //Z=15697
                const double F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);  //Z=15698
                const double arg23 = (zr+v+c+1)*atan(2.0*(x1z-x2z));  //Z=15699
                const double nen23 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);  //Z=15700
                const double arg24 = (zr+v+c+1)*atan(2.0*(x1z+x2z));  //Z=15701
                const double nen24 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);  //Z=15702
                const double arg25 = (zr+v+c)*atan(2.0*(x1z-x2z));  //Z=15703
                const double nen25 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);  //Z=15704
                const double arg26 = (zr+v+c)*atan(2.0*(x1z+x2z));  //Z=15705
                const double nen26 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);  //Z=15706
                const double arg27 = (zr+v+c-1)*atan(2.0*(x1z-x2z));  //Z=15707
                const double nen27 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);  //Z=15708
                const double arg28 = (zr+v+c-1)*atan(2.0*(x1z+x2z));  //Z=15709
                const double nen28 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);  //Z=15710

                const double a22as21z = (1/2.0)*ee0*e0*pzvc;  //Z=15712
                const double F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);  //Z=15713
                const double a22as22z = (1/2.0)*ee0*e1*(1/(2.0*x2z))*pzvc1;  //Z=15714
                const double F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);  //Z=15715
                const double a22as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1;  //Z=15716
                const double F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);  //Z=15717
                const double a22as24z = (1/2.0)*ee1*e1*(1/(2.0*x1z))*(1/(2.0*x2z))*pzvc2;  //Z=15718
                const double F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);  //Z=15719
                const double F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);  //Z=15720
                const double F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);  //Z=15721
                //F22asz = F22as1z+F22as2z;  //Z=15722
                F22 = F22as1z0+F22as2z;  //Z=15723
                //F22 = F22asz0;  //Z=15724
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
                const double F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);  //Z=15735
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);  //Z=15736

                const double arg310 = (z+v-2*a1+1)*atan(2.0*x1z);  //Z=15738
                const double nen310 = pow(1.0+4*x1z*x1z,(z+v-2*a1+1)/2.0);  //Z=15739
                const double arg320 = (z+v-2*a1)*atan(2.0*x1z);  //Z=15740
                const double nen320 = pow(1.0+4*x1z*x1z,(z+v-2*a1)/2.0);  //Z=15741
                const double F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);  //Z=15742
                const double F32as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);  //Z=15743
                const double F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);  //Z=15744

                const double arg33 = (zr+v+c+1)*atan(4.0*x1z);  //Z=15746
                const double nen33 = pow(1.0+16*x1z*x1z,(zr+v+c+1)/2.0);  //Z=15747
                const double arg34 = (zr+v+c)*atan(4.0*x1z);  //Z=15748
                const double nen34 = pow(1.0+16*x1z*x1z,(zr+v+c)/2.0);  //Z=15749
                const double arg35 = (zr+v+c-1)*atan(4.0*x1z);  //Z=15750
                const double nen35 = pow(1.0+16*x1z*x1z,(zr+v+c-1)/2.0);  //Z=15751
                const double F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);  //Z=15752
                const double F32as22z = (1/2.0)*ee0*e1*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);  //Z=15753
                const double F32as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);  //Z=15754
                const double F32as24z = (1/2.0)*ee1*e1*(1/(4.0*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);  //Z=15755
                const double F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);  //Z=15756
                const double F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);  //Z=15757
                //F32asz = F32as1z+F32as2z;  //Z=15758
                F32 = F32as1z0+F32as2z;  //Z=15759
                //F32 = F32asz0;  //Z=15760
            }/*4*/  //Z=15761

            /* ** term #4 asymptote ** */  //Z=15764
            if ( xrad>=lim4 )
            {/*4*/  //Z=15765
                const double F42as10z = preg4*preg4*pow(x22z,-2*a1);  //Z=15766
                //F42as1sumz = pva0;  //Z=15767
                //F42as1z = F42as10z*F42as1sumz;  //Z=15768
                const double F42as1z0 = F42as10z*pza;  //Z=15769

                const double arg41 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=15771
                const double nen41 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=15772
                const double arg42 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=15773
                const double nen42 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=15774
                //arg43 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=15775
                //nen43 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=15776
                const double F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=15777
                const double F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=15778
                const double F42as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);  //Z=15779
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);  //Z=15780
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);  //Z=15781
                const double F42as2z0 = F42as20z*(F42as21+F42as22);  //Z=15782

                const double F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=15784
                const double F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=15785
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);  //Z=15786
                const double F42as26 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);  //Z=15787
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);  //Z=15788
                const double F42as3z0 = F42as30z*(F42as24+F42as26);  //Z=15789

                const double F42as40z = preg3*preg3*pow(x2z*x2z,c);  //Z=15791
                const double arg44 = (zr+2*c+1)*atan(4.0*x2z);  //Z=15792
                const double nen44 = pow(1.0+16*x2z*x2z,(zr+2*c+1)/2.0);  //Z=15793
                const double arg45 = (zr+2*c)*atan(4.0*x2z);  //Z=15794
                const double nen45 = pow(1.0+16*x2z*x2z,(zr+2*c)/2.0);  //Z=15795
                const double F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);  //Z=15796
                const double F42as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=15797
                const double F42as29 = (1/2.0)*e1*e0*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=15798
                const double F42as4z = F42as40z*(F42as27+F42as28+F42as29);  //Z=15799
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;  //Z=15800
                F42 = F42as1z0+F42as2z0+F42as3z0+F42as4z;  //Z=15801
                //F42 = F42asz0;  //Z=15802
            }/*4*/  //Z=15803

            /* ** term #5 asymptote ** */  //Z=15806
            if ( xradp>=lim5 )
            {/*4*/  //Z=15807
                const double F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);  //Z=15808
                //F52as1sumz = pva0;  //Z=15809
                //F52as1z = F52as10z*F52as1sumz;  //Z=15810
                const double F52as1z0 = F52as10z*pza;  //Z=15811

                const double F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);  //Z=15813
                const double arg51 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=15814
                const double nen51 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=15815
                const double arg52 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=15816
                const double nen52 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=15817
                //arg53 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=15818
                //nen53 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=15819
                const double F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);  //Z=15820
                const double F52as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);  //Z=15821
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);  //Z=15822
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);  //Z=15823
                const double F52as2z0 = F52as20z*(F52as21+F52as22);  //Z=15824

                const double F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);  //Z=15826
                const double arg54 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=15827
                const double nen54 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=15828
                //arg55 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=15829
                //nen55 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=15830
                const double arg56 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=15831
                const double nen56 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=15832
                const double F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);  //Z=15833
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);  //Z=15834
                const double F52as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);  //Z=15835
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);  //Z=15836
                const double F52as3z0 = F52as30z*(F52as24+F52as26);  //Z=15837

                const double F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);  //Z=15839
                const double arg57 = (zr+2*c+1)*atan(2.0*(x1z-x2z));  //Z=15840
                const double nen57 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);  //Z=15841
                const double arg58 = (zr+2*c+1)*atan(2.0*(x1z+x2z));  //Z=15842
                const double nen58 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);  //Z=15843
                const double arg59 = (zr+2*c)*atan(2.0*(x1z-x2z));  //Z=15844
                const double nen59 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);  //Z=15845
                const double arg510 = (zr+2*c)*atan(2.0*(x1z+x2z));  //Z=15846
                const double nen510 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);  //Z=15847
                const double F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);  //Z=15848
                const double F52as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=15849
                const double F52as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=15850
                const double F52as4z = F52as40z*(F52as27+F52as28+F52as29);  //Z=15851
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;  //Z=15852
                F52 = F52as1z0+F52as2z0+F52as3z0+F52as4z;  //Z=15853
                //F52 = F52asz0;  //Z=15854
            }/*4*/  //Z=15855

            /* ** term #6 asymptote ** */  //Z=15857
            if ( xradp>=lim6 )
            {/*4*/  //Z=15858
                const double F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);  //Z=15859
                //F62as1sumz = pva0;  //Z=15860
                //F62as1z = F62as10z*F62as1sumz;  //Z=15861
                const double F62as1z0 = F62as10z*pza;  //Z=15862

                const double F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=15864
                const double arg61 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=15865
                const double nen61 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=15866
                const double arg62 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=15867
                const double nen62 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=15868
                //arg63 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=15869
                //nen63 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=15870
                const double F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=15871
                const double F62as22 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=15872
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);  //Z=15873
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);  //Z=15874
                const double F62as2z0 = F62as20z*(F62as21+F62as22);  //Z=15875

                const double F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=15877
                const double F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=15878
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);  //Z=15879
                const double F62as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=15880
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);  //Z=15881
                const double F62as3z0 = F62as30z*(F62as24+F62as26);  //Z=15882

                const double F62as40z = preg3*preg3*pow(x1z*x1z,c);  //Z=15884
                const double arg64 = (zr+2*c+1)*atan(4.0*x1z);  //Z=15885
                const double nen64 = pow(1.0+16*x1z*x1z,(zr+2*c+1)/2.0);  //Z=15886
                const double arg65 = (zr+2*c)*atan(4.0*x1z);  //Z=15887
                const double nen65 = pow(1.0+16*x1z*x1z,(zr+2*c)/2.0);  //Z=15888
                const double F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);  //Z=15889
                const double F62as28 = (1/2.0)*e0*e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=15890
                const double F62as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=15891
                const double F62as4z = F62as40z*(F62as27+F62as28+F62as29);  //Z=15892
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;  //Z=15893
                F62 = F62as1z0+F62as2z0+F62as3z0+F62as4z;  //Z=15894
                //F62 = F62asz0;  //Z=15895
            }/*4*/  //Z=15896

            /*formpq:=*/ return (cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;  //Z=15898

            /* formpq:=pqcoreshellin(1.0,rho,p1,1.0,0.001,alfa,radiusm,3,sigmar,q);  //Z=15901 */
        }/*3*/ /*  of inhomogeneous core/shell sphere  */  //Z=15902

        /*  myelin sphere  */  //Z=15906
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=15907

            /*  sphere parameters  */  //Z=15909
            const double v = -2;  //Z=15910
            const double e0 = 1;  //Z=15911
            const double e1 = -1;  //Z=15912
            const double preg1 = 3/4.0;  //Z=15913
            const double pz2v = 1/(zr*(zr-1)*(zr-2)*(zr-3));  //Z=15914
            const double pz2v1 = pz2v/(zr-4);  //Z=15915
            const double pz2v2 = pz2v1/(zr-5);  //Z=15916
            const double lim = 18*exp(-5*params.sigma);  //Z=15917
            const double lim1 = lim*1.4;  //Z=15918
            const double rad = params.CR->myarray[1];  //Z=15919
            const double inmax = round(params.CR->myarray[14]);  //Z=15920
            const double vvm = params.CR->myarray[15];  //Z=15921
            const double rmax = params.CR->myarray[16];  //Z=15922
            const double xmax = q*rmax;  //Z=15923

            if ( xmax<(lim1) )
            {/*4*/  //Z=15925
                /* fkv[0]:=1;  //Z=15926 */
                // Auf der GPU ist es effizienter, eine kleine Berechnung zu machen als ein grosses Array zu haben.
                double qqnn;  //Z=15927
                //for ( nser=1; nser<=120; nser++ )
                //{   //Z=15928
                //    qqn[nser] = qqn[nser-1]*q*q;  //Z=15929, muss bleiben!
                //    /* fkv[nser]:=fkv[nser-1]*nser;  //Z=15930 */
                //}   //Z=15931

                double F12sum = 0.0;  //Z=15933
                for ( int ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=15934
                    for ( int jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=15935
                        double F12 = 1.0;  //Z=15936
                        double oldF12 = 1.0;  //Z=15937
                        qqnn = 1.0;
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=15938
                            qqnn = qqnn * q*q;
                            double pqsum = 0;  //Z=15939
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=15940
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]);  //Z=15941 */
                                pqsum += pow(params.CR->carr7p[ii],2*mser)*pow(params.CR->carr7p[jj],2*(nser-mser))/(params.CR->carr6p[mser]*params.CR->carr6p[nser-mser]);  //Z=15942

                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=15944 */
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx];  //Z=15945 */
                            }/*8*/  //Z=15946
                            F12 += params.CR->carr4p[nser]*qqnn*pqsum;  //Z=15947
                            double delser = fabs((F12-oldF12)/F12);  //Z=15948
                            if ( delser<0.0001 ) break; /* goto 250; */  //Z=15949
                            oldF12 = F12;  //Z=15950
                        }/*7*/  //Z=15951
                        /*250:*/  //Z=15952
                        F12sum += params.CR->carr5p[ii]*params.CR->carr5p[jj]*F12;  //Z=15953
                    }/*6*/  //Z=15954
                }/*5*/  //Z=15955
                return /*F12 =*/ F12sum/vvm;  //Z=15956
                //F12 = F12ser;  //Z=15957
            }/*4*/  //Z=15958
            else
            {/*4*/  //Z=15959
                const double xrz = q*rad/(zr+1);  //Z=15960
                const double arg = (zr+2*v+1)*atan(2.0*xrz);  //Z=15961
                const double nen = pow(1.0+4*xrz*xrz,(zr+2*v+1)/2.0);  //Z=15962
                const double arg1 = (zr+2*v)*atan(2.0*xrz);  //Z=15963
                const double nen1 = pow(1.0+4*xrz*xrz,(zr+2*v)/2.0);  //Z=15964
                const double arg2 = (zr+2*v-1)*atan(2.0*xrz);  //Z=15965
                const double nen2 = pow(1.0+4*xrz*xrz,(zr+2*v-1)/2.0);  //Z=15966

                double F12asz = 0.0;  //Z=15968
                for ( int ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=15969
                    const double a1m = params.CR->carr5p[ii]*pow(params.CR->carr7p[ii],v);   /*  carr7p[ii]:=pp[ii];  //Z=15970 */
                    for ( int jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=15971
                        const double a2m = params.CR->carr5p[jj]*pow(params.CR->carr7p[jj],v);  //Z=15972
                        const double xijm = (params.CR->carr3p[ii]-params.CR->carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii];  //Z=15973 */
                        const double arglmz = (zr+1)*atan(xijm);  //Z=15974
                        const double nenlmz = pow(1.0+xijm*xijm,(zr+1)/2.0);  //Z=15975
                        const double xijp = (params.CR->carr3p[ii]+params.CR->carr3p[jj])*q/(zr+1);  //Z=15976
                        const double arglpz = (zr+1)*atan(xijp);  //Z=15977
                        const double nenlpz = pow(1.0+xijp*xijp,(zr+1)/2.0);  //Z=15978
                        const double F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(M_PI*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(M_PI*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=15979
                        const double F12as2z = e0*e1*(1/(params.CR->carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=15980
                        const double F12as3z = e1*e0*(1/(params.CR->carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=15981
                        const double F12as4z = e1*e1*(1/(params.CR->carr7p[ii]*params.CR->carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(M_PI*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(M_PI*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=15982

                        F12asz += a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);  //Z=15984
                    }/*6*/  //Z=15985
                }/*5*/  //Z=15986
                return /*F12 =*/ preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;  //Z=15987
                //F12 = F12asy;  //Z=15988
            }/*4*/  //Z=15989
            // /*formpq:=*/ return F12;  //Z=15990

            /* formpq:=polyliposome(llipt,radius,lliph,lin,lout,nom,sigmar,sigmal,phiax,philiph,philipt,phiin,phiout,2,q);  //Z=15992 */
            /* formpq:=polyliposome(2.0,200,1.0,3.5,3.5,1,sigmar,sigmal,0.001,-0.55,-0.7,0.001,0.001,3,q);  //Z=15993 */
            /* formpq:=pql;  //Z=15994 */
        }/*3*/ /*  of myelin sphere  */  //Z=15995

    //}/*2*/ /*  of sphere  */  //Z=15997

    return 0.0; // to avoid compiler warning
}/*1*/  //Z=18387


#endif // SC_LIB_FORMPQ_partSphere_H
