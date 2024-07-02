#ifndef SC_LIB_FORMFQ_H
#define SC_LIB_FORMFQ_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formfq( double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const
{  //Z=18395

    // Tests mit LOAD PARAMS C:\SimLab\sas-crystal\Manuskript HyperSeries (Prof.Förster)\SimulTest(Disks-gaus).ini
    // Gespräch mit Prof. Förster (05.Jun.2023): Es sollten immer die mit f hinten genutzt werden...

    double pqsum, oldpqsum, binsum, delser, argq, arglq, pqr, pqr1, pqr2, pql;  //Z=18401
    double ccc1, ccc2, ccc3, vv3;  //Z=18402
    double cc1, cc4, cc6;   //Z=18403
    double F121, F122, F123;  //Z=18406
    double dim, xrad, xradp, x1z, x12z, x2z, x22z, lim, lim1, lim4, lim6;  //Z=18408
    double a1, b1, b2, b1s, v, c, e0, e1, ee0, ee1;  //Z=18409
    double gb1s, pz2v, pz2v1, gz1, preg1, preg3, preg4;  //Z=18410
    double pzc, pzc1, pza, sumc;  //Z=18411
    double del, delc, F12, F12sez, oldF12sez, F42, F42sez, oldF42sez, F62, F62sez, oldF62sez;  //Z=18412
    double arg11, nen11, arg12, nen12, F12as1z, F12as2z, F12asz;  //Z=18413
    double arg44, nen44, arg45, nen45, F42as10z, F42as1z0, F42as40z, F42as27, F42as28, F42as4z, F42asz0;  //Z=18414
    double arg64, nen64, arg65, nen65, F62as10z, F62as1z0, F62as40z, F62as27, F62as28, F62as4z, F62asz0, FF1;  //Z=18415

    const double zl = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=18418
    const double zr = (1-sqr(params.sigma))/sqr(params.sigma);  //Z=18419
    const double radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=18420

    // TODO Unbekannte Variablen:
    double argpq, zz, pqr3, qnarg, binsum1;
    double qxn[121], qyn[121];
    double qz=1.0; // TODO: in qrombdeltac Aufruf verwendet, siehe auch bei formpq()
    double qxhklt=0,qyhklt=0,qzhklt=0,qhkl=0;

    CHECKENDTHREAD_VAL

    /* ************ */  //Z=18422
    /* ** sphere ** */  //Z=18423
    /* ************ */  //Z=18424
    if ( params.part==0 )
    {/*2*/  //Z=18425
        // TODO: hier sollte das Programm nie mehr hinkommen
        return 0.0;

#ifdef undef
        /* ** homogeneous sphere ** */  //Z=18426
        if ( params.cs==0 )
        {/*3*/  //Z=18427
            if ( q<(0.4*params.limq4f) )
            {/*4*/  //Z=18428
                double pqsum = 1.0;  //Z=18429
                double oldpqsum = 0.0;  //Z=18430
                double qqnn = 1.0;  //Z=18431
                for ( int nser=1; nser<=100; nser++ )
                {/*5*/  //Z=18432
                    qqnn = qqnn*q*q;  //Z=18433
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=18434
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18435
                    if ( delser<0.0001 ) break; /* goto 50; */  //Z=18436
                    oldpqsum = pqsum;  //Z=18437
                }/*5*/  //Z=18438
                /*50:*/  //Z=18439
                /*formfq:=*/ return pqsum;  //Z=18440
            }/*4*/  //Z=18441
            else
            {/*4*/  //Z=18442
                const double argq = q*params.radius/(zr+1);  //Z=18443
                const double pqr = (1/(zr*(zr-1))*pow(argq,-2));  //Z=18444
                const double pq1 = pqr*cos((zr-1)*atan(argq))/pow(1.0+argq*argq,(zr-1)/2.0);  //Z=18445
                const double pq2 = (pqr/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1.0+argq*argq,(zr-2)/2.0);  //Z=18446
                const double pq3 = 3*(pq2-pq1);  //Z=18447
                /*formfq:=*/ return pq3*pq3;  //Z=18448
            }/*4*/  //Z=18449
        }/*3*/ /*  of homogeneous sphere */  //Z=18450

        /* ** core/shell sphere ** */  //Z=18452
        if ( params.cs==1 )
        {/*3*/  //Z=18453
            const double ccc1 = sqr(1-params.rho)*pow(params.p1,6);  //Z=18454
            const double ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,3);  //Z=18455
            const double ccc3 = sqr(params.rho);  //Z=18456
            const double vv3 = sqr((1-params.rho)*pow(params.p1,3)+params.rho);  //Z=18457

            const double argq = q*radiusm/(zr+1);  //Z=18459
            const double argpq = q*params.radius/(zr+1);  //Z=18460

            double F121, F122, F123;    // unten wird nur F123 und nicht die F121+F122+F123 genutzt (TODO)

            /*  F121 sphere  */  //Z=18462
            if ( q<(params.limq4f/2.0) )
            {/*4*/  //Z=18463
                double qqnn = 1.0;  //Z=18464
                double pqsum = 1.0;  //Z=18465
                double oldpqsum = 0.0;  //Z=18466
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18467
                    qqnn = qqnn*q*q;  //Z=18468
                    pqsum = pqsum+qqnn*params.CR->carr4f[nser];  //Z=18469
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18470
                    if ( delser<0.0001 ) break; /* goto 51; */  //Z=18471
                    oldpqsum = pqsum;  //Z=18472
                }/*5*/  //Z=18473
                /*51:*/  //Z=18474
                F121 = ccc1*pqsum/vv3;  //Z=18475
            }/*4*/  //Z=18476
            else
            {/*4*/  //Z=18477
                const double pqr = (1/(zr*(zr-1))*pow(argpq,-2));  //Z=18478
                const double pq1 = pqr*cos((zr-1)*atan(argpq))/pow(1.0+argpq*argpq,(zr-1)/2.0);  //Z=18479
                const double pq2 = (pqr/((zr-2)*argpq))*sin((zr-2)*atan(argpq))/pow(1.0+argpq*argpq,(zr-2)/2.0);  //Z=18480
                const double pq3 = 3*(pq2-pq1);  //Z=18481
                F121 = ccc1*pq3*pq3/vv3;  //Z=18482
            }/*4*/  //Z=18483

            /*  F122 sphere  */  //Z=18485
            if ( q<(0.3*params.limq5f) )
            {/*4*/  //Z=18486
                double qqnn = 1.0;  //Z=18487
                double pqsum = 1.0;  //Z=18488
                double oldpqsum = 0.0;  //Z=18489
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18490
                    qqnn = qqnn*q*q;  //Z=18491
                    pqsum = pqsum+qqnn*params.CR->carr5f[nser];  //Z=18492
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18493
                    if ( delser<0.0001 ) break; /* goto 52; */  //Z=18494
                    oldpqsum = pqsum;  //Z=18495
                }/*5*/  //Z=18496
                /*52:*/  //Z=18497
                F122 = ccc2*pqsum/vv3;  //Z=18498
            }/*4*/  //Z=18499
            else
            {/*4*/  //Z=18500
                const double pqr1 = (1/(zr*(zr-1))*pow(argpq,-2));  //Z=18501
                const double pq1 = pqr1*cos((zr-1)*atan(argpq))/pow(1.0+argpq*argpq,(zr-1)/2.0);  //Z=18502
                const double pq2 = (pqr1/((zr-2)*argpq))*sin((zr-2)*atan(argpq))/pow(1.0+argpq*argpq,(zr-2)/2.0);  //Z=18503
                const double pq3 = 3*(pq2-pq1);  //Z=18504

                const double pqr2 = (1/(zr*(zr-1))*pow(argq,-2));  //Z=18506
                const double pq4 = pqr2*cos((zr-1)*atan(argq))/pow(1.0+argq*argq,(zr-1)/2.0);  //Z=18507
                const double pq5 = (pqr2/((zr-2)*argq))*sin((zr-2)*atan(argq))/pow(1.0+argq*argq,(zr-2)/2.0);  //Z=18508
                const double pq6 = 3*(pq5-pq4);  //Z=18509

                F122 = ccc2*pq3*pq6/vv3;  //Z=18511
            }/*4*/  //Z=18512

            /*  F123 sphere  */  //Z=18514
            if ( q<(params.limq6f/2.0) )
            {/*4*/  //Z=18515
                double qqnn = 1.0;  //Z=18516
                double pqsum = 1.0;  //Z=18517
                double oldpqsum = 0.0;  //Z=18518
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18519
                    qqnn = qqnn*q*q;  //Z=18520
                    pqsum = pqsum+qqnn*params.CR->carr6f[nser];  //Z=18521
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18522
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
            /* formfq:=F121+F122+F123;  //Z=18536 */   // TODO: welche Zeile stimmt hier?
            /*formfq:=*/ return F123;  //Z=18537
        }/*3*/  /*  of core/shell sphere  */  //Z=18538

        /* ** inhomogeneous core/shell sphere ** */  //Z=18540
        if ( params.cs==2 )
        {/*3*/  //Z=18541

            const double dim = 3;  //Z=18543
            const double delc = 0.0001;  //Z=18544
            const double zz = zr;  //Z=18545
            const double xrad = q*radiusm;  //Z=18546
            const double xradp = q*params.radius;  //Z=18547
            const double x1z = q*params.radius/(2.0*(zz+1));  //Z=18548
            const double x12z = x1z*x1z;  //Z=18549
            const double x2z = q*radiusm/(2.0*(zz+1));  //Z=18550
            const double x22z = x2z*x2z;  //Z=18551

            const double lim = 18*exp(-5*params.sigma);  //Z=18553
            const double lim1 = lim;  //Z=18554
            //lim2 = lim*0.7;  //Z=18555
            //lim3 = lim;  //Z=18556
            const double lim4 = lim;  //Z=18557
            //lim5 = lim*0.7;  //Z=18558
            const double lim6 = lim*1.2;  //Z=18559

            const double a1 = (dim-params.alphash1)/2.0;  //Z=18561
            const double b1 = dim/2.0;  //Z=18562
            const double b2 = (dim+2-params.alphash1)/2.0;  //Z=18563
            const double b1s = (dim+2)/2.0;  //Z=18564
            const double v = -b1s+1/2.0;  //Z=18565
            const double c = a1-b1-b2+1/2.0;  //Z=18566
            //d0 = 1;  //Z=18567
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=18568
            const double e0 = 1.0;  //Z=18569
            const double e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=18570
            const double ee0 = 1.0;  //Z=18571
            const double ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=18572

            const double gb1s = 3*sqrt(M_PI)/4.0;  //Z=18574
            const double pz2v = 1/(zr*(zr-1));  //Z=18575
            const double pz2v1 = pz2v/(zr-2);  //Z=18576
            //pz2v2 = pz2v1/(zr-3);  //Z=18577

            const double gz1 = gamma(zr+1);  //Z=18579
            const double preg1 = gb1s/sqrt(M_PI);  //Z=18580
            const double preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=18581
            const double preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=18582
            //pzvc = gamma(zr+1+v+c)/gz1;  //Z=18583
            //pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=18584
            //pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=18585
            //pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=18586
            //pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=18587
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=18588
            const double pzc = gamma(zr+1+c)/gz1;  //Z=18589
            const double pzc1 = gamma(zr+1+c-1)/gz1;  //Z=18590
            const double pza = gamma(zr+1-2*a1)/gz1;  //Z=18591
            //pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=18592
            //pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=18593
            //dnv0 = 1;  //Z=18594
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=18595
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=18596
            //pva0 = gamma(zr+1-2*a1)/gz1;  //Z=18597

            const double cc1 = 1/dim;  //Z=18599
            const double cc4 = params.rho/((dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=18600
            const double cc6 = -params.rho/(dim-params.alphash1);  //Z=18601
            const double sumc = cc1+cc4+cc6;  //Z=18602

            double F12, F42, F62;

            /*  term #1 series  */  //Z=18604
            if ( (xradp)<lim1 )
            {/*4*/  //Z=18605
                //z12v[0] = 1;  //Z=18606
                //b1sv[0] = 1;  //Z=18607
                //fkv[0] = 1;  //Z=18608
                double qqnn = 1.0;  //Z=18609
                F12 = 1.0;  //Z=18610
                double oldF12 = 1.0;  //Z=18611
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=18612
                    qqnn = qqnn*q*q;  //Z=18613 war: qnn[..], ist zwar definiert, aber oben wurde qqn[0]=1 gesetzt... (TODO)
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=18614
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=18615
                    //fkv[n] = fkv[n-1]*n;  //Z=18616
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]);  //Z=18617 */

                    F12 += params.CR->carr4f[n]*qqnn;  //Z=18619

                    double del = fabs((F12-oldF12)/F12);  //Z=18621
                    if ( del<delc ) break; /* goto 101; */  //Z=18622
                    oldF12 = F12;  //Z=18623
                }/*5*/  //Z=18624
                /*101:*/  //Z=18625
                //F12 = F12sez;  //Z=18626
            }/*4*/  //Z=18627

            /*  term #4 series  */  //Z=18629
            if ( xradp<lim4 )
            {/*4*/  //Z=18630
                //z12v[0] = 1;  //Z=18631
                //a1v[0] = 1;  //Z=18632
                //b1v[0] = 1;  //Z=18633
                //b2v[0] = 1;  //Z=18634
                //b1sv[0] = 1;  //Z=18635
                //fkv[0] = 1;  //Z=18636
                double qqnn = 1.0;  //Z=18637
                F42 = 1.0;  //Z=18638
                double oldF42 = 1.0;  //Z=18639
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=18640
                    qqnn = qqnn*q*q;  //Z=18641
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=18642
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=18643
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=18644
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=18645
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=18646
                    //fkv[n] = fkv[n-1]*n;  //Z=18647
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=18648 */

                    F42 += params.CR->carr5f[n]*qqnn;  //Z=18650

                    double del = fabs((F42-oldF42)/F42);  //Z=18652
                    if ( del<delc ) break; /* goto 104; */  //Z=18653
                    oldF42 = F42;  //Z=18654
                }/*5*/  //Z=18655
                /*104:*/  //Z=18656
                //F42 = F42sez;  //Z=18657
            }/*4*/  //Z=18658

            /*  term #6 series  */  //Z=18660
            if ( xradp<lim6 )
            {/*4*/  //Z=18661
                //z12v[0] = 1;  //Z=18662
                //a1v[0] = 1;  //Z=18663
                //b1v[0] = 1;  //Z=18664
                //b2v[0] = 1;  //Z=18665
                //b1sv[0] = 1;  //Z=18666
                //fkv[0] = 1;  //Z=18667
                double qqnn = 1.0;  //Z=18668
                F62 = 1.0;  //Z=18669
                double oldF62 = 1.0;  //Z=18670
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=18671
                    qqnn = qqnn*q*q;  //Z=18672
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=18673
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=18674
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=18675
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=18676
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=18677
                    //fkv[n] = fkv[n-1]*n;  //Z=18678
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=18679 */

                    F62 += F62+params.CR->carr6f[n]*qqnn;  //Z=18681

                    double del = fabs((F62-oldF62)/F62);  //Z=18683
                    if ( del<delc ) break; /* goto 106; */  //Z=18684
                    oldF62 = F62;  //Z=18685
                }/*5*/  //Z=18686
                /*106:*/  //Z=18687
                //F62 = F62sez;  //Z=18688
            }/*4*/  //Z=18689

            /* ** term #1 asymptote ** */  //Z=18691
            if ( xradp>=lim1 )
            {/*4*/  //Z=18692
                const double arg11 = (zr+v+1)*atan(2.0*x1z);  //Z=18693
                const double nen11 = pow(1.0+4*x1z*x1z,(zr+v+1)/2.0);  //Z=18694
                const double arg12 = (zr+v)*atan(2.0*x1z);  //Z=18695
                const double nen12 = pow(1.0+4*x1z*x1z,(zr+v)/2.0);  //Z=18696

                const double F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);  //Z=18698
                const double F12as2z = ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);  //Z=18699
                F12 = preg1*pow(x1z,v)*(F12as1z+F12as2z);  //Z=18700
                //F12 = F12asz;  //Z=18701
            }/*4*/  //Z=18702

            /* ** term #4 asymptote ** */  //Z=18704
            if ( xrad>=lim4 )
            {/*4*/  //Z=18705
                const double F42as10z = preg4*pow(x22z,-a1);  //Z=18706
                //F42as1sumz = pva0;  //Z=18707
                //F42as1z = F42as10z*F42as1sumz;  //Z=18708
                const double F42as1z0 = F42as10z*pza;   /* * */  //Z=18709

                const double F42as40z = preg3*pow(x2z,c);  //Z=18711
                const double arg44 = (zr+c+1)*atan(2.0*x2z);  //Z=18712
                const double nen44 = pow(1.0+4*x2z*x2z,(zr+c+1)/2.0);  //Z=18713
                const double arg45 = (zr+c)*atan(2.0*x2z);  //Z=18714
                const double nen45 = pow(1.0+4*x2z*x2z,(zr+c)/2.0);  //Z=18715
                const double F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);  //Z=18716
                const double F42as28 = e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);  //Z=18717
                const double F42as4z = F42as40z*(F42as27+F42as28);  //Z=18718
                //F42asz = F42as1z+F42as4z;  //Z=18719
                F42 = F42as1z0+F42as4z;  //Z=18720
                //F42 = F42asz0;  //Z=18721
            }/*4*/  //Z=18722

            /* ** term #6 asymptote ** */  //Z=18724
            if ( xradp>=lim6 )
            {/*4*/  //Z=18725
                const double F62as10z = preg4*pow(x12z,-a1);  //Z=18726
                //F62as1sumz = pva0;  //Z=18727
                //F62as1z = F62as10z*F62as1sumz;  //Z=18728
                const double F62as1z0 = F62as10z*pza;     /* * */  //Z=18729

                const double F62as40z = preg3*pow(x1z,c);  //Z=18731
                const double arg64 = (zr+c+1)*atan(2.0*x1z);  //Z=18732
                const double nen64 = pow(1.0+4*x1z*x1z,(zr+c+1)/2.0);  //Z=18733
                const double arg65 = (zr+c)*atan(2.0*x1z);  //Z=18734
                const double nen65 = pow(1.0+4*x1z*x1z,(zr+c)/2.0);  //Z=18735
                const double F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);  //Z=18736
                const double F62as28 = e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);  //Z=18737
                const double F62as4z = F62as40z*(F62as27+F62as28);  //Z=18738
                //F62asz = F62as1z+F62as4z;  //Z=18739
                F62 = F62as1z0+F62as4z;  //Z=18740
                //F62 = F62asz0;  //Z=18741
            }/*4*/  //Z=18742

            /* FF1:=(cc1*F12+cc4*F42+cc6*F62)/sumc;  //Z=18744 */
            FF1 = (cc1*F12)/sumc;  //Z=18745

            /*formfq:=*/ return FF1*FF1;  //Z=18747


            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,3,sigmar,q);  //Z=18750 */


        }/*3*/  /*  of inhomogeneous core/shell sphere  */  //Z=18753
#endif
    }/*2*/ /*  of sphere  */  //Z=18755


    /* ********** */  //Z=18758
    /*  cylinder  */  //Z=18759
    /* ********** */  //Z=18760
    if ( params.part==1 )
    {/*2*/  //Z=18761

        double pql, pqr;

        /* ** longitudinal part ** */  //Z=18763
        /* ** isotropic ** */  //Z=18764
        if ( ordis==7 )
        {/*3*/  //Z=18765
            if ( q<(0.6*params.limq1f) )
            {/*4*/  //Z=18766
                double pqsum = 1.0;  //Z=18767
                double oldpqsum = 0.0;  //Z=18768
                double qqnn = 1.0;  //Z=18769
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18770
                    qqnn = qqnn*q*q;  //Z=18771
                    pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=18772
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18773
                    if ( delser<0.0001 ) break; /* goto 60; */  //Z=18774
                    oldpqsum = pqsum;  //Z=18775
                }/*5*/  //Z=18776
                /*60:*/  //Z=18777
                pql = pqsum;  //Z=18778
            }/*4*/  //Z=18779
            else
            {/*4*/   /*  = P(q)  */  //Z=18780
                const double arglq = q*params.length/(zl+1);  //Z=18781
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
                    double pqsum = 1.0;  //Z=18793
                    double oldpqsum = 0.0;  //Z=18794
                    double qqnn = 1.0;  //Z=18795
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18796
                        qqnn = qqnn*(qxs+qys)*(qxs+qys);  //Z=18797
                        pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=18798
                        double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18799
                        if ( delser<0.0001 ) break; /* goto 65; */  //Z=18800
                        oldpqsum = pqsum;  //Z=18801
                    }/*6*/  //Z=18802
                    /*65:*/  //Z=18803
                    pql = pqsum;  //Z=18804
                }/*5*/  //Z=18805
                else
                {/*5*/  //Z=18806
                    const double arglq = (qxs+qys+eps9)*params.length/(zl+1);  //Z=18807
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
                    double pqsum = 1.0;  //Z=18822
                    double oldpqsum = 0.0;  //Z=18823
                    double qxn[121], qyn[121];
                    qxn[0] = 1.0;  //Z=18824
                    qyn[0] = 1.0;  //Z=18825
                    if ( params.orcase==1 )
                    {/*6*/  //Z=18826
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18827
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18828
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18829
                            double binsum = 0.0;  //Z=18830
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18831
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18832 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18833 */
                                binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qyn[mser]*qxn[nser-mser];  //Z=18834
                            }/*8*/  //Z=18835
                            pqsum = pqsum+params.CR->carr1f[nser]*binsum;  //Z=18836
                            double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18837
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=18838
                            oldpqsum = pqsum;  //Z=18839
                        }/*7*/  //Z=18840
                    }/*6*/  //Z=18841

                    if ( params.orcase==2 )
                    {/*6*/  /*  x-axis  */  //Z=18843
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18844
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18845
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18846
                            double binsum = 0.0;  //Z=18847
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18848
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18849 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18850 */
                                binsum = binsum+params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];  //Z=18851
                            }/*8*/  //Z=18852
                            pqsum = pqsum+params.CR->carr1f[nser]*binsum;  //Z=18853
                            double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18854
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=18855
                            oldpqsum = pqsum;  //Z=18856
                        }/*7*/  //Z=18857
                    }/*6*/  //Z=18858

                    if ( params.orcase==3 )
                    {/*6*/  /*  y-axis  */  //Z=18861
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18862
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18863
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18864
                            double binsum = 0.0;  //Z=18865
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18866
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18867 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18868 */
                                binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[mser]*qyn[nser-mser];  //Z=18869
                            }/*8*/  //Z=18870
                            pqsum = pqsum+params.CR->carr1f[nser]*binsum;  //Z=18871
                            double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18872
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
                double pqsum = 1.0;  //Z=18894
                double oldpqsum = 0.0;  //Z=18895
                double qqnn = 1.0;  //Z=18896
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18897
                    qqnn = qqnn*q*q;  //Z=18898
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=18899
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18900
                    if ( delser<0.0001 ) break; /* goto 61; */  //Z=18901
                    oldpqsum = pqsum;  //Z=18902
                }/*5*/  //Z=18903
                /*61:*/  //Z=18904
                pqr = pqsum;  //Z=18905
            }/*4*/  //Z=18906
            else
            {/*4*/  //Z=18907
                const double argpq = q*params.radius/(zr+1);  //Z=18908
                const double pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-1/2.0)/2.0);  //Z=18909
                const double pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-3/2.0)/2.0);  //Z=18910
                const double pqr3 = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);  //Z=18911
                pqr = pqr3*pqr3;  //Z=18912
            }/*4*/  //Z=18913
            /*formfq:=*/ return pql*pqr;  //Z=18914
            /* formfq:=pql;  //Z=18915 */
        }/*3*/ /*  of homogeneous cylinder  */  //Z=18916

        /*  homogeneous core/shell cylinder  */  //Z=18918
        if ( params.cs==1 )
        {/*3*/  //Z=18919
            const double ccc1 = sqr(1-params.rho)*pow(params.p1,4);  //Z=18920
            const double ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,2);  //Z=18921
            const double ccc3 = sqr(params.rho);  //Z=18922
            const double vv3 = sqr((1-params.rho)*pow(params.p1,2)+params.rho);  //Z=18923

            const double zz = zr;  // TODO: zz war in diesem Zweig nicht gesetzt
            const double argq = q*radiusm/(zz+1);  //Z=18925
            const double argpq = q*params.radius/(zz+1);  //Z=18926

            double F121, F122, F123;

            /*  F121 cylinder  */  //Z=18928
            if ( q<(0.7*params.limq4f) )
            {/*4*/  //Z=18929
                /* ** series expansion ** */  //Z=18930
                double pqsum = 1.0;  //Z=18931
                double oldpqsum = 0.0;  //Z=18932
                double qqnn = 1.0;  //Z=18933
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18934
                    qqnn = qqnn*q*q;  //Z=18935
                    pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=18936
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18937
                    if ( delser<0.0001 ) break; /* goto 62; */  //Z=18938
                    oldpqsum = pqsum;  //Z=18939
                }/*5*/  //Z=18940
                /*62:*/  //Z=18941
                F121 = ccc1*pqsum/vv3;  //Z=18942
            }/*4*/  //Z=18943
            else
            {/*4*/  //Z=18944
                const double pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-1/2.0)/2.0);  //Z=18945
                const double pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1.0+argpq*argpq,(zr-3/2.0)/2.0);  //Z=18946
                const double pqr = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);  //Z=18947
                //pqr = pqr3;  //Z=18948
                F121 = ccc1*pqr*pqr/vv3;  //Z=18949
            }/*4*/  //Z=18950

            /*  F122 cylinder  */  //Z=18952
            if ( q<(0.3*params.limq5f) )
            {/*4*/  //Z=18953
                /* ** series expansion ** */  //Z=18954
                double pqsum = 1.0;  //Z=18955
                double oldpqsum = 0.0;  //Z=18956
                double qqnn = 1.0;  //Z=18957
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18958
                    qqnn = qqnn*q*q;  //Z=18959
                    pqsum = pqsum+params.CR->carr5f[nser]*qqnn;  //Z=18960
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18961
                    if ( delser<0.0001 ) break; /* goto 63; */  //Z=18962
                    oldpqsum = pqsum;  //Z=18963
                }/*5*/  //Z=18964
                /*63:*/  //Z=18965
                F122 = ccc2*pqsum/vv3;  //Z=18966
            }/*4*/  //Z=18967
            else
            {/*4*/  //Z=18968
                const double pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argpq,-3/2.0)*(sin((zr-1/2.0)*atan(argpq))-cos((zr-1/2.0)*atan(argpq)))/pow(1.0+argpq*argq,(zr-1/2.0)/2.0);  //Z=18969
                const double pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argpq,-5/2.0)*(sin((zr-3/2.0)*atan(argpq))+cos((zr-3/2.0)*atan(argpq)))/pow(1.0+argpq*argq,(zr-3/2.0)/2.0);  //Z=18970
                const double pqr3 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argq,-3/2.0)*(sin((zr-1/2.0)*atan(argq))-cos((zr-1/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-1/2.0)/2.0);  //Z=18971
                const double pqr4 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argq,-5/2.0)*(sin((zr-3/2.0)*atan(argq))+cos((zr-3/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-3/2.0)/2.0);  //Z=18972
                const double pqr = (4/M_PI)*(pqr1+(9/16.0)*pqr2)*(pqr3+(9/16.0)*pqr4);  //Z=18973
                F122 = ccc2*pqr/vv3;  //Z=18974
            }/*4*/  //Z=18975

            /*  F123 cylinder  */  //Z=18977
            if ( q<(0.6*params.limq6f) )
            {/*4*/  //Z=18978
                /* ** series expansion ** */  //Z=18979
                double pqsum = 1.0;  //Z=18980
                double oldpqsum = 0.0;  //Z=18981
                double qqnn = 1.0;  //Z=18982
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=18983
                    qqnn = qqnn*q*q;  //Z=18984
                    pqsum = pqsum+params.CR->carr6f[nser]*qqnn;  //Z=18985
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18986
                    if ( delser<0.0001 ) break; /* goto 64; */  //Z=18987
                    oldpqsum = pqsum;  //Z=18988
                }/*5*/  //Z=18989
                /*64:*/  //Z=18990
                F123 = ccc3*pqsum/vv3;  //Z=18991
            }/*4*/  //Z=18992
            else
            {/*4*/  //Z=18993
                const double pqr1 = (gamma(zr-1/2.0)/gamma(zr+1))*pow(argq,-3/2.0)*(sin((zr-1/2.0)*atan(argq))-cos((zr-1/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-1/2.0)/2.0);  //Z=18994
                const double pqr2 = (gamma(zr-3/2.0)/gamma(zr+1))*pow(argq,-5/2.0)*(sin((zr-3/2.0)*atan(argq))+cos((zr-3/2.0)*atan(argq)))/pow(1.0+argq*argq,(zr-3/2.0)/2.0);  //Z=18995
                const double pqr = (2/sqrt(M_PI))*(pqr1+(9/16.0)*pqr2);  //Z=18996
                //pqr = pqr3;  //Z=18997
                F123 = ccc3*pqr*pqr/vv3;  //Z=18998
            }/*4*/  //Z=18999
            /*formfq:=*/ return pql*(F121+F122+F123);  //Z=19000
            /* formfq:=(F121+F122+F123);  //Z=19001 */
        }/*3*/ /*  of homogeneous core/shell cylinder  */  //Z=19002

        /* ** inhomogeneous core/shell cylinder ** */  //Z=19004
        if ( params.cs==2 )
        {/*3*/  //Z=19005

            const double dim = 2;  //Z=19007
            const double delc = 0.0001;  //Z=19008
            const double zz = zr;  //Z=19009
            const double xrad = q*radiusm;  //Z=19010
            const double xradp = q*params.radius;  //Z=19011
            const double x1z = q*params.radius/(2.0*(zz+1));  //Z=19012
            const double x12z = x1z*x1z;  //Z=19013
            const double x2z = q*radiusm/(2.0*(zz+1));  //Z=19014
            const double x22z = x2z*x2z;  //Z=19015

            const double lim = 18*exp(-5*params.sigma);  //Z=19017
            const double lim1 = lim;  //Z=19018
            //lim2 = lim*0.7;  //Z=19019
            //lim3 = lim;  //Z=19020
            const double lim4 = lim;  //Z=19021
            //lim5 = lim*0.7;  //Z=19022
            const double lim6 = lim*1.2;  //Z=19023

            const double a1 = (dim-params.alphash1)/2.0;  //Z=19025
            const double b1 = dim/2.0;  //Z=19026
            const double b2 = (dim+2-params.alphash1)/2.0;  //Z=19027
            const double b1s = (dim+2)/2.0;  //Z=19028
            const double v = -b1s+1/2.0;  //Z=19029
            const double c = a1-b1-b2+1/2.0;  //Z=19030
            //d0 = 1;  //Z=19031
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=19032
            const double e0 = 1.0;  //Z=19033
            const double e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=19034
            const double ee0 = 1.0;  //Z=19035
            const double ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=19036

            const double gb1s = 3*sqrt(M_PI)/4.0;  //Z=19038
            const double pz2v = 1/(zr*(zr-1));  //Z=19039
            const double pz2v1 = pz2v/(zr-2);  //Z=19040
            //pz2v2 = pz2v1/(zr-3);  //Z=19041

            const double gz1 = gamma(zr+1);  //Z=19043
            const double preg1 = gb1s/sqrt(M_PI);  //Z=19044
            const double preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=19045
            const double preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=19046
            //pzvc = gamma(zr+1+v+c)/gz1;  //Z=19047
            //pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=19048
            //pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=19049
            //pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=19050
            //pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=19051
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=19052
            const double pzc = gamma(zr+1+c)/gz1;  //Z=19053
            const double pzc1 = gamma(zr+1+c-1)/gz1;  //Z=19054
            const double pza = gamma(zr+1-2*a1)/gz1;  //Z=19055
            //pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=19056
            //pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19057
            //dnv0 = 1;  //Z=19058
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=19059
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=19060
            //pva0 = gamma(zr+1-2*a1)/gz1;  //Z=19061

            const double cc1 = 1/dim;  //Z=19063
            const double cc4 = params.rho/((dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=19064
            const double cc6 = -params.rho/(dim-params.alphash1);  //Z=19065
            const double sumc = cc1+cc4+cc6;  //Z=19066

            double F12, F42, F62;

            /*  term #1 series  */  //Z=19068
            if ( (xradp)<lim1 )
            {/*4*/  //Z=19069
                //z12v[0] = 1;  //Z=19070
                //b1sv[0] = 1;  //Z=19071
                //fkv[0] = 1;  //Z=19072
                double qqnn = 1.0;  //Z=19073
                F12 = 1.0;  //Z=19074
                double oldF12 = 1.0;  //Z=19075
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=19076
                    qqnn = qqnn*q*q;  //Z=19077
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19078
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19079
                    //fkv[n] = fkv[n-1]*n;  //Z=19080
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]/(b1sv[n]*fkv[n]);  //Z=19081 */

                    F12 += params.CR->carr4f[n]*qqnn;  //Z=19083

                    double del = fabs((F12-oldF12)/F12);  //Z=19085
                    if ( del<delc ) break; /* goto 111; */  //Z=19086
                    oldF12 = F12;  //Z=19087
                }/*5*/  //Z=19088
                /*111:*/  //Z=19089
                //F12 = F12sez;  //Z=19090
            }/*4*/  //Z=19091

            /*  term #4 series  */  //Z=19093
            if ( xradp<lim4 )
            {/*4*/  //Z=19094
                //z12v[0] = 1;  //Z=19095
                //a1v[0] = 1;  //Z=19096
                //b1v[0] = 1;  //Z=19097
                //b2v[0] = 1;  //Z=19098
                //b1sv[0] = 1;  //Z=19099
                //fkv[0] = 1;  //Z=19100
                double qqnn = 1.0;  //Z=19101
                F42 = 1.0;  //Z=19102
                double oldF42 = 1.0;  //Z=19103
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=19104
                    qqnn = qqnn*q*q;  //Z=19105
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19106
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19107
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19108
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19109
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19110
                    //fkv[n] = fkv[n-1]*n;  //Z=19111
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19112 */

                    F42 += params.CR->carr5f[n]*qqnn;  //Z=19114

                    del = fabs((F42-oldF42)/F42);  //Z=19116
                    if ( del<delc ) break; /* goto 114; */  //Z=19117
                    oldF42 = F42;  //Z=19118
                }/*5*/  //Z=19119
                /*114:*/  //Z=19120
                //F42 = F42sez;  //Z=19121
            }/*4*/  //Z=19122

            /*  term #6 series  */  //Z=19124
            if ( xradp<lim6 )
            {/*4*/  //Z=19125
                //z12v[0] = 1;  //Z=19126
                //a1v[0] = 1;  //Z=19127
                //b1v[0] = 1;  //Z=19128
                //b2v[0] = 1;  //Z=19129
                //b1sv[0] = 1;  //Z=19130
                //fkv[0] = 1;  //Z=19131
                double qqnn = 1.0;  //Z=19132
                F62 = 1.0;  //Z=19133
                double oldF62 = 1.0;  //Z=19134
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=19135
                    qqnn = qqnn*q*q;  //Z=19136
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19137
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=19138
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=19139
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=19140
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19141
                    //fkv[n] = fkv[n-1]*n;  //Z=19142
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*a1v[n]/(b1v[n]*b2v[n]*fkv[n]);  //Z=19143 */

                    F62 += params.CR->carr6f[n]*qqnn;  //Z=19145

                    double del = fabs((F62-oldF62)/F62);  //Z=19147
                    if ( del<delc ) break; /* goto 116; */  //Z=19148
                    oldF62 = F62;  //Z=19149
                }/*5*/  //Z=19150
                /*116:*/  //Z=19151
                //F62 = F62sez;  //Z=19152
            }/*4*/  //Z=19153

            /* ** term #1 asymptote ** */  //Z=19155
            if ( xradp>=lim1 )
            {/*4*/  //Z=19156
                const double arg11 = (zr+v+1)*atan(2.0*x1z);  //Z=19157
                const double nen11 = pow(1.0+4*x1z*x1z,(zr+v+1)/2.0);  //Z=19158
                const double arg12 = (zr+v)*atan(2.0*x1z);  //Z=19159
                const double nen12 = pow(1.0+4*x1z*x1z,(zr+v)/2.0);  //Z=19160

                const double F12as1z = ee0*pz2v*(cos(M_PI*v/2.0)*cos(arg11)/nen11-sin(M_PI*v/2.0)*sin(arg11)/nen11);  //Z=19162
                const double F12as2z = ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(v-1)/2.0)*sin(arg12)/nen12);  //Z=19163
                F12 = preg1*pow(x1z,v)*(F12as1z+F12as2z);  //Z=19164
                //F12 = F12asz;  //Z=19165
            }/*4*/  //Z=19166

            /* ** term #4 asymptote ** */  //Z=19168
            if ( xrad>=lim4 )
            {/*4*/  //Z=19169
                const double F42as10z = preg4*pow(x22z,-a1);  //Z=19170
                //F42as1sumz = pva0;  //Z=19171
                //F42as1z = F42as10z*F42as1sumz;  //Z=19172
                const double F42as1z0 = F42as10z*pza;   /* * */  //Z=19173

                const double F42as40z = preg3*pow(x2z,c);  //Z=19175
                const double arg44 = (zr+c+1)*atan(2.0*x2z);  //Z=19176
                const double nen44 = pow(1.0+4*x2z*x2z,(zr+c+1)/2.0);  //Z=19177
                const double arg45 = (zr+c)*atan(2.0*x2z);  //Z=19178
                const double nen45 = pow(1.0+4*x2z*x2z,(zr+c)/2.0);  //Z=19179
                const double F42as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg44)/nen44-sin(M_PI*c/2.0)*sin(arg44)/nen44);  //Z=19180
                const double F42as28 = e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(c-1)/2.0)*sin(arg45)/nen45);  //Z=19181
                const double F42as4z = F42as40z*(F42as27+F42as28);  //Z=19182
                //F42asz = F42as1z+F42as4z;  //Z=19183
                F42 = F42as1z0+F42as4z;  //Z=19184
                //F42 = F42asz0;  //Z=19185
            }/*4*/  //Z=19186

            /* ** term #6 asymptote ** */  //Z=19188
            if ( xradp>=lim6 )
            {/*4*/  //Z=19189
                const double F62as10z = preg4*pow(x12z,-a1);  //Z=19190
                //F62as1sumz = pva0;  //Z=19191
                //F62as1z = F62as10z*F62as1sumz;  //Z=19192
                const double F62as1z0 = F62as10z*pza;     /* * */  //Z=19193

                const double F62as40z = preg3*pow(x1z,c);  //Z=19195
                const double arg64 = (zr+c+1)*atan(2.0*x1z);  //Z=19196
                const double nen64 = pow(1.0+4*x1z*x1z,(zr+c+1)/2.0);  //Z=19197
                const double arg65 = (zr+c)*atan(2.0*x1z);  //Z=19198
                const double nen65 = pow(1.0+4*x1z*x1z,(zr+c)/2.0);  //Z=19199
                const double F62as27 = e0*pzc*(cos(M_PI*c/2.0)*cos(arg64)/nen64-sin(M_PI*c/2.0)*sin(arg64)/nen64);  //Z=19200
                const double F62as28 = e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(c-1)/2.0)*sin(arg65)/nen65);  //Z=19201
                const double F62as4z = F62as40z*(F62as27+F62as28);  //Z=19202
                //F62asz = F62as1z+F62as4z;  //Z=19203
                F62 = F62as1z0+F62as4z;  //Z=19204
                //F62 = F62asz0;  //Z=19205
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
                for ( int nser=1; nser<=80; nser++ )
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
                    for ( int nser=1; nser<=120; nser++ )
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
                    for ( int nser=1; nser<=120; nser++ )
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
                    for ( int nser=1; nser<=120; nser++ )
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
                    for ( int nser=1; nser<=120; nser++ )
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

                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19329
                        qqnn = qqnn*q*q;  //Z=19330
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19331
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19332

                        binsum = 0.0;  //Z=19334
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19335
                            binsum1 = 0.0;  //Z=19336
                            for ( int lser=0; lser<=mser; lser++ )
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

                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19373
                        qqnn = qqnn*q*q;  //Z=19374
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19375
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19376

                        binsum = 0.0;  //Z=19378
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19379
                            binsum1 = 0.0;  //Z=19380
                            for ( int lser=0; lser<=mser; lser++ )
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

                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=19417
                        qqnn = qqnn*q*q;  //Z=19418
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=19419
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=19420

                        binsum = 0.0;  //Z=19422
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=19423
                            binsum1 = 0.0;  //Z=19424
                            for ( int lser=0; lser<=mser; lser++ )
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
                    for ( int nser=1; nser<=120; nser++ )
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
                for ( int nser=1; nser<=120; nser++ )
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
                for ( int nser=1; nser<=120; nser++ )
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
                for ( int nser=1; nser<=120; nser++ )
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
                for ( int nser=1; nser<=120; nser++ )
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
                //z12v[0] = 1;  //Z=19650
                //b1sv[0] = 1;  //Z=19651
                //fkv[0] = 1;  //Z=19652
                double qqnn = 1.0;  //Z=19653
                F12sez = 1.0;  //Z=19654
                oldF12sez = 1.0;  //Z=19655
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=19656
                    qqnn = qqnn*q*q;  //Z=19657
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=19658
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=19659
                    //fkv[n] = fkv[n-1]*n;  //Z=19660
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
                //z12v[0] = 1;  //Z=19675
                //a1v[0] = 1;  //Z=19676
                //b1v[0] = 1;  //Z=19677
                //b2v[0] = 1;  //Z=19678
                //b1sv[0] = 1;  //Z=19679
                //fkv[0] = 1;  //Z=19680
                double qqnn = 1.0;  //Z=19681
                F42sez = 1.0;  //Z=19682
                oldF42sez = 1.0;  //Z=19683
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
                //z12v[0] = 1;  //Z=19706
                //a1v[0] = 1;  //Z=19707
                //b1v[0] = 1;  //Z=19708
                //b2v[0] = 1;  //Z=19709
                //b1sv[0] = 1;  //Z=19710
                //fkv[0] = 1;  //Z=19711
                double qqnn = 1.0;  //Z=19712
                F62sez = 1.0;  //Z=19713
                oldF62sez = 1.0;  //Z=19714
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
                for ( int nser=1; nser<=120; nser++ )
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
                qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,qxhklt,qyhklt,qzhklt,qhkl, /*!9*/
                                                                                                                                                params.ax1.length(),params.ax2.length(),params.ax3.length(),
                            params.ax1.x(),params.ax1.y(),params.ax1.z(),
                            params.ax2.x(),params.ax2.y(),params.ax2.z(),
                            params.ax3.x(),params.ax3.y(),params.ax3.z(),
                            params.sig.x(),params.sig.y(),params.sig.z(),
                            ordis,3,7,12,7,1,0,params.CR->carr1f,pql);  //Z=19818
                /*formfq:=*/ return pql/(M_PI/2.0);  //Z=19819
            }/*4*/  //Z=19820
        }/*3*/ /*  of homogeneous cube */  //Z=19821

#ifdef procnotused
        /*  core/shell cube  */  //Z=19823
        if ( params.cs==1 )
        {/*3*/  //Z=19824
            /*formfq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=19825
        }/*3*/  //Z=19826
#endif

    }/*2*/  /*  of cube  */  //Z=19828

    //qDebug() << "formfq() undef";

    return 0.0;
}  //Z=19829


#endif // SC_LIB_FORMFQ_H
