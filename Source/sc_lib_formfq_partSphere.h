#ifndef SC_LIB_FORMFQ_partSphere_H
#define SC_LIB_FORMFQ_partSphere_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formfq_partSphere( double q ) const
{  //Z=18395

    double pqr, pq1, pq2, pq3;  //Z=18401

    //const double zl = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=18418
    const double zr = (1-sqr(params.sigma))/sqr(params.sigma);  //Z=18419
    const double radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=18420

    CHECKENDTHREAD_VAL

    /* ************ */  //Z=18422
    /* ** sphere ** */  //Z=18423
    /* ************ */  //Z=18424
    //if ( params.part==0 )
    //{/*2*/  //Z=18425
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

//#define UseFormfqF121     (TODO?)
//#define UseFormfqF122
        // in dem Programm vom 01.03.2024 wird nur noch F123 verwendet und nicht mehr
        //  die vorherige Summe F121+F122+F123. Um Compiler-Warnungen zu vermeiden und
        //  die Laufzeit zu verkürzen, werden die nicht verwendeten Elemente ausgeblendet.

        if ( params.cs==1 )
        {/*3*/  //Z=18453
            const double ccc1 = sqr(1-params.rho)*pow(params.p1,6);  //Z=18454
#ifdef UseFormfqF122
            const double ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,3);  //Z=18455
#endif
            const double ccc3 = sqr(params.rho);  //Z=18456
            const double vv3 = sqr((1-params.rho)*pow(params.p1,3)+params.rho);  //Z=18457

            const double argq = q*radiusm/(zr+1);  //Z=18459
#if defined(UseFormfqF121) || defined(UseFormfqF122)
            const double argpq = q*params.radius/(zr+1);  //Z=18460
#endif

#ifdef UseFormfqF121
            double F121;
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
#endif // UseFormfqF121

#ifdef UseFormfqF122
            double F122;
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
#endif // UseFormfqF122

            double F123;
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

            /* formfq:=F121+F122+F123;  //Z=18536 */
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
            const double FF1 = (cc1*F12)/sumc;  //Z=18745

            /*formfq:=*/ return FF1*FF1;  //Z=18747


            /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,3,sigmar,q);  //Z=18750 */


        }/*3*/  /*  of inhomogeneous core/shell sphere  */  //Z=18753

    //}/*2*/ /*  of sphere  */  //Z=18755

    return 0.0;
}  //Z=19829

#endif // SC_LIB_FORMFQ_partSphere_H
