#ifndef SC_LIB_FORMPQ_partSuperball_H
#define SC_LIB_FORMPQ_partSuperball_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partSuperball(double qxs, double qys, double q, int ordis) const   /*Z=14910*/
{/*1*/  //Z=15188
    // double sigmal ist nicht zu ersetzen, da es an einer Stelle CALC.epsilon ist, sonst nur CALC.params.sigmal
    // int ordis ist nicht zu ersetzen, da es einmal eine feste Zahl ist, sonst nur CALC.ordis

    double pqsum, oldpqsum, binsum, delser, pql;  //Z=15195
    double zr, argqx, argqy, pqrx, pqry;  //Z=15196

    zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232

    double pq;
    double qxn[121], qyn[121];

    CHECKENDTHREAD_VAL

    /* ** superball ** */  //Z=18312
    //if ( params.part==8 )
    //{/*2*/  //Z=18313
        /*  homogeneous isotropic superball  */  //Z=18314
        if ( ordis==7 )
        {/*3*/  //Z=18315
            if ( params.cs==0 )
            {/*4*/  //Z=18316
                if ( q<1.1*params.limq4 )
                {/*5*/  //Z=18317
                    pqsum = 1.0;  //Z=18318
                    oldpqsum = 0.0;  //Z=18319
                    double qqnn = 1.0;  //Z=18320
                    for ( int nser=1; nser<=40; nser++ )
                    {/*6*/  //Z=18321
                        qqnn = qqnn*q*q;  //Z=18322
                        pqsum = pqsum+params.CR->carr4p[nser]*qqnn;  //Z=18323
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18324
                        if ( delser<0.0001 ) break; /* goto 274; */  //Z=18325
                        oldpqsum = pqsum;  //Z=18326
                    }/*6*/  //Z=18327
                    /*274:*/  //Z=18328
                    /*formpq:=*/ return pqsum;  //Z=18329
                }/*5*/  //Z=18330
                else
                {/*5*/  //Z=18331
                    //if ( q>=1.1*limq4 ) <<< macht keinen Sinn
                    pq = params.por/(q*q*q*q);  //Z=18332
                    /* else begin  //Z=18333 */
                    /*    qrombdeltac(length,radius,p1,sigma,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,7,14,7,0,0,carr1p,pql);  //Z=18334 */
                    /*    pq:=pql/(pi/2);  //Z=18335 */
                    /* end;  //Z=18336 */
                    /*formpq:=*/ return pq;  //Z=18337
                }/*5*/  //Z=18338
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18339

#ifdef procnotused
            /*  core/shell isotropic superball  */  //Z=18341
            if ( params.cs==1 )
            {/*4*/  //Z=18342
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18343
            }/*4*/  //Z=18344
#endif
        }/*3*/  /*  of isotropic superball  */  //Z=18345

        /*  perfectly oriented superball  */  //Z=18347
        if ( ordis==6 )
        {/*3*/  //Z=18348
            /* if (orcase=1) then begin  //Z=18349 */
            if ( 1==1 )
            {/*4*/  //Z=18350

                if ( q<(1/params.radius) )
                {/*5*/  //Z=18352
                    pqsum = 1.0;  //Z=18353
                    oldpqsum = 0.0;  //Z=18354
                    qxn[0] = 1.0;  //Z=18355
                    qyn[0] = 1.0;  //Z=18356
                    for ( int nser=1; nser<=80; nser++ )
                    {/*6*/  //Z=18357
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18358
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18359

                        binsum = 0.0;  //Z=18361
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=18362
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18363 */
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=18364 */
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=18365
                        }/*7*/  //Z=18366
                        pqsum = pqsum+binsum;  //Z=18367
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18368
                        if ( delser<0.0001 ) break; /* goto 275; */  //Z=18369
                        oldpqsum = pqsum;  //Z=18370
                    }/*6*/  //Z=18371
                    /*275:*/  //Z=18372
                    pql = pqsum;  //Z=18373
                }/*5*/  //Z=18374
                else
                {/*5*/  //Z=18375
                    argqx = qxs*params.radius/(zr+1)+eps9;  //Z=18376
                    argqy = qys*params.radius/(zr+1)+eps9;  //Z=18377
                    pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18378
                    pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18379
                    pql = pqrx*pqry;  //Z=18380
                }/*5*/  //Z=18381
                /*formpq:=*/ return pql;  //Z=18382
            }/*4*/  /*  of orcase=1  */  //Z=18383
        }/*3*/  /*  of perfect triaxial ellipsoid  */  //Z=18384
    //}/*2*/  /*  of superball  */  //Z=18385

    return 0.0;
}/*1*/  //Z=18387


#endif // SC_LIB_FORMPQ_partSuperball_H
