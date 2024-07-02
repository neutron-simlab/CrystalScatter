#ifndef SC_LIB_FORMPQ_partSuperEllips_H
#define SC_LIB_FORMPQ_partSuperEllips_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partSuperEllips(double qxs, double qys, double q, int ordis) const   /*Z=14910*/
{/*1*/  //Z=15188
    // double sigmal ist nicht zu ersetzen, da es an einer Stelle CALC.epsilon ist, sonst nur CALC.params.sigmal
    // int ordis ist nicht zu ersetzen, da es einmal eine feste Zahl ist, sonst nur CALC.ordis

    double pqsum, oldpqsum, binsum, delser, pql;  //Z=15195
    double zr, argqx, argqy, pqrx, pqry;  //Z=15196

    zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232

    double pq;
    double qxn[121], qyn[121];

    CHECKENDTHREAD_VAL

    /* ** super ellipsoid, barrel ** */  //Z=18237
    //if ( params.part==7 )
    //{/*2*/  //Z=18238
        /*  homogeneous isotropic super ellipsoid  */  //Z=18239
        if ( ordis==7 )
        {/*3*/  //Z=18240
            if ( params.cs==0 )
            {/*4*/  //Z=18241
                if ( q<0.9*params.limq4 )
                {/*5*/  //Z=18242
                    pqsum = 1.0;  //Z=18243
                    oldpqsum = 0.0;  //Z=18244
                    double qqnn = 1.0;  //Z=18245
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18246
                        qqnn = qqnn*q*q;  //Z=18247
                        pqsum = pqsum+params.CR->carr4p[nser]*qqnn;  //Z=18248
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18249
                        if ( delser<0.0001 ) break; /* goto 270; */  //Z=18250
                        oldpqsum = pqsum;  //Z=18251
                    }/*6*/  //Z=18252
                    /*270:*/  //Z=18253
                    /*formpq:=*/ return pqsum;  //Z=18254
                }/*5*/  //Z=18255
                else
                {/*5*/  //Z=18256
                    //if ( q>=0.9*limq4 ) <<< macht keinen Sinn
                    pq = params.por/(q*q*q*q);  //Z=18257
                    /* else begin  //Z=18258 */
                    /*    qrombdeltac(length,radius,p1,sigma,dbeta,theta,0,qxs,qys,qz,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,ordis,3,7,14,7,0,0,carr1p,pql);  //Z=18259 */
                    /*    pq:=pql/(pi/2);  //Z=18260 */
                    /* end;  //Z=18261 */
                    /*formpq:=*/ return pq;  //Z=18262
                }/*5*/  //Z=18263
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18264

#ifdef procnotused
            /*  core/shell isotropic ellipsoid  */  //Z=18266
            if ( params.cs==1 )
            {/*4*/  //Z=18267
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18268
            }/*4*/  //Z=18269
#endif
        }/*3*/  /*  of isotropic triaxial ellipsoid  */  //Z=18270

        /*  perfectly oriented super ellipsoid  */  //Z=18272
        if ( ordis==6 )
        {/*3*/  //Z=18273
            /* if (orcase=1) then begin  //Z=18274 */
            if ( 1==1 )
            {/*4*/  //Z=18275

                if ( q<(1/params.radius) )
                {/*5*/  //Z=18277
                    pqsum = 1.0;  //Z=18278
                    oldpqsum = 0.0;  //Z=18279
                    qxn[0] = 1.0;  //Z=18280
                    qyn[0] = 1.0;  //Z=18281
                    for ( int nser=1; nser<=80; nser++ )
                    {/*6*/  //Z=18282
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18283
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18284

                        binsum = 0.0;  //Z=18286
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=18287
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18288 */
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=18289 */
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=18290
                        }/*7*/  //Z=18291
                        pqsum = pqsum+binsum;  //Z=18292
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18293
                        if ( delser<0.0001 ) break; /* goto 271; */  //Z=18294
                        oldpqsum = pqsum;  //Z=18295
                    }/*6*/  //Z=18296
                    /*271:*/  //Z=18297
                    pql = pqsum;  //Z=18298
                }/*5*/  //Z=18299
                else
                {/*5*/  //Z=18300
                    argqx = qxs*params.radius/(zr+1)+eps9;  //Z=18301
                    argqy = qys*params.radius/(zr+1)+eps9;  //Z=18302
                    pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18303
                    pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18304
                    pql = pqrx*pqry;  //Z=18305
                }/*5*/  //Z=18306
                /*formpq:=*/ return pql;  //Z=18307
            }/*4*/  /*  of orcase=1  */  //Z=18308
        }/*3*/  /*  of perfect triaxial ellipsoid  */  //Z=18309
    //}/*2*/  /*  of super ellipsoid  */  //Z=18310

    return 0.0;
}/*1*/  //Z=18387


#endif // SC_LIB_FORMPQ_partSuperEllips_H
