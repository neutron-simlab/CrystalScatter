#ifndef SC_LIB_FORMPQ_partCylinder_H
#define SC_LIB_FORMPQ_partCylinder_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partCylinder(double qxs, double qys, double q) const   /*Z=14910*/
{/*1*/  //Z=15188

    double pqr, pql;  //Z=15195

    const double zl = params.sigmal==1.0 ? eps9 : (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=15231
    const double zr = params.sigma ==1.0 ? eps9 : (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232
    const double radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=15233

    // Die Konstante 'eps' muss durch 'eps9' ersetzt werden.
    const double zz = zr;
    const double z  = zl;

    CHECKENDTHREAD_VAL;

    /* ********** */  //Z=16072
    /*  cylinder  */  //Z=16073
    /* ********** */  //Z=16074
    //    if ( params.part==1 )
    //    {/*2*/  //Z=16075

        /* ** longitudinal part ** */  //Z=16077
        /* ** isotropic ** */  //Z=16078
        if ( params.ordis==7 )
        {/*3*/  //Z=16079
            /*  exact average  */  //Z=16080
            if ( (params.length/params.radius)<2 )
            {/*4*/  //Z=16081
                if ( q<(5*params.limq2) )
                {  //Z=16082
                    /*  double sum  */  //Z=16095
                    double qqn[101];
                    qqn[0] = 1.0;  //Z=16096  Hier ist das Array sinnvoll, aber kleiner und lokal
                    for ( int nser=1; nser<=100; nser++ ) qqn[nser] = qqn[nser-1]*q*q;  //Z=16097
                    pql = 0.0;  //Z=16098
                    double oldpqsum = -10.0;  //Z=16099
                    for ( int nser=0; nser<=100; nser++ )
                    {  //Z=16100
                        double binsum = 0.0;  //Z=16101
                        for ( int mser=0; mser<=100; mser++ ) binsum += params.CR->carr11pm[nser][mser]*qqn[mser];  //Z=16102
                        pql += params.CR->carr2p[nser]*qqn[nser]*binsum;  //Z=16103
                        double delser = fabs((pql-oldpqsum)/pql);  //Z=16104
                        if ( delser<0.0001 ) break; /* goto 601; */  //Z=16105
                        oldpqsum = pql;  //Z=16106
                    }  //Z=16107

                    /*601:*/  //Z=16109
                    //pql = pqsum;  //Z=16110
                }  //Z=16111
                else
                {  //Z=16112
                    pql = params.por/pow(q,4);  //Z=16113
                }  //Z=16114
            }/*4*/  //Z=16115

            /*  factorization  */  //Z=16117
            else
            {/*4*/  //Z=16118
                if ( q<(0.6*params.limq1) )
                {  //Z=16119
                    pql = 1.0;  //Z=16120
                    double oldpqsum = 0.0;  //Z=16121
                    double qqnn = 1.0;  //Z=16122
                    for ( int nser=1; nser<=120; nser++ )
                    {  //Z=16123
                        qqnn = qqnn*q*q;  //Z=16124
                        pql += params.CR->carr1p[nser]*qqnn;  //Z=16125
                        double delser = fabs((pql-oldpqsum)/pql);  //Z=16126
                        if ( delser<0.0001 ) break; /* goto 60; */  //Z=16127
                        oldpqsum = pql;  //Z=16128
                    }  //Z=16129
                    /*60:*/  //Z=16130
                    //pql = pqsum;  //Z=16131
                }  //Z=16132
                else
                {  //Z=16133
                    const double arglq = q*params.length/(zl+1);  //Z=16134
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2));  //Z=16135 */
                    pql = (M_PI/(2.0*zl))*(1/arglq);  //Z=16136
                    pql = pql-(1/(2.0*zl*(zl-1)*arglq*arglq))*cos((zl-1)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-1)/2.0);  //Z=16137
                }  //Z=16138
            }/*4*/  //Z=16139
        }/*3*/   /*  of isotropic  */  //Z=16140

        /*  perfect  */  //Z=16142
        if ( params.ordis==6 )
        {/*3*/  //Z=16143
            if ( params.orcase==4 )
                pql = 1.0;  //Z=16144
            else
            {/*4*/  //Z=16145
                if ( q<(0.6*params.limq1) ) //20240301 - war ( limql<(0.5*params.limq1) )
                {/*5*/  //Z=16146
                    /* if (sqrt(qx*qx*length*length+qy*qy*radius*radius+eps)<10) then begin  //Z=16147 */
                    pql = 1.0;  //Z=16148
                    double oldpqsum = 0.0;  //Z=16149
                    double qqnn = 1.0;  //Z=16150
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=16151
                        qqnn = qqnn*(qxs+qys)*(qxs+qys);     /*  (qxs,qys)=(qx,0) for x, (0,qy) for y, (qx,qy) for z  */  //Z=16152
                        pql += params.CR->carr1p[nser]*qqnn;  //Z=16153
                        const double delser = fabs((pql-oldpqsum)/pql);  //Z=16154
                        if ( delser<0.0001 ) break; /* goto 65; */  //Z=16155
                        oldpqsum = pql;  //Z=16156
                    }/*6*/  //Z=16157
                    /*65:*/  //Z=16158
                    //pql = pqsum;  //Z=16159
                }/*5*/  //Z=16160
                else
                {/*5*/  //Z=16161
                    const double arglq = (qxs+qys+eps9)*params.length/(zl+1);  //Z=16162
                    /* pql:=(1/(2*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*arctan(2*arglq))/power(1+4*arglq*arglq,(zl-1)/2));  //Z=16163 */
                    /* pql:=(pi/(2*zl))*(1/arglq);  //Z=16164 */
                    pql = (1/(2.0*zl*(zl-1)))*(1/(arglq*arglq))*(1-cos((zl-1)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-1)/2.0));  //Z=16165
                }/*5*/  //Z=16166
            }/*4*/  //Z=16167
        }/*3*/   /*  of perfect  */  //Z=16168

        /*  general  */  //Z=16170
        if ( params.ordis==0 )
        {/*3*/  //Z=16171
            if ( params.orcase==4 )
                pql = 1.0;  //Z=16172
            else
            {/*4*/  //Z=16173
                if ( q<(0.4*params.limq1) ) //20240301 - war ( limql<(2.0*params.limq1) )
                {/*5*/  //Z=16175
                    /* if (limql<0.1) then begin  //Z=16176 */
                    /* if (q<0.05) then begin  //Z=16177 */
                    pql = 1.0;  //Z=16178
                    double oldpqsum = 0.0;  //Z=16179
                    double qxn[121], qyn[121];
                    qxn[0] = 1.0;  //Z=16180
                    qyn[0] = 1.0;  //Z=16181
                    if ( params.orcase==1 )
                    {/*6*/  //Z=16182
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16183
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=16184
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=16185
                            double binsum = 0.0;  //Z=16186
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16187
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16188 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=16189 */
                                binsum += params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];  //Z=16190 //20240301
                            }/*8*/  //Z=16191
                            pql += params.CR->carr1p[nser]*binsum;  //Z=16192
                            double delser = fabs((pql-oldpqsum)/pql);  //Z=16193
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=16194
                            oldpqsum = pql;  //Z=16195
                        }/*7*/  //Z=16196
                    }/*6*/  //Z=16197

                    else if ( params.orcase==2 )
                    {/*6*/  /*  x-axis  */  //Z=16199
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16200
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=16201
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=16202
                            double binsum = 0.0;  //Z=16203
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16204
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16205 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=16206 */
                                binsum += params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];  //Z=16207 //20240301
                            }/*8*/  //Z=16208
                            pql += params.CR->carr1p[nser]*binsum;  //Z=16209
                            double delser = fabs((pql-oldpqsum)/pql);  //Z=16210
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=16211
                            oldpqsum = pql;  //Z=16212
                        }/*7*/  //Z=16213
                    }/*6*/  //Z=16214

                    else if ( params.orcase==3 )
                    {/*6*/  /*  y-axis  */  //Z=16217
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16218
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=16219
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=16220
                            double binsum = 0.0;  //Z=16221
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16222
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16223 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=16224 */
                                binsum += params.CR->carr11pm[nser-mser][mser]*qyn[mser]*qxn[nser-mser];  //Z=16225 //20240301
                            }/*8*/  //Z=16226
                            pql += params.CR->carr1p[nser]*binsum;  //Z=16227
                            double delser = fabs((pql-oldpqsum)/pql);  //Z=16228
                            if ( delser<0.0001 ) break; /* goto 66; */  //Z=16229
                            oldpqsum = pql;  //Z=16230
                        }/*7*/  //Z=16231
                    }/*6*/  //Z=16232
                    /*66:*/  //Z=16233
                    //pql = pqsum;  //Z=16234
                }/*5*/  //Z=16235
                else
                {/*5*/  //Z=16236
                    double qz=1.0;
                    qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,0,qxs,qys,qz, // 9,9,
                                9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,1,4,params.orcase,0,0,0,params.CR->carr1p,pql);  //Z=16237
                    pql = pql/params.norm;  //Z=16238
                }/*5*/  //Z=16239
            }/*4*/  //Z=16240
        }/*3*/   /*  of general  */  //Z=16241


        /*  transverse part  */  //Z=16244
        /*  homogeneous cylinder  */  //Z=16245
        if ( params.cs==0 )
        {/*3*/  //Z=16246
            /*  exact average  */  //Z=16247
            if ( (params.length/params.radius)<2 )
                pqr = 1;  //Z=16248
            /*  factorization  */  //Z=16249
            else
            {/*4*/  //Z=16250
                if ( q<(0.65*params.limq4) ) //20240301 - war 0.3
                {/*5*/  //Z=16251
                    /* if (sqrt(qx*qx*length*length+qy*qy*radius*radius+eps)<10) then begin  //Z=16252 */
                    pqr = 1.0;  //Z=16253
                    double oldpqsum = 0.0;  //Z=16254
                    double qqnn = 1.0;  //Z=16255
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=16256
                        qqnn = qqnn*q*q;  //Z=16257
                        /* qqn[nser]:=qqn[nser-1]*(qxs+qys)*(qxs+qys);     (* (qxs,qys)=(qx,0) for x, (0,qy) for y, (qx,qy) for z *)  //Z=16258 */
                        /* qqn[nser]:=qqn[nser-1]*qx*qx;       (* GISAXS *)  //Z=16259 */
                        pqr += params.CR->carr4p[nser]*qqnn;  //Z=16260
                        const double delser = fabs((pqr-oldpqsum)/pqr);  //Z=16261
                        if ( delser<0.0001 ) break; /* goto 61; */  //Z=16262
                        oldpqsum = pqr;  //Z=16263
                    }/*6*/  //Z=16264
                    /*61:*/  //Z=16265
                    //pqr = pqsum;  //Z=16266
                }/*5*/  //Z=16267
                else
                {/*5*/  //Z=16268
                    /* q:=sqrt(qxs*qxs+qys*qys+eps);       (* for perfect orientation *)  //Z=16269 */
                    /* q:=abs(qx+eps);                               (* for GISAXS *)  //Z=16270 */
                    const double argpq = q*params.radius/(zr+1);  //Z=16271
                    const double pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3);  //Z=16272
                    const double pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3)*sin((zr-2)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-2)/2.0);  //Z=16273
                    const double pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argpq,-4)*cos((zr-3)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-3)/2.0);  //Z=16274
                    pqr = (4./M_PI)*(pqr1-pqr2-(9./8.0)*pqr3);  //Z=16275
                    /*  add more terms, if necessary  */  //Z=16276
                }/*5*/  //Z=16277
            }/*4*/  //Z=16278
            /*formpq:=*/ return pql*pqr;  //Z=16279
        }/*3*/ /*  of homogeneous  */  //Z=16282

        /*  core/shell cylinder  */  //Z=16284
        if ( params.cs==1 )
        {/*3*/  //Z=16285
            //const double cc1 = sqr(params.rho);  //Z=16286
            //const double cc2 = 2*params.p1*params.rho*(1-params.rho);  //Z=16287
            //const double cc3 = sqr(1-params.rho)*sqr(params.p1);  //Z=16288
            //const double cc4 = -2*sqr(params.rho);  //Z=16289
            //const double cc5 = -2*params.p1*params.rho*(1-params.rho);  //Z=16290
            //const double cc6 = sqr(params.rho);  //Z=16291
            //cc7 = -2*params.rho*(1-params.rho);  //Z=16292
            //cc8 = -sqr(1-params.rho)*2*params.p1;  //Z=16293
            //cc9 = 2*params.rho*(1-params.rho);  //Z=16294
            //cc10 = sqr(1-params.rho);  //Z=16295

            const double ccc1 = sqr(1-params.rho)*pow(params.p1,4);  //Z=16297
            const double ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,2);  //Z=16298
            const double ccc3 = params.rho*params.rho;  //Z=16299
            const double vv3 = sqr((1-params.rho)*pow(params.p1,2)+params.rho);  //Z=16300

            const double argq = q*radiusm/(zz+1);  //Z=16302
            const double argpq = q*params.radius/(zz+1);  //Z=16303

            double F121, F122, F123;

            /*  F121 cylinder  */  //Z=16305
            if ( q<(0.7*params.limq4) )
            {/*4*/  //Z=16306
                /* ** series expansion ** */  //Z=16307
                double pqsum = 1.0;  //Z=16308
                double oldpqsum = 0.0;  //Z=16309
                double qqnn = 1.0;  //Z=16310
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=16311
                    qqnn = qqnn*q*q;  //Z=16312
                    pqsum = pqsum+params.CR->carr4p[nser]*qqnn;  //Z=16313
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16314
                    if ( delser<0.0001 ) break; /* goto 62; */  //Z=16315
                    oldpqsum = pqsum;  //Z=16316
                }/*5*/  //Z=16317
                /*62:*/  //Z=16318
                F121 = ccc1*pqsum/vv3;  //Z=16319
            }/*4*/  //Z=16320
            else
            {/*4*/  //Z=16321
                const double pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3);  //Z=16322
                const double pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argpq,-3)*sin((zr-2)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-2)/2.0);  //Z=16323
                const double pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argpq,-4)*cos((zr-3)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-3)/2.0);  //Z=16324
                const double pqr = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=16325
                F121 = ccc1*pqr/vv3;  //Z=16326
                /*  add more terms, if necessary  */  //Z=16327
            }/*4*/  //Z=16328

            /*  F122 cylinder  */  //Z=16330
            if ( q<(1.5*params.limq5) )
            {/*4*/  //Z=16331
                /* ** series expansion ** */  //Z=16332
                double pqsum = 1.0;  //Z=16333
                double oldpqsum = 0.0;  //Z=16334
                double qqnn = 1.0;  //Z=16335
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=16336
                    qqnn = qqnn*q*q;  //Z=16337
                    pqsum = pqsum+params.CR->carr5p[nser]*qqnn;  //Z=16338
                    double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16339
                    if ( delser<0.0001 ) break; /* goto 63; */  //Z=16340
                    oldpqsum = pqsum;  //Z=16341
                }/*5*/  //Z=16342
                /*63:*/  //Z=16343
                F122 = ccc2*pqsum/vv3;  //Z=16344
            }/*4*/  //Z=16345
            else
            {/*4*/  //Z=16346
                const double argbm = (zr-2)*atan(argpq-argq);  //Z=16347
                const double nenbm = pow(1.0+sqr(argpq-argq),(zr-2)/2.0);  //Z=16348
                const double argbp = (zr-2)*atan(argpq+argq);  //Z=16349
                const double nenbp = pow(1.0+sqr(argpq+argq),(zr-2)/2.0);  //Z=16350
                const double argem = (zr-3)*atan(argpq-argq);  //Z=16351
                const double nenem = pow(1.0+sqr(argpq-argq),(zr-3)/2.0);  //Z=16352
                const double argep = (zr-3)*atan(argpq+argq);  //Z=16353
                const double nenep = pow(1.0+sqr(argpq+argq),(zr-3)/2.0);  //Z=16354
                const double arggm = (zr-4)*atan(argpq-argq);  //Z=16355
                const double nengm = pow(1.0+sqr(argpq-argq),(zr-4)/2.0);  //Z=16356
                const double arggp = (zr-4)*atan(argpq+argq);  //Z=16357
                const double nengp = pow(1.0+sqr(argpq+argq),(zr-4)/2.0);  //Z=16358

                const double pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3)*(cos(argbm)/nenbm-sin(argbp)/nenbp);  //Z=16360
                const double pqr2 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4)*((1-1/params.p1)*sin(argem)/nenem-(1+1/params.p1)*cos(argep)/nenep);  //Z=16361
                const double pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)*(zr-4)))*pow(argq,-5)*(1/params.p1)*(cos(arggm)/nengm-sin(arggp)/nengp);  //Z=16362
                const double pqr = (4/M_PI)*pow(params.p1,-3/2.0)*(pqr1+(9/16.0)*pqr2+(9/16.0)*(9/16.0)*pqr3);  //Z=16363
                F122 = ccc2*pqr/vv3;  //Z=16364
            }/*4*/  //Z=16365

            /*  F123 cylinder  */  //Z=16367
            if ( q<(0.6*params.limq6) )
            {/*4*/  //Z=16368
                /* ** series expansion ** */  //Z=16369
                double pqsum = 1.0;  //Z=16370
                double oldpqsum = 0.0;  //Z=16371
                double qqnn = 1.0;  //Z=16372
                for ( int nser=1; nser<=120; nser++ )
                {/*5*/  //Z=16373
                    qqnn = qqnn*q*q;  //Z=16374
                    pqsum = pqsum+params.CR->carr6p[nser]*qqnn;  //Z=16375
                    const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=16376
                    if ( delser<0.0001 ) break; /* goto 64; */  //Z=16377
                    oldpqsum = pqsum;  //Z=16378
                }/*5*/  //Z=16379
                /*64:*/  //Z=16380
                F123 = ccc3*pqsum/vv3;  //Z=16381
            }/*4*/  //Z=16382
            else
            {/*4*/  //Z=16383
                const double pqr1 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3);  //Z=16384
                const double pqr2 = (1/(zr*(zr-1)*(zr-2)))*pow(argq,-3)*sin((zr-2)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-2)/2.0);  //Z=16385
                const double pqr3 = (1/(zr*(zr-1)*(zr-2)*(zr-3)))*pow(argq,-4)*cos((zr-3)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-3)/2.0);  //Z=16386
                const double pqr = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=16387
                F123 = ccc3*pqr/vv3;  //Z=16388
                /*  add more terms, if necessary  */  //Z=16389
            }/*4*/  //Z=16390
            /*formpq:=*/ return pql*(F121+F122+F123);  //Z=16391
            /* formpq:=F122;  //Z=16392 */
        }/*3*/ /*  of core/shell  */  //Z=16393

        /*  inhomogeneous core/shell cylinder  */  //Z=16395
        if ( params.cs==2 )
        {/*3*/  //Z=16396

            const double dim = 2;  //Z=16398
            const double delc = 0.0001;  //Z=16399
            const double xrad = q*radiusm;  //Z=16400
            const double xradp = q*params.radius;  //Z=16401
            const double x1z = q*params.radius/(2.0*(zr+1));  //Z=16402
            const double x12z = x1z*x1z;  //Z=16403
            const double x2z = q*radiusm/(2.0*(zr+1));  //Z=16404
            const double x22z = x2z*x2z;  //Z=16405

            const double lim = 18*exp(-5*params.sigma);  //Z=16407
            const double lim1 = lim;  //Z=16408
            const double lim2 = lim*0.7;  //Z=16409
            const double lim3 = lim;  //Z=16410
            const double lim4 = lim;  //Z=16411
            const double lim5 = lim*0.7;  //Z=16412
            const double lim6 = lim*1.2;  //Z=16413

            const double a1 = (dim-params.alphash1)/2.0;  //Z=16415
            const double b1 = dim/2.0;  //Z=16416
            const double b2 = (dim+2-params.alphash1)/2.0;  //Z=16417
            const double b1s = (dim+2)/2.0;  //Z=16418
            const double v = -b1s+1/2.0;  //Z=16419
            const double c = a1-b1-b2+1/2.0;  //Z=16420
            const double d0 = 1;  //Z=16421
            //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=16422
            const double e0 = 1.0;  //Z=16423
            const double e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=16424
            const double ee0 = 1.0;  //Z=16425
            const double ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=16426

            const double gb1s = 1;  //Z=16428
            const double pz2v = 1/(zr*(zr-1)*(zr-2));  //Z=16429
            const double pz2v1 = pz2v/(zr-3);  //Z=16430
            const double pz2v2 = pz2v1/(zr-4);  //Z=16431

            const double gz1 = gamma(zr+1);  //Z=16433
            const double preg1 = gb1s/sqrt(M_PI);  //Z=16434
            const double preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=16435
            const double preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=16436
            const double pzvc = gamma(zr+1+v+c)/gz1;  //Z=16437
            const double pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=16438
            const double pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=16439
            const double pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=16440
            const double pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=16441
            //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=16442
            const double pzc = gamma(zr+1+2*c)/gz1;  //Z=16443
            const double pzc1 = gamma(zr+1+2*c-1)/gz1;  //Z=16444
            const double pza = gamma(zr+1-4*a1)/gz1;  //Z=16445
            const double pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=16446
            const double pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=16447
            //dnv0 = 1;  //Z=16448
            //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=16449
            //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=16450
            //pva0 = gamma(zr+1-4*a1)/gz1;  //Z=16451

            const double cc1 = 1/(dim*dim);  //Z=16453
            const double cc2 = 2*params.rho/(dim*(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=16454
            const double cc3 = -2*params.rho/(dim*(dim-params.alphash1));  //Z=16455
            const double cc4 = sqr(params.rho)/(sqr(dim-params.alphash1)*pow(sqr(params.p1),dim-params.alphash1));  //Z=16456
            const double cc5 = -2*sqr(params.rho)/(sqr(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=16457
            const double cc6 = sqr(params.rho)/sqr(dim-params.alphash1);  //Z=16458
            const double vv3 = cc1+cc2+cc3+cc4+cc5+cc6;  //Z=16459

            double F12, F22, F32, F42, F52, F62;

            /*  term #1 series  */  //Z=16461
            if ( (xradp)<lim1 )
            {/*4*/  //Z=16462
                //z12v[0] = 1;  //Z=16463
                //b1sv[0] = 1;  //Z=16464
                //fkv[0] = 1;  //Z=16465
                F12 = 1.0;  //Z=16466
                double oldF12 = 0.0;  //Z=16467
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16468
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=16469
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16470
                    qqnn = qqnn * sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16471
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16472
                    //fkv[n] = fkv[n-1]*n;  //Z=16473
                    //sum12[n] = 0;  //Z=16474
                    /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);  //Z=16475 */
                    /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n];  //Z=16476 */

                    F12 += params.CR->carr4p[n]*qqnn; //qqn[n];  //Z=16478

                    const double del = fabs((F12-oldF12)/F12);  //Z=16480
                    if ( del<delc ) break; /* goto 211; */  //Z=16481
                    oldF12 = F12;  //Z=16482
                }/*5*/  //Z=16483
                /*211:*/  //Z=16484
                //F12 = F12sez;  //Z=16485
            }/*4*/  //Z=16486

            /*  term #2 series  */  //Z=16488
            if ( (xradp)<lim2 )
            {/*4*/  //Z=16489
                double a1v[121], b1v[121], b2v[121], b1sv[121], fkv[121];
                double z12vn = 1;  //Z=16490
                a1v[0] = 1;  //Z=16491
                b1v[0] = 1;  //Z=16492
                b2v[0] = 1;  //Z=16493
                b1sv[0] = 1;  //Z=16494
                fkv[0] = 1;  //Z=16495
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16496
                F22 = 1.0;  //Z=16497
                double oldF22 = 0.0;  //Z=16498
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=16499
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16500
                    qqnn = qqnn * sqr(q);
                    z12vn = z12vn*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16501
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16502
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16503
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16504
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16505 TODO
                    fkv[n] = fkv[n-1]*n;  //Z=16506
                    double sum22 = 0;  //Z=16507
                    for ( int m=0; m<=n; m++ ) sum22 += a1v[n-m]*pow(sqr(params.p1),m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16508
                    F22 += pow(-x22z,n)*z12vn*sum22;  //Z=16509

                    /* F22sez:=F22sez+carr5p[n]*qqn[n];  //Z=16511 */

                    const double del = fabs((F22-oldF22)/F22);  //Z=16513
                    if ( del<delc ) break; /* goto 212; */  //Z=16514
                    oldF22 = F22;  //Z=16515
                }/*5*/  //Z=16516
                /*212:*/  //Z=16517
                //F22 = F22sez;  //Z=16518
            }/*4*/  //Z=16519

            /*  term #3 series  */  //Z=16521
            if ( (xradp)<lim3 )
            {/*4*/  //Z=16522
                double a1v[121], b1v[121], b2v[121], b1sv[121], fkv[121];
                double z12vn = 1;  //Z=16523
                a1v[0] = 1;  //Z=16524
                b1v[0] = 1;  //Z=16525
                b2v[0] = 1;  //Z=16526
                b1sv[0] = 1;  //Z=16527
                fkv[0] = 1;  //Z=16528
                //qqn[0] = 1.0;  //Z=16529
                F32 = 1.0;  //Z=16530
                double oldF32 = 0.0;  //Z=16531
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=16532
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16533
                    z12vn = z12vn*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16534
                    a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16535
                    b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16536
                    b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16537
                    b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16538 TODO
                    fkv[n] = fkv[n-1]*n;  //Z=16539
                    double sum32 = 0;  //Z=16540
                    for ( int m=0; m<=n; m++ ) sum32 += a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16541
                    F32 += pow(-x12z,n)*z12vn*sum32;  //Z=16542

                    /* F32sez:=F32sez+carr6p[n]*qqn[n];  //Z=16544 */

                    const double del = fabs((F32-oldF32)/F32);  //Z=16546
                    if ( del<delc ) break; /* goto 213; */  //Z=16547
                    oldF32 = F32;  //Z=16548
                }/*5*/  //Z=16549
                /*213:*/  //Z=16550
                //F32 = F32sez;  //Z=16551
            }/*4*/  //Z=16552

            /*  term #4 series  */  //Z=16554
            if ( (xradp)<lim4 )
            {/*4*/  //Z=16555
                //z12v[0] = 1;  //Z=16556
                //a1v[0] = 1;  //Z=16557
                //b1v[0] = 1;  //Z=16558
                //b2v[0] = 1;  //Z=16559
                //b1sv[0] = 1;  //Z=16560
                //fkv[0] = 1;  //Z=16561
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16562
                F42 = 1.0;  //Z=16563
                double oldF42 = 0.0;  //Z=16564
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=16565
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16566
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16567
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16568
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16569
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16570
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16571
                    //fkv[n] = fkv[n-1]*n;  //Z=16572
                    //sum42[n] = 0;  //Z=16573
                    /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16574 */
                    /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n];  //Z=16575 */

                    F42 += params.CR->carr7p[n]*qqnn; //qqn[n];  //Z=16577

                    const double del = fabs((F42-oldF42)/F42);  //Z=16579
                    if ( del<delc ) break; /* goto 214; */  //Z=16580
                    oldF42 = F42;  //Z=16581
                }/*5*/  //Z=16582
                /*214:*/  //Z=16583
                //F42 = F42sez;  //Z=16584
            }/*4*/  //Z=16585

            /*  term #5 series  */  //Z=16587
            if ( (xradp)<lim5 )
            {/*4*/  //Z=16588
                //z12v[0] = 1;  //Z=16589
                //a1v[0] = 1;  //Z=16590
                //b1v[0] = 1;  //Z=16591
                //b2v[0] = 1;  //Z=16592
                //b1sv[0] = 1;  //Z=16593
                //fkv[0] = 1;  //Z=16594
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16595
                F52 = 1.0;  //Z=16596
                double oldF52 = 0.0;  //Z=16597
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=16598
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16599
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16600
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16601
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16602
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16603
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16604
                    //fkv[n] = fkv[n-1]*n;  //Z=16605
                    //sum52[n] = 0;  //Z=16606
                    /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16607 */
                    /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n];  //Z=16608 */

                    F52 += params.CR->carr8p[n]*qqnn; //qqn[n];  //Z=16610

                    const double del = fabs((F52-oldF52)/F52);  //Z=16612
                    if ( del<delc ) break; /* goto 215; */  //Z=16613
                    oldF52 = F52;  //Z=16614
                }/*5*/  //Z=16615
                /*215:*/  //Z=16616
                //F52 = F52sez;  //Z=16617
            }/*4*/  //Z=16618

            /*  term #6 series  */  //Z=16620
            if ( (xradp)<lim6 )
            {/*4*/  //Z=16621
                //z12v[0] = 1;  //Z=16622
                //a1v[0] = 1;  //Z=16623
                //b1v[0] = 1;  //Z=16624
                //b2v[0] = 1;  //Z=16625
                //b1sv[0] = 1;  //Z=16626
                //fkv[0] = 1;  //Z=16627
                double qqnn = 1.0; //qqn[0] = 1.0;  //Z=16628
                F62 = 1.0;  //Z=16629
                double oldF62 = 0.0;  //Z=16630
                for ( int n=1; n<=120; n++ )
                {/*5*/  //Z=16631
                    //qqn[n] = qqn[n-1]*q*q;  //Z=16632
                    qqnn *= sqr(q);
                    //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=16633
                    //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=16634
                    //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=16635
                    //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=16636
                    //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=16637
                    //fkv[n] = fkv[n-1]*n;  //Z=16638
                    //sum62[n] = 0;  //Z=16639
                    /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=16640 */
                    /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n];  //Z=16641 */

                    F62 += params.CR->carr9p[n]*qqnn; //qqn[n];  //Z=16643

                    const double del = fabs((F62-oldF62)/F62);  //Z=16645
                    if ( del<delc ) break; /* goto 216; */  //Z=16646
                    oldF62 = F62;  //Z=16647
                }/*5*/  //Z=16648
                /*216:*/  //Z=16649
                //F62 = F62sez;  //Z=16650
            }/*4*/  //Z=16651


            /* ** term #1 asymptote ** */  //Z=16654
            if ( xradp>=lim1 )
            {/*4*/  //Z=16655
                const double arg11 = (zr+2*v+1)*atan(4.0*x1z);  //Z=16656
                const double nen11 = pow(1.0+16*x1z*x1z,(zr+2*v+1)/2.0);  //Z=16657
                const double arg12 = (zr+2*v)*atan(4.0*x1z);  //Z=16658
                const double nen12 = pow(1.0+16*x1z*x1z,(zr+2*v)/2.0);  //Z=16659
                const double arg13 = (zr+2*v-1)*atan(4.0*x1z);  //Z=16660
                const double nen13 = pow(1.0+16*x1z*x1z,(zr+2*v-1)/2.0);  //Z=16661

                const double F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);  //Z=16663
                const double F12as2z = 2*ee0*ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);  //Z=16664
                const double F12as3z = ee1*ee1*(1/(4.0*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);  //Z=16665
                F12 = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);  //Z=16666
                //F12 = F12asz;  //Z=16667
            }/*4*/  //Z=16668

            /* ** term #2 asymptote ** */  //Z=16670
            if ( xradp>=lim2 )
            {/*4*/  //Z=16671
                //arg21 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=16672
                //nen21 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=16673
                //arg22 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=16674
                //nen22 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=16675
                //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);  //Z=16676
                //F22as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);  //Z=16677
                const double F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);  //Z=16678
                //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);  //Z=16679

                const double arg210 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=16681
                const double nen210 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=16682
                const double arg220 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=16683
                const double nen220 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=16684
                const double F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);  //Z=16685
                const double F22as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);  //Z=16686
                const double F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);  //Z=16687
                const double arg23 = (zr+v+c+1)*atan(2.0*(x1z-x2z));  //Z=16688
                const double nen23 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);  //Z=16689
                const double arg24 = (zr+v+c+1)*atan(2.0*(x1z+x2z));  //Z=16690
                const double nen24 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);  //Z=16691
                const double arg25 = (zr+v+c)*atan(2.0*(x1z-x2z));  //Z=16692
                const double nen25 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);  //Z=16693
                const double arg26 = (zr+v+c)*atan(2.0*(x1z+x2z));  //Z=16694
                const double nen26 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);  //Z=16695
                const double arg27 = (zr+v+c-1)*atan(2.0*(x1z-x2z));  //Z=16696
                const double nen27 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);  //Z=16697
                const double arg28 = (zr+v+c-1)*atan(2.0*(x1z+x2z));  //Z=16698
                const double nen28 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);  //Z=16699

                const double a22as21z = (1/2.0)*ee0*e0*pzvc;  //Z=16701
                const double F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);  //Z=16702
                const double a22as22z = (1/2.0)*ee0*e1*(1/(2.0*x2z))*pzvc1;  //Z=16703
                const double F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);  //Z=16704
                const double a22as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1;  //Z=16705
                const double F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);  //Z=16706
                const double a22as24z = (1/2.0)*ee1*e1*(1/(2.0*x1z))*(1/(2.0*x2z))*pzvc2;  //Z=16707
                const double F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);  //Z=16708
                const double F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);  //Z=16709
                const double F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);  //Z=16710
                //F22asz = F22as1z+F22as2z;  //Z=16711
                F22 = F22as1z0+F22as2z;  //Z=16712
                //F22 = F22asz0;  //Z=16713
            }/*4*/  //Z=16714

            /* ** term #3 asymptote ** */  //Z=16716
            if ( xradp>=lim3 )
            {/*4*/  //Z=16717
                //arg31 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=16718
                //nen31 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=16719
                //arg32 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=16720
                //nen32 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=16721
                //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);  //Z=16722
                //F32as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);  //Z=16723
                const double F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);  //Z=16724
                //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);  //Z=16725

                const double arg310 = (z+v-2*a1+1)*atan(2.0*x1z);  //Z=16727
                const double nen310 = pow(1.0+4*x1z*x1z,(z+v-2*a1+1)/2.0);  //Z=16728
                const double arg320 = (z+v-2*a1)*atan(2.0*x1z);  //Z=16729
                const double nen320 = pow(1.0+4*x1z*x1z,(z+v-2*a1)/2.0);  //Z=16730
                const double F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);  //Z=16731
                const double F32as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);  //Z=16732
                const double F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);  //Z=16733

                const double arg33 = (zr+v+c+1)*atan(4.0*x1z);  //Z=16735
                const double nen33 = pow(1.0+16*x1z*x1z,(zr+v+c+1)/2.0);  //Z=16736
                const double arg34 = (zr+v+c)*atan(4.0*x1z);  //Z=16737
                const double nen34 = pow(1.0+16*x1z*x1z,(zr+v+c)/2.0);  //Z=16738
                const double arg35 = (zr+v+c-1)*atan(4.0*x1z);  //Z=16739
                const double nen35 = pow(1.0+16*x1z*x1z,(zr+v+c-1)/2.0);  //Z=16740
                const double F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);  //Z=16741
                const double F32as22z = (1/2.0)*ee0*e1*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);  //Z=16742
                const double F32as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);  //Z=16743
                const double F32as24z = (1/2.0)*ee1*e1*(1/(4.0*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);  //Z=16744
                const double F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);  //Z=16745
                const double F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);  //Z=16746
                //F32asz = F32as1z+F32as2z;  //Z=16747
                F32 = F32as1z0+F32as2z;  //Z=16748
                //F32 = F32asz0;  //Z=16749
            }/*4*/  //Z=16750


            /* ** term #4 asymptote ** */  //Z=16753
            if ( xrad>=lim4 )
            {/*4*/  //Z=16754
                const double F42as10z = preg4*preg4*pow(x22z,-2*a1);  //Z=16755
                //F42as1sumz = pva0;  //Z=16756
                //F42as1z = F42as10z*F42as1sumz;  //Z=16757
                const double F42as1z0 = F42as10z*pza;  //Z=16758

                const double arg41 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=16760
                const double nen41 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=16761
                const double arg42 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=16762
                const double nen42 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=16763
                //arg43 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=16764
                //nen43 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=16765
                const double F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=16766
                const double F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=16767
                const double F42as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);  //Z=16768
                //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);  //Z=16769
                //F42as2z = F42as20z*(F42as21+F42as22+F42as23);  //Z=16770
                const double F42as2z0 = F42as20z*(F42as21+F42as22);  //Z=16771

                const double F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=16773
                const double F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=16774
                //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);  //Z=16775
                const double F42as26 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);  //Z=16776
                //F42as3z = F42as30z*(F42as24+F42as25+F42as26);  //Z=16777
                const double F42as3z0 = F42as30z*(F42as24+F42as26);  //Z=16778

                const double F42as40z = preg3*preg3*pow(x2z*x2z,c);  //Z=16780
                const double arg44 = (zr+2*c+1)*atan(4.0*x2z);  //Z=16781
                const double nen44 = pow(1.0+16*x2z*x2z,(zr+2*c+1)/2.0);  //Z=16782
                const double arg45 = (zr+2*c)*atan(4.0*x2z);  //Z=16783
                const double nen45 = pow(1.0+16*x2z*x2z,(zr+2*c)/2.0);  //Z=16784
                const double F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);  //Z=16785
                const double F42as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=16786
                const double F42as29 = (1/2.0)*e1*e0*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=16787
                const double F42as4z = F42as40z*(F42as27+F42as28+F42as29);  //Z=16788
                //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;  //Z=16789
                F42 = F42as1z0+F42as2z0+F42as3z0+F42as4z;  //Z=16790
                //F42 = F42asz0;  //Z=16791
            }/*4*/  //Z=16792


            /* ** term #5 asymptote ** */  //Z=16795
            if ( xradp>=lim5 )
            {/*4*/  //Z=16796
                const double F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);  //Z=16797
                //F52as1sumz = pva0;  //Z=16798
                //F52as1z = F52as10z*F52as1sumz;  //Z=16799
                const double F52as1z0 = F52as10z*pza;  //Z=16800

                const double F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);  //Z=16802
                const double arg51 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=16803
                const double nen51 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=16804
                const double arg52 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=16805
                const double nen52 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=16806
                //arg53 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=16807
                //nen53 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=16808
                const double F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);  //Z=16809
                const double F52as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);  //Z=16810
                //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);  //Z=16811
                //F52as2z = F52as20z*(F52as21+F52as22+F52as23);  //Z=16812
                const double F52as2z0 = F52as20z*(F52as21+F52as22);  //Z=16813

                const double F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);  //Z=16815
                const double arg54 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=16816
                const double nen54 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=16817
                //arg55 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=16818
                //nen55 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=16819
                const double arg56 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=16820
                const double nen56 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=16821
                const double F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);  //Z=16822
                //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);  //Z=16823
                const double F52as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);  //Z=16824
                //F52as3z = F52as30z*(F52as24+F52as25+F52as26);  //Z=16825
                const double F52as3z0 = F52as30z*(F52as24+F52as26);  //Z=16826

                const double F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);  //Z=16828
                const double arg57 = (zr+2*c+1)*atan(2.0*(x1z-x2z));  //Z=16829
                const double nen57 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);  //Z=16830
                const double arg58 = (zr+2*c+1)*atan(2.0*(x1z+x2z));  //Z=16831
                const double nen58 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);  //Z=16832
                const double arg59 = (zr+2*c)*atan(2.0*(x1z-x2z));  //Z=16833
                const double nen59 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);  //Z=16834
                const double arg510 = (zr+2*c)*atan(2.0*(x1z+x2z));  //Z=16835
                const double nen510 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);  //Z=16836
                const double F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);  //Z=16837
                const double F52as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=16838
                const double F52as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=16839
                const double F52as4z = F52as40z*(F52as27+F52as28+F52as29);  //Z=16840
                //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;  //Z=16841
                F52 = F52as1z0+F52as2z0+F52as3z0+F52as4z;  //Z=16842
                //F52 = F52asz0;  //Z=16843
            }/*4*/  //Z=16844

            /* ** term #6 asymptote ** */  //Z=16846
            if ( xradp>=lim6 )
            {/*4*/  //Z=16847
                const double F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);  //Z=16848
                //F62as1sumz = pva0;  //Z=16849
                //F62as1z = F62as10z*F62as1sumz;  //Z=16850
                const double F62as1z0 = F62as10z*pza;  //Z=16851

                const double F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=16853
                const double arg61 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=16854
                const double nen61 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=16855
                const double arg62 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=16856
                const double nen62 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=16857
                //arg63 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=16858
                //nen63 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=16859
                const double F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=16860
                const double F62as22 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=16861
                //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);  //Z=16862
                //F62as2z = F62as20z*(F62as21+F62as22+F62as23);  //Z=16863
                const double F62as2z0 = F62as20z*(F62as21+F62as22);  //Z=16864

                const double F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=16866
                const double F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=16867
                //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);  //Z=16868
                const double F62as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=16869
                //F62as3z = F62as30z*(F62as24+F62as25+F62as26);  //Z=16870
                const double F62as3z0 = F62as30z*(F62as24+F62as26);  //Z=16871

                const double F62as40z = preg3*preg3*pow(x1z*x1z,c);  //Z=16873
                const double arg64 = (zr+2*c+1)*atan(4.0*x1z);  //Z=16874
                const double nen64 = pow(1.0+16*x1z*x1z,(zr+2*c+1)/2.0);  //Z=16875
                const double arg65 = (zr+2*c)*atan(4.0*x1z);  //Z=16876
                const double nen65 = pow(1.0+16*x1z*x1z,(zr+2*c)/2.0);  //Z=16877
                const double F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);  //Z=16878
                const double F62as28 = (1/2.0)*e0*e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=16879
                const double F62as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=16880
                const double F62as4z = F62as40z*(F62as27+F62as28+F62as29);  //Z=16881
                //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;  //Z=16882
                F62 = F62as1z0+F62as2z0+F62as3z0+F62as4z;  //Z=16883
                //F62 = F62asz0;  //Z=16884
            }/*4*/  //Z=16885

            /*formpq:=*/ return pql*(cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;  //Z=16887

        }/*3*/ /*  of inhomogeneous core/shell  */  //Z=16893

        /*  myelin cylinder  */  //Z=16896
        if ( (params.cs==3) || (params.cs==4) )
        {/*3*/  //Z=16897

            /*  cylinder parameters  */  //Z=16899
            const double v = -3/2.0;  //Z=16900
            const double e0 = 1;  //Z=16901
            const double e1 = -9/16.0;  //Z=16902
            const double preg1 = 1/sqrt(M_PI);  //Z=16903
            const double pz2v = 1/(zr*(zr-1)*(zr-2));  //Z=16904
            const double pz2v1 = pz2v/(zr-3);  //Z=16905
            const double pz2v2 = pz2v1/(zr-4);  //Z=16906
            const double lim = 18*exp(-5*params.sigma);  //Z=16907
            const double lim1 = lim*1.2;  //Z=16908
            const double rad = params.CR->myarray[1];  //Z=16909
            const int    inmax = round(params.CR->myarray[14]);  //Z=16910
            const double vvm = params.CR->myarray[15];  //Z=16911
            const double rmax = params.CR->myarray[16];  //Z=16912
            const double xmax = q*rmax;  //Z=16913

            double F12;

            if ( xmax<(lim1) )
            {/*4*/  //Z=16915
                F12 = 0.0;  //Z=16923
                for ( int ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=16924
                    for ( int jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=16925
                        double F12sez = 1.0;  //Z=16926
                        double oldF12sez = 1.0;  //Z=16927
                        double qqnn = 1.0;
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=16928
                            qqnn = qqnn * q*q;
                            double pqsum = 0;  //Z=16929
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=16930
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]);  //Z=16931 */
                                pqsum += pow(params.CR->carr7p[ii],2*mser)*pow(params.CR->carr7p[jj],2*(nser-mser))/(params.CR->carr6p[mser]*params.CR->carr6p[nser-mser]);  //Z=16932

                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=16934 */
                                /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx];  //Z=16935 */
                            }/*8*/  //Z=16936
                            F12sez += params.CR->carr4p[nser]*qqnn*pqsum;  //Z=16937
                            const double delser = fabs((F12sez-oldF12sez)/F12sez);  //Z=16938
                            if ( delser<0.0001 ) break; /* goto 101; */  //Z=16939
                            oldF12sez = F12sez;  //Z=16940
                        }/*7*/  //Z=16941
                        /*101:*/  //Z=16942
                        F12 += params.CR->carr5p[ii]*params.CR->carr5p[jj]*F12sez;  //Z=16943
                    }/*6*/  //Z=16944
                }/*5*/  //Z=16945
                F12 = F12/vvm;  //Z=16946
            }/*4*/  //Z=16948
            else
            {/*4*/  //Z=16949
                const double xrz = q*rad/(zr+1);  //Z=16950
                const double arg = (zr+2*v+1)*atan(2.0*xrz);  //Z=16951
                const double nen = pow(1.0+4*xrz*xrz,(zr+2*v+1)/2.0);  //Z=16952
                const double arg1 = (zr+2*v)*atan(2.0*xrz);  //Z=16953
                const double nen1 = pow(1.0+4*xrz*xrz,(zr+2*v)/2.0);  //Z=16954
                const double arg2 = (zr+2*v-1)*atan(2.0*xrz);  //Z=16955
                const double nen2 = pow(1.0+4*xrz*xrz,(zr+2*v-1)/2.0);  //Z=16956

                F12 = 0.0;  //Z=16958
                for ( int ii=1; ii<=inmax; ii++ )
                {/*5*/  //Z=16959
                    const double a1m = params.CR->carr5p[ii]*pow(params.CR->carr7p[ii],v);   /*  carr7p[ii]:=pp[ii];  //Z=16960 */
                    for ( int jj=1; jj<=inmax; jj++ )
                    {/*6*/  //Z=16961
                        const double a2m = params.CR->carr5p[jj]*pow(params.CR->carr7p[jj],v);  //Z=16962
                        const double xijm = (params.CR->carr3p[ii]-params.CR->carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii];  //Z=16963 */
                        const double arglmz = (zr+1)*atan(xijm);  //Z=16964
                        const double nenlmz = pow(1.0+xijm*xijm,(zr+1)/2.0);  //Z=16965
                        const double xijp = (params.CR->carr3p[ii]+params.CR->carr3p[jj])*q/(zr+1);  //Z=16966
                        const double arglpz = (zr+1)*atan(xijp);  //Z=16967
                        const double nenlpz = pow(1.0+xijp*xijp,(zr+1)/2.0);  //Z=16968
                        /* F12as1z:=e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(pi*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(pi*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=16969 */
                        const double F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(0-(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=16970
                        /* F12as2z:=e0*e1*(1/(carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(pi*(2*v-1)/2)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(pi*(2*v-1)/2)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=16971 */
                        const double F12as2z = e0*e1*(1/(params.CR->carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(1*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-0)/(nen1*nenlpz));  //Z=16972
                        /* F12as3z:=e1*e0*(1/(carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(pi*(2*v-1)/2)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(pi*(2*v-1)/2)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=16973 */
                        const double F12as3z = e1*e0*(1/(params.CR->carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(1*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-0)/(nen1*nenlpz));  //Z=16974
                        /* F12as4z:=e1*e1*(1/(carr7p[ii]*carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(pi*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(pi*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=16975 */
                        const double F12as4z = e1*e1*(1/(params.CR->carr7p[ii]*params.CR->carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(0+1*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=16976

                        F12 += a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);  //Z=16978
                    }/*6*/  //Z=16979
                }/*5*/  //Z=16980
                F12 = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12/vvm;  //Z=16981
                //F12 = F12asy;  //Z=16982
            }/*4*/  //Z=16983
            /*formpq:=*/ return pql*F12;  //Z=16984

        }/*3*/ /*  of myelin  */  //Z=16990

//    }/*2*/ /*  of cylinder  */  //Z=16992

        return 0.0;
}/*1*/  //Z=18387


#endif // SC_LIB_FORMPQ_partCylinder_H
