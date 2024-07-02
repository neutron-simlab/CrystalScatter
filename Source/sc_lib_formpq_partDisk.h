#ifndef SC_LIB_FORMPQ_partDisk_H
#define SC_LIB_FORMPQ_partDisk_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partDisk(double limql, double qx, double qy, double qxs, double qys, double q, int ordis ) const   /*Z=14910*/
{  //Z=15188
    // double sigmal ist nicht zu ersetzen, da es an einer Stelle CALC.epsilon ist, sonst nur CALC.params.sigmal
    // int ordis ist nicht zu ersetzen, da es einmal eine feste Zahl ist, sonst nur CALC.ordis

    double pqr, pql;  //Z=15195

    const double zl = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=15231
    const double zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232
    const double radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=15233

    // Noch fehlende (globale) Variablen und sonstige Anpassungen:
    // An machen Stellen muss "params." eingefügt werden.
    // Die Konstante 'eps' muss durch 'eps9' ersetzt werden.
    const double zz = zr;
    const double z  = zl;
    const double qz = 1; // wird für qrombdeltac() verwendet

    CHECKENDTHREAD_VAL;

    /* ****** */  //Z=16996
    /*  disk  */  //Z=16997
    /* ****** */  //Z=16998
    //if ( params.part==2 )
    //{/*2*/  //Z=16999

    /* ** longitudinal part ** */  //Z=17001
    /* ** isotropic ** */  //Z=17002
    if ( ordis==7 )
    {/*3*/  //Z=17003
        if ( q<(0.5*params.limq1) )
        {/*4*/  //Z=17004
            pql = 1.0;  //Z=17005
            double oldpqsum = 0.0;  //Z=17006
            double qqnn = 1.0;  //Z=17007
            for ( int nser=1; nser<=80; nser++ )
            {/*5*/  //Z=17008
                qqnn = qqnn*q*q;  //Z=17009
                pql += params.CR->carr1p[nser]*qqnn;  //Z=17010
                const double delser = fabs((pql-oldpqsum)/pql);  //Z=17011
                if ( delser<0.0001 ) break; /* goto 70; */  //Z=17012
                oldpqsum = pql;  //Z=17013
            }/*5*/  //Z=17014
            /*70:*/  //Z=17015
            //pql = pqsum;  //Z=17016
        }/*4*/  //Z=17017
        else
        {/*4*/  //Z=17018
            const double arglq = q*params.length/(zl+1);  //Z=17019
            pql = (2/(zl*(zl-1)))*pow(arglq,-2);  //Z=17020
        }/*4*/  //Z=17021
    }/*3*/  /*  of isotropic  */  //Z=17022

    /*  perfect  */  //Z=17024
    if ( ordis==6 )
    {/*3*/  //Z=17025
        if ( (limql*params.length) < 5.0 ) //20240301 - war if ( limql<0.7*params.limq1 )
        {/*4*/  //Z=17026

            pql = 1.0;  //Z=17028
            double oldpqsum = 0.0;  //Z=17029
            double qqnn = 1.0;  //Z=17030
            if ( params.orcase==1 )
            {/*5*/  //Z=17031
                const double argq = qxs+qys;  //Z=17032
                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17033
                    qqnn = qqnn*q*q;  //Z=17034
                    pql += params.CR->carr1p[nser]*qqnn*pow(1.0-argq*argq,nser);  //Z=17035
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17036
                    if ( delser<0.0001 ) break; /* goto 76; */  //Z=17037
                    oldpqsum = pql;  //Z=17038
                }/*6*/  //Z=17039
            }/*5*/  //Z=17040
            if ( params.orcase==2 )
            {/*5*/  //Z=17041
                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17042
                    qqnn = qqnn*q*q;  //Z=17043
                    pql += params.CR->carr1p[nser]*qqnn*pow(1.0-qxs*qxs,nser);  //Z=17044
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17045
                    if ( delser<0.0001 ) break; /* goto 76; */  //Z=17046
                    oldpqsum = pql;  //Z=17047
                }/*6*/  //Z=17048
            }/*5*/  //Z=17049
            if ( params.orcase==3 )
            {/*5*/  //Z=17050
                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17051
                    qqnn = qqnn*q*q;  //Z=17052
                    pql += params.CR->carr1p[nser]*qqnn*pow(1.0-qys*qys,nser);  //Z=17053
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17054
                    if ( delser<0.0001 ) break; /* goto 76; */  //Z=17055
                    oldpqsum = pql;  //Z=17056
                }/*6*/  //Z=17057
            }/*5*/  //Z=17058
            if ( params.orcase==4 )
            {/*5*/  //Z=17059
                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17060
                    qqnn = qqnn*q*q;  //Z=17061
                    pql += params.CR->carr1p[nser]*qqnn;  //Z=17062
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17063
                    if ( delser<0.0001 ) break; /* goto 76; */  //Z=17064
                    oldpqsum = pql;  //Z=17065
                }/*6*/  //Z=17066
            }/*5*/  //Z=17067
            /*76:*/  //Z=17068
            //pql = pqsum;  //Z=17069
        }/*4*/  //Z=17070
        else
        {/*4*/  //Z=17071
            double arglq;
            if ( params.orcase==1 )
            {/*5*/  //Z=17072
                const double qnarg = qxs+qys;  //Z=17073
                arglq = sqrt(1.0-qnarg*qnarg)*q*params.length/(zl+1)+eps9;  //Z=17074
            }/*5*/  //Z=17075
            if ( params.orcase==2 ) arglq = sqrt(1.0-qxs*qxs)*q*params.length/(zl+1)+eps9;  //Z=17076
            if ( params.orcase==3 ) arglq = sqrt(1.0-qys*qys)*q*params.length/(zl+1)+eps9;  //Z=17077
            if ( params.orcase==4 ) arglq = q*params.length/(zl+1)+eps9;  //Z=17078

            const double pqr1 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3);  //Z=17080
            const double pqr2 = (1/(zl*(zl-1)*(zl-2)))*pow(arglq,-3)*sin((zl-2)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-2)/2.0);  //Z=17081
            const double pqr3 = (1/(zl*(zl-1)*(zl-2)*(zl-3)))*pow(arglq,-4)*cos((zl-3)*atan(2.0*arglq))/pow(1.0+4*arglq*arglq,(zl-3)/2.0);  //Z=17082
            pql = (4/M_PI)*(pqr1-pqr2-(9/8.0)*pqr3);  //Z=17083
        }/*4*/  //Z=17084
    }/*3*/   /*  of perfect  */  //Z=17085

    /*  orientational distribution  */  //Z=17087
    if ( ordis==0 )
    {/*3*/       /*  general  */  //Z=17088
        double qxn[121], qyn[121];
        if ( params.orcase==1 )
        {/*4*/  //Z=17089
            if ( q < (6.0*params.limq1) ) //20240301 - war if ( limql<(0.3*params.limq1) )
            {/*5*/  //Z=17090
                pql = 1.0;  //Z=17091
                double oldpqsum = 0.0;  //Z=17092
                double qqnn = 1.0;  //Z=17093
                qxn[0] = 1.0;  //Z=17094
                qyn[0] = 1.0;  //Z=17095

                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17097
                    qqnn = qqnn*q*q;  //Z=17098
                    qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17099
                    qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17100

                    double binsum = 0.0;  //Z=17102
                    for ( int mser=0; mser<=nser; mser++ )
                    {/*7*/  //Z=17103
                        double binsum1 = 0.0;  //Z=17104
                        for ( int lser=0; lser<=mser; lser++ )
                        {/*8*/  //Z=17105
                            /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=17106 */
                            /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=17107 */
                            binsum1 += params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=17108
                        }/*8*/  //Z=17109
                        /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17110 */
                        /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=17111 */
                        binsum += params.CR->carr11pm[nser][mser]*binsum1;  //Z=17112
                    }/*7*/  //Z=17113
                    pql += params.CR->carr1p[nser]*qqnn*binsum;  //Z=17114
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17115
                    if ( delser<0.0001 ) break; /* goto 77; */  //Z=17116
                    oldpqsum = pql;  //Z=17117
                }/*6*/  //Z=17118
                /*77:*/  //Z=17119
                //pql = pqsum;  //Z=17120
            }/*5*/  //Z=17121
            else
            {/*5*/  //Z=17122
                /*  disk: length = disk radius  */  //Z=17123
                /*  always use Bessel function approximation  */  //Z=17124
                qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17125
                pql = pql/params.norm;  //Z=17126
            }/*5*/  //Z=17127
        }/*4*/  //Z=17128

        if ( params.orcase==2 )
        {/*4*/   /*  x-axis  */  //Z=17130
            if ( q < (6000.0*params.limq1) ) //20240301 - war if ( limql<(0.9*params.limq1) )
            {/*5*/  //Z=17131
                pql = 1.0;  //Z=17132
                double oldpqsum = 0.0;  //Z=17133
                double qqnn = 1.0;  //Z=17134
                qxn[0] = 1.0;  //Z=17135
                qyn[0] = 1.0;  //Z=17136

                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17138
                    qqnn = qqnn*q*q;  //Z=17139
                    qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17140
                    qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17141

                    double binsum = 0.0;  //Z=17143
                    for ( int mser=0; mser<=nser; mser++ )
                    {/*7*/  //Z=17144
                        double binsum1 = 0.0;  //Z=17145
                        for ( int lser=0; lser<=mser; lser++ )
                        {/*8*/  //Z=17146
                            /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=17147 */
                            /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=17148 */
                            binsum1 += params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=17149
                        }/*8*/  //Z=17150
                        /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17151 */
                        /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=17152 */
                        binsum += params.CR->carr11pm[nser][mser]*binsum1;  //Z=17153
                    }/*7*/  //Z=17154
                    pql += params.CR->carr1p[nser]*qqnn*binsum;  //Z=17155
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17156
                    if ( delser<0.0001 ) break; /* goto 78; */  //Z=17157
                    oldpqsum = pql;  //Z=17158
                }/*6*/  //Z=17159
                /*78:*/  //Z=17160
                //pql = pqsum;  //Z=17161
            }/*5*/  //Z=17162
            else
            {/*5*/  //Z=17163
                /*  disk: length = disk radius  */  //Z=17164
                /*  always use Bessel function approximation  */  //Z=17165
                qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17166
                pql = pql/params.norm;  //Z=17167
                /* pql:=pq/1e-5;  //Z=17168 */
                /* pql:=0.5;  //Z=17169 */
            }/*5*/  //Z=17170
        }/*4*/  //Z=17171

        if ( params.orcase==3 )
        {/*4*/     /*  y-axis  */  //Z=17173
            if ( q < (6000.0*params.limq1) ) //20240301 - war if ( limql<(0.9*params.limq1) )
            {/*5*/  //Z=17174
                pql = 1.0;  //Z=17175
                double oldpqsum = 0.0;  //Z=17176
                double qqnn = 1.0;  //Z=17177
                qxn[0] = 1.0;  //Z=17178
                qyn[0] = 1.0;  //Z=17179

                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17181
                    qqnn = qqnn*q*q;  //Z=17182
                    qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17183
                    qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17184

                    double binsum = 0.0;  //Z=17186
                    for ( int mser=0; mser<=nser; mser++ )
                    {/*7*/  //Z=17187
                        double binsum1 = 0.0;  //Z=17188
                        for ( int lser=0; lser<=mser; lser++ )
                        {/*8*/  //Z=17189
                            /* indx:=lser+1+round(mser*(mser+1)/2);  //Z=17190 */
                            /* binsum1:=binsum1+carr2pm[indx]*qxn[lser]*qyn[mser-lser];  //Z=17191 */
                            binsum1 += params.CR->carr22pm[mser][lser]*qxn[lser]*qyn[mser-lser];  //Z=17192
                        }/*8*/  //Z=17193
                        /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17194 */
                        /* binsum:=binsum+carr1pm[indx]*binsum1;  //Z=17195 */
                        binsum += params.CR->carr11pm[nser][mser]*binsum1;  //Z=17196
                    }/*7*/  //Z=17197
                    pql += params.CR->carr1p[nser]*qqnn*binsum;  //Z=17198
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17199
                    if ( delser<0.0001 ) break; /* goto 79; */  //Z=17200
                    oldpqsum = pql;  //Z=17201
                }/*6*/  //Z=17202
                /*79:*/  //Z=17203
                //pql = pqsum;  //Z=17204
            }/*5*/  //Z=17205
            else
            {/*5*/  //Z=17206
                /*  disk: length = disk radius  */  //Z=17207
                /*  always use Bessel function approximation  */  //Z=17208
                qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17209
                pql = pql/params.norm;  //Z=17210
            }/*5*/  //Z=17211
        }/*4*/  //Z=17212

        if ( params.orcase==4 )
        {/*4*/  //Z=17214
            if ( q < (6000.0*params.limq1) ) //20240301 - war if ( limql<(0.7*params.limq1) )
            {/*5*/  //Z=17215
                pql = 1.0;  //Z=17216
                double oldpqsum = 0.0;  //Z=17217
                double qqnn = 1.0;  //Z=17218
                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17219
                    qqnn = qqnn*q*q;  //Z=17220
                    pql += params.CR->carr1p[nser]*qqnn;  //Z=17221
                    const double delser = fabs((pql-oldpqsum)/pql);  //Z=17222
                    if ( delser<0.0001 ) break; /* goto 80; */  //Z=17223
                    oldpqsum = pql;  //Z=17224
                }/*6*/  //Z=17225
                /*80:*/  //Z=17226
                //pql = pqsum;  //Z=17227
            }/*5*/  //Z=17228
            else
            {/*5*/  //Z=17229
                /*  disk: length = disk radius  */  //Z=17230
                /*  always use Bessel function approximation  */  //Z=17231
                qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,params.polPhi,qx,qy,qz,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/ordis,2,6,params.orcase+7,0,0,0,params.CR->carr2p,pql);  //Z=17232
                pql = pql/params.norm;  //Z=17233
            }/*5*/  //Z=17234
        }/*4*/  //Z=17235
    }/*3*/   /*  of orientational distribution  */  //Z=17236


    /*  transverse part  */  //Z=17239
    /*  disk: radius = disk thickness/2  */  //Z=17240
    /*  homogeneous disk  */  //Z=17241
    if ( params.cs==0 )
    {/*3*/  //Z=17242
        if ( q<(0.5*params.limq4) ) //20240301 - war 0.3
        {/*4*/  //Z=17243
            pqr = 1.0;  //Z=17244
            double oldpqsum = 0.0;  //Z=17245
            double qqnn = 1.0;  //Z=17246
            for ( int nser=1; nser<=100; nser++ )
            {/*5*/  //Z=17247
                qqnn = qqnn*q*q;  //Z=17248
                pqr += params.CR->carr4p[nser]*qqnn;  //Z=17249
                const double delser = fabs((pqr-oldpqsum)/pqr);  //Z=17250
                if ( delser<0.0001 ) break; /* goto 71; */  //Z=17251
                oldpqsum = pqr;  //Z=17252
            }/*5*/  //Z=17253
            /*71:*/  //Z=17254
            //pqr = pqsum;  //Z=17255
        }/*4*/  //Z=17256
        else
        {/*4*/  //Z=17257
            const double argpq = q*params.radius/(zr+1);  //Z=17258
            pqr = (1/(2.0*zr*(zr-1)))*pow(argpq,-2)*(1-cos((zr-1)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-1)/2.0));  //Z=17259
        }/*4*/  //Z=17260
        /*formpq:=*/ return pql*pqr;  //Z=17261
        /* formpq:=pql;  //Z=17262 */
    }/*3*/ /*  of homogeneous  */  //Z=17263

    /*  core/shell disk  */  //Z=17265
    if ( params.cs==1 )
    {/*3*/  //Z=17266
        const double ccc1 = sqr(1-params.rho)*pow(params.p1,2);  //Z=17267
        const double ccc2 = 2*params.rho*(1-params.rho)*pow(params.p1,1);  //Z=17268
        const double ccc3 = params.rho*params.rho;  //Z=17269
        const double vv3 = sqr((1-params.rho)*pow(params.p1,1)+params.rho);  //Z=17270

        const double argq = q*radiusm/(zz+1);  //Z=17272
        const double argpq = q*params.radius/(zz+1);  //Z=17273

        double F121, F122, F123;

        /*  F121 disk  */  //Z=17275
        if ( q<(0.8*params.limq4) )
        {/*4*/  //Z=17276
            /* ** series expansion ** */  //Z=17277
            double pqsum = 1.0;  //Z=17278
            double oldpqsum = 0.0;  //Z=17279
            double qqnn = 1.0;  //Z=17280
            for ( int nser=1; nser<=120; nser++ )
            {/*5*/  //Z=17281
                qqnn = qqnn*q*q;  //Z=17282
                pqsum += params.CR->carr4p[nser]*qqnn;  //Z=17283
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17284
                if ( delser<0.0001 ) break; /* goto 72; */  //Z=17285
                oldpqsum = pqsum;  //Z=17286
            }/*5*/  //Z=17287
            /*72:*/  //Z=17288
            F121 = ccc1*pqsum/vv3;  //Z=17289
        }/*4*/  //Z=17290
        else
        {/*4*/  //Z=17291
            const double pqr = (1/(2.0*zr*(zr-1)))*pow(argpq,-2)*(1-cos((zr-1)*atan(2.0*argpq))/pow(1.0+4*argpq*argpq,(zr-1)/2.0));  //Z=17292
            F121 = ccc1*pqr/vv3;  //Z=17293
        }/*4*/  //Z=17294

        /*  F122 disk  */  //Z=17296
        if ( q<(2.0*params.limq5) )
        {/*4*/  //Z=17297
            /* ** series expansion ** */  //Z=17298
            double pqsum = 1.0;  //Z=17299
            double oldpqsum = 0.0;  //Z=17300
            double qqnn = 1.0;  //Z=17301
            for ( int nser=1; nser<=120; nser++ )
            {/*5*/  //Z=17302
                qqnn = qqnn*q*q;  //Z=17303
                pqsum += params.CR->carr5p[nser]*qqnn;  //Z=17304
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17305
                if ( delser<0.0001 ) break; /* goto 73; */  //Z=17306
                oldpqsum = pqsum;  //Z=17307
            }/*5*/  //Z=17308
            /*73:*/  //Z=17309
            F122 = ccc2*pqsum/vv3;  //Z=17310
        }/*4*/  //Z=17311
        else
        {/*4*/  //Z=17312
            const double argbm = (zr+1)*atan(argpq-argq);  //Z=17313
            const double nenbm = pow(1.0+4*sqr(argpq-argq),(zr+1)/2.0);  //Z=17314
            const double argbp = (zr+1)*atan(argpq+argq);  //Z=17315
            const double nenbp = pow(1.0+4*sqr(argpq+argq),(zr+1)/2.0);  //Z=17316

            const double pqr = (1/(2.0*zr*(zr-1)*(zr-2)*argpq*argq))*(cos(argbm)/nenbm-cos(argbp)/nenbp);  //Z=17318
            F122 = ccc2*pqr/vv3;  //Z=17319
        }/*4*/  //Z=17320

        /*  F123 disk  */  //Z=17322
        if ( q<(0.3*params.limq6) )
        {/*4*/  //Z=17323
            /* ** series expansion ** */  //Z=17324
            double pqsum = 1.0;  //Z=17325
            double oldpqsum = 0.0;  //Z=17326
            double qqnn = 1.0;  //Z=17327
            for ( int nser=1; nser<=120; nser++ )
            {/*5*/  //Z=17328
                qqnn = qqnn*q*q;  //Z=17329
                pqsum += params.CR->carr6p[nser]*qqnn;  //Z=17330
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17331
                if ( delser<0.0001 ) break; /* goto 74; */  //Z=17332
                oldpqsum = pqsum;  //Z=17333
            }/*5*/  //Z=17334
            /*74:*/  //Z=17335
            F123 = ccc3*pqsum/vv3;  //Z=17336
        }/*4*/  //Z=17337
        else
        {/*4*/  //Z=17338
            const double pqr = (1/(2.0*zr*(zr-1)))*pow(argq,-2)*(1-cos((zr-1)*atan(2.0*argq))/pow(1.0+4*argq*argq,(zr-1)/2.0));  //Z=17339
            F123 = ccc3*pqr/vv3;  //Z=17340
            /*  add more terms, if necessary  */  //Z=17341
        }/*4*/  //Z=17342
        /*formpq:=*/ return pql*(F121+F122+F123);  //Z=17343
        /* formpq:=F122;  //Z=17344 */
    }/*3*/ /*  of core/shell-disk  */  //Z=17345

    /*  inhomogeneous core/shell disk  */  //Z=17347
    if ( params.cs==2 )
    {/*3*/  //Z=17348

        const double dim = 1;  //Z=17350
        const double delc = 0.0001;  //Z=17351
        const double xrad = q*radiusm;  //Z=17352
        const double xradp = q*params.radius;  //Z=17353
        const double x1z = q*params.radius/(2.0*(zr+1));  //Z=17354
        const double x12z = x1z*x1z;  //Z=17355
        const double x2z = q*radiusm/(2.0*(zr+1));  //Z=17356
        const double x22z = x2z*x2z;  //Z=17357

        const double lim = 18*exp(-5*params.sigma);  //Z=17359
        const double lim1 = lim;  //Z=17360
        const double lim2 = lim*0.7;  //Z=17361
        const double lim3 = lim;  //Z=17362
        const double lim4 = lim;  //Z=17363
        const double lim5 = lim*0.7;  //Z=17364
        const double lim6 = lim*1.2;  //Z=17365

        const double a1 = (dim-params.alphash1)/2.0;  //Z=17367
        const double b1 = dim/2.0;  //Z=17368
        const double b2 = (dim+2-params.alphash1)/2.0;  //Z=17369
        const double b1s = (dim+2)/2.0;  //Z=17370
        const double v = -b1s+1/2.0;  //Z=17371
        const double c = a1-b1-b2+1/2.0;  //Z=17372
        const double d0 = 1;  //Z=17373
        //d1 = a1*(1+a1-b1)*(1+a1-b2);  //Z=17374
        const double e0 = 1.0;  //Z=17375
        const double e1 = (3/8.0)-(b1+b2)+((b1-b2)*(b1-b2)-3*a1*a1+2*a1*(1+b1+b2))/2.0;  //Z=17376
        const double ee0 = 1.0;  //Z=17377
        const double ee1 = 3*(3-8*b1s+4*b1s*b1s)/(16.0*(1-b1s));  //Z=17378

        const double gb1s = sqrt(M_PI)/2.0;  //Z=17380
        const double pz2v = 1/(zr*(zr-1));  //Z=17381
        const double pz2v1 = pz2v/(zr-2);  //Z=17382
        const double pz2v2 = pz2v1/(zr-3);  //Z=17383

        const double gz1 = gamma(zr+1);  //Z=17385
        const double preg1 = gb1s/sqrt(M_PI);  //Z=17386
        const double preg3 = gamma(b1)*gamma(b2)/(gamma(a1)*sqrt(M_PI));  //Z=17387
        const double preg4 = gamma(b1)*gamma(b2)/(gamma(b1-a1)*gamma(b2-a1));  //Z=17388
        const double pzvc = gamma(zr+1+v+c)/gz1;  //Z=17389
        const double pzvc1 = gamma(zr+1+v+c-1)/gz1;  //Z=17390
        const double pzvc2 = gamma(zr+1+v+c-2)/gz1;  //Z=17391
        const double pzac = gamma(zr+1-2*a1+c)/gz1;  //Z=17392
        const double pzac1 = gamma(zr+1-2*a1+c-1)/gz1;  //Z=17393
        //pzac2 = gamma(zr+1-2*a1+c+2)/gz1;  //Z=17394
        const double pzc = gamma(zr+1+2*c)/gz1;  //Z=17395
        const double pzc1 = gamma(zr+1+2*c-1)/gz1;  //Z=17396
        const double pza = gamma(zr+1-4*a1)/gz1;  //Z=17397
        const double pzva = gamma(zr+1+v-2*a1)/gz1;  //Z=17398
        const double pzva1 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=17399
        //dnv0 = 1;  //Z=17400
        //pvav0 = gamma(zr+1+v-2*a1)/gz1;  //Z=17401
        //pvav10 = gamma(zr+1+v-2*a1-1)/gz1;  //Z=17402
        //pva0 = gamma(zr+1-4*a1)/gz1;  //Z=17403

        const double cc1 = 1/(dim*dim);  //Z=17405
        const double cc2 = 2*params.rho/(dim*(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=17406
        const double cc3 = -2*params.rho/(dim*(dim-params.alphash1));  //Z=17407
        const double cc4 = sqr(params.rho)/(sqr(dim-params.alphash1)*pow(sqr(params.p1),dim-params.alphash1));  //Z=17408
        const double cc5 = -2*sqr(params.rho)/(sqr(dim-params.alphash1)*pow(params.p1,dim-params.alphash1));  //Z=17409
        const double cc6 = sqr(params.rho)/sqr(dim-params.alphash1);  //Z=17410
        const double vv3 = cc1+cc2+cc3+cc4+cc5+cc6;  //Z=17411

        double F12, F22, F32, F42, F52, F62;

        /*  term #1 series  */  //Z=17413
        if ( (xradp)<lim1 )
        {/*4*/  //Z=17414
            //z12v[0] = 1;  //Z=17415
            //b1sv[0] = 1;  //Z=17416
            //fkv[0] = 1;  //Z=17417
            double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17418
            F12 = 1.0;  //Z=17419
            double oldF12 = 0.0;  //Z=17420
            for ( int n=1; n<=120; n++ )
            {/*5*/  //Z=17421
                //qqn[n] = qqn[n-1]*q*q;  //Z=17422
                qqnn *= sqr(q);
                //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17423
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17424
                //fkv[n] = fkv[n-1]*n;  //Z=17425
                //sum12[n] = 0;  //Z=17426
                /* for m:=0 to n do sum12[n]:=sum12[n]+1/(b1sv[m]*b1sv[n-m]*fkv[m]*fkv[n-m]);  //Z=17427 */
                /* F12sez:=F12sez+power(-x12z,n)*z12v[n]*sum12[n];  //Z=17428 */

                F12 += params.CR->carr4p[n]*qqnn; //qqn[n];  //Z=17430

                const double del = fabs((F12-oldF12)/F12);  //Z=17432
                if ( del<delc ) break; /* goto 221; */  //Z=17433
                oldF12 = F12;  //Z=17434
            }/*5*/  //Z=17435
            /*221:*/  //Z=17436
            //F12 = F12sez;  //Z=17437
        }/*4*/  //Z=17438

        /*  term #2 series  */  //Z=17440
        if ( (xradp)<lim2 )
        {/*4*/  //Z=17441
            //z12v[0] = 1;  //Z=17442
            //a1v[0] = 1;  //Z=17443
            //b1v[0] = 1;  //Z=17444
            //b2v[0] = 1;  //Z=17445
            //b1sv[0] = 1;  //Z=17446
            //fkv[0] = 1;  //Z=17447
            double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17448
            F22 = 1.0;  //Z=17449
            double oldF22 = 0.0;  //Z=17450
            for ( int n=1; n<=120; n++ )
            {/*5*/  //Z=17451
                //qqn[n] = qqn[n-1]*q*q;  //Z=17452
                qqnn *= sqr(q);
                //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17453
                //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17454
                //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17455
                //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17456
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17457
                //fkv[n] = fkv[n-1]*n;  //Z=17458
                //sum22[n] = 0;  //Z=17459
                /* for m:=0 to n do sum22[n]:=sum22[n]+a1v[n-m]*power(p1*p1,m)/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17460 */
                /* F22sez:=F22sez+power(-x22z,n)*z12v[n]*sum22[n];  //Z=17461 */

                F22 += params.CR->carr5p[n]*qqnn; //qqn[n];  //Z=17463

                const double del = fabs((F22-oldF22)/F22);  //Z=17465
                if ( del<delc ) break; /* goto 222; */  //Z=17466
                oldF22 = F22;  //Z=17467
            }/*5*/  //Z=17468
            /*222:*/  //Z=17469
            //F22 = F22sez;  //Z=17470
        }/*4*/  //Z=17471

        /*  term #3 series  */  //Z=17473
        if ( (xradp)<lim3 )
        {/*4*/  //Z=17474
            //z12v[0] = 1;  //Z=17475
            //a1v[0] = 1;  //Z=17476
            //b1v[0] = 1;  //Z=17477
            //b2v[0] = 1;  //Z=17478
            //b1sv[0] = 1;  //Z=17479
            //fkv[0] = 1;  //Z=17480
            double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17481
            F32 = 1.0;  //Z=17482
            double oldF32 = 0.0;  //Z=17483
            for ( int n=1; n<=120; n++ )
            {/*5*/  //Z=17484
                //qqn[n] = qqn[n-1]*q*q;  //Z=17485
                qqnn *= sqr(q);
                //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17486
                //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17487
                //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17488
                //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17489
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17490
                //fkv[n] = fkv[n-1]*n;  //Z=17491
                //sum32[n] = 0;  //Z=17492
                /* for m:=0 to n do sum32[n]:=sum32[n]+a1v[n-m]/(b1sv[m]*b1v[n-m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17493 */
                /* F32sez:=F32sez+power(-x12z,n)*z12v[n]*sum32[n];  //Z=17494 */

                F32 += params.CR->carr6p[n]*qqnn; //qqn[n];  //Z=17496

                const double del = fabs((F32-oldF32)/F32);  //Z=17498
                if ( del<delc ) break; /* goto 223; */  //Z=17499
                oldF32 = F32;  //Z=17500
            }/*5*/  //Z=17501
            /*223:*/  //Z=17502
            //F32 = F32sez;  //Z=17503
        }/*4*/  //Z=17504

        /*  term #4 series  */  //Z=17506
        if ( (xradp)<lim4 )
        {/*4*/  //Z=17507
            //z12v[0] = 1;  //Z=17508
            //a1v[0] = 1;  //Z=17509
            //b1v[0] = 1;  //Z=17510
            //b2v[0] = 1;  //Z=17511
            //b1sv[0] = 1;  //Z=17512
            //fkv[0] = 1;  //Z=17513
            double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17514
            F42 = 1.0;  //Z=17515
            double oldF42 = 0.0;  //Z=17516
            for ( int n=1; n<=120; n++ )
            {/*5*/  //Z=17517
                //qqn[n] = qqn[n-1]*q*q;  //Z=17518
                qqnn *= sqr(q);
                //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17519
                //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17520
                //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17521
                //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17522
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17523
                //fkv[n] = fkv[n-1]*n;  //Z=17524
                //sum42[n] = 0;  //Z=17525
                /* for m:=0 to n do sum42[n]:=sum42[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17526 */
                /* F42sez:=F42sez+power(-x22z,n)*z12v[n]*sum42[n];  //Z=17527 */

                F42 += params.CR->carr7p[n]*qqnn; //qqn[n];  //Z=17529

                const double del = fabs((F42-oldF42)/F42);  //Z=17531
                if ( del<delc ) break; /* goto 224; */  //Z=17532
                oldF42 = F42;  //Z=17533
            }/*5*/  //Z=17534
            /*224:*/  //Z=17535
            //F42 = F42sez;  //Z=17536
        }/*4*/  //Z=17537

        /*  term #5 series  */  //Z=17539
        if ( (xradp)<lim5 )
        {/*4*/  //Z=17540
            //z12v[0] = 1;  //Z=17541
            //a1v[0] = 1;  //Z=17542
            //b1v[0] = 1;  //Z=17543
            //b2v[0] = 1;  //Z=17544
            //b1sv[0] = 1;  //Z=17545
            //fkv[0] = 1;  //Z=17546
            double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17547
            F52 = 1.0;  //Z=17548
            double oldF52 = 0.0;  //Z=17549
            for ( int n=1; n<=120; n++ )
            {/*5*/  //Z=17550
                //qqn[n] = qqn[n-1]*q*q;  //Z=17551
                qqnn *= sqr(q);
                //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17552
                //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17553
                //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17554
                //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17555
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17556
                //fkv[n] = fkv[n-1]*n;  //Z=17557
                //sum52[n] = 0;  //Z=17558
                /* for m:=0 to n do sum52[n]:=sum52[n]+a1v[m]*a1v[n-m]*power(p1*p1,m)/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17559 */
                /* F52sez:=F52sez+power(-x22z,n)*z12v[n]*sum52[n];  //Z=17560 */

                F52 += params.CR->carr8p[n]*qqnn; //qqn[n];  //Z=17562

                const double del = fabs((F52-oldF52)/F52);  //Z=17564
                if ( del<delc ) break; /* goto 225; */  //Z=17565
                oldF52 = F52;  //Z=17566
            }/*5*/  //Z=17567
            /*225:*/  //Z=17568
            //F52 = F52sez;  //Z=17569
        }/*4*/  //Z=17570

        /*  term #6 series  */  //Z=17572
        if ( (xradp)<lim6 )
        {/*4*/  //Z=17573
            //z12v[0] = 1;  //Z=17574
            //a1v[0] = 1;  //Z=17575
            //b1v[0] = 1;  //Z=17576
            //b2v[0] = 1;  //Z=17577
            //b1sv[0] = 1;  //Z=17578
            //fkv[0] = 1;  //Z=17579
            double qqnn = 1.0; //qqn[0] = 1.0;  //Z=17580
            F62 = 1.0;  //Z=17581
            double oldF62 = 0.0;  //Z=17582
            for ( int n=1; n<=120; n++ )
            {/*5*/  //Z=17583
                //qqn[n] = qqn[n-1]*q*q;  //Z=17584
                qqnn *= sqr(q);
                //z12v[n] = z12v[n-1]*((zr+1)-2+2*n)*((zr+1)-1+2*n);  //Z=17585
                //a1v[n] = a1v[n-1]*(a1-1+n);  //Z=17586
                //b1v[n] = b1v[n-1]*(b1-1+n);  //Z=17587
                //b2v[n] = b2v[n-1]*(b2-1+n);  //Z=17588
                //b1sv[n] = b1sv[n-1]*(b1s-1+n);  //Z=17589
                //fkv[n] = fkv[n-1]*n;  //Z=17590
                //sum62[n] = 0;  //Z=17591
                /* for m:=0 to n do sum62[n]:=sum62[n]+a1v[m]*a1v[n-m]/(b1v[m]*b1v[n-m]*b2v[m]*b2v[n-m]*fkv[m]*fkv[n-m]);  //Z=17592 */
                /* F62sez:=F62sez+power(-x12z,n)*z12v[n]*sum62[n];  //Z=17593 */

                F62 += params.CR->carr9p[n]*qqnn; //qqn[n];  //Z=17595

                const double del = fabs((F62-oldF62)/F62);  //Z=17597
                if ( del<delc ) break; /* goto 226; */  //Z=17598
                oldF62 = F62;  //Z=17599
            }/*5*/  //Z=17600
            /*226:*/  //Z=17601
            //F62 = F62sez;  //Z=17602
        }/*4*/  //Z=17603


        /* ** term #1 asymptote ** */  //Z=17606
        if ( xradp>=lim1 )
        {/*4*/  //Z=17607
            const double arg11 = (zr+2*v+1)*atan(4.0*x1z);  //Z=17608
            const double nen11 = pow(1.0+16*x1z*x1z,(zr+2*v+1)/2.0);  //Z=17609
            const double arg12 = (zr+2*v)*atan(4.0*x1z);  //Z=17610
            const double nen12 = pow(1.0+16*x1z*x1z,(zr+2*v)/2.0);  //Z=17611
            const double arg13 = (zr+2*v-1)*atan(4.0*x1z);  //Z=17612
            const double nen13 = pow(1.0+16*x1z*x1z,(zr+2*v-1)/2.0);  //Z=17613

            const double F12as1z = ee0*ee0*pz2v*(1+cos(M_PI*v)*cos(arg11)/nen11-sin(M_PI*v)*sin(arg11)/nen11);  //Z=17615
            const double F12as2z = 2*ee0*ee1*(1/(2.0*x1z))*pz2v1*(cos(M_PI*(2*v-1)/2.0)*cos(arg12)/nen12-sin(M_PI*(2*v-1)/2.0)*sin(arg12)/nen12);  //Z=17616
            const double F12as3z = ee1*ee1*(1/(4.0*x1z*x1z))*pz2v2*(1+cos(M_PI*(v-1))*cos(arg13)/nen13-sin(M_PI*(v-1))*sin(arg13)/nen13);  //Z=17617
            F12 = preg1*preg1*pow(x1z,2*v)*(1/2.0)*(F12as1z+F12as2z+F12as3z);  //Z=17618
            //F12 = F12asz;  //Z=17619
        }/*4*/  //Z=17620

        /* ** term #2 asymptote ** */  //Z=17622
        if ( xradp>=lim2 )
        {/*4*/  //Z=17623
            //arg21 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=17624
            //nen21 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=17625
            //arg22 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=17626
            //nen22 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=17627
            //F22as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg21)/nen21-sin(M_PI*v/2.0)*sin(arg21)/nen21);  //Z=17628
            //F22as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg22)/nen22-sin(M_PI*(v-1)/2.0)*sin(arg22)/nen22);  //Z=17629
            const double F22as10z = preg1*preg4*pow(x1z,v)*pow(x22z,-a1);  //Z=17630
            //F22as1z = F22as10z*(F22as1sum1z+F22as1sum2z);  //Z=17631

            const double arg210 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=17633
            const double nen210 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=17634
            const double arg220 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=17635
            const double nen220 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=17636
            const double F22as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg210)/nen210-sin(M_PI*v/2.0)*sin(arg210)/nen210);  //Z=17637
            const double F22as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg220)/nen220-sin(M_PI*(v-1)/2.0)*sin(arg220)/nen220);  //Z=17638
            const double F22as1z0 = F22as10z*(F22as1sum1z0+F22as1sum2z0);  //Z=17639
            const double arg23 = (zr+v+c+1)*atan(2.0*(x1z-x2z));  //Z=17640
            const double nen23 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c+1)/2.0);  //Z=17641
            const double arg24 = (zr+v+c+1)*atan(2.0*(x1z+x2z));  //Z=17642
            const double nen24 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c+1)/2.0);  //Z=17643
            const double arg25 = (zr+v+c)*atan(2.0*(x1z-x2z));  //Z=17644
            const double nen25 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c)/2.0);  //Z=17645
            const double arg26 = (zr+v+c)*atan(2.0*(x1z+x2z));  //Z=17646
            const double nen26 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c)/2.0);  //Z=17647
            const double arg27 = (zr+v+c-1)*atan(2.0*(x1z-x2z));  //Z=17648
            const double nen27 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+v+c-1)/2.0);  //Z=17649
            const double arg28 = (zr+v+c-1)*atan(2.0*(x1z+x2z));  //Z=17650
            const double nen28 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+v+c-1)/2.0);  //Z=17651

            const double a22as21z = (1/2.0)*ee0*e0*pzvc;  //Z=17653
            const double F22as21z = a22as21z*(cos(M_PI*(v-c)/2.0)*cos(arg23)/nen23-sin(M_PI*(v-c)/2.0)*sin(arg23)/nen23+cos(M_PI*(v+c)/2.0)*cos(arg24)/nen24-sin(M_PI*(v+c)/2.0)*sin(arg24)/nen24);  //Z=17654
            const double a22as22z = (1/2.0)*ee0*e1*(1/(2.0*x2z))*pzvc1;  //Z=17655
            const double F22as22z = a22as22z*(cos(M_PI*(v-c+1)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-c+1)/2.0)*sin(arg25)/nen25+cos(M_PI*(v+c-1)/2.0)*cos(arg26)/nen26-sin(M_PI*(v+c-1)/2.0)*sin(arg26)/nen26);  //Z=17656
            const double a22as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1;  //Z=17657
            const double F22as23z = a22as23z*(cos(M_PI*(v-1-c)/2.0)*cos(arg25)/nen25-sin(M_PI*(v-1-c)/2.0)*sin(arg25)/nen25+cos(M_PI*(v-1+c)/2.0)*cos(arg26)/nen26-sin(M_PI*(v-1+c)/2.0)*sin(arg26)/nen26);  //Z=17658
            const double a22as24z = (1/2.0)*ee1*e1*(1/(2.0*x1z))*(1/(2.0*x2z))*pzvc2;  //Z=17659
            const double F22as24z = a22as24z*(cos(M_PI*(v-1-c+1)/2.0)*cos(arg27)/nen27-sin(M_PI*(v-1-c+1)/2.0)*sin(arg27)/nen27+cos(M_PI*(v-1+c-1)/2.0)*cos(arg28)/nen28-sin(M_PI*(v-1+c-1)/2.0)*sin(arg28)/nen28);  //Z=17660
            const double F22as20z = preg1*preg3*pow(x1z,v)*pow(x2z,c);  //Z=17661
            const double F22as2z = F22as20z*(F22as21z+F22as22z+F22as23z+F22as24z);  //Z=17662
            //F22asz = F22as1z+F22as2z;  //Z=17663
            F22 = F22as1z0+F22as2z;  //Z=17664
            //F22 = F22asz0;  //Z=17665
        }/*4*/  //Z=17666

        /* ** term #3 asymptote ** */  //Z=17668
        if ( xradp>=lim3 )
        {/*4*/  //Z=17669
            //arg31 = (zr+v-2*a1+1)*atan(2.0*x1z);  //Z=17670
            //nen31 = pow(1.0+4*x1z*x1z,(zr+v-2*a1+1)/2.0);  //Z=17671
            //arg32 = (zr+v-2*a1)*atan(2.0*x1z);  //Z=17672
            //nen32 = pow(1.0+4*x1z*x1z,(zr+v-2*a1)/2.0);  //Z=17673
            //F32as1sum1z = dnv0*ee0*pvav0*(cos(M_PI*v/2.0)*cos(arg31)/nen31-sin(M_PI*v/2.0)*sin(arg31)/nen31);  //Z=17674
            //F32as1sum2z = dnv0*ee1*(1/(2.0*x1z))*pvav10*(cos(M_PI*(v-1)/2.0)*cos(arg32)/nen32-sin(M_PI*(v-1)/2.0)*sin(arg32)/nen32);  //Z=17675
            const double F32as10z = preg1*preg4*pow(x1z,v)*pow(x12z,-a1);  //Z=17676
            //F32as1z = F32as10z*(F32as1sum1z+F32as1sum2z);  //Z=17677

            const double arg310 = (z+v-2*a1+1)*atan(2.0*x1z);  //Z=17679
            const double nen310 = pow(1.0+4*x1z*x1z,(z+v-2*a1+1)/2.0);  //Z=17680
            const double arg320 = (z+v-2*a1)*atan(2.0*x1z);  //Z=17681
            const double nen320 = pow(1.0+4*x1z*x1z,(z+v-2*a1)/2.0);  //Z=17682
            const double F32as1sum1z0 = ee0*pzva*(cos(M_PI*v/2.0)*cos(arg310)/nen310-sin(M_PI*v/2.0)*sin(arg310)/nen310);  //Z=17683
            const double F32as1sum2z0 = ee1*(1/(2.0*x1z))*pzva1*(cos(M_PI*(v-1)/2.0)*cos(arg320)/nen320-sin(M_PI*(v-1)/2.0)*sin(arg320)/nen320);  //Z=17684
            const double F32as1z0 = F32as10z*(F32as1sum1z0+F32as1sum2z0);  //Z=17685

            const double arg33 = (zr+v+c+1)*atan(4.0*x1z);  //Z=17687
            const double nen33 = pow(1.0+16*x1z*x1z,(zr+v+c+1)/2.0);  //Z=17688
            const double arg34 = (zr+v+c)*atan(4.0*x1z);  //Z=17689
            const double nen34 = pow(1.0+16*x1z*x1z,(zr+v+c)/2.0);  //Z=17690
            const double arg35 = (zr+v+c-1)*atan(4.0*x1z);  //Z=17691
            const double nen35 = pow(1.0+16*x1z*x1z,(zr+v+c-1)/2.0);  //Z=17692
            const double F32as21z = (1/2.0)*ee0*e0*pzvc*(cos(M_PI*(v-c)/2.0)+cos(M_PI*(v+c)/2.0)*cos(arg33)/nen33-sin(M_PI*(v+c)/2.0)*sin(arg33)/nen33);  //Z=17693
            const double F32as22z = (1/2.0)*ee0*e1*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-c+1)/2.0)+cos(M_PI*(v+c-1)/2.0)*cos(arg34)/nen34-sin(M_PI*(v+c-1)/2.0)*sin(arg34)/nen34);  //Z=17694
            const double F32as23z = (1/2.0)*ee1*e0*(1/(2.0*x1z))*pzvc1*(cos(M_PI*(v-1-c)/2.0)+cos(M_PI*(v-1+c)/2.0)*cos(arg34)/nen34-sin(M_PI*(v-1+c)/2.0)*sin(arg34)/nen34);  //Z=17695
            const double F32as24z = (1/2.0)*ee1*e1*(1/(4.0*x1z*x1z))*pzvc2*(cos(M_PI*(v-1-c+1)/2.0)+cos(M_PI*(v-1+c-1)/2.0)*cos(arg35)/nen35-sin(M_PI*(v-1+c-1)/2.0)*sin(arg35)/nen35);  //Z=17696
            const double F32as20z = preg1*preg3*pow(x1z,v)*pow(x1z,c);  //Z=17697
            const double F32as2z = F32as20z*(F32as21z+F32as22z+F32as23z+F32as24z);  //Z=17698
            //F32asz = F32as1z+F32as2z;  //Z=17699
            F32 = F32as1z0+F32as2z;  //Z=17700
            //F32 = F32asz0;  //Z=17701
        }/*4*/  //Z=17702


        /* ** term #4 asymptote ** */  //Z=17705
        if ( xrad>=lim4 )
        {/*4*/  //Z=17706
            const double F42as10z = preg4*preg4*pow(x22z,-2*a1);  //Z=17707
            //F42as1sumz = pva0;  //Z=17708
            //F42as1z = F42as10z*F42as1sumz;  //Z=17709
            const double F42as1z0 = F42as10z*pza;  //Z=17710

            const double arg41 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=17712
            const double nen41 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=17713
            const double arg42 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=17714
            const double nen42 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=17715
            //arg43 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=17716
            //nen43 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=17717
            const double F42as20z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=17718
            const double F42as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=17719
            const double F42as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c-1)/2.0)*sin(arg42)/nen42);  //Z=17720
            //F42as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg43)/nen43-sin(M_PI*c/2.0)*sin(arg43)/arg43);  //Z=17721
            //F42as2z = F42as20z*(F42as21+F42as22+F42as23);  //Z=17722
            const double F42as2z0 = F42as20z*(F42as21+F42as22);  //Z=17723

            const double F42as30z = preg4*preg3*pow(x22z,-a1)*pow(x2z,c);  //Z=17725
            const double F42as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg41)/nen41-sin(M_PI*c/2.0)*sin(arg41)/nen41);  //Z=17726
            //F42as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c-1)/2.0)*cos(arg43)/nen43-sin(M_PI*(c-1)/2.0)*sin(arg43)/nen43);  //Z=17727
            const double F42as26 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c+1)/2.0)*cos(arg42)/nen42-sin(M_PI*(c+1)/2.0)*sin(arg42)/nen42);  //Z=17728
            //F42as3z = F42as30z*(F42as24+F42as25+F42as26);  //Z=17729
            const double F42as3z0 = F42as30z*(F42as24+F42as26);  //Z=17730

            const double F42as40z = preg3*preg3*pow(x2z*x2z,c);  //Z=17732
            const double arg44 = (zr+2*c+1)*atan(4.0*x2z);  //Z=17733
            const double nen44 = pow(1.0+16*x2z*x2z,(zr+2*c+1)/2.0);  //Z=17734
            const double arg45 = (zr+2*c)*atan(4.0*x2z);  //Z=17735
            const double nen45 = pow(1.0+16*x2z*x2z,(zr+2*c)/2.0);  //Z=17736
            const double F42as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg44)/nen44-sin(M_PI*c)*sin(arg44)/nen44);  //Z=17737
            const double F42as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=17738
            const double F42as29 = (1/2.0)*e1*e0*(1/(2.0*x2z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg45)/nen45-sin(M_PI*(2*c-1)/2.0)*sin(arg45)/nen45);  //Z=17739
            const double F42as4z = F42as40z*(F42as27+F42as28+F42as29);  //Z=17740
            //F42asz = F42as1z+F42as2z+F42as3z+F42as4z;  //Z=17741
            F42 = F42as1z0+F42as2z0+F42as3z0+F42as4z;  //Z=17742
            //F42 = F42asz0;  //Z=17743
        }/*4*/  //Z=17744


        /* ** term #5 asymptote ** */  //Z=17747
        if ( xradp>=lim5 )
        {/*4*/  //Z=17748
            const double F52as10z = preg4*preg4*pow(x12z,-a1)*pow(x22z,-a1);  //Z=17749
            //F52as1sumz = pva0;  //Z=17750
            //F52as1z = F52as10z*F52as1sumz;  //Z=17751
            const double F52as1z0 = F52as10z*pza;  //Z=17752

            const double F52as20z = preg4*preg3*pow(x12z,-a1)*pow(x2z,c);  //Z=17754
            const double arg51 = (zr-2*a1+c+1)*atan(2.0*x2z);  //Z=17755
            const double nen51 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+1)/2.0);  //Z=17756
            const double arg52 = (zr-2*a1+c)*atan(2.0*x2z);  //Z=17757
            const double nen52 = pow(1.0+4*x2z*x2z,(zr-2*a1+c)/2.0);  //Z=17758
            //arg53 = (zr-2*a1+c+3)*atan(2.0*x2z);  //Z=17759
            //nen53 = pow(1.0+4*x2z*x2z,(zr-2*a1+c+3)/2.0);  //Z=17760
            const double F52as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg51)/nen51-sin(M_PI*c/2.0)*sin(arg51)/nen51);  //Z=17761
            const double F52as22 = d0*e1*pzac1*(1/(2.0*x2z))*(cos(M_PI*(c-1)/2.0)*cos(arg52)/nen52-sin(M_PI*(c-1)/2.0)*sin(arg52)/nen52);  //Z=17762
            //F52as23 = d1*e0*pzac2*(-x22z)*(cos(M_PI*c/2.0)*cos(arg53)/nen53-sin(M_PI*c/2.0)*sin(arg53)/nen53);  //Z=17763
            //F52as2z = F52as20z*(F52as21+F52as22+F52as23);  //Z=17764
            const double F52as2z0 = F52as20z*(F52as21+F52as22);  //Z=17765

            const double F52as30z = preg4*preg3*pow(x22z,-a1)*pow(x1z,c);  //Z=17767
            const double arg54 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=17768
            const double nen54 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=17769
            //arg55 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=17770
            //nen55 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=17771
            const double arg56 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=17772
            const double nen56 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=17773
            const double F52as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg54)/nen54-sin(M_PI*c/2.0)*sin(arg54)/nen54);  //Z=17774
            //F52as25 = d1*e0*pzac2*(-x22z)*(cos(M_PI*(c+1)/2.0)*cos(arg55)/nen55-sin(M_PI*(c+1)/2.0)*sin(arg55)/nen55);  //Z=17775
            const double F52as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg56)/nen56-sin(M_PI*(c-1)/2.0)*sin(arg56)/nen56);  //Z=17776
            //F52as3z = F52as30z*(F52as24+F52as25+F52as26);  //Z=17777
            const double F52as3z0 = F52as30z*(F52as24+F52as26);  //Z=17778

            const double F52as40z = preg3*preg3*pow(x1z,c)*pow(x2z,c);  //Z=17780
            const double arg57 = (zr+2*c+1)*atan(2.0*(x1z-x2z));  //Z=17781
            const double nen57 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c+1)/2.0);  //Z=17782
            const double arg58 = (zr+2*c+1)*atan(2.0*(x1z+x2z));  //Z=17783
            const double nen58 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c+1)/2.0);  //Z=17784
            const double arg59 = (zr+2*c)*atan(2.0*(x1z-x2z));  //Z=17785
            const double nen59 = pow(1.0+4*(x1z-x2z)*(x1z-x2z),(zr+2*c)/2.0);  //Z=17786
            const double arg510 = (zr+2*c)*atan(2.0*(x1z+x2z));  //Z=17787
            const double nen510 = pow(1.0+4*(x1z+x2z)*(x1z+x2z),(zr+2*c)/2.0);  //Z=17788
            const double F52as27 = (1/2.0)*e0*e0*pzc*(cos(M_PI*(c-c)/2.0)*cos(arg57)/nen57-sin(M_PI*(c-c)/2.0)*sin(arg57)/nen57+cos(M_PI*c)*cos(arg58)/nen58-sin(M_PI*c)*sin(arg58)/nen58);  //Z=17789
            const double F52as28 = (1/2.0)*e0*e1*(1/(2.0*x2z))*pzc1*(0+sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=17790
            const double F52as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(0-sin(arg59)/nen59+cos(M_PI*(2*c-1)/2.0)*cos(arg510)/nen510-sin(M_PI*(2*c-1)/2.0)*sin(arg510)/nen510);  //Z=17791
            const double F52as4z = F52as40z*(F52as27+F52as28+F52as29);  //Z=17792
            //F52asz = F52as1z+F52as2z+F52as3z+F52as4z;  //Z=17793
            F52 = F52as1z0+F52as2z0+F52as3z0+F52as4z;  //Z=17794
            //F52 = F52asz0;  //Z=17795
        }/*4*/  //Z=17796

        /* ** term #6 asymptote ** */  //Z=17798
        if ( xradp>=lim6 )
        {/*4*/  //Z=17799
            const double F62as10z = preg4*preg4*pow(x12z,-a1)*pow(x12z,-a1);  //Z=17800
            //F62as1sumz = pva0;  //Z=17801
            //F62as1z = F62as10z*F62as1sumz;  //Z=17802
            const double F62as1z0 = F62as10z*pza;  //Z=17803

            const double F62as20z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=17805
            const double arg61 = (zr-2*a1+c+1)*atan(2.0*x1z);  //Z=17806
            const double nen61 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+1)/2.0);  //Z=17807
            const double arg62 = (zr-2*a1+c)*atan(2.0*x1z);  //Z=17808
            const double nen62 = pow(1.0+4*x1z*x1z,(zr-2*a1+c)/2.0);  //Z=17809
            //arg63 = (zr-2*a1+c+3)*atan(2.0*x1z);  //Z=17810
            //nen63 = pow(1.0+4*x1z*x1z,(zr-2*a1+c+3)/2.0);  //Z=17811
            const double F62as21 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=17812
            const double F62as22 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=17813
            //F62as23 = d1*e0*pzac2*(-x12z)*(cos(M_PI*c/2.0)*cos(arg63)/nen63-sin(M_PI*c/2.0)*sin(arg63)/nen63);  //Z=17814
            //F62as2z = F62as20z*(F62as21+F62as22+F62as23);  //Z=17815
            const double F62as2z0 = F62as20z*(F62as21+F62as22);  //Z=17816

            const double F62as30z = preg4*preg3*pow(x12z,-a1)*pow(x1z,c);  //Z=17818
            const double F62as24 = d0*e0*pzac*(cos(M_PI*c/2.0)*cos(arg61)/nen61-sin(M_PI*c/2.0)*sin(arg61)/nen61);  //Z=17819
            //F62as25 = d1*e0*pzac2*(-x12z)*(cos(M_PI*(c+1)/2.0)*cos(arg63)/nen63-sin(M_PI*(c+1)/2.0)*sin(arg63)/nen63);  //Z=17820
            const double F62as26 = d0*e1*pzac1*(1/(2.0*x1z))*(cos(M_PI*(c-1)/2.0)*cos(arg62)/nen62-sin(M_PI*(c-1)/2.0)*sin(arg62)/nen62);  //Z=17821
            //F62as3z = F62as30z*(F62as24+F62as25+F62as26);  //Z=17822
            const double F62as3z0 = F62as30z*(F62as24+F62as26);  //Z=17823

            const double F62as40z = preg3*preg3*pow(x1z*x1z,c);  //Z=17825
            const double arg64 = (zr+2*c+1)*atan(4.0*x1z);  //Z=17826
            const double nen64 = pow(1.0+16*x1z*x1z,(zr+2*c+1)/2.0);  //Z=17827
            const double arg65 = (zr+2*c)*atan(4.0*x1z);  //Z=17828
            const double nen65 = pow(1.0+16*x1z*x1z,(zr+2*c)/2.0);  //Z=17829
            const double F62as27 = (1/2.0)*e0*e0*pzc*(1+cos(M_PI*c)*cos(arg64)/nen64-sin(M_PI*c)*sin(arg64)/nen64);  //Z=17830
            const double F62as28 = (1/2.0)*e0*e1*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=17831
            const double F62as29 = (1/2.0)*e1*e0*(1/(2.0*x1z))*pzc1*(cos(M_PI*(2*c-1)/2.0)*cos(arg65)/nen65-sin(M_PI*(2*c-1)/2.0)*sin(arg65)/nen65);  //Z=17832
            const double F62as4z = F62as40z*(F62as27+F62as28+F62as29);  //Z=17833
            //F62asz = F62as1z+F62as2z+F62as3z+F62as4z;  //Z=17834
            F62 = F62as1z0+F62as2z0+F62as3z0+F62as4z;  //Z=17835
            //F62 = F62asz0;  //Z=17836
        }/*4*/  //Z=17837

        /*formpq:=*/ return pql*(cc1*F12+cc2*F22+cc3*F32+cc4*F42+cc5*F52+cc6*F62)/vv3;  //Z=17839

    }/*3*/ /*  of inhomogeneous core/shell-disk  */  //Z=17845

    /*  myelin disk  */  //Z=17847
    if ( (params.cs==3) || (params.cs==4) )
    {/*3*/  //Z=17848

        /*  disk parameters  */  //Z=17850
        const double v = -1;  //Z=17851
        const double e0 = 1;  //Z=17852
        const double e1 = 0;  //Z=17853
        const double preg1 = 1/2.0;  //Z=17854
        const double pz2v = 1/(zr*(zr-1));  //Z=17855
        const double pz2v1 = pz2v/(zr-2);  //Z=17856
        const double pz2v2 = pz2v1/(zr-3);  //Z=17857
        const double lim = 18*exp(-5*params.sigma);  //Z=17858
        const double lim1 = lim*1.2;  //Z=17859
        const double rad = params.CR->myarray[1];  //Z=17860
        const int    inmax = round(params.CR->myarray[14]);  //Z=17861
        const double vvm = params.CR->myarray[15];  //Z=17862
        const double rmax = params.CR->myarray[16];  //Z=17863
        const double xmax = q*rmax;  //Z=17864

        double F12;

        if ( xmax<(lim1) )
        {/*4*/  //Z=17866
            // Auf der GPU ist es effizienter, eine kleine Berechnung zu machen als ein grosses Array zu haben.
            /* fkv[0]:=1;  //Z=17867 */
            double qqnn;  //Z=17868

            F12 = 0.0;  //Z=17874
            for ( int ii=1; ii<=inmax; ii++ )
            {/*5*/  //Z=17875
                for ( int jj=1; jj<=inmax; jj++ )
                {/*6*/  //Z=17876
                    double F12sez = 1.0;  //Z=17877
                    double oldF12sez = 1.0;  //Z=17878
                    qqnn = 1.0;
                    for ( int nser=1; nser<=120; nser++ )
                    {/*7*/  //Z=17879
                        qqnn = qqnn * q*q;
                        double pqsum = 0;  //Z=17880
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*8*/  //Z=17881
                            /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))/((mser+1)*fkv[mser]*(nser-mser+1)*fkv[nser-mser]*fkv[mser]*fkv[nser-mser]);  //Z=17882 */
                            pqsum += pow(params.CR->carr7p[ii],2*mser)*pow(params.CR->carr7p[jj],2*(nser-mser))/(params.CR->carr6p[mser]*params.CR->carr6p[nser-mser]);  //Z=17883

                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17885 */
                            /* pqsum:=pqsum+power(carr7p[ii],2*mser)*power(carr7p[jj],2*(nser-mser))*carr1pm[indx];  //Z=17886 */
                        }/*8*/  //Z=17887
                        F12sez += params.CR->carr4p[nser]*qqnn*pqsum;  //Z=17888
                        const double delser = fabs((F12sez-oldF12sez)/F12sez);  //Z=17889
                        if ( delser<0.0001 ) break; /* goto 251; */  //Z=17890
                        oldF12sez = F12sez;  //Z=17891
                    }/*7*/  //Z=17892
                    /*251:*/  //Z=17893
                    F12 += params.CR->carr5p[ii]*params.CR->carr5p[jj]*F12sez;  //Z=17894
                }/*6*/  //Z=17895
            }/*5*/  //Z=17896
            F12 = F12/vvm;  //Z=17897
            //F12 = F12ser;  //Z=17898
        }/*4*/  //Z=17899
        else
        {/*4*/  //Z=17900
            const double xrz = q*rad/(zr+1);  //Z=17901
            const double arg = (zr+2*v+1)*atan(2.0*xrz);  //Z=17902
            const double nen = pow(1.0+4*xrz*xrz,(zr+2*v+1)/2.0);  //Z=17903
            const double arg1 = (zr+2*v)*atan(2.0*xrz);  //Z=17904
            const double nen1 = pow(1.0+4*xrz*xrz,(zr+2*v)/2.0);  //Z=17905
            const double arg2 = (zr+2*v-1)*atan(2.0*xrz);  //Z=17906
            const double nen2 = pow(1.0+4*xrz*xrz,(zr+2*v-1)/2.0);  //Z=17907

            double F12asz = 0.0;  //Z=17909
            for ( int ii=1; ii<=inmax; ii++ )
            {/*5*/  //Z=17910
                const double a1m = params.CR->carr5p[ii]*pow(params.CR->carr7p[ii],v);   /*  carr7p[ii]:=pp[ii];  //Z=17911 */
                for ( int jj=1; jj<=inmax; jj++ )
                {/*6*/  //Z=17912
                    const double a2m = params.CR->carr5p[jj]*pow(params.CR->carr7p[jj],v);  //Z=17913
                    const double xijm = (params.CR->carr3p[ii]-params.CR->carr3p[jj])*q/(zr+1);      /*   carr3p[ii]:=ll[ii];  //Z=17914 */
                    const double arglmz = (zr+1)*atan(xijm);  //Z=17915
                    const double nenlmz = pow(1.0+xijm*xijm,(zr+1)/2.0);  //Z=17916
                    const double xijp = (params.CR->carr3p[ii]+params.CR->carr3p[jj])*q/(zr+1);  //Z=17917
                    const double arglpz = (zr+1)*atan(xijp);  //Z=17918
                    const double nenlpz = pow(1.0+xijp*xijp,(zr+1)/2.0);  //Z=17919
                    const double F12as1z = e0*e0*pz2v*(cos(arglmz)/nenlmz+(cos(M_PI*v)*(cos(arg)*cos(arglpz)-sin(arg)*sin(arglpz))-sin(M_PI*v)*(sin(arg)*cos(arglpz)+cos(arg)*sin(arglpz)))/(nen*nenlpz));  //Z=17920
                    const double F12as2z = e0*e1*(1/(params.CR->carr7p[jj]*xrz))*pz2v1*(-sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=17921
                    const double F12as3z = e1*e0*(1/(params.CR->carr7p[ii]*xrz))*pz2v1*(sin(arglmz)/nenlmz+(cos(M_PI*(2*v-1)/2.0)*(cos(arg1)*cos(arglpz)-sin(arg1)*sin(arglpz))-sin(M_PI*(2*v-1)/2.0)*(sin(arg1)*cos(arglpz)+cos(arg1)*sin(arglpz)))/(nen1*nenlpz));  //Z=17922
                    const double F12as4z = e1*e1*(1/(params.CR->carr7p[ii]*params.CR->carr7p[jj]*xrz*xrz))*pz2v2*(cos(arglmz)/nenlmz+(cos(M_PI*(v-1))*(cos(arg2)*cos(arglpz)-sin(arg2)*sin(arglpz))-sin(M_PI*(v-1))*(sin(arg2)*cos(arglpz)+cos(arg2)*sin(arglpz)))/(nen2*nenlpz));  //Z=17923

                    F12asz += a1m*a2m*(F12as1z+F12as2z+F12as3z+F12as4z);  //Z=17925
                }/*6*/  //Z=17926
            }/*5*/  //Z=17927
            F12 = preg1*preg1*pow(xrz/2.0,2*v)*(1/2.0)*F12asz/vvm;  //Z=17928
            //F12 = F12asy;  //Z=17929
        }/*4*/  //Z=17930
        /*formpq:=*/ return pql*F12;  //Z=17931

    }/*3*/ /*  of myelin disk  */  //Z=17937

    //}/*2*/ /*  of disk  */  //Z=17939

    return 0.0;
}  //Z=18387


#endif // SC_LIB_FORMPQ_partDisk_H
