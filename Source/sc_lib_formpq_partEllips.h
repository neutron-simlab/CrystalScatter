#ifndef SC_LIB_FORMPQ_partEllips_H
#define SC_LIB_FORMPQ_partEllips_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partEllips(double sigmal, double qx, double qy, double qxs, double qys, double q) const   /*Z=14910*/
{/*1*/  //Z=15188
    // double sigmal ist nicht zu ersetzen, da es an einer Stelle CALC.epsilon ist, sonst nur CALC.params.sigmal

    double pqsum, oldpqsum, binsum, delser, pql;  //Z=15195
    double zr;  //Z=15196
    double a1;  //Z=15207

    /*begin*/  //Z=15230
    zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232

    // Noch fehlende (globale) Variablen und sonstige Anpassungen:
    // An machen Stellen muss "params." eingefügt werden.
    // Die Konstante 'eps' muss durch 'eps9' ersetzt werden.
    double /*binsum1,*/ pq;
    double qxn[121], qyn[121];
    double qz = 1; // wird für qrombchid() verwendet
    // Dort werden auch (einmal) p11..p33, ax1[nxyz],ax2[nxyz],ax3[nxyz], sig[xyz] verwendet

    CHECKENDTHREAD_VAL;

    /* ** biaxial ellipsoid ** */  //Z=18011
    //if ( params.part==5 )
    //{/*2*/  //Z=18012
        //qDebug() << ordis << params.cs;
        /*  homogeneous isotropic ellipsoid  */  //Z=18013
        if ( params.ordis==7 )
        {/*3*/  //Z=18014
            if ( params.cs==0 )
            {/*4*/  //Z=18015
                if ( q<0.8*params.limq4 ) //20240301 - war 0.6  // 20221004 - crystal3d1.pas: 0.8
                {/*5*/  //Z=18016
                    pqsum = 1.0;  //Z=18017
                    oldpqsum = 0.0;  //Z=18018
                    double qqnn = 1.0;  //Z=18019
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18020
                        qqnn = qqnn*q*q;  //Z=18021
                        pqsum = pqsum+params.CR->carr4p[nser]*qqnn;  //Z=18022
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18023
                        if ( delser<0.0001 ) break; /* goto 260; */  //Z=18024
                        oldpqsum = pqsum;  //Z=18025
                    }/*6*/  //Z=18026
                    /*260:*/  //Z=18027
                    /*formpq:=*/ return pqsum;  //Z=18028
                }/*5*/  //Z=18029
                else
                {/*5*/  //Z=18030
                    if ( q>=1.5*params.limq4 )
                        pq = params.por/(q*q*q*q);  //Z=18031
                    else
                    {/*6*/  //Z=18032
                        //qrombchid(l,r,p1,sigma,alfa,dbeta,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2053
                        double qhkl = 1;
                        double epsi = sigmal;  // TODO: passt das?
                        //qDebug() << "EPSI:" << q << 0.6*params.limq4 << 1.5*params.limq4 << ":" << epsi;
                        qrombchid(params.length,params.radius,params.p1,params.sigma,params.alphash1,epsi,params.polTheta,params.polPhi,
                                  qx,qy,qz,
                                  params.p11,params.p12,params.p13,params.p21,params.p22,params.p23,params.p31,params.p32,params.p33,
                                  qx,qy,0,qhkl,
                                  //params.ax1.length(),params.ax2.length(),params.ax3.length(),
                                  //params.ax1.x(),params.ax1.y(),params.ax1.z(),
                                  //params.ax2.x(),params.ax2.y(),params.ax2.z(),
                                  //params.ax3.x(),params.ax3.y(),params.ax3.z(),
                                  //params.sig.x(),params.sig.y(),params.sig.z(),
                                  /*params.ordis,3,8,*/13,7,0,0,params.CR->carr1p,pql);  //Z=18033
                        pq = pql;  //Z=18034
                    }/*6*/  //Z=18035
                    /*formpq:=*/ return pq;  //Z=18036
                }/*5*/  //Z=18037
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18038

#ifdef procnotused
            /*  core/shell isotropic ellipsoid  */  //Z=18040
            if ( params.cs==1 )
            {/*4*/  //Z=18041
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18042
            }/*4*/  //Z=18043
#endif
        }/*3*/  /*  of isotropic ellipsoid  */  //Z=18044

        /*  perfectly oriented ellipsoid  */  //Z=18046
        if ( params.ordis==6 )
        {/*3*/  //Z=18047
            /* if (orcase=1) then begin  //Z=18048 */
            if ( 1==1 )
            {/*4*/  //Z=18049

                if ( sqrt(qx*qx*sqr(params.length)+qy*qy*sqr(params.radius)+eps9)<15 )
                {/*5*/  //Z=18051
                    qxn[0] = 1.0;  //Z=18052
                    qyn[0] = 1.0;  //Z=18053
                    for ( int nser=1; nser<=81; nser++ )
                    {/*6*/  //Z=18054
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18055
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18056
                    }/*6*/  //Z=18057
                    pqsum = 0.0;  //Z=18058
                    oldpqsum = 0.0;  //Z=18059
                    for ( int nser=0; nser<=80; nser++ )
                    {/*6*/  //Z=18060
                        binsum = 0.0;  //Z=18061
                        for ( int mser=0; mser<=80; mser++ ) binsum = binsum+params.CR->carr5p[mser]*params.CR->carr11pm[nser][mser]*qyn[mser];  //Z=18062
                        pqsum = pqsum+params.CR->carr4p[nser]*qxn[nser]*binsum;  //Z=18063
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18064
                        if ( delser<0.0001 ) break; /* goto 261; */  //Z=18065
                        oldpqsum = pqsum;  //Z=18066
                    }/*6*/  //Z=18067
                    /*261:*/  //Z=18068
                    pql = pqsum;  //Z=18069
                }/*5*/  //Z=18070
                else
                {/*5*/  //Z=18071
                    a1 = 9*pow(zr+1,4)/(2.0*zr*(zr-1)*(zr-2)*(zr-3));  //Z=18072
                    pql = a1/(sqr(qy*qy*sqr(params.radius)+qx*qx*sqr(params.length))+eps9);  //Z=18073
                }/*5*/  //Z=18074
                /*formpq:=*/ return pql;  //Z=18075
            }/*4*/  /*  of orcase=1  */  //Z=18076
        }/*3*/  /*  of perfect ellipsoid  */  //Z=18077

        /*  general  */  //Z=18080
        if ( params.ordis==0 )
        {/*3*/  //Z=18081
            if ( params.orcase==4 )
                pql = 1.0;  //Z=18082
            else
            {/*4*/  //Z=18083
                if ( sqrt(qx*qx*sqr(params.length)+qy*qy*sqr(params.radius)+eps9)<10 )
                {/*5*/  //Z=18084
                    double pqsum = 0.0;  //Z=18085
                    double oldpqsum = -10.0;  //Z=18086
                    qxn[0] = 1.0;  //Z=18087
                    qyn[0] = 1.0;  //Z=18088
                    // Auf der GPU ist es effizienter, eine kleine Berechnung zu machen als ein grosses Array zu haben.
                    // ==> Da hier aber zwei unterschiedliche Schleifen mit diesen Werten versorgt werden, lasse ich
                    //     im Moment mal ein lokales, in der Größe angepasstes Array stehen
                    double qqn[121];
                    qqn[0] = 1.0;  //Z=18089

                    if ( params.orcase==1 )
                    {/*6*/  //Z=18091
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18092
                            qxn[nser] = qxn[nser-1]*qxs*qxs/(q*q);  //Z=18093
                            qyn[nser] = qyn[nser-1]*qys*qys/(q*q);  //Z=18094
                            qqn[nser] = qqn[nser-1]*q*q;  //Z=18095
                        }/*7*/  //Z=18096
                        for ( int nser=0; nser<=120; nser++ )
                        {/*7*/  //Z=18097
                            double binsum1 = 0.0;  //Z=18098
                            for ( int lser=0; lser<=nser; lser++ )  //Z=18099
                                binsum1 += params.CR->carr1p[lser]*qxn[lser]*qyn[nser-lser];  //Z=18100
                            double binsum = 0.0;  //Z=18101
                            for ( int mser=0; mser<=120; mser++ )  //Z=18102
                                binsum += params.CR->carr2p[mser]*qqn[mser]*params.CR->carr11pm[nser][mser];  //Z=18103
                            pqsum += params.CR->carr3p[nser]*qqn[nser]*binsum*binsum1/pow(4.0,nser);  //Z=18104  TODO es stand: pow(4,n)
                            const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18105
                            if ( delser<0.0001 ) break; /* goto 273; */  //Z=18106
                            oldpqsum = pqsum;  //Z=18107
                        }/*7*/  //Z=18108
                    }/*6*/  //Z=18109

                    if ( params.orcase==2 )
                    {/*6*/  /*  x-axis  */  //Z=18111
                        for ( int nser=1; nser<=110; nser++ )
                        {/*7*/  //Z=18112
                            qxn[nser] = qxn[nser-1]*qxs*qxs/(q*q);  //Z=18113
                            qyn[nser] = qyn[nser-1]*qys*qys/(q*q);  //Z=18114
                            qqn[nser] = qqn[nser-1]*q*q;  //Z=18115
                        }/*7*/  //Z=18116
                        for ( int nser=0; nser<=100; nser++ )
                        {/*7*/  //Z=18117
                            double binsum1 = 0.0;  //Z=18118
                            for ( int lser=0; lser<=nser; lser++ )  //Z=18119
                                /* binsum1:=binsum1+carr4p[lser]*qxn[lser]*qyn[nser-lser];  //Z=18120 */
                                binsum1 += params.CR->carr22pm[nser][lser]*qxn[lser]*qyn[nser-lser];  //Z=18121
                            double binsum = 0.0;  //Z=18122
                            for ( int mser=0; mser<=100; mser++ )  //Z=18123
                                binsum += params.CR->carr5p[mser]*qqn[mser]*params.CR->carr11pm[nser][mser];  //Z=18124
                            pqsum = pqsum+params.CR->carr6p[nser]*qqn[nser]*binsum*binsum1/pow(4.0,nser);  //Z=18125
                            const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18126
                            if ( delser<0.0001 ) break; /* goto 273; */  //Z=18127
                            oldpqsum = pqsum;  //Z=18128
                        }/*7*/  //Z=18129
                    }/*6*/  //Z=18130

                    if ( params.orcase==3 )
                    {/*6*/  /*  y-axis  */  //Z=18132
                        for ( int nser=1; nser<=120; nser++ )
                        {/*7*/  //Z=18133
                            qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18134
                            qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18135
                            double binsum = 0.0;  //Z=18136
                            for ( int mser=0; mser<=nser; mser++ )
                            {/*8*/  //Z=18137
                                /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18138 */
                                /* binsum:=binsum+carr1pm[indx]*qxn[mser]*qyn[nser-mser];  //Z=18139 */
                                binsum += params.CR->carr11pm[nser-mser][mser]*qyn[mser]*qxn[nser-mser];  //Z=18140
                            }/*8*/  //Z=18141
                            pqsum += params.CR->carr1p[nser]*binsum;  //Z=18142
                            const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18143
                            if ( delser<0.0001 ) break; /* goto 273; */  //Z=18144
                            oldpqsum = pqsum;  //Z=18145
                        }/*7*/  //Z=18146
                    }/*6*/  //Z=18147
                    /*273:*/  //Z=18148
                    pql = pqsum;  //Z=18149
                }/*5*/   /*  of q<lim  */  //Z=18150
                else
                {/*5*/  //Z=18151
                    qrombdeltac(params.p1,sigmal,params.alphash1,params.polTheta,0,qxs,qys,qz,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,params.ordis,3,*/4,16,0,0,0,params.CR->carr1p,pql);  //Z=18152
                    pql = pql/params.norm;  //Z=18153
                }/*5*/  //Z=18154
            }/*4*/  /*  of orcase=1,2,3  */  //Z=18155
            /*formpq:=*/ return pql;  //Z=18156
        }/*3*/   /*  of general  */  //Z=18157
    //}/*2*/  /*  of biaxial ellipsoid  */  //Z=18158

    return 0.0;
}/*1*/  //Z=18387


#endif // SC_LIB_FORMPQ_partEllips_H
