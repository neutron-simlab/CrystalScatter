#ifndef SC_LIB_FORMPQ_partTriaxEllips_H
#define SC_LIB_FORMPQ_partTriaxEllips_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partTriaxEllips(double qxs, double qys, double q) const   /*Z=14910*/
{/*1*/  //Z=15188

    double pqsum, oldpqsum, binsum, delser, pql;  //Z=15195
    double zr, argqx, argqy, pqrx, pqry;  //Z=15196

    /*begin*/  //Z=15230
    zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232

    double pq;
    double qxn[121], qyn[121];
    double qz = 1;

    CHECKENDTHREAD_VAL;

    /* ** triaxial ellipsoid ** */  //Z=18161
    //if ( params.part==6 )
    //{/*2*/  //Z=18162
        /*  homogeneous isotropic triaxial ellipsoid  */  //Z=18163
        if ( params.ordis==7 )
        {/*3*/  //Z=18164
            if ( params.cs==0 )
            {/*4*/  //Z=18165
                if ( q<0.05*params.limq4 )
                {/*5*/  //Z=18166
                    pqsum = 1.0;  //Z=18167
                    oldpqsum = 0.0;  //Z=18168
                    double qqnn = 1.0;  //Z=18169
                    for ( int nser=1; nser<=120; nser++ )
                    {/*6*/  //Z=18170
                        qqnn = qqnn*q*q;  //Z=18171
                        pqsum = pqsum+params.CR->carr4p[nser]*qqnn;  //Z=18172
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18173
                        if ( delser<0.0001 ) break; /* goto 263; */  //Z=18174
                        oldpqsum = pqsum;  //Z=18175
                    }/*6*/  //Z=18176
                    /*263:*/  //Z=18177
                    /*formpq:=*/ return pqsum;  //Z=18178
                }/*5*/  //Z=18179
                else
                {/*5*/  //Z=18180
                    if ( q>2*params.limq4 )
                        pq = params.por/(q*q*q*q);  //Z=18181
                    else
                    {/*6*/  //Z=18182
                        qrombdeltac(params.p1,params.sigma,params.alphash1,params.polTheta,0,qxs,qys,qz,9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,3,7,14,7,0,0,params.CR->carr1p,pql);  //Z=18183
                        pq = pql/(M_PI/2.0);  //Z=18184
                    }/*6*/  //Z=18185
                    /*formpq:=*/ return pq;  //Z=18186
                }/*5*/  //Z=18187
            }/*4*/ /*  of homogeneous isotropic ellipsoid */  //Z=18188

#ifdef procnotused
            /*  core/shell isotropic ellipsoid  */  //Z=18190
            if ( params.cs==1 )
            {/*4*/  //Z=18191
                /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=18192
            }/*4*/  //Z=18193
#endif
        }/*3*/  /*  of isotropic triaxial ellipsoid  */  //Z=18194

        /*  perfectly oriented triaxial ellipsoid  */  //Z=18196
        if ( params.ordis==6 )
        {/*3*/  //Z=18197
            /* if (orcase=1) then begin  //Z=18198 */
            if ( 1==1 )
            {/*4*/  //Z=18199

                if ( q<(1/params.radius) )
                {/*5*/  //Z=18201
                    pqsum = 1.0;  //Z=18202
                    oldpqsum = 0.0;  //Z=18203
                    qxn[0] = 1.0;  //Z=18204
                    qyn[0] = 1.0;  //Z=18205
                    for ( int nser=1; nser<=80; nser++ )
                    {/*6*/  //Z=18206
                        qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=18207
                        qyn[nser] = qyn[nser-1]*qys*qys;  //Z=18208

                        binsum = 0.0;  //Z=18210
                        for ( int mser=0; mser<=nser; mser++ )
                        {/*7*/  //Z=18211
                            /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=18212 */
                            /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=18213 */
                            binsum = binsum+params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=18214
                        }/*7*/  //Z=18215
                        pqsum = pqsum+binsum;  //Z=18216
                        delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18217
                        if ( delser<0.0001 ) break; /* goto 264; */  //Z=18218
                        oldpqsum = pqsum;  //Z=18219
                    }/*6*/  //Z=18220
                    /*264:*/  //Z=18221
                    pql = pqsum;  //Z=18222
                }/*5*/  //Z=18223
                else
                {/*5*/  //Z=18224
                    argqx = qxs*params.radius/(zr+1)+eps9;  //Z=18225
                    argqy = qys*params.radius/(zr+1)+eps9;  //Z=18226
                    pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18227
                    pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18228
                    pql = pqrx*pqry;  //Z=18229
                }/*5*/  //Z=18230
                /*formpq:=*/ return pql;  //Z=18231
            }/*4*/  /*  of orcase=1  */  //Z=18232
        }/*3*/  /*  of perfect triaxial ellipsoid  */  //Z=18233
    //}/*2*/  /*  of triaxial ellipsoid  */  //Z=18234

    return 0.0;
}/*1*/  //Z=18387


#endif // SC_LIB_FORMPQ_partTriaxEllips_H
