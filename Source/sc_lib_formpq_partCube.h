#ifndef SC_LIB_FORMPQ_partCube_H
#define SC_LIB_FORMPQ_partCube_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formpq_partCube(double qxs, double qys, double q, int ordis) const   /*Z=14910*/
{/*1*/  //Z=15188

    double pql;  //Z=15195

    const double zr = (1-sqr(params.sigma))/(sqr(params.sigma));  //Z=15232

    double qxn[121], qyn[121];

    CHECKENDTHREAD_VAL;

    /*  cube  */  //Z=17942
    //if ( params.part==4 )
    //{/*2*/  //Z=17943

    //if ( q<0.02 )
    //    qDebug() << "formpq_partCube: ordis"<<ordis << "cs"<<params.cs << "q"<<q
    //             << "limq4"<<params.limq4 << "zr"<<zr << "por"<<params.por << "radius"<<params.radius;

    /*  homogeneous isotropic cube  */  //Z=17944
    if ( ordis==7 )
    {/*3*/  //Z=17945
        if ( params.cs==0 )
        {/*4*/  //Z=17946
            if ( q<1.2*params.limq4 ) //20240301 - war 0.7
            {/*5*/  //Z=17947
                double pqsum = 1.0;  //Z=17948
                double oldpqsum = 0.0;  //Z=17949
                double qqn = 1.0;  //Z=17950
                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=17951
                    qqn = qqn*q*q;  //Z=17952
                    pqsum += params.CR->carr4p[nser]*qqn;  //Z=17953
                    const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=17954
                    if ( delser<0.0001 ) break; /* goto 81; */  //Z=17955
                    oldpqsum = pqsum;  //Z=17956
                }/*6*/  //Z=17957
                /*81:*/  //Z=17958
                //qDebug() << "formpq A" << q << 0.7*limq4 << "=" << pqsum << nser;
                /*formpq:=*/ return pqsum;  //Z=17959
            }/*5*/  //Z=17960
            else
            {
                //if ( q < 0.75*limq4 )
                //    qDebug() << "formpq B" << q << 0.7*limq4 << "=" << por/(q*q*q*q) << "por" << por;
                /*formpq:=*/ return params.por/(q*q*q*q);  //Z=17961
            }
        }/*4*/ /*  of homogeneous isotropic cube */  //Z=17962

#ifdef procnotused
        /*  core/shell isotropic cube  */  //Z=17964
        if ( params.cs==1 )
        {/*4*/  //Z=17965
            /*formpq:=*/ return polycscube(1.0,params.rho,params.p1,1.0,0.001,0.0001,2*params.radiusi,0,params.sigma,q);  //Z=17966
        }/*4*/  //Z=17967
#endif
    }/*3*/  /*  of isotropic cube  */  //Z=17968

    /*  perfectly oriented cube  */  //Z=17970
    if ( ordis==6 )
    {/*3*/  //Z=17971
        /* if (orcase=1) then begin  //Z=17972 */
        //if ( 1==1 )
        //{/*4*/  //Z=17973

        if ( q<(1.0/params.radius) )
        {/*5*/  //Z=17975
            pql = 1.0;  //Z=17976
            double oldpqsum = 0.0;  //Z=17977
            qxn[0] = 1.0;  //Z=17978
            qyn[0] = 1.0;  //Z=17979
            for ( int nser=1; nser<=80; nser++ )
            {/*6*/  //Z=17980
                qxn[nser] = qxn[nser-1]*qxs*qxs;  //Z=17981
                qyn[nser] = qyn[nser-1]*qys*qys;  //Z=17982

                double binsum = 0.0;  //Z=17984
                for ( int mser=0; mser<=nser; mser++ )
                {/*7*/  //Z=17985
                    /* indx:=mser+1+round(nser*(nser+1)/2);  //Z=17986 */
                    /* binsum:=binsum+carr1pm[indx]*qxn[nser-mser]*qyn[mser];  //Z=17987 */
                    binsum += params.CR->carr11pm[mser][nser-mser]*qxn[nser-mser]*qyn[mser];  //Z=17988
                }/*7*/  //Z=17989
                pql += binsum;  //Z=17990
                const double delser = fabs((pql-oldpqsum)/pql);  //Z=17991
                if ( delser<0.0001 ) break; /* goto 84; */  //Z=17992
                oldpqsum = pql;  //Z=17993
            }/*6*/  //Z=17994
            /*84:*/  //Z=17995
            //pql = pqsum;  //Z=17996
        }/*5*/  //Z=17997
        else
        {/*5*/  //Z=17998
            const double argqx = qxs*params.radius/(zr+1)+eps9;  //Z=17999
            const double argqy = qys*params.radius/(zr+1)+eps9;  //Z=18000
            const double pqrx = (1/(2.0*zr*(zr-1)))*(1/(argqx*argqx))*(1-cos((zr-1)*atan(2.0*argqx))/pow(1.0+4*argqx*argqx,(zr-1)/2.0));  //Z=18001
            const double pqry = (1/(2.0*zr*(zr-1)))*(1/(argqy*argqy))*(1-cos((zr-1)*atan(2.0*argqy))/pow(1.0+4*argqy*argqy,(zr-1)/2.0));  //Z=18002
            pql = pqrx*pqry;  //Z=18003
        }/*5*/  //Z=18004
        /*formpq:=*/ return pql;  //Z=18005
        //}/*4*/  /*  of orcase=1  */  //Z=18006
    }/*3*/  /*  of perfect cube  */  //Z=18007
    //}/*2*/  /*  of cube  */  //Z=18008

    return 0.0;
}/*1*/  //Z=18387


#endif // SC_LIB_FORMPQ_partCube_H
