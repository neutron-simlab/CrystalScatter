#ifndef SC_LIB_FORMFQ_partCylinder_H
#define SC_LIB_FORMFQ_partCylinder_H


#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formfq_partCylinder( double limql, double qxs, double qys, double q ) const
{  //Z=18395

    const double zl = (1-sqr(params.sigmal))/sqr(params.sigmal);  //Z=18418
    const double zr = (1-sqr(params.sigma))/sqr(params.sigma);  //Z=18419
    const double radiusm = params.radius/params.p1;   /*  outer radius of core/shell particle  */  //Z=18420

    CHECKENDTHREAD_VAL;

    /* ********** */  //Z=18758
    /*  cylinder  */  //Z=18759
    /* ********** */  //Z=18760
    //if ( params.part==1 )
    //{/*2*/  //Z=18761

    double pql, pqr;

    /* ** longitudinal part ** */  //Z=18763
    /* ** isotropic ** */  //Z=18764
    if ( params.ordis==7 )
    {/*3*/  //Z=18765
        if ( q<(0.6*params.limq1f) )
        {/*4*/  //Z=18766
            double pqsum = 1.0;  //Z=18767
            double oldpqsum = 0.0;  //Z=18768
            double qqnn = 1.0;  //Z=18769
            for ( int nser=1; nser<=120; nser++ )
            {/*5*/  //Z=18770
                qqnn = qqnn*q*q;  //Z=18771
                pqsum += params.CR->carr1f[nser]*qqnn;  //Z=18772
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18773
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
    if ( params.ordis==6 )
    {/*3*/  //Z=18789
        if ( params.orcase==4 )
            pql = 1.0;  //Z=18790
        else
        {/*4*/  //Z=18791
            if ( limql<(0.6*params.limq1f) ) //20240301 - war 4
            {/*5*/  //Z=18792
                double pqsum = 1.0;  //Z=18793
                double oldpqsum = 0.0;  //Z=18794
                double qqnn = 1.0;  //Z=18795
                for ( int nser=1; nser<=120; nser++ )
                {/*6*/  //Z=18796
                    qqnn = qqnn*(qxs+qys)*(qxs+qys);  //Z=18797
                    pqsum = pqsum+params.CR->carr1f[nser]*qqnn;  //Z=18798
                    const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18799
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
    if ( params.ordis==0 )
    {/*3*/  //Z=18818
        if ( params.orcase==4 )
            pql = 1.0;  //Z=18819
        else
        {/*4*/  //Z=18820
            if ( limql<(0.2*params.limq1f) ) //20240301 - war 0.7
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
                            binsum += params.CR->carr11pm[mser][nser-mser]*qyn[mser]*qxn[nser-mser];  //Z=18834
                        }/*8*/  //Z=18835
                        pqsum += params.CR->carr1f[nser]*binsum;  //Z=18836
                        const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18837
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
                            binsum += params.CR->carr11pm[nser-mser][mser]*qxn[mser]*qyn[nser-mser];  //Z=18851
                        }/*8*/  //Z=18852
                        pqsum += params.CR->carr1f[nser]*binsum;  //Z=18853
                        const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18854
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
                            binsum += params.CR->carr11pm[mser][nser-mser]*qxn[mser]*qyn[nser-mser];  //Z=18869
                        }/*8*/  //Z=18870
                        pqsum += params.CR->carr1f[nser]*binsum;  //Z=18871
                        const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18872
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
                double qz=1.0; // TODO: nur hier verwendet...
                qrombdeltac(params.p1,params.sigmal,params.alphash1,params.polTheta,0,qxs,qys,qz, // 9,9,
                            9,9,9,9,/*9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,*/params.ordis,1,4,params.orcase,0,0,0,params.CR->carr1f,pql);  //Z=18884
                pql = pql/params.norm;  //Z=18885
            }/*5*/  //Z=18886
        }/*4*/  //Z=18887
    }/*3*/   /*  of general  */  //Z=18888

    /*  transverse part  */  //Z=18890
    /*  homogeneous cylinder  */  //Z=18891
    if ( params.cs==0 )
    {/*3*/  //Z=18892
        if ( q<(1.0*params.limq4f) ) //20240301 - war 1.5
        {/*4*/  //Z=18893
            double pqsum = 1.0;  //Z=18894
            double oldpqsum = 0.0;  //Z=18895
            double qqnn = 1.0;  //Z=18896
            for ( int nser=1; nser<=120; nser++ )
            {/*5*/  //Z=18897
                qqnn = qqnn*q*q;  //Z=18898
                pqsum += params.CR->carr4f[nser]*qqnn;  //Z=18899
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18900
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
                pqsum += params.CR->carr4f[nser]*qqnn;  //Z=18936
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18937
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
                pqsum += params.CR->carr5f[nser]*qqnn;  //Z=18960
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18961
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
                pqsum += params.CR->carr6f[nser]*qqnn;  //Z=18985
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=18986
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

                const double del = fabs((F42-oldF42)/F42);  //Z=19116
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

        const double FF1 = (cc1*F12+cc4*F42+cc6*F62)/sumc;  //Z=19208
        /* FF1:=(cc6*F62)/sumc;  //Z=19209 */

        /*formfq:=*/ return FF1*FF1;  //Z=19212

        /* formfq:=pqcoreshellinf(1.0,rho,p1,1.0,0.001,alfa,radiusm,2,sigmar,q);  //Z=19214 */
    }/*3*/ /*  of inhomogeneous core/shell cylinder */  //Z=19215


    //}/*2*/ /*  of cylinder  */  //Z=19218

    return 0.0;
}  //Z=19829


#endif // SC_LIB_FORMFQ_partCylinder_H
