#ifndef SC_LIB_FORMFQ_partCube_H
#define SC_LIB_FORMFQ_partCube_H

// wird anscheinend nicht verwendet

#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::formfq_partCube( /*double limql,*/ double qx, double qy, /*double qxs, double qys,*/ double q, int ordis ) const
{  //Z=18395

    CHECKENDTHREAD_VAL;

    /*  cube  */  //Z=19799
    //if ( params.part==5 )
    //{/*2*/  //Z=19800
    /*  homogeneous cube  */  //Z=19801
    if ( params.cs==0 )
    {/*3*/  //Z=19802
        if ( q<0.7*params.limq4f )
        {/*4*/  //Z=19803
            double pqsum = 1.0;  //Z=19804
            double oldpqsum = 0.0;  //Z=19805
            double qqnn = 1.0;  //Z=19806
            for ( int nser=1; nser<=120; nser++ )
            {/*5*/  //Z=19807
                qqnn = qqnn*q*q;  //Z=19808
                pqsum = pqsum+params.CR->carr4f[nser]*qqnn;  //Z=19809
                const double delser = fabs((pqsum-oldpqsum)/pqsum);  //Z=19810
                if ( delser<0.0001 ) break; /* goto 81; */  //Z=19811
                oldpqsum = pqsum;  //Z=19812
            }/*5*/  //Z=19813
            /*81:*/  //Z=19814
            /*formfq:=*/ return pqsum;  //Z=19815
        }/*4*/  //Z=19816
        else
        {/*4*/  //Z=19817
            double pql;
            double qz=1.0;                              // TODO: passt das?
            double qxhklt=0,qyhklt=0,qzhklt=0,qhkl=0;   // TODO: passt das?
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

    //}/*2*/  /*  of cube  */  //Z=19828

    return 0.0;
}  //Z=19829


#endif // SC_LIB_FORMFQ_partCube_H
