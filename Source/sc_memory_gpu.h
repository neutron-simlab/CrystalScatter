//#ifndef SC_MEMORY_GPU_H
//#define SC_MEMORY_GPU_H

//#include "sc_math.h"
//#ifndef __CUDACC__
//#include <QDebug>
//#endif


//class SC_GpuMemory
//{
//public:
//    SC_GpuMemory();

    void initMemory();

    int dbgHelper[5]; // für latpar1

public:
    void memcleanup( void *arr );

    bool gpuAvailable() { return !noGPUavailable; }

#ifdef __CUDACC__
    __host__ __device__
#endif
    inline size_t IDX( int x, int y ) const
    {
        return (-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y));
    }
    int minX() { return _xmin; }
    int maxX() { return _xmax; }
    int minY() { return _ymin; }
    int maxY() { return _ymax; }
    //double xyIntensity( int x, int y ) { return arrXYIntensity[3 + (-_xmin + (x)) + (_xmax-_xmin)*(-_ymin + (y))]; }
    double xyIntensity( int x, int y ) { return arrXYIntensity[3 + IDX(x,y)]; }
    double *data() { return arrXYIntensity+3; }

    int lastX() { return arrXYIntensity==nullptr ? -10000 : static_cast<int>(arrXYIntensity[0]); }
    int lastY() { return arrXYIntensity==nullptr ? -10000 : static_cast<int>(arrXYIntensity[1]); }
#ifndef __CUDACC__
    bool dbgFlag() const { return arrXYIntensity==nullptr ? false : (static_cast<int>(arrXYIntensity[2]) != 0); }
#endif

#ifdef COPY_FITDATA_TO_GPU
    bool setArrDataForFit( const double *data )
    {
        if ( data == nullptr )
        {
            arrDataForFitUsed = false;
            return false;
        }
        if ( arrDataForFit == nullptr ) return false;
        memcpy( arrDataForFit, data, sizeof(double)*_arrCount );
        arrDataForFitUsed = true;
        return true;
    }

#ifdef CALC_FQS_IN_GPU
    double getFQS()
    {
        if ( !arrDataForFitUsed ) return 0.0;
        double *ptr = arrFqsForFit;
        double sum = 0.0;
        for ( size_t i=0; i<_arrCount; i++ ) sum += *(ptr++);
        return sum;
    }
#endif
#endif

#ifdef FITDATA_IN_GPU  // real func decl
    bool setFitData( int sx, int sy, const double *data );
#endif

private:
    void createMemory( void **ptr, size_t lensoll, size_t &lenist, bool gpuonly, const char *dbgInfo );

    void checkArrays( int minx, int maxx, int miny, int maxy );

#ifdef __CUDACC__
    __host__ __device__
#endif
    void setXYIntensity( int x/*ihex*/, int y/*i*/, double val ) const
    {
        size_t idx = IDX(x,y); // (-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y));
        if ( idx >= _arrCount )
        {
#ifndef __CUDACC__
            qDebug() << "MEMFAIL" << "x"<<x << "y"<<y << "idx"<<idx << "cnt"<<_arrCount
                     << "xmin"<<_xmin << "xmax"<<_xmax << "ymin"<<_ymin << "ymax"<<_ymax;
#endif
            return;
        }
        arrXYIntensity[idx+3] = val;
        arrXYIntensity[0] = x;      // Merker für letztes gerechnetes Pixel, wird vom Loggig-Thread
        arrXYIntensity[1] = y;      //  genutzt, um den Fortschritt anzuzeigen.
        //arrXYIntensity[2] ist ein Debug-Flag

#ifdef COPY_FITDATA_TO_GPU
#ifdef CALC_FQS_IN_GPU
        if ( arrDataForFitUsed )
        {   // Spezielle Aktionen (Maskieren und FQS) für das Simplex-2D-Fit
            // Diese Berechnung läuft hier, damit die CPU dann nur noch aufsummieren muss
            if ( val > 0 && arrDataForFit[idx] > 0 )
            {
                arrFqsForFit[idx] = FQSVERGL( val,                     // berechnetes Pixel
                                              arrDataForFit[idx] );    // anzufittendes Pixel
            }
            else
            {
                arrFqsForFit[idx] = 0.0;
            }
        }
#endif
#endif
    }

#ifndef __CUDACC__
    void setDebugIntensity( bool val ) const
    {
        if ( arrXYIntensity == nullptr ) return;
        arrXYIntensity[2] = val ? 1 : 0;
    }
#endif

//private:
    bool   noGPUavailable;

    // Intensity-Array
    // xyintensity: array[-imax..imax,-imax..imax] of real;
    double *arrXYIntensity;     //!< array for the resulatant image
    int _xmin, _xmax,           //!< array limits for the z-values (display horizontal)
        _ymin, _ymax;           //!< array limits for the i-values (display vertical)
    size_t _arrCount;           //!< number of double values in the allocated memory
    size_t arrXYsize;

#ifdef COPY_FITDATA_TO_GPU
    double *arrDataForFit;      //!< array for fit data (same dim as arrXYIntensity), [0]=FQS
    size_t arrDataSize;
#ifdef CALC_FQS_IN_GPU
    double *arrFqsForFit;       //!< array for calculated fqs (same dim as arrXYIntensity), [0]=FQS
    size_t arrFqsSize;
#endif
    bool    arrDataForFitUsed;  //!< flag set if this feature will be used
#endif

#ifdef FITDATA_IN_GPU
    double *arrFitData;
    double *arrFitFqs;
    int    _fitWidth,
           _fitHeight;
    size_t _arrFitSize;
    bool   _fitEnabled;
    size_t arrFitDSize, arrFitFSize;
#endif

//};

//#endif // SC_MEMORY_GPU_H
