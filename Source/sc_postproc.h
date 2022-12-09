#ifndef SASCALC_POSTPROC_H
#define SASCALC_POSTPROC_H

#include <QDebug>


//  http://paulbourke.net/miscellaneous/dft/
//#define USE_COMPLEX_WEB - im .pro File wenn fftw3 nicht gefunden wird
#ifdef USE_COMPLEX_WEB
#define TRUE  1
#define FALSE 0
typedef struct
{
    double real, imag;
} COMPLEX;
int FFT2D(COMPLEX **c,int nx,int ny,int dir);
int FFT(int dir,int m,double *x,double *y);
int Powerof2(int n,int *m,int *twopm);
#endif


class SasCalc_PostProc
{
public:
    static SasCalc_PostProc *inst();

    void setLogging( bool f ) { _doLogOutput = f; }

    typedef enum { outReal, outImag, outAbs, outSpec } _outType;

    double *generateRphi(int x0, int x1, int y0, int y1,    // Source size
                         int xbs, int ybs,                  // Beamstop (=Center)
                         int s,                             // Dest. size
                         const double *d,                   // Input data
                         bool scale, bool clip , bool clip40);  // Scaling flags
    double *calculateIFFT(bool foreward,
                          int x0, int x1, int y0, int y1,   // Source size
                          const double *d,                  // Source data
                          int s,                            // Dest. size
                          _outType ot,                      // Output flag
                          bool scaled, bool clip, bool clip40, // Scaling flags
                          bool swapout );                   // Swap corners

    void scaleAndClipData(double *data/*In+Out*/, int len, bool scale, bool clip, bool clip40, bool genlog);

private:
    SasCalc_PostProc();
    static SasCalc_PostProc *_instancePtr;
    bool _doLogOutput;

    // Inputimage
    const double *_dataPtr;
    int _xmin, _xmax, _ymin, _ymax;
    int _len;
    double data( int x, int y )
    {
        if ( _dataPtr == nullptr ) return 0.0;
        //return _dataPtr[(-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y))];
        int idx = (-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y));
#ifdef CONSOLENPROG
        if ( idx < 0 || idx >= _len ) { return 0.0; }
#else
        if ( idx < 0 || idx >= _len ) { qDebug()<<"PostProc:DATA"<<x<<y<<"idx"<<idx<<"len"<<_len; return -1.0; }
#endif
        return _dataPtr[idx];
    }
#ifdef undef
    double datafft( int x, int y )
    {
        if ( _arrFFT == nullptr ) return 0.0;
        //return _dataPtr[(-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y))];
        int idx = (-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y));
        if ( idx < 0 || idx >= _len ) { qDebug()<<"PostProc:DATAFFT"<<x<<y<<"idx"<<idx<<"len"<<_len; return -1.0; }
        return _arrFFT[idx];
    }
#endif

    // Data for the FFT
    double *_arrFFT;
    int     _arrFFTdim;
#ifdef USE_COMPLEX_WEB
    COMPLEX **_complex;
#endif

    // Arrays for the (r,phi) calculation
    double *_arrPol;    // size = _arrPolDim * _arrPolDim
    int     _arrPolDim;

    // Data for the iFFT
    double *_arrIFFT;
    int calculateIndex( bool swap, int x, int y, int dim );

};

#endif // SASCALC_POSTPROC_H
