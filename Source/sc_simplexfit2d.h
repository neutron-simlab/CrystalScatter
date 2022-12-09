#ifndef SASCALC_SIMPLEXFIT2D_H
#define SASCALC_SIMPLEXFIT2D_H

#ifdef CONSOLENPROG
#include "sc_calcCons.h"
#else
#include "sc_calcgui.h"
#endif
#include "sc_globalConfig.h"



//#define USEREPETITIONS
// Wenn nicht definiert, wird die äußere Schleife (Rep) um die eigentliche Fit-Routine ignoriert


// Typedefinitions for 2D-Fit
typedef enum { fitNone, fitNumeric, fitCbs } _fitTypes;
typedef struct
{
    double min, max,    // Limits for the fit
           fitstart,    // Startvalue for the fit
           fitres,      // Resultvalue after fitting
           orgval;      // Original value before the fit wa started
    _fitTypes fitType;  // fitCbs has only integer steps
    bool used;          // true=used for the fit, false=not changed during fit
    bool fitvalid;      // Flag, if .fitres is valid
} _fitLimits;
typedef QHash< QString/*name*/, _fitLimits* > _param2fitval;


class SasCalc_SimplexFit2D
{
public:
#ifdef CONSOLENPROG
    SasCalc_SimplexFit2D( SC_CalcCons *c )
#else
    SasCalc_SimplexFit2D( SC_CalcGUI *c )
#endif
    {
        calc=c;
#ifdef FITDATA_IN_GPU
        useGpuFit = false;
#endif
    } /* Konstruktor */

    inline void setImageInfo( int x0, int x1, int y0, int y1, const double *d )
    {
        intensityForFit = d;
        _xmin = x0;
        _xmax = x1;
        _ymin = y0;
        _ymax = y1;
        _xbs  = (x0+x1)/2; // Beamstop = center
        _ybs  = (y0+y1)/2;
    }

    inline void setImageInfo( int x0, int x1, int y0, int y1, int bsx, int bsy, const double *d )
    {
        intensityForFit = d;
        _xmin = x0;
        _xmax = x1;
        _ymin = y0;
        _ymax = y1;
        _xbs  = bsx;
        _ybs  = bsy;
    }

    void doSimplexFit2D( int numThreads, // Number of Threads used in calculation (if no GPU used)
                         double stp,     // step percentage to change the parameters
                         int maxit,      // maximum number of iterations
                         double ftol,    // Tolerance
                         int borderpix,  // Number of pixels ignored at each border
                         int bstoppix,   // Number of pixels ignored around the beam stop
                         progressLogging pl,    // Logging function
                         _param2fitval *vals,   // Hash der Werte (incl. Used, Min, Max, Cur)
                         QString &retinfo );    // Error-/Info-Message

    double higResTimerElapsed;
    int repetitions;
    int numImgCalc;
    bool aborted;

    double getResiduenPixel( int ihex, int i );

private:
    progressLogging progLogging;

    const double *intensityForFit;
    int _xmin, _xmax, _ymin, _ymax, _xbs, _ybs;
    double fitData( int x, int y )
    {
        if ( intensityForFit == nullptr ) return 0.0;
        return intensityForFit[(-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y))];
    }
#ifdef COPY_FITDATA_TO_GPU
    bool useGpuForMask;
#endif
#ifdef FITDATA_IN_GPU  // glob var def
    bool useGpuFit;
#endif

#ifdef CONSOLENPROG
    SC_CalcCons *calc;
#else
    SC_CalcGUI *calc;
#endif

    double fehlerquadratsumme( int numThreads, double *params, QString &info );
    int _borderPixel, _bstoppixel;

    static const int nap = 31; // Maximum array len
    static const int map = 30; // Max Anzahl Parameter

    int    mmax;            //!< number of parameters (i.e. length of arrays)
    int    numFit;          //!< number of variable parameters (<=mmax)
    bool   enabled[nap];    //!< flags: true=fittable, false=fixed
    QStringList indVec;     //!< index vector: indVec[0] is first fittable parameter, values[indVec[0]] is the value for it
    double values[nap],     //!< vector of all current parameters (not in the calculation)
           minVals[nap],    //!< vector of the minimal values to be fitted
           maxVals[nap];    //!< vector of the maximal values to be fitted
    _param2fitval *parVals;

    double amotry(int numThreads, double p[nap][map], double *y,
                   double *sum, int ihi, double fac , QString &info );
};

#endif // SASCALC_SIMPLEXFIT2D_H
