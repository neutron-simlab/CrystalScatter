#ifndef SC_CALC_GENERIC_H
#define SC_CALC_GENERIC_H

#include "sc_calc.h"
#include "sc_calc_generic_gpu.h"


class SC_Calc_GENERIC : public SC_Calc
{
public:
    SC_Calc_GENERIC();

    QString methodName() { return "FCC Spheres"; }  // (TODO) Wenn hier geändert wird, dann müüssen auch die INI-Files angepasst werden.

    QStringList guiLayout();

    void prepareData( _dataGetter );

    void doCalculation( int numThreads, progressAndAbort pa )
    { calc->doCalculation( numThreads, pa ); }

    double doFitCalculation( int numThreads, int bstop, int border, long &cnt, long &nancnt )
    { return calc->doFitCalculation( numThreads, bstop, border, cnt, nancnt ); }

    std::string tpvPerformRandom( std::list<std::string> ids )
    { return calc->tpvPerformRandom(ids); }

    double higResTimerElapsedPrep() { return calc->higResTimerElapsedPrep; }

    double higResTimerElapsedCalc() { return calc->higResTimerElapsedCalc; }

    void endThread() { calc->endThread(); }

    void cleanup() { calc->cleanup(); }
    bool gpuAvailable() { return calc->gpuAvailable(); }
    int minX() { return calc->minX(); }
    int maxX() { return calc->maxX(); }
    int minY() { return calc->minY(); }
    int maxY() { return calc->maxY(); }
    double *data() { return calc->data(); }
    double xyIntensity( int x, int y ) { return calc->xyIntensity(x,y); }
    int lastX() { return calc->lastX(); }
    int lastY() { return calc->lastY(); }

#ifdef COPY_FITDATA_TO_GPU  // FITDATA_IN_GPU ok, real func
    bool setArrDataForFit( const double *data ) { return calc->setArrDataForFit(data); }
#ifdef CALC_FQS_IN_GPU
    double getFQS() { return calc->getFQS(); }
#endif
#endif
#ifdef FITDATA_IN_GPU  // real func
    bool setFitData( int sx, int sy, const double *data ) { return calc->setFitData(sx,sy,data); }
#endif
    void setNoFitRect( int id, int x0, int y0, int x1, int y1 )
    {
        calc->setNoFitRect( id, x0, y0, x1, y1 );
    }

private:
    SasCalc_GENERIC_calculation *calc;

};

#endif // SC_CALC_GENERIC_H
