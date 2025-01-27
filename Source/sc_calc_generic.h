#ifndef SC_CALC_GENERIC_H
#define SC_CALC_GENERIC_H

#include "sc_calc_generic_gpu.h"


/**
 * @brief data transfer helper
 */
typedef struct
{
    bool    checked;
    int     select;
    QString str;
    double  value;
    Double3 vec;
} _valueTypes;

typedef void (*_dataGetter)( QString, _valueTypes& );
typedef void (*_dataSetter)( QString, _valueTypes& );



class SC_Calc_GENERIC
{
public:
    SC_Calc_GENERIC();

    QString methodName() { return "FCC Spheres"; }
    // This name has historical reasons and might be changed in the future (then change all parameter files too)

    QStringList guiLayoutNeu();

    void prepareData( _dataGetter, bool );
    void updateOutputData( _dataSetter );

    void doCalculation( int numThreads, bool bIgnoreNewSwitch )
    { calc->doCalculation( numThreads, bIgnoreNewSwitch ); } //230512//

    double doFitCalculation( int numThreads, int bstop, int border, long &cnt, long &nancnt )
    { return calc->doFitCalculation( numThreads, bstop, border, cnt, nancnt ); }

    std::string tpvPerformRandom(std::list<std::string> ids) { return calc->tpvPerformRandom(ids); }

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
    void scaleIntensity( bool linlog ) { calc->scaleIntensity(linlog); }
    double xyIntensity( int x, int y ) { return calc->xyIntensity(x,y); }
    int lastX() { return calc->lastX(); }
    int lastY() { return calc->lastY(); }

    bool newSwitchUsed() { return calc->newSwitchUsed(); }

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
