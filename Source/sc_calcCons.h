#ifndef SC_CALCCONS_H
#define SC_CALCCONS_H

#include <QHash>

#include "sc_calc_generic.h"


// Diese Struktur sollte auch alle Daten ohne GUI enthalten können
typedef struct
{
    enum { /*undef,*/ numdbl, numint, select, toggle, outdbl } type;
    QString key;    // Name des Parameters
    bool fitparam;
    struct
    {
        double  number; // also for selections
        bool    flag;
    } value;    // immer genutzt
    double minNum, maxNum;      // Bei der GUI in der SpinBox für den Fit gebraucht
    QStringList slSelectValues; // Bei der GUI sind die in der ComboBox
} paramConsHelper;
// Ähnliche Struktur wie in sc_calcgui.h (paramHelper) nur hier ohne Qt-GUI-Elemente


/*
class calcConsHelper
{
public:
    calcConsHelper( SC_Calc_GENERIC *c );
    SC_Calc_GENERIC *subCalc;
    QHash<QString,paramConsHelper*> params;
};*/


class SC_CalcCons
{
public:
    explicit SC_CalcCons();

    void cleanup() { if ( calcGeneric != nullptr ) calcGeneric->cleanup(); }

    bool gpuAvailable(){ if ( calcGeneric != nullptr ) return calcGeneric->gpuAvailable(); else return false; }

    void loadParameter( QString fn );
    void saveParameter( QString fn );

    void prepareData(bool fromFit, bool use1d);
    // Globale Inputs für die Berechnungen
    static QHash<QString,Double3> inpVectors;
    static QHash<QString,double>  inpValues;
    static QHash<QString,double>  inpSingleValueVectors;
    static QHash<QString,paramConsHelper*> params;

    QStringList paramsForMethod( bool num, bool glob, bool fit );
    double currentParamValue( QString p );
    bool limitsOfParamValue( QString p, double &min, double &max, bool &countable );
    bool isCurrentParameterValid( QString p );

    bool updateParamValue( QString p, double v );

    void doCalculation( int numThreads, bool ignNewSwitch );
    double doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt);
    typedef enum { htimPrep, htimCalc, htimBoth } whichHigResTimer;
    double higResTimerElapsed( whichHigResTimer f );

    std::string tpvPerformRandom( std::list<std::string> ids )
    {
        if ( calcGeneric!=nullptr )
            return calcGeneric->tpvPerformRandom(ids);
        return "TPV: Error (no calc class)";
    }

    int minX() { return (calcGeneric!=nullptr) ? calcGeneric->minX() : 0; }
    int maxX() { return (calcGeneric!=nullptr) ? calcGeneric->maxX() : 0; }
    int minY() { return (calcGeneric!=nullptr) ? calcGeneric->minY() : 0; }
    int maxY() { return (calcGeneric!=nullptr) ? calcGeneric->maxY() : 0; }
    double *data() { return (calcGeneric!=nullptr) ? calcGeneric->data() : nullptr; }
    void scaleIntensity(bool linlog) { if ( calcGeneric!=nullptr ) calcGeneric->scaleIntensity(linlog); }

    SC_Calc_GENERIC *getCalcPtrWrapper() { return calcGenericWrapper; }
    SasCalc_GENERIC_calculation *getCalcPtr() { return calcGeneric; }

#ifdef COPY_FITDATA_TO_GPU
    bool setArrDataForFit( const double *data ) { return (memory!=nullptr) ? memory->setArrDataForFit(data) : false; }
#ifdef CALC_FQS_IN_GPU
    double getFQS() { return (memory!=nullptr) ? memory->getFQS() : 0.0; }
#endif
#endif
#ifdef FITDATA_IN_GPU  // real func außen
    bool setFitData( int sx, int sy, const double *data )
    { return (calcGeneric!=nullptr) ? calcGeneric->setFitData(sx,sy,data) : false; }
#endif

private:
    SC_Calc_GENERIC *calcGenericWrapper = nullptr;
    SasCalc_GENERIC_calculation *calcGeneric = nullptr;

    static void dataGetter( QString p, _valueTypes &v );

};

#endif // SC_CALCCONS_H
