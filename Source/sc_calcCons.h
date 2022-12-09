#ifndef SC_CALCCONS_H
#define SC_CALCCONS_H

#include <QHash>

#include "sc_calc.h"


// Diese Struktur sollte auch alle Daten ohne GUI enthalten können
typedef struct
{
    enum { /*text,*/ number, select, toggle } type;
    QString key;    // Name des Parameters
    bool fitparam;
    struct
    {
        //QString text;
        double  number; // also for selections
        bool    flag;
    } value;    // immer genutzt
    double minNum, maxNum;  // Bei der GUI in der SpinBox für den Fit gebraucht
    //QString str;
} paramConsHelper;
// Gleiche Struktur wie in cs_calcgui.h (paramHelper) nur hier ohne Qt-GUI-Elemente


class calcConsHelper
{
public:
    calcConsHelper( SC_Calc *c );
    SC_Calc *subCalc;
    QHash<QString,paramConsHelper*> params;
};


class SC_CalcCons
{
public:
    explicit SC_CalcCons();
    QStringList getCalcTypes();

    void cleanup() { if ( memory != nullptr ) memory->cleanup(); }

    bool gpuAvailable()
    {
        if ( memory == nullptr )
        {
            if ( methods.isEmpty() ) return false; // noch keine Methode definiert
            return methods[ methods.keys().first() ]->subCalc->gpuAvailable();
        }
        return memory->gpuAvailable();
    }

    void loadParameter( QString fn );
    void saveParameter( QString fn );

    void prepareCalculation( QString m, bool getData );
    // Globale Inputs für die Berechnungen
    static QHash<QString,Double3> inpVectors;
    static QHash<QString,double>  inpValues;
    static calcConsHelper *curMethod;
    static QHash<QString,double>  inpSingleValueVectors;

    QStringList paramsForMethod( QString m, bool num, bool glob, bool fit );
    double currentParamValue( QString m, QString p );
    bool limitsOfParamValue( QString m, QString p, double &min, double &max, bool &countable );
    bool isCurrentParameterValid( QString m, QString p );

    bool updateParamValue( QString m, QString p, double v );
    bool updateParamValueForFit( QString p, double v, bool /*dbg*/ ) { return updateParamValue("",p,v); }

    void doCalculation( int numThreads, progressAndAbort pa );
    double doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt);
    typedef enum { htimPrep, htimCalc, htimBoth } whichHigResTimer;
    double higResTimerElapsed( whichHigResTimer f );

    std::string tpvPerformRandom( std::list<std::string> ids )
    {
        if ( memory!=nullptr )
            return memory->tpvPerformRandom(ids);
        return "TPV: Error (no calc class)";
    }

    int minX() { return (memory!=nullptr) ? memory->minX() : false; }
    int maxX() { return (memory!=nullptr) ? memory->maxX() : false; }
    int minY() { return (memory!=nullptr) ? memory->minY() : false; }
    int maxY() { return (memory!=nullptr) ? memory->maxY() : false; }
    double *data() { return (memory!=nullptr) ? memory->data() : nullptr; }
    SC_Calc *getMemPtr() { return memory; }

#ifdef COPY_FITDATA_TO_GPU
    bool setArrDataForFit( const double *data ) { return (memory!=nullptr) ? memory->setArrDataForFit(data) : false; }
#ifdef CALC_FQS_IN_GPU
    double getFQS() { return (memory!=nullptr) ? memory->getFQS() : 0.0; }
#endif
#endif
#ifdef FITDATA_IN_GPU  // real func außen
    bool setFitData( int sx, int sy, const double *data )
    { return (memory!=nullptr) ? memory->setFitData(sx,sy,data) : false; }
#endif

private:
    QHash<QString,calcConsHelper*> methods;

    SC_Calc *memory;

    static bool dataGetter( QString p, _valueTypes &v );

};

#endif // SC_CALCCONS_H
