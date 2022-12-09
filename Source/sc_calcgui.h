#ifndef SC_CALCGUI_H
#define SC_CALCGUI_H

#include <QWidget>
#include <QLabel>
#include <QSlider>
#include <QLineEdit>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QGridLayout>
#include <QVector>
#include <QTabWidget>
#include <QMenu>
#include <QAction>
#include <QStatusBar>
#include <QHash>

#include "sc_calc.h"
//#include "sc_libs_gpu.h"


// Daten zum Vergleich der aktuellen Parameter mit einem File
typedef struct
{
    QString cur, par;
} _CompValues;


// Diese Struktur sollte auch alle Daten ohne GUI enthalten können
typedef struct
{
    enum { /*text,*/ number, select, toggle } type;
    QString key;    // Name des Parameters
    bool fitparam;  // true=kann gefittet werden
    bool disabled;  // true=Parameter ist gesperrt (bleibt aber in der Struktur erhalten)
    QLabel *lbl1;   // in der GUI
    union
    {
        //QLineEdit *inp;         // Text (?)
        QDoubleSpinBox *num;    // Zahl (double, int)
        QComboBox *cbs;         // Auswahl (int)
        QCheckBox *tog;         // Bool
        QWidget *w;             // vereinfacht den Zugriff
    } gui;      // nur in der GUI-Version gefüllt

    struct
    {
        //QString text;
        double  number; // also for selections
        bool    flag;
    } value;    // immer genutzt
    //QString str;
} paramHelper;


class calcHelper
{
public:
    calcHelper( SC_Calc *c, bool gui );
    SC_Calc *subCalc;
    QHash<QString,paramHelper*> params;
    QGridLayout *subLayout;
    static QStringList slDisWidgets;
};


#ifdef undef
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
#endif


class SC_CalcGUI : public QObject
{
    Q_OBJECT
public:
    explicit SC_CalcGUI();
    QStringList getCalcTypes();
    void createTabs( QTabWidget *tab, QStatusBar *sb );

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

    void saveParameter( QString fn );
    QString loadParameter(QString fn, QString onlyMethod);
    void compareParameter( QSettings &sets, QHash<QString,_CompValues*> &compWerte );

    bool setParamToAllMethods( QString par, double val );

    QString index2methodname( int i ) { if ( i<0 || i>=methodNamesSorted.size() ) return "?"; return methodNamesSorted[i]; }

    void prepareCalculation( QString m, bool getData );
    // Globale Inputs für die Berechnungen
    static QHash<QString,Double3> inpVectors;
    static QHash<QString,double>  inpValues;
    static calcHelper *curMethod;
    static QHash<QString,double>  inpSingleValueVectors;

    QStringList paramsForMethod( int im, bool num, bool glob, bool fit );
    QStringList paramsForMethod( QString m, bool num, bool glob, bool fit );

    QString currentParamValueStr( QString m, QString p, bool text );
    double  currentParamValueDbl( QString m, QString p );
    int     currentParamValueInt( QString m, QString p );
    bool    currentParamValueLog( QString m, QString p );

    bool limitsOfParamValue( QString m, QString p, double &min, double &max, bool &countable );
    bool updateParamValue( QString m, QString p, double v, QColor col, bool dbg=false );
    bool updateParamValueForFit( QString p, double v, bool dbg );
    typedef QHash<QString,double> _numericalParams;
    _numericalParams allNumericalParams( QString m );
    bool isCurrentParameterValid( QString m, QString p, bool forfit );
    void resetParamColorMarker( QColor col );

    void doCalculation( int numThreads, progressAndAbort pa );
    double doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt);
    void endThread() { memory->endThread(); }

    typedef enum { htimPrep, htimCalc, htimBoth } whichHigResTimer;
    double higResTimerElapsed( whichHigResTimer f );

    int minX() { return (memory!=nullptr) ? memory->minX() : 0; }
    int maxX() { return (memory!=nullptr) ? memory->maxX() : 0; }
    int minY() { return (memory!=nullptr) ? memory->minY() : 0; }
    int maxY() { return (memory!=nullptr) ? memory->maxY() : 0; }
    double *data() { return (memory!=nullptr) ? memory->data() : nullptr; }
    SC_Calc *getMemPtr() { return memory; }
    bool getLastXY( int &x, int &y );

    paramHelper *getParamPtr( QString m, QString p );

#ifdef COPY_FITDATA_TO_GPU  // FITDATA_IN_GPU ok, real func außen
    bool setArrDataForFit( const double *data ) { return (memory!=nullptr) ? memory->setArrDataForFit(data) : false; }
#ifdef CALC_FQS_IN_GPU
    double getFQS() { return (memory!=nullptr) ? memory->getFQS() : 0.0; }
#endif
#endif
#ifdef FITDATA_IN_GPU  // real func außen
    bool setFitData( int sx, int sy, const double *data )
    { return (memory!=nullptr) ? memory->setFitData(sx,sy,data) : false; }
#endif
    void setNoFitRect( int id, int x0, int y0, int x1, int y1 );

private slots:
    void customContextMenuRequested(const QPoint &pos);

private:
    QHash<QString,calcHelper*> methods;
    QStringList methodNamesSorted;

    //SC_GpuMemory *memory;
    SC_Calc *memory;
    //SC_Libs *libs;

    static bool dataGetter( QString p, _valueTypes &v );
    static bool dataGetterForFit( QString p, _valueTypes &v );

    QMenu   *popupMenu;
    QAction *actionSetDefault,
            *actionCopyToAll;

    QStatusBar *mainStatusBar;
};

#endif // SC_CALCGUI_H
