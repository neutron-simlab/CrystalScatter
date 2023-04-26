#ifndef SC_CALCGUI_H
#define SC_CALCGUI_H

#include <QWidget>
#include <QLabel>
#include <QSlider>
#include <QLineEdit>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QRadioButton>
#include <QGridLayout>
#include <QVector>
#include <QTabWidget>
#include <QMenu>
#include <QAction>
#include <QStatusBar>
#include <QHash>
#include <QSettings>

#include "sc_calc_generic.h"


#define SETCOL(w,c) setColHelper::setcol(w,c)
#define SETCOLMARK_TRAINTBL  Qt::yellow
#define SETCOLMARK_PARDIFF   Qt::green
#define SETCOLMARK_CHGIMG    Qt::cyan
#define SETCOLMARK_IGNORED   Qt::black
#define SETCOLMARK_CLEARED   Qt::white

class setColHelper
{
public:
    setColHelper() {}

    static void setcol(QWidget *w, QColor c)
    {
        QPalette pal=w->palette();
        if ( w->objectName().startsWith("cbs") )
        {
            if ( c == SETCOLMARK_CLEARED )
                pal.setBrush(QPalette::Text,Qt::black);
            else
                pal.setBrush(QPalette::Text,c.darker(150));
        }
        else if ( w->objectName().startsWith("tog") )
        {
            if ( c == SETCOLMARK_CLEARED )
                pal.setBrush(QPalette::Button,Qt::white); // QColor(240,240,240)); // aus Designer als Default gelesen
            else
                pal.setBrush(QPalette::Button,c.darker(150));
        }
        else // inp
            pal.setBrush(QPalette::Base,c);
        w->setPalette(pal);
    }
};

// Daten zum Vergleich der aktuellen Parameter mit einem File
typedef struct
{
    QString cur, par;
} _CompValues;


// Diese Struktur sollte auch alle Daten ohne GUI enthalten können
typedef struct
{
    enum { undef, numdbl, numint, select, toggle } type;
    QString key;        // Name des Parameters
    QCheckBox *togFit;  // !=0 wenn ein Fit möglich ist (enabled und checked ==> kommt in die Tabelle)
    union
    {
        QDoubleSpinBox *numd;   // Zahl (double)
        QSpinBox  *numi;        // Zahl (int)
        QComboBox *cbs;         // Auswahl (int)
        QCheckBox *tog;         // Bool
        QRadioButton *rad;      // Bool, nur für die spezielle Auswahl von Qmax intern verwendet
        QWidget *w;             // vereinfacht den Zugriff
    } gui;              // nur in der GUI-Version gefüllt

    struct
    {
        double  number; // also for selections
        bool    flag;
    } value;    // immer genutzt
} paramHelper;


class calcHelper
{
public:
    calcHelper( SC_Calc_GENERIC *c );
    SC_Calc_GENERIC *subCalc;
    QHash<QString,paramHelper*> params;
};


class SC_CalcGUI : public QObject
{
    Q_OBJECT
public:
    explicit SC_CalcGUI();
    QStringList getCalcTypes();
    void createTabs(QStatusBar *sb );

    void cleanup() { if ( memory != nullptr ) memory->cleanup(); }

    bool gpuAvailable()
    {
        if ( memory == nullptr )
        {
            if ( curMethod == nullptr ) return false; // noch keine Methode definiert
            return curMethod->subCalc->gpuAvailable();
        }
        return memory->gpuAvailable();
    }

    void saveParameter( QString fn );
    QString loadParameter(QString fn, QString onlyMethod, bool &hkl, bool &grid);
    void compareParameter( QSettings &sets, QHash<QString,_CompValues*> &compWerte );

    void saveFitFlags( QSettings &sets );
    void loadFitFlags( QSettings &sets );

    QString index2methodname() { return curMethod->subCalc->methodName(); /*methodNamesSorted[0];*/ }

    void prepareCalculation( bool fromFit );
    // Globale Inputs für die Berechnungen
    static QHash<QString,Double3> inpVectors;
    static QHash<QString,double>  inpValues;
    static calcHelper *curMethod;
    static QHash<QString,double>  inpSingleValueVectors;

    QStringList paramsForMethod(bool num, bool glob, bool fit );

    QString currentParamValueStr(QString p, bool text );
    double  currentParamValueDbl(QString p );
    int     currentParamValueInt(QString p );

    bool limitsOfParamValue( QString p, double &min, double &max, bool &countable );
    bool updateParamValue( QString p, double v, QColor col, bool dbg=false );
    bool updateParamValueForFit( QString p, double v, bool dbg );
    bool updateParamValueColor( QString p, QColor col );
    typedef QHash<QString,double> _numericalParams;
    _numericalParams allNumericalParams();
    bool isCurrentParameterValid( QString p, bool forfit );
    void resetParamColorMarker( QColor col );

    void doCalculation( int numThreads );
    double doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt);
    void endThread() { memory->endThread(); }

    typedef enum { htimPrep, htimCalc, htimBoth } whichHigResTimer;
    double higResTimerElapsed( whichHigResTimer f );

    int minX() { return (memory!=nullptr) ? memory->minX() : 0; }
    int maxX() { return (memory!=nullptr) ? memory->maxX() : 0; }
    int minY() { return (memory!=nullptr) ? memory->minY() : 0; }
    int maxY() { return (memory!=nullptr) ? memory->maxY() : 0; }
    double *data() { return (memory!=nullptr) ? memory->data() : nullptr; }
    SC_Calc_GENERIC *getMemPtr() { return memory; }
    bool getLastXY( int &x, int &y );

    paramHelper *getParamPtr( QString p );

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

private:
    //QHash<QString,calcHelper*> methods;
    QStringList methodNamesSorted;

    SC_Calc_GENERIC *memory;

    static bool dataGetter( QString p, _valueTypes &v );
    static bool dataGetterForFit( QString p, _valueTypes &v );

    void updateToolTipForCompare( QWidget *w, QString txt );

    QStatusBar *mainStatusBar;
};

#endif // SC_CALCGUI_H
