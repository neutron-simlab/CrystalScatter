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
#include "myguiparam.h"


#define SETCOL(w,c) setColHelper::setcol(w,c)
#define SETCOLMARK_TRAINTBL  Qt::yellow         // on_butTPVreadTrainingTable_clicked()
#define SETCOLMARK_PARDIFF   Qt::green          // compareParameter()
#define SETCOLMARK_CHGIMG    Qt::cyan           // on_actionFind_parameters_changing_image_triggered
#define SETCOLMARK_OLDPARAM  QColor(150,220,255)// mark parameters not set by reading old parameter files
#define SETCOLMARK_NOTOOLTIP QColor(250,150,100)// mark parameters without tool tips
#define SETCOLMARK_IGNORED   Qt::black          // updateParamValue(?) ignored
#define SETCOLMARK_CLEARED   Qt::white          // clear color mark

class setColHelper
{
public:
    setColHelper() {}
    static void setcol(QWidget *w, QColor c);
};

// Daten zum Vergleich der aktuellen Parameter mit einem File
typedef struct
{
    QString cur, par;
} _CompValues;

/*
// Diese Struktur sollte auch alle Daten ohne GUI enthalten können
typedef struct
{
    enum { undef, numdbl, numint, select, toggle, outdbl } type;
    QString key;        // Name des Parameters
    QString tooltip;    // Default ToolTip ohne Fitrange oder andere Zusatzinfos (aus der Liste)
    myGuiParam *gpFit;  // !=0 wenn ein Fit möglich ist (enabled und checked ==> kommt in die Tabelle)
    union
    {
        QDoubleSpinBox *numd;   // Zahl (double)
        QSpinBox  *numi;        // Zahl (int)
        QComboBox *cbs;         // Auswahl (int)
        QCheckBox *tog;         // Bool
        QRadioButton *rad;      // Bool, nur für die spezielle Auswahl von Qmax intern verwendet
        QLineEdit *out;         // Für Outputs (ReadOnly)
        QWidget *w;             // vereinfacht den Zugriff
    } gui;              // nur in der GUI-Version gefüllt
    bool  enabled;
    struct
    {
        double  number; // also for selections
        bool    flag;
    } value;    // immer genutzt
} paramHelper;
*/

class SC_CalcGUI : public QObject
{
    Q_OBJECT
public:
    explicit SC_CalcGUI();

    void cleanup() { if ( calcGeneric != nullptr ) calcGeneric->cleanup(); }

    bool gpuAvailable()
    {
        if ( calcGeneric == nullptr ) return false;
        return calcGeneric->gpuAvailable();
    }

    void saveParameter( QString fn );
    QString loadParameter(QString fn, QString onlyMethod, bool &hkl, bool &grid, bool logload);
    void compareParameter( QSettings &sets, QHash<QString,_CompValues*> &compWerte, QStringList tbign );
    bool    loadparamWithLog( QSettings &sets, QString key, bool def,    bool dolog, bool oldkey=false );
    int     loadparamWithLog( QSettings &sets, QString key, int def,     bool dolog, bool oldkey=false );
    double  loadparamWithLog( QSettings &sets, QString key, double def,  bool dolog, bool oldkey=false );
    QString loadparamWithLog( QSettings &sets, QString key, QString def, bool dolog, bool oldkey=false );

    QString loadparamWithLog( QSettings &sets, QString key, const char *def, bool dolog, bool oldkey=false )
    { return loadparamWithLog(sets,key,QString(def),dolog,oldkey); }

    void loadParamLogmessage(QSettings &sets, QString key, QString msg );

    void saveFitFlags( QSettings &sets );
    void loadFitFlags( QSettings &sets );

    void prepareCalculation( bool fromFit, bool only1d );
    // Globale Inputs für die Berechnungen
    static QHash<QString,Double3> inpVectors;
    static QHash<QString,double>  inpValues;
    static QHash<QString,double>  inpSingleValueVectors;
    static QHash<QString,paramHelper*> params;
    void updateOutputData();

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
    bool isCurrentParameterVisible(QString p, QString &dbg);

    void doCalculation(int numThreads, bool bIgnoreNewSwitch);
    double doFitCalculation(int numThreads, int bstop, int border, long &cnt, long &nancnt);
    void endThread() { calcGeneric->endThread(); }

    typedef enum { htimPrep, htimCalc, htimBoth } whichHigResTimer;
    double higResTimerElapsed( whichHigResTimer f );

    int minX() { return (calcGeneric!=nullptr) ? calcGeneric->minX() : 0; }
    int maxX() { return (calcGeneric!=nullptr) ? calcGeneric->maxX() : 0; }
    int minY() { return (calcGeneric!=nullptr) ? calcGeneric->minY() : 0; }
    int maxY() { return (calcGeneric!=nullptr) ? calcGeneric->maxY() : 0; }
    double *data() { return (calcGeneric!=nullptr) ? calcGeneric->data() : nullptr; }
    SC_Calc_GENERIC *getCalcPtr() { return calcGeneric; }
    bool getLastXY( int &x, int &y );

    paramHelper *getParamPtr( QString p );

#ifdef FITDATA_IN_GPU  // real func außen
    bool setFitData( int sx, int sy, const double *data )
    { return (calcGeneric!=nullptr) ? calcGeneric->setFitData(sx,sy,data) : false; }
#endif
    void setNoFitRect( int id, int x0, int y0, int x1, int y1 );

    void updateToolTipForCompare( QWidget *w, QString txt );

private:
    SC_Calc_GENERIC *calcGeneric;

    static void dataGetter( QString p, _valueTypes &v );
    static void dataGetterForFit( QString p, _valueTypes &v );
    static void dataSetter( QString p, _valueTypes &v );
};

#endif // SC_CALCGUI_H
