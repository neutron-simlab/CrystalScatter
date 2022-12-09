#ifndef SASMAIN_H
#define SASMAIN_H

#ifndef CONSOLENPROG // blendet alles aus

#include <QApplication>
#include <QMainWindow>
#include <QHash>
#include <QTableWidgetItem>
#include <QThread>
#include <QTimer>
#include <QElapsedTimer>
#include <QFile>
#include <QPushButton>
#include <QRadioButton>
#include <QGroupBox>
#include <QProcess>
#include "sc_calcgui.h"
#include "sc_simplexfit2d.h"
#include "widimage.h"
#include "sc_readdata.h"

// Zeilendefinitionen für ui->tblLattice3D
#define tblLattWaveLen  0
#define tblLattSaDist   1
#define tblLattPixelX   2
#define tblLattPixelY   3
#define tblLattNbRows   4
#define tblLattNbCols   5
#define tblLattBeamX    6
#define tblLattBeamY    7



class tblCheckBox : public QWidget
{
public:
    tblCheckBox( bool f ) : QWidget()
    {
        _tog = new QCheckBox();
        _lay = new QHBoxLayout(this);
        _lay->addWidget(_tog);
        _lay->setAlignment(Qt::AlignCenter);
        _lay->setContentsMargins(0,0,0,0);
        _tog->setChecked( f );
    }
    QCheckBox *tog() { return _tog; }
private:
    QCheckBox *_tog;
    QHBoxLayout *_lay;
};


class myCalcThread : public QThread
{
public:
    myCalcThread( SC_CalcGUI *cg, QObject *parent=nullptr );
    void run();
    void beenden() { _exit=true; }
    void setThreads( int n ) { numThreads=n; }
private:
    int numThreads;
    SC_CalcGUI *calcGui;
    bool _exit;
};


QT_BEGIN_NAMESPACE
namespace Ui { class SC_MainGUI; }
QT_END_NAMESPACE

class SC_MainGUI : public QMainWindow
{
    Q_OBJECT

public:
    SC_MainGUI(QWidget *parent = nullptr);
    ~SC_MainGUI();

    static SC_MainGUI *getInst() { return current; }

    void updateLogList( QString msg );

    static QString configParamsFile;

protected:
    void closeEvent(QCloseEvent *event);

private slots:
    void imageWindowClosed(QObject*);
    void on_actionExit_triggered();
    void on_tabMethods_currentChanged(int index);
    void on_tabMain_currentChanged(int index);
    void on_butCalc_clicked();
    void on_butAbbruch_clicked();
    void on_radNewImageCfg_toggled(bool checked);
    void on_radLastImageCfg_toggled(bool checked);
    void on_radNewImageCal_toggled(bool checked);
    void on_radLastImageCal_toggled(bool checked);
    void on_actionSave_all_Parameters_triggered();
    void on_actionLoad_all_Parameters_triggered();
    void on_lisDataWindows_currentTextChanged(const QString &currentText);
    void on_togIgnoreUpdates_toggled(bool checked);
    void on_butIFFT_clicked();
    void on_butFFT_clicked();
    void on_inpGridPoints_valueChanged(int arg1);
    void on_radQ1_toggled(bool checked);
    void on_radQ2_toggled(bool checked);
    void on_radQ4_toggled(bool checked);
    void on_butUseQMax_clicked();
    void on_tblLattice3DValues_cellChanged(int /*row*/, int /*column*/);
    void on_butOpenMeasFile_clicked();
    void on_butSaveAllImages_clicked();
    void on_butDataSetMask_clicked();
    void on_butDataFindCenter_clicked();
    void on_butTestGo_clicked();
    void on_cbsDefaultColTbl_activated(int index);
    void on_togOnlyNewWindow_toggled(bool checked);
    void on_togLimitRuntime_toggled(bool checked);
    void on_tblFitValues_itemChanged(QTableWidgetItem *item);
    void on_butSaveLog_clicked();
    void on_butFitStart_clicked();
    void tblFitUsed_toggled(bool checked);

    // Tab "AI"
    void on_cbsMethod_currentIndexChanged(const QString &arg1);
    void on_butUseForAI_clicked();
    void on_butSaveTable_clicked();
    void on_butLoadTable_clicked();
    void on_butClearTable_clicked();
    void on_butFileName_clicked();
    void on_butSelectDir_clicked();
    void on_butAIcheck_clicked();
    void on_butAIsaveBkgOp_clicked();
    void on_butAIstart_clicked();
    void on_butSaveImagesAI_clicked();
    void on_actionTest_read_first_AI_RUN_Line_triggered();
    void on_actionTest_10_Calculate_triggered();
    void aiBackProg_error(QProcess::ProcessError);
    void aiBackProg_finished(int,QProcess::ExitStatus);
    void aiBackProg_readyRead();

    void logThreadTimer();

    void on_togFitUseMask_toggled(bool checked);
    void on_butFitConfig_clicked();
    void on_butFitUseResult_clicked();
    void on_butFitAutomatic_clicked();

    void autoProcessingTimer();

    void on_butClearEditAutoFit_clicked();
    void on_butShowResiduen_clicked();
    void on_butFitCurSave_clicked();
    void on_butFitCurLoad_clicked();
    void on_butFitHistoShow_clicked();
    void on_butFitHistoClear_clicked();
    void on_butDataCopyScaling_clicked();
    void on_actionTest_Compare_data_files_triggered();
    void on_actionLoad_only_current_parameters_triggered();

    void automaticRecalc();

    void on_actionCompare_current_parameters_with_file_triggered();
    void on_butRemoveVar_clicked();
    void on_tblListe_cellClicked(int row, int column);

    void on_grpExtractImage_toggled(bool arg1);
    void on_inpExtractCenterX_valueChanged(double arg1);
    void on_inpExtractCenterY_valueChanged(double arg1);
    void on_inpExtractGridSize_valueChanged(int arg1);
    void on_butDoExtract_clicked();
    void on_cbsExtractFrameCol_activated(int index);
    void on_cbsExtractScale_activated(int index);

    void on_grpNoFitRegions_toggled(bool arg1);
    void on_tblNoFitRegions_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);
    void on_tblNoFitRegions_cellChanged(int row, int column);

    void extractUpdateRect(QRect);  // vom widImage

    void on_butParamSearchPath_clicked();
    void on_butParamSearchDoit_clicked();
    void on_lisParamSearchResult_itemSelectionChanged();
    void on_butParamSearchGenerate_clicked();
    void on_inpParamSearchPath_textEdited(const QString &arg1);

    void on_butTPVoutPath_clicked();
    void on_butTPVsaveOptions_clicked();
    void on_butTPVloadOptions_clicked();
    void on_butTPVstart_clicked();
    void on_togTPVaddBS_toggled(bool checked);

    void on_togFFTscaleRphi_toggled(bool checked);
    void on_togFFTclipRphi_toggled(bool checked);
    void on_togFFTscaleOut_toggled(bool checked);
    void on_togFFTclipOut_toggled(bool checked);
    void on_butFFTverifyRphi_clicked();

    void on_butReadTrainingTable_clicked();
    void on_actionFind_parameters_changing_image_triggered();

    void on_butResetColorMarker_clicked();

    void on_butDataSearchDoit_clicked();

private:
    Ui::SC_MainGUI *ui;

    SC_CalcGUI *calcGui;
    bool bIgnoreUpdates;    // GUI updates

    bool bIgnoreRecalc;     // Load Parameters etc.

    QString fnTempParamFile;

    QFile *autoProcessingFile;

    void fillDataFromInputs();
    void prepareCalculation( bool getMD, bool progbar );
    void finishCalculation(bool showtime);

    void performIFFTcalculations(int curimg, QString tit, bool foreward);

    void performOneFitLoop();
    void performFirstFitLoop( widImage *fitImg );

    void disableUnusedElements();

    void searchParameterFilesHelper( QString d, int bl, QString msk );
    bool bSearchParamsRunning;

    QVector<widImage*> images;
    int numberCounter;
    int imgPosX, imgPosX0, imgPosY, imgPosY0;
    QString dataPath;  // Save path to last image
    bool closeMainActive;  // True if the main window closes all image windows due to exit
    widImage* addImage( bool neu, int x0, int x1, int y0, int y1, double *d, QString title, bool meta );
    widImage *lastUsedImage;        // das letzte (aktuellste) gerechnete Bild
    widImage *curFitImage;          // das anzufittende Bild
    widImage *lastSelectedImage;    // das letzte selektierte Bild in der Liste der Bilder
    int noFitCurRect;

    QHash< QString/*method*/, _param2fitval* > method2fitparams;
    SasCalc_SimplexFit2D *fitClass;
    bool oneFitParamUsed, fitIsRunning, fitMinimalOutput;
    double timeForAll;
    int loopsForAll;
    int imgGenForAll;
    QString curMethod;
    // public: void updateLogList( QString msg );
    static bool myProgressLogging( char *msg );
    QFile *fileFitLog;
    QString fileFitLogName;
    double fitMeanChangePercent;
    QString lastAutoFitLine;
    QStringList savedFitParams;
    double fitOrgMin, fitOrgmax;
    bool useFixedScaling;
    double minFixedScale, maxFixedScale;

    // Daten für die History (Verlauf der Werteänderungen beim Fit)
#define UnusedValue 1e10    // Dieser Wert in der Liste bedeutet, dass er nicht berechnet wurde
#ifdef UnusedValue
    QHash< QString/*param*/, QVector<double>/*Werte*/ > fitHistory;
    void printHistory( QFile *ftex );
#endif

    static SC_MainGUI *current;
    static bool _bAbbruch;
    void updateProgBar( int val );
    static bool myProgressAndAbort(int);

    static widImage *myAddImage( int x0, int x1, int y0, int y1, double *d, QString title )
    {   return current->addImage( true, x0, x1, y0, y1, d, title, false );
    }

    QTimer *_logThreadTimer;
    int lastX, lastY, lastPrz, zzmin, zzrange, iimin, iirange, lastTime;
    QElapsedTimer calcRunTime;
    int calcMaxTime;
    myCalcThread *_calcThread;
    inline void waitForCalcThread()
    {
        while ( ! _calcThread->wait(50) )
        {
            qApp->processEvents();
            if ( _bAbbruch )
            {
                _calcThread->beenden();
                break;
            }
        }
    }

    // Tab "AI"
    typedef QMap<QString/*Variable*/, double/*value*/> _loopVariables;
    typedef QMultiMap<QString/*CLASS*/, _loopVariables > _loopDefinition;
    _loopDefinition getLoopDefinition();
    _loopDefinition getTPVLoopDefinition();
    bool globFlagTPV;
    double evaluateFormula( QString m, QString formel );
    _loopVariables calcval;  // Speicher für alle Variablen / Werte für einen Rechenschritt
    bool performSaveAIOperation( QString fn, bool doOverwrite, bool interactive, bool useTPV );
    void performStartAI_TPV( bool useTPV );
    void performSaveParamOperation( QString fn );
    QProcess *aiBackProg;
    void aiBackProgAddLog( QString );
    bool localRemoveDirectory( QString fn );

    // other
    void findBeamCenter( widImage *ikws, int &xerg, int &yerg );
    void setCurrentExtractRect();
    void setCurrentNoFitRect();

    QString local_Load_all_Parameters(QString fn, QString onlyMethod);
    bool local_OpenMeasFile( QString fn, widImage **imgout );
    void copyMetaToLatticeTable( widImage *img );
    bool copyParamsToFitTable();

    void loadScatterParameter( QString fn );
    void checkData( QString data, QString soll );
    void setValue( QString data, QString name );

    QStringList tpvParamsKeys;
    bool vergleicheCheckData( QString par, double val, const double *check, int mx );

    // Daten zum Vergleich der aktuellen Parameter mit einem File
    /*typedef struct
    {
        QString cur, par;
    } _CompValues;*/
    QHash<QString/*name*/,_CompValues*> compWerte;
    void compBool( QGroupBox *tog, QSettings &sets, QString key, QString prmt="" );
    void compBool( QCheckBox *tog, QSettings &sets, QString key, QString prmt="" );
    void compBool( QRadioButton *tog, QSettings &sets, QString key, QString prmt="" );
    void compDouble( QDoubleSpinBox *inp, QSettings &sets, QString key );
    void compInt( QSpinBox *inp, QSettings &sets, QString key );
    void compInt( QComboBox *inp, QSettings &sets, QString key, QString prmt );
    void compString( QLineEdit *inp, QSettings &sets, QString key, QString prmt="" );

};
#endif // CONSOLENPROG

#endif // SASMAIN_H
