#ifndef SASMAIN_H
#define SASMAIN_H

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
#include <QShowEvent>
#include <QClipboard>
#include "sc_calcgui.h"
#include "sc_simplexfit2d.h"
#include "widimage.h"
#include "sc_readdata.h"
#ifndef ChatbotDisabled
#include <QJsonObject>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#endif

// Zeilendefinitionen für ui->tblLattice3D
#define tblLattWaveLen  0
#define tblLattSaDist   1
#define tblLattPixelX   2
#define tblLattPixelY   3
#define tblLattNbRows   4
#define tblLattNbCols   5
#define tblLattBeamX    6
#define tblLattBeamY    7



/**
 * @brief The tblCheckBox class
 * Normal QCheckBox used in QTable-cells (2d fit)
 */
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


/**
 * @brief The myCalcThread class
 * Main calculation thread. All parameters must be set before into the calculation class variables.
 * This thread is implemented to prevent the gui from freezing.
 */
class myCalcThread : public QThread
{
public:
    myCalcThread( SC_CalcGUI *cg, QObject *parent=nullptr );
    void run();
    void beenden() { calcGuiThread->endThread(); /*_exit=true;*/ }
    void setThreads( int n ) { numThreads=n; }
    void setIgnoreNewSwitch( bool f ) { bIgnoreNewSwitch=f; }
private:
    int numThreads;
    bool bIgnoreNewSwitch;
    SC_CalcGUI *calcGuiThread;  // myCalcThread
    //bool _exit;
};


class dlgConfigAutoFit;


QT_BEGIN_NAMESPACE
namespace Ui { class SC_MainGUI; }
QT_END_NAMESPACE

/**
 * @brief The SC_MainGUI class
 * Main gui application window.
 */
class SC_MainGUI : public QMainWindow
{
    Q_OBJECT

public:
    SC_MainGUI(QWidget *parent = nullptr);
    ~SC_MainGUI();

    static SC_MainGUI *getInst() { return current; }
    SC_CalcGUI *getCalcGui() { return calcGui; }

    void updateLogList( QString msg );

    bool is1Dused() { return _is1Dused; }

protected:
    void closeEvent(QCloseEvent *event);

private slots:
    void initCbsAfterStart();
    void imageWindowClosed(QObject*);
    void on_actionExit_triggered();
    void on_tabMain_currentChanged(int index);
    void on_butCalc_clicked();      // mit speichern der Parameter im Temp-File
    void local_butCalc_clicked(bool doshow=true);   // ohne speichern der Parameter im Temp-File (wg Timer als Slot)
    void on_butAbbruch_clicked();
    void on_radNewImageCfg_toggled(bool checked);
    void on_radLastImageCfg_toggled(bool checked);
    void on_radNewImageCal_toggled(bool checked);
    void on_radLastImageCal_toggled(bool checked);

    void on_actionSave_all_Parameters_triggered();
    void on_actionLoad_all_Parameters_triggered();
    void on_butLoadParams_clicked() { on_actionLoad_all_Parameters_triggered(); }
    void on_butSaveParams_clicked() { on_actionSave_all_Parameters_triggered(); }

    void on_lisDataWindows_currentTextChanged(const QString &currentText);
    void on_togIgnoreUpdates_toggled(bool checked);
    void on_butIFFT_clicked();
    void on_butFFT_clicked();
    void on_intGridPoints_valueChanged(int arg1);
    void on_radQ1_toggled(bool checked);
    void on_radQ2_toggled(bool checked);
    void on_radQ4_toggled(bool checked);
    void on_tblHeaderData_cellChanged(int /*row*/, int /*column*/);
    void on_butOpenMeasFile_clicked();
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
    void on_actionTest_10_Calculate_triggered();
    void aiBackProg_error(QProcess::ProcessError);
    void aiBackProg_finished(int,QProcess::ExitStatus);
    void aiBackProg_readyRead();

    void logThreadTimer();

    void on_togFitUseMask_toggled(bool checked);
    void on_butFitUseResult_clicked();
    void on_butFitAutomatic_clicked();
    void fitAutomaticCreateFile(dlgConfigAutoFit*);

    void autoProcessingTimer();

    void on_butClearEditAutoFit_clicked();
    void on_butShowResiduen_clicked();
    void on_butFitCurSave_clicked();
    void on_butFitCurLoad_clicked();
    void on_butFitHistoShow_clicked();
    void on_butFitHistoClear_clicked();
    void on_butDataCopyScaling_clicked();
    void on_actionTest_Compare_data_files_triggered();

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
    void on_butTPVreadTrainingTable_clicked();
    void on_butTPVupdateAllFiles_clicked();

    void on_togFFTscaleRphi_toggled(bool checked);
    void on_togFFTclipRphi_toggled(bool checked);
    void on_togFFTscaleOut_toggled(bool checked);
    void on_togFFTclipOut_toggled(bool checked);
    void on_butFFTverifyRphi_clicked();

    void on_actionFind_parameters_changing_image_triggered();
    void on_butResetColorMarker_clicked();
    void on_butDataSearchDoit_clicked();
    void on_actionLoad_last_calculaition_parameters_triggered();

    void on_cbsLType_currentIndexChanged(int index);
    void on_cbsComboBoxParticle_currentIndexChanged(int index);
    void on_cbsOrdis_currentIndexChanged(int index);
    void on_cbsComboBoxInterior_currentIndexChanged(int index);
    void on_cbsComboBoxPeak_currentIndexChanged(int index);

    void on_actionHide_unused_values_triggered(bool checked);
    void on_butCopyHdrdata2Params_clicked();

    void on_actionAbout_triggered();
    void on_actionAbout_Qt_triggered();

    void on_actionStart_autoprocessing_file_triggered();

#ifndef ChatbotDisabled
    void on_butChatbotSearch_clicked();
    void on_butChatbotStart_clicked();
    void on_butChatbotStop_clicked();
    void chatbotBackProg_error(QProcess::ProcessError error);
    void chatbotBackProg_finished(int exitCode, QProcess::ExitStatus exitStatus);
    void chatbotBackProg_readyRead();
    void on_butChatbotLogSearch_clicked();
    void on_butChatbotTrainfileSearch_clicked();
    void on_butChatbotSaveConfig_clicked();
    void on_butChatbotReadClipboard_clicked();
    void chatbotClipboardChanged(QClipboard::Mode);
    void on_butParamSearchChatbot_clicked();
    void on_butChatbotSrvConnect_clicked();
    void on_butChatbotSrvDisconn_clicked();
    void on_tabCBCommSel_currentChanged(int index);
    void on_inpChatbotTrainfile_textEdited(const QString &arg1);
    void on_inpChatbotSrvInput_textChanged();
    void on_inpChatbotSrvInput_returnPressed();
    void on_togChatbotCalcMain_toggled(bool arg1);
    void on_butChatbotSrvClear_clicked();
    void on_butChatbotSrvSend_clicked();

    // Rest-Api
    void chatbotDoSend(QString cmd);
    void chatbotSlotError(QNetworkReply::NetworkError);
    void chatbotReplyFinished(QNetworkReply*);
#endif

    void on_togUseAdaptiveStep_toggled(bool checked);
    void on_tab1D2D_currentChanged(int index);

    void on_actionSave_current_values_as_text_triggered();
    void on_actionShow_Parameter_w_o_Tooltip_triggered();

    void on_actionGenerate_all_combinations_of_CBs_triggered();

private:
    Ui::SC_MainGUI *ui;

    SC_CalcGUI *calcGui = nullptr;
    bool bIgnoreUpdates;    // GUI updates bei performTimingTests oder per Toggle von der GUI
    bool bNoUpdateFromCBS;  // in den CBS-Callbacks werden keine anderen Daten modifiziert (wegen LoadParameter)
    bool bIgnoreRecalc;     // Load Parameters etc.

    QString fnTempParamFile;

    QFile *autoProcessingFile;
    void openAutoProcessingFile( QString fn );
    void closeAutoProcessingFile( QString reason );

    void fillDataFromInputs();
    void prepareGuiAndData();
    void finishCalculation(bool showtime);

    void performIFFTcalculations(int curimg, QString tit, bool foreward);

    void performOneFitLoop();
    void performFirstFitLoop( widImage *fitImg );

    void searchParameterFilesHelper( QString d, int bl, QString msk );
    bool bSearchParamsRunning;

    QVector<widImage*> images;
    int numberCounter;
    int imgPosX, imgPosX0, imgPosY, imgPosY0;
    QString dataPath;       // Save path to last image
    QString lastDataFile;   // das letzte geöffnete Messdatenfile
    bool closeMainActive;   // True if the main window closes all image windows due to exit
    widImage* addImage( bool neu, int x0, int x1, int y0, int y1, double *d, QString title, bool meta, bool doshow=true );
    widImage *lastUsedImage;        // das letzte (aktuellste) gerechnete Bild
    widImage *curFitImage;          // das anzufittende Bild
    widImage *lastSelectedImage;    // das letzte selektierte Bild in der Liste der Bilder
    int noFitCurRect;
    bool _is1Dused;         // set when Calculate starts (visibility of 1D tab)

    _param2fitval *fitparams;
    SasCalc_SimplexFit2D *fitClass;
    bool oneFitParamUsed, fitIsRunning, fitMinimalOutput;
    double timeForAll;
    int loopsForAll;
    int imgGenForAll;
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

    void performTimingTests(QString out);
    bool timingTestsRunning;

    // Daten für die History (Verlauf der Werteänderungen beim Fit)
#define UnusedValue 1e10    // Dieser Wert in der Liste bedeutet, dass er nicht berechnet wurde
#ifdef UnusedValue
    QHash< QString/*param*/, QVector<double>/*Werte*/ > fitHistory;
    void printHistory( QFile *ftex );
#endif

    static SC_MainGUI *current;
    static bool _bAbbruch, _bAbortTimer;
    void updateProgBar( int val );

    static widImage *myAddImage( int x0, int x1, int y0, int y1, double *d, QString title )
    {   return current->addImage( true, x0, x1, y0, y1, d, title, false );
    }

    QTimer *_logThreadTimer;
    int lastX, lastY, lastPrz, zzmin, zzrange, iimin, iirange, lastTime;
    QElapsedTimer calcRunTime;
    int calcMaxTime;
    myCalcThread *_calcThread;
    void automaticRecalcDoit(); // Wird von cbsLType aufgerufen

    // Tab "AI"
    typedef QMap<QString/*Variable*/, double/*value*/> _loopVariables;
    typedef QMultiMap<QString/*CLASS*/, _loopVariables > _loopDefinition;
    //typedef QHash<QString/*Variable*/, double/*value*/> _loopVariable;
    //typedef QMultiHash<QString/*CLASS*/, _loopVariables > _loopDefinition;
    //typedef QHash<QString/*Variable*/, double/*value*/> _loopVariable;
    //typedef QMultiHash<QString/*CLASS*/, _loopVariables > _loopDefinition;
    _loopDefinition getLoopDefinition();
    _loopDefinition getTPVLoopDefinition();
    bool globFlagTPV;
    double evaluateFormula(QString formel );
    _loopVariables calcval;  // Speicher für alle Variablen / Werte für einen Rechenschritt
    bool performSaveAIOperation( QString fn, bool doOverwrite, bool interactive, bool useTPV );
    void performStartAI_TPV( bool useTPV );
    void performSaveParamOperation( QString fn );
    QProcess *aiBackProg;
    void aiBackProgAddLog( QString );
    bool localRemoveDirectory( QString fn );
    QString getConsoleExecutable(bool domsg);

    // other
    void findBeamCenter( widImage *ikws, int &xerg, int &yerg );
    void setCurrentExtractRect();
    void setCurrentNoFitRect();

    QString local_Load_all_Parameters(QString fn);
    double loadParameter_checkOldFormat(QSettings &sets, QString oldkey, QString newkey, QString &rv, bool logLoad);
    bool local_OpenMeasFile( QString fn, widImage **imgout );
    void copyMetaToLatticeTable( widImage *img );
    bool copyParamsToFitTable();

    void loadScatterParameter( QString fn );
    void checkData( QString data, QString soll );
    void setValue( QString data, QString name );

    QStringList tpvParamsKeys;
    QStringList parChgImg_doCalculation(QStringList inpParams, QLabel *lblinfo);
    bool        parChgImg_compareData(QString par, double val, const double *check, int mx);
    QStringList parChgImg_generateParamList();

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

    void showColorMarker( QColor c );

    bool localCheckToolTip(QWidget *w, bool &used);

    bool genIniFile_CombinationOfCBs(QFile *flog, QString outpath, bool genimg, bool checkNan, bool genJson);

#ifndef ChatbotDisabled
    // Tab "ChatBot"
    QString chatbotConsProg;
    QProcess *chatbotBackProg;
    void chatbotBackProgAddLog(QString msg);
    void chatbotSaveConfig(QString fn, QString descr);
    void chatbotSaveConfigHelper(QJsonObject &jsParam, QString parkey, QString jskey);
    QString jsval2str( QString key, QJsonValue val );
    QHash<QString,QString> paramkey2jsonkey;
    void generateKeyHash();
    void loadChatbotParamFile(QString fn, bool logLoading);

    QNetworkAccessManager *chatbotManager;
    QUrl chatbotUrl;
    QNetworkReply *chatbotReply;
#endif

};

#endif // SASMAIN_H
