#include "sc_maingui.h"
#include "ui_sc_maingui.h"
#include "sc_postproc.h"
#include "dlgtimetests.h"
#include "dlgconfigautofit.h"
#include <QDebug>
#include <QHeaderView>
#include <QFileDialog>
#include <QSettings>
#include <QMessageBox>
#include <QDesktopWidget>
#include <QDate>
#include <QPainter>
#include <thread>


// Da CUDA nicht unter Windows nutzbar ist, versuche ich mal mit OpenCL
// --> https://www.qt.io/blog/2010/04/07/using-opencl-with-qt
// Ist jetzt aber erst nur ein Gedanke (noch nichts programmiert).


SC_MainGUI *SC_MainGUI::current = nullptr;
bool SC_MainGUI::_bAbbruch = false;
QString SC_MainGUI::configParamsFile = "";



SC_MainGUI::SC_MainGUI(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::SC_MainGUI)
{
    QSettings sets(SETT_APP,SETT_GUI);
    int lastTab = sets.value("LastMethod",-1).toInt();

    configParamsFile = sets.value( "ConfigParamFile", "" ).toString();

    _logThreadTimer = nullptr;
    _calcThread = nullptr;
    fitIsRunning = false;
    fitMinimalOutput = false;
    lastAutoFitLine = "";
    bIgnoreRecalc = true; // during init
    lastUsedImage = nullptr;
    curFitImage = nullptr;
    lastSelectedImage = nullptr;
    useFixedScaling = false;
    noFitCurRect = -1;
    autoProcessingFile = nullptr;

    ui->setupUi(this);
    current = this;

    ui->actionTest_read_first_AI_RUN_Line->setVisible(false); // speziell für die Doku wegnehmen
    //ui->actionFind_parameters_changing_image->setEnabled(false);

    /*
    setWindowTitle( "SasCrystal - Version Juni 2021" );
    QLabel *lblVersInfo = new QLabel("New Version - Juni 2021");
    QPalette pal = lblVersInfo->palette();
    pal.setColor( QPalette::WindowText, Qt::red );
    lblVersInfo->setPalette(pal);
    ui->tabMain->setCornerWidget(lblVersInfo);
    */

    // Get the data directory path
    QDir fileDir( qApp->applicationDirPath() );
    while ( fileDir.absolutePath().contains("debug", Qt::CaseInsensitive) ||
            fileDir.absolutePath().contains("release", Qt::CaseInsensitive) )
        fileDir.cdUp();
    fileDir.cd( "data" );
    dataPath = fileDir.absolutePath();
    numberCounter = 0;
    closeMainActive = false;

    // Automatische Sicherung der Parameter beim Calc-Button.
    fnTempParamFile = fileDir.absoluteFilePath("TempParamSave.ini");

    restoreGeometry( sets.value("GUIgeometry").toByteArray() );

    // ---------------------------
    // TAB 0 --- Configuration ---
    // ---------------------------
    // Find the number of threads possible
    int nthreads = static_cast<int>(std::thread::hardware_concurrency());
    ui->inpNumCores->setMaximum( nthreads );
    // Minimum settings later
    ui->radNewImageCfg->setChecked(true);
#ifdef Q_OS_WIN
    bIgnoreUpdates = false;
#else
    ui->togIgnoreUpdates->setChecked(true);
    bIgnoreUpdates = true;
#endif
    ui->togAutoPosit->setChecked( sets.value("ImgAutoPosit",false).toBool() );
    imgPosX  = 0;
    imgPosX0 = 0;
    imgPosY  = 0;
    imgPosY0 = 0;

    ui->togOnlyNewWindow->setChecked( sets.value("OnlyNewWindow",false).toBool() );
    ui->togLimitRuntime->setChecked( sets.value("LimitRuntimeFlag",false).toBool() );
    ui->inpLimitRuntime->setValue( sets.value("LimitRuntimeValue",60).toInt() );
    on_togLimitRuntime_toggled( ui->togLimitRuntime->isChecked() );

    QStringList sl = widImage::slColorNames();
    sl.removeOne(tblResGlo);
    sl.removeOne(tblResTemp);
    ui->cbsDefaultColTbl->addItems( sl );
    ui->cbsDefaultColTbl->setCurrentIndex( sets.value("DefColTbl",1).toInt() );
    on_cbsDefaultColTbl_activated( sets.value("DefColTbl",1).toInt() );

    ui->inpParamSearchPath->setText( sets.value("ParamSearchPath","").toString() );
    ui->lisParamSearchResult->setEnabled(false);
    ui->butParamSearchGenerate->setEnabled(false);
    bSearchParamsRunning = false;

    // ---------------------------
    // TAB 1 --- Data ------------
    // ---------------------------
    ui->butDataSetMask->setEnabled(false);
    ui->butDataCopyScaling->setEnabled(false);
    ui->grpExtractImage->setEnabled(false);
    ui->grpNoFitRegions->setEnabled(false);
    ui->butDataFindCenter->setEnabled(false);
    ui->butSaveAllImages->hide(); // setEnabled(false); TODO

    ui->cbsExtractScale->clear();
    ui->cbsExtractScale->addItem(" * 1 ", 1 );
    ui->cbsExtractScale->addItem(" * 2 ", 2 );
    ui->cbsExtractScale->addItem(" * 4 ", 4 );

    // ---------------------------
    // TAB 2 --- Calculations ----
    // ---------------------------
    ui->cbsMethod->blockSignals(true);  // AI Tab, hier schon geblockt, da tabMethods dort auch anpasst
    ui->tabMethods->blockSignals(true);
    ui->radNewImageCal->setChecked(true);
    ui->tabMethods->clear();
    calcGui = new SC_CalcGUI;
    calcGui->createTabs( ui->tabMethods, ui->statusbar );   // TabBar versteckt, wenn nur eine Methode
    disableUnusedElements(); // TODO: Per Flag schaltbar?
    if ( lastTab >= 0 && ui->tabMethods->count() > 1 )
    {   // Damit auf jeden Fall das Signal ausgelöst wird, erst auf ein anderes Tab setzen
        ui->tabMethods->setCurrentIndex( lastTab==0 ? 1 : 0 );
        ui->tabMethods->blockSignals(false);
        ui->tabMethods->setCurrentIndex(lastTab);
    }
    else
    {
        ui->tabMethods->blockSignals(false);
        on_tabMethods_currentChanged(0);
    }
    if ( ui->tabMethods->count() == 1 )
        ui->radTestAllMethods->hide();
    // Check number of threads
    if ( calcGui->gpuAvailable() )
    {
        ui->inpNumCores->setValue(0); // 0=GPU
    }
    else
    {
        ui->inpNumCores->setSpecialValueText("");
        ui->inpNumCores->setMinimum(1); // keine GPU möglich
        ui->inpNumCores->setValue( nthreads ); // immer auf Maximum als Default
    }
    ui->inpGridPoints->setValue( sets.value("GridPoints",64).toInt() );  // Damit BeamCenter gesetzt wird
    ui->butCalc->setEnabled(true);
    ui->butAbbruch->setEnabled(false);
    ui->progressBar->setEnabled(false);
    ui->radQ1->setChecked( sets.value("Quadrant",2).toInt() == 1 );
    ui->radQ2->setChecked( sets.value("Quadrant",2).toInt() == 2 );
    ui->radQ4->setChecked( sets.value("Quadrant",2).toInt() == 4 );
    ui->togExpandImage->setChecked( sets.value("ExpandImg",true).toBool() );
    ui->togExpandImage->setEnabled( sets.value("Quadrant",2).toInt() != 4 );
    ui->lblLoadPrompt->hide();
    ui->lblLoadFilename->hide();

    // Jetzt noch die Parameter für die globalen Werte setzen....
    foreach ( QString m, calcGui->getCalcTypes() )
    {
        foreach ( QString p, calcGui->paramsForMethod(m,true,true,false) ) // numerical, global, no fit
        {
            paramHelper *phlp = calcGui->getParamPtr( m, p );
            if ( phlp != nullptr )
            {
                if ( phlp->gui.w == nullptr )
                {
                    if ( phlp->key == "Editdom1" ) phlp->gui.num = ui->inpSigX;
                    if ( phlp->key == "Editdom2" ) phlp->gui.num = ui->inpSigY;
                    if ( phlp->key == "Editdom3" ) phlp->gui.num = ui->inpSigZ;
                }
                if ( phlp->lbl1 == nullptr ) phlp->lbl1 = ui->lblParamSig;
            }
        } // p
    } // m

    // Slot connections for AutomaticRecalulation
    // In der UI sind leider nur die Slots für die vordefinierten Signale zu setzen.
    connect( ui->inpGridPoints, SIGNAL(valueChanged(int)), this, SLOT(automaticRecalc()) );
    connect( ui->inpHKLmax, SIGNAL(valueChanged(int)), this, SLOT(automaticRecalc()) );
    connect( ui->inpBCenterX, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpBCenterY, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpU1, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpU2, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpU3, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpV1, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpV2, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpV3, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpN1, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpN2, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpN3, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAx1, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAy1, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAz1, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAx2, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAy2, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAz2, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAx3, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAy3, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpAz3, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpSigX, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpSigY, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );
    connect( ui->inpSigZ, SIGNAL(valueChanged(double)), this, SLOT(automaticRecalc()) );

    // ---------------------------
    // TAB 3 --- FFT -------------
    // ---------------------------
    // Grp Input Preprocessing
    ui->radFFTLinInput->setChecked(  sets.value("FFTLinInput",true).toBool() );
    ui->radFFTLogInput->setChecked( !sets.value("FFTLinInput",true).toBool() );
    // Grp (r,phi)
    ui->grpFFTuseRphi->setChecked( sets.value("FFTuseRphi",true).toBool() );
    ui->cbsFFTsizeRphi->setCurrentIndex( sets.value("FFTsizeRphi",2).toInt() );
    ui->togFFTscaleRphi->setChecked( sets.value("FFTScaleRphi",true).toBool() );
    ui->togFFTclipRphi->setChecked( sets.value("FFTclipRphi",false).toBool() );
    ui->togFFTclip40Rphi->setChecked( sets.value("FFTclip40Rphi",false).toBool() );
    ui->togFFTdispRphi->setChecked( sets.value("DispRphi",true).toBool() );
    // Grp iFFT Output
    ui->togFFTscaleOut->setChecked( sets.value("FFTScaleOutput",true).toBool() );
    ui->togFFTclipOut->setChecked( sets.value("FFTclipOutput",false).toBool() );
    ui->togFFTclip40Out->setChecked( sets.value("FFTclip40Output",false).toBool() );
    ui->togIFFTSwap->setChecked( sets.value("FFTSwapOutput",true).toBool() );
    ui->cbsFFTsizeOut->setCurrentIndex( sets.value("FFTsizeOut",2).toInt() );
    switch ( sets.value("FFToutput",0).toInt() )
    {
    case 0: // Real
        ui->radFFToutReal->setChecked(true);
        break;
    case 1: // Imag
        ui->radFFToutImag->setChecked(true);
        break;
    case 2: // Abs
        ui->radFFToutBetrag->setChecked(true);
        break;
    case 3: // Spec
        ui->radFFToutSpectrum->setChecked(true);
        break;
    }
    ui->butIFFT->setEnabled(false); // Init
    ui->butFFTverifyRphi->hide();   // war nur für mich zum Testen, daher im Normalbetrieb wegnehmen

    // ---------------------------
    // TAB 4 --- Fit -------------
    // ---------------------------
#ifndef USEREPETITIONS
    ui->inpFitRepetitions->hide();
    ui->lblFitCurRep->hide();
    ui->lblFitRepetitions->hide();
    ui->lblPrmFitCurRep->hide();
#endif
    ui->butFitStart->setEnabled(false);
    ui->butFitUseResult->setEnabled(false);
    ui->butShowResiduen->setEnabled(false);
    ui->butFitHistoShow->setEnabled(false);
#ifndef UnusedValue
    ui->butFitHistoShow->hide();
    ui->butFitHistoClear->hide();
    ui->lblFitHisto->hide();
#endif
    ui->lblFitCurRep->setText("");
    ui->lblFitCurIter->setText("");
    ui->lblAutofitInfo->setText("");

    if ( configParamsFile.size() > 55 )
        ui->lblFitConfigFile->setText(configParamsFile.left(3) + "..." + configParamsFile.right(49));
    else
        ui->lblFitConfigFile->setText(configParamsFile);
    fitClass = nullptr;
    fileFitLog = nullptr;
    fileFitLogName = "";
    sets.beginGroup( "Fit-Globals" );
    ui->inpFitBStop->setValue( sets.value("BeamStop", 0).toInt() );
    ui->inpFitBorder->setValue( sets.value( "Border", 0).toInt() );
    ui->inpFitMaxIter->setValue( sets.value( "MaxIter", 10).toInt() );
    ui->inpFitStepSize->setValue( sets.value( "StepSize", 1.0).toDouble() );
    ui->inpFitTolerance->setText( QString::number( sets.value( "Tolerance", 0.5).toDouble() ) );
    ui->inpFitRepetitions->setValue( sets.value( "Repetitions", 1).toInt() );
    ui->togFitUseMask->setChecked( sets.value( "UseMask", false).toBool() );
    sets.endGroup();
    QStringList slAllMethods = calcGui->getCalcTypes();
    oneFitParamUsed = false;
    savedFitParams.clear();
    ui->butFitCurLoad->setEnabled(false);
    foreach (QString m, slAllMethods)
    {
        method2fitparams.insert( m, new _param2fitval );
        sets.beginGroup( "Fit-"+m );
        QStringList slKeys = sets.allKeys();
        //qDebug() << m << slKeys;
        foreach (QString k, slKeys)
        {
            if ( ! calcGui->isCurrentParameterValid(m,k,true) ) continue;   // Hier nur Fitparameter
            QStringList slVal = sets.value(k,"0:0:0:0").toString().split(":");
            _fitLimits *fl = new _fitLimits;
            fl->used = slVal[0].toInt() != 0;
            fl->min  = slVal[1].toDouble();
            fl->max  = slVal[2].toDouble();
            fl->fitType = static_cast<_fitTypes>(slVal[3].toInt());
            fl->fitvalid = false;  // fitresult not stored!
            //qDebug() << "Init:" << m << k << slVal;
            method2fitparams[m]->insert(k,fl);
            oneFitParamUsed |= fl->used;
        }
        sets.endGroup();
    }

    // ---------------------------
    // TAB 5 --- AI --------------
    // ---------------------------
    // Die globalen Daten richtig übernehmen
    fillDataFromInputs();

    sets.beginGroup("AI");
    ui->inpSubDir->setText( sets.value("LastSubDir",".").toString() );
    ui->cbsAIoutputFormat->setCurrentIndex( sets.value("Grayscale",1).toInt() );
    ui->grpFileInput->setChecked( sets.value("FileInputEna",false).toBool() );
    ui->inpFileName->setText( sets.value("FileInputLast","").toString() );
    ui->inpFileClass->setText( sets.value("FileClass","{M}").toString() );
    ui->togAIUseFFTOutput->setChecked( sets.value("GenerateIFFT",false).toBool() );
    ui->radAILinOutput->setChecked(   sets.value("LinOut",true).toBool() );
    ui->radAILogOutput->setChecked( ! sets.value("LinOut",true).toBool() );
    ui->togAIScaleOutput->setChecked( sets.value("ScaleOut",false).toBool() );
    // Diese Schleifen werden immer fest im Cascadiermodus gestartet.
    ui->radLoopsCascade->setChecked( true /*   sets.value("Cascading",true).toBool()*/ );
    ui->radLoopsSingle->setChecked( false /* ! sets.value("Cascading",true).toBool()*/ );
    sets.endGroup();
    QString tabName = ui->tabMethods->tabText(lastTab);
    ui->cbsMethod->addItems( calcGui->getCalcTypes() );
    ui->cbsMethod->setCurrentIndex( -1 );
    ui->cbsMethod->blockSignals(false);
    ui->cbsMethod->setCurrentIndex( ui->cbsMethod->findText(tabName) );
    QStringList hdr;
    hdr << "Name" << "Start" << "End / Variation" << "Step / Calculation";
    ui->tblListe->setColumnCount( hdr.size() );
    ui->tblListe->setHorizontalHeaderLabels( hdr );
    ui->tblListe->resizeColumnsToContents();

    ui->butAIcheck->setEnabled(false);        // Erst freigegeben, wenn die AI-Tabelle gefüllt ist
    ui->butAIsaveBkgOp->setEnabled(false);    // -"-
    ui->butAIstart->setEnabled(false);        // -"-
    ui->butSaveImagesAI->setEnabled(false); // Nur Freigegeben, wenn images da sind
    ui->butRemoveVar->setEnabled(false);

    ui->lblAIbackProc->setEnabled(false);
    ui->lisAIbackProgOut->setEnabled(false);
    ui->lblTPVbackProc->setEnabled(false);
    ui->lisTPVbackProgOut->setEnabled(false);
    aiBackProg = nullptr;

    if ( calcGui->getCalcTypes().size() == 1 )
    {   // Bei nur einer Methode die Auswahlen bei AI sperren
        ui->cbsMethod->hide(); // setEnabled(false);
        ui->togUseAllMethods->hide(); //setEnabled(false);
        ui->lblAIMethods->hide();
        ui->hspacerMethods->changeSize(1,1);
    }

    // --------------------------------------------
    // TAB 6 --- Training Parameters Variation ----
    // --------------------------------------------
    sets.beginGroup("TPV");
    ui->togTPVio->setChecked( sets.value("ena_io",true).toBool() );
    ui->inpTPVio->setValue( sets.value("val_io",0.5).toDouble() );
    ui->inpTPVio->setEnabled(ui->togTPVio->isChecked());
    ui->togTPVbase->setChecked( sets.value("ena_base",true).toBool() );
    ui->inpTPVbase->setValue( sets.value("val_base",0.5).toDouble() );
    ui->inpTPVbase->setEnabled(ui->togTPVbase->isChecked());
    ui->togTPVradius->setChecked( sets.value("ena_radius",true).toBool() );
    ui->inpTPVradius->setValue( sets.value("val_radius",0.5).toDouble() );
    ui->inpTPVradius->setEnabled(ui->togTPVradius->isChecked());
    ui->togTPVsigmar->setChecked( sets.value("ena_sigmar",true).toBool() );
    ui->inpTPVsigmar->setValue( sets.value("val_sigmar",0.5).toDouble() );
    ui->inpTPVsigmar->setEnabled(ui->togTPVsigmar->isChecked());
    ui->togTPVlength->setChecked( sets.value("ena_length",true).toBool() );
    ui->inpTPVlength->setValue( sets.value("val_length",0.5).toDouble() );
    ui->inpTPVlength->setEnabled(ui->togTPVlength->isChecked());
    ui->togTPVsigmal->setChecked( sets.value("ena_sigmal",true).toBool() );
    ui->inpTPVsigmal->setValue( sets.value("val_sigmal",0.5).toDouble() );
    ui->inpTPVsigmal->setEnabled(ui->togTPVsigmal->isChecked());
    ui->togTPVphi->setChecked( sets.value("ena_phi",true).toBool() );
    //ui->inpTPVphi->setValue( sets.value("val_phi",0.5).toDouble() );
    //ui->inpTPVphi->setEnabled(ui->togTPVphi->isChecked());
    ui->togTPVa->setChecked( sets.value("ena_a",true).toBool() );
    ui->inpTPVa->setValue( sets.value("val_a",0.5).toDouble() );
    ui->inpTPVa->setEnabled(ui->togTPVa->isChecked());
    ui->togTPVb->setChecked( sets.value("ena_b",true).toBool() );
    ui->inpTPVb->setValue( sets.value("val_b",0.5).toDouble() );
    ui->inpTPVb->setEnabled(ui->togTPVb->isChecked());
    ui->togTPVc->setChecked( sets.value("ena_c",true).toBool() );
    ui->inpTPVc->setValue( sets.value("val_c",0.5).toDouble() );
    ui->inpTPVc->setEnabled(ui->togTPVc->isChecked());    
    ui->togTPVrho->setChecked( sets.value("ena_rho",true).toBool() );
    ui->inpTPVrho->setValue( sets.value("val_rho",0.5).toDouble() );
    ui->inpTPVrho->setEnabled(ui->togTPVrho->isChecked());
    ui->togTPVpsi->setChecked( sets.value("ena_psi",true).toBool() );
    ui->inpTPVpsi->setValue( sets.value("val_psi",0.5).toDouble() );
    ui->inpTPVpsi->setEnabled(ui->togTPVpsi->isChecked());
    ui->togTPVdbeta->setChecked( sets.value("ena_dbeta",true).toBool() );
    ui->inpTPVdbeta->setValue( sets.value("val_dbeta",0.5).toDouble() );
    ui->inpTPVdbeta->setEnabled(ui->togTPVdbeta->isChecked());
    ui->togTPVwidth->setChecked( sets.value("ena_width",true).toBool() );
    ui->inpTPVwidth->setValue( sets.value("val_width",0.5).toDouble() );
    ui->inpTPVwidth->setEnabled(ui->togTPVwidth->isChecked());
    ui->togTPVphiwidth->setChecked( sets.value("ena_phiwidth",true).toBool() );
    ui->inpTPVphiwidth->setValue( sets.value("val_phiwidth",0.5).toDouble() );
    ui->inpTPVphiwidth->setEnabled(ui->togTPVphiwidth->isChecked());
    ui->togTPVdisplacement->setChecked( sets.value("ena_displacement",true).toBool() );
    ui->inpTPVdisplacement->setValue( sets.value("val_displacement",0.5).toDouble() );
    ui->inpTPVdisplacement->setEnabled(ui->togTPVdisplacement->isChecked());
    ui->togTPVaddBS->setChecked( sets.value("ena_addBS",true).toBool() );
    ui->inpTPVbeamx0->setValue( sets.value("val_beamx0",0.5).toDouble() );
    ui->inpTPVbeamx0->setEnabled(ui->togTPVaddBS->isChecked());
    ui->inpTPVbeamy0->setValue( sets.value("val_beamy0",0.5).toDouble() );
    ui->inpTPVbeamy0->setEnabled(ui->togTPVaddBS->isChecked());
    ui->inpTPVnbeam->setValue( sets.value("val_nbeam",0.5).toDouble() );
    ui->inpTPVnbeam->setEnabled(ui->togTPVaddBS->isChecked());
    ui->inpTPVmbeam->setValue( sets.value("val_mbeam",0.5).toDouble() );
    ui->inpTPVmbeam->setEnabled(ui->togTPVaddBS->isChecked());
    // TODO: addLines
    ui->togTPVaddLines->setToolTip("Not yet implemented");
    ui->togTPVaddLines->setChecked( false /*sets.value("ena_addLines",false).toBool()*/ );
    ui->togTPVaddLines->setEnabled(false);
    ui->inpTPVaddLinesH->setValue( sets.value("val_addLinesH",1).toInt() );
    ui->inpTPVaddLinesV->setValue( sets.value("val_addLinesV",5).toInt() );
    ui->inpTPVaddLinesH->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesV->setEnabled(ui->togTPVaddLines->isChecked());
    //
    ui->togTPVaddNoise->setChecked( sets.value("ena_addNoise",true).toBool() );
    // TODO: Convolute
    ui->togTPVconvolute->setToolTip("Not yet implemented");
    ui->togTPVconvolute->setChecked( false /*sets.value("ena_convolute",true).toBool()*/ );
    ui->togTPVconvolute->setEnabled(false);
    //
    ui->togTPVcalcRphi->setChecked( sets.value("ena_calcRphi",true).toBool() );
    ui->togTPVcalcFFT->setChecked( sets.value("ena_calcFFT",true).toBool() );
    ui->inpTPVnumimg->setValue( sets.value("val_numimg",10).toInt() );
    ui->inpTPVoutPath->setText( sets.value("val_outPath",".").toString() );
    sets.endGroup();

    // ---------------------------
    // TAB 7 --- Imagewindows ----
    // ---------------------------
#ifdef IMG_ONE_WINDOW
    qDebug() << "********************* ANDROID *******************************";


#else
    ui->tabMain->removeTab(7);
#endif

    // -------------------------
    // Abschluss ---------------
    // -------------------------
    ui->tabMain->setCurrentWidget( ui->tabCalc );
    adjustSize();

    if ( qApp->arguments().size() > 1 )     // [0] = Executable file path
    {
        autoProcessingFile = new QFile( qApp->arguments().at(1) );
        if ( autoProcessingFile->open(QIODevice::ReadOnly) )
        {
            qDebug() << autoProcessingFile->fileName() << "Starting";
            QTimer::singleShot( 100, this, SLOT(autoProcessingTimer()) );
        }
        else
        {
            qDebug() << autoProcessingFile->fileName() << autoProcessingFile->errorString();
            autoProcessingFile = nullptr;
            bIgnoreRecalc = false;  // Init done AutoProc-Error
        }
    }
    else
    {
        bIgnoreRecalc = false;  // Init done w/o AutoProc
    }
}

SC_MainGUI::~SC_MainGUI()
{
    delete ui;
}

void SC_MainGUI::closeEvent(QCloseEvent *event)
{
    if ( !closeMainActive )
    {
        closeMainActive = true; // Avoid calls to imageWindowClosed()

        QSettings sets(SETT_APP,SETT_GUI);
        // Calculation & Configuration Tab (Default group)
        sets.setValue( "LastMethod", ui->tabMethods->currentIndex() );
        sets.setValue( "GridPoints", ui->inpGridPoints->value() );
        sets.setValue( "ImgAutoPosit", ui->togAutoPosit->isChecked() );
        sets.setValue( "DefColTbl", ui->cbsDefaultColTbl->currentIndex() );
        sets.setValue( "OnlyNewWindow", ui->togOnlyNewWindow->isChecked() );
        sets.setValue( "LimitRuntimeFlag", ui->togLimitRuntime->isChecked() );
        sets.setValue( "LimitRuntimeValue", ui->inpLimitRuntime->value() );
        sets.setValue( "GUIgeometry", saveGeometry() );
        if ( ui->radQ1->isChecked() ) sets.setValue("Quadrant",1);
        if ( ui->radQ2->isChecked() ) sets.setValue("Quadrant",2);
        if ( ui->radQ4->isChecked() ) sets.setValue("Quadrant",4);
        sets.setValue( "ExpandImg", ui->togExpandImage->isChecked() );
        sets.setValue( "ParamSearchPath", ui->inpParamSearchPath->text() );

        // FIT Tab (Subgroups)
        QHash< QString/*method*/, _param2fitval* >::const_iterator im = method2fitparams.constBegin();
        while ( im != method2fitparams.constEnd() )
        {
            if ( im.value() != nullptr )
            {
                sets.beginGroup( "Fit-" + im.key() );
                QHash< QString/*name*/, _fitLimits* >::const_iterator il = im.value()->constBegin();
                while ( il != im.value()->constEnd() )
                {
                    sets.setValue( il.key(), QString("%1:%2:%3:%4").arg(il.value()->used).arg(il.value()->min).arg(il.value()->max).arg(il.value()->fitType) );
                    ++il;
                }
                sets.endGroup();
            }
            ++im;
        }
        sets.beginGroup( "Fit-Globals" );
        sets.setValue( "BeamStop",    ui->inpFitBStop->value() );
        sets.setValue( "Border",      ui->inpFitBorder->value() );
        sets.setValue( "MaxIter",     ui->inpFitMaxIter->value() );
        sets.setValue( "StepSize",    ui->inpFitStepSize->value() );
        sets.setValue( "Tolerance",   ui->inpFitTolerance->text().toDouble() );
        sets.setValue( "Repetitions", ui->inpFitRepetitions->value() );
        sets.setValue( "UseMask",     ui->togFitUseMask->isChecked() );
        sets.endGroup();

        // FFT Tab (Default group)
        sets.setValue( "FFTLinInput", ui->radFFTLinInput->isChecked() );
        sets.setValue( "FFTScaleRphi", ui->togFFTscaleRphi->isChecked() );
        sets.setValue( "FFTclipRphi", ui->togFFTclipRphi->isChecked() );
        sets.setValue( "FFTclip40Rphi", ui->togFFTclip40Rphi->isChecked() );
        sets.setValue( "FFTsizeRphi", ui->cbsFFTsizeRphi->currentIndex() );
        sets.setValue( "DispRphi", ui->togFFTdispRphi->isChecked() );
        sets.setValue( "FFTuseRphi", ui->grpFFTuseRphi->isChecked() );
        sets.setValue( "FFTScaleOutput", ui->togFFTscaleOut->isChecked() );
        sets.setValue( "FFTclipOutput", ui->togFFTclipOut->isChecked() );
        sets.setValue( "FFTclip40Output", ui->togFFTclip40Out->isChecked() );
        sets.setValue( "FFTSwapOutput", ui->togIFFTSwap->isChecked() );
        sets.setValue( "FFTsizeOut", ui->cbsFFTsizeOut->currentIndex() );
        if ( ui->radFFToutReal->isChecked() )
            sets.setValue( "FFToutput", 0 );
        else if ( ui->radFFToutImag->isChecked() )
            sets.setValue( "FFToutput", 1 );
        else if ( ui->radFFToutBetrag->isChecked() )
            sets.setValue( "FFToutput", 2 );
        else if ( ui->radFFToutSpectrum->isChecked() )
            sets.setValue( "FFToutput", 3 );

        // AI definitions Tab (AI Group)
        sets.beginGroup("AI");
        sets.setValue( "LastSubDir", ui->inpSubDir->text() );
        sets.setValue( "Grayscale", ui->cbsAIoutputFormat->currentIndex() );
        sets.setValue( "FileInputEna", ui->grpFileInput->isChecked() );
        sets.setValue( "FileInputLast", ui->inpFileName->text() );
        sets.setValue( "FileClass", ui->inpFileClass->text() );
        sets.setValue( "GenerateIFFT", ui->togAIUseFFTOutput->isChecked() );
        sets.setValue( "LinOut", ui->radAILinOutput->isChecked() );
        sets.setValue( "ScaleOut", ui->togAIScaleOutput->isChecked() );
        //sets.setValue( "Cascading", ui->radLoopsCascade->isChecked() );  s.o.
        sets.endGroup();

        // TAB 6 --- Training Parameters Variation ----
        sets.beginGroup("TPV");
        sets.setValue("ena_io",ui->togTPVio->isChecked());
        sets.setValue("val_io",ui->inpTPVio->value());        
        sets.setValue("ena_base",ui->togTPVbase->isChecked());
        sets.setValue("val_base",ui->inpTPVbase->value());
        sets.setValue("ena_radius",ui->togTPVradius->isChecked());
        sets.setValue("val_radius",ui->inpTPVradius->value());
        sets.setValue("ena_sigmar",ui->togTPVsigmar->isChecked());
        sets.setValue("val_sigmar",ui->inpTPVsigmar->value());
        sets.setValue("ena_length",ui->togTPVlength->isChecked());
        sets.setValue("val_length",ui->inpTPVlength->value());
        sets.setValue("ena_sigmal",ui->togTPVsigmal->isChecked());
        sets.setValue("val_sigmal",ui->inpTPVsigmal->value());
        sets.setValue("ena_phi",ui->togTPVphi->isChecked());
        //sets.setValue("val_phi",ui->inpTPVphi->value());
        sets.setValue("ena_a",ui->togTPVa->isChecked());
        sets.setValue("val_a",ui->inpTPVa->value());
        sets.setValue("ena_b",ui->togTPVb->isChecked());
        sets.setValue("val_b",ui->inpTPVb->value());
        sets.setValue("ena_c",ui->togTPVc->isChecked());
        sets.setValue("val_c",ui->inpTPVc->value());        
        sets.setValue("ena_rho",ui->togTPVrho->isChecked());
        sets.setValue("val_rho",ui->inpTPVrho->value());
        sets.setValue("ena_psi",ui->togTPVpsi->isChecked());
        sets.setValue("val_psi",ui->inpTPVpsi->value());        
        sets.setValue("ena_dbeta",ui->togTPVdbeta->isChecked());
        sets.setValue("val_dbeta",ui->inpTPVdbeta->value());
        sets.setValue("ena_width",ui->togTPVwidth->isChecked());
        sets.setValue("val_width",ui->inpTPVwidth->value());
        sets.setValue("ena_phiwidth",ui->togTPVphiwidth->isChecked());
        sets.setValue("val_phiwidth",ui->inpTPVphiwidth->value());
        sets.setValue("ena_displacement",ui->togTPVdisplacement->isChecked());
        sets.setValue("val_displacement",ui->inpTPVdisplacement->value());
        sets.setValue("val_beamx0",ui->inpTPVbeamx0->value());
        sets.setValue("val_beamy0",ui->inpTPVbeamy0->value());
        sets.setValue("val_nbeam",ui->inpTPVnbeam->value());
        sets.setValue("val_mbeam",ui->inpTPVmbeam->value());
        sets.setValue("ena_addLines",ui->togTPVaddLines->isChecked());
        sets.setValue("val_addLinesH",ui->inpTPVaddLinesH->value());
        sets.setValue("val_addLinesV",ui->inpTPVaddLinesV->value());
        sets.setValue("ena_addBS",ui->togTPVaddBS->isChecked());
        sets.setValue("ena_addNoise",ui->togTPVaddNoise->isChecked());
        sets.setValue("ena_convolute",ui->togTPVconvolute->isChecked());
        sets.setValue("ena_calcRphi",ui->togTPVcalcRphi->isChecked());
        sets.setValue("ena_calcFFT",ui->togTPVcalcFFT->isChecked());
        sets.setValue("val_numimg",ui->inpTPVnumimg->value());
        sets.setValue("val_outPath",ui->inpTPVoutPath->text());
        sets.endGroup();

        // close all image windows
        for ( int i=0; i<images.size(); i++ )
            images.at(i)->close();

        // free memory in the GPU
        calcGui->cleanup();
    }
    if ( event )
        event->accept(); // close the main window
}

void SC_MainGUI::on_actionExit_triggered()
{
    closeEvent(nullptr);
    qApp->quit();
}

void SC_MainGUI::on_tabMethods_currentChanged(int index)
{
    if ( calcGui->getCalcTypes().size() == 1 )
    {
        curMethod = ui->tabMethods->tabText(index);
        return; // komplett ignorieren, wenn keine Tabs zur Auswahl stehen
    }
    // Aktuelle Methode speichern
    //QSettings data(SETT_APP,SETT_GUI);
    //data.setValue("LastMethod",index);
    // Den Button "Calculate" mit dem passenden Methodentext ergänzen, damit diese immer sichtbar ist
    QString hdr = ui->tabMethods->tabText(index);
    curMethod = hdr; // calcGui->index2methodname(ui->tabMethods->currentIndex());
    int pos = hdr.indexOf(" ");
    if ( pos < 0 ) pos = hdr.length();
    ui->butCalc->setText( "Calculate "+hdr.left(pos) );
    // Auf neues Image umschalten
    ui->radNewImageCal->setChecked(true);
    // Die aktuelle Methode für AI anpassen
    ui->cbsMethod->setCurrentIndex( ui->cbsMethod->findText(hdr) );
}

void SC_MainGUI::fillDataFromInputs()
{
    // Globale Inputs für die Berechnungen
    //QHash<QString,Double3> SC_CalcGUI::inpVectors;
    //QHash<QString,double>  SC_CalcGUI::inpValues;
    //QHash<QString,double>  SC_CalcGUI::inpSingleValueVectors;

    SC_CalcGUI::inpValues.insert( "EditGridPoints", ui->inpGridPoints->value() );
    SC_CalcGUI::inpValues.insert( "Edithklmax", ui->inpHKLmax->value() );
    SC_CalcGUI::inpValues.insert( "RadioButtonQ1", ui->radQ1->isChecked() );
    SC_CalcGUI::inpValues.insert( "RadioButtonQ2", ui->radQ2->isChecked() );
    SC_CalcGUI::inpValues.insert( "RadioButtonQ4", ui->radQ4->isChecked() );
    SC_CalcGUI::inpValues.insert( "ExpandImage", ui->togExpandImage->isChecked() );
    SC_CalcGUI::inpValues.insert( "BeamPosX", ui->inpBCenterX->value() );
    SC_CalcGUI::inpValues.insert( "BeamPosY", ui->inpBCenterY->value() );

    // Lattice Werte
    SC_CalcGUI::inpValues.insert( "LATTcols", ui->tblLattice3DValues->item(tblLattNbCols,0)->text().toInt() );
    SC_CalcGUI::inpValues.insert( "LATTrows", ui->tblLattice3DValues->item(tblLattNbRows,0)->text().toInt() );
    SC_CalcGUI::inpValues.insert( "LATTcenx", ui->tblLattice3DValues->item(tblLattBeamX,0)->text().toDouble() );
    SC_CalcGUI::inpValues.insert( "LATTceny", ui->tblLattice3DValues->item(tblLattBeamY,0)->text().toDouble() );
    SC_CalcGUI::inpValues.insert( "LATTwlen", ui->tblLattice3DValues->item(tblLattWaveLen,0)->text().toDouble() );
    SC_CalcGUI::inpValues.insert( "LATTdist", ui->tblLattice3DValues->item(tblLattSaDist,0)->text().toDouble() );
    SC_CalcGUI::inpValues.insert( "LATTpixx", ui->tblLattice3DValues->item(tblLattPixelX,0)->text().toDouble() );
    SC_CalcGUI::inpValues.insert( "LATTpixy", ui->tblLattice3DValues->item(tblLattPixelY,0)->text().toDouble() );

    SC_CalcGUI::inpVectors.insert( "Uvec", Double3(ui->inpU1->value(),
                                                   ui->inpU2->value(),
                                                   ui->inpU3->value()) );
    SC_CalcGUI::inpVectors.insert( "Vvec", Double3(ui->inpV1->value(),
                                                   ui->inpV2->value(),
                                                   ui->inpV3->value()) );
    SC_CalcGUI::inpVectors.insert( "Nvec", Double3(ui->inpN1->value(),
                                                   ui->inpN2->value(),
                                                   ui->inpN3->value()) );
    SC_CalcGUI::inpVectors.insert( "Ax1", Double3(ui->inpAx1->value(),
                                                  ui->inpAy1->value(),
                                                  ui->inpAz1->value()) );
    SC_CalcGUI::inpVectors.insert( "Ax2", Double3(ui->inpAx2->value(),
                                                  ui->inpAy2->value(),
                                                  ui->inpAz2->value()) );
    SC_CalcGUI::inpVectors.insert( "Ax3", Double3(ui->inpAx3->value(),
                                                  ui->inpAy3->value(),
                                                  ui->inpAz3->value()) );
    SC_CalcGUI::inpVectors.insert( "SigXYZ", Double3(ui->inpSigX->value(),
                                                     ui->inpSigY->value(),
                                                     ui->inpSigZ->value()) );

    // In normal calculation modes the vector notation is used, but in the
    // AI mode all the values are accessed as simgle values.
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis1x",  SC_CalcGUI::inpVectors["Ax1"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis1y",  SC_CalcGUI::inpVectors["Ax1"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis1z",  SC_CalcGUI::inpVectors["Ax1"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis2x",  SC_CalcGUI::inpVectors["Ax2"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis2y",  SC_CalcGUI::inpVectors["Ax2"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis2z",  SC_CalcGUI::inpVectors["Ax2"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis3x",  SC_CalcGUI::inpVectors["Ax3"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis3y",  SC_CalcGUI::inpVectors["Ax3"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "EditAxis3z",  SC_CalcGUI::inpVectors["Ax3"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editxrel",    SC_CalcGUI::inpVectors["Nvec"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Edityrel",    SC_CalcGUI::inpVectors["Nvec"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editzrel",    SC_CalcGUI::inpVectors["Nvec"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editx1rel",   SC_CalcGUI::inpVectors["Uvec"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Edity1rel",   SC_CalcGUI::inpVectors["Uvec"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editz1rel",   SC_CalcGUI::inpVectors["Uvec"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editx2rel",   SC_CalcGUI::inpVectors["Vvec"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Edity2rel",   SC_CalcGUI::inpVectors["Vvec"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editz2rel",   SC_CalcGUI::inpVectors["Vvec"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editdom1",    SC_CalcGUI::inpVectors["SigXYZ"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editdom2",    SC_CalcGUI::inpVectors["SigXYZ"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editdom3",    SC_CalcGUI::inpVectors["SigXYZ"].z() );
}

void SC_MainGUI::prepareCalculation( bool getMD, bool progbar )
{
    ui->butCalc->setEnabled(false);
    ui->butFitStart->setEnabled(false);
    ui->butIFFT->setEnabled(false); // Prepare Calc
    ui->butAbbruch->setEnabled(true);
    ui->butAbbruch->setText("Abort");
    ui->progressBar->setEnabled(progbar);
    qApp->processEvents();
    if ( getMD ) fillDataFromInputs();
    calcGui->prepareCalculation( ui->tabMethods->tabText(ui->tabMethods->currentIndex()), getMD );
}

void SC_MainGUI::finishCalculation(bool showtime)
{
    if ( showtime )
    {
        ui->lblTimeUsed->setText( QString("%1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimBoth)) );
        ui->statusbar->showMessage( QString("%1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimBoth)), 5000 );
    }
    ui->progressBar->setValue(100); // it is possible (due to roundings) that the 100% is not reached. So set it here to avoid blinking.
    ui->butCalc->setEnabled(true);
    ui->butAbbruch->setEnabled(false);
    ui->butAbbruch->setText("Abort");
    ui->progressBar->setEnabled(false);
    ui->butFitStart->setEnabled(oneFitParamUsed);
    ui->butIFFT->setEnabled(true);  // Finish Calc
    calcGui->inpValues.insert( "_CalcTime_", calcGui->higResTimerElapsed(SC_CalcGUI::htimCalc) );
    calcGui->inpValues.insert( "_PrepTime_", calcGui->higResTimerElapsed(SC_CalcGUI::htimPrep) );
}

void SC_MainGUI::on_butCalc_clicked()
{
#ifdef COPY_FITDATA_TO_GPU  // butCalc
    if ( sender() != nullptr )
        calcGui->setArrDataForFit( nullptr );
#endif
#ifdef FITDATA_IN_GPU  // butCalc
    if ( sender() != nullptr )
        calcGui->setFitData( 0, 0, nullptr );
#endif
    if ( _calcThread == nullptr )
        _calcThread = new myCalcThread(calcGui);
    _calcThread->setThreads( ui->inpNumCores->value() );

    // Da es mir schon öfter passiert ist, dass ich Parameter ausgesucht und dann das
    // Programm einfach geschlossen hatte, ohne die Parameter zu sichern. Daher werden
    // bei jedem Berechnungsstart alle Parameter in einer temporären Datei gespeichert.
    performSaveParamOperation( fnTempParamFile );

    prepareCalculation( true, true );

    if ( _logThreadTimer == nullptr )
    {
        _logThreadTimer = new QTimer;
        connect( _logThreadTimer, SIGNAL(timeout()), this, SLOT(logThreadTimer()) );
    }
    lastX = -1;
    lastY = -1;
    lastPrz = -1;
    ui->progressBar->setValue(0);
    if ( ui->radQ4->isChecked() )
    {
        zzmin = - ui->inpGridPoints->value();
        zzrange = 2 * ui->inpGridPoints->value();
    }
    else
    {
        zzmin = 0;
        zzrange = ui->inpGridPoints->value();
    }
    if ( ui->radQ1->isChecked() )
    {
        iimin = 0;
        iirange = ui->inpGridPoints->value();
    }
    else
    {
        iimin = - ui->inpGridPoints->value();
        iirange = 2 * ui->inpGridPoints->value();
    }
    lastTime = 0;

    _bAbbruch = false;
    if ( ui->togLimitRuntime->isChecked() )
        calcMaxTime = ui->inpLimitRuntime->value() * 1000;
    else
        calcMaxTime = 0;
    calcRunTime.start();
    _calcThread->start();
    _logThreadTimer->start( 100 );
    waitForCalcThread();
    //qDebug() << "MainGUI: nach Wait" << _bAbbruch;
    //_calcThread->beenden(); // noch notwendig?
    _logThreadTimer->stop();

    finishCalculation(true);
    // Generate Image window
    lastUsedImage = addImage( false,    // Die Radio-Buttons werden beachtet
                              calcGui->minX(), calcGui->maxX(),
                              calcGui->minY(), calcGui->maxY(),
                              calcGui->data(), "Image", true );

    if ( !bIgnoreRecalc && autoProcessingFile!=nullptr && autoProcessingFile->isOpen() )
    {
        QTimer::singleShot( 200, this, SLOT(autoProcessingTimer()) );
    }
    //else
    //    qDebug() << "butCalc:" << bIgnoreRecalc << ((autoProcessingFile!=nullptr)?autoProcessingFile->isOpen():false);
}

myCalcThread::myCalcThread(SC_CalcGUI *cg, QObject *parent ) : QThread(parent)
{
    _exit = false;
    calcGui = cg;
}

void myCalcThread::run()
{
    calcGui->doCalculation( numThreads, nullptr );
}


/**
 * @brief SC_MainGUI::automaticRecalc [SLOT]
 * Called if any of the parameter values are changed and calls the calculation if this
 * is enabled with the QCheckBox 'togAutoRecalc'.
 */
void SC_MainGUI::automaticRecalc()
{
    if ( bIgnoreRecalc ) return; // Werte wurden vom Programm geändert
    if ( ! ui->togAutoRecalc->isChecked() ) return; // Recalc vom User gesperrt

    on_butCalc_clicked();

    // Wenn die Rechnung länger als 1 sec dauert, dann diese Funktion abschalten
    if ( calcRunTime.elapsed() > 1000 )
        ui->togAutoRecalc->setChecked(false);
}



/**
 * @brief SC_MainGUI::on_actionTest_10_Calculate_triggered
 * Es werden 10 Calculation-Durchläufe gemacht, immer in das gleiche Fenster und mit den
 * gleichen Parametern. Dazu wird die Zeit gewertet (Min/Mean/Max) und zwar für Calc und
 * Prep und zum Schluß als Metadata zum Image gegeben.
 */
void SC_MainGUI::on_actionTest_10_Calculate_triggered()
{
    dlgTimeTests *dlg = new dlgTimeTests(dataPath,calcGui->gpuAvailable(),ui->inpNumCores->maximum(),this);
    if ( dlg->exec() != QDialog::Accepted ) return;

    int numLoops = dlg->getRepetitions();
    QVector<int> threads   = dlg->getThreads();
    QVector<int> hklmax    = dlg->getHKLmax();
    QVector<int> quadrants = dlg->getQuadrants();
    QString saveFilename   = dlg->getSaveFilename();
    QString saveBasePath   = "";
    if ( ! saveFilename.isEmpty() )
    {
        if ( saveFilename.indexOf("/") < 0 )
            saveFilename = dataPath+"/"+saveFilename;
        QFile f(saveFilename);
        if ( ! f.open(QIODevice::Append) )
            if ( ! f.open(QIODevice::WriteOnly) )
            {
                qDebug() << saveFilename << f.errorString();
                return;
            }
        f.write( qPrintable(EOL "% " + dlg->getComment() + EOL) );
        f.write( qPrintable("% Loaded parameter: "+ui->lblLoadFilename->text()+EOL) );
        f.write( " & Threads & HKLmax & Quadrants & Min/ Mean/ Max (Prep) in ms" EOL );
        f.close();
        QFileInfo fi(saveFilename);
        saveBasePath = fi.absolutePath()+"/"+fi.baseName(); // save images
    }

    double c, sumCalc, minCalc, maxCalc, p, sumPrep, minPrep, maxPrep;
    QStringList latex;
    bool savExpand = ui->togExpandImage->isChecked();
    ui->togExpandImage->setChecked(true);
    bool savQ1 = ui->radQ1->isChecked();
    bool savQ2 = ui->radQ2->isChecked();
    bool savQ4 = ui->radQ4->isChecked();
    bool savUpd = ui->togIgnoreUpdates->isChecked();
    if ( dlg->getEnaUpdates() )
    {   // GUI-Updates sind gewünscht
        ui->togIgnoreUpdates->setChecked( false );
        bIgnoreUpdates = false;
    }
    else
    {   // Keine Updates bei den Berechnungen
        ui->togIgnoreUpdates->setChecked( true );
        bIgnoreUpdates = true;
    }
    bool saveImages = dlg->getEnaSaveImages();

    foreach ( int quadr, quadrants )
    {
        if ( quadr > 0 )
        {
            ui->radQ1->setChecked( quadr == 1 );
            ui->radQ2->setChecked( quadr == 2 );
            ui->radQ4->setChecked( quadr == 4 );
        }
        for ( int h=0; h<hklmax.size(); h++ )
        {
            ui->inpHKLmax->setValue( hklmax[h] );
            for ( int t=0; t<threads.size(); t++ )
            {
                if ( threads[t] == 0 && !calcGui->gpuAvailable() ) continue;
                if ( threads[t] == 100 )
                {
                    if ( ui->inpNumCores->maximum() == 4 ) continue;
                    ui->inpNumCores->setValue( ui->inpNumCores->maximum() );
                }
                else
                {
                    ui->inpNumCores->setValue( threads[t] ); // 0=GPU
                }

                sumCalc=0;
                minCalc=10000000;
                maxCalc=0;
                sumPrep=0;
                minPrep=10000000;
                maxPrep=0;
                for ( int i=0; i<numLoops; i++ )
                {
                    statusBar()->showMessage( QString("TimingTest: Threads=%1, HKLmax=%2, Quadrants=%5, Loop %3/%4")
                                                 .arg(ui->inpNumCores->value()).arg(hklmax[h])
                                                 .arg(i+1).arg(numLoops).arg(quadr) );
                    on_butCalc_clicked();
                    if ( _bAbbruch ) break;
                    c = calcGui->higResTimerElapsed(SC_CalcGUI::htimCalc);
                    p = calcGui->higResTimerElapsed(SC_CalcGUI::htimPrep);
                    sumCalc += c;
                    sumPrep += p;
                    if ( p < minPrep ) minPrep = p;
                    else if ( p > maxPrep ) maxPrep = p;
                    if ( c < minCalc ) minCalc = c;
                    else if ( c > maxCalc ) maxCalc = c;
                }
                if ( _bAbbruch ) break;
                lastUsedImage->addMetaInfo( "_NCalcMin_", QString("%1 ms").arg(minCalc) );
                lastUsedImage->addMetaInfo( "_NCalcMax_", QString("%1 ms").arg(maxCalc) );
                lastUsedImage->addMetaInfo( "_NCalcMean_", QString("%1 ms").arg(sumCalc/numLoops) );
                lastUsedImage->addMetaInfo( "_NPrepMin_", QString("%1 ms").arg(minPrep) );
                lastUsedImage->addMetaInfo( "_NPrepMax_", QString("%1 ms").arg(maxPrep) );
                lastUsedImage->addMetaInfo( "_NPrepMean_", QString("%1 ms").arg(sumPrep/numLoops) );
                lastUsedImage->addMetaInfo( "@", "" ); // Sortieren der Meta-Tabelle
                if ( saveImages )
                {
                    lastUsedImage->saveImage( saveBasePath + QString("-t=%1-hkl=%2.png").arg(ui->inpNumCores->value()).arg(hklmax[h]) );
                }
                qDebug() << "ERG T=" << ui->inpNumCores->value() << "HKLmax=" << hklmax[h]
                         << "Q=" << quadr
                         << "Calc=" << minCalc << sumCalc/numLoops << maxCalc
                         << "Prep=" << minPrep << sumPrep/numLoops << maxPrep;
                // LaTeX Ausgabe mit den Rundungen (am Ende auf der Konsole oder in eine Datei)
                if ( saveFilename.isEmpty() )
                    latex << QString(" & %1 & %2 & %7 & %3/ %4/ %5 (%6)").arg(ui->inpNumCores->value()).arg(hklmax[h])
                                 .arg(minCalc,0,'f',3).arg(sumCalc/numLoops,0,'f',3).arg(maxCalc,0,'f',3)
                                 .arg(sumPrep/numLoops,0,'f',3).arg(quadr);
                else
                {
                    QFile f(saveFilename);
                    if ( ! f.open(QIODevice::Append) )
                        qDebug() << saveFilename << f.errorString();
                    if ( f.isOpen() )
                    {
                        f.write( QString(" & %1 & %2 & %8 & %3/ %4/ %5 (%6)%7")
                                    .arg(ui->inpNumCores->value()).arg(hklmax[h])
                                    .arg(minCalc,0,'f',3).arg(sumCalc/numLoops,0,'f',3).arg(maxCalc,0,'f',3)
                                    .arg(sumPrep/numLoops,0,'f',3).arg(EOL).arg(quadr).toLatin1() );
                        f.close();
                    }
                }
            } // for h
            if ( _bAbbruch ) break;
        } // for t
        if ( _bAbbruch ) break;
    } // for q
    if ( _bAbbruch )
    {
        qDebug() << "Aborted." << saveFilename;
        if ( ! saveFilename.isEmpty() )
        {
            QFile f(saveFilename);
            if ( ! f.open(QIODevice::Append) )
                qDebug() << saveFilename << f.errorString();
            if ( f.isOpen() )
            {
                f.write( "Aborted." EOL );
                f.close();
            }
        }
    }
    else
    {
        qDebug() << "Complete" << saveFilename;
        if ( latex.size() > 0 ) qDebug() << latex;
    }
    ui->togIgnoreUpdates->setChecked( savUpd );
    ui->togExpandImage->setChecked( savExpand );
    ui->radQ1->setChecked( savQ1 );
    ui->radQ2->setChecked( savQ2 );
    ui->radQ4->setChecked( savQ4 );
    bIgnoreUpdates = savUpd;
}

void SC_MainGUI::on_butAbbruch_clicked()
{
    _bAbbruch = true;
    ui->butAbbruch->setEnabled(false);
    ui->butAbbruch->setText("Aborting...");
    qApp->processEvents();
}

bool SC_MainGUI::myProgressAndAbort( int val )
{
    if ( current ) current->updateProgBar(val);
    return _bAbbruch;
}

void SC_MainGUI::updateProgBar( int val )
{
    if ( val < 0 )
    {   // Signals the usage of the GPU
        ui->progressBar->setFormat("-GPU-");
        ui->progressBar->setValue(100);
        qApp->processEvents();
    }
    else if ( !bIgnoreUpdates )
    {   // Normal percentage display
        if ( ui->progressBar->value() == val ) return; // optimize GUI updates
        ui->progressBar->setValue(val);
        qApp->processEvents();
    }
}


void SC_MainGUI::on_radNewImageCfg_toggled(bool checked)
{
    if ( !checked ) return;
    ui->radNewImageCal->blockSignals(true);
    ui->radNewImageCal->setChecked(true);
    ui->radNewImageCal->blockSignals(false);
}

void SC_MainGUI::on_radLastImageCfg_toggled(bool checked)
{
    if ( !checked ) return;
    ui->radLastImageCal->blockSignals(true);
    ui->radLastImageCal->setChecked(true);
    ui->radLastImageCal->blockSignals(false);
}

void SC_MainGUI::on_radNewImageCal_toggled(bool checked)
{
    if ( !checked ) return;
    ui->radNewImageCfg->blockSignals(true);
    ui->radNewImageCfg->setChecked(true);
    ui->radNewImageCfg->blockSignals(false);
}

void SC_MainGUI::on_radLastImageCal_toggled(bool checked)
{
    if ( !checked ) return;
    ui->radLastImageCfg->blockSignals(true);
    ui->radLastImageCfg->setChecked(true);
    ui->radLastImageCfg->blockSignals(false);
}


void SC_MainGUI::on_actionSave_all_Parameters_triggered()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString fn = data.value("LastParam",".").toString();
    fn = QFileDialog::getSaveFileName( this, "Save Parameter", fn, "Parameter (*.ini)", nullptr, QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".ini",Qt::CaseInsensitive) ) fn += ".ini";
    data.setValue("LastParam",fn);
    performSaveParamOperation( fn );
} /* on_actionSave_all_Parameters_triggered() */

void SC_MainGUI::performSaveParamOperation( QString fn )
{
    // Save all method specific parameters
    calcGui->saveParameter(fn);
    // Save all common parameters
    QSettings sets( fn, QSettings::IniFormat );
    sets.beginGroup( "Inputs" );
    sets.setValue( "RadioButtonQ1",  ui->radQ1->isChecked() );
    sets.setValue( "RadioButtonQ2",  ui->radQ2->isChecked() );
    sets.setValue( "RadioButtonQ4",  ui->radQ4->isChecked() );
    sets.setValue( "ExpandImage",    ui->togExpandImage->isChecked() );
    sets.setValue( "EditAxis1x",     ui->inpAx1->value() );
    sets.setValue( "EditAxis1y",     ui->inpAy1->value() );
    sets.setValue( "EditAxis1z",     ui->inpAz1->value() );
    sets.setValue( "EditAxis2x",     ui->inpAx2->value() );
    sets.setValue( "EditAxis2y",     ui->inpAy2->value() );
    sets.setValue( "EditAxis2z",     ui->inpAz2->value() );
    sets.setValue( "EditAxis3x",     ui->inpAx3->value() );
    sets.setValue( "EditAxis3y",     ui->inpAy3->value() );
    sets.setValue( "EditAxis3z",     ui->inpAz3->value() );
    sets.setValue( "Editxrel",       ui->inpN1->value() );
    sets.setValue( "Edityrel",       ui->inpN2->value() );
    sets.setValue( "Editzrel",       ui->inpN3->value() );
    sets.setValue( "Editdom1",       ui->inpSigX->value() );
    sets.setValue( "Editdom2",       ui->inpSigY->value() );
    sets.setValue( "Editdom3",       ui->inpSigZ->value() );
    sets.setValue( "Editx1rel",      ui->inpU1->value() );
    sets.setValue( "Edity1rel",      ui->inpU2->value() );
    sets.setValue( "Editz1rel",      ui->inpU3->value() );
    sets.setValue( "Editx2rel",      ui->inpV1->value() );
    sets.setValue( "Edity2rel",      ui->inpV2->value() );
    sets.setValue( "Editz2rel",      ui->inpV3->value() );
    sets.setValue( "Edithklmax",     ui->inpHKLmax->value() );
    sets.setValue( "EditGridPoints", ui->inpGridPoints->value() );
    sets.setValue( "UseLattice3D",   ui->grpLattice3D->isChecked() );
    sets.setValue( "Threads",        ui->inpNumCores->value() );
    sets.setValue( "CurMethod",      ui->tabMethods->tabText( ui->tabMethods->currentIndex() ) );
    sets.setValue( "EditCenterX",    ui->inpBCenterX->value() );
    sets.setValue( "EditCenterY",    ui->inpBCenterY->value() );
    sets.endGroup();
    sets.beginGroup( "AI" );
    sets.setValue( "LastSubDir", ui->inpSubDir->text() );
    sets.setValue( "Grayscale", ui->cbsAIoutputFormat->currentIndex() );
    sets.setValue( "FileInputEna", ui->grpFileInput->isChecked() );
    sets.setValue( "FileInputLast", ui->inpFileName->text() );
    sets.setValue( "FileClass", ui->inpFileClass->text() );
    sets.setValue( "GenerateIFFT", ui->togAIUseFFTOutput->isChecked() );
    sets.setValue( "LinOut", ui->radAILinOutput->isChecked() );
    sets.setValue( "ScaleOut", ui->togAIScaleOutput->isChecked() );
    sets.endGroup();
} /* performSaveParamOperation() */


void SC_MainGUI::on_actionLoad_all_Parameters_triggered()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString fn = data.value("LastParam",".").toString();
    fn = QFileDialog::getOpenFileName( this, "Load Parameter", fn,
                                       "Parameter (*.ini);;Scatter-Parameter (*.par)",
                                       nullptr, QFileDialog::DontUseNativeDialog | QFileDialog::DontResolveSymlinks );
    if ( fn.isEmpty() ) return;
    data.setValue("LastParam",fn);
    local_Load_all_Parameters(fn,"");
}

void SC_MainGUI::on_actionLoad_only_current_parameters_triggered()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString fn = data.value("LastParam",".").toString();
    fn = QFileDialog::getOpenFileName( this, "Load Parameter for "+curMethod, fn,
                                       "Parameter (*.ini)",
                                       nullptr, QFileDialog::DontUseNativeDialog | QFileDialog::DontResolveSymlinks );
    if ( fn.isEmpty() ) return;
    data.setValue("LastParam",fn);
    local_Load_all_Parameters(fn,curMethod);
} /* on_actionLoad_only_current_parameters_triggered() */

QString SC_MainGUI::local_Load_all_Parameters(QString fn, QString onlyMethod)
{
    bIgnoreRecalc = true;  // Load start
    if ( fn.endsWith(".par",Qt::CaseInsensitive) )
    {
        qDebug() << fn;
        loadScatterParameter(fn);
        bIgnoreRecalc = false;  // Load Scatter done
        return "";
    }
    ui->lblLoadPrompt->show();
    ui->lblLoadFilename->show();
    ui->lblLoadFilename->setText(fn);
    // Auf neues Image umschalten
    ui->radNewImageCal->setChecked(true);
    // Load all method specific parameters
    bool isloaded = false;
    QString rv="";
    if ( onlyMethod.isEmpty() )
    {
        // Wenn im aktuellen Entwicklungsschritt (nur die Methode "Generic" bekannt) alte Files geladen werden
        // sollen, die noch mehrere Methoden enthalten, dann gibt es den Schlüssel "CurMethod" in der Gruppe
        // "Inputs". Dies soll dann die Methode sein, die geladen werden soll, auch wenn dieser Name hier nicht
        // bekannt ist. Zugleich wird der Wert vom Parameter LType auf die Methode gesetzt. Somit können die
        // Parameterfiles aus dem alten Format wieder normal gelesen werden.
        QSettings sets( fn, QSettings::IniFormat );
        QStringList gr = sets.childGroups();
        gr.removeOne("Inputs");
        gr.removeOne("AI");     // TODO: neue feste Gruppen hier ausblenden
        if ( gr.size() > 1 )
        {   // Jetzt gibt es mehr als eine Methode
            sets.beginGroup( "Inputs" );
            QString m = sets.value("CurMethod","").toString();
            sets.endGroup();
            if ( ! m.isEmpty() )
            {
                QString mm = m;
                if ( mm.indexOf(" ") > 0 ) mm.truncate(mm.indexOf(" "));
                //qDebug() << "Load spec" << m << mm;
                int val = 0;
                while ( true )
                {
                    calcGui->updateParamValue( "", "LType", val, Qt::black, false );
                    QString tmp = calcGui->currentParamValueStr( "", "LType", true ); // hier in aktueller Methode suchen...
                    //qDebug() << "    " << val << tmp;
                    if ( tmp.isEmpty() || tmp[0] == '?' ) break; // Ende der Liste oder Fehlermeldung
                    if ( tmp.startsWith(mm,Qt::CaseInsensitive) ) break; // gefunden
                    val++;
                }
                rv = calcGui->loadParameter(fn,"@"+m);
                isloaded = true;
            }
        }
    } // if ( onlyMethod.isEmpty() )
    if ( !isloaded )
        rv = calcGui->loadParameter(fn,onlyMethod);

    if ( onlyMethod.isEmpty() )
    {
        // Load all common parameters
        QSettings sets( fn, QSettings::IniFormat );
        sets.beginGroup( "Inputs" );
        ui->radQ1->setChecked( sets.value( "RadioButtonQ1", false ).toBool() );
        ui->radQ2->setChecked( sets.value( "RadioButtonQ2", false ).toBool() );
        ui->radQ4->setChecked( sets.value( "RadioButtonQ4", false ).toBool() );
        ui->togExpandImage->setChecked( sets.value( "ExpandImage", false ).toBool() );
        ui->togExpandImage->setEnabled( ! ui->radQ4->isChecked() );
        ui->inpAx1->setValue( sets.value( "EditAxis1x", 0.0 ).toDouble() );
        ui->inpAy1->setValue( sets.value( "EditAxis1y", 0.0 ).toDouble() );
        ui->inpAz1->setValue( sets.value( "EditAxis1z", 0.0 ).toDouble() );
        ui->inpAx2->setValue( sets.value( "EditAxis2x", 0.0 ).toDouble() );
        ui->inpAy2->setValue( sets.value( "EditAxis2y", 0.0 ).toDouble() );
        ui->inpAz2->setValue( sets.value( "EditAxis2z", 0.0 ).toDouble() );
        ui->inpAx3->setValue( sets.value( "EditAxis3x", 0.0 ).toDouble() );
        ui->inpAy3->setValue( sets.value( "EditAxis3y", 0.0 ).toDouble() );
        ui->inpAz3->setValue( sets.value( "EditAxis3z", 0.0 ).toDouble() );
        ui->inpN1->setValue( sets.value( "Editxrel", 0.0 ).toDouble() );
        ui->inpN2->setValue( sets.value( "Edityrel", 0.0 ).toDouble() );
        ui->inpN3->setValue( sets.value( "Editzrel", 0.0 ).toDouble() );
        ui->inpSigX->setValue( sets.value( "Editdom1", 0.0 ).toDouble() );
        ui->inpSigY->setValue( sets.value( "Editdom2", 0.0 ).toDouble() );
        ui->inpSigZ->setValue( sets.value( "Editdom3", 0.0 ).toDouble() );
        ui->inpU1->setValue( sets.value( "Editx1rel", 0.0 ).toDouble() );
        ui->inpU2->setValue( sets.value( "Edity1rel", 0.0 ).toDouble() );
        ui->inpU3->setValue( sets.value( "Editz1rel", 0.0 ).toDouble() );
        ui->inpV1->setValue( sets.value( "Editx2rel", 0.0 ).toDouble() );
        ui->inpV2->setValue( sets.value( "Edity2rel", 0.0 ).toDouble() );
        ui->inpV3->setValue( sets.value( "Editz2rel", 0.0 ).toDouble() );
        ui->inpHKLmax->setValue( sets.value( "Edithklmax", 3 ).toInt() );
        ui->inpGridPoints->setValue( sets.value( "EditGridPoints", 100 ).toInt() );
        ui->grpLattice3D->setChecked( sets.value( "UseLattice3D", false ).toBool() );
        // Wenn im Parametersatz die GPU aktiviert ist, im aktuellen System diese aber nicht vorhanden ist,
        // so würde das Eingabefeld auf 1 Thread schalten. Das ist unklug und wird hier auf die maximale
        // Anzahl von Threads umgebogen. Kommt vor, wenn Parameterfiles vom Linux (mit GPU) zum Laptop
        // (ohne GPU) kopiert werden.
        int thr = sets.value( "Threads", 0 ).toInt();
        if ( thr == 0 /*GPU*/ && ui->inpNumCores->minimum() == 1 /*keine GPU vorhanden*/ )
            thr = ui->inpNumCores->maximum();
        ui->inpNumCores->setValue( thr );
        ui->inpBCenterX->setValue( sets.value( "EditCenterX", 0 ).toDouble() );
        ui->inpBCenterY->setValue( sets.value( "EditCenterY", 0 ).toDouble() );
        if ( ui->tabMethods->count() > 1 )
        {
            QString m = sets.value( "CurMethod", "" ).toString();
            if ( ! m.isEmpty() )
            {
                for ( int i=0; i<ui->tabMethods->count(); i++ )
                    if ( ui->tabMethods->tabText(i) == m )
                    {
                        ui->tabMethods->setCurrentIndex(i);
                        //on_tabMethods_currentChanged(i);
                        break;
                    }
            }
        }
        sets.endGroup();
        sets.beginGroup( "AI" );
        ui->inpSubDir->setText( sets.value( "LastSubDir", dataPath ).toString() );
        ui->cbsAIoutputFormat->setCurrentIndex( sets.value( "Grayscale", 1 ).toInt() );
        ui->grpFileInput->setChecked( sets.value("FileInputEna",false).toBool() );
        ui->inpFileName->setText( sets.value("FileInputLast","").toString() );
        ui->inpFileClass->setText( sets.value("FileClass","{M}").toString() );
        ui->togAIUseFFTOutput->setChecked( sets.value("GenerateIFFT",false).toBool() );
        ui->radAILinOutput->setChecked(   sets.value("LinOut").toBool() );
        ui->radAILogOutput->setChecked( ! sets.value("LinOut").toBool() );
        ui->togAIScaleOutput->setChecked( sets.value("ScaleOut").toBool() );
        sets.endGroup();
    } // if ( onlyMethod.isEmpty() )
    bIgnoreRecalc = false;  // Load done
    return rv;
} /* local_Load_all_Parameters() */


widImage* SC_MainGUI::addImage( bool neu, int x0, int x1, int y0, int y1, double *d, QString title, bool meta )
{
    if ( d == nullptr ) return nullptr;
    widImage *img = nullptr;
    if ( !neu && ui->radLastImageCfg->isChecked() && !ui->togOnlyNewWindow->isChecked() )
    {
        img = nullptr;
        for ( int i=images.size()-1; i>=0; i-- )
            if ( images[i]->isVisible() )
            {
                img = images[i];
                break;
            }
    }
    if ( img == nullptr )
    {
        numberCounter++;    // Global counter, no reset
#ifdef IMG_ONE_WINDOW
        QMdiSubWindow *tmp = new QMdiSubWindow(ui->mdiArea);
        img = new widImage( QString("%1 #%2").arg(title).arg(numberCounter), dataPath );
        tmp->setWidget(img);
        qApp->processEvents();
#else
        img = new widImage( QString("%1 #%2").arg(title).arg(numberCounter), dataPath );
#endif
        connect( img, SIGNAL(destroyed(QObject*)), this, SLOT(imageWindowClosed(QObject*)) );
        ui->radLastImageCfg->setChecked( !ui->togOnlyNewWindow->isChecked() );
        images.append(img);
        neu = true; // for the positioning of the new image (no last img available)
    }
    //else
    //    img->addMetaInfo("",""); // Alles löschen

    ui->butSaveImagesAI->setEnabled( images.size() > 0 );
    ui->butIFFT->setEnabled( images.size() > 0 );  // Add Img

    img->setData( x0, x1, y0, y1, d );
    if ( meta )
    {   // Meta-Daten übertragen
        // Methodenspezifische Metadaten
        img->addMetaInfo( "_Calculation_", calcGui->curMethod->subCalc->methodName() );
        QHash<QString,paramHelper*>::iterator ip = calcGui->curMethod->params.begin();
        while ( ip != calcGui->curMethod->params.end() )
        {
            img->addMetaInfo( ip.key(), //ip.value()->str );
                              calcGui->currentParamValueStr( calcGui->curMethod->subCalc->methodName(),
                                                             ip.key(), true ) );
            ++ip;
        }
        // Globale Eingabefelder
        QHash<QString,Double3>::iterator iv = calcGui->inpVectors.begin();
        while ( iv != calcGui->inpVectors.end() )
        {
            img->addMetaInfo( iv.key(), iv.value().toString() );
            ++iv;
        }
        QHash<QString,double>::iterator ii = calcGui->inpValues.begin();
        while ( ii != calcGui->inpValues.end() )
        {
            if ( ii.key().startsWith("RadioButtonQ") || ii.key().startsWith("ExpandImage") )
                img->addMetaInfo( ii.key(), ii.value() ? "True" : "False" );
            else if ( ii.key().startsWith("_CalcTime_") ||
                      ii.key().startsWith("_PrepTime_") ||
                      ii.key().startsWith("_NCalc"    ) ||
                      ii.key().startsWith("_NPrep"    ) )
                img->addMetaInfo( ii.key(), QString("%1 ms").arg(ii.value()) );
            else
                img->addMetaInfo( ii.key(), QString::number(ii.value()) );
            ++ii;
        }
        img->addMetaInfo("@",""); // Sortieren
    }
    img->show();
    qApp->processEvents();
#ifdef IMG_ONE_WINDOW
    img->parentWidget()->adjustSize();
#endif
    if ( neu && ui->togAutoPosit->isChecked() )
    {   // Posit new images automatically (try to...)
        if ( img->width() + imgPosX + 10 >= qApp->desktop()->width() )
        {
            imgPosX = imgPosX0;
            imgPosY += img->height();
            if ( imgPosY + 20 >= qApp->desktop()->height() )
            {
                imgPosX0 += 20;
                imgPosY0 += 20;
                imgPosX  = imgPosX0;
                imgPosY  = imgPosY0;
            }
        }
        //qDebug() << "MOVE addImg" << imgPosX << imgPosY << img->windowTitle();
        img->move( imgPosX, imgPosY );
        imgPosX += img->width();
    }
    // Tab Configuration: listWidget
    // Tab FFT: cbsFFTWindows
    // Tab Fit: cbsImageWindows
    QList<QListWidgetItem*> items = ui->lisDataWindows->findItems( img->windowTitle(), Qt::MatchStartsWith );
    if ( items.size() == 0 )
    {
        ui->lisDataWindows->addItem( img->windowTitle()+" ("+img->getMetaTitle()+")" );
        ui->cbsFFTWindows->addItem( img->windowTitle() );
        ui->cbsFitImageWindows->addItem( img->windowTitle() );
    }
    //ui->butSaveAllImages->setEnabled( ui->lisDataWindows->count() > 0 );
    return img;
} /* addImage() */


/**
 * @brief SasMain::imageWindowClosed
 * @param o - object to be closed (class widImage)
 * One image window will be closed, remove it from the internal list of images
 * and remove it from the selection box for the fit.
 */
void SC_MainGUI::imageWindowClosed( QObject *o )
{
    // if the main window will be closed, then ignore the rest
    if ( closeMainActive ) return;
    // sometimes the parmeter will be zero, then use the sender object
    if ( o == nullptr ) o = sender();
    // search the internal list of images
    for ( int i=0; i<images.size(); i++ )
        if ( images.at(i) == o )
        {
            // Remove the item in the list
            QList<QListWidgetItem*> items = ui->lisDataWindows->findItems( images.at(i)->windowTitle(), Qt::MatchStartsWith );
            if ( items.size() == 1 )
            {
                ui->lisDataWindows->takeItem( ui->lisDataWindows->row(items[0]) );
                delete items[0];
                //ui->butSaveAllImages->setEnabled( ui->lisDataWindows->count() > 0 );
                // and find the index in the selection boxes
                int id = ui->cbsFitImageWindows->findText( images.at(i)->windowTitle() );
                ui->cbsFitImageWindows->removeItem( id );
                ui->cbsFFTWindows->removeItem( id ); // same list
                // if this was the last image, disable the fit button
                if ( ui->cbsFitImageWindows->count() == 0 )
                {
                    ui->butFitStart->setEnabled(false);
                    ui->butIFFT->setEnabled(false); // Img Window Close
                }
                // if there was only one image left (or nothing) activate the NewImage radio button
                if ( ui->cbsFitImageWindows->count() <= 1 )
                    ui->radNewImageCfg->setChecked(true);
            }
            // remove the image window from the internal list
            images.removeAt(i);
            ui->butSaveImagesAI->setEnabled( images.size() > 0 );
            if ( images.size() == 0 )
            {   // No images, so reset the positions
                imgPosX = imgPosX0 = 0;
                imgPosY = imgPosY0 = 0;
                numberCounter = 0;
            }
            else if ( images.size() == 1 )
            {   // Only one image, adjust the positions
                imgPosX = imgPosX0 = images[0]->pos().x() + images[0]->width();
                imgPosY = imgPosY0 = images[0]->pos().y();
                numberCounter = 1;
            }
            return;
        }
} /* imageWindowsClosed() */


/**
 * @brief SC_MainGUI::on_listWidget_currentTextChanged [SLOT]
 * @param currentText - clicked text in the list
 * Searches the image window in the internal list and raise the window.
 */
void SC_MainGUI::on_lisDataWindows_currentTextChanged(const QString &currentText)
{
    ui->butDataFindCenter->setEnabled( currentText.startsWith("KWS") );
    ui->butDataSetMask->setEnabled( currentText.startsWith("KWS") && ui->lisDataWindows->count() > 1 );
    ui->butDataCopyScaling->setEnabled( ui->lisDataWindows->count() > 1 );
    ui->grpExtractImage->setEnabled( ui->lisDataWindows->count() > 0 );
    ui->grpNoFitRegions->setEnabled( ui->lisDataWindows->count() > 0 );
    if ( lastSelectedImage != nullptr )
    {
        lastSelectedImage->disconnect( SIGNAL(extractUpdateRect(QRect)) );
        lastSelectedImage->enableExtract(false);
    }
    for ( int i=0; i<images.size(); i++ )
        if ( currentText.startsWith(images.at(i)->windowTitle()) )
        {
            lastSelectedImage = images.at(i);
            connect( lastSelectedImage, SIGNAL(extractUpdateRect(QRect)), this, SLOT(extractUpdateRect(QRect)) );
            on_grpExtractImage_toggled( ui->grpExtractImage->isChecked() );
            on_grpNoFitRegions_toggled( ui->grpNoFitRegions->isChecked() );
            lastSelectedImage->raise();
            copyMetaToLatticeTable( lastSelectedImage );
            // Jetzt noch eventuelle NoFitRect-Metadaten in die GUI-Tabelle kopieren
            QRect rc;
            ui->tblNoFitRegions->blockSignals(true);
            for ( int r=0; r<4; r++ )
            {
                rc = lastSelectedImage->getNoFitRect(r);
                ui->tblNoFitRegions->setItem( r, 0, new QTableWidgetItem(QString::number(rc.left())) );
                ui->tblNoFitRegions->setItem( r, 1, new QTableWidgetItem(QString::number(rc.right())) );
                ui->tblNoFitRegions->setItem( r, 2, new QTableWidgetItem(QString::number(rc.top())) );
                ui->tblNoFitRegions->setItem( r, 3, new QTableWidgetItem(QString::number(rc.bottom())) );
            }
            ui->tblNoFitRegions->blockSignals(false);
            return;
        }
    lastSelectedImage = nullptr;
} /* on_listWidget_currentTextChanged() */


/**
 * @brief SC_MainGUI::on_togIgnoreUpdates_toggled [SLOT]
 * @param checked - current state of the check box
 * This saves the current state in a global variale.
 */
void SC_MainGUI::on_togIgnoreUpdates_toggled(bool checked)
{
    bIgnoreUpdates = checked;
}


void SC_MainGUI::on_butIFFT_clicked()
{
    QString tit = ui->cbsFFTWindows->currentText();
    int curimg = -1;
    for ( int i=0; i<images.size(); i++ )
        if ( tit.compare(images.at(i)->windowTitle()) == 0 )
        {
            curimg = i;
            break;
        }
    if ( curimg < 0 ) return;

    SasCalc_PostProc::inst()->setLogging( true );  // besser für Debug
    performIFFTcalculations(curimg,tit,false);

    if ( !bIgnoreRecalc && autoProcessingFile!=nullptr && autoProcessingFile->isOpen() )
    {
        QTimer::singleShot( 200, this, SLOT(autoProcessingTimer()) );
    }
}

void SC_MainGUI::on_butFFT_clicked()
{
    QString tit = ui->cbsFFTWindows->currentText();
    int curimg = -1;
    for ( int i=0; i<images.size(); i++ )
        if ( tit.compare(images.at(i)->windowTitle()) == 0 )
        {
            curimg = i;
            break;
        }
    if ( curimg < 0 ) return;

    SasCalc_PostProc::inst()->setLogging( true );  // besser für Debug
    performIFFTcalculations(curimg,tit,true);

    if ( !bIgnoreRecalc && autoProcessingFile!=nullptr && autoProcessingFile->isOpen() )
    {
        QTimer::singleShot( 200, this, SLOT(autoProcessingTimer()) );
    }
}

void SC_MainGUI::performIFFTcalculations(int curimg, QString tit, bool foreward)
{
    QHash<QString,QString> metaData;

    int sx = images[curimg]->xmax() - images[curimg]->xmin();
    int sy = images[curimg]->ymax() - images[curimg]->ymin();

    int bsx = images[curimg]->getFileInfos()->centerX;
    int bsy = images[curimg]->getFileInfos()->centerY;
    if ( bsx == 0 || bsy == 0 )
    {
        bsx = sx / 2;
        bsy = sy / 2;
    }

    //qDebug() << "performIFFTcalculations:"
    //         << "X" << images[curimg]->getFileInfos()->li.centerX << images[curimg]->xmin()
    //         << "Y" << images[curimg]->getFileInfos()->li.centerY << images[curimg]->ymin();

    double *data = nullptr;

    if ( ui->radFFTLogInput->isChecked() )
    {
        const double *src = images[curimg]->dataPtr();
        data = new double[sx*sy];
        //qDebug() << "r,phi new data" << sx*sy;
        double vlogmin=1e26, vlogmax=0;
        for ( int i=0, y=images[curimg]->ymin(); y<images[curimg]->ymax(); y++ )
            for ( int x=images[curimg]->xmin(); x<images[curimg]->xmax(); x++, i++ )
                if ( (images[curimg]->xmin() < 0 && x != 0 && y != 0) ||
                     (images[curimg]->xmin() == 0) )
                {   // Wenn der Nullpunkt im Bild ist, dann ist dort i.d.R. ein Extremwert (z.B. -7.4e+25 bei FCC)
                    // Bei den anderen Bildern (r,phi oder FFT) ist immer der Nullpunkt in einer Ecke
                    if ( src[i] > 0 )
                    {
                        if ( log10(src[i]) < vlogmin )
                            vlogmin = log10(src[i]);
                        else if ( log10(src[i]) > vlogmax )
                            vlogmax = log10(src[i]);
                    }
                }
        metaData.insert( "LogMin", QString::number(vlogmin) );
        metaData.insert( "LogMax", QString::number(vlogmax) );
        for ( int i=0; i<sx*sy; i++ )
            data[i] = (log10(src[i]) - vlogmin) / (vlogmax-vlogmin);
    }

    metaData.insert( "From Image", tit );
    metaData.insert( "BeamPosX", QString::number(bsx) );
    metaData.insert( "BeamPosY", QString::number(bsy) );

    int srphi = ui->cbsFFTsizeRphi->currentText().left(3).trimmed().toInt();
    int sout  = ui->cbsFFTsizeOut->currentText().left(3).trimmed().toInt();
    // TODO: was ist bei unterschiedlichen Größen?

    double *outr = nullptr;

    if ( ui->grpFFTuseRphi->isChecked() )
    {   // Jetzt soll das (r,phi) verwendet werden.

        // Das (r,phi) Bild berechnen
        auto anfr = std::chrono::high_resolution_clock::now();
        outr = SasCalc_PostProc::inst()->generateRphi(0, sx, 0, sy, bsx, bsy, srphi,
                                                      (data==nullptr)?images[curimg]->dataPtr():data,
                                                      ui->togFFTscaleRphi->isChecked(),
                                                      ui->togFFTclipRphi->isChecked(),
                                                      ui->togFFTclip40Rphi->isChecked() );
        auto endr = std::chrono::high_resolution_clock::now();
        auto timr = std::chrono::duration_cast<std::chrono::duration<float>>(endr-anfr);
        //qDebug() << "Timing (r,phi):" << timr.count()*1000.0 << "ms";


        if ( ui->togFFTdispRphi->isChecked() && outr != nullptr )
        {
            widImage* img = addImage( true, 0, srphi, 0, srphi, outr, "(r,phi)-Image", false );
            QHash<QString,QString>::iterator ii = metaData.begin();
            while ( ii != metaData.end() )
            {
                img->addMetaInfo( ii.key(), ii.value() );
                ++ii;
            }
            img->addMetaInfo( "NbRows", QString::number(srphi) );
            img->addMetaInfo( "NbCols", QString::number(srphi) );
            img->addMetaInfo( "_CalcTime_", QString("%1 ms").arg(timr.count()*1000.0) );
            img->addMetaInfo("@",""); // Sortieren
        } // if ( ui->togFFTdispRphi->isChecked() )
    } // if ( ui->grpFFTuseRphi->isChecked() )

    //if ( outr != nullptr ) return; // TODO: nur zum Test der r,phi Routine

    // und jetzt das IFFT Bild
    const double *fftinp;
    int x0, x1, y0, y1;
    if ( outr != nullptr )
    {
        fftinp = outr; // r,phi
        x0 = 0;
        x1 = srphi;
        y0 = 0;
        y1 = srphi;
    }
    else if ( data != nullptr )
    {
        fftinp = data; // Log
        x0 = images[curimg]->xmin();
        x1 = images[curimg]->xmax();
        y0 = images[curimg]->ymin();
        y1 = images[curimg]->ymax();
    }
    else
    {
        fftinp = images[curimg]->dataPtr(); // Lin
        x0 = images[curimg]->xmin();
        x1 = images[curimg]->xmax();
        y0 = images[curimg]->ymin();
        y1 = images[curimg]->ymax();
    }
    double *outi = nullptr;
    SasCalc_PostProc::_outType ot;
    if ( ui->radFFToutReal->isChecked() )
    {
        ot = SasCalc_PostProc::outReal;
        metaData.insert( "FFToutput", "Real values" );
    }
    else if ( ui->radFFToutImag->isChecked() )
    {
        ot = SasCalc_PostProc::outImag;
        metaData.insert( "FFToutput", "Imag values" );
    }
    else if ( ui->radFFToutBetrag->isChecked() )
    {
        ot = SasCalc_PostProc::outAbs;
        metaData.insert( "FFToutput", "Absolute values" );
    }
    else //if ( ui->radFFToutSpectrum->isChecked() )
    {
        ot = SasCalc_PostProc::outSpec;
        metaData.insert( "FFToutput", "Power spectrum" );
    }

    auto anfi = std::chrono::high_resolution_clock::now();
    outi = SasCalc_PostProc::inst()->calculateIFFT(foreward,
                                                   x0, x1, y0, y1, fftinp, sout, ot,
                                                   ui->togFFTscaleOut->isChecked(),
                                                   ui->togFFTclipOut->isChecked(),
                                                   ui->togFFTclip40Out->isChecked(),
                                                   ui->togIFFTSwap->isChecked() );

    auto endi = std::chrono::high_resolution_clock::now();
    auto timi = std::chrono::duration_cast<std::chrono::duration<float>>(endi-anfi);
    //qDebug() << "Timing (IFFT):" << timi.count()*1000.0 << "ms";
    if ( outi != nullptr )
    {
        widImage* img = addImage( true, 0, sout, 0, sout, outi, "iFFT-Image", false );
        QHash<QString,QString>::iterator ii = metaData.begin();
        while ( ii != metaData.end() )
        {
            img->addMetaInfo( ii.key(), ii.value() );
            ++ii;
        }
        img->addMetaInfo( "NbRows", QString::number(sout) );
        img->addMetaInfo( "NbCols", QString::number(sout) );
        img->addMetaInfo( "_CalcTime_", QString("%1 ms").arg(timi.count()*1000.0) );
        img->addMetaInfo("@",""); // Sortieren
    }
    // Next is allways a new image ...
    ui->radNewImageCfg->setChecked(true);
}

void SC_MainGUI::on_inpGridPoints_valueChanged(int arg1)
{
    Q_UNUSED(arg1)
    /* TODO
    if ( ui->radQ1->isChecked() )
    {   // Nur 1 Quadrant -> BS auf ?
    }
    else if ( ui->radQ2->isChecked() )
    {   // Nur 2 Quadranten -> BS auf ?
    }
    else // if ( ui->radQ4->isChecked() )
    {   // Alle 4 Quadranten -> BS auf ?
    }
    ui->inpBCenterX->setValue( arg1/2. );
    ui->inpBCenterY->setValue( arg1/2. );
    */
}

void SC_MainGUI::on_butUseQMax_clicked()
{
    double val = ui->inpQMax->text().mid( ui->inpQMax->text().indexOf("=")+1 ).trimmed().toDouble();
    calcGui->setParamToAllMethods( "EditQmax", val );
}

void SC_MainGUI::on_tblLattice3DValues_cellChanged(int /*row*/, int /*column*/)
{
    // Berechne qmax (siehe sascalc_fcc_gpu.cu: doIntCalc_FCC_F)
    int cols = ui->tblLattice3DValues->item(tblLattNbCols,0)->text().toInt();
    int rows = ui->tblLattice3DValues->item(tblLattNbRows,0)->text().toInt();
    double posx = ui->tblLattice3DValues->item(tblLattBeamX,0)->text().toDouble();
    double posy = ui->tblLattice3DValues->item(tblLattBeamY,0)->text().toDouble();
    double wlen_m = ui->tblLattice3DValues->item(tblLattWaveLen,0)->text().toDouble()*1.0E-9; /*nm->m*/
    double dist_m = ui->tblLattice3DValues->item(tblLattSaDist,0)->text().toDouble();
    double pixx_m = ui->tblLattice3DValues->item(tblLattPixelX,0)->text().toDouble()*1.0E-3; /*mm->m*/
    double pixy_m = ui->tblLattice3DValues->item(tblLattPixelY,0)->text().toDouble()*1.0E-3; /*mm->m*/

    int ihex_max = qMax( cols - posx, posx );
    int i_max    = qMax( rows - posy, posy );
    double qmax_cols = ( 2. * M_PI / (wlen_m * dist_m) ) * ihex_max * pixx_m * 1.0E-9; /* wieder in nm zurück */
    double qmax_rows = ( 2. * M_PI / (wlen_m * dist_m) ) * i_max    * pixy_m * 1.0E-9; /* wieder in nm zurück */

    ui->inpQMax->setText( QString("Qmax = %1").arg(qMax(qmax_cols,qmax_rows)) );
}

void SC_MainGUI::on_butOpenMeasFile_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    QString fn = data.value("LastImage",dataPath).toString();
    //qDebug() << fn;
    // Es muss der aktuelle Filter des letzten Images bestimmt werden, dann wird dieser auch selektiert
    // und alle passenden Dateien werden angezeigt.
    QStringList filter;
    filter << "TIFF-Images (*.tif *.tiff)"
           << "KWS-Data (*.dat *.data)"
           << "ESRF/Klora images (*.edf)"
           << "Spreadsheet-Format (*.spr)";
#ifdef NOHDF5
    if ( fn.endsWith(".h5") || fn.endsWith(".hdf") || fn.endsWith(".hdf5") || fn.endsWith(".nxs") )
        fn = dataPath;
    filter << "2D SANS (*.xml *.csv *.dat)";    // Hier kein NXS (ist HDF5, Nexus)
#else
    filter << "HDF5 Files (*.h5 *.hdf *.hdf5)";
    filter << "2D SANS (*.xml *.csv *.nxs *.dat)";
#endif

    filter << "SasCrystal (*.dat *.csv *.txt)";
    filter << "All files (*.*)";
    QStringList tmp = filter.filter( "." + QFileInfo(fn).suffix(), Qt::CaseInsensitive );
    QString curFilter = "";
    if ( tmp.size() > 0 )
        curFilter = tmp.at(0);
    else
        curFilter = filter.last();
    fn = QFileDialog::getOpenFileName( this, "Load Data File Image", fn, filter.join(";;"), &curFilter,
                                       QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;
    data.setValue("LastImage",fn);
    if ( ! local_OpenMeasFile(fn,nullptr) )
        qDebug() << fn << "unknown format";
}

bool SC_MainGUI::local_OpenMeasFile( QString fn, widImage **imgout )
{
    qApp->setOverrideCursor(Qt::WaitCursor);
    widImage *img = nullptr;
    if ( fn.endsWith(".dat",Qt::CaseInsensitive) ||
         fn.endsWith(".csv",Qt::CaseInsensitive) ||
         fn.endsWith(".txt",Qt::CaseInsensitive) )
    {   // Special Test to read back previously saved date from this program
        int pos = fn.lastIndexOf(".");
        QString base = fn.left(pos);
        if ( QFileInfo::exists(base+".csv") &&
             QFileInfo::exists(base+".dat") &&
             QFileInfo::exists(base+".txt") )
        {
            img = SC_ReadData::readImageSasCrystal( myAddImage, fn );
        }
    }
    if ( img == nullptr )
    {
        // Routines in sc_readdata.cpp
        if ( fn.endsWith(".tif",Qt::CaseInsensitive) ||
            fn.endsWith(".tiff",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageTiff( myAddImage, fn );
        else if ( fn.endsWith(".dat",Qt::CaseInsensitive) ||
                  fn.endsWith(".data",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageKWSData( myAddImage, fn );
        else if ( fn.endsWith(".edf",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageEDF( myAddImage, fn );
        else if ( fn.endsWith(".spr",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageSpreadsheet( myAddImage, fn );
#ifndef NOHDF5
        else if ( fn.endsWith(".h5",Qt::CaseInsensitive) ||
                  fn.endsWith(".hdf",Qt::CaseInsensitive) ||
                  fn.endsWith(".hdf5",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageHDF5( myAddImage, fn, fitIsRunning, false/*swap vert*/ );
        else if ( fn.endsWith(".nxs",Qt::CaseInsensitive) )     // 2D SANS (ILL)
            img = SC_ReadData::readImageHDF5( myAddImage, fn, fitIsRunning, true/*swap vert*/ );
#endif
        else if ( fn.endsWith(".xml",Qt::CaseInsensitive) ||    // 2D SANS (ILL)
                  fn.endsWith(".csv",Qt::CaseInsensitive) ||    // 2D SANS (ILL)
                  fn.endsWith(".txt",Qt::CaseInsensitive) )     // 2D SANS (ILL)
            img = SC_ReadData::readImage2dSans( myAddImage, fn );
    }
    qApp->restoreOverrideCursor();
    if ( img == nullptr ) return false; // Lesefehler oder unbekannte Dateiendung

    if ( imgout != nullptr ) *imgout = img;

    // Jetzt werden noch Meta-Daten in die GUI kopiert. Das mache ich bei allen Datentypen,
    //  da ich nicht wissen kann, ob diese Meta-Infos nicht doch enthalten sind.
    copyMetaToLatticeTable( img );

    // Datenfiles werden nie bei neuen Berechnungen automatisch überschrieben!
    ui->radLastImageCfg->setChecked( false );
    ui->radNewImageCfg->setChecked( true );
    return true;
}

void SC_MainGUI::copyMetaToLatticeTable( widImage *img )
{
    ui->tblLattice3DValues->setItem( tblLattWaveLen, 0, new QTableWidgetItem(img->metaInfo("wavelength",
                                     calcGui->currentParamValueStr(curMethod,"EditWavelength",true))) );
    ui->tblLattice3DValues->setItem( tblLattSaDist, 0, new QTableWidgetItem(img->metaInfo("SampleDist",
                                     calcGui->currentParamValueStr(curMethod,"EditDet",true))) );
    ui->tblLattice3DValues->setItem( tblLattPixelX, 0, new QTableWidgetItem(img->metaInfo("Pixel_X",
                                     QString::number(calcGui->currentParamValueDbl(curMethod,"EditPixelX")*1000))) );
    ui->tblLattice3DValues->setItem( tblLattPixelY, 0, new QTableWidgetItem(img->metaInfo("Pixel_Y",
                                     QString::number(calcGui->currentParamValueDbl(curMethod,"EditPixelY")*1000))) );
    ui->tblLattice3DValues->setItem( tblLattNbRows, 0, new QTableWidgetItem(img->metaInfo("NbRows")) );
    ui->tblLattice3DValues->setItem( tblLattNbCols, 0, new QTableWidgetItem(img->metaInfo("NbCols")) );
    int row = img->metaInfo("NbRows").toInt();
    int col = img->metaInfo("NbCols").toInt();
    int bsx = col/2 - ui->inpBCenterX->value();
    int bsy = row/2 - ui->inpBCenterY->value();
    ui->tblLattice3DValues->setItem( tblLattBeamX, 0, new QTableWidgetItem(img->metaInfo("BeamPosX", QString::number(bsx))) );
    ui->tblLattice3DValues->setItem( tblLattBeamY, 0, new QTableWidgetItem(img->metaInfo("BeamPosY", QString::number(bsy))) );
}

void SC_MainGUI::on_butSaveAllImages_clicked()
{
    // TODO
}





void SC_MainGUI::on_cbsMethod_currentIndexChanged(const QString &arg1)
{
    ui->cbsVariables->clear();
    ui->cbsVariables->addItems( calcGui->paramsForMethod(arg1,true,true,false) ); // nur numerische Werte, incl. globales
}

void SC_MainGUI::on_butUseForAI_clicked()
{
    QString key = ui->cbsVariables->currentText();
    for ( int r=0; r<ui->tblListe->rowCount(); r++ )
        if ( ui->tblListe->item(r,0)->text() == key )
            return;

    // Add a new row
    int row = ui->tblListe->rowCount();
    ui->tblListe->setRowCount( row+1 );
    // and set the values
    ui->tblListe->setItem( row, 0, new QTableWidgetItem(key) );
    QString val = calcGui->currentParamValueStr( ui->cbsMethod->currentText(), key, true );  // ComboBox als Text
    ui->tblListe->setItem( row, 1, new QTableWidgetItem(val) );
    ui->tblListe->setItem( row, 2, new QTableWidgetItem(val) );
    ui->tblListe->setItem( row, 3, new QTableWidgetItem("0") );

    ui->tblListe->resizeColumnToContents(0);
    //ui->tblListe->sortByColumn(0,Qt::AscendingOrder);

    ui->butAIsaveBkgOp->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butAIcheck->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butAIstart->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butRemoveVar->setEnabled(false);
}

void SC_MainGUI::on_butSaveTable_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("AI");
    QString fn = data.value("LastSaveFile",".").toString();
    fn = QFileDialog::getSaveFileName( this, "Save AI infos", fn, "Infofiles (*.sas_ai)" );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".sas_ai") ) fn += ".sas_ai";
    QFile f(fn);
    if ( !f.open(QIODevice::WriteOnly) )
    {
        qDebug() << f.errorString();
        return;
    }
    data.setValue("LastSaveFile",fn);
    data.endGroup();
    for ( int r=0; r<ui->tblListe->rowCount(); r++ )
    {
        QString line = "";
        for ( int c=0; c<4; c++ )
            line += ui->tblListe->item(r,c)->text() + "|";
        f.write( qPrintable(line+"\r\n") );
    }
    f.close();
}

void SC_MainGUI::on_butLoadTable_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("AI");
    QString fn = data.value("LastSaveFile",".").toString();
    fn = QFileDialog::getOpenFileName( this, "Load AI infos", fn, "Infofiles (*.sas_ai)" );
    if ( fn.isEmpty() ) return;
    data.setValue("LastSaveFile",fn);
    data.endGroup();
    QFile f(fn);
    if ( !f.open(QIODevice::ReadOnly) )
    {
        qDebug() << f.errorString();
        return;
    }
    ui->tblListe->clearContents();
    ui->tblListe->setRowCount(0);
    while ( !f.atEnd() )
    {
        QStringList sl = QString(f.readLine()).split("|");
        if ( sl.count() < 4 ) continue;

        int row = ui->tblListe->rowCount();
        ui->tblListe->setRowCount( row+1 );
        for ( int c=0; c<4/*5*/; c++ )
            ui->tblListe->setItem( row, c, new QTableWidgetItem(sl[c]) );
    }
    f.close();
    ui->tblListe->resizeColumnToContents(0);

    ui->butAIsaveBkgOp->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butAIcheck->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butAIstart->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butRemoveVar->setEnabled(false);
}

void SC_MainGUI::on_butClearTable_clicked()
{
    ui->tblListe->clearContents();
    ui->tblListe->setRowCount(0);

    ui->butAIsaveBkgOp->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butAIcheck->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butAIstart->setEnabled( ui->tblListe->rowCount() > 0 );
    ui->butRemoveVar->setEnabled(false);
}

void SC_MainGUI::on_tblListe_cellClicked(int row, int column)
{
    Q_UNUSED(row)
    Q_UNUSED(column)
    ui->butRemoveVar->setEnabled(true);
}

void SC_MainGUI::on_butRemoveVar_clicked()
{
    ui->tblListe->removeRow( ui->tblListe->currentRow() );
    ui->butRemoveVar->setEnabled(false);
}

void SC_MainGUI::on_butFileName_clicked()
{
    QString fn;
    fn = QFileDialog::getExistingDirectory( this, "Directory for inputfiles", ui->inpFileName->text(), QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;
    ui->inpFileName->setText( fn );
}

void SC_MainGUI::on_butSelectDir_clicked()
{
    QString fn = ui->inpSubDir->text();
    fn = QFileDialog::getExistingDirectory( this, "Save files into", fn, QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;
    ui->inpSubDir->setText( fn );
}


void SC_MainGUI::on_butAIcheck_clicked()
{
    _loopDefinition rv = getLoopDefinition();
/**/
    qDebug() << "CHECK +++ Start";
    int cnt = 0;
    foreach ( QString k, rv.uniqueKeys() )
    {
        //QMultiHash<QString, _loopVariables>::iterator i = rv.find(k);
        QMultiMap<QString, _loopVariables>::iterator i = rv.find(k);
        while ( i != rv.end() && i.key() == k )
        {
            qDebug() << cnt << k << i.value();
            ++i;
            cnt++;
        }
    }
    //qDebug() << rv;
    qDebug() << "CHECK --- Ende" << cnt;
/**/
    QMessageBox::information( this, "Loop definitions",
                              "Classes: "+rv.uniqueKeys().join(", ")
                              +QString("\n==> %1 images").arg(rv.size()),
                              QMessageBox::Ok );
}

void SC_MainGUI::on_butAIsaveBkgOp_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("AI");
    QString fn = data.value("LastSaveBkgFile",".").toString();
    fn = QFileDialog::getSaveFileName( this, "Save loops for background calculation", fn, "Runfiles (*.sas_airun)",
                                       nullptr, QFileDialog::DontConfirmOverwrite );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".sas_airun") ) fn += ".sas_airun";
    data.setValue("LastSaveBkgFile",fn);
    data.endGroup();

    bool doOverwrite = true;
    if ( QFileInfo(fn).exists() )
    {
        doOverwrite = QMessageBox::question( this, "File exists", "The file\n"+fn+"\nexists.", "Overwrite", "Append" )
                      == 0;
        // Overwrite = 0
        // Append    = 1
    }
    performSaveAIOperation( fn, doOverwrite, true, false );
}

bool SC_MainGUI::performSaveAIOperation( QString fn, bool doOverwrite, bool interactive, bool useTPV )
{
    QFile f(fn);
    if ( !f.open( doOverwrite ? QIODevice::WriteOnly : QIODevice::ReadWrite ) )
    {
        qDebug() << fn << f.errorString();
        return false;
    }
    if ( !doOverwrite ) f.seek( f.size() );

    _loopDefinition rv;

    if ( useTPV )
    {
        rv = getTPVLoopDefinition();
        f.write( qPrintable("PATH|"+ui->inpTPVoutPath->text()+EOL) );
        f.write( "COLTBL|*DAT*" EOL ); // Spezielle Kennung für das Datenformat
        // Weitere Spezial-Parameter zur Bild-Nachbearbeitung
        QString tpv = "TPV|xx"; // Damit immer 2 Einträge da sind ...
        if ( ui->togTPVaddBS->isChecked() )
            tpv += QString("|BS%1;%2;%3;%4").arg(ui->inpTPVbeamx0->value()).arg(ui->inpTPVbeamy0->value())
                       .arg(ui->inpTPVnbeam->value()).arg(ui->inpTPVmbeam->value());
        if ( ui->togTPVaddLines->isChecked() )
            tpv += QString("|LI%1;%2").arg(ui->inpTPVaddLinesH->value()).arg(ui->inpTPVaddLinesV->value());
        if ( ui->togTPVaddNoise->isChecked() )
            tpv += "|NO";
        if ( ui->togTPVconvolute->isChecked() )
            tpv += "|CO";
        if ( ui->togTPVcalcRphi->isChecked() )
            tpv += "|RP";
        f.write( qPrintable(tpv+EOL) );
        // FFT wird weiter unten gesetzt
    }
    else
    {
        rv = getLoopDefinition();
        f.write( qPrintable("PATH|"+ui->inpSubDir->text()+EOL) );
        switch ( ui->cbsAIoutputFormat->currentIndex() )
        {
        case 0:
            // Kein Eintrag, es wird die Default-Farbtabelle genutzt (oder die per Argument)
            break;
        case 1:
            f.write( qPrintable("COLTBL|"+ui->cbsDefaultColTbl->currentText()+EOL) );
            break;
        case 2:
            f.write( "COLTBL|*DAT*" EOL ); // Spezielle Kennung für das Datenformat
            break;
        }
    }

    if ( useTPV ? ui->togTPVcalcFFT->isChecked() : ui->togAIUseFFTOutput->isChecked() )
    {   // FFT-Postprocessing Options
        QString fftOutFormat = "";
        if (      ui->radFFToutReal->isChecked()     ) fftOutFormat = "OutRe";
        else if ( ui->radFFToutImag->isChecked()     ) fftOutFormat = "OutIm";
        else if ( ui->radFFToutBetrag->isChecked()   ) fftOutFormat = "OutAbs";
        else if ( ui->radFFToutSpectrum->isChecked() ) fftOutFormat = "OutSpec";
        QString rphiScale = "";
        if ( ui->togFFTscaleRphi->isChecked()  ) rphiScale += "Scale ";
        if ( ui->togFFTclipRphi->isChecked()   ) rphiScale += "Clip1 ";
        if ( ui->togFFTclip40Rphi->isChecked() ) rphiScale += "Clip4 ";
        QString fftScale = "";
        if ( ui->togFFTscaleOut->isChecked()  ) fftScale += "Scale ";
        if ( ui->togFFTclipOut->isChecked()   ) fftScale += "Clip1 ";
        if ( ui->togFFTclip40Out->isChecked() ) fftScale += "Clip4 ";
        f.write( qPrintable(QString("FFT|%1;%2;%3;%4;%5;%6;%7")
                               .arg(ui->radFFTLinInput->isChecked()?"InLin":"InLog")
                               .arg(ui->cbsFFTsizeRphi->currentText().left(3).trimmed())
                               .arg(rphiScale.trimmed())
                               .arg(ui->cbsFFTsizeOut->currentText().left(3).trimmed())
                               .arg(fftOutFormat)
                               .arg(fftScale.trimmed())
                               .arg(ui->togIFFTSwap->isChecked()?"OutSwap":"OutNoSwap")
                           +EOL) );
    }
    ui->progressBar->setEnabled(true);
    int progCount = 0;
    foreach ( QString k, rv.uniqueKeys() )
    {
        QMultiMap<QString, _loopVariables>::iterator i = rv.find(k);
        while ( i != rv.end() && i.key() == k )
        {
            updateProgBar( ++progCount * 100 / rv.size() );
            QString line = k + "|";
            QMap<QString,double>::iterator i2 = i.value().begin();
            while ( i2 != i.value().end() )
            {
                line += QString("%1=%2|").arg(i2.key()).arg(i2.value());
                ++i2;
            }
            f.write( qPrintable(line+EOL) );
            ++i;
        }
    }
    ui->progressBar->setValue(0);
    ui->progressBar->setEnabled(false);
    f.close();
    if ( interactive )
        QMessageBox::information( this, "Loop definitions",
                                  "Classes: "+rv.uniqueKeys().join(", ")+"\n==> "+fn,
                                  QMessageBox::Ok );
    return true;
}


/**
 * @brief SC_MainGUI::on_butStart_clicked
 * Used to start the Cons-Prog in the background.
 */
void SC_MainGUI::on_butAIstart_clicked()
{
    performStartAI_TPV( false );
}

void SC_MainGUI::performStartAI_TPV( bool useTPV )
{
    globFlagTPV = useTPV; // Damit die Logausgabe in die richtige Liste geht
    // Test, ob das Ausgabeverzeichnis schreibbar ist
    QDir dout( useTPV ? ui->inpTPVoutPath->text() : ui->inpSubDir->text() );
    if ( ! dout.exists() )
    {
        if ( ! dout.mkpath(dout.absolutePath()) )
        {
            QMessageBox::critical( this, "Path not found", "The output path\n"+dout.absolutePath()+"\nis not available.", QMessageBox::Ok );
            return;
        }
    }

    QString cons = "sas_scatter2Cons";
#ifdef Q_OS_WIN
    // Unter Windows muss hier das passende Exe-File gesucht werden. Unter Linux muss
    // im User-Bin Verzeichnis (im PATH) ein passender Link vorhanden sein.
    QDir d( QDir::currentPath());
    while ( d.absolutePath().contains("build-") ) d.cdUp();
    QFileInfoList fil = d.entryInfoList( QStringList()<<"*"+cons+"*" );
    //qDebug() << fil;
    // (QFileInfo(C:\SimLab\sas-crystal\build-sas_scatter2Cons-Desktop_Static-Release))
    if ( fil.size() == 0 )
    {
        QMessageBox::critical( this, "Executable", "The executable "+cons+"\nis not found.", QMessageBox::Ok );
        return;
    }
    if ( fil.size() > 1 )
    {
        QMessageBox::critical( this, "Executable", "The executable "+cons+"\nis found in more than one versin.", QMessageBox::Ok );
        return;
    }
    QString tmp = fil[0].absoluteFilePath();
    if ( tmp.endsWith("Release") )
        tmp += "/release/" + cons + ".exe";
    else if ( tmp.endsWith("Debug") )
        tmp += "/debug/" + cons + ".exe";
    else
        tmp += "/" + cons + ".exe";
    qDebug() << tmp;
    if ( ! QFileInfo::exists(tmp) )
    {
        QMessageBox::critical( this, "Executable", "The executable "+tmp+"\nis not found.", QMessageBox::Ok );
        return;
    }
    cons = tmp;
#endif
    //--> Parameter- und Logfiles in temporärem Verzeichnis <--
    //QString fnrun = QDir::tempPath()+"/tmpForAI.sas_airun";
    //QString fnpar = QDir::tempPath()+"/tmpForAI.ini";
    //QString fnlog = QDir::tempPath()+"/AIgeneration.log";
    //--> Parameter- und Logfiles in Datenausgabe Verzeichnis <--
    QString fnrun = dout.absoluteFilePath("tmpForAI.sas_airun");
    QString fnpar = dout.absoluteFilePath("tmpForAI.ini");
    QString fnlog = dout.absoluteFilePath("AIgeneration.log");
    // --> Für jeden Typ nur eines von beiden freigeben !!!

    if ( !performSaveAIOperation( fnrun, true, false, useTPV ) ) return;
    performSaveParamOperation( fnpar );

    if ( aiBackProg == nullptr )
    {
        aiBackProg = new QProcess;
        connect( aiBackProg, SIGNAL(errorOccurred(QProcess::ProcessError)),
                this, SLOT(aiBackProg_error(QProcess::ProcessError)) );
        connect( aiBackProg, SIGNAL(finished(int,QProcess::ExitStatus)),
                this, SLOT(aiBackProg_finished(int,QProcess::ExitStatus)) );
        connect( aiBackProg, SIGNAL(readyReadStandardError()),
                this, SLOT(aiBackProg_readyRead()) );
        connect( aiBackProg, SIGNAL(readyReadStandardOutput()),
                this, SLOT(aiBackProg_readyRead()) );
        /*
        void errorOccurred(QProcess::ProcessError error)
        void finished(int exitCode, QProcess::ExitStatus exitStatus)
        void readyReadStandardError()
        void readyReadStandardOutput()
        void started()
        void stateChanged(QProcess::ProcessState newState)
        */
    }
    /*
    if ( QProcess::startDetached( cons, QStringList()<<"-p"<<fnpar<<"-c"<<fnrun ) )
        QMessageBox::information( this, "Loop definitions",
                                  "Background process to generate the images is started ...",
                                  QMessageBox::Ok );
    else
        QMessageBox::critical( this, "Loop definitions",
                               "Background process executable not found.\n"+cons,
                               QMessageBox::Ok );
    */
    if ( useTPV )
    {
        ui->lblTPVbackProc->setEnabled(true);
        ui->lisTPVbackProgOut->setEnabled(true);
        ui->lisTPVbackProgOut->clear();
    }
    else
    {
        ui->lblAIbackProc->setEnabled(true);
        ui->lisAIbackProgOut->setEnabled(true);
        ui->lisAIbackProgOut->clear();
    }
    if ( ui->togRemoveFiles->isChecked() || useTPV )
    {
        aiBackProgAddLog( "Remove files in "+dout.absolutePath() );
        localRemoveDirectory( dout.absolutePath() );
    }

    QStringList params;
    params<<"-p"<<fnpar<<"-c"<<fnrun<<"-l"<<fnlog<<"-t"<<QString::number(ui->inpNumCores->value());
    qDebug() << cons << params.join(" ");

    aiBackProg->start( cons, params );
}

bool SC_MainGUI::localRemoveDirectory( QString fn )
{
    qDebug() << "Scan" << fn;
    bool mainpath = false;
    QFileInfoList sl = QDir(fn).entryInfoList( QDir::Dirs | QDir::Files | QDir::NoDotAndDotDot );
    foreach ( QFileInfo f, sl )
    {
        if ( f.isDir() )
        {
            if ( ! localRemoveDirectory(f.absoluteFilePath()) )
                return false;
        }
        else if ( !f.absoluteFilePath().contains("tmpForAI") )
        {
            //qDebug() << "Remove" << f.absoluteFilePath();
            if ( ! QFile::remove(f.absoluteFilePath()) )
            {
                qDebug() << "Remove fail" << f.absoluteFilePath();
                return false;
            }
        }
        else
            mainpath = true;
    }
    if ( !mainpath )
    {
        QDir d(fn);
        QString dn = d.dirName();
        d.cdUp();
        qDebug() << "remove" << dn << d.absolutePath();
        if ( ! d.rmdir(dn) )
        {
            qDebug() << "Remove dir fail" << d.absolutePath() << dn;
            return false;
        }
    }
    qDebug() << "Complete" << fn;
    return true;
}

void SC_MainGUI::aiBackProg_error( QProcess::ProcessError err )
{
    switch ( err )
    {
    case QProcess::FailedToStart:   // The process failed to start. Either the invoked program is missing, or you may have insufficient permissions to invoke the program.
        QMessageBox::critical( this, "Loop definitions",
                              "Background process executable not found.\n"+aiBackProg->program(),
                              QMessageBox::Ok );
        break;
    case QProcess::Crashed:         // The process crashed some time after starting successfully.
    case QProcess::Timedout:        // The last waitFor...() function timed out. The state of QProcess is unchanged, and you can try calling waitFor...() again.
        break;
    case QProcess::WriteError:      // An error occurred when attempting to write to the process. For example, the process may not be running, or it may have closed its input channel.
    case QProcess::ReadError:       // An error occurred when attempting to read from the process. For example, the process may not be running.
        QMessageBox::critical( this, "Loop definitions",
                              "Background process unable to write to output device.",
                              QMessageBox::Ok );
        break;
    case QProcess::UnknownError:    // An unknown error occurred. This is the default return value of error().
        break;
    }
}
void SC_MainGUI::aiBackProg_finished( int /*code*/, QProcess::ExitStatus sta )
{
    if ( sta == QProcess::CrashExit )
        aiBackProgAddLog("Background process crashed.");
}
void SC_MainGUI::aiBackProg_readyRead()
{
    QString str;
    str = aiBackProg->readAllStandardOutput();
    aiBackProgAddLog( str.trimmed() );
    str = aiBackProg->readAllStandardError();
    aiBackProgAddLog( str.trimmed() );
}
void SC_MainGUI::aiBackProgAddLog( QString msg )
{
    if ( msg.isEmpty() ) return;
    if ( globFlagTPV )
    {
        ui->lisTPVbackProgOut->addItem(msg);
        ui->lisTPVbackProgOut->scrollToBottom();
        while ( ui->lisTPVbackProgOut->count() > 500 ) ui->lisTPVbackProgOut->takeItem(0);
    }
    else
    {
        ui->lisAIbackProgOut->addItem(msg);
        ui->lisAIbackProgOut->scrollToBottom();
        while ( ui->lisAIbackProgOut->count() > 500 ) ui->lisAIbackProgOut->takeItem(0);
    }
}



/**
 * @brief SC_MainGUI::on_butSaveImagesAI_clicked
 */
void SC_MainGUI::on_butSaveImagesAI_clicked()
{
    int saveGrayscale = ui->cbsAIoutputFormat->currentIndex();

    QString path = ui->inpSubDir->text();
    if ( ! path.endsWith("/") ) path += "/";

    QString imgdir = path;
    // Im TensorFlow-Script werden die Trainingsdaten aus dem Ordner "xxx" gelesen, dabei sind die
    // Unterverzeichnisse die Klassen. Die Messdaten werden aus "xxxfiles" gelesen.
    if ( ! imgdir.endsWith("files/") )
    {
        imgdir = imgdir.left(imgdir.length()-1) + "files/";
        QDir tmp(imgdir);
        if ( ! tmp.exists() ) tmp.mkpath(imgdir);
    }

    // For every open image window: reduce the size of the image and save the imagefile as *.png
    for ( int i=0; i<images.size(); i++ )
    {
        QString fn = QString("File%1").arg(i,3,10,QChar('0'));
        switch ( saveGrayscale )
        {
        case 0:
            images.at(i)->saveImageGray( imgdir+fn+".png", QSize() );
            break;
        case 1:
            images.at(i)->saveImage( imgdir+fn+".png", QSize() );
            break;
        case 2:
            // TODO
            break;
        }
        qDebug() << "Save" << i << imgdir+fn+".png";
    }
}

//typedef QHash<QString/*Variable*/, double/*value*/> _loopVariable;
//typedef QMultiHash<QString/*CLASS*/, _loopVariables > _loopDefinition;
SC_MainGUI::_loopDefinition SC_MainGUI::getLoopDefinition()
{
    _loopDefinition retval;  // MultiHash für alle Rechenschritte

    // Zuerst einige globale Informationen

    QFileInfoList files;
    if ( ui->grpFileInput->isChecked() )
    {
        files = QDir(ui->inpFileName->text()).entryInfoList( QDir::Files | QDir::NoDotAndDotDot );
        qDebug() << files.size() << "Files found in" << ui->inpFileName->text();
    }

    QStringList slKeys; // Schlüssel für den Daten-Hash, damit immer gleiche Reihenfolge genutzt wird
    typedef struct { double start, end, step, current, var; int randcount; QString formel; } _werte;
    QHash<QString,_werte*> daten;
    for ( int r=0; r<ui->tblListe->rowCount(); r++ )
        if ( ! ui->tblListe->item(r,3)->text().isEmpty() /*&&
             master->metaData.contains(ui->tblListe->item(r,0)->text())*/ )
        {   // Werte aus anderen Berechnungsformeln und leere Schrittweiten werden nicht weiter verfolgt
            _werte *w = new _werte;
            w->start  = ui->tblListe->item(r,1)->text().toDouble();
            if ( ui->tblListe->item(r,2)->text().startsWith("V") )
            {   // V<num> Variation
                w->var = ui->tblListe->item(r,2)->text().mid(1).toDouble();
                w->end = w->start;
            }
            else
            {   // Normaler Endwert
                w->end = ui->tblListe->item(r,2)->text().toDouble();
                w->var = 0;
            }
            if ( ui->tblListe->item(r,3)->text().startsWith("R") )
            {   // R<num> ist Zufallszahlen
                w->randcount = ui->tblListe->item(r,3)->text().mid(1).toInt();
                w->step = 0;
                w->formel = "";
            }
            else if ( ui->tblListe->item(r,3)->text().startsWith("=") )
            {   // = ... calculation ...
                w->randcount = 0;
                w->step = 0;
                w->var = 0;
                w->formel = ui->tblListe->item(r,3)->text().mid(1).trimmed();
            }
            else
            {   // <num> ist Schrittweite
                w->step = ui->tblListe->item(r,3)->text().toDouble();
                w->randcount = 0;
                w->formel = "";
                if ( w->var > 0 )
                {
                    QMessageBox::critical( this, "AI generation error",
                                          QString("In row %1: numerical stepsize and variation given. This row is ignored.").arg(r+1),
                                          QMessageBox::Ok );
                    continue;
                }
            }
            w->current = w->start;
            daten.insert( ui->tblListe->item(r,0)->text(), w );
            slKeys.append( ui->tblListe->item(r,0)->text() );
        }

    QString curCLASS;

    // Äußerste Schleife über die Methoden
    QStringList slCalcArt = calcGui->getCalcTypes();
    for ( int iCurMethod=0; iCurMethod<slCalcArt.size(); iCurMethod++ )
    {
        if ( ui->togUseAllMethods->isChecked() )
        {   // all methods are used, so set the metaData structure
            //master->getMetaDataFromTbl(iCurMethod);
            //QString SC_CalcGUI::currentParamValue( QString m, QString p, bool text )
        }
        else if ( ui->cbsMethod->currentText() != slCalcArt[iCurMethod] )
        {   // Skip this method
            continue;
        }

        calcval.clear();

        // Neue Methode, Tabellenwerte auf Start setzen
        foreach ( QString k, slKeys )
        {
            daten[k]->current = daten[k]->start;
            if ( daten[k]->randcount > 0 )
                daten[k]->step = daten[k]->randcount;
        }

        int usedKey = 0; // Verwendeter Schlüssel, falls keine Kaskadierten Schleifen genutzt werden sollen

        // Nächste Schachtelungstiefe: die Dateien
        int iFiles = -1;
        while ( true )
        {
            if ( files.size() > 0 )
            {
                iFiles++;
                if ( iFiles >= files.size() )
                    break;
                // Aktuelle Datei einlesen und vorbereiten
                QFile finp(files.at(iFiles).absoluteFilePath());
                if ( !finp.fileName().endsWith(".inp") ) continue;
                if ( !finp.open(QIODevice::ReadOnly) )
                {
                    qDebug() << finp.fileName() << finp.errorString();
                    break;
                }
                // The first line contains the keys for the following values
                // First char '#' has to be skipped
                QStringList keys = QString(finp.readLine()).mid(1).trimmed().split("\t");
                // All other lines contains the same number of values
                QString lines = finp.readAll();
                lines.replace("\r\n","\t");
                lines.replace("\n","\t");
                QStringList values = lines.trimmed().split("\t");
                finp.close();
                for ( int k=0; k<keys.size(); k++ )
                {
                    if ( keys.at(k) == "CLASS" )
                    {   // Normalerweise ein String
                        curCLASS = values.at(k);
                    }
                    else
                    {   // Normalerweise Zahlen
                        bool ok;
                        double v = values.at(k).toDouble(&ok);
                        if ( !ok )
                            qDebug() << k << "invalid double" << values.at(k);
                        else if ( calcGui->isCurrentParameterValid( slCalcArt[iCurMethod], keys.at(k), false ) )
                            calcval.insert( keys.at(k), v );
                        else if ( calcGui->isCurrentParameterValid( slCalcArt[iCurMethod], "Edit"+keys.at(k), false ) )
                            calcval.insert( "Edit"+keys.at(k), v );
                        else
                            qDebug() << k << "not found" << keys.at(k) << values.at(k);
                    }
                }
                // Bei jedem File die Tabellenwerte auf Start setzen
                foreach ( QString k, slKeys )
                {
                    daten[k]->current = daten[k]->start;
                    if ( daten[k]->randcount > 0 )
                        daten[k]->step = daten[k]->randcount;
                }
            } // Read current file

            while ( true )
            {
                // Innerste Schleife über die Tabelle
                foreach ( QString k, slKeys )
                    if ( daten[k]->formel.isEmpty() && calcGui->isCurrentParameterValid( slCalcArt[iCurMethod], k, false ) )
                        calcval.insert( k, daten[k]->current );
                foreach ( QString k, slKeys )
                    if ( ! daten[k]->formel.isEmpty() && calcGui->isCurrentParameterValid( slCalcArt[iCurMethod], k, false ) )
                        calcval.insert( k, evaluateFormula( slCalcArt[iCurMethod], daten[k]->formel ) );
                calcval.insert( "#Calculation#", iCurMethod );
/*
                // Abfrage, ob ein File mit "default" im Namen gelesen wurde
                if ( files.size() > 0 )
                {
                    if ( files.at(iFiles).fileName().contains("default",Qt::CaseInsensitive) )
                    {
                        calcval.insert( "#DEF#", 0.0 );
                        // Wenn ja, dann wird eine spezielle Kennung in die Datei geschrieben,
                        // wodurch der Filename den Zusatz "_DEF" bekommt. Ist dann leichter
                        // zu finden.
                    }
                    else
                    {
                        calcval.remove("#DEF#");
                        // Remove the key from the list. We use no clear() to erase the whole list,
                        // the insert() do a replace if the key is available. This is for a better
                        // memory management.
                    }
                }
*/
                // Berechnung der CLASS
                QString cls;
                cls = ui->inpFileClass->text();
                cls.replace( "{M}", slCalcArt[iCurMethod] ); // Historisch...
                int pos;
                while ( (pos = cls.indexOf("{")) >= 0 )
                {
                    int p2 = cls.indexOf("}");
                    QString k = cls.mid(pos+1,p2-pos-1);
                    QString korg = k;
                    double div = 1.0;
                    int p = k.indexOf("/");
                    if ( p > 0 )
                    {
                        div = k.mid(p+1).toDouble();
                        k   = k.left(p);
                        if ( div == 0 ) div = 1;
                    }
                    if ( calcval.contains(k) )
                    {
                        if ( div != 1 )
                        {
                            int tmp = static_cast<int>(calcval[k] / div);
                            cls.replace( "{"+korg+"}", QString::number(tmp*div) );
                        }
                        else
                            cls.replace( "{"+korg+"}", QString::number(calcval[k]) );
                    }
                    else if ( calcGui->isCurrentParameterValid( slCalcArt[iCurMethod], k, false ) )
                    {
                        if ( div != 1 )
                        {
                            int tmp = calcGui->currentParamValueInt(slCalcArt[iCurMethod],k) / div;
                            cls.replace( "{"+korg+"}", QString::number(tmp*div) );
                        }
                        else
                            cls.replace( "{"+korg+"}", calcGui->currentParamValueStr(slCalcArt[iCurMethod],k,false) );
                    }
                    else
                        cls.replace( "{"+korg+"}", "unknown" );
                }

                curCLASS = cls;

                if ( ui->radLoopsSingle->isChecked() )
                {
                    QString k = slKeys[usedKey];
                    calcval.insert( "#FN_"+k, daten[k]->current );

                    if (  daten[k]->randcount == 0 &&       // Keine Zufallszahlen
                         !daten[k]->formel.isEmpty() )      // und Formel angegeben
                        curCLASS = "";                      // --> jetzt nicht speichern
                }

                // Hier sind alle Daten für einen Rechenschritt fertig ...
                if ( !curCLASS.isEmpty() )
                    retval.insert( curCLASS, calcval );

                if ( slKeys.size() > 0 )
                {
                    if ( ui->radLoopsCascade->isChecked() )
                    {
                        // Einen Wert in der Tabelle erhöhen und Test, ob alle fertig sind
                        int fertig = 0;
                        foreach (QString k, slKeys)
                        {
                            if ( daten[k]->randcount > 0 )
                            {   // Zufallswerte, count gibt die Anzahl vor
                                if ( daten[k]->var > 0 )
                                {   // Spezialfall: Variation angegeben
                                    // PAS: wert * ( 1 + varianz * ( random(200)-100)/100 )
                                    double zuf = ((static_cast<double>(rand())/RAND_MAX * 200.0) - 100.0) / 100.0; // -1 .. 1
                                    daten[k]->current = daten[k]->start * ( 1.0 + daten[k]->var * zuf );
                                }
                                else
                                {   // Normalfall: Start / Ende angegeben
                                    daten[k]->current = static_cast<double>(rand())/RAND_MAX
                                                            * ( daten[k]->end - daten[k]->start )
                                                        + daten[k]->start;
                                }
                                if ( --daten[k]->step > 0 ) break;
                                fertig++;
                                daten[k]->step = daten[k]->randcount;
                            }
                            else if ( daten[k]->formel.isEmpty() )
                            {   // Normale Schrittweite
                                daten[k]->current += daten[k]->step;
                                if ( daten[k]->step < 0 )
                                {   // Negative Schrittweite beachten
                                    if ( daten[k]->current >= daten[k]->end ) break;
                                    fertig++;
                                    daten[k]->current = daten[k]->start; // Startwert setzen
                                }
                                else if ( daten[k]->step > 0 )
                                {
                                    if ( daten[k]->current <= daten[k]->end ) break;
                                    fertig++;
                                    daten[k]->current = daten[k]->start; // Startwert setzen
                                }
                                else
                                    fertig++; // TODO: oder Fehlermeldung schon weiter oben?
                            }
                            else
                            {   // Die Formel wird oben neu berechnet... ist aber hier immer fertig
                                fertig++;
                            }
                        }
                        if ( fertig == daten.size() ) break;
                    } // if ( ui->radLoopsCascade->isChecked() )
                    else if ( ui->radLoopsSingle->isChecked() )
                    {
                        QString k = slKeys[usedKey];
                        bool fertig = false;
                        if ( daten[k]->randcount > 0 )
                        {   // Zufallswerte, count gibt die Anzahl vor
                            daten[k]->current = static_cast<double>(rand())/RAND_MAX
                                                    * ( daten[k]->end - daten[k]->start )
                                                + daten[k]->start;
                            //qDebug() << "ran" << k << daten[k]->current << daten[k]->step;
                            if ( --daten[k]->step == 0 )
                            {
                                fertig = true;
                                daten[k]->step = daten[k]->randcount;
                            }
                        }
                        else if ( daten[k]->formel.isEmpty() )
                        {   // Normale Schrittweite
                            daten[k]->current += daten[k]->step;
                            //qDebug() << "stp" << k << daten[k]->current << daten[k]->step;
                            if ( daten[k]->step < 0 )
                            {   // Negative Schrittweite beachten
                                if ( daten[k]->current < daten[k]->end )
                                {
                                    fertig = true;
                                    daten[k]->current = daten[k]->start; // Startwert setzen
                                }
                            }
                            else
                            {
                                if ( daten[k]->current > daten[k]->end )
                                {
                                    fertig = true;
                                    daten[k]->current = daten[k]->start; // Startwert setzen
                                }
                            }
                        }
                        else
                        {   // Die Formel wird oben neu berechnet... ist aber hier immer fertig
                            //qDebug() << "formel" << k << daten[k]->current << daten[k]->formel;
                            fertig = true;
                        }
                        if ( fertig )
                        {   // Jetzt sind alle Schritte der Zeile abgearbeitet
                            calcval.remove( "#FN_"+slKeys[usedKey] );
                            if ( ++usedKey >= slKeys.size() ) break;
                        }
                    } // if ( ui->radLoopsSingle->isChecked() )
                    else
                    {
                    } // if ( ui->radLoopsFixed->isChecked() )
                } // if slKeys.size() > 0
                else
                    break; // Keine Tabellen-Keys, keine Schleife hier
            } // while (true) - Schleife über die Tabelle
            if ( files.size() == 0 ) break; // Keine Files, also hier fertig.
        } // while (true) - Schleife über die Files
    } // foreach ( QString curMethod, slCalcArt ) - Schleife über die Methoden
    return retval;
}

double SC_MainGUI::evaluateFormula( QString m, QString formel )
{
    QStringList sl = formel.split(QRegExp("[+-*/ ]"),Qt::SkipEmptyParts);
    //qDebug() << "CALC" << formel << sl;
    double val = 0;
    if ( calcval.contains(sl[0]) )
        val = calcval[ sl[0] ];
    else if ( calcGui->isCurrentParameterValid(m,sl[0],false) )
        val = calcGui->currentParamValueDbl(m,sl[0]);
    if ( sl.size() < 2 ) return val;
    double val2 = sl[1].toDouble();
    formel = formel.mid(sl[0].length()).trimmed();
    //qDebug() << "CALC-" << formel << val << val2;
    switch ( formel[0].toLatin1() )
    {
    case '+': return val + val2;
    case '-': return val - val2;
    case '*': return val * val2;
    case '/': if ( val2 != 0 ) return val / val2;
        else return 0;
    }
    return 0;
}


/**
 * @brief SC_MainGUI::on_butTestGo_clicked
 * Generate Test Images ...
 */
void SC_MainGUI::on_butTestGo_clicked()
{
    prepareCalculation( true, false );

    // TODO: LogTimer starten (TestGo)

    int mx = ui->inpGridPoints->value();
    double *data = new double[4*mx*mx]; // mx is one quarter of result image

    QString title = "??";
    SasCalc_PostProc::inst()->setLogging( false );

    if ( ui->radTestImg1->isChecked() )
    {   // Test Image 1 (-)
        title = "Test-Image (-)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( x >= mx-2 && x <= mx+2 ) // Linie entlang y in der Mitte
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
    }
    else if ( ui->radTestImg2->isChecked() )
    {   // Test Image 1 (|)
        title = "Test-Image (|)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( y >= mx-2 && y <= mx+2 ) // Linie entlang x in der Mitte
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
    }
    else if ( ui->radTestImg3->isChecked() )
    {   // Test Image 1 (x)
        title = "Test-Image (x)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( abs(x-y) <= 2 || abs(x-(2*mx-y-1)) <= 2 )
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
    }
    else if ( ui->radTestImg4->isChecked() )
    {   // Test Image 4 (+)
        title = "Test-Image (+)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( (x >= mx-2 && x <= mx+2) || (y >= mx-2 && y <= mx+2) )
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
    }
    else if ( ui->radTestImg5->isChecked() )
    {   // Test Image 5 (:.)
        title = "Test-Image (:.)";
        double fac = 2.0*M_PI / (mx*2);
        for ( int y=0, idx=0; y<2*mx; y++ )
            for ( int x=0; x<2*mx; x++, idx++ )
                if ( x < mx )
                    data[idx] = (1.0 + sin( x * fac ) * sin( y * fac )) * 128.0;  // -> 0..255
                else if ( y > mx )
                    data[idx] = (1.0 + sin( x * fac ) * sin( y * fac )) * 64.0;  // -> 0..128
                else
                    data[idx] = 128;    // Quadrant rechts oben (ohne Drehung) konstant
    }
    else if ( ui->radTestImg6->isChecked() )
    {   // Test Image 6 (.)
        title = "Test-Image (.)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( x >= mx-2 && x <= mx+2 && y >= mx-2 && y <= mx+2 ) // Mittelpunkt
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
    }
    else if ( ui->radTestCurMethod->isChecked() )
    {   // Current method
        title = "Test-" + calcGui->curMethod->subCalc->methodName();
        calcGui->doCalculation( ui->inpNumCores->value(), static_cast<progressAndAbort>(&this->myProgressAndAbort) );
        memcpy( data, calcGui->data(), 4*mx*mx*sizeof(double) );
    }
    else if ( ui->radTestAllMethods->isChecked() )
    {   // All methods
        QStringList slMethods = calcGui->getCalcTypes();
        foreach( QString m, slMethods )
        {
            calcGui->prepareCalculation( m, true );
            calcGui->doCalculation( ui->inpNumCores->value(), static_cast<progressAndAbort>(&this->myProgressAndAbort) );
            title = "Test-" + calcGui->curMethod->subCalc->methodName();
            widImage* img = addImage( true, 0, 2*mx, 0, 2*mx, calcGui->data(), title, false );
            img->addMetaInfo( "From Image", title );
            // FFT
            if ( ui->togTestUseFFT->isChecked() )
            {
                performIFFTcalculations( images.size()-1, title, false );  // use the last image
            }
        }
        finishCalculation(false);
        return;
    }
    else
    {
        finishCalculation(false);
        return;
    }

    widImage* img = addImage( true, 0, 2*mx, 0, 2*mx, data, title, false );
    img->addMetaInfo( "From Image", title );

    // FFT
    if ( ui->togTestUseFFT->isChecked() )
    {
        SasCalc_PostProc::inst()->setLogging( true );  // besser für Debug
        performIFFTcalculations( images.size()-1, title, false );  // use the last image
    }

    finishCalculation(false);
}

void SC_MainGUI::on_cbsDefaultColTbl_activated(int index)
{
    if ( index < 0 ) index = 1;
    // Es kann bei falschem Settings-Eintrag mal die -1 hier vorkommen.
    // Die macht dann später einen Crash, weil sonst nicht auf eine ungültige
    // Farbtabelle abgefragt wird.
    widImage::setDefColTbl(index);
}



void SC_MainGUI::on_actionTest_read_first_AI_RUN_Line_triggered()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("AI");
    QString fn = data.value("LastSaveBkgFile",".").toString();

    QFile fcalc(fn);
    if ( !fcalc.open(QIODevice::ReadOnly) )
    {
        qDebug() << fn + " " + fcalc.errorString();
        return;
    }
    QString m = "";
    while ( !fcalc.atEnd() )
    {
        QStringList sl = QString(fcalc.readLine()).trimmed().split("|",Qt::SkipEmptyParts);
        if ( sl.size() < 3 ) continue;
        qDebug() << sl;

        for ( int i=1; i<sl.size(); i++ )
        {
            QStringList vv = sl[i].split("=");
            if ( sl[i].startsWith("#Calc") )
            {
                m = calcGui->getCalcTypes()[vv[1].toInt()];
            }
            else if ( ! sl[i].startsWith("#DEF") )
            {
                if ( calcGui->updateParamValue( m, vv[0], vv[1].toDouble(), Qt::black ) ) continue;
                //bool SC_CalcGUI::updateParamValue( QString m, QString p, double v )

                if ( vv[0] == "EditAxis1x" )
                    ui->inpAx1->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis1y" )
                    ui->inpAy1->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis1z" )
                    ui->inpAz1->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis2x" )
                    ui->inpAx2->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis2y" )
                    ui->inpAy2->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis2z" )
                    ui->inpAz2->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis3x" )
                    ui->inpAx3->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis3y" )
                    ui->inpAy3->setValue( vv[1].toDouble() );
                else if ( vv[0] == "EditAxis3z" )
                    ui->inpAz3->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editxrel" )
                    ui->inpN1->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Edityrel" )
                    ui->inpN2->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editzrel" )
                    ui->inpN3->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editdom1" )
                    ui->inpSigX->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editdom2" )
                    ui->inpSigY->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editdom3" )
                    ui->inpSigZ->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editx1rel" )
                    ui->inpU1->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Edity1rel" )
                    ui->inpU2->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editz1rel" )
                    ui->inpU3->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editx2rel" )
                    ui->inpV1->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Edity2rel" )
                    ui->inpV2->setValue( vv[1].toDouble() );
                else if ( vv[0] == "Editz2rel" )
                    ui->inpV3->setValue( vv[1].toDouble() );
                else
                    qDebug() << vv[0] << "Falscher Key???";
            }
        }

        break;
    }
}



void SC_MainGUI::loadScatterParameter( QString fn )
{
    QFile fin(fn);
    if ( ! fin.open(QIODevice::ReadOnly) )
    {
        qDebug() << fin.errorString();
        return;
    }
    QString all = fin.readAll();
    fin.close();
    QStringList data;
    if ( all.contains("\r\n") )
        data = all.split("\r\n");
    else
        data = all.split("\n");
    qDebug() << data.size();

    checkData( data[0], "Parameter file" );

    /* File Info */
    //outArray[1] = EditFilename.Text;  ==> meistens leer
    //outArray[2] = EditBackFilename.Text;
    checkData( data[3], "0" );          /* A --> nm */
    checkData( data[4], "0" );          /* s -_> q */
    checkData( data[5], "0" );          /* theta --> q */

    setValue( data[6], "EditI0" );
    //outArray[7] = EditBase.Text;
    //outArray[8] = EditSubFactor.Text;
    //outArray[9] = IntToStr(ord(CheckBoxSub.Checked));

    checkData( data[10], "1" );         /* binning */
    checkData( data[11], "1" );         /* data points */
    //outArray[12] = EditMaxInt.Text;

    /** Crystal Lattice Info **/
    //setValue( data[14], "ComboBoxLattice" ); // Berechnungsmethode
    setValue( data[15], "ComboBoxPeak" );
    setValue( data[16], "EditLattice" );
    setValue( data[17], "EditLatticeb" );
    setValue( data[18], "EditLatticec" );
    setValue( data[19], "EditDomainSize" );
    setValue( data[20], "EditAzi" );
    setValue( data[21], "EditDebyeWaller" );
    checkData( data[22], "0.05" );       /* ceff for PY */
    checkData( data[23], "10" );         /* Reff for PY */
    setValue( data[24], "EditTwratio" );
    //outArray[25] = EditAstack.Text;
    //outArray[26] = EditBstack.Text;
    setValue( data[27], "Edithklmax" );
    checkData( data[28], "0.1" );         /* qz for GISAXS */
    setValue( data[29], "EditPeakPar" );
    //outArray[30] = EditCeffCyl.Text;
    //outArray[31] = Editacrit.Text;

    /* Particle Info */
    setValue( data[32], "ComboBoxParticle" );
    setValue( data[33], "ComboBoxInterior" );
    setValue( data[34], "EditRadius" );
    setValue( data[35], "EditStdDev" );
    setValue( data[36], "EditLength" );
    setValue( data[37], "EditSigmal" );
    setValue( data[38], "EditRadiusi" );
    setValue( data[39], "EditAlpha" );
    setValue( data[40], "EditRho" );
    setValue( data[41], "EditDbeta" );

    /* Aniso-Gauss Info */
    setValue( data[42], "EditAxis1x" );
    setValue( data[43], "EditAxis1y" );
    setValue( data[44], "EditAxis1z" );
    setValue( data[45], "EditAxis2x" );
    setValue( data[46], "EditAxis2y" );
    setValue( data[47], "EditAxis2z" );
    setValue( data[48], "EditAxis3x" );
    setValue( data[49], "EditAxis3y" );
    setValue( data[50], "EditAxis3z" );
    setValue( data[51], "Editdom1" );
    setValue( data[52], "Editdom2" );
    setValue( data[53], "Editdom3" );

    /* Fit Info */
    checkData( data[54], "1" );       /* Fit: first */
    //outArray[57] = '1';       /* Fit: last */
    //outArray[58] = '10';      /* Fit: iterations */
    //outArray[59] = '0';       /* Fit: current interation */
    //outArray[60] = '0';       /* Simulated annealing */
    //outArray[61] = '1';       /* Levenberg-Marquardt */

    //outArray[62] = '0';       /* Fitted parameters ... */
    //outArray[63] = '0';
    //outArray[64] = '0';
    //outArray[65] = '0';
    //outArray[66] = '0';
    //outArray[67] = '0';
    //outArray[68] = '0';
    //outArray[69] = '0';
    //outArray[70] = '0';
    //outArray[71] = '0';
    //outArray[72] = '0';
    //outArray[73] = '0';
    //outArray[76] = '0';
    //outArray[77] = '0';
    //outArray[78] = '0';

    //outArray[79] = '0';
    //outArray[80] = '0';
    //outArray[81] = '0';
    //outArray[82] = '0';
    //outArray[83] = '0';
    //outArray[84] = '0';
    //outArray[85] = '0';
    //outArray[86] = '0';
    //outArray[87] = '0';
    //outArray[88] = '0';
    //outArray[89] = '0';

    //outArray[90] = '0';
    //outArray[91] = '0';
    //outArray[92] = '0';

    //outArray[93] = '10';
    //outArray[94] = '1e-15';
    //outArray[95] = '';
    //outArray[97] = '';

    /* Multiphase Info */
    //outArray[102] = EditPhase1.Text;
    //outArray[103] = EditPhase2.Text;
    //outArray[104] = EditPhase2.Text;
    //outArray[107] = EditMultBase.Text;
    //outArray[108] = IntToStr(ord(RadioButtonPhase1.Checked));
    //outArray[109] = IntToStr(ord(RadioButtonPhase2.Checked));
    //outArray[110] = IntToStr(ord(RadioButtonPhase3.Checked));

    /* Porod Info */
    setValue( data[111], "EditPhi" );
    //outArray[112] = EditArea.Text;
    //outArray[113] = Editlp.Text;
    checkData( data[114], "0.01" );            /* qmin */
    //outArray[115] = '1';               /* qmax */
    setValue( data[116], "EditBFactor" );
    //outArray[117] = '';                /* Porod factor */
    //outArray[118] = '';                /* Integral */
    //outArray[119] = '';                /* Porod length */
    //outArray[120] = '';                /* spec. area */

    /* Displacement Info */
    setValue( data[121], "EditDist" );
    //outArray[122] = '';                /* rel. near. neighbor dist. */
    setValue( data[123], "EditOrd" );      /* !! order parameter */
    /* outArray[124]:=EditCailleT.Text; */

    /* Plot Info */
    //outArray[125] = IntToStr(ord(RadioButtonXlog.Checked));
    //outArray[126] = IntToStr(ord(RadioButtonXlin.Checked));
    //outArray[127] = IntToStr(ord(RadioButtonZlog.Checked));
    //outArray[128] = IntToStr(ord(RadioButtonZlin.Checked));
    //outArray[129] = EditXmin.Text;
    //outArray[130] = EditXmax.Text;
    //outArray[131] = EditYmin.Text;
    //outArray[132] = EditYmax.Text;
    //outArray[133] = '1';     /* autoscale */
    //outArray[134] = '0';     /* overlay */
    //outArray[135] = '1';     /* solution */
    //outArray[136] = '0';     /* solvent */
    //outArray[137] = '0';     /* sample */
    //outArray[138] = '0';     /* fit */
    //outArray[139] = '0';     /* calculation */
    //outArray[140] = '0';     /* tau distribution */
    //outArray[141] = '0';     /* size distribution */
    //outArray[142] = '0';     /* multi phase */
    //outArray[143] = '0';     /* plot */
    //outArray[144] = IntToStr(ord(RadioButtonTickFew.Checked));
    //outArray[145] = IntToStr(ord(RadioButtonTickMany.Checked));
    //outArray[146] = '0';     /* counter */

    //outArray[147] = IntToStr(ord(RadioButtonTickNone.Checked));
    //outArray[148] = IntToStr(ord(RadioButtonZlog.Checked));
    //outArray[149] = IntToStr(ord(RadioButtonZlin.Checked));
    //outArray[150] = IntToStr(ComboBoxMap.ItemIndex);

    /* Calculation */
    //outArray[151] = EditHex2h.Text;
    //outArray[152] = EditHex2k.Text;
    //outArray[153] = EditHex2i.Text;
    //outArray[154] = EditHex2l.Text;
    setValue( data[155], "Editu1" );
    setValue( data[156], "Editu2" );
    setValue( data[157], "Editu3" );
    setValue( data[158], "Editv1" );
    setValue( data[159], "Editv2" );
    setValue( data[160], "Editv3" );
    setValue( data[161], "Rot_X" );
    setValue( data[162], "Rot_Y" );
    setValue( data[163], "Rot_Z" );
    //outArray[164] = EditDrehx.Text;
    //outArray[165] = EditDrehy.Text;
    //outArray[166] = EditDrehz.Text;
    setValue( data[167], "Rot_Angle" /*"EditRotAlpha"*/ );
    setValue( data[168], "EditDetector" );
    setValue( data[169], "EditGridPoints" );
    setValue( data[170], "EditQmax" );
    //outArray[171] = IntToStr(ord(RadioButtonFixU.Checked));
    //outArray[172] = IntToStr(ord(RadioButtonFixV.Checked));
    //outArray[173] = IntToStr(ord(RadioButtonFixNone.Checked));
    setValue( data[174], "RadioButtonQ1" );
    setValue( data[175], "RadioButtonQ2" );
    setValue( data[176], "RadioButtonQ4" );
    // TODO     ui->togExpandImage->setChecked( sets.value( "ExpandImage", false ).toBool() );
    //outArray[177] = IntToStr(ord(RadioButtonMidPoint.Checked));
    //outArray[178] = IntToStr(ord(RadioButtonBeam.Checked));
    //outArray[179] = IntToStr(ord(RadioButtonQData.Checked));
    //outArray[180] = IntToStr(ord(RadioButtonPreset.Checked));

    /* Analysis */
    //outArray[181] = EditX.Text;
    //outArray[182] = EditY.Text;
    //outArray[183] = EditZ.Text;
    setValue( data[184], "EditAnglexy" );
    //outArray[185] = EditIntMap.Text;
    //outArray[186] = EditIntTheo.Text;
    //outArray[187] = Edithout.Text;
    //outArray[188] = Editkout.Text;
    //outArray[189] = Editlout.Text;
    //outArray[190] = EditEwaldDist.Text;

    //outArray[191] = IntToStr(ord(CheckBoxLeftHoriz.Checked));
    //outArray[192] = IntToStr(ord(CheckBoxLeftVert.Checked));
    //outArray[193] = IntToStr(ord(CheckBoxLeftCircle.Checked));
    //outArray[194] = IntToStr(ord(CheckBoxRightHoriz.Checked));
    //outArray[195] = IntToStr(ord(CheckBoxRightVert.Checked));
    //outArray[196] = IntToStr(ord(CheckBoxRightCircle.Checked));
    setValue( data[197], "Tilt_Angle" /*"EditTilt"*/ );
    setValue( data[198], "Editx" );
    setValue( data[199], "Edity" );
    setValue( data[200], "EditAnglexy" );

    //outArray[201] = EditRedAngle.Text;
    //outArray[202] = EditGreenAngle.Text;

    /* Data info */
    //outArray[203] = EditTitle.Text;
    //outArray[204] = EditILLDate.Text;
    //outArray[205] = EditILLRun.Text;
    setValue( data[206], "EditWavelength" );
    setValue( data[207], "EditDet" );
    //outArray[208] = EditILLTemp.Text;
    setValue( data[209], "EditPixelX" );
    setValue( data[210], "EditPixelY" );
    setValue( data[211], "EditCenterX" );
    setValue( data[212], "EditCenterY" );
    setValue( data[213], "EditPixelNoX" );
    setValue( data[214], "EditPixelNoY" );
    setValue( data[215], "EditAlphai" );
    //outArray[216] = EditMinInt.Text;
    //outArray[217] = EditMaxInt.Text;
    setValue( data[218], "Editqamax" );
    setValue( data[219], "CheckBoxWAXS" );
    setValue( data[220], "EditWAXSangle" );
    //outArray[221] = EditWAXSvx.Text;
    //outArray[222] = EditWAXSvy.Text;
    //outArray[223] = IntToStr(ord(RadioButtonHoriz.Checked));
    //outArray[224] = IntToStr(ord(RadioButtonVertical.Checked));
}

void SC_MainGUI::checkData( QString data, QString soll )
{
    if ( data.trimmed().compare(soll,Qt::CaseInsensitive) != 0 )
        qDebug() << "checkData ungleich" << data << soll;
}

// Verwendet von loadScatterParameter()
void SC_MainGUI::setValue( QString data, QString name )
{
    data.replace(",",".");
    QString m = ui->tabMethods->tabText(ui->tabMethods->currentIndex());
    //qDebug() << "setScatterValue" << name << data;
    if ( calcGui->updateParamValue( m, name, data.toDouble(), Qt::black ) ) return;

    if ( name == "EditAxis1x" )
        ui->inpAx1->setValue( data.toDouble() );
    else if ( name == "EditAxis1y" )
        ui->inpAy1->setValue( data.toDouble() );
    else if ( name == "EditAxis1z" )
        ui->inpAz1->setValue( data.toDouble() );
    else if ( name == "EditAxis2x" )
        ui->inpAx2->setValue( data.toDouble() );
    else if ( name == "EditAxis2y" )
        ui->inpAy2->setValue( data.toDouble() );
    else if ( name == "EditAxis2z" )
        ui->inpAz2->setValue( data.toDouble() );
    else if ( name == "EditAxis3x" )
        ui->inpAx3->setValue( data.toDouble() );
    else if ( name == "EditAxis3y" )
        ui->inpAy3->setValue( data.toDouble() );
    else if ( name == "EditAxis3z" )
        ui->inpAz3->setValue( data.toDouble() );
    else if ( name == "Editxrel" )
        ui->inpN1->setValue( data.toDouble() );
    else if ( name == "Edityrel" )
        ui->inpN2->setValue( data.toDouble() );
    else if ( name == "Editzrel" )
        ui->inpN3->setValue( data.toDouble() );
    else if ( name == "Editdom1" )
        ui->inpSigX->setValue( data.toDouble() );
    else if ( name == "Editdom2" )
        ui->inpSigY->setValue( data.toDouble() );
    else if ( name == "Editdom3" )
        ui->inpSigZ->setValue( data.toDouble() );
    else if ( name == "Editx1rel" || name == "Editu1" )
        ui->inpU1->setValue( data.toDouble() );
    else if ( name == "Edity1rel" || name == "Editu2" )
        ui->inpU2->setValue( data.toDouble() );
    else if ( name == "Editz1rel" || name == "Editu3" )
        ui->inpU3->setValue( data.toDouble() );
    else if ( name == "Editx2rel" || name == "Editv1" )
        ui->inpV1->setValue( data.toDouble() );
    else if ( name == "Edity2rel" || name == "Editv2" )
        ui->inpV2->setValue( data.toDouble() );
    else if ( name == "Editz2rel" || name == "Editv3" )
        ui->inpV3->setValue( data.toDouble() );
    else if ( name == "Edithklmax" )
        ui->inpHKLmax->setValue( data.toInt() );
    else if ( name == "RadioButtonQ1" )
        ui->radQ1->setChecked( data.toInt() != 0 );
    else if ( name == "RadioButtonQ2" )
        ui->radQ2->setChecked( data.toInt() != 0 );
    else if ( name == "RadioButtonQ4" )
        ui->radQ4->setChecked( data.toInt() != 0 );
    else if ( name == "ExpandImage" )
        ui->togExpandImage->setChecked( data.toInt() != 0 );
    else if ( name == "EditGridPoints" )
        ui->inpGridPoints->setValue( data.toInt() );
    else if ( name == "EditCenterX" )
        ui->inpBCenterX->setValue( data.toInt() );
    else if ( name == "EditCenterY" )
        ui->inpBCenterY->setValue( data.toInt() );
    // Falsche / unbekannte Bezeichnungen
    else
    {
        static QStringList tmp;
        if ( tmp.size() == 0 )
            tmp = calcGui->paramsForMethod(m,false,true,false);
        if ( tmp.contains(name,Qt::CaseInsensitive) )
        {
            for ( int i=0; i<tmp.size(); i++ )
                if ( tmp[i].compare(name,Qt::CaseInsensitive) == 0 )
                {
                    qDebug() << "setScatterValue" << name << "-->" << tmp[i];
                    return;
                }
        }
        else
            qDebug() << "setScatterValue" << name << "unbekannt" << data;
    }
}


/**
 * @brief SC_MainGUI::on_togOnlyNewWindow_toggled
 * @param checked - state of the toggle
 * Enables / Disables the Use new/last Image window function.
 */
void SC_MainGUI::on_togOnlyNewWindow_toggled(bool checked)
{
    ui->radNewImageCal->setChecked(checked);
    ui->radNewImageCfg->setChecked(checked);
    ui->radNewImageCal->setDisabled(checked);
    ui->radNewImageCfg->setDisabled(checked);
    ui->radLastImageCal->setDisabled(checked);
    ui->radLastImageCfg->setDisabled(checked);
}

void SC_MainGUI::on_togLimitRuntime_toggled(bool checked)
{
    ui->inpLimitRuntime->setEnabled(checked);
}


bool SC_MainGUI::copyParamsToFitTable()
{
    QStringList slParams = calcGui->paramsForMethod( ui->tabMethods->currentIndex(), false, false, true );
    QStringList slHdr;
    slHdr << "Used" << "Min" << "Start" << "Max" << "Result";

    ui->tblFitValues->clear();
    ui->tblFitValues->setRowCount( slParams.size() );
    ui->tblFitValues->setColumnCount( slHdr.size() );
    ui->tblFitValues->setHorizontalHeaderLabels( slHdr );
    ui->tblFitValues->setVerticalHeaderLabels( slParams );

    curMethod = calcGui->index2methodname(ui->tabMethods->currentIndex());
    ui->lblFitUsedMethod->setText( curMethod );
    _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
    if ( p2f == nullptr )
    {
        p2f = new _param2fitval;
        method2fitparams.insert( curMethod, p2f );
    }
    ui->tblFitValues->blockSignals(true);
    oneFitParamUsed = false;
    foreach (QString p, slParams)
    {
        _fitLimits *fl = p2f->value(p,nullptr);
        if ( fl == nullptr )
        {
            fl = new _fitLimits;
            fl->orgval = fl->fitstart = calcGui->currentParamValueDbl( curMethod, p );
            fl->fitType  = _fitTypes::fitNone;
            fl->fitvalid = false;
            fl->min = 0;
            fl->max = 1000;
            fl->used = false;
            p2f->insert( p, fl );
            //qDebug() << "FIT: create new" << p << fl->fitstart;
        }
        else
        {
            fl->orgval = fl->fitstart = calcGui->currentParamValueDbl( curMethod, p );
            //qDebug() << "FIT: use old" << p << fl->fitstart;
        }

        bool cnt;
        double min, max;
        if ( calcGui->limitsOfParamValue( curMethod, p, min, max, cnt ) )
        {
            if ( fl->fitType == _fitTypes::fitNone )
            {   // Nur, wenn die Grenzen nicht definiert sind, diese aus der internen Struktur holen.
                // Sonst sind diese schon angepasst und beim Start aus der Registry gelesen worden.
                if ( cnt )
                    fl->fitType = _fitTypes::fitCbs;
                else
                    fl->fitType = _fitTypes::fitNumeric;
                fl->min  = min;
                fl->max  = max;
                fl->used = true;
            }
            else
            {   // Zur Sicherheit die Grenzwerte prüfen und die internen Daten nicht überschreiten
                if ( fl->min < min ) fl->min = min;
                if ( fl->max > max ) fl->max = max;
            }
        }

        int currow = slParams.indexOf(p);
        // CheckBox centered
        tblCheckBox *tog = new tblCheckBox( fl->used );
        connect( tog->tog(), SIGNAL(toggled(bool)), this, SLOT(tblFitUsed_toggled(bool)) );
        ui->tblFitValues->setCellWidget( currow, 0, tog );
        oneFitParamUsed |= fl->used;
        // other numerical values
        ui->tblFitValues->setItem( currow, 1, new QTableWidgetItem(QString::number(fl->min)) );
        ui->tblFitValues->setItem( currow, 2, new QTableWidgetItem(QString::number(fl->fitstart)) );
        ui->tblFitValues->setItem( currow, 3, new QTableWidgetItem(QString::number(fl->max)) );
        // result may be undefined and is not editable
        if ( fl->fitvalid )
            ui->tblFitValues->setItem( currow, 4, new QTableWidgetItem(QString::number(fl->fitres)) );
        else
            ui->tblFitValues->setItem( currow, 4, new QTableWidgetItem("?") );
        ui->tblFitValues->item(currow,4)->setFlags(Qt::NoItemFlags | Qt::ItemIsEnabled);
        // all textual values are centered
        for ( int c=1; c<=4; c++ )
            ui->tblFitValues->item(currow,c)->setTextAlignment(Qt::AlignCenter);
    }
    ui->tblFitValues->blockSignals(false);
    ui->tblFitValues->resizeColumnsToContents();
    //qDebug() << "copyParamsToFitTable() done" << oneFitParamUsed;
    return oneFitParamUsed;
}

/**
 * @brief SC_MainGUI::on_tabMain_currentChanged
 * @param index - index of the current visible tab
 * Updates the Fit-Data-Table if it becomes visible
 */
void SC_MainGUI::on_tabMain_currentChanged(int index)
{
    if ( index != 4 /*fit*/ ) return;

    if ( fitIsRunning ) return;

    bool fl = copyParamsToFitTable();

    ui->butFitStart->setEnabled( fl && ui->cbsFitImageWindows->count()>0 );
}

void SC_MainGUI::tblFitUsed_toggled(bool checked)
{
    QCheckBox *cbs = static_cast<QCheckBox*>(sender());
    oneFitParamUsed = false;
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        if ( static_cast<tblCheckBox*>(ui->tblFitValues->cellWidget(r,0))->tog() == cbs )
        {
            _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
            if ( p2f == nullptr ) return;
            QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
            p2f->value(p)->used = checked;
        }
        oneFitParamUsed |= static_cast<tblCheckBox*>(ui->tblFitValues->cellWidget(r,0))->tog()->isChecked();
    }
    ui->butFitStart->setEnabled( oneFitParamUsed && ui->cbsFitImageWindows->count()>0 );
}


void SC_MainGUI::on_tblFitValues_itemChanged(QTableWidgetItem *item)
{
    if ( item == nullptr ) return;
    _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
    if ( p2f == nullptr ) return;
    QString p = ui->tblFitValues->verticalHeaderItem(item->row())->text();
    switch ( item->column() )
    {
    case 1: p2f->value(p)->min = item->text().toDouble(); break;
    case 2: p2f->value(p)->fitstart = item->text().toDouble(); break;
    case 3: p2f->value(p)->max = item->text().toDouble(); break;
    }
}

bool SC_MainGUI::myProgressLogging( char *msg )
{
    if ( current ) current->updateLogList( QString(msg) );
    return _bAbbruch;
}

void SC_MainGUI::updateLogList( QString msg )
{
    if ( msg.startsWith("@REP=") )
    {
        ui->lblFitCurIter->setText( msg.mid(5) );
        if ( !bIgnoreUpdates )
            qApp->processEvents();
        return;
    }
    // Alles in das Logfile, sofern es da ist
    if ( fileFitLog )
    {
        fileFitLog->write( qPrintable(msg+EOL) );
        fileFitLog->flush();
    }

    if ( fitMinimalOutput )
    {   // Filtern bestimmter Meldungen
        if ( msg.startsWith("FQS(") ) return;
        if ( msg.startsWith("----- rtol=") ) return;
        if ( msg.startsWith("iter=") ) return;
        if ( msg.startsWith("amotry(") ) return;
        if ( msg.startsWith("===========") ) return;
        //if ( msg.contains(" -> ") ) return; // Variablen-Änderungen
        //if ( msg.startsWith("Overall:") ) return; // Anzeigen der Zeit dieser Iteration
        //qDebug() << "fitMinimalOutput" << msg;
    }
    ui->lisFitLog->addItem( msg );
    while ( ui->lisFitLog->count() > 1000 ) // Diese Listen können Probleme machen, wenn
        ui->lisFitLog->takeItem(0);         //  die Anzahl der Items zu groß wird.
    if ( !bIgnoreUpdates )
    {
        ui->lisFitLog->scrollToBottom();
        qApp->processEvents();
    }
}

void SC_MainGUI::on_butSaveLog_clicked()
{
    QString fn = QFileDialog::getSaveFileName( this, "Save Loggings", dataPath, "Logfiles (*.log)", nullptr, QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".log",Qt::CaseInsensitive) ) fn += ".log";
    if ( fileFitLogName.isEmpty() )
    {
        QFile f(fn);
        if ( f.open(QIODevice::WriteOnly) )
        {
            for ( int i=0; i<ui->lisFitLog->count(); i++ )
                f.write( qPrintable(ui->lisFitLog->item(i)->text()+EOL) );
            f.close();
            qDebug() << "Logfile saved from list";
        }
        else
            qDebug() << f.errorString();
    }
    else
    {
        QFile::remove(fn);  // copy did not overwrite
        QFile::copy( fileFitLogName, fn );
        statusBar()->showMessage( "Logfile copied to "+fn, 5000 );
    }
}

void SC_MainGUI::on_butFitStart_clicked()
{
    ui->butFitStart->setEnabled(false);
    bool savAutoEna = ui->butFitAutomatic->isEnabled();
    //ui->butFitAutomatic->setEnabled(false);
    //ui->butFitUseResult->setEnabled(false);
    //ui->butShowResiduen->setEnabled(false);

    // Automatic button enabled => the start button was pressed
    if ( savAutoEna ) fitIsRunning = true;

    // Get datapointer of selected image
    curFitImage = nullptr;
    for ( int i=0; i<images.size(); i++ )
    {
        if ( images.at(i)->windowTitle() == ui->cbsFitImageWindows->currentText() )
        {
            curFitImage = images.at(i);
            break;
        }
    }
    if ( curFitImage == nullptr )
    {
        ui->lisFitLog->addItem("No image found '"+ui->cbsFitImageWindows->currentText()+"'");
        ui->butFitStart->setEnabled(true);
        //ui->butFitAutomatic->setEnabled(savAutoEna);
        return;
    }
    ui->grpFitParams->setEnabled(false);
    ui->grpFitResult->setEnabled(false);
    ui->grpFitStart->setEnabled(false);
    //qDebug() << "Start fit: used image" << curFitImage->getMetaTitle() << curFitImage->windowTitle();

    performFirstFitLoop( curFitImage );

    ui->butFitStart->setEnabled(true);
    //ui->butFitAutomatic->setEnabled(savAutoEna);
    ui->butFitUseResult->setEnabled(true);  // Diese 3 Buttons werden beim Start gesperrt
    ui->butShowResiduen->setEnabled(true);  // und müssen daher nach derm ersten Fit-lauf
    ui->butFitHistoShow->setEnabled(true);  // explizit freigegeben werden.
    ui->grpFitParams->setEnabled(true);
    ui->grpFitResult->setEnabled(true);
    ui->grpFitStart->setEnabled(true);
    if ( savAutoEna ) fitIsRunning = false;
}

void SC_MainGUI::performFirstFitLoop( widImage *fitImg )
{
    curFitImage = fitImg;

    if ( !fitMinimalOutput ) ui->lisFitLog->clear();

    if ( fitClass == nullptr )
        fitClass = new SasCalc_SimplexFit2D( calcGui );

    fitClass->setImageInfo( fitImg->xmin(), fitImg->xmax(), fitImg->ymin(), fitImg->ymax(),
                           fitImg->getFileInfos()->centerX, fitImg->getFileInfos()->centerY,
                           fitImg->dataPtr() );
    for ( int r=0; r<4; r++ )
    {
        QRect rc = fitImg->getNoFitRect(r);
        if ( rc.width() > 0 )
            calcGui->setNoFitRect( r, rc.left(), rc.top(), rc.right(), rc.bottom() );
        else
            calcGui->setNoFitRect( r, -1, -1, -1, -1 );
    }

    performOneFitLoop();
}

void SC_MainGUI::performOneFitLoop()
{
    fileFitLogName = dataPath + "/Simplex2D-Fit-Temp.log";  // wird bei updateLogList() geschrieben
    fileFitLog = new QFile(fileFitLogName);
    if ( ! fileFitLog->open(QIODevice::WriteOnly) )
    {
        qDebug() << fileFitLogName << fileFitLog->errorString();
        fileFitLog = nullptr;
    }

    curFitImage->getVarScaling( fitOrgMin, fitOrgmax );

    _bAbbruch = false;
    prepareCalculation( true, true );

    _param2fitval *parameter = method2fitparams.value(curMethod);
    QHash<QString,_fitLimits*>::const_iterator it;

    if ( lastAutoFitLine.isEmpty() )
    {
        // Put all relevant informations into the editAutoFit TextField
        //  so the user can copy them into an automatic definition file
        QString line = "";
        it = parameter->constBegin();
        while ( it != parameter->constEnd() )
        {
            if ( it.value()->used ) line += ", "+it.key();
            ++it;
        }
        line = "\nUse: " + line.mid(2) + "\nFit: ";  // das erste ", " weglassen
#ifdef USEREPETITIONS
        line += QString("Rep=%1; ").arg(ui->inpFitRepetitions->value());
#endif
        line += QString("Stp=%1; ").arg(ui->inpFitStepSize->value());
        line += QString("Iter=%1; ").arg(ui->inpFitMaxIter->value());
        line += QString("Tol=%1; ").arg(ui->inpFitTolerance->text().toDouble());
        line += "Diff<2; Kenn=-@D";
        ui->editAutoFit->appendPlainText(line);
    }

#ifdef UnusedValue
    // Perpare History-Buffer
    it = parameter->constBegin();
    while ( it != parameter->constEnd() )
    {
        if ( ! fitHistory.contains(it.key()) )
        {
            QVector<double> tmp;
            if ( fitHistory.size() > 0 )
            {
                int fhLen = fitHistory.value(fitHistory.keys().first()).size();
                if ( fhLen > 0 )
                    tmp.fill( UnusedValue, fhLen );
            }
            fitHistory.insert( it.key(), tmp );
        }
        ++it;
    }
    //qDebug() << fitHistory;
#endif

    // Update table
    _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
    if ( p2f != nullptr )
    {
        for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
        {
            QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
            ui->tblFitValues->item( r, 4 )->setText("?");
            p2f->value(p)->fitvalid = false;
        }
    }

    timeForAll = 0;
    loopsForAll   = 0;
    imgGenForAll  = 0;
    bool ok;
    double tolerance  = ui->inpFitTolerance->text().toDouble(&ok);
    if ( !ok )
    {   // TODO: bessere Fehlermeldung
        qDebug() << "Tolerance invalid:" << ui->inpFitTolerance->text();
        tolerance = 0.001;
    }
#ifdef USEREPETITIONS
    for ( int rep=0; rep<ui->inpFitRepetitions->value(); rep++ )
    {
        ui->lblFitCurIter->setText("");
        ui->lblFitCurRep->setText( QString::number(rep+1) );
        if ( rep > 0 )
        {
            updateLogList( QString("%1 ms / %2 loops / %3 img").arg(fitClass->higResTimerElapsed)
                           .arg(fitClass->repetitions).arg(fitClass->numImgCalc) );
            it = parameter->constBegin();
            while ( it != parameter->constEnd() )
            {
                if ( it.value()->used )
                {
                    updateLogList( QString("%1: %2 -> %3").arg(it.key())
                                   .arg(it.value()->fitstart).arg(it.value()->fitres) );
                    it.value()->fitstart = it.value()->fitres;
                }
                ++it;
            }
            updateLogList("=============================");
        }
#endif
        QString retinfo;
        fitClass->doSimplexFit2D( ui->inpNumCores->value(),
                                  ui->inpFitStepSize->value() * 0.01, // wie im Pascalprog TODO ???
                                  ui->inpFitMaxIter->value(),
                                  tolerance, //ui->inpFitTolerance->text().toDouble(),
                                  ui->inpFitBorder->value(),
                                  ui->togFitUseMask->isChecked() ? -1 : ui->inpFitBStop->value(),
                                  static_cast<progressLogging>(&this->myProgressLogging),
                                  method2fitparams.value(curMethod), retinfo );
        timeForAll += fitClass->higResTimerElapsed;
        loopsForAll += fitClass->repetitions;
        imgGenForAll += fitClass->numImgCalc;
        if ( ! retinfo.isEmpty() )
        {
            QMessageBox::critical( this, "Simplex 2d Fit", retinfo );
#ifdef USEREPETITIONS
            break;
#endif
        }
#ifdef USEREPETITIONS
        if ( _bAbbruch || fitClass->aborted ) break;
    }
#endif

    finishCalculation(false);

    if ( _bAbbruch || fitClass->aborted )
        updateLogList("*** Aborted ***");

    updateLogList( QString("%1 ms / %2 loops / %3 img").arg(fitClass->higResTimerElapsed)
                   .arg(fitClass->repetitions).arg(fitClass->numImgCalc) );
    updateLogList("=============================");
    it = parameter->constBegin();
    fitMeanChangePercent = 0.0;
    int fmc = 0;
    while ( it != parameter->constEnd() )
    {
        if ( it.value()->used /*&& it.value()->fitvalid*/ )
        {
            double tmp;
            if ( it.value()->orgval == 0 )
                tmp = fabs(it.value()->orgval - it.value()->fitres);
            else
                tmp = fabs(it.value()->orgval - it.value()->fitres) * 100.0 / it.value()->orgval;
            updateLogList( QString("%1: %2 -> %3 (%4%)").arg(it.key())
                           .arg(it.value()->orgval).arg(it.value()->fitres).arg(tmp) );
            //  org     =100%
            // |org-res|= ? %
            //  sum     +=|org-res|*100/org
            fitMeanChangePercent += tmp;
            fmc++;

#ifdef UnusedValue
            // History
            fitHistory[it.key()].append(it.value()->fitres);
#endif
        }
#ifdef UnusedValue
        else
            fitHistory[it.key()].append(UnusedValue);
#endif
        ++it;
    }
    fitMeanChangePercent = fitMeanChangePercent / fmc;
    updateLogList( QString("Overall: %1 ms / %2 loops / %3 img / MeanDif %4%").arg(timeForAll)
                   .arg(loopsForAll).arg(imgGenForAll).arg(fitMeanChangePercent) );

    if ( fileFitLog )
    {
        fileFitLog->close();
        fileFitLog->deleteLater();
        fileFitLog = nullptr;
    }

    // Update table
    //_param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
    if ( p2f != nullptr )
    {
        for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
        {
            QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
            ui->tblFitValues->item( r, 4 )->setText(QString::number(p2f->value(p)->fitres));
            p2f->value(p)->fitvalid = true;
            // Restore the old values in the Calculation Tab
            calcGui->updateParamValue( "", p, p2f->value(p)->orgval, Qt::black, false/*ohne debug*/ );
        }
    }
    ui->tblFitValues->resizeColumnsToContents();
    ui->lisFitLog->scrollToBottom();    // if updates are disabled during calculations
}

void SC_MainGUI::on_butFitUseResult_clicked()
{
    _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
    if ( p2f == nullptr ) return;
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
        QString val = ui->tblFitValues->item( r, 4 )->text().trimmed();
        if ( val.isEmpty() || val.contains("?") ) continue;
        calcGui->updateParamValue( "", p, val.toDouble(), Qt::black, false/*ohne debug*/ );
        ui->tblFitValues->item( r, 2 )->setText( val );
        p2f->value(p)->fitstart = p2f->value(p)->fitres;
        p2f->value(p)->orgval   = p2f->value(p)->fitres;
    }
    ui->tblFitValues->resizeColumnsToContents();
}

void SC_MainGUI::on_radQ1_toggled(bool checked)
{
    if ( checked ) ui->togExpandImage->setEnabled(true);
}

void SC_MainGUI::on_radQ2_toggled(bool checked)
{
    if ( checked ) ui->togExpandImage->setEnabled(true);
}

void SC_MainGUI::on_radQ4_toggled(bool checked)
{
    if ( checked ) ui->togExpandImage->setEnabled(false);
}


/**
 * @brief SC_MainGUI::logThreadTimer
 * Timer function to update the progress bar and to terminate the calculation thread if the
 * maximum calculation time was expired.
 */
void SC_MainGUI::logThreadTimer()
{
    //qDebug() << "   - - - QlogThread test" << calcRunTime.elapsed();
    int x, y;
    if ( ! calcGui->getLastXY(x,y) )
    {
        qDebug() << "   - - - QlogThread ret 1";
        return;
    }
    if ( x <= -10000 )
    {
        qDebug() << "   - - - QlogThread ret 2" << calcRunTime.elapsed();
        if ( _logThreadTimer->interval() < 500 )
            _logThreadTimer->setInterval(1000);
        return; // Speicher noch nicht da
    }
    if ( calcRunTime.elapsed() < 2000 && x==0 && y==0 && lastTime == 0 )
    {
        qDebug() << "   - - - QlogThread ret 3";
        return; // Sonst kommt: x=0, y=0, v=48% ganz zu Beginn ...
    }
    if ( lastTime == 0 && (x-zzmin > zzrange/2 || y-iimin > iirange/2) )
    {   // Bei einem schnellen Test zu Beginn eines weiteren Images kann es vorkommen, dass die Prepare-Routine länger
        // dauert und somit die x,y Werte noch hoch sind. Dann startet der ProgressBar mit 99% ... das sieht blöd aus.
        //qDebug() << "   - - - QlogThread Anfang" << x << y;
        return;
    }
    if ( _bAbbruch || (calcMaxTime > 0 && calcRunTime.elapsed() > calcMaxTime) )
    {
        if ( _logThreadTimer->interval() > 1000 )
        {   // Jetzt ist hier schon einmal der Stop gemacht worden, also den Hauptthread terminieren.
            // Passiert, wenn nur 1 Thread genutzt wird (z.B. beim Debug)
            qDebug() << "   - - - QlogThread: timeout, terminating calc" << calcRunTime.elapsed();
            _calcThread->terminate();
        }
        else
        {
            if ( _bAbbruch )
                qDebug() << "   - - - QlogThread Abbruch";
            else
                qDebug() << "   - - - QlogThread: timeout" << calcRunTime.elapsed();
            calcGui->endThread();
            // Da dieser Zustand länger kommen kann (wenn die Threads nicht direkt abbrechen),
            // den Timer etwas verlangsamen für den Kontroll-Output
            _logThreadTimer->setInterval(2000);
            ui->butAbbruch->setText( _bAbbruch ? "Aborting" : "Timeout" );
            ui->butAbbruch->setEnabled(false);
        }
    }
    else
    {
        if ( _logThreadTimer->interval() > 500 )
            _logThreadTimer->setInterval(100);
    }
    if ( x != lastX || y != lastY )
    {
        //int v = ( y == iimin ) ? 0 : ((y-iimin-1)*zzrange);
        //v = 100 * ( v + (x - zzmin) ) / (1.0*zzrange*iirange);
        int v = 100.0 * (x - zzmin) / (1.0*zzrange);
        lastX = x;  // Äußere Schleife
        lastY = y;  // Innere Schleife
        if ( v > lastPrz )
        {
            lastPrz = v;
            //qDebug() << "   - - - QlogThread:" << x << y << v << "%, t=" << calcRunTime.elapsed()
            //         << "ii" << iimin << iirange << "zz" << zzmin << zzrange;
            updateProgBar( v );
            lastTime = calcRunTime.elapsed();
        }
    }
    if ( lastPrz >= 0 &&
         calcRunTime.elapsed() - lastTime > 10000 )
    {
        lastTime = calcRunTime.elapsed();
        //qDebug() << "   - - - QlogThread:" << x << y << v << "%, t=" << calcRunTime.elapsed();
    }
}



void SC_MainGUI::on_butDataSetMask_clicked()
{
    widImage *ikws = nullptr;
    for ( int i=0; i<images.size(); i++ )
        if ( images.at(i)->windowTitle().startsWith("KWS-Image") )
        {
            ikws = images.at(i);
            break;
        }
    if ( ikws == nullptr ) return;
    for ( int i=0; i<images.size(); i++ )
        if ( ! images.at(i)->windowTitle().startsWith("KWS") ) // Ignore KWS-Image and KWS-Mask
        {
            widImage *icur = images.at(i);
            if ( icur->myWidth() != ikws->myWidth() ||
                 icur->myHeight() != ikws->myHeight() ) continue;
            // Nur gleich große Bilder bearbeiten
            double *data = new double[icur->myWidth()*icur->myHeight()];
            for ( int p=0; p<icur->myWidth()*icur->myHeight(); p++ )
                if ( ikws->dataPtr()[p] > ikws->dataPtr()[0] )
                    data[p] = icur->dataPtr()[p];
                else
                    data[p] = 0;
            //widImage* img = addImage( true, 0, icur->myWidth(), 0, icur->myHeight(), data,
            //                          "KWS-Mask from '"+icur->windowTitle()+"'", false );
            widImage* img = addImage( true, icur->xmin(), icur->xmax(), icur->ymin(), icur->ymax(), data,
                                      "KWS-Mask from '"+icur->windowTitle()+"'", false );
            img->addMetaInfo( "NbRows", QString::number(icur->myHeight()) );
            img->addMetaInfo( "NbCols", QString::number(icur->myWidth()) );
            img->addMetaInfo( "BeamPosX", QString::number(ikws->getFileInfos()->centerX) );
            img->addMetaInfo( "BeamPosY", QString::number(ikws->getFileInfos()->centerY) );
            delete[] data;
        }
}

/**
 * @brief SC_MainGUI::on_butDataCopyScaling_clicked
 * Sets the min/max of the selected window as the fixed scaling for all other windows.
 * This is not only for KWS-Images but for all types of data and for all sizes.
 */
void SC_MainGUI::on_butDataCopyScaling_clicked()
{
    if ( ui->lisDataWindows->currentItem() == nullptr )
        return;
    QString t = ui->lisDataWindows->currentItem()->text();
    //qDebug() << "on_butDataCopyScaling_clicked() Titel=" << t;
    int in = -1;
    for ( int i=0; i<images.size(); i++ )
    {
        //qDebug() << "on_butDataCopyScaling_clicked()" << i << images.at(i)->windowTitle();
        if ( t.startsWith(images.at(i)->windowTitle()) )
        {
            in = i;
            //qDebug() << "on_butDataCopyScaling_clicked() Image=" << i;
            break;
        }
    }
    if ( in < 0 ) return;
    double min, max;

    if ( useFixedScaling ) // from Automatisc fit
    {
        min = minFixedScale;
        max = maxFixedScale;
    }
    else
        images.at(in)->getVarScaling( min, max );
    //qDebug() << "on_butDataCopyScaling_clicked() Scale=" << min << max << useFixedScaling;
    for ( int i=0; i<images.size(); i++ )
        if ( i != in ) // ! images.at(i)->windowTitle().startsWith(t) )
        {
            widImage *icur = images.at(i);
            //qDebug() << "on_butDataCopyScaling_clicked() Set=" << icur->windowTitle();
            icur->setFixScaling(min,max);
        }
}

void SC_MainGUI::on_butDataFindCenter_clicked()
{
    QList<QListWidgetItem*> items = ui->lisDataWindows->selectedItems();
    if ( items.count() == 0 ) return;
    widImage *ikws = nullptr;
    for ( int i=0; i<images.size(); i++ )
    {
        //qDebug() << "CHECK" << items.at(0)->text() << i << images.at(i)->windowTitle();
        if ( items.at(0)->text().startsWith( images.at(i)->windowTitle() ) )
        {
            ikws = images.at(i);
            break;
        }
    }
    if ( ikws == nullptr ) return;

    int x, y;
    SC_ReadData::findBeamCenter( ikws, x, y );

    ikws->addMetaInfo( "BeamPosX", QString::number(x) );
    ikws->addMetaInfo( "BeamPosY", QString::number(y) );

    int flg = QMessageBox::question( this, "Calculated beam center",
                                     QString("The calculated beam center is %1/%2\n"
                                             "For the simulation the values are %3/%4\n"
                                             "Should they be set in the input fields?")
                                     .arg(x).arg(y).arg(x-ikws->myWidth()/2.).arg(y-ikws->myHeight()/2.),
                                     QMessageBox::Yes, QMessageBox::No );
    if ( flg == QMessageBox::Yes )
    {
        ui->inpBCenterX->setValue(x-ikws->myWidth()/2.);
        ui->inpBCenterY->setValue(y-ikws->myHeight()/2.);
    }
}

void SC_MainGUI::on_togFitUseMask_toggled(bool checked)
{
    if ( checked )
    {
        ui->inpFitBStop->setValue(0);
        ui->inpFitBorder->setValue(0);
    }
    ui->inpFitBStop->setEnabled( !checked );
    ui->inpFitBorder->setEnabled( !checked );
}

void SC_MainGUI::on_butFitConfig_clicked()
{
    QString fn = QFileDialog::getOpenFileName( this, "Parameter configuration", configParamsFile,
                                               "Ini-Files (*.ini)", nullptr, QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;

    QSettings sets(SETT_APP,SETT_GUI);
    sets.setValue( "ConfigParamFile", fn );

    if ( fn.size() > 55 )
        ui->lblFitConfigFile->setText("("+fn.left(3) + "..." + fn.right(49)+")");
    else
        ui->lblFitConfigFile->setText("("+fn+")");

    QMessageBox::information( this, "Parameter configuration",
                              "The file\n"+fn+"\nis used at the next start of this program.",
                              QMessageBox::Ok );
}



void SC_MainGUI::autoProcessingTimer()
{
    qDebug() << "AutoProc: autoProcessingTimer()";

    while ( true )
    {
        if ( autoProcessingFile->atEnd() )
        {
            autoProcessingFile->close();
            qDebug() << "AutoProc: finished";
            bIgnoreRecalc = false;  // Init done AutoProc-Finish
            return;
        }
        QString line = autoProcessingFile->readLine().trimmed();
        if ( line.startsWith("#") )
        {
            qDebug() << "AutoProc:" << line;
            continue;   // Comment
        }
        if ( line.isEmpty() ) continue;         // Ignore empty lines

        qDebug() << "AutoProc:" << line;

        if ( line.startsWith("TAB ") )
        {
            bool found = false;
            for ( int i=0; i<ui->tabMain->count(); i++ )
                if ( ui->tabMain->tabText(i).startsWith(line.mid(4)) )
                {
                    ui->tabMain->setCurrentIndex(i);
                    found = true;
                    break;
                }
            if ( found ) continue;
            qDebug() << "AutoProc-ERROR: invalid tab name:";
            for ( int i=0; i<ui->tabMain->count(); i++ )
                qDebug() << "    :" << ui->tabMain->tabText(i);
            autoProcessingFile->close();
            bIgnoreRecalc = false;  // Init done AutoProc-Error
            return;
        }

        else if ( line.startsWith("OPEN MEAS FILE ") )
        {
            QString fn;
            if ( line.mid(15).startsWith("DEFAULT") )
            {
                QSettings data(SETT_APP,SETT_GUI);
                fn = data.value("LastImage",dataPath).toString();
            }
            else
                fn = line.mid(15);
            //qDebug() << "AutoProc: Load meas file" << fn;
            if ( local_OpenMeasFile(fn,nullptr) ) continue;
            qDebug() << "AutoProc-ERROR: meas file not found / unknown format";
            autoProcessingFile->close();
            bIgnoreRecalc = false;  // Init done AutoProc-Error
            return;
        }

        else if ( line.startsWith("LOAD PARAMS ") )
        {
            QString fn;
            if ( line.mid(12).startsWith("DEFAULT") )
            {
                QSettings data(SETT_APP,SETT_PAR);
                fn = data.value("LastParam",".").toString();
            }
            else
                fn = line.mid(12);
            //qDebug() << "AutoProc: Load params" << fn;
            if ( QFile::exists(fn) )
                local_Load_all_Parameters(fn,"");
            else
            {
                qDebug() << "AutoProc-ERROR: parameter file not found" << fn;
                autoProcessingFile->close();
                bIgnoreRecalc = false;  // Init done AutoProc-Error
                return;
            }
        }

        else if ( line.startsWith("USE QMAX FROM DATA") )
        {
            on_butUseQMax_clicked();
        }

        else if ( line.startsWith("THREADS ") )
        {
            int t;
            t = line.mid(7).trimmed().toInt();
            //qDebug() << "AutoProc: Threads" << t;
            ui->inpNumCores->setValue(t);
        }

        else if ( line.startsWith("DOCALC") )
        {
            QTimer::singleShot( 100, this, SLOT(on_butCalc_clicked()) );
            // Autoproc is restarted after calculation
            return;
        }

        else if ( line.startsWith("DOFFT") )
        {
            QTimer::singleShot( 100, this, SLOT(on_butIFFT_clicked()) );
            // Autoproc is restarted after calculation
            return;
        }

        else if ( line.startsWith("TEST ") )
        {
            switch ( line.mid(5).trimmed().toInt() )
            {
            case 1: ui->radTestImg1->setChecked(true); break;
            case 2: ui->radTestImg2->setChecked(true); break;
            case 3: ui->radTestImg3->setChecked(true); break;
            case 4: ui->radTestImg4->setChecked(true); break;
            case 5: ui->radTestImg5->setChecked(true); break;
            case 6: ui->radTestImg6->setChecked(true); break;
            default: continue; // while loop
            }
            on_butTestGo_clicked();
            ui->cbsFFTWindows->setCurrentIndex( ui->cbsFFTWindows->count()-1 );
        }

        break; // endless loop
    } // while true
    qDebug() << "AutoProc: restart timer" << bIgnoreRecalc;
    QTimer::singleShot( 200, this, SLOT(autoProcessingTimer()) );
}


void SC_MainGUI::on_butFitAutomatic_clicked()
{
    // Abfragen des Input-Files und diverser Flags
    dlgConfigAutoFit *dlgCfgAuto = new dlgConfigAutoFit( this );
    if ( dlgCfgAuto->exec() != QDialog::Accepted ) return;
    QString fnInput = dlgCfgAuto->getFilename();

    QFile finp( fnInput );
    if ( ! finp.open(QIODevice::ReadOnly) )
    {
        qDebug() << fnInput << finp.errorString();
        return;
    }

    // Basispfad für alle Outputs, immer im gleichen Verzeichnis wie das Inputfile
    QString basePath = QFileInfo(fnInput).absolutePath() + "/";

    if ( dlgCfgAuto->isDeleteFiles() )
    {   // Jetzt werden alle Files außer dem Inputfile gelöscht
        QStringList sl = QDir(basePath).entryList( QDir::Files | QDir::NoDotAndDotDot );
        sl.removeOne(QFileInfo(fnInput).fileName());
        foreach ( QString f, sl )
        {
            if ( ! QFile::remove(basePath+f) )
            {
                qDebug() << "Remove fail" << basePath+f;
                return;
            }
        }
    }
    fitMinimalOutput = dlgCfgAuto->isMinimalOutput();
    if ( fitMinimalOutput ) ui->lisFitLog->clear();

    ui->editAutoFit->appendPlainText("# Automatic start from file "
                                     +QDateTime::currentDateTime().toString("dd.MMM.yyyy hh:mm:ss")
                                     +"\n# "+fnInput);

    ui->butFitAutomatic->setEnabled(false);
    fitIsRunning = true;
    QString kenn = QDate::currentDate().toString("-yyyyMMdd-");
    bool firstScanOfFile = true;
    QString dirMask = "*";
    QStringList usedImages; // wird während des Abarbeitens geleert, daher kein Index nötig
    widImage *lastFitImage = nullptr;
    QFile *fglobLog = nullptr;
    QFile *ftex = nullptr;
    QTextStream tsTex;
    QDateTime startZeit;

    QString latexImg1, latexImg2, latexImg3;    // nebeneinander

    typedef struct
    {
        double anf, end;
    } _globValues;
    QHash<QString/*param*/,_globValues*> param2values;

    double globSummeTimeForAll = 0;
    int    globSummeLoopsForAll = 0;
    int    globSummeImgGenForAll = 0;

    QStringList slInputLines;
    QString infoText="";

    useFixedScaling = false;

    bool firstFit = true; // die Vorarbeiten werden nur einmal gemacht

    while ( true )
    {
        double summeTimeForAll = 0;
        int    summeLoopsForAll = 0;
        int    summeImgGenForAll = 0;

        latexImg1="";
        latexImg2="";
        latexImg3="";

        while ( ! finp.atEnd() )
        {
            QString line = finp.readLine().trimmed();
            if ( line.isEmpty() ) continue;         // Ignore empty lines
            if ( firstScanOfFile ) slInputLines << line;
            if ( line.startsWith("#") ) continue;   // Comment

#ifdef Q_OS_LINUX
            line.replace( "C:/SimLab/", "/home/wagener/" );
#endif
            qDebug() << "AutoFit:" << line;

            if ( line.startsWith("EOF") ) break; // Hier ist das File zu Ende
            // Das erspart beim Testen das Auskommentieren der restlichen Zeilen

            if ( firstScanOfFile )
            {   // Einige Einträge werden nur beim ersten Durchlauf interpretiert,
                // also beim ersten anzufittenden File aus einer Liste.

                if ( line.startsWith("Scale:") )
                {   // "Scale: <min>, <max>
                    // Mit diesen Werten werden alle Images mit der gleichen Skalierung gespeichert.
                    // Fehlt dieser Eintrag, wird das Orginalbild auf Min/Max skaliert gespeichert und
                    // das gefittete Bild bekommt die gleiche Skala.
                    QStringList sl = line.mid(6).trimmed().split(",");
                    if ( sl.size() != 2 ) continue;
                    bool ok1, ok2;
                    double tmp1 = sl[0].toDouble(&ok1);
                    double tmp2 = sl[1].toDouble(&ok2);
                    if ( ok1 && ok2 )
                    {
                        useFixedScaling = true;
                        minFixedScale = tmp1;
                        maxFixedScale = tmp2;
                    }
                    continue;
                }

                if ( line.startsWith("GlobLog:") )
                {   // "GlobLog: <filename>  Globales Logfile für alle Informationen
                    if ( fglobLog != nullptr ) continue;
                    fglobLog = new QFile( basePath + line.mid(8).trimmed() );
                    if ( ! fglobLog->open(QIODevice::WriteOnly) )
                    {
                        qDebug() << fglobLog->fileName() << fglobLog->errorString();
                        fglobLog->deleteLater();
                        fglobLog = nullptr;
                        continue;
                    }
                    startZeit = QDateTime::currentDateTime();
                    fglobLog->write(qPrintable("Start: "+startZeit.toString("dd.MMM.yyyy hh:MM:ss")+EOL));
                    if ( dlgCfgAuto->isLatexEnabled() )
                    {
                        ftex = new QFile( basePath + "latex-output.tex" );
                        if ( ! ftex->open(QIODevice::WriteOnly) )
                        {
                            qDebug() << ftex->fileName() << ftex->errorString();
                            ftex->deleteLater();
                            ftex = nullptr;
                        }
                        else
                        {
                            static QStringList slTextBase = { "% !TEX TS-program = pdflatex",
                                                              "% !TEX encoding = UTF-8 Unicode",
                                                              "% Author: Michael Wagener",
                                                              "\\documentclass[11pt]{article}",
                                                              "\\usepackage{tikz}",
                                                              "\\usetikzlibrary{arrows}",
                                                              "\\usepackage[utf8]{inputenc}",
                                                              "\\usepackage{geometry}",
                                                              "\\geometry{a4paper}",
                                                              "\\geometry{margin=15mm}",
                                                              "\\usepackage{booktabs}",
                                                              "\\usepackage{listings}",
                                                              "\\usepackage{color}",
                                                              "\\definecolor{rowcolor}{rgb}{0.94, 0.97, 1.0}",
                                                              "\\author{Michael Wagener, JCNS-1}",
                                                              "\\title{SAS Crystal - Simplex Fit 2D}",
                                                              "\\usepackage{fancyhdr}",
                                                              "\\pagestyle{fancy}",
                                                              "\\renewcommand{\\headrulewidth}{0pt}",
                                                              "\\lhead{}\\chead{\\textbf{SAS Crystal Simplex Fit 2D}}\\rhead{\\today}",
                                                              "\\lfoot{Michael Wagener}\\cfoot{\\thepage}\\rfoot{JCNS-1}",
                                                              "\\usepackage{sectsty}",
                                                              "\\allsectionsfont{\\sffamily\\mdseries\\upshape}",
                                                              "\\usepackage{longtable}",
                                                              "\\usepackage{multirow}",
                                                              "\\usepackage{colortbl}",
                                                              "\\usepackage{tabularx}",
                                                              "\\newcolumntype{L}[1]{>{\\raggedright\\arraybackslash}p{#1}}",
                                                              "\\newcolumntype{C}[1]{>{\\centering\\arraybackslash}p{#1}}",
                                                              "\\newcolumntype{R}[1]{>{\\raggedleft\\arraybackslash}p{#1}}",
                                                              "\\setlength{\\tabcolsep}{1mm}",
                                                              "\\begin{document}" };
                            tsTex.setDevice(ftex);
                            tsTex.setCodec("UTF-8");
                            foreach ( QString s, slTextBase )
                            {
                                tsTex << s << EOL;
                                //ftex->write(qPrintable(s+EOL));
                            }
                            //ftex->write(qPrintable("\\section{Simplex 2D Fit - "+QDate::currentDate().toString("dd. MMM. yyyy")+"}"+EOL));
                            fglobLog->write("Use LaTeX-Output." EOL);
                        }
                    } // if ( dlgCfgAuto->isLatexEnabled() )
                    continue;
                } // if ( line.startsWith("GlobLog:") )

                if ( line.startsWith("Param:") )
                {   // "Param: <filepath> Datei mit allen Parametern
                    if ( fglobLog ) fglobLog->write(qPrintable("Load Paramfile "+line.mid(6).trimmed()+EOL) );
                    local_Load_all_Parameters( line.mid(6).trimmed(), "" ); // alle, damit auch die globalen Daten
                    //fitIsRunning = false;
                    //on_tabMain_currentChanged(4);   // Noch aus den Eingabefeldern in die Fit-Tabelle übernehmen
                    //fitIsRunning = true;
                    copyParamsToFitTable(); // Fit parameter later enabled
                    continue;
                }

                if ( line.startsWith("DirMask:") )
                {   // "DirMask: *.dat" Filemask für die Verzeichnis-Suche
                    dirMask = line.mid(8).trimmed();
                    if ( fglobLog ) fglobLog->write(qPrintable("Set directory mask to "+dirMask+EOL) );
                    //qDebug() << "AutoFit: DirMask" << dirMask;
                    continue;
                }

                if ( line.startsWith("DirUp:") )
                {   // "DirUp: <dir>" Ganzes Verzeichnis in alphanumerisch aufsteigender Reihenfolge
                    QDir d(line.mid(6).trimmed());
                    QFileInfoList fil = d.entryInfoList( QStringList()<<dirMask,
                                                         QDir::Files,
                                                         QDir::Name | QDir::LocaleAware | QDir::IgnoreCase );
                    foreach ( QFileInfo fi, fil )
                        usedImages << fi.absoluteFilePath();
                    qDebug() << "AutoFit: DirUp" << usedImages.size() << usedImages.first();
                    if ( fglobLog ) fglobLog->write(qPrintable(QString("Dir (Up) %1 files found.").arg(fil.size())+EOL) );
                    continue;
                }

                if ( line.startsWith("DirDown:") )
                {   // "DirDown: <dir>" Ganzes Verzeichnis in alphanumerisch absteigender Reihenfolge
                    QDir d(line.mid(8).trimmed());
                    QFileInfoList fil = d.entryInfoList( QStringList()<<dirMask,
                                                         QDir::Files,
                                                         QDir::Name | QDir::LocaleAware | QDir::IgnoreCase
                                                         | QDir::Reversed );
                    foreach ( QFileInfo fi, fil )
                        usedImages << fi.absoluteFilePath();
                    qDebug() << "AutoFit: DirDown" << usedImages.size() << usedImages.first();
                    if ( fglobLog ) fglobLog->write(qPrintable(QString("Dir (Down) %1 files found.").arg(fil.size())+EOL) );
                    continue;
                }

                if ( line.startsWith("File:") )
                {   // "File: <filepath>" Angabe von einzelnen Dateien (mehrfacher Eintrag möglich)
                    usedImages.append( line.mid(5).trimmed() );
                    //qDebug() << "AutoFit: UsedFiles" << usedImages;
                    if ( fglobLog ) fglobLog->write(qPrintable("File added "+line.mid(5).trimmed()+EOL) );
                    continue;
                }

                if ( line.startsWith("Limits:") )
                {   // "Limits: Base; 0; 1000" Angabe von Grenzen der Parameter
                    QStringList sl = line.mid(7).trimmed().split(";");
                    for ( int i=0; i<sl.size(); i++ )
                        sl[i] = sl[i].trimmed();
                    _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
                    if ( p2f != nullptr )
                    {
                        for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
                        {
                            QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
                            if ( p == sl[0] )
                            {
                                if ( sl.size() > 1 )
                                {   // Minimum
                                    ui->tblFitValues->item(r,1)->setText(sl[1]);
                                    p2f->value(p)->min = sl[1].toDouble();
                                }
                                if ( sl.size() > 2 )
                                {   // Maximum
                                    ui->tblFitValues->item(r,3)->setText(sl[2]);
                                    p2f->value(p)->max = sl[2].toDouble();
                                }
                                break;
                            }
                        }
                    }
                    continue;
                }

            } // if ( firstScanOfFile )

            if ( line.startsWith("Info:") )
            {   // "Info: <text>" information text for user
                infoText = line.mid(5).trimmed();
                ui->lblAutofitInfo->setText(infoText);
                continue;
            }

            if ( line.startsWith("Threads:") )
            {   // "Threads: <n>"  Number of threads to use (0=GPU)
                bool ok;
                int t = line.mid(8).trimmed().toInt(&ok);
                if ( !ok )
                {
                    qDebug() << "         Invalid number.";
                    continue;
                }
                if ( t < ui->inpNumCores->minimum() || t > ui->inpNumCores->maximum() )
                {
                    qDebug() << "         Number out of range, set to max";
                    t = ui->inpNumCores->maximum();
                }
                ui->inpNumCores->setValue( t );
                continue;
            }

            if ( line.startsWith("Use:") )
            {   // "Use: Base, I0" Definiert die zu fittenden Variablen
                lastAutoFitLine = "\n" + line;
                QStringList slNames = line.mid(4).trimmed().split(",");
                for ( int i=0; i<slNames.size(); i++ )
                    slNames[i] = slNames[i].trimmed();
                if ( fglobLog ) fglobLog->write(qPrintable("Used Parameter: "+slNames.join(",")+EOL) );
                //qDebug() << slNames;
                //if ( firstScanOfFile && ftex != nullptr ) ftex->write(qPrintable("Use: "+slNames.join(", ")+EOL));
                _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
                if ( p2f != nullptr )
                {
                    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
                    {
                        QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
                        tblCheckBox *tog = static_cast<tblCheckBox*>(ui->tblFitValues->cellWidget(r,0));
                        p2f->value(p)->used = slNames.contains(p);
                        tog->tog()->setChecked(p2f->value(p)->used);
                        if ( p2f->value(p)->used )
                        {
                            if ( ! param2values.contains(p) )
                            {
                                _globValues *gv = new _globValues;
                                gv->anf = p2f->value(p)->orgval;
                                param2values.insert( p, gv );
                            }
                        }
                    }
                }
                continue;
            } // "Use:"

            // Neue Befehle (ANF) - zum (auch bedingtem) Setzen von Variablen

            if ( line.startsWith("SetVarAnf:") )
            {   // "SetVarAnf: <string> = <name> = <wert>"
                // Diese Variable wird für jede Datei deren Filename mit <string> beginnt
                // an dieser Stelle fest gesetzt.
                // Damit können z.B. verschiedene Drehungen beim Fit betrachtet werden.
                if ( usedImages.size() == 0 ) continue;  // Keine Fileliste: ignorieren
                QStringList sl = line.mid(10).split("=");
                if ( sl.size() != 3 ) continue; // Syntaxfehler: Zeile ignorieren
                if ( ! QFileInfo(usedImages.first()).baseName().startsWith(sl[0].trimmed()) )
                    continue;  // Nicht der passende Filename
                if ( calcGui->updateParamValue( "", sl[1].trimmed(), sl[2].trimmed().toDouble(), Qt::black, false ) )
                {   // Parameter aktualisiert
                    if ( fglobLog ) fglobLog->write(qPrintable("SetVar: "+sl.join(" = ")+EOL) );
                }
                else
                {   // Fehlermeldung
                    qDebug() << "         SetVarAnf FAILED." << sl;
                    if ( fglobLog ) fglobLog->write(qPrintable("SetVarAnf FAILED ("+sl[0]+")"+EOL) );
                }
                continue;
            }

            if ( line.startsWith("SetVarEnd:") )
            {   // "SetVarEnd: <string> = <name> = <wert>"
                // Diese Variable wird für jede Datei deren Filename mit <string> endet
                // an dieser Stelle fest gesetzt.
                // Damit können z.B. verschiedene Drehungen beim Fit betrachtet werden.
                if ( usedImages.size() == 0 ) continue;  // Keine Fileliste: ignorieren
                QStringList sl = line.mid(10).split("=");
                if ( sl.size() != 3 ) continue; // Syntaxfehler: Zeile ignorieren
                if ( ! QFileInfo(usedImages.first()).baseName().endsWith(sl[0].trimmed()) )
                    continue;  // Nicht der passende Filename
                if ( calcGui->updateParamValue( "", sl[1].trimmed(), sl[2].trimmed().toDouble(), Qt::black, false ) )
                {   // Parameter aktualisiert
                    if ( fglobLog ) fglobLog->write(qPrintable("SetVar: "+sl.join(" = ")+EOL) );
                }
                else
                {   // Fehlermeldung
                    qDebug() << "         SetVarEnd FAILED." << sl;
                    if ( fglobLog ) fglobLog->write(qPrintable("SetVarEnd FAILED ("+sl[0]+")"+EOL) );
                }
                continue;
            }

            if ( line.startsWith("SetVar:") )
            {   // "SetVar: <name> = <wert>"
                // Diese Variable wird für jede Datei an dieser Stelle fest gesetzt.
                // Damit können z.B. verschiedene Drehungen beim Fit betrachtet werden.
                QStringList sl = line.mid(7).split("=");
                if ( sl.size() != 2 ) continue; // Syntaxfehler: Zeile ignorieren
                if ( calcGui->updateParamValue( "", sl[0].trimmed(), sl[1].trimmed().toDouble(), Qt::black, false ) )
                {   // Parameter aktualisiert
                    if ( fglobLog ) fglobLog->write(qPrintable("SetVar: "+sl.join(" = ")+EOL) );
                }
                else
                {   // Fehlermeldung
                    qDebug() << "         SetVar FAILED." << sl;
                    if ( fglobLog ) fglobLog->write(qPrintable("SetVar FAILED ("+sl[0]+")"+EOL) );
                }
                continue;
            }

            // Neue Befehle (END)

            if ( line.startsWith("Fit:") )
            {   // "Fit: Rep=10; Stp=3.0; Iter=20; Tol=0.0010; Diff<5" Start of Fit
                //  ...Rep=?... kann ignoriert werden
                lastAutoFitLine += "\n" + line;
                ui->editAutoFit->appendPlainText(lastAutoFitLine);
                QStringList slCmd = line.mid(4).trimmed().split(";");
                for ( int i=0; i<slCmd.size(); i++ )
                    slCmd[i] = slCmd[i].trimmed();
                //qDebug() << slCmd;
                //if ( firstScanOfFile && ftex != nullptr ) ftex->write(qPrintable("Fit: "+slCmd.join(", ")+EOL));
                double maxDif=100, stp=3.0, tol=0.001;
#ifdef USEREPETITIONS
                int rep=10;
#endif
                int iter=10;
                for ( int i=0; i<slCmd.size(); i++ )
                {
#ifdef USEREPETITIONS
                    if ( slCmd[i].startsWith("Rep=") )
                        rep = slCmd[i].mid(4).toInt();
                    else
#endif
                        if ( slCmd[i].startsWith("Stp=") )
                            stp = slCmd[i].mid(4).toDouble();
                        else if ( slCmd[i].startsWith("Iter=") )
                            iter = slCmd[i].mid(5).toInt();
                        else if ( slCmd[i].startsWith("Tol=") )
                            tol = slCmd[i].mid(4).toDouble();
                        else if ( slCmd[i].startsWith("Diff<") )
                            maxDif = slCmd[i].mid(5).toDouble();
                        else if ( slCmd[i].startsWith("Kenn=") )
                        {
                            kenn = slCmd[i].mid(5);
                            kenn = kenn.replace( "@D", QDate::currentDate().toString("yyyyMMdd") );
                            if ( usedImages.size() > 0 )
                                kenn = kenn.replace("@F", QFileInfo(usedImages.first()).baseName() );
                            else
                                kenn = kenn.replace("@F", "" );
                        }
                }
#ifdef USEREPETITIONS
                ui->inpFitRepetitions->setValue(rep);
#endif
                ui->inpFitStepSize->setValue(stp);
                ui->inpFitMaxIter->setValue(iter);
                ui->inpFitTolerance->setText(QString::number(tol));

                if ( usedImages.size() > 0 )
                {   // Jetzt muss das erste Image geladen werden ...
                    if ( lastFitImage == nullptr )
                    {   // Erst das Image suchen
                        for ( int i=0; i<images.size(); i++ )
                        {
                            if ( images.at(i)->windowTitle() == ui->cbsFitImageWindows->currentText() )
                            {
                                lastFitImage = images.at(i);
                                break;
                            }
                        }
                        if ( lastFitImage == nullptr && images.size() > 0 )
                        {   // Kein Fit-Image gefunden, Abbruch
                            if ( fglobLog )
                            {
                                fglobLog->write(qPrintable("NO FIT IMAGE FOUND" EOL) );
                                fglobLog->close();
                            }
                            if ( ftex ) ftex->close();
                            ui->lblAutofitInfo->setText("NO FIT IMG");
                            ui->butFitAutomatic->setEnabled(true);
                            fitIsRunning = false;
                            return; // TODO: Fehlermeldung
                        }
                    }
                    if ( lastFitImage == nullptr || lastFitImage->getFileInfos()->filePath != usedImages.first() )
                    {   // Jetzt muss ein neues File geladen werden ...
                        if ( lastFitImage != nullptr )
                        {
                            lastFitImage->close();
                            lastFitImage->deleteLater();
                        }
                        if ( ! local_OpenMeasFile( usedImages.first(), &lastFitImage ) )
                        {   // File nicht gefunden, Abbruch
                            if ( fglobLog )
                            {
                                fglobLog->write(qPrintable("FILE OPEN ERROR: "+usedImages.first()+EOL) );
                                fglobLog->close();
                            }
                            if ( ftex ) ftex->close();
                            ui->lblAutofitInfo->setText("NO MEAS FILE");
                            ui->butFitAutomatic->setEnabled(true);
                            fitIsRunning = false;
                            return; // TODO: Fehlermeldung
                        }
                        QString tmpfn = basePath + QFileInfo(usedImages.first()).baseName() + "_org.png";
                        if ( useFixedScaling ) lastFitImage->setFixScaling( minFixedScale, maxFixedScale );
                        lastFitImage->saveImage( tmpfn );
                        if ( fglobLog ) fglobLog->write(qPrintable("Data file: "+usedImages.first()+EOL) );
                        //for ( int r=0; r<4; r++ )
                        //{
                        //    QRect rc = lastFitImage->getNoFitRect(r);
                        //    if ( rc.width() > 0 )
                        //    {
                        //        calcGui->setNoFitRect( r, rc.left(), rc.top(), rc.right(), rc.bottom() );
                        //    }
                        //}
                        if ( ftex != nullptr && dlgCfgAuto->isLatexOrgImg() )
                        {
                            latexImg1 = tmpfn;
                            // Output at the end of the loop
                        }
                    }
                    statusBar()->showMessage(QString("Use Dataset %1 (Rest=%2)").arg(usedImages.first()).arg(usedImages.size()-1));
                } // if usedImages > 0

                if ( fglobLog ) { fglobLog->write(qPrintable("Start Fit: "+slCmd.join(",")+EOL) ); fglobLog->flush(); }

                for ( int maxloop=0; maxloop<20; maxloop++ ) // keine Endlosschleife...
                {
                    ui->lblAutofitInfo->setText(infoText+QString(" (%1/20)").arg(maxloop+1));

                    if ( firstFit )
                    {
                        performFirstFitLoop( lastFitImage );
                        firstFit = false;
                    }
                    else
                        performOneFitLoop();

                    qDebug() << "AutoFit: loop" << maxloop+1 << "MeanDif =" << fitMeanChangePercent
                             << "MaxDif =" << maxDif;
                    if ( _bAbbruch ) break;

                    summeTimeForAll += timeForAll;
                    summeLoopsForAll += loopsForAll;
                    summeImgGenForAll += imgGenForAll;

                    on_butFitUseResult_clicked();
                    if ( dlgCfgAuto->isLogEnabled() && ! fileFitLogName.isEmpty() )
                    {   // Save Logging
                        QString tmpfn = basePath + kenn + QString("%1").arg(maxloop+1,2,10,QChar('0')) + "_calc.log";
                        QFile::copy( fileFitLogName, tmpfn );
                        qDebug() << "AutoFit:" << tmpfn;
                        if ( fglobLog ) fglobLog->write(qPrintable("  Logfile: "+tmpfn+EOL));
                    }
                    if ( fglobLog )
                    {
                        fglobLog->write( qPrintable( QString("  Loop %1, MeanDif=%2, MaxDif=%3:")
                                                     .arg(maxloop+1).arg(fitMeanChangePercent).arg(maxDif)+EOL ) );
                        //if ( ftex != nullptr && fitMeanChangePercent<maxDif )
                        //    ftex->write( qPrintable( QString("%  Loops %1 used, MeanDif=%2, MaxDif=%3:")
                        //                             .arg(maxloop+1).arg(fitMeanChangePercent).arg(maxDif)+EOL ) );
                        _param2fitval *parameter = method2fitparams.value(curMethod);
                        QHash<QString,_fitLimits*>::const_iterator it = parameter->constBegin();
                        while ( it != parameter->constEnd() )
                        {
                            if ( it.value()->used /*&& it.value()->fitvalid*/ )
                            {
                                fglobLog->write(qPrintable(QString("    %1: %2 -> %3").arg(it.key())
                                                           .arg(it.value()->orgval).arg(it.value()->fitres)+EOL) );
                                if ( fitMeanChangePercent<maxDif )
                                {
                                    param2values[it.key()]->end = it.value()->fitres;
                                }
                            }
                            ++it;
                        }
                        fglobLog->flush();
                    } // if fglobLog
                    if ( fitMeanChangePercent < maxDif ) break;
                }
                lastAutoFitLine = "";
                if ( _bAbbruch ) break;
            } // "Fit:"

        } // while !atEnd
        if ( _bAbbruch ) break;

        if ( usedImages.size() >= 1 )
        {
            QString tmpfn = basePath + QFileInfo(usedImages.first()).baseName() + "_out";
            QString resifn = basePath + QFileInfo(usedImages.first()).baseName() + "_resi.png";
            on_butCalc_clicked();

            if ( useFixedScaling )
                lastUsedImage->setFixScaling( minFixedScale, maxFixedScale );
            else
                lastUsedImage->setFixScaling( fitOrgMin, fitOrgmax );

            lastUsedImage->saveImage( tmpfn+".png" );
            performSaveParamOperation( tmpfn+".ini" );
            if ( ftex != nullptr && dlgCfgAuto->isLatexResiduenImg() )
            {
                double *imgData = new double[(calcGui->maxX() - calcGui->minX()) * (calcGui->maxY() - calcGui->minY())];
                for ( int ihex=calcGui->minX(); ihex<calcGui->maxX(); ihex++ )
                    for ( int i=calcGui->minY(); i<calcGui->maxY(); i++ )
                        imgData[XY2IDX(calcGui->minX(),calcGui->maxX(),calcGui->minY(),calcGui->maxY(),ihex,i)]
                                = fitClass->getResiduenPixel(ihex,i);

                widImage* img = addImage( true, calcGui->minX(),calcGui->maxX(),calcGui->minY(),calcGui->maxY(), imgData, "Residuen", false );
                img->addMetaInfo( "Residuenplot", "True" );
                img->addMetaInfo( "NbRows", QString::number(calcGui->maxY() - calcGui->minY()) );
                img->addMetaInfo( "NbCols", QString::number(calcGui->maxX() - calcGui->minX()) );
                img->saveImageColor( resifn, tblResGlo );
                latexImg3 = resifn;
                img->close();
                img->deleteLater();
            }
            lastUsedImage->close();
            lastUsedImage->deleteLater();
            lastUsedImage = nullptr;

            if ( fglobLog )
            {
                fglobLog->write( qPrintable("Save Image  to "+tmpfn+".png"+EOL) );
                fglobLog->write( qPrintable("Save Params to "+tmpfn+".ini"+EOL) );
                if ( !latexImg3.isEmpty() )
                    fglobLog->write( qPrintable("Save Residuenplot to "+resifn+EOL) );
                fglobLog->flush();
            }
            if ( ftex )
            {
                if ( dlgCfgAuto->isLatexCalcImg() ) latexImg2 = tmpfn+".png";
                if ( dlgCfgAuto->isLatexImages() )
                {
                    QStringList hdr;
                    if ( dlgCfgAuto->isLatexOrgImg()   ) hdr<<"\\bf{Datafile}";       // latexImg1 Filename
                    if ( dlgCfgAuto->isLatexCalcImg()  ) hdr<<"\\bf{Fit output}";     // latexImg2 Filename
                    if ( dlgCfgAuto->isLatexResiduenImg() ) hdr<<"\\bf{Residuen plot}";  // latexImg3 Filename
                    if ( !firstScanOfFile ) tsTex << "\\clearpage" EOL;
                    tsTex << "\\section*{\\small{" << usedImages.first() << "}}" EOL;
                    tsTex << "\\begin{longtable}{";
                    for ( int i=0; i<hdr.size(); i++ ) tsTex << "c";
                    tsTex << "}" EOL;
                    tsTex << hdr.join(" & ") << " \\\\ \\hline" EOL;
                    QString siz = "0.45";
                    if ( hdr.size() > 2 ) siz = "0.30";
                    QStringList sl;
                    if ( !latexImg1.isEmpty() ) sl << "\\includegraphics[width="+siz+"\\textwidth]{"+latexImg1+"}";
                    if ( !latexImg2.isEmpty() ) sl << "\\includegraphics[width="+siz+"\\textwidth]{"+latexImg2+"}";
                    if ( !latexImg3.isEmpty() ) sl << "\\includegraphics[width="+siz+"\\textwidth]{"+latexImg3+"}";
                    tsTex << sl.join(" & " EOL) << EOL "\\end{longtable}" EOL;
                }
                tsTex << QString("\\centerline{Statistics: %1 ms running time / %2 iterations / %3 images generated}")
                         .arg(summeTimeForAll).arg(summeLoopsForAll).arg(summeImgGenForAll) << EOL;

                globSummeTimeForAll += summeTimeForAll;
                globSummeLoopsForAll += summeLoopsForAll;
                globSummeImgGenForAll += summeImgGenForAll;
                tsTex << "\\begin{longtable}{|C{3cm}|C{4cm}|C{4cm}|}" EOL;
                tsTex << "\\hline\\rowcolor{rowcolor} {\\bf Parameter} & {\\bf Startvalue} & {\\bf Fit-Resultvalue} \\\\" EOL;
                tsTex << "\\endfirsthead" EOL;
                tsTex << "\\hline\\rowcolor{rowcolor} {\\bf Parameter} & {\\bf Startvalue} & {\\bf Fit-Resultvalue} \\\\" EOL;
                tsTex << "\\endhead" EOL;
                tsTex << "\\hline" EOL;

                QStringList gvKeys = param2values.keys();
                gvKeys.sort();  // Sortierte Schlüssel ausgeben
                foreach ( QString k, gvKeys )
                {
                    QString n=k;
                    n.replace("_","\\_");
                    tsTex << QString("%1 & %2 & %3 \\\\ \\hline")
                             .arg(n).arg(param2values[k]->anf).arg(param2values[k]->end) << EOL;
                }
                tsTex << "\\end{longtable}" EOL;
                ftex->flush();
                param2values.clear();
            }
            firstScanOfFile = false;
            finp.seek(0);
            usedImages.takeFirst();
            if ( usedImages.size() == 0 ) break;
        }
        if ( usedImages.size() == 0 ) break;
    } // while endless

    finp.close();
    if ( fglobLog )
    {
        if ( _bAbbruch ) fglobLog->write( EOL "*** ABORTED ***" EOL );
        QDateTime endZeit = QDateTime::currentDateTime();
        fglobLog->write(qPrintable("Ende: "+endZeit.toString("dd.MMM.yyyy hh:MM:ss")+EOL));
        fglobLog->write(qPrintable("Duration: "+QString::number(startZeit.secsTo(endZeit))+" sec"+EOL));
        fglobLog->close();
    }
    QTime tGlobTime(0,0,0);
    tGlobTime = tGlobTime.addMSecs(globSummeTimeForAll);
    updateLogList( QString("Sum over all images: %1 running time / %2 iterations / %3 images generated")
                   .arg(tGlobTime.toString("hh:mm:ss.zzz")).arg(globSummeLoopsForAll).arg(globSummeImgGenForAll) );
    if ( ftex )
    {
        if ( _bAbbruch ) tsTex << EOL "{\\bf CALCULATION ABORTED - LaTeX file may be corrupted.} \\newline" EOL;
        tsTex << EOL << QString("\\centerline{Sum over all images: %1 running time / %2 iterations / %3 images generated}")
                               .arg(tGlobTime.toString("hh:mm:ss.zzz")).arg(globSummeLoopsForAll).arg(globSummeImgGenForAll) << EOL;
        if ( dlgCfgAuto->isLatexTrend() )
        {
            tsTex.flush();
            printHistory( ftex );
        }
        if ( dlgCfgAuto->isLatexInput() )
        {
            tsTex << EOL "\\clearpage" EOL "\\centerline{input command file}" EOL;
            tsTex << "\\begin{lstlisting}[frame=single, xleftmargin=0.5cm, xrightmargin=0.5cm]" EOL; // % language=Python,
            foreach ( QString s, slInputLines )
            {
                if ( dlgCfgAuto->isLatexComments() || s[0] != '#' )
                {
                    if ( s.length() > 73 )
                    {   // Zeile zu lang, umbrechen
                        if ( s.startsWith("File:") || s.startsWith("Param:") )
                        {   // Spezialfall, den Filenamen am / trennen
                            int pos = s.lastIndexOf('/') + 1;
                            tsTex << s.left(pos).toUtf8() << EOL;
                            tsTex << "        " << s.mid(pos).toUtf8() << EOL;
                        }
                        else
                        {   // Alles andere am ' ' trennen oder an der festen Position
                            int pos = s.lastIndexOf(' ');
                            if ( pos < 0 )
                                pos = 73;
                            tsTex << s.left(pos).toUtf8() << EOL;
                            tsTex << "        " << s.mid(pos).toUtf8() << EOL;
                        }
                    }
                    else
                    {   // Zeile ist kurz genug
                        tsTex << s.toUtf8() << EOL;
                    }
                }
            }
            tsTex << "\\end{lstlisting}" EOL;
        }
        tsTex << EOL "\\end{document}" EOL;

        //if ( _bAbbruch ) ftex->write(EOL "CALCULATION ABORTED - LaTeX file may be corrupted." EOL);
        //ftex->write(qPrintable(EOL+QString("\\centerline{Sum over all images: %1 running time / %2 iterations / %3 images generated}")
        //                       .arg(tGlobTime.toString("hh:mm:ss.zzz")).arg(globSummeLoopsForAll).arg(globSummeImgGenForAll)+EOL) );
        //ftex->write( EOL "\\clearpage" EOL );
        //ftex->write( "\\begin{lstlisting}[frame=single, xleftmargin=0.5cm, xrightmargin=0.5cm]" EOL ); // % language=Python,
        //foreach ( QString s, slInputLines )
        //{
        //    ftex->write( (s+EOL).toUtf8() );
        //}
        //ftex->write( "\\end{lstlisting}" EOL );
        //ftex->write( EOL "\\end{document}" EOL );
        ftex->close();
#ifndef Q_OS_WIN
        QStringList args;
        args << "-output-directory=" + QFileInfo(ftex->fileName()).absolutePath();
        args << "-interaction=batchmode";
        args << "-halt-on-error";
        args << ftex->fileName();
        //qDebug() << "PDFLATEX:" << args;
        QProcess::startDetached( "pdflatex", args );
#endif
    }
    ui->lblAutofitInfo->setText("");
    ui->butFitAutomatic->setEnabled(true);
    fitIsRunning = false;
    fitMinimalOutput = false;
    ui->editAutoFit->appendPlainText("# Automatic finished"
                                     +QDateTime::currentDateTime().toString("dd.MMM.yyyy hh:mm:ss") );
    statusBar()->showMessage( "AutoFit: finished after "+tGlobTime.toString("hh:mm:ss.zzz"), 5000 );

    if ( !_bAbbruch && dlgCfgAuto->isShowResult() )
    {
        // Zum Abschluss nochmal das Bild rechnen
        if ( ui->togAutoPosit->isChecked() )
        {   // Vorher aber die Positionen so bestimmen, dass dieses neue Bild direkt
            // neben das bisherige (anzufittende) Bild kommt.
            imgPosX = curFitImage->pos().x() + curFitImage->width();
            imgPosY = curFitImage->pos().y();
            //qDebug() << "Pos vor Calc" << imgPosX << imgPosY << curFitImage->pos() << curFitImage->size();
        }
        on_butCalc_clicked();
        //qDebug() << "Pos nach Calc" << imgPosX << imgPosY;
        for ( int i=0; i<ui->lisDataWindows->count(); i++ )
        {
            if ( ui->lisDataWindows->item(i)->text().startsWith(curFitImage->windowTitle()) )
            {
                ui->lisDataWindows->setCurrentItem( ui->lisDataWindows->item(i), QItemSelectionModel::Select );
                break;
            }
        }
        on_butDataSetMask_clicked();
        //qDebug() << "Pos nach setScale" << imgPosX << imgPosY;
        on_butDataCopyScaling_clicked();
        //qDebug() << "Pos nach copy" << imgPosX << imgPosY;
        //lower();
    }
    useFixedScaling = false; // damit später keine seltsamen Effekte auftreten.
}

void SC_MainGUI::on_butClearEditAutoFit_clicked()
{
    ui->editAutoFit->clear();
}

void SC_MainGUI::on_butShowResiduen_clicked()
{
    double *imgData = new double[(calcGui->maxX() - calcGui->minX()) * (calcGui->maxY() - calcGui->minY())];
    for ( int ihex=calcGui->minX(); ihex<calcGui->maxX(); ihex++ )
        for ( int i=calcGui->minY(); i<calcGui->maxY(); i++ )
            imgData[XY2IDX(calcGui->minX(),calcGui->maxX(),calcGui->minY(),calcGui->maxY(),ihex,i)]
                    = fitClass->getResiduenPixel(ihex,i);

    widImage* img = addImage( true, calcGui->minX(),calcGui->maxX(),calcGui->minY(),calcGui->maxY(), imgData, "Residuen", false );
    img->addMetaInfo( "Residuenplot", "True" );
    img->addMetaInfo( "NbRows", QString::number(calcGui->maxY() - calcGui->minY()) );
    img->addMetaInfo( "NbCols", QString::number(calcGui->maxX() - calcGui->minX()) );
}

void SC_MainGUI::on_butFitCurSave_clicked()
{
    savedFitParams.clear();
    savedFitParams << QString("FLAGS:%1:%2:%3:%4:%5:%6:%7")
                      .arg(ui->inpFitBStop->value())
                      .arg(ui->inpFitBorder->value())
                      .arg(ui->inpFitMaxIter->value())
                      .arg(ui->inpFitStepSize->value())
                      .arg(ui->inpFitTolerance->text())
                      .arg(ui->inpFitRepetitions->value())
                      .arg(ui->togFitUseMask->isChecked()?"true":"false");
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        QStringList tmp;
        tmp << ui->tblFitValues->verticalHeaderItem(r)->text(); // Label=Bezeichnung
        tmp << (static_cast<tblCheckBox*>(ui->tblFitValues->cellWidget(r,0))->tog()->isChecked() ? "true" : "false");
        for ( int c=1; c<ui->tblFitValues->columnCount(); c++ )
            tmp << ui->tblFitValues->item(r,c)->text();
        savedFitParams << tmp.join(":");
    }
    ui->butFitCurLoad->setEnabled(true);
    //qDebug() << savedFitParams.size() << savedFitParams;
}

void SC_MainGUI::on_butFitCurLoad_clicked()
{
    _param2fitval *p2f = method2fitparams.value(curMethod,nullptr);
    if ( p2f == nullptr ) return;
    foreach (QString s, savedFitParams)
    {
        QStringList tmp = s.split(":");
        if ( tmp[0] == "FLAGS" )
        {   // Eingabefelder
            ui->inpFitBStop->setValue(tmp[1].toInt());
            ui->inpFitBorder->setValue(tmp[2].toInt());
            ui->inpFitMaxIter->setValue(tmp[3].toInt());
            ui->inpFitStepSize->setValue(tmp[4].toDouble());
            ui->inpFitTolerance->setText(tmp[5]);
            ui->inpFitRepetitions->setValue(tmp[6].toInt());
            ui->togFitUseMask->setChecked(tmp[7]=="true");
        }
        else
        {   // Tabelleneintrag
            _fitLimits *fl = p2f->value(tmp[0]);
            for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
                if ( tmp[0] == ui->tblFitValues->verticalHeaderItem(r)->text() )
                {   // Zur Sicherheit die Zeile suchen...
                    static_cast<tblCheckBox*>(ui->tblFitValues->cellWidget(r,0))->tog()->setChecked(tmp[1]=="true");
                    for ( int c=1; c<ui->tblFitValues->columnCount(); c++ )
                        ui->tblFitValues->item(r,c)->setText(tmp[c+1]);
                    fl->used     = tmp[1]=="true";
                    fl->min      = tmp[2].toDouble();
                    fl->fitstart = tmp[3].toDouble();
                    fl->max      = tmp[4].toDouble();
                    break;
                }
        }
    }
    ui->tblFitValues->resizeColumnsToContents();
}

void SC_MainGUI::on_actionTest_Compare_data_files_triggered()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString fn1 = data.value("CompareF1",".").toString();
    QString fn2 = data.value("CompareF2",".").toString();

    fn1 = QFileDialog::getOpenFileName( this, "File 1", fn1, "Data files (*.dat)" );
    if ( fn1.isEmpty() ) return;
    fn2 = QFileDialog::getOpenFileName( this, "File 2", fn2, "Data files (*.dat)" );
    if ( fn2.isEmpty() || fn1 == fn2 ) return;

    data.setValue("CompareF1",fn1);
    data.setValue("CompareF2",fn2);

    union
    {
        double e;
        char c[1];
    } val1, val2;

    QFile f1(fn1);
    QFile f2(fn2);

    if ( !f1.open(QIODevice::ReadOnly) || !f2.open(QIODevice::ReadOnly) )
    {
        qDebug() << fn1 << f1.errorString();
        qDebug() << fn2 << f2.errorString();
        return;
    }
    QString info1 = f1.read(19);
    QString info2 = f2.read(19);
    qDebug() << fn1 << info1;
    qDebug() << fn2 << info2;
    if ( info1 != info2 )
    {
        f1.close();
        f2.close();
        qDebug() << "ERROR: size differs";
        return;
    }
    int minX, maxX, minY, maxY;
    QStringList sl = info1.split(" ",Qt::SkipEmptyParts);
    qDebug() << sl;
    if ( sl.size() != 4 )
    {
        f1.close();
        f2.close();
        qDebug() << "ERROR: no 4 size values";
        return;
    }
    minX = sl[0].toInt();
    maxX = sl[1].toInt();
    minY = sl[2].toInt();
    maxY = sl[3].toInt();

    int cnt = 0;
    for ( int y=minY; y<maxY; y++ )
    {
        for ( int x=minX; x<maxX; x++ )
        {
            f1.read( val1.c, sizeof(double) );
            f2.read( val2.c, sizeof(double) );
            if ( fabs(val1.e - val2.e) > 1.0e-7 )
            {
                qDebug() << "DIFF" << x << y << val1.e << val2.e;
                cnt++;
            }
        }
    }

    qDebug() << cnt << "Differenzen";
    f1.close();
    f2.close();
}


void SC_MainGUI::on_butFitHistoShow_clicked()
{
#ifdef UnusedValue
    printHistory( nullptr );
#else
    qDebug() << "nyi";
#endif
}

#ifdef UnusedValue
void SC_MainGUI::printHistory(QFile *ftex ) // Damit per nullptr auch die Ausgabe auf die Konsole kommt
{
    qDebug() << "History:";
    QStringList usedKeys;
    QHash< QString/*param*/, QVector<double>/*Werte*/ >::const_iterator ih = fitHistory.constBegin();
    while ( ih != fitHistory.constEnd() )
    {
        for ( int i=0; i<ih.value().size(); i++ )
            if ( ih.value()[i] != UnusedValue )
            {
                usedKeys << ih.key();
                break;
            }
        //qDebug() << ih.key() << ih.value();
        ++ih;
    }
    usedKeys.sort();
    qDebug() << usedKeys;

    QTextStream tsTex;
    if ( ftex )
    {
        tsTex.setDevice( ftex );
        tsTex << EOL "\\clearpage" EOL;
        tsTex << "The following table(s) shows each parameter value after each fit run to see the trend." EOL;
    }
    QString tabWidth;
    switch ( usedKeys.size() )  // Damit die Spaltenbreite auch bei mehr als 8 Parametern immer gleich ist
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
        tabWidth = "|C{3.1cm}";   // Größere Spalten sehen nicht gut aus...
        break;
    case 6:
        tabWidth = "|C{2.7cm}";
        break;
    case 7:
        tabWidth = "|C{2.3cm}";
        break;
    default:
        tabWidth = "|C{1.9cm}";   // Kleinere Spalten sind zu klein für die Zahlen
        break;
    }
#define MAXCOL 8
    while ( usedKeys.size() > 0 )
    {
        if ( ftex )
        {
            tsTex << "\\begin{longtable}[l]{|C{0.9cm}";  // immer Linksbündig, sieht bei mehreren Tabellen besser aus
            QString hdr = "{\\bf Nr.}";
            for ( int i=0; i<usedKeys.size() && i<MAXCOL; i++ )
            {
                tsTex << tabWidth;
                if ( hdr.size() > 0 ) hdr += " & ";
                hdr += "{\\bf "+usedKeys[i]+"}";
            }
            hdr.replace("Edit",""); // Header etwas kürzen...
            tsTex << "|}" EOL;
            //tsTex << "\\caption{History of the fit data}\\\\" EOL;
            tsTex << "\\hline\\rowcolor{rowcolor} " << hdr << " \\\\" EOL;
            tsTex << "\\endfirsthead" EOL;
            tsTex << "\\hline\\rowcolor{rowcolor} " << hdr << " \\\\" EOL;
            tsTex << "\\endhead" EOL;
            tsTex << "\\hline" EOL;
        }
        else
        {
            QString hdr = usedKeys[0];
            for ( int i=1; i<usedKeys.size() && i<MAXCOL; i++ ) hdr += ", " + usedKeys[i];
            qDebug() << hdr;
        }
        QString trenn = ( ftex != nullptr ) ? " & " : ", ";
        int vi = 0;
        while ( true )
        {
            if ( fitHistory.value(usedKeys[0]).size() <= vi ) break;
            QString line = QString("%1").arg(vi);
            for ( int k=0; k<usedKeys.size() && k<MAXCOL; k++ )
            {
                /*if ( k > 0 )*/ line += trenn;
                if ( fitHistory.value(usedKeys[k])[vi] == UnusedValue )
                    line += "-/-";
                else
                    line += QString("%1").arg(fitHistory.value(usedKeys[k])[vi],0,'f',6);
            }
            if ( ftex )
                tsTex << line << " \\\\ \\hline" EOL;
            else
                qDebug() << line;
            vi++;
        }
        if ( ftex )
        {
            tsTex << "\\end{longtable}" EOL;
        }
        for ( int k=0; usedKeys.size()>0 && k<MAXCOL; k++ ) usedKeys.takeFirst();
    }
}
#endif

void SC_MainGUI::on_butFitHistoClear_clicked()
{
#ifdef UnusedValue
    QHash< QString/*param*/, QVector<double>/*Werte*/ >::iterator ih = fitHistory.begin();
    while ( ih != fitHistory.end() )
    {
        ih.value().clear();
        ++ih;
    }
    fitHistory.clear();
#else
    qDebug() << "nyi";
#endif
}



void SC_MainGUI::compString( QLineEdit *inp, QSettings &sets, QString key, QString prmt)
{
    QString cur = inp->text();
    QString par = sets.value(key,0.0).toString();
    if ( cur == par ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = cur;
    cv->par = par;
    compWerte.insert( prmt+key, cv );
}
void SC_MainGUI::compDouble( QDoubleSpinBox *inp, QSettings &sets, QString key )
{
    double cur = inp->value();
    double par = sets.value(key,0.0).toDouble();
    if ( fabs(cur-par) < 1e-6 ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = QString::number(cur);
    cv->par = QString::number(par);
    compWerte.insert( key, cv );
}
void SC_MainGUI::compInt( QSpinBox *inp, QSettings &sets, QString key )
{
    int cur = inp->value();
    int par = sets.value(key,0).toInt();
    if ( cur == par ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = QString::number(cur);
    cv->par = QString::number(par);
    compWerte.insert( key, cv );
}
void SC_MainGUI::compInt( QComboBox *inp, QSettings &sets, QString key, QString prmt )
{
    int cur = inp->currentIndex();
    int par = sets.value(key,0).toInt();
    if ( cur == par ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = QString::number(cur);
    cv->par = QString::number(par);
    compWerte.insert( prmt+key, cv );
}
void SC_MainGUI::compBool( QGroupBox *tog, QSettings &sets, QString key, QString prmt )
{
    bool cur = tog->isChecked();
    bool par = sets.value(key,false).toBool();
    if ( cur == par ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = cur ? "True" : "False";
    cv->par = par ? "True" : "False";
    compWerte.insert( prmt+key, cv );
}
void SC_MainGUI::compBool( QCheckBox *tog, QSettings &sets, QString key, QString prmt )
{
    bool cur = tog->isChecked();
    bool par = sets.value(key,false).toBool();
    if ( cur == par ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = cur ? "True" : "False";
    cv->par = par ? "True" : "False";
    compWerte.insert( prmt+key, cv );
}
void SC_MainGUI::compBool( QRadioButton *tog, QSettings &sets, QString key, QString prmt )
{
    bool cur = tog->isChecked();
    bool par = sets.value(key,false).toBool();
    if ( cur == par ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = cur ? "True" : "False";
    cv->par = par ? "True" : "False";
    compWerte.insert( prmt+key, cv );
}

void SC_MainGUI::on_actionCompare_current_parameters_with_file_triggered()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString fn = data.value("LastParam",".").toString();
    fn = QFileDialog::getOpenFileName( this, "Compare Parameter for "+curMethod, fn,
                                      "Parameter (*.ini)",
                                      nullptr, QFileDialog::DontUseNativeDialog | QFileDialog::DontResolveSymlinks );
    if ( fn.isEmpty() ) return;
    //data.setValue("LastParam",fn);

    QSettings sets( fn, QSettings::IniFormat );
    sets.beginGroup( "Inputs" );
    QString m = sets.value( "CurMethod", "" ).toString();
    if ( ! m.isEmpty() )
    {
        if ( ui->tabMethods->tabText(ui->tabMethods->currentIndex()) != m )
        {
            _CompValues *cv = new _CompValues;
            cv->cur = ui->tabMethods->tabText(ui->tabMethods->currentIndex());
            cv->par = m;
            compWerte.insert( "CurMethod", cv );
        }
    }
    compBool( ui->radQ1, sets, "RadioButtonQ1" );
    compBool( ui->radQ2, sets, "RadioButtonQ2" );
    compBool( ui->radQ4, sets, "RadioButtonQ4" );
    compBool( ui->togExpandImage, sets, "ExpandImage" );
    compDouble( ui->inpAx1, sets, "EditAxis1x" );
    compDouble( ui->inpAy1, sets, "EditAxis1y" );
    compDouble( ui->inpAz1, sets, "EditAxis1z" );
    compDouble( ui->inpAx2, sets, "EditAxis2x" );
    compDouble( ui->inpAy2, sets, "EditAxis2y" );
    compDouble( ui->inpAz2, sets, "EditAxis2z" );
    compDouble( ui->inpAx3, sets, "EditAxis3x" );
    compDouble( ui->inpAy3, sets, "EditAxis3y" );
    compDouble( ui->inpAz3, sets, "EditAxis3z" );
    compDouble( ui->inpN1, sets, "Editxrel" );
    compDouble( ui->inpN2, sets, "Edityrel" );
    compDouble( ui->inpN3, sets, "Editzrel" );
    compDouble( ui->inpSigX, sets, "Editdom1" );
    compDouble( ui->inpSigY, sets, "Editdom2" );
    compDouble( ui->inpSigZ, sets, "Editdom3" );
    compDouble( ui->inpU1, sets, "Editx1rel" );
    compDouble( ui->inpU2, sets, "Edity1rel" );
    compDouble( ui->inpU3, sets, "Editz1rel" );
    compDouble( ui->inpV1, sets, "Editx2rel" );
    compDouble( ui->inpV2, sets, "Edity2rel" );
    compDouble( ui->inpV3, sets, "Editz2rel" );
    compInt( ui->inpHKLmax, sets, "Edithklmax" );
    compInt( ui->inpGridPoints, sets, "EditGridPoints" );
    compBool( ui->grpLattice3D, sets, "UseLattice3D" );
    compInt( ui->inpNumCores, sets, "Threads" );
    compDouble( ui->inpBCenterX, sets, "EditCenterX" );
    compDouble( ui->inpBCenterY, sets, "EditCenterY" );
    sets.endGroup();
    sets.beginGroup( "AI" );
    compString( ui->inpSubDir, sets, "LastSubDir", "AI:" );
    compInt( ui->cbsAIoutputFormat, sets, "Grayscale", "AI:" );
    compBool( ui->grpFileInput, sets, "FileInputEna", "AI:" );
    compString( ui->inpFileName, sets, "FileInputLast", "AI:" );
    compString( ui->inpFileClass, sets, "FileClass", "AI:" );
    compBool( ui->togAIUseFFTOutput, sets, "GenerateIFFT", "AI:" );
    compBool( ui->radAILinOutput, sets, "LinOut", "AI:" );
    //compBool( ui->radAILogOutput, ! sets.value("LinOut").toBool() );
    compBool( ui->togAIScaleOutput, sets, "ScaleOut", "AI:" );
    sets.endGroup();

    calcGui->compareParameter( sets, compWerte );

    QFile fout( fn+".diff.csv" );
    if ( ! fout.open(QIODevice::WriteOnly) )
        qDebug() << fout.errorString();
    else
        fout.write("Key ; Current ; Parameter from file" EOL);
    qDebug() << "Key" << "Cur" << "Par";
    QStringList slKeys = compWerte.keys();
    slKeys.sort();
    foreach ( QString k, slKeys )
    {
        qDebug() << k << compWerte[k]->cur << compWerte[k]->par;
        if ( fout.isOpen() ) fout.write( qPrintable(k+" ; "+compWerte[k]->cur+" ; "+compWerte[k]->par+EOL) );
    }
    qDebug() << "Vergleich zu Ende";
    if ( fout.isOpen() ) fout.close();

    QMessageBox::information( this, "Compare parameters", "The differences are saved to\n"+fout.fileName(), QMessageBox::Ok );
}


void SC_MainGUI::disableUnusedElements()
{
#define DISPARAM(x) x->hide()
//#define DISPARAM(x) x->setEnabled(false)
    //qDebug() << "DIS anf" << calcHelper::slDisWidgets;
    if ( calcHelper::slDisWidgets.removeOne("Uvec") )
    {
        DISPARAM(ui->inpU1);
        DISPARAM(ui->inpU2);
        DISPARAM(ui->inpU3);
        DISPARAM(ui->lblParamU);
    }
    if ( calcHelper::slDisWidgets.removeOne("Vvec") )
    {
        DISPARAM(ui->inpV1);
        DISPARAM(ui->inpV2);
        DISPARAM(ui->inpV3);
        DISPARAM(ui->lblParamV);
    }
    if ( calcHelper::slDisWidgets.removeOne("Nvec") )
    {
        DISPARAM(ui->inpN1);
        DISPARAM(ui->inpN2);
        DISPARAM(ui->inpN3);
        DISPARAM(ui->lblParamN);
    }
    while ( calcHelper::slDisWidgets.size() > 0 )
    {
        QStringList sl = calcHelper::slDisWidgets.first().split(":");
        if ( sl.size() == 2 )
        {
            paramHelper *phlp = calcGui->getParamPtr( sl[0], sl[1] );
            if ( phlp != nullptr )
            {
                //qDebug() << calcHelper::slDisWidgets.first() << phlp->key;
                if ( phlp->lbl1 != nullptr )
                    DISPARAM(phlp->lbl1);
                if ( phlp->gui.w != nullptr )
                    DISPARAM(phlp->gui.w);
            }
        }
        calcHelper::slDisWidgets.removeFirst();
    }
    //qDebug() << "DIS end" << calcHelper::slDisWidgets;
#undef DISPARAM//(x)
}


// Funktionen zum Extrahieren eines Bereiches als neues Image
// Idee: im widImage wird das Zeichnen eines Rechteckes mit der Maus aktiviert
//       und hier werden die Koordinaten dieses Rechtecks als Mittelpunkt und
//       halbe Breite gespeichert. Diese können dann auch manuell geändert
//       werden und werden im widImage entsprechend aktualisiert. Das Ausschneiden
//       selber wird hier erledigt.

void SC_MainGUI::on_grpExtractImage_toggled(bool arg1)
{
    if ( lastSelectedImage == nullptr ) return;
    if ( arg1 && ui->grpNoFitRegions->isChecked() )
    {
        ui->grpNoFitRegions->setChecked(false);
    }
    lastSelectedImage->enableExtract(arg1);
    ui->grpFrameColor->setEnabled( arg1 || ui->grpNoFitRegions->isChecked() );
    setCurrentExtractRect();
}

void SC_MainGUI::on_inpExtractCenterX_valueChanged(double )
{
    setCurrentExtractRect();
}

void SC_MainGUI::on_inpExtractCenterY_valueChanged(double )
{
    setCurrentExtractRect();
}

void SC_MainGUI::on_inpExtractGridSize_valueChanged(int )
{
    setCurrentExtractRect();
}

void SC_MainGUI::on_butDoExtract_clicked()
{
    if ( lastSelectedImage == nullptr ) return;

    // Datenbereich für das neue Image
    int dstX0 = 0;
    int dstX1 = ui->inpExtractGridSize->value()*2;
    int dstY0 = 0;
    int dstY1 = ui->inpExtractGridSize->value()*2;

    // Grenzen des RoI im alten Image
    int grfac = ui->cbsExtractScale->currentData().toInt();
    int grsiz = ui->inpExtractGridSize->value() * grfac;
    int srcX0 = lastSelectedImage->xmin() + ui->inpExtractCenterX->value() - grsiz;
    int srcX1 = lastSelectedImage->xmin() + ui->inpExtractCenterX->value() + grsiz;
    int srcY0 = lastSelectedImage->ymin() + ui->inpExtractCenterY->value() - grsiz;
    int srcY1 = lastSelectedImage->ymin() + ui->inpExtractCenterY->value() + grsiz;

    //qDebug() << "SRC" << srcX0 << srcX1 << srcY0 << srcY1;
    //qDebug() << "DST" << dstX0 << dstX1 << dstY0 << dstY1;
    //qDebug() << "IMG" << lastSelectedImage->debugInfo();

    double *dst = new double[(dstX1-dstX0)*(dstY1-dstY0)];
    double *d = dst;
    for ( int sy=srcY0; sy<srcY1; sy+=grfac )
        for ( int sx=srcX0; sx<srcX1; sx+=grfac )
        {
            double sum=0;
            int cnt=0;
            for ( int yy=0; yy<grfac; yy++ )
                for ( int xx=0; xx<grfac; xx++ )
                {
                    sum += lastSelectedImage->getData(sx+xx,sy+yy);
                    cnt++;
                }
            *(d++) = sum/cnt; // lastSelectedImage->getData(sx,sy);
        }

    widImage *tmp = addImage( true, dstX0, dstX1, dstY0, dstY1, dst,
                             "Extract ("+lastSelectedImage->windowTitle()+")", false );
    tmp->addMetaInfo( "Extract from", lastSelectedImage->getMetaTitle() );
    tmp->addMetaInfo( "NbRows", QString::number(dstY1-dstY0) );
    tmp->addMetaInfo( "NbCols", QString::number(dstX1-dstX0) );
    tmp->addMetaInfo( "Org center X", QString::number(ui->inpExtractCenterX->value()) );
    tmp->addMetaInfo( "Org center Y", QString::number(ui->inpExtractCenterY->value()) );
    tmp->addMetaInfo( "Extract scale", ui->cbsExtractScale->currentText() );

    delete[] dst;
}

void SC_MainGUI::on_cbsExtractFrameCol_activated(int index)
{
    if ( lastSelectedImage == nullptr ) return;
    lastSelectedImage->setExtractFrameCol( index );
}

void SC_MainGUI::on_cbsExtractScale_activated(int index)
{
    Q_UNUSED(index)
    setCurrentExtractRect();
}

void SC_MainGUI::extractUpdateRect( QRect rc )  // SLOT vom widImage
{
    if ( lastSelectedImage == nullptr ) return;

    if ( ui->grpExtractImage->isChecked() )
    {   // Mittelpunkt und GridPoints (=halbe Breite) für den Extract
        ui->inpExtractGridSize->blockSignals(true);
        ui->inpExtractCenterX->blockSignals(true);
        ui->inpExtractCenterY->blockSignals(true);

        ui->inpExtractCenterX->setValue( rc.center().x() );
        ui->inpExtractCenterY->setValue( rc.center().y() );
        ui->inpExtractGridSize->setValue( (qMax(rc.width(), rc.height()) / 2) / ui->cbsExtractScale->currentData().toInt() );

        ui->inpExtractGridSize->blockSignals(false);
        ui->inpExtractCenterX->blockSignals(false);
        ui->inpExtractCenterY->blockSignals(false);
    }
    else if ( ui->grpNoFitRegions->isChecked() && noFitCurRect > -1 )
    {   // Echte Koordinaten für die Ausschluss-Bereiche zum Fit
        QStringList vals;
        vals << QString::number(rc.left());
        vals << QString::number(rc.right());
        vals << QString::number(rc.top());
        vals << QString::number(rc.bottom());
        QTableWidgetItem *item;
        ui->tblNoFitRegions->blockSignals(true);
        for ( int c=0; c<4; c++ )
        {
            item = ui->tblNoFitRegions->item(noFitCurRect,c);
            if ( item != nullptr )
                item->setText(vals[c]);
            else
                ui->tblNoFitRegions->setItem(noFitCurRect,c,new QTableWidgetItem(vals[c]));
        }
        ui->tblNoFitRegions->blockSignals(false);
        lastSelectedImage->setNoFitRect( noFitCurRect, rc );
    }
}

void SC_MainGUI::setCurrentExtractRect()
{
    if ( lastSelectedImage == nullptr || !ui->grpExtractImage->isChecked() ) return;
    QRect rc;
    int grsiz = ui->inpExtractGridSize->value() * ui->cbsExtractScale->currentData().toInt();
    rc.setLeft( ui->inpExtractCenterX->value()   - grsiz );
    rc.setRight( ui->inpExtractCenterX->value()  + grsiz );
    rc.setTop( ui->inpExtractCenterY->value()    - grsiz );
    rc.setBottom( ui->inpExtractCenterY->value() + grsiz );
    lastSelectedImage->setExtractRect( rc );
}

void SC_MainGUI::setCurrentNoFitRect()
{
    if ( lastSelectedImage == nullptr || !ui->grpNoFitRegions->isChecked() ) return;
    QRect rc;
    QTableWidgetItem *item;
    for ( int r=0; r<4; r++ )
    {
        rc.setRect(0,0,0,0);
        item = ui->tblNoFitRegions->item(r,0);
        if ( item != nullptr ) rc.setLeft( item->text().toInt() );
        item = ui->tblNoFitRegions->item(r,1);
        if ( item != nullptr ) rc.setRight( item->text().toInt() );
        item = ui->tblNoFitRegions->item(r,2);
        if ( item != nullptr ) rc.setTop( item->text().toInt() );
        item = ui->tblNoFitRegions->item(r,3);
        if ( item != nullptr ) rc.setBottom( item->text().toInt() );
        //if ( rc.width() > 0 && rc.height() > 0 )
        lastSelectedImage->setNoFitRect( r, rc );   // Auch leere Rechtecke zum löschen
    }
}

void SC_MainGUI::on_grpNoFitRegions_toggled(bool arg1)
{
    if ( lastSelectedImage == nullptr ) return;
    if ( arg1 && ui->grpExtractImage->isChecked() )
    {
        ui->grpExtractImage->setChecked(false);
    }
    lastSelectedImage->enableNoFit(arg1);
    ui->grpFrameColor->setEnabled( arg1 || ui->grpExtractImage->isChecked() );
    setCurrentNoFitRect();
}

// Aktivieren einer Zelle
void SC_MainGUI::on_tblNoFitRegions_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn)
{
    //qDebug() << "NoFit curCelChg" << currentRow << currentColumn << previousRow << previousColumn;
    Q_UNUSED(currentColumn)
    Q_UNUSED(previousRow)
    Q_UNUSED(previousColumn)
    noFitCurRect = currentRow;
}

// Ändern des Wertes einer Zelle
void SC_MainGUI::on_tblNoFitRegions_cellChanged(int row, int column)
{
    qDebug() << "NoFit cellChg" << row << column;
}



void SC_MainGUI::on_butParamSearchPath_clicked()
{
    if ( bSearchParamsRunning ) return;
    QString dir = ui->inpParamSearchPath->text();
    dir = QFileDialog::getExistingDirectory( this, "Parameter file directory", dir );
    if ( dir.isEmpty() ) return;
    ui->inpParamSearchPath->setText(dir);
    ui->lisParamSearchResult->clear();
    ui->lisParamSearchResult->setEnabled(false);
    ui->butParamSearchGenerate->setEnabled(false);
    ui->lblParamSearchCount->setText(" ");
}

void SC_MainGUI::on_inpParamSearchPath_textEdited(const QString &arg1)
{
    if ( bSearchParamsRunning ) return;
    Q_UNUSED(arg1)
    if ( !ui->lisParamSearchResult->isEnabled() ) return;
    ui->lisParamSearchResult->clear();
    ui->lisParamSearchResult->setEnabled(false);
    ui->butParamSearchGenerate->setEnabled(false);
    ui->lblParamSearchCount->setText(" ");
}

void SC_MainGUI::searchParameterFilesHelper( QString d, int bl, QString msk )
{
    //qDebug() << "searchParameterFilesHelper" << d << bl;
    QDir dir(d);
    QFileInfoList fild = dir.entryInfoList( QStringList()<<"*", QDir::AllDirs | QDir::NoDotAndDotDot );
    QFileInfoList fili = dir.entryInfoList( QStringList()<<msk, QDir::Files );
    bool all = msk.endsWith("*");
    foreach ( QFileInfo fi, fili )
    {
        QString fn = fi.absoluteFilePath();
        if ( all )
        {
            if ( fn.endsWith(".ini") ) continue;
            if ( fn.endsWith(".png") ) continue;
            if ( fn.endsWith("ParamSearch.log") ) continue;
            if ( fn.contains("training table") ) continue;
            if ( fn.endsWith(".spr") )
            {   // Spezielle Abfrage, damit von den vielen *.spr Files nicht alle verwendet werden
                if ( fn.lastIndexOf('.') - fn.lastIndexOf('_') > 2 ) continue;
            }
        }
        else
        {
            if ( fn.endsWith("_out.ini") ) continue;
        }
        ui->lisParamSearchResult->addItem( fn.mid(bl) );
    }
    foreach ( QFileInfo fi, fild )
    {
        searchParameterFilesHelper( fi.absoluteFilePath(), bl, msk );
    }
}

void SC_MainGUI::on_butParamSearchDoit_clicked()
{
    if ( bSearchParamsRunning ) return;
    ui->lisParamSearchResult->setEnabled(true);
    ui->lisParamSearchResult->clear();
    int bl = ui->inpParamSearchPath->text().length();
    if ( !ui->inpParamSearchPath->text().endsWith("/") && !ui->inpParamSearchPath->text().endsWith("\\") ) bl++;
    qApp->setOverrideCursor(Qt::WaitCursor);
    searchParameterFilesHelper(ui->inpParamSearchPath->text(),bl,"*.ini");
    qApp->restoreOverrideCursor();
    ui->butParamSearchGenerate->setEnabled(false);
    ui->lblParamSearchCount->setText( QString::number(ui->lisParamSearchResult->count()) );
}

void SC_MainGUI::on_butDataSearchDoit_clicked()
{
    if ( bSearchParamsRunning ) return;
    ui->lisParamSearchResult->setEnabled(true);
    ui->lisParamSearchResult->clear();
    int bl = ui->inpParamSearchPath->text().length();
    if ( !ui->inpParamSearchPath->text().endsWith("/") && !ui->inpParamSearchPath->text().endsWith("\\") ) bl++;
    qApp->setOverrideCursor(Qt::WaitCursor);
    searchParameterFilesHelper(ui->inpParamSearchPath->text(),bl,"*.*");
    qApp->restoreOverrideCursor();
    ui->butParamSearchGenerate->setEnabled(false);
    ui->lblParamSearchCount->setText( QString::number(ui->lisParamSearchResult->count()) );
}

void SC_MainGUI::on_lisParamSearchResult_itemSelectionChanged()
{
    if ( bSearchParamsRunning ) return;
    bool ena = ui->lisParamSearchResult->count() > 0 &&
               ui->lisParamSearchResult->selectedItems().count() > 0;
    ui->butParamSearchGenerate->setEnabled( ena );
}

void SC_MainGUI::on_butParamSearchGenerate_clicked()
{
    ui->butParamSearchGenerate->setEnabled(false);
    QString base = ui->inpParamSearchPath->text();
    if ( !base.endsWith("/") ) base += "/";
    // Speichern der aktuellen Parameterwerte
    QString myLocPar = fnTempParamFile;
    myLocPar.replace(".ini","_psg.ini");
    performSaveParamOperation( myLocPar );
    // Logfile
    QFile flog(base+"ParamSearch.log");
    if ( !flog.open(QIODevice::WriteOnly) )
    {
        QMessageBox::critical( this, "Error", flog.fileName()+":"+flog.errorString(), QMessageBox::Ok );
        return;
    }
    bSearchParamsRunning = true;
    // Gang durch alle selektierten Files
    foreach ( QListWidgetItem *item, ui->lisParamSearchResult->selectedItems() )
    {
        QString fn = base + item->text();
        //item->setSelected(false);
        ui->lisParamSearchResult->setCurrentItem(item,QItemSelectionModel::Deselect);
        qApp->processEvents();
        //qDebug() << fn;
        if ( fn.endsWith(".ini") )
        {   // Parameterfiles
            if ( fn.endsWith("ConfigParams.ini") ) continue;
            flog.write( qPrintable("Read "+fn+EOL) );
            QString rv = local_Load_all_Parameters(fn,"");
            if ( !rv.isEmpty() ) flog.write( qPrintable("  "+rv) ); // in 'rv' wird immer ein EOL am Ende gesetzt.
            flog.flush();
            on_butCalc_clicked();
            if ( _bAbbruch )
            {
                flog.write("Aborted by user" EOL);
                break;
                // Debug: "C:/SimLab/sas-crystal/20220608 - 2D Fits und Fiber-Pattern/TestParamsCylinder.ini"
            }
            fn.replace(".ini",".png");
            flog.write( qPrintable(QString("  PrepTime %1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimPrep),12,'f',3)+EOL) );
            flog.write( qPrintable(QString("  CalcTime %1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimCalc),12,'f',3)+EOL) );
        }
        else
        {   // Datenfiles
            if ( ! local_OpenMeasFile( fn, &lastUsedImage ) ) continue;
            flog.write( qPrintable("Read "+fn+EOL) );
            if ( fn.endsWith(".spr") ) lastUsedImage->setLogScaling( false );
            fn.append(".png");
        }
        flog.write( qPrintable("Save "+fn+EOL+EOL) );
        lastUsedImage->saveImage(fn);
        lastUsedImage->close();
    }
    flog.close();
    bSearchParamsRunning = false;
    // Wiederherstellen der vorherigen Parameterwerte
    local_Load_all_Parameters(myLocPar,"");
}


//
// Definitions for the generation of simulation datasets
//


void SC_MainGUI::on_togTPVaddBS_toggled(bool checked)
{
    ui->inpTPVbeamx0->setEnabled( checked );
    ui->inpTPVbeamy0->setEnabled( checked );
    ui->inpTPVnbeam->setEnabled( checked );
    ui->inpTPVmbeam->setEnabled( checked );
}

void SC_MainGUI::on_butTPVoutPath_clicked()
{
    QString fn = QFileDialog::getExistingDirectory( this, "Save files in", ui->inpTPVoutPath->text() );
    if ( fn.isEmpty() ) return;
    ui->inpTPVoutPath->setText(fn);
}

void SC_MainGUI::on_butTPVsaveOptions_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("TPV");
    QString fn = data.value("LastSaveFile",".").toString();
    fn = QFileDialog::getSaveFileName( this, "Save TPV infos", fn, "Infofiles (*.sas_tpv)" );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".sas_tpv") ) fn += ".sas_tpv";
    QFile f(fn);
    if ( !f.open(QIODevice::WriteOnly) )
    {
        qDebug() << f.errorString();
        return;
    }
    f.close();
    data.setValue("LastSaveFile",fn);
    data.endGroup();
    QSettings sets(fn,QSettings::IniFormat);
    // Group [General]
    sets.setValue("ena_io",ui->togTPVio->isChecked());
    sets.setValue("val_io",ui->inpTPVio->value());    
    sets.setValue("ena_base",ui->togTPVbase->isChecked());
    sets.setValue("val_base",ui->inpTPVbase->value());
    sets.setValue("ena_radius",ui->togTPVradius->isChecked());
    sets.setValue("val_radius",ui->inpTPVradius->value());
    sets.setValue("ena_sigmar",ui->togTPVsigmar->isChecked());
    sets.setValue("val_sigmar",ui->inpTPVsigmar->value());
    sets.setValue("ena_length",ui->togTPVlength->isChecked());
    sets.setValue("val_length",ui->inpTPVlength->value());
    sets.setValue("ena_sigmal",ui->togTPVsigmal->isChecked());
    sets.setValue("val_sigmal",ui->inpTPVsigmal->value());
    sets.setValue("ena_phi",ui->togTPVphi->isChecked());
    //sets.setValue("val_phi",ui->inpTPVphi->value());
    sets.setValue("ena_a",ui->togTPVa->isChecked());
    sets.setValue("val_a",ui->inpTPVa->value());
    sets.setValue("ena_b",ui->togTPVb->isChecked());
    sets.setValue("val_b",ui->inpTPVb->value());
    sets.setValue("ena_c",ui->togTPVc->isChecked());
    sets.setValue("val_c",ui->inpTPVc->value());    
    sets.setValue("ena_rho",ui->togTPVrho->isChecked());
    sets.setValue("val_rho",ui->inpTPVrho->value());
    sets.setValue("ena_psi",ui->togTPVpsi->isChecked());
    sets.setValue("val_psi",ui->inpTPVpsi->value());
    sets.setValue("ena_dbeta",ui->togTPVdbeta->isChecked());
    sets.setValue("val_dbeta",ui->inpTPVdbeta->value());
    sets.setValue("ena_width",ui->togTPVwidth->isChecked());
    sets.setValue("val_width",ui->inpTPVwidth->value());
    sets.setValue("ena_phiwidth",ui->togTPVphiwidth->isChecked());
    sets.setValue("val_phiwidth",ui->inpTPVphiwidth->value());
    sets.setValue("ena_displacement",ui->togTPVdisplacement->isChecked());
    sets.setValue("val_displacement",ui->inpTPVdisplacement->value());
    sets.setValue("val_beamx0",ui->inpTPVbeamx0->value());
    sets.setValue("val_beamy0",ui->inpTPVbeamy0->value());
    sets.setValue("val_nbeam",ui->inpTPVnbeam->value());
    sets.setValue("val_mbeam",ui->inpTPVmbeam->value());
    sets.setValue("ena_addLines",ui->togTPVaddLines->isChecked());
    sets.setValue("val_addLinesH",ui->inpTPVaddLinesH->value());
    sets.setValue("val_addLinesV",ui->inpTPVaddLinesV->value());
    sets.setValue("ena_addBS",ui->togTPVaddBS->isChecked());
    sets.setValue("ena_addNoise",ui->togTPVaddNoise->isChecked());
    sets.setValue("ena_convolute",ui->togTPVconvolute->isChecked());
    sets.setValue("ena_calcRphi",ui->togTPVcalcRphi->isChecked());
    sets.setValue("ena_calcFFT",ui->togTPVcalcFFT->isChecked());
    sets.setValue("val_numimg",ui->inpTPVnumimg->value());
    sets.setValue("val_outPath",ui->inpTPVoutPath->text());
    sets.sync();
    // In dieses File kommen jetzt auch alle Daten
    // Groups [FCC%20Spheres] und [Inputs]
    performSaveParamOperation( fn );
    // FFT Tab
    sets.beginGroup("FFT");
    sets.setValue( "FFTLinInput", ui->radFFTLinInput->isChecked() );
    sets.setValue( "FFTScaleRphi", ui->togFFTscaleRphi->isChecked() );
    sets.setValue( "FFTclipRphi", ui->togFFTclipRphi->isChecked() );
    sets.setValue( "FFTclip40Rphi", ui->togFFTclip40Rphi->isChecked() );
    sets.setValue( "FFTsizeRphi", ui->cbsFFTsizeRphi->currentIndex() );
    sets.setValue( "DispRphi", ui->togFFTdispRphi->isChecked() );
    sets.setValue( "FFTuseRphi", ui->grpFFTuseRphi->isChecked() );
    sets.setValue( "FFTScaleOutput", ui->togFFTscaleOut->isChecked() );
    sets.setValue( "FFTclipOutput", ui->togFFTclipOut->isChecked() );
    sets.setValue( "FFTclip40Output", ui->togFFTclip40Out->isChecked() );
    sets.setValue( "FFTSwapOutput", ui->togIFFTSwap->isChecked() );
    sets.setValue( "FFTsizeOut", ui->cbsFFTsizeOut->currentIndex() );
    if ( ui->radFFToutReal->isChecked() )
        sets.setValue( "FFToutput", 0 );
    else if ( ui->radFFToutImag->isChecked() )
        sets.setValue( "FFToutput", 1 );
    else if ( ui->radFFToutBetrag->isChecked() )
        sets.setValue( "FFToutput", 2 );
    else if ( ui->radFFToutSpectrum->isChecked() )
        sets.setValue( "FFToutput", 3 );
    sets.endGroup();
}

void SC_MainGUI::on_butTPVloadOptions_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("TPV");
    QString fn = data.value("LastSaveFile",".").toString();
    fn = QFileDialog::getOpenFileName( this, "Open TPV infos", fn, "Infofiles (*.sas_tpv)" );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".sas_tpv") ) fn += ".sas_tpv";
    data.setValue("LastSaveFile",fn);
    data.endGroup();
    QSettings sets(fn,QSettings::IniFormat);
    ui->togTPVio->setChecked( sets.value("ena_io",true).toBool() );
    ui->inpTPVio->setValue( sets.value("val_io",0.5).toDouble() );
    ui->inpTPVio->setEnabled(ui->togTPVio->isChecked());    
    ui->togTPVbase->setChecked( sets.value("ena_base",true).toBool() );
    ui->inpTPVbase->setValue( sets.value("val_base",0.5).toDouble() );
    ui->inpTPVbase->setEnabled(ui->togTPVbase->isChecked());
    ui->togTPVradius->setChecked( sets.value("ena_radius",true).toBool() );
    ui->inpTPVradius->setValue( sets.value("val_radius",0.5).toDouble() );
    ui->inpTPVradius->setEnabled(ui->togTPVradius->isChecked());
    ui->togTPVsigmar->setChecked( sets.value("ena_sigmar",true).toBool() );
    ui->inpTPVsigmar->setValue( sets.value("val_sigmar",0.5).toDouble() );
    ui->inpTPVsigmar->setEnabled(ui->togTPVsigmar->isChecked());
    ui->togTPVlength->setChecked( sets.value("ena_length",true).toBool() );
    ui->inpTPVlength->setValue( sets.value("val_length",0.5).toDouble() );
    ui->inpTPVlength->setEnabled(ui->togTPVlength->isChecked());
    ui->togTPVsigmal->setChecked( sets.value("ena_sigmal",true).toBool() );
    ui->inpTPVsigmal->setValue( sets.value("val_sigmal",0.5).toDouble() );
    ui->inpTPVsigmal->setEnabled(ui->togTPVsigmal->isChecked());
    ui->togTPVphi->setChecked( sets.value("ena_phi",true).toBool() );
    //ui->inpTPVphi->setValue( sets.value("val_phi",0.5).toDouble() );
    ui->togTPVa->setChecked( sets.value("ena_a",true).toBool() );
    ui->inpTPVa->setValue( sets.value("val_a",0.5).toDouble() );
    ui->inpTPVa->setEnabled(ui->togTPVa->isChecked());
    ui->togTPVb->setChecked( sets.value("ena_b",true).toBool() );
    ui->inpTPVb->setValue( sets.value("val_b",0.5).toDouble() );
    ui->inpTPVb->setEnabled(ui->togTPVb->isChecked());
    ui->togTPVc->setChecked( sets.value("ena_c",true).toBool() );
    ui->inpTPVc->setValue( sets.value("val_c",0.5).toDouble() );
    ui->inpTPVc->setEnabled(ui->togTPVc->isChecked());    
    ui->togTPVrho->setChecked( sets.value("ena_rho",true).toBool() );
    ui->inpTPVrho->setValue( sets.value("val_rho",0.5).toDouble() );
    ui->inpTPVrho->setEnabled(ui->togTPVrho->isChecked());
    ui->togTPVpsi->setChecked( sets.value("ena_psi",true).toBool() );
    ui->inpTPVpsi->setValue( sets.value("val_psi",0.5).toDouble() );
    ui->inpTPVpsi->setEnabled(ui->togTPVpsi->isChecked());    
    ui->togTPVdbeta->setChecked( sets.value("ena_dbeta",true).toBool() );
    ui->inpTPVdbeta->setValue( sets.value("val_dbeta",0.5).toDouble() );
    ui->inpTPVdbeta->setEnabled(ui->togTPVdbeta->isChecked());
    ui->togTPVwidth->setChecked( sets.value("ena_width",true).toBool() );
    ui->inpTPVwidth->setValue( sets.value("val_width",0.5).toDouble() );
    ui->inpTPVwidth->setEnabled(ui->togTPVwidth->isChecked());
    ui->togTPVphiwidth->setChecked( sets.value("ena_phiwidth",true).toBool() );
    ui->inpTPVphiwidth->setValue( sets.value("val_phiwidth",0.5).toDouble() );
    ui->inpTPVphiwidth->setEnabled(ui->togTPVphiwidth->isChecked());
    ui->togTPVdisplacement->setChecked( sets.value("ena_displacement",true).toBool() );
    ui->inpTPVdisplacement->setValue( sets.value("val_displacement",0.5).toDouble() );
    ui->inpTPVdisplacement->setEnabled(ui->togTPVdisplacement->isChecked());
    ui->togTPVaddBS->setChecked( sets.value("ena_addBS",true).toBool() );
    ui->inpTPVbeamx0->setValue( sets.value("val_beamx0",0.5).toDouble() );
    ui->inpTPVbeamx0->setEnabled(ui->togTPVaddBS->isChecked());
    ui->inpTPVbeamy0->setValue( sets.value("val_beamy0",0.5).toDouble() );
    ui->inpTPVbeamy0->setEnabled(ui->togTPVaddBS->isChecked());
    ui->inpTPVnbeam->setValue( sets.value("val_nbeam",0.5).toDouble() );
    ui->inpTPVnbeam->setEnabled(ui->togTPVaddBS->isChecked());
    ui->inpTPVmbeam->setValue( sets.value("val_mbeam",0.5).toDouble() );
    ui->inpTPVmbeam->setEnabled(ui->togTPVaddBS->isChecked());
    ui->togTPVaddLines->setChecked( sets.value("ena_addLines",false).toBool() );
    ui->inpTPVaddLinesH->setValue( sets.value("val_addLinesH",1).toInt() );
    ui->inpTPVaddLinesV->setValue( sets.value("val_addLinesV",5).toInt() );
    ui->inpTPVaddLinesH->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesV->setEnabled(ui->togTPVaddLines->isChecked());
    ui->togTPVaddNoise->setChecked( sets.value("ena_addNoise",true).toBool() );
    ui->togTPVconvolute->setChecked( sets.value("ena_convolute",true).toBool() );
    ui->togTPVcalcRphi->setChecked( sets.value("ena_calcRphi",true).toBool() );
    ui->togTPVcalcFFT->setChecked( sets.value("ena_calcFFT",true).toBool() );
    ui->inpTPVnumimg->setValue( sets.value("val_numimg",10).toInt() );
    ui->inpTPVoutPath->setText( sets.value("val_outPath",".").toString() );
    local_Load_all_Parameters(fn,"");
}

void SC_MainGUI::on_butTPVstart_clicked()
{
    performStartAI_TPV( true );
}


//typedef QHash<QString/*Variable*/, double/*value*/> _loopVariable;
//typedef QMultiHash<QString/*CLASS*/, _loopVariables > _loopDefinition;
SC_MainGUI::_loopDefinition SC_MainGUI::getTPVLoopDefinition()
{
    _loopDefinition retval;  // MultiHash für alle Rechenschritte
    double tmp;
    //for ( int iNum=0; iNum<ui->inpTPVnumimg->value(); iNum++ )
    {
        calcval.clear();
        calcval.insert( "#Calculation#", 0 );
        calcval.insert( "#FACTOR#", ui->inpTPVnumimg->value() );
#define ADDCALCVAL(tog,inp,nam) if ( tog->isChecked() ) calcval.insert( nam, (tmp=inp->value()) )
        ADDCALCVAL( ui->togTPVio,           ui->inpTPVio,           "I0"              );
        ADDCALCVAL( ui->togTPVbase,         ui->inpTPVbase,         "Base"            );
        ADDCALCVAL( ui->togTPVradius,       ui->inpTPVradius,       "EditRadius"      );
        ADDCALCVAL( ui->togTPVsigmar,       ui->inpTPVsigmar,       "EditSigma"       );
        ADDCALCVAL( ui->togTPVlength,       ui->inpTPVlength,       "Length"          );
        ADDCALCVAL( ui->togTPVsigmal,       ui->inpTPVsigmal,       "SigmaL"          );
        ADDCALCVAL( ui->togTPVa,            ui->inpTPVa,            "uca"             );
        if ( ui->togTPVb->isChecked() ) calcval.insert( "ucb", tmp ); // =uca
        if ( ui->togTPVc->isChecked() ) calcval.insert( "ucc", tmp ); // =uca
        ADDCALCVAL( ui->togTPVrho,          ui->inpTPVrho,          "EditRho"         );
        ADDCALCVAL( ui->togTPVpsi,          ui->inpTPVpsi,          "ucpsi"           );
        ADDCALCVAL( ui->togTPVdbeta,        ui->inpTPVdbeta,        "EditDbeta"       );
        ADDCALCVAL( ui->togTPVwidth,        ui->inpTPVwidth,        "EditDomainSize"  ); // =4/EditDomainSize
        ADDCALCVAL( ui->togTPVphiwidth,     ui->inpTPVphiwidth,     "EditAzi"         ); // =4/EditAzi
        ADDCALCVAL( ui->togTPVdisplacement, ui->inpTPVdisplacement, "EditDebyeWaller" ); // =EditDebyeWaller*EditDebyeWaller/3
        //ADDCALCVAL( ui->togTPVphi,          ui->inpTPVphi,          "phi"             );
        if ( ui->togTPVphi->isChecked() ) calcval.insert( "phi", 1 );
        //retval.insert( QString("%1").arg(iNum+1,6,10,QChar('0')), calcval );  -> scat_000001.spr
        //retval.insert( QString("%1").arg(iNum), calcval ); // -> scat_0.spr
        retval.insert( "0", calcval ); // Faktor in der Zeile
    }
    return retval;
}


void SC_MainGUI::on_togFFTscaleRphi_toggled(bool checked)
{
    if ( checked ) ui->togFFTclipRphi->setChecked(false);
}

void SC_MainGUI::on_togFFTclipRphi_toggled(bool checked)
{
    if ( checked ) ui->togFFTscaleRphi->setChecked(false);
}

void SC_MainGUI::on_togFFTscaleOut_toggled(bool checked)
{
    if ( checked ) ui->togFFTclipOut->setChecked(false);
}

void SC_MainGUI::on_togFFTclipOut_toggled(bool checked)
{
    if ( checked ) ui->togFFTscaleOut->setChecked(false);
}


void SC_MainGUI::on_butFFTverifyRphi_clicked()
{
    QString dn = QFileDialog::getExistingDirectory(this, "Save Verification files",
                                                   "C:/SimLab/sas-crystal/20221004/VerifyRphi" );
    if ( dn.isEmpty() ) return;

    QDir destdir(dn);

    QPixmap pix   = QPixmap(256,256);
    double data[256][256];

    QPainter pnt(&pix);
    QPen pen( Qt::white, 4 );
    pnt.setPen( pen );
    for ( float ang=0.; ang<360.; ang+=15. )
    {
        pix.fill(Qt::black);
        QLineF line( QPointF(128,128), QPointF(128,256) );
        line.setAngle(ang);
        pnt.drawLine(line);

        pix.save( destdir.absoluteFilePath(QString("%1_pix.png").arg(ang)) );

        QImage img = pix.toImage();
        for ( int y=0; y<256; y++ )
            for ( int x=0; x<256; x++ )
                data[y][x] = img.pixel(x,y);

        //widImage *wid = addImage( true, 0, 256, 0, 256, &data[0][0], "Linie", false );
        //wid->saveImage( destdir.absoluteFilePath(QString("%1_wid.png").arg(ang)) );
        //wid->close();
        //wid->deleteLater();

        double *outr = SasCalc_PostProc::inst()->generateRphi(0, 256, 0, 256, 128, 128, 256,
                                                              &data[0][0],
                                                              true,     // scale
                                                              false,    // clip
                                                              false );  // clip40

        widImage *wid = addImage( true, 0, 256, 0, 256, outr, "R,phi", false );
        wid->saveImage( destdir.absoluteFilePath(QString("%1_rphi.png").arg(ang)) );
        wid->close();
        wid->deleteLater();

    }
    qDebug() << "Done";
}


/*
"training table.spr":
    25    10 Start pixel = (       1       1)
  intensity radius length sigma_r sigma_l rho a_cell psi phi width aziwidth dw base x0 y0 nbeam mbeam linex liney points holder hold.ang. bs.shape det.shape
  8.01E+03 1.14E+01 1.04E+01 8.50E-02 1.34E-01 1.29E-01 4.33E+01 0.00E+00 2.80E+01 1.03E-02 3.28E-02 8.33E-01 1.72E-10 4.26E+01 4.45E+01 6.00E+00 5.00E+00 -1.46E+02 1.22E+02 6.40E+01 true 2.64E+02 circ sq

"training table.txt":
    28    1000 Start pixel = (       1       1)
  run_no.   intensity radius    radius_m  length    sigma_r   sigma_l   rho       a_cell    b_cell    c_cell   psi      phi      width    aziwidth dw       dbeta    base     x0       y0       nbeam    mbeam    linex    liney    points   holder   hold_ang bs_shape det_shape
  0 1.61E+04 1.57E+01 2.68E+00 2.86E+01 6.39E-02 1.10E-01 8.40E-02 4.77E+01 4.77E+01 2.62E+01 0.00E+00 8.10E+01 1.44E-02 1.11E-02 1.00E+00 4.26E-01 9.34E-03 9.22E+01 5.95E+01 3.00E+00 7.00E+00 1.38E+02 -6.40E+01 6.40E+01 false 1.46E+02 rect sq
*/
void SC_MainGUI::on_butReadTrainingTable_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("TPV");
    QString fn = data.value("LastTrainTable",".").toString();
    fn = QFileDialog::getOpenFileName( this, "Open Training Table", fn, "training table (*.spr *.txt)" );
    if ( fn.isEmpty() ) return;
    data.setValue("LastTrainTable",fn);
    data.endGroup();

    QFile fin(fn);
    if ( ! fin.open(QIODevice::ReadOnly) )
    {
        qDebug() << fn << fin.errorString();
        return;
    }

    QString line;
    /*line =*/ fin.readLine();  // 1. Zeile, die Werte dort brauche ich nicht
    line = fin.readLine();  // 2. Zeile: Schlüsselworte
    QStringList keys = line.trimmed().split(" ",Qt::SkipEmptyParts);
    qDebug() << keys << keys.size();
    // Umsetzen der Schlüsselworte in die Kennungen der Variablen, die hier verwendet werden
    // <schlüssel in datei> | <kennung hier>     Normale Datenwerte zum Übernehmen
    // <schlüssel in datei> | -                  Hier nicht verwendet
    // <schlüssel in datei> | @<kennung hier>    Globale Werte, nicht über den Namen setzbar
    static QStringList keys2kenn = {
        "intensity|I0",
        "radius|EditRadius",
        "length|Length",
        "sigma_r|EditSigma",
        "sigma_l|SigmaL",
        "rho|EditRho",
        "a_cell|uca",
        "psi|ucpsi",
        "phi|phi",
        "width|-",
        "aziwidth|EditAzi",
        "dw|EditDebyeWaller",
        "base|Base",
        "x0|@BSX",
        "y0|@BSY",
        "nbeam|-",
        "mbeam|-",
        "linex|-",
        "liney|-",
        "points|@GP",
        "holder|-",
        "hold.ang.|-",
        "bs.shape|-",
        "det.shape|-",
        "run_no.|-",
        "radius_m|EditRadiusi",
        "b_cell|ucb",
        "c_cell|ucc",
        "dbeta|EditDbeta",
        "hold_ang|-",
        "bs_shape|-",
        "det_shape|-",
        "pixx|EditPixelX",
        "pixy|EditPixelY"
    };
    QStringList usedKeys;
    tpvParamsKeys.clear();
    for ( int i=0; i<keys.size(); i++ )
    {
        bool found = false;
        for ( int k=0; k<keys2kenn.size(); k++ )
            if ( keys2kenn[k].startsWith(keys[i]+"|") )
            {
                usedKeys << keys2kenn[k].mid(keys2kenn[k].indexOf("|")+1);
                if ( usedKeys.last()[0] != '-' && usedKeys.last()[0] != '@' )
                    tpvParamsKeys << usedKeys.last();
                found = true;
                break;
            }
        if ( !found )
            usedKeys << keys[i];
    }
    qDebug() << usedKeys;
    ui->actionFind_parameters_changing_image->setEnabled( tpvParamsKeys.size() > 0 );
    bool pixelSizeChanged = tpvParamsKeys.contains("EditPixelX");
    double defPixX = calcGui->currentParamValueDbl("","EditPixelX");
    double defPixY = calcGui->currentParamValueDbl("","EditPixelY");
    int zmaxt = ui->inpGridPoints->value();
#define SETCOL(w,c) {QPalette pal=w->palette(); pal.setBrush(QPalette::Base,c); w->setPalette(pal);}
    while ( !fin.atEnd() )
    {
        line = fin.readLine();
        if ( line.isEmpty() ) continue;
        QStringList svals = line.trimmed().split(" ",Qt::SkipEmptyParts);
        if ( svals.size() < keys.size() ) continue;     // Weniger Werte als Schlüssel werden ignoiert
        for ( int i=0; i<keys.size(); i++ )
        {
            double d = svals[i].toDouble();
            if ( usedKeys[i][0] == '@' )
            {
                if ( usedKeys[i] == "@BSX" )
                {
                    ui->inpBCenterX->setValue( d );
                    SETCOL( ui->inpBCenterX, Qt::yellow );
                }
                else if ( usedKeys[i] == "@BSY" )
                {
                    ui->inpBCenterY->setValue( d );
                    SETCOL( ui->inpBCenterY, Qt::yellow );
                }
                else if ( usedKeys[i] == "@GP" )
                {
                    ui->inpGridPoints->setValue( d );
                    SETCOL( ui->inpGridPoints, Qt::yellow );
                }
            }
            else if ( usedKeys[i][0] != '-' )
                calcGui->updateParamValue( "", usedKeys[i], d, Qt::yellow, false );
        }
        if ( ui->inpBCenterX->value() > 0 || ui->inpBCenterY->value() > 0 )
        {
            ui->inpBCenterX->setValue( ui->inpBCenterX->value() - ui->inpGridPoints->value() );
            ui->inpBCenterY->setValue( ui->inpBCenterY->value() - ui->inpGridPoints->value() );
            SETCOL( ui->inpBCenterX, Qt::yellow );
            SETCOL( ui->inpBCenterY, Qt::yellow );
        }
        // Es wird nur die erste Zeile mit Datenwerten genutzt.
        break;
#undef SETCOL//(w,c)
    }
    fin.close();

    if ( !pixelSizeChanged )
    {
        int zmax = ui->inpGridPoints->value();

        calcGui->updateParamValue( "", "EditPixelNoX", 2*zmax+1, Qt::yellow );
        calcGui->updateParamValue( "", "EditPixelNoY", 2*zmax+1, Qt::yellow );

        calcGui->updateParamValue( "", "EditPixelX", defPixX*(2*zmaxt+1)/(2*zmax+1), Qt::yellow );
        calcGui->updateParamValue( "", "EditPixelY", defPixY*(2*zmaxt+1)/(2*zmax+1), Qt::yellow );
    }

    // Jetzt wird noch der erste "scat" Datensatz aus dem verzeichnis mit der tabelle geladen
    QDir d(QFileInfo(fn).path());
    qDebug() << fn;
    d.cd("scat");
    qDebug() << d.absoluteFilePath("scat_0.spr");
    if ( ! local_OpenMeasFile(d.absoluteFilePath("scat_0.spr"),nullptr) )
        qDebug() << "Open failed.";
}

void SC_MainGUI::on_actionFind_parameters_changing_image_triggered()
{
    //qDebug() << "TPV" << tpvParamsKeys;
    QStringList all = calcGui->paramsForMethod(0,true,false,false);

    foreach ( QString s, tpvParamsKeys )
        all.removeOne(s);
    all.removeOne("RadioButtonQ1");
    all.removeOne("RadioButtonQ2");
    all.removeOne("RadioButtonQ4");
    all.removeOne("ExpandImage");
    all.removeOne("Edithklmax");
    all.removeOne("EditGridPoints");
    all.removeOne("LATTrows");
    all.removeOne("LATTcols");
    all.removeOne("LATTceny");
    all.removeOne("LATTcenx");
    all.removeOne("BeamPosX");
    all.removeOne("BeamPosY");

    all.removeOne("EditDet");       // Detektor-Abstand lassen wir hier mal weg
    all.removeOne("EditPixelX");    // Pixelsizes lassen wir hier auch weg
    all.removeOne("EditPixelY");

    // Jetzt sind in 'all' nur noch Parameter enthalten, die nicht vom obigen Lesen gesetzt wurden.
    // Diese werden jetzt durchgespielt, ob diese den Datensatz ändern.

    // Zum Test schauen, ob der Datensatz nicht zu groß wird, das würde die Rechenzeit bei dem
    //  Durchgehen der Parameter unnötig erhöhen.
    if ( ui->inpGridPoints->value() > 100 )
    {
        ui->inpGridPoints->setValue(64);
        // Dann aber auch den Beamstop mittig setzen
        ui->inpBCenterX->setValue(0);
        ui->inpBCenterY->setValue(0);
    }

    // Erzeugen des aktuellen Datensatzes
    on_butCalc_clicked();
    int mx = ui->inpGridPoints->value();
    if ( ui->radQ1->isChecked() )
        mx = 4*mx*mx;
    else if ( ui->radQ2->isChecked() )
        mx = 2*mx*mx;
    else if ( ui->radQ4->isChecked() )
        mx = mx*mx;
    double *check = new double[mx];
    memcpy( check, calcGui->data(), mx*sizeof(double) );

    // Variationen über alle Parameter
    QStringList modParams;
    foreach ( QString s, all )
    {
        double curval = calcGui->currentParamValueDbl("",s);
        double min, max;
        bool count; // uninteressant
        if ( ! calcGui->limitsOfParamValue( "", s, min, max, count ) )
        {
            qDebug() << "***" << s;
            continue;
        }
        double v1, v2;

        if ( curval == 0 )
        {
            if ( min == -10000 && max == +10000 )
            {   // Default
                v1 = 1;
                v2 = 2;
            }
            else
            {   // Spezialisiert
                v1 = (max-min) * 0.05;
                v2 = (max-min) * 0.10;
            }
        }
        else
        {
            if ( min == -10000 && max == +10000 )
            {   // Default
                v1 = curval * 0.95;
                v2 = curval * 1.05;
            }
            else
            {   // Spezialisiert
                v1 = curval * 0.95;
                if ( v1 < min ) v1 = min;
                v2 = curval * 1.05;
                if ( v2 > max ) v2 = max;
            }
        }

        //qDebug() << s << curval << min << max << "/" << v1 << v2;

        if ( vergleicheCheckData( s, v1, check, mx ) )
        {
            //qDebug() << "TRUE-1";
            //break;
            modParams << s;
        }
        if ( _bAbbruch ) break;
        if ( vergleicheCheckData( s, v2, check, mx ) )
        {
            //qDebug() << "TRUE-2";
            //break;
            if ( !modParams.contains(s) ) modParams << s;
        }
        calcGui->updateParamValue( "", s, curval, Qt::black );
        if ( _bAbbruch ) break;
    }
    qDebug() << modParams;
}

bool SC_MainGUI::vergleicheCheckData( QString par, double val, const double *check, int mx )
{
    calcGui->updateParamValue( "", par, val, Qt::black );

    on_butCalc_clicked();
    if ( _bAbbruch ) return false; // Keine Änderungsmarkierung bei Abbruch

    const double *tmp = calcGui->data();

    for ( int i=0; i<mx; i++ )
        if ( fabs( *(tmp++) - *(check++) ) > eps4 )
            return true;

    return false;
}


void SC_MainGUI::on_butResetColorMarker_clicked()
{
    calcGui->resetParamColorMarker( Qt::white );
    QPalette pal;
#define SETCOL(w,c) pal=w->palette(); pal.setBrush(QPalette::Base,c); w->setPalette(pal);
    SETCOL( ui->inpBCenterX, Qt::white );
    SETCOL( ui->inpBCenterY, Qt::white );
    SETCOL( ui->inpGridPoints, Qt::white );
#undef SETCOL//(w,c)
}
