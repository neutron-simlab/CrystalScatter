/**
 * Program: sas_scatter2
 *  Author: Michael Wagener, FZJ, JCNS-1
 *          Stephan Förster, FZJ, JCNS-1 wrote the original Scatter in Pascal.
 *    Date: March 2023
 *    File: sc_maingui.cpp
 *          Contains all routines for the main GUI compnents.
 */


#include "sc_maingui.h"
#include "ui_sc_maingui.h"
#include "sc_postproc.h"
#include "dlgtimetests.h"
#include "dlgconfigautofit.h"
#include "myguiparam.h"
#include <QDebug>
#include <QHeaderView>
#include <QFileDialog>
#include <QSettings>
#include <QMessageBox>
#include <QDesktopWidget>
#include <QDate>
#include <QPainter>
#include <thread>
#ifdef Q_OS_LINUX
#include <QStyleFactory>
#endif


//#define NOCBSCALLBACK
// Wenn definiert, werden keine Anpassungen der Labels bei den CBS-Callbacks durchgeführt. Somit
//  bleiben die Labels so, wie die Parameternamen lauten. Hilfreich für erklärende Screenshots.


//#define ChatbotIgnoreImages
// Wenn definiert, werden die im Chatbot-Tab generierten Images nicht automatisch angezeigt.
// Es werden ja nur PNG-Images erzeugt und keine normalen Datenfiles gesichert, die im Programm
//  weiterverwendet werden könnten.



// Global (static) variables
SC_MainGUI *SC_MainGUI::current = nullptr;
bool SC_MainGUI::_bAbbruch = false;
bool SC_MainGUI::_bAbortTimer = false;


/**
 * @brief SC_MainGUI::SC_MainGUI
 * @param parent - Qt parent object, normally null
 * Construktor for the main GUI application.
 */
SC_MainGUI::SC_MainGUI(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::SC_MainGUI)
{
    QSettings sets(SETT_APP,SETT_GUI);

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
    fitparams = nullptr;
    timingTestsRunning = false;

    ui->setupUi(this);
    current = this;

#ifdef Q_OS_LINUX
    // To draw frames around the QGrupBox not visible in the default linux style
    QStyle *winstyle = QStyleFactory::create("Windows");
    ui->grpDataFile->setStyle(winstyle);
    ui->grpLattice->setStyle(winstyle);
    ui->grpParticle->setStyle(winstyle);
    ui->grpPeakShape->setStyle(winstyle);
    ui->grpAniso->setStyle(winstyle);
    ui->grpOrientation->setStyle(winstyle);
    ui->grpInformations->setStyle(winstyle);
    ui->grpColorMarker->setStyle(winstyle);
    ui->grpControl->setStyle(winstyle);
    ui->grpPixelManipul->setStyle(winstyle);
    ui->grpCalculation->setStyle(winstyle);
    ui->grpCalcDestCalc->setStyle(winstyle);
    ui->grpEditQmaxGroup->setStyle(winstyle);
    ui->grpQuadrants->setStyle(winstyle);
    ui->grpExperimant->setStyle(winstyle);
    ui->grpWindows->setStyle(winstyle);
    ui->grpExtractImage->setStyle(winstyle);
    ui->grpNoFitRegions->setStyle(winstyle);
    ui->grpFrameColor->setStyle(winstyle);
    ui->grpFFToutput->setStyle(winstyle);
    ui->grpFFTuseRphi->setStyle(winstyle);
    ui->grpFFTinput->setStyle(winstyle);
    ui->grpFitParams->setStyle(winstyle);
    ui->grpFitResult->setStyle(winstyle);
    ui->grpFitStart->setStyle(winstyle);
    ui->grpTPVvariation->setStyle(winstyle);
    ui->grpTPVoutPath->setStyle(winstyle);
    ui->grpTPVothers->setStyle(winstyle);
    ui->grpVariables->setStyle(winstyle);
    ui->grpAIoptions->setStyle(winstyle);
    ui->grpFileInput->setStyle(winstyle);
    ui->grpClass->setStyle(winstyle);
    ui->grpAIoutdir->setStyle(winstyle);
    ui->grpCalcDestConf->setStyle(winstyle);
    ui->grpLastFitActions->setStyle(winstyle);
    ui->grpTimingInfo->setStyle(winstyle);
    ui->grpDefColTbl->setStyle(winstyle);
    ui->grpValChg->setStyle(winstyle);
    ui->grpCalcTests->setStyle(winstyle);
    ui->grpScanFiles->setStyle(winstyle);
#endif

    // Set menu flags to defaults
    ui->actionHide_unused_values->setChecked(true);  // true: hide unused values, false: show unused values with grey labels

    // Small text in the upper rigth corner with the current version string
    QLabel *lblVersInfo = new QLabel("Version " MYVERSION);
    QPalette palv = lblVersInfo->palette();
    palv.setColor( QPalette::WindowText, Qt::gray );
    lblVersInfo->setPalette(palv);
    lblVersInfo->setToolTip("(C) Michael Wagener, Stephan Förster, JCNS-1\nForschungszentrum Jülich GmbH");
    ui->tabMain->setCornerWidget(lblVersInfo);

    // Get the data directory path (defaults to 'data' near the program executable)
    QDir fileDir( qApp->applicationDirPath() );
    while ( fileDir.absolutePath().contains("debug", Qt::CaseInsensitive) ||
            fileDir.absolutePath().contains("release", Qt::CaseInsensitive) )
        fileDir.cdUp();
    fileDir.cd( "data" );
    dataPath = fileDir.absolutePath();
    lastDataFile = "";
    numberCounter = 0;
    closeMainActive = false;

    // Set the filename for the automatic parameter saving before the calculation starts
    fnTempParamFile = fileDir.absoluteFilePath("TempParamSave.ini");

    restoreGeometry( sets.value("GUIgeometry").toByteArray() );

    // --------------------------
    // --- TAB  Configuration ---
    // --------------------------
    // Find the number of threads possible
    int nthreads = static_cast<int>(std::thread::hardware_concurrency());
    ui->inpNumCores->setMaximum( nthreads );
    // Minimum settings later (0: with GPU, 1: without GPU)
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
    foreach ( QString s, sl )
        if ( s.startsWith("Residuen") )
            sl.removeOne(s);
    ui->cbsDefaultColTbl->addItems( sl );
    ui->cbsDefaultColTbl->setCurrentIndex( sets.value("DefColTbl",1).toInt() );
    on_cbsDefaultColTbl_activated( sets.value("DefColTbl",1).toInt() );

    ui->inpParamSearchPath->setText( sets.value("ParamSearchPath",dataPath).toString() );
    ui->lisParamSearchResult->setEnabled(false);
    ui->butParamSearchGenerate->setEnabled(false);
    ui->lblParamSearchLog->hide();
    ui->lblParamSearchLog->setText(" ");
    bSearchParamsRunning = false;

    // --------------------------
    // --- TAB  Data ------------
    // --------------------------
    ui->butDataSetMask->setEnabled(false);
    ui->butDataCopyScaling->setEnabled(false);
    ui->grpExtractImage->setEnabled(false);
    ui->grpNoFitRegions->setEnabled(false);
    ui->butDataFindCenter->setEnabled(false);

    ui->cbsExtractScale->clear();
    ui->cbsExtractScale->addItem(" * 1 ", 1 );
    ui->cbsExtractScale->addItem(" * 2 ", 2 );
    ui->cbsExtractScale->addItem(" * 4 ", 4 );

    // --------------------------
    // --- TAB  Calculations ----
    // --------------------------
    ui->radNewImageCal->setChecked(true);
    calcGui = new SC_CalcGUI;
    // Check number of threads
    if ( calcGui->gpuAvailable() )
    {
        ui->inpNumCores->setValue(0);       // 0=GPU, =SpecialText, =Default
    }
    else
    {
        ui->inpNumCores->setSpecialValueText("");
        ui->inpNumCores->setMinimum(1);     // no GPU possible
        ui->inpNumCores->setValue( nthreads ); // set max num of threads as default
    }
    ui->intGridPoints->setValue( sets.value("GridPoints",64).toInt() );  // calculate BeamCenter
    ui->butCalc->setEnabled(true);
    ui->butAbbruch->setEnabled(false);
    ui->progressBar->setEnabled(false);
    ui->radQ1->setChecked( sets.value("Quadrant",2).toInt() == 1 );
    ui->radQ2->setChecked( sets.value("Quadrant",2).toInt() == 2 );
    ui->radQ4->setChecked( sets.value("Quadrant",2).toInt() == 4 );
    ui->togExpandImage->setChecked( sets.value("ExpandImg",true).toBool() );
    ui->togExpandImage->setEnabled( sets.value("Quadrant",2).toInt() != 4 );

    ui->radCenterBeam->setChecked( sets.value("UseCenterBeam",false).toBool() );
    ui->radCenterMidpoint->setChecked( ! ui->radCenterBeam->isChecked() );

    ui->lblLoadPrompt->hide();
    ui->lblLoadFilename->hide();

    // Initialize the color marker labels
    QPalette pal=ui->lblColMarkTrainTbl->palette();
    pal.setBrush(QPalette::Window,SETCOLMARK_TRAINTBL); ui->lblColMarkTrainTbl->setPalette(pal);
    pal.setBrush(QPalette::Window,SETCOLMARK_PARDIFF);  ui->lblColMarkParDiff->setPalette(pal);
    pal.setBrush(QPalette::Window,SETCOLMARK_CHGIMG);   ui->lblColMarkChgImg->setPalette(pal);
    pal.setBrush(QPalette::Window,SETCOLMARK_OLDPARAM); ui->lblColMarkOldPar->setPalette(pal);
    showColorMarker( SETCOLMARK_CLEARED ); // Alles wegnehmen

    myGuiParam::setAllGuiParams(ui->tabCalc);
    //myGuiParam::debugGuiParams();
    calcGui->loadFitFlags(sets);

    QStringList slverthdr;
    slverthdr << "Wavelength [nm]" << "Sample-Det Dist [m]" << "Pix Width [mm]"
              << "Pix Height [mm]" << "Pix Rows" << "Pix Cols" << "Beam Cent X"
              << "Beam Cent Y";
    ui->tblHeaderData->blockSignals(true);
    ui->tblHeaderData->setRowCount(slverthdr.size());
    ui->tblHeaderData->setColumnCount(1);
    ui->tblHeaderData->setVerticalHeaderLabels(slverthdr);
    ui->tblHeaderData->setHorizontalHeaderLabels(QStringList()<<"Value");
    int wv = ui->tblHeaderData->fontMetrics().horizontalAdvance("  Value  ");
    int wl = ui->tblHeaderData->fontMetrics().horizontalAdvance(slverthdr[1]);
    ui->tblHeaderData->setColumnWidth(0,wv);
    ui->tblHeaderData->setMinimumWidth( wl+wv+40 );
    // Defaultwerte setzen
    ui->tblHeaderData->setItem(tblLattNbCols,0, new QTableWidgetItem("200") );
    ui->tblHeaderData->setItem(tblLattNbRows,0, new QTableWidgetItem("200") );
    ui->tblHeaderData->setItem(tblLattBeamX,0, new QTableWidgetItem("100") );
    ui->tblHeaderData->setItem(tblLattBeamY,0, new QTableWidgetItem("100") );
    ui->tblHeaderData->setItem(tblLattWaveLen,0, new QTableWidgetItem("0.5") );
    ui->tblHeaderData->setItem(tblLattSaDist,0, new QTableWidgetItem("4.0") );
    ui->tblHeaderData->setItem(tblLattPixelX,0, new QTableWidgetItem("5.3") );
    ui->tblHeaderData->setItem(tblLattPixelY,0, new QTableWidgetItem("5.3") );
    ui->tblHeaderData->blockSignals(false);
    on_tblHeaderData_cellChanged(0,0);  // damit die Qmax Berechnung startet
    ui->radEditQmaxPreset->setChecked(true);

    // --------------------------
    // --- TAB  FFT -------------
    // --------------------------
    // Grp Input Preprocessing
    ui->radFFTLinInput->setChecked(  sets.value("FFTLinInput",true).toBool() );
    ui->radFFTLogInput->setChecked( !sets.value("FFTLinInput",true).toBool() );
    ui->togFFTscaleImg->setChecked( sets.value("FFTScaleImg",true).toBool() );
    ui->togFFTclipImg->setChecked( sets.value("FFTclipImg",false).toBool() );
    ui->togFFTclip40Img->setChecked( sets.value("FFTclip40Img",false).toBool() );
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
    ui->inpFFTcutX->setValue( sets.value("FFTcutOutX",64).toInt() );
    ui->inpFFTcutY->setValue( sets.value("FFTcutOutY",64).toInt() );
    ui->togFFTcutEnabled->setChecked( sets.value("FFTcutOutEna",false).toBool() );
    ui->inpFFTcutX->setEnabled( ui->togFFTcutEnabled->isChecked() );
    ui->inpFFTcutY->setEnabled( ui->togFFTcutEnabled->isChecked() );
    switch ( sets.value("FFToutput",0).toInt() )
    {
    case 0: ui->radFFToutReal->setChecked(true);     break;
    case 1: ui->radFFToutImag->setChecked(true);     break;
    case 2: ui->radFFToutBetrag->setChecked(true);   break;
    case 3: ui->radFFToutSpectrum->setChecked(true); break;
    }
    ui->butIFFT->setEnabled(false); // Init
    ui->butFFT->setEnabled(false);  // Init
    ui->butFFTverifyRphi->hide();   // war nur für mich zum Testen, daher im Normalbetrieb wegnehmen

    // -------------------------------------
    // --- TAB  Simplex 2D Fit -------------
    // -------------------------------------
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

    // --------------------------
    // --- TAB  AI --------------
    // --------------------------
    // Die globalen Daten richtig übernehmen
    fillDataFromInputs();

    sets.beginGroup("AI");
    ui->togAILoadFromParam->setChecked( sets.value("LoadFromParam").toBool() );
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
    QStringList hdr;
    hdr << "Name" << "Start" << "End / Variation" << "Step / Calculation";
    ui->tblListe->setColumnCount( hdr.size() );
    ui->tblListe->setHorizontalHeaderLabels( hdr );
    ui->tblListe->resizeColumnsToContents();

    ui->butAIcheck->setEnabled(false);        // Erst freigegeben, wenn die AI-Tabelle gefüllt ist
    ui->butAIsaveBkgOp->setEnabled(false);    // -"-
    ui->butAIstart->setEnabled(false);        // -"-
    ui->butSaveImagesAI->setEnabled(false);   // Nur Freigegeben, wenn images da sind
    ui->butRemoveVar->setEnabled(false);

    ui->lblAIbackProc->setEnabled(false);
    ui->lisAIbackProgOut->setEnabled(false);
    ui->lblTPVbackProc->setEnabled(false);
    ui->lisTPVbackProgOut->setEnabled(false);
    aiBackProg = nullptr;

    ui->cbsVariables->clear();
    ui->cbsVariables->addItems( calcGui->paramsForMethod(true,true,false) ); // nur numerische Werte, incl. globales

    // -------------------------------------------
    // --- TAB  Training Parameters Variation ----
    // -------------------------------------------
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
    ui->togTPVaddLines->setChecked( sets.value("ena_addLines",false).toBool() );
    ui->inpTPVaddLinesH->setValue( sets.value("val_addLinesH",1).toInt() );
    ui->inpTPVaddLinesV->setValue( sets.value("val_addLinesV",1).toInt() );
    ui->inpTPVaddLinesHwidth->setValue( sets.value("val_addLinesHwidth",2).toInt() );
    ui->inpTPVaddLinesVwidth->setValue( sets.value("val_addLinesVwidth",2).toInt() );
    ui->inpTPVaddLinesH->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesV->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesHwidth->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesVwidth->setEnabled(ui->togTPVaddLines->isChecked());
    ui->lblTPVaddLines->setEnabled(ui->togTPVaddLines->isChecked());

    ui->togTPVaddNoise->setChecked( sets.value("ena_addNoise",true).toBool() );
    // TODO: Convolute
    ui->togTPVconvolute->setToolTip("Not yet implemented");
    ui->togTPVconvolute->setChecked( false /*sets.value("ena_convolute",true).toBool()*/ );
    ui->togTPVconvolute->setEnabled(false);
    //
    ui->togTPVcalcRphi->setChecked( sets.value("ena_calcRphi",true).toBool() );
    ui->togTPVcalcFFT->setChecked( sets.value("ena_calcFFT",true).toBool() );
    ui->togTPVsaveBaseExtra->setChecked( sets.value("ena_saveExtra",false).toBool() );
    ui->togTPVgeneratePNG->setChecked( sets.value("ena_generatePNG",true).toBool() );
    ui->togTPVscaleScat->setChecked( sets.value("ena_scaleScat",false).toBool() );
    ui->inpTPVnumimg->setValue( sets.value("val_numimg",10).toInt() );
    ui->inpTPVoutPath->setText( sets.value("val_outPath",".").toString() );
    sets.endGroup();

    // --------------------------
    // --- TAB  ChatBot ---------
    // --------------------------
    //ui->tabMain->removeTab( ui->tabMain->indexOf(ui->tabChatBot) ); // Zum Verstecken der kompletten Funktion
#ifdef ChatbotIgnoreImages
    ui->togChatbotClpShowImg->hide(); // Sofern die PNG-Files
#endif

    connect( qApp->clipboard(), SIGNAL(changed(QClipboard::Mode)),
             this, SLOT(chatbotClipboardChanged(QClipboard::Mode)) );
    ui->butChatbotReadClipboard->setEnabled(false);
    chatbotConsProg = getConsoleExecutable(false); // Sonst irritiert die Meldung beim Start den User
    chatbotBackProg = nullptr;
    sets.beginGroup("ChatBot");
    ui->inpChatbotFile->setText( sets.value("LastFile","").toString() );
    ui->inpChatbotTrainfile->setText( sets.value("LastTrainFile","").toString() );
    ui->cbsChatbotColor->setCurrentIndex( sets.value("ImgColor",0).toInt() );
    ui->cbsChatbotRotation->setCurrentIndex( sets.value("ImgRotate",0).toInt() );
    ui->cbsChatbotZoom->setCurrentIndex( sets.value("ImgZoom",0).toInt() );
    ui->inpChatbotLogfile->setText( sets.value("LogFile","").toString() );
    ui->grpChatbotLogfile->setChecked( sets.value("LogEnabled",true).toBool() );
    ui->togChatbotHor->setChecked( sets.value("ImgHor",false).toBool() );
    ui->togChatbotVert->setChecked( sets.value("ImgVert",false).toBool() );
    sets.endGroup();
    ui->butChatbotStop->setEnabled(false);
    if ( chatbotConsProg.isEmpty() )
    {
        ui->lisChatbotLogOut->clear();
        ui->butChatbotStart->setEnabled(false);
        chatbotBackProgAddLog("Background process exec not found but needed for this function.");
    }

    // --------------------------
    // --- TAB  Imagewindows ----  Nur unter Android sichtbar und verwendet
    // --------------------------
#ifdef IMG_ONE_WINDOW
    qDebug() << "********************* ANDROID *******************************";


#else
    ui->tabMain->removeTab( ui->tabMain->indexOf(ui->tabImages) );
#endif

    // -------------------------
    // Abschluss ---------------
    // -------------------------
    ui->tabMain->setCurrentWidget( ui->tabCalc );
    //DT( qDebug() << "... adjustSize ..." );
    //adjustSize();

    QTimer::singleShot( 300, this, SLOT(initCbsAfterStart()) ); // Konstruktor für cbs Trigger und adjustSize()
}


/**
 * @brief SC_MainGUI::~SC_MainGUI
 * Default destructor.
 */
SC_MainGUI::~SC_MainGUI()
{
    delete ui;
}


/**
 * @brief SC_MainGUI::closeEvent [SLOT]
 * @param event
 * Called if the window is about to close, saves all global settings read
 * in the constructor.
 */
void SC_MainGUI::closeEvent(QCloseEvent *event)
{
    if ( !closeMainActive )
    {
        closeMainActive = true; // Avoid calls to imageWindowClosed()

        if ( chatbotBackProg != nullptr && chatbotBackProg->state() != QProcess::NotRunning )
            chatbotBackProg->kill();

        QSettings sets(SETT_APP,SETT_GUI);
        // Calculation & Configuration Tab (Default group)
        calcGui->saveFitFlags(sets);
        sets.setValue( "LastMethod", "FCC Spheres" /*ui->tabMethods->currentIndex()*/ );
        sets.setValue( "GridPoints", ui->intGridPoints->value() );
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
        sets.setValue( "UseCenterBeam", ui->radCenterBeam->isChecked() );

        if ( fitparams != nullptr )
        {
            // FIT Tab (Subgroups)
            sets.beginGroup( "Fit-Limits" );
            QHash< QString/*name*/, _fitLimits* >::const_iterator il = fitparams->constBegin();
            while ( il != fitparams->constEnd() )
            {
                sets.setValue( il.key(), QString("%1:%2:%3:%4").arg(il.value()->used).arg(il.value()->min).arg(il.value()->max).arg(il.value()->fitType) );
                ++il;
            }
            sets.endGroup();
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
        sets.setValue( "FFTScaleImg", ui->togFFTscaleImg->isChecked() );
        sets.setValue( "FFTclipImg", ui->togFFTclipImg->isChecked() );
        sets.setValue( "FFTclip40Img", ui->togFFTclip40Img->isChecked() );
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
        sets.setValue( "FFTcutOutX", ui->inpFFTcutX->value() );
        sets.setValue( "FFTcutOutY", ui->inpFFTcutY->value() );
        sets.setValue( "FFTcutOutEna", ui->togFFTcutEnabled->isChecked() );
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
        sets.setValue( "LoadFromParam", ui->togAILoadFromParam->isChecked() );
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
        sets.setValue("val_addLinesHwidth",ui->inpTPVaddLinesHwidth->value());
        sets.setValue("val_addLinesVwidth",ui->inpTPVaddLinesVwidth->value());
        sets.setValue("ena_addBS",ui->togTPVaddBS->isChecked());
        sets.setValue("ena_addNoise",ui->togTPVaddNoise->isChecked());
        sets.setValue("ena_convolute",ui->togTPVconvolute->isChecked());
        sets.setValue("ena_calcRphi",ui->togTPVcalcRphi->isChecked());
        sets.setValue("ena_calcFFT",ui->togTPVcalcFFT->isChecked());
        sets.setValue("ena_saveExtra",ui->togTPVsaveBaseExtra->isChecked());
        sets.setValue("ena_generatePNG",ui->togTPVgeneratePNG->isChecked());
        sets.setValue("ena_scaleScat",ui->togTPVscaleScat->isChecked());
        sets.setValue("val_numimg",ui->inpTPVnumimg->value());
        sets.setValue("val_outPath",ui->inpTPVoutPath->text());
        sets.endGroup();

        sets.beginGroup("ChatBot");
        sets.setValue("LastFile",ui->inpChatbotFile->text());
        sets.setValue("LastTrainFile",ui->inpChatbotTrainfile->text());
        sets.setValue("ImgColor",ui->cbsChatbotColor->currentIndex());
        sets.setValue("ImgRotate",ui->cbsChatbotRotation->currentIndex());
        sets.setValue("ImgZoom",ui->cbsChatbotZoom->currentIndex());
        sets.setValue("LogFile",ui->inpChatbotLogfile->text());
        sets.setValue("LogEnabled",ui->grpChatbotLogfile->isChecked());
        sets.setValue("ImgHor",ui->togChatbotHor->isChecked());
        sets.setValue("ImgVert",ui->togChatbotVert->isChecked());
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


/**
 * @brief SC_MainGUI::on_actionExit_triggered [SLOT]
 * Exit function from the menu bar.
 */
void SC_MainGUI::on_actionExit_triggered()
{
    closeEvent(nullptr);
    qApp->quit();
}


/**
 * @brief SC_MainGUI::fillDataFromInputs
 * Copies the values from the input fields into the internal Hash-Structures.
 */
void SC_MainGUI::fillDataFromInputs()
{
    DT( qDebug() << "fillDataFromInputs()" );
    // Globale Inputs für die Berechnungen
    //QHash<QString,Double3> SC_CalcGUI::inpVectors;
    //QHash<QString,double>  SC_CalcGUI::inpValues;
    //QHash<QString,double>  SC_CalcGUI::inpSingleValueVectors;

    SC_CalcGUI::inpValues.insert( "GridPoints",    ui->intGridPoints->value() );
    SC_CalcGUI::inpValues.insert( "HKLmax",        ui->intHKLmax->value() );
    SC_CalcGUI::inpValues.insert( "RadioButtonQ1", ui->radQ1->isChecked() );
    SC_CalcGUI::inpValues.insert( "RadioButtonQ2", ui->radQ2->isChecked() );
    SC_CalcGUI::inpValues.insert( "RadioButtonQ4", ui->radQ4->isChecked() );
    SC_CalcGUI::inpValues.insert( "ExpandImage",   ui->togExpandImage->isChecked() );
    SC_CalcGUI::inpValues.insert( "BeamPosX",      ui->inpBeamPosX->value() );
    SC_CalcGUI::inpValues.insert( "BeamPosY",      ui->inpBeamPosY_BeamPosX->value() );

    SC_CalcGUI::inpVectors.insert( "Ax1", Double3(ui->inpVAx1->value(),
                                                  ui->inpAy1->value(),
                                                  ui->inpAz1->value()) );
    SC_CalcGUI::inpVectors.insert( "Ax2", Double3(ui->inpVAx2_VAx1->value(),
                                                  ui->inpAy2_Ay1->value(),
                                                  ui->inpAz2_Az1->value()) );
    SC_CalcGUI::inpVectors.insert( "Ax3", Double3(ui->inpVAx3_VAx1->value(),
                                                  ui->inpAy3_Ay1->value(),
                                                  ui->inpAz3_Az1->value()) );
    SC_CalcGUI::inpVectors.insert( "SigXYZ", Double3(ui->inpSigX->value(),
                                                     ui->inpSigY_SigX->value(),
                                                     ui->inpSigZ_SigX->value()) );

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
    /*
    SC_CalcGUI::inpSingleValueVectors.insert( "Editxrel",    SC_CalcGUI::inpVectors["Nvec"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Edityrel",    SC_CalcGUI::inpVectors["Nvec"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editzrel",    SC_CalcGUI::inpVectors["Nvec"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editx1rel",   SC_CalcGUI::inpVectors["Uvec"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Edity1rel",   SC_CalcGUI::inpVectors["Uvec"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editz1rel",   SC_CalcGUI::inpVectors["Uvec"].z() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editx2rel",   SC_CalcGUI::inpVectors["Vvec"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Edity2rel",   SC_CalcGUI::inpVectors["Vvec"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editz2rel",   SC_CalcGUI::inpVectors["Vvec"].z() );
    */
    SC_CalcGUI::inpSingleValueVectors.insert( "Editdom1",    SC_CalcGUI::inpVectors["SigXYZ"].x() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editdom2",    SC_CalcGUI::inpVectors["SigXYZ"].y() );
    SC_CalcGUI::inpSingleValueVectors.insert( "Editdom3",    SC_CalcGUI::inpVectors["SigXYZ"].z() );
}


/**
 * @brief SC_MainGUI::prepareCalculation
 * @param progbar - enable state of the progressBar
 * This procedure disables the gui elements during the calculation and starts it.
 */
void SC_MainGUI::prepareCalculation( bool progbar )
{
    ui->butCalc->setEnabled(false);
    ui->butFitStart->setEnabled(false);
    ui->butIFFT->setEnabled(false); // Prepare Calc
    ui->butFFT->setEnabled(false); // Prepare Calc
    ui->butAbbruch->setEnabled(true);
    ui->butAbbruch->setText("Abort");
    ui->progressBar->setEnabled(progbar);
    qApp->processEvents();
    fillDataFromInputs();
    calcGui->prepareCalculation( false );
}


/**
 * @brief SC_MainGUI::finishCalculation
 * @param showtime - if true, the calculation time was shown in the gui.
 * This procedure is called to reenable all gui elements.
 */
void SC_MainGUI::finishCalculation(bool showtime)
{
    if ( showtime )
    {
        ui->lblTimeUsed->setText( QString("%1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimBoth)) );
        ui->statusbar->showMessage( QString("%1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimBoth)), 5000 );
    }
    calcGui->updateOutputData();
    ui->progressBar->setValue(100); // it is possible (due to roundings) that the 100% is not reached. So set it here to avoid blinking.
    ui->butCalc->setEnabled(true);
    ui->butAbbruch->setEnabled(false);
    ui->butAbbruch->setText("Abort");
    ui->progressBar->setEnabled(false);
    ui->butFitStart->setEnabled(oneFitParamUsed);
    ui->butIFFT->setEnabled(true);  // Finish Calc
    ui->butFFT->setEnabled(true);  // Finish Calc
    calcGui->inpValues.insert( "_CalcTime_", calcGui->higResTimerElapsed(SC_CalcGUI::htimCalc) );
    calcGui->inpValues.insert( "_PrepTime_", calcGui->higResTimerElapsed(SC_CalcGUI::htimPrep) );
}


/**
 * @brief SC_MainGUI::on_butCalc_clicked
 * Button "Calculate", save the current values to a temporary parameter file before.
 */
void SC_MainGUI::on_butCalc_clicked()
{
#ifdef FITDATA_IN_GPU  // butCalc
    if ( sender() != nullptr )
        calcGui->setFitData( 0, 0, nullptr );
#endif

    // Da es mir schon öfter passiert ist, dass ich Parameter ausgesucht und dann das
    // Programm einfach geschlossen hatte, ohne die Parameter zu sichern, werden bei
    // jedem Berechnungsstart alle Parameter in einer temporären Datei gespeichert.
    performSaveParamOperation( fnTempParamFile );

    local_butCalc_clicked();
}


/**
 * @brief SC_MainGUI::local_butCalc_clicked
 * Start the calculation in a thread without saving the current values.
 */
void SC_MainGUI::local_butCalc_clicked()
{
    if ( _calcThread == nullptr )
        _calcThread = new myCalcThread(calcGui);
    _calcThread->setThreads( ui->inpNumCores->value() );
    _calcThread->setIgnoreNewSwitch( ui->togIgnNewSwitch->isChecked() );

    prepareCalculation( true );

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
        zzmin = - ui->intGridPoints->value();
        zzrange = 2 * ui->intGridPoints->value();
    }
    else
    {
        zzmin = 0;
        zzrange = ui->intGridPoints->value();
    }
    if ( ui->radQ1->isChecked() )
    {
        iimin = 0;
        iirange = ui->intGridPoints->value();
    }
    else
    {
        iimin = - ui->intGridPoints->value();
        iirange = 2 * ui->intGridPoints->value();
    }
    lastTime = 0;

    _bAbbruch = false;
    _bAbortTimer = false;
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

    if ( !timingTestsRunning && !bIgnoreRecalc && autoProcessingFile!=nullptr && autoProcessingFile->isOpen() )
    {
        QTimer::singleShot( 200, this, SLOT(autoProcessingTimer()) );
    }
}


/**
 * @brief myCalcThread::myCalcThread
 * @param cg     - CalcGUI object to be accessible by the thread
 * @param parent - normal parent object for the thread
 * Constructs the calculation thread
 */
myCalcThread::myCalcThread(SC_CalcGUI *cg, QObject *parent ) : QThread(parent)
{
    //_exit = false;
    calcGuiThread = cg;
}

/**
 * @brief myCalcThread::run
 * Thread run function, perform the real calculation
 */
void myCalcThread::run()
{
    calcGuiThread->doCalculation( numThreads, bIgnoreNewSwitch );
}


/**
 * @brief SC_MainGUI::automaticRecalc [SLOT]
 * Called if any of the parameter values are changed and calls the calculation if this
 * is enabled with the QCheckBox 'togAutoRecalc'.
 */
void SC_MainGUI::automaticRecalc()
{
    DT(QObject *oo = sender();
       if ( oo != nullptr )
         qDebug() << "automaticRecalc(), ignore:" << bIgnoreRecalc << (static_cast<QWidget*>(oo))->objectName();
       else
         qDebug() << "automaticRecalc(), ignore:" << bIgnoreRecalc;
       )
    if ( bIgnoreRecalc ) return; // Werte wurden vom Programm geändert
    if ( fitIsRunning ) return; // Bei einem Fit auch nichts automatisch neu berechnen
    if ( timingTestsRunning ) return;   // Timing Tets
    QObject *o = sender();
    if ( o != nullptr )
    {
        QWidget *w = static_cast<QWidget*>(o);
        if ( w->objectName().contains("LType",Qt::CaseInsensitive) )
            return;
        if ( w->objectName().contains("Ordis",Qt::CaseInsensitive) )
            return;
    }
    automaticRecalcDoit();
}

/**
 * @brief SC_MainGUI::automaticRecalcDoit
 */
void SC_MainGUI::automaticRecalcDoit()
{
    DT( qDebug() << "automaticRecalcDoit(), toggle:" << ui->togAutoRecalc->isChecked() );
    if ( fitIsRunning ) return; // Bei einem Fit auch nichts automatisch neu berechnen
    if ( bIgnoreRecalc ) return; // Werte wurden vom Programm geändert
    if ( timingTestsRunning ) return;
    if ( ! ui->togAutoRecalc->isChecked() ) return; // Recalc vom User gesperrt
    if ( _calcThread != nullptr && _calcThread->isRunning() ) return; // Läuft schon

    on_butCalc_clicked(); // hier mit(!) speichern der Parameter im Temp-File

    // Wenn die Rechnung länger als 1,5 sec dauert, dann diese Funktion abschalten
    if ( calcRunTime.elapsed() > 1500 )
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
    performTimingTests("");
}

void SC_MainGUI::performTimingTests(QString out)
{
    dlgTimeTests *dlg = new dlgTimeTests(dataPath,calcGui->gpuAvailable(),ui->inpNumCores->maximum(),this);
    if ( out.isEmpty() /*aus menu*/ )
    {
        dlg->setHKLmaxUsed( myGuiParam::isEnabled("HKLmax") );
        if ( dlg->exec() != QDialog::Accepted ) return;
    }

    int numLoops = dlg->getRepetitions();
    QVector<int> threads   = dlg->getThreads();
    QVector<int> hklmax    = dlg->getHKLmax();
    QVector<int> quadrants = dlg->getQuadrants();
    dlgTimeTests::_newSwitch nswitch = dlg->getNewSwitchFlag();
    bool nswitchModFN      = out.isEmpty() ? dlg->getNewSwitchModFN() : false;
    QString saveFilename   = out.isEmpty() ? dlg->getSaveFilename() : out;
    QString saveBasePath   = "";    // for Images to save
    if ( ! saveFilename.isEmpty() )
    {
        if ( saveFilename.indexOf("/") < 0 )
            saveFilename = dataPath+"/"+saveFilename;
        QFileInfo fi(saveFilename);
        saveBasePath = fi.absolutePath()+"/"+fi.baseName(); // save images
        if ( nswitchModFN )
        {
            saveBasePath += "_" + dlg->newSwitch2Str(nswitch);
            saveFilename = saveBasePath + "." + fi.suffix();
            //qDebug() << saveFilename;
            //return;
        }
        QFile f(saveFilename);
        if ( ! f.open(QIODevice::Append) )
            if ( ! f.open(QIODevice::WriteOnly) )
            {
                QMessageBox::critical( this, "Timing tests", "Error opening the File\n "+saveFilename+"\n"+f.errorString(),
                                      QMessageBox::Ok );
                qDebug() << saveFilename << f.errorString();
                return;
            }
        if ( autoProcessingFile != nullptr /*sender()->objectName().isEmpty()*/ )
            f.write( qPrintable("% Generated from Autoprocessingfile: " + autoProcessingFile->fileName() + EOL) );
        else
            f.write( qPrintable("% " + dlg->getComment() + EOL) );
        f.write( qPrintable("% Loaded parameter: "+ui->lblLoadFilename->text()+EOL) );
        f.write( "Threads & HKLmax & Quadrants & Switch & Min/ Mean/ Max (Prep) in ms" EOL );
        f.close();
    }

    double c, sumCalc, minCalc, maxCalc, p, sumPrep, minPrep, maxPrep;
    QStringList latex;
    dlgTimeTests::_newSwitch nswcurr = nswitch;
    bool savExpand = ui->togExpandImage->isChecked();
    ui->togExpandImage->setChecked(true);
    bool savQ1 = ui->radQ1->isChecked();
    bool savQ2 = ui->radQ2->isChecked();
    bool savQ4 = ui->radQ4->isChecked();
    bool savUpd = ui->togIgnoreUpdates->isChecked();
    bool savNSw = ui->togIgnNewSwitch->isChecked();
    if ( out.isEmpty() && dlg->getEnaUpdates() )
    {   // GUI-Updates sind gewünscht
        ui->togIgnoreUpdates->setChecked( false );
        bIgnoreUpdates = false;
    }
    else
    {   // Keine Updates bei den Berechnungen
        ui->togIgnoreUpdates->setChecked( true );
        bIgnoreUpdates = true;
    }
    bool saveImages = out.isEmpty() ? dlg->getEnaSaveImages() : false;
    timingTestsRunning = true;

    // ***** Loop 1 ***** Quadrants
    foreach ( int quadr, quadrants )
    {
        if ( quadr > 0 )
        {
            ui->radQ1->setChecked( quadr == 1 );
            ui->radQ2->setChecked( quadr == 2 );
            ui->radQ4->setChecked( quadr == 4 );
        }

        // ***** Loop 2 ***** HKLmax
        for ( int h=0; h<hklmax.size(); h++ )
        {
            ui->intHKLmax->setValue( hklmax[h] );

            // ***** Loop 3 ***** Threads
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

                // ***** Loop 4 ***** NewSwitch
                if ( nswitch/*Soll*/ == dlgTimeTests::swBoth )
                    nswcurr = dlgTimeTests::swOld;
                ui->togIgnNewSwitch->setChecked( nswcurr == dlgTimeTests::swOld );
                for ( int sw=0; sw<2; sw++ )
                {   // Es können maximal 2 Durchläufe sein...

                    sumCalc=0;
                    minCalc=10000000;
                    maxCalc=0;
                    sumPrep=0;
                    minPrep=10000000;
                    maxPrep=0;
                    for ( int i=0; i<numLoops; i++ )
                    {
                        statusBar()->showMessage( QString("TimingTest: Threads=%1, HKLmax=%2, Quadrants=%5, Switch=%6, Loop %3/%4")
                                                     .arg(ui->inpNumCores->value()).arg(hklmax[h])
                                                     .arg(i+1).arg(numLoops).arg(quadr).arg(dlg->newSwitch2Str(nswcurr)) );
                        local_butCalc_clicked(); // ohne speichern der Parameter im Temp-File
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
                        lastUsedImage->saveImage( saveBasePath + QString("-t=%1-hkl=%2-q=%3.png").arg(ui->inpNumCores->value()).arg(hklmax[h]).arg(quadr) );
                    }
                    qDebug() << "ERG T=" << ui->inpNumCores->value() << "HKLmax=" << hklmax[h]
                             << "Q=" << quadr << "SW=" << dlg->newSwitch2Str(nswcurr)
                             << "Calc=" << minCalc << sumCalc/numLoops << maxCalc
                             << "Prep=" << minPrep << sumPrep/numLoops << maxPrep;
                    // LaTeX Ausgabe mit den Rundungen (am Ende auf der Konsole oder in eine Datei)
                    if ( saveFilename.isEmpty() )
                        latex << QString("%1 & %2 & %7 & %8 & %3/ %4/ %5 (%6)").arg(ui->inpNumCores->value()).arg(hklmax[h])
                                     .arg(minCalc,0,'f',3).arg(sumCalc/numLoops,0,'f',3).arg(maxCalc,0,'f',3)
                                     .arg(sumPrep/numLoops,0,'f',3).arg(quadr).arg(dlg->newSwitch2Str(nswcurr));
                    else
                    {
                        QFile f(saveFilename);
                        if ( ! f.open(QIODevice::Append) )
                        {
                            QMessageBox::critical( this, "Timing tests", "Error opening the File\n "+saveFilename+"\n"+f.errorString(),
                                                  QMessageBox::Ok );
                            qDebug() << saveFilename << f.errorString();
                            saveFilename = ""; // Damit unten nicht die gleiche Meldung schon wieder kommt
                            _bAbbruch = true;
                            break;
                        }
                        if ( f.isOpen() )
                        {
                            f.write( QString("%1 & %2 & %8 & %9 & %3/ %4/ %5 (%6)%7")
                                        .arg(ui->inpNumCores->value()).arg(hklmax[h])
                                        .arg(minCalc,0,'f',3).arg(sumCalc/numLoops,0,'f',3).arg(maxCalc,0,'f',3)
                                        .arg(sumPrep/numLoops,0,'f',3).arg(EOL).arg(quadr).arg(dlg->newSwitch2Str(nswcurr)).toLatin1() );
                            f.close();
                        }
                    }
                    if ( _bAbbruch ) break;

                    //typedef enum { swBoth, swNew, swOld } _newSwitch;
                    if (nswitch/*Soll*/ == dlgTimeTests::swBoth &&
                        nswcurr/*Ist */ == dlgTimeTests::swOld)
                    {   // Jetzt muss die Schleife nochmals wiederholt werden
                        nswcurr = dlgTimeTests::swNew;
                        ui->togIgnNewSwitch->setChecked( nswcurr == dlgTimeTests::swOld );
                    }
                    else
                        sw = 10;
                } // for sw
                if ( _bAbbruch ) break;
            } // for t
            if ( _bAbbruch ) break;
        } // for h
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
    ui->togIgnNewSwitch->setChecked( savNSw );
    ui->togIgnoreUpdates->setChecked( savUpd );
    ui->togExpandImage->setChecked( savExpand );
    ui->radQ1->setChecked( savQ1 );
    ui->radQ2->setChecked( savQ2 );
    ui->radQ4->setChecked( savQ4 );
    bIgnoreUpdates = savUpd;
    timingTestsRunning = false;
    if ( !bIgnoreRecalc && autoProcessingFile!=nullptr && autoProcessingFile->isOpen() )
    {
        QTimer::singleShot( 200, this, SLOT(autoProcessingTimer()) );
    }
}

/**
 * @brief SC_MainGUI::on_butAbbruch_clicked
 * The current calculation should be aborted.
 * TODO: das Programm kann crashen, wenn die Berechnung abgebrochen wird.
 *       Das passiert leider im Debugger nie, sodass eine genaue Analyse
 *       fast unmöglich wird.
 */
void SC_MainGUI::on_butAbbruch_clicked()
{
    _bAbbruch = true;
    ui->butAbbruch->setEnabled(false);
    ui->butAbbruch->setText("Aborting...");
    qApp->processEvents();
}

/**
 * @brief SC_MainGUI::updateProgBar
 * @param val - percentage of the progress bar or -1 to denote GPU-Usage
 */
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
    fn = QFileDialog::getSaveFileName( this, "Save Parameter", fn, "Parameter (*.ini)" );
                                      //nullptr, QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".ini",Qt::CaseInsensitive) ) fn += ".ini";
    data.setValue("LastParam",fn);
    performSaveParamOperation( fn );
    ui->lblLoadPrompt->show();
    ui->lblLoadFilename->show();
    ui->lblLoadFilename->setText(fn);
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
    sets.setValue( "EditAxis1x",     ui->inpVAx1->value() );
    sets.setValue( "EditAxis1y",     ui->inpAy1->value() );
    sets.setValue( "EditAxis1z",     ui->inpAz1->value() );
    sets.setValue( "EditAxis2x",     ui->inpVAx2_VAx1->value() );
    sets.setValue( "EditAxis2y",     ui->inpAy2_Ay1->value() );
    sets.setValue( "EditAxis2z",     ui->inpAz2_Az1->value() );
    sets.setValue( "EditAxis3x",     ui->inpVAx3_VAx1->value() );
    sets.setValue( "EditAxis3y",     ui->inpAy3_Ay1->value() );
    sets.setValue( "EditAxis3z",     ui->inpAz3_Az1->value() );
    //sets.setValue( "Editxrel",       ui->inpN1->value() );
    //sets.setValue( "Edityrel",       ui->inpN2->value() );
    //sets.setValue( "Editzrel",       ui->inpN3->value() );
    sets.setValue( "Editdom1",       ui->inpSigX->value() );
    sets.setValue( "Editdom2",       ui->inpSigY_SigX->value() );
    sets.setValue( "Editdom3",       ui->inpSigZ_SigX->value() );
    //sets.setValue( "Editx1rel",      ui->inpU1->value() );
    //sets.setValue( "Edity1rel",      ui->inpU2->value() );
    //sets.setValue( "Editz1rel",      ui->inpU3->value() );
    //sets.setValue( "Editx2rel",      ui->inpV1->value() );
    //sets.setValue( "Edity2rel",      ui->inpV2->value() );
    //sets.setValue( "Editz2rel",      ui->inpV3->value() );
    sets.setValue( "HKLmax",     ui->intHKLmax->value() );
    sets.setValue( "GridPoints", ui->intGridPoints->value() );
    sets.setValue( "Threads",        ui->inpNumCores->value() );
    //sets.setValue( "CurMethod",      ui->tabMethods_o->tabText( ui->tabMethods_o->currentIndex() ) );
    sets.setValue( "EditCenterX",    ui->inpBeamPosX->value() );
    sets.setValue( "EditCenterY",    ui->inpBeamPosY_BeamPosX->value() );
    sets.setValue( "Comment", ui->txtComment->text() );
    sets.endGroup();
    // TODO
    sets.beginGroup( "AI" );
    //sets.setValue( "LoadFromParam", ui->togAILoadFromParam->isChecked() );
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
                                       "Parameter (*.ini);;Scatter-Parameter (*.par)" //);
                                      , nullptr, QFileDialog::DontUseNativeDialog | QFileDialog::DontResolveSymlinks );
    if ( fn.isEmpty() ) return;
    DT( QElapsedTimer et;   et.start() );
    data.setValue("LastParam",fn);
    data.sync();
    DT( qDebug() << "on_actionLoad_all_Parameters_triggered()" << et.elapsed() );
    on_butResetColorMarker_clicked();
    DT( qDebug() << "vor load" << et.elapsed() );
    QString rv = local_Load_all_Parameters(fn);
    if ( ! rv.isEmpty() )
    {
        //qDebug() << rv;
        if ( rv.contains("Error:") )
        {
            ui->lblLoadPrompt->hide();
            ui->lblLoadFilename->hide();
            ui->lblLoadFilename->setText("");
            QMessageBox::critical( this, "Load parameter ERROR", fn+EOL+rv, QMessageBox::Ok );
        }
        else
            QMessageBox::warning( this, "Load parameter WARNING", fn+EOL+rv, QMessageBox::Ok );
    }
    DT( qDebug() << "LOAD TIMER" << et.elapsed() );
}

void SC_MainGUI::on_actionLoad_last_calculaition_parameters_triggered()
{
    local_Load_all_Parameters(fnTempParamFile); // Hier kann kein Fehler auftreten
}


/**
 * @brief SC_MainGUI::loadParameter_checkOldFormat
 * @param cur - current value from parameter file
 * @param old - key from old file
 * @param rv  - message string, remove the key
 * @return if the key is in the rv then this value is returned, else return the cur
 */
double SC_MainGUI::loadParameter_checkOldFormat(double cur, QString key, QString &rv )
{
    int p1 = rv.indexOf(key+"=");
    if ( p1 < 0 ) return cur;
    int p2 = rv.indexOf(",",p1);
    if ( p2 < 0 ) p2 = rv.length();
    int p3 = p1 + key.length()+1;
    double val = rv.midRef(p3,p2-p3).toDouble();
    rv = rv.left(p1) + rv.mid(p2+2);
    calcGui->updateParamValueColor( key, SETCOLMARK_CLEARED );
    return val;
}

QString SC_MainGUI::local_Load_all_Parameters(QString fn)
{
    bIgnoreRecalc = true;  // Load start
    DT( qDebug() << "local_Load_all_Parameters()" << fn );
    if ( fn.endsWith(".par",Qt::CaseInsensitive) )
    {
        qDebug() << "loadScatterParameter" << fn;
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
    on_butResetColorMarker_clicked();   // Parameters not loaded are marked

    //QList<QWidget*> wid = ui->tabCalc->findChildren<QWidget*>();
    foreach ( QWidget *w, ui->tabCalc->findChildren<QWidget*>() ) w->blockSignals(true);

    bool hklloaded, gridloaded;

    // Wenn im aktuellen Entwicklungsschritt (nur die Methode "Generic" bekannt) alte Files geladen werden
    // sollen, die noch mehrere Methoden enthalten, dann gibt es den Schlüssel "CurMethod" in der Gruppe
    // "Inputs". Dies soll dann die Methode sein, die geladen werden soll, auch wenn dieser Name hier nicht
    // bekannt ist. Zugleich wird der Wert vom Parameter LType auf die Methode gesetzt. Somit können die
    // Parameterfiles aus dem alten Format wieder normal gelesen werden.
    QSettings sets( fn, QSettings::IniFormat );
    QStringList gr = sets.childGroups();
    gr.removeOne("Inputs");
    gr.removeOne("FFT");
    gr.removeOne("AI");     // TODO: neue feste Gruppen hier ausblenden
    //qDebug() << "Groups:" << gr;
    if ( gr.size() > 1 )
    {   // Jetzt gibt es mehr als eine Methode
        sets.beginGroup( "Inputs" );
        QString m = sets.value("CurMethod","").toString();
        sets.endGroup();
        if ( ! m.isEmpty() )
        {
            QString mm = m;
            if ( mm.indexOf(" ") > 0 ) mm.truncate(mm.indexOf(" "));
            qDebug() << "Load spec" << m << mm;

            // Suchen des richtigen LType-Wertes durch scan der Combobox-Texte
            for ( int val=0; val<ui->cbsLType->count(); val++ )
                if ( ui->cbsLType->itemText(val).startsWith(mm,Qt::CaseInsensitive) )
                {
                    calcGui->updateParamValue( "LType", val, SETCOLMARK_IGNORED, false );
                    break;
                }
            /*
            // Suchen des richtigen LType-Wertes mittels setzen und rücklesen
            int val = 0;
            while ( true )
            {
                calcGui->updateParamValue( "LType", val, SETCOLMARK_IGNORED, false );
                QString tmp = calcGui->currentParamValueStr( "LType", true ); // hier in aktueller Methode suchen...
                //qDebug() << "    " << val << tmp;
                if ( tmp.isEmpty() || tmp[0] == '?' ) break; // Ende der Liste oder Fehlermeldung
                if ( tmp.startsWith(mm,Qt::CaseInsensitive) ) break; // gefunden
                val++;
            }
            */
            rv = calcGui->loadParameter(fn,"@"+m,hklloaded,gridloaded);
            isloaded = true;
        }
    }

    if ( !isloaded )
        rv = calcGui->loadParameter(fn,"",hklloaded,gridloaded);

    // Load all common parameters
    sets.beginGroup( "Inputs" );
    ui->radQ1->setChecked( sets.value( "RadioButtonQ1", false ).toBool() );
    ui->radQ2->setChecked( sets.value( "RadioButtonQ2", false ).toBool() );
    ui->radQ4->setChecked( sets.value( "RadioButtonQ4", false ).toBool() );
    ui->togExpandImage->setChecked( sets.value( "ExpandImage", false ).toBool() );
    ui->togExpandImage->setEnabled( ! ui->radQ4->isChecked() );
    ui->inpVAx1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis1x",0.0).toDouble(), "VAx1", rv ) );
    ui->inpAy1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis1y",0.0).toDouble(), "Ay1", rv ) );
    ui->inpAz1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis1z",0.0).toDouble(), "Az1", rv ) );
    ui->inpVAx2_VAx1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis2x",0.0).toDouble(), "VAx2", rv ) );
    ui->inpAy2_Ay1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis2y",0.0).toDouble(), "Ay2", rv ) );
    ui->inpAz2_Az1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis2z",0.0).toDouble(), "Az2", rv ) );
    ui->inpVAx3_VAx1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis3x",0.0).toDouble(), "VAx3", rv ) );
    ui->inpAy3_Ay1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis3y",0.0).toDouble(), "Ay3", rv ) );
    ui->inpAz3_Az1->setValue( loadParameter_checkOldFormat( sets.value("EditAxis3z",0.0).toDouble(), "Az3", rv ) );
    ui->inpSigX->setValue( loadParameter_checkOldFormat( sets.value("Editdom1",0.0).toDouble(), "SigX", rv ) );
    ui->inpSigY_SigX->setValue( loadParameter_checkOldFormat( sets.value("Editdom2",0.0).toDouble(), "SigY", rv ) );
    ui->inpSigZ_SigX->setValue( loadParameter_checkOldFormat( sets.value("Editdom3",0.0).toDouble(), "SigZ", rv ) );
    ui->inpBeamPosX->setValue( sets.value( "EditCenterX", 0 ).toDouble() );
    ui->inpBeamPosY_BeamPosX->setValue( sets.value( "EditCenterY", 0 ).toDouble() );
    //qDebug() << "Load CenterBeam:" << calcGui->currentParamValueInt("CenterBeam");
    if ( calcGui->currentParamValueInt("CenterBeam")>0 )
        ui->radCenterBeam->setChecked( true );
    else
        ui->radCenterMidpoint->setChecked( true );

    // EditHKLmax und EditGridPoints sind jetzt in die normalen Parameter gewandert...
    //qDebug() << "*** LOAD ***" << hklloaded << gridloaded;
    if ( !hklloaded )
        ui->intHKLmax->setValue( loadParameter_checkOldFormat( sets.value("HKLmax",3).toInt(), "HKLmax", rv ) );
    if ( !gridloaded )
        ui->intGridPoints->setValue( loadParameter_checkOldFormat( sets.value("GridPoints",100).toInt(), "GridPoints", rv ) );

    // Werte aus dem 'rv' und den ColorMarkern entfernen, die im alten Format (noch) keinen Sinn gemacht hatten.
    loadParameter_checkOldFormat( 0.0, "CalcQmax", rv );
    loadParameter_checkOldFormat( 0.0, "CenterBeam", rv );

    if ( rv.endsWith(", ") ) rv.chop(2);
    //qDebug() << "LoadParameters:" << rv;
    if ( rv.contains("Unknown parameters in file") )
        showColorMarker( SETCOLMARK_OLDPARAM );

    // Wenn im Parametersatz die GPU aktiviert ist, im aktuellen System diese aber nicht vorhanden ist,
    // so würde das Eingabefeld auf 1 Thread schalten. Das ist unklug und wird hier auf die maximale
    // Anzahl von Threads umgebogen. Kommt vor, wenn Parameterfiles vom Linux (mit GPU) zum Laptop
    // (ohne GPU) kopiert werden.
    int thr = sets.value( "Threads", 0 ).toInt();
    if ( thr == 0 /*GPU*/ && ui->inpNumCores->minimum() == 1 /*keine GPU vorhanden*/ )
        thr = ui->inpNumCores->maximum();
    // Andersherum wird bei einem Thread-Parameter > 1 und vorhandener GPU auf die GPU gestellt.
    else if ( thr > 0 && ui->inpNumCores->minimum() == 0 /*GPU ist vorhanden*/ )
        thr = 0;
    ui->inpNumCores->setValue( thr );

    ui->txtComment->setText( sets.value("Comment","").toString() );
    if ( ui->txtComment->text().isEmpty() )
        ui->txtComment->setText( QFileInfo(fn).baseName() );

    sets.endGroup();
    if ( ui->togAILoadFromParam->isChecked() )
    {
        sets.beginGroup( "AI" );
        //ui->togAILoadFromParam->setChecked( sets.value("LoadFromParam").toBool() );
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
    }
    on_cbsLType_currentIndexChanged( ui->cbsLType->currentIndex() );
    on_cbsComboBoxParticle_currentIndexChanged( ui->cbsComboBoxParticle->currentIndex() );
    on_cbsOrdis_currentIndexChanged( ui->cbsOrdis->currentIndex() );
    on_cbsComboBoxInterior_currentIndexChanged( ui->cbsComboBoxInterior->currentIndex() );
    on_cbsComboBoxPeak_currentIndexChanged( ui->cbsComboBoxPeak->currentIndex() );

    foreach ( QWidget *w, ui->tabCalc->findChildren<QWidget*>() ) w->blockSignals(false);

    bIgnoreRecalc = false;  // Load done
    DT( qDebug() << "... local load::adjustSize ..." );
    adjustSize();   // local load
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
    ui->butFFT->setEnabled( images.size() > 0 );  // Add Img

    int cntNaN = img->setData( x0, x1, y0, y1, d );
    if ( meta )
    {   // Meta-Daten übertragen
        // Kommentar wird als Info über dem Image verwendet
        img->addMetaInfo( "_Calculation_", ui->txtComment->text() ); // calcGui->curMethod->subCalc->methodName() );
        // Das Flag für den neuen Switch wird ebenfalls kopiert
        img->addMetaInfo( "_OptimizedSwitch_", calcGui->getCalcPtr()->newSwitchUsed() ? "True" : "False" );
        if ( cntNaN > 0 )
            img->addMetaInfo( "_NaNvalues_", QString::number(cntNaN) );
        // Methodenspezifische Metadaten
        QHash<QString,paramHelper*>::iterator ip = calcGui->params.begin();
        while ( ip != calcGui->params.end() )
        {
            img->addMetaInfo( ip.key(), //ip.value()->str );
                              calcGui->currentParamValueStr( ip.key(), true ) );
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
                // and find the index in the selection boxes
                int id = ui->cbsFitImageWindows->findText( images.at(i)->windowTitle() );
                ui->cbsFitImageWindows->removeItem( id );
                ui->cbsFFTWindows->removeItem( id ); // same list
                // if this was the last image, disable the fit button
                if ( ui->cbsFitImageWindows->count() == 0 )
                {
                    ui->butFitStart->setEnabled(false);
                    ui->butIFFT->setEnabled(false); // Img Window Close
                    ui->butFFT->setEnabled(false); // Img Window Close
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
    DT( qDebug() << "on_lisDataWindows_currentTextChanged()" );
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

    if ( !timingTestsRunning && !bIgnoreRecalc && autoProcessingFile!=nullptr && autoProcessingFile->isOpen() )
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

    // Größe des Quellbildes
    int sx = images[curimg]->xmax() - images[curimg]->xmin();
    int sy = images[curimg]->ymax() - images[curimg]->ymin();

    // Beamstop
    int bsx = images[curimg]->getFileInfos()->centerX;
    int bsy = images[curimg]->getFileInfos()->centerY;
    if ( bsx == 0 || bsy == 0 )
    {
        bsx = sx / 2;
        bsy = sy / 2;
    }
    else if ( images[curimg]->xmin() < 0 || images[curimg]->ymin() < 0 )
    {   // Falls das Image um den Nullpunkt berechnet wird, so sollte der
        // Beamstop verschoben werden, da in der weiteren Berechnung das
        // Image immer im 1. Quadranten betrachtet wird.
        bsx -= images[curimg]->xmin();
        bsy -= images[curimg]->ymin();
    }

    //qDebug() << "performIFFTcalculations:"
    //         << "X" << images[curimg]->getFileInfos()->li.centerX << images[curimg]->xmin()
    //         << "Y" << images[curimg]->getFileInfos()->li.centerY << images[curimg]->ymin();

/*
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
        {
            if ( src[i] > 0 )
                data[i] = (log10(src[i]) - vlogmin) / (vlogmax-vlogmin);
            else
                data[i] = 0;
        }
    }
*/

    double *data = new double[sx*sy];
    memcpy( data, images[curimg]->dataPtr(), sx*sy*sizeof(double) );

    SasCalc_PostProc::inst()->scaleAndClipData( data/*in and out*/, sx*sy,
                                                ui->togFFTscaleImg->isChecked(),
                                                ui->togFFTclipImg->isChecked(),
                                                ui->togFFTclip40Img->isChecked(),
                                                ui->radFFTLogInput->isChecked() );

    metaData.insert( "From Image", tit );

    int srphi = ui->cbsFFTsizeRphi->currentText().left(3).trimmed().toInt();
    int sout  = ui->cbsFFTsizeOut->currentText().left(3).trimmed().toInt();

    double *outr = nullptr;

    if ( ui->grpFFTuseRphi->isChecked() )
    {   // Jetzt soll das (r,phi) verwendet werden.

        // Das (r,phi) Bild berechnen
        auto anfr = std::chrono::high_resolution_clock::now();
        outr = SasCalc_PostProc::inst()->generateRphi(0, sx, 0, sy, bsx, bsy, srphi, data/*in*/,
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
        //qDebug() << "FFT: rphi";
    }
    else //if ( data != nullptr )
    {
        fftinp = data; // Image lin/log
        x0 = images[curimg]->xmin();
        x1 = images[curimg]->xmax();
        y0 = images[curimg]->ymin();
        y1 = images[curimg]->ymax();
        //qDebug() << "FFT: data log";
    }
    /*else
    {
        fftinp = images[curimg]->dataPtr(); // Lin
        x0 = images[curimg]->xmin();
        x1 = images[curimg]->xmax();
        y0 = images[curimg]->ymin();
        y1 = images[curimg]->ymax();
        qDebug() << "FFT: data lin";
    }*/

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
                                                   ui->togIFFTSwap->isChecked(),
                                                   ui->togFFTcutEnabled->isChecked()?ui->inpFFTcutX->value():0,
                                                   ui->togFFTcutEnabled->isChecked()?ui->inpFFTcutY->value():0 );

    auto endi = std::chrono::high_resolution_clock::now();
    auto timi = std::chrono::duration_cast<std::chrono::duration<float>>(endi-anfi);
    //qDebug() << "Timing (IFFT):" << timi.count()*1000.0 << "ms";
    if ( outi != nullptr )
    {
        // (Jan2023): in der Mitte einen Teil rausschneiden, da dort die Informationen sind ==> calculateIFFT()
        int soutx=sout, souty=sout;
        if ( ui->togFFTcutEnabled->isChecked() )
        {
            soutx = ui->inpFFTcutX->value();
            souty = ui->inpFFTcutY->value();
        }
        widImage* img = addImage( true, 0, soutx, 0, souty, outi, "iFFT-Image", false );
        QHash<QString,QString>::iterator ii = metaData.begin();
        while ( ii != metaData.end() )
        {
            img->addMetaInfo( ii.key(), ii.value() );
            ++ii;
        }
        img->addMetaInfo( "NbRows", QString::number(souty) );
        img->addMetaInfo( "NbCols", QString::number(soutx) );
        img->addMetaInfo( "_CalcTime_", QString("%1 ms").arg(timi.count()*1000.0) );
        img->addMetaInfo("@",""); // Sortieren
    }
    // Next is allways a new image ...
    ui->radNewImageCfg->setChecked(true);
}

void SC_MainGUI::on_intGridPoints_valueChanged(int arg1)
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


void SC_MainGUI::on_tblHeaderData_cellChanged(int DT(row), int /*column*/)
{
    DT( qDebug() << "on_tblHeaderData_cellChanged()" << row );
    // Berechne qmax (siehe sascalc_fcc_gpu.cu: doIntCalc_FCC_F)
    // tblLattice3DValues --> tblHeaderData
    int cols = ui->tblHeaderData->item(tblLattNbCols,0)->text().toInt();
    int rows = ui->tblHeaderData->item(tblLattNbRows,0)->text().toInt();
    double posx = ui->tblHeaderData->item(tblLattBeamX,0)->text().toDouble();
    double posy = ui->tblHeaderData->item(tblLattBeamY,0)->text().toDouble();
    double wlen_m = ui->tblHeaderData->item(tblLattWaveLen,0)->text().toDouble()*1.0E-9; /*nm->m*/
    double dist_m = ui->tblHeaderData->item(tblLattSaDist,0)->text().toDouble();
    double pixx_m = ui->tblHeaderData->item(tblLattPixelX,0)->text().toDouble()*1.0E-3; /*mm->m*/
    double pixy_m = ui->tblHeaderData->item(tblLattPixelY,0)->text().toDouble()*1.0E-3; /*mm->m*/

    int ihex_max = qMax( cols - posx, posx );
    int i_max    = qMax( rows - posy, posy );
    double qmax_cols = ( 2. * M_PI / (wlen_m * dist_m) ) * ihex_max * pixx_m * 1.0E-9; /* wieder in nm zurück */
    double qmax_rows = ( 2. * M_PI / (wlen_m * dist_m) ) * i_max    * pixy_m * 1.0E-9; /* wieder in nm zurück */

    calcGui->updateParamValue( "CalcQmax", qMax(qmax_cols,qmax_rows), SETCOLMARK_IGNORED );
    // Damit beim LoadParameter auch der Wert erhalten bleibt.
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
    fn = QFileDialog::getOpenFileName( this, "Load Data File Image", fn, filter.join(";;"), &curFilter );
                                       //QFileDialog::DontUseNativeDialog );
    if ( fn.isEmpty() ) return;
    data.setValue("LastImage",fn);
    lastDataFile = fn;
    if ( ! local_OpenMeasFile(fn,nullptr) )
    {
        qDebug() << fn << "unknown format";
        QMessageBox::critical( this, "Open measurement file", "Unknown file format in\n"+fn, QMessageBox::Ok );
    }
}

bool SC_MainGUI::local_OpenMeasFile( QString fn, widImage **imgout )
{
    QElapsedTimer tt;
    tt.start();
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
    //qDebug() << "ReadData (1)" << tt.elapsed();
    if ( img == nullptr ) return false; // Lesefehler oder unbekannte Dateiendung

    if ( imgout != nullptr ) *imgout = img;

    // Jetzt werden noch Meta-Daten in die GUI kopiert. Das mache ich bei allen Datentypen,
    //  da ich nicht wissen kann, ob diese Meta-Infos nicht doch enthalten sind.
    copyMetaToLatticeTable( img );
    //qDebug() << "ReadData (2)" << tt.elapsed();

    // Datenfiles werden nie bei neuen Berechnungen automatisch überschrieben!
    ui->radLastImageCfg->setChecked( false );
    ui->radNewImageCfg->setChecked( true );
    img->show();
    return true;
}

void SC_MainGUI::copyMetaToLatticeTable( widImage *img )
{
    DT( qDebug() << "copyMetaToLatticeTable()" );
    bool savign = bIgnoreRecalc;
    bIgnoreRecalc = true;
    // tblLattice3DValues --> tblHeaderData
    ui->tblHeaderData->setItem( tblLattWaveLen, 0, new QTableWidgetItem(img->metaInfo("wavelength",
                                     calcGui->currentParamValueStr("EditWavelength",true))) );
    ui->tblHeaderData->setItem( tblLattSaDist, 0, new QTableWidgetItem(img->metaInfo("SampleDist",
                                     calcGui->currentParamValueStr("EditDet",true))) );
    ui->tblHeaderData->setItem( tblLattPixelX, 0, new QTableWidgetItem(img->metaInfo("Pixel_X",
                                     QString::number(calcGui->currentParamValueDbl("EditPixelX") ))) );  // *1000.
    ui->tblHeaderData->setItem( tblLattPixelY, 0, new QTableWidgetItem(img->metaInfo("Pixel_Y",
                                     QString::number(calcGui->currentParamValueDbl("EditPixelY") ))) );  // *1000.
    ui->tblHeaderData->setItem( tblLattNbRows, 0, new QTableWidgetItem(img->metaInfo("NbRows")) );
    ui->tblHeaderData->setItem( tblLattNbCols, 0, new QTableWidgetItem(img->metaInfo("NbCols")) );
    int row = img->metaInfo("NbRows").toInt();
    int col = img->metaInfo("NbCols").toInt();
    int bsx = col/2 - ui->inpBeamPosX->value();
    int bsy = row/2 - ui->inpBeamPosY_BeamPosX->value();
    ui->tblHeaderData->setItem( tblLattBeamX, 0, new QTableWidgetItem(img->metaInfo("BeamPosX", QString::number(bsx))) );
    ui->tblHeaderData->setItem( tblLattBeamY, 0, new QTableWidgetItem(img->metaInfo("BeamPosY", QString::number(bsy))) );
    bIgnoreRecalc = savign;
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
    QString val = calcGui->currentParamValueStr( key, true );
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
        qDebug() << fn << f.errorString();
        QMessageBox::critical( this, "Save AI infos", fn+"\n"+f.errorString(), QMessageBox::Ok );
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
        qDebug() << fn << f.errorString();
        QMessageBox::critical( this, "Load AI infos", fn+"\n"+f.errorString(), QMessageBox::Ok );
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
    fn = QFileDialog::getExistingDirectory( this, "Directory for inputfiles", ui->inpFileName->text(),
                                           QFileDialog::ShowDirsOnly /*| QFileDialog::DontUseNativeDialog*/ );
    if ( fn.isEmpty() ) return;
    ui->inpFileName->setText( fn );
}

void SC_MainGUI::on_butSelectDir_clicked()
{
    QString fn = ui->inpSubDir->text();
    fn = QFileDialog::getExistingDirectory( this, "Save files into", fn, QFileDialog::ShowDirsOnly /*| QFileDialog::DontUseNativeDialog*/ );
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
        QMessageBox::critical( this, "Save AI data", fn+"\n"+f.errorString(), QMessageBox::Ok );
        return false;
    }
    if ( !doOverwrite ) f.seek( f.size() );

    _loopDefinition rv;

    if ( useTPV )
    {
        rv = getTPVLoopDefinition();
        f.write( qPrintable("PATH|"+ui->inpTPVoutPath->text()+EOL) );
        f.write( "COLTBL|*DAT*" EOL ); // Spezielle Kennung für das Datenformat
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

    if ( useTPV && !ui->togTPVcalcFFT->isChecked() && ui->togTPVcalcRphi->isChecked() )
    {   // FFT weglassen, aber r,phi soll bestimmt werden
        QString rphiScale = "";
        if ( ui->togFFTscaleRphi->isChecked()  ) rphiScale += "Scale ";
        if ( ui->togFFTclipRphi->isChecked()   ) rphiScale += "Clip1 ";
        if ( ui->togFFTclip40Rphi->isChecked() ) rphiScale += "Clip4 ";
        f.write( qPrintable(QString("RPHI|%1;%2")
                               .arg(ui->cbsFFTsizeRphi->currentText().left(3).trimmed(),rphiScale.trimmed())
                           +EOL) );
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
                           //ui->inpFFTcutX->setValue( sets.value("FFTcutOutX",64).toInt() );
                           //ui->inpFFTcutY->setValue( sets.value("FFTcutOutY",64).toInt() );
                           +EOL) );
    }
    if ( useTPV )
    {   // Weitere Spezial-Parameter zur Bild-Nachbearbeitung
        // Diese erst nach der FFT / RPHI Definition schrieben, damit der Logoutput später die richtigen Daten anzeigt
        QString tpv = "TPV|xx"; // Damit immer 2 Einträge da sind ...
        if ( ui->togTPVaddBS->isChecked() )
            tpv += QString("|BS%1;%2;%3;%4").arg(ui->inpTPVbeamx0->value()).arg(ui->inpTPVbeamy0->value())
                       .arg(ui->inpTPVnbeam->value()).arg(ui->inpTPVmbeam->value());
        if ( ui->togTPVaddLines->isChecked() )
            tpv += QString("|LI%1;%2;%3;%4").arg(ui->inpTPVaddLinesH->value()).arg(ui->inpTPVaddLinesV->value())
                       .arg(ui->inpTPVaddLinesHwidth->value()).arg(ui->inpTPVaddLinesVwidth->value());
        if ( ui->togTPVaddNoise->isChecked() )
            tpv += "|NO";
        if ( ui->togTPVconvolute->isChecked() )
            tpv += "|CO";
        if ( ui->togTPVcalcRphi->isChecked() )
            tpv += "|RP";
        if ( ui->togTPVsaveBaseExtra->isChecked() )
            tpv += "|IX";
        if ( ui->togTPVgeneratePNG->isChecked() )
            tpv += "|GP";
        if ( ui->togTPVscaleScat->isChecked() )
            tpv += "|SC";
        f.write( qPrintable(tpv+EOL) );
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


QString SC_MainGUI::getConsoleExecutable(bool domsg)
{
    QString cons = "sas_scatter2Cons";
#ifdef Q_OS_WIN
    // Unter Windows muss hier das passende Exe-File gesucht werden. Unter Linux muss
    // im User-Bin Verzeichnis (im PATH) ein passender Link vorhanden sein.
    QDir d( QDir::currentPath());
    QString dBase = d.absolutePath();
    while ( d.absolutePath().contains("build-") ) d.cdUp();
    QFileInfoList fil = d.entryInfoList( QStringList()<<"*"+cons+"*" );
    //qDebug() << "FIL" << fil;
    //qDebug() << "CUR" << QDir::currentPath();
    //Debug: FIL (QFileInfo(C:\SimLab\sas-crystal\build-sas_scatter2Cons-Desktop_Qt_5_15_2_MinGW_64_bit-Debug),
    //            QFileInfo(C:\SimLab\sas-crystal\build-sas_scatter2Cons-Desktop_Static-Release))
    //Debug: CUR           "C:/SimLab/sas-crystal/build-sas_scatter2-Desktop_Qt_5_15_2_MinGW_64_bit-Debug"

    while ( fil.size() > 1 )
    {
        if ( fil[0].absoluteFilePath().right(30) != dBase.right(30) )
            fil.takeFirst();
        else if ( fil.size() > 1 )
            fil.takeLast();
    }
    if ( fil.size() == 0 )
    {
        if ( domsg ) QMessageBox::critical( this, "Executable", "The executable "+cons+"\nis not found.", QMessageBox::Ok );
        return "";
    }
    /*if ( fil.size() != 1 )
    {
        if ( domsg ) QMessageBox::critical( this, "Executable", "The executable "+cons+"\nis found in more than one version.", QMessageBox::Ok );
        return;
    }*/
    QString tmp = fil[0].absoluteFilePath();
    if ( tmp.endsWith("Release") )
        tmp += "/release/" + cons + ".exe";
    else if ( tmp.endsWith("Debug") )
        tmp += "/debug/" + cons + ".exe";
    else
        tmp += "/" + cons + ".exe";
    //qDebug() << tmp;
    if ( ! QFileInfo::exists(tmp) )
    {
        if ( domsg ) QMessageBox::critical( this, "Executable", "The executable "+tmp+"\nnot found.", QMessageBox::Ok );
        return "";
    }
    cons = tmp;
#else
    Q_UNUSED(domsg) // Avoid compiler warning
#endif
    return cons;
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

    QString cons = getConsoleExecutable(true);
    if ( cons.isEmpty() ) return;

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
    //qDebug() << cons << params.join(" ");

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
    // Aufbau der Tabelle pro Zeile:
    // [0]=Variablenname
    // [1]=Startwert:   <zahl> oder '*'<zahl>
    // [2]=Endwert:     <zahl> oder '*'<zahl> oder 'V'<zahl>
    // [3]=Schrittweite <zahl> oder 'R'<zahl> oder '='<formel> oder '/'<zahl>
    for ( int r=0; r<ui->tblListe->rowCount(); r++ )
        if ( ! ui->tblListe->item(r,3)->text().isEmpty() /*&&
             master->metaData.contains(ui->tblListe->item(r,0)->text())*/ )
        {   // Werte aus anderen Berechnungsformeln und leere Schrittweiten werden nicht weiter verfolgt
            _werte *w = new _werte;
            w->start = 0;
            w->end = 0;
            w->step = 0;
            //w->current = 0;
            w->var = 0;
            w->randcount = 0;
            w->formel = "";

            if ( ui->tblListe->item(r,1)->text().startsWith("*") )
            {   // Start berechnet aus aktuellem Parameter
                w->start = calcGui->currentParamValueDbl( ui->tblListe->item(r,0)->text() )
                           * ui->tblListe->item(r,1)->text().midRef(1).toDouble();
                qDebug() << "START" << w->start;
            }
            else
            {   // Start ist normale Zahl
                w->start  = ui->tblListe->item(r,1)->text().toDouble();
            }

            if ( ui->tblListe->item(r,2)->text().startsWith("V") )
            {   // V<num> Variation
                w->var = ui->tblListe->item(r,2)->text().midRef(1).toDouble();
                w->end = w->start;
            }
            else if ( ui->tblListe->item(r,2)->text().startsWith("*") )
            {   // Ende berechnet aus aktuellem Parameter
                w->end = calcGui->currentParamValueDbl( ui->tblListe->item(r,0)->text() )
                         * ui->tblListe->item(r,2)->text().midRef(1).toDouble();
                //w->var = 0;
                qDebug() << "ENDE" << w->end;
            }
            else
            {   // Normaler Endwert
                w->end = ui->tblListe->item(r,2)->text().toDouble();
                //w->var = 0;
            }

            if ( ui->tblListe->item(r,3)->text().startsWith("R") )
            {   // R<num> ist Zufallszahlen
                w->randcount = ui->tblListe->item(r,3)->text().midRef(1).toInt();
            }
            else if ( ui->tblListe->item(r,3)->text().startsWith("=") )
            {   // = ... calculation ...
                w->formel = ui->tblListe->item(r,3)->text().mid(1).trimmed();
            }
            else if ( ui->tblListe->item(r,3)->text().startsWith("/") )
            {   // /<zahl> Schrittweite berechnet
                w->step = (w->end - w->start) / ui->tblListe->item(r,3)->text().midRef(1).toDouble();
            }
            else
            {   // <num> ist Schrittweite
                w->step = ui->tblListe->item(r,3)->text().toDouble();
                //w->randcount = 0;
                //w->formel = "";
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

    // Äußerste Schleife über die Methoden ==> es gibt nur noch eine ...
    //QStringList slCalcArt = calcGui->getCalcTypes();
    //for ( int iCurMethod=0; iCurMethod<slCalcArt.size(); iCurMethod++ )
    //{
        calcval.clear();

        // Tabellenwerte auf Start setzen
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
                    QMessageBox::critical( this, "Check AI loop file", finp.fileName()+"\n"+finp.errorString(), QMessageBox::Ok );
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
                        else if ( calcGui->isCurrentParameterValid( keys.at(k), false ) )
                            calcval.insert( keys.at(k), v );
                        else if ( calcGui->isCurrentParameterValid( "Edit"+keys.at(k), false ) )
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
                    if ( daten[k]->formel.isEmpty() && calcGui->isCurrentParameterValid( k, false ) )
                        calcval.insert( k, daten[k]->current );
                foreach ( QString k, slKeys )
                    if ( ! daten[k]->formel.isEmpty() && calcGui->isCurrentParameterValid( k, false ) )
                        calcval.insert( k, evaluateFormula( daten[k]->formel ) );
                calcval.insert( "#Calculation#", 0/*iCurMethod*/ );
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
                cls.replace( "{M}", "DEFAULT" /*calcGui->memory->methodName()*/ ); // Historisch...
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
                    else if ( calcGui->isCurrentParameterValid( k, false ) )
                    {
                        if ( div != 1 )
                        {
                            int tmp = calcGui->currentParamValueInt(k) / div;
                            cls.replace( "{"+korg+"}", QString::number(tmp*div) );
                        }
                        else
                            cls.replace( "{"+korg+"}", calcGui->currentParamValueStr(k,false) );
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
    //} // foreach ( QString curMethod, slCalcArt ) - Schleife über die Methoden
    return retval;
}

double SC_MainGUI::evaluateFormula( QString formel )
{
    QStringList sl = formel.split(QRegExp("[+-*/ ]"),Qt::SkipEmptyParts);
    //qDebug() << "CALC" << formel << sl;
    double val = 0;
    if ( calcval.contains(sl[0]) )
        val = calcval[ sl[0] ];
    else if ( calcGui->isCurrentParameterValid(sl[0],false) )
        val = calcGui->currentParamValueDbl(sl[0]);
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
    prepareCalculation( true );

    // TODO: LogTimer starten (TestGo)

    int mx = ui->intGridPoints->value();
    double *data = new double[4*mx*mx]; // mx is one quarter of result image

    QString title = "??";
    SasCalc_PostProc::inst()->setLogging( false );
    double fac;

    switch ( ui->cbsTestImgSelect->currentIndex() )
    {
    case 0: // Test Image 1 (-)
        title = "Test-Image (-)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( x >= mx-2 && x <= mx+2 ) // Linie entlang y in der Mitte
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
        break;
    case 1: // Test Image 2 (|)
        title = "Test-Image (|)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( y >= mx-2 && y <= mx+2 ) // Linie entlang x in der Mitte
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
        break;
    case 2: // Test Image 3 (x)
        title = "Test-Image (x)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( abs(x-y) <= 2 || abs(x-(2*mx-y-1)) <= 2 )
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
        break;
    case 3: // Test Image 4 (+)
        title = "Test-Image (+)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( (x >= mx-2 && x <= mx+2) || (y >= mx-2 && y <= mx+2) )
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
        break;
    case 4: // Test Image 5 (:.)
        title = "Test-Image (:.)";
        fac = 2.0*M_PI / (mx*2);
        for ( int y=0, idx=0; y<2*mx; y++ )
            for ( int x=0; x<2*mx; x++, idx++ )
                if ( x < mx )
                    data[idx] = (1.0 + sin( x * fac ) * sin( y * fac )) * 128.0;  // -> 0..255
                else if ( y > mx )
                    data[idx] = (1.0 + sin( x * fac ) * sin( y * fac )) * 64.0;  // -> 0..128
                else
                    data[idx] = 128;    // Quadrant rechts oben (ohne Drehung) konstant
        break;
    case 5: // Test Image 6 (.)
        title = "Test-Image (.)";
        for ( int x=0, idx=0; x<2*mx; x++ )
            for ( int y=0; y<2*mx; y++, idx++ )
            {
                if ( x >= mx-2 && x <= mx+2 && y >= mx-2 && y <= mx+2 ) // Mittelpunkt
                    data[idx] = 100;
                else
                    data[idx] = 0.1;
            }
        break;
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


void SC_MainGUI::loadScatterParameter( QString fn )
{
    QFile fin(fn);
    if ( ! fin.open(QIODevice::ReadOnly) )
    {
        qDebug() << fin.errorString();
        QMessageBox::critical( this, "Load Scatter Parameter", fn+"\n"+fin.errorString(), QMessageBox::Ok );
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
    setValue( data[27], "HKLmax" );
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
    setValue( data[169], "GridPoints" );
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
    //QString m = ""; // "FCC Spheres"; // ui->tabMethods->tabText(ui->tabMethods->currentIndex());
    //qDebug() << "setScatterValue" << name << data;
    if ( calcGui->updateParamValue( name, data.toDouble(), SETCOLMARK_IGNORED ) ) return;

    if ( name == "EditAxis1x" )
        ui->inpVAx1->setValue( data.toDouble() );
    else if ( name == "EditAxis1y" )
        ui->inpAy1->setValue( data.toDouble() );
    else if ( name == "EditAxis1z" )
        ui->inpAz1->setValue( data.toDouble() );
    else if ( name == "EditAxis2x" )
        ui->inpVAx2_VAx1->setValue( data.toDouble() );
    else if ( name == "EditAxis2y" )
        ui->inpAy2_Ay1->setValue( data.toDouble() );
    else if ( name == "EditAxis2z" )
        ui->inpAz2_Az1->setValue( data.toDouble() );
    else if ( name == "EditAxis3x" )
        ui->inpVAx3_VAx1->setValue( data.toDouble() );
    else if ( name == "EditAxis3y" )
        ui->inpAy3_Ay1->setValue( data.toDouble() );
    else if ( name == "EditAxis3z" )
        ui->inpAz3_Az1->setValue( data.toDouble() );
    //else if ( name == "Editxrel" )
    //    ui->inpN1->setValue( data.toDouble() );
    //else if ( name == "Edityrel" )
    //    ui->inpN2->setValue( data.toDouble() );
    //else if ( name == "Editzrel" )
    //    ui->inpN3->setValue( data.toDouble() );
    else if ( name == "Editdom1" )
        ui->inpSigX->setValue( data.toDouble() );
    else if ( name == "Editdom2" )
        ui->inpSigY_SigX->setValue( data.toDouble() );
    else if ( name == "Editdom3" )
        ui->inpSigZ_SigX->setValue( data.toDouble() );
    //else if ( name == "Editx1rel" || name == "Editu1" )
    //    ui->inpU1->setValue( data.toDouble() );
    //else if ( name == "Edity1rel" || name == "Editu2" )
    //    ui->inpU2->setValue( data.toDouble() );
    //else if ( name == "Editz1rel" || name == "Editu3" )
    //    ui->inpU3->setValue( data.toDouble() );
    //else if ( name == "Editx2rel" || name == "Editv1" )
    //    ui->inpV1->setValue( data.toDouble() );
    //else if ( name == "Edity2rel" || name == "Editv2" )
    //    ui->inpV2->setValue( data.toDouble() );
    //else if ( name == "Editz2rel" || name == "Editv3" )
    //    ui->inpV3->setValue( data.toDouble() );
    else if ( name == "Edithklmax" )
        ui->intHKLmax->setValue( data.toInt() );
    else if ( name == "RadioButtonQ1" )
        ui->radQ1->setChecked( data.toInt() != 0 );
    else if ( name == "RadioButtonQ2" )
        ui->radQ2->setChecked( data.toInt() != 0 );
    else if ( name == "RadioButtonQ4" )
        ui->radQ4->setChecked( data.toInt() != 0 );
    else if ( name == "ExpandImage" )
        ui->togExpandImage->setChecked( data.toInt() != 0 );
    else if ( name == "GridPoints" )
        ui->intGridPoints->setValue( data.toInt() );
    else if ( name == "EditCenterX" )
        ui->inpBeamPosX->setValue( data.toInt() );
    else if ( name == "EditCenterY" )
        ui->inpBeamPosY_BeamPosX->setValue( data.toInt() );
    // Falsche / unbekannte Bezeichnungen
    else
    {
        static QStringList tmp;
        if ( tmp.size() == 0 )
            tmp = calcGui->paramsForMethod(false,true,false);
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
    DT( qDebug() << "copyParamsToFitTable()" );
    QStringList slParams = calcGui->paramsForMethod( false, false, true );
    QStringList slHdr;
    slHdr << "Used" << "Min" << "Start" << "Max" << "Result";

    ui->tblFitValues->clear();
    ui->tblFitValues->setRowCount( slParams.size() );
    ui->tblFitValues->setColumnCount( slHdr.size() );
    ui->tblFitValues->setHorizontalHeaderLabels( slHdr );
    ui->tblFitValues->setVerticalHeaderLabels( slParams );

    if ( fitparams == nullptr ) fitparams = new _param2fitval;

    //qDebug() << "Soll:" << slParams;
    //qDebug() << " Ist:" << fitparams->keys();

    ui->tblFitValues->blockSignals(true);
    oneFitParamUsed = false;

    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup("Fit-Limits");

    // Alle aktiven Parameter in die Fit-Liste packen
    foreach (QString p, slParams)
    {
        _fitLimits *fl = fitparams->value(p,nullptr);
        if ( fl == nullptr )
        {
            QStringList slVal = sets.value(p,"0:0:0:0").toString().split(":",SPLIT_KEEP_EMPTY_PARTS);
            while ( slVal.size() < 4 ) slVal << "0";
            fl = new _fitLimits;
            fl->fitvalid = false;
            fl->used     = slVal[0].toInt() != 0;
            fl->min      = slVal[1].toDouble();
            fl->max      = slVal[2].toDouble();
            if ( fl->min >= fl->max || fl->max == 0 )
                fl->fitType = _fitTypes::fitNone;       // Jetzt stimmen die Grenzen nicht, nehme die interne Struktur
            else
                fl->fitType  = static_cast<_fitTypes>(slVal[3].toInt());
            fitparams->insert( p, fl );
        }
        fl->orgval = fl->fitstart = calcGui->currentParamValueDbl( p );

        bool cnt;
        double min, max;
        if ( calcGui->limitsOfParamValue( p, min, max, cnt ) )
        {
            if ( fl->fitType == _fitTypes::fitNone )
            {   // Nur, wenn die Grenzen nicht definiert sind, diese aus der internen Struktur holen.
                // Sonst sind diese schon angepasst und oben aus der Registry gelesen worden.
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
    sets.endGroup();

    // Zur Sicherheit bei allen Fit-Daten in der Struktur das Use-Flag löschen, wenn es nicht verwendet wird
    QHash< QString/*name*/, _fitLimits* >::const_iterator il = fitparams->constBegin();
    while ( il != fitparams->constEnd() )
    {
        if ( slParams.indexOf(il.key()) == -1 ) il.value()->used = false;
        ++il;
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
    DT( qDebug() << "on_tabMain_currentChanged()" );
    if ( index != ui->tabMain->indexOf(ui->tabFit) ) return;

    if ( fitIsRunning ) return;

    bool fl = copyParamsToFitTable();

    ui->butFitStart->setEnabled( fl && ui->cbsFitImageWindows->count()>0 );
}

void SC_MainGUI::tblFitUsed_toggled(bool checked)
{
    DT( qDebug() << "tblFitUsed_toggled()" );
    QCheckBox *cbs = static_cast<QCheckBox*>(sender());
    oneFitParamUsed = false;
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        if ( static_cast<tblCheckBox*>(ui->tblFitValues->cellWidget(r,0))->tog() == cbs )
        {
            _param2fitval *p2f = fitparams;
            //if ( p2f == nullptr ) return; TODO
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
    _param2fitval *p2f = fitparams;
    //if ( p2f == nullptr ) return; TODO
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
    QString fn = QFileDialog::getSaveFileName( this, "Save Loggings", dataPath, "Logfiles (*.log)" ); //, nullptr, QFileDialog::DontUseNativeDialog );
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
        {
            qDebug() << f.errorString();
            QMessageBox::critical( this, "Save Loggings", fn+"\n"+f.errorString(), QMessageBox::Ok );
        }
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
    if ( curFitImage->myHeight() != curFitImage->myWidth() )
    {
        ui->lisFitLog->addItem("None square images not working '"+ui->cbsFitImageWindows->currentText()+"'");
        ui->butFitStart->setEnabled(true);
        //ui->butFitAutomatic->setEnabled(savAutoEna);
        return;
    }

    // Automatic button enabled => the start button was pressed
    if ( savAutoEna ) fitIsRunning = true;

    // Jetzt erst noch alle (geänderten) Startwerte aus der Tabelle auslesen
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
        double val = ui->tblFitValues->item( r, 2 )->text().trimmed().toDouble();
        fitparams->value(p)->fitstart = val;
        fitparams->value(p)->orgval   = val;
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
        QMessageBox::critical( this, "Fit Loop", fileFitLogName+"\n"+fileFitLog->errorString()+"\n\nNo logging will be saved",
                              QMessageBox::Ok );
        fileFitLog = nullptr;
    }

    curFitImage->getVarScaling( fitOrgMin, fitOrgmax );

    _bAbbruch = false;
    prepareCalculation( true );

    _param2fitval *parameter = fitparams;
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
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
        ui->tblFitValues->item( r, 4 )->setText("?");
        fitparams->value(p)->fitvalid = false;
    }

    timeForAll = 0;
    loopsForAll   = 0;
    imgGenForAll  = 0;
    bool ok;
    double tolerance  = ui->inpFitTolerance->text().toDouble(&ok);
    if ( !ok || (tolerance<=0) || (tolerance>100.0) )
    {   // TODO: bessere Fehlermeldung
        qDebug() << "Tolerance invalid:" << ui->inpFitTolerance->text();
        QMessageBox::critical( this, "Fit input parameter check", "Tolerance value invalid: "+ui->inpFitTolerance->text()+"\n\nSet to 0.001",
                              QMessageBox::Ok );
        tolerance = 0.001;
        ui->inpFitTolerance->setText("0.001");
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
                                  fitparams, retinfo );
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
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
        ui->tblFitValues->item( r, 4 )->setText(QString::number(fitparams->value(p)->fitres));
        fitparams->value(p)->fitvalid = true;
        // Restore the old values in the Calculation Tab
        calcGui->updateParamValue( p, fitparams->value(p)->orgval, SETCOLMARK_IGNORED, false/*ohne debug*/ );
    }
    ui->tblFitValues->resizeColumnsToContents();
    ui->lisFitLog->scrollToBottom();    // if updates are disabled during calculations
}

void SC_MainGUI::on_butFitUseResult_clicked()
{
    bool savign = bIgnoreRecalc;
    bIgnoreRecalc = true;
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
        QString val = ui->tblFitValues->item( r, 4 )->text().trimmed();
        if ( val.isEmpty() || val.contains("?") ) continue;
        calcGui->updateParamValue( p, val.toDouble(), SETCOLMARK_IGNORED, false/*ohne debug*/ );
        ui->tblFitValues->item( r, 2 )->setText( val );
        fitparams->value(p)->fitstart = fitparams->value(p)->fitres;
        fitparams->value(p)->orgval   = fitparams->value(p)->fitres;
    }
    ui->tblFitValues->resizeColumnsToContents();
    // Zur Sicherheit testen, ob bei einem Calculate auf das anzufittende Image ausgegeben werden soll.
    if ( ui->radLastImageCal->isChecked() )
    {
        if ( lastUsedImage != nullptr )
            if ( ui->cbsFitImageWindows->currentText() == lastUsedImage->windowTitle() )
            {
                ui->radNewImageCal->setChecked(true);
            }
    }
    bIgnoreRecalc = savign;
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
            {
                _bAbortTimer = true;
                qDebug() << "   - - - QlogThread: timeout" << calcRunTime.elapsed();
            }
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
        ui->inpBeamPosX->setValue(x-ikws->myWidth()/2.);
        ui->inpBeamPosY_BeamPosX->setValue(y-ikws->myHeight()/2.);
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


/**
 * @brief SC_MainGUI::initCbsAfterStart
 * This function is called with a short delay after the Constructor of the GUI is finished
 * to perform last actions and to make all inputs visible according to the selections of
 * the ComboBoxes.
 * After this it will start the AutoProc timer if an argument is given to this program.
 */
void SC_MainGUI::initCbsAfterStart()
{
    // Jetzt die jeweiligen ComboBoxen triggern, damit die Ena/Dis für die anderen Elemente gemacht werden
    bIgnoreRecalc = true;  // initCbsAfterStart
    on_cbsLType_currentIndexChanged( ui->cbsLType->currentIndex() );
    on_cbsComboBoxParticle_currentIndexChanged( ui->cbsComboBoxParticle->currentIndex() );
    on_cbsOrdis_currentIndexChanged( ui->cbsOrdis->currentIndex() );
    on_cbsComboBoxInterior_currentIndexChanged( ui->cbsComboBoxInterior->currentIndex() );
    on_cbsComboBoxPeak_currentIndexChanged( ui->cbsComboBoxPeak->currentIndex() );
    bIgnoreRecalc = false;  // initCbsAfterStart

    if ( qApp->arguments().size() > 1 )     // [0] = Executable file path
    {
        autoProcessingFile = new QFile( qApp->arguments().at(1) );
        if ( autoProcessingFile->open(QIODevice::ReadOnly) )
        {
            qDebug() << autoProcessingFile->fileName() << "Starting";
            QTimer::singleShot( 100, this, SLOT(autoProcessingTimer()) );   // Start
        }
        else
        {
            qDebug() << autoProcessingFile->fileName() << autoProcessingFile->errorString();
            QMessageBox::critical( this, "Auto processing file (first parameter)",
                                  autoProcessingFile->fileName()+"\n"+autoProcessingFile->errorString(), QMessageBox::Ok );
            autoProcessingFile = nullptr;
            bIgnoreRecalc = false;  // initCbsAfterStart AutoProc-Error
        }
    }
}


void SC_MainGUI::on_actionStart_autoprocessing_file_triggered()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString fn = data.value("LastManAutoProc",".").toString();
    fn = QFileDialog::getOpenFileName(this,"Autoprocessing file", fn );
    if ( fn.isEmpty() ) return;
    data.setValue("LastManAutoProc",fn);
    data.sync();

    autoProcessingFile = new QFile( fn );
    if ( autoProcessingFile->open(QIODevice::ReadOnly) )
    {
        qDebug() << autoProcessingFile->fileName() << "Starting";
        QTimer::singleShot( 100, this, SLOT(autoProcessingTimer()) );   // Start
    }
    else
    {
        qDebug() << autoProcessingFile->fileName() << autoProcessingFile->errorString();
        QMessageBox::critical( this, "Auto processing file",
                              autoProcessingFile->fileName()+"\n"+autoProcessingFile->errorString(), QMessageBox::Ok );
        autoProcessingFile = nullptr;
        bIgnoreRecalc = false;  // manual AutoProc-Error
    }
}


/**
 * @brief SC_MainGUI::autoProcessingTimer
 * This will execute all commands from a text file given by the command line or via menu.
 * The available commands are:
 *  #                          starts a comment line (also empty lines are ignored)
 *  TAB <name>                 select the given tab page in the gui
 *  OPEN MEAS FILE DEFAULT     opens the last used measurement file
 *  OPEN MEAS FILE <filename>  opens the given measurement file
 *  LOAD PARAMS DEFAULT        opens the last used parameter set
 *  LOAD PARAMS <filename>     opens the given parameter set
 *  USE QMAX FROM DATA         use the qmax value calculated from the data (if possible)
 *  USE QMAX PRESET            use the qmax value given by the user
 *  THREADS <n>                sets the number of threads to be used to <n> (0=GPU)
 *  DOCALC                     starts a normal calculation (after this the processings continues)
 *  DOFFT                      starts a FFT calculation with the parameters given in the gui
 *  DOIFFT                     starts a iFFT calculation with the parameters given in the gui
 *  TEST <n>                   starts a test calculation as shown in the Configuration tab
 *  TIMETEST <outfilename>     starts the TimeTest procedure specified in the Tools -> Timing tests dialog
 */
void SC_MainGUI::autoProcessingTimer()
{
    qDebug() << "AutoProc: autoProcessingTimer()";

    if ( autoProcessingFile == nullptr ) return;
    if ( ! autoProcessingFile->isOpen() ) return;

    while ( true )
    {
        if ( autoProcessingFile->atEnd() )
        {
            autoProcessingFile->close();
            qDebug() << "AutoProc: finished";
            bIgnoreRecalc = false;  // AutoProc-Finish
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
                if ( ui->tabMain->tabText(i).startsWith(line.mid(4),Qt::CaseInsensitive) )
                {
                    ui->tabMain->setCurrentIndex(i);
                    found = true;
                    break;
                }
            if ( found ) continue;
            qDebug() << "AutoProc-ERROR: invalid tab name:" << line.mid(4);
            for ( int i=0; i<ui->tabMain->count(); i++ )
                qDebug() << "    :" << ui->tabMain->tabText(i);
            autoProcessingFile->close();
            bIgnoreRecalc = false;  // Init done AutoProc-Error
            return;
        }

        else if ( line.startsWith("OPEN MEAS FILE ") )
        {
            QString fn;
            if ( line.midRef(15).startsWith("DEFAULT") )
            {
                QSettings data(SETT_APP,SETT_GUI);
                fn = data.value("LastImage",dataPath).toString();
            }
            else
                fn = line.mid(15);
            //qDebug() << "AutoProc: Load meas file" << fn;
            if ( local_OpenMeasFile(fn,nullptr) ) continue;
            qDebug() << "AutoProc-ERROR: meas file not found / unknown format:" << fn;
            QMessageBox::critical( this, "Autoprocessing ERROR", "Measurement file not found / unknown format:\n"+fn, QMessageBox::Ok );
            autoProcessingFile->close();
            bIgnoreRecalc = false;  // "OPEN MEAS FILE " AutoProc-Error
            return;
        }

        else if ( line.startsWith("LOAD PARAMS ") )
        {
            QString fn;
            if ( line.midRef(12).startsWith("DEFAULT") )
            {
                QSettings data(SETT_APP,SETT_PAR);
                fn = data.value("LastParam",".").toString();
            }
            else
                fn = line.mid(12);
            //qDebug() << "AutoProc: Load params" << fn;
            if ( QFile::exists(fn) )
            {
                QString rv = local_Load_all_Parameters(fn);
                if ( ! rv.isEmpty() )
                {
                    qDebug() << "AutoProc-ERROR: Parameterfile" << fn;
                    qDebug() << rv;
                    if ( rv.contains("Error:") )
                    {
                        QMessageBox::critical( this, "Autoprocessing ERROR", fn+EOL+rv, QMessageBox::Ok );
                        autoProcessingFile->close();
                        bIgnoreRecalc = false;  // Init done AutoProc-Error
                        return;
                    }
                    else
                        QMessageBox::warning( this, "Autoprocessing WARNING", fn+EOL+rv, QMessageBox::Ok );
                }
            }
            else
            {
                qDebug() << "AutoProc-ERROR: Parameterfile not found" << fn;
                QMessageBox::critical( this, "Autoprocessing ERROR", "Parameterfile not found:\n"+fn, QMessageBox::Ok );
                autoProcessingFile->close();
                bIgnoreRecalc = false;  // Init done AutoProc-Error
                return;
            }
        }

        else if ( line.startsWith("USE QMAX FROM DATA") )
        {
            ui->radEditQmaxData->setChecked(true);
        }
        else if ( line.startsWith("USE QMAX PRESET") )
        {
            ui->radEditQmaxPreset->setChecked(true);
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
            QTimer::singleShot( 100, this, SLOT(local_butCalc_clicked()) );  // ohne speichern der Parameter im Temp-File
            // Autoproc is restarted after calculation
            return;
        }

        else if ( line.startsWith("DOFFT") )
        {
            QTimer::singleShot( 100, this, SLOT(on_butFFT_clicked()) );
            // Autoproc is restarted after calculation
            return;
        }
        else if ( line.startsWith("DOIFFT") )
        {
            QTimer::singleShot( 100, this, SLOT(on_butIFFT_clicked()) );
            // Autoproc is restarted after calculation
            return;
        }

        else if ( line.startsWith("TEST ") )
        {
            ui->cbsTestImgSelect->setCurrentIndex( line.mid(5).trimmed().toInt() );
            on_butTestGo_clicked();
            ui->cbsFFTWindows->setCurrentIndex( ui->cbsFFTWindows->count()-1 );
        }

        else if ( line.startsWith("TIMETEST") )
        {   // TIMETEST <outfilename>     starts the TimeTest procedure specified in the Tools -> Timing tests dialog
            QString fn = line.mid(8).trimmed();
            if ( fn.isEmpty() )
            {
                qDebug() << "AutoProc-ERROR: timetest file invalid";
                QMessageBox::critical( this, "Autoprocessing ERROR", "Timetest outputfile not specified.", QMessageBox::Ok );
                autoProcessingFile->close();
                bIgnoreRecalc = false;  // "TIMETEST" AutoProc-Error
                return;
            }
            performTimingTests(fn);
            // Autoproc is restarted after calculation
            return;
        }

        break; // endless loop
    } // while true
    qDebug() << "AutoProc: restart timer" << bIgnoreRecalc;
    QTimer::singleShot( 200, this, SLOT(autoProcessingTimer()) );
}


void SC_MainGUI::fitAutomaticCreateFile(dlgConfigAutoFit *dlg)
{
    if ( lastDataFile.isEmpty() )
    {
        QMessageBox::information(dlg,"Create automatic file","No datafile loaded, no automatic fit file is written.");
        return;
    }

    QFile f(dlg->getFilename());
    if ( !f.open(QIODevice::WriteOnly) )
    {
        QMessageBox::critical(dlg,"Open file error",f.fileName()+":\n"+f.errorString());
        return;
    }
    f.write("# This file is generated from the current settings." EOL);
    f.write("# Uncomment lines to activate the commands." EOL);
    f.write(EOL);
    f.write("# -- The following commands are interpreted during the first scan of the file:" EOL);
    f.write(EOL);
    f.write("# Scale: <min>, <max>  # this can be used to scale all images to the same range." EOL);
    f.write("#                      otherwise all images are scaled to the min/max of the original image." EOL);
    f.write(EOL);
    QFileInfo fi(f.fileName());
    f.write(qPrintable("GlobLog: "+fi.completeBaseName()+".log" EOL));
    f.write("# The given filename (without path) is used for the global logfile." EOL);
    f.write(EOL);
    f.write(qPrintable("Param: "+ui->lblLoadFilename->text()+EOL));
    f.write("# The parameter input file to be used during the fit as a first start." EOL);
    f.write(EOL);
    f.write("# DirMask: *.dat   # Filemask used for the DirUp / DirDown commands." EOL);
    f.write("# DirUp: <dir>     # Directory to be scanned in alphanumeric ascending order." EOL);
    f.write("# DirDown: <dir>   # Directory to be scanned in alphanumeric descending order." EOL);
    f.write(qPrintable("File: "+lastDataFile+"  # Single file to be used (multiple occurances possible)." EOL));
    f.write(EOL);
    QStringList slUse, slRest;
    for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
    {
        QString tmp = ui->tblFitValues->verticalHeaderItem(r)->text();
        tblCheckBox *tog = static_cast<tblCheckBox*>(ui->tblFitValues->cellWidget(r,0));
        if ( tog->tog()->isChecked() ) slUse << tmp; else slRest << tmp;
        tmp = "Limits: "+tmp+QString("; %1; %2  # Min; Max")
                                     .arg(ui->tblFitValues->item(r,1)->text(),
                                          ui->tblFitValues->item(r,3)->text())+EOL;
        f.write(qPrintable(tmp));
    }
    f.write(EOL);
    f.write("# -- The following commands are interpreted for every input file." EOL);
    f.write(EOL);
    f.write("# Info: <text>  # This text (up to 20 chars) will be shown in the GUI." EOL);
    f.write(EOL);
    f.write("Threads: 0  # Number of threads to use (0=GPU if available)" EOL);
    f.write(EOL);
    f.write(qPrintable(QString("Border: %1; %2; %3" EOL)
                           .arg(ui->inpFitBorder->value())
                           .arg(ui->inpFitBStop->value())
                           .arg(ui->togFitUseMask->isChecked()?1:0)));
    f.write("# Border size; Beamstop half-size; Flag(0/1) to ignore pixel less equal corner pixel" EOL);
    f.write(EOL);
    f.write(qPrintable("Use: "+slUse.join(", ")+EOL));
    f.write("# List of parameternames used during the next fit step." EOL);
    f.write(qPrintable("# Rest of parameters selected for fit: "+slRest.join(", ")+EOL));
    f.write(EOL);
    f.write(qPrintable(QString("Fit: Stp=%1; Iter=%2; Tol=%3; Diff<5; Kenn=A-@D-@F" EOL)
                           .arg(ui->inpFitStepSize->value(),0,'f',1)
                           .arg(ui->inpFitMaxIter->value())
                           .arg(ui->inpFitTolerance->text()) ));
    f.write("# Perform one fit operation. The values of Stp= (Step size), Iter= (Maximum iterations)," EOL);
    f.write("# Tol= (Tolerance) are the same as in the GUI used for this fit run. After each fit run" EOL);
    f.write("# the percentual change of each parameter is summed up and if this sum is less than the" EOL);
    f.write("# Diff= value the fit step is finished. If not, this fit run is repeated up to 20 times" EOL);
    f.write("# to avoid endless loops. Kenn= is used for the logfiles of each run, @D is the current" EOL);
    f.write("# date, @F is the current input file to fit, a runnumber is appended automatically." EOL);
    f.close();
}


void SC_MainGUI::on_butFitAutomatic_clicked()
{
    // Abfragen des Input-Files und diverser Flags
    dlgConfigAutoFit *dlgCfgAuto = new dlgConfigAutoFit( this );
    connect( dlgCfgAuto, SIGNAL(fitAutoCreateFile(dlgConfigAutoFit*)),
             this, SLOT(fitAutomaticCreateFile(dlgConfigAutoFit*)) );
    if ( dlgCfgAuto->exec() != QDialog::Accepted ) return;
    QString fnInput = dlgCfgAuto->getFilename();

    QFile finp( fnInput );
    if ( ! finp.open(QIODevice::ReadOnly) )
    {
        qDebug() << fnInput << finp.errorString();
        QMessageBox::critical( this, "Automatic Fit", finp.fileName()+"\n"+finp.errorString(), QMessageBox::Ok );
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
    QString kenn = QDate::currentDate().toString("-yyyyMMdd-");     // Als Default
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
            int posc = line.indexOf("#");
            if ( posc >= 0 ) { line.truncate(posc); line=line.trimmed(); }
            if ( line.isEmpty() ) continue;         // Ignore empty lines
            if ( firstScanOfFile ) slInputLines << line;

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
                        QMessageBox::critical( this, "Automatic Fit",
                                              "Error open global logfile:\n"+fglobLog->fileName()+"\n"+fglobLog->errorString(),
                                              QMessageBox::Ok );
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
                            QMessageBox::critical( this, "Automatic Fit",
                                                  "Error open LaTeX output file:\n"+ftex->fileName()+"\n"+ftex->errorString(),
                                                  QMessageBox::Ok );
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
                    QString rv = local_Load_all_Parameters( line.mid(6).trimmed() );
                    if ( ! rv.isEmpty() )
                    {
                        if ( fglobLog ) fglobLog->write(qPrintable(rv));
                        if ( rv.contains("Error:") )
                        {
                            QMessageBox::critical( this, "Parameterfile invalid", line.mid(6).trimmed()+EOL+rv, QMessageBox::Ok );
                            _bAbbruch = true;
                            break; // while ( ! finp.atEnd() )
                        }
                    }
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
                    if ( sl.size() != 3 ) continue;     // Es müssen immer genau "Name;Min;Max" angegeben werden
                    for ( int i=0; i<sl.size(); i++ )
                        sl[i] = sl[i].trimmed();
                    _param2fitval *p2f = fitparams;
                    if ( p2f != nullptr )
                    {
                        for ( int r=0; r<ui->tblFitValues->rowCount(); r++ )
                        {
                            QString p = ui->tblFitValues->verticalHeaderItem(r)->text();
                            if ( p == sl[0] )
                            {
                                double min = sl[1].toDouble();
                                double max = sl[2].toDouble();
                                if ( min >= max ) break;
                                ui->tblFitValues->item(r,1)->setText(sl[1]);
                                p2f->value(p)->min = min;
                                ui->tblFitValues->item(r,3)->setText(sl[2]);
                                p2f->value(p)->max = max;
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
                _param2fitval *p2f = fitparams;
                if ( p2f != nullptr ) // TODO?
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
            }

            if ( line.startsWith("Border:") )  // NEU 25.07.2024
            {   // "Border: 0; 0; 0"
                // Border size; Beamstop half-size; ignore pixel
                lastAutoFitLine += "\n" + line;
                ui->editAutoFit->appendPlainText(lastAutoFitLine);
                QStringList slCmd = line.mid(7).trimmed().split(";");
                for ( int i=0; i<slCmd.size(); i++ )
                    slCmd[i] = slCmd[i].trimmed();
                //qDebug() << slCmd;
                while ( slCmd.size() < 3 ) slCmd << "0";
                ui->inpFitBorder->setValue(slCmd[0].toInt());
                ui->inpFitBStop->setValue(slCmd[1].toInt());
                ui->togFitUseMask->setChecked(slCmd[2].toInt()!=0);
                continue;
            }

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
                if ( calcGui->updateParamValue( sl[1].trimmed(), sl[2].trimmed().toDouble(), SETCOLMARK_IGNORED, false ) )
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
                if ( calcGui->updateParamValue( sl[1].trimmed(), sl[2].trimmed().toDouble(), SETCOLMARK_IGNORED, false ) )
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
                if ( calcGui->updateParamValue( sl[0].trimmed(), sl[1].trimmed().toDouble(), SETCOLMARK_IGNORED, false ) )
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
                            stp = slCmd[i].midRef(4).toDouble();
                        else if ( slCmd[i].startsWith("Iter=") )
                            iter = slCmd[i].midRef(5).toInt();
                        else if ( slCmd[i].startsWith("Tol=") )
                            tol = slCmd[i].midRef(4).toDouble();
                        else if ( slCmd[i].startsWith("Diff<") )
                            maxDif = slCmd[i].midRef(5).toDouble();
                        else if ( slCmd[i].startsWith("Kenn=") )
                        {
                            kenn = slCmd[i].mid(5);
                            kenn = kenn.replace( "@D", QDate::currentDate().toString("yyyyMMdd") );
                            if ( usedImages.size() > 0 )
                                kenn = kenn.replace("@F", QFileInfo(usedImages.first()).baseName() );
                            else
                                kenn = kenn.replace("@F", "" );
                            if ( !kenn.endsWith("-") ) kenn += "-";  // Damit man den Zähler erkennt.
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
                                // Jetzt auch noch als _org abspeichern
                                QString tmpfn = basePath + QFileInfo(usedImages.first()).baseName() + "_org.png";
                                if ( useFixedScaling ) lastFitImage->setFixScaling( minFixedScale, maxFixedScale );
                                lastFitImage->saveImage( tmpfn );
                                if ( fglobLog ) fglobLog->write(qPrintable("Data file: "+usedImages.first()+EOL) );
                                if ( ftex != nullptr && dlgCfgAuto->isLatexOrgImg() )
                                    latexImg1 = tmpfn;
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
                            QMessageBox::critical( this, "Fit image", "No fit image found", QMessageBox::Ok );
                            return;
                        }
                        if ( lastFitImage != nullptr && lastFitImage->myHeight() != lastFitImage->myWidth() )
                        {
                            if ( fglobLog )
                            {
                                fglobLog->write(qPrintable("None square images not supported: "+usedImages.first()+EOL) );
                                fglobLog->close();
                            }
                            if ( ftex ) ftex->close();
                            ui->lblAutofitInfo->setText("None square data");
                            ui->butFitAutomatic->setEnabled(true);
                            fitIsRunning = false;
                            QMessageBox::critical( this, "Fit image", "None square images not supported" EOL+usedImages.first(), QMessageBox::Ok );
                            return;
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
                            QMessageBox::critical( this, "Fit image", "File open error" EOL+usedImages.first(), QMessageBox::Ok );
                            return;
                        }
                        if ( lastFitImage->myHeight() != lastFitImage->myWidth() )
                        {
                            if ( fglobLog )
                            {
                                fglobLog->write(qPrintable("None square images not supported: "+usedImages.first()+EOL) );
                                fglobLog->close();
                            }
                            if ( ftex ) ftex->close();
                            ui->lblAutofitInfo->setText("None square data");
                            ui->butFitAutomatic->setEnabled(true);
                            fitIsRunning = false;
                            QMessageBox::critical( this, "Fit image", "None square images not supported" EOL+usedImages.first(), QMessageBox::Ok );
                            return;
                        }
                        QString tmpfn = basePath + QFileInfo(usedImages.first()).baseName() + "_org.png";
                        if ( useFixedScaling ) lastFitImage->setFixScaling( minFixedScale, maxFixedScale );
                        lastFitImage->saveImage( tmpfn );
                        if ( fglobLog ) fglobLog->write(qPrintable("Data file: "+usedImages.first()+EOL) );

                        // Damit beim butFitStart auch das richtige Image gefunden wird, muss dieses in der
                        //  ComboBox der FitWindows als Current stehen. Anscheinend ist bei der neueren
                        //  Qt-Version etwas im Handling anders, sodass ich es hier selber setze.
                        ui->cbsFitImageWindows->setCurrentText( lastFitImage->windowTitle() );

                        // TODO: Wenn das gerade geladene Bild von den Dimensionen nicht passt,
                        // könnte der Fit-Lauf abgebrochen werden ...

                        // TODO: Aktivieren der Rechtecke zum Ausblenden des Fits...
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
                        QHash<QString,_fitLimits*>::const_iterator it = fitparams->constBegin();
                        while ( it != fitparams->constEnd() )
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
            local_butCalc_clicked(); // ohne speichern der Parameter im Temp-File

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
                        imgData[XY2IDX(calcGui->minX(),calcGui->maxX(),calcGui->minY(),calcGui->maxY(),ihex,i)] // fitAutomatic
                                = fitClass->getResiduenPixel(ihex,i);

                widImage* img = addImage( true, calcGui->minX(),calcGui->maxX(),calcGui->minY(),calcGui->maxY(), imgData, "Residuen", false );
                img->addMetaInfo( "Residuenplot", "True" );
                img->addMetaInfo( "NbRows", QString::number(calcGui->maxY() - calcGui->minY()) );
                img->addMetaInfo( "NbCols", QString::number(calcGui->maxX() - calcGui->minX()) );
                QString dbgsave;
                img->saveImageColor( resifn, tblResGlo, dbgsave );  // TODO: Returnvalue false -> save failed
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
        if ( ui->togAutoPosit->isChecked() && curFitImage != nullptr )
        {   // Vorher aber die Positionen so bestimmen, dass dieses neue Bild direkt
            // neben das bisherige (anzufittende) Bild kommt.
            imgPosX = curFitImage->pos().x() + curFitImage->width();
            imgPosY = curFitImage->pos().y();
            //qDebug() << "Pos vor Calc" << imgPosX << imgPosY << curFitImage->pos() << curFitImage->size();
        }
        local_butCalc_clicked(); // ohne speichern der Parameter im Temp-File
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
            imgData[XY2IDX(calcGui->minX(),calcGui->maxX(),calcGui->minY(),calcGui->maxY(),ihex,i)] // butShowResiduen
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
    //_param2fitval *p2f = fitparams;
    if ( fitparams == nullptr ) fitparams = new _param2fitval;
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
            _fitLimits *fl = fitparams->value(tmp[0],nullptr);
            if ( fl == nullptr )
            {
                fl = new _fitLimits;
                fitparams->insert( tmp[0], fl );
            }
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
        qDebug() << fn1 << f1.errorString();  // TODO: MessageBox::critical...
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

    qDebug() << cnt << "Differences";
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
    SETCOL( inp, SETCOLMARK_PARDIFF );
    calcGui->updateToolTipForCompare( inp, cv->par );
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
    SETCOL( inp, SETCOLMARK_PARDIFF );
    calcGui->updateToolTipForCompare( inp, cv->par );
}
void SC_MainGUI::compInt( QSpinBox *inp, QSettings &sets, QString key )
{
    int cur = inp->value();
    int par = sets.value(key,0).toInt();
    //qDebug() << "compInt" << key << cur << par << inp->toolTip();
    if ( cur == par ) return;
    _CompValues *cv = new _CompValues;
    cv->cur = QString::number(cur);
    cv->par = QString::number(par);
    compWerte.insert( key, cv );
    SETCOL( inp, SETCOLMARK_PARDIFF );
    calcGui->updateToolTipForCompare( inp, cv->par );
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
    SETCOL( inp, SETCOLMARK_PARDIFF );
    calcGui->updateToolTipForCompare( inp, cv->par );
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
    SETCOL( tog, SETCOLMARK_PARDIFF );
    calcGui->updateToolTipForCompare( tog, cv->par );
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
    SETCOL( tog, SETCOLMARK_PARDIFF );
    calcGui->updateToolTipForCompare( tog, cv->par );
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
    SETCOL( tog, SETCOLMARK_PARDIFF );
    calcGui->updateToolTipForCompare( tog, cv->par );
}

void SC_MainGUI::on_actionCompare_current_parameters_with_file_triggered()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString fn = data.value("LastParam",".").toString();
    fn = data.value("LastCompareParam",fn).toString();
    fn = QFileDialog::getOpenFileName( this, "Compare current Parameter", fn,
                                      "Parameter (*.ini *.sas_tpv)" );
                                      //nullptr, QFileDialog::DontUseNativeDialog | QFileDialog::DontResolveSymlinks );
    if ( fn.isEmpty() ) return;
    data.setValue("LastCompareParam",fn);

    on_butResetColorMarker_clicked();

    QSettings sets( fn, QSettings::IniFormat );

    // Die folgenden ObjectNames der GUI-Elemente werden in der Schleife über alle Elemente ignoriert.
    // Der Grund: diese Werte sind im Paremeterfile in einer anderen Gruppe gespeichert.
    QStringList tobeignored;
    tobeignored << ui->radQ1->objectName()
                << ui->radQ2->objectName()
                << ui->radQ4->objectName()
                << ui->togExpandImage->objectName()
                << ui->inpVAx1->objectName()
                << ui->inpAy1->objectName()
                << ui->inpAz1->objectName()
                << ui->inpVAx2_VAx1->objectName()
                << ui->inpAy2_Ay1->objectName()
                << ui->inpAz2_Az1->objectName()
                << ui->inpVAx3_VAx1->objectName()
                << ui->inpAy3_Ay1->objectName()
                << ui->inpAz3_Az1->objectName()
                << ui->inpSigX->objectName()
                << ui->inpSigY_SigX->objectName()
                << ui->inpSigZ_SigX->objectName()
                //<< ui->intHKLmax->objectName()
                //<< ui->intGridPoints->objectName()
                //<< ui->inpNumCores->objectName()
                << ui->inpBeamPosX->objectName()
                << ui->inpBeamPosY_BeamPosX->objectName();
    //sets.beginGroup( "AI" );
    //compString( ui->inpSubDir, sets, "LastSubDir", "AI:" );
    //compInt( ui->cbsAIoutputFormat, sets, "Grayscale", "AI:" );
    //compBool( ui->grpFileInput, sets, "FileInputEna", "AI:" );
    //compString( ui->inpFileName, sets, "FileInputLast", "AI:" );
    //compString( ui->inpFileClass, sets, "FileClass", "AI:" );
    //compBool( ui->togAIUseFFTOutput, sets, "GenerateIFFT", "AI:" );
    //compBool( ui->radAILinOutput, sets, "LinOut", "AI:" );
    //compBool( ui->togAIScaleOutput, sets, "ScaleOut", "AI:" );
    //sets.endGroup();

    calcGui->compareParameter( sets, compWerte, tobeignored );

    sets.beginGroup( "Inputs" );
    compBool( ui->radQ1, sets, "RadioButtonQ1" );
    compBool( ui->radQ2, sets, "RadioButtonQ2" );
    compBool( ui->radQ4, sets, "RadioButtonQ4" );
    compBool( ui->togExpandImage, sets, "ExpandImage" );
    compDouble( ui->inpVAx1, sets, "EditAxis1x" );
    compDouble( ui->inpAy1, sets, "EditAxis1y" );
    compDouble( ui->inpAz1, sets, "EditAxis1z" );
    compDouble( ui->inpVAx2_VAx1, sets, "EditAxis2x" );
    compDouble( ui->inpAy2_Ay1, sets, "EditAxis2y" );
    compDouble( ui->inpAz2_Az1, sets, "EditAxis2z" );
    compDouble( ui->inpVAx3_VAx1, sets, "EditAxis3x" );
    compDouble( ui->inpAy3_Ay1, sets, "EditAxis3y" );
    compDouble( ui->inpAz3_Az1, sets, "EditAxis3z" );
    compDouble( ui->inpSigX, sets, "Editdom1" );
    compDouble( ui->inpSigY_SigX, sets, "Editdom2" );
    compDouble( ui->inpSigZ_SigX, sets, "Editdom3" );
    compInt( ui->intHKLmax, sets, "HKLmax" );               // TODO
    compInt( ui->intGridPoints, sets, "GridPoints" );       // TODO
    //compInt( ui->inpNumCores, sets, "Threads" ); - macht hier weniger Sinn...
    compDouble( ui->inpBeamPosX, sets, "EditCenterX" );
    compDouble( ui->inpBeamPosY_BeamPosX, sets, "EditCenterY" );
    sets.endGroup();
    //sets.beginGroup( "AI" );
    //compString( ui->inpSubDir, sets, "LastSubDir", "AI:" );
    //compInt( ui->cbsAIoutputFormat, sets, "Grayscale", "AI:" );
    //compBool( ui->grpFileInput, sets, "FileInputEna", "AI:" );
    //compString( ui->inpFileName, sets, "FileInputLast", "AI:" );
    //compString( ui->inpFileClass, sets, "FileClass", "AI:" );
    //compBool( ui->togAIUseFFTOutput, sets, "GenerateIFFT", "AI:" );
    //compBool( ui->radAILinOutput, sets, "LinOut", "AI:" );
    //compBool( ui->togAIScaleOutput, sets, "ScaleOut", "AI:" );
    //sets.endGroup();

    showColorMarker( SETCOLMARK_PARDIFF );
/*
    if ( compWerte.size() > 0 )
    {
        QFile fout( fn+".diff.csv" );
        if ( ! fout.open(QIODevice::WriteOnly) )
            qDebug() << fout.fileName() << fout.errorString();
        else
            fout.write(qPrintable("Key ; Current ; Parameter from file "+fn+EOL));
        // Auch bei einem Fehler beim Öffnen der Datei werden die Ergebnisse auf der Console angezeigt
        qDebug() << "Key" << "Cur" << fn;
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
    else
        QMessageBox::information( this, "Compare parameters", "The given file contains no other values.", QMessageBox::Ok );
*/
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
    ui->lblParamSearchLog->hide();
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
    ui->lblParamSearchLog->hide();
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
        if ( fn.endsWith("ParamSearch.log") ) continue;
        if ( fn.contains("tmpForAI") ) continue;
        if ( fn.contains("TempParamSave") ) continue;
        if ( all )
        {
            if ( fn.endsWith(".ini") ) continue;
            if ( fn.endsWith(".png") ) continue;
            if ( fn.contains("training table") ) continue;
            /*if ( fn.endsWith(".spr") )
            {   // Spezielle Abfrage, damit von den vielen *.spr Files nicht alle verwendet werden
                if ( fn.lastIndexOf('.') - fn.lastIndexOf('_') > 3 ) continue;
            }*/
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
    ui->lblParamSearchLog->hide();
    ui->butParamSearchGenerate->setEnabled(false);
    qApp->setOverrideCursor(Qt::WaitCursor);
    searchParameterFilesHelper(ui->inpParamSearchPath->text(),bl,"*.ini");
    searchParameterFilesHelper(ui->inpParamSearchPath->text(),bl,"*.sas_tpv");
    qApp->restoreOverrideCursor();
    ui->lblParamSearchCount->setText( QString::number(ui->lisParamSearchResult->count()) );
}

void SC_MainGUI::on_butDataSearchDoit_clicked()
{
    if ( bSearchParamsRunning ) return;
    ui->lisParamSearchResult->setEnabled(true);
    ui->lisParamSearchResult->clear();
    int bl = ui->inpParamSearchPath->text().length();
    if ( !ui->inpParamSearchPath->text().endsWith("/") && !ui->inpParamSearchPath->text().endsWith("\\") ) bl++;
    ui->lblParamSearchLog->hide();
    ui->butParamSearchGenerate->setEnabled(false);
    qApp->setOverrideCursor(Qt::WaitCursor);
    searchParameterFilesHelper(ui->inpParamSearchPath->text(),bl,"*.*");
    qApp->restoreOverrideCursor();
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
    ui->lblParamSearchLog->setText("See <a href=\"file:///"+base+"ParamSearch.log\">"+base+"ParamSearch.log</a>");
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
        if ( fn.endsWith(".ini") || fn.endsWith(".sas_tpv") )
        {   // Parameterfiles
            if ( fn.endsWith("ConfigParams.ini") ) continue;
            flog.write( qPrintable("Read "+fn+EOL) );
            QString rv = local_Load_all_Parameters(fn);
            if ( !rv.isEmpty() )
            {
                QStringList tmp = rv.split(EOL,Qt::SkipEmptyParts);
                foreach ( QString s, tmp )
                    flog.write( qPrintable("  "+s+EOL) );
                if ( rv.contains("Error:") ) continue; // Jetzt kein Bild rechnen (käme nur Blödsinn raus)
            }
            flog.flush();
            local_butCalc_clicked(); // ohne speichern der Parameter im Temp-File
            if ( _bAbbruch )
            {
                flog.write("Aborted by user" EOL);
                break;
                // Debug: "C:/SimLab/sas-crystal/20220608 - 2D Fits und Fiber-Pattern/TestParamsCylinder.ini"
            }
            fn.replace(".ini",".png");
            fn.replace(".sas_tpv",".png");
            flog.write( qPrintable(QString("  PrepTime %1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimPrep),12,'f',3)+EOL) );
            flog.write( qPrintable(QString("  CalcTime %1 ms").arg(calcGui->higResTimerElapsed(SC_CalcGUI::htimCalc),12,'f',3)+EOL) );
            // Loggen der bislang definierten Unterscheidungsmerkmale der Kernel
            static QString part2str[] = {"partSphere", "partCylinder", "partDisk", "partVesicle", "partCube", "partEllips",
                                         "partTriaxEllips", "partSuperEllips", "partSuperball", "partChain", "partKPChain"};
            static QString peak2str[] = {"peakLorentzian", "peakGauss", "peakMod1Lor", "peakMod2Lor", "peakPseudoVoigt",
                                         "peakPearsonVII", "peakGamma", "peakAnisoGauss"};
            static QString latt2str[] = {"lattLamellae", "lattHexPackCyl", "lattSquarePackCyl", "lattRectCentCyl", "lattBCC",
                                         "lattFCC", "lattHCP", "lattSC", "lattBCT", "lattGyroid", "lattOBDD", "lattPlumNightmare",
                                         "lattNone", "lattCPLayers", "latt2DHexG", "latt2DSquareG", "latt1DLamG", "lattFd3m",
                                         "lattOrthorombic", "lattLDQ12", "lattPercYevick", "lattTeubnerStrey", "lattPm3n",
                                         "lattP42mnm", "lattFdddNetwork"};
            static QString ordi2str[] = {"ordisGauss", "ordisExponent", "ordisOnsager", "ordisMaierSaupe", "ordisCutOff",
                                         "ordisLaguerre", "ordisZDir", "ordisIsotropic", "ordisMirrorGauss", "ordisMirrorExponent",
                                         "ordisMirrorOnsager", "ordisMirrorMaierSaupe", "ordisMirrorCutOff", "ordisFiberPattern"};
            QString latt = calcGui->params["LType"]->value.number < 0 ? "lattINVALID" : latt2str[(int)(calcGui->params["LType"]->value.number)];
            QString ordi = calcGui->params["Ordis"]->value.number < 0 ? "ordisINVALID" : ordi2str[(int)(calcGui->params["Ordis"]->value.number)];

            flog.write( qPrintable(QString("  Kernel = %1_%2_%3_%4")
                                      .arg(part2str[(int)(calcGui->params["ComboBoxParticle"]->value.number)])
                                      .arg(peak2str[(int)(calcGui->params["ComboBoxPeak"]->value.number)])
                                      .arg(latt)
                                      .arg(ordi)
                                  +EOL) );
            //Debug: ("LType", "ComboBoxPeak", "ucn2", "CenterBeam", "EditSigma", "RadioButtonPara", "EditQmaxData",
            //      "RotAlpha", "EditAzi", "EditPeakPar", "EditWavelength", "CalcQmax", "SigY", "SigmaL", "CheckBoxWAXS",
            //      "ucalpha", "EditRadiusi", "CheckBoxTwinned", "ucgamma", "VAx1", "SigX", "EditPixelNoY", "EditDebyeWaller",
            //      "Alpha", "ucn3", "VAx2", "EditPixelY", "Base", "ComboBoxParticle", "ucc", "iso", "ucb", "ucn1",
            //      "CenterMidpoint", "EditPixelNoX", "Ordis", "EditQmax", "Az2", "EditRadius", "ifluc", "EditDbeta",
            //      "ucbeta", "EditCeffcyl", "SigZ", "acpl", "theta", "EditPixelX", "EditBFactor", "Az3", "bcpl", "I0",
            //      "EditRho", "reff", "rotPhi", "ComboBoxInterior", "VAx3", "Ay1", "rotTheta", "uca", "rfluc", "EditDet",
            //      "P1", "RadButDebyeScherrer", "phi", "ucpsi", "Length", "EditQmaxPreset", "EditRelDis", "BeamPosY",
            //      "Ay2", "GridPoints", "Az1", "BeamPosX", "Ay3", "HKLmax", "EditDomainSize", "EditDist", "EditCeff")
        }
        else
        {   // Datenfiles
            if ( ! local_OpenMeasFile( fn, &lastUsedImage ) ) continue;
            flog.write( qPrintable("Read "+fn+EOL) );
            if ( fn.endsWith(".spr") ) lastUsedImage->setLogScaling( false );
            fn.append(".png");
        }
        if ( _bAbortTimer )
        {
            flog.write("Aborted due to timeout in calculation" EOL EOL);
        }
        else
        {
            flog.write( qPrintable("Save "+fn+EOL+EOL) );
            lastUsedImage->saveImage(fn);
            lastUsedImage->close();
        }
        flog.flush();
    }
    flog.close();
    ui->lblParamSearchLog->show();
    bSearchParamsRunning = false;
    local_Load_all_Parameters(myLocPar);    // Wiederherstellen der vorherigen Parameterwerte
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
        QMessageBox::critical( this, "TPV save file", f.fileName()+"\n"+f.errorString(), QMessageBox::Ok );
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
    sets.setValue("val_addLinesHwidth",ui->inpTPVaddLinesHwidth->value());
    sets.setValue("val_addLinesVwidth",ui->inpTPVaddLinesVwidth->value());
    sets.setValue("ena_addBS",ui->togTPVaddBS->isChecked());
    sets.setValue("ena_addNoise",ui->togTPVaddNoise->isChecked());
    sets.setValue("ena_convolute",ui->togTPVconvolute->isChecked());
    sets.setValue("ena_calcRphi",ui->togTPVcalcRphi->isChecked());
    sets.setValue("ena_calcFFT",ui->togTPVcalcFFT->isChecked());
    sets.setValue("ena_saveExtra",ui->togTPVsaveBaseExtra->isChecked());
    sets.setValue("ena_generatePNG",ui->togTPVgeneratePNG->isChecked());
    sets.setValue("ena_scaleScat",ui->togTPVscaleScat->isChecked());
    sets.setValue("val_numimg",ui->inpTPVnumimg->value());
    sets.setValue("val_outPath",ui->inpTPVoutPath->text());
    sets.sync();
    // In dieses File kommen jetzt auch alle Daten
    // Groups [FCC%20Spheres] und [Inputs]
    performSaveParamOperation( fn );
    // FFT Tab
    sets.beginGroup("FFT");
    sets.setValue( "FFTLinInput", ui->radFFTLinInput->isChecked() );
    sets.setValue( "FFTScaleImg", ui->togFFTscaleImg->isChecked() );
    sets.setValue( "FFTclipImg", ui->togFFTclipImg->isChecked() );
    sets.setValue( "FFTclip40Img", ui->togFFTclip40Img->isChecked() );
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
    sets.setValue( "FFTcutOutX", ui->inpFFTcutX->value() );
    sets.setValue( "FFTcutOutY", ui->inpFFTcutY->value() );
    sets.setValue( "FFTcutOutEna", ui->togFFTcutEnabled->isChecked() );
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
    on_butResetColorMarker_clicked();
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
    ui->inpTPVaddLinesHwidth->setValue( sets.value("val_addLinesHwidth",2).toInt() );
    ui->inpTPVaddLinesVwidth->setValue( sets.value("val_addLinesVwidth",2).toInt() );
    ui->inpTPVaddLinesH->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesV->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesHwidth->setEnabled(ui->togTPVaddLines->isChecked());
    ui->inpTPVaddLinesVwidth->setEnabled(ui->togTPVaddLines->isChecked());
    ui->lblTPVaddLines->setEnabled(ui->togTPVaddLines->isChecked());
    ui->togTPVaddNoise->setChecked( sets.value("ena_addNoise",true).toBool() );
    ui->togTPVconvolute->setChecked( sets.value("ena_convolute",true).toBool() );
    ui->togTPVcalcRphi->setChecked( sets.value("ena_calcRphi",true).toBool() );
    ui->togTPVcalcFFT->setChecked( sets.value("ena_calcFFT",true).toBool() );
    ui->togTPVsaveBaseExtra->setChecked( sets.value("ena_saveExtra",false).toBool() );
    ui->togTPVgeneratePNG->setChecked( sets.value("ena_generatePNG",true).toBool() );
    ui->togTPVscaleScat->setChecked( sets.value("ena_scaleScat",false).toBool() );
    ui->inpTPVnumimg->setValue( sets.value("val_numimg",10).toInt() );
    ui->inpTPVoutPath->setText( sets.value("val_outPath",".").toString() );
    // Alle Berechnungdaten
    QString rv = local_Load_all_Parameters(fn);
    if ( ! rv.isEmpty() )
    {
        if ( rv.contains("Error:") )
            QMessageBox::critical( this, "Parameterfile ERROR", fn+EOL+rv, QMessageBox::Ok );
        else
            QMessageBox::warning( this, "Parameterfile WARNING", fn+EOL+rv, QMessageBox::Ok );
        // Trotzdem den Rest noch laden
    }
    // FFT Tab
    sets.beginGroup("FFT");
    ui->radFFTLinInput->setChecked(  sets.value("FFTLinInput",true).toBool() );
    ui->radFFTLogInput->setChecked( !sets.value("FFTLinInput",true).toBool() );
    ui->togFFTscaleImg->setChecked( sets.value("FFTScaleImg",true).toBool() );
    ui->togFFTclipImg->setChecked( sets.value("FFTclipImg",false).toBool() );
    ui->togFFTclip40Img->setChecked( sets.value("FFTclip40Img",false).toBool() );
    ui->grpFFTuseRphi->setChecked( sets.value("FFTuseRphi",true).toBool() );
    ui->cbsFFTsizeRphi->setCurrentIndex( sets.value("FFTsizeRphi",2).toInt() );
    ui->togFFTscaleRphi->setChecked( sets.value("FFTScaleRphi",true).toBool() );
    ui->togFFTclipRphi->setChecked( sets.value("FFTclipRphi",false).toBool() );
    ui->togFFTclip40Rphi->setChecked( sets.value("FFTclip40Rphi",false).toBool() );
    ui->togFFTdispRphi->setChecked( sets.value("DispRphi",true).toBool() );
    ui->togFFTscaleOut->setChecked( sets.value("FFTScaleOutput",true).toBool() );
    ui->togFFTclipOut->setChecked( sets.value("FFTclipOutput",false).toBool() );
    ui->togFFTclip40Out->setChecked( sets.value("FFTclip40Output",false).toBool() );
    ui->togIFFTSwap->setChecked( sets.value("FFTSwapOutput",true).toBool() );
    ui->cbsFFTsizeOut->setCurrentIndex( sets.value("FFTsizeOut",2).toInt() );
    ui->inpFFTcutX->setValue( sets.value("FFTcutOutX",64).toInt() );
    ui->inpFFTcutY->setValue( sets.value("FFTcutOutY",64).toInt() );
    ui->togFFTcutEnabled->setChecked( sets.value("FFTcutOutEna",false).toBool() );
    ui->inpFFTcutX->setEnabled( ui->togFFTcutEnabled->isChecked() );
    ui->inpFFTcutY->setEnabled( ui->togFFTcutEnabled->isChecked() );
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
    sets.endGroup();
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
"training table.txt":
    28    1000 Start pixel = (       1       1)
  run_no.   intensity radius    radius_m  length    sigma_r   sigma_l   rho       a_cell    b_cell    c_cell   psi      phi      width    aziwidth dw       dbeta    base     x0       y0       nbeam    mbeam    linex    liney    points   holder   hold_ang bs_shape det_shape
  0 1.61E+04 1.57E+01 2.68E+00 2.86E+01 6.39E-02 1.10E-01 8.40E-02 4.77E+01 4.77E+01 2.62E+01 0.00E+00 8.10E+01 1.44E-02 1.11E-02 1.00E+00 4.26E-01 9.34E-03 9.22E+01 5.95E+01 3.00E+00 7.00E+00 1.38E+02 -6.40E+01 6.40E+01 false 1.46E+02 rect sq
*/
void SC_MainGUI::on_butTPVreadTrainingTable_clicked()
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
        QMessageBox::critical( this, "TPV Read training table", fin.fileName()+"\n"+fin.errorString(), QMessageBox::Ok );
        return;
    }
    on_butResetColorMarker_clicked();

    QString line;
    /*line =*/ fin.readLine();  // 1. Zeile, die Werte dort brauche ich nicht
    line = fin.readLine();  // 2. Zeile: Schlüsselworte
    QStringList keys = line.trimmed().split(" ",Qt::SkipEmptyParts);
    //qDebug() << keys << keys.size();
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

    // Jetzt einmal die maximale Run-No bestimmen und den User fragen, welchen Run er laden möchte
    int maxrun = 0;
    while ( !fin.atEnd() )
    {
        line = fin.readLine();
        if ( line.isEmpty() ) continue;
        QStringList svals = line.trimmed().split(" ",Qt::SkipEmptyParts);
        if ( svals.size() < 2 ) continue;
        if ( maxrun < svals[0].toInt() ) maxrun = svals[0].toInt();
    }
    fin.seek(0);
    /*line =*/ fin.readLine();  // 1. Zeile (Dimensionen), die Werte dort brauche ich nicht
    /*line =*/ fin.readLine();  // 2. Zeile (Schlüssel), die Werte dort brauche ich jetzt nicht mehr

    // Um einen besseren Startwert zu haben, versuche ich den LType anhand des Verzeichnis-Namens zu ermitteln
    QString mm = QFileInfo(fin).absoluteDir().dirName();
    int pos = mm.indexOf(QRegExp("[_ ]"));
    //qDebug() << "TPV base" << fn << mm << pos;
    if ( pos > 0 ) mm.truncate(pos);
    // Suchen des richtigen LType-Wertes durch scan der Combobox-Texte
    QString ltypeval = "";
    int ltypeind = -1;
    for ( int val=0; val<ui->cbsLType->count(); val++ )
        if ( ui->cbsLType->itemText(val).startsWith(mm,Qt::CaseInsensitive) )
        {
            ltypeval = ui->cbsLType->itemText(val);
            ltypeind = val;
            break;
        }

    QDialog *dlg = new QDialog(this);
    QCheckBox *tog = nullptr;
    if ( ltypeind >= 0 )
    {
        tog = new QCheckBox("Set LType to: "+ltypeval);
        tog->setChecked(true);
    }
    QLabel *lbl = new QLabel(QString("Enter run no to load (0 - %1):").arg(maxrun));
    QSpinBox *inp = new QSpinBox();
    inp->setRange( 0, maxrun );
    QPushButton *but = new QPushButton("Load");
    QVBoxLayout *lay = new QVBoxLayout();
    lay->addWidget(lbl);
    lay->addWidget(inp);
    if ( ltypeind >= 0 )
        lay->addWidget(tog);
    lay->addWidget(but);
    dlg->setLayout(lay);
    connect( but, SIGNAL(clicked()), dlg, SLOT(accept()) );
    dlg->exec();

    int usedRun = inp->value();
    //qDebug() << usedRun;
    if ( ltypeind >= 0 )
        if ( tog->isChecked() )
            calcGui->updateParamValue( "LType", ltypeind, SETCOLMARK_IGNORED, false );

    bool pixelSizeChanged = tpvParamsKeys.contains("EditPixelX");
    double defPixX = calcGui->currentParamValueDbl("EditPixelX");
    double defPixY = calcGui->currentParamValueDbl("EditPixelY");
    int zmaxt = ui->intGridPoints->value();
    bool updFlag = bIgnoreRecalc;   // Während Parameterupdates macht hier ein Auto-Recalc keinen Sinn
    bIgnoreRecalc = true;
    while ( !fin.atEnd() )
    {
        line = fin.readLine();
        if ( line.isEmpty() ) continue;
        QStringList svals = line.trimmed().split(" ",Qt::SkipEmptyParts);
        if ( svals.size() < keys.size() ) continue;     // Weniger Werte als Schlüssel werden ignoiert
        if ( usedRun != svals[0].toInt() ) continue;
        for ( int i=0; i<usedKeys.size(); i++ )
        {
            double d = svals[i].toDouble();
            if ( usedKeys[i][0] == '@' )
            {
                if ( usedKeys[i] == "@BSX" )
                {
                    ui->inpBeamPosX->setValue( d );
                    SETCOL( ui->inpBeamPosX, SETCOLMARK_TRAINTBL );
                }
                else if ( usedKeys[i] == "@BSY" )
                {
                    ui->inpBeamPosY_BeamPosX->setValue( d );
                    SETCOL( ui->inpBeamPosY_BeamPosX, SETCOLMARK_TRAINTBL );
                }
                else if ( usedKeys[i] == "@GP" )
                {
                    ui->intGridPoints->setValue( d );
                    SETCOL( ui->intGridPoints, SETCOLMARK_TRAINTBL );
                }
            }
            else if ( usedKeys[i][0] != '-' )
                calcGui->updateParamValue( usedKeys[i], d, SETCOLMARK_TRAINTBL, false );
        }
        if ( ui->inpBeamPosX->value() > 0 || ui->inpBeamPosY_BeamPosX->value() > 0 )
        {
            ui->inpBeamPosX->setValue( ui->inpBeamPosX->value() - ui->intGridPoints->value() );
            ui->inpBeamPosY_BeamPosX->setValue( ui->inpBeamPosY_BeamPosX->value() - ui->intGridPoints->value() );
            SETCOL( ui->inpBeamPosX, SETCOLMARK_TRAINTBL );
            SETCOL( ui->inpBeamPosY_BeamPosX, SETCOLMARK_TRAINTBL );
        }
        // Es wird nur eine Zeile mit Datenwerten genutzt.
        break;
    }
    fin.close();

    if ( !pixelSizeChanged )
    {
        int zmax = ui->intGridPoints->value();

        calcGui->updateParamValue( "EditPixelNoX", 2*zmax+1, SETCOLMARK_TRAINTBL );
        calcGui->updateParamValue( "EditPixelNoY", 2*zmax+1, SETCOLMARK_TRAINTBL );

        if ( fabs( (2*zmaxt+1)/(2*zmax+1) - 1.0 ) > 0.001 )
        {
            calcGui->updateParamValue( "EditPixelX", defPixX*(2*zmaxt+1)/(2*zmax+1), SETCOLMARK_TRAINTBL );
            calcGui->updateParamValue( "EditPixelY", defPixY*(2*zmaxt+1)/(2*zmax+1), SETCOLMARK_TRAINTBL );
        }
    }
    showColorMarker( SETCOLMARK_TRAINTBL );

    // Jetzt wird noch der gewünschte "scat" Datensatz aus dem verzeichnis mit der Tabelle geladen.
    QDir d(QFileInfo(fn).path());
    //qDebug() << fn;
    d.cd("scat");
    //qDebug() << d.absoluteFilePath(QString("scat_%1.spr").arg(usedRun));
    widImage *img;
    if ( ! local_OpenMeasFile(d.absoluteFilePath(QString("scat_%1.spr").arg(usedRun)),&img) )
    {
        qDebug() << "Open failed.";
        QMessageBox::critical( this, "TPV Read training table", "Unable to open "+d.absoluteFilePath(QString("scat_%1.spr").arg(usedRun)),
                              QMessageBox::Ok );
    }
    else
    {
        img->setLogScaling(false);
        QHash<QString,paramHelper*>::iterator ip = calcGui->params.begin();
        while ( ip != calcGui->params.end() )
        {
            if ( usedKeys.contains(ip.key()) )
                img->addMetaInfo( ip.key(), calcGui->currentParamValueStr( ip.key(), true ) );
            ++ip;
        }
        img->addMetaInfo( QString("from TrainTbl[%1]").arg(usedRun), fn );
        img->addMetaInfo( "@", "" );    // Sortieren
    }
    bIgnoreRecalc = updFlag;
}

void SC_MainGUI::on_actionFind_parameters_changing_image_triggered()
{
    //qDebug() << "TPV" << tpvParamsKeys;
    QStringList all = calcGui->paramsForMethod(true,false,false);

    foreach ( QString s, tpvParamsKeys )
        all.removeOne(s);
    all.removeOne("RadioButtonQ1");
    all.removeOne("RadioButtonQ2");
    all.removeOne("RadioButtonQ4");
    all.removeOne("ExpandImage");
    all.removeOne("Edithklmax");
    all.removeOne("HKLmax");
    all.removeOne("GridPoints");
    all.removeOne("LATTrows");
    all.removeOne("LATTcols");
    all.removeOne("LATTceny");
    all.removeOne("LATTcenx");
    all.removeOne("BeamPosX");
    all.removeOne("BeamPosY");
    all.removeOne("BeamcenterX");
    all.removeOne("BeamcenterY");

    all.removeOne("EditDet");       // Detektor-Abstand lassen wir hier mal weg
    all.removeOne("EditPixelX");    // Pixelsizes lassen wir hier auch weg
    all.removeOne("EditPixelY");

    all.removeOne("iso");           // Diese Werte ändern auch immer die Daten
    all.removeOne("I0");
    all.removeOne("Base");

    // Jetzt sind in 'all' nur noch Parameter enthalten, die nicht vom obigen Lesen gesetzt wurden.
    // Diese werden jetzt durchgespielt, ob diese den Datensatz ändern.

    //on_butResetColorMarker_clicked();

    // Zum Test schauen, ob der Datensatz nicht zu groß wird, das würde die Rechenzeit bei dem
    //  Durchgehen der Parameter unnötig erhöhen.
    /*if ( ui->inpGridPoints->value() > 100 )
    {
        ui->inpGridPoints->setValue(64);
        // Dann aber auch den Beamstop mittig setzen
        ui->inpBCenterX->setValue(0);
        ui->inpBCenterY->setValue(0);
    }*/

    // Erzeugen des aktuellen Datensatzes
    local_butCalc_clicked(); // ohne speichern der Parameter im Temp-File
    int mx = ui->intGridPoints->value();
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
        double curval = calcGui->currentParamValueDbl(s);
        double min, max;
        bool count; // uninteressant
        if ( ! calcGui->limitsOfParamValue( s, min, max, count ) )
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
        calcGui->updateParamValue( s, curval, SETCOLMARK_IGNORED );
        if ( _bAbbruch ) break;
    }
    foreach ( QString s, modParams )
        calcGui->updateParamValueColor( s, SETCOLMARK_CHGIMG );
    showColorMarker( SETCOLMARK_CHGIMG );

    qDebug() << modParams;
}

bool SC_MainGUI::vergleicheCheckData( QString par, double val, const double *check, int mx )
{
    calcGui->updateParamValue( par, val, SETCOLMARK_IGNORED );

    local_butCalc_clicked(); // ohne speichern der Parameter im Temp-File
    if ( _bAbbruch ) return false; // Keine Änderungsmarkierung bei Abbruch

    const double *tmp = calcGui->data();

    for ( int i=0; i<mx; i++ )
        if ( fabs( *(tmp++) - *(check++) ) > eps4 )
            return true;

    return false;
}


void SC_MainGUI::on_butResetColorMarker_clicked()
{
    DT( qDebug() << "on_butResetColorMarker_clicked()" );
    showColorMarker( SETCOLMARK_CLEARED );
    calcGui->resetParamColorMarker( SETCOLMARK_CLEARED );
    SETCOL( ui->radQ1,             SETCOLMARK_CLEARED );
    SETCOL( ui->radQ2,             SETCOLMARK_CLEARED );
    SETCOL( ui->radQ4,             SETCOLMARK_CLEARED );
    SETCOL( ui->togExpandImage,    SETCOLMARK_CLEARED );
    SETCOL( ui->inpVAx1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpAy1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpAz1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpVAx2_VAx1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpAy2_Ay1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpAz2_Az1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpVAx3_VAx1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpAy3_Ay1,            SETCOLMARK_CLEARED );
    SETCOL( ui->inpAz3_Az1,            SETCOLMARK_CLEARED );
    //SETCOL( ui->inpN1,             SETCOLMARK_CLEARED );
    //SETCOL( ui->inpN2,             SETCOLMARK_CLEARED );
    //SETCOL( ui->inpN3,             SETCOLMARK_CLEARED );
    SETCOL( ui->inpSigX,           SETCOLMARK_CLEARED );
    SETCOL( ui->inpSigY_SigX,           SETCOLMARK_CLEARED );
    SETCOL( ui->inpSigZ_SigX,           SETCOLMARK_CLEARED );
    //SETCOL( ui->inpU1,             SETCOLMARK_CLEARED );
    //SETCOL( ui->inpU2,             SETCOLMARK_CLEARED );
    //SETCOL( ui->inpU3,             SETCOLMARK_CLEARED );
    //SETCOL( ui->inpV1,             SETCOLMARK_CLEARED );
    //SETCOL( ui->inpV2,             SETCOLMARK_CLEARED );
    //SETCOL( ui->inpV3,             SETCOLMARK_CLEARED );
    SETCOL( ui->intHKLmax,         SETCOLMARK_CLEARED );
    SETCOL( ui->intGridPoints,     SETCOLMARK_CLEARED );
    SETCOL( ui->inpNumCores,       SETCOLMARK_CLEARED );
    SETCOL( ui->inpBeamPosX,       SETCOLMARK_CLEARED );
    SETCOL( ui->inpBeamPosY_BeamPosX,       SETCOLMARK_CLEARED );
    SETCOL( ui->inpSubDir,         SETCOLMARK_CLEARED );
    SETCOL( ui->cbsAIoutputFormat, SETCOLMARK_CLEARED );
    SETCOL( ui->grpFileInput,      SETCOLMARK_CLEARED );
    SETCOL( ui->inpFileName,       SETCOLMARK_CLEARED );
    SETCOL( ui->inpFileClass,      SETCOLMARK_CLEARED );
    SETCOL( ui->togAIUseFFTOutput, SETCOLMARK_CLEARED );
    SETCOL( ui->radAILinOutput,    SETCOLMARK_CLEARED );
    SETCOL( ui->radAILogOutput,    SETCOLMARK_CLEARED );
    SETCOL( ui->togAIScaleOutput,  SETCOLMARK_CLEARED );
}

void SC_MainGUI::on_butTPVupdateAllFiles_clicked()
{
    QSettings data(SETT_APP,SETT_GUI);
    data.beginGroup("TPV");
    QString fn = data.value("LastSaveFile",".").toString();
    data.endGroup();
    fn = QFileDialog::getExistingDirectory( this, "Update TPV infos", fn );
    if ( fn.isEmpty() ) return;

    QDir d(fn);
    QFileInfoList fil = d.entryInfoList( QStringList()<<"*.sas_tpv" );

    foreach ( QFileInfo fi, fil )
    {
        qDebug() << "Update" << fi.absoluteFilePath();

        QSettings sets(fi.absoluteFilePath(),QSettings::IniFormat);
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
        sets.setValue("val_addLinesHwidth",ui->inpTPVaddLinesHwidth->value());
        sets.setValue("val_addLinesVwidth",ui->inpTPVaddLinesVwidth->value());
        sets.setValue("ena_addBS",ui->togTPVaddBS->isChecked());
        sets.setValue("ena_addNoise",ui->togTPVaddNoise->isChecked());
        sets.setValue("ena_convolute",ui->togTPVconvolute->isChecked());
        sets.setValue("ena_calcRphi",ui->togTPVcalcRphi->isChecked());
        sets.setValue("ena_calcFFT",ui->togTPVcalcFFT->isChecked());
        sets.setValue("ena_saveExtra",ui->togTPVsaveBaseExtra->isChecked());
        sets.setValue("ena_generatePNG",ui->togTPVgeneratePNG->isChecked());
        sets.setValue("ena_scaleScat",ui->togTPVscaleScat->isChecked());
        sets.setValue("val_numimg",ui->inpTPVnumimg->value());
        // Nur der Pfad bleibt unberührt !
        //sets.setValue("val_outPath",ui->inpTPVoutPath->text());
        sets.sync();
    }
}


void SC_MainGUI::showColorMarker( QColor c )
{
    DT( qDebug() << "showColorMarker()" << c );
    if ( c == SETCOLMARK_CLEARED )
    {
        ui->lblColMarkChgImg->setVisible(false);
        ui->lblColMarkParDiff->setVisible(false);
        ui->lblColMarkTrainTbl->setVisible(false);
        ui->lblColMarkOldPar->setVisible(false);
        ui->grpColorMarker->setVisible(false);
        return;
    }
    if ( c == SETCOLMARK_CHGIMG   ) ui->lblColMarkChgImg->setVisible(true);
    if ( c == SETCOLMARK_PARDIFF  ) ui->lblColMarkParDiff->setVisible(true);
    if ( c == SETCOLMARK_TRAINTBL ) ui->lblColMarkTrainTbl->setVisible(true);
    if ( c == SETCOLMARK_OLDPARAM ) ui->lblColMarkOldPar->setVisible(true);
    ui->grpColorMarker->setVisible(true);
}


// Im Vorgriff auf die komplette Neustrukturierungen (nur noch die GENERAL-Klasse)
// werden hier mal Signale verwendet.

void SC_MainGUI::on_cbsLType_currentIndexChanged(int index)
{
#ifndef NOCBSCALLBACK
    DT( qDebug() << "on_cbsLType_currentIndexChanged()" << index << bIgnoreRecalc );
    bool updFlag = bIgnoreRecalc;
    bIgnoreRecalc = true;
    //Z=36496 = ComboBoxLatticeChange
    myGuiParam *gp;
    switch ( index )
    {
    case  0:  /*  Lamellae  */  //Z=36497
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 2 );
        gp=myGuiParam::getGuiParam("Length");
        if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSys1DClick(Sender);  //Z=43426
        /*
            ComboBoxCubic.ItemIndex = -1;  //Z=41593
            ComboBoxHexagonal.ItemIndex = -1;  //Z=41594
            ComboBoxTrigonal.ItemIndex = -1;  //Z=41595
            ComboBoxTetragonal.ItemIndex = -1;  //Z=41596
            ComboBoxOrthorhombic.ItemIndex = -1;  //Z=41597
            ComboBoxMonoclinic.ItemIndex = -1;  //Z=41598
            ComboBoxTriclinic.ItemIndex = -1;  //Z=41599
            ComboBox2D.ItemIndex = -1;  //Z=41600
            ComboBox1D.ItemIndex = 0;  //Z=41601

            TrackBarCellA.visible = true;  //Z=41603
            EditCellA.visible = true;  //Z=41604
            EditAmax.visible = true;  //Z=41605
            TrackBarCellB.visible = false;  //Z=41606
            EditCellB.visible = false;  //Z=41607
            EditBmax.visible = false;  //Z=41608
            TrackBarCellC.visible = false;  //Z=41609
            EditCellC.visible = false;  //Z=41610
            EditCmax.visible = false;  //Z=41611
            TrackBarCellAlpha.visible = false;  //Z=41612
            EditCellAlpha.visible = false;  //Z=41613
            EditCellAlpha.Text = '90';  //Z=41614
            TrackBarCellBeta.visible = false;  //Z=41615
            EditCellBeta.visible = false;  //Z=41616
            EditCellBeta.Text = '90';  //Z=41617
            TrackBarCellGamma.visible = true;  //Z=41618
            EditCellGamma.visible = true;  //Z=41619
            EditCellGamma.Text = '90';  //Z=41620
        */
        break;

    case  1:  /*  hexagonal cylinders  */  //Z=36544
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 1 );
        gp=myGuiParam::getGuiParam("Length");
        if(gp)gp->inp()->setValue( 500 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true );
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSys2DClick(Sender);  //Z=43431
        break;

    case  2:  /*  square cylinders  */  //Z=36592
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 1 );
        gp=myGuiParam::getGuiParam("Length");
        if(gp)gp->inp()->setValue( 500 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSys2DClick(Sender);  //Z=43436
        break;

    case  3:  /*  rectangular centered cylinders  */  //Z=36640
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", true  );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 1 );
        gp=myGuiParam::getGuiParam("Length");
        if(gp)gp->inp()->setValue( 500 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSys2DClick(Sender);  //Z=43441
        break;

    case  4:  /*  bcc  */  //Z=36688
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", true );                  /*  twinned  */  //Z=36511
        //gp->tog()->setChecked(true);
        //??EditTwRatio.Visible = true;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysCubicClick(Sender);  //Z=43446
        break;

    case  5:  /*  fcc  */  //Z=36735
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", true );                  /*  twinned  */  //Z=36511
        //gp->tog()->setChecked(true);
        //??EditTwRatio.Visible = true;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = true;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysCubicClick(Sender);  //Z=43451
        break;

    case  6:  /*  hcp  */  //Z=36782
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", true );                  /*  twinned  */  //Z=36511
        //gp->tog()->setChecked(true);
        //??EditTwRatio.Visible = true;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysHexClick(Sender);  //Z=43456
        break;

    case  7:  /*  sc  */  //Z=36829
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysCubicClick(Sender);  //Z=43461
        break;

    case  8:  /*  tetragonal centered spheres  */  //Z=36876
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", true );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysTetragonalClick(Sender);  //Z=43466
        break;

    case 17:  /*  fd3m  */  //Z=36923
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysCubicClick(Sender);  //Z=43486
        break;

    case  9:  /*  gyroid  */  //Z=37018
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 3 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysCubicClick(Sender);  //Z=43471
        break;

    case 10:  /*  OBDD  */  //Z=37066
        //ButtonPn3mVesClick(Sender);  //Z=37067
        {/*1*/  //Z=45878
            const double b = 0.5484;  //Z=45879
            const double c = 0.6200;  //Z=45880
            const double n = 2;  //Z=45881
            double a, Ri, Ro, phi; //, area, lp;  //Z=45883
            gp=myGuiParam::getGuiParam("uca");  //Z=45886
            a = gp->inp()->value();
            Ro = c*a;  //Z=45887
            Ri = b*Ro;  //Z=45888
            phi = n*4*M_PI*c*c*c*(1-b*b*b)/3.0;  //Z=45889
            //area = 8*M_PI*c*c*(1+b*b)/a;  //Z=45890
            //lp = 4*phi*(1-phi)/area;  //Z=45891
            //s.u. ComboBoxParticle.ItemIndex = 3;  //Z=45892
            gp=myGuiParam::getGuiParam("EditRadius");
            if(gp)gp->inp()->setValue( Ri );  //Z=45893
            gp=myGuiParam::getGuiParam("EditLength");
            if(gp)gp->inp()->setValue( Ro );  //Z=45894
            gp=myGuiParam::getGuiParam("EditPhi");
            if(gp)gp->inp()->setValue( phi );  //Z=45895
            //EditArea.Text = FloatToStr(area);  //Z=45896
            //Editlp.Text = FloatToStr(lp);  //Z=45897
        }/*1*/  //Z=45898
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 3 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysCubicClick(Sender);  //Z=43476
        break;

    case 11:  /*  Im3m  */  //Z=37114
        //ButtonIm3mVesClick(Sender);  //Z=37115
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 3 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysCubicClick(Sender);  //Z=43481
        break;

    case 12:  /*  none  */  //Z=37162
        myGuiParam::setEnabled( "uca", false );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, false );
        myGuiParam::setEnabled( "EditDomainSize", false );
        myGuiParam::setEnabled( "EditAzi", false );
        myGuiParam::setEnabled( "EditDebyeWaller", false );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", false );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", false );
        myGuiParam::setEnabled( "EditCeffcyl", true );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = true;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = true;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->nichts
        break;

    case 13:  /*  cpl  */  //Z=37206
        myGuiParam::setEnabled( "uca", true  );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", false );
        myGuiParam::setEnabled( "EditAzi", false );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", true );     // Editastack
        myGuiParam::setEnabled( "bcpl", true );     // Editbstack
        gp=myGuiParam::setEnabled( "ComboBoxPeak", true );
        if(gp)gp->cbs()->setCurrentIndex( 7 );   //Z=37242
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //GroupBoxAniso.Visible = true;  //Z=37238
        gp=myGuiParam::getGuiParam("SigX");
        if(gp)gp->inp()->setValue(200);  //Editdom1.Text = '200';  //Z=37239
        if(gp)gp->inp2()->setValue(200); //Editdom2.Text = '200';  //Z=37240
        if(gp)gp->inp3()->setValue(200); //Editdom3.Text = '200';  //Z=37241
        //ButtonHKLClick(Sender);
        //-->nichts
        break;

    case 14:  /*  2D-Hex-GiSAXS  */  //Z=37256
        myGuiParam::setEnabled( "uca", true );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", false );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = true;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = true;
        //??EditAbsorb.Visible = true;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->nichts
        //EditHex2h.Text = '0';  //Z=37292
        //EditHex2k.Text = '1';  //Z=37293
        //EditHex2l.Text = '0';  //Z=37294
        //RadioButtonFixU.checked = true;  //Z=37301
        //GroupBoxDataInfo.visible = true;  //Z=37306
        break;

    case 15:  /*  2D-Square-GiSAXS  */  //Z=37314
        myGuiParam::setEnabled( "uca", true );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", false );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = true;
        myGuiParam::setEnabled( "iso", false );
        //??Editacrit.visible = true;
        //??EditAbsorb.Visible = true;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->nichts
        //EditHex2h.Text = '0';  //Z=37350
        //EditHex2k.Text = '1';  //Z=37351
        //EditHex2l.Text = '0';  //Z=37352
        //RadioButtonFixU.checked = true;  //Z=37359
        //GroupBoxDataInfo.visible = true;  //Z=37364
        break;

    case 16:  /*  1D-Lam-GiSAXS  */  //Z=37372
        myGuiParam::setEnabled( "uca", true );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", false );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", true );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 2 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = true;
        myGuiParam::setEnabled( "iso", false );
        //??Editacrit.visible = true;
        //??EditAbsorb.Visible = true;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->nichts
        //EditHex2h.Text = '0';  //Z=37408
        //EditHex2k.Text = '1';  //Z=37409
        //EditHex2l.Text = '0';  //Z=37410
        //RadioButtonFixU.checked = true;  //Z=37417
        //GroupBoxDataInfo.visible = true;  //Z=37422
        break;

    case 18:  /*  orthogonal centered spheres  */  //Z=37430
        myGuiParam::setEnabled( "uca", true );
        myGuiParam::setEnabled( "ucb", true );
        myGuiParam::setEnabled( "ucc", true );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", true );
        myGuiParam::setEnabled( "EditAzi", true );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", false );
        myGuiParam::setEnabled( "bcpl", false );
        myGuiParam::setEnabled( "ComboBoxPeak", true );
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = true;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = false;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //ButtonHKLClick(Sender);
        //-->RadioButtonSysOrthoClick(Sender);  //Z=43491
        break;

    case 19:  /*  quasicrystal  */  //Z=37477
        myGuiParam::setEnabled( "uca", true );
        myGuiParam::setEnabled( "ucb", false );
        myGuiParam::setEnabled( "ucc", false );
        myGuiParam::setEnabled( ui->lblucabc, true );
        myGuiParam::setEnabled( "EditDomainSize", false );
        myGuiParam::setEnabled( "EditAzi", false );
        myGuiParam::setEnabled( "EditDebyeWaller", true );
        gp=myGuiParam::setEnabled( "CheckBoxTwinned", false );                  /*  twinned  */  //Z=36511
        gp->tog()->setChecked(false);
        //??EditTwRatio.Visible = false;
        myGuiParam::setEnabled( "HKLmax", true );
        myGuiParam::setEnabled( "acpl", true );     // Editastack
        myGuiParam::setEnabled( "bcpl", true );     // Editbstack
        gp=myGuiParam::setEnabled( "ComboBoxPeak", true );
        if(gp)gp->cbs()->setCurrentIndex( 7 );   //Z=37513
        myGuiParam::setEnabled( "EditCeffcyl", false );
        gp=myGuiParam::getGuiParam("ComboboxParticle");
        if(gp)gp->cbs()->setCurrentIndex( 0 );
        //gp=myGuiParam::getGuiParam("Length");
        //if(gp)gp->inp()->setValue( 250 );
        //??ButtonBiCont.Visible = false;
        //myGuiParam::setEnabled(ui->grpOrientation, true);
        //??EditAlphai.Visible = false;
        myGuiParam::setEnabled( "iso", true );
        //??Editacrit.visible = false;
        //??EditAbsorb.Visible = true;
        myGuiParam::setEnabled( "phi", true );
        myGuiParam::setEnabled( "theta", true );
        //??GroupBoxAniso.Visible = true;  //Z=37509
        gp=myGuiParam::getGuiParam("SigX");
        if(gp)gp->inp()->setValue(200);  //Editdom1.Text = '200';  //Z=37510
        if(gp)gp->inp2()->setValue(200); //Editdom2.Text = '200';  //Z=37511
        if(gp)gp->inp3()->setValue(200); //Editdom3.Text = '200';  //Z=37512
        //ButtonHKLClick(Sender);
        //-->nichts
        break;

    } // switch index

    on_cbsComboBoxPeak_currentIndexChanged( ui->cbsComboBoxPeak->currentIndex() );
    on_cbsComboBoxParticle_currentIndexChanged( ui->cbsComboBoxParticle->currentIndex() );

    if ( !updFlag && sender() != nullptr )
    {
        DT( qDebug() << "... cbsLType::adjustSize ..." );
        adjustSize();   // cbs LType
    }
    bIgnoreRecalc = updFlag;
    if ( !updFlag && ui->togAutoRecalc->isChecked() ) automaticRecalcDoit();
#endif
} /* on_cbsLType_currentIndexChanged() */



void SC_MainGUI::on_cbsComboBoxParticle_currentIndexChanged(int index)
{
#ifndef NOCBSCALLBACK
    DT( qDebug() << "on_cbsComboBoxParticle_currentIndexChanged()" << index << bIgnoreRecalc );
    bool updFlag = bIgnoreRecalc;
    bIgnoreRecalc = true;
    switch ( index )
    {
    case 0:   /*  sphere  */  //Z=37531
        myGuiParam::setEnabled( "EditRadius",  true,  "Rc [nm]:"       );   /*  radius  */  //Z=37532
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(R):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      false, "L [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      false, "sigma(L):"      );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 1:   /*  cylinder  */  //Z=37552
        myGuiParam::setEnabled( "EditRadius",  true,  "Rc [nm]:"       );   /*  radius  */  //Z=37553
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(R):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "L [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(L):"      );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 2:    /*  disk  */  //Z=37573
        myGuiParam::setEnabled( "EditRadius",  true,  "d/2 [nm]:"      );   /*  radius  */  //Z=37574
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(d):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "R [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(R):"      );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 3:   /*  vesicle  */  //Z=37594
        myGuiParam::setEnabled( "EditRadius",  true,  "Ri [nm]:"       );   /*  radius  */
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(R):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "Ro [nm]:"       );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(d):"      );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", false, "@0" );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 4:  /*  cube  */  //Z=37617
        myGuiParam::setEnabled( "EditRadius",  true,  "a [nm]:"        );   /*  radius  */  //Z=37618
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(a):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      false, "L [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      false, "sigma(L):"      );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", false, "@0" );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 5:   /*  ellipsoid  */  //Z=37640
        myGuiParam::setEnabled( "EditRadius",  true,  "a [nm]:"        );   /*  radius  */  //Z=37641
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(a):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "c [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(c):"      );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 6:   /*  triaxial ellipsoid  */  //Z=37661
        myGuiParam::setEnabled( "EditRadius",  true,  "a [nm]:"        );   /*  radius  */  //Z=37662
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(a):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "b [nm]:"        );   /*  length, semiaxis b  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(c):"      );   /*  sigma L  */
        myGuiParam::setEnabled( "EditRadiusi", true,  "c [nm]:"        );   /*  semiaxis c  */  //Z=37675 <<< NEU <<<<<<<<<<<<<<<<<<<<
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 7:   /*  super ellipsoid  */  //Z=37685
        myGuiParam::setEnabled( "EditRadius",  true,  "R [nm]:"        );   /*  radius  */  //Z=37686
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(R):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "L [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(L):"      );   /*  sigma L  */
        myGuiParam::setEnabled( "EditRadiusi", true,  "c [nm]:"        );   /*  semiaxis c  */  // <<< NEU <<<<<<<<<<<<<<<<<<<<
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "Alpha",       true,  "k:"             );   /*  k-parameter  */  //Z=37704 <<< NEU <<<<<<<<<<<<<<<<<<<<
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 8:   /*  superball  */  //Z=37712
        myGuiParam::setEnabled( "EditRadius",  true,  "R [nm]:"        );   /*  radius  */  //Z=37713
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(R):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "L [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(L):"      );   /*  sigma L  */
        myGuiParam::setEnabled( "EditRadiusi", true,  "c [nm]:"        );   /*  semiaxis c  */  // <<< NEU <<<<<<<<<<<<<<<<<<<<
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "Alpha",       true,  "k:"             );   /*  k-parameter  */  //Z=37704 <<< NEU <<<<<<<<<<<<<<<<<<<<
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 9:   /*  excluded volume chain  */  //Z=37739
        myGuiParam::setEnabled( "EditRadius",  true,  "Rc [nm]:"       );   /*  radius  */  //Z=37740
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(R):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "L [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "sigma(L):"      );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    case 10:   /*  Kratky Porod chain  */  //Z=37760
        myGuiParam::setEnabled( "EditRadius",  true,  "Rc [nm]:"       );   /*  radius  */  //Z=37761
        myGuiParam::setEnabled( "EditSigma",   true,  "sigma(R):"      );   /*  sigma R  */
        myGuiParam::setEnabled( "Length",      true,  "L [nm]:"        );   /*  length  */
        myGuiParam::setEnabled( "SigmaL",      true,  "Lp:"            );   /*  sigma L  */
        //myGuiParam::setEnabled( "ShellNo",     false, "no. of shells:" );   /*  no. of shells  */
        myGuiParam::setEnabled( "EditRadiusi", false, "EditRadiusi:"   );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false, "Alpha:"         );   /*  DEFAULT  */
        myGuiParam::setEnabled( "ComboBoxInterior", true );
        //myGuiParam::setEnabled( ui->grpOrientation, true );
        break;
    } // switch index
    if ( !updFlag && sender() != nullptr )
    {
        DT( qDebug() << "... cbsParticle::adjustSize ..." );
        adjustSize();   // cbs Particle
    }
    bIgnoreRecalc = updFlag;
    if ( !updFlag && ui->togAutoRecalc->isChecked() ) automaticRecalcDoit();
#endif
}


void SC_MainGUI::on_cbsOrdis_currentIndexChanged(int index)
{
#ifndef NOCBSCALLBACK
    DT( qDebug() << "on_cbsOrdis_currentIndexChanged()" << index << bIgnoreRecalc );
    bool updFlag = bIgnoreRecalc;
    bIgnoreRecalc = true;
    myGuiParam::setEnabled( "ifluc", index==6/*z-dir*/  );
    myGuiParam::setEnabled( "rfluc", index==6/*z-dir*/ );
    if ( !updFlag && sender() != nullptr )
    {
        DT( qDebug() << "... cbsOrdis::adjustSize ..." );
        adjustSize();   // cbs Ordis
    }
    bIgnoreRecalc = updFlag;
    if ( !updFlag && ui->togAutoRecalc->isChecked() ) automaticRecalcDoit();
#endif
}


void SC_MainGUI::on_cbsComboBoxInterior_currentIndexChanged(int index)
{
#ifndef NOCBSCALLBACK
    DT( qDebug() << "on_cbsComboBoxInterior_currentIndexChanged()" << index << bIgnoreRecalc );
    bool updFlag = bIgnoreRecalc;
    bIgnoreRecalc = true;
    switch ( index )
    {
    case 0:   /*  homogeneous core */  //Z=41005
        myGuiParam::setEnabled( "EditRadius",  true  );   /*  Rm     */  //Z=41006
        myGuiParam::setEnabled( "EditRadiusi", false );   /*  DEFAULT  */
        myGuiParam::setEnabled( "Alpha",       false );   /*  alpha  */  //Z=41008 - unklar, ob es der richtige Parameter ist (EditAlpha)
        myGuiParam::setEnabled( "EditRho",     false );   /*  rho    */  //Z=41010
        break;

    case 1:   /*  homogeneous core/shell */  //Z=41014
        myGuiParam::setEnabled( "EditRadius",  true  );   /*  Rm     */  //Z=41015
        myGuiParam::setEnabled( "EditRadiusi", true  );   //Z=41017  <<<NEU<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        myGuiParam::setEnabled( "Alpha",       false );   /*  alpha  */  //Z=41018 - unklar, ob es der richtige Parameter ist (EditAlpha)
        myGuiParam::setEnabled( "EditRho",     true  );   /*  rho    */  //Z=41020
        break;

    case 2:   /*  inhomogeneous core/shell */  //Z=41024
        myGuiParam::setEnabled( "EditRadius",  true  );   /*  Rm     */  //Z=41025
        myGuiParam::setEnabled( "EditRadiusi", true  );   //Z=41027  <<<NEU<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        myGuiParam::setEnabled( "Alpha",       true  );   /*  alpha  */  //Z=41028 - unklar, ob es der richtige Parameter ist (EditAlpha)
        myGuiParam::setEnabled( "EditRho",     true  );   /*  rho    */  //Z=41030
        break;

    case 3:   // if ( (ComboBoxInterior.ItemIndex == 3) || (ComboBoxInterior.ItemIndex == 4) )
    case 4:   /*  myelin  */  //Z=41034
#ifdef undef
        EditRadiusi.Visible = true;  //Z=41035
        EditRadiusi.Text = '1';  //Z=41036
        label19.visible = true;             /*  Rm  */  //Z=41037
        label19.caption = 'n:';  //Z=41038
        EditRadius.Visible = true;  //Z=41039
        EditRadius.Text = '200';  //Z=41040
        EditLength.Text = '5000';  //Z=41041
        label20.caption = 'l_ex(nm):';  //Z=41042
        label20.visible = true;             /*  alpha  */  //Z=41043
        EditAlpha.Visible = true;  //Z=41044
        EditAlpha.Text = '3.5';  //Z=41045
        label21.caption = 'l_in(nm):';  //Z=41046
        label21.visible = true;             /*  rho  */  //Z=41047
        Editrho.visible = true;  //Z=41048
        Editrho.Text = '3.5';  //Z=41049
        label2.Visible = true;  //Z=41050
        label2.Caption = 'head';  //Z=41051
        label88.Visible = true;  //Z=41052
        label88.Caption = 'tail';  //Z=41053
        label22.Visible = true;  //Z=41054
        label22.Caption = 'axon';  //Z=41055
        label50.Visible = true;  //Z=41056
        EditLattice.text = '0.001';  //Z=41057
        label50.Caption = 'intra';  //Z=41058
        label52.Visible = true;  //Z=41059
        EditLatticeb.Text = '0.001';  //Z=41060
        label52.Caption = 'extra';  //Z=41061
        EditLatticec.Text = '0.001';  //Z=41062
        EditDomainSize.visible = true;  //Z=41063
        EditDomainSize.Text = '-0.55';  //Z=41064
        EditAzi.Visible = true;  //Z=41065
        EditAzi.Text = '-0.7';  //Z=41066
        Label95.Visible = true;  //Z=41067
        Label95.Caption = 'l_h(nm):';  //Z=41068
        EditAstack.Visible = true;  //Z=41069
        EditAstack.Text = '1.0';  //Z=41070
        Label143.Visible = true;  //Z=41071
        Label143.Caption = 'l_t(nm):';  //Z=41072
        EditBstack.Visible = true;  //Z=41073
        EditBstack.Text = '2.0';  //Z=41074

        Editlattice.visible = true;  //Z=41076
        Editlatticeb.visible = true;  //Z=41077
        Editlatticec.visible = true;  //Z=41078
#endif
        break;
    }
    if ( !updFlag && sender() != nullptr )
    {
        DT( qDebug() << "... cbsInterior::adjustSize ..." );
        adjustSize();   // cbs Interior
    }
    bIgnoreRecalc = updFlag;
    if ( !updFlag && ui->togAutoRecalc->isChecked() ) automaticRecalcDoit();
#endif
}


void SC_MainGUI::on_cbsComboBoxPeak_currentIndexChanged(int index)
{
#ifndef NOCBSCALLBACK
    DT( qDebug() << "on_cbsComboBoxPeak_currentIndexChanged()" << index << bIgnoreRecalc );
    bool updFlag = bIgnoreRecalc;
    bIgnoreRecalc = true;
    switch ( index )
    {
    case 0:     /*  Gaussian  */  //Z=37787
        myGuiParam::setEnabled( "EditPeakPar", false );     /*  peak par  */  //Z=37788
        myGuiParam::setEnabled( "SigX", false );
        myGuiParam::setEnabled( ui->grpAniso, false );  //Z=37790
        break;
    case 1:    /*  Lorentzian  */  //Z=37793
        myGuiParam::setEnabled( "EditPeakPar", false );     /*  peak par  */  //Z=37794
        myGuiParam::setEnabled( "SigX", false );
        myGuiParam::setEnabled( ui->grpAniso, false );  //Z=37796
        break;
    case 2:     /*  Lorentzian I  */  //Z=37799
        myGuiParam::setEnabled( "EditPeakPar", false );     /*  peak par  */  //Z=37800
        myGuiParam::setEnabled( "SigX", false );
        myGuiParam::setEnabled( ui->grpAniso, false );  //Z=37802
        break;
    case 3:     /*  Lorentzian II  */  //Z=37805
        myGuiParam::setEnabled( "EditPeakPar", false );     /*  peak par  */  //Z=37806
        myGuiParam::setEnabled( "SigX", false );
        myGuiParam::setEnabled( ui->grpAniso, false );  //Z=37808
        break;
    case 4:     /*  Pseudo-Voigt  */  //Z=37811
        myGuiParam::setEnabled( "EditPeakPar", true  );     /*  peak par  */  //Z=37812
        myGuiParam::setEnabled( "SigX", false );
        myGuiParam::setEnabled( ui->grpAniso, false );  //Z=37814
        break;
    case 5:     /*  Pearson  */  //Z=37817
        myGuiParam::setEnabled( "EditPeakPar", true  );     /*  peak par  */  //Z=37818
        myGuiParam::setEnabled( "SigX", false );
        myGuiParam::setEnabled( ui->grpAniso, false );  //Z=37820
        break;
    case 6:    /*  gamma  */  //Z=37823
        myGuiParam::setEnabled( "EditPeakPar", true  );     /*  peak par  */  //Z=37824
        myGuiParam::setEnabled( "SigX", false );
        myGuiParam::setEnabled( ui->grpAniso, false );  //Z=37826
        break;
    case 7:     /*  anisotropic  */  //Z=37829
        myGuiParam::setEnabled( "EditPeakPar", true  );     /*  peak par  */  //Z=37830
        myGuiParam::setEnabled( "SigX", myGuiParam::isEnabled("ComboBoxPeak") );
        myGuiParam::setEnabled( ui->grpAniso, myGuiParam::isEnabled("ComboBoxPeak") );  //Z=37832
        break;
    }
    if ( !updFlag && sender() != nullptr )
    {
        DT( qDebug() << "... cbsPeak::adjustSize ..." );
        adjustSize();   // cbs Peak
    }
    bIgnoreRecalc = updFlag;
    if ( !updFlag && ui->togAutoRecalc->isChecked() ) automaticRecalcDoit();
#endif
}


/**
 * @brief SC_MainGUI::on_actionHide_unused_values_triggered
 * @param checked - Flag from the checkable menu entry
 */
void SC_MainGUI::on_actionHide_unused_values_triggered(bool checked)
{
    DT( qDebug() << "on_actionHide_unused_values_triggered()" << bIgnoreRecalc );
    bool updFlag = bIgnoreRecalc;
    bIgnoreRecalc = true;
    // Zuerst alle UI-Elemente ohne Namen sichtbar schalten, da sonst diese auf das Flag reagieren
    myGuiParam::setEnabled( ui->grpAniso, true );
    //myGuiParam::setEnabled(ui->grpOrientation, true);   wird nirgend auf false gesetzt...
    // erst jetzt das globale Flag setzen und alle Elemente mit Namen schalten
    myGuiParam::setHideFlag(checked);
    // Jetzt zur Sicherheit noch alle CBS Callbacks aufrufen
    on_cbsLType_currentIndexChanged( ui->cbsLType->currentIndex() );
    on_cbsComboBoxParticle_currentIndexChanged( ui->cbsComboBoxParticle->currentIndex() );
    on_cbsOrdis_currentIndexChanged( ui->cbsOrdis->currentIndex() );
    on_cbsComboBoxInterior_currentIndexChanged( ui->cbsComboBoxInterior->currentIndex() );
    on_cbsComboBoxPeak_currentIndexChanged( ui->cbsComboBoxPeak->currentIndex() );
    if ( !updFlag && sender() != nullptr )
    {
        DT( qDebug() << "... actionHide::adjustSize ..." );
        adjustSize();   // actionHide
    }
    bIgnoreRecalc = updFlag;
    if ( !updFlag && ui->togAutoRecalc->isChecked() ) automaticRecalcDoit();
}


/**
 * @brief SC_MainGUI::on_butCopyHdrdata2Params_clicked
 * Copies the data from the headertable into the calculation parameters.
 */
void SC_MainGUI::on_butCopyHdrdata2Params_clicked()
{
    DT( qDebug() << "on_butCopyHdrdata2Params_clicked()" );
    int c = ui->tblHeaderData->item(tblLattNbCols,0)->text().toInt();
    int r = ui->tblHeaderData->item(tblLattNbRows,0)->text().toInt();
    int gp = qMax(r,c) / 2;
    calcGui->updateParamValue( "GridPoints", gp, SETCOLMARK_IGNORED );

    double bsx = ui->tblHeaderData->item(tblLattBeamX,  0)->text().toDouble();
    double bsy = ui->tblHeaderData->item(tblLattBeamY,  0)->text().toDouble();
    if ( bsx > 0 || bsy > 0 )
    {
        bsx = bsx - gp;
        bsy = bsy - gp;
    }
    calcGui->updateParamValue( "BeamPosX",       bsx, SETCOLMARK_IGNORED );
    calcGui->updateParamValue( "BeamPosY",       bsy, SETCOLMARK_IGNORED );
    calcGui->updateParamValue( "EditWavelength", ui->tblHeaderData->item(tblLattWaveLen,0)->text().toDouble(), SETCOLMARK_IGNORED );
    calcGui->updateParamValue( "EditDet",        ui->tblHeaderData->item(tblLattSaDist, 0)->text().toDouble(), SETCOLMARK_IGNORED );
    calcGui->updateParamValue( "EditPixelX",     ui->tblHeaderData->item(tblLattPixelX, 0)->text().toDouble(), SETCOLMARK_IGNORED );  //  /1000.
    calcGui->updateParamValue( "EditPixelY",     ui->tblHeaderData->item(tblLattPixelY, 0)->text().toDouble(), SETCOLMARK_IGNORED );  //  /1000.
    // In tblHeaderData steht die Pixelsize in [mm], im Eingabefeld in [m] ...
}


/**
 * @brief SC_MainGUI::on_actionAbout_triggered
 * Shows the authors and links to the documentation.
 */
void SC_MainGUI::on_actionAbout_triggered()
{
    QMessageBox::about(this, "Scatter",
                       "Original <b>Scatter</b> in <i>Pascal</i> from <a href=\"mailto:s.foerster@fz-juelich.de\">Prof. S. Förster</a><br>"
                       "Adapted to <i>C++</i> with <i>Qt</i> and <i>Cuda</i> from <a href=\"mailto:m.wagener@fz-juelich.de\">M. Wagener</a><br>"
                       "Forschungszentrum Jülich GmbH, JCNS-1<br><br>"
                       "The idea behind this is available as a <a href=\"https://www.nature.com/articles/s41598-023-27558-8\">Nature Scientific Report</a><br><br>"
                       "Software and documentation is available on <a href=\"https://github.com/neutron-simlab/CrystalScatter\">github</a>"
                       );
}


/**
 * @brief SC_MainGUI::on_actionAbout_Qt_triggered
 * Shows the Qt version information
 */
void SC_MainGUI::on_actionAbout_Qt_triggered()
{
    QMessageBox::aboutQt(this);
}


/*
    Informationen zum Testen der Software (ein automatischer Test kommt später):

    1. Normale Berechnungen:
        Parameter aus "C:\SimLab\sas-crystal\Manuskript HyperSeries (Prof.Förster)" laden

    2. FFT / r,phi:
        Es können die gleichen Parameter wie oben verwendet werden

    3. Simplex 2D-Fit:
    3a die GUI-Version:
        einen schnell zu berechnenden Parametersatz laden und dann einfach den Fit starten

    3b die Consolenversion:
        "C:\SimLab\sas-crystal\20211206 - Fit Rezept\AutoFit-2023\autoFitCommands.txt" laden und starten

    4. Training Param Variation (testet die Consolenversion):
        Konfigurationsfiles aus "C:\SimLab\sas-crystal\20220909" nutzen, oder das Bat/Shell-Script dort

    5. AI Definitions (werden z.Zt. nicht verwendet, kommt aber sicher später wieder)

 */


void SC_MainGUI::on_butChatbotSearch_clicked()
{
    QString fn = QFileDialog::getSaveFileName( this, "Open ChatBot config file", ui->inpChatbotFile->text(), "config file (*.txt)",
                                              nullptr, QFileDialog::DontConfirmOverwrite );
    if ( fn.isEmpty() ) return;
    ui->inpChatbotFile->setText(fn);
}

void SC_MainGUI::on_butChatbotLogSearch_clicked()
{
    QString fn = QFileDialog::getSaveFileName( this, "Open ChatBot logfile", ui->inpChatbotLogfile->text(), "logfile (*.log)" );
    if ( fn.isEmpty() ) return;
    ui->inpChatbotLogfile->setText(fn);
}

void SC_MainGUI::on_butChatbotStart_clicked()
{
    ui->butChatbotStart->setEnabled(false);
    ui->butChatbotStop->setEnabled(true);

    QString fnpar = QDir::tempPath()+"/tmpForChatbot.ini";
    performSaveParamOperation( fnpar );

    if ( chatbotBackProg == nullptr )
    {
        chatbotBackProg = new QProcess;
        connect( chatbotBackProg, SIGNAL(errorOccurred(QProcess::ProcessError)),
                this, SLOT(chatbotBackProg_error(QProcess::ProcessError)) );
        connect( chatbotBackProg, SIGNAL(finished(int,QProcess::ExitStatus)),
                this, SLOT(chatbotBackProg_finished(int,QProcess::ExitStatus)) );
        connect( chatbotBackProg, SIGNAL(readyReadStandardError()),
                this, SLOT(chatbotBackProg_readyRead()) );
        connect( chatbotBackProg, SIGNAL(readyReadStandardOutput()),
                this, SLOT(chatbotBackProg_readyRead()) );
    }
    ui->lisChatbotLogOut->clear();

    QStringList params;
    params << "--chatbot" << ui->inpChatbotFile->text();
    if ( ui->togChatbotDebug->isChecked() ) params << "--cbkeep";
    params << "-p" << fnpar;
    params << "-t" << QString::number(ui->inpNumCores->value());
    if ( ui->cbsChatbotColor->currentIndex() != 3 )     // 0=grey, 1=glow, 2=earth, 3,def=temp
        params << "--color" << ui->cbsChatbotColor->currentText();
    if ( ui->cbsChatbotRotation->currentIndex() > 0 )
        params << "--rotate" << QString::number(ui->cbsChatbotRotation->currentIndex());
    if ( ui->cbsChatbotZoom->currentIndex() > 0 )
        params << "--zoom" << ui->cbsChatbotZoom->currentText();
    if ( ui->togChatbotHor->isChecked() || ui->togChatbotVert->isChecked() )
        params << "--swap" << (QString(ui->togChatbotHor->isChecked()?"H":"") + QString(ui->togChatbotVert->isChecked()?"V":""));
    if ( ui->grpChatbotLogfile->isChecked() && ! ui->inpChatbotLogfile->text().isEmpty() )
        params << "--logfile" << ui->inpChatbotLogfile->text();

    chatbotBackProgAddLog( "Params: "+params.join(" ") );
    chatbotBackProgAddLog( "Start: "+chatbotConsProg );
    chatbotBackProgAddLog( "----------" );

    chatbotBackProg->start( chatbotConsProg, params );
}

void SC_MainGUI::on_butChatbotStop_clicked()
{
    if ( chatbotBackProg != nullptr && chatbotBackProg->state() != QProcess::NotRunning )
        chatbotBackProg->kill();
    else
    {
        ui->butChatbotStart->setEnabled(true);
        ui->butChatbotStop->setEnabled(false);
    }
}

void SC_MainGUI::chatbotBackProg_error( QProcess::ProcessError err )
{
    switch ( err )
    {
    case QProcess::FailedToStart:   // The process failed to start. Either the invoked program is missing, or you may have insufficient permissions to invoke the program.
        QMessageBox::critical( this, "Chatbot input",
                              "Background process executable not found.\n"+chatbotBackProg->program(),
                              QMessageBox::Ok );
        break;
    case QProcess::Crashed:         // The process crashed some time after starting successfully.
    case QProcess::Timedout:        // The last waitFor...() function timed out. The state of QProcess is unchanged, and you can try calling waitFor...() again.
        break;
    case QProcess::WriteError:      // An error occurred when attempting to write to the process. For example, the process may not be running, or it may have closed its input channel.
    case QProcess::ReadError:       // An error occurred when attempting to read from the process. For example, the process may not be running.
        QMessageBox::critical( this, "Chatbot input",
                              "Background process unable to write to output device.",
                              QMessageBox::Ok );
        break;
    case QProcess::UnknownError:    // An unknown error occurred. This is the default return value of error().
        break;
    }
}

void SC_MainGUI::chatbotBackProg_finished( int /*code*/, QProcess::ExitStatus sta )
{
    if ( sta == QProcess::CrashExit )
        chatbotBackProgAddLog("Background process exited.");
    ui->butChatbotStart->setEnabled(true);
    ui->butChatbotStop->setEnabled(false);
}

void SC_MainGUI::chatbotBackProg_readyRead()
{
    QString str;
    str = chatbotBackProg->readAllStandardOutput();
    chatbotBackProgAddLog( str.trimmed() );
    str = chatbotBackProg->readAllStandardError();
    chatbotBackProgAddLog( str.trimmed() );
}

void SC_MainGUI::chatbotBackProgAddLog( QString msg )
{
    if ( msg.isEmpty() ) return;
    QStringList sl = msg.split(EOL);
    ui->lisChatbotLogOut->addItems(sl);
    ui->lisChatbotLogOut->scrollToBottom();
    while ( ui->lisChatbotLogOut->count() > 500 ) ui->lisChatbotLogOut->takeItem(0);
#ifndef ChatbotIgnoreImages
    if ( ui->togChatbotClpShowImg->isChecked() )
    {   // Wenn das Bild erzeugt wurde, soll es angezeigt werden.
        // Die Meldung lautet dann: "IMG:<filename> -> <zeit>ms"
        int pos = msg.indexOf("IMG:");
        qDebug() << "##" << pos << msg;
        if ( pos < 0 ) return;
        QString fn = msg.mid(pos+4);
        pos = fn.indexOf(".png");
        fn.truncate(pos+4);
        qDebug() << "CBFILE:" << fn << "  Reading not yet implemented";

    }
#endif
}

void SC_MainGUI::on_butChatbotTrainfileSearch_clicked()
{
    QString fn = QFileDialog::getSaveFileName( this, "Save ChatBot trainfile", ui->inpChatbotTrainfile->text(),
                                              "Train Config (*.txt)", nullptr, QFileDialog::DontConfirmOverwrite );
    if ( fn.isEmpty() ) return;
    ui->inpChatbotTrainfile->setText(fn);
}

void SC_MainGUI::chatbotSaveConfigHelper( QFile &fTrain, QString key, bool showena, bool isena )
{
    if ( key.contains("BeamPos") ) qDebug() << key;
    QString dbg="";
    QString prtkey;
    if ( key.startsWith("VAx")  ) prtkey = key.mid(1); // Das 'V' stört bei der Ausgabe
    else if ( key.startsWith("Edit") ) prtkey = key.mid(4); // Das 'Edit' stört auch
    else if ( key.startsWith("CheckBox") ) prtkey = key.mid(8);
    else if ( key.startsWith("ComboBox") ) prtkey = "CB"+key.mid(8);
    else if ( key.startsWith("RadBut") ) prtkey = key.mid(6);
    else if ( key.startsWith("RadioButton") ) prtkey = "RB"+key.mid(11);
    else prtkey = key;
    QString val = calcGui->currentParamValueStr(key, true/*text*/ );
    if ( showena || isena )
    {
        paramHelper *par = calcGui->params.value(key,nullptr);
        if ( par != nullptr && ! par->tooltip.isEmpty() && prtkey!=par->tooltip )
            fTrain.write(qPrintable("# "+prtkey+"  ("+par->tooltip+")"+EOL));
        else
            fTrain.write(qPrintable("# Enable/Disable "+prtkey+EOL));
        if ( isena )
            fTrain.write(qPrintable("ena_"+prtkey+"=true"+EOL));
        else
        {
            fTrain.write(qPrintable("ena_"+prtkey+"="+(calcGui->isCurrentParameterVisible(key,dbg)?"true":"false")+EOL));
            if ( dbg == "Norm" ) dbg=""; else dbg="   # "+dbg;
        }
    }
    if ( val.contains("{") )
    {   // Bei den ComboBoxen sollten nur Zahlen als Wert und der Text als Kommentar geschrieben werden
        int p1 = val.indexOf("{");
        int p2 = val.indexOf("}");
        dbg = "   # "+val.left(p1-1).trimmed();
        val = val.mid(p1+1,p2-p1-1);
    }
    fTrain.write(qPrintable("val_"+prtkey+"="+val+dbg+EOL+EOL));
}

void SC_MainGUI::on_butChatbotSaveConfig_clicked()
{
    QFile fTrain(ui->inpChatbotTrainfile->text());
    bool isappended = true;
    //if ( ! fTrain.open(QIODevice::Append) ) // ist zum Testen einfacher, wenn die Datei immer überschrieben wird
    {
        isappended = false;
        if ( ! fTrain.open(QIODevice::WriteOnly) )
        {
            QMessageBox::critical(this,"Saving training configuration",fTrain.errorString(),QMessageBox::Ok);
            return;
        }
    }

    fTrain.write("# Please add the name of the experiment" EOL);
    fTrain.write(qPrintable("Name of the experiment = "+ui->txtComment->text()+EOL+EOL));

    QStringList allkeys = calcGui->paramsForMethod(false/*num*/, true/*glob*/, false/*fit*/ );
    while ( allkeys.size() > 0 )
    {
        QString key = allkeys.first();
        if ( key.startsWith("EditAxis") || key.startsWith("Editdom") ||
             key.startsWith("ExpandImage") || key.startsWith("RadioButtonQ") ||
             key.contains("BeamPos") || key.contains("CalcQmax") ||
             key.contains("GridPoints") || key.contains("CenterBeam") ||
             key.contains("CenterMidpoint") || key.startsWith("Qmax") ||
             key.contains("EditRelDis") || key.contains("EditDist")
            )
        {
            allkeys.removeAll(key);
            continue;
        }
        if ( calcGui->isCurrentParameterValid(key, false/*forfit*/) )
        {
            chatbotSaveConfigHelper( fTrain, key, true, false );
        }
        allkeys.removeAll(key); // damit auch doppelte rausgehen
    }

    chatbotSaveConfigHelper( fTrain, "GridPoints", true, true );

    fTrain.write("# Enable/Disable BeamPos" EOL);
    fTrain.write("ena_BeamPos=true" EOL);
    if ( calcGui->currentParamValueInt("CenterBeam") )
    {   // angegebene Koordinaten
        fTrain.write(qPrintable("val_BeamPosX="+calcGui->currentParamValueStr("BeamPosX",true/*text*/)+EOL));
        fTrain.write(qPrintable("val_BeamPosY="+calcGui->currentParamValueStr("BeamPosY",true/*text*/)+EOL+EOL));
    }
    else
    {   // Mittelpunkt ist 0 wenn es von -GridPoints bis +GridPoints geht
        fTrain.write("val_BeamPosX=0  # -GridPoints .. +GridPoints" EOL);
        fTrain.write("val_BeamPosY=0  # -GridPoints .. +GridPoints" EOL EOL);
    }

    //fTrain.write("# Enable/Disable Generate PNG" EOL);
    //fTrain.write("ena_generatePNG=true" EOL EOL);

    //fTrain.write("# Number of Images" EOL);     // Beispiel wohl von TPV
    //fTrain.write("val_numimg=1" EOL EOL);

    fTrain.write("# Output Path" EOL);
    fTrain.write("val_outPath=." EOL EOL);  // Outputfile kommt dann bei der Berechnung... (TODO)

    fTrain.write(EOL EOL EOL EOL);  // Trennung zwischen den Blöcken
    fTrain.close();
    ui->statusbar->showMessage("File "+fTrain.fileName()+(isappended?" appended.":" generated."),5000);
}


void SC_MainGUI::on_butChatbotReadClipboard_clicked()
{
    QFile f(ui->inpChatbotFile->text());
    if ( ! f.open(QIODevice::WriteOnly) )
    {
        ui->statusbar->showMessage( f.fileName()+" "+f.errorString(), 5000 );
        ui->butChatbotReadClipboard->setEnabled(false);
        return;
    }
    f.write( qPrintable(qApp->clipboard()->text()) );
    f.write(EOL);
    f.close();
    qApp->clipboard()->clear();
    if ( chatbotBackProg == nullptr || ! chatbotBackProg->isOpen() )
    {   // Jetzt läuft der Hintergrundprozess NICHT ...
        on_butChatbotStart_clicked(); // Macht sofort die erste Berechnung
    }
    // Ansonsten wird der Prozess die Datei finden.
}

void SC_MainGUI::chatbotClipboardChanged(QClipboard::Mode)
{
    if ( chatbotConsProg.isEmpty() ) return;
    ui->butChatbotReadClipboard->setEnabled( ! qApp->clipboard()->text().isEmpty() );
}

void SC_MainGUI::on_togUseAdaptiveStep_toggled(bool checked)
{
    myGuiParam::updateAdaptiveSteps(checked);
}
