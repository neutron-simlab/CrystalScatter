#include "sc_calcCons.h"
#include <QCoreApplication>
#include <QCommandLineParser>
#include <QFile>
#include <QDir>
#include <QDateTime>
#include <QStringList>
#include <thread>
#include <iostream>
#include "sc_postproc.h"
#include "widimage.h"
#include "sc_readdata.h"
#include "sc_simplexfit2d.h"

#include <QDebug> // nur zum Test

#define USEGENCOLORTBL


//#define USE_INTERNAL_DEFAULTS
// Versuch, aber noch nicht verwendet...
#define INTDEF_PARAMFILE "sas_scatter2/Test B-0-6 (BCT).ini"
#define INTDEF_IMGFILE   "20211206 - Fit Rezept/Daten/46424-AWM3-OX-SF01-CDCl3-B-2-95.dat"
#define INTDEF_THREADS   "4"
#define INTDEF_LOGFILE   "20211206 - Fit Rezept/AutoFit-Multifile/AutoFit-Test" // -Win.log oder -Gpu.log
#define INTDEF_CSVFILE   "20211206 - Fit Rezept/AutoFit-Test" // -Win.csv oder -Gpu.csv
#define INTDEF_AUTOFIT   "20211206 - Fit Rezept/AutoFit-Multifile/autoFitCommands.txt"
#define INTDEF_AFOUTPUT  "20211206 - Fit Rezept/AutoFit-Multifile"


#ifdef Q_OS_LINUX
#define EOL "\n"
#else
#define EOL "\r\n"
#endif



#ifdef WIN32
//#include "debughandler.h"
#else
#include <QLoggingCategory>
#endif


//QStringList slCalcArt;
int iCurCalcArt;

QString category;

QFile *flog=nullptr, *fcsv=nullptr;

QString imgColorTbl;
bool    imgSwapH, imgSwapV;
int     imgRot, imgZoom;
bool    noConsolOutput;

bool    doPrintParameter;

QStringList strColTblNames = {"grey", // 0="Greyscale"
                              "glow", // 1="Glowing colors"
                              "earth",  // 2="Earth colors"
                              "temp", // 3="Temperature",
                              "scol",  /* 4*/"Scatter-Color", /* 5*/"Scatter-Colorwheel", /* 6*/"Scatter-Temperature",  /* 7*/"Scatter-Geo",
// /* 8*/"Scatter-Desy",  /* 9*/"Scatter-Grey",       /*10*/"Scatter-Inverse grey", /*11*/"Scatter-Jet",
// /*12*/"Residuen Glow", /*13*/"Residuen Temp" };
};



static widImage *addImage( int x0, int x1, int y0, int y1, double *d, QString /*title*/ )
{
    widImage *img = new widImage();
    img->setConfig(imgColorTbl, imgSwapH, imgSwapV, imgRot, imgZoom, false);
    img->setData( x0, x1, y0, y1, d );
    return img;
}





int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);

    QCoreApplication::setApplicationName("sas-scatter2");
    QCoreApplication::setApplicationVersion("Version " MYVERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription("Calculation of simulated SAXS / SANS using a theoretical model of the structure.");
    parser.addHelpOption();
    parser.addVersionOption();

    QCommandLineOption paramFileOption( QStringList() << "p" << "paramfile",
                                        "Filepath loaded as the parameter set for all calculations (*.ini).",
                                        "paramfile" );
    parser.addOption(paramFileOption);

    QCommandLineOption aiFileOption( QStringList() << "c" << "seqfile",
                                     "Filepath loaded as the AI definition sequence (*.sas_airun).",
                                     "aifile" );
    parser.addOption(aiFileOption);

    QCommandLineOption imgOutputOption( QStringList() << "i" << "img",
                                        "Filepath and name to save the image (*.png).\nOr the image to load for automatic fit (*.dat)",
                                        "imgfile" );
    parser.addOption(imgOutputOption);

    QCommandLineOption threadsOption( QStringList() << "t" << "threads",
                                      "Number of threads used (default=max cores), set to 0 to use GPU if available.",
                                      "threads" );
    parser.addOption(threadsOption);

    QCommandLineOption loggingOption( QStringList() << "l" << "logfile",
                                      "Filepath and name of the logfile for informations of the process (*.log).",
                                      "logfile" );
    parser.addOption(loggingOption);

    QCommandLineOption csvOption( "csv",
                                  "Filepath and name of the csvfile for some statistic outputs during 2D-Fit (*.csv).",
                                  "csvfile" );
    parser.addOption(csvOption);

    QCommandLineOption autofitOption( QStringList() << "f" << "autofit",
                                      "Filepath and name for the automatic fit routine (*.txt)."
                                      " The parameterfile and the dataset are taken from this fit file."
                                      " Or use the -p for the parameterfile and the -i for the dataset."
                                      " Use -l for the global logfile.",
                                      "autofit" );
    parser.addOption(autofitOption);

    QCommandLineOption outputOption( QStringList() << "o" << "output",
                                     "Filepath and name of the single output image of the automatic fit routine or the path to save the multifile outputs.",
                                     "output" );
    parser.addOption(outputOption);

    QCommandLineOption nofftOption( "nofft",
                                    "If set, the FFT calculation during AI generation is skipped." );
    parser.addOption(nofftOption);

    QCommandLineOption noConsolOption( "noconsole",
                                       "If set, no output is printed to the console during auto fit (except statistic info)." );
    parser.addOption(noConsolOption);

    QCommandLineOption imgSwapOption( "swap",
                                      "Set the image swapping and can contain the characters 'H' and/or 'V'.",
                                      "swap", "" );
    parser.addOption(imgSwapOption);

    QCommandLineOption imgRotateOption( "rotate",
                                        "Set the image ratoation to 0(=0°), 1=(90°), 2=(180°), 3=(270°).",
                                        "rotate", "0" );
    parser.addOption(imgRotateOption);

    QCommandLineOption imgZoomOption( "zoom",
                                      "Set the image zoom factor to 1, 2, 4.",
                                      "zoom", "1" );
    parser.addOption(imgZoomOption);

    QCommandLineOption imgColorOption( "color",
                                       "Set the used colortable to one of: grey, glow, earth, temp.",
                                       "color" /*, "grey"*/ );
    parser.addOption(imgColorOption);

    QCommandLineOption tpvOption( "tpv",
                                  "Inputfile (*.sas_tpv) for the train parameter variation section to generate a number of random modified datasets (incl. r,phi and fft images)."
                                  " This parameter can be given more than once. If this is the last parameter, a wildcard for the filename can be used too.",
                                  "tpv" );
    parser.addOption(tpvOption);

    QCommandLineOption tpvFactorOption( "tpvimg",
                                       "Number of images to generate, will overwrite the factor value in the *.sas_tpv files.",
                                       "tpvimg" );
    parser.addOption(tpvFactorOption);

#ifdef USEGENCOLORTBL
    QCommandLineOption genColTblOption( "gencolortable",
                                        "Special option to convert an image file to a usable color table.",
                                        "file" );
    parser.addOption(genColTblOption);
#endif // USEGENCOLORTBL

    QCommandLineOption timetestOption( "timetest",
                                       "Start timing test. No output image is generated. Use -o to set the"
                                       " generated output file with the statistics and -p for the parametr"
                                       " file for the calculation. The given config file can be generated"
                                       " with the GUI and contains the test definitions.",
                                       "timetest" );
    parser.addOption(timetestOption);

    QCommandLineOption chatbotOption( "chatbot",
                                      "The program waits for the given input file (generated from a chatbot),"
                                      " calculates and saves the image and removes the input file.",
                                      "file" );
    parser.addOption(chatbotOption);

    QCommandLineOption cbkeepOption( "cbkeep",
                                     "If set, the inputfile is kept and the program ends after one calculation."
                                     " This is useful during debugging, not needed in normal operations." );
    parser.addOption(cbkeepOption);

    QCommandLineOption printParOption( "prtpar",
                                    "If set, the parameter are printed before calculation (internal)." );
    printParOption.setFlags(QCommandLineOption::HiddenFromHelp);
    parser.addOption(printParOption);

    // Process the actual command line arguments given by the user
    parser.process(app);

#ifdef USEGENCOLORTBL
    QString tmp = parser.value(genColTblOption);
    if ( ! tmp.isEmpty() )
    {
        std::cerr << "GenColTbl: " << qPrintable(tmp) << std::endl;
        QImage img(tmp);
        if ( img.isNull() )
        {
            std::cerr << "Image load error" << std::endl;
            return 1;
        }
        std::cerr << "Image size " << img.width() << "x" << img.height() << std::endl;
        if ( img.width() < img.height() )
        {
            QTransform mat = QTransform(); // Identity
            img = img.transformed( mat.rotate( 90 ) );
        }
        if ( img.width() != 256 )
        {
            img = img.scaledToWidth( 256 );
        }
        int y = img.height() / 2;
        //static int rEarth[10] = { 0,  85,   0,  32,   0,   0,   0, 124, 255, 255 };
        //static int gEarth[10] = { 0,  26,   0, 178, 139, 100, 139, 252, 255, 255 };
        //static int bEarth[10] = { 0, 139, 128, 170,  69,   0,   0,   0,   0, 255 };
        QStringList r,g,b;
        for ( int x=0; x<img.width(); x++ )
        {
            r << QString("%1").arg(QColor(img.pixel(x,y)).red(),4);
            g << QString("%1").arg(QColor(img.pixel(x,y)).green(),4);
            b << QString("%1").arg(QColor(img.pixel(x,y)).blue(),4);
        }
        std::cerr << "static int r[256] = {" << qPrintable(r.join(",")) << "};" << std::endl;
        std::cerr << "static int g[256] = {" << qPrintable(g.join(",")) << "};" << std::endl;
        std::cerr << "static int b[256] = {" << qPrintable(b.join(",")) << "};" << std::endl;
        return  0;
    }
#endif // USEGENCOLORTBL

    QString paramfile = parser.value(paramFileOption);
    QString aifile    = parser.value(aiFileOption);
    QString imgfile   = parser.value(imgOutputOption);
    QString threadpar = parser.value(threadsOption);
    QString logfile   = parser.value(loggingOption);
    QString csvfile   = parser.value(csvOption);
    QString autofit   = parser.value(autofitOption);
    QString af_output = parser.value(outputOption);
    bool nofft        = parser.isSet(nofftOption);
    noConsolOutput    = parser.isSet(noConsolOption);
    QString timetest  = parser.value(timetestOption);
    QString chatbot   = parser.value(chatbotOption);
    bool cbkeepfile   = parser.isSet(cbkeepOption);
    doPrintParameter  = parser.isSet(printParOption);

    // Die TPV-Files werden hinter der Option "--tpv" geschrieben. Wenn der User
    // jetzt Files angibt (z.B. mit *) dann sollen die Positions-Argumente nicht
    // ignoriert werden.
    QStringList checkWildcard(QString);
    QStringList tpvfiles = parser.values(tpvOption);
    int tpvImagefactor = 0;
    if ( tpvfiles.size() > 0 )
    {
        std::cerr << "----------" << std::endl;
        foreach ( QString s, tpvfiles )
        {
            tpvfiles.removeOne(s);
            QStringList tmp = checkWildcard(s);
            foreach ( QString s2, tmp )
                if ( s2.endsWith(".sas_tpv") )
                    tpvfiles << s2;
        }
        QStringList tpvAddon = parser.positionalArguments();
        foreach ( QString s, tpvAddon )
        {
            QStringList tmp = checkWildcard(s);
            foreach ( QString s2, tmp )
                if ( s2.endsWith(".sas_tpv") )
                    tpvfiles << s2;
        }
        if ( tpvfiles.size() > 0 )
        {
            QString tf = parser.value(tpvFactorOption);
            bool ok;
            int tmp = tf.toInt(&ok);
            if ( ok && tmp > 0 ) tpvImagefactor = tmp;
        }
        std::cerr << "----------" << std::endl;
    }

    bool ok;
    imgSwapH = parser.value(imgSwapOption).contains("H",Qt::CaseInsensitive);
    imgSwapV = parser.value(imgSwapOption).contains("V",Qt::CaseInsensitive);
    imgRot   = parser.value(imgRotateOption).toInt(&ok);
    if ( !ok || (imgRot < 0) || (imgRot > 3) ) imgRot = 0;
    imgZoom  = parser.value(imgZoomOption).toInt(&ok);
    if ( !ok || (imgZoom!=1 && imgZoom!=2 && imgZoom!=4) ) imgZoom = 1;
    if ( parser.value(imgColorOption).startsWith("grey",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames().at(0);
    else if ( parser.value(imgColorOption).startsWith("glow",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames().at(1);
    else if ( parser.value(imgColorOption).startsWith("earth",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames().at(2);
    else if ( parser.value(imgColorOption).startsWith("temp",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames().at(6);
    else
        imgColorTbl = widImage::slColorNames().at(0); // tblGrey;

    // /* 0*/"Greyscale",     /* 1*/"Glowing colors",     /* 2*/"Earth colors",         /* 3*/"Temperature",
    // /* 4*/"Scatter-Color", /* 5*/"Scatter-Colorwheel", /* 6*/"Scatter-Temperature",  /* 7*/"Scatter-Geo",
    // /* 8*/"Scatter-Desy",  /* 9*/"Scatter-Grey",       /*10*/"Scatter-Inverse grey", /*11*/"Scatter-Jet",
    // /*12*/"Residuen Glow", /*13*/"Residuen Temp" };


#ifdef USE_INTERNAL_DEFAULTS
#ifdef WIN32
    // Autofit-Test
//    if ( paramfile.isEmpty() )
//        paramfile = "C:/SimLab/sas-crystal/sas_scatter2/Test B-0-6 (BCT).ini";
//    if ( imgfile.isEmpty() )
//        imgfile = "C:/SimLab/sas-crystal/20211206 - Fit Rezept/Daten/46424-AWM3-OX-SF01-CDCl3-B-2-95.dat";
    if ( threadpar.isEmpty() )
        threadpar = "4";
    if ( logfile.isEmpty() )
        logfile = "C:/SimLab/sas-crystal/20211206 - Fit Rezept/AutoFit-Multifile/AutoFit-Test-Win.log";
//    if ( csvfile.isEmpty() )
//        csvfile = "C:/SimLab/sas-crystal/20211206 - Fit Rezept/AutoFit-Test-Win.csv";
    if ( autofit.isEmpty() )
        autofit = "C:/SimLab/sas-crystal/20211206 - Fit Rezept/AutoFit-Multifile/autoFitCommands.txt";
    if ( af_output.isEmpty() )
        af_output = "C:/SimLab/sas-crystal/20211206 - Fit Rezept/AutoFit-Multifile"; // AutoFit-Test-Win.png";
    nofft = true;
#else
    // Autofit-Test
    //if ( paramfile.isEmpty() )
    //    paramfile = "/home/wagener/sas-crystal/sas_scatter2/Test B-0-6 (BCT).ini";
    //if ( imgfile.isEmpty() )
    //    imgfile = "/home/wagener/sas-crystal/20211206 - Fit Rezept/Daten/46424-AWM3-OX-SF01-CDCl3-B-2-95.dat";
    if ( threadpar.isEmpty() )
        threadpar = "6";
    if ( logfile.isEmpty() )
        logfile = "/home/wagener/sas-crystal/20211206 - Fit Rezept/AutoFit-Multifile/AutoFit-Test-CPU.log";
    //if ( csvfile.isEmpty() )
    //    csvfile = "/home/wagener/sas-crystal/20211206 - Fit Rezept/AutoFit-Multifile/AutoFit-Test-CPU.csv";
    if ( autofit.isEmpty() )
        autofit = "/home/wagener/sas-crystal/20211206 - Fit Rezept/AutoFit-Multifile/AutoFit/autoFitCommands.txt";
    if ( af_output.isEmpty() )
        af_output = "/home/wagener/sas-crystal/20211206 - Fit Rezept/AutoFit-Multifile"; // AutoFit-Test-CPU.png";
    nofft = true;
#endif
#endif // USE_INTERNAL_DEFAULTS

    if ( parser.value(imgColorOption).isEmpty() )
        imgColorTbl = widImage::slColorNames()[6]; // tblGRtemp;

#ifdef WIN32
    // To avoid empty lines printed every second after one qDebug() output
    //qInstallMessageHandler(debugWinMsgHandler);
#else
    // To enable qDebug() under Linux (https://forum.qt.io/topic/81430/unable-to-see-the-qdebug-messages-on-console/15)
    QLoggingCategory::defaultCategory()->setEnabled(QtDebugMsg, true);
#endif

    if ( !timetest.isEmpty() )
    {
        // --timetest "C:\SimLab\sas-crystal\20240115-Neuer_Switch_Ablauf\timetest.ini"
        //   -p "C:\SimLab\sas-crystal\20240115-Neuer_Switch_Ablauf\partSphere_peakGauss_lattFCC_ordisGauss.ini"
        //   -o "C:\SimLab\sas-crystal\20240115-Neuer_Switch_Ablauf\partSphere_peakGauss_lattFCC_ordisGauss_CONStimetest.txt"
        std::cerr << "  TIMETEST: " << qPrintable(timetest) << std::endl;
        std::cerr << "    output: " << qPrintable(af_output) << std::endl;
        std::cerr << " parameter: " << qPrintable(paramfile) << std::endl;
        void performTimeTest( QString par, QString cfg, QString out );
        if ( af_output.isEmpty() || paramfile.isEmpty() )
            std::cerr << "No outputfile and no parameterfile given. Aborted." << std::endl;
        else
            performTimeTest( paramfile, timetest, af_output );
        return 0;
    }

    if ( !autofit.isEmpty() )
    {   // Bei automatischem Fit wird die Eingabe von AI ignoriert
        aifile = "";
        tpvfiles.clear();
    }

    //------------------------------------------------------------------

    bool onlyCalc = aifile.isEmpty();

    int nthreads = static_cast<int>(std::thread::hardware_concurrency());
    if ( !threadpar.isEmpty() )
    {
        bool ok;
        int n = threadpar.toInt(&ok);
        if ( n >= 0 && n <= nthreads && ok ) nthreads = n;
    }

    if ( !csvfile.isEmpty() )
    {
        fcsv = new QFile(csvfile);
        if ( !fcsv->open(QIODevice::WriteOnly) )
        {
            std::cerr << qPrintable(csvfile + " " + fcsv->errorString()) << std::endl;
            return 1;
        }
    }

    if ( !logfile.isEmpty() )
    {
        flog = new QFile(logfile);
        if ( !flog->open(QIODevice::WriteOnly) )
        {
            std::cerr << qPrintable(logfile + " " + flog->errorString()) << std::endl;
            return 1;
        }
        if ( onlyCalc && autofit.isEmpty() && chatbot.isEmpty() )
            flog->write(qPrintable(QString("Only one image calculated.")+EOL));

        flog->write( qPrintable(        "AI-Calculationfile: "+aifile+EOL) );
        flog->write( qPrintable(        "     Parameterfile: "+paramfile+EOL) );
        flog->write( qPrintable(        "         Imagefile: "+imgfile+EOL) );
        flog->write( qPrintable(QString("      Threads used: %1").arg(nthreads)+EOL) );
        flog->write( qPrintable(        "     Automatic fit: "+autofit+EOL) );
        flog->write( qPrintable(        "   Auto fit output: "+af_output+EOL) );
        flog->write( qPrintable(        "    CSV-Outputfile: "+csvfile+EOL) );
        flog->write( qPrintable(        "    TPV-Inputfiles: "+QString::number(tpvfiles.size())+EOL) );
        flog->write( qPrintable(QString("       No FFT Flag: %1").arg(nofft)+EOL) );
        flog->write( qPrintable(QString("             Image:%1%2 Rot=%3deg, Zoom=%4, Color=%5")
                                   .arg((imgSwapH?" SwapHor,":""),(imgSwapV?" SwapVert,":""))
                                   .arg(imgRot*90).arg(imgZoom).arg(imgColorTbl)+EOL) );
        flog->write( qPrintable(        "      Chatbot file: "+chatbot+EOL) );
        flog->write( qPrintable(QString(" Chatbot debug mode: %1").arg(cbkeepfile)+EOL) );
    }
    std::cerr << "   CALC: " << qPrintable(aifile)    << std::endl;
    std::cerr << " PARAMS: " << qPrintable(paramfile) << std::endl;
    std::cerr << "  IMAGE: " << qPrintable(imgfile)   << std::endl;
    std::cerr << "THREADS: " << nthreads              << std::endl;
    std::cerr << "AUTOFIT: " << qPrintable(autofit)   << std::endl;
    std::cerr << " OUTPUT: " << qPrintable(af_output) << std::endl;
    std::cerr << "LOGFILE: " << qPrintable(logfile)   << std::endl;
    std::cerr << "    CSV: " << qPrintable(csvfile)   << std::endl;
    std::cerr << "    TPV: " << qPrintable(tpvfiles.join(", ")) << std::endl;
    std::cerr << " no FFT: " << nofft                 << std::endl;
    std::cerr << qPrintable(QString("  Image:%1%2 Rot=%3deg, Zoom=%4, Color=%5")
                                .arg((imgSwapH?" SwapHor,":""),(imgSwapV?" SwapVert,":""))
                                .arg(imgRot*90).arg(imgZoom).arg(imgColorTbl)) << std::endl;
    std::cerr << "CHATBOT: " << qPrintable(chatbot)    << std::endl;
    std::cerr << "CBDebug: " << cbkeepfile << std::endl;
    if ( doPrintParameter )
        std::cerr << " PrtPar: yes (internal)" << std::endl;

    QDateTime dtStart = QDateTime::currentDateTime();

    SC_CalcCons *calc = new SC_CalcCons;

    //slCalcArt = calc->getCalcTypes();
    if ( !paramfile.isEmpty() ) calc->loadParameter( paramfile );

    if ( ! chatbot.isEmpty() )
    {
        void waitForChatbot( SC_CalcCons *calc, QString cbfile, /*QString param,*/ int nthreads, bool keep );
        if ( flog )
        {
            flog->write("*****************************************" EOL);
            flog->write("* Chatbot mode                          *" EOL);
            flog->write("*****************************************" EOL);
            if ( fcsv ) flog->write("--> CSV file will be ignored." EOL);
            flog->close(); // Wird bei jedem Eintrag unten wieder geöffnet
        }
        if ( fcsv ) fcsv->close();
        waitForChatbot( calc, chatbot, /*paramfile,*/ nthreads, cbkeepfile );
        // Endlos-Schleife. Wird nur durch ein Terminate unterbrochen!
        return 0;
    }
    else if ( tpvfiles.size() > 0 )
    {   // Generate AI files with TPV from one file ...
        void generateAIfiles( SC_CalcCons *calc, QString aifile, bool nofft, int nthreads );
        if ( flog )
        {
            flog->write("*****************************************" EOL);
            flog->write("* Generate AI files from TPV filelist   *" EOL);
            flog->write("*****************************************" EOL);
        }
        foreach ( QString f, tpvfiles )
        {
            // In diesem File stehen auch die Rechen-Parameter und können direkt verarbeitet werden.
            calc->loadParameter( f );
            // Jetzt muss aus den vorhandenen TPV-Flags und FFT-Einstellungen das AI File als
            //  passende Eingabe für die generateAIfiles() Routine umgesetzt werden.
            // Im Parameter 'aifile' wird der Fileinhalt als String mit dem Trennzeichen '@'
            //  und beginnend mit '@' übergeben:
            //    PATH|C:/SimLab/sas-crystal/20220909/Testfiles/SphLargMono
            //    COLTBL|*DAT*
            //    TPV|xx|BS0.5;0.5;0.5;0.5|RP
            //    FFT|InLog;128;;128;OutAbs;;OutSwap
            //    0|#Calculation#=0|#FACTOR#=5|EditDomainSize=0.1|EditSigma=0.1|I0=0.1|
            QString aifile;
            QSettings sets(f,QSettings::IniFormat);
            if ( sets.value("val_outPath","").toString().isEmpty() )
            {
                if ( flog )
                    flog->write(qPrintable(EOL+QString("ERROR: file not found: ")+f+EOL));
                std::cerr << std::endl << "ERROR: file not found: " << qPrintable(f) << std::endl << std::flush;
                return 0;
            }
            aifile  = "@PATH|" + sets.value("val_outPath").toString();
            aifile += "@COLTBL|*DAT*";
            aifile += "@TPV|xx"; // Damit immer 2 Einträge da sind ...
            if ( sets.value("ena_addBS",false).toBool() )
                aifile += QString("|BS%1;%2;%3;%4").arg(sets.value("val_beamx0",0).toDouble())
                           .arg(sets.value("val_beamy0",0).toDouble())
                           .arg(sets.value("val_nbeam",0).toDouble())
                           .arg(sets.value("val_mbeam",0).toDouble());
            if ( sets.value("ena_addLines",false).toBool() )
                aifile += QString("|LI%1;%2;%3;%4").arg(sets.value("val_addLinesH",1).toInt())
                              .arg(sets.value("val_addLinesV",1).toInt())
                              .arg(sets.value("val_addLinesHwidth",2).toInt())
                              .arg(sets.value("val_addLinesVwidth",2).toInt());
            if ( sets.value("ena_addNoise",false).toBool() )
                aifile += "|NO";
            if ( sets.value("ena_convolute",false).toBool() )
                aifile += "|CO";
            bool haveRphi = sets.value("ena_calcRphi",false).toBool();
            if ( haveRphi )
                aifile += "|RP";
            if ( sets.value("ena_saveExtra",false).toBool() )
                aifile += "|IX";
            if ( sets.value("ena_generatePNG",false).toBool() )
                aifile += "|GP";
            if ( sets.value("ena_scaleScat",false).toBool() )
                aifile += "|SC";
            bool haveFFT = sets.value("ena_calcFFT",false).toBool();
            if ( haveFFT || haveRphi )
            {
                sets.beginGroup("FFT");
                static const QString fftOutFormat[4] = { "OutRe", "OutIm", "OutAbs", "OutSpec" };
                static const QString fftSize[5] = { "32", "64", "128", "256", "512" };
                QString rphiScale = "";
                if ( haveRphi )
                {
                    if ( sets.value("FFTScaleRphi",false).toBool()  ) rphiScale += "Scale ";
                    if ( sets.value("FFTclipRphi",false).toBool()   ) rphiScale += "Clip1 ";
                    if ( sets.value("FFTclip40Rphi",false).toBool() ) rphiScale += "Clip4 ";
                }
                QString fftScale = "";
                if ( haveFFT )
                {
                    if ( sets.value("FFTScaleOutput",false).toBool()  ) fftScale += "Scale ";
                    if ( sets.value("FFTclipOutput",false).toBool()   ) fftScale += "Clip1 ";
                    if ( sets.value("FFTclip40Output",false).toBool() ) fftScale += "Clip4 ";
                }
                aifile += QString("@FFT|%1;%2;%3;%4;%5;%6;%7;%8;%9")
                              .arg(sets.value("FFTLinInput",false).toBool()?"InLin":"InLog")
                              .arg(haveRphi?fftSize[sets.value("FFTsizeRphi",2).toInt()]:0)
                              .arg(rphiScale.trimmed())
                              .arg(haveFFT?fftSize[sets.value("FFTsizeOut",2).toInt()]:0)
                              .arg(fftOutFormat[sets.value("FFToutput",0).toInt()])
                              .arg(fftScale.trimmed())
                              .arg(sets.value("FFTSwapOutput",false).toBool()?"OutSwap":"OutNoSwap")
                              .arg(sets.value("FFTcutOutEna",false).toBool()?sets.value("FFTcutOutX",64).toInt():0)
                              .arg(sets.value("FFTcutOutEna",false).toBool()?sets.value("FFTcutOutY",64).toInt():0);
                sets.endGroup();
            }

            aifile += "@0|#Calculation#=0|#FACTOR#=";
            if ( tpvImagefactor > 0 )
                aifile += QString::number(tpvImagefactor);
            else
                aifile += QString::number(sets.value("val_numimg",1).toInt());
            if ( sets.value("ena_io",false).toBool() )
                aifile += "|I0=" + sets.value("val_io","0.5").toString();
            if ( sets.value("ena_base",false).toBool() )
                aifile += "|Base=" + sets.value("val_base","0.5").toString();
            if ( sets.value("ena_radius",false).toBool() )
                aifile += "|EditRadius=" + sets.value("val_radius","0.5").toString();
            if ( sets.value("ena_sigmar",false).toBool() )
                aifile += "|EditSigma=" + sets.value("val_sigmar","0.5").toString();
            if ( sets.value("ena_length",false).toBool() )
                aifile += "|Length=" + sets.value("val_length","0.5").toString();
            if ( sets.value("ena_sigmal",false).toBool() )
                aifile += "|SigmaL=" + sets.value("val_sigmal","0.5").toString();
            if ( sets.value("ena_phi",false).toBool() )
                aifile += "|phi=1"; // [2..87]
            if ( sets.value("ena_a",false).toBool() )
                aifile += "|uca=" + sets.value("val_a","0.5").toString(); // ucb, ucc automatisch
            if ( sets.value("ena_rho",false).toBool() )
                aifile += "|EditRho=" + sets.value("val_rho","0.5").toString();
            if ( sets.value("ena_psi",false).toBool() )
                aifile += "|ucpsi=" + sets.value("val_psi","0.5").toString();
            if ( sets.value("ena_dbeta",false).toBool() )
                aifile += "|EditDbeta=" + sets.value("val_dbeta","0.5").toString();
            if ( sets.value("ena_width",false).toBool() )
                aifile += "|EditDomainSize=" + sets.value("val_width","0.5").toString();
            if ( sets.value("ena_phiwidth",false).toBool() )
                aifile += "|EditAzi=" + sets.value("val_phiwidth","0.5").toString();
            if ( sets.value("ena_displacement",false).toBool() )
                aifile += "|EditDebyeWaller=" + sets.value("val_displacement","0.5").toString();

            std::cerr << "TPV-File=" << qPrintable(f) << std::endl;
            std::cerr << "TPV-Line: " << qPrintable(aifile) << std::endl << std::endl;

            generateAIfiles( calc, aifile, nofft, nthreads );
        }
    }
    else if ( ! aifile.isEmpty() )
    {   // Generate AI files from definition file
        void generateAIfiles( SC_CalcCons *calc, QString aifile, bool nofft, int nthreads );
        if ( flog )
        {
            flog->write("*****************************************" EOL);
            flog->write("* Generate AI files                     *" EOL);
            flog->write("*****************************************" EOL);
        }
        generateAIfiles( calc, aifile, nofft, nthreads );
    }
    else if ( ! autofit.isEmpty() )
    {   // Calculate 2D Fits
        void automaticFit( SC_CalcCons *calc, QString imgfile, QString autofit, int nthreads,
                           QString fnpng );
        if ( flog )
        {
            flog->write("*****************************************" EOL);
            flog->write("* Perform automatic fit                 *" EOL);
            flog->write("*****************************************" EOL);
        }
        automaticFit( calc, imgfile, autofit, nthreads, af_output );
    }
    else if ( ! imgfile.isEmpty() )
    {   // Generate only one image
        void generateSingleImage( SC_CalcCons *calc, QString fnpng, int nthreads );
        if ( flog )
        {
            flog->write("*****************************************" EOL);
            flog->write("* Generate single image file            *" EOL);
            flog->write("*****************************************" EOL);
        }
        generateSingleImage( calc, imgfile, nthreads );
    }

    QDateTime dtEnd = QDateTime::currentDateTime();
    //std::cerr << "DBG: dtStart=" << qPrintable(dtStart.toString()) << std::endl;
    //std::cerr << "DBG: dtEnd  =" << qPrintable(dtEnd.toString()) << std::endl;
    QTime tGlobTime(0,0,0);
    tGlobTime = tGlobTime.addSecs(dtStart.secsTo(dtEnd));
    std::cerr << "Done in " << dtStart.secsTo(dtEnd) << " sec ="
              << qPrintable(tGlobTime.toString("hh:mm:ss")) << std::endl;
    if ( flog )
    {
        flog->write( qPrintable(EOL+QString("Done in %1 sec = %2.").arg(dtStart.secsTo(dtEnd)).arg(tGlobTime.toString("hh:mm:ss.zzz"))+EOL) );
        flog->close();
    }
    if ( fcsv ) fcsv->close();
    return 0; // app.exec();
} /* main() */


QStringList checkWildcard(QString fn)
{
    if ( ! fn.contains("*") && ! fn.contains("%") )
        return QStringList() << fn;
    QStringList rv;
    QFileInfo fi(fn);
    QDir d(fi.absolutePath());
    std::cerr << "checkWildcard(" << qPrintable(fn) << "):" << std::endl;
    std::cerr << "     dir=" << qPrintable(d.absolutePath()) << std::endl;
    std::cerr << "      fn=" << qPrintable(fi.fileName()) << std::endl;
    QFileInfoList fil = d.entryInfoList( QStringList()<<fi.fileName(), QDir::Files );
    foreach ( QFileInfo f, fil )
    {
        rv << f.absoluteFilePath();
        std::cerr << "   found=" << qPrintable(f.absoluteFilePath()) << std::endl;
    }
    return rv;
}


double zuf()
{
    return ((static_cast<double>(rand())/RAND_MAX * 200.0) - 100.0) / 100.0; // -1 .. 1
}


void generateAIfiles( SC_CalcCons *calc, QString aifile, bool nofft, int nthreads )
{
    void saveAIfile( QString fn, double *src, int sx, int sy, int linflg ); // Hilfsroutine zum Speichern

    // das File mit den Berechnungsinfos öffnen
    QFile fcalc;
    QStringList fcalcLines;
    if ( aifile[0] == '@' )
    {   // Spezielle Syntax: Übergabe des Fileinhalts im String
        fcalcLines = aifile.split("@");
    }
    else
    {   // Normaler Filename zum Öffnen
        fcalc.setFileName(aifile);
        if ( !fcalc.open(QIODevice::ReadOnly) )
        {
            std::cerr << qPrintable(aifile + " " + fcalc.errorString()) << std::endl;
            if ( flog )
            {
                flog->write(qPrintable(aifile + " " + fcalc.errorString() + EOL));
            }
            return;
        }
    }

    QString path;

    int imgCount = 0;

    SasCalc_PostProc::inst()->setLogging(false);
    bool fftInputLin=false,
         fftRphiScale=false, fftRphiClip=false, fftRphiClip40=false,    // to avoid compiler warnings
         fftOutScale=false,  fftOutClip=false,  fftOutClip40=false,
         fftOutputSwap=false;
    int  fftRphiSize=0, fftOutSize=0;
#define fftOutRe   0
#define fftOutIm   1
#define fftOutAbs  2
#define fftOutSpec 3
    int fftOutFormat=0;
    int fftOutCutX=0, fftOutCutY=0;

    QString fn, fnpng;
    bool isDef;
    QString dbg;
    // Merker für die "Training Parameter Variation" (TPV)
    bool useTPV = false;
    double tpvBSx0=-1, tpvBSy0=-1, tpvBSxs=-1, tpvBSys=-1;
    int tpvLinesX=0, tpvLinesY=0, tpvLinesXwidth=1, tpvLinesYwidth=1;
    bool tpvAddNoise=false, tpvConvolute=false, tpvAddRphi=false,
         tpvSaveExtra=false, tpvScaleScat=false, tpvGenPNG=false;

    bool factorInLine = false;
    int  factorInLineCount = 0;

    QStringList sl;

    while ( true )
    {
        // Ende-Bedingungen testen
        if ( factorInLine && factorInLineCount==0 ) break;
        if ( aifile[0] == '@' )
        {   // Spezielle Syntax: Übergabe des Fileinhalts im String
            if ( !(factorInLine && factorInLineCount>0) && (fcalcLines.size() == 0) ) break;
        }
        else
        {
            if ( !(factorInLine && factorInLineCount>0) && fcalc.atEnd() ) break;
        }
        // Bearbeiten
        if ( !factorInLine )
        {
            if ( aifile[0] == '@' )
                dbg = fcalcLines.takeFirst();
            else
                dbg = fcalc.readLine().trimmed();
            sl = dbg.split("|",Qt::SkipEmptyParts);
            if ( sl.size() < 2 ) continue;

            if ( sl.size() == 2 )
            {
                std::cerr << "SL: (" << qPrintable(sl.join(",")) << ") '" << qPrintable(dbg) << "'" << std::endl;
                if ( sl[0] == "PATH" )
                {
                    path = sl[1];  // Das ist der Basis-Pfad
#ifdef Q_OS_WIN
                    if ( path.startsWith("/home/wagener") )
                        path.replace("/home/wagener","C:/SimLab");
#else
                    if ( path.startsWith("C:/SimLab") )
                        path.replace("C:/SimLab","/home/wagener");
                    if ( path.startsWith("c:/SimLab") )
                        path.replace("c:/SimLab","/home/wagener");
#endif
                    if ( !path.endsWith("/") ) path += "/";
                    std::cerr << "imgPATH: " << qPrintable(path) << std::endl;
                    continue;
                }
                if ( sl[0] == "COLTBL" )
                {
                    imgColorTbl = sl[1];
                    continue;
                }
                if ( sl[0] == "FFT" && !nofft )
                {
                    static QStringList fftOut = { "OutRe", "OutIm", "OutAbs", "OutSpec" };
                    QStringList fft = sl[1].split(";");
                    while ( fft.size() < 9 ) fft << "0";
                    // FFT|InLog;256;Scale;256;OutAbs;Clip1 Clip4;OutSwap;64;64
                    //      0     1   2     3   4      5           6       7  8
                    fftInputLin   = fft[0] == "InLin";     // InLog
                    fftRphiSize   = fft[1].toInt();
                    fftRphiScale  = fft[2].contains("Scale");
                    fftRphiClip   = fft[2].contains("Clip1");
                    fftRphiClip40 = fft[2].contains("Clip4");
                    fftOutSize    = fft[3].toInt();
                    fftOutFormat  = fftOut.indexOf(fft[4]);
                    fftOutScale   = fft[5].contains("Scale");
                    fftOutClip    = fft[5].contains("Clip1");
                    fftOutClip40  = fft[5].contains("Clip4");
                    fftOutputSwap = fft[6] == "OutSwap";   // OutNoSwap
                    fftOutCutX    = fft[7].toInt();
                    fftOutCutY    = fft[8].toInt();
                    continue;
                }
                if ( sl[0] == "RPHI" )
                {
                    QStringList fft = sl[1].split(";");
                    while ( fft.size() < 2 ) fft << "0";
                    // RPHI|256;Scale
                    //       0   1
                    fftRphiSize   = fft[0].toInt();
                    fftRphiScale  = fft[1].contains("Scale");
                    fftRphiClip   = fft[1].contains("Clip1");
                    fftRphiClip40 = fft[1].contains("Clip4");
                    continue;
                }
            } // if ( sl.size() == 2 )
            if ( sl[0] == "TPV" )
            {
                useTPV = true;
                for ( int i=1; i<sl.size(); i++ )
                {
                    if ( sl[i].startsWith("BS") )
                    {
                        QStringList val = sl[i].mid(2).split(";");
                        if ( val.size() != 4 ) continue;
                        tpvBSx0 = val[0].toDouble();
                        tpvBSy0 = val[1].toDouble();
                        tpvBSxs = val[2].toDouble();
                        tpvBSys = val[3].toDouble();
                        continue;
                    }
                    if ( sl[i].startsWith("LI") )
                    {
                        QStringList val = sl[i].mid(2).split(";");
                        if ( val.size() >= 2 )
                        {
                            tpvLinesX = val[0].toInt();
                            tpvLinesY = val[1].toInt();
                        }
                        if ( val.size() >= 4 )
                        {
                            tpvLinesXwidth = val[2].toInt();
                            tpvLinesYwidth = val[3].toInt();
                        }
                        continue;
                    }
                    if ( sl[i].startsWith("NO") )
                    {
                        tpvAddNoise = true;
                        continue;
                    }
                    if ( sl[i].startsWith("CO") )
                    {
                        tpvConvolute = true;
                        continue;
                    }
                    if ( sl[i].startsWith("RP") )
                    {
                        tpvAddRphi = true;
                        continue;
                    }
                    if ( sl[i].startsWith("IX") )
                    {
                        tpvSaveExtra = true;
                        continue;
                    }
                    if ( sl[i].startsWith("GP") )
                    {
                        tpvGenPNG = true;
                        continue;
                    }
                    if ( sl[i].startsWith("SC") )
                    {
                        tpvScaleScat = true;
                        continue;
                    }
                }
                QStringList parstr;
                if ( tpvBSx0 > 0 )
                    parstr << QString("BS=%1x%2 (%3x%4px)").arg(tpvBSx0).arg(tpvBSy0).arg(tpvBSxs).arg(tpvBSys);
                if ( tpvLinesX > 0 || tpvLinesY > 0 )
                    parstr << QString("Lines: X=%1(%2px) / Y=%3(%4px)").arg(tpvLinesX).arg(tpvLinesXwidth).arg(tpvLinesY).arg(tpvLinesYwidth);
                if ( tpvAddNoise ) parstr << "Noise";
                if ( tpvConvolute ) parstr << "Convolute(TODO)";
                if ( tpvAddRphi ) parstr << QString("Save(r,phi) %1px").arg(fftRphiSize);
                if ( tpvSaveExtra ) parstr << "SaveUnmodified";
                if ( tpvGenPNG ) parstr << "GeneratePNG";
                if ( tpvScaleScat ) parstr << "ScaleScat";
                std::cerr << "TPV params: " << qPrintable(parstr.join("; ")) << std::endl;
                continue;
            } // if ( sl[0] == "TPV" )

            if ( tpvScaleScat ) fftInputLin = true;
            // If alll output images are scaled to [0..1] then the fft gets this input directly
            //  and must not rescale the log values.

        } // if ( !factorInLine )

        isDef = false;
        std::list<std::string> tpvDataStd;
        if ( useTPV )
        {
            if ( factorInLine )
                fn = "_scat_" + QString::number(factorInLineCount-1);
            else
                fn = "_scat_" + sl[0];
            category = "scat";
            iCurCalcArt = 0;
            tpvDataStd.clear();
            for ( int i=1; i<sl.size(); i++ )
            {
                if ( !factorInLine )
                {
                    if ( sl[i].startsWith("#FAC") )
                    {
                        QStringList vv = sl[i].split("=");
                        factorInLine = true;
                        factorInLineCount = vv[1].toInt();
                        fn = "_scat_" + QString::number(factorInLineCount-1);
                        if ( flog )
                            flog->write(qPrintable(QString("TPV Factor %1").arg(factorInLineCount)+EOL));
                    }
                }
                if ( sl[i].startsWith("#") ) continue;
                //QStringList vv = sl[i].split("=");
                //calc->updateParamValue( slCalcArt[iCurCalcArt], vv[0], vv[1].toDouble() );
                tpvDataStd.push_back( sl[i].toStdString() );
                std::cerr << "  TPVDATA: " << sl[i].toStdString() << std::endl;
            }
        }
        else
        {
            calc->tpvPerformRandom(std::list<std::string>()); // Löschen der internen Daten
            // FCC_Latt=10|#Calculation#=0|EditLattice=12.6295|Editx1rel=0.485339|Editx2rel=-0.823865|Editxrel=-0.292733|Edity1rel=-0.602674|Edity2rel=-0.0726763|Edityrel=-0.794671|Editz1rel=0.633427|Editz2rel=0.562107|Editzrel=-0.531795|
            category = sl[0];  // Das wird das Unterverzeichnis

            fn = "";
            bool fnfixed = false;
            for ( int i=1; i<sl.size(); i++ )
            {
                QStringList vv = sl[i].split("=");
                if ( sl[i].startsWith("#Calc") )
                    iCurCalcArt = vv[1].toInt();
                else if ( sl[i].startsWith("#DEF") )
                    isDef = true;
                else if ( sl[i].startsWith("#FN_") )
                {
                    fnfixed = true;
                    fn = vv[0].mid(3) + QString("=%1").arg(vv[1].toDouble(),0,'f',3);
                }
                else
                {
                    calc->updateParamValue( vv[0], vv[1].toDouble() );
                    if ( !fnfixed )
                        fn += QString("_%1=%2").arg(vv[0]).arg(vv[1].toDouble(),0,'f',3); // nur 3 Digits müssen reichen, sonst wird der Name zu lang
                }
            }
            fn = fn.replace("Edit","",Qt::CaseInsensitive);  // Viele Variablen beginnen mit "Edit", das macht den Filenamen lang
            // Jetzt noch die 0 in den Nachkommastellen wegnehmen, macht den Filenamen noch kürzer,
            //  allerdings könnte es sein, dass die Filenamen nicht mehr gleich lang sind
            while ( fn.contains("0_") ) fn=fn.replace("0_","_");
            fn = fn.replace("._","_");
            while ( fn.endsWith("0") ) fn.chop(1);
            if ( fn.endsWith(".") ) fn.chop(1);
            if ( isDef ) fn += "_DEF";  // Marker for default orientation
        }
        fnpng = path + category + "/" + fn.mid(1);  // fn startet mit "_"
        if ( ! QDir().exists(path + category) )
        {
            QDir().mkpath( path + category );
            if ( useTPV )
            {
                // category = "scat" in diesem Fall
                if ( tpvGenPNG )
                    QDir().mkpath( path + "scat_png" );
                if ( tpvAddRphi )
                {
                    QDir().mkpath( path + "rphi" );
                    if ( tpvGenPNG )
                        QDir().mkpath( path + "rphi_png" );
                }
                if ( (fftOutSize > 0 /*|| fftRphiSize > 0*/) && !nofft )
                {
                    QDir().mkpath( path + "fft"  );
                    if ( tpvGenPNG )
                        QDir().mkpath( path + "fft_png"  );
                }
                if ( tpvSaveExtra )
                {
                    QDir().mkpath( path + "scat_org" );
                    if ( tpvGenPNG )
                        QDir().mkpath( path + "scat_org_png" );
                }
            }
        }

        // Berechnen
        calc->prepareCalculation( false/*fromfit*/, false/*use1d*/ );
        if ( useTPV )
        {
            if ( (factorInLine && factorInLineCount > 1) || !factorInLine )
            {
                std::string rv = calc->tpvPerformRandom(tpvDataStd);
                //QString str = QString::fromStdString(rv);
                //std::cerr << "TPVinfo: " << rv << std::endl;
                if ( flog )
                    flog->write((EOL+rv).c_str());
                    //flog->write(qPrintable(EOL+str));
            }
            else if ( flog )
            {
                flog->write(EOL "TPV: no random values for this image." EOL);
                //std::cerr << "TPVinfo: no random values for this image." << std::endl;
            }
        }
        calc->doCalculation( nthreads, false /*ignNewSwitch*/ );

        std::cerr << qPrintable(fn) << " -> " << calc->higResTimerElapsed(SC_CalcCons::htimBoth) << std::endl;
        if ( flog )
            flog->write(qPrintable(QString("Calctime: %1ms").arg(calc->higResTimerElapsed(SC_CalcCons::htimBoth))+EOL));

        double *src = calc->data();
        int sx = calc->maxX() - calc->minX();
        int sy = calc->maxY() - calc->minY();
        int len = sx * sy;
        int bsx0 = calc->currentParamValue("BeamPosX");   // Position
        int bsy0 = calc->currentParamValue("BeamPosY");
        int bsdx = 2; // TODO  Größe
        int bsdy = 2;

        if ( useTPV )
        {   // Nachbearbeitung
            QString str = "";
            if ( tpvScaleScat )
            {   // Auch das Orginalbild soll auf [0 .. 1] skaliert werden
                calc->scaleIntensity(true); // In diesem Fall immer auf LOG skalieren
            }
            if ( tpvSaveExtra )
            {   // Bild auch ohne Nachbearbeitung speichern
                saveAIfile( fnpng, src, sx, sy, tpvGenPNG?13:3 );    // Orginalbild vor der Manipulation
            }
            if ( tpvConvolute )
            {   // Convolution
                // TODO
            }
            if ( tpvAddNoise )
            {   // Noise
                double *d = src;
                double minz=200, maxz=-200;
                for ( int i=0; i<len; i++, d++ )
                {
                    double z = zuf();
                    if ( z < minz ) minz = z;
                    if ( z > maxz ) maxz = z;
                    if ( tpvScaleScat )                 // if the data values are in the range [0 .. 1]
                        *d = *d + (z/100.0) * sqrt(*d); // then a random number [-1 .. 1] is too big.
                    else
                        *d = *d + z * sqrt(*d);
                }
                str += QString("AddNoise [%1,%2]").arg(minz).arg(maxz) + EOL;
            }
            if ( tpvBSx0 > 0 )
            {   // Beamstop
                bsx0 = bsx0 * ( 1.0 + tpvBSx0 * zuf() );
                bsy0 = bsy0 * ( 1.0 + tpvBSy0 * zuf() );
                bsdx = bsdx * ( 1.0 + tpvBSxs * zuf() );
                bsdy = bsdy * ( 1.0 + tpvBSys * zuf() );
                for ( int x=calc->minX(); x<calc->maxX(); x++ )
                    if ( x >= bsx0-bsdx && x <= bsx0+bsdx )
                        for ( int y=calc->minY(); y<calc->maxY(); y++ )
                            if ( y >= bsy0-bsdy && y <= bsy0+bsdy )
                            {
                                int ii = (-calc->minX() + x) + sx*(-calc->minY() + y);
                                if ( ii >= 0 && ii < len ) src[ii] = 0; // 1e-20;
                            }
                str += QString("AddBeamstop Pos=(%1,%2) Size=(%3,%4)").arg(bsx0).arg(bsy0).arg(bsdx).arg(bsdy) + EOL;
            }
            if ( tpvLinesX > 0 )
            {
                int xc = (calc->maxX() - calc->minX() + 1) / tpvLinesX;
                // zuf() = -1 .. +1
                str += "AddLinesX at";
                for ( int i=0; i<tpvLinesX; i++ )
                {
                    int xp = calc->minX() + i * xc + fabs(zuf()) * xc;
                    str += QString(" %1").arg(xp);
                    for ( int x=xp; x<xp+tpvLinesXwidth; x++ )
                        for ( int y=calc->minY(); y<calc->maxY(); y++ )
                        {
                            int ii = (-calc->minX() + x) + sx*(-calc->minY() + y);
                            if ( ii >= 0 && ii < len ) src[ii] = 0; // 1e-20;
                        }
                }
                str += EOL;
            }
            if ( tpvLinesY > 0 )
            {
                int yc = (calc->maxY() - calc->minY() + 1) / tpvLinesY;
                // zuf() = -1 .. +1
                str += "AddLinesY at";
                for ( int i=0; i<tpvLinesY; i++ )
                {
                    int yp = calc->minY() + i * yc + fabs(zuf()) * yc;
                    str += QString(" %1").arg(yp);
                    for ( int y=yp; y<yp+tpvLinesYwidth; y++ )
                        for ( int x=calc->minX(); x<calc->maxX(); x++ )
                        {
                            int ii = (-calc->minX() + x) + sx*(-calc->minY() + y);
                            if ( ii >= 0 && ii < len ) src[ii] = 0; // 1e-20;
                        }
                }
                str += EOL;
            }
            if ( flog )
                flog->write(qPrintable(str));
        }

        saveAIfile( fnpng, src, sx, sy, tpvGenPNG?10:0 );    // Berechnetes Bild

        if ( (fftOutSize > 0 && !nofft) || fftRphiSize > 0 )
        {   // noch die FFT berechnen ...
            static double *data = nullptr;

            if ( ! fftInputLin )
            {   // Das Bild wird logarithmisch verarbeitet

                const double *lsrc = calc->data();
                if ( data == nullptr )
                    data = new double[sx*sy];   // Pro Lauf sind die Bilder immer gleich
                double vlogmin=1e26, vlogmax=0;
                for ( int i=0, y=calc->minY(); y<calc->maxY(); y++ )
                    for ( int x=calc->minX(); x<calc->maxX(); x++, i++ )
                        if ( (calc->minX() < 0 && x != 0 && y != 0) ||
                             (calc->minX() == 0) )
                        {   // Wenn der Nullpunkt im Bild ist, dann ist dort i.d.R. ein Extremwert (z.B. -7.4e+25 bei FCC)
                            // Bei den anderen Bildern (r,phi oder FFT) ist immer der Nullpunkt in einer Ecke
                            if ( lsrc[i] > 0 )
                            {
                                if ( log10(lsrc[i]) < vlogmin )
                                    vlogmin = log10(lsrc[i]);
                                else if ( log10(lsrc[i]) > vlogmax )
                                    vlogmax = log10(lsrc[i]);
                            }
                        }
                //qDebug() << "FFT scale log" << vlogmin << vlogmax;
                for ( int i=0; i<sx*sy; i++ )
                    if ( lsrc[i] > 0 )
                        data[i] = (log10(lsrc[i]) - vlogmin) / (vlogmax-vlogmin);
                    else
                        data[i] = 0;
            }

            double *rphi = nullptr;
            if ( fftRphiSize > 0 )
            {
                rphi = SasCalc_PostProc::inst()->generateRphi(calc->minX(), calc->maxX(),
                                                              calc->minY(), calc->maxY(),
                                                              bsx0, bsy0, fftRphiSize,
                                                              (data==nullptr)?calc->data():data,
                                                              fftRphiScale,    // scale [0,1]
                                                              fftRphiClip,     // clip [1e-6, 1]
                                                              fftRphiClip40 ); // clip40 [40%,100%] -> [0,1]

                SasCalc_PostProc::inst()->scaleAndClipData( rphi, fftRphiSize*fftRphiSize, false, false, false, true/*Log*/ );

                saveAIfile( fnpng, rphi, fftRphiSize, fftRphiSize, tpvGenPNG?11:1 ); // r,phi speichern
            }

            if ( fftOutSize > 0 )
            {
                double *fftindata;
                int x0, x1, y0, y1;
                if ( rphi != nullptr )
                {
                    fftindata = rphi; // r,phi
                    x0 = 0;
                    x1 = fftRphiSize;
                    y0 = 0;
                    y1 = fftRphiSize;
                }
                else if ( data != nullptr )
                {
                    fftindata = data; // Log
                    x0 = calc->minX();
                    x1 = calc->maxX();
                    y0 = calc->minY();
                    y1 = calc->maxY();
                }
                else
                {
                    fftindata = calc->data(); // Lin
                    x0 = calc->minX();
                    x1 = calc->maxX();
                    y0 = calc->minY();
                    y1 = calc->maxY();
                }
                src = SasCalc_PostProc::inst()->calculateIFFT(false, // Backward
                                                              x0, x1, y0, y1,
                                                              fftindata,
                                                              fftOutSize,
                                                              static_cast<SasCalc_PostProc::_outType>(fftOutFormat),
                                                              fftOutScale,    // scale [0,1]
                                                              fftOutClip,     // clip [1e-6, 1]
                                                              fftOutClip40,   // clip40 [40%,100%] -> [0,1]
                                                              fftOutputSwap,
                                                              fftOutCutX,fftOutCutY);

                if ( src != nullptr )
                    saveAIfile( fnpng, src, fftOutSize, fftOutSize, tpvGenPNG?12:2 ); // fft speichern
            }
        } // if ( fftImgSize > 0 && !nofft )

        if ( (++imgCount) % 1000 == 0 )
        {   // Flush der Logs, damit ich extern etwas sehen kann
            std::cerr << std::flush;
            if ( flog ) flog->flush();
        }

        if ( useTPV && factorInLine )
        {
            std::cerr << "   TPV & Factor " << factorInLineCount << std::endl;
            if ( --factorInLineCount <= 0 ) break;
        }

        // TEST
        //if ( imgCount > 1 ) break;

    } // while ( !fcalc.atEnd() )
    return;
} /* generateAIfiles() */

void saveAIfile( QString fn, double *src, int sx, int sy, int linflg )
{   // linflg: 0|10=png, 1|11=r/phi, 2|12=fft, 3|13=Scat_Org
    //         1x=PNG zusätzlich

    static widImage *img = nullptr;
    if ( imgColorTbl[0] == '*' )
    {   // Spezielle Kennung zum Speichern in Textformat (AIgen)
        static const int prec = 6;
        //static const int digi = 2;
        if ( linflg == 1 || linflg == 11 )
            fn.replace("scat","rphi");
        else if ( linflg == 2 || linflg == 12 )
            fn.replace("scat","fft");
        else if ( linflg == 3 || linflg == 13 )
            fn.replace("scat/","scat_org/"); // Only directory, not filename!
        QFile fo(fn+".spr");
        if ( fo.open(QIODevice::WriteOnly) )
        {
            //xystring[0]:=concat('    ',IntToStr(points),'    ',IntToStr(points),' Start pixel = (       1       1)');
            //for i:=1 to points do begin     (* rows *)
            //   xystring[i]:=' ';
            //   for j:=0 to points-1 do      (* columns *)
            //   xystring[i]:=concat(xystring[i],' ',FloatToStrF(calcIma[i,j],form,prec,digi));
            //end;
            fo.write( qPrintable(QString("    %1    %2 Start pixel = (       1       1)").arg(sx).arg(sy)+EOL) );
            for ( int y=0; y<sy; y++ )
            {
                QString line = " ";
                for ( int x=0; x<sx; x++ )
                    line += QString(" %1").arg(src[x+sx*y],prec+7,'e',prec);    // -1.123456e+01
                fo.write( qPrintable(line+EOL) );
            }
            fo.close();
            if ( flog )
                flog->write(qPrintable(QString("Save %1").arg(fn)+EOL));
            // Jetzt den Filenamen (bzw. den Pfad) so anpassen, dass das Bild (unten)
            // in einem anderen Unterpfad gespeichert wird.
            int i = fn.lastIndexOf("/");
            if ( i < 0 ) i = fn.lastIndexOf("\\");
            fn.insert(i,"_png");
            std::cerr << qPrintable(fn) << std::endl;
        }
        else if ( flog )
            flog->write(qPrintable(QString("ERROR %1 (%2)").arg(fo.errorString()).arg(fo.fileName())+EOL));
        if ( linflg < 10 ) return;
    }
    //else  Bei TPV auch das farbige Image speichern (macht das Testen einfacher)
    {   // Normale Speicherung als PNG
        if ( img == nullptr ) img = new widImage();
        bool useLin = false;
        if ( linflg == 1 || linflg == 11 )      // rphi
            useLin = true;
        else if ( linflg == 2 || linflg == 12 ) // fft
            useLin = true;
        //else if ( linflg == 3 || linflg == 13 ) // orgimg
        img->setConfig(imgColorTbl, imgSwapH, imgSwapV, imgRot, imgZoom, useLin);
        //img->setData( calc->minX(), calc->maxX(), calc->minY(), calc->maxY(), src );
        img->setData( 0, sx, 0, sy, src );
        img->saveImage( fn+".png" );    // auch in Farbe speichern möglich
        if ( flog )
            flog->write(qPrintable(QString("Save %1.png").arg(fn)+EOL));
    }

    /*if ( flog )       TODO
        {
            double min, max;
            img->getVarScaling( min, max );
            std::cerr << "  MinVal=" << min << ",  MaxVal=" << max << std::endl;
            flog->write(qPrintable(QString("  MinVal=%1,  MaxVal=%2").arg(min).arg(max)));
        }*/
}


void generateSingleImage( SC_CalcCons *calc, QString fnpng, int nthreads )
{
    if ( !fnpng.endsWith(".png",Qt::CaseInsensitive) ) fnpng += ".png";
    iCurCalcArt = 0;

    // Berechnen
    calc->prepareCalculation( false/*fromfit*/, false/*use1d*/ );
    calc->doCalculation( nthreads, false /*ignNewSwitch*/ );
    std::cerr << qPrintable(fnpng) << " -> " << calc->higResTimerElapsed(SC_CalcCons::htimBoth) << std::endl;
    if ( flog )
        flog->write(qPrintable(QString("%1 -> %2ms").arg(fnpng).arg(calc->higResTimerElapsed(SC_CalcCons::htimBoth))+EOL));

    widImage *img = new widImage();
    img->setConfig(imgColorTbl, imgSwapH, imgSwapV, imgRot, imgZoom, false);
    img->setData( calc->minX(), calc->maxX(), calc->minY(), calc->maxY(), calc->data() );
    img->saveImageGray( fnpng );
} /* generateSingleImage() */



static SasCalc_SimplexFit2D *fitClass = nullptr;
static _param2fitval p2f;  // Merker für alle zu fittenden Variablen

void updateLogList( QString msg )
{
    if ( msg.startsWith("@REP=") ) return;
    if ( flog )
    {
        flog->write(qPrintable(msg+EOL));
        if ( msg.startsWith("iter=") ) flog->flush();
    }
    if ( fcsv!=nullptr && msg.startsWith("----- rtol=") )
    {   // ----- rtol=1.12315 <= ftol=1e-015
        QString tmp = msg.mid(11);
        tmp.truncate(tmp.indexOf(" <="));
        tmp.replace(".",",");
        fcsv->write(qPrintable(" ; "+tmp+EOL));
    }
    // Damit etwas auf die Konsole ausgegeben wird, aber nicht zuviel:
    if ( msg.startsWith("FQS") ) return;
    if ( msg.startsWith("amotry") ) return;
    if ( !noConsolOutput ) std::cerr << qPrintable(msg) << std::endl; // to avoid " " around the output
}

bool myProgressLogging( char *msg )
{
    updateLogList( QString(msg) );
    return false; // keine Abbruch-Möglichkeit, da keine GUI
}


void printAllParameter( SC_CalcCons *calc, QFile *f )
{
    if ( f == nullptr )
    {
        std::cerr << "printAllParameter - no file" << std::endl;
        return;
    }
    std::cerr << "printAllParameter - " << qPrintable(f->fileName()) << std::endl;
    f->write("------------------------------ Parameters Start" EOL);

    QStringList meta = calc->getCalcPtr()->guiLayoutNeu();
    for ( int m=0; m<meta.size(); m++ )
    {
        QStringList sl = meta[m].split(";");
        while ( sl.size() < 5 ) sl << " ";  // tooltip;default sind optional

        // TODO: Neue Struktur des guiLayout()
        // Each element must be in the form "type;kenn;typespec;tooltip;default" with:
        //  type     = C:Selection, N:Double, I:Integer, T:Toggle, O:DoubleOut
        //             Zweites Zeichen ist 'F' für Fittable oder '-' für nicht fittable
        //  kenn     = internal parameter name to connect to correct gui element
        //  typespec = C:Selection : "...|...|..."  (required)
        //             N:Double    : "frac|min|max|unit"  (optional, default: "2|-10000|+10000|")
        //             I:Integer   : "min|max|unit"  (optional, default: "-10000|+10000|")
        //             T:Toggle    : (empty)
        //             O:DoubleOut : (empty)
        //  tooltip  = this is the tooltip set to both prompt label and inputfield (optional)
        //  default  = the default value (optional)

        QString isfit  = sl[0][1] == 'F' ? "Fittable" : "";
        QString key = sl[1].trimmed();
        QStringList typval = sl[2].split("|",Qt::SkipEmptyParts);
        paramConsHelper *par = calc->params[key];
        if ( par == nullptr )
        {
            f->write(qPrintable(QString("%1 = unknown" EOL).arg(key)));
            continue;
        }
        switch ( sl[0][0].toLatin1() )
        {
        case 'C':   // ComboBox  : "...|...|..." mit den jeweiligen Werten
            f->write(qPrintable(QString("%1 = %2  ").arg(par->key).arg(par->value.number)));
            if ( typval.size() > par->value.number )
                f->write(qPrintable(typval[par->value.number]));
            f->write(EOL);
            break;
        case 'N':   // Zahlenwert: "frac|min|max|unit" mit Nachkommastellen und Grenzwerten (optional)
        case 'I':   // Integer-Zahlenwert: "min|max|unit"
            f->write(qPrintable(QString("%1 = %2  %3" EOL).arg(par->key).arg(par->value.number).arg(isfit)));
            break;
        case 'T':   // CheckBox  : "tog"
            f->write(qPrintable(QString("%1 = %2" EOL).arg(par->key,par->value.flag?"true":"false")));
            break;
        }
    }
    f->write("------------------------------ Parameters End" EOL);
}



double doFitStart( SC_CalcCons *calc, widImage *img, int nthreads,
                   double inpFitStepSize, double tolerance,
                   int inpFitRepetitions, int inpFitMaxIter,
                   int inpFitBorder, int inpFitBStop,
                   double &timeForAll, int &loopsForAll, int &imgGenForAll )
{
    if ( fitClass == nullptr )
        fitClass = new SasCalc_SimplexFit2D( calc );

    fitClass->setImageInfo( img->xmin(), img->xmax(), img->ymin(), img->ymax(),
                            img->getFileInfos()->centerX, img->getFileInfos()->centerY,
                            img->dataPtr() );

    calc->prepareCalculation( true/*fromfit*/, false/*use1d*/ );

    QHash<QString,_fitLimits*>::const_iterator it;

    timeForAll   = 0;
    loopsForAll  = 0;
    imgGenForAll = 0;

    for ( int rep=0; rep<inpFitRepetitions; rep++ )
    {
        if ( rep > 0 )
        {
            updateLogList( QString("rep=%4: %1 ms / %2 loops / %3 img").arg(fitClass->higResTimerElapsed)
                           .arg(fitClass->repetitions).arg(fitClass->numImgCalc).arg(rep) );
            it = p2f.constBegin();
            while ( it != p2f.constEnd() )
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
        QString retinfo;
        fitClass->doSimplexFit2D( nthreads,
                                  inpFitStepSize * 0.01, // wie im Pascalprog TODO ???
                                  inpFitMaxIter,
                                  tolerance,
                                  inpFitBorder,
                                  inpFitBStop,  // or -1 to use the masking
                                  static_cast<progressLogging>(&myProgressLogging),
                                  &p2f, retinfo );
        timeForAll += fitClass->higResTimerElapsed;
        loopsForAll += fitClass->repetitions;
        imgGenForAll += fitClass->numImgCalc;
        if ( ! retinfo.isEmpty() && flog!=nullptr )
        {
            flog->write(qPrintable("Simplex 2d Fit - critical:"+retinfo+EOL));
            return -1;
        }
        if ( flog ) flog->flush();
        if ( fcsv ) fcsv->flush();
        if ( fitClass->aborted ) break;
        // Jetzt das bisherige Ergebnis als neuen Startwert nehmen, orgval bleibt aber wie vorher
        it = p2f.constBegin();
        while ( it != p2f.constEnd() )
        {
            if ( it.value()->used )
                it.value()->fitstart = it.value()->fitres;
            ++it;
        }

    }

    if ( fitClass->aborted )
        updateLogList("*** Aborted ***");

    updateLogList( QString("rep=%4: %1 ms / %2 loops / %3 img").arg(fitClass->higResTimerElapsed)
                   .arg(fitClass->repetitions).arg(fitClass->numImgCalc).arg(inpFitRepetitions) );
    updateLogList("=============================");
    it = p2f.constBegin();
    double fitMeanChangePercent = 0.0;
    int fmc = 0;
    while ( it != p2f.constEnd() )
    {
        if ( it.value()->used )
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

            it.value()->fitstart = it.value()->fitres;
            it.value()->orgval   = it.value()->fitres;
        }
        ++it;
    }
    fitMeanChangePercent = fitMeanChangePercent / fmc;
    if ( fcsv )
    {
        QString tmp = QString(" ;  ; %1").arg(fitMeanChangePercent);
        tmp.replace(".",",");
        fcsv->write(qPrintable(tmp+EOL));
    }
    updateLogList( QString("Overall: %1 ms / %2 loops / %3 img / MeanDif %4%").arg(timeForAll)
                   .arg(loopsForAll).arg(imgGenForAll).arg(fitMeanChangePercent) );
    return fitMeanChangePercent;
}


widImage *local_OpenMeasFile(QString fn, bool &iskws)
{
    iskws = false;
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
            img = SC_ReadData::readImageSasCrystal( addImage, fn );
        }
    }
    if ( img == nullptr )
    {
        // Routines in sc_readdata.cpp
        if ( fn.endsWith(".tif",Qt::CaseInsensitive) ||
            fn.endsWith(".tiff",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageTiff( addImage, fn );
        else if ( fn.endsWith(".dat",Qt::CaseInsensitive) ||
                 fn.endsWith(".data",Qt::CaseInsensitive) )
        {
            img = SC_ReadData::readImageKWSData( addImage, fn );
            iskws = true;
        }
        else if ( fn.endsWith(".edf",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageEDF( addImage, fn );
        else if ( fn.endsWith(".spr",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageSpreadsheet( addImage, fn );
#ifndef NOHDF5
        else if ( fn.endsWith(".h5",Qt::CaseInsensitive) ||
                 fn.endsWith(".hdf",Qt::CaseInsensitive) ||
                 fn.endsWith(".hdf5",Qt::CaseInsensitive) )
            img = SC_ReadData::readImageHDF5( addImage, fn, true/*onlyOneImage*/, false/*swap vert*/ );
        else if ( fn.endsWith(".nxs",Qt::CaseInsensitive) )     // 2D SANS (ILL)
            img = SC_ReadData::readImageHDF5( addImage, fn, true/*onlyOneImage*/, true/*swap vert*/ );
#endif
        else if ( fn.endsWith(".xml",Qt::CaseInsensitive) ||    // 2D SANS (ILL)
                 fn.endsWith(".csv",Qt::CaseInsensitive) ||    // 2D SANS (ILL)
                 fn.endsWith(".txt",Qt::CaseInsensitive) )     // 2D SANS (ILL)
            img = SC_ReadData::readImage2dSans( addImage, fn );
    }
    return img;
}



/**
 * @brief automaticFit - performs the automatic 2d simplex fit with recipe
 * @param calc      = pointer to calculation class
 * @param imgfile   = input imagefilepath/name if only one image is fitted
 *                    empty if multiple images are fitted (given in recipe)
 * @param autofit   = filepath/name with recipe
 * @param nthreads  = number of threads to be used (0=GPU if available)
 * @param fnpng     = output image filepath/name if only one image is fitted
 *                    filepath for all other output file if multiple images are fitted
 */
void automaticFit(SC_CalcCons *calc, QString imgfile, QString autofit, int nthreads, QString fnpng)
{
    QFile finp( autofit );
    if ( ! finp.open(QIODevice::ReadOnly) )
    {
        if ( flog ) flog->write(qPrintable(autofit+": "+finp.errorString()+EOL));
        std::cerr << qPrintable(autofit) << ": " << qPrintable(finp.errorString()) << std::endl;
        return;
    }

    // Basispfad für alle Outputs, immer im gleichen Verzeichnis wie das Inputfile
    QString basePath;
    if ( fnpng.isEmpty() )
        basePath = QFileInfo(autofit).absolutePath() + "/";
    else
    {
        basePath = fnpng;
        if ( !basePath.endsWith("/") ) basePath += "/";
        QDir().mkpath(basePath);
    }
    std::cerr << "Used basepath: " << qPrintable(basePath) << std::endl;

    // Jetzt werden alle Files außer dem Inputfile gelöscht
    QStringList sl = QDir(basePath).entryList( QDir::Files | QDir::NoDotAndDotDot );
    sl.removeOne(QFileInfo(autofit).fileName());
    if ( flog ) sl.removeOne(QFileInfo(*flog).fileName());  // Das Logfile ist geöffnet und kann nicht gelöscht werden
    std::cerr << "Remove: " << qPrintable(sl.join(", ")) << std::endl;
    foreach ( QString f, sl )
    {
        if ( ! QFile::remove(basePath+f) )
            std::cerr << "Remove fail: " << qPrintable(basePath+f) << std::endl;
            // Beim Fehler nur eine Meldung, kein Abbruch!
    }

    widImage *imgZiel = nullptr;
    int bs_x=0, bs_y=0;

    if ( !imgfile.isEmpty() )
    {
        bool iskws;
        imgZiel = local_OpenMeasFile(imgfile, iskws);
        if ( imgZiel == nullptr )
        {
            if ( flog ) flog->write(qPrintable("unable to open "+imgfile+EOL));
            std::cerr << "Unable to open data file: " << qPrintable(imgfile) << std::endl;
            return;
        }
        if ( iskws )
        {
            SC_ReadData::findBeamCenter( imgZiel, bs_x, bs_y );
            imgZiel->addMetaInfo( "BeamPosX", QString::number(bs_x) );
            imgZiel->addMetaInfo( "BeamPosY", QString::number(bs_y) );
            calc->updateParamValue( "BeamPosX", bs_x - imgZiel->myWidth()/2. );
            calc->updateParamValue( "BeamPosY", bs_y - imgZiel->myHeight()/2. );
            if ( flog )
                flog->write(qPrintable(QString("AutoBeamstop: in data (%1 / %2), in calc (%3 / %4)").arg(bs_x).arg(bs_y)
                                           .arg(bs_x - imgZiel->myWidth()/2.).arg(bs_y - imgZiel->myHeight()/2.)+EOL));
            std::cerr << "ImageInfos: X=" << imgZiel->xmin() << " .. " << imgZiel->xmax()
                      << ", Y=" << imgZiel->ymin() << " .. " << imgZiel->ymax()
                      << ", BS=" << imgZiel->getFileInfos()->centerX << " / " << imgZiel->getFileInfos()->centerY
                      << std::endl;
        }
        else
            std::cerr << "ImageInfos: X=" << imgZiel->xmin() << " .. " << imgZiel->xmax()
                      << ", Y=" << imgZiel->ymin() << " .. " << imgZiel->ymax()
                      << ", BS not searched" << std::endl;
    }

    if ( fcsv ) fcsv->write("Variables ; rtol ; MeanDiffPercentage" EOL);

    // Die globale Variable 'flog' wird auch bei der Berechnung des Fits mit verwendet. Beim Aufruf hier ist das
    // eine globale Logdatei, die per Parameter geschaltet werden kann. Diese nehme ich hier lokal in 'globLog'
    // und nutze die 'flog' um für jedes Bild der Serie eine eigene Logdatei zu schreiben, die bleibt dann kleiner.
    QFile *globLog = flog;
    flog = nullptr;

    // Jetzt die Struktur p2f füllen mit allen fitbaren Variablen aus den Definitionen der angegebenen Methode
    QStringList slParams = calc->paramsForMethod( false, false, true );
    foreach (QString p, slParams)
    {
        bool cnt;
        double min, max;
        //qDebug() << "FIT: create new" << p;
        _fitLimits *fl = new _fitLimits;
        fl->orgval = fl->fitstart = calc->currentParamValue( p );
        fl->fitType  = _fitTypes::fitNone;
        fl->fitvalid = false;
        fl->used = false;
        if ( calc->limitsOfParamValue( p, min, max, cnt ) )
        {
            if ( cnt )
                fl->fitType = _fitTypes::fitCbs;
            else
                fl->fitType = _fitTypes::fitNumeric;
            fl->min  = min;
            fl->max  = max;
            //qDebug() << "DEF:" << p << min << max;
        }
        else
        {
            //qDebug() << "DEF:" << p << "no limits ***********************************";
            fl->fitType = _fitTypes::fitNumeric;
            fl->min  = 0;
            fl->max  = 1000;
        }
        p2f.insert( p, fl );
    } // foreach (QString p, slParams)

    // Jetzt kommen die Min/Max/Start Werte Modifikationen aus den Settings
    // (wenn das Programm im manuellen Modus auf diesem Rechner schon mal lief)
    // Diese Daten werdn aus den globalen Settings geholt, nicht aus dem Parameterfile
    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup( "Fit-Limits" );
    QStringList slKeys = sets.allKeys();
    foreach (QString k, slKeys)
    {
        if ( ! calc->isCurrentParameterValid(k) ) continue;
        _fitLimits *fl = p2f.value(k,nullptr);
        if ( fl == nullptr ) continue;
        QStringList slVal = sets.value(k,"0:0:0:0").toString().split(":");
        fl->used = slVal[0].toInt() != 0;
        fl->min  = slVal[1].toDouble();
        fl->max  = slVal[2].toDouble();
        fl->fitType = static_cast<_fitTypes>(slVal[3].toInt());
        fl->fitvalid = false;  // fitresult not stored!
        //std::cerr << "REG: " << qPrintable(k) << "=" << fl->min << " .. " << fl->max << std::endl;
    }

    QString kenn = QDate::currentDate().toString("-yyyyMMdd-");     // Als Default
    bool firstScanOfFile = true;
    QString dirMask = "*";
    QStringList usedImages; // wird während des Abarbeitens geleert, daher kein Index nötig
    QFile *ftex = nullptr;
    QDateTime startZeit;
    int fitBorder=0;    // keine Pixel am Rand ignorieren
    int fitBeamstop=0;  // 0=keinen BS in der Mitte, -1=Eckpixel maskiert, >0=BS in der Mitte ausblenden

    typedef struct
    {
        double anf, end;
    } _globValues;
    QHash<QString/*param*/,_globValues*> param2values;

    double globTimeForAllFits=0;    // Globale Summen über alle Dateien (falls Multifile)
    int    globLoopsForAllFits=0;
    int    globImgGenForAllFits=0;

    QStringList slInputLines;

    while ( true )
    {   // Multifile Fit Loop

        double timeForAllFits=0;        // Globale Summen über alle Fits eines Files
        int    loopsForAllFits=0;
        int    imgGenForAllFits=0;

        QString curImgZiel = "";

        while ( ! finp.atEnd() )
        {
            QString line = finp.readLine().trimmed();
            int pos = line.indexOf("#");
            if ( pos >= 0 ) { line.truncate(pos); line=line.trimmed(); }
            if ( line.isEmpty() ) continue;         // Ignore empty lines

            // Keine Kommentare in das LaTeX File schreiben lassen, dann gibt es auch keine
            //  Probleme mit den Umlauten (bzw. UTF-8)
            if ( firstScanOfFile ) slInputLines << line;

            std::cerr << "AutoFit: " << qPrintable(line) << std::endl;
            if ( globLog != nullptr && firstScanOfFile ) globLog->write(qPrintable("Autofit: "+line+EOL));

            if ( line.startsWith("EOF") ) break; // Hier ist das File zu Ende
            // Das erspart beim Testen das Auskommentieren der restlichen Zeilen

            if ( firstScanOfFile )
            {   // Einige Einträge werden nur beim ersten Durchlauf interpretiert

                // Scale: TODO

                if ( line.startsWith("GlobLog:") )
                {   // "GlobLog: <filename>  Globales Logfile für alle Informationen

                    // Hier ist 'globLog' ein globales Logfile von gesamten Konsolenprogramm.
                    // 'flog' dagegen ist die lokale Logdatei für die Berechnungen.
                    /*
                    if ( flog != nullptr ) flog->close();
                    flog = new QFile( basePath + line.mid(8).trimmed() );
                    if ( ! flog->open(QIODevice::WriteOnly) )
                    {
                        if ( globLog )
                            globLog->write(qPrintable(flog->fileName()+": "+flog->errorString()+EOL));
                        flog->deleteLater();
                        flog = nullptr;
                    }
                    else if ( globLog )
                        globLog->write(qPrintable("Write to "+flog->fileName()+EOL));
                    */

                    // Hier verwendet für den LaTeX Output
                    startZeit = QDateTime::currentDateTime();
                    if ( globLog ) globLog->write(qPrintable("Start: "+startZeit.toString("dd.MMM.yyyy hh:MM:ss")+EOL));
                    ftex = new QFile( basePath + "latex-output.tex" );
                    if ( ! ftex->open(QIODevice::WriteOnly) )
                    {
                        std::cerr << qPrintable(ftex->fileName()) << ": " << qPrintable(ftex->errorString()) << std::endl;
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
                                                          "\\begin{document}" EOL };
                        foreach ( QString s, slTextBase )
                        {
                            ftex->write(qPrintable(s+EOL));
                        }
                        //ftex->write(qPrintable("\\section{Simplex 2D Fit - "+QDate::currentDate().toString("dd. MMM. yyyy")+"}"+EOL));
                        if ( globLog ) globLog->write(qPrintable("Use LaTeX-Output "+ftex->fileName()+" as a test." EOL));
                    }
                    continue;
                }

                if ( line.startsWith("Param:") )
                {   // "Param: <filepath>
                    if ( globLog ) globLog->write(qPrintable("Load Paramfile "+line.mid(6).trimmed()+EOL) );
                    calc->loadParameter( line.mid(6).trimmed() );

                    if ( doPrintParameter ) printAllParameter( calc, globLog );

                    continue;
                }

                if ( line.startsWith("DirMask:") )
                {   // "DirMask: *.dat" Filemask für die Verzeichnis-Suche
                    dirMask = line.mid(8).trimmed();
                    if ( globLog ) globLog->write(qPrintable("Set directory mask to "+dirMask+EOL) );
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
                    std::cerr << "AutoFit: DirUp " << usedImages.size() << " Files found" << std::endl; // << usedImages.first();
                    if ( globLog ) globLog->write(qPrintable(QString("Dir (Up) %1 files found.").arg(fil.size())+EOL) );
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
                    std::cerr << "AutoFit: DirDown " << usedImages.size() << " Files found" << std::endl; // << usedImages.first();
                    if ( globLog ) globLog->write(qPrintable(QString("Dir (Down) %1 files found.").arg(fil.size())+EOL) );
                    continue;
                }

                if ( line.startsWith("File:") )
                {   // "File: <filepath>" Angabe von einzelnen Dateien (mehrfacher Eintrag möglich)
                    usedImages.append( line.mid(5).trimmed() );
                    //qDebug() << "AutoFit: UsedFiles" << usedImages;
                    if ( globLog ) globLog->write(qPrintable("File added "+line.mid(5).trimmed()+EOL) );
                    continue;
                }

                if ( line.startsWith("Limits:") )
                {   // "Limits: Editdom2; 0; 500" - Angabe von Min/Max Werten für spezielle Variablen
                    QStringList slNames = line.mid(7).trimmed().split(";");
                    for ( int i=0; i<slNames.size(); i++ ) slNames[i] = slNames[i].trimmed();
                    QHash<QString,_fitLimits*>::const_iterator it = p2f.constBegin();
                    while ( it != p2f.constEnd() )
                    {
                        if ( slNames[0].contains(it.key(),Qt::CaseInsensitive) )
                        {
                            if ( slNames.size() >= 2 ) it.value()->min = slNames[1].toDouble();
                            if ( slNames.size() >= 3 ) it.value()->max = slNames[2].toDouble();
                            //qDebug() << "INP:" << it.key() << it.value()->min << it.value()->max;
                            break;
                        }
                        ++it;
                    }
                    continue;
                }

            } // if ( firstScanOfFile )

            if ( line.startsWith("Info:") )
            {   // "Info: <text>" information text for user
                if ( globLog != nullptr && !firstScanOfFile ) globLog->write(qPrintable("Autofit: "+line+EOL));
                continue;
            }

            // Threads: TODO

            if ( line.startsWith("Use:") )
            {   // "Use: Base, I0"
                bool oneUsed = false;
                if ( fcsv ) fcsv->write(qPrintable(line.mid(4)+EOL));
                QStringList slNames = line.mid(4).trimmed().split(",");
                for ( int i=0; i<slNames.size(); i++ ) slNames[i] = slNames[i].trimmed();
                QHash<QString,_fitLimits*>::const_iterator it = p2f.constBegin();
                while ( it != p2f.constEnd() )
                {

                    it.value()->used = slNames.contains(it.key(),Qt::CaseInsensitive);
                    it.value()->orgval = it.value()->fitstart = calc->currentParamValue( it.key() );
                    oneUsed |= it.value()->used;
                    //qDebug() << it.key() << it.value()->used;
                    if ( ! param2values.contains(it.key()) && it.value()->used )
                    {   // Nur die verwendeten Parameter werden aufgelistet
                        _globValues *gv = new _globValues;
                        gv->anf = it.value()->orgval;
                        param2values.insert( it.key(), gv );
                    }
                    ++it;
                }
                if ( !oneUsed )
                {
                    if ( globLog ) globLog->write(qPrintable("ERROR no parameter used: "+line+EOL));
                    if ( flog ) flog->close();
                    flog = globLog;
                    if ( ftex ) ftex->close();
                    qDebug() << "ERROR" << slNames << p2f.keys();
                    finp.close();
                    return;
                }
                if ( globLog )
                {
                    globLog->write(qPrintable("USE Names="+slNames.join(", ")+EOL));
                    globLog->write(qPrintable("USE p2fkeys="+p2f.keys().join(", ")+EOL));
                    globLog->write(qPrintable("USE param2values="+param2values.keys().join(", ")+EOL));
                    globLog->flush();
                }
                continue;
            } // if "Use:..."

            if ( line.startsWith("Border:") )  // NEU 25.07.2024
            {   // "Border: 0; 0; 0"
                // Border size; Beamstop half-size; ignore pixel
                QStringList slCmd = line.mid(7).trimmed().split(";");
                for ( int i=0; i<slCmd.size(); i++ )
                    slCmd[i] = slCmd[i].trimmed();
                while ( slCmd.size() < 3 ) slCmd << "0";
                fitBorder   = slCmd[0].toInt(); // Pixel am Rand ignorieren
                fitBeamstop = slCmd[1].toInt(); // 0=keinen BS in der Mitte, >0=BS in der Mitte ausblenden
                if ( slCmd[2].toInt() != 0 ) fitBeamstop = -1; // -1=Eckpixel maskiert
                continue;
            }

            // SetVarAnf: - noch experimentell und nur in der GUI
            // SetVarEnd: - noch experimentell und nur in der GUI
            // SetVar:    - noch experimentell und nur in der GUI

            else if ( line.startsWith("Fit:") )
            {   // "Fit: Rep=10; Stp=3.0; Iter=20; Tol=0.0010; Diff<5; Kenn=..."
                QStringList slCmd = line.mid(4).trimmed().split(";");
                for ( int i=0; i<slCmd.size(); i++ )
                    slCmd[i] = slCmd[i].trimmed();
                //qDebug() << slCmd;
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
                    {   // Diese Kennung war zum speichern der Logfiles in der interaktiven Session.
                        kenn = slCmd[i].mid(5);
                        kenn = kenn.replace( "@D", QDate::currentDate().toString("yyyyMMdd") );
                        if ( usedImages.size() > 0 )
                            kenn = kenn.replace("@F", QFileInfo(usedImages.first()).baseName() );
                        else
                            kenn = kenn.replace("@F", "" );
                        if ( !kenn.endsWith("-") ) kenn += "-";  // Damit man den Zähler erkennt.
                    }
                }
#ifndef USEREPETITIONS
                int rep = 1;    // ist so am einfachsten...
#endif
                double timeForAll=0;        // Summen über alle Loops
                int    loopsForAll=0;
                int    imgGenForAll=0;

                if ( imgfile.isEmpty() &&
                     usedImages.size() > 0 &&
                     curImgZiel != usedImages.first() )
                {   // Jetzt muss das erste Image geladen werden ...

                    // TODO: Zudem die richtigen Metadaten neben dem File lesen

                    bool iskws;
                    imgZiel = local_OpenMeasFile(usedImages.first(), iskws);
                    //imgZiel = SC_ReadData::readImageKWSData( addImage, usedImages.first() );
                    if ( imgZiel == nullptr )
                    {
                        if ( globLog ) globLog->write(qPrintable("unable to open "+usedImages.first()+EOL));
                        if ( flog ) flog->close();
                        flog = globLog;
                        if ( ftex ) ftex->close();
                        return;
                    }
                    curImgZiel = usedImages.first();
                    if ( globLog )
                        globLog->write(qPrintable("Imgfile "+curImgZiel+EOL));
/*
                    if ( flog ) flog->close();
                    flog = new QFile( basePath+QFileInfo(curImgZiel).baseName()+"_calc.log");
                    if ( ! flog->open(QIODevice::WriteOnly) )
                    {
                        if ( globLog ) globLog->write(qPrintable("unable to open logfile ("+flog->fileName()+"). This run is without logfile." EOL));
                        flog->deleteLater();
                        flog = nullptr;
                    }
                    else if ( globLog )
                        globLog->write(qPrintable("Logfile "+flog->fileName()+EOL));
*/
                    if ( iskws )
                    {
                        SC_ReadData::findBeamCenter( imgZiel, bs_x, bs_y );
                        imgZiel->addMetaInfo( "BeamPosX", QString::number(bs_x) );
                        imgZiel->addMetaInfo( "BeamPosY", QString::number(bs_y) );
                        calc->updateParamValue( "BeamPosX", bs_x - imgZiel->myWidth()/2. );
                        calc->updateParamValue( "BeamPosY", bs_y - imgZiel->myHeight()/2. );
                        if ( globLog )
                            globLog->write(qPrintable(QString("AutoBeamstop: in data (%1 / %2), in calc (%3 / %4)").arg(bs_x).arg(bs_y)
                                                       .arg(bs_x - imgZiel->myWidth()/2.).arg(bs_y - imgZiel->myHeight()/2.)+EOL));
                    }
                    if ( ftex )
                    {
                        QString tmpfn = basePath + QFileInfo(usedImages.first()).baseName() + "_org.png";
                        imgZiel->saveImage( tmpfn );
                        if ( !firstScanOfFile ) ftex->write( EOL "\\clearpage" EOL );
                        ftex->write( qPrintable("\\section*{\\small{"+usedImages.first()+"}}"+EOL) );
                        ftex->write( "\\begin{longtable}{cc}" EOL );
                        ftex->write( "\\bf{Datafile} & \\bf{Fit output} \\\\ \\hline" EOL );
                        ftex->write( qPrintable("\\includegraphics[width=0.45\\textwidth]{"+
                                                QFileInfo(tmpfn).fileName()+"} &"+EOL) );
                        ftex->flush();
                    }

                    std::cerr << "ImageInfos: X=" << imgZiel->xmin() << " .. " << imgZiel->xmax()
                              << ", Y=" << imgZiel->ymin() << " .. " << imgZiel->ymax()
                              << ", BS=" << imgZiel->getFileInfos()->centerX << " / " << imgZiel->getFileInfos()->centerY
                              << ", iskws=" << iskws << std::endl;
                    if ( globLog ) globLog->flush();
                } // if usedImages > 0

                for ( int maxloop=0; maxloop<20; maxloop++ ) // keine Endlosschleife...
                {
                    if ( flog ) flog->close();
                    flog = new QFile(basePath + kenn + QString("%1").arg(maxloop+1,2,10,QChar('0')) + "_calc.log");
                    if ( ! flog->open(QIODevice::WriteOnly) )
                    {
                        if ( globLog ) globLog->write(qPrintable("unable to open logfile ("+flog->fileName()+"). This run is without logfile." EOL));
                        flog->deleteLater();
                        flog = nullptr;
                    }
                    else if ( globLog )
                    {
                        globLog->write(qPrintable("Logfile "+flog->fileName()+EOL));
                        globLog->flush();
                    }

                    double timeForOne;      // Summen über alle Repetitions/Iterations
                    int    loopsForOne;
                    int    imgGenForOne;
                    double fitMeanChangePercent = doFitStart( calc, imgZiel, nthreads,
                                                              stp, tol, rep, iter,
                                                              fitBorder, fitBeamstop,
                                                              timeForOne, loopsForOne, imgGenForOne );
                    timeForAll += timeForOne;
                    loopsForAll += loopsForOne;
                    imgGenForAll += imgGenForOne;
                    timeForAllFits += timeForOne;
                    loopsForAllFits += loopsForOne;
                    imgGenForAllFits += imgGenForOne;
                    globTimeForAllFits += timeForOne;
                    globLoopsForAllFits += loopsForOne;
                    globImgGenForAllFits += imgGenForOne;
                    std::cerr << "AutoFit: loop " << maxloop+1 << ", MeanDif=" << fitMeanChangePercent
                             << ", MaxDif =" << maxDif << std::endl;
                    if ( fitMeanChangePercent < 0 )
                    {   // Errormessage
                        finp.close();
                        if ( globLog ) globLog->write("Autofit MeanChangePercent < 0" EOL);
                        if ( flog ) flog->close();
                        flog = globLog;
                        if ( ftex ) ftex->close();
                        return;
                    }
                    if ( flog )
                    {
                        flog->write(qPrintable(EOL+QString("AutoFit: loop=%1, MeanDiff=%2, MaxDiff=%3").arg(maxloop+1).arg(fitMeanChangePercent).arg(maxDif)+EOL+EOL));
                        flog->flush();
                    }
                    if ( fitMeanChangePercent < maxDif )
                    {
                        QHash<QString,_fitLimits*>::const_iterator it = p2f.constBegin();
                        while ( it != p2f.constEnd() )
                        {   // Nur für die verwendeten Parameter!
                            if ( it.value()->used )
                                param2values[it.key()]->end = it.value()->fitres;
                            ++it;
                        }
                        break;
                    }
                }
                //if ( _bAbbruch ) break;
                continue;
            } // if "Fit:..."

        } // while !atEnd

        updateLogList( QString("AutoFit: Sum over all fits: %1 sec / %2 loops / %3 img").arg(timeForAllFits/1000)
                       .arg(loopsForAllFits).arg(imgGenForAllFits) );

        if ( imgfile.isEmpty() &&
             usedImages.size() > 0 )
        {
            QString tmpfn = basePath + QFileInfo(usedImages.first()).baseName() + "_out";
            calc->prepareCalculation( true/*fromfit*/, false/*use1d*/ );
            calc->doCalculation( nthreads, false /*ignNewSwitch*/ );
            widImage *img = new widImage();
            img->setConfig(imgColorTbl, imgSwapH, imgSwapV, imgRot, imgZoom, false);
            img->setData( calc->minX(), calc->maxX(), calc->minY(), calc->maxY(), calc->data() );
            img->saveImage( tmpfn+".png" );
            calc->saveParameter( tmpfn+".ini" );
            if ( flog )
            {
                flog->write( qPrintable("Save Image  to "+tmpfn+".png"+EOL) );
                flog->write( qPrintable("Save Params to "+tmpfn+".ini"+EOL+EOL) );
                flog->flush();
            }
            if ( globLog )
            {
                globLog->write( qPrintable("Save Image to "+tmpfn+".png and params to .ini"+EOL) );
                globLog->flush();
            }
            if ( ftex )
            {
                ftex->write( qPrintable("\\includegraphics[width=0.45\\textwidth]{"+
                                        QFileInfo(tmpfn).fileName()+".png} \\\\"+EOL) );
                ftex->write( "\\end{longtable}" EOL );
                ftex->write(qPrintable(QString("\\centerline{Statistics: %1 ms running time / %2 iterations / %3 images generated}")
                                       .arg(timeForAllFits).arg(loopsForAllFits).arg(imgGenForAllFits)+EOL) );
                ftex->write( "\\begin{longtable}{|C{3cm}|C{4cm}|C{4cm}|}" EOL );
                ftex->write( "\\hline\\rowcolor{rowcolor} {\\bf Parameter} & {\\bf Startvalue} & {\\bf Fit-Resultvalue} \\\\" EOL );
                ftex->write( "\\endfirsthead" EOL );
                ftex->write( "\\hline\\rowcolor{rowcolor} {\\bf Parameter} & {\\bf Startvalue} & {\\bf Fit-Resultvalue} \\\\" EOL );
                ftex->write( "\\endhead" EOL );
                ftex->write( "\\hline" EOL );
                QStringList gvKeys = param2values.keys();
                gvKeys.sort();  // Sortierte Schlüssel ausgeben
                foreach ( QString k, gvKeys )
                {
                    QString n=k;
                    n.replace("_","\\_");
                    ftex->write(qPrintable(QString("%1 & %2 & %3 \\\\ \\hline")
                                           .arg(n).arg(param2values[k]->anf).arg(param2values[k]->end)+EOL));
                }
                ftex->write( "\\end{longtable}" EOL );
                ftex->flush();
                param2values.clear();
            }
            firstScanOfFile = false;
            finp.seek(0);
            usedImages.takeFirst();
            if ( usedImages.size() == 0 ) break;
        } // if ( imgfile.isEmpty() && usedImages.size() > 0 )
        else
            break;

    } // while endless

    if ( flog ) flog->close();
    flog = globLog;

    finp.close();

    QTime tGlobTime(0,0,0);
    tGlobTime = tGlobTime.addMSecs(globTimeForAllFits);
    if ( flog )
    {
        QDateTime endZeit = QDateTime::currentDateTime();
        flog->write(qPrintable("Ende: "+endZeit.toString("dd.MMM.yyyy hh:MM:ss")+EOL));
    }
    if ( ftex )
    {
        ftex->write(qPrintable(EOL+QString("\\centerline{Sum over all images: %1 running time / %2 iterations / %3 images generated}")
                               .arg(tGlobTime.toString("hh:mm:ss.zzz")).arg(globLoopsForAllFits).arg(globImgGenForAllFits)+EOL) );
        ftex->write( EOL "\\clearpage" EOL );
        ftex->write( "\\begin{lstlisting}[frame=single, xleftmargin=0.5cm, xrightmargin=0.5cm]" EOL ); // % language=Python,
        foreach ( QString s, slInputLines )
        {
            ftex->write( qPrintable(s+EOL) );
        }
        ftex->write( "\\end{lstlisting}" EOL );
        ftex->write( EOL "\\end{document}" EOL );
        ftex->close();
    }

    if ( ! imgfile.isEmpty() )
    {   // Einzelfile, also jetzt erst das generierte Bild speichern
        calc->prepareCalculation( true/*fromfit*/, false/*use1d*/ );
        calc->doCalculation( nthreads, false /*ignNewSwitch*/ );
        std::cerr << qPrintable(fnpng) << " -> " << calc->higResTimerElapsed(SC_CalcCons::htimBoth) << std::endl;
        if ( flog )
        {
            flog->write(qPrintable(QString("%1 -> %2ms").arg(fnpng).arg(calc->higResTimerElapsed(SC_CalcCons::htimBoth))));
            flog->write(EOL EOL "All Parameters:" EOL);
            QStringList par = p2f.keys();
            par.sort();                 // Sortiert ist besser zu vergleichen
            foreach (QString k, par)
            {
                if ( p2f.value(k)->used )
                {
                    QString k1 = k;
                    k1.replace("_","\\_");
                    flog->write(qPrintable(QString("%1 = %2").arg(k1).arg(p2f.value(k)->fitres)+EOL) );
                }
            }
        }

        widImage *img = new widImage();
        img->setConfig(imgColorTbl, imgSwapH, imgSwapV, imgRot, imgZoom, false);
        img->setData( calc->minX(), calc->maxX(), calc->minY(), calc->maxY(), calc->data() );
        img->saveImageGray( fnpng );

        fnpng.replace(".png",".ini");
        calc->saveParameter( fnpng );
    }
} /* automaticFit() */



void performTimeTest( QString par, QString cfg, QString out )
{
    // Outputfile vorbereiten. Falls das schief geht, gehen wir hier schon raus.
    QFile fout(out);
    if ( ! fout.open(QIODevice::Append) )
        if ( ! fout.open(QIODevice::WriteOnly) )
        {
            std::cerr << "Error open file " << qPrintable(out+": "+fout.errorString()) << std::endl;
            return;
        }
#ifdef Q_OS_WIN
    QString hname = qEnvironmentVariable("COMPUTERNAME");
#else
    QString hname = qEnvironmentVariable("HOSTNAME");
#endif
    fout.write( qPrintable("% Generated with sas_scatter2Cons "+QDateTime::currentDateTime().toString("dd.MMM.yyyy hh:mm:ss")+" on "+hname+EOL) );
    fout.write( qPrintable("%       Parameterfile: "+par+EOL) );
    fout.write( qPrintable("% TimeTest configfile: "+cfg+EOL) );
    fout.write( "Threads & HKLmax & Quadrants & Switch & Min/ Mean/ Max (Prep) in ms" EOL );
    fout.close();
    // Einlesen der TimeTest-Konfiguration
    QSettings sets(cfg,QSettings::IniFormat);
    int nthreads = static_cast<int>(std::thread::hardware_concurrency());
    QVector<int> threads;
    if ( sets.value("thread0",true).toBool() ) threads << 0;    // TODO: ist die GPU da?
    if ( sets.value("thread1",true).toBool() ) threads << 1;
    if ( sets.value("thread4",true).toBool() ) threads << 4;
    if ( sets.value("threadM",true).toBool() && nthreads>4 ) threads << nthreads;
    QVector<int> hklmax;
    if ( sets.value("hkl2",false).toBool() ) hklmax << sets.value("hkl2val",2).toInt();
    if ( sets.value("hkl3",true ).toBool() ) hklmax << sets.value("hkl3val",3).toInt();
    if ( sets.value("hkl4",false).toBool() ) hklmax << sets.value("hkl4val",4).toInt();
    if ( sets.value("hkl5",true ).toBool() ) hklmax << sets.value("hkl5val",5).toInt();
    QVector<int> quadrants;
    if ( sets.value("usequadrants",true).toBool() )
    {
        if ( sets.value("useq1",false).toBool() ) quadrants << 1;
        if ( sets.value("useq2",false).toBool() ) quadrants << 2;
        if ( sets.value("useq4",true ).toBool() ) quadrants << 4;
    }
    else
        quadrants << 0;
    typedef enum { swBoth, swNew, swOld } _newSwitch;
    _newSwitch nswitch = static_cast<_newSwitch>(sets.value("newswitch",swBoth).toInt());
    _newSwitch nswcurr;
    //ui->togSwitchModFN->setChecked(sets.value("switchmodfn",true).toBool());
    //ui->togEnaUpdates->setChecked(sets.value("enaupdates",false).toBool());
    //ui->togSaveImages->setChecked(sets.value("saveimages",false).toBool());
    //ui->togSaveFile->setChecked(sets.value("savefile",false).toBool());
    //on_togSaveFile_toggled(ui->togSaveFile->isChecked()); // Sicherheitshalber, falls das set keine Änderung macht
    int numLoops = sets.value("loops",10).toInt();
    //ui->inpSaveFilename->setText(sets.value("filename","").toString());
    //ui->inpComment->setText(sets.value("comment","").toString());

    SC_CalcCons *calc = new SC_CalcCons;
    calc->loadParameter( par );

    if ( calc->currentParamValue("HKLmax") == 1 )
    {   // Kennung, dass der Wert hier nicht verwendet wird
        hklmax.clear();
        hklmax << 1;
    }

    // ***** Loop 1 ***** Quadrants
    foreach ( int quadr, quadrants )
    {
        if ( quadr > 0 )
        {
            calc->updateParamValue( "RadioButtonQ1", quadr == 1 );
            calc->updateParamValue( "RadioButtonQ2", quadr == 2 );
            calc->updateParamValue( "RadioButtonQ4", quadr == 4 );
        }

        // ***** Loop 2 ***** HKLmax
        for ( int h=0; h<hklmax.size(); h++ )
        {
            calc->updateParamValue( "HKLmax", hklmax[h] );

            // ***** Loop 3 ***** Threads
            for ( int t=0; t<threads.size(); t++ )
            {
                if ( threads[t] == 0 && !calc->gpuAvailable() ) continue;

                // ***** Loop 4 ***** NewSwitch
                if ( nswitch/*Soll*/ == swBoth )
                    nswcurr = swOld;
                else
                    nswcurr = nswitch;
                for ( int sw=0; sw<2; sw++ )
                {   // Es können maximal 2 Durchläufe sein...

                    double sumCalc=0;
                    double minCalc=10000000;
                    double maxCalc=0;
                    double sumPrep=0;
                    double minPrep=10000000;
                    double maxPrep=0;
                    std::cerr << "TimingTest: Threads=" << threads[t] << ", HKLmax=" << hklmax[h]
                              << ", Quadrants=" << quadr << ", Switch=" << nswcurr << ", Loops="
                              << numLoops << std::endl;
                    for ( int i=0; i<numLoops; i++ )
                    {
                        // Berechnen
                        calc->prepareCalculation( false/*fromfit*/, false/*use1d*/ );
                        calc->doCalculation( threads[t], nswcurr == swOld );

                        double c = calc->higResTimerElapsed(SC_CalcCons::htimCalc);
                        double p = calc->higResTimerElapsed(SC_CalcCons::htimPrep);
                        sumCalc += c;
                        sumPrep += p;
                        if ( p < minPrep ) minPrep = p;
                        else if ( p > maxPrep ) maxPrep = p;
                        if ( c < minCalc ) minCalc = c;
                        else if ( c > maxCalc ) maxCalc = c;
                    } // i < numLoops
                    std::cerr << "ERG Calc=" << minCalc << " / " << sumCalc/numLoops << " / " << maxCalc
                              << ", Prep=" << minPrep << " / " << sumPrep/numLoops << " / " << maxPrep << std::endl;
                    // LaTeX Ausgabe mit den Rundungen (am Ende auf der Konsole oder in eine Datei)
                    QFile f(out);
                    if ( ! f.open(QIODevice::Append) )
                    {
                        std::cerr << "Error open file " << qPrintable(out+": "+f.errorString()) << std::endl;
                        return;
                    }
                    f.write( QString("%1 & %2 & %8 & %9 & %3/ %4/ %5 (%6)%7")
                                .arg(threads[t]).arg(hklmax[h])
                                .arg(minCalc,0,'f',3).arg(sumCalc/numLoops,0,'f',3).arg(maxCalc,0,'f',3)
                                .arg(sumPrep/numLoops,0,'f',3).arg(EOL).arg(quadr).arg(nswcurr==swOld?"old":"new").toLatin1() );
                    f.close();

                    //typedef enum { swBoth, swNew, swOld } _newSwitch;
                    if (nswitch/*Soll*/ == swBoth &&
                        nswcurr/*Ist */ == swOld)
                    {   // Jetzt muss die Schleife nochmals wiederholt werden
                        nswcurr = swNew;
                    }
                    else
                        sw = 10;
                } // for sw
            } // for t
        } // for h
    } // for q
} /* performTimeTest() */


void waitForChatbot(SC_CalcCons *calc, QString cbfile, /*QString param,*/ int nthreads, bool keep)
{
    iCurCalcArt = 0;
    QString basepath = QFileInfo(cbfile).absolutePath();

    while ( true )
    {   // Endlosschleife. Wird vom rufenden Programm via Terminate beendet
        if ( ! QFile::exists(cbfile) )
        {   // Warten...
            QThread::sleep( 2 /*sec*/ );
            continue;
        }
        // default-Parameter laden - wurde schon vor dem Aufruf gemacht.
        //if ( !param.isEmpty() ) calc->loadParameter( param );
        // cbfile einlesen ...
        QString fnpng="";
        QFile fcb(cbfile);
        if ( fcb.open(QIODevice::ReadOnly) )
        {
            // TODO: im cbfile sollte ein Filename für das Image sein
            while ( true )
            {
                QString line = fcb.readLine().trimmed();
                if ( line.isEmpty() )
                {
                    if ( fcb.atEnd() ) break;
                    continue;
                }
                int pos;
                if ( (pos=line.indexOf('#')) >= 0 )
                {   // Kommentare sind erlaubt
                    line.truncate(pos);
                    line = line.trimmed();
                    if ( line.isEmpty() ) continue;
                }

                QStringList sl = line.split("=");
                if ( sl.size() != 2 ) continue;         // Kein '=' bedeutet falsche Syntax, ignorieren

                if ( sl[0].startsWith("ena_") ) continue;
                // Diese Flags brauche ich hier nicht (mehr)

                if ( ! sl[0].startsWith("val_") ) continue; // Falsche Syntax

                sl[1].replace("\"", "");

                if ( keep )
                    std::cerr << "READ '" << qPrintable(sl[0]) << "' = '" << qPrintable(sl[1]) << "'" << std::endl;

                if ( sl[0].startsWith("val_numimg",Qt::CaseInsensitive) ) continue;
                // Erstmal ignorieren ...

                if ( sl[0].contains("Dist",Qt::CaseInsensitive) ) continue; // war von früher

                if ( sl[0].startsWith("val_outPath",Qt::CaseInsensitive) )
                {   // Jetzt wird der Pfad angegeben, eventuell auch das File.

                    fnpng = sl[1].trimmed();
                    if ( fnpng.length() < 4 )
                    {   // Jetzt ist wohl nur ein '.' angegeben
                        fnpng = basepath + "/" + QFileInfo(cbfile).baseName() + ".png";
                    }
                    if ( fnpng.contains("[path]") )
                    {
                        fnpng.replace( "[path]", QFileInfo(cbfile).baseName()+".png" );
                    }
                    if ( fnpng.startsWith("/home/user/") )
                    {
                        fnpng = basepath + "/" + QFileInfo(cbfile).baseName() + ".png";
                    }
                    std::cerr << "OUTfile: (" << QFileInfo(fnpng).isRelative() << ") "
                              << qPrintable(QFileInfo(fnpng).absolutePath()) << std::endl;
                    if ( QFileInfo(fnpng).isRelative() )
                    {
                        fnpng = QFileInfo(basepath+"/"+fnpng).absoluteFilePath();
                        //std::cerr << "OUTfile: neu " << qPrintable(fnpng) << std::endl;
                    }
                    std::cerr << "OUTfile: " << qPrintable(fnpng) << std::endl;
                }
                else if ( sl[1][0] != '[' )
                {   // Wenn im Chatbot nichts angegeben war, dann wird der Wert "[value]" geschrieben,
                    // dann sollte der Default bleiben, sonst gibt es kein Image.
                    // Die Keys wurden zm Schreiben etwas modifiziert, das muss wieder rückgängig gemacht werden.
                    //+ if ( key.startsWith("VAx")  ) prtkey = key.mid(1); // Das 'V' stört bei der Ausgabe
                    //- else if ( key.startsWith("Edit") ) prtkey = key.mid(4); // Das 'Edit' stört auch
                    //- else if ( key.startsWith("CheckBox") ) prtkey = key.mid(8);
                    //+ else if ( key.startsWith("ComboBox") ) prtkey = "CB"+key.mid(8);
                    //- else if ( key.startsWith("RadBut") ) prtkey = key.mid(6);
                    //+ else if ( key.startsWith("RadioButton") ) prtkey = "RB"+key.mid(11);

                    double val = sl[1].trimmed().toDouble();
                    if ( sl[1].contains("false",Qt::CaseInsensitive) ) val = 0;
                    else if ( sl[1].contains("true",Qt::CaseInsensitive) ) val = 1;

                    QString key = sl[0].mid(4).trimmed();
                    if ( key.startsWith("Ax") ) key = "V"+key;
                    else if ( key.startsWith("CB") ) key = "ComboBox"+key.mid(2);
                    else if ( key.startsWith("RB") ) key = "RadioButton"+key.mid(2);
                    bool rv = calc->updateParamValue( key, val );
                    if ( !rv ) rv = calc->updateParamValue( "Edit"+key, val );
                    if ( !rv ) rv = calc->updateParamValue( "CheckBox"+key, val );
                    if ( !rv ) rv = calc->updateParamValue( "RadBut"+key, val );
                    // Bei return false ist der Parameter unbekannt
                    if ( !rv && flog != nullptr )
                    {
                        flog->open(QIODevice::Append);
                        flog->write(qPrintable("Parameter "+sl[0].mid(4)+" unknown."+EOL));
                        flog->close();
                    }
                }
            } // while true -> bis eof lesen
            fcb.close();
            if ( fnpng.isEmpty() )
            {
                fnpng = cbfile + ".png"; // Default, falls nichts angegeben wird
                if ( flog )
                {
                    flog->open(QIODevice::Append);
                    flog->write(qPrintable("No outputfile given, use same name as inputfile: "+fnpng+EOL));
                    flog->close();
                }
            }
        } // if open ok
        else
        {   // Lesefehler, obwohl Datei vorhanden ist...
            std::cerr << qPrintable(cbfile) << " -> " << qPrintable(fcb.errorString()) << std::endl;
            std::cerr << "Program exits." << std::endl;
            if ( flog )
            {
                flog->open(QIODevice::Append);
                flog->write(qPrintable(cbfile+" -> "+fcb.errorString()+EOL));
                flog->write("Program exits." EOL);
                flog->close();
            }
            if ( !keep ) fcb.remove(); // Versuch, wird aber wahrscheinlich scheitern ....
            break; // Jetzt nicht mehr weiter machen
        }
        // Berechnen
        calc->prepareCalculation( false/*fromfit*/, false/*use1d*/ );
        calc->doCalculation( nthreads, false /*ignNewSwitch*/ );
        if ( flog )
        {
            flog->open(QIODevice::Append);
            flog->write(qPrintable(QString("%1 -> %2ms").arg(fnpng).arg(calc->higResTimerElapsed(SC_CalcCons::htimBoth))+EOL));
            flog->close();
        }

        widImage *img = new widImage();
        img->setConfig(imgColorTbl, imgSwapH, imgSwapV, imgRot, imgZoom, false);
        img->setData( calc->minX(), calc->maxX(), calc->minY(), calc->maxY(), calc->data() );
        QString dbgsave;
        if ( ! img->saveImageColor( fnpng, imgColorTbl, dbgsave ) )
        {   // fnpng enthält schon den Pfad (wird oben sichergestellt)
            std::cerr << "Error saving image " << qPrintable(dbgsave) << std::endl;
            if ( flog )
            {
                flog->open(QIODevice::Append);
                flog->write(qPrintable(QString("Error saving image: ")+dbgsave+EOL));
                flog->close();
            }
        }
        else
        {
            std::cerr << "IMG: " << qPrintable(fnpng) << " -> " << calc->higResTimerElapsed(SC_CalcCons::htimBoth) << "ms" << std::endl;
            if ( flog )
            {
                flog->open(QIODevice::Append);
                flog->write(qPrintable(QString("Save image: ")+dbgsave+EOL));
                flog->close();
            }
        }

        if ( keep ) break;

        if ( ! fcb.remove() )
        {
            std::cerr << "Error removing inputfile. Program exits." << std::endl;
            if ( flog )
            {
                flog->open(QIODevice::Append);
                flog->write("Error removing file. Program exits." EOL);
                flog->close();
            }
            break;
        }
    } // endless
} /* waitForChatbot() */


/*
TODO:
- äußere Schleife (Rep) erstmal weglassen oder Sicherstellen, dass es eine Erweiterung der Iterationen ist
- Bilderserien zum anfitten verfügbar machen (Wildcarts/sortierung/feste Namen)
- Optimierung: Berechnungsmaske erstellen und übergeben.
- Editsigma muss noch gefittet werden können

Sequenz mit Förster, kam sehr gut an:

Use: Base, I0
Fit: Rep=2; Stp=10; Iter=10; Tol=1e-15; Diff<2; Kenn=-@D

Use: uca, Base, ucc, ucb, I0
Fit: Rep=2; Stp=10; Iter=10; Tol=1e-15; Diff<2; Kenn=-@D

Use: Base, EditRadius, I0
Fit: Rep=2; Stp=10; Iter=10; Tol=1e-15; Diff<2; Kenn=-@D

Use: Base, I0, EditDebyeWaller
Fit: Rep=2; Stp=10; Iter=10; Tol=1e-15; Diff<2; Kenn=-@D

EditDebyeWaller=1, da der Wert doch weggelaufen ist

Use: Editdom3, Base, I0, EditDebyeWaller, Editdom2, Editdom1
Fit: Rep=2; Stp=10; Iter=10; Tol=1e-15; Diff<2; Kenn=-@D

Use: Editdom3, Base, I0, EditDebyeWaller, Editdom2, Editdom1 + Editsigma
Fit: Rep=2; Stp=10; Iter=10; Tol=1e-15; Diff<2; Kenn=-@D

*/
