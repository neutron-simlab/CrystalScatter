#include "sc_calcCons.h"
#include <QCoreApplication>
#include <QCommandLineParser>
#include <QPixmap>
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
#include "sc_maingui.h"  // mittels #ifndef CONSOLENPROG alles bis auf Definitionen ausgeblendet

#include <QDebug> // nur zum Test

//#define USEGENCOLORTBL


//#define USE_INTERNAL_DEFAULTS
// Versuch, aber noch nicht verwendet...
#define INTDEF_PARAMFILE "sas_scatter2/Test B-0-6 (BCT).ini"
#define INTDEF_IMGFILE   "20211206 - Fit Rezept/Daten/46424-AWM3-OX-SF01-CDCl3-B-2-95.dat"
#define INTDEF_THREADS   "4"
#define INTDEF_LOGFILE   "20211206 - Fit Rezept/AutoFit-Multifile/AutoFit-Test" // -Win.log oder -Gpu.log
#define INTDEF_CSVFILE   "20211206 - Fit Rezept/AutoFit-Test" // -Win.csv oder -Gpu.csv
#define INTDEF_AUTOFIT   "20211206 - Fit Rezept/AutoFit-Multifile/autoFitCommands.txt"
#define INTDEF_AFOUTPUT  "20211206 - Fit Rezept/AutoFit-Multifile"
#define INTDEF_METHOD    "bct"


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


QStringList slCalcArt;
int iCurCalcArt;

QString category;

QFile *flog=nullptr, *fcsv=nullptr;

QString imgColorTbl;
bool    imgSwapH, imgSwapV;
int     imgRot, imgZoom;
bool    noConsolOutput;

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

    QCoreApplication::setApplicationName("sas-scatter");
    QCoreApplication::setApplicationVersion("1.0");

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
                                      "Filepath and name for the automatic fit routine (*.txt). If this is given, the AI-File is ignored and the imagefile is an input dataset.",
                                      "autofit" );
    parser.addOption(autofitOption);

    QCommandLineOption outputOption( QStringList() << "o" << "output",
                                     "Filepath and name of the single output image of the automatic fit routine or the path to save the multifile outputs.",
                                     "output" );
    parser.addOption(outputOption);

    QCommandLineOption methodOption( QStringList() << "m" << "method",
                                     "Method to be used for single image calculation and 2D-Fit.",
                                     "method" );
    parser.addOption(methodOption);

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
                                  " This parameter can be given more than once.",
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
    QString parMethod = parser.value(methodOption);
    bool nofft        = parser.isSet(nofftOption);
    noConsolOutput    = parser.isSet(noConsolOption);

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
    if ( !ok || imgRot < 0 || imgRot > 3 ) imgRot = 0;
    imgZoom  = parser.value(imgZoomOption).toInt(&ok);
    if ( !ok || (imgZoom!=1 && imgZoom!=2 && imgZoom!=4) ) imgZoom = 1;
    if ( parser.value(imgColorOption).startsWith("grey",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames()[0];
    else if ( parser.value(imgColorOption).startsWith("glow",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames()[1];
    else if ( parser.value(imgColorOption).startsWith("earth",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames()[2];
    else if ( parser.value(imgColorOption).startsWith("temp",Qt::CaseInsensitive) )
        imgColorTbl = widImage::slColorNames()[6];
    else
        imgColorTbl = widImage::slColorNames()[0]; // tblGrey;

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
    if ( parMethod.isEmpty() )
        parMethod = "bct";
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
    if ( parMethod.isEmpty() )
        parMethod = "bct";
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


    if ( !autofit.isEmpty() )
    {   // Bei automatischem Fit wird die Eingabe von AI ignoriert
        aifile = "";
        tpvfiles.clear();
        //if ( parMethod.isEmpty() ) parMethod = "bct"; - jetzt nur noch Generic
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
        if ( onlyCalc && autofit.isEmpty() )
            flog->write(qPrintable(QString("Only one image calculated.")+EOL));

        flog->write( qPrintable(        "AI-Calculationfile: "+aifile+EOL) );
        flog->write( qPrintable(        "     Parameterfile: "+paramfile+EOL) );
        flog->write( qPrintable(        "         Imagefile: "+imgfile+EOL) );
        flog->write( qPrintable(QString("      Threads used: %1").arg(nthreads)+EOL) );
        flog->write( qPrintable(        "     Automatic fit: "+autofit+EOL) );
        flog->write( qPrintable(        "   Auto fit output: "+af_output+EOL) );
        flog->write( qPrintable(        "    CSV-Outputfile: "+csvfile+EOL) );
        flog->write( qPrintable(        "    TPV-Inputfiles: "+QString::number(tpvfiles.size())+EOL) );
        flog->write( qPrintable(        "   Selected method: "+parMethod+EOL) );
        flog->write( qPrintable(QString("       No FFT Flag: %1").arg(nofft)+EOL) );
        flog->write( qPrintable(QString("             Image:%1%2 Rot=%3deg, Zoom=%4, Color=%5")
                               .arg(imgSwapH?" SwapHor,":"").arg(imgSwapV?" SwapVert,":"")
                               .arg(imgRot*90).arg(imgZoom).arg(imgColorTbl)+EOL) );
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
    std::cerr << " METHOD: " << qPrintable(parMethod) << std::endl;
    std::cerr << " no FFT: " << nofft                 << std::endl;
    std::cerr << qPrintable(QString("  Image:%1%2 Rot=%3deg, Zoom=%4, Color=%5")
                           .arg(imgSwapH?" SwapHor,":"").arg(imgSwapV?" SwapVert,":"")
                           .arg(imgRot*90).arg(imgZoom).arg(imgColorTbl)) << std::endl;

    QDateTime dtStart = QDateTime::currentDateTime();

    SC_CalcCons *calc = new SC_CalcCons;

    slCalcArt = calc->getCalcTypes();
    if ( !paramfile.isEmpty() ) calc->loadParameter( paramfile );

    if ( tpvfiles.size() > 0 )
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
                aifile += QString("|LI%1;%2").arg(sets.value("val_addLinesH",1).toInt())
                           .arg(sets.value("val_addLinesV",1).toInt());
            if ( sets.value("ena_addNoise",false).toBool() )
                aifile += "|NO";
            if ( sets.value("ena_convolute",false).toBool() )
                aifile += "|CO";
            if ( sets.value("ena_calcRphi",false).toBool() )
                aifile += "|RP";

            if ( sets.value("ena_calcFFT",false).toBool() )
            {
                sets.beginGroup("FFT");
                static const QString fftOutFormat[4] = { "OutRe", "OutIm", "OutAbs", "OutSpec" };
                static const QString fftSize[5] = { "32", "64", "128", "256", "512" };
                QString rphiScale = "";
                if ( sets.value("FFTScaleRphi",false).toBool()  ) rphiScale += "Scale ";
                if ( sets.value("FFTclipRphi",false).toBool()   ) rphiScale += "Clip1 ";
                if ( sets.value("FFTclip40Rphi",false).toBool() ) rphiScale += "Clip4 ";
                QString fftScale = "";
                if ( sets.value("FFTScaleOutput",false).toBool()  ) fftScale += "Scale ";
                if ( sets.value("FFTclipOutput",false).toBool()   ) fftScale += "Clip1 ";
                if ( sets.value("FFTclip40Output",false).toBool() ) fftScale += "Clip4 ";
                aifile += QString("@FFT|%1;%2;%3;%4;%5;%6;%7")
                              .arg(sets.value("FFTLinInput",false).toBool()?"InLin":"InLog")
                              .arg(fftSize[sets.value("FFTsizeRphi",2).toInt()])
                              .arg(rphiScale.trimmed())
                              .arg(fftSize[sets.value("FFTsizeOut",2).toInt()])
                              .arg(fftOutFormat[sets.value("FFToutput",0).toInt()])
                              .arg(fftScale.trimmed())
                              .arg(sets.value("FFTSwapOutput",false).toBool()?"OutSwap":"OutNoSwap");
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
                           QString parMethod, QString fnpng );
        if ( flog )
        {
            flog->write("*****************************************" EOL);
            flog->write("* Perform automatic fit                 *" EOL);
            flog->write("*****************************************" EOL);
        }
        automaticFit( calc, imgfile, autofit, nthreads, parMethod, af_output );
    }
    else
    {   // Generate only one image
        void generateSingleImage( SC_CalcCons *calc, QString fnpng, int nthreads, QString parMethod );
        if ( flog )
        {
            flog->write("*****************************************" EOL);
            flog->write("* Generate single image file            *" EOL);
            flog->write("*****************************************" EOL);
        }
        generateSingleImage( calc, imgfile, nthreads, parMethod );
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
    bool fftInputLin,
         fftRphiScale, fftRphiClip, fftRphiClip40,
         fftOutScale,  fftOutClip,  fftOutClip40,
         fftOutputSwap;
    int  fftRphiSize=0, fftOutSize=0;
#define fftOutRe   0
#define fftOutIm   1
#define fftOutAbs  2
#define fftOutSpec 3
    int fftOutFormat;

    QString fn, fnpng;
    bool isDef;
    QString dbg;
    // Merker für die "Training Parameter Variation" (TPV)
    bool useTPV = false;
    double tpvBSx0=-1, tpvBSy0=-1, tpvBSxs=-1, tpvBSys=-1;
    int tpvLinesX=0, tpvLinesY=0;
    bool tpvAddNoise=false, tpvConvolute=false, tpvAddRphi=false;

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
                    /*
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
                    */
                    static QStringList fftOut = { "OutRe", "OutIm", "OutAbs", "OutSpec" };
                    QStringList fft = sl[1].split(";");
                    while ( fft.size() < 7 ) fft << "0";
                    // FFT|InLog;256;Scale;256;OutAbs;Clip1 Clip4;OutSwap
                    //      0     1   2     3   4      5           6
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
                        if ( val.size() != 2 ) continue;
                        tpvLinesX = val[0].toInt();
                        tpvLinesY = val[1].toInt();
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
                }
                std::cerr << "TPV params:";
                if ( tpvBSx0 > 0 ) std::cerr << " BS: " << tpvBSx0 << "/" << tpvBSy0
                              << " (" << tpvBSxs << "x" << tpvBSys << ")";
                if ( tpvLinesX > 0 ) std::cerr << " Lines: " << tpvLinesX << "/" << tpvLinesY;
                if ( tpvAddNoise ) std::cerr << " Noise";
                if ( tpvConvolute ) std::cerr << " Convolute";
                if ( tpvAddRphi ) std::cerr << " Save(r,phi)";
                std::cerr << std::endl;
                continue;
            } // if ( sl[0] == "TPV" )

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
                    calc->updateParamValue( slCalcArt[iCurCalcArt], vv[0], vv[1].toDouble() );
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
        QDir().mkpath( path + category );
        if ( useTPV )
        {
            // category = "scat" in diesem Fall
            QDir().mkpath( path + "rphi" );
            QDir().mkpath( path + "fft"  );
            QDir().mkpath( path + "scat_png" );
            QDir().mkpath( path + "rphi_png" );
            QDir().mkpath( path + "fft_png"  );
        }

        // Berechnen
        calc->prepareCalculation( slCalcArt[iCurCalcArt], true );
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
        calc->doCalculation( nthreads, nullptr );

        double *src = calc->data();
        int sx = calc->maxX() - calc->minX();
        int sy = calc->maxY() - calc->minY();
        int len = sx * sy;
        int bsx0 = calc->currentParamValue("","BeamPosX");   // Position
        int bsy0 = calc->currentParamValue("","BeamPosY");
        int bsdx = 2; // TODO  Größe
        int bsdy = 2;

        if ( useTPV )
        {   // Nachbearbeitung
            QString str = "";
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

            }
            if ( tpvLinesY > 0 )
            {

            }
            if ( flog )
                flog->write(qPrintable(str));
        }

        std::cerr << qPrintable(fn) << " -> " << calc->higResTimerElapsed(SC_CalcCons::htimBoth) << std::endl;
        if ( flog )
            flog->write(qPrintable(QString("Calctime: %1ms").arg(calc->higResTimerElapsed(SC_CalcCons::htimBoth))+EOL));

        saveAIfile( fnpng, src, sx, sy, 0 );    // Berechnetes Bild

        if ( (fftOutSize > 0 || fftRphiSize > 0) && !nofft )
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

                saveAIfile( fnpng, rphi, fftRphiSize, fftRphiSize, 1 ); // r,phi speichern
            }

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
                                                          fftOutputSwap );

            if ( src != nullptr )
                saveAIfile( fnpng, src, fftOutSize, fftOutSize, 2 );
        } // if ( fftImgSize > 0 && !nofft )

        imgCount++;

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
{   // linflg: 0=png, 1=r/phi, 2=fft

    static widImage *img = nullptr;
    if ( imgColorTbl[0] == '*' )
    {   // Spezielle Kennung zum Speichern in Textformat (AIgen)
        static const int prec = 6;
        //static const int digi = 2;
        if ( linflg == 1 )
            fn.replace("scat","rphi");
        else if ( linflg == 2 )
            fn.replace("scat","fft");
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
    }
    //else  Bei TPV auch das farbige Image speichern (macht das Testen einfacher)
    {   // Normale Speicherung als PNG
        if ( img == nullptr ) img = new widImage();
        img->setConfig(imgColorTbl, imgSwapH, imgSwapV, imgRot, imgZoom, linflg!=0);
        //img->setData( calc->minX(), calc->maxX(), calc->minY(), calc->maxY(), src );
        img->setData( 0, sx, 0, sy, src );
        img->saveImage( fn+".png" );    // auch in Farbe speichern möglich
    }

    /*if ( flog )       TODO
        {
            double min, max;
            img->getVarScaling( min, max );
            std::cerr << "  MinVal=" << min << ",  MaxVal=" << max << std::endl;
            flog->write(qPrintable(QString("  MinVal=%1,  MaxVal=%2").arg(min).arg(max)));
        }*/
}


void generateSingleImage( SC_CalcCons *calc, QString fnpng, int nthreads, QString parMethod )
{
    if ( !fnpng.endsWith(".png",Qt::CaseInsensitive) ) fnpng += ".png";
    if ( parMethod.isEmpty() )
        iCurCalcArt = 0;
    else
    {
        iCurCalcArt = 0;
        foreach (QString s, slCalcArt)
        {
            if ( s.startsWith(parMethod,Qt::CaseInsensitive) )
                break;
            iCurCalcArt++;
        }
    }

    // Berechnen
    calc->prepareCalculation( slCalcArt[iCurCalcArt], true );
    calc->doCalculation( nthreads, nullptr );
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
    if ( flog ) flog->write(qPrintable(msg+EOL));
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


double doFitStart( SC_CalcCons *calc, widImage *img, int nthreads, QString parMethod,
                   double inpFitStepSize, double tolerance,
                   int inpFitRepetitions, int inpFitMaxIter,
                   int inpFitBorder, int inpFitBStop/*-1 for mask*/,
                   double &timeForAll, int &loopsForAll, int &imgGenForAll )
{
    if ( fitClass == nullptr )
        fitClass = new SasCalc_SimplexFit2D( calc );

    fitClass->setImageInfo( img->xmin(), img->xmax(), img->ymin(), img->ymax(),
                            img->getFileInfos()->centerX, img->getFileInfos()->centerY,
                            img->dataPtr() );

    calc->prepareCalculation( parMethod, true );

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


/**
 * @brief automaticFit - performs the automatic 2d simplex fit with recipe
 * @param calc      = pointer to calculation class
 * @param imgfile   = input imagefilepath/name if only one image is fitted
 *                    empty if multiple images are fitted (given in recipe)
 * @param autofit   = filepath/name with recipe
 * @param nthreads  = number of threads to be used (0=GPU if available)
 * @param parMethod = calculation method to be used
 * @param fnpng     = output image filepath/name if only one image is fitted
 *                    filepath for all other output file if multiple images are fitted
 */
void automaticFit( SC_CalcCons *calc, QString imgfile, QString autofit, int nthreads,
                   QString parMethod, QString fnpng )
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

    widImage *imgZiel = nullptr;
    int bs_x=0, bs_y=0;

    if ( !imgfile.isEmpty() )
    {
        imgZiel = SC_ReadData::readImageKWSData( addImage, imgfile );
        if ( imgZiel == nullptr )
        {
            if ( flog ) flog->write(qPrintable("unable to open "+imgfile+EOL));
            std::cerr << "Unable to open data file: " << qPrintable(imgfile) << std::endl;
            return;
        }
        SC_ReadData::findBeamCenter( imgZiel, bs_x, bs_y );
        imgZiel->addMetaInfo( "BeamPosX", QString::number(bs_x) );
        imgZiel->addMetaInfo( "BeamPosY", QString::number(bs_y) );
        calc->updateParamValue( parMethod, "BeamPosX", bs_x - imgZiel->myWidth()/2. );
        calc->updateParamValue( parMethod, "BeamPosY", bs_y - imgZiel->myHeight()/2. );
        if ( flog )
            flog->write(qPrintable(QString("AutoBeamstop: in data (%1 / %2), in calc (%3 / %4)").arg(bs_x).arg(bs_y)
                                   .arg(bs_x - imgZiel->myWidth()/2.).arg(bs_y - imgZiel->myHeight()/2.)+EOL));
        std::cerr << "ImageInfos: X=" << imgZiel->xmin() << " .. " << imgZiel->xmax()
                  << ", Y=" << imgZiel->ymin() << " .. " << imgZiel->ymax()
                  << ", BS=" << imgZiel->getFileInfos()->centerX << " / " << imgZiel->getFileInfos()->centerY
                  << std::endl;
    }

    if ( fcsv ) fcsv->write("Variables ; rtol ; MeanDiffPercentage" EOL);

    // Die globale Variable 'flog' wird auch bei der Berechnung des Fits mit verwendet. Beim Aufruf hier ist das
    // eine globale Logdatei, die per Parameter geschaltet werden kann. Diese nehme ich hier lokal in 'globLog'
    // und nutze die 'flog' um für jedes Bild der Serie eine eigene Logdatei zu schreiben, die bleibt dann kleiner.
    QFile *globLog = flog;
    flog = nullptr;

    // Da die Eingaben z.T. in anderer Groß/Kleinschreibung kommt als in der Definition
    // festgelegt, wird hier ein Angleich ausgeführt, damit die folgenden Routinen laufen.
    QStringList types = calc->getCalcTypes();
    foreach (QString s, types)
    {
        if ( s.startsWith(parMethod,Qt::CaseInsensitive) )
        {
            parMethod = s;
            break;
        }
    }

    // Jetzt die Struktur p2f füllen mit allen fitbaren Variablen aus den Definitionne der angegebenen Methode
    QStringList slParams = calc->paramsForMethod( parMethod, false, false, true );
    foreach (QString p, slParams)
    {
        bool cnt;
        double min, max;
        //qDebug() << "FIT: create new" << p;
        _fitLimits *fl = new _fitLimits;
        fl->orgval = fl->fitstart = calc->currentParamValue( parMethod, p );
        fl->fitType  = _fitTypes::fitNone;
        fl->fitvalid = false;
        fl->used = false;
        if ( calc->limitsOfParamValue( parMethod, p, min, max, cnt ) )
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
    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup( "Fit-"+parMethod );
    QStringList slKeys = sets.allKeys();
    //qDebug() << m << slKeys;
    foreach (QString k, slKeys)
    {
        if ( ! calc->isCurrentParameterValid(parMethod,k) ) continue;
        _fitLimits *fl = p2f.value(k,nullptr);
        if ( fl == nullptr ) continue;
        QStringList slVal = sets.value(k,"0:0:0:0").toString().split(":");
        fl->used = slVal[0].toInt() != 0;
        fl->min  = slVal[1].toDouble();
        fl->max  = slVal[2].toDouble();
        fl->fitType = static_cast<_fitTypes>(slVal[3].toInt());
        fl->fitvalid = false;  // fitresult not stored!
        //qDebug() << "REG:" << k << fl->min << fl->max;
    }

    bool firstScanOfFile = true;
    QString dirMask = "*";
    QStringList usedImages; // wird während des Abarbeitens geleert, daher kein Index nötig
    QFile *ftex = nullptr;
    QDateTime startZeit;

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
            if ( line.isEmpty() ) continue;         // Ignore empty lines
            if ( line.startsWith("#") ) continue;   // Comment

            // Keine Kommentare in das LaTeX File schreiben lassen, dann gibt es auch keine
            //  Probleme mit den Umlauten (bzw. UTF-8)
            if ( firstScanOfFile ) slInputLines << line;

            std::cerr << "AutoFit:" << qPrintable(line) << std::endl;
            if ( globLog != nullptr && firstScanOfFile ) globLog->write(qPrintable("Autofit: "+line+EOL));

            if ( firstScanOfFile )
            {   // Einige Einträge werden nur beim ersten Durchlauf interpretiert

                if ( line.startsWith("GlobLog:") )
                {   // "GlobLog: <filename>  Globales Logfile für alle Informationen
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
                    it.value()->orgval = it.value()->fitstart = calc->currentParamValue( parMethod, it.key() );
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
                    qDebug() << "ERROR" << slNames << p2f.keys();
                    finp.close();
                    return;
                }
                continue;
            } // if "Use:..."

            else if ( line.startsWith("Fit:") )
            {   // "Fit: Rep=10; Stp=3.0; Iter=20; Tol=0.0010; Diff<5; Kenn=..."
                QStringList slCmd = line.mid(4).trimmed().split(";");
                for ( int i=0; i<slCmd.size(); i++ )
                    slCmd[i] = slCmd[i].trimmed();
                //qDebug() << slCmd;
                double maxDif=100, stp=3.0, tol=0.001;
                int rep=10;
                int iter=10;
                for ( int i=0; i<slCmd.size(); i++ )
                {
                    if ( slCmd[i].startsWith("Rep=") )
                        rep = slCmd[i].mid(4).toInt();
                    else if ( slCmd[i].startsWith("Stp=") )
                        stp = slCmd[i].mid(4).toDouble();
                    else if ( slCmd[i].startsWith("Iter=") )
                        iter = slCmd[i].mid(5).toInt();
                    else if ( slCmd[i].startsWith("Tol=") )
                        tol = slCmd[i].mid(4).toDouble();
                    else if ( slCmd[i].startsWith("Diff<") )
                        maxDif = slCmd[i].mid(5).toDouble();
                    else if ( slCmd[i].startsWith("Kenn=") )
                    {   // Diese Kennung war zum speichern der Logfiles in der interaktiven Session.
                        //kenn = slCmd[i].mid(5);
                        //kenn = kenn.replace( "@D", QDate::currentDate().toString("yyyyMMdd") );
                    }
                }
#ifndef USEREPETITIONS
                rep = 1;    // ist so am einfachsten...
#endif
                double timeForAll=0;        // Summen über alle Loops
                int    loopsForAll=0;
                int    imgGenForAll=0;

                if ( imgfile.isEmpty() &&
                     usedImages.size() > 0 &&
                     curImgZiel != usedImages.first() )
                {   // Jetzt muss das erste Image geladen werden ...

                    // TODO: jedes Imageformat lesen (ConsFit)
                    // TODO: Zudem die richtigen Metadaten neben dem File lesen

                    imgZiel = SC_ReadData::readImageKWSData( addImage, usedImages.first() );
                    if ( imgZiel == nullptr )
                    {
                        if ( globLog ) globLog->write(qPrintable("unable to open "+usedImages.first()+EOL));
                        if ( flog ) flog->close();
                        flog = globLog;
                        if ( ftex ) ftex->close();
                        return;
                    }
                    curImgZiel = usedImages.first();
                    if ( flog ) flog->close();
                    flog = new QFile( basePath+QFileInfo(curImgZiel).baseName()+"_calc.log");
                    if ( ! flog->open(QIODevice::WriteOnly) )
                    {
                        if ( globLog ) globLog->write(qPrintable("unable to open logfile ("+flog->fileName()+"). This run is without logfile." EOL));
                        flog->deleteLater();
                        flog = nullptr;
                    }
                    else if ( globLog )
                        globLog->write(qPrintable("Use "+flog->fileName()+EOL));
                    SC_ReadData::findBeamCenter( imgZiel, bs_x, bs_y );
                    imgZiel->addMetaInfo( "BeamPosX", QString::number(bs_x) );
                    imgZiel->addMetaInfo( "BeamPosY", QString::number(bs_y) );
                    calc->updateParamValue( parMethod, "BeamPosX", bs_x - imgZiel->myWidth()/2. );
                    calc->updateParamValue( parMethod, "BeamPosY", bs_y - imgZiel->myHeight()/2. );
                    if ( flog )
                        flog->write(qPrintable(QString("AutoBeamstop: in data (%1 / %2), in calc (%3 / %4)").arg(bs_x).arg(bs_y)
                                                  .arg(bs_x - imgZiel->myWidth()/2.).arg(bs_y - imgZiel->myHeight()/2.)+EOL));
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
                              << std::endl;
                } // if usedImages > 0

                for ( int maxloop=0; maxloop<20; maxloop++ ) // keine Endlosschleife...
                {
                    double timeForOne;      // Summen über alle Repetitions/Iterations
                    int    loopsForOne;
                    int    imgGenForOne;
                    double fitMeanChangePercent = doFitStart( calc, imgZiel, nthreads, parMethod,
                                                              stp, tol, rep, iter,
                                                              0, -1/* for mask*/,
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
            calc->prepareCalculation( parMethod, false );
            calc->doCalculation( nthreads, nullptr );
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
        calc->prepareCalculation( parMethod, false );
        calc->doCalculation( nthreads, nullptr );
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
