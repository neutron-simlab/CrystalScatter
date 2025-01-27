#include "widimage.h"
#ifndef CONSOLENPROG
#include "ui_widimage.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QToolTip>
#include <QPainter>
#include <QStyleFactory>

#ifndef NOQWT
#include <QwtLogScaleEngine>
#include <QwtPlotRenderer>
#include <QwtText>
#endif

#endif
#include <QFile>
#include <QDebug>
#include <math.h>
#include <iostream>


#define NOADJUST
// Wenn definiert, werden alle adjustSize() deaktiviert. Dadurch sollte die Größe
// des Fensters stabil bleiben. Durch die ScrollArea um das Image-Label ist die
// MinimumSize des Labels nicht mehr relevant für das Fenster.


//#define UseAutoNoBorder
// Wenn definiert, gibt es eine zusätzliche Skalierungsauswahl:
// 'Automatische Skalierung ohne Berücksichtigung der Ränder'
// Ich habe einmal gesehen, dass an den Rändern unklare Artefakte waren, durch die
// die automatische Skalierung unbrauchbar wurde.


#define CountNanValues
// Wenn definiert, wird bei setData() die Anzahl der NaN Werte gezählt und mit qDebug() ausgegeben



#define SETT_APP   "JCNS-1-SasCrystal"  // auch in sc_maingui.h
#define SETT_IMG   "ImageScaling"

#define HISTO_SY  60


QString widImage::lastSavedFilePath = "";
int widImage::lastColorTbl = 1;

#ifndef CONSOLENPROG
static float zoom2pix[] = {  1,  2,  4, 0.5, 0.25, 0      };
//                          No, *2, *4,  /2,   /4, 640*640
#endif

QStringList widImage::_slColornames = { /* 0*/"Greyscale",     /* 1*/"Glowing colors",     /* 2*/"Earth colors",         /* 3*/"Temperature",
                                        /* 4*/"Scatter-Color", /* 5*/"Scatter-Colorwheel", /* 6*/"Scatter-Temperature",  /* 7*/"Scatter-Geo",
                                        /* 8*/"Scatter-Desy",  /* 9*/"Scatter-Grey",       /*10*/"Scatter-Inverse grey", /*11*/"Scatter-Jet",
                                        /*12*/"Residuen Glow", /*13*/"Residuen Temp",      /*14*/"Mathematica",          /*15*/"Mathematica-Inv",
                                        /*16*/"SasView" };


widImage::widImage() :
    noGUI(true)
{
#ifndef CONSOLENPROG
    getConfig();
#endif
}


#ifndef CONSOLENPROG
widImage::widImage( QString tit, QString dp, QWidget *parent ) :
    QWidget(parent),
    ui(new Ui::widImage),
    noGUI(false)
{
    useImageFileOnly = false;
    histoThread = nullptr;
    histogramValid = false;
    useResiduenColors = false;
    enaExtract = false;
    enaNoFit   = false;
    firstView  = true;
    extractDrawByMouse = false;
#ifndef NOQWT
    plot1d  = nullptr;
    curve1d = nullptr;
    grid1d  = nullptr;
#endif
    if ( lastSavedFilePath.isEmpty() && !dp.isEmpty() ) lastSavedFilePath = dp;
    getConfig();
    ui->setupUi(this);
    ui->scrollArea1d->hide();
#ifdef Q_OS_LINUX
    // Sonst haben die GroupBox keinen Rahmen...
    QStyle *winstyle = QStyleFactory::create("Windows");
    ui->grpImage->setStyle(winstyle);
    ui->grpScaling->setStyle(winstyle);
#endif
    setWindowTitle( tit );
    QStringList sl;
    sl << "Key" << "Value";
    ui->tblMetaData->setColumnCount(sl.size());
    ui->tblMetaData->setHorizontalHeaderLabels(sl);
    ui->tblMetaData->hide();
    ui->lblHistogram->hide();
    ui->cbsColorTbl->addItems( _slColornames );
    ui->cbsColorTbl->setCurrentIndex(lastColorTbl);
    on_cbsColorTbl_textActivated( ui->cbsColorTbl->currentText() );
    ui->lblFilename->setText( "not saved yet" );
    ui->lblFilename->setEnabled(false);
    ui->cbsZoom->setCurrentIndex( zoomFactor );
    ui->cbsRotation->setCurrentIndex( rotIndex );
    ui->togMirrorH->setChecked( swapHor );
    ui->togMirrorV->setChecked( swapVert );
    ui->lblHistogram->setMaximumHeight( 10 + HISTO_SY );
    ui->lblHistogram->setMinimumHeight( 10 + HISTO_SY );
#ifndef UseAutoNoBorder
    ui->radScaleAutoNoBorder->hide();
#endif
    setMouseTracking(true);
    if ( /*tit.startsWith("(r,phi)") ||*/ tit.startsWith("iFFT") )
        ui->radShowLin->setChecked( true );
    bRekursiv = false;
}

void widImage::getConfig()
{
    QSettings set(SETT_APP,SETT_IMG);
    swapHor  = set.value("SwapHor",false).toBool();
    swapVert = set.value("SwapVert",false).toBool();
    rotIndex = set.value("RotIndex",0).toInt();
    zoomFactor = set.value("ZoomFactor",0).toInt();
}

widImage::~widImage()
{
    if ( noGUI ) return; // Destructor
    delete ui;
}

void widImage::closeEvent(QCloseEvent *event)
{
    QSettings set(SETT_APP,SETT_IMG);
    set.setValue("SwapHor",swapHor);
    set.setValue("SwapVert",swapVert);
    set.setValue("RotIndex",rotIndex);
    set.setValue("ZoomFactor",zoomFactor);
    emit destroyed(this);
    event->accept();
}

void widImage::saveNoFitRects()
{
    if ( fileInfos.filePath.isEmpty() ) return;
    if ( !QFile::exists(fileInfos.filePath) ) return;
    QFile fout(fileInfos.filePath+".nofit");
    if ( fout.open(QIODevice::WriteOnly) )
    {
        for ( int i=0; i<4; i++ )
            fout.write( qPrintable(QString("%1 %2 %3 %4\n")
                                      .arg(noFitRect[i].left()).arg(noFitRect[i].top())
                                      .arg(noFitRect[i].right()).arg(noFitRect[i].bottom())) );
        fout.close();
    }
}
#endif // !CONSOLENPROG

QString widImage::getMetaTitle()
{
#ifndef CONSOLENPROG
    return ui->grpImage->title().remove("Image ");
#else
    return "?";
#endif
}

void widImage::addMetaInfo( QString key, QString val )
{
#ifndef CONSOLENPROG
    if ( noGUI ) return; // addMetaInfo
    if ( key.isEmpty() )
    {   // Alles löschen
        ui->tblMetaData->clearContents();
        ui->tblMetaData->setRowCount(0);
        //fileInfos.isValid = 0; wird bei setData() gemacht
        return;
    }
    if ( key == "@" )
    {   // Abschluß, Sortieren
        ui->lblFilename->setText( "not saved yet" );
        ui->lblFilename->setEnabled(false);
        ui->tblMetaData->sortItems(0);
        ui->tblMetaData->resizeColumnsToContents();
        return;
    }
    if ( val == "?2" )
        return;     // ReadOnly Parameter, werden im Laufe der Berechnungen bestimmt
    if ( key == "_Calculation_" )
        ui->grpImage->setTitle( "Image "+val );
    else if ( key == "From Image" )
        ui->grpImage->setTitle( val );
    //--- Merker für den genutzten Mittelpunkt: das Bild oder die Beamstop-Werte
    else if ( key == "CenterMidpoint" )
        return;
    else if ( key == "CenterBeam" )
    {
        useCenterBeamPos = val[0] == 'T';
        //qDebug() << key << val;
        //Debug: "CenterBeam" "False"
        //Debug: "CenterBeam" "True"
    }
    //--- Merker für den genutzten QMax-Wert, wird nur beim xy2q() verwendet
    /*else if ( key.contains("qmax",Qt::CaseInsensitive) )
    {
        //qmax = val.toDouble();
        //qDebug() << key << val;
        //Debug: "EditQmaxData" "False"  --> wenn True, dann CalcQmax nutzen
        //Debug: "CalcQmax" "1.665"
        //Debug: "EditQmax" "2"
        //Debug: "EditQmaxPreset" "True" --> wenn True, dann EditQmax nutzen
    }*/
    else if ( key == "EditQmax" && !data1D )
    {
        qmaxset = false;
        editQmax = val.toDouble();
    }
    else if ( key == "CalcQmax" && !data1D )
    {
        qmaxset = false;
        calcQmax = val.toDouble();
    }
    else if ( key == "EditQmaxPreset" && !data1D )
    {
        qmaxedit = val[0] == 'T';
    }
    else if ( key == "EditQmaxData" && !data1D )
    {
        qmaxcalc = val[0] == 'T';
    }
    //---
    else if ( key == "From File" )
    {
        fileInfos.filePath = val;
        QString fn = val;
        QFile fin;
        if ( fn.endsWith(".csv",Qt::CaseInsensitive) )
            fin.setFileName( fn.left(fn.length()-4)+".dat.nofit" );
        else
            fin.setFileName( fn+".nofit" );
        if ( fin.open(QIODevice::ReadOnly) )
        {
            for ( int i=0; i<4; i++ )
            {
                QStringList sl = QString(fin.readLine()).simplified().split(" ");
                while ( sl.size() < 4 ) sl << "-1";
                noFitRect[i].setLeft(   sl[0].toInt() );
                noFitRect[i].setTop(    sl[1].toInt() );
                noFitRect[i].setRight(  sl[2].toInt() );
                noFitRect[i].setBottom( sl[3].toInt() );
            }
            fin.close();
        }
        if ( fn.size() > 40 )
            fn = fn.left(3) + "..." + fn.right(34);
        ui->grpImage->setTitle( fn );
    }
    // Bestimmte Werte werden beim 1D bzw. 2D ausgeblendet
    if ( data1D )
    {
        if ( key == "GridPoints"    ) return;
        if ( key == "RadioButtonQ1" ) return;
        if ( key == "RadioButtonQ2" ) return;
        if ( key == "RadioButtonQ4" ) return;
        if ( key == "ExpandImage"   ) return;
    }
    else
    {
        if ( key == "EditQmin"   ) return;
        if ( key == "EditQsteps" ) return;
    }
    QList<QTableWidgetItem*> items = ui->tblMetaData->findItems( key, Qt::MatchFixedString );
    if ( items.size() == 0 )
    {
        int row = ui->tblMetaData->rowCount();
        ui->tblMetaData->setRowCount( row + 1 );
        ui->tblMetaData->setItem( row, 0, new QTableWidgetItem(key) );
        ui->tblMetaData->setItem( row, 1, new QTableWidgetItem(val) );
    }
    else
        ui->tblMetaData->item(items[0]->row(),1)->setText(val);
#endif
    // Spezielle Werte noch lokal speichern
    QString val2 = val.left(val.indexOf(' ')).trimmed();
    if ( val2.isEmpty() ) val2 = val.trimmed();
    if ( key == "BeamPosX" )
    {
        fileInfos.centerX = val2.toDouble();
        fileInfos.isValid |= fivalidCenX;
    }
    else if ( key == "BeamPosY" )
    {
        fileInfos.centerY = val2.toDouble();
        fileInfos.isValid |= fivalidCenY;
    }
    else if ( key.contains("wave",Qt::CaseInsensitive) )
    {
        fileInfos.wavelen = val2.toDouble();
        fileInfos.isValid |= fivalidWave;
    }
    else if ( key == "EditDet" || key == "SampleDist" )
    {
        fileInfos.distance = val2.toDouble();
        fileInfos.isValid |= fivalidDist;
    }
    else if ( key == "EditPixelX" || key == "Pixel_X" )
    {
        fileInfos.pixWidthX = val2.toDouble();
        fileInfos.isValid |= fivalidWidX;
    }
    else if ( key == "EditPixelY" || key == "Pixel_Y" )
    {
        fileInfos.pixWidthY = val2.toDouble();
        fileInfos.isValid |= fivalidWidY;
    }
    else if ( key.startsWith("NoFitRect_") )
    {
        QStringList sl = val.simplified().split(" ");
        while ( sl.size() < 4 ) sl << "-1";
        int i = key.right(1).toInt();
        if ( i >= 0 && i < 4 )
        {
            noFitRect[i].setLeft(   sl[0].toInt() );
            noFitRect[i].setTop(    sl[1].toInt() );
            noFitRect[i].setRight(  sl[2].toInt() );
            noFitRect[i].setBottom( sl[3].toInt() );
        }
    }
    else if ( key == "Residuenplot" )
    {
#ifndef CONSOLENPROG
        ui->cbsColorTbl->setCurrentIndex( ui->cbsColorTbl->findText(_slColornames[12]) ); // tblResGlo
#endif
        on_cbsColorTbl_textActivated(_slColornames[12]); // tblResGlo
#ifndef CONSOLENPROG
        on_butUpdate_clicked();
#else
        imgNormal = generateImage();
#endif
    }
#ifndef CONSOLENPROG
    if ( fileInfos.isValid == fivalid_All )
    {
        fileInfos.isValid |= fivalidDone;
        /*
        int minx, miny;
        double q, qx, qy, qz, qmin;
        qmin = 1e6;
        //for ( int x=0; x<fileInfos.pixelX; x++ )
        //    for ( int y=0; y<fileInfos.pixelY; y++ )
        for ( int x=minX; x<maxX; x++ )
            for ( int y=minY; y<maxY; y++ )
            {
                q = xy2q( x, y, qx, qy, qz );
                if ( q < qmin )
                {
                    qmin = q;
                    minx = x;
                    miny = y;
                }
            }
        qDebug() << "Calc minQ:" << minx << miny << qmin << "pixel" << fileInfos.pixelX << fileInfos.pixelY;
        qDebug() << "          pixWid" << fileInfos.pixWidthX << fileInfos.pixWidthY
                 << "center" << fileInfos.centerX << fileInfos.centerY
                 << "dist" << fileInfos.distance << "wave" << fileInfos.wavelen;
        qDebug() << "          min,max|x,y" << minX << maxX << minY << maxY;
        */
    }
#endif
}

QString widImage::metaInfo( QString key, QString def )
{
#ifndef CONSOLENPROG
    if ( noGUI ) return def; // metaInfo
    QList<QTableWidgetItem*> items = ui->tblMetaData->findItems( key, Qt::MatchFixedString );
    if ( items.size() == 0 ) return def;
    QString tmp = ui->tblMetaData->item(items[0]->row(),1)->text();
    tmp.replace(" nm","");
    tmp.replace(" mm","");
    tmp.replace(" m","");
    tmp.replace(" deg","");
    return tmp;
#else
    Q_UNUSED(key)
    return def;
#endif
}


#if !defined(CONSOLENPROG) && !defined(NOQWT)
// Die Grenzwerte der X-Achse muss hier schon bekannt sein.
int widImage::setData1D(int x1, double q0, double q1, double *d)
{
    data1D = true;

    switchGuiFor1d();
    // Preset internal structure
    fileInfos.isValid = fivalidPixX | fivalidPixY; //  | fivalidCenX | fivalidCenY
    fileInfos.filePath = "<calc1d>";
    fileInfos.pixelX = x1;
    fileInfos.pixelY = 1;
    fileInfos.centerX = x1 / 2;
    fileInfos.centerY = 0;

    qmax = q1;
    qmaxset = true; // Verhindert ein ändern bei den Metadaten
    qmin_1d = q0;
    minX = 0;
    maxX = x1;
    minY = 0;
    maxY = 1;
    //data.clear();
    int len = x1;
    //qDebug() << "setData1D" << qmin_1d << qmax << "(1d)";
#ifdef CountNanValues
    int cntNan=0;
#endif
    data.resize( len );    // TODO brauche ich die hier?
    points1d.resize(len);
    double qstp = (qmax - qmin_1d) / (double)len;
    for ( int i=0; i<len; i++ )
    {
        QPointF pnt;
        pnt.setX( i*qstp );
        if ( isnan(*d) )
        {
#ifdef CountNanValues
            cntNan++;
#endif
            d++;
            data[i] = 0.0;
            pnt.setY(0.0);
        }
        else
        {
            pnt.setY(*d);
            data[i] = *(d++);
        }
        points1d[i] = pnt;
    }
    loadImageFile("");  // Reset disabled gui elements

#ifdef CountNanValues
    if ( cntNan > 0 ) qDebug() << "widImage::setData  len=" << len << ", NaN=" << cntNan;
#endif

    histogramValid = false;
    if ( !noGUI )   // setData, ena cbs
    {
        ui->butUpdate->setEnabled(true);
        ui->cbsZoom->setEnabled(true);
    }
    on_butUpdate_clicked(); // Redraw Image
#ifdef CountNanValues
    if ( noGUI ) return cntNan;    // setData
#else
    if ( noGUI ) return -1;    // setData
#endif
//#ifndef NOADJUST
// Hier noch ein adjustSize() zulassen, sonst ist das Fenster zu groß
#ifdef IMG_ONE_WINDOW
    parentWidget()->adjustSize();
#else
    if ( firstView ) adjustSize();
    firstView = false;
#endif
    //#endif
    raise(); // show Window on top of other windows
#ifdef CountNanValues
    return cntNan;
#else
    return -1;
#endif
}
#else // #if !defined(CONSOLENPROG) && !defined(NOQWT)
int widImage::setData1D(int x1, double q0, double q1, double *d) { Q_UNUSED(x1) Q_UNUSED(q0) Q_UNUSED(q1) Q_UNUSED(d) return -1; }
#endif


/**
 * @brief widImage::setData
 * @param x0 - X min (1D=0)
 * @param x1 - X max (1D=num points)
 * @param y0 - Y min (1D=0)
 * @param y1 - Y max (1D=1)
 * @param d  - data values
 */
int widImage::setData( int x0, int x1, int y0, int y1, double *d )
{
    //qDebug() << "widImage::setData(" << x0 << x1 << y0 << y1 << d << " )";
    if ( d == nullptr )
    {
#ifndef CONSOLENPROG
        QMessageBox::critical( this, "No data", "No image data given", QMessageBox::Ok );
#endif
        return -2;
    }
#ifndef CONSOLENPROG
    setEnabled(true);   // Komplettes Fenster wird gesperrt wenn falsche Daten da waren,
                        // zur Sicherheit hier wieder freigeben.
#endif

    if ( y1-y0 == 1 )
        return -1;
#ifndef NOQWT
    data1D = false;
#endif
#ifndef CONSOLENPROG
    switchGuiFor1d();
#endif

    // Preset internal structure
    fileInfos.isValid = fivalidPixX | fivalidPixY; //  | fivalidCenX | fivalidCenY
    fileInfos.filePath = "<calc>";
    fileInfos.pixelX = x1 - x0;
    fileInfos.pixelY = y1 - y0;
    fileInfos.centerX = (x1 + x0) / 2;
    fileInfos.centerY = (y1 + y0) / 2;
    //qDebug() << "                  " << fileInfos.pixelX << fileInfos.pixelY
    //         << fileInfos.centerX << fileInfos.centerY;

    minX = x0;
    maxX = x1;
    minY = y0;
    maxY = y1;
    data.clear();
    int len = (x1-x0)*(y1-y0);
    //qDebug() << len << minX << maxX << minY << maxY;
#ifdef CountNanValues
    int cntNan=0;
#endif
    data.reserve( len );
    for ( int i=0; i<len; i++ )
    {
        if ( isnan(*d) )
        {
#ifdef CountNanValues
            cntNan++;
#endif
            d++;
            data.append(0.0);
        }
        else
            data.append(*(d++));
    }
#ifdef CONSOLENPROG
    imgNoGUI = generateImage(); // Zum Abspeichern wichtig...
#else
    loadImageFile("");  // Reset disabled gui elements

#ifdef CountNanValues
    if ( cntNan > 0 ) qDebug() << "widImage::setData  len=" << len << ", NaN=" << cntNan;
#endif

    histogramValid = false;
    if ( !noGUI )   // setData, ena cbs
    {
        ui->butUpdate->setEnabled(true);
        ui->cbsZoom->setEnabled(true);
    }
    on_butUpdate_clicked(); // Redraw Image
#ifdef CountNanValues
    if ( noGUI ) return cntNan;    // setData
#else
    if ( noGUI ) return -1;    // setData
#endif
//#ifndef NOADJUST
// Hier noch ein adjustSize() zulassen, sonst ist das Fenster zu groß
#ifdef IMG_ONE_WINDOW
    parentWidget()->adjustSize();
#else
    if ( firstView ) adjustSize();
    firstView = false;
#endif
//#endif
    raise(); // show Window on top of other windows
#endif
#ifdef CountNanValues
    return cntNan;
#else
    return -1;
#endif
}


#ifndef CONSOLENPROG
/**
 * @brief widImage::loadImageFile
 * @param fn - imagefilename to load (if empty reset the disabled gui elements)
 * @return true if successful, false if failed
 * Load the given imagefile, disables all image manipulations in the gui.
 */
bool widImage::loadImageFile( QString fn )
{
    if ( fn.isEmpty() )
    {
        useImageFileOnly = false;
        ui->grpScaling->setEnabled(true);
        return true;
    }
    bool rv;
    if ( noGUI )
        rv = imgNoGUI.load(fn);
    else
        rv = imgNormal.load(fn);
    if ( !rv )
    {
        qDebug() << "widImage::loadImageFile - error loading image" << fn;
        useImageFileOnly = false;
        ui->grpScaling->setEnabled(true);
        return false;
    }

    useImageFileOnly = true;
#ifndef NOQWT
    data1D = false;
#endif
    ui->grpScaling->setEnabled(false);
    on_butUpdate_clicked();
    return true;
}
#endif


void widImage::saveImage( QString fn )
{
    if ( noGUI )
        imgNoGUI.save(fn);
    else
        imgNormal.save(fn);
}

bool widImage::saveImageColor(QString fn, QString coltbl, QString &dbg)
{
#ifndef CONSOLENPROG
    if ( useImageFileOnly ) return false;
    QString savColTbl = ui->cbsColorTbl->currentText();
    //qDebug() << "saveImageColor" << fn << coltbl << "Old:" << savColTbl;
#endif
    on_cbsColorTbl_textActivated(coltbl);
    QImage tmp = generateImage();
    dbg = QString("FN=%1, COL=%2, IMG=%3*%4").arg(fn,coltbl).arg(tmp.width()).arg(tmp.height());
    bool rv = tmp.save(fn);
#ifndef CONSOLENPROG
    on_cbsColorTbl_textActivated(savColTbl);
#endif
    return rv;
}


void widImage::saveImageGray( QString fn )
{
    QImage imgGray;
    if ( noGUI )
        imgGray = imgNoGUI.convertToFormat( QImage::Format_Grayscale8  );
    else
        imgGray = imgNormal.convertToFormat( QImage::Format_Grayscale8  );
    if ( ! imgGray.save( fn ) )
    {
        qDebug() << "*** Image Write Error ***";
        qDebug() << "*** imgGray:" << imgGray.size() << ", imgNoGUI:" << imgNoGUI.size()
                 << ", imgNormal:" << imgNormal.size() << ", Flag:" << noGUI;
    }
}

void widImage::saveImageBinary( QString fn )
{
    if ( !noGUI ) return;
    QImage imgGray = imgNoGUI.convertToFormat( QImage::Format_Grayscale8  );
    QFile f(fn);
    if ( ! f.open(QIODevice::WriteOnly) ) return;
    for ( int y=0; y<imgGray.height(); y++ )
    {
        const uchar *bits = imgGray.constScanLine(y);
        for ( int x=0; x<imgGray.bytesPerLine(); x++ )
            f.write( qPrintable(QString(" %1").arg(*(bits+x))) );
        f.write("\n");
    }
    f.close();
}

void widImage::saveImage( QString fn, QSize siz )
{
    if ( noGUI ) return;
    if ( siz.isNull() )
        imgNormal.save(fn);
    else
        imgNormal.scaled(siz).save(fn);
}

void widImage::saveImageGray( QString fn, QSize siz )
{
    QImage imgGray;
    if ( siz.isNull() )
    {
        if ( noGUI )
            imgGray = imgNoGUI.convertToFormat( QImage::Format_Grayscale8  );
        else
            imgGray = imgNormal.convertToFormat( QImage::Format_Grayscale8  );
    }
    else
    {
        if ( noGUI )
            imgGray = imgNoGUI.scaled(siz).convertToFormat( QImage::Format_Grayscale8  );
        else
            imgGray = imgNormal.scaled(siz).convertToFormat( QImage::Format_Grayscale8  );
    }
    imgGray.convertToFormat(QImage::Format_RGB888).save( fn );
    // Diese Funktion wird nur verwendet, wenn Messdaten für TensorFlow exportiert
    // werden sollen. Diese werden später nicht als .png eingelesen sondern werden
    // in ein numpy-Array umgewandelt. Dazu muss es aber als RGB gespeichert sein.
    // Das Format_Grayscale8 ist daher an der Stelle nicht richtig, wandelt das Bild
    // aber schnell in ein Graubild um....
}

void widImage::saveImageBinary( QString fn, QSize siz )
{
    if ( noGUI ) return;
    QImage imgGray = imgNormal.convertToFormat( QImage::Format_Grayscale8  );
    imgGray = imgGray.scaled(siz);
    QFile f(fn);
    if ( ! f.open(QIODevice::WriteOnly) ) return;
    for ( int y=0; y<imgGray.height(); y++ )
    {
        const uchar *bits = imgGray.constScanLine(y);
        for ( int x=0; x<imgGray.bytesPerLine(); x++ )
            f.write( qPrintable(QString(" %1").arg(*(bits+x))) );
        f.write("\n");
    }
    f.close();
}


void widImage::on_cbsColorTbl_textActivated(const QString &arg1)
{
    // Glühfarben von Stahl
    static int rGlow[10] = { 0, 102, 129, 153, 204, 254, 255, 255, 255, 255 };
    static int gGlow[10] = { 0,  50,   0,   1,   0,   0, 102, 204, 255, 255 };
    static int bGlow[10] = { 0,   0,   1,   0,   1,   0,   0,   0,   0, 203 };
    // Erdfarben
    static int rEarth[10] = { 0,  85,   0,  32,   0,   0,   0, 124, 255, 255 };
    static int gEarth[10] = { 0,  26,   0, 178, 139, 100, 139, 252, 255, 255 };
    static int bEarth[10] = { 0, 139, 128, 170,  69,   0,   0,   0,   0, 255 };
    // GR: temperature.png (to avoid too long initializers split them into two arrays)
    static int r0[512] = {255,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  2,  3,  4,  5,  6,  7,  8,  8,  9, 10, 11, 11, 12, 13, 14, 14, 15, 16, 17, 18, 19, 20, 21, 21, 22, 23, 24, 24, 25, 26, 27, 27, 28, 29, 30, 30, 31, 32, 33, 34, 35, 36, 37, 37, 38, 39, 40, 40, 41, 42, 43, 43, 44, 45, 46, 47, 48, 49, 50, 50, 51, 52, 53, 53, 54, 55, 56, 56, 57, 58, 59, 59, 60, 61, 62, 63, 64, 65, 66, 66, 67, 68, 69, 69, 70, 71, 72, 72, 73, 74, 75, 76, 77, 78, 78, 79, 80, 81, 81, 82, 83, 84, 84, 85, 86, 87, 88, 89, 90, 91, 91, 92, 93, 94, 94, 95, 96, 97, 97, 98, 99,100,100,101,102,103,104,105,106,107,107,108,109,110,110,111,112,113,113,114,115,116,117,118,119,120,120,121,122,123};
    static int g0[512] = {255,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  2,  2,  3,  4,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 34, 35, 36, 37, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 66, 67, 68, 69, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 98, 99,100,101,102,104,105,106,107,108,109,110,111,112,113,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,133,134,135,136,137,138,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,165,166,167,168,169,170,171,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,231,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,253,254,254,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255};
    static int b0[512] = {255,128,129,130,131,133,134,135,136,137,138,139,140,142,143,144,145,146,147,148,149,151,152,153,154,155,156,157,158,160,161,162,163,164,165,166,167,169,170,171,172,174,175,176,177,178,179,180,181,183,184,185,186,187,188,189,190,192,193,194,195,196,197,198,199,201,202,203,204,205,206,207,208,209,211,212,213,214,215,216,217,218,220,221,222,223,224,225,226,227,229,230,231,232,234,235,236,237,238,239,240,241,243,244,245,246,247,248,249,250,252,253,254,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,254,254,253,252,251,251,250,249,248,248,247,246,245,244,243,242,241,241,240,239,238,238,237,236,235,235,234,233,232,231,230,229,228,228,227,226,225,225,224,223,222,222,221,220,219,219,218,217,216,215,214,213,212,212,211,210,209,209,208,207,206,206,205,204,203,202,201,200,199,199,198,197,196,196,195,194,193,193,192,191,190,190,189,188,187,186,185,184,183,183,182,181,180,180,179,178,177,177,176,175,174,173,172,171,170,170,169,168,168,167,166,165,165,164,163,162,161,160,159,158,158,157,156,155,155,154,153,152,152,151,150,149,149,148,147,146,145,144,143,142,142,141,140,139,139,138,137,136,136,135,134,133,132,131,130,129,129,128,127,126,126,125,124,123};
    static int r1[512] = {123,124,125,126,126,127,128,129,129,130,131,132,133,134,135,136,136,137,138,139,139,140,141,142,142,143,144,145,146,147,148,149,149,150,151,152,152,153,154,155,155,156,157,158,158,159,160,161,162,163,164,165,165,166,167,168,168,169,170,170,171,172,173,174,175,176,177,177,178,179,180,180,181,182,183,183,184,185,186,187,188,189,190,190,191,192,193,193,194,195,196,196,197,198,199,199,200,201,202,203,204,205,206,206,207,208,209,209,210,211,212,212,213,214,215,216,217,218,219,219,220,221,222,222,223,224,225,225,226,227,228,228,229,230,231,232,233,234,235,235,236,237,238,238,239,240,241,241,242,243,244,245,246,247,248,248,249,250,251,251,252,253,254,254,254,254,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,253,252,250,249,248,247,246,245,244,243,241,240,239,238,237,236,235,234,232,231,230,229,228,227,226,225,223,222,221,220,218,217,216,215,214,213,212,211,209,208,207,206,205,204,203,202,201,199,198,197,196,195,194,193,192,190,189,188,187,186,185,184,183,181,180,179,178,177,176,175,174,172,171,170,169,167,166,165,164,163,162,161,160,158,157,156,155,154,153,152,151,149,148,147,146,145,144,143,142,140,139,138,137,136,135,134,133,131,130,129,128,128};
    static int g1[512] = {255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,253,252,251,250,249,248,247,247,246,245,244,243,242,241,240,239,238,237,236,236,235,234,233,232,231,230,230,229,228,227,226,225,224,223,222,221,220,219,219,218,217,216,215,214,213,212,211,210,209,208,208,207,206,205,204,203,202,201,200,199,198,197,196,195,194,193,193,192,191,190,189,188,187,186,185,184,183,182,182,181,180,179,178,177,176,175,174,173,172,171,171,170,169,168,167,166,165,164,163,162,161,160,159,158,157,156,156,155,154,153,152,151,150,149,148,147,146,145,145,144,143,142,141,140,139,138,137,136,135,134,134,133,132,131,130,129,128,127,126,125,124,123,122,121,120,120,119,118,117,116,115,114,113,112,111,110,109,109,108,107,106,105,104,103,102,101,100, 99, 98, 97, 96, 95, 94, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
    static int b1[512] = {123,122,121,120,120,119,118,117,116,115,114,113,113,112,111,110,110,109,108,107,107,106,105,104,103,102,101,100,100, 99, 98, 97, 97, 96, 95, 94, 94, 93, 92, 91, 91, 90, 89, 88, 87, 86, 85, 84, 84, 83, 82, 81, 81, 80, 79, 78, 78, 77, 76, 75, 74, 73, 72, 72, 71, 70, 69, 69, 68, 67, 66, 66, 65, 64, 63, 62, 61, 60, 59, 59, 58, 57, 56, 56, 55, 54, 53, 53, 52, 51, 50, 50, 49, 48, 47, 46, 45, 44, 43, 43, 42, 41, 40, 40, 39, 38, 37, 37, 36, 35, 34, 33, 32, 31, 30, 30, 29, 28, 27, 27, 26, 25, 24, 24, 23, 22, 21, 21, 20, 19, 18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10,  9,  8,  8,  7,  6,  5,  4,  3,  2,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
    // Aus den Mathematica-Bildern aus dem Veröffentlichung
    static int rMat[256] = { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 254, 254, 254, 254, 253, 253,253, 253, 252, 252, 252, 252, 251, 251, 251, 251, 250, 250, 250, 250, 249, 249, 249, 249, 248, 248, 248, 247, 247, 247, 247, 247, 247, 246, 246, 246, 246, 245, 245, 245, 245, 244, 244, 244, 243, 243, 243, 243, 242, 242, 242, 242, 241, 241, 241, 241, 240, 240, 240, 240, 239, 239, 239, 239, 239, 238, 238, 238, 238, 237, 237, 237, 237, 236, 236, 236, 236, 235, 235, 235, 235, 234, 234, 234, 234, 233, 233, 233, 233, 232, 232, 232, 232, 231, 231, 231, 231, 231, 230, 228, 228, 226, 224, 222, 222, 220, 217,215, 213, 213, 211, 209, 207, 207, 205, 203, 201, 201, 199, 197, 195, 195, 193, 191, 189, 187, 187, 185, 183, 181, 181, 179, 177, 175, 175, 173, 171, 169, 167, 167, 165, 163, 161, 161, 159, 157, 155, 155, 153, 151, 149, 149, 147, 145, 143, 141, 141, 139, 137, 135, 135, 133, 131, 129, 129, 127, 125, 123, 123, 121, 119, 117, 115, 115, 112, 110, 108, 108, 106, 104, 101, 101,  99,  97,  95,  95,  93,  91,  89,  87,  87,  85,  83,  81,  81,  79,  77,  75,  75,  73,  71,  69,  67,  67,  66,  64,  63,  63, 61,  59,  58,  58,  56,  55,  53,  53,  51,  50,  48,  46,  46,  45,  43,  42,  42,  40};
    static int gMat[256] = { 242, 241, 241, 240, 239, 238, 238, 237, 236, 235, 233, 233, 232, 231, 230, 230, 229, 228, 227, 227, 226, 225, 224, 224, 223, 222, 221, 221, 221, 220, 219, 218, 218, 217, 217, 216,216, 215, 214, 213, 213, 213, 212, 211, 210, 210, 209, 209, 208, 208, 207, 206, 205, 205, 205, 204, 203, 202, 202, 201, 201, 200, 200, 199, 198, 197, 197, 197, 196, 195, 195, 194, 193, 193, 192, 192, 191, 190, 189, 189, 189, 188, 187, 187, 186, 185, 185, 185, 184, 183, 182, 181, 181, 181, 180, 179, 179, 178, 177, 177, 177, 176, 175, 174, 173, 173, 173, 172, 171, 171, 170, 169, 169, 169, 168, 167, 166, 166, 165, 165, 164, 163, 163, 162, 161, 161, 161, 160, 159, 158, 158, 158, 157, 156, 156, 155, 154,153, 153, 153, 152, 151, 150, 150, 150, 149, 148, 148, 148, 147, 146, 146, 146, 145, 144, 143, 143, 143, 142, 141, 141, 141, 140, 139, 139, 139, 138, 137, 136, 136, 136, 135, 134, 134, 134, 133, 132, 132, 131, 131, 130, 130, 129, 129, 128, 127, 127, 127, 126, 125, 125, 124, 124, 123, 123, 122, 122, 121, 121, 120, 119, 119, 118, 118, 117, 116, 115, 115, 115, 114, 113, 113, 112, 112, 111, 111, 110, 110, 109, 108, 108, 108, 107, 106, 106, 105, 105, 104, 104, 103, 103, 102, 101, 101, 100,  99,  98,  98, 97,  97,  96,  96,  95,  94,  93,  93,  92,  91,  90,  89,  89,  88,  88,  87,  87,  86};
    static int bMat[256] = { 190, 188, 188, 186, 184, 182, 182, 180, 178, 176, 174, 174, 172, 170, 168, 168, 167, 165, 163, 163, 161, 159, 157, 157, 155, 153, 152, 151, 151, 150, 149, 148, 148, 147, 146, 145,145, 144, 143, 142, 141, 141, 140, 139, 138, 138, 137, 136, 135, 135, 134, 133, 132, 132, 131, 130, 129, 128, 128, 127, 126, 125, 125, 124, 123, 122, 122, 121, 120, 119, 119, 118, 117, 116, 115, 115, 114, 113, 112, 112, 111, 110, 109, 109, 108, 107, 106, 106, 105, 104, 103, 103, 103, 102, 101, 100, 100,  99,  98,  97,  97,  96,  95,  94,  93,  93,  92,  91,  90,  90,  89,  88,  87,  87,  86,  85,  84,  84,  83,  82,  81,  80,  80,  79,  78,  77,  77,  76,  75,  75,  75,  76,  77,  78,  78,  80,  82, 83,  84,  84,  85,  87,  88,  88,  89,  90,  91,  91,  92,  94,  95,  95,  96,  97,  98,  99,  99, 101, 102, 103, 103, 104, 105, 107, 107, 108, 109, 110, 111, 111, 112, 114, 115, 115, 116, 117, 118, 118, 119, 121, 122, 122, 123, 124, 125, 127, 127, 128, 129, 130, 130, 131, 132, 134, 134, 135, 136, 137, 137, 138, 140, 141, 142, 142, 144, 145, 146, 146, 148, 149, 150, 150, 151, 152, 154, 154, 155, 156, 157, 158, 158, 160, 161, 162, 162, 163, 164, 165, 165, 167, 168, 167, 166, 166, 164, 163, 161, 161,160, 158, 157, 157, 155, 154, 152, 152, 151, 149, 148, 146, 146, 145, 143, 142, 142, 140};
    // Aus SasView - UMDREHEN
    static int rSasV[256] = { 120, 132, 141, 141, 146, 155, 155, 159, 168, 168, 173, 182, 182, 187, 191, 200, 200, 205, 214, 214, 218, 228, 228, 232, 241, 241, 246, 250, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 254, 251, 248, 244, 241, 238, 235, 231, 228, 225, 222, 219, 215, 212, 209, 206, 202, 199, 196, 193, 190, 186, 183, 180, 177, 173, 170, 167, 164, 160, 157, 154, 151, 148, 144, 141, 138, 135, 131, 128, 125, 122, 119, 115, 112, 109, 106, 102,  99,  96,  93,  90,  86,  83,  80,  77,  73,  70,  67,  64,  60,  57,  54,  51,  48,  44,  41,  38,  35,  31,  28,  25,  22,  19,  15,  12,   9,   6,   2,   0,   0,   0,   0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};
    static int gSasV[256] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   8,   8,  11,  15,  19,  22,  26,  34,  34,  37,  45,  45,  48,  56,  56, 59,  63,  67,  71,  74,  78,  82,  85,  93,  93,  96, 104, 104, 108, 115, 115, 119, 122, 126, 130, 134, 137, 141, 145, 152, 152, 156, 163, 163, 167, 171, 174, 178, 182, 185, 189, 193, 196, 200, 204, 211, 211, 215, 219, 222, 226, 230, 234, 237, 241, 245, 248, 252, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 253, 249, 245, 241, 237, 232, 229, 225, 221, 217,213, 209, 205, 200, 197, 193, 189, 185, 181, 177, 173, 168, 165, 161, 157, 153, 149, 145, 141, 141, 133, 129, 125, 121, 116, 113, 109, 105, 100,  97,  93,  89,  89,  81,  77,  77,  68,  65,  61,  57,  52,  49,  45,  41,  36,  33,  29,  25,  25,  17,  13,  13,   4,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};
    static int bSasV[256] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   6,   9,  12,  15,  19,  22,  25,  28,  31,  35,  38,  41,  44,  48,  51,  54,  57,  60,  64,  67,  70,  73,  77,  80,  83,  86,  90,  93,  96,  99, 102, 106, 109, 112, 115, 119, 122, 125, 128, 131, 135, 138, 141, 144, 148, 151, 154, 157, 160, 164, 167, 170, 173, 177, 180, 183, 186, 190, 193, 196, 199, 202, 206, 209, 212, 215, 219, 222, 225, 228, 231, 235, 238, 241, 244, 248, 251, 254, 255,255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 250, 246, 241, 241, 232, 227, 227, 218, 214, 214, 205, 200, 196, 191, 187, 182, 182, 173, 168, 168, 159, 155, 155, 146, 141, 141, 132, 120};

    // aus dem Programm von Prof. Förster (gespeicherte Farbtabelle)
    //static int rScatter[256] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   7,  14,  20,  33,  39,  45,  52,  58,  64,  71,  77,  83,  90, 102, 109, 115, 121, 127, 134, 140, 146, 153, 159, 172, 178, 184, 191, 197, 203, 210, 216, 222, 229, 241, 248, 254, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
    //static int gScatter[256] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,   0,   0,   5,  11,  17,  24,  30,  36,  43,  49,  61,  68,  74,  80,  87,  93,  99, 106, 112, 118, 131, 137, 144, 150, 156, 163, 169, 175, 182, 188, 201, 207, 213, 220, 226, 232, 239, 245, 251, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 250, 243, 237, 231, 224, 218, 212, 199, 193, 186, 180, 174, 167, 161, 155, 148, 136, 129, 123, 117, 110, 104,  98,  91,  85,  79,  66,  60,  53,  47,  41,  34,  28,  22,  15,   9,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  8,  14,  27,  33,  40,  46,  52,  59,  65,  71,  78,  90,  97, 103, 109, 116, 122, 128, 135, 141, 147, 160, 166, 173, 179, 185, 192, 198, 204, 211, 217, 230, 236, 242, 249, 255};
    //static int bScatter[256] = {   6,  13,  19,  25,  32,  44,  51,  57,  63,  70,  76,  82,  89,  95, 101, 114, 120, 127, 133, 139, 146, 152, 158, 165, 171, 184, 190, 196, 203, 209, 215, 222, 228,234, 247, 253, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 252, 240, 233, 227, 221, 214, 208, 202, 195, 189, 183, 170, 164, 157, 151, 145, 138, 132, 126, 119, 113, 100,  94,  88,  81,  75,  69,  62,  56,  50,  37,  31,  24,  18,  12,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,  10,  16,  23,  29,  35,  42,  48,  54,  61,  73,  80,  86,  92,  99, 105, 111, 118, 124, 130, 143, 149, 156, 162, 168, 175, 181, 187, 194, 200, 213, 219, 225, 231, 238, 244, 250, 255,255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};

#ifndef CONSOLENPROG
    if ( !noGUI ) lastColorTbl = ui->cbsColorTbl->currentIndex();
#endif

    useResiduenColors = false;
    curColorSelected = _slColornames.indexOf(arg1);
    if ( curColorSelected < 0 ) curColorSelected = 4; // auch bei TPV
    switch ( curColorSelected )
    {
    // 0"Greyscale", 1"Glowing colors", 2"Earth colors", 3"Temperature", 4"Scatter-Color",
    // 5"Scatter-Colorwheel", 6"Scatter-Temperature", 7"Scatter-Geo", 8"Scatter-Desy",
    // 9"Scatter-Grey", 10"Scatter-Inverse grey", 11"Scatter-Jet",
    // 12"Residuen Glow", 13"Residuen Temp", 14"Mathematica", 15"Mathematica-Inv" };
    case  0: //if ( arg1 == tblGrey )
        farbtabelle.clear();
        for ( int c=0; c<256; c++ )
            farbtabelle.append( QColor(c,c,c) );
        break;
    case  1: //else if ( arg1 == tblGlow )
        updateColortable( rGlow, gGlow, bGlow );
        break;
    case 2: //else if ( arg1 == tblEarth )
        updateColortable( rEarth, gEarth, bEarth );
        break;
    case  3: //else if ( arg1 == tblGRtemp )
        farbtabelle.clear();
        for ( int i=0; i<512; i++ )
            farbtabelle.append( QColor(r0[i],g0[i],b0[i]) );
        for ( int i=0; i<512; i++ )
            farbtabelle.append( QColor(r1[i],g1[i],b1[i]) );
        break;
    case  4: // 4"Scatter-Color",
        farbtabelle.clear();
        for ( int i=0; i<=1785; i++ )
            farbtabelle.append( IntToColor( i/1785., 0, 1 ) );
        break;
    case  5: // 5"Scatter-Colorwheel",
        farbtabelle.clear();
        for ( int i=0; i<=1275; i++ )
            farbtabelle.append( IntToColor( i/1275., 0, 2 ) );
        break;
    case  6: // 6"Scatter-Temperature",
        farbtabelle.clear();
        for ( int i=0; i<=733; i++ )
            farbtabelle.append( IntToColor( i/733., 0, 3 ) );
        break;
    case  7: // 7"Scatter-Geo",
        farbtabelle.clear();
        for ( int i=0; i<=1337; i++ )
            farbtabelle.append( IntToColor( i/1337., 0, 4 ) );
        break;
    case  8: // 8"Scatter-Desy",
        farbtabelle.clear();
        for ( int i=0; i<=1404; i++ )
            farbtabelle.append( IntToColor( i/1404., 0, 5 ) );
        break;
    case  9: // 9"Scatter-Grey",
        farbtabelle.clear();
        for ( int i=0; i<=256; i++ )
            farbtabelle.append( IntToColor( i/256., 0, 6 ) );
        break;
    case 10: // 10"Scatter-Inverse grey",
        farbtabelle.clear();
        for ( int i=0; i<=256; i++ )
            farbtabelle.append( IntToColor( i/256., 0, 7 ) );
        break;
    case 11: // 11"Scatter-Jet",
        farbtabelle.clear();
        for ( int i=0; i<=1132; i++ )
            farbtabelle.append( IntToColor( i/1132., 0, 8 ) );
        break;
    case 12: //else if ( arg1 == tblResGlo )
        useResiduenColors = true;
        updateColortable( rGlow, gGlow, bGlow );
        break;
    case 13: //else if ( arg1 == tblResTem )
        useResiduenColors = true;
        farbtabelle.clear();
        for ( int i=0; i<512; i++ )
            farbtabelle.append( QColor(r0[i],g0[i],b0[i]) );
        for ( int i=0; i<512; i++ )
            farbtabelle.append( QColor(r1[i],g1[i],b1[i]) );
        break;
    case 14: // Mathematica
        farbtabelle.clear();
        for ( int i=0; i<256; i++ )
            farbtabelle.append( QColor(rMat[i],gMat[i],bMat[i]) );
        break;
    case 15: // Mathematica-Inv
        farbtabelle.clear();
        for ( int i=255; i>=0; i-- )
            farbtabelle.append( QColor(rMat[i],gMat[i],bMat[i]) );
        break;
    case 16: // SasView - UMDREHEN
        farbtabelle.clear();
        for ( int i=255; i>=0; i-- )
            farbtabelle.append( QColor(rSasV[i],gSasV[i],bSasV[i]) );
        break;
    }
#ifndef CONSOLENPROG
    if ( noGUI ) return;

    // Anzeigen der Farbtabelle in einem schmalen Bild. Damit es bei 1024 Farben nicht zu breit wird,
    //  wird es dann geviertelt.
    float off;
//    if ( farbtabelle.size() > 1500 ) off = 8;
//    else if ( farbtabelle.size() > 1100 ) off = 6;
//    else if ( farbtabelle.size() >  500 ) off = 4;
//    else off = 1;
    QImage img( 256 /*farbtabelle.size()/off*/, 10, QImage::Format_RGB32 );
    off = farbtabelle.size() / 256.0;
    if ( useResiduenColors )
    {
        int x, xc;
        for ( x=0, xc=255; x<img.width()/2; xc-=(img.width()/(256/2)),x++ )
            for ( int y=0; y<10; y++ )
                img.setPixelColor( x, y, QColor(xc,xc,xc) );
        for ( xc=0; xc<farbtabelle.size(); xc+=2*off,x++ )
            for ( int y=0; y<10; y++ )
                img.setPixelColor( x, y, farbtabelle.at(xc) );
    }
    else
    {
        float xc;
        int x;
        //qDebug() << farbtabelle.size() << off << img.width();
        for ( x=0, xc=0; xc<farbtabelle.size(); xc+=off,x++ )
            for ( int y=0; y<10; y++ )
                img.setPixelColor( x, y, farbtabelle.at(static_cast<int>(xc)) );
    }
    ui->lblColorTable->setPixmap( QPixmap::fromImage(img) );
    histogramValid = false;
    //if ( ui->tblMetaData->isVisible() )
    //    on_butMetaData_toggled(true);
#endif
}

void widImage::updateColortable( int *r, int *g, int *b )
{
    farbtabelle.clear();
    float rr=r[0], gg=g[0], bb=b[0];
    float stp = 256 / 9.0f + 1;
    for ( int i=1; i<10; i++ )
    {
        float dr = (r[i] - rr) / stp;
        float dg = (g[i] - gg) / stp;
        float db = (b[i] - bb) / stp;
        for ( int c=0; c<stp && farbtabelle.size()<256; c++ )
        {
            farbtabelle.append( QColor(static_cast<int>(rr),static_cast<int>(gg),static_cast<int>(bb)) );
            rr += dr; if (rr<0) rr=0; else if (rr>255) rr=255;
            gg += dg; if (gg<0) gg=0; else if (gg>255) gg=255;
            bb += db; if (bb<0) bb=0; else if (bb>255) bb=255;
        }
    }
}

#ifndef CONSOLENPROG
void widImage::on_butUpdate_clicked()
{
    if ( data1D )
    {
        imgNormal = generateImage();
        return;
    }
    qApp->setOverrideCursor(Qt::WaitCursor);
    if ( !useImageFileOnly )
        imgNormal = generateImage();
    if ( imgNormal.isNull() )
    {
        if ( !noGUI )
        {
            ui->lblImage->setText( "Invalid min/max" );
            setEnabled(false); // Komplettes Fenster gesperrt
        }
        qApp->restoreOverrideCursor();
        return;
    }

    if ( noGUI )
        imgNoGUI = imgNormal;
    else if ( useImageFileOnly )
    {
        ui->lblImage->setPixmap( QPixmap::fromImage(imgNormal) );
        //qDebug() << "img resize" << imgNormal.size();
        ui->lblImage->setMinimumSize( imgNormal.size() );

        QSize siz = imgNormal.size();
        siz += QSize(10,10);
        ui->scrollArea2d->setMinimumSize( siz );
        adjustSize();
    }
    else
    {   // Display the image
        if ( ui->cbsZoom->currentIndex() > 0 )
        {   // Displays the image with the correct zoom factor
            on_cbsZoom_activated( ui->cbsZoom->currentIndex() );
        }
        else
        {
            if ( enaExtract || enaNoFit )
            {
                extractUpdateImage(false);
            }
            else
            {
                ui->lblImage->setPixmap( QPixmap::fromImage(imgNormal) );
            }
            //qDebug() << "img resize" << imgNormal.size();
            ui->lblImage->setMinimumSize( imgNormal.size() );

            QSize siz = imgNormal.size();
            siz += QSize(10,10);
            ui->scrollArea2d->setMinimumSize( siz );
            adjustSize();
        }
#ifndef NOADJUST
#ifdef IMG_ONE_WINDOW
        parentWidget()->adjustSize();
#else
        adjustSize();
#endif
#endif
        if ( ui->tblMetaData->isVisible() )
            on_butMetaData_toggled(true);
    }
    qApp->restoreOverrideCursor();
}

QColor widImage::extractFrameColor()
{
    switch ( extractFrameCol )
    {
    case 0: return farbtabelle.first();
    case 1: return farbtabelle.last();
    case 2: return Qt::red;
    case 3:
    default: break;
    }
    return Qt::white; // default
}

void widImage::extractUpdateImage(bool flg)
{
    if ( enaExtract || enaNoFit )
    {
        QImage tmp = imgNormal;
        QPainter pnt(&tmp);
        pnt.setBrush( Qt::NoBrush );
        if ( zoom2pix[ui->cbsZoom->currentIndex()] == 0 /*640*640*/ )
            pnt.setPen( extractFrameColor() );
        else if ( zoom2pix[ui->cbsZoom->currentIndex()] < 1 )
            pnt.setPen( QPen( extractFrameColor(), 1.0/zoom2pix[ui->cbsZoom->currentIndex()]) );
        else
            pnt.setPen( extractFrameColor() );
        if ( enaExtract )
            pnt.drawRect( extractRect );
        else
        {
            for ( int i=0; i<4; i++ )
                if ( ! noFitRect[i].isEmpty() )
                    pnt.drawRect( noFitRect[i] );
        }

        if ( ui->cbsZoom->currentIndex() > 0 )
        {
            if ( zoom2pix[ui->cbsZoom->currentIndex()] == 0 /*640*640*/ )
                ui->lblImage->setPixmap( QPixmap::fromImage(tmp).scaled(640,640) );
            else
                ui->lblImage->setPixmap( QPixmap::fromImage(tmp).scaled(tmp.size()*zoom2pix[ui->cbsZoom->currentIndex()]) );
        }
        else
        {
            ui->lblImage->setPixmap( QPixmap::fromImage(tmp) );
        }
    }
    else if ( flg )
    {
        if ( ui->cbsZoom->currentIndex() > 0 )
        {
            if ( zoom2pix[ui->cbsZoom->currentIndex()] == 0 /*640*640*/ )
                ui->lblImage->setPixmap( QPixmap::fromImage(imgNormal).scaled(640,640) );
            else
                ui->lblImage->setPixmap( QPixmap::fromImage(imgNormal).scaled(imgNormal.size()*zoom2pix[ui->cbsZoom->currentIndex()]) );
        }
        else
        {
            ui->lblImage->setPixmap( QPixmap::fromImage(imgNormal) );
        }
    }
}


void widImage::switchGuiFor1d()
{
    if ( data1D )
    {   // für 1d nicht verwendete Elemente wegnehmen
        ui->cbsZoom->hide(); // eventuell doch???
        ui->cbsRotation->hide();
        ui->togMirrorH->hide();
        ui->togMirrorV->hide();
        ui->lblColorTable->hide();
        ui->lblCbsColorTbl->hide();
        ui->cbsColorTbl->hide();
        ui->scrollArea2d->hide();  //qDebug()<<"img hide";
        ui->scrollArea1d->show();
    }
    else
    {   // für 2d wieder zeigen
        ui->cbsZoom->show();
        ui->cbsRotation->show();
        ui->togMirrorH->show();
        ui->togMirrorV->show();
        ui->lblColorTable->show();
        ui->lblCbsColorTbl->show();
        ui->cbsColorTbl->show();
        ui->scrollArea1d->hide();
        ui->scrollArea2d->show();  //qDebug()<<"img show";
    }
    adjustSize();
}

#endif // !CONSOLENPROG

#ifndef NOQWT
void widImage::plot1dCurve(double vmin, double vmax, bool isLog)
{
    //qDebug() << "plot1dCurve()" << vmin << vmax << isLog << "pnts:" << points1d.size();
    bool doresize = false;
    if ( plot1d == nullptr )
    {
        doresize = true;
        plot1d = new QwtPlot( ui->widPlot );
        plot1d->setCanvasBackground( Qt::white );
        // Ein Grid zur Zeichenfläche hinzufügen
        grid1d = new QwtPlotGrid();
        grid1d->setMajorPen(QPen(Qt::DotLine));
        grid1d->attach( plot1d );
        // Achsenbeschriftungen im Default Font lassen
        plot1d->setAxisTitle( QwtPlot::xBottom, "q [nm-1]" );
        plot1d->setAxisTitle( QwtPlot::yLeft, "I(q)" );
    }
    if ( curve1d == nullptr )
    {
        curve1d = new QwtPlotCurve();
        curve1d->setPen( Qt::black, 2 ); // color and thickness in pixels
        curve1d->setRenderHint( QwtPlotItem::RenderAntialiased, true ); // use antialiasing
        curve1d->attach( plot1d );
    }
    // Titel über dem Plot aus dem Widget nehmen (ist vorher gesetzt worden)
    QString tit = ui->grpImage->title();
    if ( tit.startsWith("Image") ) tit = tit.mid(6).trimmed();
    if ( tit.isEmpty() ) tit = "Curve";
    // Den Titel etwas kleiner machen (Bold ist wohl immer)
    QFont font;
    font.setPointSize(10);
    QwtText title;
    title.setText(tit);
    title.setFont(font);
    plot1d->setTitle( title );

    // Die Punkte hinzufügen
    curve1d->setSamples( points1d );
    // Anscheinend hat Qt6 Probleme mit der Qwt-Lib wenn diese mit Qt5 erstellt wurde, diese Routine wird nicht gefunden.

    // Jetzt noch die Y-Achse logarithmisch sofern gewünscht
    plot1d->setAxisAutoScale(QwtPlot::yLeft, false);
    if ( isLog )
    {
        QwtLogScaleEngine *ylogScale = new QwtLogScaleEngine();
        double x1=vmin, x2=vmax;
        double stp=0;
        if ( vmin == 0 && vmax > 10 ) x1 = 0.0002;
        //ylogScale->autoScale( 5 /*maxNumSteps*/, x1, x2, stp );
        QwtScaleDiv yDiv = ylogScale->divideScale(x1, x2, 2/*maxMajorSteps*/, 10/*maxMinorSteps*/, stp/*stepSize*/);
        plot1d->setAxisScaleEngine(QwtPlot::yLeft, ylogScale);
        plot1d->setAxisScaleDiv(QwtPlot::yLeft, yDiv);
        qDebug() << "yLog:" << vmin << vmax << "->" << x1 << x2 << ":" << stp;
    }
    else
    {
        QwtLinearScaleEngine *ylinScale = new QwtLinearScaleEngine();
        QwtScaleDiv yDiv = ylinScale->divideScale(vmin, vmax, 10, 10, 0.0); // to calculate
        plot1d->setAxisScaleEngine(QwtPlot::yLeft, ylinScale);
        plot1d->setAxisScaleDiv(QwtPlot::yLeft, yDiv);
    }

    if ( doresize )
    {
        plot1d->adjustSize();
        QSize siz = plot1d->sizeHint();
        siz += QSize(10,10);
        //qDebug() << "plt size" << plot1d->sizeHint() << siz;
        ui->widPlot->setMinimumSize( plot1d->sizeHint() );
        ui->scrollArea1d->setMinimumSize( siz );
        adjustSize();
    }
    else
        plot1d->replot();
}
#endif


QImage widImage::generateImage()
{
    double vmin, vmin2, vmax, vlogmin, vlogmax;
#ifndef CONSOLENPROG
#ifdef UseAutoNoBorder
    bool useBorders = !noGUI && !ui->radScaleAutoNoBorder->isChecked();
#endif
    if ( !noGUI && ui->radShowLin->isChecked() )
#else
    if ( useLinOut )
#endif
    {   // Lin
        vlogmin = vlogmax = 0;
        vmin = vmin2 = vmax = data.first();
        QVector<double>::const_iterator i;
#ifdef UseAutoNoBorder
        for ( i = data.constBegin(); i != data.constEnd(); )
            for ( int xx=minX; xx<maxX; xx++ )
                for ( int yy=minY; yy<maxY; yy++, ++i )
                    if ( useBorders || ( xx>minX && xx<maxX-1 && yy>minY && yy<maxY-1 ) )
#else
        for ( i = data.constBegin(); i != data.constEnd(); ++i )
#endif
                    {
                        if ( *i > 0 && *i < vmin2 )
                            vmin2 = *i;
                        if ( *i < vmin )
                            vmin = *i;
                        if ( *i > vmax )
                            vmax = *i;
                    }
    }
    else
    {   // Log
        vlogmin=1e26; vlogmax=0;
        vmin = vmin2 = vmax = data.first();
        QVector<double>::const_iterator i;
#ifdef UseAutoNoBorder
        for ( i = data.constBegin(); i != data.constEnd(); )
            for ( int xx=minX; xx<maxX; xx++ )
                for ( int yy=minY; yy<maxY; yy++, ++i )
                    if ( useBorders || ( xx>minX && xx<maxX-1 && yy>minY && yy<maxY-1 ) )
#else
        for ( i = data.constBegin(); i != data.constEnd(); ++i )
#endif
                        if ( *i > -1.0e10 )
                        {
                            if ( *i > 0 && *i < vmin2 )
                                vmin2 = *i;
                            if ( *i < vmin )
                                vmin = *i;
                            if ( *i > vmax )
                                vmax = *i;
                            if ( *i > 0 )
                            {
                                if ( log10(*i) < vlogmin )
                                    vlogmin = log10(*i);
                                if ( log10(*i) > vlogmax )
                                    vlogmax = log10(*i);
                            }
                        }
    }
    if ( fabs(vmin-vmax) < 1e-6 )
    {   // Bei lin und log ...
        vmin *= 0.998;
        vmax *= 1.002;
        if ( vmax == 0 ) vmax = 0.2;
        //qDebug() << "change vmin/vmax" << vmin << vmax;
    }
    //qDebug() << "generateImage(): vmin/vmax" << vmin << vmax;
#ifndef CONSOLENPROG
    if ( !noGUI )
    {
        if ( vmin == 0 )
            ui->lblDispMinLin->setText( QString("(0)  %1").arg(vmin2) );
        else
            ui->lblDispMinLin->setText( QString::number(vmin) );
        ui->lblDispMaxLin->setText( QString::number(vmax) );
        ui->lblDispMinLog->setText( QString::number(vlogmin) );
        ui->lblDispMaxLog->setText( QString::number(vlogmax) );
    }
    if ( !noGUI && ui->radScaleFixed->isChecked() )
    {
        bool ok;
        double tmp;
        //if ( ui->radShowLin->isChecked() )
        //{
            tmp = input2val( ui->inpFixedMax, &ok ); if ( ok ) vmax = tmp;
            tmp = input2val( ui->inpFixedMin, &ok ); if ( ok ) vmin = tmp;
        //}
        //else
        if ( !ui->radShowLin->isChecked() )
        {
            tmp = input2val( ui->inpFixedMax, &ok ); if ( ok && tmp > 0 ) vlogmax = log10(tmp);
            tmp = input2val( ui->inpFixedMin, &ok ); if ( ok && tmp > 0 ) vlogmin = log10(tmp);
        }
    }
#endif
#ifndef NOQWT
    if ( data1D )
    {
        if ( !noGUI && ui->radShowLin->isChecked() )
        {   // Lin
            plot1dCurve( vmin, vmax, false );
        }
        else
        {   // Log
            plot1dCurve( vmin, vmax, true );    // nicht vlogmin/vlogmax!
        }
        plot1d->replot();
        QRect rc = plot1d->geometry();
        QImage img( rc.width(), rc.height(), QImage::Format_RGB32 );
        QPainter pnt(&img);
        //pnt.fillRect( rc, Qt::white );
        //qDebug() << rc;
        QwtPlotRenderer rend;
        rend.render(plot1d,&pnt,rc);
        pnt.end();

        return img;
    }
#endif // !NOQWT

    if ( vmin >= vmax ) return QImage();

    if ( farbtabelle.size() < 256 )
#ifndef CONSOLENPROG
        on_cbsColorTbl_textActivated( noGUI ? _slColornames[1]/*tblGlow*/ : ui->cbsColorTbl->currentText() );
#else
        on_cbsColorTbl_textActivated( lastColorName );
#endif

    QImage img( maxX - minX, maxY - minY, QImage::Format_RGB32 );
    //qDebug() << "GenImg" << img.size();

#ifndef CONSOLENPROG
    if ( !noGUI && ui->radShowLin->isChecked() )
#else
    if ( useLinOut )
#endif
    {   // Lin
        for ( int y=minY, yp=0; y<maxY; y++, yp++ )
            for ( int x=minX, xp=0; x<maxX; x++, xp++ )
            {
                img.setPixelColor( yp, xp, data2pixelcol(getData(x,y),vmin,vmax,vmax,false) );
                //QColor widImage::data2pixelcol( double val, double min, double max, bool islog )
            }
    }
    else
    {   // Log
        for ( int y=minY, yp=0; y<maxY; y++, yp++ )
            for ( int x=minX, xp=0; x<maxX; x++, xp++ )
            {
                img.setPixelColor( xp, yp, data2pixelcol(getData(x,y),vlogmin,vlogmax,vmax,true) );
                //QColor widImage::data2pixelcol( double val, double min, double max, bool islog )
            }
    }

    // TEST:
    //for ( int y=minY+1, yp=1; y<minY+21; y++, yp++ )
    //    img.setPixelColor( 1, yp, farbtabelle.last() );
    //for ( int x=minX+1, xp=1; x<minX+11; x++, xp++ )
    //    img.setPixelColor( xp, 1, farbtabelle.last() );

    // Perform Mirroring
    if ( swapHor || swapVert )
        img = img.mirrored( swapHor, swapVert );

    // Perform rotation
    // ACHTUNG: im Scatter ist +Y nach oben und +X nach rechts.
    //          hier ist bei mat.rotate(0) +Y nach links und +X nach unten.
    //          Daher wird das Bild immer um 270° gedreht, dann passen beide wieder.
    int locrot = (rotIndex + 3) % 4;
    if ( locrot == 0 /*0°*/ ) return img;     // to speed up a little bit
    QTransform mat = QTransform(); // Identity
    return img.transformed( mat.rotate( locrot * 90 ) );
}


#ifndef CONSOLENPROG
double widImage::input2val( QLineEdit *inp, bool *ok )
{
    // TODO: Alternative: bei allen Inputfeldern eine passende Maske setzen, damit kein , eingegeben werden kann
    QString val = inp->text();
    val = val.replace(",",".");
    return val.toDouble(ok);
}
#endif


QColor widImage::data2pixelcol( double val, double min, double max, double datamax, bool islog )
{
    if ( useResiduenColors )
    {   // Spezielle Anzeigeform
        if ( val >= 0 )
        {   // Glüh- oder Temp-Farben, sprich: aktuelle Farbtabelle
            double relval = val / max;     // Min ist ja 0
            int pixval = static_cast<int>( relval * farbtabelle.size() );
            if ( pixval < 0 ) pixval = 0;
            else if ( pixval >= farbtabelle.size() ) pixval = farbtabelle.size()-1;
            return farbtabelle.at(pixval);
        }
        else
        {   // Grauwerte (invertiert, d.h. weiß ist Minimum)
            double relval = (val - min) / (0-min);  // Max ist 0
            int pixval = 256 - static_cast<int>( relval * 256 );
            if ( pixval < 0 ) pixval = 0;
            else if ( pixval >= 255 ) pixval = 255;
            return QColor(pixval,pixval,pixval);
        }
    }

    if ( curColorSelected >= 4 && curColorSelected <= 11 )
    {   // Spezielle Scatter-Farben
        if ( isnan(val) ) return QColor(0,0,0);
        if ( islog )
        {
            double logdecade = max - min;
            return IntToColor( val/datamax, logdecade, curColorSelected-3 );
        }
        return IntToColor( val/datamax, 0, curColorSelected-3 );
    }

    double relval;
    if ( islog )
    {
        if ( val > 0 )
            relval = (log10(val) - min) / (max-min);
        else
            relval = 0;
    }
    else
        relval = (val - min) / (max-min);
    int pixval = static_cast<int>( relval * farbtabelle.size() );
    if ( pixval < 0 ) pixval = 0;
    else if ( pixval >= farbtabelle.size() ) pixval = farbtabelle.size()-1;
    return farbtabelle.at(pixval);
}


#ifndef CONSOLENPROG
void widImage::saveParams( QSettings &sets )
{
    sets.setValue( "ScaleFix",   ui->radScaleFixed->isChecked() ); // radScaleAuto = !radScaleFixed
    sets.setValue( "DispLin",    ui->radShowLin->isChecked() ); // radShowLog = !radShowLin
    sets.setValue( "FixedMin",   ui->inpFixedMin->text() );
    sets.setValue( "FixedMax",   ui->inpFixedMax->text() );
    sets.setValue( "ColorTable", ui->cbsColorTbl->currentText() );
}

void widImage::loadParams( QSettings &sets )
{
    ui->radScaleFixed->setChecked( sets.value("ScaleFix").toBool() ); // radScaleAuto = !radScaleFixed
    ui->radShowLin->setChecked( sets.value("DispLin").toBool() ); // radShowLog = !radShowLin
    ui->inpFixedMin->setText( sets.value("FixedMin").toString() );
    ui->inpFixedMax->setText( sets.value("FixedMax").toString() );
    ui->cbsColorTbl->setCurrentIndex( ui->cbsColorTbl->findText( sets.value("ColorTable").toString() ) );
}

void widImage::on_butSave_clicked()
{
    QString fn;
    if ( data1D )
    {
        fn = QFileDialog::getSaveFileName( this, "Save", lastSavedFilePath, "1d images (*.png)" );
        if ( fn.isEmpty() ) return;
        if ( !fn.endsWith(".png",Qt::CaseInsensitive) ) fn += ".png";
    }
    else
    {
        fn = QFileDialog::getSaveFileName( this, "Save", lastSavedFilePath, "Datafiles (*.dat)" );
        if ( fn.isEmpty() ) return;
        if ( !fn.endsWith(".dat",Qt::CaseInsensitive) ) fn += ".dat";
    }
    lastSavedFilePath = fn;
    QString fnbase = fn.left(fn.length()-4);
    QFile f(fn);
    ui->lblFilename->setText( fn );
    ui->lblFilename->setEnabled(true);

    if ( ! data1D )
    {
        // ----- Save the image in binary format
        union
        {
            double e;
            char c[1];
        } val;
        if ( f.open(QIODevice::WriteOnly) )
        {
            // (1) 19 Bytes as text with the current dimensions of the image.
            QString tmp = QString("%1 %2 %3 %4").arg(minX,4).arg(maxX,4).arg(minY,4).arg(maxY,4);
            f.write( qPrintable(tmp) );
            // (2) then all values in binary format.
            for ( int y=minY; y<maxY; y++ )
            {
                for ( int x=minX; x<maxX; x++ )
                {
                    val.e = getData(x,y);
                    f.write( val.c, sizeof(double) );
                }
            }
            f.close();
        }
        else
            qDebug() << f.fileName() << f.errorString();

        // ----- Save the image in textual format
        fn = fnbase + ".csv";
        f.setFileName(fn);
        if ( f.open(QIODevice::WriteOnly) )
        {
            QTextStream ts(&f);
            // (1) the first two lines contains the dimensions.
            ts << "X:;" << minX << ";" << maxX << EOL;
            ts << "Y:;" << minY << ";" << maxY << EOL;
            // (2) then all values line by line separated with a ; so this can be read by Excel.
            for ( int y=minY; y<maxY; y++ )
            {
                ts << "Y=" << y << ":";
                for ( int x=minX; x<maxX; x++ )
                {
                    ts << ";" << getData(x,y);
                }
                ts << EOL;
            }
            f.close();
        }
        else
            qDebug() << f.fileName() << f.errorString();
    } // !data1D

    // ----- Save the image in PNG format
    saveImage(fnbase+".png");

    // ----- Save the MetaData in textual format
    fn = fnbase + ".txt";
    f.setFileName(fn);
    if ( f.open(QIODevice::WriteOnly) )
    {
        QTextStream ts(&f);
        for ( int r=0; r<ui->tblMetaData->rowCount(); r++ )
        {
            ts << ui->tblMetaData->item(r,0)->text() << "\t" << ui->tblMetaData->item(r,1)->text() << EOL;
        }
        f.close();
    }
    else
        qDebug() << f.fileName() << f.errorString();

    if ( !data1D )
    {
        // ----- Save NoFitRect
        fileInfos.filePath = lastSavedFilePath;
        saveNoFitRects();
    }
}

void widImage::on_butMetaData_toggled(bool checked)
{
    ui->tblMetaData->setVisible(checked);
    ui->lblHistogram->setVisible(checked);
    ui->butMetaData->setText( checked ? "<<<" : "Meta >>>" );
    if ( checked )
    {
        ui->tblMetaData->resizeColumnsToContents();
        if ( data1D )
            ui->lblHistogram->hide();
        else
        {
            ui->lblHistogram->show();
            if ( histoThread == nullptr )
            {
                histoThread = new calcHistorgramm(this);
                connect( histoThread, SIGNAL(setHistoPixmap(QPixmap)), this, SLOT(setHistogram(QPixmap)) );
            }
            if ( ! histogramValid )
            {   // Muss ja nur neu berechnet werden, wenn neue Daten da sind.
                ui->lblHistogram->setPixmap(QPixmap());
                ui->lblHistogram->setText("generating histogram ...");
                histoThread->start();   // Bei EDF-Dateien (1920*1920) dauert es doch schon bis zu 10 sec.
            }
        }
    }
    else
        qApp->processEvents();
//#ifndef NOADJUST
#ifdef IMG_ONE_WINDOW
    parentWidget()->adjustSize();
#else
    adjustSize();
#endif
//#endif
}

QPixmap widImage::getColorTable()
{
    // Qt-Versionen:
    //  Büro   qmake     5.9.7  (Anaconda)
    //         QtCreator 5.15.2 (/usr/bin/qmake-qt5)
    //  Laptop qmake     <nicht im Pfad>
    //         QtCreator 5.14.2
#if (QT_VERSION <= QT_VERSION_CHECK(5, 14, 2))
    return *ui->lblColorTable->pixmap();
#else
    return ui->lblColorTable->pixmap(Qt::ReturnByValue);
#endif
}

QImage widImage::getCurrentImage()
{
    if ( noGUI )
        return imgNoGUI;
    return imgNormal;
}

calcHistorgramm::calcHistorgramm( widImage *i ) // QThread
{
    img = i;
}

void calcHistorgramm::run() // QThread
{
    QPixmap coltbl = img->getColorTable();
    QImage src = img->getCurrentImage();
    QVector<QColor> farbtabelle = img->getFarbtabelle();
    QPixmap dst( coltbl.width(), coltbl.height() + HISTO_SY );
    QPainter pnt(&dst);
    pnt.drawPixmap( 0, HISTO_SY, coltbl );
    pnt.fillRect( 0, 0, coltbl.width(), HISTO_SY, Qt::white );
    float facx = farbtabelle.size() / coltbl.width();
    int cnt[static_cast<int>(farbtabelle.size()/facx+1)];
    memset( cnt, 0, static_cast<int>(farbtabelle.size()/facx+1)*sizeof(int) );
    int cntmax = 0;
    //int indmax = -1;
    for ( int x=0; x<src.width(); x++ )
        for ( int y=0; y<src.height(); y++ )
        {
            int p = farbtabelle.indexOf( QColor(src.pixel(x,y)) );
            if ( p >= 0 )
            {
                int pi = p / facx;
                cnt[pi]++;
                if ( cnt[pi] > cntmax ) { cntmax = cnt[pi]; /*indmax = p/facx;*/ }
            }
        }
    int max2 = 0;
    for ( int x=0; x<farbtabelle.size()/facx; x++ )
        if ( cnt[x] > max2 && cnt[x] < cntmax ) max2 = cnt[x];
    float facy = (1.0 * HISTO_SY ) / max2;
    pnt.setPen( Qt::black );
    for ( int x=0; x<farbtabelle.size()/facx; x++ )
    {
        int y = cnt[x] * facy;
        if ( y < 0 ) y=0;
        else if ( y >= HISTO_SY ) y=HISTO_SY-1;
        else if ( y == 0 && cnt[x] > 0 ) y++;
        if ( y > 0 ) pnt.drawLine( x, HISTO_SY-y, x, HISTO_SY );
    }
    emit setHistoPixmap(dst);
}

void widImage::setHistogram( QPixmap dst )
{
    ui->lblHistogram->setPixmap( dst );
    histogramValid = true;
}


/**
 * @brief widImage::on_cbsZoom_activated [SLOT]
 * @param index - Zoomfactor as an index (0=1*, 1=2*, 2=4*) selected from the combobox
 *                and Factor set directly from the UpdateButton Slot (-1=2*, -2=4*) - NOT USED
 */
void widImage::on_cbsZoom_activated(int index)
{
    if ( bRekursiv ) return;
    zoomFactor = index;
#ifndef NOADJUST
    QPoint pnt = pos();
#endif
    //qDebug() << "img zoom";
    QSize siz;
    if ( zoom2pix[index] == 0 /*640*640*/ )
    {
        ui->lblImage->setMinimumSize( 640, 640 );
        siz = QSize(640,640);
    }
    else
    {
        ui->lblImage->setMinimumSize( imgNormal.size()*zoom2pix[index] );
        siz = imgNormal.size()*zoom2pix[index];
    }
    siz += QSize(10,10);
    if ( siz.width() > 1000 ) siz = QSize(1000,1000);
    ui->scrollArea2d->setMinimumSize( siz );
    //adjustSize();


    if ( enaExtract || enaNoFit )
    {
        extractUpdateImage(false);
    }
    else if ( ! imgNormal.isNull() )
    {
        QPixmap pix = QPixmap::fromImage(imgNormal);
        if ( zoom2pix[index] == 0 /*640*640*/ )
            ui->lblImage->setPixmap( pix.scaled(640,640) );
        else
            ui->lblImage->setPixmap( pix.scaled(pix.size()*zoom2pix[index]) );
    }
    qApp->processEvents();
#ifndef NOADJUST
#ifdef IMG_ONE_WINDOW
    parentWidget()->adjustSize();
    move( pnt );
    parentWidget()->adjustSize();
#else
    adjustSize();
    //qDebug() << "MOVE cbsZoom" << pnt << windowTitle() << sender();
    if ( sender() != nullptr )
    {
        move( pnt );
        adjustSize();
        qApp->processEvents();
    }
#endif
#endif
}

void widImage::on_cbsRotation_activated( int index )
{
    if ( bRekursiv ) return;
    rotIndex = index;
    // Die Drehung wird in generateImage() berücksichtigt
    on_cbsZoom_activated( ui->cbsZoom->currentIndex() );
}

void widImage::on_togMirrorH_toggled( bool checked )
{
    if ( bRekursiv ) return;
    swapHor = checked;
    // Die Spiegelung wird in generateImage() berücksichtigt
    on_cbsZoom_activated( ui->cbsZoom->currentIndex() );
}

void widImage::on_togMirrorV_toggled( bool checked )
{
    if ( bRekursiv ) return;
    swapVert = checked;
    // Die Spiegelung wird in generateImage() berücksichtigt
    on_cbsZoom_activated( ui->cbsZoom->currentIndex() );
}

bool widImage::mouse2data( QPoint pos, int &x, int &y, QPoint &pnt )
{
#if (QT_VERSION <= QT_VERSION_CHECK(5, 14, 2))
    QRect rp = ui->lblImage->pixmap()->rect();  // Rechteck des Images
#else
    QRect rp = ui->lblImage->pixmap(Qt::ReturnByValue).rect();  // Rechteck des Images
#endif
    QSize si = ui->lblImage->size();            // Pixelgröße des Labels komplett
    //QSize si = ui->scrollArea->size();            // Pixelgröße des Labels komplett
    int dx = si.width() / 2 - rp.width() / 2;
    int dy = si.height() / 2 - rp.height() / 2;
    QRect rpi = rp.adjusted( dx, dy, dx, dy );    // Imagerechteck in Labelkoordinaten
    if ( ! rpi.contains(pos) ) return false;

    QLineF line( rpi.center(), pos );             // Linie vom Mittelpunkt des Rechtecks zum Mauspunkt

    // In <pnt> werden die Koordinaten für die Extraktionsauswahl mit der Maus
    //  gespeichert und sollten nicht gedreht werden.
    if ( zoom2pix[zoomFactor] == 0 /*640*640*/ )
    {
        pnt.setX( (pos.x() - rpi.left()) );
        pnt.setY( (pos.y() - rpi.top() ) );
        if ( ! imgNormal.isNull() )
        {
            pnt.setX( pnt.x() * imgNormal.width()  / 640.0 );
            pnt.setY( pnt.y() * imgNormal.height() / 640.0 );
        }
    }
    else
    {
        pnt.setX( (pos.x() - rpi.left()) / zoom2pix[zoomFactor] );
        pnt.setY( (pos.y() - rpi.top() ) / zoom2pix[zoomFactor] );
    }

    // The image is first mirrored, then rotated. So we must first rotate back then mirror back...
    // The rotation is performed around the center of the image!

    // Perform rotation
    line.setAngle( line.angle() + 90 * rotIndex );

    // Perform Mirroring
    if ( swapHor )
        line.setAngle( 180 - line.angle() );
    if ( swapVert )
        line.setAngle( - line.angle() );

    if ( zoom2pix[zoomFactor] == 0 /*640*640*/ )
    {
        if ( imgNormal.isNull() )
        {
            x = minX + (line.p2().x() - rpi.left());
            y = minY + (line.p2().y() - rpi.top() );
        }
        else
        {
            x = minX + (line.p2().x() - rpi.left()) * imgNormal.width()  / 640.0;
            y = minY + (line.p2().y() - rpi.top() ) * imgNormal.height() / 640.0;
        }
    }
    else
    {
        x = minX + (line.p2().x() - rpi.left()) / zoom2pix[zoomFactor];
        y = minY + (line.p2().y() - rpi.top())  / zoom2pix[zoomFactor];
    }

    return true;
}

void widImage::mousePressEvent(QMouseEvent *event)
{
    if ( !enaExtract && !enaNoFit ) return;
    if ( event->button() != Qt::LeftButton ) return;

#if (QT_VERSION <= QT_VERSION_CHECK(6,0,0))
    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPos() );
#else
    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPosition() ).toPoint();
#endif
    int x, y;
    QPoint pnt;
    if ( !mouse2data( pos, x, y, pnt ) )
    {
        extractDrawByMouse = false;
        return;
    }
    extractDrawByMouse = true;
    extractRect.setTopLeft( pnt );
    extractRect.setBottomRight( pnt );
    emit extractUpdateRect(extractRect);
    extractMouseStart = pnt;
    extractUpdateImage(false);
}

void widImage::mouseReleaseEvent(QMouseEvent *event)
{
    if ( !enaExtract && !enaNoFit ) return;
    if ( event->button() != Qt::LeftButton || !extractDrawByMouse ) return;

#if (QT_VERSION <= QT_VERSION_CHECK(6,0,0))
    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPos() );
#else
    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPosition() ).toPoint();
#endif
    int x, y;
    QPoint pnt;
    if ( !mouse2data( pos, x, y, pnt ) )
    {
        extractDrawByMouse = false;
        return;
    }
    extractDrawByMouse = false;
    extractRect.setLeft(   qMin(extractMouseStart.x(), pnt.x()) );
    extractRect.setRight(  qMax(extractMouseStart.x(), pnt.x()) );
    extractRect.setTop(    qMin(extractMouseStart.y(), pnt.y()) );
    extractRect.setBottom( qMax(extractMouseStart.y(), pnt.y()) );
    emit extractUpdateRect(extractRect);
    extractUpdateImage(false);
}

void widImage::mouseMoveEvent(QMouseEvent *event)
{
    if ( (event->buttons() & Qt::LeftButton) == 0 ) return;

#if (QT_VERSION <= QT_VERSION_CHECK(6,0,0))
    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPos() );
#else
    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPosition() ).toPoint();
#endif

    int x, y;
    QPoint pnt;
    if ( !mouse2data( pos, x, y, pnt ) )
    {
        //QToolTip::hideText();
        return;
    }

    if ( (enaExtract || enaNoFit) &&
         (event->buttons() & Qt::LeftButton) && extractDrawByMouse )
    {
        extractRect.setLeft(   qMin(extractMouseStart.x(), pnt.x()) );
        extractRect.setRight(  qMax(extractMouseStart.x(), pnt.x()) );
        extractRect.setTop(    qMin(extractMouseStart.y(), pnt.y()) );
        extractRect.setBottom( qMax(extractMouseStart.y(), pnt.y()) );
        emit extractUpdateRect(extractRect);
        extractUpdateImage(false);
        return; // Zeichnen des Mausrahmens und der Tooltip geht nicht gut zusammen.
    }

    int idx = XY2IDX(minX,maxX,minY,maxY,x,y);
    if ( idx < 0 || idx >= data.size() ) return;

    double qx, qy, qz;
    double q = xy2q( x, y, qx, qy, qz );
/*
    if ( isinf(q) )
        qDebug() << "P" << fileInfos.pixWidthX << fileInfos.pixWidthY
                 << "C" << fileInfos.centerX   << fileInfos.centerY
                 << "D" << fileInfos.distance
                 << "W" << fileInfos.wavelen;
*/
    //if ( ui->radShowLin->isChecked() )
    //{   // Lin
#if (QT_VERSION <= QT_VERSION_CHECK(6,0,0))
    QToolTip::showText( event->globalPos(),
#else
    QToolTip::showText( event->globalPosition().toPoint(),
#endif
                        QString("[%1,%2]=%3 (%8,%9)\nqx=%4,qy=%5,qz=%6)\nq=%7  %10").arg(x).arg(y).arg(data.at(idx))
                           .arg(qx).arg(qy).arg(qz).arg(q).arg(pos.x()).arg(pos.y())
                           .arg(useCenterBeamPos?"BS":"Mid"));
    //}
    //else
    //{   // Log
    //    QToolTip::showText( event->globalPos(),
    //                        QString("[%1,%2]=%3").arg(x).arg(y).arg(myLog10(data.at(idx))) );
    //}
/*
If mouse tracking is disabled (the default), the widget only receives mouse
move events when at least one mouse button is pressed while the mouse is
being moved.
If mouse tracking is enabled, the widget receives mouse move events even
if no buttons are pressed.
Access functions:
bool hasMouseTracking() const
void setMouseTracking(bool enable)

If you want to show a tooltip immediately, while the mouse is moving
(e.g., to get the mouse coordinates with QMouseEvent::pos() and show
them as a tooltip), you must first enable mouse tracking as described
above. Then, to ensure that the tooltip is updated immediately, you
must call QToolTip::showText() instead of setToolTip() in your
implementation of mouseMoveEvent().
*/
}


double widImage::xy2q( int x, int y, double &qx, double &qy, double &qz )
{
    if ( useCenterBeamPos )
    {   // Einrechnen des Beamstops (d.h. Verschiebung des Zentrums)
        double mdet = (y/*+CALC.zmax+1*/)*fileInfos.pixelY/(double)(maxY-minY);      /* mth pixel */
        double ndet = (x/*+CALC.zmax+1*/)*fileInfos.pixelX/(double)(maxX-minX);      /* nth pixel */
        // Im Pascal-Programm ist der Beamstop für Null in der Ecke.
        // Hier ist der Beamstop für Null in der Mitte gerechnet.
        double xdet = fileInfos.pixWidthX * (ndet - fileInfos.centerX);
        double ydet = fileInfos.pixWidthY * (mdet - fileInfos.centerY);
        double rdet = sqrt(xdet*xdet+ydet*ydet);
        double phidet = atan2(ydet,xdet);
        double thetadet = atan2(rdet,fileInfos.distance);
        qx = 2*M_PI*cos(phidet)*sin(thetadet)/fileInfos.wavelen;
        qy = 2*M_PI*sin(phidet)*sin(thetadet)/fileInfos.wavelen;
        qz = 2*M_PI*(1-cos(thetadet))/fileInfos.wavelen;
    }
    else
    {   // Den Mittelpunkt nutzen (ihex=0,i=0)
        if ( !qmaxset )
        {
            qmaxset = true;
            if ( qmaxcalc ) qmax = calcQmax;
            if ( qmaxedit ) qmax = editQmax;
        }
        qx = qmax * x / ((maxX-minX)/2.0); // lamu
        qy = qmax * y / ((maxY-minY)/2.0);  // lamv
        qz = 1e-20;
    }
    return sqrt(qx*qx+qy*qy+qz*qz)+eps9;
}


/**
 * @brief widImage::on_lblImage_customContextMenuRequested
 *        Click on right mouse button inside the label (image) window.
 *        This is used to mark the beam center position.
 * @param pos - position of the click in local window coordinates.
 */
void widImage::on_lblImage_customContextMenuRequested(const QPoint &pos)
{
    int centerX, centerY;
    QPoint pnt;
    if ( !mouse2data( pos, centerX, centerY, pnt ) ) return;

    fileInfos.centerX = centerX;
    fileInfos.centerY = centerY;

    // Update Metadatatable
    bool metafound = false;
    for ( int r=0; r<ui->tblMetaData->rowCount(); r++ )
    {
        if ( ui->tblMetaData->item(r,0)->text() == "BeamPosX" )
        {
            metafound = true;
            ui->tblMetaData->item(r,1)->setText( QString::number(fileInfos.centerX) );
        }
        else if ( ui->tblMetaData->item(r,0)->text() == "BeamPosY" )
        {
            metafound = true;
            ui->tblMetaData->item(r,1)->setText( QString::number(fileInfos.centerY) );
        }
    }
    if ( ! metafound )
    {
        int row = ui->tblMetaData->rowCount();
        ui->tblMetaData->setRowCount( row + 1 );
        ui->tblMetaData->setItem( row, 0, new QTableWidgetItem("BeamPosX") );
        ui->tblMetaData->setItem( row, 1, new QTableWidgetItem(QString::number(fileInfos.centerX)) );
        row = ui->tblMetaData->rowCount();
        ui->tblMetaData->setRowCount( row + 1 );
        ui->tblMetaData->setItem( row, 0, new QTableWidgetItem("BeamPosY") );
        ui->tblMetaData->setItem( row, 1, new QTableWidgetItem(QString::number(fileInfos.centerY)) );
    }
}

void widImage::getVarScaling( double &min, double &max )
{
    if ( ui->lblDispMinLin->text().startsWith("(0)") )
        min = ui->lblDispMinLin->text().mid(4).toDouble();
    else
        min = ui->lblDispMinLin->text().toDouble();
    max = ui->lblDispMaxLin->text().toDouble();
}

void widImage::setFixScaling( double min, double max )
{
    ui->inpFixedMin->setText( QString::number(min) );
    ui->inpFixedMax->setText( QString::number(max) );
    ui->radScaleFixed->setChecked(true);
    on_butUpdate_clicked();
}

void widImage::setLogScaling( bool on )
{
    ui->radShowLog->setChecked( on);
    ui->radShowLin->setChecked(!on);
    on_butUpdate_clicked();
}


#else

void widImage::setConfig( QString col, bool sh, bool sv, int r, int z, bool lin )
{
    lastColorName = col;
    swapHor  = sh;
    swapVert = sv;
    rotIndex = r;
    zoomFactor = z;
    useLinOut = lin;
}
#endif


/**
 * @brief widImage::IntToColor
 * @param val    - Datenwert
 * @param decade - 0=Linear, >0=Log-Decade
 * @param tbl    - Farbtabelle: 1=Color, 2=Color Wheel, 3=Temperature, 4=Geo, 5=Desy, 6=Grey, 7=Inverse grey, 8=Jet
 * @return passender Farbwert des Datenwertes
 */
QColor widImage::IntToColor( double ywert, double decade, int tbl )
{
    int farbwert,greyvalue;
    if ( decade < 0 || isnan(ywert) ) return QColor(0,0,0);

    switch ( tbl )
    {
    case 1:     /* color */
        if ( decade > 0 )
        {
            if ( ywert <= exp(-decade*log(10.)) )
                return QColor(0,0,0);
            else
                farbwert=trunc((1785/decade)*(decade+log(ywert)/log(10.)));
        }
        else
            farbwert=round(ywert*1785);
        if ( farbwert < 0 ) return QColor(0,0,0);
        if (farbwert<=255                  ) return QColor( 0, 0, farbwert );               /*** 0*R+0*G+i*B schwarz->blau  ***/
        if (farbwert> 255 && farbwert<= 510) return QColor( 0, farbwert-255, 255 );         /*** 0*R+i*G+  B blau->türkis   ***/
        if (farbwert> 510 && farbwert<= 765) return QColor( 0, 255, 255-(farbwert-510) );   /*** 0*R+  G-i*B türkis>grün   ***/
        if (farbwert> 765 && farbwert<=1020) return QColor( farbwert-765, 255, 0 );         /*** i*R+  G+0*B grün->gelb   ***/
        if (farbwert>1020 && farbwert<=1275) return QColor( 255, 255-(farbwert-1020), 0 );  /***   R-i*G+0*B gelb->rot   ***/
        if (farbwert>1275 && farbwert<=1530) return QColor( 255, 0, farbwert-1275 );        /***   R+0*G+i*B rot->rosa  ***/
        if (farbwert>1530 && farbwert<=1785) return QColor( 255, farbwert-1530, 255 );      /***   R+i*G+  B rosa->weiss ***/
        if (farbwert>1785                  ) return QColor( 255, 255, 255 );                /* weiss */
        break; // to avoid compiler warning

    case 2:     /* color wheel */
        if ( decade>0 )
        {
            if (ywert<=exp(-decade*log(10.)))
                return QColor(0,0,0);
            else
                farbwert=trunc((1275/decade)*(decade+log(ywert)/log(10.)));
        }
        else
            farbwert=round(ywert*1275);
        if ( farbwert < 0 ) return QColor(0,0,0);
        if (farbwert<=255                  ) return QColor( 255-farbwert, 0, 255 );         /*** i*R+0*G+  B violett->blau  ***/
        if (farbwert> 255 && farbwert<= 510) return QColor( 0, farbwert-255, 255 );         /*** 0*R+i*G+  B blau->türkis   ***/
        if (farbwert> 510 && farbwert<= 765) return QColor( 0, 255, 255-(farbwert-510) );   /*** 0*R+  G-i*B türkis>grün   ***/
        if (farbwert> 765 && farbwert<=1020) return QColor( farbwert-765, 255, 0 );         /*** i*R+  G+0*B grün->gelb   ***/
        if (farbwert>1020 && farbwert<=1275) return QColor( 255, 255-(farbwert-1020), 0 );  /***   R-i*G+0*B gelb->rot   ***/
        if (farbwert>1275                  ) return QColor( 255, 255, 255 );
        break; // to avoid compiler warning

    case 3:     /* temperature */
        if ( decade>0 )
        {
            if (ywert<=exp(-decade*log(10.)))
                return QColor(0,0,0);
            else
                farbwert=trunc((733/decade)*(decade+log(ywert)/log(10.)));
        }
        else
            farbwert=round(ywert*733);
        if ( farbwert < 0 ) return QColor(0,0,0);
        if (farbwert<=223                 ) return QColor( farbwert, 0, 0 );
        if (farbwert> 223 && farbwert<=350) return QColor( 223+round((farbwert-223)/4.), farbwert-223, round((farbwert-223)/4.) );
        if (farbwert> 350 && farbwert<=478) return QColor( 255, farbwert-223, round((478-farbwert)/4.) );
        if (farbwert> 478 && farbwert<=733) return QColor( 255, 255, farbwert-478 );
        if (farbwert> 733                 ) return QColor( 255, 255, 255 );
        break; // to avoid compiler warning

    case 4:     /* geo */
        if ( decade>0 )
        {
            if (ywert<=exp(-decade*log(10.)))
                return QColor(0,0,0);
            else
                farbwert=trunc((1337/decade)*(decade+log(ywert)/log(10.)));
        }
        else
            farbwert=round(ywert*1337);
        if ( farbwert < 0 ) return QColor(0,0,0);
        if (farbwert<=255                  ) return QColor( 0, 0, farbwert );
        if (farbwert> 255 && farbwert<= 446) return QColor( farbwert-255, farbwert-255, 255 );
        if (farbwert> 446 && farbwert<= 637) return QColor( 191-round((farbwert-446)/3.), 191+round((farbwert-446)/6.), 255-(farbwert-446) );
        if (farbwert> 637 && farbwert<= 764) return QColor( 127+(farbwert-637), 223+round((farbwert-637)/4.), round(63+(farbwert-637)/2.) );
        if (farbwert> 764 && farbwert<= 891) return QColor( 255, 255-round((farbwert-764)/2.), 127-(farbwert-764) );
        if (farbwert> 891 && farbwert<=1146) return QColor( 255-round((farbwert-891)/8.), 191-round((farbwert-891)/2.), farbwert-891 );
        if (farbwert>1146 && farbwert<=1337) return QColor( 223+round((farbwert-1146)/6.), 63+(farbwert-1146), 255 );
        if (farbwert>1337                  ) return QColor( 255, 255, 255 );
        break; // to avoid compiler warning

    case 5:     /* DESY */
        if ( decade>0 )
        {
            if (ywert<=exp(-decade*log(10.)))
                return QColor(0,0,0);
            else
                farbwert=trunc((1404/decade)*(decade+log(ywert)/log(10.)));
        }
        else
            farbwert=round(ywert*1404);
        if ( farbwert < 0 ) return QColor(0,0,0);
        if (farbwert<=192                  ) return QColor( 255-round(farbwert/6), 255-farbwert, 255 );
        if (farbwert> 192 && farbwert<= 447) return QColor( 223-round((farbwert-192)/8), round(63+(farbwert-192)/4), 255+192-farbwert );
        if (farbwert> 447 && farbwert<= 511) return QColor( 191+farbwert-447, 127+farbwert-447, 0 );
        if (farbwert> 511 && farbwert<= 638) return QColor( 255, 191+round((farbwert-511)/2), farbwert-511 );
        if (farbwert> 638 && farbwert<= 766) return QColor( 255-farbwert+638, 255-round((farbwert-638)/4), 127-round((farbwert-638)/2) );
        if (farbwert> 766 && farbwert<= 958) return QColor( 127+round((farbwert-766)/3), 223-round((farbwert-766)/6), 63+farbwert-766 );
        if (farbwert> 958 && farbwert<=1149) return QColor( 191-farbwert+958, 191-farbwert+958, 255 );
        if (farbwert>1149 && farbwert<=1404) return QColor( 0, 0, 255-farbwert+1149 );
        if (farbwert>1404                  ) return QColor( 0, 0, 0 );
        break; // to avoid compiler warning

    case 6:     /* grey nur 0..1 */
        if ( ywert>1 )
            return QColor( 255, 255, 255 );
        if ( ywert<=0 )
            return QColor(0,0,0);
        if ( decade>0 )
            greyvalue=trunc(255/decade*(decade+log(ywert)/log(10.)));
        else
            greyvalue=round(255*ywert);
        return QColor( greyvalue, greyvalue, greyvalue );

    case 7:     /* inverse grey nur 0..1 */
        if ( ywert>1 )
            return QColor(0,0,0);
        if ( ywert<=0 )
            return QColor( 255, 255, 255 );
        if ( decade>0 )
            greyvalue=trunc(255/decade*(decade+log(ywert)/log(10.)));
        else
            greyvalue=round(255*ywert);
        return QColor( 255-greyvalue, 255-greyvalue, 255-greyvalue );

    case 8:     /* jet */
        if ( decade>0 )
        {
            if (ywert<=exp(-decade*log(10.)))
                return QColor(0,0,0);
            else
                farbwert=trunc((1132/decade)*(decade+log(ywert)/log(10.)));
        }
        else
            farbwert=round(ywert*1132);
        if ( farbwert < 0 ) return QColor(0,0,0);
        if (farbwert<=112                  ) return QColor( 0, 0, 143+farbwert );           /*** dunkelblau > blau  ***/
        if (farbwert> 112 && farbwert<= 367) return QColor( 0, farbwert-112, 255 );         /*** blau + grün   ***/
        if (farbwert> 367 && farbwert<= 622) return QColor( 0, 255, 255-(farbwert-367) );   /*** -blau ***/
        if (farbwert> 622 && farbwert<= 877) return QColor( farbwert-622, 255, 0 );         /*** + rot   ***/
        if (farbwert> 877 && farbwert<=1132) return QColor( 255, 255-(farbwert-877), 0 );   /***   -grün   ***/
        if (farbwert>1132                  ) return QColor( 255, 255, 255 );
        break; // to avoid compiler warning

    } // switch tbl
    return QColor(0,0,0);
} /* IntToColor() */
