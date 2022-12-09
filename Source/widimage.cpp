#include "widimage.h"
#ifndef CONSOLENPROG
#include "ui_widimage.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QToolTip>
#include <QPainter>
#endif
#include <QFile>
#include <QDebug>
#include <math.h>
#include <iostream>


#define NOADJUST
// Wenn definiert, werden alle adjustSize() deaktiviert. Dadurch sollte die Größe
// des Fensters stabil bleiben. Durch die ScrollArea um das Image-Label ist die
// MinimumSize des Labels nicht mehr relevant für das Fenster.


#define SETT_APP   "JCNS-1-SasCrystal"  // auch in sc_maingui.h
#define SETT_IMG   "ImageScaling"

#define HISTO_SY  60


QString widImage::lastSavedFilePath = "";
int widImage::lastColorTbl = 1;

#ifndef CONSOLENPROG
static float zoom2pix[] = { 1, 2, 4, 0.5, 0.25 };
#endif

QStringList widImage::_slColornames = { /* 0*/"Greyscale",     /* 1*/"Glowing colors",     /* 2*/"Earth colors",         /* 3*/"Temperature",
                                        /* 4*/"Scatter-Color", /* 5*/"Scatter-Colorwheel", /* 6*/"Scatter-Temperature",  /* 7*/"Scatter-Geo",
                                        /* 8*/"Scatter-Desy",  /* 9*/"Scatter-Grey",       /*10*/"Scatter-Inverse grey", /*11*/"Scatter-Jet",
                                        /*12*/"Residuen Glow", /*13*/"Residuen Temp" };


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
    histoThread = nullptr;
    histogramValid = false;
    useResiduenColors = false;
    enaExtract = false;
    enaNoFit   = false;
    firstView  = true;
    extractDrawByMouse = false;
    if ( lastSavedFilePath.isEmpty() ) lastSavedFilePath = dp;
    getConfig();
    ui->setupUi(this);
    setWindowTitle( tit );
    QStringList sl;
    sl << "Key" << "Value";
    ui->tblMetaData->setColumnCount(sl.size());
    ui->tblMetaData->setHorizontalHeaderLabels(sl);
    ui->tblMetaData->hide();
    ui->lblHistogram->hide();
    ui->cbsColorTbl->addItems( _slColornames );
    ui->cbsColorTbl->setCurrentIndex(lastColorTbl);
    on_cbsColorTbl_activated( ui->cbsColorTbl->currentText() );
    ui->lblFilename->setText( "not saved yet" );
    ui->lblFilename->setEnabled(false);
    ui->cbsZoom->setCurrentIndex( zoomFactor );
    ui->cbsRotation->setCurrentIndex( rotIndex );
    ui->togMirrorH->setChecked( swapHor );
    ui->togMirrorV->setChecked( swapVert );
    ui->lblHistogram->setMaximumHeight( 10 + HISTO_SY );
    ui->lblHistogram->setMinimumHeight( 10 + HISTO_SY );
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
    if ( key == "_Calculation_" )
        ui->grpImage->setTitle( "Image "+val );
    else if ( key == "From Image" )
        ui->grpImage->setTitle( val );
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
        fileInfos.centerX = val2.toInt();
    else if ( key == "BeamPosY" )
        fileInfos.centerY = val2.toInt();
    else if ( key.contains("wave",Qt::CaseInsensitive) )
        fileInfos.wavelen = val2.toDouble();
    else if ( key == "EditDet" || key == "SampleDist" )
        fileInfos.distance = val2.toDouble();
    else if ( key == "EditPixelX" || key == "Pixel_X" )
        fileInfos.pixWidthX = val2.toDouble();
    else if ( key == "EditPixelY" || key == "Pixel_Y" )
        fileInfos.pixWidthY = val2.toDouble();
    else if ( key.startsWith("NoFitRect_") )
    {
        QStringList sl = val.simplified().split(" ");
        while ( sl.size() < 4 ) sl << "-1";
        int i = key.rightRef(1).toInt();
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
        on_cbsColorTbl_activated(_slColornames[12]); // tblResGlo
#ifndef CONSOLENPROG
        on_butUpdate_clicked();
#else
        imgNormal = generateImage();
#endif
    }
    //qDebug() << "META" << key << val;
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

/**
 * @brief widImage::setData
 * @param x0 - X min
 * @param x1 - X max
 * @param y0 - Y min
 * @param y1 - Y max
 * @param md - Metadata Hash
 * @param d  - data values
 */
void widImage::setData( int x0, int x1, int y0, int y1, double *d )
{
    //qDebug() << "widImage::setData(" << x0 << x1 << y0 << y1 << d << " )";
    if ( d == nullptr )
    {
#ifndef CONSOLENPROG
        QMessageBox::critical( this, "No data", "No image data given", QMessageBox::Ok );
#endif
        return;
    }

    // Preset internal structure
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
    //qDebug() << len;
    data.reserve( len );
    for ( int i=0; i<len; i++ )
    {
        if ( isnan(*d) )
        {
            d++;
            data.append(0.0);
        }
        else
            data.append(*(d++));
    }
#ifdef CONSOLENPROG
    imgNoGUI = generateImage(); // Zum Abspeichern wichtig...
#else
    histogramValid = false;
    if ( !noGUI )   // setData, ena cbs
    {
        ui->butUpdate->setEnabled(true);
        ui->cbsZoom->setEnabled(true);
    }
    on_butUpdate_clicked(); // Redraw Image
    if ( noGUI ) return;    // setData
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
}

void widImage::saveImage( QString fn )
{
    if ( noGUI )
        imgNoGUI.save(fn);
    else
        imgNormal.save(fn);
}

void widImage::saveImageColor( QString fn, QString coltbl )
{
#ifndef CONSOLENPROG
    QString savColTbl = ui->cbsColorTbl->currentText();
    qDebug() << "saveImageColor" << fn << coltbl << "Old:" << savColTbl;
#endif
    on_cbsColorTbl_activated(coltbl);
    QImage tmp = generateImage();
    if ( ! tmp.save(fn) )
        qDebug() << "saveImageColor: error on save";
#ifndef CONSOLENPROG
    on_cbsColorTbl_activated(savColTbl);
#endif
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


void widImage::on_cbsColorTbl_activated(const QString &arg1)
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

    // aus den Screenshots von Prof. Förster
    //static int rScatter[256] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  11,  11,  11,  29,  29,  45,  45,  62,  62,  62,  78,  78,  95,  95,  95, 111, 111, 127, 127, 127, 144, 144, 160, 160, 160, 177, 177, 193, 193, 210, 210, 210, 226, 226, 244, 244, 244, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
    //static int gScatter[256] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  10,  10,  27,  27,  44,  44,  44,  60,  60,  76,  76,  76,  93,  93, 109, 109, 109, 126, 126, 142, 142, 142, 159, 159, 175, 175, 192, 192, 192, 208, 208, 225, 225, 225, 242, 242, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 251, 251, 234, 234, 234, 217, 217, 201, 201, 185, 185, 185, 168, 168, 151, 151, 151, 135, 135, 118, 118, 118, 102, 102,  85,  85,  85,  69,  69,  52,  52,  36,  36,  36,  19,  19,   3,   3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,   6,   6,  22,  22,  39,  39,  39,  56,  56,  72,  72,  72,  89,  89, 105, 105, 122, 122, 122, 138, 138, 155, 155, 155, 171, 171, 188, 188, 188, 204, 204};
    //static int bScatter[256] = {   0,   0,  34,  34,  34,  51,  51,  67,  67,  67,  84,  84, 100, 100, 100, 117, 117, 133, 133, 133, 150, 150, 166, 166, 183, 183, 183, 199, 199, 216, 216, 216, 233, 233, 249, 249, 249, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 252, 252, 252, 236, 236, 219, 219, 203, 203, 203, 186, 186, 170, 170, 170, 153, 153, 137, 137, 137, 120, 120, 104, 104, 104,  87,  87,  70,  70,  54,  54,  54,  37,  37,  21,  21,  21,   4,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  13,  13,  30,  30,  30,  47,  47,  63,  63,  80,  80,  80,  96,  96, 113, 113, 113, 129, 129, 146, 146, 146, 162, 162, 179, 179, 179, 195, 195, 212, 212, 228, 228, 228, 245, 245, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};


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
    // 12"Residuen Glow", 13"Residuen Temp" };
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
    qApp->setOverrideCursor(Qt::WaitCursor);
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
            ui->lblImage->setMinimumSize( imgNormal.size() );
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
        if ( zoom2pix[ui->cbsZoom->currentIndex()] < 1 )
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
            ui->lblImage->setPixmap( QPixmap::fromImage(imgNormal).scaled(imgNormal.size()*zoom2pix[ui->cbsZoom->currentIndex()]) );
        }
        else
        {
            ui->lblImage->setPixmap( QPixmap::fromImage(imgNormal) );
        }
    }
}

#endif

QImage widImage::generateImage()
{
    // ----- Read the image and displays it
    QImage img( maxX - minX, maxY - minY, QImage::Format_RGB32 );
    //img.fill( Qt::black );
    double vmin=1e26, vmin2=1e26, vmax=0, vlogmin=1e26, vlogmax=0;
#ifndef CONSOLENPROG
    if ( !noGUI && ui->radShowLin->isChecked() )
#else
    if ( useLinOut )
#endif
    {   // Lin
        vlogmin = vlogmax = 0;
        for ( int y=minY; y<maxY; y++ )
            for ( int x=minX; x<maxX; x++ )
                if ( x != 0 && y != 0 )
                {
                    if ( getData(x,y) > 0 && getData(x,y) < vmin2 )
                        vmin2 = getData(x,y);
                    if ( getData(x,y) < vmin )
                        vmin = getData(x,y);
                    else if ( getData(x,y) > vmax )
                        vmax = getData(x,y);
                }
    }
    else
    {   // Log
        for ( int y=minY; y<maxY; y++ )
            for ( int x=minX; x<maxX; x++ )
                if ( ( (minX < 0 && x != 0 && y != 0) ||
                       (minX == 0 && x != fileInfos.centerX && y != fileInfos.centerY) )
                     && getData(x,y) > -1.0e10 )
                {   // Im Mittelpunkt ist i.d.R. ein Extremwert (z.B. -7.4e+25 bei FCC)
                    if ( getData(x,y) > 0 && getData(x,y) < vmin2 )
                        vmin2 = getData(x,y);
                    if ( getData(x,y) < vmin )
                        vmin = getData(x,y);
                    else if ( getData(x,y) > vmax )
                        vmax = getData(x,y);
                    if ( getData(x,y) > 0 )
                    {
                        if ( log10(getData(x,y)) < vlogmin )
                            vlogmin = log10(getData(x,y));
                        else if ( log10(getData(x,y)) > vlogmax )
                            vlogmax = log10(getData(x,y));
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
        if ( ui->radShowLin->isChecked() )
        {
            tmp = input2val( ui->inpFixedMax, &ok ); if ( ok ) vmax = tmp;
            tmp = input2val( ui->inpFixedMin, &ok ); if ( ok ) vmin = tmp;
        }
        else
        {
            tmp = input2val( ui->inpFixedMax, &ok ); if ( ok && tmp > 0 ) vlogmax = log10(tmp);
            tmp = input2val( ui->inpFixedMin, &ok ); if ( ok && tmp > 0 ) vlogmin = log10(tmp);
        }
    }
#endif
    if ( vmin >= vmax ) return QImage();

    if ( farbtabelle.size() < 256 )
#ifndef CONSOLENPROG
        on_cbsColorTbl_activated( noGUI ? _slColornames[1]/*tblGlow*/ : ui->cbsColorTbl->currentText() );
#else
        on_cbsColorTbl_activated( lastColorName );
#endif

#ifndef CONSOLENPROG
    if ( !noGUI && ui->radShowLin->isChecked() )
#else
    if ( useLinOut )
#endif
    {   // Lin
        for ( int y=minY, yp=0; y<maxY; y++, yp++ )
            for ( int x=minX, xp=0; x<maxX; x++, xp++ )
            {
                img.setPixelColor( xp, yp, data2pixelcol(getData(x,y),vmin,vmax,vmax,false) );
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
    /*
    int yp=0, xp=0;
    for ( int c=0; c<farbtabelle.size(); c++ )
    {
        img.setPixelColor( xp, yp, farbtabelle.at(c) );
        if ( ++xp >= img.width() )
        {
            xp = 0;
            yp++;
        }
    }
    */
    /*
    // TEST:
    for ( int y=minY+1, yp=1; y<minY+21; y++, yp++ )
        img.setPixelColor( yp, 1, farbtabelle.last() );
    for ( int x=minX+1, xp=1; x<minX+11; x++, xp++ )
        img.setPixelColor( 1, xp, farbtabelle.last() );
    */

    // Perform Mirroring
    if ( swapHor || swapVert )
        img = img.mirrored( swapHor, swapVert );

    // Perform rotation
    if ( rotIndex == 0 /*0°*/ ) return img;     // to speed up a little bit
    QTransform mat = QTransform(); // Identity
    return img.transformed( mat.rotate( rotIndex * 90 ) );
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
    QString fn = QFileDialog::getSaveFileName( this, "Save", lastSavedFilePath, "Datafiles (*.dat)" );
    if ( fn.isEmpty() ) return;
    if ( !fn.endsWith(".dat",Qt::CaseInsensitive) ) fn += ".dat";
    lastSavedFilePath = fn;
    QFile f(fn);
    ui->lblFilename->setText( fn );
    ui->lblFilename->setEnabled(true);

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
        qDebug() << f.errorString();

    // ----- Save the image in PNG format
    fn = fn.replace(".dat",".png");
    saveImage(fn);

    // ----- Save the image in textual format
    fn = fn.replace(".png",".csv");
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
        qDebug() << f.errorString();

    // ----- Save the MetaData in textual format
    fn = fn.replace(".csv",".txt");
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
        qDebug() << f.errorString();

    // ----- Save NoFitRect
    fileInfos.filePath = lastSavedFilePath;
    saveNoFitRects();
}

void widImage::on_butMetaData_toggled(bool checked)
{
    ui->tblMetaData->setVisible(checked);
    ui->lblHistogram->setVisible(checked);
    ui->butMetaData->setText( checked ? "<<<" : "Meta >>>" );
    if ( checked )
    {
        ui->tblMetaData->resizeColumnsToContents();
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
    //qDebug() << "Histo" << cntmax << max2 << indmax << "/" << facx << facy;
    pnt.setPen( Qt::black );
    for ( int x=0; x<farbtabelle.size()/facx; x++ )
    {
        int y = cnt[x] * facy;
        if ( y < 0 ) { /*qDebug()<<"Histo"<<x<<y<<"<"<<cnt[x];*/ y=0; }
        else if ( y >= HISTO_SY ) { /*qDebug()<<"Histo"<<x<<y<<">"<<cnt[x];*/ y=HISTO_SY-1; }
        else if ( y == 0 && cnt[x] > 0 ) y++;
        //else qDebug()<< "Histo" << x << y << cnt[x];
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
    ui->lblImage->setMinimumSize( imgNormal.size()*zoom2pix[index] );
    if ( enaExtract || enaNoFit )
    {
        extractUpdateImage(false);
    }
    else if ( ! imgNormal.isNull() )
    {
        //qDebug() << "cbsZoom";
        QPixmap pix = QPixmap::fromImage(imgNormal);
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
    pnt.setX( (pos.x() - rpi.left()) / zoom2pix[zoomFactor] );
    pnt.setY( (pos.y() - rpi.top() ) / zoom2pix[zoomFactor] );

    // The image is first mirrored, then rotated. So we must first rotate back then mirror back...
    // The rotation is performed around the center of the image!

    // Perform rotation
    line.setAngle( line.angle() + 90 * rotIndex );

    // Perform Mirroring
    if ( swapHor )
        line.setAngle( 180 - line.angle() );
    if ( swapVert )
        line.setAngle( - line.angle() );

    x = minX + (line.p2().x() - rpi.left()) / zoom2pix[zoomFactor];
    y = minY + (line.p2().y() - rpi.top())  / zoom2pix[zoomFactor];

    return true;
}

void widImage::mousePressEvent(QMouseEvent *event)
{
    if ( !enaExtract && !enaNoFit ) return;
    if ( event->button() != Qt::LeftButton ) return;

    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPos() );
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
    //qDebug() << "Press" << pnt << extractRect;
    extractUpdateImage(false);
}

void widImage::mouseReleaseEvent(QMouseEvent *event)
{
    if ( !enaExtract && !enaNoFit ) return;
    if ( event->button() != Qt::LeftButton || !extractDrawByMouse ) return;

    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPos() );
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
    //qDebug() << "Release" << pnt << extractRect;
    extractUpdateImage(false);
}

void widImage::mouseMoveEvent(QMouseEvent *event)
{
    if ( (event->buttons() & Qt::LeftButton) == 0 ) return;

    QPoint pos = ui->lblImage->mapFromGlobal( event->globalPos() );

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
        //qDebug() << swapVert << swapHor << rotIndex;
        extractRect.setLeft(   qMin(extractMouseStart.x(), pnt.x()) );
        extractRect.setRight(  qMax(extractMouseStart.x(), pnt.x()) );
        extractRect.setTop(    qMin(extractMouseStart.y(), pnt.y()) );
        extractRect.setBottom( qMax(extractMouseStart.y(), pnt.y()) );
        emit extractUpdateRect(extractRect);
        //qDebug() << "move" << pnt << extractRect;
        extractUpdateImage(false);
        return; // Zeichnen des Mausrahmens und der Tooltip geht nicht gut zusammen.
    }

    int idx = XY2IDX(minX,maxX,minY,maxY,x,y);
    if ( idx < 0 || idx >= data.size() ) return;

    //int ix = x - minX;
    //int iy = y - minY;

    double xdet = fileInfos.pixWidthX * (x + fileInfos.centerX);
    double ydet = fileInfos.pixWidthY * (y + fileInfos.centerY);
    double rdet = sqrt(xdet*xdet+ydet*ydet);
    double phidet = atan2(ydet,xdet);
    double thetadet = atan2(rdet,fileInfos.distance);
    double qx = 2*M_PI*cos(phidet)*sin(thetadet)/fileInfos.wavelen;
    double qy = 2*M_PI*sin(phidet)*sin(thetadet)/fileInfos.wavelen;
    double qz = 2*M_PI*(1-cos(thetadet))/fileInfos.wavelen;
    double q = sqrt(qx*qx+qy*qy+qz*qz)+eps9;
/*
    if ( isinf(q) )
        qDebug() << "P" << fileInfos.pixWidthX << fileInfos.pixWidthY
                 << "C" << fileInfos.centerX   << fileInfos.centerY
                 << "D" << fileInfos.distance
                 << "W" << fileInfos.wavelen;
*/
    //if ( ui->radShowLin->isChecked() )
    //{   // Lin
    QToolTip::showText( event->globalPos(),
                        QString("[%1,%2]=%3\nqx=%4,qy=%5,qz=%6)\nq=%7").arg(x).arg(y).arg(data.at(idx)).arg(qx).arg(qy).arg(qz).arg(q) );
                        //QString("[%1,%2]=%3\n%4,%5=%6").arg(x).arg(y).arg(data.at(idx))
                        //.arg(ix).arg(iy).arg(imgNormal.pixel(ix,iy),0,16) );
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


/**
 * @brief widImage::on_lblImage_customContextMenuRequested
 *        Click on right mouse button inside the label (image) window.
 *        This is used to mark the beam center position.
 * @param pos - position of the click in local window coordinates.
 */
void widImage::on_lblImage_customContextMenuRequested(const QPoint &pos)
{
    //qDebug() << "customContextMenuRequested" << pos;
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

    //std::cerr << "IntToColor( " << ywert << ", " << decade << ", " << tbl << " )" << std::endl;
    //IntToColor( 2.57273e+181, 0, 1 )
    //QImage::setPixelColor: color is invalid

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
