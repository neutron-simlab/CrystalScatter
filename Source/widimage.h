#ifndef WIDIMAGE_H
#define WIDIMAGE_H

// CONSOLENPROG wird definiert, wenn diese Klasse beim Consolen-Programm ohne GUI verwendet wird

#ifndef CONSOLENPROG
#include <QWidget>
#include <QCloseEvent>
#include <QLineEdit>
#ifdef IMG_ONE_WINDOW
#include <QMdiSubWindow>
#endif // IMG_ONE_WINDOW
#endif // CONSOLENPROG
#include <QVector>
#include <QColor>
#include <QSettings>
#include <QHash>
#include <QDebug>
#include <QImage>
#include <QThread>
#include "sc_math.h"


#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))

#ifndef CONSOLENPROG
#define tblResGlo  widImage::slColorNames()[12]
#define tblResTemp widImage::slColorNames()[13]
#endif



// Qt-Versionen:
//  Büro   qmake     5.9.7  (Anaconda)
//         QtCreator 5.15.2 (/usr/bin/qmake-qt5)
//  Laptop qmake     <nicht im Pfad>
//         QtCreator 5.14.2

//#if (QT_VERSION <= QT_VERSION_CHECK(5, 12, 11))
//  ...
//#else
//  ...
//#endif


#ifndef CONSOLENPROG
namespace Ui {
class widImage;
}


class widImage;
class calcHistorgramm : public QThread
{
    Q_OBJECT
public:
    calcHistorgramm( widImage *i );
    void run() override;
signals:
    void setHistoPixmap( QPixmap );
private:
    widImage *img;
};



class widImage : public QWidget
{
    Q_OBJECT

public:
    widImage();
    widImage( QString tit, QString dp, QWidget *parent = nullptr);
    ~widImage();
    void getConfig();

    // Funktionen für die Extract-Funktion:
    void enableExtract( bool f ) { enaExtract=f; extractUpdateImage(true); }
    void setExtractRect( QRect rc ) { extractRect=rc; extractUpdateImage(false); }
    void setExtractFrameCol( int id ) {  extractFrameCol=id; extractUpdateImage(false); }

    void enableNoFit( bool f ) { enaNoFit=f; extractUpdateImage(true); }
    void setNoFitRect( int r, QRect rc ) { noFitRect[r]=rc; extractUpdateImage(false); }
    QRect getNoFitRect( int r ) { return noFitRect[r]; }
    void saveNoFitRects();

#else

class widImage
{
public:
    widImage();
    void setConfig(QString col, bool sh, bool sv, int r, int z, bool lin);
#endif
    static QStringList slColorNames() { return _slColornames; }

    static void setDefColTbl( int idx ) { lastColorTbl=idx; }

    void saveParams( QSettings &sets );
    void loadParams( QSettings &sets );

    void addMetaInfo( QString key, QString val );
    QString getMetaTitle();
    QString metaInfo( QString key, QString def="" );

    void getVarScaling( double &min, double &max );
    void setFixScaling( double min, double max );
    void setLogScaling( bool on );

    void setData(int x0, int x1, int y0, int y1, /*_metaData &md,*/ double *d);
    const double *dataPtr() { return data.constData(); }
    double xmin() { return minX; }
    double xmax() { return maxX; }
    double ymin() { return minY; }
    double ymax() { return maxY; }
    int myWidth() { return maxX - minX; }
    int myHeight() { return maxY - minY; }
    inline double getData(int x,int y) { return data[(tmpidx=XY2IDX(minX,maxX,minY,maxY,x,y))]; }

    QString debugInfo() { return QString("X:%1 .. %2 / Y:%3 .. %4").arg(minX).arg(maxX).arg(minY).arg(maxY); }

    typedef struct
    {
        double  wavelen;
        double  distance;
        double  pixWidthX, pixWidthY;
        int     pixelX, pixelY;
        int     centerX, centerY;
        QString filePath;
    } _fileInfos;
    inline _fileInfos *getFileInfos() { return &fileInfos; }

    void saveImage( QString fn );
//#ifndef CONSOLENPROG
    void saveImageColor( QString fn, QString coltbl );
//#endif
    void saveImageGray( QString fn );
    void saveImageBinary( QString fn );
    void saveImage( QString fn, QSize siz );
    void saveImageGray( QString fn, QSize siz );
    void saveImageBinary( QString fn, QSize siz );

    void calculateHistogram( QPixmap coltbl, QImage src );

    QPixmap getColorTable();
    QImage  getCurrentImage();
    QVector<QColor> getFarbtabelle() { return farbtabelle; }
    QColor  IntToColor( double ywert, double decade, int tbl ); // aus Pascal-programm

#ifndef CONSOLENPROG
protected:
    void closeEvent(QCloseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

private slots:
    void on_cbsColorTbl_activated(const QString &arg1);
    void on_butUpdate_clicked();
    void on_butSave_clicked();
    void on_butMetaData_toggled(bool checked);
    void on_cbsZoom_activated(int index);
    void on_cbsRotation_activated(int index);
    void on_togMirrorH_toggled(bool checked);
    void on_togMirrorV_toggled(bool checked);
    void on_lblImage_customContextMenuRequested(const QPoint &pos);

    void setHistogram( QPixmap dst );

private:
    Ui::widImage *ui;
    bool mouse2data(QPoint pos, int &x, int &y, QPoint &pnt);
#else
private:
    void on_cbsColorTbl_activated(const QString &arg1); // wird für die Farbtabelle gebraucht
    inline void on_butUpdate_clicked() { imgNormal = generateImage(); }
    QString lastColorName;
#endif

private:
    bool noGUI;
    QImage imgNoGUI;
    QImage imgNormal;

    QImage generateImage();

    static QString lastSavedFilePath;
    static int lastColorTbl;

    bool swapHor, swapVert;
    int  rotIndex;
    int  zoomFactor;
    bool useLinOut;

    bool firstView; // zum Sperren weiterer adjustSize() bei setData(...)
    // TODO: Wenn sich Bildgröße oder Zoom geändert hat, dann einmal adjustSize() wieder freigeben.

    _fileInfos fileInfos;
    bool bRekursiv;

    int minX, maxX, minY, maxY;
    QVector<double> data;
    int tmpidx;

    static QStringList _slColornames;
    int curColorSelected;
    QVector<QColor> farbtabelle;
    void updateColortable( int *r, int *g, int *b );
    bool useResiduenColors; // value<0 --> greyscale / value>=0 --> Glow or Temp
    QColor data2pixelcol( double val, double min, double max, double datamax, bool islog );

    QRect  noFitRect[4];

#ifndef CONSOLENPROG
    double input2val(QLineEdit *inp, bool *ok);

    calcHistorgramm *histoThread;
    bool histogramValid;

    bool   enaExtract;
    QRect  extractRect;
    int    extractFrameCol;
    QColor extractFrameColor();
    void   extractUpdateImage(bool flg);
    bool   extractDrawByMouse;
    QPoint extractMouseStart;
    bool   enaNoFit;
signals:
    void extractUpdateRect(QRect);

#endif
};

#endif // WIDIMAGE_H
