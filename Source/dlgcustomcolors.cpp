#include "dlgcustomcolors.h"
#include "ui_dlgcustomcolors.h"
#include <QSettings>
#include <QDebug>


dlgCustomColors::dlgCustomColors(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::dlgCustomColors)
{
    ui->setupUi(this);

    blockUpdate = true;
    colHelper.append( new myColorHelper(BGKeyLattice,ui->lblLattice,ui->togLattice,ui->inpLatticeR,ui->inpLatticeG,ui->inpLatticeB) );
    colHelper.append( new myColorHelper(BGKeyParticle,ui->lblParticle,ui->togParticle,ui->inpParticleR,ui->inpParticleG,ui->inpParticleB) );
    colHelper.append( new myColorHelper(BGKeyPeakShape,ui->lblPeak,ui->togPeak,ui->inpPeakR,ui->inpPeakG,ui->inpPeakB) );
    colHelper.append( new myColorHelper(BGKeyOrientation,ui->lblOrientation,ui->togOrientation,ui->inpOrientationR,ui->inpOrientationG,ui->inpOrientationB) );
    colHelper.append( new myColorHelper(BGKeyExperiment,ui->lblExperiment,ui->togExperiment,ui->inpExperimentR,ui->inpExperimentG,ui->inpExperimentB) );
    colHelper.append( new myColorHelper(BGKeyPixelMan,ui->lblPixel,ui->togPixel,ui->inpPixelR,ui->inpPixelG,ui->inpPixelB) );

    colHelper.append( new myColorHelper(BGKeyDisText,ui->lblDisText,nullptr,ui->inpDisTextR,ui->inpDisTextG,ui->inpDisTextB) );
    colHelper.append( new myColorHelper(BGKeyOnly1D,ui->lblOnly1D,nullptr,ui->inpOnly1DR,ui->inpOnly1DG,ui->inpOnly1DB) );
    colHelper.append( new myColorHelper(BGKEYOnly2D,ui->lblOnly2D,nullptr,ui->inpOnly2DR,ui->inpOnly2DG,ui->inpOnly2DB) );
    colHelper.append( new myColorHelper(BGKey1Dand2D,ui->lblBothD,nullptr,ui->inpBothDR,ui->inpBothDG,ui->inpBothDB) );

    blockUpdate = false;

    adjustSize();;
}

dlgCustomColors::~dlgCustomColors()
{
    delete ui;
}

void dlgCustomColors::closeEvent(QCloseEvent *event)
{
    //qDebug() << "dlgCustomColors:closeEvent";
    blockUpdate = true;
    foreach (myColorHelper *ch, colHelper)
        delete ch;
    colHelper.clear();
    event->accept();
}

void dlgCustomColors::updateColor(QString key)
{
    if ( blockUpdate ) return;
    foreach (myColorHelper *ch, colHelper)
    {
        if ( ch->key == key )
        {
            ch->update();
            return;
        }
    }
}

void dlgCustomColors::on_butExit_clicked()
{
    //qDebug() << "dlgCustomColors:on_butExit_clicked";
    accept();
}

void dlgCustomColors::on_butSave_clicked()
{
    //qDebug() << "dlgCustomColors:on_butSave_clicked";
    foreach (myColorHelper *ch, colHelper)
        ch->save();
    accept();
}




myColorHelper::myColorHelper( QString k, QLabel *l, QCheckBox *t, QSpinBox *r, QSpinBox *g, QSpinBox *b )
{
    key=k;
    lbl=l;
    tog=t;
    red=r;
    grn=g;
    blu=b;
    //qDebug() << "myColorHelper::create()" << lbl->text() << k;

    // Get the current settings from the config file
    QSettings data(SETT_APP,SETT_PAR);
    if ( tog != nullptr )
        tog->setChecked( data.value( key+"_tog", false ).toBool() );
    red->setValue( data.value( key+"_red", 100 ).toInt() );
    grn->setValue( data.value( key+"_grn", 200 ).toInt() );
    blu->setValue( data.value( key+"_blu", 100 ).toInt() );
    update();
}

void myColorHelper::save()
{
    QSettings data(SETT_APP,SETT_PAR);
    if ( tog != nullptr )
        data.setValue( key+"_tog", tog->isChecked() );
    data.setValue( key+"_red", red->value() );
    data.setValue( key+"_grn", grn->value() );
    data.setValue( key+"_blu", blu->value() );
    //qDebug() << "myColorHelper::save" << key << SETT_APP << SETT_PAR
    //         << red->value() << grn->value() << blu->value();
}

void myColorHelper::update()
{
    //qDebug() << "myColorHelper::update()" << key;
    if ( tog == nullptr || tog->isChecked() )
    {
        QPalette pal = lbl->palette();
        QColor col = QColor(red->value(),grn->value(),blu->value());
        pal.setBrush( QPalette::Window, col );
        //qDebug() << lbl->text() << col << qGray(col.rgb());
        if ( qGray(col.rgb()) < 150 )
            pal.setColor( QPalette::WindowText, Qt::white );
        else
            pal.setColor( QPalette::WindowText, Qt::black );
        lbl->setPalette(pal);
        lbl->setAutoFillBackground(true);
    }
    else
        lbl->setAutoFillBackground(false);
}
