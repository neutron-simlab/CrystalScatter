#include "dlgtimetests.h"
#include "ui_dlgtimetests.h"
#include <QSettings>
#include <QFileDialog>

#define SETT_APP   "JCNS-1-SasCrystal"  // aus sc_maingui.h
#define SETT_PAR   "TimeTests"


dlgTimeTests::dlgTimeTests( QString path, bool gpu, int maxt, QWidget *parent ) :
    QDialog(parent),
    ui(new Ui::dlgTimeTests)
{
    ui->setupUi(this);
    dataPath = path;
    QSettings sets(SETT_APP,SETT_PAR);
    ui->togThread0->setChecked(sets.value("thread0",true).toBool() && gpu);
    ui->togThread1->setChecked(sets.value("thread1",true).toBool());
    ui->togThread4->setChecked(sets.value("thread4",true).toBool());
    ui->togThreadM->setChecked(sets.value("threadM",true).toBool() && (maxt>4));
    ui->togHKL2->setChecked(sets.value("hkl2",false).toBool());
    ui->togHKL3->setChecked(sets.value("hkl3",true ).toBool());
    ui->togHKL4->setChecked(sets.value("hkl4",false).toBool());
    ui->togHKL5->setChecked(sets.value("hkl5",true ).toBool());
    ui->inpHKL2->setValue(sets.value("hkl2val",2).toInt());
    ui->inpHKL3->setValue(sets.value("hkl3val",3).toInt());
    ui->inpHKL4->setValue(sets.value("hkl4val",4).toInt());
    ui->inpHKL5->setValue(sets.value("hkl5val",5).toInt());
    if ( maxt <= 4 ) ui->togThreadM->setEnabled(false);
    if ( maxt <  4 ) ui->togThread4->setEnabled(false);
    if ( ! gpu     ) ui->togThread0->setEnabled(false);
    ui->togQ1->setChecked(sets.value("useq1",false).toBool());
    ui->togQ2->setChecked(sets.value("useq2",false).toBool());
    ui->togQ4->setChecked(sets.value("useq4",true ).toBool());
    ui->grpQuadrants->setChecked(sets.value("usequadrants",true).toBool());
    onTogToggled();
    ui->togEnaUpdates->setChecked(sets.value("enaupdates",false).toBool());
    ui->togSaveImages->setChecked(sets.value("saveimages",false).toBool());
    ui->togSaveFile->setChecked(sets.value("savefile",false).toBool());
    on_togSaveFile_toggled(ui->togSaveFile->isChecked()); // Sicherheitshalber, falls das set keine Änderung macht
    ui->inpNumLoops->setValue(sets.value("loops",10).toInt());
    ui->inpSaveFilename->setText(sets.value("filename","").toString());
    ui->inpComment->setText(sets.value("comment","").toString());
}

dlgTimeTests::~dlgTimeTests()
{
    delete ui;
}

void dlgTimeTests::onTogToggled()
{
    int cntThreads, cntHKL, cntQ;
    cntThreads = 0;
    cntHKL = 0;
    cntQ = 0;
    if ( ui->togThread0->isChecked() ) cntThreads++;
    if ( ui->togThread1->isChecked() ) cntThreads++;
    if ( ui->togThread4->isChecked() ) cntThreads++;
    if ( ui->togThreadM->isChecked() ) cntThreads++;
    if ( ui->togHKL2->isChecked() ) cntHKL++;
    if ( ui->togHKL3->isChecked() ) cntHKL++;
    if ( ui->togHKL4->isChecked() ) cntHKL++;
    if ( ui->togHKL5->isChecked() ) cntHKL++;
    if ( ui->grpQuadrants->isChecked() )
    {
        if ( ui->togQ1->isChecked() ) cntQ++;
        if ( ui->togQ2->isChecked() ) cntQ++;
        if ( ui->togQ4->isChecked() ) cntQ++;
    }
    else
        cntQ = 1;
    ui->lblNumCalc->setNum( cntThreads * cntHKL * cntQ );
}

void dlgTimeTests::on_butStart_clicked()
{
    QSettings sets(SETT_APP,SETT_PAR);
    sets.setValue("thread0",ui->togThread0->isChecked());
    sets.setValue("thread1",ui->togThread1->isChecked());
    sets.setValue("thread4",ui->togThread4->isChecked());
    sets.setValue("threadM",ui->togThreadM->isChecked());
    sets.setValue("hkl2",ui->togHKL2->isChecked());
    sets.setValue("hkl3",ui->togHKL3->isChecked());
    sets.setValue("hkl4",ui->togHKL4->isChecked());
    sets.setValue("hkl5",ui->togHKL5->isChecked());
    sets.setValue("hkl2val",ui->inpHKL2->value());
    sets.setValue("hkl3val",ui->inpHKL3->value());
    sets.setValue("hkl4val",ui->inpHKL4->value());
    sets.setValue("hkl5val",ui->inpHKL5->value());
    sets.setValue("useq1",ui->togQ1->isChecked());
    sets.setValue("useq2",ui->togQ2->isChecked());
    sets.setValue("useq4",ui->togQ4->isChecked());
    sets.setValue("usequadrants",ui->grpQuadrants->isChecked());
    sets.setValue("enaupdates",ui->togEnaUpdates->isChecked());
    sets.setValue("saveimages",ui->togSaveImages->isChecked());
    sets.setValue("savefile",ui->togSaveFile->isChecked());
    sets.setValue("loops",ui->inpNumLoops->value());
    sets.setValue("filename",ui->inpSaveFilename->text());
    sets.setValue("comment",ui->inpComment->text());
    accept();
}

void dlgTimeTests::on_butCancel_clicked()
{
    reject();
}

QVector<int> dlgTimeTests::getThreads()
{
    QVector<int> rv;
    if ( ui->togThread0->isChecked() ) rv.append( 0 );
    if ( ui->togThread1->isChecked() ) rv.append( 1 );
    if ( ui->togThread4->isChecked() ) rv.append( 4 );
    if ( ui->togThreadM->isChecked() ) rv.append( 100 );
    return rv;
}

QVector<int> dlgTimeTests::getHKLmax()
{
    QVector<int> rv;
    if ( ui->togHKL2->isChecked() ) rv.append( ui->inpHKL2->value() );
    if ( ui->togHKL3->isChecked() ) rv.append( ui->inpHKL3->value() );
    if ( ui->togHKL4->isChecked() ) rv.append( ui->inpHKL4->value() );
    if ( ui->togHKL5->isChecked() ) rv.append( ui->inpHKL5->value() );
    return rv;
}

QVector<int> dlgTimeTests::getQuadrants()
{
    QVector<int> rv;
    if ( ui->grpQuadrants->isChecked() )
    {
        if ( ui->togQ1->isChecked() ) rv.append( 1 );
        if ( ui->togQ2->isChecked() ) rv.append( 2 );
        if ( ui->togQ4->isChecked() ) rv.append( 4 );
    }
    else
        rv.append( 0 ); // Kennung für keine Änderung
    return rv;
}

int dlgTimeTests::getRepetitions()
{
    return ui->inpNumLoops->value();
}

QString dlgTimeTests::getSaveFilename()
{
    if ( ui->togSaveFile->isChecked() )
        return ui->inpSaveFilename->text();
    return "";
}

QString dlgTimeTests::getComment()
{
    if ( ui->togSaveFile->isChecked() )
        return ui->inpComment->text();
    return "";
}

bool dlgTimeTests::getEnaUpdates()
{
    return ui->togEnaUpdates->isChecked();
}

bool dlgTimeTests::getEnaSaveImages()
{
    return ui->togSaveImages->isChecked() && ui->togSaveFile->isChecked() && !ui->inpSaveFilename->text().isEmpty();
}

void dlgTimeTests::on_togSaveFile_toggled(bool checked)
{
    ui->inpSaveFilename->setEnabled(checked);
    ui->togSaveImages->setEnabled(checked);
    ui->inpComment->setEnabled(checked);
    ui->lblComment->setEnabled(checked);
}

void dlgTimeTests::on_butSaveFilename_clicked()
{
    QString fn = ui->inpSaveFilename->text();
    if ( fn.isEmpty() ) fn = dataPath;
    fn = QFileDialog::getSaveFileName( this, "Select file to save timings", fn );
    if ( fn.isEmpty() ) return;
    ui->inpSaveFilename->setText(fn);
}
