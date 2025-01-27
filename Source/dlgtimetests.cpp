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
    ui->tog1d->setChecked(sets.value("use1d",false).toBool());
    ui->radSwitchBoth->setChecked(sets.value("newswitch",swBoth).toInt() == swBoth);
    ui->radSwitchNew->setChecked(sets.value("newswitch",swBoth).toInt() == swNew);
    ui->radSwitchOld->setChecked(sets.value("newswitch",swBoth).toInt() == swOld);
    ui->togSwitchModFN->setChecked(sets.value("switchmodfn",true).toBool());
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
    if ( ui->grpHKLmax->isEnabled() )
    {
        if ( ui->togHKL2->isChecked() ) cntHKL++;
        if ( ui->togHKL3->isChecked() ) cntHKL++;
        if ( ui->togHKL4->isChecked() ) cntHKL++;
        if ( ui->togHKL5->isChecked() ) cntHKL++;
    }
    else
        cntHKL = 1;
    if ( ui->tog1d->isChecked() ) cntQ++;
    if ( ui->togQ1->isChecked() ) cntQ++;
    if ( ui->togQ2->isChecked() ) cntQ++;
    if ( ui->togQ4->isChecked() ) cntQ++;
    if ( ui->radSwitchBoth->isChecked() ) cntThreads *= 2;
    ui->lblNumCalc->setNum( cntThreads * cntHKL * cntQ );
    ui->butStart->setEnabled( cntThreads * cntHKL * cntQ > 0 );
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
    sets.setValue("use1d",ui->tog1d->isChecked());
    sets.setValue("enaupdates",ui->togEnaUpdates->isChecked());
    sets.setValue("saveimages",ui->togSaveImages->isChecked());
    sets.setValue("savefile",ui->togSaveFile->isChecked());
    sets.setValue("loops",ui->inpNumLoops->value());
    sets.setValue("filename",ui->inpSaveFilename->text());
    sets.setValue("comment",ui->inpComment->text());
    if ( ui->radSwitchBoth->isChecked() ) sets.setValue("newswitch",swBoth);
    if ( ui->radSwitchNew ->isChecked() ) sets.setValue("newswitch",swNew);
    if ( ui->radSwitchOld ->isChecked() ) sets.setValue("newswitch",swOld);
    sets.setValue("switchmodfn",ui->togSwitchModFN->isChecked());
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

void dlgTimeTests::setHKLmaxUsed(bool f)
{
    ui->grpHKLmax->setEnabled(f);
    onTogToggled();
}

QVector<int> dlgTimeTests::getHKLmax()
{
    QVector<int> rv;
    if ( ui->grpHKLmax->isEnabled() )
    {
        if ( ui->togHKL2->isChecked() ) rv.append( ui->inpHKL2->value() );
        if ( ui->togHKL3->isChecked() ) rv.append( ui->inpHKL3->value() );
        if ( ui->togHKL4->isChecked() ) rv.append( ui->inpHKL4->value() );
        if ( ui->togHKL5->isChecked() ) rv.append( ui->inpHKL5->value() );
    }
    else
        rv.append(1);
    return rv;
}

QVector<int> dlgTimeTests::getQuadrants()
{
    QVector<int> rv;
    if ( ui->tog1d->isChecked() ) rv.append( 0 );  // Kennung für 1D
    if ( ui->togQ1->isChecked() ) rv.append( 1 );
    if ( ui->togQ2->isChecked() ) rv.append( 2 );
    if ( ui->togQ4->isChecked() ) rv.append( 4 );
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

dlgTimeTests::_newSwitch dlgTimeTests::getNewSwitchFlag()
{
    // typedef enum { swBoth, swNew, swOld } _newSwitch;
    if ( ui->radSwitchBoth->isChecked() ) return swBoth;
    if ( ui->radSwitchNew->isChecked() ) return swNew;
    if ( ui->radSwitchOld->isChecked() ) return swOld;
    return swBoth;
}

bool dlgTimeTests::getNewSwitchModFN()
{
    return ui->togSwitchModFN->isChecked();
}

QString dlgTimeTests::newSwitch2Str(_newSwitch nsw)
{
    switch ( nsw )
    {
    case swBoth: return "both";
    case swNew:  return "new";
    case swOld:  return "old";
    }
    return "both";
}

void dlgTimeTests::on_butSaveTests_clicked()
{
    QSettings sets(SETT_APP,SETT_PAR);
    QString fn = sets.value("saveparams",".").toString();
    fn = QFileDialog::getSaveFileName(this,"Save Test Params",fn,"Parameterfiles (*.ini)");
    if ( fn.isEmpty() ) return;
    sets.setValue("saveparams",fn);
    sets.sync();

    QSettings set2(fn,QSettings::IniFormat);
    set2.setValue("thread0",ui->togThread0->isChecked());
    set2.setValue("thread1",ui->togThread1->isChecked());
    set2.setValue("thread4",ui->togThread4->isChecked());
    set2.setValue("threadM",ui->togThreadM->isChecked());
    set2.setValue("hkl2",ui->togHKL2->isChecked());
    set2.setValue("hkl3",ui->togHKL3->isChecked());
    set2.setValue("hkl4",ui->togHKL4->isChecked());
    set2.setValue("hkl5",ui->togHKL5->isChecked());
    set2.setValue("hkl2val",ui->inpHKL2->value());
    set2.setValue("hkl3val",ui->inpHKL3->value());
    set2.setValue("hkl4val",ui->inpHKL4->value());
    set2.setValue("hkl5val",ui->inpHKL5->value());
    set2.setValue("useq1",ui->togQ1->isChecked());
    set2.setValue("useq2",ui->togQ2->isChecked());
    set2.setValue("useq4",ui->togQ4->isChecked());
    set2.setValue("use1d",ui->tog1d->isChecked());
    set2.setValue("enaupdates",ui->togEnaUpdates->isChecked());
    set2.setValue("saveimages",ui->togSaveImages->isChecked());
    set2.setValue("savefile",ui->togSaveFile->isChecked());
    set2.setValue("loops",ui->inpNumLoops->value());
    set2.setValue("filename",ui->inpSaveFilename->text());
    set2.setValue("comment",ui->inpComment->text());
    if ( ui->radSwitchBoth->isChecked() ) set2.setValue("newswitch",swBoth);
    if ( ui->radSwitchNew ->isChecked() ) set2.setValue("newswitch",swNew);
    if ( ui->radSwitchOld ->isChecked() ) set2.setValue("newswitch",swOld);
    set2.setValue("switchmodfn",ui->togSwitchModFN->isChecked());
    set2.sync();
}

#ifdef undef
void dlgTimeTests::on_togQ1_toggled(bool checked)
{
    if ( checked )
    {
        ui->tog1d->blockSignals(true); ui->tog1d->setChecked(false); ui->tog1d->blockSignals(false);
    }
    ui->tog1d->setDisabled(checked);
    onTogToggled(); // Update counter
}
void dlgTimeTests::on_togQ2_toggled(bool checked)
{
    if ( checked )
    {
        ui->tog1d->blockSignals(true); ui->tog1d->setChecked(false); ui->tog1d->blockSignals(false);
    }
    ui->tog1d->setDisabled(checked);
    onTogToggled(); // Update counter
}
void dlgTimeTests::on_togQ4_toggled(bool checked)
{
    if ( checked )
    {
        ui->tog1d->blockSignals(true); ui->tog1d->setChecked(false); ui->tog1d->blockSignals(false);
    }
    ui->tog1d->setDisabled(checked);
    onTogToggled(); // Update counter
}
void dlgTimeTests::on_tog1d_toggled(bool checked)
{
    if ( checked )
    {
        ui->togQ1->blockSignals(true); ui->togQ1->setChecked(false); ui->togQ1->blockSignals(false);
        ui->togQ2->blockSignals(true); ui->togQ2->setChecked(false); ui->togQ2->blockSignals(false);
        ui->togQ4->blockSignals(true); ui->togQ4->setChecked(false); ui->togQ4->blockSignals(false);
    }
    ui->togQ1->setDisabled(checked);
    ui->togQ2->setDisabled(checked);
    ui->togQ4->setDisabled(checked);
    onTogToggled(); // Update counter
}
#endif
