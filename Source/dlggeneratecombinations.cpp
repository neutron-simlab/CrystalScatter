#include "dlggeneratecombinations.h"
#include "ui_dlggeneratecombinations.h"
#include "sc_globalConfig.h"
#include <QCheckBox>
#include <QSettings>
#include <QFileDialog>
#include <QDebug>

#define SETT_GENCOMB "GenComb"



dlgGenerateCombinations::dlgGenerateCombinations(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::dlgGenerateCombinations)
{
    ui->setupUi(this);
    //ui->togDoCalc->hide(); // Erstmal nur die INI-Dateien erzeugen
    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup(SETT_GENCOMB);
    ui->inpOutPath->setText(sets.value("OutPath",".").toString());
    ui->togDoCalc->setChecked(sets.value("DoCalcImg",false).toBool());
    ui->togCheckNaN->setChecked(sets.value("DoCheckNaN",false).toBool());
    ui->togCheckNaN->setEnabled( ui->togDoCalc->isChecked() );
#ifndef ChatbotDisabled
    ui->togGenJson->setChecked(sets.value("DoCheckJson",false).toBool());
    ui->togGenJson->setEnabled( ui->togDoCalc->isChecked() );
#else
    ui->togGenJson->hide();
#endif
    ui->togCalcAll->setChecked(sets.value("DoCalcAll",false).toBool());
    // An dieser Stelle wird die Zeitbegrenzung nicht aus dem Parameterfile geholt,
    //  sie wird direkt im Hauptprogramm aus der globalen Konfiguration genommen.
}

dlgGenerateCombinations::~dlgGenerateCombinations()
{
    delete ui;
}

void dlgGenerateCombinations::on_butDone_clicked()
{
    reject();
}

int dlgGenerateCombinations::setValues( QString hdr, QStringList val )
{
    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup(SETT_GENCOMB+hdr);
    int col = ui->tblValues->columnCount();
    ui->tblValues->setColumnCount(col+1);
    ui->tblValues->setHorizontalHeaderItem( col, new QTableWidgetItem(hdr) );
    if ( ui->tblValues->rowCount() < val.size() )
        ui->tblValues->setRowCount(val.size());
    for ( int r=0; r<val.size(); r++ )
    {
        QCheckBox *cb = new QCheckBox(val.at(r));
        cb->setChecked(sets.value(val.at(r),true).toBool());
        ui->tblValues->setCellWidget( r, col, cb );
    }
    ui->tblValues->resizeColumnToContents(col);
    return ui->tblValues->columnWidth(col);
}

void dlgGenerateCombinations::on_butStartGen_clicked()
{
    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup(SETT_GENCOMB);
    sets.setValue("OutPath",ui->inpOutPath->text());
    sets.setValue("DoCalcImg",ui->togDoCalc->isChecked());
    sets.setValue("DoCheckNaN",ui->togCheckNaN->isChecked());
#ifndef ChatbotDisabled
    sets.setValue("DoCheckJson",ui->togGenJson->isChecked());
#endif
    sets.setValue("DoCalcAll",ui->togCalcAll->isChecked());
    accept();
}

QStringList dlgGenerateCombinations::getValues( QString hdr )
{
    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup(SETT_GENCOMB+hdr);
    QStringList rv;
    for ( int c=0; c<ui->tblValues->columnCount(); c++ )
        if ( ui->tblValues->horizontalHeaderItem(c)->text() == hdr )
        {
            for ( int r=0; r<ui->tblValues->rowCount(); r++ )
            {
                QCheckBox *cb = static_cast<QCheckBox*>(ui->tblValues->cellWidget(r,c));
                if ( cb == nullptr ) break;
                sets.setValue(cb->text(),cb->isChecked());
                if ( cb->isChecked() ) rv << cb->text();
            }
        }
    return rv;
}

QList<int> dlgGenerateCombinations::getValueIDs( QString hdr )
{
    QSettings sets(SETT_APP,SETT_GUI);
    sets.beginGroup(SETT_GENCOMB+hdr);
    QList<int> rv;
    for ( int c=0; c<ui->tblValues->columnCount(); c++ )
        if ( ui->tblValues->horizontalHeaderItem(c)->text() == hdr )
        {
            for ( int r=0; r<ui->tblValues->rowCount(); r++ )
            {
                QCheckBox *cb = static_cast<QCheckBox*>(ui->tblValues->cellWidget(r,c));
                if ( cb == nullptr ) break;
                sets.setValue(cb->text(),cb->isChecked());
                if ( cb->isChecked() ) rv << r;
            }
        }
    return rv;
}

void dlgGenerateCombinations::setCalcTimeLimit( int secs )
{
    if ( secs <= 0 )
    {
        ui->togLimitCalc->setChecked(false);
        ui->inpLimitCalc->setEnabled(false);
    }
    else
    {
        ui->togLimitCalc->setChecked(true);
        ui->inpLimitCalc->setEnabled(true);
        ui->inpLimitCalc->setValue(secs);
    }
}

int dlgGenerateCombinations::getCalcTimeLimit()
{
    if ( ui->togLimitCalc->isChecked() )
        return ui->inpLimitCalc->value();
    return 0;
}

void dlgGenerateCombinations::on_butOutPath_clicked()
{
    QString fn = QFileDialog::getExistingDirectory( this, "Directory to save Combinations", ui->inpOutPath->text() );
    if ( fn.isEmpty() ) return;
    ui->inpOutPath->setText(fn);
}

QString dlgGenerateCombinations::getPath()
{
    return ui->inpOutPath->text();
}

bool dlgGenerateCombinations::getCalcDirect()
{
    return ui->togDoCalc->isChecked();
}

bool dlgGenerateCombinations::getCalcAll()
{
    return ui->togCalcAll->isChecked();
}

bool dlgGenerateCombinations::getCalcNaNCheck()
{
    return ui->togCheckNaN->isChecked();
}

bool dlgGenerateCombinations::getCalcGenJson()
{
#ifndef ChatbotDisabled
    return ui->togGenJson->isChecked();
#else
    return false;
#endif
}

void dlgGenerateCombinations::on_butSaveSettings_clicked()
{
    QSettings gsets(SETT_APP,SETT_GUI);
    gsets.beginGroup(SETT_GENCOMB);
    QString fn = gsets.value("SettingsFile",ui->inpOutPath->text()+"/GenerateCombinations.ini").toString();
    fn = QFileDialog::getSaveFileName( this, "Save settings for console prog", fn, "Setting files (*.ini)" );
    if ( fn.isEmpty() ) return;
    gsets.setValue("SettingsFile",fn);

    QSettings sets(fn, QSettings::IniFormat);
    // Globale Flags in Generic
    sets.setValue("OutPath",ui->inpOutPath->text());
    sets.setValue("DoCalcImg",ui->togDoCalc->isChecked());
    sets.setValue("DoCheckNaN",ui->togCheckNaN->isChecked());
#ifndef ChatbotDisabled
    sets.setValue("DoCheckJson",ui->togGenJson->isChecked());
#endif
    sets.setValue("DoCalcAll",ui->togCalcAll->isChecked());
    sets.setValue("CalcLimitFlag",ui->togLimitCalc->isChecked());
    sets.setValue("CalcLimitVal",ui->inpLimitCalc->value());

    // Einstellungen in Gruppen
    for ( int c=0; c<ui->tblValues->columnCount(); c++ )
    {
        QString hdr = ui->tblValues->horizontalHeaderItem(c)->text();
        sets.beginGroup(hdr);
        for ( int r=0; r<ui->tblValues->rowCount(); r++ )
        {
            QCheckBox *cb = static_cast<QCheckBox*>(ui->tblValues->cellWidget(r,c));
            if ( cb == nullptr ) break;
            sets.setValue(cb->text(),cb->isChecked());
        }
        sets.endGroup();
    }
    qDebug() << fn;
}

void dlgGenerateCombinations::on_butLoadSettings_clicked()
{
    QSettings gsets(SETT_APP,SETT_GUI);
    gsets.beginGroup(SETT_GENCOMB);
    QString fn = gsets.value("SettingsFile",ui->inpOutPath->text()+"/GenerateCombinations.ini").toString();
    fn = QFileDialog::getOpenFileName( this, "Load settings for console prog", fn, "Setting files (*.ini)" );
    if ( fn.isEmpty() ) return;
    gsets.setValue("SettingsFile",fn);

    QSettings sets(fn, QSettings::IniFormat);
    // Globale Flags in Generic
    ui->inpOutPath->setText( sets.value("OutPath",".").toString() );
    ui->togDoCalc->setChecked( sets.value("DoCalcImg",false).toBool() );
    ui->togCheckNaN->setChecked( sets.value("DoCheckNaN",false).toBool() );
#ifndef ChatbotDisabled
    ui->togGenJson->setChecked( sets.value("DoCheckJson",false).toBool() );
#endif
    ui->togCalcAll->setChecked( sets.value("DoCalcAll",false).toBool() );
    ui->togLimitCalc->setChecked( sets.value("CalcLimitFlag",true).toBool() );
    ui->inpLimitCalc->setValue( sets.value("CalcLimitVal",5).toInt() );

    // Einstellungen in Gruppen
    for ( int c=0; c<ui->tblValues->columnCount(); c++ )
    {
        QString hdr = ui->tblValues->horizontalHeaderItem(c)->text();
        sets.beginGroup(hdr);
        for ( int r=0; r<ui->tblValues->rowCount(); r++ )
        {
            QCheckBox *cb = static_cast<QCheckBox*>(ui->tblValues->cellWidget(r,c));
            if ( cb == nullptr ) break;
            cb->setChecked( sets.value(cb->text(),false).toBool() );
        }
        sets.endGroup();
    }
}

