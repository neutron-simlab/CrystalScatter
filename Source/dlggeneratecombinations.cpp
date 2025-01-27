#include "dlggeneratecombinations.h"
#include "ui_dlggeneratecombinations.h"
#include "sc_globalConfig.h"
#include <QCheckBox>
#include <QSettings>
#include <QFileDialog>

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
