#include "dlgconfigautofit.h"
#include "ui_dlgconfigautofit.h"
#include "sc_globalConfig.h"
#include <QFileDialog>
#include <QSettings>
#include <QProcess>


dlgConfigAutoFit::dlgConfigAutoFit(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::dlgConfigAutoFit)
{
    ui->setupUi(this);

    QSettings data(SETT_APP,SETT_PAR);
    ui->inpInputfile->setText( data.value("LastAutoFit",".").toString() );
    ui->togRunLogs->setChecked( data.value("AutoFitRunLogs",true).toBool() );
    ui->grpLatexEnabled->setChecked( data.value("AutoFitLatexEna",true).toBool() );
    ui->grpImages->setChecked( data.value("AutoFitLatexImages",true).toBool() );
    ui->togOrgImgStep->setChecked( data.value("AutoFitLatexOrgImgStep",true).toBool() );
    ui->togCalcStep->setChecked( data.value("AutoFitLatexCalcImgStep",true).toBool() );
    ui->togResiStep->setChecked( data.value("AutoFitLatexResiduenStep",false).toBool() );
    ui->togLatexTrend->setChecked( data.value("AutoFitLatexTrend",true).toBool() );
    ui->togLatexInput->setChecked( data.value("AutoFitLatexInput",true).toBool() );
    ui->togLatexComments->setChecked( data.value("AutoFitLatexComment",false).toBool() );
    ui->togLatexComments->setEnabled( ui->togLatexInput->isChecked() );
    ui->togDelFiles->setChecked( data.value("AutoFitDelFiles",false).toBool() );
    ui->togMinimalOut->setChecked( data.value("AutoFitMinimalOutput",false).toBool() );
    ui->togShowResult->setChecked( data.value("AutoFitShowResult",false).toBool() );

    adjustSize();
}

dlgConfigAutoFit::~dlgConfigAutoFit()
{
    delete ui;
}

void dlgConfigAutoFit::on_butCancel_clicked()
{
    reject();
}

void dlgConfigAutoFit::on_butStart_clicked()
{
    QSettings data(SETT_APP,SETT_PAR);
    data.setValue( "LastAutoFit",              ui->inpInputfile->text() );
    data.setValue( "AutoFitRunLogs",           ui->togRunLogs->isChecked() );
    data.setValue( "AutoFitLatexEna",          ui->grpLatexEnabled->isChecked() );
    data.setValue( "AutoFitLatexImages",       ui->grpImages->isChecked() );
    data.setValue( "AutoFitLatexOrgImgStep",   ui->togOrgImgStep->isChecked() );
    data.setValue( "AutoFitLatexCalcImgStep",  ui->togCalcStep->isChecked() );
    data.setValue( "AutoFitLatexResiduenStep", ui->togResiStep->isChecked() );
    data.setValue( "AutoFitLatexTrend",        ui->togLatexTrend->isChecked() );
    data.setValue( "AutoFitLatexInput",        ui->togLatexInput->isChecked() );
    data.setValue( "AutoFitLatexComment",      ui->togLatexComments->isChecked() );
    data.setValue( "AutoFitDelFiles",          ui->togDelFiles->isChecked() );
    data.setValue( "AutoFitMinimalOutput",     ui->togMinimalOut->isChecked() );
    data.setValue( "AutoFitShowResult",        ui->togShowResult->isChecked() );
    accept();
}

void dlgConfigAutoFit::on_butInputfile_clicked()
{
    QString fn = QFileDialog::getOpenFileName(this,"Automatic sequence",ui->inpInputfile->text());
    if ( fn.isEmpty() ) return;
    ui->inpInputfile->setText(fn);
}


QString dlgConfigAutoFit::getFilename()  { return ui->inpInputfile->text(); }
bool dlgConfigAutoFit::isLogEnabled()    { return ui->togRunLogs->isChecked(); }
bool dlgConfigAutoFit::isLatexEnabled()  { return ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isLatexImages()   { return ui->grpImages->isChecked()        && ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isLatexOrgImg()   { return ui->togOrgImgStep->isChecked()    && ui->grpImages->isChecked() && ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isLatexCalcImg()  { return ui->togCalcStep->isChecked()      && ui->grpImages->isChecked() && ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isLatexResiduenImg() { return ui->togResiStep->isChecked()   && ui->grpImages->isChecked() && ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isLatexTrend()    { return ui->togLatexTrend->isChecked()    && ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isLatexInput()    { return ui->togLatexInput->isChecked()    && ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isLatexComments() { return ui->togLatexComments->isChecked() && ui->togLatexInput->isChecked() && ui->grpLatexEnabled->isChecked(); }
bool dlgConfigAutoFit::isDeleteFiles()   { return ui->togDelFiles->isChecked(); }
bool dlgConfigAutoFit::isMinimalOutput() { return ui->togMinimalOut->isChecked(); }
bool dlgConfigAutoFit::isShowResult()    { return ui->togShowResult->isChecked(); }


void dlgConfigAutoFit::on_butEditFile_clicked()
{
    QSettings data(SETT_APP,SETT_PAR);
    QString editor = data.value("AutoFitEditor","").toString();
    if ( editor.isEmpty() )
    {
        editor = QFileDialog::getOpenFileName( this, "Select Editor", "." );
        if ( editor.isEmpty() ) return;
        data.setValue("AutoFitEditor",editor);
    }
    QProcess::startDetached( editor, QStringList()<<ui->inpInputfile->text() );
}
