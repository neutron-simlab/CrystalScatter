#ifndef DLGCONFIGAUTOFIT_H
#define DLGCONFIGAUTOFIT_H

#include <QDialog>

namespace Ui {
class dlgConfigAutoFit;
}

class dlgConfigAutoFit : public QDialog
{
    Q_OBJECT

public:
    explicit dlgConfigAutoFit( QWidget *parent = nullptr);
    ~dlgConfigAutoFit();

    QString getFilename();
    bool isLogEnabled();
    bool isLatexEnabled();
    bool isLatexImages();    
    bool isLatexOrgImg();
    bool isLatexCalcImg();
    bool isLatexResiduenImg();
    bool isLatexTrend();
    bool isLatexInput();
    bool isLatexComments();
    bool isDeleteFiles();
    bool isMinimalOutput();
    bool isShowResult();

private slots:
    void on_butCancel_clicked();
    void on_butStart_clicked();
    void on_butInputfile_clicked();
    void on_butEditFile_clicked();

private:
    Ui::dlgConfigAutoFit *ui;
};

#endif // DLGCONFIGAUTOFIT_H
