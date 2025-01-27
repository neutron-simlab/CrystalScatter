#ifndef DLGGENERATECOMBINATIONS_H
#define DLGGENERATECOMBINATIONS_H

#include <QDialog>

namespace Ui {
class dlgGenerateCombinations;
}

class dlgGenerateCombinations : public QDialog
{
    Q_OBJECT

public:
    explicit dlgGenerateCombinations(QWidget *parent = nullptr);
    ~dlgGenerateCombinations();

    int setValues( QString hdr, QStringList val );
    QStringList getValues( QString hdr );
    QList<int> getValueIDs( QString hdr );

    QString getPath();
    bool getCalcDirect();
    bool getCalcAll();

private slots:
    void on_butDone_clicked();
    void on_butStartGen_clicked();
    void on_butOutPath_clicked();

private:
    Ui::dlgGenerateCombinations *ui;
};

#endif // DLGGENERATECOMBINATIONS_H
