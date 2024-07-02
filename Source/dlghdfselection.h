#ifndef DLGHDFSELECTION_H
#define DLGHDFSELECTION_H

#include <QDialog>
#include <QCheckBox>
#include <QHash>


class dlgHdfSelection : public QDialog
{
    Q_OBJECT

public:
    explicit dlgHdfSelection( QString fn, QStringList names, QWidget *parent = nullptr);
    bool isChecked( QString k );
    bool makeSquare();

private slots:
    void onButCloseClicked();

private:
    QHash<QString,QCheckBox*> map;
    QCheckBox *togSquare;
};

#endif // DLGHDFSELECTION_H
