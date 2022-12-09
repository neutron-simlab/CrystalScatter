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

private slots:
    void on_butClose_clicked();

private:
    QHash<QString,QCheckBox*> map;
};

#endif // DLGHDFSELECTION_H
