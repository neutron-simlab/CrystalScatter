#include "dlghdfselection.h"
#include <QVBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QApplication>


dlgHdfSelection::dlgHdfSelection(QString fn, QStringList names, QWidget *parent ) :
    QDialog(parent)
{
    setWindowTitle(fn);
    QVBoxLayout *lay = new QVBoxLayout(this);
    QLabel *lbl1 = new QLabel("The following datasets found in the HDF5 input file:");
    lbl1->setWordWrap(true);
    lay->addWidget(lbl1);
    foreach (QString k, names)
    {
        QCheckBox *b = new QCheckBox( k );
        b->setChecked(true);
        lay->addWidget(b);
        map.insert( k, b );
    }
    QLabel *lbl2 = new QLabel("Uncheck all datasets you do not want to load.");
    lbl2->setWordWrap(true);
    lay->addWidget(lbl2);
    QPushButton *butClose = new QPushButton("Finished and continue");
    connect( butClose, SIGNAL(clicked()), this, SLOT(on_butClose_clicked()) );
    lay->addWidget( butClose );
    setLayout(lay);
    show();
    qApp->restoreOverrideCursor();  // Damit der Cursor wieder zum normalen Pfeil wird
}

void dlgHdfSelection::on_butClose_clicked()
{
    accept();
}

bool dlgHdfSelection::isChecked( QString k )
{
    if ( map.contains(k) )
        return map.value(k)->isChecked();
    return false;
}
