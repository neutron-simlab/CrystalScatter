#include "dlghdfselection.h"
#include <QVBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QApplication>
#include <QPalette>
#include <QFileInfo>
#include <QFont>


dlgHdfSelection::dlgHdfSelection(QString fn, QStringList names, QWidget *parent ) :
    QDialog(parent)
{
    setWindowTitle(QFileInfo(fn).fileName());
    QVBoxLayout *lay = new QVBoxLayout(this);
    QLabel *lbl1 = new QLabel("The following datasets found in the HDF5 input file:");
    lbl1->setWordWrap(true);
    lay->addWidget(lbl1);
    bool defcheck = true;
    foreach (QString k, names)
    {
        QCheckBox *b = new QCheckBox( k );
        b->setChecked(defcheck);
        lay->addWidget(b);
        map.insert( k, b );
        defcheck = false; // nur das erste Image per Default einlesen
    }
    QLabel *lbl2 = new QLabel("Uncheck all datasets you do not want to load.");
    lbl2->setWordWrap(true);
    QFont fnt = lbl2->font();
    fnt.setItalic(true);
    lbl2->setFont(fnt);
    lay->addWidget(lbl2);

    togSquare = new QCheckBox("Make all images square if dimensions differ");
    togSquare->setChecked(true);
    QPalette pal = togSquare->palette();
    pal.setBrush( QPalette::WindowText, Qt::darkBlue );
    togSquare->setPalette(pal);
    lay->addWidget(togSquare);

    QPushButton *butClose = new QPushButton("Finished and continue");
    connect( butClose, SIGNAL(clicked()), this, SLOT(onButCloseClicked()) );
    lay->addWidget( butClose );
    setLayout(lay);
    show();
    qApp->restoreOverrideCursor();  // Damit der Cursor wieder zum normalen Pfeil wird
}

void dlgHdfSelection::onButCloseClicked()
{
    accept();
}

bool dlgHdfSelection::isChecked( QString k )
{
    if ( map.contains(k) )
        return map.value(k)->isChecked();
    return false;
}

bool dlgHdfSelection::makeSquare()
{
    return togSquare->isChecked();
}
