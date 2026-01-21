#ifndef DLGCUSTOMCOLORS_H
#define DLGCUSTOMCOLORS_H

#include <QDialog>
#include <QLabel>
#include <QCheckBox>
#include <QSpinBox>
#include <QList>
#include <QCloseEvent>
#include "sc_globalConfig.h"


namespace Ui {
class dlgCustomColors;
}


class myColorHelper
{
public:
    myColorHelper( QString k, QLabel *l, QCheckBox *t, QSpinBox *r, QSpinBox *g, QSpinBox *b );

    void update();
    QString key;
    void save();

private:
    QLabel *lbl;
    QSpinBox *red, *grn, *blu;
    QCheckBox *tog;
};


class dlgCustomColors : public QDialog
{
    Q_OBJECT

public:
    explicit dlgCustomColors(QWidget *parent = nullptr);
    ~dlgCustomColors();

private slots:
    void on_butExit_clicked();
    void on_butSave_clicked();

    void on_togLattice_toggled(bool) { updateColor(BGKeyLattice); }
    void on_inpLatticeR_valueChanged(int) { updateColor(BGKeyLattice); }
    void on_inpLatticeG_valueChanged(int) { updateColor(BGKeyLattice); }
    void on_inpLatticeB_valueChanged(int) { updateColor(BGKeyLattice); }

    void on_togParticle_toggled(bool) { updateColor(BGKeyParticle); }
    void on_inpParticleR_valueChanged(int) { updateColor(BGKeyParticle); }
    void on_inpParticleG_valueChanged(int) { updateColor(BGKeyParticle); }
    void on_inpParticleB_valueChanged(int) { updateColor(BGKeyParticle); }

    void on_togPeak_toggled(bool) { updateColor(BGKeyPeakShape); }
    void on_inpPeakR_valueChanged(int) { updateColor(BGKeyPeakShape); }
    void on_inpPeakG_valueChanged(int) { updateColor(BGKeyPeakShape); }
    void on_inpPeakB_valueChanged(int) { updateColor(BGKeyPeakShape); }

    void on_togOrientation_toggled(bool) { updateColor(BGKeyOrientation); }
    void on_inpOrientationR_valueChanged(int) { updateColor(BGKeyOrientation); }
    void on_inpOrientationG_valueChanged(int) { updateColor(BGKeyOrientation); }
    void on_inpOrientationB_valueChanged(int) { updateColor(BGKeyOrientation); }

    void on_togExperiment_toggled(bool) { updateColor(BGKeyExperiment); }
    void on_inpExperimentR_valueChanged(int) { updateColor(BGKeyExperiment); }
    void on_inpExperimentG_valueChanged(int) { updateColor(BGKeyExperiment); }
    void on_inpExperimentB_valueChanged(int) { updateColor(BGKeyExperiment); }

    void on_togPixel_toggled(bool) { updateColor(BGKeyPixelMan); }
    void on_inpPixelR_valueChanged(int) { updateColor(BGKeyPixelMan); }
    void on_inpPixelG_valueChanged(int) { updateColor(BGKeyPixelMan); }
    void on_inpPixelB_valueChanged(int) { updateColor(BGKeyPixelMan); }

    void on_inpDisTextR_valueChanged(int) { updateColor(BGKeyDisText); }
    void on_inpDisTextG_valueChanged(int) { updateColor(BGKeyDisText); }
    void on_inpDisTextB_valueChanged(int) { updateColor(BGKeyDisText); }

    void on_inpOnly1DR_valueChanged(int) { updateColor(BGKeyOnly1D); }
    void on_inpOnly1DG_valueChanged(int) { updateColor(BGKeyOnly1D); }
    void on_inpOnly1DB_valueChanged(int) { updateColor(BGKeyOnly1D); }

    void on_inpOnly2DR_valueChanged(int) { updateColor(BGKEYOnly2D); }
    void on_inpOnly2DG_valueChanged(int) { updateColor(BGKEYOnly2D); }
    void on_inpOnly2DB_valueChanged(int) { updateColor(BGKEYOnly2D); }

    void on_inpBothDR_valueChanged(int) { updateColor(BGKey1Dand2D); }
    void on_inpBothDG_valueChanged(int) { updateColor(BGKey1Dand2D); }
    void on_inpBothDB_valueChanged(int) { updateColor(BGKey1Dand2D); }

protected:
    void closeEvent(QCloseEvent *event);

private:
    Ui::dlgCustomColors *ui;

    QList<myColorHelper*> colHelper;
    void updateColor(QString key);
    bool blockUpdate;

};

#endif // DLGCUSTOMCOLORS_H
