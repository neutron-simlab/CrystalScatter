#ifndef DLGTIMETESTS_H
#define DLGTIMETESTS_H

#include <QDialog>
#include <QVector>


namespace Ui {
class dlgTimeTests;
}

class dlgTimeTests : public QDialog
{
    Q_OBJECT

public:
    explicit dlgTimeTests( QString path, bool gpu, int maxt, QWidget *parent = nullptr);
    ~dlgTimeTests();
    QVector<int> getThreads();
    QVector<int> getHKLmax();
    QVector<int> getQuadrants();
    int getRepetitions();
    QString getSaveFilename();
    QString getComment();
    bool getEnaUpdates();
    bool getEnaSaveImages();

private slots:
    void on_togThread0_toggled(bool) { onTogToggled(); }
    void on_togThread1_toggled(bool) { onTogToggled(); }
    void on_togThread4_toggled(bool) { onTogToggled(); }
    void on_togThreadM_toggled(bool) { onTogToggled(); }
    void on_togHKL2_toggled(bool) { onTogToggled(); }
    void on_togHKL3_toggled(bool) { onTogToggled(); }
    void on_togHKL4_toggled(bool) { onTogToggled(); }
    void on_togHKL5_toggled(bool) { onTogToggled(); }
    void on_grpQuadrants_toggled(bool) { onTogToggled(); }
    void on_togQ1_toggled(bool) { onTogToggled(); }
    void on_togQ2_toggled(bool) { onTogToggled(); }
    void on_togQ4_toggled(bool) { onTogToggled(); }
    void on_butStart_clicked();
    void on_butCancel_clicked();
    void on_togSaveFile_toggled(bool checked);
    void on_butSaveFilename_clicked();

private:
    Ui::dlgTimeTests *ui;
    QString dataPath;

    void onTogToggled();
};

#endif // DLGTIMETESTS_H
