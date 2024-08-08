#ifndef MYGUIPARAM_H
#define MYGUIPARAM_H

#include <QWidget>
#include <QCheckBox>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLayout>
#include <QVector>
#include <QFontMetricsF>
#include <QRadioButton>


class myGuiParam;

// Diese Struktur sollte auch alle Daten ohne GUI enthalten können
typedef struct
{
    enum { undef, numdbl, numint, select, toggle, outdbl } type;
    QString key;        // Name des Parameters
    QString tooltip;    // Default ToolTip ohne Fitrange oder andere Zusatzinfos (aus der Liste)
    myGuiParam *gpFit;  // !=0 wenn ein Fit möglich ist (enabled und checked ==> kommt in die Tabelle)
    union
    {
        QDoubleSpinBox *numd;   // Zahl (double)
        QSpinBox  *numi;        // Zahl (int)
        QComboBox *cbs;         // Auswahl (int)
        QCheckBox *tog;         // Bool
        QRadioButton *rad;      // Bool, nur für die spezielle Auswahl von Qmax intern verwendet
        QLineEdit *out;         // Für Outputs (ReadOnly)
        QWidget *w;             // vereinfacht den Zugriff
    } gui;              // nur in der GUI-Version gefüllt
    bool  enabled;
    struct
    {
        double  number; // also for selections
        bool    flag;
    } value;    // immer genutzt
} paramHelper;



class myGuiParam : public QObject
{
    Q_OBJECT
public:
    explicit myGuiParam(QObject *parent = nullptr);

    static void setHideFlag(bool f);

    static void        setAllGuiParams(QWidget *wtab);
    static void        setLabel(QWidget *w);
    static myGuiParam *setSelect(QWidget *w, QStringList sl, QString def);
    static myGuiParam *setInputInt(QWidget *w, int min, int max, int def, QString pref="", QString suff="");
    static myGuiParam *setInputDbl(QWidget *w, double min, double max, int prec, double def, QString pref="", QString suff="");
    static myGuiParam *setInputTog(QWidget *w);
    static myGuiParam *setFitToggle(QWidget *w);
    static void        setFitToggleHelper(QSpinBox *w);
    static void        setFitToggleHelper(QDoubleSpinBox *w);
    static myGuiParam *setOutputDbl(QWidget *w);

    static QComboBox      *getSelect(QString key);
    static QDoubleSpinBox *getInputDbl(QString key);
    static QSpinBox       *getInputInt(QString key);
    static QCheckBox      *getToggle(QString key);

    static myGuiParam *getGuiParam(QString key);
    static myGuiParam *setEnabled(QString key, bool ena, QString lbl="");
    static void        setEnabled(QWidget *w, bool ena);
    static void        setLabel(QString key, QString txt);
    static bool        isEnabled(QString key);
    //static bool        isEnabled(QWidget *w);

    static void updateAdaptiveSteps(bool flg);

    void setLocLabel(QWidget *w, QString on);
    void setLocSelect(QWidget *w, QString on, QStringList sl, QString def);
    void setLocInputInt(QWidget *w, QString on, int min, int max, int def, QString pref="", QString suff="");
    void setLocInputDbl(QWidget *w, QString on, double min, double max, int prec, double def, QString pref="", QString suff="");
    void setLocInputTog(QWidget *w, QString on);
    void setLocFitToggle(QWidget *w, QString on);
    void setLocOutputDbl(QWidget *w, QString on);

    bool isFitEnabled() { return _fit != nullptr; }
    bool isFitUsed();
    void setFitCheck(bool f); // { if ( _fit != nullptr ) { _fit->setChecked(f); _fitUsed=f; } }

    void setValueSel( int v );
    void setValueSel( QString v );
    void setValueInt( int v );
    void setValueDbl( double v );
    void setValue2Dbl( double v1, double v2 );
    void setValue3Dbl( double v1, double v2, double v3 );
    void setValueTog( bool v );

    int     valueSel();
    QString valueSelStr();
    int     valueInt();
    double  valueDbl();
    bool    valueTog();
    void    value2Dbl( double &v1, double &v2 );
    void    value3Dbl( double &v1, double &v2, double &v3 );

    QString debug();
    static void debugGuiParams();

    QString keyName()  { return _keyName;  }
    QString keyName2() { return _keyName2; }
    QString keyName3() { return _keyName3; }

    inline QComboBox *cbs() { return _cbs; }
    inline QSpinBox  *intinp() { return _intinp; }
    inline QDoubleSpinBox *inp() { return _inp; }
    inline QDoubleSpinBox *inp2() { return _inp2; }
    inline QDoubleSpinBox *inp3() { return _inp3; }
    inline QCheckBox *tog() { return _tog; }
    inline QCheckBox *fit() { return _fit; }
    inline QLabel    *lbl() { return _lbl; }
    inline QLineEdit *out() { return _out; }

private slots:
    void fitToggled(bool);

signals:
    void valueChanged();

private:
    QLabel *_lbl;
    QComboBox *_cbs;
    QSpinBox *_intinp;
    QDoubleSpinBox *_inp, *_inp2, *_inp3;
    QCheckBox *_tog;
    QCheckBox *_fit;
    QLineEdit *_out;

    QString _keyName, _keyName2, _keyName3;

    bool _fitUsed;

    paramHelper *_par;

    static myGuiParam *searchParam(QWidget *w, QString &on);
    static QString searchSubKey(QString key);
    static QVector<myGuiParam*> allParams;
    static bool hideValues, oldHideFlag;
    static int  spinBoxMinH;
};

#endif // MYGUIPARAM_H
