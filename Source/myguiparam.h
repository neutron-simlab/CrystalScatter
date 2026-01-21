#ifndef MYGUIPARAM_H
#define MYGUIPARAM_H

#include <QWidget>
#include <QCheckBox>
#include <QComboBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QGroupBox>
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
    typedef enum { only1d, only2d, both } _use1d2d;
    _use1d2d use1d2d;
    QString key;        // Name des Parameters
    QString tooltip;    // Default ToolTip ohne Fitrange oder andere Zusatzinfos (aus der Liste)
    QString deflabel;   // Default prompt label
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
    QString bgkey;       // optional für die Hintergrundfarbe
} paramHelper;



class myGuiParam : public QObject
{
    Q_OBJECT
public:
    explicit myGuiParam(QObject *parent = nullptr);

    static void setHideFlag(bool f);
    static inline bool hideFlag() { return hideValues; }

    typedef enum { showColNone, showCol1d2d, showColParams } _showColorLabels;
    static void setColBkgFlag(_showColorLabels f, bool force=false);
    static QHash<QString,QColor> colorForBkg;
    static _showColorLabels showColBkg;
    static QColor ColBkg1D, ColBkg2D, ColBkgBoth, ColFgDis, ColFgNorm;

    static void setTooltipFlag(bool l) { longerTooltips=l; }
    static bool getTooltipFlag() { return longerTooltips; }
    static void updateTooltip(myGuiParam *gp, QString spec="");
    static void localUpdateTooltipHelper(QWidget *w, QString ttdef, QString fitrg, QString spec);
    static void updateLabelTooltips();

    static void setParamsDimUsage(paramHelper::_use1d2d f );

    static void        setAllGuiParams(QWidget *wtab);
    static void        setLabel(QWidget *w);
    static myGuiParam *setSelect(QWidget *w, QStringList sl, QString def);
    static myGuiParam *setInputInt(QWidget *w, int min, int max, int def, QString pref="", QString suff="");
    static myGuiParam *setInputDbl(QWidget *w, double min, double max, int prec, double def, QString pref="", QString suff="");
    static myGuiParam *setInputTog(QWidget *w);
    static myGuiParam *setFitToggle(QWidget *w);
    static myGuiParam *setOutputDbl(QWidget *w);
    static void        setInputFieldSize( QDoubleSpinBox *box, double max, int prec, QString pref, QString suff );

    static QComboBox      *getSelect(QString key);
    static QDoubleSpinBox *getInputDbl(QString key);
    static QSpinBox       *getInputInt(QString key);
    static QCheckBox      *getToggle(QString key);

    static myGuiParam *getGuiParam(QString key);
    static myGuiParam *setEnabled(QString key, bool ena, QString lbl="", QString tt="");
    static myGuiParam *setEnabled(myGuiParam *gp, bool ena, QString lbl="", QString tt="");
    static void        setEnabled(QWidget *w, bool ena);
    static void        setLabel(QString key, QString txt);
    static bool        isEnabled(QString key);
    static bool        isEnabled(myGuiParam *gp);
    //static bool        isEnabled(QWidget *w);

    bool isParEnabled() { if (_par==nullptr) return false; return _par->enabled; }

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
    void setFitCheck(bool f, QString dbg); // { if ( _fit != nullptr ) { _fit->setChecked(f); _fitUsed=f; } }

    QString debug(bool nurvis=false);
    static void debugGuiParams(bool nurvis=false, QString fn="");

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

    static myGuiParam *searchParam(QWidget *w, QString &on, bool search=true);
    static QString searchSubKey(QString key);
    static QVector<myGuiParam*> allParams;
    static bool hideValues;
    static bool longerTooltips;
    static paramHelper::_use1d2d use12Flag;
    static int  spinBoxMinH;
    static bool localSetLongerTooltip(QWidget *w); // in myguiparam_ToolTipsLong.cpp

    static QHash<QString,QList<myGuiParam*>> grpLabels;
    static QHash<QString,QGroupBox*> grpBoxes;

    static void generateTestFile();

    static void setLabelColor( QLabel *lbl, QColor fg, QColor bg=QColor() );
    static void setLabelColor( myGuiParam *gp, bool ena );
};

#endif // MYGUIPARAM_H
