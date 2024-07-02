#include "myguiparam.h"
#include <QDebug>
#include "sc_maingui.h"

#define TT(x) //x // Untersuchungen zum ToolTip ...
#define FF(x) //x // Untersuchungen bzgl. Fit-Toggle



/*static*/ QVector<myGuiParam*> myGuiParam::allParams;
/*static*/ bool myGuiParam::hideValues  = true;
/*static*/ bool myGuiParam::oldHideFlag = true;
/*static*/ int  myGuiParam::spinBoxMinH = -1;



myGuiParam::myGuiParam(QObject *parent) : QObject(parent)
{
    _lbl      = nullptr;
    _cbs      = nullptr;
    _intinp   = nullptr;
    _inp      = nullptr;
    _inp2     = nullptr;
    _inp3     = nullptr;
    _tog      = nullptr;
    _fit      = nullptr;
    _out      = nullptr;
    _keyName  = "??";
    _keyName2 = "??";
    _keyName3 = "??";
    _par      = nullptr;
    allParams.append(this);
}


/*static*/ myGuiParam *myGuiParam::searchParam(QWidget *w, QString &on)
{
    QStringList sl = w->objectName().mid(3).split("_");
    if ( sl.size() == 1 )
        sl << sl[0];
    else if ( sl[1].at(0).isDigit() )
        sl[1] = sl[0];
    //TT( if ( w->objectName().contains("Pixel") ) qDebug() << "  search" << w->objectName() << sl );
    on = sl[0];
    foreach ( myGuiParam*p, allParams )
    {
        if ( sl[1].compare(p->keyName()) == 0 )
        {
            return p;
        }
    }
    return nullptr;
}

/*static*/void myGuiParam::setLabel( QWidget *w )
{
    QString on;
    myGuiParam *p = searchParam(w,on);
    if ( p == nullptr ) p = new myGuiParam;
    p->setLocLabel(w,on);
}
void myGuiParam::setLocLabel( QWidget *w, QString on )
{
    _lbl = static_cast<QLabel*>(w);
    //_lbl->setAlignment(Qt::AlignRight|Qt::AlignVCenter);
    _keyName = on;
    //w->setEnabled(false);
}


/*static*/myGuiParam *myGuiParam::setSelect( QWidget *w, QStringList sl, QString def )
{
    QString on;
    myGuiParam *p = searchParam(w,on);
    if ( p == nullptr ) p = new myGuiParam;
    p->setLocSelect(w,on,sl,def);
    return p;
}
void myGuiParam::setLocSelect( QWidget *w, QString on, QStringList sl, QString def )
{
    _cbs = static_cast<QComboBox*>(w);
    _cbs->clear();
    _cbs->addItems(sl);
#ifdef CALC_INPUT_MAXWIDTH
    _cbs->setMaximumWidth(CALC_INPUT_MAXWIDTH);
#endif
    // Calculate width of ComboBox-List to make all item readable.
    QFontMetrics fm = _cbs->fontMetrics();
    int ws = 0;
    foreach ( QString s, sl )
    {
        if ( fm.horizontalAdvance(s) > ws )
            ws = fm.horizontalAdvance(s);
    }
    _cbs->view()->setMinimumWidth(ws*1.1);
    _cbs->setMaxVisibleItems(sl.size()); // -> no scroll bar
    _keyName = on;
    if ( !def.isEmpty() )
    {   // Defaultwert setzen
        _cbs->setCurrentIndex( _cbs->findText(def,Qt::MatchStartsWith) );
    }
    _cbs->connect( _cbs, SIGNAL(currentIndexChanged(int)),
                   SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
    DT( qDebug() << "myGuiParam::Cbs  connect" << on );
    //w->setEnabled(false);
}


/*sttaic*/myGuiParam *myGuiParam::setInputInt( QWidget *w, int min, int max, int def, QString pref, QString suff )
{
    QString on;
    myGuiParam *p = searchParam(w,on);
    if ( p == nullptr ) p = new myGuiParam;
    p->setLocInputInt(w,on,min,max,def,pref,suff);
    return p;
}
void myGuiParam::setLocInputInt(QWidget *w, QString on, int min, int max, int def, QString pref, QString suff )
{
    _intinp = static_cast<QSpinBox*>(w);
    _intinp->setRange( min, max );
    if ( !suff.isEmpty() ) _intinp->setSuffix(suff);
    if ( !pref.isEmpty() ) _intinp->setPrefix(pref);
    _intinp->setValue(def);
#ifdef CALC_INPUT_MAXWIDTH
    //_intinp->setMaximumWidth(CALC_INPUT_MAXWIDTH);
#endif
    if ( spinBoxMinH < 0 )
    {
        spinBoxMinH = _intinp->minimumSizeHint().height();
        //qDebug() << "spinBoxMinH" << spinBoxMinH;
    }
    else
        _intinp->setMinimumHeight(spinBoxMinH);
    _keyName = on;
    //w->setEnabled(false);
    _intinp->connect( _intinp, SIGNAL(valueChanged(int)),
                     SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
    DT( qDebug() << "myGuiParam::Int  connect" << on );
    //QString tt = QString("Fitrange: %1 .. %2").arg(min).arg(max);
    //if ( w->toolTip().isEmpty() )
    //    w->setToolTip(tt);
    //else
    //    w->setToolTip( w->toolTip()+"\n"+tt );
}


/*static*/ myGuiParam *myGuiParam::setOutputDbl(QWidget *w)
{
    QString on;
    myGuiParam *p = searchParam(w,on);
    if ( p == nullptr ) p = new myGuiParam;
    p->setLocOutputDbl(w,on);
    return p;
}
void myGuiParam::setLocOutputDbl(QWidget *w, QString on)
{
    if ( _out == nullptr )
    {
        _out = static_cast<QLineEdit*>(w);
#ifdef CALC_INPUT_MAXWIDTH
        //_inp->setMaximumWidth(CALC_INPUT_MAXWIDTH);
#endif
        _keyName = on;
        //w->setEnabled(false);
        DT( qDebug() << "myGuiParam::OutDbl  connect" << on );
    }
}


/*static*/myGuiParam *myGuiParam::setInputDbl(QWidget *w, double min, double max, int prec, double def, QString pref, QString suff )
{
    QString on;
    myGuiParam *p = searchParam(w,on);
    if ( p == nullptr ) p = new myGuiParam;
    p->setLocInputDbl(w,on,min,max,prec,def,pref,suff);
    return p;
}
void myGuiParam::setLocInputDbl(QWidget *w, QString on, double min, double max, int prec, double def, QString pref, QString suff )
{
    if ( _inp == nullptr )
    {
        _inp = static_cast<QDoubleSpinBox*>(w);
        _inp->setRange( min, max );
        _inp->setDecimals(prec);
        if ( prec > 1 ) _inp->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
        if ( !suff.isEmpty() ) _inp->setSuffix(suff);
        if ( !pref.isEmpty() ) _inp->setPrefix(pref);
        _inp->setValue(def);
#ifdef CALC_INPUT_MAXWIDTH
        //_inp->setMaximumWidth(CALC_INPUT_MAXWIDTH);
#endif
        if ( spinBoxMinH < 0 )
        {
            spinBoxMinH = _inp->minimumSizeHint().height();
            //qDebug() << "spinBoxMinH" << spinBoxMinH;
        }
        else
            _inp->setMinimumHeight(spinBoxMinH);
        _keyName = on;
        //w->setEnabled(false);
        if ( suff=="°" && min==0 && max==360 )
             _inp->setWrapping(true); // Bei Winkeln aktiviert
        _inp->connect( _inp, SIGNAL(valueChanged(double)),
                       SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
        DT( qDebug() << "myGuiParam::Dbl  connect" << on );
    }
    else if ( _inp2 == nullptr )
    {
        _inp2 = static_cast<QDoubleSpinBox*>(w);
        _inp2->setRange( min, max );
        _inp2->setDecimals(prec);
        if ( prec > 1 ) _inp2->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
        if ( !suff.isEmpty() ) _inp2->setSuffix(suff);
        if ( !pref.isEmpty() ) _inp2->setPrefix(pref);
        _inp2->setValue(def);
#ifdef CALC_INPUT_MAXWIDTH
        //_inp2->setMaximumWidth(CALC_INPUT_MAXWIDTH);
#endif
        if ( spinBoxMinH < 0 )
        {
             spinBoxMinH = _inp2->minimumSizeHint().height();
             //qDebug() << "spinBoxMinH" << spinBoxMinH;
        }
        else
             _inp2->setMinimumHeight(spinBoxMinH);
        _keyName2 = on;
        //w->setEnabled(false);
        _inp2->connect( _inp2, SIGNAL(valueChanged(double)),
                        SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
        DT( qDebug() << "myGuiParam::Dbl2 connect" << on );
    }
    else if ( _inp3 == nullptr )
    {
        _inp3 = static_cast<QDoubleSpinBox*>(w);
        _inp3->setRange( min, max );
        _inp3->setDecimals(prec);
        if ( prec > 1 ) _inp3->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
        if ( !suff.isEmpty() ) _inp3->setSuffix(suff);
        if ( !pref.isEmpty() ) _inp3->setPrefix(pref);
        _inp3->setValue(def);
#ifdef CALC_INPUT_MAXWIDTH
        //_inp3->setMaximumWidth(CALC_INPUT_MAXWIDTH);
#endif
        if ( spinBoxMinH < 0 )
        {
             spinBoxMinH = _inp3->minimumSizeHint().height();
             //qDebug() << "spinBoxMinH" << spinBoxMinH;
        }
        else
             _inp3->setMinimumHeight(spinBoxMinH);
        _keyName3 = on;
        //w->setEnabled(false);
        _inp3->connect( _inp3, SIGNAL(valueChanged(double)),
                        SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
        DT( qDebug() << "myGuiParam::Dbl3 connect" << on );
    }
    else
        return;
    TT( if ( w->objectName().contains("Pixel") ) qDebug() << "  set" << w->objectName() << debug() );
    //QString tt = "TEST "+on+"\n" + w->objectName(); // QString("Fitrange: %1 .. %2").arg(min).arg(max);
    //if ( w->toolTip().isEmpty() )
    //    w->setToolTip(tt);
    //else
    //    w->setToolTip( w->toolTip()+"\n"+tt );
}


/*static*/myGuiParam *myGuiParam::setInputTog( QWidget *w )
{
    QString on;
    myGuiParam *p = searchParam(w,on);
    if ( p == nullptr )
    {
        if ( on.startsWith("Expand") || on.startsWith("Auto") ) return nullptr;
        p = new myGuiParam;
    }
    p->setLocInputTog(w,on);
    return p;
}
void myGuiParam::setLocInputTog( QWidget *w, QString on )
{
    _tog = static_cast<QCheckBox*>(w);
    _keyName = on;
    //w->setEnabled(false);
    // Es werden für QMax spezielle Toggles (Radiobuttons) verwendet, um
    // das richtige Input-Feld zu wählen... Diese sollten aber kein
    // Recalc auslösen, das käme nämlich zweimal...
    // Das gleiche gilt für die BeamCenter Auswahl.
    // Daher wird hier einer der beiden Radiobuttons ausgeblendet.
    if ( ! w->objectName().startsWith("radEditQmaxPreset") &&
         ! w->objectName().startsWith("radCenterMidpoint") )
    {
        _tog->connect( _tog, SIGNAL(toggled(bool)),
                       SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
        // für die Markierungen mit class setColHelper muss der Hintergrund gezeichnet werden
        _tog->setAutoFillBackground(true);
        QPalette pal = _tog->palette();
        pal.setBrush(QPalette::Button,Qt::white); // ... und der Default gesetzt werden.
        _tog->setPalette(pal);
        DT( qDebug() << "myGuiParam::Tog  connect" << on );
    }
    /*else
    {
        qDebug() << "setLocInputTog::radEditQmaxPreset / radCenterMidpoint";
    }*/
}


/*static*/myGuiParam *myGuiParam::setFitToggle( QWidget *w )
{
    QString on;
    myGuiParam *p = searchParam(w,on);
    if ( p == nullptr ) p = new myGuiParam;
    p->setLocFitToggle(w,on);
    return p;
}
/*static*/ void myGuiParam::setFitToggleHelper(QSpinBox *w)
{
    if ( w == nullptr ) return;
    if ( w->toolTip().contains("Fitrange") ) return;
    QString tt = QString("Fitrange: %1 .. %2").arg(w->minimum()).arg(w->maximum());
    if ( w->toolTip().isEmpty() )
        w->setToolTip(tt);
    else
        w->setToolTip( w->toolTip()+"\n"+tt );
    TT( if ( w->objectName().contains("Pixel") ) qDebug() << "   tt" << w->objectName() << tt );
}
/*static*/ void myGuiParam::setFitToggleHelper(QDoubleSpinBox *w)
{
    if ( w == nullptr ) return;
    if ( w->toolTip().contains("Fitrange") ) return;
    QString tt = QString("Fitrange: %1 .. %2").arg(w->minimum()).arg(w->maximum());
    if ( w->toolTip().isEmpty() )
        w->setToolTip(tt);
    else
        w->setToolTip( w->toolTip()+"\n"+tt );
    TT( if ( w->objectName().contains("Pixel") ) qDebug() << "   tt" << w->objectName() << tt );
}
void myGuiParam::setLocFitToggle( QWidget *w, QString on )
{
    _fit = static_cast<QCheckBox*>(w);
    _fit->setToolTip("If checked, use this variable in the Simplex 2D Fit");
    _keyName = on;
    _fit->connect( _fit, SIGNAL(toggled(bool)), this, SLOT(fitToggled(bool)) );
    _fitUsed = false;
}
void myGuiParam::fitToggled(bool f)
{
    FF( if ( _fit->objectName().startsWith("fitSig") ) qDebug() << "fitToggled" << f << _fit->objectName() );
    _fitUsed = f;
}


bool myGuiParam::isFitUsed()
{
    if ( _fit == nullptr ) return false;
    FF( if ( _fit->objectName().startsWith("fitSig") ) qDebug() << "isFitUsed" << _fit->isEnabled() << _fit->isChecked() /*<< _fitUsed*/ << _fit->objectName() );
    return _fit->isEnabled() && _fit->isChecked(); // && _fitUsed;
}

void myGuiParam::setFitCheck(bool f)
{
    if ( _fit == nullptr ) return;
    FF( if ( _fit->objectName().startsWith("fitSig") ) qDebug() << "setFitCheck" << f << _fit->objectName() );
    _fit->setChecked(f);
    _fitUsed=f;
}


void myGuiParam::setValueSel( int v )
{
    Q_UNUSED(v)
}

void myGuiParam::setValueSel( QString v )
{
    Q_UNUSED(v)
}

void myGuiParam::setValueInt( int v )
{
    Q_UNUSED(v)
}

void myGuiParam::setValueDbl( double v )
{
    Q_UNUSED(v)
}

void myGuiParam::setValue2Dbl( double v1, double v2 )
{
    Q_UNUSED(v1)
    Q_UNUSED(v2)
}

void myGuiParam::setValue3Dbl( double v1, double v2, double v3 )
{
    Q_UNUSED(v1)
    Q_UNUSED(v2)
    Q_UNUSED(v3)
}

void myGuiParam::setValueTog( bool v )
{
    Q_UNUSED(v)
}


int myGuiParam::valueSel()
{
    return 0;
}

QString myGuiParam::valueSelStr()
{
    return "";
}

int myGuiParam::valueInt()
{
    return 0;
}

double myGuiParam::valueDbl()
{
    return 0;
}

bool myGuiParam::valueTog()
{
    return false;
}

void myGuiParam::value2Dbl( double &v1, double &v2 )
{
    Q_UNUSED(v1)
    Q_UNUSED(v2)
}

void myGuiParam::value3Dbl( double &v1, double &v2, double &v3 )
{
    Q_UNUSED(v1)
    Q_UNUSED(v2)
    Q_UNUSED(v3)
}


QString myGuiParam::debug()
{
    QString str = "";
    if ( _lbl != nullptr )
        str += _lbl->text();
    if ( _cbs != nullptr )
        str += QString(" Select{%1:%2} (%3)").arg(_cbs->count()).arg(_cbs->currentText(),_keyName);
    if ( _intinp != nullptr )
        str += " IntInp ("+_keyName+")";
    if ( _inp != nullptr )
        str += " Input ("+_keyName+")";
    if ( _inp2 != nullptr )
        str += " Input2 ("+_keyName2+")";
    if ( _inp3 != nullptr )
        str += " Input3 ("+_keyName3+")";
    if ( _tog != nullptr )
        str += " Toggle ("+_keyName+")";
    if ( _fit != nullptr )
        str += " FitEnabled";
    return str;
}


/*static*/ void myGuiParam::debugGuiParams()
{
    if ( allParams.size() == 0 ) return;
    qDebug() << "myGuiParam::debugGuiParams() - Start";
    qDebug() << "Anzahl:" << allParams.size();
    foreach ( myGuiParam *p, allParams )
    {
        qDebug() << p->debug();
    }
    qDebug() << "myGuiParam::debugGuiParams() - Finished";
}


/*static*/ void myGuiParam::setAllGuiParams(QWidget *wtab)
{
    DT( qDebug() << "setAllGuiParams()" );
    if ( SC_MainGUI::getInst()->getCalcGui() == nullptr )
    {
        qDebug() << "ERR1";
        return;
    }
    if ( SC_MainGUI::getInst()->getCalcGui()->getCalcPtr() == nullptr )
    {
        qDebug() << "ERR2";
        return;
    }
    // Im ersten Schritt die Namensliste der Parameter aufbauen
    QStringList names;
    QStringList meta = SC_MainGUI::getInst()->getCalcGui()->getCalcPtr()->guiLayoutNeu();
    foreach ( QString mm, meta )
    {
        QStringList sl = mm.split(";");
        names << sl[1].trimmed();
    }

    //  /*name=QString()*/, Qt::FindChildOptions options = Qt::FindChildrenRecursively) const
    QList<QWidget*> wid = wtab->findChildren<QWidget*>();
    for ( int us=0; us<2; us++ )
    {
        foreach ( QWidget *w, wid )
        {
            if ( us == 0 &&   w->objectName().contains("_") ) continue;
            if ( us == 1 && ! w->objectName().contains("_") ) continue;
            // Unten werden für die Eingabefelder die Tooltips gesetzt. Aus Tests habe ich gelernt,
            // dass die Widgetnamen mit '_' (d.h. inp2, inp3) immer im Nachgang bearbeitet werden
            // müssen, damit immer die richtigen Gruppen angelegt werden.

            QString on;
            myGuiParam *gp = searchParam(w, on);
            if ( ! names.contains(on) )
            {
                //qDebug() << "no name" << on;
                continue;
            }

            QStringList sl;
            foreach ( QString mm, meta )
            {
                if ( mm.indexOf(";"+on+";") > 0 )
                {
                    sl = mm.split(";");
                    break;
                }
            }
            while ( sl.size() < 5 ) sl << " ";  // tooltip;default sind optional

            // Each element must be in the form "type;kenn;typespec;tooltip;default" with:
            //  [0]type     = C:Selection, N:Double, I:Integer, T:Toggle, O:DoubleOut
            //  [1]kenn     = internal parameter name to connect to correct gui element
            //  [2]typespec = C:Selection : "...|...|..."  (required)
            //                N:Double    : "frac|min|max|unit"  (optional, default: "2|-10000|+10000|")
            //                I:Integer   : "min|max|unit"  (optional, default: "-10000|+10000|")
            //                T:Toggle    : (empty)
            //                O:DoubleOut : (empty)
            //  [3]tooltip  = this is the tooltip set to both prompt label and inputfield (optional)
            //  [4]default  = the default value (optional)

            if ( ! sl[3].trimmed().isEmpty() )
                w->setToolTip(sl[3]);

            QStringList typval = sl[2].split("|",SPLIT_SKIP_EMPTY_PARTS);

            paramHelper *par = SC_MainGUI::getInst()->getCalcGui()->params.value(on,nullptr);
            if ( par == nullptr )
            {
                qDebug() << "getpar" << on << "FAIL";
                continue;
            }
            //else
            //    qDebug() << "getpar" << on << par->type;
            par->gpFit = nullptr;
            par->tooltip = sl[3].trimmed();

            if ( w->objectName().startsWith("lbl") )
            {
                myGuiParam::setLabel(w);
            }
            else if ( w->objectName().startsWith("cbs") )
            {   //  [2]typespec = C:Selection : "...|...|..."  (required)
                // Weil die ComboBoxen in der neuen Version mit Signalen versehen werden um andere Elemente freizugeben,
                // sollten die Signale beim Füllen geblockt werden.
                w->blockSignals(true);
                for ( int i=0; i<typval.size(); i++ ) typval[i] = typval[i]+QString(" {%1}").arg(i);
                myGuiParam::setSelect(w, typval, sl[4]);
                //if ( par->gui.w != nullptr ) { par->gui.w->deleteLater(); }
                par->gui.w = w;
                w->blockSignals(false);
            }
            else if ( w->objectName().startsWith("int") )
            {   // [2]typespec = I:Integer   : "min|max|unit"  (optional, default: "-10000|+10000|")
                int min = ( typval.size() > 0 ) ? typval[0].toInt() : 0;
                int max = ( typval.size() > 1 ) ? typval[1].toInt() : 1000;
                QString suff = ( typval.size() > 2 ) ? typval[2] : "";
                int def = ( sl[4].isEmpty() ) ? 0 : sl[4].toInt();
                myGuiParam::setInputInt(w, min, max, def, "", suff);
                //if ( par->gui.w != nullptr ) { par->gui.w->deleteLater(); }
                par->gui.w = w;
            }
            else if ( w->objectName().startsWith("inp") || w->objectName() == "outCalcQmax" )
            {   // [2]typespec = N:Double    : "frac|min|max|unit"  (optional, default: "2|-10000|+10000|")
                int frac = ( typval.size() > 0 ) ? typval[0].toInt() : 2;
                double min = ( typval.size() > 1 ) ? typval[1].toDouble() : -10000.0;
                double max = ( typval.size() > 2 ) ? typval[2].toDouble() : +10000.0;
                QString suff = ( typval.size() > 3 ) ? typval[3] : "";
                double def = ( sl[4].isEmpty() ) ? 0 : sl[4].toDouble();
                myGuiParam::setInputDbl(w, min, max, frac, def, "", suff);
                //if ( par->gui.w != nullptr ) { par->gui.w->deleteLater(); }
                par->gui.w = w;
            }
            else if ( w->objectName().startsWith("fit") )
            {
                myGuiParam::setFitToggle(w);
            }
            else if ( w->objectName().startsWith("tog") ||  w->objectName().startsWith("rad") )
            {   // [2]typespec = T:Toggle    : (empty)
                myGuiParam::setInputTog(w);
                //if ( par->gui.w != nullptr ) { par->gui.w->deleteLater(); }
                par->gui.w = w;
            }
            //else if ( w->objectName().startsWith("grp") )
            //    qDebug() << "Group " << w->objectName();
            //else if ( w->objectName().startsWith("rad") )
            //    qDebug() << "Radio " << w->objectName();
            else if ( w->objectName().startsWith("out") )
            {   // [2]typespec = O:DoubleOut : (empty)
                myGuiParam::setOutputDbl(w);
                par->gui.w = w;
            }
            else
                qDebug() << "???" << w->objectName() << on;
        }
    }
    // Da es sein kann, dass der Fit-Toggle erst nach dem passenden Input gefunden und eingetragen wird,
    // muss hier ein zweiter Durchlauf gemacht weren, wobei nur bei den Inputs der Fit-Toggle gesucht wird.
    foreach ( QWidget *w, wid )
    {
        QString on;
        myGuiParam *hlpPar = searchParam(w, on);
        if ( ! names.contains(on) ) continue;
        paramHelper *par = SC_MainGUI::getInst()->getCalcGui()->params[on];
        if ( par == nullptr ) continue;
        par->gpFit = hlpPar; // ->fit();
        hlpPar->_par = par;
        //qDebug() << "setAllGuiParams" << on << "-> paramHelper";
        if ( hlpPar->fit() != nullptr )
        {
            TT( if ( w->objectName().contains("Pixel") ) qDebug() << "ToolTip:" << w->objectName()
                                                               << hlpPar->debug() );
            setFitToggleHelper(hlpPar->intinp());
            setFitToggleHelper(hlpPar->inp());
            setFitToggleHelper(hlpPar->inp2());
            setFitToggleHelper(hlpPar->inp3());
        }
    }
}


/*static*/ QComboBox *myGuiParam::getSelect(QString key)
{
    foreach ( myGuiParam *p, allParams )
    {
        if ( p->keyName() == key )
            return p->cbs();
    }
    return nullptr;
}

/*static*/ QDoubleSpinBox *myGuiParam::getInputDbl(QString key)
{
    foreach ( myGuiParam *p, allParams )
    {
        if ( p->keyName() == key )
            return p->inp();
        if ( p->keyName2() == key )
            return p->inp2();
        if ( p->keyName3() == key )
            return p->inp3();
    }
    return nullptr;
}

/*static*/ QSpinBox *myGuiParam::getInputInt(QString key)
{
    foreach ( myGuiParam *p, allParams )
    {
        if ( p->keyName() == key )
            return p->intinp();
    }
    return nullptr;
}

/*static*/ QCheckBox *myGuiParam::getToggle(QString key)
{
    myGuiParam *gp = getGuiParam(key);
    if ( gp != nullptr ) return gp->tog();
    return nullptr;
}

/*static*/ myGuiParam *myGuiParam::getGuiParam(QString key)
{
    foreach ( myGuiParam *p, allParams )
    {
        if ( p->keyName() == key )
            return p;
        if ( p->keyName2() == key )
            return p;
        if ( p->keyName3() == key )
            return p;
    }
    return nullptr;
}

/*static*/ void myGuiParam::setHideFlag(bool f)
{
    if ( hideValues == f ) return;
    hideValues = f;
    foreach ( myGuiParam *p, allParams )
    {
        // Hier werden alle Werte aus setEnabled() zurückgesetzt, damit sie dann in setEnabled() richtig genutzt werden
        //if ( hideValues )
        {   // Jetzt sollen alle Werte versteckt werden (visible)
            if ( p->cbs()    != nullptr ) { p->cbs()->setVisible(true);    p->cbs()->setEnabled(true); }
            if ( p->intinp() != nullptr ) { p->intinp()->setVisible(true); p->intinp()->setEnabled(true); }
            if ( p->inp()    != nullptr ) { p->inp()->setVisible(true);    p->inp()->setEnabled(true); }
            if ( p->inp2()   != nullptr ) { p->inp2()->setVisible(true);   p->inp2()->setEnabled(true); }
            if ( p->inp3()   != nullptr ) { p->inp3()->setVisible(true);   p->inp3()->setEnabled(true); }
            if ( p->tog()    != nullptr ) { p->tog()->setVisible(true);    p->tog()->setEnabled(true); }
            if ( p->fit()    != nullptr ) { p->fit()->setVisible(true);    p->fit()->setEnabled(true);   /*p->setFitCheck(true);*/ }
            if ( p->lbl()    != nullptr ) { p->lbl()->setVisible(true);    p->lbl()->setEnabled(true); }
        }
        /*else
        {   // Jetzt werden die Werte nur gesperrt, aber die Eingaben bleiben verwendbar
            //if ( p->cbs()    != nullptr ) p->cbs()->setEnabled(true);
            //if ( p->intinp() != nullptr ) p->intinp()->setEnabled(true);
            //if ( p->inp()    != nullptr ) p->inp()->setEnabled(true);
            //if ( p->inp2()   != nullptr ) p->inp2()->setEnabled(true);
            //if ( p->inp3()   != nullptr ) p->inp3()->setEnabled(true);
            //if ( p->tog()    != nullptr ) p->tog()->setEnabled(true);
            if ( p->fit()    != nullptr ) p->fit()->setEnabled(true);
            if ( p->lbl()    != nullptr ) p->lbl()->setEnabled(true);
        }*/
    }
}

/*static*/ QString myGuiParam::searchSubKey(QString key)
{
    foreach ( myGuiParam *p, allParams )
    {
        if ( key.compare(p->keyName(),Qt::CaseInsensitive)  == 0 ) return "0:" + p->keyName();
        if ( key.compare(p->keyName2(),Qt::CaseInsensitive) == 0 ) return "2:" + p->keyName2();
        if ( key.compare(p->keyName3(),Qt::CaseInsensitive) == 0 ) return "3:" + p->keyName3();
    }
    return "?";
}

/*static*/ myGuiParam *myGuiParam::setEnabled(QString key, bool ena, QString lbl)
{
    myGuiParam *gp = getGuiParam(key);
    if ( gp == nullptr )
    {
        qDebug() << "myGuiParam::setEnabled()  unknown key" << key << searchSubKey(key);
        return nullptr;
    }
    if ( gp->_par == nullptr )
    {
        paramHelper *par;
        par = SC_MainGUI::getInst()->getCalcGui()->params.value(gp->keyName(),nullptr);
        if ( par != nullptr )
            gp->_par = par;
        else if ( ! gp->keyName2().isEmpty() )
        {
            par = SC_MainGUI::getInst()->getCalcGui()->params.value(gp->keyName2(),nullptr);
            if ( par != nullptr )
                gp->_par = par;
            else if ( ! gp->keyName3().isEmpty() )
            {
                par = SC_MainGUI::getInst()->getCalcGui()->params.value(gp->keyName3(),nullptr);
                if ( par != nullptr ) gp->_par = par;
            }
        }
    }
    if ( gp->_par != nullptr )
    {
        //qDebug() << "setEnabled(" << key << ") =" << ena;
        gp->_par->enabled = ena;
    }
    else
        qDebug() << "setEnabled(" << key << ") _par unbekannt";
    if ( hideValues )
    {   // Jetzt sollen alle Werte versteckt werden (visible)
        if ( gp->cbs()    != nullptr ) { gp->cbs()->setVisible(ena);    gp->cbs()->setEnabled(ena); }
        if ( gp->intinp() != nullptr ) { gp->intinp()->setVisible(ena); gp->intinp()->setEnabled(ena); }
        if ( gp->inp()    != nullptr ) { gp->inp()->setVisible(ena);    gp->inp()->setEnabled(ena); }
        if ( gp->inp2()   != nullptr ) { gp->inp2()->setVisible(ena);   gp->inp2()->setEnabled(ena); }
        if ( gp->inp3()   != nullptr ) { gp->inp3()->setVisible(ena);   gp->inp3()->setEnabled(ena); }
        if ( gp->tog()    != nullptr ) { gp->tog()->setVisible(ena);    gp->tog()->setEnabled(ena); }
        if ( gp->fit()    != nullptr ) { gp->fit()->setVisible(ena);    gp->fit()->setEnabled(ena);  gp->setFitCheck(ena); }
        if ( gp->lbl()    != nullptr ) { gp->lbl()->setVisible(ena);    gp->lbl()->setEnabled(ena); }
    }
    else
    {   // Jetzt werden die Werte nur gesperrt, aber die Eingaben bleiben verwendbar
        //if ( gp->cbs()    != nullptr ) gp->cbs()->setEnabled(ena);
        //if ( gp->intinp() != nullptr ) gp->intinp()->setEnabled(ena);
        //if ( gp->inp()    != nullptr ) gp->inp()->setEnabled(ena);
        //if ( gp->inp2()   != nullptr ) gp->inp2()->setEnabled(ena);
        //if ( gp->inp3()   != nullptr ) gp->inp3()->setEnabled(ena);
        //if ( gp->tog()    != nullptr ) gp->tog()->setEnabled(ena);
        if ( gp->fit()    != nullptr ) gp->fit()->setEnabled(ena);
        if ( gp->lbl()    != nullptr ) gp->lbl()->setEnabled(ena);
    }
    if ( !lbl.isEmpty() )
    {
        if ( lbl[0] == '@' && gp->cbs() != nullptr )
            gp->cbs()->setCurrentIndex( lbl.midRef(1).toInt() );
        else if ( gp->lbl() != nullptr )
            gp->lbl()->setText(lbl);
    }
    return gp;
}

/*static*/ void myGuiParam::setEnabled(QWidget *w, bool ena)
{
    if ( hideValues )
        w->setVisible(ena);
    else
        w->setEnabled(ena);
}

/*static*/ bool myGuiParam::isEnabled(QString key)
{
    myGuiParam *gp = getGuiParam(key);
    if ( gp == nullptr )
    {
        qDebug() << "myGuiParam::isEnabled()  unknown key" << key << searchSubKey(key);
        return false;
    }
    if ( gp->_par != nullptr )
        return gp->_par->enabled;
    qDebug() << "myGuiParam::isEnabled("<<key<<")  no _par set.";

    if ( hideValues )
    {   // Jetzt sollen alle Werte versteckt werden (visible)
        if ( gp->cbs()    != nullptr ) return gp->cbs()->isVisible();
        if ( gp->intinp() != nullptr ) return gp->intinp()->isVisible();
        if ( gp->inp()    != nullptr ) return gp->inp()->isVisible();
        if ( gp->inp2()   != nullptr ) return gp->inp2()->isVisible();
        if ( gp->inp3()   != nullptr ) return gp->inp3()->isVisible();
        if ( gp->tog()    != nullptr ) return gp->tog()->isVisible();
        if ( gp->fit()    != nullptr ) return gp->isFitEnabled(); // fit()->isVisible();  eigentlich nicht verwendet
        if ( gp->lbl()    != nullptr ) return gp->lbl()->isVisible();
    }
    else
    {   // Jetzt werden die Werte nur gesperrt, aber die Eingaben bleiben verwendbar
        //if ( gp->cbs()    != nullptr ) return gp->cbs()->isEnabled();
        //if ( gp->intinp() != nullptr ) return gp->intinp()->isEnabled();
        //if ( gp->inp()    != nullptr ) return gp->inp()->isEnabled();
        //if ( gp->inp2()   != nullptr ) return gp->inp2()->isEnabled();
        //if ( gp->inp3()   != nullptr ) return gp->inp3()->isEnabled();
        //if ( gp->tog()    != nullptr ) return gp->tog()->isEnabled();
        if ( gp->fit()    != nullptr ) return gp->fit()->isEnabled();
        if ( gp->lbl()    != nullptr ) return gp->lbl()->isEnabled();
    }
    return true;
}

#ifdef undef
/*static*/ bool myGuiParam::isEnabled(QWidget *w)
{
    if ( hideValues )
        return w->isVisible();
    else
        return w->isEnabled();
}
#endif

/*static*/ void myGuiParam::setLabel(QString key, QString txt)
{
    myGuiParam *gp = getGuiParam(key);
    if ( gp == nullptr )
    {
        qDebug() << "myGuiParam::setLabel()  unknown key" << key;
        return;
    }
    if ( gp->lbl() != nullptr ) gp->lbl()->setText(txt);
}
