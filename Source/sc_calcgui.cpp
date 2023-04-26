#include "sc_calcgui.h"

#include <QDebug>
#include <QWidget>
#include <QSettings>
#include <QFrame>
#include <QAbstractItemView>
#include "sc_maingui.h"
#include "sc_calc_generic.h"


calcHelper *SC_CalcGUI::curMethod = nullptr;
QHash<QString,Double3> SC_CalcGUI::inpVectors;
QHash<QString,double>  SC_CalcGUI::inpValues;
QHash<QString,double>  SC_CalcGUI::inpSingleValueVectors;


#define D(x)  //x // Debuging...


/**
 * @brief SC_CalcGUI::SC_CalcGUI
 * Constructor, it instantiate all calculation subclasses.
 */
SC_CalcGUI::SC_CalcGUI() : QObject()
{
    memory = nullptr;

    curMethod = new calcHelper( new SC_Calc_GENERIC() );

    if ( methodNamesSorted.size() == 0 )
        methodNamesSorted << curMethod->subCalc->methodName();
}

/**
 * @brief SC_CalcGUI::createTabs
 * @param tab - TabWidget witch will be filled with pages for each method.
 */
void SC_CalcGUI::createTabs( /*QTabWidget *tab,*/ QStatusBar *sb )
{
    mainStatusBar = sb;
    // Immer alphabetisch sortiert
    QStringList names = getCalcTypes();
    /*
    for ( int i=0; i<names.size(); i++ )
    {
        QWidget *wid = new QWidget;
        QFrame *frm = new QFrame(wid);
        frm->setObjectName(names[i]);
        frm->setFrameStyle( QFrame::Box );
        frm->setLayout( methods[names[i]]->subLayout );
        tab->addTab( frm, names[i] );
    }
    */
    // Reihenfolge eher zufällig
    /*
    QHash<QString,calcHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        QWidget *wid = new QWidget;
        wid->setLayout( im.value()->subLayout );
        tab->addTab( wid, im.key() );
        ++im;
    }
    */
    /*if ( names.size() == 1 )
    {
        //tab->tabBar()->hide();
        curMethod = methods.begin().value();
    }*/
}

/**
 * @brief SC_CalcGUI::getCalcTypes
 * @return the names of the methods (aka calculation types)
 */
QStringList SC_CalcGUI::getCalcTypes()
{
    return methodNamesSorted;
}

/**
 * @brief SC_CalcGUI::saveParameter
 * @param fn - filename (ini)
 * This will save all parameters of all methods into the given file.
 */
void SC_CalcGUI::saveParameter( QString fn )
{
    QSettings sets( fn, QSettings::IniFormat );
    //QHash<QString,calcHelper*>::iterator im = methods.begin();
    //while ( im != methods.end() )
    //{
        sets.beginGroup( curMethod->subCalc->methodName() );
        sets.remove(""); // Remove all previous keys in this group to save only the current settings
        QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
        while ( ip != curMethod->params.end() )
        {
            paramHelper *par = ip.value();
            switch ( par->type )
            {
            case paramHelper::undef:
                break;
            case paramHelper::numdbl:
                if ( par->gui.w )
                    par->value.number = par->gui.numd->value();
                sets.setValue( par->key, par->value.number /*double*/ );
                break;
            case paramHelper::numint:
                if ( par->gui.w )
                    par->value.number = par->gui.numi->value();
                sets.setValue( par->key, par->value.number /*double*/ );
                break;
            case paramHelper::select:
                if ( par->gui.w )
                    par->value.number = par->gui.cbs->currentIndex();
                sets.setValue( par->key, par->value.number /*double*/ );
                break;
            case paramHelper::toggle:
                if ( par->gui.w )
                    par->value.flag = par->gui.tog->isChecked();
                sets.setValue( par->key, par->value.flag /*bool*/ );
                break;
            }
            if ( par->gui.w != nullptr && par->togFit != nullptr )
            {
                QString s = "00";
                if ( par->togFit->isChecked() )
                {
                    s[0] = '1';
                    //s[1] = par->used4fit ? '1' : '0';
                }
                sets.setValue( "FITFLAG_"+par->key, s );
            }
            ++ip;
        }
        sets.endGroup();
        //++im;
    //}
}


/**
 * @brief SC_CalcGUI::loadParameter
 * @param fn         - filename (ini)
 * @param onlyMethod - method name to use
 * This will load from the given file all parameters:
 *   - of all known methods (onlyMethod=="")
 *   - for one specific known method (onlyMethod!="")
 *   - from a specific method in the file to the current method (onlyMetod="@...").
 * There are no default values if a key is not found.
 * If in GUI mode, the input fields are updated.
 */
QString SC_CalcGUI::loadParameter(QString fn, QString onlyMethod , bool &hkl, bool &grid)
{
    DT( qDebug() << "loadParameter()" << fn );
    QString rv="";
    hkl  = false;
    grid = false;
    QSettings sets( fn, QSettings::IniFormat );
    //QHash<QString,calcHelper*>::iterator im = methods.begin();
    //while ( im != methods.end() )
    //{
        if ( onlyMethod.isEmpty() /*|| im.key() == onlyMethod*/ || onlyMethod[0]=='@' )
        {
            if ( !onlyMethod.isEmpty() && onlyMethod[0] == '@' )
                sets.beginGroup( onlyMethod.mid(1) );
            else
                sets.beginGroup( curMethod->subCalc->methodName() );
            QStringList slKeys = sets.allKeys();
            QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
            while ( ip != curMethod->params.end() )
            {
                paramHelper *par = ip.value();
                slKeys.removeOne( par->key );

                //qDebug() << "LOAD" << par->key;
                if ( par->key == "EditQmaxPreset" || par->key == "EditQmaxData" )
                {   // Die RadioButtons für die Qmax Auswahl sind auch in dieser Liste....
                    ++ip;
                    continue;
                }
                if ( par->key.startsWith("FITFLAG_") )
                {   // Das Flag zum Fitten wird hier auch gespeichert...
                    ++ip;
                    continue;
                }

                //par->fitparam = false;
                switch ( par->type )
                {
                case paramHelper::numdbl:
                    // Damit bei unbekannten Schlüsseln die Daten nicht auf 0 gesetzt werden
                    par->value.number = sets.value( par->key, par->value.number ).toDouble();
                    if ( par->gui.w )
                        par->gui.numd->setValue( par->value.number );
                    break;
                case paramHelper::numint:
                    // Damit bei unbekannten Schlüsseln die Daten nicht auf 0 gesetzt werden
                    par->value.number = sets.value( par->key, par->value.number ).toDouble();
                    if ( par->gui.w )
                        par->gui.numi->setValue( par->value.number );
                    //if ( par->key.contains("HKLmex",Qt::CaseInsensitive) ) hkl=true;      TODO
                    //if ( par->key.contains("GridPoints",Qt::CaseInsensitive) ) grid=true;
                    break;
                case paramHelper::select:
                    par->value.number = sets.value( par->key, par->value.number ).toDouble();
                    if ( par->gui.w )
                        par->gui.cbs->setCurrentIndex( par->value.number );
                    //if ( par->key.contains("LType",Qt::CaseInsensitive) )
                    //    qDebug() << par->key << par->value.number;
                    break;
                case paramHelper::toggle:
                    par->value.flag = sets.value( par->key, par->value.flag ).toBool();
                    if ( par->gui.w )
                        par->gui.tog->setChecked( par->value.flag );
                    break;
                default:
                    break;
                }
                if ( par->gui.w != nullptr && par->togFit != nullptr )
                {
                    QString s = sets.value( "FITFLAG_"+par->key, "?" ).toString();
                    // Wenn kein Wert im INI-File vorhanden ist, dann bleibt der Status des Fit-
                    //  Toggles in der GUI erhalten.
                    if ( s[0] == '?' )
                        s = par->togFit->isChecked() ? "11" : "00";
                    par->togFit->setChecked(s[0]!='0');
                    //par->fitparam = s[1]!='0';

                }
                ++ip;
            }
            if ( slKeys.size() > 0 )
                rv += /*"Method: "+im.key()+*/"Unknown keys in file: "+slKeys.join(", ")+EOL;
            sets.endGroup();
            //if ( !onlyMethod.isEmpty() && onlyMethod[0] == '@' ) break;
        } // if method found
        //++im;
    //} // while methods
    return rv;
}


void SC_CalcGUI::saveFitFlags( QSettings &sets )
{
    //QHash<QString,calcHelper*>::iterator im = methods.begin();
    //while ( im != methods.end() )
    //{
        sets.beginGroup( curMethod->subCalc->methodName() );
        QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
        while ( ip != curMethod->params.end() )
        {
            paramHelper *par = ip.value();
            if ( par->gui.w != nullptr && par->togFit != nullptr )
            {
                QString s = "00";
                if ( par->togFit->isChecked() )
                {
                    s[0] = '1';
                    //s[1] = par->used4fit ? '1' : '0';
                }
                sets.setValue( "FITFLAG_"+par->key, s );
            }
            ++ip;
        }
        sets.endGroup();
        //++im;
    //}
}

void SC_CalcGUI::loadFitFlags( QSettings &sets )
{
    DT( qDebug() << "loadFitFlags()" );
    //QHash<QString,calcHelper*>::iterator im = methods.begin();
    //while ( im != methods.end() )
    //{
        sets.beginGroup( curMethod->subCalc->methodName() );
        QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
        while ( ip != curMethod->params.end() )
        {
            paramHelper *par = ip.value();
            if ( par->key.startsWith("FITFLAG_") )
            {   // Das Flag zum Fitten wird hier auch gespeichert...
                ++ip;
                continue;
            }
            if ( par->gui.w != nullptr && par->togFit != nullptr )
            {
                QString s = sets.value( "FITFLAG_"+par->key, "00" ).toString();
                par->togFit->setChecked(s[0]!='0');
                //par->fitparam = s[1]!='0';
            }
            ++ip;
        }
        sets.endGroup();
        //++im;
    //} // while methods
}


void SC_CalcGUI::updateToolTipForCompare( QWidget *w, QString txt )
{
    QString tt = w->toolTip();
    int p = tt.indexOf("File:");
    if ( p > 1 )
    {
        tt.truncate(p+6);
        w->setToolTip( tt + txt );
    }
    else if ( tt.isEmpty() || p == 0 )
        w->setToolTip( "File: " + txt );
    else
        w->setToolTip( tt + "\nFile: " + txt );
}

void SC_CalcGUI::compareParameter( QSettings &sets, QHash<QString,_CompValues*> &compWerte )
{
    //QHash<QString,calcHelper*>::iterator im = methods.begin();
    //while ( im != methods.end() )
    //{
        sets.beginGroup( curMethod->subCalc->methodName() );
        QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
        while ( ip != curMethod->params.end() )
        {
            paramHelper *par = ip.value();
            switch ( par->type )
            {
            case paramHelper::numdbl:
            {
                double c = par->gui.numd->value();
                double p = sets.value( par->key ).toDouble();
                if ( fabs(c-p) > 1e-6 )
                {
                    _CompValues *cv = new _CompValues;
                    cv->cur = QString::number(c);
                    cv->par = QString::number(p);
                    compWerte.insert( ip.key(), cv );
                    if ( par->gui.w != nullptr )
                    {
                        SETCOL( par->gui.w, SETCOLMARK_PARDIFF );
                        updateToolTipForCompare( par->gui.w, cv->par );
                    }
                }
                break;
            }
            case paramHelper::numint:
            {
                int c = par->gui.numi->value();
                int p = sets.value( par->key ).toInt();
                if ( c != p )
                {
                    _CompValues *cv = new _CompValues;
                    cv->cur = QString::number(c);
                    cv->par = QString::number(p);
                    compWerte.insert( ip.key(), cv );
                    if ( par->gui.w != nullptr )
                    {
                        SETCOL( par->gui.w, SETCOLMARK_PARDIFF );
                        updateToolTipForCompare( par->gui.w, cv->par );
                    }
                }
                break;
            }
            case paramHelper::select:
            {
                int c = par->gui.cbs->currentIndex();
                int p = sets.value( par->key ).toInt();
                if ( c != p )
                {
                    _CompValues *cv = new _CompValues;
                    cv->cur = QString::number(c)+" ("+par->gui.cbs->currentText()+")";
                    cv->par = QString::number(p);
                    compWerte.insert( ip.key(), cv );
                    if ( par->gui.w != nullptr )
                    {
                        SETCOL( par->gui.w, SETCOLMARK_PARDIFF );
                        updateToolTipForCompare( par->gui.w, cv->par );
                    }
                }
                break;
            }
            case paramHelper::toggle:
            {
                bool c = par->gui.tog->isChecked();
                bool p = sets.value( par->key ).toBool();
                if ( c != p )
                {
                    _CompValues *cv = new _CompValues;
                    cv->cur = c ? "True" : "False";
                    cv->par = p ? "True" : "False";
                    compWerte.insert( ip.key(), cv );
                    if ( par->gui.w != nullptr )
                    {
                        SETCOL( par->gui.w, SETCOLMARK_PARDIFF );
                        updateToolTipForCompare( par->gui.w, cv->par );
                    }
                }
                break;
            }
            default:
                break;
            }
            ++ip;
        }
        sets.endGroup();
        //++im;
    //}
}


/**
 * @brief SC_CalcGUI::paramsForMethod
 * @param num  - if true, returns only numerical parameter names (no toggles, no selections)
 * @param glob - if true, add all global variable names to the return list
 * @param fit  - if true, return only fittable parameter names
 * @return Stringlist of all parameter names of this method (no global things)
 */
QStringList SC_CalcGUI::paramsForMethod(bool num, bool glob, bool fit )
{
    DT( qDebug() << "paramsForMethod()" << num << glob << fit );
    QString m = curMethod->subCalc->methodName();
    QStringList rv;
    QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
    while ( ip != curMethod->params.end() )
    {
        if ( fit )
        {   // Fit Parameter gehen hier vor
            if ( ip.value()->togFit != nullptr && ip.value()->togFit->isChecked() && ip.value()->togFit->isEnabled() )
            {
                //ip.value()->fitparam = true;
                rv << ip.value()->key;
            }
            //else
            //    ip.value()->fitparam = false;
        }
        else if ( !num || (ip.value()->type == paramHelper::numdbl) || (ip.value()->type == paramHelper::numint) )
            rv << ip.value()->key;
        ++ip;
    }
    if ( fit )
    {   // Add global fittable parameter (at the end)
        rv.sort();
        //-- QStringList tmp = inpSingleValueVectors.keys();
        //-- tmp.sort();
        //-- rv << tmp;  --> werden zuviele
        //rv << "Editdom1" << "Editdom2" << "Editdom3";  // sind ja schon in der Parameter-Liste...
        return rv;
    }
    if ( glob )
    {   // Add global parameters
        rv << inpValues.keys();
        rv << inpSingleValueVectors.keys();
    }
    rv.sort();
    return rv;
}

/**
 * @brief SC_CalcGUI::currentParamValue
 * @param m    - current metod to use
 * @param p    - parameter name to read
 * @param text - if true and in GUI mode then the ComboBox values are returned
 *               as text, else as index
 * @return String of the current value of this parameter, "?" if not known.
 */
QString SC_CalcGUI::currentParamValueStr(QString p, bool text )
{
    DT( qDebug() << "currentParamValueStr()" << p );
    QString m = curMethod->subCalc->methodName();
    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
            break;
        case paramHelper::numdbl:
            if ( par->gui.w )
                par->value.number = par->gui.numd->value();
            return QString::number(par->value.number);
        case paramHelper::numint:
            if ( par->gui.w )
                par->value.number = par->gui.numi->value();
            return QString::number(par->value.number);
        case paramHelper::select:
            if ( par->gui.w )
            {
                par->value.number = par->gui.cbs->currentIndex();
                if ( text )
                    return par->gui.cbs->currentText();
            }
            return QString::number(par->value.number);
        case paramHelper::toggle:
            if ( par->gui.w )
                par->value.flag = par->gui.tog->isChecked();
            return par->value.flag ? "True" : "False";
        }
        return "?2";
    }
    // Globaler Wert?
    if ( inpValues.contains(p) )
        return QString::number(inpValues[p]);
    // Globaler Vektor-Wert als Einzelwert?
    if ( inpSingleValueVectors.contains(p) )
        return QString::number(inpSingleValueVectors[p]);
    return "?3";
}
double SC_CalcGUI::currentParamValueDbl( QString p )
{
    DT( qDebug() << "currentParamValueDbl()" << p );
    // Bzgl. der Eingabe von QMax gibt es zwei Wertefelder:
    //  inpEditQmax (EditQmax) - kommt vom User
    //  outCalcQmax (CalcQmax) - kommt aus der HeaderData-Tabelle
    // Unterschieden werden die beiden Eingaben über die Radiobuttons
    //  radEditQmaxData (EditQmaxData)
    //  radEditQmaxPreset (EditQmaxPreset)

    QString m = curMethod->subCalc->methodName();

    //if ( p.contains("qmax",Qt::CaseInsensitive) )
    //    qDebug() << "curParDbl" << p << methods[m]->params.contains(p)
    //             << inpValues.contains(p) << inpSingleValueVectors.contains(p);

    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
            break;
        case paramHelper::numdbl:
            if ( par->gui.w )
                par->value.number = par->gui.numd->value();
            return par->value.number;
        case paramHelper::numint:
            if ( par->gui.w )
                par->value.number = par->gui.numi->value();
            return par->value.number;
        case paramHelper::select:
            if ( par->gui.w )
                par->value.number = par->gui.cbs->currentIndex();
            return par->value.number;
        case paramHelper::toggle:
            if ( par->gui.w )
                par->value.flag = par->gui.tog->isChecked();
            return par->value.flag;
        }
        return 0.0;
    }
    // Globaler Wert?
    if ( inpValues.contains(p) )
        return inpValues[p];
    // Globaler Vektor-Wert als Einzelwert?
    if ( inpSingleValueVectors.contains(p) )
        return inpSingleValueVectors[p];
    return 0.0;
}
int SC_CalcGUI::currentParamValueInt( QString p )
{
    DT( qDebug() << "currentParamValueInt()" << p );
    QString m = curMethod->subCalc->methodName();
    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
            break;
        case paramHelper::numdbl:
            if ( par->gui.w )
                par->value.number = par->gui.numd->value();
            return par->value.number;
        case paramHelper::numint:
            if ( par->gui.w )
                par->value.number = par->gui.numi->value();
            return par->value.number;
        case paramHelper::select:
            if ( par->gui.w )
                par->value.number = par->gui.cbs->currentIndex();
            return par->value.number;
        case paramHelper::toggle:
            if ( par->gui.w )
                par->value.flag = par->gui.tog->isChecked();
            return par->value.flag;
        }
        return 0;
    }
    // Globaler Wert?
    if ( inpValues.contains(p) )
        return inpValues[p];
    // Globaler Vektor-Wert als Einzelwert?
    if ( inpSingleValueVectors.contains(p) )
        return inpSingleValueVectors[p];
    return 0;
}

bool SC_CalcGUI::limitsOfParamValue(QString p, double &min, double &max, bool &countable )
{
    DT( qDebug() << "limitsOfParamValue()" << p );
    min = 0;
    max = 0;
    countable = false;
    QString m = curMethod->subCalc->methodName();
    //if ( ! methods.contains(m) ) return false;
    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
            return false;
        case paramHelper::numdbl:
            if ( par->gui.w )
            {
                min = par->gui.numd->minimum();
                max = par->gui.numd->maximum();
                return true;
            }
            return false;
        case paramHelper::numint:
            if ( par->gui.w )
            {
                min = par->gui.numi->minimum();
                max = par->gui.numi->maximum();
                return true;
            }
            return false;
        case paramHelper::select:
            if ( par->gui.w )
            {
                min = 0;
                max = par->gui.cbs->count();
                countable = true;
                return true;
            }
            return false;
        case paramHelper::toggle:
            min = 0;
            max = 1;
            countable = true;
            return true;
        }
        return false;
    }
    return false;
}

bool SC_CalcGUI::updateParamValue(QString p, double v, QColor col, bool dbg )
{
    DT( qDebug() << "updateParamValue()" << p );
    if ( p.contains("Grid",Qt::CaseInsensitive) ) dbg=true;
    if ( p.contains("BeamPos",Qt::CaseInsensitive) ) dbg=true;
    if ( p.contains("Beamcenter",Qt::CaseInsensitive) ) dbg=true;
    if ( p.contains("Debye",Qt::CaseInsensitive) ) dbg=true;

    QString m = curMethod->subCalc->methodName();
    //if ( ! methods.contains(m) )
    //{
    //    if ( dbg ) qDebug() << "updateParamValue: unknown method" << m;
    //    return false;
    //}

    if ( curMethod->params.contains(p) )
    {
        if ( dbg ) qDebug() << "updateParamValue: update" << m << p << v;
        paramHelper *par = curMethod->params[p];

        if ( par->gui.w != nullptr && col != SETCOLMARK_IGNORED )
        {
            SETCOL( par->gui.w, col );
            //QPalette pal = par->gui.num->palette();
            //pal.setBrush( QPalette::Base, col );
            //par->gui.num->setPalette(pal);
        }
        switch ( par->type )
        {
        case paramHelper::numdbl:
            par->value.number = v;
            if ( par->gui.w )
                par->gui.numd->setValue( par->value.number );
            return true;
        case paramHelper::numint:
            par->value.number = v;
            if ( par->gui.w )
                par->gui.numi->setValue( par->value.number );
            return true;
        case paramHelper::select:
            par->value.number = v;
            if ( par->gui.w )
                par->gui.cbs->setCurrentIndex( par->value.number );
            return true;
        case paramHelper::toggle:
            par->value.flag = v != 0;
            if ( par->gui.w )
                par->gui.tog->setChecked( par->value.flag );
            return true;
        default:
            if ( dbg ) qDebug() << "updateParamValue: unknown type" << m << p;
            break;
        }
    }
    if ( dbg ) qDebug() << "updateParamValue: unknown parameter" << m << p;
    return false;
}

bool SC_CalcGUI::updateParamValueColor(QString p, QColor col )
{
    DT( qDebug() << "updateParamValueColor()" << p );
    QString m = curMethod->subCalc->methodName();
    //if ( ! methods.contains(m) ) return false;

    if ( ! curMethod->params.contains(p) ) return false;

    paramHelper *par = curMethod->params[p];
    if ( par->gui.w != nullptr && col != Qt::black )
    {
        SETCOL( par->gui.w, col );
        return true;
    }
    return false;
}

bool SC_CalcGUI::updateParamValueForFit( QString p, double v, bool dbg )
{
    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params[p];
        switch ( par->type )
        {
        case paramHelper::numdbl:
        case paramHelper::numint:
            if ( dbg && par->value.number != v )
                qDebug() << "updateParamValueForFit:Number" << p << v << "old:" << par->value.number;
            par->value.number = v;
            return true;
        case paramHelper::select:
            if ( dbg && par->value.number != v )
                qDebug() << "updateParamValueForFit:Select" << p << v << "old:" << par->value.number;
            par->value.number = v;
            return true;
        case paramHelper::toggle:
            if ( dbg ) qDebug() << "updateParamValueForFit:Toggle" << p << v << "old:" << par->value.flag;
            par->value.flag = v != 0;
            return true;
        default:
            break;
        }
    }
    qDebug() << "updateParamValueForFit: unknown param" << p;
    return false;
}

void SC_CalcGUI::resetParamColorMarker( QColor col )
{
    DT( qDebug() << "resetParamColorMarker()" << col );
    //QHash<QString,calcHelper*>::iterator im = methods.begin();
    //while ( im != methods.end() )
    //{
        //qDebug() << "resetParamColorMarker im" << im.key();
        QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
        while ( ip != curMethod->params.end() )
        {
            //qDebug() << "resetParamColorMarker ip" << ip.key();
            if ( ip.value()->gui.w != nullptr )
            {
                SETCOL(ip.value()->gui.w, col);
                if ( col == SETCOLMARK_CLEARED )
                {
                    QString tt = ip.value()->gui.w->toolTip();
                    int p = tt.indexOf("File:");
                    if ( p >= 0 ) tt.truncate(p);
                    ip.value()->gui.w->setToolTip( tt.trimmed() );
                }
            }
            ++ip;
        }
        //++im;
    //}
    //qDebug() << "resetParamColorMarker ENDE";
}

/**
 * @brief SC_CalcGUI::allNumericalParams
 * @return a QHash with all numerical parameter values of this method.
 * If in GUI mode the internal values are updated before returning.
 */
SC_CalcGUI::_numericalParams SC_CalcGUI::allNumericalParams()
{
    _numericalParams rv;
    QHash<QString,paramHelper*>::iterator ip = curMethod->params.begin();
    while ( ip != curMethod->params.end() )
    {
        paramHelper *par = ip.value();
        switch ( par->type )
        {
        case paramHelper::undef:
            break;
        case paramHelper::toggle:
            if ( par->gui.w )
                par->value.flag = par->gui.tog->isChecked();
            rv.insert( par->key, par->value.flag );
            break;
        case paramHelper::numdbl:
            if ( par->gui.w )
                par->value.number = par->gui.numd->value();
            rv.insert( par->key, par->value.number );
            break;
        case paramHelper::numint:
            if ( par->gui.w )
                par->value.number = par->gui.numi->value();
            rv.insert( par->key, par->value.number );
            break;
        case paramHelper::select:
            if ( par->gui.w )
                par->value.number = par->gui.cbs->currentIndex();
            rv.insert( par->key, par->value.number );
            break;
        }
        ++ip;
    }

    // Globale Werte
    QHash<QString,double>::iterator id = inpValues.begin();
    while ( id != inpValues.end() )
    {
        rv.insert( id.key(), id.value() );
        ++id;
    }

    // Globale Vektor-Werte als Einzelwerte
    QHash<QString,double>::iterator iv = inpSingleValueVectors.begin();
    while ( iv != inpSingleValueVectors.end() )
    {
        rv.insert( iv.key(), iv.value() );
        ++iv;
    }

    return rv;
}

bool SC_CalcGUI::isCurrentParameterValid(QString p, bool forfit )
{
    // Methodenspezifischer Parameter?
    if ( curMethod->params.contains(p) )
    {
        if ( forfit )
        {   // Die Fit-Parameter anders untersuchen
            if ( curMethod->params[p]->togFit == nullptr ) return false;
            return curMethod->params[p]->togFit->isEnabled() && curMethod->params[p]->togFit->isChecked();
        }
        return true;
    }
    // Globaler Wert?
    if ( inpValues.contains(p) ) return true;
    // Globaler Vektor-Wert?
    if ( inpVectors.contains(p) ) return true;
    // Globaler Vektor-Wert als Einzelwert?
    if ( inpSingleValueVectors.contains(p) ) return true;
    // Nicht gefunden
    return false;
}


/**
 * @brief SC_CalcGUI::prepareCalculation
 * @param fromFit - if true called from the 2d-fit (use other dataGetter)
 * @param getData - if true call the method specific prepare function
 */
void SC_CalcGUI::prepareCalculation(bool fromFit)
{
    if ( fromFit )
    {   // Special call during 2D-Fit to update the parameters
        curMethod->subCalc->prepareData( &dataGetterForFit );
        return;
    }
    memory = curMethod->subCalc;
    curMethod->subCalc->prepareData( &dataGetter );
}

/**
 * @brief SC_CalcGUI::dataGetter [static]
 * @param p - parameter name to retrieve (method is set globally)
 * @param v - value structure to return
 * @return true if the parameter was found, false in case of error
 * If in GUI mode the internal values are updated before returning.
 */
bool SC_CalcGUI::dataGetter( QString p, _valueTypes &v )
{
    DT( qDebug() << "dataGetter()" << p );
    if ( curMethod == nullptr ) return false;

    if ( p.contains("qmax",Qt::CaseInsensitive) )
    {
        // Bzgl. der Eingabe von QMax gibt es zwei Wertefelder:
        // "G;0;EditQmax;", -> ist oben schon definiert
        // "G;0;CalcQmax;inp|4;Qmax calculated from data header above;",
        // Unterschieden werden die beiden Eingaben über die Radiobuttons
        //"G;0;EditQmaxData;tog;;",   // auch wenn das Radiobuttons sind
        //"G;0;EditQmaxPreset;tog;;", // -"-
        paramHelper *par = curMethod->params["EditQmaxData"];
        if ( par->gui.w )
            par->value.flag = par->gui.rad->isChecked();  // RADIOBUTTON !!!
        if ( par->value.flag )
            p = "CalcQmax";
        else
            p = "EditQmax";
    }
    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params[p];
        switch ( par->type )
        {
        case paramHelper::numdbl:
            if ( par->gui.w )
                par->value.number = par->gui.numd->value();
            v.value = par->value.number;
            D(qDebug() << "dataGetter() paramHelper::numdbl" << p << v.value;)
            break;
        case paramHelper::numint:
            if ( par->gui.w )
                par->value.number = par->gui.numi->value();
            v.value = par->value.number;
            D(qDebug() << "dataGetter() paramHelper::numint" << p << v.value;)
            break;
        case paramHelper::select:
            if ( par->gui.w )
                par->value.number = par->gui.cbs->currentIndex();
            v.select = par->value.number;
            D(qDebug() << "dataGetter() paramHelper::select" << p << v.select;)
            break;
        case paramHelper::toggle:
            if ( par->gui.w )
                par->value.flag = par->gui.tog->isChecked();
            v.checked = par->value.flag;
            D(qDebug() << "dataGetter() paramHelper::toggle" << p << v.checked;)
            break;
        default:
            return false;
        }
        return true;
    }
    if ( inpValues.contains(p) )
    {
        v.value = inpValues[p];
        D(qDebug() << "dataGetter() inpValues" << p << v.value;)
        return true;
    }
    if ( inpVectors.contains(p) )
    {
        v.vec = inpVectors[p];
        D(qDebug() << "dataGetter() inpVectors" << p << v.vec.toString();)
        return true;
    }
    qDebug() << "dataGetter:" << p << "not found";
    return false;
}

/**
 * @brief SC_CalcGUI::dataGetterForFit [static]
 * @param p - parameter name to retrieve (method is set globally)
 * @param v - value structure to return
 * @return true if the parameter was found, false in case of error
 * Same function as dataGetter except no GUI access (for 2D-Fit)
 */
bool SC_CalcGUI::dataGetterForFit( QString p, _valueTypes &v )
{
    if ( curMethod == nullptr ) return false;

    if ( p.contains("qmax",Qt::CaseInsensitive) )
    {
        //qDebug() << "dataGetter" << p;
        // Bzgl. der Eingabe von QMax gibt es zwei Wertefelder:
        // "G;0;EditQmax;", -> ist oben schon definiert
        // "G;0;CalcQmax;inp|4;Qmax calculated from data header above;",
        // Unterschieden werden die beiden Eingaben über die Radiobuttons
        //"G;0;EditQmaxData;tog;;",   // auch wenn das Radiobuttons sind
        //"G;0;EditQmaxPreset;tog;;", // -"-
        paramHelper *par = curMethod->params["EditQmaxData"];
        //if ( par->gui.w )
        //    par->value.flag = par->gui.rad->isChecked();  // RADIOBUTTON !!!
        if ( par->value.flag )
            p = "CalcQmax";
        else
            p = "EditQmax";
        //qDebug() << "dataGetterForFit" << p << par->value.flag;
    }

    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params[p];
        switch ( par->type )
        {
        case paramHelper::numdbl:
        case paramHelper::numint:
            v.value = par->value.number;
            D(qDebug() << "dataGetterForFit() paramHelper::number" << p << v.value;)
            return true;
        case paramHelper::select:
            v.select = par->value.number;
            D(qDebug() << "dataGetterForFit() paramHelper::select" << p << v.select;)
            return true;
        case paramHelper::toggle:
            v.checked = par->value.flag;
            D(qDebug() << "dataGetterForFit() paramHelper::toggle" << p << v.checked;)
            return true;
        default:
            return false;
        }
    }
    if ( inpValues.contains(p) )
    {
        v.value = inpValues[p];
        D(qDebug() << "dataGetterForFit() inpValues" << p << v.value;)
        return true;
    }
    if ( inpVectors.contains(p) )
    {
        v.vec = inpVectors[p];
        D(qDebug() << "dataGetterForFit() inpVectors" << p << v.vec.toString();)
        return true;
    }
    return false;
}


/**
 * @brief SC_CalcGUI::doCalculation
 * @param numThreads - number of used threads (ignored in GPU mode)
 * @param pa         - function to show the progress and get the abort flag
 * Starts the calculation of the current method.
 */
void SC_CalcGUI::doCalculation( int numThreads )
{
    if ( curMethod == nullptr ) return;
    curMethod->subCalc->doCalculation( numThreads );
}

/**
 * @brief SC_CalcGUI::doCalculation
 * @param numThreads - number of used threads (ignored in GPU mode)
 * @param pa         - function to show the progress and get the abort flag
 * Starts the calculation of the current method.
 */
double SC_CalcGUI::doFitCalculation( int numThreads, int bstop, int border, long &cnt, long &nancnt )
{
    if ( curMethod == nullptr ) return 0.0;
    return curMethod->subCalc->doFitCalculation( numThreads, bstop, border, cnt, nancnt );
}


void SC_CalcGUI::setNoFitRect( int id, int x0, int y0, int x1, int y1 )
{
    if ( curMethod == nullptr ) return;
    if ( x0>=0 || y0>=0 || x1>0 || y1>0 )
        qDebug() << " SC_CalcGUI::setNoFitRect" << id << x0 << y0 << x1 << y1;
    curMethod->subCalc->setNoFitRect( id, x0, y0, x1, y1 );
}


/**
 * @brief SC_CalcGUI::higResTimerElapsed
 * @return duration of the last calculation in milliseconds
 */
double SC_CalcGUI::higResTimerElapsed( whichHigResTimer f )
{
    if ( curMethod == nullptr ) return 0;
    switch ( f )
    {
    case htimPrep:
        return curMethod->subCalc->higResTimerElapsedPrep();
    case htimCalc:
        return curMethod->subCalc->higResTimerElapsedCalc();
    case htimBoth:
        return curMethod->subCalc->higResTimerElapsedPrep() +
               curMethod->subCalc->higResTimerElapsedCalc();
    }
    return 0;
}

bool SC_CalcGUI::getLastXY( int &x, int &y )
{
    if ( memory == nullptr ) return false;
    x = memory->lastX();
    y = memory->lastY();
    return true;
}

paramHelper *SC_CalcGUI::getParamPtr( QString p )
{
    DT( qDebug() << "getParamPtr()" << p );
    return curMethod->params.value(p,nullptr);
}





/**
 * @brief calcHelper::calcHelper
 * @param c   - Calculation class pointer
 */
calcHelper::calcHelper(SC_Calc_GENERIC *c)
{
    DT( qDebug() << "calcHelper()" << c->methodName() );
    subCalc = c;
    QStringList meta = c->guiLayout();
    for ( int m=0; m<meta.size(); m++ )
    {
        QStringList sl = meta[m].split(";");
        if ( sl.size() < 4 )
        {   // x;y;prompt;type müssen immer da sein
            qDebug() << "calcHelper: Ignore" << sl;
            continue;
        }
        //while ( sl.size() < 6 ) sl << " ";  // tooltip;default sind optional
        // "x;y;prompt;type;tooltip;default" with:  ==> sc_calc.h for details
        //  x;y     = index in the grid (0,1,2,....)
        //  prompt  = prompting text label left of the inputfield
        //  type    = Selection : "cbs|...|...|..."
        //            Numericals: "inp|frac|min|max|unit"
        //            CheckBox  : "tog"
        //            Infolabel : "lbl"
        //  tooltip = this is the tooltip set to both prompt label and inputfield (optional)
        //  default = the default value (optional)
        paramHelper *e = new paramHelper;
        e->key = sl[2].trimmed();
        e->value.number = 0;
        e->value.flag   = false;
        e->gui.w = nullptr;
        switch ( sl[0][0].toLatin1() )
        {
        case 'N':
            e->type = paramHelper::numdbl;
            break;
        case 'I':
            e->type = paramHelper::numint;
            break;
        case 'T':
            e->type = paramHelper::toggle;
            break;
        case 'C':
            e->type = paramHelper::select;
            break;
        }
        params.insert( e->key, e );
    } // for m
}
