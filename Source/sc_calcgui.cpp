#include "sc_calcgui.h"

#include <QDebug>
#include <QWidget>
#include <QSettings>
#include <QFrame>
#include <QAbstractItemView>
#include "sc_maingui.h"
#include "sc_calc_generic.h"


QHash<QString,Double3> SC_CalcGUI::inpVectors;
QHash<QString,double>  SC_CalcGUI::inpValues;
QHash<QString,double>  SC_CalcGUI::inpSingleValueVectors;
QHash<QString,paramHelper*> SC_CalcGUI::params;


#define D(x)  //x // Debuging...


/**
 * @brief SC_CalcGUI::SC_CalcGUI
 * Constructor, it instantiate all calculation subclasses.
 */
SC_CalcGUI::SC_CalcGUI() : QObject()
{
    calcGeneric = new SC_Calc_GENERIC();

    QStringList meta = calcGeneric->guiLayoutNeu();
    for ( int m=0; m<meta.size(); m++ )
    {
        // Each element must be in the form "type;kenn;typespec;tooltip;default" with:
        //  type     = C:Selection, N:Double, I:Integer, T:Toggle, O:DoubleOut
        //  kenn     = internal parameter name to connect to correct gui element
        //  typespec = C:Selection : "...|...|..."  (required)
        //             N:Double    : "frac|min|max|unit"  (optional, default: "2|-10000|+10000|")
        //             I:Integer   : "min|max|unit"  (optional, default: "-10000|+10000|")
        //             T:Toggle    : (empty)
        //             O:DoubleOut : (empty)
        //  tooltip  = this is the tooltip set to both prompt label and inputfield (optional)
        //  default  = the default value (optional)
        QStringList sl = meta[m].split(";");
        paramHelper *e = new paramHelper;
        e->key = sl[1].trimmed();
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
        case 'O':
            e->type = paramHelper::outdbl;
            break;
        }
        params.insert( e->key, e );
    } // for m
}


/**
 * @brief SC_CalcGUI::saveParameter
 * @param fn - filename (ini)
 * This will save all parameters of all methods into the given file.
 */
void SC_CalcGUI::saveParameter( QString fn )
{
    QSettings sets( fn, QSettings::IniFormat );
    sets.beginGroup( calcGeneric->methodName() );
    sets.remove(""); // Remove all previous keys in this group to save only the current settings
    QHash<QString,paramHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
    {
        paramHelper *par = ip.value();
        switch ( par->type )
        {
        case paramHelper::undef:
        case paramHelper::outdbl:
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
QString SC_CalcGUI::loadParameter(QString fn, QString onlyMethod, bool &hkl, bool &grid)
{
    DT( qDebug() << "loadParameter()" << fn );
    QString rv="";
    hkl  = false;
    grid = false;
    QSettings sets( fn, QSettings::IniFormat );
    if ( !onlyMethod.isEmpty() && onlyMethod[0] == '@' )
        sets.beginGroup( onlyMethod.mid(1) );
    else
        sets.beginGroup( calcGeneric->methodName() );
    QStringList slKeys = sets.allKeys();
    QHash<QString,paramHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
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
    return rv;
}


void SC_CalcGUI::saveFitFlags( QSettings &sets )
{
    sets.beginGroup( calcGeneric->methodName() );
    QHash<QString,paramHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
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
}

void SC_CalcGUI::loadFitFlags( QSettings &sets )
{
    DT( qDebug() << "loadFitFlags()" );
        sets.beginGroup( calcGeneric->methodName() );
        QHash<QString,paramHelper*>::iterator ip = params.begin();
        while ( ip != params.end() )
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
    sets.beginGroup( calcGeneric->methodName() );
    QHash<QString,paramHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
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
    QStringList rv;
    QHash<QString,paramHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
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
    if ( params.contains(p) )
    {
        paramHelper *par = params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
        case paramHelper::outdbl:
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

    //QString m = curMethod->subCalc->methodName();

    //if ( p.contains("qmax",Qt::CaseInsensitive) )
    //    qDebug() << "curParDbl" << p << methods[m]->params.contains(p)
    //             << inpValues.contains(p) << inpSingleValueVectors.contains(p);

    if ( params.contains(p) )
    {
        paramHelper *par = params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
        case paramHelper::outdbl:
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
    if ( params.contains(p) )
    {
        paramHelper *par = params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
        case paramHelper::outdbl:
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
    if ( params.contains(p) )
    {
        paramHelper *par = params.value(p);
        switch ( par->type )
        {
        case paramHelper::undef:
        case paramHelper::outdbl:
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

    if ( params.contains(p) )
    {
        if ( dbg ) qDebug() << "updateParamValue: update" << p << v;
        paramHelper *par = params[p];

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
            if ( dbg ) qDebug() << "updateParamValue: unknown type" << p;
            break;
        }
    }
    if ( dbg ) qDebug() << "updateParamValue: unknown parameter" << p;
    return false;
}

bool SC_CalcGUI::updateParamValueColor(QString p, QColor col )
{
    DT( qDebug() << "updateParamValueColor()" << p );
    if ( ! params.contains(p) ) return false;

    paramHelper *par = params[p];
    if ( par->gui.w != nullptr && col != Qt::black )
    {
        SETCOL( par->gui.w, col );
        return true;
    }
    return false;
}

bool SC_CalcGUI::updateParamValueForFit( QString p, double v, bool dbg )
{
    if ( params.contains(p) )
    {
        paramHelper *par = params[p];
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
    QHash<QString,paramHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
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
}

/**
 * @brief SC_CalcGUI::allNumericalParams
 * @return a QHash with all numerical parameter values of this method.
 * If in GUI mode the internal values are updated before returning.
 */
SC_CalcGUI::_numericalParams SC_CalcGUI::allNumericalParams()
{
    _numericalParams rv;
    QHash<QString,paramHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
    {
        paramHelper *par = ip.value();
        switch ( par->type )
        {
        case paramHelper::undef:
        case paramHelper::outdbl:
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
    if ( params.contains(p) )
    {
        if ( forfit )
        {   // Die Fit-Parameter anders untersuchen
            if ( params[p]->togFit == nullptr ) return false;
            return params[p]->togFit->isEnabled() && params[p]->togFit->isChecked();
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
        calcGeneric->prepareData( &dataGetterForFit );
        return;
    }
    calcGeneric->prepareData( &dataGetter );
}

void SC_CalcGUI::updateOutputData()
{
    calcGeneric->updateOutputData( &dataSetter );
}


/**
 * @brief SC_CalcGUI::dataGetter [static]
 * @param p - parameter name to retrieve
 * @param v - value structure to return
 * @return true if the parameter was found, false in case of error
 * If in GUI mode the internal values are updated before returning.
 */
bool SC_CalcGUI::dataGetter( QString p, _valueTypes &v )
{
    DT( qDebug() << "dataGetter()" << p );

    if ( p.contains("qmax",Qt::CaseInsensitive) )
    {
        // Bzgl. der Eingabe von QMax gibt es zwei Wertefelder:
        // "G;0;EditQmax;", -> ist oben schon definiert
        // "G;0;CalcQmax;inp|4;Qmax calculated from data header above;",
        // Unterschieden werden die beiden Eingaben über die Radiobuttons
        //"G;0;EditQmaxData;tog;;",   // auch wenn das Radiobuttons sind
        //"G;0;EditQmaxPreset;tog;;", // -"-
        paramHelper *par = params["EditQmaxData"];
        if ( par->gui.w )
            par->value.flag = par->gui.rad->isChecked();  // RADIOBUTTON !!!
        if ( par->value.flag )
            p = "CalcQmax";
        else
            p = "EditQmax";
    }
    if ( params.contains(p) )
    {
        paramHelper *par = params[p];
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

bool SC_CalcGUI::dataSetter( QString p, _valueTypes &v )
{
    DT( qDebug() << "dataSetter()" << p );
    if ( ! params.contains(p) ) return false;

    // Ein Update von Werte zur Anzeige ist nur für bestimmte Daten verfügbar.

    paramHelper *par = params[p];
    if ( par->gui.w == nullptr ) return false;
    if ( par->type == paramHelper::outdbl )
    {
        //par->gui.out->setText( QString("%1").arg(v.value,0,'f',4) );
        par->gui.out->setText( QString("%1").arg(v.value) );
        return true;
    }
    qDebug() << "dataSetter:" << p << "not found";
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
    if ( p.contains("qmax",Qt::CaseInsensitive) )
    {
        //qDebug() << "dataGetter" << p;
        // Bzgl. der Eingabe von QMax gibt es zwei Wertefelder:
        // "G;0;EditQmax;", -> ist oben schon definiert
        // "G;0;CalcQmax;inp|4;Qmax calculated from data header above;",
        // Unterschieden werden die beiden Eingaben über die Radiobuttons
        //"G;0;EditQmaxData;tog;;",   // auch wenn das Radiobuttons sind
        //"G;0;EditQmaxPreset;tog;;", // -"-
        paramHelper *par = params["EditQmaxData"];
        //if ( par->gui.w )
        //    par->value.flag = par->gui.rad->isChecked();  // RADIOBUTTON !!!
        if ( par->value.flag )
            p = "CalcQmax";
        else
            p = "EditQmax";
        //qDebug() << "dataGetterForFit" << p << par->value.flag;
    }

    if ( params.contains(p) )
    {
        paramHelper *par = params[p];
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
 * Starts the calculation of the current method.
 */
void SC_CalcGUI::doCalculation( int numThreads )
{
    calcGeneric->doCalculation( numThreads );
}

/**
 * @brief SC_CalcGUI::doCalculation
 * @param numThreads - number of used threads (ignored in GPU mode)
 * Starts the calculation of the current method.
 */
double SC_CalcGUI::doFitCalculation( int numThreads, int bstop, int border, long &cnt, long &nancnt )
{
    return calcGeneric->doFitCalculation( numThreads, bstop, border, cnt, nancnt );
}


void SC_CalcGUI::setNoFitRect( int id, int x0, int y0, int x1, int y1 )
{
    if ( x0>=0 || y0>=0 || x1>0 || y1>0 )
        qDebug() << " SC_CalcGUI::setNoFitRect" << id << x0 << y0 << x1 << y1;
    calcGeneric->setNoFitRect( id, x0, y0, x1, y1 );
}


/**
 * @brief SC_CalcGUI::higResTimerElapsed
 * @return duration of the last calculation in milliseconds
 */
double SC_CalcGUI::higResTimerElapsed( whichHigResTimer f )
{
    switch ( f )
    {
    case htimPrep:
        return calcGeneric->higResTimerElapsedPrep();
    case htimCalc:
        return calcGeneric->higResTimerElapsedCalc();
    case htimBoth:
        return calcGeneric->higResTimerElapsedPrep() +
               calcGeneric->higResTimerElapsedCalc();
    }
    return 0;
}

bool SC_CalcGUI::getLastXY( int &x, int &y )
{
    if ( calcGeneric == nullptr ) return false;
    x = calcGeneric->lastX();
    y = calcGeneric->lastY();
    return true;
}

paramHelper *SC_CalcGUI::getParamPtr( QString p )
{
    DT( qDebug() << "getParamPtr()" << p );
    return params.value(p,nullptr);
}
