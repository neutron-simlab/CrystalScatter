#include "sc_calcCons.h"

//#include <QDebug>  funktioniert nicht
#include <QSettings>
#include <iostream>

QHash<QString,Double3> SC_CalcCons::inpVectors;
QHash<QString,double>  SC_CalcCons::inpValues;
QHash<QString,double>  SC_CalcCons::inpSingleValueVectors;
calcConsHelper *SC_CalcCons::curMethod;


// ADD NEW METHOD: First, add the include here, ...
#include "sc_calc_generic.h"
//#include "sc_calc_fcc.h"
//#include "sc_calc_bcc.h"  ist schon angefangen...
//#include "sc_calc_bct.h"
//#include "sc_calc_sc.h"
//#include "sc_calc_hcp.h"
//#include "sc_calc_fd3m.h"

#define D(x) // x // Debuging... ACHTUNG: das qDebug() sollte durch std::cerr ersetzt werden!



/**
 * @brief SC_CalcCons::SC_CalcCons
 * Constructor, it instantiate all calculation subclasses.
 */
SC_CalcCons::SC_CalcCons()
{
    memory = nullptr;

    // ADD NEW METHOD: ... and add the constructor call here - Thats all!
    calcConsHelper *c;
    c = new calcConsHelper( new SC_Calc_GENERIC );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcConsHelper( new SC_Calc_FCC  );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcConsHelper( new SC_Calc_BCC  );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcConsHelper( new SC_Calc_BCT  );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcConsHelper( new SC_Calc_SC   );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcConsHelper( new SC_Calc_HCP  );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcConsHelper( new SC_Calc_FD3M );   methods.insert( c->subCalc->methodName(), c );
}


/**
 * @brief SC_CalcCons::getCalcTypes
 * @return the names of the methods (aka calculation types)
 */
QStringList SC_CalcCons::getCalcTypes()
{
    QStringList rv = methods.keys();
    rv.sort();
    return rv;
}


/**
 * @brief SC_CalcCons::loadParameter
 * @param fn - filename (ini)
 * This will load all parameters of all methods from the given file.
 * There are no default values if a key is not found.
 * If in GUI mode, the input fields are updated.
 */
void SC_CalcCons::loadParameter( QString fn )
{
    QSettings sets( fn, QSettings::IniFormat );
    // Erst alle Methodenspezifischen Parameter
    QHash<QString,calcConsHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        sets.beginGroup( im.key() );
        QHash<QString,paramConsHelper*>::iterator ip = im.value()->params.begin();
        while ( ip != im.value()->params.end() )
        {
            paramConsHelper *par = ip.value();
            switch ( par->type )
            {
            //case paramConsHelper::text:
            //    par->value.text = sets.value( par->key ).toString();
            //    break;
            case paramConsHelper::number:
                par->value.number = sets.value( par->key ).toDouble();
                break;
            case paramConsHelper::select:
                par->value.number = sets.value( par->key ).toDouble();
                break;
            case paramConsHelper::toggle:
                par->value.flag = sets.value( par->key ).toBool();
                break;
            default:
                break;
            }
            ++ip;
        }
        sets.endGroup();
        ++im;
    }
    // Globale Inputs für die Berechnungen
    //QHash<QString,Double3> inpVectors;
    //QHash<QString,double>  inpValues;
    //QHash<QString,double>  inpSingleValueVectors;
    // dann noch alle gemeinsamen Parameter
    sets.beginGroup( "Inputs" );
    inpValues.insert( "EditGridPoints", sets.value( "EditGridPoints", 100 ).toInt() );
    inpValues.insert( "Edithklmax",     sets.value( "Edithklmax", 3 ).toInt() );
    inpValues.insert( "RadioButtonQ1",  sets.value( "RadioButtonQ1", false ).toBool() );
    inpValues.insert( "RadioButtonQ2",  sets.value( "RadioButtonQ2", false ).toBool() );
    inpValues.insert( "RadioButtonQ4",  sets.value( "RadioButtonQ4", false ).toBool() );
    // sets.value( "Threads", 0 ) -> per CmdLine Parameter
    // sets.value( "CurMethod", "" ) -> per CmdLine Parameter
    inpValues.insert( "ExpandImage", sets.value( "ExpandImage", false ).toBool() );
    inpValues.insert( "BeamPosX", sets.value( "EditCenterX", 0 ).toInt() );
    inpValues.insert( "BeamPosY", sets.value( "EditCenterY", 0 ).toInt() );

    // Lattice Werte
    inpValues.insert( "LATTcols", sets.value( "LATTcols", 0 ).toInt() );
    inpValues.insert( "LATTrows", sets.value( "LATTrows", 0 ).toInt() );
    inpValues.insert( "LATTcenx", sets.value( "LATTcenx", 0 ).toDouble() );
    inpValues.insert( "LATTceny", sets.value( "LATTceny", 0 ).toDouble() );
    inpValues.insert( "LATTwlen", sets.value( "LATTwlen", 0 ).toDouble() );
    inpValues.insert( "LATTdist", sets.value( "LATTdist", 0 ).toDouble() );
    inpValues.insert( "LATTpixx", sets.value( "LATTpixx", 0 ).toDouble() );
    inpValues.insert( "LATTpixy", sets.value( "LATTpixy", 0 ).toDouble() );

    inpVectors.insert( "Uvec", Double3(sets.value( "Editx1rel", 0.0 ).toDouble(),
                                       sets.value( "Edity1rel", 0.0 ).toDouble(),
                                       sets.value( "Editz1rel", 0.0 ).toDouble()) );
    inpVectors.insert( "Vvec", Double3(sets.value( "Editx2rel", 0.0 ).toDouble(),
                                       sets.value( "Edity2rel", 0.0 ).toDouble(),
                                       sets.value( "Editz2rel", 0.0 ).toDouble()) );
    inpVectors.insert( "Nvec", Double3(sets.value( "Editxrel", 0.0 ).toDouble(),
                                       sets.value( "Edityrel", 0.0 ).toDouble(),
                                       sets.value( "Editzrel", 0.0 ).toDouble()) );
    inpVectors.insert( "Ax1", Double3(sets.value( "EditAxis1x", 0.0 ).toDouble(),
                                      sets.value( "EditAxis1y", 0.0 ).toDouble(),
                                      sets.value( "EditAxis1z", 0.0 ).toDouble()) );
    inpVectors.insert( "Ax2", Double3(sets.value( "EditAxis2x", 0.0 ).toDouble(),
                                      sets.value( "EditAxis2y", 0.0 ).toDouble(),
                                      sets.value( "EditAxis2z", 0.0 ).toDouble()) );
    inpVectors.insert( "Ax3", Double3(sets.value( "EditAxis3x", 0.0 ).toDouble(),
                                      sets.value( "EditAxis3y", 0.0 ).toDouble(),
                                      sets.value( "EditAxis3z", 0.0 ).toDouble()) );
    inpVectors.insert( "SigXYZ", Double3(sets.value( "Editdom1", 0.0 ).toDouble(),
                                         sets.value( "Editdom2", 0.0 ).toDouble(),
                                         sets.value( "Editdom3", 0.0 ).toDouble()) );
    sets.endGroup();

    // In normal calculation modes the vector notation is used, but in the
    // AI mode all the values are accessed as simgle values.
    inpSingleValueVectors.insert( "EditAxis1x",  inpVectors["Ax1"].x() );
    inpSingleValueVectors.insert( "EditAxis1y",  inpVectors["Ax1"].y() );
    inpSingleValueVectors.insert( "EditAxis1z",  inpVectors["Ax1"].z() );
    inpSingleValueVectors.insert( "EditAxis2x",  inpVectors["Ax2"].x() );
    inpSingleValueVectors.insert( "EditAxis2y",  inpVectors["Ax2"].y() );
    inpSingleValueVectors.insert( "EditAxis2z",  inpVectors["Ax2"].z() );
    inpSingleValueVectors.insert( "EditAxis3x",  inpVectors["Ax3"].x() );
    inpSingleValueVectors.insert( "EditAxis3y",  inpVectors["Ax3"].y() );
    inpSingleValueVectors.insert( "EditAxis3z",  inpVectors["Ax3"].z() );
    inpSingleValueVectors.insert( "Editxrel",    inpVectors["Nvec"].x() );
    inpSingleValueVectors.insert( "Edityrel",    inpVectors["Nvec"].y() );
    inpSingleValueVectors.insert( "Editzrel",    inpVectors["Nvec"].z() );
    inpSingleValueVectors.insert( "Editx1rel",   inpVectors["Uvec"].x() );
    inpSingleValueVectors.insert( "Edity1rel",   inpVectors["Uvec"].y() );
    inpSingleValueVectors.insert( "Editz1rel",   inpVectors["Uvec"].z() );
    inpSingleValueVectors.insert( "Editx2rel",   inpVectors["Vvec"].x() );
    inpSingleValueVectors.insert( "Edity2rel",   inpVectors["Vvec"].y() );
    inpSingleValueVectors.insert( "Editz2rel",   inpVectors["Vvec"].z() );
    inpSingleValueVectors.insert( "Editdom1",    inpVectors["SigXYZ"].x() );
    inpSingleValueVectors.insert( "Editdom2",    inpVectors["SigXYZ"].y() );
    inpSingleValueVectors.insert( "Editdom3",    inpVectors["SigXYZ"].z() );
}


void SC_CalcCons::saveParameter( QString fn )
{
    QSettings sets( fn, QSettings::IniFormat );
    QHash<QString,calcConsHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        sets.beginGroup( im.key() );
        sets.remove(""); // Remove all previous keys in this group to save only the current settings
        QHash<QString,paramConsHelper*>::iterator ip = im.value()->params.begin();
        while ( ip != im.value()->params.end() )
        {
            paramConsHelper *par = ip.value();
            switch ( par->type )
            {
            //case paramConsHelper::text:
            //    sets.setValue( par->key, par->value.text /*string*/ );
            //    break;
            case paramConsHelper::number:
                sets.setValue( par->key, par->value.number /*double*/ );
                break;
            case paramConsHelper::select:
                sets.setValue( par->key, par->value.number /*double*/ );
                break;
            case paramConsHelper::toggle:
                sets.setValue( par->key, par->value.flag /*bool*/ );
                break;
            }
            ++ip;
        }
        sets.endGroup();
        ++im;
    }
    // Globale Inputs für die Berechnungen
    //QHash<QString,Double3> inpVectors;
    //QHash<QString,double>  inpValues;
    //QHash<QString,double>  inpSingleValueVectors;
    // dann noch alle gemeinsamen Parameter
    sets.beginGroup( "Inputs" );
    sets.setValue( "EditGridPoints", inpValues["EditGridPoints"] );
    sets.setValue( "Edithklmax",     inpValues["Edithklmax"] );
    sets.setValue( "RadioButtonQ1",  inpValues["RadioButtonQ1"] );
    sets.setValue( "RadioButtonQ2",  inpValues["RadioButtonQ2"] );
    sets.setValue( "RadioButtonQ4",  inpValues["RadioButtonQ4"] );
    sets.setValue( "CurMethod",      curMethod->subCalc->methodName() );
    sets.setValue( "ExpandImage",    inpValues["ExpandImage"] );
    sets.setValue( "EditCenterX",    inpValues["BeamPosX"] );
    sets.setValue( "EditCenterY",    inpValues["BeamPosY"] );

    sets.setValue( "Editx1rel",  inpVectors["Uvec"].x() );
    sets.setValue( "Edity1rel",  inpVectors["Uvec"].y() );
    sets.setValue( "Editz1rel",  inpVectors["Uvec"].z() );
    sets.setValue( "Editx2rel",  inpVectors["Vvec"].x() );
    sets.setValue( "Edity2rel",  inpVectors["Vvec"].y() );
    sets.setValue( "Editz2rel",  inpVectors["Vvec"].z() );
    sets.setValue( "Editxrel",   inpVectors["Nvec"].x() );
    sets.setValue( "Edityrel",   inpVectors["Nvec"].y() );
    sets.setValue( "Editzrel",   inpVectors["Nvec"].z() );
    sets.setValue( "EditAxis1x", inpVectors["Ax1"].x() );
    sets.setValue( "EditAxis1y", inpVectors["Ax1"].y() );
    sets.setValue( "EditAxis1z", inpVectors["Ax1"].z() );
    sets.setValue( "EditAxis2x", inpVectors["Ax2"].x() );
    sets.setValue( "EditAxis2y", inpVectors["Ax2"].y() );
    sets.setValue( "EditAxis2z", inpVectors["Ax2"].z() );
    sets.setValue( "EditAxis3x", inpVectors["Ax3"].x() );
    sets.setValue( "EditAxis3y", inpVectors["Ax3"].y() );
    sets.setValue( "EditAxis3z", inpVectors["Ax3"].z() );
    sets.setValue( "Editdom1",   inpVectors["SigXYZ"].x() );
    sets.setValue( "Editdom2",   inpVectors["SigXYZ"].y() );
    sets.setValue( "Editdom3",   inpVectors["SigXYZ"].z() );
    sets.endGroup();
}



QStringList SC_CalcCons::paramsForMethod( QString m, bool num, bool glob, bool fit )
{
    QStringList rv;
    if ( ! methods.contains(m) ) return rv;
    QHash<QString,paramConsHelper*>::iterator ip = methods[m]->params.begin();
    while ( ip != methods[m]->params.end() )
    {
        if ( fit )
        {   // Fit Parameter gehen hier vor
            if ( ip.value()->fitparam )
                rv << ip.value()->key;
        }
        else if ( !num || (ip.value()->type == paramConsHelper::number) )
            rv << ip.value()->key;
        ++ip;
    }
    if ( fit )
    {   // Add global fittable parameter (at the end)
        rv.sort();
        //-- QStringList tmp = inpSingleValueVectors.keys();
        //-- tmp.sort();
        //-- rv << tmp;  --> werden zuviele
        rv << "Editdom1" << "Editdom2" << "Editdom3";  // und mehr nicht.
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


double SC_CalcCons::currentParamValue( QString m, QString p )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) ) return 0;
    if ( methods[m]->params.contains(p) )
    {
        paramConsHelper *par = methods[m]->params.value(p);
        switch ( par->type )
        {
        //case paramConsHelper::text:
        //    return par->value.text;
        case paramConsHelper::number:
            return par->value.number;
        case paramConsHelper::select:
            return par->value.number;
        case paramConsHelper::toggle:
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


bool SC_CalcCons::limitsOfParamValue( QString m, QString p, double &min, double &max, bool &countable )
{
    min = 0;
    max = 0;
    countable = false;
    if ( ! methods.contains(m) ) return false;
    if ( methods[m]->params.contains(p) )
    {
        paramConsHelper *par = methods[m]->params.value(p);
        switch ( par->type )
        {
        //case paramConsHelper::text:
        //    return false;
        case paramConsHelper::number:
            min = par->minNum;
            max = par->maxNum;
            return true;
        case paramConsHelper::select:
            min = 0;
            max = par->maxNum;
            countable = true;
            return true;
        case paramConsHelper::toggle:
            min = 0;
            max = 1;
            countable = true;
            return true;
        }
        return false;
    }
    return false;
}


bool SC_CalcCons::updateParamValue( QString m, QString p, double v )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) ) return false;
    if ( methods[m]->params.contains(p) )
    {
        paramConsHelper *par = methods[m]->params.value(p);
        switch ( par->type )
        {
        //case paramConsHelper::text:
        //    par->value.text = v;
        //    return true;
        case paramConsHelper::number:
            par->value.number = v;
            return true;
        case paramConsHelper::select:
            par->value.number = v;
            return true;
        case paramConsHelper::toggle:
            par->value.flag = v != 0;
            return true;
        }
        return false;
    }
    // Globaler Wert?
    /*if ( inpValues.contains(p) )
    {
        inpValues[p] = Double3::
        return true;
    }*/
    // Globaler Vektor-Wert als Einzelwert?
    if ( inpSingleValueVectors.contains(p) )
    {
        inpSingleValueVectors[p] = v;
        return true;
    }
    return false;
}


bool SC_CalcCons::isCurrentParameterValid( QString m, QString p )
{
    // Methode gültig?
    if ( ! methods.contains(m) ) return false;
    // Methodenspezifischer Parameter?
    if ( methods[m]->params.contains(p) ) return true;
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
 * @brief SC_CalcCons::prepareCalculation
 * @param m       - current method to use ("*" during init)
 * @param getData - if true call the method specific prepare function (not in fit)
 */
void SC_CalcCons::prepareCalculation( QString m, bool getData )
{
    if ( m == "F" )
    {   // Special call during 2D-Fit to update the parameters
        curMethod->subCalc->prepareData( &dataGetter );
        return;
    }
    curMethod = methods.value(m,nullptr);
    memory = (curMethod!=nullptr) ? curMethod->subCalc : nullptr;
    if ( curMethod == nullptr )
        std::cerr << "ERROR " << qPrintable(m) << " " << qPrintable(methods.keys().join(","))
                  << "*********************" << std::endl;
    if ( getData )
        curMethod->subCalc->prepareData( &dataGetter );
}


/**
 * @brief SC_CalcCons::dataGetter [static]
 * @param p - parameter name to retrieve (method is set globally)
 * @param v - value structure to return
 * @return true if the parameter was found, false in case of error
 * If in GUI mode the internal values are updated before returning.
 */
bool SC_CalcCons::dataGetter( QString p, _valueTypes &v )
{
    if ( curMethod == nullptr ) return false;
    if ( curMethod->params.contains(p) )
    {
        paramConsHelper *par = curMethod->params[p];
        switch ( par->type )
        {
        //case paramConsHelper::text:
        //    v.str = par->value.text;
        //    D(qDebug() << "paramHelper::text" << p << v.str;)
        //    break;
        case paramConsHelper::number:
            v.value = par->value.number;
            D(qDebug() << "paramHelper::number" << p << v.value;)
            break;
        case paramConsHelper::select:
            v.select = par->value.number;
            D(qDebug() << "paramHelper::select" << p << v.select;)
            break;
        case paramConsHelper::toggle:
            v.checked = par->value.flag;
            D(qDebug() << "paramHelper::toggle" << p << v.checked;)
            break;
        default:
            return false;
        }
        return true;
    }
    if ( inpValues.contains(p) )
    {
        v.value = inpValues[p];
        D(qDebug() << "inpValues" << p << v.value;)
        return true;
    }
    if ( inpVectors.contains(p) )
    {
        v.vec = inpVectors[p];
        D(qDebug() << "inpVectors" << p << v.vec.toString();)
        return true;
    }
    std::cerr << "dataGetter: '" << qPrintable(p) << "' not found" << std::endl;
    return false;
}


/**
 * @brief SC_CalcCons::doCalculation
 * @param numThreads - number of used threads (ignored in GPU mode)
 * @param pa         - function to show the progress and get the abort flag
 * Starts the calculation of the current method.
 */
void SC_CalcCons::doCalculation( int numThreads, progressAndAbort pa )
{
    if ( curMethod == nullptr ) return;
    curMethod->subCalc->doCalculation( numThreads, pa );
}

double SC_CalcCons::doFitCalculation( int numThreads, int bstop, int border, long &cnt, long &nancnt )
{
    if ( curMethod == nullptr ) return 0.0;
    return curMethod->subCalc->doFitCalculation( numThreads, bstop, border, cnt, nancnt );
}



/**
 * @brief SC_CalcCons::higResTimerElapsed
 * @return duration of the last calculation in milliseconds
 */
double SC_CalcCons::higResTimerElapsed( whichHigResTimer f )
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






/**
 * @brief calcHelper::calcHelper
 * @param c   - Calculation class pointer
 */
calcConsHelper::calcConsHelper( SC_Calc *c )
{
    /* Daten aus dem Config-File zum Überschreiben der Werte aus dem Code:
    # [ <methode> ]
    #
    # Wenn der Parameter fittable ist gilt:
    #
    # <paramname> = <a> : <b> : <c>
    #    <a> ist der Defaultwert
    #    <b> ist der Minimalwert beim Fit
    #    <c> ist der Maximalwert beim Fit
    #
    # Wenn der Parameter nicht zum Fit geeignet ist gilt:
    #
    # <paramname> = <a>
    #    <a> ist der Defaultwert
    */

    // TODO

    subCalc = c;
    QStringList meta = c->guiLayout();
    for ( int m=0; m<meta.size(); m++ )
    {
        QStringList sl = meta[m].split(";");
        if ( sl.size() < 4 ) continue;  // x;y;prompt;type müssen immer da sein
        while ( sl.size() < 6 ) sl << " ";  // tooltip;default sind optional
        // "x;y;prompt;type;tooltip;default" with:  ==> sc_calc.h for details
        //  x;y     = index in the grid (0,1,2,....)
        //  prompt  = prompting text label left of the inputfield
        //  type    = Selection : "cbs|...|...|..."
        //             Fittable : "cbsfit|...|...|..."
        //obsolete:            Textinput : "txt|len"
        //            Numericals: "inp|frac|min|max|unit"
        //              Fittable: "inpfit|frac|min|max|unit"
        //            CheckBox  : "tog"
        //            Infolabel : "lbl"
        //  tooltip = this is the tooltip set to both prompt label and inputfield (optional)
        //  default = the default value (optional)
        paramConsHelper *e = new paramConsHelper;
        e->fitparam = sl[3].mid(3,3) == "fit";
        e->key = sl[2].trimmed();
        //e->value.text   = "";
        e->value.number = 0;
        e->value.flag   = false;
        e->minNum       = 0;
        e->maxNum       = 0;
        QStringList typval = sl[3].mid(e->fitparam?7:4).split("|",Qt::SkipEmptyParts);
        if ( sl[3].startsWith("cbs") )
        {   // ComboBox  : "cbs|...|...|..." mit den jeweiligen Werten
            e->type = paramConsHelper::select;
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.number = typval.indexOf(sl[5]);
            }
            e->maxNum = typval.size();
        }
        /*else if ( sl[3].startsWith("txt") )
        {   // Textinput : "txt|len" mit maximaler Textlänge (optional)
            e->type = paramConsHelper::text;
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.text = sl[5];
            }
        }*/
        else if ( sl[3].startsWith("inp") )
        {   // Zahlenwert: "inp|frac|min|max" mit Nachkommastellen und Grenzwerten (optional)
            e->type = paramConsHelper::number;
            if ( typval.size() > 1 )
                e->minNum = typval[1].toDouble();
            else
                e->minNum = -10000.0;
            if ( typval.size() > 2 )
                e->maxNum = typval[2].toDouble();
            else
                e->maxNum =  10000.0;
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.number = sl[5].toDouble();
            }
        }
        else if ( sl[3].startsWith("tog") )
        {   // CheckBox  : "tog"
            e->type = paramConsHelper::toggle;
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.flag = sl[5].toInt() != 0;
            }
            e->maxNum = 1;
        }
        else if ( sl[3].startsWith("lbl") )
        {   // Infolabel : "lbl" (hier wird der Prompt-Text über beide Spalten gezogen)
            // und es wird nicht in der Parameter-Struktur gespeichert.
            continue;
        }
        else
            continue;
        params.insert( e->key, e );
    }
}
