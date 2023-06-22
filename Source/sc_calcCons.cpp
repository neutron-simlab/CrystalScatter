#include "sc_calcCons.h"

//#include <QDebug>  funktioniert nicht
#include <QSettings>
#include <iostream>

QHash<QString,Double3> SC_CalcCons::inpVectors;
QHash<QString,double>  SC_CalcCons::inpValues;
QHash<QString,double>  SC_CalcCons::inpSingleValueVectors;
QHash<QString,paramConsHelper*> SC_CalcCons::params;


#include "sc_calc_generic.h"


#define D(x) // x // Debuging... ACHTUNG: das qDebug() sollte durch std::cerr ersetzt werden!



/**
 * @brief SC_CalcCons::SC_CalcCons
 * Constructor, it instantiate calculation subclass.
 */
SC_CalcCons::SC_CalcCons()
{
    calcGeneric = new SC_Calc_GENERIC;
    QStringList meta = calcGeneric->guiLayoutNeu();
    for ( int m=0; m<meta.size(); m++ )
    {
        QStringList sl = meta[m].split(";");
        while ( sl.size() < 5 ) sl << " ";  // tooltip;default sind optional

        // TODO: Neue Struktur des guiLayout()

        paramConsHelper *e = new paramConsHelper;
        e->fitparam = true; // TODO ??? sl[3].mid(3,3) == "fit";
        e->key = sl[1].trimmed();
        e->value.number = 0;
        e->value.flag   = false;
        e->minNum       = 0;
        e->maxNum       = 0;
        QStringList typval = sl[2].split("|",Qt::SkipEmptyParts);
        switch ( sl[0][0].toLatin1() )
        {
        case 'C':   // ComboBox  : "...|...|..." mit den jeweiligen Werten
            e->type = paramConsHelper::select;
            if ( !sl[4].isEmpty() )
            {   // set default in local variable
                e->value.number = typval.indexOf(sl[4]);
            }
            e->maxNum = typval.size();
            break;
        case 'N':   // Zahlenwert: "frac|min|max|unit" mit Nachkommastellen und Grenzwerten (optional)
            e->type = paramConsHelper::number;
            e->minNum = ( typval.size() > 1 ) ? typval[1].toDouble() : -10000.0;
            e->maxNum = ( typval.size() > 2 ) ? typval[2].toDouble() :  10000.0;
            if ( !sl[4].isEmpty() )
            {   // set default in local variable
                e->value.number = sl[4].toDouble();
            }
            break;
        case 'I':   // Integer-Zahlenwert: "min|max|unit"
            e->type = paramConsHelper::number;
            e->minNum = ( typval.size() > 0 ) ? typval[0].toDouble() : -10000.0;
            e->maxNum = ( typval.size() > 1 ) ? typval[1].toDouble() :  10000.0;
            if ( !sl[4].isEmpty() )
            {   // set default in local variable
                e->value.number = sl[4].toDouble();
            }
            break;
        case 'T':   // CheckBox  : "tog"
            e->type = paramConsHelper::toggle;
            if ( !sl[4].isEmpty() )
            {   // set default in local variable
                e->value.flag = sl[4].toInt() != 0;
            }
            e->maxNum = 1;
            break;
        default:
            // Outputs werden hier ignoriert
            e->key = "";
            break;
        }
        if ( !e->key.isEmpty() ) params.insert( e->key, e );
    }
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
    sets.beginGroup( calcGeneric->methodName() );
    QHash<QString,paramConsHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
    {
        paramConsHelper *par = ip.value();
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

        switch ( par->type )
        {
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
    // Globale Inputs für die Berechnungen
    //QHash<QString,Double3> inpVectors;
    //QHash<QString,double>  inpValues;
    //QHash<QString,double>  inpSingleValueVectors;
    // dann noch alle gemeinsamen Parameter
    sets.beginGroup( "Inputs" );
    int gp = sets.value( "EditGridPoints", -1 ).toInt();
    if ( gp < 0 )
        inpValues.insert( "GridPoints", sets.value( "GridPoints", 100 ).toInt() );
    else
        inpValues.insert( "GridPoints", gp );
    int hkl = sets.value( "Edithklmax", -1 ).toInt();
    if ( hkl < 0 )
        inpValues.insert( "HKLmax", sets.value( "HKLmax", 3 ).toInt() );
    else
        inpValues.insert( "HKLmax", sets.value( "Edithklmax", 3 ).toInt() );
    inpValues.insert( "RadioButtonQ1",  sets.value( "RadioButtonQ1", false ).toBool() );
    inpValues.insert( "RadioButtonQ2",  sets.value( "RadioButtonQ2", false ).toBool() );
    if ( sets.value( "RadioButtonQ4", false ).toBool() )
    {
        inpValues.insert( "RadioButtonQ4",  true );
        inpValues.insert( "ExpandImage", false );
    }
    else
    {
        inpValues.insert( "RadioButtonQ4",  false );
        inpValues.insert( "ExpandImage", sets.value( "ExpandImage", false ).toBool() );
    }
    // sets.value( "Threads", 0 ) -> per CmdLine Parameter
    inpValues.insert( "BeamPosX", sets.value( "EditCenterX", 0 ).toInt() );
    inpValues.insert( "BeamPosY", sets.value( "EditCenterY", 0 ).toInt() );

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
    inpSingleValueVectors.insert( "Editdom1",    inpVectors["SigXYZ"].x() );
    inpSingleValueVectors.insert( "Editdom2",    inpVectors["SigXYZ"].y() );
    inpSingleValueVectors.insert( "Editdom3",    inpVectors["SigXYZ"].z() );
}


void SC_CalcCons::saveParameter( QString fn )
{
    QSettings sets( fn, QSettings::IniFormat );
    sets.beginGroup( calcGeneric->methodName() );
    sets.remove(""); // Remove all previous keys in this group to save only the current settings
    QHash<QString,paramConsHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
    {
        paramConsHelper *par = ip.value();
        switch ( par->type )
        {
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
    // Globale Inputs für die Berechnungen
    //QHash<QString,Double3> inpVectors;
    //QHash<QString,double>  inpValues;
    //QHash<QString,double>  inpSingleValueVectors;
    // dann noch alle gemeinsamen Parameter
    sets.beginGroup( "Inputs" );
    sets.setValue( "GridPoints", inpValues["GridPoints"] );
    sets.setValue( "HKLmax",     inpValues["HKLmax"] );
    sets.setValue( "RadioButtonQ1",  inpValues["RadioButtonQ1"] );
    sets.setValue( "RadioButtonQ2",  inpValues["RadioButtonQ2"] );
    sets.setValue( "RadioButtonQ4",  inpValues["RadioButtonQ4"] );
    sets.setValue( "ExpandImage",    inpValues["ExpandImage"] );
    sets.setValue( "EditCenterX",    inpValues["BeamPosX"] );
    sets.setValue( "EditCenterY",    inpValues["BeamPosY"] );
    //sets.setValue( "Threads",        ui->inpNumCores->value() );  TODO?

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



QStringList SC_CalcCons::paramsForMethod( bool num, bool glob, bool fit )
{
    QStringList rv;
    QHash<QString,paramConsHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
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


double SC_CalcCons::currentParamValue( QString p )
{
    if ( params.contains(p) )
    {
        paramConsHelper *par = params.value(p);
        switch ( par->type )
        {
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


bool SC_CalcCons::limitsOfParamValue( QString p, double &min, double &max, bool &countable )
{
    min = 0;
    max = 0;
    countable = false;
    if ( params.contains(p) )
    {
        paramConsHelper *par = params.value(p);
        switch ( par->type )
        {
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


bool SC_CalcCons::updateParamValue( QString p, double v )
{
    if ( params.contains(p) )
    {
        paramConsHelper *par = params.value(p);
        switch ( par->type )
        {
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


bool SC_CalcCons::isCurrentParameterValid( QString p )
{
    // Methodenspezifischer Parameter?
    if ( params.contains(p) ) return true;
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
 * @param getData - if true call the method specific prepare function (not in fit)
 */
void SC_CalcCons::prepareCalculation( bool fromFit )
{
    if ( fromFit )
    {   // Special call during 2D-Fit to update the parameters
        calcGeneric->prepareData( &dataGetter );
        return;
    }
    calcGeneric->prepareData( &dataGetter );
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
    if ( params.contains(p) )
    {
        paramConsHelper *par = params[p];
        switch ( par->type )
        {
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
void SC_CalcCons::doCalculation( int numThreads )
{
    calcGeneric->doCalculation( numThreads );
}

double SC_CalcCons::doFitCalculation( int numThreads, int bstop, int border, long &cnt, long &nancnt )
{
    return calcGeneric->doFitCalculation( numThreads, bstop, border, cnt, nancnt );
}



/**
 * @brief SC_CalcCons::higResTimerElapsed
 * @return duration of the last calculation in milliseconds
 */
double SC_CalcCons::higResTimerElapsed( whichHigResTimer f )
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





#ifdef undef
/**
 * @brief calcHelper::calcHelper
 * @param c   - Calculation class pointer
 */
calcConsHelper::calcConsHelper(SC_Calc_GENERIC *c )
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
#endif
