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
        // Each element must be in the form "type;kenn;typespec;tooltip;default" with:
        //  type     = C:Selection, N:Double, I:Integer, T:Toggle, O:DoubleOut
        //             Zweites Zeichen ist 'F' für Fittable oder '-' für nicht fittable
        //  kenn     = internal parameter name to connect to correct gui element
        //  typespec = C:Selection : "...|...|..."  (required)
        //             N:Double    : "frac|min|max|unit"  (optional, default: "2|-10000|+10000|")
        //             I:Integer   : "min|max|unit"  (optional, default: "-10000|+10000|")
        //             T:Toggle    : (empty)
        //             O:DoubleOut : (empty)
        //  tooltip  = this is the tooltip set to both prompt label and inputfield (optional)
        //  default  = the default value (optional)

        paramConsHelper *e = new paramConsHelper;
        e->fitparam = sl[0][1] == 'F';      // Dieses Flag wird in der GUI anders gespeichert, aber hier lassen wir es in dieser Struktur
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
            e->type = paramConsHelper::numdbl;
            e->minNum = ( typval.size() > 1 ) ? typval[1].toDouble() : -10000.0;
            e->maxNum = ( typval.size() > 2 ) ? typval[2].toDouble() :  10000.0;
            if ( !sl[4].isEmpty() )
            {   // set default in local variable
                e->value.number = sl[4].toDouble();
            }
            break;
        case 'I':   // Integer-Zahlenwert: "min|max|unit"
            e->type = paramConsHelper::numint;
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
    // Jetzt bessere Kopie aus der GUI
    QStringList gr = sets.childGroups();
    gr.removeOne("Inputs");
    gr.removeOne("FFT");
    gr.removeOne("AI");     // TODO: neue feste Gruppen hier ausblenden
    QString usedGroup = calcGeneric->methodName();
    if ( gr.size() > 1 )
    {   // Jetzt gibt es mehr als eine Methode (Alte Dateien)
        sets.beginGroup( "Inputs" );
        QString m = sets.value("CurMethod","").toString();
        sets.endGroup();
        if ( ! m.isEmpty() )
        {
            QString mm = m;
            if ( mm.indexOf(" ") > 0 ) mm.truncate(mm.indexOf(" "));
            foreach ( QString g, gr )
            {
                if ( g.startsWith(mm) )
                {
                    usedGroup = g;
                    break;
                }
            }
            std::cerr << "Load spec" << qPrintable(m) << std::endl;
        }
    }
    // Erst alle Methodenspezifischen Parameter
    sets.beginGroup( usedGroup );
    std::cerr << "Load param '" << qPrintable(usedGroup+"' from: "+fn) << std::endl;
    QStringList slKeys = sets.allKeys();
    QHash<QString,paramConsHelper*>::iterator ip = params.begin();
    while ( ip != params.end() )
    {
        paramConsHelper *par = ip.value();
        //std::cerr << "LoadParam: " << qPrintable(par->key+"="+sets.value(par->key).toString()) << std::endl;
        slKeys.removeOne(par->key);
        slKeys.removeOne("FITFLAG_"+par->key);
        if (par->key == "EditQmaxPreset" ||     // Die RadioButtons für die Qmax Auswahl, nur einen nutzen...
            par->key == "CenterMidpoint" )      // Radiobuttons für BeamCenter, nur einen nutzen...
        {
            ++ip;
            continue;
        }
        QString fitflag;
        switch ( par->type )
        {
        case paramConsHelper::numdbl:
            par->value.number = sets.value( par->key ).toDouble();
            if ( par->key == "EditPixelX" || par->key == "EditPixelY" )
            {   // In den "alten" Datensätzen waren diese Werte in Metern angegeben (also 0.0001m für 1mm)
                // daher sollten diese Werte jetzt angepasst werden, da hier Millimeter verwendet werden
                if ( par->value.number < 0.01 )
                    par->value.number *= 1000.0;
            }
            fitflag = sets.value("FITFLAG_"+par->key,"??").toString();
            break;
        case paramConsHelper::numint:
            par->value.number = sets.value( par->key ).toDouble();
            fitflag = sets.value("FITFLAG_"+par->key,"??").toString();
            break;
        case paramConsHelper::select:
            par->value.number = sets.value( par->key ).toDouble();
            fitflag = "??";
            break;
        case paramConsHelper::toggle:
            par->value.flag = sets.value( par->key ).toBool();
            fitflag = "??";
            break;
        default:
            break;
        }
        if ( fitflag != "??" )
        {
            par->fitparam = fitflag[0] != '0';
            slKeys.removeOne( "FITFLAG_"+par->key );
            //std::cerr << "LoadParam: " << qPrintable(par->key+" -> "+fitflag) << std::endl;
        }
        ++ip;
    }
    if ( slKeys.size() > 0 )
    {
        std::cerr << "LoadParams: ini keys not used: " << qPrintable(slKeys.join(", ")) << std::endl;
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

    // In der GUI werden Eingabefelder ein-/ausgeblendet, wenn die ComboBoxen geändert werden (per Callback).
    // Diese Funktion ist in der Konsolenanwendung nicht notwendig, da gesperrte Eingaben nicht von den
    // Berechnungen verwendet werden.
    // ABER: bei bestimmten Werten vom LType werden andere Parameter manipuliert. Und das wird hier gemacht.
    switch ( static_cast<int>(currentParamValue("LType")) )
    {
    case  0:  /*  Lamellae  */  //Z=36497
        updateParamValue( "CheckBoxTwinned", 0 );
        updateParamValue( "ComboboxParticle", 2 );
        updateParamValue( "Length", 250 );
        break;

    case  1:  /*  hexagonal cylinders  */  //Z=36544
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 1 );
        updateParamValue( "Length", 500 );
        break;

    case  2:  /*  square cylinders  */  //Z=36592
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 1 );
        updateParamValue( "Length", 500 );
        break;

    case  3:  /*  rectangular centered cylinders  */  //Z=36640
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 1 );
        updateParamValue( "Length", 500 );
        break;

    case  4:  /*  bcc  */  //Z=36688
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case  5:  /*  fcc  */  //Z=36735
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case  6:  /*  hcp  */  //Z=36782
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case  7:  /*  sc  */  //Z=36829
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case  8:  /*  tetragonal centered spheres  */  //Z=36876
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case 17:  /*  fd3m  */  //Z=36923
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case  9:  /*  gyroid  */  //Z=37018
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 3 );
        break;

    case 10:  /*  OBDD  */  //Z=37066
        {/*1*/  //Z=45878
        const double b = 0.5484;  //Z=45879
        const double c = 0.6200;  //Z=45880
        const double n = 2;  //Z=45881
        double a, Ri, Ro, phi; //, area, lp;  //Z=45883
        a = currentParamValue("uca");
        Ro = c*a;  //Z=45887
        Ri = b*Ro;  //Z=45888
        phi = n*4*M_PI*c*c*c*(1-b*b*b)/3.0;  //Z=45889
        updateParamValue( "EditRadius", Ri );  //Z=45893
        updateParamValue( "EditLength", Ro );  //Z=45894
        updateParamValue( "EditPhi", phi );  //Z=45895
        }/*1*/  //Z=45898
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 3 );
        break;

    case 11:  /*  Im3m  */  //Z=37114
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 3 );
        break;

    case 12:  /*  none  */  //Z=37162
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case 13:  /*  cpl  */  //Z=37206
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboBoxPeak", 7 );
        updateParamValue( "ComboboxParticle", 0 );
        updateParamValue( "Editdom1", 200 );  //Editdom1.Text = '200';
        updateParamValue( "Editdom2", 200 ); //Editdom2.Text = '200';
        updateParamValue( "Editdom3", 200 ); //Editdom3.Text = '200';
        inpVectors["SigXYZ"].setX( 200 );
        inpVectors["SigXYZ"].setY( 200 );
        inpVectors["SigXYZ"].setZ( 200 );
        break;

    case 14:  /*  2D-Hex-GiSAXS  */  //Z=37256
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case 15:  /*  2D-Square-GiSAXS  */  //Z=37314
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case 16:  /*  1D-Lam-GiSAXS  */  //Z=37372
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 2 );
        break;

    case 18:  /*  orthogonal centered spheres  */  //Z=37430
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboboxParticle", 0 );
        break;

    case 19:  /*  quasicrystal  */  //Z=37477
        updateParamValue( "CheckBoxTwinned", false );
        updateParamValue( "ComboBoxPeak", 7 );
        updateParamValue( "ComboboxParticle", 0 );
        updateParamValue( "Editdom1", 200 );  //Editdom1.Text = '200';
        updateParamValue( "Editdom2", 200 ); //Editdom2.Text = '200';
        updateParamValue( "Editdom3", 200 ); //Editdom3.Text = '200';
        inpVectors["SigXYZ"].setX( 200 );
        inpVectors["SigXYZ"].setY( 200 );
        inpVectors["SigXYZ"].setZ( 200 );
        break;

    } // switch currentParamValue("LType")
}


void SC_CalcCons::saveParameter( QString fn )
{
    QSettings sets( fn, QSettings::IniFormat );
    sets.beginGroup( calcGeneric->methodName() );
    sets.remove(""); // Remove all previous keys in this group to save only the current settings
    QHash<QString,paramConsHelper*>::iterator ip = params.begin();
    //std::cerr << "saveParameter: " << qPrintable(params.keys().join(", ")) << std::endl;
    while ( ip != params.end() )
    {
        paramConsHelper *par = ip.value();
        if ( par == nullptr )
        {
            //saveParameter: EditRelDis empty
            //saveParameter: EditDist empty
            //std::cerr << "saveParameter: " << qPrintable(ip.key()) << " empty" << std::endl;
            ++ip;
            continue;
        }
        switch ( par->type )
        {
        case paramConsHelper::numdbl:
        case paramConsHelper::numint:
            sets.setValue( par->key, par->value.number /*double*/ );
            break;
        case paramConsHelper::select:
            sets.setValue( par->key, par->value.number /*double*/ );
            break;
        case paramConsHelper::toggle:
            sets.setValue( par->key, par->value.flag /*bool*/ );
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
        else if ( !num ||
                 (ip.value()->type == paramConsHelper::numdbl) ||
                 (ip.value()->type == paramConsHelper::numint) )
            rv << ip.value()->key;
        ++ip;
    }
    if ( fit )
    {   // Add global fittable parameter (at the end)
        //-- QStringList tmp = inpSingleValueVectors.keys();
        //-- tmp.sort();
        //-- rv << tmp;  --> werden zuviele
        rv << "Editdom1" << "Editdom2" << "Editdom3";  // und mehr nicht.
        rv.sort();
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
        case paramConsHelper::numdbl:
        case paramConsHelper::numint:
            return par->value.number;
        case paramConsHelper::select:
            return par->value.number;
        case paramConsHelper::toggle:
            return par->value.flag;
        default:
            break;
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
        case paramConsHelper::numdbl:
        case paramConsHelper::numint:
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
        default:
            break;
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
        case paramConsHelper::numdbl:
        case paramConsHelper::numint:
            par->value.number = v;
            return true;
        case paramConsHelper::select:
            par->value.number = v;
            return true;
        case paramConsHelper::toggle:
            par->value.flag = v != 0;
            return true;
        default:
            break;
        }
        return false;
    }
    // Globaler Wert?
    if ( inpValues.contains(p) )
    {
        std::cerr << "UPDATE PARAM (3) " << qPrintable(p) << std::endl;
        //inpValues[p] = Double3::
        return true;
    }
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
void SC_CalcCons::dataGetter( QString p, _valueTypes &v )
{
    if ( params.contains(p) )
    {
        paramConsHelper *par = params[p];
        switch ( par->type )
        {
        case paramConsHelper::numdbl:
        case paramConsHelper::numint:
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
            break;
        }
        return;
    }
    if ( inpValues.contains(p) )
    {
        v.value = inpValues[p];
        D(qDebug() << "inpValues" << p << v.value;)
        return;
    }
    if ( inpVectors.contains(p) )
    {
        v.vec = inpVectors[p];
        D(qDebug() << "inpVectors" << p << v.vec.toString();)
        return;
    }
    std::cerr << "dataGetter: '" << qPrintable(p) << "' not found" << std::endl;
}


/**
 * @brief SC_CalcCons::doCalculation
 * @param numThreads - number of used threads (ignored in GPU mode)
 * @param pa         - function to show the progress and get the abort flag
 * Starts the calculation of the current method.
 */
void SC_CalcCons::doCalculation(int numThreads, bool ignNewSwitch)
{
    calcGeneric->doCalculation( numThreads, ignNewSwitch );
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
