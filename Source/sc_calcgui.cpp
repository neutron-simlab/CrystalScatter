#include "sc_calcgui.h"

#include <QDebug>
#include <QWidget>
#include <QSettings>
#include <QFrame>
#include <QAbstractItemView>
#include "sc_maingui.h"

calcHelper *SC_CalcGUI::curMethod = nullptr;
QHash<QString,Double3> SC_CalcGUI::inpVectors;
QHash<QString,double>  SC_CalcGUI::inpValues;
QHash<QString,double>  SC_CalcGUI::inpSingleValueVectors;


// ADD NEW METHOD: First, add the include here, ...
#include "sc_calc_generic.h"
//#include "sc_calc_fcc.h"
//#include "sc_calc_bcc.h"  ist schon angefangen...
//#include "sc_calc_bct.h"
//#include "sc_calc_sc.h"
//#include "sc_calc_hcp.h"
//#include "sc_calc_fd3m.h"

#define D(x)  //x // Debuging...

/* Es wird im Endausbau die folgenden Methoden geben, die Zeilennummern sind aus crystal3d1.pas(16.6./18.8.):
 Ab Zeile X/11204 wird aus 'ComboBoxLattice.ItemIndex' der passende Wert für 'ltype' bestimmt. Aus den
 dortigen Kommentaren und den Kommentaren bei den jeweiligen Startzeile erkenne ich die Methodenbezeichnung.
 10842/12023: Berechnungen für Lamellae(0). . . . . . . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 11159/12340: Berechnungen für No Structure(12). Hier per Testmode unterschiedliche Berechnungen. . ==> ?
 11749/12923: Berechnungen für HEX cylinders(1) . . . . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 12047/13221: Berechnungen für Square Cylinders(2). . . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 12339/13513: Berechnungen für Centered Rectangular Cylinders(3). . . . . . . . . . . . . . . . . . ==> ?
 12658/13832: Berechnungen für CP-Layers(13). . . . . . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 12962/14136: Berechnungen für Quasicrystal, dodecagonal(19). . . . . . . . . . . . . . . . . . . . ==> ?
 13279/14453: Berechnungen für BCC Spheres(4) . . . . . . . . . . . . . . . . . . . . . . . . . . . ==> i.A.
 13700/14874: Berechnungen für FCC Spheres(5) wenn kein Testlauf(=neue Berechnungen) aktiviert ist. ==> ok
 14391/15564: Berechnungen für Generic(30), ist auch beim Testlauf(=neue Berechnungen FCC)  . . . . ==> ok
 14703/16755: Berechnungen für SC Spheres(7). . . . . . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 15025/17077: Berechnungen für HCP Spheres(6) . . . . . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 15459/17510: Berechnungen für BCT Tetragonal Spheres(8) wenn kein Generic  . . . . . . . . . . . . ==> i.A.
 15787/17839: Berechnungen für Orthorhombic Spheres(18) wenn kein Generic . . . . . . . . . . . . . ==> ?
 16126/18178: Berechnungen für Fd3m Spheres = diamond lattice(17) wenn kein Generic . . . . . . . . ==> ?
 16608/18660: Berechnungen für Ia3d(9) wenn kein Generic  . . . . . . . . . . . . . . . . . . . . . ==> ?
 17235/19287: Berechnungen für Pn3m(10) wenn kein Generic . . . . . . . . . . . . . . . . . . . . . ==> ?
 17638/19690: Berechnungen für Im3m(11) wenn kein Generic . . . . . . . . . . . . . . . . . . . . . ==> ?
 18046/20098: Berechnungen für 2D-Hex, GISAXS(14) . . . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 18276/20327: Berechnungen für 2D-Square, GISAXS(15). . . . . . . . . . . . . . . . . . . . . . . . ==> ?
 18471/20523: Berechnungen für 1D-Lam, GISAXS(16) . . . . . . . . . . . . . . . . . . . . . . . . . ==> ?

 Achtung: bei 0 wird das lat1d-Flag gesetzt,
          bei 1,2,3 wird das lat2d-Flag gesetzt,
          sonst ist das lat3d-Flag gesetzt.
          Durch diese Flags gibt es andere Vorberechnungen!
          (ab Zeile 11078)
*/

QStringList calcHelper::slDisWidgets;


/**
 * @brief SC_CalcGUI::SC_CalcGUI
 * Constructor, it instantiate all calculation subclasses.
 */
SC_CalcGUI::SC_CalcGUI() : QObject()
{
    popupMenu = nullptr;
    memory = nullptr;

    // ADD NEW METHOD: ... and add the constructor call here - Thats all!
    calcHelper *c;
    c = new calcHelper( new SC_Calc_GENERIC(),  true );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcHelper( new SC_Calc_FCC(),  true );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcHelper( new SC_Calc_BCC(),  true );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcHelper( new SC_Calc_BCT(),  true );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcHelper( new SC_Calc_SC,   true );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcHelper( new SC_Calc_HCP,  true );   methods.insert( c->subCalc->methodName(), c );
    //c = new calcHelper( new SC_Calc_FD3M, true );   methods.insert( c->subCalc->methodName(), c );

    methodNamesSorted = methods.keys();
    methodNamesSorted.sort();
}

/**
 * @brief SC_CalcGUI::createTabs
 * @param tab - TabWidget witch will be filled with pages for each method.
 */
void SC_CalcGUI::createTabs( QTabWidget *tab, QStatusBar *sb )
{
    mainStatusBar = sb;
    // Immer alphabetisch sortiert
    QStringList names = getCalcTypes();
    for ( int i=0; i<names.size(); i++ )
    {
        QWidget *wid = new QWidget;
        QFrame *frm = new QFrame(wid);
        frm->setObjectName(names[i]);
        frm->setFrameStyle( QFrame::Box );
        frm->setLayout( methods[names[i]]->subLayout );
        frm->setContextMenuPolicy( Qt::CustomContextMenu );
        connect( frm, SIGNAL(customContextMenuRequested(QPoint)),
                 this, SLOT(customContextMenuRequested(QPoint)) );
        tab->addTab( frm, names[i] );
    }
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
    if ( names.size() == 1 )
    {
        tab->tabBar()->hide();
        curMethod = methods.begin().value();
    }
}

void SC_CalcGUI::customContextMenuRequested(const QPoint &pos)
{
    if ( popupMenu == nullptr )
    {
        popupMenu = new QMenu();
        actionSetDefault = popupMenu->addAction( "Set to default" );
        actionCopyToAll  = popupMenu->addAction( "Copy to all methods" );
    }
    QWidget *w = static_cast<QWidget*>(sender());
    QLabel *lbl = static_cast<QLabel*>(w->childAt(pos));
    if ( lbl == nullptr ) return;   // Nicht auf einem Label geklickt
    QAction *a = popupMenu->exec( w->mapToGlobal(pos) );
    if ( a == nullptr ) return;     // ESC oder daneben geklickt
    // w ist der Frame und w->objectName() ist der Methodenname
    calcHelper *c = methods[w->objectName()];
    if ( a == actionSetDefault )
    {
        if ( c->params.contains(lbl->text()) )
        {   // Da der Default-Wert nicht gespeichert wird, muss ich diesen hier suchen
            QStringList sl = c->subCalc->guiLayout();
            for ( int i=0; i<sl.size(); i++ )
            {
                if ( sl[i].indexOf(";"+lbl->text()+";") > 0 )
                {
                    QStringList v = sl[i].split(";");
                    while ( v.size() < 6 ) v << "";
                    paramHelper *p = c->params[lbl->text()];
                    switch ( p->type )
                    {
                    /*case paramHelper::text:
                        // Default ist leerer Text
                        p->gui.inp->setText( v[5] );
                        break;*/
                    case paramHelper::number:
                        if ( v[5].isEmpty() )
                            mainStatusBar->showMessage( "No default value known", 5000 );
                        else
                            p->gui.num->setValue( v[5].toDouble() );
                        break;
                    case paramHelper::select:
                        if ( v[5].isEmpty() )
                            p->gui.cbs->setCurrentIndex(0);
                        else
                            p->gui.cbs->setCurrentIndex( p->gui.cbs->findText(v[5]) );
                        break;
                    case paramHelper::toggle:
                        if ( v[5].isEmpty() )
                            p->gui.tog->setChecked(false);
                        else
                            p->gui.tog->setChecked( v[5].toInt() != 0 );
                        break;
                    }
                    return;
                }
            }
        }
    }
    else if ( a == actionCopyToAll )
    {
        QString val = currentParamValueStr( w->objectName(), lbl->text(), false );
        if ( ! setParamToAllMethods( lbl->text(), val.toDouble() ) )
            mainStatusBar->showMessage( "No numerical parameter, no copy allowed", 5000 );
    }
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
    QHash<QString,calcHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        sets.beginGroup( im.key() );
        sets.remove(""); // Remove all previous keys in this group to save only the current settings
        QHash<QString,paramHelper*>::iterator ip = im.value()->params.begin();
        while ( ip != im.value()->params.end() )
        {
            paramHelper *par = ip.value();
            switch ( par->type )
            {
            /*case paramHelper::text:
                if ( par->gui.w )
                    par->value.text = par->gui.inp->text();
                sets.setValue( par->key, par->value.text );
                break;*/
            case paramHelper::number:
                if ( par->gui.w )
                    par->value.number = par->gui.num->value();
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
            ++ip;
        }
        sets.endGroup();
        ++im;
    }
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
QString SC_CalcGUI::loadParameter( QString fn, QString onlyMethod )
{
    QString rv="";
    QSettings sets( fn, QSettings::IniFormat );
    QHash<QString,calcHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        if ( onlyMethod.isEmpty() || im.key() == onlyMethod || onlyMethod[0]=='@' )
        {
            if ( !onlyMethod.isEmpty() && onlyMethod[0] == '@' )
                sets.beginGroup( onlyMethod.mid(1) );
            else
                sets.beginGroup( im.key() );
            QStringList slKeys = sets.allKeys();
            QHash<QString,paramHelper*>::iterator ip = im.value()->params.begin();
            while ( ip != im.value()->params.end() )
            {
                paramHelper *par = ip.value();
                slKeys.removeOne( par->key );
                switch ( par->type )
                {
                /*case paramHelper::text:
                    par->value.text = sets.value( par->key ).toString();
                    if ( par->gui.w )
                        par->gui.inp->setText( par->value.text );
                    break;*/
                case paramHelper::number:
                    // Damit bei unbekannten Schlüsseln die Daten nicht auf 0 gesetzt werden
                    par->value.number = sets.value( par->key, par->value.number ).toDouble();
                    if ( par->gui.w )
                        par->gui.num->setValue( par->value.number );
                    break;
                case paramHelper::select:
                    par->value.number = sets.value( par->key, par->value.number ).toDouble();
                    if ( par->gui.w )
                        par->gui.cbs->setCurrentIndex( par->value.number );
                    if ( par->key.contains("LType",Qt::CaseInsensitive) )
                        qDebug() << par->key << par->value.number;
                    break;
                case paramHelper::toggle:
                    par->value.flag = sets.value( par->key, par->value.flag ).toBool();
                    if ( par->gui.w )
                        par->gui.tog->setChecked( par->value.flag );
                    break;
                default:
                    break;
                }
                ++ip;
            }
            if ( slKeys.size() > 0 )
                rv += "Method: "+im.key()+" Unknown keys in file: "+slKeys.join(", ")+EOL;
            sets.endGroup();
            if ( !onlyMethod.isEmpty() && onlyMethod[0] == '@' ) break;
        } // if method found
        ++im;
    } // while methods
    return rv;
}


void SC_CalcGUI::compareParameter( QSettings &sets, QHash<QString,_CompValues*> &compWerte )
{
    QHash<QString,calcHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        sets.beginGroup( im.key() );
        QHash<QString,paramHelper*>::iterator ip = im.value()->params.begin();
        while ( ip != im.value()->params.end() )
        {
            paramHelper *par = ip.value();
            switch ( par->type )
            {
            case paramHelper::number:
            {
                double c = par->gui.num->value();
                double p = sets.value( par->key ).toDouble();
                if ( fabs(c-p) > 1e-6 )
                {
                    _CompValues *cv = new _CompValues;
                    cv->cur = QString::number(c);
                    cv->par = QString::number(p);
                    compWerte.insert( im.key()+":"+ip.key(), cv );
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
                    compWerte.insert( im.key()+":"+ip.key(), cv );
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
                    compWerte.insert( im.key()+":"+ip.key(), cv );
                }
                break;
            }
            default:
                break;
            }
            ++ip;
        }
        sets.endGroup();
        ++im;
    }
}


/**
 * @brief SC_CalcGUI::setParamToAllMethods
 * @param key - parameter name
 * @param val - value to set
 * @return true if one or more parameters are set, false if nothing set
 * This will set the value to the given parameter name in all methods, if found.
 * If in GUI mode, the input fields are updated.
 */
bool SC_CalcGUI::setParamToAllMethods( QString key, double val )
{
    bool retval = false;
    QHash<QString,calcHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        if ( im.value()->params.contains(key) )
        {
            paramHelper *par = im.value()->params[key];
            if ( par->type == paramHelper::number )
            {
                par->value.number = val;
                if ( par->gui.w )
                    par->gui.num->setValue( val );
                retval = true;
            }
        }
        ++im;
    }
    return retval;
}

QStringList SC_CalcGUI::paramsForMethod( int im, bool num, bool glob, bool fit )
{
    return paramsForMethod( methodNamesSorted[im], num, glob, fit );
}

/**
 * @brief SC_CalcGUI::paramsForMethod
 * @param m    - current method to use
 * @param num  - if true, returns only numerical parameter names (no toggles, no selections)
 * @param glob - if true, add all global variable names to the return list
 * @param fit  - if true, return only fittable parameter names
 * @return Stringlist of all parameter names of this method (no global things)
 */
QStringList SC_CalcGUI::paramsForMethod( QString m, bool num, bool glob, bool fit )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    QStringList rv;
    if ( ! methods.contains(m) ) return rv;
    QHash<QString,paramHelper*>::iterator ip = methods[m]->params.begin();
    while ( ip != methods[m]->params.end() )
    {
        if ( fit )
        {   // Fit Parameter gehen hier vor
            if ( ip.value()->fitparam )
                rv << ip.value()->key;
        }
        else if ( !num || (ip.value()->type == paramHelper::number) )
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
QString SC_CalcGUI::currentParamValueStr( QString m, QString p, bool text )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) ) return "?1";
    if ( methods[m]->params.contains(p) )
    {
        paramHelper *par = methods[m]->params.value(p);
        switch ( par->type )
        {
        /*case paramHelper::text:
            if ( par->gui.w )
                par->value.text = par->gui.inp->text();
            return par->value.text;*/
        case paramHelper::number:
            if ( par->gui.w )
                par->value.number = par->gui.num->value();
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
double SC_CalcGUI::currentParamValueDbl( QString m, QString p )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) ) return 0.0;
    if ( methods[m]->params.contains(p) )
    {
        paramHelper *par = methods[m]->params.value(p);
        switch ( par->type )
        {
        case paramHelper::number:
            if ( par->gui.w )
                par->value.number = par->gui.num->value();
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
int SC_CalcGUI::currentParamValueInt( QString m, QString p )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) ) return 0;
    if ( methods[m]->params.contains(p) )
    {
        paramHelper *par = methods[m]->params.value(p);
        switch ( par->type )
        {
        case paramHelper::number:
            if ( par->gui.w )
                par->value.number = par->gui.num->value();
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
bool SC_CalcGUI::currentParamValueLog( QString m, QString p )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) ) return false;
    if ( methods[m]->params.contains(p) )
    {
        paramHelper *par = methods[m]->params.value(p);
        if ( par->type != paramHelper::toggle ) return false;
        if ( par->gui.w )
            par->value.flag = par->gui.tog->isChecked();
        return par->value.flag;
    }
    // Keine globalen Werte, da dort keine bool sind
    return false;
}

bool SC_CalcGUI::limitsOfParamValue( QString m, QString p, double &min, double &max, bool &countable )
{
    min = 0;
    max = 0;
    countable = false;
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) ) return false;
    if ( methods[m]->params.contains(p) )
    {
        paramHelper *par = methods[m]->params.value(p);
        switch ( par->type )
        {
        //case paramHelper::text:
        //    return false;
        case paramHelper::number:
            if ( par->gui.w )
            {
                min = par->gui.num->minimum();
                max = par->gui.num->maximum();
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

bool SC_CalcGUI::updateParamValue( QString m, QString p, double v, QColor col, bool dbg )
{
    if ( m.isEmpty() && curMethod != nullptr ) m = curMethod->subCalc->methodName();
    if ( ! methods.contains(m) )
    {
        if ( dbg ) qDebug() << "updateParamValue: unknown method" << m;
        return false;
    }

    if ( methods[m]->params.contains(p) )
    {
        if ( dbg ) qDebug() << "updateParamValue: update" << m << p << v;
        paramHelper *par = methods[m]->params[p];

        if ( par->gui.w != nullptr && col != Qt::black )
        {
            QPalette pal = par->gui.num->palette();
            pal.setBrush( QPalette::Base, col );
            par->gui.num->setPalette(pal);
        }
        switch ( par->type )
        {
        /*case paramHelper::text:
            par->value.text = v;
            if ( par->gui.w )
                par->gui.inp->setText( par->value.text );
            return true;*/
        case paramHelper::number:
            par->value.number = v;
            if ( par->gui.w )
                par->gui.num->setValue( par->value.number );
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

bool SC_CalcGUI::updateParamValueForFit( QString p, double v, bool dbg )
{
    QString m;
    if ( curMethod != nullptr )
        m = curMethod->subCalc->methodName();
    else
    {
        qDebug() << "updateParamValueForFit: method undef.";
        return false;
    }
    if ( ! methods.contains(m) )
    {
        qDebug() << "updateParamValueForFit: inv methos" << m;
        return false;
    }

    if ( methods[m]->params.contains(p) )
    {
        paramHelper *par = methods[m]->params[p];
        switch ( par->type )
        {
        /*case paramHelper::text:
            if ( dbg ) qDebug() << "updateParamValueForFit:Text" << p << v << "old:" << par->value.text;
            par->value.text = v;
            return true;*/
        case paramHelper::number:
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
    QHash<QString,calcHelper*>::iterator im = methods.begin();
    while ( im != methods.end() )
    {
        QHash<QString,paramHelper*>::iterator ip = methods[im.key()]->params.begin();
        while ( ip != methods[im.key()]->params.end() )
        {
            if ( ip.value()->gui.w != nullptr )
            {
                QPalette pal = ip.value()->gui.num->palette();
                pal.setBrush( QPalette::Base, col );
                ip.value()->gui.num->setPalette(pal);
            }
            ++ip;
        }
        ++im;
    }
}

/**
 * @brief SC_CalcGUI::allNumericalParams
 * @param m - current method to use
 * @return a QHash with all numerical parameter values of this method.
 * If in GUI mode the internal values are updated before returning.
 */
SC_CalcGUI::_numericalParams SC_CalcGUI::allNumericalParams( QString m )
{
    _numericalParams rv;
    if ( ! methods.contains(m) ) return rv;
    QHash<QString,paramHelper*>::iterator ip = methods[m]->params.begin();
    while ( ip != methods[m]->params.end() )
    {
        paramHelper *par = ip.value();
        switch ( par->type )
        {
        //case paramHelper::text:
        default:
            break; // TODO: text umwandeln in eine Zahl? - wohl eher nicht bzw. Textfelder kommen kaum vor
        case paramHelper::toggle:
            if ( par->gui.w )
                par->value.flag = par->gui.tog->isChecked();
            rv.insert( par->key, par->value.flag );
            break;
        case paramHelper::number:
            if ( par->gui.w )
                par->value.number = par->gui.num->value();
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

bool SC_CalcGUI::isCurrentParameterValid( QString m, QString p, bool forfit )
{
    // Methode gültig?
    if ( ! methods.contains(m) ) return false;
    // Methodenspezifischer Parameter?
    if ( methods[m]->params.contains(p) )
    {
        if ( methods[m]->params[p]->disabled ) return false;
        if ( forfit && !methods[m]->params[p]->fitparam ) return false;
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
 * @param m       - current method to use ("*" during init, "F" during 2D-Fit)
 * @param getData - if true call the method specific prepare function
 */
void SC_CalcGUI::prepareCalculation( QString m, bool getData )
{
    if ( m == "F" )
    {   // Special call during 2D-Fit to update the parameters
        curMethod->subCalc->prepareData( &dataGetterForFit );
        return;
    }
    curMethod = methods[m];
    if ( curMethod == nullptr )
    {
        memory = nullptr;
        return;
    }
    memory = curMethod->subCalc;
    if ( getData )
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
    if ( curMethod == nullptr ) return false;
    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params[p];
        switch ( par->type )
        {
        /*case paramHelper::text:
            if ( par->gui.w )
                par->value.text = par->gui.inp->text();
            v.str = par->value.text;
            D(qDebug() << "dataGetter() paramHelper::text" << p << v.str;)
            break;*/
        case paramHelper::number:
            if ( par->gui.w )
                par->value.number = par->gui.num->value();
            v.value = par->value.number;
            D(qDebug() << "dataGetter() paramHelper::number" << p << v.value;)
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
    if ( curMethod->params.contains(p) )
    {
        paramHelper *par = curMethod->params[p];
        switch ( par->type )
        {
        /*case paramHelper::text:
            v.str = par->value.text;
            return true;*/
        case paramHelper::number:
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
void SC_CalcGUI::doCalculation( int numThreads, progressAndAbort pa )
{
    if ( curMethod == nullptr ) return;
    //curMethod->subCalc->tpvPerformRandom(ids); // TODO
    curMethod->subCalc->doCalculation( numThreads, pa );
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

paramHelper *SC_CalcGUI::getParamPtr( QString m, QString p )
{
    if ( ! methods.contains(m) ) return nullptr;
    if ( methods[m]->params.contains(p) )
    {
        return methods[m]->params.value(p);
    }
    return nullptr;
}





/**
 * @brief calcHelper::calcHelper
 * @param c   - Calculation class pointer
 * @param gui - if true, create all GUI components
 */
calcHelper::calcHelper( SC_Calc *c, bool gui )
{
    //qDebug() << "calcHelper" << c->methodName() << gui;
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
    QSettings globsets(SETT_APP,SETT_GUI);
    QSettings *cfg = nullptr;
    if ( ! SC_MainGUI::configParamsFile.isEmpty() )
    {
        cfg = new QSettings( SC_MainGUI::configParamsFile, QSettings::IniFormat );
        cfg->beginGroup( c->methodName() );
        if ( cfg->allKeys().size() > 0 ) qDebug() << "calcHelper (ext config): " << cfg->group() << cfg->allKeys();
        globsets.beginGroup( "Fit-"+cfg->group() );
    }
    subCalc = c;
    QStringList meta = c->guiLayout();
    subLayout = new QGridLayout;
    subLayout->setVerticalSpacing(3);
    for ( int m=0; m<meta.size(); m++ )
    {
        QStringList sl = meta[m].split(";");
        bool parDis = false;
        if ( sl.size() < 4 )
        {   // x;y;prompt;type müssen immer da sein
            qDebug() << "calcHelper: Ignore" << sl;
            continue;
        }
        if ( sl[0].at(0) == 'X' )
        {   // Gesperrtes Feld
            if ( sl[1].at(0).isDigit() )
            {   // Normales Feld
                slDisWidgets << c->methodName()+":"+sl[2];
                sl[0] = sl[0].mid(1); // bleibt aber erhalten in der Struktur
                //qDebug() << "DISa" << sl;
                parDis = true;
            }
            else
            {   // Spezielle Einträge um globale Parameter sperren zu können
                sl.takeFirst();
                //qDebug() << "DISb" << sl;
                slDisWidgets << sl;
                continue;
            }
        }
        while ( sl.size() < 6 ) sl << " ";  // tooltip;default sind optional
        // "x;y;prompt;type;tooltip;default" with:  ==> sc_calc.h for details
        //  x;y     = index in the grid (0,1,2,....)
        //  prompt  = prompting text label left of the inputfield
        //  type    = Selection : "cbs|...|...|..."
        //             Fittable : "cbsfit|...|...|..."
        //            Textinput : "txt|len"
        //            Numericals: "inp|frac|min|max|unit"
        //              Fittable: "inpfit|frac|min|max|unit"
        //            CheckBox  : "tog"
        //            Infolabel : "lbl"
        //  tooltip = this is the tooltip set to both prompt label and inputfield (optional)
        //  default = the default value (optional)
        paramHelper *e = new paramHelper;
        e->fitparam = sl[3].mid(3,3) == "fit";
        e->disabled = parDis;
        e->key = sl[2].trimmed();
        //e->value.text   = "";
        e->value.number = 0;
        e->value.flag   = false;
        e->gui.w = nullptr;
        QStringList typval = sl[3].mid( e->fitparam ? 7 : 4 ).split("|",Qt::SkipEmptyParts);
        if ( sl[0] == "G" )
        {   // Spezielle Einträge für globale Werte zum Fitten
            e->lbl1 = nullptr;  // wird im Programm-Konstruktor gefüllt (ebenso wie .gui.w)
            e->type = paramHelper::number; // nur damit möglich
            params.insert( e->key, e );
            continue;
        }
        if ( cfg != nullptr && cfg->contains(sl[2]) )
        {
            globsets.remove( sl[2] );
            QStringList tmp = cfg->value(sl[2],"").toString().split(":");
            sl[5] = tmp[0];     // Default
            if ( tmp.size() > 1 )
            {   // Fittable
                e->fitparam = true;
                if ( typval.size() >= 3 )
                {
                    typval[2] = tmp[1]; // Min
                    typval[3] = tmp[2]; // Max
                }
                else
                {
                    if ( typval.size() == 0 ) typval << "3"; // Fraction
                    if ( typval.size() == 1 ) typval << tmp[1]; // Min
                    if ( typval.size() == 2 ) typval << tmp[2]; // Max
                }
            }
            else
            {   // Nicht fittable
                e->fitparam = false;
            }
        }
        if ( gui )
        {
            e->lbl1 = new QLabel( e->key );
            e->lbl1->setObjectName( "lbl_"+e->key );
            e->lbl1->setAlignment( Qt::AlignVCenter | Qt::AlignRight );
            e->lbl1->setContextMenuPolicy( Qt::NoContextMenu );
            if ( e->fitparam )
            {
                QFont fnt = e->lbl1->font();
                fnt.setBold(true);
                e->lbl1->setFont(fnt);
            }
        }
        if ( sl[3].startsWith("cbs") )
        {   // ComboBox  : "cbs|...|...|..." mit den jeweiligen Werten
            e->type = paramHelper::select;
            if ( gui )
            {
                e->gui.cbs = new QComboBox;
#ifdef CALC_INPUT_MAXWIDTH
                e->gui.cbs->setMaximumWidth(CALC_INPUT_MAXWIDTH);
#endif
                e->gui.cbs->view()->setMinimumWidth(200);
                /**/
                for ( int i=0; i<typval.size(); i++ ) typval[i] = typval[i]+QString(" {%1}").arg(i);
                /**/
                e->gui.cbs->addItems( typval );
                if ( !sl[5].isEmpty() )
                {   // Defaultwert setzen
                    e->gui.cbs->setCurrentIndex( e->gui.cbs->findText(sl[5],Qt::MatchStartsWith) );
                }
                e->gui.cbs->connect( e->gui.cbs, SIGNAL(currentIndexChanged(int)),
                                     SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
            }
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.number = typval.indexOf(sl[5]);
            }
        }
        /*else if ( sl[3].startsWith("txt") )
        {   // Textinput : "txt|len" mit maximaler Textlänge (optional)
            e->type = paramHelper::text;
            if ( gui )
            {
                e->gui.inp = new QLineEdit;
                if ( typval.size() > 0 )
                    e->gui.inp->setMaxLength( typval[0].toInt() );
                if ( !sl[5].isEmpty() )
                {   // Defaultwert setzen
                    e->gui.inp->setText( sl[5] );
                }
            }
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.text = sl[5];
            }
        }*/
        else if ( sl[3].startsWith("inp") )
        {   // Zahlenwert: "inp|frac|min|max|unit" mit Nachkommastellen, Grenzwerten und Einheit (optional)
            e->type = paramHelper::number;
            if ( gui )
            {
                e->gui.num = new QDoubleSpinBox;
#ifdef CALC_INPUT_MAXWIDTH
                e->gui.num->setMaximumWidth( CALC_INPUT_MAXWIDTH );
#endif
                //e->gui.num->setContentsMargins( 0, 0, 0, 0 ); geht nicht wirklich (TODO)
                if ( typval.size() > 0 )
                    e->gui.num->setDecimals( typval[0].toInt() );
                if ( typval.size() > 1 )
                    e->gui.num->setMinimum( typval[1].toDouble() );
                else
                    e->gui.num->setMinimum( -10000.0 );
                if ( typval.size() > 2 )
                    e->gui.num->setMaximum( typval[2].toDouble() );
                else
                    e->gui.num->setMaximum( 10000.0 );
                if ( typval.size() > 3 )
                {
                    e->gui.num->setSuffix( typval[3] );
                    if ( typval[3]=="°" && e->gui.num->minimum()==0 && e->gui.num->maximum()==360 )
                        e->gui.num->setWrapping(true); // Bei Winkeln aktiviert
                }
                if ( !sl[5].isEmpty() )
                {   // Defaultwert setzen
                    e->gui.num->setValue( sl[5].toDouble() );
                }
                e->gui.num->connect( e->gui.num, SIGNAL(valueChanged(double)),
                                     SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
            }
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.number = sl[5].toDouble();
            }
        }
        else if ( sl[3].startsWith("tog") )
        {   // CheckBox  : "tog"
            e->type = paramHelper::toggle;
            if ( gui )
            {
                e->gui.tog = new QCheckBox;
                e->gui.tog->connect( e->gui.tog, SIGNAL(stateChanged(int)),
                                     SC_MainGUI::getInst(), SLOT(automaticRecalc()) );
                if ( !sl[5].isEmpty() )
                {   // Defaultwert setzen
                    e->gui.tog->setChecked( sl[5].toInt() != 0 );
                }
            }
            if ( !sl[5].isEmpty() )
            {   // set default in local variable
                e->value.flag = sl[5].toInt() != 0;
            }
        }
        else if ( sl[3].startsWith("lbl") && gui )
        {   // Infolabel : "lbl" (hier wird der Prompt-Text über beide Spalten gezogen)
            // und es wird nicht in der Parameter-Struktur gespeichert.
            subLayout->addWidget( e->lbl1, sl[1].toInt(), sl[0].toInt()*3+0, 1, 2 );
            continue;
        }
        else
        {
            qDebug() << "calcHelper: Unknown" << sl;
            continue;
        }
        if ( !sl[4].trimmed().isEmpty() )
        {
            //e->lbl1->setToolTip( sl[4].trimmed() ); -> nur Info bzgl. Kontext-Menü
            e->gui.w->setToolTip( sl[4].trimmed() );
            //qDebug() << "TOOLTIP:" << c->methodName() << e->lbl1->text() << sl[4];
        }
        params.insert( e->key, e );
        //qDebug() << "calcHelper: Use" << sl;
        subLayout->addWidget( e->lbl1,  sl[1].toInt(), sl[0].toInt()*3+0 );
        subLayout->addWidget( e->gui.w, sl[1].toInt(), sl[0].toInt()*3+1 );
    }
    int cmax = subLayout->columnCount();
    // Die letzte Zeile/Spalte als leeren Stretchbereich einfügen
    subLayout->setColumnStretch( subLayout->columnCount(), 1 );
    subLayout->setRowStretch( subLayout->rowCount(), 1 );
    // Zwischen den Spalten eine vertikale Linie
    for ( int c=2; c<cmax; c+=3 )
    {
        QFrame *frm = new QFrame;
        frm->setFrameShape( QFrame::VLine );
        subLayout->addWidget( frm, 0, c, subLayout->rowCount(), 1 );
    }
    if ( cfg ) globsets.endGroup();
}
