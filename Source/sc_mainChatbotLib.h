#include "sc_maingui.h"
#include "ui_sc_maingui.h"


#ifndef ChatbotDisabled

#include <QFileDialog>
#include <QMessageBox>
#include <QJsonDocument>
#include <QJsonArray>


/**
 * Sammlung von Unterprogrammen aus der MainGui für das Chatbot-Interface.
 * Hier wird keine neue Klasse definiert, sondern die Routinen nur aus dem
 * Sourcefile sc_mainGui.cpp ausgelagert, der Übersichtlichkeit wegen.
 */



/**
 * @brief SC_MainGUI::on_butChatbotSearch_clicked [Button Slot]
 * Auswahl des Chatbot-Files, welches vom Chatbot geschrieben werden kann und vom
 * Hintergrundprozess gelesen und interpretiert wird.
 */
void SC_MainGUI::on_butChatbotSearch_clicked()
{
    QString fn = QFileDialog::getSaveFileName( this, "Open ChatBot config file", ui->inpChatbotFile->text(),
                                              "config file (*.json)", nullptr, QFileDialog::DontConfirmOverwrite );
    if ( fn.isEmpty() ) return;
    ui->inpChatbotFile->setText(fn);
}


/**
 * @brief SC_MainGUI::on_butChatbotLogSearch_clicked [Button Slot]
 * Auswahl des Logfiles für die Berechnungen des Hintergrundprozesses.
 */
void SC_MainGUI::on_butChatbotLogSearch_clicked()
{
    QString fn = QFileDialog::getSaveFileName( this, "Open ChatBot logfile", ui->inpChatbotLogfile->text(),
                                              "logfile (*.log)" );
    if ( fn.isEmpty() ) return;
    ui->inpChatbotLogfile->setText(fn);
    ui->lblChatbotOpenLogfile->setText("<a href=\"file:///"+ui->inpChatbotLogfile->text()+"\">Open</a>");
    ui->lblChatbotOpenLogfile->setToolTip("Open "+ui->inpChatbotLogfile->text());
}


/**
 * @brief SC_MainGUI::on_butChatbotStart_clicked [Button Slot]
 * Start des Hintergrundprozesses. Alle Werte im Calculation-Tab werden hier als
 * Defaultwerte angenommen falls im Chatbotfile kein Wert zu dem Parameter steht.
 */
void SC_MainGUI::on_butChatbotStart_clicked()
{
    ui->butChatbotStart->setEnabled(false);
    ui->butChatbotStop->setEnabled(true);

    QString fnpar = QDir::tempPath()+"/tmpForChatbot.ini";
    performSaveParamOperation( fnpar );

    if ( chatbotBackProg == nullptr )
    {
        chatbotBackProg = new QProcess;
        connect( chatbotBackProg, SIGNAL(errorOccurred(QProcess::ProcessError)),
                this, SLOT(chatbotBackProg_error(QProcess::ProcessError)) );
        connect( chatbotBackProg, SIGNAL(finished(int,QProcess::ExitStatus)),
                this, SLOT(chatbotBackProg_finished(int,QProcess::ExitStatus)) );
        connect( chatbotBackProg, SIGNAL(readyReadStandardError()),
                this, SLOT(chatbotBackProg_readyRead()) );
        connect( chatbotBackProg, SIGNAL(readyReadStandardOutput()),
                this, SLOT(chatbotBackProg_readyRead()) );
    }
    ui->lisChatbotLogOut->clear();

    QStringList params;
    params << "--chatbot" << ui->inpChatbotFile->text();
    if ( ui->togChatbotDebug->isChecked() ) params << "--cbkeep";
    params << "-p" << fnpar;
    params << "-t" << QString::number(ui->inpNumCores->value());
    if ( ui->cbsChatbotColor->currentIndex() != 3 )     // 0=grey, 1=glow, 2=earth, 3,def=temp
        params << "--color" << ui->cbsChatbotColor->currentText();
    if ( ui->cbsChatbotRotation->currentIndex() > 0 )
        params << "--rotate" << QString::number(ui->cbsChatbotRotation->currentIndex());
    if ( ui->cbsChatbotZoom->currentIndex() > 0 )
        params << "--zoom" << ui->cbsChatbotZoom->currentText();
    if ( ui->togChatbotHor->isChecked() || ui->togChatbotVert->isChecked() )
        params << "--swap" << (QString(ui->togChatbotHor->isChecked()?"H":"") + QString(ui->togChatbotVert->isChecked()?"V":""));
    if ( ui->grpChatbotLogfile->isChecked() && ! ui->inpChatbotLogfile->text().isEmpty() )
        params << "--logfile" << ui->inpChatbotLogfile->text();

    chatbotBackProgAddLog( "Params: "+params.join(" ") );
    chatbotBackProgAddLog( "Start: "+chatbotConsProg );
    chatbotBackProgAddLog( "----------" );

    chatbotBackProg->start( chatbotConsProg, params );
}


/**
 * @brief SC_MainGUI::on_butChatbotStop_clicked [Button Slot]
 * Der Hintergrundprozess wird gestoppt.
 */
void SC_MainGUI::on_butChatbotStop_clicked()
{
    if ( chatbotBackProg != nullptr && chatbotBackProg->state() != QProcess::NotRunning )
        chatbotBackProg->kill();
    else
    {
        ui->butChatbotStart->setEnabled(true);
        ui->butChatbotStop->setEnabled(false);
    }
}


/**
 * @brief SC_MainGUI::chatbotBackProg_error [QProcess Slot]
 * @param err - Fehlercode vom System
 * Fehlermeldung vom Hintergrundprozess
 */
void SC_MainGUI::chatbotBackProg_error( QProcess::ProcessError err )
{
    switch ( err )
    {
    case QProcess::FailedToStart:   // The process failed to start. Either the invoked program is missing, or you may have insufficient permissions to invoke the program.
        QMessageBox::critical( this, "Chatbot input",
                              "Background process executable not found.\n"+chatbotBackProg->program(),
                              QMessageBox::Ok );
        break;
    case QProcess::Crashed:         // The process crashed some time after starting successfully.
    case QProcess::Timedout:        // The last waitFor...() function timed out. The state of QProcess is unchanged, and you can try calling waitFor...() again.
        break;
    case QProcess::WriteError:      // An error occurred when attempting to write to the process. For example, the process may not be running, or it may have closed its input channel.
    case QProcess::ReadError:       // An error occurred when attempting to read from the process. For example, the process may not be running.
        QMessageBox::critical( this, "Chatbot input",
                              "Background process unable to write to output device.",
                              QMessageBox::Ok );
        break;
    case QProcess::UnknownError:    // An unknown error occurred. This is the default return value of error().
        break;
    }
}


/**
 * @brief SC_MainGUI::chatbotBackProg_finished [QProcess Slot]
 * @param sta - Status vom System
 * Ende-Meldung vom Hintergrundprozess
 */
void SC_MainGUI::chatbotBackProg_finished( int /*code*/, QProcess::ExitStatus sta )
{
    if ( sta == QProcess::CrashExit )
        chatbotBackProgAddLog("Background process exited.");
    ui->butChatbotStart->setEnabled(true);
    ui->butChatbotStop->setEnabled(false);
}


/**
 * @brief SC_MainGUI::chatbotBackProg_readyRead [QProcess Slot]
 * Meldung vom Hintergrundprozess, dass eine Ausgabe verfügbar ist, die dann
 * in die Liste in der GUI geschrieben wird.
 */
void SC_MainGUI::chatbotBackProg_readyRead()
{
    QString str;
    str = chatbotBackProg->readAllStandardOutput();
    chatbotBackProgAddLog( str.trimmed() );
    str = chatbotBackProg->readAllStandardError();
    chatbotBackProgAddLog( str.trimmed() );
}


/**
 * @brief SC_MainGUI::chatbotBackProgAddLog
 * @param msg - Daten von stdout / stderr
 * Hilfsfunktion zum Auswerten der Meldungsstrings vom Hintergrundprozess und einfügen
 * in die Liste in der GUI. Dort werden maximal 500 Zeilen gespeichert.
 * Hier wird auch das Image gelesen und angezeigt, sofern das gewünscht wird.
 */
void SC_MainGUI::chatbotBackProgAddLog( QString msg )
{
    if ( msg.isEmpty() ) return;
    QStringList sl = msg.split(EOL);
    ui->lisChatbotLogOut->addItems(sl);
    ui->lisChatbotLogOut->scrollToBottom();
    while ( ui->lisChatbotLogOut->count() > 500 ) ui->lisChatbotLogOut->takeItem(0);
#ifndef ChatbotIgnoreImages
    if ( ui->togChatbotClpShowImg->isChecked() )
    {   // Wenn das Bild erzeugt wurde, soll es angezeigt werden.
        // Die Meldung lautet dann: "IMG: <filename> -> <zeit>ms"
        int pos = msg.indexOf("IMG:");
        if ( pos < 0 ) return;
        QString fn = msg.mid(pos+4).trimmed();
        pos = fn.indexOf(".png");
        fn.truncate(pos+4);
        widImage *img = new widImage("Chatbot","");
        connect( img, SIGNAL(destroyed(QObject*)), this, SLOT(imageWindowClosed(QObject*)) );
        images.append(img);
        img->loadImageFile(fn);
        img->show();
        QList<QListWidgetItem*> items = ui->lisDataWindows->findItems( img->windowTitle(), Qt::MatchStartsWith );
        if ( items.size() == 0 )
        {
            ui->lisDataWindows->addItem( img->windowTitle()+" ("+img->getMetaTitle()+")" );
        }
    }
#endif
}


/**
 * @brief SC_MainGUI::on_butChatbotTrainfileSearch_clicked [Button SLot]
 * Auswahl des Trainingsfiles für den Chatbot.
 */
void SC_MainGUI::on_butChatbotTrainfileSearch_clicked()
{
    QString fn = QFileDialog::getSaveFileName( this, "Save ChatBot trainfile", ui->inpChatbotTrainfile->text(),
                                              "Train Config (*.json)", nullptr, QFileDialog::DontConfirmOverwrite );
    if ( fn.isEmpty() ) return;
    ui->inpChatbotTrainfile->setText(fn);
    ui->lblChatbotOpenConfig->setText("<a href=\"file:///"+ui->inpChatbotTrainfile->text()+"\">Open</a>");
    ui->lblChatbotOpenConfig->setToolTip("Open "+ui->inpChatbotTrainfile->text());
}


void SC_MainGUI::generateKeyHash()
{
    // ACHTUNG: auch in sc_mainCons.cpp:waitForChatbot()
    paramkey2jsonkey.clear();
    QStringList allkeys = calcGui->paramsForMethod(false/*num*/, true/*glob*/, false/*fit*/ );
    while ( allkeys.size() > 0 )
    {
        QString key = allkeys.first();
        if (key.startsWith("EditAxis")     || key.startsWith("Editdom")        ||
            key.startsWith("ExpandImage")  || key.startsWith("RadioButtonQ")   ||
            key.contains("EditQsteps")     || key.contains("EditQmin")         ||
            key.startsWith("P1")           || key.contains("EditBFactor")      ||
            key.contains("CheckBoxWAXS")   || key.contains("CalcQmax")         ||
            key.startsWith("EditQmaxData") || key.startsWith("EditQmaxPreset") ||
            key.contains("CenterBeam")     || key.contains("CenterMidpoint")   ||
            key.contains("EditRelDis")     || key.contains("EditDist")
            )
        {
            //if ( key.contains("Qmax") ) qDebug() << "REM" << key;
            allkeys.removeAll(key);
            continue;
        }
        if ( calcGui->isCurrentParameterValid(key, false/*forfit*/) )
        {
            //if ( key.contains("Qmax") ) qDebug() << "use" << key;
            QString prtkey;
            if ( key.startsWith("VAx")  ) prtkey = key.mid(1); // Das 'V' stört bei der Ausgabe
            else if ( key.startsWith("Edit") ) prtkey = key.mid(4); // Das 'Edit' stört auch
            else if ( key.startsWith("CheckBox") ) prtkey = key.mid(8);
            else if ( key.startsWith("ComboBox") ) prtkey = "CB"+key.mid(8);
            else if ( key.startsWith("RadBut") ) prtkey = key.mid(6);
            else if ( key.startsWith("RadioButton") ) prtkey = "RB"+key.mid(11);
            else prtkey = key;
            paramkey2jsonkey.insert( key, prtkey.toLower() );
        }
        else
            qDebug() << "generateKeyHash() param not valid:" << key;
        allkeys.removeAll(key); // damit auch doppelte rausgehen
    }
}

/**
 * @brief SC_MainGUI::chatbotSaveConfigHelper
 * @param fTrain
 * @param parkey
 * @param showena
 * @param isena
 * Hilfsfunktion zum Schreiben des Trainingsfiles (für einen Parameter)
 */
void SC_MainGUI::chatbotSaveConfigHelper(QJsonObject &jsParam, QString parkey, QString jskey )
{
    // Versteckte Parameter werden nicht ausgegeben
    if ( !calcGui->isCurrentParameterVisible(parkey) ) return;

    // Unbekannte Parameter werden ignoriert
    double min, max;
    bool countable;
    if ( !calcGui->limitsOfParamValue( parkey, min, max, countable ) ) return;

    // Daten für das Json-Objekt
    QString jsDesc   = "";      // "description": "Unit cell parameter a.",
    QString jsDefVal = "";      // "default_value": None,
    QString jsRange  = "";      // "range": [0.1, 1000.0],
    QString jsUnit   = "";      // "unit": "Å"

    paramHelper *par = calcGui->params.value(parkey,nullptr);
    if ( countable )
    {
        if ( max == 1 )
            jsRange = "[false, true]";
        else
        {
            QStringList sl;
            static QRegularExpression rex("[,({]");
            for ( int i=0; i<par->gui.cbs->count(); i++ )
            {
                QString item = par->gui.cbs->itemText(i);
                int pos = item.indexOf( rex );
                if ( pos > 0 ) item.truncate(pos);
                sl << "\"" + item.trimmed() + "\"";
            }
            jsRange = "[" + sl.join(", ") + "]";
        }
    }
    else if ( min < max && min > -10000 )
        jsRange = QString("[%1, %2]").arg(min).arg(max);

    if ( par->gui.w->objectName().startsWith("inp") )
        jsUnit = par->gui.numd->suffix().trimmed();

    if ( ! par->tooltip.isEmpty() && jskey!=par->tooltip )
        jsDesc = par->tooltip;

    // Hier eine Liste, welche Parameter NICHT mit einem Default versehen werden sollten,
    // damit der User nicht zuviele Werte eingeben muss.
    QStringList slNoDef = {"EditRadius", "EditRadiusi", "EditSigma", "EditDbeta",
                           "Length", "SigmaL", "HKLmax", "EditQmax"};
    QString jsDefValNum = "";
    if ( ! slNoDef.contains(parkey) )
    {
        jsDefVal = calcGui->currentParamValueStr(parkey, true/*text*/ );
        if ( jsDefVal.contains("{") )
        {
            // Bei den ComboBoxen wird der Text in das JSon File übernommen, mit den Zahlen kann der Chatbot-Anwender weniger anfangen.
            // Die Zahl wird als Zusatz-Wert eingefügt. Ob der Chatbot etwas damit macht, bleibt abzuwarten.
            int p1 = jsDefVal.indexOf("{");
            int p2 = jsDefVal.indexOf("}");
            jsDefValNum = jsDefVal.mid(p1+1,p2-p1-1).trimmed();
            jsDefVal.truncate(p1);
            jsDefVal = jsDefVal.trimmed();
        }
    }
    else
        jsDefValNum = calcGui->currentParamValueStr(parkey, true/*text*/ ); // Bei def val = null wird der aktuelle Wert mit angegeben

    QJsonObject jsObj;
    jsObj.insert( "description", jsDesc.isEmpty() ? QJsonValue::Null : QJsonValue(jsDesc) );
    jsObj.insert( "default_value", jsDefVal.isEmpty() ? QJsonValue::Null : QJsonValue(jsDefVal) );
    if ( !jsDefValNum.isEmpty() )
        jsObj.insert( "default_value_number", jsDefValNum );
    jsObj.insert( "range", jsRange.isEmpty() ? QJsonValue::Null : QJsonValue(jsRange) );
    jsObj.insert( "unit", jsUnit.isEmpty() ? QJsonValue::Null : QJsonValue(jsUnit) );

    //qDebug() << "JSW" << jskey << "-->" << (jsDefVal.isEmpty() ? "NULL" : jsDefVal);

    jsParam.insert( jskey, jsObj );

}


/**
 * @brief SC_MainGUI::on_butChatbotSaveConfig_clicked [Button Slot]
 * Schreiben (append) des Trainingsfiles für den Chatbot.
 */
void SC_MainGUI::on_butChatbotSaveConfig_clicked()
{
    QString descr = ui->inpChatbotExpDescr->text();
    if ( descr.isEmpty() ) descr = ui->txtComment->text();
    chatbotSaveConfig(ui->inpChatbotTrainfile->text().trimmed(), descr);
}


/**
 * @brief SC_MainGUI::chatbotSaveConfig
 * @param fn    - filename
 * @param descr - Beschreibung des Experiments
 * Schreiben (append) des Trainingsfiles für den Chatbot. Vom Button und vom Autoprocessing genutzt.
 */
void SC_MainGUI::chatbotSaveConfig(QString fn, QString descr)
{
    QString str;

    QJsonObject jsObjExp;
    jsObjExp.insert( "description", descr.isEmpty() ? QJsonValue::Null : QJsonValue(descr) );
    jsObjExp.insert( "default_value", QJsonValue::Null );
    jsObjExp.insert( "range", "string" );
    jsObjExp.insert( "unit", "text" );

    QJsonObject jsObjParams;
    jsObjParams.insert( "Experiment_name", jsObjExp );
    // Das QJsonObject.toJson() gibt alles alphabetisch sortiert aus. Daher dieses hier
    //  mit einem großen 'E' und alle andere Keys werden in Kleinbuchstaben gespeichert.

    generateKeyHash();
    QHash<QString,QString>::const_iterator ik = paramkey2jsonkey.constBegin();
    while ( ik != paramkey2jsonkey.constEnd() )
    {
        chatbotSaveConfigHelper( jsObjParams, ik.key(), ik.value() );
        ++ik;
    }

    QJsonObject jsObjPart;
    str = ui->cbsComboBoxParticle->currentText();
    int pos = str.indexOf("{");
    if ( pos > 0 ) str.truncate(pos);
    jsObjPart.insert( "name", str.trimmed() );
    jsObjPart.insert( "description", ui->txtComment->text() );
    jsObjPart.insert( "parameters", jsObjParams );

    QJsonObject jsObjPartType;
    jsObjPartType.insert( "particle_type", jsObjPart );

    QJsonDocument jsTrain;
    jsTrain.setObject(jsObjPartType);

    QFile fTrain(fn);
    bool isappended = true;
    //if ( ! fTrain.open(QIODevice::Append) ) // ist zum Testen einfacher, wenn die Datei immer überschrieben wird
    {
        isappended = false;
        if ( ! fTrain.open(QIODevice::WriteOnly) )
        {
            QMessageBox::critical(this,"Saving training configuration",fTrain.errorString(),QMessageBox::Ok);
            return;
        }
    }

    fTrain.write(jsTrain.toJson());

    ui->statusbar->showMessage("File "+fTrain.fileName()+(isappended?" appended.":" generated."),5000);
}


/**
 * @brief SC_MainGUI::loadChatbotParamFile
 * @param fn - Filename
 * Lese alle Daten aus dem JSON-File in der Art, wie diese oben gespeichert wurden.
 */
void SC_MainGUI::loadChatbotParamFile(QString fn, bool logLoading)
{
    QFile fin(fn);
    if ( ! fin.open(QIODevice::ReadOnly) )
    {
        qDebug() << "ERROR" << fin.fileName() << fin.errorString();
        return;
    }
    QJsonDocument jsInputFile = QJsonDocument::fromJson( fin.readAll() );
    fin.close();

    // ACHTUNG: auch in sc_mainCons.cpp:waitForChatbot()

    QJsonObject jsAllObj = jsInputFile.object();

    QJsonObject jsObjParticle = jsAllObj.take("particle_type").toObject();  // jetzt ist jsAllObj leer!
    QString expDescr = jsval2str( "description", jsObjParticle.value("description") );

    QJsonObject jsObjParams = jsObjParticle.take("parameters").toObject();
    QJsonObject jsObjExpNam = jsObjParams.take("Experiment_name").toObject();
    if ( expDescr == "NULL" )
        expDescr = jsval2str( "description", jsObjExpNam.value("description") );
    if ( logLoading ) calcGui->loadParamLogmessage("JSON",expDescr,"description","");
    ui->inpChatbotExpDescr->setText(expDescr);
    ui->txtComment->setText(expDescr);

    generateKeyHash();
    showColorMarker(SETCOLMARK_CLEARED);

    QJsonObject::const_iterator it;
    for ( it=jsObjParams.constBegin(); it != jsObjParams.constEnd(); ++it )
    {
        QJsonObject jobj = it.value().toObject();
        QJsonValue joval = jobj.value("value");
        if ( joval.isUndefined() )
            joval = jobj.value("default_value");
        QString val = jsval2str( it.key(), joval );
        //qDebug() << it.key() << "==>" << val << paramkey2jsonkey.key(it.key(),"?");
        // val: <zahl oder -> / "NULL" / "False" / "True" / combobox-wert

        QString pkey = paramkey2jsonkey.key(it.key(),"?");

        if ( val.at(0).isNumber() || val.at(0) == '-' || val.at(0) == '+' )
        {
            //qDebug() << it.key() << "==>" << val << pkey << "NUMBER";
            calcGui->updateParamValue( pkey, val.toDouble(), SETCOLMARK_IGNORED );
            if ( logLoading ) calcGui->loadParamLogmessage("JSON",val,pkey,"");
        }
        else if ( val.at(0) == 'F' || val.at(0) == 'T' )
        {
            //qDebug() << it.key() << "==>" << val << pkey << "BOOL";
            calcGui->updateParamValue( pkey, val.at(0)=='T', SETCOLMARK_IGNORED );
            if ( logLoading ) calcGui->loadParamLogmessage("JSON",val,pkey,"");
        }
        else if ( val == "NULL" )
        {
            /*
            Debug: "dbeta" ==> "NULL" "EditDbeta" -NULL-
            Debug: "gridpoints" ==> "NULL" "GridPoints" -NULL-
            Debug: "radius" ==> "NULL" "EditRadius" -NULL-
            Debug: "sigma" ==> "NULL" "EditSigma" -NULL-
            */
            val = jsval2str( it.key(), jobj.value("default_value_number") );    // Mehr für mich ohne reale Daten vom Chatbot
            calcGui->updateParamValue( pkey, val.toDouble(), SETCOLMARK_IGNORED );
            //qDebug() << it.key() << "==>" << val << pkey << "-NULL-";
            if ( logLoading ) calcGui->loadParamLogmessage("JSON",val,it.key(),"");
        }
        else
        {
            //qDebug() << it.key() << "==>" << val << pkey << "CBS";
            myGuiParam *gp = myGuiParam::getGuiParam(pkey);
            if ( gp == nullptr || gp->cbs() == nullptr )
            {
                qDebug() << "*** UNKNOWN ***";
                if ( logLoading ) calcGui->loadParamLogmessage("JSON",val,it.key(),"UNKNOWN KEY");
            }
            else
            {
                int p = gp->cbs()->findText( val, Qt::MatchStartsWith );
                if ( p >= 0 )
                {
                    gp->cbs()->setCurrentIndex(p);
                    //SETCOL( gp->cbs(), SETCOLMARK_IGNORED );
                    if ( logLoading ) calcGui->loadParamLogmessage("JSON",val,it.key(),gp->cbs()->currentText());
                }
                else
                    if ( logLoading ) calcGui->loadParamLogmessage("JSON",val,it.key(),"UNKOWN VALUE");
                // Die Callbacks der Comboboxen werden zum Schluss nach dem Laden vom Hauptprogramm aufgerufen, damit die Eingaben geschaltet werden
            }
        }
    } // for it

    // Set some defaults
    if ( logLoading ) calcGui->loadParamLogmessage("","","","[JSON] Set default quadrants" );
    ui->radQ1->setChecked( false );
    ui->radQ2->setChecked( false );
    ui->radQ4->setChecked( true  );
    ui->togExpandImage->setChecked( false );
    ui->togExpandImage->setEnabled( false );
    ui->radEditQmaxData->setChecked( false );
    ui->radEditQmaxPreset->setChecked( true );
}


QString SC_MainGUI::jsval2str( QString key, QJsonValue val )
{
    switch ( val.type() )
    {
    case QJsonValue::Null:
        return "NULL";
    case QJsonValue::Bool:
        return val.toBool() ? "true" : "false";
    case QJsonValue::Double:
        return QString::number(val.toDouble());
    case QJsonValue::String:
        return val.toString();
    case QJsonValue::Array:
    {
        QJsonArray arr = val.toArray();
        qDebug() << "JSR" << key << "Array" << arr.size();
        QString rv = "[ " + jsval2str(key+"[0]", arr[0]);
        for ( int i=1; i<arr.size(); i++ )
        {
            rv += ", " + jsval2str(key+QString("[%1]").arg(i), arr[i]);
        }
        return rv+" ]";
    }
    case QJsonValue::Object:
    {
        /*
        QJsonObject jobj = val.toObject();
        QJsonObject::const_iterator it;
        for ( it=jobj.constBegin(); it != jobj.constEnd(); ++it )
        {
            QString str = jsval2str( it.key(), it.value() );
            if ( it.key().contains("value",Qt::CaseInsensitive) && ! it.key().contains("number",Qt::CaseInsensitive) )
            {
                qDebug() << "JSR" << key << "-->" << str;
            }
        }
        */
        return "obj";
    }
    case QJsonValue::Undefined:
        return "NULL";
    }
    qDebug() << "jsr" << key << "?" << val;
    return "";
}


/**
 * @brief SC_MainGUI::on_butChatbotReadClipboard_clicked [Button SLot]
 * Lesen des Clipboards, schreibe das File und der Hintergrundprozess rechnet.
 */
void SC_MainGUI::on_butChatbotReadClipboard_clicked()
{
    QFile f(ui->inpChatbotFile->text());
    if ( ! f.open(QIODevice::WriteOnly) )
    {
        ui->statusbar->showMessage( f.fileName()+" "+f.errorString(), 5000 );
        ui->butChatbotReadClipboard->setEnabled(false);
        return;
    }
    f.write( qPrintable(qApp->clipboard()->text()) );
    f.write(EOL);
    f.close();
    qApp->clipboard()->clear();
    if ( chatbotBackProg == nullptr || ! chatbotBackProg->isOpen() )
    {   // Jetzt läuft der Hintergrundprozess NICHT ...
        on_butChatbotStart_clicked(); // Macht sofort die erste Berechnung
    }
    // Ansonsten wird der Prozess die Datei finden.
}


/**
 * @brief SC_MainGUI::chatbotClipboardChanged [Slot]
 * Meldung, wenn etwas im Clipboard geändert wurde. Hier wird dann der Button freigegeben.
 */
void SC_MainGUI::chatbotClipboardChanged(QClipboard::Mode)
{
    if ( chatbotConsProg.isEmpty() ) return;
    qDebug() << "chatbotClipboardChanged - Check (Start)";
    ui->butChatbotReadClipboard->setEnabled( ! qApp->clipboard()->text().isEmpty() );
    qDebug() << "chatbotClipboardChanged - Check" << ui->butChatbotReadClipboard->isEnabled();
}

#endif // ChatbotDisabled
