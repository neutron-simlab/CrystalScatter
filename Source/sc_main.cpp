#include "sc_maingui.h"
#include <QApplication>

#ifdef WIN32
#include "debughandler.h"
#else
#include <QLoggingCategory>
#endif


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

#ifdef WIN32
    // To avoid empty lines printed every second after one qDebug() output
    qInstallMessageHandler(debugWinMsgHandler);
#else
    // To enable qDebug() under Linux (https://forum.qt.io/topic/81430/unable-to-see-the-qdebug-messages-on-console/15)
    QLoggingCategory::defaultCategory()->setEnabled(QtDebugMsg, true);
#endif

    SC_MainGUI w;
    w.show();
    return a.exec();
}
