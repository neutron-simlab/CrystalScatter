#include "sc_maingui.h"
#include <QApplication>

#ifdef WIN32
#include "debughandler.h"
#else
#include <QLoggingCategory>
#endif

#include <QDebug>
#include <QStyleFactory>


int main(int argc, char *argv[])
{
    //QStringList styl = QStyleFactory::keys();
    //qDebug() << styl;
    // Qt6.7.2@Win11: QList("windows11", "windowsvista", "Windows", "Fusion")
    //                       ^ Crash      ^ Screenshot    ^ Screenshot
    //                                                               ^ noch am ehesten wie Win10 aber dafür sehr hoch...

    // Qt5.15.2@Linux: QList("Adwaita-Dark", "Adwaita", "Windows", "Fusion")
    //                        ^ Dunkel,       ^ Default unter Linux
    //                                                    ^ Used    ^ GroupBox mit Hintergrund u.a.
#ifdef WIN32
#if (QT_VERSION >= QT_VERSION_CHECK(6,0,0))
    //QStyle *style = QStyleFactory::create("windows");         // ==> Macht Probleme beim Debugger (TODO)
    //if ( style != nullptr ) QApplication::setStyle(style);
#endif
#endif
// #ifdef Q_OS_LINUX
//    QApplication::setStyle(QStyleFactory::create("Windows"));
// #endif

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
