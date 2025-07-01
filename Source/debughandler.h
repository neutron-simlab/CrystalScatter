#ifndef DEBUGHANDLER_H
#define DEBUGHANDLER_H

//#define MitSourceInfos

#ifdef USELOGWINDOW

#include "debugwindow.h"

void debugWinMsgHandler( QtMsgType type, const QMessageLogContext &context, const QString &msg )
{
    DebugWindow* deb = DebugWindow::self();
    deb->show();
    switch (type) {
    case QtDebugMsg:
        deb->append( QString("Debug: %1 (%2:%3, %4)").arg(msg).arg(context.file).arg(context.line).arg(context.function));
        break;
    case QtInfoMsg:
        deb->append( QString("Info: %1 (%2:%3, %4)").arg(msg).arg(context.file).arg(context.line).arg(context.function));
        break;
    case QtWarningMsg:
        deb->append( QString("Warning: %1 (%2:%3, %4)").arg(msg).arg(context.file).arg(context.line).arg(context.function));
        break;
    case QtCriticalMsg:
        deb->append( QString("Critical: %1 (%2:%3, %4)").arg(msg).arg(context.file).arg(context.line).arg(context.function));
        break;
    case QtFatalMsg:
        // we need to create a modal dialogue here as everything else would crash with us
        QMessageBox::critical(nullptr, "Debug - Fatal", msg);
        abort();
    }
}

#else

#include <QMessageBox>

void debugWinMsgHandler( QtMsgType type, const QMessageLogContext &context, const QString &msg )
{
#ifndef MitSourceInfos
    Q_UNUSED(context);
#endif
    switch (type) {
    case QtDebugMsg:
#ifdef MitSourceInfos
        fprintf( stderr, "Debug: %s (%s:%d, %s)\n", qPrintable(msg), context.file, context.line, context.function );
#else
        if ( msg.startsWith("onSelectionChange (") ) break;
        fprintf( stderr, "Debug: %s\n", qPrintable(msg) );
#endif
        break;
    case QtInfoMsg:
#ifdef MitSourceInfos
        fprintf( stderr, "Info: %s (%s:%d, %s)\n", qPrintable(msg), context.file, context.line, context.function );
#else
        fprintf( stderr, "Info: %s\n", qPrintable(msg) );
#endif
        break;
    case QtWarningMsg:
        // Warning: setGeometry: Unable to set geometry 481x546+0+39 (frame: 499x593-9+1) on QWidgetWindow/"widImageWindow" on "\\.\DISPLAY1". Resulting geometry: 479x538+1+46 (frame: 497x585-8+8) margins: 9, 38, 9, 9 minimum size: 374x437 MINMAXINFO(maxSize=POINT(x=0, y=0), maxpos=POINT(x=0, y=0), maxtrack=POINT(x=0, y=0), mintrack=POINT(x=486, y=593)))
        if ( msg.contains("setGeometry: Unable to set geometry") ) break;
        // Diese Meldung kommt, weil Qt das Image-Fenster anscheinend nicht richtig anzeigen kann.
        // Meiner Menung nach ist das aber nicht so entscheidend.
#ifdef MitSourceInfos
        fprintf( stderr, "Warning: %s (%s:%d, %s)\n", qPrintable(msg), context.file, context.line, context.function );
#else
        fprintf( stderr, "Warning: %s\n", qPrintable(msg) );
#endif
        break;
    case QtCriticalMsg:
#ifdef MitSourceInfos
        fprintf( stderr, "Critical: %s (%s:%d, %s)\n", qPrintable(msg), context.file, context.line, context.function );
#else
        fprintf( stderr, "Critical: %s\n", qPrintable(msg) );
#endif
        break;
    case QtFatalMsg:
        // we need to create a modal dialogue here as everything else would crash with us
        QMessageBox::critical(nullptr, "Debug - Fatal", msg);
        abort();
    }
#ifdef WIN32
    // under windows the debug output has to be flushed to be readable in the development environment
    fflush(stderr);
#endif
}

#endif

#endif
