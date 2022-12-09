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
