#ifndef SC_GLOBALCONFIG_H
#define SC_GLOBALCONFIG_H

#ifndef __CUDACC__
#include <QtGlobal>
#endif


/* Hier stehen alle Definitionen (kein Code), welche an verschiedenen Stellen des Programms
 * bedingte Kompilationen auslösen.
 */


// This version information is displayed in the top right corner in the GUI and with -v in the console
#define MYVERSION "2.1.1  (Jun 2024)"


#define FITDATA_IN_GPU  // Neu im Juli 2022
// Wenn definiert, wird das anzufittende Bild beim Start in einen Speicher in der GPU kopiert
// und dort wird die FQS berechnet. Dabei kann das anzufittende Bild eine andere Größe haben,
// der Algorithmus geht pixelweise durch und mittels der Werte in _latticeForFit (in sc_libs.h)
// werden die passenden qx,qy,qz Werte berechnet und dann der simulierte Intensitätswert, der
// dann sofort mit dem Datenpixel verarbeitet werden kann. Die jeweiligen FQS-Werte werden
// in einem weiteren Array gespeichert und später von der CPU addiert. Das ist einfacher, als
// bei jedem Zugriff das Locking zu machen. (Vielleicht wird es später mal anders probiert)


// Vergleich verschiedener Berechnungen für die FQS:
#define FQSVERGL(A,B) sqr( log(A)/log(10.) - log(B)/log(10.) )  // log(x) is e-based logarithm
//--#define FQSVERGL(A,B) sqr( log(A) - log(B) )
//--#define FQSVERGL(A,B) sqr( log10(A) - log10(B) )



// Definitionen für die QSettings
#define SETT_APP   "JCNS-1-SasCrystal"  // auch in widimage.cpp und dlgtimetests.cpp
#define SETT_GUI   "GUISettings"
#define SETT_PAR   "Parameter"
//#define SETT_AIGEN "AI"



// Definition der Breite der Eingabeelemente (cbs,inp) im Calculation-Tab, wird in myguiparam.cpp verwendet
// Fehlt der Eintrag, wird der Qt-Default genommen
#ifdef Q_OS_WIN
#define CALC_INPUT_MAXWIDTH  90
#else
#define CALC_INPUT_MAXWIDTH  100
#endif



// Da es beim Laden von Parameter-Files u.U. zu unklaren Verzögerungen kommt, können mit dieser
// Definition zu Beginn jeder Routine die Probleme eventuell erkannt werden...
// Vor dem Upload nach Github sollte das allerdings wieder raus.
#define DT(x) //x


#ifndef __CUDACC__

// Syntaxänderungen bei neueren Qt-Versionen
#if (QT_VERSION <= QT_VERSION_CHECK(5, 12, 11))
#define SPLIT_SKIP_EMPTY_PARTS QString::SkipEmptyParts
#define SPLIT_KEEP_EMPTY_PARTS QString::KeepEmptyParts
#define FONTMETRIC_WIDTH(s) width(s)
#else
#define SPLIT_SKIP_EMPTY_PARTS Qt::SkipEmptyParts
#define SPLIT_KEEP_EMPTY_PARTS Qt::KeepEmptyParts
#define FONTMETRIC_WIDTH(s) horizontalAdvance(s)
#endif

#endif // CUDACC


// Da in der Routine ButtonHKLClick() u.U. bestimmte Variablen modifiziert werden, sollten diese auch
// im Umfeld mit aktualisiert werden. Dazu gibt es diesen Returnwert.
typedef enum
{
    chgNone,    // Nichts geändert
    // params.* wurden geändert
    chgB,       // ucb
    chgBC,      // ucb, ucc
    chgBCabc    // ucb, ucc, alpha_deg, beta_deg, gamma_deg
} rvButtonHKLClick;

#endif // SC_GLOBALCONFIG_H
