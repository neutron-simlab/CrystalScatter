#ifndef SC_GLOBALCONFIG_H
#define SC_GLOBALCONFIG_H

#ifndef __CUDACC__
#include <QtGlobal>
#endif


/* Hier stehen alle Definitionen (kein Code), welche an verschiedenen Stellen des Programms
 * bedingte Kompilationen auslösen.
 */


// This version information is displayed in the top right corner in the GUI and with -v in the console
#define MYVERSION "2.1.3  (Jan 2025)"


#define FITDATA_IN_GPU  // Neu im Juli 2022 - ACHTUNG unten
// Wenn definiert, wird das anzufittende Bild beim Start in einen Speicher in der GPU kopiert
// und dort wird die FQS berechnet. Dabei kann das anzufittende Bild eine andere Größe haben,
// der Algorithmus geht pixelweise durch und mittels der Werte in _latticeForFit (in sc_libs.h)
// werden die passenden qx,qy,qz Werte berechnet und dann der simulierte Intensitätswert, der
// dann sofort mit dem Datenpixel verarbeitet werden kann. Die jeweiligen FQS-Werte werden
// in einem weiteren Array gespeichert und später von der CPU addiert. Das ist einfacher, als
// bei jedem Zugriff das Locking zu machen. (Vielleicht wird es später mal anders probiert)
// ACHTUNG: wenn dieses Symbol deaktiviert wird, ist das Programm nicht mehr compilierbar,
//   da einige Programmstücke im Laufe der Zeit entwickelt wurden und diese Definition dabei
//   nicht mehr bedacht wurde.


// Vergleich verschiedener Berechnungen für die FQS:
#define FQSVERGL(A,B) sqr( log(A)/log(10.) - log(B)/log(10.) )  // log(x) is e-based logarithm
//--#define FQSVERGL(A,B) sqr( log(A) - log(B) )
//--#define FQSVERGL(A,B) sqr( log10(A) - log10(B) )



//#define NOCBSCALLBACK
// If defined, no ComboBox Callbacks are performed and all labels show the parameter names
// and not some special informations. Helpful for screenshots for the documentation.



//#define ChatbotDisabled
// If defined, the complete Chatbot function will not be visible to the user and no functions
// will be active in the background.

//#define ChatbotIgnoreImages  // Not yet implemented
// If defined, the generated images from the Chatbot function will not be shown automatically.
// This function only generates PNG files and not normal datafiles with meta informations.
// Only valid if ChatbotDisabled is not defined.



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
    // Qt 5.12.11 oder früher
#define SPLIT_SKIP_EMPTY_PARTS QString::SkipEmptyParts
#define SPLIT_KEEP_EMPTY_PARTS QString::KeepEmptyParts
#define FONTMETRIC_WIDTH(s) width(s)
#else
    // Qt 5.12.12 oder neuer
#define SPLIT_SKIP_EMPTY_PARTS Qt::SkipEmptyParts
#define SPLIT_KEEP_EMPTY_PARTS Qt::KeepEmptyParts
#define FONTMETRIC_WIDTH(s) horizontalAdvance(s)
#endif

#if (QT_VERSION <= QT_VERSION_CHECK(6,0,0))
    // Qt 5 und früher
#define XMLCOMP(x,s) ((x) == s)
#else
    // Qt 6 und neuer
#define XMLCOMP(x,s) ((x).compare(s) == 0)
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
