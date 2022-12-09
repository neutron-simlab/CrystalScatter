#ifndef SC_GLOBALCONFIG_H
#define SC_GLOBALCONFIG_H

/* Hier stehen alle Definitionen (kein Code), welche an verschiedenen Stellen des Programms
 * bedingte Kompilationen auslösen.
 */


//#define COPY_FITDATA_TO_GPU
// Wenn definiert, wird das anzufittende Bild beim Start in einen Speicher in die GPU kopiert
// und dort werden die Vergleiche bzgl. Corner-Pixel-Mask durchgeführt (wird besonders auf dem
// PC wohl schneller).
// Die FQS wird dort nicht berechnet, da es recht aufwendig ist, unter CUDA das Thread-Locking zu machen.
// Vielleicht später, wenn ich den Lehrgang durchgemacht habe.

//#define CALC_FQS_IN_GPU
// Wenn COPY_FITDATA_TO_GPU definiert ist, wird hiermit gleichzeitig die Fehlerquadratsumme in
// einem weiteren Array pixelweise festgehalten und in der Simplex-Routine kann dann die Summe
// ohne weitere Berechnungen schnell bestimmt werden.


//#ifndef CONSOLENPROG
#define FITDATA_IN_GPU  // Neu im Juli 2022
// Wenn definiert, wird das anzufittende Bild beim Start in einen Speicher in die GPU kopiert
// und dort wird die FQS berechnet. Dabei kann das anzufittende Bild eine andere Größe haben,
// der Algorithmus geht pixelweise durch und mittels der Werte in _latticeForFit (in sc_libs.h)
// werden die passenden qx,qy,qz Werte berechnet und dann der simulierte Intensitätswert, der
// dann sofort mit dem Datenpixel verarbeitet werden kann. Die jeweiligen FQS-Werte werden
// in einem weiteren Array gespeichert und später von der CPU addiert. Das ist einfacher, als
// bei jedem Zugriff das Locking zu machen. (Vielleicht wird es später mal anders probiert)
//#endif


// Vergleich verschiedener Berechnungen für die FQS:
#define FQSVERGL(A,B) sqr( log(A)/log(10.) - log(B)/log(10.) )  // log(x) is e-based logarithm
//--#define FQSVERGL(A,B) sqr( log(A) - log(B) )
//--#define FQSVERGL(A,B) sqr( log10(A) - log10(B) )



// Definitionen für die QSettings
#define SETT_APP   "JCNS-1-SasCrystal"  // auch in widimage.cpp und dlgtimetests.cpp
#define SETT_GUI   "GUISettings"
#define SETT_PAR   "Parameter"
//#define SETT_AIGEN "AI"



// Definition der Breite der Eingabeelemente (cbs,inp) im Calculation-Tab, wird in sc_calcgui.cpp verwendet
// Fehlt der Eintrag, wird der Qt-Default genommen
#ifdef Q_OS_WIN
#define CALC_INPUT_MAXWIDTH  85
#else
#define CALC_INPUT_MAXWIDTH  180
#endif


#endif // SC_GLOBALCONFIG_H
