#include "sc_postproc.h"
#include <QDebug>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef USE_COMPLEX_WEB
#include "fftw3.h"
#endif
#include "sc_math.h"


#define FROM_FOERSTER


SasCalc_PostProc *SasCalc_PostProc::_instancePtr = nullptr;


SasCalc_PostProc *SasCalc_PostProc::inst()
{
    if ( _instancePtr == nullptr )
        _instancePtr = new SasCalc_PostProc;
    return _instancePtr;
}

SasCalc_PostProc::SasCalc_PostProc()
{
    _arrFFT    = nullptr;
    _arrIFFT   = nullptr;
    _arrPol    = nullptr;
    _arrPolDim = 0; // size = dim*dim
    _arrFFTdim = 0;
    _len       = 0;
    _doLogOutput = true; // per Call abschaltbar...
}


/**
 * @brief SasCalc_PostProc::generateRphi
 * @param x0      - source min x
 * @param x1      - source max x
 * @param y0      - source min y
 * @param y1      - source max y
 * @param Xcenter - source beam center x
 * @param Ycenter - source beam center y
 * @param s       - output size (each dimension)
 * @param d       - source data
 * @param scale   - flag scaling [0,1] lin
 * @param clip    - flag clipping [1e-10,1] log
 * @param clip40  - flag clipping [40%,100%] -> [0,1]
 * @return pointer to output data
 */
double *SasCalc_PostProc::generateRphi(int x0, int x1, int y0, int y1,
                                       int Xcenter, int Ycenter,
                                       int s, const double *d,
                                       bool scale, bool clip, bool clip40 )
{
    // Die 'alten' Routinen sind via GIT abrufbar (vor dem 21.09.2022)
    // Jetzt versuche ich mal den Algorithmus vom Henrich Frielinghaus,
    // den er als F90 Programm am 20.09.2022 geschickt hat.
    // Das 'Problem': der Algorithmus vom Henrich liefert ein 200*150 Array,
    // meine weiteren Verarbeitungen verlangen ein <s>*<s> Array.

    // Das Array wird hier jetzt lokal angelegt (nicht mehr zusammen mit den FFT-Daten)
    if ( _arrPolDim != s && _arrPol != nullptr )
    {
        delete _arrPol;
        _arrPol = nullptr;
    }
    if ( _arrPol == nullptr )
    {
        _arrPolDim = s;
        _arrPol    = new double[s*s];
    }
    // Dieser Algorithmus kommt nicht klar, wenn x0 bzw. y0 < 0 sind.
    if ( x0 < 0 )
    {
        x1      -= x0;
        Xcenter -= x0;
        x0       = 0;
    }
    if ( y0 < 0 )
    {
        y1      -= y0;
        Ycenter -= y0;
        y0       = 0;
    }
    // Parameter in globale Variabln übernehmen
    _dataPtr = d;
    _xmin = x0;
    _xmax = x1;
    _ymin = y0;
    _ymax = y1;
    _len  = (_xmax - _xmin) * (_ymax - _ymin);
    int srcsx = _xmax - _xmin;
    int srcsy = _ymax - _ymin;
    int arrlen = _arrPolDim * _arrPolDim;

    if ( _doLogOutput )
        qDebug() << "POSTPROC:generateRphi x"<<x0<<x1<<srcsx << "y"<<y0<<y1<<srcsy << "bs"<<Xcenter<<Ycenter << "len"<<_len << "dim"<<s;

    //------------------------------------
    // transformation to polar coordinates
    //------------------------------------

    // Temporäre Daten zum Arbeiten
    //double polint[s/*200*/][s/*150*/];    // Intensit"aten in Polarkoordinaten
    //double polier[s/*200*/][s/*150*/];    // Fehler
    //int    polmsk[s/*200*/][s/*150*/];    // Maske

    int       phisteps = s/*200*/;   // range 2..200 gives segments per full circle 2pi
    const int phisymm  = 0;     // 1=on, 0=off  use symmetry
    // Wenn auf 1 gesetzt, dann sollte unten noch der Code angepasst werden

    const int plaincor = 0;     // 1=on, 0=off  correction for planar detector
    // Wenn auf 1 gesetzt, dann müssen <detdist> und <detelem> aus den Daten kommen

    double detdist  = 400.0;    // detector distance
    double detelem  = 0.5;      // detector element
    //double lambda   = 7.0;      // wavelength AA

    //memset( polint, 0, sizeof(polint) );
    //memset( polier, 0, sizeof(polier) );
    //memset( polmsk, 0, sizeof(polmsk) );
    //do i=1,200                 ! reset des Datenfeldes
    // do j=1,150
    //  polint(i,j) = 0.0
    //  polier(i,j) = 0.0
    //  polmsk(i,j) = 0
    // enddo
    //enddo

    // TODO if ( phisymm >= 1 ) phisteps = int( (phisteps+1)*0.5 ) * 2;     // even steps

    double phiber;             // Segment in phi
    phiber   = 2.*M_PI/phisteps;  // angular range of a single segment

    double rmax;               // maximaler Radius in r

#define itsty(x,y)  (data(x,y))
//#define errint(x,y) (data(x,y)*0.01)

    rmax =           sqrt(sqr(       1. -Xcenter) + sqr(       1. -Ycenter));
    rmax = qMax(rmax,sqrt(sqr((srcsx-1.)-Xcenter) + sqr(       1. -Ycenter)));
    rmax = qMax(rmax,sqrt(sqr((srcsx-1.)-Xcenter) + sqr((srcsy-1.)-Ycenter)));
    rmax = qMax(rmax,sqrt(sqr(       1. -Xcenter) + sqr((srcsy-1.)-Ycenter)));

    //qDebug() << "r,phi:" << Xcenter << Ycenter << s << rmax;

    //open (7,file='pol_coord',status='unknown')
    //open (8,file='pol_data', status='unknown')

    int irmax = rmax + 0.5;
    float irstp = float(irmax) / float(s-1);
    //qDebug() << "irmax" << irmax << irstp;
    for ( float ir=0; ir<=irmax; ir+=irstp )
    {
        //double qr = 2.*M_PI/lambda * ir/detdist;
        double rad;                // tats"achlicher Radius
        int    umf;                // Zahl der Winkel in phiber
        double phi;                // tats"achlicher Winkel

        if ( plaincor == 0 )
            rad = ir;
        else
            rad = detdist / detelem * tan( 2.0*asin(0.5*ir*detelem/detdist) );

        for ( int i=0; i<=phisteps/(1+phisymm); i++ )
        {
            double twt = 0.0;      // totales Gewicht    >>>RESET<<<
            double avg = 0.0;      // Wert
            //double avge= 0.0;      // Fehler
            //int    anz = 0;        // Anzahl der Punkte (vielleicht nicht n"otig)

            umf = qMax(1024-int(1024.-phiber*rad),1);

            for ( int j=1; j<=umf; j++ )
            {
                phi = (i-1.5)*phiber + (j-0.5)*phiber/umf;

                //if (ir.eq.2) write(7,*) phi

                double xs  = Xcenter + rad * cos(phi);
                double ys  = Ycenter + rad * sin(phi);

                int ix1 = int(xs);         // Eckpunkte
                int iy1 = int(ys);
                int ix2 = ix1 + 1;
                int iy2 = iy1 + 1;

                double rex = xs - ix1;        // Reste f"ur Gewichtung sp"ater
                double rey = ys - iy1;
                double wt;

                int msk = 0;               // linker unterer Punkt
                //if (ix1>=1 && ix1<=128 && iy1>=1 && iy1<=128)
                if (ix1>=0 && ix1<srcsx && iy1>=0 && iy1<srcsy)
                {
                    msk = 1; // mask(ix1,iy1)  TODO: Bislang habe ich keine Maske für die Daten
                    wt  = msk*(1.0-rex)*(1.0-rey);
                    twt = twt + wt;
                    avg = avg + wt   *itsty(ix1,iy1);
                    //avge= avge+ wt*wt*sqr(errint(ix1,iy1));
                    //anz = anz + msk;
                }

                msk = 0;               // rechter unterer Punkt
                //if (ix2>=1 && ix2<=128 && iy1>=1 && iy1<=128)
                if (ix2>=0 && ix2<srcsx && iy1>=0 && iy1<srcsy)
                {
                    msk = 1; // s.o. mask(ix2,iy1)
                    wt  = msk*(    rex)*(1.0-rey);
                    twt = twt + wt;
                    avg = avg + wt   *itsty(ix2,iy1);
                    //avge= avge+ wt*wt*sqr(errint(ix2,iy1));
                    //anz = anz + msk;
                }

                msk = 0;               // linker oberer Punkt
                //if (ix1>=1 && ix1<=128 && iy2>=1 && iy2<=128)
                if (ix1>=0 && ix1<srcsx && iy2>=0 && iy2<srcsy)
                {
                    msk = 1; // s.o. mask(ix1,iy2)
                    wt  = msk*(1.0-rex)*(    rey);
                    twt = twt + wt;
                    avg = avg + wt   *itsty(ix1,iy2);
                    //avge= avge+ wt*wt*sqr(errint(ix1,iy2));
                    //anz = anz + msk;
                }

                msk = 0;               // rechter oberer Punkt
                //if (ix2>=1 && ix2<=128 && iy2>=1 && iy2<=128)
                if (ix2>=0 && ix2<srcsx && iy2>=0 && iy2<srcsy)
                {
                    msk = 1; // s.o. mask(ix2,iy2)
                    wt  = msk*(    rex)*(    rey);
                    twt = twt + wt;
                    avg = avg + wt   *itsty(ix2,iy2);
                    //avge= avge+ wt*wt*sqr(errint(ix2,iy2));
                    //anz = anz + msk;
                }

            } // for ( int j=1; j<=umf; j++ )

/*
            if (phisymm.ge.1) then

                    do j=1,umf

                        phi = pi + (i-1.5)*phiber + (j-0.5)*phiber/umf

                                                        xs  = Xcenter + rad * cos(phi)
                                        ys  = Ycenter + rad * sin(phi)

                                    ix1 = int(xs)         ! Eckpunkte
                        iy1 = int(ys)
                        ix2 = ix1 + 1
                      iy2 = iy1 + 1

                      rex = xs - ix1        ! Reste f"ur Gewichtung sp"ater
                          rey = ys - iy1

                          msk = 0               ! linker unterer Punkt

                    if (ix1.ge.1 .and. ix1.le.128 .and. iy1.ge.1 .and. iy1.le.128) then
                        msk = mask(ix1,iy1)
                        wt  = msk*(1.0-rex)*(1.0-rey)
                          twt = twt + wt
                          avg = avg + wt   *itsty(ix1,iy1)
                                    avge= avge+ wt*wt*errint(ix1,iy1)**2
                                 anz = anz + msk
                          endif

                              msk = 0               ! rechter unterer Punkt

                    if (ix2.ge.1 .and. ix2.le.128 .and. iy1.ge.1 .and. iy1.le.128) then
                        msk = mask(ix2,iy1)
                        wt  = msk*(    rex)*(1.0-rey)
                          twt = twt + wt
                          avg = avg + wt   *itsty(ix2,iy1)
                                    avge= avge+ wt*wt*errint(ix2,iy1)**2
                                 anz = anz + msk
                          endif

                              msk = 0               ! linker oberer Punkt

                    if (ix1.ge.1 .and. ix1.le.128 .and. iy2.ge.1 .and. iy2.le.128) then
                        msk = mask(ix1,iy2)
                        wt  = msk*(1.0-rex)*(    rey)
                          twt = twt + wt
                          avg = avg + wt   *itsty(ix1,iy2)
                                    avge= avge+ wt*wt*errint(ix1,iy2)**2
                                 anz = anz + msk
                          endif

                              msk = 0               ! rechter oberer Punkt

                    if (ix2.ge.1 .and. ix2.le.128 .and. iy2.ge.1 .and. iy2.le.128) then
                        msk = mask(ix2,iy2)
                        wt  = msk*(    rex)*(    rey)
                          twt = twt + wt
                          avg = avg + wt   *itsty(ix2,iy2)
                                    avge= avge+ wt*wt*errint(ix2,iy2)**2
                                 anz = anz + msk
                          endif

                              enddo

                                  endif
*/

            if (twt>0.0)
            {
                avg = avg / twt;
                //avge= sqrt(avge) / twt;
            }

            int idx = i * s + (ir/irstp);
            if ( idx < arrlen )
                _arrPol[idx] = (avg < 0 ) ? 0 : avg;
            // Im Bild von unten nach oben der Winkel gegen den Uhrzeigersinn
            // (0°=unten ist im ScrBild vom Beamstop nach rechts)

            //polint[i][ir] = avg;
            //polier[i][ir] = avge;
            //if (twt>0.0) polmsk[i][ir] = 1;

        } // for ( int i=1; i<=phisteps/(1+phisymm); i++ )

    } // for ( float ir=1; ir<=irmax && ir<s; ir+=irstp )

    scaleAndClipData( _arrPol, s*s, scale, clip, clip40, false );

    return _arrPol;
}


/**
 * @brief SasCalc_PostProc::calculateIFFT
 * @param foreward - flag foreward(true) or backward(false)
 * @param x0       - source min x
 * @param x1       - source max x
 * @param y0       - source min y
 * @param y1       - source max y
 * @param d        - source data
 * @param s        - output size (each dimension)
 * @param ot       - output type (Re, Im, ...)
 * @param scaled   - flag scaling [0,1] lin
 * @param clip     - flag clipping [1e-10,1] log
 * @param clip40   - flag clipping [40%,100%] -> [0,1]
 * @param swapout  - flag swap output (0° in center)
 * @return pointer to output data
 */
double *SasCalc_PostProc::calculateIFFT(bool foreward,
                                        int x0, int x1, int y0, int y1, // Source size
                                        const double *d,                // Source data
                                        int s,                          // Dest. size
                                        _outType ot,
                                        bool scaled, bool clip, bool clip40,
                                        bool swapout )
{
    // Parameter in globale Variabeln übernehmen
    _dataPtr = d;
    // TODO: Was wenn die Größen (input<->output) unterschiedlich sind ???
    _xmin = x0;
    _xmax = x1;
    _ymin = y0;
    _ymax = y1;
    _len  = (_xmax - _xmin) * (_ymax - _ymin);
    int sx = _xmax - _xmin;
    int sy = _ymax - _ymin;

    if ( _doLogOutput )
        qDebug() << "POSTPROC:calculateIFFT: Input X:" << _xmin << _xmax << sx
                 << "Y:" << _ymin << _ymax << sy
                 << "Outsize:" << s
                 << "foreward:" << foreward;

    // Find Max
    double maxinput=0;
    const double *tmp=_dataPtr;
    for ( int i=0; i<_len; i++, tmp++ )
        if ( *tmp > maxinput )
            maxinput = *tmp;
    if ( _doLogOutput )
        qDebug() << "POSTPROC:calculateIFFT: maxinput" << maxinput;

    if ( maxinput == 0 ) return nullptr;

#ifdef USE_COMPLEX_WEB
    // Jetzt wird der Algorithmus aus  http://paulbourke.net/miscellaneous/dft/   verwendet

    if ( _arrFFTdim != s && _arrFFT != nullptr )
    {
        delete _arrFFT;
        delete _complex;
        _arrFFT = nullptr;
    }
    if ( _arrFFT == nullptr )
    {
        _arrFFTdim = s;
        _arrFFT    = new double[_len];
        _complex   = new COMPLEX*[sy];
        for ( int y=0; y<sy; y++ )
            _complex[y] = new COMPLEX[sx];
    }

    //TODO: Input-Scaling wie unten

    for ( int y=0; y<sy; y++ )
        for ( int x=0; x<sx; x++ )
        {
            _complex[y][x].real = data(x,y) / maxinput;  // Input skalieren auf 0..1
            _complex[y][x].imag = 0.0;
        }

    int rv = FFT2D( _complex, sx, sy, foreward?1:-1 );
    // The direction dir, 1 for forward, -1 for reverse
    if ( !rv )
    {
        qDebug() << "FAIL";
        // TODO ?
    }

    // Calculate the output
    int i;
    switch ( ot )
    {
    case outReal:
        for ( int y=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrFFT[i] = _complex[x][y].real;
            }
        break;
    case outImag:
        for ( int y=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrFFT[i] = _complex[x][y].imag;
            }
        break;
    case outAbs:
        for ( int y=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrFFT[i] = sqrt( sqr(_complex[x][y].real) + sqr(_complex[x][y].imag) );
            }
        break;
    case outSpec:
        for ( int y=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrFFT[i] = (_complex[x][y].real+_complex[x][y].imag) * (_complex[x][y].real-_complex[x][y].imag);
            }
        break;
    }

    scaleAndClipData( _arrFFT, _arrFFTdim*_arrFFTdim, scaled, clip, clip40, false );

    return _arrFFT;

#else

    // Das Array wird hier jetzt lokal angelegt
    if ( _arrFFTdim != s && _arrIFFT != nullptr )
    {
        delete _arrIFFT;
        _arrIFFT = nullptr;
    }
    if ( _arrIFFT == nullptr )
    {
        _arrFFTdim = s;
        _arrIFFT   = new double[s*s];
    }

    fftw_complex *in  = fftw_alloc_complex(_arrFFTdim*_arrFFTdim);
    fftw_complex *out = fftw_alloc_complex(_arrFFTdim*_arrFFTdim);        /* double [2] */

    // Create the calculaton pan for the lib (in,out are modified)
    fftw_plan plan = fftw_plan_dft_2d(_arrFFTdim, _arrFFTdim, in, out,
                                      (foreward?FFTW_FORWARD:FFTW_BACKWARD),
                                      FFTW_ESTIMATE );

    memset( in,  0, _arrFFTdim*_arrFFTdim*2*sizeof(double) );
    memset( out, 0, _arrFFTdim*_arrFFTdim*2*sizeof(double) );

    // Scale the input to the range of [0..1] - and dimensions
    if ( sx < s && sy < s )
    {   // Data smaller than FFT
        maxinput = 0; // wird neu bestimmt
        float fx = float(sx) / float(s);
        float fy = float(sy) / float(s);
        int no = 0;
        for ( int dx=0; dx<s; dx++ )
            for ( int dy=0; dy<s; dy++ )
            {
                int i = dx * s + dy;
                int x = dx * fx;
                int y = dy * fy;
                //qDebug() << dx << dy << "=" << i << "/" << x << y;
                if ( i < 0 || i >= _arrFFTdim*_arrFFTdim ) { no++; continue; }
                if ( x >= sx || y >= sy ) { no++; continue; }
                in[i][0] += data(x,y);
                if ( in[i][0] > maxinput ) maxinput = in[i][0];
            }
        //qDebug() << "FFT: src < dst:" << fx << fy << "max:" << maxinput << "IgnIdx:" << no;
        for ( int i=0; i<_arrFFTdim*_arrFFTdim; i++ )
            in[i][0] /= maxinput;
    }
    else if ( sx > s && sy > s )
    {   // Data larger than FFT
        maxinput = 0; // wird neu bestimmt
        float dx = float(s) / float(sx);
        float dy = float(s) / float(sy);
        int no = 0;
        for ( int x=_xmin; x<_xmax; x++ )
            for ( int y=_ymin; y<_ymax; y++ )
            {
                int i = (int((x-_xmin)*dx) * _arrFFTdim) + int((y-_ymin)*dy);
                //qDebug() << x << y << "/" << x*dx << y*dy << "=" << i << "/" << dx << dy;
                if ( i < 0 || i >= _arrFFTdim*_arrFFTdim ) { no++; continue; }
                in[i][0] += data(x,y);
                if ( in[i][0] > maxinput ) maxinput = in[i][0];
            }
        //qDebug() << "FFT: src > dst:" << dx << dy << "max:" << maxinput << "IgnIdx:" << no;
        for ( int i=0; i<_arrFFTdim*_arrFFTdim; i++ )
            in[i][0] /= maxinput;
    }
    else if ( sx != s || sy != s )
    {   // Data size don't match FFT size - IGNORED
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
        return nullptr;
    }
    else
    {   // Data same size as FFT
        for ( int i=0, x=_xmin; x<_xmax && i<_arrFFTdim*_arrFFTdim; x++ )
            for ( int y=_ymin; y<_ymax && i<_arrFFTdim*_arrFFTdim; y++, i++ )
            {
                in[i][0] = data(x,y) / maxinput;
            }
        //qDebug() << "FFT: src = dst:" << maxinput;
    }

    // Calculate the IFFT
    fftw_execute(plan);

    // Calculate the output
    int i;
    switch ( ot )
    {
    case outReal:
        for ( int y=0, k=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++, k++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrIFFT[i] = out[k][0];
            }
        break;
    case outImag:
        for ( int y=0, k=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++, k++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrIFFT[i] = out[k][1];
            }
        break;
    case outAbs:
        for ( int y=0, k=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++, k++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrIFFT[i] = sqrt(out[k][0]*out[k][0]+out[k][1]*out[k][1]);
            }
        break;
    case outSpec:
        for ( int y=0, k=0; y<_arrFFTdim; y++ )
            for ( int x=0; x<_arrFFTdim; x++, k++ )
            {
                i = calculateIndex( swapout, x, y, _arrFFTdim );
                _arrIFFT[i] = (out[k][0]+out[k][1]) * (out[k][0]-out[k][1]);
            }
        break;
    }

    scaleAndClipData( _arrIFFT, _arrFFTdim*_arrFFTdim, scaled, clip, clip40, false );

    fftw_destroy_plan(plan);

    fftw_free(in);
    fftw_free(out);

    fftw_cleanup();

    return _arrIFFT;
#endif // !USE_COMPLEX_WEB
}

/**
 * @brief SasCalc_PostProc::calculateIndex
 * @param swap - flag swap output (0° in center)
 * @param x    - input x index
 * @param y    - input y index
 * @return calculated index in data array
 */
int SasCalc_PostProc::calculateIndex(bool swap, int x, int y, int dim)
{
    if ( swap )
    {
        if ( y < dim/2 && x < dim/2 )
        {   // Links unten
            return (dim/2 - x) + dim*(dim/2 - y);
        }
        if ( y >= dim/2 && x < dim/2 )
        {   // Links oben
            return (dim/2 - x) + dim*(dim-1 - y + dim/2);
        }
        if ( y < dim/2 && x >= dim/2 )
        {   // Rechts unten
            return (dim-1 - x + dim/2) + dim*(dim/2 - y);
        }
        // if ( y >= dim/2 && x >= dim/2 )
        {   // Rechts oben
            return (dim-1 - x + dim/2) + dim*(dim-1 - y + dim/2);
        }
    }
    return x + dim*y; // Default ohne Swap
}



/**
 * @brief SasCalc_PostProc::scaleAndClipData
 * @param data   - data pointer (in/out)
 * @param len    - length of data in values
 * @param scale  - flag scaling [0,1] lin
 * @param clip   - flag clipping [1e-10,1] log
 * @param clip40 - flag clipping [40%,100%] -> [0,1]
 * @param genlog - flag generate log data
 */
void SasCalc_PostProc::scaleAndClipData( double *data, int len, bool scale, bool clip, bool clip40, bool genlog )
{
    if ( scale )
    {   // Scale to max
        double max = 0;
        double *ptr = data;
        for ( int i=0; i<len; i++, ptr++)
            if ( max < *ptr )
                max = *ptr;
        ptr = data;
        for ( int i=0; i<len; i++, ptr++ )
            *ptr = *ptr / max;
    }
    if ( clip )
    {   // Clip auf [1e-6 .. 1]
        double *ptr = data;
        for ( int i=0; i<len; i++, ptr++)
        {
            *ptr = ( ( log( abs(*ptr) / 1e-20 ) / log(10) ) + 6.0 ) / 6.0;
            if ( *ptr < 0 ) *ptr = 0;
        }
    }
    if ( clip40 )
    {   // Clip auf 40% .. 100%
        double *ptr = data;
        for ( int i=0; i<len; i++, ptr++ )
        {
            *ptr = 1.0 - ( 1.0 - *ptr ) / (1.0 - 0.4);
            if ( *ptr < 0 ) *ptr = 0;
        }
    }
    if ( genlog )
    {   //
        double *ptr = data;
        double vlogmin=1e26, vlogmax=0;
        for ( int i=0; i<len; i++, ptr++)
            if ( *ptr > 0 )
            {
                if ( log10(*ptr) < vlogmin )
                    vlogmin = log10(*ptr);
                else if ( log10(*ptr) > vlogmax )
                    vlogmax = log10(*ptr);
            }
        ptr = data;
        for ( int i=0; i<len; i++, ptr++)
            if ( *ptr > 0 )
                *ptr = (log10(*ptr) - vlogmin) / (vlogmax-vlogmin);
            else
                *ptr = 0;
    }
}



#ifdef USE_COMPLEX_WEB
/*#############################################################*/
/*   http://paulbourke.net/miscellaneous/dft/                  */
/*#############################################################*/

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
int FFT2D(COMPLEX **c,int nx,int ny,int dir)
{
   int i,j;
   int m,twopm;
   double *real,*imag;

   /* Transform the rows */
   real = (double *)malloc(nx * sizeof(double));
   imag = (double *)malloc(nx * sizeof(double));
   if (real == NULL || imag == NULL)
      return(FALSE);
   if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return(FALSE);
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         real[i] = c[i][j].real;
         imag[i] = c[i][j].imag;
      }
      FFT(dir,m,real,imag);
      for (i=0;i<nx;i++) {
         c[i][j].real = real[i];
         c[i][j].imag = imag[i];
      }
   }
   free(real);
   free(imag);

   /* Transform the columns */
   real = (double *)malloc(ny * sizeof(double));
   imag = (double *)malloc(ny * sizeof(double));
   if (real == NULL || imag == NULL)
      return(FALSE);
   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return(FALSE);
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         real[j] = c[i][j].real;
         imag[j] = c[i][j].imag;
      }
      FFT(dir,m,real,imag);
      for (j=0;j<ny;j++) {
         c[i][j].real = real[j];
         c[i][j].imag = imag[j];
      }
   }
   free(real);
   free(imag);

   return(TRUE);
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/
int FFT(int dir,int m,double *x,double *y)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

   return(TRUE);
}

/*-------------------------------------------------------------------------
   Calculate the closest but lower power of two of a number
   twopm = 2**m <= n
   Return TRUE if 2**m == n
*/
int Powerof2(int n,int *m,int *twopm)
{
   if (n <= 1) {
      *m = 0;
      *twopm = 1;
      return(FALSE);
   }

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
      return(FALSE);
   else
      return(TRUE);
}
#endif
