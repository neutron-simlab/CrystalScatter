#include "sc_simplexfit2d.h"
#include <QDebug>

// Code von Herrn Prof. Förster
// Erklärungen auch unter https://de.wikipedia.org/wiki/Downhill-Simplex-Verfahren


#define DBGMODIFY(x) //x
// Ich habe bei Tests zufällig gesehen, dass bei einem Lauf durch FQS die Werte von allen
// Parametern geändert werden. Ich weiß nicht, ob das wirklich so gewollt ist. Mit dieser
// Definition umklammere ich alle Outputs zu diesem Thema. Es sollte dann auch der Parameter
// bei calc->updateParamValueForFit() auf true gesetzt werden, um die Werteänderunge zu sehen.



double SasCalc_SimplexFit2D::fehlerquadratsumme( int numThreads, double *params, QString &info )
{
    // Im Pascalprogramm sind alle Parameter in lokalen Variablen, daher hier die Werte aus dem Parameter nehmen
    // und somit keine globalen Daten.

    // (* Randbedingungen für Parameter *)
    for ( int j=0; j<numFit; j++ )
    {
        if ( enabled[j] )
        {
            // Versuch: Winkelangaben (0..360) hier nicht prüfen sondern durchlassen.
            if ( minVals[j] == 0.0 && maxVals[j] == 360.0 )
                continue;

            if ( params[j] < minVals[j] )   // vorher <= (wenn Val=0 und Min=0 folgte immer ein Fehler)
            {
                info = QString("MIN Index:%1=%2<%3").arg(indVec[j]).arg(params[j]).arg(minVals[j]);
                return 1e20;
            }
            if ( params[j] > maxVals[j] )   // vorher >= (mal sehen, ob es klappt)
            {
                info = QString("MAX Index:%1=%2>%3").arg(indVec[j]).arg(params[j]).arg(minVals[j]);
                return 1e20;
            }
        }
    }
    for ( int j=0; j<mmax && j<indVec.size(); j++ ) // all params
    {
        calc->updateParamValueForFit( indVec[j], params[j], false/*ohne Ausgaben, wenn Werte geändert werden*/ );
    }
    calc->prepareCalculation( "F", false );     // to transfer parameter values into the calculation class

    // Berechne das Image...
#ifdef FITDATA_IN_GPU  // fqs
    if ( useGpuFit )
    {   // Neuer Algorithmus (Juli 2022), damit die FQS auf der GPU berechnet wird
        long   cnt = 0;     // Pixelcounter for Debug
        long   nancnt = 0;  // Counter for nan values
        double summe = calc->doFitCalculation(numThreads,_bstoppixel,_borderPixel,cnt,nancnt);
        numImgCalc++;
        info = QString("(%1px/%2px)").arg(cnt).arg((_xmax-_xmin)*(_ymax-_ymin));
        if ( nancnt > 0 ) info += QString(" NaN:%1").arg(nancnt);
        return summe;
    }
#endif
    calc->doCalculation(numThreads,nullptr);    // calculation without progressBar
    numImgCalc++;
    // TODO qApp->processEvents();

#ifdef COPY_FITDATA_TO_GPU
#ifdef CALC_FQS_IN_GPU
    //double fqs_gpu = -10.0;
    if ( useGpuForMask )
    {
        info = QString("(gpu/%1px)").arg((_xmax-_xmin)*(_ymax-_ymin));
        //fqs_gpu = calc->getFQS();
        return calc->getFQS();
    }
#endif
#endif

    double summe = 0.0;
    long   cnt = 0;     // Pixelcounter for Debug
    long   nancnt = 0;  // Counter for nan values
    long   nullcnt = 0; // Counter for null values after beamstop check

    if ( _borderPixel > 0 || _bstoppixel > 0 )
    {   // Ausblendungen über Pixelangaben an Rand und Mitte
        for ( int ihex=calc->minX(), ihexd=_xmin; ihex<calc->maxX() && ihexd<_xmax; ihex++, ihexd++ )  // (***   z-loop  ***)
        {
            // Border
            //if ( ihex  < calc->minX()+_borderPixel || ihex  >= calc->maxX()-_borderPixel ) continue;
            if ( ihexd < _xmin       +_borderPixel || ihexd >= _xmax       -_borderPixel ) continue;
            // loop over y
            for ( int i=calc->minY(), id=_ymin; i<calc->maxY() && id<_ymax; i++, id++ )
            {
                // Border
                //if ( i  < calc->minY()+_borderPixel || i  >= calc->maxY()-_borderPixel ) continue;
                if ( id < _ymin       +_borderPixel || id >= _ymax       -_borderPixel ) continue;
                // Beamstop (Center)
                if ( ihexd >= _xbs-_bstoppixel && ihexd < _xbs+_bstoppixel &&
                     id    >= _ybs-_bstoppixel && id    < _ybs+_bstoppixel ) continue;
                // Test: Werte <= 0.0 werden ignoriert, da der log10(0.0) keinen Sinn ergibt
                if ( calc->getMemPtr()->xyIntensity(ihex,i) <= 0.0 ||
                     fitData(ihexd,id) <= 0.0 )
                {
                    nullcnt++;
                    continue;
                }
                // Pixelcounter for Debug
                cnt++;
                summe += FQSVERGL( calc->getMemPtr()->xyIntensity(ihex,i),   // oben berechnetes Pixel
                                   fitData(ihexd,id) );                      // anzufittendes Pixel
                if ( std::isnan(calc->getMemPtr()->xyIntensity(ihex,i)) )
                    nancnt++;
            } // for i
        } // for ihex
    }
    else if ( _bstoppixel == -1 )
    {   // Ausblendungen per Eck-Pixel-Wert
        double ecke = fitData(_xmin,_ymin);
        for ( int ihex=calc->minX(), ihexd=_xmin; ihex<calc->maxX() && ihexd<_xmax; ihex++, ihexd++ )  // (***   z-loop  ***)
        {
            for ( int i=calc->minY(), id=_ymin; i<calc->maxY() && id<_ymax; i++, id++ )
            {
                // Special check: ignore pixel if source is less equal corner pixel
                if ( fitData(ihexd,id) <= ecke ) continue;
                // Test: Werte <= 0.0 werden ignoriert, da der log10(0.0) keinen Sinn ergibt
                if ( calc->getMemPtr()->xyIntensity(ihex,i) <= 0.0 ||
                     fitData(ihexd,id) <= 0.0 )
                {
                    nullcnt++;
                    continue;
                }
                // Pixelcounter for Debug
                cnt++;
                summe += FQSVERGL( calc->getMemPtr()->xyIntensity(ihex,i),   // oben berechnetes Pixel
                                   fitData(ihexd,id) );                      // anzufittendes Pixel
/*
                if ( numImgCalc == 2 ) // && ihex==-30 && i==-30 )
                    qDebug() << "cpu" << ihex << i
                             << "idx:" << (-_xmin + (ihexd)) + (_xmax-_xmin)*(-_ymin + (id))
                             << "calc:" << calc->getMemPtr()->xyIntensity(ihex,i)
                             << "fit:" << fitData(ihexd,id);
                             //<< "x:" << _xmin << _xmax << "y:" << _ymin << _ymax
                             //<< "#" << numImgCalc;
*/
                if ( std::isnan(calc->getMemPtr()->xyIntensity(ihex,i)) )
                    nancnt++;
            } // for i
        } // for ihex
    }
    else
    {   // Keine Ausblenungen
        for ( int ihex=calc->minX(), ihexd=_xmin; ihex<calc->maxX() && ihexd<_xmax; ihex++, ihexd++ )  // (***   z-loop  ***)
        {
            for ( int i=calc->minY(), id=_ymin; i<calc->maxY() && id<_ymax; i++, id++ )
            {
                // Test: Werte <= 0.0 werden ignoriert, da der log10(0.0) keinen Sinn ergibt
                if ( calc->getMemPtr()->xyIntensity(ihex,i) <= 0.0 ||
                     fitData(ihexd,id) <= 0.0 )
                {
                    nullcnt++;
                    continue;
                }
                // Pixelcounter for Debug
                cnt++;
                summe += FQSVERGL( calc->getMemPtr()->xyIntensity(ihex,i),   // oben berechnetes Pixel
                                   fitData(ihexd,id) );                      // anzufittendes Pixel
                if ( std::isnan(calc->getMemPtr()->xyIntensity(ihex,i)) )
                    nancnt++;
            } // for i
        } // for ihex
    }

    info = QString("(%1px/%2px)").arg(cnt).arg((_xmax-_xmin)*(_ymax-_ymin));
    if ( nancnt > 0 ) info += QString(" NaN:%1").arg(nancnt);
    if ( nullcnt > 0 ) info += QString(" NULL:%1").arg(nullcnt);
//#ifdef COPY_FITDATA_TO_GPU
//    if ( fqs_gpu >= 0 ) info += QString(" GPU=%1").arg(fqs_gpu);
//#endif
    return summe;
}


double SasCalc_SimplexFit2D::getResiduenPixel( int ihex, int i )
{
    int ihexd = ihex - calc->minX() + _xmin;
    int id    = i    - calc->minY() + _ymin;

    // Border
    if ( ihexd < _xmin       +_borderPixel || ihexd >= _xmax       -_borderPixel ) return 0;
    if ( id    < _ymin       +_borderPixel || id    >= _ymax       -_borderPixel ) return 0;
    // Beamstop (Center)
    if ( _bstoppixel > 0 &&
         ihexd >= _xbs-_bstoppixel && ihexd < _xbs+_bstoppixel &&
         id    >= _ybs-_bstoppixel && id    < _ybs+_bstoppixel ) return 0;
    // Special check: ignore pixel if source is less equal corner pixel
    if ( _bstoppixel == -1 &&
         fitData(ihexd,id) <= fitData(_xmin,_ymin) ) return 0;
    // Test: Werte <= 0.0 werden ignoriert, da der log(0.0) keinen Sinn ergibt
    if ( calc->getMemPtr()->xyIntensity(ihex,i) <= 0.0 ||
         fitData(ihexd,id) <= 0.0 ) return 0;
    // Log(I(Rechnung) – Log(I(Experiment)
    return  log(calc->getMemPtr()->xyIntensity(ihex,i)) - log(fitData(ihexd,id));
}


// (* function amotry (Simplex) *)
double SasCalc_SimplexFit2D::amotry( int numThreads,
                                     double p[nap][map],
                                     double *y,
                                     double *sum,
                                     int ihi, double fac,
                                     QString &info )
{
    double ptry[map];
    for ( int j=0; j<mmax; j++ ) ptry[j]=sum[j]; // Alle kopieren

    double fac1=(1.0-fac)/mmax; // (* Faktor mit fac = -alfa, beta, -gamma; o.k. *)
    double fac2=fac1-fac;       // (* Faktor mit center anstatt centroid; o.k. *)
    for ( int j=0; j<numFit; j++ ) // nur die variablen Werte
    {
        DBGMODIFY( qDebug() << "AMOTRY" << j << sum[j] << fac1 << ihi << p[ihi][j] << fac2; )
        ptry[j] = sum[j] * fac1 - p[ihi][j] * fac2;   // (* Berechnung des Testvertex, o.k. *)
    }
    double ytry=fehlerquadratsumme(numThreads,ptry,info);  // (* Berechnung der Fehlerquadratsumme des Testvertex; o.k. *)
    info += QString(" %1<y[%2]=%3?").arg(ytry).arg(ihi).arg(y[ihi]);
    if ( ytry < y[ihi] )
    {                // (* Testvertex besser als schlechtester *)
        y[ihi]=ytry; // (* dann erstmal den schlechtesten ersetzen *)
        for ( int j=0; j<numFit; j++ )
        {
            sum[j] += ptry[j]-p[ihi][j]; // (* Centrum schon mal umrechnen *)
            p[ihi][j]=ptry[j];           // (* Vertex durch Testvertex ersetzen *)
            DBGMODIFY( qDebug() << "amotry set" << ihi << j << "=" << p[ihi][j]; )
        }
    }
    return ytry; // (* gibt y = FQS des Testvertex zurück *)
}


void SasCalc_SimplexFit2D::doSimplexFit2D(int numThreads, double stp, int maxit, double ftol,
                                          int borderpix, int bstoppix, progressLogging pl,
                                          _param2fitval *vals, QString &retinfo )
{
    aborted = false;
    retinfo = "";
#ifdef COPY_FITDATA_TO_GPU
    useGpuForMask = false;
#endif

    parVals = vals;
    //typedef enum { fitNone, fitNumeric, fitCbs } _fitTypes;
    //typedef struct
    //{
    //    double min, max,    // Limits for the fit
    //           fitstart,    // Startvalue for the fit
    //           fitres;      // Resultvalue after fitting
    //    _fitTypes fitType;  // fitCbs has only integer steps
    //    bool used;          // true=used for the fit, false=not changed during fit
    //    bool fitvalid;      // Flag, if .fitres is valid
    //} _fitLimits;
    //typedef QHash< QString/*name*/, _fitLimits* > _param2fitval;

    progLogging  = pl;
    numImgCalc   = 0;
    _borderPixel = borderpix;
    _bstoppixel  = bstoppix;

    int i=0;
    indVec.clear();
    QStringList keys = parVals->keys();
    keys.sort();
    // Dadurch bleiben die Parameternamen immer in der gleichen Reihenfolge.
    foreach ( QString k, keys )
    {
        _fitLimits *fl = parVals->value(k);
        if ( fl->used )
        {
            indVec << k; // it.key();
            enabled[i] = true;
            values[i]  = fl->fitstart; // it.value()->fitstart;
            minVals[i] = fl->min;      // it.value()->min;
            maxVals[i] = fl->max;      // it.value()->max;
            //qDebug() << "Var" << i << k << values[i] << minVals[i] << maxVals[i];
            if ( ++i >= nap ) break;
        }
    }
    numFit = indVec.size();   // Number of fitable parameters
    DBGMODIFY( qDebug() << "Fittable" << numFit << indVec; )

    // Jetzt noch die festen Parameter übernehmen
    foreach ( QString k, keys )
    {
        _fitLimits *fl = parVals->value(k);
        if ( ! fl->used )
        {
            indVec << k; // it.key();
            enabled[i] = false;
            values[i]  = fl->fitstart; // it.value()->fitstart;
            minVals[i] = fl->min;      // it.value()->min;
            maxVals[i] = fl->max;      // it.value()->max;
            //qDebug() << "Fix" << i << k << values[i] << minVals[i] << maxVals[i];
            if ( ++i >= nap ) break;
        }
    }
    mmax = indVec.size();   // Number of parameters in total
    DBGMODIFY( qDebug() << "All Params" << mmax << indVec; )
    if ( mmax >= nap )
    {
        retinfo = QString("The maximum array length (%1) was exeeded by the current number of fitable parameters (%2)")
                  .arg(nap).arg(mmax);
        return;
    }
    //qDebug() << "Fit: indVec=" << indVec << "mmax,numFit:" << mmax << numFit;

#ifndef CONSOLENPROG
#ifdef FITDATA_IN_GPU  // init fit
    useGpuFit = calc->setFitData( _xmax - _xmin, _ymax - _ymin, intensityForFit );
    qDebug() << "Fit: useGpu =" << useGpuFit;
#endif
#endif

    QString info;   // For informations from fehlerquadratsumme()

    auto start   = std::chrono::high_resolution_clock::now();

    const double alfa  = 1.0; // (* simplex reflection *)
    const double beta  = 0.5; // (* simplex contraction *)
    const double gamma = 2.0; // (* simplex expansion *)

    double y[nap];     // y: ArrayN;  (* ArrayN = array[1..nap=31] of extended; *)
    double pp[map];    // pp: arrayM;
    double q[map];     // q: arrayM;
    double step[map];  // step: arrayM;
    double psum[map];  // psum: ^ArrayM; (* ArrayM = RealArrayMA = array[1..map=30] of extended; *)
    double nmax = mmax + 1.0;
    int    ndim = mmax;
    int    mpts = numFit /*ndim*/ + 1;
    // p[1,i] = pars[i] = globalParams.values[i]
    double p[nap][map];
    char logbuffer[300];

    double rootn1 = sqrt(nmax) - 1.0; // (* Schrittweite *)
    double mroot2 = 1.0 / (mmax*sqrt(2.0));
    //double m1 = 1.0/mmax;
    for ( int i=0; i<mmax; i++ ) // Hier alle Werte nutzen
    {
        psum[i] = values[i];
        p[0][i] = values[i]; // (* erster Vertex p[1,...] *)
        if ( i < numFit )
        {
            if ( values[i] == 0 )
            {   // Startwert Null macht wenig Sinn
                if ( minVals[i] == 0 )
                    step[i] = stp;  // Min==0 nehmen wir einfach einmal die Schrittweite
                else
                    step[i] = minVals[i] * stp;
            }
            else
                step[i] = values[i] * stp;
        }
        else
            step[i] = 0.0;
        pp[i]   = step[i] * (rootn1+mmax)*mroot2;
        q[i]    = step[i] * rootn1*mroot2;
        //if ( enabled[i] )
        //    qDebug() << indVec[i] << "Val" << values[i] << "Step" << step[i]
        //                << "pp,q,psum" << pp[i] << q[i] << psum[i];
    }
    for ( int i=mmax; i<map; i++ )
    {
        psum[i] = 0.0;
        p[0][i] = 0.0;
        step[i] = 0.0;
        pp[i]   = 0.0;
        q[i]    = 0.0;
    }
    for ( int i=0; i<nap; i++ )
    {
        y[i] = 0.0;
        for ( int j=0; i<map; i++ )
            p[i][j] = 0.0;
    }

    DBGMODIFY( qDebug() << "Fit: 0 fqs"; )
    y[0] = fehlerquadratsumme(numThreads,values,info);

    //return;

    if ( progLogging )
    {
        if ( y[0] > 1e10 )
            sprintf( logbuffer, "FQS(0) = %s", qPrintable(info) );
        else
            sprintf( logbuffer, "FQS(0) = %lf  %s", y[0], qPrintable(info) );
        if ( progLogging(logbuffer) ) return;
        DBGMODIFY( qDebug() << logbuffer; )
    }

#ifdef COPY_FITDATA_TO_GPU
    // Da der GPU-Speicher (auch intensityForFit) erst bei der ersten Berechnung
    //  angelegt wird, erst danach testen, ob diese Option für uns passt.
    if ( _bstoppixel == -1 )
    {
        if ( !useGpuForMask )
        {
            useGpuForMask = calc->setArrDataForFit( intensityForFit );
#ifndef __CUDACC__
            qDebug() << "ArrDataForFit:" << useGpuForMask;
#endif
        }
    }
    else
        useGpuForMask = false;
#endif

    for ( int i=1; i<nmax; i++ ) // (* restliche Vertices p[i,...] *)
    {
        for ( int j=0; j<mmax; j++ )
            p[i][j] = p[0][j] + q[j];
        // predi:=pred(i);  -->  predi = i-1
        p[i][i-1] = p[i][i-1] + pp[i-1];
    }

    for ( int i=1; i<nmax; i++ )
    {
        for ( int j=0; j<numFit; j++ )
            values[j] = p[i][j];
        if ( i > numFit+1 )       // This speeds up a little bit with GPU
            y[i] = y[numFit+1];   // but it saves more time without GPU
        else
        {
            DBGMODIFY( qDebug() << "Fit:" << i << "fqs"; )
            y[i] = fehlerquadratsumme(numThreads,values,info);
        }
        if ( progLogging )
        {
            if ( y[i] > 1e10 )
                sprintf( logbuffer, "FQS(%d) = %s", i, qPrintable(info) );
            else
                sprintf( logbuffer, "FQS(%d) = %lf  %s", i, y[i], qPrintable(info) );
            if ( progLogging(logbuffer) ) return;
            DBGMODIFY( qDebug() << logbuffer; )
            /*if ( info.contains("NULL:") )
            {
                retinfo = info;
                return;
            }*/
        }
    }

    for ( int j=0; j<ndim; j++ ) // (* center *)
    {
        double sum = 0.0;
        for ( int i=0; i<mpts; i++ ) sum += p[i][j];
        psum[j] = sum;
    }

    //(*********************  Hauptschleife zur Suche *****************)
    //(****************************************************************)

    int ilo=0, ihi, inhi;
    for ( repetitions=0; repetitions<maxit; repetitions++ )
    {
        if ( progLogging )
        {   // Update the label in GUI
            sprintf( logbuffer, "@REP=%d", repetitions+1 );
            if ( progLogging(logbuffer) ) break; //return;
            DBGMODIFY( qDebug() << logbuffer; )
        }

        ilo = 0;
        if ( y[0] > y[1] )
        {
            ihi=0;
            inhi=1;
        }
        else
        {
            ihi=1;
            inhi=0;
        }
        for ( int i=0; i<mpts; i++ )
        {
            if ( y[i] < y[ilo] ) ilo=i;
            if ( y[i] > y[ihi] )
            {
                inhi=ihi;
                ihi=i;
            }
            else if ( y[i] > y[inhi] )
            {
                if ( i != ihi ) inhi=i;
            }
        }

        if ( fabs(y[ihi]) + fabs(y[ilo]) == 0 )
        {   // ==> rtol=nan, makes no sense...
            aborted = true;
            qDebug() << "ABORT rtol=NaN";
            break;
        }

        // rtol := 2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo]));
        double rtol = 2.0 * fabs( y[ihi] - y[ilo] ) / ( fabs(y[ihi]) + fabs(y[ilo]) );

        if ( isnan(rtol) )
            qDebug() << "rtol=NAN" << ihi << ilo << y[ihi] << y[ilo];

        if ( progLogging )
        {
            sprintf( logbuffer, "----- rtol=%lg <= ftol=%lg", rtol, ftol );
            if ( progLogging(logbuffer) ) break; //return;
            DBGMODIFY( qDebug() << logbuffer; )
        }
        if ( rtol <= ftol ) break;

        if ( progLogging )
        {
            sprintf( logbuffer, "iter=%d, ilo=%d, ihi=%d, inhi=%d", repetitions, ilo, ihi, inhi );
            if ( progLogging(logbuffer) ) break; //return;
            DBGMODIFY( qDebug() << logbuffer; )
        }

        DBGMODIFY( qDebug() << "... amotry -alpha"; )
        //     ytry = amotry(ma,p,y,psum^,diff,xydata,xybdata,ndim,ihi,-alfa);
        double ytry = amotry( numThreads, p, y, psum, ihi, -alfa, info );
        if ( progLogging )
        {
            if ( ytry > 1e10 )
                sprintf( logbuffer, "amotry(simplex reflection) = %s", qPrintable(info) );
            else
                sprintf( logbuffer, "amotry(simplex reflection) = %lf  %s", ytry, qPrintable(info) );
            if ( progLogging(logbuffer) ) break; //return;
            DBGMODIFY( qDebug() << logbuffer; )
        }

        if ( ytry <= y[ilo] )
        {
            DBGMODIFY( qDebug() << "... amotry gamma"; )
            //ytry = amotry(ma,p,y,psum^,diff,xydata,xybdata,ndim,ihi,gamma);
            ytry = amotry( numThreads, p, y, psum, ihi, gamma, info );
            if ( progLogging )
            {
                if ( ytry > 1e10 )
                    sprintf( logbuffer, "amotry(simplex expansion) = %s", qPrintable(info) );
                else
                    sprintf( logbuffer, "amotry(simplex expansion) = %lf  %s", ytry, qPrintable(info) );
                if ( progLogging(logbuffer) ) break; //return;
                DBGMODIFY( qDebug() << logbuffer; )
            }
        }
        else if ( ytry >= y[inhi] )
        {
            double ysave=y[ihi];
            DBGMODIFY( qDebug() << "... amotry beta"; )
            //ytry = amotry(ma,p,y,psum^,diff,xydata,xybdata,ndim,ihi,beta);
            ytry = amotry( numThreads, p, y, psum, ihi, beta, info );
            if ( progLogging )
            {
                if ( ytry > 1e10 )
                    sprintf( logbuffer, "amotry(simplex contraction) = %s", qPrintable(info) );
                else
                    sprintf( logbuffer, "amotry(simplex contraction) = %lf  %s", ytry, qPrintable(info) );
                if ( progLogging(logbuffer) ) break; //return;
                DBGMODIFY( qDebug() << logbuffer; )
            }
            if ( ytry >= ysave )
            {
                for ( int i=0; i<mpts; i++ )    // Anzahl der Fittable Parameter
                    if ( i != ilo )
                    {
                        for ( int j=0; j<ndim; j++ )    // ALLE Parameter
                        {
                            psum[j] = 0.5*(p[i][j]+p[ilo][j]);
                            p[i][j] = psum[j];
                        }
                        DBGMODIFY( qDebug() << "..." << i << "fqs" << mpts << ndim; )
                        //y[i] = fqs(ndata,xydata,xybdata,psum^,diff);
                        y[i] = fehlerquadratsumme(numThreads,psum,info);
                        if ( progLogging )
                        {
                            if ( y[i] > 1e10 )
                                sprintf( logbuffer, "FQS(%d) = %s", i, qPrintable(info) );
                            else
                                sprintf( logbuffer, "FQS(%d) = %lf  %s", i, y[i], qPrintable(info) );
                            if ( progLogging(logbuffer) ) { aborted=true; break; } //return;
                        }
                    }
                if ( aborted ) break; // for repetitions
                for ( int j=0; j<ndim; j++ )    // Hier werden ALLE Parameter manipuliert ???
                {
                    double sum=0.0;
                    for ( int i=0; i<mpts; i++ )
                        sum += p[i][j];
                    DBGMODIFY( qDebug() << "..." << j << "update psum[]" << ndim << mpts; )
                    psum[j]=sum;
                }
            }
        }
    } // for rep

    for ( int i=0; i<mmax; i++ )
    {
        // Ergebis zurückgeben
        parVals->value( indVec[i] )->fitres = p[ilo][i];
        //calc->updateParamValue( "", indVec[i], QString::number(p[ilo][i]) ); -> nicht die aktuellen Daten übernehmen
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast<std::chrono::duration<float>>(end-start);
    higResTimerElapsed = calc_time.count()*1000.0;
}
