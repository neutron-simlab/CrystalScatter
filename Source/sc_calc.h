#ifndef SC_CALC_H
#define SC_CALC_H

#include <QString>
#include <QSettings>
#include "sc_math.h"
//#include "sc_memory_gpu.h"
#include "sc_globalConfig.h"
#include <string>
#include <list>



/**
 * @brief data transfer helper
 */
typedef struct
{
    bool    checked;
    int     select;
    QString str;
    double  value;
    Double3 vec;
} _valueTypes;

typedef bool (*_dataGetter)( QString, _valueTypes& );



/**
 * @brief The SC_Calc class is the baseclass for all calculation methods.
 */
class SC_Calc
{
public:
    SC_Calc() {}

    /**
     * @brief methodName
     * @return the method name to be displayed in the GUI
     */
    virtual QString methodName() = 0;

    /**
     * @brief guiLayout
     * @return a list of definitions for the parameter GUI
     * Each element must be in the form "x;y;prompt;type;tooltip;default" with:
     *  x;y     = index in the grid (0,1,2,....)
     *  prompt  = prompting text label left of the inputfield
     *  type    = Selection : "cbs|...|...|..."
     *             Fittable : "cbsfit|...|...|..."
     *                        with the list of possible values
     *            Textinput : "txt|len"
     *                        with the maximum character count (optional, default: 16k)
     *            Numericals: "inp|frac|min|max|unit"
     *              Fittable: "inpfit|frac|min|max|unit"
     *                        with the fraction count and the limits (optional, default: 2, -10000, +10000)
     *            CheckBox  : "tog"
     *                        only true/false possible
     *            Infolabel : "lbl"
     *                        this expands the prompt over both columns and has no value
     *  tooltip = this is the tooltip set to the inputfield (optional)
     *  default = the default value (optional)
     * Space characters at the beginning and end of each field are discarded before use.
     */
    virtual QStringList guiLayout() = 0;

    /**
     * @brief prepareData
     * @param _dataGetter - function pointer to get all data values
     * Called during data preparation before the calculation can start and let the
     * method function retrieve the data from the GUI calss.
     */
    virtual void prepareData( _dataGetter ) = 0;

    /**
     * @brief doCalculation
     * @param numThreads   - is the number of used threads if no GPU is found
     * @param pa           - pointer to global function to update the progess bar and to ask the Abort button
     * Starts the calculation for this method.
     */
    virtual void doCalculation( int, progressAndAbort ) = 0;
    virtual double doFitCalculation( int, int, int, long&, long& ) = 0;

    virtual std::string tpvPerformRandom( std::list<std::string> ids ) = 0;

    /**
     * @brief higResTimerElapsed
     * @return the elapsed time in milliseconds
     */
    virtual double higResTimerElapsedPrep() = 0;

    /**
     * @brief higResTimerElapsed
     * @return the elapsed time in milliseconds
     */
    virtual double higResTimerElapsedCalc() = 0;

    virtual void endThread() = 0;

    virtual void cleanup() = 0;
    virtual bool gpuAvailable() = 0;
    virtual int minX() = 0;
    virtual int maxX() = 0;
    virtual int minY() = 0;
    virtual int maxY() = 0;
    virtual double *data() = 0;
    virtual double xyIntensity( int x, int y ) = 0;
    virtual int lastX() = 0;
    virtual int lastY() = 0;
#ifdef COPY_FITDATA_TO_GPU  // FITDATA_IN_GPU ok, virtual func def
    virtual bool setArrDataForFit( const double *data ) = 0;
#ifdef CALC_FQS_IN_GPU
    virtual double getFQS() = 0;
#endif
#endif
#ifdef FITDATA_IN_GPU  // virtual func def
    virtual bool setFitData( int sx, int sy, const double *data ) = 0;
#endif
    virtual void setNoFitRect( int id, int x0, int y0, int x1, int y1 ) = 0;

};

#endif // SC_CALC_H
