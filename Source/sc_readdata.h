#ifndef SC_READDATA_H
#define SC_READDATA_H

#include "widimage.h"
#ifndef NOHDF5
#include "H5Cpp.h"
using namespace H5;
#endif



typedef enum {
    dftFromExt,         // Aus der Extension erkennen
    dft2dTiff,          // "2d:TIFF-Images (*.tif *.tiff)"
    dft2dKwsData,       // "2d:KWS-Data (*.dat *.data)"
    dft2dEsrf,          // "2d:ESRF/Klora images (*.edf)"
    dft2dSpreadsheet,   // "2d:Spreadsheet-Format (*.spr)";
    dft2dSans,          // "2d:2D SANS (*.xml *.csv *.dat)";    // Hier kein NXS (ist HDF5, Nexus)
#ifndef NOHDF5
    dft2dHdf5,          // "2d:HDF5 Files (*.h5 *.hdf *.hdf5)";
    dft2dSansNexus,     // "2d:SANS (*.xml *.csv *.nxs *.dat)";
#endif
    dft2dSasCrystal,    // "2d:SasCrystal (*.dat *.csv *.txt)";
    dft1dTextual,       // "1d:Textual data (*.txt *.dat)"; // filetype:=1  TODO
    dft1dCsv,           // "1d:CSV data (*.csv)"; // filetype:=2
    dft1dHeimoPdh,      // "1d:Heimo-PDH-Format (*.pri *.mul *.sub *.fun *.pdh)"; // filetype:=3
    dft1dXRay,          // "1d:X-Ray-Format (*.asc)"; // filetype:=4
    dft1dAcXrd,         // "1d:AC-XRD-Format (*.uxd)"; // filetype:=5
    dft1dFit2dOut,      // "1d:Fit2D output (*.chi)"; // filetype:=6
    dft1dEsrf,          // "1d:ESRF-data (*.out)"; // filetype:=7
    dft1d2DinData,      // "1d:2dim-data (*.csw)"; // filetype:=8
    dftAllFiles         // "All files (*.*)";
} _dataFileTypes;



typedef widImage* (*myAddImage)( int x0, int x1, int y0, int y1, double *d, QString title );


class SC_ReadData
{
public:
    static widImage *readImageTiff( myAddImage add, QString fn );
    static void swapBytes( char *c, unsigned int len );

    static widImage *readImageSasCrystal( myAddImage add, QString fn );

    static widImage *readImageKWSData( myAddImage add, QString fn );
    static void findBeamCenter( widImage *ikws, int &xerg, int &yerg );

    static widImage *readImage2dSans( myAddImage add, QString fn );

#ifndef NOHDF5
    static widImage *readImageHDF5( myAddImage add, QString fn, bool onlyOneImg, bool swapvert );
    static QString hdfFileName;
    static int iterate( QString file_name, QString grp_path, const hid_t loc_id );
    static QStringList slHdfPathAll, slHdfPathImages;
#endif

    static widImage *readImageEDF( myAddImage add, QString fn );
    //static inline int genIDX( int pixCnt, int NbRows, int NbCols )
    //{   // Helper routine to calculate Index
    //    /*  data is written line by line */
    //    int rv = ((NbRows-(pixCnt / NbCols)-1) + NbCols*(pixCnt % NbCols));
    //    if ( rv < 0 || rv >= NbRows*NbCols ) return 0;
    //    return rv;
    //}

    static widImage *readImageSpreadsheet( myAddImage add, QString fn );

};

#endif // SC_READDATA_H
