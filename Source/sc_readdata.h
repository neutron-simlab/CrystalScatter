#ifndef SC_READDATA_H
#define SC_READDATA_H

#include "widimage.h"
#ifndef NOHDF5
#include "H5Cpp.h"
using namespace H5;
#endif

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
