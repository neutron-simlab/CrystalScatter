#include "sc_readdata.h"
#include <QDebug>
#include <QFile>
//#include <QDomDocument>
#include <QXmlStreamReader>

#ifndef NOHDF5
#include "dlghdfselection.h"

#include "H5Cpp.h"
using namespace H5;

QString SC_ReadData::hdfFileName;
QStringList SC_ReadData::slHdfPathAll, SC_ReadData::slHdfPathImages;
#endif


widImage *SC_ReadData::readImageTiff( myAddImage add, QString fn )
{
    /* (vor 16.06.2021) crystal3d1-n.pas line 28523 */
    /* (18.08.2021) crystal3d1.pas line 7740 */
    /* According to http://www.rayonix.com/support/data.htm, the file structure is:
       |-- 1024 bytes TIFF HEADER -------------|
       |-- 3072 bytes frame_header structure --|
       |-- nfast*nslow*depth bytes image ------|
       (Website not available in Mar. 2020)
    */

    QFile fimg(fn);
    if ( ! fimg.open(QIODevice::ReadOnly) )
    {
        qDebug() << fimg.fileName() << fimg.errorString();
        return nullptr;
    }

    // Read the header informations (126 informations)
    // <?xml version="1.0" ?><SAXSLABparameters><param name="det_pixel_size">172e-6 172e-6</param> ...
    // ... </param></SAXSLABparameters>\x00\x00
    // Chars 518 to 5214 in line 15 (Notepad++)
    //
    union
    {
        char     c[20000];
        int16_t  w[10];
        uint16_t uw[10];
        int32_t  i[10];
        uint32_t ui[10];
        float    f[10];
        double   d[10];
    } buf;
    // Part 1: Finden des richtigen Headers. Dieser beginnt mit <?xml version="1.0" ?>
    //         Mal ab Pos 709 oder 643 oder ???
    long long RecRead = fimg.read( buf.c, 20000 );
    // Da auch viele Null-Bytes vorkommen, kann nicht mit String-Suche gearbeitet werden.
    int anf = 0;
    for ( int i=5; i<RecRead; i++ )
    {
        if ( strncmp( &buf.c[i], "<?xml version=", 14 ) == 0 )
        {
            anf = i;
            break;
        }
    }
    qDebug() << "Anfang" << anf;
    // Jetzt ab dem Anfang erneut lesen und die XML-Daten extrahieren
    fimg.seek(anf);
    RecRead = fimg.read( buf.c, 20000 );
    if ( RecRead != 20000 )
    {
        qDebug() << fn << "read" << RecRead;
        fimg.close();
        return nullptr;
    }
/*
    Debug: PARAM ("1.0", "det_pixel_size", "det_thickness", "det_exposure_time", "det_exposure_period", "det_tau", "det_count_cutoff", "det_threshold_setting", "det_n_excluded_pixels", "det_excluded_pixels", "det_flat_field", "det_trim_directory", "datatype", "detectortype", "detector_function", "detector_sn", "meastype", "start_timestamp", "end_timestamp", "save_timestamp", "realtime", "livetime", "pixelsize", "beamcenter_nominal", "beamcenter_actual", "WAXSdet_conf", "data_mean", "data_min", "data_max", "data_rms", "data_p10", "data_p90", "calibrationtype", "kcal", "pixelcal", "koffset", "wavelength", "detector_dist", "saxsconf_r1", "saxsconf_r2", "saxsconf_r3", "saxsconf_l1", "saxsconf_l2", "saxsconf_l3", "detector_dist", "saxsconf_wavelength", "saxsconf_dwavelength", "saxsconf_Imon", "saxsconf_Ieff", "saxsconf_Izero", "saxsconf_det_offx", "saxsconf_det_offy", "saxsconf_det_rotx", "saxsconf_det_roty", "saxsconf_det_pixsizez", "saxsconf_det_pixsizey", "saxsconf_det_resx_0", "saxsconf_det_resy_0", "saxsconf_abs_int_fact", "sample_transfact", "sample_thickness", "sample_ypos", "sample_zpos", "sample_angle1", "sample_angle2", "sample_angle3", "sample_temp", "sample_pressure", "sample_strain", "sample_stress", "sample_shear_rate", "sample_concentration", "sample_buffer", "sample_ID", "hg1", "hp1", "vg1", "vp1", "hg2", "hp2", "vg2", "vp2", "xgon", "zgon", "ygon", "thsam", "ysam", "xsam", "chisam", "phisam", "detx", "dety", "detz", "bstop", "pd", "source_type", "source_runningtime", "source_kV", "source_ma", "xaxis", "xaxisfull", "yaxis", "error_norm_fact", "xaxisbintype", "log", "reduction_type", "reduction_state", "raw_filename", "bsmask_configuration", "mask_filename", "flatfield_filename", "empty_filename", "solvent_filename", "darkcurrent_filename", "readoutnoise_filename", "zinger_removal", "data_added_constant", "data_multiplied_constant", "Img.Class", "Img.MonitorMethod", "Img.ImgType", "Img.Site", "Img.Group", "Img.Researcher", "Img.Operator", "Img.Administrator", "Meas.Description") 127
    Debug:   NUM (" ?", "172e-6 172e-6", "0.000320", "7180.0", "7251.8", "383.8e-09", "1093518", "4024", "16", "badpix_mask.tif", "(nil)", "p300k0140_E8048_T4024_vrf_m0p15.bin", "tiff", "Pilatus", "saxs", "dec383", "", "Mon Mar 02 20:09:38 2020", "", "", "", "7200.00", "0.172 0.172", "372.00    199.00", "371.94    198.97", "", "", "", "", "", "", "", "geom", "", "", "", "1.3414", "536.7868", "0.5000", "0.1500", "0.5000", "725", "400", "200", "536.7868", "1.3414", "0.004", "73171053", "1.14100", "83488172", "0", "0", "0", "0", "0.172", "0.172", "", "", "", "0.06668", "0.20000", "4.700", "1.000", "-71.200", "", "", "", "", "", "", "", "", "", "", "1.000000", "0.072485", "1.000000", "-0.000000", "0.300000", "-0.021621", "0.300000", "-0.008054", "-48.000000", "64.000000", "77.500000", "-71.200007", "4.700000", "-9.949000", "0.000000", "0.850000", "350.000000", "1.041875", "4.926437", "-17.880000", "10.000000", "GENIX3D", "", "0.70", "5.10", "", "", "", "1", "lin", "log", "s", "", "", "", "", "", "", "", "", "", "0", "0", "1", "", "", "2D", "TUM", "", "", "", "", "FeOx_LF7, in conf 12, for 7200 seconds with  0.07 transmission at a Temperature of 22.27C") 127
*/
    QStringList param, num;

    param.clear();
    num.clear();
    bool pstart=false;
    bool nstart=false;
    for ( int ii=1; ii<20000; ii++ )
    {   // This is not the best way to do but it works
        if ( pstart && (buf.c[ii] != '"') ) param.last() += buf.c[ii];
        if ( nstart && (buf.c[ii] != '>') && (buf.c[ii] != '<') ) num.last() += buf.c[ii];
        if ( nstart && (buf.c[ii] == '<') ) nstart=false;
        if ( buf.c[ii] == '"' )
        {
            if ( pstart == false )
            {
                pstart=true;
                param << "";
                num << "";
                //qDebug() << "PARAM: new element" << ii;
            }
            else
            {
                pstart=false;
                nstart=true;
                //qDebug() << "PARAM: fertig" << ii << param.last() << num.last();
                if ( param.last().at(0) == 0x00 )
                {
                    qDebug() << "No parameters found in file.";
                    param.clear();
                    num.clear();
                    break;
                }
            }
        }
    }
    if ( param.size() > 0 )
    {
        param.takeFirst();
        num.takeFirst();
    }
    //qDebug() << "PARAM" << param << param.size();
    //qDebug() << "  NUM" << num << num.size();
    /*28568*/
    bool LittleEndian;
    // Read the first three bytes of the file and check if it is the correct TIF(F) format
    fimg.seek(0);
    fimg.read( buf.c, 8 );  // Bytes 0,1,2=Tiffcode; 4,5,6,7=Offset
    if ( (buf.c[0] == 'I' ) && (buf.c[1] == 'I') )
        LittleEndian = true;
    else if ( (buf.c[0] == 'M') && (buf.c[1] == 'M') )
        LittleEndian = false;
    else
    {
        qDebug() << fn << "ERR: not a tiff file";
        fimg.close();
        return nullptr;
    }
    if ( buf.c[2] != 42 /* '*' */ )
    {
        qDebug() << fn << "ERR: not a tiff file";
        fimg.close();
        return nullptr;
    }
    /*28596*/
    int OffsetIFD = buf.i[1];  // Bytes 4 to 7
    qDebug() << "OffsetIFD" << OffsetIFD;
    fimg.seek( OffsetIFD );
    /*28602*/
    fimg.read( buf.c, 2 );  // BlockRead(MarTifFile,buff2B,1,RecRead);
    int NbTags = buf.w[0];
    //qDebug() << "Offset" << OffsetIFD << "Tags" << NbTags;
    bool isMarTif = false;
    bool isRGBtiff = false;

    int NbCols=0, NbRows=0, BitsPerSample=0, StripOffsets=0;

    for ( int ii=0; ii<NbTags; ii++ )
    {
        fimg.read( buf.c, 12 );  // 2Bytes=TagCode + 2Bytes=TagDataType + 4Bytes=NumValues + 4Bytes=TagData/Offset
        //BlockRead(MarTifFile,buff2B,1,RecRead); //read TagCode
        unsigned int TagCode = buf.uw[0];
        //BlockRead(MarTifFile,buff2B,1,RecRead); //read TagDataType
        int TagDataType = buf.w[1];     // ==> nur für Debug im default:
        //  (* TagDataTypes: 1: byte; 2: ASCII; 3: short; 4: long; 5: rational;
        //     TIFF 6.0 data types: 6: sbyte; 7: undefined; 8: sshort; 9: slong;
        //                     10: srational; 11: float; 12: double
        //     Note: sbyte/sshort/slong = signed byte/short/long *)
        //BlockRead(MarTifFile,buff4B,2,RecRead); //read Nb of Values
        int NbVal = buf.i[1];
        //BlockRead(MarTifFile,buff4B,2,RecRead); //read TagData or offset to it
        /*28643*/
        switch ( TagCode )
        {
        /*
        case 254:Debug: ---TagCode: 254 TagDataType: 4 NbVal: 1 "buf[2]: 0=0x0"
        case 272:Debug: ---TagCode: 272 TagDataType: 2 NbVal: 94 "buf[2]: 996=0x3e4"
        case 278:Debug: ---TagCode: 278 TagDataType: 4 NbVal: 1 "buf[2]: 200=0xc8"
        case 282:Debug: ---TagCode: 282 TagDataType: 5 NbVal: 1 "buf[2]: 1090=0x442"
        case 283:Debug: ---TagCode: 283 TagDataType: 5 NbVal: 1 "buf[2]: 1098=0x44a"
        case 284:Debug: ---TagCode: 284 TagDataType: 3 NbVal: 1 "buf[2]: 1=0x1"
        case 306:Debug: ---TagCode: 306 TagDataType: 2 NbVal: 62 "buf[2]: 1144=0x478"
        case 315:Debug: ---TagCode: 315 TagDataType: 2 NbVal: 20 "buf[2]: 1206=0x4b6"
        case 339:Debug: ---TagCode: 339 TagDataType: 3 NbVal: 1 "buf[2]: 2=0x2"
        case 36865:Debug: ---TagCode: 36865 TagDataType: 4 NbVal: 1 "buf[2]: 0=0x0"
        case 36866:Debug: ---TagCode: 36866 TagDataType: 4 NbVal: 1 "buf[2]: 0=0x0"
        case 36867:Debug: ---TagCode: 36867 TagDataType: 11 NbVal: 1 "buf[2]: 1097859072=0x41700000"
        case 36868:Debug: ---TagCode: 36868 TagDataType: 11 NbVal: 1 "buf[2]: 0=0x0"
        case 36876:Debug: ---TagCode: 36876 TagDataType: 11 NbVal: 1 "buf[2]: 0=0x0"
        case 36877:Debug: ---TagCode: 36877 TagDataType: 11 NbVal: 1 "buf[2]: 0=0x0"
        case 36878:Debug: ---TagCode: 36878 TagDataType: 11 NbVal: 1 "buf[2]: 0=0x0"
        case 36879:Debug: ---TagCode: 36879 TagDataType: 11 NbVal: 1 "buf[2]: 0=0x0"
        case 36880:Debug: ---TagCode: 36880 TagDataType: 11 NbVal: 1 "buf[2]: 0=0x0"
        case 37120:Debug: ---TagCode: 37120 TagDataType: 4 NbVal: 20 "buf[2]: 1226=0x4ca"
            break; // Nicht ausgewertet
        */

        case 256: //ImageWidth (long)
            if ( NbVal == 1 )
            {
                NbCols = buf.i[2];
                qDebug() << "NbCols =" << NbCols;
                param << "NbCols";
                num   << QString::number(NbCols);
            }
            break;
        case 257: //ImageLength (long)
            if ( NbVal == 1 )
            {
                NbRows = buf.i[2];
                qDebug() << "NbRows =" << NbRows;
                param << "NbRows";
                num   << QString::number(NbRows);
            }
            break;
        case 258: //BitsPerSample (short; uses only 2 Bytes)
            if ( NbVal == 1 )
            {
                BitsPerSample = buf.i[2];
                qDebug() << "BitsPerSample =" << BitsPerSample;
                param << "BitsPerSample";
                num   << QString::number(BitsPerSample);
            }
            break;
        case 259: //Compression (short; uses only 2 Bytes)
            if ( NbVal == 1 )
            {
                if ( buf.i[2] != 1 )
                {
                    qDebug() << "ERROR: compressed data. Not supported";
                    fimg.close();
#ifdef WIN32
                    fflush(stderr);
#endif
                    return nullptr;
                }
                qDebug() << "No compression used.";
            }
            break;
        case 262: //PhotometricInterpretation (short; uses only 2 Bytes)
            if ( NbVal == 1 )
            {
                if ( buf.i[2] != 0 && buf.i[2] != 1 )
                {
                    qDebug() << "ERROR: PhotometricInterpretation=" << buf.i[2] << " Not supported";
                    fimg.close();
#ifdef WIN32
                    fflush(stderr);
#endif
                    return nullptr;
                    // isRGBtiff = true;
                }
                qDebug() << "PhotomInterp =" << buf.i[2];
                param << "PhotomInterp";
                num   << QString::number(buf.i[2]);
            }
            break;
        case 269: //DocumentName (ASCII). NbVal can be > 4
            qDebug() << "Document name ignored."; // TODO
            break;
        case 270: //ImageDescription (ASCII). NbVal can be > 4
            qDebug() << "Image description ignored."; // TODO
            break; //ignore ImageDescription
        case 273: //StripOffsets (long)
            if ( NbVal == 1 )
                StripOffsets = buf.i[2];
            else
            {
                qint64 RecPos = fimg.pos();
                fimg.seek( buf.i[2] );
                fimg.read( buf.c, 4 );
                StripOffsets = buf.i[0];
                fimg.seek( RecPos );
            }
            qDebug() << "StripOffsets=" << StripOffsets;
            param << "StripOffsets";
            num   << QString::number(StripOffsets);
            break;
        case 277: //SamplesPerPixel (short; uses only 2 Bytes)
            if ( NbVal == 1 )
            {
                if ( buf.i[2] != 1 )
                {
                    qDebug() << "ERROR: SamplesPerPixel=" << buf.i[2] << " Not supported";
                    fimg.close();
#ifdef WIN32
                    fflush(stderr);
#endif
                    return nullptr;
                }
                qDebug() << "SamplesPerPixel =" << buf.i[2];
                param << "SamplesPerPixel";
                num   << QString::number(buf.i[2]);
            }
            break;
        case 279: //StripByteCounts (short or long). NbVal can be > 1
            qDebug() << "Strip Byte Count ignored."; // TODO
            break;
        case 305: //Software (ASCII). NbVal can be > 4
            qDebug() << "Software ignored."; // TODO
            break;
        case 34710: //Mar specific tag
            if ( NbVal == 1 )
            {
                if ( buf.i[2] == 1024 )
                    isMarTif = true;
                qDebug() << "Mar TagData =" << buf.i[2] << isMarTif;
            }
            break;
        default:
            qDebug() << "---TagCode:" << TagCode << "TagDataType:" << TagDataType << "NbVal:" << NbVal
                     << QString("buf[2]: %1=0x%2").arg(buf.i[2]).arg(buf.i[2],0,16);
            break;
        } // switch
        /*28731*/
    } // for ii<NbTags

    if ( isRGBtiff )
    {
        fimg.close();
        // Tried to check/uncheck CheckBoxSeries before/after Tiff16bitClick, but
        // the actions after the call are not executed anymore.
        //CheckBoxSeries.Checked := True; //to avoid another OpenDialog
        //CheckBoxViewSeries.Checked := True; //to display the image and update EditYmax & EditYmin
        //DataInfoValid := False; //values in the fields (if any) are false
        qDebug() << "*** Tiff16bitClick(Sender); ***"; //LoadTiffFromFile(fileName,totoBmp)
        //CheckBoxSeries.Checked := False;
        //CheckBoxViewSeries.Checked := False;
        return nullptr;
    }
    /*28748*/

    //Reset(MarTifFile,4); // go back to the beginning. RecSize = 4 Bytes.
    //(* The Mar frame_header starts at Byte 1024 *)
    fimg.seek(256*4);
    fimg.read( buf.c, 4 );
    int HeadType = buf.i[0];
    buf.c[4] = 0;
    qDebug() << "HeadType:" << HeadType << "MarTifFlag:" << isMarTif << "LittleEndian:" << LittleEndian << buf.c;
    if ( isMarTif && (HeadType == 2) )
    {
        /*28756*/
        // TODO: diesen Bereich lasse ich weg, da das
        //       angegebene Beispiel nicht diesen Header-Typ hat
        //Seek(MarTifFile,257); //we are already here
        // ...
        /*28903*/
    }
    // Use only NbCols, NbRows, BitsPerSample, StripOffsets to read the image
    int BytesPerPix = BitsPerSample / 8; //BitsPerSample is a multiple of 8
    //Reset(MarTifFile,BytesPerPix); // open MarTifFile for read. RecSize = BytesPerPix
    // Mar and Pilatus tiff images start at StripOffsets=4096
    fimg.seek( StripOffsets );
    int pixCnt, ii=0, jj;
    //bool bkgrd = false; // Flag if background or original image. TODO: get ist from GUI? Then uncomment all code lines with it

    double *origIma = new double[NbRows*NbCols];
//    origIma,zoomIma,origIma1,origIma2: Array[0..maxdim,0..maxdim] of real; //was: double
    double val;

    qDebug() << "READ" << BytesPerPix << LittleEndian << "Start at" << StripOffsets << NbRows*NbCols;

    switch ( BytesPerPix ) // Demo-Picture = 4
    {
    case 1:
        pixCnt = 0;
        while ( ! fimg.atEnd() && pixCnt < NbRows*NbCols )
        {
            //BlockRead(MarTifFile,dataBuff,Trunc(2048/BytesPerPix),RecRead);
            RecRead = fimg.read( buf.c, 2048 );
            for ( int nn=0; nn<RecRead; nn++ )
            {
                val = buf.c[nn];
                ii = pixCnt / NbCols; // data is written line by line
                jj = pixCnt % NbCols;
                //if ( bkgrd )
                //    origIma1[ii+1,jj+1] = val; //or: NbRows-ii if y-axis is inverted
                //else
                    origIma[ii+NbCols*jj] = val;
                pixCnt++;
            } // for nn
        } //while
        break; // 1
    case 2:
        if (LittleEndian)
        {
            pixCnt = 0;
            while ( ! fimg.atEnd() && pixCnt < NbRows*NbCols )
            {
                //BlockRead(MarTifFile,dataBuff,Trunc(2048/BytesPerPix),RecRead);
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<RecRead/2; nn++ )
                {
                    val = buf.w[nn];
                    ii = pixCnt / NbCols; // data is written line by line
                    jj = pixCnt % NbCols;
                    //if ( bkgrd )
                    //    origIma1[ii+1,jj+1] = val; //or: NbRows-ii if y-axis is inverted
                    //else
                        origIma[ii+NbCols*jj] = val;
                    pixCnt++;
                } // for nn
            } //while
        }
        else
        {   //i.e. Big Endian
            pixCnt = 0;
            while ( ! fimg.atEnd() && pixCnt < NbRows*NbCols )
            {
                //BlockRead(MarTifFile,dataBuff,Trunc(2048/BytesPerPix),RecRead);
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<RecRead/2; nn++ )
                {
                    swapBytes( &buf.c[2*nn], 2 );
                    val = buf.w[nn];
                    ii = pixCnt / NbCols; // data is written line by line
                    jj = pixCnt % NbCols;
                    //if ( bkgrd )
                    //    origIma1[ii+1,jj+1] = val; //or: NbRows-ii if y-axis is inverted
                    //else
                        origIma[ii+NbCols*jj] = val;
                    pixCnt++;
                } // for nn
            } //while
        } //if
        break; //Case BytesPerPix=2
    case 4:
/*
        if (LittleEndian=True) then begin
          pixCnt := 0;
          while ((not Eof(MarTifFile)) and (pixCnt < NbRows*NbCols)) do begin
            BlockRead(MarTifFile,dataBuff,Trunc(2048/BytesPerPix),RecRead);
            for nn:= 0 to RecRead-1 do begin
              Ptr2LongWord := Addr(dataBuff[4*nn]);
              ii := pixCnt div NbCols; // data is written line by line
              jj := pixCnt mod NbCols;
               if bkgrd then origIma1[ii+1,jj+1] := Ptr2LongWord^ //or: NbRows-ii if y-axis is inverted
               else origIma[ii+1,jj+1] := Ptr2LongWord^;
              if (Ptr2LongWord^ > maxVal) then
                maxVal := Ptr2LongWord^;
              if ((Ptr2LongWord^ > 0) and (Ptr2LongWord^ < minVal))
                then minVal := Ptr2LongWord^;
              Inc(pixCnt);
            end;
            ProgressBar1.Position:=ii;
          end; //while
        end
*/
        if (LittleEndian)
        {
            int nn=0;
            pixCnt = 0;
            while ( ! fimg.atEnd() && pixCnt < NbRows*NbCols )
            {
                //BlockRead(MarTifFile,dataBuff,Trunc(2048/BytesPerPix),RecRead);
                RecRead = fimg.read( buf.c, 2048 );
                for ( nn=0; nn<RecRead/4 && pixCnt < NbRows*NbCols; nn++ )
                {
                    val = buf.i[nn];
                    ii = pixCnt / NbCols; // data is written line by line
                    jj = pixCnt % NbCols;
                    //if ( bkgrd )
                    //    origIma1[ii+1,jj+1] = val; //or: NbRows-ii if y-axis is inverted
                    //else
                    //origIma[pixCnt] = val;
                    origIma[jj+NbCols*(NbRows-1-ii)] = val;
                    //if ( (ii < 10 && jj < 10) || (ii >= NbRows-3 && jj >= NbCols-3) )
                    //    qDebug() << ii << jj << jj+NbCols*(NbRows-1-ii);
                    pixCnt++;
                } // for nn
            } //while
            qDebug() << "Pixelcount:" << pixCnt;
        }
        else
        {   //i.e. Big Endian
            pixCnt = 0;
            while ( ! fimg.atEnd() && pixCnt < NbRows*NbCols )
            {
                //BlockRead(MarTifFile,dataBuff,Trunc(2048/BytesPerPix),RecRead);
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<RecRead/4; nn++ )
                {
                    swapBytes( &buf.c[4*nn], 4 );
                    val = buf.i[nn];
                    ii = pixCnt / NbCols; // data is written line by line
                    jj = pixCnt % NbCols;
                    //if ( bkgrd )
                    //    origIma1[ii+1,jj+1] = val; //or: NbRows-ii if y-axis is inverted
                    //else
                    origIma[ii+NbCols*jj] = val;
                    pixCnt++;
                } // for nn
            } //while
        } //if
        break; //Case BytesPerPix=4
    default:
        break;
    }
    fimg.close();
/*29025*/

    widImage* img = add( 0, NbCols, 0, NbRows, origIma, "TIF-Image" );
    img->addMetaInfo( "From File", fn );    // TIFF
    for ( int i=0; i<param.size(); i++ )
        img->addMetaInfo( param.at(i), num.at(i) );
    img->addMetaInfo("@",""); // Sortieren
    return img;
}

void SC_ReadData::swapBytes( char *c, unsigned int len )
{
    char buf[8];
    for ( unsigned int k=0; k<len; k++ )
        buf[k] = c[len-k-1];
    memcpy( buf, c, len );
}





widImage *SC_ReadData::readImageSasCrystal( myAddImage add, QString fnpar )
{
    int pos = fnpar.lastIndexOf(".");

    QString fncsv = fnpar.left(pos)+".csv"; // Einlesen für die Daten
    /* (1) the first two lines contains the dimensions.
       ts << "X:;" << minX << ";" << maxX << EOL;
       ts << "Y:;" << minY << ";" << maxY << EOL;
       (2) then all values line by line separated with a ; so this can be read by Excel.
       for ( int y=minY; y<maxY; y++ )
       { ts << "Y=" << y << ":";
         for ( int x=minX; x<maxX; x++ ) ts << ";" << getData(x,y);
         ts << EOL;
       }
    */

    QString fndat = fnpar.left(pos)+".dat"; // könnten auch 2D-SANS(ILL) Daten sein
    /* union { double e; char c[1]; } val;
       (1) 19 Bytes as text with the current dimensions of the image.
       QString tmp = QString("%1 %2 %3 %4").arg(minX,4).arg(maxX,4).arg(minY,4).arg(maxY,4);
       f.write( qPrintable(tmp) );
       (2) then all values in binary format.
       for ( int y=minY; y<maxY; y++ )
         for ( int x=minX; x<maxX; x++ )
         { val.e = getData(x,y);
           f.write( val.c, sizeof(double) );
         }
    */

    QString fntxt = fnpar.left(pos)+".txt"; // EInlesen für die Meta-Informationen
    /* for ( int r=0; r<ui->tblMetaData->rowCount(); r++ )
         ts << ui->tblMetaData->item(r,0)->text() << "\t"
            << ui->tblMetaData->item(r,1)->text() << EOL;
    */

    QFile fimg(fncsv);
    if ( ! fimg.open(QIODevice::ReadOnly) )
    {
        qDebug() << fimg.fileName() << fimg.errorString();
        return nullptr;
    }

    QVector<double> daten;
    QStringList sl;
    int minX=0, maxX=0, minY=0, maxY=0, NbRows=0, NbCols=0;
    while ( !fimg.atEnd() )
    {
        sl = QString(fimg.readLine()).trimmed().split(";",Qt::SkipEmptyParts);
        if ( sl.size() < 2 ) continue;
        if ( sl[0] == "X:" )
        {   // ts << "X:;" << minX << ";" << maxX << EOL;
            minX = sl[1].toInt();
            maxX = sl[2].toInt();
            NbCols = maxX - minX;
            //qDebug() << "cols" << minX << maxX << NbCols;
        }
        else if ( sl[0] == "Y:" )
        {   // ts << "Y:;" << minY << ";" << maxY << EOL;
            minY = sl[1].toInt();
            maxY = sl[2].toInt();
            NbRows = maxY - minY;
            //qDebug() << "rows" << minY << maxY << NbRows;
        }
        else if ( sl[0].startsWith("Y=") && sl.size() > NbCols )
        {   // ts << "Y=" << y << ":";
            // for ( int x=minX; x<maxX; x++ ) ts << ";" << getData(x,y);
            if ( daten.size() == 0 )
                daten.reserve( NbRows * NbCols );
            //qDebug() << sl.size() << sl;
            for ( int i=1; i<sl.size(); i++ )
                daten.append( sl[i].toDouble() );
        }
    }
    fimg.close();
    //qDebug() << daten;
    //qDebug() << daten.size();

    // Das Image muss noch um 90° gedreht werden, dann entspricht es dem Ausduck
    // Die Y-Achse muss gespiegelt werden, da die Grafik mit Y=0 oben arbeitet
    // Dann wird auch das Simulationsbild mit den gleichen BeamStop-Koordinaten
    // passend dazu berechnet.
    double *imgData = new double[daten.size()];
    for ( int r=0, d=0; r<NbRows; r++ )     // y
        for ( int c=0; c<NbCols; c++, d++ ) // x
            imgData[XY2IDX(0,NbCols,0,NbRows,c,r)] = daten[d];
    //#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))

    QFile ftxt(fntxt);
    if ( ! ftxt.open(QIODevice::ReadOnly) )
    {
        qDebug() << ftxt.fileName() << ftxt.errorString();
        return nullptr;
    }
    widImage* img = add( 0, NbCols, 0, NbRows, imgData, "SasCrystal-Image" );
    img->addMetaInfo( "From File", fncsv ); // SasCrystalImage
    img->addMetaInfo( "NbRows", QString::number(NbRows) );
    img->addMetaInfo( "NbCols", QString::number(NbCols) );
    while ( ! ftxt.atEnd() )
    {
        sl = QString(ftxt.readLine()).trimmed().split("\t",Qt::SkipEmptyParts);
        if ( sl.size() < 2 ) continue;
        if ( sl[0] == "From File" ) continue;
        if ( sl[0] == "NbRows" ) continue;
        if ( sl[0] == "NbCols" ) continue;
        img->addMetaInfo( sl[0], sl[1] );
    }
    ftxt.close();
    return img;
}






widImage *SC_ReadData::readImageKWSData( myAddImage add, QString fn )
{
    QFile fimg(fn);
    if ( ! fimg.open(QIODevice::ReadOnly) )
    {
        qDebug() << fimg.fileName() << fimg.errorString();
        return nullptr;
    }

    // Sondertest: die Endung ".dat" kann auch von 2D-SANS (ILL) Daten stammen. Dann kommt in den ersten 30 Zeilen
    //             auf jeden Fall der Text "Data columns" vor. Wenn das der Fall ist, muss die Routine
    //             readImage2dSans() stattdessen aufgerufen werden.
    for ( int i=0; i<40 && !fimg.atEnd(); i++ )
    {
        if ( QString(fimg.readLine()).trimmed().contains("Data columns") )
        {
            qDebug() << "KWS -> 2D-SANS";
            fimg.close();
            return readImage2dSans(add,fn);
        }
    }
    qDebug() << "Real KWS";
    fimg.seek(0);
    // Sondertest ende

    // Fileformat (Leerzeilen entfernt):
    //****************************************************************
    // KWS1_MEASUREMENT KFA-IFF spectrum_data file version V-5.00
    //>> Jetzt stehen die Detektordaten als Integerwerte im File (Rohdaten) sonst sind sie in Float
    // wood2_3152_Stan_C4_S2_D0.DAT
    // Standard_Sample measurement started by Mr(s). kws1 at 21-Nov-2014 17:42:18.00
    // (* Statistics *)
    // Measurement produced 0 Informations 0 Warnings 0 Errors 0 Fatals
    // Cyclus_Number Reduce_Data Date_field from to
    //         3152          NO               0  0
    // (* Comment *)
    // kws1 | n.szekely@fz-juelich.de
    // wood decomposition IL 100 oC | lambda=4.73A
    //>>                              Wellenlänge
    // wood IL old - 20m
    // (* Collimation discription *)
    // Coll_Position Wind(1)_Pos Beamwindow_X Beamwindow_Y
    //          [m]         [m]         [mm]         [mm]
    //            4           4           30           30
    //            4           4           30           30
    //>>                                 Beamsize
    // (* Detector Discription *)
    // Li6 Detector is in normal mode. Angle: 0.00 grd
    // Offset Z_Position X_Position Y_Position
    //   [m]        [m]       [mm]       [mm]
    //  -0.38       1.50      -4.10     9.00
    //  -0.38       1.50      -4.10     9.00
    //>>           DetPos
    // (* Sample discription *)
    // Sample_Nr Sample_Pos Thicknes Beamwindow_X Beamwindow_Y Time_Factor
    //                [mm]     [mm]         [mm]         [mm]           *
    //        2     224.10         2.00            8            8        1.00
    //        2     318.00        0            16           17
    // (* Temperature discription *)
    // Temperature dummy line
    // (* Data_field and Time per Data_Field *)
    // Data_Field     Time Time_Factor  Repetition
    //         0    5 min        1.00           1
    // (* Selector and Monitor Counter *)
    // Selector Monitor_1 Monitor_2 Monitor_3
    //  141481  17186386   2216543       565
    // (* Measurement stop state *)
    // measurement STOPPED by USER command
    // (* Calculate measurement time = Time per Data_field*Time_Factor*Repetition *)
    //   5 min
    // (* Real measurement time for detector data *)
    //     300 sec
    // (* Detector Data Sum (Sum_Data_Field = 32 values = 6 * 5 and 2) *)
    // @
    // 2.876773E+07 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
    // 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
    // 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
    // 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
    // 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
    // 0.000000E+00 0.000000E+00
    // (* Detector Data (Data_Field = 16384 values = 8 * 2048) *)
    // $
    //>> Start of data
    //       0        0        0        0        0        0        0        0
    //
    //****************************************************************
    // Anderes Fileformat:
    //****************************************************************
    // $
    //>> Start of data
    // 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00

    bool    isHdr=true,
            lookForLambda=false,
            //lookForColli=false,
            lookForDetpos=false, detHdr=true;
    QVector<double> daten;
    QString line;
    QString wavelen = "0.5 nm";
    QString detpos = "4.0 m";
    while ( !fimg.atEnd() )
    {
        line = QString(fimg.readLine()).trimmed();
        if ( isHdr )
        {
            if ( line.startsWith("(* Comment *)") )
                lookForLambda = true;
            else if ( lookForLambda && line.contains("lambda=") )
            {
                lookForLambda = false;
                int p = line.indexOf("lambda=") + 7;
                line = line.mid(p);
                double wl = line.left(line.length()-1).trimmed().toDouble();
                if ( line.endsWith("A") )
                    wl = wl / 10;
                wavelen = QString("%1 nm").arg(wl);
                //qDebug() << "WAVELEN:" << line << wavelen;
            }
            /*if ( line.startsWith("(* Collimation discription *)") )
                lookForColli = true;
            else if ( lookForColli )
            {
                if ( line.startsWith("(*") )
                    lookForColli = false;
                else
                    qDebug() << "COLLI:" << line;
            }*/
            if ( line.startsWith("(* Detector Discription *)") )
                lookForDetpos = true;
            else if ( lookForDetpos )
            {
                if ( line.startsWith("(*") )
                    lookForDetpos = false;
                else if ( detHdr && line[0] == '[' )
                    detHdr = false;
                else if ( !detHdr )
                {
                    QStringList sl = line.split(" ",Qt::SkipEmptyParts);
                    if ( sl.size() == 4 )
                        detpos = sl[1] + " m";
                    //qDebug() << "DETPOS:" << line << detpos;
                }
            }
            if ( line.startsWith("$") )
                isHdr = false;
            continue;
        } // if isHdr
        QStringList sl = line.split(" ",Qt::SkipEmptyParts);
        if ( sl.size() < 8 ) continue; // evtl. Leerzeilen am Ende
        for ( int i=0; i<sl.size(); i++ )
        {
            double val = sl.at(i).toDouble();
            if ( val < 0 ) { /*qDebug() << "KWS data" << val;*/ val = 0; }
            daten.append(val);
        }
    }
    fimg.close();
    //qDebug() << daten;
    //qDebug() << daten.size();

    int NbRows, NbCols, BeamPosX, BeamPosY;
    if ( daten.size() == 128*128 )
    {
        NbRows = NbCols = 128;
        BeamPosX = BeamPosY = 64;
    }
    else if ( daten.size() == 144*144 )
    {
        NbRows = NbCols = 144;
        BeamPosX = BeamPosY = 72;
    }
    else
    {
        qDebug() << daten.size() << "undefined format";
        return nullptr;
    }

    // Das Image muss noch um 90° gedreht werden, dann entspricht es dem Ausduck
    // Die Y-Achse muss gespiegelt werden, da die Grafik mit Y=0 oben arbeitet
    // Dann wird auch das Simulationsbild mit den gleichen BeamStop-Koordinaten
    // passend dazu berechnet.
    double *imgData = new double[daten.size()];
    for ( int r=0, d=0; r<NbRows; r++ )     // y
        for ( int c=0; c<NbCols; c++, d++ ) // x
            imgData[XY2IDX(0,NbCols,0,NbRows,c,r)] = daten[d];
    //#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))
    widImage* img = add( 0, NbCols, 0, NbRows, imgData, "KWS-Image" );
    img->addMetaInfo( "From File", fn );    // KWS-Data
    img->addMetaInfo( "NbRows", QString::number(NbRows) );
    img->addMetaInfo( "NbCols", QString::number(NbCols) );
    img->addMetaInfo( "BeamPosX", QString::number(BeamPosX) );
    img->addMetaInfo( "BeamPosY", QString::number(BeamPosY) );
    img->addMetaInfo( "wavelength", wavelen );
    img->addMetaInfo( "EditWAXSangle", "0.0 deg" );
    img->addMetaInfo( "SampleDist", detpos );
    img->addMetaInfo( "Pixel_X", QString::number(5.3) );
    img->addMetaInfo( "Pixel_Y", QString::number(5.3) );
    // bei den wenigen nicht nötig img->addMetaInfo("@",""); // Sortieren

    // Jetzt noch das automatische Suchen nach dem BeamCenter. Dazu muss ein Wert vorgegeben sein
    // (oben als MetaInfo) und das Ergebnis wird dort überschrieben.
    //findBeamCenter( img, BeamPosX, BeamPosY );
    //img->addMetaInfo( "BeamPosX", QString::number(BeamPosX) );
    //img->addMetaInfo( "BeamPosY", QString::number(BeamPosY) );
    return img;
}


void SC_ReadData::findBeamCenter( widImage *ikws, int &xerg, int &yerg )
{
    double ref = ikws->dataPtr()[0];
    if ( ref > 0.001 )
    {
        //qDebug() << "findBeamCenter" << ref;
        ref = 0;
    }
    int x0 = ikws->getFileInfos()->centerX;
    int y0 = ikws->getFileInfos()->centerY;
    int bsxmin=x0, bsxmax=x0, bsymin=y0, bsymax=y0;
    if ( ikws->getData(x0,y0) > ref )
    {   // Echte Suche starten ...
        bsxmin = ikws->xmin();
        bsxmax = ikws->xmax()-1;
        bsymin = ikws->ymin();
        bsymax = ikws->ymax()-1;
        int maxbsx = 2;
        int maxbsy = 2;
        for ( int y=ikws->ymin(); y<ikws->ymax(); y++ )
        {
            // Für jede Zeile:
            // 1. Suche die Grenzen von außen
            int xx0=ikws->xmin(), xx1=ikws->xmax()-1;
            while ( xx0 < xx1 && ikws->getData(xx0,y) <= ref ) xx0++;
            while ( xx0 < xx1 && ikws->getData(xx1,y) <= ref ) xx1--;
            // 2. Suche innerhalb des Unterbereiches den größten Bereich 'außen' (> 3 px)
            bool inner = false;
            int tmpx;
            for ( int x=xx0; x<=xx1; x++ )
            {
                if ( ikws->getData(x,y) <= ref )
                {   // BS?
                    if ( inner )
                    {   // War BS
                        continue;
                    }
                    else
                    {   // War Daten
                        tmpx = x;
                        inner = true;
                    }
                }
                else
                {   // Ist Daten
                    if ( inner )
                    {   // War BS
                        if ( x - tmpx > maxbsx )
                        {
                            bsxmin = tmpx;
                            bsxmax = x-1;
                            maxbsx = x - tmpx;
                        }
                        inner = false;
                    }
                    else
                    {   // War Daten
                        continue;
                    }
                }
            }
        }

        for ( int x=ikws->xmin(); x<ikws->xmax(); x++ )
        {
            // Für jede Spalte:
            // 1. Suche die Grenzen von außen
            int yy0=ikws->ymin(), yy1=ikws->ymax()-1;
            while ( yy0 < yy1 && ikws->getData(x,yy0) <= ref ) yy0++;
            while ( yy0 < yy1 && ikws->getData(x,yy1) <= ref ) yy1--;
            // 2. Suche innerhalb des Unterbereiches den größten Bereich 'außen' (> 3 px)
            bool inner = false;
            int tmpy;
            for ( int y=yy0; y<=yy1; y++ )
            {
                if ( ikws->getData(x,y) <= ref )
                {   // BS?
                    if ( inner )
                    {   // War BS
                        continue;
                    }
                    else
                    {   // War Daten
                        tmpy = y;
                        inner = true;
                    }
                }
                else
                {   // Ist Daten
                    if ( inner )
                    {   // War BS
                        if ( y - tmpy > maxbsy )
                        {
                            bsymin = tmpy;
                            bsymax = y-1;
                            maxbsy = y - tmpy;
                        }
                        inner = false;
                    }
                    else
                    {   // War Daten
                        continue;
                    }
                }
            }
        }
        xerg = (bsxmin + bsxmax) / 2;
        yerg = (bsymin + bsymax) / 2;
    }
    else
    {
        while ( bsxmin >= ikws->xmin() && ikws->getData(bsxmin,y0) <= ref ) bsxmin--;
        while ( bsxmax <  ikws->xmax() && ikws->getData(bsxmax,y0) <= ref ) bsxmax++;
        while ( bsymin >= ikws->ymin() && ikws->getData(x0,bsymin) <= ref ) bsymin--;
        while ( bsymax <  ikws->ymax() && ikws->getData(x0,bsymax) <= ref ) bsymax++;
        xerg = (bsxmin + bsxmax) / 2;
        yerg = (bsymin + bsymax) / 2;
    }
}



#ifndef NOHDF5

typedef struct
{
    char *name;
    int type;
} name_type_t;

// Callback functions for HDF5 routines
#ifdef __cplusplus
extern "C" {
#endif

static herr_t count_objects_cb( hid_t, const char *, const H5L_info_t *, void *_op_data )
{
    hsize_t *op_data = (hsize_t *)_op_data;
    (*op_data)++;
    return H5_ITER_CONT;
}

static herr_t get_name_type_cb( hid_t loc_id, const char *name, const H5L_info_t *, void *op_data )
{
    H5G_stat_t stat_buf;
    H5Gget_objinfo(loc_id, name, 0, &stat_buf); // < 0)
    ((name_type_t *)op_data)->type = stat_buf.type;
    ((name_type_t *)op_data)->name = (char *)strdup(name);
    // define H5_ITER_STOP for return. This will cause the iterator to stop
    return H5_ITER_STOP;
}

#ifdef __cplusplus
}
#endif


widImage *SC_ReadData::readImageHDF5( myAddImage add, QString fn, bool onlyOneImg, bool swapvert )
{
    QList<widImage*> images;    // Wenn mehrere Datensätze in einem File sind, dann werden
                                //  auch mehrere Images angezeigt. Das kann später noch
                                //  geändert werden (TODO).

    H5File file( qPrintable(fn), H5F_ACC_RDONLY );
    //if ( file.isValid(...) )

    hid_t fid = file.getId();
    //if ( (fid = H5Fopen( qPrintable(fn), H5F_ACC_RDONLY, H5P_DEFAULT ) ) < 0 )
    if ( fid < 0 )
    {
        qDebug() << "Unable to open" << fn;
        hdfFileName = "";
        return nullptr;
    }
    hdfFileName = fn;

    /*
    //---------- Globale Informationen ausgeben (nur zum Test)
    H5F_info_t finfo;
    if ( H5Fget_info( fid, &finfo ) >= 0 )
    {
        qDebug() << "H5F-INFO:" << fn;
        qDebug() << "  SUPER:" << finfo.super.version << finfo.super.super_size << finfo.super.super_ext_size;
        qDebug() << "   FREE:" << finfo.free.version << finfo.free.tot_space << finfo.free.meta_size;
        qDebug() << "   SOHM:" << finfo.sohm.version << finfo.sohm.hdr_size
                 << finfo.sohm.msgs_info.heap_size<< finfo.sohm.msgs_info.index_size;
    }
    */

    //---------- Iterieren über alle Elemente im File
    slHdfPathAll.clear();
    slHdfPathImages.clear();
    iterate( fn, "/", fid );

    qDebug() << "ALL:" << slHdfPathAll;
    //qDebug() << "IMG:" << slHdfPathImages;

    if ( slHdfPathImages.size() > 0 )
    {   // Datensätze zu lesen macht nur Sinn, wenn auch welche da sind.

        dlgHdfSelection *dlg = nullptr;
        if ( slHdfPathImages.size() > 1 )
        {   // Nur bei mehr als einem Datensatz den User fragen...
            dlg = new dlgHdfSelection( fn, slHdfPathImages );
            dlg->exec();
        }

        foreach ( QString dsname, slHdfPathImages )
        {
            if ( dlg != nullptr &&  ! dlg->isChecked(dsname) ) continue;
            QString dbg = "DS=" + dsname;
            DataSet dataset = file.openDataSet( qPrintable(dsname) );

            DataSpace dataspace = dataset.getSpace();
            int rank = dataspace.getSimpleExtentNdims();
            hsize_t dims_out[rank];
            int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
            dbg += QString(", rank=%1, ndims=%2, ").arg(rank).arg(ndims);
            size_t dimsize = 1;
            for ( int i=0; i<rank; i++ )
            {
                dimsize *= dims_out[i];
                dbg += QString(" [%1]=%2").arg(i).arg(dims_out[i]);
            }
            qDebug() << dbg;    // Nur eine Zeile pro Image

            H5T_class_t type_class = dataset.getTypeClass();
            if ( type_class == H5T_INTEGER )
            {
                qDebug() << "**INT**" << dimsize;  // ist im Beispiel vorhanden

                int buffer[dimsize];
                dataset.read( buffer, PredType::NATIVE_INT ); // , memspace, dataspace );

                double dblbuf[dimsize];
                if ( swapvert )
                {   // Vertikal umdrehen
                    for ( hsize_t r=0; r<dims_out[0]; r++ )     // y
                        for ( hsize_t c=0; c<dims_out[1]; c++ ) // x
                            dblbuf[XY2IDX(0,dims_out[1],0,dims_out[0],c,r)] = buffer[XY2IDX(0,dims_out[1],0,dims_out[0],c,dims_out[0]-r-1)];
                    //#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))

                }
                else
                {   // Daten nur Formatmäßig kopieren
                    int    *pi = buffer;
                    double *pd = dblbuf;
                    for ( size_t i=0; i<dimsize; i++ ) *(pd++) = *(pi++);
                }
                widImage* img = add( 0, dims_out[1], 0, dims_out[0], dblbuf, "HDF5-Image" );
                img->addMetaInfo( "Dataset", dsname );
                img->addMetaInfo( "From File", fn );    // HDF5 (Int)
                img->addMetaInfo( "NbRows", QString::number(dims_out[0]) );
                img->addMetaInfo( "NbCols", QString::number(dims_out[1]) );
                //img->addMetaInfo( "BeamPosX", QString::number(BeamPosX) );
                //img->addMetaInfo( "BeamPosY", QString::number(BeamPosY) );
                //img->addMetaInfo( "wavelength", wavelen );
                //img->addMetaInfo( "EditWAXSangle", "0.0 deg" );
                //img->addMetaInfo( "SampleDist", detpos );
                //img->addMetaInfo( "Pixel_X", QString::number(5.3) );
                //img->addMetaInfo( "Pixel_Y", QString::number(5.3) );
                images.append(img);
            }
            else if ( type_class == H5T_FLOAT )
            {
                qDebug() << "**FLOAT**" << dimsize;

                float buffer[dimsize];
                dataset.read( buffer, PredType::NATIVE_FLOAT ); // , memspace, dataspace );

                double dblbuf[dimsize];
                if ( swapvert )
                {   // Vertikal umdrehen
                    for ( hsize_t r=0; r<dims_out[0]; r++ )     // y
                        for ( hsize_t c=0; c<dims_out[1]; c++ ) // x
                            dblbuf[XY2IDX(0,dims_out[1],0,dims_out[0],c,r)] = buffer[XY2IDX(0,dims_out[1],0,dims_out[0],c,dims_out[0]-r-1)];
                    //#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))

                }
                else
                {   // Daten nur Formatmäßig kopieren
                    float  *pi = buffer;
                    double *pd = dblbuf;
                    for ( size_t i=0; i<dimsize; i++ ) *(pd++) = *(pi++);
                }
                widImage* img = add( 0, dims_out[1], 0, dims_out[0], dblbuf, "HDF5-Image" );
                img->addMetaInfo( "Dataset", dsname );
                img->addMetaInfo( "From File", fn );    // HDF5 (Float)
                img->addMetaInfo( "NbRows", QString::number(dims_out[0]) );
                img->addMetaInfo( "NbCols", QString::number(dims_out[1]) );
                //img->addMetaInfo( "BeamPosX", QString::number(BeamPosX) );
                //img->addMetaInfo( "BeamPosY", QString::number(BeamPosY) );
                //img->addMetaInfo( "wavelength", wavelen );
                //img->addMetaInfo( "EditWAXSangle", "0.0 deg" );
                //img->addMetaInfo( "SampleDist", detpos );
                //img->addMetaInfo( "Pixel_X", QString::number(5.3) );
                //img->addMetaInfo( "Pixel_Y", QString::number(5.3) );
                images.append(img);
            }
            else
            {
                qDebug() << "Invalid data type";
            }
            dataset.close();
            if ( onlyOneImg ) break; // Automatic fit sollte nur ein Image lesen
        } // foreach ( QString dsname, slHdfPathImages )

    }

    H5Fclose(fid);

    if ( images.size() == 0 ) return nullptr;
    return images.first();
}

int SC_ReadData::iterate( QString file_name, QString grp_path, const hid_t loc_id )
{
    hsize_t nbr_objects = 0;
    hsize_t index;
    name_type_t info;
    QString path;
    bool do_iterate;

    H5Literate( loc_id, H5_INDEX_NAME, H5_ITER_INC, NULL, count_objects_cb, &nbr_objects );
    //qDebug() << "ITERATE:" << grp_path << nbr_objects;

    for ( hsize_t idx_obj=0; idx_obj < nbr_objects; idx_obj++ )
    {
        index = idx_obj;
        H5Literate( loc_id, H5_INDEX_NAME, H5_ITER_INC, &index, get_name_type_cb, &info );
        // initialize path
        path = grp_path;
        if ( grp_path != "/" )
            path += "/";
        path += info.name;

        switch(info.type)
        {
        ///////////////////////////////////////////////////////////////////////////////////////
        //H5G_GROUP
        //////////////////////////////////////////////////////////////////////////////////////
        case H5G_GROUP:
        {
            hid_t gid;
            H5O_info_t oinfo_buf;

            gid = H5Gopen2(loc_id, info.name, H5P_DEFAULT);
            do_iterate = true;
            //get object info
            H5Oget_info(gid, &oinfo_buf,0);
            //qDebug() << "group" << info.name;

            //iterate in sub group passing QTreeWidgetItem for group as parent
            if ( do_iterate ) iterate( file_name, path, gid );

            // TODO (separate Routine): get_attributes(file_name, path, gid, item_grp);

            H5Gclose(gid);
            free(info.name);
            break;
        }

        ///////////////////////////////////////////////////////////////////////////////////////
        //H5G_DATASET
        //////////////////////////////////////////////////////////////////////////////////////
        case H5G_DATASET:
        {
            hid_t did, sid, ftid, mtid, attr;
            size_t datatype_size;
            H5T_sign_t datatype_sign;
            H5T_class_t datatype_class;
            hsize_t dims[H5S_MAX_RANK];
            herr_t ret;
            int rank;

            did = H5Dopen2(loc_id, info.name, H5P_DEFAULT);
            sid = H5Dget_space(did);
            rank = H5Sget_simple_extent_dims(sid, dims, NULL);
            ftid = H5Dget_type(did);
            mtid = H5Tget_native_type(ftid, H5T_DIR_DEFAULT);

            ///////////////////////////////////////////////////////////////////////////////////////
            //store datatype sizes and metadata needed to display HDF5 buffer data
            ///////////////////////////////////////////////////////////////////////////////////////
            datatype_size = H5Tget_size(mtid); // ) == 0)
            datatype_sign = H5Tget_sign(mtid); // ) < 0)
            datatype_class = H5Tget_class(mtid); // ) < 0)
            H5Sclose(sid);
            H5Tclose(ftid);
            H5Tclose(mtid);
            static QString dtc2str[] = { "Integer", "Float", "Time", "String", "Bitfield",
                                         "Opaque", "Compound", "Reference", "Enum", "VLen",
                                         "Array" };
            if ( datatype_class == H5T_NO_CLASS ||
                 datatype_class >= static_cast<int>(sizeof(dtc2str)/sizeof(dtc2str[0])) )
                qDebug() << "dataset" << path << datatype_size << datatype_sign << "*DataTypeClass ERROR*";
            else
            {
                QString tmp = "Dims:";
                bool allDimsOne = true;
                for ( int i=0; i<rank && i<H5S_MAX_RANK; i++ )
                {
                    tmp += QString(" %1").arg(dims[i]);
                    if ( dims[i] > 1 ) allDimsOne = false;
                }
                QString sign;
                switch ( datatype_sign )
                {
                case H5T_SGN_NONE:  sign = "Unsigned "; break;
                case H5T_SGN_2:     sign = "Signed ";   break;
                case H5T_SGN_ERROR:
                default:            sign = "";          break;
                }
                if ( rank < 2 || allDimsOne )
                {
                    switch ( datatype_class )
                    {
                    default:
                    case H5T_NO_CLASS:   //= -1, /**< error                                   */
                        break;
                    case H5T_INTEGER:    //= 0,  /**< integer types                           */
                    {
                        //int val;
                        qDebug() << "dataset" << grp_path << info.name << datatype_size << "Bytes" << sign + dtc2str[datatype_class] << tmp;
                        //attr = H5Aopen_by_name(did, qPrintable(grp_path), info.name, H5P_DEFAULT, H5P_DEFAULT);
                        //ret = H5Aread(attr, H5T_NATIVE_INT, &val);
                        //ret = H5Aclose(attr);
                        //qDebug() << "        =" << val;
                        break;
                    }
                    case H5T_FLOAT:      //= 1,  /**< floating-point types                    */

                        break;
                    case H5T_TIME:       //= 2,  /**< date and time types                     */
                        break;
                    case H5T_STRING:     //= 3,  /**< character string types                  */

                        break;
                    case H5T_BITFIELD:   //= 4,  /**< bit field types                         */
                        break;
                    case H5T_OPAQUE:     //= 5,  /**< opaque types                            */
                        break;
                    case H5T_COMPOUND:   //= 6,  /**< compound types                          */
                        break;
                    case H5T_REFERENCE:  //= 7,  /**< reference types                         */
                        break;
                    case H5T_ENUM:       //= 8,  /**< enumeration types                       */
                        break;
                    case H5T_VLEN:       //= 9,  /**< variable-Length types                   */
                        break;
                    case H5T_ARRAY:      //= 10, /**< array types                             */
                        break;
                    }
                }
                //qDebug() << "dataset" << path << datatype_size << "Bytes" << sign + dtc2str[datatype_class] << tmp;
                slHdfPathAll << path + "  (" + sign + dtc2str[datatype_class] + ", " + tmp + ")";
                if ( rank >= 2 && !allDimsOne /*&& dims[rank-1] >= 32 && dims[rank-2] >= 32*/ )    // TODO: sind Bilder kleiner 32*32 sinnvoll?
                    slHdfPathImages << path; // + "  (" + sign + dtc2str[datatype_class] + ", " + tmp + ")";
            }

            H5Dclose(did);
            break;
        } // case H5G_DATASET:
        } // switch(info.type)
    }
    return 0;
}
#endif // !NOHDF5






widImage *SC_ReadData::readImage2dSans( myAddImage add, QString fn )
{
    // TODO: es fehlen noch Dialoge, um das Format genauer zu spezifizieren

    QFile fimg(fn);
    if ( ! fimg.open(QIODevice::ReadOnly) )
    {
        qDebug() << fimg.fileName() << fimg.errorString();
        return nullptr;
    }

    int NbRows=0, NbCols=0;
    QVector<double> daten;
    QHash<QString,QString> metadata;

    if ( fn.endsWith(".xml",Qt::CaseInsensitive) )
    {
        bool firstrow = true;
        QStringList path;
        QString tmp;
        QXmlStreamReader xml(&fimg);
        while ( ! xml.atEnd() )
        {
            xml.readNext();
            switch ( xml.tokenType() )
            {
            case QXmlStreamReader::NoToken:
                // The reader has not yet read anything.
                break;
            case QXmlStreamReader::Invalid:
                // An error has occurred, reported in error() and errorString().
                break;
            case QXmlStreamReader::StartDocument:
                // The reader reports the XML version number in documentVersion(),
                // and the encoding as specified in the XML document in documentEncoding().
                // If the document is declared standalone, isStandaloneDocument() returns true;
                // otherwise it returns false.
                //qDebug() << "START doc" << xml.documentVersion() << xml.documentEncoding();
                path.clear();
                break;
            case QXmlStreamReader::EndDocument:
                // The reader reports the end of the document.
                qDebug() << "END doc" << path;
                break;
            case QXmlStreamReader::StartElement:
                // The reader reports the start of an element with namespaceUri() and name().
                // Empty elements are also reported as StartElement, followed directly by
                // EndElement. The convenience function readElementText() can be called to
                // concatenate all content until the corresponding EndElement. Attributes
                // are reported in attributes(), namespace declarations in namespaceDeclarations().
                if ( xml.name() == "SASdata" )
                {
                    NbRows++;
                    if ( NbRows > 1 ) firstrow = false;
                }
                else if ( xml.name() == "Idata" && firstrow )
                {
                    NbCols++;
                }
                else if ( xml.name() == "I" )
                {
                    QString val = xml.readElementText();
                    daten.append( val.toDouble() );
                }
                else if ( xml.name() == "Q" ||
                          xml.name() == "Idev" ||
                          xml.name() == "Qdev" )
                {
                }
                else if ( xml.name() == "ID" ||
                          xml.name() == "name" ||
                          xml.name() == "radiation" ||
                          xml.name() == "date" )
                {
                    tmp = xml.readElementText();
                    qDebug() << path << xml.name() << tmp;
                    if ( ! tmp.isEmpty() )
                        metadata.insert( path.join("/")+"/"+xml.name(), tmp );
                }
                else if ( xml.name() != "Idata" &&
                          xml.name() != "SASroot" &&
                          xml.name() != "SASentry" )
                {
                    path.append( xml.name().toString() );
                    //qDebug() << "  Start" << xml.name();
                }

                break;
            case QXmlStreamReader::EndElement:
                // The reader reports the end of an element with namespaceUri() and name().
                tmp = xml.name().toString();
                if ( path.contains(tmp) )
                {
                    while ( path.contains(tmp) ) path.takeLast();
                    //qDebug() << "  End" << xml.name() << path;
                }
                break;
            case QXmlStreamReader::Characters:
                // The reader reports characters in text(). If the characters are all white-space,
                // isWhitespace() returns true. If the characters stem from a CDATA section,
                // isCDATA() returns true.
            case QXmlStreamReader::Comment:
                // The reader reports a comment in text().
            case QXmlStreamReader::DTD:
                // The reader reports a DTD in text(), notation declarations in notationDeclarations(),
                // and entity declarations in entityDeclarations(). Details of the DTD declaration
                // are reported in in dtdName(), dtdPublicId(), and dtdSystemId().
            case QXmlStreamReader::EntityReference:
                // The reader reports an entity reference that could not be resolved. The name
                // of the reference is reported in name(), the replacement text in text().
            case QXmlStreamReader::ProcessingInstruction:
                // The reader reports a processing instruction in processingInstructionTarget()
                // and processingInstructionData().
                break;
            }
        }
        if ( xml.hasError() )
        {
            qDebug() << "Read Error:" << xml.errorString()
                     << "Line" << xml.lineNumber(); // , columnNumber(), and characterOffset() ;
            fimg.close();
            return nullptr;
        }
        fimg.close();
        //qDebug() << NbRows << NbCols;
        //return nullptr;
    } // xml

    else if ( fn.endsWith(".csv",Qt::CaseInsensitive) )
    {
        // A       -0.03375,       -0.03325,       -0.03275, ...
        // 0            nan,            nan,            nan, ...
        // 1            nan,            nan,            nan, ...
        // ...
        // 135            nan,            nan,            nan, ...
        //
        // ERRORS
        // 0            nan,            nan,            nan, ...
        // 1            nan,            nan,            nan, ...
        // ...
        // 135            nan,            nan,            nan, ...

        QString line;
        bool ok;
        while ( !fimg.atEnd() )
        {
            line = QString(fimg.readLine()).trimmed();
            if ( line.length() < 2 ) break; // nur der erste Block wird verwendet
            if ( line.startsWith("A") ) continue; // erste Zeile ignorieren
            QStringList sl = line.mid(3).split(",",Qt::SkipEmptyParts);
            if ( sl.size() < 8 ) continue; // evtl. Leerzeilen am Ende
            int i = line.left(3).trimmed().toInt(&ok);
            if ( ok ) NbRows = i+1;
            if ( NbCols == 0 ) NbCols = sl.size();
            for ( i=0; i<sl.size(); i++ )
            {
                if ( sl.at(i).startsWith("nan") )
                    daten.append(0.0);
                else
                {
                    double val = sl.at(i).trimmed().toDouble();
                    //if ( val < 0 ) { /*qDebug() << "2D SANS data" << val;*/ val = 0; }
                    daten.append(val);
                }
            }
        }
        //qDebug() << daten;
        qDebug() << daten.size();

        fimg.close();
    } // csv

    else if ( fn.endsWith(".dat",Qt::CaseInsensitive) )
    {
        // FORMAT -A-
        //  Data columns Qx - Qy - I(Qx,Qy) - err(I)
        //  ASCII data
        //  -0.149  -0.149  0.77225  0.218756
        //  ...
        //
        // FORMAT -B-
        //  FILE: SILIC010.SA3_SRK_S110   CREATED: 22-JAN-2008 02:45:55
        //  LABEL: Silica 2 pct in D2O Scatt 4M
        //  MON CNT   LAMBDA (A)  DET_OFF(cm)   DET_DIST(m)   TRANS   THICK(cm)
        //  4.4712e+06  6       0     4     0.72452     0.2
        //  BCENT(X,Y)   A1(mm)   A2(mm)   A1A2DIST(m)   DL/L   BSTOP(mm)   DET_TYP
        //  68.15  64.79  50    3    8.5505    0.142    50.8    ORNL
        //  SAM: RAW Data File: SILIC010.SA3_SRK_S110
        //  BGD: none
        //  EMP: none
        //  DIV: none
        //  MASK: none
        //  ABS Parameters (3-6): none
        //  Average Choices: none
        //
        //  *** Data written from RAW folder and may not be a fully corrected data file ***
        //  Data columns are Qx - Qy - I(Qx,Qy) - Qz - SigmaQx - SigmaQy - fSubS(beam stop shadow)
        //
        //  ASCII data created Tue, Mar 09, 2010 3:40:01 PM
        //
        //  -0.08731621	-0.08294715	7	0.006948369	0.005297193	0.00507565	1
        //  ...

        QStringList cols, sl;
        QString line;
        bool ok, hdr=true;
        int index=-1;
        while ( !fimg.atEnd() )
        {
            line = QString(fimg.readLine()).trimmed();
            if ( line.length() < 2 ) continue; // Leerzeilen werden ignoriert
            if ( hdr )
            {   // Header
                qDebug() << line;
                if ( line.startsWith("Data columns") )
                {
                    qDebug() << "COLS:" << line;
                    line.replace("Data columns ","");
                    line.replace("are ","");
                    cols = line.split(" - ");
                    for ( int i=0; i<cols.size(); i++ )
                        if ( cols[i].startsWith("I") )
                        {
                            index = i;
                            break;
                        }
                    qDebug() << "     " << cols << index;
                }
                else if ( line.startsWith("ASCII") )
                {
                    hdr = false;
                    if ( index < 0 )
                    {
                        fimg.close();
                        return nullptr;
                    }
                }
                else if ( line.contains(":") )
                {
                    int pos = line.indexOf(":"); // es können mehrere : vorkommen (Uhrzeit)
                    metadata.insert( line.left(pos), line.mid(pos+1).trimmed() );
                }
                // Andere Headerdaten werden erstmal ignoriert (TODO)
            }
            else
            {   // Daten
                sl = line.split( QRegExp("[\t ]"), Qt::SkipEmptyParts );
                //qDebug() << line << sl;
                //break;
                double val = sl[index].toDouble(&ok);
                if ( ok )
                    daten.append( val );
                else
                    daten.append( 0.0 );
            }
        }
        NbRows = sqrt( daten.size() );
        if ( NbRows < 2 ) NbRows = 2;
        NbCols = daten.size() / NbRows;
        qDebug() << daten.size() << NbRows << NbCols;
        fimg.close();
    }
    else
    {
        fimg.close();
        return nullptr;
    } // txt

    // Bild vertikal invertieren, nicht nur visuell
    double *imgData = new double[daten.size()];
    for ( int r=0; r<NbRows; r++ )     // y
        for ( int c=0; c<NbCols; c++ ) // x
            imgData[XY2IDX(0,NbCols,0,NbRows,c,r)] = daten[XY2IDX(0,NbCols,0,NbRows,c,NbRows-r-1)];
    //#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))
    widImage* img = add( 0, NbCols, 0, NbRows, imgData, "2D SANS-Image" );
//    widImage* img = add( 0, NbCols, 0, NbRows, daten.data(), "2D SANS-Image" );
    img->addMetaInfo( "From File", fn );    // 2d SANS
    img->addMetaInfo( "NbRows", QString::number(NbRows) );
    img->addMetaInfo( "NbCols", QString::number(NbCols) );

    QHash<QString, QString>::const_iterator i = metadata.constBegin();
    while (i != metadata.constEnd())
    {
        img->addMetaInfo( i.key(), i.value() );
        ++i;
    }
    //img->addMetaInfo("@",""); // Sortieren

    return img;
}




widImage *SC_ReadData::readImageEDF( myAddImage add, QString fn )
{
    //{NV}X/5003 procedure TFormCrystal3d.EdfDataClick(Sender: TObject);
    QFile fimg(fn);
    if ( ! fimg.open(QIODevice::ReadOnly) )
    {
        qDebug() << fimg.fileName() << fimg.errorString();
        return nullptr;
    }
    QString line;
    while ( !fimg.atEnd() )
    {
        line = fimg.readLine().trimmed();
        if ( line.startsWith("{") ) break;
    }

    QString fileType;
    QStringList keywords;
    int EDF_BinarySize, headSize, NbCols, NbRows;
    QString ByteOrder, DataType, Title;
    double ExposureTime, BeamPosX, BeamPosY, pixel_x, pixel_y, SampleDist, Wavelength, DetectorPos;

    line = fimg.readLine().trimmed();
    if ( line.startsWith("EDF_DataBlockID") )
    {
        // EDF_DataFormatVersion (Klora, several DataBlocks e.g. image & variance)
        // EDF_DataBlockID (Klora, 1 DataBlock per file)
        // HeaderID (edf, 1 DataBlock per file (backward compatibility)

        fileType = "K1"; // Klora, 1 DataBlock per file
        keywords.clear();
        keywords << "EDF_BinarySize" << "EDF_HeaderSize" << "ByteOrder" << "DataType" << "Dim_1"
                 << "Dim_2" << "Title" << "ExposureTime" << "Center_1" << "Center_2" << "PSize_1"
                 << "PSize_2" << "SampleDistance" << "WaveLength" << "DetectorPosition";
        // In principle the keywords appear always in this order,
        // but it should work anyway thank to the StringToCaseSelect function.
        while ( ! fimg.atEnd() )
        {
            line = fimg.readLine().trimmed();
            if ( line.contains("}") ) break;
            QStringList tmp = line.split(" = ");
            if ( tmp.size() < 2 ) continue;
            if ( tmp[1].endsWith(";") ) tmp[1] = tmp[1].left(tmp[1].length()-1).trimmed();
            //NOTE: standard Case doesn't work with Strings
            switch ( keywords.indexOf(tmp[0]) )
            {
            case  0 : EDF_BinarySize = tmp[1].toInt(); break;
            case  1 : headSize = tmp[1].toInt(); break;
            case  2 : ByteOrder = tmp[1]; break;
            case  3 : DataType = tmp[1]; break;
            case  4 : NbCols = tmp[1].toInt(); break;
            case  5 : NbRows = tmp[1].toInt(); break;
            case  6 : Title = tmp[1]; break;
            case  7 : ExposureTime = tmp[1].toDouble(); break;
            case  8 : BeamPosX = tmp[1].toDouble(); break; //[pixels]. Was: StrToInt
            case  9 : BeamPosY = tmp[1].toDouble(); break; //[pixels]. Was: StrToInt
            case 10 : pixel_x  = 1e3 * tmp[1].toDouble(); break; //[mm]
            case 11 : pixel_y  = 1e3 * tmp[1].toDouble(); break; //[mm]
            case 12 : SampleDist = tmp[1].toDouble(); break; //[m]
            case 13 : Wavelength = 1e9 * tmp[1].toDouble(); break; //[nm]
            case 14 : DetectorPos = tmp[1].toDouble(); break;
            } // switch
        }
    } // if ( line.startsWith("EDF_DataBlockID") )

    if ( line.startsWith("EDF_DataFormatVersion") )
    {
        qDebug() << "ATTN: more than 1 DataBlock per file. FileType 'K2' is not supported."; // TODO
        fimg.close();
        return nullptr;
    }

    if ( line.startsWith("HeaderID") )
    {
        fileType = "E1"; // edf, 1 DataBlock per file
        keywords.clear();
        keywords << "ByteOrder" << "DataType" << "Dim_1" << "Dim_2" << "Size";
        //NOTE: one can also find  in the header: count_time, col_bin, row_bin.
        //      As they are not indispensable, we can ignore them.
        while ( ! fimg.atEnd() )
        {
            line = fimg.readLine().trimmed();
            if ( line.contains("}") ) break;
            QStringList tmp = line.split(" = ");
            if ( tmp.size() < 2 ) continue;
            if ( tmp[1].endsWith(";") ) tmp[1] = tmp[1].left(tmp[1].length()-1).trimmed();
            //NOTE: standard Case doesn't work with Strings
            switch ( keywords.indexOf(tmp[0]) )
            {
            case 0 : ByteOrder = tmp[1]; break;
            case 1 : DataType = tmp[1]; break;
            case 2 : NbCols = tmp[1].toInt(); break;
            case 3 : NbRows = tmp[1].toInt(); break;
            case 4 : EDF_BinarySize = tmp[1].toInt(); break;
            case 5 : if ( ! tmp[1].contains("na",Qt::CaseInsensitive) ) ExposureTime = tmp[1].toDouble();
                     break;
            } // switch
        }
        BeamPosX   = round(NbCols/2.);
        BeamPosY   = round(NbRows/2.);
        pixel_x    = 0.260;      //[mm]
        pixel_y    = 0.260;      //[mm]
        SampleDist = 5.005;   //[m] Range: [1.5 to 8]
        Wavelength = 1.24e-1; //[nm]
    } //if ( line.startsWith("HeaderID") )

    int RecSize = EDF_BinarySize / NbCols / NbRows; //record size (in Bytes)

    // Theoretisch sollte das File schon richtig stehen (Nach der Zeile mit "}" ...
    // Seek(edfBinFile, Trunc(headSize/RecSize));

    // ATTN: y-axis is inverted. Data is written line by line
    double maxVal = 0.0;
    double minVal = 1e20;
    int pixCnt;

    union
    {
        char     c[20000];
        int16_t  w[10];
        uint16_t uw[10];
        int32_t  i[10];
        uint32_t ui[10];
        float    f[10];
        double   d[10];
    } buf;
    long long RecRead;

    double *origIma = new double[NbRows*NbCols];
    double val;

    //    /*  data is written line by line */
#define genIDX(pc,nr,nc) ((nr-(pc/nc)-1) + nc*(pc%nc))

    //qDebug() << RecSize << DataType << ByteOrder << NbRows << NbCols;
    //Debug: 4 "FloatValue" "LowByteFirst" 1920 1920

    if ( (RecSize == 4) && (DataType=="FloatValue" || DataType=="FloatIEEE32") )
    {
        if ( ByteOrder == "LowByteFirst" )
        {
            // Testdatensatz hier!
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<RecRead/4; nn++ )
                {
                    val = buf.f[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        }
        else
        { /* i.e. if ByteOrder = HighByteFirst */
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<=RecRead/4; nn++ )
                {
                    swapBytes( &buf.c[4*nn], 4 );
                    val = buf.f[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        } /* else if (ByteOrder) */
    } /* if (RecSize = 4) and AnsiMatchStr(DataType,['FloatValue'...]) */

    else if ( (RecSize == 4) && (DataType=="SignedInteger" || DataType=="Signed32" || DataType=="SignedLong") )
    {
        if ( ByteOrder == "LowByteFirst" )
        {
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<=RecRead/4; nn++ )
                {
                    val = buf.i[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        }
        else
        { /* i.e. if ByteOrder = HighByteFirst */
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<=RecRead/4; nn++ )
                {
                    swapBytes( &buf.c[4*nn], 4 );
                    val = buf.i[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        } /* else if (ByteOrder) */
    } /* if (RecSize = 4) and AnsiMatchStr(DataType,['SignedInteger'...]) */

    else if ( (RecSize == 2) && (DataType=="UnsignedShort" || DataType=="Unsigned16") )
    {/*3*/
        if ( ByteOrder == "LowByteFirst" )
        {
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<=RecRead/2; nn++ )
                {
                    val = buf.uw[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        }
        else
        { /* i.e. if ByteOrder = HighByteFirst */
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<=RecRead/2; nn++ )
                {
                    swapBytes( &buf.c[2*nn], 2 );
                    val = buf.uw[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        } /* else if (ByteOrder) */
    } /* if (RecSize = 2) and AnsiMatchStr(DataType,['UnsignedShort'...]) */

    else if ( (RecSize == 8) && (DataType=="DoubleValue" || DataType=="DoubleIEEE64") )
    {/*1*/
        if ( ByteOrder == "LowByteFirst" )
        {
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<=RecRead/8; nn++ )
                {
                    val = buf.d[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        }
        else
        { /* i.e. if ByteOrder = HighByteFirst */
            pixCnt  =  0;
            while ( !fimg.atEnd() )
            {
                RecRead = fimg.read( buf.c, 2048 );
                for ( int nn=0; nn<=RecRead/8; nn++ )
                {
                    swapBytes( &buf.c[8*nn], 8 );
                    val = buf.d[nn];
                    origIma[genIDX(pixCnt,NbRows,NbCols)] = val;
                    if ( val > maxVal )
                        maxVal = val;
                    if ( val > 0 && val < minVal )
                        minVal = val;
                    pixCnt++;
                }
            } /* while */
        } /* else if (ByteOrder) */
    } /* if (RecSize = 8) and AnsiMatchStr(DataType,['DoubleValue'...]) */

    else
    {
        qDebug() << "EdfDataClick: data format not supported";
        fimg.close();
        return nullptr;
    }
    fimg.close();

/*
    // Das Image muss noch um 90° gedreht werden, dann entspricht es dem Ausduck
    // Die Y-Achse muss gespiegelt werden, da die Grafik mit Y=0 oben arbeitet
    // Dann wird auch das Simulationsbild mit den gleichen BeamStop-Koordinaten
    // passend dazu berechnet.
    double *imgData = new double[daten.size()];
    for ( int r=0, d=0; r<NbRows; r++ )     // y
        for ( int c=0; c<NbCols; c++, d++ ) // x
            imgData[XY2IDX(0,NbCols,0,NbRows,c,r)] = daten[d];
    //#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))
*/
    widImage* img = add( 0, NbCols, 0, NbRows, origIma, "EDF-Image" );
    img->addMetaInfo( "From File", fn );    // EDF
    img->addMetaInfo( "NbRows", QString::number(NbRows) );
    img->addMetaInfo( "NbCols", QString::number(NbCols) );
    img->addMetaInfo( "BeamPosX", QString::number(BeamPosX) );
    img->addMetaInfo( "BeamPosY", QString::number(BeamPosY) );
    img->addMetaInfo( "FileType", fileType );
    img->addMetaInfo( "Pixel_x", QString::number(pixel_x) );
    img->addMetaInfo( "Pixel_y", QString::number(pixel_y) );
    img->addMetaInfo( "SampleDist", QString::number(SampleDist) );
    img->addMetaInfo( "wavelength", QString::number(Wavelength) );
    img->addMetaInfo( "EditWAXSangle", "0.0 deg" );
    img->addMetaInfo( "EDF_BinarySize", QString::number(EDF_BinarySize) );
    img->addMetaInfo( "headSize", QString::number(headSize) );
    img->addMetaInfo( "ByteOrder", ByteOrder );
    img->addMetaInfo( "DataType", DataType );
    img->addMetaInfo( "Title", Title );
    img->addMetaInfo( "ExposureTime", QString::number(ExposureTime) );
    img->addMetaInfo( "DetectorPos", QString::number(DetectorPos) );
    img->addMetaInfo("@",""); // Sortieren

    return img;
}



widImage *SC_ReadData::readImageSpreadsheet( myAddImage add, QString fn )
{
    QFile fimg(fn);
    if ( ! fimg.open(QIODevice::ReadOnly) )
    {
        qDebug() << fimg.fileName() << fimg.errorString();
        return nullptr;
    }

    // Fileformat:
    //
    // 161    161 Start pixel = (       1       1)
    //    3.27986E-01 2.55160E-01 3.56233E-01 3.37742E-01 3.65466E-01 .....
    //    ....

/*
    fo.write( qPrintable(QString("    %1    %2 Start pixel = (       1       1)").arg(sx).arg(sy)+EOL) );
    for ( int y=0; y<sy; y++ )
    {
        QString line = " ";
        for ( int x=0; x<sx; x++ )
            line += QString(" %1").arg(src[x+sx*y],prec+7,'e',prec);    // -1.123456e+01
        fo.write( qPrintable(line+EOL) );
    }
*/

    QVector<double> daten;
    QString line;
    int pos;
    line = QString(fimg.readLine()).trimmed();
    if ( (pos=line.indexOf("Start pixel")) < 0 )
    {
        fimg.close();
        return nullptr;
    }
    QStringList sl = line.left(pos).split(" ",Qt::SkipEmptyParts);
    int NbRows, NbCols, BeamPosX, BeamPosY;
    NbRows = sl[0].toInt();
    NbCols = sl[1].toInt();
    BeamPosX = NbCols/2;  // Da Meta-Infos im File fehlen
    BeamPosY = NbRows/2;

    for ( int z=0; z<NbRows; z++ )
    {
        if ( fimg.atEnd() ) break;

        line = QString(fimg.readLine()).trimmed();
        sl = line.split(" ",Qt::SkipEmptyParts);
        if ( sl.size() < NbCols ) continue; // evtl. Leerzeilen am Ende
        for ( int i=0; i<sl.size(); i++ )
        {
            double val = sl.at(i).toDouble();
            daten.append(val);
        }
    }
    fimg.close();
    //qDebug() << daten;
    //qDebug() << daten.size();

    // Das Image muss noch um 90° gedreht werden, dann entspricht es dem Ausduck
    // Die Y-Achse muss gespiegelt werden, da die Grafik mit Y=0 oben arbeitet
    // Dann wird auch das Simulationsbild mit den gleichen BeamStop-Koordinaten
    // passend dazu berechnet.
    double *imgData = new double[daten.size()];
    for ( int r=0, d=0; r<NbRows; r++ )     // y
        for ( int c=0; c<NbCols; c++, d++ ) // x
            imgData[XY2IDX(0,NbCols,0,NbRows,c,r)] = daten[d];
    //#define XY2IDX(X0,X1,Y0,Y1,x,y) ((-X0 + (x)) + (X1-X0)*(-Y0 + (y)))
    widImage* img = add( 0, NbCols, 0, NbRows, imgData, "SPR-Image" );
    img->addMetaInfo( "From File", fn );
    img->addMetaInfo( "NbRows", QString::number(NbRows) );
    img->addMetaInfo( "NbCols", QString::number(NbCols) );
    img->addMetaInfo( "BeamPosX", QString::number(BeamPosX) );
    img->addMetaInfo( "BeamPosY", QString::number(BeamPosY) );
    // bei den wenigen nicht nötig img->addMetaInfo("@",""); // Sortieren

    // Jetzt noch das automatische Suchen nach dem BeamCenter. Dazu muss ein Wert vorgegeben sein
    // (oben als MetaInfo) und das Ergebnis wird dort überschrieben.
    findBeamCenter( img, BeamPosX, BeamPosY );
    img->addMetaInfo( "BeamPosX", QString::number(BeamPosX) );
    img->addMetaInfo( "BeamPosY", QString::number(BeamPosY) );
    return img;
}


/* Notizen zum Filetype *.NXS (ist HDF5 oder Nexus)
 * File: C:\SimLab\sas-crystal\20220221 - 2D SANS\028114.nxs
 *
 * "/entry0/D11/Beam/center_x  (Float, Dims: 1)"
 * "/entry0/D11/Beam/center_y  (Float, Dims: 1)"
 * "/entry0/D11/Beam/sample_ap_x_or_diam  (Float, Dims: 1)"
 * "/entry0/D11/Beam/sample_ap_y  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/data  (Signed Integer, Dims: 256 192 1)" ----------------
 * "/entry0/D11/Detector 1/description  (String, Dims: 1)"
 * "/entry0/D11/Detector 1/det_actual  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/det_calc  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/det_offset  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/det_requested  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/detrate  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/detsum  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/gaz_pressure  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/pixel_size_x  (Float, Dims: 1)"
 * "/entry0/D11/Detector 1/pixel_size_y  (Float, Dims: 1)"
 * "/entry0/D11/Detector 2/data  (Signed Integer, Dims: 32 256 1)"  ----------------
 * "/entry0/D11/Detector 2/description  (String, Dims: 1)"
 * "/entry0/D11/Detector 2/det_actual  (Float, Dims: 1)"
 * "/entry0/D11/Detector 2/det_calc  (Float, Dims: 1)"
 * "/entry0/D11/Detector 2/detrate  (Float, Dims: 1)"
 * "/entry0/D11/Detector 2/detsum  (Float, Dims: 1)"
 * "/entry0/D11/Detector 2/gaz_pressure  (Float, Dims: 1)"
 * "/entry0/D11/Detector 2/pixel_size_x  (Float, Dims: 1)"
 * "/entry0/D11/Detector 2/pixel_size_y  (Float, Dims: 1)"
 * "/entry0/D11/Detector 3/data  (Signed Integer, Dims: 32 256 1)"  ----------------
 * "/entry0/D11/Detector 3/description  (String, Dims: 1)"
 * "/entry0/D11/Detector 3/det_actual  (Float, Dims: 1)"
 * "/entry0/D11/Detector 3/det_calc  (Float, Dims: 1)"
 * "/entry0/D11/Detector 3/detrate  (Float, Dims: 1)"
 * "/entry0/D11/Detector 3/detsum  (Float, Dims: 1)"
 * "/entry0/D11/Detector 3/gaz_pressure  (Float, Dims: 1)"
 * "/entry0/D11/Detector 3/pixel_size_x  (Float, Dims: 1)"
 * "/entry0/D11/Detector 3/pixel_size_y  (Float, Dims: 1)"
 * "/entry0/D11/Distance/Gap  (Float, Dims: 1)"
 * "/entry0/D11/Distance/det1_actual  (Float, Dims: 1)"
 * "/entry0/D11/Distance/det1_calc  (Float, Dims: 1)"
 * "/entry0/D11/Distance/det2_calc  (Float, Dims: 1)"
 * "/entry0/D11/Distance/det3_calc  (Float, Dims: 1)"
 * "/entry0/D11/Distance/detectorTankWindow  (String, Dims: 1)"
 * "/entry0/D11/IN6Beam/final_wavelength  (Float, Dims: 1)"
 * "/entry0/D11/MeasurementType/SampleProposalComment  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/SampleProposalName  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/SampleProposalid  (Signed Integer, Dims: 1)"
 * "/entry0/D11/MeasurementType/subType  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/subTypeID  (Signed Integer, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeEnvEch  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeEnvEchID  (Signed Integer, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeGroupInst  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeGroupInstID  (Signed Integer, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeOfMeasure  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeOfMeasureID  (Signed Integer, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeOfSample  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/typeOfSampleID  (Signed Integer, Dims: 1)"
 * "/entry0/D11/MeasurementType/typePorteEch  (String, Dims: 1)"
 * "/entry0/D11/MeasurementType/typePorteEchID  (Signed Integer, Dims: 1)"
 * "/entry0/D11/attenuator/attenuation_coefficient  (Float, Dims: 1)"
 * "/entry0/D11/attenuator/attenuation_state  (String, Dims: 1)"
 * "/entry0/D11/attenuator/attenuation_transmission  (Float, Dims: 1)"
 * "/entry0/D11/attenuator/attenuation_value  (String, Dims: 1)"
 * "/entry0/D11/attenuator/attenuator_distance  (String, Dims: 1)"
 * "/entry0/D11/attenuator/attenuator_type  (String, Dims: 1)"
 * "/entry0/D11/attenuator/position  (Signed Integer, Dims: 1)"
 * "/entry0/D11/beamstop/actual_beamstop_number  (Signed Integer, Dims: 1)"
 * "/entry0/D11/beamstop/bx1_actual  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx1_offset  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx1_requested  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx2_actual  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx2_offset  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx2_requested  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx_actual  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx_offset  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/bx_requested  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/by_actual  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/by_offset  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/by_requested  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/height  (Float, Dims: 1)"
 * "/entry0/D11/beamstop/width  (Float, Dims: 1)"
 * "/entry0/D11/collimation/Disk1  (Float, Dims: 1)"
 * "/entry0/D11/collimation/Disk2  (Float, Dims: 1)"
 * "/entry0/D11/collimation/Disk3  (Float, Dims: 1)"
 * "/entry0/D11/collimation/Disk4  (Float, Dims: 1)"
 * "/entry0/D11/collimation/Disk5  (Float, Dims: 1)"
 * "/entry0/D11/collimation/actual_position  (Float, Dims: 1)"
 * "/entry0/D11/collimation/col10_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col11_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col12_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col13_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col1_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col2_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col3_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col4_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col5_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col6_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col7_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col8_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/col9_actual_state  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/diaphragm1_position  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/diaphragm2_position  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/diaphragm3_position  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/diaphragm4_position  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/diaphragm5_position  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/diaphragm6_position  (Signed Integer, Dims: 1)"
 * "/entry0/D11/collimation/guide_exit_cross_section_height  (Float, Dims: 1)"
 * "/entry0/D11/collimation/guide_exit_cross_section_width  (Float, Dims: 1)"
 * "/entry0/D11/collimation/sourceAperture  (String, Dims: 1)"
 * "/entry0/D11/collimation/sourceDistance  (Float, Dims: 1)"
 * "/entry0/D11/name  (String, Dims: 1)"
 * "/entry0/D11/selector/rotation_speed  (Float, Dims: 1)"
 * "/entry0/D11/selector/selector_speed_average  (Float, Dims: 1)"
 * "/entry0/D11/selector/selector_speed_sigma  (Float, Dims: 1)"
 * "/entry0/D11/selector/selrot_actual  (Float, Dims: 1)"
 * "/entry0/D11/selector/selrot_offset  (Float, Dims: 1)"
 * "/entry0/D11/selector/selrot_requested  (Float, Dims: 1)"
 * "/entry0/D11/selector/wavelength  (Float, Dims: 1)"
 * "/entry0/D11/selector/wavelength_res  (Float, Dims: 1)"
 * "/entry0/D11/version  (String, Dims: 1)"
 * "/entry0/acquisition_mode  (Signed Integer, Dims: 1)"
 * "/entry0/data1/MultiDet1_data  (Signed Integer, Dims: 256 192 1)"    ----------------
 * "/entry0/data2/MultiDet2_data  (Signed Integer, Dims: 32 256 1)"     ----------------
 * "/entry0/data3/MultiDet3_data  (Signed Integer, Dims: 32 256 1)"     ----------------
 * "/entry0/duration  (Float, Dims: 1)"
 * "/entry0/end_time  (String, Dims: 1)"
 * "/entry0/experiment_identifier  (String, Dims: 1)"
 * "/entry0/experiment_title  (String, Dims: 1)"
 * "/entry0/instrument_name  (String, Dims: 1)"
 * "/entry0/mode  (Float, Dims: 1)"
 * "/entry0/modestring  (String, Dims: 1)"
 * "/entry0/monitor1/data  (Signed Integer, Dims: 1 1 1)"
 * "/entry0/monitor1/mode  (String, Dims: 1)"
 * "/entry0/monitor1/monrate  (Float, Dims: 1)"
 * "/entry0/monitor1/monsum  (Float, Dims: 1)"
 * "/entry0/monitor1/preset  (Float, Dims: 1)"
 * "/entry0/monitor2/data  (Signed Integer, Dims: 1 1 1)"
 * "/entry0/monitor2/mode  (String, Dims: 1)"
 * "/entry0/monitor2/monrate  (Float, Dims: 1)"
 * "/entry0/monitor2/monsum  (Float, Dims: 1)"
 * "/entry0/monitor2/preset  (Float, Dims: 1)"
 * "/entry0/reactor_power  (Float, Dims: 1)"
 * "/entry0/run_number  (Signed Integer, Dims: 1)"
 * "/entry0/sample/Actual bath  (String, Dims: 1)"
 * "/entry0/sample/Eurotherm818_1_actualTemperature  (Float, Dims: 1)"
 * "/entry0/sample/Eurotherm818_1_wantedTemperature  (Float, Dims: 1)"
 * "/entry0/sample/Eurotherm818_2_actualTemperature  (Float, Dims: 1)"
 * "/entry0/sample/Eurotherm818_2_wantedTemperature  (Float, Dims: 1)"
 * "/entry0/sample/Pressure1  (Float, Dims: 1)"
 * "/entry0/sample/Pressure2  (Float, Dims: 1)"
 * "/entry0/sample/Pressure3  (Float, Dims: 1)"
 * "/entry0/sample/Pressure4  (Float, Dims: 1)"
 * "/entry0/sample/Pressure5  (Float, Dims: 1)"
 * "/entry0/sample/Pressure6  (Float, Dims: 1)"
 * "/entry0/sample/Pressure7  (Float, Dims: 1)"
 * "/entry0/sample/Pressure8  (Float, Dims: 1)"
 * "/entry0/sample/Wavegenerator_frequency_requested  (Float, Dims: 1)"
 * "/entry0/sample/Wavegenerator_voltage_requested  (Float, Dims: 1)"
 * "/entry0/sample/Wavegenerator_waveform_requested  (String, Dims: 1)"
 * "/entry0/sample/air_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath1_regulation_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath1_setpoint_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath1_temperature_requested  (Float, Dims: 1)"
 * "/entry0/sample/bath2_regulation_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath2_setpoint_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath2_temperature_requested  (Float, Dims: 1)"
 * "/entry0/sample/bath_reservoir_regulation_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath_reservoir_requested  (Float, Dims: 1)"
 * "/entry0/sample/bath_reservoir_setpoint_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath_sample_humidity  (Float, Dims: 1)"
 * "/entry0/sample/bath_sample_regulation_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath_sample_setpoint_temperature  (Float, Dims: 1)"
 * "/entry0/sample/bath_sample_temperature_requested  (Float, Dims: 1)"
 * "/entry0/sample/chemicalFormula  (String, Dims: 1)"
 * "/entry0/sample/cold_valve_pressure  (Float, Dims: 1)"
 * "/entry0/sample/consistance  (String, Dims: 1)"
 * "/entry0/sample/density  (String, Dims: 1)"
 * "/entry0/sample/field_actual  (Float, Dims: 1)"
 * "/entry0/sample/field_requested  (Float, Dims: 1)"
 * "/entry0/sample/hp_temperature2  (Float, Dims: 1)"
 * "/entry0/sample/keithley_voltage  (Float, Dims: 1)"
 * "/entry0/sample/mass  (String, Dims: 1)"
 * "/entry0/sample/omega_actual  (Float, Dims: 1)"
 * "/entry0/sample/omega_offset  (Float, Dims: 1)"
 * "/entry0/sample/omega_requested  (Float, Dims: 1)"
 * "/entry0/sample/ph  (Float, Dims: 1)"
 * "/entry0/sample/ph_voltage  (Float, Dims: 1)"
 * "/entry0/sample/phi_actual  (Float, Dims: 1)"
 * "/entry0/sample/phi_offset  (Float, Dims: 1)"
 * "/entry0/sample/phi_requested  (Float, Dims: 1)"
 * "/entry0/sample/ps1_current  (Float, Dims: 1)"
 * "/entry0/sample/ps1_voltage  (Float, Dims: 1)"
 * "/entry0/sample/ps2_current  (Float, Dims: 1)"
 * "/entry0/sample/ps2_voltage  (Float, Dims: 1)"
 * "/entry0/sample/pump1_speed  (Float, Dims: 1)"
 * "/entry0/sample/pump2_speed  (Float, Dims: 1)"
 * "/entry0/sample/pump3_speed  (Float, Dims: 1)"
 * "/entry0/sample/quantum_sample_temp  (Float, Dims: 1)"
 * "/entry0/sample/quantum_wanted_temp  (Float, Dims: 1)"
 * "/entry0/sample/rack_temperature  (Float, Dims: 1)"
 * "/entry0/sample/regulation_temperature  (Float, Dims: 1)"
 * "/entry0/sample/sampleId  (String, Dims: 1)"
 * "/entry0/sample/sample_changer_nickname  (String, Dims: 1)"
 * "/entry0/sample/sample_changer_slot_value  (Signed Integer, Dims: 1)"
 * "/entry0/sample/sample_changer_value  (Signed Integer, Dims: 1)"
 * "/entry0/sample/sample_pressure  (Float, Dims: 1)"
 * "/entry0/sample/sample_pressure2  (Float, Dims: 1)"
 * "/entry0/sample/sample_temperature2  (Float, Dims: 1)"
 * "/entry0/sample/san_actual  (Float, Dims: 1)"
 * "/entry0/sample/san_offset  (Float, Dims: 1)"
 * "/entry0/sample/san_requested  (Float, Dims: 1)"
 * "/entry0/sample/scatteringLengthDensity  (String, Dims: 1)"
 * "/entry0/sample/sdi_actual  (Float, Dims: 1)"
 * "/entry0/sample/sdi_offset  (Float, Dims: 1)"
 * "/entry0/sample/sdi_requested  (Float, Dims: 1)"
 * "/entry0/sample/setpoint_temperature  (Float, Dims: 1)"
 * "/entry0/sample/sht_actual  (Float, Dims: 1)"
 * "/entry0/sample/sht_offset  (Float, Dims: 1)"
 * "/entry0/sample/sht_requested  (Float, Dims: 1)"
 * "/entry0/sample/size  (String, Dims: 1)"
 * "/entry0/sample/str_actual  (Float, Dims: 1)"
 * "/entry0/sample/str_offset  (Float, Dims: 1)"
 * "/entry0/sample/str_requested  (Float, Dims: 1)"
 * "/entry0/sample/surfaceArea  (String, Dims: 1)"
 * "/entry0/sample/temperature  (Float, Dims: 1)"
 * "/entry0/sample/thermo_temperature1  (Float, Dims: 1)"
 * "/entry0/sample/thermo_temperature2  (Float, Dims: 1)"
 * "/entry0/sample/thickness  (Float, Dims: 1)"
 * "/entry0/sample/transmission  (Float, Dims: 1)"
 * "/entry0/sample/trs_actual  (Float, Dims: 1)"
 * "/entry0/sample/trs_offset  (Float, Dims: 1)"
 * "/entry0/sample/trs_requested  (Float, Dims: 1)"
 * "/entry0/sample/turbidity_voltage  (Float, Dims: 1)"
 * "/entry0/sample/typeInstall  (String, Dims: 1)"
 * "/entry0/sample/unitCellClass  (String, Dims: 1)"
 * "/entry0/sample/volumeUnitCell  (String, Dims: 1)"
 * "/entry0/sample_description  (String, Dims: 1)"
 * "/entry0/start_time  (String, Dims: 1)"
 * "/entry0/useListMode  (Signed Integer, Dims: 1)"
 * "/entry0/user/name  (String, Dims: 1)"
 * "/entry0/user/namelocalcontact  (String, Dims: 1)"
 * "/entry0/user/proposal  (String, Dims: 1)"
 * (237 Werte)
 */
