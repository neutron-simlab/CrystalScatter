/** Crystal calculation.
  *
  * A class with in interface for all parameters an the result and a calculation function
  * to calculate the image.
  * This was originally written in Pascal and translated into C++ by m.wagener@fz-juelich.de
  * During the translation, the variable and type names are preserved.
  *
  * ONLY GPU routines here!
  */


#include "sc_calc_generic_gpu.h"
#include <stdlib.h>
#include <string.h>

#define DBGSPEC(x) //x
// Diese Debugausgaben sind innerhalb der Schleifen und produzieren viel Output.
//  Werden aber sowohl für die CPU als auch für die GPU genutzt.

#include <iostream>
#include <chrono>
#include <unistd.h>
#include <signal.h>

#ifdef __CUDACC__
#include <cuda.h>
#endif


#define myacos(x) ((x)<-1)?0.0:(((x)>+1)?M_PI:acos((x)))



//#############################################################################################################
// Attention for CUDA version 12 on Linux:
//
//  The Lazy Loading feature (introduced in CUDA 11.7) is now enabled by default on Linux with the 535 driver.
//  To disable this feature on Linux, set the environment variable CUDA_MODULE_LOADING=EAGER before launch.
//  Default enablement for Windows will happen in a future CUDA driver release. To enable this feature on
//  Windows, set the environment variable CUDA_MODULE_LOADING=LAZY before launch.
//
// With this Lazy Loading activated, this program will not run (August 2023).
//#############################################################################################################


#ifdef __CUDACC__
#define KERNEL_PREWORK  int ihex_x, i_y; \
                        if ( dofit ) \
                        { \
                            ihex_x = blockIdx.x * blockDim.x + threadIdx.x; \
                            if ( ihex_x >= CALC._fitWidth ) return; \
                            i_y    = blockIdx.y * blockDim.y + threadIdx.y; \
                            if ( i_y    >= CALC._fitHeight ) return; \
                            if ( !CALC.checkFitRanges( CALC, ihex_x, i_y ) ) return; \
                        } else { \
                            ihex_x = blockIdx.x * blockDim.x + threadIdx.x + CALC.zzmin; \
                            if ( ihex_x >= CALC.zzmax ) return; \
                            i_y    = blockIdx.y * blockDim.y + threadIdx.y + CALC.iimin; \
                            if ( i_y    >= CALC.iimax ) return; \
                        } \
                        double qx, qy, qz; \
                        CALC.convert_idx2q( CALC, ihex_x, i_y, qx, qy, qz );

#define KERNEL_POSTWORK if ( CALC._endThread ) return; \
                        if ( dofit ) \
                        { \
                            size_t idx = CALC.IDX(ihex_x,i_y); \
                            if ( pixval > 0 && CALC.arrFitData[idx] > 0 ) \
                                CALC.arrFitFqs[idx] = FQSVERGL( pixval, CALC.arrFitData[idx] ); \
                            else \
                                CALC.arrFitFqs[idx] = 0.0; \
                        } else \
                            CALC.setXYIntensity( ihex_x, i_y, pixval );

#endif


#include "sc_gpu_generic.h"
#include "sc_gpu_pCyl_lNon_oGau.h"
#include "sc_gpu_pCyl_lNon_oZDi.h"
#include "sc_gpu_pCyl_pGau_lBCT_oZDi.h"
#include "sc_gpu_pSph_lNon_oIso.h"
#include "sc_gpu_pSph_pGau_lFCC_oGau.h"
#include "sc_gpu_pCub_lNon_oIso.h"


void SasCalc_GENERIC_calculation::kernel_selection( int Nx, int Ny, bool dofit )
{
#ifdef __CUDACC__
    _endThread = false;

    // GPU Programming with CUDA @ JSC (Kurs Apr.2022)
    // GPU 01 Introduction.pdf - Seite 52
    dim3 blockDimDef(16, 16);
    dim3 gridDimDef((Nx % blockDimDef.x == 0) ? (Nx / blockDimDef.x) : (Nx / blockDimDef.x + 1),
                    (Ny % blockDimDef.y == 0) ? (Ny / blockDimDef.y) : (Ny / blockDimDef.y + 1));
    dim3 blockDimSpec(16, 16);
    dim3 gridDimSpec((Nx % blockDimSpec.x == 0) ? (Nx / blockDimSpec.x) : (Nx / blockDimSpec.x + 1),
                     (Ny % blockDimSpec.y == 0) ? (Ny / blockDimSpec.y) : (Ny / blockDimSpec.y + 1));

    //std::cerr << "Nxy: " << Nx << "/" << Ny << ", gridDef: " << gridDimDef.x << "/" << gridDimDef.y
    //          << ", gridSpec: " << gridDimSpec.x << "/" << gridDimSpec.y << std::endl;

    bool done=false;
    if ( !bIgnoreNewSwitch )
    {
        done = true;
        // ComboBoxParticle mit 11 Varianten
        // ComboBoxPeak mit 8 Varianten (shp) nicht verwendet bei LatticeType=12
        // LatticeType mit 25 Varianten
        // Ordis-Parameter mit 14 Varianten
        // Eine 4-Byte-Zahl berechnen und dann nur ein Switch. Ist übersichtlicher im Code.
        switch ( (ComboBoxParticle << 24) | (ltype==12?0:(shp << 16)) | (ltype << 8) | params.ordis )
        {
        case 0x00000C07:  // cbpartSphere(0), shp(X), ltype(None=12), ordis_Isotropic(7)
            kernel_partSphere_lattNone_ordisIsotropic<<<gridDimSpec,blockDimSpec>>>(*this,dofit);
            break;
        case 0x00020500:  // cbpartSphere(0), cbpeakGaussian(2), ltype(FCC=5), ordis_Gaussian(0)
            kernel_partSphere_peakGauss_lattFCC_ordisGauss<<<gridDimSpec,blockDimSpec>>>(*this,dofit);
            break;
        case 0x01000C00:  // cbpartCylinder(1), shp(X), ltype(None=12), ordis_Gaussian(0):
            kernel_partCylinder_lattNone_ordisGauss<<<gridDimSpec,blockDimSpec>>>(*this,dofit);
            break;
        case 0x01000C06:  // cbpartCylinder(1), shp(X), ltype(None=12), ordis_ZDir(6):
            kernel_partCylinder_lattNone_ordisZDir<<<gridDimSpec,blockDimSpec>>>(*this,dofit);
            break;
        case 0x01020806:  // cbpartCylinder(1), cbpeakGaussian(2), ltype(BCT=8), ordis_ZDir(6):
            kernel_partCylinder_peakGauss_lattBCT_ordisZDir<<<gridDimSpec,blockDimSpec>>>(*this,dofit);
            break;
        case 0x04000C07:  // cbpartCube(4), shp(X), ltype(None=12), ordis_Isotropic(7)
            kernel_partCube_lattNone_ordisIsotropic<<<gridDimSpec,blockDimSpec>>>(*this,dofit);
            break;
        default:          // Alle anderen Kombinationen
            done = false;
            break;
        } // switch (code)
    } // if ( !bIgnoreNewSwitch )
    if ( ! done )
    {
        //std::cerr << "Use kernel_GENERIC " << (dofit?"Fit":"Calc") << std::endl;
        kernel_GENERIC<<<gridDimDef,blockDimDef>>>(*this,dofit);
    }
    setNewSwitchUsed(done);

    cudaError_t err = cudaGetLastError();
    if ( err != cudaSuccess )
        std::cerr << cudaGetErrorString(err) << std::endl;

    cudaDeviceSynchronize();
#else // __CUDACC__
    Q_UNUSED(Nx)
    Q_UNUSED(Ny)
    Q_UNUSED(dofit)
    std::cerr << "NO CUDA" << std::endl;
#endif // __CUDACC__
} /* kernel_selection() */


/**
 * @brief SasCalc_gpu::doIntCalc_GENERIC_F
 * @param CALC   - reference to class with parameters and subfunctions
 * @param dofit  - true:2d-Fit, false: normal image generation
 * @param ihex_x - horizontal pixel index
 * @param i_y    - vertical pixel index
 * This function is called ONLY from the CPU-Threads to calculate the pixelvalue and store it
 *  if no specialized function is defined (see top of file).
 */
void SasCalc_GENERIC_calculation::calcCPU_selection( const SasCalc_GENERIC_calculation& CALC,
                                                      bool dofit, int ihex_x, int i_y )
{
    if ( dofit )
    {
        if ( !CALC.checkFitRanges( CALC, ihex_x, i_y ) ) return;
    }

    // Bestimmen von qx,qy,qz
    double qx, qy, qz;
    convert_idx2q( CALC, ihex_x, i_y, qx, qy, qz );

    // Call q(x,y,z) Routine (also from 2dFit)
    double pixval;
    bool done=false;

    // Crashtest wegen 60sec Zeitüberschreitung (TODO) mit
    // "C:\SimLab\sas-crystal\20241206 - CombinationTests\Comb_lt00_pa02_or00_in00_pe07.ini"

    if ( !CALC.bIgnoreNewSwitch )
    {
        done = true;
        // 0xPP?????? ComboBoxParticle mit 11 Varianten
        // 0x??SS???? ComboBoxPeak mit 8 Varianten (shp), nicht verwendet bei LatticeType=12
        // 0x????LL?? LatticeType mit 25 Varianten
        // 0x??????OO Ordis-Parameter mit 14 Varianten
        // Eine 4-Byte-Zahl berechnen und dann nur ein Switch. Ist übersichtlicher im Code und wahrscheinlich auch
        //  etwas effizienter im Code als viele geschachtelte switch.
        switch ( (CALC.ComboBoxParticle << 24) | (ltype==12?0:(CALC.shp << 16)) | (ltype << 8) | params.ordis )
        {
        case 0x00000C07:  // cbpartSphere(0), shp(X), ltype(None=12), ordis_Isotropic(7)
            pixval = CALC.calc_partSphere_lattNone_ordisIsotropic( CALC, qx, qy, qz );
            break;
        case 0x00020500:  // cbpartSphere(0), cbpeakGaussian(2), ltype(FCC=5), ordis_Gaussian(0)
            pixval = CALC.calc_partSphere_peakGauss_lattFCC_ordisGauss( CALC, qx, qy, qz );
            break;
        case 0x01000C00:  // cbpartCylinder(1), shp(X), ltype(None=12), ordis_Gaussian(0):
            pixval = CALC.calc_partCylinder_lattNone_ordisGauss( CALC, qx, qy, qz );
            break;
        case 0x01000C06:  // cbpartCylinder(1), shp(X), ltype(None=12), ordis_ZDir(6):
            pixval = CALC.calc_partCylinder_lattNone_ordisZDir( CALC, qx, qy, qz );
            break;
        case 0x01020806:  // cbpartCylinder(1), cbpeakGaussian(2), ltype(BCT=8), ordis_ZDir(6):
            pixval = CALC.calc_partCylinder_peakGauss_lattBCT_ordisZDir( CALC, qx, qy, qz );
            break;
        case 0x04000C07:  // cbpartCube(4), shp(X), ltype(None=12), ordis_Isotropic(7)
            pixval = CALC.calc_partCube_lattNone_ordisIsotropic( CALC, qx, qy, qz );
            break;
        default:          // Alle anderen Kombinationen
            done = false;
            break;
        } // switch (code)
    } // if ( !CALC.bIgnoreNewSwitch )
    if ( ! done )
    {
        pixval = CALC.calc_GENERIC( CALC, qx, qy, qz );
    }

    if ( CALC._endThread ) return;
    CALC.setNewSwitchUsed(done);

    if ( dofit )
    {
        size_t idx = CALC.IDX(ihex_x,i_y); // (-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y));

        if ( pixval > 0 && CALC.arrFitData[idx] > 0 )
            CALC.arrFitFqs[idx] = FQSVERGL( pixval, CALC.arrFitData[idx] );
        else
            CALC.arrFitFqs[idx] = 0.0;
    }
    else
    {
        CALC.setXYIntensity( ihex_x, i_y, pixval );
    }
} /* calcCPU_selection() */



#ifdef __CUDACC__
__host__ __device__
#endif
/*static*/ void SasCalc_GENERIC_calculation::convert_idx2q( const SasCalc_GENERIC_calculation& CALC,
                                                        int ihex, int i, double &qx, double &qy, double &qz )
{
    if ( CALC.use1d )
    {
        qx = (CALC.qmax - CALC.qmin) / (double)(CALC.zmax) * (double)ihex;
        qy = 1e-20;
        qz = 1e-20;
        //qDebug() << "(1d) ihex/i" << ihex << i << "zmax" << CALC.zmax << "qx" << qx;
    }
    else if ( CALC.useBeamStop )
    {   // Einrechnen des Beamstops (d.h. Verschiebung des Zentrums)
        // Hier wird die Pixelsize aber nicht der qmax Wert beachtet.
        double mdet = (ihex/*+CALC.zmax+1*/)*CALC.pixnoy/(2.0*CALC.zmax);      /* mth pixel */
        double ndet = (i   /*+CALC.zmax+1*/)*CALC.pixnox/(2.0*CALC.zmax);      /* nth pixel */
        /*Z=24635*/
        // Im Pascal-Programm ist der Beamstop für Null in der Ecke.
        // Hier ist der Beamstop für Null in der Mitte gerechnet.
        double xdet = CALC.pixx_m * (ndet - CALC.beamX0);
        double ydet = CALC.pixy_m * (mdet - CALC.beamY0);
        double rdet = sqrt(xdet*xdet+ydet*ydet);
        double phidet = atan2(ydet,xdet);
        double thetadet = atan2(rdet,CALC.det);
        qx = 2*M_PI*cos(phidet)*sin(thetadet)/CALC.wave;
        qy = 2*M_PI*sin(phidet)*sin(thetadet)/CALC.wave;
        qz = 2*M_PI*(1-cos(thetadet))/CALC.wave;
        //if ( ihex == 10 && i >= 10 && i <= 20 )
        //    qDebug() << "ihex/i" << ihex << i << "BSx/y" << CALC.beamX0 << CALC.beamY0 << "zmax" << CALC.zmax << "m/ndet" << mdet << ndet << "qxy" << qx << qy;
    }
    else
    {   // Den Mittelpunkt nutzen (ihex=0,i=0)
        // hier wird nur der qmax Wert beachtet und die Pixelsize außer Acht gelassen.
        qx = CALC.qmax * i    / (double)(CALC.zmax); // lamu
        qy = CALC.qmax * ihex / (double)(CALC.zmax); // lamv
        qz = 1e-20;
        //if ( ihex == 10 && i >= 10 && i <= 20 )
        //    qDebug() << "ihex/i" << ihex << i << "BSx/y" << CALC.beamX0 << CALC.beamY0 << "zmax" << CALC.zmax << "qxy" << qx << qy;
    }
}


#ifdef __CUDACC__
__host__ __device__
#endif
/*static*/ bool SasCalc_GENERIC_calculation::checkFitRanges( const SasCalc_GENERIC_calculation& CALC,
                                                int &ihex_x, int &i_y )
{
    // There are up to 4 rectangle regions to be ignored during the fit.
    for ( int i=0; i<4; i++ )
    {
        if ( CALC.noFitX0[i] < 0 ) continue;
        if ( ihex_x < CALC.noFitX0[i] ) continue;
        if ( CALC.noFitX1[i] < ihex_x ) continue;
        if ( i_y < CALC.noFitY0[i] ) continue;
        if ( CALC.noFitY1[i] < i_y ) continue;
        size_t idx = ihex_x + (CALC._fitWidth * i_y);
        CALC.arrFitFqs[idx] = 0.0;
        return false;
    }

    if ( CALC.fitBorderPixel > 0 || CALC.fitBStopPixel > 0 )
    {   // Ausblendungen über Pixelangaben an Rand und Mitte
        if (ihex_x < CALC._xmin + CALC.fitBorderPixel || ihex_x >= CALC._xmax - CALC.fitBorderPixel ||
            i_y < CALC._ymin + CALC.fitBorderPixel || i_y >= CALC._ymax - CALC.fitBorderPixel)
        {
            size_t idx = ihex_x + (CALC._fitWidth * i_y);
            CALC.arrFitFqs[idx] = 0.0;
            return false;
        }
        if (ihex_x >= CALC.beamX0 - CALC.fitBStopPixel &&
            ihex_x <  CALC.beamX0 + CALC.fitBStopPixel &&
            i_y    >= CALC.beamY0 - CALC.fitBStopPixel &&
            i_y    <  CALC.beamY0 + CALC.fitBStopPixel )
        {
            size_t idx = ihex_x + (CALC._fitWidth * i_y);
            CALC.arrFitFqs[idx] = 0.0;
            return false;
        }
    }
    else if ( CALC.fitBStopPixel == -1 )
    {   // Ausblendungen per Eck-Pixel-Wert
        size_t idx = ihex_x + (CALC._fitWidth * i_y);
        if ( CALC.arrFitData[idx] <= CALC.arrFitData[0] )
        {
            CALC.arrFitFqs[idx] = 0.0;
            return false;
        }
    }
    //else
    //{   // Keine Ausblenungen
    //}

    // int i    = x+CALC.iimin; // - CALC.beamX0;
    // int ihex = y+CALC.iimin; // - CALC.beamY0;

    ihex_x += CALC.iimin;
    i_y    += CALC.iimin;
    return true;
}


#ifdef __CUDACC__
__host__ __device__
#endif
void SasCalc_GENERIC_calculation::setXYIntensity( int x/*ihex*/, int y/*i*/, double val ) const
{
    size_t idx = IDX(x,y); // (-_xmin + (x)) + (_xmax-_xmin/*+1*/)*(-_ymin + (y));
    if ( idx >= _arrCount )
    {
#ifndef __CUDACC__
        qDebug() << "MEMFAIL" << "x"<<x << "y"<<y << "idx"<<idx << "cnt"<<_arrCount
                 << "xmin"<<_xmin << "xmax"<<_xmax << "ymin"<<_ymin << "ymax"<<_ymax;
#endif
        arrXYIntensity[0] = -1;
        arrXYIntensity[1] = idx;
        return;
    }
    arrXYIntensity[idx+3] = val;
    arrXYIntensity[0] = x;      // Merker für letztes gerechnetes Pixel, wird vom Loggig-Thread
    arrXYIntensity[1] = y;      //  genutzt, um den Fortschritt anzuzeigen.
    //arrXYIntensity[2] ist ein Debug-Flag
}


// Memory-Functions need some Cuda-Code!
void SasCalc_GENERIC_calculation::initMemory()
{
    arrXYIntensity = nullptr;
    arrXYsize      = 0;
#ifdef FITDATA_IN_GPU
    arrFitData  = nullptr;
    arrFitFqs   = nullptr;
    _fitWidth   = 0;
    _fitHeight  = 0;
    _arrFitSize = 0;
    _fitEnabled = false;
#endif
    _xmin=0;
    _xmax=0;
    _ymin=0;
    _ymax=0;
#ifdef __CUDACC__
    std::cerr << "vor cudaDeviceReset()" << std::endl;
    cudaDeviceReset();
    int deviceCount = 0;
    std::cerr << "vor cudaGetDeviceCount()" << std::endl;
    cudaGetDeviceCount(&deviceCount);
    std::cerr << "SasCalc_GENERIC_calculation::initMemory(): GPU device count: " << deviceCount << std::endl;
    for (int device = 0; device < deviceCount; ++device)
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);
        std::cerr << "  Device " << device << " (" << deviceProp.name << ") has compute capability " << deviceProp.major << "." << deviceProp.minor << std::endl;
        std::cerr << "                    canMapHostMemory: " << deviceProp.canMapHostMemory << std::endl;
        std::cerr << "                      totalGlobalMem: " << deviceProp.totalGlobalMem << " bytes" << std::endl;
        std::cerr << "                  maxThreadsPerBlock: " << deviceProp.maxThreadsPerBlock << std::endl;
        std::cerr << "                       maxThreadsDim: " << deviceProp.maxThreadsDim[0] << "," << deviceProp.maxThreadsDim[1] << "," << deviceProp.maxThreadsDim[2] << std::endl;
        std::cerr << "                         maxGridSize: " << deviceProp.maxGridSize[0] << "," << deviceProp.maxGridSize[1] << "," << deviceProp.maxGridSize[2] << std::endl;
        std::cerr << "                        regsPerBlock: " << deviceProp.regsPerBlock << std::endl; // is the maximum number of 32-bit registers available to a thread block;
        std::cerr << "                 multiProcessorCount: " << deviceProp.multiProcessorCount << std::endl;
        std::cerr << "         maxThreadsPerMultiProcessor: " << deviceProp.maxThreadsPerMultiProcessor << std::endl;
        std::cerr << "                         l2CacheSize: " << deviceProp.l2CacheSize << " bytes" << std::endl;
        std::cerr << "            persistingL2CacheMaxSize: " << deviceProp.persistingL2CacheMaxSize << " bytes" << std::endl;
        std::cerr << "      directManagedMemAccessFromHost: " << deviceProp.directManagedMemAccessFromHost << std::endl;

        /* - \ref ::cudaDeviceProp::uuid "uuid" is a 16-byte unique identifier.
        * - \ref ::cudaDeviceProp::sharedMemPerBlock "sharedMemPerBlock" is the maximum amount of shared memory available to a thread block in bytes;
        * - \ref ::cudaDeviceProp::warpSize "warpSize" is the warp size in threads;
        * - \ref ::cudaDeviceProp::memPitch "memPitch" is the maximum pitch in bytes allowed by the memory copy functions that involve memory regions allocated through ::cudaMallocPitch();
        * - \ref ::cudaDeviceProp::clockRate "clockRate" is the clock frequency in kilohertz;
        * - \ref ::cudaDeviceProp::totalConstMem "totalConstMem" is the total amount of constant memory available on the device in bytes;
        * - \ref ::cudaDeviceProp::kernelExecTimeoutEnabled "kernelExecTimeoutEnabled" is 1 if there is a run time limit for kernels executed on the device, or 0 if not.
        * - \ref ::cudaDeviceProp::integrated "integrated" is 1 if the device is an integrated (motherboard) GPU and 0 if it is a discrete (card) component.
        * - \ref ::cudaDeviceProp::computeMode "computeMode" is the compute mode that the device is currently in. Available modes are as follows:
        *   - cudaComputeModeDefault: Default mode - Device is not restricted and multiple threads can use ::cudaSetDevice() with this device.
        *   - cudaComputeModeExclusive: Compute-exclusive mode - Only one thread will be able to use ::cudaSetDevice() with this device.
        *   - cudaComputeModeProhibited: Compute-prohibited mode - No threads can use ::cudaSetDevice() with this device.
        *   - cudaComputeModeExclusiveProcess: Compute-exclusive-process mode - Many threads in one process will be able to use ::cudaSetDevice() with this device.
        *   If ::cudaSetDevice() is called on an already occupied \p device with computeMode ::cudaComputeModeExclusive, ::cudaErrorDeviceAlreadyInUse
        *   will be immediately returned indicating the device cannot be used. When an occupied exclusive mode device is chosen with ::cudaSetDevice,
        *   all subsequent non-device management runtime functions will return ::cudaErrorDevicesUnavailable.
        * - \ref ::cudaDeviceProp::concurrentKernels "concurrentKernels" is 1 if the device supports executing multiple kernels within the same context simultaneously, or 0 if not.
                   It is not guaranteed that multiple kernels will be resident on the device concurrently so this feature should not be relied upon for correctness;
        * - \ref ::cudaDeviceProp::ECCEnabled "ECCEnabled" is 1 if the device has ECC support turned on, or 0 if not.
        * - \ref ::cudaDeviceProp::pciBusID "pciBusID" is the PCI bus identifier of the device.
        * - \ref ::cudaDeviceProp::pciDeviceID "pciDeviceID" is the PCI device (sometimes called slot) identifier of the device.
        * - \ref ::cudaDeviceProp::pciDomainID "pciDomainID" is the PCI domain identifier of the device.
        * - \ref ::cudaDeviceProp::tccDriver "tccDriver" is 1 if the device is using a TCC driver or 0 if not.
        * - \ref ::cudaDeviceProp::asyncEngineCount "asyncEngineCount" is 1 when the device can concurrently copy memory between host and device while executing a kernel.
                   It is 2 when the device can concurrently copy memory between host and device in both directions and execute a kernel at the same time. It is 0 if neither of these is supported.
        * - \ref ::cudaDeviceProp::unifiedAddressing "unifiedAddressing" is 1 if the device shares a unified address space with the host and 0 otherwise.
        * - \ref ::cudaDeviceProp::memoryClockRate "memoryClockRate" is the peak memory clock frequency in kilohertz.
        * - \ref ::cudaDeviceProp::memoryBusWidth "memoryBusWidth" is the memory bus width in bits.
        * - \ref ::cudaDeviceProp::streamPrioritiesSupported "streamPrioritiesSupported" is 1 if the device supports stream priorities, or 0 if it is not supported.
        * - \ref ::cudaDeviceProp::globalL1CacheSupported "globalL1CacheSupported" is 1 if the device supports caching of globals in L1 cache, or 0 if it is not supported.
        * - \ref ::cudaDeviceProp::localL1CacheSupported "localL1CacheSupported" is 1 if the device supports caching of locals in L1 cache, or 0 if it is not supported.
        * - \ref ::cudaDeviceProp::sharedMemPerMultiprocessor "sharedMemPerMultiprocessor" is the maximum amount of shared memory available to a multiprocessor in bytes; this amount is shared by all thread blocks simultaneously resident on a multiprocessor;
        * - \ref ::cudaDeviceProp::regsPerMultiprocessor "regsPerMultiprocessor" is the maximum number of 32-bit registers available to a multiprocessor; this number is shared by all thread blocks simultaneously resident on a multiprocessor;
        * - \ref ::cudaDeviceProp::managedMemory "managedMemory" is 1 if the device supports allocating managed memory on this system, or 0 if it is not supported.
        * - \ref ::cudaDeviceProp::isMultiGpuBoard "isMultiGpuBoard" is 1 if the device is on a multi-GPU board (e.g. Gemini cards), and 0 if not;
        * - \ref ::cudaDeviceProp::multiGpuBoardGroupID "multiGpuBoardGroupID" is a unique identifier for a group of devices associated with the same board. Devices on the same multi-GPU board will share the same identifier;
        * - \ref ::cudaDeviceProp::singleToDoublePrecisionPerfRatio "singleToDoublePrecisionPerfRatio" is the ratio of single precision performance (in floating-point operations per second) to double precision performance.
        * - \ref ::cudaDeviceProp::pageableMemoryAccess "pageableMemoryAccess" is 1 if the device supports coherently accessing pageable memory without calling cudaHostRegister on it, and 0 otherwise.
        * - \ref ::cudaDeviceProp::concurrentManagedAccess "concurrentManagedAccess" is 1 if the device can coherently access managed memory concurrently with the CPU, and 0 otherwise.
        * - \ref ::cudaDeviceProp::computePreemptionSupported "computePreemptionSupported" is 1 if the device supports Compute Preemption, and 0 otherwise.
        * - \ref ::cudaDeviceProp::canUseHostPointerForRegisteredMem "canUseHostPointerForRegisteredMem" is 1 if the device can access host registered memory at the same virtual address as the CPU, and 0 otherwise.
        * - \ref ::cudaDeviceProp::cooperativeLaunch "cooperativeLaunch" is 1 if the device supports launching cooperative kernels via ::cudaLaunchCooperativeKernel, and 0 otherwise.
        * - \ref ::cudaDeviceProp::cooperativeMultiDeviceLaunch "cooperativeMultiDeviceLaunch" is 1 if the device supports launching cooperative kernels via ::cudaLaunchCooperativeKernelMultiDevice, and 0 otherwise.
        * - \ref ::cudaDeviceProp::pageableMemoryAccessUsesHostPageTables "pageableMemoryAccessUsesHostPageTables" is 1 if the device accesses pageable memory via the host's page tables, and 0 otherwise.
        * - \ref ::cudaDeviceProp::directManagedMemAccessFromHost "directManagedMemAccessFromHost" is 1 if the host can directly access managed memory on the device without migration, and 0 otherwise.
        * - \ref ::cudaDeviceProp::maxBlocksPerMultiProcessor "maxBlocksPerMultiProcessor" is the maximum number of thread blocks that can reside on a multiprocessor.
        * - \ref ::cudaDeviceProp::accessPolicyMaxWindowSize "accessPolicyMaxWindowSize" is the maximum value of ::cudaAccessPolicyWindow::num_bytes.
        */

        // IFF1585:
        //  Device 0 (GeForce RTX 2070 SUPER) has compute capability 7.5
        //                    canMapHostMemory: 1
        //                      totalGlobalMem: 8366915584 bytes (=8GB)
        //                  maxThreadsPerBlock: 1024
        //                       maxThreadsDim: 1024,1024,64
        //                         maxGridSize: 2147483647,65535,65535
        //                 multiProcessorCount: 40
        //         maxThreadsPerMultiProcessor: 1024
        //                         l2CacheSize: 4194304 bytes
        //            persistingL2CacheMaxSize: 0 bytes
        //      directManagedMemAccessFromHost: 0

        //cudaCtxResetPersistingL2Cache();
    }

    //size_t s;
    //cudaDeviceGetLimit( &s, cudaLimitPrintfFifoSize );
    //std::cerr << "                    Printf Fifo Size: " << s << std::endl;

    //size_t memFree, memTotal;
    //cudaMemGetInfo( &memFree, &memTotal );
    //std::cerr << "Memory: free=" << memFree << " of " << memTotal << std::endl;

    noGPUavailable = deviceCount == 0;
#else
    noGPUavailable = true;
    std::cerr << "SasCalc_GENERIC_calculation::initMemory(): no nvcc used" << std::endl;
#endif
}

// Memory-Functions need some Cuda-Code!
void SasCalc_GENERIC_calculation::createMemory( void **ptr, size_t lensoll, size_t &lenist, bool gpuonly, const char *dbgInfo )
{
#ifndef __CUDACC__
    //Q_UNUSED(dbgInfo);
    Q_UNUSED(gpuonly);

#else

//#define MEMALLOC(p,s) cudaMallocManaged(p,s)      //ist wohl doch etwas langsamer...
//#define MEMFREE(p)    cudaFree(p)                 //und nach einem Crash wird nicht richtig aufgeräumt

#define MEMALLOC(p,s) cudaMallocHost(p,s)
#define MEMFREE(p)    cudaFreeHost(p)

#endif
    if ( *ptr != nullptr )
    {
        if ( lensoll <= lenist ) return;    // nur, wenn wirklich mehr gebraucht wird.
        // realloc new size
        if ( noGPUavailable )
        {
            std::cerr << "re-allocate " << lensoll << " Bytes CPU-Memory (ist=" << lenist << ")";
            if ( dbgInfo != nullptr ) std::cerr << " (" << dbgInfo << ")";
            std::cerr << std::endl;
            *ptr = realloc( *ptr, lensoll );
        }
#ifdef __CUDACC__
        else
        {
            std::cerr << "re-allocate " << lensoll << " Bytes GPU-Memory (ist=" << lenist << ")";
            if ( dbgInfo != nullptr ) std::cerr << " (" << dbgInfo << ")";
            std::cerr << std::endl;
            // Da der Speicher immer vor den Berechnungen angelegt wird und die alten Daten
            // nicht verwendet werden (hier nur das Image), muss nichts kopiert werden.
            // Das klappte nämlich nicht wirklich - zumindest bei den ersten Versuchen
            cudaError_t e;
            e = MEMFREE( *ptr );
            if ( e != cudaSuccess )
            {
                std::cerr << "GPU free memory failed:" << cudaGetErrorString(e) << std::endl;
            }
            e = MEMALLOC( ptr, lensoll );
            if ( *ptr == nullptr )
            {
                std::cerr << "GPU allocate memory failed, ignore GPU features" << std::endl;
                std::cerr << "Error is: " << cudaGetErrorString(e) << std::endl;
                *ptr = malloc( lensoll );
                noGPUavailable = true;
            }
        }
#endif
    }
    else
    {   // alloc memory
        if ( noGPUavailable )
        {
            std::cerr << "allocate " << lensoll << " Bytes CPU-Memory";
            if ( dbgInfo != nullptr ) std::cerr << " (" << dbgInfo << ")";
            std::cerr << std::endl;
            *ptr = malloc( lensoll );
        }
#ifdef __CUDACC__
        else
        {
            std::cerr << "allocate " << lensoll << " Bytes GPU-Memory";
            if ( dbgInfo != nullptr ) std::cerr << " (" << dbgInfo << ")";
            std::cerr << std::endl;
            cudaError_t e = MEMALLOC( ptr, lensoll );
            if ( *ptr == nullptr )
            {
                std::cerr << "GPU allocate memory failed, ignore GPU features" << std::endl;
                std::cerr << "Error is: " << cudaGetErrorString(e) << std::endl;
                *ptr = malloc( lensoll );
                noGPUavailable = true;
            }
            if ( e != cudaSuccess )
                std::cerr << cudaGetErrorString(e) << std::endl;
        }
#endif
    }
    if ( *ptr != nullptr )
    {
        memset( *ptr, 0, lensoll );
        lenist = lensoll;
    }
    else
    {
        std::cerr << "   allocate FAILED" << std::endl;
        lenist = 0;
    }

}

// Memory-Functions need some Cuda-Code!
void SasCalc_GENERIC_calculation::memcleanup( void *arr )
{
    if ( arr != nullptr )
    {
        if ( noGPUavailable )
            free( static_cast<void*>( arr ) );
#ifdef __CUDACC__
        else
            MEMFREE( arr );
#endif
    }
}



//(************************** Gamma function ***************************)
/**
 * @brief SasCalc_GENERIC_calculation::gamma
 * @param z
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::gamma(double z) const          //{NV}
{
    // Die einzelnen Zahlen werden in Arrays umkopiert. Daher hier direkt die Arrys fertigstellen.
    // Dabei ist der erste Wert immer 0, da die Pascal-Arrays bei 1 beginnen und ich nicht alle Indizes umsetzen will.
    const double bm[] = { 0.0,
                         -0.577191652, 0.988205891, -0.897056937, 0.918206857, -0.756704078, 0.482199394, -0.193527818, 0.035868343 };
    const double cm[] = {  0.0,
                         1.0000000000000000,  0.5772156649015329, -0.6558780715202538,  -0.0420026350340952,  0.1665386113822915,
                         -0.0421977345555443, -0.0096219715278770,  0.00721894322466630, -0.0011651675918591, -0.0002152416741149,
                         0.0001280502823882, -0.0000201348547807, -0.0000012504934821,   0.0000011330272320, -0.0000002056338417,
                         0.0000000061160950,  0.0000000050020075, -0.0000000011812746,   0.0000000001043427,  0.0000000000077823,
                         -0.0000000000036968,  0.0000000000005100, -0.0000000000000206,  -0.0000000000000054,  0.0000000000000014,
                         0.0000000000000001 };

    bool negative;

    if ( fabs(z) < eps9 ) z = z + eps9;
    if ( z > 0 )
        negative = false;
    else //if ( z < 0 )
    {
        negative = true;
        if ( fabs(frac(z)) < eps9 ) z = z+eps9;
    }
    const double absz = fabs(z);

    if ( absz > 0.0 && absz <= 1 )
    {
        double fct = 0;
        for ( int i=1; i<=26; i++ ) fct += cm[i]*pow(absz,i);
        if ( negative ) return -M_PI*fct/(absz*sin(M_PI*absz));
        return 1./fct;
    }
    if ( absz > 1 && absz <= 1000 )
    {
        if ( fabs(frac(absz)) < eps9 )  // (* integer z *)
        {
            double fct = 1;
            for ( int i=1; i<=trunc(absz)-1; i++ ) fct = fct*i;
            if ( negative ) return -M_PI/(absz*fct*sin(M_PI*absz));
            return fct;
        }
        const double y = absz-1;
        const double di = trunc(y);
        const double x = absz-trunc(absz);
        double fak = y;
        if ( di == 0 ) fak = 1;
        else for ( int i=1; i<di; i++ ) fak = fak*(y-i);
        double f = 1;
        for ( int i=1; i<=8; i++ ) f = f+bm[i]*pow(x,i);
        if ( negative ) return -M_PI/(absz*fak*f*sin(M_PI*absz));
        return fak*f;
    }
    //if ( absz > 1000 )
    //{
    double fct = exp(-absz)*exp((absz-0.5)*log(absz))*sqrt(2*M_PI);
    if ( negative ) return -M_PI/(absz*fct*sin(M_PI*absz));
    return fct;
    //}
    //return 0.0;
}


/********************* gamma(z+a)/gamma(z+b) for large z *******************/
double SasCalc_GENERIC_calculation::gammaratio( double a, double b, double x) const
{
   const double c1=pow(1.+a/x, x+a-0.5);
   const double c2=pow(1.+b/x, x+b-0.5);
   return exp(-a+b)*exp((a-b)*log(x))*c1/c2;
}


//(* *************************** vesicles **************************** *)
/**
 * @brief SasCalc_GENERIC_calculation::polyvesicle
 * @param ro     = SasCalc_GENERIC_calculation::length
 * @param ri     = SasCalc_GENERIC_calculation::radius
 * @param sigmar = SasCalc_GENERIC_calculation::sigma
 * @param sigmal = SasCalc_GENERIC_calculation::sigmal
 * @param q      = local var from calling routine
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::polyvesicle(double q) const
{
    const double r = (params.length+params.radius)/2.;
    if ( q*r < 0.001 )
        return 1;

    const double d = (params.length-params.radius)/2.;
    const double zr = (1.-sqr(params.sigma))/sqr(params.sigma);
    const double zd = (1.-sqr(params.sigmal))/sqr(params.sigmal);
    const double yr = q*r/(zr+1.);
    const double yd = q*d/(zd+1.);
    const double c = 1.+d*d/(12.*r*r);

    const double gr1 = 1./(2.*zr*(zr-1.)*yr*yr);
    const double gr2 = gr1/((zr-2.)*yr);
    const double gr3 = gr2/((zr-3.)*yr);
    const double gd1 = 1./(2.*zd*(zd-1.)*yd*yd);
    const double gd2 = 1./(2.*zd*yd);
    const double gd3 = 1./2.;

    const double hr1 = gr1*(1.-cos((zr-2.+1.)*atan(2.*yr))/exp((zr-2.+1.)*log(1.+4.*yr*yr)/2.));    //(* <sin^2(yr)yr^-2> *)
    const double hr2 = gr2*(   sin((zr-3.+1.)*atan(2.*yr))/exp((zr-3.+1.)*log(1.+4.*yr*yr)/2.));    //(* <sin(yr)cos(yr)yr^-3> *)
    const double hr3 = gr3*(1.+cos((zr-4.+1.)*atan(2.*yr))/exp((zr-4.+1.)*log(1.+4.*yr*yr)/2.));    //(* <cos^2(yr)yr^-4> *)
    const double hd1 = gd1*(1.-cos((zd-2.+1.)*atan(2.*yd))/exp((zd-2.+1.)*log(1.+4.*yd*yd)/2.));    //(* <sin^2(yd)yd^-2> *)
    const double hd2 = gd2*(   sin((zd-1.+1.)*atan(2.*yd))/exp((zd-1.+1.)*log(1.+4.*yd*yd)/2.));    //(* <sin(yd)cos(yd)yd^-1> *)
    const double hd3 = gd3*(1.+cos((zd+0.+1.)*atan(2.*yd))/exp((zd+0.+1.)*log(1.+4.*yd*yd)/2.));    //(* <cos^2(yd)yd^0> *)

    return (hr1*hd1-2.*hr2*hd2+2.*hr2*hd1+hr3*hd3-2.*hr3*hd2+hr3*hd1)/(c*c);
}


/**
 * @brief SasCalc_GENERIC_calculation::f2dpolyvesicle
 * @param ro     = SasCalc_GENERIC_calculation::length
 * @param ri     = SasCalc_GENERIC_calculation::radius
 * @param sigmar = SasCalc_GENERIC_calculation::sigma
 * @param sigmad = ???  SasCalc_GENERIC_calculation::sigmal
 * @param q      = local var from calling routine
 * @return
 */
double SasCalc_GENERIC_calculation::f2dpolyvesicle(double q) const
{
    /* *************************** vesicles **************************** */
    double r,d,yr,yd,zr,zd,c,pq;
    double gr1,gr2,/*gr3,*/gd1,gd2,/*gd3,*/hr1,hr2,hr3,hd1,hd2,hd3;

    const double ro = params.length;
    const double ri = params.radius;
    const double sigmar = params.sigma;
    const double sigmad = params.sigmal;

    r=(ro+ri)/2;
    d=(ro-ri)/2;
    zr=(1-sqr(sigmar))/sqr(sigmar);
    zd=(1-sqr(sigmad))/sqr(sigmad);
    yr=q*r/(zr+1);
    yd=q*d/(zd+1);
    c=1+d*d/(12*r*r);

    if ( q*r<0.001 )
        pq=1;
    else
    {
        if ( sigmar<0.05 )
        {
            gr1=gammaratio(-1+1,1,zr)*exp(-1*log(yr));
            gr2=gammaratio(-2+1,1,zr)*exp(-2*log(yr));
            //gr3=gr2;
        }
        else
        {
            gr1=gamma(zr-1+1)*exp(-1*log(yr))/gamma(zr+1);       /* <sin(yr)yr^-1> */
            gr2=gamma(zr-2+1)*exp(-2*log(yr))/gamma(zr+1);       /* <cos(yr)yr^-2> */
            //gr3=gr2;                                            /* <cos(yr)yr^-2> */
        }
        if ( sigmad<0.05 )
        {
            gd1=gammaratio(-1+1,1,zd)*exp(-1*log(yd));
            gd2=gammaratio(-0+1,1,zd);
            //gd3=gd1;
        }
        else
        {
            gd1=gamma(zd-1+1)*exp(-1*log(yd))/gamma(zd+1);       /* <sin(yd)yd^-1> */
            gd2=gamma(zd-0+1)/gamma(zd+1);                      /* <cos(yd)yd^0>  */
            //gd3=gd1;                                            /* <sin(yd)yd^-1> */
        }

        hr1=gr1*sin((zr-1+1)*atan(yr))/exp((zr-1+1)*log(1+yr*yr)/2);            /* <sin(yr)yr^-1> */
        hr2=gr2*cos((zr-2+1)*atan(yr))/exp((zr-2+1)*log(1+yr*yr)/2);            /* <cos(yr)yr^-2> */
        hr3=hr2;                                                                /* <cos(yr)yr^-2> */
        hd1=gd1*sin((zd-1+1)*atan(yd))/exp((zd-1+1)*log(1+yd*yd)/2);            /* <sin(yd)yd^-1> */
        hd2=gd2*cos((zd-0+1)*atan(yd))/exp((zd-0+1)*log(1+yd*yd)/2);            /* <cos(yd)yd^0> */
        hd3=hd1;                                                                /* <sin(yd)yd^-1> */

        pq=(hr1*hd1-hr2*hd2+hr3*hd3)/c;
    }
    return pq*pq;
}



//(* ************************** 3d-integral over lorentz(x)*sin(x) ********************************* *)
/**
 * @brief SasCalc_GENERIC_calculation::lorentznorm3
 * @param a
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::lorentznorm3(double a) const
{
    if ( a > 0.01 )
    {
        const double b0 = log(1+a*M_PI*M_PI)/2.;
        const double b1 = (a*M_PI*M_PI-log(1+a*M_PI*M_PI))/12.;
        const double b2 = (-2*a*M_PI*M_PI+a*a*M_PI*M_PI*M_PI*M_PI+2*log(1+a*M_PI*M_PI))/480.;
        const double b3 = (a*M_PI*M_PI*(6-3*a*M_PI*M_PI+2*a*a*M_PI*M_PI*M_PI*M_PI)-6*log(1+a*M_PI*M_PI))/60480.;
        return b0/a-b1/(a*a)+b2/(a*a*a)-b3/(a*a*a*a);
    }
    //if a<=0.01 then begin
    const double a0 = 2;
    const double a1 = (M_PI*M_PI-4);
    const double a2 = (18-12*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI)/2.;
    const double a3 = (-1440+360*M_PI*M_PI-30*M_PI*M_PI*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI*M_PI*M_PI)/6.;
    return a0-a1*a+a2*a*a-a3*a*a*a;
}



//(* ************************** 3d-integral over gauss(x)*sin(x) ********************************* *)
/**
 * @brief SasCalc_GENERIC_calculation::gaussnorm3
 * @param a
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::gaussnorm3(double a) const         //{NV}
{
    if ( a > 0.01 )
    {
        const double b0 = (1-exp(-a*M_PI*M_PI))/2.;
        const double b1 = (1-exp(-a*M_PI*M_PI)*(1+a*M_PI*M_PI))/12.;
        const double b2 = (2-exp(-a*M_PI*M_PI)*(2+2*a*M_PI*M_PI+a*a*M_PI*M_PI*M_PI*M_PI))/240.;
        const double b3 = (6-exp(-a*M_PI*M_PI)*(6+6*a*M_PI*M_PI+3*a*a*M_PI*M_PI*M_PI*M_PI+a*a*a*M_PI*M_PI*M_PI*M_PI*M_PI*M_PI))/10080.;
        return b0/a-b1/(a*a)+b2/(a*a*a)-b3/(a*a*a*a);
    }
    //if a<=0.01 then begin
    const double a0 = 2;
    const double a1 = (M_PI*M_PI-4);
    const double a2 = (18-12*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI)/2.;
    const double a3 = (-1440+360*M_PI*M_PI-30*M_PI*M_PI*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI*M_PI*M_PI)/6.;
    return a0-a1*a+a2*a*a-a3*a*a*a;
}



//(* ************************** 3d-integral over pearson(x)*sin(x) ********************************* *)
/**
 * @brief SasCalc_GENERIC_calculation::pearsonnorm3
 * @param a
 * @param b
 * @return
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::pearsonnorm3(double a, double b) const         //{NV}
{
    b += 0.000001;
    if ( a > 20. )
    {
        const double b0 = (1-exp((1-b)*log(1+a*M_PI*M_PI)))/(2*(b-1));
        //b1 = (-1+a*a*M_PI*M_PI-b*a*M_PI*M_PI*(1+a*M_PI*M_PI)+exp(b*log(1+a*M_PI*M_PI)))/(12*(b-2)*(b-1)*exp(b*log(1+a*M_PI*M_PI)));
        //b2 = (-2*a*M_PI*M_PI+a*a*M_PI*M_PI*M_PI*M_PI+2*log(1+a*M_PI*M_PI))/480;
        //b3 = (a*M_PI*M_PI*(6-3*a*M_PI*M_PI+2*a*a*M_PI*M_PI*M_PI*M_PI)-6*log(1+a*M_PI*M_PI))/60480;
        return b0/a;
    }
    if ( a <= 0.01 )
    {
        const double a0 = 2.;
        const double a1 = (M_PI*M_PI-4)*b;
        const double a2 = (18-12*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI)*b*(b+1)/2.;
        const double a3 = (-1440+360*M_PI*M_PI-30*M_PI*M_PI*M_PI*M_PI+M_PI*M_PI*M_PI*M_PI*M_PI*M_PI)*b*(b+1)*(b+2)/6.;
        return a0-a1*a+a2*a*a-a3*a*a*a;
    }

    const int    maxit = 10;
    int    trapzditx;
    double inttf=0;
    pearsonintegral3(0,M_PI,a,b,inttf,1,trapzditx);
    double oldinttf=inttf;
    for ( int itf=2; itf<=maxit; itf++ )
    {
        pearsonintegral3(0,M_PI,a,b,inttf,itf,trapzditx);
        if ( fabs(1-inttf/oldinttf) < eps3 ) break;
        oldinttf = inttf;
    }
    return inttf;
}



//(* ********************* integration procedure for pearson(x)*sin(x) ***************************** *)
/**
 * @brief SasCalc_GENERIC_calculation::pearsonintegral3
 * @param at
 * @param bt
 * @param a
 * @param b
 * @param sx
 * @param nx
 * @param trapzditx
 */
#ifdef __CUDACC__
__host__ __device__
#endif
void SasCalc_GENERIC_calculation::pearsonintegral3(double at, double bt, double a, double b,
                                                  double &sx, int nx, int &trapzditx) const         //{NV}
{
    if ( nx == 1 )
    {
        const double fa = sin(at)/exp(b*log(1+a*at*at));
        const double fb = sin(bt)/exp(b*log(1+a*bt*bt));
        sx = 0.5*(bt-at)*(fa+fb);
        trapzditx = 1;
    }
    else
    {
        const double tnmx = trapzditx;
        const double delx = (bt-at)/tnmx;
        double xt = at+0.5*delx;
        double sumx = 0.0;
        for ( int jx=1; jx<=trapzditx; jx++ )
        {
            sumx += sin(xt)/exp(b*log(1+a*xt*xt));
            //sumx += fx;
            xt = xt+delx;
        }
        sx = 0.5*(sx+(bt-at)*sumx/tnmx);
        trapzditx = 2*trapzditx;
    }
}





//(* *********************** Romberg integration ****************************** *)
//(* returns integral in the limits a and b *)

//(*** integration routine use trapezoidal rule ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void SasCalc_GENERIC_calculation::trapzddeltac( double a, double b, double l, double r, double p1, double sigma, double alfa,
                                              double dbeta, double theta, double phi,
                                              double qx/*1*/, double qy/*1*/, double qz/*1*/,
                                              double p11, double p12, double p13, double p21, double p22, double p23,
                                              double p31, double p32, double p33,
                                              double qxn/*9*/, double qyn/*9*/, double qzn/*9*/, double qhkl/*9*/,
                                              //double ax1n, double ax2n, double ax3n,
                                              //double ax1x, double ax1y, double ax1z,
                                              //double ax2x, double ax2y, double ax2z,
                                              //double ax3x, double ax3y, double ax3z,
                                              //double sigx, double sigy, double sigz,
                                              //int ordis, int dim,
                                              int i0, int i1, int i2, int i3, int i4,
                                              double *carr1, double &pq, int n, int &trapzddeltac_cnt ) const
{
    int j;  //Zupq1=2045
    double x, tnm, sump, /*sumf, sumn, sums,*/ del;  //Zupq1=2046
    double /*fa, fb, fx,*/ pa, pb, px; //, na, nb, nx, sa, sb, sx;  //Zupq1=2047

    CHECKENDTHREAD_RET;

    if ( n==1 )
    {/*2*/  //Zupq1=2050
        switch ( i0 )
        {
        case 1:   /*  delta/chi integration  */  //Zupq1=2051
#ifdef nichtVerwendet
            switch ( i2 )   // TODO: nur case 0 wird verwendet
            {
            case 5: /*  delta and chi integration  */  //Zupq1=2052
                qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,/*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,*/ordis,dim,i0,i1,i2,i3,i4,carr1,pa);  //Zupq1=2053
                qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,/*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,*/ordis,dim,i0,i1,i2,i3,i4,carr1,pb);  //Zupq1=2054
                pa = pa/(2.0*M_PI);  //Zupq1=2055
                pb = pb/(2.0*M_PI);  //Zupq1=2056
                break;  //Zupq1=2057
            default: /*  just delta integration  */  //Zupq1=2058
                pa = pow(cos(a),i3)*pow(sin(a),i4);  //Zupq1=2059
                pb = pow(cos(b),i3)*pow(sin(b),i4);  //Zupq1=2060
                break;  //Zupq1=2061
            case 0: /*  Gauss  */  //Zupq1=2062
#endif
                pa = sin(a)*exp(-a*a/(dbeta*dbeta)); // *pa;  //Zupq1=2063
                pb = sin(b)*exp(-b*b/(dbeta*dbeta)); // *pb;  //Zupq1=2064
#ifdef nichtVerwendet
                break;  //Zupq1=2065
            case 1: /*  Exponential  */  //Zupq1=2066
                pa = sin(a)*exp(-a/dbeta)*pa;  //Zupq1=2067
                pb = sin(b)*exp(-b/dbeta)*pb;  //Zupq1=2068
                break;  //Zupq1=2069
            case 2: /*  Onsager  */  //Zupq1=2070
                pa = sin(a)*exp(-sin(a)/dbeta)*pa;  //Zupq1=2071
                pb = sin(b)*exp(-sin(b)/dbeta)*pb;  //Zupq1=2072
                break;  //Zupq1=2073
            }  //Zupq1=2074
#endif
            break;
        case 2:   /*  norm  */  //Zupq1=2075
#ifdef nichtVerwendet
            switch ( i2 )   // TODO: nur case 0  wird verwendet
            {
            case 0: /*  Gauss  */  //Zupq1=2078
#endif
                pa = sin(a)*exp(-a*a/(dbeta*dbeta));  //Zupq1=2079
                pb = sin(b)*exp(-b*b/(dbeta*dbeta));  //Zupq1=2080
#ifdef nichtVerwendet
                break;  //Zupq1=2081
            case 1: /*  Exponential  */  //Zupq1=2082
                pa = sin(a)*exp(-a/dbeta*dbeta);  //Zupq1=2083
                pb = sin(b)*exp(-b/dbeta*dbeta);  //Zupq1=2084
                break;  //Zupq1=2085
            case 2: /*  Onsager  */  //Zupq1=2086
                pa = sin(a)*exp(-sin(a)/dbeta*dbeta);  //Zupq1=2087
                pb = sin(b)*exp(-sin(b)/dbeta*dbeta);  //Zupq1=2088
                break;  //Zupq1=2089
            }  //Zupq1=2092
#endif
            break;
        case 3:   /*  order parameter  */  //Zupq1=2093
#ifdef nichtVerwendet
            switch ( i2 )   // TODO: nur case 0 wird verwendet
            {
            case 0: /*  Gauss  */  //Zupq1=2096
#endif
                pa = sin(a)*exp(-a*a/(dbeta*dbeta))*(3*cos(a)*cos(a)-1)/2.0;  //Zupq1=2097
                pb = sin(b)*exp(-b*b/(dbeta*dbeta))*(3*cos(b)*cos(b)-1)/2.0;  //Zupq1=2098
#ifdef nichtVerwendet
                break;;  //Zupq1=2099
            case 1: /*  Exponential  */  //Zupq1=2100
                pa = sin(a)*exp(-a/dbeta)*(3*cos(a)*cos(a)-1)/2.0;  //Zupq1=2101
                pb = sin(b)*exp(-b/dbeta)*(3*cos(b)*cos(b)-1)/2.0;  //Zupq1=2102
                break;  //Zupq1=2103
            case 2: /*  Onsager  */  //Zupq1=2104
                pa = sin(a)*exp(-sin(a)/dbeta)*(3*cos(a)*cos(a)-1)/2.0;  //Zupq1=2105
                pb = sin(b)*exp(-sin(b)/dbeta)*(3*cos(b)*cos(b)-1)/2.0;  //Zupq1=2106
                break;  //Zupq1=2107
            }  //Zupq1=2110
#endif
            break;
        case 4:   /*  cylinder formfactor  */  //Zupq1=2111
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pa);  //Zupq1=2112
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pb);  //Zupq1=2113
#ifdef nichtVerwendet
            switch ( i2 )   // TODO: nur case 0 wird verwendet
            {
            case 0: /*  Gauss  */  //Zupq1=2114
#endif
                pa = sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2.0*M_PI);  //Zupq1=2115
                pb = sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2.0*M_PI);  //Zupq1=2116
#ifdef nichtVerwendet
                break;  //Zupq1=2117
            case 1: /*  Exponential  */  //Zupq1=2118
                pa = sin(a)*exp(-a/dbeta)*pa/(2.0*M_PI);  //Zupq1=2119
                pb = sin(b)*exp(-b/dbeta)*pb/(2.0*M_PI);  //Zupq1=2120
                break;  //Zupq1=2121
            case 2: /*  Onsager  */  //Zupq1=2122
                pa = sin(a)*exp(-sin(a)/dbeta)*pa/(2.0*M_PI);  //Zupq1=2123
                pb = sin(b)*exp(-sin(b)/dbeta)*pb/(2.0*M_PI);  //Zupq1=2124
                break;  //Zupq1=2125
            }  //Zupq1=2126
#endif
            break;
        case 5:   /*  unit cell rotation  */  //Zupq1=2127
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pa);  //Zupq1=2128
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pb);  //Zupq1=2129
#ifdef nichtVerwendet
            switch ( i2 )   // TODO: nur case 0 wird verwendet
            {
            case 0: /*  Gauss  */  //Zupq1=2130
#endif
                pa = 8*sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  //Zupq1=2131
                pb = 8*sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  //Zupq1=2132
#ifdef nichtVerwendet
                break;  //Zupq1=2133
            case 1: /*  Exponential  */  //Zupq1=2134
                pa = 8*sin(a)*exp(-a/dbeta)*pa/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  //Zupq1=2135
                pb = 8*sin(b)*exp(-b/dbeta)*pb/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  //Zupq1=2136
                break;  //Zupq1=2137
            case 2: /*  Onsager  */  //Zupq1=2138
                pa = 8*sin(a)*exp(-sin(a)/dbeta)*pa/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  //Zupq1=2139
                pb = 8*sin(b)*exp(-sin(b)/dbeta)*pb/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  //Zupq1=2140
                break;  //Zupq1=2141
            }  //Zupq1=2142
#endif
            break;
        case 6:   /*  disk formfactor  */  //Zupq1=2143
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pa);  //Zupq1=2144
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pb);  //Zupq1=2145
#ifdef nichtVerwendet
            switch ( i2 )   // TODO: nur case 0 wird verwendet
            {
            case 0: /*  Gauss  */  //Zupq1=2146
#endif
                pa = sin(a)*exp(-a*a/(dbeta*dbeta))*pa/(2.0*M_PI);  //Zupq1=2147
                pb = sin(b)*exp(-b*b/(dbeta*dbeta))*pb/(2.0*M_PI);  //Zupq1=2148
#ifdef nichtVerwendet
                break;  //Zupq1=2149
            case 1: /*  Exponential  */  //Zupq1=2150
                pa = sin(a)*exp(-a/dbeta)*pa/(2.0*M_PI);  //Zupq1=2151
                pb = sin(b)*exp(-b/dbeta)*pb/(2.0*M_PI);  //Zupq1=2152
                break;  //Zupq1=2153
            case 2: /*  Onsager  */  //Zupq1=2154
                pa = sin(a)*exp(-sin(a)/dbeta)*pa/(2.0*M_PI);  //Zupq1=2155
                pb = sin(b)*exp(-sin(b)/dbeta)*pb/(2.0*M_PI);  //Zupq1=2156
                break;  //Zupq1=2157
                /* pa:=sin(a)*exp(-(a-pi/2)*(a-pi/2)/(dbeta*dbeta))*pa/(2*pi);  //Zupq1=2158 */
                /* pb:=sin(b)*exp(-(b-pi/2)*(b-pi/2)/(dbeta*dbeta))*pb/(2*pi);  //Zupq1=2159 */
            }  //Zupq1=2160
#endif
            break;
        case 7:    /*  cube-, triaxial ellipsoid-integration  */  //Zupq1=2161
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pa);  //Zupq1=2162
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pb);  //Zupq1=2163
            pa = pa*sin(a);  //Zupq1=2164
            pb = pb*sin(b);  //Zupq1=2165
            break;  //Zupq1=2166
        case 8:   /*  superball integration  */  //Zupq1=2167
            qrombchid(l,r,p1,sigma,alfa,a,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pa);  //Zupq1=2168
            qrombchid(l,r,p1,sigma,alfa,b,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                      /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,pb);  //Zupq1=2169
            /* pa:=1.05;  //Zupq1=2170 */
            /* pb:=1.0;  //Zupq1=2171 */
            break;  //Zupq1=2172
        } // switch i0
        pq = 0.5*(b-a)*(pa+pb);  //Zupq1=2173
        trapzddeltac_cnt = 1;  //Zupq1=2174
    }

    else    // ==========================================

    {/*2*/  //Zupq1=2176
        tnm = trapzddeltac_cnt;  //Zupq1=2177
        del = (b-a)/tnm;  //Zupq1=2178
        x = a+0.5*del;  //Zupq1=2179
        sump = 0.0;  //Zupq1=2180
        for ( j=1; j<=trapzddeltac_cnt; j++ )
        {/*3*/  //Zupq1=2181
            CHECKENDTHREAD_RET;
            switch ( i0 )
            {
            case 1:  /*  delta integration  */  //Zupq1=2182
#ifdef nichtVerwendet
                switch ( i2 )
                {
                case 5:  /*  delta and chi integration  */  //Zupq1=2184
                    qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,/*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,*/ordis,dim,i0,i1,i2,i3,i4,carr1,px);  //Zupq1=2185
                    px = px/(2.0*M_PI);  //Zupq1=2186
                    break;  //Zupq1=2187
                    /*  just delta integration  */  //Zupq1=2188
                default:
                    px = pow(cos(x),i3)*pow(sin(x),i4);  //Zupq1=2189
                    break;
                case 0:
#endif
                    px = sin(x)*exp(-x*x/(dbeta*dbeta)); // *px;  /*  Gauss  */  //Zupq1=2190
#ifdef nichtVerwendet
                    break;
                case 1:
                    px = sin(x)*exp(-x/dbeta)*px;  /*  Exponential  */  //Zupq1=2191
                    break;
                case 2:
                    px = sin(x)*exp(-sin(x)/dbeta)*px;  /*  Onsager  */  //Zupq1=2192
                    break;
                }/*4*/  //Zupq1=2193
#endif
                break;

            case 2:  /*  norm  */  //Zupq1=2195
#ifdef nichtVerwendet
                if ( i2==0 )
#endif
                    px = sin(x)*exp(-x*x/(dbeta*dbeta));  /*  Gauss  */  //Zupq1=2197
#ifdef nichtVerwendet
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta);  /*  Exponential  */  //Zupq1=2198
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta);  /*  Onsager  */  //Zupq1=2199
#endif
                break;  //Zupq1=2200

            case 3:  /*  order parameter  */  //Zupq1=2202
#ifdef nichtVerwendet
                if ( i2==0 )
#endif
                    px = sin(x)*exp(-x*x/(dbeta*dbeta))*(3*cos(x)*cos(x)-1)/2.0;  /*  Gauss  */  //Zupq1=2204
#ifdef nichtVerwendet
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta)*(3*cos(x)*cos(x)-1)/2.0;  /*  Exponential  */  //Zupq1=2205
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta)*(3*cos(x)*cos(x)-1)/2.0;  /*  Onsager  */  //Zupq1=2206
#endif
                break;  //Zupq1=2207
                /* if i0=3 then px:=sin(x)*exp(-(x-pi/2)*(x-pi/2)/(dbeta*dbeta))*(3*cos(x)*cos(x)-1)/2;  //Zupq1=2208 */

            case 4:  /*  cylinder formfactor  */  //Zupq1=2210
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,px);  //Zupq1=2211
#ifdef nichtVerwendet
                if ( i2==0 )
#endif
                    px = sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2.0*M_PI);   /*  Gauss  */  //Zupq1=2212
#ifdef nichtVerwendet
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta)*px/(2.0*M_PI);   /*  Exponential  */  //Zupq1=2213
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta)*px/(2.0*M_PI);   /*  Onsager  */  //Zupq1=2214
#endif
                break;  //Zupq1=2215

            case 5:  /*  unit cell rotation  */  //Zupq1=2216
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,px);  //Zupq1=2217
#ifdef nichtVerwendet
                if ( i2==0 )
#endif
                    px = 8*sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  /*  Gauss  */  //Zupq1=2218
#ifdef nichtVerwendet
                if ( i2==1 ) px = 8*sin(x)*exp(-x/dbeta)*px/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  /*  Exponential  */  //Zupq1=2219
                if ( i2==2 ) px = 8*sin(x)*exp(-sin(x)/dbeta)*px/(2.0*M_PI*pow(M_PI,3/2.0)*params.sig.x()*params.sig.y()*params.sig.z());  /*  Onsager  */  //Zupq1=2220
#endif
                break;  //Zupq1=2221

            case 6:  /*  disk formfactor  */  //Zupq1=2222
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,px);  //Zupq1=2223
#ifdef nichtVerwendet
                if ( i2==0 )
#endif
                    px = sin(x)*exp(-x*x/(dbeta*dbeta))*px/(2.0*M_PI);  /*  Gauss  */  //Zupq1=2224
#ifdef nichtVerwendet
                if ( i2==1 ) px = sin(x)*exp(-x/dbeta)*px/(2.0*M_PI);  /*  Exponential  */  //Zupq1=2225
                if ( i2==2 ) px = sin(x)*exp(-sin(x)/dbeta)*px/(2.0*M_PI);  /*  Onsager  */  //Zupq1=2226
#endif
                /* px:=sin(x)*exp(-(x-pi/2)*(x-pi/2)/(dbeta*dbeta))*px/(2*pi);  //Zupq1=2227 */
                break;  //Zupq1=2228

            case 7:  /*  cube-, triaxial ellipsoid-integration  */  //Zupq1=2229
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,px);  //Zupq1=2230
                px = px*sin(x);  //Zupq1=2231
                break;  //Zupq1=2232

            case 8:  /*  superball integration  */  //Zupq1=2233
                qrombchid(l,r,p1,sigma,alfa,x,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,qxn,qyn,qzn,qhkl,
                          /*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,ordis,dim,i0,*/i1,i2,i3,i4,carr1,px);  //Zupq1=2234
                /* px:=0.99;  //Zupq1=2235 */
                break;  //Zupq1=2236
            } // switch i0
            sump = sump+px;  //Zupq1=2237
            x = x+del;  //Zupq1=2238
        }/* for j */  //Zupq1=2239
        pq = 0.5*(pq+(b-a)*sump/tnm);  //Zupq1=2240
        trapzddeltac_cnt = 2*trapzddeltac_cnt;  //Zupq1=2241
    }/* n==1 else */  //Zupq1=2242
}  //Zupq1=2243



//(*** interpolation routine ***)
#ifdef __CUDACC__
__host__ __device__
#endif
void SasCalc_GENERIC_calculation::polint( double *xa/*RealArrayNP*/, double *ya/*RealArrayNP*/,
                                        int n, double /*x*/, double &y, double &dy, const char * /* DSM(dbg) */ ) const
{
    // Parameter x ist bei allen Aufrufen immer 0.0 - daher kann der auch weggelassen werden!
    // Parameter xa= 1 0.25 0.0625 0.015625 0.00390625 bei allen Aufrufen in qrombdeltac() und qrombchid()
    int ns,m,i;
    double w,hp,ho,dift,dif,den;
    RealArrayNP c, d; // ^RealArrayNP;     =array[1..np=5] of extended;

    //D8L( qDebug() << "   polint Anf xa=" <<xa[1]<<xa[2]<<xa[3]<<xa[4]<<xa[5] << "ya=" <<ya[1]<<ya[2]<<ya[3]<<ya[4]<<ya[5] << dbg );
    ns=1;
    dif=fabs(xa[1]);           // fabs(x-xa[1])
    for ( i=1/*?? 2*/; i<=n; i++ )
    {
        dift=fabs(xa[i]);      // fabs(x-xa[i])
        if ( dift < dif )
        {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    //D8L( qDebug() << "   polint Inner" << dif << ns << ya[ns] );
    y=ya[ns];
    ns=ns-1;
    for ( m=1; m<=n-1; m++ )
    {
        for ( i=1; i<=n-m; i++ )
        {
            ho=xa[i];       // -x
            hp=xa[i+m];     // -x
            w=c[i+1]-d[i];
            den=ho-hp;
            //if ( fabs(den) > eps5 ) den=w/den;
            if ( den != 0 ) den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        } // for i
        if ( 2*ns < n-m )
            dy=c[ns+1];
        else
        {
            dy=d[ns];
            ns=ns-1;
        }
        y=y+dy;
        //qDebug() << m << "c=" << arr2str(c,n) << "d=" << arr2str(d,n) << y << dy;
        //D8L( qDebug() << "   polint m" << m << y << dy );
    } // for m

    DSM( if ( fabs(dy) < 1.0e-200 && fabs(dy) > 0 )
            qDebug() << "polint ya=" <<ya[1]<<ya[2]<<ya[3]<<ya[4]<<ya[5] << dbg
            << y << dy; )

}




#ifdef __CUDACC__
__host__ __device__
#endif
//void SasCalc_GENERIC_calculation::qrombdeltac(
//    double l,     = params.length   -> raus
//    double r,     = params.radius   -> raus
//    double p1,    = params.p1 | params.radiusi(nur einmal)   ==> bleibt
//    double sigma, = params.sigma | sigmal(formpq-Parameter) | params.sigmal   ==> bleibt
//    double alfa,  = params.alphash | params.alpha   ==> bleibt
//    double dbeta, = params.dbeta   -> raus
//    double theta, = params.polTheta | theta(loc.Var)   ==> bleibt
//    double phi,   = params.polPhi | phi(loc.Var)   ==> bleibt
//    double qx,    = qx(loc.Var) | qxs(Param) | 1   ==> bleibt
//    double qy,    = qy(loc.Var) | qys(Param) | 1   ==> bleibt
//    double qz,    = qz(loc.Var) | 1   ==> bleibt
//    double qxn,   = 9 | qxhkl(loc.Var)   ==> bleibt
//    double qyn,   = 9 | qyhkl(loc.Var)   ==> bleibt
//    double qzn,   = 9 | qzhkl(loc.Var)   ==> bleibt
//    double qhkl,  = 9 | qhkl(loc.Var)   ==> bleibt
//    int ordis,    = ordis(Param in formpq) | params.ordis   ==> bleibt
//    int dim,      = 2 | 3 u.a.   ==> bleibt
//    int i0,       = 6 | 5 u.a.   ==> bleibt
//    int i1,       = params.orcase+7 | 6 u.a.   ==> bleibt
//    int i2,       = 0 u.a.   ==> bleibt
//    int i3,       = 0 u.a.   ==> bleibt
//    int i4,       = 0 u.a.   ==> bleibt
//    double *carr1,= params.CR->carr2p | params.CR->carr1p   ==> bleibt
//    double &pq    = (Returnwert)
//   ) const

void SasCalc_GENERIC_calculation::qrombdeltac( double p1, double sigma, double alfa,
                                             double theta, double phi, double qx/*1*/, double qy/*1*/, double qz/*1*/,
                                             double qxn/*9*/, double qyn/*9*/, double qzn/*9*/, double qhkl/*9*/,
                                             //double ax1n, double ax2n, double ax3n, double ax1x, double ax1y, double ax1z,
                                             //double ax2x, double ax2y, double ax2z, double ax3x, double ax3y, double ax3z,
                                             //double sigx, double sigy, double sigz,
                                             //int ordis, int dim,
                                             int i0, int i1, int i2, int i3, int i4,
                                             double *carr1, double &pq ) const
{
    /* label 99; */  //Zupq1=1968
    const double eps = 1.0e-4;  //Zupq1=1970
    //const double eps1 = 1.0e-8;  //Zupq1=1971
    const int    jmax = 20;  //Zupq1=1972
    //const int    jmaxp = 21;  //Zupq1=1973 global
    const int    k = 5;  //Zupq1=1974
    //const int    np = 5;  //Zupq1=1975 global
    const double prmin =  0.00001;       /*  minimum pr(delta)-value  */  //Zupq1=1976

    //double ax1n, ax2n, ax3n, ax1x, ax1y, ax1z, ax2x, ax2y, ax2z, ax3x, ax3y, ax3z, sigx, sigy, sigz;

    DTL( static time_t startzeit = 0;
        static long   dbgCntDiff = 0 );

    //typedef double RealArrayJMAXP[jmaxp+1];  //Zupq1=1978 global
    //typedef double RealArrayNP[np+1];  //Zupq1=1979 global

    int i, j;  //Zupq1=1981
    //int trapzdit, midpntit;  //Zupq1=1982
    double dssp, delmax=0; //, norm;  //Zupq1=1983
    RealArrayJMAXP hp, sp;  //Zupq1=1984
    RealArrayNP    cp, dp;  //Zupq1=1985
    float alim, blim, p11, p12, p13, p21, p22, p23, p31, p32, p33;       //Zupq1=1986

    int trapzddeltac_cnt;

    DTL( if ( startzeit == 0 ) startzeit = time(nullptr) );

    //D8L( qDebug() << "qrombdeltac(" << r << theta << phi << "qx/y/z" << qx << qy << qz << "q?n" << qxn << qyn << qzn << qhkl << "ax?n" << ax1n << ax2n << ax3n
    //                               << "ax1?" << ax1x << ax1y << ax1z << "ax2?" << ax2x << ax2y << ax2z << "ax3?" << ax3x << ax3y << ax3z
    //                               << "sig?" << sigx << sigy << sigz << "ordis" << ordis << dim << "i?" << i0 << i1 << i2 << i3 << i4 << "&pq )" << params.dbeta );

    //memset( hp, 0, sizeof(hp) );
    //memset( sp, 0, sizeof(sp) );
    //memset( cp, 0, sizeof(cp) );
    //memset( dp, 0, sizeof(dp) );

    // params.dbeta ist die globale Variable in Grad.
    double dbeta = params.dbeta*M_PI/180.0;  //Zupq1=2247
    theta = theta*M_PI/180.0;  //Zupq1=2248
    phi = phi*M_PI/180.0;  //Zupq1=2249

    /*  search for maximum integration angle  */  //Zupq1=2250
    // ==> i2=0 or i2=7, nothing else!
    if ( (i2==0) /*|| (i2==2) || (i2==3) || (i2==5)*/ ) delmax = dbeta*sqrt(log(1.0/prmin));  /*  Gaussian, Onsager, Maier-Saupe, Laguerre  */  //Zupq1=2251
    //if ( i2==1 ) delmax = dbeta*log(1.0/prmin);        /*  Exponential  */  //Zupq1=2252
    /* if i2=2 then delmax:=arcsin(dbeta*ln(1/prmin));  //Zupq1=2253 */
    //if ( i2==4 ) delmax = dbeta;                    /*  cut-off  */  //Zupq1=2254
    if ( (i2==7) /*|| (i2==8) || (i2==9) || (i2==10) || (i2==11) || (i2==12)*/ ) delmax = M_PI/2.0; /*  isotropic, mirrored distributions  */  //Zupq1=2255

    if ( delmax>M_PI/2.0 ) delmax = M_PI/2.0;  //Zupq1=2257
    alim = 0.0;  //Zupq1=2258
    blim = delmax;  //Zupq1=2259

    if ( i1==12 )
    {   /*  for cubes  */  //Zupq1=2261
        alim = eps;  //Zupq1=2262
        blim = M_PI/2.0-eps;  //Zupq1=2263
    }  //Zupq1=2264

    if ( i1==17 )
    {   /*  for superballs  */  //Zupq1=2266
        alim = 0.0001;  //Zupq1=2267
        /* blim:=r-0.0001;  //Zupq1=2268 */
        blim = M_PI/2.0-0.0001;  //Zupq1=2269
    }  //Zupq1=2270

    /*  for disks: mirrored Gaussian  */  //Zupq1=2273
    /* if ((qzn=2) or (i0=6)) then begin  //Zupq1=2274 */
    /*    alim:=pi/2-delmax;  //Zupq1=2275 */
    /*    blim:=pi/2;  //Zupq1=2276 */
    /* end;  //Zupq1=2277 */

    p11 = -cos(phi)*cos(theta);  //Zupq1=2279
    p12 = sin(phi);  //Zupq1=2280
    p13 = cos(phi)*sin(theta);  //Zupq1=2281
    p21 = -cos(phi);  //Zupq1=2282
    p22 = -sin(phi)*cos(theta);  //Zupq1=2283
    p23 = sin(phi)*sin(theta);  //Zupq1=2284
    p31 = -sin(theta);  //Zupq1=2285
    p32 = 0;  //Zupq1=2286
    p33 = -cos(theta);  //Zupq1=2287

    hp[1] = 1.0;  //Zupq1=2293
    for ( j=1; j<=jmax; j++ )
    {  //Zupq1=2294
        CHECKENDTHREAD_RET;
        trapzddeltac(alim,blim,params.length,params.radius,p1,sigma,alfa,dbeta,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,p31,p32,p33,
                     qxn,qyn,qzn,qhkl,/*ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,*/
                     /*ordis,dim,*/i0,i1,i2,i3,i4,carr1,sp[j],j, trapzddeltac_cnt );  //Zupq1=2295
        if ( j>=k )
        {  //Zupq1=2296
            for ( i=1; i<=k; i++ )
            {  //Zupq1=2297
                cp[i] = hp[j-k+i];  //Zupq1=2298
                dp[i] = sp[j-k+i];  //Zupq1=2299
            }  //Zupq1=2300
            polint(cp,dp,k,0.0,pq,dssp,"qrombdeltac");  //Zupq1=2301
            D8L( qDebug() << "qrombdeltac inner" << j << dssp << pq );
            if ( abs(dssp)<(eps*abs(pq)) ) break; /* goto 99 */  //Zupq1=2302
            if ( fabs(pq) < 1.0e-100 ) break;  // nicht im Pascal-Programm
        }  //Zupq1=2303
        sp[j+1] = sp[j];  //Zupq1=2304
        hp[j+1] = 0.25*hp[j];  //Zupq1=2305
    }  //Zupq1=2306

    D8L( qDebug() << "qrombdeltac: trapzddeltac_cnt =" << trapzddeltac_cnt );

}



#ifdef __CUDACC__
__host__ __device__
#endif
void SasCalc_GENERIC_calculation::qrombchid( double l, double r, double p1, double sigma, double alfa, double delta,
                                           double theta, double phi,
                                           double qx, double qy, double qz,
                                           double p11, double p12, double p13, double p21, double p22, double p23,
                                           double p31, double p32, double p33,
                                           double qxn, double qyn, double qzn, double qhkl,
                                           //double ax1n, double ax2n, double ax3n,
                                           //double ax1x, double ax1y, double ax1z,
                                           //double ax2x, double ax2y, double ax2z,
                                           //double ax3x, double ax3y, double ax3z,
                                           //double sigx, double sigy, double sigz,
                                           //int ordis, int dim, int i0,
                                           int i1, int i2, int i3, int i4,
                                           double *carr1, double &pq ) const
{/*1*/  //Zupq1=53

    const double eps = 1.0e-4;  //Zupq1=57
    const int    jmax = 20;  //Zupq1=60
    const int    k = 5;  //Zupq1=62

    double dssp;
    RealArrayJMAXP hp, sp;  //Zupq1=72
    RealArrayNP    cp, dp;  //Zupq1=73
    float alim, blim;       //Zupq1=74

    int trapzdchid_cnt;

    //qDebug() << "qrombchid(" << r << sigma << dbeta << delta << theta << phi << qx << qy << qz << p11 << p12 << p13 << p21
    //         << p22 << p23 << p31 << p32 << p33 << qxn << qyn << qzn << qhkl << ax1n << ax2n << ax3n << ax1x << ax1y
    //         << ax1z << ax2x << ax2y << ax2z << ax3x << ax3y << ax3z << sigx << sigy << sigz << "ordis" << ordis << dim
    //         << i0 << i1 << i2 << i3 << i4 << "&pq" << " );";

    if ( i1==12 )
    {   /*  for cubes  */  //Zupq1=1424
        alim = eps;  //Zupq1=1425
        blim = M_PI/2.0-eps;  //Zupq1=1426
    }   //Zupq1=1427
    else if ( (i1==13) || (i1==14) )
    {   /*  for ellipsoids  */  //Zupq1=1429
        alim = 0;  //Zupq1=1430
        blim = M_PI/2.0;  //Zupq1=1431
    }   //Zupq1=1432
    else if ( i1==15 )
    {   /*  for area integration  */  //Zupq1=1434
        alim = 0.0;  //Zupq1=1435
        blim = l;  //Zupq1=1436
    }   //Zupq1=1437
    else if ( i1==17 )
    {   /*  superball integration  */  //Zupq1=1439
        alim = 0;  //Zupq1=1440
        /* blim:=r*power(1-power(delta/r,p1),1/p1)-0.0005;  //Zupq1=1441 */
        blim = M_PI/2.0;  //Zupq1=1442
    }   //Zupq1=1443
    else
    {   // Default
        alim = 0.0;  //Zupq1=1421
        blim = 2*M_PI;  //Zupq1=1422
    }

    hp[1] = 1.0;  //Zupq1=1449
    for ( int j=1; j<=jmax; j++ )
    {  //Zupq1=1450
        CHECKENDTHREAD_RET;
        /* trapzdchid(alim,blim,l,r,p1,rho,alfa,sigmal,sigma,cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,dbeta,phi,theta,qxn,qyn,qzn,q,maxit,i0,i1,i2,i3,i4,sp^[j],j);  //Zupq1=1451 */
        trapzdchid(alim,blim,l,r,p1,sigma,alfa,delta,theta,phi,qx,qy,qz,p11,p12,p13,p21,p22,p23,
                   p31,p32,p33,qxn,qyn,qzn,qhkl,
                   //ax1n,ax2n,ax3n,ax1x,ax1y,ax1z,ax2x,ax2y,ax2z,ax3x,ax3y,ax3z,sigx,sigy,sigz,
                   //ordis,dim,i0,
                   i1,i2,i3,i4,carr1,sp[j],j, trapzdchid_cnt );  //Zupq1=1452
        if ( j>=k )
        {  //Zupq1=1453
            for ( int i=1; i<=k; i++ )
            {  //Zupq1=1454
                cp[i] = hp[j-k+i];  //Zupq1=1455
                dp[i] = sp[j-k+i];  //Zupq1=1456
            }  //Zupq1=1457
            polint(cp,dp,k,0.0,pq,dssp,"qrombchid");  //Zupq1=1458
            if ( fabs(dssp)<(eps*fabs(pq)) ) break; /* goto 99*/  //Zupq1=1459
            if ( fabs(pq) < 1.0e-300 )
            {
                //qDebug() << " ... qrombchid" << dssp << pq;
                break; // TODO ?
            }
        }  //Zupq1=1460
        sp[j+1] = sp[j];  //Zupq1=1461
        hp[j+1] = 0.25*hp[j];  //Zupq1=1462
    }  //Zupq1=1463
    /*99:*/  //Zupq1=1464
} /* qrombchid() */



#ifdef __CUDACC__
__host__ __device__
#endif
void SasCalc_GENERIC_calculation::trapzdchid( double a, double b, double l, double r, double p1, double sigma, double alfa,
                                            double delta, double theta, double phi,
                                            double qx, double qy, double qz,
                                            double p11, double p12, double p13, double p21, double p22, double p23,
                                            double p31, double /*p32*/, double p33,
                                            double qxn, double qyn, double qzn, double qhkl,
                                            //double ax1n, double ax2n, double ax3n,
                                            //double ax1x, double ax1y, double ax1z,
                                            //double ax2x, double ax2y, double ax2z,
                                            //double ax3x, double ax3y, double ax3z,
                                            //double sigx, double sigy, double sigz,
                                            //int /*ordis*/, int /*dim*/, int /*i0*/,
                                            int i1, int /*i2*/, int i3, int /*i4*/,
                                            double *carr1, double &pq, int n,
                                            int &trapzdchid_cnt ) const
{
    // Aufgerufen von qrombchid()

    const double eps =  0.0000000001;  //Zupq1=133
    /*
    ax1n = params.ax1.length();
    ax2n = params.ax2.length();
    ax3n = params.ax3.length();
    ax1x = params.ax1.x();
    ax1y = params.ax1.y();
    ax1z = params.ax1.z();
    ax2x = params.ax2.x();
    ax2y = params.ax2.y();
    ax2z = params.ax2.z();
    ax3x = params.ax3.x();
    ax3y = params.ax3.y();
    ax3z = params.ax3.z();
    sigx = params.sig.x();
    sigy = params.sig.y();
    sigz = params.sig.z();
    */

    CHECKENDTHREAD_RET;

    // Die Randbedingungen der weiteren Abfragen steht nur im oberen Teil (n=1) ist aber im unteren Teil identisch
    if ( n==1 )
    {/*2*/  //Zupq1=163
        double pa=0, pb=0;
        if ( i1==1 )    // nur bei orcase==1, immer i3=0
        {/*3*/    /*  cylinder, general case  */  //Zupq1=164
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=165
            //const double a1 = 1/(2.0*z*(z-1));  //Zupq1=166
            const double mla1 = p11*cos(a)*sin(delta)+p12*sin(a)*sin(delta)+p13*cos(delta);  //Zupq1=167
            const double mla2 = p21*sin(a)*sin(delta)+p22*cos(a)*sin(delta)+p23*cos(delta);  //Zupq1=168
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta);  //Zupq1=169 */
            const double mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=170
            const double mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=171
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta);  //Zupq1=172 */
            const double arga = (qx*mla1+qy*mla2+eps)*l/(z+1);  //Zupq1=173
            const double argb = (qx*mlb1+qy*mlb2+eps)*l/(z+1);  //Zupq1=174
#ifdef nichtVerwendet
            if ( i3==0 )
#endif
            {/*4*/ /*  P(q)  */  //Zupq1=175
                const double a1 = 1/(2.0*z*(z-1));  //Zupq1=176
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=177
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=178
            }/*4*/  //Zupq1=179
#ifdef nichtVerwendet
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=180
                const double pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=181
                pa = pa1*pa1;  //Zupq1=182
                const double pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=183
                pb = pb1*pb1;  //Zupq1=184
            }/*4*/  //Zupq1=185
#endif
        }/*3*/  //Zupq1=186

        if ( i1==2 )    // nur bei orcase==2, immer i3=0
        {/*3*/   /*  cylinder, x-axis  */  //Zupq1=187
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=188
            const double arga = (qx*cos(delta)-qy*sin(a)*sin(delta)+eps)*l/(z+1);  //Zupq1=189
            const double argb = (qx*cos(delta)-qy*sin(b)*sin(delta)+eps)*l/(z+1);  //Zupq1=190
#ifdef nichtVerwendet
            if ( i3==0 )
#endif
            {/*4*/ /*  P(q)  */  //Zupq1=191
                const double a1 = 1/(2.0*z*(z-1));  //Zupq1=192
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=193
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=194
            }/*4*/  //Zupq1=195
#ifdef nichtVerwendet
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=196
                const double pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=197
                pa = pa1*pa1;  //Zupq1=198
                const double pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=199
                pb = pb1*pb1;  //Zupq1=200
            }/*4*/  //Zupq1=201
#endif
        }/*3*/  //Zupq1=202

        if ( i1==3 )    // nur bei orcase==3, immer i3=0
        {/*3*/   /*  cylinder, y-axis  */  //Zupq1=203
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=204
            const double arga = (qx*sin(a)*sin(delta)+qy*cos(delta)+eps)*l/(z+1);  //Zupq1=205
            const double argb = (qx*sin(b)*sin(delta)+qy*cos(delta)+eps)*l/(z+1);  //Zupq1=206
#ifdef nichtVerwendet
            if ( i3==0 )
#endif
            {/*4*/ /*  P(q)  */  //Zupq1=207
                const double a1 = 1/(2.0*z*(z-1));  //Zupq1=208
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=209
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=210
            }/*4*/  //Zupq1=211
#ifdef nichtVerwendet
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=212
                const double pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=213
                pa = pa1*pa1;  //Zupq1=214
                const double pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=215
                pb = pb1*pb1;  //Zupq1=216
            }/*4*/  //Zupq1=217
#endif
        }/*3*/  //Zupq1=218

        if ( i1==4 )    // nur bei orcase==4, immer i3=0
        {/*3*/   /*  cylinder, -z-axis  */  //Zupq1=219
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=220
            const double arga = (qx*sin(a)*sin(delta)-qy*cos(a)*sin(delta)+eps)*l/(z+1);  //Zupq1=221
            const double argb = (qx*sin(b)*sin(delta)-qy*cos(b)*sin(delta)+eps)*l/(z+1);  //Zupq1=222
#ifdef nichtVerwendet
            if ( i3==0 )
#endif
            {/*4*/ /*  P(q)  */  //Zupq1=223
                const double a1 = 1/(2.0*z*(z-1));  //Zupq1=224
                pa = (a1/(arga*arga))*(1-cos((z-1)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-1)/2.0));  //Zupq1=225
                pb = (a1/(argb*argb))*(1-cos((z-1)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-1)/2.0));  //Zupq1=226
            }/*4*/  //Zupq1=227
#ifdef nichtVerwendet
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=228
                const double pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=229
                pa = pa1*pa1;  //Zupq1=230
                const double pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=231
                pb = pb1*pb1;  //Zupq1=232
            }/*4*/  //Zupq1=233
#endif
        }/*3*/  //Zupq1=234

#ifdef nichtVerwendet
        if ( i1==5 )
        {/*3*/    /*  general series expansion  */  //Zupq1=235
            const double mla1 = p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=236
            const double mla2 = p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=237
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta);  //Zupq1=238 */
            const double mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=239
            const double mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=240
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta);  //Zupq1=241 */
            pa = pow(mla1,i3)*pow(mla2,i4);  //Zupq1=242
            pb = pow(mlb1,i3)*pow(mlb2,i4);  //Zupq1=243
        }/*3*/  //Zupq1=244
#endif

        if ( i1==6 )
        {/*3*/    /*  unit cell rotation  */  //Zupq1=245
            const double dqxa = qx-qhkl*(p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));  //Zupq1=246
            const double dqya = qy-qhkl*(p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));  //Zupq1=247
            const double dqza = qz-qhkl*(p31*cos(a)*sin(delta)+p33*cos(delta));  //Zupq1=248
            const double dqxb = qx-qhkl*(p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta));  //Zupq1=249
            const double dqyb = qy-qhkl*(p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta));  //Zupq1=250
            const double dqzb = qz-qhkl*(p31*cos(b)*sin(delta)+p33*cos(delta));  //Zupq1=251
            const double dqs1a = (dqxa*params.ax1.x()+dqya*params.ax1.y()+dqza*params.ax1.z())/(params.ax1.length()/*ax1n*/*params.sig.x());;  //Zupq1=252
            const double dqs2a = (dqxa*params.ax2.x()+dqya*params.ax2.y()+dqza*params.ax2.z())/(params.ax2.length()/*ax2n*/*params.sig.y());  //Zupq1=253
            const double dqs3a = (dqxa*params.ax3.x()+dqya*params.ax3.y()+dqza*params.ax3.z())/(params.ax3.length()/*ax3n*/*params.sig.z());  //Zupq1=254
            const double dqs1b = (dqxb*params.ax1.x()+dqyb*params.ax1.y()+dqzb*params.ax1.z())/(params.ax1.length()/*ax1n*/*params.sig.x());  //Zupq1=255
            const double dqs2b = (dqxb*params.ax2.x()+dqyb*params.ax2.y()+dqzb*params.ax2.z())/(params.ax2.length()/*ax2n*/*params.sig.y());  //Zupq1=256
            const double dqs3b = (dqxb*params.ax3.x()+dqyb*params.ax3.y()+dqzb*params.ax3.z())/(params.ax3.length()/*ax3n*/*params.sig.z());  //Zupq1=257
            pa = exp(-4*(dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a)/M_PI);  //Zupq1=258
            pb = exp(-4*(dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b)/M_PI);  //Zupq1=259
        }/*3*/  //Zupq1=260

        if ( i1==7 )    // nur hier: i2=h, i3=k, i4=l, werden aber nicht verwendet
        {/*3*/    /*  fiber rotation  */  //Zupq1=261
            /*  rotation axis, director  */  //Zupq1=262
            const double l1 = sin(theta)*cos(phi);  //Zupq1=263
            const double l2 = sin(theta)*sin(phi);  //Zupq1=264
            const double l3 = -cos(theta);  //Zupq1=265

            /*  rotation matrix Ra  */  //Zupq1=267
            const double r11a = cos(a)+(1-cos(a))*l1*l1;  //Zupq1=268
            const double r12a = -l3*sin(a)+(1-cos(a))*l1*l2;  //Zupq1=269
            const double r13a = -l2*sin(a)+(1-cos(a))*l1*l3;  //Zupq1=270
            const double r21a = l3*sin(a)+(1-cos(a))*l1*l2;  //Zupq1=271
            const double r22a = cos(a)+(1-cos(a))*l2*l2;  //Zupq1=272
            const double r23a = l1*sin(a)+(1-cos(a))*l2*l3;  //Zupq1=273
            const double r31a = l2*sin(a)+(1-cos(a))*l1*l3;  //Zupq1=274
            const double r32a = -l1*sin(a)+(1-cos(a))*l2*l3;  //Zupq1=275
            const double r33a = cos(a)+(1-cos(a))*l3*l3;  //Zupq1=276

            /*  rotation matrix Rb  */  //Zupq1=278
            const double r11b = cos(b)+(1-cos(b))*l1*l1;  //Zupq1=279
            const double r12b = -l3*sin(b)+(1-cos(b))*l1*l2;  //Zupq1=280
            const double r13b = -l2*sin(b)+(1-cos(b))*l1*l3;  //Zupq1=281
            const double r21b = l3*sin(b)+(1-cos(b))*l1*l2;  //Zupq1=282
            const double r22b = cos(b)+(1-cos(b))*l2*l2;  //Zupq1=283
            const double r23b = l1*sin(b)+(1-cos(b))*l2*l3;  //Zupq1=284
            const double r31b = l2*sin(b)+(1-cos(b))*l1*l3;  //Zupq1=285
            const double r32b = -l1*sin(b)+(1-cos(b))*l2*l3;  //Zupq1=286
            const double r33b = cos(b)+(1-cos(b))*l3*l3;  //Zupq1=287

            /*  rotate scattering vector  */  //Zupq1=294
            const double qxhkla = r11a*qxn+r12a*qyn+r13a*qzn;  //Zupq1=295
            const double qyhkla = r21a*qxn+r22a*qyn+r23a*qzn;  //Zupq1=296
            const double qzhkla = r31a*qxn+r32a*qyn+r33a*qzn;  //Zupq1=297
            const double qxhklb = r11b*qxn+r12b*qyn+r13b*qzn;  //Zupq1=298
            const double qyhklb = r21b*qxn+r22b*qyn+r23b*qzn;  //Zupq1=299
            const double qzhklb = r31b*qxn+r32b*qyn+r33b*qzn;  //Zupq1=300

            const double dqxa = qx-qxhkla;  //Zupq1=378
            const double dqya = qy-qyhkla;  //Zupq1=379
            const double dqza = qz-qzhkla;  //Zupq1=380

            const double dqxb = qx-qxhklb;  //Zupq1=382
            const double dqyb = qy-qyhklb;  //Zupq1=383
            const double dqzb = qz-qzhklb;  //Zupq1=384

            const double ax1xa = (r11a*params.ax1.x()+r12a*params.ax1.y()+r13a*params.ax1.z());  //Zupq1=393
            const double ax1ya = (r21a*params.ax1.x()+r22a*params.ax1.y()+r23a*params.ax1.z());  //Zupq1=394
            const double ax1za = (r31a*params.ax1.x()+r32a*params.ax1.y()+r33a*params.ax1.z());  //Zupq1=395
            const double ax1xb = (r11b*params.ax1.x()+r12b*params.ax1.y()+r13b*params.ax1.z());  //Zupq1=396
            const double ax1yb = (r21b*params.ax1.x()+r22b*params.ax1.y()+r23b*params.ax1.z());  //Zupq1=397
            const double ax1zb = (r31b*params.ax1.x()+r32b*params.ax1.y()+r33b*params.ax1.z());  //Zupq1=398
            const double ax2xa = (r11a*params.ax2.x()+r12a*params.ax2.y()+r13a*params.ax2.z());  //Zupq1=399
            const double ax2ya = (r21a*params.ax2.x()+r22a*params.ax2.y()+r23a*params.ax2.z());  //Zupq1=400
            const double ax2za = (r31a*params.ax2.x()+r32a*params.ax2.y()+r33a*params.ax2.z());  //Zupq1=401
            const double ax2xb = (r11b*params.ax2.x()+r12b*params.ax2.y()+r13b*params.ax2.z());  //Zupq1=402
            const double ax2yb = (r21b*params.ax2.x()+r22b*params.ax2.y()+r23b*params.ax2.z());  //Zupq1=403
            const double ax2zb = (r31b*params.ax2.x()+r32b*params.ax2.y()+r33b*params.ax2.z());  //Zupq1=404
            const double ax3xa = (r11a*params.ax3.x()+r12a*params.ax3.y()+r13a*params.ax3.z());  //Zupq1=405
            const double ax3ya = (r21a*params.ax3.x()+r22a*params.ax3.y()+r23a*params.ax3.z());  //Zupq1=406
            const double ax3za = (r31a*params.ax3.x()+r32a*params.ax3.y()+r33a*params.ax3.z());  //Zupq1=407
            const double ax3xb = (r11b*params.ax3.x()+r12b*params.ax3.y()+r13b*params.ax3.z());  //Zupq1=408
            const double ax3yb = (r21b*params.ax3.x()+r22b*params.ax3.y()+r23b*params.ax3.z());  //Zupq1=409
            const double ax3zb = (r31b*params.ax3.x()+r32b*params.ax3.y()+r33b*params.ax3.z());  //Zupq1=410

            const double dqs1a = (dqxa*ax1xa+dqya*ax1ya+dqza*ax1za)/(params.ax1.length()/*ax1n*/*params.sig.x());  //Zupq1=412
            const double dqs2a = (dqxa*ax2xa+dqya*ax2ya+dqza*ax2za)/(params.ax2.length()/*ax2n*/*params.sig.y());  //Zupq1=413
            const double dqs3a = (dqxa*ax3xa+dqya*ax3ya+dqza*ax3za)/(params.ax3.length()/*ax3n*/*params.sig.z());  //Zupq1=414
            const double dqs1b = (dqxb*ax1xb+dqyb*ax1yb+dqzb*ax1zb)/(params.ax1.length()/*ax1n*/*params.sig.x());  //Zupq1=415
            const double dqs2b = (dqxb*ax2xb+dqyb*ax2yb+dqzb*ax2zb)/(params.ax2.length()/*ax2n*/*params.sig.y());  //Zupq1=416
            const double dqs3b = (dqxb*ax3xb+dqyb*ax3yb+dqzb*ax3zb)/(params.ax3.length()/*ax3n*/*params.sig.z());  //Zupq1=417

            const double arga = dqs1a*dqs1a+dqs2a*dqs2a+dqs3a*dqs3a;  //Zupq1=426
            const double argb = dqs1b*dqs1b+dqs2b*dqs2b+dqs3b*dqs3b;  //Zupq1=427
            pa = exp(-4*arga/M_PI);  //Zupq1=428
            pb = exp(-4*argb/M_PI);  //Zupq1=429
            /* if (arga>2) then pa:=eps  //Zupq1=430 */
            /*    else pa:=exp(-4*arga/pi);  //Zupq1=431 */
            /* if (argb>2) then pb:=eps  //Zupq1=432 */
            /*    else pb:=exp(-4*argb/pi);  //Zupq1=433 */
        }/*3*/  //Zupq1=434

        if ( i1==8 )    // nur bei orcase==1, immer i3=0
        {/*3*/    /*  disk, general case  */  //Zupq1=435
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=436
            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=437
            const double qxl = qx/qn;  //Zupq1=438
            const double qyl = qy/qn;  //Zupq1=439
            const double mla1 = p11*cos(a)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=440
            const double mla2 = p21*sin(a)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=441
            /* mla3:=p31*cos(a)*sin(delta)+p33*cos(delta);  //Zupq1=442 */
            const double mlb1 = p11*cos(b)*sin(delta)+p12*sin(b)*sin(delta)+p13*cos(delta);  //Zupq1=443
            const double mlb2 = p21*sin(b)*sin(delta)+p22*cos(b)*sin(delta)+p23*cos(delta);  //Zupq1=444
            /* mlb3:=p31*cos(b)*sin(delta)+p33*cos(delta);  //Zupq1=445 */
            const double qnna = (qxl*mla1+qyl*mla2);  //Zupq1=446
            const double qnnb = (qxl*mlb1+qyl*mlb2);  //Zupq1=447
            const double arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=448
            const double argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=449

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=451
                if ( arga<0.015 )
                {/*5*/  //Zupq1=452
                    pa = 1;  //Zupq1=453
                    double oldpa = 0;  //Zupq1=454
                    double argser = 1;  //Zupq1=455
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=456
                        argser = argser*arga*arga/4.0;  //Zupq1=457
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=458 */
                        pa = pa+carr1[i]*argser;  //Zupq1=459
                        const double delser = fabs((pa-oldpa)/pa);  //Zupq1=460
                        if ( delser<0.0001 ) break; /* goto 12; */  //Zupq1=461
                        oldpa = pa;  //Zupq1=462
                    }/*6*/  //Zupq1=463
                    /*12:*/  //Zupq1=464
#ifdef nichtVerwendet
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=465
#endif
                }/*5*/  //Zupq1=466
                else
                {/*5*/  //Zupq1=467
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/ /*  P(q)  */  //Zupq1=468
                        const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=469
                        const double pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=470
                        const double pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=471
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=472
                    }/*6*/  //Zupq1=473
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=474
                        const double pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=475
                        const double pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=476
                        const double pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=477
                        pa = pa3*pa3;  //Zupq1=478
                    }/*6*/  //Zupq1=479
#endif
                }/*5*/  //Zupq1=480
                if ( argb<0.015 )
                {/*5*/  //Zupq1=481
                    pb = 1;  //Zupq1=482
                    double oldpb = 0;  //Zupq1=483
                    double argser = 1;  //Zupq1=484
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=485
                        argser = argser*argb*argb/4.0;  //Zupq1=486
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=487 */
                        pb = pb+carr1[i]*argser;  //Zupq1=488
                        const double delser = fabs((pb-oldpb)/pb);  //Zupq1=489
                        if ( delser<0.0001 ) break; /* goto 13; */  //Zupq1=490
                        oldpb = pb;  //Zupq1=491
                    }/*6*/  //Zupq1=492
                    /*13:*/  //Zupq1=493
#ifdef nichtVerwendet
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=494
#endif
                }/*5*/  //Zupq1=495
                else
                {/*5*/  //Zupq1=496
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/  /*  P(q)  */  //Zupq1=497
                        const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=498
                        const double pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=499
                        const double pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=500
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=501
                    }/*6*/  //Zupq1=502
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=503
                        const double pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=504
                        const double pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=505
                        const double pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=506
                        pb = pb3*pb3;  //Zupq1=507
                    }/*6*/  //Zupq1=508
#endif
                }/*5*/  //Zupq1=509
            }/*4*/  //Zupq1=510
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=511
                const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=512
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=513
                const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=514
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=515
            }/*4*/  //Zupq1=516
        }/*3*/  //Zupq1=517

        if ( i1==9 )    // nur bei orcase==2, immer i3=0
        {/*3*/   /*  disk, x-axis  */  //Zupq1=519
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=520
            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=521
            const double qxl = qx/qn;  //Zupq1=522
            const double qyl = qy/qn;  //Zupq1=523
            const double qnna = qxl*cos(delta)-qyl*sin(a)*sin(delta);  //Zupq1=524
            const double qnnb = qxl*cos(delta)-qyl*sin(b)*sin(delta);  //Zupq1=525
            const double arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=526
            const double argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=527

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=529
                if ( arga<0.015 )
                {/*5*/  //Zupq1=530
                    pa = 1;  //Zupq1=531
                    double oldpa = 0;  //Zupq1=532
                    double argser = 1;  //Zupq1=533
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=534
                        argser = argser*arga*arga/4.0;  //Zupq1=535
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=536 */
                        pa = pa+carr1[i]*argser;  //Zupq1=537
                        const double delser = fabs((pa-oldpa)/pa);  //Zupq1=538
                        if ( delser<0.0001 ) break; /* goto 15; */  //Zupq1=539
                        oldpa = pa;  //Zupq1=540
                    }/*6*/  //Zupq1=541
                    /*15:*/  //Zupq1=542
#ifdef nichtVerwendet
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=543
#endif
                }/*5*/  //Zupq1=544
                else
                {/*5*/  //Zupq1=545
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/ /*  P(q)  */  //Zupq1=546
                        const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=547
                        const double pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=548
                        const double pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=549
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=550
                    }/*6*/  //Zupq1=551
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=552
                        const double pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=553
                        const double pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=554
                        const double pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=555
                        pa = pa3*pa3;  //Zupq1=556
                    }/*6*/  //Zupq1=557
#endif
                }/*5*/  //Zupq1=558
                if ( argb<0.015 )
                {/*5*/  //Zupq1=559
                    pb = 1;  //Zupq1=560
                    double oldpb = 0;  //Zupq1=561
                    double argser = 1;  //Zupq1=562
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=563
                        argser = argser*argb*argb/4.0;  //Zupq1=564
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=565 */
                        pb = pb+carr1[i]*argser;  //Zupq1=566
                        const double delser = fabs((pb-oldpb)/pb);  //Zupq1=567
                        if ( delser<0.0001 ) break; /* goto 16; */  //Zupq1=568
                        oldpb = pb;  //Zupq1=569
                    }/*6*/  //Zupq1=570
                    /*16:*/  //Zupq1=571
#ifdef nichtVerwendet
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=572
#endif
                }/*5*/  //Zupq1=573
                else
                {/*5*/  //Zupq1=574
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/  /*  P(q)  */  //Zupq1=575
                        const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=576
                        const double pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=577
                        const double pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=578
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=579
                    }/*6*/  //Zupq1=580
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=581
                        const double pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=582
                        const double pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=583
                        const double pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=584
                        pb = pb3*pb3;  //Zupq1=585
                    }/*6*/  //Zupq1=586
#endif
                }/*5*/  //Zupq1=587
            }/*4*/  //Zupq1=588
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=589
                const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=590
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=591
                const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=592
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=593
            }/*4*/  //Zupq1=594
        }/*3*/  //Zupq1=595

        if ( i1==10 )   // nur bei orcase==3, immer i3=0
        {/*3*/   /*  disk, y-axis  */  //Zupq1=597
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=598
            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=599
            const double qxl = qx/qn;  //Zupq1=600
            const double qyl = qy/qn;  //Zupq1=601
            const double qnna = qxl*sin(a)*sin(delta)+qyl*cos(delta);  //Zupq1=602
            const double qnnb = qxl*sin(b)*sin(delta)+qyl*cos(delta);  //Zupq1=603
            const double arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=604
            const double argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=605

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=607
                if ( arga<0.015 )
                {/*5*/  //Zupq1=608
                    pa = 1;  //Zupq1=609
                    double oldpa = 0;  //Zupq1=610
                    double argser = 1;  //Zupq1=611
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=612
                        argser = argser*arga*arga/4.0;  //Zupq1=613
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=614 */
                        pa = pa+carr1[i]*argser;  //Zupq1=615
                        const double delser = fabs((pa-oldpa)/pa);  //Zupq1=616
                        if ( delser<0.0001 ) break; /* goto 18; */  //Zupq1=617
                        oldpa = pa;  //Zupq1=618
                    }/*6*/  //Zupq1=619
                    /*18:*/  //Zupq1=620
#ifdef nichtVerwendet
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=621
#endif
                }/*5*/  //Zupq1=622
                else
                {/*5*/  //Zupq1=623
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/ /*  P(q)  */  //Zupq1=624
                        const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=625
                        const double pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=626
                        const double pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=627
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=628
                    }/*6*/  //Zupq1=629
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=630
                        const double pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=631
                        const double pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=632
                        const double pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=633
                        pa = pa3*pa3;  //Zupq1=634
                    }/*6*/  //Zupq1=635
#endif
                }/*5*/  //Zupq1=636
                if ( argb<0.015 )
                {/*5*/  //Zupq1=637
                    pb = 1;  //Zupq1=638
                    double oldpb = 0;  //Zupq1=639
                    double argser = 1;  //Zupq1=640
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=641
                        argser = argser*argb*argb/4.0;  //Zupq1=642
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=643 */
                        pb = pb+carr1[i]*argser;  //Zupq1=644
                        const double delser = fabs((pb-oldpb)/pb);  //Zupq1=645
                        if ( delser<0.0001 ) break; /* goto 19; */  //Zupq1=646
                        oldpb = pb;  //Zupq1=647
                    }/*6*/  //Zupq1=648
                    /*19:*/  //Zupq1=649
#ifdef nichtVerwendet
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=650
#endif
                }/*5*/  //Zupq1=651
                else
                {/*5*/  //Zupq1=652
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/  /*  P(q)  */  //Zupq1=653
                        const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=654
                        const double pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=655
                        const double pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=656
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=657
                    }/*6*/  //Zupq1=658
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=659
                        const double pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=660
                        const double pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=661
                        const double pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=662
                        pb = pb3*pb3;  //Zupq1=663
                    }/*6*/  //Zupq1=664
#endif
                }/*5*/  //Zupq1=665
            }/*4*/  //Zupq1=666
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=667
                const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=668
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=669
                const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=670
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=671
            }/*4*/  //Zupq1=672
        }/*3*/  //Zupq1=673

        if ( i1==11 )   // nur bei orcase==4, immer i3=0
        {/*3*/   /*  disk, -z-axis  */  //Zupq1=675
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=676
            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=677
            const double qxl = qx/qn;  //Zupq1=678
            const double qyl = qy/qn;  //Zupq1=679
            const double qnna = qxl*sin(a)*sin(delta)-qyl*cos(a)*sin(delta);  //Zupq1=680
            const double qnnb = qxl*sin(b)*sin(delta)-qyl*cos(b)*sin(delta);  //Zupq1=681
            const double arga = sqrt(1.0-qnna*qnna+eps)*l*qn/(z+1);  //Zupq1=682
            const double argb = sqrt(1.0-qnnb*qnnb+eps)*l*qn/(z+1);  //Zupq1=683

            if ( sigma<0.15 )
            {/*4*/  /*  series expansion/asymptote  */  //Zupq1=685
                if ( arga<0.015 )
                {/*5*/  //Zupq1=686
                    pa = 1;  //Zupq1=687
                    double oldpa = 0;  //Zupq1=688
                    double argser = 1;  //Zupq1=689
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=690
                        argser = argser*arga*arga/4.0;  //Zupq1=691
                        /* pa:=pa+carr1[i]*power(arga/2,2*i);  //Zupq1=692 */
                        pa = pa+carr1[i]*argser;  //Zupq1=693
                        const double delser = fabs((pa-oldpa)/pa);  //Zupq1=694
                        if ( delser<0.0001 ) break; /* goto 21; */  //Zupq1=695
                        oldpa = pa;  //Zupq1=696
                    }/*6*/  //Zupq1=697
                    /*21:*/  //Zupq1=698
#ifdef nichtVerwendet
                    if ( i3==1 ) pa = pa*pa;  //Zupq1=699
#endif
                }/*5*/  //Zupq1=700
                else
                {/*5*/  //Zupq1=701
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/ /*  P(q)  */  //Zupq1=702
                        const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=703
                        const double pa2 = (1/(z*(z-1)*(z-2)))*pow(arga,-3)*sin((z-2)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-2)/2.0);  //Zupq1=704
                        const double pa3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(arga,-4)*cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0);  //Zupq1=705
                        pa = (4/M_PI)*(pa1-pa2-(9/8.0)*pa3);  //Zupq1=706
                    }/*6*/  //Zupq1=707
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=708
                        const double pa1 = (gamma(z-1/2.0)/gamma(z+1))*pow(arga,-3/2.0)*(sin((z-1/2.0)*atan(arga))-cos((z-1/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-1/2.0)/2.0);  //Zupq1=709
                        const double pa2 = (gamma(z-3/2.0)/gamma(z+1))*pow(arga,-5/2.0)*(sin((z-3/2.0)*atan(arga))+cos((z-3/2.0)*atan(arga)))/pow(1.0+arga*arga,(z-3/2.0)/2.0);  //Zupq1=710
                        const double pa3 = (2/sqrt(M_PI))*(pa1+(9/16.0)*pa2);  //Zupq1=711
                        pa = pa3*pa3;  //Zupq1=712
                    }/*6*/  //Zupq1=713
#endif
                }/*5*/  //Zupq1=714
                if ( argb<0.015 )
                {/*5*/  //Zupq1=715
                    pb = 1;  //Zupq1=716
                    double oldpb = 0;  //Zupq1=717
                    double argser = 1;  //Zupq1=718
                    for ( int i=1; i<=50; i++ )
                    {/*6*/  //Zupq1=719
                        argser = argser*argb*argb/4.0;  //Zupq1=720
                        /* pb:=pb+carr1[i]*power(argb/2,2*i);  //Zupq1=721 */
                        pb = pb+carr1[i]*argser;  //Zupq1=722
                        const double delser = fabs((pb-oldpb)/pb);  //Zupq1=723
                        if ( delser<0.0001 ) break; /* goto 22; */  //Zupq1=724
                        oldpb = pb;  //Zupq1=725
                    }/*6*/  //Zupq1=726
                    /*22:*/  //Zupq1=727
#ifdef nichtVerwendet
                    if ( i3==1 ) pb = pb*pb;  //Zupq1=728
#endif
                }/*5*/  //Zupq1=729
                else
                {/*5*/  //Zupq1=730
#ifdef nichtVerwendet
                    if ( i3==0 )
#endif
                    {/*6*/  /*  P(q)  */  //Zupq1=731
                        const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=732
                        const double pb2 = (1/(z*(z-1)*(z-2)))*pow(argb,-3)*sin((z-2)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-2)/2.0);  //Zupq1=733
                        const double pb3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argb,-4)*cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0);  //Zupq1=734
                        pb = (4/M_PI)*(pb1-pb2-(9/8.0)*pb3);  //Zupq1=735
                    }/*6*/  //Zupq1=736
#ifdef nichtVerwendet
                    if ( i3==1 )
                    {/*6*/  /*  F(q)  */  //Zupq1=737
                        const double pb1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argb,-3/2.0)*(sin((z-1/2.0)*atan(argb))-cos((z-1/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-1/2.0)/2.0);  //Zupq1=738
                        const double pb2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argb,-5/2.0)*(sin((z-3/2.0)*atan(argb))+cos((z-3/2.0)*atan(argb)))/pow(1.0+argb*argb,(z-3/2.0)/2.0);  //Zupq1=739
                        const double pb3 = (2/sqrt(M_PI))*(pb1+(9/16.0)*pb2);  //Zupq1=740
                        pb = pb3*pb3;  //Zupq1=741
                    }/*6*/  //Zupq1=742
#endif
                }/*5*/  //Zupq1=743
            }/*4*/  //Zupq1=744
            else
            {/*4*/  /*  OZ-type  */  //Zupq1=745
                const double pa1 = (1/(z*(z-1)*(z-2)))*pow(arga,-3);  //Zupq1=746
                pa = 1/(1.0+M_PI/(4.0*pa1));  //Zupq1=747
                const double pb1 = (1/(z*(z-1)*(z-2)))*pow(argb,-3);  //Zupq1=748
                pb = 1/(1.0+M_PI/(4.0*pb1));  //Zupq1=749
            }/*4*/  //Zupq1=750
        }/*3*/  //Zupq1=751

        if ( i1==12 )   // immer i3=1
        {/*3*/   /*  isotropic cube  */  //Zupq1=753
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=754

            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=787
            const double arga1 = qn*sin(delta)*cos(a)*r/(z+1)+eps;  //Zupq1=788
            const double arga2 = qn*sin(delta)*sin(a)*r/(z+1)+eps;  //Zupq1=789
            const double arga3 = qn*cos(delta)*r/(z+1)+eps;  //Zupq1=790
#ifdef nichtVerwendet
            if ( i3==0 )
            {/*4*/  /*  P(q)  */  //Zupq1=791
                const double a1 = 1/(2.0*z*(z-1));  //Zupq1=792
                double pa1, pa2, pa3;
                if ( arga1 < 0.001 )
                    pa1 = 1;
                else  //Zupq1=793
                    pa1 = (a1/(arga1*arga1+eps))*(1-cos((z-1)*atan(2.0*arga1))/pow(1.0+4*arga1*arga1,(z-1)/2.0));  //Zupq1=794
                if ( arga2 < 0.001 )
                    pa2 = 1;
                else  //Zupq1=795
                    pa2 = (a1/(arga2*arga2+eps))*(1-cos((z-1)*atan(2.0*arga2))/pow(1.0+4*arga2*arga2,(z-1)/2.0));  //Zupq1=796
                if ( arga3 < 0.001 )
                    pa3 = 1;
                else  //Zupq1=797
                    pa3 = (a1/(arga3*arga3+eps))*(1-cos((z-1)*atan(2.0*arga3))/pow(1.0+4*arga3*arga3,(z-1)/2.0));  //Zupq1=798
                pa = pa1*pa2*pa3;  //Zupq1=799
            }/*4*/  //Zupq1=800
#endif
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=801
                const double a1 = 1/z;  //Zupq1=802
                double pa1, pa2, pa3;
                if ( arga1 < 0.001 ) pa1 = 1; else  //Zupq1=803
                    pa1 = (a1/(arga1+eps))*sin(z*atan(arga1))/pow(1.0+arga1*arga1,z/2.0);  //Zupq1=804
                if ( arga2 < 0.001 ) pa2 = 1; else  //Zupq1=805
                    pa2 = (a1/(arga2+eps))*sin(z*atan(arga2))/pow(1.0+arga2*arga2,z/2.0);  //Zupq1=806
                if ( arga3 < 0.001 ) pa3 = 1; else  //Zupq1=807
                    pa3 = (a1/(arga3+eps))*sin(z*atan(arga3))/pow(1.0+arga3*arga3,z/2.0);  //Zupq1=808
                pa = pa1*pa1*pa2*pa2*pa3*pa3;  //Zupq1=809
            }/*4*/  //Zupq1=810
            const double argb1 = qn*sin(delta)*cos(b)*r/(z+1);  //Zupq1=811
            const double argb2 = qn*sin(delta)*sin(b)*r/(z+1);  //Zupq1=812
            const double argb3 = qn*cos(delta)*r/(z+1);  //Zupq1=813
#ifdef nichtVerwendet
            if ( i3==0 )
            {/*4*/   /*  P(q)  */  //Zupq1=814
                const double a1 = 1/(2.0*z*(z-1));  //Zupq1=815
                double pb1, pb2, pb3;
                if ( argb1 < 0.001 ) pb1 = 1; else  //Zupq1=816
                    pb1 = (a1/(argb1*argb1+eps))*(1-cos((z-1)*atan(2.0*argb1))/pow(1.0+4*argb1*argb1,(z-1)/2.0));  //Zupq1=817
                if ( argb2 < 0.001 ) pb2 = 1; else  //Zupq1=818
                    pb2 = (a1/(argb2*argb2+eps))*(1-cos((z-1)*atan(2.0*argb2))/pow(1.0+4*argb2*argb2,(z-1)/2.0));  //Zupq1=819
                if ( argb3 < 0.001 ) pb3 = 1; else  //Zupq1=820
                    pb3 = (a1/(argb3*argb3+eps))*(1-cos((z-1)*atan(2.0*argb3))/pow(1.0+4*argb3*argb3,(z-1)/2.0));  //Zupq1=821
                pb = pb1*pb2*pb3;  //Zupq1=822
            }/*4*/  //Zupq1=823
#endif
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=824
                const double a1 = 1/z;  //Zupq1=825
                double pb1, pb2, pb3;
                if ( argb1 < 0.001 ) pb1 = 1; else  //Zupq1=826
                    pb1 = (a1/(argb1+eps))*sin(z*atan(argb1))/pow(1.0+argb1*argb1,z/2.0);  //Zupq1=827
                if ( argb2 < 0.001 ) pb2 = 1; else  //Zupq1=828
                    pb2 = (a1/(argb2+eps))*sin(z*atan(argb2))/pow(1.0+argb2*argb2,z/2.0);  //Zupq1=829
                if ( arga3 < 0.001 ) pb3 = 1; else  //Zupq1=830
                    pb3 = (a1/(argb3+eps))*sin(z*atan(argb3))/pow(1.0+argb3*argb3,z/2.0);  //Zupq1=831
                pb = pb1*pb1*pb2*pb2*pb3*pb3;  //Zupq1=832
            }/*4*/  //Zupq1=833
        }/*3*/  //Zupq1=834

        if ( i1==13 )
        {/*3*/   /*  biaxial ellipsoid, isotropic  */  //Zupq1=835
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=836
            const double epsi = l/r;  //Zupq1=837
            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=838
            const double arga = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(a)*cos(a))/(z+1);  //Zupq1=839
            const double a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=840
            const double pa1 = a1*pow(arga,-4)*(1+cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0));  //Zupq1=841
            const double pa2 = (a1/(z-4))*pow(arga,-5)*sin((z-4)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-4)/2.0);  //Zupq1=842
            const double pa3 = (a1/((z-4)*(z-5)))*pow(arga,-6)*(1-cos((z-5)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-5)/2.0));  //Zupq1=843
            pa = 9*(pa1-2*pa2+pa3)*sin(a);  //Zupq1=844
            const double argb = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(b)*cos(b))/(z+1);  //Zupq1=845
            const double pb1 = a1*pow(argb,-4)*(1+cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0));  //Zupq1=846
            const double pb2 = (a1/(z-4))*pow(argb,-5)*sin((z-4)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-4)/2.0);  //Zupq1=847
            const double pb3 = (a1/((z-4)*(z-5)))*pow(argb,-6)*(1-cos((z-5)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-5)/2.0));  //Zupq1=848
            pb = 9*(pb1-2*pb2+pb3)*sin(b);  //Zupq1=849
        }/*3*/  /*  of biaxial ellipsoid  */  //Zupq1=850

        if ( i1==14 )
        {/*3*/   /*  triaxial ellipsoid, isotropic  */  //Zupq1=852
            const double ella = r;  //Zupq1=853
            const double ellb = l;  //Zupq1=854
            const double ellc = r/p1;  //Zupq1=855
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=856
            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=857
            const double arga = qn*sqrt(pow(ella*cos(a)*sin(delta),2)+pow(ellb*sin(a)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);  //Zupq1=858
            const double a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=859
            const double pa1 = a1*pow(arga,-4)*(1+cos((z-3)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-3)/2.0));  //Zupq1=860
            const double pa2 = (a1/(z-4))*pow(arga,-5)*sin((z-4)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-4)/2.0);  //Zupq1=861
            const double pa3 = (a1/((z-4)*(z-5)))*pow(arga,-6)*(1-cos((z-5)*atan(2.0*arga))/pow(1.0+4*arga*arga,(z-5)/2.0));  //Zupq1=862
            pa = 9*(pa1-2*pa2+pa3);  //Zupq1=863
            /* pa:=power(3*(sin(qn*r)-qn*r*cos(qn*r))/(qn*qn*qn*r*r*r+eps),2);   }  //Zupq1=864 */
            /* pa:=1.05;  //Zupq1=865 */
            const double argb = qn*sqrt(pow(ella*cos(b)*sin(delta),2)+pow(ellb*sin(b)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);  //Zupq1=866
            const double pb1 = a1*pow(argb,-4)*(1+cos((z-3)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-3)/2.0));  //Zupq1=867
            const double pb2 = (a1/(z-4))*pow(argb,-5)*sin((z-4)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-4)/2.0);  //Zupq1=868
            const double pb3 = (a1/((z-4)*(z-5)))*pow(argb,-6)*(1-cos((z-5)*atan(2.0*argb))/pow(1.0+4*argb*argb,(z-5)/2.0));  //Zupq1=869
            pb = 9*(pb1-2*pb2+pb3);  //Zupq1=870
            /* pb:=power(3*(sin(qn*r*1.04)-qn*r*cos(qn*r))/(qn*qn*qn*r*r*r+eps),2);   }  //Zupq1=871 */
            /* pb:=1.04;  //Zupq1=872 */
        }/*3*/  /*  of ellipsoid  */  //Zupq1=873

        if ( i1==15 )
        {/*3*/   /*  barrel area integration  */  //Zupq1=875
            const double pa1 = r*pow(1.0-pow(a/l,p1),1/p1);  //Zupq1=876
            const double arga = -r*pow(a/l,p1-1)*pow(1.0-pow(a/l,p1),(1/p1)-1)/l;  //Zupq1=877
            pa = pa1*sqrt(1.0+arga*arga);  //Zupq1=878
            const double pb1 = r*pow(1.0-pow(b/l,p1),1/p1);  //Zupq1=879
            const double argb = -r*pow(b/l,p1-1)*pow(1.0-pow(b/l,p1),(1/p1)-1)/l;  //Zupq1=880
            pb = pb1*sqrt(1.0+argb*argb);  //Zupq1=881
        }/*3*/  //Zupq1=882

        if ( i1==16 )   // immer i3=0
        {/*3*/   /*  barrel, x-axis  */  //Zupq1=884
            const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=885
            const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=886
            const double arga = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(a)*sin(delta)+eps));  //Zupq1=887
            const double argb = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(b)*sin(delta)+eps));  //Zupq1=888
#ifdef nichtVerwendet
            if ( i3==0 )
#endif
            {/*4*/ /*  P(q)  */  //Zupq1=889
                const double a1 = 9*pow(z+1,4)/(2.0*z*(z-1)*(z-2)*(z-3));  //Zupq1=890
                pa = (a1/(arga*arga));  //Zupq1=891
                pb = (a1/(argb*argb));  //Zupq1=892
            }/*4*/  //Zupq1=893
#ifdef nichtVerwendet
            if ( i3==1 )
            {/*4*/ /*  F(q)  */  //Zupq1=894
                const double pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=895
                pa = pa1*pa1;  //Zupq1=896
                const double pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=897
                pb = pb1*pb1;  //Zupq1=898
            }/*4*/  //Zupq1=899
#endif
        }/*3*/  //Zupq1=900

        if ( i1==17 )
        {/*3*/   /*  superball integration  */  //Zupq1=902
            const double aa = r;  //Zupq1=903
            const double bb = p1;  //Zupq1=904
            const double cc = l;  //Zupq1=905
            /* arga1:=-(l/r)*power(delta/r,p1-1)*power(1-power(delta/r,p1)-power(a/r,p1),(1/p1)-1);  //Zupq1=906 */
            /* arga2:=-(l/r)*power(a/r,p1-1)*power(1-power(delta/r,p1)-power(a/r,p1),(1/p1)-1);  //Zupq1=907 */
            /* pa:=sqrt(1+arga1*arga1+arga2*arga2);  //Zupq1=908 */
            /* pa:=power(r,4)*power(l,2)*(power(power(r*r*cos(delta),p1)+(power(cos(a),p1)+power(sin(a),p1))*power(r*l*sin(delta),p1),-(2+p1)/p1))*  //Zupq1=909 */
            /*       sqrt(power(r,4*p1)*power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(a),2*p1-2)+power(sin(a),2*p1-2))*power(r*l*sin(delta),2*p1));  //Zupq1=910 */
            pa = aa*aa*bb*bb*cc*cc*(pow(pow(aa*bb*cos(delta),alfa)+(pow(bb*cc*cos(a),alfa)+pow(aa*cc*sin(a),alfa))*pow(sin(delta),alfa),-(2+alfa)/alfa))*  //Zupq1=911
                 sqrt(pow(aa*bb,2*alfa)*pow(cos(delta),2*alfa-2)*pow(sin(delta),2)+(pow(bb*cc,2*alfa)*pow(cos(a),2*alfa-2)+pow(aa*cc,2*alfa)*pow(sin(a),2*alfa-2))*pow(sin(delta),2*alfa));  //Zupq1=912
            /* pa:=r*r*(power(power(cos(delta),p1)+(power(cos(a),p1)+power(sin(a),p1))*power(sin(delta),p1),-(2+p1)/p1))*  //Zupq1=913 */
            /*       sqrt(power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(a),2*p1-2)+power(sin(a),2*p1-2))*power(sin(delta),2*p1));  //Zupq1=914 */
            /* argb1:=-(l/r)*power(delta/r,p1-1)*power(1-power(delta/r,p1)-power(b/r,p1),(1/p1)-1);  //Zupq1=915 */
            /* argb2:=-(l/r)*power(b/r,p1-1)*power(1-power(delta/r,p1)-power(b/r,p1),(1/p1)-1);  //Zupq1=916 */
            /* pb:=sqrt(1+argb1*argb1+argb2*argb2);  //Zupq1=917 */
            /* pb:=power(r,4)*power(l,2)*(power(power(r*r*cos(delta),p1)+(power(cos(b),p1)+power(sin(b),p1))*power(r*l*sin(delta),p1),-(2+p1)/p1))*  //Zupq1=918 */
            /*       sqrt(power(r,4*p1)*power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(b),2*p1-2)+power(sin(b),2*p1-2))*power(r*l*sin(delta),2*p1));  //Zupq1=919 */
            pb = aa*aa*bb*bb*cc*cc*(pow(pow(aa*bb*cos(delta),alfa)+(pow(bb*cc*cos(b),alfa)+pow(aa*cc*sin(b),alfa))*pow(sin(delta),alfa),-(2+alfa)/alfa))*  //Zupq1=920
                 sqrt(pow(aa*bb,2*alfa)*pow(cos(delta),2*alfa-2)*pow(sin(delta),2)+(pow(bb*cc,2*alfa)*pow(cos(b),2*alfa-2)+pow(aa*cc,2*alfa)*pow(sin(b),2*alfa-2))*pow(sin(delta),2*alfa));  //Zupq1=921
            /* pb:=r*r*(power(power(cos(delta),p1)+(power(cos(b),p1)+power(sin(b),p1))*power(sin(delta),p1),-(2+p1)/p1))*  //Zupq1=922 */
            /*       sqrt(power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(b),2*p1-2)+power(sin(b),2*p1-2))*power(sin(delta),2*p1));  //Zupq1=923 */
        }/*3*/  //Zupq1=924

        pq = 0.5*(b-a)*(pa+pb);  //Zupq1=926
        trapzdchid_cnt = 1;  //Zupq1=927
    }/*2 if n==1 */  //Zupq1=928

    else    // =======================================================

    {/*2*/  //Zupq1=929
        double px=0;
        const double tnm = trapzdchid_cnt;  //Zupq1=930
        const double del = (b-a)/tnm;  //Zupq1=931
        double x = a+0.5*del;  //Zupq1=932
        double sump = 0.0;  //Zupq1=933
        for ( int j=1; j<=trapzdchid_cnt; j++ )
        {/*3*/  //Zupq1=934
            if ( i1==1 )
            {/*4*/  /*  cylinder, general case  */  //Zupq1=935
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=936
                const double mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);  //Zupq1=937
                const double mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);  //Zupq1=938
                /* mlx3:=p31*cos(x)*sin(delta)+p33*cos(delta);  //Zupq1=939 */
                const double argx = (qx*mlx1+qy*mlx2+eps)*l/(z+1);  //Zupq1=940
#ifdef nichtVerwendet
                if ( i3==0 )
#endif
                {   /*  P(q)  */  //Zupq1=941
                    const double a1 = 1/(2.0*z*(z-1));  //Zupq1=942
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=943
                }   //Zupq1=944
#ifdef nichtVerwendet
                if ( i3==1 )
                {   /*  F(q)  */  //Zupq1=945
                    const double px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=946
                    px = px1*px1;  //Zupq1=947
                }   //Zupq1=948
#endif
            }/*4*/  //Zupq1=949

            if ( i1==2 )
            {/*4*/   /*  cylinder, x-axis  */  //Zupq1=950
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=951
                const double argx = (qx*cos(delta)-qy*sin(x)*sin(delta)+eps)*l/(z+1);  //Zupq1=952
#ifdef nichtVerwendet
                if ( i3==0 )
#endif
                {   /*  P(q)  */  //Zupq1=953
                    const double a1 = 1/(2.0*z*(z-1));  //Zupq1=954
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=955
                }   //Zupq1=956
#ifdef nichtVerwendet
                if ( i3==1 )
                {   /*  F(q)  */  //Zupq1=957
                    const double px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=958
                    px = px1*px1;  //Zupq1=959
                }   //Zupq1=960
#endif
            }/*4*/  //Zupq1=961

            if ( i1==3 )
            {/*4*/   /*  cylinder, y-axis  */  //Zupq1=962
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=963
                //const double a1 = 1/(2.0*z*(z-1));  //Zupq1=964
                const double argx = (qx*sin(x)*sin(delta)+qy*cos(delta)+eps)*l/(z+1);  //Zupq1=965
#ifdef nichtVerwendet
                if ( i3==0 )
#endif
                {   /*  P(q)  */  //Zupq1=966
                    const double a1 = 1/(2.0*z*(z-1));  //Zupq1=967
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=968
                }   //Zupq1=969
#ifdef nichtVerwendet
                if ( i3==1 )
                {   /*  F(q)  */  //Zupq1=970
                    const double px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=971
                    px = px1*px1;  //Zupq1=972
                }   //Zupq1=973
#endif
            }/*4*/  //Zupq1=974

            if ( i1==4 )
            {/*4*/   /*  cylinder, -z-axis  */  //Zupq1=975
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=976
                //const double a1 = 1/(2.0*z*(z-1));  //Zupq1=977
                const double argx = (qx*sin(x)*sin(delta)-qy*cos(x)*sin(delta)+eps)*l/(z+1);  //Zupq1=978
#ifdef nichtVerwendet
                if ( i3==0 )
#endif
                {   /*  P(q)  */  //Zupq1=979
                    const double a1 = 1/(2.0*z*(z-1));  //Zupq1=980
                    px = (a1/(argx*argx))*(1-cos((z-1)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-1)/2.0));  //Zupq1=981
                }   //Zupq1=982
#ifdef nichtVerwendet
                if ( i3==1 )
                {   /*  F(q)  */  //Zupq1=983
                    const double px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=984
                    px = px1*px1;  //Zupq1=985
                }   //Zupq1=986
#endif
            }/*4*/  //Zupq1=987

#ifdef nichtVerwendet
            if ( i1==5 )
            {/*4*/   /*  general series expansion  */  //Zupq1=988
                const double mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);  //Zupq1=989
                const double mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);  //Zupq1=990
                //mlx3 = p31*cos(x)*sin(delta)+p33*cos(delta);  //Zupq1=991
                px = pow(mlx1,i3)*pow(mlx2,i4);  //Zupq1=992
            }/*4*/  //Zupq1=993
#endif
            if ( i1==6 )
            {/*4*/    /*  unit cell rotation  */  //Zupq1=994
                const double dqxx = qx-qhkl*(p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta));  //Zupq1=995
                const double dqyx = qy-qhkl*(p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta));  //Zupq1=996
                const double dqzx = qz-qhkl*(p31*cos(x)*sin(delta)+p33*cos(delta));  //Zupq1=997
                const double dqs1x = (dqxx*params.ax1.x()+dqyx*params.ax1.y()+dqzx*params.ax1.z())/(params.ax1.length()/*ax1n*/*params.sig.x());  //Zupq1=998
                const double dqs2x = (dqxx*params.ax2.x()+dqyx*params.ax2.y()+dqzx*params.ax2.z())/(params.ax2.length()/*ax2n*/*params.sig.y());  //Zupq1=999
                const double dqs3x = (dqxx*params.ax3.x()+dqyx*params.ax3.y()+dqzx*params.ax3.z())/(params.ax3.length()/*ax3n*/*params.sig.z());  //Zupq1=1000
                px = exp(-4*(dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x)/M_PI);  //Zupq1=1001
            }/*4*/  //Zupq1=1002

            if ( i1==7 )
            {/*4*/    /*  fiber unit cell rotation  */  //Zupq1=1003
                const double l1 = sin(theta)*cos(phi);  //Zupq1=1004
                const double l2 = sin(theta)*sin(phi);  //Zupq1=1005
                const double l3 = -cos(theta);  //Zupq1=1006

                const double r11x = cos(x)+(1-cos(x))*l1*l1;  //Zupq1=1008
                const double r12x = -l3*sin(x)+(1-cos(x))*l1*l2;  //Zupq1=1009
                const double r13x = -l2*sin(x)+(1-cos(x))*l1*l3;  //Zupq1=1010
                const double r21x = l3*sin(x)+(1-cos(x))*l1*l2;  //Zupq1=1011
                const double r22x = cos(x)+(1-cos(x))*l2*l2;  //Zupq1=1012
                const double r23x = l1*sin(x)+(1-cos(x))*l2*l3;  //Zupq1=1013
                const double r31x = l2*sin(x)+(1-cos(x))*l1*l3;  //Zupq1=1014
                const double r32x = -l1*sin(x)+(1-cos(x))*l2*l3;  //Zupq1=1015
                const double r33x = cos(x)+(1-cos(x))*l3*l3;  //Zupq1=1016

                /*  rotate this scattering vector  */  //Zupq1=1023
                const double qxhklx = r11x*qxn+r12x*qyn+r13x*qzn;  //Zupq1=1024
                const double qyhklx = r21x*qxn+r22x*qyn+r23x*qzn;  //Zupq1=1025
                const double qzhklx = r31x*qxn+r32x*qyn+r33x*qzn;  //Zupq1=1026

                const double dqxx = qx-qxhklx;  //Zupq1=1093
                const double dqyx = qy-qyhklx;  //Zupq1=1094
                const double dqzx = qz-qzhklx;  //Zupq1=1095

                const double ax1xx = (r11x*params.ax1.x()+r12x*params.ax1.y()+r13x*params.ax1.z());  //Zupq1=1097
                const double ax1yx = (r21x*params.ax1.x()+r22x*params.ax1.y()+r23x*params.ax1.z());  //Zupq1=1098
                const double ax1zx = (r31x*params.ax1.x()+r32x*params.ax1.y()+r33x*params.ax1.z());  //Zupq1=1099
                const double ax2xx = (r11x*params.ax2.x()+r12x*params.ax2.y()+r13x*params.ax2.z());  //Zupq1=1100
                const double ax2yx = (r21x*params.ax2.x()+r22x*params.ax2.y()+r23x*params.ax2.z());  //Zupq1=1101
                const double ax2zx = (r31x*params.ax2.x()+r32x*params.ax2.y()+r33x*params.ax2.z());  //Zupq1=1102
                const double ax3xx = (r11x*params.ax3.x()+r12x*params.ax3.y()+r13x*params.ax3.z());  //Zupq1=1103
                const double ax3yx = (r21x*params.ax3.x()+r22x*params.ax3.y()+r23x*params.ax3.z());  //Zupq1=1104
                const double ax3zx = (r31x*params.ax3.x()+r32x*params.ax3.y()+r33x*params.ax3.z());  //Zupq1=1105

                const double dqs1x = (dqxx*ax1xx+dqyx*ax1yx+dqzx*ax1zx)/(params.ax1.length()/*ax1n*/*params.sig.x());  //Zupq1=1107
                const double dqs2x = (dqxx*ax2xx+dqyx*ax2yx+dqzx*ax2zx)/(params.ax2.length()/*ax2n*/*params.sig.y());  //Zupq1=1108
                const double dqs3x = (dqxx*ax3xx+dqyx*ax3yx+dqzx*ax3zx)/(params.ax3.length()/*ax3n*/*params.sig.z());  //Zupq1=1109

                const double argx = dqs1x*dqs1x+dqs2x*dqs2x+dqs3x*dqs3x;  //Zupq1=1115
                px = exp(-4*argx/M_PI);  //Zupq1=1116
                /* if (argx>2) then px:=eps  //Zupq1=1117 */
                /*    else px:=exp(-4*argx/pi);  //Zupq1=1118 */
            }/*4*/  //Zupq1=1119

            if ( i1==8 )
            {/*4*/  /*  disk, general case  */  //Zupq1=1120
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1121
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1122
                const double qxl = qx/qn;  //Zupq1=1123
                const double qyl = qy/qn;  //Zupq1=1124
                const double mlx1 = p11*cos(x)*sin(delta)+p12*sin(x)*sin(delta)+p13*cos(delta);  //Zupq1=1125
                const double mlx2 = p21*sin(x)*sin(delta)+p22*cos(x)*sin(delta)+p23*cos(delta);  //Zupq1=1126
                /* mlx3:=p31*cos(x)*sin(delta)+p33*cos(delta);  //Zupq1=1127 */
                const double qnnx = qxl*mlx1+qyl*mlx2;  //Zupq1=1128
                const double argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1129

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1131
                    if ( argx<0.015 )
                    {/*6*/  //Zupq1=1132
                        px = 1;  //Zupq1=1133
                        double oldpx = 0;  //Zupq1=1134
                        double argser = 1;  //Zupq1=1135
                        for ( int i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1136
                            argser = argser*argx*argx/4.0;  //Zupq1=1137
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1138 */
                            px = px+carr1[i]*argser;  //Zupq1=1139
                            const double delser = fabs((px-oldpx)/px);  //Zupq1=1140
                            if ( delser<0.0001 ) break; /* goto 14; */  //Zupq1=1141
                            oldpx = px;  //Zupq1=1142
                        }/*7*/  //Zupq1=1143
                        /*14:*/  //Zupq1=1144
#ifdef nichtVerwendet
                        if ( i3==1 ) px = px*px;  //Zupq1=1145
#endif
                    }/*6*/  //Zupq1=1146
                    else
                    {/*6*/  //Zupq1=1147
#ifdef nichtVerwendet
                        if ( i3==0 )
#endif
                        {/*7*/  /*  P(q)  */  //Zupq1=1148
                            const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1149
                            const double px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1150
                            const double px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1151
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1152
                        }/*7*/  //Zupq1=1153
#ifdef nichtVerwendet
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1154
                            const double px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1155
                            const double px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1156
                            const double px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1157
                            px = px3*px3;  //Zupq1=1158
                        }/*7*/  //Zupq1=1159
#endif
                    }/*6*/  //Zupq1=1160
                }/*5*/  //Zupq1=1161
                else
                {/*5*/  //Zupq1=1162
                    const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1163
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1164
                }/*5*/  //Zupq1=1165
            }/*4*/  //Zupq1=1166

            if ( i1==9 )
            {/*4*/   /*  disk, x-axis  */  //Zupq1=1168
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1169
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1170
                const double qxl = qx/qn;  //Zupq1=1171
                const double qyl = qy/qn;  //Zupq1=1172
                const double qnnx = qxl*cos(delta)-qyl*sin(x)*sin(delta);  //Zupq1=1173
                const double argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1174

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1176
                    if ( argx<0.015 )
                    {/*6*/  /*  series expansion  */  //Zupq1=1177
                        px = 1;  //Zupq1=1178
                        double oldpx = 0;  //Zupq1=1179
                        double argser = 1;  //Zupq1=1180
                        for ( int i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1181
                            argser = argser*argx*argx/4.0;  //Zupq1=1182
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1183 */
                            px = px+carr1[i]*argser;  //Zupq1=1184
                            const double delser = fabs((px-oldpx)/px);  //Zupq1=1185
                            if ( delser<0.0001 ) break; /* goto 17; */  //Zupq1=1186
                            oldpx = px;  //Zupq1=1187
                        }/*7*/  //Zupq1=1188
                        /*17:*/  //Zupq1=1189
#ifdef nichtVerwendet
                        if ( i3==1 ) px = px*px;  //Zupq1=1190
#endif
                    }/*6*/  //Zupq1=1191
                    else
                    {/*6*/  //Zupq1=1192
#ifdef nichtVerwendet
                        if ( i3==0 )
#endif
                        {/*7*/  /*  P(q)  */  //Zupq1=1193
                            const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1194
                            const double px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1195
                            const double px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1196
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1197
                        }/*7*/  //Zupq1=1198
#ifdef nichtVerwendet
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1199
                            const double px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1200
                            const double px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1201
                            const double px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1202
                            px = px3*px3;  //Zupq1=1203
                        }/*7*/  //Zupq1=1204
#endif
                    }/*6*/  //Zupq1=1205
                }/*5*/  //Zupq1=1206
                else
                {/*5*/  //Zupq1=1207
                    const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1208 TODO war px
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1209
                }/*5*/  //Zupq1=1210
            }/*4*/  //Zupq1=1211

            if ( i1==10 )
            {/*4*/   /*  disk, y-axis  */  //Zupq1=1213
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1214
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1215
                const double qxl = qx/qn;  //Zupq1=1216
                const double qyl = qy/qn;  //Zupq1=1217
                const double qnnx = qxl*sin(x)*sin(delta)+qyl*cos(delta);  //Zupq1=1218
                const double argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1219

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1221
                    if ( argx<0.015 )
                    {/*6*/  //Zupq1=1222
                        px = 1;  //Zupq1=1223
                        double oldpx = 0;  //Zupq1=1224
                        double argser = 1;  //Zupq1=1225
                        for ( int i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1226
                            argser = argser*argx*argx/4.0;  //Zupq1=1227
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1228 */
                            px = px+carr1[i]*argser;  //Zupq1=1229
                            const double delser = fabs((px-oldpx)/px);  //Zupq1=1230
                            if ( delser<0.0001 ) break; /* goto 20; */  //Zupq1=1231
                            oldpx = px;  //Zupq1=1232
                        }/*7*/  //Zupq1=1233
                        /*20:*/  //Zupq1=1234
#ifdef nichtVerwendet
                        if ( i3==1 ) px = px*px;  //Zupq1=1235
#endif
                    }/*6*/  //Zupq1=1236
                    else
                    {/*6*/  //Zupq1=1237
#ifdef nichtVerwendet
                        if ( i3==0 )
#endif
                        {/*7*/  /*  P(q)  */  //Zupq1=1238
                            const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1239
                            const double px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1240
                            const double px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1241
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1242
                        }/*7*/  //Zupq1=1243
#ifdef nichtVerwendet
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1244
                            const double px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1245
                            const double px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1246
                            const double px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1247
                            px = px3*px3;  //Zupq1=1248
                        }/*7*/  //Zupq1=1249
#endif
                    }/*6*/  //Zupq1=1250
                }/*5*/  //Zupq1=1251
                else
                {/*5*/  //Zupq1=1252
                    const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1253
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1254
                }/*5*/  //Zupq1=1255
            }/*4*/  //Zupq1=1256

            if ( i1==11 )
            {/*4*/   /*  disk, z-axis  */  //Zupq1=1258
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1259
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1260
                const double qxl = qx/qn;  //Zupq1=1261
                const double qyl = qy/qn;  //Zupq1=1262
                const double qnnx = qxl*sin(x)*sin(delta)-qyl*cos(x)*sin(delta);  //Zupq1=1263
                const double argx = sqrt(1.0-qnnx*qnnx+eps)*l*qn/(z+1);  //Zupq1=1264

                if ( sigma<0.15 )
                {/*5*/  //Zupq1=1266
                    if ( argx<0.015 )
                    {/*6*/  //Zupq1=1267
                        px = 1;  //Zupq1=1268
                        double oldpx = 0;  //Zupq1=1269
                        double argser = 1;  //Zupq1=1270
                        for ( int i=1; i<=50; i++ )
                        {/*7*/  //Zupq1=1271
                            argser = argser*argx*argx/4.0;  //Zupq1=1272
                            /* px:=px+carr1[i]*power(argx/2,2*i);  //Zupq1=1273 */
                            px = px+carr1[i]*argser;  //Zupq1=1274
                            const double delser = fabs((px-oldpx)/px);  //Zupq1=1275
                            if ( delser<0.0001 ) break; /* goto 23; */  //Zupq1=1276
                            oldpx = px;  //Zupq1=1277
                        }/*7*/  //Zupq1=1278
                        /*23:*/  //Zupq1=1279
#ifdef nichtVerwendet
                        if ( i3==1 ) px = px*px;  //Zupq1=1280
#endif
                    }/*6*/  //Zupq1=1281
                    else
                    {/*6*/  //Zupq1=1282
#ifdef nichtVerwendet
                        if ( i3==0 )
#endif
                        {/*7*/  /*  P(q)  */  //Zupq1=1283
                            const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1284
                            const double px2 = (1/(z*(z-1)*(z-2)))*pow(argx,-3)*sin((z-2)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-2)/2.0);  //Zupq1=1285
                            const double px3 = (1/(z*(z-1)*(z-2)*(z-3)))*pow(argx,-4)*cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0);  //Zupq1=1286
                            px = (4/M_PI)*(px1-px2-(9/8.0)*px3);  //Zupq1=1287
                        }/*7*/  //Zupq1=1288
#ifdef nichtVerwendet
                        if ( i3==1 )
                        {/*7*/  /*  F(q)  */  //Zupq1=1289
                            const double px1 = (gamma(z-1/2.0)/gamma(z+1))*pow(argx,-3/2.0)*(sin((z-1/2.0)*atan(argx))-cos((z-1/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-1/2.0)/2.0);  //Zupq1=1290
                            const double px2 = (gamma(z-3/2.0)/gamma(z+1))*pow(argx,-5/2.0)*(sin((z-3/2.0)*atan(argx))+cos((z-3/2.0)*atan(argx)))/pow(1.0+argx*argx,(z-3/2.0)/2.0);  //Zupq1=1291
                            const double px3 = (2/sqrt(M_PI))*(px1+(9/16.0)*px2);  //Zupq1=1292
                            px = px3*px3;  //Zupq1=1293
                        }/*7*/  //Zupq1=1294
#endif
                    }/*6*/  //Zupq1=1295
                }/*5*/  //Zupq1=1296
                else
                {/*5*/  //Zupq1=1297
                    const double px1 = (1/(z*(z-1)*(z-2)))*pow(argx,-3);  //Zupq1=1298
                    px = 1/(1.0+M_PI/(4.0*px1));  //Zupq1=1299
                }/*5*/  //Zupq1=1300
            }/*4*/  //Zupq1=1301

            if ( i1==12 )
            {/*4*/   /*  isotropic cube  */  //Zupq1=1303
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1322
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1323
                const double argx1 = qn*sin(delta)*cos(x)*r/(z+1);  //Zupq1=1324
                const double argx2 = qn*sin(delta)*sin(x)*r/(z+1);  //Zupq1=1325
                const double argx3 = qn*cos(delta)*l/(z+1);  //Zupq1=1326
#ifdef nichtVerwendet
                if ( i3==0 )
                {/*5*/   /*  P(q)  */  //Zupq1=1327
                    const double a1 = 1/(2.0*z*(z-1));  //Zupq1=1328
                    double px1, px2, px3;
                    if ( argx1 < 0.001 ) px1 = 1; else  //Zupq1=1329
                        px1 = (a1/(argx1*argx1+eps))*(1-cos((z-1)*atan(2.0*argx1))/pow(1.0+4*argx1*argx1,(z-1)/2.0));  //Zupq1=1330
                    if ( argx2 < 0.001 ) px2 = 1; else  //Zupq1=1331
                        px2 = (a1/(argx2*argx2+eps))*(1-cos((z-1)*atan(2.0*argx2))/pow(1.0+4*argx2*argx2,(z-1)/2.0));  //Zupq1=1332
                    if ( argx3 < 0.001 ) px3 = 1; else  //Zupq1=1333
                        px3 = (a1/(argx3*argx3+eps))*(1-cos((z-1)*atan(2.0*argx3))/pow(1.0+4*argx3*argx3,(z-1)/2.0));  //Zupq1=1334
                    px = px1*px2*px3;  //Zupq1=1335
                }/*5*/  //Zupq1=1336
#endif
                if ( i3==1 )
                {/*5*/ /*  F(q)  */  //Zupq1=1337
                    const double a1 = 1/z;  //Zupq1=1338
                    double px1, px2, px3;
                    if ( argx1 < 0.001 ) px1 = 1; else  //Zupq1=1339
                        px1 = (a1/(argx1+eps))*sin(z*atan(argx1))/pow(1.0+argx1*argx1,z/2.0);  //Zupq1=1340
                    if ( argx2 < 0.001 ) px2 = 1; else  //Zupq1=1341
                        px2 = (a1/(argx2+eps))*sin(z*atan(argx2))/pow(1.0+argx2*argx2,z/2.0);  //Zupq1=1342
                    if ( argx3 < 0.001 ) px3 = 1; else  //Zupq1=1343 TODO: es stand arga3 hier...
                        px3 = (a1/(argx3+eps))*sin(z*atan(argx3))/pow(1.0+argx3*argx3,z/2.0);  //Zupq1=1344
                    px = px1*px1*px2*px2*px3*px3;  //Zupq1=1345
                }/*5*/  //Zupq1=1346
            }/*4*/  //Zupq1=1347

            if ( i1==13 )
            {/*4*/   /*  biaxial ellipsoid, isotropic  */  //Zupq1=1349
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1350
                const double epsi = l/r;  //Zupq1=1351
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1352
                const double argx = qn*r*sqrt(1.0+(epsi*epsi-1)*cos(x)*cos(x))/(z+1);  //Zupq1=1353
                const double a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=1354
                const double px1 = a1*pow(argx,-4)*(1+cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0));  //Zupq1=1355
                const double px2 = (a1/(z-4))*pow(argx,-5)*sin((z-4)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-4)/2.0);  //Zupq1=1356
                const double px3 = (a1/((z-4)*(z-5)))*pow(argx,-6)*(1-cos((z-5)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-5)/2.0));  //Zupq1=1357
                px = 9*(px1-2*px2+px3)*sin(x);  //Zupq1=1358
            }/*4*/  /*  of biaxial ellipsoid  */  //Zupq1=1359

            if ( i1==14 )
            {/*4*/   /*  triaxial ellipsoid, isotropic  */  //Zupq1=1361
                const double ella = r;  //Zupq1=1362
                const double ellb = l;  //Zupq1=1363
                const double ellc = r/p1;  //Zupq1=1364
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1365
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1366
                const double argx = qn*sqrt(pow(ella*cos(x)*sin(delta),2)+pow(ellb*sin(x)*sin(delta),2)+pow(ellc*cos(delta),2))/(z+1);  //Zupq1=1367
                const double a1 = (1/(2.0*z*(z-1)*(z-2)*(z-3)));  //Zupq1=1368
                const double px1 = a1*pow(argx,-4)*(1+cos((z-3)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-3)/2.0));  //Zupq1=1369
                const double px2 = (a1/(z-4))*pow(argx,-5)*sin((z-4)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-4)/2.0);  //Zupq1=1370
                const double px3 = (a1/((z-4)*(z-5)))*pow(argx,-6)*(1-cos((z-5)*atan(2.0*argx))/pow(1.0+4*argx*argx,(z-5)/2.0));  //Zupq1=1371
                px = 9*(px1-2*px2+px3);  //Zupq1=1372
            }/*4*/  /*  of triaxial ellipsoid  */  //Zupq1=1373

            if ( i1==15 )
            {/*4*/   /*  barrel area integration  */  //Zupq1=1375
                const double px1 = r*pow(1.0-pow(x/l,p1),1/p1);  //Zupq1=1376
                const double argx = -r*pow(x/l,p1-1)*pow(1.0-pow(x/l,p1),(1/p1)-1)/l;  //Zupq1=1377
                px = px1*sqrt(1.0+argx*argx);  //Zupq1=1378
            }/*4*/  //Zupq1=1379

            if ( i1==16 )
            {/*4*/   /*  barrel, x-axis  */  //Zupq1=1381
                const double z = (1-sigma*sigma)/(sigma*sigma);  //Zupq1=1382
                const double qn = sqrt(qx*qx+qy*qy+eps);  //Zupq1=1383
                const double argx = sqr(qn*r)+sqr(qn*l*((qx/qn)*cos(delta)-(qy/qn)*sin(x)*sin(delta)+eps));  //Zupq1=1384
#ifdef nichtVerwendet
                if ( i3==0 )
#endif
                {   /*  P(q)  */  //Zupq1=1385
                    const double a1 = 9*pow(z+1,4)/(2.0*z*(z-1)*(z-2)*(z-3));  //Zupq1=1386
                    px = (a1/(argx*argx));  //Zupq1=1387
                }   //Zupq1=1388
#ifdef nichtVerwendet
                if ( i3==1 )
                {   /*  F(q)  */  //Zupq1=1389
                    // pa und pb werden hier nicht mehr berechnet (Copy/Paste-Fehler)
                    //pa1 = (1/z)*(1/arga)*sin(z*atan(arga))/pow(1.0+arga*arga,z/2.0);  //Zupq1=1390
                    //pa = pa1*pa1;  //Zupq1=1391
                    //pb1 = (1/z)*(1/argb)*sin(z*atan(argb))/pow(1.0+argb*argb,z/2.0);  //Zupq1=1392
                    //pb = pb1*pb1;  //Zupq1=1393
                    const double px1 = (1/z)*(1/argx)*sin(z*atan(argx))/pow(1.0+argx*argx,z/2.0);  //Zupq1=1390
                    px = px1*px1;  //Zupq1=1391
                }   //Zupq1=1394
#endif
            }/*4*/  //Zupq1=1395

            if ( i1==17 )
            {/*4*/   /*  superball integration  */  //Zupq1=1397
                const double aa = r;     // von oben kopiert... (TODO?)
                const double bb = p1;
                const double cc = l;
                /*  argx1:=-(l/r)*power(delta/r,p1-1)*power(1-power(delta/r,p1)-power(x/r,p1),(1/p1)-1);  //Zupq1=1398 */
                /*  argx2:=-(l/r)*power(x/r,p1-1)*power(1-power(delta/r,p1)-power(x/r,p1),(1/p1)-1);  //Zupq1=1399 */
                /*  px:=sqrt(1+argx1*argx1+argx2*argx2);  //Zupq1=1400 */
                /* px:=power(r,4)*power(l,2)*(power(power(r*r*cos(delta),p1)+(power(cos(x),p1)+power(sin(x),p1))*power(r*l*sin(delta),p1),-(2+p1)/p1))*  //Zupq1=1401 */
                /*    sqrt(power(r,4*p1)*power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(x),2*p1-2)+power(sin(x),2*p1-2))*power(r*l*sin(delta),2*p1));  //Zupq1=1402 */
                px = aa*aa*bb*bb*cc*cc*(pow(pow(aa*bb*cos(delta),alfa)+(pow(bb*cc*cos(x),alfa)+pow(aa*cc*sin(x),alfa))*pow(sin(delta),alfa),-(2+alfa)/alfa))*  //Zupq1=1403
                     sqrt(pow(aa*bb,2*alfa)*pow(cos(delta),2*alfa-2)*pow(sin(delta),2)+(pow(bb*cc,2*alfa)*pow(cos(x),2*alfa-2)+pow(aa*cc,2*alfa)*pow(sin(x),2*alfa-2))*pow(sin(delta),2*alfa));  //Zupq1=1404
                /* px:=r*r*(power(power(cos(delta),p1)+(power(cos(x),p1)+power(sin(x),p1))*power(sin(delta),p1),-(2+p1)/p1))*  //Zupq1=1405 */
                /*    sqrt(power(cos(delta),2*p1-2)*power(sin(delta),2)+(power(cos(x),2*p1-2)+power(sin(x),2*p1-2))*power(sin(delta),2*p1));  //Zupq1=1406 */
            }/*4*/  //Zupq1=1407

            sump += px;   //Zupq1=1409
            x    += del;  //Zupq1=1410
        }/*3*/  //Zupq1=1411
        pq = 0.5*(pq+(b-a)*sump/tnm);  //Zupq1=1412
        trapzdchid_cnt = 2*trapzdchid_cnt;  //Zupq1=1413
    }/*2*/  //Zupq1=1414
}/*1*/  //Zupq1=1415


#ifdef procnotused
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::polycscube(double /*rho1*/, double /*rho2*/, double /*p1*/, double /*p2*/,
                                            double /*alf1*/, double /*alf2*/, double /*rn*/, double /*pf*/,
                                            double /*sigma*/, double /*q*/) const
{
    // Gespräch mit Prof. Förster (05.Jun.2023): Diese Routine wird noch nicht verwendet (bei Erweiterungen?).
    return 0;
}
#endif



#ifdef undef
/* ****************** Percus-Yevick S(q) ************************ */
// aus 20210628 - Upq.pas
#ifdef __CUDACC__
__host__ __device__
#endif
double SasCalc_GENERIC_calculation::spy( double q ) const
{
    double d,v,c,a0,a1,a2,cdr,i1,i2,i3,s1,s2,s3;

    d = 2*params.radius;
    v = 4.*M_PI*params.radius*params.radius*params.radius/3.;
    c = params.dbeta/v;
    a0 = -1*sqr(1+2*params.dbeta)/sqr(sqr(1-params.dbeta));
    a1 = 6.*params.dbeta*sqr(1+params.dbeta/2.)/sqr(sqr(1-params.dbeta));
    a2 = params.dbeta*a0/2.;
    i1 = sin(d*q)-d*q*cos(d*q);
    i2 = 2*cos(d*q)-sqr(d*q)*cos(d*q)+2*d*q*sin(d*q)-2;
    i3 = 6*sqr(params.radius*q)*cos(d*q)-3*cos(d*q)-2*sqr(sqr(params.radius*q))*cos(d*q)
         -3*d*q*sin(d*q)+4*params.radius*q*sqr(params.radius*q)*sin(d*q)+3;
    s1 = 4*M_PI*a0*i1/(q*q*q);
    s2 = 4*M_PI*a1*i2/(d*sqr(sqr(q)));
    s3 = 4*M_PI*a2*8*i3/(d*d*d*sqr(q)*sqr(q)*sqr(q));
    cdr = s1+s2+s3;
    return 1./(1.-c*cdr);
}
#endif


#include "sc_lib_formpq_partSphere.h"       //  0, cbpartSphere
#include "sc_lib_formpq_partCylinder.h"     //  1, cbpartCylinder
#include "sc_lib_formpq_partDisk.h"         //  2, cbpartDisk
//                                              3, cbpartVesicle
#include "sc_lib_formpq_partCube.h"         //  4, cbpartCube
#include "sc_lib_formpq_partEllips.h"       //  5, cbpartEllipsoide
#include "sc_lib_formpq_partTriaxEllips.h"  //  6, cbpartTriaxEllips
#include "sc_lib_formpq_partSuperEllips.h"  //  7, cbpartSuperEllips
#include "sc_lib_formpq_partSuperball.h"    //  8, cbpartSuperball
//                                              9, cbpartChain
//                                             10, cbpartkpchain

#include "sc_lib_formfq_partSphere.h"       //  0, cbpartSphere
#include "sc_lib_formfq_partCylinder.h"     //  1, cbpartCylinder
#include "sc_lib_formfq_partDisk.h"         //  2, cbpartDisk
//                                              3, cbpartVesicle
//#include "sc_lib_formfq_partCube.h"           4, cbpartCube - wird nicht verwendet
//                                              5, cbpartEllipsoide - müsste implementiert werden, wenn ltype!=None
//                                              6, cbpartTriaxEllips
//                                              7, cbpartSuperEllips
//                                              8, cbpartSuperball
//                                              9, cbpartChain
//                                             10, cbpartkpchain
