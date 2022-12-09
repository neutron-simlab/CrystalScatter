/** 
  */

//#include "sc_memory_gpu.h"
//#include <stdlib.h>
//#include <string.h>

//#include <iostream>
//#include <chrono>

//#ifndef __CUDACC__
//#include <QDebug>
//#endif

//#ifdef __CUDACC__
//#include <cuda.h>
//#endif

void CLASSLIB::initMemory()
{
    arrXYIntensity = nullptr;
    arrXYsize      = 0;
#ifdef COPY_FITDATA_TO_GPU
    arrDataForFit  = nullptr;
    arrDataSize    = 0;
#ifdef CALC_FQS_IN_GPU
    arrFqsForFit   = nullptr;
    arrFqsSize     = 0;
#endif
    arrDataForFitUsed = false;
#endif
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
    cudaDeviceReset();
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    std::cerr << "CLASSLIB::initMemory(): GPU device count: " << deviceCount << std::endl;
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
        std::cerr << "                 multiProcessorCount: " << deviceProp.multiProcessorCount << std::endl;
        std::cerr << "         maxThreadsPerMultiProcessor: " << deviceProp.maxThreadsPerMultiProcessor << std::endl;
        std::cerr << "                         l2CacheSize: " << deviceProp.l2CacheSize << " bytes" << std::endl;
        std::cerr << "            persistingL2CacheMaxSize: " << deviceProp.persistingL2CacheMaxSize << " bytes" << std::endl;
        std::cerr << "      directManagedMemAccessFromHost: " << deviceProp.directManagedMemAccessFromHost << std::endl;

        /* - \ref ::cudaDeviceProp::uuid "uuid" is a 16-byte unique identifier.
        * - \ref ::cudaDeviceProp::sharedMemPerBlock "sharedMemPerBlock" is the maximum amount of shared memory available to a thread block in bytes;
        * - \ref ::cudaDeviceProp::regsPerBlock "regsPerBlock" is the maximum number of 32-bit registers available to a thread block;
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
    //cudaDeviceGetLimit( &s, cudaLimitPrintfFifoSize );  diese Routine braucht sehr lange ...
    //std::cerr << "Printf Fifo Size: " << s << std::endl;

    //size_t memFree, memTotal;
    //cudaMemGetInfo( &memFree, &memTotal );
    //std::cerr << "Memory: free=" << memFree << " of " << memTotal << std::endl;

    noGPUavailable = deviceCount == 0;
#else
    noGPUavailable = true;
    std::cerr << "CLASSLIB::initMemory(): no nvcc used" << std::endl;
#endif
}



void CLASSLIB::createMemory( void **ptr, size_t lensoll, size_t &lenist, bool gpuonly, const char *dbgInfo )
{
#ifndef __CUDACC__
    //Q_UNUSED(dbgInfo);

#else

//#define MEMALLOC(p,s) cudaMallocManaged(p,s)      //ist wohl doch etwas langsamer...
//#define MEMFREE(p)    cudaFree(p)                 //und nach einem Crash wird nicht richtig aufgeräumt

#define MEMALLOC(p,s) cudaMallocHost(p,s)
#define MEMFREE(p)    cudaFreeHost(p)

#endif
    if ( *ptr != nullptr )
    {
        if ( lensoll == lenist ) return;
        // realloc new size
        if ( noGPUavailable )
        {
            std::cerr << "re-allocate " << lensoll << " Bytes CPU-Memory";
            if ( dbgInfo != nullptr ) std::cerr << " (" << dbgInfo << ")";
            std::cerr << std::endl;
            *ptr = realloc( *ptr, lensoll );
        }
#ifdef __CUDACC__
        else
        {
            std::cerr << "re-allocate " << lensoll << " Bytes GPU-Memory";
            if ( dbgInfo != nullptr ) std::cerr << " (" << dbgInfo << ")";
            std::cerr << std::endl;
            // Da der Speicher immer vor den Berechnungen angelegt wird und die alten Daten
            // nicht verwendet werden (hier nur das Image), muss nichts kopiert werden.
            // Das klappte nämlich nicht wirklich - zumindest bei den ersten Versuchen
            cudaError_t e;
            e = MEMFREE( *ptr );
            if ( e != CUDA_SUCCESS )
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



/**
 * @brief SC_GpuMemory::checkXYarray - helper function to initialize the dynamic result array
 * @param minx - minimal horizontal index (incl.)
 * @param maxx - maximal horizontal index (incl.)
 * @param miny - minimal vertical index (incl.)
 * @param maxy - maximal vertical index (incl.)
 */
void CLASSLIB::checkArrays( int minx, int maxx, int miny, int maxy )
{
    _arrCount = (maxx-minx) * (maxy-miny) + 1;
    size_t fullsize = sizeof(double) * (_arrCount + 3); // Es werden die letzten (x,y) Indizes und ein Debugflag mit gespeichert

    if ( arrXYIntensity == nullptr || (_xmax-_xmin) != (maxx-minx) || (_ymax-_ymin) != (maxy-miny) )
    {
        createMemory( (void **)(&arrXYIntensity), fullsize, arrXYsize, false, "xyIntensity" );
#ifdef COPY_FITDATA_TO_GPU
        createMemory( (void **)(&arrDataForFit), fullsize, arrDataSize, false, "arrDataForFit" );
#ifdef CALC_FQS_IN_GPU
        createMemory( (void **)(&arrFqsForFit),  fullsize, arrFqsSize, false, "arrFqsForFit" );
#endif
#endif
        //qDebug() << "checkArrays(" << minx << maxx << miny << maxy << ")" << _arrCount << fullsize;
    }
    // save geometry
    _xmin = minx;
    _xmax = maxx;
    _ymin = miny;
    _ymax = maxy;
}

void CLASSLIB::memcleanup( void *arr )
{
    //if ( arr == nullptr )
    //    arr = arrXYIntensity;
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


#ifdef FITDATA_IN_GPU  // real func prog
bool CLASSLIB::setFitData( int sx, int sy, const double *data )
{
    _fitWidth  = sx;
    _fitHeight = sy;
    _arrFitSize = sx * sy;
    if ( data == nullptr || _arrFitSize == 0 )
    {
        if ( arrFitData != nullptr ) { memcleanup( arrFitData );   arrFitData = nullptr; }
        if ( arrFitFqs  != nullptr ) { memcleanup( arrFitFqs  );   arrFitFqs  = nullptr; }
        _fitEnabled = false;
        return false;
    }
    createMemory( (void **)(&arrFitData), _arrFitSize * sizeof(double), arrFitDSize, false, "arrFitData" );
    createMemory( (void **)(&arrFitFqs),  _arrFitSize * sizeof(double), arrFitFSize, false, "arrFitFqs" );
    _fitEnabled = arrFitData != nullptr && arrFitFqs != nullptr;
    if ( _fitEnabled )
        memcpy( arrFitData, data, arrFitDSize );

    std::cerr << "Fit-Craete: " << sx << "x" << sy << "=" << _arrFitSize << std::endl;

    return _fitEnabled;
}
#endif
