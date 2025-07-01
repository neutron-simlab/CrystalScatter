# +++++ Settings of all path variables used for the specialized libraries:

unix:{
    CUDA_DIR = /usr/local/cuda-12.2   # Path to cuda sdk/toolkit install
    CUDA_EXE = $$CUDA_DIR/bin/nvcc
    CUDA_ARCH = sm_75  # native #  compute_52      # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'

    FFTW3_PATH = /usr/local/include
    FFTW3_LIBS = /usr/local/lib

    HDF5_BASE = /usr/local/hdf5

    QWT_BASE = /usr/local/qwt
}

win32:{
    CUDA_DIR = "C:/SimLab/CUDA/v12.2"  # link to "C:/Program\ Files/NVIDIA\ GPU\ Computing\ Toolkit/CUDA/v12.2"
    CUDA_EXE = "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.2/bin/nvcc.exeX"  # append X to disable it
    CUDA_ARCH = sm_75  # native #  compute_52      # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'

    FFTW3_PATH = ../fftw-3.3.5-dll64    # not used in static configuration
    FFTW3_LIBS = ../fftw-3.3.5-dll64

    HDF5_VERSIONS = hdf5-1.12.1 hdf5_1_14_5 # list of possible installed versions
    # resultant path: ../../CMake-{version}/{version}/build-{version}-Dekstop_Qt_{qtversion}_MinGW_64_bin-[Debug|MinSizeRel]
    # ZLIB_DIR must be set globally:  ZLIB_DIR = $$HDF5_BASE/build-ZLib-Desktop_Qt_$${QTVERSSTR}_MinGW_64_bit-Debug

    QWT_BASE = "C:/qwt-6.3.0"
}

macx|macos:{
    CUDA_DIR = /usr/local/cuda-12.2   # Path to cuda sdk/toolkit install
    CUDA_EXE = $$CUDA_DIR/bin/nvcc
    CUDA_ARCH = sm_75  # native #  compute_52      # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'

    FFTW3_PATH = /usr/local/include
    FFTW3_LIBS = /usr/local/lib

    HDF5_BASE = /usr/local/hdf5

    QWT_BASE = /usr/local/qwt
}
