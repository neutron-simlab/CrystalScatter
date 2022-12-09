QT       += core
#QT       -= gui

TARGET = sas_scatter2Cons
TEMPLATE = app

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

DEFINES += CONSOLENPROG

SOURCES += \
    sc_mainCons.cpp \
    sc_calcCons.cpp \
    sc_postproc.cpp \
    sc_calc_generic.cpp \
    sc_readdata.cpp \
    sc_simplexfit2d.cpp \
    widimage.cpp

    #sc_calc_fcc.cpp \
    #sc_calc_bct.cpp \
    #sc_calc_bcc.cpp \
    #sc_calc_fd3m.cpp \
    #sc_calc_sc.cpp \
    #sc_calc_hcp.cpp \

HEADERS += \
    sc_calcCons.h \
    sc_postproc.h \
    sc_calc_generic.h \
    sc_calc_generic_gpu.h \
    sc_math.h \
    sc_readdata.h \
    sc_simplexfit2d.h \
    widimage.h

    #sc_calc_fcc.h \
    #sc_calc_fcc_gpu.h \
    #sc_calc_bct.h \
    #sc_calc_bct_gpu.h \
    #sc_calc_bcc.h \
    #sc_calc_bcc_gpu.h \
    #sc_calc_fd3m.h \
    #sc_calc_fd3m_gpu.h \
    #sc_calc_sc.h \
    #sc_calc_sc_gpu.h \
    #sc_calc_hcp.h \
    #sc_calc_hcp_gpu.h \


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/sas_scatter2/bin
!isEmpty(target.path): INSTALLS += target


# CUDA definitions
CUDA_DIR = /usr/local/cuda   # Path to cuda sdk/toolkit install
CUDA_SOURCES += sc_calc_generic_gpu.cu
                #sc_calc_fcc_gpu.cu
                #sc_calc_bct_gpu.cu
                #sc_calc_bcc_gpu.cu \
                #sc_calc_sc_gpu.cu \
                #sc_calc_hcp_gpu.cu \
                #sc_calc_fd3m_gpu.cu

!exists($$CUDA_DIR/bin/nvcc) {
    # if no cuda support is installed, treat this as a normal c++ file
    SOURCES += $$CUDA_SOURCES
    QMAKE_CFLAGS += -x c++ -std=gnu++11
} else {
    # CUDA settings <-- may change depending on your system
    SYSTEM_TYPE = 64            # '32' or '64', depending on your system
    CUDA_ARCH = compute_52          # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'
    NVCC_OPTIONS = --use_fast_math -std=c++11

    # Add the necessary libraries
    LIBS += -L$$CUDA_DIR/lib64 -lcuda -lcudart

    # Add include path
    CUDA_INC = -I$$CUDA_DIR/include

    # Configuration of the Cuda compiler
    CONFIG(debug, debug|release) {
        # Debug mode
        cuda_d.input = CUDA_SOURCES
        cuda_d.output = ${QMAKE_FILE_BASE}_cuda.o
        cuda_d.commands = $$CUDA_DIR/bin/nvcc -D_DEBUG $$NVCC_OPTIONS $$CUDA_INC $$LIBS --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
        cuda_d.dependency_type = TYPE_C
        QMAKE_EXTRA_COMPILERS += cuda_d
    }
    else {
        # Release mode
        cuda.input = CUDA_SOURCES
        cuda.output = ${QMAKE_FILE_BASE}_cuda.o
        cuda.commands = $$CUDA_DIR/bin/nvcc $$NVCC_OPTIONS $$CUDA_INC $$LIBS --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
        cuda.dependency_type = TYPE_C
        QMAKE_EXTRA_COMPILERS += cuda
    }
}


win32: {
    CONFIG(static) {
        DEFINES += USE_COMPLEX_WEB
    }
    else {
        # FFTW3 Library (Laptop)
        FFTW3_PATH = ../fftw-3.3.5-dll64
        exists($$FFTW3_PATH/fftw3.h) {
            LIBS += -L$$FFTW3_PATH -lfftw3-3
            INCLUDEPATH += $$FFTW3_PATH
            HEADERS += $$FFTW3_PATH/fftw3.h
        }
        else {
            DEFINES += USE_COMPLEX_WEB
        }
    }

    # HDF5 Library
    #HDF5_BASE = ..\..\CMake-hdf5-1.12.0
    #HDF5_GEN  = $$HDF5_BASE/build-hdf5-1.12.0-Desktop_Qt_5_14_2_MinGW_64_bit-MinSizeRel
    #exists($$HDF5_GEN/bin) {
    #    INCLUDEPATH += $$HDF5_BASE/hdf5-1.12.0/src $$HDF5_BASE/hdf5-1.12.0/c++/src $$HDF5_GEN
    #    LIBS += -L$$HDF5_GEN/bin -lhdf5_cpp -lhdf5_hl -lhdf5
    #}
    #else {
        DEFINES += NOHDF5
    #}
}
unix: { # exists(../fftw-3.3.9) {
    # FFTW3 Library (Büro-PC)
    LIBS += -L/usr/local/lib -lfftw3
    INCLUDEPATH += /usr/local/include
    HEADERS += /usr/local/include/fftw3.h

    DEFINES += NOHDF5

}

DISTFILES += \
    ConfigParams.ini \
    sc_libs_gpu.h \
    sc_memory_gpu.h \
    sc_memory_gpu.cu \
    sc_libs_gpu.cu
