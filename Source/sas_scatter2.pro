QT       += core gui xml

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = sas_scatter2
TEMPLATE = app


#DEFINES += IMG_ONE_WINDOW
# Test für ein anderes Layout, damit unter Android die Images im Hauptfenster bleiben


# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
    dlgconfigautofit.cpp \
    dlgtimetests.cpp \
    myguiparam.cpp \
    sc_calc_generic.cpp \
    sc_calcgui.cpp \
    sc_main.cpp \
    sc_maingui.cpp \
    sc_postproc.cpp \
    sc_readdata.cpp \
    sc_simplexfit2d.cpp \
    widimage.cpp

HEADERS += \
    debughandler.h \
    dlgconfigautofit.h \
    dlgtimetests.h \
    myguiparam.h \
    sc_calc_generic.h \
    sc_calc_generic_gpu.h \
    sc_libs_gpu.h \
    sc_memory_gpu.h \
    sc_calcgui.h \
    sc_globalConfig.h \
    sc_maingui.h \
    sc_math.h \
    sc_postproc.h \
    sc_readdata.h \
    sc_simplexfit2d.h \
    widimage.h

FORMS += \
    dlgconfigautofit.ui \
    dlgtimetests.ui \
    sc_maingui.ui \
    widimage.ui


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/sas_scatter2/bin
!isEmpty(target.path): INSTALLS += target


# CUDA definitions
unix:{
    CUDA_DIR = /usr/local/cuda-12.2   # Path to cuda sdk/toolkit install
    CUDA_EXE = $$CUDA_DIR/bin/nvcc
}
win32:{
    #CUDA_DIR = "C:/Program\ Files/NVIDIA\ GPU\ Computing\ Toolkit/CUDA/v11.0"
    # -> nvcc fatal   : A single input file is required for a non-link phase when an outputfile is specified
    #                   -> nvcc has problems with spaces in filenames...
    CUDA_DIR = "C:/SimLab/CUDA/v12.2"  # link to above dir

    # 'X' hinten zum Ausblenden dieser Funktion.
    CUDA_EXE = "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.2/bin/nvcc.exeX"
}
# nvcc can work with multiple input files but then they must be linked with a special cuda linker
#  and this is not so easy to include in this Qt-Project. So only a single input is used and the
#  other .cu files are included there.
CUDA_SOURCES += sc_calc_generic_gpu.cu

!exists($$CUDA_EXE) {
    # if no cuda support is installed, treat this as a normal c++ file
    SOURCES += $$CUDA_SOURCES
    QMAKE_CFLAGS += -x c++ -std=gnu++11
} else {
    # CUDA settings <-- may change depending on your system
    SYSTEM_NAME = Win64         # Depending on your system either 'Win32', 'x64', or 'Win64'
    SYSTEM_TYPE = 64            # '32' or '64', depending on your system
    CUDA_ARCH = compute_52      # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'
    NVCC_OPTIONS = --use_fast_math -std=c++11  # todo: -warning 550 ...
    #NVCC_OPTIONS += --ptxas-options=--verbose  # Specify options directly to ptxas, the PTX optimizing assembler.

# The following library conflicts with something in Cuda
    QMAKE_LFLAGS_RELEASE = /NODEFAULTLIB:msvcrt.lib         # TODO: Unter Linux: file not found...
    QMAKE_LFLAGS_DEBUG   = /NODEFAULTLIB:msvcrtd.lib

    message("Use CUDA")

    # Add the necessary libraries and Add include path
    win32:{
        # c++  cpp  g++  gcc  x86_64-w64-mingw32-c++  x86_64-w64-mingw32-g++  x86_64-w64-mingw32-gcc
        NVCC_OPTIONS += --dont-use-profile \
            --allow-unsupported-compiler \
            -ccbin=C:/Qt/Tools/mingw810_64/bin/x86_64-w64-mingw32-g++.exe \
            --forward-unknown-to-host-compiler \
            --forward-unknown-to-host-linker \
            --drive-prefix="/" --verbose
        LIBS += "-l$$CUDA_DIR/lib/x64/cuda"
        LIBS += "-l$$CUDA_DIR/lib/x64/cudart"
        CUDA_INC = "$$CUDA_DIR/include"
    } else {
        LIBS += -L$$CUDA_DIR/lib64 -lcuda -lcudart
        CUDA_INC = -I$$CUDA_DIR/include
        CUDA_INC += -I/usr/include/qt5 -I/usr/include/qt5/QtCore
    }

    # Configuration of the Cuda compiler
    CONFIG(debug, debug|release) {
        # Debug mode
        cuda_d.input = CUDA_SOURCES
        cuda_d.output = ${QMAKE_FILE_BASE}_cuda.o
        cuda_d.commands = $$CUDA_EXE -D_DEBUG $$NVCC_OPTIONS -I $$CUDA_INC --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
        cuda_d.dependency_type = TYPE_C
        QMAKE_EXTRA_COMPILERS += cuda_d
    }
    else {
        # Release mode
        cuda.input = CUDA_SOURCES
        cuda.output = ${QMAKE_FILE_BASE}_cuda.o
        cuda.commands = $$CUDA_EXE $$NVCC_OPTIONS -I $$CUDA_INC $$LIBS --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
        cuda.dependency_type = TYPE_C
        QMAKE_EXTRA_COMPILERS += cuda
    }
}

win32: {
    CONFIG(static) {
        message("No FFTW-Lib (static)")
        DEFINES += USE_COMPLEX_WEB
    }
    else {
        # FFTW3 Library (Laptop)
        FFTW3_PATH = ../fftw-3.3.5-dll64
        exists($$FFTW3_PATH/fftw3.h) {
            LIBS += -L$$FFTW3_PATH -lfftw3-3
            INCLUDEPATH += $$FFTW3_PATH
            HEADERS += $$FFTW3_PATH/fftw3.h
            message("Use FFTW-Lib in " $$FFTW3_PATH)
        }
        else {
            DEFINES += USE_COMPLEX_WEB
            message("No FFTW-Lib (not found) " $$FFTW3_PATH)
        }
    }

    # HDF5 Library
    HDF5_BASE = ..\..\CMake-hdf5-1.12.1
    HDF5_SRC  = $$HDF5_BASE/hdf5-1.12.1
    CONFIG(static) {
        # C:\SimLab\CMake-hdf5-1.12.1\build-ZLib-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\bin\libzlib_D.a
        # C:\Windows\System32>echo %ZLIB_DIR%
        # C:\SimLab\CMake-hdf5-1.12.1\build-ZLib-Desktop_Qt_5_14_2_MinGW_64_bit-Debug

        HDF5_GEN  = $$HDF5_BASE/build-hdf5-1.12.1-Desktop_Qt_5_14_2_MinGW_64_bit-MinSizeRel
        HDF5_LIBS = -L$$HDF5_GEN/bin -lhdf5_cpp -lhdf5_hl_cpp -lhdf5 -lhdf5_hl \
                    -L$$(ZLIB_DIR)/bin -lzlib_D
    }
    else {
        HDF5_GEN  = $$HDF5_BASE/build-hdf5-1.12.1-Desktop_Qt_5_14_2_MinGW_64_bit-Debug
        HDF5_LIBS = -L$$HDF5_GEN/bin -lhdf5_cpp_D -lhdf5_hl_cpp_D -lhdf5_D -lhdf5_hl_D
    }
    exists($$HDF5_GEN/bin) {
        INCLUDEPATH += $$HDF5_SRC/src $$HDF5_SRC/c++/src $$HDF5_GEN $$HDF5_GEN/src
        LIBS += $$HDF5_LIBS
        SOURCES += dlghdfselection.cpp
        HEADERS += dlghdfselection.h
    }
    else {
        DEFINES += NOHDF5
        message("No HDF5 found in " $$HDF5_GEN)
    }
}

unix: {     # no static !

    # FFTW3 Library
    FFTW3_PATH = /usr/local/include
    exists($$FFTW3_PATH/fftw3.h) {
        LIBS += -L/usr/local/lib -lfftw3
        INCLUDEPATH += /usr/local/include
        HEADERS += /usr/local/include/fftw3.h
        message("Use FFTW-Lib in " $$FFTW3_PATH)
    }
    else {
        DEFINES += USE_COMPLEX_WEB
        message("No FFTW-Lib (not found) " $$FFTW3_PATH)
    }

    # HDF5 Library
    HDF5_BASE = /usr/local/hdf5
    exists($$HDF5_BASE/bin) {
        INCLUDEPATH += $$HDF5_BASE/include
        LIBS += -L$$HDF5_BASE/lib -lhdf5_cpp -lhdf5_hl_cpp -lhdf5 -lhdf5_hl
        SOURCES += dlghdfselection.cpp
        HEADERS += dlghdfselection.h
    }
    else {
        DEFINES += NOHDF5
        message("No HDF5 found in " $$HDF5_BASE)
    }

}

DISTFILES += \
    "../Pascal-Sourcecodes/20221004 - crystal3d1.pas" \
    "../Pascal-Sourcecodes/20230420 - crystal3d1.pas" \
    sc_memory_gpu.cu \
    sc_libs_gpu.cu
    #sc_libs_gpu.h \
    #sc_memory_gpu.h \

    #"../20210616 - Neue Routinen/FHKL_routine_full.txt" \
    #"../20210616 - Neue Routinen/HKL_routine_basic.txt" \
    #"../20210616 - Neue Routinen/coordinate_rotation_routine_full.txt" \
    #"../20210616 - Neue Routinen/formfactor_coefficient_routine_full.txt" \
    #"../20210616 - Neue Routinen/generic_routine_basic.txt" \
    #"../20210616 - Neue Routinen/generic_routine_full.txt" \
    #"../20210616 - Neue Routinen/main_routine_basic.txt" \
    #"../20210616 - Neue Routinen/main_routine_full.txt" \
    #"../20210616 - Neue Routinen/numerical_integration_routines_full.txt" \
    #"../20210616 - Neue Routinen/cpu_gpu_lattice.txt" \
