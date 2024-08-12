# +++++ Settings of all path variables used for the specialized libraries:

unix:{
    CUDA_DIR = /usr/local/cuda-12.2   # Path to cuda sdk/toolkit install
    CUDA_EXE = $$CUDA_DIR/bin/nvcc

    FFTW3_PATH = /usr/local/include
    FFTW3_LIBS = /usr/local/lib

    HDF5_BASE = /usr/local/hdf5
}

win32:{
    CUDA_DIR = "C:/SimLab/CUDA/v12.2"  # link to "C:/Program\ Files/NVIDIA\ GPU\ Computing\ Toolkit/CUDA/v12.2"
    CUDA_EXE = "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.2/bin/nvcc.exeX"  # append X to disable ist

    FFTW3_PATH = ../fftw-3.3.5-dll64    # not used in static configuration
    FFTW3_LIBS = ../fftw-3.3.5-dll64

    HDF5_BASE = ../../CMake-hdf5-1.12.1
    HDF5_SRC  = $$HDF5_BASE/hdf5-1.12.1
    CONFIG(static) {
        HDF5_GEN  = $$HDF5_BASE/build-hdf5-1.12.1-Desktop_Qt_5_14_2_MinGW_64_bit-MinSizeRel
    }
    else {
        HDF5_GEN  = $$HDF5_BASE/build-hdf5-1.12.1-Desktop_Qt_5_14_2_MinGW_64_bit-Debug
    }
    # ZLIB_DIR must be set
}

# ----- User configuration ends here.


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
    sc_calc_generic_cpu.cpp \
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
    sc_calcgui.h \
    sc_globalConfig.h \
    sc_gpu_generic.h \
    sc_gpu_pCub_lNon_oIso.h \
    sc_gpu_pCyl_lNon_oGau.h \
    sc_gpu_pSph_lNon_oIso.h \
    sc_gpu_pSph_pGau_lFCC_oGau.h \
    sc_gpu_pCyl_lNon_oZDi.h \
    sc_gpu_pCyl_pGau_lBCT_oZDi.h \
    sc_lib_formfq_partCylinder.h \
    sc_lib_formfq_partDisk.h \
    sc_lib_formfq_partSphere.h \
    sc_lib_formpq_partCube.h \
    sc_lib_formpq_partCylinder.h \
    sc_lib_formpq_partDisk.h \
    sc_lib_formpq_partEllips.h \
    sc_lib_formpq_partSphere.h \
    sc_lib_formpq_partSuperEllips.h \
    sc_lib_formpq_partSuperball.h \
    sc_lib_formpq_partTriaxEllips.h \
    sc_maingui.h \
    sc_math.h \
    sc_postproc.h \
    sc_readdata.h \
    sc_simplexfit2d.h \
    widimage.h
    #sc_lib_formfq_partCube.h \

FORMS += \
    dlgconfigautofit.ui \
    dlgtimetests.ui \
    sc_maingui.ui \
    widimage.ui


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/sas_scatter2/bin
!isEmpty(target.path): INSTALLS += target


# OpenCL definitions - not used, only for future tests
#win32:{
#    DEFINES += CL_HPP_TARGET_OPENCL_VERSION=300
#    OCL_DIR = "../../OpenCL-SDK/install"
#}


# nvcc can work with multiple input files but then they must be linked with a special cuda linker
#  and this is not so easy to include in this Qt-Project. So only a single input is used and the
#  other .cu files are included there.
CUDA_SOURCES += sc_calc_generic_gpu.cu

!exists($$CUDA_EXE) {
    # if no cuda support is installed, treat this as a normal c++ file
    SOURCES += $$CUDA_SOURCES
    QMAKE_CFLAGS += -x c++ -std=gnu++11

    #exists($$OCL_DIR) {
    #    # OpenCL usage
    #    INCLUDEPATH += $$OCL_DIR/include
    #    HEADERS += $$OCL_DIR/include/CL/opencl.hpp
    #    LIBS += -L$$OCL_DIR/lib -lOpenCL
    #    DEFINES += __OPENCL__
    #    message("Use OpenCL in " $$OCL_DIR)
    #}
} else {
    # CUDA settings <-- may change depending on your system
    SYSTEM_NAME = Win64         # Depending on your system either 'Win32', 'x64', or 'Win64'
    CUDA_ARCH = sm_75  # native #  compute_52      # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'
    NVCC_OPTIONS = --use_fast_math -std=c++11 # --verbose --keep
    #NVCC_OPTIONS += --ptxas-options=--verbose  # Specify options directly to ptxas, the PTX optimizing assembler.
    NVCC_OPTIONS += -diag-suppress 550         # variable "xx" set but never used

    win32:{
        # The following library conflicts with something in Cuda
        QMAKE_LFLAGS_RELEASE = /NODEFAULTLIB:msvcrt.lib
        QMAKE_LFLAGS_DEBUG   = /NODEFAULTLIB:msvcrtd.lib
    }

    message("Use CUDA")

    # Add the necessary libraries and Add include path
    win32:{
        # c++  cpp  g++  gcc  x86_64-w64-mingw32-c++  x86_64-w64-mingw32-g++  x86_64-w64-mingw32-gcc
        NVCC_OPTIONS += --dont-use-profile \
            --allow-unsupported-compiler --use-local-env \
            -ccbin=C:/Qt/Tools/mingw810_64/bin/x86_64-w64-mingw32-c++.exe \ # not working on my Laptop (Win10) :(
            --forward-unknown-to-host-compiler \
            --forward-unknown-to-host-linker \
            --verbose --drive-prefix="/"
        LIBS += "-l$$CUDA_DIR/lib/x64/cuda"
        LIBS += "-l$$CUDA_DIR/lib/x64/cudart"
        CUDA_INC = "$$CUDA_DIR/include"
        #--compiler-bindir <path>                        (-ccbin)
        #        Specify the directory in which the host compiler executable resides.  The
        #        host compiler executable name can be also specified to ensure that the correct
        #        host compiler is selected.  In addition, driver prefix options ('--input-drive-prefix',
        #        '--dependency-drive-prefix', or '--drive-prefix') may need to be specified,
        #        if nvcc is executed in a Cygwin shell or a MinGW shell on Windows.
    } else {
        LIBS += -L$$CUDA_DIR/lib64 -lcuda -lcudart
        CUDA_INC = -I$$CUDA_DIR/include
        #CUDA_INC += -I/usr/include/qt5 -I/usr/include/qt5/QtCore
    }

    # Configuration of the Cuda compiler
    #   --machine 64  is the only allowed and default
    CONFIG(debug, debug|release) {
        # Debug mode
        cuda_d.input = CUDA_SOURCES
        cuda_d.output = ${QMAKE_FILE_BASE}_cuda.o
        cuda_d.commands = $$CUDA_EXE -D_DEBUG $$NVCC_OPTIONS -I $$CUDA_INC -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
        cuda_d.dependency_type = TYPE_C
        QMAKE_EXTRA_COMPILERS += cuda_d
    }
    else {
        # Release mode
        cuda.input = CUDA_SOURCES
        cuda.output = ${QMAKE_FILE_BASE}_cuda.o
        cuda.commands = $$CUDA_EXE $$NVCC_OPTIONS -I $$CUDA_INC $$LIBS -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
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
        exists($$FFTW3_PATH/fftw3.h) {
            LIBS += -L$$FFTW3_LIBS -lfftw3-3
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
    CONFIG(static) {
        # C:\SimLab\CMake-hdf5-1.12.1\build-ZLib-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\bin\libzlib_D.a
        # C:\Windows\System32>echo %ZLIB_DIR%
        # C:\SimLab\CMake-hdf5-1.12.1\build-ZLib-Desktop_Qt_5_14_2_MinGW_64_bit-Debug
        HDF5_LIBS = -L$$HDF5_GEN/bin -lhdf5_cpp -lhdf5_hl_cpp -lhdf5 -lhdf5_hl \
                    -L$$(ZLIB_DIR)/bin -lzlib_D
    }
    else {
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
    exists($$FFTW3_PATH/fftw3.h) {
        LIBS += -L$$FFTW3_LIBS -lfftw3
        INCLUDEPATH += $$FFTW3_PATH
        HEADERS += $$FFTW3_PATH/fftw3.h
        message("Use FFTW-Lib in " $$FFTW3_PATH)
    }
    else {
        DEFINES += USE_COMPLEX_WEB
        message("No FFTW-Lib (not found) " $$FFTW3_PATH)
    }

    # HDF5 Library
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

win32:{
    # Macht den Zugriff im QtCreator einfacher
    DISTFILES += \
        "../Pascal-Sourcecodes/20230420 - crystal3d1.pas" \
        "../Pascal-Sourcecodes/20220730 - Upq1.pas" \
        "../Pascal-Sourcecodes/20240301 - crystal3d1.pas"

        #"../Pascal-Sourcecodes/20221004 - crystal3d1.pas" \
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
}
