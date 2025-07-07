# +++++ Settings of all path variables used for the specialized libraries:
include(sas_scatter2.pri)
# ----- User configuration ends here.

FFTW3_PATH = ""
FFTW3_LIBS = ""


QT       += core
#QT       -= gui  wird für QImage und QColor gebraucht um die Images zu erzeugen

TARGET = sas_scatter2Cons
TEMPLATE = app


# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

greaterThan(QT_MAJOR_VERSION, 5) {
    CONFIG += c++17
    QMAKE_CFLAGS += -x c++ -std=gnu++17
    #message("C++17")
} else {
    CONFIG += c++11
    QMAKE_CFLAGS += -x c++ -std=gnu++11
    #message("C++11")
}

DEFINES += CONSOLENPROG
DEFINES += NOQWT

SOURCES += \
    sc_mainCons.cpp \
    sc_calcCons.cpp \
    sc_postproc.cpp \
    sc_calc_generic.cpp \
    sc_calc_generic_cpu.cpp \
    sc_readdata.cpp \
    sc_simplexfit2d.cpp \
    widimage.cpp

HEADERS += \
    sc_calcCons.h \
    sc_postproc.h \
    sc_calc_generic.h \
    sc_calc_generic_gpu.h \
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
    sc_math.h \
    sc_readdata.h \
    sc_simplexfit2d.h \
    widimage.h


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/sas_scatter2/bin
!isEmpty(target.path): INSTALLS += target


# nvcc can work with multiple input files but then they must be linked with a special cuda linker
#  and this is not so easy to include in this Qt-Project. So only a single input is used and the
#  other .cu files are included there.
CUDA_SOURCES += sc_calc_generic_gpu.cu

!exists($$CUDA_EXE) {
    # if no cuda support is installed, treat this as a normal c++ file
    SOURCES += $$CUDA_SOURCES
    CONFIG += -nocudalib

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
    #SYSTEM_NAME = Win64         # Depending on your system either 'Win32', 'x64', or 'Win64'
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
            -ccbin=C:/Qt/Tools/mingw810_64/bin/x86_64-w64-mingw32-c++.exe \
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
        DEFINES += USE_COMPLEX_WEB
    }
    else {
        # FFTW3 Library (Laptop)
        exists($$FFTW3_PATH/fftw3.h) {
            LIBS += -L$$FFTW3_LIBS -lfftw3-3
            INCLUDEPATH += $$FFTW3_PATH
            HEADERS += $$FFTW3_PATH/fftw3.h
        }
        else {
            DEFINES += USE_COMPLEX_WEB
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
    }
    else {
        DEFINES += NOHDF5
        message("No HDF5 found in " $$HDF5_BASE)
    }

}
