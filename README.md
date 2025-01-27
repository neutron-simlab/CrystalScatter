# CrystalScatter
Development is in local iffgit in Forschungszentrum Jülich.

- Base pascal program by Stephan Förster (s.foerster@fz-juelich.de).
- Conversion to C++ with Cuda by Michael Wagener (m.wagener@fz-juelich.de).

# Current state
The software is still under development but the source here is functional.
Some input parameter settings may not produce the desired output.

There is no Cuda Windows library available for the currently used compiler.
Therefore no GPU support for the windows platform available. The program
will use only all found CPU-Cores for the maximum number of threads.

# Documentation
The documentation provided here might have older screenshots. Except the
SasScatter2-Internals.tex and .pdf, this document is up to date and includes
code examples in C++ and Python to calculate a sphere.

Not all implemented features are documented in the SasScatter2-Userdoc.tex
and .pdf at the moment. Inline source documentation will be added.

The document HyperSeries_SI_rev.pdf is the Supporting Information for the
Scientific Report (https://www.nature.com/articles/s41598-023-27558-8).

# License
This software is developed under Open Source License.
© (2023-2025) Forschungszentrum Jülich GmbH, JCNS-1.

# Used libraries
The given version numbers are used during the development and are working.
Other library versions may work too. To configure the used path of the
external libraries edit the sas_scatter2.pri file in the source directory.

### mandatory
- Qt 5.14.x (first development)
- Qt 6.7.2 (current version)

### optional / Windows
- Qwt-6.3.0
- fftw-3.3.5-dll64
- hdf5-1.12.1 (CMake build)

### optional / Linux
- cuda-rhel8-11-0-local
- qwt-qt5-devel.x86_64 (6.1.5-5.el8)
- fftw-devel.x86_64 (3.3.5-11.el8)
- fftw-libs.x86_64 (3.3.5-11.el8)
- hdf-devel.x86_64 (4.2.14-5.el8)
- or download and install manually from hdf5 website
  (checked with 1.12.1 and 1.14.4-2)

If you want to use a cuda version newer than 11.7 then you have
to set the environment variable CUDA_MODULE_LOADING=EAGER before
compile and launch to disable the lazy loading feature.
