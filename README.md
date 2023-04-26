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
SasScatter2-Internals.tex and .pdf, this document is up to date.

Not all implemented features are documented in the SasScatter2-Userdoc.tex
and .pdf at the moment. Inline source documentation will be added.

The document HyperSeries_SI_rev.pdf is the Supporting Information for the
Scientific Report (https://www.nature.com/articles/s41598-023-27558-8).

# License
This software is developed under Open Source License.
© (2023) Forschungszentrum Jülich GmbH, JCNS-1.

# Used libraries
The given version numbers are used during the development and are working.
Other library versions may work too.

### mandatory
- Qt 5.14.x (no Qt6)

### optional / Windows
- fftw-3.3.5-dll64
- hdf5-1.12.1 (CMake build)

### optional / Linux
- cuda-rhel8-11-0-local
- fftw-devel.x86_64 (3.3.5-11.el8)
- fftw-libs.x86_64 (3.3.5-11.el8)
- hdf-devel.x86_64 (4.2.14-5.el8)
