# CrystalScatter
Development is in local iffgit in Forschungszentrum Jülich.

- Base pascal program by Stephan Förster (s.foerster@fz-juelich.de).
- Conversion to C++ with Cuda by Michael Wagener (m.wagener@fz-juelich.de).

# Current state
The software is still under development but the source here is functional
but some input parameter settings may not produce the desired output.

# Documentation
The documentation provided here might have older screenshots.
Not all implemented features are documented at the moment.

Inline source documentation will be added.

# Used libraries
### mandatory
- Qt 5.14.x (no Qt6)

### optional / Windows
- Cuda 11.0
- fftw-3.3.5-dll64
- hdf5-1.12.1 (CMake build)

### optional / Linux
- cuda-rhel8-11-0-local
- fftw-devel.x86_64 (3.3.5-11.el8)
- fftw-libs.x86_64 (3.3.5-11.el8)
- hdf-devel.x86_64 (4.2.14-5.el8)
