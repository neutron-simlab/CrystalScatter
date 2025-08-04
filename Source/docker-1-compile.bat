@echo on

rem *** simply compile the source files
docker run --rm  -v "%cd%:/home/user/project" g76r/qt6-builder:qt-6.8.3-debug sh -c "cd /home/user/project; rm -rf build; mkdir build; cd build; qmake ../sas_scatter2Cons.pro; make"

pause
