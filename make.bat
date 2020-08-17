@echo off

IF [%1] == [clean] GOTO :CLEAN

cl /EHa /O2 /openmp /Qpar /arch:AVX /std:c++17 genplate.cpp || EXIT /B 1
cl /EHa /O2 /openmp /Qpar /arch:AVX /std:c++17 gencube.cpp || EXIT /B 1

CALL :RUN plate genplate 100 2000 100 100
CALL :RUN cube16 gencube 16 2000 100 16
CALL :RUN cube32 gencube 32 2000 100 16
CALL :RUN cube gencube 100 2000 100 16

EXIT /B 0

:RUN
    IF NOT EXIST %1 mkdir %1
    pushd %1
    del /Q *.txt *.csv *.tec
    ..\%2 %3 %4 %5 %6
    popd
    EXIT /B 0

:CLEAN
    del /q /s *.exe *.obj *.pdb *.ilk plate cube16 cube32 cube
    rmdir /q /s plate cube16 cube32 cube
    EXIT /B 0
