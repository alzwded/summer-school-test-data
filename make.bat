@echo off

cl /EHa /O2 /openmp /Qpar /arch:AVX /std:c++17 genplate.cpp
cl /EHa /O2 /openmp /Qpar /arch:AVX /std:c++17 gencube.cpp

CALL :RUN plate genplate 100 2000 100 100
CALL :RUN cube16 gencube 16 2000 100 16
CALL :RUN cube32 gencube 32 2000 100 16
CALL :RUN cube gencube 100 2000 100 16

EXIT /B 0

:RUN
    IF NOT EXIST %1 mkdir %1
    pushd %1
    del *.txt *.csv *.tec
    ..\%2 %3 %4 %5 %6
    popd %1
    EXIT /B 0
