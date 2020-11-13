This directory contains the compiled mex files for https://github.com/pyomeca/ezc3d version 1.1.0.

# Linux

I had to use the CMAKE_POSITION_INDEPENDENT_CODE option to revolve an issue as explained here https://stackoverflow.com/questions/38296756/what-is-the-idiomatic-way-in-cmake-to-add-the-fpic-compiler-option/38297422

Initially I compiled the mex files with R2019a; however these would not work with R016b. Compiling with R2016b and R2017b failed due to an error with mxGetDoubles. I recompiled the linux mex files with R2018b. These work up to R2020b, but not with older versions.

roboos@mentat001> module load cmake/3.11.3
roboos@mentat001> cmake -fPIC -D BUILD_SHARED_LIBS=FALSE -D BINDER_MATLAB=ON -D Matlab_ROOT_DIR=/opt/matlab/R2019a -D CMAKE_POSITION_INDEPENDENT_CODE=ON .
roboos@mentat001> make

# MacOS

roboos@mac036> cmake -D BUILD_SHARED_LIBS=FALSE -D BINDER_MATLAB=ON -D Matlab_ROOT_DIR=/Applications/MATLAB_R2019a.app/ .
roboos@mac036> make

# Windows

I executed the following in the "x64 Native Tools Command Prompt for VS 2017"

cmake -D CMAKE_BUILD_TYPE=Release -D BUILD_SHARED_LIBS=TRUE -D BINDER_MATLAB=ON -D Matlab_ROOT_DIR="c:\Program Files\MATLAB\R2019a" -D Matlab_MEX_LIBRARY="C:\Program Files\MATLAB\R2019a\extern\lib\win64\microsoft\libmex.lib" -D Matlab_MX_LIBRARY="C:\Program Files\MATLAB\R2019a\extern\lib\win64\microsoft\libmx.lib" -G"NMake Makefiles" .
nmake

Subsequently I had to copy the ezc3d.dll file from the top level directory to the location with the mex files. Compiling without shared libs did not work.
