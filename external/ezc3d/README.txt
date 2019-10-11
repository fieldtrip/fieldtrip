This directory contains the compiled mex files for https://github.com/pyomeca/ezc3d version 1.1.0.

# Linux

I had to use the CMAKE_POSITION_INDEPENDENT_CODE option to revolve an issue as explained here https://stackoverflow.com/questions/38296756/what-is-the-idiomatic-way-in-cmake-to-add-the-fpic-compiler-option/38297422

roboos@mentat001> module load cmake/3.11.3
roboos@mentat001> cmake -fPIC -D BUILD_SHARED_LIBS=FALSE -D BINDER_MATLAB=ON -D Matlab_ROOT_DIR=/opt/matlab/R2019a -D CMAKE_POSITION_INDEPENDENT_CODE=ON .
roboos@mentat001> make

# MacOS

roboos@mac036> cmake -D BUILD_SHARED_LIBS=FALSE -D BINDER_MATLAB=ON -D Matlab_ROOT_DIR=/Applications/MATLAB_R2019a.app/ .
roboos@mac036> make

# Windows

I do not have access to a windows machine at the moment, so cannot compile the binaries.

