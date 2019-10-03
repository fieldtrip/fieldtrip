disp('Compiling mex files');

current_path = cd;

addpath([current_path,'/data/synth'])

cd private

mex gc_aux_mex.cpp

cd(current_path)
