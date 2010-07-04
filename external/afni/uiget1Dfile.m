function [ff, pp, fi] = uiget1Dfile()
% open a gui for 1D file selection
% returns the same output as uigetfile

[ff, pp, fi] = uigetfile( {'*.1D;*.1D.dset' , ...
                         'AFNI 1D Files (*.1D, *.1D.dset)';},...
                         'Select 1D file' );
return
