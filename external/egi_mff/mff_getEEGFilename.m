%% mff_getEEGFilename.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012
%  Copyright 2012 EGI. All rights reserved.
%  Support routine for MFF Matlab code. Not intended to be called directly.
%  
%  Gets the filename of the bin file with the EEG in it. 
%  todo?: Replace with new call that returns EEG bin file. 
%%
function EEGFile = mff_getEEGFilename(mfffileObj)
binfiles = mfffileObj.getSignalResourceList(false);
EEGBinInd = 0; 
%EEGFile = binfiles.elementAt(EEGBinInd);
EEGFile = binfiles.get(EEGBinInd);
