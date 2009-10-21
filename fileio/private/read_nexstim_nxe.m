function [dat] = read_nexstim_nxe(filename, begsample, endsample, chanindx)

% READ_NEXSTIM_NXE reads specified samples from a NXE continous datafile
%
% Use as
%   [hdr] = read_nexstim_nxe(filename)
% where
%    filename        name of the datafile, including the .bdf extension
% This returns a header structure with the following elements
%   hdr.Fs           sampling frequency
%   hdr.nChans       number of channels
%   hdr.nSamples     number of samples per trial
%   hdr.nSamplesPre  number of pre-trigger samples in each trial
%   hdr.nTrials      number of trials
%   hdr.label        cell-array with labels of each channel
%
% Or use as
%   [dat] = read_nexstim_nxe(filename, begsample, endsample, chanindx)
% where
%    filename        name of the datafile, including the .nxe extension
%    begsample       index of the first sample to read
%    endsample       index of the last sample to read
%    chanindx	     index of channels to read (optional, default is all)
% This returns a Nchans X Nsamples data matrix

% Written by Vladimir Litvak based on functions provided by Nexstim
%
% Copyright (C) 2007, Vladimir Litvak
%
% $Log: read_nexstim_nxe.m,v $
% Revision 1.2  2009/02/27 12:42:03  vlalit
% Fix for a bug spotted by Miriam Klein
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.3  2007/12/17 13:03:52  roboos
% Vladimir found and fixed some bugs pertaining to the nexstim_nxe format
%
% Revision 1.2  2007/12/17 08:42:34  roboos
% fixed typo
%
% Revision 1.1  2007/12/17 08:24:28  roboos
% added support for nexstim_nxe, thanks to Vladimir
% the low-level code has not been tested by myself
% 

if nargin==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the header
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hdr.Fs          = 1450;
    hdr.nChans      = 64;
    hdr.label  = cell(64,1);
    hdr.label(1:4)= {'GATE', 'TRIG1', 'TRIG2','EOG'};
    hdr.label(5:64) = {
        'Fp1'
        'Fpz'
        'Fp2'
        'AF1'
        'AFz'
        'AF2'
        'F7'
        'F3'
        'F1'
        'Fz'
        'F2'
        'F4'
        'F8'
        'FT9'
        'FT7'
        'FC5'
        'FC3'
        'FC1'
        'FCz'
        'FC2'
        'FC4'
        'FC6'
        'FT8'
        'FT10'
        'T7'
        'C5'
        'C3'
        'C1'
        'Cz'
        'C2'
        'C4'
        'C6'
        'T8'
        'TP9'
        'TP7'
        'CP5'
        'CP3'
        'CP1'
        'CPz'
        'CP2'
        'CP4'
        'CP6'
        'TP8'
        'TP10'
        'P9'
        'P7'
        'P3'
        'P1'
        'Pz'
        'P2'
        'P4'
        'P8'
        'P10'
        'PO3'
        'POz'
        'PO4'
        'O1'
        'Oz'
        'O2'
        'Iz'};

    % it is continuous data, therefore append all records in one trial
    hdr.nTrials     = 1;

    fid=fopen(filename,'r','l');
    fseek(fid,0,'eof');
    numBytes = ftell(fid);
    hdr.nSamples = (numBytes/2)/hdr.nChans;

    hdr.nSamplesPre = 0;

    fclose(fid);

    % return the header
    dat = hdr;

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sfEEG = (1/2000) * (10/65535) * 1000000;
    sfEOG = (1/400) * (10/65535) * 1000000;
    sfTRIG = 10 * (10/65535);

    numChannels = 64;

    fid = fopen(filename,'r','l');
    fseek(fid, 2*numChannels*(begsample-1),'bof');
    data = fread(fid,[numChannels endsample-begsample+1],'short');
    fclose(fid);

    data(1:3,:) = sfTRIG.*data(1:3,:);
    data(4,:) = sfEOG.*data(4,:);
    data(5:64,:) = sfEEG.*data(5:64,:);

    if nargin<4
        chanindx = 1:numChannels;
    end

    dat = data(chanindx,:);
end

