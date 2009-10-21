function [backpanel, frontpanel] = read_ctf_trigger(dataset)

% READ_CTF_TRIGGER reads the STIM channel from a dataset and detects
% the trigger moments and values
%
% [backpanel, frontpanel] = read_ctf_trigger(dataset)
% 
% This returns all samples of the STIM channel, converted to backpanel
% and frontpanel trigger values. Triggers are placed at the rising flank
% of the STIM channel.
%
% Triggers should be at least 9 samples long (for 1200Hz samplerate) and
% should not overlap each other.
%
% See also READ_CTF_MEG4, READ_CTF_RES4

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_ctf_trigger.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.8  2006/02/06 09:02:35  roboos
% added support for triggers in UPPT001, USPT001, UTRG001 (thanks to Tom)
%
% Revision 1.7  2004/03/12 11:21:19  roberto
% The function outputs backpanel and frontpanel seperately. However
% it performs the upflank/kernel-convolution-correction of the
% trigger-channel on the the full 32-bit data. Usually the frontpanel
% is used to record responses while the backpanel is used for recording
% stimulus-trigger. While the timecourse of backpanel-triggers can
% be fully controlled the temporal overlap of these triggers with
% subject-responses can not. Therefore I think read_ctf_trigger should
% do the correction on the splitted data as done in the attached
% version of the function. This version works as before for nonoverlapping
% cases but can also handle backpanel-frontpanel overlap.
%
% Revision 1.6  2004/02/18 13:55:39  roberto
% fixed most significant bit, which was messed up by reading the data as signed int
%
% Revision 1.5  2004/02/12 09:55:35  roberto
% fixed bug, length of trigger was trigshift samples too short
%
% Revision 1.4  2003/10/14 12:37:17  roberto
% made teh trigshift adaptive to sampling frequency
%
% Revision 1.3  2003/09/12 09:27:11  roberto
% added comment to help about how triggers are supposed to behave
%
% Revision 1.2  2003/05/21 10:59:52  roberto
% fixed bugs in front/backpanel
%
% Revision 1.1  2003/05/19 14:40:32  roberto
% new implementation, replaces inline code in framework/preprocessing
%

[path, file, ext] = fileparts(dataset);
datafile   = fullfile(dataset, [file '.meg4']);
headerfile = fullfile(dataset, [file '.res4']);

% read the header from the raw CTF data
hdr = read_ctf_res4(headerfile);

% number of samples to shift the assesment of the trigger value
% this is needed because it takes some time for the rising flank to get to the correct value
trigshift = fix(hdr.Fs * 9/1200);

% read the stimulus channel from raw CTF data: with the new electronics
% control console used at NIH, there is now no longer any channel
% named "STIM". Instead, there are a variety of channel names. The
% electronics used at FCDC still uses the "STIM" channel.
stimindx = find(strcmp(hdr.label, 'UPPT001'));
if isempty(stimindx)
  stimindx = find(strcmp(hdr.label, 'USPT001'));
end
if isempty(stimindx)
  stimindx = find(strcmp(hdr.label, 'UTRG001'));
end
if isempty(stimindx)
  stimindx = find(strcmp(hdr.label, 'STIM'));
end
stim =  read_ctf_meg4(datafile, hdr, 1, hdr.nTrials*hdr.nSamples, stimindx);

% correct for reading stimulus channel as signed integer, whereas it should be an unsigned int
stim(find(stim<0)) = stim(find(stim<0)) + 2^32;

% split backpanel and frontpanel data
bpstim = fix(stim / 2^16);
fpstim = double(bitand(uint32(stim), 2^16-1));

% determine the precise timing of the triggers
bpupflank = [0 (diff(bpstim)>0 & bpstim(1:(end-1))==0)];
backpanel = bpupflank(1:(end-trigshift)).*bpstim((1+trigshift):end);
fpupflank = [0 (diff(fpstim)>0 & fpstim(1:(end-1))==0)];
frontpanel = fpupflank(1:(end-trigshift)).*fpstim((1+trigshift):end);

% pad with zeros to ensure same length as data
backpanel(length(stim)) = 0;
frontpanel(length(stim)) = 0;

