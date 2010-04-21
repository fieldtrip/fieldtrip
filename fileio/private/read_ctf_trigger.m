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
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

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

