function [hdr] = fetch_header(data)

% FETCH_HEADER mimics the behaviour of READ_HEADER, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   [event] = fetch_header(data)
%
% See also READ_HEADER, FETCH_DATA, FETCH_EVENT

% Copyright (C) 2008, Esther Meeuwissen
%
% $Log: fetch_header.m,v $
% Revision 1.2  2009/10/01 12:36:01  jansch
% workaround if input data has been resampled, so that the comparison with
% trl matrix does not lead to a crash
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.2  2008/09/29 21:12:39  roboos
% cleaned up the code from Esther, added copyrights, updated documentation
%

% check whether input is data
data = checkdata(data, 'datatype', 'raw');

% get trial definition according to original data file
trl    = findcfg(data.cfg, 'trl');
trlnum = length(data.trial);
trllen = zeros(trlnum,1);
for trllop=1:trlnum
  trllen(trllop) = size(data.trial{trllop},2);
end

% check whether data.trial is consistent with trl
if size(trl,1)~=length(data.trial)
  error('trial definition is not internally consistent')
elseif any(trllen~=(trl(:,2)-trl(:,1)+1)) && ~isempty(findcfg(data.cfg, 'resamplefs')) && ~isempty(findcfg(data.cfg,'resampletrl')),
  warning('the data have been resampled along the way, the trl-definition is in the original sampling rate, attempt to adjust for this may introduce some timing inaccuracies');
  trlold = trl;
  trl    = findcfg(data.cfg, 'resampletrl');   
end

%this has to be done again
if any(trllen~=(trl(:,2)-trl(:,1)+1))
  error('trial definition is not internally consistent')
end

% fill in hdr.nChans 
hdr.nChans = length(data.label);

% fill in hdr.label 
hdr.label = data.label;

% fill in hdr.Fs (sample frequency)
hdr.Fs = data.fsample;

% determine hdr.nSamples, hdr.nSamplesPre, hdr.nTrials
% always pretend that it is continuous data
hdr.nSamples    = max(trl(:,2));
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;

% fill in hdr.grad or hdr.elec
if isfield(data, 'grad')
  hdr.grad=data.grad;
elseif isfield(data, 'elec')
  hdr.elec=data.elec;
end

