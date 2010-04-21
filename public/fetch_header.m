function [hdr] = fetch_header(data)

% FETCH_HEADER mimics the behaviour of READ_HEADER, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   [hdr] = fetch_header(data)
%
% See also READ_HEADER, FETCH_DATA, FETCH_EVENT

% Copyright (C) 2008, Esther Meeuwissen
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

