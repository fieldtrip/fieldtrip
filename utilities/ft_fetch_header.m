function [hdr] = ft_fetch_header(data)

% FT_FETCH_HEADER mimics the behavior of FT_READ_HEADER, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   hdr = ft_fetch_header(data)
%
% See also FT_READ_HEADER, FT_FETCH_DATA, FT_FETCH_EVENT

% Copyright (C) 2008, Esther Meeuwissen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');

trlnum = length(data.trial);
trllen = zeros(trlnum,1);
for trllop=1:trlnum
  trllen(trllop) = size(data.trial{trllop},2);
end

if isfield(data, 'sampleinfo')
  % construct the trial definition according to thew samples from the original data file
  trl = data.sampleinfo;
  trl(:,3) = 0;
else
  % construct the trial definition as if it is a continuous piece of data
  trl = [1 sum(trllen) 0];
end

% fill in some header details
hdr.Fs     = data.fsample;
hdr.label  = data.label(:);
hdr.nChans = numel(data.label);

% fill in the channel type
if isfield(data, 'chantype')
  hdr.chantype = data.chantype(:);
elseif isfield(data, 'hdr') && isfield(data.hdr, 'chantype')
  % keep them ordered according to the FieldTrip data structure, which might differ from the original header
  [datindx, hdrindx] = match_str(data.label, data.hdr.label);
  hdr.chantype = repmat({'unknown'}, hdr.nChans, 1);
  hdr.chantype(datindx) = data.hdr.chantype(hdrindx);
end

% fill in the channel unit
if isfield(data, 'chanunit')
  hdr.chanunit = data.chanunit(:);
elseif isfield(data, 'hdr') && isfield(data.hdr, 'chanunit')
  % keep them ordered according to the FieldTrip data structure, which might differ from the original header
  [datindx, hdrindx] = match_str(data.label, data.hdr.label);
  hdr.chanunit = repmat({'unknown'}, hdr.nChans, 1);
  hdr.chanunit(datindx) = data.hdr.chanunit(hdrindx);
end

% try to determine them on the basis of heuristics, when already present they will stay the same
hdr.chantype = ft_chantype(hdr);
hdr.chanunit = ft_chanunit(hdr);

% determine hdr.nSamples, hdr.nSamplesPre, hdr.nTrials
% always pretend that it is continuous data
hdr.nSamples    = max(trl(:,2));
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;

% retrieve the gradiometer and/or electrode and/or optode information
if isfield(data, 'grad')
  hdr.grad = data.grad;
elseif isfield(data, 'hdr') && isfield(data.hdr, 'grad')
  hdr.grad = data.hdr.grad;
end
if isfield(data, 'elec')
  hdr.elec = data.elec;
elseif isfield(data, 'hdr') && isfield(data.hdr, 'elec')
  hdr.elec = data.hdr.elec;
end
if isfield(data, 'opto')
  hdr.opto = data.opto;
elseif isfield(data, 'hdr') && isfield(data.hdr, 'opto')
  hdr.opto = data.hdr.opto;
end

% retrieve the synchronization information
if isfield(data, 'hdr') && isfield(data.hdr, 'FirstTimeStamp')
  hdr.FirstTimeStamp = data.hdr.FirstTimeStamp;
end
if isfield(data, 'hdr') && isfield(data.hdr, 'TimeStampPerSample')
  hdr.TimeStampPerSample = data.hdr.TimeStampPerSample;
end
