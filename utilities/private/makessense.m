function status = makessense(data, field)

% MAKESSENSE determines whether a some specific fields in a FieldTrip data structure
% make sense.
%
% Use as
%   status = makessense(data, field)
%
% See also GETDIMORD, GETDIMSIZ, GETDATFIELD

% Copyright (C) 2018-2019, Robert Oostenveld
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

%% determine the size of the data
nrpt = nan;
nsmp = nan;

% it should check all possible data fields that are well-behaved
datfield = setdiff(fieldnames(data), ignorefields('makessense'));
datfield = datfield(~endsWith(datfield, 'dimord'));

for i=1:numel(datfield)
  if isfield(data, datfield{i})
    dimord = getdimord(data, datfield{i});
    dimtok = tokenize(dimord, '_');
    dimsiz = getdimsiz(data, datfield{i}, numel(dimtok));
    if strcmp(dimord, '{rpt}_chan_time') || strcmp(dimord, '{subj}_chan_time')
      nrpt  = dimsiz(1);
      nsmp  = cellfun(@(x)size(x,2), data.(datfield{i})(:)); % it should be a column array
      break;
    elseif all(ismember({'rpt', 'time'}, dimtok))
      nrpt = dimsiz(strcmp(dimtok, 'rpt'));
      nsmp = repmat(dimsiz(strcmp(dimtok, 'time')), [nrpt 1]);
      break;
    elseif all(ismember({'subj', 'time'}, dimtok))
      nrpt = dimsiz(strcmp(dimtok, 'subj'));
      nsmp = repmat(dimsiz(strcmp(dimtok, 'time')), [nrpt 1]);
      break;
    elseif all(ismember({'rpt', 'freq'}, dimtok))
      nrpt = dimsiz(strcmp(dimtok, 'rpt'));
      nsmp = nan; % cannot determine number of time points
      break;
    end
  end
end

%% determine whether the specific field makes sense

switch field
  case 'sampleinfo'
    if isfield(data, 'sampleinfo')
      % check whether the existing field makes sense
      if size(data.sampleinfo,1)~=nrpt
        status = false;
      elseif ~all(data.sampleinfo(:,2)-data.sampleinfo(:,1)+1 == nsmp)
        status = false;
      else
        status = true;
      end
    else
      % this is the default
      status = false;
    end
    
  case 'trialinfo'
    if isfield(data, 'trialinfo')
      % check whether the existing field makes sense
      if size(data.trialinfo,1)~=nrpt
        status = false;
      else
        status = true;
      end
    else
      % this is the default
      status = false;
    end
    
  case 'cumtapcnt'
    if isfield(data, 'cumtapcnt')
      % check whether the existing field makes sense
      if size(data.cumtapcnt,1)~=nrpt
        status = false;
      else
        status = true;
      end
    else
      % this is the default
      status = false;
    end
    
  otherwise
    % the default is to assume that it does not make sense
    status = false;
end

if isfield(data, field) && ~status
  % this should be dealt with elsewhere
  ft_warning('inconsistent %s', field);
end
