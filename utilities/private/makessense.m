function status = makessense(data, field)

% MAKESSENSE determines whether an existing field in the data structure makes sense
%
% Use as
%   status = makessense(data, field)
%
% See also GETDIMSIZ, GETDATFIELD

% Copyright (C) 2018, Robert Oostenveld
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

% it should not only check trial and avg, but also other possible data fields
datfield = setdiff(fieldnames(data), ignorefields('makessense'));

for i=1:numel(datfield)
  if isfield(data, datfield{i})
    dimord = getdimord(data, datfield{i});
    dimtok = tokenize(dimord, '_');
    dimsiz = getdimsiz(data, datfield{i});
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
    
    if isfield(data, 'sampleinfo') && ~status
      % this should be dealt with elsewhere
      ft_warning('inconsistent sampleinfo');
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
    
    if isfield(data, 'trialinfo') && ~status
      % this should be dealt with elsewhere
      ft_warning('inconsistent trialinfo');
    end
    
  otherwise
    % the default is to assume that it does not make sense
    status = false;
end
