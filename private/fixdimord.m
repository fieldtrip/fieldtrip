function [data] = fixdimord(data)

% FIXDIMORD ensures consistency between the dimord string and the axes
% that describe the data dimensions. The main purpose of this function
% is to ensure backward compatibility of all functions with data that has
% been processed by older FieldTrip versions.
%
% Use as
%   [data] = fixdimord(data)
% This will modify the data.dimord field to ensure consistency.
% The name of the axis is the same as the name of the dimord, i.e. if
% dimord='freq_time', then data.freq and data.time should be present.
%
% The default dimensions in the data are described by
%  'time'
%  'freq'
%  'chan'
%  'chancmb'
%  'refchan'
%  'subj'
%  'rpt'
%  'rpttap'
%  'pos'
%  'ori'
%  'rgb'
%  'comp'
%  'voxel'

% Copyright (C) 2009-2014, Robert Oostenveld, Jan-Mathijs Schoffelen
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

% if nargin<2, keepsourcedimord = 0; end
%
% if any(ft_datatype(data, {'source', 'volume'})) && isfield(data, 'dimord') && ~keepsourcedimord
%   % the old source data representation does not have a dimord, whereas the new source data representation does have a dimord
%   ft_warning(sprintf('removing dimord "%s" from source representation data', data.dimord));
%   data = rmfield(data, 'dimord');
%   return
% else
%   % it is ok
%   return
% end

if ~isfield(data, 'dimord')
  if ft_datatype(data, 'raw')
    % it is raw data, which does not have a dimord -> this is ok
    return
  elseif ft_datatype(data, 'comp')
    % it is component data, which resembles raw data -> this is ok
    return
  elseif ft_datatype(data, 'volume')
    % it is volume data, which does not have a dimord -> this is ok
    return
  elseif ft_datatype(data, 'source') || ft_datatype(data, 'parcellation')
    % it is old-style source data -> this is ok
    return
  else
    % find the XXXdimord fields
    fn = fieldnames(data);
    sel = true(size(fn));
    for i=1:length(fn)
      sel(i) = ~isempty(strfind(fn{i}, 'dimord'));
    end
    df = fn(sel);
    % use this function recursively on the XXXdimord fields
    for i=1:length(df)
      data.dimord = data.(df{i});
      data = fixdimord(data);
      data.(df{i}) = data.dimord;
      data = rmfield(data, 'dimord');
    end
    % after the recursive call it should be ok
    return
  end
end % if no dimord

if strcmp(data.dimord, 'voxel')
  % this means that it is position
  data.dimord = 'pos';
end

dimtok = tokenize(data.dimord, '_');
if strncmp('{pos_pos}', data.dimord, 9)
  % keep these together for bivariate source structures
  dimtok = {'{pos_pos}', dimtok{3:end}};
end

for i=1:length(dimtok)
  switch dimtok{i}
    case {'tim' 'time' 'toi' 'latency'}
      dimtok{i} = 'time';
      
    case {'frq' 'freq' 'foi' 'frequency'}
      dimtok{i} = 'freq';
      
    case {'sgn' 'label' 'chan'}
      dimtok{i} = 'chan';
      
    case {'rpt' 'trial'}
      dimtok{i} = 'rpt';
      
    case {'subj' 'subject'}
      dimtok{i} = 'subj';
      
    case {'comp'}
      % don't change, it is ok
      
    case {'sgncmb' 'labelcmb' 'chancmb'}
      dimtok{i} = 'chancmb';
      
    case {'rpttap'}
      % this is a 2D field, coding trials and tapers along the same dimension
      % don't change, it is ok
      
    case {'refchan'}
      % don't change, it is ok
      
    case {'ori'}
      % don't change, it is ok
      
    case {'rgb'}
      % don't change, it is ok
      
    case {'voxel' 'vox' 'repl' 'wcond'}
      % these are used in some FieldTrip functions, but are not considered standard
      ft_warning('unexpected dimord "%s"', data.dimord);
      
    case {'pos'}
      % this is for source data on a 3D grid, a cortical sheet, or unstructured positions
      
    case {'{pos}'}
      % this is for source data on a 3D grid, a cortical sheet, or unstructured positions
      % the data itself is represented in a cell-array, e.g. source.mom or source.leadfield
      
    case {'{pos_pos}'}
      % this is for bivariate source data on a 3D grid, a cortical sheet, or unstructured positions
      
    otherwise
      ft_error('unexpected dimord "%s"', data.dimord);
      
  end % switch dimtok
end % for length dimtok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, 'tim'),         data.time      = data.tim         ; data = rmfield(data, 'tim')        ; end
if isfield(data, 'toi'),         data.time      = data.toi         ; data = rmfield(data, 'toi')        ; end
if isfield(data, 'latency'),     data.time      = data.latency     ; data = rmfield(data, 'latency')    ; end
if isfield(data, 'frq'),         data.freq      = data.frq         ; data = rmfield(data, 'frq')        ; end
if isfield(data, 'foi'),         data.freq      = data.foi         ; data = rmfield(data, 'foi')        ; end
if isfield(data, 'frequency'),   data.freq      = data.frequency   ; data = rmfield(data, 'frequency')  ; end
if isfield(data, 'sgn'),         data.label     = data.sgn         ; data = rmfield(data, 'sgn')        ; end
if isfield(data, 'chan'),        data.label     = data.chan        ; data = rmfield(data, 'chan')       ; end
% if isfield(data, 'trial'),         data.rpt     = data.trial         ; data = rmfield(data, 'trial')        ; end  % DO NOT CONVERT -> this is an exception
if isfield(data, 'subject'),     data.subj      = data.subject     ; data = rmfield(data, 'subject')    ; end
if isfield(data, 'sgncmb'),      data.labelcmb  = data.sgncmb      ; data = rmfield(data, 'sgncmb')     ; end
if isfield(data, 'chancmb'),     data.labelcmb  = data.chancmb     ; data = rmfield(data, 'chancmb')    ; end

% ensure that it is a column
if isfield(data, 'label')
  data.label = data.label(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if isfield(data, 'trial')
%   mat = data.trial;
% elseif isfield(data, 'individual')
%   mat = data.individual;
% elseif isfield(data, 'avg')
%   mat = data.avg;
% elseif isfield(data, 'crsspctrm')
%   mat = data.crsspctrm;
% elseif isfield(data, 'powspctrm')
%   mat = data.powspctrm;
% elseif isfield(data, 'fourierspctrm')
%   mat = data.fourierspctrm;
% end
%
% add the descriptive axis for each dimension
% for i=1:length(dimtok)
%   if isfield(data, dimtok{i})
%     % the dimension is already described with its own axis
%     % data = setfield(data, dimtok{i}, getfield(data, dimtok{i}));
%   else
%     % add an axis to the output data
%     data = setfield(data, dimtok{i}, 1:size(mat,i));
%   end
% end

% undo the tokenization
data.dimord = dimtok{1};
for i=2:length(dimtok)
  data.dimord = [data.dimord '_' dimtok{i}];
end
