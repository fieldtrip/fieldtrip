function [select] = parameterselection(param, data)

% PARAMETERSELECTION selects the parameters that are present as a volume in the data
% add that have a dimension that is compatible with the specified dimensions of the
% volume, i.e. either as a vector or as a 3D volume.
%
% Use as
%   [select] = parameterselection(param, data)
% where
%   param    cell-array, or single string, can be 'all'
%   data     structure with anatomical or functional data
%   select   returns the selected parameters as a cell-array

% Copyright (C) 2005-2008, Robert oostenveld
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

if ischar(param)
  param = {param};   % it should be a cell-array
elseif isempty(param)
  param = {};        % even being empty, it should be a cell-array
end

sel = find(strcmp(param, 'all'));
if ~isempty(sel)
  % the old default was a list of all known volume parameters
  % the new default is to try all fields present in the data
  allparam = fieldnames(data);
  % fields can be nested in source.avg
  if isfield(data, 'avg') && isstruct(data.avg)
    tmp = fieldnames(data.avg);
    for i=1:length(tmp)
      tmp{i} = ['avg.' tmp{i}];
    end
    allparam = cat(1, allparam, tmp);
  end
  % fields can be nested in source.trial
  if isfield(data, 'trial') && isstruct(data.trial)
    tmp = fieldnames(data.trial);
    for i=1:length(tmp)
      tmp{i} = ['trial.' tmp{i}];
    end
    allparam = cat(1, allparam, tmp);
  end
  param(sel) = [];                          % remove the 'all'
  param      = [param(:)' allparam(:)'];    % add the list of all possible parameters, these will be tested later
else
  % check all specified parameters and give support for some parameters like 'pow' and 'coh'
  % which most often will indicate 'avg.pow' and 'avg.coh'
  for i=1:length(param)
    if ~issubfield(data, param{i}) && issubfield(data, ['avg.' param{i}])
      % replace the parameter xxx by avg.xxx
      param{i} = ['avg.' param{i}];
    end
  end
end

% remove empty fields
param(cellfun('isempty', param)) = [];

% ensure that there are no double entries
param = unique(param);

select = {};
for i=1:length(param)
  if issubfield(data, param{i})
    % the field is present, check whether the dimension is correct
    dim = size(getsubfield(data, param{i}));
    if isfield(data, 'dim') && isequal(dim(:), data.dim(:))
      select{end+1} = param{i};
    elseif isfield(data, 'dim') && prod(dim)==prod(data.dim)
      select{end+1} = param{i};
    elseif isfield(data, 'dim') && numel(dim)==3 && isequal(dim(1:3)', data.dim(:))
      select{end+1} = param{i};
    elseif isfield(data, 'pos') && (prod(dim)==size(data.pos, 1) || dim(1)==size(data.pos,1))
      select{end+1} = param{i};
    elseif isfield(data, 'dimord') && (isfield(data, 'pos') || isfield(data, 'transform')),
      dimtok = tokenize(data.dimord, '_');
      nels   = 1;
      for k=1:numel(dimtok)
        if strcmp(dimtok{k}, 'rpt') || strcmp(dimtok{k}, 'rpttap')
          nels = nels*dim(k);
        elseif strcmp(dimtok{k}, 'pos') && isfield(data, 'pos')
          nels = nels*size(data.pos,1);
        elseif strcmp(dimtok{k}, '{pos}') && isfield(data, 'pos')
          nels = nels*size(data.pos,1);
        elseif isfield(data, dimtok{k})
          nels = nels*numel(getfield(data, dimtok{k}));
        end
      end
      if nels==prod(dim),
        select{end+1} = param{i};
      end
    elseif isfield(data, 'dim') && numel(dim)>3 && isequal(dim(1:3), data.dim(1:3))
      select{end+1} = param{i};
    end
  end
end

