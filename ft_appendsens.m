function [sens] = ft_appendsens(cfg, varargin)

% FT_APPENDSENS concatenates multiple sensor definitions that have been processed
% separately.
%
% Use as
%   combined = ft_appendsens(cfg, sens1, sens2, ...)
%
% A call to FT_APPENDSENS results in the label, pos and ori fields to be
% concatenated, and the tra matrix to be merged. Any duplicate electrodes
% will be removed. The labelold and chanposold fields are kept under the
% condition that they are identical across the inputs.
%
% See also FT_ELECTRODEPLACEMENT, FT_ELECTRODEREALIGN, FT_DATAYPE_SENS,
% FT_APPENDDATA, FT_APPENDTIMELOCK, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2017-2018, Arjen Stolk
% Copyright (C) 2023, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    varargin
ft_preamble provenance varargin

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check that the input data is valid for this function
% and ensure it is up to the latest standards
for i=1:length(varargin)
  assert(ft_datatype(varargin{i}, 'sens'), 'incorrect input, should be a grad/elec/opto structure');
  varargin{i} = ft_datatype_sens(varargin{i});
end

% do a basic check whether the senstype, units, and coordinate systems match
for i=2:length(varargin)
  if ~strcmp(ft_senstype(varargin{1}), ft_senstype(varargin{i}))
    ft_error('the senstype does not match');
  end
end

if isfield(varargin{1}, 'unit')
  for i=2:length(varargin)
    if ~strcmp(varargin{1}.unit, varargin{i}.unit)
      ft_error('the unit does not match');
    end
  end
end

if isfield(varargin{1}, 'coordsys')
  for i=2:length(varargin)
    if ~strcmp(varargin{1}.coordsys, varargin{i}.coordsys)
      ft_error('the coordsys does not match');
    end
  end
end

% keep these fields (when present) in the output
sens = keepfields(varargin{1}, {'type', 'unit', 'coordsys'});

% these must be present
label      = {};
chanpos    = {};
% these are optional
chantype   = {};
chanunit   = {};
chanori    = {};
elecpos    = {};
optopos    = {};
coilpos    = {};
coilori    = {};
tra        = {};
labelold   = {};
chanposold = {};

for i=1:length(varargin)
  % the following fields should be present in any sens structure
  if isfield(varargin{i}, 'label')
    label{i} = varargin{i}.label(:); % ensure column orientation
  end
  if isfield(varargin{i}, 'chanpos')
    chanpos{i} = varargin{i}.chanpos;
  end

  % some the following fields are likely present in a sens structure
  if isfield(varargin{i}, 'chantype')
    chantype{i} = varargin{i}.chantype(:); % ensure column orientation
  end
  if isfield(varargin{i}, 'chanunit')
    chanunit{i} = varargin{i}.chanunit(:); % ensure column orientation
  end
  if isfield(varargin{i}, 'chanori')
    chanori{i} = varargin{i}.chanori;
  end
  if isfield(varargin{i}, 'elecpos') % EEG
    elecpos{i} = varargin{i}.elecpos;
  end
  if isfield(varargin{i}, 'optopos') % NIRS
    optopos{i} = varargin{i}.optopos;
  end
  if isfield(varargin{i}, 'coilpos') % MEG
    coilpos{i} = varargin{i}.coilpos;
  end
  if isfield(varargin{i}, 'coilori') % MEG
    coilori{i} = varargin{i}.coilori;
  end
  if isfield(varargin{i}, 'tra') % tra
    tra{i} = varargin{i}.tra;
  end

  % the following fields might be present in a sens structure
  if isfield(varargin{i}, 'labelold')
    labelold{i} = varargin{i}.labelold(:); % ensure column orientation
  end
  if isfield(varargin{i}, 'chanposold')
    chanposold{i} = varargin{i}.chanposold;
  end
end

% concatenate all fields
if ~isempty(label)
  sens.label = cat(1,label{:});
end

if ~isempty(chantype)
  sens.chantype = cat(1,chantype{:});
end

if ~isempty(chanunit)
  sens.chanunit = cat(1,chanunit{:});
end

if ~isempty(chanpos)
  sens.chanpos = cat(1,chanpos{:});
end

if ~isempty(chanori)
  sens.chanori = cat(1,chanori{:});
end

if ~isempty(tra) && ~any(cellfun(@isempty, tra))
  sens.tra = [];
  % this needs to be blockwise shifted
  for t = 1:numel(tra)
    trarow = [1:size(tra{t},1)] + size(sens.tra,1);
    tracol = [1:size(tra{t},2)] + size(sens.tra,2);
    sens.tra(trarow, tracol) = tra{t};
  end
end

if ~isempty(elecpos)
  sens.elecpos = cat(1,elecpos{:});
end

if ~isempty(optopos)
  sens.optopos = cat(1,optopos{:});
end

if ~isempty(coilpos)
  sens.coilpos = cat(1,coilpos{:});
end

if ~isempty(coilori)
  sens.coilori = cat(1,coilori{:});
end

% keep the following fields only when identical across inputs
if ~isempty(labelold) && all(isequal(labelold{1}, labelold{:})) % labeloldmatch
  sens.labelold = labelold{1};
end

if ~isempty(chanposold) && all(isequal(chanposold{1}, chanposold{:})) % chanposoldmatch
  sens.chanposold = chanposold{1};
end

% remove duplicate channels
[lab, labidx] = unique(sens.label, 'stable');
if ~isequal(lab, sens.label)
  ft_warning('removing channels with duplicate labels')

  % when electrodes, coils or optodes are removed, the columns of the tra matrix also needs to be pruned
  prunecolumns = false;

  % these need to be checked before pruning chanpos and chanori
  if isfield(sens, 'elecpos') && isequal(sens.chanpos, sens.elecpos)
    sens.elecpos = sens.elecpos(labidx,:);
    prunecolumns = true;
  end
  if isfield(sens, 'optopos') && isequal(sens.chanpos, sens.optopos)
    sens.optopos = sens.optopos(labidx,:);
    prunecolumns = true;
  end
  if isfield(sens, 'coilpos') && isequal(sens.chanpos, sens.coilpos)
    sens.coilpos = sens.coilpos(labidx,:);
    prunecolumns = true;
  end
  if isfield(sens, 'coilori') && isequal(sens.chanori, sens.coilori)
    sens.coilori = sens.coilori(labidx,:);
    prunecolumns = true;
  end

  if isfield(sens, 'tra')
    if prunecolumns
      sens.tra = sens.tra(labidx,labidx);
    else
      sens.tra = sens.tra(labidx,:);
    end
  end

  sens.label = sens.label(labidx);
  sens.chanpos = sens.chanpos(labidx,:);

  if isfield(sens, 'chantype')
    sens.chantype = sens.chantype(labidx);
  end
  if isfield(sens, 'chanunit')
    sens.chanunit = sens.chanunit(labidx);
  end
  if isfield(sens, 'chanori')
    sens.chanori = sens.chanori(labidx,:);
  end

end % if duplicate channels

% remove duplicate sensors
if isfield(sens, 'elecpos') && ~isequal(sens.chanpos, sens.elecpos)
  ft_warning('removing duplicate electrodes')
  [pos, index, swap] = unique(sens.elecpos, 'rows', 'stable');
  if isfield(sens, 'tra')
    % reconstruct the tra matrix
    nchan = length(sens.label);
    nsens = length(index);
    tra = zeros(nchan, nsens);
    for i=1:nchan
      for j=1:nsens
        tra(i,j) = sum(sens.tra(i,swap==index(j)));
      end
    end
    sens.tra = tra;
  end
  sens.elecpos = pos;

elseif isfield(sens, 'optopos') && ~isequal(sens.chanpos, sens.optopos)
  ft_warning('removing duplicate optodes')
  [pos, index, swap] = unique(sens.optopos, 'rows', 'stable');
  if isfield(sens, 'tra')
    % reconstruct the tra matrix
    nchan = length(sens.label);
    nsens = length(index);
    tra = zeros(nchan, nsens);
    for i=1:nchan
      for j=1:nsens
        tra(i,j) = sum(sens.tra(i,swap==index(j)));
      end
    end
    sens.tra = tra;
  end
  sens.optopos = pos;

elseif isfield(sens, 'coilpos') && ~isequal(sens.chanpos, sens.coilpos)
  ft_warning('removing duplicate magnetometer or gradiometer coils')
  [posori, index, swap] = unique([sens.coilpos sens.coilori], 'rows', 'stable');
  if isfield(sens, 'tra')
    % reconstruct the tra matrix
    nchan = length(sens.label);
    nsens = length(index);
    tra = zeros(nchan, nsens);
    for i=1:nchan
      for j=1:nsens
        tra(i,j) = sum(sens.tra(i,swap==index(j)));
      end
    end
    sens.tra = tra;
  end
  sens.coilpos = posori(:,1:3);
  sens.coilori = posori(:,4:6);
  
end % if duplicate sensors

% ensure up-to-date and consistent output sensor description
sens = ft_datatype_sens(sens);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous varargin
ft_postamble provenance sens
ft_postamble history sens
ft_postamble savevar sens
