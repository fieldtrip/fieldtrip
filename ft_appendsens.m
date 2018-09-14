function [sens] = ft_appendsens(cfg, varargin)

% FT_APPENDSENS concatenates multiple sensor definitions that have been processed
% separately.
%
% Use as
%   combined = ft_appendsens(cfg, sens1, sens2, ...)
%
% A call to FT_APPENDSENS results in the label, pos and ori fields to be
% concatenated, and the tra matrix to be merged. Any duplicates will be removed.
% The labelold and chanposold fields are kept under the condition that they
% are identical across the inputs.
%
% See also FT_ELECTRODEPLACEMENT, FT_ELECTRODEREALIGN, FT_DATAYPE_SENS,
% FT_APPENDDATA, FT_APPENDTIMELOCK, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2017-2018, Arjen Stolk
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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
% and ensure it is up to the latest standards
for i=1:length(varargin)
  assert(ft_datatype(varargin{i}, 'sens'), 'incorrect input, should be a grad/elec/opto structure');
  varargin{i} = ft_datatype_sens(varargin{i});
end

% do a basic check whether the senstype, units, and coordinate systems match
senstype = cell(1,length(varargin));
for i=1:length(varargin)
  senstype{i} = ft_senstype(varargin{i});
end
typematch = all(strcmp(senstype{1}, senstype));

unitmatch = 1;
if isfield(varargin{1}, 'unit')
  unit = cell(1,length(varargin));
  for i=1:length(varargin)
    unit{i} = varargin{i}.unit;
  end
  unitmatch = all(strcmp(unit{1}, unit));
end

coordsysmatch = 1;
if isfield(varargin{1}, 'coordsys')
  coordsys = cell(1,length(varargin));
  for i=1:length(varargin)
    coordsys{i} = varargin{i}.coordsys;
  end
  coordsysmatch = all(strcmp(coordsys{1}, coordsys));
end

if ~typematch || ~unitmatch || ~coordsysmatch
  ft_error('the senstype, units, or coordinate systems of the inputs do not match');
end

% keep these fields (when present) in the output
sens = keepfields(varargin{1}, {'type', 'unit', 'coordsys'});

% make inventory
haselecpos = 0;
hascoilpos = 0;
hascoilori = 0;
haschanori = 0;
hasoptopos = 0;
hastra     = 0;
haslabelold = 0;
haschanposold = 0;
for i=1:length(varargin)
  % the following fields should be present in any sens structure
  if isfield(varargin{i}, 'label')
    label{i} = varargin{i}.label(:); % ensure column orientation
  end
  if isfield(varargin{i}, 'chanpos')
    chanpos{i} = varargin{i}.chanpos;
  end
  
  % some the following fields are likely present in a sens structure
  if isfield(varargin{i}, 'elecpos') % EEG
    elecpos{i} = varargin{i}.elecpos;
    haselecpos = 1;
  end
  if isfield(varargin{i}, 'coilpos') % MEG
    coilpos{i} = varargin{i}.coilpos;
    hascoilpos = 1;
  end
  if isfield(varargin{i}, 'coilori') % MEG
    coilori{i} = varargin{i}.coilori;
    hascoilori = 1;
  end
  if isfield(varargin{i}, 'chanori') % MEG
    chanori{i} = varargin{i}.chanori;
    haschanori = 1;
  end
  if isfield(varargin{i}, 'optopos') % NIRS
    optopos{i} = varargin{i}.optopos;
    hasoptopos = 1;
  end
  if isfield(varargin{i}, 'tra') % tra
    tra{i} = varargin{i}.tra;
    hastra = 1;
  end
  
  % the following fields might be present in a sens structure
  if isfield(varargin{i}, 'labelold')
    labelold{i} = varargin{i}.labelold(:); % ensure column orientation
    haslabelold = 1;
  end
  if isfield(varargin{i}, 'chanposold')
    chanposold{i} = varargin{i}.chanposold;
    haschanposold = 1;
  end
end

% concatenate the main fields and remove duplicates
sens.label = cat(1,label{:});
[~, labidx] = unique(sens.label, 'stable');
if ~isequal(numel(labidx), numel(sens.label))
  fprintf('removing duplicate labels\n')
  sens.label = sens.label(labidx);
end

sens.chanpos = cat(1,chanpos{:});
[~, chanidx] = unique(sens.chanpos, 'rows', 'stable');
if ~isequal(numel(chanidx), size(sens.chanpos,1))
  fprintf('removing duplicate channels\n')
  sens.chanpos = sens.chanpos(chanidx,:);
  if ~isequal(labidx, chanidx) % check for matching order
    ft_error('inconsistent order or number of channel labels and positions')
  end
end
if ~isequal(numel(sens.label), size(sens.chanpos,1)) % check for matching number
  ft_error('inconsistent number of channel labels and positions')
end
if haschanori
  sens.chanori = cat(1,chanori{:});
  if ~isequal(numel(chanidx), size(sens.chanpos,1))
    sens.chanori = sens.chanori(chanidx,:); % chanori should match chanpos
  end
end

if hastra && ~any(cellfun(@isempty, tra))
  sens.tra = [];
  for t = 1:numel(tra)
    trarow = [1:size(tra{t},1)]+size(sens.tra,1);
    tracol = [1:size(tra{t},2)]+size(sens.tra,2);
    sens.tra(trarow, tracol) = tra{t};
  end
end

if haselecpos
  sens.elecpos = cat(1,elecpos{:});
  [~, elecidx, elecidx2] = unique(sens.elecpos, 'rows', 'stable');
  if ~isequal(numel(elecidx), size(sens.elecpos,1))
    fprintf('removing duplicate electrodes\n')
    sens.elecpos = sens.elecpos(elecidx,:);
  end
  if isfield(sens, 'tra')
    % shape duplicates into a single column, if necessary
    for idx = 1:numel(elecidx)
      tmp(:, idx) = sum(sens.tra(chanidx, find(idx==elecidx2)),2); % find and take sum over duplicates
    end
    sens.tra = tmp;
    % check for expected size and non-empty rows
    if ~isequal(size(sens.tra,1), size(sens.chanpos,1)) || ~isequal(size(sens.tra,2), size(sens.elecpos,1)) ...
        || any(any(sens.tra,2)==0) % all channels need to be linked to their origins
      fprintf('removing inconsistent tra matrix\n')
      sens = rmfield(sens, 'tra');
    end
  end
end

if hasoptopos
  sens.optopos = cat(1,optopos{:});
  [~, optoidx, optoidx2] = unique(sens.optopos, 'rows', 'stable');
  if ~isequal(numel(optoidx), size(sens.optopos,1))
    fprintf('removing duplicate optodes\n')
    sens.optopos = sens.optopos(optoidx,:);
  end
  if isfield(sens, 'tra')
    % shape duplicates into a single column, if necessary
    for idx = 1:numel(optoidx)
      tmp(:, idx) = sum(sens.tra(chanidx, find(idx==optoidx2)),2); % find and take sum over duplicates
    end
    sens.tra = tmp;
    % check for expected size and non-empty rows
    if ~isequal(size(sens.tra,1), size(sens.chanpos,1)) || ~isequal(size(sens.tra,2), size(sens.optopos,1)) ...
        || any(any(sens.tra,2)==0) % all channels need to be linked to their origins
      fprintf('removing inconsistent tra matrix\n')
      sens = rmfield(sens, 'tra');
    end
  end
end

if hascoilpos
  sens.coilpos = cat(1,coilpos{:});
  [~, coilidx, coilidx2] = unique(sens.coilpos, 'rows', 'stable');
  if ~isequal(numel(coilidx), size(sens.coilpos,1))
    fprintf('removing duplicate coils\n')
    sens.coilpos = sens.coilpos(coilidx,:);
  end
  if isfield(sens, 'tra')
    % shape duplicates into a single column, if necessary
    for idx = 1:numel(coilidx)
      tmp(:, idx) = sum(sens.tra(chanidx, find(idx==coilidx2)),2); % find and take sum over duplicates
    end
    sens.tra = tmp;
    % check for expected size and non-empty rows
    if ~isequal(size(sens.tra,1), size(sens.chanpos,1)) || ~isequal(size(sens.tra,2), size(sens.coilpos,1)) ...
        || any(any(sens.tra,2)==0) % all channels need to be linked to their origins
      fprintf('removing inconsistent tra matrix\n')
      sens = rmfield(sens, 'tra');
    end
  end
end
if hascoilori
  sens.coilori = cat(1,coilori{:});
  if ~isequal(numel(coilidx), size(sens.coilori,1))
    sens.coilori = sens.coilori(coilidx,:); % coilori should match coilpos
  end
end

% keep the following fields only when identical across inputs
if haslabelold && all(isequal(labelold{1}, labelold{:})) % labeloldmatch
  sens.labelold = labelold{1};
end
if haschanposold && all(isequal(chanposold{1}, chanposold{:})) % chanposoldmatch
  sens.chanposold = chanposold{1};
end

% ensure up-to-date and consistent output sensor description
sens = ft_datatype_sens(sens);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance sens
ft_postamble history sens
ft_postamble savevar sens
