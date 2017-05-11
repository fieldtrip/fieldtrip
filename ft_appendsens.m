function [sens] = ft_appendsens(cfg, varargin)

% FT_APPENDSENS concatenates multiple sensor definitions that have been processed
% separately. This is specifically designed for multiple intracranial electrode
% descriptions, e.g. to combine ECoG and sEEG.
%
% Use as
%   combined = ft_appendsens(cfg, sens1, sens2, ...)
%
% A call to FT_APPENDSENS results in the label, chanpos, and when applicable, 
% the chanori fields to be concatenated. Any duplicates will be removed. 
% Moreover, the elecpos, labelold, and chanposold fields are copied over
% from the first input under the condition that these fields are identical 
% across the inputs. If applicable, the tra matrices are also merged.
%
% See also FT_ELECTRODEPLACEMENT, FT_ELECTRODEREALIGN, FT_DATAYPE_SENS,
% FT_APPENDDATA, FT_APPENDTIMELOCK, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2017, Arjen Stolk
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
  error('the senstype, units, or coordinate systems of the inputs are not equal');
end

% keep these fields (when present) in the output
sens = keepfields(varargin{1}, {'type', 'unit', 'coordsys'});

% make inventory
haslabelold = 0;
haschanposold = 0;
haselecpos = 0;
hascoilpos = 0;
hascoilori = 0;
haschanori = 0;
hasoptopos = 0;
for i=1:length(varargin)
  % the following fields should be present in any sens structure
  if isfield(varargin{i}, 'label')
    label{i} = varargin{i}.label;
  end
  if isfield(varargin{i}, 'chanpos')
    chanpos{i} = varargin{i}.chanpos;
  end
  
  % the following fields may be present in a subset of sens structures
  if isfield(varargin{i}, 'labelold')
    labelold{i} = varargin{i}.labelold;
    haslabelold = 1;
  end
  if isfield(varargin{i}, 'chanposold')
    chanposold{i} = varargin{i}.chanposold;
    haschanposold = 1;
  end
  
  % the following fields might be present in a sens structure
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
end

% concatenate channels - see test_pull393.m for a test script
sens.label = cat(1,label{:});
sens.chanpos = cat(1,chanpos{:});
if haschanori
  sens.chanori = cat(1,chanori{:});
end

% remove duplicate channels
[~, idx] = unique(sens.label);
if ~isequal(numel(idx), numel(sens.label))
  warning('duplicate channel labels found and removed, assuming that the chanpos fields match')
  idx = sort(idx);
  sens.label = sens.label(idx);
end


% copy sensor information when identical across inputs (thus likely to have the same origin)
if haselecpos && all(isequal(elecpos{1}, elecpos{:})) % elecposmatch
  sens.elecpos = elecpos{1};
  if hastra && ~any(cellfun(@isempty, tra))
    sens.tra = cat(1,tra{:});
  end
end
if hasoptopos && all(isequal(optopos{1}, optopos{:})) % optoposmatch
  sens.optopos = optopos{1};
  if hastra && ~any(cellfun(@isempty, tra))
    sens.tra = cat(1,tra{:});
  end
end
if hascoilpos && all(isequal(coilpos{1}, coilpos{:})) % coilposmatch
  sens.coilpos = coilpos{1};
  if hastra && ~any(cellfun(@isempty, tra))
    sens.tra = cat(1,tra{:});
  end
end
if hascoilori && all(isequal(coilori{1}, coilori{:})) % coilorimatch
  sens.coilori = coilori{1};
end
if haslabelold && all(strcmp(labelold{1}, labelold)) % labeloldmatch
  sens.labelold = labelold{1};
end
if haschanposold && all(isequal(chanposold{1}, chanposold{:})) % chanposoldmatch
  sens.chanposold = chanposold{1};
end

% ensure up-to-date output sensor description
sens = ft_datatype_sens(sens);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance sens
ft_postamble history sens
ft_postamble savevar sens
