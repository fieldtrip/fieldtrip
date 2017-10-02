function [grandavg] = ft_sourcegrandaverage(cfg, varargin)

% FT_SOURCEGRANDAVERAGE averages source reconstructions over either multiple
% subjects or conditions. It computes the average and variance for all
% known source parameters. The output can be used in FT_SOURCESTATISTICS
% with the method 'parametric'.
%
% Alternatively, it can construct an average for multiple input source
% reconstructions in two conditions after randomly reassigning the
% input data over the two conditions. The output then can be used in
% FT_SOURCESTATISTICS with the method 'randomization' or 'randcluster'.
%
% The input source structures should be spatially alligned to each other
% and should have the same positions for the source grid.
%
% Use as
%  [grandavg] = ft_sourcegrandaverage(cfg, source1, source2, ...)
%
% where the source structures are obtained from FT_SOURCEANALYSIS or
% from FT_VOLUMENORMALISE, and the configuration can contain the
% following fields:
%   cfg.parameter          = string, describing the functional data to be processed, e.g. 'pow', 'nai' or 'coh'
%   cfg.keepindividual     = 'no' or 'yes'
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. For this particular function, the input data
% should be structured as a single cell array.
%
% See also FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES, FT_SOURCESTATISTICS, FT_MATH

% Undocumented local options
%  You can also use FT_SOURCEGRANDAVERAGE to compute averages after
% randomizing the assignment of the functional data over two conditions.
% The resulting output can then be used in a statistical test just like
% the randomized single-subject source reconstruction that results from
% randomization in FT_SOURCEANALYSIS. This involves the following options
%   cfg.randomization      = 'no' or 'yes'
%   cfg.permutation        = 'no' or 'yes'
%   cfg.numrandomization   = number, e.g. 500
%   cfg.numpermutation     = number, e.g. 500 or 'all'
%   cfg.c1                 = list with subjects belonging to condition 1 (or A)
%   cfg.c2                 = list with subjects belonging to condition 2 (or B)

% Copyright (C) 2005, Robert Oostenveld
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
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'source'}, 'feedback', 'no', 'inside', 'logical');
  varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'concatenate', 'randomization', 'permutation', 'c1', 'c2'});

% set the defaults
cfg.keepindividual  = ft_getopt(cfg, 'keepindividual', 'no');
cfg.parameter       = ft_getopt(cfg, 'parameter', 'pow');

if strncmp(cfg.parameter, 'avg.', 4)
  cfg.parameter = cfg.parameter(5:end); % remove the 'avg.' part
end
for i=1:length(varargin)
  assert(isfield(varargin{i}, cfg.parameter), 'data does not contain parameter "%s"', cfg.parameter);
end

% check that these fields are identical for each input source
checkfields = {'pos' 'dim' 'xgrid' 'ygrid' 'zgrid' 'transform' 'inside' 'outside'};
for k = 1:numel(checkfields)
  tmpstr = checkfields{k};
  if isfield(varargin{1}, tmpstr)
    tmpvar1 = varargin{1}.(tmpstr);
    for i=2:length(varargin)
      tmpvar2 = varargin{i}.(tmpstr);
      if any(size(tmpvar1)~=size(tmpvar2)) || any(tmpvar1(:)~=tmpvar2(:))
        ft_error('the input sources vary in the field %s', tmpstr);
      end
    end
  end
end

% ensure a consistent selection of the data over all inputs
tmpcfg = keepfields(cfg, {'parameter', 'trials', 'latency', 'frequency', 'foilim', 'showcallinfo'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

% start with an empty output structure
grandavg = [];

if iscell(varargin{1}.(cfg.parameter))

  % collect the data
  dat = cellfun(@getfield, varargin, repmat({cfg.parameter}, size(varargin)), 'UniformOutput', false);

  npos = numel(dat{1});
  nrpt = numel(dat);
  dat  = cat(2, dat{:}); % make it {pos_rpt}

  if isfield(varargin{1}, 'inside')
    % it is logically indexed, take the first inside source location
    probe = find(varargin{1}.inside, 1, 'first');
  else
    % just take the first source position
    probe = 1;
  end

  olddim = size(dat{probe,1});
  newdim = [1 olddim];

  if strcmp(cfg.keepindividual, 'yes')
    dat = cellfun(@reshape,  dat, repmat({newdim}, size(dat)), 'UniformOutput', false);
    for i=1:npos
      dat{i,1} = cat(1, dat{i,:}); % concatenate them into the first one
    end
    grandavg.(cfg.parameter) = dat(:,1);

    if ~isequal(size(grandavg.(cfg.parameter)), size(varargin{1}.(cfg.parameter)))
      % this is a bit unexpected, but let's reshape it back into the original size
      grandavg.(cfg.parameter) = reshape(grandavg.(cfg.parameter), size(varargin{1}.(cfg.parameter)));
    end

    if isfield(varargin{1}, [cfg.parameter 'dimord'])
      dimord = varargin{1}.([cfg.parameter 'dimord']);
      dimtok = tokenize(dimord, '_');
      dimtok = {dimtok{1} 'rpt' dimtok{2:end}};
      dimord = sprintf('%s_', dimtok{:});
      dimord = dimord(1:end-1); % remove the trailing '_'
      grandavg.([cfg.parameter 'dimord']) = dimord;

    elseif isfield(varargin{1}, 'dimord')
      dimord = varargin{1}.dimord;
      dimtok = tokenize(dimord, '_');
      dimtok = {dimtok{1} 'rpt' dimtok{2:end}};
      dimord = sprintf('%s_', dimtok{:});
      dimord = dimord(1:end-1); % remove the trailing '_'
      grandavg.dimord = dimord;
    end

  else
    for i=1:npos
      for j=2:nrpt
        dat{i,1} = dat{i,1} + dat{i,j}; % add them all to the first one
      end
      dat{i,1} = dat{i,1}/nrpt;
    end
    grandavg.(cfg.parameter) = dat(:,1);

    if isfield(varargin{1}, [cfg.parameter 'dimord'])
      grandavg.([cfg.parameter 'dimord']) = varargin{1}.([cfg.parameter 'dimord']);
    elseif isfield(varargin{1}, 'dimord')
      grandavg.dimord = varargin{1}.dimord;
    end

  end % if keepindividual
  clear dat

else
  % determine the dimensions, include the new repetition dimension
  olddim = size(varargin{1}.(cfg.parameter));
  newdim = [1 olddim];

  % collect and reshape the data
  dat = cellfun(@getfield, varargin, repmat({cfg.parameter}, size(varargin)), 'UniformOutput', false);

  if strcmp(cfg.keepindividual, 'yes')
    % concatenate the data into a single array
    dat = cellfun(@reshape,  dat, repmat({newdim}, size(dat)), 'UniformOutput', false);
    grandavg.(cfg.parameter) = cat(1, dat{:});
    if isfield(varargin{1}, [cfg.parameter 'dimord'])
      grandavg.([cfg.parameter 'dimord']) = ['rpt_' varargin{1}.([cfg.parameter 'dimord'])];
    elseif isfield(varargin{1}, 'dimord')
      grandavg.dimord = ['rpt_' varargin{1}.dimord];
    end

  else
    % sum the data in a single array
    for i=2:length(dat)
      dat{1} = dat{1} + dat{i};
    end
    grandavg.(cfg.parameter) = dat{1}/length(varargin);
    if isfield(varargin{1}, [cfg.parameter 'dimord'])
      grandavg.([cfg.parameter 'dimord']) = varargin{1}.([cfg.parameter 'dimord']);
    elseif isfield(varargin{1}, 'dimord')
      grandavg.dimord = varargin{1}.dimord;
    end
  end % if keepindividual
  clear dat

end % if iscell

% the fields that describe the actual data need to be copied over from the input to the output
grandavg = copyfields(varargin{1}, grandavg, {'pos', 'time', 'freq', 'dim', 'transform', 'inside', 'outside', 'unit', 'coordsys'});

% these fields might not be needed
dimord = getdimord(grandavg, cfg.parameter);
if isempty(strfind(dimord, 'time')), grandavg = removefields(grandavg, 'time'); end
if isempty(strfind(dimord, 'freq')), grandavg = removefields(grandavg, 'freq'); end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance grandavg
ft_postamble history    grandavg
ft_postamble savevar    grandavg
