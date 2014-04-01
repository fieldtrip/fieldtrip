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
% See also FT_SOURCEANALYSIS, FT_VOLUMENORMALISE, FT_SOURCESTATISTICS

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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar varargin

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'source'}, 'feedback', 'no', 'sourcerepresentation', 'new');
  varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'concatenate', 'randomization', 'permutation', 'c1', 'c2'});

% set the defaults
if ~isfield(cfg, 'parameter'),      cfg.parameter = 'pow';     end
if ~isfield(cfg, 'keepindividual'), cfg.keepindividual = 'no'; end

% check that these fields are identical for each input source
checkfields = {'pos' 'dim' 'xgrid' 'ygrid' 'zgrid' 'transform' 'inside' 'outside'};
for k = 1:numel(checkfields)
  tmpstr = checkfields{k};
  if isfield(varargin{1}, tmpstr)
    tmpvar1 = varargin{1}.(tmpstr);
    for i=2:length(varargin)
      tmpvar2 = varargin{i}.(tmpstr);
      if any(size(tmpvar1)~=size(tmpvar2)) || any(tmpvar1(:)~=tmpvar2(:))
        error('the input sources vary in the field %s', tmpstr);
      end
    end
    grandavg.(tmpstr) = varargin{1}.(tmpstr);
  end
end

% ensure a consistent selection of the data over all inputs
[varargin{:}] = ft_selectdata(cfg, varargin{:});

% deal with the fields that are always present
grandavg.pos = varargin{1}.pos;
if isfield(varargin{1}, 'inside')
  grandavg.inside = varargin{1}.inside;
end
if isfield(varargin{1}, 'outside')
  grandavg.outside = varargin{1}.outside;
end
if isfield(varargin{1}, 'dim')
  grandavg.dim = varargin{1}.dim;
end

if strcmp(cfg.keepindividual, 'yes')
  olddim = size(varargin{i}.(cfg.parameter));
  newdim = [1 olddim];
  dat = cell(size(varargin));
  for i=1:length(varargin)
    dat{i} = reshape(varargin{i}.(cfg.parameter), newdim);
  end
  grandavg.(cfg.parameter) = cat(1, dat{:});
  if isfield(varargin{1}, [cfg.parameter 'dimord'])
    grandavg.([cfg.parameter 'dimord']) = ['rpt_' varargin{1}.([cfg.parameter 'dimord'])];
  elseif isfield(varargin{1}, 'dimord')
    grandavg.dimord = ['rpt_' varargin{1}.dimord];
  end
  
else
  dat = varargin{1}.(cfg.parameter);
  for i=2:length(varargin)
    dat = dat + varargin{i}.(cfg.parameter);
  end
  grandavg.(cfg.parameter) = dat/length(varargin);
  if isfield(varargin{1}, [cfg.parameter 'dimord'])
    grandavg.([cfg.parameter 'dimord']) = varargin{1}.([cfg.parameter 'dimord']);
  elseif isfield(varargin{1}, 'dimord')
    grandavg.dimord = varargin{1}.('dimord');
  end
end % if keepindividual

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
ft_postamble history grandavg
ft_postamble savevar grandavg
