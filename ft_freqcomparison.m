function [freq] = ft_freqcomparison(cfg, varargin)

% FT_FREQCOMPARISON performs a comparison between two (rpt-/subj-)
% chan-freq datasets.
%
% Differently from ft_freqbaseline, ft_freqcomparison requires two datasets
% for input arguments (n==2) where the first dataset is considered the norm.
% The output contains the difference of the second dataset to this norm,
% expressed in units as determined by cfg.comparisontype.
%
% Use as
%   [freq] = ft_freqcomparison(cfg, dataset1, dataset2);
%
% where the freq data comes from FT_SPECEST_MTMFFT and the configuration
% should contain
%
%   cfg.comparisontype = 'absolute' 'relchange' 'relative' (default =
%   'absolute')
%
% See also FT_FREQBASELINE, FT_FREQANALYSIS, FT_TIMELOCKBASELINE

% FIXME add support for cohspctrm

% Copyright (C) 2010-2011, Arjen Stolk, DCCN, Donders Institute
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
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig

% nargin check
if length(varargin)~=2
  error('two conditions required as input');
end

% check if the input data is valid for this function
varargin{1} = ft_checkdata(varargin{1}, 'datatype', 'freq', 'feedback', 'yes');
varargin{2} = ft_checkdata(varargin{2}, 'datatype', 'freq', 'feedback', 'yes');

% set the defaults
if ~isfield(cfg, 'comparisontype'), cfg.comparisontype = 'absolute';     end

% we have to ensure that we don't end up with an inconsistent dataset
% remove cross-spectral densities since coherence cannot be computed any more
if isfield(varargin{1}, 'crsspctrm')
  varargin{1} = rmfield(varargin{1}, 'crsspctrm');
end
if isfield(varargin{1}, 'cohspctrmsem')
  varargin{1} = rmfield(varargin{1}, 'cohspctrmsem');
end
if isfield(varargin{1}, 'powspctrmsem')
  varargin{1} = rmfield(varargin{1}, 'powspctrmsem');
end
if isfield(varargin{2}, 'crsspctrm')
  varargin{2} = rmfield(varargin{2}, 'crsspctrm');
end
if isfield(varargin{2}, 'cohspctrmsem')
  varargin{2} = rmfield(varargin{2}, 'cohspctrmsem');
end
if isfield(varargin{2}, 'powspctrmsem')
  varargin{2} = rmfield(varargin{2}, 'powspctrmsem');
end

% initiate output variable
freq = varargin{1};

% frequency comparison for multiple trials/subjects
if strcmp(varargin{1}.dimord, 'rpt_chan_freq') || strcmp(varargin{1}.dimord, 'subj_chan_freq')
  
  if size(varargin{1}.powspctrm,3) ~= size(varargin{2}.powspctrm,3)
    error('input conditions have different sizes');
  end
  
  if strcmp(cfg.comparisontype, 'absolute')
    for j = 1:size(varargin{2}.powspctrm,1)
      freq.powspctrm(j,:,:) = varargin{2}.powspctrm(j,:,:) - mean(varargin{1}.powspctrm,1);
    end
  elseif strcmp(cfg.comparisontype, 'relchange')
    for j = 1:size(varargin{2}.powspctrm,1)
      freq.powspctrm(j,:,:) = (varargin{2}.powspctrm(j,:,:) - mean(varargin{1}.powspctrm,1)) ./ mean(varargin{1}.powspctrm,1);
    end
  elseif strcmp(cfg.comparisontype, 'relative')
    for j = 1:size(varargin{2}.powspctrm,1)
      freq.powspctrm(j,:,:) = varargin{2}.powspctrm(j,:,:) ./ mean(varargin{1}.powspctrm,1);
    end
  else
    error('unsupported comparisontype');
  end
  
  % frequency comparison for averages
elseif strcmp(varargin{1}.dimord, 'chan_freq')
  
  if size(varargin{1}.powspctrm,2) ~= size(varargin{2}.powspctrm,2)
    error('input conditions have different sizes');
  end
  
  if strcmp(cfg.comparisontype, 'absolute')
    freq.powspctrm = varargin{2}.powspctrm - varargin{1}.powspctrm;
  elseif strcmp(cfg.comparisontype, 'relchange')
    freq.powspctrm = (varargin{2}.powspctrm - varargin{1}.powspctrm) ./ varargin{1}.powspctrm;
  elseif strcmp(cfg.comparisontype, 'relative')
    freq.powspctrm = varargin{2}.powspctrm ./ varargin{1}.powspctrm;
  else
    error('unsupported comparisontype');
  end
end

% user update
fprintf('performing %s comparison \n', cfg.comparisontype);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance

