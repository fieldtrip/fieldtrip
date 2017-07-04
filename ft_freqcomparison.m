function [freq] = ft_freqcomparison(cfg, varargin)

% FT_FREQCOMPARISON is deprecated, please use FT_MATH instead.
%
% See also FT_MATH, FT_FREQBASELINE

% Copyright (C) 2010-2011, Arjen Stolk, DCCN, Donders Institute
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

% DEPRECATED by roboos on 30 July 2013
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2222 for more details
% support for this functionality can be removed at the end of 2013
warning('FT_FREQCOMPARISON is deprecated, please use FT_MATH instead.')

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% nargin check
if length(varargin)~=2
  ft_error('two conditions required as input');
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

if strcmp(varargin{1}.dimord, 'rpt_chan_freq') || strcmp(varargin{1}.dimord, 'subj_chan_freq')
  % frequency comparison for multiple trials/subjects

  if size(varargin{1}.powspctrm,3) ~= size(varargin{2}.powspctrm,3)
    ft_error('input conditions have different sizes');
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
    ft_error('unsupported comparisontype');
  end

elseif strcmp(varargin{1}.dimord, 'chan_freq') || strcmp(varargin{1}.dimord, 'chan_freq_time')
  % frequency comparison for averages

  if size(varargin{1}.powspctrm,2) ~= size(varargin{2}.powspctrm,2)
    ft_error('input conditions have different sizes');
  end

  if strcmp(cfg.comparisontype, 'absolute')
    freq.powspctrm = varargin{2}.powspctrm - varargin{1}.powspctrm;
  elseif strcmp(cfg.comparisontype, 'relchange')
    freq.powspctrm = (varargin{2}.powspctrm - varargin{1}.powspctrm) ./ varargin{1}.powspctrm;
  elseif strcmp(cfg.comparisontype, 'relative')
    freq.powspctrm = varargin{2}.powspctrm ./ varargin{1}.powspctrm;
  else
    ft_error('unsupported comparisontype');
  end

else
  ft_error('unsupported dimord')
end

% user update
fprintf('performing %s comparison \n', cfg.comparisontype);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance freq
ft_postamble history    freq
