function [freq] = ft_freqcomparison(cfg, varargin)

% FT_FREQCOMPARISON performs a baseline-like comparison for channel-frequency data
%
% Use as
%   [freq] = ft_freqcomparison(cfg, dataset1, dataset2);
%
% where the freq data comes from FREQANALYSIS_MTMFFT and the configuration
% should contain
%   cfg.baselinetype = 'absolute' 'relchange' 'relative' (default =
%   'absolute')
%
% where the first dataset is seen as baseline
%
% See also FT_FREQBASELINE, FT_FREQANALYSIS, FT_TIMELOCKBASELINE
%
% Copyright (C) 2010-2011, Arjen Stolk, DCCN, Donders Institute
% A.Stolk@donders.ru.nl
% 
% not yet supported: cohspctrm

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

ft_defaults

% nargin check
if nargin ~= 3
    error('two conditions required for input');
end

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
varargin{1} = ft_checkdata(varargin{1}, 'datatype', 'freq', 'feedback', 'yes');
varargin{2} = ft_checkdata(varargin{2}, 'datatype', 'freq', 'feedback', 'yes');

% set the defaults
if ~isfield(cfg, 'baselinetype'), cfg.baselinetype = 'absolute';     end

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

% frequency comparison for multiple trials/subjects
if strcmp(varargin{1}.dimord, 'rpt_chan_freq') || strcmp(varargin{1}.dimord, 'subj_chan_freq')
    
    if size(varargin{1}.powspctrm,3) ~= size(varargin{2}.powspctrm,3)
        error('input conditions have different sizes');
    end
    
    if strcmp(cfg.baselinetype, 'absolute')
        for j = 1:size(varargin{2}.powspctrm,1)
            freq.powspctrm(j,:,:) = varargin{2}.powspctrm(j,:,:) - mean(varargin{1}.powspctrm,1);
        end
    elseif strcmp(cfg.baselinetype, 'relchange')
        for j = 1:size(varargin{2}.powspctrm,1)
            freq.powspctrm(j,:,:) = (varargin{2}.powspctrm(j,:,:) - mean(varargin{1}.powspctrm,1)) ./ mean(varargin{1}.powspctrm,1);
        end
    elseif strcmp(cfg.baselinetype, 'relative')
        for j = 1:size(varargin{2}.powspctrm,1)
            freq.powspctrm(j,:,:) = varargin{2}.powspctrm(j,:,:) ./ mean(varargin{1}.powspctrm,1);
        end
    else
        error('unsupported baselinetype');
    end
  
% frequency comparison for averages 
elseif strcmp(varargin{1}.dimord, 'chan_freq')
    
    if size(varargin{1}.powspctrm,2) ~= size(varargin{2}.powspctrm,2)
        error('input conditions have different sizes');
    end
    
    if strcmp(cfg.baselinetype, 'absolute')
        freq.powspctrm = varargin{2}.powspctrm - varargin{1}.powspctrm;
    elseif strcmp(cfg.baselinetype, 'relchange')
        freq.powspctrm = (varargin{2}.powspctrm - varargin{1}.powspctrm) ./ varargin{1}.powspctrm;
    elseif strcmp(cfg.baselinetype, 'relative')
        freq.powspctrm = varargin{2}.powspctrm ./ varargin{1}.powspctrm;
    else
        error('unsupported baselinetype');
    end
end

% user update
fprintf('performing %s comparison \n', cfg.baselinetype);

% add the remaining parameters
if isfield(varargin{1}, 'label')
    freq.label = varargin{1}.label;
end
if isfield(varargin{1}, 'dimord')
    freq.dimord = varargin{1}.dimord;
end
if isfield(varargin{1}, 'freq')
    freq.freq = varargin{1}.freq;
end
if isfield(varargin{1}, 'cumsumcnt')
    freq.cumsumcnt = varargin{1}.cumsumcnt;
end
if isfield(varargin{1}, 'cumtapcnt')
    freq.cumtapcnt = varargin{1}.cumtapcnt;
end
if isfield(varargin{1}, 'grad')
    freq.grad = varargin{1}.grad;
end
if isfield(varargin{1}, 'cfg')
    freq.cfg = varargin{1}.cfg;
end
