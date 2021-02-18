function cfg = ft_badchannel(cfg, data)

% FT_BADCHANNEL tries to identify bad channels in a MEG or EEG dataset. Different
% methods are implemented to identify bad channels, these are described in detail
% below.
%
% STD - identifies channels with a large standard deviation.
%
% RANGE - identifies channels with a large range, i.e. difference between the minimal
% and maximal value observed in that channel.
%
% NEIGHCORREL - identifies channels that have low correlation with their neighbours.
% The rationale is that "bad" channel have inherent noise that is uncorrelated with
% other sensors. The average correlation of MEG magnetometer channels in a BTI248 MEG
% system with their neighbours within a distance of 4 cm is between 0.6 and 0.8.
%
% NEIGHSTDRATIO - identifies channels that have a standard deviation which is very
% different from their neighbours. This computes the difference in the standard
% deviation of each channel to each of its neighbours, relative to that of the
% neighbours. A value of 0.5 has empirically been determined as adequate for
% identifying bad channels in a BTI248 MEG system.
%
% Use as
%   cfg = ft_badchannel(cfg, data)
% where the input data corresponds to the output from FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.method         = 'std', 'range', 'neighcorrel', 'neighstdratio'
%   cfg.threshold      = scalar
%   cfg.neighbours     = neighbourhood structure, see FT_PREPARE_NEIGHBOURS for details
%
% The 'neighcorrel' and 'neighstdratio' methods implement the bad channel detection
% according to the paper "Adding dynamics to the Human Connectome Project with MEG",
% Larson-Prior et al. https://doi.org/10.1016/j.neuroimage.2013.05.056
%
% See also FT_CHANNELREPAIR

% Copyright (C) 2021, Robert Oostenveld
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
ft_preamble loadvar data
ft_preamble provenance
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% store the original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% select channels and trials of interest
tmpcfg = keepfields(cfg, {'channel', 'trials', 'tolerance', 'showcallinfo'});
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

ft_info('concatenating data');
tim = cat(2, data.time{:});
dat = cat(2, data.trial{:});
lab = data.label;
nchan = length(lab);

switch cfg.method
  case 'range'
    % compute the range for each channel
    chanrange = nan(size(lab));
    for i=1:nchan
      chanrange(i) = range(dat(i,:));
    end
    
    % give some feedback
    ft_info('min(range)    = %g\n', min(chanrange));
    ft_info('max(range)    = %g\n', max(chanrange));
    ft_info('mean(range)   = %g\n', mean(chanrange));
    ft_info('median(range) = %g\n', median(chanrange));
    ft_info('mad(range)    = %g\n', mad(chanrange));
    
    % find channels above the threshold
    bad = chanrange>cfg.threshold;
    
  case 'std'
    % compute the standard deviation for each channel
    chanstd = nan(size(lab));
    for i=1:nchan
      chanstd(i) = std(dat(i,:));
    end
    
    % give some feedback
    ft_info('min(range)    = %g\n', min(chanstd));
    ft_info('max(range)    = %g\n', max(chanstd));
    ft_info('mean(range)   = %g\n', mean(chanstd));
    ft_info('median(range) = %g\n', median(chanstd));
    ft_info('mad(range)    = %g\n', mad(chanstd));
    
    % find channels above the threshold
    bad = chanstd>cfg.threshold;
    
  case 'neighstdratio'
    % compute the standard deviation for each channel
    chanstd = nan(size(lab));
    for i=1:nchan
      chanstd(i) = std(dat(i,:));
    end
    
    chanstdratio = nan(length(lab), length(lab));
    for i=1:nchan
      selthis = find(strcmp({cfg.neighbours.label}, lab{i}));
      selnb   = find(ismember(lab, cfg.neighbours(selthis).neighblabel));
      
      labnb = sprintf('%s, ', lab{selnb});
      labnb = labnb(1:end-2);
      ft_debug('channel %s has neighbours %s\n', lab{i}, labnb);
      
      for j=selnb(:)'
        chanstdratio(i,j) = abs(chanstd(i)-chanstd(j)) / chanstd(j);
      end
    end
    
    % give some feedback
    tmp = chanstdratio(~isnan(chanstdratio(:)));
    ft_info('min(range)    = %g\n', min(tmp));
    ft_info('max(range)    = %g\n', max(tmp));
    ft_info('mean(range)   = %g\n', mean(tmp));
    ft_info('median(range) = %g\n', median(tmp));
    ft_info('mad(range)    = %g\n', mad(tmp));
    figure; histogram(tmp, round(length(tmp)/10));
    figure; imagesc(chanstdratio);
    
    % find channels with a median neighbour ABOVE the threshold
    bad = nanmedian(chanstdratio,2)>cfg.threshold;
    
  case 'neighcorrel'
    channeighcorrel = nan(length(lab), length(lab));
    for i=1:nchan
      selthis = find(strcmp({cfg.neighbours.label}, lab{i}));
      selnb   = find(ismember(lab, cfg.neighbours(selthis).neighblabel));
      
      labnb = sprintf('%s, ', lab{selnb});
      labnb = labnb(1:end-2);
      ft_debug('channel %s has neighbours %s\n', lab{i}, labnb);
      
      for j=selnb(:)'
        channeighcorrel(i,j) = corr(dat(i,:)', dat(j,:)');
      end
    end
    
    % give some feedback
    tmp = channeighcorrel(~isnan(channeighcorrel(:)));
    ft_info('min(range)    = %g\n', min(tmp));
    ft_info('max(range)    = %g\n', max(tmp));
    ft_info('mean(range)   = %g\n', mean(tmp));
    ft_info('median(range) = %g\n', median(tmp));
    ft_info('mad(range)    = %g\n', mad(tmp));
    figure; histogram(tmp, round(length(tmp)/10));
    figure; imagesc(channeighcorrel);
    
    % find channels with a median neighbour value BELOW the threshold
    bad = nanmedian(channeighcorrel,2)<cfg.threshold;
    
  otherwise
    ft_error('unsupported method');
end


ft_info('identified %d bad channels from %d channels\n', sum(bad), length(bad));
cfg.badchannel = lab(bad);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance
ft_postamble hastoolbox
