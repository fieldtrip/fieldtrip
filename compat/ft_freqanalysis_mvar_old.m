function [freq] = ft_freqanalysis_mvar(cfg, data)

% FT_FREQANALYSIS_MVAR performs frequency analysis on
% mvar data.
%
% Use as
%   [freq] = ft_freqanalysis(cfg, data)

% Copyright (C) 2009, Jan-Mathijs Schoffelen
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
% $Id: ft_freqanalysis_mvar.m 1980 2010-10-27 10:45:10Z jansch $

if ~isfield(cfg, 'channel'),    cfg.channel    = 'all';          end
if ~isfield(cfg, 'channelcmb'), cfg.channelcmb = {'all' 'all'};  end
if ~isfield(cfg, 'foi'),        cfg.foi        = 'all';          end
if ~isfield(cfg, 'keeptrials'), cfg.keeptrials = 'no';           end
if ~isfield(cfg, 'jackknife'),  cfg.jackknife  = 'no';           end
if ~isfield(cfg, 'keeptapers'), cfg.keeptapers = 'yes';          end
if ~isfield(cfg, 'feedback'),   cfg.feedback   = 'none';         end

if strcmp(cfg.foi, 'all'),
  cfg.foi = (0:1:data.fsampleorig./2);
end

cfg.channel    = ft_channelselection(cfg.channel,      data.label);
%cfg.channelcmb = channelcombination(cfg.channelcmb, data.label);

%keeprpt  = strcmp(cfg.keeptrials, 'yes');
%keeptap  = strcmp(cfg.keeptapers, 'yes');
%dojack   = strcmp(cfg.jackknife,  'yes');
%dozscore = strcmp(cfg.zscore,     'yes');

%if ~keeptap, error('not keeping tapers is not possible yet'); end
%if dojack && keeprpt, error('you cannot simultaneously keep trials and do jackknifing'); end

nfoi     = length(cfg.foi);
if isfield(data, 'time')
  ntoi = numel(data.time);
else
  ntoi = 1;
end
nlag     = size(data.coeffs,3); %change in due course
chanindx = match_str(data.label, cfg.channel);
nchan    = length(chanindx);
label    = data.label(chanindx);

%---allocate memory
h         = complex(zeros(nchan, nchan,  nfoi, ntoi), zeros(nchan, nchan,  nfoi, ntoi));
a         = complex(zeros(nchan, nchan,  nfoi, ntoi), zeros(nchan, nchan,  nfoi, ntoi));
crsspctrm = complex(zeros(nchan, nchan,  nfoi, ntoi), zeros(nchan, nchan,  nfoi, ntoi));

%FIXME build in repetitions

%---loop over the tois
ft_progress('init', cfg.feedback, 'computing MAR-model based TFR');
for j = 1:ntoi
  ft_progress(j/ntoi, 'processing timewindow %d from %d\n', j, ntoi);
 
  %---compute transfer function
  ar = reshape(data.coeffs(:,:,:,j), [nchan nchan*nlag]);
  [h(:,:,:,j), a(:,:,:,j)] = ar2h(ar, cfg.foi, data.fsampleorig);

  %---compute cross-spectra
  nc = data.noisecov(:,:,j);
  for k = 1:nfoi
    tmph               = h(:,:,k,j);
    crsspctrm(:,:,k,j) = tmph*nc*tmph';
  end
end  
ft_progress('close');

%---create output-structure
freq          = [];
freq.label    = label;
freq.freq     = cfg.foi;
%freq.cumtapcnt= ones(ntrl, 1)*ntap;
if ntoi>1
  freq.time   = data.time;
  freq.dimord = 'chan_chan_freq_time';
else
  freq.dimord = 'chan_chan_freq';
end
freq.transfer  = h;
freq.itransfer = a;
freq.noisecov  = data.noisecov;
freq.crsspctrm = crsspctrm;
freq.dof       = data.dof;

try
  cfg.previous = data.cfg;
end
freq.cfg     = cfg; 

%---SUBFUNCTION to compute transfer-function from ar-parameters
function [h, zar] = ar2h(ar, foi, fsample)

nchan = size(ar,1);
ncmb  = nchan*nchan;
nfoi  = length(foi);

%---z-transform frequency axis
zfoi  = exp(-2.*pi.*i.*(foi./fsample));

%---reorganize the ar-parameters
ar  = reshape(ar, [ncmb size(ar,2)./nchan]);
ar  = fliplr([reshape(eye(nchan), [ncmb 1]) -ar]);

zar = complex(zeros(ncmb, nfoi), zeros(ncmb, nfoi));
for k = 1:ncmb
  zar(k,:) = polyval(ar(k,:),zfoi);
end
zar = reshape(zar, [nchan nchan nfoi]);
h   = zeros(size(zar));
for k = 1:nfoi
  h(:,:,k) = inv(zar(:,:,k)); 
end
h   = sqrt(2).*h; %account for the negative frequencies, normalization necessary for
%comparison with non-parametric (fft based) results in fieldtrip
%FIXME probably the normalization for the zero Hz bin is incorrect
zar = zar./sqrt(2);
