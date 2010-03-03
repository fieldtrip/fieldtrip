function [freq] = ft_freqanalysis_mvar(cfg, data)

% FT_FREQANALYSIS_MVAR performs frequency analysis on
% mvar data.
%
% Use as
%   [freq] = ft_freqanalysis(cfg, data)

% Copyright (C) 2009, Jan-Mathijs Schoffelen
% $Log: freqanalysis_mvar.m,v $
% Revision 1.1  2009/10/02 13:55:48  jansch
% first implementation in fieldtrip
%

if ~isfield(cfg, 'channel'),    cfg.channel    = 'all';          end
if ~isfield(cfg, 'channelcmb'), cfg.channelcmb = {'all' 'all'};  end
if ~isfield(cfg, 'foi'),        cfg.foi        = 'all';          end
if ~isfield(cfg, 'keeptrials'), cfg.keeptrials = 'no';           end
if ~isfield(cfg, 'jackknife'),  cfg.jackknife  = 'no';           end
if ~isfield(cfg, 'keeptapers'), cfg.keeptapers = 'yes';          end
if ~isfield(cfg, 'feedback'),   cfg.feedback   = 'none';         end

if strcmp(cfg.foi, 'all'),
  cfg.foi = [0:1:data.fsampleorig./2];
end

cfg.channel    = channelselection(cfg.channel,      data.label);
%cfg.channelcmb = channelcombination(cfg.channelcmb, data.label);

%keeprpt  = strcmp(cfg.keeptrials, 'yes');
%keeptap  = strcmp(cfg.keeptapers, 'yes');
%dojack   = strcmp(cfg.jackknife,  'yes');
%dozscore = strcmp(cfg.zscore,     'yes');

%if ~keeptap, error('not keeping tapers is not possible yet'); end
%if dojack && keeprpt, error('you cannot simultaneously keep trials and do jackknifing'); end

nfoi     = length(cfg.foi);
ntoi     = length(data.time);
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
progress('init', cfg.feedback, 'computing MAR-model based TFR');
for j = 1:ntoi
  progress(j/ntoi, 'processing timewindow %d from %d\n', j, ntoi);
 
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
progress('close');

%---create output-structure
freq          = [];
freq.label    = label;
freq.freq     = cfg.foi;
freq.time     = data.time;
%freq.cumtapcnt= ones(ntrl, 1)*ntap;
freq.dimord    = 'chan_chan_freq_time';
freq.transfer  = h;
freq.itransfer = a;
freq.noisecov  = data.noisecov;
freq.crsspctrm = crsspctrm;
freq.dof       = data.dof;

try,
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
for k = 1:nfoi
  h(:,:,k) = inv(zar(:,:,k)); 
end
h   = sqrt(2).*h; %account for the negative frequencies, normalization necessary for
%comparison with non-parametric (fft based) results in fieldtrip
zar = zar./sqrt(2);
