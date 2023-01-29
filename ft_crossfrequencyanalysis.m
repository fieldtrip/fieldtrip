function crossfreq = ft_crossfrequencyanalysis(cfg, freqlow, freqhigh)

% FT_CROSSFREQUENCYANALYSIS performs cross-frequency analysis
%
% Use as
%   crossfreq = ft_crossfrequencyanalysis(cfg, freq)
%   crossfreq = ft_crossfrequencyanalysis(cfg, freqlo, freqhi)
%
% The input data should be organised in a structure as obtained from the
% FT_FREQANALYSIS function. The configuration should be according to
%
%   cfg.freqlow    = scalar or vector, selection of frequencies for the low frequency data
%   cfg.freqhigh   = scalar or vector, selection of frequencies for the high frequency data
%
% Channel selection can be specified according to whether one wants to perform within- or
% cross-channel analysis.
%
% For within-channel analysis (default), you should specifies only a single channel selection:
%   cfg.channel    = cell-array with selection of channels, see FT_CHANNELSELECTION
% In this case, the output "dimord" will be "chan_freqlow_freqhigh"
%
% For cross-channel analysis, you should specifies two channel selections:
%   cfg.chanlow    = cell-array with selection of channels for the phase providing channels from the
%                    freqlow data argument, with wildcards allowed, see FT_CHANNELSELECTION
%   cfg.chanhigh   = cell-array with selection of channels for the amplitude providing channels from the
%                    freqhigh data argument, with wildcards allowed, see FT_CHANNELSELECTION
% In this case, the output "dimord" will be "chancmb_freqlow_freqhigh" and "label"
% field will be replaced with "labelcmb" (corresponding to the dimension "chancmb")
% describing the pairs of channel combinations as
%   {'chanlow01' 'chanhigh01'
%    'chanlow01' 'chanhigh02'
%    ...
%    'chanlow02' 'chanhigh01'
%    'chanlow02' 'chanhigh02'
%    ...
%    }
% N.B.: The order of channels corresponds to their order in the original "label" field
%
% Various metrics for cross-frequency coupling have been introduced in a number of
% scientific publications, but these do not use a consistent method naming scheme,
% nor implement it in exactly the same way. The particular implementation in this
% code tries to follow the most common format, generalizing where possible. If you
% want details about the algorithms, please look into the code.
%   cfg.method     = string, can be
%                     'coh' - coherence
%                     'plv' - phase locking value
%                     'mvl' - mean vector length
%                     'mi'  - modulation index
%                     'pac' - phase amplitude coupling
%
% The modulation index and phase amplitude coupling implement
%   Tort A. B. L., Komorowski R., Eichenbaum H., Kopell N. (2010). Measuring Phase-Amplitude
%   Coupling Between Neuronal Oscillations of Different Frequencies. J Neurophysiol 104:
%   1195?1210. doi:10.1152/jn.00106.2010
%
% cfg.keeptrials = string, can be 'yes' or 'no'
%
% See also FT_FREQANALYSIS, FT_CONNECTIVITYANALYSIS

% Copyright (C) 2014-2017, Donders Centre for Cognitive Neuroimaging
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
ft_preamble loadvar freqlow freqhigh
ft_preamble provenance freqlow freqhi

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

if nargin<3
  % use the same data for the low and high frequencies
  freqhigh = freqlow;
end

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was read from
% an old *.mat file
freqlow  = ft_checkdata(freqlow,  'datatype', 'freq', 'feedback', 'yes');
freqhigh = ft_checkdata(freqhigh, 'datatype', 'freq', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729

% FIXME the below is a bit hacky but it does the trick
if isfield(cfg, 'chanlow') && isfield(cfg, 'chanhigh')
  docrosschan   = true;
  cfg.chanlow   = ft_channelselection(cfg.chanlow, freqlow.label);
  cfg.chanhigh  = ft_channelselection(cfg.chanhigh, freqhigh.label);
  labelcmb = ft_channelcombination({cfg.chanlow,cfg.chanhigh},union(freqlow.label, freqhigh.label),1,2);
  %labelcmb(ismember(labelcmb(:,1),cfg.chanhigh)&strcmp(labelcmb(:,1),labelcmb(:,2)),:) = [];
elseif ~isfield(cfg, 'chanlow') && ~isfield(cfg, 'chanhigh')  % within-channel analysis (default)
  docrosschan = false;
  % ensure that we are working on the intersection of the channels
  cfg.channel  = ft_getopt(cfg, 'channel',  'all');
  cfg.channel  = ft_channelselection(cfg.channel, intersect(freqlow.label, freqhigh.label));
  cfg.chanlow  = cfg.channel;
  cfg.chanhigh = cfg.channel;
  labelcmb = horzcat(cfg.channel,cfg.channel);
else
  ft_error('you should either specify both cfg.chanlow and cfg.chanhigh, or none of these options');
end

% get the defaults
cfg.freqlow    = ft_getopt(cfg, 'freqlow',  'all');
cfg.freqhigh   = ft_getopt(cfg, 'freqhigh', 'all');
cfg.nphase     = ft_getopt(cfg, 'nphase', 20);
cfg.keeptrials = ft_getopt(cfg, 'keeptrials');

% make selection of frequencies and channels
tmpcfg = [];
tmpcfg.channel   = unique(labelcmb(:,1));
tmpcfg.frequency = cfg.freqlow;
freqlow = ft_selectdata(tmpcfg, freqlow);
[tmpcfg, freqlow] = rollback_provenance(cfg, freqlow);
try, cfg.freqlow = tmpcfg.frequency; end

% make selection of frequencies and channels
tmpcfg = [];
tmpcfg.channel   = unique(labelcmb(:,2));
tmpcfg.frequency = cfg.freqhigh;
freqhigh = ft_selectdata(tmpcfg, freqhigh);
[tmpcfg, freqhigh] = rollback_provenance(cfg, freqhigh);
try, cfg.freqhigh = tmpcfg.frequency; end

LF = freqlow.freq;
HF = freqhigh.freq;
ntrial = size(freqlow.fourierspctrm,1); % FIXME the dimord might be different
nchan  = size(labelcmb,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.method
  
  case 'coh'
    % coherence
    cohdatas = zeros(ntrial,nchan,numel(LF),numel(HF));
    for  i =1:nchan
      chandataLF = freqlow.fourierspctrm(:,strcmp(freqlow.label,labelcmb{i,1}),:,:);
      chandataHF = freqhigh.fourierspctrm(:,strcmp(freqhigh.label,labelcmb{i,2}),:,:);
      for j = 1:ntrial
        cohdatas(j,i,:,:) = data2coh(squeeze(chandataLF(j,:,:,:)),squeeze(chandataHF(j,:,:,:)));
      end
    end
    cfcdata = cohdatas;
    
  case 'plv'
    % phase locking value
    plvdatas = zeros(ntrial,nchan,numel(LF),numel(HF));
    for  i =1:nchan
      chandataLF = freqlow.fourierspctrm(:,strcmp(freqlow.label,labelcmb{i,1}),:,:);
      chandataHF = freqhigh.fourierspctrm(:,strcmp(freqhigh.label,labelcmb{i,2}),:,:);
      for j = 1:ntrial
        plvdatas(j,i,:,:) = data2plv(squeeze(chandataLF(j,:,:,:)),squeeze(chandataHF(j,:,:,:)));
      end
    end
    cfcdata = plvdatas;
    
  case  'mvl'
    % mean vector length
    mvldatas = zeros(ntrial,nchan,numel(LF),numel(HF));
    for  i =1:nchan
      chandataLF = freqlow.fourierspctrm(:,strcmp(freqlow.label,labelcmb{i,1}),:,:);
      chandataHF = freqhigh.fourierspctrm(:,strcmp(freqhigh.label,labelcmb{i,2}),:,:);
      for j = 1:ntrial
        mvldatas(j,i,:,:) = data2mvl(squeeze(chandataLF(j,:,:,:)),squeeze(chandataHF(j,:,:,:)));
      end
    end
    cfcdata = mvldatas;
    
  case  {'mi','pac'}
    % modulation index
    pacdatas   = zeros(ntrial,nchan,numel(LF),numel(HF),cfg.nphase);
    for  i =1:nchan
      chandataLF = freqlow.fourierspctrm(:,strcmp(freqlow.label,labelcmb{i,1}),:,:);
      chandataHF = freqhigh.fourierspctrm(:,strcmp(freqhigh.label,labelcmb{i,2}),:,:);
      for j = 1:ntrial
        [pacdatas(j,i,:,:,:), phasebins] = data2pac(squeeze(chandataLF(j,:,:,:)),squeeze(chandataHF(j,:,:,:)),cfg.nphase);
      end
    end
    cfcdata = pacdatas;
    
end % switch method for data preparation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the actual computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch cfg.method
  
  case 'coh'
    [ntrial,nchan,nlf,nhf] = size(cfcdata);
    if strcmp(cfg.keeptrials, 'no')
      crsspctrm = reshape(abs(mean(cfcdata,1)), [nchan, nlf, nhf]);
      dimord = 'chan_freqlow_freqhigh' ;
    else
      crsspctrm = abs(cfcdata);
      dimord = 'rpt_chan_freqlow_freqhigh' ;
    end
    
  case 'plv'
    [ntrial,nchan,nlf,nhf] = size(cfcdata);
    if strcmp(cfg.keeptrials, 'no')
      crsspctrm = reshape(abs(mean(cfcdata,1)), [nchan, nlf, nhf]);
      dimord = 'chan_freqlow_freqhigh' ;
    else
      crsspctrm = abs(cfcdata);
      dimord = 'rpt_chan_freqlow_freqhigh' ;
    end
    
  case  'mvl'
    [ntrial,nchan,nlf,nhf] = size(cfcdata);
    if strcmp(cfg.keeptrials, 'no')
      crsspctrm = reshape(abs(mean(cfcdata,1)), [nchan, nlf, nhf]);
      dimord = 'chan_freqlow_freqhigh' ;
    else
      crsspctrm = abs(cfcdata);
      dimord = 'rpt_chan_freqlow_freqhigh' ;
    end
    
  case  'mi'
    [ntrial,nchan,nlf,nhf,nbin] = size(cfcdata);
    
    if strcmp(cfg.keeptrials, 'yes')
      dimord = 'rpt_chan_freqlow_freqhigh' ;
      crsspctrm = zeros(ntrial,nchan,nlf,nhf);
      for k =1:ntrial
        for n=1:nchan
          pac = squeeze(cfcdata(k,n,:,:,:));
          Q =ones(nbin,1)/nbin;                             % uniform distribution
          mi = zeros(nlf,nhf);
          
          for i=1:nlf
            for j=1:nhf
              P = squeeze(pac(i,j,:))/ nansum(pac(i,j,:));  % normalized distribution
              % KL distance
              mi(i,j) = nansum(P.* log2(P./Q))./log2(nbin);
            end
          end
          crsspctrm(k,n,:,:) = mi;
          
        end
      end
      
    else
      dimord = 'chan_freqlow_freqhigh' ;
      crsspctrm = zeros(nchan,nlf,nhf);
      cfcdatamean = reshape(mean(cfcdata,1),[nchan nlf nhf nbin 1]);
      
      for k =1:nchan
        pac = squeeze(cfcdatamean(k,:,:,:));
        Q =ones(nbin,1)/nbin;                             % uniform distribution
        mi = zeros(nlf,nhf);
        
        for i=1:nlf
          for j=1:nhf
            P = squeeze(pac(i,j,:))/ nansum(pac(i,j,:));  % normalized distribution
            % KL distance
            mi(i,j) = nansum(P.* log2(P./Q))./log2(nbin);
          end
        end
        crsspctrm(k,:,:) = mi;
      end
      
    end % if keeptrials
    
  case 'pac'
    [ntrial,nchan,nlf,nhf,nbin] = size(cfcdata);
    
    if strcmp(cfg.keeptrials, 'yes')
      dimord = 'rpt_chan_freqlow_freqhigh_phase' ;
      crsspctrm = cfcdata;
      
    else
      dimord = 'chan_freqlow_freqhigh_phase' ;
      crsspctrm = reshape(mean(cfcdata,1),[nchan nlf nhf nbin 1]);
      crsspctrm(isnan(crsspctrm)) = 0;
      
    end % if keeptrials
    
end % switch method for actual computation

crossfreq.crsspctrm  = crsspctrm;
crossfreq.dimord     = dimord;
crossfreq.freqlow    = LF;
crossfreq.freqhigh   = HF;

if any(strcmp(strsplit(dimord,'_'),'phase'))
  crossfreq.phase = phasebins;
end

if docrosschan
  crossfreq.labelcmb = labelcmb;
  crossfreq.dimord   = strrep(crossfreq.dimord,'chan','chancmb');
else
  crossfreq.label    = cfg.channel;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   freqlow freqhigh
ft_postamble provenance crossfreq
ft_postamble history    crossfreq
ft_postamble savevar    crossfreq

end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cohdata] = data2coh(LFsigtemp,HFsigtemp)

HFamp    = abs(HFsigtemp);
HFamp(isnan(HFamp(:))) = 0;                              % replace nan with 0
HFphas   = angle(hilbert(HFamp'))';
HFsig    = HFamp .* exp(sqrt(-1)*HFphas);

LFsig = LFsigtemp;
LFsig(isnan(LFsig(:))) = 0;                              % replace nan with 0

cohdata = zeros(size(LFsig,1),size(HFsig,1));
for i = 1:size(LFsig,1)
  for j = 1:size(HFsig,1)
    Nx  = sum(~isnan(LFsigtemp(i,:) .* LFsigtemp(i,:)));
    Ny  = sum(~isnan(HFsigtemp(j,:) .* HFsigtemp(j,:)));
    Nxy = sum(~isnan(LFsigtemp(i,:) .* HFsigtemp(j,:)));
    
    Px  = LFsig(i,:) * ctranspose(LFsig(i,:)) ./ Nx;
    Py  = HFsig(j,:) * ctranspose(HFsig(j,:)) ./ Ny;
    Cxy = LFsig(i,:) * ctranspose(HFsig(j,:)) ./ Nxy;
    
    cohdata(i,j) = Cxy / sqrt(Px * Py);
  end
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plvdata] = data2plv(LFsigtemp,HFsigtemp)

LFphas   = angle(LFsigtemp);
HFamp    = abs(HFsigtemp);
HFamp(isnan(HFamp(:))) = 0;                              % replace nan with 0
HFphas   = angle(hilbert(HFamp'))';
plvdata  = zeros(size(LFsigtemp,1),size(HFsigtemp,1));   % phase locking value

for i = 1:size(LFsigtemp,1)
  for j = 1:size(HFsigtemp,1)
    plvdata(i,j) = nanmean(exp(1i*(LFphas(i,:)-HFphas(j,:))));
  end
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mvldata] = data2mvl(LFsigtemp,HFsigtemp)
% calculate  mean vector length (complex value) per trial
% mvldata dim: LF*HF

LFphas   = angle(LFsigtemp);
HFamp    = abs(HFsigtemp);
mvldata  = zeros(size(LFsigtemp,1),size(HFsigtemp,1));    % mean vector length

for i = 1:size(LFsigtemp,1)
  for j = 1:size(HFsigtemp,1)
    mvldata(i,j) = nanmean(HFamp(j,:).*exp(1i*LFphas(i,:)));
  end
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pacdata, phasebins] = data2pac(LFsigtemp,HFsigtemp,nbin)
% calculate phase amplitude distribution per trial
% pacdata dim: LF*HF*Phasebin

pacdata = zeros(size(LFsigtemp,1),size(HFsigtemp,1),nbin);

Ang  = angle(LFsigtemp);
Amp  = abs(HFsigtemp);
phasebins = linspace(-pi,pi,nbin);

% histc takes the edges rather than the centres of the bins
phasebinedges = (2*pi)/(nbin-1)/2;
phasebinedges = linspace(-pi-phasebinedges,pi+phasebinedges,nbin+1);

[dum,bin] = histc(Ang, phasebinedges);  % binned low frequency phase
binamp = zeros (size(HFsigtemp,1),nbin);      % binned amplitude

for i = 1:size(Ang,1)
  for k = 1:nbin
    idx = (bin(i,:)==k);
    binamp(:,k) = mean(Amp(:,idx),2);
  end
  pacdata(i,:,:) = binamp;
end

end % function
