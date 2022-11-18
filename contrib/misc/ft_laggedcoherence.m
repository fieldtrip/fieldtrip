function [dataout] = ft_laggedcoherence(cfg, datain)

% FT_LAGGEDCOHERENCE calculates the lagged coherence for a set of
% channels pairs, a number of frequencies, and a number of lags. The channel
% pairs are the complete set of pairs that can be formed from the channels
% in datain, including the auto-pairs (one channel for which the lagged
% coherence is calculated at different lags).
%
% Use as
%   outdata = ft_laggedcoherence(cfg, indata)
% where cfg is a configuration structure (see below) and indata is the
% output of FT_PREPROCESSING.
%
% The configuration structure should contain
%   cfg.foi       = vector 1 x numfoi, frequencies of interest
%   cfg.loi       = vector 1 x numloi, lags of interest, this must be a vector of
%                   integers with starting value 0 or higher
%   cfg.numcycles = integer, number of cycles of the Fourier basis functions that
%                   are used to calculate the Fourier coefficients that are the
%                   basis for calculating lagged coherence
%
%
% When using the results of this function in a publication, please cite:
%   Fransen, A. M., van Ede, F., & Maris, E. (2015). Identifying neuronal
%   oscillations using rhythmicity. Neuroimage, 118, 256-267.
%
% See also FT_PREPROCESSING, FT_FREQANALYSIS, FT_CONNECTIVITYANALYSIS

% Copyright (C) 2019-2020, DCC, Eric Maris & Anne Fransen; DCCN, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
% FieldTrip is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% FieldTrip is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% agreement = {
%   'When using the results of this function in a publication, please cite:'
%   ''
%   'Fransen, A. M., van Ede, F., & Maris, E. (2015). Identifying neuronal'
%   'oscillations using rhythmicity. Neuroimage, 118, 256-267.'
%   };
%
% if ~strcmp(questdlg(agreement, 'User agreement', 'Yes', 'Cancel', 'Cancel'), 'Yes')
%   return
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin = nargin;
ft_nargout = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar datain
ft_preamble provenance datain


% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the required options are present
cfg.channel = ft_getopt(cfg, 'channel', 'all');
cfg.trials = ft_getopt(cfg, 'trials', 'all');

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
datain = ft_selectdata(tmpcfg, datain);
% restore the provenance information
[cfg, datain] = rollback_provenance(cfg, datain);

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
datain = ft_checkdata(datain, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'foi', 'loi', 'numcycles'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataout = [];

% Run a time-resolved frequency analysis (method = mtmconvol) to produce
% Fourier coefficients that will later be used for calculating lagged
% coherence. We begin by calculating the toi, foi and t_ftimwin configuration
% fields for the call to ft_freqanalysis.

Fs = datain.fsample;
minfoi = min(cfg.foi);
nsamplespercycleforminfoi = ceil(Fs/minfoi);
nsamplespercycleall = 2:nsamplespercycleforminfoi;
freqall = Fs*ones(size(nsamplespercycleall))./nsamplespercycleall;
nsamplespertimwinall = nsamplespercycleall*cfg.numcycles;
% remove all elements in nsamplespertimwinall that are even (because only
% an odd number of samples is consistent with a cfg.toi element that has
% (nsamplespertimwinall-1)/2 samples at each side.
oddnsamples = mod(nsamplespertimwinall,2) == 1;
nsamplespertimwinall = nsamplespertimwinall(oddnsamples);
freqall = freqall(oddnsamples);
timepertimwinall = nsamplespertimwinall/Fs;
% Find the best matching frequencies
freq = [];
t_ftimwin = [];
nsamplespertimwin = [];
for foiind = 1:length(cfg.foi)
  absdiff = abs(freqall-cfg.foi(foiind));
  [minval,minvalind] = min(absdiff);
  if isempty(freq) || ~any(freq == freqall(minvalind))
    freq(end+1) = freqall(minvalind);
    t_ftimwin(end+1) = timepertimwinall(minvalind);
    nsamplespertimwin(end+1) = nsamplespertimwinall(minvalind);
  end
end

% Construct a frequency-specific toi vector for every element in freq,
% which will later be passed as an argument to ft_freqanalysis.
% For the moment, we assume that all trials in datain have the same time
% vector. This may have to be generalized.
timevec = datain.time{1};
ok = true;
for k = 1:numel(datain.trial)
  ok = ok && isequal(datain.time{k},timevec);
end
if ~ok
  ft_error('the input data should have the same time axis on each trial');
end

% remove the frequencies with a nsamplespertimwin that is larger than the trial length.
remfreq = false(size(freq));
for freqindx = 1:length(freq)
  remfreq(freqindx) = nsamplespertimwin(freqindx)>numel(timevec);
end
freq=freq(~remfreq);
nsamplespertimwin=nsamplespertimwin(~remfreq);
t_ftimwin=t_ftimwin(~remfreq);


% construct the toi vectors for each of the frequencies
toicell = cell(size(freq));
for freqindx = 1:length(freq)
  nsegments = floor(length(timevec)/nsamplespertimwin(freqindx));
  nsamplesnext2toi = ceil(nsamplespertimwin(freqindx)./2); %(nsamplespertimwin(freqindx)-1)/2;
  toisamples = (nsamplesnext2toi+1):nsamplespertimwin(freqindx):nsegments*nsamplespertimwin(freqindx);
  toicell{freqindx} = timevec(toisamples);
end

% Build the cfg for ft_freqanalysis
cfg_freq = [];
cfg_freq.method = 'mtmconvol';
cfg_freq.output = 'fourier';
cfg_freq.keeptrials = 'yes';
cfg_freq.taper = 'hanning';

cfg_lcoh.method = 'laggedcoherence';
cfg_lcoh.channelcmb = ft_channelcombination({'all' 'all'},datain.label);

% Loop over the frequencies
freqout = cell(1,numel(freq));
cohout = cell(1,numel(freq));
for freqindx = 1:length(freq)
  cfg_freq.foi = freq(freqindx);
  cfg_freq.t_ftimwin = t_ftimwin(freqindx);
  cfg_freq.toi = toicell{freqindx};
  cfg_freq.pad = ceil(numel(timevec)./Fs.*cfg_freq.foi)./cfg_freq.foi;
%
% Code for checking whether the requested cfg_freq.t_ftimwin (calculated as t_ftimwin) allows for the
% requested cfg_freq.foi (calculated as freq) using the requested number of cycles (cfg.numcycles).
%
% cyclelengths=t_ftimwin/cfg.numcycles
% correspondingfreqs=ones(size(cyclelengths))./cyclelengths
% freq % compare with correspondingfreqs

  freqout{freqindx} = ft_freqanalysis(cfg_freq,datain);
  
  cfg_lcoh.laggedcoherence.lags = cfg.numcycles.*cfg.loi./freqout{freqindx}.freq;
  cohout{freqindx} = ft_connectivityanalysis(cfg_lcoh,freqout{freqindx});
  if freqindx == 1
    labelcmb = cohout{1}.labelcmb;
    label = cell(size(cohout{1}.labelcmb,1),1);
    for m = 1:numel(label)
      label{m} = [cohout{1}.labelcmb{m,1} '_' cohout{1}.labelcmb{m,2}];
    end
  end
  cohout{freqindx}.label = label;
  cohout{freqindx}.dimord = 'chan_freq_time';
  cohout{freqindx} = rmfield(cohout{freqindx},'labelcmb');
  cohout{freqindx}.time = round(cohout{freqindx}.time.*cohout{freqindx}.freq./cfg.numcycles);
end

cfg_append = [];
cfg_append.parameter = 'lcohspctrm';
dataout = ft_appendfreq(cfg_append,cohout{:});
dataout.labelcmb = labelcmb;
dataout.dimord = 'chancmb_freq_lag';
dataout.lag = dataout.time;
dataout = removefields(dataout,{'label','time'});

% this might involve more active checking of whether the input options
% are consistent with the data and with each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug

ft_postamble previous datain
ft_postamble provenance dataout
ft_postamble history dataout
ft_postamble savevar dataout
