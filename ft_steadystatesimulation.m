function [data] = ft_steadystatesimulation(cfg)

% FT_STEADYSTATESIMULATION creates a simulated EEG/MEG dataset. This function
% allows to simulate the effect of several independent stimulus trains. These can
% be presented as a periodic sequence, or as single (or few) transient stimuli.
% This function creates a single block of data. You can call it repeatedly and use
% FT_APPENDDATA to combine different blocks.
%
% Use as
%   data = ft_steadystatesimulation(cfg)
% where cfg is a configuration structure that should contain
%   cfg.fsample   = scalar, sampling frequency in Hz (default = 512)
%   cfg.duration  = scalar, trial length in seconds (default = 4.56)
%   cfg.baseline  = scalar, baseline length in seconds (default = 0)
%   cfg.ntrials   = integer N, number of trials (default = 320)
%   cfg.iti       = scalar, inter-trial interval in seconds (default = 1)
%
% Each trial can contain multiple nested experimental manipulations
%   cfg.level1.condition = scalar, or vector of length L1 (default = 1)
%   cfg.level1.gain      = scalar, or vector of length L1 (default = 1)
%   cfg.level2.condition = scalar, or vector of length L2 (default = 1)
%   cfg.level2.gain      = scalar, or vector of length L2 (default = 1)
%   cfg.level3.condition = scalar, or vector of length L3 (default = 1)
%   cfg.level3.gain      = scalar, or vector of length L3 (default = 1)
% If you don't need level 2 and up, specify the condition and gain as empty.
% Idem for level 3 and up.
%
% Stimuli are created at the lowest experimental level, and are modulated according to the product of the gain of all levels.
% Each trial can contain one or multiple stimuli.
% The behavior of each stimuli is specified with
%   cfg.stimulus1.mode = 'periodic', 'transient' or 'off' (default = 'periodic')
%   cfg.stimulus2.mode = 'periodic', 'transient' or 'off' (default = 'transient')
%
% If the stimulus is periodic (below as example for stimulus1), the following options apply
%   cfg.stimulus1.number          = does not apply for periodic stimuli
%   cfg.stimulus1.onset           = in seconds, first stimulus relative to the start of the trial (default = 0)
%   cfg.stimulus1.onsetjitter     = in seconds, max jitter that is added to the onset (default = 0)
%   cfg.stimulus1.isi             = in seconds, i.e. for 10Hz you would specify 0.1 seconds as the interstimulus interval (default = 0.1176)
%   cfg.stimulus1.isijitter       = in seconds, max jitter relative to the previous stimulus (default = 0)
%   cfg.stimulus2.condition       = does not apply for periodic stimuli
%   cfg.stimulus2.gain            = does not apply for periodic stimuli
%   cfg.stimulus1.kernelshape     = 'sine'
%   cfg.stimulus1.kernelduration  = in seconds (default = isi)
%
% If the stimulus is transient (below as example for stimulus2), the following options apply
%   cfg.stimulus2.number          = scalar M, how many transients are to be presented per trial (default = 4)
%   cfg.stimulus2.onset           = in seconds, first stimulus relative to the start of the trial (default = 0.7)
%   cfg.stimulus2.onsetjitter     = in seconds, max jitter that is added to the onset (default = 0.2)
%   cfg.stimulus2.isi             = in seconds as the interstimulus interval (default = 0.7)
%   cfg.stimulus2.isijitter       = in seconds, max jitter relative to the previous stimulus ((default = 0.2)
%   cfg.stimulus2.condition       = 1xM vector with condition codes for each transient within a trial (default = [1 1 2 2])
%   cfg.stimulus2.gain            = 1xM vector with gain for each condition for each transient within a trial(default = [1 1 1 1])
%   cfg.stimulus2.kernelshape     = 'hanning'
%   cfg.stimulus2.kernelduration  = in seconds (default = 0.75*isi)
%
% RANDOMIZATIONS:
% - The onsetjitter is randomized between 0 and the value given, and is always added to the onset.
% - The isijitter is randomized between 0 and the value given, and is always added to the interstimulus interval (isi).
% - For periodic stimuli, which are constant within a trial, the condition code and gain are shuffled over all trials.
% - For transient stimuli, the condition code and gain are shuffled within each trial.
%
% Using the default settings, we model a peripherally presented flickering stimulus
% that appears at different excentricities together with a centrally presented
% transient stimulus that appears 4x per trial. To simulate the experiment described
% at , you have to call this 4 times with a different cfg.configuration and
% cfg.gain to model the task load and use FT_APPENDDATA to concatenate the trials. In
% this case cfg.condition models the factor "task load" (2 levels, low and high),
% cfg.stimulus1.condition models the factor "excentricity" (4 levels), and
% cfg.stimulation2.condition models the factor "stimulus type" (2 levels, non-target
% or target).
%
% See also FT_DIPOLESIMULATION, FT_TIMELOCKSIMULATION, FT_FREQSIMULATION,
% FT_CONNECTIVITYSIMULATION, FT_APPENDDATA

% Copyright (C) 2017, Stefan Wiens, Malina Szychowska, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% get the options or set the defaults
cfg.level1 = ft_getopt(cfg, 'level1');
cfg.level2 = ft_getopt(cfg, 'level2');
cfg.level3 = ft_getopt(cfg, 'level3');
cfg.level1.condition = ft_getopt(cfg.level1, 'condition', 1);
cfg.level2.condition = ft_getopt(cfg.level2, 'condition', []);
cfg.level3.condition = ft_getopt(cfg.level3, 'condition', []);
cfg.level1.gain = ft_getopt(cfg.level1, 'gain', 1);
cfg.level2.gain = ft_getopt(cfg.level2, 'gain', []);
cfg.level3.gain = ft_getopt(cfg.level3, 'gain', []);

cfg.ntrials   = ft_getopt(cfg, 'ntrials', 320);
cfg.fsample   = ft_getopt(cfg, 'fsample', 512);
cfg.duration  = ft_getopt(cfg, 'duration', 4.5);
cfg.baseline  = ft_getopt(cfg, 'baseline', 0);
cfg.iti       = ft_getopt(cfg, 'iti', 1);

cfg.stimulus1 = ft_getopt(cfg, 'stimulus1');
cfg.stimulus2 = ft_getopt(cfg, 'stimulus2');

cfg.stimulus1.mode            = ft_getopt(cfg.stimulus1, 'mode', 'periodic');
cfg.stimulus1.number          = ft_getopt(cfg.stimulus1, 'number', 0);
cfg.stimulus1.onset           = ft_getopt(cfg.stimulus1, 'onset', 0);
cfg.stimulus1.onsetjitter     = ft_getopt(cfg.stimulus1, 'onsetjitter', 0);
cfg.stimulus1.isi             = ft_getopt(cfg.stimulus1, 'isi', 0.1176);
cfg.stimulus1.isijitter       = ft_getopt(cfg.stimulus1, 'isijitter', 0);
cfg.stimulus1.condition       = ft_getopt(cfg.stimulus1, 'condition');
cfg.stimulus1.gain            = ft_getopt(cfg.stimulus1, 'gain');
cfg.stimulus1.kernelshape     = ft_getopt(cfg.stimulus1, 'kernelshape', 'sine');
cfg.stimulus1.kernelduration  = ft_getopt(cfg.stimulus1, 'kernelduration', cfg.stimulus1.isi);

cfg.stimulus2.mode            = ft_getopt(cfg.stimulus2, 'mode', 'transient');
cfg.stimulus2.number          = ft_getopt(cfg.stimulus2, 'number', 4);
cfg.stimulus2.onset           = ft_getopt(cfg.stimulus2, 'onset', 0.7);
cfg.stimulus2.onsetjitter     = ft_getopt(cfg.stimulus2, 'onsetjitter', 0);
cfg.stimulus2.isi             = ft_getopt(cfg.stimulus2, 'isi', 0.7);
cfg.stimulus2.isijitter       = ft_getopt(cfg.stimulus2, 'isijitter', 0);
cfg.stimulus2.condition       = ft_getopt(cfg.stimulus2, 'condition', [1 2 3 4]);
cfg.stimulus2.gain            = ft_getopt(cfg.stimulus2, 'gain', [1 1 1 1]);
cfg.stimulus2.kernelshape     = ft_getopt(cfg.stimulus2, 'kernelshape', 'hanning');
cfg.stimulus2.kernelduration  = ft_getopt(cfg.stimulus2, 'kernelduration', 0.75*cfg.stimulus2.isi);

ncondition1 = length(cfg.level1.condition);
ncondition2 = length(cfg.level2.condition);
ncondition3 = length(cfg.level3.condition);

switch cfg.stimulus1.mode
  case 'off'
    cfg.stimulus1.number = 0;
    cfg.stimulus1.condition = [];
    cfg.stimulus1.gain = [];
  case 'transient'
    assert(numel(cfg.stimulus1.condition)==cfg.stimulus1.number);
    assert(numel(cfg.stimulus1.gain)==cfg.stimulus1.number);
end

switch cfg.stimulus2.mode
  case 'off'
    cfg.stimulus2.number = 0;
    cfg.stimulus2.condition = [];
    cfg.stimulus2.gain = [];
  case 'transient'
    assert(numel(cfg.stimulus2.condition)==cfg.stimulus2.number);
    assert(numel(cfg.stimulus2.gain)==cfg.stimulus2.number);
end

% they cannot be the same, as that would confuse the trialinfo
assert(~isequal(cfg.stimulus1.mode, cfg.stimulus2.mode));

% determine which levels are being used
level1 = ncondition1>0;
level2 = ncondition2>0;
level3 = ncondition3>0;

% count the number of blocks and number of trials per block
nblock = max(ncondition1,1) * max(ncondition2,1) * max(ncondition3,1);
ntrial = cfg.ntrials / nblock;

assert(ntrial==round(ntrial), 'number of trials inconsistent');

% the time (followinf onset) and baseline are the same for each trial
baseline = (-cfg.baseline):(1/cfg.fsample):0;
baseline(end) = []; % remove the last sample
time = 0:(1/cfg.fsample):cfg.duration;

% construct a sine-wave kernel, feval fails on this one
sine = @(n) sin(2*pi*((1:n)-1)/n);

% construct the kernels for the response to each stimulus
if strcmp(cfg.stimulus1.kernelshape, 'sine')
  kernel1 = sine(round(cfg.stimulus1.kernelduration*cfg.fsample));
else
  kernel1 = feval(cfg.stimulus1.kernelshape, round(cfg.stimulus1.kernelduration*cfg.fsample));
end
if strcmp(cfg.stimulus2.kernelshape, 'sine')
  kernel2 = sine(round(cfg.stimulus2.kernelduration*cfg.fsample));
else
  kernel2 = feval(cfg.stimulus2.kernelshape, round(cfg.stimulus2.kernelduration*cfg.fsample));
end
% these should be row-vectors
kernel1 = kernel1(:)';
kernel2 = kernel2(:)';

% start with empty data
data = [];

% this keeps track of the trials
k = 1;

shuffle = randperm(ncondition3);
tmpcondition3 = cfg.level3.condition(shuffle);
tmpgain3      = cfg.level3.gain(shuffle);
for cond3=1:max(ncondition3, 1)
  if level3
    trialgain = tmpgain3(cond3);
  else
    trialgain = 1;
  end
  
  shuffle = randperm(ncondition2);
  tmpcondition2 = cfg.level2.condition(shuffle);
  tmpgain2      = cfg.level2.gain(shuffle);
  for cond2=1:max(ncondition2, 1)
    if level2
      trialgain = trialgain * tmpgain2(cond2);
    end
    
    shuffle = randperm(ncondition1);
    tmpcondition1 = cfg.level1.condition(shuffle);
    tmpgain1      = cfg.level1.gain(shuffle);
    for cond1=1:max(ncondition1, 1)
      if level1
        trialgain = trialgain * tmpgain1(cond1);
      end
      
      for trial=1:ntrial
        trialinfo = k; % start with the trial number
        varname = {'trial'};
        
        if level1
          trialinfo = [trialinfo tmpcondition1(cond1)];
          varname(end+1) = {'l1_cond'};
        end
        if level2>0
          trialinfo = [trialinfo tmpcondition2(cond2)];
          varname(end+1) = {'l2_cond'};
        end
        if level3>0
          trialinfo = [trialinfo tmpcondition3(cond3)];
          varname(end+1) = {'l3_cond'};
        end
        
        dat = [];
        
        stimcfg = cfg.stimulus1;
        signal  = zeros(size(time)); % the baseline will be added later
        
        switch stimcfg.mode
          case 'periodic'
            sample = nearest(time, stimcfg.onset + rand(1)*stimcfg.onsetjitter);
            for i=1:ceil((cfg.duration-stimcfg.onset)/stimcfg.isi)
              sample = sample + round(cfg.fsample*(stimcfg.isi + stimcfg.isijitter*rand(1)));
              if sample>length(signal)
                % this periodic stimulus falls outside the trial
                sample = nan;
              else
                signal(sample) = 1;
              end
              trialinfo = [trialinfo sample+length(baseline)];
              varname(end+1) = {sprintf('p%d_sample', i)};
            end % while
            
          case 'transient'
            shuffle = randperm(numel(stimcfg.gain));
            sample  = nearest(time, stimcfg.onset + rand(1)*stimcfg.onsetjitter);
            transientgain = stimcfg.gain(shuffle);
            transientcondition = stimcfg.condition(shuffle);
            for i=1:stimcfg.number
              sample = sample + round(cfg.fsample*(stimcfg.isi + stimcfg.isijitter*rand(1)));
              if sample>length(signal)
                warning('transient stimulus falls outside of trial');
                sample = nan;
              else
                signal(sample) = transientgain(i);
              end
              trialinfo = [trialinfo transientcondition(i)];
              varname(end+1) = {sprintf('t%d_cond', i)};
              trialinfo = [trialinfo sample+length(baseline)];
              varname(end+1) = {sprintf('t%d_sample', i)};
              
            end % while
        end % switch mode
        
        dat(2,:) = [zeros(size(baseline))             signal*trialgain          ];
        dat(3,:) = [zeros(size(baseline)) erpconvolve(signal*trialgain, kernel1)];
        
        stimcfg  = cfg.stimulus2;
        signal   = zeros(size(time));  % the baseline will be added later
        
        switch stimcfg.mode
          case 'periodic'
            sample = nearest(time, stimcfg.onset + rand(1)*stimcfg.onsetjitter);
            for i=1:ceil((cfg.duration-stimcfg.onset)/stimcfg.isi)
              sample = sample + round(cfg.fsample*(stimcfg.isi + stimcfg.isijitter*rand(1)));
              if sample>length(signal)
                % this periodic stimulus falls outside the trial
                sample = nan;
              else
                signal(sample) = 1;
              end
              trialinfo = [trialinfo sample+length(baseline)];
              varname(end+1) = {sprintf('p%d_sample', i)};
            end % while
            
          case 'transient'
            shuffle = randperm(numel(stimcfg.gain));
            sample  = nearest(time, stimcfg.onset + rand(1)*stimcfg.onsetjitter);
            transientgain = stimcfg.gain(shuffle);
            transientcondition = stimcfg.condition(shuffle);
            for i=1:stimcfg.number
              sample = sample + round(cfg.fsample*(stimcfg.isi + stimcfg.isijitter*rand(1)));
              if sample>length(signal)
                warning('transient stimulus falls outside of trial');
                sample = nan;
              else
                signal(sample) = transientgain(i);
              end
              trialinfo = [trialinfo transientcondition(i)];
              varname(end+1) = {sprintf('t%d_cond', i)};
              trialinfo = [trialinfo sample+length(baseline)];
              varname(end+1) = {sprintf('t%d_sample', i)};
            end % while
            
        end % switch mode
        
        dat(4,:) = [zeros(size(baseline))             signal*trialgain          ];
        dat(5,:) = [zeros(size(baseline)) erpconvolve(signal*trialgain, kernel2)];
        
        % the first channel contains the sum of the two convolved signals
        dat(1,:) = sum(dat([3 5],:), 1);
        
        data.trial{k} = dat;
        data.time{k}  = [baseline time];
        data.trialinfo(k,:) = trialinfo;
        
        k = k+1;
      end % for each trial in each block
    end % for all conditions in level 1
  end % for all conditions in level 2
end % for all conditions in level 3

data.trialinfo = array2table(data.trialinfo, 'VariableNames', varname);
data.fsample = cfg.fsample;
data.label = {
  'combined'
  'stickfunction1'
  'convolved1'
  'stickfunction2'
  'convolved2'
  };

% construct the sample information, start with consecutive segments
nsample = length(baseline) + length(time);
data.sampleinfo(:,1) = ((1:ntrial)-1)*nsample + 1;
data.sampleinfo(:,2) = ((1:ntrial)  )*nsample;
% add the inter-trial-interval
data.sampleinfo(:,1) = data.sampleinfo(:,1) + transpose((1:ntrial)-1)*round(cfg.iti*cfg.fsample);
data.sampleinfo(:,2) = data.sampleinfo(:,2) + transpose((1:ntrial)-1)*round(cfg.iti*cfg.fsample);

ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dat = erpconvolve(stickfunction, kernel)
dat = convn(stickfunction, kernel, 'full');
dat = dat(1:length(stickfunction));
