function [data] = ft_freqsimulation(cfg)

% FT_FREQSIMULATION simulates channel-level time-series data . The data is built up
% from different frequencies and can contain a signal in which the different
% frequencies interact (i.e. cross-frequency coherent). Different methods are
% possible to make data with specific properties.
%
% Use as
%   [data] = ft_freqsimulation(cfg)
% which will return a raw data structure that resembles the output of
% FT_PREPROCESSING.
%
% The configuration options can include
%   cfg.method     = The methods are explained in more detail below, but they can be
%                     'superimposed'    simply add the contribution of the different frequencies
%                     'broadband'       create a single broadband signal component
%                     'phalow_amphigh'  phase of low freq correlated with amplitude of high freq
%                     'amplow_amphigh'  amplitude of low freq correlated with amplithude of high freq
%                     'phalow_freqhigh' phase of low freq correlated with frequency of high signal
%                     'asymmetric'      single signal component with asymmetric positive/negative deflections
%   cfg.output     = which channels should be in the output data, can be 'mixed' or 'all' (default = 'all')
%   cfg.randomseed = 'yes' or a number or vector with the seed value (default = 'yes')
%
% The number of trials and the time axes of the trials can be specified by
%   cfg.fsample    = simulated sample frequency (default = 1200)
%   cfg.trllen     = length of simulated trials in seconds (default = 1)
%   cfg.numtrl     = number of simulated trials (default = 1)
%   cfg.baseline   = number (default = 0)
% or by
%   cfg.time       = cell-array with one time axis per trial, which are for example obtained from an existing dataset
%
% For each of the methods default parameters are configured to generate
% example data, including noise. To get full control over the generated
% data you should explicitely set all parameters involved in the method
% of your choise. The interpretation of the following signal components
% depends on the specified method:
%
% cfg.s1.freq     = frequency of signal 1
% cfg.s1.phase    = phase (in rad) relative to cosine of signal 1  (default depends on method)
%                 = number or 'random'
% cfg.s1.ampl     = amplitude of signal 1
% cfg.s2.freq     = frequency of signal 2
% cfg.s2.phase    = phase (in rad) relative to cosine of signal 1  (default depends on method)
%                 = number or 'random'
% cfg.s2.ampl     = amplitude of signal 2
% cfg.s3.freq     = frequency of signal 3
% cfg.s3.phase    = phase (in rad) relative to cosine of signal 1  (default depends on method)
%                 = number or 'random'
% cfg.s3.ampl     = amplitude of signal 3
% cfg.s4.freq     = frequency of signal 4
% cfg.s4.phase    = phase (in rad) relative to cosine of signal 1  (default depends on method)
%                 = number or 'random'
% cfg.s4.ampl     = amplitude of signal 4
%
% cfg.n1.ampl     = root-mean-square amplitude of wide-band signal prior to filtering
% cfg.n1.bpfreq   = [Flow Fhigh]
% cfg.n2.ampl     = root-mean-square amplitude of wide-band signal prior to filtering
% cfg.n2.bpfreq   = [Flow Fhigh]
%
% cfg.asymmetry   = amount of asymmetry (default = 0, which is none)
% cfg.noise.ampl  = amplitude of noise
%
%
% In the method 'superimposed' the signal contains just the sum of the different frequency contributions:
%     s1: first frequency
%     s2: second frequency
%     s3: third frequency
% and the output consists of the following channels:
%     1st channel: mixed signal = s1 + s2 + s3 + noise
%     2nd channel: s1
%     3rd channel: s2
%     4th channel: s3
%     5th channel: noise
%
% In the method 'broadband' the signal contains a the superposition of two
% broadband signal components, which are created by bandpass filtering a
% Gaussian noise signal:
%     n1: first broadband signal
%     n2: second broadband signal
% and the output consists of the following channels:
%     1st channel: mixed signal = n1 + n2 + noise
%     2nd channel: n1
%     3rd channel: n2
%     4th channel: noise
%
% In the method 'phalow_amphigh' the signal is build up of 4 components; s1, s2, s3 and noise:
%     s1: amplitude modulation (AM), frequency of this signal should be lower than s2
%     s2: second frequency, frequncy that becomes amplitude modulated
%     s3: DC shift of s1, should have frequency of 0
% and the output consists of the following channels:
%     1st channel: mixed signal = (s1 + s3)*s2 + noise,
%     2nd channel: s1
%     3rd channel: s2
%     4th channel: s3
%     5th channel: noise
%
% In the method 'amplow_amphigh' the signal is build up of 5 components; s1, s2, s3, s4 and noise.
%     s1: first frequency
%     s2: second frequency
%     s3: DC shift of s1 and s2, should have frequency of 0
%     s4: amplitude modulation (AM), frequency of this signal should be lower than s1 and s2
% and the output consists of the following channels:
%     1st channel: mixed signal = (s4 + s3)*s1 + (s4 + s3)*s2 + noise,
%     2nd channel: s1
%     3rd channel: s2
%     4th channel: s3
%     5th channel: noise
%     6th channel: s4
%     7th channel: mixed part 1: (s4 + s3)*s1
%     8th channel: mixed part 2: (s4 + s3)*s2
%
% In the method 'phalow_freqhigh' a frequency modulated signal is created.
%   signal is build up of 3 components; s1, s2 and noise.
%     s1: represents the base signal that will be modulated
%     s2: signal that will be used for the frequency modulation
% and the output consists of the following channels:
%     1st channel: mixed signal = s1.ampl * cos(ins_pha) + noise
%     2nd channel: s1
%     3rd channel: s2
%     4th channel: noise
%     5th channel: inst_pha_base   instantaneous phase of the high (=base) frequency signal s1
%     6th channel: inst_pha_mod    low frequency phase modulation, this is equal to s2
%     7th channel: inst_pha        instantaneous phase, i.e. inst_pha_base + inst_pha_mod
%
% In the method 'asymmetric' there is only one periodic signal, but that
% signal is more peaked for the positive than for the negative deflections.
% The average of the signal over time is zero.
%     s1: represents the frequency of the base signal
% and the output consists of the following channels:
%     1st channel: mixed signal = asymmetric signal + noise
%     2nd channel: sine wave with base frequency and phase, i.e. s1
%     3rd channel: asymmetric signal
%     4th channel: noise
%
% See also FT_FREQANALYSIS, FT_TIMELOCKSIMULATION, FT_DIPOLESIMULATION,
% FT_CONNECTIVITYSIMULATION

% Copyright (C) 2007-2008, Ingrid Nieuwenhuis & Robert Oostenveld
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
ft_preamble provenance
ft_preamble randomseed

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% return immediately after distributed execution
if ~isempty(ft_getopt(cfg, 'distribute'))
  return
end

% set defaults
if ~isfield(cfg, 'method'),        cfg.method = 'phalow_amphigh';         end
if ~isfield(cfg, 'output'),        cfg.output = 'all';                    end
if ~isfield(cfg, 'time'),          cfg.time = [];                         end

if isempty(cfg.time)
  cfg.fsample   = ft_getopt(cfg, 'fsample', 1200);
  cfg.trllen    = ft_getopt(cfg, 'trllen', 1);
  cfg.numtrl    = ft_getopt(cfg, 'numtrl', 1);
  cfg.baseline  = ft_getopt(cfg, 'baseline', 0);
else
  cfg.trllen  = [];                         % can be variable
  cfg.fsample = 1/mean(diff(cfg.time{1}));  % determine from time-axis
  cfg.numtrl  = length(cfg.time);
end

if strcmp(cfg.method, 'superimposed')
  if ~isfield(cfg, 's1'),           cfg.s1 = [];                          end
  if ~isfield(cfg.s1, 'freq'),      cfg.s1.freq  = 10;                    end
  if ~isfield(cfg.s1, 'phase'),     cfg.s1.phase = 0;                     end
  if ~isfield(cfg.s1, 'ampl'),      cfg.s1.ampl  = 1;                     end
  if ~isfield(cfg, 's2'),           cfg.s2 = [];                          end
  if ~isfield(cfg.s2, 'freq'),      cfg.s2.freq  = 20;                    end
  if ~isfield(cfg.s2, 'phase'),     cfg.s2.phase = 0;                     end
  if ~isfield(cfg.s2, 'ampl'),      cfg.s2.ampl  = 0;                     end
  if ~isfield(cfg, 's3'),           cfg.s3 = [];                          end
  if ~isfield(cfg.s3, 'freq'),      cfg.s3.freq  = 30;                    end
  if ~isfield(cfg.s3, 'phase'),     cfg.s3.phase = 0;                     end
  if ~isfield(cfg.s3, 'ampl'),      cfg.s3.ampl  = 0;                     end
end

if strcmp(cfg.method, 'broadband')
  if ~isfield(cfg, 'n1'),           cfg.n1 = [];                          end
  if ~isfield(cfg.n1, 'ampl'),      cfg.n1.ampl  = 1;                     end
  if ~isfield(cfg.n1, 'bpfreq'),    cfg.n1.bpfreq  = [30 50];             end
  if ~isfield(cfg, 'n2'),           cfg.n2 = [];                          end
  if ~isfield(cfg.n2, 'ampl'),      cfg.n2.ampl  = 1;                     end
  if ~isfield(cfg.n2, 'bpfreq'),    cfg.n2.bpfreq  = [80 120];            end
end

if strcmp(cfg.method, 'phalow_amphigh')
  if ~isfield(cfg, 's1'),           cfg.s1 = [];                          end
  if ~isfield(cfg.s1, 'freq'),      cfg.s1.freq = 3;                      end
  if ~isfield(cfg.s1, 'phase'),     cfg.s1.phase = -1*pi;                 end
  if ~isfield(cfg.s1, 'ampl'),      cfg.s1.ampl = 1;                      end
  if ~isfield(cfg, 's2'),           cfg.s2 = [];                          end
  if ~isfield(cfg.s2, 'freq'),      cfg.s2.freq = 20;                     end
  if ~isfield(cfg.s2, 'phase'),     cfg.s2.phase = 0;                     end
  if ~isfield(cfg.s2, 'ampl'),      cfg.s2.ampl = 1;                      end
  if ~isfield(cfg, 's3'),           cfg.s3 = [];                          end
  if ~isfield(cfg.s3, 'freq'),      cfg.s3.freq = 0;                      end
  if ~isfield(cfg.s3, 'phase'),     cfg.s3.phase = 0;                     end
  if ~isfield(cfg.s3, 'ampl'),      cfg.s3.ampl = cfg.s1.ampl;            end
end

if strcmp(cfg.method, 'amplow_amphigh')
  if ~isfield(cfg, 's1'),           cfg.s1 = [];                          end
  if ~isfield(cfg.s1, 'freq'),      cfg.s1.freq = 6;                      end
  if ~isfield(cfg.s1, 'phase'),     cfg.s1.phase = 0;                     end
  if ~isfield(cfg.s1, 'ampl'),      cfg.s1.ampl = 1;                      end
  if ~isfield(cfg, 's2'),           cfg.s2 = [];                          end
  if ~isfield(cfg.s2, 'freq'),      cfg.s2.freq = 20;                     end
  if ~isfield(cfg.s2, 'phase'),     cfg.s2.phase = 0;                     end
  if ~isfield(cfg.s2, 'ampl'),      cfg.s2.ampl = 1;                      end
  if ~isfield(cfg, 's4'),           cfg.s4 = [];                          end
  if ~isfield(cfg.s4, 'freq'),      cfg.s4.freq = 1;                      end
  if ~isfield(cfg.s4, 'phase'),     cfg.s4.phase = -1*pi;                 end
  if ~isfield(cfg.s4, 'ampl'),      cfg.s4.ampl = 1;                      end
  if ~isfield(cfg, 's3'),           cfg.s3 = [];                          end
  if ~isfield(cfg.s3, 'freq'),      cfg.s3.freq = 0;                      end
  if ~isfield(cfg.s3, 'phase'),     cfg.s3.phase = 0;                     end
  if ~isfield(cfg.s3, 'ampl'),      cfg.s3.ampl = cfg.s4.ampl;            end
end

if strcmp(cfg.method, 'phalow_freqhigh')
  if ~isfield(cfg, 's1'),           cfg.s1 = [];                          end
  if ~isfield(cfg.s1, 'freq'),      cfg.s1.freq = 20;                     end
  if ~isfield(cfg.s1, 'phase'),     cfg.s1.phase = 0;                     end
  if ~isfield(cfg.s1, 'ampl'),      cfg.s1.ampl = 1;                      end
  if ~isfield(cfg, 's2'),           cfg.s2 = [];                          end
  if ~isfield(cfg.s2, 'freq'),      cfg.s2.freq = 2;                      end
  if ~isfield(cfg.s2, 'phase'),     cfg.s2.phase = -0.5 * pi;             end %then base freq at t=0
  if ~isfield(cfg.s2, 'ampl'),      cfg.s2.ampl = pi;                     end
end

if strcmp(cfg.method, 'asymmetric')
  if ~isfield(cfg, 's1'),           cfg.s1 = [];                          end
  if ~isfield(cfg.s1, 'freq'),      cfg.s1.freq = 6;                      end
  if ~isfield(cfg.s1, 'phase'),     cfg.s1.phase = 0;                     end
  if ~isfield(cfg.s1, 'ampl'),      cfg.s1.ampl = 1;                      end
  if ~isfield(cfg, 'noise'),        cfg.noise = [];                       end
  if ~isfield(cfg.noise, 'ampl'),   cfg.noise.ampl = 0.1;                 end % default should not be too high
end

if ~isfield(cfg, 'noise'),         cfg.noise = [];                        end
if ~isfield(cfg.noise, 'ampl'),    cfg.noise.ampl = 1;                    end

if ~isempty(cfg.time)
  % use the user-supplied time vectors
  timevec = cfg.time;
else
  % give the user some feedback
  ft_debug('using %f as samping frequency', cfg.fsample);
  ft_debug('using %d trials of %f seconds long', cfg.numtrl, cfg.trllen);
  nsample = round(cfg.trllen*cfg.fsample);
  timevec = cell(1, cfg.numtrl);
  for iTr = 1:cfg.numtrl
    timevec{iTr} = (((1:nsample)-1)/cfg.fsample) - cfg.baseline;
  end
end

% give the user some feedback
ft_info('simulating data using %s method', cfg.method);

%%%%%%% SUPERIMPOSED, SIMPLY ADD THE SIGNALS %%%%%%%%%
if strcmp(cfg.method, 'superimposed')
  
  % make data
  for iTr = 1:length(timevec)
    if ischar(cfg.s1.phase); phase_s1 = rand * 2 * pi; else phase_s1 = cfg.s1.phase; end
    if ischar(cfg.s2.phase); phase_s2 = rand * 2 * pi; else phase_s2 = cfg.s2.phase; end
    if ischar(cfg.s3.phase); phase_s3 = rand * 2 * pi; else phase_s3 = cfg.s3.phase; end
    
    s1    = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_s1);
    s2    = cfg.s2.ampl*cos(2*pi*cfg.s2.freq*timevec{iTr} + phase_s2);
    s3    = cfg.s3.ampl*cos(2*pi*cfg.s3.freq*timevec{iTr} + phase_s3);
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix   = s1 + s2 + s3 + noise;
    
    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output, 'all')
      data.trial{iTr}(2,:) = s1;
      data.trial{iTr}(3,:) = s2;
      data.trial{iTr}(4,:) = s3;
      data.trial{iTr}(5,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr
  
  data.label{1} = 'mix';
  if strcmp(cfg.output, 'all')
    data.label{2} = 's1';
    data.label{3} = 's2';
    data.label{4} = 's3';
    data.label{5} = 'noise';
  end
  data.fsample = cfg.fsample;
  
  %%%%%%% SUPERIMPOSED BROADBAND SIGNAL %%%%%%%%%
elseif strcmp(cfg.method, 'broadband')
  
  % make data
  for iTr = 1:length(timevec)
    n1    = ft_preproc_bandpassfilter(cfg.n1.ampl*randn(size(timevec{iTr})), cfg.fsample, cfg.n1.bpfreq);
    n2    = ft_preproc_bandpassfilter(cfg.n2.ampl*randn(size(timevec{iTr})), cfg.fsample, cfg.n2.bpfreq);
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix   = n1 + n2 + noise;
    
    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output, 'all')
      data.trial{iTr}(2,:) = n1;
      data.trial{iTr}(3,:) = n2;
      data.trial{iTr}(4,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr
  
  data.label{1} = 'mix';
  if strcmp(cfg.output, 'all')
    data.label{2} = 'n1';
    data.label{3} = 'n2';
    data.label{4} = 'noise';
  end
  data.fsample = cfg.fsample;
  
  %%%%%%% PHASE TO AMPLITUDE CORRELATION %%%%%%%%%
elseif strcmp(cfg.method, 'phalow_amphigh')
  
  % sanity checks
  if cfg.s2.freq < cfg.s1.freq
    ft_error('with method is phalow_amphigh freq s2 should be higher than freq s1')
  end
  if cfg.s2.freq > cfg.fsample/2
    ft_error('you cannot have a frequency higher than the sample frequency/2')
  end
  if cfg.s3.freq ~= 0 || cfg.s3.phase ~= 0
    ft_warning('for method phalow_amphigh s3 is DC and therefore expect freq and phase to be zero but they are not')
  end
  if cfg.s3.ampl < cfg.s1.ampl
    ft_warning('expect amplitude s3 (=DC) not to be smaller than amplitude s1 (=low frequency)')
  end
  
  % make data
  for iTr = 1:length(timevec)
    
    if ischar(cfg.s1.phase); phase_AM   = rand * 2 * pi; else phase_AM   = cfg.s1.phase;  end
    if ischar(cfg.s2.phase); phase_high = rand * 2 * pi; else phase_high = cfg.s2.phase; end
    if ischar(cfg.s3.phase); phase_DC   = rand * 2 * pi; else phase_DC   = cfg.s3.phase;   end
    high  = cfg.s2.ampl*cos(2*pi*cfg.s2.freq*timevec{iTr} + phase_high);
    AM    = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_AM);
    DC    = cfg.s3.ampl*cos(2*pi*0*timevec{iTr} + phase_DC);
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix   = ((AM + DC) .* high) + noise;
    
    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output, 'all')
      data.trial{iTr}(2,:) = AM;
      data.trial{iTr}(3,:) = high;
      data.trial{iTr}(4,:) = DC;
      data.trial{iTr}(5,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr
  
  data.label{1} = 'mix';
  if strcmp(cfg.output, 'all')
    data.label{2} = 's1 (AM)';
    data.label{3} = 's2 (high)';
    data.label{4} = 's3 (DC)';
    data.label{5} = 'noise';
  end
  data.fsample = cfg.fsample;
  
  %%%%%%% POWER TO POWER CORRELATION %%%%%%%%%
elseif strcmp(cfg.method, 'amplow_amphigh')
  
  % sanity checks
  if cfg.s2.freq < cfg.s1.freq || cfg.s1.freq < cfg.s4.freq
    ft_error('with method is powlow_powhigh freq s4 < s1 < s2')
  end
  if cfg.s2.freq > cfg.fsample/2
    ft_error('you cannot have a frequency higher than the sample frequency/2')
  end
  if cfg.s3.freq ~= 0 || cfg.s3.phase ~= 0
    ft_warning('for method powlow_powhigh s3 is DC and therefore expect freq and phase to be zero but they are not')
  end
  if cfg.s3.ampl < cfg.s4.ampl
    ft_warning('expect amplitude s3 (=DC) not to be smaller than amplitude s4 (= AM frequency)')
  end
  
  % make data
  for iTr = 1:length(timevec)
    
    if ischar(cfg.s1.phase); phase_low  = rand * 2 * pi; else phase_low = cfg.s1.phase;    end
    if ischar(cfg.s2.phase); phase_high = rand * 2 * pi; else phase_high = cfg.s2.phase;   end
    if ischar(cfg.s3.phase); phase_DC   = rand * 2 * pi; else phase_DC = cfg.s3.phase;     end
    if ischar(cfg.s4.phase); phase_AM   = rand * 2 * pi; else phase_AM = cfg.s4.phase; end
    high     = cfg.s2.ampl*cos(2*pi*cfg.s2.freq*timevec{iTr} + phase_high);
    low      = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_low);
    AM       = cfg.s4.ampl*cos(2*pi*cfg.s4.freq*timevec{iTr} + phase_AM);
    DC       = cfg.s3.ampl*cos(2*pi*0*timevec{iTr} + phase_DC);
    noise    = cfg.noise.ampl*randn(size(timevec{iTr}));
    lowmix  = ((AM + DC) .* low);
    highmix = ((AM + DC) .* high);
    mix     = lowmix + highmix + noise;
    
    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output, 'all')
      data.trial{iTr}(2,:) = low;
      data.trial{iTr}(3,:) = high;
      data.trial{iTr}(4,:) = DC;
      data.trial{iTr}(5,:) = noise;
      data.trial{iTr}(6,:) = AM;
      data.trial{iTr}(7,:) = lowmix;
      data.trial{iTr}(8,:) = highmix;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr
  
  data.label{1} = 'mix';
  if strcmp(cfg.output, 'all')
    data.label{2} = 's1 (low)';
    data.label{3} = 's2 (high)';
    data.label{4} = 's3 (DC)';
    data.label{5} = 'noise';
    data.label{6} = 's4 (AM)';
    data.label{7} = 'mixlow';
    data.label{8} = 'mixhigh';
  end
  data.fsample = cfg.fsample;
  
  %%%%%%% PHASE TO FREQUENCY CORRELATION %%%%%%%%%
elseif strcmp(cfg.method, 'phalow_freqhigh')
  
  % sanity checks
  if cfg.s1.freq > cfg.fsample/2 || cfg.s2.freq > cfg.fsample/2
    ft_error('you cannot have a frequency higher than the sample frequency/2')
  end
  
  % make data
  for iTr = 1:length(timevec)
    
    if ischar(cfg.s1.phase); phase_s1 = rand * 2 * pi; else phase_s1 = cfg.s1.phase;    end
    if ischar(cfg.s2.phase); phase_s2 = rand * 2 * pi; else phase_s2= cfg.s2.phase;    end
    s1            = cfg.s1.ampl .* cos(2*pi*cfg.s1.freq * timevec{iTr} + phase_s1); % to be modulated signal
    s2            = cfg.s2.ampl .* cos(2*pi*cfg.s2.freq * timevec{iTr} + phase_s2); % modulation of instantaneous phase
    inst_pha_base = 2*pi*cfg.s1.freq * timevec{iTr} + phase_s1; % unmodulated instantaneous phase s1 (linear)
    inst_pha_mod  = s2;                                    % modulation of instantaneous phase
    inst_pha      = inst_pha_base + inst_pha_mod;
    noise         = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix           = cfg.s1.ampl .* cos(inst_pha) + noise;
    
    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output, 'all')
      data.trial{iTr}(2,:) = s1;
      data.trial{iTr}(3,:) = s2;
      data.trial{iTr}(4,:) = noise;
      data.trial{iTr}(5,:) = inst_pha_base;
      data.trial{iTr}(6,:) = inst_pha_mod;
      data.trial{iTr}(7,:) = inst_pha;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr
  
  data.label{1} = 'mix';
  if strcmp(cfg.output, 'all')
    data.label{2} = 's1';
    data.label{3} = 's2';
    data.label{4} = 'noise';
    data.label{5} = 'inst phase base';
    data.label{6} = 'inst phase modulation (=s2)';
    data.label{7} = 'inst phase';
  end
  data.fsample = cfg.fsample;
  
  %%%%%%% ASYMETRIC POSITIVE AND NEGATIVE PEAKS %%%%%%%%%
elseif strcmp(cfg.method, 'asymmetric')
  
  % make data
  for iTr = 1:length(timevec)
    if ischar(cfg.s1.phase); phase_s1 = rand * 2 *pi; else phase_s1 = cfg.s1.phase; end
    
    s1    = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_s1);
    tmp   = cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_s1);  % same signal but with unit amplitude
    tmp   = (tmp+1)/2;                                 % scaled and shifted between 0 and 1
    tmp   = tmp.^(cfg.asymmetry+1);                    % made asymmetric
    tmp   = (tmp - mean(tmp))*2*cfg.s1.ampl;           % rescale
    s2    = tmp;
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix   = s2 + noise;
    
    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output, 'all')
      data.trial{iTr}(2,:) = s1;
      data.trial{iTr}(3,:) = s2;
      data.trial{iTr}(4,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr
  
  data.label{1} = 'mix';
  if strcmp(cfg.output, 'all')
    data.label{2} = 's1';
    data.label{3} = 's2';
    data.label{4} = 'noise';
  end
  data.fsample = cfg.fsample;
  
else
  ft_error('unknown method specified')
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble randomseed
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
