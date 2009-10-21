function [data] = freqsimulation(cfg)

% FREQSIMULATION makes simulated data in FieldTrip format. The data is built
% up from fifferent frequencies and can contain a signal in which the
% different frequencies interact (i.e. cross-frequency coherent). Different
% methods are possible to make data with special properties.
%
% Use as
%   [data] = freqsimulation(cfg)
%
% The configuration options include
%   cfg.method      = The methods are explained in more detail below, but they can be
%                     'superimposed'    simply add the contribution of the different frequencies
%                     'broadband'       create a single broadband signal component
%                     'phalow_amphigh'  phase of low freq correlated with amplitude of high freq
%                     'amplow_amphigh'  amplitude of low freq correlated with amplithude of high freq
%                     'phalow_freqhigh' phase of low freq correlated with frequency of high signal
%                     'asymmetric'      single signal component with asymmetric positive/negative deflections
%   cfg.output      = which channels should be in the output data, can be 'mixed' or 'all' (default = 'all')
%
% The number of trials and the time axes of the trials can be specified by
%   cfg.fsample     = simulated sample frequency
%   cfg.trllen      = length of simulated trials in seconds
%   cfg.numtrl      = number of simulated trials
% or by 
%   cfg.time        = cell-array with one time axis per trial (i.e. from another dataset)
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

% Copyright (C) 2007-2008, Ingrid Nieuwenhuis & Robert Oostenveld, F.C. Donders Centre
%
% $Log: freqsimulation.m,v $
% Revision 1.15  2009/10/12 14:15:06  jansch
% some typo fixes in comments
%
% Revision 1.14  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.13  2008/07/10 08:39:36  roboos
% updated documentation
%
% Revision 1.12  2008/06/17 16:13:18  sashae
% now using preproc_modules
%
% Revision 1.11  2008/02/26 14:39:25  roboos
% fixed bug in cfg.time
%
% Revision 1.10  2008/01/31 09:21:54  roboos
% look at length of timevec instead of cfg.numtrl
%
% Revision 1.9  2008/01/24 12:18:56  roboos
% allow user to specify his own time axes
%
% Revision 1.8  2008/01/24 11:55:45  roboos
% added method 'asymmetric'
%
% Revision 1.7  2007/11/08 12:57:26  roboos
% renamed method amplow_amphigh_1chan into amplow_amphigh
%
% Revision 1.6  2007/11/07 09:14:28  ingnie
% added method 'phalow_freqhigh'
%
% Revision 1.5  2007/11/06 14:47:53  ingnie
% changed naming from power to amplitude, changed documentation methods phalow_amphigh and amplow_amphigh, changed some variable names to match documentation
%
% Revision 1.4  2007/11/05 11:29:54  roboos
% added cfg.method=broadband
% add version details to output cfg
%
% Revision 1.3  2007/11/05 10:36:32  roboos
% added method superimposed, changed whitespace, moved some code out of for-loops
%
% Revision 1.2  2007/10/23 15:34:42  ingnie
% added option powlow_powhigh_1chan
%
% Revision 1.1  2007/08/08 06:33:10  roboos
% renamed simluatedata to freqsimuation
%
% Revision 1.1  2007/08/07 10:45:33  ingnie
% first implementation
%

fieldtripdefs

% set defaults
if ~isfield(cfg, 'method'),        cfg.method = 'phalow_amphigh';         end
if ~isfield(cfg, 'output'),        cfg.output = 'all';                    end
if ~isfield(cfg, 'time'),          cfg.time = [];                         end

if isempty(cfg.time)
  if ~isfield(cfg, 'fsample'),       cfg.fsample = 1200;                    end
  if ~isfield(cfg, 'trllen'),        cfg.trllen = 1;                        end
  if ~isfield(cfg, 'numtrl'),        cfg.numtrl = 1;                        end
else
  cfg.trllen = [];                                    % can be variable
  cfg.fsample = 1/(cfg.time{1}(2) - cfg.time{1}(1));  % determine from time-axis
  cfg.numtrl = length(cfg.time);
end

if strcmp(cfg.method,'superimposed')
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
  if ~isfield(cfg.s3, 'ampl'),    	cfg.s3.ampl  = 0;                     end
end

if strcmp(cfg.method,'broadband')
  if ~isfield(cfg, 'n1'),           cfg.n1 = [];                          end
  if ~isfield(cfg.n1, 'ampl'),      cfg.n1.ampl  = 1;                     end
  if ~isfield(cfg.n1, 'bpfreq'),    cfg.n1.bpfreq  = [30 50];             end
  if ~isfield(cfg, 'n2'),           cfg.n2 = [];                          end
  if ~isfield(cfg.n2, 'ampl'),      cfg.n2.ampl  = 1;                     end
  if ~isfield(cfg.n2, 'bpfreq'),    cfg.n2.bpfreq  = [80 120];            end
end

if strcmp(cfg.method,'phalow_amphigh')
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
  if ~isfield(cfg.s3, 'ampl'),    	cfg.s3.ampl = cfg.s1.ampl;            end
end

if strcmp(cfg.method,'amplow_amphigh')
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
  if ~isfield(cfg.s3, 'ampl'),    	cfg.s3.ampl = cfg.s4.ampl;            end
end

if strcmp(cfg.method,'phalow_freqhigh')
  if ~isfield(cfg, 's1'),           cfg.s1 = [];                          end
  if ~isfield(cfg.s1, 'freq'),      cfg.s1.freq = 20;                     end
  if ~isfield(cfg.s1, 'phase'),     cfg.s1.phase = 0;                     end
  if ~isfield(cfg.s1, 'ampl'),      cfg.s1.ampl = 1;                      end
  if ~isfield(cfg, 's2'),           cfg.s2 = [];                          end
  if ~isfield(cfg.s2, 'freq'),      cfg.s2.freq = 2;                      end
  if ~isfield(cfg.s2, 'phase'),     cfg.s2.phase = -0.5 * pi;             end %then base freq at t=0
  if ~isfield(cfg.s2, 'ampl'),      cfg.s2.ampl = pi;                     end
end

if strcmp(cfg.method,'asymmetric')
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
  Nsamp_tr = cfg.fsample * cfg.trllen;
  for iTr = 1 : cfg.numtrl
    timevec{iTr} = (1:Nsamp_tr)/cfg.fsample;
  end
end


%%%%%%% SUPERIMPOSED, SIMPLY ADD THE SIGNALS %%%%%%%%%
if strcmp(cfg.method,'superimposed')

  % make data
  for iTr = 1 : length(timevec)
    if isstr(cfg.s1.phase); phase_s1 = rand * 2 *pi; else phase_s1 = cfg.s1.phase; end
    if isstr(cfg.s2.phase); phase_s2 = rand * 2 *pi; else phase_s2 = cfg.s2.phase; end
    if isstr(cfg.s3.phase); phase_s3 = rand * 2 *pi; else phase_s3 = cfg.s3.phase; end

    s1    = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_s1);
    s2    = cfg.s2.ampl*cos(2*pi*cfg.s2.freq*timevec{iTr} + phase_s2);
    s3    = cfg.s3.ampl*cos(2*pi*cfg.s3.freq*timevec{iTr} + phase_s3);
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix   = s1 + s2 + s3 + noise;

    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output,'all')
      data.trial{iTr}(2,:) = s1;
      data.trial{iTr}(3,:) = s2;
      data.trial{iTr}(4,:) = s3;
      data.trial{iTr}(5,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr

  data.label{1} = 'mix';
  if strcmp(cfg.output,'all')
    data.label{2} = 's1';
    data.label{3} = 's2';
    data.label{4} = 's3';
    data.label{5} = 'noise';
  end
  data.fsample = cfg.fsample;

  %%%%%%% SUPERIMPOSED BROADBAND SIGNAL %%%%%%%%%
elseif strcmp(cfg.method,'broadband')

  % make data
  for iTr = 1 : length(timevec)
    n1    = preproc_bandpassfilter(cfg.n1.ampl*randn(size(timevec{iTr})), cfg.fsample, cfg.n1.bpfreq);
    n2    = preproc_bandpassfilter(cfg.n2.ampl*randn(size(timevec{iTr})), cfg.fsample, cfg.n2.bpfreq);
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix   = n1 + n2 + noise;

    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output,'all')
      data.trial{iTr}(2,:) = n1;
      data.trial{iTr}(3,:) = n2;
      data.trial{iTr}(4,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr

  data.label{1} = 'mix';
  if strcmp(cfg.output,'all')
    data.label{2} = 'n1';
    data.label{3} = 'n2';
    data.label{4} = 'noise';
  end
  data.fsample = cfg.fsample;

  %%%%%%% PHASE TO AMPLITUDE CORRELATION %%%%%%%%%
elseif strcmp(cfg.method,'phalow_amphigh')

  % sanity checks
  if cfg.s2.freq < cfg.s1.freq
    error('with method is phalow_amphigh freq s2 should be higher than freq s1')
  end
  if cfg.s2.freq > cfg.fsample/2
    error('you cannot have a frequency higher than the sample frequency/2')
  end
  if cfg.s3.freq ~= 0 || cfg.s3.phase ~= 0
    warning('for method phalow_amphigh s3 is DC and therefore expect freq and phase to be zero but they are not')
  end
  if cfg.s3.ampl < cfg.s1.ampl
    warning('expect amplitude s3 (=DC) not to be smaller than amplitude s1 (=low frequency)')
  end

  % make data
  for iTr = 1 : length(timevec)

    if isstr(cfg.s1.phase); phase_AM = rand * 2 *pi;  else phase_AM = cfg.s1.phase;  end
    if isstr(cfg.s2.phase); phase_high = rand * 2 *pi; else phase_high = cfg.s2.phase; end
    if isstr(cfg.s3.phase); phase_DC = rand * 2 *pi;   else phase_DC = cfg.s3.phase;   end
    high  = cfg.s2.ampl*cos(2*pi*cfg.s2.freq*timevec{iTr} + phase_high);
    AM   = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_AM);
    DC    = cfg.s3.ampl*cos(2*pi*0*timevec{iTr} + phase_DC);
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix = ((AM + DC) .* high) + noise;

    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output,'all')
      data.trial{iTr}(2,:) = AM;
      data.trial{iTr}(3,:) = high;
      data.trial{iTr}(4,:) = DC;
      data.trial{iTr}(5,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr

  data.label{1} = 'mix';
  if strcmp(cfg.output,'all')
    data.label{2} = 's1 (AM)';
    data.label{3} = 's2 (high)';
    data.label{4} = 's3 (DC)';
    data.label{5} = 'noise';
  end
  data.fsample = cfg.fsample;

  %%%%%%% POWER TO POWER CORRELATION %%%%%%%%%
elseif strcmp(cfg.method,'amplow_amphigh')

  % sanity checks
  if cfg.s2.freq < cfg.s1.freq || cfg.s1.freq < cfg.s4.freq
    error('with method is powlow_powhigh freq s4 < s1 < s2')
  end
  if cfg.s2.freq > cfg.fsample/2
    error('you cannot have a frequency higher than the sample frequency/2')
  end
  if cfg.s3.freq ~= 0 || cfg.s3.phase ~= 0
    warning('for method powlow_powhigh s3 is DC and therefore expect freq and phase to be zero but they are not')
  end
  if cfg.s3.ampl < cfg.s4.ampl
    warning('expect amplitude s3 (=DC) not to be smaller than amplitude s4 (= AM frequency)')
  end

  % make data
  for iTr = 1 : length(timevec)

    if isstr(cfg.s1.phase); phase_low = rand * 2 *pi;    else phase_low = cfg.s1.phase;    end
    if isstr(cfg.s2.phase); phase_high = rand * 2 *pi;   else phase_high = cfg.s2.phase;   end
    if isstr(cfg.s3.phase); phase_DC = rand * 2 *pi;     else phase_DC = cfg.s3.phase;     end
    if isstr(cfg.s4.phase); phase_AM = rand * 2 *pi;     else phase_AM = cfg.s4.phase; end
    high     = cfg.s2.ampl*cos(2*pi*cfg.s2.freq*timevec{iTr} + phase_high);
    low      = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_low);
    AM       = cfg.s4.ampl*cos(2*pi*cfg.s4.freq*timevec{iTr} + phase_AM);
    DC       = cfg.s3.ampl*cos(2*pi*0*timevec{iTr} + phase_DC);
    noise    = cfg.noise.ampl*randn(size(timevec{iTr}));
    lowmix  = ((AM + DC) .* low);
    highmix = ((AM + DC) .* high);
    mix     = lowmix + highmix + noise;

    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output,'all')
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
  if strcmp(cfg.output,'all')
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
elseif strcmp(cfg.method,'phalow_freqhigh')

  % sanity checks
  if cfg.s1.freq > cfg.fsample/2 || cfg.s2.freq > cfg.fsample/2
    error('you cannot have a frequency higher than the sample frequency/2')
  end

  % make data
  for iTr = 1 : length(timevec)

    if isstr(cfg.s1.phase); phase_s1  = rand * 2 *pi;    else phase_s1 = cfg.s1.phase;    end
    if isstr(cfg.s2.phase); phase_s2 = rand * 2 *pi;     else phase_s2= cfg.s2.phase;    end
    s1            = cfg.s1.ampl .* cos(2*pi*cfg.s1.freq * timevec{iTr} + phase_s1); % to be modulated signal
    s2            = cfg.s2.ampl .* cos(2*pi*cfg.s2.freq * timevec{iTr} + phase_s2); % modulation of instantaneous phase
    inst_pha_base = 2*pi*cfg.s1.freq * timevec{iTr} + phase_s1; % unmodulated instantaneous phase s1 (linear)
    inst_pha_mod  = s2;                                    % modulation of instantaneous phase
    inst_pha      = inst_pha_base + inst_pha_mod;
    noise         = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix           = cfg.s1.ampl .* cos(inst_pha) + noise;

    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output,'all')
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
  if strcmp(cfg.output,'all')
    data.label{2} = 's1';
    data.label{3} = 's2';
    data.label{4} = 'noise';
    data.label{5} = 'inst phase base';
    data.label{6} = 'inst phase modulation (=s2)';
    data.label{7} = 'inst phase';
  end
  data.fsample = cfg.fsample;

  %%%%%%% ASYMETRIC POSITIVE AND NEGATIVE PEAKS %%%%%%%%%
elseif strcmp(cfg.method,'asymmetric')

  % make data
  for iTr = 1 : length(timevec)
    if isstr(cfg.s1.phase); phase_s1 = rand * 2 *pi; else phase_s1 = cfg.s1.phase; end

    s1    = cfg.s1.ampl*cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_s1);
    tmp   = cos(2*pi*cfg.s1.freq*timevec{iTr} + phase_s1);  % same signal but with unit amplitude
    tmp   = (tmp+1)/2;                                 % scaled and shifted between 0 and 1
    tmp   = tmp.^(cfg.asymmetry+1);                    % made asymmetric
    tmp   = (tmp - mean(tmp))*2*cfg.s1.ampl;           % rescale
    s2    = tmp;
    noise = cfg.noise.ampl*randn(size(timevec{iTr}));
    mix   = s2 + noise;

    data.trial{iTr}(1,:) = mix;
    if strcmp(cfg.output,'all')
      data.trial{iTr}(2,:) = s1;
      data.trial{iTr}(3,:) = s2;
      data.trial{iTr}(4,:) = noise;
    end
    data.time{iTr} = timevec{iTr};
  end % for iTr

  data.label{1} = 'mix';
  if strcmp(cfg.output,'all')
    data.label{2} = 's1';
    data.label{3} = 's2';
    data.label{4} = 'noise';
  end
  data.fsample = cfg.fsample;

else
  error('unknown method specified')
end

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: freqsimulation.m,v 1.15 2009/10/12 14:15:06 jansch Exp $';

% remember the exact configuration details in the output
data.cfg = cfg;
