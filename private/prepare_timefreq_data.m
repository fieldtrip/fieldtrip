function [cfg, data] = prepare_timefreq_data(cfg, varargin);

% PREPARE_TIMEFREQ_DATA collects the overlapping data from multiple ERPs
% or ERFs according to the channel/time/freq-selection in the configuration
% and returns it with a dimord of 'repl_chan_freq_time'.
%
% Supported input data is from TIMELOCKANALYSIS, TIMELOCKGRANDAVERAGE,
% FREQANALYSIS or FREQDESCRIPTIVES.
%
% Use as
%   [cfg, data] = prepare_timefreq_data(cfg, varargin)
% where the configuration can contain
%   cfg.channel            = cell-array of size Nx1, see CHANNELSELECTION, or
%   cfg.channelcmb         = cell-array of size Nx2, see CHANNELCOMBINATION
%   cfg.latency            = [begin end], can be 'all' (default)
%   cfg.frequency          = [begin end], can be 'all' (default)
%   cfg.avgoverchan        = 'yes' or 'no' (default)
%   cfg.avgovertime        = 'yes' or 'no' (default)
%   cfg.avgoverfreq        = 'yes' or 'no' (default)
%   cfg.datarepresentation = 'cell-array' or 'concatenated'
%   cfg.precision          = 'single' or 'double' (default)
%
% If the data contains both power- and cross-spectra, they will be concatenated
% and treated together. In that case you should specify cfg.channelcmb instead
% of cfg.channel.
%
% This function can serve as a subfunction for TIMEFREQRANDOMIZATION,
% TIMEFREQCLUSTER and other statistical functions, so that those functions
% do not have to do the bookkeeping themselves.

% KNOWN BUGS AND LIMITATIONS
%   copying of data is not memory efficient for cell-array output
%   cfg.precision is not used for cell-array output

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: prepare_timefreq_data.m,v $
% Revision 1.28  2009/01/26 12:44:38  marvger
% added control condition when computing output.dof; computed
% incorrectly for single trial data. probably needs a more elegant solution.
%
% Revision 1.27  2006/09/05 10:42:25  roboos
% added fixdimord for backward compatibility with old data (toi->time etc)
%
% Revision 1.26  2006/05/08 09:02:02  roboos
% fixed bug in assesment of the dimensionality and size of the data, relevant vor averaging over time or frequencies
%
% Revision 1.25  2006/03/21 15:03:43  roboos
% added support for fourier data, not yet keeping track of tapers and trials
%
% Revision 1.24  2006/02/28 12:18:16  erimar
% Added handling of degrees of freedom, which may be passed by freqanalysis.
%
% Revision 1.23  2006/02/23 10:28:17  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.22  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.21  2005/04/05 06:54:08  roboos
% fixed bug in data selectino which occurred when the order of the channels was different in the different input datasets
%
% Revision 1.20  2005/03/23 08:01:27  roboos
% fixed typing mistake, causing invalid assignment of variable size
%
% Revision 1.19  2005/03/23 07:25:18  roboos
% extended support for single precision to ERPs and power
% fixed bug with singleton dimensions for frequency and time-selection
%
% Revision 1.18  2005/01/19 12:39:56  erimar
% Added pack for efficient memory use.
%
% Revision 1.17  2004/11/22 08:45:25  roboos
% fixed bug in data.label for avgoverchan (was string instead of cell-array)
%
% Revision 1.16  2004/11/19 09:12:40  roboos
% fixed bug in sprintf with cell-array
%
% Revision 1.15  2004/11/10 12:57:54  erimar
% correct wrong Nfreq and Ntime values (made them dataset-dependent)
%
% Revision 1.14  2004/11/09 15:53:17  roboos
% fixed bug in the separate handling of chansel for cross- and auto-spectra
%
% Revision 1.13  2004/11/09 13:59:17  roboos
% fixed bug in checking whether all chan/freq/time bins were selected
%
% Revision 1.12  2004/11/08 13:06:53  roboos
% made changes to improve memory efficiency, mainly for concatenated output and cross-spectra
% implemented single-precision output datatype (cfg.precision)
% split handling of auto- and cross-spectra in forcedimord
% made avgoverdim optional, do not call if avgdim is empty
% made dataselection optional, do not explicitely select if the selection contains everything
%
% Revision 1.11  2004/11/02 11:33:19  roboos
% fixed double concatenation of power and cross-spectra
% improved output of time/freq axis if time/freq is not present in data
%
% Revision 1.10  2004/11/01 15:55:43  roboos
% fixed hascrsspctrm in forcedimord for ERP input data
%
% Revision 1.9  2004/11/01 13:51:03  roboos
% fixed bug in labelcmb v.s. channelcmb
% fixed bug in default channelcmb
% print either the string "channels" or "channelcombinations" in fprintf feedback
% cleaned up internal handling of label and labelcombinatinos using two subfunctions that convert between the two representations
%
% Revision 1.8  2004/10/29 11:34:16  roboos
% cleaned up the code on many places
% renamed subfunction in forcedimord
% made permute and reshape in forcedimord conditional
%
% Revision 1.7  2004/10/29 09:45:54  roboos
% complete redesign with numerous changes, with the most important ones being:
% added support for data in matlab files
% attempt to make memory efficient
% output dimord is fixed repl_chan_freq_time
% added output design matrix (in case of concatenated)
%
% Revision 1.6  2004/10/22 16:26:04  roboos
% added default for cfg.channelcmb
% added backward compatibility support for ERPs without dimord
% fixed some small bugs
%
% Revision 1.5  2004/10/22 10:48:28  roboos
% implemented support for cross-spectra and selec tion of channel combinations
%
% Revision 1.4  2004/10/13 14:11:30  roboos
% changed cfg.previous, now consistent over functions with try-statement
% and will also work if the input already had a cfg.previous
%
% Revision 1.3  2004/10/13 11:24:03  roboos
% added support for cfg.datarepresentation = 'cell-array' or 'concatenated'
%
% Revision 1.2  2004/10/12 13:51:19  roboos
% implemented avgoverchan/freq/time
% many general improvements
%
% Revision 1.1  2004/10/06 15:12:13  roboos
% initial implementation to serve as subfunction for main statistical functions
%

% set the defaults
if ~isfield(cfg, 'channel'),              cfg.channel = 'all';                     end
if ~isfield(cfg, 'channelcmb'),           cfg.channelcmb = [];                     end
if ~isfield(cfg, 'latency'),              cfg.latency = 'all';                     end
if ~isfield(cfg, 'frequency'),            cfg.frequency = 'all';                   end
if ~isfield(cfg, 'avgoverchan'),          cfg.avgoverchan = 'no';                  end
if ~isfield(cfg, 'avgovertime'),          cfg.avgovertime = 'no';                  end
if ~isfield(cfg, 'avgoverfreq'),          cfg.avgoverfreq = 'no';                  end
if ~isfield(cfg, 'datarepresentation'),   cfg.datarepresentation = 'cell-array';   end
if ~isfield(cfg, 'precision'),            cfg.precision = 'double';                end

% initialize
containsdof = false;

% count the number of input datasets
Nvarargin = length(varargin);
remember = {};

% initialize
swapmemfile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect information over all data sets required for a consistent selection
% this is the first time that the data will be swapped from file into memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c=1:Nvarargin
  % swap the data into memory
  varargin{c} = swapmemfile(varargin{c});
  % for backward compatibility with old data structures
  varargin{c} = fixdimord(varargin{c});
  % reformat the data, and combine cross and power spectra
  [remember{c}, hascrsspctrm] = forcedimord(varargin{c});
  % keep a small part of the data in memory
  if hascrsspctrm
    remember{c} = rmfield(remember{c}, 'dat');
    remember{c} = rmfield(remember{c}, 'crs');
  else
    remember{c} = rmfield(remember{c}, 'dat');
  end
  % swap the data out of memory
  varargin{c} = swapmemfile(varargin{c});
end

if hascrsspctrm
  % by default this should be empty
  cfg.channel = [];
  if isempty(cfg.channelcmb)
    % default is to copy the channelcombinations from the first dataset
    % and these have been concatenated by forcedimord -> tear them apart
    cfg.channelcmb = label2cmb(remember{1}.label);
  end
  % instead of matching channelcombinations, we will match channels that
  % consist of the concatenated labels of the channelcombinations
  cfg.channel = cmb2label(cfg.channelcmb);
  cfg.channelcmb = []; % this will be reconstructed at the end
else
  % by default this should be empty
  cfg.channelcmb = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust the selection of channels, latency and frequency so that it
% matches with all datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c=1:Nvarargin
  % match the channel selection with this dataset
  cfg.channel = channelselection(cfg.channel, remember{c}.label);

  % match the selected latency with this dataset
  if isempty(findstr(remember{c}.dimord, 'time')) || all(isnan(remember{c}.time))
    cfg.latency = [];
  elseif ischar(cfg.latency) && strcmp(cfg.latency, 'all')
    cfg.latency = [];
    cfg.latency(1) = min(remember{c}.time);
    cfg.latency(2) = max(remember{c}.time);
  else
    cfg.latency(1) = max(cfg.latency(1), min(remember{c}.time));
    cfg.latency(2) = min(cfg.latency(2), max(remember{c}.time));
  end

  % match the selected frequency with this dataset
  if isempty(findstr(remember{c}.dimord, 'freq')) || all(isnan(remember{c}.freq))
    cfg.frequency = [];
  elseif ischar(cfg.frequency) && strcmp(cfg.frequency, 'all')
    cfg.frequency = [];
    cfg.frequency(1) = min(remember{c}.freq);
    cfg.frequency(2) = max(remember{c}.freq);
  else
    cfg.frequency(1) = max(cfg.frequency(1), min(remember{c}.freq));
    cfg.frequency(2) = min(cfg.frequency(2), max(remember{c}.freq));
  end

  % keep track of the number of replications
  Nrepl(c) = length(remember{c}.repl);
end

% count the number of channels, time bins and frequency bins
dimord = remember{c}.dimord;
if ~isempty(cfg.latency) && ~all(isnan(cfg.latency))
  timesel = nearest(remember{1}.time, cfg.latency(1)):nearest(remember{1}.time, cfg.latency(2));
else
  timesel = 1;
end
if ~isempty(cfg.frequency) && ~all(isnan(cfg.frequency))
  freqsel = nearest(remember{1}.freq, cfg.frequency(1)):nearest(remember{1}.freq, cfg.frequency(2));
else
  freqsel = 1;
end
Nchan = length(cfg.channel);
Ntime = length(timesel);
Nfreq = length(freqsel);

if ~hascrsspctrm
  % print the information for channels
  if Nchan<1
    error('channel selection is empty')
  elseif strcmp(cfg.avgoverchan, 'yes')
    fprintf('averaging over %d channels\n', Nchan);
    Nchan = 1; % after averaging
  else
    fprintf('selected %d channels\n', Nchan);
  end
else
  % print the (same) information for channelcombinations
  if Nchan<1
    error('channelcombination selection is empty')
  elseif strcmp(cfg.avgoverchan, 'yes')
    fprintf('averaging over %d channelcombinations\n', Nchan);
    Nchan = 1; % after averaging
  else
    fprintf('selected %d channelcombinations\n', Nchan);
  end
end

if Ntime<1
  error('latency selection is empty')
elseif strcmp(cfg.avgovertime, 'yes')
  fprintf('averaging over %d time bins\n', Ntime);
  Ntime = 1; % after averaging
else
  fprintf('selected %d time bins\n', Ntime);
end

if Nfreq<1
  error('latency selection is empty')
elseif strcmp(cfg.avgoverfreq, 'yes')
  fprintf('averaging over %d frequency bins\n', Nfreq);
  Nfreq = 1; % after averaging
else
  fprintf('selected %d frequency bins\n', Nfreq);
end

% determine whether, and over which dimensions the data should be averaged
tok = tokenize(dimord, '_');
chandim = find(strcmp(tok, 'chan'));
timedim = find(strcmp(tok, 'time'));
freqdim = find(strcmp(tok, 'freq'));
avgdim = [];
if strcmp(cfg.avgoverchan, 'yes') && ~isempty(chandim)
  avgdim = [avgdim chandim];
end
if strcmp(cfg.avgovertime, 'yes') && ~isempty(timedim)
  avgdim = [avgdim timedim];
end
if strcmp(cfg.avgoverfreq, 'yes') && ~isempty(freqdim)
  avgdim = [avgdim freqdim];
end

pack;
% this will hold the output data
if strcmp(cfg.datarepresentation, 'cell-array')
  dat = {};
  dof = {};
elseif strcmp(cfg.datarepresentation, 'concatenated')
  % preallocate memory to hold all the data
  if hascrsspctrm
    if str2num(version('-release'))>=14
      dat = complex(zeros(sum(Nrepl), Nchan, Nfreq, Ntime, cfg.precision));
    else
      % use double precision for default and Matlab versions prior to release 14
      if ~strcmp(cfg.precision, 'double')
        warning('Matlab versions prior to release 14 only support double precision');
      end
      dat = complex(zeros(sum(Nrepl), Nchan, Nfreq, Ntime));
    end
  else
    if str2num(version('-release'))>=14
      dat = zeros(sum(Nrepl), Nchan, Nfreq, Ntime, cfg.precision);
    else
      % use double precision for default and Matlab versions prior to release 14
      if ~strcmp(cfg.precision, 'double')
        warning('Matlab versions prior to release 14 only support double precision');
      end
      dat = zeros(sum(Nrepl), Nchan, Nfreq, Ntime);
    end
  end
  if isempty(avgdim)
    dof = zeros(sum(Nrepl), Nfreq, Ntime);
  end;
else
  error('unsupported output data representation');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the desired channels, latency and frequency bins from each dataset
% this is the second time that the data will be swapped from file into memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c=1:Nvarargin
  % swap the data into memory
  varargin{c} = swapmemfile(varargin{c});
  % reformat the data, and combine cross and power spectra
  [varargin{c}, hascrsspctrm] = forcedimord(varargin{c});

  % select the channels, sorted according to the order in the configuration
  [dum, chansel] = match_str(cfg.channel, varargin{c}.label);
  if ~isempty(cfg.latency)
    timesel = nearest(varargin{c}.time, cfg.latency(1)):nearest(varargin{c}.time, cfg.latency(2));
  end
  if ~isempty(cfg.frequency)
    freqsel = nearest(varargin{c}.freq, cfg.frequency(1)):nearest(varargin{c}.freq, cfg.frequency(2));
  end

  if hascrsspctrm
    Nfreq = size(varargin{c}.crs, 3);
    Ntime = size(varargin{c}.crs, 4);
    % separate selection for the auto-spectra and cross-spectra
    Ndat = size(varargin{c}.dat, 2);                                        % count the total number of auto-spectra
    Ncrs = size(varargin{c}.crs, 2);                                        % count the total number of cross-spectra
    [dum, chansel_dat] = ismember(chansel, 1:Ndat);                         % selection of auto-spectra
    chansel_dat=chansel_dat(chansel_dat~=0)';
    [dum, chansel_crs] = ismember(chansel, (1:Ncrs)+Ndat);                  % selection of cross-spectra
    chansel_crs=chansel_crs(chansel_crs~=0)';
    Ndatsel = length(chansel_dat);                                          % count the number of selected auto-spectra
    Ncrssel = length(chansel_crs);                                          % count the number of selected cross-spectra
  else
    Nfreq = size(varargin{c}.dat, 3);
    Ntime = size(varargin{c}.dat, 4);
  end

  % collect the data, average if required
  if strcmp(cfg.datarepresentation, 'cell-array')
    if hascrsspctrm
      dat{c} = cat(2, avgoverdim(varargin{c}.dat(:, chansel_dat, freqsel, timesel), avgdim), avgoverdim(varargin{c}.crs(:, chansel_crs, freqsel, timesel), avgdim));
    else
      dat{c} = avgoverdim(varargin{c}.dat(:, chansel, freqsel, timesel), avgdim);
    end
    if isfield(varargin{c},'dof')
      containsdof = isempty(avgdim);
      if containsdof
        dof{c} = varargin{c}.dof(:,freqsel,timesel);
      else
         warning('Degrees of freedom cannot be calculated if averaging over channels, frequencies, or time bins is requested.'); 
      end;
    end;

  elseif strcmp(cfg.datarepresentation, 'concatenated')
    if c==1
      trialsel = 1:Nrepl(1);
    else
      trialsel = (1+sum(Nrepl(1:(c-1)))):sum(Nrepl(1:c));
    end

    % THIS IS IMPLEMENTED WITH STRICT MEMORY CONSIDERATIONS
    % handle cross-spectrum and auto-spectrum separately
    % first copy complex part (crs) then real part (pow), otherwise dat first will be made real and then complex again
    % copy selection of "all" as complete matrix and not as selection
    % copy selection with for loops (not yet implemented below)

    if isempty(avgdim) && hascrsspctrm
      if length(chansel_crs)==Ncrs && length(freqsel)==Nfreq && length(timesel)==Ntime && all(chansel_crs==(1:Ncrs)) && all(freqsel==(1:Nfreq)) && all(timesel==(1:Ntime))
        % copy everything into the output data array
        % this is more memory efficient than selecting everything
        dat(trialsel,(1:Ncrssel)+Ndatsel,:,:) = varargin{c}.crs;
      else
        % copy only a selection into the output data array
        % this could be made more memory efficient using for-loops
        dat(trialsel,(1:Ncrssel)+Ndatsel,:,:) = varargin{c}.crs(:, chansel_crs, freqsel, timesel);
      end
      dat(trialsel,1:Ndatsel,:,:) = varargin{c}.dat(:, chansel_dat, freqsel, timesel);

    elseif ~isempty(avgdim) && hascrsspctrm
      if length(chansel_crs)==Ncrs && length(freqsel)==Nfreq && length(timesel)==Ntime && all(chansel_crs==(1:Ncrs)) && all(freqsel==(1:Nfreq)) && all(timesel==(1:Ntime))
        % copy everything into the output data array
        % this is more memory efficient than selecting everything
        dat(trialsel,(1:Ncrssel)+Ndatsel,:,:) = avgoverdim(varargin{c}.crs, avgdim);
      else
        % copy only a selection into the output data array
        % this could be made more memory efficient using for-loops
        dat(trialsel,(1:Ncrssel)+Ndatsel,:,:) = avgoverdim(varargin{c}.crs(:, chansel_crs, freqsel, timesel), avgdim);
      end
      dat(trialsel,1:Ndatsel,:,:) = avgoverdim(varargin{c}.dat(:, chansel_dat, freqsel, timesel), avgdim);

    elseif isempty(avgdim) && ~hascrsspctrm
      dat(trialsel,:,:,:) = varargin{c}.dat(:, chansel, freqsel, timesel);

    elseif ~isempty(avgdim) && ~hascrsspctrm
      dat(trialsel,:,:,:) = avgoverdim(varargin{c}.dat(:, chansel, freqsel, timesel), avgdim);
    end

    if isfield(varargin{c},'dof')
      containsdof = isempty(avgdim);
      if containsdof
        dof(trialsel,:,:) = varargin{c}.dof(:,freqsel,timesel);
      else
         warning('Degrees of freedom cannot be calculated if averaging over channels, frequencies, or time bins is requested.'); 
      end;
    end;

  end
  % swap the data out of memory
  varargin{c} = swapmemfile(varargin{c});
end

% collect the output data
data.biol   = dat;
data.dimord = dimord;
if containsdof
    data.dof = dof;
end;

if strcmp(cfg.datarepresentation, 'concatenated')
  % all trials are combined together, construct a design matrix to ensure
  % that they can be identified afterwards
  data.design(:,1) = 1:sum(Nrepl);
  for c=1:Nvarargin
    if c==1
      trialsel = 1:Nrepl(1);
    else
      trialsel = (1+sum(Nrepl(1:(c-1)))):sum(Nrepl(1:c));
    end
    data.design(trialsel,2) = 1:length(trialsel);
    data.design(trialsel,3) = c;
  end
end

if ~hascrsspctrm
  if strcmp(cfg.avgoverchan, 'yes')
    concatlabel = sprintf('%s+', remember{end}.label{chansel});   % concatenate the labels with a '+' in between
    concatlabel(end) = [];                                        % remove the last '+'
    data.label = {concatlabel};                                   % this should be a cell-array
  else
    data.label  = remember{end}.label(chansel);
  end
else
  % specify the correct dimord and channelcombinations
  labelcmb = label2cmb(remember{end}.label(chansel));
  if strcmp(cfg.avgoverchan, 'yes')
    % the cross-spectra have been averaged, therefore there is only one combination left
    str1 = sprintf('%s+', labelcmb{:,1}); str1(end) = [];
    str2 = sprintf('%s+', labelcmb{:,2}); str2(end) = [];
    data.labelcmb = {str1, str2};
  else
    data.labelcmb = labelcmb;
  end
  % the dimord should specify chancmb instead of chan
  s = strfind(dimord, 'chan');
  data.dimord = [data.dimord(1:(s-1)) 'chancmb' data.dimord((s+4):end)];
end

if ~all(isnan(cfg.latency))
  if strcmp(cfg.avgovertime, 'yes')
    data.time = mean(remember{end}.time(timesel));
  else
    data.time = remember{end}.time(timesel);
  end
else
  % the output data contains a time dimension, but the input data did not contain a time axis
  data.time = nan;
  cfg.latency = [];
end

if ~all(isnan(cfg.frequency))
  if strcmp(cfg.avgoverfreq, 'yes')
    data.freq = mean(remember{end}.freq(freqsel));
  else
    data.freq = remember{end}.freq(freqsel);
  end
else
  % the output data contains a freq dimension, but the input data did not contain a freq axis
  data.freq = nan;
  cfg.frequency = [];
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i).name;
end
cfg.version.id = '$Id: prepare_timefreq_data.m,v 1.28 2009/01/26 12:44:38 marvger Exp $';

cfg.previous = [];
% remember the configuration details from the input data
for c=1:Nvarargin
  try, cfg.previous{c} = remember{c}.cfg; end
end

% remember the full configuration details
data.cfg    = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to help with the averaging over multiple dimensions at the
% same time. Furthermore, the output data will have the same NUMBER of
% dimensions as the input data (i.e. no squeezing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dat] = avgoverdim(dat, dim);
siz = size(dat);
siz((end+1):4) = 1;           % this should work at least up to 4 dimensions, also if not present in  the data
for i=dim
  siz(i) = 1;
  dat = mean(dat, i);         % compute the mean over this dimension
  dat = reshape(dat, siz);    % ensure that the number of dimenstions remains the same
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that reformats the data into a consistent struct with a fixed
% dimord and that combines power and cross-spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output, hascrsspctrm] = forcedimord(input);

% for backward compatibility with old timelockanalysis data that does not contain a dimord
if ~isfield(input, 'dimord') && isfield(input, 'avg') && isfield(input, 'var') && isfield(input, 'trial')
  input.dimord = 'repl_chan_time';
elseif ~isfield(input, 'dimord') && isfield(input, 'avg') && isfield(input, 'var') && ~isfield(input, 'trial')
  input.dimord = 'chan_time';
end

outputdim = tokenize('repl_chan_freq_time', '_');
inputdim  = tokenize(input.dimord, '_');
% ensure consistent names of the dimensions
try, inputdim{strcmp(inputdim, 'rpt')}    = 'repl'; end;
try, inputdim{strcmp(inputdim, 'rpttap')} = 'repl'; end;
try, inputdim{strcmp(inputdim, 'trial')}  = 'repl'; end;
try, inputdim{strcmp(inputdim, 'subj')}   = 'repl'; end;
try, inputdim{strcmp(inputdim, 'sgncmb')} = 'chan'; end;
try, inputdim{strcmp(inputdim, 'sgn')}    = 'chan'; end;
try, inputdim{strcmp(inputdim, 'frq')}    = 'freq'; end;
try, inputdim{strcmp(inputdim, 'tim')}    = 'time'; end;

try, output.cfg  = input.cfg;        end; % general
try, output.time = input.time;       end; % for ERPs and TFRs
try, output.freq = input.freq;       end; % for frequency data

if isfield(input, 'powspctrm') && isfield(input, 'crsspctrm')
  hascrsspctrm = 1;
  chandim = find(strcmp(inputdim, 'chan'));
  % concatenate powspctrm and crsspctrm
  % output.dat = cat(chandim, input.powspctrm, input.crsspctrm);
  % seperate output for auto and cross-spectra
  output.dat = input.powspctrm;
  output.crs = input.crsspctrm;
  % concatenate the labels of the channelcombinations into one Nx1 cell-array
  output.label = cmb2label([[input.label(:) input.label(:)]; input.labelcmb]);
elseif isfield(input, 'powspctrm')
  hascrsspctrm = 0;
  % output only power spectrum
  output.dat    = input.powspctrm;
  output.label  = input.label;
elseif isfield(input, 'fourierspctrm')
  hascrsspctrm = 0;
  % output only fourier spectrum
  output.dat    = input.fourierspctrm;
  output.label  = input.label;
else
  hascrsspctrm = 0;
  try, output.dat  = input.avg;        end; % for average ERP
  try, output.dat  = input.trial;      end; % for average ERP with keeptrial
  try, output.dat  = input.individual; end; % for grandaverage ERP with keepsubject
  output.label  = input.label;
end

% determine the size of each dimension
repldim = find(strcmp(inputdim, 'repl'));
chandim = find(strcmp(inputdim, 'chan'));
freqdim = find(strcmp(inputdim, 'freq'));
timedim = find(strcmp(inputdim, 'time'));
if isempty(repldim)
  Nrepl = 1;
else
  Nrepl = size(output.dat, repldim);
end
if isempty(chandim)
  error('data must contain channels');
else
  if hascrsspctrm
    Nchan = size(output.dat, chandim) + size(output.crs, chandim);
    Ncrs  = size(output.crs, chandim);
  else
    Nchan = size(output.dat, chandim);
  end
end
if isempty(freqdim)
  Nfreq = 1;
else
  Nfreq = size(output.dat, freqdim);
end
if isempty(timedim)
  Ntime = 1;
else
  Ntime = size(output.dat, timedim);
end

% find the overlapping dimensions in the input and output
[dum, inputorder, outputorder] = intersect(inputdim, outputdim);
% sort them according to the specified output dimensions
[outputorder, indx] = sort(outputorder);
inputorder = inputorder(indx);

% rearrange the input dimensions in the same order as the desired output
if ~all(inputorder==sort(inputorder))
  if hascrsspctrm
    output.dat = permute(output.dat, inputorder);
    output.crs = permute(output.crs, inputorder);
  else
    output.dat = permute(output.dat, inputorder);
  end
end

% reshape the data so that it contains the (optional) singleton dimensions
if ~(size(output.dat,1)==Nrepl && size(output.dat,2)==Nchan && size(output.dat,3)==Nfreq && size(output.dat,4)==Ntime)
  if hascrsspctrm
    output.dat  = reshape(output.dat, Nrepl, Nchan-Ncrs, Nfreq, Ntime);
    output.crs  = reshape(output.crs, Nrepl, Ncrs      , Nfreq, Ntime);
  else
    output.dat  = reshape(output.dat, Nrepl, Nchan, Nfreq, Ntime);
  end
end

if isfield(input,'dof') && numel(input.dof)==Nrepl*Nfreq*Ntime
    output.dof=reshape(input.dof, Nrepl, Nfreq, Ntime);
end;

% add the fields that describe the axes of the data
if ~isfield(output, 'time')
  output.time = nan;
end
if ~isfield(output, 'freq')
  output.freq = nan;
end
% create replication axis, analogous to time and frequency axis
output.repl   = 1:Nrepl;
output.dimord = 'repl_chan_freq_time';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION convert Nx2 cell array with channel combinations into Nx1 cell array with labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label = cmb2label(labelcmb);
label = {};
for i=1:size(labelcmb)
  label{i,1} = [labelcmb{i,1} '&' labelcmb{i,2}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION  convert Nx1 cell array with labels into N2x cell array with channel combinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labelcmb = label2cmb(label);
labelcmb = {};
for i=1:length(label)
  labelcmb(i,1:2) = tokenize(label{i}, '&');
end
