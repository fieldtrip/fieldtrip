function [data] = besa2fieldtrip(filename);

% BESA2FIELDTRIP reads and converts various BESA datafiles into a FieldTrip
% data structure, which subsequently can be used for statistical analysis
% or other analysis methods implemented in Fieldtrip.
%
% Use as
%   [data] = besa2fieldtrip(filename)
% where the filename should point to a BESA datafile (or data that
% is exported by BESA). The output is a Matlab structure that is
% compatible with FieldTrip.
%
% The format of the output structure depends on the type of datafile:
%   *.avr is converted to a structure similar to the output of TIMELOCKANALYSIS
%   *.mul is converted to a structure similar to the output of TIMELOCKANALYSIS
%   *.swf is converted to a structure similar to the output of TIMELOCKANALYSIS (*)
%   *.tfc is converted to a structure similar to the output of FREQANALYSIS     (*)
%   *.dat is converted to a structure similar to the output of SOURCANALYSIS
%   *.dat combined with a *.gen or *.generic is converted to a structure similar to the output of PREPROCESSING
%
% Note (*): If the BESA toolbox by Karsten Hochstatter is found on your
% Matlab path, the readBESAxxx functions will be used (where xxx=tfc/swf),
% alternatively the private functions from FieldTrip will be used.
%
% See also EEGLAB2FIELDTRIP

% Undocumented local options:
% cfg.filename
% cfg.version

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: besa2fieldtrip.m,v $
% Revision 1.16  2008/10/21 09:34:12  roboos
% fixed sampling frequency (factor 1000 wrgon), thanks to Stephan Moratti
%
% Revision 1.15  2008/09/22 21:23:11  roboos
% always add besa path (now as external/besa)
%
% Revision 1.14  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.13  2007/04/24 09:45:56  roboos
% added support for besa simple binary, requires besa toolbox
%
% Revision 1.12  2006/10/05 08:54:14  roboos
% convert time to seconds for TFR data
%
% Revision 1.11  2006/06/27 13:26:51  roboos
% use readBESAswf from besa toolbox if available, otherwise use fieldtrip/private function
%
% Revision 1.10  2006/06/07 09:34:19  roboos
% changed checktoolbox into hastoolbox
%
% Revision 1.9  2006/05/15 13:18:23  roboos
% implemented fixlabels subfunction, use it to convert between char and cell-array
%
% Revision 1.8  2006/04/26 11:36:52  roboos
% detect the presence of the BESA toolbox and added support for reading
% *.tfc using Karstens readBESAtfc function
%
% Revision 1.7  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.6  2006/03/16 17:36:13  roboos
% added besa_swf
%
% Revision 1.5  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.4  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.3  2005/09/01 09:25:33  roboos
% added beamformer source reconstruction
%
% Revision 1.2  2005/07/29 13:23:58  roboos
% various improvements and fixes
%
% Revision 1.1  2005/07/28 15:08:58  roboos
% new implementation, see documentation on website
%

fieldtripdefs

% This function can either use the reading functions included in FieldTrip
% (with contributions from Karsten, Vladimir and Robert), or the official
% released functions by Karsten Hoechstetter from BESA. The functions in the
% official toolbox have precedence.
hasbesa = hastoolbox('besa',1, 1);

type = filetype(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type, 'besa_avr')
  fprintf('reading ERP/ERF\n');
  % this should be similar to the output of TIMELOCKANALYSIS
  tmp = read_besa_avr(filename);
  % convert into a TIMELOCKANALYSIS compatible data structure
  data = [];
  data.label   = fixlabels(tmp.label);
  data.avg     = tmp.data;
  data.time    = (0:(tmp.npnt-1)) * tmp.di + tmp.tsb;
  data.time    = data.time / 1000; % convert to seconds
  data.fsample = 1000/tmp.di;
  data.dimord  = 'chan_time';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_mul')
  fprintf('reading ERP/ERF\n');
  % this should be similar to the output of TIMELOCKANALYSIS
  tmp = read_besa_mul(filename);
  % convert into a TIMELOCKANALYSIS compatible data structure
  data = [];
  data.label   = tmp.label(:);
  data.avg     = tmp.data;
  data.time    = (0:(tmp.TimePoints-1)) * tmp.SamplingInterval_ms_ + tmp.BeginSweep_ms_;
  data.time    = data.time / 1000; % convert to seconds
  data.fsample = 1000/tmp.SamplingInterval_ms_;
  data.dimord  = 'chan_time';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_sb') 
  if hasbesa
    fprintf('reading preprocessed channel data using BESA toolbox\n');
  else
    error('this data format requires the BESA toolbox');
  end
  [p, f, x] = fileparts(filename);
  filename = fullfile(p, [f '.dat']);
  [time,buf,ntrial] = readBESAsb(filename);
  time  = time/1000;                 % convert from ms to sec
  nchan = size(buf,1);
  ntime = size(buf,3);

  % convert into a PREPROCESSING compatible data structure
  data       = [];
  data.trial = {};
  data.time  = {};
  for i=1:ntrial
    data.trial{i} = reshape(buf(:,i,:), [nchan, ntime]);
    data.time{i} = time;
  end
  data.label   = {};
  for i=1:size(buf,1)
    data.label{i,1} = sprintf('chan%03d', i);
  end
  data.fsample = 1/(time(2)-time(1));  % time is already in seconds

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_tfc') && hasbesa
  fprintf('reading time-frequency representation using BESA toolbox\n');
  % this should be similar to the output of FREQANALYSIS
  tfc = readBESAtfc(filename);
  Nchan = size(tfc.ChannelLabels,1);
  % convert into a FREQANALYSIS compatible data structure
  data = [];
  data.time = tfc.Time(:)';
  data.freq = tfc.Frequency(:)';
  if isfield(tfc, 'DataType') && strcmp(tfc.DataType, 'COHERENCE_SQUARED')
    % it contains coherence between channel pairs
    fprintf('reading coherence between %d channel pairs\n', Nchan);
    for i=1:Nchan
      tmp = tokenize(deblank(tfc.ChannelLabels(i,:)), '-');
      data.labelcmb{i,1} = tmp{1};
      data.labelcmb{i,2} = tmp{2};
    end
    data.cohspctrm = permute(tfc.Data, [1 3 2]);
  else
    % it contains power on channels
    fprintf('reading power on %d channels\n', Nchan);
    for i=1:Nchan
      data.label{i,1} = deblank(tfc.ChannelLabels(i,:));
    end
    data.powspctrm = permute(tfc.Data, [1 3 2]);
  end
  data.dimord    = 'chan_freq_time';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_tfc') && ~hasbesa
  fprintf('reading time-frequency representation\n');
  % this should be similar to the output of FREQANALYSIS
  [ChannelLabels, Time, Frequency, Data, Info] = read_besa_tfc(filename);
  Nchan = size(ChannelLabels,1);
  % convert into a FREQANALYSIS compatible data structure
  data = [];
  data.time = Time * 1e-3; % convert to seconds;
  data.freq = Frequency;
  if isfield(Info, 'DataType') && strcmp(Info.DataType, 'COHERENCE_SQUARED')
    % it contains coherence between channel pairs
    fprintf('reading coherence between %d channel pairs\n', Nchan);
    for i=1:Nchan
      tmp = tokenize(deblank(ChannelLabels(i,:)), '-');
      data.labelcmb{i,1} = tmp{1};
      data.labelcmb{i,2} = tmp{2};
    end
    data.cohspctrm = permute(Data, [1 3 2]);
  else
    % it contains power on channels
    fprintf('reading power on %d channels\n', Nchan);
    for i=1:Nchan
      data.label{i} = deblank(ChannelLabels(i,:));
    end
    data.powspctrm = permute(Data, [1 3 2]);
  end
  data.dimord    = 'chan_freq_time';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_swf') && hasbesa
  fprintf('reading source waveform using BESA toolbox\n');
  swf = readBESAswf(filename);
  % convert into a TIMELOCKANALYSIS compatible data structure
  data         = [];
  data.label   = fixlabels(swf.waveName);
  data.avg     = swf.data;
  data.time    = swf.Time * 1e-3; % convert to seconds
  data.fsample = 1/(data.time(2)-data.time(1));
  data.dimord  = 'chan_time';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_swf') && ~hasbesa
  fprintf('reading source waveform\n');
  % hmm, I guess that this should be similar to the output of TIMELOCKANALYSIS
  tmp = read_besa_swf(filename);
  % convert into a TIMELOCKANALYSIS compatible data structure
  data = [];
  data.label   = fixlabels(tmp.label);
  data.avg     = tmp.data;
  data.time    = (0:(tmp.npnt-1)) * tmp.di + tmp.tsb;
  data.time    = data.time / 1000; % convert to seconds
  data.fsample = 1000/tmp.di;
  data.dimord  = 'chan_time';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_src')
  src = read_besa_src(filename);
  data.xgrid = linspace(src.X(1), src.X(2), src.X(3));
  data.ygrid = linspace(src.Y(1), src.Y(2), src.Y(3));
  data.zgrid = linspace(src.Z(1), src.Z(2), src.Z(3));
  data.avg.pow = src.vol;
  data.dim = size(src.vol);
  [X, Y, Z] = ndgrid(data.xgrid, data.ygrid, data.zgrid);
  data.pos = [X(:) Y(:) Z(:)];
  % cannot determine which voxels are inside the brain volume
  data.inside = 1:prod(data.dim);
  data.outside = [];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'besa_pdg')
  % hmmm, I have to think about this one...
  error('sorry, pdg is not yet supported');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error('unrecognized file format for importing BESA data');
end

% construct and add a configuration to the output
cfg = [];
cfg.filename = filename;
% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: besa2fieldtrip.m,v 1.16 2008/10/21 09:34:12 roboos Exp $';
data.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that fixes the channel labels, should be a cell-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newlabels] = fixlabels(labels);
if iscell(labels) && length(labels)>1
  % seems to be ok
  newlabels = labels;
elseif iscell(labels) && length(labels)==1
  % could be a cell with a single long string in it
  if length(tokenize(labels{1}, ' '))>1
    % seems like a long string that accidentaly ended up in a single
    cell
    newlabels = tokenize(labels{1}, ' ');
  else
    % seems to be ok
    newlabels = labels;
  end
elseif ischar(labels) && any(size(labels)==1)
  newlabels = tokenize(labels(:)', ' '); % also ensure that it is a
  row-string
elseif ischar(labels) && ~any(size(labels)==1)
  for i=1:size(labels)
    newlabels{i} = fliplr(deblank(fliplr(deblank(labels(i,:)))));
  end
end
% convert to column
newlabels = newlabels(:);
