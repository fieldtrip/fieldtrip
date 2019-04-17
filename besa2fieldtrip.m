function data = besa2fieldtrip(input)

% BESA2FIELDTRIP reads and converts various BESA datafiles into a FieldTrip
% data structure, which subsequently can be used for statistical analysis
% or other analysis methods implemented in Fieldtrip.
%
% Use as
%   [data] = besa2fieldtrip(filename)
% where the filename should point to a BESA datafile (or data that
% was exported by BESA). The output is a MATLAB structure that is
% compatible with FieldTrip.
%
% The format of the output structure depends on the type of datafile:
%   *.avr is converted to a structure similar to the output of FT_TIMELOCKANALYSIS
%   *.mul is converted to a structure similar to the output of FT_TIMELOCKANALYSIS
%   *.swf is converted to a structure similar to the output of FT_TIMELOCKANALYSIS (*)
%   *.tfc is converted to a structure similar to the output of FT_FREQANALYSIS     (*)
%   *.dat is converted to a structure similar to the output of FT_SOURCANALYSIS
%   *.dat combined with a *.gen or *.generic is converted to a structure similar to the output of FT_PREPROCESSING
%
% (*) If the BESA toolbox by Karsten Hochstatter is found on your
% MATLAB path, the readBESAxxx functions will be used (where xxx=tfc/swf),
% alternatively the private functions from FieldTrip will be used.
%
% See also EEGLAB2FIELDTRIP, SPM2FIELDTRIP

% Copyright (C) 2005-2010, Robert Oostenveld
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

if isstruct(input) && numel(input)>1
  % use a recursive call to convert multiple inputs
  data = cell(size(input));
  for i=1:numel(input)
    data{i} = besa2fieldtrip(input(i));
  end
  return
end

if isstruct(input)
  fprintf('besa2fieldtrip: converting structure\n');

  %---------------------TFC-------------------------------------------------%
  if strcmp(input.structtype, 'besa_tfc')
    %fprintf('BESA tfc\n');

    data.time   = input.latencies;
    data.freq   = input.frequencies;
    temp_chans  = char(input.channellabels');
    Nchan       = size(temp_chans,1);
    %{
    if strcmp(input.type, 'COHERENCE_SQUARED')
         % it contains coherence between channel pairs
         fprintf('reading coherence between %d channel pairs\n', Nchan);
         for i=1:Nchan
             tmp = tokenize(deblank(temp_chans(i,:)), '-');
             data.labelcmb{i,1} = deblank(tmp{1});
             data.labelcmb{i,2} = deblank(tmp{2});
             data.label{i,1} = deblank(temp_chans(i,:));
         end
         data.cohspctrm = input.data;
    else
    %}
    % it contains power on channels
    fprintf('reading power on %d channels\n', Nchan);
    for i=1:Nchan
      data.label{i,1} = deblank(temp_chans(i,:));
    end
    data.powspctrm = input.data;
    data.dimord    = 'chan_freq_time';
    data.condition = input.condition; %not original Fieldtrip fieldname

    %end

    clear temp;

    %--------------------Image------------------------------------------------%
  elseif strcmp(input.structtype, 'besa_image')
    %fprintf('BESA image\n');
    data.avg.pow  = input.data;
    xTemp         = input.xcoordinates;
    yTemp         = input.ycoordinates;
    zTemp         = input.zcoordinates;
    data.xgrid    = xTemp;
    data.ygrid    = yTemp;
    data.zgrid    = zTemp;
    nx            = size(data.xgrid,2);
    ny            = size(data.ygrid,2);
    nz            = size(data.zgrid,2);
    % Number of points in each dimension
    data.dim      = [nx ny nz];
    % Array with all possible positions (x,y,z)
    data.pos      = WritePosArray(xTemp,yTemp,zTemp,nx,ny,nz);
    data.inside   = 1:prod(data.dim); %as in Fieldtrip - not correct
    data.outside  = [];

    %--------------------Source Waveform--------------------------------------%
  elseif strcmp(input.structtype, 'besa_sourcewaveforms')
    %fprintf('BESA source waveforms\n');
    data.label         = input.labels'; %not the same as Fieldtrip!
    data.dimord        = 'chan_time';
    data.fsample       = input.samplingrate;
    data.time          = input.latencies / 1000.0;
    data.avg           = input.waveforms';
    data.cfg.filename  = input.datafile;

    %--------------------Data Export------------------------------------------%
  elseif strcmp(input.structtype, 'besa_channels')
    %fprintf('BESA data export\n');

    if isfield(input, 'datatype')
      switch input.ft_datatype
        case {'Raw_Data', 'Epoched_Data', 'Segment'}
          data.fsample    = input.samplingrate;
          data.label      = input.channellabels';
          for k=1:size(input.data,2)
            data.time{1,k}  = input.data(k).latencies / 1000.0';
            data.trial{1,k} = input.data(k).amplitudes';
          end
        otherwise
          fprintf('ft_datatype other than Raw_Data, Epoched or Segment');
      end
    else
      fprintf('workspace created with earlier MATLAB version');
    end

    %--------------------else-------------------------------------------------%
  else
    ft_error('unrecognized format of the input structure');
  end

elseif ischar(input)
  fprintf('besa2fieldtrip: reading from file\n');

  % This function can either use the reading functions included in FieldTrip
  % (with contributions from Karsten, Vladimir and Robert), or the official
  % released functions by Karsten Hoechstetter from BESA. The functions in the
  % official toolbox have precedence.
  hasbesa = ft_hastoolbox('besa',1, 1);

  type = ft_filetype(input);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(type, 'besa_avr') && hasbesa
    fprintf('reading ERP/ERF\n');
    % this should be similar to the output of TIMELOCKANALYSIS
    tmp = readBESAavr(input);
    % convert into a TIMELOCKANALYSIS compatible data structure
    data = [];
    data.label = [];
    if isfield(tmp, 'ChannelLabels')
      data.label = fixlabels(tmp.ChannelLabels);
    end
    data.avg     = tmp.Data;
    data.time    = tmp.Time / 1000; % convert to seconds
    data.fsample = 1000/tmp.DI;
    data.dimord  = 'chan_time';
  elseif strcmp(type, 'besa_avr') && ~hasbesa
    fprintf('reading ERP/ERF\n');
    % this should be similar to the output of TIMELOCKANALYSIS
    tmp = read_besa_avr(input);
    % convert into a TIMELOCKANALYSIS compatible data structure
    data = [];
    data.label   = fixlabels(tmp.label);
    data.avg     = tmp.data;
    data.time    = (0:(tmp.npnt-1)) * tmp.di + tmp.tsb;
    data.time    = data.time / 1000; % convert to seconds
    data.fsample = 1000/tmp.di;
    data.dimord  = 'chan_time';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(type, 'besa_mul') && hasbesa
    fprintf('reading ERP/ERF\n');
    % this should be similar to the output of TIMELOCKANALYSIS
    tmp = readBESAmul(input);
    % convert into a TIMELOCKANALYSIS compatible data structure
    data = [];
    data.label    = tmp.ChannelLabels(:);
    data.avg      = tmp.data';
    data.time     = (0:(tmp.Npts-1)) * tmp.DI + tmp.TSB;
    data.time     = data.time / 1000; %convert to seconds
    data.fsample  = 1000/tmp.DI;
    data.dimord   = 'chan_time';
  elseif strcmp(type, 'besa_mul') && ~hasbesa
    fprintf('reading ERP/ERF\n');
    % this should be similar to the output of TIMELOCKANALYSIS
    tmp = read_besa_mul(input);
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
      ft_error('this data format requires the BESA toolbox');
    end
    [p, f, x] = fileparts(input);
    input = fullfile(p, [f '.dat']);
    [time,buf,ntrial] = readBESAsb(input);
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
    data.fsample = 1./mean(diff(time));  % time is already in seconds

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(type, 'besa_tfc') && hasbesa
    fprintf('reading time-frequency representation using BESA toolbox\n');
    % this should be similar to the output of FREQANALYSIS
    tfc = readBESAtfc(input);
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
    data.condition = tfc.ConditionName;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(type, 'besa_tfc') && ~hasbesa
    fprintf('reading time-frequency representation\n');
    % this should be similar to the output of FREQANALYSIS
    [ChannelLabels, Time, Frequency, Data, Info] = read_besa_tfc(input);
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
    swf = readBESAswf(input);
    % convert into a TIMELOCKANALYSIS compatible data structure
    data         = [];
    data.label   = fixlabels(swf.waveName);
    data.avg     = swf.data;
    data.time    = swf.Time * 1e-3; % convert to seconds
    data.fsample = 1/mean(diff(data.time));
    data.dimord  = 'chan_time';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(type, 'besa_swf') && ~hasbesa
    fprintf('reading source waveform\n');
    % hmm, I guess that this should be similar to the output of TIMELOCKANALYSIS
    tmp = read_besa_swf(input);
    % convert into a TIMELOCKANALYSIS compatible data structure
    data = [];
    data.label   = fixlabels(tmp.label);
    data.avg     = tmp.data;
    data.time    = (0:(tmp.npnt-1)) * tmp.di + tmp.tsb;
    data.time    = data.time / 1000; % convert to seconds
    data.fsample = 1000/tmp.di;
    data.dimord  = 'chan_time';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(type, 'besa_src') && hasbesa
    src = readBESAimage(input);
    data.xgrid = src.Coordinates.X;
    data.ygrid = src.Coordinates.Y;
    data.zgrid = src.Coordinates.Z;
    data.avg.pow = src.Data;
    data.dim = size(src.Data);
    [X, Y, Z] = ndgrid(data.xgrid, data.ygrid, data.zgrid);
    data.pos = [X(:) Y(:) Z(:)];
    % cannot determine which voxels are inside the brain volume
    data.inside = 1:prod(data.dim);
    data.outside = [];
  elseif strcmp(type, 'besa_src') && ~hasbesa
    src = read_besa_src(input);
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
    ft_error('sorry, pdg is not yet supported');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else
    ft_error('unrecognized file format for importing BESA data');
  end

end % isstruct || ischar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that fixes the channel labels, should be a cell-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newlabels] = fixlabels(labels)
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
  newlabels = tokenize(labels(:)', ' '); % also ensure that it is a row-string
elseif ischar(labels) && ~any(size(labels)==1)
  for i=1:size(labels)
    newlabels{i} = strtrim(labels(i,:));
  end
end
% convert to column
newlabels = newlabels(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PArray] = WritePosArray(x,y,z,mx,my,mz)
A1  = repmat(x,1,my*mz);
A21 = repmat(y,mx,mz);
A2  = reshape(A21,1,mx*my*mz);
A31 = repmat(z,mx*my,1);
A3  = reshape(A31,1,mx*my*mz);
PArray = [A1;A2;A3]';
