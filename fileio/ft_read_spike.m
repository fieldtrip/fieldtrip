function [spike] = ft_read_spike(filename, varargin);

% FT_READ_SPIKE reads spike timestamps and waveforms from various data
% formats.
%
% Use as
%  [spike] = ft_read_spike(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'spikeformat'
%
% The output spike structure contains
%   spike.label     = 1xNchans cell-array, with channel labels
%   spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nsamples X Nspikes)
%   spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
%   spike.unit      = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
%
% See also FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT

% Copyright (C) 2007-2010 Robert Oostenveld
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
% $Id$

% get the options
spikeformat   = keyval('spikeformat',   varargin);

% determine the filetype
if isempty(spikeformat)
  spikeformat = ft_filetype(filename);
end

switch spikeformat
  case {'neuralynx_ncs' 'plexon_ddt'}
    % these files only contain continuous data
    error('file does not contain spike timestamps or waveforms');

  case 'matlab'
    % plain matlab file with a single variable in it
    load(filename, 'spike');

  case 'mclust_t'
    fp = fopen(filename, 'rb', 'ieee-le');
    H = ReadHeader(fp);
    fclose(fp);
    % read only from one file
    S = LoadSpikes({filename});
    spike.hdr = H(:);
    spike.timestamp = S;
    [p, f, x] = fileparts(filename);
    spike.label     = {f};  % use the filename as label for the spike channel
    spike.waveform  = {};   % this is unknown
    spike.unit      = {};   % this is unknown

  case 'neuralynx_nse'
    % single channel file, read all records
    nse = read_neuralynx_nse(filename);
    if isfield(nse.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nse.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nse.hdr.AcqEntName};
    end
    spike.timestamp = {nse.TimeStamp};
    spike.waveform  = {nse.dat};
    spike.unit      = {nse.CellNumber};
    spike.hdr       = nse.hdr;

  case 'neuralynx_nst'
    % single channel stereotrode file, read all records
    nst = read_neuralynx_nst(filename, 1, inf);
    if isfield(nst.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nst.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nst.hdr.AcqEntName};
    end
    spike.timestamp = {nst.TimeStamp};
    spike.waveform  = {nst.dat};
    spike.unit      = {nst.CellNumber};
    spike.hdr       = nst.hdr;

  case 'neuralynx_ntt'
    % single channel stereotrode file, read all records
    ntt = read_neuralynx_ntt(filename);
    if isfield(ntt.hdr, 'NLX_Base_Class_Name')
      spike.label   = {ntt.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {ntt.hdr.AcqEntName};
    end
    spike.timestamp = {ntt.TimeStamp};
    spike.waveform  = {ntt.dat};
    spike.unit      = {ntt.CellNumber};
    spike.hdr       = ntt.hdr;

  case 'neuralynx_nts'
    % single channel file, read all records
    nts = read_neuralynx_nts(filename);
    if isfield(nte.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nts.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nts.hdr.AcqEntName};
    end
    spike.timestamp = {nts.TimeStamp(:)'};
    spike.waveform  = {zeros(0,length(nts.TimeStamp))};  % does not contain waveforms
    spike.unit      = {zeros(0,length(nts.TimeStamp))};  % does not contain units
    spike.hdr       = nts.hdr;

  case 'plexon_nex'
    % a single file can contain multiple channels of different types
    hdr  = read_plexon_nex(filename);
    typ  = [hdr.VarHeader.Type];
    chan = 0;

    spike.label     = {};
    spike.waveform  = {};
    spike.unit      = {};
    spike.timestamp = {};

    for i=1:length(typ)
      if typ(i)==0
        % neurons, only timestamps
        nex = read_plexon_nex(filename, 'channel', i);
        nspike = length(nex.ts);
        chan = chan + 1;
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = zeros(0, nspike);
        spike.unit{chan}      = nan*ones(1,nspike);
        spike.timestamp{chan} = nex.ts;
      elseif typ(i)==3
        % neurons, timestamps and waveforms
        nex = read_plexon_nex(filename, 'channel', i);
        chan = chan + 1;
        nspike = length(nex.ts);
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = nex.dat;
        spike.unit{chan}      = nan*ones(1,nspike);
        spike.timestamp{chan} = nex.ts;
      end
    end
    spike.hdr = hdr;

  case 'plexon_plx'
    % read the header information
    hdr   = read_plexon_plx(filename);
    nchan = length(hdr.ChannelHeader);
    typ   = [hdr.DataBlockHeader.Type];
    unit  = [hdr.DataBlockHeader.Unit];
    chan  = [hdr.DataBlockHeader.Channel];

    for i=1:nchan
      % select the data blocks that contain spike waveforms and that belong to this channel
      sel = (typ==1 & chan==hdr.ChannelHeader(i).Channel);

      if any(sel)
        % get the timestamps that correspond with this spike channel
        tsl = [hdr.DataBlockHeader(sel).TimeStamp];
        tsh = [hdr.DataBlockHeader(sel).UpperByteOf5ByteTimestamp];
        % convert the 16 bit high timestamp into a 32 bit integer
        ts = timestamp_plexon(tsl, tsh);
        spike.timestamp{i} = ts;
        spike.unit{i}      = unit(sel);
      else
        % this spike channel is empty
        spike.timestamp{i} = [];
        spike.unit{i}      = [];
      end
    end
    for i=1:nchan
      spike.label{i}    = deblank(hdr.ChannelHeader(i).Name);
      spike.waveform{i} = read_plexon_plx(filename, 'ChannelIndex', i, 'header', hdr);
    end
    spike.hdr = hdr;
    
  case 'neuroshare' % NOTE: still under development
    % check that the required neuroshare toolbox is available
    ft_hastoolbox('neuroshare', 1);

    tmp = read_neuroshare(filename, 'readspike', 'yes');
    spike.label = {tmp.hdr.entityinfo(tmp.list.segment).EntityLabel};
    for i=1:length(spike.label)
      spike.waveform{i}  = tmp.spikew.data(:,:,i);
      spike.timestamp{i} = tmp.spikew.timestamp(:,i)';
      spike.unit{i}      = tmp.spikew.unitID(:,i)';
    end

  otherwise
    error('unsupported data format');
end

