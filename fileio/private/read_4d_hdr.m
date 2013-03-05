function [header] = read_4d_hdr(datafile, configfile)

% hdr=READ_4D_HDR(datafile, configfile)
% Collects the required Fieldtrip header data from the data file 'filename'
% and the associated 'config' file for that data.
%
% Adapted from the MSI>>Matlab code written by Eugene Kronberg

% Copyright (C) 2008-2009, Centre for Cognitive Neuroimaging, Glasgow, Gavin Paterson & J.M.Schoffelen
% Copyright (C) 2010-2011, Donders Institute for Brain, Cognition and Behavior, J.M.Schoffelen
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

%read header
if nargin ~= 2
  [path, file, ext] = fileparts(datafile);
  configfile        = fullfile(path, 'config');
end

if ~isempty(datafile),
  %always big endian  
  fid = fopen(datafile, 'r', 'b');
  
  if fid == -1
    error('Cannot open file %s', datafile);
  end
  
  fseek(fid, 0, 'eof');
  header_end = ftell(fid);
  %last 8 bytes of the pdf is header offset
  fseek(fid, -8, 'eof');
  header_offset = fread(fid,1,'uint64');
  
  %first byte of the header
  fseek(fid, header_offset, 'bof');
  
  % read header data
  align_file_pointer(fid)
  header.header_data.FileType  = fread(fid, 1, 'uint16=>uint16');
  file_type                    = char(fread(fid, 5, 'uchar'))';
  header.header_data.file_type = file_type(file_type>0);
  fseek(fid, 1, 'cof');
  format = fread(fid, 1, 'int16=>int16');
  switch format
      case 1
          header.header_data.Format = 'SHORT';
      case 2
          header.header_data.Format = 'LONG';
      case 3
          header.header_data.Format = 'FLOAT';
      case 4
          header.header_data.Format ='DOUBLE';
  end
  header.header_data.acq_mode           = fread(fid, 1,  'uint16=>uint16');
  header.header_data.TotalEpochs        = fread(fid, 1,  'uint32=>double');
  header.header_data.input_epochs       = fread(fid, 1,  'uint32=>uint32');
  header.header_data.TotalEvents        = fread(fid, 1,  'uint32=>uint32');
  header.header_data.total_fixed_events = fread(fid, 1,  'uint32=>uint32');
  header.header_data.SamplePeriod       = fread(fid, 1,  'float32=>float64');
  header.header_data.SampleFrequency    = 1/header.header_data.SamplePeriod;
  xaxis_label                           = char(fread(fid, 16, 'uchar'))';
  header.header_data.xaxis_label        = xaxis_label(xaxis_label>0);
  header.header_data.total_processes    = fread(fid, 1,  'uint32=>uint32');
  header.header_data.TotalChannels      = fread(fid, 1,  'uint16=>double');
  fseek(fid, 2, 'cof');
  header.header_data.checksum           = fread(fid, 1,  'int32=>int32');
  header.header_data.total_ed_classes   = fread(fid, 1,  'uint32=>uint32');
  header.header_data.total_associated_files = fread(fid, 1, 'uint16=>uint16');
  header.header_data.last_file_index    = fread(fid, 1,  'uint16=>uint16');
  header.header_data.timestamp          = fread(fid, 1,  'uint32=>uint32');
  header.header_data.reserved           = fread(fid, 20, 'uchar')';
  fseek(fid, 4, 'cof');
   
  %read epoch_data
  for epoch = 1:header.header_data.TotalEpochs;
    align_file_pointer(fid)
  
    header.epoch_data(epoch).pts_in_epoch     = fread(fid, 1,  'uint32=>uint32');
    header.epoch_data(epoch).epoch_duration   = fread(fid, 1,  'float32=>float32');
    header.epoch_data(epoch).expected_iti     = fread(fid, 1,  'float32=>float32');
    header.epoch_data(epoch).actual_iti       = fread(fid, 1,  'float32=>float32');
    header.epoch_data(epoch).total_var_events = fread(fid, 1,  'uint32=>uint32');
    header.epoch_data(epoch).checksum         = fread(fid, 1,  'int32=>int32');
    header.epoch_data(epoch).epoch_timestamp  = fread(fid, 1,  'int32=>int32');
    header.epoch_data(epoch).reserved         = fread(fid, 28, 'uchar')';
    header.header_data.SlicesPerEpoch         = double(header.epoch_data(1).pts_in_epoch);
  
    %read event data (var_events)
    for event = 1:header.epoch_data(epoch).total_var_events
      align_file_pointer(fid)
  
      event_name = char(fread(fid, 16, 'uchar'))';
      header.epoch_data(epoch).var_event{event}.event_name  = event_name(event_name>0);
      header.epoch_data(epoch).var_event{event}.start_lat   = fread(fid, 1, 'float32=>float32');
      header.epoch_data(epoch).var_event{event}.end_lat     = fread(fid, 1, 'float32=>float32');
      header.epoch_data(epoch).var_event{event}.step_size   = fread(fid, 1, 'float32=>float32');
      header.epoch_data(epoch).var_event{event}.fixed_event = fread(fid, 1, 'uint16=>uint16');
      fseek(fid, 2, 'cof');
      header.epoch_data(epoch).var_event{event}.checksum    = fread(fid, 1, 'int32=>int32');
      header.epoch_data(epoch).var_event{event}.reserved    = fread(fid, 32, 'uchar')';
      fseek(fid, 4, 'cof');
    end
  end
  
  %read  channel ref data
  for channel = 1:header.header_data.TotalChannels
    align_file_pointer(fid)
 
    chan_label                                 = (fread(fid, 16, 'uint8=>char'))';
    header.channel_data(channel).chan_label    = chan_label(chan_label>0);
    header.channel_data(channel).chan_no       = fread(fid, 1, 'uint16=>uint16');
    header.channel_data(channel).attributes    = fread(fid, 1, 'uint16=>uint16');
    header.channel_data(channel).scale         = fread(fid, 1, 'float32=>float32');
    yaxis_label                                = char(fread(fid, 16, 'uint8=>char'))';
    header.channel_data(channel).yaxis_label   = yaxis_label(yaxis_label>0);
    header.channel_data(channel).valid_min_max = fread(fid, 1, 'uint16=>uint16');
    fseek(fid, 6, 'cof');
    header.channel_data(channel).ymin          = fread(fid, 1,  'float64');
    header.channel_data(channel).ymax          = fread(fid, 1,  'float64');
    header.channel_data(channel).index         = fread(fid, 1,  'uint32=>uint32');
    header.channel_data(channel).checksum      = fread(fid, 1,  'int32=>int32');
    header.channel_data(channel).whatisit      = char(fread(fid, 4, 'uint8=>char'))';
    header.channel_data(channel).reserved      = fread(fid, 28, 'uint8')';
  end
  
  %read event data
  for event = 1:header.header_data.total_fixed_events
    align_file_pointer(fid)
    event_name                           = char(fread(fid, 16, 'uchar'))';
    header.event_data(event).event_name  = event_name(event_name>0);
    header.event_data(event).start_lat   = fread(fid, 1,  'float32=>float32');
    header.event_data(event).end_lat     = fread(fid, 1,  'float32=>float32');
    header.event_data(event).step_size   = fread(fid, 1,  'float32=>float32');
    header.event_data(event).fixed_event = fread(fid, 1,  'uint16=>uint16');
    fseek(fid, 2, 'cof');
    header.event_data(event).checksum    = fread(fid, 1,  'int32=>int32');
    header.event_data(event).reserved    = fread(fid, 32, 'uchar')';
    fseek(fid, 4, 'cof');
  end
  header.header_data.FirstLatency = double(header.event_data(1).start_lat);
  
  %experimental: read process information
  for np = 1:header.header_data.total_processes
    align_file_pointer(fid)
    nbytes                          = fread(fid, 1,  'uint32=>uint32');
    fp                              = ftell(fid);
    header.process(np).hdr.nbytes   = nbytes;
    type                            = char(fread(fid, 20, 'uchar'))';
    header.process(np).hdr.type     = type(type>0);
    header.process(np).hdr.checksum = fread(fid, 1,  'int32=>int32'); 
    user                            = char(fread(fid, 32, 'uchar'))';
    header.process(np).user         = user(user>0);
    header.process(np).timestamp    = fread(fid, 1,  'uint32=>uint32');
    fname                           = char(fread(fid, 32, 'uchar'))';
    header.process(np).filename     = fname(fname>0);
    fseek(fid, 28*8, 'cof'); %dont know
    header.process(np).totalsteps   = fread(fid, 1,  'uint32=>uint32');
    header.process(np).checksum     = fread(fid, 1,  'int32=>int32');
    header.process(np).reserved     = fread(fid, 32, 'uchar')'; 
    for ns = 1:header.process(np).totalsteps
      align_file_pointer(fid)
      nbytes2                                   = fread(fid, 1, 'uint32=>uint32');
      header.process(np).step(ns).hdr.nbytes    = nbytes2;
      type                                      = char(fread(fid, 20, 'uchar'))';
      header.process(np).step(ns).hdr.type      = type(type>0); %dont know how to interpret the first two
      header.process(np).step(ns).hdr.checksum  = fread(fid, 1, 'int32=>int32');
      userblocksize                             = fread(fid, 1, 'int32=>int32'); %we are at 32 bytes here
      header.process(np).step(ns).userblocksize = userblocksize;
      fseek(fid, nbytes2 - 32, 'cof');
      
      if strcmp(header.process(np).step(ns).hdr.type, 'PDF_Weight_Table'),
        warning('reading in weight table: no warranty that this is correct. it seems to work for the Glasgow 248-magnetometer system. if you have some code yourself, and/or would like to test it on your own data, please contact Jan-Mathijs');
        tmpfp = ftell(fid);
        tmp   = fread(fid, 1, 'uint8');
        Nchan = fread(fid, 1, 'uint32');
        Nref  = fread(fid, 1, 'uint32');
        for k = 1:Nref
          name = fread(fid, 17, 'uchar'); %strange number, but seems to be true
          header.process(np).step(ns).RefChan{k,1} = char(name(name>0))';
        end
        fseek(fid, 152, 'cof');
        for k = 1:Nchan
          name = fread(fid, 17, 'uchar');
          header.process(np).step(ns).Chan{k,1}   = char(name(name>0))';
        end
        %fseek(fid, 20, 'cof');
        %fseek(fid, 4216, 'cof');
        header.process(np).step(ns).stuff1  = fread(fid, 4236, 'uint8');
        name                                = fread(fid, 16, 'uchar');
        header.process(np).step(ns).Creator = char(name(name>0))';
        %some stuff I don't understand yet
        %fseek(fid, 136, 'cof');
        header.process(np).step(ns).stuff2  = fread(fid, 136, 'uint8');
        %now something strange is going to happen: the weights are probably little-endian encoded.
        %here we go: check whether this applies to the whole PDF weight table
        fp = ftell(fid);
        fclose(fid);
        fid = fopen(datafile, 'r', 'l');
        fseek(fid, fp, 'bof');
        for k = 1:Nchan
          header.process(np).step(ns).Weights(k,:) = fread(fid, 23, 'float32=>float32')';
          fseek(fid, 36, 'cof');
        end
      else
        if userblocksize < 1e6,
          %for one reason or another userblocksize can assume strangely high values
          fseek(fid, userblocksize, 'cof');
        end
      end    
    end
  end
  fclose(fid);
end 
%end read header

%read config file
fid = fopen(configfile, 'r', 'b');

if fid == -1
  error('Cannot open config file');
end

header.config_data.version           = fread(fid, 1, 'uint16=>uint16');
site_name                            = char(fread(fid, 32, 'uchar'))';
header.config_data.site_name         = site_name(site_name>0);
dap_hostname                         = char(fread(fid, 16, 'uchar'))';
header.config_data.dap_hostname      = dap_hostname(dap_hostname>0);
header.config_data.sys_type          = fread(fid, 1, 'uint16=>uint16');
header.config_data.sys_options       = fread(fid, 1, 'uint32=>uint32');
header.config_data.supply_freq       = fread(fid, 1, 'uint16=>uint16');
header.config_data.total_chans       = fread(fid, 1, 'uint16=>uint16');
header.config_data.system_fixed_gain = fread(fid, 1, 'float32=>float32');
header.config_data.volts_per_bit     = fread(fid, 1, 'float32=>float32');
header.config_data.total_sensors     = fread(fid, 1, 'uint16=>uint16');
header.config_data.total_user_blocks = fread(fid, 1, 'uint16=>uint16');
header.config_data.next_derived_channel_number = fread(fid, 1, 'uint16=>uint16');
fseek(fid, 2, 'cof');
header.config_data.checksum          = fread(fid, 1, 'int32=>int32');
header.config_data.reserved          = fread(fid, 32, 'uchar=>uchar')';

header.config.Xfm = fread(fid, [4 4], 'double');

%user blocks
for ub = 1:header.config_data.total_user_blocks
  align_file_pointer(fid)
  header.user_block_data{ub}.hdr.nbytes   = fread(fid, 1, 'uint32=>uint32');
  type                                    = char(fread(fid, 20, 'uchar'))';
  header.user_block_data{ub}.hdr.type     = type(type>0);
  header.user_block_data{ub}.hdr.checksum = fread(fid, 1, 'int32=>int32');
  user                                    = char(fread(fid, 32, 'uchar'))';
  header.user_block_data{ub}.user         = user(user>0);
  header.user_block_data{ub}.timestamp    = fread(fid, 1, 'uint32=>uint32');
  header.user_block_data{ub}.user_space_size = fread(fid, 1, 'uint32=>uint32');
  header.user_block_data{ub}.reserved     = fread(fid, 32, 'uchar=>uchar')';
  fseek(fid, 4, 'cof');
  user_space_size                         = double(header.user_block_data{ub}.user_space_size);
  if strcmp(type(type>0), 'B_weights_used'), 
    %warning('reading in weight table: no warranty that this is correct. it seems to work for the Glasgow 248-magnetometer system. if you have some code yourself, and/or would like to test it on your own data, please contact Jan-Mathijs');
    tmpfp = ftell(fid);
    %read user_block_data weights
    %there is information in the 4th and 8th byte, these might be related to the settings?
    version  = fread(fid, 1, 'uint32');
    header.user_block_data{ub}.version = version;
    if version==1,
      Nbytes   = fread(fid,1,'uint32');
      Nchan    = fread(fid,1,'uint32');
      Position = fread(fid, 32, 'uchar');
      header.user_block_data{ub}.position = char(Position(Position>0))';
      fseek(fid,tmpfp+user_space_size - Nbytes*Nchan, 'bof');
      Ndigital = floor((Nbytes - 4*2) / 4);
      Nanalog  = 3; %lucky guess?
      % how to know number of analog weights vs digital weights???
      for ch = 1:Nchan
        % for Konstanz -- comment for others?
        header.user_block_data{ub}.aweights(ch,:) = fread(fid, [1 Nanalog],  'int16')'; 
        fseek(fid,2,'cof'); % alignment
        header.user_block_data{ub}.dweights(ch,:) = fread(fid, [1 Ndigital], 'single=>double')';
      end
      fseek(fid, tmpfp, 'bof');
      %there is no information with respect to the channels here.
      %the best guess would be to assume the order identical to the order in header.config.channel_data
      %for the digital weights it would be the order of the references in that list
      %for the analog weights I would not know
    elseif version==2,
      unknown2 = fread(fid, 1, 'uint32');
      Nchan    = fread(fid, 1, 'uint32');
      Position = fread(fid, 32, 'uchar');
      header.user_block_data{ub}.position = char(Position(Position>0))';
      fseek(fid, tmpfp+124, 'bof');
      Nanalog  = fread(fid, 1, 'uint32');
      Ndigital = fread(fid, 1, 'uint32');
      fseek(fid, tmpfp+204, 'bof');
      for k = 1:Nchan
        Name     = fread(fid, 16, 'uchar');
        header.user_block_data{ub}.channames{k,1} = char(Name(Name>0))';
      end
      for k = 1:Nanalog
        Name     = fread(fid, 16, 'uchar');
        header.user_block_data{ub}.arefnames{k,1} = char(Name(Name>0))';
      end
      for k = 1:Ndigital
        Name     = fread(fid, 16, 'uchar');
        header.user_block_data{ub}.drefnames{k,1} = char(Name(Name>0))';
      end

      header.user_block_data{ub}.dweights = fread(fid, [Ndigital Nchan], 'single=>double')';
      header.user_block_data{ub}.aweights = fread(fid, [Nanalog  Nchan],  'int16')'; 
      fseek(fid, tmpfp, 'bof');
    end
  elseif strcmp(type(type>0), 'B_E_table_used'),
    %warning('reading in weight table: no warranty that this is correct');
    %tmpfp = ftell(fid);
    %fseek(fid, 4, 'cof'); %there's info here dont know how to interpret
    %Nx    = fread(fid, 1, 'uint32');
    %Nchan = fread(fid, 1, 'uint32');
    %type  = fread(fid, 32, 'uchar'); %don't know whether correct
    %header.user_block_data{ub}.type = char(type(type>0))';
    %fseek(fid, 16, 'cof');
    %for k = 1:Nchan
    %  name                                 = fread(fid, 16, 'uchar');
    %  header.user_block_data{ub}.name{k,1} = char(name(name>0))';
    %end
  elseif strcmp(type(type>0), 'B_COH_Points'),
    tmpfp = ftell(fid);
    Ncoil = fread(fid, 1,         'uint32');
    N     = fread(fid, 1,         'uint32');
    coils = fread(fid, [7 Ncoil], 'double');

    header.user_block_data{ub}.pnt   = coils(1:3,:)';
    header.user_block_data{ub}.ori   = coils(4:6,:)';
    header.user_block_data{ub}.Ncoil = Ncoil;
    header.user_block_data{ub}.N     = N;
    tmp = fread(fid, (904-288)/8, 'double');
    header.user_block_data{ub}.tmp   = tmp; %FIXME try to find out what these bytes mean
    fseek(fid, tmpfp, 'bof');
  elseif strcmp(type(type>0), 'b_ccp_xfm_block'),
    tmpfp = ftell(fid);
    tmp1 = fread(fid, 1, 'uint32');
    %tmp = fread(fid, [4 4], 'double');
    %tmp = fread(fid, [4 4], 'double');
    %the next part seems to be in little endian format (at least when I tried) 
    tmp = fread(fid, 128, 'uint8');
    tmp = uint8(reshape(tmp, [8 16])');
    xfm = zeros(4,4);
    for k = 1:size(tmp,1)
      xfm(k) = typecast(tmp(k,:), 'double');
      if abs(xfm(k))<1e-10 || abs(xfm(k))>1e10, xfm(k) = typecast(fliplr(tmp(k,:)), 'double');end
    end
    fseek(fid, tmpfp, 'bof'); %FIXME try to find out why this looks so strange
  elseif strcmp(type(type>0), 'b_eeg_elec_locs'),
    %this block contains the digitized coil and electrode positions
    tmpfp   = ftell(fid);
    Npoints = user_space_size./40;
    for k = 1:Npoints
      tmp      = fread(fid, 16, 'uchar');
      %tmplabel = char(tmp(tmp>47 & tmp<128)'); %stick to plain ASCII
      
      % store up until the first space
      tmplabel = char(tmp(1:max(1,(find(tmp==0,1,'first')-1)))'); %stick to plain ASCII
      
      %if strmatch('Coil', tmplabel), 
      %  label{k} = tmplabel(1:5);
      %elseif ismember(tmplabel(1), {'L' 'R' 'C' 'N' 'I'}),
      %  label{k} = tmplabel(1);
      %else
      %  label{k} = '';
      %end
      label{k} = tmplabel;
      tmp      = fread(fid, 3, 'double');
      pnt(k,:) = tmp(:)';
    end

    % post-processing of the labels
    % it seems the following can happen
    %  - a sequence of L R N C I, i.e. the coordinate system defining landmarks
    for k = 1:numel(label)
      firstletter(k) = label{k}(1);
    end
    sel = strfind(firstletter, 'LRNCI');
    if ~isempty(sel)
      label{sel}   = label{sel}(1);
      label{sel+1} = label{sel+1}(1);
      label{sel+2} = label{sel+2}(1);
      label{sel+3} = label{sel+3}(1);
      label{sel+4} = label{sel+4}(1);
    end
    %  - a sequence of coil1...coil5 i.e. the localization coils
    for k = 1:numel(label)
       if strncmpi(label{k},'coil',4)
         label{k} = label{k}(1:5);
       end
    end
    %  - something else: EEG electrodes?
    header.user_block_data{ub}.label = label(:);
    header.user_block_data{ub}.pnt   = pnt;
    fseek(fid, tmpfp, 'bof');
  end
  fseek(fid, user_space_size, 'cof');
end

%channels
for ch = 1:header.config_data.total_chans
  align_file_pointer(fid)
  name                                       = char(fread(fid, 16, 'uchar'))';
  header.config.channel_data(ch).name        = name(name>0);
  %FIXME this is a very dirty fix to get the reading in of continuous headlocalization
  %correct. At the moment, the numbering of the hmt related channels seems to start with 1000
  %which I don't understand, but seems rather nonsensical.
  chan_no                                    = fread(fid, 1, 'uint16=>uint16');
  if chan_no > header.config_data.total_chans,
    
    %FIXME fix the number in header.channel_data as well
    sel     = find([header.channel_data.chan_no]== chan_no);
    if ~isempty(sel),
      chan_no = ch;
      header.channel_data(sel).chan_no    = chan_no;
      header.channel_data(sel).chan_label = header.config.channel_data(ch).name;
    else
      %does not matter
    end
  end
  header.config.channel_data(ch).chan_no     = chan_no;
  header.config.channel_data(ch).type        = fread(fid, 1, 'uint16=>uint16');
  header.config.channel_data(ch).sensor_no   = fread(fid, 1, 'int16=>int16');
  fseek(fid, 2, 'cof');
  header.config.channel_data(ch).gain        = fread(fid, 1, 'float32=>float32');
  header.config.channel_data(ch).units_per_bit = fread(fid, 1, 'float32=>float32');
  yaxis_label                                = char(fread(fid, 16, 'uchar'))';
  header.config.channel_data(ch).yaxis_label = yaxis_label(yaxis_label>0);
  header.config.channel_data(ch).aar_val     = fread(fid, 1, 'double');
  header.config.channel_data(ch).checksum    = fread(fid, 1, 'int32=>int32');
  header.config.channel_data(ch).reserved    = fread(fid, 32, 'uchar=>uchar')';
  fseek(fid, 4, 'cof');

  align_file_pointer(fid)
  header.config.channel_data(ch).device_data.hdr.size     = fread(fid, 1, 'uint32=>uint32');
  header.config.channel_data(ch).device_data.hdr.checksum = fread(fid, 1, 'int32=>int32');
  header.config.channel_data(ch).device_data.hdr.reserved = fread(fid, 32, 'uchar=>uchar')';

  switch header.config.channel_data(ch).type
    case {1, 3}%meg/ref

      header.config.channel_data(ch).device_data.inductance  = fread(fid, 1, 'float32=>float32');
      fseek(fid, 4, 'cof');
      header.config.channel_data(ch).device_data.Xfm         = fread(fid, [4 4], 'double');
      header.config.channel_data(ch).device_data.xform_flag  = fread(fid, 1, 'uint16=>uint16');
      header.config.channel_data(ch).device_data.total_loops = fread(fid, 1, 'uint16=>uint16');
      header.config.channel_data(ch).device_data.reserved    = fread(fid, 32, 'uchar=>uchar')';
      fseek(fid, 4, 'cof');

      for loop = 1:header.config.channel_data(ch).device_data.total_loops
        align_file_pointer(fid)
        header.config.channel_data(ch).device_data.loop_data(loop).position    = fread(fid, 3, 'double');
        header.config.channel_data(ch).device_data.loop_data(loop).direction   = fread(fid, 3, 'double');
        header.config.channel_data(ch).device_data.loop_data(loop).radius      = fread(fid, 1, 'double');
        header.config.channel_data(ch).device_data.loop_data(loop).wire_radius = fread(fid, 1, 'double');
        header.config.channel_data(ch).device_data.loop_data(loop).turns       = fread(fid, 1, 'uint16=>uint16');
        fseek(fid, 2, 'cof');
        header.config.channel_data(ch).device_data.loop_data(loop).checksum    = fread(fid, 1, 'int32=>int32');
        header.config.channel_data(ch).device_data.loop_data(loop).reserved    = fread(fid, 32, 'uchar=>uchar')';
      end
    case 2%eeg
      header.config.channel_data(ch).device_data.impedance       = fread(fid, 1, 'float32=>float32');
      fseek(fid, 4, 'cof');
      header.config.channel_data(ch).device_data.Xfm             = fread(fid, [4 4], 'double');
      header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
    case 4%external
      header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
      header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
      fseek(fid, 4, 'cof');
    case 5%TRIGGER
      header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
      header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
      fseek(fid, 4, 'cof');
    case 6%utility
      header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
      header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
      fseek(fid, 4, 'cof');
    case 7%derived
      header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
      header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
      fseek(fid, 4, 'cof');
    case 8%shorted
      header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
    otherwise
      error('Unknown device type: %d\n', header.config.channel_data(ch).type);
  end
end

fclose(fid);
%end read config file

header.header_data.FileDescriptor = 0; %no obvious field to take this from
header.header_data.Events         = 1;%no obvious field to take this from
header.header_data.EventCodes     = 0;%no obvious field to take this from

if isfield(header, 'channel_data'),
  header.ChannelGain        = double([header.config.channel_data([header.channel_data.chan_no]).gain]');
  header.ChannelUnitsPerBit = double([header.config.channel_data([header.channel_data.chan_no]).units_per_bit]');
  header.Channel            = {header.config.channel_data([header.channel_data.chan_no]).name}';
  header.ChannelType        = double([header.config.channel_data([header.channel_data.chan_no]).type]');
  %header.Channel            = {header.channel_data.chan_label}';
  %header.Channel            = {header.channel_data([header.channel_data.index]+1).chan_label}';
  header.Format             = header.header_data.Format;
  
  % take the EEG labels from the channel_data, and the rest of the labels
  % from the config.channel_data. Some systems have overloaded MEG channel
  % labels, which clash with the labels in the grad-structure. This will
  % lead to problems in forward/inverse modelling. Use the following
  % convention: keep the ones from the config.
  % Some systems have overloaded EEG channel
  % labels, rather than Exxx have a human interpretable form. Use these,
  % to prevent a clash with the elec-structure, if present. This is a
  % bit clunky (because EEG is treated different from MEG), but inherent is
  % inherent in how the header information is organised.
  header.Channel(header.ChannelType==2) = {header.channel_data(header.ChannelType==2).chan_label}';
  
end

function align_file_pointer(fid)
current_position = ftell(fid);
if mod(current_position, 8) ~= 0
    offset = 8 - mod(current_position,8);
    fseek(fid, offset, 'cof');
end
