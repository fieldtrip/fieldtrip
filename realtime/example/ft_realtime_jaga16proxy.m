function ft_realtime_jaga16proxy(cfg)

% FT_REALTIME_JAGA16PROXY reads continuous EEG data from a Jinga-Hi JAGA16 system
% through the UDP network interface and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_jaga16proxy(cfg)
%
% The configuration should contain
%   cfg.port                 = number, UDP port to listen on (default = 55000)
%   cfg.channel              = cell-array with channel names, see FT_CHANNELSELECTION
%   cfg.blocksize            = number, in seconds (default = 0.5)
%   cfg.decimate             = integer number (default = 1)
%   cfg.calibration          = number, in uV per bit (default = 1)
%   cfg.feedback             = 'yes' or 'no' (default = 'no')
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2015-2016, Robert Oostenveld
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

cfg = ft_checkconfig(cfg);

% set the defaults
cfg.port               = ft_getopt(cfg, 'port', 55000);
cfg.channel            = ft_getopt(cfg, 'channel', 'all');
cfg.precision          = ft_getopt(cfg, 'precision', 'single');
cfg.demean             = ft_getopt(cfg, 'demean', 'yes');
cfg.decimate           = ft_getopt(cfg, 'decimate', 1);
cfg.blocksize          = ft_getopt(cfg, 'blocksize', 0.5);    % seconds, this is only approximate
cfg.timeout            = ft_getopt(cfg, 'timeout', 5);        % seconds
cfg.calibration        = ft_getopt(cfg, 'calibration', 0.2);  % uv/bit
cfg.feedback           = ft_getopt(cfg, 'feedback', 'yes');
cfg.target             = ft_getopt(cfg, 'target', []);
cfg.target.datafile    = ft_getopt(cfg.target, 'datafile', 'buffer://localhost:1972');
cfg.target.dataformat  = ft_getopt(cfg.target, 'dataformat', []); % default is to use autodetection of the output format

% in principle this could be implemented using the MATLAB UDP object,
% but that requires the Mathwoks Instrument Control Toolbox

% this requires an external toolbox for the TCP communication
ft_hastoolbox('tcp_udp_ip', 1);

% ensure that the persistent variables inside these functions are reinitialized
clear pnet
clear buffer

% start listening on the UDP port
con = pnet('udpsocket', cfg.port);

if (con<0)
  error('unable to establish connection with host');
else
  fprintf('listening on UDP port %d\n', cfg.port);
end

% I am not sure whether this has any effect
pnet(con, 'setreadtimeout', cfg.timeout)

% part of this is hard coded for the Jinga-Hi JAGA16
hdr = [];
hdr.Fs          = 1000;       % sampling frequency, is updated further down
hdr.nChans      = 16;         % number of channels, is updated further down
hdr.nSamples    = 0;          % number of samples
hdr.nSamplesPre = 0;          % number of pre-trigger samples, not used
hdr.nTrials     = 1;          % number of trials, represent as a continuous strea,
hdr.label       = cell(16,1); % Nx1 cell-array with the label of each channel
for i=1:hdr.nChans
  hdr.label{i} = sprintf('%d', i);
end

% determine the selection of channels to be transmitted
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);

% update the header according to the selection
hdr.nChans = numel(chanindx);
hdr.label  = hdr.label(chanindx);

count = 0;
dat   = cast([], cfg.precision);

while (true)
  % read a packet
  
  buf = [];
  while isempty(buf)
    siz = pnet(con,'readpacket');
    buf = pnet(con, 'read');
    if isempty(buf)
      pause(43/1023);
    end
  end
  packet = jaga16_packet(buf, false);
  
  if strcmp(cfg.feedback, 'yes')
    fprintf('received packet\n');
  end
  
  % do some sanity checks on the hard-coded parameters
  assert(packet.nchan==16);
  assert(packet.nbit==16);
  % update the sampling frequency
  hdr.Fs = packet.fsample;
  
  % get the data and do some minimal conditioning
  buf = double(packet.dat(chanindx,:));
  if strcmp(cfg.demean, 'yes')
    % the data is uint16 and hovers around 32768
    buf = buf - double(intmax('uint16')/2);
  end
  if cfg.calibration~=1
    buf = cfg.calibration * buf;
  end
  buf = cast(buf, cfg.precision);
  
  % concatenate the data fragments until we have a complete block
  dat = cat(2, dat, buf);
  if size(dat,2) < hdr.Fs*cfg.blocksize
    continue
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % from here onward it is specific to writing the data to another stream
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % only write every Nth sample
  % this is a kludge, as it drops occasional samples and ignores aliassing
  if cfg.decimate~=1
    hdr.Fs = hdr.Fs/cfg.decimate;     % update the sampling frequency
    dat = dat(:,1:cfg.decimate:end);
  end
  
  count = count + 1;
  if strcmp(cfg.feedback, 'yes')
    fprintf('writing %d channels, %d samples\n', size(dat,1), size(dat,2));
  end
  if count==1
    % flush the file, write the header and subsequently write the data segment
    ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', false);
  else
    % write the data segment
    ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', true);
  end
  
  if cfg.decimate~=1
    % revert the updated sampling frequency
    hdr.Fs = hdr.Fs*cfg.decimate;
  end
  
  % start with an empty data block for the next iteration
  dat = cast([], cfg.precision);
  
end % while true

% FIXME close the connection, this should be handled in some try-catch statement
pnet(con, 'close');
