function ft_realtime_modeegproxy(cfg)

% FT_REALTIME_MODEEGPROXY reads continuous data from a modeeg EEG acquisition system
% through the serial port or through BlueTooth and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_modeegproxy(cfg)
%
% The configuration should contain
%   cfg.filename             = string, name of the serial port (default = '/dev/tty.FireFly-B106-SPP')
%   cfg.feedback             = 'yes' or 'no' (default = 'no')
%   cfg.blocksize            = number, in seconds (default = 0.125)
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2012, Robert Oostenveld
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
if ~isfield(cfg, 'filename'),           cfg.filename = '/dev/tty.FireFly-B106-SPP';       end
if ~isfield(cfg, 'channel'),            cfg.channel = 'all';                              end
if ~isfield(cfg, 'feedback'),           cfg.feedback = 'no';                              end
if ~isfield(cfg, 'blocksize'),          cfg.blocksize = 0.125;                            end
if ~isfield(cfg, 'target'),             cfg.target = [];                                  end
if ~isfield(cfg.target, 'datafile'),    cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg.target, 'dataformat'),  cfg.target.dataformat = [];                       end % default is to use autodetection of the output format

% make a connection to the serial port
fid = fopen(cfg.filename, 'r');
if fid<0
  ft_error('cannot open %s', cfg.filename);
else
  fprintf('opened %s\n', cfg.filename);
  c = onCleanup(@()fclose(fid));
end

rem   = [];   % this will contain the remaining bytes following each conversion
bytes = 0;    % number of bytes read in total
count = 0;    % number of blocks read
blocksize = round(17*256*cfg.blocksize); % 17 bytes per sample, 256 Hz sampling rate
feedback  = strcmp(cfg.feedback, 'yes'); % convert to boolean value

% create a fixed header
hdr.Fs          = 256;
hdr.nChans      = 6;
ndr.nSamples    = inf;
hdr.nSamplesPre = 0;
hdr.label       = cellfun(@num2str, {1, 2, 3, 4, 5, 6}, 'uniformoutput', false);
chanindx        = 1:2;

while (true)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read some data from the serial stream
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  raw = fread(fid, [1 blocksize], 'uint8=>uint8');
  bytes = bytes + numel(raw);
  if feedback
    fprintf('%d bytes read\n', bytes);
  end
  
  raw = cat(2, raw, rem);
  [dat, rem] = decode_modeeg(raw);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % from here onward it is specific to writing the data to another stream
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ~isempty(dat)
    count = count + 1;
    fprintf('writing %d channels, %d samples\n', size(dat,1), size(dat,2));
    if count==1
      % flush the file, write the header and subsequently write the data segment
      ft_write_data(cfg.target.datafile, dat(chanindx,:), 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', false);
    else
      % write the data segment
      ft_write_data(cfg.target.datafile, dat(chanindx,:), 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', true);
    end
  end
  
  
end % while true
