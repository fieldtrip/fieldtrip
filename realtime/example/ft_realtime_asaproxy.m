function ft_realtime_asaproxy(cfg)

% FT_REALTIME_ASAPROXY reads continuous data from the ASA acquisition system and
% writes it to a FieldTrip buffer. This function uses the NeuroSDK software, which
% can be obtained from ANT.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_asaproxy(cfg)
%
% The configuration should contain
%   cfg.channel              = cell-array, see FT_CHANNELSELECTION (default = 'all')
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2009, Robert Oostenveld
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

% set the defaults
if ~isfield(cfg, 'target'),               cfg.target = [];                                  end
if ~isfield(cfg.target, 'headerformat'),  cfg.target.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.target, 'dataformat'),    cfg.target.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'datafile'),      cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg, 'blocksize'),            cfg.blocksize = 1;                                end % in seconds
% if ~isfield(cfg, 'minblocksize'),         cfg.minblocksize = 0;                             end % in seconds
% if ~isfield(cfg, 'maxblocksize'),         cfg.maxblocksize = 1;                             end % in seconds
% if ~isfield(cfg, 'channel'),              cfg.channel = 'all';                              end

sHost               = 'asa';       % change this to 'asa', 'eemagine', 'tmsi' or 'test' other work with other hosts
bDirectAmp          = true;

% connect to amplifier server
if (strcmpi(sHost, 'asa') || strcmpi(sHost, 'eemagine'))
  bDirectAmp   = false;
  oAcquisition =  actxserver('Eemagine.Acquisition');
end

oDevice = device(sHost);
oDevice = connect(oDevice);

try
  pause(1);
  count = 0;

  while true

    % get new EEG data
    dat = getEEG(oDevice, cfg.blocksize);
    dat = dat';
    
    % FIXME this data is the last N seconds, which means that there will be overlap 
    % between subsequent blocks that are processed by this function 

    count = count + 1;
    nchans = size(dat,1);
    nsamples = size(dat,2);
    fprintf('processing segment %d, %d nchans, %d nsamples\n', count, nchans, nsamples);

    hdr = [];
    hdr.Fs = (nsamples-1)/cfg.blocksize;
    hdr.nChans = nchans;
    hdr.nSamples = nsamples;
    % hdr.nSamplesPre = 0;
    % hdr.label = {};
    % for i=1:nchans
    %   hdr.label{i} = sprintf('%d', i);
    % end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to writing the data to another stream
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if count==1
      % flush the file, write the header and subsequently write the data segment
      ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
    else
      % write the data segment
      ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', true);
    end

    % pause for some time
    pause(cfg.blocksize);
  end
catch
  % disconnect the device
  oDevice = disconnect(oDevice);
  disp(lasterr);
end

