function rt_asaproxy(cfg)

% RT_ASAPROXY reads continuous data from the ASA acquisition system
% and writes it to a FieldTrip buffer. This function uses the
% NeuroSDK software, which can be obtained from ANT.
%
% The FieldTrip buffer is a network transparent server that allows the
% acquisition client to stream data to it. An analysis client can connect
% to read the data upon request. Multiple clients can connect simultaneously,
% each analyzing a specific aspect of the data concurrently.
%
% Use as
%   rt_asaproxy(cfg)
%
% The configuration should contain
%   cfg.channel              = cell-array, see CHANNELSELECTION (default = 'all')
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: rt_asaproxy.m,v $
% Revision 1.6  2009/06/17 13:50:20  roboos
% wait for the specified amount of time instead of !))ms
%
% Revision 1.5  2009/05/01 08:04:19  roboos
% updated documentation
%
% Revision 1.4  2009/01/21 20:59:06  roboos
% added cfg.blocksize, fixed sampling frequency
%
% Revision 1.3  2009/01/20 16:12:18  roboos
% added correct sampling rate to header
% fixed display error
%
% Revision 1.2  2009/01/20 15:43:49  roboos
% first proper implementation
%
% Revision 1.1  2009/01/20 13:14:13  roboos
% added two (still empty) proxy functions for ASA and Plexon
%

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
      write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
    else
      % write the data segment
      write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', true);
    end

    % pause for some time
    pause(cfg.blocksize);
  end
catch
  % disconnect the device
  oDevice = disconnect(oDevice);
  disp(lasterr);
end

