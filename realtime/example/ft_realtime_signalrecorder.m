function ft_realtime_signalrecorder(cfg)

% FT_REALTIME_SIGNALRECORDER is an example realtime application for recording of data
% that is streaming to the buffer in real-time. It should work both for EEG and MEG.
%
% Use as
%   ft_realtime_signalrecorder(cfg)
% with the following configuration options
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'last')
%   cfg.jumptoeof  = whether to skip to the end of the stream/file at startup (default = 'yes')
%
% The source of the data, i.e. where it comes from, is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% The target for the data, i.e. where it goes to, is configured as
%   cfg.export.dataset    = string with the output file name
%   cfg.export.dataformat = string describing the output file format, see FT_WRITE_DATA
%
% Some notes about skipping data and catching up with the data stream:
%
% cfg.jumptoeof='yes' causes the realtime function to jump to the end
% when the function _starts_. It causes all data acquired prior to
% starting the realtime function to be skipped.
%
% cfg.bufferdata='last' causes the realtime function to jump to the last
% available data while _running_. If the realtime loop is not fast enough,
% it causes some data to be dropped.
%
% If you want to skip all data that was acquired before you start the
% RT function, but don't want to miss any data that was acquired while
% the realtime function is started, then you should use jumptoeof=yes and
% bufferdata='first'. If you want to analyse data from a file, then you
% should use cfg.jumptoeof='no' and cfg.bufferdata='first'.
%
% To stop this realtime function, you will have have to press Ctrl-C.

% Copyright (C) 2012, Robert Oostenveld
%
% $Id$

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0;          end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last
if ~isfield(cfg, 'jumptoeof'),      cfg.jumptoeof = 'no';     end % jump to end of file at initialization

if ~isfield(cfg, 'dataset') && ~isfield(cfg, 'header') && ~isfield(cfg, 'datafile')
  cfg.dataset = 'buffer://localhost:1972';
end

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true, 'retry', true);

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  ft_error('no channels were selected');
end

% make a copy of the header that will be passed to the writing function, update with the channel selection
writehdr          = hdr;
writehdr.nChans   = length(chanindx);
writehdr.label    = writehdr.label(chanindx);
writehdr.chantype = writehdr.chantype(chanindx);
writehdr.chanunit = writehdr.chanunit(chanindx);

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);

if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample = 0;
end
count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true
  
  % determine number of samples available in buffer
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);
  
  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);
  
  if newsamples>=blocksize
    
    % determine the samples to process
    if strcmp(cfg.bufferdata, 'last')
      begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
      endsample  = hdr.nSamples*hdr.nTrials;
    elseif strcmp(cfg.bufferdata, 'first')
      begsample  = prevSample+1;
      endsample  = prevSample+blocksize ;
    else
      ft_error('unsupported value for cfg.bufferdata');
    end
    
    % this allows overlapping data segments
    if overlap && (begsample>overlap)
      begsample = begsample - overlap;
      endsample = endsample - overlap;
    end
    
    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
    
    % read data segment from buffer
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to saving of the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if count==1
      ft_write_data(cfg.export.dataset, dat, 'header', writehdr, 'dataformat', cfg.export.dataformat, 'append', false);
    else
      ft_write_data(cfg.export.dataset, dat, 'header', writehdr, 'dataformat', cfg.export.dataformat, 'append', true);
    end
    
  end % if enough new samples
end % while true
