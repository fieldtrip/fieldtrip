function rt_fileproxy(cfg)

% RT_FILEPROXY reads continuous data from an EEG/MEG file and writes it to
% a FieldTrip buffer. This works for any file format that is supported by
% FieldTrip.
%
% The FieldTrip buffer is a network transparent server that allows the
% acquisition client to stream data to it. An analysis client can connect
% to read the data upon request. Multiple clients can connect simultaneously,
% each analyzing a specific aspect of the data concurrently.
%
% Use as
%   rt_fileproxy(cfg)
% with the following configuration options
%   cfg.minblocksize         = number, in seconds (default = 0)
%   cfg.maxblocksize         = number, in seconds (default = 1)
%   cfg.channel              = cell-array, see CHANNELSELECTION (default = 'all')
%   cfg.jumptoeof            = jump to end of file at initialization (default = 'no')
%   cfg.readevent            = whether or not to copy events (default = 'no')
%   cfg.speed                = relative speed at which data is written (default = inf)
%
% The source of the data is configured as
%   cfg.source.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.source.datafile      = string
%   cfg.source.headerfile    = string
%   cfg.source.eventfile     = string
%   cfg.source.dataformat    = string, default is determined automatic
%   cfg.source.headerformat  = string, default is determined automatic
%   cfg.source.eventformat   = string, default is determined automatic
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: rt_fileproxy.m,v $
% Revision 1.20  2009/07/29 13:32:07  marvger
% source->target
%
% Revision 1.19  2009/07/22 12:12:08  marvger
% also changed cfg.target.eventfile and added defaults to keep our online
% BCI from breaking
%
% Revision 1.18  2009/07/20 09:45:09  roboos
% small change in comments, nothing functional
%
% Revision 1.17  2009/07/20 09:17:11  roboos
% read events from cfg.source.eventfile and not datafile
%
% Revision 1.16  2009/05/01 15:41:21  roboos
% consistent handling of cfg.speed
%
% Revision 1.15  2009/04/23 12:50:49  marvger
% update of the BCI realtime code
%
% Revision 1.14  2009/02/10 12:42:54  marvger
% maxblocksize = 1
%
% Revision 1.13  2009/02/10 10:52:58  marvger
% default blocksize set to 0.1 (100 ms blocks as in the realtime protocol)
%
% Revision 1.12  2009/02/10 10:43:37  marvger
% small changes
%
% Revision 1.11  2009/02/04 09:08:07  roboos
% ensure that the persistent variables related to header caching are cleared
% this is needed when switching the headerformat (from ctf_res4 to ctf_old) while continuing on the same file
%
% Revision 1.10  2009/02/03 20:25:06  marvger
% renamed rt_timer to rt_packettimer
%
% Revision 1.9  2009/02/02 08:15:12  marvger
% placed data reading/writing in a try-catch block
%
% Revision 1.8  2009/01/29 09:14:50  marvger
% experienced problems with read/write event; needs fixing!
%
% Revision 1.7  2009/01/15 12:09:06  marvger
% restored while loop
%
% Revision 1.6  2009/01/15 11:22:24  marvger
% changed reading/writing of events
%
% Revision 1.5  2009/01/14 21:16:52  marvger
% changes related to realtime processing
%
% Revision 1.4  2008/12/01 14:48:57  roboos
% merged in the changes made in Lyon, general cleanup
%
% Revision 1.3  2008/11/14 16:23:41  roboos
% numerous changes to make the rt_xxx functions more similar
%
% Revision 1.2  2008/10/28 14:05:16  roboos
% updated documentation and defaults
%
% Revision 1.1  2008/10/24 08:51:38  roboos
% new implementation
%

% set the defaults
if ~isfield(cfg, 'source'),               cfg.source = [];                                  end
if ~isfield(cfg, 'target'),               cfg.target = [];                                  end
if ~isfield(cfg.source, 'headerformat'),  cfg.source.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.source, 'dataformat'),    cfg.source.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'headerformat'),  cfg.target.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.target, 'dataformat'),    cfg.target.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'datafile'),      cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg, 'minblocksize'),         cfg.minblocksize = 0;                             end % in seconds
if ~isfield(cfg, 'maxblocksize'),         cfg.maxblocksize = 1;                             end % in seconds
if ~isfield(cfg, 'channel'),              cfg.channel = 'all';                              end
if ~isfield(cfg, 'jumptoeof'),            cfg.jumptoeof = 'no';                             end % jump to end of file at initialization
if ~isfield(cfg, 'readevent'),            cfg.readevent = 'no';                             end % capture events?
if ~isfield(cfg, 'speed'),                cfg.speed = inf ;                                 end % inf -> run as fast as possible

% translate dataset into datafile+headerfile
cfg.source = checkconfig(cfg.source, 'dataset2files', 'yes');
cfg.target = checkconfig(cfg.target, 'dataset2files', 'yes');
checkconfig(cfg.source, 'required', {'datafile' 'headerfile'});
checkconfig(cfg.target, 'required', {'datafile' 'headerfile'});

if ~isfield(cfg.source,'eventfile') || isempty(cfg.source.eventfile)
  cfg.source.eventfile = cfg.source.datafile;
end

if ~isfield(cfg.target,'eventfile') || isempty(cfg.target.eventfile)
  cfg.target.eventfile = cfg.target.datafile;
end

% ensure that the persistent variables related to caching are cleared
clear read_header
% read the header for the first time
hdr = read_header(cfg.source.headerfile);
fprintf('updating the header information, %d samples available\n', hdr.nSamples*hdr.nTrials);

% define a subset of channels for reading
cfg.channel = channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
end

minblocksmp = round(cfg.minblocksize*hdr.Fs);
minblocksmp = max(minblocksmp, 1);
maxblocksmp = round(cfg.maxblocksize*hdr.Fs);
count       = 0;

if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample  = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evt = [];
while true

  % determine number of samples available in buffer
  hdr = read_header(cfg.source.headerfile, 'cache', true);

  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

  if newsamples>=minblocksmp

    begsample  = prevSample+1;
    endsample  = prevSample+min(newsamples,maxblocksmp);

    % remember up to where the data was read
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

    % read data segment
    dat = read_data(cfg.source.datafile,'header', hdr, 'dataformat', cfg.source.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    % it only makes sense to read those events associated with the currently processed data
    if strcmp(cfg.readevent,'yes')
      evt = read_event(cfg.source.eventfile, 'header', hdr, 'minsample', begsample, 'maxsample', endsample);
    end

    prevSample  = endsample;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to writing the data to another stream
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if count==1
      % flush the file, write the header and subsequently write the data segment
      write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', false);
      if strcmp(cfg.readevent,'yes')
        write_event(cfg.target.eventfile,evt,'append',false);
      end
    else
      % write the data segment
      write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', true);
      if strcmp(cfg.readevent,'yes')
        write_event(cfg.target.eventfile,evt,'append',true);
      end
    end % if count==1

    % wait for a realistic amount of time
    pause(((endsample-begsample+1)/hdr.Fs)/cfg.speed);

  end % if enough new samples
end % while true
