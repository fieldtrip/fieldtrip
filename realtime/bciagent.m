function [cfg] = bciagent(cfg)

% BCIAGENT starts the core processing loop of a brain-computer interface (bci)
% which handles the reading and writing of events and data around calls to the
% real-time application function that us supplied by the user.
%
% Use as
%   bciagent(cfg)
%
% The input configuration cfg is a structure whose fields specify the
% sources and targets of the flow of information (events and data matrices)
% in the bci system, as well as user defined and specified processing that
% should be applied to the data. These fields include:
%   cfg.source.datafile  = (string) name of stream for incoming data (required)
%   cfg.source.eventfile = (string) name of stream for incoming (e.g. trigger) events
%   cfg.target.datafile  = (string) name of stream for outgoing data (optional)
%   cfg.target.eventfile = (string) name of stream for outgoing (e.g. control) events
%   cfg.trialfun         = (string) name of the trial function, see below
%   cfg.bcifun           = (string) name of the bci function, see below
%   cfg.bci              = configuration structure, for storing the parameters of cfg.bcifun
%   cfg.channel          = determines which channels to read from the data source
%   cfg.continuous       = 'yes' or 'no' whether the file contains continuous data (default 
%                          is determined automatic)
%   cfg.delay            = specifies the time to wait at the end of each iteration (default = 0)
%   cfg.delaymethod      = 'fixed', 'adapt' or 'smart' specifies how the
%                          delay is implemented (default = 'fixed')
%   cfg.bufferdata       = 'last', 'all', 'each' and 'first' specifies how and which
%                          data should be passed to cfg.bcifun (default = 'last')
%   cfg.readheader       = 'once' or 'each' (default = 'once')
%
% The function specified in cfg.trialfun detects the data segments that
% have to be processed. It should have the interface
%    [trl, event, filter] = trialfun(cfg)
% where the cfg passed to trialfun corresponds to cfg of the main
% application. Trialfuns are explained in more detail in DEFINETRIAL. The
% cfg passed to trialfun will contain the field cfg.event, i.e. the
% trialfun is not required or supposed to read the events (in contrast
% to trialfuns used in DEFINETRIAL). Furthermore, the cfg passed to the
% trialfun will contain the field cfg.fsample and the field cfg.header,
% i.e. the trialfun is also not required/supposed to read the header
% (in contrast to trialfuns used in DEFINETRIAL).
%
% The function specified in cfg.bcifun processes the data segments and
% should have the interface
%    [event, data] = bcifun(cfg, data, event)
% where the cfg passed to bcifun corresponds to cfg.bci of the main
% application. The input data is a single preprocessed segment of data in
% a format that corresponds with the output of PREPROCESSING.
%
% See also READ_EVENT, WRITE_EVENT, READ_DATA, WRITE_DATA, PREPROCESSING, CHANNELSELECTION

% Undocumented options
%  cfg.source.dataset
%  cfg.source.headerfile
%  cfg.source.dataformat
%  cfg.source.headerformat
%  cfg.source.eventformat
%  cfg.target.dataformat
%  cfg.target.headerformat
%  cfg.target.eventformat
%  cfg.filter
%  TODO implement cfg.target.dataset,  see dataset2files
%  TODO get rid of cfg.preproc, instead the bcifun should call preprocessing

% Copyright (C) 2007-2008, Christian Hesse & Robert Oostenveld
% F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% $Log: bciagent.m,v $
% Revision 1.2  2008/12/16 20:13:53  roboos
% only change in some comments, no functional change
%
% Revision 1.1  2008/11/13 22:32:11  roboos
% moved to realtime directory
%
% Revision 1.27  2008/10/14 10:17:01  sashae
% replaced call to dataset2files with checkconfig
%
% Revision 1.26  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.25  2008/09/22 19:41:43  roboos
% updated documentation
%
% Revision 1.24  2008/01/30 10:41:32  roboos
% moved appendeevnt to seperate function/file, shared with other fileio stuff
%
% Revision 1.23  2008/01/27 20:49:25  roboos
% clear sdata if isempty(trl)
% removed call to preproc and removed rdata
% use numel instead of min(size(..))
% added subfunction appendevent and cleaned up the relevant section of code
%
% Revision 1.22  2008/01/15 10:17:26  roboos
% ensure that filter is not updated in case trialfun fails with an error
%
% Revision 1.21  2008/01/14 20:08:37  roboos
% many small stylistic changes, no functional change
%
% Revision 1.20  2008/01/10 21:26:46  roboos
% added timing information for debugging (trl and bci loop)
% added cfg.readheader, default is once (similar as before) but can also be each for each bci loop
% added hdr to the tmpcfg that is passed to the trialfun
% changed an if-elseif-elseif ladder into a swith-case-case ladder
% some formatting changes
%
% Revision 1.19  2007/12/20 19:03:21  roboos
% quite some changes, here are the most important ones:
% - reimplemented the buffering of data and bcifunskip, it did not work reliably for 'each'
% - removed the auto-filter update for timestamps (was not used)
% - changed the filter update, the trialfun can now return a third output argument which should contain the filter to be used the next time, the default is still the same as before
%
% Revision 1.18  2007/10/29 10:54:15  chrhes
% added some flexibility w.r.t. calculation of the loop delay using the new
% option cfg.delaymethod. changed the functionality associated with the data
% buffering option.
%
% Revision 1.17  2007/10/16 12:36:41  roboos
% added event filter update for mintimestamp
% made start with implementing bufferdata=first/last/all/each
%
% Revision 1.16  2007/09/22 13:41:29  chrhes
% moved the processing of special source events to the end of the bci loop so
% that the BCIFUN can also act on these (e.g., save its current state to file
% on a "stop" event). Moreover, source events of type "stop" are sent to the
% target event queue, thereby propagating throughout the whole bci system.
%
% Revision 1.15  2007/09/15 14:32:50  chrhes
% simplified the code for reading the source events and updating the event
% filter settings so that this is now preceeds the call to the TRIALFUN which
% does additional parsing of the source events
%
% Revision 1.14  2007/09/15 12:04:12  chrhes
% added event filter updates immediately after event reading in the trialfun
% case BEFORE sevent is generated, since the filter settings don't generally
% refer to sevents generated by the trialfun
%
% Revision 1.13  2007/08/22 22:49:05  chrhes
% changed the cfg.filter.minsample update from using endsample to begsample
% as a temporary measure (the entire approach to updating of the event filter
% after each iteration will have to be revised sooner rather than later)
%
% Revision 1.12  2007/08/21 17:13:53  chrhes
% implemented option cfg.bufferdata which, when enabled, causes all avaiable
% trials to be buffered in sdata before bcifun is called; made a few cosmetic
% changes here and there
%
% Revision 1.11  2007/08/21 14:16:35  chrhes
% fixed some typos in documentation
%
% Revision 1.10  2007/08/21 14:12:35  chrhes
% implemented option cfg.source.dataset
% additional checks on function handles for cfg.trialfun and cfg.bcifun
% reading of source events also works without specifying a cfg.trialfun
% the loop does not require source data or source events to run
% response to sevent.type = 'debug' brings up the command prompt
%
% Revision 1.9  2007/08/16 14:25:21  chrhes
% put the writing of target events and target data in try-catches and added
% functionality for appending these after calls to bcifun prior to the next
% writing attempt
%
% Revision 1.8  2007/08/15 15:33:50  chrhes
% updated some documentation
%
% Revision 1.7  2007/08/15 07:38:59  chrhes
% made the display of error messages in a few try-catch statements dependent
% on the global feedback ('fb') variable
%
% Revision 1.6  2007/08/01 16:24:50  roboos
% added fsample to tmpcfg for trialfun
% changed format of cfg.datasource, now souce.datafile/headerfile/eventfile etc.
% idem for other source/target cfgs
% give error if no channels selected (otherwise real all channels by default)
% specify headerformat when reading hdr
% some minor cleanup and comment changed
%
% Revision 1.5  2007/07/29 09:06:24  roboos
% added default for cfg.channel
% added and implemented cfg.continuous (for checking boundary when reading)
% minor changes to handling of some cfg default settings
% fixed bug in defaults for event filters, should be part of cfg
% instead of processing only the first trial from trialfun, loop over all trials
% updated documentation, explain trialfun and bcifun
%
% Revision 1.4  2007/07/27 12:14:47  roboos
% set checkboundary=0 for reading
% filter inside the read_event, and do not call seperate filter_event function
% pass the header to teh read_even function (prevents reading multiple times)
% some cosmetic changes
%
% Revision 1.3  2007/06/13 13:31:36  roboos
% convert raw data matrix into data structure
% added some debugging fprintf's
%
% Revision 1.2  2007/06/11 15:20:29  chrhes
% implemented changes which ensure that the bcifun receives data in a format
% resembling the standard raw data structure as produced by PREPROCESSING
%
% Revision 1.1  2007/06/06 12:38:18  roboos
% renamed bci functions
%
% Revision 1.4  2007/06/06 12:36:03  roboos
% changed the event filtering, implemented the use of a trialfun to determine
% begin and end sample
%
% Revision 1.3  2007/06/06 07:15:54  roboos
% added some serious filtering in external function
% changed data selection prior to reading
%
% Revision 1.2  2007/06/05 13:54:53  roboos
% changed some aspects of the logic of the loop
% implemented writing of target data
% added (optional) cfg parameters for source and target formats (only needed
% since filetype is not yet up-to-date) some coding style changes
%
% Revision 1.1  2007/06/04 19:31:15  chrhes
% initial version added to CVS repository needs further debugging and testing
%

fieldtripdefs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare and initialize global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global fb
if isempty(fb) || ~islogical(fb), fb = true; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin~=1), error('incorrect number of input arguments'); end
% cfg must be a structure
if ~isstruct(cfg), error('argument cfg must be a structure'); end
% set default configuration options
if ~isfield(cfg, 'source'),              cfg.source = [];                   end
if ~isfield(cfg, 'target'),              cfg.target = [];                   end
if ~isfield(cfg.source, 'dataset'),      cfg.source.dataset = [];           end
if ~isfield(cfg.source, 'datafile'),     cfg.source.datafile = [];          end
if ~isfield(cfg.source, 'dataformat'),   cfg.source.dataformat = [];        end
if ~isfield(cfg.source, 'headerfile'),   cfg.source.headerfile = [];        end
if ~isfield(cfg.source, 'headerformat'), cfg.source.headerformat = [];      end
if ~isfield(cfg.source, 'eventfile'),    cfg.source.eventfile = [];         end
if ~isfield(cfg.source, 'eventformat'),  cfg.source.eventformat = [];       end
if ~isfield(cfg.target, 'datafile'),     cfg.target.datafile = [];          end
if ~isfield(cfg.target, 'dataformat'),   cfg.target.dataformat = [];        end
if ~isfield(cfg.target, 'eventfile'),    cfg.target.eventfile = [];         end
if ~isfield(cfg.target, 'eventformat'),  cfg.target.eventformat = [];       end
if ~isfield(cfg, 'trialfun'),            cfg.trialfun = [];                 end
if ~isfield(cfg, 'preproc'),             cfg.preproc = [];                  end
if ~isfield(cfg, 'bci'),                 cfg.bci = [];                      end
if ~isfield(cfg, 'bcifun'),              cfg.bcifun = [];                   end
if ~isfield(cfg, 'delay'),               cfg.delay = 0;                     end
if ~isfield(cfg, 'delaymethod'),         cfg.delaymethod = 'fixed';         end
if ~isfield(cfg, 'channel'),             cfg.channel = 'all';               end
if ~isfield(cfg, 'bufferdata'),          cfg.bufferdata = 'last';           end
if ~isfield(cfg, 'readheader'),          cfg.readheader = 'once';           end

% set the default filters for the incoming events
if ~isfield(cfg, 'filter')
  cfg.filter = [];
  cfg.filter.minsample = 1;
  cfg.filter.type = {'newdata','stop','debug','trigger'};
end

% check whether the source data has been specified as a dataset
if isstr(cfg.source.dataset)
  cfg.source = checkconfig(cfg.source, 'dataset2files', {'yes'});
end

if ~isempty(cfg.source.datafile) && ~isstr(cfg.source.datafile)
  error('field cfg.source.datafile must be a string, empty or omitted');
end

if ~isempty(cfg.source.eventfile) && ~isstr(cfg.source.eventfile)
  warning('field cfg.source.eventfile must be a string, empty or omitted');
end

if ~isempty(cfg.target.datafile) && ~isstr(cfg.target.datafile)
  error('field cfg.target.datafile must be a string, empty or omitted');
end

if ~isempty(cfg.target.eventfile) && ~isstr(cfg.target.eventfile)
  warning('field cfg.target.eventfile must be a string, empty or omitted');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set/check fields to do with channel selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.source.headerfile)
  hdr         = read_header(cfg.source.headerfile, 'headerformat', cfg.source.headerformat);
  cfg.channel = channelselection(cfg.channel, hdr.label);
  chanindx    = match_str(hdr.label, cfg.channel);
else
  hdr      = [];
  chanindx = [];
end

if isempty(chanindx)
  if fb, disp('no channels selected'); end;
end

if fb
  fprintf('nchan      = %d\n', length(chanindx));
end

% this option relates to reading over trial boundaries in a pseudo-continuous dataset
if ~isfield(cfg, 'continuous') && isstruct(hdr)
  if isfield(cfg, 'datatype') && strcmp(cfg.datatype, 'continuous')
    cfg.continuous = 'yes';
  elseif hdr.nTrials==1
    cfg.continuous = 'yes';
  else
    cfg.continuous = 'no';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set/check function handles for cfg.trialfun and cfg.bcifun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a function handle for cfg.trialfun
if isstr(cfg.trialfun)
  trialfun = str2func(cfg.trialfun);
  if ~isa(trialfun,'function_handle')
    error('could not generate a valid function handle from field cfg.trialfun');
  end
else
  trialfun = [];
  if fb,
    disp('field cfg.trialfun not specified');
  end
end

% generate a function handle for cfg.bcifun
if isstr(cfg.bcifun)
  bcifun = str2func(cfg.bcifun);
  if ~isa(bcifun,'function_handle')
    error('could not generate a valid function handle from field cfg.bcifun');
  end
else
  bcifun = [];
  if fb,
    disp('field cfg.bcifun not specified');
  end
end

% NOTE: should the convention be to initialize the BCIFUN with a call???
% cfg.bci = bcifun(cfg.bci);

% flag indicating whether to skip calls to bcifun etc during data buffering
% is always initialized to FALSE
bcifunskip = false;

% check validity of cfg.delay
if ~isscalar(cfg.delay) || ~isreal(cfg.delay) || (cfg.delay<0)
  error('invalid specification of cfg.delay');
end

% check cfg.delaymethod and change into a number to speed up processing
switch cfg.delaymethod
  case 'fixed', delaymethod = 0;
  case 'adapt', delaymethod = 1;
  case 'smart', delaymethod = 2;
  otherwise
    error('invalid setting for cfg.delaymethod');
end % switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start the loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = 0; % used for controlling the delay
t2 = 0; % used for controlling the delay
dt = 0; % used for controlling the delay
tevent = [];
tdata  = [];

again  = 1;
while again
  % NOTE: there is currently no way out of this loop other than "ctrl + C" or
  % through an externtal event with event.type = 'stop'

  % this is for the timing
  t_bciloop = clock;
  
  % check/record time at start of iteration
  if  (cfg.delay>0) && (delaymethod~=0), t1 = clock; end;

  % show progress on screen
  if fb, fprintf('--------------------------------------------------\n'); end;
  if fb, fprintf('running bci loop\n'); end;

  switch cfg.readheader
  case 'once'
    % reuse the header that was read prior to starting the bci loop
  case 'each'
   % rereading it is usefull if hdr.nSamples increases when chasing the end of file, but it can also slow down the bci loop
    hdr = read_header(cfg.source.headerfile, 'headerformat', cfg.source.headerformat);
  otherwise
    error('incorrect specification of cfg.readheader');
  end

  % read (and filter) the source events
  % apply the latest event filter settings
  try
    filter = cfg2keyval(cfg.filter);
    sevent = read_event(cfg.source.eventfile, 'header', hdr, 'eventformat', cfg.source.eventformat, filter{:});
  catch
    if fb, disp(lasterr); end;
    sevent = [];
  end

  % use call to trial function to determine the data segments of interest for further processing
  if isa(trialfun,'function_handle')
    try
      tmpcfg         = cfg;
      tmpcfg.fsample = hdr.Fs;  % it is permitted to use this inside the trialfun
      tmpcfg.hdr     = hdr;     % it is permitted to use this inside the trialfun
      tmpcfg.event   = sevent;  % it is permitted to use this inside the trialfun
      % call trial function (this should only change/remove trial related events)
      switch nargout(trialfun)
      case 1
        [trl] = trialfun(tmpcfg);
      case 2
        [trl, sevent] = trialfun(tmpcfg);
      case 3
        [trl, sevent, filter] = trialfun(tmpcfg);
      end
    catch
      if fb, disp(lasterr); end;
      trl = [];
      filter = [];
    end
  else
    trl = [];
    filter = [];
  end

  if nargout(trialfun)<3
    % update the event filter settings using the default rules (this is still somewhat ad hoc)
    if isstruct(sevent) && numel(sevent)>0 && isfield(cfg.filter,'minsample')
      tmp1 = [sevent(:).sample];
      tmp1 = tmp1(:);
      tmp0 = find(isfinite(tmp1));
      if ~isempty(tmp0)
        cfg.filter.minsample = max([cfg.filter.minsample, max(tmp1(tmp0))]) + 1;
      end
    end
  elseif ~isempty(filter)
    % use the updated filters as returned by the trialfun
    cfg.filter = filter;
  end

  % determine how to continue if multiple trials were defined
  if ~isempty(trl)
    switch cfg.bufferdata
    case 'last'
      trl = trl(end,:);
      bcifunskip = false;
    case 'first'
      trl = trl(1,:);
      bcifunskip = false;
    case 'all'
      % put a temporary lock on calls to bcifun while data is being buffered
      bcifunskip = true;
    case 'each'
      % call bcifun for each trial directly after the data has been read
      bcifunskip = false;
    otherwise
      error('incorrect specification of cfg.bufferdata');
    end
  end

  if fb
    fprintf('nevent     = %d\n', length(sevent));
    fprintf('ntrl       = %d\n', size(trl,1));
  end

  % this is for debugging only
  % if ~isempty(trl), toc, tic, end;

  if isempty(trl)
    % this is to call the bcifun once when trl is empty, it will then only get the events as input
    trlsel = 0;
  else
    % execute the loop below for each trial 
    trlsel=1:size(trl,1);
  end
  
  for trllop=trlsel

    % this is for the timing
    t_trlloop = clock;

    if trllop==0
      % no trials were defined, call the bcifun with only the events
      bcifunskip = false; % make sure bcifun gets called
      sdata = [];

    else
      % one or multiple trials were defined, read the data for each of them
      begsample = trl(trllop,1);
      endsample = trl(trllop,2);
      offset    = trl(trllop,3);

      if fb
        fprintf('trial      = %d\n', trllop);
        fprintf('begsample  = %d\n', begsample);
        fprintf('endsample  = %d\n', endsample);
        fprintf('offset     = %d\n', offset);
      end

      if trllop==size(trl,1)
        % the last trial will be read on this iteration
        bcifunskip = false; % make sure bcifun gets called
      end

      try
        if strcmp(cfg.bufferdata,'each') || trllop==1
          % read the data from the stream and
          % make a new raw data structure containing a single trial
          sdata          = [];
          sdata.trial{1} = read_data(cfg.source.datafile, 'header', hdr, 'dataformat', cfg.source.dataformat, 'chanindx', chanindx, 'begsample', begsample, 'endsample', endsample, 'checkboundary', strcmp(cfg.continuous, 'no'));
          sdata.time{1}  = offset2time(offset, hdr.Fs, endsample-begsample+1);
          sdata.cfg.trl  = trl(trllop,:);       % for consistency with preprocessing and e.g. redefinetrial, the bcifun is allowed to know where this trial originated from
          sdata.label    = hdr.label(chanindx); % add the labels of the selected channels
          sdata.hdr      = hdr;                 % add the header details
          sdata.fsample  = hdr.Fs;              % add the sampling frequency
          if isfield(hdr, 'grad')
            sdata.grad   = hdr.grad;            % add the gradiometer system in head coordinates
          end
        elseif strcmp(cfg.bufferdata,'all')
          % read the data from the stream and
          % add the new trial to the existing raw data structure
          sdata.trial{trllop} = read_data(cfg.source.datafile, 'header', hdr, 'dataformat', cfg.source.dataformat, 'chanindx', chanindx, 'begsample', begsample, 'endsample', endsample, 'checkboundary', strcmp(cfg.continuous, 'no'));
          sdata.time{trllop}  = offset2time(offset, hdr.Fs, endsample-begsample+1);
          sdata.cfg.trl       = trl;            % for consistency with preprocessing and e.g. redefinetrial, the bcifun is allowed to know where the trials originated from
        end
      catch
        if fb, disp(lasterr); end;
        sdata = [];
      end % try reading and preprocessing

    end % if trllop

    % call bcifun, write_event and write_data once all data is available
    if ~bcifunskip

      % do something with sdata and/or sevent, this produces a target event and optionally more data
      [tempevent, tempdata] = bcifun(cfg.bci, sdata, sevent);

      % append tevent
      if isempty(tevent)
        tevent = tempevent;
      elseif ~isempty(tempevent)
        tevent = appendevent(tevent,tempevent);
      end
      clear tempevent;

      % append tdata
      if isempty(tdata)
        tdata = tempdata;
      elseif ~isempty(tempdata)
        tdata = appenddata([],tdata,tempdata);
      end
      clear tempdata;

      % write outgoing event
      if isempty(tevent) && ~isempty(cfg.target.eventfile)
        try
          write_event(cfg.target.eventfile, tevent, 'eventformat', cfg.target.eventformat);
          tevent = [];
        catch
          if fb, disp(lasterr); end;
          % current tevent will be appended on next iteration
        end
      end

      % write outgoing data
      if ~isempty(tdata) && ~isempty(cfg.target.datafile)
        try
          write_data(cfg.target.datafile, tdata, 'header', hdr, 'chanindx', chanindx, 'dataformat', cfg.target.dataformat);
          tdata = [];
        catch
          if fb, disp(lasterr); end;
          % current tdata will be appended on next iteration
        end
      end

    end % if ~bcifunskip

    if fb
      fprintf('time elapsed in trl loop = %f seconds\n', etime(clock, t_trlloop));
    end

  end % for trllop

  % process any "special" events
  if isstruct(sevent) && isfield(sevent,'type') 
    if any(strcmp({sevent.type}, 'stop'))
      % send the "stop" event to the target event queue
      if ~isempty(cfg.target.eventfile)
        ii = find(strcmp({sevent.type}, 'stop'));
        try
          write_event(cfg.target.eventfile, sevent(ii), 'eventformat', cfg.target.eventformat);
        catch
          if fb, disp(lasterr); end;
        end
      end
      % break out of the while-loop
      again = 0;
      continue;
    end
    if any(strcmp({sevent.type}, 'debug'))
      % FIXME also send the "debug" event to the target event queue?
      % go to the command prompt
      keyboard
    end
  end

  % wait for a specified delay period
  if (cfg.delay>0)
    switch delaymethod
    case 0 % fixed
      pause(cfg.delay);
    case 1 % adapt
      t2 = clock;
      dt = cfg.delay - etime(t2,t1);
      if (dt>=0.001)
        pause(dt);
      end
    case 2 % smart
      t2 = clock;
      dt = dt + (cfg.delay - etime(t2,t1));
      if (dt>=0.001)
        pause(dt);
        dt = 0;
      end
    otherwise
      error('invalid setting for cfg.delaymethod');
    end % switch delaymethod
  end % if cfg.delay

  if fb
    fprintf('time elapsed in bci loop = %f seconds\n', etime(clock, t_bciloop));
  end
  
end % while again

disp('bci loop stopped');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

