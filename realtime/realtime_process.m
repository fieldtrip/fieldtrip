function realtime_process(cfg)

% REALTIME_PROCESS is used to process the data according to a specified
% trialfun and bcifun. This function can be considered as the general
% workhorse for realtime processing.
%
% Use as
%   realtime_process(cfg)
% with the following configuration options
%   cfg.channel    = cell-array, see CHANNELSELECTION (default = 'all')
%   cfg.trialfun   = string with the trial function
%   cfg.bcifun     = string with the realtime analysis function
%   cfg.readevent  = {'yes' 'no'}; no reading is useful in asychronous
%                       setups where we analyze data continuously (default = 'no')
%   cfg.nexamples  = the number of data segments to be processed (default = inf)
%   cfg.ostream    = the output stream that is used to send a command via
%                     write_event (default = []). E.g., use
%                     'tcp://presentation011:1976' for tcp output to the
%                     presentation machine
%
% The source of the data is configured as
%   cfg.dataset       = string (default = 'buffer://odin:1972')
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% To stop the realtime function, you have to press Ctrl-C

% Copyright (C) 2009, Marcel van Gerven, Robert Oostenveld 
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

  % to ensure clean output
  warning off all;

  % set the default configuration options
  if ~isfield(cfg, 'trialfun'),       cfg.trialfun = @trialfun_realtime;  end
  if ~isfield(cfg, 'bcifun'),         cfg.bcifun = @bcifun_latidx;        end % example neurofeedback bcifun
  if ~isfield(cfg, 'dataset'),        cfg.dataset = 'buffer://odin:1972'; end % default for FCDC CTF system
  if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];                end % default is detected automatically
  if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];              end % default is detected automatically
  if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];               end % default is detected automatically
  if ~isfield(cfg, 'channel'),        cfg.channel = 'all';                end
  if ~isfield(cfg, 'readevent'),      cfg.readevent = 'no';               end % asynchronous/synchronous mode
  if ~isfield(cfg, 'nexamples'),      cfg.nexamples = inf;                end
  if ~isfield(cfg, 'ostream'),        cfg.ostream = [];                   end % output stream
  
  if strcmp(cfg.readevent,'no')
    cfg.event = []; % dummy representation of event in case no triggers are read 
  end
  
  % translate dataset into datafile+headerfile
  cfg = checkconfig(cfg, 'dataset2files', 'yes');
  cfg = checkconfig(cfg, 'required', {'datafile' 'headerfile'});

  % start by reading the header from the realtime buffer
  hdr = read_header(cfg.headerfile, 'cache', true,'retry',true);
  
  % bug fix converts channel names from CTF to fcdc_buffer format
  % FIXME: fcdc_buffer ideally retains the original channel names
  % Note that when fcdc_buffer is called we could always use all channels
  % since channel selection can be assumed to be performed at an earlier
  % stage (on the side of the acquisition machine)
  if isequal(cfg.dataformat,'fcdc_buffer') && ...
      ~any(strcmp(cfg.channel,'all')) && ...
      any(cellfun(@(x)(isnan(str2double(x))),cfg.channel))
    
    % assuming the standard CTF system
    [t,t,t,t,labels] = read_lay('CTF275.lay');
  
    cfg.channel = channelselection(cfg.channel, labels);
  
    x = find(ismember(labels,cfg.channel));
    for j=1:length(x)
      cfg.channel{j} = num2str(x(j),'%.3d'); %  FIXME: this currently fails since fcdc_buffer now represents channel names differently 
    end
  end
  
  % define a subset of channels for reading
  cfg.channel = channelselection(cfg.channel, hdr.label);
  chanindx    = match_str(hdr.label, cfg.channel);
  nchan       = length(chanindx);
  if nchan==0
    error('no channels were selected');
  end
  
  % set the constant fields of data
  data.label     = hdr.label(chanindx);
  data.fsample   = hdr.Fs;
  
  cfg.minsample = 0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % this is the general BCI loop where realtime incoming data is handled
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg.count = 0; % current segment; can be used in bcifun
  while cfg.count < cfg.nexamples
 
    % determine latest header and event information
    if strcmp(cfg.readevent,'yes')
      cfg.event = read_event(cfg.dataset, 'minsample', cfg.minsample+1);  % only consider events that are later than the data processed sofar
    end
    cfg.hdr = read_header(cfg.dataset, 'cache', true); % the trialfun might want to use this, but it is not required
    
    % evaluate the trialfun, note that the trialfun should not re-read the events and header
    %fprintf('evaluating ''%s'' based on %d events\n', func2str(cfg.trialfun), length(cfg.event));
    trl = cfg.trialfun(cfg);

    if size(trl,1) 
      fprintf('processing %d trials using %s\n', size(trl,1),func2str(cfg.bcifun));
    end
    
    for trlidx=1:size(trl,1)

      begsample = trl(trlidx,1);
      endsample = trl(trlidx,2);
      offset    = trl(trlidx,3);
      condition = trl(trlidx,4);

      % remember up to where the data was read
      cfg.minsample  = endsample;
      cfg.count      = cfg.count + 1;
      fprintf('processing segment %d from sample %d to %d, condition = %d\n', cfg.count, begsample, endsample, condition);

      % read data segment from buffer     
      dat = read_data(cfg.datafile, 'dataformat', cfg.dataformat, 'header', cfg.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

      % put the data in a fieldtrip-like raw structure
      data.trial{1}  = dat;
      data.time{1}   = offset2time(offset, hdr.Fs, endsample-begsample+1);
      data.hdr       = hdr;
      data.count     = cfg.count;
      % beg and endsample in data.cfg.trl?
      data.cfg.trl = trl(trlidx,:);

      cmd = cfg.bcifun(cfg, data);
            
      fprintf('%s\n',num2str(cmd));
      
      if ~isempty(cfg.ostream)

          % send command
          cmdevent.type = 'uint';
          cmdevent.offset = [];
          cmdevent.duration = [];
          cmdevent.sample = abs(data.time{1}(1)*data.fsample);
          cmdevent.timestamp = data.time{1}(1);
          cmdevent.value = cmd;
          
          write_event(cfg.ostream,cmdevent);
          
      end     
      
      if cfg.count > cfg.nexamples, break; end

    end % looping over new trials
  end % while true
end
