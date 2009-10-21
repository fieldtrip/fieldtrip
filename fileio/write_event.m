function write_event(filename, event, varargin)

% WRITE_EVENT writes an event structure to a file, a message daemon
% listening on a network socked, or to another computer connected through
% the serial port.
%
% Use as
%   write_event(filename, event, ...)
%
% The first argument is a string containing the filename. The second
% argument is a structure with the event. Multiple events can be
% represented as a structure array.
%
% Events are represented as
%   event.type      string
%   event.sample    expressed in samples, the first sample of a recording is 1
%   event.value     number or string
%   event.offset    expressed in samples
%   event.duration  expressed in samples
%   event.timestamp expressed in timestamp units, which vary over systems (optional)
%
% Events can also be written to special communication streams
% by specifying the target as URI instead of a filename. Supported are
%   buffer://<host>:<port>
%   fifo://<filename>
%   tcp://<host>:<port>
%   udp://<host>:<port>
%   mysql://<user>:<password>@<host>:<port>
%   rfb://<password>@<host>:<port>
%   serial:<port>?key1=value1&key2=value2&...
%   rfb://<password>@<host>:<port>
%
% See also READ_HEADER, READ_DATA, READ_EVENT, WRITE_DATA

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: write_event.m,v $
% Revision 1.35  2009/05/22 09:02:29  marvger
% changed tcp handling
%
% Revision 1.34  2009/04/28 08:33:05  marvger
% small changes
%
% Revision 1.33  2009/01/22 15:31:59  marvger
% updated catch handling
%
% Revision 1.32  2009/01/21 11:34:44  marvger
% update in hostname detection for fcdc_buffer
%
% Revision 1.31  2009/01/20 08:56:51  marvger
% fixed catch me bug (illegal syntax in older matlab versions)
%
% Revision 1.30  2009/01/16 11:38:51  marvger
% update tcp/udp
%
% Revision 1.29  2009/01/14 21:16:52  marvger
% changes related to realtime processing
%
% Revision 1.28  2008/12/19 14:39:25  marvger
% added support for udp, tcp and fifo
%
% Revision 1.27  2008/06/19 19:31:44  roboos
% made fcdc_buffer more robust
%
% Revision 1.26  2008/06/18 06:22:24  roboos
% added support for fcdc_buffer
%
% Revision 1.25  2008/01/30 10:40:54  roboos
% moved catevent to seperate function and renamed to appendevent
%
% Revision 1.24  2007/12/12 11:29:07  roboos
% chedk for presence of timestamp prior to trying to concatenate events
%
% Revision 1.23  2007/11/07 11:03:58  roboos
% added line of documentation
%
% Revision 1.22  2007/11/07 10:49:07  roboos
% cleaned up the reading and writing from/to mysql database, using db_xxx helper functions (see mysql directory)
%
% Revision 1.21  2007/11/05 17:04:03  roboos
% moved insert_query to seperate function
% ved declarationof persistent and global variables to the beginning of the function
% some cosmetic changes
%
% Revision 1.20  2007/10/30 20:47:04  roboos
% added support for remote frame buffer (rfb), i.e. VNC server
%
% Revision 1.19  2007/10/16 12:34:43  roboos
% use recursion to write to multiple event sources
% implemented fcdc_global
%
% Revision 1.18  2007/09/15 14:36:26  chrhes
% changed the default mode for writing to a serial port to synchronous
%
% Revision 1.17  2007/08/16 14:28:59  chrhes
% removed the try-catches introduced in previous revision: the relevant error
% is now dealt with elsewhere
%
% Revision 1.16  2007/08/16 14:03:20  chrhes
% put some try-catches around the save command in the default option
%
% Revision 1.15  2007/06/19 11:11:28  chrhes
% changed the implementation of how the target serial port is located if it
% exists
%
% Revision 1.14  2007/06/19 10:11:37  chrhes
% restricted the functionality of serial port writing such that only the
% content of the field event.value is written as a character array
%
% Revision 1.13  2007/06/13 14:46:35  roboos
% removed type/subtype, added type/value/sample for mysql
%
% Revision 1.12  2007/06/13 08:06:45  roboos
% updated help
%
% Revision 1.11  2007/06/12 19:37:28  roboos
% added support for mysql port specification (default is 3306)
% added support for writing multiple events to mysql, each requires a single query
%
% Revision 1.10  2007/06/12 16:34:49  roboos
% first implementation of writing to a mysql database: approx. 30ms per event
%
% Revision 1.9  2007/06/07 12:43:38  chrhes
% added an option for specifying the maximum length of the event queue for
% the case where the events are being written (appended) to a .mat file
%
% Revision 1.8  2007/06/06 20:14:49  chrhes
% fixed a small bug to do with sting comparison
%
% Revision 1.7  2007/06/06 16:00:31  chrhes
% updated some documentation
%
% Revision 1.6  2007/06/06 15:55:30  chrhes
% extended functionality for serial port writing so that the code recognises
% when to write the whole event structure or only a single control character
%
% Revision 1.5  2007/06/06 15:45:53  chrhes
% implemented option for writing events to the serial port
%
% Revision 1.4  2007/06/06 12:38:57  roboos
% write events to uncompressed matlab v6 file
%
% Revision 1.3  2007/06/06 07:14:29  roboos
% switched to using filetype_check_uri for detection and parsing of filename
% switched to using external struct2msg function
% implemented fcdc_fifo
% added optinoal input arguments (key-value pairs)
%
% Revision 1.2  2007/05/31 09:54:05  roboos
% implemented writing events to a plain matlab file, the default is to append them
%
% Revision 1.1  2007/05/31 09:14:34  roboos
% initial implementation, sofar only tcpsocket
%

global event_queue   % for fcdc_global
global db_blob       % for fcdc_mysql
if isempty(db_blob)
  db_blob = 0;
end

if iscell(filename)
  % use recursion to write to multiple event targets
  for i=1:numel(filename)
    write_event(filename{i}, event, varargin);
  end
  return
end

% set the defaults
eventformat = keyval('eventformat', varargin); if isempty(eventformat), eventformat = filetype(filename); end
swapping    = keyval('swapping',    varargin); if isempty(swapping),    swapping = 'native';              end
append      = keyval('append',      varargin); if isempty(append),      append = 'yes';                   end
maxqlength  = keyval('maxqlength',  varargin); if isempty(maxqlength),  maxqlength = Inf;                 end

switch eventformat
  case 'disp'
    % display it on screen, this is only for debugging
    disp(event);

  case 'fcdc_global'
    % store it in a global variable, this is only for debugging
    if isempty(event_queue) || ~isstruct(event_queue)
      event_queue = event;
    else
      event_queue = appendevent(event_queue, event);
    end

  case 'fcdc_rfb'
   % remote frame buffer, i.e. VNC server on another computer
   [password, host, port] = filetype_check_uri(filename);
   rfbevent(sprintf('%s:%d', host, port), password, event.type, event.value);
   
  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    [host, port] = filetype_check_uri(filename);

    type = {
      'char'
      'uint8'
      'uint16'
      'uint32'
      'uint64'
      'int8'
      'int16'
      'int32'
      'int64'
      'single'
      'double'
      };

    wordsize = {
      1 % 'char'
      1 % 'uint8'
      2 % 'uint16'
      4 % 'uint32'
      8 % 'uint64'
      1 % 'int8'
      2 % 'int16'
      4 % 'int32'
      8 % 'int64'
      4 % 'single'
      8 % 'double'
      };

    for i=1:length(event)
      evt = [];
      buf = [];
      bufsize = 0;
      
      % convert the field "type" into the message representation
      this_type = class(event(i).type);
      this_size = numel(event(i).type);
      switch this_type
        case 'char'
          buf = cat(2, buf, event(i).type);
        otherwise
          buf = cat(2, buf, typecast(event(i).type, 'uint8'));
      end
      bufsize = bufsize + wordsize{strcmp(type, this_type)}*this_size;
      evt.type_type  = find(strcmp(type, this_type))-1; % zero-offset
      evt.type_numel = this_size;
      
      % convert the field "value" into the message representation
      this_type = class(event(i).value);
      this_size = numel(event(i).value);
      event(i).value = event(i).value(:)';  % it must be represented as a vector
      switch this_type
        case 'char'
          buf = cat(2, buf, event(i).value);
        otherwise
          buf = cat(2, buf, typecast(event(i).value, 'uint8'));
      end
      bufsize = bufsize + wordsize{strcmp(type, this_type)}*this_size;
      evt.value_type  = find(strcmp(type, this_type))-1; % zero-offset
      evt.value_numel = this_size;
      
      % although sample, offset and duration do not play a role in BCI2000,
      % they must exist for the buffer
      if ~isfield(event(i), 'sample') || isempty(event(i).sample)
        event(i).sample = 0;
      end
      if ~isfield(event(i), 'offset') || isempty(event(i).offset)
        event(i).offset = 0;
      end
      if ~isfield(event(i), 'duration') || isempty(event(i).duration)
        event(i).duration = 0;
      end
      
      % the other fields are simple, because they have a fixed type and only a single elements
      evt.sample   = int32(event(i).sample);
      evt.offset   = int32(event(i).offset);
      evt.duration = int32(event(i).duration);

      evt.bufsize = bufsize;
      evt.buf     = uint8(buf);
      
      % check if buffer is open and keep trying to fill it...
      try
        
        if strcmp(append,'no')        
          buffer('flush_evt', [], host, port);  % flush event
        end

        buffer('put_evt', evt, host, port);  % indices should be zero-offset
        
      catch
        
        % retrieve hostname
        [ret, hname] = system('hostname');
        if ret ~= 0,
          if ispc
            hname = getenv('COMPUTERNAME');
          else
            hname = getenv('HOSTNAME');
          end
        end

        if strcmpi(host,'localhost') || strcmpi(host,hname)
          
          warning('starting fieldtrip buffer on localhost');

          % try starting a local buffer
          buffer('tcpserver', 'init', host, port);
          pause(1);

          % write packet until succeed
          bhdr = false;
          while ~bhdr
            try
              bhdr = true;

              % try writing a dummy header
              dumhdr.fsample   = 0;
              dumhdr.nchans    = 0;
              dumhdr.nsamples  = 0;
              dumhdr.nevents   = 0;
              dumhdr.data_type = 0;

              buffer('put_hdr', dumhdr, host, port);
              buffer('put_evt', evt, host, port);  % indices should be zero-offset

            catch
              bhdr = false;
            end
          end
        end
      end
    end

  case 'fcdc_serial'
    % serial port on windows or linux platform
    s = [];
    [port, opt] = filetype_check_uri(filename);
    % determine whether any serial port objects are already associated with the
    % target serial port
    temp = instrfind;
    if isa(temp,'instrument')
      % find all serial ports
      i1 = strcmp('serial',{temp(:).Type});
      if any(i1)
        % find all serial ports whose name matches that of the specified port
        i2 = strmatch(lower(['Serial-',port]),lower({temp(find(i1)==1).Name}));
        % set s to the (first) matching port if present (and open if necessary)
        if ~isempty(i2)
          s = temp(i2(1));
          if ~strcmp(s.Status,'open'), fopen(s); end;
        end
      end
    end
    % create, configure a serial port object if necessary and open the port
    if ~isa(s,'serial')
      s = serial(port);
      if ~isempty(opt) && iscell(opt), s = set(s,opt); end;
      fopen(s);
    end

    %     % convert the event structure into an appropriate message
    %     if isfield(event,'type') && strcmp(event.type,'ctrlchar')
    %       % use only a single control character
    %       msg = char(event.value(1));
    %     else
    %       % convert the entire event structure into a message
    %       msg = struct2msg(event);
    %     end

    % write the contents of the field event.value to the serial port as a string
    if isfield(event,'value') && ~isempty(event.value)
      msg = char(event.value);
    else
      msg = [];
    end
    % write the message to the serial port
    if ~isempty(msg) && isa(s,'serial') && strcmp(s.Status,'open')
      %fprintf(s,msg,'async');
      fprintf(s,msg);
    else
      error('could not write event to serial port');
    end

  case 'fcdc_mysql'
    % write to a MySQL server listening somewhere else on the network
    db_open(filename);
    for i=1:length(event)
      if db_blob
        % insert the structure into the database table as a binary blob
        db_insert_blob('fieldtrip.event', 'msg', event(i));
      else
        % make a structure with the same elements as the fields in the database table
        s = struct;
        % these fields also exist as elements in the table and as such can be used for filtering
        if isa(event(i).type, 'char')
          s.type = event(i).type;
        end
        if isa(event(i).value, 'numeric') && numel(event(i).value)==1
          s.value = event(i).value;
        end
        if isa(event(i).sample, 'numeric') && numel(event(i).sample)==1
          s.sample = event(i).sample;
        end
        if isa(event(i).sample, 'numeric') && numel(event(i).offset)==1
          s.offset = event(i).offset;
        end
        if isa(event(i).sample, 'numeric') && numel(event(i).duration)==1
          s.duration = event(i).duration;
        end
        % insert the structure into the database table
        db_insert('fieldtrip.event', s);
      end
    end

   case 'fcdc_fifo'
    
      % these are opened in blocking mode, i.e. reading/writing will block until boths sides are connected
      fifo = filetype_check_uri(filename);
      
      if ~exist(fifo,'file')
          warning('the FIFO %s does not exist; attempting to create it', fifo);          
          system(sprintf('mkfifo -m 0666 %s',fifo));          
      end

      fid = fopen(fifo, 'w');
      for i=1:length(event)

        try
          % convert the event into a network message
          msg = mxSerialize(event(i));
          num = fwrite(fid, msg, 'uint8');
        catch
          warning(lasterr);
        end

        if num~=length(msg)
          error('problem writing to FIFO %s', fifo);
        end
      end
      fclose(fid);
      
    case 'fcdc_tcp'

        % TCP network socket
        [host, port] = filetype_check_uri(filename);

        con=pnet('tcpconnect',host,port);

        pnet(con,'setwritetimeout',1);

        if con~=-1,

            try % Failsafe

                for i=1:length(event)

                    % convert the event into a network message
                    msg = mxSerialize(event(i));

                    % tell the message daemon that a message will be sent, and send it
                    pnet(con,'printf',num2str(msg));
                    pnet(con,'printf','\n');
                end
%            catch             
%                warning(lasterr);
            end
            
            pnet(con,'close');
        end

    case 'fcdc_udp'

        % UDP network socket
        
        [host, port] = filetype_check_uri(filename);
        udp=pnet('udpsocket',port);

        if udp~=-1,
            try % Failsafe

                for i=1:length(event)

                    % convert the event into a network message
                    msg = mxSerialize(event(i));

                    % tell the message daemon that a message will be sent, and send it
                    pnet(udp,'write',uint8(msg),1000);
                    pnet(udp,'writepacket',host,port);   % Send buffer as UDP packet to host
                end

            catch
              warning(lasterr);
            end
            pnet(udp,'close');
        end


    otherwise
        % assume that it is a file. Since the file probably does not yet
        % exist, determine its type by only looking at the extension
        if filetype_check_extension(filename, '.mat')
            % write the events to a matlab file
            if exist(filename,'file') && strcmp(append, 'yes')
                try
                    tmp = load(filename, 'event');
                    event = cat(1, tmp.event(:), event(:));
                catch
                    event = event(:);
                end
                % optionally restric the length of the event queue to flush old events
                if isfinite(maxqlength) && isreal(maxqlength) && (maxqlength>0) && (length(event)>maxqlength)
                    event = event(end-maxqlength+1:end);
                    % NOTE: this could be done using the filter event function, but
                    % then this is just a temporary solution that will probably be
                    % removed in a future versions of the code
                end
                save(filename, 'event', '-append', '-v6');
                % NOTE: the -append option in this call to the save function does
                % not actually do anything useful w.r.t. the event variable since the
                % events are being appended in the code above and the the save function
                % will just overwrite the existing event variable in the file.
                % However, if there are other variables in the file (whatever) then the
                % append option preservs them
            else
                save(filename, 'event', '-v6');
            end
        else
            error('unsupported file type')
        end
end
