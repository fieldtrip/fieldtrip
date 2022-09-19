function ft_write_event(filename, event, varargin)

% FT_WRITE_EVENT writes an event structure to a file, a message daemon listening on a
% network socked, or to another computer connected through the serial port. Note that
% this function is mostly for real-time streaming of events. For most data files on
% disk the writing of events is done simultaneously with the header and data in
% FT_WRITE_DATA.
%
% Use as
%   ft_write_event(filename, event, ...)
%
% The first argument is a string containing the filename. The second argument is a
% structure with the event. Multiple events can be represented as a structure array.
% Events are represented in the same format as those returned by FT_READ_EVENT.
%   event.type      = string
%   event.sample    = expressed in samples, the first sample of a recording is 1
%   event.value     = number or string
%   event.offset    = expressed in samples
%   event.duration  = expressed in samples
%   event.timestamp = expressed in timestamp units, which vary over systems (optional)
%
% Additional options should be specified in key-value pairs and can be
%   'eventformat'  = string, see below
%   'append'       = boolean, not supported for all formats
%
% Events can be written to special communication streams by specifying the target as
% URI instead of a filename. Supported are
%   buffer://<host>:<port>
%   fifo://<filename>
%   tcp://<host>:<port>
%   udp://<host>:<port>
%   mysql://<user>:<password>@<host>:<port>
%   rfb://<password>@<host>:<port>
%   serial:<port>?key1=value1&key2=value2&...
%   rfb://<password>@<host>:<port>
%
% See also FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, FT_WRITE_DATA

% Copyright (C) 2007-2020, Robert Oostenveld
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

global event_queue   % for fcdc_global
global db_blob       % for fcdc_mysql
if isempty(db_blob)
  db_blob = 0;
end

if iscell(filename)
  % use recursion to write to multiple event targets
  for i=1:numel(filename)
    ft_write_event(filename{i}, event, varargin);
  end
  return
end

% set the defaults
eventformat = ft_getopt(varargin, 'eventformat', ft_filetype(filename));
append      = ft_getopt(varargin, 'append', 'yes');
maxqlength  = ft_getopt(varargin, 'maxqlength', inf);

switch eventformat
  
  case 'empty'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % just pretend that we are writing the events, this is only for debugging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ft_info('Pretending to write %i events...\n', length(event));
    
  case 'fcdc_global'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store it in a global variable, this is only for debugging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(event_queue) || ~isstruct(event_queue)
      event_queue = event;
    else
      event_queue = appendstruct(event_queue, event);
    end
    
  case 'fcdc_rfb'
    % remote frame buffer, i.e. VNC server on another computer
    [password, host, port] = filetype_check_uri(filename);
    rfbevent(sprintf('%s:%d', host, port), password, event.type, event.value);
    
  case 'fcdc_buffer'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to a network transparent buffer for realtime analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [host, port] = filetype_check_uri(filename);
    
    if strcmp(append, 'no')
      buffer('flush_evt', [], host, port); % flush the existing events
    end
    
    % the MEX file now can handle various MATLAB types directly and respects the fields
    % sample, offset, duration
    %   -- these must all be numeric and non-empty (only first element is of interest)
    % type, value
    %   -- these can be strings or any numeric type (double, single, [u]int[8-64])
    %      will be transmitted as if vectorised
    
    % hack to fix empty event.offset and event.duration fields
    % Maybe should track down why they're empty? But this works for now
    % ES, 10-may-2019
    
    if ~isempty(event)
      for k = 1:numel(event)
        if isfield(event(k), 'offset') && isempty(event(k).offset)
          event(k).offset = 0;
        end
        if isfield(event(k), 'duration') && isempty(event(k).duration)
          event(k).duration = 0;
        end
      end
    end
    
    buffer('put_evt', event, host, port); % append the new events
    
    % SK: There was some code here for firing up a FieldTrip buffer locally,
    % but this is very likely to be senseless because we have no proper header
    % information here. Please explicitly use ft_create_buffer instead.
    
  case 'fcdc_serial'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to a serial port
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_serial_event(filename, event);
    
  case 'fcdc_mysql'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to a MySQL server listening somewhere else on the network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check that the required low-level toolbox is available
    ft_hastoolbox('mysql', 1);
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to a FIFO special file, which is a named pipe
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % these are opened in blocking mode, i.e. reading/writing will block until boths sides are connected
    fifo = filetype_check_uri(filename);
    
    if ~exist(fifo,'file')
      ft_warning('the FIFO %s does not exist; attempting to create it', fifo);
      system(sprintf('mkfifo -m 0666 %s',fifo));
    end
    
    fid = fopen_or_error(fifo, 'w');
    for i=1:length(event)
      
      try
        % convert the event into a network message
        msg = mxSerialize(event(i));
        num = fwrite(fid, msg, 'uint8');
      catch
        ft_warning(lasterr);
      end
      
      if num~=length(msg)
        ft_error('problem writing to FIFO %s', fifo);
      end
    end
    fclose(fid);
    
  case 'fcdc_tcp'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % serialize the events and write the message over a TCP network socket
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [host, port] = filetype_check_uri(filename);
    con = pnet('tcpconnect', host, port);
    
    pnet(con,'setwritetimeout',1);
    
    if con~=-1
      try % failsafe, don't give error when there is no connection
        for i=1:length(event)
          % convert the event into a network message
          msg = mxSerialize(event(i));
          
          % tell the message daemon that a message will be sent, and send it
          pnet(con,'printf',num2str(msg));
          pnet(con,'printf','\n');
        end
      catch
        ft_warning(lasterr);
      end % try
      
      pnet(con, 'close');
    end % if connected
    
  case 'fcdc_udp'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % serialize the events and write the message over a UDP network socket
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [host, port] = filetype_check_uri(filename);
    udp = pnet('udpsocket', port);
    
    if udp~=-1
      try % failsafe, don't give error when there is no connection
        for i=1:length(event)
          % convert the event into a network message
          msg = mxSerialize(event(i));
          
          % tell the message daemon that a message will be sent, and send it
          pnet(udp,'write',uint8(msg),1000);
          pnet(udp,'writepacket',host,port);   % Send buffer as UDP packet to host
        end
        
      catch
        ft_warning(lasterr);
      end % try
      pnet(udp,'close');
    end % if connected
    
    
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % since the file probably does not yet exist, determine its type by only looking at the extension
    if filetype_check_extension(filename, '.mat')
      % write the events to a MATLAB file
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
      if ~isempty(evt)
        ft_error('writing events to this fileformat is not supported here, please see FT_WRITE_DATA');
      end
      
    end % if file extension
    
end % switch eventformat
