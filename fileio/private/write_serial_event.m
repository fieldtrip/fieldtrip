function write_serial_event(filename, event)

% WRITE_SERIAL_EVENT
%
% changed A.Hadjipapas 2010
%
% write to phyiscal serial port
% serial port on windows or linux platform

% Copyright (C) 2007, Christian Hesse
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

[port, opt] = filetype_check_uri(filename);
temp=instrfind; %% find isntruments including serial ports

%% check if a serial port with the specified name is already available
if isempty(temp)||~any(strcmpi(temp.Port,port))
  %% if not make one
  s=serial(port);
  %% if user has specified Baudrate, apply that, %%  FIX THIS: currently no
  %% automatic implementation for more options
  if ~isempty(opt)
    s.Baudrate =str2double(opt{2}); %% set Baudrate using optional parameters
  end
  fopen(s); %% open serial port
end

%% write the contents of the field event.value to the serial port as a string
if isfield(event,'value') && ~isempty(event.value)
  msg = event.value;
else
  msg = [];
end

%% write the message to the serial port if it is open
if strcmp(s.Status,'open')
  %% if numeric use fwrite otherwise fprintf
  if isnumeric(msg)
    fwrite(s,msg);
  end

  if ischar(msg)
    fprintf(s,'%s',msg);
  end
else
  error('could not write event to serial port');
end

%% this is the old code
% 
%     % serial port on windows or linux platform
%     s = [];
%     [port, opt] = filetype_check_uri(filename);
% 
%     % determine whether any serial port objects are already associated with the target serial port
%     temp = instrfind;
%     if isa(temp,'instrument')
%       % find all serial ports
%       i1 = strcmp('serial',{temp(:).Type});
%       if any(i1)
%         % find all serial ports whose name matches that of the specified port
%         i2 = strmatch(lower(['Serial-',port]),lower({temp(find(i1)==1).Name}));
%         % set s to the (first) matching port if present (and open if necessary)
%         if ~isempty(i2)
%           s = temp(i2(1));
%           if ~strcmp(s.Status,'open'), fopen(s); end;
%         end
%       end
%     end
% 
%     % create, configure a serial port object if necessary and open the port
%     if ~isa(s,'serial')
%       s = serial(port);
%       if ~isempty(opt) && iscell(opt), s = set(s,opt); end;
%       fopen(s);
%     end
% 
%     %     % convert the event structure into an appropriate message
%     %     if isfield(event,'type') && strcmp(event.type,'ctrlchar')
%     %       % use only a single control character
%     %       msg = char(event.value(1));
%     %     else
%     %       % convert the entire event structure into a message
%     %       msg = struct2msg(event);
%     %     end
% 
%     % write the message to the serial port
%     if isa(s,'serial') && strcmp(s.Status,'open')
%       fwrite(s,char(event.value));
%     else
%       error('could not write event to serial port');
%     end
