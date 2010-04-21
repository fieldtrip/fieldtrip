function event = read_serial_event(filename)

% READ_SERIAL_EVENT
%
% changed A.Hadjipapas 2010
%
% The only thing transmitted is the event.value (no info about sample) but it works

% Copyright (C) 2007, Christian Hesse
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

%% serial port on windows or linux platform
[port, opt] = filetype_check_uri(filename);
% determine whether any serial port objects are already associated with the
% target serial port
temp = instrfind;
if isempty(temp)||~any(strcmpi(temp.Port,port))
  s=serial(port);
  s.Baudrate =str2double(opt{2});
  fopen(s);
end

%% try to read a message from the serial port
if ~isempty(s.BytesAvailable) && s.BytesAvailable~=0
  %% NOTE: Here I only use fread,which delivers numerical outputs. I
  %% then check whether the message is longer than 3 and if it is
  %% I assume it was a string to begin with and then deocde it
  %% back to char
  msg = fread(s,s.BytesAvailable);
  if length(msg)>3

    msg=char(msg');
    event.type='char';
  else
    event.type='uint';
  end
  %% The only thing transmitted is the event.value
  event.value=msg;
  %% FIX THIS: here sample is always at a fixed value (i.e true sample (e.g. write_event) is not
  %% transmitted)
  event.sample=1;
  event.duration=[];
end;


%% this is the old code
%
%     % serial port on windows or linux platform
%     [port, opt] = filetype_check_uri(filename);
%     % determine whether any serial port objects are already associated with the
%     % target serial port
%     s = [];
%     temp = instrfind;
%     if isa(temp,'instrument')
%       % find all serial ports
%       i1 = strcmpi({temp(:).Type},'serial');
%       if any(i1)
%         % find all serial ports whose name matches that of the specified port
%         i2 = strmatch(lower(port),lower({temp(find(i1)).Name}));
%         % set s to the (first) matching port if present (and open if necessary)
%         if ~isempty(i2)
%           s = temp(i2(1));
%           if ~strcmp(s.Status,'open'), fopen(s); end;
%         end
%       end
%     end
%     % create, configure a serial port object if necessary and open the port
%     if ~isa(s,'serial')
%       s = serial(port);
%       if ~isempty(opt) && iscell(opt), s = set(s,opt); end;
%       fopen(s);
%     end
%     % try to read a message from the serial port
%     msg = [];
%     % FIXME: this currently assumes that all messages are terminated by the
%     % "newline" character (ascii character 10)
%     try
%       msg = fscanf(s,'%s\n');
%     end;
%     % convert message to event structure
%     event = msg2struct(msg);
