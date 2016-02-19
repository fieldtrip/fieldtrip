function varargout = midiIn(varargin)

% midiIn -- connect to a MIDI input device and record notes or other events
%
% >> S = midiIn('L')
% Returns a list of MIDI devices as a structure array. Fields are 'index',
% 'name', 'input' and 'output'. The last two are flags (1 or 0) that determine
% whether the device can be used as an input, or output, or both.
%
% >> midiIn('O', index)
% Opens a MIDI device by its index (as returned by aforementioned 'L' call).
% If the device is already opened, nothing will happen. If another device is
% opened, that one will be closed first.
%
% >> midiIn('C')
% Closes the current connection.
%
% >> A = midiIn('G')
% Get all notes that were recorded since starting or since the previous call. This
% returns a matrix with 4 columns representing the channel, note, velocity and
% timestamp (in miliseconds).
%
% >> midiIn('F')
% Flush all recorded notes.
%
% >> midiIn('V', value)
% Toggle verbose mode, can be 0 or 1.
%
% For example, receive notes from the piano
% >> midiIn('O', 1)
% % play some melody
% >> a = midiIn('G')
%
% Subsequently play the sequence of notes on the piano
% >> midiOut('O', 2)
% >> note     = a(:,2);
% >> velocity = a(:,3);
% >> delay    = diff([a(1,4); a(:,4)])/1000;
% >> for i=1:length(note); pause(delay(i)); midiOut('+', 1, note(i), velocity(i)); end
%
% See also MIDIOUT

% Copyright (C) 2015, Robert Oostenveld
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

error('the %s mex file is not available for your platform (%s)', mfilename, mexext);
