function varargout = midiOut(varargin)

% midiOut -- connect to a MIDI output device and send notes
%
% >> S = midiOut('L')
% Returns a list of MIDI devices as a structure array. Fields are 'index',
% 'name', 'input' and 'output'. The last two are flags (1 or 0) that determine
% whether the device can be used as an input, or output, or both.
%
% >> midiOut('O', index)
% Opens a MIDI device by its index (as returned by aforementioned 'L' call).
% If the device is already opened, nothing will happen. If another device is
% opened, that one will be closed first.
%
% >> midiOut('C')
% Closes the current connection.
%
% >> midiOut('+', channel, notes, velocities)
% Send one or more 'note on' messages at the given channel (1..16).
% 'notes' and 'velocities' must be given as double-precision arrays with
% equal number of elements. Messages will be transmitted in order.
%
% >> midiOut('-', channel, notes, velocities)
% Send one or more 'note off' messages, same parameters as for 'note on'.
%
% >> midiOut('.', channel)
% Send an 'all notes off' message at the given channel.
%
% >> midiOut('P', channel, program)
% Send a 'program change' to the given channel.
%
% >> midiOut(msg)
% Send raw data. 'msg' must be of type 'uint8' and have a multiple of 3 elements.
%
% For example, play the note C4 (261.63 Hz) on the piano
% >> midiOut('O', 2)
% >> midiOut('+', 1, 60, 127)
%
% See also MIDIIN

% Copyright (C) 2010, Stefan Klanke
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
