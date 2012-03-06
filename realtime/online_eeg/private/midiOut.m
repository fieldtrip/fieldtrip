function midiOut(varargin)
%midiOut  --  connect to a MIDI output device and send notes
%
% Usage:
% >> S = midiOut('L')   or without arguments
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

% (C) 2010 Stefan Klanke

error 'M-file called: It seems you have not compiled the midiOut MEX file.'



