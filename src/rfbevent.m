function [varargout] = rfbevent(varargin)

% RFBEVENT sends a keyboard or mouse event to a VNC server
%
% RFB ("remote frame buffer") is a simple protocol for remote access to
% graphical user interfaces. Because it works at the framebuffer level it
% is applicable to all windowing systems and applications, including X11,
% Windows and Macintosh. RFB is the protocol used in VNC (Virtual Network
% Computing).
%
% The remote endpoint where the user sits (i.e. the display plus keyboard
% and/or pointer) is called the RFB client or viewer. The endpoint where
% changes to the framebuffer originate (i.e. the windowing system and
% applications) is known as the RFB server.
%
% Use as
%   rfbevent(display, passwd, eventtype, eventvalue, ...)
%
% Some examples
%   rfbevent('vncserver:5901', 'yourpasswd', 'Text',   'xclock')          % type multiple characters
%   rfbevent('vncserver:5901', 'yourpasswd', 'Button', 'Return')          % single key event, press and release
%   rfbevent('vncserver:5901', 'yourpasswd', 'Button', 'Return',  0)      % single key event, press and release
%   rfbevent('vncserver:5901', 'yourpasswd', 'Button', 'Return',  1)      % single key event, press only
%   rfbevent('vncserver:5901', 'yourpasswd', 'Button', 'Return', -1)      % single key event, release only
%   rfbevent('vncserver:5901', 'yourpasswd', 'Pointer', [20 100])         % only mouse position
%   rfbevent('vncserver:5901', 'yourpasswd', 'Pointer', [20 100 1])       % mouse position and button 1, press and release
%   rfbevent('vncserver:5901', 'yourpasswd', 'Pointer', [20 100 1],  0)   % mouse position and button 1, press and release
%   rfbevent('vncserver:5901', 'yourpasswd', 'Pointer', [20 100 1],  1)   % mouse position and button 1, press only
%   rfbevent('vncserver:5901', 'yourpasswd', 'Pointer', [20 100 1], -1)   % mouse position and button 1, release only
%
% Note that the password has to be represented as plain text in the matlab
% script/function that is using RFBEVENT, which poses a potential security
% problem. The password is sent over the network to the VNC server after
% being encrypted.
%
% This implements the KeyEvent and PointerEvent messages according to
% "The RFB Protocol" by Tristan Richardson (RealVNC Ltd, formerly of
% Olivetti Research Ltd / AT&T Labs Cambridge) Version 3.8 (Last updated 8
% June 2007), http://www.realvnc.com/docs/rfbproto.pdf

% Copyright (C) 2008, Robert Oostenveld
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

% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);

  % the following is specific for this particular mex file
  mex rfbevent.c d3des.c

  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end

