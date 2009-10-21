function [varargout] = buffer(varargin)

% BUFFER manages and accesses the realtime data acquisition buffer
% This function is implented as mex file.
%
% Use as
%   retval = buffer(cmd, detail, host, port)
%
% To read data from a buffer server over the network
%   hdr = buffer('get_hdr', [],     host, port)
%   dat = buffer('get_dat', datsel, host, port)
%   evt = buffer('get_evt', evtsel, host, port)
%
% The selection for data and events should be zero-offset and contain
%   datsel = [begsample endsample]
%   evtsel = [begevent  endevent]
%
% To write data to a buffer server over the network
%   buffer('put_hdr', hdr, host, port)
%   buffer('put_dat', dat, host, port)
%   buffer('put_evt', evt, host, port)
%
% To implement a local buffer server and have other clients
% connect to it at a specified network port
%   buffer('tcpserver', 'init', [], port)
%   buffer('tcpserver', 'exit', [], port)
%
% To implement a local acquisition client and have it buffer the
% data locally
%   buffer(acqclient, 'init')
%   buffer(acqclient, 'exit')
%
% To implement a local acquisition client and have it send the data
% through the network to a buffer server
%   buffer(acqclient, 'init', host, port)
%   buffer(acqclient, 'exit', host, port)
%
% Supported acquisition clients should be
%  'sinewave'
%  'ctf'
%  'brainamp'
%  'biosemi'
%  'tmsi'
%  'eldith'

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: buffer.m,v $
% Revision 1.7  2009/01/21 19:58:07  roboos
% also call the function after compilation
%
% Revision 1.6  2009/01/20 15:38:37  roboos
% try to compile the mex file on the fly
%
% Revision 1.5  2008/10/24 07:32:05  roboos
% should be a function, not a script
%
% Revision 1.4  2008/07/08 20:31:34  roboos
% extended the list of acquisition threads in the mex file, also added all of them to thread stopping
%
% Revision 1.3  2008/05/22 09:23:09  roboos
% updated documentation
%
% Revision 1.2  2008/03/10 10:41:58  roboos
% updated documentation
%
% Revision 1.1  2008/03/09 22:53:04  roboos
% new function, only placeholder for the documentation since the actual implementation is in a mex file
%

% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);

  % the following is specific for this particular mex file
  mex -I../src -I../../pthreads-win32/include -c buffer_gethdr.c
  mex -I../src -I../../pthreads-win32/include -c buffer_getdat.c
  mex -I../src -I../../pthreads-win32/include -c buffer_getevt.c
  mex -I../src -I../../pthreads-win32/include -c buffer_getprp.c
  mex -I../src -I../../pthreads-win32/include -c buffer_puthdr.c
  mex -I../src -I../../pthreads-win32/include -c buffer_putdat.c
  mex -I../src -I../../pthreads-win32/include -c buffer_putevt.c
  mex -I../src -I../../pthreads-win32/include -c buffer_putprp.c
  mex -I../src -I../../pthreads-win32/include -c buffer_flushhdr.c
  mex -I../src -I../../pthreads-win32/include -c buffer_flushdat.c
  mex -I../src -I../../pthreads-win32/include -c buffer_flushevt.c

  % this is for the demo acquisition threads
  mex -I../src -I../../pthreads-win32/include -c ../test/sinewave.c
  mex -I../src -I../../pthreads-win32/include -c ../test/event.c

  % link all the objects together into a mex file
  if ispc
    % this is needed for Borland on windows, but not for other platforms
    mex -I../src -I../../pthreads-win32/include -c ../src/gettimeofday.c -I../src -I../../pthreads-win32/include
    % link all the objects together into a mex file
    mex -I../src -I../../pthreads-win32/include -L../src -L../../pthreads-win32/lib buffer.c event.obj sinewave.obj gettimeofday.obj buffer_gethdr.obj buffer_getdat.obj buffer_getevt.obj buffer_getprp.obj buffer_flushhdr.obj buffer_flushdat.obj buffer_flushevt.obj buffer_puthdr.obj buffer_putdat.obj buffer_putevt.obj buffer_putprp.obj -lbuffer -lpthreadVC2.bcb
  else
    % link all the objects together into a mex file
    mex -I../src -L../src buffer.c event.o sinewave.o buffer_gethdr.o buffer_getdat.o buffer_getevt.o buffer_getprp.o buffer_flushhdr.o buffer_flushdat.o buffer_flushevt.o buffer_puthdr.o buffer_putdat.o buffer_putevt.o buffer_putprp.o -lbuffer
  end

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
