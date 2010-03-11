function peercompile

% PEERCOMPILE compile the low-level mex file implements that implements
% the TCP/IP, announce and discover services.

% -----------------------------------------------------------------------
% Copyright (C) 2010, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
% -----------------------------------------------------------------------

if ispc

  olddir = pwd;
  try
    % go into the directory containing the source code
    newdir = fullfile(fileparts(which(mfilename)), 'src');
    cd(newdir);

     mex -I. -I../pthreads-win32/include -c announce.c
     mex -I. -I../pthreads-win32/include -c discover.c
     mex -I. -I../pthreads-win32/include -c expire.c
     mex -I. -I../pthreads-win32/include -c extern.c
     mex -I. -I../pthreads-win32/include -c fairshare.c
     mex -I. -I../pthreads-win32/include -c peerinit.c
     mex -I. -I../pthreads-win32/include -c util.c
     mex -I. -I../pthreads-win32/include -c tcpserver.c
     mex -I. -I../pthreads-win32/include -c tcpsocket.c
     mex -I. -I../pthreads-win32/include -c security.c

     % link all the objects together into a mex file
     mex -I. -I../pthreads-win32/include -L../pthreads-win32/lib -output ../peer peer.c tcpserver.obj tcpsocket.obj discover.obj announce.obj expire.obj peerinit.obj util.obj extern.obj fairshare.obj security.obj -lpthreadVC2.bcb

    % also compile the memory profiler
    mex -I. -I../pthreads-win32/include -L../pthreads-win32/lib -output ../memprofile memprofile.c -lpthreadVC2.bcb

    % also compile the delayed exit
    mex -I. -I../pthreads-win32/include -L../pthreads-win32/lib -output ../delayedexit delayedexit.c -lpthreadVC2.bcb

    % return to the original directory
    cd(olddir);

  catch
    % return to the original directory
    cd(olddir);
    % report the error
    rethrow(lasterror);
  end


elseif isunix || ismac

  olddir = pwd;
  try
    % go into the directory containing the source code
    newdir = fullfile(fileparts(which(mfilename)), 'src');
    cd(newdir);

    % these are general
    mex  -I. -c announce.c
    mex  -I. -c discover.c
    mex  -I. -c expire.c
    mex  -I. -c extern.c
    mex  -I. -c fairshare.c
    mex  -I. -c peerinit.c
    mex  -I. -c util.c
    mex  -I. -c tcpserver.c
    mex  -I. -c tcpsocket.c
    mex  -I. -c security.c

    % link all the objects together into a mex file
    mex -I. -output ../peer peer.c tcpserver.o tcpsocket.o discover.o announce.o expire.o peerinit.o util.o extern.o fairshare.o security.o -lpthread

    % also compile the memory profiler
    mex -I. -output ../memprofile memprofile.c

    % also compile the delayed exit
    mex -I. -output ../delayedexit delayedexit.c

    % return to the original directory
    cd(olddir);

  catch
    % return to the original directory
    cd(olddir);
    % report the error
    rethrow(lasterror);
  end

end

