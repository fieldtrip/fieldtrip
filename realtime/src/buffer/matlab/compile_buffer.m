function compile_buffer(cc, libpath)

% COMPILE is used for compiling and linking the 'buffer' MEX file.
%
% On Linux and MacOS X, you should just be able to call this function without arguments.
%
% On Windows, you can select a compiler using one of the following options:
%   compile('lcc')   - LCC compiler (shipped with Matlab on 32bit platforms)
%   compile('bcb')   - Borland C++ Builder
%   compile('bcc55') - Borland C++ 5.5 (free command line tools)
%   compile('mingw') - MinGW (GCC port without cygwin dependency)
%   compile('vc')    - Visual Studio C++ 2005 or 2008
%
% Please note that this script does NOT set up your MEX environment for you, so in case
% you haven't selected the C compiler on Windows yet, you need to type 'mex -setup' first
% to choose either the Borland or Microsoft compiler. If you want to use MinGW, you also
% need to install Gnumex (http://gnumex.sourceforget.net), which comes with its own
% procedure for setting up the MEX environment.

% Copyright (C) 2010-2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% You can tweak this a bit for setting platform-independent options, e.g for optimisation
% or debugging. The default is to enable debug information (-g).
% The -I../src is for including the header files of the buffer C library.
cflags = '-I../src -g';

if ispc
  % On Windows, we need to link against different libraries
  % If you are having problems compiling the buffer, you can try to tweak
  % the variables 'extra_cflags' and 'ldflags' specific to your compiler,
  % or you can add your own compiler flags.
  
  if strcmp(computer,'PCWIN64')
    extra_cflags = '-I../../external/pthreads-win64/include';
    amd64 = true;
  else
    extra_cflags = '-I../../external/pthreads-win32/include';
    amd64 = false;
  end
  suffix = 'obj';
  
  if nargin<1
    cc = 'LCC';
    if amd64
      error('Running on 64-bit Windows. Please setup your compiler environment for either Visual C++ or MinGW64, and then call compile(''vc'') or compile(''mingw'').');
    end
  end
  switch upper(cc)
    case 'BCB'
      ldflags = '-L../../external/pthreads-win32/lib -lpthreadVC2.bcb';
    case 'BCC55'
      ldflags = '-L../../external/pthreads-win32/lib -lpthreadVC2_bcc55';
    case 'MINGW'
      % For MinGW/Gnumex, it seems to be easier to just directly refer to the archives, since
      % the MEX tools expect libraries to end with .lib, whereas MinGW uses the .a suffix.
      if amd64
        ldflags = '../../external/pthreads-win64/lib/libpthread.a';
        ws2lib = 'C:/mingw64/x86_64-w64-mingw32/lib/libws2_32.a';
      else
        ldflags = '../../external/pthreads-win32/lib/libpthreadGC2.a';
        ws2lib = 'C:/mingw/lib/libws2_32.a';
      end
      if nargin<2
        if ~exist(ws2lib,'file')
          fprintf(1,'Library file WS2_32 does not exist in guessed location %s.\n', ws2lib);
          fprintf(1,'Please re-run this script with a second argument pointing to your MinGW\n');
          fprintf(1,'installation''s LIB folder, e.g.,   compile(''mingw'',''D:\MinGW\lib''\n');
          error('unsuccesful');
        end
      else
        ws2lib = [libpath '/libws2_32.a'];
        if ~exist(ws2lib,'file')
          fprintf(1,'Library file %s does not exist.\n', ws2lib);
          error('unsuccesful');
        end
      end
      ldflags = [ldflags ' ' ws2lib];
    case 'VC'
      if amd64
        ldflags = '-L../../external/pthreads-win64/lib -lpthreadVC2 ws2_32.lib';
      else
        ldflags = '-L../../external/pthreads-win32/lib -lpthreadVC2 ws2_32.lib ';
      end
    case 'LCC'
      ldflags = '-L../../external/pthreads-win32/lib  -lpthreadGC2_lcc';
      ldflags = [ldflags ' "' matlabroot '\sys\lcc\lib\wsock32.lib"'];
      %ldflags = [ldflags ' "' matlabroot '\sys\lcc\lib\kernel32.lib"'];
    otherwise
      error 'Unsupported compiler - select one of LCC, VC, MINGW, BCC55, BCB';
  end
else
  % On POSIX systems such as MacOS X and Linux, the following should work without tweaking
  ldflags = '-lpthread';
  extra_cflags = '';
  suffix = 'o';
end

% This is the list of C files from the low-level buffer implementation that we need
libfuncs = {
  'util'
  'printstruct'
  'clientrequest'
  'tcprequest'
  'tcpserver'
  'tcpsocket'
  'extern'
  'dmarequest'
  'endianutil'
  'cleanup'
  'clock_gettime'
  };

% If you want to add a new helper function to the MEX file, you should just add
% its name here (without the .c suffix)
helpers = {'buffer_gethdr'
  'buffer_getdat'
  'buffer_getevt'
  'buffer_puthdr'
  'buffer_putdat'
  'buffer_putevt'
  'buffer_flushhdr'
  'buffer_flushdat'
  'buffer_flushevt'
  'buffer_waitdat'
  'buffer_mxutils'
  };


%%%% Please do not change anything below this line %%%%

% this will become the list of objects files for inclusion during linking
allObjects = '';

% This is for compiling all the library functions (no linking yet).
for i=1:length(libfuncs)
  fprintf(1,'Compiling library functions in %s.c ...\n', libfuncs{i});
  cmd = sprintf('mex -c %s %s ../src/%s.c', cflags, extra_cflags, libfuncs{i});
  eval(cmd);
  
  % append newly created object file to the list of files we need to link
  allObjects = sprintf('%s %s.%s', allObjects, libfuncs{i}, suffix);
end

fprintf(1,'\n');

% This is for compiling all the helper functions (no linking yet).
for i=1:length(helpers)
  fprintf(1,'Compiling helper functions in %s.c ...\n', helpers{i});
  cmd = sprintf('mex -c %s %s %s.c', cflags, extra_cflags, helpers{i});
  eval(cmd);
  
  % append newly created object file to the list of files we need to link
  allObjects = sprintf('%s %s.%s', allObjects, helpers{i}, suffix);
end

% This is for compiling buffer.c and linking everything into the MEX file.
fprintf(1,'Compiling and linking MEX file:\n');
cmd = sprintf('mex %s %s buffer.c %s %s',cflags,extra_cflags,allObjects,ldflags);
disp(cmd)
eval(cmd);
