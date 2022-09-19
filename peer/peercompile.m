function peercompile(cc)

% PEERCOMPILE  This script/function is used for compiling and linking the 'peer' MEX file.
%
% On Linux and macOS, you should just be able to call this function without arguments.
%
% On Windows, you need to select a compiler using one of the following options:
%   compile('bcb')   - Borland C++ Builder
%   compile('bcc55') - Borland C++ 5.5 (free command line tools)
%   compile('mingw') - MinGW (GCC port without cygwin dependency)
%   compile('vc')    - Visual Studio C++ 2005 or 2008
%
% Please note that this script does NOT set up your MEX environment for you, so in case
% you haven't selected a C compiler on Windows yet, you need to type 'mex -setup' first
% to choose either the Borland or Microsoft compiler. If you want to use MinGW, you also
% need to install Gnumex (http://gnumex.sourceforget.net), which comes with its own
% procedure for setting up the MEX environment.

% -----------------------------------------------------------------------
% Copyright (C) 2010, Stefan Klanke
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
%
% $Id$
% -----------------------------------------------------------------------

% You can tweak this a bit for setting platform-independent options, e.g for optimisation
% or debugging. The default is to enable debug information (-g).
% The -I. is for including the header files.
cflags = '-I. -g';

% change to the directory containing the source code
cwd = pwd;
src = fullfile(fileparts(which(mfilename)),'src');
cd(src);

if ispc
    % On Windows, we need to link against different libraries
    % If you are having problems compiling the peer mex file, you can try to tweak
    % the variables 'extra_cflags' and 'ldflags' specific to your compiler,
    % or you can add your own compiler flags.
    
    if strcmp(computer, 'PCWIN')
        extra_cflags = '-I../pthreads-win32/include -DSYSLOG=3';
    elseif strcmp(computer, 'PCWIN64')
        extra_cflags = '-I../pthreads-win64/include -DSYSLOG=3';
    end
    suffix = 'obj';
    
    if nargin<1	        % this is just to make sure the switch statement works
        cc = 'NONE';    % but you can also put a default entry here
    end
    switch upper(cc)
        case 'BCB'
            ldflags = '-L../src -L../pthreads-win32/lib -lpthreadVC2.bcb';
        case 'BCC55'
            ldflags = '-L../src -L../pthreads-win32/lib -lpthreadVC2_bcc55';
        case 'MINGW'
            %ldflags = '-L../src -L../pthreads-win32/lib -lpthreadGC2';
            % For MinGW/Gnumex, it seems to be easier to just directly refer to the archives
            if strcmp(computer, 'PCWIN')
                ldflags = '../pthreads-win32/lib/libpthreadGC2.a C:/mingw/lib/libws2_32.a';
            elseif strcmp(computer, 'PCWIN64')
                ldflags = '../pthreads-win64/lib/libpthread.a C:/MinGW64/x86_64-w64-mingw32/lib/libws2_32.a';
            end
        case 'VC'
            if strcmp(computer, 'PCWIN')
                ldflags = '-L../src -L../pthreads-win32/lib -lpthreadVC2';
            elseif strcmp(computer, 'PCWIN64')
                ldflags = '-L../src -L../pthreads-win64/lib -lpthreadVC64';
            end
        otherwise
            error 'On a PC, you must call this function with a string argument to select a compiler';
    end
else
    % On POSIX systems such as MacOS X and Linux, the following should work without tweaking
    ldflags = '';
    extra_cflags = '-DSYSLOG=3 -g';
    suffix = 'o';
end

% If you want to add a new helper function to the MEX file, you should just add
% its name here (without the .c suffix)
helpers = {'announce' 'discover' 'expire' 'extern' ...
    'peerinit' 'util' 'udsserver' 'tcpserver' 'tcpsocket' ...
    'security' 'localhost' 'smartshare' 'smartmem' 'smartcpu' 'connect'};

headers = {'compiler' 'peer' 'platform_includes' 'extern' 'platform' ...
    'swapbytes'};

%%%% Please do not change anything below this line %%%%

% this will become the list of objects files for inclusion during linking
allObjects = '';

% get the timestamp of each file
for i=1:length(helpers)
    helperinfo(i) = dir([helpers{i} '.c']);
end
helperdate = [helperinfo.datenum];
for i=1:length(headers)
    headerinfo(i) = dir([headers{i}, '.h']);
end
headerdate = [headerinfo.datenum];

% This is for compiling all the helper functions (no linking yet).
for i=1:length(helpers)

    % get the timestamp of the object file
    objinfo = dir([helpers{i}, '.' suffix]);
    if isempty(objinfo)
        objdate = -inf;
    else
        objdate = objinfo.datenum;
    end
    recompile = false;
    recompile = recompile | helperdate(i)>objdate;
    recompile = recompile | any(headerdate>objdate);
    
    if ~recompile
        fprintf(1,'%s.%s is up to date\n', helpers{i}, suffix);
    else
        fprintf(1,'Compiling helper functions in %s ...\n', helpers{i});
        cmd = sprintf('mex -c %s %s %s.c', cflags, extra_cflags, helpers{i});
        disp(cmd);
        eval(cmd);
    end
    % append newly created object file to the list of files we need to link
    allObjects = sprintf('%s %s.%s', allObjects, helpers{i}, suffix);
end

% This is for compiling peer.c and linking everything into the MEX file.
fprintf(1,'Compiling and linking MEX file:\n');
cmd = sprintf('mex -outdir ../private %s %s peer.c %s %s',cflags,extra_cflags,allObjects,ldflags);
disp(cmd);
eval(cmd);

cmd = sprintf('mex -outdir ../private %s %s memprofile.c %s %s',cflags,extra_cflags,allObjects,ldflags);
disp(cmd);
eval(cmd);

cmd = sprintf('mex -outdir ../private %s %s ft_getopt.c %s %s',cflags,extra_cflags,allObjects,ldflags);
disp(cmd);
eval(cmd);

cmd = sprintf('mex -outdir ../private %s %s watchdog.c %s %s',cflags,extra_cflags,allObjects,ldflags);
disp(cmd);
eval(cmd);

cmd = sprintf('mex -outdir ../private %s %s time.c %s',cflags,extra_cflags,ldflags);
disp(cmd);
eval(cmd);

% change back to the original directory
cd(cwd);
