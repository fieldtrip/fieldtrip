function compile(cc)

% COMPILE  This script/function is used for compiling and linking the 'buffer' MEX file.
%
% On Linux and MacOS X, you should just be able to call this function without arguments.
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
%
% Please also note that you can only use this function, AFTER you have compiled the
% 'buffer' library from '../src' with the SAME compiler you use here.

% You can tweak this a bit for setting platform-independent options, e.g for optimisation 
% or debugging. The default is to enable debug information (-g).
% The -I../src is for including the header files of the buffer C library.
cflags = '-I../src -g'; 

if ispc
	% On Windows, we need to link against different libraries
	% If you are having problems compiling the buffer, you can try to tweak
	% the variables 'extra_cflags' and 'ldflags' specific to your compiler,
	% or you can add your own compiler flags.

	extra_cflags = '-I../pthreads-win32/include';
	suffix = 'obj';

	if nargin<1	        % this is just to make sure the switch statement works
		cc = 'NONE';    % but you can also put a default entry here
	end
	switch upper(cc)
		case 'BCB'
			ldflags = '-L../src -L../pthreads-win32/lib -lbuffer -lpthreadVC2.bcb';
		case 'BCC55'
			ldflags = '-L../src -L../pthreads-win32/lib -lbuffer -lpthreadVC2_bcc55';
		case 'MINGW'
			%ldflags = '-L../src -L../pthreads-win32/lib -lbuffer -lpthreadGC2';
			% For MinGW/Gnumex, it seems to be easier to just directly refer to the archives
			ldflags = '../src/libbuffer.a ../pthreads-win32/lib/libpthreadGC2.a C:/mingw/lib/libws2_32.a';
		case 'VC'
			ldflags = '-L../src -L../pthreads-win32/lib -lbuffer -lpthreadVC2 ';
		otherwise
			error 'On a PC, you must call this function with a string argument to select a compiler';
	end
else
	% On POSIX systems such as MacOS X and Linux, the following should work without tweaking
	ldflags = '-lbuffer -L../src';
	extra_cflags = '';
	suffix = 'o';
end

% If you want to add a new helper function to the MEX file, you should just add 
% its name here (without the .c suffix)
helpers = {'buffer_gethdr', 'buffer_getdat', 'buffer_getevt', 'buffer_getprp', ...
		   'buffer_puthdr', 'buffer_putdat', 'buffer_putevt', 'buffer_putprp', ...
		   'buffer_flushhdr', 'buffer_flushdat', 'buffer_flushevt'};


%%%% Please do not change anything below this line %%%%

% this will become the list of objects files for inclusion during linking
allObjects = '';

% This is for compiling all the helper functions (no linking yet).
for i=1:length(helpers)
	fprintf(1,'Compiling helper function %s ...\n', helpers{i});
	cmd = sprintf('mex -c %s %s %s.c', cflags, extra_cflags, helpers{i});
	eval(cmd);

	% append newly created object file to the list of files we need to link
	allObjects = sprintf('%s %s.%s', allObjects, helpers{i}, suffix);
end

% This is for compiling buffer.c and linking everything into the MEX file.
fprintf(1,'Compiling and linking MEX file:\n');
cmd = sprintf('mex -outdir ../../../fileio/private %s %s buffer.c %s %s',cflags,extra_cflags,allObjects,ldflags);
disp(cmd)
eval(cmd);
