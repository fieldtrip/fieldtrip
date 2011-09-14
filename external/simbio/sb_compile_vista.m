function sb_compile_vista(force)

% FT_COMPILE_VISTA is used for compiling most of the Vista MEX files that are used in FieldTrip
%
% Please note that this script does NOT set up your MEX environment for you, so in case
% you haven't selected the C compiler on Windows yet, you need to type 'mex -setup' first
% to choose either the Borland or Microsoft compiler. If you want to use MinGW, you also
% need to install Gnumex (http://gnumex.sourceforget.net), which comes with its own
% procedure for setting up the MEX environment.

% The logic in this script is to first build a list of files that actually need compilation for the
% particular platform that Matlab is running on, and then to go through that list.
% Functions are added to the list by giving their destination directory and (relative to that) the 
% name of the source file (without the .c). Optionally, you can specify a list of platform this
% file needs to be compiled on only, and a list of platforms where you don't compile it on.
% Finally, you can give extra arguments to the MEX command, e.g., for including other c-sources or
% giving compiler flags.

% Copyright (C) 2010, Stefan Klanke, 2011 Cristiano Micheli
%
% $Log$

% ft_compile_vista compiles the files belonging to the fileio Vista format library
% 
% vista/libvista.a
% vista/read_vista_mesh.cpp
% vista/write_vista_mesh.cpp
% vista/write_vista_vol.cpp
% vista/vistaprimitive.h
% vista/vistaprimitive.cpp
% vista/write_vista_vol.m

if nargin<1
   force=false;
end

% Possible COMPUTER types
% GLNX86
% GLNXA64
% PCWIN
% PCWIN64
% MAC
% MACI
% MACI64

% At the moment this works only for Linux 32-64 bit
% Vista folder has to be in the same folder as sb_compile_vista

% Step 1
% compile the libvista.a shared library for UNIX/LINUX
L = [];
L = add_mex_source(L,'vista','FIL_Vista_Attr.c',{'GLNX86', 'GLNXA64'},[],'-Ivista -c -o FIL_Vista_Attr.o');
L = add_mex_source(L,'vista','FIL_Vista_Basic.c',{'GLNX86', 'GLNXA64'},[],'-Ivista -c -o FIL_Vista_Basic.o');
L = add_mex_source(L,'vista','FIL_Vista_File.c',{'GLNX86', 'GLNXA64'},[],'-Ivista -c -o FIL_Vista_File.o');
L = add_mex_source(L,'vista','FIL_Vista_Graph.c',{'GLNX86', 'GLNXA64'},[],'-Ivista -c -o FIL_Vista_Graph.o');
L = add_mex_source(L,'vista','FIL_Vista_Image.c',{'GLNX86', 'GLNXA64'},[],'-Ivista -c -o FIL_Vista_Image.o');

oldDir = pwd;
[baseDir, myName] = fileparts(mfilename('fullpath'));
try
  compile_mex_list_vista(L, baseDir, force); 
catch
  me = lasterror;
  cd(oldDir);
  rethrow(me);
end
cd(baseDir);
% FIXME: the following should be done inside a function called add_shared_library
if isunix
  system('ar rcs libvista.a vista/FIL_Vista_Attr.o vista/FIL_Vista_Basic.o vista/FIL_Vista_File.o vista/FIL_Vista_Graph.o vista/FIL_Vista_Image.o');
end

% Step 2
% compile the mesh/vol functions
L = [];
L = add_mex_source(L,'.','vistaprimitive.cpp',{'GLNX86', 'GLNXA64'},[],'-Ivista -c -o vistaprimitive.o');
L = add_mex_source(L,'.','read_vista_mesh.cpp',{'GLNX86', 'GLNXA64'},[],'-Ivista vistaprimitive.o libvista.a');
L = add_mex_source(L,'.','write_vista_mesh.cpp',{'GLNX86', 'GLNXA64'},[],'-Ivista libvista.a');
L = add_mex_source(L,'.','write_vista_vol.cpp',{'GLNX86', 'GLNXA64'},[],'-Ivista libvista.a');

try
  compile_mex_list_vista(L, baseDir, force); 
catch
  me = lasterror;
  cd(oldDir);
  rethrow(me);
end
cd(oldDir);

function compile_mex_list_vista(L, baseDir, force)
% function compile_mex_list(L, baseDir)
%
% Compile a list of MEX files as determined by the input argument L.
% The second argument 'baseDir' is the common base directory for the
% files listed in L. The third argument is a flag that determines
% whether to force (re-)compilation even if the MEX file is up-to-date.
%
% See also ft_compile_mex, add_mex_source.

% (C) 2010 S. Klanke, 2011, Cristiano Micheli

for i=1:length(L)
   [~, name] = fileparts(L(i).relName);

   sfname = [baseDir filesep L(i).dir filesep L(i).relName ];
   SF = dir(sfname);
   if numel(SF)<1
      fprintf(1,'Error: source file %s cannot be found.', sfname);
      continue;
   end
   
   if ~force
      mfname = [baseDir filesep L(i).dir filesep name '.' mexext];
      MF = dir(mfname);
      if numel(MF)==1 && datenum(SF.date) <= datenum(MF.date)
         fprintf(1,'Skipping up-to-date MEX file %s/%s\n', L(i).dir, name);
         continue;
      end 
   end
   fprintf(1,'Compiling MEX file %s/%s ...\n', L(i).dir, name);
   cd([baseDir '/' L(i).dir]);
   cmd = sprintf('mex %s %s', L(i). relName, L(i).extras);
   eval(cmd);
end

function L = add_mex_source(L, directory, relName, matchPlatform, excludePlatform, extras)
% function L = add_mex_source(L, directory, relName, matchPlatform, excludePlatform, extras)
%
% Input + output argument L is a structure array of directory names, source file names,
% and extra arguments required for the compilation of MEX files. This function will
% create a new element of this structure and append it to L.
%
% Further inputs:
%   directory
%      target directory of the mex-file
%   relName
%      source file relative to 'directory'
%   matchPlatform
%      list of platforms this MEX file should only be compiled for.
%      use an empty matrix [] to compile for all platforms
%   excludePlatform
%      list of platforms this MEX file should NOT be compiled for.
%   extras
%      extra arguments to the MEX command, e.g. additional source files

% (C) 2010 S. Klanke

% Check if this file only needs compilation on certain platforms (including this one)
if nargin>3 && ~isempty(matchPlatform) 
   ok = false;
   for k=1:numel(matchPlatform)
       if strcmp(matchPlatform{k}, computer)
		  ok = true;
		  break;
	   end
	end
	if ~ok
	   return
	end
end

% Check if this file cannot be compiled on certain platforms (including this one)
if nargin>4 && ~isempty(excludePlatform) 
   ok = true;
   for k=1:numel(excludePlatform)
      if strcmp(excludePlatform{k}, computer)
         return;
      end
   end
end

L(end+1).dir   = directory;
L(end).relName = relName;
if nargin>5
  L(end).extras = extras;
end
