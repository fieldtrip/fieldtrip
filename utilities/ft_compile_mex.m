function ft_compile_mex(force)

% FT_COMPILE_MEX can be used for compiling most of the FieldTrip MEX files
% Note that this function does not put the MEX files in the correct
% location in the private folders, this is managed by SVN autosync. In case
% you are not working with SVN and you want to recompile the mex files for
% your platform, you can find all mex files for your platform and move them
% to a backup directory that is not on your MATLAB path. Subsequently you
% can rtun this function to recompile it on your platform with your
% compiler settings
%
% The standards procedure for compiling mex files is detailled on
% http://fieldtriptoolbox.org/development/guidelines/code#compiling_mex_files
%
% Please note that this script does NOT set up your MEX environment for
% you, so in case you haven't selected the C compiler on Windows yet, you
% need to type 'mex -setup' first to choose either the LCC, Borland or
% Microsoft compiler. If you want to use MinGW, you also need to install
% Gnumex (http://gnumex.sourceforget.net), which comes with its own
% procedure for setting up the MEX environment.
%
% The logic in this script is to first build a list of files that actually
% need compilation for the particular platform that MATLAB is running on,
% and then to go through that list. Functions are added to the list by
% giving their destination directory and (relative to that) the name of the
% source file (without the .c). Optionally, you can specify a list of
% platform this file needs to be compiled on only, and a list of platforms
% where you don't compile it on. Finally, you can give extra arguments to
% the MEX command, e.g., for including other c-sources or giving compiler
% flags.
%
% See also MEX

% Copyright (C) 2010, Stefan Klanke
%
% $Id$

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

[ftver, ftpath] = ft_version;

L = [];
L = add_mex_source(L,'fileio/@uint64','abs');
L = add_mex_source(L,'fileio/@uint64','min');
L = add_mex_source(L,'fileio/@uint64','max');
L = add_mex_source(L,'fileio/@uint64','plus');
L = add_mex_source(L,'fileio/@uint64','minus');
L = add_mex_source(L,'fileio/@uint64','times');
L = add_mex_source(L,'fileio/@uint64','rdivide');

L = add_mex_source(L,'@config/private','deepcopy');
L = add_mex_source(L,'@config/private','increment');
L = add_mex_source(L,'@config/private','setzero');

L = add_mex_source(L,'realtime/online_mri','ft_omri_smooth_volume');
L = add_mex_source(L,'realtime/src/acquisition/siemens/src', 'sap2matlab', [], [], 'siemensap.c -I../include');

L = add_mex_source(L,'src','ft_getopt');
L = add_mex_source(L,'src','read_16bit');
L = add_mex_source(L,'src','read_24bit');
L = add_mex_source(L,'src','read_ctf_shm', {'GLNX86'});  % only compile on GLNX86
L = add_mex_source(L,'src','write_ctf_shm', {'GLNX86'}); % only compile on GLNX86

L = add_mex_source(L,'src','lmoutr',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','ltrisect',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','plinproj',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','ptriproj',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','routlm',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','solid_angle',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','rfbevent',[],{'PCWIN', 'PCWIN64'},'d3des.c -I.'); % do not compile on WIN32 and WIN64
L = add_mex_source(L,'src','meg_leadfield1');
L = add_mex_source(L,'src','splint_gh');
L = add_mex_source(L,'src','plgndr');

L = add_mex_source(L,'src','ft_spike_sub_crossx');
L = add_mex_source(L,'src','rename');
L = add_mex_source(L,'src','getpid');

L = add_mex_source(L,'src','nanmean');
L = add_mex_source(L,'src','nanstd');
L = add_mex_source(L,'src','nanvar');
L = add_mex_source(L,'src','nansum');
L = add_mex_source(L,'src','nanstd');
L = add_mex_source(L,'src','det2x2');
L = add_mex_source(L,'src','inv2x2');
L = add_mex_source(L,'src','mtimes2x2');
L = add_mex_source(L,'src','sandwich2x2');
L = add_mex_source(L,'src','combineClusters');

% this one is located elsewhere
L = add_mex_source(L,'external/fileexchange','CalcMD5',[],[],'CFLAGS=''-std=c99 -fPIC''');

% this one depends on the MATLAB version
if ft_platform_supports('libmx_c_interface')
  % use the C interface
  L = add_mex_source(L,'src','mxSerialize_c');
  L = add_mex_source(L,'src','mxDeserialize_c');
else
  % use the C++ interface
  L = add_mex_source(L,'src','mxSerialize_cpp');
  L = add_mex_source(L,'src','mxDeserialize_cpp');
end

oldDir = pwd;
try
  compile_mex_list(L, ftpath, force);
catch
  % the "catch me" syntax is broken on MATLAB74, this fixes it
  me = lasterror;
  cd(oldDir);
  rethrow(me);
end
cd(oldDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%
% Use as
%   list = add_mex_source(list, directory, relName, includePlatform, excludePlatform, extras)
%
% The list is a structure array of directory names, source file names, and
% extra arguments required for the compilation of MEX files. This function will
% create a new element of this structure and append it to L.
%
%   directory       = target directory of the mex-file
%   relName         = source file relative to 'directory'
%   includePlatform   = list of platforms this MEX file should only be compiled for.
%                     use an empty matrix [] to compile for all platforms
%   excludePlatform = list of platforms this MEX file should NOT be compiled for.
%   extras          = extra arguments to the MEX command, e.g. additional source files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = add_mex_source(L, directory, relName, includePlatform, excludePlatform, extras)

% Check if this file only needs compilation on certain platforms (including this one)
if nargin>3 && ~isempty(includePlatform)
  ok = false;
  for k=1:numel(includePlatform)
    if strcmp(includePlatform{k}, computer)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%
% Use as
%  compile_mex_list(L, baseDir)
%
% Compile a list of MEX files as determined by the input argument L.
% The second argument 'baseDir' is the common base directory for the
% files listed in L. The third argument is a flag that determines
% whether to force (re-)compilation even if the MEX file is up-to-date.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compile_mex_list(L, baseDir, force)

for i=1:length(L)
  [relDir, name] = fileparts(L(i).relName);
  sfname1 = fullfile(baseDir, L(i).dir, [L(i).relName '.c']);
  sfname2 = fullfile(baseDir, L(i).dir, [L(i).relName '.cpp']);
  if exist(sfname1, 'file')
    sfname = sfname1;
    L(i).ext = 'c';
  elseif exist(sfname2, 'file')
    sfname = sfname2;
    L(i).ext = 'cpp';
  else
    sfname = '';
    L(i).ext = '';
    fprintf(1,'Error: source file for %s cannot be found.\n', L(i).relName);
    continue;
  end
  SF = dir(sfname);
  
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
  cmd = sprintf('mex %s.%s %s', L(i).relName, L(i).ext, L(i).extras);
  eval(cmd);
end
