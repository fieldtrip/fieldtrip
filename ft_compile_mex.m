function ft_compile_mex(force)

% FT_COMPILE_MEX is used for compiling most of the MEX files that are used in FieldTrip
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

% Copyright (C) 2010, Stefan Klanke
%
% $Log$

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
L = add_mex_source(L,'@config/private','reset');

L = add_mex_source(L,'src','ft_getopt');
L = add_mex_source(L,'src','keyval');
L = add_mex_source(L,'src','read_16bit');
L = add_mex_source(L,'src','read_24bit');
L = add_mex_source(L,'src','mxSerialize');
L = add_mex_source(L,'src','mxDeserialize');
L = add_mex_source(L,'src','read_ctf_shm', {'GLNX86'});  % only compile on GLNX86
L = add_mex_source(L,'src','write_ctf_shm', {'GLNX86'}); % only compile on GLNX86

L = add_mex_source(L,'src','lmoutr',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','ltrisect',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','plinproj',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','ptriproj',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','routlm',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','solid_angle',[],[],'geometry.c -I.');
L = add_mex_source(L,'src','rfbevent',[],[],'d3des.c -I.');
L = add_mex_source(L,'src','meg_leadfield1');
L = add_mex_source(L,'src','splint_gh');
L = add_mex_source(L,'src','plgndr');
L = add_mex_source(L,'src','ft_spike_sub_crossx');

L = add_mex_source(L,'realtime/online_mri','ft_omri_smooth_volume');
L = add_mex_source(L,'realtime/acquisition/siemens', 'sap2matlab',[],[],'siemensap.c -I.');


oldDir = pwd;
[baseDir, myName] = fileparts(mfilename('fullpath'));
try
  compile_mex_list(L, baseDir, force);
catch
  % the "catch me" syntax is broken on MATLAB74, this fixes it
  me = lasterror;
  cd(oldDir);
  rethrow(me);
end
cd(oldDir);

