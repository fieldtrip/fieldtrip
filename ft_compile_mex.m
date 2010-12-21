function ft_compile_mex(force)
% FT_COMPILE_MEX  This script/function is used for compiling the various MEX files used in FieldTrip
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
% (C) 2010 S. Klanke

if nargin<1
   force=false;
end

%Possible COMPUTER types
%GLNX86
%GLNXA64
%PCWIN
%PCWIN64
%MAC
%MACI
%MACI64

L = [];
L = add_mex_source(L,'fileio/@uint64','abs');
L = add_mex_source(L,'fileio/@uint64','min');
L = add_mex_source(L,'fileio/@uint64','max');
L = add_mex_source(L,'fileio/@uint64','plus');
L = add_mex_source(L,'fileio/@uint64','minus');
L = add_mex_source(L,'fileio/@uint64','times');
L = add_mex_source(L,'fileio/@uint64','rdivide');

L = add_mex_source(L,'fileio/private','read_16bit');
L = add_mex_source(L,'fileio/private','read_24bit');
L = add_mex_source(L,'fileio/private','mxSerialize');
L = add_mex_source(L,'fileio/private','mxDeserialize');
L = add_mex_source(L,'fileio/private','read_ctf_shm', {'GLNX86'});  % only compile on GLNX86
L = add_mex_source(L,'fileio/private','write_ctf_shm', {'GLNX86'}); % only compile on GLNX86
L = add_mex_source(L,'fileio/private','../../realtime/datasource/siemens/sap2matlab',[],[],'../../realtime/datasource/siemens/siemensap.c -I../../realtime/datasource/siemens/');

L = add_mex_source(L,'forward/private','plgndr');
L = add_mex_source(L,'forward/private','meg_leadfield1');
L = add_mex_source(L,'forward/private','lmoutr',[],[],'geometry.c -I.');
L = add_mex_source(L,'forward/private','solid_angle',[],[],'geometry.c -I.');
L = add_mex_source(L,'forward/private','ptriproj',[],[],'geometry.c -I.');

L = add_mex_source(L,'private','splint_gh');
L = add_mex_source(L,'private','ltrisect',[],[],'geometry.c -I.');
L = add_mex_source(L,'private','routlm',[],[],'geometry.c -I.');
L = add_mex_source(L,'private','plinproj',[],[],'geometry.c -I.');

L = add_mex_source(L,'realtime/online_mri','ft_omri_smooth_volume');


oldDir = pwd;
[baseDir, myName] = fileparts(mfilename('fullpath'));
try
  compile_mex_list(L, baseDir, force);
catch me
  cd(oldDir);
  rethrow(me);
end
cd(oldDir);

