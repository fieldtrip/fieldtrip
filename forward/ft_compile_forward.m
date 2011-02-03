function ft_compile_forward(force)
% FT_COMPILE_FORWARD  This script/function is used for compiling the various MEX files used in FieldTrip 'forward'
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
L = add_mex_source(L,'private','plgndr');
L = add_mex_source(L,'private','meg_leadfield1');
L = add_mex_source(L,'private','lmoutr',[],[],'geometry.c -I.');
L = add_mex_source(L,'private','solid_angle',[],[],'geometry.c -I.');
L = add_mex_source(L,'private','ptriproj',[],[],'geometry.c -I.');

oldDir = pwd;
[baseDir, myName] = fileparts(mfilename('fullpath'));
try
  compile_mex_list(L, baseDir, force);
catch me
  cd(oldDir);
  rethrow(me);
end
cd(oldDir);


