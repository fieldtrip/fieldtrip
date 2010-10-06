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
L = addSource(L,'fileio/@uint64','abs');
L = addSource(L,'fileio/@uint64','min');
L = addSource(L,'fileio/@uint64','max');
L = addSource(L,'fileio/@uint64','plus');
L = addSource(L,'fileio/@uint64','minus');
L = addSource(L,'fileio/@uint64','times');
L = addSource(L,'fileio/@uint64','rdivide');

L = addSource(L,'fileio/private','read_16bit');
L = addSource(L,'fileio/private','read_24bit');
L = addSource(L,'fileio/private','mxSerialize');
L = addSource(L,'fileio/private','mxDeserialize');
L = addSource(L,'fileio/private','read_ctf_shm', {'GLNX86'});  % only compile on GLNX86
L = addSource(L,'fileio/private','write_ctf_shm', {'GLNX86'}); % only compile on GLNX86
L = addSource(L,'fileio/private','../../realtime/datasource/siemens/sap2matlab',[],[],'../../realtime/datasource/siemens/siemensap.c -I../../realtime/datasource/siemens/');

L = addSource(L,'forward/private','plgndr');
L = addSource(L,'forward/private','meg_leadfield1');
L = addSource(L,'forward/private','lmoutr',[],[],'geometry.c -I.');
L = addSource(L,'forward/private','solid_angle',[],[],'geometry.c -I.');
L = addSource(L,'forward/private','ptriproj',[],[],'geometry.c -I.');

L = addSource(L,'private','splint_gh');
L = addSource(L,'private','ltrisect',[],[],'geometry.c -I.');
L = addSource(L,'private','routlm',[],[],'geometry.c -I.');
L = addSource(L,'private','plinproj',[],[],'geometry.c -I.');

L = addSource(L,'realtime/online_mri','ft_omri_smooth_volume');


oldDir = pwd;
[baseDir, myName] = fileparts(mfilename('fullpath'));

for i=1:length(L)
   [relDir, name] = fileparts(L(i).relName);

   sfname = [baseDir filesep L(i).dir filesep L(i).relName '.c'];
   SF = dir(sfname);
   if numel(SF)<1
      fprintf(1,'Error: source file %s cannot be found.', sfname);
      continue;
   end

   if ~force
      mfname = [baseDir filesep L(i).dir filesep name '.' mexext];
      MF = dir(mfname);
      if numel(MF)==1 && SF.datenum <= MF.datenum
         fprintf(1,'Skipping up-to-date MEX file %s/%s\n', L(i).dir, name);
         continue;
      end 
   end
   fprintf(1,'Compiling MEX file %s/%s ...\n', L(i).dir, name);
   cd([baseDir '/' L(i).dir]);
   cmd = sprintf('mex %s.c %s', L(i). relName, L(i).extras);
   eval(cmd);
end

cd(oldDir);



function L = addSource(L, directory, relName, matchPlatform, excludePlatform, extras)

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
