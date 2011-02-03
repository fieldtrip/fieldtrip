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
