function items = ft_dependencies(cfg)

% FT_RECOMPILE supports FieldTrip developers in compilation by finding outdated
% targets (e.g. compiled mex files, executables), and by recompiling mex files.
% NOTE: WORK IN PROGRESS!
% 
% See also FT_COMPILE_MEX.
 
% Use as outdata = ft_recompile(cfg, indata) where indata is <<describe
% the type of data or where it comes from>> and cfg is a configuratioun
% structure that should contain:
%
%  cfg.option1 = value, explain the value here (default = something)
%  cfg.option2 = value, describe the value here and if needed continue here
%  to allow automatic parsing of the help
%
% Target are specified as follows:
%   item.target_dir      = target directory of the MEX file
%   item.relname         = source file relative to target_dir
%   item.platforms_white = list of platforms for this target. See @@.
%   item.platforms_black = list of blacklisted platforms for this target.
%   item.mex_arg = extra argument for the MEX command
%
% The configuration can optionally contain cfg.option3 = value, explain it
% here (default is automatic).
%
%
% See also FT_COMPILE_MEX.

% Copyright (c) 2010, Stefan Klanke 
%
% $Id: ft_compile_mex.m 5057 2011-12-20 14:03:21Z jansch $


% Mex files can be found using:
% `find . -name "*.c" | xargs grep -l mexFunction | sort`
% TODO: what should we do with external mex-files?

warning('This function is not finished yet!');

items        = mex_item('@config/private','deepcopy');
items(end+1) = mex_item('@config/private','increment');
items(end+1) = mex_item('@config/private','reset');
items(end+1) = mex_item('fileio/@uint64','abs');
items(end+1) = mex_item('fileio/@uint64','max');
items(end+1) = mex_item('fileio/@uint64','min');
items(end+1) = mex_item('fileio/@uint64','minus');
items(end+1) = mex_item('fileio/@uint64','plus');
items(end+1) = mex_item('fileio/@uint64','rdivide');
items(end+1) = mex_item('fileio/@uint64','times');
items(end+1) = mex_item('src','det2x2');
items(end+1) = mex_item('src','ft_getopt');
items(end+1) = mex_item('src','ft_spike_sub_crossx');
items(end+1) = mex_item('src','getpid');
items(end+1) = mex_item('src','inv2x2');
items(end+1) = mex_item('src','meg_leadfield1');
items(end+1) = mex_item('src','mtimes2x2');
items(end+1) = mex_item('src','mxDeserialize');
items(end+1) = mex_item('src','mxSerialize');
items(end+1) = mex_item('src','nanmean');
items(end+1) = mex_item('src','nanstd');
items(end+1) = mex_item('src','nanvar');
items(end+1) = mex_item('src','nansum');
items(end+1) = mex_item('src','plgndr');
items(end+1) = mex_item('src','read_16bit');
items(end+1) = mex_item('src','read_24bit');
items(end+1) = mex_item('src','rename');
items(end+1) = mex_item('src','sandwich2x2');
items(end+1) = mex_item('src','splint_gh');
items(end+1) = mex_item('src','lmoutr'  , 'extras', 'geometry.c -I.');
items(end+1) = mex_item('src','ltrisect', 'extras', 'geometry.c -I.');
items(end+1) = mex_item('src','plinproj', 'extras' ,'geometry.c -I.');
items(end+1) = mex_item('src','ptriproj', 'extras' ,'geometry.c -I.');
items(end+1) = mex_item('src','routlm', 'extras','geometry.c -I.');
items(end+1) = mex_item('src','solid_angle', 'extras','geometry.c -I.');
items(end+1) = mex_item('src','rfbevent', 'excludePlatform',{'PCWIN', 'PCWIN64'}, 'extras', 'd3des.c -I.');
items(end+1) = mex_item('src','read_ctf_shm', 'matchPlatform', {'GitemsNX86'});
items(end+1) = mex_item('src','write_ctf_shm', 'matchPlatform', {'GitemsNX86'});


function item = mex_item(directory, relName, varargin)
p = inputParser;
p.addRequired('directory'         , @ischar);
p.addRequired('relName'           , @ischar);
p.addParamValue('matchPlatform'   , []        , @iscellstr);
p.addParamValue('excludePlatform' , []        , @iscellstr);
p.addParamValue('extras'          , ''        , @ischar);
p.parse(directory, relName, varargin{:});

item.target_dir         = p.Results.directory;
item.rel_source         = p.Results.relName;
item.platform_whitelist = p.Results.matchPlatform;
item.platform_blacklist = p.Results.excludePlatform;
item.mex_flags          = p.Results.extras;
