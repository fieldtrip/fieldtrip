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
items = [items mex_item('@config/private','increment')];
items = [items mex_item('@config/private','reset')];
items = [items mex_item('fileio/@uint64','abs')];
items = [items mex_item('fileio/@uint64','max')];
items = [items mex_item('fileio/@uint64','min')];
items = [items mex_item('fileio/@uint64','minus')];
items = [items mex_item('fileio/@uint64','plus')];
items = [items mex_item('fileio/@uint64','rdivide')];
items = [items mex_item('fileio/@uint64','times')];
items = [items mex_item('src','det2x2')];
items = [items mex_item('src','ft_getopt')];
items = [items mex_item('src','ft_spike_sub_crossx')];
items = [items mex_item('src','getpid')];
items = [items mex_item('src','inv2x2')];
items = [items mex_item('src','meg_leadfield1')];
items = [items mex_item('src','mtimes2x2')];
items = [items mex_item('src','mxDeserialize')];
items = [items mex_item('src','mxSerialize')];
items = [items mex_item('src','nanmean')];
items = [items mex_item('src','nanstd')];
items = [items mex_item('src','nanvar')];
items = [items mex_item('src','nansum')];
items = [items mex_item('src','plgndr')];
items = [items mex_item('src','read_16bit')];
items = [items mex_item('src','read_24bit')];
items = [items mex_item('src','rename')]; 
items = [items mex_item('src','sandwich2x2')];
items = [items mex_item('src','splint_gh')];
items = [items mex_item('src', 'lmoutr'     , 'depends', {'geometry.c'}, 'mexFlags' , '-I.')];
items = [items mex_item('src', 'ltrisect'   , 'depends', {'geometry.c'}, 'mexFlags' , '-I.')];
items = [items mex_item('src', 'plinproj'   , 'depends', {'geometry.c'}, 'mexFlags' , '-I.')];
items = [items mex_item('src', 'ptriproj'   , 'depends', {'geometry.c'}, 'mexFlags' , '-I.')];
items = [items mex_item('src', 'routlm'     , 'depends', {'geometry.c'}, 'mexFlags' , ' -I.')];
items = [items mex_item('src', 'solid_angle', 'depends', {'geometry.c'}, 'mexFlags' , '-I.')];
items = [items mex_item('src','rfbevent', 'excludePlatform',{'win32', 'win64'}, 'extras', 'd3des.c -I.')];
items = [items mex_item('src','read_ctf_shm', 'matchPlatform', {'glnx86'})];
items = [items mex_item('src','write_ctf_shm', 'matchPlatform', {'glnx86'})];


function item = mex_item(directory, relName, varargin)
p = inputParser;
p.addRequired('directory'         , @ischar);
p.addRequired('relName'           , @ischar);
p.addParamValue('matchPlatform'   , []        , @iscellstr);
p.addParamValue('excludePlatform' , []        , @iscellstr);
p.addParamValue('extras'          , ''        , @ischar);
p.addParamValue('mexFlags'        , ''        , @ischar);
p.addParamValue('depends'         , []        , @iscellstr);
p.parse(directory                 , relName   , varargin{:});

MEXEXT = mexext('all');  % FIXME: probably we can cache this func-call.
for i = 1:length(MEXEXT)

  % Test if dependency needs to be generated for this platform:
  if ~isempty(p.Results.matchPlatform)
    if ~any(strncmp(MEXEXT(i).arch, p.Results.matchPlatform, ...
      length(MEXEXT(i).arch)))
      continue
    end
  end
  if ~isempty(p.Results.excludePlatform)
    if any(strncmp(MEXEXT(i).arch, p.Results.excludePlatform, ...
      length(MEXEXT(i).arch)))
      continue
    end
  end

  item(i).mex_flags  = p.Results.extras;
  item(i).arch       = MEXEXT(i).arch;
  item(i).target     = [p.Results.directory, '/', p.Results.relName, '.', MEXEXT(i).arch];
  item(i).source     = [[p.Results.directory, '/', p.Results.relName, '.c'], ...
    p.Results.depends];
  
end
