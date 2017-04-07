function test_bug1114

% MEM 1500mb
% WALLTIME 00:10:00

% This function parses all FieldTrip main and module functions and determines
% whether there are any dependencies on fieldtrip/compat or any other
% compat directory. If so, the files are printed with an error.

[ftver, ftpath] = ft_version;

% ensure that the compat directories are on the path
addpath(fullfile(ftpath, 'compat'))
addpath(fullfile(ftpath, 'fileio/compat'))
addpath(fullfile(ftpath, 'forward/compat'))
addpath(fullfile(ftpath, 'plotting/compat'))
addpath(fullfile(ftpath, 'preproc/compat'))
addpath(fullfile(ftpath, 'utilities/compat'))

dirlist = {
  ftpath
  fullfile(ftpath, 'fileio')
  fullfile(ftpath, 'forward')
  fullfile(ftpath, 'inverse')
  fullfile(ftpath, 'plotting')
  fullfile(ftpath, 'connectivity')
  fullfile(ftpath, 'specest')
  fullfile(ftpath, 'trialfun')
  fullfile(ftpath, 'statfun')
  fullfile(ftpath, 'utilities')
  fullfile(ftpath, 'private')
  };

for dirindex=1:length(dirlist)
  functionlist = dir(fullfile(dirlist{dirindex}, '*.m'));
  functionlist = {functionlist.name};
  
  fprintf('==== processing directory %s ====\n', dirlist{dirindex});
  
  % find the dependencies
  [outlist, depmat] = mydepfun(functionlist);
  
  compat = false(size(outlist));
  for i=1:length(outlist)
    compat(i) = ~isempty(regexp(outlist{i}, '/compat', 'once'));
  end
  % switch to list indices
  compat = find(compat);
  
  if ~isempty(compat)
    % report on the problems
    fprintf('\nThe compat functions ...\n');
    disp(outlist(compat));
    fprintf('... are being used by\n');
    disp(functionlist(any(depmat(:,compat),2))');
    error('some of the FT functions depend on compat functions');
    % warning('some of the FT functions depend on compat functions');
  end
end % for dirlist

