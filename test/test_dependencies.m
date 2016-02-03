function [funlist, deplist, depmat] = test_dependencies

% WALLTIME 0:20:00
% MEM 1gb

% TEST_DEPENDENCIES checks the dependencies on the backward compatibility functions
% and on the external toolboxes.

ft_defaults

p = fileparts(which(mfilename));

f1 = dir(fullfile(p, '*.m'));
f1 = {f1.name}';

f2 = dir(fullfile(p, 'utilities', '*.m'));
f2 = {f2.name}';

f3 = dir(fullfile(p, 'preproc', '*.m'));
f3 = {f3.name}';

f4 = dir(fullfile(p, 'fileio', '*.m'));
f4 = {f4.name}';

f5 = dir(fullfile(p, 'forward', '*.m'));
f5 = {f5.name}';

f6 = dir(fullfile(p, 'inverse', '*.m'));
f6 = {f6.name}';

f7 = dir(fullfile(p, 'realtime', '*.m'));
f7 = {f7.name}';

f8 = dir(fullfile(p, 'realtime', 'datasource', '*.m'));
f8 = {f8.name}';

f9 = dir(fullfile(p, 'peer', '*.m'));
f9 = {f9.name}';

f10 = dir(fullfile(p, 'plotting', '*.m'));
f10 = {f10.name}';

f11 = dir(fullfile(p, 'statfun', '*.m'));
f11 = {f11.name}';

f12 = dir(fullfile(p, 'specest', '*.m'));
f12 = {f12.name}';

f13 = dir(fullfile(p, 'connectivity', '*.m'));
f13 = {f13.name}';

f14 = dir(fullfile(p, 'contrib', '*.m'));
f14 = {f14.name}';

f15 = dir(fullfile(p, 'qsub', '*.m'));
f15 = {f15.name}';

f16 = dir(fullfile(p, 'contrib', 'spike', '*.m'));
f16 = {f16.name}';

f17 = dir(fullfile(p, 'contrib', 'trentool', '*.m'));
f17 = {f17.name}';

f18 = dir(fullfile(p, 'utility', '*.m'));
f18 = {f18.name}';

funlist = cat(1, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18);

for i=1:length(funlist)
  [p, funlist{i}, x] = fileparts(funlist{i});
end

% determine the dependencies, toponly, including matlab
[deplist, depmat] = mydepfun(funlist, true, true);

compatfun = ~cellfun(@isempty, regexp(deplist, [filesep 'compat' filesep]));
if any(compatfun)
  fprintf('\nThere is a dependency from function "%s" onto a compat function\n', funlist{any(depmat(:,compatfun),2)});
  fprintf('\nThere is a dependency onto the compat function "%s"\n', deplist{compatfun});
  error('dependencies on compat functions are not allowed');
end

toolboxfun = find(~cellfun(@isempty, regexp(deplist, matlabroot)));
toolbox = cell(size(toolboxfun));
t0 = length(tokenize(matlabroot, filesep));
for i=1:length(toolboxfun)
  t = tokenize(deplist{toolboxfun(i)}, filesep);
  toolbox{i} = t{t0+2};
end
toolbox = unique(toolbox);
toolbox = setdiff(toolbox, {'matlab'});
for i=1:length(toolbox)
  fprintf('\nThe following functions depend on the Mathworks "%s" toolbox:\n', toolbox{i});
  toolboxfun = ~cellfun(@isempty, regexp(deplist, [matlabroot filesep 'toolbox' filesep toolbox{i}]));
  fprintf('\t%s\n', funlist{any(depmat(:,toolboxfun),2)});
end

if ~nargout
  % do not return anything
  clear all
end

