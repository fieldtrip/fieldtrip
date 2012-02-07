% function test_bugXXX

% This function parses all "fieldtrip/ft_*.m" main functions and determines
% whether there are any dependencies on fieldtrip/compat or any other
% compat directory. If so, the files are printed with an error.

% TODO
%  crimic should be explained how to use compat
%  fixed a spike function (channelselection), should be committed
%  fieldtrip/compat/openmeeg.m should be removed
%  fieldtrip/ft_headmodelplot.m should be removed
%  compat/ft_prepare_bemmodel.m and ft_prepare_bemmodel.m should be merged, the compat one should then be removed

fieldtripdir = fileparts(which('ft_defaults'));

% ensure that the compat directories are on the path
addpath(fullfile(fieldtripdir, 'compat'))
addpath(fullfile(fieldtripdir, 'fileio/compat'))
addpath(fullfile(fieldtripdir, 'forward/compat'))
addpath(fullfile(fieldtripdir, 'plotting/compat'))
addpath(fullfile(fieldtripdir, 'preproc/compat'))
addpath(fullfile(fieldtripdir, 'utilities/compat'))

mainfunctionlist = dir(fullfile(fieldtripdir, 'ft_*.m'));
mainfunctionlist = {mainfunctionlist.name};

% find the dependencies
[outlist, depmat] = mydepfun(mainfunctionlist);

compat = false(size(outlist));
for i=1:length(outlist)
  compat(i) = ~isempty(regexp(outlist{i}, '/compat', 'once'));
end
% switch to list indices
compat = find(compat);

if ~isempty(compat)
  % report on the problems
  disp(outlist(compat));
  disp(mainfunctionlist(any(depmat(:,compat),2))');
  error('some of the FT main functions depend on compat functions');
end
