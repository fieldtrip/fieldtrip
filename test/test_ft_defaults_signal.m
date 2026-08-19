function test_ft_defaults_signal

% MEM 1gb
% WALLTIME 00:03:00
% DEPENDENCY ft_hastoolbox ft_defaults
% DATA no

% external/signal must stay on the path when the Signal Processing Toolbox is
% licensed but not installed, across repeated re-initialisation of ft_defaults

origpath = path;
cleanup  = onCleanup(@() path(origpath)); %#ok<NASGU>

% simulate a toolbox that is licensed but not installed, as on CI runners
sigroot = fullfile(matlabroot, 'toolbox', 'signal');
p = strsplit(path, pathsep);
p = p(~strncmpi(p, sigroot, numel(sigroot)));
path(strjoin(p, pathsep));

external_signal = fullfile(fileparts(which('ft_defaults')), 'external', 'signal');

for k = 1:4
  clear ft_defaults
  ft_defaults;

  if ~contains(path, external_signal)
    error('external/signal is not on the path after %d ft_defaults call(s)', k);
  end
  if isempty(which('butter'))
    error('butter is undefined after %d ft_defaults call(s)', k);
  end
end
