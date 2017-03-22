function test_bug1916

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_scalingfactor

if isempty(which('ft_scalingfactor'))
  % it is in the fieldtrip/private directory
  [ftver, ftpath] = ft_version;
  cd(fullfile(ftpath, 'private'));
end

old = 'mT/mm';
new = 'fT/km';

% the first time for any conversion takes 150 ms
stopwatch = tic;
factor = ft_scalingfactor(old, new);
toc(stopwatch)

n = 300;
old = repmat({old}, 1, n);
new = repmat({new}, 1, n);

% this used to work fine, because all inputs are the same
stopwatch = tic;
factor = cellfun(@ft_scalingfactor, old, new);
toc(stopwatch)

n = 100;
old = repmat({'T/cm' 'T/cm' 'T'}, 1, n);
new = repmat({'fT/m' 'fT/m' 'fT'}, 1, n);

% this used to be very slow, because of the A-A-B pattern
stopwatch = tic;
factor = cellfun(@ft_scalingfactor, old, new);
toc(stopwatch)


