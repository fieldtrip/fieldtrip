function test_external_signal

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY hanning

filelist = {'barthannwin'
            'blackmanharris'
            'bohmanwin'
            'boxcar'
            'butter'
            'filtfilt'
            'flattopwin'
            'gausswin'
            'hann'
            'hanning'
            'hilbert'
            'kaiser'
            'nuttallwin'
            'parzenwin'
            'rectwin'
            'tukeywin'
            'triang'
            'window'};

[ftver, ftpath] = ft_version;
restoredefaultpath

% ensure that the the 'compat' (i.e. external/signal) is used
global ft_default
ft_default.toolbox.signal = 'compat';
addpath(ftpath);
ft_defaults;

for k = 1:numel(filelist)
  assert(exist(filelist{k}, 'file'));

  fprintf('testing the functionality of %s\n', filelist{k});
  switch filelist{k}
    case 'barthannwin'
      assert(isequal(barthannwin(1), 1));
      assert(isequal(barthannwin (2), zeros(2, 1)));
    case 'blackmanharris'
      assert(isequal(blackmanharris(1),  1));
      assert(isalmostequal(blackmanharris(2),  0.00006 * ones(2, 1), 'abstol', eps));
      assert(isalmostequal(blackmanharris(15), flipud(blackmanharris(15)), 'abstol', 10*eps));
      assert(isalmostequal(blackmanharris(16), flipud(blackmanharris(16)), 'abstol', 10*eps));
      assert(isequal(blackmanharris(15), blackmanharris(15, 'symmetric')));
      tmp = blackmanharris(16);
      assert(isequal(tmp(1:15), blackmanharris(15, 'periodic')));
    case 'bohmanwin'
      assert(isequal(bohmanwin(1), 1));
      assert(isequal(bohmanwin(2), zeros(2, 1)));
    case 'boxcar'
      assert(isequal(boxcar(1), 1));
      assert(isequal(boxcar(2), ones(2, 1)));
      assert(isequal(boxcar(100), ones(100, 1)));
    case 'butter'
      % FIXME the octave tests are either too detailed or only check for
      % output orientation
    case 'filtfilt'
      r = randn(1,200);
      [b,a] = butter(10, [.2, .25]);
      yfb = filtfilt(b, a, r);
      assert(isequal(size(r), size(yfb)));
      assert(mean(abs(yfb)) < 1e3);
      assert(mean(abs(yfb)) < mean(abs(r)));
      ybf = fliplr(filtfilt(b, a, fliplr(r)));
      assert(mean(abs(ybf)) < 1e3);
      assert(mean(abs(ybf)) < mean(abs(r)));

      r = randn(1,1000);
      s = 10 * sin(pi * 4e-2 * (1:length(r)));
      [b,a] = butter(2, [4e-4 8e-2]);
      y = filtfilt(b, a, [r.' s.']);
      yr = filtfilt(b, a, r);
      ys = filtfilt(b, a, s);
      assert(isequal(y, [yr.' ys.']));
      y2 = filtfilt(b.', a.', [r.' s.']);
      assert(isequal(y, y2));
    case 'flattopwin'
      assert(isequal(flattopwin(1), 1));
      assert(isalmostequal(flattopwin(2), -0.421051e-3 * ones (2, 1), 'abstol', eps)); % slightly different coeffs.
      assert(isalmostequal(flattopwin(15), flipud(flattopwin(15)), 'abstol', 10*eps));
      assert(isalmostequal(flattopwin(16), flipud(flattopwin(16)), 'abstol', 10*eps));
      assert(isequal(flattopwin(15), flattopwin(15, 'symmetric')));
      tmp = flattopwin(16);
      assert(isequal(tmp(1:15), flattopwin(15, 'periodic')));
    case 'gausswin'
      assert(isequal(gausswin(1), 1));
      assert(isequal(gausswin(2), [exp(-3.125); exp(-3.125)]));
      assert(isequal(gausswin(3), [exp(-3.125); 1; exp(-3.125)]));
    case 'hann'
      assert(isequal(hann(1), 1));
      assert(isalmostequal(hann(2), [0;0], 'abstol', eps)); 
      assert(isalmostequal(hann(16), flipud(hann(16)), 'abstol', 10*eps));
      assert(isalmostequal(hann(15), flipud(hann(15)), 'abstol', 10*eps));
      assert(isequal(hann(15), hann(15, 'symmetric')));
      tmp = hann(16, 'symmetric');
      assert(isequal(tmp(1:15), hann(15, 'periodic')));
    case 'hanning'
      % this is more or less tested with 'hann'
    case 'hilbert'
      % FIXME
    case 'kaiser'
      assert(isequal(kaiser(1), 1));
    case 'nuttallwin'
      assert(isequal(nuttallwin(1), 1));
      assert(isalmostequal(nuttallwin(2), 0.3628e-3*ones(2,1), 'abstol', eps));
      assert(isalmostequal(nuttallwin(15), flipud(nuttallwin(15)), 'abstol', 10*eps));
      assert(isalmostequal(nuttallwin(16), flipud(nuttallwin(16)), 'abstol', 10*eps));
      assert(isequal(nuttallwin(15), nuttallwin(15, 'symmetric')));
      tmp = nuttallwin(16);
      assert(isequal(tmp(1:15), nuttallwin(15, 'periodic')));
    case 'parzenwin'
      assert(isequal(parzenwin(1), 1));
      assert(isequal(parzenwin(2), 0.25 * ones (2, 1)));
    case 'rectwin'
      assert(isequal(rectwin(1),   1));
      assert(isequal(rectwin(2),   ones(2, 1)));
      assert(isequal(rectwin(100), ones(100, 1)));
    case 'triang'
      assert(isequal(triang(1), 1));
      assert(isequal(triang(2), [1; 1]/2));
      assert(isequal(triang(3), [1; 2; 1]/2));
      assert(isequal(triang(4), [1; 3; 3; 1]/4));
    case 'tukeywin'
      assert(isequal(tukeywin(1), 1));
      assert(isequal(tukeywin(2), zeros (2, 1)));
      assert(isequal(tukeywin(3), [0; 1; 0]));
      assert(isequal(tukeywin(16, 0), rectwin (16)));
      assert(isequal(tukeywin(16, 1), hanning (16)));
    case 'window'
      assert(isequal(window(@hanning, 16), window('hanning', 16)));
      assert(isequal(window(@triang,  16), window('triang', 16)));
    otherwise
      ft_error('function %s is not part of the official external/signal directory', filelist{k});
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare the external/signal output with the matlab version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:numel(filelist)
  try
    % this only works for the functions that create a window
    funhandle = str2func(filelist{k});
    output1.(filelist{k}) = funhandle(50);
  end
end

restoredefaultpath

% ensure that the the matlab version is used
ft_default.toolbox.signal = 'matlab';
addpath(ftpath);
ft_defaults;

for k = 1:numel(filelist)
  try
    % this only works for the functions that create a window
    funhandle = str2func(filelist{k});
    output2.(filelist{k}) = funhandle(50);
  end
end

fn1 = fieldnames(output1);
fn2 = fieldnames(output2);
assert(isequal(sort(fn1), sort(fn2)));
for k = 1:numel(fn1)
  assert(isalmostequal(output1.(fn1{k}), output2.(fn1{k}), 'abstol', 10*eps));
end
