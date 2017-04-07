function test_bug2033

% WALLTIME 00:10:00
% MEM 1gb


type = {
  'boxcar'
  'barthannwin'
  'blackmanharris'
  'bohmanwin'
  % 'chebwin'
  'flattopwin'
  'gausswin'
  'hann'
  'kaiser'
  'nuttallwin'
  'parzenwin'
  'rectwin'
  'triang'
  'tukeywin'
  };

for i=1:length(type)
  w = window(type{i}, 1000);
end

[a, b] = butter(10, 0.5);
h = hilbert(randn(1,10));

