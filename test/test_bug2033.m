% function test_bug2033

% TEST test_bug2033

type = {
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

