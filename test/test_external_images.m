function test_external_images

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY hsv2rgb

filelist = {'hsv2rgb'};

% NOTE: hsv2rgb was part of standard matlab in older matlab versions

[ftver, ftpath] = ft_version;
restoredefaultpath

% ensure that the the 'compat' (i.e. external/signal) is used
global ft_default
ft_default.toolbox.signal = 'compat';
addpath(ftpath);
ft_defaults;

for k = 1:numel(filelist)
  assert(exist(filelist{k}, 'file')==2);

  fprintf('testing the functionality of %s\n', filelist{k});
  switch filelist{k}
    case 'hsv2rgb'
      assert(isequal(hsv2rgb([0 1 1]),[1 0 0]));
      assert(isequal(hsv2rgb([1 0 1]),[1 1 1]));
    otherwise
      ft_error('function %s is not part of the official external/images directory', filelist{k});
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare the external/images output with the matlab version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = rand(10,10,3);
for k = 1:numel(filelist)
  try
    % this only works for the functions that create a window
    funhandle = str2func(filelist{k});
    output1.(filelist{k}) = funhandle(tmp);
  end
end

restoredefaultpath

% ensure that the the matlab version is used
ft_default.toolbox.imagesc = 'matlab';
addpath(ftpath);
ft_defaults;

for k = 1:numel(filelist)
  try
    funhandle = str2func(filelist{k});
    output2.(filelist{k}) = funhandle(tmp);
  end
end

fn1 = fieldnames(output1);
fn2 = fieldnames(output2);
assert(isequal(sort(fn1), sort(fn2)));
for k = 1:numel(fn1)
  assert(isalmostequal(output1.(fn1{k}), output2.(fn1{k}), 'abstol', 10*eps));
end
