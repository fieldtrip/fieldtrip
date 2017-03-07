function test_bug2570

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_apply_montage ft_scalingfactor

montage = [];
montage.tra = 1e6;
montage.labelold = {'Cz'};
montage.labelnew = {'Cz'};
montage.chanunitold = {'V'};
montage.chanunitnew = {'uV'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = ft_apply_montage(montage, montage);
assert(isequal(output, montage)); % it should not have changed

output = ft_apply_montage(montage, montage, 'inverse', 1);
assert(isequal(output.tra, 1)); % this should be an identity transform
assert(isequal(output.chanunitold, output.chanunitnew));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_V = [];
data_V.label = {'Cz'};
data_V.time = {1:10};
data_V.trial = {ones(1,10)};
data_V.chanunit = {'V'};

data_mV = [];
data_mV.label = {'Cz'};
data_mV.time = {1:10};
data_mV.trial = {ones(1,10)*1000};
data_mV.chanunit = {'mV'};

data_uV = [];
data_uV.label = {'Cz'};
data_uV.time = {1:10};
data_uV.trial = {ones(1,10)*1000000};
data_uV.chanunit = {'uV'};

output1 = ft_apply_montage(data_V,  montage);
output2 = ft_apply_montage(data_mV, montage);
output3 = ft_apply_montage(data_uV, montage);

assert(isequal(output1, output2));
assert(isequal(output1, output3));
assert(isequal(output1, data_uV)); % they should all be in uV

output1 = ft_apply_montage(data_V,  montage, 'inverse', 'yes');
output2 = ft_apply_montage(data_mV, montage, 'inverse', 'yes');
output3 = ft_apply_montage(data_uV, montage, 'inverse', 'yes');

assert(isequal(output1, output2));
assert(isequal(output1, output3));
assert(isequal(output1, data_V)); % they should all be in V

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_T = [];
data_T.label = {'Cz'};
data_T.time = {1:10};
data_T.trial = {ones(1,10)};
data_T.chanunit = {'T'};

try
  ft_apply_montage(data_T,  montage);
  ok = false;
catch
  ok = true;
end
if ~ok
  error('this should have given an error, T cannot be converted to V');
end

