function failed_neuromag_units

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_neuromag_units
% TEST ft_read_header ft_read_data

% http://bugzilla.fcdonders.nl/show_bug.cgi?id=953
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=963

% this is a test dataset from Rik Henson that contains both MEG and EEG
dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/run_01_raw.fif');

hdr = ft_read_header(dataset);

chanunit1 = hdr.chanunit; % this is in T, T/m and in V
chanunit2 = chanunit1;
chanunit2(strcmp(chanunit2, 'T/m')) = {'T/cm'};
chanunit3 = chanunit1;
chanunit3(strcmp(chanunit3, 'T')) = {'fT'};
chanunit4 = chanunit1;
chanunit4(strcmp(chanunit4, 'V')) = {'uV'};
chanunit5 = chanunit1;
chanunit5(strcmp(chanunit5, 'T')) = {'fT'}; chanunit5(strcmp(chanunit5, 'T/m')) = {'fT/cm'};

begsample = 1;
endsample = hdr.Fs;

dat0 = ft_read_data(dataset,                        'begsample', begsample, 'endsample', endsample);
dat1 = ft_read_data(dataset, 'chanunit', chanunit1, 'begsample', begsample, 'endsample', endsample);
dat2 = ft_read_data(dataset, 'chanunit', chanunit2, 'begsample', begsample, 'endsample', endsample);
dat3 = ft_read_data(dataset, 'chanunit', chanunit3, 'begsample', begsample, 'endsample', endsample);
dat4 = ft_read_data(dataset, 'chanunit', chanunit4, 'begsample', begsample, 'endsample', endsample);
dat5 = ft_read_data(dataset, 'chanunit', chanunit5, 'begsample', begsample, 'endsample', endsample);

% assert that the scaling is correct
% channel 1 and 2 are planar gradiometers, channel 3 is a magnetometer, channel 307 is eeg
assert_almost_equal     = @(x, y) assert(all(abs(x-y)/((x+y)/2)<1e3*eps));
assert_reasonably_equal = @(x, y) assert(all(abs(x-y)/((x+y)/2)<1e9*eps));

assert_almost_equal(dat1(1,1)./dat0(1,1), 1)
assert_almost_equal(dat1(3,1)./dat0(3,1), 1)
assert_almost_equal(dat1(307,1)./dat0(307,1), 1)

assert_almost_equal(dat2(1,1)./dat0(1,1), 0.01) % /cm instead /m
assert_almost_equal(dat2(3,1)./dat0(3,1), 1)
assert_almost_equal(dat2(307,1)./dat0(307,1), 1)

assert_almost_equal(dat3(1,1)./dat0(1,1), 1)
assert_almost_equal(dat3(3,1)./dat0(3,1), 1e15) % fT instead of T
assert_almost_equal(dat3(307,1)./dat0(307,1), 1)

assert_almost_equal(dat4(1,1)./dat0(1,1), 1)
assert_almost_equal(dat4(3,1)./dat0(3,1), 1)
assert_almost_equal(dat4(307,1)./dat0(307,1), 1e6) % uV instead of V

assert_almost_equal(dat5(1,1)./dat0(1,1), 1e13) % fT/cm instead of T/cm
assert_almost_equal(dat5(3,1)./dat0(3,1), 1e15) % fT instead of T
assert_almost_equal(dat5(307,1)./dat0(307,1), 1)

vol1 = [];
vol1.type = 'infinite';
vol1.unit = 'm';

vol2 = ft_convert_units(vol1, 'dm');
vol3 = ft_convert_units(vol1, 'cm');
vol4 = ft_convert_units(vol1, 'mm');

grad0 = hdr.grad; % % this is in T, T/m and in V
grad1 = ft_convert_units(grad0, 'm');
grad2 = ft_convert_units(grad0, 'dm');
grad3 = ft_convert_units(grad0, 'cm');
grad4 = ft_convert_units(grad0, 'mm');

assert_almost_equal(grad2.coilpos(:)\grad1.coilpos(:), 1e1\1); % note the use of the backslash
assert_almost_equal(grad2.chanpos(:)\grad1.chanpos(:), 1e1\1); % note the use of the backslash

assert_almost_equal(grad3.coilpos(:)\grad1.coilpos(:), 1e2\1); % note the use of the backslash
assert_almost_equal(grad3.chanpos(:)\grad1.chanpos(:), 1e2\1); % note the use of the backslash

assert_almost_equal(grad4.coilpos(:)\grad1.coilpos(:), 1e3\1); % note the use of the backslash
assert_almost_equal(grad4.chanpos(:)\grad1.chanpos(:), 1e3\1); % note the use of the backslash

pos1 = grad1.chanpos(1,:)/2; % in m
pos2 = grad2.chanpos(1,:)/2; % in dm
pos3 = grad3.chanpos(1,:)/2; % in cm
pos4 = grad4.chanpos(1,:)/2; % in mm

lf11 = ft_compute_leadfield(pos1, grad1, vol1); % it assumes that pos1 has the same units as grad1 and vol1
lf22 = ft_compute_leadfield(pos2, grad2, vol2); % it assumes that pos2 has the same units as grad3 and vol3
lf33 = ft_compute_leadfield(pos3, grad3, vol3); % it assumes that pos3 has the same units as grad3 and vol3
lf44 = ft_compute_leadfield(pos4, grad4, vol4); % it assumes that pos1 has the same units as grad1 and vol1

% assuming arbitrary units, the dipole is 10x further away with each step
assert_almost_equal(lf22(1,1)./lf11(1,1), 1e-3); % it is   10 dm away instead of 1 m
assert_almost_equal(lf33(1,1)./lf11(1,1), 1e-6); % it is  100 cm away instead of 1 m
assert_almost_equal(lf44(1,1)./lf11(1,1), 1e-9); % it is 1000 mm away instead of 1 m

distance      = {'m' 'dm' 'cm' 'mm'};
amplitude = {'T' 'mT' 'uT' 'nT' 'pT' 'fT'};
for i=1:4
  for j=1:6
    fprintf('distance = %s, amplitude = %s\n', distance{i}, amplitude{j});
    grad1a = ft_convert_grad(grad1, amplitude{j}, distance{i}, 'amplitude');
    grad2a = ft_convert_grad(grad2, amplitude{j}, distance{i}, 'amplitude');
    grad3a = ft_convert_grad(grad3, amplitude{j}, distance{i}, 'amplitude');
    grad4a = ft_convert_grad(grad4, amplitude{j}, distance{i}, 'amplitude');
    grad1b = ft_convert_grad(grad1, amplitude{j}, distance{i}, 'amplitude/distance');
    grad2b = ft_convert_grad(grad2, amplitude{j}, distance{i}, 'amplitude/distance');
    grad3b = ft_convert_grad(grad3, amplitude{j}, distance{i}, 'amplitude/distance');
    grad4b = ft_convert_grad(grad4, amplitude{j}, distance{i}, 'amplitude/distance');
    grad1c = ft_convert_grad(grad1b, amplitude{j}, distance{i}, 'amplitude'); % should again be the same as "a"
    grad2c = ft_convert_grad(grad2b, amplitude{j}, distance{i}, 'amplitude'); % should again be the same as "a"
    grad3c = ft_convert_grad(grad3b, amplitude{j}, distance{i}, 'amplitude'); % should again be the same as "a"
    grad4c = ft_convert_grad(grad4b, amplitude{j}, distance{i}, 'amplitude'); % should again be the same as "a"
    assert(isalmostequal(grad1c, grad1a, 'reltol', 100*eps));
    assert(isalmostequal(grad2c, grad2a, 'reltol', 100*eps));
    assert(isalmostequal(grad3c, grad3a, 'reltol', 100*eps));
    assert(isalmostequal(grad4c, grad4a, 'reltol', 100*eps));
  end
end

try
  haserror = true;
  ft_compute_leadfield(pos1, grad1, vol2); % should give an error
  haserror = false;
end % try
if ~haserror
  error('the invalid input should have resulted in an error')
end

try
  haserror = true;
  ft_compute_leadfield(pos1, grad1, vol3); % should give an error
  haserror = false;
end % try
if ~haserror
  error('the invalid input should have resulted in an error')
end

lf11 = ft_compute_leadfield(pos1, grad1, vol1); % these are in T and m
lfsi = ft_compute_leadfield(pos1, grad1, vol1, 'chanunit', grad1.chanunit);
assert_almost_equal(lf11(1,1)./lfsi(1,1), 1);

grad1a = ft_convert_grad(grad1, 'fT',  'm', 'amplitude');           vol1a = vol1;
grad1b = ft_convert_grad(grad1, 'fT', 'cm', 'amplitude');           vol1b = ft_convert_units(vol1, 'cm');
grad1c = ft_convert_grad(grad1, 'fT',  'm', 'amplitude/distance');  vol1c = vol1;
grad1d = ft_convert_grad(grad1, 'fT', 'cm', 'amplitude/distance');  vol1d = ft_convert_units(vol1, 'cm');
lf1a = ft_compute_leadfield(pos1    , grad1a, vol1a);
lf1b = ft_compute_leadfield(pos1*100, grad1b, vol1b); % in cm
lf1c = ft_compute_leadfield(pos1    , grad1c, vol1c);
lf1d = ft_compute_leadfield(pos1*100, grad1d, vol1d); % in cm

assert_almost_equal(lf1a(1,1)./lfsi(1,1), 1e15);
assert_almost_equal(lf1c(3,1)./lfsi(3,1), 1e15);

assert_almost_equal(lf1b(1,1)./lfsi(1,1), 1e15/(100^3));
assert_almost_equal(lf1b(3,1)./lfsi(3,1), 1e15/(100^3));

% assert_almost_equal(lf1c(1,1)./lfsi(1,1), 1e15); % this does not hold, as it is scaled with the baseline
assert_almost_equal(lf1c(3,1)./lfsi(3,1), 1e15);

% assert_almost_equal(lf1d(1,1)./lfsi(1,1), 1e15/(100^3)); % this does not hold, as it is scaled with the baseline
assert_almost_equal(lf1d(3,1)./lfsi(3,1), 1e15/(100^3));

% these should differ with one-over-the-baseline
% note that the numerical errors accumulate
assert_reasonably_equal(lf1c(1,1)./lf1a(1,1), 1/0.0168); % baseline is 0.0168 m
assert_reasonably_equal(lf1d(1,1)./lf1b(1,1), 1/1.6800); % baseline is 1.6800 cm

