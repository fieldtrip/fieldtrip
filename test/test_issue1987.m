function test_issue1987

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY read_brainvision_eeg

%% Unit test for update to read_brainvision_eeg.m
%
% requires test data "seg_noxfo1" in the same folder as this function
% dat:  https://osf.io/mvuzt/
% vhdr: https://osf.io/c6yq7/
% vmrk: https://osf.io/6qu4y/
%
% 1. read header             (ft_read_header)
% 2. read data in MICROVOLTS (ft_read_data)
% 3. compare 1st and 100th sample from all channels against stored values
%
% Stored values are in VOLTS from mne.io.read_raw_brainvision
% <https://mne.tools/stable/generated/mne.io.read_raw_brainvision.html>
%
% Expect: data * 10^-6 == stored values (volts to microvolts)
% with error <10^-13

% go to where the test data is located
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1987'))

hdr = ft_read_header('seg_noxfo1.vhdr');

% read raw data in MICROVOLTS
raw = ft_read_data('seg_noxfo1.dat');

% GOLD STANDARD values from mne.io.read_raw_brainvision in VOLTS
first_sample = [7.87292492e-05, 9.45078444e-05, 1.08382639e-04, 9.06084607e-05, ...
  9.41521391e-05, 1.27865892e-04, 1.84153825e-05, 7.93829051e-05, ...
  8.76488965e-05, 6.53165902e-05, 9.09365019e-05, 4.19619186e-05, ...
  9.03568586e-05, 8.79602216e-05, 6.22964448e-05, 6.64892360e-05, ...
  1.14737543e-04, 4.16034933e-05, 1.98621509e-04, 6.31107836e-05, ...
  1.21482065e-04, 6.33632289e-05, 6.20244227e-05, 5.14866151e-05, ...
  6.52697839e-05, 3.93936240e-05];
hundr_sample = [7.10163812e-05, 9.96401458e-05, 9.13076874e-05, 8.18402495e-05, ...
  8.09276441e-05, 1.17282465e-04, 1.69571098e-05, 8.36039747e-05, ...
  1.02514001e-04, 9.44410642e-05, 1.00473428e-04, 3.28999904e-05, ...
  8.70283902e-05, 8.66390318e-05, 6.32346431e-05, 6.71046917e-05, ...
  1.11633760e-04, 6.67019968e-05, 2.17220386e-04, 8.57964489e-05, ...
  1.26016558e-04, 8.72791991e-05, 3.02727266e-05, 3.69033323e-05, ...
  4.44250992e-05, 2.83058476e-05];

% test against MNE data, scaled to VOLTS
first_err = (raw(:,1) .* 10^-6)  - first_sample';
hundr_err = (raw(:,100) .* 10^-6) - hundr_sample';

max_err = max([first_err; hundr_err]);

if max_err > 10^-12
  error('FAIL: Data does not match expected values (error: %g)\n', max_err);
else
  fprintf('SUCCESS: Data are consistent with expected values (error: %g)\n', max_err);
end
