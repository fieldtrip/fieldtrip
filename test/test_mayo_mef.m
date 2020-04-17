% TEST_MAYO_MEF an example to read MEF 3.0/2.1 dataset into FieldTrip
% 
% Note:
%   `test_mayo_mef.m` should be in the folder /root/of/fieldtrip/test.

%% Test reading MEF 3.0 dataset
clear

% look for session path of data
p = mfilename('fullpath');
sesspath = fullfile(fileparts(fileparts(p)), 'external', 'mayo_mef', 'sample_mef',...
    'mef_3p0.mefd');

% set the password
password = struct('Level1Password', 'password1', 'Level2Password',...
    'password2', 'AccessLevel', 2);


%% read data with ft_read_data() but specifying time interval using seconds
% -------------------------------------------------------------------------
% Let's import 1.5 seconds data at the beginning of the recording
hdr = ft_read_header(sesspath, 'password', password);
fs  = hdr.samplingrate;
dat = ft_read_data(sesspath,...
    'begsample', 1,...
    'endsample', 10000,...
    'header', hdr,...
    'password', password,...
    'chanindx', [4 1 2 3]); % the order of read data can be decided by
                            % key-value 'chanindx'

figure
plot(t, dat')
xlim([0 1.5]+be_second(1))
xlabel('Time (s)')
legend(hdr.label{4}, hdr.label{1}, hdr.label{2}, hdr.label{3})


%% read data with ft_preprocessing()
% ---------------------------------
% Let's import 5 trials/epochs. Each trial is 1.50 second long.  The
% trigger time of the 5 trials are at 0.5, 2.0, 3.5, 5.0 and 6.5 seconds,
% with the pre-stimulus length of 0.5 second.
% 
% setup trial information
trig = [.5, 2, 3.5, 5, 6.5]; % in seconds
n_trig = length(trig); % number of triggers
trigger = mef_ft.SessionUnitConvert(trig, in_unit, out_unit)';
prestim = mef_ft.SessionUnitConvert(.5, in_unit, out_unit)*ones(n_trig, 1);
poststim = mef_ft.SessionUnitConvert(1., in_unit, out_unit)*ones(n_trig, 1);
trl = [trigger-prestim+1, trigger+poststim, -prestim];

% read the data
cfg = [];
cfg.dataset = sesspath;
cfg.password = password;
cfg.header = hdr;
cfg.trl = trl;
dat_ieeg = ft_preprocessing(cfg);

% plot the data
cfg.viewmode = 'vertical';
brwview = ft_databrowser(cfg, dat_ieeg);

%% Copyright (c) 2020 MNL group
% Created on Fri 04/17/2020 12:21:46.992 PM
% Revision: 0.1  Date: Fri 04/17/2020 12:21:46.992 PM
%
% Multimodal Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)