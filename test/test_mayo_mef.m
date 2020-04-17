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
% Let's import 1.5 seconds data at the beginning of the recording, where
% the trigger offset is -0.5 second
hdr = ft_read_header(sesspath, 'password', password);
fs  = hdr.Fs;
start_sample = 1;
end_sample = round(1.5*fs);
dat = ft_read_data(sesspath,...
    'begsample', start_sample,...
    'endsample', end_sample,...
    'header', hdr,...
    'password', password,...
    'chanindx', [4 1 2 3]); % the order of read data can be decided by
                            % key-value 'chanindx'

t = linspace(0, 1.5, end_sample-start_sample+1)-0.5;                            
figure
plot(t, dat')
line(gca, [0 0], ylim, 'Color', 'k')
xlabel('Time (s)')
ylabel('Signal amplitude')
legend(hdr.label{4}, hdr.label{1}, hdr.label{2}, hdr.label{3})


%% Copyright (c) 2020 MNL group
% Created on Fri 04/17/2020 12:21:46.992 PM
% Revision: 0.1  Date: Fri 04/17/2020 12:21:46.992 PM
%
% Multimodal Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)