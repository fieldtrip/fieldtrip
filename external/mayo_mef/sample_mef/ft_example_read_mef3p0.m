% FT_EXAMPLE_READ_MEF3P0 an example to read MEF 3.0 data set into FieldTrip

% clear all
clearvars

% look for session path of data
p = mfilename('fullpath');
sesspath = fullfile(fileparts(p), 'mef_3p0.mefd');

% set the password
password = struct('Level1Password', 'password1', 'Level2Password',...
    'password2', 'AccessLevel', 2);

% get the object to deal with FieldTrip
mef_ft = MEFFieldTrip_3p0(sesspath, password);

% data may be imported with 'mef_ft', but let's demonstrate the interaction
% with FieldTrip. Let's import 10 seconds data at the beginning of the
% recording with an assumption that the trigger was at 0 second.
in_unit = 'second';
be_second = [0, 10, 0]; % 10-second time of data from the start
out_unit = 'index';
be_sample = mef_ft.SessionUnitConvert(be_second, in_unit, out_unit);

% read data with ft_read_data()
dat = ft_read_data(sesspath, 'begsample', be_sample(1), 'endsample', be_sample(2),...
    'password', password);

% read data with ft_preprocessing
cfg = [];
cfg.dataset = sesspath;
cfg.password = password;
cfg.trl = be_sample;
dat_ieeg = ft_preprocessing(cfg);

% plot the data
cfg.viewmode = 'vertical';
brwview = ft_databrowser(cfg, dat_ieeg);

%% Copyright 2020 Richard J. Cui
% Created: Sun 03/22/2020  9:03:27.318 PM;
% Revision: 0.1  Date: Sun 03/22/2020  9:03:27.336 PM
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca