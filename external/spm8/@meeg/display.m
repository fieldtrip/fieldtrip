function str = display(this)
% Method for displaying information about an meeg object
% FORMAT display(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: display.m 3219 2009-06-25 10:15:36Z vladimir $

str = ['SPM M/EEG data object\n'...
    'Type: ' type(this) '\n'...
    'Transform: ' transformtype(this) '\n'...
    num2str(nconditions(this)), ' conditions\n'...
    num2str(nchannels(this)), ' channels\n'
    ];

if strncmpi(transformtype(this),'TF',2)
    str = [str num2str(nfrequencies(this)), ' frequencies\n'];
end

str = [str ...
    num2str(nsamples(this)), ' samples/trial\n'...
    num2str(ntrials(this)), ' trials\n'...
    'Sampling frequency: ' num2str(fsample(this)) ' Hz\n'...
    'Loaded from file  %s\n\n'
    ];

if strncmpi(transformtype(this),'TF',2)
    str = [str  'Use the syntax D(channels, frequencies, samples, trials) to access the data\n'];
else
     str = [str  'Use the syntax D(channels, samples, trials) to access the data\n'];
end

str = [str  'Type "methods(''meeg'')" for the list of methods performing other operations with the object\n'...
    'Type "help meeg/method_name" to get help about methods\n'];

str = sprintf(str, fullfile(this.path, this.fname));

disp(str);