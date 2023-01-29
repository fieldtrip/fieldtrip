function test_mayo_mef()

    % WALLTIME 00:10:00
    % MEM 3gb
    % DEPENDENCY read_mayo_mef21 read_mayo_mef30

    % Copyright 2020-2023 Richard J. Cui. Created: Sat 03/21/2020 10:35:23.147 PM
    % $Revision: 0.4 $  $Date: Thu 01/26/2023 10:35:01.233 PM $
    %
    % Mayo Foundation for Medical Education and Research
    % Mayo Clinic St. Mary Campus
    % Rochester, MN 55905
    %
    % Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)

    % =========================================================================
    % Note of sample dataset
    % =========================================================================
    % Go to https://github.com/jiecui/mef_reader_fieldtrip
    % and follow the instructions to download the sample dataset. Set the path
    % to the sample dataset below.

    %% set path to sample dataset
    % please set the path to sample dataset
    mef21_data = '/path/to/sample_mef/mef_2p1';
    mef30_data = '/path/to/sample_mef/mef_3p0.mefd';

    %% install toolbox and build mex files
    disp('-----------------------------------')
    disp('install toolbox and build mex files')
    disp('-----------------------------------')
    ft_hastoolbox('mayo_mef', 1);
    setup_mayo_mex()

    %% MEF version 2.1
    disp(' ')
    disp('--------------------------')
    disp('testing MEF version 2.1...')
    disp('--------------------------')

    % the password of MEF 2.1 sample dataset
    mef21_pw = struct('Subject', 'erlichda', 'Session', 'sieve', 'Data', '');

    % low-level
    hdr = ft_read_header(mef21_data, 'password', mef21_pw);
    dat = ft_read_data(mef21_data, 'password', mef21_pw);
    event = ft_read_event(mef21_data, 'password', mef21_pw);

    % high-level
    cfg = [];
    cfg.dataset = mef21_data;
    cfg.password = mef21_pw;
    data = ft_preprocessing(cfg);

    %% MEF version 3.0
    disp(' ')
    disp('--------------------------')
    disp('testing MEF version 3.0...')
    disp('--------------------------')

    % the password of MEF 3.0 sample dataset
    mef30_pw = struct('Level1Password', 'password1', 'Level2Password', ...
    'password2', 'AccessLevel', 2);

    % low-level
    hdr = ft_read_header(mef30_data, 'password', mef30_pw);
    dat = ft_read_data(mef30_data, 'password', mef30_pw);
    event = ft_read_event(mef30_data, 'password', mef30_pw);

    % high-level
    cfg = [];
    cfg.dataset = mef30_data;
    cfg.password = mef30_pw;
    data = ft_preprocessing(cfg);

    disp('********************')
    disp('* You are all set! *')
    disp('********************')

end % function

% [EOF]
