function test_mayo_mef(datapath)

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY read_mayo_mef21 read_mayo_mef30

% Copyright 2020-2023 Richard J. Cui. Created: Sat 03/21/2020 10:35:23.147 PM
% $Revision: 0.6 $  $Date: Wed 10/11/2023 12:44:00.481 AM $
%
% Mayo Foundation for Medical Education and Research
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)

% ======================================================================
% Note of sample dataset, when outside the DCCN
% ======================================================================
% Follow the instructions
% (https://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path/)
% to add FieldTrip to your MATLAB path. Then, download the sample
% dataset from https://github.com/jiecui/mef_reader_fieldtrip and follow
% the instructions to download the sample dataset. The input argument to
% this function should contain the path to where the data are located. Also
% the mef_reader_fieldtrip toolbox should be available on the path, and
% properly set up

if nargin<1
  datapath = dccnpath('/home/common/matlab/fieldtrip/data/test/original/lfp/mef');
end
mef21_data = fullfile(datapath, 'mef_2p1');
mef30_data = fullfile(datapath, 'mef_3p0.mefd');
med10_data = fullfile(datapath, 'med_1p0');

ft_hastoolbox('mayo_mef');

try
  % MEF version 2.1
  % ---------------
  disp(' ')
  disp('--------------------------')
  disp('testing MEF version 2.1...')
  disp('--------------------------')

  % the password of MEF 2.1 sample dataset
  mef21_pw = struct('Subject', 'erlichda', 'Session', 'sieve', 'Data', '');

  % low-level testing
  fprintf('low-level testing...\n')
  hdr = ft_read_header(mef21_data, 'password', mef21_pw);

  if isempty(hdr)
    warning('failed to read header of MEF 2.1 sample dataset')
  end

  dat = ft_read_data(mef21_data, 'password', mef21_pw);

  if isempty(dat)
    warning('failed to read data of MEF 2.1 sample dataset')
  end

  evt = ft_read_event(mef21_data, 'password', mef21_pw);

  if isempty(evt)
    warning('failed to read event of MEF 2.1 sample dataset')
  end

  % high-level testing
  fprintf('high-level testing...\n')
  cfg = [];
  cfg.dataset = mef21_data;
  cfg.password = mef21_pw;
  data = ft_preprocessing(cfg);

  if isempty(data)
    warning('failed to read data of MEF 2.1 sample dataset in high-level testing')
  end
catch
  ft_warning('failed to read data of MEF 2.1 sample data');
end

try
  % MEF version 3.0
  % ---------------
  disp(' ')
  disp('--------------------------')
  disp('testing MEF version 3.0...')
  disp('--------------------------')

  % the password of MEF 3.0 sample dataset
  mef30_pw = struct('Level1Password', 'password1', 'Level2Password', ...
    'password2', 'AccessLevel', 2);

  % low-level testing
  fprintf('low-level testing...\n')
  hdr = ft_read_header(mef30_data, 'password', mef30_pw);

  if isempty(hdr)
    warning('failed to read header of MEF 3.0 sample dataset')
  end

  dat = ft_read_data(mef30_data, 'password', mef30_pw);

  if isempty(dat)
    warning('failed to read data of MEF 3.0 sample dataset')
  end

  evt = ft_read_event(mef30_data, 'password', mef30_pw);

  if isempty(evt)
    warning('failed to read event of MEF 3.0 sample dataset')
  end

  % high-level testing
  fprintf('high-level testing...\n')
  cfg = [];
  cfg.dataset = mef30_data;
  cfg.password = mef30_pw;
  data = ft_preprocessing(cfg);

  if isempty(data)
    warning('failed to read data of MEF 3.0 sample dataset in high-level testing')
  end
catch
  ft_warning('failed to read data of MEF 3.0 sample data');
end

try
  % JM comment, this is likely to fail, since the compiled binaries
  % downloaded from darkhorseneuro.com (which are necessary for this
  % functionality) may have been compiled against libraries that are not
  % available on the test computer -> at least I did not manage to get the
  % below code to run on the DCCN cluster or my local Macbook pro.

  % MED version 1.0
  % ---------------
  disp(' ')
  disp('--------------------------')
  disp('testing MED version 1.0...')
  disp('--------------------------')

  % the password of MED 1.0 sample dataset
  med10_pw = struct('Level1Password', 'L1_password', 'Level2Password', ...
    'L2_password', 'AccessLevel', 2);


  % low-level testing
  fprintf('low-level testing...\n')
  hdr = ft_read_header(med10_data, 'password', med10_pw);

  if isempty(hdr)
    warning('failed to read header of MED 1.0 sample dataset')
  end

  dat = ft_read_data(med10_data, 'password', med10_pw);

  if isempty(dat)
    warning('failed to read data of MED 1.0 sample dataset')
  end

  evt = ft_read_event(med10_data, 'password', med10_pw);

  if isempty(evt)
    warning('failed to read event of MED 1.0 sample dataset')
  end

  % high-level testing
  fprintf('high-level testing...\n')
  cfg = [];
  cfg.dataset = med10_data;
  cfg.password = med10_pw;
  data = ft_preprocessing(cfg);

  if isempty(data)
    warning('failed to read data of MED 1.0 sample dataset in high-level testing')
  end
catch
  ft_warning('failed to read data of MED 1.0 sample data');
end
