function test_dccnpath

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY dccnpath
% DATA public

global ft_default
ft_default.dccnpath = [];

% remember the original directory
origdir = pwd;

hostname = gethostname();
if startsWith(hostname, 'DCCN') || startsWith(hostname, 'dccn')

  %% Alternative0: Finds the right path for users with access to DCCN storage, either on the compute cluster or on DCCN desktops

  filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf');
  if ispc
    if ~strcmp(filename, 'H:\common\matlab\fieldtrip\data\ftp\test\edf\testAlphaIR20170321-0.edf')
      error('dccnpath() does not find the right DCCN path')
    end
  else
    if ~strcmp(filename, '/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf')
      error('dccnpath() does not find the right DCCN path')
    end
  end

else

  %% Alternative1: Test with local files in the present working directory

  tmpdir = tempname;
  mkdir(fullfile(tmpdir, 'ftp', 'test', 'edf'));
  cd(fullfile(tmpdir, 'ftp', 'test', 'edf'));
  websave('testAlphaIR20170321-0.edf', 'https://download.fieldtriptoolbox.org/test/edf/testAlphaIR20170321-0.edf');

  filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf');
  if ~startsWith(filename, tmpdir)
    error('dccnpath() does not find the right working directory');
  end

  % return to the original directory
  cd(origdir);

  %% Alternative2: It downloads test data to the local computer, if test data is not already downloaded. It has 6 different tests that are listed below

  %% Test1: When ft_default.dccnpath is given by the user, dccnpath() should always select alternative2

  % ... continue with the temporary copy from the previous test
  ft_default.dccnpath = tmpdir;

  filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf');
  if ~startsWith(filename, tmpdir)
    error('dccnpath() does not find the right working directory');
  end

  if startsWith(ft_default.dccnpath, tempdir)
    rmdir(ft_default.dccnpath, 's');
  end

  %% Test2: Automatically downloading a file from the HTTPS download server

  tmpdir = tempname;
  ft_default.dccnpath = tmpdir;

  filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls');

  if ~exist(filename, 'file')
    error('File is not downloaded from the HTTPS FieldTrip server');
  end

  %% Test3: File exists in the local copy and it doesn't get downloaded again

  % ... continue with the temporary copy from the previous test

  filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls'); %'/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls' has already been downloaded by test2
  if ~exist(filename, 'file')
    error('File exists in the local copy, but dccnpath() can not find it');
  end

  if startsWith(ft_default.dccnpath, tempdir)
    rmdir(ft_default.dccnpath, 's');
  end

  %% Test4: Downloading a folder from the HTTPS download server

  tmpdir = tempname;
  ft_default.dccnpath = tmpdir;

  foldername = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds');
  if ~exist(foldername, 'dir') % Here I check only if the main dir exists and not if this dir has the necessary contents.
    error('Folder is not downloaded from the HTTPS FieldTrip server');
  end

  %% Test5: Folder already exists in the local copy and it doesn't get downloaded again

  % ... continue with the temporary copy from the previous test

  foldername = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds'); %'/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds/' has already been downloaded by test4
  if ~exist(foldername, 'dir')
    error('Folder exists in the local copy, but dccnpath() can not find it');
  end

  if startsWith(ft_default.dccnpath, tempdir)
    rmdir(ft_default.dccnpath, 's');
  end

  %% Test6: When ft_default.dccnpath is not specified, data should be saved automatically to temporary directory

  % Download public data
  if isfield(ft_default, 'dccnpath')
    ft_default = rmfield(ft_default, 'dccnpath');
  end

  hostname = gethostname();
  if ~startsWith(hostname, 'DCCN') && ~startsWith(hostname, 'dccn') % this test applies to external contributors and not for users in the DCCN cluster
    filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls');
    if ~exist(filename, 'file') || ~contains(filename, tempdir)
      error('Data were not automatically downloaded to tempdir');
    end
  end

  % Do not download private data
  if isfield(ft_default, 'dccnpath')
    ft_default = rmfield(ft_default, 'dccnpath');
  end

  hostname = gethostname();
  if ~startsWith(hostname, 'DCCN') && ~startsWith(hostname, 'dccn') % this test applies to external contributors and not for users in the DCCN cluster
    try
      filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/bug1027.mat');
      ok = true;
    catch
      ok = false;
    end
    assert(~ok, 'Private data should not have been downloaded')
  end

end % if inside or outside the DCCN

