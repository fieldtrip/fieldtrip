function test_dccnpath

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY dccnpath
% DATA public

global ft_default
ft_default.dccnpath = [];

% remember the original directory
origdir = pwd;

%% Alternative0: Finds the right path for users that have access to the DCCN computer cluster. It will replace '/home' by 'H:' and will replace forward by backward slashes.

try
  filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf');
  if ispc
    if ~strcmp(filename, 'H:\common\matlab\fieldtrip\data\ftp\test\edf\testAlphaIR20170321-0.edf')
      error('Alternative0 does not work');
    end
  else
    if ~strcmp(filename, '/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf')
      error('Alternative0 does not work');
    end
  end
catch
  hostname = gethostname();
  if startsWith(hostname, 'DCCN') || startsWith(hostname, 'dccn')
    error('dccnpath() does not find the right DCCN path')
  end
end

%% Alternative1: Test with local files in the present working directory
hostname = gethostname();
if ~startsWith(hostname, 'DCCN') && ~startsWith(hostname, 'dccn') % this test applies to external contributors and not for users in the DCCN cluster
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
end

%% Alternative2: It downloads test data to the local computer, if test data is not already downloaded. It has 6 different tests that are listed below

%% Test1: When ft_default.dccnpath is given by the user, dccnpath() should always select alternative2 (this is not done now)

hostname = gethostname();
if ~startsWith(hostname, 'DCCN') && ~startsWith(hostname, 'dccn') % this test applies to external contributors and not for users in the DCCN cluster
    % ... continue with the temporary copy from the previous test
    ft_default.dccnpath = tmpdir;
    
    filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf');
    if ~startsWith(filename, tmpdir)
      error('dccnpath() does not find the right working directory');
    end

    rmdir(ft_default.dccnpath, 's')
end

%% Test2: Downloading a file from the HTTPS download server

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

if exist(ft_default.dccnpath,'dir')
    rmdir(ft_default.dccnpath, 's');
end

%% Test4: Downloading a folder from the HTTPS download server

tmpdir = tempname;
ft_default.dccnpath = tmpdir;

foldername = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds');
if ~exist(foldername, 'dir') % Here I check only if the main dir exists and not if this dir has the necessary contents.
  error('Folder is not downloaded from the HTTPS FieldTrip server');
end

%% Test5: Folder exists in the local copy and it doesn't get downloaded again

% ... continue with the temporary copy from the previous test

foldername = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds'); %'/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds/' has already been downloaded by test4
if ~exist(foldername, 'dir')
  error('Folder exists in the local copy, but dccnpath() can not find it');
end

if exist(ft_default.dccnpath,'dir')
    rmdir(ft_default.dccnpath, 's');
end

%% Test6: When ft_default.dccnpath is not specified, then data should be saved automatically to tempdir

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
