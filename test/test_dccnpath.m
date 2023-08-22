function test_dccnpath

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY dccnpath
% DATA public

%% Alternative0: Finds the right path for users that have access to the DCCN computer cluster. It will replace '/home' by 'H:' and will replace forward by backward slashes.

try
filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf'); 
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
    hostname=gethostname();
    if startsWith(hostname, 'DCCN')
        error('dccnpath() does not find the right DCCN path')
    end
end


%% Alternative1: Allows to test with local files in the present working directory

ft_defaults % '/home/common/matlab/fieldtrip/data/ftp/test/fileio' now belongs to the present working directory
filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/fileio');
if ~contains(filename, 'fileio')
        error('dccnpath() does not find the right working directory');
end


%% Alternative2: It downloads test data to the local computer, if test data is not already downloaded. It has 6 different tests that are listed below

%% Test1: When ft_default.dccnpath is given by the user, dccnpath() should always select alternative2 (this is not done now)

% if ~isempty(ft_default.dccnpath)
%     filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf');
%     if ~contains(filename, ft_default.dccnpath)
%             error('The path defined by ft_default.dccnpath was not used, while it should have');
%     end
% end
% 
% delete(filename) 

    
%% Test2: Downloading a file from the HTTPS download server
testdata=fullfile(tempdir,'testdata');
mkdir(testdata);

global ft_default
ft_default.dccnpath=testdata;

filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls');
if ~exist(filename, 'file')
    error('File is not downloaded from the HTTPS FieldTrip server');
end


%% Test3: File exists in the local copy and it doesn't get downloaded
ft_default.dccnpath=fullfile(testdata,'ClassFile.cls');

filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls'); %'/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls' has already been downloaded by test2
if ~exist(filename, 'file')
    error('File exists in the local copy, but dccnpath() can not find it');
end
 
delete(ft_default.dccnpath);

%% Test4: Downloading a folder from the HTTPS download server 

ft_default.dccnpath=testdata;

foldername=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds/');
if ~exist(foldername, 'dir') % Here I check only if the main dir exists and not if this dir has the necessary contents. 
    error('Folder is not downloaded from the HTTPS FieldTrip server');
end


%% Test5: Folder exists in the local copy and it doesn't get downloaded 

foldername=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds/'); %'/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds/' has already been downloaded by test4
if ~exist(foldername, 'dir')
    error('Folder exists in the local copy, but dccnpath() can not find it');
end


rmdir(testdata, 's');

%% Test6: When ft_default.dccnpath is not specified then data should be saved automatically to tempdir

% Download
ft_default=rmfield(ft_default,'dccnpath');

delete(fullfile(tempdir,'data/ftp/test/ctf/Subject01.ds/ClassFile.cls'));

filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls'); 
if ~exist(filename, 'file') || ~contains(filename,tempdir)
    error('Data are not automatically downloaded to tempdir when ft_default.dccnpath is not specified');
end


% Do not download
if isfield(ft_default,'dccnpath')
    ft_default=rmfield(ft_default,'dccnpath');
end 

filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/ClassFile.cls'); 
if ~exist(filename, 'file') || ~contains(filename,tempdir)
    error('Data exist in the tempdir, but dccnpath() can not find it');
end

