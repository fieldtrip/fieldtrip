function test_dccnpath

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY
% DATA public private

%% Alternative0

global ft_default;

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


%% Alternative1

ft_defaults
filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/fileio');
if ~strcmp(filename, 'fileio')
        error('Alternative0 does not work');
end


%% Alternative2

% Create the right paths inside tempdir in order for alternative2 to work
global ft_default;
ft_default.dccnpath = tempdir;
ft_default.dccnpath = ft_default.dccnpath(1:end-1); % we do not want the temppdir to end with a '/' or '\'


if ispc
    path=strcat(ft_default.dccnpath,'\data\ftp');
else
    path=strcat(ft_default.dccnpath,'/data/ftp');
end
if ~isfolder(path) % make the directories data/ftp, if they do not exist in the temp path
    mkdir(path);
end


%% Test1: When ft_default.dccnpath is given by the user, dccnpath() should select alternative2

if ~isempty(ft_default.dccnpath)
    filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/edf/testAlphaIR20170321-0.edf');
    if ~contains(filename, ft_default.dccnpath)
            error('The alternative2 was not selected while it should have');
    end
end

delete(filename) 

    
%% Test2: Downloading a file from the HTTPS download server
 
filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/BadChannels');
if ~exist(filename, 'file')
    error('Alternative2 does not work for test2');
end


%% Test3: File exists in the local copy and it doesn't get downloaded
filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds/BadChannels');
if ~exist(filename, 'file')
    error('Alternative2 does not work for test3');
end

delete(filename)  


%% Test4: Downloading a folder from the HTTPS download server 

filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds/');
if ~exist(filename, 'dir') % Here I check only if the main dir exists and not if this dir has the contents that had to be downloaded 
    error('Alternative2 does not work for test4');
end


%% Test5: Folder exists in the local copy and it doesn't get downloaded 

filename=dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject02.ds/');
if ~exist(filename, 'dir')
    error('Alternative2 does not work for test5');
end

rmdir(filename,'s'); 


% Todo: FTP connection , file
% Todo: FTP connection, folder

% Todo: delete the folder data/ftp from the temp path. Then probably 'rmdir' and 'delete' are not
% needed.