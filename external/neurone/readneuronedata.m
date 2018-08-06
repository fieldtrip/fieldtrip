% readneuronedata() -  Read NeurOne Binary files
%
% Usage:  >> data = readneuronedata(dataFiles,nChannels,chans)
%
% ===================================================================
% Inputs: 
%         dataFiles      - A cell structure containing full paths for data 
%                          files to be read. This is input argument is
%                          required.
%         nChannels      - Another required input argument:
%                          the total number of channels.
%         chans          - A vector containing the channel numbers to be
%                          read. They have to be arranged in an ascending
%                          order. This input argument is optional.
% ====================================================================
% Output:
%         data           - A matrix where all numerical data is stored.
%                          The data is arranged so that the data for each
%                          channel is stored in a row according
%                          to its number. 
%
% ========================================================================
% NOTE:
% This file is part of the NeurOne data import plugin for EEGLAB.
% ========================================================================
% 
% Current version: 1.0.3.4 (2016-06-17)
% Author: Mega Electronics

function data = readneuronedata(dataFiles,nChannels,chans)
%% Parse input arguments
p = inputParser;
p.addRequired('dataFiles', @iscellstr);
p.addRequired('nChannels', @isnumeric);
p.addOptional('chans', 0, @isnumeric);

p.parse(dataFiles, nChannels, chans);
arglist=p.Results;

%% Prepare reading data
nDataFiles = numel(arglist.dataFiles);
chans = arglist.chans;
nChannels = arglist.nChannels;

if chans==0
    chans=1:nChannels;
end

% Get all file sizes
fSize=zeros(1,nDataFiles);
tmp={};

for k=1:nDataFiles
    tmp{k,1}=dir(dataFiles{k,1});
    fSize(1,k)=tmp{k,1}.bytes;
end

fSizes=cumsum(fSize);

% Total size of the data
totalSize=sum(fSize);
% Total number of data points (per channel)
dataPntsTotal=(totalSize/4)/nChannels;
% Number of data points in each binary file (per channel)
dataPnts=(fSize./4)./nChannels;

%% Read binary data

fprintf('Allocating memory...\n')
data=zeros(numel(chans),dataPntsTotal);


% Option 1: Read all channels at once (faster method)
if numel(chans)==nChannels
    fprintf('Reading data...');
    data=data(:);
    readidx=[0 (fSizes/4)];
    for n=1:nDataFiles
        fid = fopen([dataFiles{n}], 'rb');
        data((1+readidx(n)):readidx(n+1),1)=fread(fid, ...
            fSize(n)/4,'int32');
        fclose(fid);
    end
    data=reshape(data,nChannels,dataPntsTotal);
    fprintf('%s','Done');
else
    % Option 2: Read only specific channels to save memory
    fprintf('Reading data... 0%%');
    readidx=[0 (fSizes/4)./nChannels];
    
    % Read channels one by one
    for k=1:numel(chans)
        for n=1:nDataFiles
        
            fid = fopen([dataFiles{n}], 'rb');
            fseek(fid, 4*(chans(k) - 1), 'bof');
            data(k,(1+readidx(n)):readidx(n+1)) = fread(fid, ...
                dataPnts(n),'int32',4*(nChannels-1))';
            fclose(fid);
        end
        
    % Print progress to the command window
    ii = floor(k/numel(chans)*100);
        if ii<10
            ii = [num2str(ii) '%'];
            fprintf(1,'\b\b%s',ii);
        elseif ii==100;
            ii = 'Done';
            fprintf(1,'\b\b\b%s',ii);
        else
            ii = [num2str(ii) '%'];
            fprintf(1,'\b\b\b%s',ii);
                
        end
    end                
end

fprintf('\n');

%% Convert data from nanovolts to microvolts
data=data./1000;

end
