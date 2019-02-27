function [err,info] = sqdwrite(source,destination,varargin)
% [ERR,INFO] = SQDWRITE(SOURCE,DEST,'ACTION',<'APPEND','OVERWRITE'>,...
%           'CHANNELS',[C1 C2 C3,...],'SAMPLES',[S1 SN],'DATA',DATA_ARRAY); 
% Writes new Channel Data to .sqd file given by string DEST. The template 
% (Acquisition and Setup Information) is chosen from SOURCE sqd-file. It assumes
% that the channeldata is in a format (shown below) and appropriately modifies
% the information regarding the data written to the sqd-file. It returns ERR
% which is the number of values successfully written. It returns -1 if
% there has been a complete failure. 
% Motive behind this function is to able to write data (obtained from sqdread)
% which has been modified (filtered,etc.) in MATLAB back to a sqd-file.
% MEG data should be in units of picotesla (pT).
%
% Sometimes the data array DATA_ARRAY must be given for all 192 channels.
%
% Usage:
% err = sqdwrite(source,destination,data);
%   Write new data to sqd-file. Data should be arranged as following:
%	|Chan(0)Sample(1)   Chan(1)Sample(1)	... Chan(end)Sample(1)  |
%	|Chan(0)Sample(2)   Chan(1)Sample(2)	...	        .	        |
%	|   .			        .			                .	        |
%	|   .			        .			                .           |
%	|Chan(0)Sample(end) Chan(1)Sample(end)	... Chan(end)Sample(end)|
% PUTDATA will take the channel and sample assignments from the data
% entered.
% 
% err = sqdwrite(source,destination,...
%   'ACTION',<'Append','Overwrite'>,...
%   'CHANNELS',[CH1 CH2 .. CHn],...
%   'SAMPLES',[SAMP1 SAMPN],...
%	'DATA',data);	
%   Write data to an already existing sqd-file. You may specify the Action
%   to be 'Append' or 'Overwrite'.  
%   APPEND merely appends to an existing file. It assumes the data is in the above
%   format (ie. (Samplenum,Channelnum)).The input arguments - 'CHANNELS' and
%   'SAMPLES' will be ignored if given.
%   OVERWRITE requires destination file to already exist. It writes the raw data 
%   writes over the data at the location specified by Channelnumber and 
%   sample number. It returns an error if there is no data at the specified 
%   location. If neither Channelnumber nor Samplesnumber is specified, it
%   is assumed that overwrite will take place for entire existing data.
%   Similar logic applies if only ChannelNumber or only Samplenumber is
%   specified. Channels must be contiguous.
%
% Examples:
% If the path for the sqd-toolbox has been set properly and all of the files present 
% (see README), there should be a sample raw sqd-file "testsqd.sqd" for testing. The
% sqd-file "testsqd.sqd" has data for 192 channels, 10 seconds of data and sampled
% at 1000Hz. The data is sinusoids of varying frequencies for the first 157 channels
% and random triggers for the rest. For the examples below, "testsqd.sqd" is used as
% the original sqd-file. There is also a file "sqdtestdata.mat" which has a 2-D array
% "data" which is 1000x192. "data" will used as the test data to be written to the
% new sqd-file.
% 1. To make a new sqd file "new.sqd" using an array "data" which has data for 192 
%    channels and 1000 samples:
%   load sqdtestdata.mat;           % Load "data" into workspace. This is the data
%                                   % you wish to write to the new sqd-file
%   size(data)                      % Make sure data is provided for all channels
%                                   % Here size(data) = (1000 x 192) ie (NumSamples x NumChannels)
%   newsqdfile     = 'new.sqd';     % Name of the new sqdfile
%   templatesqdfile= 'testsqd.sqd'; % Template sqd-file. This is the original
%                                   % sqd-file from which all the header (acquisition
%                                   % information) will be copied into the newfile.
%                                   % sqdwrite will not work without a valid template-sqd-file
%   sqdwrite(templatesqdfile,newsqdfile,data); % Call sqdwrite
%   ls -l new.sqd                   % To confirm that the file was made
%   % To check if the creation has worked, read data from the new.sqd using sqdread and plot
%   checkdata = sqdread(newsqdfile,'Channels',0); % read in channel# 0
%   plot(checkdata);                % should be a sinusoid
%
% 2. To add more data to an already existing sqd-file:
%    Here we will be appending data (500 samples x 192 channels) to the file created
%    from example #1.
%   newdata = rand(500,192)*10^4;   % Get the data that you wish to append to the 
%                                   % existing sqd-file
%   size(newdata)                   % Make sure data is provided for all channels
%                                   % Here size(data) = (500 x 192) ie (NumSamples x NumChannels)
%   existingsqdfile = 'new.sqd';    % The existing sqd-file
%   templatesqdfile= existingsqdfile; % Since file already exists, it can act as its own template
%   sqdwrite(templatesqdfile,existingsqdfile,...
%       'Action','Append','Data',newdata); % Append to file
%   % To check if the append has worked, read data from the new.sqd using sqdread and plot
%   checkdata = sqdread(existingsqdfile,'Channels',0); % read in channel# 0
%   plot(checkdata);                % should be a sinusoid with noise at the end
%
% 3. To overwrite part of the data in an already existing file:
%    Here we are going to overwrite the data for channel N, from samples S1 to SN.
%   N  = 0;                         % Channel numner = 0;
%   S1 = 101;                       % Start sample
%   SN = 600;                       % End sample => NumSamples = (600-101+1) = 500
%   newdata = rand(500,1)*10^4;     % Get the data that you wish to overwrite to the 
%                                   % existing sqd-file
%   size(newdata)                   % Make sure data is provided for all specified channels
%                                   % Here size(data) = (500 x 1) ie (NumSamples x NumChannels)
%   existingsqdfile = 'new.sqd';    % The existing sqd-file
%   templatesqdfile= existingsqdfile; % Since file already exists, it can act as its own template
%   sqdwrite(templatesqdfile,existingsqdfile,...
%       'Action','Overwrite','Channels',N,...
%       'Samples',[S1 SN],'Data',newdata); % Overwrite
%   % To check if the append has worked, read data from the new.sqd using sqdread and plot
%   checkdata = sqdread(existingsqdfile,'Channels',N); % read in channel# 0
%   plot(checkdata);                % should be a sinusoid with noise at the end 
%                                   % and noise in between (101 to 600)
%
%   See also, 
%   sqdhandle/putdata,getdata,sqdhandle,get,set

% Note:
% This function is merely a wrapper around putdata method of SQDHANDLE.
%
% To be included in the coming versions:
% #Action - INSERT allows you to insert data into an already existing sqdfile. You 
%   must specify the location of insertion by  specifying 'Samples'
% #Variable numebr of channels instead of assuming all 192 channels.


if nargin<3,
    error('Atleast 3 input arguments required');
end;

if ischar(source)
    info = sqdhandle(source);           % Get sqdhandle_obj from template sqd-file
    if nargin==3,   % Ex. sqdwrite(sourcefile,destfile,data);
        err = putdata(info,destination,'Data',varargin{1});
    else            % Ex. sqdwrite(sourcefile,destfile,'Action','Append','Data',data)
        err = putdata(info,destination,varargin{:});
    end;
    if nargout>=2
        info = sqdhandle(destination,'info');
    end;
else
   error('First argument must be a valid filename'); 
end;
