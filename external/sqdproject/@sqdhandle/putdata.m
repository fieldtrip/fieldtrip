function err = putdata(t,filename,varargin)
% ERR = PUTDATA(SQDHANDLE,FILENAME,'DATA',channeldata); 
% Writes new Channel Data to .sqd file given by 'filename'. The template 
% (Acquisition and Setup Information) is chosen from SQDHANDLE. It assumes
% that the channeldata is in a format (shown below) and appropriately modifies
% the information regarding the data written to the sqd-file. It returns ERR
% which is the number of values successfully written. It returns -1 if
% there hs been a complete failure. INFO is the sqdhandle to the new
% sqd-file.
% Assumption: Data is given for all 192 channels.
%
% Usage:
% err = putdata(t,filename,'Data',data);
%   Write new data to sqd-file. Data should be arranged as following:
%	|Chan(0)Sample(1)   Chan(1)Sample(1)	... Chan(end)Sample(1)  |
%	|Chan(0)Sample(2)   Chan(1)Sample(2)	...	        .	        |
%	|   .			        .			                .	        |
%	|   .			        .			                .           |
%	|Chan(0)Sample(end) Chan(1)Sample(end)	... Chan(end)Sample(end)|
% PUTDATA will take the channel and sample assignments from the data
% entered.
% 
% err = putdata(t,filename,...
%   'ACTION',<'Append','Overwrite'>,...
%   'CHANNELS',[CH1 CH2 .. CHn],...
%   'SAMPLES',[SAMP1 SAMPN],...
%	'DATA',data);	
%   Write data to an already existing sqd-file. You may specify the Action
%   to be 'Append' or 'Overwrite'.  
%   APPEND merely appends to an existing file. It assumes the data is in the above
%   format (Samplenum,Channelnum).The input arguments - 'CHANNELS' and
%   'SAMPLES' will be ignored.
%   OVERWRITE requires destination file to already exist. It writes the raw data 
%   writes over the data at the location specified by Channelnumber and 
%   sample number. It returns an error if there is no data at the specified 
%   location. If neither Channelnumber nor Samplesnumber is specified, it
%   is assumed that overwrite will take place for entire existing data.
%   Similar logic applies if only ChannelNumber or only Samplenumber is
%   specified. 
%
%   See also, 
%   sqdhandle/getdata,sqdhandle,get,set

% To be included in the coming versions:
% #Action - INSERT allows you to insert data into an already existing sqdfile. You 
%   must specify the location of insertion by  specifying 'Samples'
% #Variable numebr of channels instead of assuming all 192 channels.


if nargin<4
    error('Atleast 4 arguments required');
end;
if ~isa(t,'sqdhandle')
    error('First argument must belong to class - sqdhandle');
end;
if rem(nargin-2,2)
    error('Property and range must come in pairs');
end;
err = 0;
rawoffset = get(t,'RawOffset');


% Initialize default values for actions:
action = 'append';
format	 = get(t,'Datatype');
switch lower(format)
case 'int16'
    numbytes = 16/8;
case 'double'
    numbytes = 64/8;
end;

if ~rem(nargin-2,2)&fix((nargin-2)/2) 
    inargs = varargin;
    while length(inargs)>=2
        prop = inargs{1};
        range = inargs{2};
        inargs = inargs(3:end);
        switch lower(prop)
        case 'channels'
            channels = sort(range);
        case 'samples'
            samples = sort(range);
            if length(samples)>2
                error('Sorry - Please input only start and stop of the range')
            end;
            if (samples<=0)
                error('Value of Sample-Range must be >0');
            end;
        case 'data'
            data = range;
        case 'action'
            action = lower(range);
        otherwise
            error('Sorry - Incorrect Property used');
        end;
    end;
end;

if ~exist(filename,'file')
    fid = fopen(get(t,'FileName'),'rb','l'); % Open sqd-read-file to obtain template    
    if (fid==-1)
        disp('Error opening files');
        err = -1;
        return;
    end;
    template = fread(fid,rawoffset,'uchar'); % Read in template as bytes
    fclose(fid);
    fid = fopen(filename,'w','l');
    err = fwrite(fid,template,'uchar');  % Write template to file
    fclose(fid);
    fid = fopen(filename,'r+','l');
    updateinfo(t,fid,0);
    fclose(fid);
end;
% Switch on action:
switch action
case 'append'
    rawdata = convert2raw(t,data,format);
    fid = fopen(filename,'a+','l');
    [samplenum,channum] = size(data);      % Get number of chanels to be added
    err = err+fwrite(fid,rawdata,format);  % Write data and update output err
    fclose(fid);
    fid = fopen(filename,'r+','l');
    updateinfo(t,fid,samplenum);           % update info in file
    fclose(fid);
case 'insert'
    disp('sorry cant do it yet');       
    return;
case 'overwrite'
    fid = fopen(filename,'r','l');
%     fid = fopen(filename,'r');
    % Do checks so that size of data in the file, size of input data and
    % value of input arguments are in agreement
    % 1. Check for data present in file
    fseek(fid,0,1);                             % Go to end of the file
    if ftell(fid)<rawoffset                     % If new sqd-file
        error('No data exists in sqd-file - cannot overwrite');
    end;
    fclose(fid);
    [numsamples,numchannels] = size(data);
    numsamplesold = get(t,'SamplesAvailable');
    channelcount = get(t,'ChannelCount');
    if ~exist('samples','var')      
        samples = [1,numsamplesold];    % Default for number of samples
    end;
    if ~exist('channels','var')
        channels = 0:channelcount-1;    % Default for number of channels
    end;
    
    numsamples1 = samples(2)-samples(1)+1;
    numchannels1 = length(channels);
    % 2. Check for input arguments agreeing with each other
    if (numsamples~=numsamples1)|(numchannels~=numchannels1),
        disp('Error: Number of channels and number of samples given in rawdata');
        disp('       doesnot match those given as input arguments.');
        disp('       Check if input data has been given in the right format');
        error('quiting');
    end;
    % 3. Check if overwrite is legal ie there is enough data in file to
    % overwrite.
    if (numsamples>numsamplesold)|(numchannels>channelcount),
        error('Not enough data in file to overwrite');
    end;
    % End of checks
    % Convert to raw
    rawdata = convert2raw(t,data,format,channels);
    % Open file in overwrite mode
    [fid,message] = fopen(filename,'r+','l');    
%     [fid,message] = fopen(filename,'r+');
    % Begin writing data
    if (numsamples == numsamplesold) & (numchannels == channelcount)
        fseek(fid,rawoffset,-1);
        err = err+fwrite(fid,rawdata,format);
        % changing conditions here 11/20/08 dbhertz jzsimon
        % If all channels
    elseif (length(channels)==channelcount)
        offset = channelcount*(samples(1)-1)*numbytes;
        fseek(fid, rawoffset+offset,-1);
        err = err+fwrite(fid,...                                    % File ID
            rawdata,...                                             % DATA
            [num2str(numchannels) '*' format]);                     % PRECISION

        
    elseif (channels(2:end)-channels(1:end-1))==ones(1,numchannels-1)
        % If contiguous channels
        offset = (channels(1)+channelcount*(samples(1)-1))*numbytes-...
            (channelcount - channels(end))*numbytes;
        fseek(fid,rawoffset+offset,-1);
        err = err+fwrite(fid,...                                    % File ID
            rawdata,...                                             % DATA
            [num2str(numchannels) '*' format],...                   % PRECISION
            (channelcount - channels(end)+channels(1)-1)*numbytes); % SKIP
    else % If non-contiguous channels
        disp('Error: Channels must be contiguous');
        disp('       Please try again.');
        error('quiting');
%         for i = 1:length(channels)
%             chan_offset = channels(i);
%             sample_offset = channelcount*(samples(1)-1);
%             fseek(fid,rawoffset+(sample_offset+chan_offset)*numbytes-(channelcount-1)*numbytes,-1);
%             err = err + fwrite(fid,...              % FILE ID
%                 rawdata(:,i),...                    % DATA
%                 format,...                          % PRECISION
%                 (channelcount-1)*numbytes);         % SKIP
%         end;
    end;
    fclose(fid);
end;
% fclose('all');
return;


function rawdata = convert2raw(t,data,format,channels)
% % Converts data from 'double' to format
% switch format
%     case 'double'
%         gain = get(t,'SensitivityGains');
%         if nargin<4,
%             channels = 1:get(t,'ChannelCount');
%         end;
%         gain = gain(channels)';
%         inversegain = (ones(size(gain))./gain);
%         convfactor = ones(size(data,1),1)*inversegain*(2^12)/50;
%         rawdata = data';
%         rawdata = rawdata';
%         rawdata = rawdata(:);
%     case {'int16','short'}
persistent FileName;
persistent Gain;
if nargin<4,
    channels = [1:get(t,'ChannelCount')]-1;
end;
%
%gain = get(t.Channels,channels+1,'Sensitivity.Gain');
%gain = gain';
%
if isempty(FileName)
    FileName = 'new';    
end;

if isempty(Gain)|(~strcmp(FileName,t.FileName))
    Gain = get(t,'SensitivityGain');
end;

chgain = Gain(channels+1);
% ampgain = get(t,'OutputGain')/get(t,'InputGain');
% inversegain = (ones(size(chgain))./chgain);
% convfactor = ones(size(data,1),1)*inversegain*(2^12)/ampgain;
ampgain = get(t,'OutputGain')*get(t,'InputGain');
inversegain = (ones(size(chgain))./chgain);
% see comments in getdata.m for what these different factors mean:
convfactor = ones(size(data,1),1)*ampgain*inversegain*(2^12)/10;
%    convfactor = ones(size(data,1),1)*10/(2^12)*gain/ampgain;
rawdata = convfactor.*data;
rawdata = rawdata';
rawdata = rawdata(:);

function [err,t] = updateinfo(t,fid,samplenum)
% Updates information regarding sample number

% fseek to appropriate location
SEEK_SET = -1;
fseek( fid, 128, SEEK_SET );
acq_offset = fread(fid,1,'long');
fseek( fid, acq_offset+(32+64)/8, SEEK_SET );
switch get(t,'AcquisitionType')
    case 1  
        % Read acquisition parameters
        fseek(fid,(32/8),0);
        if ~samplenum
            err = fwrite(fid,round(samplenum),'long');
        else
            samplenum = samplenum + fread(fid,1,'long');
            fseek(fid,-(32/8),0);
            err = fwrite(fid,round(samplenum),'long');
        end;
        t = set(t,'ActSamplesAcquired',round(samplenum));
        t = set(t,'SamplesAcquired',round(samplenum));
    case {2,3}
        % Read acquisition parameters
        if ~samplenum
            err = fwrite(fid,round(samplenum),'long');
        else
            samplenum = samplenum + fread(fid,1,'long');
            fseek(fid,-(32/8),0);
            err = fwrite(fid,round(samplenum),'long');
        end;
        t = set(t,'FrameLength',round(samplenum));
end;
return;
