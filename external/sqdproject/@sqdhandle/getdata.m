function data = getdata(t,varargin)
% GETDATA Get Channel data pointed by handle 't' using specifications from
% other options like 'Samples', 'Channels','Format'
% Usage:
% data = getdata(t);
%   Get all the raw data pointed by 't'. Data is arranged as following:
%	|Chan(0)Sample(1)   Chan(1)Sample(1)	... Chan(end)Sample(1)  |
%	|Chan(0)Sample(2)   Chan(1)Sample(2)	...	        .	        |
%	|   .			        .			                .	        |
%	|   .			        .			                .           |
%	|Chan(0)Sample(end) Chan(1)Sample(end)	... Chan(end)Sample(end)|
% data = getdata(t,'Channels',[0 1 2 3 4],'Samples',[2001 3000],...
%	'Format','double');	
%   Get data specified by Channels,Samples, and Format fields of inputs. Data
%   will be arranged as shown above in a increasing order

if nargin<1
    error('Atleast one input argument is required');
end;
if ~isa(t,'sqdhandle')
    error('First argument must belong to class - sqdhandle');
end;
if rem(nargin-1,2)
    error('Property and range must come in pairs');
end;

% Init Parameters required to read data
raw_offset      = get(t,'RawOffset');         % Location of start for reading data
framesamples    = get(t,'SamplesAvailable');  % Number of samples available er channel    
channelcount    = get(t,'ChannelCount');  
channels        = [0:channelcount-1];              
numsamples      = framesamples;               % Number of samples to be read
format	        = get(t,'Datatype');           
oldformat       = get(t,'Datatype');          % Format in which was written
newformat       = 'double';                   % Format to which data read will be converted
acquisitiontype = get(t,'AcquisitionType');   % need to know if epoched data (type 3)
if acquisitiontype==3                         % Need to get all epochs
    actaveragecount = get(get(t,'acqparam'),'ActAverageCount'); % number of epochs
    numsamples = numsamples * actaveragecount;
end
samples         = [1 numsamples];

% Check user input for channels,samples and format
if ~rem(nargin-1,2)&fix((nargin-1)/2) 
    fields = {'Channels','Samples','Format'};
    inargs = varargin;
    while length(inargs)>=2
        prop = inargs{1};
        range = inargs{2};
        inargs = inargs(3:end);
        switch lower(prop)
        case 'channels'
            channels = sort(range);
            if (find(channels>channelcount-1))|(find(channels<0))
                error('Sorry incorrect choice for channel numbers');
            end;
        case 'samples'
            samples = sort(range);
            if length(samples)>2
                error('Sorry - Please input only start and stop of the range')
            end;
            if (samples(1)<=0)|(samples(2)>numsamples)
                error('Samples chosen are out of range');
            end;
        case 'format'
            if ~strcmpi(range,'double')
                format = [format '=>' range];
                newformat = range;
            end;
        otherwise
            error('Sorry - Incorrect Property used');
        end;
    end;
end;
fid = fopen(t.FileName,'rb','l');
numsamples  = samples(2)-samples(1)+1;
switch lower(get(t,'Datatype'))
case 'int16'
    numbytes = 16/8;
case 'double'
    numbytes = 64/8;
end;
numchannels = length(channels);
fseek(fid,raw_offset,-1);

if (numsamples == framesamples) & (numchannels == channelcount) 
    % If reading all the data that is there in the file
    data = zeros(framesamples*channelcount,1);
    data = fread(fid,framesamples*channelcount,format);
    data = reshape(data,channelcount,framesamples)';
elseif all(diff(channels)==1)
    % if the channels specified are contiguous Ex. 0,1,2 and not 1,3,5
    data = zeros(numsamples,numchannels);
    offset = (channels(1)+channelcount*(samples(1)-1))*numbytes; % Goto 1st channel 1st Sample
    fseek(fid,offset,0);
    data = fread(fid,...                                        % FILE ID
        numchannels*numsamples,...                              % DATASAMPLES
        [num2str(numchannels) '*' format],...                   % PRECISION
        (channelcount - channels(end)+channels(1)-1)*numbytes); % SKIP
    data = reshape(data,numchannels,numsamples)';       % Reshape
else % If non-contiguous channels	
    data = zeros(numsamples,numchannels);
    for i = 1:length(channels)
        chan_offset = channels(i);
        sample_offset = channelcount*(samples(1)-1);
        fseek(fid,raw_offset+(sample_offset+chan_offset)*numbytes,-1);
        data(:,i) = fread(fid,...                               % FILE ID
            numsamples,...                                      % Num data Samples
            format,...                                          % Precision
            (channelcount-1)*numbytes);                         % SKIP
    end;
end;
fclose(fid);

if strcmpi(newformat,'double')
    gain = get(t,'SensitivityGain');
    gain = gain(channels+1);
%     ampgain = get(t,'OutputGain')/get(t,'InputGain');
%     switch lower(oldformat)
%     case 'int16'        
%         convfactor = ones(size(data,1),1)*ampgain*gain/(2^12);
%     case 'double'
%         convfactor = ones(size(data,1),1)*ampgain*gain/(2^12);
%     end;
    ampgain = get(t,'OutputGain')*get(t,'InputGain');
    % 10 is NI AD card ranve of volts
    % 2^12 is conversion of int16 to doubles
    % gain is conversion of plain double to picoTesla on an individual channel basis
    % ampgain undoes the hardware gain of created by the InputGain and the OutputGain.
    convfactor = ones(size(data,1),1)*10/(2^12)*gain/ampgain;
    data = convfactor.*data;
end;

if (acquisitiontype==3) & (size(data,1)==actaveragecount*framesamples)  % Fold epochs if all data was read
    % that is, turn Matrix of size (actaveragecount*framesamples,nchannels) into
    % Matrix of size (framesamples, nchannels, actaveragecount)
    data=permute(reshape(data.',[size(data,2),size(data,1)/actaveragecount,actaveragecount]),[2,1,3]);
end

