function out = checkvalid(t)
% OUT = CHECKVALID(SQDHANDLE_OBJ); % Checks to see if the data in the sqd file
% is valid and returns corresponding error message.
% If the results of any test are fine,Then OUT(i) = 1
% If the results of any test could not be ddetermined then OUT(i) = 
% Tests:
%   1. If sqd-file doesnot exist - OUT(1) = -1
%   2. If sqd-file has more data than specified in the header - OUT(2) = -1
%      If sqd-file has more data than specified in the header - OUT(2) = -2

% Init output
out = [0,0];

% 1. Check to see if file exists
fid = fopen(get(t,'FileName'),'r','l');
if fid<0
    out(1) = -1;
    disp('Error: Could not find sqd file');
else
    out(1) = 1;
end;

numbytes = 0;   % Init number of bytes available for reading
if out(1)>0
    % 2. Check number of samples available in the sqdfile and 
    %    compare to number of samples specified in the header    
    
    % Get number of samples avaialble from the file directly
    fseek(fid,0,1);                 % Go to end of the file
    numbytes = ftell(fid)-get(t,'RawOffset');
    
    switch lower(get(t,'Datatype')) % Number of bytes for each sample
    case 'int16'
        szsamp = 2;
    case 'double'
        szsamp = 8;
    end;
    
    sampavailable = numbytes/szsamp;
    
    % Get number of samples from the header
    sampheader    = get(t,'SamplesAvailable')*get(t,'ChannelCount');
    
    if sampavailable>sampheader
        out(2) = -1;
        disp('Error: The SQD file has more data than can be accounted');
        disp('      for from the header');
    elseif sampavailable<sampheader
        out(2) = -1;
        disp('Error: The SQD file has insufficient amount of data');
    else
        out(2) = 1;
    end;
end;

% Print out a general message if there is an error
if min(out)<=0
    disp('If you think the above message is caused by a bug in the program');
    disp('Please email to shantanu@isr.umd.edu, with the following data');
    disp(['FileName     :' get(t,'FileName')]);
    disp(['Channels     :' num2str(get(t,'ChannelCount'))]);
    disp(['Datatype     :' get(t,'Datatype')]);
    disp(['AcqType      :' num2str(get(t,'AcquisitionType'))]);
    disp(['NumBytes     :' num2str(numbytes)]);
    disp(['SampHeader   :' num2str(get(t,'SamplesAvailable'))]);
    disp(['ActSampHd    :' num2str(get(t,'ActSamplesAcquired'))]);
    disp(['FrmlengthHd  :' num2str(get(t,'FrameLength'))]);
    disp(['Error        :' num2str(out)]);
end;