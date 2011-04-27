function varargout = sqdread(filename,varargin)
% [DATA,INFO] = SQDREAD(FILENAME,'CHANNELS',[C1 C2 C3,...],'SAMPLES',[S1 SN],'FORMAT','double');    % READ SQD-FILE
%
% SQDREAD Get Channel data and info saved in sqd file 'filename' using 
% specifications from other options like 'Info','Samples', 'Channels',
% 'Format'. When the format is 'double', the MEG data has units of picotesla (pT).
% Usage:
% info = sqdread(filename,'Info');
%   where 'info' belongs to class sqdhandle. Use GET and SET methods to review or
%   change object properties. It contains experiment parameter values and
%   data acquisition parameters.
%
% [data,info] = sqdread(filename);
%   Get all the raw data saved in 'filename'. Data is arranged as following:
%	|Chan(0)Sample(1)   Chan(1)Sample(1)	... Chan(end-1)Sample(1)  |
%	|Chan(0)Sample(2)   Chan(1)Sample(2)	...	        .	        |
%	|   .			        .			                .	        |
%	|   .			        .			                .           |
%	|Chan(0)Sample(end) Chan(1)Sample(end)	... Chan(end-1)Sample(end)|
%
% [data,info] = sqdread(filename,'Channels',[0 1 2 3 4],'Samples',[2001 3000],...
%	'Format','int16');	
%   Get data specified by Channels,Samples, and Format fields of inputs, ie
%   in the above example, data is returned for Channel numbers [0,1,2,3,4] and Samples
%   from 2001 to 3000. The data is returned as 'INT16'. Data will be arranged
%   as shown above in a increasing order. You may choose above options in any
%   order. If option is not specified, default values are assumed. 
%   Default values for :
%       'Channels' - all available channels
%       'Samples'  - all available samples
%       'Format'   - 'double'
% 
% Examples:
% If the path for the sqd-toolbox has been set properly and all of the files present 
% (see README), there should be a sample raw sqd-file "testsqd.sqd" for testing. The
% sqd-file "testsqd.sqd" has data for 192 channels, 10 seconds of data and sampled
% at 1000Hz. The data is sinusoids of varying frequencies for the first 157 channels
% and random triggers for the rest. For the examples below, "testsqd.sqd" is used as
% the test sqd-file.
%
% 1. To read only the information regarding the data acquisition parameters:
%   a. info = sqdread('testsqd.sqd','info');
%       This will return you "info" which contains all the acquisition information.
%       Please make sure that when you specify the filename, it either includes the
%       filepath or you are in the directory where the file exists.
%   b. get(info);
%       This will tell you the parameters included in "info"
%   c. To get value of a parameter "fieldname" of "info".Choose fieldname
%      from parameters returned from get(info).
%      Ex. 
%       info.PatientInfo        % this returns the PatientInfo = Matt Cheely,...
%       info.InputGain          % this returns the Amplifier InputGain = 2
%       info.SamplesAvailable   % this returns number of samples available per channel = 10000
%      In general, to get value for parameter "fieldname"
%       info.fieldname
%       get(info,'fieldname')
%
% 2. To read all the data:
%   data = sqdread('testsqd.sqd');
%       This will return you a 2-Dimensional matrix array "data" whose size will
%   be 10000 x 192 ie (Number of Samples x Number of Channels)
%
% 3. To read data and info together:
%   [data,info] = sqdread('testsqd.sqd');
%
% 4. To read in data for only specific channels and/or samples
%   a. To read in all the data for only channel N
%      N = 0;                                       % First channel
%      data = sqdread('testsqd.sqd','Channels',N);  % Read in the data
%      plot(data);                                  % To look at the signal
%      size(data)                                   % should be 10000 x 1 
%                                                   % ie  (NumSamples x NumChannels)
%   b. To read in the data for all channels but only from Sample1 to SampleN
%      Sample1 = 1;                                 % First Sample 
%      SampleN = 1000;                              % 1000th sample = 1second of data
%      data = sqdread('testsqd.sqd','Samples',[Sample1 SampleN]); % Read in the data
%      size(data)                                   % should be 1000 x 192 
%                                                   % ie (NumSamples x NumChannels)
%   c. To read in data for channels N1,N2 and samples from Sample1 to SampleN
%      N1 = 0;                                      % First Channel
%      N2 = 99;                                     % 100th channel
%      Sample1 = 1;                                 % First Sample
%      SampleN = 1000;                              % 1000th sample = 1second of data
%      % Read in the data
%      data = sqdread('testsqd.sqd','Channels',[N1,N2],'Samples',[Sample1 SampleN]);
%      size(data)                                   % should be 1000 x 2 
%                                                   % ie (NumSamples x NumChannels)
%
%
%   SEE ALSO,
%   SQDHANDLE,GETDATA,GET,SET


% This function is just a wrapper around the methods of CLASS SQDHANDLE
% It uses the filename to make a sqdhandle_obj and then calls its methods
% to return the information asked for. Very little error correction is done
% here. It is assumed that the methods of SQDHANDLE will take care of it.

% Get filename
% If cannot find file, ask user input
if nargin<1,
    filenameqn = input(...
        'Do you wish to choose an sqd-file to open? [Return,Yes|No]',...
        's');
    if isempty(filenameqn)|strncmpi(filenameqn,'y',1),
        [f,p] = uigetfile({'*.sqd;*.SQD'},'Get SQD-FILE');
        if ((~isstr(f)) | ~min(size(f))), return;end;  % If Load cancelled
        filename = [p f];
    else
        error('Valid sqd-filename required');
    end;
end;

if ischar(filename)
    info = sqdhandle(filename);                     % Get sqdhandle_obj from file
    switch nargin-1                                 
        case 0                                      % Ex. data = sqdread(file.sqd);
            varargout{1} = getdata(info);           % Get data from sqdhandle_obj
            if nargout>=2                           % Ex. [data,info] = sqdread(file.sqd);
                varargout{2} = info;                % Return sqdhandle_obj as info
            end;
        case 1
            if strcmpi(varargin{1},'info')          % Ex. info = sqdread(file.sqd,'info');
                varargout{1} = info;                % Return sqdhandle_obj as info
            else
                error('Incorrect inputs')
            end;            
        otherwise        
            % Ex. data = sqdread(file.sqd,'Channels',[1 2],'Samples',[1 1000]);
            varargout{1} = getdata(info,varargin{:});
            if nargout>=2
                varargout{2} = info;
            end;
    end;
else
   error('First argument must be a valid filename'); 
end;
