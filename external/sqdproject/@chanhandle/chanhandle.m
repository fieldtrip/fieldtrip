function t = chanhandle(varargin)
% T = CHANHANDLE; Constructor for class CHANHANDLE 
% Chanhandle class contains channel-specific information only. 
% An object belonging to chanhandle has the following properties:
%   ChannelNumber
%   Sensitivity
%   SensorInfo
% Usage:
% t = chanhandle; Default values
% t = chanhandle(t1); where t1 is a chanhandle
% t = chanhandle(fname,channelnum); where fname is name of the sqd file
% See also,
% @chanhandle/get.m,set.m,writechanloc.m

% Note:
% Object structure:
% Chanhandle can be a handle to only one channel, or several channels.
% To accomodate an object array, the object is structured as follows:
% Chanhandle_obj has only one field - "handle". Subfields of "handle" are
% the above named properties. ie
% ChanHandle_Obj ->
%       handle   ->
%           ChannelNumber
%           Sensitivity
%           SensorInfo
% When the ChanHandle_obj contains several channels, "handle" is an array of
% structures.
%
% Also, note that:
% The ChannelNumber subfield indicates the hardware channel number. Whereas, the
% index of handle is software index to it. The index of handle is monotonically
% increasing from 1. So Chanhandle_obj contains HwChannels = 0,1,3
% then  handle(1) -> HwChannel = 0
%       handle(2) -> HwChannel = 1
%       handle(3) -> HwChannel = 3
% ChanHandle_obj is be indexed as the "handle" substructure.
%
% For internal use only, an additional way of initialising object:
%   t = chanhandle(handle); where "handle" is a structure with the channel's info
%       Fields for "handle": ChannelNumber,Sensitivity,SensorInfo


% initialize the fields of the object to default values
t.handle.ChannelNumber             = 0; % hardware channel number (0-191)
t.handle.Sensitivity               = []; % Sensitivity Gain/Offset
t.handle.SensorInfo                = []; % type,position,size,etc. of sensor

% Define class definition
t = class(t,'chanhandle');
switch nargin
case 0
    % Default object, return
case 1
    if isa(varargin{1},'chanhandle'), % if argin = chanhandle_obj
        t = varargin{1};
    elseif ischar(varargin{1})        % if argin = sqd-filename
        fname = varargin{1};
        t.handle = readsqdinfo(t.handle,varargin{1});
    elseif isstruct(varargin{1})      % if argin = "handle" substructure
        t.handle = varargin{1};
    else
        error('Incorrect input to chanhandle.m');
    end;
case 2                                  % if argin = sqd-filename,channelnumber
    if ischar(varargin{1})&isnumeric(varargin{2})
        for i = 1:length(varargin{2})
            t.handle(i).ChannelNumber = varargin{2}(i);
            if ~isempty(varargin{1})    % In case no file name is passed - default values
                t.handle(i) = readsqdinfo(t.handle(i),varargin{1});
            end;
        end;
    else
        error('Incorrect inputs to chanhandle.m');
    end;
otherwise
    error('Incorrect inputs to chanhandle.m');
end;

function t = readsqdinfo(t,fname)
% T = READSQDINFO(T,FNAME); Update values of CHANHANDLE object 'T' 
% to those given in sqd-file FNAME

seekset=-1; % Start point for fseek: -1 => beginning of the file

if nargin<2|~ischar(fname),error('Incorrect arguments');end;

% open file as binary and little endian
fid = fopen(fname,'rb','l');
if fid==-1,
    disp('File not found');
    [f,p] = uigetfile('*.sqd','MEG160 files');
    if ((~isstr(f)) | ~min(size(f))), 
        t = [];
        return; 
    end;
    fname = [p f];
    t.FileName = [p f];
    fid = fopen(fname,'rb','l');
end;

%Read in info regarding sensitivity of sensors
t.Sensitivity   = sensitivityinfo(fid,seekset,t.ChannelNumber);
% Read in sensor location
t.SensorInfo    = sensorinfo(fid,seekset,t.ChannelNumber);

% Close file
fclose(fid);