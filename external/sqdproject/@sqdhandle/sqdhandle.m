function t = sqdhandle(varargin)
% T = SQDHANDLE(SQD_FILENAME); SQDHANDLE Constructor
% 
% Usage:
% t = sqdhandle; Default values
% t = sqdhandle(t1); where t1 is a sqdhandle
% t = sqdhandle(filename); where filename is name of the sqd-file
%
% SQDHANDLE saves all the information from the sqd-file and saves it 
% in a convenient object.
% Please see other methods of SQDHANDLE: 
%   GET.M, SET.M, GETDATA.M, PUTDATA.M

% Note: 
% Object Structure:
%   Version
%   Revision
%   Systemname
%   ModelName
%   Comment
%   FileName
%   PatientInfo - a structure with following fields:
%       id
%       Name
%       birthdate
%       handedness
%       gender
%   Amplifier - a structure with the foll. fields
%	    InputGainMask
%	    InputGainBit
%	    InputGain
%	    OutputGainMask
%	    OutputGainBit
%	    OutputGain
%   Channels
%   ChannelCount
% Inherited Properties
%   SQDHANDLE inherits from CLASS ACQPARAM
% Aggregation
%   SQDHANDLE aggregates CLASS CHANHANDLE into field 'Channels'

% init object structure
t.Version       = 1;
t.Revision      = 1;
t.SystemID      = '';
t.SystemName    = '';
t.ModelName     = '';
t.FileName      = '';
t.Amplifier     = [];
t.ChannelCount	= 192;
t.Channels      = chanhandle(t.FileName,0);
t.Comment       = '';
t.PatientInfo   = [];
acqtype	        = acqparam;
t = class(t,'sqdhandle',acqtype);

switch nargin
    case 0
        % Default object, return
    case 1
	if isa(varargin{1},'sqdhandle'),% argin = sqdhandle_obj
	    t = varargin{1};
	elseif ischar(varargin{1})      % argin = sqd-filename
	    t.FileName = varargin{1};
	    t = readsqdinfo(t,varargin{1});
	else
	    error('Incorrect input');
	end;
end;
