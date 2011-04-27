function t = acqparam(varargin)
% T = ACQPARAM(SQD_FILENAME); ACQPARAM CLASS CONSTRUCTOR
% ACQPARAM initializes acqparam object which holds only information pertaining
% to the data-acquisition by MEG160 which is common to each channel
% Usage:
%   t = acqparam; Default values
%   t = acqparam(t1); where t1 is an object of class acqparam
%   t = acqparam(filename,acq_offset,raw_offset); where 'filename' is the 
%       name of the sqd file and acq_offset and raw_offset are offsets
%       into the file to read acquisition info 
%   t = acqparam(fid,acq_offset,raw_offset); where fid is the file pointer
%       of a sqd file
% See also,
% @ACQPARAM/get.m,set.m

% Note: 
% This class was made more for convenience so that there would be
% independent functions which would take care of the different
% acquisition types.
%
% Object structure
% The object acqparam has the following fields:
%       SampleRate          Data acquisition sampling rate
%       AcquisitionType     non-triggered/triggered/averaged
%       Offsets             Location of the start of the data
%       Datatype            dataformat in file
%       SampleInfo          Number of samples,etc
%           Depending on type of acquisition:
%            - SamplesAcquired,ActSamplesAcquired,SamplesAvailable
%            - FrameLength,PreTrigger,AverageCount,ActAverageCount
%       Private fields:
%           Offset2Raw
%           Datatypes

% Initialize default values 
t.SampleRate          = 1000;
t.AcquisitionType     = 1;
t.Private.Offset2Raw  = [144;160;144]; % Offsets to read raw_offsets
t.RawOffset           = 0;         % Actual offsets to info [acq;raw]
t.Private.Datatypes   = {'int16';'double';'int16'};
t.Datatype            = t.Private.Datatypes{t.AcquisitionType};
t.SampleInfo          = [];

% Class definition
t = class(t,'acqparam');

switch nargin
case 0
case 1
    if isa(varargin{1},'acqparam')  % argin = acqparam object
        t = varargin{1};        
    else      % argin = filename or filepointer
        t = readsqdinfo(t,varargin{1});
    end;
case {2,3}                          % argins = fid,acq_offset,raw_offset
    t = readsqdinfo(t,varargin{:}); %        = filename,acq_offset,raw_offset
end;