function D = meeg(varargin)
% Function for creating meeg objects.
% FORMAT m = meeg(varargin)
%
% SPM8 format consists of a binary memory - mapped data file and a header 
% mat file.
%
% The header file will contain an object of class meeg called D. All 
% information other than data is contained in this object and access to the
% data is via methods of the object. Also, arbitrary data can be stored 
% inside the object if their field names do not conflict with existing 
% methods' names.
%
% The following is a description of the internal implementation of meeg.
%
% Fields of meeg:
% .type - type of data in the file: 'continuous', 'single', 'evoked'
% .Fsample - sampling rate
% .data - substruct containing the information related to the bimary data 
% file and memory mapping.
%
%   Subfields of .data
%       .y - reference to the data file_array
%       .fnamedat - name of the .dat file
%       .datatype - type of the data for file_array
%
%
% .Nsamples - length of the trial (whole data if the file is continuous).
% .timeOnset - the peri-stimulus time of the first sample in the trial (in sec)
%
% .fname, .path - strings added by spm_eeg_load to keep track of where a 
%                 header struct was loaded from.
%
% .trials - this describes the segments of the epoched file and is also a 
%           structure array.
%
%   Subfields of .trials
%
%       .label - user-specified string for the condition
%       .onset - time of the first sample in seconds in terms of the 
%                original file
%       .bad - 0 or 1 flag to allow rejection of trials.
%       .repl - for epochs that are averages - number of replications used 
%               for the average.
%       .events - this is a structure array describing events related to 
%                 each trial.
%
%           Subfields of .events
%
%           .type - string (e.g. 'front panel trigger')
%           .value - number or string, can be empty (e.g. 'Trig 1').
%           .time - in seconds in terms of the original file
%           .duration - in seconds
%
% .channels - This is a structure array which is a field of meeg.
%             length(channels) should equal size(.data.y, 1) and the order 
%             must correspond to the order of channels in the data.
%
%   Subfields of .channels
%
%       .label - channel label which is always a string
%       .type - a string, possible values - 'MEG', 'EEG', 'VEOG', 'HEOG', 
%               'EMG' ,'LFP' etc.
%       .units - units of the data in the channel.
%       .bad - 0 or 1 flag to mark channels as bad.
%       .X_plot2D, .Y_plot2D - positions on 2D plane (formerly in ctf). NaN
%                              for channels that should not be plotted.
%
% .sensors
%
%
%   Subfields of .sensors (optional)
%       .meg - struct with sensor positions for MEG (subfields: .pnt .ori .tra .label)
%       .eeg - struct with sensor positions for MEG (subfields: .pnt .tra .label)
%
% .fiducials - headshape and fiducials for coregistration with sMRI
%      
%   Subfiels of .fiducials (optional)
%       .pnt - headshape points
%       .fid.pnt - fiducial points
%       .fid.label - fiducial labels
%
% .artifacts - structure array with fields .start and .stop expressed in 
%              seconds in original file time.
%
% .history - structure array describing commands that modified the file.
%
%   Subfields of .history:
%
%       .function - string, the function name
%       .arguments - cell array, the function arguments
%       .time - when function call was made
%
% .other - structure used to store other information bits, not fitting the
%          object structure at the moment,
%       for example:
%       .inv - structure array corresponding to the forw/inv problem in MEEG.
%       .ival - index of the 'inv' solution currently used.
%
% ______________________________________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: meeg.m 3940 2010-06-21 12:51:19Z christophe $

if nargin == 1
    if isstruct(varargin{1})
        [OK D] = checkmeeg(varargin{1}, 'basic');
        if OK
            D  = class(D, 'meeg');
        else
            error('Struct not fit for conversion to meeg');
        end
        return;
    elseif isa(varargin{1},'meeg'),
        D = varargin{1};
        return;
    end;
end;

D          = [];
D.Nsamples = 0;
[res D]    = checkmeeg(D);
D          = class(D, 'meeg');
