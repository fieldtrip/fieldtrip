function event_table = getMAFEvent(this, maf_file)
% GETMAFEVENT Get event table from MAF file
%   
% Syntax:
%   event_table = getMAFEvent(this, maf_file)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
%   maf_file        - [str] filepath + filename of MAF file
% 
% Output(s)
%   event_table     - [table] event table including Type (string) and
%                     Latency (in uUTC)
% 
% See also .

% Copyright 2019 Richard J. Cui. Created: Mon 05/27/2019  5:08:20.075 PM
% $Revision: 0.5 $  $Date: Sun 06/02/2019  2:51:12.672 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% main
% =========================================================================
% parse inputs
% ------------
q = parseInputs(this, maf_file);
maf_file = q.maf_file;
[fp, fn, ext] = fileparts(maf_file);
this.MafFilePath = fp;
this.MafFileName = [fn, ext];

maf = xml2struct(maf_file);
this.MAF = maf;

% get events
% ----------
maf_episode = maf.XREDE.Dataset.Subject.Episode;
time_unit = maf_episode.Attributes.time_units;
num_mafevents = numel(maf_episode.Event);
event = {};
var_names = {'Type', 'Latency'};
for k = 1:num_mafevents
    maf_event_k = maf_episode.Event{k};
    
    maf_eventtype_k = maf_event_k.Attributes.type;
    switch lower(maf_eventtype_k)
        case 'seizure'
            % seizure onset
            event_type = 'SeizureOn';
            event_latency = str2double(maf_event_k.Timestamp{1}.Attributes.onset);
            event_k = seizureEvent(this, event_type, event_latency, time_unit);
            event = cat(1, event, event_k);
            
            % seizure offset
            event_type = 'SeizureOff';
            event_latency = str2double(maf_event_k.Timestamp{1}.Attributes.offset);
            event_k = seizureEvent(this, event_type, event_latency, time_unit);
            event = cat(1, event, event_k); 
            
        otherwise
            fprintf('Unknown event type %s.\n', maf_eventtype_k)
            
    end % switch
    
end % for
event_table = cell2table(event, 'VariableNames', var_names);

% sort event_table according to latency
% -------------------------------------
[~, lat_index] = sort(event_table.Latency);
event_table = event_table(lat_index, :);

end

% =========================================================================
% subroutines
% =========================================================================
function seizure_e = seizureEvent(this, event_type, event_latency, time_unit)

if ~strcmpi(time_unit, 'uutc')
    event_latency_ind = this.SampleTime2Index(event_latency, time_unit);
    event_latency = this.SampleIndex2Time(event_latency_ind, 'uutc');
end % if
seizure_e = {event_type, event_latency};

end % function

function q = parseInputs(varargin)

% defaults

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('maf_file', @isstr);

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]