function dc_event = findDiscontEvent(this, start_end, unit)
% MEFFIELDTRIP.FINDDISCONTEVENT process dicountinuity events
% 
% Syntax:
%   dc_event = findDiscontEvent(this)
%   dc_event = findDiscontEvent(__, start_end)
%   dc_event = findDiscontEvent(__, start_end, unit)
% 
% Input(s):
%   this            - [obj] MEFEEGLab_2p1 object
%   start_end       - [1 x 2 array] (optional) [start time/index, end time/index] of 
%                     the signal to be extracted from the file (default:
%                     the entire signal)
%   unit            - [str] (optional) unit of start_end: 'Index', 'uUTC',
%                     'Second', 'Minute', 'Hour', and 'Day' (default = '')
% Output(s):
%   dc_event        - [struct] discountinuity event structure
% 
% Note:
% 
% See also .

% Copyright 2020 Richard J. Cui. Created: Sat 03/21/2020 10:35:23.147 PM
% $Revision: 0.1 $  $Date: Sat 03/21/2020 10:35:23.147 PM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, start_end, unit);
start_end = q.start_end;
if isempty(start_end)
    start_end = this.StartEnd;
end % if
unit = q.unit;
if isempty(unit)
    unit = this.SEUnit;
end % if

% =========================================================================
% main process
% =========================================================================
% converte start_end to index if not
% if strcmpi(unit, 'index')
%     se_ind = start_end;
% else
%     se_ind = this.SampleTime2Index(start_end, unit);
% end % if
se_ind = this.SessionUnitCovert(start_end, unit, 'index');

% find the continuity blocks
seg_cont = this.SessionContinuity;
cont_ind = se_ind(1) <= seg_cont.SampleIndexEnd...
    & se_ind(2) >= seg_cont.SampleIndexStart;
dc_start = seg_cont.SampleIndexStart(cont_ind); % discont start in index

% find the relative index of start of discontinuity
rel_dc = dc_start - se_ind(1);
rel_dc(rel_dc < 0) = []; % get rid of index < 0

% construct the event of EEGLAB
num_event = numel(rel_dc);
t = table('Size', [num_event, 3], 'VariableTypes', {'string', 'double', 'double'},...
    'VariableNames', {'type', 'latency', 'urevent'});
t.type = repmat('Discont', num_event, 1);
t.latency = rel_dc(:);
t.urevent = (1:num_event)';

% output
dc_event = table2struct(t);

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
defaultSE = [];
defaultUnit = '';
expectedUnit = {'index', 'uutc', 'second', 'minute', 'hour', 'day'};

% parse rules
p = inputParser;
p.addRequired('this', @(x) isobject(x) || strcmpi(class(x), 'MEFEEGLab_2p1'));
p.addOptional('start_end', defaultSE,...
    @(x) isnumeric(x) & numel(x) == 2 & x(1) <= x(2));
p.addOptional('unit', defaultUnit,...
    @(x) any(validatestring(x, expectedUnit)));

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]