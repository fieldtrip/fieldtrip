function timestamps = edgeDetect(signal, threshold, type)

% edgeDetect
% 
% Detects rising or falling edges of a signal whenever the signal crosses
% the threshold.
%
%
% Use timestamps = edgeDetect(signal, threshold, type)
% 
% Variables fname and threshold are optional.
%
%   signal:       Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file. 
%
%   threshold:    The value used for threshold crossings. This function
%                 will detect when the signal crosses this value.
%                 DEFAULT: It will automatically choose the threshold at
%                 80% of the maximum value in the signal.
%
%
%   type:         Contains the points in the signal where the threshold
%                 crossing occurs.
%                 DEFAULT: It will find the rising edges.
%
%   Example 1: 
%   timestamps = edgeDetect(signal, 1000);
%
%   In the example above, the points where signal crosses value 1000 
%   (rising edge) is detected and returned in variable timestamps.
%
%   timestamps = edgeDetect(signal, 1000, 'falling');
%
%   In the example above, the points where signal crosses value 1000 
%   (falling edge) is detected and returned in variable timestamps.
%
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.1.0.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.1.0.0:
%   - Added automatic edge detection.
%   - Updated help.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validations

% Validating input arguments.
if nargin < 1 || nargin > 3
    disp('Invalid number of input arguments. Use ''help edgeDetect'' for more information.');
    return;
end

if ~exist('threshold', 'var')
    threshold = 0.8 * max(signal);
    disp(['Threshold was not provided. It wss automatically calculated and set at ' num2str(threshold) '.']);
end

% Validating variable 'type'
if ~exist('type', 'var')
    type = 'rising';
end

% Validating type and determining the threshold crossing points
if strcmpi(type, 'rising')
    timestamps = signal>threshold;
elseif strcmpi(type, 'falling')
    timestamps = signal<threshold;
else
    disp('Type does not exist. Please type ''help edgeDetect'' to see all available types.');
    timestamps = 0;
    return;
end

% Finding all the points where the signal crosses the threshold
timestamps = diff(timestamps);
timestamps = find(timestamps==1);