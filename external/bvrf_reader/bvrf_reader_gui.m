% bvrf_reader_gui() - BVRF Reader UI
% Takes an structure with fields specified below and launch a GUI for the
% reader. 
%
% Usage:
%   >> [cfg, result] = bvrf_reader_gui(info);
%
% Inputs:
%  info      - Structure with following fields
%              Where info:
%              info.bvrfVersion           -> Version of BVRF file
%              info.dataType              -> Data type of data         
%              info.samples               -> Length on samples of the recording
%              info.srate                 -> Sampling rate of the recording
%              info.participantNames      -> Number of participants in the recording
%              info.participantChans      -> Number of channels per participants
%              info.hasImpedances         -> Flag if impedance file is present
%
% Outputs:
%   cfg             - Structure with user selection turned into options
%   cfg.participantId           -> Participant Id
%   cfg.sampleInterval          -> Samples interval
%   cfg.channelIndx             -> Indices of channels to import
%   cfg.flagImportMarkers       -> Boolean. True: import Markers information 
%   cfg.flagImportImpedances    -> Boolean. True: import Impedance information
%
% Author: Ramon Martinez-Cancino, Brain Products GmbH, 2025
%
% Copyright (C) 2025 Brain Products GmbH
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function cfg = bvrf_reader_gui(info)

% ---------- safety / defaults ----------
if nargin < 1
    error('bvrf_reader_gui: info struct required');
end

nPart = numel(info.participantNames);
partNames = info.participantNames;
partChans = info.participantChans;
% end

% ---------- strings for header ----------
line1_left  = sprintf('BVRF: %s', info.bvrfVersion);
line1_right = sprintf('Type: %s', info.dataType);
line2_left  = sprintf('Samples: %d', info.samples);
line2_right = sprintf('Sampling rate: %.4g Hz', info.srate);
line3       = sprintf('Participants: %d', nPart);

% "table" of participants (Subject / channel count)
partListStr = cell(nPart,1);
for k = 1:nPart
    if k <= numel(partChans)
        ch = partChans(k);
    else
        ch = NaN;
    end
    partListStr{k} = sprintf('%s    %g%s', partNames{k}, ch, '(Channels)');
end

% popup strings
popupPartStr  = ['All Participants|' strjoin(partNames, '|')];
popupUnitsStr = 'Samples|Seconds';

% GEOMETRY 
% =========================================================
geometry = { ...
    [1]           ... % 1  "Dataset Info"
    [1 1]         ... % 2  BVRF / Type
    [1 1]         ... % 3  Samples / Fs
    [1]           ... % 4  Participants count
    [1]           ... % 5  Participants listbox
    [1]           ... % 5  Spacer
    [1]           ... %     Reading options
    [0.7 1.3]     ... % 6  Participants label + popup
    [1]           ... % 7  "Interval"
    [0.8 0.8 0.8] ... % 8  Start / End / Units labels
    [0.8 0.8 0.8] ... % 9  Start / End / Units edits / popup
    [1]           ... %10  "Channels"
    [1]           ... %11  Channels edit
    [1]           ... %12  "Additional data to import"
    [0.9 0.9]     ... %13  Markers / Impedances
    [1]           ... %14  Use Poly
    };

vertval =  1;
geomvert = [ ...
    vertval ...     % title
    vertval ...     % BVRF / Type
    vertval ...     % Samples / Fs
    vertval ...     % Participants count
    3 ...           % listbox
    vertval ...     % Spacer
    vertval ...     % Reading Options
    vertval+0.5...  % Participants popup row
    vertval+0.5...  % "Interval"
    vertval ...     % Start/End/Units labels
    vertval ...     % Start/End/Units controls
    vertval+0.5 ... % Channels label
    vertval ...     % Channels edit
    vertval+0.5 ... % "Additional data" label
    vertval ...     % checkboxes
    vertval ...     % checkbox poly
    ];

% UILIST  (order matters)
% =========================================================
uilist = {};

% 1: Dataset info title
uilist{end+1} = {'style','text','string','Dataset Info','fontweight','bold'};

% 2: BVRF / Type
uilist{end+1} = {'style','text','string',line1_left};
uilist{end+1} = {'style','text','string',line1_right};

% 3: Samples / Fs
uilist{end+1} = {'style','text','string',line2_left};
uilist{end+1} = {'style','text','string',line2_right};

% 4: Participants count
uilist{end+1} = {'style','text','string',line3};

% 5: Participants listbox
uilist{end+1} = {'style','listbox','string',partListStr, ...
                 'min',1,'max',1,'tag','participants_table'};
% 6: Separator
uilist{end+1} = {};

% 7: REading Options
uilist{end+1} = {'style','text','string','Reading options','fontweight','bold'};

% 8: Participants label line
uilist{end+1} = {'style','text','string','Participants','fontweight','bold'};

% 9: Participants popup line
uilist{end+1} = {'style','popupmenu','string',popupPartStr,'value',1, ...
                 'tag','popup_participants'};

% 10: Interval label
uilist{end+1} = {'style','text','string','Interval','fontweight','bold'};

% 11: Start / End / Units labels
uilist{end+1} = {'style','text','string','Start'};
uilist{end+1} = {'style','text','string','End'};
uilist{end+1} = {'style','text','string','Units'};

%12: Start / End / Units controls
uilist{end+1} = {'style','edit','string','','tag','edit_start'};
uilist{end+1} = {'style','edit','string','','tag','edit_end'};
uilist{end+1} = {'style','popupmenu','string',popupUnitsStr,'value',1, ...
                 'tag','popup_units'};

%13: Channels label
uilist{end+1} = {'style','text','string','Channels','fontweight','bold'};

%14: Channels edit
uilist{end+1} = {'style','edit','string','','tag','edit_channels'};

%15: Additional data label
uilist{end+1} = {'style','text','string','Additional data and calibrations','fontweight','bold'};

%16: Markers / Impedances checkboxes
uilist{end+1} = {'style','checkbox','string','Markers', ...
                 'value',info.hasMarkers,'tag','cb_markers'};
uilist{end+1} = {'style','checkbox','string','Impedances', ...
                 'value',info.hasImpedances,'tag','cb_impedances'};

uilist{end-1}{end+1} = 'enable';  uilist{end-1}{end+1} = 'on';

if ~info.hasImpedances
    uilist{end}{end+1}   = 'enable';  uilist{end}{end+1}   = 'off'; % disable box if data not present
end

%17: Use Poly checkboxes
uilist{end+1} = {'style','checkbox','string','Apply sensor calibration', ...
                 'value',true,'tag','cb_usepoly'};


% Call inputgui
% =========================================================
[result, ~, ~, ~] = inputgui( ...
    'geometry', geometry, ...
    'uilist',   uilist, ...
    'geomvert', geomvert, ...
    'title',    'BVRRF Reader', ...
    'helpbut', 'Help', ...
    'helpcom',  'pophelp(''pop_loadbvrf'')',...
    'minwidth', 1);

if isempty(result)   % user cancelled
    cfg = [];
    return;
end

% ---------- Map results to options ----------
cfg = struct();

% Particpant
if result{2} == 1
    cfg.participantId = [];
else
    cfg.participantId = partNames{result{2}-1};
end

% Units
unitList = strsplit(popupUnitsStr, '|');
intervalUnit = unitList{result{5}};

% Range Interval
if isempty(result{3}) || isempty(result{4})
    cfg.sampleInterval   = [];
else
    range = [str2num(result{3})  str2num(result{4})];
    if strcmp(intervalUnit, 'Samples')
        cfg.sampleInterval = range;
    else
        cfg.sampleInterval  = range*info.srate;
    end
end

% Channels
if isempty(result{6})
    cfg.channelIndx  = [];
else
    cfg.channelIndx  = str2num(result{6});
end

% Markers and Impedances
cfg.flagImportMarkers    = logical(result{7});
cfg.flagImportImpedances = logical(result{8}) && info.hasImpedances;

cfg.usePoly = logical(result{9});
end