% eegplugin_bvrfimport() - EEGLAB plugin to import BrainVision BVRF files
%
% Usage:
%   >> eegplugin_bvrfimport(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
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

function vers = eegplugin_bvrfimport(fig, trystrs, catchstrs)

vers = 'bvrfimport v1.0.0';

if nargin < 3
    error('eegplugin_bvrfimport requires 3 input arguments');
end

% Find the "Import data" menu in EEGLAB
hImport = findobj(fig, 'tag', 'import data');

% Add our menu item under "File -> Import data"
% The callback:
%   - Calls [EEGarray, com] = pop_loadbvrf;
%   - Expects EEGarray to be a struct array with one EEG per dataset;
%   - Stores each dataset in ALLEEG;
%   - Redraws EEGLAB and sets LASTCOM for history.

uimenu(hImport, 'Label', 'From BrainVision (BVRF)...', ...
    'Callback', [ ...
    trystrs.no_check ...
    '[EEGarray, com] = pop_loadbvrf; ' ...                           % call main function
    'if ~isempty(EEGarray),' ...
    'if ~iscell(EEGarray), EEGarray = num2cell(EEGarray); end;' ... % ensure cell of EEG structs
    'for k = 1:numel(EEGarray),' ...
    'EEG = EEGarray{k};' ...                                        % take one EEG at a time
    '[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);' ...
    'end;' ...
    'eeglab(''redraw'');' ...
    'end;' ...
    'LASTCOM = com;' ...
    catchstrs.add_to_hist ...
    ]);
end
