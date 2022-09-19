function filename = bst_fullfile(filename, varargin)
% BST_FILEPARTS: Same as Matlab's fullfile, but with much less tests.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2013

% Look for existing '\' in the filename
if any(filename == '\')
    sep = '\';
    otherSep = '/';
else
    sep = '/';
    otherSep = '\';
end
% Remove trailing separator
if isempty(filename)
    filename = '';
elseif (filename(end) == sep)
    filename = filename(1:end-1);
end
% Add all the inputs
for i = 1:length(varargin)
    p = varargin{i};
    % Empty: skip
    if isempty(p)
        continue;
    end
    % Replace with the correct separator
    p(p == otherSep) = sep;
    % If the new element starts with the separator: do not add an extra one
    if (p(1) == sep) || isempty(filename)
        filename = [filename p];
    % Add a separator between the base path and the new element
    else
        filename = [filename sep p];
    end
end
