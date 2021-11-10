function varargout = bst_fileparts(filename, isFolder)
% BST_FILEPARTS: Same as Matlab's fileparts, but consider '\' as file separators on Linux/MacOS systems.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2018 University of Southern California & McGill University
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
% Authors: Francois Tadel, 2010

% Parse inputs
if (nargin < 2) || isempty(isFolder)
    isFolder = 0;
end
% MacOS and Linux: Replace \ with /
if ~ispc
    filename(filename == '\') = '/';
end
% Check if link
if (length(filename) > 5) && strcmpi(filename(1:5), 'link|')
    if (nargout >= 1)
        varargout{1} = 'link';
    end
    if (nargout >= 2)
        varargout{2} = filename;
    end
    if (nargout >= 3)
        varargout{3} = [];
    end
% Folder: split only with str_split
elseif isFolder
    % Replace dots with invalid characters, to avoid splitting the folder names that include dots
    filename = strrep(filename, '.', '>');
    % Split filename
    [varargout{1:nargout}] = fileparts(filename);
    % Restore the dots
    for i = 1:length(varargout)
        varargout{i} = strrep(varargout{i}, '>', '.');
    end
% Regular file
else 
    [varargout{1:nargout}] = fileparts(filename);
end




