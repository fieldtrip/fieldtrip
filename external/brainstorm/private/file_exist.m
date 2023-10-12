function isExist = file_exist(filename)
% FILE_EXIST: Test the existence of an absolute file path on the disk.
%             Replacement for the Matlab call "exist(..., 'file')" that gives unexpected results
%             when the file is somewhere in the Matlab path.

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
% Author: Francois Tadel, 2011-2013

% Empty variable
if isempty(filename)
    isExist = 0;
% File is a root folder: let's consider it exists
% elseif strcmp(filename, bst_fileparts(filename))
elseif isequal(filename, '/') || isequal(filename, '\')
    isExist = 1;
% File does not exist: no ambiguity
elseif isempty(dir(filename))
    isExist = 0;
% File exists: we have to make sure it is NOT relatively to the current path
else
    % Check the file relative to the current Matlab folder
    pwdPath = pwd;
    pwdFile = [pwdPath '/' filename];
    % If it does NOT exist: no ambiguity, file exists
    if isempty(dir(pwdFile))
        isExist = 1;
    % Else: we consider that the file exists only in the current folder is a root folder
    % (=> folder = folder's parent)
    elseif strcmp(pwdPath, bst_fileparts(pwdPath, 1))
        isExist = 1;
    else
        isExist = 0;
    end
end


