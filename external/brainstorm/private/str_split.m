function splStr = str_split( str, delimiters, isCollapse )
% STR_SPLIT: Split string.
%
% USAGE:  str_split( str, delimiters, isCollapse=1 ) : delimiters in an array of char delimiters
%         str_split( str )             : default are file delimiters ('\' and '/')
%
% INPUT: 
%    - str        : String to split
%    - delimiters : String that contains all the characters used to split, default = '/\'
%    - isCollapse : If 1, remove all the empty entries
% 
% OUTPUT: 
%    - splStr : cell array of blocks found between separators

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
% Authors: Francois Tadel, 2008-2014

% Default delimiters: file delimiters ('\', '/')
if (nargin < 3) || isempty(isCollapse)
    isCollapse = 1;
end
if (nargin < 2) || isempty(delimiters)
    delimiters = '/\';
end
% Empty input
if isempty(str)
    splStr = {};
    return
end

% Find all delimiters in string
iDelim = [];
for i=1:length(delimiters)
    iDelim = [iDelim strfind(str, delimiters(i))];
end
iDelim = unique(iDelim);

% If no delimiter: return the whole string
if isempty(iDelim)
    splStr = {str};
    return
end

% Allocates the split array
splStr = cell(1, length(iDelim)+1);
    
% First part (before first delimiter)
if (iDelim(1) ~= 1)
    iSplitStr = 1;
    splStr{iSplitStr} = str(1:iDelim(1)-1);
else
    iSplitStr = 0;
end
    
% Loop over all other delimiters
for i = 2:length(iDelim)
    if (isCollapse  && (iDelim(i) - iDelim(i-1) > 1)) || ...
       (~isCollapse && (iDelim(i) - iDelim(i-1) >= 1))
        iSplitStr = iSplitStr + 1;
        splStr{iSplitStr} = str(iDelim(i-1)+1:iDelim(i)-1);
    end
end
    
% Last part (after last delimiter)
if (iDelim(end) ~= length(str))
    iSplitStr = iSplitStr + 1;
    splStr{iSplitStr} = str(iDelim(end)+1:end);
end
    
% Remove all the unused entries
if (iSplitStr < length(splStr))
    splStr(iSplitStr+1:end) = [];
end






    