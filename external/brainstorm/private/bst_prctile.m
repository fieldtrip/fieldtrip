function value = bst_prctile(vector, percentile)
% BST_PRCTILE: Returns the percentile value in vector
%
% USAGE: value = bst_prctile(vector, percentile)

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
% Authors: Martin Cousineau, 2020

% Try to use toolbox function
try
    value = prctile(vector, percentile);
    return;
catch
end

if ~isvector(vector)
    error('Only vectors supported.');
end

% Custom implementation
vector = sort(vector);
rank   = percentile / 100 * (length(vector) + 1);
lowerRank = floor(rank);
upperRank = ceil(rank);
fraction  = rank - lowerRank;

if fraction == 0
    value = vector(rank);
else
    value = fraction * (vector(upperRank) - vector(lowerRank)) + vector(lowerRank);
end
