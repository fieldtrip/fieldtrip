function F = in_fread_manscan(sFile, sfid, iEpoch, SamplesBounds)
% IN_FREAD_MANSCAN:  Read a block of recordings from a MANSCAN file
%
% USAGE:  F = in_fread_manscan(sFile, sfid, iEpoch, SamplesBounds)  : Read all channels
%         F = in_fread_manscan(sFile, sfid)                         : Read all channels, all the times

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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
% Authors: Francois Tadel, 2012

if (nargin < 3) || isempty(iEpoch)
    iEpoch = 1;
end
if (nargin < 4) || isempty(SamplesBounds)
    if ~isempty(sFile.epochs)
        SamplesBounds = sFile.epochs(iEpoch).samples;
    else
        SamplesBounds = sFile.prop.samples;
    end
end

% Read data block
nChannels  = length(sFile.channelmat.Channel);
nTimes = SamplesBounds(2) - SamplesBounds(1) + 1;
% Epoch offset
epochOffset = sFile.header.epoch(iEpoch).StartData1;
% Time offset
timeOffset = 2 * SamplesBounds(1) * nChannels;
% Total offset
totalOffset = epochOffset + timeOffset;

% Set position in file
fseek(sfid, totalOffset, 'bof');
% Read value
F = fread(sfid, [nChannels,nTimes], 'int16');

% Apply gains
F = bst_bsxfun(@rdivide, double(F), double(sFile.header.Gains));


function C = bst_bsxfun(fun, A, B)
% BST_BSXFUN:  Compatible version of bsxfun function.
% 
% DESCRIPTION:
%     Matlab function bsxfun() is a useful, fast, and memory efficient 
%     way to apply element-by-element operations on huge matrices. 
%     The problem is that this function only exists in Matlab versions >= 7.4
%     This function check Matlab version, and use bsxfun if possible, if not
%     it finds another way to perform the same operation.

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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

% Old Matlab version: do it old school
if ~exist('bsxfun', 'builtin')
    sA = [size(A,1), size(A,2), size(A,3)];
    sB = [size(B,1), size(B,2), size(B,3)];
    % If arguments were not provided in the correct order
    if all(sA == [1 1 1]) || all(sB == [1 1 1])
        C = fun(A, B);
    elseif all(sA == sB)
        C = fun(A, B);
    % Dim 1
    elseif (sB(1) == 1) && (sA(2) == sB(2)) && (sA(3) == sB(3))
        C = fun(A, repmat(B, [sA(1), 1, 1]));
    elseif (sA(1) == 1) && (sA(2) == sB(2)) && (sA(3) == sB(3))
        C = fun(repmat(A, [sB(1), 1, 1]), B);
    % Dim 2
    elseif (sA(1) == sB(1)) && (sB(2) == 1) && (sA(3) == sB(3))
        C = fun(A, repmat(B, [1, sA(2), 1]));
    elseif (sA(1) == sB(1)) && (sA(2) == 1) && (sA(3) == sB(3))
        C = fun(repmat(A, [1, sB(2), 1]), B);
    % Dim 3
    elseif (sA(1) == sB(1)) && (sA(2) == sB(2)) && (sB(3) == 1)
        C = fun(A, repmat(B, [1, 1, sA(3)]));
    elseif (sA(1) == sB(1)) && (sA(2) == sB(2)) && (sA(3) == 1)
        C = fun(repmat(A, [1, 1, sB(3)]), B);
    else
        error('A and B must have enough common dimensions.');
    end
% New Matlab version: use bsxfun
else
   C = bsxfun(fun, A, B);
end







