function F = in_fread_manscan(sFile, sfid, iEpoch, SamplesBounds)
% IN_FREAD_MANSCAN:  Read a block of recordings from a MANSCAN file
%
% USAGE:  F = in_fread_manscan(sFile, sfid, iEpoch, SamplesBounds)  : Read all channels
%         F = in_fread_manscan(sFile, sfid)                         : Read all channels, all the times

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
% Authors: Francois Tadel, 2012-2014

if (nargin < 3) || isempty(iEpoch)
    iEpoch = 1;
end
if (nargin < 4) || isempty(SamplesBounds)
    if ~isempty(sFile.epochs)
        SamplesBounds = round(sFile.epochs(iEpoch).times .* sFile.prop.sfreq);
    else
        SamplesBounds = round(sFile.prop.times .* sFile.prop.sfreq);
    end
end

% Read data block
nChannels  = length(sFile.header.epoch(1).ChannelOrder);
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





