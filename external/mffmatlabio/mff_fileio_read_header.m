% mff_fileio_read_header() - import MFF header
%
% Usage:
%   >> header = read_eeglabheader(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   header   - FILEIO toolbox type structure

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function header = mff_fileio_read_header(mffFile)

if nargin < 1
    help mff_fileio_read_header;
    return;
end

% import basic information
header.orig = mff_import(mffFile);
header.Fs          = header.orig.srate;
header.nChans      = header.orig.nbchan;
header.nSamples    = header.orig.pnts;

% Import epoch information if any
header.nTrials     = header.orig.trials;
header.nSamplesPre = -header.orig.xmin*header.orig.srate;
header.nSamples    = header.orig.pnts;

% import coordinate layout in EEGLAB format
chanlocs = header.orig.chanlocs;
header.elec.pnt   = zeros(length( chanlocs ), 3);
for ind = 1:length( chanlocs )
    header.elec.label{ind} = chanlocs(ind).labels;
    if ~isempty(chanlocs(ind).X)
        header.elec.pnt(ind,1) = chanlocs(ind).X;
        header.elec.pnt(ind,2) = chanlocs(ind).Y;
        header.elec.pnt(ind,3) = chanlocs(ind).Z;
    else
        header.elec.pnt(ind,:) = [0 0 0];
    end
end
if isfield(header, 'elec') && isfield(header.elec, 'pnt') && isempty(header.elec.pnt)
    header = rmfield(header, 'elec');
end

if isfield(header, 'elec') && isfield(header.elec, 'label')
    header.label = header.elec.label;
else
    header.label = {};
end
