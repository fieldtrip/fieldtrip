% read_erplabheader() - import ERPLAB dataset files
%
% Usage:
%   >> header = read_erplabheader(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   header   - FILEIO toolbox type structure
%
% Modified from read_eeglabheader

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function header = read_erplabheader(filename)

if nargin < 1
  help read_erplabheader;
  return;
end;

if ~isstruct(filename)
  load('-mat', filename);
else
  ERP = filename;
end;

header.Fs          = ERP.srate;
header.nChans      = ERP.nchan;
header.nSamples    = ERP.pnts;
header.nSamplesPre = -ERP.xmin*ERP.srate;
header.nTrials     = ERP.nbin;
try
  header.label       = { ERP.chanlocs.labels }';
catch
  ft_warning('creating default channel names');
  for i=1:header.nChans
    header.label{i} = sprintf('chan%03d', i);
  end
end
ind = 1;
for i = 1:length( ERP.chanlocs )
    if isfield(ERP.chanlocs(i), 'X') && ~isempty(ERP.chanlocs(i).X)
        header.elec.label{ind, 1} = ERP.chanlocs(i).labels;
        % this channel has a position
        header.elec.pnt(ind,1) = ERP.chanlocs(i).X;
        header.elec.pnt(ind,2) = ERP.chanlocs(i).Y;
        header.elec.pnt(ind,3) = ERP.chanlocs(i).Z;
        ind = ind+1;
    end;
end;

header.orig = ERP;

