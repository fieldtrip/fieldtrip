% read_eeglabheader() - import EEGLAB dataset files
%
% Usage:
%   >> header = read_eeglabheader(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   header   - FILEIO toolbox type structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2008-

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

% $Log: read_eeglabheader.m,v $
% Revision 1.6  2009/08/08 04:07:02  josdie
% Bug in Joe's brain fixed.  Change to header.nSamplesPre calculation changed back.
%
% Revision 1.5  2009/08/08 03:17:26  josdie
% Fixed bug that was causing hdr.label to have as many labels as there are time points rather than matching the number of channels.
%
% Revision 1.4  2009/08/08 03:05:29  josdie
% Fixed bug in calculation of header.nSamplesPre.
%
% Revision 1.3  2009/07/01 16:08:21  vlalit
% Fixing a bug in converting channel locations to elec struct (reproted by Jakib Scherer)
%
% Revision 1.2  2009/01/23 15:35:46  roboos
% create default channel names if EEG.chanlocs.labels is missing
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.2  2008/04/21 18:45:59  roboos
% fixed bug, ori should be orig
%
% Revision 1.1  2008/04/18 14:04:48  roboos
% new implementation by Arno, shoudl be tested
%

function header = read_eeglabheader(filename)

if nargin < 1
  help read_eeglabheader;
  return;
end;

if ~isstruct(filename)
  load('-mat', filename);
else
  EEG = filename;
end;

header.Fs          = EEG.srate;
header.nChans      = EEG.nbchan;
header.nSamples    = EEG.pnts;
header.nSamplesPre = -EEG.xmin*EEG.srate;
header.nTrials     = EEG.trials;
try
  header.label       = { EEG.chanlocs.labels }';
catch
  warning('creating default channel names');
  for i=1:header.nChans
    header.label{i} = sprintf('chan%03d', i);
  end
end
ind = 1;
for i = 1:length( EEG.chanlocs )
    if ~isempty(EEG.chanlocs(i).X)
        header.elec.label{ind, 1} = EEG.chanlocs(i).labels;
        % this channel has a position
        header.elec.pnt(ind,1) = EEG.chanlocs(i).X;
        header.elec.pnt(ind,2) = EEG.chanlocs(i).Y;
        header.elec.pnt(ind,3) = EEG.chanlocs(i).Z;
        ind = ind+1;
    end;
end;

% remove data
% -----------
%if isfield(EEG, 'datfile')
%    if ~isempty(EEG.datfile)
%        EEG.data = EEG.datfile;
%    end;
%else
%    EEG.data = 'in set file';
%end;
EEG.icaact = [];

header.orig = EEG;
