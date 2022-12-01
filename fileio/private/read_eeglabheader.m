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

function header = read_eeglabheader(filename)

if nargin < 1
  help read_eeglabheader
  return
end

if ~isstruct(filename)
  s = whos('-file', filename);
  if any(strcmp({s.name}', 'EEG'))
    load('-mat', filename, 'EEG');
  else % EEGLAB > 2021.0 saves content of EEG as default
    EEG = load('-mat', filename);
  end
else
  EEG = filename;
end

header.Fs          =  EEG.srate;
header.nChans      =  EEG.nbchan;
header.nSamples    =  EEG.pnts;
header.nSamplesPre = -EEG.xmin*EEG.srate;
header.nTrials     =  EEG.trials;
try
  header.label     = { EEG.chanlocs.labels }';
catch
  ft_warning('creating default channel names');
  for i=1:header.nChans
    header.label{i} = sprintf('chan%03d', i);
  end
end
ind = 1;
for i = 1:length( EEG.chanlocs )
  if isfield(EEG.chanlocs(i), 'X') && ~isempty(EEG.chanlocs(i).X)
    header.elec.label{ind, 1} = EEG.chanlocs(i).labels;
    % this channel has a position
    header.elec.elecpos(ind,1) = EEG.chanlocs(i).X;
    header.elec.elecpos(ind,2) = EEG.chanlocs(i).Y;
    header.elec.elecpos(ind,3) = EEG.chanlocs(i).Z;
    ind = ind+1;
  end
end

if ind>1 && isfield(EEG, 'chaninfo') && isfield(EEG.chaninfo, 'nosedir') && ~isempty(EEG.chaninfo.nosedir)
  switch EEG.chaninfo.nosedir
    case '+X'
      header.elec.coordsys = 'ctf';
    case '-X'
      header.elec.coordsys = 'pls';
    case '+Y'
      header.elec.coordsys = 'ras';
    case '-Y'
      header.elec.coordsys = 'lpi';
    otherwise
  end
end

% remove data
% -----------
%if isfield(EEG, 'datfile')
%    if ~isempty(EEG.datfile)
%        EEG.data = EEG.datfile;
%    end
%else
%    EEG.data = 'in set file';
%end;
EEG.icaact = [];

header.orig = EEG;
