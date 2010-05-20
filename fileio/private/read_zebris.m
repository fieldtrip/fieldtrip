function [fid, sens, label, sens_label] = read_zebris(Fname_zeb,skip)

% Reads Zebris files:
%   fiducials locations, and
%   either sensor file or headshape file or both
%
% FORMAT [fid, sens, label] = read_zebris(Fname_zeb,skip)
% Input:
% Fname_zeb  - Zebris ASCII file containing sensor locations (mm)
%             (headshape can also be considered here instead of sensors)
% skip       - first channels to skip
%
% Output:
% fid        - fiducial         locations (mm) in rows
% sens       - sensor/headshape locations (mm) in rows
% label      - labels of the fiducials
% sens_label - labels of the surface points, electrodes + headshape
%
% IMPORTANT: Note that Zebris data files should be -ASCII files with
% extension .sfp
% It is assumed that the .sfp file contains the location (mm) of fiducials
% (possibly twice), possibly followed by some additional named points for
% the electrodes, and then so more named location starting with 'sfl' for
% headshape locations.
% In some instances the first few channel locations may pertain to
% reference channels; the skip variable allows these to be skipped if
% necessary.
% The fiducial locations are flaged with the strings 'fidt9','fidnz' and
% 'fidt10'; indicating the leaft ear, nasion, and right ear, respectively.
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id$


% checks and assigments
%--------------------------------------------------------------------------
try, skip;    catch, skip    = 0; end

[pth,nam,ext] = fileparts(Fname_zeb);
if ~strcmp(ext,'.sfp')
  warndlg(sprintf('Wrong input file format\n'));
  return
end


% --- READ Zebris Sensor + fiducial locations ---
%==========================================================================
try
  file = textread(Fname_zeb,'%s');
catch
  file = textread(fullfile(pwd,[nam ext]),'%s');
end
% remove zeros at the end
bool = 0;
while bool == 0
  if strcmp(file{end},'0')
    file = file(1:end-1);
  else
    bool = 1;
  end
end

% read fiducials
%--------------------------------------------------------------------------
NZ   = [];
LE   = [];
RE   = [];
temp = 0;
nl   = 1;

while temp == 0
  if strcmp(file{nl},'fidnz')
    NZ = [NZ ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
    nl = nl + 4;
  elseif strcmp(file{nl},'fidt9')
    LE = [LE ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
    nl = nl + 4;
  elseif strcmp(file{nl},'fidt10')
    RE = [RE ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
    nl = nl + 4;
  else
    temp = 1;
  end
end

% convert from cm to mm
%--------------------------------------------------------------------------
NZ    = mean(NZ,1); LE = mean(LE,1); RE = mean(RE,1);
fid   = [NZ; LE; RE];

label = [{'nas', 'lpa', 'rpa'}];

% read sensor locations or headshape locations
%--------------------------------------------------------------------------
sens  = [];
sens_label = [];

start = nl + skip*3;
for ii = start:4:length(file)
  sens_label = [sens_label, file(ii)];
  sens = [sens; ...
    str2num(file{ii+1}) str2num(file{ii+2}) str2num(file{ii+3})];
end
