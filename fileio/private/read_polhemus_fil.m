function [fid, sens, label] = read_polhemus_fil(Fname_pol,skip)

% Reads Polhemus files:
%   either sensor file or headshape file or both
%
% FORMAT [fid, sens, label] = read_polhemus_fil(Fname_pol,skip)
% Input:
% Fname_pol - Polhemus ASCII file containing sensor locations (cm)
%             (headshape can also be considered here instead of sensors)
% skip      - first channels to skip
%
% Output:
% fid       - fiducial         locations (mm) in rows
% sens      - sensor/headshape locations (mm) in rows
% label - labels of the fiducials
%
% IMPORTANT: Note that Polhemus data files should be -ASCII files with
% extension .pol
% It is assumed that the .pol file contains the location (cm) of fiducials
% (sampled twice), possibly followed by some additional named points and 
% then unnamed location of the sensors.  In some instances the first
% few channel locations may pertain to reference channels; the skip 
% variable allows these to be skipped if necessary. The fiducial locations
% are flaged with the strings 'NZ','LE' and 'RE'; indicating the Nasion,
% left and right eare respectively.
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id$

% checks and assigments
%--------------------------------------------------------------------------
try, skip;    catch, skip    = 0; end

[pth,nam,ext] = fileparts(Fname_pol);
if ~strcmp(ext,'.pol')
    warndlg(sprintf('Wrong input file format\n'));
    return
end

% --- READ Polhemus Sensor + fiducial locations ---
%==========================================================================
try
    file = textread(Fname_pol,'%s');
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
label = {};
fid = [];
temp = 0;
nl   = 1;

while temp == 0
    if strcmp(file{nl},'NZ')
        NZ = [NZ ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
        nl = nl + 4;
    elseif strcmp(file{nl},'LE') || strcmp(file{nl},'OG')
        LE = [LE ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
        nl = nl + 4;
    elseif strcmp(file{nl},'RE') || strcmp(file{nl},'OD')
        RE = [RE ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
        nl = nl + 4;
    elseif isempty(str2num(file{nl})) % Add possible other named points
        label = [label, file(nl)];
        fid = [fid; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
        nl = nl + 4;
    else
        temp = 1;
    end
end

% convert from cm to mm
%--------------------------------------------------------------------------
NZ    = mean(NZ,1); LE = mean(LE,1); RE = mean(RE,1);
fid   = [NZ; LE; RE; fid]; 

label = [{'nas', 'lpa', 'rpa'}, label];

% read sensor locations or headshape locations
%--------------------------------------------------------------------------
sens  = [];
start = nl + skip*3;
for i = start:3:length(file)
    sens = [sens; str2num(file{i}) str2num(file{i+1}) str2num(file{i+2})];
end
