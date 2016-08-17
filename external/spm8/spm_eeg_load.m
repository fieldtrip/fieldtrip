function D = spm_eeg_load(P)
% Load an M/EEG file in SPM format
% FORMAT D = spm_eeg_load(P)
%
% P         - filename of M/EEG file
% D         - MEEG object 
%__________________________________________________________________________
% 
% spm_eeg_load loads an M/EEG file using the SPM8 format. Importantly, the
% data array is memory-mapped and the struct is converted to MEEG object.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_load.m 3879 2010-05-07 20:03:34Z vladimir $

% bypass if the input is already an MEEG object
%--------------------------------------------------------------------------
if nargin && isa(P, 'meeg')
    D = P;
    return;
end

% get filename
%--------------------------------------------------------------------------
try
    P = deblank(P);
catch
    [P, sts] = spm_select(1, 'mat', 'Select SPM M/EEG file');
    if ~sts, D = []; return; end
end

Ppath = fileparts(P);
if isempty(Ppath)
    Ppath = pwd;
end

% load MAT file
%--------------------------------------------------------------------------
try
    load(P);
catch    
    error('Trouble reading file %s', P);
end

% check whether there is a struct D
%--------------------------------------------------------------------------
if ~exist('D','var')
    error('%s doesn''t contain SPM M/EEG data', P);
end

% This is for the case when people save the object in a file
%--------------------------------------------------------------------------
if ~isa(D, 'struct')
    try
        D = struct(D);
    catch
        error('The file should contain an SPM M/EEG struct named D');
    end
end

% save path in structure and make sure 'type' is set
%--------------------------------------------------------------------------
D.path = Ppath;

try
    D.type;
catch
    D.type = 'other';
end

% return an MEEG object
%--------------------------------------------------------------------------
D = meeg(D);
