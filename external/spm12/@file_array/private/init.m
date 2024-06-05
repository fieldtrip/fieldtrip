function init(fname, nbytes, opts)
% Initialise binary file on disk
% FORMAT init(fname, nbytes[, opts])
% fname   - filename
% nbytes  - data size {bytes}
% opts    - optional structure with fields:
%   .offset   - file offset {bytes} [default: 0]
%   .wipe     - overwrite existing values with 0 [default: false]
%   .truncate - truncate file if larger than requested size [default: true]
%
% This function is normally called by file_array/initialise
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
error('init.c not compiled - see Makefile');
