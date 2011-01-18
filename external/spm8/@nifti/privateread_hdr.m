function vol = read_hdr(fname)
% Get a variety of information from a NIFTI-1 header.
% FORMAT vol = read_hdr(fname)
% fname - filename of image
% vol   - various bits of information
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: read_hdr.m 1143 2008-02-07 19:33:33Z spm $


persistent d
if isempty(d), d = getdict; end;

[pth,nam,ext,num] = spm_fileparts(fname);
switch ext
case '.hdr'
    hname = fullfile(pth,[nam '.hdr']);
case '.HDR'
    hname = fullfile(pth,[nam '.HDR']);
case '.img'
    hname = fullfile(pth,[nam '.hdr']);
case '.IMG'
    hname = fullfile(pth,[nam '.HDR']);
case '.nii'
    hname = fullfile(pth,[nam '.nii']);
case '.NII'
    hname = fullfile(pth,[nam '.NII']);
otherwise
    hname = fullfile(pth,[nam ext]);
end
[hdr,be] = read_hdr_raw(hname);
if isempty(hdr)
    error(['Error reading header file "' hname '"']);
end;

if ~isfield(hdr,'magic'),
    % A patch for some possible SPM2 datatypes
    switch hdr.datatype,
    case 130, hdr.datatype = 256; %  int8
    case 132, hdr.datatype = 512; % uint16
    case 136, hdr.datatype = 768; % uint32
    end;
end;

dt = [];
for i=1:length(d.dtype)
    if hdr.datatype == d.dtype(i).code
        dt = d.dtype(i);
        break;
    end;
end;
if isempty(dt)
    error(['Unrecognised datatype (' num2str(double(hdr.datatype)) ') for "' fname '.'] );
end
if isfield(hdr,'magic')
    switch deblank(hdr.magic)
    case {'n+1'}
        iname = hname;
        if hdr.vox_offset < hdr.sizeof_hdr
            error(['Bad vox_offset (' num2str(double(hdr.vox_offset)) ') for "' fname '.'] );
        end
    case {'ni1'}
        if strcmp(ext,lower(ext)),
            iname = fullfile(pth,[nam '.img']);
        else
            iname = fullfile(pth,[nam '.IMG']);
        end;
    otherwise
        error(['Bad magic (' hdr.magic ') for "' fname '.'] );
    end
else
    if strcmp(ext,lower(ext)),
        iname = fullfile(pth,[nam '.img']);
    else
        iname = fullfile(pth,[nam '.IMG']);
    end;
end
if rem(double(hdr.vox_offset),dt.size)
   error(['Bad alignment of voxels (' num2str(double(hdr.vox_offset)) '/' num2str(double(dt.size)) ') for "' fname '.'] );
end;

vol      = struct('hdr',hdr,'be',be,'hname',hname,'iname',iname);
return
