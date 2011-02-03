function [hdr,be] = read_hdr_raw(fname)
% Read a NIFTI-1 hdr file
% FORMAT [hdr,be] = read_hdr_raw(fname)
% fname - filename of image
% hdr   - a structure containing hdr info
% be    - whether big-endian or not
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id$


hdr = [];
be  = [];
ok  = true;

% Name of header file
[pth,nam,ext] = fileparts(fname);
switch ext
case {'.img','.hdr'}
    hname = fullfile(pth,[nam '.hdr']);
case {'.nii'}
    hname = fullfile(pth,[nam '.nii']);
otherwise
    hname = fullfile(pth,[nam '.hdr']);
end;

% Open it if possible
fp  = fopen(hname,'r','native');
if fp==-1
    hdr = [];
    return;
end;

% Sort out endian issues of header
[unused,unused,mach] = fopen(fp);
if strcmp(mach,'ieee-be') || strcmp(mach,'ieee-be.l64')
    be = true;
else
    be = false;
end;
fseek(fp,0,'bof');
fseek(fp,40,'bof');
nd = fread(fp,1,'int16')';
if isempty(nd),
    fclose(fp);
    hdr = [];
    return;
elseif nd<1 || nd>7
    be = ~be;
    fclose(fp);
    if be, mach = 'ieee-be';
    else   mach = 'ieee-le';
    end;
    fp = fopen(hname,'r',mach);
    if fp==-1
        hdr = [];
        return;
    end;
end;

% Is it NIFTI or not
fseek(fp,0,'bof');
fseek(fp,344,'bof');
mgc = deblank(char(fread(fp,4,'uint8')'));
switch mgc
case {'ni1','n+1'}
    org = niftistruc;
otherwise
    org = mayostruc;
end;
fseek(fp,0,'bof');
% Read the fields
for i=1:length(org)
    tmp = fread(fp,org(i).len,['*' org(i).dtype.prec])';
    if length(tmp) ~= org(i).len
disp([length(tmp) org(i).len]);
        tmp = org(i).def;
        ok  = false;
    end;
    tmp = feval(org(i).dtype.conv,tmp);
    hdr.(org(i).label) = tmp;
end;

fclose(fp);
if ~ok,
     fprintf('There was a problem reading the header of\n');
     fprintf('"%s".\n', fname);
     fprintf('It may be corrupted in some way.');
end;
return;

