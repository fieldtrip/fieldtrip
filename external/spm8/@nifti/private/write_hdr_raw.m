function ok = write_hdr_raw(fname,hdr,be)
% Write a NIFTI-1 .hdr file.
% FORMAT ok = write_hdr_raw(fname,hdr,be)
% fname - filename of image
% hdr   - a structure containing hdr info
% be    - whether big-endian or not
% ok    - status (1=good, 0=bad)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id$


[pth,nam,ext] = fileparts(fname);
if isempty(pth), pth = pwd; end

if isfield(hdr,'magic')
    org = niftistruc;
    switch deblank(hdr.magic)
    case {'ni1'}
        hname = fullfile(pth,[nam '.hdr']);
    case {'n+1'}
        hname = fullfile(pth,[nam '.nii']);
    otherwise
        error('Bad header.');
    end;
else
    org   = mayostruc;
    hname = fullfile(pth,[nam '.hdr']);
end;

if nargin >=3
    if be, mach = 'ieee-be';
    else   mach = 'ieee-le';
    end;
else       mach = 'native';
end;

ok = true;
if spm_existfile(hname),
    fp = fopen(hname,'r+',mach);
else
    fp = fopen(hname,'w+',mach);
end
if fp == -1,
    ok = false;
    return;
end

for i=1:length(org)
    if isfield(hdr,org(i).label),
        dat = hdr.(org(i).label);
        if length(dat) ~= org(i).len,
            if length(dat)< org(i).len,
                dat = [dat(:) ; zeros(org(i).len-length(dat),1)];
            else
                dat = dat(1:org(i).len);
            end;
        end;
    else
        dat = org(i).def;
    end;
    % fprintf('%s=\n',org(i).label)
    % disp(dat)
    len = fwrite(fp,dat,org(i).dtype.prec);
    if len ~= org(i).len,
        ok = false;
    end;
end;
fclose(fp);
if ~ok,
     fprintf('There was a problem writing to the header of\n');
     fprintf('"%s"\n', fname);
end;
return;

