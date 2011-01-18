function extras = write_extras(fname,extras)
% Write extra bits of information
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: write_extras.m 1143 2008-02-07 19:33:33Z spm $


[pth,nam,ext] = fileparts(fname);
switch ext
case {'.hdr','.img','.nii'}
    mname = fullfile(pth,[nam '.mat']);
case {'.HDR','.IMG','.NII'}
    mname = fullfile(pth,[nam '.MAT']);
otherwise
    mname = fullfile(pth,[nam '.mat']);
end
if isstruct(extras) && ~isempty(fieldnames(extras)),
    savefields(mname,extras);
end;

function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
for i_=1:length(fn),
    eval([fn{i_} '= p.' fn{i_} ';']);
end;
if str2num(version('-release'))>=14,
    fn = {'-V6',fn{:}};
end;
if numel(fn)>0,
    save(fnam,fn{:});
end;
return;

