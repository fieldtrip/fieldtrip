function t = structn(obj)
% Convert a NIFTI-1 object into a form of struct
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: structn.m 1143 2008-02-07 19:33:33Z spm $


if numel(obj)~=1,
    error('Too many elements to convert');
end;
fn = fieldnames(obj);
for i=1:length(fn)
    tmp = subsref(obj,struct('type','.','subs',fn{i}));
    if ~isempty(tmp)
        t.(fn{i}) = tmp;
    end;
end;
return;
