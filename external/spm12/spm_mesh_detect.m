function s = spm_mesh_detect(F)
% True for valid representation of a mesh
% FORMAT s = spm_mesh_detect(F)
% F        - variable to query: filename, vol structure, patch structure
% s        - true if F corresponds to a mesh, and false otherwise
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


s = false;

if iscellstr(F) || ischar(F)
    F = cellstr(F);
    [p,n,e] = spm_fileparts(F{1});
    if strcmpi(e,'.gii')
        s = true;
    end
elseif isstruct(F)
    if isfield(F,'vertices') && isfield(F,'faces')
        s = true;
    elseif isfield(F,'cdata')
        s = true;
    elseif isfield(F,'private') && isa(F(1).private,'gifti')
        s = true;
    end
elseif isa(F,'gifti')
    s = true;
end
