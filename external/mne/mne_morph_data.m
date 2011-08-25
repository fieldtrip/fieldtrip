function [stcs] = mne_morph_data(from,to,stcs,grade)
% MNE_MORPH_DATA   Returns data morphed to a new subject.
%
%   SYNTAX
%       [STCS] = MNE_MORPH_DATA(FROM, TO, STCS, GRADE)
%
%   from : name of origin subject
%   to : name of destination subject
%   stcs : stc data to morph
%   grade : (optional) resolution of the icosahedral mesh (typically 5)
%
% Note : The functions requires to set MNE_ROOT and SUBJECTS_DIR variables.
%
% Example:
%  from = 'sample';
%  to = 'fsaverage';
%  stcs_morph = mne_morph_data(from,to,stcs,5);
%
%

%
%   Author : Alexandre Gramfort, MGH Martinos Center
%            Matti Hamalainen
%   License : BSD 3-clause
%

if nargin < 5
    grade = 5;
end

sphere_from = [getenv('SUBJECTS_DIR'),'/',from,'/surf/lh.sphere.reg'];
[tmp, tris{1}] = mne_read_surface(sphere_from);
sphere_from = [getenv('SUBJECTS_DIR'),'/',from,'/surf/rh.sphere.reg'];
[tmp, tris{2}] = mne_read_surface(sphere_from);

[map{1}, map{2}] = mne_read_morph_map(from,to);

for hemi = 1:2
    e = mne_mesh_edges(tris{hemi});
    e = e==2;
    n_vertices = length(e);
    e = e + speye(n_vertices, n_vertices);
    idx_use = stcs{hemi}.vertices + 1;
    n_iter = 100; % nb of smoothing iterations
    for k = 1:n_iter
        data1 = e(:,idx_use) * ones(length(idx_use),1);
        stcs{hemi}.data  = e(:,idx_use)*stcs{hemi}.data;
        idx_use = find(data1);
        fprintf(1,'%d/%d ',k,length(idx_use));
        if ( k == n_iter ) || ( length(idx_use) >= n_vertices )
            stcs{hemi}.data(idx_use,:) = bsxfun(@rdivide, stcs{hemi}.data(idx_use,:), data1(idx_use));
            break;
        else
            stcs{hemi}.data = bsxfun(@rdivide, stcs{hemi}.data(idx_use,:), data1(idx_use));
        end
    end
    fprintf(1,'\n');
    stcs{hemi}.data = map{hemi}*stcs{hemi}.data;
end

ico_file_name = [getenv('MNE_ROOT'),'/share/mne/icos.fif'];

s = mne_read_bem_surfaces(ico_file_name);

for k = 1:length(s)
    if (s(k).id == 9000 + grade)
        ico = s(k);
        break;
    end
end

sphere = [getenv('SUBJECTS_DIR'),'/',to,'/surf/lh.sphere.reg'];
lhs = mne_read_surface(sphere);
sphere = [getenv('SUBJECTS_DIR'),'/',to,'/surf/rh.sphere.reg'];
rhs = mne_read_surface(sphere);

nearest{1} = zeros(ico.np, 1);
nearest{2} = zeros(ico.np, 1);

lhs =  bsxfun(@rdivide, lhs, sqrt(sum(lhs.^2, 2)));
rhs =  bsxfun(@rdivide, rhs, sqrt(sum(rhs.^2, 2)));

for k = 1:size(ico.rr,1)
    dots = ico.rr(k,:)*lhs';
    [tmp, nearest{1}(k)] = max(dots);
    dots = ico.rr(k,:)*rhs';
    [tmp, nearest{2}(k)] = max(dots);
end

stcs{1}.data = stcs{1}.data(nearest{1},:);
stcs{2}.data = stcs{2}.data(nearest{2},:);
stcs{1}.vertices = nearest{1}-1;
stcs{2}.vertices = nearest{2}-1;

fprintf(1, '\n');
