function nmt_spm8_plot(cfg)
global st

if(~exist('cfg','var'))
    cfg = st.nmt.cfg;
end

% if there is no time dimension, set cfg.time to index first and only column
if(~isfield(cfg,'time'))
    cfg.time = [1 1];
end

if(cfg.time(1) == cfg.time(2)) % single time point selected
    fun = st.nmt.fun(:,cfg.time(1));
    
    scalemax = max(abs(st.nmt.fun(:)));
    scalemin = -max(abs(st.nmt.fun(:)));
    
    % set colorscale: anatomical (first 64) + functional (second 64)
    % grayscale for MRI; jet for painted activation
    set(st.fig,'Colormap',[gray(64);jet(64)]);
else % time interval selected
    fun = nmt_ts_intervalpower(st.nmt.fun(:,cfg.time(1):cfg.time(2)));
    
    scalemax = max(fun(:));
    scalemin = min(fun(:));
    
    % set colorscale: anatomical (first 64) + functional (second 64)
    % grayscale for MRI; jet for painted activation
    set(st.fig,'Colormap',[gray(64);hot(64)]);
end




%% figure out voxelsize
% in case of nonstandard nifti orientation
voxelspacing = abs(st.nmt.pos(2,:) - st.nmt.pos(1,:));
voxeldir = find(voxelspacing); % determine which dimension contains a non-zero value
if(length(voxeldir) > 1)
    error('something is strange... is your voxel grid uniform and aligned to your coordinate system?')
end
voxelsize = st.nmt.pos(2,voxeldir) - st.nmt.pos(1,voxeldir);

%% pos in MRI space
blob2mri_tfm = [voxelsize          0          0  min(st.nmt.pos(:,1))-voxelsize
                        0  voxelsize          0  min(st.nmt.pos(:,2))-voxelsize
                        0          0  voxelsize  min(st.nmt.pos(:,3))-voxelsize
                        0          0          0                               1 ];
%blob = ft_transform_geometry(inv(blob2mri_tfm),functional);
blob.pos = nmt_transform_coord(inv(blob2mri_tfm),st.nmt.pos);

%blob.transform = blob2mri_tfm;




if(~isfield(st.vols{1},'blobs'))
    % create new "hot" opaque blob, with colorbar
    spm_orthviews('addblobs',1,blob.pos',fun',blob2mri_tfm);
else
    % blob already exists, just update values
    st.vols{1}.blobs{1}.vol = reshape(fun',size(st.vols{1}.blobs{1}.vol));
end
st.vols{1}.blobs{1}.max = scalemax;
st.vols{1}.blobs{1}.min = scalemin;

% ft_image
spm_orthviews('redraw'); % blob doesn't always show until redraw is forced


function fun = nmt_ts_intervalpower(fun);
% when time interval is selected, calculate desired representation of interval power
switch(2)
    case 1 % mean power
        fun = mean(fun.^2,2);
    case 2 % RMS power
        fun = sqrt(mean(fun.^2,2));
    case 3  % simple average
        fun = mean(fun,2);
end