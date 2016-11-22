function nmt_spm_plot(cfg,axsel)
% nmt_spm_plot(cfg)
% plots desired activation on SPM8 viewer
% designed for use via nmt_sourceplot, but possibly usable independently

global st

if(~exist('cfg','var') || isempty(cfg))
    cfg = st.nmt.cfg;
end

if(~exist('axsel','var'))
    axsel = 1;
end

if(~isfield(cfg,'colormap'))
    cfg.colormap = jet(64);
end


% if there is no time dimension, set cfg.time_idx to index first and only column
if(~isfield(cfg,'time_idx'))
    cfg.time_idx = [1 1];
end

% if there is no frequency dimension, set cfg.time_idx to index first and only column
if(~isfield(cfg,'freq_idx'))
    cfg.freq_idx = [1 1];
end


if(cfg.time_idx(1) == cfg.time_idx(2) && cfg.freq_idx(1) == cfg.freq_idx(2)) % single time point selected
    fun = st.nmt.fun{axsel}(:,cfg.time_idx(1),cfg.freq_idx(1));
    % set colorscale: anatomical (first 64) + functional (second 64)
    % grayscale for MRI; jet for painted activation
    set(st.fig,'Colormap',[gray(64);cfg.colormap]);
    
    scalemax = max(abs(st.nmt.fun{axsel}(:)));
    scalemin = -max(abs(st.nmt.fun{axsel}(:)));
else % time interval selected
    switch(cfg.plottype)
        case 'tf'
            fun = nmt_ts_intervalpower(st.nmt.fun{axsel}(:,cfg.time_idx(1):cfg.time_idx(2),cfg.freq_idx(1):cfg.freq_idx(2)),'mean');
            % set colorscale: anatomical (first 64) + functional (second 64)
            % grayscale for MRI; jet for painted activation
            set(st.fig,'Colormap',[gray(64);cfg.colormap]);
        otherwise
            fun = nmt_ts_intervalpower(st.nmt.fun{axsel}(:,cfg.time_idx(1):cfg.time_idx(2),cfg.freq_idx(1):cfg.freq_idx(2)),'rms');
            % set colorscale: anatomical (first 64) + functional (second 64)
            % grayscale for MRI; jet for painted activation
            set(st.fig,'Colormap',[gray(64);hot(64)]);
    end
    
    scalemax = max(fun(:));
    scalemin = min(fun(:));
    
end


%% apply mask
msk = st.nmt.msk;
msk(msk==0) = NaN;
maskedfun = (msk(:,cfg.time_idx(1),cfg.freq_idx(1)).*fun)';

%% figure out voxelsize
% voxelsize assumed to be most common non-zero difference value
diffcoords = abs(diff(st.nmt.pos));
diffcoords(diffcoords==0)=NaN; % replace zeroes with NaN
voxelsize = mode(diffcoords(:));

%% pos in MRI space
blob2mri_tfm = [voxelsize          0          0  min(st.nmt.pos(:,1))-voxelsize
                        0  voxelsize          0  min(st.nmt.pos(:,2))-voxelsize
                        0          0  voxelsize  min(st.nmt.pos(:,3))-voxelsize
                        0          0          0                               1 ];
blob.pos = nmt_transform_coord(inv(blob2mri_tfm),st.nmt.pos);


if(~isfield(st.vols{1},'blobs'))
    % create new "hot" opaque blob, with colorbar
    spm_orthviews('addblobs',1,blob.pos',maskedfun,blob2mri_tfm);
else
    % blob already exists, just update values
    st.vols{1}.blobs{1}.vol = reshape(maskedfun,size(st.vols{1}.blobs{1}.vol));
end

if(isfield(cfg,'zlim'))
    st.vols{1}.blobs{1}.max = cfg.zlim(2);
    st.vols{1}.blobs{1}.min = cfg.zlim(1);
else
    st.vols{1}.blobs{1}.max = scalemax;
    st.vols{1}.blobs{1}.min = scalemin;
end

% ft_image
spm_orthviews('redraw'); % blob doesn't always show until redraw is forced

%% refresh current voxel and functional quantity display
posmrimm = spm_orthviews('pos')';
% set(st.nmt.gui.megp,'String',sprintf('%.1f %.1f %.1f',posmeg));

blobidx = nmt_transform_coord(inv(st.vols{1}.blobs{1}.mat),posmrimm);
blobidx = round(blobidx);
blobdim = size(st.vols{1}.blobs{1}.vol);
if(all(blobidx <= blobdim & blobidx > 0))
    st.nmt.cfg.vox_idx  = sub2ind(size(st.vols{1}.blobs{1}.vol),blobidx(1),blobidx(2),blobidx(3));
    funval = fun(st.nmt.cfg.vox_idx);
else % if out of blob's bounds
    funval = NaN;
    st.nmt.cfg.vox_idx = NaN;
end

set(st.nmt.gui.beamin,'String',sprintf('%+g',funval));


%nmt_update_panel(axsel)


