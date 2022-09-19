function [vox_idx, t_idx] = nmt_peaksearch(cfg)
% [vox_idx, t_idx] = nmt_peaksearch(cfg)
% cfg.time = single time point, or time range, or 'current'; if unspecified, search over time at specified voxel
% cfg.vox = find peak time at specified voxel, or 'current'; if unspecified, search over voxels at specified time
% cfg.peaktype = 'mag' (max magnitude, default) or 'max' or 'min'
% cfg.searchradius = minimum and maximum distance to search for peak

global st

if(~isfield(cfg,'peaktype'))
    cfg.peaktype = 'mag';
end

if(~isfield('cfg','axsel'))
    cfg.axsel = 1;
end

currvox = st.nmt.cfg.vox_idx; % currently selected voxel
currtime = st.nmt.cfg.time_idx; % currently selected time
currfreq = st.nmt.cfg.freq_idx; % currently selected freq

if(isfield(cfg,'vox'))
    switch(class(cfg.vox))
        case 'char'
            switch(cfg.vox)
                case 'current'
                    cfg.vox = st.nmt.vox_idx; % current voxel
                otherwise
                    error('either ''current'' or voxel coords/index expected.')
            end
        case 'double'
            if(length(cfg.vox)==3) % coordinate rather than index
                cfg.vox = dsearchn(st.nmt.pos,cfg.vox);
            end
        otherwise
            error('something unexpected happened here...')
    end
end
                               
if(isfield(cfg,'time'))
    switch(class(cfg.time))
        case 'char'
            cfg.time = st.nmt.cfg.time_idx; % current time point
        case 'double'
            if(length(cfg.time)==3)  % search for index corresponding to requested time
                cfg.time = (dsearchn(st.nmt.time',cfg.time'))';
            end
        otherwise
            error('something unexpected happened here...')
    end
end

%% define search area
if(isfield(cfg,'searchradius'))
    distances = nmt_rownorm(nmt_coord_diff(st.nmt.pos,st.nmt.pos(currvox,:)));
    excludeindices = find(distances<cfg.searchradius(1) | distances>cfg.searchradius(2));
else
    excludeindices = [];
end

fun = st.nmt.fun{cfg.axsel}(:,:,st.nmt.cfg.freq_idx(1));
fun(excludeindices,:,:,:) = NaN; % NaN out values outside search range

%%
if(~isfield(cfg,'time') && ~isfield(cfg,'vox'))
    switch(cfg.peaktype)
        case 'max'
            [dum,peakind] = max(fun(:));
            [peakvoxind,peaktimeind] = ind2sub(size(fun),peakind);
        case 'mag'
            [dum,peakind] = max(abs(fun(:)));
            [peakvoxind,peaktimeind] = ind2sub(size(fun),peakind);
        case 'min'
            [dum,peakind] = min(fun(:));
            [peakvoxind,peaktimeind] = ind2sub(size(fun),peakind);
    end
    t_idx = peaktimeind;
    vox_idx = peakvoxind;
elseif(~isfield(cfg,'time') && isfield(cfg,'vox'))
    switch(cfg.peaktype)
        case 'max'
            [peakval,t_idx_peak] = max(fun(cfg.vox,:));
        case 'mag'
            [peakval,t_idx_peak] = max(abs(fun(cfg.vox,:)));
        case 'min'
            [peakval,t_idx_peak] = min(fun(cfg.vox,:));
    end
    t_idx = t_idx_peak;
    vox_idx = cfg.vox; % searchign over time, so voxel doesn't change
elseif(isfield(cfg,'time') && ~isfield(cfg,'vox'))
    t_idx = cfg.time;
    if(t_idx(2)~=t_idx(1))
        funval = nmt_ts_intervalpower(fun(:,t_idx(1):t_idx(2)));
    else
        funval = fun(:,t_idx(1));
    end
    switch(cfg.peaktype)
        case 'max'
            [peakval,vox_idx_peak] = max(funval);
        case 'mag'
            [peakval,vox_idx_peak] = max(abs(funval));
        case 'min'
            [peakval,vox_idx_peak] = min(funval);
    end
    vox_idx = vox_idx_peak;
    t_idx = [t_idx(1) t_idx(end)]; %  searching over space, so time doesn't change
end

