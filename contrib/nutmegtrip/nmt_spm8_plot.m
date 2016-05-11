function nmt_spm8_plot(cfg)
% nmt_spm8_plot(cfg)
% plots desired activation on SPM8 viewer
% designed for use via nmt_sourceplot_spm8, but possibly usable independently

global st

if(~exist('cfg','var'))
    cfg = st.nmt.cfg;
end

if(~isfield('cfg','axsel'))
    cfg.axsel = 1;
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
    fun = st.nmt.fun{cfg.axsel}(:,cfg.time_idx(1),cfg.freq_idx(1));
    % set colorscale: anatomical (first 64) + functional (second 64)
    % grayscale for MRI; jet for painted activation
    set(st.fig,'Colormap',[gray(64);jet(64)]);
    
    scalemax = max(abs(st.nmt.fun{cfg.axsel}(:)));
    scalemin = -max(abs(st.nmt.fun{cfg.axsel}(:)));
else % time interval selected
    switch(cfg.plottype)
        case 'tf'
            fun = nmt_ts_intervalpower(st.nmt.fun{cfg.axsel}(:,cfg.time_idx(1):cfg.time_idx(2),cfg.freq_idx(1):cfg.freq_idx(2)),'mean');
            % set colorscale: anatomical (first 64) + functional (second 64)
            % grayscale for MRI; jet for painted activation
            set(st.fig,'Colormap',[gray(64);jet(64)]);
        otherwise
            fun = nmt_ts_intervalpower(st.nmt.fun{cfg.axsel}(:,cfg.time_idx(1):cfg.time_idx(2),cfg.freq_idx(1):cfg.freq_idx(2)),'rms');
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



%% update time series, if applicable
if(isfield(st.nmt,'time')) %& ~isfield(st.nmt,'freq'))
    set(st.nmt.gui.timeguih,'Visible','On'); % ensure plot is visible
    switch(st.nmt.cfg.plottype)
        case 'tf'
            set(st.nmt.gui.freqguih,'Visible','On'); % ensure plot is visible
            nmt_tfplot(st.nmt.gui.ax_ts(cfg.axsel),st.nmt.time,st.nmt.freq,squeeze(st.nmt.fun{cfg.axsel}(st.nmt.cfg.vox_idx,:,:)),@nmt_repos_start);
        case 'ts'
            if(isfinite(st.nmt.cfg.vox_idx))
                plot(st.nmt.gui.ax_ts(cfg.axsel),st.nmt.time,squeeze(st.nmt.fun{cfg.axsel}(st.nmt.cfg.vox_idx,:,:)));
                grid(st.nmt.gui.ax_ts(cfg.axsel),'on');
            else
                plot(st.nmt.gui.ax_ts(cfg.axsel),st.nmt.time,nan(length(st.nmt.time),1));
            end
            
            % set y-axis range based on whole volume range (single time point), or
            % selected timeslice range
            if(st.nmt.cfg.time_idx(1)==st.nmt.cfg.time_idx(2))
                ymin = min(st.nmt.fun{cfg.axsel}(:));
                ymax = max(st.nmt.fun{cfg.axsel}(:));
            else
                ymin =  min(min(st.nmt.fun{cfg.axsel}(:,st.nmt.cfg.time_idx(1):st.nmt.cfg.time_idx(2))));
                ymax =  max(max(st.nmt.fun{cfg.axsel}(:,st.nmt.cfg.time_idx(1):st.nmt.cfg.time_idx(2))));
            end
            set(st.nmt.gui.ax_ts(cfg.axsel),'YLim',[ymin ymax]);
            
    end
    set(st.nmt.gui.ax_ts(cfg.axsel),'XLim',st.nmt.time([1 end]));
    xlabel(st.nmt.gui.ax_ts(cfg.axsel),['Time (s)']);
    
    
    %% plot vertical line indicating selected time point
    axisvalues = axis(st.nmt.gui.ax_ts(cfg.axsel));
    line_max=axisvalues(4);
    line_min=axisvalues(3);
    ts_xhair_color = 'red';
    ts_xhair_width = 1;
    
    switch(st.nmt.cfg.plottype)
        case 'ts'
            ylim=get(st.nmt.gui.ax_ts(cfg.axsel),'YLim');
        case 'tf'
            ylim = [st.nmt.freq(st.nmt.cfg.freq_idx(1),1) st.nmt.freq(st.nmt.cfg.freq_idx(2),2)]
    end
    axes(st.nmt.gui.ax_ts(cfg.axsel));
    h=patch([st.nmt.time(st.nmt.cfg.time_idx(1)) st.nmt.time(st.nmt.cfg.time_idx(2)) st.nmt.time(st.nmt.cfg.time_idx(2)) st.nmt.time(st.nmt.cfg.time_idx(1))]',[ylim(1) ylim(1) ylim(2) ylim(2)]',[400 400 400 400],[1 0.4 0.4],'EdgeColor','red');

    %% update GUI textboxes
    set(st.nmt.gui.t1,'String',num2str(st.nmt.time(st.nmt.cfg.time_idx(1))));
    set(st.nmt.gui.t2,'String',num2str(st.nmt.time(st.nmt.cfg.time_idx(2))));

    if(1)
        set(st.nmt.gui.f1,'Value',st.nmt.cfg.freq_idx(1));
        set(st.nmt.gui.f2,'Value',st.nmt.cfg.freq_idx(2));
    else
        set(st.nmt.gui.f1,'String',num2str(st.nmt.freq(st.nmt.cfg.freq_idx(1),1)));
        set(st.nmt.gui.f2,'String',num2str(st.nmt.freq(st.nmt.cfg.freq_idx(2),2)));
    end
        
    set(st.nmt.gui.ax_ts(cfg.axsel),'ButtonDownFcn',@nmt_repos_start);
    
    %% optionally add topoplot
    switch(st.nmt.cfg.topoplot)
        case 'timelock'
            cfgplot.xlim = st.nmt.time(st.nmt.cfg.time_idx);
            cfgplot.zlim = 'maxabs';
            nmt_addtopo(cfgplot,st.nmt.timelock);
        case {'spatialfilter','leadfield','leadfieldX','leadfieldY','leadfieldZ','leadfieldori'}
            % FIXME: this isn't so elegant :-)
            % currently steals some information from timelock structure
            topo.label = st.nmt.timelock.label;
            topo.grad = st.nmt.timelock.grad;
            topo.dimord = 'chan_time';
            topo.time = [0 1 2 3]; % not really time, but selection of magnitude or x/y/z components
            switch(cfg.topoplot)
                case 'spatialfilter'
                    topo.time = 0; % fake time (take vector norm)
                    cfgplot.xlim = [topo.time topo.time];
                    cfgplot.zlim = 'maxabs';
                    topo.avg = st.nmt.spatialfilter{st.nmt.cfg.vox_idx}';
                case {'leadfield','leadfieldX','leadfieldY','leadfieldZ'}
                    switch(cfg.topoplot)
                        case 'leadfield'
                            topo.time = 0; % fake time (take vector norm)
                            topo.avg = nut_rownorm(st.nmt.grid.leadfield{st.nmt.cfg.vox_idx});  % TODO: remove dependency on nut_rownorm (from NUTMEG)
                        case 'leadfieldX'
                            topo.time = 1;
                            topo.avg = st.nmt.grid.leadfield{st.nmt.cfg.vox_idx};
                        case 'leadfieldY'
                            topo.time = 2;
                            topo.avg = st.nmt.grid.leadfield{st.nmt.cfg.vox_idx};
                        case 'leadfieldZ'
                            topo.time = 3;
                            topo.avg = st.nmt.grid.leadfield{st.nmt.cfg.vox_idx};
                    end
                    cfgplot.xlim = [topo.time topo.time];
                    cfgplot.zlim = 'maxabs';
%                 case 'leadfieldX'
%                     topo.time = 1; % fake time (corresponds to x-orientation)
%                     cfgplot.xlim = [topo.time topo.time];
%                     cfgplot.zlim = 'maxabs';
%                     topo.avg = st.nmt.grid.leadfield{st.nmt.cfg.vox_idx};
%                 case 'leadfieldY'
%                     topo.time = 2; % fake time (corresponds to y-orientation)
%                     cfgplot.xlim = [topo.time topo.time];
%                     cfgplot.zlim = 'maxabs';
%                     topo.avg = st.nmt.grid.leadfield{st.nmt.cfg.vox_idx};
%                 case 'leadfieldZ'
%                     topo.time = 3; % fake time (corresponds to z-orientation)
%                     cfgplot.xlim = [topo.time topo.time];
%                     cfgplot.zlim = 'maxabs';
%                     topo.avg = nut_rownorm(st.nmt.grid.leadfield{st.nmt.cfg.vox_idx});
%                 case 'leadfieldori'
%                     error('TODO: not yet implemented');
            end
            
            nmt_addtopo(cfgplot,topo);
        case {'no',''}
            % do nothing
        otherwise
            error('requested cfg.topoplot parameter is not supported.');
    end
end




function nmt_repos_start(varargin)
global st
set(gcbf,'windowbuttonmotionfcn',@nmt_repos_move, 'windowbuttonupfcn',@nmt_repos_end);
nmt_timeselect('ts');
%_______________________________________________________________________
%_______________________________________________________________________
function nmt_repos_move(varargin)
nmt_timeselect('ts');
%_______________________________________________________________________
%_______________________________________________________________________
function nmt_repos_end(varargin)
set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');

