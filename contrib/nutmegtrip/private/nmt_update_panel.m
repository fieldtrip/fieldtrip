function nmt_update_panel(axsel)

global st
%% update time series, if applicable
if(isfield(st.nmt,'time')) %& ~isfield(st.nmt,'freq'))
    set(st.nmt.gui.timeguih,'Visible','On'); % ensure plot is visible
    switch(st.nmt.cfg.plottype)
        case 'tf'
            clim = max(abs(st.nmt.fun{axsel}(:))) * [-1 1];
            set(st.nmt.gui.freqguih,'Visible','On'); % ensure plot is visible
            if(all(st.nmt.freq(1,:) == [0 0])) % if the first row contains evoked data
                st.nmt.gui.h_tf = nmt_tfplot(st.nmt.gui.ax_ts(axsel),st.nmt.time,st.nmt.freq(2:end,:),squeeze(st.nmt.fun{axsel}(st.nmt.cfg.vox_idx,:,2:end)),clim,@nmt_repos_start);
            else
                st.nmt.gui.h_tf = nmt_tfplot(st.nmt.gui.ax_ts(axsel),st.nmt.time,st.nmt.freq,squeeze(st.nmt.fun{axsel}(st.nmt.cfg.vox_idx,:,:)),clim,@nmt_repos_start);
            end
        case 'ts'
            if(isfinite(st.nmt.cfg.vox_idx))
                plot(st.nmt.gui.ax_ts(axsel),st.nmt.time,squeeze(st.nmt.fun{axsel}(st.nmt.cfg.vox_idx,:,:)));
                grid(st.nmt.gui.ax_ts(axsel),'on');
            else
                plot(st.nmt.gui.ax_ts(axsel),st.nmt.time,nan(length(st.nmt.time),1));
            end
            
            % set y-axis range based on whole volume range (single time point), or
            % selected timeslice range
            if(st.nmt.cfg.time_idx(1)==st.nmt.cfg.time_idx(2))
                ymin = min(st.nmt.fun{axsel}(:));
                ymax = max(st.nmt.fun{axsel}(:));
            else
                ymin =  min(min(st.nmt.fun{axsel}(:,st.nmt.cfg.time_idx(1):st.nmt.cfg.time_idx(2))));
                ymax =  max(max(st.nmt.fun{axsel}(:,st.nmt.cfg.time_idx(1):st.nmt.cfg.time_idx(2))));
            end
            set(st.nmt.gui.ax_ts(axsel),'YLim',[ymin ymax]);
            
    end
    set(st.nmt.gui.ax_ts(axsel),'XLim',st.nmt.time([1 end]));
    xlabel(st.nmt.gui.ax_ts(axsel),['Time (s)']);
    
    
    %% plot vertical line indicating selected time point
    switch(st.nmt.cfg.plottype)
        case 'ts'
            ylim=get(st.nmt.gui.ax_ts(axsel),'YLim');
            facealpha = 1; % make patch opaque
        case 'tf'
            ylim = [st.nmt.freq(st.nmt.cfg.freq_idx(1),1) st.nmt.freq(st.nmt.cfg.freq_idx(2),2)];
            facealpha = 0; % make patch see-through
    end
    axes(st.nmt.gui.ax_ts(axsel));
    zlim = get(st.nmt.gui.ax_ts(axsel),'ZLim'); % selection patch needs to have a Z higher than the TF range so that it's not hidden by the TF plot
    h=patch([st.nmt.time(st.nmt.cfg.time_idx(1)) st.nmt.time(st.nmt.cfg.time_idx(2)) st.nmt.time(st.nmt.cfg.time_idx(2)) st.nmt.time(st.nmt.cfg.time_idx(1))]',[ylim(1) ylim(1) ylim(2) ylim(2)]',[zlim(2) zlim(2) zlim(2) zlim(2)],[1 0.4 0.4],'EdgeColor','red','FaceAlpha',facealpha);

    set(st.nmt.gui.ax_ts(axsel),'ButtonDownFcn',@nmt_repos_start);
    
    
    switch('disabled')  % TODO: hook for overlaying time series on TF plot;
        %       perhaps better solution is independent nmt_spm_plot call for each funparameter
        case 'tf'
            if(st.nmt.cfg.evokedoverlay)
                % trick to overlay time series taken from plotyy.m
                ylim=get(st.nmt.gui.ax_ts(axsel,2),'YLim');
                ts=squeeze(st.nmt.fun{axsel}(st.nmt.cfg.vox_idx,:,1));
                
                axes(st.nmt.gui.ax_ts(axsel,2));
                plot(st.nmt.gui.ax_ts(axsel,2),st.nmt.time,ts);
                set(st.nmt.gui.ax_ts(axsel,2),'YAxisLocation','right','Color','none', ...
                    'XGrid','off','YGrid','off','Box','off', ...
                    'HitTest','off','Visible','on');
            end
    end

    %% update GUI textboxes
    set(st.nmt.gui.t1,'String',num2str(st.nmt.time(st.nmt.cfg.time_idx(1))));
    set(st.nmt.gui.t2,'String',num2str(st.nmt.time(st.nmt.cfg.time_idx(2))));

    if(1) % this method could theoretically allow selection of multiple bands
        set(st.nmt.gui.f1,'Value',st.nmt.cfg.freq_idx(1));
        set(st.nmt.gui.f2,'Value',st.nmt.cfg.freq_idx(2));
    else % this method allows only single bands
        set(st.nmt.gui.f1,'String',num2str(st.nmt.freq(st.nmt.cfg.freq_idx(1),1)));
        set(st.nmt.gui.f2,'String',num2str(st.nmt.freq(st.nmt.cfg.freq_idx(2),2)));
    end
        
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

