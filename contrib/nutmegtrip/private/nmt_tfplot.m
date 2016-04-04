function nmt_tfplot(axh,twin,fwin,tf)

if(0)
    handles = guidata(rivets.fig);
    rivets.tf_ax = handles.nut_ts_axes;
    
    axes(rivets.tf_ax);
    tmp=unique(beam.bands(:));
    set(rivets.tf_ax,'YTick',tmp);
    if ~isinteger(tmp)
        tmp=round(tmp*10)/10;
    end
    set(rivets.tf_ax,'YTickLabel',num2str(tmp));
    if get(handles.nut_check_logscale,'Value')
        set(rivets.tf_ax,'YScale','log')           %,'YTick',[4 12 30 60 120 180 220]);
    else
        set(rivets.tf_ax,'YScale','linear')
    end
    
    tfh = zeros(rivets.timeselect,rivets.freqselect);
    
    if(size(beam.s{1},1)==1)
        beam.s{1} = beam.s{1}';
    end
    
    cla(rivets.tf_ax);
end





%%
% reformat beam.bands and corresponding tf info to include gaps
% (NaN?)
t = unique(twin(:));
f = unique(fwin(:));

[fnogaps,fnogaps_idx] = intersect(f,fwin(:,1)); % find gaps in frequency

% create TF matrix; gaps remain filled with "k", to code points with no data
% surf breaks with Nan or Inf so we pick realmax
k = realmax;
z = k*ones(length(f),length(t));
z(fnogaps_idx,:) = tf';

global st
st.nmt.gui.h_tf = mesh(axh,t,f,z);
view(axh,2); % '2-D view' of spectrogram
%set(st.nmt.gui.h_tf,'LineStyle','none'); % useful if plot made with 'surf'
caxis(axh,[st.vols{1}.blobs{1}.min st.vols{1}.blobs{1}.max*33/32]);
colormap(axh,[jet(128); 1 1 1]);
% set(axh,'YScale','log');
% xlabel(axh,'Time');
ylabel(axh,'Frequency (Hz)');


%%




if(0)        
    
    time = str2double(get(handles.nut_time_text,'String'));
    rivets.timeselect = dsearchn(beam.timepts,time);
    
    %% hack to bring selected patch to foreground with complete red box
    delete(tfh(rivets.timeselect,rivets.freqselect));
    for timebin=rivets.timeselect
        for freqbin=rivets.freqselect
            tfh(timebin,freqbin) = patch([beam.timewindow(timebin,1) beam.timewindow(timebin,1) beam.timewindow(timebin,2) beam.timewindow(timebin,2)],[beam.bands(freqbin,[1 2]) beam.bands(freqbin,[2 1])],rivets.s(MEGvoxelindex,timebin,freqbin),'Parent',rivets.tf_ax,'HitTest','off');
        end
    end
    set(tfh(rivets.timeselect,rivets.freqselect),'EdgeColor',[1 0 0],'LineStyle','-','LineWidth',2);
    %%
    
    axis(rivets.tf_ax,[min(beam.timewindow(:)) max(beam.timewindow(:)) 0 max(beam.bands(:))]);
    caxis(rivets.tf_ax,[rivets.scalemin rivets.scalemax]);
    rivets.cbar = colorbar('peer',rivets.tf_ax);
    if ( isfield(beam,'labels') && isfield(beam.labels,'colorbar') )
        cblab = beam.labels.colorbar;
    else
        cblab = 'dB';
    end
    % dbH = title(rivets.cbar,cblab,'FontSize',14,'FontWeight','bold');
    set(rivets.tf_ax,'YDir','normal');
    set(rivets.tf_ax,'ylim',[beam.bands(1,1) beam.bands(end,end)])
    drawnow;        % This avoids stupid matlab bug in version 2003a
    axes(rivets.tf_ax);
    if isfield(beam,'labels')
        xlabel(beam.labels.xaxis)
        ylabel(beam.labels.yaxis)
    else
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
    end
end