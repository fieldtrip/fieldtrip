function nmt_tfplot(axh,twin,fwin,tf)
% nutmegtrip uses this helper function to plot time-frequency data

% reformat beam.bands and corresponding tf info to include gaps
% (NaN?)
t = unique(twin(:));
f = fwin';
f = f(:);

jj = 1;
for ii=1:size(tf,2)
    z(jj,:) = tf(:,ii);
    z(jj+1,:) = tf(:,ii);
    jj = jj+2;
end

global st
st.nmt.gui.h_tf = mesh(axh,t,f,z);
view(axh,2); % '2-D view' of spectrogram
%set(st.nmt.gui.h_tf,'LineStyle','none'); % useful if plot made with 'surf'

if(verLessThan('matlab','8.4')) % necessary to preserve colormap on functional image for Matlab R2014a and earlier
    for ii=1:3
        freezeColors(st.vols{1}.ax{ii}.ax);
    end
end

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