function nmt_addtopo(cfgplot,erp)
% configures and display topoplot on nutmegtrip viewer
% (not intended for use outside of nutmegtrip plot functions)

global st
figure(st.fig);
if(verLessThan('matlab','8.4')) % necessary to preserve colormap on functional image for Matlab R2014a and earlier
    for ii=1:3
        freezeColors(st.vols{1}.ax{ii}.ax);
    end
end

newway = true
if(newway)
    axes(st.nmt.gui.ax_ts(2,1));
    set(st.nmt.gui.ax_ts(2,1),'Visible');
    cfgplot.interactive = 'no';
    % cfgplot.colorbar = 'yes';  % how can we transfer colorbar over?
    cfgplot.comment = 'no';
    ft_topoplotER(cfgplot,erp);
    
    set(st.nmt.gui.ax_ts(2,1),'Units','pixels','Position',st.nmt.gui.ax_topo_pos);
    colormap(st.nmt.gui.ax_ts(2,1),jet)

else
    if(isfield(st.nmt.gui,'ax_topo'))
        delete(st.nmt.gui.ax_topo);
    end
    topofig = figure('Visible','off');
    cfgplot.interactive = 'no';
    % cfgplot.colorbar = 'yes';  % how can we transfer colorbar over?
    cfgplot.comment = 'no';
    ft_topoplotER(cfgplot,erp);
    
    st.nmt.gui.ax_topo = copyobj(gca,st.fig);
    delete(topofig);
    
    set(st.nmt.gui.ax_topo,'Units','pixels','Position',st.nmt.gui.ax_topo_pos);
    colormap(st.nmt.gui.ax_topo,jet)
end
end