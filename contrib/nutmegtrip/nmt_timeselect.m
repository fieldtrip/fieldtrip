function nmt_reposition(op)
global st

if(~exist('op','var'))
    op = '';
end

switch(op)
    case 'ts'
        % if this was an interactive button-press in time series
        currpt = get(st.nmt.gui.ax_ts,'CurrentPoint');
        currpt = currpt(1,1);
        
        if strcmpi(get(gcbf,'SelectionType'),'alt') % right-click
            st.nmt.cfg.time_idx(2) = dsearchn(st.nmt.time',currpt);
        else % left-click
            st.nmt.cfg.time_idx(1) = dsearchn(st.nmt.time',currpt);
            st.nmt.cfg.time_idx(2) = dsearchn(st.nmt.time',currpt);
        end
    case 'textbox'
        t1 = str2num(get(st.nmt.gui.t1,'String'));
        t2 = str2num(get(st.nmt.gui.t2,'String'));
        
        st.nmt.cfg.time_idx(1) = dsearchn(st.nmt.time',t1);
        st.nmt.cfg.time_idx(2) = dsearchn(st.nmt.time',t2);
    otherwise
        % otherwise, time interval was already specified, so nothing to do
end

if(st.nmt.cfg.time_idx(2) < st.nmt.cfg.time_idx(1))
    % enforce that time1 < time2
    st.nmt.cfg.time_idx = st.nmt.cfg.time_idx([2 1]);
end


nmt_spm8_plot

