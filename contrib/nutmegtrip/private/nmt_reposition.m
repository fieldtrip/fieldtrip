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
    otherwise
        % otherwise, time interval was already specified, so nothing to do
end

if(st.nmt.cfg.time_idx(2) < st.nmt.cfg.time_idx(1))
    % enforce that time1 < time2
    st.nmt.cfg.time_idx = st.nmt.cfg.time_idx([2 1]);
end


nmt_spm8_plot

