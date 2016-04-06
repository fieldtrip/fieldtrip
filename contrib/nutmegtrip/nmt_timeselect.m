function nmt_timeselect(op)
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

        if(isfield(st.nmt,'freq'))
            f1 = str2num(get(st.nmt.gui.f1,'String'));
            f2 = str2num(get(st.nmt.gui.f2,'String'));
            st.nmt.cfg.freq_idx(1) = dsearchn(st.nmt.freq(:,1),f1);
            st.nmt.cfg.freq_idx(2) = dsearchn(st.nmt.freq(:,1),f2);
        end
    otherwise
        % otherwise, time interval was already specified, so nothing to do
end

if(st.nmt.cfg.time_idx(2) < st.nmt.cfg.time_idx(1))
    % enforce that time1 < time2
    st.nmt.cfg.time_idx = st.nmt.cfg.time_idx([2 1]);
end


nmt_spm8_plot

