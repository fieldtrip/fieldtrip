function nmt_timeselect(op)
global st

if(~exist('op','var'))
    op = '';
end

%axsel = st.nmt.gui.ax_ts(1);
axsel = gca; % for multi panel support

switch(op)
    case 'ts'
        % if this was an interactive button-press in time series
        currpt = get(axsel,'CurrentPoint');
        currtimept = currpt(1,1);
        currfreqpt = currpt(1,2);
        
        if strcmpi(get(gcbf,'SelectionType'),'alt') % right-click
            st.nmt.cfg.time_idx(2) = dsearchn(st.nmt.time',currtimept);
        else % left-click
            st.nmt.cfg.time_idx(1) = dsearchn(st.nmt.time',currtimept);
            st.nmt.cfg.time_idx(2) = dsearchn(st.nmt.time',currtimept);
            
            st.nmt.cfg.freq_idx = find(diff(st.nmt.freq(:,:)>currfreqpt,[],2));
            if(isempty(st.nmt.cfg.freq_idx))
                st.nmt.cfg.freq_idx = [1 1];
            else
                st.nmt.cfg.freq_idx(2) = st.nmt.cfg.freq_idx;
            end
        end
    case 'textbox'
        t1 = str2num(get(st.nmt.gui.t1,'String'));
        t2 = str2num(get(st.nmt.gui.t2,'String'));
        st.nmt.cfg.time_idx(1) = dsearchn(st.nmt.time',t1);
        st.nmt.cfg.time_idx(2) = dsearchn(st.nmt.time',t2);

        if(isfield(st.nmt,'freq'))
            if(1)
                st.nmt.cfg.freq_idx(1) = get(st.nmt.gui.f1,'Value');
                st.nmt.cfg.freq_idx(2) = get(st.nmt.gui.f2,'Value');
                
                if(st.nmt.cfg.freq_idx(1) > st.nmt.cfg.freq_idx(2))
                    st.nmt.cfg.freq_idx(1) = st.nmt.cfg.freq_idx(2);
                    set(st.nmt.gui.f1,'Value',st.nmt.cfg.freq_idx(1));
                end
            else
                f1 = str2num(get(st.nmt.gui.f1,'String'));
                f2 = str2num(get(st.nmt.gui.f2,'String'));
                st.nmt.cfg.freq_idx(1) = dsearchn(st.nmt.freq(:,1),f1);
                st.nmt.cfg.freq_idx(2) = dsearchn(st.nmt.freq(:,2),f2);
            end
        end
    otherwise
        % otherwise, time interval was already specified, so nothing to do
end

if(st.nmt.cfg.time_idx(2) < st.nmt.cfg.time_idx(1))
    % enforce that time1 < time2
    st.nmt.cfg.time_idx = st.nmt.cfg.time_idx([2 1]);
end

nmt_spm_plot
for funidx=1:length(st.nmt.fun);
    nmt_update_panel(funidx)
end