function nmt_peaksearch_helper
% simplifies Callback for peak search functionality in GUI
global st

cfg=[];
cfg.searchradius = [str2num(get(st.nmt.gui.searchradius1,'String')) str2num(get(st.nmt.gui.searchradius2,'String'))];



cfg.peaktype=get(st.nmt.gui.peaktype,'string');
cfg.peaktype=cfg.peaktype{get(st.nmt.gui.peaktype,'Value')};

peakdomain=get(st.nmt.gui.peakdomain,'string');
peakdomain=peakdomain{get(st.nmt.gui.peakdomain,'Value')};

switch(peakdomain)
    case 'spatial'
        cfg.time = st.nmt.cfg.time_idx;
    case 'temporal'
        cfg.vox = st.nmt.cfg.vox_idx;
    case 'spatiotemporal'
        % nothing to do, this is default behavior
    otherwise
        error('well this is unexpected...')
end

[v,t]=nmt_peaksearch(cfg);

nmt_repos(v,t);